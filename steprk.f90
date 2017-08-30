      SUBROUTINE STEPRK(deltat)
!
!      Runge-Kutta-Gill method
!       (Ralston & Wilf Vol I, pag. 117)
!
use para_derivnd
use deriva ! (icon,npd,dxden,dyden,dzden,pderx,pdery,pderz,llap,xlap)
use field  ! (pot4,hpsi,uext,limp)
use grid   ! (nx,ny,nz,nxyz,dxyz)
use gridk  ! (px,py,pz)
use he4    ! (h2o2m4)
! use impur  !
use classicimp
use rho    ! (psi, psiold & hpsiold)
use util1  ! (vdt,nn,mmx,iw)
use work1  ! (temporal storage)
use rkpc   ! (Storage for Steprk & Steppc rutines)

implicit none

real (kind=8) :: arun(4),brun(4),crun(4)

integer (kind=4) :: ix,iy,iz,jrun,i
real    (kind=8) :: deltat
real    (kind=8) :: aux2r(3)
complex (kind=8) :: auxc(10)
complex (kind=8) :: aux1c,aux2c
complex (kind=8) :: ci=cmplx(0.0d0,1.0d0)
Logical   :: Lprint=.false.

arun(1)=0.5d0
arun(2)=1.0d0-1.d0/dsqrt(2.d0)
arun(3)=1.0d0+1.d0/dsqrt(2.d0)
arun(4)=1.d0/6.d0

brun(1)=2.0d0
brun(2)=1.0d0
brun(3)=1.0d0
brun(4)=2.0d0

 crun(1)=0.5d0
 crun(2)=1.0d0-1.d0/dsqrt(2.d0)
 crun(3)=1.0d0+1.d0/dsqrt(2.d0)
 crun(4)=0.5d0


!........................................
!.. Laplacian of H^{iorder-1}�Psi (0) ...
!........................................
!
!   icon   =  0 ! Take the derivative.
!   icon   = 20 ! Take the derivative. Use Wigner-Seitz conditions.
!   icon   =  8 ! Take the derivative. Use periodic conditions.
!

do jrun=1,4

  call derivnd(2,nn,hx,1,psi,sto1c,icon)
  call derivnd(2,nn,hy,2,psi,sto2c,icon)
  call derivnd(2,nn,hz,3,psi,sto3c,icon)

  !call pdergc(2,npd,nn,hx,psi,sto1c,3,1,mmx,iw,icon)
  !call pdergc(2,npd,nn,hy,psi,sto2c,3,2,mmx,iw,icon)
  !call pdergc(2,npd,nn,hz,psi,sto3c,3,3,mmx,iw,icon)



!
!   We compute H�Psi
!
  do iz=1,nz
    do iy=1,ny
      do ix=1,nx
        aux1c = timec(ix,iy,iz)*((sto1c(ix,iy,iz)+sto2c(ix,iy,iz)+sto3c(ix,iy,iz))*h2o2m4 &
                    - pot4(ix,iy,iz)*psi(ix,iy,iz)) - ci*uimp(ix,iy,iz)*psi(ix,iy,iz)
        aux2c = arun(jrun)*(aux1c - brun(jrun)*q(ix,iy,iz))
        q(ix,iy,iz) = q(ix,iy,iz) + 3.*aux2c - crun(jrun)*aux1c
        if(jrun.eq.1)then
          hpsiold(ix,iy,iz,2) = hpsiold(ix,iy,iz,1)
          hpsiold(ix,iy,iz,1) = aux1c
          psiold(ix,iy,iz,3) = psiold(ix,iy,iz,2)
          psiold(ix,iy,iz,2) = psiold(ix,iy,iz,1)
          psiold(ix,iy,iz,1) = psi(ix,iy,iz)
        endif
        psi(ix,iy,iz) = psi(ix,iy,iz) + deltat*aux2c
        den(ix,iy,iz) = Abs(psi(ix,iy,iz))**2
      enddo
    enddo
  enddo

!...............................................................test
! print*,'inside steprk, jrun=',jrun
 Lprint=.false.
 do iz=1,nz ; do iy=1,ny ; do ix=1,nx 
   if(.not.(den(ix,iy,iz).gt.0))then
     write(6,*)ix,iy,iz,den(ix,iy,iz)
     Lprint=.true.
     Go to 10
   endif
 end do ; enddo ; enddo
10 Continue
 If(Lprint)Then
    write(6,*)'inside steprk, jrun=',jrun
    write(6,'("rimp(1,2,3)..",1p,3e15.6)')rimp
    write(6,'("vimp(1,2,3)..",1p,3e15.6)')vimp
    write(6,'("Invar.....",/,1p,(2E15.8))')(invar(i),i=1,ninvar)
    write(6,'("Hinvar....",/,1p,(2E15.8))')(Hinvar(i),i=1,ninvar)
    call respar(x,y,z,nx,ny,nz,2,'uimp','den',uimp,den)
    STOP 'Steprk 01'
 EndIf
!...............................................................test

!
! Impurity evolution if it is necessary
!
!.................!
!... positions ...!
!.................!
  aux2r(:) = arun(jrun)*(vimp(:) - brun(jrun)*qr(:))
     qr(:) = qr(:) + 3.*aux2r(:) - crun(jrun)*vimp(:)
  if(jrun.eq.1)then
!    vpold(:,:,2) = vpold(:,:,1) 
!    vpold(:,:,1) =    vp(:,:)            
   rimpold(:,3) = rimpold(:,2) 
   rimpold(:,2) = rimpold(:,1) 
   rimpold(:,1) = rimp(:) 
  endif

  rimp(:) = rimp(:) + deltat*aux2r(:)

!..................!
!... velocities ...!
!..................!
  aux2r(:) = arun(jrun)*(aimp(:) - brun(jrun)*qv(:))
  qv(:) = qv(:) + 3.*aux2r(:) - crun(jrun)*aimp(:)
  if(jrun.eq.1)then
   aimpold(:,2) = aimpold(:,1) 
   aimpold(:,1) =    aimp(:)            
   vimpold(:,3) = vimpold(:,2) 
   vimpold(:,2) = vimpold(:,1) 
   vimpold(:,1) =    vimp(:) 
  endif
  vimp(:) = vimp(:) + deltat*aux2r(:)

If(.Not.Lfix_lambda)Then
!..........................!
!... Internal Variables ...!
!..........................!
  auxc(:) = arun(jrun)*(Hinvar(:) - brun(jrun)*qiv(:))
     qiv(:) = qiv(:) + 3.*auxc(:) - crun(jrun)*Hinvar(:)
  if(jrun.eq.1)then
    Hinvarold(:,2) = Hinvarold(:,1)
    Hinvarold(:,1) = Hinvar(:)
    invarold(:,3) = invarold(:,2)
    invarold(:,2) = invarold(:,1)
    invarold(:,1) = invar(:)
  endif

  invar(:) = invar(:) + deltat*auxc(:)
Endif

! ...........................................................

  if(jrun.le.3)then
!    Lfirst=.true.
!    Lprint_invar=.true.
    call potenimp(rimp,invar)
    call poten()
    call forceimp(rimp,aimp)
    aimp(:) = aimp(:)/mAg_u
!    Lfirst=.true.
!    Lprint_invar=.true.
    call potinvar(rimp,invar,Hinvar)
  endif
enddo

do ix=1,3
  ioldp(ix)=ix
  ioldr(ix)=ix
  ioldv(ix)=ix
  ioldiv(ix)=ix
enddo

do ix=1,2
  ioldh(ix)=ix
  iolda(ix)=ix
  ioldhiv(ix)=ix
enddo

return
end
