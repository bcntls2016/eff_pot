      SUBROUTINE STEPPC(deltat,errHe,errimp,errvimp,erriv)
!
!      Predictor-Modifier-Corrector method
!       (Ralston & Wilf Vol I, pag. 99)
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
use rkpc   ! (Storage for Steprk & Steppc)

implicit none

real (kind=8) :: c112=112.d0/121.d0
real (kind=8) :: c9  =9.d0/121.d0
real (kind=8) :: c1o3=1.d0/3.d0
real (kind=8) :: c4o3=4.d0/3.d0
real (kind=8) :: c5o3=5.d0/3.d0

integer (kind=4) :: ix,iy,iz,iaux
real    (kind=8) :: deltat
real    (kind=8) :: errHe, errimp,erriv, errvimp
real    (kind=8) :: auxr(3)
complex (kind=8) :: auxc(10)
complex (kind=8) :: aux1c,aux2c,aux3c,aux4c
complex (kind=8) :: ci=cmplx(0.0d0,1.0d0)

call derivnd(2,nn,hx,1,psi,sto1c,icon)
call derivnd(2,nn,hy,2,psi,sto2c,icon)
call derivnd(2,nn,hz,3,psi,sto3c,icon)


!
!   We compute H�Psi
!
!-------------------------------------------------------
!-------------------------------------------------------
!write(6,*)" Control setppc(1): Ioldp(1-3)",Ioldp(1),ioldp(2),ioldp(3)
!write(6,*)" Control setppc(1): Ioldh(1-2)",Ioldh(1),ioldh(2)
Call Flush(6)
!-------------------------------------------------------
!-------------------------------------------------------
!
!      Predictor:
!
      sto4c = timec*((sto1c+sto2c+sto3c)*h2o2m4 - pot4*psi) - ci*uimp*psi
      sto1c = psiold(:,:,:,ioldp(3))                               &
            + c4o3*deltat*(2.d0*sto4c - hpsiold(:,:,:,ioldh(1))    &
            + 2.d0*hpsiold(:,:,:,ioldh(2)))
!
!      Modificador:
!
      psiold(:,:,:,ioldp(3)) = psi
      psi = sto1c - c112*pc
      pc  = sto1c
      hpsiold(:,:,:,ioldh(2)) = sto4c
      den = Abs(psi)**2
!
! Aqui reubicamos los indices para no tener que mover las fuciones
!
iaux=ioldh(2)
ioldh(2)=ioldh(1)
ioldh(1)=iaux

!................!
!... position ...!
!................!
! Predictor
 auxr(:) = rimpold(:,ioldr(3)) + c4o3*deltat*(2.d0*vimp(:) - vimpold(:,ioldv(1)) + 2.d0*vimpold(:,ioldv(2)))
! Modificador
 rimpold(:,ioldr(3)) = rimp
 rimp = auxr(:) - c112*pcr(:)
 pcr = auxr

!..................!
!... velocities ...!
!..................!
! Predictor
 auxr(:) = vimpold(:,ioldv(3)) + c4o3*deltat*(2.d0*aimp(:) - aimpold(:,iolda(1)) + 2.d0*aimpold(:,iolda(2)))
! Modificador
 vimpold(:,ioldv(3)) =  auxr(:) - c112*pcv(:)
 pcv = auxr

If(.Not.Lfix_lambda)Then
  !..........................!
  !... Internal variables ...!
  !..........................!
! Predictor
   auxc(:) = invarold(:,ioldiv(3)) + c4o3*deltat*(2.d0*hinvar(:) - hinvarold(:,ioldhiv(1)) + 2.d0*hinvarold(:,ioldhiv(2)))
!  Modificador
   invarold(:,ioldiv(3)) = invar
   invar = auxc(:) - c112*pciv(:)
   pciv = (0.d0, 0.d0)
   pciv = auxc
Endif

 aimpold(:,iolda(2)) = aimp(:)
 Hinvarold(:,ioldhiv(2)) = hinvar(:)
  ! Reubicacion indices
 iaux=iolda(2)  ; iolda(2)=iolda(1)   ; iolda(1)=iaux
 iaux=ioldhiv(2); ioldhiv(2)=ioldhiv(1); ioldhiv(1)=iaux


!........................................
    call potenimp(rimp,invar)
    call poten()
    call forceimp(rimp,aimp)
    aimp(:) = aimp(:)/mAg_u
    call potinvar(rimp,invar,hinvar)
!........................................

call derivnd(2,nn,hx,1,psi,sto1c,icon)
call derivnd(2,nn,hy,2,psi,sto2c,icon)
call derivnd(2,nn,hz,3,psi,sto3c,icon)


errHe=0.0d0
!-------------------------------------------------------
!-------------------------------------------------------
!write(6,*)" Control setppc(2): Ioldp(1-3)",Ioldp(1),ioldp(2),ioldp(3)
!write(6,*)" Control setppc(2): Ioldh(1-2)",Ioldh(1),ioldh(2)
Call Flush(6)
!-------------------------------------------------------
!-------------------------------------------------------
!
!      Seguim amb el modificador:
!
      sto4c = timec*((sto1c+sto2c+sto3c)*h2o2m4 &
                  - pot4*psi) - ci*uimp*psi
!
!     Corrector:
!
      sto5c = 0.125d0*( 9.d0*psiold(:,:,:,ioldp(3)) - psiold(:,:,:,ioldp(2))       &
                       +3.d0*deltat*(sto4c + 2.d0*hpsiold(:,:,:,ioldh(1))          &
                                                - hpsiold(:,:,:,ioldh(2))  ))
      pc = pc - sto5c

!
!     Valor final:
!
      psi = sto5c + c9*pc
      errHe = c9*sum(Abs(pc))
      den = Abs(psi)**2
errHe=errHe/nxyz

!
! Aqui reubicamos los indices para no tener que mover las fuciones
!
      iaux=ioldp(3)
      ioldp(3)=ioldp(2)
      ioldp(2)=ioldp(1)
      ioldp(1)=iaux
!-------------------------------------------------------
!-------------------------------------------------------
!write(6,*)" Control setppc(3): Ioldp(1-3)",Ioldp(1),ioldp(2),ioldp(3)
!write(6,*)" Control setppc(3): Ioldh(1-2)",Ioldh(1),ioldh(2)
! Call Flush(6)
!-------------------------------------------------------
!-------------------------------------------------------

!.................!
!... positions ...!
!.................!
! Corrector:
auxr(:) = 0.125d0*( 9.d0*rimpold(:,ioldr(3)) - rimpold(:,ioldr(2))     &
                 +3.d0*deltat*(vimpold(:,ioldv(3)) + 2.d0*vimp(:) - vimpold(:,ioldv(1)) ))
! vpold3 is actually the v_temporal just computed, so it is the 'newest'.
! The combination of v and vold is different because THE INDEXS HAVE NOT BEEN REALLOCATED YET.
pcr = pcr -auxr
! Valor final:
  rimp = auxr + c9*pcr
errimp = sum(Abs(c9*pcr))*0.3333333333d0
! Reubicacion
iaux=ioldr(3) ; ioldr(3)=ioldr(2) ; ioldr(2)=ioldr(1) ; ioldr(1)=iaux

!..................!
!... velocities ...!
!..................!
! Corrector:
auxr(:) = 0.125d0*( 9.d0*vimp(:) - vimpold(:,ioldv(2))     &
                 +3.d0*deltat*(aimp(:) + 2.d0*aimpold(:,iolda(1)) - aimpold(:,iolda(2)) ))
pcv = pcv -auxr
vimpold(:,ioldv(3)) = vimp(:)
! Valor final:
  vimp = auxr + c9*pcv
 errvimp = sum(Abs(c9*pcv))*0.3333333333d0

If(.Not.Lfix_lambda)Then
  !..........................!
  !... Internal variables ...!
  !..........................!
  ! Corrector:
  auxc(:) = 0.125d0*( 9.d0*invarold(:,ioldiv(3)) - invarold(:,ioldiv(2))     &
                 +3.d0*deltat*(hinvar(:) + 2.d0*hinvarold(:,ioldhiv(1)) - hinvarold(:,ioldhiv(2)) ))
  If(LState.Eq.'P')Then
    pciv(1:6) = pciv(1:6) -auxc(1:6)
  Else
    pciv = pciv -auxc
  Endif
  ! Valor final:
  invar = (0.d0, 0.d0)
  If(LState.Eq.'P')Then
    invar(1:6) = auxc(1:6) + c9*pciv(1:6)
  Else
    invar = auxc + c9*pciv
  Endif
  If(LState.Eq.'P')Then
    erriv = sum(Abs(c9*pciv(1:6)))*0.1666666666666d0
  Else
    erriv = sum(Abs(c9*pciv))*0.1d0
  Endif
Else
  erriv = 0.d0
Endif


! Reubicacion
iaux=ioldv(3) ; ioldv(3)=ioldv(2) ; ioldv(2)=ioldv(1) ; ioldv(1)=iaux
iaux=ioldiv(3) ; ioldiv(3)=ioldiv(2) ; ioldiv(2)=ioldiv(1) ; ioldiv(1)=iaux


return
end
