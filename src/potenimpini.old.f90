subroutine potenimpini
use interpol !, only:DelInter,potpi,potdel,npot,vpi,delta
use grid
use classicimp, only: rimp
implicit none
integer (kind=4) :: i
real    (kind=8) :: r,rmax,V_Pi,V_Sig
real    (kind=8) :: r1,r2


npot =100*max(nx,ny,nz)+1
rmax = 4.d0*max(xmax,ymax,zmax)

DelInter = rmax/dfloat(npot-1)

!.....................................
! Interpolation for V_Pi

! allocate(rpot(npot))
allocate(potpi(0:npot+1))
!allocate(potpi(npot))
allocate(potdel(0:npot+1))
!allocate(potdel(npot))
allocate(  Vpi(nx,ny,nz))
allocate(Delta(nx,ny,nz))

do i=1,npot
 r = dfloat(i-1)*DelInter
!  rpot(i) = r
 potpi(i) = V_Pi(r)
! if(r.Eq.0.0d0)then
   potdel(i) = V_Sig(r) - V_Pi(r)
!   potdel(i) = (V_Sig(r) - V_Pi(r))/r**2
! else
!   potdel(i)=0.d0
! endif
enddo

potpi(0)=potpi(1); potpi(npot+1)=potpi(npot)
potdel(0)=potdel(1); potdel(npot+1)=potdel(npot)

!.....................................


! ! ! ! ! ! ! ! ! ! ! TEST MODEdd
! do r=0.01d0,40.d0,0.01d0
!    r1   =  potpi(int(r/DelInter)+1)*(1.d0-mod(r,DelInter)/DelInter) +  potpi(int(r/DelInter)+2)*mod(r,DelInter)/DelInter
!    r2 = potdel(int(r/DelInter)+1)*(1.d0-mod(r,DelInter)) + potdel(int(r/DelInter)+2)*mod(r,DelInter)
! !    write(*,*)r,r1,r2,V_pi(r),(V_Sig(r)-V_Pi(r))/r**2
!    write(*,*)r,r1,V_pi(r)
! enddo
! STOP


lastr = (/1000.,1000.,1000./)

call updatepoten(rimp)

end subroutine potenimpini



! suboroutine potinter







