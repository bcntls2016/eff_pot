!............................................................
!...                      Subroutine energy               ...
!............................................................

!
! In order to have actualized arrays. It is very convient
! that before to call this subroutine calls the POTEN
! subroutine (for delj4)  and derden (for derivadores of 
! the helium4 density)subroutines


subroutine energy()

use para_derivnd
use alphasterm
use deriva
use energies
use field
! use impur
use classicimp
use lenard4
use grid
use gridk
use he4
use rho
use util1
use work1

implicit none

real (kind=8)    ::  aux1,aux2,aux3,aux4,aux5 ! Auxiliar variables
complex (kind=8) :: invars(10), vlsaux
integer (kind=4) :: ix,iy,iz,i,j, np=6, nd=10



!................................ Use derivatives for (grad(den))**2
!
!  icon   =  0 ! Take the derivative.
!  icon   = 20 ! Take the derivative. Use Wigner-Seitz conditions.
!  icon   =  8 ! Take the derivative. Use periodic conditions.

call derivnd(1,nn,hx,1,psi,sto1c,icon)
call derivnd(1,nn,hy,2,psi,sto2c,icon)
call derivnd(1,nn,hz,3,psi,sto3c,icon)

!call pdergc(1,npd,nn,hx,psi,sto1c,3,1,mmx,iw,icon)   ! derivative respect X
!call pdergc(1,npd,nn,hy,psi,sto2c,3,2,mmx,iw,icon)   ! derivative respect Y
!call pdergc(1,npd,nn,hz,psi,sto3c,3,3,mmx,iw,icon)   ! derivative respect Z


aux1   = 0.5d0*cp4
aux2   = (1.0d0/3.0d0)*cpp4

!..................................................................
!... Calculate the density of energy (kinetic, and correlation) ...
!..................................................................

ekin4 = 0.0d0
elj4  = 0.0d0
ecor4 = 0.0d0

do iz=1,nz
   do iy=1,ny
      do ix=1,nx
          aux3  = den(ix,iy,iz)
          aux4  = dencg(ix,iy,iz)
          aux5  = aux3*aux4**2*(aux1+aux2*aux4)
          ekin4 = ekin4 + abs(sto1c(ix,iy,iz))**2+abs(sto2c(ix,iy,iz))**2+abs(sto3c(ix,iy,iz))**2
          elj4  = elj4  + delj4(ix,iy,iz)*aux3   
          ecor4 = ecor4 + aux5
      end do
   end do
end do


!......................................................
!... Calculate the density of energy (alpha_s term) ...
!......................................................

select case(core4)
   case('OTE','OTC')
     ealphas = sum(falfs*( dxden*intxalf + dyden*intyalf + dzden*intzalf) )
     ealphas = -h2o2m4*0.5d0*alphas*ealphas*dxyz
   case default
     continue
end select

!....................!
!... penalty term ...!
!....................!
esolid = C*sum(den*(1.d0+dtanh(beta*(den-den_m))))*dxyz

ekin4 = h2o2m4*ekin4*dxyz           ! TOTAL Kinetic energy for 4He
elj4  = 0.5d0 *elj4 *dxyz           ! TOTAL Lennard-Jones energy
ecor4 =        ecor4*dxyz           ! TOTAL Correlation energy for 4He
etot4 = ekin4+elj4+ecor4            ! TOTAL ENERGY without impurity

select case(core4)
   case('OTE','OTC')
     etot4 = etot4+ealphas          ! TOTAL ENERGY including Alpha_s term
   case default
     continue
end select

etot4    =  etot4 + sum(uext*den)*dxyz

! Classic vecotrial particle energy:
 ekinx = 0.5d0*mAg_u*sum(vimp*vimp)
 eHeX = sum(uimp*den)*dxyz

!
! Spin-Orbit contribution:
!
If(Lstate.Eq.'P')Then
 invars(:)=conjg(invar(:))
 eso = Also2*( 2.d0*real(invars(1)*invar(6)-invars(2)*invar(5))                                     + & 
               2.d0*aimag(invars(1)*invar(3) + invars(4)*invar(2) + invars(3)*invar(6) + invars(4)*invar(5) ) )
Else
!
!  La matriu SOD s'omple dintre de la subrutina Instates
!
  do i = 1,nd
    invars(i) = (0.d0,0.d0)
    do j = 1,nd
      invars(i) = invars(i) + SOD(i,j)*invar(j)
    end do
  end do
  vlsaux=(0.d0,0.d0)
  do i = 1,nd
    vlsaux = vlsaux + invars(i)*conjg(invar(i))
  end do  
  vlsaux=vlsaux*0.5d0 
  eso=vlsaux*Als
EndIf

! Change while in Linz: The order of the aimag items is wrong (3 1 instead of 1 3 etc).
! New Change in Barcelona: it was ok the first time!

eimpu = ekinx + eHeX + eso

etot   =  etot4 + eimpu


return

end
