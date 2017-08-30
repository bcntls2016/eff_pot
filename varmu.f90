!............................................................
!...                      Subroutine varmu                ...
!............................................................
!
! It is supposed that the density is already normalized.
! Warning!!!
! After the use of 'evolo' what actually is inthe array 'hpsi'
! is H*Psi-Mu*Psi


subroutine varmu(npart,mu,errmu)

use field
use grid
use rho

implicit none
real    (kind=8), intent(IN)  :: mu     ! (IN)  Chemical Potential
real    (kind=8), intent(OUT) :: errmu  ! (OUT) Error
integer (kind=4), intent(IN)  :: npart  ! (IN)  Number of particles

!errmu =  sqrt(sum((psi*(hpsi-mu*psi))**2)) / npart
!errmu =  sqrt(sum((hpsi-mu*psi)**2)*dxyz)/ npart

!errmu = sum( hpsi**2 ) * dxyz
!errmu = errmu/(nxyz*npart) - mu**2
!errmu = sqrt(abs(errmu))

errmu =  sqrt( (sum(hpsi**2)*dxyz))/  npart 

return
end
