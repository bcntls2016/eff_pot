!............................................................
!...                      Subroutine vareps               ...
!............................................................

! Number of particules in the impurity = 1
! It is supposed that the density is already normalized....

!WARNING:
! After the use of 'evolox' what acually is in the array
! hpsix is 'H*PSI-eps*PSI'



subroutine vareps(epsx,errepsx)

use field
use impur
use grid

implicit none
real    (kind=8), intent(IN)  :: epsx    ! (IN) Autovalue
real    (kind=8), intent(OUT) :: errepsx ! (OUT) Error...


!errepsx =  sqrt(sum((psix*(hpsix-epsx*psix))**2))
!errepsx =  sqrt(sum( (hpsix-epsx*psix)**2 ) *dxyz )

!errepsx  = sum( hpsix**2 ) * dxyz
!errepsx  = (errepsx/nxyz) - epsx**2
!errepsx  = sqrt(abs(errepsx))
errepsx  = sqrt( sum( hpsix**2 ) * dxyz  )

return

end
