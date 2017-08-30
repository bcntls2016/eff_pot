!------------------------------------------------------------------
!---                    Subroutine readenc                      ---
!------------------------------------------------------------------

subroutine readenc(npart,densat,fileden,fileimp,mode,rimpur,r_clust)

! INPUT QUANTITIES:
!
! From list:
!
! npart   ----> Number of particles
! densat  ----> Density of saturation (useful for fermi distributions)
! fileden ----> Name of the file with the density
! mode    ----> Select the density:
!                       0 = continue a calculation
!                       1 = Build from scratch a new density
!                       2 = Build from scratch a new density with impurity
! rimpur  ----> Radius of the impurity
!
!
! OUTPUT QUANTITIES:
!
!  psi    ----> He Wave function)
!  den    ----> Density (module rho)
!  r_clust ---> Size of the cluster
!  psix ------> Wave function for the impurity
!
!
! NOTE:
! This subroutine check the consistency of namelist input data when an 
! external density is used as a input. The quantities to check consistency
! came from modules
!-------------------------------------------------------------------------------

use rho
use field
use grid
! use impur
use classicimp
use util1
use interpol , only:Delta

implicit none
integer   (kind=4), intent(in)    :: npart
real      (kind=8)                :: nph     ! Num. part in hole
real      (kind=8), intent(in)    :: densat
character (len=60), intent(in)    :: fileden,fileimp
integer   (kind=4), intent(in)    :: mode
real      (kind=8), intent(out)   :: r_clust
real      (kind=8), intent(in)    :: rimpur

real      (kind=8) :: xmaxp,ymaxp,zmaxp,xcp,ycp,zcp,hxp,hyp,hzp
real      (kind=8) :: aux1, aux2, aux3,aux
real      (kind=8) :: aux1b,aux2b,aux3b,vls
! real      (kind=8) :: ximp,yimp,zimp

integer   (kind=4) :: nxp,nyp,nzp,nninvar0
integer   (kind=4) :: ix,iy,iz,isalto
logical :: limp
real      (kind=8) :: rr

real      (kind=8) :: Uzr ! Convolution of rho with (V_zz - V_rr)
real      (kind=8) :: AUU ! = Als + 2(Uz-Ur)
real      (kind=8) :: lambda,norm

real      (kind=8) , allocatable  :: z2r2(:,:,:) ! (z-zimp)**2- (x-ximp)**2
Complex   (kind=8), Allocatable :: iinvar0(:)

! real      (kind=8)    :: Oneosq3,Oneosq2! = dsqrt(1.d0/3.d0)
! real      (kind=8)    :: tre,fiv,norm,rr
! real      (kind=8) , parameter :: p1=0.573585,p2=0.584807
! parameter(Oneosq3 = dsqrt(1.d0/3.d0))
! parameter(Oneosq2 = dsqrt(1.d0/2.d0))
! parameter(tre=0.297044d0)
! parameter(fiv=0.495074d0)
!...........................................
!With 'mode' select the kind of density...
!...........................................
!

r_clust=0.d0

select case(mode)
!-------------------------------------------------------------------
  case(0)   ! Continue a calculation (read a previous density)
!-------------------------------------------------------------------

     open(unit=1,file=fileden,status='old')
     call titols(1,cchar,isalto)
     read(1,*) xmaxp,ymaxp,zmaxp,hxp,hyp,hzp,nxp,nyp,nzp,limp,ximp,yimp,zimp
     read(1,*) den
     close(1)
 
      rimp = (/ximp,yimp,zimp/)


!---------------------------------------------------------------------------------------
 case(1)   ! Continue a calculation (read a previous wave functon for non-excited state)
!---------------------------------------------------------------------------------------
   open(unit=1,file=fileden,status='old')
     call titols(1,cchar,isalto)
     read(1,*) xmaxp,ymaxp,zmaxp,hxp,hyp,hzp,nxp,nyp,nzp!,ximp,yimp,zimp,vximp,vyimp,vzimp
     read(1,*) rimp
     read(1,*) vimp
     !read(1,*) ninvar
     !read(1,*) invar
     read(1,*) psi
     den=Conjg(psi)*psi
     ximp = rimp(1)
     yimp = rimp(2)
     zimp = rimp(3)
   close(1)



!---------------------------------------------------------------
 case(2)   ! Continue a calculation (read a previous wave functon)
!---------------------------------------------------------------
   open(unit=1,file=fileden,status='old')
     call titols(1,cchar,isalto)
     read(1,*) xmaxp,ymaxp,zmaxp,hxp,hyp,hzp,nxp,nyp,nzp!,ximp,yimp,zimp,vximp,vyimp,vzimp
     read(1,*) rimp
     read(1,*) vimp
     read(1,*) nninvar0
     Allocate (iinvar0(nninvar0))
     read(1,*) iinvar0
     read(1,*) psi
     den=Conjg(psi)*psi
    close(1)
    ximp = rimp(1)
    yimp = rimp(2)
    zimp = rimp(3)
    Call Updatepoten(rimp)
    Lprint_invar=.false.
!
!  Per emplenar les matrius del Spi-Orbita
!
    call instates()

     If(nninvar0.Lt.10.Or.LState.Eq.'P')Then
       invar = (0.d0, 0.0d0)
       Do ix=1, 6
         invar(ix) = iinvar0(ix)
       Enddo
     Else
       invar = iinvar0
     Endif
     Deallocate(iinvar0)
    Return
!---------------------------------------------------------------
 case(3)   ! Continue a calculation (read a previous wave functon)
         !---------------------------------------------------------------
            open(unit=1,file=fileden,status='old')
            call titols(1,cchar,isalto)
            read(1,*) xmaxp,ymaxp,zmaxp,hxp,hyp,hzp,nxp,nyp,nzp!,ximp,yimp,zimp,vximp,vyimp,vzimp
            read(1,*) rimp
            read(1,*) vimp
            read(1,*) nninvar0
            Allocate (iinvar0(nninvar0))
            read(1,*) iinvar0
            read(1,*) psi
            den=Conjg(psi)*psi
            close(1)
            ximp = rimp(1)
            yimp = rimp(2)
            zimp = rimp(3)
            Call Updatepoten(rimp)
     Deallocate(iinvar0)

!
!-------------------------------------------------
 case default ! Non programmed situation
!-------------------------------------------------
!
      stop 'ERROR READEN: This mode is still not programmed'
end select


Lprint_invar=.true.

!  PREPARE INITIAL STATE FOR INVAR
if(instate>=0)then
  

  Call Updatepoten(rimp)

  select case(instate)
!-------------------------------------------------
  case(0) ! Pi 1/2 1/2 (<V_SO> approx -Als)  If(State='P')
          ! D  3/2 1/2                       If(State='D')
!-------------------------------------------------
    
      call instates()

!-------------------------------------------------
  case(1) ! Pi 3/2 3/2 (<V_SO> approx +Als/2) If(State='P')
          ! D  3/2 3/2                        If(State='D')
!-------------------------------------------------
      call instates()

!-------------------------------------------------
  case(2) ! Pi 3/2 1/2 (<V_SO> equals +Als/2) If(State='P')
          ! D  5/2 1/2                        If(State='D')
!-------------------------------------------------

      call instates()

!-------------------------------------------------
  case(3) ! Pi 3/2 1/2 (<V_SO> equals +Als/2) If(State='P')
          ! D  5/2 3/2                        If(State='D')
!-------------------------------------------------

      call instates()

!-------------------------------------------------
  case(4) ! Pi 3/2 1/2 (<V_SO> equals +Als/2) If(State='P')
          ! D  5/2 5/2                        If(State='D')
!-------------------------------------------------

      call instates()

!-------------------------------------------------
  case(5) ! Read from invar.dat
!-------------------------------------------------
   open(unit=1,file='invar.dat',status='OLD')
     read(1,*) invar(1) 
     read(1,*) invar(2) 
     read(1,*) invar(3) 
     read(1,*) invar(4) 
     read(1,*) invar(5) 
     read(1,*) invar(6) 
     read(1,*) invar(7)
     read(1,*) invar(8)
     read(1,*) invar(9)
     read(1,*) invar(10)
   close(1)

!-------------------------------------------------
  case default ! other
!-------------------------------------------------
    write(*,*)'instate value invalid',instate
    STOP
  end select
   
else
  write(*,*)'instate less than 0, still not programmed',instate
  STOP
endif

 Lprint_invar=.false.

return
end
