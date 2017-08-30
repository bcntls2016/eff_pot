!------------------------------------------------------------------
!---                    Subroutine printout                     ---
!------------------------------------------------------------------

!    iprint = 1 ---> Last printout before exit.
!    iprint = 2 ---> Print due a change of Paflov parameter.
!    iprint = 3 ---> Security copy.
!
!  If we use negative values the output files will be written
!  in binary format.


subroutine printoutc(time,iprint,namefile,namefile1,             &
                    psi,nx,ny,nz,hx,hy,hz,       &
                    xmax,ymax,zmax,rimp,vimp,psix,    &
                    paflv,iter,invar,ninvar)

implicit none
! logical                        :: limp
integer   (kind=4), intent(in) :: iter,iprint
integer   (kind=4), intent(in) :: nx,ny,nz,ninvar
complex   (kind=8), intent(in) :: psi(nx,ny,nz)
complex   (kind=8), intent(in) :: psix(nx,ny,nz)
real      (kind=8), intent(in) :: xmax,ymax,zmax,hx,hy,hz
real      (kind=8), intent(in) :: paflv
real      (kind=8), intent(in) :: time
real      (kind=8), intent(in) :: rimp(3)
real      (kind=8), intent(in) :: vimp(3)
! character (len=3)              :: elem
character (len=60)             :: namefile,namefile1
complex   (kind=8) :: invar(ninvar)
!---------------------------------------------------------------
!---------------------------------------------------------------
!---------------------------------------------------------------
real      (kind=8) :: xnorm,xcm,ycm,zcm,auxx,auxy,auxz,sto
integer   (kind=4) :: ix, iy, iz
!
! Before we print the w.f., we will compute the helium C.M.
!
xcm = 0.d0
ycm = 0.d0
zcm = 0.d0
xnorm = 0.d0
Do iz = 1, nz
  auxz = -zmax + (iz-1)*hz
  Do iy = 1, ny
    auxy = -ymax + (iy-1)*hy
    Do ix = 1, nx
      sto = psi(ix,iy,iz)*conjg(psi(ix,iy,iz))
      xnorm = xnorm + sto
      auxx = -xmax + (ix-1)*hx
      xcm = xcm + auxx*sto 
      ycm = ycm + auxy*sto 
      zcm = zcm + auxz*sto
    EndDo  
  EndDo  
EndDo  
xcm = xcm/xnorm
ycm = ycm/xnorm
zcm = zcm/xnorm
!---------------------------------------------------------------
!---------------------------------------------------------------
!---------------------------------------------------------------

!-----------------------------------------------
! If iprint < 0 then print binary files...   ---
!-----------------------------------------------
!
if(iprint.lt.0) then
   open(10,file=namefile,form='UNFORMATTED',BUFFERED='yes')
   write(10) xmax,ymax,zmax,hx,hy,hz,nx,ny,nz
   write(10) rimp
   write(10) vimp
   write(10) ninvar
   write(10) invar
   write(10) psi
   close(10)
!    if(limp) then
!       open(11,file=namefile1,form='UNFORMATTED')
!       write(11) xmax,ymax,zmax,hx,hy,hz,nx,ny,nz,limp,ximp,yimp,zimp
!       write(11) psix
!       close(11)
!    endif
else
!------------------------------------------------------------
! If iprint > 0 then print 'normal files'                 ---
!               the first lines of the file are comments  ---
!               telling you the kind of backup file       ---
!------------------------------------------------------------
!
  open(10,file=namefile,form='FORMATTED',BUFFERED='yes')
!   if(limp) then
!     write(10,1010) iter,paflv,elem
!   else
    write(10,1020) iter,paflv
    Write(10,'("#  Total evolution time(ps)...:",1p,E20.12)')time
    Write(10,'("#  xcm,ycm,zcm:",1p,3E16.8)')xcm,ycm,zcm
!   end if
  select case(iprint)
    case(1)  !........................................ last printout
      write(10,1030)
    case(2)  !............ print due of a change of Paflov parameter
      write(10,1040)
    case(3)  !........................... print due to a backup copy
      write(10,1050)
  end select
!...(Print the density of helium)
  write(10,*) xmax,ymax,zmax,hx,hy,hz,nx,ny,nz
  write(10,*) rimp
  write(10,*) vimp
  write(10,*) ninvar
  write(10,*) invar
  write(10,*) psi
  close(10)

!... If (limp.eq.true) then print the wave function
!   if(limp) then
!      open(10,file=namefile1,form='FORMATTED')
!      write(10,1010) iter,paflv,elem
!      select case(iprint)
!        case(1)  !....................................... last printout
!          write(10,1030)
!        case(2)  ! print due of a change of Paflov parameter
!          write(10,1040)
!        case(3)  ! print due to a backup copy
!          write(10,1050)
!      end select
!      write(10,*) xmax,ymax,zmax,hx,hy,hz,nx,ny,nz,limp,ximp,yimp,zimp
!      write(10,*) psix
!      close(10)
!   end if
end if

return

1010 format('#  Density after ',I15,' iterations.',/,   &
            '#  Paflov factor ',0P,F6.3,/,             &
            '#  Impurity      ',A)
1020 format('#  Density after ',I15,' iterations.',/,   &
            '#  Actual delatatps....:',1P,E20.10)
1030 format('#  Last printout of the run.')
1040 format('#  Change of Paflov parameter.')
1050 format('#  Backup copy.')
end
