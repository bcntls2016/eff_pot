program DFT4HeImpd

!
!                         ******************
!                         *** DFT4HeImpd ***
!                         ******************
!
!

!
! This program computes time dependent Helium_4 drops with/without impurities using
! Density Functional Theory.
! The interation can be Orsay-Paris or Orsay-Trento
! The code is fully 3-dimensional
!
!-----------------------------------------------------------------------------
! Version 0  (alpha)   Barcelona April-19,2004   R. Mayol & M. Pi
! Version 99 (Epsilon) Barcelona April-19,2006   R. Mayol & M. Pi
!-----------------------------------------------------------------------------
!07oct2004/11:30AM  
!06oct2004/11:30AM
!23sep2004/11:30AM
!Feb2005
! (LE FALTA: Ficheros binarios)
!-----------------------------------------------------------------------------

use seleccio_de_potencial
use Modificacio_De_Select_Pot
use para_derivnd
use alphasterm
use deriva
use energies
use rho
use field
use fftmodule
use grid 
use gridk
use classicimp
use rkpc
! use impur
use lenard4
use he4
use util1
use work1

implicit none

real (kind=4) :: t0,t1,t2,t3,t4,t5,t6   ! Variables used to printout run time execution.

logical              :: lfilepv         ! T-> print save file when change Paflo parameter.
logical      :: L_v_aditional=.false.   ! To add and aditional velocity when we continue a calculation form a previous w.f.
logical              :: lpaflv=.false.  ! T-> allows change of Paflov coeffient
logical              :: lrkpc=.true.    ! T-> allows to use diferent evolution procedures
logical              :: lrk=.false.     ! T-> allows to only Runge-Kutta method
integer    (kind=4)  :: naux            ! Auxiliar variable
integer    (kind=4)  :: nstepp=1        ! Number of 'Paflov parameter'.
integer    (kind=4)  :: ix ,iy ,iz      ! Used as indexs of arrays
integer    (kind=4)  :: ipx,ipy,ipz     ! Used as indexs of arrays
integer    (kind=4)  :: iter,niter      ! Control the number of iterations
integer    (kind=4)  :: pener=50        ! Computes energy each 'pener' iterations
integer    (kind=4)  :: pdenpar=50      ! Writes partial densities each 'pdenpar' iterations
integer    (kind=4)  :: pchem=50        ! Writes partial chemical potential each 'pchem' iter.
integer    (kind=4)  :: ppot=1          ! Computes the potential each 'ppot' iterations
integer    (kind=4)  :: tstgrid
integer    (kind=4)  :: norder          ! Order of Taylor expansion for time evolution
integer    (kind=4)  :: nv_iter         ! Number of iterationsfor the first time evolution
real       (kind=8)  :: n4real          ! Number of He atoms calculated
real       (kind=8)  :: norma           ! Use for normalization
real       (kind=8)  :: mimpur          ! Impurity mass in uma
real       (kind=8)  :: denxmax         ! Maximum value for density of the impurity
real       (kind=8)  :: r_clust         ! Radius of the helium cluster
real       (kind=8)  :: rimpur=5.5      ! Radius of the impurity
real       (kind=8)  :: p,px2,py2,pz2   ! Temporary variables for momentum values
real       (kind=8)  :: mu4,mu4err      ! Value of Chemical potential and associated error
real       (kind=8)  :: epsx,epsxerr    ! Value of autovalue and associated error
real       (kind=8)  :: errmu4          ! Relative change betwen iteration for chemical potential
real       (kind=8)  :: deltat=0.d0     ! Time step
real       (kind=8)  :: deltatps=0.d0   ! Time step in picoseconds
real       (kind=8)  :: xlamda=0.0d0    ! Monopol costraint
real       (kind=8)  :: Ktops=7.63822291d0
real       (kind=8)  :: pstoK=1.d0/7.63822291d0
real       (kind=8)  :: xlamdax=0.0d0    ! He Dipol costraint x direction
real       (kind=8)  :: xlamday=0.0d0    ! He Dipol costraint y direction
real       (kind=8)  :: xlamdaz=0.0d0    ! He Dipol costraint z direction
real       (kind=8)  :: xlamdax_x=0.0d0    ! Impurity Dipol costraint x direction
real       (kind=8)  :: xlamdax_y=0.0d0    ! Impurity Dipol costraint y direction
real       (kind=8)  :: xlamdax_z=0.0d0    ! Impurity Dipol costraint z direction
real       (kind=8)  :: eold            ! Auxiliar variables
real       (kind=8)  :: temps           ! Auxiliar variables
real       (kind=8)  :: aux,aux1,aux2,aux3,time0=-1             ! Auxiliar variables
real       (kind=8)  :: aux4,aux5,aux6                 ! Auxiliar variables
integer    (kind=4),allocatable  :: nitera(:) ! In wich iteration change to the
real       (kind=8),allocatable  :: paflv(:)  !    corresponging Paflov coeffiecient
real       (kind=8)  :: cnorm           !    Normalization constant
real       (kind=8)  :: pmod1        !    Work variable for controlling p-values
real       (kind=8)  :: pmod2        !    Work variable for controlling p-values
real       (kind=8)  :: xcm4,ycm4,zcm4,xcmx,ycmx,zcmx ! Center of mass Drop and Impurity
real       (kind=8)  :: distx,disty,distz  ! Distance between center of masses
real       (kind=8)  :: errHe, errrimp, errvimp, erriv     ! Error evolution (only form Predictor-Corrector-Modificator)
real       (kind=8)  :: Zsurf = -25.d0, Morse_HeTiO2_1D, Morse_HeTiO2_3D, auxn4 ! Position of the surface
real       (kind=8)  :: vimp0(3)    = 0.d0 ! Auxilar vector for adding aditional velocity to the impurity
real       (kind=8)  :: xlambda0(3) = 0.d0 ! Auxilar vector for adding aditional velocity to the drop

integer    (kind=4)  :: n4=300         ! Number of helium_4 atoms
integer    (kind=4)  :: mode=0          ! Way to start the program (see readen subroutine)
integer    (kind=4)  :: iter0=0         ! Starting point for iteration procedure
integer    (kind=4)  :: ncurr,pcurr
integer    (kind=4)  :: icurr = 0
integer    (kind=4)  :: ndmax = 2

character  (len=40)  :: title         = 'Helium4 - 3dim.  '
character  (len=60)  :: fileout       = 'DFT4He3d.res'
character  (len=60)  :: filedenin     = 'he4.input.density'
character  (len=60)  :: filedenout    = 'he4.output.density'
character  (len=60)  :: fileimpin     = 'X.input.wf'
character  (len=60)  :: fileimpout    = 'X.output.wf'
character  (len=60)  :: filerimp      = 'rimp.out'
character  (len=60)  :: filevimp      = 'vimp.out'
character  (len=60)  :: fileaimp      = 'aimp.out'
character  (len=60)  :: filelambda    = 'lambda.out'
character  (len=60)  :: namefile,namefile1
character  (len=23)  :: curvfile
character  (len=3)   :: chariter
logical              :: lsurf=.false.        ! include TiO2 surface or not
logical              :: lsurf3D=.false.        ! include TiO2 surface or not
real       (kind=8)  :: Lambdah,tzmean,tzsurf,rt
integer              :: i
integer (kind=4)	 :: nthread   ! Number of threads
COMPLEX (kind=8)     :: uim = (0.d0,1.d0)

! interface
!   double precision function v_alka(d,elem)
!        character (len=3),  intent(in) :: elem
!        real      (kind=8), intent(in) :: d
!   end function v_alka
! end interface



!....................Variables read in a NAMELIST statement ..............................

namelist /input/title,fftwplan,nthread,nsfiles,                         &
                fileout,filerimp,filevimp,fileaimp,                     &
                filedenin,filedenout,                                   &
                fileimpin,fileimpout,                                   &
                n4,mode,instate,                                        &
                nx,ny,nz,xmax,ymax,zmax,                                &
                xc,yc,zc,afermi,                                        &
                eps4,sigma4,core4,l,                                    &
                cp4,cpp4,den4c,alphas,h2o2m4,                           &
                denmin,psimin,npd,icon,                                 &
                niter,printpot,pchem,irespar,                           &
                pdenpar,pener,ppot,iron,ironx,                          &
                ximp,yimp,zimp,rimpur,mimpur,                           &
!                 gwf,leepot,nq,tol,rmin,vxpot,rinfi,                     &
                norder,nv_iter,deltat,lrkpc,lrk,xlamda,                 &
                xlamdax,xlamday,xlamdaz,deltatps,                       &
                xlamdax_x,xlamdax_y,xlamdax_z,iter0,time0,Zsurf,        &
                pcurr,icurr,lsurf,lsurf3D,Lambdah,tzmean,tzsurf,        &
                vximp,vyimp,vzimp,ndmax,ninvar,                         &
                selec_gs,selec_pi,selec_sigma,selec_delta,              &
                r_cutoff_gs,r_cutoff_pi,r_cutoff_sigma,r_cutoff_delta,  &
                umax_gs,umax_pi,umax_sigma,umax_delta,Ev0,              &
                Lexcite_state,instate,Als_P,Als_D,                      &
                Lexcite_state_fix,Lexciplex_axis,                       &
                Quita_C4_Ba_plus_gs_fix_C4,                             &
                Quita_C4_Ba_plus_pi_fix_C4,                             &
                Quita_C4_Ba_plus_sigma_fix_C4,Lstate,Ldiag_jz,Ljz,      &
                Laverage_P_value,Lexcite_state_external,                &
                Exciplex, Lexciplex_state_fix,r_exc, Lfix_lambda,       &
                Lmixture,a_pi_3o2,a_sig_1o2,L_v_aditional,              &
                Lfilter_exciplex_force, z_exciplex_exclusion,           &
                z_exciplex_position
 
!................................ Start main Program ..............................
call timer(t0)

!.............................................
!... Inicializate some numerical constants ...
!.............................................

pi     = 4.0d0*datan(1.0d0) ! Initialization of pi
twopi  = 2.0d0*pi
fourpi = 4.0d0*pi
piq    = pi*pi


!...............................
!... Read  master-input file ...
!...............................
!Lambdah = 0.d0
!tzmean   = 0.d0
!tzsurf   = 1.d0

read(5,nml=input,end=999)
open(10,file="DFT4He3d.namelist.read")
write(10,nml=input)
call flush(10)

nn(1)  = nx ; nn(2)  = ny ; nn(3)  = nz;                ! Initialize things for PDERG
mmx(1) = nx ; mmx(2) = ny ; mmx(3) = nx ; mmx(4) = ny   ! (NO NOT MOVE THAT!!!!!!!)
Also2 = 0.5d0*Als

!.............................................................
!.. Check if the size of the grid is among the valid values ..
!.............................................................

if(tstgrid(nx).ne.0 ) stop 'SEVERE ERROR: NX is not correct'
if(tstgrid(ny).ne.0 ) stop 'SEVERE ERROR: NY is not correct'
if(tstgrid(nz).ne.0 ) stop 'SEVERE ERROR: NZ is not correct'

!...................................................
!.. Controls Paflov parameters (read and storage ...
!...................................................

      Call Init_deriv_p(npd,ndmax,nthread)

close(5)
close(10)

!.........................................
!.. Some consistency check on Delta t  ...
!.........................................
If(deltat.eq.0.d0)then
 if(deltatps.eq.0.d0)then    ! deltat ==  0, deltatps ==  0
  print*,'You must specify either Deltat (in kelvin) or Deltatps (in picosecond)'
  STOP
 else                        ! deltat ==  0, deltatps =/= 0
  deltat = deltatps/Ktops
 endif
Else
 if(deltatps.eq.0.d0)then    ! deltat =/= 0, deltatps ==  0
  deltatps = deltat*Ktops
 else                        ! deltat =/= 0, deltatps =/= 0
  if(deltatps.ne.(deltat*Ktops))then
   print *,'Inconsistent deltat - deltatps'
   STOP
  endif
 endif
Endif
write(*,*)'Time step is ',deltat,' kelvins or ',deltatps,' picoseconds.'


!................................................
!.. Some consistency check on input variables ...
!................................................

nthread=abs(nthread)

! if(mode.eq.2) then
! !   if(.not.limp) then
!      write(6,*) ' '
!      write(6,*) ' Inconsistency input error.'
!      write(6,*) ' '
!      write(6,*) ' If mode=2. MUST be limp=.true.'
!      write(6,*) ' '
!      write(6,*) ' ACTION: ABORT EXECUTION'
!      write(6,*) ' '
!      STOP
! !   end if
! end if

mAg_u = mimpur*mp_u

hx    = 2.0d0*abs(xmax)/(nx)  ! Step in x-grid
hy    = 2.0d0*abs(ymax)/(ny)  ! Step in y-grid
hz    = 2.0d0*abs(zmax)/(nz)  ! Step in z-grid

dxyz  = hx*hy*hz              ! Element of volum in real space
nxyz  = nx*ny*nz              ! Total number of points
dVomAg = dxyz/mAg_u

hpx   = 1.0d0/(nx*hx)         ! Step in px-grid
hpy   = 1.0d0/(ny*hy)         ! Step in py-grid
hpz   = 1.0d0/(nz*hz)         ! Step in pz-grid

pmaxx = 1.0d0/(2.0d0*hx)      ! Maximum 'frequency' in X-grid
pmaxy = 1.0d0/(2.0d0*hy)      ! Maximum 'frequency' in Y-grid
pmaxz = 1.0d0/(2.0d0*hz)      ! Maximum 'frequency' in Z-grid


!...............................
!.. Dimensionate main ARRAYS ...
!...............................

call dimen()

!.....................................
!.. Initial value of rimp and vimp ...
!.....................................
rimp = (/ximp,yimp,zimp/)
vimp = (/vximp,vyimp,vzimp/)

!If(L_v_aditional)Then
!  Write(6,'(" We add this velocity (Ang./ps) to the impurity...:",1p, 3E15.6)')vimp
!Endif

xlambda0 = -vimp*mimpur/(4.d0*n4)

vimp = vimp*Ktops ! Because it is given in A/ps, we have to transform to A*K

vimp0 = vimp

!................................
!... Build grid in real space ...
!................................

do ix=1,nx  !.................... Grid X
 x(ix) = -xmax+hx*(ix-1)
end do
do iy=1,ny  !.................... Grid Y
 y(iy) = -ymax+hy*(iy-1)
end do
do iz=1,nz  !.................... Grid  Z
 z(iz) = -zmax+hz*(iz-1)
end do

!....................................
!... Build grid in momentum space ...
!....................................

!.... Build p-grid. In order to use FFTW the grid must
!     start from frequency zero to the maximum and then continue from
!     (-maximum) to zero (negative).

!............................................ grid Px
do ipx=1,nx/2+1
   px(ipx) =        hpx*(ipx-1)
end do
do ipx=nx/2+2,nx
   px(ipx) = -pmaxx+hpx*(ipx-(nx/2+1))
end do

!............................................ grid Py
do ipy=1,ny/2+1
   py(ipy) =        hpy*(ipy-1)
end do
do ipy=ny/2+2,ny
   py(ipy) = -pmaxy+hpy*(ipy-(ny/2+1))
end do

!............................................ grid Pz
do ipz=1,nz/2+1
   pz(ipz) =        hpz*(ipz-1)
end do
do ipz=nz/2+2,nz
   pz(ipz) = -pmaxz+hpz*(ipz-(nz/2+1))
end do

!............................................ Compule modulus of p
do ipz=1,nz
  pz2=pz(ipz)**2
  do ipy=1,ny
    py2=py(ipy)**2
    do ipx=1,nx/2+1
      px2               = px(ipx)**2
      pmod(ipx,ipy,ipz) = sqrt(px2+py2+pz2)
    end do
  end do
end do

pmod1=maxval(pmod)
pmod2=sqrt(pmaxx**2+pmaxy**2+pmaxz**2)


!write(6,*) '    Initialize Linear Interpolation for V_Pi and V_Sig'
call flush(8)
call potenimpini() ! interpolation + first call to updatepoten

!................................
!... read density or build-it ...
!................................

call readenc(n4,densat4,filedenin,fileimpin,mode,instate,rimpur,r_clust)

If(L_v_aditional)Then
  vimp = vimp + vimp0
  !write(6,'(" We impose linear momentum conservation..:",1p,3E15.6)')xlambda0
  xlambda0 = xlambda0*Ktops/(2.d0*h2o2m4)
   do iz=1,nz
     do iy=1,ny
       do ix=1,nx
         aux=x(ix)*xlambda0(1)       &
            +y(iy)*xlambda0(2)       &
            +z(iz)*xlambda0(3)
         psi(ix,iy,iz) = psi(ix,iy,iz) &
                       * cmplx(cos(aux),sin(aux))
       end do
     end do
   end do
Endif        

!...............................................................test
 do iz=1,nz ; do iy=1,ny ; do ix=1,nx 
  if(.not.(den(ix,iy,iz).gt.0))print*,ix,iy,iz,den(ix,iy,iz)
 end do ; enddo ; enddo
!...............................................................test


! if(irespar.ne.0) then
!    if(limp) then
!       call respar(x,y,z,nx,ny,nz,2,'hedenini','impurdenini',den,denx)
!    else
!      call respar(x,y,z,nx,ny,nz,1,'hedenini','hedenini',den,den)
!    end if
! end if

!....................................
!.. Print-out initial quantities  ...
!....................................

!open(6,file=fileout)

write(6,6010) title

select case(mode)
   case(0) !................................... Start a dynamic form static calcultaion
!       if(limp) then
!          write(6,6111) filedenin,filedenout,fileimpin,fileimpout
!       else
         !write(6,6011) filedenin,filedenout
!       end if
   case(1) !................................... Continue a dynamic calculation from a prevous one (non-excited)
      !write(6,6013) filedenout,fileimpout
   case(2) !................................... Continue a dynamic calculation from a prevous one
      !write(6,6013) filedenout,fileimpout
   case(3) !................................... Continue a dynamic calculation from a prevous one
      !write(6,6013) filedenout,fileimpout
      !write(6,*) ' '
   case default !.............................. Start still not programed.
      !write(6,*) ' The variable mode has no acceptable value.'
      !write(6,*) ' '
      stop
end select

!...............................................................
!................................ Write the input parameters ...
!...............................................................

!write(6,6018) nthread,niter
!if(mode.ne.0) then
!  write(6,6020) n4,r_clust
!else
!  write(6,6025) n4
!end if
!write(6,6030) nx,ny,nz,hx, hy, hz, x(1) ,y(1) ,z(1) ,x(nx),y(ny),z(nz)
!write(6,6035)          hpx,hpz,hpz,px(1),py(1),pz(1),pmaxx,pmaxy,pmaxz,&
!                       pmod1,pmod2
!write(6,6037) cp4,cpp4,den4c,alphas,l,den0s,h2o2m4


!...............
!.. Impurity ...
!...............


!...................................................................
!... Compute the FFT of Lennard-Jones                            ...
!... Prepara \alpha_s term in case of Orsay-Trento Interaction.  ...
!...................................................................

select case(core4)
   case('OP ')
     h4=h4op
     !write(6,*) '    Use Orsay-Paris-Barcelona Interaction.'
     !write(6,6040) core4,h4,eps4,sigma4,b4
   case('OT ')
     h4=h4ot
     !write(6,*) '    Use Orsay-Trento Interaction. (ONLY CORE)'
     !write(6,*) '    Do not calculate Alpha_s in the field neither energy.'
     !write(6,6040) core4,h4,eps4,sigma4,b4
   case('OTC')
     h4=h4ot
     !write(6,*) '    Use Orsay-Trento Interaction.'
     !write(6,*) '    Full Orsay-Trento calculation. (Field and Energy)'
     !write(6,6040) core4,h4,eps4,sigma4,b4
!     allocate( denalf(nx  ,ny,nz))                                                            
     allocate(  falfs(nx  ,ny,nz))
     allocate(kalfs(nx/2+1,ny,nz))
!     allocate(intxalf(nx  ,ny,nz))
!     allocate(intyalf(nx  ,ny,nz))
!     allocate(intzalf(nx  ,ny,nz))
!     allocate(ualphas(nx  ,ny,nz))
   case default
     print *,' ***************** WARNING ************************'
     print *,' '
     print *,' I do not know how to work with the core: ',core4
     print *,' '
     print *,' **************************************************'
     print *,' '
     STOP ' ... The program stops due to a severe error.'
end select

     allocate( denalf(nx  ,ny,nz))                                                            
     allocate(intxalf(nx  ,ny,nz))
     allocate(intyalf(nx  ,ny,nz))
     allocate(intzalf(nx  ,ny,nz))
     allocate(ualphas(nx  ,ny,nz))


!...............................
!... Prepare plans for FFTWs ...
!...............................
write(6,*) '    Initialize Plans for FFT.'
!call fftini(nx,ny,nz)
call fftini(nx,ny,nz,nthread)

!...........................................................
!... Form  factor for Lennard-Jones and for the impurity ...
!...........................................................

write(6,*) '    Compute the FFT of the kernel of Lennard-Jones integrals.'

call fforma(core4,b4,eps4,sigma4,h4,nx,ny,nz,pmod,fvlj4)


!........................................
!.. Initialize coarse-graining kernel ...
!........................................

write(6,*) '    Initialize Coarse-graining kernel.'
call initcg(h4,wcgk)


if(core4.eq.'OTC') then
   forall(ix=1:nx/2+1,iy=1:ny,iz=1:nz)
      kalfs(ix,iy,iz) = exp(-(pi*l*pmod(ix,iy,iz))**2)
   end forall
end if

!
!    We compute the perturbate wave function
!
! Initial velocity : Altough xlamdax_i is in K units, xlamda
! is introduced in Angstrom/picosecond for the sake of lazyness.
auxn4 = sum(den)*dxyz

if(xlamda.ne.0.d0)then
write(*,*)'xlambda not equal zero'
xlamdaz = xlamda*Ktops/(2.d0*h2o2m4)
! xlamdaz = xlamda*0.630247988
endif

If(mode.eq.0)then
   do iz=1,nz
     do iy=1,ny
       do ix=1,nx
         aux=x(ix)*xlamdax       &
            +y(iy)*xlamday       &
            +z(iz)*xlamdaz
         psi(ix,iy,iz) = sqrt(den(ix,iy,iz)) &
                       * cmplx(cos(aux),sin(aux))
       end do
     end do
   end do
Endif
!.................................
!.. First call to total energy ...
!.................................

Lfirst = .true.
Lprint_invar=.true.

den=Conjg(psi)*psi
call potenimp(rimp,invar)
call poten()              ! First Potential  (for Lagrange Equation)
call energy()             ! Calculate energies

Lfirst = .false.
Lprint_invar=.false.

If(Ev0.Eq.0.0d0)Ev0 = eHeX + eso

!Write(6,'("Energia que fem servir per restar una fase global",/,                 &
!          "per la evolució dels estats interns. Ev0........:",1p,E15.6)')Ev0
If(Lfilter_exciplex_force)Then
!  Write(6,'("From potenimpini: We will use ",(A)," for He-He interaction")')Trim(selec_gs)

  ipx = 1.5 +(rimp(1)+xmax)/hx
  ipy = 1.5 +(rimp(2)+ymax)/hy
  ipz = 1.5 +(rimp(3)+zmax)/hz
!
!  Open(Unit = 1, File='uimp-x.0.dat')
!  Write(1,'("#x,uimp( x, rimp(2), rimp(3) )")')
!  Do ix=1, nx
!    Write(1,'(1p,2E15.6)')x(ix),Uimp(ix,ipy,ipz)
!  EndDo
!  Close(Unit=1)
!
!  Open(Unit = 1, File='uimp-y.0.dat')
!  Write(1,'("#y,uimp( rimp(1), y, rimp(3) )")')
!  Do iy=1, ny
!    Write(1,'(1p,2E15.6)')y(iy),Uimp(ipx,iy,ipz)
!  EndDo
!  Close(Unit=1)
!
!  Open(Unit = 1, File='uHe_He-z.0.dat')
!  Write(1,'("#z, uHe_He( rimp(1), rimp(2), z )")')
!  Do iz=1, nz
!    Write(1,'(1p,2E15.6)')z(iz),UHe_He(ipx,ipy,iz)
!  EndDo
!  Close(Unit=1)
EndiF

call forceimp(rimp,aimp)
aimp(:) = aimp(:)/mAg_u
call potinvar(rimp,invar,Hinvar)


!Write(6,'("Variables internas....:",/,1p,(2E15.6))')invar

!    Write(6,'("From Main(1)...")')
!    write(6,'("rimp(1,2,3)..",1p,3e15.6)')rimp
!    write(6,'("vimp(1,2,3)..",1p,3e15.6)')vimp
!    write(6,'("Invar.....",/,1p,(2E15.8))')(invar(i),i=1,ninvar)
!    write(6,'("Hinvar....",/,1p,(2E15.8))')(Hinvar(i),i=1,ninvar)


! TEST: Treu el valor de uimp
!call respar(x,y,z,nx,ny,nz,1,'uimp','den',uimp,den)

!
! Escrivim el valor de Uimp pels aixos que passen per l'impureça
!

ipx = 1.5 +(rimp(1)+xmax)/hx
ipy = 1.5 +(rimp(2)+ymax)/hy
ipz = 1.5 +(rimp(3)+zmax)/hz

!ipx = nx/2 + 1
!ipy = ny/2 + 1
!ipz = nz/2 + 1

!Open(Unit = 1, File='uimp-x.0.dat')
!Write(1,'("#x,uimp( x, rimp(2), rimp(3) )")')
!Do ix=1, nx
!  Write(1,'(1p,2E15.6)')x(ix),Uimp(ix,ipy,ipz)
!EndDo
!Close(Unit=1)
!
!Open(Unit = 1, File='uimp-y.0.dat')
!Write(1,'("#y,uimp( rimp(1), y, rimp(3) )")')
!Do iy=1, ny
!  Write(1,'(1p,2E15.6)')y(iy),Uimp(ipx,iy,ipz)
!EndDo
!Close(Unit=1)
!
!Open(Unit = 1, File='uimp-z.0.dat')
!Write(1,'("#z, uimp( rimp(1), rimp(2), z )")')
!Do iz=1, nz
!  Write(1,'(1p,2E15.6)')z(iz),Uimp(ipx,ipy,iz)
!EndDo
!Close(Unit=1)

!write(6,'("Number of He4 atoms",1P,E15.6)')auxn4

!write(6,6050) auxn4,etot4,etot4/auxn4,ekin4,elj4,ealphas,esolid,ecor4
write(6,6060) eimpu,ekinx,eHeX,eso,etot
write(6,6065) rimp(1),rimp(2),rimp(3)
do iz=1,100
	rimp(3)=rimp(3)+hz
	call potenimp(rimp,invar)
	call energy()
!	write(6,6050) auxn4,etot4,etot4/auxn4,ekin4,elj4,ealphas,esolid,ecor4
	write(6,6060) eimpu,ekinx,eHeX,eso,etot
	write(6,6065) rimp(1),rimp(2),rimp(3)
end do



!write(6,6065) rimp(1),rimp(2),rimp(3)

!eold = etot4
!call flush(6)
!
!!!-------------------------------------------------------------------------------
!!!---                            Iterative procedure                           --
!!!-------------------------------------------------------------------------------
!!
!!! TIME CONSTANT !
!!! This time it's a cylinder.
!!do iz=1,nz
!! do iy=1,ny
!!  do ix=1,nx
!!!    rt = dsqrt(x(ix)*x(ix)+y(iy)*y(iy)+z(iz)*z(iz))
!!   rt = dsqrt(x(ix)*x(ix)+y(iy)*y(iy))
!!!    timec(ix,iy,iz)=cmplx(Lambdah*(1.d0+tanh((rt-tmean)/tsurf)),1.d0)
!!   !timec(ix,iy,iz)=cmplx(Lambdah*(1.d0+tanh((z(iz)-tzmean)/tzsurf)),1.d0)
!!   timec(ix,iy,iz)=cmplx(Lambdah*(1.d0+tanh((abs(z(iz))-tzmean)/tzsurf)),1.d0)
!!  enddo
!! enddo
!!enddo
!!
!!
!!!plot it
!!! call respar(x,y,z,nx,ny,nz,1,'timec','den',timec,den)
!!open(unit=32,file='timec.dat')
!!do iz=1,nz
!! do ix=1,nx
!!  write(32,*)x(ix),z(iz),real(timec(ix,ny/2+1,iz))
!! enddo
!! write(32,*)''
!!enddo
!!close(32)
!!
!!open(81,file=filerimp)
!!open(82,file=filevimp)
!!open(83,file=fileaimp)
!!open(84,file=filelambda)
!!
!!
!!! If(limp)Then
!!!   write(81,1381)iter0*deltat,aux1,aux2,aux3,aux4,aux5,aux6
!!! Else
!!
!!! Endif
!!call flush(81)
!!
!!lfilepv  = .false.
!!
!!! write(6,7000)
!!call timer(t5)
!!
!!uext = 0.0d0
!!
!!!       print*,rimp
!!!       print*,vimp
!!qr(:) = 0.d0
!!qv(:) = 0.d0
!!! q(:,:,:)=0.d0
!!rimpold(:,:) = 0.d0
!!vimpold(:,:) = 0.d0
!!aimpold(:,:) = 0.d0
!!
!!!
!!! Eigenvectors of SO matrix (P states)
!!!
!!
!!!  ev1pSO = (/(0.d0,0.d0) , (0.d0,0.d0) , uim        , (0.d0,0.d0), (0.d0,0.d0), (1.d0,0.d0)/)
!!!  ev2pSO = (/(0.d0,0.d0) , (0.d0,0.d0) , (0.d0,0.d0), -uim       , (-1.d0,0.d0), (0.d0,0.d0)/)
!!!  ev3pSO = (/-uim        , (0.d0,0.d0) , (0.d0,0.d0), (0.d0,0.d0), (0.d0,0.d0), uim        /)
!!!  ev4pSO = (/(0.d0,0.d0) , uim	      , (0.d0,0.d0), (0.d0,0.d0), uim        , (0.d0,0.d0)/)
!!!  ev5pSO = (/(0.d0,0.d0) , (-1.d0,0.d0), (0.d0,0.d0), -uim       , (0.d0,0.d0), (0.d0,0.d0)/)
!!!  ev6pSO = (/(1.d0,0.d0) , (0.d0,0.d0) , -uim       , (0.d0,0.d0), (0.d0,0.d0), (0.d0,0.d0)/)
!!! 
!!!  do i = 1,ninvar
!!!     evpSO (i,1) = ev1pSO(i)
!!!     evpSO (i,2) = ev2pSO(i)
!!!     evpSO (i,3) = ev3pSO(i)
!!!     evpSO (i,4) = ev4pSO(i)
!!!     evpSO (i,5) = ev5pSO(i)
!!!     evpSO (i,6) = ev6pSO(i)
!!!  end do
!!
!! evpSO = (0.d0,0.d0)
!! 
!! evpSO(1,3)=-uim
!! evpSO(2,4)=uim
!! evpSO(2,5)=(-1.d0,0.d0)
!! evpSO(1,6)=(1.d0,0.d0)
!! evpSO(4,5)=-uim
!! evpSO(3,6)=-uim
!! 
!! evpSO(3,1)=uim
!! evpSO(4,2)=-uim
!! evpSO(5,2)=(-1.d0,0.d0)
!! evpSO(6,1)=(1.d0,0.d0)
!! evpSO(5,4)=uim
!! evpSO(6,3)=uim
!!
!!if(time0 .LT. 0.d0) time0 = iter0*deltatps
!!write(81,'("# Tiempo(ps), x(AA), y(AA), z(AA)")')
!!write(81,'(1x,1p,E15.6,3E18.10)')time0, rimp(1), rimp(2), rimp(3)
!!write(82,'("# Tiempo(ps), Vx(AA/ps), Vy(AA/ps), Vz(AA/ps)")')
!!write(82,'(1x,1p,E15.6,3E18.10)')time0, vimp(1)*pstoK, vimp(2)*pstoK, vimp(3)*pstoK
!!write(83,'("# Tiempo(ps), Ax(AA/ps**2), Ay(AA/ps**2), Az(AA/ps**2)")')
!!write(83,'(1x,1p,E15.6,3E18.10)')time0, aimp(1)*pstoK*pstoK, aimp(2)*pstoK*pstoK, aimp(3)*pstoK*pstoK
!!write(84,'(I10,1p,20E18.10)')iter0, invar
!!
!!open(11,file='projections.dat')
!!
!!      prodLambda = 0.d0
!!      do i=1,ninvar
!!	prodLambda(i) = DOT_PRODUCT(conjg(invar),evpSO(i,:))
!!      end do 
!!	aux1 = DOT_PRODUCT(invar,prodLambda)*0.5
!!  
!!      auxn4 = sum(den)*dxyz
!!   
!!      write(11,'("# Iteracion, <L·S>, |<Lamb0|Lambt>|**2")')
!!      write(11,'(I10,1p,2E16.8)') iter0,aux1,ABS(DOT_PRODUCT(invar0,invar))**2
!!
!!iter0=iter0+1
!!
!!
!!do iter=iter0,niter  ! <--------- Iterative procedure starts here.
!!
!!!    Write(6,'("From Main(2)...")')
!!!    write(6,'("rimp(1,2,3)..",1p,3e15.6)')rimp
!!!    write(6,'("vimp(1,2,3)..",1p,3e15.6)')vimp
!!!    write(6,'("Invar.....",/,1p,(2E15.8))')(invar(i),i=1,ninvar)
!!!    write(6,'("Hinvar....",/,1p,(2E15.8))')(Hinvar(i),i=1,ninvar)
!!
!!    if(iter.le.(2+iter0).Or.lrk)then
!!      call steprk(deltat)
!!
!!    else
!!      call steppc(deltat,ErrHe,Errrimp,Errvimp,Erriv)
!!
!!      write(6,'(" Steppc: ErrHe, Errrimp, Errvimp, Erriv...",1p,4E15.6)')ErrHe, Errrimp, Errvimp, Erriv
!!
!!    endif
!!
!!
!!    call potenimp(rimp,invar)
!!    call poten()
!!    call forceimp(rimp,aimp)
!!    aimp(:) = aimp(:)/mAg_u
!!    call potinvar(rimp,invar,Hinvar)
!!
!!    aux1 = time0 + (iter-iter0+1)*deltatps
!!    temps = aux1
!!
!!write(81,'(1x,1p,E15.6,3E18.10)')aux1, rimp(1), rimp(2), rimp(3)
!!write(82,'(1x,1p,E15.6,3E18.10)')aux1, vimp(1)*pstoK, vimp(2)*pstoK, vimp(3)*pstoK
!!write(83,'(1x,1p,E15.6,3E18.10)')aux1, aimp(1)*pstoK*pstoK, aimp(2)*pstoK*pstoK, aimp(3)*pstoK*pstoK
!!
!!call flush(81)
!!call flush(82)
!!call flush(83)
!!!
!!! Proyección < lambda | lambda_SO >
!!!
!!
!!   if(mod(iter,pener).eq.0) then          ! Compute New energy and max of density
!!   
!!      prodLambda = 0.d0
!!      do i=1,ninvar
!!	prodLambda(i) = DOT_PRODUCT(invar,evpSO(:,i))
!!      end do 
!!	aux1 = DOT_PRODUCT(conjg(invar),prodLambda)*0.5
!!  
!!      auxn4 = sum(den)*dxyz
!!      
!!      write(11,'(I10,1p,2E16.8)') iter,aux1,ABS(DOT_PRODUCT(invar0,invar))**2
!!   
!!      call energy()
!!      print*," "
!!      print*,"Iteration:", iter
!!      print*," "
!!      write(*,*)'... Internal Variables ...'
!!!       do ix=1,ninvar
!!!        write(*,'("invar ",I1,1X,2E13.5)')ix,real(invar(ix)),aimag(invar(ix))
!!!       enddo
!!      write(*,*)'..........................'
!!      write(*,*)'Mod of internal vars:',sum(abs(invar)**2)
!!      
!!      write(84,'(I10,1p,20E18.10)')iter, invar
!!      call flush(84)
!!      write(6,'("Re, Img, Mod: 1,2...",/,(0p,I4,1p,3E16.8))') (ix,real(invar(ix)), aimag(invar(ix)), abs(invar(ix)**2),ix=1,ninvar)
!!      
!!      write(6,7010) auxn4,etot4,(etot4-eold),etot4/auxn4,ekin4,elj4,ealphas,esolid,ecor4
!!!      call varmu(n4,mu4,errmu4)           ! Error in mu4
!!      write(6,7017) mu4,errmu4
!!
!!           write(6,7015) eimpu,ekinx,eHeX,eso,etot
!!
!!      eold = etot4
!!
!!        call r_cm(den,n4,xcm4,ycm4,zcm4)    ! Center of mass of 4He Drop
!!        write(6,7100) xcm4,ycm4,zcm4
!!        write(*,*)'Impurity position:',rimp
!!        write(*,*)'Impurity velocity:',vimp
!!
!!   end if
!!
!!!..............................................................................
!!
!!   if(mod(iter,pcurr).eq.0) then        ! Save wavefunction for current
!!
!!     ncurr = (iter-iter0+1)/pcurr + icurr
!!      select case (ncurr)
!!       case(1:9)
!!      write(chariter,8011)ncurr
!!       case(10:99)
!!      write(chariter,8012)ncurr
!!       case(100:999)
!!      write(chariter,8013)ncurr
!!      end select
!!      namefile='density.'//chariter//'.dat'
!!!        call printoutc(3,namefile,namefile1,psi,elem,nx,ny,nz,hx,hy,hz,limp, &
!!!                     xmax,ymax,zmax,ximp,yimp,zimp,psi,         &
!!!                     deltat,iter)
!!       call printoutc(temps,3,namefile,namefile1,psi,nx,ny,nz,hx,hy,hz, &
!!                    xmax,ymax,zmax,rimp,vimp,psi,         &
!!                    deltatps,iter,invar,ninvar)
!!
!!   endif
!!
!!!    if(mod(iter,pdenpar).eq.0) then        ! Save partial densities
!!!      nsfaux = nsfaux+1
!!!      nsfaux2= nsfaux2+1 ! The same as nsfaux but not in mod(nsfiles)
!!!      if(nsfaux.gt.nsfiles) nsfaux=mod(nsfaux,nsfiles)
!!!      select case(nsfaux)
!!!          case(1:9)
!!!            write(namefile, 8010) nsfaux
!!!            write(namefile1,8015) nsfaux
!!!          case(10:99)
!!!            write(namefile, 8020) nsfaux
!!!            write(namefile1,8025) nsfaux
!!!          case(100:999)
!!!            write(namefile, 8030) nsfaux
!!!            write(namefile1,8035) nsfaux
!!!      end select
!!! 
!!! ! Contour plots
!!!       select case (nsfaux2)
!!!        case(1:9)
!!!       write(chariter,8011)nsfaux2
!!!        case(10:99)
!!!       write(chariter,8012)nsfaux2
!!!        case(100:999)
!!!       write(chariter,8013)nsfaux2
!!!       end select
!!! 
!!!      curvfile='curvasnivel.y=0.'//chariter//'.dat'
!!!      open(74,file=curvfile)
!!!      do ix=1,nx
!!!        do iz=1,nz
!!!          write(74,7109) x(ix),z(iz),den(ix,ny/2+1,iz)!,dene(ix,ny/2+1,iz)
!!!        end do
!!!        write(74,7119)
!!!      Enddo
!!!      close(74)
!!! 7109 format(T5,0P,2F15.5,3x,1P,1E15.5)
!!! 7119 format()
!!! 
!!!      if(limp) then
!!!        if(irespar.ne.0) call respar(x,y,z,nx,ny,nz,2,'folding','foldingx',potx4,upotx)
!!!        call printoutc(3,namefile,namefile1,psi,elem,nx,ny,nz,hx,hy,hz,limp, &
!!!                     xmax,ymax,zmax,ximp,yimp,zimp,psix,         &
!!!                     deltat,iter)
!!!      else
!!!        call respar(x,y,z,nx,ny,nz,1,'den','den',den,den)
!!!        call printoutc(3,namefile,namefile1,psi,elem,nx,ny,nz,hx,hy,hz,limp, &
!!!                     xmax,ymax,zmax,ximp,yimp,zimp,psi,         &
!!!                     deltat,iter)
!!!      end if
!!!    end if
!!
!!!..............................................................................
!!
!!   if(lfilepv) then
!!     select case(iter)
!!         case(1:9)
!!           write(namefile ,5010) iter
!!           write(namefile1,5015) iter
!!         case(10:99)
!!           write(namefile ,5020) iter
!!           write(namefile1,5025) iter
!!         case(100:999)
!!           write(namefile ,5030) iter
!!           write(namefile1,5035) iter
!!         case(1000:9999)
!!           write(namefile ,5040) iter
!!           write(namefile1,5045) iter
!!         case(10000:99999)
!!           write(namefile ,5050) iter
!!           write(namefile1,5055) iter
!!         case(100000:999999)
!!           write(namefile ,5060) iter
!!           write(namefile1,5065) iter
!!     end select
!!
!!
!!!      if(limp) then
!!!         call printoutc(temps,2,namefile,namefile1,psi,elem,nx,ny,nz,hx,hy,hz,limp, &
!!!                       xmax,ymax,zmax,ximp,yimp,zimp,psix,                   &
!!!                       deltat,iter)
!!!      else
!!        call printoutc(temps,2,namefile,namefile1,psi,nx,ny,nz,hx,hy,hz, &
!!                      xmax,ymax,zmax,rimp,vimp,psi,                   &
!!                      deltatps,iter,invar,ninvar)
!!!      end if
!!
!!
!!   end if
!!!..............................................................................
!!
!!   call timer(t6)                         ! Compute use time
!!
!!!    if(mod(iter,pchem).eq.0) then          ! Save partial densities
!!! !      if(limp) then
!!! !        write(6,7001)
!!! !        write(6,7035) iter,mu4,mu4err,epsx,epsxerr,t6,(t6-t5)
!!! !      else
!!!        write(6,7000)
!!!        write(6,7030) iter,mu4,mu4err,t6,(t6-t5)
!!! !      end if
!!!    end if
!!   t5=t6
!!end do     ! <---- END OF ITERATIVE PROCEDURE
!!
!!close(11)
!
!
!! if(limp) then
!!    call printoutc(1,filedenout,fileimpout,psi,elem,nx,ny,nz,hx,hy,hz,limp, &
!!                  xmax,ymax,zmax,ximp,yimp,zimp,psix,         &
!!                  deltat,iter)
!! else
!!call printoutc(temps,1,filedenout,fileimpout,psi,nx,ny,nz,hx,hy,hz, &
!!				xmax,ymax,zmax,rimp,vimp,psi, &
!!				deltatps,iter,invar,ninvar)
!! end if
!
!call timer(t4)
!print *,' Total  ',t4-t0

stop
999 stop 'DFT3He3d. Error in input master file. Too short'

!...............
!... Formats ...
!...............

3100 format(3x,0P,f9.4,2x,1P,E13.5)

3156 format(10E13.5)

6010 format(//,&
T10,'   ######  ####### ####### #       #     #          #####          ',/,  &
T10,'   #     # #          #    #    #  #     #  ###### #     #  #####  ',/,  &
T10,'   #     # #          #    #    #  #     #  #            #  #    # ',/,  &
T10,'   #     # #####      #    #    #  #######  #####   #####   #    # ',/,  &
T10,'   #     # #          #    ####### #     #  #            #  #    # ',/,  &
T10,'   #     # #          #         #  #     #  #      #     #  #    # ',/,  &
T10,'   ######  #          #         #  #     #  ######  #####   #####  ',//, &
T6,'Title of the run: ',A)

6011 format(//,T6,'CONTINUE a calculation:',//,&
               T6,'Input  densitity file: ',A,/,&
               T6,'Output densitity file: ',A)

6111 format(//,T6,'CONTINUE a calculation:',//,&
               T6,'Input file with helium densitity    : ',A,/,&
               T6,'Input file with impurity wave func. : ',A,/,&
               T6,'Output file with helium densitity   : ',A,/,&
               T6,'Output file with impurity wave func.: ',A,/,' ')

6012 format(//,T6,'Start a new calculation:',//,&
               T6,'Output densitity file: ',A)
6013 format(//,T6,'Start a new calculation with an impurity:',//,&
               T6,'Output file for Helium density: ',A,/,        &
               T6,'Output file for the impurity wave function: ',A)
6018 format(//,T6,'Number of threads:    ',I6,/,&
               T6,'Number of iterations: ',i6)
6020 format(//,T6,'Number of particles:    ',0P,I10,/,&
               T6,'Radius of the cluster : ',F10.3,' A')
6025 format(//,T6,'Number of particles:    ',0P,I10,/,' ')
6030 format(//,T6,'+-------------------+----------------------------------------+',/,&
               T6,'| REAL GRID         |     X-grid       Y-grid       Z-grid   |',/,&
               T6,'+-------------------+----------------------------------------+',/,&
               T6,'| Number of points  |',0P,T32,I4,T45,I4,T58,I4,T66,' |',/,&
               T6,'| Step              |',1P,T27,E12.5,T40,E12.5,T53,E12.5,T66,' |',/,&
               T6,'| Min value         |',1P,T27,E12.5,T40,E12.5,T53,E12.5,T66,' |',/,&
               T6,'| Max value         |',1P,T27,E12.5,T40,E12.5,T53,E12.5,T66,' |',/,&
               T6,'+-------------------+----------------------------------------+')
6035 format(//,T6,'+-------------------+----------------------------------------+',/,&
               T6,'| MOMEMTUM GRID     |    Px-grid      Py-grid      Pz-grid   |',/,&
               T6,'+-------------------+----------------------------------------+',/,&
               T6,'| Step              |',1P,T27,E12.5,T40,E12.5,T53,E12.5,T66,' |',/,&
               T6,'| Min value         |',1P,T27,E12.5,T40,E12.5,T53,E12.5,T66,' |',/,&
               T6,'| Max value         |',1P,T27,E12.5,T40,E12.5,T53,E12.5,T66,' |',/,&
               T6,'+-------------------+----------------------------------------+' ,//,&
               T6,'Maximum modulus of p ',0P,F11.6,3x,F11.6,/,' ')
6037 format(//,T6,'Parameters for the funcional:',//,          &
               T8,'cp4 ...... ',1P,E13.5 ,' K \AA**6'      ,/, &
               T8,'cpp4 ..... ',1P,E13.5 ,' K \AA**9'      ,/, &
               T8,'den4c .... ',0P,F13.3 ,' \AA**{-3}'     ,/, &
               T8,'Alphas ... ',0P,F13.3 ,' K ^-1 \AA**3'  ,/, &
               T8,'L ........ ',0P,F13.2 ,' \AA'           ,/, &
               T8,'den0s .... ',0P,F13.2 ,' \AA**-3'       ,/, &
               T8,'h2o2m4 ... ',0P,F14.11,' hbar**2 / (2 m_4)' )
6138 format(' ',T8,'h2o2mx ... ',0P,F14.11,' hbar**2 / (2 m_x)' )
6038 format(//,T6,'Change of Paflov parameter allowed: ',//, &
     T18,'From     to      iter     Factor',/, &
     T19,'------  ------  ------  -----------')

6039 format(1x,0P,T17,i6,T25,I6,T33,i6,t42,f11.7)

6150 format(//,T6,'Pavlov parameter fixed for all the run to: ',F8.4)

6040 format( /,T6,'Lennard-Jones parameters:',//,&
               T10,'Core    ',A3,/,&
               T10,'h     ',F11.7,' A',/,&
               T10,'eps   ',F11.7,' K',/,&
               T10,'sigma ',F11.7,' A',/,&
               T10,'b     ',F11.3,' K A**3 '//,' ')
6050 format(//,T5,'FIRST ENERGY BALANCE: ',                    //     &
              ,T5,'Number of particles ..........: ',F18.6,'  ',/,    &
               T5,'TOTAL   energy (He) ..........: ',F18.6,' K',/,    &
               T5,'Energy per particle (He) .....: ',F18.6,' K',/,    &
               T5,'Kinetic energy (He) ..........: ',F18.6,' K',/,    &
               T5,'Lennard-Jones energy (He) ....: ',F18.6,' K',/,    &
               T5,'Alpha_s term  energy (He) ....: ',F18.6,' K',/,    &
               T5,'Solid energy (He)  ...........: ',F18.6,' K',/,    &
               T5,'Correlation energy   (He) ....: ',F18.6,' K')
6060 format(1x,T5,'Impurity energy (X) ..........: ',F18.6,' K',/,    &
               T5,'Kinetic energy (X) ...........: ',F18.6,' K',/,    &
               T5,'Interaction energy (X-He) ....: ',F18.6,' K',/,    &
               T5,'Spin-Orbit energy (X) ........: ',F18.6,' K',/,    &
               T5,'TOTAL ENERGY (He+X) ..........: ',F18.6,' K',/)
6065 format(1x,T5,'Impurity location  (x-axix) ..: ',F18.6,' A',/,    &
               T5,'                   (y-axis) ..: ',F18.6,' A',/,    &
               T5,'                   (z-axis) ..: ',F18.6,' A',/)
7000 format(//,1x,T2,'Iter     Mu(K)      Err(Mu)    Ttime  / Lap Time',/,&
 '--------------------------------------------------')
7001 format(//,1x,T2, &
     'Iter     Mu(K)      Err(Mu)    Autovalue(K)   err(K)   ETtime  / Lap Time',&
     /,74('-'))

7010 format(//,T5,'ITERATIVE PROCEDURE ',                                    //  &
               T5,'Number of particles ....... ',0P,F18.6,'  ',/,                &
               T5,'Total Energy (He).......... ',0P,F18.6,' K +- ',1P,e12.4,' K',&
             /,T5,'Energy per particle (He)... ',0P,F18.6,' K',/,                &
             /,T5,'Kinetic Energy (He)........ ',0P,F18.6,' K',                  &
             /,T5,'Lennard-Jones Energy (He).. ',0P,F18.6,' K',                  &
             /,T5,'Alpha_s term  energy (He).. ',0P,F18.6,' K',                  &
             /,T5,'Solid energy (He) ......... ',0P,F18.6,' K',                  &             
             /,T5,'Correlation Energy  (He)... ',0P,F18.6,' K')
7015 format(   T5,'Impurity energy (X->He) ... ',0P,F18.6,' K',/,&
               T5,'Kinetic energy (X) ........ ',0P,F18.6,' K',/,&
               T5,'Interaction energy (X-He) . ',0P,F18.6,' K',/,    &
               T5,'Spin-Orbit energy (X) ..... ',0P,F18.6,' K',/,    &
               T5,'TOTAL energy (He+X) ....... ',0P,F18.6,' K',/,' ')

7016 format(1x,T5,'Impurity location  (x-axix) ..: ',F18.6,' K',/,    &
               T5,'                   (y-axis) ..: ',F18.6,' K',/,    &
               T5,'                   (z-axis) ..: ',F18.6,' K',/,' ')

7017 format(   T5,'Chemical Potential ........ ',0P,F18.6,' K +- ',1P,e12.4,'K')
7018 format(   T5,'Autovalue (impurity) ...... ',0P,F18.6,' K +- ',1P,e12.4,'K')

! 7100 format(/1x,T5,'Center of Mass of the Helium ...(', &
!                      0P,F10.6,',',F10.6,',',F10.6,') A')
7100 format(/1x,T5,'Center of Mass of the Helium ...(', &
                    0P,F9.4,',',F9.4,',',F9.4,') A')
7110 format(/1x,T5,'Center of Mass of the Helium ........(', &
                    0P,F9.4,',',F9.4,',',F9.4,') A',/,       &
             1x,T5,'Center of Mass of the Impurity ......(', &
                    0P,F9.4,',',F9.4,',',F9.4,') A',/,       &
             1x,T5,'Distances between centers of mass ...(', &
                    0P,F9.4,',',F9.4,',',F9.4,') A',/,' ')

7020 format(1x,1P,2E14.5)
7030 format(0P,I5,T7,F13.7,T21,1P,E9.2,T32,0P,F8.0,'/',F8.2)
7035 format(0P,I5,T7,F13.7,T21,1P,E9.2,T32,0P,F13.7,T46,1P,E9.2,T57,0P,F8.0,'/',F8.2)


8010 format('partial.density' ,i1)
8015 format('partial.densityx',i1)
8020 format('partial.density' ,i2)
8025 format('partial.densityx',i2)
8030 format('partial.density' ,i3)
8035 format('partial.densityx',i3)

8011 format('00',i1)
8012 format('0',i2)
8013 format(i3)

5010 format('density.',SS,i1,'.out')
5020 format('density.',SS,i2,'.out')
5030 format('density.',SS,i3,'.out')
5040 format('density.',SS,i4,'.out')
5050 format('density.',SS,i5,'.out')
5060 format('density.',SS,i6,'.out')

5015 format('densityx.',SS,i1,'.out')
5025 format('densityx.',SS,i2,'.out')
5035 format('densityx.',SS,i3,'.out')
5045 format('densityx.',SS,i4,'.out')
5055 format('densityx.',SS,i5,'.out')
5065 format('densityx.',SS,i6,'.out')
!         1         2         3         4         5         6         7         8
!|2345678901234567890123456789012345678901234567890123456789012345678901234567890

end program 
