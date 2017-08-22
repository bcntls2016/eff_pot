PROGRAM Eff_pot

IMPLICIT NONE

character (len = 40)	::	fileden = "den.inp"		! Input file containing the wavefunction
character (len = 1) 	:: cchar = "#"			! Character to be skipped on routine titols
character (len = 80)	:: selec = "Rb_plus_Fausto"	! Selects the potential
logical (kind = 4)	:: OMP_Dynamic_Enable = .false.
integer (kind = 4)	:: nx, ny, nz			! Number of points of the grid in each direction
integer (kind = 4)	:: ix, iy, iz			! Counters for grid loops
integer (kind = 4)	:: ninvar			! Number of components of lambda
integer (kind = 4)	:: isalto			! Number of lines to skip at reading
integer (kind = 4)	:: mode				! Calculation mode
real (kind = 8)		:: xmax, ymax, zmax	! Maximum values of x, y and z
real (kind = 8)		:: hx, hy, hz		! x, y and z step for the grid
real (kind = 8)		:: dxyz            	! Grid-element volume. (dxyz=hx*hy*hz)
real (kind = 8)		:: rimp(3)			! Vector storing the position of the impurity
real (kind = 8)		:: rXden(3)			! Vector from density element to impurity 
real (kind = 8)		:: vimp(3)			! Vector storing the velocity of the impurity
real (kind = 8)		:: zimp, cdm(20)	! Impurity location (Z)
real (kind = 8)		:: sumConvolution	! Result of the convolution
real (kind = 8)		:: r				! Distance to the grid origin
real (kind = 8)		:: v				! Module of the velocity
real (kind = 8)		:: vasymp			! Asymptotic velocity
real (kind = 8)		:: Ekin				! Kinetic Energy
real (kind = 8)		:: Select_Pot		! Function to select the potential
real (kind = 8)		:: umax				! Value of the potential at r=r_cutoff
real (kind = 8)		:: r_cutoff			! Distance for the cutoff
real (kind = 8)		:: mimpur			! Impurity mass (au)
real (kind = 8)		:: tevo				! Amount of time evolution needed in ps		
real (kind = 8)		:: dt				! Time step in ps between recorded wave functions		
real (kind = 8)		:: dz				! Impurity position change corresponding to a
										! 	z-distance traveled in dt at velocity v(3)
real (kind = 8)		:: zcontinue = 0	! z-position to continue from. Usually zero	
real (kind = 8)		:: zfin				! Total z-distances traveled after tevo	
real (kind = 8) , parameter	:: mp_u = 0.0207569277d0  	! proton mass in Angs**-2 * Kelvin**-1, go figure!
real (kind=8)	, parameter	::	Ktops=7.63822291d0
real (kind = 8), ALLOCATABLE	:: x(:), y(:), z(:)		! Values in X, Y and Z
real (kind = 8), ALLOCATABLE	:: den(:,:,:)			! Helium density
complex(kind = 8), ALLOCATABLE	:: invar(:)			! Lambda (internal variables)
complex(kind = 8), ALLOCATABLE	:: psi(:,:,:)			! Wave Function. Density = MOD(psi)**2

!$ CALL OMP_SET_DYNAMIC(OMP_Dynamic_Enable)

40 format(5e26.15)
namelist /input/ selec, fileden, mode, zcontinue, tevo, dt, r_cutoff, umax, mimpur                             
read(5,nml=input)
mimpur = mimpur * mp_u

open (unit=1, file=fileden)
call titols(1,cchar,isalto)
read(1,*) xmax,ymax,zmax,hx,hy,hz,nx,ny,nz,rimp,vimp,ninvar
ALLOCATE (x(nx))
ALLOCATE (y(ny))
ALLOCATE (z(nz))
ALLOCATE (invar(ninvar))
ALLOCATE (psi(nx,ny,nz))
ALLOCATE (den(nx,ny,nz))
read(1,*) invar
read(1,*) psi
close(1)

den = Conjg(psi) * psi
hx = 2.0d0*abs(xmax)/(nx)
hy = 2.0d0*abs(ymax)/(ny)
hz = 2.0d0*abs(zmax)/(nz)    
dxyz = hx*hy*hz
    
!$OMP PARALLEL
!$OMP DO PRIVATE(ix)    
DO ix=1,nx  !.................... Grid X
	x(ix) = -xmax + hx * (ix - 1)
END DO
!$OMP END DO NOWAIT

!$OMP DO PRIVATE(iy)    
DO iy=1,ny  !.................... Grid Y
	y(iy) = -ymax + hy * (iy - 1)
END DO
!$OMP END DO NOWAIT

!$OMP DO PRIVATE(iz)    
DO iz=1,nz  !.................... Grid  Z
	z(iz) = -zmax + hz * (iz - 1)
END DO
!$OMP END DO
!$OMP END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF (mode == 0) THEN		! Uses all the densities
	!$OMP PARALLEL PRIVATE(ix,iy,iz,rXden,r)
	!$OMP DO REDUCTION(+:sumConvolution)  
	DO iz=1,nz
		DO iy = 1,ny
			DO ix = 1,nx	    
				rXden(1) = rimp(1) - x(ix)
				rXden(2) = rimp(2) - y(iy)
				rXden(3) = rimp(3) - z(iz)
				r = SQRT(sum(rXden * rXden))
				sumConvolution = sumConvolution + (den(ix,iy,iz) * Select_Pot(selec,r,r_cutoff,umax))
			END DO
		END DO
	END DO
	!$OMP END DO
	!$OMP END PARALLEL    
	sumConvolution = sumConvolution * dxyz
	v = SQRT(sum(vimp * vimp))
	Ekin = 0.5 * mimpur * v * v
	vasymp = SQRT(2 * (Ekin + sumConvolution) / mimpur)
	WRITE(6,40) rimp(3) + zimp, v / Ktops, Ekin, sumConvolution, vasymp / Ktops
	CALL FLUSH(6)
END IF 
  
If (mode == 1) THEN !Mode == 1, uses only one density and changes the z-axis position of the impurity by hand
	dz = vimp(3) / Ktops * dt
	zfin = tevo / dt * dz + zcontinue
	OPEN(unit = 10, file = "status.dat")
	WRITE(6,*) "# rimp(3)+zimp [A], v* [A/ps], Ekin [K], sumConvolution [K], v+ [A/ps]"
	CALL FLUSH(6)
	DO zimp = zcontinue + dz, zfin, dz    
		sumConvolution = 0.d0
		!$OMP PARALLEL PRIVATE(ix,iy,iz,rXden,r)
		!$OMP DO REDUCTION(+:sumConvolution)
		DO iz=1,nz
			DO iy = 1,ny
				DO ix = 1,nx
					rXden(1) = rimp(1) - x(ix)
					rXden(2) = rimp(2) - y(iy)
					rXden(3) = rimp(3) - z(iz) + zimp 	      
					r = SQRT(sum(rXden * rXden))
					sumConvolution = sumConvolution + (den(ix,iy,iz) * Select_Pot(selec,r,r_cutoff,umax))
				END DO
			END DO
		END DO
		!$OMP END DO
		!$OMP END PARALLEL
		sumConvolution = sumConvolution * dxyz
		v = SQRT(sum(vimp * vimp))
		Ekin = 0.5 * mimpur * v * v
		vasymp = SQRT(2 * (Ekin + sumConvolution) / mimpur)
		WRITE(6,40) rimp(3) + zimp, v / Ktops, Ekin, sumConvolution, vasymp / Ktops
		CALL FLUSH(6)
		WRITE(10,*), zimp, " - DONE"
	END DO
	
	CLOSE(10)
END IF

END PROGRAM
