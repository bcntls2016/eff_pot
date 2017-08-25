subroutine potenimpini
use interpol !, only:DelInter,potpi,potdel,npot,vpi,delta
use grid
use classicimp  !, only: rimp
implicit none
integer (kind=4) :: i
real    (kind=8) :: r,rmax,V_Pi,V_Sigma, V_Delta, Aux, V_gs
real    (kind=8) :: r1,r2


npot =100*max(nx,ny,nz)+1
rmax = 4.d0*max(xmax,ymax,zmax)

DelInter = rmax/dfloat(npot-1)

!.....................................
! Interpolation 
!.....................................

allocate(potdel(0:npot+1))
If(Lfilter_exciplex_force)Then
  allocate(potHe_He(0:npot+1)); potHe_He = 1.d6
Endif        
allocate(Delta(nx,ny,nz))
If(LState.Eq.'P')Then
  Als=Als_P
  allocate(potpi(0:npot+1))
  allocate(  Vpi(nx,ny,nz))
Else
  Als=Als_D
  allocate(Pi_Del(nx,ny,nz)); allocate(Sig_Del(nx,ny,nz))
  allocate(potpidel(0:npot+1)); allocate(potsigdel(0:npot+1))
EndIf

do i=1,npot
  r  = dfloat(i-1)*DelInter
  If(LState.Eq.'P')Then
    Aux           = V_Pi(r)
    potpi(i)      = Aux
    potdel(i)     = V_Sigma(r) - Aux
  Else
    Aux           = V_Delta(r)
    potdel(i)     = Aux
    potpidel(i)   = V_Pi(r)    - Aux
    potsigdel(i)  = V_Sigma(r) - Aux
  Endif
  If(Lfilter_exciplex_force)Then
    potHe_He(i)   = V_gs(r)         ! We must define selec_gs='LJ_OT' or 'Aziz_He', and thier corresponding: r_cutoff_gs and umax_gs 
  Endif        
enddo

potdel(0)=potdel(1); potdel(npot+1)=potdel(npot)
If(LState.Eq.'P')Then
  potpi(0)=potpi(1); potpi(npot+1)=potpi(npot)
Else
  potpidel(0)=potpidel(1); potpidel(npot+1)=potpidel(npot)
  potsigdel(0)=potsigdel(1); potsigdel(npot+1)=potsigdel(npot)
EndIf
If(Lfilter_exciplex_force)Then
  potHe_He(0)   = potHe_He(1); potHe_He(npot+1)   = potHe_He(npot) 
!  Open(Unit=11, File='Pot_He_he.out')
!  Do i=1, npot
!    r  = dfloat(i-1)*DelInter
!    Write(11,'(1p,2E16.6)')r,potHe_He(i)
!  EndDo
Endif  

!.....................................

call updatepoten(rimp)

end subroutine potenimpini
