subroutine instates()

use grid
use rho
use classicimp
use energies
use field
use interpol

implicit none
Real (kind=8), Allocatable :: r_He_exc(:,:)
real    (kind=8) :: ar(10,10), ai(10,10), w(10), zr(10,10), zi(10,10), fv1(10), fv2(10), fm1(2,10)
complex (kind=8) :: Mat(10,10), SO(10,10), MatSO(10,10)
REAL (kind=8) :: d0, d1, d2, d3, d4, d5, sd0, sd1, sd2, sd3, sd4, sd5, dV1, dV2
REAL (kind=8) :: A,B,C,Vx,Vy,Vz,Vxy,Vxz,Vyz,dVpi,dVsig,dVdel, dum2,Dv,vls,Ep1,Ep2,Ep3,V0gs,U_ext
REAL (kind=8) :: armonicd(6), pi, twopi, fpi, dx, dy, dz, zx, yx, xx, yy, zz, r, r2, c1osqr2
Integer (kind=4)  :: ierr,i,j,ix,iy,iz,il,nm,np,nd,il1, k, kk, l, ll, njz, ki, N_He_exc,iaxis
Integer (kind=4),save  :: iaux
REAL (kind=8) :: sqr3, sqr3o2, insqr3, in3, difd3d5, sumd3d5, cm1_to_K = 1.439d0, Uext_new
REAL (kind=8) :: sumVdel(10,10), MP(10,10), MS(10,10), MatD(10,10), U_total
REAL (kind=8) :: Sto, Aux, V_pi, V_Sigma
COMPLEX (kind=8) ::  MatDSO(10,10), eigv(10,10), testv(10,10), invars(10), vlsaux, SO_aux(10), invar1(10)
complex (kind=8) :: auxEVEN, auxODD, auxc
complex   (kind=8) :: xxxr(3)
complex   (kind=8) :: uim = (0.d0,1.d0)

Data nm/10/, np/6/, nd/10/, iaux/0/

pi    = 4.0D0*DATAN(1.0D0)
twopi = 2.d0*pi
fpi   = 4.d0*pi
c1osqr2=1.0d0/Dsqrt(2.0d0)

invar = (0.0d0, 0.0d0)

If (Lstate=='P') then

  Als = Als_P
  Also2 = Als*0.5d0

!......................................
!... Dimensionamos algunas matrices ...
!......................................

  rimp(1)=ximp; rimp(2)=yimp; rimp(3)=zimp

  Vx   = 0.d0
  Vy   = 0.d0
  Vz   = 0.d0
  Vxy  = 0.d0
  Vxz  = 0.d0
  Vyz  = 0.d0

  If(Lexciplex_state_fix)Then
    If(Trim(Exciplex).Eq.'Ring')then
      N_He_exc=6
    ElseIf(Trim(Exciplex).Eq.'Linear')Then
      N_He_exc=2
    Else
      Write(6,'("From Instates: Exciplex variable is not correct...",A)')Exciplex
      Stop 'From instates  001'
    Endif
    Allocate(r_He_exc(N_He_exc,3))
    iaxis = 0
    If(Lexciplex_axis.Eq.'X')iaxis=1
    If(Lexciplex_axis.Eq.'Y')iaxis=2
    If(Lexciplex_axis.Eq.'Z')iaxis=3
    If(iaxis.Eq.0)Then
      Write(6,'("From Instates: El valor de iaxis no pot ser 0")')
      Stop 'Instates ---> Iaxis = 0'
    Endif
    If(iaxis.Eq.3)Then
      ix=1
      iy=2
    ElseIf(iaxis.Eq.2)Then
      ix=1
      iy=3
    ElseIf(iaxis.Eq.1)Then
      ix=2
      iy=3
    Endif
    r_He_exc = 0.0d0
    If(N_He_exc.Eq.2)Then
      r_He_exc(1,iaxis)=r_exc
      r_He_exc(2,iaxis)=-r_exc
    Else
      Sto = 2.0d0*Pi/Dfloat(N_He_exc)
      Do i = 1, N_He_exc
         Aux =(i-1)*Sto
         r_He_exc(i,ix) = r_exc*Dcos(Aux)
         r_He_exc(i,iy) = r_exc*Dsin(Aux)
      EndDo
    EndIf
    Do i = 1, N_He_exc
!      xx=r_He_exc(i,1) - rimp(1)
!      yy=r_He_exc(i,2) - rimp(2)
!      zz=r_He_exc(i,3) - rimp(3)
      xx=r_He_exc(i,1) 
      yy=r_He_exc(i,2) 
      zz=r_He_exc(i,3) 
      r2=xx**2+yy**2+zz**2
      r=dsqrt(r2)
      aux=V_pi(r)
      dVpi = aux
      dV = V_Sigma(r) - aux
      Vx  = Vx  + (dVpi+dV*(xx/r)**2)
      Vy  = Vy  + (dVpi+dV*(yy/r)**2)
      Vz  = Vz  + (dVpi+dV*(zz/r)**2)
      Vxy = Vxy + dV*xx*yy/r2
      Vxz = Vxz + dV*xx*zz/r2
      Vyz = Vyz + dV*yy*zz/r2
    EndDo
    DeAllocate(r_He_exc)

  Else

    dx   = hx
    dy   = hy
    dz   = hz
    A = 0.0d0
    do iz=1,nz
      zz = z(iz) - rimp(3)
      do iy=1,ny
        yy = y(iy) - rimp(2)
        do ix=1,nx
          xx = x(ix) - rimp(1)
            dVpi = Vpi(ix,iy,iz)
            dum2=den(ix,iy,iz)
            dV=Delta(ix,iy,iz)     ! Delta ja esta dividi per r**2
            A = A + dum2
            Vx  = Vx  + dum2*(dVpi+dV*xx**2)
            Vy  = Vy  + dum2*(dVpi+dV*yy**2)
            Vz  = Vz  + dum2*(dVpi+dV*zz**2)
            Vxy = Vxy + dum2*dV*xx*yy
            Vxz = Vxz + dum2*dV*xx*zz
            Vyz = Vyz + dum2*dV*yy*zz
        end do
      end do
    end do
    Vx  = Vx*dxyz
    Vy  = Vy*dxyz
    Vz  = Vz*dxyz
    Vxy = Vxy*dxyz
    Vxz = Vxz*dxyz
    Vyz = Vyz*dxyz
    A   = A*dxyz
    If(Lprint_invar)write(6,'(" Atomos de helio",1p,E18.10)') A
  EndIf
!
!  Matriz V^jiss' doble
!
  Mat=(0.d0,0.d0)
  
  Mat(1,1) = Vx
  Mat(2,2) = Vx
  Mat(3,3) = Vy
  Mat(4,4) = Vy
  Mat(5,5) = Vz
  Mat(6,6) = Vz
  Mat(1,3) = Vxy
  Mat(2,4) = Vxy
  Mat(1,5) = Vxz
  Mat(2,6) = Vxz
  Mat(3,5) = Vyz
  Mat(4,6) = Vyz
  
  do i=1,np
    do j=1,np
      Mat(j,i)=Mat(i,j)
    end do
   end do
!
!Matriz SO
!
  SO=(0.d0,0.d0)
 
  SO(1,3)=-uim
  SO(2,4)=uim
  SO(2,5)=(-1.d0,0.d0)
  SO(1,6)=(1.d0,0.d0)
  SO(4,5)=-uim
  SO(3,6)=-uim
 
  SO(3,1)=uim
  SO(4,2)=-uim
  SO(5,2)=(-1.d0,0.d0)
  SO(6,1)=(1.d0,0.d0)
  SO(5,4)=uim
  SO(6,3)=uim
!
!Suma matriz V y matriz SO. Llamada a rutina de diagonalización
!
  MatSO=Mat+(Also2*SO)
  ar=real(MatSO)
  ai=aimag(MatSO)
  call ch(nm,np,ar,ai,w,1,zr,zi,fv1,fv2,fm1,ierr)
  If(Lprint_invar)write(6,'(" Error",1p,E18.10)') ierr
  do il = 1,np
    invar = cmplx(0.d0, 0.d0)
    If(Lprint_invar)write(6,'(" Eigenvalues",I5,1p,E18.10)') il,w(il)
    do i = 1,np
      invar(i) = cmplx(zr(i,il),zi(i,il))
    end do
!
! We store all the eigenvectors for diagonalize J_z if it is needed
!
    eigv(il,:) = invar(:)

    invars = conjg(invar)
    vls = ( 2.d0*real(invars(1)*invar(6)-invars(2)*invar(5))                                     + &
               2.d0*aimag(invars(1)*invar(3) + invars(4)*invar(2) + invars(3)*invar(6) + invars(4)*invar(5) ))
    vls=vls*0.5d0
               
    If(Lprint_invar)write(6,'("Valor esperado de L·S, ESO", 1p, 2E15.6)') vls, (vls*Als)
    do i = 1,np
      invars(i) = (0.d0,0.d0)
      do j = 1,np
        invars(i) = invars(i) + SO(i,j)*invar(j)
      end do
    end do
    vlsaux=(0.d0,0.d0)
    do i = 1,np
      vlsaux = vlsaux + invars(i)*conjg(invar(i))
    end do
    vlsaux=vlsaux*0.5d0
    If(Lprint_invar)write(6,'("Valor esperado de L·S matrices, ESO", 1p, 4E15.6)') vlsaux, (vlsaux*Als)
  end do
  If(Ldiag_jz)Then
    If(Instate.Eq.0)Then
      ki=0; njz=2
    Else
      ki=2; njz=4
    Endif
      MatSO = cmplx(0.0d0, 0.0d0)
      Do kk=1, njz
        k = kk + ki
        Do ll=1, njz
          l = ll + ki
          MatSO(kk,ll) = 0.5d0*conjg(eigv(k,1))*eigv(l,1)          &
                       - uim*conjg(eigv(k,1))*eigv(l,3)            &
                       - 0.5d0*conjg(eigv(k,2))*eigv(l,2)          &
                       - uim*conjg(eigv(k,2))*eigv(l,4)            &
                       + uim*conjg(eigv(k,3))*eigv(l,1)            &
                       + 0.5d0*conjg(eigv(k,3))*eigv(l,3)          &
                       + uim*conjg(eigv(k,4))*eigv(l,2)            &
                       - 0.5d0*conjg(eigv(k,4))*eigv(l,4)          &
                       + 0.5d0*conjg(eigv(k,5))*eigv(l,5)          &
                       - 0.5d0*conjg(eigv(k,6))*eigv(l,6)          
        EndDo
      EndDo
      ar=real(MatSO)
      ai=aimag(MatSO)
      call ch(nm,njz,ar,ai,w,1,zr,zi,fv1,fv2,fm1,ierr)
      Do il = 1, njz
        If(Lprint_invar)write(6,'(" Eigenvalues",I5,1p,E18.10)') il,w(il)
      EndDo
!
!  Where, we select the state that give an approiate <jz> value
!
      il=0
      If(Instate.Eq.0)Then
        If(Ljz.Eq.'')Then
          Ljz='-1/2'
        Endif
        If(Ljz.Eq.'-1/2')il=1
        If(Ljz.Eq.'+1/2')il=2
      Else
        If(Ljz.Eq.'')Then
          If(Instate.Eq.1)Ljz='-3/2'
          If(Instate.Eq.2)Ljz='-1/2'
        Endif
        If(Ljz.Eq.'-3/2')il=1
        If(Ljz.Eq.'-1/2')il=2
        If(Ljz.Eq.'+1/2')il=3
        If(Ljz.Eq.'+3/2')il=4
      Endif
      If(il.Eq.0)Then
        Write(6,'("Wrong Ljz value...:",A)')Ljz
        Stop 'Wrong Ljz value..'
      Endif
      Do i=1,np
        auxc = (0.d0, 0.d0)
        Do kk=1, njz
          k= kk + ki
          auxc = auxc + Cmplx(zr(kk,il),zi(kk,il))*eigv(k,i)
        EndDo
        invar(i) = auxc
      Enddo
!
!  We compute < invar | Jz | invar >
!
      Auxc  =  0.5d0*Abs(Invar(1))**2 - 0.5d0*Abs(Invar(2))**2             &   !
             + 0.5d0*Abs(Invar(3))**2 - 0.5d0*Abs(Invar(4))**2             &   !  Termes diagonals
             + 0.5d0*Abs(Invar(5))**2 - 0.5d0*Abs(Invar(6))**2             &   !
!
             - uim*Conjg(Invar(1))*Invar(3) - uim*Conjg(Invar(2))*Invar(4)  &   !
             + uim*Conjg(Invar(3))*Invar(1) + uim*Conjg(Invar(4))*Invar(2)      !  Termes no diagonals
                                                                                !

      If(Lprint_invar)Write(6,'(" Control for selected <Jz>...:",A,"...>:",1p,2E15.6)')Ljz,Auxc
  ElseIf(Lexcite_state_fix)Then
!
!   Aqui fijamos el estado sin diagonalizar
!
    Invar = (0.0d0, 0.0d0)
    If(Instate.Eq.0)Then
!
!   Seleccionamos |x up>
!
      invar(1) = (1.0d0,0.0d0)
    ElseIf(Instate.Eq.1)Then
!
!   Seleccionamos 1/Sqrt(2)*(|x up> +|y up>)
!
      invar(1) = (1.0d0,0.0d0)/Dsqrt(2.0d0)
      invar(3) = (1.0d0,0.0d0)/Dsqrt(2.0d0)
    Else
!
!   Seleccionamos 1/Sqrt(3)*(|x up> +|y up> +|z up>)
!
      invar(1) = (1.0d0,0.0d0)/Dsqrt(3.0d0)
      invar(3) = (1.0d0,0.0d0)/Dsqrt(3.0d0)
      invar(5) = (1.0d0,0.0d0)/Dsqrt(3.0d0)
      If(Lprint_invar)Then
        Write(6,'("Hemos seleccionado un estado sin diagonalizar")')
      Endif
    Endif
  ElseIf(Laverage_P_value)Then
    Aux=1.0d0/Dsqrt(Dfloat(np))
    invar = (0.0d0, 0.0d0)
    Do il = 1, np
      do i = 1,np
        invar(i) = invar(i) + Aux*cmplx(zr(i,il),zi(i,il))
      end do
    end do
  ElseIf(Lmixture)Then
    !      
    !  Instate =0 corresponds to Pi_1o2
    !  Instate =1 corresponds to Pi_3o2
    !  Instate =2 corresponds to Sigma_1o2
    !    
    !  il  = 1  + 2*instate
    !
    il  = 1  + 2
    il1 = il + 1
      do i = 1,np
        invar(i) = DSqrt(a_Pi_3o2)*(c1osqr2*cmplx(zr(i,il),zi(i,il)) + c1osqr2*cmplx(zr(i,il1),zi(i,il1)))
      end do
    il  = 1  + 4
    il1 = il + 1
      do i = 1,np
        invar(i) = invar(i) + Dsqrt(a_Sig_1o2)*(c1osqr2*cmplx(zr(i,il),zi(i,il)) + c1osqr2*cmplx(zr(i,il1),zi(i,il1)))
      end do
    il  = 1  + 0
    il1 = il + 1
      do i = 1,np
        invar(i) = invar(i) + Dsqrt(Max((1.d0-a_Pi_3o2-a_Sig_1o2),0.d0))*(c1osqr2*cmplx(zr(i,il),zi(i,il)) + c1osqr2*cmplx(zr(i,il1),zi(i,il1)))
      end do
  Else
    il  = 1  + 2*instate
    il1 = il + 1
      do i = 1,np
        invar(i) = c1osqr2*cmplx(zr(i,il),zi(i,il)) + c1osqr2*cmplx(zr(i,il1),zi(i,il1))
      end do
  Endif
  If(Lprint_invar)Then
    Write(6,'("Norma del estado sleccionado....:",1p,2E15.6)')Sum(invar*conjg(invar))
  Endif
  invar1 = invar
!
!! polinomio caracteristico para encontrar los valores propios:
!! t^3 + C t^2 + B t + A == 0
!
  A = Vxz**2*Vy-2*Vxy*Vxz*Vyz+Vx*Vyz**2+Vxy**2*Vz-Vx*Vy*Vz
  B = -Vxy**2-Vxz**2-Vyz**2+Vx*Vy+Vx*Vz+Vy*Vz
  C = -(Vx+Vy+Vz)
!
!! adicion del termino de spin-orbita
!! Delta = 3/2 Als
  A = A + 0.25d0*(Als**3-C*Als**2)
  B = B - 0.75d0*Als**2
  C = C
!! valores propios:
  call cubic(1.d0,C,B,A,xxxr)
  Ep1=real(xxxr(1))
  Ep2=real(xxxr(2))
  Ep3=real(xxxr(3))
  If(Lprint_invar)write(6,'(" Ep1,EP2,Ep3,V0gs",1p,4E18.10)') Ep1,Ep2,Ep3,V0gs
!
!   Comprobació del estat selecionat
!
  If(Lprint_invar)Write(6,'(//,"Valor seleccionat",//)')
  invars = conjg(invar)
  vls = ( 2.d0*real(invars(1)*invar(6)-invars(2)*invar(5))                                     + & 
               2.d0*aimag(invars(1)*invar(3) + invars(4)*invar(2) + invars(3)*invar(6) + invars(4)*invar(5) ))
  vls=vls*0.5d0 
  If(Lprint_invar)write(6,'("Valor esperado de L·S, ESO", 1p, 2E15.6)') vls, (vls*Als)
  do i = 1,np
    invars(i) = (0.d0,0.d0)
    do j = 1,np
      invars(i) = invars(i) + SO(i,j)*invar(j)
    end do
  end do
  vlsaux=(0.d0,0.d0)
  do i = 1,np
    vlsaux = vlsaux + invars(i)*conjg(invar(i))
  end do  
  vlsaux=vlsaux*0.5d0 
  If(Lprint_invar)write(6,'("Valor esperado de L·S matrices, ESO", 1p, 4E15.6)') vlsaux, (vlsaux*Als)
!
!  Calculate uext using internal state
!
  do iz=1,nz
   zz = z(iz)-rimp(3)
   do iy=1,ny
     yy = y(iy)-rimp(2)
     do ix=1,nx
      xx = x(ix)-rimp(1)
!      r2=xx**xx+yy*yy+zz*zz
!      r=sqrt(r2)
      auxODD =invar(1)*xx + invar(3)*yy +invar(5)*zz
      auxEVEN=invar(2)*xx + invar(4)*yy +invar(6)*zz
      Uext_new = VPi(ix,iy,iz) + Delta(ix,iy,iz)*(auxODD*conjg(auxODD) + auxEVEN*conjg(auxEVEN))
      uimp(ix,iy,iz) = Uext_new
    enddo
   enddo
  enddo
!  Lfirst=.false.
!
!  Calculate Spin-Orbit energy contribution
!
  invars(:)=conjg(invar(:))
  eso = Also2*( 2.d0*real(invars(1)*invar(6)-invars(2)*invar(5))                                     + &
               2.d0*aimag(invars(1)*invar(3) + invars(4)*invar(2) + invars(3)*invar(6) + invars(4)*invar(5) ) )


  If(Lprint_invar)write(6,'("Spin-Orbit energy contribution", 1p, 4E15.6)') eso


Else
! !
! ! Calculo de la matriz D
! !
  Als   = Als_D
  Also2 = Als*0.5d0
!
! Constantes e inicialización de matrices
!
  sqr3=sqrt(3.d0)
  sqr3o2=sqrt(3.d0)/2.d0
  insqr3=1.d0/sqrt(3.d0)
  in3=1.d0/3.d0
  MP=0.d0
  MS=0.d0
  sumVdel = 0.d0
!
! Inicio bucle que recorre la malla de cálculo
!
  do iz=1,nz
    zz=z(iz) - rimp(3)
    do iy=1,ny
      yy=y(iy) - rimp(2)
      do ix=1,nx
        xx=x(ix) - rimp(1)
!        r=dsqrt(xx**2+yy**2+zz**2)
!        if(r.Eq.0.d0) Cycle
        dum2=den(ix,iy,iz)
        dVdel= Delta(ix,iy,iz)*dum2
    !dV1=(dVsig-dVdel)
        dV1=Pi_Del(ix,iy,iz)*dum2     !   Lleva un r**4 dividiendo
    !dV2=(dVpi-dVdel)                 !
        dV2=Sig_Del(ix,iy,iz)*dum2    !     "   "   "       "
!
!  Armónicos Esféricos
!
        d0=xx*xx+yy*yy+zz*zz
        d1=sqr3*xx*yy
        d2=sqr3*yy*zz
        d3=(-xx*xx-yy*yy+2*zz*zz)/2.
        d4=sqr3*xx*zz
        d5=sqr3o2*(xx*xx-yy*yy)
        armonicd = (/ d1, d2, d3, d4, d5, d0 /)
!
!   Definición de variables para facilitar escritura
!
        sd0=d0*d0
        sd1=d1*d1
        sd2=d2*d2
        sd3=d3*d3
        sd4=d4*d4
        sd5=d5*d5
        sumd3d5 = d3 + insqr3*d5
        difd3d5 = d3 - insqr3*d5
!
!Definición de los elementos del triángulo superior la matriz (5x5)
!convolucionados con la densidad (rho4)
!
        MP(1,1) = MP(1,1) + (in3*(sd2+sd4+4.d0*sd5))*dV2
        MP(2,2) = MP(2,2) + (in3*(sd1+sd4) + sumd3d5*sumd3d5)*dV2
        MP(3,3) = MP(3,3) + (sd4 + sd2)*dV2
        MP(4,4) = MP(4,4) + (in3*(sd1+sd2) + difd3d5*difd3d5)*dV2
        MP(5,5) = MP(5,5) + (in3*(sd2+sd4+4.d0*sd1))*dV2
!
        MP(1,2) = MP(1,2) + (-4.d0*in3*d1*d2 + insqr3*d4*d0)*dV2
        MP(1,3) = MP(1,3) + (-2.d0*insqr3*d2*d4)*dV2
        MP(1,4) = MP(1,4) + (-4.d0*in3*d1*d4 + insqr3*d2*d0)*dV2
        MP(1,5) = MP(1,5) + (-4.d0*in3*d1*d5)*dV2
!
        MP(2,3) = MP(2,3) + (insqr3*d1*d4 - d2*sumd3d5)*dV2
        MP(2,4) = MP(2,4) + (-4.d0*in3*d2*d4 + insqr3*d1*d0)*dV2
        MP(2,5) = MP(2,5) + (-d1*d4 - insqr3*d2*sumd3d5)*dV2
!
        MP(3,4) = MP(3,4) + (insqr3*d1*d2 - d4*difd3d5)*dV2
        MP(3,5) = MP(3,5) + (insqr3*(sd2-sd4))*dV2
!
        MP(4,5) = MP(4,5) + (d1*d2 + insqr3*d4*difd3d5)*dV2
        do i=1,nd/2
          do j=i,nd/2
            MS(i,j) = MS(i,j) + (armonicd(i)*armonicd(j))*dV1
            IF (i==j) sumVdel(i,j) = sumVdel(i,j) + dVdel
          end do
        end do
     end do
   end do
  end do
!
!Construcción de las matrices completas (5x5)
!
  do i=1,nd/2
    do j=i,nd/2
      MP(i,j) = MP(i,j)*dxyz
      MP(j,i) = MP(i,j)
      MS(i,j) = MS(i,j)*dxyz
      MS(j,i) = MS(i,j)
      sumVdel(i,j) = sumVdel(i,j)*dxyz
    end do
  end do
!
! V_delta + (V_sigma - V_delta)*MS + (V_pi - V_delta)*MP
!
  MP = sumVdel + MS + MP
!
!Producto de Kronecker de las matriz anterior (5x5)
!por la matriz unidad 2x2, creando una matriz (10x10)
!
  MatD = 0.d0
  do i=1,nd/2
    do j=1,nd/2
      MatD(2*i-1,2*j-1) = MP(i,j)
      MatD(2*i,2*j) = MP(i,j)
    end do
  end do
!
! Definición elementos matriz Spin-órbita (10x10)
!
  SOD = (0.d0,0.d0)

  SOD(1,4)  = (1.d0, 0.d0)
  SOD(1,8)  = -uim
  SOD(1,9)  = 2.d0*uim
  SOD(2,3)  = (-1.d0, 0.d0)
  SOD(2,7)  = -uim
  SOD(2,10) = -2.d0*uim
  SOD(3,6)  = -uim*sqr3
  SOD(3,7)  = uim
  SOD(3,10) = -uim
  SOD(4,5)  = -uim*sqr3
  SOD(4,8)  = -uim
  SOD(4,9)  = -uim
  SOD(5,8)  = -sqr3
  SOD(6,7)  = sqr3
  SOD(7,10) = (-1.d0, 0.d0)
  SOD(8,9)  = (1.d0, 0.d0)
  do i=1,nd
    do j=i,nd
      SOD(j,i) = conjg(SOD(i,j))
    enddo
  enddo
!
!Suma matriz MatD con la de Spin-órbita (10x10) y llamada
!a rutina de diagonalización 
!
  MatDSO = MatD + (Also2*SOD)
  ar=real(MatDSO)
  ai=aimag(MatDSO)
 
  call ch(nm,nd,ar,ai,w,1,zr,zi,fv1,fv2,fm1,ierr)

  do il = 1,nd
    If(Lprint_invar)write(6,'(" Eigenvalues",I5,1p,E18.10)') il,w(il)
    do i = 1,nd
      invar(i) = cmplx(zr(i,il),zi(i,il))
    end do
!
! We store all the eigenvectors for diagonalize J_z if it is needed
!
      eigv(il,:) = invar(:)

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
    eso=(vlsaux*Als)
   If(Lprint_invar)write(6,'("Valor esperado de L·S matrices, ESO", 1p, 4E15.6)') vlsaux, eso
  end do
  If(Ldiag_jz)Then
    If(Instate.Eq.0.Or.Instate.Eq.1)Then
      njz=4    !
      ki=0     ! Corresponde a los estado de j=3/2
               !
!
!  Where, we select the state that give an approiate <jz> value
!
      il=0
      If(Ljz.Eq.'')Then
        If(Instate.Eq.0)Ljz='-3/2'
        If(Instate.Eq.1)Ljz='-1/2'
      Endif
      If(Ljz.Eq.'-3/2')il=1
      If(Ljz.Eq.'-1/2')il=2
      If(Ljz.Eq.'+1/2')il=3
      If(Ljz.Eq.'+3/2')il=4
      If(il.Eq.0)Then
        Write(6,'("You have selected a wrong state..:",A,"Instate=",I4)')Ljz,instate
        Stop 'From instates:-(001)-Worng selection of Ljz'
      Endif
    Else
      njz=6    !
      ki=4     ! Corresponde a los estado de j=5/2
               !
!
!  Where, we select the state that give an approiate <jz> value
!
      il=0
      If(Ljz.Eq.'')Then
        If(Instate.Eq.2)Ljz='-5/2'
        If(Instate.Eq.3)Ljz='-3/2'
        If(Instate.Eq.4)Ljz='-1/2'
      Endif
      If(Ljz.Eq.'-5/2')il=1
      If(Ljz.Eq.'-3/2')il=2
      If(Ljz.Eq.'-1/2')il=3
      If(Ljz.Eq.'+1/2')il=4
      If(Ljz.Eq.'+3/2')il=5
      If(Ljz.Eq.'+5/2')il=6
      If(il.Eq.0)Then
        Write(6,'("You have selected a wrong state..:",A,"Instate=",I4)')Ljz,instate
        Stop 'From instates:-(002)-Worng selection of Ljz'
      Endif
    Endif
    MatSO = cmplx(0.0d0, 0.0d0)
    Do kk=1,njz
      k=kk+ki
      Do ll=1, njz
        l=ll+ki
        MatSO(kk,ll) = 0.5d0*conjg(eigv(k,9))*eigv(l,9)          &
                     - 2.0d0*uim*conjg(eigv(k,9))*eigv(l,1)      &
                     - 0.5d0*conjg(eigv(k,10))*eigv(l,10)        &
                     - 2.0d0*uim*conjg(eigv(k,10))*eigv(l,2)     &
                     + 0.5d0*conjg(eigv(k,3))*eigv(l,3)          &
                     + uim*conjg(eigv(k,3))*eigv(l,7)            &
                     - 0.5d0*conjg(eigv(k,4))*eigv(l,4)          &
                     + uim*conjg(eigv(k,4))*eigv(l,8)            &
                     + 0.5d0*conjg(eigv(k,5))*eigv(l,5)          &
                     - 0.5d0*conjg(eigv(k,6))*eigv(l,6)          &
                     - uim*conjg(eigv(k,7))*eigv(l,3)            &
                     + 0.5d0*conjg(eigv(k,7))*eigv(l,7)          &
                     - uim*conjg(eigv(k,8))*eigv(l,4)            &
                     - 0.5d0*conjg(eigv(k,8))*eigv(l,8)          &
                     + 2.0d0*uim*conjg(eigv(k,1))*eigv(l,9)      &
                     + 0.5d0*conjg(eigv(k,1))*eigv(l,1)          &
                     + 2.0d0*uim*conjg(eigv(k,2))*eigv(l,10)     &
                     - 0.5d0*conjg(eigv(k,2))*eigv(l,2)
      EndDo
    EndDo
    ar=real(MatSO)
    ai=aimag(MatSO)
    call ch(nm,njz,ar,ai,w,1,zr,zi,fv1,fv2,fm1,ierr)
    Do i = 1, njz
      If(Lprint_invar)write(6,'(" Eigenvalues",I5,1p,E18.10)') i,w(i)
    EndDo
    Do i=1,nd
      auxc = (0.d0, 0.d0)
      Do k=1, njz
          kk=k+ki
          auxc = auxc + Cmplx(zr(k,il),zi(k,il))*eigv(kk,i)
      EndDo
      invar(i) = auxc
    Enddo
!
!  We compute < invar | Jz | invar >
!
    Auxc = 0.5d0*conjg(invar(9))*invar(9)          &
         - 2.0d0*uim*conjg(invar(9))*invar(1)      &
         - 0.5d0*conjg(invar(10))*invar(10)        &
         - 2.0d0*uim*conjg(invar(10))*invar(2)     &
         + 0.5d0*conjg(invar(3))*invar(3)          &
         + uim*conjg(invar(3))*invar(7)            &
         - 0.5d0*conjg(invar(4))*invar(4)          &
         + uim*conjg(invar(4))*invar(8)            &
         + 0.5d0*conjg(invar(5))*invar(5)          &
         - 0.5d0*conjg(invar(6))*invar(6)          &
         - uim*conjg(invar(7))*invar(3)            &
         + 0.5d0*conjg(invar(7))*invar(7)          &
         - uim*conjg(invar(8))*invar(4)            &
         - 0.5d0*conjg(invar(8))*invar(8)          &
         + 2.0d0*uim*conjg(invar(1))*invar(9)      &
         + 0.5d0*conjg(invar(1))*invar(1)          &
         + 2.0d0*uim*conjg(invar(2))*invar(10)     &
         - 0.5d0*conjg(invar(2))*invar(2)

    If(Lprint_invar)Write(6,'(" Control for selected <Jz>...:",A,"...>:",1p,2E15.6)')Ljz,Auxc
!    Stop 'Temporal stop'
  Else
    il  =  1 + 2*instate
    il1 = il + 1
    do i = 1,nd
      invar(i) = c1osqr2*cmplx(zr(i,il),zi(i,il)) + c1osqr2*cmplx(zr(i,il1),zi(i,il1))
    end do
  Endif
!
!
  invar1 = invar   ! Guardem l'estat selecionat
!
!
  if(Lprint_invar) then    
    SO_aux=(0.d0,0.d0)
    do il = 1,nd
      do i = 1,nd
        invar(i) = cmplx(zr(i,il),zi(i,il))
      end do
      do i = 1,nd
        invars(i) = (0.d0,0.d0)
        Do j = 1,nd
          invars(i) = invars(i) + SOD(i,j)*invar(j)
        End do
      end do
      do i = 1,nd
        SO_aux(il)= SO_aux(il) + invars(i)*conjg(invar1(i))
      end do   
    end do
     write(6,'("Vector Li*SO*L_selec")')
    do i = 1,nd
     write(6,'(1p,2E15.6)') (SO_aux(i)*Also2)
    end do
  end if
!
!
  invar = invar1 ! Recuperem l'estat selecionat
!
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
  eso=(vlsaux*Als)
  If(Lprint_invar) then
    Write(6,'(//," Valor seleccionat",//)')
    write(6,'("Valor esperado de L·S matrices, ESO", 1p, 4E15.6)') vlsaux, eso
  end if
!!!!!!!!!!!!!!!!!!!!!!!!
!! Cálculo de Uext_D  !!
!!!!!!!!!!!!!!!!!!!!!!!!
  U_ext = 0.d0
  uimp=0.d0
  do iz=1,nz
    zz=z(iz) - rimp(3)
    do iy=1,ny
      yy=y(iy) - rimp(2)
      do ix=1,nx
        xx=x(ix) - rimp(1)
        MP=0.d0
        MS=0.d0
        MatD=0.d0
        sumVdel=0.d0
        r=dsqrt(xx**2+yy**2+zz**2)
!        if(r.eq.0.d0)r=1.d-8
        dVdel= Delta(ix,iy,iz)
    !dV1=(dVsig-dVdel)
        dV1=Pi_Del(ix,iy,iz)     !  Lleva un r**4 dividiendo
    !dV2=(dVpi-dVdel)            !
        dV2=Sig_Del(ix,iy,iz)    !    "    "  "       "
!
!Armónicos Esféricos
!
        d0=xx*xx+yy*yy+zz*zz
        d1=sqr3*xx*yy
        d2=sqr3*yy*zz
        d3=(-xx*xx-yy*yy+2*zz*zz)/2.
        d4=sqr3*xx*zz
        d5=sqr3o2*(xx*xx-yy*yy)
        armonicd = (/ d1, d2, d3, d4, d5, d0 /)
!
!Definición de variables para facilitar escritura
!
        sd0=d0*d0
        sd1=d1*d1
        sd2=d2*d2
        sd3=d3*d3
        sd4=d4*d4
        sd5=d5*d5
        sumd3d5 = d3 + insqr3*d5
        difd3d5 = d3 - insqr3*d5
!
!Definición de los elementos del triángulo superior la matriz (5x5)
!convolucionados con la densidad (rho4)
!
        MP(1,1) = (in3*(sd2+sd4+4.d0*sd5))*dV2
        MP(2,2) = (in3*(sd1+sd4) + sumd3d5*sumd3d5)*dV2
        MP(3,3) = (sd4 + sd2)*dV2
        MP(4,4) = (in3*(sd1+sd2) + difd3d5*difd3d5)*dV2
        MP(5,5) = (in3*(sd2+sd4+4.d0*sd1))*dV2
!
        MP(1,2) = (-4.d0*in3*d1*d2 + insqr3*d4*d0)*dV2
        MP(1,3) = (-2.d0*insqr3*d2*d4)*dV2
        MP(1,4) = (-4.d0*in3*d1*d4 + insqr3*d2*d0)*dV2
        MP(1,5) = (-4.d0*in3*d1*d5)*dV2
!
        MP(2,3) = (insqr3*d1*d4 - d2*sumd3d5)*dV2
        MP(2,4) = (-4.d0*in3*d2*d4 + insqr3*d1*d0)*dV2
        MP(2,5) = (-d1*d4 - insqr3*d2*sumd3d5)*dV2
!
        MP(3,4) = (insqr3*d1*d2 - d4*difd3d5)*dV2
        MP(3,5) = (insqr3*(sd2-sd4))*dV2
!
        MP(4,5) = (d1*d2 + insqr3*d4*difd3d5)*dV2
        do i=1,nd/2
          do j=i,nd/2
            MS(i,j) = (armonicd(i)*armonicd(j))*dV1
            IF (i==j) sumVdel(i,j) = dVdel
          end do
        end do
!
!Construcción de las matrices completas (5x5)
!
        do i=1,nd/2
          do j=i,nd/2
            MP(j,i) = MP(i,j)
            MS(j,i) = MS(i,j)
          end do
        end do
!
! V_delta + (V_sigma - V_delta)*MS + (V_pi - V_delta)*MP
!
        MP = sumVdel + MS + MP
!
!Producto de Kronecker de las matriz anterior (5x5)
!por la matriz unidad 2x2, creando una matriz (10x10)
!
        do i=1,nd/2
          do j=1,nd/2
            MatD(2*i-1,2*j-1) = MP(i,j)
            MatD(2*i,2*j) = MP(i,j)
          end do
        end do
!
!   Ja hem selecionat l'estat previament
!
!        il = 1 + 2*instate
!        do i = 1,nd
!          invar(i) = cmplx(zr(i,il),zi(i,il))
!        end do
        do i = 1,nd
          invars(i) = (0.d0,0.d0)
          do j = 1,nd
            invars(i) = invars(i) + MatD(i,j)*invar(j)
          end do
        end do
        A=0.d0
        do i = 1,nd
          A = A + invars(i)*conjg(invar(i))
        end do
        uimp(ix,iy,iz) = A
        If(Lprint_invar) then
          dum2=den(ix,iy,iz)
          U_ext = U_ext + A * dum2
        Endif
      enddo
    enddo
  enddo
  If(Lprint_invar) then
    write(6,'("ESO", 1p, E15.6)') eso
    U_ext = U_ext*dxyz
    write(6,'("Calculo a partir de U_ext", 1p, E15.6)') U_ext
    U_total = U_ext + eso
    write(6,'("U_total", 1p, E15.6)') U_total
  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Fin Cálculo de Uext_D  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Lfirst=.false.
End if
return
end



subroutine cubic(a,b,c,d,x)

implicit double precision (a-h,o-z)
complex(kind=8) ::  x(3)
real   (kind=8) :: a,b,c,d

pi=4.0D0*DATAN(1.0D0)   

! Step 0: If a is 0 use the quadratic formula. -------------------------
IF(a.eq.0.d0)THEN

   if(b.eq.0.d0)then
       if(c.eq.0.d0)then  ! We have a non-equation
          nroot = 0
       else               ! We have a linear equation with 1 root
          nroot = 1
          x(1) = cmplx(-d/c, 0.)
       endif
   else                   ! We have  quadratic equation
       nroot = 2
       DD = c*c-4.*b*d
       if(DD .ge. 0.d0)then
          x(1) = cmplx((-c+sqrt(DD))/2./b, 0.)
          x(2) = cmplx((-c-sqrt(DD))/2./b, 0.)
       else
          x(1) = cmplx(-c/2./b, +sqrt(-DD)/2./b)
          x(2) = cmplx(-c/2./b, -sqrt(-DD)/2./b)
       endif
   endif

ELSE ! Cubic equation with 3 root  (OJO SOLUCION NO GENERAL!!!!!)
   nroot = 3
   p  = c/a - b*b/a/a/3.d0
   q  = (2.d0*b*b*b/a/a/a - 9.d0*b*c/a/a + 27.d0*d/a) / 27.d0
   DD = p*p*p/27.d0 + q*q/4.d0   ! Calculate the discriminant

   phi    = acos( -(q/dabs(q)) * min( 1.d0, dabs(q)/2.d0/dsqrt(dabs(p*p*p)/27.d0) ) )
   temp1  = 2.d0*dsqrt(dabs(p)/3.d0)

   y1     =  temp1*dcos(phi/3.d0)
   y2     = -temp1*dcos((phi+pi)/3.d0)
   y3     = -temp1*dcos((phi-pi)/3.d0)
   temp1  = b/a/3.d0
   y1     = y1-temp1
   y2     = y2-temp1
   y3     = y3-temp1
   x(1)   = dcmplx( y1,  0.)
   x(2)   = dcmplx( y2,  0.)
   x(3)   = dcmplx( y3,  0.)
ENDIF

return
end


! *************Fin paquete de subrutinas****************************





