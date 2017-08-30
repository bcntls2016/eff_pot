subroutine potenimp(rrimp,iv)
! Esta rutina calcula el potencial sentido por
! el He debido a la impureza. Sabiendo su posicion,
! evalua el potencial de interacion en la malla de trabajo.
! use field
use classicimp !, only: uimp!,Also2!,pairpot
use grid
use energies
use interpol !, only: potdel,potpi,DelInter
use rho
implicit none
integer (kind=4) :: ix,iy,iz,i,j,il,il1, k, kk, l, ll, njz, ki,N_He_exc
integer (kind=4) :: ir, np=6, nd=10
!real    (kind=8) :: zt,yt,r,r2,r4,c1osqr2,Aux
real    (kind=8) :: Aux
real    (kind=8) :: xx,yy,zz
real    (kind=8) , intent(in) :: rrimp(3)
complex (kind=8) , intent(in) :: iv(10)
complex (kind=8)              :: ivs(10)
complex (kind=8)              :: auxEVEN,auxODD
real    (kind=8)              :: rmod
REAL (kind=8) :: sumVdel(10,10), MP(10,10), MS(10,10), MatD(10,10)
REAL (kind=8) :: armonicd(6)
REAL (kind=8) :: sqr3, sqr3o2, insqr3, in3, difd3d5, sumd3d5
REAL (kind=8) :: d0, d1, d2, d3, d4, d5, sd0, sd1, sd2, sd3, sd4, sd5, dV1, dV2, dVDel
Real (Kind=8) :: dum2, U_total, U_Ext
!integer (kind=4), save :: iprint=1
!...................................
!... First, update uimp (for He) ...
!...................................

call updatepoten(rrimp)

If(Lstate.Eq.'P')Then
  do iz=1,nz
   zz = z(iz)-rrimp(3)
   do iy=1,ny
     yy = y(iy)-rrimp(2)
     do ix=1,nx
       xx = x(ix)-rrimp(1)
       auxODD =iv(1)*xx + iv(3)*yy +iv(5)*zz
       auxEVEN=iv(2)*xx + iv(4)*yy +iv(6)*zz
       uimp(ix,iy,iz) = VPi(ix,iy,iz) + Delta(ix,iy,iz)*(auxODD*conjg(auxODD) + auxEVEN*conjg(auxEVEN))
     enddo
   enddo
  enddo
Else
  Als = Als_D
!
! Constantes e inicialización de matrices
!
  sqr3=sqrt(3.d0)
  sqr3o2=sqrt(3.d0)/2.d0
  insqr3=1.d0/sqrt(3.d0)
  in3=1.d0/3.d0
  U_ext = 0.d0
  uimp=0.d0
  do iz=1,nz
    zz=z(iz) - rrimp(3)
    do iy=1,ny
      yy=y(iy) - rrimp(2)
      do ix=1,nx
        xx=x(ix) - rrimp(1)
        MP=0.d0
        MS=0.d0
        MatD=0.d0
        sumVdel=0.d0
!        r=dsqrt(xx**2+yy**2+zz**2)
!        r2 = r*r
!        r4 = r2*r2
!        if(r.eq.0.d0)Cycle
        dVdel= Delta(ix,iy,iz)
        dV1=Pi_Del(ix,iy,iz)         ! Ya esta dividido por r**4
        dV2=Sig_Del(ix,iy,iz)        ! "    "       "    "   "
!
!Armónicos Esféricos
!
        d0=xx*xx+yy*yy+zz*zz
        d1=sqr3*xx*yy
        d2=sqr3*yy*zz
        d3=(-xx*xx-yy*yy+2*zz*zz)/2
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
        do i = 1,nd
          ivs(i) = (0.d0,0.d0)
          do j = 1,nd
            ivs(i) = ivs(i) + MatD(i,j)*iv(j)
          end do
        end do
        Aux=0.d0
        do i = 1,nd
          Aux = Aux + ivs(i)*conjg(iv(i))
        end do
        uimp(ix,iy,iz) = Aux
        If(Lprint_invar) then
          dum2=den(ix,iy,iz)
          U_ext = U_ext + Aux * dum2
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
  Endif
  Lfirst=.false.
!    Write(6,'("From Potenimp(iprint)...",I5)')iprint
!    write(6,'("rimp(1,2,3)..",1p,3e15.6)')rimp
!    write(6,'("vimp(1,2,3)..",1p,3e15.6)')vimp
!    write(6,'("Invar.....",/,1p,(2E15.8))')(invar(i),i=1,ninvar)
!    write(6,'("Hinvar....",/,1p,(2E15.8))')(Hinvar(i),i=1,ninvar)
!    iprint = iprint +1
EndIf

! Capamiento

forall(ix=1:nx,iy=1:ny,iz=1:nz)
 uimp(ix,iy,iz) = min(uimp(ix,iy,iz),1000.d0)
end forall


end subroutine potenimp

!---------------------------------------------------------------------------!

subroutine forceimp(rimp,F)
! Esta rutina calcula el potencial sentido por
! el He debido a la impureza. Sabiendo su posicion,
! evalua el potencial de interacion en la malla de trabajo.
use classicimp , only: uimp, Lfilter_exciplex_force, z_exciplex_exclusion, uHe_He
use deriva
! use field
use grid
use rho
implicit none
integer (kind=4) :: ix,iy,iz, iz0, nz0
real    (kind=8) :: zt,yt,zt2,yt2,r,d
real    (kind=8) :: aux1,aux2,aux3,drV_HeXor,dzV_XTiO,dV
real    (kind=8) , intent(in) :: rimp(3)
real    (kind=8) , intent(out):: F(3)

!.................................!
!... First, the term due to He ...!
!.................................!


F(1) = -sum(dxden*uimp)*dxyz
F(2) = -sum(dyden*uimp)*dxyz
If(Lfilter_exciplex_force)Then
  If(rimp(3).Lt.0.d0)Then
    iz0=(rimp(3) - z(1) + z_exciplex_exclusion)/hz + 1.5
    iz0=max(1,iz0)
    aux1 = 0.d0
    aux2 = 0.d0
    Do iz=iz0,nz
      Do iy=1,ny
        Do ix=1,nx
          aux1 =aux1 + dzden(ix,iy,iz)*uimp(ix,iy,iz)
!          aux2 =aux2 + dzden(ix,iy,iz)*uHe_He(ix,iy,iz)
        EndDo  
      EndDo  
    EndDo
  Else  
    nz0=(rimp(3) - z(1) - z_exciplex_exclusion)/hz + 1.5
    nz0=min(nz,nz0)
    aux1 = 0.d0
    aux2 = 0.d0
    Do iz=1,nz0
      Do iy=1,ny
        Do ix=1,nx
          aux1 =aux1 + dzden(ix,iy,iz)*uimp(ix,iy,iz)
 !         aux2 =aux2 + dzden(ix,iy,iz)*uHe_He(ix,iy,iz)
        EndDo  
      EndDo  
    EndDo
  Endif
  aux1 = -aux1*dxyz
  aux2 = -aux2*dxyz
!  Write(6,'("From Forceimp: dominant term & correction term:",1p,2E15.6)')aux1, aux2
  F(3) = aux1+aux2  
Else        
  F(3) = -sum(dzden*uimp)*dxyz
Endif
!...........................................!
!... Second, the term due to the surface ...!
!...........................................!

!F(3) = F(3) - dzV_XTiO(rimp(3))


end subroutine forceimp


subroutine potinvar(rrimp,iv,Hiv)
!..................................
!!!!!!... This routine gives -i*H.iv      ...
!... This routine gives -i*(H.iv - Ev0.iv)...
!..................................
use grid
! use field
use rho, only: den
use classicimp !, only:Also2, Ev0
use interpol !, only: potdel,potpi,DelInter

implicit none
real    (kind=8) , intent(in) :: rrimp(3)
complex (kind=8) , intent(in) :: iv(10)
complex (kind=8) , intent(out):: Hiv(10)
! complex (kind=8)              :: auxz, auxz1, auxz2, auxz3
! complex (kind=8)              :: cxx,cxy,cxz,cyy,cyz,czz
real    (kind=8)              :: uxx,uxy,uxz,uyy,uyz,uzz
real    (kind=8)              :: xx,yy,zz,yt,zt,r,V_pi,V_sig,aux1, aux2, aux3
integer (kind=4) :: ix,iy,iz,ir, i, j, np=6, nd=10
complex (kind=8)              :: uim = (0.d0,1.d0)
real    (kind=8)              :: rmod
complex (kind=8)              :: MatDSO(10,10)
REAL (kind=8) :: sumVdel(10,10), MP(10,10), MS(10,10), MatD(10,10)
REAL (kind=8) :: armonicd(6)
REAL (kind=8) :: sqr3, sqr3o2, insqr3, in3, difd3d5, sumd3d5
REAL (kind=8) :: d0, d1, d2, d3, d4, d5, sd0, sd1, sd2, sd3, sd4, sd5, dV1, dV2, dVDel
Real (Kind=8) :: dum2, eso, U_total, U_Ext, aux

If(.Not.Lfirst)Then
  Lprint_invar = .false.
Else
  Lprint_invar = .true.
Endif
uxx = 0.d0 ; uxy = 0.d0 ; uxz = 0.d0
uyy = 0.d0 ; uyz = 0.d0 ; uzz = 0.d0

If(Lstate.Eq.'P')Then
!
! Para los estados P
!
  do iz=1,nz
    zz = z(iz)-rrimp(3)
   do iy=1,ny
     yy = y(iy)-rrimp(2)
     do ix=1,nx
       xx = x(ix)-rrimp(1)
       aux1 = den(ix,iy,iz)
       aux2 = delta(ix,iy,iz)
       aux3 = VPi(ix,iy,iz)
       uxx = uxx + aux1 * (aux2*xx*xx + aux3)
       uxy = uxy + aux1 *  aux2*xx*yy
       uxz = uxz + aux1 *  aux2*xx*zz
       uyy = uyy + aux1 * (aux2*yy*yy + aux3)
       uyz = uyz + aux1 *  aux2*yy*zz
       uzz = uzz + aux1 * (aux2*zz*zz + aux3)
     enddo
   enddo
  enddo

  Also2=Als_P/2.0d0
  Hiv(1) =  (uxx*iv(1) + uxy*iv(3) + uxz*iv(5))*dxyz    + Also2*( -uim*iv(3) +     iv(6) )
  Hiv(2) =  (uxx*iv(2) + uxy*iv(4) + uxz*iv(6))*dxyz    + Also2*(  uim*iv(4) -     iv(5) )
  Hiv(3) =  (uxy*iv(1) + uyy*iv(3) + uyz*iv(5))*dxyz    + Also2*(  uim*iv(1) - uim*iv(6) )
  Hiv(4) =  (uxy*iv(2) + uyy*iv(4) + uyz*iv(6))*dxyz    + Also2*( -uim*iv(2) - uim*iv(5) )
  Hiv(5) =  (uxz*iv(1) + uyz*iv(3) + uzz*iv(5))*dxyz    + Also2*(     -iv(2) + uim*iv(4) )
  Hiv(6) =  (uxz*iv(2) + uyz*iv(4) + uzz*iv(6))*dxyz    + Also2*(      iv(1) + uim*iv(3) )
!
! Hiv = int(rho*(Vs-Vp)*(x_j*iv^j)*x_i) + VSO_{ij}*iv^j
!
Else
!
! Para los estados D
!
  Als   = Als_D
  Also2 = Als*0.5d0
!
! Constantes e inicialización de matrices
!
  sqr3=sqrt(3.d0)
  sqr3o2=sqrt(3.d0)/2.d0
  insqr3=1.d0/sqrt(3.d0)
  in3=1.d0/3.d0
        MP=0.d0        !
        MS=0.d0        !  Estas matrices las utilizamos para calcular
        MatD=0.d0      !  el hamiltoniano de | \lambda >, por tanto
        sumVdel=0.d0   !  hemos de integrar lo potenciales con la densidad
  do iz=1,nz
    zz=z(iz) - rrimp(3)
    do iy=1,ny
      yy=y(iy) - rrimp(2)
      do ix=1,nx
        xx=x(ix) - rrimp(1)
!        r=dsqrt(xx**2+yy**2+zz**2) !  Las r4 estan dentro de los potencialse
!        if(r.eq.0.d0)r=1.d-8       !
        Aux = den(ix,iy,iz)
        dVdel= Delta(ix,iy,iz)*Aux
        dV1=Pi_Del(ix,iy,iz)*Aux
        dV2=Sig_Del(ix,iy,iz)*Aux
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
        MP(1,1) = MP(1,1) + (in3*(sd2+sd4+4.d0*sd5))*dV2
        MP(2,2) = MP(2,2) + (in3*(sd1+sd4) + sumd3d5*sumd3d5)*dV2
        MP(3,3) = MP(3,3) + (sd4 + sd2)*dV2
        MP(4,4) = MP(4,4) + (in3*(sd1+sd2) + difd3d5*difd3d5)*dV2
        MP(5,5) = MP(5,5) +  (in3*(sd2+sd4+4.d0*sd1))*dV2
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
      enddo
    enddo
  enddo
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
!  Multiplicamos por el elemento de volumen para acabar la integracion

  MP = (sumVdel + MS + MP)*dxyz
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
!Suma matriz MatD con la de Spin-órbita (10x10) 
!
  Also2=Als_D/2.0d0
  MatDSO = MatD + Also2*SOD
  do i = 1,nd
    Hiv(i) = (0.d0,0.d0)
    do j = 1,nd
      Hiv(i) = Hiv(i) + MatDSO(i,j)*iv(j)
     end do
  end do
Endif
If(Lprint_invar) then
  Aux1=0.d0
  do i = 1,nd
    Aux1 = Aux1 + Hiv(i)*conjg(iv(i))
  end do
  Write(6,'("From potinvar: Total Impurity Energy...",1p,E15.6)')Aux1
Endif

!
! Hiv = {int(rho*[Vd*\delta{i,j}+(Vp-Vd)*Mp{i,j}+(Vs-Vd)*Ms{i,j}]) + VSO_{ij}}*iv^j
!

Hiv = -uim*(Hiv-Ev0*iv)    ! Aqui restem una constant per restar una fase global

!Lfirst = .false.

end subroutine potinvar

double precision function V_HeX(r)
implicit none
real (kind=8) :: r,ff, V_gs
V_HeX = V_gs(r)
end function

double precision function V_Sig(r)
implicit none
real (kind=8) :: r,ff,V_sigma
V_sig = V_sigma(r)

end function

double precision function dzV_XTiO(z)
implicit none
real (kind=4) :: z

dzV_XTiO = 0.d0

end function

subroutine updatepoten(rrimp)
use classicimp !, only: uimp!,Also2!,pairpot
use seleccio_de_potencial
use grid
use interpol
implicit none
real    (kind=8) , intent(in) :: rrimp(3)
real    (kind=8)              :: dist(3)
real    (kind=8)              :: r,yt,zt,rmod,update,rmax,r2,r4
!real    (kind=8) , parameter  :: tolerance=1.d-40
real    (kind=8) , parameter  :: zeror=1.d-3
integer (kind=8)              :: ix,iy,iz,ir
!Logical :: Lparab =.false.
Logical :: Lparab =.true.

rmax=(npot-1)*DelInter
IF(Lstate.Eq.'P')Then
  do iz=1,nz
    zt = (z(iz)-rrimp(3))**2
    do iy=1,ny
      yt = (y(iy)-rrimp(2))**2 + zt
      do ix=1,nx
        r2 = (x(ix)-rrimp(1))**2 + yt
        r  = dsqrt(r2)
        Vpi(ix,iy,iz)=update(r,Lparab,Npot,rmax,potpi)
        If(r.Gt.zeror)Then
          Delta(ix,iy,iz)=update(r,Lparab,Npot,rmax,potdel)/r2
        Else
          Delta(ix,iy,iz)=umax_sigma
        Endif
      enddo
    enddo
  enddo
Else
  do iz=1,nz
    zt = (z(iz)-rrimp(3))**2
    do iy=1,ny
      yt = (y(iy)-rrimp(2))**2 + zt
      do ix=1,nx
        r2 = (x(ix)-rrimp(1))**2 + yt
        r4 = r2*r2
        r  = dsqrt(r2)
        Delta(ix,iy,iz) = update(r,Lparab,Npot,rmax,potdel)
        If(r.Gt.zeror)Then
          Pi_Del(ix,iy,iz)  = update(r,Lparab,Npot,rmax,potpidel)/r4
          Sig_Del(ix,iy,iz) = update(r,Lparab,Npot,rmax,potsigdel)/r4
        Else
          Pi_Del(ix,iy,iz)  = umax_pi
          Sig_Del(ix,iy,iz) = umax_sigma
!          Pi_Del(ix,iy,iz)  = 0.d0
!          Sig_Del(ix,iy,iz) = 0.d0
        Endif
      enddo
    enddo
  enddo
Endif
If(Lfilter_exciplex_force)Then
  do iz=1,nz
    zt = (z(iz) - rrimp(3) - z_exciplex_position)**2
    do iy=1,ny
      yt = (y(iy)-rrimp(2))**2 + zt
      do ix=1,nx
        r2 = (x(ix)-rrimp(1))**2 + yt
        r  = dsqrt(r2)
        uHe_He(ix,iy,iz) = update(r,Lparab,Npot,rmax,pothe_he)
      enddo
    enddo
  enddo
Endif

end subroutine updatepoten


Double Precision Function update(r,Lparab,Nint,Rmax,f)

implicit none
integer (kind=4)   :: i0,i1,i2
real    (kind=8)   :: h,f0,f1,f2,r0,r1,r2,aa,bb,cc,c1oh,c1o2h2
Logical , intent(in) ::  Lparab
integer (kind=4) , intent(in) :: nint
real    (kind=8) , intent(in) :: r,rmax,f(0:Nint+1)

    h=rmax/(Nint-1)
    c1oh=1.0d0/h
    c1o2h2=c1oh*c1oh*0.5d0

    i0 = 1.5 + r/h
    r0 = (i0-1)*h
    f0 = f(i0)
    If(r.Lt.r0)Then
      i1 = i0 - 1
      i2 = i0 + 1
    Else
      i1 = i0 + 1
      i2 = i0 - 1
    Endif
    f1 = f(i1)
    If(Lparab)Then
!
! Interpolació parabòlica
!
      f2 = f(i2)
      aa = (f2-2.0d0*f0+f1)*c1o2h2
      If(r.Lt.r0)Then
        bb = (f2-f0)*c1oh - aa*(2.0d0*r0+h)
      Else
        bb = (f1-f0)*c1oh - aa*(2.0d0*r0+h)
      Endif
      cc = f0 - r0*(bb + r0*aa)
      Update = r*(aa*r+bb)+cc
    Else
!
! Interpolació lineal
!
      Update = f0 - (f0-f1)*c1oh*Abs(r-r0)
    Endif
Return
End


