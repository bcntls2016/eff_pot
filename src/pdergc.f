      SUBROUTINE PDERGC(K,M,NN,H,F,DKF,NV,IX,NMX,IW,ICON)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 F,DKF,SUMA,SUMB
      LOGICAL COND
      COMMON/SDERA/ A(25792)
      DIMENSION F(*),DKF(*),NN(*),NMX(*),IW(*)
!
!       K:  Ordre de la derivada
!       M:  Numero de punts farem servir per trovar la derivada
!      NN:  Vector on posarem las dimensions reals de les variables
!           1 fins a NV. Diemensio >= NV
!       H:  Pas per la variables IX
!       F:  Matriu de la funcio multivariable
!     DKF:  Matriu de la deivada parcial de F respecte la variable IX
!      NV:  Numero de variables de la funcio
!      IX:  Index corresponent a la variable de la que volem calcular
!           la derivada parcial
!     NMX:  Vector on posarem les dimensions maximas dels NV-1 primers
!           index, dimensio >= 2*NV-2:
!                 1  fins a   NV-1 per la matriu F
!                 NV fins a 2*NV-2 per la matriu DKF
!      IW:  Vector auxiliar de dimensio >=4*NV-1
!
!     Funciones especiales:
!      ICON:
!         0 Calcula la derivada deseada
!
!         2 Es com si afagisim un punt mes de la funcio a la
!           dreta del ultim punt del grid, amb valor 0.
!
!         3 El mateix que per ICON=2 i a fa la mateixa hipotesi
!           a l'esquerra del primer punt del grid.
!
!         4 Fa una reflexio parell de la funcio per els primers
!           punts, per tal de tenir els punts de la derivada
!           centrats, per aquests primers punts
!
!         5 El mateix que per ICON=4 pero amb reflexio sanar
!
!         6 El mateix que ICON=4 i ICON=2
!
!         7 El mateix que ICON=5 i ICON=2
!
!         8 Condicions periodicas: ( I<1 <--> I+N, I>N <--> I-N )
!
!        12 Soposem la funcio nulla a la dreta del ultim punt
!           del grid per tal que sempre fem servir els coeficients
!           centrats
!
!        13 El mateix que per ICON=12 i fa la mateixa hipotesi
!           a l'esquerra del primer punt del grid.
!
!        16 El mateix que ICON=4 i ICON=12
!
!        17 El mateix que ICON=5 i ICON=12
!
!        20 Condicions periodiques de Wigner_Seitz ( N+i <--> N-i )
!
!        24 El mateix que ICON= 4 i ICON=20
!
!        25 El mateix que ICON= 5 i ICON=20
!
!
!     Restriccions:
!       M:  Te que esser sanas, mes gran o igual a 3 
!           y mes petit o igual que 25
!       K:  Pot valdre desde 1 fins M-1
!  NN(IX):  Te que esser mes gran o igual que M
!
!
      N=NN(IX)
      IF(M.LT.3.OR.M.GT.25.OR.M.GT.N.OR.MOD(M,2).EQ.0
     &  .OR.K.LT.0.OR.K.GT.M-1)THEN
        WRITE(*,1100)N,M,K
        RETURN
      ENDIF
      SIG=1.0D0-2.0D0*MOD(K,2)
      MM1=M-1
      COND=K.LT.MM1
      HK=H**K
      JI=1
      M2=M*M
      MM1O2=MM1/2
      M21O2=(M2+1)/2
      FACT=1.0D0/HK
      DO I=3,M-2,2
        JI=JI+(I*I+1)*(I-2)/2+I
      ENDDO
      JI=JI+M21O2*(K-1)
      IF(COND)THEN
        JF=JI+M21O2-1
      ELSE
        JF=JI+MM1
      ENDIF
      N0=MM1O2+1
      NF=N-MM1O2
      NRV=1
      NVM1=NV-1
      NVM2=NV-2
      I0=0
      IND0=NV-1
      IPF0=2*NV-1
      IPD0=3*NV-1
      IPF1=1+IPF0
      IPD1=1+IPD0
      IW(IPF1)=1
      IW(IPD1)=1
      DO IV=1,NVM1
        IF(IV.EQ.IX)I0=1
        I1=IV+I0
        IV1=IV+1
        NRV=NRV*NN(I1)
        IW(IV)=I1
        IND1=IND0+I1
        IW(IND1)=1
        IPFV=IV+IPF0
        IPDV=IV+IPD0
        IPFV1=IV1+IPF0
        IPDV1=IV1+IPD0
        IW(IPFV1)=IW(IPFV)*NMX(IV)
        IW(IPDV1)=IW(IPDV)*NMX(IV+NVM1)
      ENDDO
      IOV=IW(1)
      INDV=IOV+IND0
      IW(INDV)=0
      DO IT=1,NRV
        IJKF0=1
        IJKD0=1
        IOV=IW(1)
        INDV=IOV+IND0
        IW(INDV)=IW(INDV)+1
        DO IV=1,NVM2
          IOV=IW(IV)
          II=NN(IOV)
          INDV=IOV+IND0
          JJ=IW(INDV)
          IF(JJ.GT.II)THEN
            IV1=IV+1
            IW(INDV)=1
            IOV1=IW(IV1)
            INDV1=IOV1+IND0
            IW(INDV1)=IW(INDV1)+1
          ENDIF
          INDV=IOV+IND0
          IPFV=IOV+IPF0
          IPDV=IOV+IPD0
          IJKF0=IJKF0+(IW(INDV)-1)*IW(IPFV)
          IJKD0=IJKD0+(IW(INDV)-1)*IW(IPDV)
        ENDDO
        INDV=IW(NVM1)+IND0
        IPFV=IW(NVM1)+IPF0
        IPDV=IW(NVM1)+IPD0
        IJKF0=IJKF0+(IW(INDV)-1)*IW(IPFV)
        IJKD0=IJKD0+(IW(INDV)-1)*IW(IPDV)
        IPFX=IX+IPF0
        IPDX=IX+IPD0
        IF(ICON.EQ.0)THEN
!
!          Calcula la derivada desitjada
!
          IF(COND)THEN
            DO I=1,MM1O2
              L0=JI+(I-1)*M
              L1=L0+MM1
              SUMA=0.0D0
              L=L0-1
              DO J=1,M
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            ENDDO
            L0=JI+MM1O2*M
            L1=L0+MM1O2
            DO I=N0,NF
              SUMA=0.0D0
              L=L0-1
              DO J=I-MM1O2,I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              L=L1
              SUMB=0.0D0
              DO J=I+1,MM1O2+I
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=(SUMA+SUMB*SIG)*FACT
            ENDDO
            DO II=MM1O2,1,-1
              L1=JI+(II-1)*M
              L0=L1+MM1
              SUMA=0.0D0
              L=L0+1
              I=NF+MM1O2-II+1
              DO J=NF-MM1O2,N
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*SIG*FACT
            ENDDO
          ELSE
            DO I=N0,NF
              SUMA=0.0D0
              L=JI-1
              DO J=I-MM1O2,I+MM1O2
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            ENDDO
            IJKN0=IJKD0+(N0-1)*IW(IPDX)
            IJKNF=IJKD0+(NF-1)*IW(IPDX)
            DO II=1,MM1O2
              I=N0-II
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=DKF(IJKN0)
              I=NF+II
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=DKF(IJKNF)
            ENDDO
          ENDIF
        ELSE IF(ICON.EQ.2)THEN
!
!          Es com si afagisim un punt mes de la funcio a la
!          dreta del ultim punt del grid, amb valor 0.
!
          IF(COND)THEN
            DO I=1,MM1O2
              L0=JI+(I-1)*M
              L1=L0+MM1
              SUMA=0.0D0
              L=L0-1
              DO J=1,M
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            ENDDO
            L0=JI+MM1O2*M
            L1=L0+MM1O2
            DO I=N0,NF
              SUMA=0.0D0
              L=L0-1
              DO J=I-MM1O2,I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              L=L1
              SUMB=0.0D0
              DO J=I+1,MM1O2+I
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=(SUMA+SUMB*SIG)*FACT
            ENDDO
              I=NF+1
              SUMA=0.0D0
              L=L0-1
              DO J=I-MM1O2,I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              L=L1
              SUMB=0.0D0
              DO J=I+1,MM1O2+I-1
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=(SUMA+SUMB*SIG)*FACT
            DO II=MM1O2-1,1,-1
              L1=JI+II*M
              L0=L1+MM1
              SUMA=0.0D0
              L=L0+1
              I=NF+MM1O2-II+1
              DO J=NF-MM1O2+1,N
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*SIG*FACT
            ENDDO
          ELSE
            DO I=N0,NF
              SUMA=0.0D0
              L=JI-1
              DO J=I-MM1O2,I+MM1O2
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            ENDDO
            IJKN0=IJKD0+(N0-1)*IW(IPDX)
            DO II=1,MM1O2
              I=N0-II
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=DKF(IJKN0)
            ENDDO
              I=NF+1
              SUMA=0.0D0
              L=JI-1
              DO J=I-MM1O2,N
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
              IJKNF=IJK
            DO I=NF+2,N
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=DKF(IJKNF)
            ENDDO
          ENDIF
        ELSE IF(ICON.EQ.3)THEN
!
!          El mateix que per ICON=2 i a fa la mateixa hipotesi
!          a l'esquerra del primer punt del grid.
!
          IF(COND)THEN
            DO I=1,MM1O2-1
              L0=JI+I*M
              L1=L0+MM1
              SUMA=0.0D0
              L=L0
              DO J=1,MM1
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            ENDDO
            L0=JI+MM1O2*M
            L1=L0+MM1O2
              I=MM1O2
              SUMA=0.0D0
              L=L0
              DO J=I-MM1O2+1,I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              L=L1
              SUMB=0.0D0
              DO J=I+1,MM1O2+I
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=(SUMA+SUMB*SIG)*FACT
            DO I=N0,NF
              SUMA=0.0D0
              L=L0-1
              DO J=I-MM1O2,I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              L=L1
              SUMB=0.0D0
              DO J=I+1,MM1O2+I
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=(SUMA+SUMB*SIG)*FACT
            ENDDO
              I=NF+1
              SUMA=0.0D0
              L=L0-1
              DO J=I-MM1O2,I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              L=L1
              SUMB=0.0D0
              DO J=I+1,MM1O2+I-1
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=(SUMA+SUMB*SIG)*FACT
            DO II=MM1O2-1,1,-1
              L1=JI+II*M
              L0=L1+MM1
              SUMA=0.0D0
              L=L0+1
              I=NF+MM1O2-II+1
              DO J=NF-MM1O2+1,N
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*SIG*FACT
            ENDDO
          ELSE
              I=1
              SUMA=0.0D0
              L=JI
              DO J=1,MM1
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            DO I=2,MM1O2
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=DKF(IJKD0)
            ENDDO
            DO I=N0,NF
              SUMA=0.0D0
              L=JI-1
              DO J=I-MM1O2,I+MM1O2
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            ENDDO
              I=NF+1
              SUMA=0.0D0
              L=JI-1
              DO J=I-MM1O2,N
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            IJKNF=IJK
            DO I=NF+2,N
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=DKF(IJKNF)
            ENDDO
          ENDIF
        ELSE IF(ICON.EQ.4)THEN
!
!           Fa una reflexio parell de la funcio per els primers
!           punts, per tal de tenir els punts de la derivada
!           centrats, per aquests primers punts
!
          IF(COND)THEN
            L0=JI+MM1O2*M
            L1=L0+MM1O2
            DO I=1,MM1O2
              SUMA=0.0D0
              L=L0-1
              DO JJ=MM1O2+1-I,1,-1
                J=JJ+1
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              DO J=1,I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              L=L1
              SUMB=0.0D0
              DO J=I+1,MM1O2+I
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=(SUMA+SUMB*SIG)*FACT
            ENDDO
            DO I=N0,NF
              SUMA=0.0D0
              L=L0-1
              DO J=I-MM1O2,I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              L=L1
              SUMB=0.0D0
              DO J=I+1,MM1O2+I
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=(SUMA+SUMB*SIG)*FACT
            ENDDO
            DO II=MM1O2,1,-1
              L1=JI+(II-1)*M
              L0=L1+MM1
              SUMA=0.0D0
              L=L0+1
              I=NF+MM1O2-II+1
              DO J=NF-MM1O2,N
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*SIG*FACT
            ENDDO
          ELSE
            DO I=1,MM1O2
              L=JI-1
              SUMA=0.0D0
              DO JJ=MM1O2+1-I,1,-1
                J=JJ+1
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              DO J=1,MM1O2+I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            ENDDO
            DO I=N0,NF
              SUMA=0.0D0
              L=JI-1
              DO J=I-MM1O2,I+MM1O2
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            ENDDO
            IJKNF=IJKD0+(NF-1)*IW(IPDX)
            DO II=1,MM1O2
              I=NF+II
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=DKF(IJKNF)
            ENDDO
          ENDIF
        ELSE IF(ICON.EQ.5)THEN
!
!          Fa una reflexio sanar de la funcio per els primers
!          punts, per tal de tenir els punts de la derivada
!          sempre centrats.
!
          IF(COND)THEN
            L0=JI+MM1O2*M
            L1=L0+MM1O2
            DO I=1,MM1O2
              SUMA=0.0D0
              L=L0-1
              DO JJ=MM1O2+1-I,1,-1
                J=JJ+1
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA-A(L)*F(IJK)
              ENDDO
              DO J=1,I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              SUMB=0.0D0
              L=L1
              DO J=I+1,MM1O2+I
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=(SUMA+SUMB*SIG)*FACT
            ENDDO
            DO I=N0,NF
              SUMA=0.0D0
              L=L0-1
              DO J=I-MM1O2,I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              SUMB=0.0D0
              L=L1
              DO J=I+1,MM1O2+I
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=(SUMA+SUMB*SIG)*FACT
            ENDDO
            DO II=MM1O2,1,-1
              L1=JI+(II-1)*M
              L0=L1+MM1
              SUMA=0.0D0
              L=L0+1
              I=NF+MM1O2-II+1
              DO J=NF-MM1O2,N
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*SIG*FACT
            ENDDO
          ELSE
            DO I=1,MM1O2
              L=JI-1
              SUMA=0.0D0
              DO JJ=MM1O2+1-I,1,-1
                J=JJ+1
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA-A(L)*F(IJK)
              ENDDO
              DO J=1,MM1O2+I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            ENDDO
            DO I=N0,NF
              SUMA=0.0D0
              L=JI-1
              DO J=I-MM1O2,I+MM1O2
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            ENDDO
            IJKNF=IJKD0+(NF-1)*IW(IPDX)
            DO II=1,MM1O2
              I=NF+II
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=DKF(IJKNF)
            ENDDO
          ENDIF
        ELSE IF(ICON.EQ.6)THEN
!
!          El mateix que per ICON=4 i ICON=2
!
          IF(COND)THEN
            L0=JI+MM1O2*M
            L1=L0+MM1O2
            DO I=1,MM1O2
              SUMA=0.0D0
              L=L0-1
              DO JJ=MM1O2+1-I,1,-1
                J=JJ+1
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              DO J=1,I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              L=L1
              SUMB=0.0D0
              DO J=I+1,MM1O2+I
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=(SUMA+SUMB*SIG)*FACT
            ENDDO
            DO I=N0,NF
              SUMA=0.0D0
              L=L0-1
              DO J=I-MM1O2,I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              L=L1
              SUMB=0.0D0
              DO J=I+1,MM1O2+I
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=(SUMA+SUMB*SIG)*FACT
            ENDDO
              I=NF+1
              SUMA=0.0D0
              L=L0-1
              DO J=I-MM1O2,I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              L=L1
              SUMB=0.0D0
              DO J=I+1,MM1O2+I-1
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=(SUMA+SUMB*SIG)*FACT
            DO II=MM1O2-1,1,-1
              L1=JI+II*M
              L0=L1+MM1
              SUMA=0.0D0
              L=L0+1
              I=NF+MM1O2-II+1
              DO J=NF-MM1O2+1,N
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*SIG*FACT
            ENDDO
          ELSE
            DO I=1,MM1O2
              L=JI-1
              SUMA=0.0D0
              DO JJ=MM1O2+1-I,1,-1
                J=JJ+1
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              DO J=1,MM1O2+I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            ENDDO
            DO I=N0,NF
              SUMA=0.0D0
              L=JI-1
              DO J=I-MM1O2,I+MM1O2
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            ENDDO
              I=NF+1
              SUMA=0.0D0
              L=JI-1
              DO J=I-MM1O2,N
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            IJKNF=IJK
            DO I=NF+2,N
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=DKF(IJKNF)
            ENDDO
          ENDIF
        ELSE IF(ICON.EQ.7)THEN
!
!          El mateix que per ICON=5 i ICON=2
!
          IF(COND)THEN
            L0=JI+MM1O2*M
            L1=L0+MM1O2
            DO I=1,MM1O2
              SUMA=0.0D0
              L=L0-1
              DO JJ=MM1O2+1-I,1,-1
                J=JJ+1
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA-A(L)*F(IJK)
              ENDDO
              DO J=1,I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              SUMB=0.0D0
              L=L1
              DO J=I+1,MM1O2+I
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=(SUMA+SUMB*SIG)*FACT
            ENDDO
            DO I=N0,NF
              SUMA=0.0D0
              L=L0-1
              DO J=I-MM1O2,I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              SUMB=0.0D0
              L=L1
              DO J=I+1,MM1O2+I
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=(SUMA+SUMB*SIG)*FACT
            ENDDO
              I=NF+1
              SUMA=0.0D0
              L=L0-1
              DO J=I-MM1O2,I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              L=L1
              SUMB=0.0D0
              DO J=I+1,MM1O2+I-1
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=(SUMA+SUMB*SIG)*FACT
            DO II=MM1O2-1,1,-1
              L1=JI+II*M
              L0=L1+MM1
              SUMA=0.0D0
              L=L0+1
              I=NF+MM1O2-II+1
              DO J=NF-MM1O2+1,N
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*SIG*FACT
            ENDDO
          ELSE
            DO I=1,MM1O2
              L=JI-1
              SUMA=0.0D0
              DO JJ=MM1O2+1-I,1,-1
                J=JJ+1
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA-A(L)*F(IJK)
              ENDDO
              DO J=1,MM1O2+I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            ENDDO
            DO I=N0,NF
              SUMA=0.0D0
              L=JI-1
              DO J=I-MM1O2,I+MM1O2
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            ENDDO
              I=NF+1
              SUMA=0.0D0
              L=JI-1
              DO J=I-MM1O2,N
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            IJKNF=IJK
            DO I=NF+2,N
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=DKF(IJKNF)
            ENDDO
          ENDIF
        ELSE IF(ICON.EQ.8)THEN
!
!          Condicions periodicas: ( I<1 <--> I+N, I>N <--> I-N )
!
          IF(COND)THEN
            L0=JI+MM1O2*M
            L1=L0+MM1O2
            DO I=1,MM1O2
              SUMA=0.0D0
              L=L0-1
              DO JJ=I-MM1O2,0
                L=L+1
                J=JJ+N
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              DO J=1,I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              L=L1
              SUMB=0.0D0
              DO J=I+1,MM1O2+I
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=(SUMA+SUMB*SIG)*FACT
            ENDDO
            DO I=N0,NF
              SUMA=0.0D0
              L=L0-1
              DO J=I-MM1O2,I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              L=L1
              SUMB=0.0D0
              DO J=I+1,MM1O2+I
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
               SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=(SUMA+SUMB*SIG)*FACT
            ENDDO
            DO I=NF+1,N
              SUMA=0.0D0
              L=L0-1
              DO J=I-MM1O2,I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              L=L1
              SUMB=0.0D0
              DO J=I+1,N
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              DO JJ=N+1,MM1O2+I
                L=L-1
                J=JJ-N
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=(SUMA+SUMB*SIG)*FACT
            ENDDO
          ELSE
            DO I=1,MM1O2
              SUMA=0.0D0
              L=JI-1
              DO JJ=I-MM1O2,0
                L=L+1
                J=JJ+N
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              DO J=1,I+MM1O2
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            ENDDO
            DO I=N0,NF
              SUMA=0.0D0
              L=JI-1
              DO J=I-MM1O2,I+MM1O2
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            ENDDO
            DO I=NF+1,N
              SUMA=0.0D0
              L=JI-1
              DO J=I-MM1O2,N
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              DO JJ=N+1,I+MM1O2
                L=L+1
                J=JJ-N
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            ENDDO
          ENDIF
        ELSE IF(ICON.EQ.12)THEN
!
!           Soposem la funcio nulla a la dreta del ultim punt
!           del grid per tal que sempre fem servir els coeficients
!           centrats
!
          IF(COND)THEN
            DO I=1,MM1O2
              L0=JI+(I-1)*M
              L1=L0+MM1
              SUMA=0.0D0
              L=L0-1
              DO J=1,M
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            ENDDO
            L0=JI+MM1O2*M
            L1=L0+MM1O2
            DO I=N0,NF
              SUMA=0.0D0
              L=L0-1
              DO J=I-MM1O2,I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              SUMB=0.0D0
              L=L1
              DO J=I+1,MM1O2+I
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=(SUMA+SUMB*SIG)*FACT
            ENDDO
            DO I=NF+1,N
              SUMA=0.0D0
              L=L0-1
              DO J=I-MM1O2,I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              SUMB=0.0D0
              L=L1
              DO J=I+1,N
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=(SUMA+SUMB*SIG)*FACT
            ENDDO
          ELSE
            DO I=N0,NF
              SUMA=0.0D0
              L=JI-1
              DO J=I-MM1O2,I+MM1O2
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            ENDDO
            IJKN0=IJKD0+(N0-1)*IW(IPDX)
            DO II=1,MM1O2
              I=N0-II
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=DKF(IJKN0)
            ENDDO
            DO I=NF+1,N
              SUMA=0.0D0
              L=JI-1
              DO J=I-MM1O2,N
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            ENDDO
          ENDIF
        ELSE IF(ICON.EQ.13)THEN
!
!         El mateix que per ICON=12 i fa la mateixa hipotesi
!         a l'esquerra del primer punt del grid.
!
          IF(COND)THEN
            L0=JI+MM1O2*M
            L1=L0+MM1O2
            DO I=1,MM1O2
              SUMA=0.0D0
              L=L0-I+MM1O2
              DO J=1,I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              SUMB=0.0D0
              L=L1
              DO J=I+1,MM1O2+I
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=(SUMA+SUMB*SIG)*FACT
            ENDDO
            DO I=N0,NF
              SUMA=0.0D0
              L=L0-1
              DO J=I-MM1O2,I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              SUMB=0.0D0
              L=L1
              DO J=I+1,MM1O2+I
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=(SUMA+SUMB*SIG)*FACT
            ENDDO
            DO I=NF+1,N
              SUMA=0.0D0
              L=L0-1
              DO J=I-MM1O2,I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              SUMB=0.0D0
              L=L1
              DO J=I+1,N
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=(SUMA+SUMB*SIG)*FACT
            ENDDO
          ELSE
            DO I=1,MM1O2
              SUMA=0.0D0
              L=JI-I+MM1O2
              DO J=1,I+MM1O2
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            ENDDO
            DO I=N0,NF
              SUMA=0.0D0
              L=JI-1
              DO J=I-MM1O2,I+MM1O2
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            ENDDO
            DO I=NF+1,N
              SUMA=0.0D0
              L=JI-1
              DO J=I-MM1O2,N
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            ENDDO
          ENDIF
        ELSE IF(ICON.EQ.16)THEN
!
!          El mateix que ICON= 4 i ICON=12
!
          IF(COND)THEN
            L0=JI+MM1O2*M
            L1=L0+MM1O2
            DO I=1,MM1O2
              SUMA=0.0D0
              L=L0-1
              DO JJ=MM1O2+1-I,1,-1
                J=JJ+1
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              DO J=1,I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              SUMB=0.0D0
              L=L1
              DO J=I+1,MM1O2+I
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=(SUMA+SUMB*SIG)*FACT
            ENDDO
            DO I=N0,NF
              SUMA=0.0D0
              L=L0-1
              DO J=I-MM1O2,I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              SUMB=0.0D0
              L=L1
              DO J=I+1,MM1O2+I
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=(SUMA+SUMB*SIG)*FACT
            ENDDO
            DO I=NF+1,N
              SUMA=0.0D0
              L=L0-1
              DO J=I-MM1O2,I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              SUMB=0.0D0
              L=L1
              DO J=I+1,N
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=(SUMA+SUMB*SIG)*FACT
            ENDDO
          ELSE
            DO I=1,MM1O2
              L=JI-1
              SUMA=0.0D0
              DO JJ=MM1O2+1-I,1,-1
                J=JJ+1
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              DO J=1,MM1O2+I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            ENDDO
            DO I=N0,NF
              SUMA=0.0D0
              L=JI-1
              DO J=I-MM1O2,I+MM1O2
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            ENDDO
            DO I=NF+1,N
              SUMA=0.0D0
              L=JI-1
              DO J=I-MM1O2,N
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            ENDDO
          ENDIF
        ELSE IF(ICON.EQ.17)THEN
!
!          El mateix que ICON= 5 i ICON=12
!
!
          IF(COND)THEN
            L0=JI+MM1O2*M
            L1=L0+MM1O2
            DO I=1,MM1O2
              SUMA=0.0D0
              L=L0-1
              DO JJ=MM1O2+1-I,1,-1
                J=JJ+1
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA-A(L)*F(IJK)
              ENDDO
              DO J=1,I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              SUMB=0.0D0
              L=L1
              DO J=I+1,MM1O2+I
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=(SUMA+SUMB*SIG)*FACT
            ENDDO
            DO I=N0,NF
              SUMA=0.0D0
              L=L0-1
              DO J=I-MM1O2,I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              SUMB=0.0D0
              L=L1
              DO J=I+1,MM1O2+I
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=(SUMA+SUMB*SIG)*FACT
            ENDDO
            DO I=NF+1,N
              SUMA=0.0D0
              L=L0-1
              DO J=I-MM1O2,I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              SUMB=0.0D0
              L=L1
              DO J=I+1,N
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=(SUMA+SUMB*SIG)*FACT
            ENDDO
          ELSE
            DO I=1,MM1O2
              L=JI-1
              SUMA=0.0D0
              DO JJ=MM1O2+1-I,1,-1
                J=JJ+1
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA-A(L)*F(IJK)
              ENDDO
              DO J=1,MM1O2+I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            ENDDO
            DO I=N0,NF
              SUMA=0.0D0
              L=JI-1
              DO J=I-MM1O2,I+MM1O2
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            ENDDO
            DO I=NF+1,N
              SUMA=0.0D0
              L=JI-1
              DO J=I-MM1O2,N
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            ENDDO
          ENDIF
        ELSE IF(ICON.EQ.20)THEN
!
!        Condicions periodiques de Wigner_Seitz ( N+i <--> N-i )
!
          IF(COND)THEN
            DO I=1,MM1O2
              L0=JI+(I-1)*M
              L1=L0+MM1
              SUMA=0.0D0
              L=L0-1
              DO J=1,M
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            ENDDO
            L0=JI+MM1O2*M
            L1=L0+MM1O2
            DO I=N0,NF
              SUMA=0.0D0
              L=L0-1
              DO J=I-MM1O2,I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              L=L1
              SUMB=0.0D0
              DO J=I+1,MM1O2+I
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=(SUMA+SUMB*SIG)*FACT
            ENDDO
            DO I=NF+1,N
              SUMA=0.0D0
              L=L0-1
              DO J=I-MM1O2,I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              L=L1
              SUMB=0.0D0
              DO J=I+1,N
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              DO JJ=N+1,MM1O2+I
                L=L-1
                J=2*N-JJ
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=(SUMA+SUMB*SIG)*FACT
            ENDDO
          ELSE
            DO I=N0,NF
              SUMA=0.0D0
              L=JI-1
              DO J=I-MM1O2,I+MM1O2
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            ENDDO
            IJKN0=IJKD0+(N0-1)*IW(IPDX)
            DO II=1,MM1O2
              I=N0-II
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=DKF(IJKN0)
            ENDDO
            DO I=NF+1,N
              SUMA=0.0D0
              L=JI-1
              DO J=I-MM1O2,N
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              DO JJ=N+1,I+MM1O2
                L=L+1
                J=2*N-JJ
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            ENDDO
          ENDIF
        ELSE IF(ICON.EQ.24)THEN
!
!        El mateix que ICON=4 i ICON=20
!
          IF(COND)THEN
            L0=JI+MM1O2*M
            L1=L0+MM1O2
            DO I=1,MM1O2
              SUMA=0.0D0
              L=L0-1
              DO JJ=MM1O2+1-I,1,-1
                J=JJ+1
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              DO J=1,I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              L=L1
              SUMB=0.0D0
              DO J=I+1,MM1O2+I
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=(SUMA+SUMB*SIG)*FACT
            ENDDO
            DO I=N0,NF
              SUMA=0.0D0
              L=L0-1
              DO J=I-MM1O2,I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              L=L1
              SUMB=0.0D0
              DO J=I+1,MM1O2+I
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=(SUMA+SUMB*SIG)*FACT
            ENDDO
            DO I=NF+1,N
              SUMA=0.0D0
              L=L0-1
              DO J=I-MM1O2,I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              L=L1
              SUMB=0.0D0
              DO J=I+1,N
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              DO JJ=N+1,MM1O2+I
                L=L-1
                J=2*N-JJ
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=(SUMA+SUMB*SIG)*FACT
            ENDDO
          ELSE
            DO I=1,MM1O2
              L=JI-1
              SUMA=0.0D0
              DO JJ=MM1O2+1-I,1,-1
                J=JJ+1
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              DO J=1,MM1O2+I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            ENDDO
            DO I=N0,NF
              SUMA=0.0D0
              L=JI-1
              DO J=I-MM1O2,I+MM1O2
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            ENDDO
            DO I=NF+1,N
              SUMA=0.0D0
              L=JI-1
              DO J=I-MM1O2,N
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              DO JJ=N+1,I+MM1O2
                L=L+1
                J=2*N-JJ
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            ENDDO
          ENDIF
        ELSE IF(ICON.EQ.25)THEN
!
!        El mateix que ICON=5 i ICON=20
!
          IF(COND)THEN
            L0=JI+MM1O2*M
            L1=L0+MM1O2
            DO I=1,MM1O2
              SUMA=0.0D0
              L=L0-1
              DO JJ=MM1O2+1-I,1,-1
                J=JJ+1
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA-A(L)*F(IJK)
              ENDDO
              DO J=1,I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              SUMB=0.0D0
              L=L1
              DO J=I+1,MM1O2+I
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=(SUMA+SUMB*SIG)*FACT
            ENDDO
            DO I=N0,NF
              SUMA=0.0D0
              L=L0-1
              DO J=I-MM1O2,I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              L=L1
              SUMB=0.0D0
              DO J=I+1,MM1O2+I
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=(SUMA+SUMB*SIG)*FACT
            ENDDO
            DO I=NF+1,N
              SUMA=0.0D0
              L=L0-1
              DO J=I-MM1O2,I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              L=L1
              SUMB=0.0D0
              DO J=I+1,N
                L=L-1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              DO JJ=N+1,MM1O2+I
                L=L-1
                J=2*N-JJ
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMB=SUMB+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=(SUMA+SUMB*SIG)*FACT
            ENDDO
          ELSE
            DO I=1,MM1O2
              L=JI-1
              SUMA=0.0D0
              DO JJ=MM1O2+1-I,1,-1
                J=JJ+1
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA-A(L)*F(IJK)
              ENDDO
              DO J=1,MM1O2+I
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            ENDDO
            DO I=N0,NF
              SUMA=0.0D0
              L=JI-1
              DO J=I-MM1O2,I+MM1O2
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            ENDDO
            DO I=NF+1,N
              SUMA=0.0D0
              L=JI-1
              DO J=I-MM1O2,N
                L=L+1
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              DO JJ=N+1,I+MM1O2
                L=L+1
                J=2*N-JJ
                IJK=IJKF0+(J-1)*IW(IPFX)
                SUMA=SUMA+A(L)*F(IJK)
              ENDDO
              IJK=IJKD0+(I-1)*IW(IPDX)
              DKF(IJK)=SUMA*FACT
            ENDDO
          ENDIF
        ELSE
          WRITE(6,*) ' Subroutine pderg mal utilitzada: ICON-->',ICON
          RETURN
        ENDIF
      ENDDO
      RETURN
1100  FORMAT(' Subroutine Pderg mal utilitzada: N,M,K------->',3I5)
      END
