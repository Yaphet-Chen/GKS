c----program for 2D mach reflection problem
cc---  f77 -r8 -O2 
      INTEGER  NP1,NP2
      PARAMETER ( NP1=484,NP2=124)
      
      DIMENSION DENS(NP1,NP2), XMOM(NP1,NP2), YMOM(NP1,NP2),
     *     ENER(NP1,NP2),
     *     U(NP1,NP2), SIG(NP1,NP2), 
     *     Q(NP1,NP2), 
     *     ADENS(NP1),AXMOM(NP1),AYMOM(NP1),AENER(NP1),
     *     ASDENS(NP1),ASXMOM(NP1),ASYMOM(NP1),ASENER(NP1),
     *     TFM(NP1),TFP(NP1),TFPX(NP1),TFE(NP1),
     *     FM(NP1),FP(NP1),FPX(NP1),FE(NP1),AAD(4),AAX(4),AAY(4),
     *     AAE(4),DB(2),XB(2),YB(2),EB(2),
     *     AAAD(NP1),AAAX(NP1),AAAY(NP1),AAAE(NP1)


      DX=1.0
      DY=1.0
      REY=10.0
      EP=1.0E-7
      AA=1.0


      PI=4.0*atan(1.0)
      GAM = 1.4
      RGAMH = SQRT(GAM/2.0)
      PROT1  = 1.0E-10
      CK=3.0
      EPS=0.4
      TIME   = 0.0
      TAUMAX = 10.0
      UBD    = 0.0
      IT     = 1

      OPEN(UNIT=1,FILE='data1.out',STATUS='unknown')
      OPEN(UNIT=2,FILE='data2.out',STATUS='unknown')

C     INITIAL CONDITIONS:
      NPIX=3
      NPX=482
      NPIY=3
      NPY=122
      DO 102 IY=3,NPY
         DO 101 IX = 1, NPX
            H=SQRT(3.0)*(IX-22.0)+3.0
            IF(IY.GE.H) THEN
               DENS(IX,IY) = 8.0 
               XMOM(IX,IY) = 57.1597
               YMOM(IX,IY) = -33.0012
               ENER(IX,IY) = 563.544
               GOTO 101
            ENDIF
            DENS(IX,IY) = 1.40
            XMOM(IX,IY) = 0.0
            YMOM(IX,IY) = 0.0
            ENER(IX,IY) = 2.5
 101     CONTINUE
 102  CONTINUE

      VIS=3.0*40.0/REY
      
C     COURANT CONDITIONS:

      XSTOP=22.0+200.0*SQRT(3.0)

 601  DO 201 IY = NPIY,NPY
         DO 202 IX=NPIX,NPX
            U(IX,IY)=SQRT(XMOM(IX,IY)*XMOM(IX,IY)+
     *           YMOM(IX,IY)*YMOM(IX,IY))/DENS(IX,IY)
            ELAM=0.25*(CK + 2.0)*DENS(IX,IY)/(ENER(IX,IY) 
     *           -0.5*XMOM(IX,IY)*XMOM(IX,IY)/DENS(IX,IY)
     *           -0.5 *YMOM(IX,IY)*YMOM(IX,IY)/DENS(IX,IY))
            IF(ELAM.LE.0.0) THEN 
               print *, 'ELAM NEGATIVE',DENS(IX,IY),XMOM(IX,IY),
     *              YMOM(IX,IY),ENER(IX,IY), IX,IY
            ENDIF
            SIG(IX,IY) = 1.0/(SQRT(ELAM))
            Q(IX,IY)   = U(IX,IY)/SIG(IX,IY)
            UABS = ABS(U(IX,IY))
            S    = UABS + RGAMH * SIG(IX,IY)
            UBD  = MAX(UBD,S)
            
 202     CONTINUE
 201  CONTINUE
      TBD  = 1.0/(UBD + PROT1)
      TAU =(1.0 - EPS) * MIN(TAUMAX,TBD)
      UBD = 0.0
      print *, ' IT, TIME, TAU: ',IT,TIME,TAU

      XEF= 22.0 + 40.0*SQRT(3.0)+ SQRT(3.0)*20.0*TIME/3.0

      
      IF(XEF.GE.XSTOP) THEN

         DO IY=3,122
            DO IX=3,362

               WRITE(1,10) FLOAT(IX)-2.0, FLOAT(IY)-2.0, DENS(IX,IY)
 10            FORMAT(E10.4,2X,E10.4,2X,E10.4)

               AP=0.4*(ENER(IX,IY)
     *              -0.5*(XMOM(IX,IY)**2+YMOM(IX,IY)**2)/
     *              DENS(IX,IY))
               WRITE(2,11) FLOAT(IX)-2.0, FLOAT(IY)-2.0, AP
 11            FORMAT(E10.4,2X,E10.4,2X,E10.4)

            enddo
         enddo

         stop

      ENDIF
      
      TIME=TIME+TAU
      IT=IT+1
C     THE FOLLOWING IS THE METHOD FOR CACULATING FLUXES.


      DO 161 IY=3,NPY
         IXEND=NPX
         DO 162 IX=3,IXEND
            ADENS(IX)=DENS(IX,IY)
            AXMOM(IX)=XMOM(IX,IY)
            AYMOM(IX)=YMOM(IX,IY)
            AENER(IX)=ENER(IX,IY)
 162     CONTINUE
         
          ADENS(1)= 8.0
          AXMOM(1)= 57.1597
          AYMOM(1)= -33.0012
          AENER(1)= 563.544
          ADENS(2)= 8.0
          AXMOM(2)= 57.1597
          AYMOM(2)= -33.0012
          AENER(2)= 563.544

          ADENS(IXEND+1)=ADENS(IXEND)
          AXMOM(IXEND+1)=AXMOM(IXEND)
          AYMOM(IXEND+1)=AYMOM(IXEND)
          AENER(IXEND+1)=AENER(IXEND)
          ADENS(IXEND+2)=ADENS(IXEND-1)
          AXMOM(IXEND+2)=AXMOM(IXEND-1)
          AYMOM(IXEND+2)=AYMOM(IXEND-1)
          AENER(IXEND+2)=AENER(IXEND-1)
          

         DO II=2,IXEND+1

            DFD1=ADENS(II)-ADENS(II-1)
            DFD2=ADENS(II+1)-ADENS(II)

            DFX1=AXMOM(II)-AXMOM(II-1)
            DFX2=AXMOM(II+1)-AXMOM(II)

            DFY1=AYMOM(II)-AYMOM(II-1)
            DFY2=AYMOM(II+1)-AYMOM(II)

            DFE1=AENER(II)-AENER(II-1)
            DFE2=AENER(II+1)-AENER(II)
            
            AAAD(II)=(SIGN(AA,DFD1)+SIGN(AA,DFD2))
     *           *(ABS(DFD1)*ABS(DFD2))/(ABS(DFD1)+ABS(DFD2)+EP)
            AAAX(II)=(SIGN(AA,DFX1)+SIGN(AA,DFX2))
     *           *(ABS(DFX1)*ABS(DFX2))/(ABS(DFX1)+ABS(DFX2)+EP)
            AAAY(II)=(SIGN(AA,DFY1)+SIGN(AA,DFY2))
     *           *(ABS(DFY1)*ABS(DFY2))/(ABS(DFY1)+ABS(DFY2)+EP)
            AAAE(II)=(SIGN(AA,DFE1)+SIGN(AA,DFE2))
     *           *(ABS(DFE1)*ABS(DFE2))/(ABS(DFE1)+ABS(DFE2)+EP)

         ENDDO


          DO 999 I= 2,IXEND

               I0=I-1
               I1=I+1
               I2=I+2

               AAD(1)=AAAD(I)/DX
               AAD(2)=ADENS(I)
               AAD(3)=ADENS(I1)
               AAD(4)=AAAD(I1)/DX
               
               AAX(1)=AAAX(I)/DX
               AAX(2)=AXMOM(I)
               AAX(3)=AXMOM(I1)
               AAX(4)=AAAX(I1)/DX

               AAY(1)=AAAY(I)/DX
               AAY(2)=AYMOM(I)
               AAY(3)=AYMOM(I1)
               AAY(4)=AAAY(I1)/DX
               
               AAE(1)=AAAE(I)/DX
               AAE(2)=AENER(I)
               AAE(3)=AENER(I1)
               AAE(4)=AAAE(I1)/DX

               CALL DXEE(AAD,AAX,AAY,AAE,DX,TAU,VIS,AFM,AFP,AFX,AFE)

               FM(I)=AFM
               FP(I)=AFP
               FPX(I)=AFX
               FE(I)=AFE
 999     CONTINUE

         DO 36 IP=3,IXEND
            ASDENS(IP)=ADENS(IP)+FM(IP-1)-FM(IP)
            ASXMOM(IP)=AXMOM(IP)+FP(IP-1)-FP(IP)
            ASYMOM(IP)=AYMOM(IP)+FPX(IP-1)-FPX(IP)
            ASENER(IP)=AENER(IP)+FE(IP-1)-FE(IP)
 36      CONTINUE


         DO 502 IX=3,IXEND
            DENS(IX,IY) = ASDENS(IX)
            XMOM(IX,IY) = ASXMOM(IX)
            YMOM(IX,IY) = ASYMOM(IX)
            ENER(IX,IY) = ASENER(IX)
 502     CONTINUE
 161  CONTINUE


      DO 141 IX=3, NPX
         INI=3
         IINI=2
         IEND=NPY
         DO 1309 IY=INI,IEND
            ADENS(IY)=DENS(IX,IY)
            AYMOM(IY)=-XMOM(IX,IY)
            AXMOM(IY)=YMOM(IX,IY)
            AENER(IY)=ENER(IX,IY)
 1309    CONTINUE
         ADENS(2)= ADENS(3)
         AXMOM(2)= -AXMOM(3)
         AYMOM(2)= AYMOM(3)
         AENER(2)= AENER(3)
         ADENS(1)= ADENS(4)
         AXMOM(1)= -AXMOM(4)
         AYMOM(1)= AYMOM(4)
         AENER(1)= AENER(4)
         
         IF(IX.LE.22) THEN
            ADENS(2)= 8.0
            AXMOM(2)= -33.0012
            AYMOM(2)= -57.1597
            AENER(2)= 563.544
            ADENS(1)= 8.0
            AXMOM(1)= -33.0012
            AYMOM(1)= -57.1597
            AENER(1)= 563.544
         ENDIF

         ADENS(IEND+1)= 1.40
         AXMOM(IEND+1)= 0.0
         AYMOM(IEND+1)= 0.0
         AENER(IEND+1)= 2.50
         ADENS(IEND+2)= 1.40
         AXMOM(IEND+2)= 0.0
         AYMOM(IEND+2)= 0.0 
         AENER(IEND+2)= 2.50

         IF(IX.LE.XEF) THEN
            ADENS(IEND+1)= 8.0
            AXMOM(IEND+1)= -33.0012
            AYMOM(IEND+1)= -57.1597
            AENER(IEND+1)= 563.544
            ADENS(IEND+2)= 8.0
            AXMOM(IEND+2)= -33.0012
            AYMOM(IEND+2)= -57.1597
            AENER(IEND+2)= 563.544
         ENDIF

          DO II=IINI,IEND+1
            DFD1=ADENS(II)-ADENS(II-1)
            DFD2=ADENS(II+1)-ADENS(II)

            DFX1=AXMOM(II)-AXMOM(II-1)
            DFX2=AXMOM(II+1)-AXMOM(II)

            DFY1=AYMOM(II)-AYMOM(II-1)
            DFY2=AYMOM(II+1)-AYMOM(II)

            DFE1=AENER(II)-AENER(II-1)
            DFE2=AENER(II+1)-AENER(II)
            
            AAAD(II)=(SIGN(AA,DFD1)+SIGN(AA,DFD2))
     *           *(ABS(DFD1)*ABS(DFD2))/(ABS(DFD1)+ABS(DFD2)+EP)
            AAAX(II)=(SIGN(AA,DFX1)+SIGN(AA,DFX2))
     *           *(ABS(DFX1)*ABS(DFX2))/(ABS(DFX1)+ABS(DFX2)+EP)
            AAAY(II)=(SIGN(AA,DFY1)+SIGN(AA,DFY2))
     *           *(ABS(DFY1)*ABS(DFY2))/(ABS(DFY1)+ABS(DFY2)+EP)
            AAAE(II)=(SIGN(AA,DFE1)+SIGN(AA,DFE2))
     *           *(ABS(DFE1)*ABS(DFE2))/(ABS(DFE1)+ABS(DFE2)+EP)
               
         ENDDO



         DO 3991  I =  IINI , IEND

                I0=I-1
                I1=I+1
                I2=I+2

            AAD(1)=AAAD(I)/DY
            AAD(2)=ADENS(I)
            AAD(3)=ADENS(I1)
            AAD(4)=AAAD(I1)/DY

            AAX(1)=AAAX(I)/DY
            AAX(2)=AXMOM(I)
            AAX(3)=AXMOM(I1)
            AAX(4)=AAAX(I1)/DY

            AAY(1)=AAAY(I)/DY
            AAY(2)=AYMOM(I)
            AAY(3)=AYMOM(I1)
            AAY(4)=AAAY(I1)/DY

            AAE(1)=AAAE(I)/DY
            AAE(2)=AENER(I)
            AAE(3)=AENER(I1)
            AAE(4)=AAAE(I1)/DY

            CALL DXEE(AAD,AAX,AAY,AAE,DY,TAU,VIS,AFM,AFP,AFX,AFE)



                FM(I)=AFM
                FP(I)=AFP
                FPX(I)=AFX
                FE(I)=AFE

 3991        CONTINUE

             DO 331 IP=INI,IEND
                ASDENS(IP)=ADENS(IP)+FM(IP-1)-FM(IP)
                ASXMOM(IP)=AXMOM(IP)+FP(IP-1)-FP(IP)
                ASYMOM(IP)=AYMOM(IP)+FPX(IP-1)-FPX(IP)
                ASENER(IP)=AENER(IP)+FE(IP-1)-FE(IP)
 331      CONTINUE

           
            
          DO 3196 IQ=INI,IEND
             DENS(IX,IQ)=ASDENS(IQ)
             XMOM(IX,IQ)=-ASYMOM(IQ)
             YMOM(IX,IQ)=ASXMOM(IQ)
             ENER(IX,IQ)=ASENER(IQ)
 3196     CONTINUE

 141   CONTINUE
      
       GOTO 601
 965   STOP
       END



      
      SUBROUTINE DXE(AXU,AYU,AL,A,B1,B2,C,XQ,Y1Q,Y2Q,ZQ)
      CC=2.0*C-(AXU*AXU+AYU*AYU+2.5/AL)*A
      BB=B2-AYU*A
      AA=B1-AXU*A
      ZQ=0.400*AL*AL*(CC-2.0*AXU*AA-2.0*AYU*BB)
      Y1Q=2.0*AL*(AA-AXU*ZQ/AL)
      Y2Q=2.0*AL*(BB-AYU*ZQ/AL)
      XQ=A-Y1Q*AXU-Y2Q*AYU-ZQ*(AXU*AXU+AYU*AYU+2.5/AL)
      RETURN
      END

      SUBROUTINE DXEE(AD,AX,AY,AE,DD,TAU,VIS,AFM,AFP,AFPX,AFE)
      DIMENSION AD(4),AX(4),AY(4),AE(4)
      CK = 3.0
      PI=4.0*atan(1.0)
      
      ADE1=AD(2)+0.5*DD*AD(1)
      AXM1=AX(2)+0.5*DD*AX(1)
      AYM1=AY(2)+0.5*DD*AY(1)
      AEN1=AE(2)+0.5*DD*AE(1)
      
      ADE2=AD(3)-0.5*DD*AD(4)
      AXM2=AX(3)-0.5*DD*AX(4)
      AYM2=AY(3)-0.5*DD*AY(4)
      AEN2=AE(3)-0.5*DD*AE(4)

      RADE1=1.0/ADE1
      RADE2=1.0/ADE2

      AE1=0.25*(CK+2)*ADE1/(AEN1
     *     -0.50*AXM1*AXM1/ADE1-0.50*AYM1*AYM1/ADE1)
      AXU1=AXM1*RADE1
      AYU1=AYM1*RADE1

      AE2=0.25*(CK+2)*ADE2/(AEN2
     *     -0.50*AXM2*AXM2/ADE2-0.50*AYM2*AYM2/ADE2)
      AXU2=AXM2*RADE2
      AYU2=AYM2*RADE2

      AW11=AD(1)
      BXW11=AX(1)
      BYW11=AY(1)
      CW11=AE(1)

      AW22=AD(4)
      BXW22=AX(4)
      BYW22=AY(4)
      CW22=AE(4)



      IF(AE1.LE.0.0.or.AE2.LE.0.0) THEN
         ADE1=AD(2)
         AXM1=AX(2)
         AYM1=AY(2)
         AEN1=AE(2)
         AW11=0.0
         BXW11=0.0
         BYW11=0.0
         CW11=0.0
         RADE1=1.0/ADE1
         AE1=0.25*(CK+2)*ADE1/(AEN1
     *        -0.50*AXM1*AXM1/ADE1
     *        -0.50*AYM1*AYM1/ADE1)
         AXU1=AXM1/ADE1
         AYU1=AYM1/ADE1

         ADE2=AD(3)
         AXM2=AX(3)
         AYM2=AY(3)
         AEN2=AE(3)
         AW22=0.0
         BXW22=0.0
         BYW22=0.0
         CW22=0.0
         RADE2=1.0/ADE2

         AE2=0.25*(CK+2)*ADE2/(AEN2
     *        -0.50*AXM2*AXM2/ADE2
     *        -0.50*AYM2*AYM2/ADE2)
         AXU2=AXM2/ADE2
         AYU2=AYM2/ADE2
      ENDIF

      CALL DXE(AXU1,AYU1,AE1,AW11,BXW11,BYW11,CW11,X1,Y1,YX1,Z1)
      CALL DXE(AXU2,AYU2,AE2,AW22,BXW22,BYW22,CW22,X2,Y2,YX2,Z2)
      
            
      RAE1=1.0/AE1
      RAE2=1.0/AE2
      
      TEU0=0.5*ERFC(-AXU1*SQRT(AE1))
      TEU1=AXU1*TEU0+0.50*EXP(-AE1*AXU1*AXU1)/SQRT(AE1*PI)
      TEU2=AXU1*TEU1+0.5*TEU0*RAE1
      TEU3=AXU1*TEU2+1.0*TEU1*RAE1
      TEU4=AXU1*TEU3+1.5*TEU2*RAE1
      TEU5=AXU1*TEU4+2.0*TEU3*RAE1
      TEU6=AXU1*TEU5+2.5*TEU4*RAE1
          
      TGU0=0.5*(ERFC(AXU2*SQRT(AE2)))
      TGU1=AXU2*TGU0-0.50*EXP(-AE2*AXU2*AXU2)/SQRT(AE2*PI)
      TGU2=AXU2*TGU1+0.5*TGU0*RAE2
      TGU3=AXU2*TGU2+1.0*TGU1*RAE2
      TGU4=AXU2*TGU3+1.5*TGU2*RAE2
      TGU5=AXU2*TGU4+2.0*TGU3*RAE2
      TGU6=AXU2*TGU5+2.5*TGU4*RAE2


      EY0=1.0          
      EY1=AYU1
      EY2=AYU1*EY1+0.5*EY0*RAE1
      EY3=AYU1*EY2+1.0*EY1*RAE1
      EY4=AYU1*EY3+1.5*EY2*RAE1
      EY5=AYU1*EY4+2.0*EY3*RAE1
      EI2=0.5*CK*RAE1
      EI4=(0.75*CK+0.25*CK*(CK-1.0))*RAE1*RAE1


      GY0=1.0
      GY1=AYU2
      GY2=AYU2*GY1+0.5*GY0*RAE2
      GY3=AYU2*GY2+1.0*GY1*RAE2
      GY4=AYU2*GY3+1.5*GY2*RAE2
      GY5=AYU2*GY4+2.0*GY3*RAE2
      GI2=0.5*CK*RAE2
      GI4=(0.75*CK+0.25*CK*(CK-1.0))*RAE2*RAE2
            

      ADE=ADE1*TEU0+ADE2*TGU0
      AXM=ADE1*TEU1+ADE2*TGU1
      AYM=ADE1*TEU0*EY1+ADE2*TGU0*GY1
      AEN=0.5*(ADE1*(TEU2+TEU0*EY2+TEU0*EI2)+
     *     ADE2*(TGU2+TGU0*GY2+TGU0*GI2))

      
      RADE=1.0/ADE
      AE0=0.25*(CK+2)*ADE/(AEN
     *     -0.50*AXM*AXM/ADE-0.50*AYM*AYM/ADE)
      AXU0=AXM*RADE
      AYU0=AYM*RADE
      RAE0=1.0/AE0


      AA1=ADE1/AE1
      AA2=ADE2/AE2
      
c      TE=2.0*AE0*VIS+TAU*ABS(AA1-AA2)/(AA1+AA2)

      TE=0.01*TAU+TAU*ABS(AA1-AA2)/(AA1+AA2)
      
      SP=-TAU/TE
      SW=EXP(SP)
      SE1=TE*(1.0-SW)
      SET0=TE*(-TAU+SE1)
      SET1=-TE*TAU*SW+TE*TE*(1.0-SW)
      SET2=TE*(-TAU+SE1)+SET1
      

      FMI1=SE1*ADE1*TEU0-SET1*
     *     (X1*TEU1+Y1*TEU2+YX1*TEU1*EY1+
     *     Z1*(TEU3+TEU1*EY2+TEU1*EI2))
      FPI1=SE1*ADE1*TEU1-SET1*
     *     (X1*TEU2+Y1*TEU3+YX1*TEU2*EY1+
     *     Z1*(TEU4+TEU2*EY2+TEU2*EI2))
      FPIX1=SE1*ADE1*TEU0*EY1-SET1*
     *     (X1*TEU1*EY1+Y1*TEU2*EY1
     *     +YX1*TEU1*EY2+
     *     Z1*(TEU3*EY1+TEU1*EY3+TEU1*EY1*EI2))
      FEI1=0.5*SE1*ADE1*(TEU2+TEU0*EY2+TEU0*EI2)
     *     -0.5*SET1*
     *     (X1*(TEU3+TEU1*EY2+TEU1*EI2)+
     *     Y1*(TEU4+TEU2*EY2+TEU2*EI2)+
     *     YX1*(TEU3*EY1+TEU1*EY3+TEU1*EY1*EI2)+
     *     Z1*(TEU5+TEU1*EY4+TEU1*EI4
     *     +2.0*TEU3*EY2+2.0*TEU3*EI2+2.0*TEU1*EY2*EI2))

      
      FMI2=SE1*ADE2*TGU0-SET1*
     *     (X2*TGU1+Y2*TGU2+YX2*TGU1*GY1+
     *     Z2*(TGU3+TGU1*GY2+TGU1*GI2))
      FPI2=SE1*ADE2*TGU1-SET1*
     *     (X2*TGU2+Y2*TGU3+YX2*TGU2*GY1+
     *     Z2*(TGU4+TGU2*GY2+TGU2*GI2))
      FPIX2=SE1*ADE2*TGU0*GY1-SET1*
     *     (X2*TGU1*GY1+Y2*TGU2*GY1
     *     +YX2*TGU1*GY2+
     *     Z2*(TGU3*GY1+TGU1*GY3+TGU1*GY1*GI2))
      FEI2=0.5*SE1*ADE2*(TGU2+TGU0*GY2+TGU0*GI2)
     *     -0.5*SET1*
     *     (X2*(TGU3+TGU1*GY2+TGU1*GI2)+
     *     Y2*(TGU4+TGU2*GY2+TGU2*GI2)+
     *     YX2*(TGU3*GY1+TGU1*GY3+TGU1*GY1*GI2)+
     *     Z2*(TGU5+TGU1*GY4+TGU1*GI4
     *     +2.0*TGU3*GY2+2.0*TGU3*GI2+2.0*TGU1*GY2*GI2))
      
      AMI=FMI1+FMI2
      API=FPI1+FPI2
      APIX=FPIX1+FPIX2
      AEI=FEI1+FEI2


      SETV=SET1

      FM1=ADE1*(SE1*TEU1)-(SETV)*
     *     (X1*TEU2+Y1*TEU3+YX1*TEU2*EY1
     *     +Z1*(TEU4+TEU2*EY2+TEU2*EI2))
      FP1=ADE1*(SE1*TEU2)-(SETV)*
     *     (X1*TEU3+Y1*TEU4+YX1*TEU3*EY1
     *     +Z1*(TEU5+TEU3*EY2+TEU3*EI2))
      FPX1=ADE1*(SE1*TEU1*EY1)-(SETV)*
     *     (X1*TEU2*EY1+Y1*TEU3*EY1
     *     +YX1*TEU2*EY2
     *     +Z1*(TEU4*EY1+TEU2*EY3+TEU2*EY1*EI2))
      FE1=ADE1*(0.5*SE1*(TEU3+TEU1*EI2+TEU1*EY2))
     *     -0.5*(SETV)*
     *     (X1*(TEU4+TEU2*EY2+TEU2*EI2)+
     *     Y1*(TEU5+TEU3*EY2+TEU3*EI2)+
     *     YX1*(TEU4*EY1+TEU2*EY3+TEU2*EY1*EI2)+
     *     Z1*(TEU6+TEU2*EY4+TEU2*EI4+2.0*TEU4*EY2
     *     +2.0*TEU4*EI2+2.0*TEU2*EY2*EI2))
      
      FM2=ADE2*(SE1*TGU1)-(SETV)*
     *     (X2*TGU2+Y2*TGU3+YX2*TGU2*GY1
     *     +Z2*(TGU4+TGU2*GY2+TGU2*GI2))
      FP2=ADE2*(SE1*TGU2)-(SETV)*
     *     (X2*TGU3+Y2*TGU4+YX2*TGU3*GY1
     *     +Z2*(TGU5+TGU3*GY2+TGU3*GI2))
      FPX2=ADE2*(SE1*TGU1*GY1)-(SETV)*
     *     (X2*TGU2*GY1+Y2*TGU3*GY1
     *     +YX2*TGU2*GY2
     *     +Z2*(TGU4*GY1+TGU2*GY3+TGU2*GY1*GI2))
      FE2=ADE2*(0.5*SE1*(TGU3+TGU1*GI2+TGU1*GY2))
     *     -0.5*(SETV)*
     *     (X2*(TGU4+TGU2*GY2+TGU2*GI2)+
     *     Y2*(TGU5+TGU3*GY2+TGU3*GI2)+
     *     YX2*(TGU4*GY1+TGU2*GY3+TGU2*GY1*GI2)+
     *     Z2*(TGU6+TGU2*GY4+TGU2*GI4+2.0*TGU4*GY2+2.0*TGU4*GI2+
     *     2.0*TGU2*GY2*GI2))
            
cc--- in this program, a single slope is used for g, which is 
cc--- consistent with the 93 paper. The VKI lecture describes
cc--- the schemes with two slopes for g, which is more advanced.

      AWL=(AD(3)-AD(2))/DD
      BXWL=(AX(3)-AX(2))/DD
      BYWL=(AY(3)-AY(2))/DD
      CWL=(AE(3)-AE(2))/DD

      
      CALL DXE(AXU0,AYU0,AE0,AWL,BXWL,BYWL,CWL,XL,YL,YXL,ZL)

      T0=1.0
      T1=AXU0
      T2=AXU0*T1+0.5*T0*RAE0
      T3=AXU0*T2+1.0*T1*RAE0
      T4=AXU0*T3+1.5*T2*RAE0
      T5=AXU0*T4+2.0*T3*RAE0
      T6=AXU0*T5+2.5*T4*RAE0

 
      Y0=1.0          
      Y1=AYU0
      Y2=AYU0*Y1+0.5*Y0*RAE0
      Y3=AYU0*Y2+1.0*Y1*RAE0
      Y4=AYU0*Y3+1.5*Y2*RAE0
      Y5=AYU0*Y4+2.0*Y3*RAE0
      EI2=0.5*CK*RAE0
      EI4=(0.75*CK+0.25*CK*(CK-1.0))*RAE0*RAE0

            
      TRIU=(XL*T1+YL*T2+YXL*T1*Y1+
     *     ZL*(T3+T1*Y2+T1*EI2))
      TRIU2=(XL*T2+YL*T3+YXL*T2*Y1+
     *     ZL*(T4+T2*Y2+T2*EI2))
      TRIUY1=(XL*T1*Y1+YL*T2*Y1
     *     +YXL*T1*Y2+
     *     ZL*(T3*Y1+T1*Y3+T1*Y1*EI2))
      TRIEU=0.5*
     *     (XL*(T3+T1*Y2+T1*EI2)+
     *     YL*(T4+T2*Y2+T2*EI2)+
     *     YXL*(T3*Y1+T1*Y3+T1*Y1*EI2)+
     *     ZL*(T5+T1*Y4+T1*EI4
     *     +2.0*T3*Y2+2.0*T3*EI2+2.0*T1*Y2*EI2))

      TRIU3=(XL*T3+YL*T4+YXL*T3*Y1
     *     +ZL*(T5+T3*Y2+T3*EI2))
      TRIU2Y1=(XL*T2*Y1+YL*T3*Y1
     *     +YXL*T2*Y2
     *     +ZL*(T4*Y1+T2*Y3+T2*Y1*EI2))
      TRIEU2=0.5*
     *     (XL*(T4+T2*Y2+T2*EI2)+
     *     YL*(T5+T3*Y2+T3*EI2)+
     *     YXL*(T4*Y1+T2*Y3+T2*Y1*EI2)+
     *     ZL*(T6+T2*Y4+T2*EI4+2.0*T4*Y2+2.0*T4*EI2+
     *     2.0*T2*Y2*EI2))




ccc---------------------------End two slopes for g


      AS=(-(TRIU)*SET2+SE1*ADE-AMI)/SET0
      BS=(-(TRIU2)*SET2+SE1*AXM-API)/SET0
      BSX=(-(TRIUY1)*SET2+SE1*AYM-APIX)/SET0
      CS=(-(TRIEU)*SET2+SE1*AEN-AEI)/SET0

      CALL DXE(AXU0,AYU0,AE0,AS,BS,BSX,CS,A1,B1,C1,D1)
            
CCFC----------------------------------------------------

      ATRIU=A1*T1+B1*T2+C1*T1*Y1+D1*(T3+T1*EI2+T1*Y2)
      ATRIUY1=A1*T1*Y1+B1*T2*Y1+C1*T1*Y2
     *     +D1*(T3*Y1+T1*EI2*Y1+T1*Y3)
      ATRIU2=A1*T2+B1*T3+C1*T2*Y1+D1*(T4+T2*EI2+T2*Y2)
      ATRIEU=0.5*(
     *     A1*(T3+T1*EI2+T1*Y2)+
     *     B1*(T4+T2*EI2+T2*Y2)+
     *     C1*(T3*Y1+T1*EI2*Y1+T1*Y3)+
     *     D1*(T5+T1*Y4+T1*EI4+2.0*T3*Y2
     *     +2.0*T3*EI2+2.0*T1*EI2*Y2))
      
      A=ADE*(TAU-SE1)
      B=SET2
      D=(0.5*TAU*TAU
     *     -TE*TAU+TE*TE*(1.0-SW))
      
      AFM=A*T1+B*TRIU2+D*ATRIU+(FM1+FM2)
      AFP=A*T2+B*TRIU3+D*ATRIU2+(FP1+FP2)
      AFPX=A*T1*Y1+B*TRIU2Y1+D*ATRIUY1
     *     +(FPX1+FPX2)
      AFE=0.5*A*(T3+T1*EI2+T1*Y2)+B*TRIEU2+
     *     D*ATRIEU+(FE1+FE2)

      RETURN
      END




