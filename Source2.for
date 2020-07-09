       PROGRAM CART3D

C ------------- 04.05.18
CC          11.05.18
CCC            10.06.18 --------
cv      ---- 14.10.18

        PARAMETER (NC=75000,NR=300,NR2=100,ICMAX=100)
        IMPLICIT REAL*8 (A-H,O-Z)
        REAL*8 M1,M2,M3,MSP     
C        ,KAP2,M3A
        REAL*8 XCL(NC),YCL(NC),
C         ZCL(NC),
     ,         UCL(NC),VCL(NC),
C         WCL(NC),
     ,         DD2(ICMAX),CUT(NR),DVR(ICMAX),DDU(ICMAX),
     ,         DDV(ICMAX),DDX(ICMAX),DDY(ICMAX),CUTB(NR)
        
         INTEGER NM2(NC),NUM2(ICMAX)
c        NCP(NC),

C        REAL*8 BINT0(NR2),BINT1(NR2) ! precomputed integrals
      

       DOUBLE PRECISION
     ,     DSIN,DCOS,DBLE,DEXP,DTAN,DASIN,DABS,DSQRT,
     ,     DATAN2,DLOG

     
        ROmega(xx)=(260.d0*dexp(-xx/150.d0-(12.96d0 /(xx*xx)))+
     +  360.d0*dexp(-xx/3.3d0-0.1d0/xx))

       character(len=*), parameter::path=
     =  "C:\solution\Final\"




      OPEN (10,FILE=path//'INPdec.DAT')
      OPEN (11,FILE=path//'GALAC.DAT')
      OPEN (15,FILE=path//'LASTN.DAT',
     ,            ACCESS='DIRECT',RECL=64)
      OPEN (16,FILE=path//'UVVEL.DAT',
     ,            ACCESS='DIRECT',RECL=16*(NC))
cv     ,            ACCESS='DIRECT',RECL=4*(NC))

      OPEN (25,FILE=path//'COORD.DAT',
     ,           ACCESS='DIRECT',RECL=16*(NC))
cv     ,           ACCESS='DIRECT',RECL=4*(NC))

      OPEN (26,FILE=path//'CUT.DAT',
     ,           ACCESS='DIRECT',RECL=8*NR)


      open (45,file=path//'TESTCoorD.dat')
      open (55,file=path//'test2.dat')
      open (56,file=path//'testClump.dat')
      OPEN (57,FILE=path//'MOMENTUM_FORCE.DAT')
      OPEN (58,FILE=path//'E_KIN.DAT')    ! kinetic energy
      OPEN (59,FILE=path//'FULL_IMP_MOM.DAT')    ! full impulse momentum
c      OPEN (60,FILE=path//'\MOTION.DAT',
c     ,            ACCESS='DIRECT',RECL=32*(NC))
      OPEN (60,FILE=path//'bartest.DAT')

      open (27,file=path//'XY.DAT')
      open (28,file=path//'NUMCOLL.DAT')
      
C      OPEN (29,FILE=PATH//'BARPOT.DAT')
C      OPEN (30,FILE=PATH//'SPIRPOT.DAT')
      OPEN (31,FILE=PATH//'POT.DAT')
      
cv      OPEN (62,file=path//'test_coll\8hTEST_COLL1.DAT')
cv      OPEN (63,file=path//'test_coll\8hTEST_COLL2.DAT')
      
C      OPEN (98,file=path//'TEST_COLL1.DAT')
C      OPEN (99,file=path//'TEST_COLL2.DAT')      

        READ(10,1) NPUSK,NCOL,IDSTR,NST3,CR3,NCC,RTH,
     ,             SCALE,EPSR,NITM,RGI,
     ,             R0S,TCO,RG,RMIN,R0,OMP,
     ,             M1,M2,M3,
     ,             B1,A2,B2,A3,
     ,             ALF,PA,DISP,RGO,
     ,             IJR,KLR,IJT,KLT,IJV,KLV,
     ,             IJZ,KLZ,IJC,KLC,MSP,CSEC,
     ,             BT,KT,KREC,KSCR,DRAN,
     ,             GMO,COFT,DCL,
     ,             DSPCH,GRD,MSP4,A4A2,HI44,HI22,
     ,             rcl0,tcl0,
     ,             ABAR,CBAR,BARM,
     ,             HZCL,HZMAX,ELAST
1     FORMAT(I1/I1/I1/I1/F5.2/I5/F5.2/
     /         F5.2/F6.4/I3/F5.2/
     /         F5.2/F10.5/F5.2/F5.2/F5.2/F5.2/
     /         3(E9.4/),
     ,         F7.4/F4.2/F5.1/F6.1/
     /         F6.3/F6.2/F5.1/F5.1/
     ,         10(I7/),F4.1/F4.2/
     /         F5.2/I7/I5/I5/F5.2/
     /         D7.1/F3.1/F8.5/
     /         F5.3/F6.3/I1/F5.2/F7.1/F7.1/
     /         D7.1/D7.1/
     /         F5.2/F5.2/D8.1/
     /         F5.3/F8.3/F5.3)

        WRITE (*,2) NPUSK,NCOL,IDSTR,NST3,CR3,NCC,RTH,
     ,             SCALE,EPSR,NITM,RGI,R0S,
     ,             TCO,RG,RMIN,R0,OMP,
     ,             M1,M2,M3,B1,A2,B2,A3,
     ,             ALF,PA,DISP,RGO,
     ,             IJR,KLR,IJT,KLT,IJV,KLV,IJZ,KLZ,IJC,KLC,
     ,             MSP,CSEC,BT,KT,KREC,KSCR,DRAN,
     ,             GMO,COFT,DCL,
     ,             DSPCH,GRD,MSP4,A4A2,HI44,HI22,
     ,             rcl0,tcl0,
     ,             ABAR,CBAR,BARM,
     ,             HZCL,HZMAX,ELAST
2      FORMAT(1X,'NPUSK=',I1,2X,'NCOL=',I1,2X,'IDSTR=',I1,2X,'NST3=',I4/
     ,         1X,'CR3=',F5.2,2X,'NCC=',I5,2X,'RTH=',F5.2/
     /         1X,'SCALE=',F5.2,2X,'EPSR=',F6.4,2X,'NITM=',I3,
     ,                2X,'RGI=',F5.2/
     /         1X,'R0S=',F5.2,2X,'TCO=',F8.2,2X,'RG=',F5.2/
     ,         1X,'RMIN=',F5.2,2X,'R0=',F5.2,2X,'OMP=',F5.2/
     ,         1X,'M1=',E10.2,2X,'M2=',E10.2,2X,'M3=',E10.2/
     ,         1X,'B1=',E10.2,2X,'A2=',E10.2,2X,'B2=',E10.2,
     ,         1X,'A3=',D10.2/
     /         1X,'ALF=',F7.4,2X,'PA=',F5.1,2X,'DISP=',F6.1,
     ,             2X,'RGO=',F5.1/
     /         1X,'IJR=',I7,2X,'KLR=',I7,2X,'IJT=',I7/
     /         1X,'KLT=',I7,2X,'IJV=',I7,2X,'KLV=',I7/
     /         1X,'IJZ=',I7,2X,'KLZ=',I7,1X',IJC=',I7,
     ,                     1X,'KLC=',I7/
     ,         1X,'MSP=',F3.0,2X,'CSEC=',F7.3,2X,'BT=',F8.3/
     /         1X,'KT=',I10,2X,'KREC=',I5,2X,'KSCR=',I5,2X,
     ,                'DRAN=',F5.2/
     / 1X,'GMO=',D7.1,' M_sun',2X,'COFT=',F3.1,2X,'DCL=',F6.4/
     /    1X,'DSPCH=',F5.3,2X,'GRD=',F6.3/
     ,    1X,'MSP4=',I1,2X,'A4A2=',F5.2,2X,'HI44=',F7.1,
     ,      2X,'HI22=',F7.1/
     ,     1x,'rcl0=',f5.2,2x,'tcl0=',f6.2/
     /     1X,'ABAR=',F5.2,1X,'CBAR=',F5.2,2X,'BARM=',D8.1/
     /     1X,'HZCL=',F5.3,2X,'HZMAX=',F8.3,2X,
     /        2X,'ELAST=',F5.3)
     

     
C======CONSTANTS=======

        G=4.446666666667D-6
        PI=2.D0*DASIN(1.D0)
        PID2=0.5D0*PI
        PI2=2.D0*PI
       PA=PA/180.D0*PI

       TGI=DTAN(PA)
       HI44=HI44/180.*PI
       HI22=HI22/180.*PI

       CTGI=1.D0/TGI
       TGI=TGI*COFT


         COTM=MSP*CTGI

        OPE2=0.5D0*(1.D0+ELAST)

         GMO=G*GMO
          OMPS=OMP*OMP
           OMP2=2.D0*OMP

           R=R0

            OM0=ROmega(r)/r
             OM2 = OM0*OM0

            PER=PI2/OM0
      
      READ (15,REC=1) NREC,RLIN,RC,RLOT,TIME,DRG,TAU,NL
      
      PRINT *, '---------'
      PRINT *, 'RLIN = ',RLIN
      PRINT *, 'RC = ',RC
      PRINT *, 'RLOT = ',RLOT
      PRINT *, 'DRG = ',DRG
      PRINT *, '---------'
c //////////////////////////////////////////////
c      TAU=1.D-4
       write(*,*) "TAU = ", TAU


C      WRITE (*,*) 'LAST NREC=',NREC
      WRITE (*,*) 'START NREC=',NREC   
      WRITE (*,*) 'TIME=',TIME
C      WRITE (*,*) 'INPUT NREC'
cv      READ  (*,*)  NREC
C      WRITE(*,*) 'NREC = 1'
C      WRITE (*,*) 'NREC=',NREC
      WRITE (*,*) '  NREC GOOD ? --> ENTER'
C      READ (*,*)


      RZIL=RLIN-CSEC
      RZIR=RLIN+CSEC

        RZOL=RLOT-CSEC
        RZOR=RLOT+CSEC
      
      PRINT *, 'RZIL = ',RZIL
      PRINT *, 'RZOR = ',RZOR
      
C ---  COMMENT BY CCC ----> FOR K=FUNCTION OF R
CCC   WN0=MSP*CTGI/R0S*DSQRT((RKO-R0S)/(R0S-RKI))
CCC       AMF=ALF*OM2*R0S/WN0
CCC        AMFM=AMF*MSP
CCC      AMFW=AMF*WN0


C  AMPLITUDE BELOW FOR CONSTANT PITCH ANGLE --->
            WNAM=MSP*DABS(CTGI)/R0
            AMF=ALF*OM2*R0/WNAM
            AMFM=AMF*MSP
      

          READ (25,REC=NREC) XCL,YCL !,ZCL
          READ (16,REC=NREC) UCL,VCL !,WCL
          READ (26,REC=1) CUT
          
c ------------ BAR CONSTANTS ----------->
        CUTB = 1.D0 - CUT

        CHI0 = PID2-0.3D0
        
        THB0 = -CHI0+CTGI*DLOG(RLIN/R0S)
        RB = RLIN
        
         RB3 = RB*RB*RB
          AMPB = -AMFM
           OMB = OMP
C            THB0 = 0.D0 ! 0.25D0*PI     ! BAR INITIAL PHASE
             THB02 = THB0*2.D0

        THB = THB0
         SNTHB = DSIN(THB)
          CSTHB = DCOS(THB) 
        
        PRINT *, 'RB = ', RB  
        PRINT *, 'THB = ',THB0

        DRWR = 0.1D0
        NI = IDINT(RG/DRWR)
        print*, 'NI = ', ni
C POTENTIAL WRITE --->
        CHI0 = 0.D0
         AS = AMFM
          AMPB = -AMFM      ! /MSP
C          AMPB = 0.D0                     ! TEST NO BAR
           PRINT *, 'AS = ',AS
            PRINT *, 'AB = ',AMPB

            
      GOTO 207
        DO 203 J = -NI,NI
         DO 204 I = -NI,NI
            X = DRWR*DBLE(I)
             Y = DRWR*DBLE(J)
              X2 = X*X
               Y2 = Y*Y
                R = DSQRT(X2+Y2)
                 TH = DATAN2(Y,X)

C SPIR --->            
            IF (R.GT.RZIL .AND. R.LE.RZOR) THEN
               II=IDINT(R/DRG)
               CUTI=CUT(II)
        CUTI=CUTI+(CUT(II+1)-CUTI)*(R/DRG-DBLE(II))
            ELSE
               CUTI=0.D0
            ENDIF
            
            AS = AMFM*CUTI
            
            IF (R .LT. 1.D-5) THEN
              AS = 0.D0
              GOTO 205
            ENDIF
            
            CHIS = CTGI*DLOG(R/R0S)-TH
             CHIS2 = 2.D0*CHIS
205           PHIS = AS*DCOS(CHIS2)
            
c205         WRITE (30,30) X,Y,PHIS

C BAR --->
         IF (R .GT. RB) THEN
              RBR = RB/R
               RBR3 = RBR*RBR*RBR
                ABAR = -1.D0*AMPB*RBR3
            ELSE
              RBR = R/RB
               RBR3 = RBR*RBR*RBR
                ABAR = AMPB*(RBR3 - 2.D0)
            ENDIF

             CHIB = -TH-THB0
              CHIB = 2.D0*CHIB
               PHIB = ABAR*DCOS(CHIB)

            IF (R.LE.RZIR) THEN
               II = IDINT(R/DRG)
               CUTI=CUTB(II)
      CUTI=CUTI+(CUTB(II+1)-CUTI)*(R/DRG-DBLE(II))         
            ELSE
               CUTI = 0.D0
            ENDIF

            PHIB = PHIB*CUTI
C            WRITE (29,30) X,Y,PHIB
            
            SUMPOT = PHIB+PHIS
            
            WRITE (31,30) X,Y,PHIS/AMFM,PHIB/AMFM,SUMPOT/AMFM
30           FORMAT (1X,5(F12.5,1X))
            
204      CONTINUE
203     CONTINUE
207   CONTINUE
C       STOP 'BAR AND SPIR POT'

C ------------ CONSTANTS FOR COLLISIONS ------->

       RCL=DCL*0.5d0

         DEF=DCL 
c  +0.001D0

           DEF2=DEF*DEF

         NCOLL = 0
C  R_MIN & R_MAX --- FOR COLLIDING PARTICLES

C              RMIN=RZIL
             RMIN = 0.1D0
              RMIN2=RMIN*RMIN

C   ----------  BELOW "5" IS RATHER ARBITRARY
c             RMAX=RG-5.D0*DCL     !!!!

              RMAX=RG-0.1D0

               RMAX2=RMAX*RMAX      ! CHECK

C ====== TEST PARTICLES MOTION =====
c       XCL(1) = 5.D0
c        YCL(1) = 3.D0
c         X = XCL(1)
c          Y = YCL(1)
c           R = DSQRT(X*X+Y*Y)
c            OM = ROMEGA(R)
c       UCL(1) = OM*Y/R
c       VCL(1) = -OM*X/R
cc       AMFM = 0.D0
c       AMPB = 0.D0
c       PRINT *, 'R = ', R
CCC TEST WRITE        
c           DO 346 I=1,NC
c            XC = XCL(I)
c            YC = YCL(I)
c            WRITE (27,27) XC,YC
cC 27           FORMAT (1X,2(F8.3,2X))
c346       CONTINUE
CCC END TEST WRITE

C ========  TRAJECTORY ===============

         TAU2=TAU*0.5D0
          TAU6=TAU/6.D0
           KWRT=0
            KPRR=0

C ==========  TIME CYCLE =============

      WRITE (*,*) '### START TIME CYCLE ###'

      DO 300 KS=1,KT

        CAMPL1=1.D0-DEXP(-BT*OM0*TIME)
        CAMPL2=1.D0-DEXP(-BT*OM0*(TIME+TAU2))
        CAMPL3=1.D0-DEXP(-BT*OM0*(TIME+TAU))
c        CAMPL1 = 1.0
c        CAMPL2 = 1.0
c        CAMPL3 = 1.0


           OMPT=OMP*TIME
             KPRR=KPRR+1
           IF (KPRR.EQ.KSCR) THEN
             WRITE (*,4) KS,CAMPL1,TIME
4             FORMAT (1X,'KS==',I7,2X,'CAMPL1 ='
     ,           ,F6.3,3X,'TIME =',F6.3)
               KPRR=0
           END IF

      AMT1=CAMPL1*AMFM
       AMT2=CAMPL2*AMFM
        AMT3=CAMPL3*AMFM
      
C --------- COLLISIONS ---------------->

       IF (NCOL.EQ.0) GOTO 303

CCCCCCCCC        NM1=0 !!!!!!!!!!!!!!  NOT NEED!!!!!!!!!

             NM2=0
              NUMCOL = 0

      
        DO 350 N1=1,NC-1
        
          IF(NM2(N1).GT.0) GOTO 350

          X1=XCL(N1)
          Y1=YCL(N1)
          R11=X1*X1+Y1*Y1

          IF (R11.GT.RMAX2) GOTO 350

          U1=UCL(N1)
          V1=VCL(N1)

          I1 = IDINT(X1/RCL)
          IF (X1.LT.0) I1=I1-1 

          J1 = IDINT(Y1/RCL)
          IF (Y1.LT.0) J1=J1-1 

          IC=0  !  NUMBER OF PARTIC 2 SUSPECTED FOR COLL
      
          DO 332 N2=N1+1,NC
         
            IF(NM2(N2).GT.0) GOTO  332

            X2=XCL(N2)

            I2 = IDINT(X2/RCL)
            IF (X2.LT.0) I2=I2-1  

            IF (I1.LT.I2-2 .OR. I1.GT.I2+2) GOTO 332

            Y2=YCL(N2)
      
            J2 = IDINT(Y2/RCL)
            IF (Y2.LT.0) J2=J2-1

            IF (J1.LT.J2-2 .OR .J1.GT.J2+2) GOTO 332
            R22 = X2*X2+Y2*Y2
            IF (R22 .GT. RMAX2) GOTO 332
      
            DX = X1 - X2
            DY = Y1 - Y2
      
            DU=U1-UCL(N2)
            DV=V1-VCL(N2)
         

            DVDR=DU*DX+DV*DY
                
            IF (DVDR.GE.0.D0) THEN        !!  NO COLL !!!!!
               GOTO 332
            ELSE
               IC=IC+1
        IF (IC.GE.ICMAX) STOP ' == IC > IC max == 100 '
cv         PRINT *, 'IC = ', IC
             DRM2=DX*DX+DY*DY
              NUM2(IC)=N2
               DD2(IC)=DRM2
                DVR(IC)=DVDR
                 DDU(IC)=DU
                  DDV(IC)=DV
                   DDX(IC)=DX
                    DDY(IC)=DY
      
C ****** ASIGN DD X&Y !!!!!!!!!!!!!!!!!
          !     GOTO 332
           
            ENDIF

332       CONTINUE

C ------ IF IC = 0 cloud N1 does not collide WITH ANY OTHERS !!!!!!!!!!!!!!!!

          IF (IC.EQ.0) GOTO 350

          DRM2=DD2(1)
           NMIN = 1

CHECK         BELOW IF **************************
          IF (IC.EQ.1 .AND. DRM2.LE.DEF2) THEN
c            NMIN=1
            GOTO 315
          ENDIF

          DO 333 IN=2,IC
          
            D12=DD2(IN)  
              
            IF (DRM2.GT.D12) THEN     
              DRM2=D12
               NMIN=IN 
            ENDIF
333       CONTINUE
         
          IF (DRM2.GE.DEF2) GOTO 350     ! DIST > DIAM -> NO COLL
      
315       DU=DDU(NMIN)
           DV=DDV(NMIN)
            DVDR=DVR(NMIN)

          DX = DDX(NMIN)
           DY = DDY(NMIN)

c         DRM2=DD2(NMIN)     ! probably not need

          DVM2=DU*DU+DV*DV
           DVM = DSQRT(DVM2)
            DRM = DSQRT(DRM2)
            
          N2 = NUM2(NMIN)
           U2 = UCL(N2)
            V2 = VCL(N2)
               
          COSAL=DVDR/(DVM*DRM)

C ------ TO PREVENT COS BE < -1

          IF (COSAL.LE.-1.D0) THEN
            
                UCL(N1)=U1-OPE2*DU
                VCL(N1)=V1-OPE2*DV
c                 WCL(NIJ)=W1-OPE2*DW
                                          
                UCL(N2)=U2+OPE2*DU
                VCL(N2)=V2+OPE2*DV
c                 WCL(N2)=W1+OPE2*DW

          ELSE

                SINAL=DSQRT(1.D0-COSAL*COSAL)

                DIMP=DRM*SINAL

                CSHI=2.D0*DIMP*DIMP/DEF2-1.D0
                SNHI=DSQRT(1.D0-CSHI*CSHI)


                SNHIAL=SNHI/SINAL
                BR1=ELAST*(CSHI-COSAL*SNHIAL)-1.D0
                BR2=ELAST*SNHIAL*DVM/DRM

                DVX=0.5D0*(DU*BR1+DX*BR2)
                DVY=0.5D0*(DV*BR1+DY*BR2)
                
c                 DVZ=0.5D0*(DW*BR1+DZ*BR2)

                UCL(N1)=U1+DVX
                VCL(N1)=V1+DVY
c                 WCL(NIJ)=W1+DVZ
               
                UCL(N2)=U2-DVX
                VCL(N2)=V2-DVY

          ENDIF

          NUMCOL=NUMCOL+1
          NCOLL = NCOLL+1

          NM2(N2)=1

C            R22 = XCL(N2)*XCL(N2)+YCL(N2)*YCL(N2)
C            WRITE (28,*) N1,N2,R22,DRM,DVM,COSAL 
350     CONTINUE

C TEST WRITE COLL
C        WRITE (28,28) KS, NUMCOL
C28      FORMAT (1X, 2(I7,2X))
C        PRINT *,NUMCOL,RMAX2
        WRITE (28,28) KS, NUMCOL, NCOLL
28      FORMAT (1X, 3(I7,2X))
C      STOP 'END COLL'
C ================ END  COLL =====>>>>>

303   CONTINUE        ! AFTER COLL

C CALCULATIONS BEFORE RUNGE-KUTTA ------>
         OMBT = -OMB*TIME
          OMBT2 = -OMB*(TIME+TAU2)
           OMBT4 = -OMB*(TIME+TAU)
         
         HIB = MSP*(OMBT-THB0)
          HIB2 = MSP*(OMBT2-THB0)
           HIB4 = MSP*(OMBT4-THB0)
         
         SNHIB = DSIN(HIB)
          SNHIB2 = DSIN(HIB2)
           SNHIB4 = DSIN(HIB4)
         
         CSHIB = DCOS(HIB)
          CSHIB2 = DCOS(HIB2)
           CSHIB4 = DCOS(HIB4)
         
         OMPT = -OMP*TIME
          OMPT2 = -OMP*(TIME+TAU2)
           OMPT4 = -OMP*(TIME+TAU)
           
C RUNGE-KUTTA-4 METHOD =================>>>

      DO 600 N = 1,NC

C VARIABLES AT PREVIOUS TIMESTEP 
C (WITH COLLISIONS TAKEN INTO ACCOUNT)
        XS = XCL(N)
         YS = YCL(N)
        
        US = UCL(N)
         VS = VCL(N)
         
C STEP 1 ------->
        X = XS
         Y = YS
         
        U = US
         V = VS
        
        X2 = X*X
         Y2 = Y*Y
          R2 = X2+Y2
           R = DSQRT(R2)
            R3 = R2*R

        PK1 = U
         PL1 = V
        
        IF (R .LE. 1.D-8) THEN      ! SKIP CALC WITH 1/R (FORCE = 0)
           PM1 = 0.D0
            PN1 = 0.D0
           GOTO 611
        ENDIF
        
C CENTRAL --->
        OM = ROMEGA(R)/R
         OM2 = OM*OM
         
         PM1 = -OM2*X
          PN1 = -OM2*Y
C BAR --->

        IF (R.LE.RZIR) THEN

          I = IDINT(R/DRG)
          
          IF (I.EQ.0) THEN
               I=1
          ENDIF
          
           CUTI=CUTB(I)
       CUTI=CUTI+(CUTB(I+1)-CUTI)*(R/DRG-DBLE(I))
          
          XR = X/R
           YR = Y/R
        
          XR2 = X*X/R2
           YR2 = Y*Y/R2
            XRYR2 = 2.D0*X*Y/R2
         
          IF (R .GT. RB) THEN                 ! OUTER REGION
            ABAR = -AMPB*CUTI*RB3/R3
            DABAR = -3.D0*ABAR/R
          ELSE                                ! INNER REGION
            ABAR = AMPB*CUTI*(R3/RB3-2.D0)
            DABAR = 3.D0*(ABAR+2.D0)/R
          ENDIF
        
          BR0 = 1.D0 - 2.D0*YR2
           BR1 = BR0*CSHIB + XRYR2*SNHIB
            BR2 = BR0*SNHIB - XRYR2*CSHIB
         
          DABAR = -DABAR*BR1
           ABAR = ABAR*MSP*BR2/R
        
          FBX = DABAR*XR + ABAR*YR
           FBY = DABAR*YR - ABAR*XR

          FBX = FBX*CAMPL1
           FBY = FBY*CAMPL1
        
          PM1 = PM1 + FBX
           PN1 = PN1 + FBY
        ENDIF                 ! END BAR
C SPIRAL --->
        IF (R.GT.RZIL.AND.R.LE.RZOR) THEN
           I=IDINT(R/DRG)
           
           IF (I.EQ.0) THEN
               I=1
           ENDIF
           
            CUTI=CUT(I)
          CUTI=CUTI+(CUT(I+1)-CUTI)*(R/DRG-DBLE(I))
        
            AMT=AMT1*CUTI/(R2*R2)
        
C BELOW FOR CONSTANT PITCH ANGLE ---->

            HI=MSP*(CTGI*DLOG(R/R0S)+OMPT)
             SNHI=DSIN(HI)
              CSHI=DCOS(HI)
        
            FS2=AMT*((X2-Y2)*SNHI-2.D0*X*Y*CSHI)
        
            FSX=FS2*(CTGI*X+Y)
             FSY=FS2*(CTGI*Y-X)
        
            PM1=PM1 + FSX
             PN1=PN1 + FSY
      ENDIF                 ! END SPIRAL
        
611     CONTINUE

C STEP 2 ------->
        X = XS + PK1*TAU2
         Y = YS + PL1*TAU2
          U = US + PM1*TAU2
           V = VS + PN1*TAU2
        
        X2 = X*X
         Y2 = Y*Y
          R2 = X2+Y2
           R = DSQRT(R2)
            R3 = R2*R

        PK2 = U
         PL2 = V
        
        IF (R .LE. 1.D-8) THEN      ! SKIP CALC WITH 1/R (FORCE = 0)
           PM2 = 0.D0
            PN2 = 0.D0
           GOTO 612
        ENDIF
        
C CENTRAL --->
        OM = ROMEGA(R)/R
         OM2 = OM*OM
        
         PM2 = -OM2*X
          PN2 = -OM2*Y
C BAR --->
        IF (R.LE.RZIR) THEN

          I = IDINT(R/DRG)
          
          IF (I.EQ.0) THEN
               I=1
          ENDIF
          
           CUTI=CUTB(I)
       CUTI=CUTI+(CUTB(I+1)-CUTI)*(R/DRG-DBLE(I))
       
          XR = X/R
           YR = Y/R
        
          XR2 = X*X/R2
           YR2 = Y*Y/R2
            XRYR2 = 2.D0*X*Y/R2
         
          IF (R .GT. RB) THEN                 ! OUTER REGION
            ABAR = -AMPB*CUTI*RB3/R3
             DABAR = -3.D0*ABAR/R
          ELSE                                ! INNER REGION
            ABAR = AMPB*CUTI*(R3/RB3-2.D0)
             DABAR = 3.D0*(ABAR+2.D0)/R
          ENDIF
        
          BR0 = 1.D0 - 2.D0*YR2
           BR1 = BR0*CSHIB2 + XRYR2*SNHIB2
            BR2 = BR0*SNHIB2 - XRYR2*CSHIB2
        
          DABAR = -DABAR*BR1
           ABAR = ABAR*MSP*BR2/R
        
          FBX = DABAR*XR + ABAR*YR
           FBY = DABAR*YR - ABAR*XR
         
          FBX = FBX*CAMPL2
           FBY = FBY*CAMPL2
        
          PM2 = PM2 + FBX
           PN2 = PN2 + FBY
        ENDIF                  ! END BAR
C SPIRAL --->
        IF (R.GT.RZIL.AND.R.LE.RZOR) THEN
           I=IDINT(R/DRG)
           
           IF (I.EQ.0) THEN
               I=1
           ENDIF
           
            CUTI=CUT(I)
          CUTI=CUTI+(CUT(I+1)-CUTI)*(R/DRG-DBLE(I))
        
        AMT=AMT2*CUTI/(R2*R2)
        
C BELOW FOR CONSTANT PITCH ANGLE ---->

        HI=MSP*(CTGI*DLOG(R/R0S)+OMPT2)
         SNHI=DSIN(HI)
          CSHI=DCOS(HI)
        
        FS2=AMT*((X2-Y2)*SNHI-2.D0*X*Y*CSHI)
        
        FSX=FS2*(CTGI*X+Y)
         FSY=FS2*(CTGI*Y-X)
        
        PM2=PM2 + FSX
         PN2=PN2 + FSY
       ENDIF                   ! END SPIRAL
         
612     CONTINUE

C STEP 3 ------->
        X = XS + PK2*TAU2
         Y = YS + PL2*TAU2
          U = US + PM2*TAU2
           V = VS + PN2*TAU2
        
        X2 = X*X
         Y2 = Y*Y
          R2 = X2+Y2
           R = DSQRT(R2)
            R3 = R2*R

        PK3 = U
         PL3 = V
        
        IF (R .LE. 1.D-8) THEN      ! SKIP CALC WITH 1/R (FORCE = 0)
           PM3 = 0.D0
            PN3 = 0.D0
           GOTO 613
        ENDIF
        
C CENTRAL --->
        OM = ROMEGA(R)/R
         OM2 = OM*OM
         
         PM3 = -OM2*X
          PN3 = -OM2*Y
C BAR --->
        IF (R.LE.RZIR) THEN

          I = IDINT(R/DRG)
          
          IF (I.EQ.0) THEN
               I=1
          ENDIF
          
           CUTI=CUTB(I)
       CUTI=CUTI+(CUTB(I+1)-CUTI)*(R/DRG-DBLE(I))
        
          XR = X/R
           YR = Y/R
        
          XR2 = X*X/R2
           YR2 = Y*Y/R2
            XRYR2 = 2.D0*X*Y/R2
         
          IF (R .GT. RB) THEN                 ! OUTER REGION
            ABAR = -AMPB*CUTI*RB3/R3
             DABAR = -3.D0*ABAR/R
          ELSE                                ! INNER REGION
            ABAR = AMPB*CUTI*(R3/RB3-2.D0)
             DABAR = 3.D0*(ABAR+2.D0)/R
          ENDIF
        
          BR0 = 1.D0 - 2.D0*YR2
           BR1 = BR0*CSHIB2 + XRYR2*SNHIB2
            BR2 = BR0*SNHIB2 - XRYR2*CSHIB2
        
          DABAR = -DABAR*BR1
           ABAR = ABAR*MSP*BR2/R
        
          FBX = DABAR*XR + ABAR*YR
           FBY = DABAR*YR - ABAR*XR

          FBX = FBX*CAMPL2
           FBY = FBY*CAMPL2
        
          PM3 = PM3 + FBX
           PN3 = PN3 + FBY
        ENDIF                    ! END BAR
C SPIRAL --->
        IF (R.GT.RZIL.AND.R.LE.RZOR) THEN
           I=IDINT(R/DRG)
           
           IF (I.EQ.0) THEN
               I=1
           ENDIF
           
            CUTI=CUT(I)
          CUTI=CUTI+(CUT(I+1)-CUTI)*(R/DRG-DBLE(I))
                
        AMT=AMT2*CUTI/(R2*R2)
        
C BELOW FOR CONSTANT PITCH ANGLE ---->

        HI=MSP*(CTGI*DLOG(R/R0S)+OMPT2)
         SNHI=DSIN(HI)
          CSHI=DCOS(HI)
        
        FS2=AMT*((X2-Y2)*SNHI-2.D0*X*Y*CSHI)
        
        FSX=FS2*(CTGI*X+Y)
         FSY=FS2*(CTGI*Y-X)
        
        PM3=PM3 + FSX
         PN3=PN3 + FSY
        ENDIF                 ! END SPIRAL
        
        
613     CONTINUE

C STEP 4 ------->
        X = XS + PK3*TAU
         Y = YS + PL3*TAU
          U = US + PM3*TAU
           V = VS + PN3*TAU
        
        X2 = X*X
         Y2 = Y*Y
          R2 = X2+Y2
           R = DSQRT(R2)
            R3 = R2*R

        PK4 = U
         PL4 = V
        
        IF (R .LE. 1.D-8) THEN      ! SKIP CALC WITH 1/R (FORCE = 0)
           PM4 = 0.D0
            PN4 = 0.D0
           GOTO 614
        ENDIF
        
C CENTRAL --->
        OM = ROMEGA(R)/R
         OM2 = OM*OM
        
         PM4 = -OM2*X
          PN4 = -OM2*Y
C BAR --->
        IF (R.LE.RZIR) THEN

          I = IDINT(R/DRG)
          
          IF (I.EQ.0) THEN
               I=1
          ENDIF
          
           CUTI=CUTB(I)
       CUTI=CUTI+(CUTB(I+1)-CUTI)*(R/DRG-DBLE(I))
       
          XR = X/R
           YR = Y/R
         
          XR2 = X*X/R2
           YR2 = Y*Y/R2
            XRYR2 = 2.D0*X*Y/R2
         
          IF (R .GT. RB) THEN                 ! OUTER REGION
            ABAR = -AMPB*CUTI*RB3/R3
             DABAR = -3.D0*ABAR/R
          ELSE                                ! INNER REGION
            ABAR = AMPB*CUTI*(R3/RB3-2.D0)
             DABAR = 3.D0*(ABAR+2.D0)/R
          ENDIF
        
          BR0 = 1.D0 - 2.D0*YR2
           BR1 = BR0*CSHIB4 + XRYR2*SNHIB4
            BR2 = BR0*SNHIB4 - XRYR2*CSHIB4
        
          DABAR = -DABAR*BR1
           ABAR = ABAR*MSP*BR2/R
        
          FBX = DABAR*XR + ABAR*YR
           FBY = DABAR*YR - ABAR*XR
         
          FBX = FBX*CAMPL3
           FBY = FBY*CAMPL3
        
          PM4 = PM4 + FBX
           PN4 = PN4 + FBY
        ENDIF                ! END BAR
C SPIRAL --->
        IF (R.GT.RZIL.AND.R.LE.RZOR) THEN
           I=IDINT(R/DRG)
           
           IF (I.EQ.0) THEN
               I=1
           ENDIF
           
            CUTI=CUT(I)
          CUTI=CUTI+(CUT(I+1)-CUTI)*(R/DRG-DBLE(I))
                
        AMT=AMT3*CUTI/(R2*R2)
        
C BELOW FOR CONSTANT PITCH ANGLE ---->

        HI=MSP*(CTGI*DLOG(R/R0S)+OMPT4)
         SNHI=DSIN(HI)
          CSHI=DCOS(HI)
        
        FS2=AMT*((X2-Y2)*SNHI-2.D0*X*Y*CSHI)
        
        FSX=FS2*(CTGI*X+Y)
         FSY=FS2*(CTGI*Y-X)
        
        PM4=PM4 + FSX
         PN4=PN4+ FSY
        ENDIF                 ! END SPIRAL
        
614     CONTINUE

C SUMMARY --------------------->
        XCL(N)=XS+TAU6*(PK1+2.D0*(PK2+PK3)+PK4)
        YCL(N)=YS+TAU6*(PL1+2.D0*(PL2+PL3)+PL4)
        
        UCL(N)=US+TAU6*(PM1+2.D0*(PM2+PM3)+PM4)
        VCL(N)=VS+TAU6*(PN1+2.D0*(PN2+PN3)+PN4)
        
600   CONTINUE
C END RUNGE-KUTTA <<<====================

        TIME = TIME + TAU


C WRITE BLOCK FOR CLOUDS =====>        
        KWRT = KWRT + 1
        
        IF (KWRT.EQ.KREC) THEN
           KWRT = 0
            NREC = NREC + 1
           
           WRITE (25,REC=NREC) XCL,YCL
           
           WRITE (16,REC=NREC) UCL,VCL
           
CCC TEST WRITE        
c           DO 345 I=1,NC
c            XC = XCL(I)
c            YC = YCL(I)
c            WRITE (27,27) XC,YC
c27           FORMAT (1X,2(F8.3,2X))
c345       CONTINUE
CCC END TEST WRITE
        ENDIF
C END WRITE BLOCK <<<----------

300   CONTINUE
C END TIME CYCLE
        WRITE (28,'(I7,A)') NCOLL, ' TOTAL NUM COLL'
      
        WRITE (*,51) NREC,TIME/PER,TIME
51      FORMAT (1X,'LAST NREC=',I5,2X,'TIME/PER=', F6.2,2X,
     ,         'TIME=',F8.4,1x,'Gyrs')
      
      WRITE (*,*) '#####  FINISH  #####'
      READ (*,*)
      STOP
      END
