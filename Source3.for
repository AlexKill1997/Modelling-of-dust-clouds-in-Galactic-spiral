       PROGRAM D3_AN_8

C ------------- 04.05.2018
C            10.05.2018

        PARAMETER (NC=75000,NR=300,NAN=350, LC=1, NV=400, ICC=200,
     ,     NUMCEL = 30, NUMDR = 50, NUMDIS = 50) !cc - cells count

        IMPLICIT REAL*8 (A-H,O-Z)
        REAL*8 M1,M2,M3,M3A,KAP2,MSP
        REAL*8 RCL(NC),TCL(NC),XCL(NC),YCL(NC),
     ,         UCL(NC),VCL(NC),CUT(NR),
     ,         ZCL(NC),WCL(NC),
     ,         DD(NC)
C     ,         R1KS(0:NAN),X1(0:NAN),Y1(0:NAN),
C     ,         X2(0:NAN),Y2(0:NAN),TET2(0:NAN),
C     ,       ZZ(1000),ZDEN(-100:100),ZDN(-100:100),
C     ,       VFR(0:140,0:450),VFT(0:140,0:450)

         INTEGER IJK(-NUMCEL-2:NUMCEL+1,-NUMCEL-2:NUMCEL+1)
         INTEGER CLIST(NC)
         INTEGER NCLR(NUMDR)
         INTEGER NDD(NUMDR,NUMDIS)
         REAL*8 SCLR(NUMDR)
     
C          INTEGER NIZ(0:NAN),NIJ(0:140,0:450)

       DOUBLE PRECISION
     ,     DSIN,DCOS,DBLE,DEXP,DTAN,DASIN,DABS,DSQRT,
     ,     DLOG,DACOS,DSIGN,DATAN

C        character(len=*), 
C     ,   parameter::path="H:\Astro\CLOUDS\program\data"
        
C       character(len=*), parameter::path="D:\_Zolotarev\2Ddata"
       character(len=*), parameter::path=
     = "C:\solution\Final"
        !"D:\Astro\clouds\22-07-2018_21-50-18\2Ddata\"


        ROmega(xx)=(260.d0*dexp(-xx/150.d0-(12.96d0 /(xx*xx)))+
     +  360.d0*dexp(-xx/3.3d0-0.1d0/xx))



      OPEN (10,FILE=path//'\INPdec.dat')
      OPEN (11,FILE=path//'\GALAC.DAT')
      OPEN (15,FILE=path//'\LASTN.DAT',
     ,            ACCESS='DIRECT',RECL=64)
      OPEN (16,FILE=path//'\UVVEL.DAT',
     ,            ACCESS='DIRECT',RECL=16*(NC))
cv     ,            ACCESS='DIRECT',RECL=4*(NC))
C      OPEN (17,FILE='c:\users\YURI\twarog\WVEL.DAT',
C     ,            ACCESS='DIRECT',RECL=16*NC)
      OPEN (25,FILE=path//'\COORD.DAT',
     ,           ACCESS='DIRECT',RECL=16*(NC))
cv     ,           ACCESS='DIRECT',RECL=4*(NC))

      OPEN (27,FILE=path//'\DISTANCESD.DAT',
     ,           ACCESS='DIRECT',RECL=8*(NC))
cv     ,           ACCESS='DIRECT',RECL=4*(NC))

      OPEN (28,FILE=path//'\DISTANCES.DAT')

      OPEN (29,FILE=PATH//'\DISTAN_DISTRIB.DAT')

      OPEN (37,FILE=PATH//'\R-DIST-DISTRIB.DAT')

      OPEN (26,FILE=path//'\CUT.DAT',
     ,           ACCESS='DIRECT',RECL=8*NR)

C      OPEN (30,FILE='c:\users\YURI\3DTRAJ\CHEM.DAT',
C     ,           ACCESS='DIRECT',RECL=24*NC+8)

C      open (31,file='c:\users\YURI\3DTRAJ\tst_pol.dat')
C     open (32,file='c:\users\YURI\3DTRAJ\traj_JAC.dat')

      OPEN (33,FILE=path//'\X_Y.DAT')
C      OPEN (34,FILE=path//'\TET_Z.DAT')
      OPEN (35,FILE=path//'\DENSR.DAT')
      OPEN (36,FILE=PATH//'\VELOCR.DAT')

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
     ,             RCL0,TCL0,
     ,             ABAR,CBAR,BARM,
     ,             HZCL,HZMAX,ELAST
2     FORMAT(1X,'NPUSK=',I1,2X,'NCOL=',I1,2X,'IDSTR=',I1,2X,'NST3=',I4/
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
         ODSI=1.D0/DSQRT(1.D0+CTGI*CTGI)

       TGI=TGI*COFT

         COTM=MSP*CTGI
           CTGI4=CTGI*0.5
             TGI4=TGI*2.D0
              COTM4=COTM

C         TST0=TST0*PI/180.D0

C           FRV=ELAST*0.5D0
        OPE2=0.5D0*(1.D0+ELAST)

         GMO=G*GMO
         OMPS=OMP*OMP
         OMP2=2.D0*OMP


           R=R0

            OM0=ROmega(r)/r
            OM2 = OM0*OM0

            PER=PI2/OM0
        TAU=PER/TCO
c     TAU=1.D-4

c      READ (15,REC=1) NREC,RLIN,RC,RLOT,CAMPL1,TIME,DRG,TAU,NL
       READ (15,REC=1) NREC,RLIN,RC,RLOT,TIME,DRG,TAU,NL


            hi02=hi22
            hi04=hi44



      WRITE (*,*) 'LAST NREC=',NREC
      WRITE (*,*) 'TIME=',TIME,' CAMPL1=',CAMPL1
      WRITE (*,*) 'INPUT NREC'
      READ  (*,*)  NREC
      WRITE (*,*) 'NREC=',NREC,'  WAIT'



C  AMPLITUDE BELOW FOR CONSTANT PITCH ANGLE --->
            WNAM=MSP*DABS(CTGI)/R0
            AMF=ALF*OM2*R0/WNAM
            AMFM=AMF*MSP


C ------------->>  BELOW RCL & TCL == XCL & YCL





          READ (25,REC=NREC) RCL,TCL !, ZCL

          READ (25,REC=NREC) XCL,YCL

          READ (16,REC=NREC) UCL,VCL !, WCL

            ! MAKE THEM INVERSE FOR TEMPERATURE BRIGHTNESS CALC
          !UCL=-UCL
              !VCL=-VCL
          ! END

          READ (26,REC=1) CUT

c------------2D--------
      ZCL = 0.d0
      WCL = 0.d0
c-----------/2D--------

C      GOTO 195   !100 ! 399

C ------------- PICTURE ON GALACTIC PLANE  ==========>>>>

C ----------> TRANSFORM TO ASTRON COORD SYSTEM ---->>>>
CT       TCL=PI-TCL

C  ---------  NPRINT===>> PRINT EACH "NPRINT" POINT
         NPRINT=0

            NCUR=0
       DO 200 N=1,NC
           NCUR=NCUR+1
         IF (NCUR.LT.NPRINT) THEN
            GOTO 200
         ENDIF

C *********************
c           WRITE (*,*) 'N= ',N,'  NCUR=',NCUR


        X=RCL(N)
        Y=TCL(N)

            WRITE (33,33) X,Y
33    FORMAT (1X,F8.3,2X,F8.3)

c/////////  removed this one
c           WRITE (33,33) -X,-Y


            NCUR=0

200    CONTINUE

C  ---- AVERAGE DENSITY VS R ======>>>>>>

        R1=0.1D0
        R2=20.D0
         NDR=IDINT((R2-R1)/DRAN)

      WRITE (35,35) R1,0.D0

       DO 201 I=0,NDR-1
         RI=R1+DRAN*I
         R2=RI+DRAN
         SQ=PI*(R2*R2-RI*RI)
         DEN=0.D0

         DO 202 N=1,NC
           X=XCL(N)
           Y=YCL(N)
           R = DSQRT(X*X+Y*Y)
           IF (R.GE.RI .AND. R.LT.R2) DEN=DEN+1.D0
202      CONTINUE

         DEN=DEN/SQ

         WRITE (35,35) RI,DEN
         WRITE (35,35) R2,DEN
35       FORMAT (1X,F10.3,2X,F16.4)

201    CONTINUE

       DEN=0.D0
      WRITE (35,35) R2,DEN
      WRITE (*,*) '----- SURF DENS VS R DONE ------'

c      STOP 'SURF DENS'
      
      
C ------- AVERAGE R-VELOC VS R >>>

      DELTR = 0.2D0
      RST = 2.0D0
      NDR = IDINT((RG - RST)/DELTR)

      WRITE (36,36) RST, 0.D0
36    FORMAT (2(1X,F8.3))
      DO 210 I = 0,NDR-1

          R1 = RST + DELTR*DBLE(I)
          R2 = R1 + DELTR

          NCL = 0
          SUMVR = 0.D0

          DO 211 N = 1,NC
             X = XCL(N)
             Y = YCL(N)
             VX = UCL(N)
             VY = VCL(N)
             
             R = DSQRT(X*X+Y*Y)

             IF (R.GE.R1 .AND. R.LT.R2) THEN
                VR = (VX*X + VY*Y)/R
                SUMVR = SUMVR + VR
                NCL = NCL+1
             ENDIF

211       CONTINUE

          IF (NCL.NEQV.0) THEN
          
              VR = SUMVR/DBLE(NCL)
          ENDIF    
          
          WRITE (36,36) R1,VR
          WRITE (36,36) R2,VR
          
210   CONTINUE      

      WRITE (36,36) R2, 0.D0
      WRITE (*,*) '----- VELOC R DONE ------'
      
      
      GOTO 100
C      stop '---------- PICTURE DISTRIB ---------'
C      GOTO 195
C DISTANCE BETWEEN PARTICLES (SIMPLE METHOD)

399   CONTINUE
      RG2 = RG*RG
      DO 199 I = 1,NC
        DR2 = RG2

        IF (MOD(I,10000).EQ.0) PRINT *,I
        
        DO 299 J = 1,NC

          IF (I.EQ.J) GOTO 299

          DX = XCL(I)-XCL(J)
          DY = YCL(I)-YCL(J)

          R2 = DX*DX + DY*DY

          IF (R2 .LT. DR2) DR2 = R2          


299     CONTINUE

        DD(I) = DR2
        
199   CONTINUE

      WRITE (27,REC=1) DD

C      DO 198 I = 1,NC
C        WRITE (28,*) DSQRT(DD(I))
C198   CONTINUE

      WRITE (*,*)' --- DISTANCES DONE --- '
      
C DISTANCES --> METHOD 2
100   CONTINUE
      WRITE (*,*)' --- START TO CALC DISTANCES --- '
      
      IJK = 0
      CLIST = 0
      RG2 = RG*RG
      DELT = RG/DBLE(NUMCEL+1)
      
      DO 101 N = 1,NC
         X = XCL(N)
         Y = YCL(N)
         R2 = X*X+Y*Y

         IF (R2 .GT. RG2) GOTO 101

         I = IDINT(X/DELT)
         J = IDINT(Y/DELT)

         IF (X .LT. 0) I = I-1
         IF (Y .LT. 0) J = J-1

         NIJK = IJK(I,J)

         IF (NIJK .GT. 0) CLIST(N) = NIJK

         IJK(I,J) = N
         
101   CONTINUE

      DO 102 N = 1,NC
         X = XCL(N)
         Y = YCL(N)
         R2 = X*X+Y*Y

         IF (R2 .GT. RG2) GOTO 102

         I1 = IDINT(X/DELT)
         J1 = IDINT(Y/DELT)

         IF (X .LT. 0) I1 = I1-1
         IF (Y .LT. 0) J1 = J1-1

         IS1 = I1-1
         IF1 = I1+1
            
         JS1 = J1-1
         JF1 = J1+1

         DRMIN2 = RG2

         DO 341 I = IS1,IF1

           DO 342 J = JS1,JF1

              N2 = IJK(I,J)

343           CONTINUE

                 IF (N2 .EQ. 0) GOTO 342

                 IF (N .EQ. N2) THEN
                    N2 = CLIST(N2)
                    GOTO 343
                 ENDIF

                 DX = X - XCL(N2)
                 DY = Y - YCL(N2)
                  
                 DR2 = DX*DX+DY*DY

                 IF (DR2 .LT. DRMIN2) THEN
                    DRMIN2  = DR2
                 ENDIF

                 N2 = CLIST(N2)
              GOTO 343
                     
                  
342            CONTINUE
341         CONTINUE

         DD(N) = DRMIN2
         IF (MOD(N,100000).EQ.0) PRINT *,N
102   CONTINUE

      WRITE (*,*)' --- DISTANCES DONE --- '

      WRITE (27,REC=1) DD

C      DO 197 I = 1,NC
C        WRITE (28,*) DSQRT(DD(I))
C197   CONTINUE

C SORT DD TO MAKE A HISTOGRAM

      DMIN = 0.D0
      DMAX = 0.1D0
      DR = 0.001D0
      NDR = IDINT((DMAX-DMIN)/DR)

      WRITE (28,28) DMIN,0
      
      DO 150 I = 0, NDR-1
         R1 = DBLE(I)*DR
         R2 = R1+DR
         R11 = R1*R1
         R22 = R2*R2
         NI = 0
         DO 151 N = 1,NC
            DDN = DD(N)
            IF (DDN.GT.R11 .AND. DDN.LE.R22) NI = NI+1
151      CONTINUE
         WRITE (28,28) R1,NI
         WRITE (28,28) R2,NI
28       FORMAT (1X,F10.5,2X,I6)
150   CONTINUE
      WRITE (28,28) R2,0

      WRITE (*,*) '-- DISTANCES HIST DONE --'      

      STOP 'HIST D'
      
195   CONTINUE

C CALC AVERAGE DISTANCES DISTRIB VS R
C (FOR NON-UNIFORM DENSITY DISTRIBUTION)

      IJK = 0
      CLIST = 0
      RG2 = RG*RG
      DELT = RG/DBLE(NUMCEL+1)
      DR = RG/DBLE(NUMDR)
      DMIN = 0.D0
      DMAX = 0.05D0
      DDIS = (DMAX - DMIN)/DBLE(NUMDIS)

      NUMB = 0
      
      SCLR = 0.D0
      NCLR = 0
      WRITE (29,29) 0.D0,0.D0
      
      WRITE (*,*)' --- START TO CALC DISTANCES DISTRIB VS R--- '
      
      DO 120 N = 1,NC
         X = XCL(N)
         Y = YCL(N)
         R2 = X*X+Y*Y

         IF (R2 .GT. RG2) GOTO 120

         I = IDINT(X/DELT)
         J = IDINT(Y/DELT)

         IF (X .LT. 0) I = I-1
         IF (Y .LT. 0) J = J-1

         NIJK = IJK(I,J)

         IF (NIJK .GT. 0) CLIST(N) = NIJK

         IJK(I,J) = N
         
120   CONTINUE

      DO 125 N = 1,NC
         X = XCL(N)
         Y = YCL(N)
         R2 = X*X+Y*Y

         IF (R2 .GT. RG2) GOTO 125

         I1 = IDINT(X/DELT)
         J1 = IDINT(Y/DELT)

         IF (X .LT. 0) I1 = I1-1
         IF (Y .LT. 0) J1 = J1-1

         IS1 = I1-1
         IF1 = I1+1
            
         JS1 = J1-1
         JF1 = J1+1

         DRMIN2 = RG2

         DO 541 I = IS1,IF1

           DO 542 J = JS1,JF1

              N2 = IJK(I,J)

543           CONTINUE

                 IF (N2 .EQ. 0) GOTO 542

                 IF (N .EQ. N2) THEN
                    N2 = CLIST(N2)
                    GOTO 543
                 ENDIF

                 DX = X - XCL(N2)
                 DY = Y - YCL(N2)
                  
                 DR2 = DX*DX+DY*DY

                 IF (DR2 .LT. DRMIN2) THEN
                    DRMIN2  = DR2
                 ENDIF

                 N2 = CLIST(N2)
              GOTO 543
                     
                  
542            CONTINUE
541         CONTINUE

         DD(N) = DRMIN2

         R = DSQRT(R2)
         DRMIN = DSQRT(DRMIN2)
         IR = IDINT(R/DR) + 1
         SCLR(IR) = SCLR(IR) + DRMIN
         NCLR(IR) = NCLR(IR) + 1

         IDD = IDINT(DRMIN/DDIS) + 1

         IF (IDD .LE. NUMDIS) THEN
            NDD(IR,IDD) = NDD(IR,IDD)+1
         ENDIF

         NUMB = NUMB + 1
         IF (MOD(N,100000).EQ.0) PRINT *,N
125   CONTINUE

      WRITE (*,*) 'NUMB = ',NUMB
      WRITE (*,*)' --- DISTANCES DISTRIB VS R DONE --- '

      SCLR = SCLR/DBLE(NCLR)  !? DOES IT WORKS?

      DO 540 I = 1,NUMDR
         R1 = DBLE(I-1)*DR
         R2 = R1 + DR
         WRITE (29,29) R1,SCLR(I)
         WRITE (29,29) R2,SCLR(I)
29       FORMAT (F10.5,2X,F10.5)         
540   CONTINUE
      WRITE (29,29) R2,0.D0

      DR2 = 0.5D0*DR
      DDIS2 = 0.5D0*DDIS
      DO 551 I = 1,NUMDR
         R = DBLE(I-1)*DR + DR2
         DO 552 J = 1,NUMDIS
            DIS = DBLE(J-1)*DDIS + DDIS2
            WRITE (37,37) R,DIS,NDD(I,J)
37          FORMAT (1X,2(F10.5,2X),I6)
552      CONTINUE
551   CONTINUE

      STOP '##### FINISH ######'
      
      END

