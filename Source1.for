       PROGRAM Dec2_ST

C ------------- 04.05.18 ====>
C ---------10.06.18
cv   26.08.2018

        PARAMETER (NC=75000,NR=300, NBR=100,NBTH=360)
        IMPLICIT REAL*8 (A-H,O-Z)
        REAL*8 M1,M2,M3,M3A,MSP ! ,KAP2
        REAL*8 XCL(NC),YCL(NC),
C       ZCL(NC),
     ,       UCL(NC),VCL(NC),
C        WCL(NC),
     ,       CUT(NR)
C     ,       AR0(1:NR),AR2(1:NR),AR4(1:NR),
C     ,       BR2(1:NR),BR4(1:NR)

       DOUBLE PRECISION URR(97),URT(97),URV(97), ! URZ(97),
     ,          URC(97),
     ,          CRR,CDR,CMR,CRT,CDT,CMT,CRV,CDV,CMV,
     ,     DSIN,DCOS,DBLE,DEXP,DTAN,DASIN,DABS,DSQRT,
     ,     DLOG

cv       REAL*8 BARINT1(NBR,NBTH),BARINT2(NBR,NBTH) ! precomputed integrals

       character(len=*), parameter::path=
     = "C:\solution\Final"


       COMMON /ROTC/ OMP,MSP,RGO,M1,M2,M3,
     ,               B1,A2,B2,A3,PI,PI2,
     ,               B1S,A2S,AB2B,M3A,B2S,A3M,AM32,
     ,               RC,RLIN,RLOT,RIN4,ROT4
       COMMON /BARF/ ABAR,CBAR,BARM,CTGI,RZIL,RZIR,RZOL,RZOR,
     ,               RKI,RKO,DROI,WN0,AMF,CSEC

        ROmega(xx)=260.d0*dexp(-xx/150.d0-(3.6d0 / xx)*(3.6d0 / xx))+
     +             360.d0*dexp(-xx/3.3d0-0.1d0/xx)

C         DVF(XX)=260.D0*dexp(-xx/150.d0-(3.6d0 / xx)*(3.6d0 / xx))*
C     *            ((3.6D0*3.6D0)*2/(XX*XX*XX)-1.D0/150.D0)+
C     +            360.d0*dexp(-XX/3.3d0-0.1d0/XX)*
C     *            (0.1D0/(XX*XX)-1.D0/3.3D0)

C         DKAPS(XX)=2.D0*OMG*(OMG+DVDR)


      OPEN (10,FILE=path//'\INPdec.DAT')

      OPEN (15,FILE=path//'\LASTN.DAT',
     ,             ACCESS='DIRECT',RECL=64)
      OPEN (16,FILE=path//'\UVVEL.DAT',
     ,            ACCESS='DIRECT',RECL=16*(NC))
cv     ,            ACCESS='DIRECT',RECL=4*(NC))

C      OPEN (21,FILE='d:\_YURI\3DTRAJ\X_Y.DAT')

      OPEN (25,FILE=path//'\COORD.DAT',
     ,           ACCESS='DIRECT',RECL=16*(NC))
cv     ,           ACCESS='DIRECT',RECL=4*(NC))
      OPEN (26,FILE=path//'\CUT.DAT',
     ,           ACCESS='DIRECT',RECL=8*NR)

      open (31,file=path//'\CONC.DAT')
C      open (32,file=path//'\barint.dat')
C      open (33,file=path//'\bartest.dat')
c      open (27,file=path//'XY.DAT')



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
     ,             RCL0,TCL0,
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
     /     1X,'ABAR=',F5.2,1X,'CBAR=',F5.2,1X,'BARM=',D8.1/
     /     1X,'HZCL=',F5.3,2X,'HZMAX=',F8.3,2X,'ELAST=',F5.3)

cv        PAUSE 'INP GOOD? ---> ENTER'
         write(*,*) 'INP GOOD? ---> ENTER'
         read(*,*)

C     ND2=ND/2
c       DO 300 I=0,ND
c        X=DCL*DFLOAT(ND2-I)
c         Y=DSQRT(RG*RG-X*X)
c          J=IDINT(Y/DCL)
c           JJ(I)=J
c300   CONTINUE

C======CONSTANTS=======
        CRR=  362436.D0/16777216.D0
        CDR= 7654321.D0/16777216.D0
        CMR=16777213.D0/16777216.D0
        IR97=97
        JR97=33

        CRT=  362436.D0/16777216.D0
        CDT= 7654321.D0/16777216.D0
        CMT=16777213.D0/16777216.D0
        IT97=97
        JT97=33

        CRV=  362436.D0/16777216.D0
        CDV= 7654321.D0/16777216.D0
        CMV=16777213.D0/16777216.D0
        IV97=97
        JV97=33

        CRZ=  362436.D0/16777216.D0
        CDZ= 7654321.D0/16777216.D0
        CMZ=16777213.D0/16777216.D0
        IZ97=97
        JZ97=33

        CRC=  362436.D0/16777216.D0
        CDC= 7654321.D0/16777216.D0
        CMC=16777213.D0/16777216.D0
        IC97=97
        JC97=33

        G=4.446666666667D-6

        PI=2.D0*DASIN(1.D0)
        PI2=2.D0*PI
        PID2=0.5D0*PI

       PA=PA/180.D0*PI
       TGI=DTAN(PA)
       CTGI=1.D0/TGI
C           TGI=TGI*COFT

         COTM=MSP*CTGI
           CTGI4=CTGI*0.5d0
             TGI4=TGI*2.D0
              COTM4=COTM

       HI44=HI44/180.d0*PI
       HI22=HI22/180.d0*PI

cv         TST0=TST0*PI/180.D0

         OMPS=OMP*OMP
         OMP2=2.D0*OMP


c/////////////////////////// Precompute bar potential
C           BARINT1=0.
C            BARINT2=0.
            !call BARPRE(BARINT1,BARINT2)
            !write(*,*) 'Precomputed integrals done'

            !------ test
            !------------------
c            r1=0.01
c            r2=8.
c
c            dr=(r2-r1)/dble(NBR)
c            r=r1
c
c            OM2 = (ROmega(r0s)/r0s)**2.d0
c            AMP=0.5*ALF*OM2*R0*R0/DABS(CTGI)
c
c
c                WNAM=MSP*DABS(CTGI)/R0
c                AMF=ALF*OM2*R0/WNAM
c                AMFM=AMF*MSP
c            write(*,*) 'AMFM=', AMFM
c
c
c            th = 0.d0
c            do i=1,NBR-1
c                r=r1 + dble(i-1)*dr
c                th=CTGI*DLOG(r/R0S)
c
c                ir=int(r/dr)+1
c                barphi=BARINT1(ir,1)+BARINT2(ir,1)
c               barforce= -((BARINT1(ir+1,1)-BARINT1(ir,1))/dr+
c     +          (BARINT2(ir+1,1)-BARINT2(ir,1))/dr)
c
c                xs = r*dcos(th)
c                ys = r*dsin(th)
c                x22=xs**2
c                y22=ys**2
c                HI=0.
c                SNHI=DSIN(HI)
c                CSHI=DCOS(HI)
c
c
c                FS2=AMFM*((X22-Y22)*SNHI-2.D0*XS*YS*CSHI)/r**4
c                FSX=FS2*(CTGI*XS+YS)
c                FSY=FS2*(CTGI*YS-XS)
c
c                rn=r1 + dble(i)*dr
c                thn=CTGI*DLOG(rn/R0S)
c                xf=rn*cos(thn)
c                yf=rn*sin(thn)
c                rf=dsqrt(xf**2+yf**2)
c                spf=(FSX*(xf-xs)/rf+FSY*(yf-ys)/rf)
c
c                barforce=barforce*(xf-xs)/rf
c
c                write(33,*) r,barforce,spf
c            enddo
c
c            write(*,*) 'done bar and spiral pot test'


c            r0=0.1

c            dr=8./dble(NBR)
c            dtheta=2.*PI/dble(NBTH)

c            do i=1,NBR

c              r=r0 + dble(i-1)*dr
c              do j=1,NBTH


c                 th=dble(j-1)*dtheta

c                 x=r*dcos(th)
c                 y=r*dsin(th)

c                 phi=BARINT1(i,j)+BARINT2(i,j)
c                 write(32,*) x,y,phi

c                if(phi>0.) then
c                    write(*,*) r, th, phi
c                endif

c              enddo
c            enddo

c            pause 'pot test done'



c////////////////////////////////////////////////////////

C ------- ROTATION DATA + RESONANCES    ---------------->

          CALL ROT

        WRITE (*,*) 'RC=',RC,' RLIN=',RLIN,'  RLOT=',RLOT
        WRITE (*,*) 'RIN4=',RIN4,'  ROT4=',ROT4

c               stop 'stop resonance ?????'

          R=R0S
            OM0 = (ROmega(r)/r)

            PER=PI2/OM0

        TAU=PER/TCO

c//////////////// set explisitly

        TAU = 5.d-4     ! 5.D-4

        !!!!!!!!!!!!!!
cv        tau = 1.d-6
        write(*,*) 'TAU = ', TAU

            HI02=HI22
            HI04=HI44


C        DROI = ?????  RKI,RKO  FOR BAR ???????

       DRG=RG/DBLE(NR)
      DROI=RLOT-RLIN+2.D0*DRG
      RKI=RLIN-DRG
      RKO=RLOT+DRG



            WN0=COTM/R0S


           AMF=ALF*OM0*R0S/WN0
        WRITE(*,107) RC,AMF,TAU,OM0
107   FORMAT(1X,'RC=',F10.7,3X,'AMF=',F14.8/
     /         1X,2X,'TAU=',f14.11,1X,'GYR',2X,'OM0=',F15.10)

C -------- BAR QUADRUPOLS

C       CALL BAR (NR,NT,RG,AR0,AR2,AR4,BR2,BR4,MSP,RLIN,RLOT)
C         CLC=DSQRT(ABAR*ABAR-CBAR*CBAR)
C         CGRAV=3.D0*G*BARM/(CLC*CLC*CLC)



      RZIL=RLIN-CSEC
      RZIR=RLIN+CSEC


      RZOL=RLOT-CSEC
      RZOR=RLOT+CSEC

C *********************
          WRITE (*,*) 'RZIL=',RZIL,'  RZIR=',RZIR
          WRITE (*,*) 'RZOL=',RZOL,'  RZOR=',RZOR


C --------- CUT-OFF FACTOR ------>

            CUT=0.D0

      DO 190 I=1,NR

        RAI=DBLE(I)*DRG
      ! for test
      dd = 0.d0

      IF (RAI.LE.RZIL .OR. RAI.GE.RZOR+dd) GOTO 190

                  CUT(I)=1.D0

      IF (RAI.GT.RZIL .AND. RAI.LT.RZIR) THEN

        CUT(I)=0.5D0*
     *        (1.D0+DSIN(PID2*(RAI-RLIN)/CSEC))

       ENDIF


      IF (RAI.GT.RZOL+dd .AND. RAI.LT.RZOR+dd) THEN

      CUT(I)=0.5D0*
     *        (1.D0-DSIN(PID2*(RAI-RLOT-dd)/CSEC))

      ENDIF


190   CONTINUE

      ! test no cut
      !CUT = 1.d0
      !end test no cut

C      do i=1,NR
C          write(31,*) DBLE(I)*DRG, cut(i)
C      enddo
C         write (*,*) '== force fin ======'



           NREC=1

C +++++ INIT CLOUD COORD +++++++

            CALL RMARIN (IJR,KLR,URR)
            CALL RMARIN (IJT,KLT,URT)
C            CALL RMARIN (IJZ,KLZ,URZ)
            CALL RMARIN (IJC,KLC,URC)

C -- TIPE DISTR OVER R: IDSTR=1- UNIF, 2-EXPON, 3-MOD EXP, 4-R*EXP(-R)


C       IDSTR = 5 FOR MY PAPER AT 2000

        GOTO (501,502,503,504,505,506) IDSTR
C --  START UNIFORM DENSITY DISTR ---
501     R1=RGI

        NUMB=0
       DO 303 J=1,NC/NCC
         R2=DSQRT(R1*R1+(RG*RG-RGI*RGI)*NCC/NC)

           DO 301 I=1,NCC

              IJ=NCC*(J-1)+I
                CALL RRAND (UNI,IR97,JR97,URR,CRR,CDR,CMR)

              RCLIJ=DSQRT((R2*R2-R1*R1)*UNI+R1*R1)

c//////////////////// for quad test
c                d_size = 1.d0*RG + dcl*5.d0
c                r_x = UNI*2.d0*d_size-d_size


C ------------- EXCLUDING THE CENTRAL PART ----->>>>>>>
C           IF (RCLIJ.LE.3.D0) GOTO 301

               NUMB=NUMB+1
                IJ=NUMB

               CALL RRAND (UNIT,IT97,JT97,URT,CRT,CDT,CMT)

               TCL=PI2*UNIT

               XCL(IJ)=RCLIJ*DCOS(TCL)
                YCL(IJ)=RCLIJ*DSIN(TCL)

!21    FORMAT (1X,F8.3,20X,F8.3)
301       CONTINUE
          R1=R2
303    CONTINUE

        GOTO 410
C ------  END UNIF DENS DISTR ---



C ==== EXPONENT DENS DISTR ===
502     NUMB=0
        CNORM=(1.+RGI/SCALE)*DEXP(-RGI/SCALE)-
     -         (1.+RG/SCALE)*DEXP(-RG/SCALE)
       DO 400 I=1,NC
         CALL RRAND (UNI,IR97,JR97,URR,CRR,CDR,CMR)
         ET=1.
         NIT=1
401      ET1=ET-(DEXP(ET)*((1.+RGI/SCALE)*DEXP(-RGI/SCALE)-
     -          UNI*CNORM)-1.D0-ET)/ET
         IF (DABS(ET1-ET).GT.EPSR) THEN
           NIT=NIT+1
             IF (NIT.GT.NITM) THEN
               WRITE (*,*) 'R-ITER-STOP:  I=',I
               STOP
              END IF
           ET=ET1
           GOTO 401
         END IF
        RCLI=SCALE*ET1

C ------------- EXCLUDING OF THE CENTRAL PART ----->>>>>>>
C     IF (RCLI.LE.3.D0) GOTO 400

            NUMB=NUMB+1

        CALL RRAND (UNIT,IT97,JT97,URT,CRT,CDT,CMT)

        TCL=PI2*UNIT


           XCL(NUMB)=RCLI*DCOS(TCL)
           YCL(NUMB)=RCLI*DSIN(TCL)



400    CONTINUE

        GOTO 410
C ==================== END EXP DISTR


C ===  MOD EXPON =====
503     CONTINUE
       NUMB=0
       CALL SIMDST (RG,NST3,RTH,CR3,CNRM3)
         DO 421 I=1,NC
         CALL RRAND (UNI,IR97,JR97,URR,CRR,CDR,CMR)
CC         ET=RG/2.
      ET=RG*UNI
         NIT=1
422      CALL SIMDST (ET,NST3,RTH,CR3,CINT)
         ET1=ET+(CNRM3*UNI-CINT)*(1.+DEXP((ET-RTH)/CR3))/ET
        IF (DABS(ET1-ET).GT.EPSR) THEN
          NIT=NIT+1
          IF (NIT.GT.NITM) THEN
            WRITE (*,*) 'R-ITER-STOP IN DISTR3; I=',I
            STOP
          END IF
           ET=ET1
           GOTO 422
        END IF
         RCLI=ET1

CCC ------------- EXCLUDING OF THE CENTRAL PART ----->>>>>>>
      IF (RCLI.LE.3.D0) GOTO 421
            NUMB=NUMB+1
        RCL0=RCLI
        CALL RRAND (UNIT,IT97,JT97,URT,CRT,CDT,CMT)
          TCL0=2.d0*PI*UNIT
           XCL(I)=RCLI*DCOS(TCL0)
           YCL(I)=RCLI*DSIN(TCL0)
C     WRITE (21,21) XC0,YC0
C        XC0(I)=RCL(I)*DCOS(TCL(I))
C        YC0(I)=RCL(I)*DSIN(TCL(I))

421      CONTINUE
           GOTO 410
C ---------------------
CC504    NUMB=0
CC       CNORM=2.D0-DEXP(-RG/SCALE)*(2.D0+2.D0*(RG/SCALE)+
CC     +    (RG/SCALE)**2)
CC       DO 450 I=1,NC
CC         CALL RRAND (UNI,IR97,JR97,URR,CRR,CDR,CMR)
CC         ET=1.
CC         NIT=1
CC451    ET1=ET-(DEXP(ET)*(2.D0-UNI*CNORM)-(2.D0+2.D0*ET+ET*ET))/
CC     /             (ET*ET)
CC         IF (DABS(ET1-ET).GT.EPSR) THEN
CC           NIT=NIT+1
CC             IF (NIT.GT.NITM) THEN
CC               WRITE (*,*) 'R-ITER-STOP:  I=',I
CC               STOP
CC             END IF
CC           ET=ET1
CC           GOTO 451
CC         END IF
CC        RI=SCALE*ET1
C ------------- EXCLUDING OF THE CENTRAL PART ----->>>>>>>
C     IF (RCLI.LE.3.D0) GOTO 450
CC
CC          NUMB=NUMB+1
CC           RCL(NUMB)=RI
CC        CALL RRAND (UNIT,IT97,JT97,URT,CRT,CDT,CMT)
CC        TI=PI*UNIT
CC        TCL(NUMB)=TI
CC450    CONTINUE


504   CONTINUE

      NUMB=0

C      R1 = RGI
C          x0 = 12.d0
C          x1 = 20.d0
C          ii=NC
C          do while(ii>0)
C              R2=SQRT(R1*R1+(RG*RG-RGI*RGI)/NC)
C              CALL RRAND (UNI,IR97,JR97,URR,CRR,CDR,CMR)
C              CALL RRAND (UNIT,IT97,JT97,URT,CRT,CDT,CMT)
C              x = DSQRT((R2*R2-R1*R1)*UNI+R1*R1)
C              y = UNIT
C              distr = DistribFunc(x, x0, x1)
C              if(y<=distr) then
C                  CALL RRAND (UNIT,IT97,JT97,URT,CRT,CDT,CMT)
C                  RCLIJ=DSQRT((R2*R2-R1*R1)*UNI+R1*R1)
C                  TCL=PI2*UNIT
C                  XCL(NC-ii+1)=RCLIJ*DCOS(TCL)
C                  YCL(NC-ii+1)=RCLIJ*DSIN(TCL)
C                  ii=ii-1
C                  NUMB=NUMB+1
C              endif
C              R1 = R2
C          enddo
          GOTO 410

          
505       CONTINUE

C --------- FOR DISTR IN PAPER OF 2000 YEAR

c ***********************************
c           initial x&y

            XCL = 30.d0
             YCL = 30.d0

c ****

          NUMB = 0
          DR = 1.0D-6

          REF = 3.d0
          SCALE = 4.d0

          PRINT '(A,F6.3)', 'REF   = ',REF
          PRINT '(A,F6.3)', 'SCALE = ',SCALE

! INTEGRAL TO FIND n0 (SIMPSON FORMULA)

         NSTEP = IDINT((RG - RGI)/DR)
         DR6 = DR/6.D0
         DR2 = DR*0.5D0

         SUM1 = 0.D0
         SUM2 = 0.D0

         DO 199 I = 0,NSTEP-1
           R = RGI + DBLE(I)*DR
           FN = 1.D0+DEXP((R-REF)/SCALE)
           SUM1 = SUM1 + R/FN

           R = R + DR2
           FN = 1.D0+DEXP((R-REF)/SCALE)
           SUM2 = SUM2 + R/FN

           R = R + DR
           FN = 1.D0+DEXP((R-REF)/SCALE)
           SUM1 = SUM1 + R/FN
199      CONTINUE

         SUM2 = SUM2*4.D0
         SUM1 = SUM1+SUM2
         SUM1 = SUM1*DR6

         CN0 = DBLE(NC)/SUM1/PI2

         WRITE (*,*) '-------------------------'
         WRITE (*,*) 'n0 = ', CN0
C ------------------------------------------------

          DR = 0.1d0
          DR2 = DR*0.5D0
          DRS = 1.0D-6
          DRS2 = 0.5D0*DRS
          DRS6 = DRS/6.D0
          NRING = IDINT((RG-RGI)/DR)
          NUMB = 0

          WRITE (*,*) 'NUM RINGS = ', NRING

          DO 111 I = 0,NRING-1
            R1 = RGI + DBLE(I)*DR
            R2 = R1+DR
            RI2 = R1*R1
            SQ = R2*R2 - RI2

C INT n(r) OVER THE RING USING SIMPSON FORMULAE
            NSTEP = IDINT((R2-R1)/DRS)
            SUM1 = 0.D0
            SUM2 = 0.D0

            DO 198 J = 0,NSTEP-1
              R = R1 + DBLE(J)*DRS
              FN = 1.D0+DEXP((R-REF)/SCALE)
              SUM1 = SUM1 + R/FN

              R = R + DRS2
              FN = 1.D0+DEXP((R-REF)/SCALE)
              SUM2 = SUM2 + R/FN

              R = R + DRS
              FN = 1.D0+DEXP((R-REF)/SCALE)
              SUM1 = SUM1 + R/FN
198         CONTINUE
              SUM2 = SUM2*4.D0
              SUM1 = SUM1+SUM2
              SUM1 = SUM1*DRS6

            DNI = PI2*CN0*SUM1
C --------------------------------------
            NI = IDINT(DNI)
C <<<<-----------------------------------
            DO 121 J = 1,NI
              CALL RRAND (UNI,IR97,JR97,URR,CRR,CDR,CMR)
              RCL2 = UNI*SQ + RI2
              RCL = DSQRT(RCL2)

              CALL RRAND (UNIT,IT97,JT97,URT,CRT,CDT,CMT)
              TCL = PI2*UNIT

              NUMB = NUMB + 1

              XCL(NUMB) = RCL*DCOS(TCL)
              YCL(NUMB) = RCL*DSIN(TCL)
                
121         CONTINUE

111       CONTINUE

          WRITE (*,*) 'NUMB = ',NUMB
          WRITE (*,*) 'NC - NUMB = ',NC - NUMB

        GOTO 410

c          nDR = idint((RG-RGI)/dr)
c          r = RGI


c          tmp = 0.d0    ! integral FOR CONCENTRATION

c         do i=1,nDR

c             r = rgi+(dble(i)-0.5d0)*dr
c             tmp = tmp + r/(1.d0+dexp((r-Ref)/SCALE))
c                 IF (R.LE.Ref) TMPIN=TMP

c          enddo

c        TMPIN=TMPIN*DR*PI2
c          TMP=TMP*DR*pi2

c      WRITE (*,*) 'TMP=',TMP
ccv        PAUSE 'TMP ?????'
c          write(*,*) 'TMP ?????'
c          read(*,*)

c             DNCIN=DBLE(NC)*TMPIN/TMP

c          dn0 = dble(NC)/(tmp)

c      write (*,*)

c      WRITE (*,*) 'NUMB PART WITHIN Ref= ',DNCIN

c////////////////////           test curve n
c          r=RGI
c          dr = 0.1d0
c          n = iDInt((RG-RGI)/dr)

c          do i=1,n
c              !write(31,*) r, dn0/(1.d0+exp((r-a)/b))
c              r =RGI + i*dr
c          enddo
c/////////       end  test curve n0

C          write(*,*) 'n0 = ', DN0


c          R1 = RGI
c          R2 = R1+DR
c          k=1
c          sumni = 0

cc*************************


c            DO 1100 i=1,n+1

c             ! compute ni
c             step=dr*0.0001d0
c             m = iDInt((r2-r1)/step)
c             rr = r1
c             tmp=0.d0

c             do j=1,m
c                fa=rr/(1.d0+dexp((rr-Ref)/SCALE))
c                fb=(rr+step)/(1.d0+dexp(((rr+step)-Ref)/SCALE))
c                tmp = tmp + (fa+fb)
c                rr=rr+step
c             enddo

c             sni = PI2*DN0*STEP*0.5D0*tmp
             
c             ni = idint(sni)
cc ***************************
cc            write (*,*) 'ni=',ni

c             ! compute coords for clouds in [r1..r2]

c               sumni=sumni+sni

c       do 1501 j=1,ni

c                CALL RRAND (UNI,IR97,JR97,URR,CRR,CDR,CMR)
c                CALL RRAND (UNIT,IT97,JT97,URT,CRT,CDT,CMT)

c                r = DSQRT((r2*r2-r1*r1)*UNI+r1*r1)
c                theta = PI2*UNIT
c                XCL(k)=r*dcos(theta)
c                YCL(k)=r*dsin(theta)

c                k=k+1
c                NUMB = NUMB + 1  !!!!!!!!!  ???????????????????
c                if(k.EQ.NC+1) then
c                   write(*,*) '... sum ni = ', sumni
c           GOTO 410
c                endif



c1501     CONTINUE

c             r1=r2
c             r2=r2+dr
c1100     CONTINUE


cc***************************
c        write (*,*) 'r1=',r1,'   r2=',r2

c      GOTO 410


C ____________________________

506   CONTINUE


C      do i=1,NC
C          CALL RRAND (UNI,IR97,JR97,URR,CRR,CDR,CMR)
C          UNI=dsqrt(UNI)
C          p1=0.
C          m=-1
C          do j=1,n
C              p2=p1+MDIST(j)
C              if(UNI>p1 .and. UNI<=p2) then
C                  m=j
C                  exit  ! break do cycle
C              endif
C              p1=p2
C          enddo

C          CALL RRAND (UNI,IR97,JR97,URR,CRR,CDR,CMR)
C          CALL RRAND (UNIT,IT97,JT97,URT,CRT,CDT,CMT)

C          r1=RGI + dble(m)*dr
C          r2=r1+dr
C          r = DSQRT((r2*r2-r1*r1)*UNI+r1*r1)

C          th=UNIT*2.*PI

C          XCL(i)=r*dcos(th)
C          YCL(i)=r*dsin(th)
C
C
C          NUMB=NUMB+1
C      enddo
C
C                   k=1
C      do i=1,n*0
C          r1=RGI+dble(i-1)*dr
C          r2=r1+dr
C          do j=1,NDIST(i)
C
C                CALL RRAND (UNI,IR97,JR97,URR,CRR,CDR,CMR)
C                CALL RRAND (UNIT,IT97,JT97,URT,CRT,CDT,CMT)
C
C                r = DSQRT((r2*r2-r1*r1)*0.5+r1*r1)
C                theta = 2.*PI*UNIT
C                XCL(k)=r*dcos(theta)
C                YCL(k)=r*dsin(theta)
C                k=k+1
C          enddo
C      enddo
C
C      GOTO 410




c/////////////////////////////////////

410       CONTINUE

C       WRITE (*,*) '===  NUMB ==',NUMB, ' sumni=',sumni
C      write (*,*) 'k=  ',k
c     PAUSE 'NUMB PART === GOOD? ----> ENTER'


C ----- Z - coordinate --------->


C       DO 471 I=1,NUMB
C472      CALL RRAND(UNI,IZ97,JZ97,URZ,CRZ,CDZ,CMZ)
C          Z1=UNI*2.D0-1.D0
C         CALL RRAND(UNI,IZ97,JZ97,URZ,CRZ,CDZ,CMZ)
C          Z2=UNI*2.D0-1.D0
C           S=Z1*Z1+Z2*Z2
C             IF (S.GT.1.D0) GOTO 472
C          ZX=DSQRT(-2.D0*DLOG(S)/S)



C --========== HZCL = H_DISK --------------------->>>>>>>>
C  CHECK HZCL:  IS IT Z-DISPERSION OR WIDTH ????????


C          ZCL(I)=Z1*ZX*HZCL
C          WCL(I)=Z2*ZX*DISP

c ===  for flat  disk ==
C                    zcl(i)=0.d0
C                    wcl(i)=0.d0
c ******
C471    CONTINUE



C -- CLOUD VELOC BLOCK IN GALAC PLANE ----->>>>>>
C !! VELOC DISP IZ EQUAL IN ALL DIRECTIONS !!!!

        CALL RMARIN (IJV,KLV,URV)


       DO 302 I=1,NUMB

        XI=XCL(I)
         YI=YCL(I)
          R2=XI*XI+YI*YI
           R=DSQRT(R2)

            OMG = (ROmega(r)/r)

333     CONTINUE
           CALL RRAND (UNI,IV97,JV97,URV,CRV,CDV,CMV)
          V1U=UNI*2.D0-1.D0
           CALL RRAND (UNI,IV97,JV97,URV,CRV,CDV,CMV)
          V2V=UNI*2.D0-1.D0
         S=V1U*V1U+V2V*V2V
         IF (S.GT.1.D0) GOTO 333

C ---------------- EXPONENT VELOC DISPERSION
C         UX=DISP*DEXP(0.5D0*(R0S-RI)/SCALE)*
C     *          DSQRT(-2.D0*DLOG(S)/S)

C --- CONSTANT DISP:  C_RAD=C_THETA= DISP AT THE SUN ----
         UX=-2.D0*DLOG(S)/S
         UX=DSQRT(UX)

         U1X=V1U*DISP*UX
         V1Y=V2V*DISP*UX

C         V1Y=V2V*UX*0.6
c  -- above coeff 0.6 follows from observations FOR STARS


c******************** no chaos vel for test
            !u1x = 0
            !v1y = 0
c**********

         UCL(I)=OMG*YI+U1X
          VCL(I)=-OMG*XI+V1Y

 302    CONTINUE

C *********************

C          WRITE (*,*) 'N=',N

       DO 222 I=1,200
         R=DBLE(I)*0.1D0
           CONC=CN0/(1.D0+DEXP((R-REF)/SCALE))
              WRITE (31,31) R,CONC
31              FORMAT (1X,F5.2,2X,F8.2)

222    CONTINUE


        TIME=0.


c          write(*,*) 'Before'
c          write(*,*) 'XCL=',XCL,'YCL=',YCL
          !test for collision. 2 particles. DCL = 0.02
c          XCL(1) = 10.0
c          YCL(1) = 0.0 + 0.01
c          XCL(2) = 10.0 + 0.01
c          YCL(2) = 0.0

c          UCL(1) = 1.0
c          VCL(1) = 0.0
c          UCL(2) = 0.0
c          VCL(2) = 0.0
          !end collision test
c          write(*,*) 'XCL=',XCL,'YCL=',YCL
c          write(*,*) 'UCL=',UCL,'VCL=',VCL

c         for movement test
          !XCL(1) = 8.0d0
          !YCL(1) = 0.0d0

            !WRITE(*,*) ucl(1), vcl(1)
c          UCL(1) = 0.0d0
c          VCL(1) = OMP*XCL(1)

          NL = 1  ! for no wrning ?????
       
c          WRITE (15,REC=1) NREC,RLIN,RC,RLOT,TIME,DRG,TAU
          WRITE (15,REC=1) NREC,RLIN,RC,RLOT,TIME,DRG,TAU,NL
          WRITE (16,REC=1) UCL,VCL !,WCL
          WRITE (25,REC=1) XCL,YCL !,ZCL
          WRITE (26,REC=1) CUT

        !  write (*,*) 'NREC,RLIN,RC,RLOT,TIME,DRG,TAU,NL'
        !  write (*,*) NREC,RLIN,RC,RLOT,TIME,DRG,TAU,NL



        WRITE (*,51) NREC,TIME/PER,TIME
51      FORMAT (1X,'LAST NREC=',I5,2X,'TIME/PER=', F6.2,2X,
     ,         'TIME=',F8.4,1x,'Gyrs')


CCC TEST WRITE        
c        do 345 i=1,nc

c          Xw = XCL(I)
c           Yw = YCL(I)
c            uw = UCL(i)
c             vw = vcl(i)
c          WRITE (27,27) Xw,Yw,uw,vw
cc          print*, i
c27         FORMAT (1X,4(F8.3,2X)) 
c345     continue

CCC END TEST WRITE


        WRITE (*,*) '###### FINISH #######'
c        read(*,*)
        STOP
        END


c/////////////////////////////////////////////////////////

        SUBROUTINE RRAND (UNI,I97,J97,UU,CC,CCD,CCM)
        DOUBLE PRECISION UU(97),CC,CCD,CCM,UNI

          UNI=UU(I97)-UU(J97)
          IF (UNI.LT.0.D0) UNI=UNI+1.D0
          UU(I97)=UNI
          I97=I97-1
          IF (I97.EQ.0) I97=97
          J97=J97-1
          IF (J97.EQ.0) J97=97
          CC=CC-CCD
          IF (CC.LT.0.D0) CC=CC+CCM
          UNI=UNI-CC
          IF (UNI.LT.0.D0) UNI=UNI+1.D0

          RETURN
          END

        SUBROUTINE RMARIN (IJ,KL,U)
        DOUBLE PRECISION U(97),S,T

          I=MOD(IJ/177,177)+2
          J=MOD(IJ,177)+2
          K=MOD(KL/169,178)+1
          L=MOD(KL,169)

        DO 2 II=1,97
          S=0.D0
          T=0.5D0
         DO 3 JJ=1,24
           M=MOD(MOD(I*J,179)*K,179)
           I=J
           J=K
           K=M
           L=MOD(53*L+1,169)
         IF (MOD(L*M,64).GE.32) THEN
             S=S+T
         END IF
           T=0.5D0*T
3        CONTINUE
         U(II)=S
2      CONTINUE

       RETURN
       END

C =====================================================
       SUBROUTINE SIMDST (RO3,NS3,RTH3,CRC3,CN3)
         IMPLICIT REAL*8 (A-H,O-Z)
         REAL*8 DBLE,DEXP
         CN3=0.
         STR3=RO3/DBLE(NS3)
         ST2=0.5*STR3
         ABCS=0.
         SY0=0.
       DO 420 I=1,NS3
         ABCS=ABCS+ST2
         SY1=ABCS/(1.+DEXP((ABCS-RTH3)/CRC3))
         ABCS=ABCS+ST2
         SY2=ABCS/(1.+DEXP((ABCS-RTH3)/CRC3))
         CN3=CN3+SY0+4.*SY1+SY2
         SY0=SY2
420    CONTINUE
         CN3=ST2/3.*CN3

         RETURN
        END

C ===========================================================

       SUBROUTINE ROT
        IMPLICIT REAL*8 (A-H,O-Z)
        REAL*8 M1,M2,M3,M3A,KAP2,MSP

       COMMON /ROTC/ OMP,MSP,RGO,M1,M2,M3,
     ,               B1,A2,B2,A3,PI,PI2,
     ,               B1S,A2S,AB2B,M3A,B2S,A3M,AM32,
     ,               RC,RLIN,RLOT,RIN4,ROT4

       DOUBLE PRECISION DEXP, DSQRT !,DSIN ,DASIN,DCOS,DBLE,DTAN,
cv     ,                   DLOG,DABS
     
       character(len=*), parameter::path=
     = "C:\solution\Final"

        !"/home/roman/Astro/clouds/2Ddata/"

        ROmega(xx)=260.d0*dexp(-xx/150.d0-(3.6d0 / xx)*(3.6d0 / xx))+
     +             360.d0*dexp(-xx/3.3d0-0.1d0/xx)

         DVF(XX)=260.D0*dexp(-xx/150.d0-(3.6d0 / xx)*(3.6d0 / xx))*
     *            ((3.6D0*3.6D0)*2/(XX*XX*XX)-1.D0/150.D0)+
     +            360.d0*dexp(-XX/3.3d0-0.1d0/XX)*
     *            (0.1D0/(XX*XX)-1.D0/3.3D0)

         DKAPS(XX)=2.D0*OMG*(OMG+DVDR)

      OPEN (11,FILE=path//'\GALAC.DAT')

C======CONSTANTS=======

         OMPS=OMP*OMP
         OMP2=2.D0*OMP


         WRITE (11,42)
42       FORMAT(3X,'R',7X,'OMG',7X,'OMMK',6X,'OMPK')

            IF (OMP.LT.0.1) THEN
                WRITE (*,*) '   OMP TOO SMALL ?????????'
                GOTO 133
            END IF

         R=0.5
         VMAX=0.
144     CONTINUE

           OMG = (ROmega(r)/r)
            DVDR=DVF(R)

           KAP2=DKAPS(R)

c           write(*,*) 'KAP2 = ',KAP2

            VROT=ROMEGA(R)
         IF (VROT.GT.VMAX) VMAX=VROT

         SQKAP=DSQRT(KAP2)/MSP
         OMMK=OMG-SQKAP
         OMPK=OMG+SQKAP

          WRITE (11,40) R,OMG,OMMK,OMPK
40          format (1X,F6.3,2X,F8.3,2X,F8.3,2X,F8.3)

        IF (R.LE.RGO) THEN
             R=R+0.1
             GOTO 144
        END IF
C  -- end of omega ---

C ----- RESONANCES --- MSP = 2   ====>>>>>>>>>>>>>
        DO 166 I=1,3
        K=1
        R=1.
          DR=0.1

           OMG = (ROmega(r)/r)

         DOM=OMG-OMP
                  DVDR=DVF(R)
        IF (I.EQ.1) Y1=DOM

       IF (I.EQ.2) Y1=DOM-DSQRT(DKAPS(R))/MSP

        IF (I.EQ.3) Y1=DOM+DSQRT(DKAPS(R))/MSP

100     R=R+DR

          IF (R.GT.RGO) GOTO 102
          OMG = (ROmega(r)/r)
                  DVDR=DVF(R)

         DOM=OMG-OMP
        IF (I.EQ.1) Y2=DOM
        IF (I.EQ.2) Y2=DOM-DSQRT(DKAPS(R))/MSP
        IF (I.EQ.3) Y2=DOM+DSQRT(DKAPS(R))/MSP
        IF(K.EQ.3) GOTO 102
        YM=Y1*Y2
        Y1=Y2
        IF(YM.GT.0.0) THEN
               GOTO 100
        ELSE
               DR=-0.1*DR
               K=K+1
               GOTO 100
        ENDIF
102   CONTINUE
         IF (I.EQ.1) RC=-0.5*DR+R
         IF (I.EQ.2) RLIN=-0.5*DR+R
         IF (I.EQ.3) RLOT=-0.5*DR+R
166    CONTINUE

         RIN=RLIN
         ROUT=RLOT
C ---- END MSP = 2--------

C ----RESONANCES FOR MSP = 4 =====>>>>>>
        DO 266 I=1,2
        K=1
        R=1.
          DR=0.1
          OMG=ROMEGA(R)/R
                  DVDR=DVF(R)

         DOM=OMG-OMP
        IF (I.EQ.1) Y1=DOM-DSQRT(DKAPS(R))/MSP/2.
        IF (I.EQ.2) Y1=DOM+DSQRT(DKAPS(R))/MSP/2.
190     R=R+DR
          IF (R.GT.RGO) GOTO 192
       OMG=ROMEGA(R)/R

                  DVDR=DVF(R)
         DOM=OMG-OMP
        IF (I.EQ.1) Y2=DOM-DSQRT(DKAPS(R))/MSP/2.
        IF (I.EQ.2) Y2=DOM+DSQRT(DKAPS(R))/MSP/2.
        IF(K.EQ.2) GOTO 192
        YM=Y1*Y2
        Y1=Y2
        IF(YM.GT.0.0) THEN
               GOTO 190
        ELSE
               DR=-0.1*DR
               K=K+1
               GOTO 190
        ENDIF
  192   CONTINUE
         IF (I.EQ.1) RLIN4=-0.5*DR+R
         IF (I.EQ.2) RLOT4=-0.5*DR+R
266    CONTINUE

         RIN4=RLIN4
         ROT4=RLOT4
          WRITE (*,*) 'RIN_2=',RLIN,'  ROT_2=',RLOT
          WRITE (*,*) 'RIN_4=',RLIN4,'  ROT_4=',RLOT4


C ----- 'COROTATION CIRCLE DO NOT WRITE '-----
C       DO 320 I=1,360
C         TEI=PI2/360.*DBLE(I-1)
C         XC=RC*DCOS(TEI)
C         YC=RC*DSIN(TEI)
C320   CONTINUE

C ----- SPIRAL PATTERN ==== DO NOT WRITE-->>>

C        DO 170 I=1,100
C         R=2.+(I-1)*0.1
C         AZ=CTGI*DLOG(R/R0S)
C         XSP1=R*DCOS(AZ)
C         YSP1=R*DSIN(AZ)
C
C         XSP2=R*DCOS(AZ+PI2/MSP)
C         YSP2=R*DSIN(AZ+PI2/MSP)
C
C         XSP3=R*DCOS(AZ+2.*PI2/MSP)
C         YSP3=R*DSIN(AZ+2.*PI2/MSP)
C
C         XSP4=R*DCOS(AZ+3.*PI2/MSP)
C         YSP4=R*DSIN(AZ+3.*PI2/MSP)
C170    CONTINUE

133        WRITE (*,*) '--- FINISH ROTATION DAT ----------'
        RETURN
      END


