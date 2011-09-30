

      double precision function x3j(j1,j2,j3,m1,m2,m3)
c uses [theolib.bibl]angalg.for
      implicit double precision (a-h,o-z)
      h=cleb(j1,m1,j2,m2,j3,-m3)*(-1)**((j1-j2-m3)/2)/dsqrt(j3+1d0)
      x3j=h
      return
      end
      double precision function x6j(l1,l2,l3,j1,j2,j3)
c uses [theolib.bibl]angalg.for
      implicit double precision (a-h,o-z)
      h=racah(l1,l2,j2,j1,l3,j3)
      x6j=h*(-1)**((l1+l2+j1+j2)/2)
      return
      end
      double precision function x9j(j1,j2,j3,j4,j5,j6,j7,j8,j9)
c uses [theolib.bibl]angalg.for
      implicit double precision (a-h,o-z)
      h=coef9j(j1,j2,j4,j5,j3,j6,j7,j8,j9)
      x9j=h
      return
      end
      SUBROUTINE TH_FACINIT
      implicit double precision (a-h,o-z)
C ...  SET UP LOG OF FACTORIALS
      PARAMETER (LFACTC=200)
      LOGICAL FIRST
      COMMON / LOGFAC / FIRST,LFACT,FACLOG(LFACTC)
      DATA FIRST/.TRUE./,LFACT/LFACTC/
      FIRST=.FALSE.
      FACLOG(1)=0d00
      FACLOG(2)=0d00
      FN=1d00
      DO 10 I=3,LFACTC
      FN=FN+1d00
      FACLOG(I)=FACLOG(I-1)+DLOG(FN)
   10 CONTINUE
      RETURN
      END
C
C ***********************  CLEB, CLEBI, CLEB2I ************************
C
      double precision FUNCTION CLEBR(A,B,C,D,E,F)
      implicit double precision (a-h,o-z)
C
C      ARGUMENTS ARE REAL AND OF TRUE VALUE; J1,M1,J2,M2,J3,M3
C
      COMMON / LOGFAC / FIRST,LFACT,FACLOG(1)
      LOGICAL FIRST
CCCCCC      IA=2*J1,ID=2*M1 ETC.(J1 IS OF TRUE VALUE)
      IA=NINT(2d0*A)
      IB=NINT(2d0*C)
      IC=NINT(2d0*E)
      ID=NINT(2d0*B)
      IE=NINT(2d0*D)
      IF=NINT(2d0*F)
      GOTO 7000
C ...............  CLEBI  ...........................
C
C      ARGUMENTS ARE INTEGER AND OF TRUE VALUE
C
      ENTRY CLEBI(LL1,LM1,LL2,LM2,LL3,LM3)
      IA=2*LL1
      IB=2*LL2
      IC=2*LL3
      ID=2*LM1
      IE=2*LM2
      IF=2*LM3
      GOTO 7000
C ..............  CLEB  ..............................
C
C      ARGUMENTS ARE INTEGER AND REPRESENT TWICE THE REAL VALUE
C
      ENTRY CLEB(I2J1,I2M1,I2J2,I2M2,I2J3,I2M3)
      IA=I2J1
      IB=I2J2
      IC=I2J3
      ID=I2M1
      IE=I2M2
      IF=I2M3
 7000 IF (FIRST) CALL TH_FACINIT
      RAC=0d00
      IF(ID+IE-IF) 1000,105,1000
  105 K1=IA+IB+IC
      IF((-1)**K1) 1000,110,110
  110 K1=IA+IB-IC
      K2=IC+IA-IB
      K3=IB+IC-IA
      K4=IA-IABS (IB-IC)
      K5=IB-IABS (IC-IA)
      K6=IC-IABS (IA-IB)
      K7= MIN0 (K1,K2,K3,K4,K5,K6)
      IF(K7) 1000,120,120
  120 IF((-1)**(IA+ID)) 1000,1000,130
  130 IF((-1)**(IB+IE)) 1000,1000,140
  140 IF((-1)**(IC+IF)) 1000,1000,150
  150 IF(IA-IABS (ID)) 1000,152,152
  152 IF(IB-IABS (IE)) 1000,154,154
  154 IF(IC-IABS (IF)) 1000,160,160
  160 SIGNFC=1d00
      IAM=IA
      IBM=IB
      ICM=IC
      IDM=ID
      IEM=IE
      IFM=IF
      IF(IA-IB) 210,220,220
  210 IF(IA-IC) 215,225,225
  215 IT=IA
      IA=IB
      IB=IT
      IT=ID
      ID=IE
      IE=IT
      SIGNFC=(-1d00)**((IA+IB-IC)/2)
      GO TO 235
  220 IF(IC-IB) 225,235,235
  225 IT=IC
      IC=IB
      IB=IT
      IT=IF
      IF=-IE
      IE=-IT
      FIBM=IBM+1
      FICM=ICM+1
      SIGNFC=(-1d0)**((IAM-IDM)/2)*DSQRT (FICM/FIBM)
  235 IF(IB) 237,236,237
  236 RAC=SIGNFC
      GO TO 900
  237 IF(IE) 250,250,240
  240 SIGNFC=SIGNFC*((-1d00)**((IA+IB-IC)/2))
      ID=-ID
      IE=-IE
      IF=-IF
  250 FC2=IC+1
      IABCP=(IA+IB+IC)/2+1
      IABC=IABCP-IC
      ICAB=IABCP-IB
      IBCA=IABCP-IA
      IAPD=(IA+ID)/2+1
      IAMD=IAPD-ID
      IBPE=(IB+IE)/2+1
      IBME=IBPE-IE
      ICPF=(IC+IF)/2+1
      ICMF=ICPF-IF
      SQFCLG=0.5d0*(dLOG(FC2)-FACLOG(IABCP+1)
     1      +FACLOG(IABC)+FACLOG(ICAB)+FACLOG(IBCA)
     2      +FACLOG(IAPD)+FACLOG(IAMD)+FACLOG(IBPE)
     3      +FACLOG(IBME)+FACLOG(ICPF)+FACLOG(ICMF))
      NZMIC2=(IB-IC-ID)/2
      NZMIC3=(IA-IC+IE)/2
      NZMI= MAX0 (0,NZMIC2,NZMIC3)+1
      NZMX= MIN0 (IABC,IAMD,IBPE)
      IF(NZMI-NZMX) 310,310,900
  310 SS=0d00
      S1=(-1d00)**(NZMI-1)
      DO 400 NZ=NZMI,NZMX
      NZM1=NZ-1
      NZT1=IABC-NZM1
      NZT2=IAMD-NZM1
      NZT3=IBPE-NZM1
      NZT4=NZ-NZMIC2
      NZT5=NZ-NZMIC3
      TERMLG=SQFCLG-FACLOG(NZ)-FACLOG(NZT1)-FACLOG(NZT2)
     1           -FACLOG(NZT3)-FACLOG(NZT4)-FACLOG(NZT5)
      SSTERM=S1*EXP (TERMLG)
      SS=SS+SSTERM
  400 S1=-S1
      RAC=SIGNFC*SS
  900 IA=IAM
      IB=IBM
      IC=ICM
      ID=IDM
      IE=IEM
      IF=IFM
 1000 CLEB=RAC
      RETURN
      END
      double precision FUNCTION RACAH(JAD,JBD,JCD,JDD,JED,JFD)
C
C        CALCULATES RACAH COEFFICIENTS
C
C        ORIGINAL SOURCE : UNKNOWN
C        SOURCE : IBA_PROGRAM LIBRARY
C        MODIFIED : MARCH 1982 , OLAF SCHOLTEN
C              RUN TIME OPTIMIZED FOR VAX780 MACHINE
C
C        ENTRIES : RACAH , RACAHI , RACAHR
C            RACAH  : INTEGER ARGUMENTS = 2*J
C            RACAHI : ARGUMENTS = TRUE INTEGER VALUE
C            RACAHR : ARGUMENTS = TRUE REAL VALUE
C        EXTERNAL : TH_FACINIT , GENERATES FACTORIAL TABLE
C
      implicit double precision (a-h,o-z)
      DIMENSION I(16)
      LOGICAL FIRST
      COMMON / LOGFAC / FIRST,LFACT,G(1)
      EQUIVALENCE(I(1),I1),(I(2),I2),(I(3),I3),(I(4),I4),(I(5),I5),
     1 (I(6),I6),(I(7),I7),(I(8),I8),(I(9),I9),(I(10),I10),(I(11),I11),
     2 (I(12),I12),(I(13),I13),(I(14),I14),(I(15),I15),(I(16),I16)
C        MAKE USEFULL COMBINATIONS
      K=JAD+JBD-JED+2
      I1=K/2
      IF((2*I1).NE.K) GOTO 300
      K=JCD+JDD-JED+2
      I4=K/2
      IF((2*I4).NE.K) GOTO 300
      K=JAD+JCD-JFD+2
      I7=K/2
      IF((2*I7).NE.K) GOTO 300
      K=JBD+JDD-JFD+2
      I10=K/2
      IF((2*I10).NE.K) GOTO 300
      I13=I1+JED
      I14=I4+JED
      I15=I7+JFD
      I16=I10+JFD
      I2=I13-JAD
      I3=I13-JBD
      I5=I14-JCD
      I6=I14-JDD
      I8=I15-JAD
      I9=I15-JCD
      I11=I16-JBD
      I12=I16-JDD
C       CHECK TRIANGULAR INEQUALITIES,FIND NO. OF TERMS IN SUM
      N=MIN(I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12)-1
      IF(N) 300,2,2
C       FIND MINIMUM VALUE OF SUMMATION INDEX
    2 IL=MAX(I13,I14,I15,I16)
      IF(MIN(JAD,JBD,JCD,JDD,JED,JFD)) 300,20,1
C    ..............
      ENTRY RACAHI(JA1,JB1,JC1,JD1,JE1,JF1)
C        MAKE USEFULL COMBINATIONS
      I13=JA1+JB1+JE1+1
      I14=JC1+JD1+JE1+1
      I15=JA1+JC1+JF1+1
      I16=JB1+JD1+JF1+1
      I1=I13-JE1*2
      I2=I13-JA1*2
      I3=I13-JB1*2
      I4=I14-JE1*2
      I5=I14-JC1*2
      I6=I14-JD1*2
      I7=I15-JF1*2
      I8=I15-JA1*2
      I9=I15-JC1*2
      I10=I16-JF1*2
      I11=I16-JB1*2
      I12=I16-JD1*2
C       CHECK TRIANGULAR INEQUALITIES,FIND NO. OF TERMS IN SUM
      N=MIN(I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12)-1
      IF(N) 300,4,4
C       FIND MINIMUM VALUE OF SUMMATION INDEX
    4 IL=MAX(I13,I14,I15,I16)
      LMIN=MIN(JA1,JB1,JC1,JD1,JE1,JF1)
      IF(LMIN)300,20,1
C     ............
      ENTRY RACAHR(A,B,C,D,E,F)
C     CONVERT ARGUMENTS TO INTEGER 
      JA=NINT(2d0*A)
      JB=NINT(2d0*B)
      JC=NINT(2d0*C)
      JD=NINT(2d0*D)
      JE=NINT(2d0*E)
      JF=NINT(2d0*F)
C        MAKE USEFULL COMBINATIONS
      K=JA+JB-JE+2
      I1=K/2
      IF((2*I1-K).NE.0) GOTO 300
      K=JC+JD-JE+2
      I4=K/2
      IF((2*I4-K).NE.0) GOTO 300
      K=JA+JC-JF+2
      I7=K/2
      IF((2*I7-K).NE.0) GOTO 300
      K=JB+JD-JF+2
      I10=K/2
      IF((2*I10-K).NE.0) GOTO 300
      I13=I1+JE
      I14=I4+JE
      I15=I7+JF
      I16=I10+JF
      I2=I13-JA
      I3=I13-JB
      I5=I14-JC
      I6=I14-JD
      I8=I15-JA
      I9=I15-JC
      I11=I16-JB
      I12=I16-JD
C       CHECK TRIANGULAR INEQUALITIES,FIND NO. OF TERMS IN SUM
      N=MIN(I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12)-1
      IF(N) 300,3,3
C       FIND MINIMUM VALUE OF SUMMATION INDEX
    3 IL=MAX(I13,I14,I15,I16)
      LMIN=MIN(JA,JB,JC,JD,JE,JF)
      IF(LMIN)300,20,1
C      ------------
    1 IF(FIRST) CALL TH_FACINIT
      IF(IL.GE.LFACT) STOP 'RACAH: LENGTH FACTORIAL TABLE INSUFFICIENT'
      J1=IL-I13+1 
      J2=IL-I14+1 
      J3=IL-I15+1 
      J4=IL-I16+1
      J5=I13+I4-IL 
      J6=I15+I5-IL 
      J7=I16+I6-IL
      PH=1d0
      IF(2*(J5/2).EQ.J5) PH=-1d0
      H=PH*EXP ((G(I1)+G(I2)+G(I3)-G(I13+1)+G(I4)+G(I5)+G(I6)-
     1G(I14+1)+G(I7)+G(I8)+G(I9)-G(I15+1)+G(I10)+G(I11)+G(I12)-G(I16+1))
     2*.5d0+G(IL+1)-G(J1)-G(J2)-G(J3)-G(J4)-G(J5)-G(J6)-G(J7))
      IF(N)300,110,120
C
  110 RACAH=H 
      RETURN
C
  120 S=1d0
      K=N-1
      KL=IL+1
      J5=J5-1
      J6=J6-1
      J7=J7-1
      DO 130 J=1,N   ! K=N-J
      S=1d0-((KL+K)*(J5-K)*(J6-K)*(J7-K))
     & *S/((J1+K)*(J2+K)*(J3+K)*(J4+K))
      K=K-1
  130 CONTINUE  
      RACAH=H*S
      RETURN
C
C      ONE OF THE ARGUMENTS =0
   20 IAD=IL
      IBD=IL
      DO 21 J=13,16
      IF(IAD.LT.I(J)) GOTO 22
      IF(IAD.LT.IBD) IBD=IAD
      IAD=I(J)
      GOTO 21
   22 IF(IBD.GT.I(J)) IBD=I(J)
   21 CONTINUE
      J5=I13+I4-IL 
      PH=1d0
      IF(2*(J5/2).EQ.J5) PH=-1d0
      RACAH=PH/DSQRT(DFLOAT(IAD*IBD))
      RETURN
C
C      IMPOSSIBLE COMBINATION OF ARGUMENTS
  300 RACAH=0d0
      RETURN
      END
C
C **************   C O E F 9 J    ***********************************
C
      double precision FUNCTION COEF9J(J1,J2,J3,J4,J5,J6,J7,J8,J9)

      implicit double precision (a-h,o-z)
C
CCCC  TAMURAS NUMBERING CONVENTION IS USED HERE.
CCCCC THE ARGUMENTS OF COEF9J, WHEN NUMBERED SEQUENTIALLY 1 THROUGH 9,
CCCC    CORRESPOND TO THE ARRAY
CCCC                               1  2  5
CCCC                               3  4  6
CCCC                               7  8  9
C
      DIMENSION LT(9)
CCCCCC
CCCCCC      ALL L9 MUST BE TWICE AS LARGE AS TRUE ARGUMENTS
CCCCCC
C
C               CHANGED FOR THEORY LIBRARY 4/8/82   HK
C
      U9=0d00
      LT(1)=J1
      LT(2)=J2
      LT(3)=J3
      LT(4)=J4
      LT(5)=J5
      LT(6)=J6
      LT(7)=J7
      LT(8)=J8
      LT(9)=J9
      LMIN=LT(1)
      IMIN=1
      DO 20 I=2,9
      IF(LT(I)-LMIN) 15,20,20
   15 LMIN=LT(I)
      IMIN=I
   20 CONTINUE
      KEX=0
      GO TO (110,110,110,110,150,150,170,170,190),IMIN
  110 MM=(IMIN-1)/2+1
      M1=MM+MM-1
      M2=M1+1
      M3=MM+4
      L1=LT(7)
      LT(7)=LT(M1)
      LT(M1)=L1
      L1=LT(8)
      LT(8)=LT(M2)
      LT(M2)=L1
      L1=LT(9)
      LT(9)=LT(M3)
      LT(M3)=L1
      IMIN=IMIN+(7-M1)
      GO TO 175
  150 KEX=1
      M1=7
      M2=8
      M3=IMIN+IMIN-9
      M4=M3+1
      GO TO 180
  170 KEX=1
  175 M1=5
      M2=6
      M3=IMIN-6
      M4=M3+2
  180 L1=LT(M1)
      L1=LT(M1)
      LT(M1)=LT(M3)
      LT(M3)=L1
      L1=LT(M2)
      LT(M2)=LT(M4)
      LT(M4)=L1
      L1=LT(9)
      LT(9)=LT(IMIN)
      LT(IMIN)=L1
  190 IF(LT(9)) 200,200,300
  200 IF(LT(5)-LT(6)) 1000,210,1000
  210 IF(LT(7)-LT(8)) 1000,220,1000
  220 RT=(LT(5)+1)*(LT(7)+1)
      K=(LT(5)+LT(7)-LT(1)-LT(4))/2
      RAC= RACAH(LT(1),LT(2),LT(3),LT(4),LT(5),LT(7))
      PH=1d0
      IF (2*(K/2) .NE. K) PH=-1d0
      U9=(RAC/DSQRT(RT))*PH
      GO TO 370
  300 K1=IABS(LT(2)-LT(7))
      K2=IABS(LT(3)-LT(5))
      K3=IABS(LT(4)-LT(9))
      NMIN=MAX0(K1,K2,K3)
      K1=LT(2)+LT(7)
      K2=LT(3)+LT(5)
      K3=LT(4)+LT(9)
      NMAX=MIN0(K1,K2,K3)
      IF (NMIN-NMAX) 320, 320, 1000
  320 DO 350 N=NMIN,NMAX,2
      W1=N+1
      RAC= RACAH(LT(2),LT(5),LT(7),LT(3),LT(1),N)
      IF (RAC) 321, 350, 321
  321 W1=W1*RAC
      RAC= RACAH(LT(2),LT(4),LT(7),LT(9),LT(8),N)
      IF (RAC) 322, 350, 322
  322 W1=W1*RAC
      RAC= RACAH(LT(3),LT(4),LT(5),LT(9),LT(6),N)
      IF (RAC) 323, 350, 323
  323 U9=U9+W1*RAC
  350 CONTINUE
  370 IF(KEX) 400,1000,400
  400 KP=0
      DO 410 I=1,9
  410 KP=KP+LT(I)
      K=KP/2
      PH=1d0
      IF (2*(K/2) .NE. K) PH=-1d0
      U9=U9*PH
 1000 COEF9J=U9
      RETURN
      END
