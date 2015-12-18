C     *****************************LEGAUS*******************************
C
      SUBROUTINE LEGAUS(XS,XL,N,X,W)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N),W(N)
      ZERO=0.0D0
      HALF=0.5D0
       EIN=1.0D0
      ZWEI=2.0D0
      IF(N) 10,10,20
   10 print*,'N > 0 required'
      RETURN
   20 IF(N-2) 30,40,40
   30 X(1)=ZERO
      W(1)=HALF
      GO TO 140
   40 I=1
      G=-EIN
      IC=(N+1)/2
   50 S=G
      T=EIN
      U=EIN
      V=ZERO
      DO 60 K=2,N
      A=K
      FACT1=(ZWEI*A-EIN)/A
      FACT2=(A-EIN)/A
      P=FACT1*G*S-FACT2*T
      DP=FACT1*(S+G*U)-FACT2*V
      T=S
      S=P
      V=U
   60 U=DP
      SUM=ZERO
      IF(I-1) 90,90,70
   70 IM1=I-1
      DO 80 K=1,IM1
   80 SUM=SUM+EIN/(G-X(K))
   90 TEST=G
      G=G-P/(DP-P*SUM)
      R=DABS(TEST-G)
      IF(R.LT.1.0D-13) GO TO 100
      GO TO 50
  100 R=N
      X(I)=G
      W(I)=ZWEI/R/T/DP
      IF(IC-I) 120,120,110
  110 FIM1=IM1
      G=G-(DP-P*SUM)/((ZWEI*G*DP-A*(A+EIN)*P)/(EIN-G*G)-ZWEI*DP*SUM
     1 -P*SUM**2+FIM1*P)
      I=I+1
      GO TO 50
  120 K0=2*IC-N+2*(N/2)+1
      IC=IC+1
      DO 130 I=IC,N
      K=K0-I
      X(I)=-X(K)
  130 W(I)=W(K)
  140 FACT1=(XL-XS)/ZWEI
      FACT2=(XL+XS)/ZWEI
      DO 150 I=1,N
      W(I)=W(I)*FACT1
  150 X(I)=X(I)*FACT1+FACT2
      RETURN
      END
