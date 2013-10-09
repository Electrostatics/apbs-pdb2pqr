
SUBROUTINE EZFFTF(N, R, AZERO, A, B, WSAVE)
      IMPLICIT REAL*8(A-H,O-Z)
      
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER  :: N
      REAL*8 , INTENT(OUT) :: AZERO
      REAL*8 , INTENT(IN) :: R(*)
      REAL*8 , INTENT(OUT) :: A(*)
      REAL*8 , INTENT(OUT) :: B(*)
      REAL*8  :: WSAVE(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, NS2, NS2M
      REAL*8  :: CF, CFM
!-----------------------------------------------
!
      IF (N - 2 <= 0) THEN
         IF (N - 2 /= 0) THEN
            AZERO = R(1)
            RETURN 
         ENDIF
         AZERO = 0.5D0*(R(1)+R(2))
         A(1)  = 0.5D0*(R(1)-R(2))
         RETURN 
      ENDIF
      WSAVE(:N) = R(:N)
      CALL RFFTF (N, WSAVE, WSAVE(N+1))
      CF    = 2.D0/FLOAT(N)
      CFM   = -CF
      AZERO = 0.5D0*CF*WSAVE(1)
      NS2   = (N + 1)/2
      NS2M  = NS2 - 1
      A(:NS2M) = CF*WSAVE(2:NS2M*2:2)
      B(:NS2M) = CFM*WSAVE(3:NS2M*2+1:2)
      IF (MOD(N,2) == 1) RETURN 
      A(NS2) = 0.5D0*CF*WSAVE(N)
      B(NS2) = 0.D0
      
      RETURN 
      
END SUBROUTINE EZFFTF


SUBROUTINE EZFFTB(N, R, AZERO, A, B, WSAVE)
      IMPLICIT REAL*8(A-H,O-Z)
      
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      REAL*8 , INTENT(IN) :: AZERO
      REAL*8  :: R(*)
      REAL*8 , INTENT(IN) :: A(*)
      REAL*8 , INTENT(IN) :: B(*)
      REAL*8  :: WSAVE(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NS2, I
!-----------------------------------------------
!
      IF (N - 2 <= 0) THEN
         IF (N - 2 /= 0) THEN
            R(1) = AZERO
            RETURN 
         ENDIF
         R(1) = AZERO + A(1)
         R(2) = AZERO - A(1)
         RETURN 
      ENDIF
      NS2 = (N - 1)/2
      R(2:NS2*2:2) = 0.5D0*A(:NS2)
      R(3:NS2*2+1:2) = -0.5D0*B(:NS2)
      R(1) = AZERO
      IF (MOD(N,2) == 0) R(N) = A(NS2+1)
      CALL RFFTB (N, R, WSAVE(N+1))
      
      RETURN 
      
END SUBROUTINE EZFFTB


SUBROUTINE EZFFTI(N, WSAVE)
      
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      REAL*8  :: WSAVE(*)
      
!-----------------------------------------------
!
      IF (N == 1) RETURN 
      CALL EZFFT1 (N, WSAVE(2*N+1), WSAVE(3*N+1))
      RETURN 
END SUBROUTINE EZFFTI


SUBROUTINE EZFFT1(N, WA, IFAC)
      IMPLICIT REAL*8(A-H,O-Z)
      
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(INOUT) :: IFAC(*)
      REAL*8 , INTENT(INOUT) :: WA(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(4) :: NTRYH
      INTEGER :: NL,NF,J,NTRY,NQ,NR,I,IB,IS,NFM1,L1,K1,IP,L2,IDO,IPM,II
      REAL*8  :: TPI, DUM, ARGH, ARG1, CH1, SH1, DCH1, DSH1, CH1H
!-----------------------------------------------
      DATA NTRYH(1), NTRYH(2), NTRYH(3), NTRYH(4)/ 4, 2, 3, 5/ 
      TPI = 8.D0*ATAN(1.D0)
      NL  = N
      NF  = 0
      J   = 0
  101 CONTINUE
      J = J + 1
      IF (J - 4 <= 0) THEN
         NTRY = NTRYH(J)
      ELSE
         NTRY = NTRY + 2
      ENDIF
  104 CONTINUE
      NQ = NL/NTRY
      NR = NL - NTRY*NQ
      IF (NR /= 0) GO TO 101
      NF = NF + 1
      IFAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY == 2) THEN
         IF (NF /= 1) THEN
            IFAC(NF+2:4:(-1)) = IFAC(NF+1:3:(-1))
            IFAC(3) = 2
         ENDIF
      ENDIF
      IF (NL /= 1) GO TO 104
      IFAC(1) = N
      IFAC(2) = NF
      ARGH = TPI/FLOAT(N)
      IS   = 0
      NFM1 = NF - 1
      L1   = 1
      IF (NFM1 == 0) RETURN 
      DO K1 = 1, NFM1
         IP = IFAC(K1+2)
         L2 = L1*IP
         IDO = N/L2
         IPM = IP - 1
         ARG1 = FLOAT(L1)*ARGH
         CH1 = 1.D0
         SH1 = 0.D0
         DCH1 = COS(ARG1)
         DSH1 = SIN(ARG1)
         DO J = 1, IPM
            CH1H = DCH1*CH1 - DSH1*SH1
            SH1 = DCH1*SH1 + DSH1*CH1
            CH1 = CH1H
            I = IS + 2
            WA(I-1) = CH1
            WA(I) = SH1
            IF (IDO >= 5) THEN
               DO II = 5, IDO, 2
                  I = I + 2
                  WA(I-1) = CH1*WA(I-3) - SH1*WA(I-2)
                  WA(I) = CH1*WA(I-2) + SH1*WA(I-3)
               END DO
            ENDIF
            IS = IS + IDO
         END DO
         L1 = L2
      END DO
      RETURN 
END SUBROUTINE EZFFT1


SUBROUTINE COSTI(N, WSAVE)
       IMPLICIT REAL*8(A-H,O-Z)
       
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      REAL*8  :: WSAVE(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NM1, NP1, NS2, K, KC
      REAL*8  :: PI, DUM, DT, FK
!-----------------------------------------------
!
      PI = 4.D0*ATAN(1.0)
      IF (N <= 3) RETURN 
      NM1 = N - 1
      NP1 = N + 1
      NS2 = N/2
      DT  = PI/FLOAT(NM1)
      FK  = 0.D0
      DO K = 2, NS2
         KC = NP1 - K
         FK = FK + 1.D0
         WSAVE(K) = 2.D0*SIN(FK*DT)
         WSAVE(KC) = 2.D0*COS(FK*DT)
      END DO
      CALL RFFTI (NM1, WSAVE(N+1))
      RETURN 
      
END SUBROUTINE COSTI


SUBROUTINE COST(N, X, WSAVE)
       IMPLICIT REAL*8(A-H,O-Z)
       
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      REAL*8  :: X(*)
      REAL*8  :: WSAVE(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NM1, NP1, NS2, K, KC, MODN, I
      REAL*8  :: X1H, X1P3, TX2, C1, T1, T2, XIM2, XI
!-----------------------------------------------
!
      NM1 = N - 1
      NP1 = N + 1
      NS2 = N/2
      IF (N - 2 >= 0) THEN
         IF (N - 2 <= 0) THEN
            X1H = X(1) + X(2)
            X(2) = X(1) - X(2)
            X(1) = X1H
            RETURN 
         ENDIF
         IF (N <= 3) THEN
            X1P3 = X(1) + X(3)
            TX2 = X(2) + X(2)
            X(2) = X(1) - X(3)
            X(1) = X1P3 + TX2
            X(3) = X1P3 - TX2
            RETURN 
         ENDIF
         C1 = X(1) - X(N)
         X(1) = X(1) + X(N)
         DO K = 2, NS2
            KC = NP1 - K
            T1 = X(K) + X(KC)
            T2 = X(K) - X(KC)
            C1 = C1 + WSAVE(KC)*T2
            T2 = WSAVE(K)*T2
            X(K) = T1 - T2
            X(KC) = T1 + T2
         END DO
         MODN = MOD(N,2)
         IF (MODN /= 0) X(NS2+1) = X(NS2+1) + X(NS2+1)
         CALL RFFTF (NM1, X, WSAVE(N+1))
         XIM2 = X(2)
         X(2) = C1
         DO I = 4, N, 2
            XI = X(I)
            X(I) = X(I-2) - X(I-1)
            X(I-1) = XIM2
            XIM2 = XI
         END DO
         IF (MODN /= 0) X(N) = XIM2
      ENDIF
      RETURN 
      
END SUBROUTINE COST


SUBROUTINE SINTI(N, WSAVE)
       IMPLICIT REAL*8(A-H,O-Z)
       
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      REAL*8  :: WSAVE(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NS2, NP1, K
      REAL*8  :: PI, DUM, DT
!-----------------------------------------------
!
      PI = 4.D0*ATAN(1.D0)
      IF (N <= 1) RETURN 
      NS2 = N/2
      NP1 = N + 1
      DT  = PI/FLOAT(NP1)
      DO K = 1, NS2
         WSAVE(K) = 2.D0*SIN(K*DT)
      END DO
      CALL RFFTI (NP1, WSAVE(NS2+1))
      RETURN 
      
END SUBROUTINE SINTI


SUBROUTINE SINT(N, X, WSAVE)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      REAL*8   :: X(*)
      REAL*8   :: WSAVE(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NP1, IW1, IW2, IW3
!-----------------------------------------------
!
      NP1 = N + 1
      IW1 = N/2 + 1
      IW2 = IW1 + NP1
      IW3 = IW2 + NP1
      CALL SINT1 (N, X, WSAVE, WSAVE(IW1), WSAVE(IW2), WSAVE(IW3))
      RETURN 
      
END SUBROUTINE SINT


SUBROUTINE SINT1(N, WAR, WAS, XH, X, IFAC)
       IMPLICIT REAL*8(A-H,O-Z)
       
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      INTEGER  :: IFAC(*)
      REAL*8   :: WAR(*)
      REAL*8, INTENT(IN) :: WAS(*)
      REAL*8   :: XH(*)
      REAL*8   :: X(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, NP1, NS2, K, KC, MODN
      REAL*8  :: SQRT3, XHOLD, T1, T2
      
!-----------------------------------------------
      DATA SQRT3/ 1.73205080756888D0/ 
      
      XH(:N)  = WAR(:N)
      WAR(:N) = X(:N)
      IF (N - 2 <= 0) THEN
         IF (N - 2 /= 0) THEN
            XH(1) = XH(1) + XH(1)
            GO TO 106
         ENDIF
         XHOLD = SQRT3*(XH(1)+XH(2))
         XH(2) = SQRT3*(XH(1)-XH(2))
         XH(1) = XHOLD
         GO TO 106
      ENDIF
      NP1 = N + 1
      NS2 = N/2
      X(1) = 0.D0
      DO K = 1, NS2
         KC = NP1 - K
         T1 = XH(K) - XH(KC)
         T2 = WAS(K)*(XH(K)+XH(KC))
         X(K+1)  = T1 + T2
         X(KC+1) = T2 - T1
      END DO
      MODN = MOD(N,2)
      IF (MODN /= 0) X(NS2+2) = 4.*XH(NS2+1)
      CALL RFFTF1 (NP1, X, XH, WAR, IFAC)
      XH(1) = 0.5D0*X(1)
      DO I = 3, N, 2
         XH(I-1) = -X(I)
         XH(I) = XH(I-2) + X(I-1)
      END DO
      IF (MODN == 0) XH(N) = -X(N+1)
  106 CONTINUE
      X(:N) = WAR(:N)
      WAR(:N) = XH(:N)
      RETURN 
      
END SUBROUTINE SINT1


SUBROUTINE COSQI(N, WSAVE)
        
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      REAL*8   :: WSAVE(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K
      REAL*8  :: PIH, DUM, DT, FK
!-----------------------------------------------
!
      PIH = 2.D0*ATAN(1.D0)
      DT = PIH/FLOAT(N)
      FK = 0.D0
      DO K = 1, N
         FK = FK + 1.D0
         WSAVE(K) = COS(FK*DT)
      END DO
      CALL RFFTI (N, WSAVE(N+1))
      RETURN 
      
END SUBROUTINE COSQI


SUBROUTINE COSQF(N, X, WSAVE)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER :: N
      REAL*8  :: X(*)
      REAL*8  :: WSAVE(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL*8 :: SQRT2, TSQX
!-----------------------------------------------
      DATA SQRT2/ 1.414213562373095D0/ 
!
      IF (N - 2 >= 0) THEN
         IF (N - 2 > 0) GO TO 103
         TSQX = SQRT2*X(2)
         X(2) = X(1) - TSQX
         X(1) = X(1) + TSQX
      ENDIF
      RETURN 
  103 CONTINUE
      CALL COSQF1 (N, X, WSAVE, WSAVE(N+1))
      RETURN 
      
END SUBROUTINE COSQF


SUBROUTINE COSQF1(N, X, W, XH)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER :: N
      REAL*8  :: X(*)
      REAL*8 , INTENT(IN) :: W(*)
      REAL*8  :: XH(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NS2, NP2, K, KC, MODN, I
      REAL*8  :: XIM1
!-----------------------------------------------

      NS2 = (N + 1)/2
      NP2 = N + 2
      DO K = 2, NS2
         KC = NP2 - K
         XH(K) = X(K) + X(KC)
         XH(KC) = X(K) - X(KC)
      END DO
      MODN = MOD(N,2)
      IF (MODN == 0) XH(NS2+1) = X(NS2+1) + X(NS2+1)
      DO K = 2, NS2
         KC = NP2 - K
         X(K) = W(K-1)*XH(KC) + W(KC-1)*XH(K)
         X(KC) = W(K-1)*XH(K) - W(KC-1)*XH(KC)
      END DO
      IF (MODN == 0) X(NS2+1) = W(NS2)*XH(NS2+1)
      CALL RFFTF (N, X, XH)
      DO I = 3, N, 2
         XIM1 = X(I-1) - X(I)
         X(I) = X(I-1) + X(I)
         X(I-1) = XIM1
      END DO
      RETURN
       
END SUBROUTINE COSQF1


SUBROUTINE COSQB(N, X, WSAVE)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER :: N
      REAL*8  :: X(*)
      REAL*8  :: WSAVE(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL*8 :: TSQRT2, X1
!-----------------------------------------------
      DATA TSQRT2/ 2.82842712474619D0/ 
!
      IF (N - 2 <= 0) THEN
         IF (N - 2 /= 0) THEN
            X(1) = 4.D0*X(1)
            RETURN 
         ENDIF
         X1 = 4.D0*(X(1)+X(2))
         X(2) = TSQRT2*(X(1)-X(2))
         X(1) = X1
         RETURN 
      ENDIF
      CALL COSQB1 (N, X, WSAVE, WSAVE(N+1))
      RETURN 
END SUBROUTINE COSQB


SUBROUTINE COSQB1(N, X, W, XH)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: N
      REAL*8  :: X(*)
      REAL*8 , INTENT(IN) :: W(*)
      REAL*8  :: XH(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NS2, NP2, I, MODN, K, KC
      REAL*8  :: XIM1
      
      NS2 = (N + 1)/2
      NP2 = N + 2
      DO I = 3, N, 2
         XIM1 = X(I-1) + X(I)
         X(I) = X(I) - X(I-1)
         X(I-1) = XIM1
      END DO
      X(1) = X(1) + X(1)
      
      MODN = MOD(N,2)
      IF (MODN == 0) X(N) = X(N) + X(N)
      CALL RFFTB (N, X, XH)
      DO K = 2, NS2
         KC = NP2 - K
         XH(K) = W(K-1)*X(KC) + W(KC-1)*X(K)
         XH(KC) = W(K-1)*X(K) - W(KC-1)*X(KC)
      END DO
      IF (MODN == 0) X(NS2+1) = W(NS2)*(X(NS2+1)+X(NS2+1))
      DO K = 2, NS2
         KC = NP2 - K
         X(K) = XH(K) + XH(KC)
         X(KC) = XH(K) - XH(KC)
      END DO
      X(1) = X(1) + X(1)
      RETURN 
END SUBROUTINE COSQB1


SUBROUTINE SINQI(N, WSAVE)

      INTEGER :: N
      REAL*8  :: WSAVE(*)
      
      CALL COSQI (N, WSAVE)
      RETURN 
END SUBROUTINE SINQI


SUBROUTINE SINQF(N, X, WSAVE)

      INTEGER :: N
      REAL*8  :: X(*)
      REAL*8  :: WSAVE(*)

      INTEGER :: NS2, K, KC
      REAL*8  :: XHOLD
!
      IF (N == 1) RETURN 
      NS2 = N/2
      DO K = 1, NS2
         KC = N - K
         XHOLD = X(K)
         X(K) = X(KC+1)
         X(KC+1) = XHOLD
      END DO
      CALL COSQF (N, X, WSAVE)
      X(2:N:2) = -X(2:N:2)
      RETURN 
END SUBROUTINE SINQF


SUBROUTINE SINQB(N, X, WSAVE)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER :: N
      REAL*8  :: X(*)
      REAL*8  :: WSAVE(*)
      
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: NS2, K, KC
      REAL*8  :: XHOLD
!-----------------------------------------------
!
      IF (N <= 1) THEN
         X(1) = 4.D0*X(1)
         RETURN 
      ENDIF
      NS2 = N/2
      X(2:N:2) = -X(2:N:2)
      CALL COSQB (N, X, WSAVE)
      DO K = 1, NS2
         KC = N - K
         XHOLD = X(K)
         X(K) = X(KC+1)
         X(KC+1) = XHOLD
      END DO
      RETURN 
END SUBROUTINE SINQB


SUBROUTINE CFFTI(N, WSAVE)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: N
      REAL*8  :: WSAVE(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IW1, IW2
!-----------------------------------------------
!
      IF (N == 1) RETURN 
      IW1 = N + N + 1
      IW2 = IW1 + N + N
      CALL CFFTI1 (N, WSAVE(IW1), WSAVE(IW2))
      RETURN 
END SUBROUTINE CFFTI


SUBROUTINE CFFTI1(N, WA, IFAC)
      IMPLICIT REAL*8(A-H,O-Z)
      
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(INOUT) :: IFAC(*)
      REAL*8 , INTENT(INOUT) :: WA(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(4) :: NTRYH
      INTEGER :: NL, NF, J, NTRY, NQ, NR, I, IB, L1, K1, IP, LD, L2, IDO, &
                 IDOT, IPM, I1, II
      REAL*8  :: TPI, DUM, ARGH, FI, ARGLD, ARG
!-----------------------------------------------
      DATA NTRYH(1), NTRYH(2), NTRYH(3), NTRYH(4)/ 3, 4, 2, 5/ 
      NL = N
      NF = 0
      J  = 0
  101 CONTINUE
      J = J + 1
      IF (J - 4 <= 0) THEN
         NTRY = NTRYH(J)
      ELSE
         NTRY = NTRY + 2
      ENDIF
  104 CONTINUE
      NQ = NL/NTRY
      NR = NL - NTRY*NQ
      IF (NR /= 0) GO TO 101
      NF = NF + 1
      IFAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY == 2) THEN
         IF (NF /= 1) THEN
            IFAC(NF+2:4:(-1)) = IFAC(NF+1:3:(-1))
            IFAC(3) = 2
         ENDIF
      ENDIF
      IF (NL /= 1) GO TO 104
      IFAC(1) = N
      IFAC(2) = NF
      TPI = 8.D0*ATAN(1.D0)
      ARGH = TPI/FLOAT(N)
      I = 2
      L1 = 1
      DO K1 = 1, NF
         IP = IFAC(K1+2)
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IDOT = IDO + IDO + 2
         IPM = IP - 1
         DO J = 1, IPM
            I1 = I
            WA(I-1) = 1.D0
            WA(I) = 0.D0
            LD = LD + L1
            FI = 0.D0
            ARGLD = FLOAT(LD)*ARGH
            DO II = 4, IDOT, 2
               I = I + 2
               FI = FI + 1.D0
               ARG = FI*ARGLD
               WA(I-1) = COS(ARG)
               WA(I) = SIN(ARG)
            END DO
            IF (IP <= 5) CYCLE 
            WA(I1-1) = WA(I-1)
            WA(I1) = WA(I)
         END DO
         L1 = L2
      END DO
      RETURN 
END SUBROUTINE CFFTI1


SUBROUTINE CFFTB(N, C, WSAVE)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      REAL*8  :: C(*)
      REAL*8  :: WSAVE(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IW1, IW2
!-----------------------------------------------
!
      IF (N == 1) RETURN 
      IW1 = N + N + 1
      IW2 = IW1 + N + N
      CALL CFFTB1 (N, C, WSAVE, WSAVE(IW1), WSAVE(IW2))
      RETURN 
      
END SUBROUTINE CFFTB


SUBROUTINE CFFTB1(N, C, CH, WA, IFAC)
      IMPLICIT REAL*8(A-H,O-Z)
      
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: IFAC(*)
      REAL*8  :: C(*)
      REAL*8  :: CH(*)
      REAL*8  :: WA(*)
!----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NF,NA,L1,IW,K1,IP,L2,IDO,IDOT,IDL1,IX2,IX3,IX4,NAC,N2,I
      
!-----------------------------------------------

      NF = IFAC(2)
      NA = 0
      L1 = 1
      IW = 1
      DO K1 = 1, NF
         IP = IFAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDOT = IDO + IDO
         IDL1 = IDOT*L1
         IF (IP == 4) THEN
            IX2 = IW + IDOT
            IX3 = IX2 + IDOT
            IF (NA == 0) THEN
               CALL PASSB4 (IDOT, L1, C, CH, WA(IW), WA(IX2), WA(IX3))
            ELSE
               CALL PASSB4 (IDOT, L1, CH, C, WA(IW), WA(IX2), WA(IX3))
            ENDIF
            NA = 1 - NA
         ELSE
            IF (IP == 2) THEN
               IF (NA == 0) THEN
                  CALL PASSB2 (IDOT, L1, C, CH, WA(IW))
               ELSE
                  CALL PASSB2 (IDOT, L1, CH, C, WA(IW))
               ENDIF
               NA = 1 - NA
            ELSE
               IF (IP == 3) THEN
                  IX2 = IW + IDOT
                  IF (NA == 0) THEN
                     CALL PASSB3 (IDOT, L1, C, CH, WA(IW), WA(IX2))
                  ELSE
                     CALL PASSB3 (IDOT, L1, CH, C, WA(IW), WA(IX2))
                  ENDIF
                  NA = 1 - NA
               ELSE
                  IF (IP == 5) THEN
                     IX2 = IW + IDOT
                     IX3 = IX2 + IDOT
                     IX4 = IX3 + IDOT
                     IF (NA == 0) THEN
                        CALL PASSB5 (IDOT, L1, C, CH, WA(IW), WA(IX2), & 
                                     WA(IX3), WA(IX4))
                     ELSE
                        CALL PASSB5 (IDOT, L1, CH, C, WA(IW), WA(IX2), & 
                                     WA(IX3), WA(IX4))
                     ENDIF
                     NA = 1 - NA
                  ELSE
                     IF (NA == 0) THEN
                        CALL PASSB (NAC, IDOT, IP, L1, IDL1, C, C, C, CH, &
                                    CH, WA(IW))
                     ELSE
                        CALL PASSB (NAC, IDOT, IP, L1, IDL1, CH, CH, CH, &
                                    C, C, WA(IW))
                     ENDIF
                     IF (NAC /= 0) NA = 1 - NA
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         L1 = L2
         IW = IW + (IP - 1)*IDOT
      END DO
      IF (NA == 0) RETURN 
      N2 = N + N
      C(:N2) = CH(:N2)
      RETURN 
      
END SUBROUTINE CFFTB1


SUBROUTINE PASSB2(IDO, L1, CC, CH, WA1)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL*8 , INTENT(IN)  :: CC(IDO,2,L1)
      REAL*8 , INTENT(OUT) :: CH(IDO,L1,2)
      REAL*8 , INTENT(IN)  :: WA1(1)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, I
      REAL*8  :: TR2, TI2
!-----------------------------------------------
      IF (IDO <= 2) THEN
         CH(1,:,1) = CC(1,1,:) + CC(1,2,:)
         CH(1,:,2) = CC(1,1,:) - CC(1,2,:)
         CH(2,:,1) = CC(2,1,:) + CC(2,2,:)
         CH(2,:,2) = CC(2,1,:) - CC(2,2,:)
         RETURN 
      ENDIF
      
      DO K = 1, L1
         DO I = 2, IDO, 2
            CH(I-1,K,1) = CC(I-1,1,K) + CC(I-1,2,K)
            TR2 = CC(I-1,1,K) - CC(I-1,2,K)
            CH(I,K,1) = CC(I,1,K) + CC(I,2,K)
            TI2 = CC(I,1,K) - CC(I,2,K)
            CH(I,K,2) = WA1(I-1)*TI2 + WA1(I)*TR2
            CH(I-1,K,2) = WA1(I-1)*TR2 - WA1(I)*TI2
         END DO
      END DO
      RETURN 
      
END SUBROUTINE PASSB2


SUBROUTINE PASSB3(IDO, L1, CC, CH, WA1, WA2)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL*8 , INTENT(IN) :: CC(IDO,3,L1)
      REAL*8 , INTENT(OUT) :: CH(IDO,L1,3)
      REAL*8 , INTENT(IN) :: WA1(*)
      REAL*8 , INTENT(IN) :: WA2(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, I
      REAL*8  :: TAUR,TAUI,TR2,CR2,TI2,CI2,CR3,CI3,DR2,DR3,DI2,DI3
!-----------------------------------------------

      DATA TAUR, TAUI/ -.5D0, 0.866025403784439D0/ 
      
      IF (IDO == 2) THEN
         DO K = 1, L1
            TR2 = CC(1,2,K) + CC(1,3,K)
            CR2 = CC(1,1,K) + TAUR*TR2
            CH(1,K,1) = CC(1,1,K) + TR2
            TI2 = CC(2,2,K) + CC(2,3,K)
            CI2 = CC(2,1,K) + TAUR*TI2
            CH(2,K,1) = CC(2,1,K) + TI2
            CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))
            CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))
            CH(1,K,2) = CR2 - CI3
            CH(1,K,3) = CR2 + CI3
            CH(2,K,2) = CI2 + CR3
            CH(2,K,3) = CI2 - CR3
         END DO
         RETURN 
      ENDIF
      DO K = 1, L1
         DO I = 2, IDO, 2
            TR2 = CC(I-1,2,K) + CC(I-1,3,K)
            CR2 = CC(I-1,1,K) + TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K) + TR2
            TI2 = CC(I,2,K) + CC(I,3,K)
            CI2 = CC(I,1,K) + TAUR*TI2
            CH(I,K,1) = CC(I,1,K) + TI2
            CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
            CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
            DR2 = CR2 - CI3
            DR3 = CR2 + CI3
            DI2 = CI2 + CR3
            DI3 = CI2 - CR3
            CH(I,K,2) = WA1(I-1)*DI2 + WA1(I)*DR2
            CH(I-1,K,2) = WA1(I-1)*DR2 - WA1(I)*DI2
            CH(I,K,3) = WA2(I-1)*DI3 + WA2(I)*DR3
            CH(I-1,K,3) = WA2(I-1)*DR3 - WA2(I)*DI3
         END DO
      END DO
      RETURN 
      
END SUBROUTINE PASSB3


SUBROUTINE PASSB4(IDO, L1, CC, CH, WA1, WA2, WA3)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL*8 , INTENT(IN) :: CC(IDO,4,L1)
      REAL*8, INTENT(OUT) :: CH(IDO,L1,4)
      REAL*8, INTENT(IN) :: WA1(*)
      REAL*8, INTENT(IN) :: WA2(*)
      REAL*8, INTENT(IN) :: WA3(*)
      
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, I
      REAL*8  :: TI1,TI2,TR4,TI3,TR1,TR2,TI4,TR3,CR3,CI3,CR2,CR4,CI2,CI4
      
!-----------------------------------------------

      IF (IDO == 2) THEN
         DO K = 1, L1
            TI1 = CC(2,1,K) - CC(2,3,K)
            TI2 = CC(2,1,K) + CC(2,3,K)
            TR4 = CC(2,4,K) - CC(2,2,K)
            TI3 = CC(2,2,K) + CC(2,4,K)
            TR1 = CC(1,1,K) - CC(1,3,K)
            TR2 = CC(1,1,K) + CC(1,3,K)
            TI4 = CC(1,2,K) - CC(1,4,K)
            TR3 = CC(1,2,K) + CC(1,4,K)
            CH(1,K,1) = TR2 + TR3
            CH(1,K,3) = TR2 - TR3
            CH(2,K,1) = TI2 + TI3
            CH(2,K,3) = TI2 - TI3
            CH(1,K,2) = TR1 + TR4
            CH(1,K,4) = TR1 - TR4
            CH(2,K,2) = TI1 + TI4
            CH(2,K,4) = TI1 - TI4
         END DO
         RETURN 
      ENDIF
      DO K = 1, L1
         DO I = 2, IDO, 2
            TI1 = CC(I,1,K) - CC(I,3,K)
            TI2 = CC(I,1,K) + CC(I,3,K)
            TI3 = CC(I,2,K) + CC(I,4,K)
            TR4 = CC(I,4,K) - CC(I,2,K)
            TR1 = CC(I-1,1,K) - CC(I-1,3,K)
            TR2 = CC(I-1,1,K) + CC(I-1,3,K)
            TI4 = CC(I-1,2,K) - CC(I-1,4,K)
            TR3 = CC(I-1,2,K) + CC(I-1,4,K)
            CH(I-1,K,1) = TR2 + TR3
            CR3 = TR2 - TR3
            CH(I,K,1) = TI2 + TI3
            CI3 = TI2 - TI3
            CR2 = TR1 + TR4
            CR4 = TR1 - TR4
            CI2 = TI1 + TI4
            CI4 = TI1 - TI4
            CH(I-1,K,2) = WA1(I-1)*CR2 - WA1(I)*CI2
            CH(I,K,2) = WA1(I-1)*CI2 + WA1(I)*CR2
            CH(I-1,K,3) = WA2(I-1)*CR3 - WA2(I)*CI3
            CH(I,K,3) = WA2(I-1)*CI3 + WA2(I)*CR3
            CH(I-1,K,4) = WA3(I-1)*CR4 - WA3(I)*CI4
            CH(I,K,4) = WA3(I-1)*CI4 + WA3(I)*CR4
         END DO
      END DO
      RETURN 
END SUBROUTINE PASSB4


SUBROUTINE PASSB5(IDO, L1, CC, CH, WA1, WA2, WA3, WA4)
      IMPLICIT REAL*8(A-H,O-Z)
      
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL*8, INTENT(IN) :: CC(IDO,5,L1)
      REAL*8, INTENT(OUT) :: CH(IDO,L1,5)
      REAL*8, INTENT(IN) :: WA1(*)
      REAL*8, INTENT(IN) :: WA2(*)
      REAL*8, INTENT(IN) :: WA3(*)
      REAL*8, INTENT(IN) :: WA4(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: K, I
      REAL*8  :: TR11, TI11, TR12, TI12, TI5, TI2, TI4, TI3, TR5, TR2, TR4, &
                 TR3, CR2, CI2, CR3, CI3, CR5, CI5, CR4, CI4, DR3, DR4, DI3, &
                 DI4, DR5, DR2, DI5, DI2
                 
!-----------------------------------------------

      DATA TR11, TI11, TR12, TI12/ 0.309016994374947D0, 0.951056516295154D0, &
                                   -.809016994374947D0, 0.587785252292473D0/ 
      IF (IDO == 2) THEN
         DO K = 1, L1
            TI5 = CC(2,2,K) - CC(2,5,K)
            TI2 = CC(2,2,K) + CC(2,5,K)
            TI4 = CC(2,3,K) - CC(2,4,K)
            TI3 = CC(2,3,K) + CC(2,4,K)
            TR5 = CC(1,2,K) - CC(1,5,K)
            TR2 = CC(1,2,K) + CC(1,5,K)
            TR4 = CC(1,3,K) - CC(1,4,K)
            TR3 = CC(1,3,K) + CC(1,4,K)
            CH(1,K,1) = CC(1,1,K) + TR2 + TR3
            CH(2,K,1) = CC(2,1,K) + TI2 + TI3
            CR2 = CC(1,1,K) + TR11*TR2 + TR12*TR3
            CI2 = CC(2,1,K) + TR11*TI2 + TR12*TI3
            CR3 = CC(1,1,K) + TR12*TR2 + TR11*TR3
            CI3 = CC(2,1,K) + TR12*TI2 + TR11*TI3
            CR5 = TI11*TR5 + TI12*TR4
            CI5 = TI11*TI5 + TI12*TI4
            CR4 = TI12*TR5 - TI11*TR4
            CI4 = TI12*TI5 - TI11*TI4
            CH(1,K,2) = CR2 - CI5
            CH(1,K,5) = CR2 + CI5
            CH(2,K,2) = CI2 + CR5
            CH(2,K,3) = CI3 + CR4
            CH(1,K,3) = CR3 - CI4
            CH(1,K,4) = CR3 + CI4
            CH(2,K,4) = CI3 - CR4
            CH(2,K,5) = CI2 - CR5
         END DO
         RETURN 
      ENDIF
      DO K = 1, L1
         DO I = 2, IDO, 2
            TI5 = CC(I,2,K) - CC(I,5,K)
            TI2 = CC(I,2,K) + CC(I,5,K)
            TI4 = CC(I,3,K) - CC(I,4,K)
            TI3 = CC(I,3,K) + CC(I,4,K)
            TR5 = CC(I-1,2,K) - CC(I-1,5,K)
            TR2 = CC(I-1,2,K) + CC(I-1,5,K)
            TR4 = CC(I-1,3,K) - CC(I-1,4,K)
            TR3 = CC(I-1,3,K) + CC(I-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K) + TR2 + TR3
            CH(I,K,1) = CC(I,1,K) + TI2 + TI3
            CR2 = CC(I-1,1,K) + TR11*TR2 + TR12*TR3
            CI2 = CC(I,1,K) + TR11*TI2 + TR12*TI3
            CR3 = CC(I-1,1,K) + TR12*TR2 + TR11*TR3
            CI3 = CC(I,1,K) + TR12*TI2 + TR11*TI3
            CR5 = TI11*TR5 + TI12*TR4
            CI5 = TI11*TI5 + TI12*TI4
            CR4 = TI12*TR5 - TI11*TR4
            CI4 = TI12*TI5 - TI11*TI4
            DR3 = CR3 - CI4
            DR4 = CR3 + CI4
            DI3 = CI3 + CR4
            DI4 = CI3 - CR4
            DR5 = CR2 + CI5
            DR2 = CR2 - CI5
            DI5 = CI2 - CR5
            DI2 = CI2 + CR5
            CH(I-1,K,2) = WA1(I-1)*DR2 - WA1(I)*DI2
            CH(I,K,2) = WA1(I-1)*DI2 + WA1(I)*DR2
            CH(I-1,K,3) = WA2(I-1)*DR3 - WA2(I)*DI3
            CH(I,K,3) = WA2(I-1)*DI3 + WA2(I)*DR3
            CH(I-1,K,4) = WA3(I-1)*DR4 - WA3(I)*DI4
            CH(I,K,4) = WA3(I-1)*DI4 + WA3(I)*DR4
            CH(I-1,K,5) = WA4(I-1)*DR5 - WA4(I)*DI5
            CH(I,K,5) = WA4(I-1)*DI5 + WA4(I)*DR5
         END DO
      END DO
      RETURN 
      
END SUBROUTINE PASSB5


SUBROUTINE PASSB(NAC, IDO, IP, L1, IDL1, CC, C1, C2, CH, CH2, WA)
      IMPLICIT REAL*8(A-H,O-Z)
      
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER , INTENT(OUT) :: NAC
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: IP
      INTEGER , INTENT(IN) :: L1
      INTEGER , INTENT(IN) :: IDL1
      REAL*8, INTENT(IN) :: CC(IDO,IP,L1)
      REAL*8, INTENT(OUT) :: C1(IDO,L1,IP)
      REAL*8, INTENT(INOUT) :: C2(IDL1,IP)
      REAL*8, INTENT(INOUT) :: CH(IDO,L1,IP)
      REAL*8, INTENT(INOUT) :: CH2(IDL1,IP)
      REAL*8, INTENT(IN) :: WA(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: IDOT, NT, IPP2, IPPH, IDP, J, JC, K, I, IDL, INC, L, LC, &
                 IK, IDLJ, IDIJ, IDJ
      REAL*8  :: WAR, WAI
      
!-----------------------------------------------

      IDOT = IDO/2
      NT   = IP*IDL1
      IPP2 = IP + 2
      IPPH = (IP + 1)/2
      IDP  = IP*IDO
!
      IF (IDO >= L1) THEN
         DO J = 2, IPPH
            JC = IPP2 - J
            CH(:,:,J) = CC(:,J,:) + CC(:,JC,:)
            CH(:,:,JC) = CC(:,J,:) - CC(:,JC,:)
         END DO
         CH(:,:,1) = CC(:,1,:)
      ELSE
         DO J = 2, IPPH
            JC = IPP2 - J
            CH(:,:,J) = CC(:,J,:) + CC(:,JC,:)
            CH(:,:,JC) = CC(:,J,:) - CC(:,JC,:)
         END DO
         CH(:,:,1) = CC(:,1,:)
      ENDIF
      IDL = 2 - IDO
      INC = 0
      DO L = 2, IPPH
         LC = IPP2 - L
         IDL = IDL + IDO
         C2(:,L) = CH2(:,1) + WA(IDL-1)*CH2(:,2)
         C2(:,LC) = WA(IDL)*CH2(:,IP)
         IDLJ = IDL
         INC = INC + IDO
         DO J = 3, IPPH
            JC = IPP2 - J
            IDLJ = IDLJ + INC
            IF (IDLJ > IDP) IDLJ = IDLJ - IDP
            WAR = WA(IDLJ-1)
            WAI = WA(IDLJ)
            C2(:,L) = C2(:,L) + WAR*CH2(:,J)
            C2(:,LC) = C2(:,LC) + WAI*CH2(:,JC)
         END DO
      END DO
      DO J = 2, IPPH
         CH2(:,1) = CH2(:,1) + CH2(:,J)
      END DO
      DO J = 2, IPPH
         JC = IPP2 - J
         CH2(:IDL1-1:2,J) = C2(:IDL1-1:2,J) - C2(2:IDL1:2,JC)
         CH2(:IDL1-1:2,JC) = C2(:IDL1-1:2,J) + C2(2:IDL1:2,JC)
         CH2(2:IDL1:2,J) = C2(2:IDL1:2,J) + C2(:IDL1-1:2,JC)
         CH2(2:IDL1:2,JC) = C2(2:IDL1:2,J) - C2(:IDL1-1:2,JC)
      END DO
      NAC = 1
      IF (IDO == 2) RETURN 
      NAC = 0
      C2(:,1) = CH2(:,1)
      C1(1,:,2:IP) = CH(1,:,2:IP)
      C1(2,:,2:IP) = CH(2,:,2:IP)
      IF (IDOT <= L1) THEN
         IDIJ = 0
         DO J = 2, IP
            IDIJ = IDIJ + 2
            DO I = 4, IDO, 2
               IDIJ = IDIJ + 2
               C1(I-1,:,J) = WA(IDIJ-1)*CH(I-1,:,J) - WA(IDIJ)*CH(I,:,J)
               C1(I,:,J) = WA(IDIJ-1)*CH(I,:,J) + WA(IDIJ)*CH(I-1,:,J)
            END DO
         END DO
         RETURN 
      ENDIF
      IDJ = 2 - IDO
      DO J = 2, IP
         IDJ = IDJ + IDO
         DO K = 1, L1
            IDIJ = IDJ
            C1(3:IDO-1:2,K,J) = WA(IDIJ+1:IDO-3+IDIJ:2)*CH(3:IDO-1:2,K,J) - &
                                WA(IDIJ+2:IDO-2+IDIJ:2)*CH(4:IDO:2,K,J)
            C1(4:IDO:2,K,J) = WA(IDIJ+1:IDO-3+IDIJ:2)*CH(4:IDO:2,K,J) + &
                              WA(IDIJ+2:IDO-2+IDIJ:2)*CH(3:IDO-1:2,K,J)
         END DO
      END DO
      RETURN 
END SUBROUTINE PASSB


SUBROUTINE CFFTF(N, C, WSAVE)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER  :: N
      REAL*8 :: C(*)
      REAL*8 :: WSAVE(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IW1, IW2
!-----------------------------------------------
!
      IF (N == 1) RETURN 
      IW1 = N + N + 1
      IW2 = IW1 + N + N
      CALL CFFTF1 (N, C, WSAVE, WSAVE(IW1), WSAVE(IW2))
      RETURN 
      
END SUBROUTINE CFFTF


SUBROUTINE CFFTF1(N, C, CH, WA, IFAC)
      IMPLICIT REAL*8(A-H,O-Z)
      
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: IFAC(*)
      REAL*8 :: C(*)
      REAL*8 :: CH(*)
      REAL*8 :: WA(*)
      
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: NF,NA,L1,IW,K1,IP,L2,IDO,IDOT,IDL1,IX2,IX3,IX4,NAC,N2,I
      
!-----------------------------------------------
      NF = IFAC(2)
      NA = 0
      L1 = 1
      IW = 1
      DO K1 = 1, NF
         IP = IFAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDOT = IDO + IDO
         IDL1 = IDOT*L1
         IF (IP == 4) THEN
            IX2 = IW + IDOT
            IX3 = IX2 + IDOT
            IF (NA == 0) THEN
               CALL PASSF4 (IDOT, L1, C, CH, WA(IW), WA(IX2), WA(IX3))
            ELSE
               CALL PASSF4 (IDOT, L1, CH, C, WA(IW), WA(IX2), WA(IX3))
            ENDIF
            NA = 1 - NA
         ELSE
            IF (IP == 2) THEN
               IF (NA == 0) THEN
                  CALL PASSF2 (IDOT, L1, C, CH, WA(IW))
               ELSE
                  CALL PASSF2 (IDOT, L1, CH, C, WA(IW))
               ENDIF
               NA = 1 - NA
            ELSE
               IF (IP == 3) THEN
                  IX2 = IW + IDOT
                  IF (NA == 0) THEN
                     CALL PASSF3 (IDOT, L1, C, CH, WA(IW), WA(IX2))
                  ELSE
                     CALL PASSF3 (IDOT, L1, CH, C, WA(IW), WA(IX2))
                  ENDIF
                  NA = 1 - NA
               ELSE
                  IF (IP == 5) THEN
                     IX2 = IW + IDOT
                     IX3 = IX2 + IDOT
                     IX4 = IX3 + IDOT
                     IF (NA == 0) THEN
                        CALL PASSF5 (IDOT, L1, C, CH, WA(IW), WA(IX2),WA(IX3), WA(IX4))
                     ELSE
                        CALL PASSF5 (IDOT, L1, CH, C, WA(IW), WA(IX2),WA(IX3), WA(IX4))
                     ENDIF
                     NA = 1 - NA
                  ELSE
                     IF (NA == 0) THEN
                        CALL PASSF (NAC, IDOT, IP, L1, IDL1, C, C, C, CH, CH, WA(IW))
                     ELSE
                        CALL PASSF (NAC, IDOT, IP, L1, IDL1, CH, CH, CH, C, C, WA(IW))
                     ENDIF
                     IF (NAC /= 0) NA = 1 - NA
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         L1 = L2
         IW = IW + (IP - 1)*IDOT
      END DO
      IF (NA == 0) RETURN 
      N2 = N + N
      C(:N2) = CH(:N2)
      RETURN 
END SUBROUTINE CFFTF1


SUBROUTINE PASSF2(IDO, L1, CC, CH, WA1)
      IMPLICIT REAL*8(A-H,O-Z)
       
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL*8, INTENT(IN) :: CC(IDO,2,L1)
      REAL*8, INTENT(OUT) :: CH(IDO,L1,2)
      REAL*8, INTENT(IN) :: WA1(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, I
      REAL*8  :: TR2, TI2
!-----------------------------------------------
      IF (IDO <= 2) THEN
         CH(1,:,1) = CC(1,1,:) + CC(1,2,:)
         CH(1,:,2) = CC(1,1,:) - CC(1,2,:)
         CH(2,:,1) = CC(2,1,:) + CC(2,2,:)
         CH(2,:,2) = CC(2,1,:) - CC(2,2,:)
         RETURN 
      ENDIF
      DO K = 1, L1
         DO I = 2, IDO, 2
            CH(I-1,K,1) = CC(I-1,1,K) + CC(I-1,2,K)
            TR2 = CC(I-1,1,K) - CC(I-1,2,K)
            CH(I,K,1) = CC(I,1,K) + CC(I,2,K)
            TI2 = CC(I,1,K) - CC(I,2,K)
            CH(I,K,2) = WA1(I-1)*TI2 - WA1(I)*TR2
            CH(I-1,K,2) = WA1(I-1)*TR2 + WA1(I)*TI2
         END DO
      END DO
      RETURN 
      
END SUBROUTINE PASSF2


SUBROUTINE PASSF3(IDO, L1, CC, CH, WA1, WA2)
      IMPLICIT REAL*8(A-H,O-Z)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL*8, INTENT(IN)   :: CC(IDO,3,L1)
      REAL*8, INTENT(OUT)  :: CH(IDO,L1,3)
      REAL*8, INTENT(IN)   :: WA1(*)
      REAL*8, INTENT(IN)   :: WA2(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, I
      REAL*8  :: TAUR,TAUI,TR2,CR2,TI2,CI2,CR3,CI3,DR2,DR3,DI2,DI3
      
!-----------------------------------------------

      DATA TAUR, TAUI/ -.5D0,  - 0.866025403784439D0/ 
      IF (IDO == 2) THEN
         DO K = 1, L1
            TR2 = CC(1,2,K) + CC(1,3,K)
            CR2 = CC(1,1,K) + TAUR*TR2
            CH(1,K,1) = CC(1,1,K) + TR2
            TI2 = CC(2,2,K) + CC(2,3,K)
            CI2 = CC(2,1,K) + TAUR*TI2
            CH(2,K,1) = CC(2,1,K) + TI2
            CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))
            CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))
            CH(1,K,2) = CR2 - CI3
            CH(1,K,3) = CR2 + CI3
            CH(2,K,2) = CI2 + CR3
            CH(2,K,3) = CI2 - CR3
         END DO
         RETURN 
      ENDIF
      DO K = 1, L1
         DO I = 2, IDO, 2
            TR2 = CC(I-1,2,K) + CC(I-1,3,K)
            CR2 = CC(I-1,1,K) + TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K) + TR2
            TI2 = CC(I,2,K) + CC(I,3,K)
            CI2 = CC(I,1,K) + TAUR*TI2
            CH(I,K,1) = CC(I,1,K) + TI2
            CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
            CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
            DR2 = CR2 - CI3
            DR3 = CR2 + CI3
            DI2 = CI2 + CR3
            DI3 = CI2 - CR3
            CH(I,K,2) = WA1(I-1)*DI2 - WA1(I)*DR2
            CH(I-1,K,2) = WA1(I-1)*DR2 + WA1(I)*DI2
            CH(I,K,3) = WA2(I-1)*DI3 - WA2(I)*DR3
            CH(I-1,K,3) = WA2(I-1)*DR3 + WA2(I)*DI3
         END DO
      END DO
      RETURN 
      
END SUBROUTINE PASSF3


SUBROUTINE PASSF4(IDO, L1, CC, CH, WA1, WA2, WA3)
      IMPLICIT REAL*8(A-H,O-Z)
      
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL*8, INTENT(IN) :: CC(IDO,4,L1)
      REAL*8, INTENT(OUT) :: CH(IDO,L1,4)
      REAL*8, INTENT(IN) :: WA1(*)
      REAL*8, INTENT(IN) :: WA2(*)
      REAL*8, INTENT(IN) :: WA3(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: K, I
      REAL*8  :: TI1,TI2,TR4,TI3,TR1,TR2,TI4,TR3,CR3,CI3,CR2,CR4,CI2,CI4
!-----------------------------------------------
      IF (IDO == 2) THEN
         DO K = 1, L1
            TI1 = CC(2,1,K) - CC(2,3,K)
            TI2 = CC(2,1,K) + CC(2,3,K)
            TR4 = CC(2,2,K) - CC(2,4,K)
            TI3 = CC(2,2,K) + CC(2,4,K)
            TR1 = CC(1,1,K) - CC(1,3,K)
            TR2 = CC(1,1,K) + CC(1,3,K)
            TI4 = CC(1,4,K) - CC(1,2,K)
            TR3 = CC(1,2,K) + CC(1,4,K)
            CH(1,K,1) = TR2 + TR3
            CH(1,K,3) = TR2 - TR3
            CH(2,K,1) = TI2 + TI3
            CH(2,K,3) = TI2 - TI3
            CH(1,K,2) = TR1 + TR4
            CH(1,K,4) = TR1 - TR4
            CH(2,K,2) = TI1 + TI4
            CH(2,K,4) = TI1 - TI4
         END DO
         RETURN 
      ENDIF
      DO K = 1, L1
         DO I = 2, IDO, 2
            TI1 = CC(I,1,K) - CC(I,3,K)
            TI2 = CC(I,1,K) + CC(I,3,K)
            TI3 = CC(I,2,K) + CC(I,4,K)
            TR4 = CC(I,2,K) - CC(I,4,K)
            TR1 = CC(I-1,1,K) - CC(I-1,3,K)
            TR2 = CC(I-1,1,K) + CC(I-1,3,K)
            TI4 = CC(I-1,4,K) - CC(I-1,2,K)
            TR3 = CC(I-1,2,K) + CC(I-1,4,K)
            CH(I-1,K,1) = TR2 + TR3
            CR3 = TR2 - TR3
            CH(I,K,1) = TI2 + TI3
            CI3 = TI2 - TI3
            CR2 = TR1 + TR4
            CR4 = TR1 - TR4
            CI2 = TI1 + TI4
            CI4 = TI1 - TI4
            CH(I-1,K,2) = WA1(I-1)*CR2 + WA1(I)*CI2
            CH(I,K,2) = WA1(I-1)*CI2 - WA1(I)*CR2
            CH(I-1,K,3) = WA2(I-1)*CR3 + WA2(I)*CI3
            CH(I,K,3) = WA2(I-1)*CI3 - WA2(I)*CR3
            CH(I-1,K,4) = WA3(I-1)*CR4 + WA3(I)*CI4
            CH(I,K,4) = WA3(I-1)*CI4 - WA3(I)*CR4
         END DO
      END DO
      RETURN 
END SUBROUTINE PASSF4


SUBROUTINE PASSF5(IDO, L1, CC, CH, WA1, WA2, WA3, WA4)
      IMPLICIT REAL*8(A-H,O-Z)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL*8 , INTENT(IN) :: CC(IDO,5,L1)
      REAL*8 , INTENT(OUT) :: CH(IDO,L1,5)
      REAL*8 , INTENT(IN) :: WA1(*)
      REAL*8 , INTENT(IN) :: WA2(*)
      REAL*8 , INTENT(IN) :: WA3(*)
      REAL*8 , INTENT(IN) :: WA4(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, I
      REAL*8  :: TR11, TI11, TR12, TI12, TI5, TI2, TI4, TI3, TR5, TR2, TR4, &
                 TR3, CR2, CI2, CR3, CI3, CR5, CI5, CR4, CI4, DR3, DR4, DI3, & 
                 DI4, DR5, DR2, DI5, DI2
!-----------------------------------------------

      DATA TR11, TI11, TR12, TI12/ 0.309016994374947D0, -.951056516295154D0, &
                                   -.809016994374947D0, -.587785252292473D0/ 
      IF (IDO == 2) THEN
         DO K = 1, L1
            TI5 = CC(2,2,K) - CC(2,5,K)
            TI2 = CC(2,2,K) + CC(2,5,K)
            TI4 = CC(2,3,K) - CC(2,4,K)
            TI3 = CC(2,3,K) + CC(2,4,K)
            TR5 = CC(1,2,K) - CC(1,5,K)
            TR2 = CC(1,2,K) + CC(1,5,K)
            TR4 = CC(1,3,K) - CC(1,4,K)
            TR3 = CC(1,3,K) + CC(1,4,K)
            CH(1,K,1) = CC(1,1,K) + TR2 + TR3
            CH(2,K,1) = CC(2,1,K) + TI2 + TI3
            CR2 = CC(1,1,K) + TR11*TR2 + TR12*TR3
            CI2 = CC(2,1,K) + TR11*TI2 + TR12*TI3
            CR3 = CC(1,1,K) + TR12*TR2 + TR11*TR3
            CI3 = CC(2,1,K) + TR12*TI2 + TR11*TI3
            CR5 = TI11*TR5 + TI12*TR4
            CI5 = TI11*TI5 + TI12*TI4
            CR4 = TI12*TR5 - TI11*TR4
            CI4 = TI12*TI5 - TI11*TI4
            CH(1,K,2) = CR2 - CI5
            CH(1,K,5) = CR2 + CI5
            CH(2,K,2) = CI2 + CR5
            CH(2,K,3) = CI3 + CR4
            CH(1,K,3) = CR3 - CI4
            CH(1,K,4) = CR3 + CI4
            CH(2,K,4) = CI3 - CR4
            CH(2,K,5) = CI2 - CR5
         END DO
         RETURN 
      ENDIF
      DO K = 1, L1
         DO I = 2, IDO, 2
            TI5 = CC(I,2,K) - CC(I,5,K)
            TI2 = CC(I,2,K) + CC(I,5,K)
            TI4 = CC(I,3,K) - CC(I,4,K)
            TI3 = CC(I,3,K) + CC(I,4,K)
            TR5 = CC(I-1,2,K) - CC(I-1,5,K)
            TR2 = CC(I-1,2,K) + CC(I-1,5,K)
            TR4 = CC(I-1,3,K) - CC(I-1,4,K)
            TR3 = CC(I-1,3,K) + CC(I-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K) + TR2 + TR3
            CH(I,K,1) = CC(I,1,K) + TI2 + TI3
            CR2 = CC(I-1,1,K) + TR11*TR2 + TR12*TR3
            CI2 = CC(I,1,K) + TR11*TI2 + TR12*TI3
            CR3 = CC(I-1,1,K) + TR12*TR2 + TR11*TR3
            CI3 = CC(I,1,K) + TR12*TI2 + TR11*TI3
            CR5 = TI11*TR5 + TI12*TR4
            CI5 = TI11*TI5 + TI12*TI4
            CR4 = TI12*TR5 - TI11*TR4
            CI4 = TI12*TI5 - TI11*TI4
            DR3 = CR3 - CI4
            DR4 = CR3 + CI4
            DI3 = CI3 + CR4
            DI4 = CI3 - CR4
            DR5 = CR2 + CI5
            DR2 = CR2 - CI5
            DI5 = CI2 - CR5
            DI2 = CI2 + CR5
            CH(I-1,K,2) = WA1(I-1)*DR2 + WA1(I)*DI2
            CH(I,K,2) = WA1(I-1)*DI2 - WA1(I)*DR2
            CH(I-1,K,3) = WA2(I-1)*DR3 + WA2(I)*DI3
            CH(I,K,3) = WA2(I-1)*DI3 - WA2(I)*DR3
            CH(I-1,K,4) = WA3(I-1)*DR4 + WA3(I)*DI4
            CH(I,K,4) = WA3(I-1)*DI4 - WA3(I)*DR4
            CH(I-1,K,5) = WA4(I-1)*DR5 + WA4(I)*DI5
            CH(I,K,5) = WA4(I-1)*DI5 - WA4(I)*DR5
         END DO
      END DO
      RETURN 
      
END SUBROUTINE PASSF5


SUBROUTINE PASSF(NAC, IDO, IP, L1, IDL1, CC, C1, C2, CH, CH2, WA)
      IMPLICIT REAL*8(A-H,O-Z)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(OUT) :: NAC
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: IP
      INTEGER , INTENT(IN) :: L1
      INTEGER , INTENT(IN) :: IDL1
      REAL*8 , INTENT(IN) :: CC(IDO,IP,L1)
      REAL*8 , INTENT(OUT) :: C1(IDO,L1,IP)
      REAL*8 , INTENT(INOUT) :: C2(IDL1,IP)
      REAL*8 , INTENT(INOUT) :: CH(IDO,L1,IP)
      REAL*8 , INTENT(INOUT) :: CH2(IDL1,IP)
      REAL*8 , INTENT(IN) :: WA(*)
      
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: IDOT, NT, IPP2, IPPH, IDP, J, JC, K, I, IDL, INC, L, LC, &
                 IK, IDLJ, IDIJ, IDJ
      REAL*8  :: WAR, WAI
!-----------------------------------------------
      IDOT = IDO/2
      NT   = IP*IDL1
      IPP2 = IP + 2
      IPPH = (IP + 1)/2
      IDP  = IP*IDO
!
      IF (IDO >= L1) THEN
         DO J = 2, IPPH
            JC = IPP2 - J
            CH(:,:,J) = CC(:,J,:) + CC(:,JC,:)
            CH(:,:,JC) = CC(:,J,:) - CC(:,JC,:)
         END DO
         CH(:,:,1) = CC(:,1,:)
      ELSE
         DO J = 2, IPPH
            JC = IPP2 - J
            CH(:,:,J) = CC(:,J,:) + CC(:,JC,:)
            CH(:,:,JC) = CC(:,J,:) - CC(:,JC,:)
         END DO
         CH(:,:,1) = CC(:,1,:)
      ENDIF
      IDL = 2 - IDO
      INC = 0
      DO L = 2, IPPH
         LC = IPP2 - L
         IDL = IDL + IDO
         C2(:,L) = CH2(:,1) + WA(IDL-1)*CH2(:,2)
         C2(:,LC) = -WA(IDL)*CH2(:,IP)
         IDLJ = IDL
         INC = INC + IDO
         DO J = 3, IPPH
            JC = IPP2 - J
            IDLJ = IDLJ + INC
            IF (IDLJ > IDP) IDLJ = IDLJ - IDP
            WAR = WA(IDLJ-1)
            WAI = WA(IDLJ)
            C2(:,L) = C2(:,L) + WAR*CH2(:,J)
            C2(:,LC) = C2(:,LC) - WAI*CH2(:,JC)
         END DO
      END DO
      DO J = 2, IPPH
         CH2(:,1) = CH2(:,1) + CH2(:,J)
      END DO
      DO J = 2, IPPH
         JC = IPP2 - J
         CH2(:IDL1-1:2,J) = C2(:IDL1-1:2,J) - C2(2:IDL1:2,JC)
         CH2(:IDL1-1:2,JC) = C2(:IDL1-1:2,J) + C2(2:IDL1:2,JC)
         CH2(2:IDL1:2,J) = C2(2:IDL1:2,J) + C2(:IDL1-1:2,JC)
         CH2(2:IDL1:2,JC) = C2(2:IDL1:2,J) - C2(:IDL1-1:2,JC)
      END DO
      NAC = 1
      IF (IDO == 2) RETURN 
      NAC = 0
      C2(:,1) = CH2(:,1)
      C1(1,:,2:IP) = CH(1,:,2:IP)
      C1(2,:,2:IP) = CH(2,:,2:IP)
      IF (IDOT <= L1) THEN
         IDIJ = 0
         DO J = 2, IP
            IDIJ = IDIJ + 2
            DO I = 4, IDO, 2
               IDIJ = IDIJ + 2
               C1(I-1,:,J) = WA(IDIJ-1)*CH(I-1,:,J) + WA(IDIJ)*CH(I,:,J)
               C1(I,:,J) = WA(IDIJ-1)*CH(I,:,J) - WA(IDIJ)*CH(I-1,:,J)
            END DO
         END DO
         RETURN 
      ENDIF
      IDJ = 2 - IDO
      DO J = 2, IP
         IDJ = IDJ + IDO
         DO K = 1, L1
            IDIJ = IDJ
            C1(3:IDO-1:2,K,J) = WA(IDIJ+1:IDO-3+IDIJ:2)*CH(3:IDO-1:2,K,J) + &
                                WA(IDIJ+2:IDO-2+IDIJ:2)*CH(4:IDO:2,K,J)
            C1(4:IDO:2,K,J) = WA(IDIJ+1:IDO-3+IDIJ:2)*CH(4:IDO:2,K,J) - &
                              WA(IDIJ+2:IDO-2+IDIJ:2)*CH(3:IDO-1:2,K,J)
         END DO
      END DO
      RETURN 
      
END SUBROUTINE PASSF


SUBROUTINE RFFTI(N, WSAVE)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      REAL*8 :: WSAVE(*)
!-----------------------------------------------
!
      IF (N == 1) RETURN 
      CALL RFFTI1 (N, WSAVE(N+1), WSAVE(2*N+1))
      RETURN 
      
END SUBROUTINE RFFTI


SUBROUTINE RFFTI1(N, WA, IFAC)
      IMPLICIT REAL*8(A-H,O-Z)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(INOUT) :: IFAC(*)
      REAL*8, INTENT(OUT) :: WA(*)
      
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER , DIMENSION(4) :: NTRYH
      INTEGER :: NL, NF, J, NTRY, NQ, NR, I, IB, IS, NFM1, L1, K1, IP, & 
                 LD, L2, IDO, IPM, II
      REAL*8  :: TPI, DUM, ARGH, ARGLD, FI, ARG
!-----------------------------------------------
      DATA NTRYH(1), NTRYH(2), NTRYH(3), NTRYH(4)/ 4, 2, 3, 5/ 
      
      NL = N
      NF = 0
      J = 0
  101 CONTINUE
      J = J + 1
      IF (J - 4 <= 0) THEN
         NTRY = NTRYH(J)
      ELSE
         NTRY = NTRY + 2
      ENDIF
  104 CONTINUE
      NQ = NL/NTRY
      NR = NL - NTRY*NQ
      IF (NR /= 0) GO TO 101
      NF = NF + 1
      IFAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY == 2) THEN
         IF (NF /= 1) THEN
            IFAC(NF+2:4:(-1)) = IFAC(NF+1:3:(-1))
            IFAC(3) = 2
         ENDIF
      ENDIF
      IF (NL /= 1) GO TO 104
      IFAC(1) = N
      IFAC(2) = NF
      TPI = 8.D0*ATAN(1.0)
      ARGH = TPI/FLOAT(N)
      IS = 0
      NFM1 = NF - 1
      L1 = 1
      IF (NFM1 == 0) RETURN 
      DO K1 = 1, NFM1
         IP = IFAC(K1+2)
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IPM = IP - 1
         DO J = 1, IPM
            LD = LD + L1
            I = IS
            ARGLD = FLOAT(LD)*ARGH
            FI = 0.D0
            DO II = 3, IDO, 2
               I = I + 2
               FI = FI + 1.D0
               ARG = FI*ARGLD
               WA(I-1) = COS(ARG)
               WA(I) = SIN(ARG)
            END DO
            IS = IS + IDO
         END DO
         L1 = L2
      END DO
      RETURN 
      
END SUBROUTINE RFFTI1


SUBROUTINE RFFTB(N, R, WSAVE)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      REAL*8 :: R(*)
      REAL*8 :: WSAVE(*)
!-----------------------------------------------
!
      IF (N == 1) RETURN 
      CALL RFFTB1 (N, R, WSAVE, WSAVE(N+1), WSAVE(2*N+1))
      RETURN 
END SUBROUTINE RFFTB


SUBROUTINE RFFTB1(N, C, CH, WA, IFAC)
      IMPLICIT REAL*8(A-H,O-Z)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: IFAC(*)
      REAL*8  :: C(*)
      REAL*8  :: CH(*)
      REAL*8  :: WA(*)
      
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NF, NA, L1, IW, K1, IP, L2, IDO, IDL1, IX2, IX3, IX4, I
      
!-----------------------------------------------

      NF = IFAC(2)
      NA = 0
      L1 = 1
      IW = 1
      DO K1 = 1, NF
         IP = IFAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDL1 = IDO*L1
         IF (IP == 4) THEN
            IX2 = IW + IDO
            IX3 = IX2 + IDO
            IF (NA == 0) THEN
               CALL RADB4 (IDO, L1, C, CH, WA(IW), WA(IX2), WA(IX3))
            ELSE
               CALL RADB4 (IDO, L1, CH, C, WA(IW), WA(IX2), WA(IX3))
            ENDIF
            NA = 1 - NA
         ELSE
            IF (IP == 2) THEN
               IF (NA == 0) THEN
                  CALL RADB2 (IDO, L1, C, CH, WA(IW))
               ELSE
                  CALL RADB2 (IDO, L1, CH, C, WA(IW))
               ENDIF
               NA = 1 - NA
            ELSE
               IF (IP == 3) THEN
                  IX2 = IW + IDO
                  IF (NA == 0) THEN
                     CALL RADB3 (IDO, L1, C, CH, WA(IW), WA(IX2))
                  ELSE
                     CALL RADB3 (IDO, L1, CH, C, WA(IW), WA(IX2))
                  ENDIF
                  NA = 1 - NA
               ELSE
                  IF (IP == 5) THEN
                     IX2 = IW + IDO
                     IX3 = IX2 + IDO
                     IX4 = IX3 + IDO
                     IF (NA == 0) THEN
                        CALL RADB5 (IDO, L1, C, CH, WA(IW), WA(IX2), WA(IX3), WA(IX4))
                     ELSE
                        CALL RADB5 (IDO, L1, CH, C, WA(IW), WA(IX2), WA(IX3), WA(IX4))
                     ENDIF
                     NA = 1 - NA
                  ELSE
                     IF (NA == 0) THEN
                        CALL RADBG(IDO,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
                     ELSE
                        CALL RADBG(IDO,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
                     ENDIF
                     IF (IDO == 1) NA = 1 - NA
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         L1 = L2
         IW = IW + (IP - 1)*IDO
      END DO
      IF (NA == 0) RETURN 
      C(:N) = CH(:N)
      RETURN 
      
END SUBROUTINE RFFTB1


SUBROUTINE RADB2(IDO, L1, CC, CH, WA1)
      IMPLICIT REAL*8(A-H,O-Z)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL*8 , INTENT(IN) :: CC(IDO,2,L1)
      REAL*8 , INTENT(OUT) :: CH(IDO,L1,2)
      REAL*8 , INTENT(IN) :: WA1(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, IDP2, I, IC
      REAL*8  :: TR2, TI2
!-----------------------------------------------
      CH(1,:,1) = CC(1,1,:) + CC(IDO,2,:)
      CH(1,:,2) = CC(1,1,:) - CC(IDO,2,:)
      IF (IDO - 2 >= 0) THEN
         IF (IDO - 2 /= 0) THEN
            IDP2 = IDO + 2
            DO K = 1, L1
               DO I = 3, IDO, 2
                  IC = IDP2 - I
                  CH(I-1,K,1) = CC(I-1,1,K) + CC(IC-1,2,K)
                  TR2 = CC(I-1,1,K) - CC(IC-1,2,K)
                  CH(I,K,1) = CC(I,1,K) - CC(IC,2,K)
                  TI2 = CC(I,1,K) + CC(IC,2,K)
                  CH(I-1,K,2) = WA1(I-2)*TR2 - WA1(I-1)*TI2
                  CH(I,K,2) = WA1(I-2)*TI2 + WA1(I-1)*TR2
               END DO
            END DO
            IF (MOD(IDO,2) == 1) RETURN 
         ENDIF
         CH(IDO,:,1) = CC(IDO,1,:) + CC(IDO,1,:)
         CH(IDO,:,2) = -(CC(1,2,:)+CC(1,2,:))
      ENDIF
      RETURN 
      
END SUBROUTINE RADB2


SUBROUTINE RADB3(IDO, L1, CC, CH, WA1, WA2)
      IMPLICIT REAL*8(A-H,O-Z)
      
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL*8 , INTENT(IN)  :: CC(IDO,3,L1)
      REAL*8 , INTENT(OUT) :: CH(IDO,L1,3)
      REAL*8 , INTENT(IN)  :: WA1(*)
      REAL*8 , INTENT(IN)  :: WA2(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, IDP2, I, IC
      REAL*8  :: TAUR,TAUI,TR2,CR2,CI3,TI2,CI2,CR3,DR2,DR3,DI2,DI3
      
!-----------------------------------------------

      DATA TAUR, TAUI/ -.5D0, 0.866025403784439D0/ 
      
      DO K = 1, L1
         TR2 = CC(IDO,2,K) + CC(IDO,2,K)
         CR2 = CC(1,1,K) + TAUR*TR2
         CH(1,K,1) = CC(1,1,K) + TR2
         CI3 = TAUI*(CC(1,3,K)+CC(1,3,K))
         CH(1,K,2) = CR2 - CI3
         CH(1,K,3) = CR2 + CI3
      END DO
      IF (IDO == 1) RETURN 
      IDP2 = IDO + 2
      DO K = 1, L1
         DO I = 3, IDO, 2
            IC = IDP2 - I
            TR2 = CC(I-1,3,K) + CC(IC-1,2,K)
            CR2 = CC(I-1,1,K) + TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K) + TR2
            TI2 = CC(I,3,K) - CC(IC,2,K)
            CI2 = CC(I,1,K) + TAUR*TI2
            CH(I,K,1) = CC(I,1,K) + TI2
            CR3 = TAUI*(CC(I-1,3,K)-CC(IC-1,2,K))
            CI3 = TAUI*(CC(I,3,K)+CC(IC,2,K))
            DR2 = CR2 - CI3
            DR3 = CR2 + CI3
            DI2 = CI2 + CR3
            DI3 = CI2 - CR3
            CH(I-1,K,2) = WA1(I-2)*DR2 - WA1(I-1)*DI2
            CH(I,K,2) = WA1(I-2)*DI2 + WA1(I-1)*DR2
            CH(I-1,K,3) = WA2(I-2)*DR3 - WA2(I-1)*DI3
            CH(I,K,3) = WA2(I-2)*DI3 + WA2(I-1)*DR3
         END DO
      END DO
      RETURN 
      
END SUBROUTINE RADB3


SUBROUTINE RADB4(IDO, L1, CC, CH, WA1, WA2, WA3)

      IMPLICIT REAL*8(A-H,O-Z)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL*8 , INTENT(IN) :: CC(IDO,4,L1)
      REAL*8 , INTENT(OUT) :: CH(IDO,L1,4)
      REAL*8 , INTENT(IN) :: WA1(*)
      REAL*8 , INTENT(IN) :: WA2(*)
      REAL*8 , INTENT(IN) :: WA3(*)
      
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: K, IDP2, I, IC
      REAL*8  :: SQRT2, TR1, TR2, TR3, TR4, TI1, TI2, TI3, TI4, CR3, CI3, &
                 CR2, CR4, CI2, CI4
     
!-----------------------------------------------

      DATA SQRT2/ 1.414213562373095D0/ 
      
      DO K = 1, L1
         TR1 = CC(1,1,K) - CC(IDO,4,K)
         TR2 = CC(1,1,K) + CC(IDO,4,K)
         TR3 = CC(IDO,2,K) + CC(IDO,2,K)
         TR4 = CC(1,3,K) + CC(1,3,K)
         CH(1,K,1) = TR2 + TR3
         CH(1,K,2) = TR1 - TR4
         CH(1,K,3) = TR2 - TR3
         CH(1,K,4) = TR1 + TR4
      END DO
      IF (IDO - 2 >= 0) THEN
         IF (IDO - 2 /= 0) THEN
            IDP2 = IDO + 2
            DO K = 1, L1
               DO I = 3, IDO, 2
                  IC = IDP2 - I
                  TI1 = CC(I,1,K) + CC(IC,4,K)
                  TI2 = CC(I,1,K) - CC(IC,4,K)
                  TI3 = CC(I,3,K) - CC(IC,2,K)
                  TR4 = CC(I,3,K) + CC(IC,2,K)
                  TR1 = CC(I-1,1,K) - CC(IC-1,4,K)
                  TR2 = CC(I-1,1,K) + CC(IC-1,4,K)
                  TI4 = CC(I-1,3,K) - CC(IC-1,2,K)
                  TR3 = CC(I-1,3,K) + CC(IC-1,2,K)
                  CH(I-1,K,1) = TR2 + TR3
                  CR3 = TR2 - TR3
                  CH(I,K,1) = TI2 + TI3
                  CI3 = TI2 - TI3
                  CR2 = TR1 - TR4
                  CR4 = TR1 + TR4
                  CI2 = TI1 + TI4
                  CI4 = TI1 - TI4
                  CH(I-1,K,2) = WA1(I-2)*CR2 - WA1(I-1)*CI2
                  CH(I,K,2) = WA1(I-2)*CI2 + WA1(I-1)*CR2
                  CH(I-1,K,3) = WA2(I-2)*CR3 - WA2(I-1)*CI3
                  CH(I,K,3) = WA2(I-2)*CI3 + WA2(I-1)*CR3
                  CH(I-1,K,4) = WA3(I-2)*CR4 - WA3(I-1)*CI4
                  CH(I,K,4) = WA3(I-2)*CI4 + WA3(I-1)*CR4
               END DO
            END DO
            IF (MOD(IDO,2) == 1) RETURN 
         ENDIF
         DO K = 1, L1
            TI1 = CC(1,2,K) + CC(1,4,K)
            TI2 = CC(1,4,K) - CC(1,2,K)
            TR1 = CC(IDO,1,K) - CC(IDO,3,K)
            TR2 = CC(IDO,1,K) + CC(IDO,3,K)
            CH(IDO,K,1) = TR2 + TR2
            CH(IDO,K,2) = SQRT2*(TR1 - TI1)
            CH(IDO,K,3) = TI2 + TI2
            CH(IDO,K,4) = -SQRT2*(TR1 + TI1)
         END DO
      ENDIF
      RETURN 
      
END SUBROUTINE RADB4


SUBROUTINE RADB5(IDO, L1, CC, CH, WA1, WA2, WA3, WA4)
      IMPLICIT REAL*8(A-H,O-Z)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL*8 , INTENT(IN) :: CC(IDO,5,L1)
      REAL*8 , INTENT(OUT) :: CH(IDO,L1,5)
      REAL*8 , INTENT(IN) :: WA1(*)
      REAL*8 , INTENT(IN) :: WA2(*)
      REAL*8 , INTENT(IN) :: WA3(*)
      REAL*8 , INTENT(IN) :: WA4(*)
      
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: K, IDP2, I, IC
      REAL*8  :: TR11, TI11, TR12, TI12, TI5, TI4, TR2, TR3, CR2, CR3, CI5, &
                 CI4, TI2, TI3, TR5, TR4, CI2, CI3, CR5, CR4, DR3, DR4, DI3, &
                 DI4, DR5, DR2, DI5, DI2
                 
!-----------------------------------------------

      DATA TR11, TI11, TR12, TI12/ 0.309016994374947D0, 0.951056516295154D0, &
                                   -.809016994374947D0, 0.587785252292473D0/ 
      DO K = 1, L1
         TI5 = CC(1,3,K) + CC(1,3,K)
         TI4 = CC(1,5,K) + CC(1,5,K)
         TR2 = CC(IDO,2,K) + CC(IDO,2,K)
         TR3 = CC(IDO,4,K) + CC(IDO,4,K)
         CH(1,K,1) = CC(1,1,K) + TR2 + TR3
         CR2 = CC(1,1,K) + TR11*TR2 + TR12*TR3
         CR3 = CC(1,1,K) + TR12*TR2 + TR11*TR3
         CI5 = TI11*TI5 + TI12*TI4
         CI4 = TI12*TI5 - TI11*TI4
         CH(1,K,2) = CR2 - CI5
         CH(1,K,3) = CR3 - CI4
         CH(1,K,4) = CR3 + CI4
         CH(1,K,5) = CR2 + CI5
      END DO
      IF (IDO == 1) RETURN 
      IDP2 = IDO + 2
      DO K = 1, L1
         DO I = 3, IDO, 2
            IC = IDP2 - I
            TI5 = CC(I,3,K) + CC(IC,2,K)
            TI2 = CC(I,3,K) - CC(IC,2,K)
            TI4 = CC(I,5,K) + CC(IC,4,K)
            TI3 = CC(I,5,K) - CC(IC,4,K)
            TR5 = CC(I-1,3,K) - CC(IC-1,2,K)
            TR2 = CC(I-1,3,K) + CC(IC-1,2,K)
            TR4 = CC(I-1,5,K) - CC(IC-1,4,K)
            TR3 = CC(I-1,5,K) + CC(IC-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K) + TR2 + TR3
            CH(I,K,1) = CC(I,1,K) + TI2 + TI3
            CR2 = CC(I-1,1,K) + TR11*TR2 + TR12*TR3
            CI2 = CC(I,1,K) + TR11*TI2 + TR12*TI3
            CR3 = CC(I-1,1,K) + TR12*TR2 + TR11*TR3
            CI3 = CC(I,1,K) + TR12*TI2 + TR11*TI3
            CR5 = TI11*TR5 + TI12*TR4
            CI5 = TI11*TI5 + TI12*TI4
            CR4 = TI12*TR5 - TI11*TR4
            CI4 = TI12*TI5 - TI11*TI4
            DR3 = CR3 - CI4
            DR4 = CR3 + CI4
            DI3 = CI3 + CR4
            DI4 = CI3 - CR4
            DR5 = CR2 + CI5
            DR2 = CR2 - CI5
            DI5 = CI2 - CR5
            DI2 = CI2 + CR5
            CH(I-1,K,2) = WA1(I-2)*DR2 - WA1(I-1)*DI2
            CH(I,K,2) = WA1(I-2)*DI2 + WA1(I-1)*DR2
            CH(I-1,K,3) = WA2(I-2)*DR3 - WA2(I-1)*DI3
            CH(I,K,3) = WA2(I-2)*DI3 + WA2(I-1)*DR3
            CH(I-1,K,4) = WA3(I-2)*DR4 - WA3(I-1)*DI4
            CH(I,K,4) = WA3(I-2)*DI4 + WA3(I-1)*DR4
            CH(I-1,K,5) = WA4(I-2)*DR5 - WA4(I-1)*DI5
            CH(I,K,5) = WA4(I-2)*DI5 + WA4(I-1)*DR5
         END DO
      END DO
      RETURN 
      
END SUBROUTINE RADB5


SUBROUTINE RADBG(IDO, IP, L1, IDL1, CC, C1, C2, CH, CH2, WA)
      IMPLICIT REAL*8(A-H,O-Z)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: IP
      INTEGER , INTENT(IN) :: L1
      INTEGER , INTENT(IN) :: IDL1
      REAL*8 , INTENT(IN) :: CC(IDO,IP,L1)
      REAL*8 , INTENT(INOUT) :: C1(IDO,L1,IP)
      REAL*8 , INTENT(INOUT) :: C2(IDL1,IP)
      REAL*8 , INTENT(INOUT) :: CH(IDO,L1,IP)
      REAL*8 , INTENT(INOUT) :: CH2(IDL1,IP)
      REAL*8 , INTENT(IN) :: WA(*)
      
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: IDP2,NBD,IPP2,IPPH,K,I,J,JC,J2,IC,L,LC,IK,IS,IDIJ
      REAL*8  :: TPI,DUM,ARG,DCP,DSP,AR1,AI1,AR1H,DC2,DS2,AR2,AI2,AR2H
      
!-----------------------------------------------

      TPI = 8.D0*ATAN(1.D0)
      ARG = TPI/FLOAT(IP)
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IDP2 = IDO + 2
      NBD = (IDO - 1)/2
      IPP2 = IP + 2
      IPPH = (IP + 1)/2
      IF (IDO >= L1) THEN
         CH(:,:,1) = CC(:,1,:)
      ELSE
         CH(:,:,1) = CC(:,1,:)
      ENDIF
      DO J = 2, IPPH
         JC = IPP2 - J
         J2 = J + J
         CH(1,:,J) = CC(IDO,J2-2,:) + CC(IDO,J2-2,:)
         CH(1,:,JC) = CC(1,J2-1,:) + CC(1,J2-1,:)
      END DO
      IF (IDO /= 1) THEN
         IF (NBD >= L1) THEN
            DO J = 2, IPPH
               JC = IPP2 - J
               CH(2:IDO-1:2,:,J) = CC(2:IDO-1:2,2*J-1,:) + &
                                   CC(IDP2-4:IDP2-1-IDO:(-2),2*J-2,:)
               CH(2:IDO-1:2,:,JC) = CC(2:IDO-1:2,2*J-1,:) - &
                                    CC(IDP2-4:IDP2-1-IDO:(-2),2*J-2,:)
               CH(3:IDO:2,:,J) = CC(3:IDO:2,2*J-1,:) - &
                                 CC(IDP2-3:IDP2-IDO:(-2),2*J-2,:)
               CH(3:IDO:2,:,JC) = CC(3:IDO:2,2*J-1,:) + &
                                  CC(IDP2-3:IDP2-IDO:(-2),2*J-2,:)
            END DO
         ELSE
            DO J = 2, IPPH
               JC = IPP2 - J
               CH(2:IDO-1:2,:,J) = CC(2:IDO-1:2,2*J-1,:) + &
                                   CC(IDP2-4:IDP2-1-IDO:(-2),2*J-2,:)
               CH(2:IDO-1:2,:,JC) = CC(2:IDO-1:2,2*J-1,:) - &
                                    CC(IDP2-4:IDP2-1-IDO:(-2),2*J-2,:)
               CH(3:IDO:2,:,J) = CC(3:IDO:2,2*J-1,:) - &
                                 CC(IDP2-3:IDP2-IDO:(-2),2*J-2,:)
               CH(3:IDO:2,:,JC) = CC(3:IDO:2,2*J-1,:) + &
                                  CC(IDP2-3:IDP2-IDO:(-2),2*J-2,:)
            END DO
         ENDIF
      ENDIF
      AR1 = 1.D0
      AI1 = 0.D0
      DO L = 2, IPPH
         LC = IPP2 - L
         AR1H = DCP*AR1 - DSP*AI1
         AI1 = DCP*AI1 + DSP*AR1
         AR1 = AR1H
         C2(:,L) = CH2(:,1) + AR1*CH2(:,2)
         C2(:,LC) = AI1*CH2(:,IP)
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
         DO J = 3, IPPH
            JC = IPP2 - J
            AR2H = DC2*AR2 - DS2*AI2
            AI2 = DC2*AI2 + DS2*AR2
            AR2 = AR2H
            C2(:,L) = C2(:,L) + AR2*CH2(:,J)
            C2(:,LC) = C2(:,LC) + AI2*CH2(:,JC)
         END DO
      END DO
      DO J = 2, IPPH
         CH2(:,1) = CH2(:,1) + CH2(:,J)
      END DO
      DO J = 2, IPPH
         JC = IPP2 - J
         CH(1,:,J) = C1(1,:,J) - C1(1,:,JC)
         CH(1,:,JC) = C1(1,:,J) + C1(1,:,JC)
      END DO
      IF (IDO /= 1) THEN
         IF (NBD >= L1) THEN
            DO J = 2, IPPH
               JC = IPP2 - J
               CH(2:IDO-1:2,:,J) = C1(2:IDO-1:2,:,J) - C1(3:IDO:2,:,JC)
               CH(2:IDO-1:2,:,JC) = C1(2:IDO-1:2,:,J) + C1(3:IDO:2,:,JC)
               CH(3:IDO:2,:,J) = C1(3:IDO:2,:,J) + C1(2:IDO-1:2,:,JC)
               CH(3:IDO:2,:,JC) = C1(3:IDO:2,:,J) - C1(2:IDO-1:2,:,JC)
            END DO
         ELSE
            DO J = 2, IPPH
               JC = IPP2 - J
               CH(2:IDO-1:2,:,J) = C1(2:IDO-1:2,:,J) - C1(3:IDO:2,:,JC)
               CH(2:IDO-1:2,:,JC) = C1(2:IDO-1:2,:,J) + C1(3:IDO:2,:,JC)
               CH(3:IDO:2,:,J) = C1(3:IDO:2,:,J) + C1(2:IDO-1:2,:,JC)
               CH(3:IDO:2,:,JC) = C1(3:IDO:2,:,J) - C1(2:IDO-1:2,:,JC)
            END DO
         ENDIF
      ENDIF
      IF (IDO == 1) RETURN 
      C2(:,1) = CH2(:,1)
      C1(1,:,2:IP) = CH(1,:,2:IP)
      IF (NBD <= L1) THEN
         IS = -IDO
         DO J = 2, IP
            IS = IS + IDO
            IDIJ = IS
            DO I = 3, IDO, 2
               IDIJ = IDIJ + 2
               C1(I-1,:,J) = WA(IDIJ-1)*CH(I-1,:,J) - WA(IDIJ)*CH(I,:,J)
               C1(I,:,J) = WA(IDIJ-1)*CH(I,:,J) + WA(IDIJ)*CH(I-1,:,J)
            END DO
         END DO
      ELSE
         IS = -IDO
         DO J = 2, IP
            IS = IS + IDO
            DO K = 1, L1
               IDIJ = IS
               C1(2:IDO-1:2,K,J) = WA(IDIJ+1:IDO-2+IDIJ:2)*CH(2:IDO-1:2,K,J) - &
                                   WA(IDIJ+2:IDO-1+IDIJ:2)*CH(3:IDO:2,K,J)
               C1(3:IDO:2,K,J) = WA(IDIJ+1:IDO-2+IDIJ:2)*CH(3:IDO:2,K,J) + &
                                 WA(IDIJ+2:IDO-1+IDIJ:2)*CH(2:IDO-1:2,K,J)
            END DO
         END DO
      ENDIF
      RETURN 
      
END SUBROUTINE RADBG
     
SUBROUTINE RFFTF(N, R, WSAVE)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      REAL*8 :: R(*)
      REAL*8 :: WSAVE(*)
!-----------------------------------------------
!
      IF (N == 1) RETURN 
      CALL RFFTF1 (N, R, WSAVE, WSAVE(N+1), WSAVE(2*N+1))
      
      RETURN 
      
END SUBROUTINE RFFTF


SUBROUTINE RFFTF1(N, C, CH, WA, IFAC)
      IMPLICIT REAL*8(A-H,O-Z)
      
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: IFAC(*)
      REAL*8  :: C(*)
      REAL*8  :: CH(*)
      REAL*8  :: WA(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: NF,NA,L2,IW,K1,KH,IP,L1,IDO,IDL1,IX2,IX3,IX4,I
      
!-----------------------------------------------
      NF = IFAC(2)
      NA = 1
      L2 = N
      IW = N
      DO K1 = 1, NF
         KH = NF - K1
         IP = IFAC(KH+3)
         L1 = L2/IP
         IDO = N/L2
         IDL1 = IDO*L1
         IW = IW - (IP - 1)*IDO
         NA = 1 - NA
         IF (IP == 4) THEN
            IX2 = IW + IDO
            IX3 = IX2 + IDO
            IF (NA == 0) THEN
               CALL RADF4 (IDO, L1, C, CH, WA(IW), WA(IX2), WA(IX3))
               GO TO 110
            ENDIF
            CALL RADF4 (IDO, L1, CH, C, WA(IW), WA(IX2), WA(IX3))
            GO TO 110
         ENDIF
         IF (IP == 2) THEN
            IF (NA == 0) THEN
               CALL RADF2 (IDO, L1, C, CH, WA(IW))
               GO TO 110
            ENDIF
            CALL RADF2 (IDO, L1, CH, C, WA(IW))
            GO TO 110
         ENDIF
  104    CONTINUE
         IF (IP == 3) THEN
            IX2 = IW + IDO
            IF (NA == 0) THEN
               CALL RADF3 (IDO, L1, C, CH, WA(IW), WA(IX2))
               GO TO 110
            ENDIF
            CALL RADF3 (IDO, L1, CH, C, WA(IW), WA(IX2))
            GO TO 110
         ENDIF
  106    CONTINUE
         IF (IP == 5) THEN
            IX2 = IW + IDO
            IX3 = IX2 + IDO
            IX4 = IX3 + IDO
            IF (NA == 0) THEN
               CALL RADF5(IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
               GO TO 110
            ENDIF
            CALL RADF5(IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
            GO TO 110
         ENDIF
  108    CONTINUE
         IF (IDO == 1) NA = 1 - NA
         IF (NA == 0) THEN
            CALL RADFG (IDO, IP, L1, IDL1, C, C, C, CH, CH, WA(IW))
            NA = 1
         ELSE
            CALL RADFG (IDO, IP, L1, IDL1, CH, CH, CH, C, C, WA(IW))
            NA = 0
         ENDIF
  110    CONTINUE
         L2 = L1
      END DO
      IF (NA == 1) RETURN 
      C(:N) = CH(:N)
      
      RETURN 
      
END SUBROUTINE RFFTF1


SUBROUTINE RADF2(IDO, L1, CC, CH, WA1)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL*8 , INTENT(IN) :: CC(IDO,L1,2)
      REAL*8 , INTENT(OUT) :: CH(IDO,2,L1)
      REAL*8 , INTENT(IN) :: WA1(*)
      
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: K, IDP2, I, IC
      REAL*8  :: TR2, TI2
      
!-----------------------------------------------

      CH(1,1,:) = CC(1,:,1) + CC(1,:,2)
      CH(IDO,2,:) = CC(1,:,1) - CC(1,:,2)
      IF (IDO - 2 >= 0) THEN
         IF (IDO - 2 /= 0) THEN
            IDP2 = IDO + 2
            DO K = 1, L1
               DO I = 3, IDO, 2
                  IC = IDP2 - I
                  TR2 = WA1(I-2)*CC(I-1,K,2) + WA1(I-1)*CC(I,K,2)
                  TI2 = WA1(I-2)*CC(I,K,2) - WA1(I-1)*CC(I-1,K,2)
                  CH(I,1,K) = CC(I,K,1) + TI2
                  CH(IC,2,K) = TI2 - CC(I,K,1)
                  CH(I-1,1,K) = CC(I-1,K,1) + TR2
                  CH(IC-1,2,K) = CC(I-1,K,1) - TR2
               END DO
            END DO
            IF (MOD(IDO,2) == 1) RETURN 
         ENDIF
         CH(1,2,:)   =-CC(IDO,:,2)
         CH(IDO,1,:) = CC(IDO,:,1)
      ENDIF
      RETURN 
      
END SUBROUTINE RADF2


SUBROUTINE RADF3(IDO, L1, CC, CH, WA1, WA2)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL*8 , INTENT(IN)  :: CC(IDO,L1,3)
      REAL*8 , INTENT(OUT) :: CH(IDO,3,L1)
      REAL*8 , INTENT(IN)  :: WA1(*)
      REAL*8 , INTENT(IN)  :: WA2(*)
      
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: K, IDP2, I, IC
      REAL*8  :: TAUR,TAUI,CR2,DR2,DI2,DR3,DI3,CI2,TR2,TI2,TR3,TI3
      
!-----------------------------------------------

      DATA TAUR, TAUI/ -0.5D0, 0.866025403784439D0/ 
      
      DO K = 1, L1
         CR2 = CC(1,K,2) + CC(1,K,3)
         CH(1,1,K) = CC(1,K,1) + CR2
         CH(1,3,K) = TAUI*(CC(1,K,3)-CC(1,K,2))
         CH(IDO,2,K) = CC(1,K,1) + TAUR*CR2
      END DO
      IF (IDO == 1) RETURN 
      IDP2 = IDO + 2
      DO K = 1, L1
         DO I = 3, IDO, 2
            IC = IDP2 - I
            DR2 = WA1(I-2)*CC(I-1,K,2) + WA1(I-1)*CC(I,K,2)
            DI2 = WA1(I-2)*CC(I,K,2) - WA1(I-1)*CC(I-1,K,2)
            DR3 = WA2(I-2)*CC(I-1,K,3) + WA2(I-1)*CC(I,K,3)
            DI3 = WA2(I-2)*CC(I,K,3) - WA2(I-1)*CC(I-1,K,3)
            CR2 = DR2 + DR3
            CI2 = DI2 + DI3
            CH(I-1,1,K) = CC(I-1,K,1) + CR2
            CH(I,1,K) = CC(I,K,1) + CI2
            TR2 = CC(I-1,K,1) + TAUR*CR2
            TI2 = CC(I,K,1) + TAUR*CI2
            TR3 = TAUI*(DI2 - DI3)
            TI3 = TAUI*(DR3 - DR2)
            CH(I-1,3,K) = TR2 + TR3
            CH(IC-1,2,K) = TR2 - TR3
            CH(I,3,K) = TI2 + TI3
            CH(IC,2,K) = TI3 - TI2
         END DO
      END DO
      RETURN 
      END SUBROUTINE RADF3


SUBROUTINE RADF4(IDO, L1, CC, CH, WA1, WA2, WA3)
       IMPLICIT REAL*8(A-H,O-Z)
       
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL*8 , INTENT(IN)  :: CC(IDO,L1,4)
      REAL*8 , INTENT(OUT) :: CH(IDO,4,L1)
      REAL*8 , INTENT(IN)  :: WA1(*)
      REAL*8 , INTENT(IN)  :: WA2(*)
      REAL*8 , INTENT(IN)  :: WA3(*)
      
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, IDP2, I, IC
      REAL*8  :: HSQT2, TR1, TR2, CR2, CI2, CR3, CI3, CR4, CI4, TR4, TI1, &
                 TI4, TI2, TI3, TR3
     
!-----------------------------------------------
      DATA HSQT2/ 0.7071067811865475D0/ 
      
      DO K = 1, L1
         TR1 = CC(1,K,2) + CC(1,K,4)
         TR2 = CC(1,K,1) + CC(1,K,3)
         CH(1,1,K) = TR1 + TR2
         CH(IDO,4,K) = TR2 - TR1
         CH(IDO,2,K) = CC(1,K,1) - CC(1,K,3)
         CH(1,3,K)   = CC(1,K,4) - CC(1,K,2)
      END DO
      IF (IDO - 2 >= 0) THEN
         IF (IDO - 2 /= 0) THEN
            IDP2 = IDO + 2
            DO K = 1, L1
               DO I = 3, IDO, 2
                  IC = IDP2 - I
                  CR2 = WA1(I-2)*CC(I-1,K,2) + WA1(I-1)*CC(I,K,2)
                  CI2 = WA1(I-2)*CC(I,K,2) - WA1(I-1)*CC(I-1,K,2)
                  CR3 = WA2(I-2)*CC(I-1,K,3) + WA2(I-1)*CC(I,K,3)
                  CI3 = WA2(I-2)*CC(I,K,3) - WA2(I-1)*CC(I-1,K,3)
                  CR4 = WA3(I-2)*CC(I-1,K,4) + WA3(I-1)*CC(I,K,4)
                  CI4 = WA3(I-2)*CC(I,K,4) - WA3(I-1)*CC(I-1,K,4)
                  TR1 = CR2 + CR4
                  TR4 = CR4 - CR2
                  TI1 = CI2 + CI4
                  TI4 = CI2 - CI4
                  TI2 = CC(I,K,1) + CI3
                  TI3 = CC(I,K,1) - CI3
                  TR2 = CC(I-1,K,1) + CR3
                  TR3 = CC(I-1,K,1) - CR3
                  CH(I-1,1,K)  = TR1 + TR2
                  CH(IC-1,4,K) = TR2 - TR1
                  CH(I,1,K)    = TI1 + TI2
                  CH(IC,4,K)   = TI1 - TI2
                  CH(I-1,3,K)  = TI4 + TR3
                  CH(IC-1,2,K) = TR3 - TI4
                  CH(I,3,K)    = TR4 + TI3
                  CH(IC,2,K)   = TR4 - TI3
               END DO
            END DO
            IF (MOD(IDO,2) == 1) RETURN 
         ENDIF
         DO K = 1, L1
            TI1 = -HSQT2*(CC(IDO,K,2)+CC(IDO,K,4))
            TR1 = HSQT2*(CC(IDO,K,2)-CC(IDO,K,4))
            CH(IDO,1,K) = TR1 + CC(IDO,K,1)
            CH(IDO,3,K) = CC(IDO,K,1) - TR1
            CH(1,2,K)   = TI1 - CC(IDO,K,3)
            CH(1,4,K)   = TI1 + CC(IDO,K,3)
         END DO
      ENDIF
      RETURN 
      
END SUBROUTINE RADF4


SUBROUTINE RADF5(IDO, L1, CC, CH, WA1, WA2, WA3, WA4)
      IMPLICIT REAL*8(A-H,O-Z)
      
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL*8 , INTENT(IN)  :: CC(IDO,L1,5)
      REAL*8 , INTENT(OUT) :: CH(IDO,5,L1)
      REAL*8 , INTENT(IN)  :: WA1(*)
      REAL*8 , INTENT(IN)  :: WA2(*)
      REAL*8 , INTENT(IN)  :: WA3(*)
      REAL*8 , INTENT(IN)  :: WA4(*)
      
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: K, IDP2, I, IC
      REAL*8  :: TR11, TI11, TR12, TI12, CR2, CI5, CR3, CI4, DR2, DI2, DR3, &
                 DI3, DR4, DI4, DR5, DI5, CR5, CI2, CR4, CI3, TR2, TI2, TR3, &
                 TI3, TR5, TI5, TR4, TI4
!-----------------------------------------------
      DATA TR11, TI11, TR12, TI12/ 0.309016994374947D0, 0.951056516295154D0, &
                                   -.809016994374947D0, 0.587785252292473D0/ 
      DO K = 1, L1
         CR2 = CC(1,K,5) + CC(1,K,2)
         CI5 = CC(1,K,5) - CC(1,K,2)
         CR3 = CC(1,K,4) + CC(1,K,3)
         CI4 = CC(1,K,4) - CC(1,K,3)
         CH(1,1,K) = CC(1,K,1) + CR2 + CR3
         CH(IDO,2,K) = CC(1,K,1) + TR11*CR2 + TR12*CR3
         CH(1,3,K) = TI11*CI5 + TI12*CI4
         CH(IDO,4,K) = CC(1,K,1) + TR12*CR2 + TR11*CR3
         CH(1,5,K) = TI12*CI5 - TI11*CI4
      END DO
      IF (IDO == 1) RETURN 
      IDP2 = IDO + 2
      DO K = 1, L1
         DO I = 3, IDO, 2
            IC = IDP2 - I
            DR2 = WA1(I-2)*CC(I-1,K,2) + WA1(I-1)*CC(I,K,2)
            DI2 = WA1(I-2)*CC(I,K,2) - WA1(I-1)*CC(I-1,K,2)
            DR3 = WA2(I-2)*CC(I-1,K,3) + WA2(I-1)*CC(I,K,3)
            DI3 = WA2(I-2)*CC(I,K,3) - WA2(I-1)*CC(I-1,K,3)
            DR4 = WA3(I-2)*CC(I-1,K,4) + WA3(I-1)*CC(I,K,4)
            DI4 = WA3(I-2)*CC(I,K,4) - WA3(I-1)*CC(I-1,K,4)
            DR5 = WA4(I-2)*CC(I-1,K,5) + WA4(I-1)*CC(I,K,5)
            DI5 = WA4(I-2)*CC(I,K,5) - WA4(I-1)*CC(I-1,K,5)
            CR2 = DR2 + DR5
            CI5 = DR5 - DR2
            CR5 = DI2 - DI5
            CI2 = DI2 + DI5
            CR3 = DR3 + DR4
            CI4 = DR4 - DR3
            CR4 = DI3 - DI4
            CI3 = DI3 + DI4
            CH(I-1,1,K) = CC(I-1,K,1) + CR2 + CR3
            CH(I,1,K) = CC(I,K,1) + CI2 + CI3
            TR2 = CC(I-1,K,1) + TR11*CR2 + TR12*CR3
            TI2 = CC(I,K,1) + TR11*CI2 + TR12*CI3
            TR3 = CC(I-1,K,1) + TR12*CR2 + TR11*CR3
            TI3 = CC(I,K,1) + TR12*CI2 + TR11*CI3
            TR5 = TI11*CR5 + TI12*CR4
            TI5 = TI11*CI5 + TI12*CI4
            TR4 = TI12*CR5 - TI11*CR4
            TI4 = TI12*CI5 - TI11*CI4
            CH(I-1,3,K) = TR2 + TR5
            CH(IC-1,2,K) = TR2 - TR5
            CH(I,3,K) = TI2 + TI5
            CH(IC,2,K) = TI5 - TI2
            CH(I-1,5,K) = TR3 + TR4
            CH(IC-1,4,K) = TR3 - TR4
            CH(I,5,K) = TI3 + TI4
            CH(IC,4,K) = TI4 - TI3
         END DO
      END DO
      RETURN 
      
END SUBROUTINE RADF5


SUBROUTINE RADFG(IDO, IP, L1, IDL1, CC, C1, C2, CH, CH2, WA)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: IP
      INTEGER , INTENT(IN) :: L1
      INTEGER , INTENT(IN) :: IDL1
      REAL*8 , INTENT(OUT) :: CC(IDO,IP,L1)
      REAL*8 , INTENT(INOUT) :: C1(IDO,L1,IP)
      REAL*8 , INTENT(INOUT) :: C2(IDL1,IP)
      REAL*8 , INTENT(INOUT) :: CH(IDO,L1,IP)
      REAL*8 , INTENT(INOUT) :: CH2(IDL1,IP)
      REAL*8 , INTENT(IN) :: WA(*)
      
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: IPPH,IPP2,IDP2,NBD,IK,J,K,IS,IDIJ,I,JC,L,LC,J2,IC
      REAL*8  :: TPI,DUM,ARG,DCP,DSP,AR1,AI1,AR1H,DC2,DS2,AR2,AI2,AR2H
!-----------------------------------------------

      TPI = 8.D0*ATAN(1.D0)
      ARG = TPI/FLOAT(IP)
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IPPH = (IP + 1)/2
      IPP2 = IP + 2
      IDP2 = IDO + 2
      NBD = (IDO - 1)/2
      IF (IDO /= 1) THEN
         CH2(:,1) = C2(:,1)
         CH(1,:,2:IP) = C1(1,:,2:IP)
         IF (NBD <= L1) THEN
            IS = -IDO
            DO J = 2, IP
               IS = IS + IDO
               IDIJ = IS
               DO I = 3, IDO, 2
                  IDIJ = IDIJ + 2
                  CH(I-1,:,J) = WA(IDIJ-1)*C1(I-1,:,J)+WA(IDIJ)*C1(I,:,J)
                  CH(I,:,J)   = WA(IDIJ-1)*C1(I,:,J)-WA(IDIJ)*C1(I-1,:,J)
               END DO
            END DO
         ELSE
            IS = -IDO
            DO J = 2, IP
               IS = IS + IDO
               DO K = 1, L1
                  IDIJ = IS
                  CH(2:IDO-1:2,K,J) = WA(IDIJ+1:IDO-2+IDIJ:2)*C1(2:IDO-1:2,K,J) + &
                                      WA(IDIJ+2:IDO-1+IDIJ:2)*C1(3:IDO:2,K,J)
                  CH(3:IDO:2,K,J)  = WA(IDIJ+1:IDO-2+IDIJ:2)*C1(3:IDO:2,K,J) &
                                   - WA(IDIJ+2:IDO-1+IDIJ:2)*C1(2:IDO-1:2,K,J)
               END DO
            END DO
         ENDIF
         IF (NBD >= L1) THEN
            DO J = 2, IPPH
               JC = IPP2 - J
               C1(2:IDO-1:2,:,J)  =CH(2:IDO-1:2,:,J)+CH(2:IDO-1:2,:,JC)
               C1(2:IDO-1:2,:,JC) = CH(3:IDO:2,:,J) - CH(3:IDO:2,:,JC)
               C1(3:IDO:2,:,J)    = CH(3:IDO:2,:,J) + CH(3:IDO:2,:,JC)
               C1(3:IDO:2,:,JC)   = CH(2:IDO-1:2,:,JC) - CH(2:IDO-1:2,:,J)
            END DO
            GO TO 121
         ENDIF
         DO J = 2, IPPH
            JC = IPP2 - J
            C1(2:IDO-1:2,:,J)  = CH(2:IDO-1:2,:,J) + CH(2:IDO-1:2,:,JC)
            C1(2:IDO-1:2,:,JC) = CH(3:IDO:2,:,J) - CH(3:IDO:2,:,JC)
            C1(3:IDO:2,:,J)    = CH(3:IDO:2,:,J) + CH(3:IDO:2,:,JC)
            C1(3:IDO:2,:,JC)   = CH(2:IDO-1:2,:,JC) - CH(2:IDO-1:2,:,J)
         END DO
         GO TO 121
      ENDIF
      C2(:,1) = CH2(:,1)
  121 CONTINUE
      DO J = 2, IPPH
         JC = IPP2 - J
         C1(1,:,J) = CH(1,:,J) + CH(1,:,JC)
         C1(1,:,JC) = CH(1,:,JC) - CH(1,:,J)
      END DO
!
      AR1 = 1.D0
      AI1 = 0.D0
      DO L = 2, IPPH
         LC = IPP2 - L
         AR1H = DCP*AR1 - DSP*AI1
         AI1 = DCP*AI1 + DSP*AR1
         AR1 = AR1H
         CH2(:,L) = C2(:,1) + AR1*C2(:,2)
         CH2(:,LC) = AI1*C2(:,IP)
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
         DO J = 3, IPPH
            JC = IPP2 - J
            AR2H = DC2*AR2 - DS2*AI2
            AI2 = DC2*AI2 + DS2*AR2
            AR2 = AR2H
            CH2(:,L) = CH2(:,L) + AR2*C2(:,J)
            CH2(:,LC) = CH2(:,LC) + AI2*C2(:,JC)
         END DO
      END DO
      DO J = 2, IPPH
         CH2(:,1) = CH2(:,1) + C2(:,J)
      END DO
!
      IF (IDO >= L1) THEN
         CC(:,1,:) = CH(:,:,1)
      ELSE
         CC(:,1,:) = CH(:,:,1)
      ENDIF
      CC(IDO,2:(IPPH-1)*2:2,:) = TRANSPOSE(CH(1,:,2:IPPH))
      CC(1,3:IPPH*2-1:2,:) = TRANSPOSE(CH(1,:,IPP2-2:IPP2-IPPH:(-1)))
      IF (IDO == 1) RETURN 
      IF (NBD >= L1) THEN
         CC(2:IDO-1:2,3:IPPH*2-1:2,:) = RESHAPE(SOURCE = CH(2:IDO-1:2,:,2:IPPH) + &
         CH(2:IDO-1:2,:,IPP2-2:IPP2-IPPH:(-1)),SHAPE = (/(IDO-1)/2,IPPH-1,L1/), &
         ORDER = (/1,3,2/))
         CC(IDP2-4:IDP2-1-IDO:(-2),2:(IPPH-1)*2:2,:) = RESHAPE(SOURCE = &
          CH(2:IDO-1:2,:,2:IPPH)-CH(2:IDO-1:2,:,IPP2-2:IPP2-IPPH:(-1)), &
          SHAPE = (/(IDO-1)/2,IPPH-1,L1/),ORDER = (/1,3,2/))
         CC(3:IDO:2,3:IPPH*2-1:2,:) = RESHAPE(SOURCE = CH(3:IDO:2,:,2:IPPH) + &
          CH(3:IDO:2,:,IPP2-2:IPP2-IPPH:(-1)),SHAPE = (/(IDO-1)/2,IPPH-1,L1/), &
          ORDER = (/1,3,2/))
         CC(IDP2-3:IDP2-IDO:(-2),2:(IPPH-1)*2:2,:) = RESHAPE(SOURCE = &
          CH(3:IDO:2,:,IPP2-2:IPP2-IPPH:(-1))-CH(3:IDO:2,:,2:IPPH),SHAPE = &
          (/(IDO-1)/2,IPPH-1,L1/),ORDER = (/1,3,2/))
         RETURN 
      ENDIF
      CC(2:IDO-1:2,3:IPPH*2-1:2,:) = RESHAPE(SOURCE = CH(2:IDO-1:2,:,2:IPPH)+ &
               CH(2:IDO-1:2,:,IPP2-2:IPP2-IPPH:(-1)),SHAPE = (/(IDO-1)/2,IPPH-1,L1/), &
               ORDER = (/1,3,2/))
      CC(IDP2-4:IDP2-1-IDO:(-2),2:(IPPH-1)*2:2,:) = RESHAPE(SOURCE = &
        CH(2:IDO-1:2,:,2:IPPH)-CH(2:IDO-1:2,:,IPP2-2:IPP2-IPPH:(-1)),SHAPE = &
        (/(IDO-1)/2,IPPH-1,L1/),ORDER = (/1,3,2/))
      CC(3:IDO:2,3:IPPH*2-1:2,:) = RESHAPE(SOURCE = CH(3:IDO:2,:,2:IPPH) &
         +CH(3:IDO:2,:,IPP2-2:IPP2-IPPH:(-1)),SHAPE = (/(IDO-1)/2,IPPH-1,L1/), &
         ORDER = (/1,3,2/))
      CC(IDP2-3:IDP2-IDO:(-2),2:(IPPH-1)*2:2,:) = RESHAPE(SOURCE = &
       CH(3:IDO:2,:,IPP2-2:IPP2-IPPH:(-1))-CH(3:IDO:2,:,2:IPPH),SHAPE = &
       (/(IDO-1)/2,IPPH-1,L1/),ORDER = (/1,3,2/))
      RETURN 
      
END SUBROUTINE RADFG

!     this function is define in the file comf.f
! SEPTEMBER 1973    VERSION 1
! APRIL     1976    VERSION 2
! JANUARY   1978    VERSION 3
! DECEMBER  1979    VERSION 3.1
! FEBRUARY  1985    DOCUMENTATION UPGRADE
! NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
! June 2004 2004    fortran 90 updates