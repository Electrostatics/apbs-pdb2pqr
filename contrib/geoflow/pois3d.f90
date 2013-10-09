SUBROUTINE POIS3D(LPEROD, L, C1, MPEROD, M, C2, NPEROD, N, A, B, C, &
                 LDIMF, MDIMF, F, IERROR)
      USE fish
      IMPLICIT REAL*8(A-H,O-Z)
      TYPE (fishworkspace) :: w
      
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER  :: LPEROD
      INTEGER  :: L
      INTEGER  :: MPEROD
      INTEGER  :: M
      INTEGER  :: NPEROD
      INTEGER  :: N
      INTEGER  :: LDIMF
      INTEGER  :: MDIMF
      INTEGER  :: IERROR
      REAL*8  :: C1
      REAL*8  :: C2
      REAL*8  :: A(*)
      REAL*8  :: B(*)
      REAL*8  :: C(*)
      REAL*8  :: F(LDIMF,MDIMF,*)
      
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: LP, MP, NP, K, IRWK, ICWK
      REAL*8, DIMENSION(6) :: SAVE
      
      LP = LPEROD + 1
      MP = MPEROD + 1
      NP = NPEROD + 1
!
!     CHECK FOR INVALID INPUT.
!
      IERROR = 0
      IF (LP<1 .OR. LP>5) IERROR = 1
      IF (L < 3)          IERROR = 2
      IF (MP<1 .OR. MP>5) IERROR = 3
      IF (M < 3)          IERROR = 4
      IF (NP<1 .OR. NP>2) IERROR = 5
      IF (N < 3)          IERROR = 6
      IF (LDIMF < L)      IERROR = 7
      IF (MDIMF < M)      IERROR = 8
      IF (NP == 1) THEN
         DO K = 1, N
            IF (A(K) /= C(1)) GO TO 102
            IF (C(K) /= C(1)) GO TO 102
            IF (B(K) /= B(1)) GO TO 102
         END DO
         GO TO 104
  102    CONTINUE
         IERROR = 9
      ENDIF
      
      IF (NPEROD==1 .AND. (A(1)/=0. .OR. C(N)/=0.)) IERROR = 10
      
  104 CONTINUE
  
      IF (IERROR /= 0) RETURN 
      
!     allocate required work space length (generous estimate)
      IRWK=30+L+M+2*N+MAX0(L,M,N)+7*(INT((L+1)/2)+INT((M+1)/2))
      ICWK = 0
      CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
!     check that allocation was successful
      IF (IERROR == 20) RETURN 
      call pois3dd(LPEROD,L,C1,MPEROD,M,C2,NPEROD,N,A,B,C,LDIMF, &
                   MDIMF,F,IERROR,w%rew)
!     release work space
      CALL FISHFIN (W)
      RETURN 
      
END SUBROUTINE POIS3D


 
SUBROUTINE POIS3DD(LPEROD, L, C1, MPEROD, M, C2, NPEROD, N, A, B, &
                  C, LDIMF, MDIMF, F, IERROR, W)
      IMPLICIT REAL*8(A-H,O-Z)
      
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER , INTENT(IN) :: LPEROD
      INTEGER  :: L
      INTEGER , INTENT(IN) :: MPEROD
      INTEGER  :: M
      INTEGER , INTENT(IN) :: NPEROD
      INTEGER  :: N
      INTEGER  :: LDIMF
      INTEGER  :: MDIMF
      INTEGER  :: IERROR
      REAL*8  :: C1
      REAL*8  :: C2
      REAL*8  :: A(*)
      REAL*8  :: B(*)
      REAL*8  :: C(*)
      REAL*8  :: F(LDIMF,MDIMF,*)
      REAL*8  :: W(*)
      
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: LP, MP, NP, IWYRT, IWT, IWD, IWBB, IWX, IWY, NH, NHM1, & 
                 NODD, I, J, K, NHPK, NHMK
      REAL*8, DIMENSION(6) :: SAVE
      

      LP = LPEROD + 1
      MP = MPEROD + 1
      NP = NPEROD + 1
      IWYRT = L + 1
      IWT = IWYRT + M
      IWD = IWT + MAX0(L,M,N) + 1
      IWBB = IWD + N
      IWX = IWBB + N
      IWY = IWX + 7*((L + 1)/2) + 15
      GO TO (105,114) NP
!
!     REORDER UNKNOWNS WHEN NPEROD = 0.
!
  105 CONTINUE
      NH = (N + 1)/2
      NHM1 = NH - 1
      NODD = 1
      IF (2*NH == N) NODD = 2
      DO I = 1, L
         DO J = 1, M
            DO K = 1, NHM1
               W(K) = F(I,J,NH-K) - F(I,J,K+NH)
               W(K+NH) = F(I,J,NH-K) + F(I,J,K+NH)
            END DO
            W(NH) = 2.D0*F(I,J,NH)
            GO TO (108,107) NODD
  107       CONTINUE
            W(N) = 2.D0*F(I,J,N)
  108       CONTINUE
            F(I,J,:N) = W(:N)
         END DO
      END DO
      SAVE(1) = C(NHM1)
      SAVE(2) = A(NH)
      SAVE(3) = C(NH)
      SAVE(4) = B(NHM1)
      SAVE(5) = B(N)
      SAVE(6) = A(N)
      C(NHM1) = 0.D0
      A(NH) = 0.D0
      C(NH) = 2.D0*C(NH)
      SELECT CASE (NODD) 
      CASE DEFAULT
         B(NHM1) = B(NHM1) - A(NH-1)
         B(N) = B(N) + A(N)
      CASE (2) 
         A(N) = C(NH)
      END SELECT
  114 CONTINUE
      CALL POS3D1 (LP, L, MP, M, N, A, B, C, LDIMF, MDIMF, F, W, W(IWYRT), &
                  W(IWT), W(IWD), W(IWX), W(IWY), C1, C2, W(IWBB))
      GO TO (115,122) NP
  115 CONTINUE
      DO I = 1, L
         DO J = 1, M
            W(NH-1:NH-NHM1:(-1))=0.5D0*(F(I,J,NH+1:NHM1+NH)+F(I,J,:NHM1))
            W(NH+1:NHM1+NH) = 0.5D0*(F(I,J,NH+1:NHM1+NH)-F(I,J,:NHM1))
            W(NH) = 0.5D0*F(I,J,NH)
            GO TO (118,117) NODD
  117       CONTINUE
            W(N) = 0.5D0*F(I,J,N)
  118       CONTINUE
            F(I,J,:N) = W(:N)
         END DO
      END DO
      C(NHM1) = SAVE(1)
      A(NH)   = SAVE(2)
      C(NH)   = SAVE(3)
      B(NHM1) = SAVE(4)
      B(N)    = SAVE(5)
      A(N)    = SAVE(6)
  122 CONTINUE
      RETURN 
      
END SUBROUTINE POIS3DD


SUBROUTINE POS3D1(LP, L, MP, M, N, A, B, C, LDIMF, MDIMF, F, XRT, &
                 YRT, T, D, WX, WY, C1, C2, BB)
      IMPLICIT REAL*8(A-H,O-Z)
      
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: LP
      INTEGER , INTENT(IN) :: L
      INTEGER , INTENT(IN) :: MP
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: LDIMF
      INTEGER , INTENT(IN) :: MDIMF
      REAL*8 , INTENT(IN) :: C1
      REAL*8 , INTENT(IN) :: C2
      REAL*8  :: A(*)
      REAL*8 , INTENT(IN) :: B(*)
      REAL*8  :: C(*)
      REAL*8 , INTENT(INOUT) :: F(LDIMF,MDIMF,1)
      REAL*8 , INTENT(INOUT) :: XRT(*)
      REAL*8 , INTENT(INOUT) :: YRT(*)
      REAL*8  :: T(*)
      REAL*8  :: D(*)
      REAL*8  :: WX(*)
      REAL*8  :: WY(*)
      REAL*8  :: BB(*)
      
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: LR, MR, NR, LRDEL, I, MRDEL, J, IFWRD, IS, K
      REAL*8  :: PI, DUM, SCALX, DX, DI, SCALY, DY, DJ
      
!-----------------------------------------------

      PI = 4.D0*ATAN(1.D0)
      LR = L
      MR = M
      NR = N
!
!     GENERATE TRANSFORM ROOTS
!
      LRDEL = ((LP - 1)*(LP - 3)*(LP - 5))/3
      SCALX = LR + LRDEL
      DX = PI/(2.D0*SCALX)
      GO TO (108,103,101,102,101) LP
  101 CONTINUE
      DI = 0.5D0
      SCALX = 2.D0*SCALX
      GO TO 104
  102 CONTINUE
      DI = 1.D0
      GO TO 104
  103 CONTINUE
      DI = 0.D0
  104 CONTINUE
      DO I = 1, LR
         XRT(I) = -4.D0*C1*SIN((FLOAT(I) - DI)*DX)**2
      END DO
      SCALX = 2.D0*SCALX
      GO TO (112,106,110,107,111) LP
  106 CONTINUE
      CALL SINTI (LR, WX)
      GO TO 112
  107 CONTINUE
      CALL COSTI (LR, WX)
      GO TO 112
  108 CONTINUE
      XRT(1) = 0.D0
      XRT(LR) = -4.D0*C1
      DO I = 3, LR, 2
         XRT(I-1) = -4.D0*C1*SIN(FLOAT(I - 1)*DX)**2
         XRT(I) = XRT(I-1)
      END DO
      CALL RFFTI (LR, WX)
      GO TO 112
  110 CONTINUE
      CALL SINQI (LR, WX)
      GO TO 112
  111 CONTINUE
      CALL COSQI (LR, WX)
  112 CONTINUE
      MRDEL = ((MP - 1)*(MP - 3)*(MP - 5))/3
      SCALY = MR + MRDEL
      DY = PI/(2.D0*SCALY)
      GO TO (120,115,113,114,113) MP
  113 CONTINUE
      DJ = 0.5
      SCALY = 2.D0*SCALY
      GO TO 116
  114 CONTINUE
      DJ = 1.D0
      GO TO 116
  115 CONTINUE
      DJ = 0.D0
  116 CONTINUE
      DO J = 1, MR
         YRT(J) = -4.D0*C2*SIN((FLOAT(J) - DJ)*DY)**2
      END DO
      SCALY = 2.D0*SCALY
      GO TO (124,118,122,119,123) MP
  118 CONTINUE
      CALL SINTI (MR, WY)
      GO TO 124
  119 CONTINUE
      CALL COSTI (MR, WY)
      GO TO 124
  120 CONTINUE
      YRT(1) = 0.D0
      YRT(MR) = -4.D0*C2
      DO J = 3, MR, 2
         YRT(J-1) = -4.D0*C2*SIN(FLOAT(J - 1)*DY)**2
         YRT(J) = YRT(J-1)
      END DO
      CALL RFFTI (MR, WY)
      GO TO 124
  122 CONTINUE
      CALL SINQI (MR, WY)
      GO TO 124
  123 CONTINUE
      CALL COSQI (MR, WY)
  124 CONTINUE
      IFWRD = 1
      IS = 1
  125 CONTINUE
!
!     TRANSFORM X
!
      DO J=1,MR
	 DO K=1,NR
	    DO I=1,LR
               T(I) = F(I,J,K)
	    END DO
            GO TO (127,130,131,134,135),LP
  127       GO TO (128,129),IFWRD
  128       CALL RFFTF (LR,T,WX)
            GO TO 138
  129       CALL RFFTB (LR,T,WX)
            GO TO 138
  130       CALL SINT (LR,T,WX)
            GO TO 138
  131       GO TO (132,133),IFWRD
  132       CALL SINQF (LR,T,WX)
            GO TO 138
  133       CALL SINQB (LR,T,WX)
            GO TO 138
  134       CALL COST (LR,T,WX)
            GO TO 138
  135       GO TO (136,137),IFWRD
  136       CALL COSQF (LR,T,WX)
            GO TO 138
  137       CALL COSQB (LR,T,WX)
  138       CONTINUE
	    DO I=1,LR
               F(I,J,K) = T(I)
	    END DO
	 END DO
      END DO
      GO TO (142,164) IFWRD
!
!     TRANSFORM Y
!
  142 CONTINUE
      DO I=1,LR
	 DO K=1,NR
	    DO J=1,MR
               T(J) = F(I,J,K)
	    END DO
            GO TO (144,147,148,151,152),MP
  144       GO TO (145,146),IFWRD
  145       CALL RFFTF (MR,T,WY)
            GO TO 155
  146       CALL RFFTB (MR,T,WY)
            GO TO 155
  147       CALL SINT (MR,T,WY)
            GO TO 155
  148       GO TO (149,150),IFWRD
  149       CALL SINQF (MR,T,WY)
            GO TO 155
  150       CALL SINQB (MR,T,WY)
            GO TO 155
  151       CALL COST (MR,T,WY)
            GO TO 155
  152       GO TO (153,154),IFWRD
  153       CALL COSQF (MR,T,WY)
            GO TO 155
  154       CALL COSQB (MR,T,WY)
  155       CONTINUE
	    DO J=1,MR
               F(I,J,K) = T(J)
	    END DO
	 END DO
      END DO
      GO TO (159,125) IFWRD
  159 CONTINUE
      DO I = 1, LR
         DO J = 1, MR
            BB(:NR) = B(:NR) + XRT(I) + YRT(J)
            T(:NR) = F(I,J,:NR)
            CALL TRID (NR, A, BB, C, T, D)
            F(I,J,:NR) = T(:NR)
         END DO
      END DO
      IFWRD = 2
      IS = -1
      GO TO 142
  164 CONTINUE
      F(:LR,:MR,:NR) = F(:LR,:MR,:NR)/(SCALX*SCALY)
      RETURN 
      
END SUBROUTINE POS3D1


SUBROUTINE TRID(MR, A, B, C, Y, D)

      IMPLICIT REAL*8(A-H,O-Z)
      
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER , INTENT(IN) :: MR
      REAL*8 , INTENT(IN) :: A(*)
      REAL*8 , INTENT(IN) :: B(*)
      REAL*8 , INTENT(IN) :: C(*)
      REAL*8 , INTENT(INOUT) :: Y(*)
      REAL*8 , INTENT(INOUT) :: D(*)
      
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: M, MM1, I, IP      
      REAL*8  :: Z
      
!-----------------------------------------------

      M = MR
      MM1 = M - 1
      Z = 1.D0/B(1)
      D(1) = C(1)*Z
      Y(1) = Y(1)*Z
      DO I = 2, MM1
         Z = 1.D0/(B(I)-A(I)*D(I-1))
         D(I) = C(I)*Z
         Y(I) = (Y(I)-A(I)*Y(I-1))*Z
      END DO
      Z = B(M) - A(M)*D(MM1)
      IF (Z == 0.D0) THEN
         Y(M) = 0.D0
      ELSE
         Y(M) = (Y(M)-A(M)*Y(MM1))/Z
      ENDIF
      DO IP = 1, MM1
         I = M - IP
         Y(I) = Y(I) - D(I)*Y(I+1)
      END DO
      RETURN 
!
! REVISION HISTORY---
!
! SEPTEMBER 1973    VERSION 1
! APRIL     1976    VERSION 2
! JANUARY   1978    VERSION 3
! DECEMBER  1979    VERSION 3.1
! FEBRUARY  1985    DOCUMENTATION UPGRADE
! NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
! June      2004    Version 5.0, Fortran 90 Changes
!-----------------------------------------------------------------------
END SUBROUTINE TRID
