
!----------------------------------------------------------------------------------
!
SUBROUTINE POISOLVE3D(B1,X,NX,NY,NZ,XLEFT,XRIGHT,YLEFT,YRIGHT,ZLEFT,ZRIGHT)
     IMPLICIT REAL*8(A-H,O-Z)
     
     REAL*8  :: X(NX,NY,NZ)
     REAL*8  :: XLEFT,XRIGHT,YLEFT,YRIGHT,ZLEFT,ZRIGHT
     real*8  :: B1(nx*ny*nz)

     INTEGER :: NX, NY, NZ
     
     REAL*8  :: ELMBDA, PETERB, BDXS, BDXF, BDYS, BDYF, BDZS, BDZF
     INTEGER :: LBDCND, MBDCND, NBDCND, LDIMF, MDIMF
     
     LBDCND = 1
     MBDCND = 1
     NBDCND = 1
     
     ELMBDA = 0.D0
     
     LDIMF  = NX
     MDIMF  = NY
 
     BDXS = 0.D0
     BDXF = 0.D0

     BDYS = 0.D0
     BDYF = 0.D0

     BDZS = 0.D0
     BDZF = 0.D0

     do i=1,nx
	     do j=1,ny
		     do k=1,nz
			     ijk=(i-1)*ny*nz+(j-1)*nz+k
				 x(i,j,k)=b1(ijk)
			 enddo
		 enddo
	 enddo
	 
     CALL HW3CRT(XLEFT, XRIGHT, NX-1, LBDCND, BDXS, BDXF, &
                 YLEFT, YRIGHT, NY-1, MBDCND, BDYS, BDYF, &
                 ZLEFT, ZRIGHT, NZ-1, NBDCND, BDZS, BDZF, &
                 ELMBDA, LDIMF, MDIMF, X, PETERB, IERROR)

END SUBROUTINE 

!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!     SUBROUTINE HW3CRT (XS,XF,L,LBDCND,BDXS,BDXF,YS,YF,M,MBDCND,BDYS,
!    +                   BDYF,ZS,ZF,N,NBDCND,BDZS,BDZF,ELMBDA,LDIMF,
!    +                   MDIMF,F,PERTRB,IERROR)
!
!
! DIMENSION OF           BDXS(MDIMF,N+1),    BDXF(MDIMF,N+1),
! ARGUMENTS              BDYS(LDIMF,N+1),    BDYF(LDIMF,N+1),
!                        BDZS(LDIMF,M+1),    BDZF(LDIMF,M+1),
!                        F(LDIMF,MDIMF,N+1)
!
! LATEST REVISION        June 2004
!
! PURPOSE                SOLVES THE STANDARD FIVE-POINT FINITE
!                        DIFFERENCE APPROXIMATION TO THE HELMHOLTZ
!                        EQUATION IN CARTESIAN COORDINATES.  THIS
!                        EQUATION IS
!
!                          (D/DX)(DU/DX) + (D/DY)(DU/DY) +
!                          (D/DZ)(DU/DZ) + LAMBDA*U = F(X,Y,Z) .
!
! USAGE                  CALL HW3CRT (XS,XF,L,LBDCND,BDXS,BDXF,YS,YF,M,
!                                     MBDCND,BDYS,BDYF,ZS,ZF,N,NBDCND,
!                                     BDZS,BDZF,ELMBDA,LDIMF,MDIMF,F,
!                                     PERTRB,IERROR)
!
! ARGUMENTS
!
! ON INPUT               XS,XF
!
!                          THE RANGE OF X, I.E. XS .LE. X .LE. XF .
!                          XS MUST BE LESS THAN XF.
!
!                        L
!                          THE NUMBER OF PANELS INTO WHICH THE
!                          INTERVAL (XS,XF) IS SUBDIVIDED.
!                          HENCE, THERE WILL BE L+1 GRID POINTS
!                          IN THE X-DIRECTION GIVEN BY
!                          X(I) = XS+(I-1)DX FOR I=1,2,...,L+1,
!                          WHERE DX = (XF-XS)/L IS THE PANEL WIDTH.
!                          L MUST BE AT LEAST 5.
!
!                        LBDCND
!                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
!                          AT X = XS AND X = XF.
!
!                          = 0  IF THE SOLUTION IS PERIODIC IN X,
!                               I.E. U(L+I,J,K) = U(I,J,K).
!                          = 1  IF THE SOLUTION IS SPECIFIED AT
!                               X = XS AND X = XF.
!                          = 2  IF THE SOLUTION IS SPECIFIED AT
!                               X = XS AND THE DERIVATIVE OF THE
!                               SOLUTION WITH RESPECT TO X IS
!                               SPECIFIED AT X = XF.
!                          = 3  IF THE DERIVATIVE OF THE SOLUTION
!                               WITH RESPECT TO X IS SPECIFIED AT
!                               X = XS AND X = XF.
!                          = 4  IF THE DERIVATIVE OF THE SOLUTION
!                               WITH RESPECT TO X IS SPECIFIED AT
!                               X = XS AND THE SOLUTION IS SPECIFIED
!                               AT X=XF.
!
!                        BDXS
!                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
!                          VALUES OF THE DERIVATIVE OF THE SOLUTION
!                          WITH RESPECT TO X AT X = XS.
!
!                          WHEN LBDCND = 3 OR 4,
!
!                            BDXS(J,K) = (D/DX)U(XS,Y(J),Z(K)),
!                            J=1,2,...,M+1,      K=1,2,...,N+1.
!
!                          WHEN LBDCND HAS ANY OTHER VALUE, BDXS
!                          IS A DUMMY VARIABLE. BDXS MUST BE
!                          DIMENSIONED AT LEAST (M+1)*(N+1).
!
!                        BDXF
!                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
!                          VALUES OF THE DERIVATIVE OF THE SOLUTION
!                          WITH RESPECT TO X AT X = XF.
!
!                          WHEN LBDCND = 2 OR 3,
!
!                            BDXF(J,K) = (D/DX)U(XF,Y(J),Z(K)),
!                            J=1,2,...,M+1,      K=1,2,...,N+1.
!
!                          WHEN LBDCND HAS ANY OTHER VALUE, BDXF IS
!                          A DUMMY VARIABLE.  BDXF MUST BE
!                          DIMENSIONED AT LEAST (M+1)*(N+1).
!
!                        YS,YF
!                          THE RANGE OF Y, I.E. YS .LE. Y .LE. YF.
!                          YS MUST BE LESS THAN YF.
!
!                        M
!                          THE NUMBER OF PANELS INTO WHICH THE
!                          INTERVAL (YS,YF) IS SUBDIVIDED.
!                          HENCE, THERE WILL BE M+1 GRID POINTS IN
!                          THE Y-DIRECTION GIVEN BY Y(J) = YS+(J-1)DY
!                          FOR J=1,2,...,M+1,
!                          WHERE DY = (YF-YS)/M IS THE PANEL WIDTH.
!                          M MUST BE AT LEAST 5.
!
!                        MBDCND
!                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
!                          AT Y = YS AND Y = YF.
!
!                          = 0  IF THE SOLUTION IS PERIODIC IN Y, I.E.
!                               U(I,M+J,K) = U(I,J,K).
!                          = 1  IF THE SOLUTION IS SPECIFIED AT
!                               Y = YS AND Y = YF.
!                          = 2  IF THE SOLUTION IS SPECIFIED AT
!                               Y = YS AND THE DERIVATIVE OF THE
!                               SOLUTION WITH RESPECT TO Y IS
!                               SPECIFIED AT Y = YF.
!                          = 3  IF THE DERIVATIVE OF THE SOLUTION
!                               WITH RESPECT TO Y IS SPECIFIED AT
!                               Y = YS AND Y = YF.
!                          = 4  IF THE DERIVATIVE OF THE SOLUTION
!                               WITH RESPECT TO Y IS SPECIFIED AT
!                               AT Y = YS AND THE SOLUTION IS
!                               SPECIFIED AT Y=YF.
!
!                        BDYS
!                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES
!                          THE VALUES OF THE DERIVATIVE OF THE
!                          SOLUTION WITH RESPECT TO Y AT Y = YS.
!
!                          WHEN MBDCND = 3 OR 4,
!
!                            BDYS(I,K) = (D/DY)U(X(I),YS,Z(K)),
!                            I=1,2,...,L+1,      K=1,2,...,N+1.
!
!                          WHEN MBDCND HAS ANY OTHER VALUE, BDYS
!                          IS A DUMMY VARIABLE. BDYS MUST BE
!                          DIMENSIONED AT LEAST (L+1)*(N+1).
!
!                        BDYF
!                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES
!                          THE VALUES OF THE DERIVATIVE OF THE
!                          SOLUTION WITH RESPECT TO Y AT Y = YF.
!
!                          WHEN MBDCND = 2 OR 3,
!
!                            BDYF(I,K) = (D/DY)U(X(I),YF,Z(K)),
!                            I=1,2,...,L+1,      K=1,2,...,N+1.
!
!                          WHEN MBDCND HAS ANY OTHER VALUE, BDYF
!                          IS A DUMMY VARIABLE. BDYF MUST BE
!                          DIMENSIONED AT LEAST (L+1)*(N+1).
!
!                        ZS,ZF
!                          THE RANGE OF Z, I.E. ZS .LE. Z .LE. ZF.
!                          ZS MUST BE LESS THAN ZF.
!
!                        N
!                          THE NUMBER OF PANELS INTO WHICH THE
!                          INTERVAL (ZS,ZF) IS SUBDIVIDED.
!                          HENCE, THERE WILL BE N+1 GRID POINTS
!                          IN THE Z-DIRECTION GIVEN BY
!                          Z(K) = ZS+(K-1)DZ FOR K=1,2,...,N+1,
!                          WHERE DZ = (ZF-ZS)/N IS THE PANEL WIDTH.
!                          N MUST BE AT LEAST 5.
!
!                        NBDCND
!                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
!                          AT Z = ZS AND Z = ZF.
!
!                          = 0  IF THE SOLUTION IS PERIODIC IN Z, I.E.
!                               U(I,J,N+K) = U(I,J,K).
!                          = 1  IF THE SOLUTION IS SPECIFIED AT
!                               Z = ZS AND Z = ZF.
!                          = 2  IF THE SOLUTION IS SPECIFIED AT
!                               Z = ZS AND THE DERIVATIVE OF THE
!                               SOLUTION WITH RESPECT TO Z IS
!                               SPECIFIED AT Z = ZF.
!                          = 3  IF THE DERIVATIVE OF THE SOLUTION
!                               WITH RESPECT TO Z IS SPECIFIED AT
!                               Z = ZS AND Z = ZF.
!                          = 4  IF THE DERIVATIVE OF THE SOLUTION
!                               WITH RESPECT TO Z IS SPECIFIED AT
!                               Z = ZS AND THE SOLUTION IS SPECIFIED
!                               AT Z=ZF.
!
!                        BDZS
!                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES
!                          THE VALUES OF THE DERIVATIVE OF THE
!                          SOLUTION WITH RESPECT TO Z AT Z = ZS.
!
!                          WHEN NBDCND = 3 OR 4,
!
!                            BDZS(I,J) = (D/DZ)U(X(I),Y(J),ZS),
!                            I=1,2,...,L+1,      J=1,2,...,M+1.
!
!                          WHEN NBDCND HAS ANY OTHER VALUE, BDZS
!                          IS A DUMMY VARIABLE. BDZS MUST BE
!                          DIMENSIONED AT LEAST (L+1)*(M+1).
!
!                        BDZF
!                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES
!                          THE VALUES OF THE DERIVATIVE OF THE
!                          SOLUTION WITH RESPECT TO Z AT Z = ZF.
!
!                          WHEN NBDCND = 2 OR 3,
!
!                            BDZF(I,J) = (D/DZ)U(X(I),Y(J),ZF),
!                            I=1,2,...,L+1,      J=1,2,...,M+1.
!
!                          WHEN NBDCND HAS ANY OTHER VALUE, BDZF
!                          IS A DUMMY VARIABLE. BDZF MUST BE
!                          DIMENSIONED AT LEAST (L+1)*(M+1).
!
!                        ELMBDA
!                          THE CONSTANT LAMBDA IN THE HELMHOLTZ
!                          EQUATION. IF LAMBDA .GT. 0, A SOLUTION
!                          MAY NOT EXIST.  HOWEVER, HW3CRT WILL
!                          ATTEMPT TO FIND A SOLUTION.
!
!                        LDIMF
!                          THE ROW (OR FIRST) DIMENSION OF THE
!                          ARRAYS F,BDYS,BDYF,BDZS,AND BDZF AS IT
!                          APPEARS IN THE PROGRAM CALLING HW3CRT.
!                          THIS PARAMETER IS USED TO SPECIFY THE
!                          VARIABLE DIMENSION OF THESE ARRAYS.
!                          LDIMF MUST BE AT LEAST L+1.
!
!                        MDIMF
!                          THE COLUMN (OR SECOND) DIMENSION OF THE
!                          ARRAY F AND THE ROW (OR FIRST) DIMENSION
!                          OF THE ARRAYS BDXS AND BDXF AS IT APPEARS
!                          IN THE PROGRAM CALLING HW3CRT.  THIS
!                          PARAMETER IS USED TO SPECIFY THE VARIABLE
!                          DIMENSION OF THESE ARRAYS.
!                          MDIMF MUST BE AT LEAST M+1.
!
!                        F
!                          A THREE-DIMENSIONAL ARRAY OF DIMENSION AT
!                          AT LEAST (L+1)*(M+1)*(N+1), SPECIFYING THE
!                          VALUES OF THE RIGHT SIDE OF THE HELMHOLZ
!                          EQUATION AND BOUNDARY VALUES (IF ANY).
!
!                          ON THE INTERIOR, F IS DEFINED AS FOLLOWS:
!                          FOR I=2,3,...,L,  J=2,3,...,M,
!                          AND K=2,3,...,N
!                          F(I,J,K) = F(X(I),Y(J),Z(K)).
!
!                          ON THE BOUNDARIES, F IS DEFINED AS FOLLOWS:
!                          FOR J=1,2,...,M+1,  K=1,2,...,N+1,
!                          AND I=1,2,...,L+1
!
!                          LBDCND      F(1,J,K)         F(L+1,J,K)
!                          ------   ---------------   ---------------
!
!                            0      F(XS,Y(J),Z(K))   F(XS,Y(J),Z(K))
!                            1      U(XS,Y(J),Z(K))   U(XF,Y(J),Z(K))
!                            2      U(XS,Y(J),Z(K))   F(XF,Y(J),Z(K))
!                            3      F(XS,Y(J),Z(K))   F(XF,Y(J),Z(K))
!                            4      F(XS,Y(J),Z(K))   U(XF,Y(J),Z(K))
!
!                          MBDCND      F(I,1,K)         F(I,M+1,K)
!                          ------   ---------------   ---------------
!
!                            0      F(X(I),YS,Z(K))   F(X(I),YS,Z(K))
!                            1      U(X(I),YS,Z(K))   U(X(I),YF,Z(K))
!                            2      U(X(I),YS,Z(K))   F(X(I),YF,Z(K))
!                            3      F(X(I),YS,Z(K))   F(X(I),YF,Z(K))
!                            4      F(X(I),YS,Z(K))   U(X(I),YF,Z(K))
!
!                          NBDCND      F(I,J,1)         F(I,J,N+1)
!                          ------   ---------------   ---------------
!
!                            0      F(X(I),Y(J),ZS)   F(X(I),Y(J),ZS)
!                            1      U(X(I),Y(J),ZS)   U(X(I),Y(J),ZF)
!                            2      U(X(I),Y(J),ZS)   F(X(I),Y(J),ZF)
!                            3      F(X(I),Y(J),ZS)   F(X(I),Y(J),ZF)
!                            4      F(X(I),Y(J),ZS)   U(X(I),Y(J),ZF)
!
!                          NOTE:
!                          IF THE TABLE CALLS FOR BOTH THE SOLUTION
!                          U AND THE RIGHT SIDE F ON A BOUNDARY,
!                          THEN THE SOLUTION MUST BE SPECIFIED.
!
!
! ON OUTPUT              F
!                          CONTAINS THE SOLUTION U(I,J,K) OF THE
!                          FINITE DIFFERENCE APPROXIMATION FOR THE
!                          GRID POINT (X(I),Y(J),Z(K)) FOR
!                          I=1,2,...,L+1, J=1,2,...,M+1,
!                          AND K=1,2,...,N+1.
!
!                        PERTRB
!                          IF A COMBINATION OF PERIODIC OR DERIVATIVE
!                          BOUNDARY CONDITIONS IS SPECIFIED FOR A
!                          POISSON EQUATION (LAMBDA = 0), A SOLUTION
!                          MAY NOT EXIST.  PERTRB IS A CONSTANT,
!                          CALCULATED AND SUBTRACTED FROM F, WHICH
!                          ENSURES THAT A SOLUTION EXISTS.  PWSCRT
!                         THEN COMPUTES THIS SOLUTION, WHICH IS A
!                          LEAST SQUARES SOLUTION TO THE ORIGINAL
!                          APPROXIMATION.  THIS SOLUTION IS NOT
!                          UNIQUE AND IS UNNORMALIZED.  THE VALUE OF
!                          PERTRB SHOULD BE SMALL COMPARED TO THE
!                          THE RIGHT SIDE F.  OTHERWISE, A SOLUTION
!                          IS OBTAINED TO AN ESSENTIALLY DIFFERENT
!                          PROBLEM.  THIS COMPARISON SHOULD ALWAYS
!                          BE MADE TO INSURE THAT A MEANINGFUL
!                          SOLUTION HAS BEEN OBTAINED.
!
!                        IERROR
!                          AN ERROR FLAG THAT INDICATES INVALID INPUT
!                          PARAMETERS.  EXCEPT FOR NUMBERS 0 AND 12,
!                          A SOLUTION IS NOT ATTEMPTED.
!
!                          =  0  NO ERROR
!                          =  1  XS .GE. XF
!                          =  2  L .LT. 5
!                          =  3  LBDCND .LT. 0 .OR. LBDCND .GT. 4
!                          =  4  YS .GE. YF
!                          =  5  M .LT. 5
!                          =  6  MBDCND .LT. 0 .OR. MBDCND .GT. 4
!                          =  7  ZS .GE. ZF
!                          =  8  N .LT. 5
!                          =  9  NBDCND .LT. 0 .OR. NBDCND .GT. 4
!                          = 10  LDIMF .LT. L+1
!                          = 11  MDIMF .LT. M+1
!                          = 12  LAMBDA .GT. 0
!                          = 20 If the dynamic allocation of real and
!                               complex work space required for solution
!                               fails (for example if N,M are too large
!                               for your computer)
!
!                          SINCE THIS IS THE ONLY MEANS OF INDICATING
!                          A POSSIBLY INCORRECT CALL TO HW3CRT, THE
!                          USER SHOULD TEST IERROR AFTER THE CALL.
!
! SPECIAL CONDITIONS     NONE
!
! I/O                    NONE
!
! PRECISION              SINGLE
!
! REQUIRED Files         fish.f,pois3d.f,fftpack.f,comf.f
!
! LANGUAGE               FORTRAN 90
!
! HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN THE LATE
!                        1970'S.  RELEASED ON NCAR'S PUBLIC SOFTWARE
!                        LIBRARIES IN JANUARY 1980.
!                        Revised in June 2004 by John Adams using
!                        Fortran 90 dynamically allocated work space.
!
! PORTABILITY            FORTRAN 90
!
! ALGORITHM              THIS SUBROUTINE DEFINES THE FINITE DIFFERENCE
!                        EQUATIONS, INCORPORATES BOUNDARY DATA, AND
!                        ADJUSTS THE RIGHT SIDE OF SINGULAR SYSTEMS AND
!                        THEN CALLS POIS3D TO SOLVE THE SYSTEM.
!
! TIMING                 FOR LARGE L, M AND N, THE OPERATION COUNT
!                        IS ROUGHLY PROPORTIONAL TO
!                          L*M*N*(LOG2(L)+LOG2(M)+5),
!                        BUT ALSO DEPENDS ON INPUT PARAMETERS LBDCND
!                        AND MBDCND.
!
! ACCURACY               THE SOLUTION PROCESS EMPLOYED RESULTS IN
!                        A LOSS OF NO MORE THAN FOUR SIGNIFICANT
!                        DIGITS FOR L, M AND N AS LARGE AS 32.
!                        MORE DETAILED INFORMATION ABOUT ACCURACY
!                        CAN BE FOUND IN THE DOCUMENTATION FOR
!                        ROUTINE POIS3D WHICH IS THE ROUTINE THAT
!                        ACTUALLY SOLVES THE FINITE DIFFERENCE
!                        EQUATIONS.
!
! REFERENCES             NONE
!***********************************************************************

SUBROUTINE HW3CRT(XS, XF, L, LBDCND, BDXS, BDXF, YS, YF, M, MBDCND, &
                 BDYS, BDYF, ZS, ZF, N, NBDCND, BDZS, BDZF, ELMBDA, LDIMF, &
                 MDIMF, F, PERTRB, IERROR)
      USE fish
      IMPLICIT REAL*8(A-H,O-Z)
      
      TYPE (fishworkspace) :: w
      
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER  :: L
      INTEGER  :: LBDCND
      INTEGER  :: M
      INTEGER  :: MBDCND
      INTEGER  :: N
      INTEGER  :: NBDCND
      INTEGER  :: LDIMF
      INTEGER  :: MDIMF
      INTEGER  :: IERROR
      REAL*8  :: XS
      REAL*8  :: XF
      REAL*8  :: YS
      REAL*8  :: YF
      REAL*8  :: ZS
      REAL*8  :: ZF
      REAL*8  :: ELMBDA
      REAL*8  :: PERTRB
      REAL*8  :: BDXS(MDIMF,*)
      REAL*8  :: BDXF(MDIMF,*)
      REAL*8  :: BDYS(LDIMF,*)
      REAL*8  :: BDYF(LDIMF,*)
      REAL*8  :: BDZS(LDIMF,*)
      REAL*8  :: BDZF(LDIMF,*)
      REAL*8  :: F(LDIMF,MDIMF,*)
      
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: IRWK, ICWK
!-----------------------------------------------
!
!     CHECK FOR INVALID INPUT.
!
      IERROR = 0
      IF (XF <= XS)               IERROR = 1
      IF (L < 5)                  IERROR = 2
      IF (LBDCND<0 .OR. LBDCND>4) IERROR = 3
      IF (YF <= YS)               IERROR = 4
      IF (M < 5)                  IERROR = 5
      IF (MBDCND<0 .OR. MBDCND>4) IERROR = 6
      IF (ZF <= ZS)               IERROR = 7
      IF (N < 5)                  IERROR = 8
      IF (NBDCND<0 .OR. NBDCND>4) IERROR = 9
      IF (LDIMF < L + 1)          IERROR = 10
      IF (MDIMF < M + 1)          IERROR = 11
      
      IF (IERROR /= 0) RETURN 
      
!     allocate required work space length (generous estimate)
      IRWK=30+L+M+5*N+MAX0(L,M,N)+7*(INT((L+1)/2)+INT((M+1)/2))
      ICWK = 0
      CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
!     check that allocation was successful
      IF (IERROR == 20) RETURN 
      call hw3crtt(XS,XF,L,LBDCND,BDXS,BDXF,YS,YF,M,MBDCND,BDYS, &
                   BDYF,ZS,ZF,N,NBDCND,BDZS,BDZF,ELMBDA,LDIMF, &
                   MDIMF,F,PERTRB,IERROR,w%rew)
!     release allocated work space
      CALL FISHFIN (W)
      RETURN 
END SUBROUTINE HW3CRT


 
SUBROUTINE HW3CRTT(XS, XF, L, LBDCND, BDXS, BDXF, YS, YF, M, &
                   MBDCND, BDYS, BDYF, ZS, ZF, N, NBDCND, BDZS, BDZF, ELMBDA,  &
                   LDIMF, MDIMF, F, PERTRB, IERROR, W)
      IMPLICIT REAL*8(A-H,O-Z)
      
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER , INTENT(IN) :: L
      INTEGER  :: LBDCND
      INTEGER , INTENT(IN) :: M
      INTEGER  :: MBDCND
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: NBDCND
      INTEGER  :: LDIMF
      INTEGER  :: MDIMF
      INTEGER , INTENT(OUT) :: IERROR
      REAL*8 , INTENT(IN) :: XS
      REAL*8 , INTENT(IN) :: XF
      REAL*8 , INTENT(IN) :: YS
      REAL*8 , INTENT(IN) :: YF
      REAL*8 , INTENT(IN) :: ZS
      REAL*8 , INTENT(IN) :: ZF
      REAL*8 , INTENT(IN) :: ELMBDA
      REAL*8 , INTENT(OUT) :: PERTRB
      REAL*8 , INTENT(IN) :: BDXS(MDIMF,*)
      REAL*8 , INTENT(IN) :: BDXF(MDIMF,*)
      REAL*8 , INTENT(IN) :: BDYS(LDIMF,*)
      REAL*8 , INTENT(IN) :: BDYF(LDIMF,*)
      REAL*8 , INTENT(IN) :: BDZS(LDIMF,*)
      REAL*8 , INTENT(IN) :: BDZF(LDIMF,*)
      REAL*8  :: F(LDIMF,MDIMF,*)
      REAL*8  :: W(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: MSTART, MSTOP, MP1, MP, MUNK, NP, NP1, NSTART, NSTOP, &
                 NUNK, LP1, LP, LSTART, LSTOP, J, K, LUNK, I, IWB, IWC, IWW, &
                 MSTPM1, LSTPM1, NSTPM1, NPEROD, IR
      REAL*8  :: DY,TWBYDY,C2,DZ,TWBYDZ,C3,DX,C1,TWBYDX,XLP,YLP,ZLP,S1,S2,S
 
      DY     = (YF - YS)/FLOAT(M)
      TWBYDY = 2.D0/DY
      C2     = 1.D0/DY**2
      MSTART = 1
      MSTOP  = M
      MP1    = M + 1
      MP     = MBDCND + 1
      
      GO TO (104,101,101,102,102) MP
  101 CONTINUE
      MSTART = 2
  102 CONTINUE
      GO TO (104,104,103,103,104) MP
  103 CONTINUE
      MSTOP = MP1
  104 CONTINUE
  
      MUNK   = MSTOP - MSTART + 1
      DZ     = (ZF - ZS)/FLOAT(N)
      TWBYDZ = 2.D0/DZ
      NP     = NBDCND + 1
      C3     = 1.D0/DZ**2
      NP1    = N + 1
      NSTART = 1
      NSTOP  = N
      
      GO TO (108,105,105,106,106) NP
  105 CONTINUE
      NSTART = 2
  106 CONTINUE
      GO TO (108,108,107,107,108) NP
  107 CONTINUE
      NSTOP = NP1
  108 CONTINUE
  
      NUNK   = NSTOP - NSTART + 1
      LP1    = L + 1
      DX     = (XF - XS)/FLOAT(L)
      C1     = 1.D0/DX**2
      TWBYDX = 2.D0/DX
      LP     = LBDCND + 1
      LSTART = 1
      LSTOP  = L
!
!     ENTER BOUNDARY DATA FOR X-BOUNDARIES.
!
      GO TO (122,109,109,112,112) LP
  109 CONTINUE
      LSTART = 2
      F(2,MSTART:MSTOP,NSTART:NSTOP) = F(2,MSTART:MSTOP,NSTART:NSTOP) - & 
                                       C1*F(1,MSTART:MSTOP,NSTART:NSTOP)
      GO TO 115
  112 CONTINUE
      F(1,MSTART:MSTOP,NSTART:NSTOP) = F(1,MSTART:MSTOP,NSTART:NSTOP) + &
                                       TWBYDX*BDXS(MSTART:MSTOP,NSTART:NSTOP)
  115 CONTINUE
      GO TO (122,116,119,119,116) LP
  116 CONTINUE
      F(L,MSTART:MSTOP,NSTART:NSTOP) = F(L,MSTART:MSTOP,NSTART:NSTOP) - &
                                       C1*F(LP1,MSTART:MSTOP,NSTART:NSTOP)
      GO TO 122
  119 CONTINUE
      LSTOP = LP1
      F(LP1,MSTART:MSTOP,NSTART:NSTOP) = F(LP1,MSTART:MSTOP,NSTART:NSTOP) - &
                                         TWBYDX*BDXF(MSTART:MSTOP,NSTART:NSTOP)
  122 CONTINUE
      LUNK = LSTOP - LSTART + 1
!
!     ENTER BOUNDARY DATA FOR Y-BOUNDARIES.
!
      GO TO (136,123,123,126,126) MP
  123 CONTINUE
      F(LSTART:LSTOP,2,NSTART:NSTOP) = F(LSTART:LSTOP,2,NSTART:NSTOP) - &
                                       C2*F(LSTART:LSTOP,1,NSTART:NSTOP)
      GO TO 129
  126 CONTINUE
      F(LSTART:LSTOP,1,NSTART:NSTOP) = F(LSTART:LSTOP,1,NSTART:NSTOP) + &
                                       TWBYDY*BDYS(LSTART:LSTOP,NSTART:NSTOP)
  129 CONTINUE
      GO TO (136,130,133,133,130) MP
  130 CONTINUE
      F(LSTART:LSTOP,M,NSTART:NSTOP) = F(LSTART:LSTOP,M,NSTART:NSTOP) - & 
                                       C2*F(LSTART:LSTOP,MP1,NSTART:NSTOP)
      GO TO 136
  133 CONTINUE
      F(LSTART:LSTOP,MP1,NSTART:NSTOP) = F(LSTART:LSTOP,MP1,NSTART:NSTOP) - &
                                         TWBYDY*BDYF(LSTART:LSTOP,NSTART:NSTOP)
  136 CONTINUE
      GO TO (150,137,137,140,140) NP
  137 CONTINUE
      F(LSTART:LSTOP,MSTART:MSTOP,2) = F(LSTART:LSTOP,MSTART:MSTOP,2) - &
                                       C3*F(LSTART:LSTOP,MSTART:MSTOP,1)
      GO TO 143
  140 CONTINUE
      F(LSTART:LSTOP,MSTART:MSTOP,1) = F(LSTART:LSTOP,MSTART:MSTOP,1) + &
                                       TWBYDZ*BDZS(LSTART:LSTOP,MSTART:MSTOP)
  143 CONTINUE
      GO TO (150,144,147,147,144) NP
  144 CONTINUE
      F(LSTART:LSTOP,MSTART:MSTOP,N) = F(LSTART:LSTOP,MSTART:MSTOP,N) - &
                                       C3*F(LSTART:LSTOP,MSTART:MSTOP,NP1)
      GO TO 150
  147 CONTINUE
      F(LSTART:LSTOP,MSTART:MSTOP,NP1) = F(LSTART:LSTOP,MSTART:MSTOP,NP1) - &
                                         TWBYDZ*BDZF(LSTART:LSTOP,MSTART:MSTOP)
     
!
!     DEFINE A,B,C COEFFICIENTS IN W-ARRAY.
!

  150 CONTINUE
      IWB = NUNK + 1
      IWC = IWB + NUNK
      IWW = IWC + NUNK
      W(:NUNK) = C3
      W(IWC:NUNK-1+IWC) = C3
      W(IWB:NUNK-1+IWB) = (-2.D0*C3) + ELMBDA
      GO TO (155,155,153,152,152) NP
  152 CONTINUE
      W(IWC) = 2.D0*C3
  153 CONTINUE
      GO TO (155,155,154,154,155) NP
  154 CONTINUE
      W(IWB-1) = 2.D0*C3
  155 CONTINUE
      PERTRB = 0.D0
!
!     FOR SINGULAR PROBLEMS ADJUST DATA TO INSURE A SOLUTION WILL EXIST.
!
      GO TO (156,172,172,156,172) LP
  156 CONTINUE
      GO TO (157,172,172,157,172) MP
  157 CONTINUE
      GO TO (158,172,172,158,172) NP
  158 CONTINUE
      IF (ELMBDA >= 0.D0) THEN
         IF (ELMBDA /= 0.D0) THEN
            IERROR = 12
         ELSE
            MSTPM1 = MSTOP - 1
            LSTPM1 = LSTOP - 1
            NSTPM1 = NSTOP - 1
            XLP = (2 + LP)/3
            YLP = (2 + MP)/3
            ZLP = (2 + NP)/3
            S1 = 0.D0
            DO K = 2, NSTPM1
               DO J = 2, MSTPM1
                  S1 = S1 + SUM(F(2:LSTPM1,J,K))
                  S1 = S1 + (F(1,J,K)+F(LSTOP,J,K))/XLP
               END DO
               S2 = SUM(F(2:LSTPM1,1,K)+F(2:LSTPM1,MSTOP,K))
               S2 = (S2 + (F(1,1,K) + F(1,MSTOP,K) + F(LSTOP,1,K) + &
                           F(LSTOP,MSTOP,K))/XLP)/YLP
               S1 = S1 + S2
            END DO
            S = (F(1,1,1) + F(LSTOP,1,1) + F(1,1,NSTOP) + F(LSTOP,1,NSTOP) + &
                 F(1,MSTOP,1) + F(LSTOP,MSTOP,1) + F(1,MSTOP,NSTOP) + &
                 F(LSTOP,MSTOP,NSTOP))/(XLP*YLP)
            DO J = 2, MSTPM1
               S = S + SUM(F(2:LSTPM1,J,1)+F(2:LSTPM1,J,NSTOP))
            END DO
            S2 = 0.D0
            S2 = SUM(F(2:LSTPM1,1,1) + F(2:LSTPM1,1,NSTOP) + F(2:LSTPM1,MSTOP,1) + &
                     F(2:LSTPM1,MSTOP,NSTOP))
            S = S2/YLP + S
            S2 = 0.D0
            S2 = SUM(F(1,2:MSTPM1,1) + F(1,2:MSTPM1,NSTOP) + F(LSTOP,2:MSTPM1,1) + &
                     F(LSTOP,2:MSTPM1,NSTOP))
            S = S2/XLP + S
            PERTRB = (S/ZLP + S1)/((FLOAT(LUNK + 1) - XLP)*(FLOAT(MUNK + 1) &
                     - YLP)*(FLOAT(NUNK + 1) - ZLP))
            F(:LUNK,:MUNK,:NUNK) = F(:LUNK,:MUNK,:NUNK) - PERTRB
         ENDIF
      ENDIF
  172 CONTINUE
      NPEROD = 0
      IF (NBDCND /= 0) THEN
         NPEROD = 1
         W(1) = 0.D0
         W(IWW-1) = 0.D0
      ENDIF
      CALL POIS3DD (LBDCND, LUNK, C1, MBDCND, MUNK, C2, NPEROD, NUNK, W, &
                    W(IWB), W(IWC), LDIMF, MDIMF, F(LSTART,MSTART,NSTART), &
                    IR, W(IWW))
!
!     FILL IN SIDES FOR PERIODIC BOUNDARY CONDITIONS.
!
      IF (LP == 1) THEN
         IF (MP == 1) THEN
            F(1,MP1,NSTART:NSTOP) = F(1,1,NSTART:NSTOP)
            MSTOP = MP1
         ENDIF
         IF (NP == 1) THEN
            F(1,MSTART:MSTOP,NP1) = F(1,MSTART:MSTOP,1)
            NSTOP = NP1
         ENDIF
         F(LP1,MSTART:MSTOP,NSTART:NSTOP) = F(1,MSTART:MSTOP,NSTART:NSTOP)
      ENDIF
      IF (MP == 1) THEN
         IF (NP == 1) THEN
            F(LSTART:LSTOP,1,NP1) = F(LSTART:LSTOP,1,1)
            NSTOP = NP1
         ENDIF
         F(LSTART:LSTOP,MP1,NSTART:NSTOP) = F(LSTART:LSTOP,1,NSTART:NSTOP)
      ENDIF
      IF (NP == 1) THEN
         F(LSTART:LSTOP,MSTART:MSTOP,NP1) = F(LSTART:LSTOP,MSTART:MSTOP,1)
      ENDIF
      
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
! June      2004    Version 5.0, Fortran 90 changes
!-----------------------------------------------------------------------

END SUBROUTINE HW3CRTT
