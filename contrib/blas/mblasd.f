c* ///////////////////////////////////////////////////////////////////////////
c* PMG -- Parallel algebraic MultiGrid
c* Copyright (c) 1994-2008.  Michael Holst.
c*
c* Michael Holst <mholst@math.ucsd.edu>
c* University of California, San Diego
c* Department of Mathematics, 5739 AP&M
c* 9500 Gilman Drive, Dept. 0112
c* La Jolla, CA 92093-0112 USA                                                 
c* http://math.ucsd.edu/~mholst
c*
c* This file is part of PMG.
c*
c* This library is free software; you can redistribute it and/or
c* modify it under the terms of the GNU Lesser General Public
c* License as published by the Free Software Foundation; either
c* version 2.1 of the License, or (at your option) any later version.
c*
c* This library is distributed in the hope that it will be useful,
c* but WITHOUT ANY WARRANTY; without even the implied warranty of
c* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
c* Lesser General Public License for more details.

c* You should have received a copy of the GNU Lesser General Public
c* License along with this library; if not, write to the Free Software
c* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
c* ///////////////////////////////////////////////////////////////////////////

      subroutine daxpy(n,alpha,x,istep,y,jstep)
c* *********************************************************************
c* purpose:
c*
c*    parallel/vector vector saxpy.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      double precision x(*),y(*),alpha
      integer          n,istep,jstep,i,ii,nproc,ipara,ivect
      parameter        (nproc=1)
c*
cmdir 0 0
c*
c*    *** find parallel loops (ipara), remainder (ivect) ***
      ipara = n / nproc
      ivect = mod(n,nproc)
c*
c*    *** do parallel loops ***
cmdir 2 1
      do 10 ii = 1, nproc
cmdir 2 2
         do 11 i = 1+(ipara*(ii-1)), ipara*ii
            y(i) = y(i) + alpha * x(i)
 11      continue
 10   continue
c*
c*    *** do vector loops ***
cmdir 1 1
      do 20 i = ipara*nproc+1, n
         y(i) = y(i) + alpha * x(i)
 20   continue
c*
c*    *** return and end ***
      return
      end
      subroutine dcopy(n,x,istep,y,jstep)
c* *********************************************************************
c* purpose:
c*
c*    parallel/vector vector copy.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      double precision x(*),y(*)
      integer          n,istep,jstep,i,ii,nproc,ipara,ivect
      parameter        (nproc=1)
c*
cmdir 0 0
c*
c*    *** find parallel loops (ipara), remainder (ivect) ***
      ipara = n / nproc
      ivect = mod(n,nproc)
c*
c*    *** do parallel loops ***
cmdir 2 1
      do 10 ii = 1, nproc
cmdir 2 2
         do 11 i = 1+(ipara*(ii-1)), ipara*ii
            y(i) = x(i) 
 11      continue
 10   continue
c*
c*    *** do vector loops ***
cmdir 1 1
      do 20 i = ipara*nproc+1, n
         y(i) = x(i) 
 20   continue
c*
c*    *** return and end ***
      return
      end
      function dasum(n,x,istep)
c* *********************************************************************
c* purpose:
c*
c*    parallel/vector sum of vector components.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      double precision x(*)
      double precision dasum,dnrm1
      integer          n,istep
c*
cmdir 0 0
c*
c*    *** call the max-norm function ***
      dasum = dnrm1(n,x,istep)
c*
c*    *** return and end ***
      return
      end
      function dnrm1(n,x,istep)
c* *********************************************************************
c* purpose:
c*
c*    parallel/vector vector norm.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      double precision x(*)
      double precision tmp,dnrm1
      integer          n,istep,i,ii,nproc,ipara,ivect
      parameter        (nproc=1)
c*
cmdir 0 0
c*
c*    *** initialize ***
      tmp = 0.0d0
c*
c*    *** find parallel loops (ipara), remainder (ivect) ***
      ipara = n / nproc
      ivect = mod(n,nproc)
c*
c*    *** do parallel loops ***
cmdir 2 1
      do 10 ii = 1, nproc
cmdir 2 2
         do 11 i = 1+(ipara*(ii-1)), ipara*ii
            tmp = tmp + dabs(x(i))
 11      continue
 10   continue
c*
c*    *** do vector loops ***
cmdir 1 1
      do 20 i = ipara*nproc+1, n
         tmp = tmp + dabs(x(i))
 20   continue
c*
c*    *** finish up ***
      dnrm1 = tmp
c*
c*    *** return and end ***
      return
      end
      function dnrm2(n,x,istep)
c* *********************************************************************
c* purpose:
c*
c*    parallel/vector vector norm.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      double precision x(*)
      double precision tmp,dnrm2
      integer          n,istep,i,ii,nproc,ipara,ivect
      parameter        (nproc=1)
c*
cmdir 0 0
c*
c*    *** initialize ***
      tmp = 0.0d0
c*
c*    *** find parallel loops (ipara), remainder (ivect) ***
      ipara = n / nproc
      ivect = mod(n,nproc)
c*
c*    *** do parallel loops ***
cmdir 2 1
      do 10 ii = 1, nproc
cmdir 2 2
         do 11 i = 1+(ipara*(ii-1)), ipara*ii
            tmp = tmp + x(i)*x(i)
 11      continue
 10   continue
c*
c*    *** do vector loops ***
cmdir 1 1
      do 20 i = ipara*nproc+1, n
         tmp = tmp + x(i)*x(i)
 20   continue
c*
c*    *** finish up ***
      dnrm2 = dsqrt(tmp)
c*
c*    *** return and end ***
      return
      end
      function dnrm8(n,x,istep)
c* *********************************************************************
c* purpose:
c*
c*    parallel/vector vector norm.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      double precision x(*)
      double precision tmp,dnrm8
      integer          n,istep,i,ii,nproc,ipara,ivect
      parameter        (nproc=1)
c*
cmdir 0 0
c*
c*    *** initialize ***
      tmp = 0.0d0
c*
c*    *** find parallel loops (ipara), remainder (ivect) ***
      ipara = n / nproc
      ivect = mod(n,nproc)
c*
c*    *** do parallel loops ***
cmdir 2 1
      do 10 ii = 1, nproc
cmdir 2 2
         do 11 i = 1+(ipara*(ii-1)), ipara*ii
            tmp = dmax1( tmp , dabs(x(i)) )
 11      continue
 10   continue
c*
c*    *** do vector loops ***
cmdir 1 1
      do 20 i = ipara*nproc+1, n
         tmp = dmax1( tmp , dabs(x(i)) )
 20   continue
c*
c*    *** finish up ***
      dnrm8 = tmp
c*
c*    *** return and end ***
      return
      end
      subroutine dscal(n,fac,x,istep)
c* *********************************************************************
c* purpose:
c*
c*    parallel/vector vector scale.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      double precision x(*),fac
      integer          n,istep,i,ii,nproc,ipara,ivect
      parameter        (nproc=1)
c*
cmdir 0 0
c*
c*    *** find parallel loops (ipara), remainder (ivect) ***
      ipara = n / nproc
      ivect = mod(n,nproc)
c*
c*    *** do parallel loops ***
cmdir 2 1
      do 10 ii = 1, nproc
cmdir 2 2
         do 11 i = 1+(ipara*(ii-1)), ipara*ii
            x(i) = fac*x(i)
 11      continue
 10   continue
c*
c*    *** do vector loops ***
cmdir 1 1
      do 20 i = ipara*nproc+1, n
         x(i) = fac*x(i)
 20   continue
c*
c*    *** return and end ***
      return
      end
      function ddot(n,x,istep,y,jstep)
c* *********************************************************************
c* purpose:
c*
c*    parallel/vector dot product.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      double precision x(*),y(*),ddot
      integer          n,istep,jstep,i,ii,nproc,ipara,ivect
      parameter        (nproc=1)
c*
cmdir 0 0
c*
c*    *** initialize ***
      ddot = 0.0d0
c*
c*    *** find parallel loops (ipara), remainder (ivect) ***
      ipara = n / nproc
      ivect = mod(n,nproc)
c*
c*    *** do parallel loops ***
cmdir 2 1
      do 10 ii = 1, nproc
cmdir 2 2
         do 11 i = 1+(ipara*(ii-1)), ipara*ii
            ddot = ddot + x(i)*y(i)
 11      continue
 10   continue
c*
c*    *** do vector loops ***
cmdir 1 1
      do 20 i = ipara*nproc+1, n
         ddot = ddot + x(i)*y(i)
 20   continue
c*
c*    *** return and end ***
      return
      end
      function idamax(n,sx,incx)
c* *********************************************************************
c* purpose:
c*
c*    parallel/vector icamax.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          idamax,n,ns,incx,ii,i
      double precision smax,xmag
      double precision sx(*)
c*
cmdir 0 0
c*
      idamax = 0
      if(n.le.0) return
      idamax = 1
      if(n.le.1)return
      if(incx.eq.1)goto 20
c*
c*       code for increments not equal to 1.
c*
      smax = dabs(sx(1))
      ns = n*incx
      ii = 1
      do 10 i=1,ns,incx
          xmag = dabs(sx(i))
          if(xmag.le.smax) go to 40
          idamax = ii
          smax = xmag
 40       ii = ii + 1
 10   continue
      return
c*
c*    code for increments equal to 1.
c*
 20   smax = dabs(sx(1))
      do 30 i = 2,n
         xmag = dabs(sx(i))
         if(xmag.le.smax) go to 30
         idamax = i
         smax = xmag
 30   continue
      return
      end

      SUBROUTINE DTRSV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
c* *********************************************************************

*     .. Scalar Arguments ..
      INTEGER            INCX, LDA, N
      CHARACTER*1        DIAG, TRANS, UPLO
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  DTRSV  solves one of the systems of equations
*
*     A*x = b,   or   A'*x = b,
*
*  where b and x are n element vectors and A is an n by n unit, or
*  non-unit, upper or lower triangular matrix.
*
*  No test for singularity or near-singularity is included in this
*  routine. Such tests must be performed before calling this routine.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the equations to be solved as
*           follows:
*
*              TRANS = 'N' or 'n'   A*x = b.
*
*              TRANS = 'T' or 't'   A'*x = b.
*
*              TRANS = 'C' or 'c'   A'*x = b.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit
*           triangular as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U' or 'u', the leading n by n
*           upper triangular part of the array A must contain the upper
*           triangular matrix and the strictly lower triangular part of
*           A is not referenced.
*           Before entry with UPLO = 'L' or 'l', the leading n by n
*           lower triangular part of the array A must contain the lower
*           triangular matrix and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u', the diagonal elements of
*           A are not referenced either, but are assumed to be unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, n ).
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element right-hand side vector b. On exit, X is overwritten
*           with the solution vector x.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
*     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JX, KX
      LOGICAL            NOUNIT
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( UPLO , 'U' ).AND.
     $         .NOT.LSAME( UPLO , 'L' )      )THEN
         INFO = 1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 2
      ELSE IF( .NOT.LSAME( DIAG , 'U' ).AND.
     $         .NOT.LSAME( DIAG , 'N' )      )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTRSV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
      NOUNIT = LSAME( DIAG, 'N' )
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX  too small for descending loops.
*
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  x := inv( A )*x.
*
         IF( LSAME( UPLO, 'U' ) )THEN
            IF( INCX.EQ.1 )THEN
               DO 20, J = N, 1, -1
                  IF( X( J ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( J ) = X( J )/A( J, J )
                     TEMP = X( J )
                     DO 10, I = J - 1, 1, -1
                        X( I ) = X( I ) - TEMP*A( I, J )
   10                CONTINUE
                  END IF
   20          CONTINUE
            ELSE
               JX = KX + ( N - 1 )*INCX
               DO 40, J = N, 1, -1
                  IF( X( JX ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )/A( J, J )
                     TEMP = X( JX )
                     IX   = JX
                     DO 30, I = J - 1, 1, -1
                        IX      = IX      - INCX
                        X( IX ) = X( IX ) - TEMP*A( I, J )
   30                CONTINUE
                  END IF
                  JX = JX - INCX
   40          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 60, J = 1, N
                  IF( X( J ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( J ) = X( J )/A( J, J )
                     TEMP = X( J )
                     DO 50, I = J + 1, N
                        X( I ) = X( I ) - TEMP*A( I, J )
   50                CONTINUE
                  END IF
   60          CONTINUE
            ELSE
               JX = KX
               DO 80, J = 1, N
                  IF( X( JX ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )/A( J, J )
                     TEMP = X( JX )
                     IX   = JX
                     DO 70, I = J + 1, N
                        IX      = IX      + INCX
                        X( IX ) = X( IX ) - TEMP*A( I, J )
   70                CONTINUE
                  END IF
                  JX = JX + INCX
   80          CONTINUE
            END IF
         END IF
      ELSE
*
*        Form  x := inv( A' )*x.
*
         IF( LSAME( UPLO, 'U' ) )THEN
            IF( INCX.EQ.1 )THEN
               DO 100, J = 1, N
                  TEMP = X( J )
                  DO 90, I = 1, J - 1
                     TEMP = TEMP - A( I, J )*X( I )
   90             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( J, J )
                  X( J ) = TEMP
  100          CONTINUE
            ELSE
               JX = KX
               DO 120, J = 1, N
                  TEMP = X( JX )
                  IX   = KX
                  DO 110, I = 1, J - 1
                     TEMP = TEMP - A( I, J )*X( IX )
                     IX   = IX   + INCX
  110             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( J, J )
                  X( JX ) = TEMP
                  JX      = JX   + INCX
  120          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 140, J = N, 1, -1
                  TEMP = X( J )
                  DO 130, I = N, J + 1, -1
                     TEMP = TEMP - A( I, J )*X( I )
  130             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( J, J )
                  X( J ) = TEMP
  140          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 160, J = N, 1, -1
                  TEMP = X( JX )
                  IX   = KX
                  DO 150, I = N, J + 1, -1
                     TEMP = TEMP - A( I, J )*X( IX )
                     IX   = IX   - INCX
  150             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( J, J )
                  X( JX ) = TEMP
                  JX      = JX   - INCX
  160          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of DTRSV .
*
      END

      SUBROUTINE XERBLA( SRNAME, INFO )
c* *********************************************************************
*
*  -- LAPACK auxiliary routine (preliminary version) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER*6        SRNAME
      INTEGER            INFO
*     ..
*
*  Purpose
*  =======
*
*  XERBLA  is an error handler for the LAPACK routines.
*  It is called by an LAPACK routine if an input parameter has an
*  invalid value.  A message is printed and execution stops.
*
*  Installers may consider modifying the STOP statement in order to
*  call system-specific exception-handling facilities.
*
*  Arguments
*  =========
*
*  SRNAME  (input) CHARACTER*6
*          The name of the routine which called XERBLA.
*
*  INFO    (input) INTEGER
*          The position of the invalid parameter in the parameter list
*          of the calling routine.
*
*
      WRITE( *, FMT = 9999 )SRNAME, INFO
*
      STOP
*
 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ',
     $      'an illegal value' )
*
*     End of XERBLA
*
      END

      LOGICAL          FUNCTION LSAME( CA, CB )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     January 31, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          CA, CB
*     ..
*
*  Purpose
*  =======
*
*  LSAME returns .TRUE. if CA is the same letter as CB regardless of
*  case.
*
*  Arguments
*  =========
*
*  CA      (input) CHARACTER*1
*  CB      (input) CHARACTER*1
*          CA and CB specify the single characters to be compared.
*
* =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
*     ..
*     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
*     ..
*     .. Executable Statements ..
*
*     Test if the characters are equal
*
      LSAME = CA.EQ.CB
      IF( LSAME )
     $   RETURN
*
*     Now test for equivalence if both characters are alphabetic.
*
      ZCODE = ICHAR( 'Z' )
*
*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
*     machines, on which ICHAR returns a value with bit 8 set.
*     ICHAR('A') on Prime machines returns 193 which is the same as
*     ICHAR('A') on an EBCDIC machine.
*
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
*
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
*
*        ASCII is assumed - ZCODE is the ASCII code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
*
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
*
*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.
     $       INTA.GE.145 .AND. INTA.LE.153 .OR.
     $       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.
     $       INTB.GE.145 .AND. INTB.LE.153 .OR.
     $       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
*
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
*
*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
*        plus 128 of either lower or upper case 'Z'.
*
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
*
*     RETURN
*
*     End of LSAME
*
      END
