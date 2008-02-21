c* ///////////////////////////////////////////////////////////////////////////
c* MG/XMG -- Multilevel nonlinear scalar elliptic PDE solver and X interface
c* Copyright (C) 1995  Michael Holst
c*
c* This program is free software; you can redistribute it and/or modify
c* it under the terms of the GNU General Public License as published by
c* the Free Software Foundation; either version 2 of the License, or
c* (at your option) any later version.
c*
c* This program is distributed in the hope that it will be useful,
c* but WITHOUT ANY WARRANTY; without even the implied warranty of
c* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c* GNU General Public License for more details.
c*
c* You should have received a copy of the GNU General Public License
c* along with this program; if not, write to the Free Software
c* Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
c*
c* MG/XMG was developed by:
c*
c*    Michael Holst                TELE:  (858) 534-4899
c*    Department of Mathematics    FAX:   (858) 534-5273
c*    UC San Diego, AP&M 5739      EMAIL: mholst@math.ucsd.edu
c*    La Jolla, CA 92093 USA       WEB:   http://www.scicomp.ucsd.edu/~mholst
c*
c* See the file "maind.f" for more information and pointers to papers.
c*
c* IMPORTANT: If you intend to use or modify this code, make sure you 
c* understand your responsibilities under the GNU license.
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

