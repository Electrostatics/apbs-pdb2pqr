c**************************************************************************
c*  @file    readhb.f
c*  @author  Nathan Baker
c*  @brief   Routines for reading Harwell-Boeing format matrices
c*  @version $Id$
c*  @attention
c*
c**************************************************************************
      subroutine geteigs(nrow  , ncol , nnzero, values, colptr, 
     1                   rowind, nev   , ncv  , key   , v     , 
     2                   workl, workd , d     , resid , ax    ,
     3                   select)
c**************************************************************************
c*  Author:  Nathan Baker (based on ARPACK examples)
c*  Arguments:  nrow   Number of rows in matrix (input)
c*              ncol   Number of columns in matrix (input)
c*              nnzero Number of nonzeros in matrix (input)
c*              values Nonzeros of matrix (input)
c*              colptr Column pointers for matrix (input)
c*              rowind Row pointers for matrix (input)
c*              nev    Number of eigenvalues to calculate (input)
c*              ncv    Number of Lanczos basis vectors to use (input)
c*              key    Which eigenvalues to calculate: (input)
c*                       0 => Largest
c*                       1 => Smallest
c*                       2 => Both ends of spectrum
c*              v      Lanzcos vector or eigenvector array (output) 
c*                     Requires nrow*ncv double storage
c*              workl  Work array (output)
c*                     Requires ncv*(ncv+8) double storage
c*              workd  Work array (output)
c*                     Requires 3*ncol double storage
c*              d      First column contains eigenvalues (output)
c*                     Requires 2*ncv double storage
c*              resid  Residual values (output)
c*                     Requires nrow double storage
c*              ax     Used to compute residual norm (output)
c*                     Requires nrow double storage
c*              select Dunno (output)
c*                     Requires ncv integer storage
c*********************************************************************


      implicit none

c*  Arguments
      integer          nrow, ncol, nnzero, key, ncv, nev,
     1                 colptr(ncol+1), rowind(nnzero)
      double precision values(nnzero), v(ncol,ncv), 
     1                 workl(ncv*(ncv+8)), workd(3*nrow),
     2                 d(ncv,2), resid(ncol), ax(ncol)
      logical          select(ncol)

c* Local variables
      integer          iparam(11), ipntr(11), ido, n, lworkl, info, 
     1                 ierr, i, j, maxitr, mode, ishfts, nconv
      character        bmat*1, which*2
      logical          rvec
      double precision tol, sigma, mat

c* BLAS & LAPACK
      double precision dnrm2
c*    external         dnrm2, daxpy

c* Here we go!

      if (nev.gt.nrow) then
          print *, 'ARDRIV: ERROR -- nev > nrow!'
          stop 'ARDRIV'
      endif
      if (nrow .ne. ncol) then
          print *, 'ARDRIV: only works with square matrices!'
          stop 'ARDRIV'
      endif
      if (ncv .lt. (nev + 1)) then
          print *, 'ARDRIV: ERROR -- ncv < (nev+1)!!!'
          stop 'ARDRIV'
      endif
      n = nrow

c* Do a standard eigenvalue calculation...
      bmat = 'I'

c* Check which eigenvalues we want
      if (key .eq. 0) then
          which = 'LM'
      elseif (key .eq. 1) then
          which = 'SM'
      elseif (key .eq. 2) then
          which = 'BE'
      else
          print *, 'ARDRIV: ERROR:  INVALID KEY: ', key
          stop 'ARDRIV'
      endif 

c* Set some values
      lworkl = ncv*(ncv+8)
      tol = 1.0D-9
      info = 0
      ido = 0

c* Set shift strategy (1 = exact)
      iparam(1) = 1

c* Set maximum number of Arnoldi iterations allowed
      maxitr = 999
      iparam(3) = maxitr

c* Set Mode 1 calculation (standard eigenvalue problem; see DSAUPD
c* in ARPACK documentationa)
      iparam(7) = 1

c* OOOH -- my favorite FORTRAN coding practice:  GOTO loops
 10   continue

c* Call DSAUPD and do as instructed by ido parameter
      call dsaupd( ido, bmat, n, which, nev, tol, resid,
     1             ncv, v, n, iparam, ipntr, workd, workl,
     2             lworkl, info )
      if ((ido .eq. -1) .or. (ido .eq. 1)) then

c* Do matrix-vector multiply
          call av(nrow, ncol, nnzero, values, colptr, rowind, 
     1            workd(ipntr(1)), workd(ipntr(2)))

c* OOOH -- my favorite FORTRAN coding practice:  GOTO loops
          go to 10
      endif

c* Check for convergence
      if (info .lt. 0) then
          print *, 'ARDRIV: DSAUPD ERROR, INFO = ', info
c* No errors => post-process
      else
          rvec = .true.
          call dseupd ( rvec, 'All', select, d, v, n, sigma,
     1                  bmat, n, which, nev, tol, resid, ncv, v, n,
     2                  iparam, ipntr, workd, workl, lworkl, ierr )

          if (ierr .ne. 0) then
              print *, 'ARDRIV: DSEUPD ERROR, IERR = ', ierr
c* Eigenvalues in the first column of d, eigenvectors in the first 
c* iparam(5) columns of v
          else         
              nconv = iparam(5)
              do 20 j=1, nconv
c* Compute residual norm
                  call av(nrow, ncol, nnzero, values, colptr, rowind, 
     1                    v(1,j), ax)
                  call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                  d(j,2) = dnrm2(n, ax, 1)
                  d(j,2) = d(j,2) / abs(d(j,1))
 20           continue
c* Print results
              print *, 'ARDRIV:  Printing eigenvalues'
              call dmout(6, nconv, 2, d, ncv, -6,
     1            'Ritz values and relative residuals')
          end if

c* Other info
          if (info .eq. 1) then
              print *, 'ARDRIV: WARNING!  Max iterations reached!'
              print *, 'ARDRIV: Max iterations = ', maxitr
          else if (info .eq. 3) then
              print *, 'ARDRIV: WARNING!  No shifts applied during'
              print *, 'ARDRIV: implicit Arnoldi update.  Increase'
              print *, 'ARDRIV: NCV'
          endif
          print *, 'ARDRIV:  # eigenvalues calculated = ', nev
          print *, 'ARDRIV:  # eigenvectors calculated = ', ncv
          print *, 'ARDRIV:  Portion of spectrum = ', which
          print *, 'ARDRIV:  # of iterations = ', iparam(3)
          print *, 'ARDRIV:  Convergence criterion = ', tol
      endif

      return 
      end

      subroutine av(nrow , ncol , nnzero, values, colptr, rowind,
     1              v    , w)
c**************************************************************************
c*  Author:  Nathan Baker (based on ARPACK examples)
c*  Arguments:  nrow   Number of rows in matrix (input)
c*              ncol   Number of columns in matrix (input)
c*              nnzero Number of nonzeros in matrix (input)
c*              values Nonzeros of matrix (input)
c*              colptr Column pointers for matrix (input)
c*              rowind Row pointers for matrix (input)
c*              v      Vector to muliply matrix by (input)
c*              w      Result of matrix-vector multiply (output)
c**************************************************************************

      implicit none

c*  Arguments
      integer          nrow     , ncol     , nnzero  , key,
     1                 colptr(ncol+1), rowind(nnzero)
      double precision values(nnzero), v(nrow), w(nrow), mat

c* Local variables
      integer          i, icol, irow

c* Here we go!

      if (nrow .ne. ncol) then
          print *, 'AV:  Only works for symmetric, square matrices!'
          stop 'AV1'
      endif

      do 10 i = 1, nrow
          w(i) = 0.0
 10   continue

      do 20 icol = 1, ncol
          do 30 i = colptr(icol), (colptr(icol+1)-1)
              irow = rowind(i)
              mat = values(i)
              w(irow) = w(irow) + mat*v(icol)
 30       continue
 20   continue

      return 
      end
