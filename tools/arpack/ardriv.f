c**************************************************************************
c*  @file    readhb.f
c*  @author  Nathan Baker
c*  @brief   Routines for reading Harwell-Boeing format matrices
c*  @version $Id$
c*  @attention
c*  @verbatim
c*
c* APBS -- Adaptive Poisson-Boltzmann Solver
c*
c* Nathan A. Baker (baker@biochem.wustl.edu)
c* Dept. of Biochemistry and Molecular Biophysics
c* Center for Computational Biology
c* Washington University in St. Louis
c*
c* Additional contributing authors listed in the code documentation.
c*
c* Copyright (c) 2002-2006.  Washington University in St. Louis.
c* All Rights Reserved.
c* Portions Copyright (c) 1999-2002.  The Regents of the University 
c* of California.  
c* Portions Copyright (c) 1995.  Michael Holst.
c*
c* This file is part of APBS.
c*
c* APBS is free software; you can redistribute it and/or modify
c* it under the terms of the GNU General Public License as published
c* by c* the Free Software Foundation; either version 2 of the 
c* License, or (at your option) any later version.
c*
c* APBS is distributed in the hope that it will be useful,
c* but WITHOUT ANY WARRANTY; without even the implied warranty of
c* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c* GNU General Public License for more details.
c*
c* You should have received a copy of the GNU General Public License
c* along with APBS; if not, write to the Free Software
c* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 
c* 02111-1307  USA
c* * Linking APBS statically or dynamically with other modules is making a
c* combined work based on APBS. Thus, the terms and conditions of the GNU
c* General Public License cover the whole combination.
c* 
c* SPECIAL GPL EXCEPTION
c* In addition, as a special exception, the copyright holders of APBS
c* give you permission to combine the APBS program with free software
c* programs and libraries that are released under the GNU LGPL or with
c* code included in releases of ISIM, Ion Simulator Interface, PMV, PyMOL
c* SMOL, VMD, and Vision. Such combined software may be linked with APBS and 
c* redistributed together in original or modified form as mere aggregation
c* without requirement that the entire work be under the scope of the GNU 
c* General Public License. This special exception permission is also extended
c* to any software listed in the SPECIAL GPL EXCEPTION clauses by the PMG,
c* FEtk, MC, or MALOC libraries.
c* 
c* Note that people who make modified versions of APBS are not obligated
c* to grant this special exception for their modified versions; it is
c* their choice whether to do so. The GNU General Public License gives
c* permission to release a modified version without this exception; this
c* exception also makes it possible to release a modified version which
c* carries forward this exception.
c*
c* @endverbatim
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
