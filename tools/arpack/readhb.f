c**************************************************************************
c*  @file    arstuff.f
c*  @author  Nathan Baker
c*  @brief   Support routines for ARPACK
c*  @version $Id$
c*  @attention
c*  @verbatim
c*
c* APBS -- Adaptive Poisson-Boltzmann Solver
c*
c* Nathan A. Baker (nbaker@wasabi.ucsd.edu)
c* Dept. of Chemistry and Biochemistry
c* University of California, San Diego
c*
c* Additional contributing authors listed in the code documentation.
c*
c* Copyright (c) 1999-2002.  Nathan A. Baker
c*
c* Permission to use, copy, modify, and distribute this software and its
c* documentation for educational, research, and not-for-profit purposes,
c* without fee and without a signed licensing agreement, is hereby granted,
c* provided that the above copyright notice, this paragraph and the
c* following two paragraphs appear in all copies, modifications, and
c* distributions.
c*
c* IN NO EVENT SHALL THE AUTHOR BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
c* SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
c* ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
c* THE AUTHOR HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
c*
c* REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
c* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
c* PARTICULAR PURPOSE.  THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF
c* ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS".  THE AUTHOR HAS NO OBLIGATION
c* TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
c* MODIFICATIONS.
c*
c* @endverbatim
c**************************************************************************
      subroutine readhbhead(nrow, ncol, nnzero, path)
c**************************************************************************
c* Routine:  readhbhead
c* Purpose:  Read the header to a Harwell-Boeing format matrix file so you can
c*           allocate the necessary array space
c* Args:     nrow   -- Set to number of rows
c*           ncol   -- Set to number of columns
c*           nnzero -- Set to number of non-zeros
c* Author:   Nathan Baker (mostly ripped off from Harwell-Boeing format
c*           documentation)
c**************************************************************************
      implicit none
      character title*72 , path*256 , key*8    , mxtype*3 ,
     1          ptrfmt*16, indfmt*16, valfmt*20, rhsfmt*20
      integer   totcrd   , ptrcrd   , indcrd   , valcrd   ,
     1          nrow     , ncol     , nnzero  , neltvl,
     2          i        , myunit   , openstat, closestat,
     3          rhscrd  

c* Open the file for reading
      myunit = 12
      open(UNIT=myunit, FILE=path, STATUS="unknown")

c* Read in header block
      read(myunit, 1000) title , key   ,
     1                   totcrd, ptrcrd, indcrd, valcrd, rhscrd,
     2                   mxtype, nrow  , ncol  , nnzero, neltvl,
     3                   ptrfmt, indfmt, valfmt, rhsfmt
 1000 format(a72, a8 / 5i14 / a3, 11x, 4i14 / 2a16, 2a20)

c* Close the file
      close (UNIT=myunit)

c* Enough of this FORTRAN crap!
      return
      end

      subroutine readhb(nrow, ncol, nnzero, values, rowind, colptr,
     2  path, title, key, mxtype)
c**************************************************************************
c* Routine:  readhbhead
c* Purpose:  Read a a Harwell-Boeing format matrix file 
c* Args:     nrow   -- Set to number of rows
c*           ncol   -- Set to number of columns
c*           nnzero -- Set to number of non-zeros
c* Author:   Nathan Baker (mostly ripped off from Harwell-Boeing format
c*           documentation)
c**************************************************************************
      implicit none
      character        title*72 , path*256 , key*8    , mxtype*3 ,
     1                 ptrfmt*16, indfmt*16, valfmt*20, rhsfmt*20
      integer          totcrd   , ptrcrd   , indcrd  , valcrd,
     1                 nrow     , ncol     , nnzero  , neltvl,
     2                 i        , myunit   , openstat, closestat,
     3                 colptr(ncol+1), rowind(nnzero), rhscrd
      double precision values(nnzero)

c* Open the file for reading
      myunit = 12
      open(UNIT=myunit, FILE=path, STATUS="unknown")

c* Read in header block
      read(myunit, 1000) title , key   ,
     1                   totcrd, ptrcrd, indcrd, valcrd, rhscrd,
     2                   mxtype, nrow  , ncol  , nnzero, neltvl,
     3                   ptrfmt, indfmt, valfmt, rhsfmt
 1000 format(a72, a8 / 5i14 / a3, 11x, 4i14 / 2a16, 2a20)

c* Read matrix structure
      read(myunit, ptrfmt) (colptr(i), i = 1, ncol+1)
      read(myunit, indfmt) (rowind(i), i = 1, nnzero)

c* Read matrix values
      if (valcrd.gt.0) then
          read(myunit, valfmt) (values(i), i = 1, nnzero)
      endif

      return 
      end

