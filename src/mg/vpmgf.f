c**************************************************************************
c*  @file    vpmgf.f
c*  @ingroup Vpmg
c*  @author  Nathan Baker
c*  @brief   Class Vpmg FORTRAN routines
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
c* Copyright (c) 1999-2002. The Regents of the University of California
c*                          (Regents).  All Rights Reserved.
c*
c* Permission to use, copy, modify, and distribute this software and its
c* documentation for educational, research, and not-for-profit purposes,
c* without fee and without a signed licensing agreement, is hereby granted,
c* provided that the above copyright notice, this paragraph and the
c* following two paragraphs appear in all copies, modifications, and
c* distributions.
c*
c* IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
c* SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
c* ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
c* REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
c*
c* REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
c* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
c* PARTICULAR PURPOSE.  THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF
c* ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS".  REGENTS HAS NO OBLIGATION
c* TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
c* MODIFICATIONS.
c*
c* @endverbatim
c**************************************************************************
      subroutine pcolcomp(nrow, ncol, nnzero, values, rowind, colptr, 
     2  path, title, mxtype)
c**************************************************************************
c* Routine:  pcolcomp
c* Purpose:  Print a column-compressed matrix in Harwell-Boeing format
c* Author:   Nathan Baker (mostly ripped off from Harwell-Boeing format
c*           documentation)
c**************************************************************************
      implicit none
      character title*72 , path*72  , key*8    , mxtype*3 ,
     1          ptrfmt*16, indfmt*16, valfmt*20, rhsfmt*20
      integer   totcrd, ptrcrd, indcrd, valcrd, rhscrd,
     1          nrow  , ncol  , nnzero, neltvl 
      integer   i, myunit, openstat, closestat
      integer   colptr(*), rowind(*)
      real      values(*)

c* Open the file for reading
      myunit = 12
      open (UNIT=myunit, FILE=path, STATUS="unknown")

c* Set some default values
      key = "key"
      ptrcrd = int(ncol/13 + 1) - 1
      indcrd = int(nnzero/13 + 1) - 1
      valcrd = int(nnzero/13 + 1) - 1
      totcrd = ptrcrd + indcrd + valcrd
      rhscrd = 0
      neltvl = 0
      ptrfmt = "(13I6)"
      indfmt = "(13I6)"
      valfmt = "(5E15.8)"
      rhsfmt = "(5E15.8)"
  
c* Print the header
      write (myunit, 1000) title , key   ,
     1                     totcrd, ptrcrd, indcrd, valcrd, rhscrd,
     2                     mxtype, nrow  , ncol  , nnzero, neltvl,
     3                     ptrfmt, indfmt, valfmt, rhsfmt
 1000 format ( A72, A8 / 5I14 / A3, 11X, 4I14 / 2A16, 2A20 )

c* Write the matrix structure
      write (myunit, ptrfmt) (colptr (i), i = 1, ncol+1)
      write (myunit, indfmt) (rowind (i), i = 1, nnzero)

c* Write out the values
       if (valcrd .gt. 0) then
           write (myunit, valfmt) (values (i), i = 1, nnzero)
       endif

c* Close the file
      close (UNIT=myunit)

c* Enough of this FORTRAN crap!
      return
      end
