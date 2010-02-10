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
c* Nathan A. Baker (baker@biochem.wustl.edu)
c* Dept. of Biochemistry and Molecular Biophysics
c* Center for Computational Biology
c* Washington University in St. Louis
c*
c* Additional contributing authors listed in the code documentation.
c*
c* Copyright (c) 2002-2010, Washington University in St. Louis.
c* Portions Copyright (c) 2002-2010.  Nathan A. Baker
c* Portions Copyright (c) 1999-2002.  The Regents of the University of California.
c* Portions Copyright (c) 1995.  Michael Holst
c*
c* All rights reserved.
c* 
c* Redistribution and use in source and binary forms, with or without
c* modification, are permitted provided that the following conditions are met: 
c*
c* -  Redistributions of source code must retain the above copyright notice, this
c* list of conditions and the following disclaimer.  
c* 
c* - Redistributions in binary form must reproduce the above copyright notice,
c* this list of conditions and the following disclaimer in the documentation
c* and/or other materials provided with the distribution.
c* 
c* - Neither the name of Washington University in St. Louis nor the names of its
c* contributors may be used to endorse or promote products derived from this
c* software without specific prior written permission.
c* 
c* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
c* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
c* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
c* A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
c* CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
c* EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
c* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
c* PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
c* LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
c* NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
c* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

