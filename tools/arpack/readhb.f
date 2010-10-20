c**************************************************************************
c*  @file    arstuff.f
c*  @author  Nathan Baker
c*  @brief   Support routines for ARPACK
c*  @version $Id$
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

