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
c* Nathan A. Baker (nathan.baker@pnl.gov)
c* Pacific Northwest National Laboratory
c*
c* Additional contributing authors listed in the code documentation.
c*
c* Copyright (c) 2010.  Pacific Northwest National Laboratory.
c* Portions Copyright (c) 2002-2010.  Washington University in St. Louis.
c* Portions Copyright (c) 1999-2002.  The Regents of the University 
c* of California.  
c* Portions Copyright (c) 1995.  Michael Holst.
c*
c* All rights reserved.
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
c*
c* @endverbatim
c**************************************************************************
      subroutine bcolcomp(iparm, rparm, iwork, rwork, values, rowind, 
     2  colptr, flag)
c**************************************************************************
c* Routine:  bcolcomp
c* Purpose:  Build a column-compressed matrix in Harwell-Boeing format
c* Args:     flag   0 ==> Use Poisson operator only
c*                  1 ==> Use linearization of full operator around current
c*                        solution
c* Author:   Nathan Baker (mostly ripped off from Harwell-Boeing format
c*           documentation)
c**************************************************************************
      implicit none
      integer          nrow     , ncol     , nnzero   , i        , flag
      integer          iparm(*) , iwork(*) , colptr(*), rowind(*)
      double precision rparm(*) , rwork(*) , values(*)
      integer          nxc,nyc,nzc,nf,nc,narr,narrc,n_rpc
      integer          n_iz,n_ipc,iretot,iintot
      integer          nrwk,niwk,nx,ny,nz,nlev,ierror,maxlev,mxlv
      integer          mgcoar,mgdisc,mgsolv
      integer          k_iz
      integer          k_ipc,k_rpc,k_ac,k_cc,k_fc,k_pc

c*    *** decode some parameters ***
      nrwk   = iparm(1)
      niwk   = iparm(2)
      nx     = iparm(3)
      ny     = iparm(4)
      nz     = iparm(5)
      nlev   = iparm(6)

c*    *** some checks on input ***
      mxlv = maxlev(nx,ny,nz)

c*    *** basic grid sizes, etc. ***
      mgcoar = iparm(18)
      mgdisc = iparm(19)
      mgsolv = iparm(21)
      call mgsz(mgcoar,mgdisc,mgsolv,nx,ny,nz,nlev,nxc,nyc,nzc,
     2   nf,nc,narr,narrc,n_rpc,n_iz,n_ipc,iretot,iintot)

c*    *** split up the integer work array ***
      k_iz  = 1
      k_ipc = k_iz + n_iz

c*    *** split up the real work array ***
      k_rpc = 1
      k_cc  = k_rpc + n_rpc
      k_fc  = k_cc  + narr
      k_pc  = k_fc  + narr
      k_ac  = k_pc  + 27*narrc

      call bcolcomp2(iparm, rparm, nx, ny, nz, iwork(k_iz), 
     2  iwork(k_ipc), rwork(k_rpc), rwork(k_ac), rwork(k_cc), 
     3  values, rowind, colptr, flag)


      return 
      end

      subroutine bcolcomp2(iparm, rparm, nx, ny, nz, iz,
     2  ipc, rpc, ac, cc, values, rowind, colptr, flag)
c**************************************************************************
c* Routine:  bcolcomp2
c* Purpose:  Build a column-compressed matrix in Harwell-Boeing format
c* Args:     flag   0 ==> Use Poisson operator only
c*                  1 ==> Use linearization of full operator around current
c*                        solution
c* Author:   Nathan Baker (mostly ripped off from Harwell-Boeing format
c*           documentation)
c**************************************************************************
      implicit none

      integer          iparm(*) , ipc(*)   , iz(50,*), rowind(*) 
      integer          colptr(*), flag
      double precision rparm(*) , rpc(*)   , ac(*)    , cc(*)   
      double precision values(*) 
      integer          nlev, lev, nx, ny, nz

c*    *** decode the iparm array ***
      nlev   = iparm(6)

c*    *** build the multigrid data structure in iz ***
c*    THIS MAY HAVE BEEN DONE ALREADY, BUT IT'S OK TO DO IT AGAIN,
c*    RIGHT?
c*    call buildstr (nx,ny,nz,nlev,iz)

c*    We're interested in the finest level
      lev = 1

      call bcolcomp3(nx, ny, nz, 
     1  ipc(iz(5,lev)), rpc(iz(6,lev)), ac(iz(7,lev)),
     2  cc(iz(1,lev)), values, rowind, colptr, flag)

      return 
      end

      subroutine bcolcomp3(nx, ny, nz, 
     2  ipc, rpc, ac, cc, values, rowind, colptr, flag)
c**************************************************************************
c* Routine:  bcolcomp3
c* Purpose:  Build a column-compressed matrix in Harwell-Boeing format
c* Args:     flag   0 ==> Use Poisson operator only
c*                  1 ==> Use linearization of full operator around current
c*                        solution
c* Author:   Nathan Baker (mostly ripped off from Harwell-Boeing format
c*           documentation)
c**************************************************************************
      implicit none

      integer          ipc(*), rowind(*), colptr(*),
     1                 flag, nx, ny, nz
      double precision rpc(*), ac(nx*ny*nz,*), cc(*)
      double precision values(*)

      call bcolcomp4(nx, ny, nz, ipc, rpc, ac(1,1), cc, 
     1  ac(1,2), ac(1,3), ac(1,4), values, rowind, colptr, flag)

      return
      end

      subroutine bcolcomp4(nx, ny, nz, ipc, rpc, oC, cc, oE, oN,
     1  uC, values, rowind, colptr, flag)
c**************************************************************************
c* Routine:  bcolcomp4
c* Purpose:  Build a column-compressed matrix in Harwell-Boeing format
c* Args:     flag   0 ==> Use Poisson operator only
c*                  1 ==> Use linearization of full operator around current
c*                        solution
c* Author:   Nathan Baker (mostly ripped off from Harwell-Boeing format
c*           documentation)
c**************************************************************************
      implicit none

      integer          ipc(*), rowind(*), colptr(*),
     1                 flag, nx, ny, nz, nxm2, nym2, nzm2, ii,jj,kk,ll,
     2                 i, j, k, l, inonz, iirow, nn, nrow, ncol, nonz,
     3                 irow, n
      double precision oE(nx,ny,nz), oN(nx,ny,nz), uC(nx,ny,nz),
     1                 cc(nx,ny,nz), oC(nx,ny,nz), rpc(*), values(*)
      logical          doit
    


c*    Get some column, row, and nonzero information
      n = nx*ny*nz
      nxm2 = nx - 2
      nym2 = ny - 2
      nzm2 = nz - 2
      nn   = nxm2*nym2*nzm2
      ncol = nn
      nrow = nn
      nonz = 7*nn - 2*nxm2*nym2 - 2*nxm2 - 2

c*    Intialize some pointers
      inonz = 1

c*    Run over the dimensions of the matrix (non-zero only in the interior 
c*    of the mesh
      do 10 k=2, nz-1
c*      Offset the index to the output grid index
        kk = k - 1
        do 11 j=2, ny-1
c*        Offset the index to the output grid index
          jj = j - 1
          do 12 i=2, nx-1
c*          Offset the index to the output grid index
            ii = i - 1

c*          Get the output (i,j,k) row number in natural ordering
            ll = (kk-1)*nxm2*nym2 + (jj-1)*nxm2 + (ii-1) + 1
             l = (k -1)*nx  *ny   + (j -1)*nx   + (i -1) + 1 

c*          Store where this column starts
            colptr(ll) = inonz

c*          SUB-DIAGONAL 3
            iirow = ll - nxm2 * nym2
            irow = l - nx * ny
            doit = (iirow.ge.1).and.(iirow.le.nn)
            doit = doit.and.(irow.ge.1).and.(irow.le.n)
            if (doit) then
              values(inonz) = -uC(i,j,k-1)
c*            if (uC(i,j,k-1).eq.0) then
c*                write (*,*) 'SUB3 uC(', i, ' ', j, ' ', k-1, ') ZERO!'
c*            endif
              rowind(inonz) = iirow
              inonz = inonz + 1
            endif

c*          SUB-DIAGONAL 2 
            iirow = ll - nxm2
            irow = l - nx
            doit = (iirow.ge.1).and.(iirow.le.nn)
            doit = doit.and.(irow.ge.1).and.(irow.le.n)
            if (doit) then
              values(inonz) = -oN(i,j-1,k)
c*            if (oN(i,j-1,k).eq.0) then
c*                write (*,*) 'SUB2 oN(', i, ' ', j-1, ' ', k, ') ZERO!'
c*            endif
              rowind(inonz) = iirow
              inonz = inonz + 1
            endif 

c*          SUB-DIAGONAL 1 
            iirow = ll - 1
            irow = l - 1
            doit = (iirow.ge.1).and.(iirow.le.nn)
            doit = doit.and.(irow.ge.1).and.(irow.le.n)
            if (doit) then
              values(inonz) = -oE(i-1,j,k)
c*            if (oE(i-1,j,k).eq.0) then
c*                write (*,*) 'SUB1 oE(', i-1, ' ', j, ' ', k, ') ZERO!'
c*            endif
              rowind(inonz) = iirow
              inonz = inonz + 1
            endif 

c*          DIAGONAL 
            iirow = ll
            irow = l
c*          if (oC(i,j,k).eq.0) then
c*              write (*,*) 'oC(', i, ' ', j, ' ', k, ') ZERO!'
c*          endif
            if (flag.eq.0) then
              values(inonz) = oC(i,j,k)
            elseif (flag.eq.1) then
              values(inonz) = oC(i,j,k) + cc(i,j,k)
            else
              stop 'PMGF1'
            endif
            rowind(inonz) = iirow
            inonz = inonz + 1
  
c*          SUPER-DIAGONAL 1
            iirow = ll + 1
            irow = l + 1
            doit = (iirow.ge.1).and.(iirow.le.nn)
            doit = doit.and.(irow.ge.1).and.(irow.le.n)
            if (doit) then
              values(inonz) = -oE(i,j,k)
c*            if (oE(i,j,k).eq.0) then
c*                write (*,*) 'SUP1 oE(', i, ' ', j, ' ', k, ') ZERO!'
c*            endif
              rowind(inonz) = iirow
              inonz = inonz + 1
            endif

c*          SUPER-DIAGONAL 2
            iirow = ll + nxm2
            irow = l + nx
            doit = (iirow.ge.1).and.(iirow.le.nn)
            doit = doit.and.(irow.ge.1).and.(irow.le.n)
            if (doit) then
              values(inonz) = -oN(i,j,k)
c*            if (oN(i,j,k).eq.0) then
c*                write (*,*) 'SUP2 oN(', i, ' ', j, ' ', k, ') ZERO!'
c*            endif
              rowind(inonz) = iirow
              inonz = inonz + 1
            endif 
 
c*          SUPER-DIAGONAL 3
            iirow = ll + nxm2*nym2
            irow = l + nx*ny
            doit = (iirow.ge.1).and.(iirow.le.nn)
            doit = doit.and.(irow.ge.1).and.(irow.le.n)
            if (doit) then
              values(inonz) = -uC(i,j,k)
c*            if (uC(i,j,k).eq.0) then
c*                write (*,*) 'SUP3 uC(', i, ' ', j, ' ', k, ') ZERO!'
c*            endif
              rowind(inonz) = iirow
              inonz = inonz + 1
            endif

 12       continue
 11     continue
 10   continue

      colptr(ncol+1) = inonz

      if (inonz.ne.(nonz+1)) then
        write (*,*) 'BCOLCOMP4:  ERROR -- INONZ = ', inonz
        write (*,*) 'BCOLCOMP4:  ERROR -- NONZ = ', nonz
        stop 'PMGF2'
      endif 

      return 
      end


      subroutine pcolcomp(nrow, ncol, nnzero, values, rowind, colptr, 
     2  path, title, mxtype)
c**************************************************************************
c* Routine:  pcolcomp
c* Purpose:  Print a column-compressed matrix in Harwell-Boeing format
c* Author:   Nathan Baker (mostly ripped off from Harwell-Boeing format
c*           documentation)
c**************************************************************************
      implicit none
      character        title*72 , path*72  , key*8    , mxtype*3 ,
     1                 ptrfmt*16, indfmt*16, valfmt*20, rhsfmt*20
      integer          totcrd   , ptrcrd   , indcrd  , valcrd,
     1                 nrow     , ncol     , nnzero  , neltvl,
     2                 i        , myunit   , openstat, closestat,
     3                 colptr(ncol+1), rowind(nnzero), rhscrd
      double precision values(nnzero)

c* Open the file for reading
      myunit = 12
      open (UNIT=myunit, FILE=path, STATUS="unknown")

c* Set some default values
      key = "key"
      ptrcrd = int(ncol/10 + 1) - 1
      indcrd = int(nnzero/10 + 1) - 1
      valcrd = int(nnzero/10 + 1) - 1
      totcrd = ptrcrd + indcrd + valcrd
      rhscrd = 0
      neltvl = 0
      ptrfmt = "(10I8)"
      indfmt = "(10I8)"
      valfmt = "(5E15.8)"
      rhsfmt = "(5E15.8)"
  
c* Print the header
      write (myunit, 1000) title , key   ,
     1                     totcrd, ptrcrd, indcrd, valcrd, rhscrd,
     2                     mxtype, nrow  , ncol  , nnzero, neltvl,
     3                     ptrfmt, indfmt, valfmt, rhsfmt
 1000 format ( A72, A8 / 5I14 / A3, 11X, 4I14 / 2A16, 2A20 )

c* Write the matrix structure
      write (myunit, ptrfmt) (colptr(i), i = 1, ncol+1)
      write (myunit, indfmt) (rowind(i), i = 1, nnzero)

c* Write out the values
      if (valcrd .gt. 0) then
          write (myunit, valfmt) (values(i), i = 1, nnzero)
      endif
 
c  Do it again in a loop 
c      do 10 i=1, nnzero
c          if (values(i).ne.0.0) then 
c            write (*,*) 'values(', i, ') = ', values(i)
c          endif
c 10   continue


c* Close the file
      close (UNIT=myunit)

c* Enough of this FORTRAN crap!
      return
      end

