c* ///////////////////////////////////////////////////////////////////////////
c* @file    gsd.f
c* @author  Michael Holst
c* @brief   Gauss-Seidel iteration.
c* @version $Id: gsd.f 868 2006-04-05 18:32:05Z todd_dolinsky $
c* @attention
c* @verbatim
c*
c* PMG -- Parallel algebraic MultiGrid
c* Copyright (c) 1994-2006.  Michael Holst.
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
c* PMG is free software; you can redistribute it and/or modify
c* it under the terms of the GNU General Public License as published by
c* the Free Software Foundation; either version 2 of the License, or
c* (at your option) any later version.
c*
c* PMG is distributed in the hope that it will be useful,
c* but WITHOUT ANY WARRANTY; without even the implied warranty of
c* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c* GNU General Public License for more details.
c*
c* You should have received a copy of the GNU General Public License
c* along with PMG; if not, write to the Free Software
c* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
c*
c* Linking PMG statically or dynamically with other modules is making a
c* combined work based on PMG. Thus, the terms and conditions of the GNU
c* General Public License cover the whole combination.
c* 
c* SPECIAL GPL EXCEPTION
c* In addition, as a special exception, the copyright holders of PMG
c* give you permission to combine the PMG program with free software
c* programs and libraries that are released under the GNU LGPL or with
c* code included in releases of ISIM, Ion Simulator Interface, PMV, PyMOL
c* SMOL, VMD, and Vision. Such combined software may be linked with PMG and 
c* redistributed together in original or modified form as mere aggregation
c* without requirement that the entire work be under the scope of the GNU 
c* General Public License. This special exception permission is also extended
c* to any software listed in the SPECIAL GPL EXCEPTION clauses by the FEtk and
c* APBS libraries.
c* 
c* Note that people who make modified versions of PMG are not obligated
c* to grant this special exception for their modified versions; it is
c* their choice whether to do so. The GNU General Public License gives
c* permission to release a modified version without this exception; this
c* exception also makes it possible to release a modified version which
c* carries forward this exception.
c*
c* @endverbatim
c* ///////////////////////////////////////////////////////////////////////////

      subroutine gsrb(nx,ny,nz,ipc,rpc,ac,cc,fc,x,w1,w2,r,
     2   itmax,iters,errtol,omega,iresid,iadjoint)
c* *********************************************************************
c* purpose:
c*
c*    call the fast diagonal iterative method.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),itmax,iters,iresid,iadjoint,nx,ny,nz
      integer          numdia
      double precision omega,errtol
      double precision rpc(*),ac(nx*ny*nz,*),cc(nx,ny,nz),fc(nx,ny,nz)
      double precision x(nx,ny,nz),w1(nx,ny,nz),w2(nx,ny,nz),r(nx,ny,nz)
c*
cmdir 0 0
c*
c*    *** do in one step ***
      numdia = ipc(11)
      if (numdia .eq. 7) then
         call gsrb7(nx,ny,nz,ipc,rpc,ac(1,1),cc,fc,
     2      ac(1,2),ac(1,3),ac(1,4),
     3      x,w1,w2,r,
     4      itmax,iters,errtol,omega,iresid,iadjoint)
      elseif (numdia .eq. 27) then
         call gsrb27(nx,ny,nz,ipc,rpc,ac(1,1),cc,fc,
     2      ac(1,2),ac(1,3),ac(1,4),ac(1,5),ac(1,6),ac(1,7),ac(1,8),
     3      ac(1,9),ac(1,10),ac(1,11),ac(1,12),ac(1,13),ac(1,14),
     4      x,w1,w2,r,
     4      itmax,iters,errtol,omega,iresid,iadjoint)
      else
         call vnmprt(2,'% GSRB: invalid stencil type given...',37)
      endif
c*
c*    *** return and end ***
      return
      end
      subroutine gsrb7(nx,ny,nz,ipc,rpc,
     2   oC,cc,fc,
     3   oE,oN,uC,
     4   x,w1,w2,r,
     5   itmax,iters,errtol,omega,iresid,iadjoint)
c* *********************************************************************
c* purpose: 
c*
c*    7 diagonal gauss-seidel routine.
c*
c*    this routine applies the gauss-seidel operator or its
c*    adjoint depending on the flag iadjoint.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),itmax,iters,iresid,iadjoint,nx,ny,nz
      integer          i,j,k,i1,i2,j1,j2,k1,k2,istep
      double precision omega,errtol
      double precision rpc(*),oE(nx,ny,nz),oN(nx,ny,nz),uC(nx,ny,nz)
      double precision fc(nx,ny,nz),oC(nx,ny,nz),cc(nx,ny,nz)
      double precision x(nx,ny,nz),w1(nx,ny,nz),w2(nx,ny,nz)
      double precision r(nx,ny,nz)
c*
cmdir 0 0
c*
c*    *** do the gauss-seidel iteration itmax times ***
      i1 = (1-iadjoint) * 2 + iadjoint     * (nx-1)
      i2 = iadjoint     * 2 + (1-iadjoint) * (nx-1)
      j1 = (1-iadjoint) * 2 + iadjoint     * (ny-1)
      j2 = iadjoint     * 2 + (1-iadjoint) * (ny-1)
      k1 = (1-iadjoint) * 2 + iadjoint     * (nz-1)
      k2 = iadjoint     * 2 + (1-iadjoint) * (nz-1)
      istep = iadjoint*(-1) + (1-iadjoint)*(1)
      do 30 iters = 1, itmax
c*
c*       *** do the red points ***
cmdir 3 1
         do 10 k=k1,k2,istep
cmdir 3 2
            do 11 j=j1,j2,istep
cmdir 3 3
               do 12 i=i1,i2,istep
                  x(i,j,k) = (fc(i,j,k)+(
     2               +  oN(i,j,k)        * x(i,j+1,k)
     3               +  oN(i,j-1,k)      * x(i,j-1,k)
     4               +  oE(i,j,k)        * x(i+1,j,k)
     5               +  oE(i-1,j,k)      * x(i-1,j,k)
     6               +  uC(i,j,k-1)      * x(i,j,k-1)
     7               +  uC(i,j,k)        * x(i,j,k+1)
     8               )) /  (oC(i,j,k) + cc(i,j,k))
 12            continue
 11         continue
 10      continue
c*
c*       *** main loop ***
 30   continue
c*
c*    *** if specified, return the new residual as well ***
      if (iresid .eq. 1) then
         call mresid7_1s(nx,ny,nz,ipc,rpc,oC,cc,fc,oE,oN,uC,x,r)
      endif
c*
c*    *** return and end ***
      return
      end
      subroutine gsrb27(nx,ny,nz,ipc,rpc,
     2   oC,cc,fc,
     3   oE,oN,uC,oNE,oNW,uE,uW,uN,uS,uNE,uNW,uSE,uSW,
     4   x,w1,w2,r,
     5   itmax,iters,errtol,omega,iresid,iadjoint)
c* *********************************************************************
c* purpose: 
c*
c*    27 diagonal gauss-seidel routine.
c*
c*    this routine applies the gauss-seidel operator or its
c*    adjoint depending on the flag iadjoint.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),itmax,iters,iresid,iadjoint,nx,ny,nz
      integer          i,j,k,i1,i2,j1,j2,k1,k2,istep
      double precision omega,errtol
      double precision rpc(*),oE(nx,ny,nz),oN(nx,ny,nz),uC(nx,ny,nz)
      double precision oNE(nx,ny,nz),oNW(nx,ny,nz),uE(nx,ny,nz)
      double precision uW(nx,ny,nz),uN(nx,ny,nz),uS(nx,ny,nz)
      double precision uNE(nx,ny,nz),uNW(nx,ny,nz),uSE(nx,ny,nz)
      double precision uSW(nx,ny,nz)
      double precision fc(nx,ny,nz),oC(nx,ny,nz),cc(nx,ny,nz)
      double precision x(nx,ny,nz),w1(nx,ny,nz),w2(nx,ny,nz)
      double precision r(nx,ny,nz)
      double precision tmpO,tmpU,tmpD
c*
cmdir 0 0
c*
c*    *** do the gauss-seidel iteration itmax times ***
      i1 = (1-iadjoint) * 2 + iadjoint     * (nx-1)
      i2 = iadjoint     * 2 + (1-iadjoint) * (nx-1)
      j1 = (1-iadjoint) * 2 + iadjoint     * (ny-1)
      j2 = iadjoint     * 2 + (1-iadjoint) * (ny-1)
      k1 = (1-iadjoint) * 2 + iadjoint     * (nz-1)
      k2 = iadjoint     * 2 + (1-iadjoint) * (nz-1)
      istep = iadjoint*(-1) + (1-iadjoint)*(1)
      do 30 iters = 1, itmax
c*
c*       *** do all of the points ***
cmdir 3 1
         do 10 k=k1,k2,istep
cmdir 3 2
            do 11 j=j1,j2,istep
cmdir 3 3
               do 12 i=i1,i2,istep
                  tmpO =
     2               +  oN(i,j,k)        * x(i,j+1,k)
     3               +  oN(i,j-1,k)      * x(i,j-1,k)
     4               +  oE(i,j,k)        * x(i+1,j,k)
     5               +  oE(i-1,j,k)      * x(i-1,j,k)
     6               +  oNE(i,j,k)       * x(i+1,j+1,k)
     7               +  oNW(i,j,k)       * x(i-1,j+1,k)
     8               +  oNW(i+1,j-1,k)   * x(i+1,j-1,k)
     9               +  oNE(i-1,j-1,k)   * x(i-1,j-1,k)
                  tmpU =
     2               +  uC(i,j,k)        * x(i,j,k+1)
     3               +  uN(i,j,k)        * x(i,j+1,k+1)
     4               +  uS(i,j,k)        * x(i,j-1,k+1)
     5               +  uE(i,j,k)        * x(i+1,j,k+1)
     6               +  uW(i,j,k)        * x(i-1,j,k+1)
     7               +  uNE(i,j,k)       * x(i+1,j+1,k+1)
     8               +  uNW(i,j,k)       * x(i-1,j+1,k+1)
     9               +  uSE(i,j,k)       * x(i+1,j-1,k+1)
     9               +  uSW(i,j,k)       * x(i-1,j-1,k+1)
                  tmpD =
     2               +  uC(i,j,k-1)      * x(i,j,k-1)
     3               +  uS(i,j+1,k-1)    * x(i,j+1,k-1)
     4               +  uN(i,j-1,k-1)    * x(i,j-1,k-1)
     5               +  uW(i+1,j,k-1)    * x(i+1,j,k-1)
     6               +  uE(i-1,j,k-1)    * x(i-1,j,k-1)
     7               +  uSW(i+1,j+1,k-1) * x(i+1,j+1,k-1)
     8               +  uSE(i-1,j+1,k-1) * x(i-1,j+1,k-1)
     9               +  uNW(i+1,j-1,k-1) * x(i+1,j-1,k-1)
     9               +  uNE(i-1,j-1,k-1) * x(i-1,j-1,k-1)
                  x(i,j,k) = (fc(i,j,k)+(tmpO + tmpU + tmpD
     9               )) /  (oC(i,j,k) + cc(i,j,k))
 12            continue
 11         continue
 10      continue
c*
c*       *** main loop ***
 30   continue
c*
c*    *** if specified, return the new residual as well ***
      if (iresid .eq. 1) then
         call mresid27_1s(nx,ny,nz,ipc,rpc,oC,cc,fc,
     2      oE,oN,uC,oNE,oNW,uE,uW,uN,uS,uNE,uNW,uSE,uSW,
     3      x,r)
      endif
c*
c*    *** return and end ***
      return
      end
      subroutine gsrb7X(nx,ny,nz,ipc,rpc,
     2   oC,cc,fc,
     3   oE,oN,uC,
     4   x,w1,w2,r,
     5   itmax,iters,errtol,omega,iresid,iadjoint)
c* *********************************************************************
c* purpose: 
c*
c*    fast 7 diagonal red/black gauss-seidel routine.
c*
c*    this routine applies the red/black gauss-seidel operator or its
c*    adjoint depending on the flag iadjoint.  note that the adjoint
c*    in this case is simply doing the red or black points first.
c*
c*    note that the interior grid points from 2,...,nx-1, etc.
c*    we then begin coloring with a red point, and do a point red/black
c*    coloring.
c*
c*    the red points are:  
c*
c*       if ((j even) and (k even)) then 
c*          begin row at first point, i=2 (i even), or ioff = 0
c*       else if ((j odd) and (k even)) then
c*          begin row at second point, i=3 (i odd), or ioff = 1
c*       else if ((j even) and (k odd)) then
c*          begin row at second point, i=3 (i odd), or ioff = 1
c*       else if ((j odd) and (k odd)) then
c*          begin row at first point, i=2 (i even), or ioff = 0
c*       endif
c*       then: begin row at:  i=2+ioff
c*
c*    the appropriate ioff function for the red points is then:
c*         ioff = dabs(mod(j,2)-mod(k,2))
c*
c*    the appropriate ioff function for the black points is then:
c*         ioff = 1 - dabs(mod(j,2)-mod(k,2))
c*
c*
c*    alternatively, the red points are:
c*
c*       those whose indices add up to an even number.
c*       to see this, consider that all surrounding points are
c*       only one index different, hence the sum will differ by one.
c*       thus, if a given point has an even sum, then the surrounding
c*       points will have an odd sum.
c*
c*    thus, the black points are:
c*
c*       therefore those whose indices add up to an odd number.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),itmax,iters,iresid,iadjoint,nx,ny,nz
      integer          i,j,k,ioff
      double precision omega,errtol
      double precision rpc(*),oE(nx,ny,nz),oN(nx,ny,nz),uC(nx,ny,nz)
      double precision fc(nx,ny,nz),oC(nx,ny,nz),cc(nx,ny,nz)
      double precision x(nx,ny,nz),w1(nx,ny,nz),w2(nx,ny,nz)
      double precision r(nx,ny,nz)
c*
cmdir 0 0
c*
c*    *** do the gauss-seidel iteration itmax times ***
      do 30 iters = 1, itmax
c*
c*       *** do the red points ***
cmdir 3 1
         do 10 k=2,nz-1
cmdir 3 2
            do 11 j=2,ny-1
CZZZ           ioff = mod((j+k+2),2) 
               ioff = (1-iadjoint)*mod((j+k+2),2) 
     2              + iadjoint    *(1-mod((j+k+2),2))
cmdir 3 3
               do 12 i=2+ioff,nx-1, 2
                  x(i,j,k) = (fc(i,j,k)+(
     2               +  oN(i,j,k)        * x(i,j+1,k)
     3               +  oN(i,j-1,k)      * x(i,j-1,k)
     4               +  oE(i,j,k)        * x(i+1,j,k)
     5               +  oE(i-1,j,k)      * x(i-1,j,k)
     6               +  uC(i,j,k-1)      * x(i,j,k-1)
     7               +  uC(i,j,k)        * x(i,j,k+1)
     8               )) /  (oC(i,j,k) + cc(i,j,k))
 12            continue
 11         continue
 10      continue
c*
c*       *** do the black points ***
cmdir 3 1
         do 20 k=2,nz-1
cmdir 3 2
            do 21 j=2,ny-1
CZZZ           ioff = 1-mod((j+k+2),2)
               ioff = iadjoint    *mod((j+k+2),2) 
     2              + (1-iadjoint)*(1-mod((j+k+2),2))
cmdir 3 3
               do 22 i=2+ioff,nx-1, 2
                  x(i,j,k) = (fc(i,j,k)+(
     2               +  oN(i,j,k)        * x(i,j+1,k)
     3               +  oN(i,j-1,k)      * x(i,j-1,k)
     4               +  oE(i,j,k)        * x(i+1,j,k)
     5               +  oE(i-1,j,k)      * x(i-1,j,k)
     6               +  uC(i,j,k-1)      * x(i,j,k-1)
     7               +  uC(i,j,k)        * x(i,j,k+1)
     8               )) /  (oC(i,j,k) + cc(i,j,k))
 22            continue
 21         continue
 20      continue
c*
c*       *** main loop ***
 30   continue
c*
c*    *** if specified, return the new residual as well ***
      if (iresid .eq. 1) then
         call mresid7_1s(nx,ny,nz,ipc,rpc,oC,cc,fc,oE,oN,uC,x,r)
      endif
c*
c*    *** return and end ***
      return
      end
