c* ///////////////////////////////////////////////////////////////////////////
c* @file    nrchd.f
c* @author  Michael Holst
c* @brief   Nonlinear version of Richardson iteration.
c* @version $Id: nrchd.f 868 2006-04-05 18:32:05Z todd_dolinsky $
c* @attention
c* @verbatim
c*
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
c*
c* @endverbatim
c* ///////////////////////////////////////////////////////////////////////////

      subroutine nrich(nx,ny,nz,ipc,rpc,ac,cc,fc,x,w1,w2,r,
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
      double precision x(nx,ny,nz),w1(nx,ny,nz),w2(nx,ny,nz)
      double precision r(nx,ny,nz)
c*
cmdir 0 0
c*
c*    *** do in one step ***
      numdia = ipc(11)
      if (numdia .eq. 7) then
         call nrich7(nx,ny,nz,ipc,rpc,ac(1,1),cc,fc,
     2      ac(1,2),ac(1,3),ac(1,4),
     3      x,w1,w2,r,
     4      itmax,iters,errtol,omega,iresid,iadjoint)
      elseif (numdia .eq. 27) then
         call nrich27(nx,ny,nz,ipc,rpc,ac(1,1),cc,fc,
     2      ac(1,2),ac(1,3),ac(1,4),ac(1,5),ac(1,6),ac(1,7),ac(1,8),
     3      ac(1,9),ac(1,10),ac(1,11),ac(1,12),ac(1,13),ac(1,14),
     4      x,w1,w2,r,
     5      itmax,iters,errtol,omega,iresid,iadjoint)
      else
         call vnmprt(2,'% NRICH: invalid stencil type given...',38)
      endif
c*
c*    *** return and end ***
      return
      end
      subroutine nrich7(nx,ny,nz,ipc,rpc,
     2   oC,cc,fc,
     3   oE,oN,uC,
     4   x,w1,w2,r,
     5   itmax,iters,errtol,omega,iresid,iadjoint)
c* *********************************************************************
c* purpose: 
c*
c*   nonlinear richardson's iteration (hackbusch, page 185)
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),itmax,iters,iresid,iadjoint,nx,ny,nz
      integer          i,j,k
      double precision omega,errtol
      double precision rpc(*),oE(nx,ny,nz),oN(nx,ny,nz),uC(nx,ny,nz)
      double precision fc(nx,ny,nz),oC(nx,ny,nz),cc(nx,ny,nz)
      double precision x(nx,ny,nz),w1(nx,ny,nz),w2(nx,ny,nz)
      double precision r(nx,ny,nz)
c*
cmdir 0 0
c*
c*    *** do the richardson's iteration itmax times ***
      do 30 iters = 1, itmax
c*
c*       *** first get the residual ***
         call nmresid7_1s(nx,ny,nz,ipc,rpc,oC,cc,fc,oE,oN,uC,x,r,w1)
c*
c*       *** do it ***
cmdir 3 1
         do 10 k=2,nz-1
cmdir 3 2
            do 11 j=2,ny-1
cmdir 3 3
               do 12 i=2,nx-1
                  x(i,j,k) = x(i,j,k) + omega*r(i,j,k)
 12            continue
 11         continue
 10      continue
c*
c*       *** main loop ***
 30   continue
c*
c*    *** if specified, return the new residual as well ***
      if (iresid .eq. 1) then
         call nmresid7_1s(nx,ny,nz,ipc,rpc,oC,cc,fc,oE,oN,uC,x,r,w1)
      endif
c*
c*    *** return and end ***
      return
      end
      subroutine nrich27(nx,ny,nz,ipc,rpc,
     2   oC,cc,fc,
     3   oE,oN,uC,oNE,oNW,uE,uW,uN,uS,uNE,uNW,uSE,uSW,
     4   x,w1,w2,r,
     5   itmax,iters,errtol,omega,iresid,iadjoint)
c* *********************************************************************
c* purpose: 
c*
c*   nonlinear richardson's iteration (hackbusch, page 185)
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),itmax,iters,iresid,iadjoint,nx,ny,nz
      integer          i,j,k
      double precision omega,errtol
      double precision rpc(*),oE(nx,ny,nz),oN(nx,ny,nz),uC(nx,ny,nz)
      double precision oNE(nx,ny,nz),oNW(nx,ny,nz),uE(nx,ny,nz)
      double precision uW(nx,ny,nz),uN(nx,ny,nz),uS(nx,ny,nz)
      double precision uNE(nx,ny,nz),uNW(nx,ny,nz),uSE(nx,ny,nz)
      double precision uSW(nx,ny,nz)
      double precision fc(nx,ny,nz),oC(nx,ny,nz),cc(nx,ny,nz)
      double precision x(nx,ny,nz),w1(nx,ny,nz),w2(nx,ny,nz)
      double precision r(nx,ny,nz)
c*
cmdir 0 0
c*
c*    *** do the richardson's iteration itmax times ***
      do 30 iters = 1, itmax
c*
c*       *** first get the residual ***
         call mresid27_1s(nx,ny,nz,ipc,rpc,oC,cc,fc,
     2      oE,oN,uC,oNE,oNW,uE,uW,uN,uS,uNE,uNW,uSE,uSW,
     3      x,r)
c*
c*       *** do it ***
cmdir 3 1
         do 10 k=2,nz-1
cmdir 3 2
            do 11 j=2,ny-1
cmdir 3 3
               do 12 i=2,nx-1
                  x(i,j,k) = x(i,j,k) + omega*r(i,j,k)
 12            continue
 11         continue
 10      continue
c*
c*       *** main loop ***
 30   continue
c*
c*    *** if specified, return the new residual as well ***
      if (iresid .eq. 1) then
         call nmresid27_1s(nx,ny,nz,ipc,rpc,oC,cc,fc,
     2      oE,oN,uC,oNE,oNW,uE,uW,uN,uS,uNE,uNW,uSE,uSW,
     3      x,r,w1)
      endif
c*
c*    *** return and end ***
      return
      end

