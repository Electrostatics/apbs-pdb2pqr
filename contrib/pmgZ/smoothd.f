c* ///////////////////////////////////////////////////////////////////////////
c* @file    smoothd.f
c* @author  Michael Holst
c* @brief   Generic linear smoothing iteration.
c* @version $Id: smoothd.f 868 2006-04-05 18:32:05Z todd_dolinsky $
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
      subroutine smooth(nx,ny,nz,ipc,rpc,ac,cc,fc,x,w1,w2,r,
     2   itmax,iters,errtol,omega,iresid,iadjoint,meth)
c* *********************************************************************
c* purpose: 
c*
c*    call the appropriate linear smoothing routine.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),itmax,iters,iresid,iadjoint,nx,ny,nz,meth
      double precision omega,errtol
      double precision rpc(*),ac(nx*ny*nz,*),cc(nx,ny,nz),fc(nx,ny,nz)
      double precision x(nx,ny,nz),w1(nx,ny,nz),w2(nx,ny,nz)
      double precision r(nx,ny,nz)
c*
cmdir 0 0
c*
c*    *** do in one step ***
      if (meth .eq. 0) then
         call wjac(nx,ny,nz,ipc,rpc,ac,cc,fc,x,w1,w2,r,
     2      itmax,iters,errtol,omega,iresid,iadjoint)
      elseif (meth .eq. 1) then
         call gsrb(nx,ny,nz,ipc,rpc,ac,cc,fc,x,w1,w2,r,
     2      itmax,iters,errtol,omega,iresid,iadjoint)
      elseif (meth .eq. 2) then
         call sor(nx,ny,nz,ipc,rpc,ac,cc,fc,x,w1,w2,r,
     2      itmax,iters,errtol,omega,iresid,iadjoint)
      elseif (meth .eq. 3) then
         call rich(nx,ny,nz,ipc,rpc,ac,cc,fc,x,w1,w2,r,
     2      itmax,iters,errtol,omega,iresid,iadjoint)
      elseif (meth .eq. 4) then
         call cghs(nx,ny,nz,ipc,rpc,ac,cc,fc,x,w1,w2,r,
     2      itmax,iters,errtol,omega,iresid,iadjoint)
      else
         call vnmpri(2, '% SMOOTH: bad smoothing routine specified = ', 
     2      42, meth)
      endif
c*
c*    *** return and end ***
      return
      end

