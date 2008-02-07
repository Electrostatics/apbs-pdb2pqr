c* ///////////////////////////////////////////////////////////////////////////
c* @file    ncgd.f
c* @author  Michael Holst
c* @brief   Nonlinear (Fletcher-Reeves) variation of (Hestenes-Stiefel) CG.
c* @version $Id: ncgd.f 868 2006-04-05 18:32:05Z todd_dolinsky $
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

      subroutine ncghs(nx,ny,nz,ipc,rpc,ac,cc,fc,x,p,ap,r,zk,zkp1,tmp,
     2   itmax,iters,errtol,omega,iresid,iadjoint)
c* *********************************************************************
c* purpose:
c*
c*    solves the a(x)=b using nonlinar conjugate gradients (fletcher-reeves).
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),itmax,iters,iresid,iadjoint,nx,ny,nz
      double precision omega,errtol,rsnrm
      double precision rpc(*),ac(nx,ny,nz,*),fc(nx,ny,nz),cc(nx,ny,nz)
      double precision x(nx,ny,nz),r(nx,ny,nz),p(nx,ny,nz),ap(nx,ny,nz)
      double precision zk(nx,ny,nz),zkp1(nx,ny,nz),tmp(nx,ny,nz)
      double precision rhok1,rhok2,alpha,beta,xdot
c*
cmdir 0 0
c*
c*    *** setup for the looping ***
      iters = 0
      if ((iters .ge. itmax) .and. (iresid .eq. 0)) goto 99
      call nmresid(nx,ny,nz,ipc,rpc,ac,cc,fc,x,r,zk)
      if (iters .ge. itmax) goto 99
 30   continue
c*
c*       *** compute/check the current stopping test ***
         rhok2 = xdot(nx,ny,nz,r,r)
         rsnrm = dsqrt(rhok2)
         if (rsnrm .le. errtol) goto 99
         if (iters .ge. itmax) goto 99
c*
c*       *** form new direction vector from old one and residual ***
         if (iters .eq. 0) then
            call xcopy(nx,ny,nz,r,p)
         else
            beta = rhok2 / rhok1
            call xaxpy(nx,ny,nz,((1.0d0)/beta),r,p)
            call xscal(nx,ny,nz,beta,p)
         endif
c*
c*       *** nonlinear case: do a line search ***
c*       *** (note: "ap,zk,zkp1" passed back from line search as desired) ***
         call xcopy(nx,ny,nz,r,tmp)
         call linesearch(nx,ny,nz,alpha,
     2      ipc,rpc,ac,cc,fc,p,x,tmp,ap,zk,zkp1)
c*
c*       *** save rhok2 for next iteration ***
         rhok1 = rhok2
c*
c*       *** update solution in direction p of length alpha ***
         call xaxpy(nx,ny,nz,alpha,p,x)
c*
c*       *** update residual ***
         call xaxpy(nx,ny,nz,(-alpha),ap,r)
         call xaxpy(nx,ny,nz,(1.0d0),zk,r)
         call xaxpy(nx,ny,nz,(-1.0d0),zkp1,r)
c*
c*       *** some bookkeeping ***
         iters = iters + 1
      goto 30
c*
c*    *** return and end ***
 99   continue
      return
      end
