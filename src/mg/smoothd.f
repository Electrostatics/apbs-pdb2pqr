c* ///////////////////////////////////////////////////////////////////////////
c* MG/XMG -- Multilevel nonlinear scalar elliptic PDE solver and X interface
c* Copyright (C) 1995  Michael Holst
c*
c* This program is free software; you can redistribute it and/or modify
c* it under the terms of the GNU General Public License as published by
c* the Free Software Foundation; either version 2 of the License, or
c* (at your option) any later version.
c*
c* This program is distributed in the hope that it will be useful,
c* but WITHOUT ANY WARRANTY; without even the implied warranty of
c* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c* GNU General Public License for more details.
c*
c* You should have received a copy of the GNU General Public License
c* along with this program; if not, write to the Free Software
c* Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
c*
c* MG/XMG was developed by:
c*
c*    Michael Holst                TELE:  (858) 534-4899
c*    Department of Mathematics    FAX:   (858) 534-5273
c*    UC San Diego, AP&M 5739      EMAIL: mholst@math.ucsd.edu
c*    La Jolla, CA 92093 USA       WEB:   http://www.scicomp.ucsd.edu/~mholst
c*
c* See the file "maind.f" for more information and pointers to papers.
c*
c* IMPORTANT: If you intend to use or modify this code, make sure you 
c* understand your responsibilities under the GNU license.
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
         print*,'% SMOOTH: bad smoothing routine specified...'
      endif
c*
c*    *** return and end ***
      return
      end

