c* ///////////////////////////////////////////////////////////////////////////
c* @file    dlphsubd.f
c* @author  Michael Holst
c* @brief   Stub file to simulate interaction with DelPhi electrostatics.
c* @version $Id: dlphsubd.f 868 2006-04-05 18:32:05Z todd_dolinsky $
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

      subroutine delget(nx,ny,nz,
     2   xmin,xmax,ymin,ymax,zmin,zmax,
     3   epsin,epsout,rionst,temper,ncrgpt,
     4   iepsmap,idebmap,icrgpos,crg,phi)
c* **********************************************************************
c* DOUBLE PRECISION VERSION
c*
c* simple stub that doesn't do anything...
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          nx,ny,nz
      double precision xmin,xmax,ymin,ymax,zmin,zmax
      double precision epsin,epsout,rionst,temper
      integer          iepsmap(*),idebmap(*)
      integer          icrgpos(*),ncrgpt
      double precision phi(*),crg(*)
c*
c*    *** i/o ***
      call vnmprt(2,'% DELGET: stub routine...no function...',39)
c*
c*    *** end it ***
      return
      end

