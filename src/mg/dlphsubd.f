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
      print*,'% DELGET: stub routine...no function...'
c*
c*    *** end it ***
      return
      end

