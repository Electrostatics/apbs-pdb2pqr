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

      subroutine daxpy(n,alpha,x,istep,y,jstep)
c* *********************************************************************
c* purpose:
c*
c*    parallel/vector vector saxpy.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      double precision x(*),y(*),alpha
      integer          n,istep,jstep,i,ii,nproc,ipara,ivect
      parameter        (nproc=1)
c*
cmdir 0 0
c*
c*    *** find parallel loops (ipara), remainder (ivect) ***
      ipara = n / nproc
      ivect = mod(n,nproc)
c*
c*    *** do parallel loops ***
cmdir 2 1
      do 10 ii = 1, nproc
cmdir 2 2
         do 11 i = 1+(ipara*(ii-1)), ipara*ii
            y(i) = y(i) + alpha * x(i)
 11      continue
 10   continue
c*
c*    *** do vector loops ***
cmdir 1 1
      do 20 i = ipara*nproc+1, n
         y(i) = y(i) + alpha * x(i)
 20   continue
c*
c*    *** return and end ***
      return
      end
      subroutine dcopy(n,x,istep,y,jstep)
c* *********************************************************************
c* purpose:
c*
c*    parallel/vector vector copy.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      double precision x(*),y(*)
      integer          n,istep,jstep,i,ii,nproc,ipara,ivect
      parameter        (nproc=1)
c*
cmdir 0 0
c*
c*    *** find parallel loops (ipara), remainder (ivect) ***
      ipara = n / nproc
      ivect = mod(n,nproc)
c*
c*    *** do parallel loops ***
cmdir 2 1
      do 10 ii = 1, nproc
cmdir 2 2
         do 11 i = 1+(ipara*(ii-1)), ipara*ii
            y(i) = x(i) 
 11      continue
 10   continue
c*
c*    *** do vector loops ***
cmdir 1 1
      do 20 i = ipara*nproc+1, n
         y(i) = x(i) 
 20   continue
c*
c*    *** return and end ***
      return
      end
      function dasum(n,x,istep)
c* *********************************************************************
c* purpose:
c*
c*    parallel/vector sum of vector components.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      double precision x(*)
      double precision dasum,dnrm1
      integer          n,istep
c*
cmdir 0 0
c*
c*    *** call the max-norm function ***
      dasum = dnrm1(n,x,istep)
c*
c*    *** return and end ***
      return
      end
      function dnrm1(n,x,istep)
c* *********************************************************************
c* purpose:
c*
c*    parallel/vector vector norm.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      double precision x(*)
      double precision tmp,dnrm1
      integer          n,istep,i,ii,nproc,ipara,ivect
      parameter        (nproc=1)
c*
cmdir 0 0
c*
c*    *** initialize ***
      tmp = 0.0d0
c*
c*    *** find parallel loops (ipara), remainder (ivect) ***
      ipara = n / nproc
      ivect = mod(n,nproc)
c*
c*    *** do parallel loops ***
cmdir 2 1
      do 10 ii = 1, nproc
cmdir 2 2
         do 11 i = 1+(ipara*(ii-1)), ipara*ii
            tmp = tmp + dabs(x(i))
 11      continue
 10   continue
c*
c*    *** do vector loops ***
cmdir 1 1
      do 20 i = ipara*nproc+1, n
         tmp = tmp + dabs(x(i))
 20   continue
c*
c*    *** finish up ***
      dnrm1 = tmp
c*
c*    *** return and end ***
      return
      end
      function dnrm2(n,x,istep)
c* *********************************************************************
c* purpose:
c*
c*    parallel/vector vector norm.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      double precision x(*)
      double precision tmp,dnrm2
      integer          n,istep,i,ii,nproc,ipara,ivect
      parameter        (nproc=1)
c*
cmdir 0 0
c*
c*    *** initialize ***
      tmp = 0.0d0
c*
c*    *** find parallel loops (ipara), remainder (ivect) ***
      ipara = n / nproc
      ivect = mod(n,nproc)
c*
c*    *** do parallel loops ***
cmdir 2 1
      do 10 ii = 1, nproc
cmdir 2 2
         do 11 i = 1+(ipara*(ii-1)), ipara*ii
            tmp = tmp + x(i)*x(i)
 11      continue
 10   continue
c*
c*    *** do vector loops ***
cmdir 1 1
      do 20 i = ipara*nproc+1, n
         tmp = tmp + x(i)*x(i)
 20   continue
c*
c*    *** finish up ***
      dnrm2 = dsqrt(tmp)
c*
c*    *** return and end ***
      return
      end
      function dnrm8(n,x,istep)
c* *********************************************************************
c* purpose:
c*
c*    parallel/vector vector norm.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      double precision x(*)
      double precision tmp,dnrm8
      integer          n,istep,i,ii,nproc,ipara,ivect
      parameter        (nproc=1)
c*
cmdir 0 0
c*
c*    *** initialize ***
      tmp = 0.0d0
c*
c*    *** find parallel loops (ipara), remainder (ivect) ***
      ipara = n / nproc
      ivect = mod(n,nproc)
c*
c*    *** do parallel loops ***
cmdir 2 1
      do 10 ii = 1, nproc
cmdir 2 2
         do 11 i = 1+(ipara*(ii-1)), ipara*ii
            tmp = dmax1( tmp , dabs(x(i)) )
 11      continue
 10   continue
c*
c*    *** do vector loops ***
cmdir 1 1
      do 20 i = ipara*nproc+1, n
         tmp = dmax1( tmp , dabs(x(i)) )
 20   continue
c*
c*    *** finish up ***
      dnrm8 = tmp
c*
c*    *** return and end ***
      return
      end
      subroutine dscal(n,fac,x,istep)
c* *********************************************************************
c* purpose:
c*
c*    parallel/vector vector scale.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      double precision x(*),fac
      integer          n,istep,i,ii,nproc,ipara,ivect
      parameter        (nproc=1)
c*
cmdir 0 0
c*
c*    *** find parallel loops (ipara), remainder (ivect) ***
      ipara = n / nproc
      ivect = mod(n,nproc)
c*
c*    *** do parallel loops ***
cmdir 2 1
      do 10 ii = 1, nproc
cmdir 2 2
         do 11 i = 1+(ipara*(ii-1)), ipara*ii
            x(i) = fac*x(i)
 11      continue
 10   continue
c*
c*    *** do vector loops ***
cmdir 1 1
      do 20 i = ipara*nproc+1, n
         x(i) = fac*x(i)
 20   continue
c*
c*    *** return and end ***
      return
      end
      function ddot(n,x,istep,y,jstep)
c* *********************************************************************
c* purpose:
c*
c*    parallel/vector dot product.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      double precision x(*),y(*),ddot
      integer          n,istep,jstep,i,ii,nproc,ipara,ivect
      parameter        (nproc=1)
c*
cmdir 0 0
c*
c*    *** initialize ***
      ddot = 0.0d0
c*
c*    *** find parallel loops (ipara), remainder (ivect) ***
      ipara = n / nproc
      ivect = mod(n,nproc)
c*
c*    *** do parallel loops ***
cmdir 2 1
      do 10 ii = 1, nproc
cmdir 2 2
         do 11 i = 1+(ipara*(ii-1)), ipara*ii
            ddot = ddot + x(i)*y(i)
 11      continue
 10   continue
c*
c*    *** do vector loops ***
cmdir 1 1
      do 20 i = ipara*nproc+1, n
         ddot = ddot + x(i)*y(i)
 20   continue
c*
c*    *** return and end ***
      return
      end
      function idamax(n,sx,incx)
c* *********************************************************************
c* purpose:
c*
c*    parallel/vector icamax.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          idamax,n,ns,incx,ii,i
      double precision smax,xmag
      double precision sx(*)
c*
cmdir 0 0
c*
      idamax = 0
      if(n.le.0) return
      idamax = 1
      if(n.le.1)return
      if(incx.eq.1)goto 20
c*
c*       code for increments not equal to 1.
c*
      smax = dabs(sx(1))
      ns = n*incx
      ii = 1
      do 10 i=1,ns,incx
          xmag = dabs(sx(i))
          if(xmag.le.smax) go to 40
          idamax = ii
          smax = xmag
 40       ii = ii + 1
 10   continue
      return
c*
c*    code for increments equal to 1.
c*
 20   smax = dabs(sx(1))
      do 30 i = 2,n
         xmag = dabs(sx(i))
         if(xmag.le.smax) go to 30
         idamax = i
         smax = xmag
 30   continue
      return
      end

