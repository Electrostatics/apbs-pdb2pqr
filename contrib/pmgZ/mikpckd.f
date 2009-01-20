c* ///////////////////////////////////////////////////////////////////////////
c* @file    mikpckd.f
c* @author  Michael Holst
c* @brief   A collection of useful low-level routines (timing, etc).
c* @version $Id: mikpckd.f 868 2006-04-05 18:32:05Z todd_dolinsky $
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

      function epsmac (idum)
c* *********************************************************************
c* purpose:
c*
c*    this routine computes the unit roundoff of the machine in single
c*    precision.  this is defined as the smallest positive machine
c*    number u such that  1.0d0 + u .ne. 1.0d0 (in single precision).
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          idum
      double precision epsmac
c* ***double precision u, comp
c*
c* ****** determine the machine epsilon ***
c* ***u = 1.0d0
c* ***10 continue
c* ***   u = u*0.50d0
c* ***   comp = 1.0d0 + u
c* ***if (comp .ne. 1.0d0) goto 10
c* ***epsmac = u*2.0d0
c*
c*    *** hardcoded ***
      epsmac = 1.0d-9
c*
c*    *** return and end ***
      return
      end
      subroutine xaxpy(nx,ny,nz,alpha,x,y)
c* *********************************************************************
c* purpose:
c*
c*    saxpy operation for a grid function with boundary values.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          nx,ny,nz,i,j,k
      double precision x(nx,ny,nz),y(nx,ny,nz),alpha
c*
cmdir 0 0
c*
c*    *** do it ***
cmdir 3 1
      do 10 k = 2, nz-1
cmdir 3 2
         do 11 j = 2, ny-1
cmdir 3 3
            do 12 i = 2, nx-1
               y(i,j,k) = y(i,j,k) + alpha * x(i,j,k)
 12         continue
 11      continue
 10   continue
c*
c*    *** return and end ***
      return
      end
      subroutine xcopy(nx,ny,nz,x,y)
c* *********************************************************************
c* purpose:
c*
c*    copy operation for a grid function with boundary values.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          nx,ny,nz,i,j,k
      double precision x(nx,ny,nz),y(nx,ny,nz)
c*
cmdir 0 0
c*
c*    *** do it ***
!$OMP PARALLEL DO private(i,j,k)
      do k = 2, nz-1
         do j = 2, ny-1
            do i = 2, nx-1
               y(i,j,k) = x(i,j,k)
			end do
		 end do
	  end do
!$OMP END PARALLEL DO
c*
c*    *** return and end ***
      return
      end

      function xnrm1(nx,ny,nz,x)
c* *********************************************************************
c* purpose:
c*
c*    norm operation for a grid function with boundary values.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          nx,ny,nz,i,j,k
      double precision x(nx,ny,nz)
      double precision tmp,xnrm1
c*
cmdir 0 0
c*
c*    *** initialize ***
      tmp = 0.0d0
c*
c*    *** do it ***
cmdir 3 1
      do 10 k = 2, nz-1
cmdir 3 2
         do 11 j = 2, ny-1
cmdir 3 3
            do 12 i = 2, nx-1
               tmp = tmp + dabs(x(i,j,k))
 12         continue
 11      continue
 10   continue
c*
c*    *** finish up ***
      xnrm1 = tmp
c*
c*    *** return and end ***
      return
      end
      function xnrm2(nx,ny,nz,x)
c* *********************************************************************
c* purpose:
c*
c*    norm operation for a grid function with boundary values.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          nx,ny,nz,i,j,k
      double precision x(nx,ny,nz)
      double precision tmp,xnrm2
c*
cmdir 0 0
c*
c*    *** initialize ***
      tmp = 0.0d0
c*
c*    *** do it ***
cmdir 3 1
      do 10 k = 2, nz-1
cmdir 3 2
         do 11 j = 2, ny-1
cmdir 3 3
            do 12 i = 2, nx-1
               tmp = tmp + x(i,j,k)*x(i,j,k) 
 12         continue
 11      continue
 10   continue
c*
c*    *** finish up ***
      xnrm2 = dsqrt(tmp)
c*
c*    *** return and end ***
      return
      end
      function xnrm8(nx,ny,nz,x)
c* *********************************************************************
c* purpose:
c*
c*    norm operation for a grid function with boundary values.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          nx,ny,nz,i,j,k
      double precision x(nx,ny,nz)
      double precision tmp,xnrm8
c*
cmdir 0 0
c*
c*    *** initialize ***
      tmp = 0.0d0
c*
c*    *** do it ***
cmdir 3 1
      do 10 k = 2, nz-1
cmdir 3 2
         do 11 j = 2, ny-1
cmdir 3 3
            do 12 i = 2, nx-1
               tmp = dmax1( tmp , dabs(x(i,j,k)) )
 12         continue
 11      continue
 10   continue
c*
c*    *** finish up ***
      xnrm8 = tmp
c*
c*    *** return and end ***
      return
      end
      subroutine xscal(nx,ny,nz,fac,x)
c* *********************************************************************
c* purpose:
c*
c*    scale operation for a grid function with boundary values.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          nx,ny,nz,i,j,k
      double precision x(nx,ny,nz),fac
c*
cmdir 0 0
c*
c*    *** do it ***
cmdir 3 1
      do 10 k = 2, nz-1
cmdir 3 2
         do 11 j = 2, ny-1
cmdir 3 3
            do 12 i = 2, nx-1
               x(i,j,k) = fac*x(i,j,k)
 12         continue
 11      continue
 10   continue
c*
c*    *** return and end ***
      return
      end
      function xdot(nx,ny,nz,x,y)
c* *********************************************************************
c* purpose:
c*
c*    inner product operation for a grid function with boundary values.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          nx,ny,nz,i,j,k
      double precision x(nx,ny,nz),y(nx,ny,nz),xdot
c*
cmdir 0 0
c*
c*    *** initialize ***
      xdot = 0.0d0
c*
c*    *** do it ***
cmdir 3 1
      do 10 k = 2, nz-1
cmdir 3 2
         do 11 j = 2, ny-1
cmdir 3 3
            do 12 i = 2, nx-1
               xdot = xdot + x(i,j,k)*y(i,j,k)
 12         continue
 11      continue
 10   continue
c*
c*    *** return and end ***
      return
      end
      function xdot3(nx,ny,nz,x,y)
c* *********************************************************************
c* purpose:
c*
c*    necessary mess.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          nx,ny,nz,i,j,k
      double precision x(nx,ny,nz),y(nx,ny,nz),xdot3
c*
cmdir 0 0
c*
c*    *** initialize ***
      xdot3 = 0.0d0
c*
c*    *** do it ***
cmdir 3 1
      do 10 k = 2, nz-1
cmdir 3 2
         do 11 j = 2, ny-1
cmdir 3 3
            do 12 i = 2, nx-1
               xdot3 = xdot3 + x(i,j,k) * y(i,j,k) * y(i,j,k)
 12         continue
 11      continue
 10   continue
c*
c*    *** return and end ***
      return
      end
      subroutine zeros(nx,ny,nz,x)
c* *********************************************************************
c* purpose:
c*
c*    zero out operation for a grid function with boundary values.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          nx,ny,nz,i,j,k
      double precision x(nx,ny,nz)
c*
cmdir 0 0
c*
c*    *** do it ***
cmdir 3 1
      do 10 k = 2, nz-1
cmdir 3 2
         do 11 j = 2, ny-1
cmdir 3 3
            do 12 i = 2, nx-1
               x(i,j,k) = 0.0d0
 12         continue
 11      continue
 10   continue
c*
c*    *** return and end ***
      return
      end
      subroutine xrand(nx,ny,nz,x)
c* *********************************************************************
c* purpose:
c*
c*    fill grid function with random values.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          nx,ny,nz,i,j,k,iflag
      double precision x(nx,ny,nz)
      double precision xdum
      integer vrnd
c*
cmdir 0 0
c*
c*    *** do it ***
      iflag = 1
      xdum  = dble(vrnd())
cmdir 3 1
      do 10 k = 2, nz-1
cmdir 3 2
         do 11 j = 2, ny-1
cmdir 3 3
            do 12 i = 2, nx-1
               x(i,j,k) = dble(vrnd())
 12         continue
 11      continue
 10   continue
c*
c*    *** return and end ***
      return
      end
      subroutine cinit(nx,ny,nz,x,value)
c* *********************************************************************
c* purpose:
c*
c*    initialize a grid function to a specific constant, without
c*    changing the boundary values.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          nx,ny,nz,i,j,k
      double precision x(nx,ny,nz),value
c*
cmdir 0 0
c*
c*    *** do it ***
cmdir 3 1
      do 10 k = 2, nz-1
cmdir 3 2
         do 11 j = 2, ny-1
cmdir 3 3
            do 12 i = 2, nx-1
               x(i,j,k) = value
 12         continue
 11      continue
 10   continue
c*
c*    *** return and end ***
      return
      end
      subroutine azeros(nx,ny,nz,x)
c* *********************************************************************
c* purpose:
c*
c*    zero out operation for a grid function, including boundary values.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          nx,ny,nz
      double precision x(*)
      integer          n,i,ii,nproc,ipara,ivect
      parameter        (nproc=1)
c*
cmdir 0 0
c*
c*    *** find parallel loops (ipara), remainder (ivect) ***
      n     = nx * ny * nz

!$OMP PARALLEL DO private(i)
      do i = 1, n
       x(i) = 0.0d0
	  end do
!$OMP END PARALLEL DO

c*
c*    *** return and end ***
      return
      end
      subroutine axrand(nx,ny,nz,x)
c* *********************************************************************
c* purpose:
c*
c*    fill grid function with random values, including boundary values.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          nx,ny,nz,iflag
      double precision x(*)
      double precision xdum
      integer          n,i,ii,nproc,ipara,ivect
      parameter        (nproc=1)
      integer vrnd
c*
cmdir 0 0
c*
c*    *** find parallel loops (ipara), remainder (ivect) ***
      n     = nx * ny * nz
      ipara = n / nproc
      ivect = mod(n,nproc)
      iflag = 1
      xdum  = dble(vrnd())
c*
c*    *** do parallel loops ***
cmdir 2 1
      do 10 ii = 1, nproc
cmdir 2 2
         do 11 i = 1+(ipara*(ii-1)), ipara*ii
            x(i) = dble(vrnd())
 11      continue
 10   continue
c*
c*    *** do vector loops ***
cmdir 1 1
      do 20 i = ipara*nproc+1, n
         x(i) = dble(vrnd())
 20   continue
c*
c*    *** return and end ***
      return
      end
      subroutine xcopy_small(nx,ny,nz,x,y)
c* *********************************************************************
c* purpose:
c*
c*    copy operation for a grid function with boundary values.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          nx,ny,nz,i,j,k
      double precision x(nx,ny,nz),y(nx-2,ny-2,nz-2)
c*
cmdir 0 0
c*
c*    *** do it ***
cmdir 3 1
      do 10 k = 2, nz-1
cmdir 3 2
         do 11 j = 2, ny-1
cmdir 3 3
            do 12 i = 2, nx-1
               y(i-1,j-1,k-1) = x(i,j,k)
 12         continue
 11      continue
 10   continue
c*
c*    *** return and end ***
      return
      end
      subroutine xcopy_large(nx,ny,nz,x,y)
c* *********************************************************************
c* purpose:
c*
c*    copy operation for a grid function with boundary values.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          nx,ny,nz,i,j,k
      double precision x(nx-2,ny-2,nz-2),y(nx,ny,nz)
c*
cmdir 0 0
c*
c*    *** do it ***
cmdir 3 1
      do 10 k = 2, nz-1
cmdir 3 2
         do 11 j = 2, ny-1
cmdir 3 3
            do 12 i = 2, nx-1
               y(i,j,k) = x(i-1,j-1,k-1)
 12         continue
 11      continue
 10   continue
c*
c*    *** return and end ***
      return
      end
      subroutine fboundPMG(ibound,nx,ny,nz,x,gxc,gyc,gzc)
c* *********************************************************************
c* purpose:
c*
c*    initialize a grid function to have a certain boundary value,
c*    given in the arrays gxc,gyc,gzc
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ibound,nx,ny,nz,i,j,k
      double precision x(nx,ny,nz)
      double precision gxc(ny,nz,2),gyc(nx,nz,2),gzc(nx,ny,2)
c*
cmdir 0 0
c*
c*    *** zero dirichlet? ***
      if (ibound .eq. 0) then 
         call fboundPMG00(nx,ny,nz,x)
c*
c*    *** nonzero dirichlet ***
      else
c*       *** the (i=1) and (i=nx) boundaries ***
cmdir 2 1
         do 10 k = 1, nz
cmdir 2 2
            do 11 j = 1, ny
               x(1 ,j,k) = gxc(j,k,1)
               x(nx,j,k) = gxc(j,k,2)
 11         continue
 10      continue
c*
c*       *** the (j=1) and (j=ny) boundaries ***
cmdir 2 1
         do 20 k = 1, nz
cmdir 2 2
            do 21 i = 1, nx
               x(i,1 ,k) = gyc(i,k,1)
               x(i,ny,k) = gyc(i,k,2)
 21         continue
 20      continue
c*
c*       *** the (k=1) and (k=nz) boundaries ***
cmdir 2 1
         do 30 j = 1, ny
cmdir 2 2
            do 31 i = 1, nx
               x(i,j,1 ) = gzc(i,j,1)
               x(i,j,nz) = gzc(i,j,2)
 31         continue
 30      continue
      endif
c*
c*    *** return and end ***
      return
      end
      subroutine fboundPMG00(nx,ny,nz,x)
c* *********************************************************************
c* purpose:
c*
c*    initialize a grid function to have a zero boundary value,
c*    given in the arrays gxc,gyc,gzc
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          nx,ny,nz,i,j,k
      double precision x(nx,ny,nz)
c*
cmdir 0 0
c*
c*    *** the (i=1) and (i=nx) boundaries ***
cmdir 2 1
      do 10 k = 1, nz
cmdir 2 2
         do 11 j = 1, ny
            x(1 ,j,k) = 0.0d0
            x(nx,j,k) = 0.0d0
 11      continue
 10   continue
c*
c*    *** the (j=1) and (j=ny) boundaries ***
cmdir 2 1
      do 20 k = 1, nz
cmdir 2 2
         do 21 i = 1, nx
            x(i,1 ,k) = 0.0d0
            x(i,ny,k) = 0.0d0
 21      continue
 20   continue
c*
c*    *** the (k=1) and (k=nz) boundaries ***
cmdir 2 1
      do 30 j = 1, ny
cmdir 2 2
         do 31 i = 1, nx
            x(i,j,1 ) = 0.0d0
            x(i,j,nz) = 0.0d0
 31      continue
 30   continue
c*
c*    *** return and end ***
      return
      end
      subroutine xprint(title,nx,ny,nz,x)
c* *********************************************************************
c* purpose:
c*
c*    this routine prints a grid function.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          nx,ny,nz,i,j,k
      double precision x(nx,ny,nz)
      character*10     title
c*
c*    *** doit ***
      call vnmprt(2,title,3)
      do 10 k = 2, nz-1
         do 11 j = 2, ny-1
              write(6,100) (x(i,j,k),i=2,nx-1)
 100          format(100(1x,1pe8.1))
 11      continue
 10   continue
c*
c*    *** return and end ***
      return
      end
      subroutine prtmatd(nx,ny,nz,ipc,rpc,ac)
c* *********************************************************************
c* purpose:
c*
c*    interface to routine to print sparse diagonal-format matrix.
c*
c* note:
c*
c*    the matrix is symmetric, therefore we print only the main and
c*    upper diagonals.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),nx,ny,nz,numdia
      double precision rpc(*),ac(nx*ny*nz,*)
c*
c*    *** do the printing ***
      numdia = ipc(11)
      if (numdia .eq. 7) then
         call prtmatd7(nx,ny,nz,ipc,rpc,
     2      ac(1,1),ac(1,2),ac(1,3),ac(1,4))
      elseif (numdia .eq. 27) then
         call prtmatd27(nx,ny,nz,ipc,rpc,
     2      ac(1,1),ac(1,2),ac(1,3),ac(1,4),
     3      ac(1,5),ac(1,6),ac(1,7),ac(1,8),ac(1,9),
     4      ac(1,10),ac(1,11),ac(1,12),ac(1,13),ac(1,14))
      else
         call vnmprt(2,'% PRTMATD: invalid stencil type given...',40)
      endif
c*
c*    *** go home ***
      return
      end
      subroutine prtmatd7(nx,ny,nz,ipc,rpc,oC,oE,oN,uC)
c* *********************************************************************
c* purpose:
c*
c*    routine to print sparse 7-diagonal-format matrix.
c*
c* note:
c*
c*    the matrix is symmetric, therefore we print only the main and
c*    upper diagonals.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),nx,ny,nz,n,i,j,k
      double precision rpc(*),oC(nx,ny,nz)
      double precision oE(nx,ny,nz),oN(nx,ny,nz),uC(nx,ny,nz)
c*
c*    *** recover matrix dimension ***
      n = (nx-2) * (ny-2) * (nz-2)
      call vnmprt(2,'% ',2)
      call vnmpri(2,'% dimension of matrix = ',24,n)
      call vnmprt(2,'% *begin diagonal matrix* ',26)
c*
c*    *** handle first block ***
      do 10 k = 2, nz-1
         do 11 j = 2, ny-1
            do 12 i = 2, nx-1
               write(*,104) oC(i,j,k),oE(i,j,k),oN(i,j,k),uC(i,j,k)
 12         continue
 11      continue
 10   continue
c*
c*    *** finish up ***
      call vnmprt(2,'% *end diagonal matrix* ',24)
c*
c*    *** format statements ***
 104  format (4(1pe8.1,1x))
c*
c*    *** return and end ***
      return
      end
      subroutine prtmatd27(nx,ny,nz,ipc,rpc,oC,oE,oN,uC,
     2   oNE,oNW,uE,uW,uN,uS,uNE,uNW,uSE,uSW)
c* *********************************************************************
c* purpose:
c*
c*    routine to print sparse 27-diagonal-format matrix.
c*
c* note:
c*
c*    the matrix is symmetric, therefore we print only the main and
c*    upper diagonals.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),nx,ny,nz,n,i,j,k
      double precision rpc(*),oC(nx,ny,nz)
      double precision oE(nx,ny,nz),oN(nx,ny,nz),uC(nx,ny,nz)
      double precision oNE(nx,ny,nz),oNW(nx,ny,nz),uE(nx,ny,nz)
      double precision uW(nx,ny,nz),uN(nx,ny,nz),uS(nx,ny,nz)
      double precision uNE(nx,ny,nz),uNW(nx,ny,nz),uSE(nx,ny,nz)
      double precision uSW(nx,ny,nz)
c*
c*    *** recover matrix dimension ***
      n = (nx-2) * (ny-2) * (nz-2)
      call vnmprt(2,'% ',2)
      call vnmpri(2,'% dimension of matrix = ',24,n)
      call vnmprt(2,'% *begin diagonal matrix* ',26)
c*
c*    *** handle first block ***
      do 10 k = 2, nz-1
         do 11 j = 2, ny-1
            do 12 i = 2, nx-1
               write(*,114) oC(i,j,k),oE(i,j,k),
     2            oNW(i,j,k),oN(i,j,k),oNE(i,j,k),
     3            uSW(i,j,k),uS(i,j,k),uSE(i,j,k),
     4            uW(i,j,k),uC(i,j,k),uE(i,j,k),
     5            uNW(i,j,k),uN(i,j,k),uNE(i,j,k)
 12         continue
 11      continue
 10   continue
c*
c*    *** finish up ***
      call vnmprt(2,'% *end diagonal matrix* ',24)
c*
c*    *** format statements ***
 114  format (14(1pe8.1,1x))
c*
c*    *** return and end ***
      return
      end
      subroutine prtmatb(a,n,m,lda)
c***********************************************************************
c* purpose:
c*
c*    routine to print sparse banded-format matrix.
c*
c* note:
c*
c*    the matrix is symmetric, therefore we print only the main and
c*    upper diagonal bands.
c*
c***********************************************************************
      implicit         none
      integer          n,m,lda,i,j
      double precision a(lda,*)
c
c     *** recover matrix dimension ***
      call vnmprt(2,' ',1)
      call vnmpri(2,' dimension of matrix = ',23,n)
      call vnmprt(2,' *begin banded matrix* ',23)
c
c     *** do the printing ***
      do 10 j = 1, n
         write(*,100) (a(i,j),i=lda,1,-1)
 10   continue
c
c     *** finish up ***
      call vnmprt(2,' *end banded matrix* ',21)
c
c     *** format statements ***
 100  format (30(1x,1pe8.1))
c
c     *** go home ***
      return
      end
      subroutine linesearch(nx,ny,nz,alpha,
     2   ipc,rpc,ac,cc,fc,p,x,r,ap,zk,zkp1)
c* *********************************************************************
c* purpose:
c*
c*    line search for a minimizing steplength given a descent direction.
c*
c* input:
c*
c*    nx,ny,nz = mesh dimensions
c*    ipc,rpc  = operator/level information
c*    ac,cc,fc = operator/level entries
c*    p        = the direction 
c*    x        = the previous solution
c*    r        = the previous nonlinear residual, r=r(x)
c*
c* temporary vectors:
c*
c*    ap,zk,zkp1 = used as temp vectors, passed back as needed for the
c*                 nonlinear fletcher-reeves conjugate gradient algorithm
c*                 (can be ignored otherwise)
c*    r          = is definitely overwritten
c*
c* output:
c*
c*    alpha = steplength for minimizing the 1d problem in direction p.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
c*
c*    *** other declarations ***
      integer          ipc(*),nx,ny,nz,ipkey
      double precision rpc(*),ac(*),cc(*),fc(*)
      double precision p(*),x(*),r(*),ap(*),zk(*),zkp1(*)
      double precision alpha,xdot
      double precision ALP,pAp,pr,val,xdot3,ALP_old
      double precision zNN,zntol,DzNN
      integer          niters,nitmax,ifail_tol
c*
c*    *** nonlinear iteration tolerance and itmax ***
      nitmax    = 20
      zntol     = 1.0e-5
      ifail_tol = 0
      ipkey     = ipc(10)
c*
c* *********************************************************************
c* *** begin linear case for choosing alpha
c* *********************************************************************
c*
c*    *** compute the hackbusch/reusken damping parameter ***
c*    *** which is equivalent to the standard linear cg steplength ***
      call matvec(nx,ny,nz,ipc,rpc,ac,cc,p,ap)
      pr  = xdot(nx,ny,nz,p,r)
      pAp = xdot(nx,ny,nz,p,ap)
      alpha = pr / pAp
c*
c* *********************************************************************
c* *** begin nonlinear case for choosing alpha
c* *********************************************************************
c*
c*    *** initialize ALP (use linear approx) and zkp1 ***
      ALP = alpha
      call xcopy(nx,ny,nz,x,r)
      call xaxpy(nx,ny,nz,ALP,p,r)
      call c_vec(cc,x,zk,nx,ny,nz,ipkey)
      call c_vec(cc,r,zkp1,nx,ny,nz,ipkey)
c*
c*    *** evaluate residual of 1d system
      call xcopy(nx,ny,nz,zkp1,r)
      call xaxpy(nx,ny,nz,(-1.0d0),zk,r)
      val = xdot(nx,ny,nz,p,r)
      zNN = ALP * pAp - pr + val
c*
c*    *** setup for newton's method
      niters = 0
c*
c*    *** the 1d newton's method ***
      if (niters .gt. nitmax) goto 16
      ifail_tol = 0
      ALP_old = ALP
 15   continue
         niters = niters + 1
c*
c*       *** some i/o ***
c* **    write(6,110)
c* **2      'LINESEARCH: 1d newton iterations: niters, alpha = ',
c* **3      niters,alpha
c* **       format('% ',a,i5,1pe10.2)
c*
c*       *** construct jacobian matrix of NN ***
         call dc_vec(cc,r,zkp1,nx,ny,nz,ipkey)
         val = xdot3(nx,ny,nz,zkp1,p)
         DzNN = pAp + val
c*
c*       *** update the solution ***
         ALP_old = ALP
         ALP = ALP - zNN / DzNN
c*
c*       *** evaluate residual of 1d system ***
         call xcopy(nx,ny,nz,x,r)
         call xaxpy(nx,ny,nz,ALP,p,r)
         call c_vec(cc,r,zkp1,nx,ny,nz,ipkey)
         call xcopy(nx,ny,nz,zkp1,r)
         call xaxpy(nx,ny,nz,(-1.0d0),zk,r)
         val = xdot(nx,ny,nz,p,r)
         zNN = ALP * pAp - pr + val
c*
c*       *** tolerance and itmax check ***
         if ((dabs(zNN) .lt. zntol) .and. 
     2       (dabs((ALP-ALP_old)/ALP) .lt. zntol) ) goto 17
         if (niters .gt. nitmax) goto 16
c*
c*    *** loop ***
      goto 15
c*
c*    *** tolerance not reached ***
 16   continue
      ifail_tol = ifail_tol + 1
c*
c*    *** tolerance reached ***
 17   continue
c*
c*    *** newton's method complete -- update solution value ***
      alpha = ALP
c*
c*    *** messages ***
      if (ifail_tol .gt. 0) then
         call vnmpri(2,'% LINESEARCH: exceeded 1d newton itmax: = ',
     2      42, nitmax)
      endif
c*
c*    *** return and end ***
      return
      end
