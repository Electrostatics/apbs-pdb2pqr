c* ///////////////////////////////////////////////////////////////////////////
c* @file    matvecd.f
c* @author  Michael Holst
c* @brief   Matrix-vector operations.
c* @version $Id: matvecd.f 868 2006-04-05 18:32:05Z todd_dolinsky $
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

      subroutine matvec(nx,ny,nz,ipc,rpc,ac,cc,x,y)
c* *********************************************************************
c* purpose:
c*
c*    break the matrix data-struCture into diagoNals
c*    and then call the matrix-vector routine.
c*    
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),nx,ny,nz,numdia
      double precision rpc(*),ac(nx*ny*nz,*),cc(*),x(*),y(*)
c*
cmdir 0 0
c*
c*    *** do in oNe step ***
      numdia = ipc(11)
      if (numdia .eq. 7) then
         call matvec7(nx,ny,nz,ipc,rpc,ac,cc,x,y)
      elseif (numdia .eq. 27) then
         call matvec27(nx,ny,nz,ipc,rpc,ac,cc,x,y)
      else
         call vnmprt(2,'% MATVEC: invalid stencil type given...',39)
      endif
c*
c*    *** return and end ***
      return
      end
      subroutine matvec7(nx,ny,nz,ipc,rpc,ac,cc,x,y)
c* *********************************************************************
c* purpose:
c*
c*    break the matrix data-struCture into diagoNals
c*    and then call the matrix-vector routine.
c*    
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),nx,ny,nz
      double precision rpc(*),ac(nx*ny*nz,*),cc(*),x(*),y(*)
c*
cmdir 0 0
c*
c*    *** do in oNe step ***
      call matvec7_1s(nx,ny,nz,ipc,rpc,ac(1,1),cc,
     2   ac(1,2),ac(1,3),ac(1,4),
     3   x,y)
c*
c*    *** return and end ***
      return
      end
      subroutine matvec7_1s(nx,ny,nz,ipc,rpc,oC,cc,oE,oN,uC,x,y)
c* *********************************************************************
c* purpose:
c*
c*    compute the operator times a vector.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),nx,ny,nz,i,j,k
      double precision oE(nx,ny,nz),oN(nx,ny,nz),uC(nx,ny,nz)
      double precision cc(nx,ny,nz)
      double precision oC(nx,ny,nz)
      double precision x(nx,ny,nz),y(nx,ny,nz),rpc(*)
c*
cmdir 0 0
c*
c*    *** do it ***
cmdir 3 1
      do 10 k=2,nz-1
cmdir 3 2
         do 11 j=2,ny-1
cmdir 3 3
            do 12 i=2,nx-1
               y(i,j,k) = 
     2               -  oN(i,j,k)        * x(i,j+1,k)
     3               -  oN(i,j-1,k)      * x(i,j-1,k)
     4               -  oE(i,j,k)        * x(i+1,j,k)
     5               -  oE(i-1,j,k)      * x(i-1,j,k)
     6               -  uC(i,j,k-1)      * x(i,j,k-1)
     7               -  uC(i,j,k)        * x(i,j,k+1)
     8               +  (oC(i,j,k) + cc(i,j,k)) * x(i,j,k)
 12         continue
 11      continue
 10   continue
c*
c*    *** return and end ***
      return
      end
      subroutine matvec27(nx,ny,nz,ipc,rpc,ac,cc,x,y)
c* *********************************************************************
c* purpose:
c*
c*    break the matrix data-struCture into diagoNals
c*    and then call the matrix-vector routine.
c*    
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),nx,ny,nz
      double precision rpc(*),ac(nx*ny*nz,*),cc(*),x(*),y(*)
c*
cmdir 0 0
c*
c*    *** do in oNe step ***
      call matvec27_1s(nx,ny,nz,ipc,rpc,ac(1,1),cc,
     2   ac(1,2),ac(1,3),ac(1,4),ac(1,5),ac(1,6),ac(1,7),ac(1,8),
     3   ac(1,9),ac(1,10),ac(1,11),ac(1,12),ac(1,13),ac(1,14),
     4   x,y)
c*
c*    *** return and end ***
      return
      end
      subroutine matvec27_1s(nx,ny,nz,ipc,rpc,oC,cc,
     2   oE,oN,uC,oNE,oNW,uE,uW,uN,uS,uNE,uNW,uSE,uSW,
     3   x,y)
c* *********************************************************************
c* purpose:
c*
c*    compute the operator times a vector.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),nx,ny,nz,i,j,k
      double precision oE(nx,ny,nz),oN(nx,ny,nz),uC(nx,ny,nz)
      double precision oNE(nx,ny,nz),oNW(nx,ny,nz),uE(nx,ny,nz)
      double precision uW(nx,ny,nz),uN(nx,ny,nz),uS(nx,ny,nz)
      double precision uNE(nx,ny,nz),uNW(nx,ny,nz),uSE(nx,ny,nz)
      double precision uSW(nx,ny,nz)
      double precision cc(nx,ny,nz)
      double precision oC(nx,ny,nz)
      double precision x(nx,ny,nz),y(nx,ny,nz),rpc(*)
      double precision tmpO,tmpU,tmpD
c*
cmdir 0 0
c*
c*    *** do it ***
cmdir 3 1
      do 10 k=2,nz-1
cmdir 3 2
         do 11 j=2,ny-1
cmdir 3 3
            do 12 i=2,nx-1
                  tmpO =
     2               -  oN(i,j,k)        * x(i,j+1,k)
     3               -  oN(i,j-1,k)      * x(i,j-1,k)
     4               -  oE(i,j,k)        * x(i+1,j,k)
     5               -  oE(i-1,j,k)      * x(i-1,j,k)
     6               -  oNE(i,j,k)       * x(i+1,j+1,k)
     7               -  oNW(i,j,k)       * x(i-1,j+1,k)
     8               -  oNW(i+1,j-1,k)   * x(i+1,j-1,k)
     9               -  oNE(i-1,j-1,k)   * x(i-1,j-1,k)
                  tmpU =
     2               -  uC(i,j,k)        * x(i,j,k+1)
     3               -  uN(i,j,k)        * x(i,j+1,k+1)
     4               -  uS(i,j,k)        * x(i,j-1,k+1)
     5               -  uE(i,j,k)        * x(i+1,j,k+1)
     6               -  uW(i,j,k)        * x(i-1,j,k+1)
     7               -  uNE(i,j,k)       * x(i+1,j+1,k+1)
     8               -  uNW(i,j,k)       * x(i-1,j+1,k+1)
     9               -  uSE(i,j,k)       * x(i+1,j-1,k+1)
     9               -  uSW(i,j,k)       * x(i-1,j-1,k+1)
                  tmpD =
     2               -  uC(i,j,k-1)      * x(i,j,k-1)
     3               -  uS(i,j+1,k-1)    * x(i,j+1,k-1)
     4               -  uN(i,j-1,k-1)    * x(i,j-1,k-1)
     5               -  uW(i+1,j,k-1)    * x(i+1,j,k-1)
     6               -  uE(i-1,j,k-1)    * x(i-1,j,k-1)
     7               -  uSW(i+1,j+1,k-1) * x(i+1,j+1,k-1)
     8               -  uSE(i-1,j+1,k-1) * x(i-1,j+1,k-1)
     9               -  uNW(i+1,j-1,k-1) * x(i+1,j-1,k-1)
     9               -  uNE(i-1,j-1,k-1) * x(i-1,j-1,k-1)
               y(i,j,k) = tmpO + tmpU + tmpD
     2               +  (oC(i,j,k) + cc(i,j,k)) * x(i,j,k)
 12         continue
 11      continue
 10   continue
c*
c*    *** return and end ***
      return
      end
      subroutine mresid(nx,ny,nz,ipc,rpc,ac,cc,fc,x,r)
c* *********************************************************************
c* purpose:
c*
c*    break the matrix data-struCture into diagoNals
c*    and then call the residual routine.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),nx,ny,nz,numdia
      double precision rpc(*),ac(nx*ny*nz,*),cc(*),fc(*),x(*),r(*)
c*
cmdir 0 0
c*
c*    *** do in oNe step ***
      numdia = ipc(11)
      if (numdia .eq. 7) then
         call mresid7(nx,ny,nz,ipc,rpc,ac,cc,fc,x,r)
      elseif (numdia .eq. 27) then
         call mresid27(nx,ny,nz,ipc,rpc,ac,cc,fc,x,r)
      else
         call vnmprt(2,'% MRESID: invalid stencil type given...',39)
      endif
c*
c*    *** return and end ***
      return
      end
      subroutine mresid7(nx,ny,nz,ipc,rpc,ac,cc,fc,x,r)
c* *********************************************************************
c* purpose:
c*
c*    break the matrix data-struCture into diagoNals
c*    and then call the residual routine.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),nx,ny,nz
      double precision rpc(*),ac(nx*ny*nz,*),cc(*),fc(*),x(*),r(*)
c*
cmdir 0 0
c*
c*    *** do in oNe step ***
      call mresid7_1s(nx,ny,nz,ipc,rpc,ac(1,1),cc,fc,
     2   ac(1,2),ac(1,3),ac(1,4),
     3   x,r)
c*
c*    *** return and end ***
      return
      end
	  
      subroutine mresid7_1s(nx,ny,nz,ipc,rpc,oC,cc,fc,oE,oN,uC,x,r)
c* *********************************************************************
c* purpose:
c*
c*    compute the residual.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),nx,ny,nz,i,j,k
      double precision oE(nx,ny,nz),oN(nx,ny,nz),uC(nx,ny,nz)
      double precision cc(nx,ny,nz),fc(nx,ny,nz)
      double precision oC(nx,ny,nz)
      double precision x(nx,ny,nz),r(nx,ny,nz),rpc(*)
c*
cmdir 0 0
c*
c*    *** do it ***
!$OMP PARALLEL DO default(shared) private(i,j,k)
       do k=2,nz-1
         do j=2,ny-1
            do i=2,nx-1
               r(i,j,k) = fc(i,j,k)
     2               +  oN(i,j,k)        * x(i,j+1,k)
     3               +  oN(i,j-1,k)      * x(i,j-1,k)
     4               +  oE(i,j,k)        * x(i+1,j,k)
     5               +  oE(i-1,j,k)      * x(i-1,j,k)
     6               +  uC(i,j,k-1)      * x(i,j,k-1)
     7               +  uC(i,j,k)        * x(i,j,k+1)
     8               -  (oC(i,j,k) + cc(i,j,k)) * x(i,j,k)
			end do
		end do
	  end do
!$OMP END PARALLEL DO
c*
c*    *** return and end ***
      return
      end
	  
      subroutine mresid27(nx,ny,nz,ipc,rpc,ac,cc,fc,x,r)
c* *********************************************************************
c* purpose:
c*
c*    break the matrix data-struCture into diagoNals
c*    and then call the residual routine.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),nx,ny,nz
      double precision rpc(*),ac(nx*ny*nz,*),cc(*),fc(*),x(*),r(*)
c*
cmdir 0 0
c*
c*    *** do in oNe step ***
      call mresid27_1s(nx,ny,nz,ipc,rpc,ac(1,1),cc,fc,
     2   ac(1,2),ac(1,3),ac(1,4),ac(1,5),ac(1,6),ac(1,7),ac(1,8),
     3   ac(1,9),ac(1,10),ac(1,11),ac(1,12),ac(1,13),ac(1,14),
     4   x,r)
c*
c*    *** return and end ***
      return
      end
      subroutine mresid27_1s(nx,ny,nz,ipc,rpc,oC,cc,fc,
     2   oE,oN,uC,oNE,oNW,uE,uW,uN,uS,uNE,uNW,uSE,uSW,
     3   x,r)
c* *********************************************************************
c* purpose:
c*
c*    compute the residual.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),nx,ny,nz,i,j,k
      double precision oE(nx,ny,nz),oN(nx,ny,nz),uC(nx,ny,nz)
      double precision oNE(nx,ny,nz),oNW(nx,ny,nz),uE(nx,ny,nz)
      double precision uW(nx,ny,nz),uN(nx,ny,nz),uS(nx,ny,nz)
      double precision uNE(nx,ny,nz),uNW(nx,ny,nz),uSE(nx,ny,nz)
      double precision uSW(nx,ny,nz)
      double precision cc(nx,ny,nz),fc(nx,ny,nz)
      double precision oC(nx,ny,nz)
      double precision x(nx,ny,nz),r(nx,ny,nz),rpc(*)
      double precision tmpO,tmpU,tmpD

!$OMP PARALLEL default(shared) private(i,j,k,tmpO,tmpU,tmpD)
!$OMP DO
      do 10 k=2,nz-1
         do 11 j=2,ny-1
            do 12 i=2,nx-1
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
               r(i,j,k) = fc(i,j,k) + tmpO + tmpU + tmpD
     2               -  (oC(i,j,k) + cc(i,j,k)) * x(i,j,k)
 12         continue
 11      continue
 10   continue
!$OMP END DO
!$OMP END PARALLEL

      return
      end
      subroutine nmatvec(nx,ny,nz,ipc,rpc,ac,cc,x,y,w1)
c* *********************************************************************
c* purpose:
c*
c*    break the matrix data-struCture into diagoNals
c*    and then call the matrix-vector routine.
c*    
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),nx,ny,nz,numdia
      double precision rpc(*),ac(nx*ny*nz,*),cc(*),x(*),y(*)
      double precision w1(*)
c*
cmdir 0 0
c*
c*    *** do in oNe step ***
      numdia = ipc(11)
      if (numdia .eq. 7) then
         call nmatvec7(nx,ny,nz,ipc,rpc,ac,cc,x,y,w1)
      elseif (numdia .eq. 27) then
         call nmatvec27(nx,ny,nz,ipc,rpc,ac,cc,x,y,w1)
      else
         call vnmprt(2,'% NMATVEC: invalid stencil type given...',40)
      endif
c*
c*    *** return and end ***
      return
      end
      subroutine nmatvec7(nx,ny,nz,ipc,rpc,ac,cc,x,y,w1)
c* *********************************************************************
c* purpose:
c*
c*    break the matrix data-struCture into diagoNals
c*    and then call the matrix-vector routine.
c*    
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),nx,ny,nz
      double precision rpc(*),ac(nx*ny*nz,*),cc(*),x(*),y(*)
      double precision w1(*)
c*
cmdir 0 0
c*
c*    *** do in oNe step ***
      call nmatvecd7_1s(nx,ny,nz,ipc,rpc,ac(1,1),cc,
     2   ac(1,2),ac(1,3),ac(1,4),
     3   x,y,w1)
c*
c*    *** return and end ***
      return
      end
      subroutine nmatvecd7_1s(nx,ny,nz,ipc,rpc,
     2   oC,cc,
     3   oE,oN,uC,x,y,w1)
c* *********************************************************************
c* purpose:
c*
c*    compute the operator times a vector.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),nx,ny,nz,i,j,k,ipkey
      double precision oE(nx,ny,nz),oN(nx,ny,nz),uC(nx,ny,nz)
      double precision cc(nx,ny,nz)
      double precision oC(nx,ny,nz)
      double precision x(nx,ny,nz),y(nx,ny,nz),rpc(*),w1(nx,ny,nz)
c*
cmdir 0 0
c*
c*    *** first get vector noNlinear term to avoid subroutine calls ***
      ipkey = ipc(10)
      call c_vec(cc,x,w1,nx,ny,nz,ipkey)
c*
c*    *** the operator ***
cmdir 3 1
      do 10 k=2,nz-1
cmdir 3 2
         do 11 j=2,ny-1
cmdir 3 3
            do 12 i=2,nx-1
               y(i,j,k) = 
     2               -  oN(i,j,k)        * x(i,j+1,k)
     3               -  oN(i,j-1,k)      * x(i,j-1,k)
     4               -  oE(i,j,k)        * x(i+1,j,k)
     5               -  oE(i-1,j,k)      * x(i-1,j,k)
     6               -  uC(i,j,k-1)      * x(i,j,k-1)
     7               -  uC(i,j,k)        * x(i,j,k+1)
     8               +  oC(i,j,k) * x(i,j,k) +  w1(i,j,k)
 12         continue
 11      continue
 10   continue
c*
c*    *** return and end ***
      return
      end
      subroutine nmatvec27(nx,ny,nz,ipc,rpc,ac,cc,x,y,w1)
c* *********************************************************************
c* purpose:
c*
c*    break the matrix data-struCture into diagoNals
c*    and then call the matrix-vector routine.
c*    
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),nx,ny,nz
      double precision rpc(*),ac(nx*ny*nz,*),cc(*),x(*),y(*)
      double precision w1(*)
c*
cmdir 0 0
c*
c*    *** do in oNe step ***
      call nmatvecd27_1s(nx,ny,nz,ipc,rpc,ac(1,1),cc,
     2   ac(1,2),ac(1,3),ac(1,4),ac(1,5),ac(1,6),ac(1,7),ac(1,8),
     3   ac(1,9),ac(1,10),ac(1,11),ac(1,12),ac(1,13),ac(1,14),
     4   x,y,w1)
c*
c*    *** return and end ***
      return
      end
      subroutine nmatvecd27_1s(nx,ny,nz,ipc,rpc,
     2   oC,cc,
     3   oE,oN,uC,oNE,oNW,uE,uW,uN,uS,uNE,uNW,uSE,uSW,
     4   x,y,w1)
c* *********************************************************************
c* purpose:
c*
c*    compute the operator times a vector.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),nx,ny,nz,i,j,k,ipkey
      double precision oE(nx,ny,nz),oN(nx,ny,nz),uC(nx,ny,nz)
      double precision oNE(nx,ny,nz),oNW(nx,ny,nz),uE(nx,ny,nz)
      double precision uW(nx,ny,nz),uN(nx,ny,nz),uS(nx,ny,nz)
      double precision uNE(nx,ny,nz),uNW(nx,ny,nz),uSE(nx,ny,nz)
      double precision uSW(nx,ny,nz)
      double precision cc(nx,ny,nz)
      double precision oC(nx,ny,nz)
      double precision x(nx,ny,nz),y(nx,ny,nz),rpc(*),w1(nx,ny,nz)
      double precision tmpO,tmpU,tmpD
c*
cmdir 0 0
c*
c*    *** first get vector noNlinear term to avoid subroutine calls ***
      ipkey = ipc(10)
      call c_vec(cc,x,w1,nx,ny,nz,ipkey)
c*
c*    *** the operator ***
cmdir 3 1
      do 10 k=2,nz-1
cmdir 3 2
         do 11 j=2,ny-1
cmdir 3 3
            do 12 i=2,nx-1
                  tmpO =
     2               -  oN(i,j,k)        * x(i,j+1,k)
     3               -  oN(i,j-1,k)      * x(i,j-1,k)
     4               -  oE(i,j,k)        * x(i+1,j,k)
     5               -  oE(i-1,j,k)      * x(i-1,j,k)
     6               -  oNE(i,j,k)       * x(i+1,j+1,k)
     7               -  oNW(i,j,k)       * x(i-1,j+1,k)
     8               -  oNW(i+1,j-1,k)   * x(i+1,j-1,k)
     9               -  oNE(i-1,j-1,k)   * x(i-1,j-1,k)
                  tmpU =
     2               -  uC(i,j,k)        * x(i,j,k+1)
     3               -  uN(i,j,k)        * x(i,j+1,k+1)
     4               -  uS(i,j,k)        * x(i,j-1,k+1)
     5               -  uE(i,j,k)        * x(i+1,j,k+1)
     6               -  uW(i,j,k)        * x(i-1,j,k+1)
     7               -  uNE(i,j,k)       * x(i+1,j+1,k+1)
     8               -  uNW(i,j,k)       * x(i-1,j+1,k+1)
     9               -  uSE(i,j,k)       * x(i+1,j-1,k+1)
     9               -  uSW(i,j,k)       * x(i-1,j-1,k+1)
                  tmpD =
     2               -  uC(i,j,k-1)      * x(i,j,k-1)
     3               -  uS(i,j+1,k-1)    * x(i,j+1,k-1)
     4               -  uN(i,j-1,k-1)    * x(i,j-1,k-1)
     5               -  uW(i+1,j,k-1)    * x(i+1,j,k-1)
     6               -  uE(i-1,j,k-1)    * x(i-1,j,k-1)
     7               -  uSW(i+1,j+1,k-1) * x(i+1,j+1,k-1)
     8               -  uSE(i-1,j+1,k-1) * x(i-1,j+1,k-1)
     9               -  uNW(i+1,j-1,k-1) * x(i+1,j-1,k-1)
     9               -  uNE(i-1,j-1,k-1) * x(i-1,j-1,k-1)
               y(i,j,k) = tmpO + tmpU + tmpD
     2               +  oC(i,j,k) * x(i,j,k) +  w1(i,j,k)
 12         continue
 11      continue
 10   continue
c*
c*    *** return and end ***
      return
      end
      subroutine nmresid(nx,ny,nz,ipc,rpc,ac,cc,fc,x,r,w1)
c* *********************************************************************
c* purpose:
c*
c*    break the matrix data-struCture into diagoNals
c*    and then call the residual routine.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),nx,ny,nz,numdia
      double precision rpc(*),ac(nx*ny*nz,*),cc(*),fc(*),x(*),r(*)
      double precision w1(*)
c*
cmdir 0 0
c*
c*    *** do in oNe step ***
      numdia = ipc(11)
      if (numdia .eq. 7) then
         call nmresid7(nx,ny,nz,ipc,rpc,ac,cc,fc,x,r,w1)
      elseif (numdia .eq. 27) then
         call nmresid27(nx,ny,nz,ipc,rpc,ac,cc,fc,x,r,w1)
      else
         call vnmprt(2,'% NMRESID: invalid stencil type given...',40)
      endif
c*
c*    *** return and end ***
      return
      end
      subroutine nmresid7(nx,ny,nz,ipc,rpc,ac,cc,fc,x,r,w1)
c* *********************************************************************
c* purpose:
c*
c*    break the matrix data-struCture into diagoNals
c*    and then call the residual routine.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),nx,ny,nz
      double precision rpc(*),ac(nx*ny*nz,*),cc(*),fc(*),x(*),r(*)
      double precision w1(*)
c*
cmdir 0 0
c*
c*    *** do in oNe step ***
      call nmresid7_1s(nx,ny,nz,ipc,rpc,ac(1,1),cc,fc,
     2   ac(1,2),ac(1,3),ac(1,4),
     3   x,r,w1)
c*
c*    *** return and end ***
      return
      end
      subroutine nmresid7_1s(nx,ny,nz,ipc,rpc,
     2   oC,cc,fc,
     3   oE,oN,uC,x,r,w1)
c* *********************************************************************
c* purpose:
c*
c*    compute the residual.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),nx,ny,nz,i,j,k,ipkey
      double precision oE(nx,ny,nz),oN(nx,ny,nz),uC(nx,ny,nz)
      double precision cc(nx,ny,nz),fc(nx,ny,nz)
      double precision oC(nx,ny,nz)
      double precision x(nx,ny,nz),r(nx,ny,nz),rpc(*),w1(nx,ny,nz)
c*
cmdir 0 0
c*
c*    *** first get vector noNlinear term to avoid subroutine calls ***
      ipkey = ipc(10)
      call c_vec(cc,x,w1,nx,ny,nz,ipkey)
c*
c*    *** the residual ***
cmdir 3 1
      do 10 k=2,nz-1
cmdir 3 2
         do 11 j=2,ny-1
cmdir 3 3
            do 12 i=2,nx-1
               r(i,j,k) = fc(i,j,k)
     2               +  oN(i,j,k)        * x(i,j+1,k)
     3               +  oN(i,j-1,k)      * x(i,j-1,k)
     4               +  oE(i,j,k)        * x(i+1,j,k)
     5               +  oE(i-1,j,k)      * x(i-1,j,k)
     6               +  uC(i,j,k-1)      * x(i,j,k-1)
     7               +  uC(i,j,k)        * x(i,j,k+1)
     8               -  oC(i,j,k) * x(i,j,k) -  w1(i,j,k)
 12         continue
 11      continue
 10   continue
c*
c*    *** return and end ***
      return
      end
      subroutine nmresid27(nx,ny,nz,ipc,rpc,ac,cc,fc,x,r,w1)
c* *********************************************************************
c* purpose:
c*
c*    break the matrix data-struCture into diagoNals
c*    and then call the residual routine.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),nx,ny,nz
      double precision rpc(*),ac(nx*ny*nz,*),cc(*),fc(*),x(*),r(*)
      double precision w1(*)
c*
cmdir 0 0
c*
c*    *** do in oNe step ***
      call nmresid27_1s(nx,ny,nz,ipc,rpc,ac(1,1),cc,fc,
     2   ac(1,2),ac(1,3),ac(1,4),ac(1,5),ac(1,6),ac(1,7),ac(1,8),
     3   ac(1,9),ac(1,10),ac(1,11),ac(1,12),ac(1,13),ac(1,14),
     4   x,r,w1)
c*
c*    *** return and end ***
      return
      end
      subroutine nmresid27_1s(nx,ny,nz,ipc,rpc,
     2   oC,cc,fc,
     3   oE,oN,uC,oNE,oNW,uE,uW,uN,uS,uNE,uNW,uSE,uSW,
     4   x,r,w1)
c* *********************************************************************
c* purpose:
c*
c*    compute the residual.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),nx,ny,nz,i,j,k,ipkey
      double precision oE(nx,ny,nz),oN(nx,ny,nz),uC(nx,ny,nz)
      double precision oNE(nx,ny,nz),oNW(nx,ny,nz),uE(nx,ny,nz)
      double precision uW(nx,ny,nz),uN(nx,ny,nz),uS(nx,ny,nz)
      double precision uNE(nx,ny,nz),uNW(nx,ny,nz),uSE(nx,ny,nz)
      double precision uSW(nx,ny,nz)
      double precision cc(nx,ny,nz),fc(nx,ny,nz)
      double precision oC(nx,ny,nz)
      double precision x(nx,ny,nz),r(nx,ny,nz),rpc(*),w1(nx,ny,nz)
      double precision tmpO,tmpU,tmpD
c*
cmdir 0 0
c*
c*    *** first get vector noNlinear term to avoid subroutine calls ***
      ipkey = ipc(10)
      call c_vec(cc,x,w1,nx,ny,nz,ipkey)
c*
c*    *** the residual ***
cmdir 3 1
      do 10 k=2,nz-1
cmdir 3 2
         do 11 j=2,ny-1
cmdir 3 3
            do 12 i=2,nx-1
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
               r(i,j,k) = fc(i,j,k) + tmpO + tmpU + tmpD
     2               -  oC(i,j,k) * x(i,j,k) -  w1(i,j,k)
 12         continue
 11      continue
 10   continue
c*
c*    *** return and end ***
      return
      end
      subroutine interpPMG(nxc,nyc,nzc,nxf,nyf,nzf,xin,xout,pc)
c* *********************************************************************
c* purpose:
c*
c*    apply the prolongation operator.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          nxf,nyf,nzf,nxc,nyc,nzc
      double precision xin(*),xout(*),pc(nxc*nyc*nzc,*)
c*
c*    *** doit ***
CZZZZ call interpPMGZ(nxc,nyc,nzc,nxf,nyf,nzf,xin,xout)
      call interpPMG2(nxc,nyc,nzc,nxf,nyf,nzf,xin,xout,
     2   pc(1,1),pc(1,2),pc(1,3),pc(1,4),pc(1,5),pc(1,6),pc(1,7),
     3   pc(1,8),pc(1,9),pc(1,10),pc(1,11),pc(1,12),pc(1,13),pc(1,14),
     4   pc(1,15),pc(1,16),pc(1,17),pc(1,18),pc(1,19),pc(1,20),pc(1,21),
     5   pc(1,22),pc(1,23),pc(1,24),pc(1,25),pc(1,26),pc(1,27))
c*
c*    *** return and end ***
      return
      end
      subroutine interpPMG2(nxc,nyc,nzc,nxf,nyf,nzf,xin,xout,
     2   oPC,oPN,oPS,oPE,oPW,oPNE,oPNW,oPSE,oPSW,
     3   uPC,uPN,uPS,uPE,uPW,uPNE,uPNW,uPSE,uPSW,
     4   dPC,dPN,dPS,dPE,dPW,dPNE,dPNW,dPSE,dPSW)
c* *********************************************************************
c* purpose:
c*
c*    apply the prolongation operator.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          nxc,nyc,nzc,nxf,nyf,nzf,i,j,k,jj,ii,kk
      double precision xin(nxc,nyc,nzc),xout(nxf,nyf,nzf)
      double precision oPC(nxc,nyc,nzc),oPN(nxc,nyc,nzc)
      double precision oPS(nxc,nyc,nzc),oPE(nxc,nyc,nzc)
      double precision oPW(nxc,nyc,nzc),oPNE(nxc,nyc,nzc)
      double precision oPNW(nxc,nyc,nzc),oPSE(nxc,nyc,nzc)
      double precision oPSW(nxc,nyc,nzc),uPC(nxc,nyc,nzc)
      double precision uPN(nxc,nyc,nzc),uPS(nxc,nyc,nzc)
      double precision uPE(nxc,nyc,nzc),uPW(nxc,nyc,nzc)
      double precision uPNE(nxc,nyc,nzc),uPNW(nxc,nyc,nzc)
      double precision uPSE(nxc,nyc,nzc),uPSW(nxc,nyc,nzc)
      double precision dPC(nxc,nyc,nzc),dPN(nxc,nyc,nzc)
      double precision dPS(nxc,nyc,nzc),dPE(nxc,nyc,nzc)
      double precision dPW(nxc,nyc,nzc),dPNE(nxc,nyc,nzc)
      double precision dPNW(nxc,nyc,nzc),dPSE(nxc,nyc,nzc)
      double precision dPSW(nxc,nyc,nzc)
c*
c* *********************************************************************
c* setup 
c* *********************************************************************
c*
cmdir 0 0
c*
c*    *** verify correctness of the input boundary points ***
      call fboundPMG00(nxc,nyc,nzc,xin)
c*
c*    *** doit ***
cmdir 3 1
      do 10 k = 1, nzf-2, 2
         kk = (k-1) / 2 + 1
cmdir 3 2
         do 11 j = 1, nyf-2, 2
            jj = (j-1) / 2 + 1
cmdir 3 3
            do 12 i = 1, nxf-2, 2
               ii = (i-1) / 2 + 1
c*
c* *********************************************************************
c* type 1 -- fine grid points common to a coarse grid point
c* *********************************************************************
c*
c*             *** copy coinciding points from coarse grid to fine grid ***
               xout(i,j,k) = xin(ii,jj,kk)
c*
c* *********************************************************************
c* type 2 -- fine grid points common to a coarse grid plane
c* *********************************************************************
c*
c*             *** fine grid pts common only to y-z planes on coarse grid ***
c*             *** (intermediate pts between 2 grid points on x-row) ***
               xout(i+1,j,k) = oPE(ii,jj,kk)   * xin(ii,jj,kk) 
     2                       + oPW(ii+1,jj,kk) * xin(ii+1,jj,kk)
c*
c*             *** fine grid pts common only to x-z planes on coarse grid ***
c*             *** (intermediate pts between 2 grid points on a y-row) ***
               xout(i,j+1,k) = oPN(ii,jj,kk)   * xin(ii,jj,kk)
     2                       + oPS(ii,jj+1,kk) * xin(ii,jj+1,kk) 
c*
c*             *** fine grid pts common only to x-y planes on coarse grid ***
c*             *** (intermediate pts between 2 grid points on a z-row) ***
               xout(i,j,k+1) = uPC(ii,jj,kk)   * xin(ii,jj,kk)
     2                       + dPC(ii,jj,kk+1) * xin(ii,jj,kk+1) 
c*
c* *********************************************************************
c* type 3 -- fine grid points common to a coarse grid line
c* *********************************************************************
c*
c*             *** fine grid pts common only to z planes on coarse grid ***
c*             *** (intermediate pts between 4 grid pts on the xy-plane***
               xout(i+1,j+1,k) = oPNE(ii,jj,kk)     * xin(ii,jj,kk)   
     2                         + oPNW(ii+1,jj,kk)   * xin(ii+1,jj,kk)
     3                         + oPSE(ii,jj+1,kk)   * xin(ii,jj+1,kk)
     4                         + oPSW(ii+1,jj+1,kk) * xin(ii+1,jj+1,kk)
c*
c*             *** fine grid pts common only to y planes on coarse grid ***
c*             *** (intermediate pts between 4 grid pts on the xz-plane***
               xout(i+1,j,k+1) = uPE(ii,jj,kk)     * xin(ii,jj,kk)
     2                         + uPW(ii+1,jj,kk)   * xin(ii+1,jj,kk)
     3                         + dPE(ii,jj,kk+1)   * xin(ii,jj,kk+1)
     4                         + dPW(ii+1,jj,kk+1) * xin(ii+1,jj,kk+1)
c*
c*             *** fine grid pts common only to x planes on coarse grid ***
c*             *** (intermediate pts between 4 grid pts on the yz-plane***
               xout(i,j+1,k+1) = uPN(ii,jj,kk)     * xin(ii,jj,kk)
     2                         + uPS(ii,jj+1,kk)   * xin(ii,jj+1,kk)
     3                         + dPN(ii,jj,kk+1)   * xin(ii,jj,kk+1)
     4                         + dPS(ii,jj+1,kk+1) * xin(ii,jj+1,kk+1)
c*
c* *********************************************************************
c* type 4 -- fine grid points not common to coarse grid pts/lines/planes
c* *********************************************************************
c*
c*             *** completely interior points ***
               xout(i+1,j+1,k+1) = 
     2            + uPNE(ii,jj,kk)       * xin(ii,jj,kk)
     3            + uPNW(ii+1,jj,kk)     * xin(ii+1,jj,kk)
     4            + uPSE(ii,jj+1,kk)     * xin(ii,jj+1,kk)
     5            + uPSW(ii+1,jj+1,kk)   * xin(ii+1,jj+1,kk)
     6            + dPNE(ii,jj,kk+1)     * xin(ii,jj,kk+1)
     7            + dPNW(ii+1,jj,kk+1)   * xin(ii+1,jj,kk+1)
     8            + dPSE(ii,jj+1,kk+1)   * xin(ii,jj+1,kk+1)
     9            + dPSW(ii+1,jj+1,kk+1) * xin(ii+1,jj+1,kk+1)
 12         continue
 11      continue
 10   continue
c*
c*    *** verify correctness of the output boundary points ***
      call fboundPMG00(nxf,nyf,nzf,xout)
c*
c*    *** return and end ***
      return
      end
      subroutine restrc(nxf,nyf,nzf,nxc,nyc,nzc,xin,xout,pc)
c* *********************************************************************
c* purpose:
c*
c*    apply the restriction operator.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          nxf,nyf,nzf,nxc,nyc,nzc
      double precision xin(*),xout(*),pc(nxc*nyc*nzc,*)
c*
c*    *** doit ***
CZZZ  call restrcZ(nxf,nyf,nzf,nxc,nyc,nzc,xin,xout)
      call restrc2(nxf,nyf,nzf,nxc,nyc,nzc,xin,xout,
     2   pc(1,1),pc(1,2),pc(1,3),pc(1,4),pc(1,5),pc(1,6),pc(1,7),
     3   pc(1,8),pc(1,9),pc(1,10),pc(1,11),pc(1,12),pc(1,13),pc(1,14),
     4   pc(1,15),pc(1,16),pc(1,17),pc(1,18),pc(1,19),pc(1,20),pc(1,21),
     5   pc(1,22),pc(1,23),pc(1,24),pc(1,25),pc(1,26),pc(1,27))
c*
c*    *** return and end ***
      return
      end
	  
      subroutine restrc2(nxf,nyf,nzf,nxc,nyc,nzc,xin,xout,
     2   oPC,oPN,oPS,oPE,oPW,oPNE,oPNW,oPSE,oPSW,
     3   uPC,uPN,uPS,uPE,uPW,uPNE,uPNW,uPSE,uPSW,
     4   dPC,dPN,dPS,dPE,dPW,dPNE,dPNW,dPSE,dPSW)
c* *********************************************************************
c* purpose:
c*
c*    generic restriction of a fine grid function to coarse grid.
c*    this is the adjoint of some prolongation:
c*
c*         R = 2**dim P^T, 
c*
c*    where dim=3 is the number of spatial dimensions, and P is the
c*    prolongation operator.
c*
c*    with the inner-product weighed by 2**dim, where here dim=3,
c*    this operator is the adjoint of P.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          idimenshun
      parameter        (idimenshun = 3)
      integer          nxf,nyf,nzf,nxc,nyc,nzc,i,j,k,ii,jj,kk
      double precision dimfac
      double precision xin(nxf,nyf,nzf),xout(nxc,nyc,nzc)
      double precision oPC(nxc,nyc,nzc),oPN(nxc,nyc,nzc)
      double precision oPS(nxc,nyc,nzc),oPE(nxc,nyc,nzc)
      double precision oPW(nxc,nyc,nzc),oPNE(nxc,nyc,nzc)
      double precision oPNW(nxc,nyc,nzc),oPSE(nxc,nyc,nzc)
      double precision oPSW(nxc,nyc,nzc),uPC(nxc,nyc,nzc)
      double precision uPN(nxc,nyc,nzc),uPS(nxc,nyc,nzc)
      double precision uPE(nxc,nyc,nzc),uPW(nxc,nyc,nzc)
      double precision uPNE(nxc,nyc,nzc),uPNW(nxc,nyc,nzc)
      double precision uPSE(nxc,nyc,nzc),uPSW(nxc,nyc,nzc)
      double precision dPC(nxc,nyc,nzc),dPN(nxc,nyc,nzc)
      double precision dPS(nxc,nyc,nzc),dPE(nxc,nyc,nzc)
      double precision dPW(nxc,nyc,nzc),dPNE(nxc,nyc,nzc)
      double precision dPNW(nxc,nyc,nzc),dPSE(nxc,nyc,nzc)
      double precision dPSW(nxc,nyc,nzc)
      double precision tmpO,tmpU,tmpD
c*
cmdir 0 0
c*
c*    *** verify correctness of the input boundary points ***
      call fboundPMG00(nxf,nyf,nzf,xin)
c*
c*    *** determine dimension factor ***
      dimfac  = 2.**idimenshun
c*
c*    *** handle the interior points as average of 5 finer grid pts ***
!$OMP PARALLEL default(shared) private(k,kk,j,jj,i,ii,tmpO,tmpU,tmpD)
!$OMP DO
      do k = 2, nzc-1
         kk = (k - 1) * 2 + 1
         do j = 2, nyc-1
            jj = (j - 1) * 2 + 1
            do i = 2, nxc-1
               ii = (i - 1) * 2 + 1
c*
c*             *** compute the restriction ***
               tmpO =
     2            +  oPC(i,j,k)        * xin(ii,jj,kk)
     3            +  oPN(i,j,k)        * xin(ii,jj+1,kk)
     4            +  oPS(i,j,k)        * xin(ii,jj-1,kk)
     5            +  oPE(i,j,k)        * xin(ii+1,jj,kk)
     6            +  oPW(i,j,k)        * xin(ii-1,jj,kk)
     7            +  oPNE(i,j,k)       * xin(ii+1,jj+1,kk)
     8            +  oPNW(i,j,k)       * xin(ii-1,jj+1,kk)
     9            +  oPSE(i,j,k)       * xin(ii+1,jj-1,kk)
     9            +  oPSW(i,j,k)       * xin(ii-1,jj-1,kk)
               tmpU =
     2            +  uPC(i,j,k)        * xin(ii,jj,kk+1)
     3            +  uPN(i,j,k)        * xin(ii,jj+1,kk+1)
     4            +  uPS(i,j,k)        * xin(ii,jj-1,kk+1)
     5            +  uPE(i,j,k)        * xin(ii+1,jj,kk+1)
     6            +  uPW(i,j,k)        * xin(ii-1,jj,kk+1)
     7            +  uPNE(i,j,k)       * xin(ii+1,jj+1,kk+1)
     8            +  uPNW(i,j,k)       * xin(ii-1,jj+1,kk+1)
     9            +  uPSE(i,j,k)       * xin(ii+1,jj-1,kk+1)
     9            +  uPSW(i,j,k)       * xin(ii-1,jj-1,kk+1)
               tmpD =
     2            +  dPC(i,j,k)        * xin(ii,jj,kk-1)
     3            +  dPN(i,j,k)        * xin(ii,jj+1,kk-1)
     4            +  dPS(i,j,k)        * xin(ii,jj-1,kk-1)
     5            +  dPE(i,j,k)        * xin(ii+1,jj,kk-1)
     6            +  dPW(i,j,k)        * xin(ii-1,jj,kk-1)
     7            +  dPNE(i,j,k)       * xin(ii+1,jj+1,kk-1)
     8            +  dPNW(i,j,k)       * xin(ii-1,jj+1,kk-1)
     9            +  dPSE(i,j,k)       * xin(ii+1,jj-1,kk-1)
     9            +  dPSW(i,j,k)       * xin(ii-1,jj-1,kk-1)
               xout(i,j,k) = tmpO + tmpU + tmpD
		  end do
		end do
	   end do
!$OMP END DO
!$OMP END PARALLEL
c*
c*    *** verify correctness of the output boundary points ***
      call fboundPMG00(nxc,nyc,nzc,xout)
c*
c*    *** return and end ***
      return
      end
	  
      subroutine extrac(nxf,nyf,nzf,nxc,nyc,nzc,xin,xout)
c* *********************************************************************
c* purpose:
c*
c*    simple injection of a fine grid function into coarse grid.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          nxf,nyf,nzf,nxc,nyc,nzc,i,j,k,ii,jj,kk
      double precision xin(nxf,nyf,nzf),xout(nxc,nyc,nzc)
c*
cmdir 0 0
c*
c*    *** verify correctness of the input boundary points ***
      call fboundPMG00(nxf,nyf,nzf,xin)
c*
c*    *** doit ***
cmdir 3 1
      do 10 k = 2, nzc-1
         kk = (k - 1) * 2 + 1
cmdir 3 2
         do 11 j = 2, nyc-1
            jj = (j - 1) * 2 + 1
cmdir 3 3
            do 12 i = 2, nxc-1
               ii = (i - 1) * 2 + 1
c*
c*             *** compute the restriction ***
               xout(i,j,k) = xin(ii,jj,kk)
 12         continue
 11      continue
 10   continue
c*
c*    *** verify correctness of the output boundary points ***
      call fboundPMG00(nxc,nyc,nzc,xout)
c*
c*    *** return and end ***
      return
      end
      subroutine interpPMGZ(nxc,nyc,nzc,nxf,nyf,nzf,xin,xout)
c* *********************************************************************
c* purpose:
c*
c*    fast trilinear interpolation.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          nxc,nyc,nzc,nxf,nyf,nzf,i,j,k,jj,ii,kk
      double precision xin(nxc,nyc,nzc),xout(nxf,nyf,nzf)
      double precision c0,c1,c2,c3
c*
c* *********************************************************************
c* setup 
c* *********************************************************************
c*
cmdir 0 0
c*
c*    *** verify correctness of the input boundary points ***
      call fboundPMG00(nxc,nyc,nzc,xin)
c*
c*    *** set interpolation weights ***
      c0 = 1.0d0
      c1 = 0.5e0
      c2 = 0.25e0
      c3 = 0.125e0
c*
c*    *** doit ***
cmdir 3 1
      do 10 k = 1, nzf-2, 2
         kk = (k-1) / 2 + 1
cmdir 3 2
         do 11 j = 1, nyf-2, 2
            jj = (j-1) / 2 + 1
cmdir 3 3
            do 12 i = 1, nxf-2, 2
               ii = (i-1) / 2 + 1
c*
c* *********************************************************************
c* type 1 -- fine grid points common to a coarse grid point
c* *********************************************************************
c*
c*             *** copy coinciding points from coarse grid to fine grid ***
               xout(i,j,k) = xin(ii,jj,kk)
c*
c* *********************************************************************
c* type 2 -- fine grid points common to a coarse grid plane
c* *********************************************************************
c*
c*             *** fine grid pts common only to y-z planes on coarse grid ***
c*             *** (intermediate pts between 2 grid points on x-row) ***
               xout(i+1,j,k) = c1*(xin(ii,jj,kk) + xin(ii+1,jj,kk))
c*
c*             *** fine grid pts common only to x-z planes on coarse grid ***
c*             *** (intermediate pts between 2 grid points on a y-row) ***
               xout(i,j+1,k) = c1*(xin(ii,jj,kk) + xin(ii,jj+1,kk))
c*
c*             *** fine grid pts common only to x-y planes on coarse grid ***
c*             *** (intermediate pts between 2 grid points on a z-row) ***
               xout(i,j,k+1) = c1*(xin(ii,jj,kk) + xin(ii,jj,kk+1))
c*
c* *********************************************************************
c* type 3 -- fine grid points common to a coarse grid line
c* *********************************************************************
c*
c*             *** fine grid pts common only to z planes on coarse grid ***
c*             *** (intermediate pts between 4 grid pts on the xy-plane***
               xout(i+1,j+1,k) = c2*(xin(ii,jj,kk)  +xin(ii+1,jj,kk)
     2                              +xin(ii,jj+1,kk)+xin(ii+1,jj+1,kk))
c*
c*             *** fine grid pts common only to y planes on coarse grid ***
c*             *** (intermediate pts between 4 grid pts on the xz-plane***
               xout(i+1,j,k+1) = c2*(xin(ii,jj,kk)  +xin(ii+1,jj,kk)
     2                              +xin(ii,jj,kk+1)+xin(ii+1,jj,kk+1))
c*
c*             *** fine grid pts common only to x planes on coarse grid ***
c*             *** (intermediate pts between 4 grid pts on the yz-plane***
               xout(i,j+1,k+1) = c2*(xin(ii,jj,kk)  +xin(ii,jj+1,kk)
     2                              +xin(ii,jj,kk+1)+xin(ii,jj+1,kk+1))
c*
c* *********************************************************************
c* type 4 -- fine grid points not common to coarse grid pts/lines/planes
c* *********************************************************************
c*
c*             *** completely interior points ***
               xout(i+1,j+1,k+1) = c3*(
     2            + xin(ii,jj,kk)     + xin(ii+1,jj,kk)
     3            + xin(ii,jj+1,kk)   + xin(ii+1,jj+1,kk)
     4            + xin(ii,jj,kk+1)   + xin(ii+1,jj,kk+1)
     5            + xin(ii,jj+1,kk+1) + xin(ii+1,jj+1,kk+1) )
 12         continue
 11      continue
 10   continue
c*
c*    *** verify correctness of the output boundary points ***
      call fboundPMG00(nxf,nyf,nzf,xout)
c*
c*    *** return and end ***
      return
      end
      subroutine restrcZ(nxf,nyf,nzf,nxc,nyc,nzc,xin,xout)
c* *********************************************************************
c* purpose:
c*
c*    fast full-weighting restriction.
c*
c*    trilinear restriction of a fine grid function to coarse grid.
c*    this is the adjoint of trilinear prolongation:
c*
c*         R = 2**dim P^T, 
c*
c*    where dim=3 is the number of spatial dimensions, and P is the
c*    prolongation operator.
c*
c*    with the inner-product weighed by 2**dim, where here dim=3,
c*    this operator is the adjoint of P.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          idimenshun
      parameter        (idimenshun = 3)
      integer          nxf,nyf,nzf,nxc,nyc,nzc,i,j,k,ii,jj,kk
      double precision xin(nxf,nyf,nzf),xout(nxc,nyc,nzc)
      double precision c0,c1,c2,c3,dimfac
c*
cmdir 0 0
c*
c*    *** verify correctness of the input boundary points ***
      call fboundPMG00(nxf,nyf,nzf,xin)
c*
c*    *** determine dimension factor ***
      dimfac  = 2.**idimenshun
c*
c*    *** set restriction weights ***
      c0 = 1.0d0
      c1 = 0.5e0
      c2 = 0.25e0
      c3 = 0.125e0
c*
c*    *** handle the interior points as average of 5 finer grid pts ***
cmdir 3 1
      do 10 k = 2, nzc-1
         kk = (k - 1) * 2 + 1
cmdir 3 2
         do 11 j = 2, nyc-1
            jj = (j - 1) * 2 + 1
cmdir 3 3
            do 12 i = 2, nxc-1
               ii = (i - 1) * 2 + 1
c*
c*             *** compute the restriction ***
               xout(i,j,k) = c0 * xin(ii,jj,kk)
     2            + c1 * (xin(ii-1,jj,kk) + xin(ii+1,jj,kk)
     3                   +xin(ii,jj-1,kk) + xin(ii,jj+1,kk)   
     4                   +xin(ii,jj,kk-1) + xin(ii,jj,kk+1))
     5            + c2 * (xin(ii-1,jj-1,kk) + xin(ii+1,jj-1,kk)  
     6                   +xin(ii-1,jj+1,kk) + xin(ii+1,jj+1,kk)
     7                   +xin(ii-1,jj,kk-1) + xin(ii+1,jj,kk-1) 
     8                   +xin(ii-1,jj,kk+1) + xin(ii+1,jj,kk+1)
     9                   +xin(ii,jj-1,kk-1) + xin(ii,jj+1,kk-1) 
     9                   +xin(ii,jj-1,kk+1) + xin(ii,jj+1,kk+1))
     9            + c3 * (xin(ii-1,jj-1,kk-1) + xin(ii+1,jj-1,kk-1)
     9                   +xin(ii-1,jj+1,kk-1) + xin(ii+1,jj+1,kk-1)
     9                   +xin(ii-1,jj-1,kk+1) + xin(ii+1,jj-1,kk+1)
     9                   +xin(ii-1,jj+1,kk+1) + xin(ii+1,jj+1,kk+1))
 12       continue
 11     continue
 10   continue
c*
c*    *** verify correctness of the output boundary points ***
      call fboundPMG00(nxc,nyc,nzc,xout)
c*
c*    *** return and end ***
      return
      end

