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

      subroutine buildband(key,nx,ny,nz,ipc,rpc,ac,ipcB,rpcB,acB)
c* *********************************************************************
c* purpose:
c*
c*    build and factor a banded matrix given a matrix in diagonal form.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),ipcB(*),nx,ny,nz,numdia,n,m,lda,info,key
      double precision rpc(*),ac(nx*ny*nz,*),rpcB(*),acB(*)
CZZZZ double precision rcond
c*
cmdir 0 0
c*
c*    *** do in one step ***
      numdia = ipc(11)
      if (numdia .eq. 7) then
         print*,'% BUILDBAND: building 7-pt banded coarse operator...'
         n   = (nx-2)*(ny-2)*(nz-2)
         m   = (nx-2)*(ny-2)
         lda = m+1
         call buildband1_7(nx,ny,nz,ipc,rpc,ac(1,1),
     2      ac(1,2),ac(1,3),ac(1,4),
     3      ipcB,rpcB,acB,n,m,lda)
      elseif (numdia .eq. 27) then
         print*,'% BUILDBAND: building 27-pt banded coarse operator...'
         n   = (nx-2)*(ny-2)*(nz-2)
         m   = (nx-2)*(ny-2) + (nx-2) + 1
         lda = m+1
         call buildband1_27(nx,ny,nz,ipc,rpc,ac(1,1),
     2      ac(1,2),ac(1,3),ac(1,4),
     3      ac(1,5),ac(1,6),ac(1,7),ac(1,8),ac(1,9),
     4      ac(1,10),ac(1,11),ac(1,12),ac(1,13),ac(1,14),
     5      ipcB,rpcB,acB,n,m,lda)
      else
         print*,'% BUILDBAND: invalid stencil type given...'
      endif
c*
c*    *** print matrices ***
CZZZZ call prtmatd(nx,ny,nz,ipc,rpc,ac)
CZZZZ call prtmatb(acB,n,m,lda)
c*
c*    *** factor the system ***
      key  = 0
      info = 0
      print*,'% BUILDBAND: dpbfa factoring banded coarse system...'
      call dpbfa(acB,lda,n,m,info)
CZZZZ call dpbco(acB,lda,n,m,rcond,w1,info)
      ipcB(4) = 1
      if (info .ne. 0) then
CZZZZ    print*,'% BUILDBAND: dpbco rcond:   ',rcond
         print*,'% BUILDBAND: dpbfa problem: ',info
         print*,'% BUILDBAND: leading principle minor not PD...'
         key = 1
      endif
      print*,'% BUILDBAND: dpbfa banded factorization complete.'
c*
c*    *** return and end ***
      return
      end
      subroutine buildband1_7 (nx,ny,nz,ipc,rpc,
     2   oC,oE,oN,uC,
     3   ipcB,rpcB,acB,n,m,lda)
c* *********************************************************************
c* purpose: 
c*
c*    build the operator in banded form given the 7-diagonal form.
c*
c* notes:
c*
c*    pseudo-code from the banded linpack routines:
c*
c*    jj = 0
c*    do 10 k = 2, nz-1
c*       do 10 j = 2, ny-1
c*          do 10 i = 2, nx-1
c*             jj = jj + 1
c*             i1 = max0(1, jj-m)
c*             do 20 ii = i1, jj
c*                kk = ii - jj + m + 1
c*                acB(kk,jj) = 0.0e0
c*                if ((jj-ii).eq.0) acB(kk,jj)=oC(i,j,k)
c*                if ((jj-ii).eq.1) acB(kk,jj)=-oE(i-1,j,k)
c*                if ((jj-ii).eq.(nx-2)) acB(kk,jj)=-oN(i,j-1,k)
c*                if ((jj-ii).eq.(nx-2)*(ny-2)) acB(kk,jj)=-uC(i,j,k-1)
c* 20          continue
c* 10   continue
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),ipcB(*),nx,ny,nz,lda,m,n
      integer          i,j,k,ii,jj,kk
      double precision rpc(*),rpcB(*),acB(lda,*),oC(nx,ny,nz)
      double precision oE(nx,ny,nz),oN(nx,ny,nz),uC(nx,ny,nz)
c*
cmdir 0 0
c*
c*    *** doit ***
      ipcB(1)  = n
      ipcB(2)  = m
      ipcB(3)  = lda
      ipcB(4)  = 0
      jj = 0
cmdir 3 1
      do 10 k = 2, nz-1
cmdir 3 2
         do 11 j = 2, ny-1
cmdir 3 3
            do 12 i = 2, nx-1
               jj = jj + 1
c*
c*             *** diagonal term ***
               ii = jj
               kk = ii - jj + m + 1
               acB(kk,jj) = oC(i,j,k)
c*
c*             *** east neighbor ***
               ii = jj-1
               kk = ii - jj + m + 1
               acB(kk,jj) = - oE(i-1,j,k)
c*
c*             *** north neighbor ***
               ii = jj-(nx-2)
               kk = ii - jj + m + 1
               acB(kk,jj) = - oN(i,j-1,k)
c*
c*             *** up neighbor ***
               ii = jj-(nx-2)*(ny-2)
               kk = ii - jj + m + 1
               acB(kk,jj) = - uC(i,j,k-1)
 12         continue
 11      continue
 10   continue
c*
c*    *** return and end ***
      return
      end
      subroutine buildband1_27 (nx,ny,nz,ipc,rpc,
     2   oC,oE,oN,uC,
     3   oNE,oNW,uE,uW,uN,uS,uNE,uNW,uSE,uSW,
     4   ipcB,rpcB,acB,n,m,lda)
c* *********************************************************************
c* purpose: 
c*
c*    build the operator in banded form given the 27-diagonal form.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),ipcB(*),nx,ny,nz,lda,m,n
      integer          i,j,k,ii,jj,kk
      double precision rpc(*),rpcB(*),acB(lda,*),oC(nx,ny,nz)
      double precision oE(nx,ny,nz),oN(nx,ny,nz),uC(nx,ny,nz)
      double precision oNE(nx,ny,nz),oNW(nx,ny,nz),uE(nx,ny,nz)
      double precision uW(nx,ny,nz),uN(nx,ny,nz),uS(nx,ny,nz)
      double precision uNE(nx,ny,nz),uNW(nx,ny,nz),uSE(nx,ny,nz)
      double precision uSW(nx,ny,nz)
c*
cmdir 0 0
c*
c*    *** doit ***
      ipcB(1)  = n
      ipcB(2)  = m
      ipcB(3)  = lda
      ipcB(4)  = 0
      jj = 0
cmdir 3 1
      do 10 k = 2, nz-1
cmdir 3 2
         do 11 j = 2, ny-1
cmdir 3 3
            do 12 i = 2, nx-1
               jj = jj + 1
c*
c*             *** diagonal term ***
               ii = jj
               kk = ii - jj + m + 1
               acB(kk,jj) = oC(i,j,k)
c*
c*             *** east neighbor ***
               ii = jj-1
               kk = ii - jj + m + 1
               acB(kk,jj) = - oE(i-1,j,k)
c*
c*             *** north neighbor ***
               ii = jj-(nx-2)
               kk = ii - jj + m + 1
               acB(kk,jj) = - oN(i,j-1,k)
c*
c*             *** north-east neighbor ***
               ii = jj-(nx-2)+1
               kk = ii - jj + m + 1
               acB(kk,jj) = - oNE(i,j-1,k)
c*
c*             *** north-west neighbor ***
               ii = jj-(nx-2)-1
               kk = ii - jj + m + 1
               acB(kk,jj) = - oNW(i,j-1,k)
c*
c*             *** up neighbor ***
               ii = jj-(nx-2)*(ny-2)
               kk = ii - jj + m + 1
               acB(kk,jj) = - uC(i,j,k-1)
c*
c*             *** up-east neighbor ***
               ii = jj-(nx-2)*(ny-2)+1
               kk = ii - jj + m + 1
               acB(kk,jj) = - uE(i,j,k-1)
c*
c*             *** up-west neighbor ***
               ii = jj-(nx-2)*(ny-2)-1
               kk = ii - jj + m + 1
               acB(kk,jj) = - uW(i,j,k-1)
c*
c*             *** up-north neighbor ***
               ii = jj-(nx-2)*(ny-2)+(nx-2)
               kk = ii - jj + m + 1
               acB(kk,jj) = - uN(i,j,k-1)
c*
c*             *** up-south neighbor ***
               ii = jj-(nx-2)*(ny-2)-(nx-2)
               kk = ii - jj + m + 1
               acB(kk,jj) = - uS(i,j,k-1)
c*
c*             *** up-north-east neighbor ***
               ii = jj-(nx-2)*(ny-2)+(nx-2)+1
               kk = ii - jj + m + 1
               acB(kk,jj) = - uNE(i,j,k-1)
c*
c*             *** up-north-west neighbor ***
               ii = jj-(nx-2)*(ny-2)+(nx-2)-1
               kk = ii - jj + m + 1
               acB(kk,jj) = - uNW(i,j,k-1)
c*
c*             *** up-south-east neighbor ***
               ii = jj-(nx-2)*(ny-2)-(nx-2)+1
               kk = ii - jj + m + 1
               acB(kk,jj) = - uSE(i,j,k-1)
c*
c*             *** up-south-west neighbor ***
               ii = jj-(nx-2)*(ny-2)-(nx-2)-1
               kk = ii - jj + m + 1
               acB(kk,jj) = - uSW(i,j,k-1)
 12         continue
 11      continue
 10   continue
c*
c*    *** return and end ***
      return
      end

