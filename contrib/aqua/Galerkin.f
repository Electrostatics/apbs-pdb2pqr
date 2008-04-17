c ==========================================================================================
c Aqua 
c Copyright (C) 2007 Patrice Koehl
c
c This library is free software; you can redistribute it and/or
c modify it under the terms of the GNU Lesser General Public
c License as published by the Free Software Foundation; either
c version 2.1 of the License, or (at your option) any later version.
c
c This library is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
c Lesser General Public License for more details.
c
c You should have received a copy of the GNU Lesser General Public
c License along with this library; if not, write to the Free Software
c Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
c
c ==========================================================================================
c Works in this file have been derived or copied from PMG. The original copyright,
c license and warranty are provided below.
c
c PMG -- Parallel algebraic MultiGrid
c Copyright (c) 1994-2008.  Michael Holst.
c 
c Michael Holst <mholst@math.ucsd.edu>
c University of California, San Diego
c Department of Mathematics, 5739 AP&M
c 9500 Gilman Drive, Dept. 0112
c La Jolla, CA 92093-0112 USA                                                 
c http://math.ucsd.edu/~mholst
c 
c This library is free software; you can redistribute it and/or
c modify it under the terms of the GNU Lesser General Public
c License as published by the Free Software Foundation; either
c version 2.1 of the License, or (at your option) any later version.
c 
c This library is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
c Lesser General Public License for more details.
c
c You should have received a copy of the GNU Lesser General Public
c License along with this library; if not, write to the Free Software
c Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
c ==========================================================================================

c ==========================================================================================
c	Galerkin.f
c       (Galerkin operator, to restrict Laplacian operator for fine grid to coarse grid)
c
c	Version 1:      7/4/07
c
c	Author: Patrice Koehl (in collaboration with Marc Delarue)
c
c	This file contains a series of subroutines for solving modified Poisson-Boltzmann
c	equations. It is inspired form Michael Holst MG code, which is available at:
c	http://www.fetk.org/codes/pmg/index.html
c
c
c ==========================================================================================
c ==========================================================================================
c	Build_Galerkin.f
c ==========================================================================================
c ==========================================================================================
c
c	Form the Galerkin coarse grid system
c
c	(note that the coarse grid is always 27 diagonal (with only 14 stored, as
c	the matrix is symmetric)
c
	subroutine build_Galerkin(nxf,nyf,nzf,nxc,nyc,nzc,ipkey,numdia,
     1			ipcF,rpcF,acF,ccF,fcF,ipc,rpc,ac,cc,fc)
c
c
c ==========================================================================================
c	Declare all variables
c ==========================================================================================
c
c	Input/Output
c ==========================================================================================
c
	integer	numdia,numdia_local
	integer	nxf,nyf,nzf,nxc,nyc,nzc,ipkey
c
	integer	ipcF(*),ipc(*)
c
	real*8	rpcF(*),acF(*),ccF(*),fcF(*)
	real*8	rpc(*),ac(*),cc(*),fc(*)
c
c ==========================================================================================
c	Call Galerkin routine
c ==========================================================================================
c
	numdia_local = ipcF(11)
	call build_G(nxf,nyf,nzf,nxc,nyc,nzc,numdia_local,acF,ac)
c
	ipc(11)	= 27
	numdia  = 14
c
	ipc(10) = ipkey
c
	call restrict(nxf,nyf,nzf,nxc,nyc,nzc,ccF,cc)
	call restrict(nxf,nyf,nzf,nxc,nyc,nzc,fcF,fc)
c
	return
	end
c
c ==========================================================================================
c ==========================================================================================
c	Build_G.f
c ==========================================================================================
c ==========================================================================================
c
c	Form the Galerkin coarse grid system from the fine grid and the prolongation
c	operator
c
c	Only consider 7 diagonal fine grid matrices
c
	subroutine build_G(nxf,nyf,nzf,nxc,nyc,nzc,numdia,acF,ac)
c
c ==========================================================================================
c	Declare all variables
c ==========================================================================================
c
c	Input/Output
c ==========================================================================================
c
	integer	numdia
	integer	nxf,nyf,nzf,nxc,nyc,nzc
c
	real*8	acF(nxf*nyf*nzf,*),ac(nxc*nyc*nzc,*)
c
c ==========================================================================================
c	Call the build routine
c ==========================================================================================
c
	if(numdia.eq.7) then
c
		call build_G7(nxf,nyf,nzf,nxc,nyc,nzc,
     7		acF(1,1),acF(1,2),acF(1,3),acF(1,4),
     8		ac(1,1),ac(1,2),ac(1,3),ac(1,4),ac(1,5),ac(1,6),
     9		ac(1,7),ac(1,8),ac(1,9),ac(1,10),ac(1,11),ac(1,12),
     1		ac(1,13),ac(1,14))
c
	elseif(numdia.eq.27) then
c
		call build_G27(nxf,nyf,nzf,nxc,nyc,nzc,
     7		acF(1,1),acF(1,2),acF(1,3),acF(1,4),
     8		acF(1,5),acF(1,6),acF(1,7),acF(1,8),acF(1,9),
     9		acF(1,10),acF(1,11),acF(1,12),acF(1,13),acF(1,14),
     8		ac(1,1),ac(1,2),ac(1,3),ac(1,4),ac(1,5),ac(1,6),
     9		ac(1,7),ac(1,8),ac(1,9),ac(1,10),ac(1,11),ac(1,12),
     1		ac(1,13),ac(1,14))
c
	else
c
		write(6,*) 'Galerkin projections:'
		write(6,*) 'Only 7diag and 27diag matrices implemented'
		stop
c
	endif
c
	return
	end
c
c ==========================================================================================
c ==========================================================================================
c	Build_G7.f
c ==========================================================================================
c ==========================================================================================
c
c	Compute a 27-point Galerkin coarse grid matrix from a 7-point fine grid matrix
c
c
	subroutine build_G7(nxf,nyf,nzf,nx,ny,nz,
     4		oC,oE,oN,uC,
     5		XoC,XoE,XoN,XuC,XoNE,XoNW,XuE,XuW,XuN,XuS,XuNE,XuNW,
     6		XuSE,XuSW)
c
c ==========================================================================================
c	Declare all variables
c ==========================================================================================
c
c	Input/Output
c ==========================================================================================
c
	integer	nxf,nyf,nzf,nx,ny,nz
c
	real*8	oC(nxf,nyf,nzf),oE(nxf,nyf,nzf)
	real*8	oN(nxf,nyf,nzf),uC(nxf,nyf,nzf)
c
	real*8	XoC(nx,ny,nz),XoE(nx,ny,nz),XoN(nx,ny,nz)
	real*8	XuC(nx,ny,nz),XuE(nx,ny,nz),XuN(nx,ny,nz)
	real*8	XoNE(nx,ny,nz),XoNW(nx,ny,nz),XuW(nx,ny,nz)
	real*8	XuS(nx,ny,nz),XuNE(nx,ny,nz),XuNW(nx,ny,nz)
	real*8	XuSE(nx,ny,nz),XuSW(nx,ny,nz)
c
c	Local variables
c ==========================================================================================
c
	integer	i,j,k,ii,jj,kk
	integer	im1,ip1,im2,ip2,jm1,jp1,jm2,jp2,km1,kp1,km2,kp2
	integer	nxm1,nym1,nzm1
c
	real*8	whole,half,quarter,eighth
c
	real*8	tmp1_XoC,tmp2_XoC,tmp3_XoC,tmp4_XoC
	real*8	tmp5_XoC,tmp6_XoC,tmp7_XoC,tmp8_XoC
	real*8	tmp9_XoC
	real*8	tmp1_XoE,tmp2_XoE,tmp3_XoE,tmp4_XoE
	real*8	tmp1_XoN,tmp2_XoN,tmp3_XoN,tmp4_XoN
	real*8	tmp1_XuC,tmp2_XuC,tmp3_XuC,tmp4_XuC
c
c ==========================================================================================
c	Build the operator
c ==========================================================================================
c
c	Setup
c
	whole = 1.d0
	half = whole/2.d0
	quarter = half/2.d0
	eighth = quarter/2.d0
c
	nxm1 = nx - 1
	nym1 = ny - 1
	nzm1 = nz - 1
c
	do 300 kk = 2,nzm1
	  k = 2*kk-1
c
	  do 200 jj =2,nym1
	    j = 2*jj-1
c
	    do 100 ii = 2,nxm1
		i = 2*ii - 1
c
c		Set indicies
c              =============================================================================
c
		im1 = i-1
		ip1 = i+1
		im2 = i-2
		ip2 = i+2
		jm1 = j-1
		jp1 = j+1
		jm2 = j-2
		jp2 = j+2
		km1 = k-1
		kp1 = k+1
		km2 = k-2
		kp2 = k+2
c
c		XoC(ii,jj,kk)
c		=============================================================================
c
		tmp1_XoC =
     1			eighth * (oC(im1,jm1,km1)*eighth
     2					- uC(im1,jm1,km1)*quarter
     3					- oN(im1,jm1,km1)*quarter
     4					- oE(im1,jm1,km1)*quarter)
c
     5		+	quarter * (-oE(i,jp1,k)*half
     6					- oN(ip1,j,k)*half
     7					- uC(ip1,jp1,km1)*eighth
     8					+ oC(ip1,jp1,k)*quarter
     9					- uC(ip1,jp1,k)*eighth)
c
     5		+	eighth * (-oE(i,jp1,km1)*quarter
     6					- oN(ip1,j,km1)*quarter
     8					+ oC(ip1,jp1,km1)*eighth
     7					- uC(ip1,jp1,km1)*quarter)
c
c
		tmp2_XoC =
     1			eighth * (-oE(i,jm1,km1)*quarter
     2					+ oC(ip1,jm1,km1)*eighth
     3					- uC(ip1,jm1,km1)*quarter
     4					- oN(ip1,jm1,km1)*quarter)
c
     5		+	quarter * (-oE(i,j,kp1)*half
     6					- oN(ip1,jm1,kp1)*eighth
     7					- uC(ip1,j,k)*half
     8					+ oC(ip1,j,kp1)*quarter
     9					- oN(ip1,j,kp1)*eighth)
c
     5		+	half * (-oE(i,j,k)*whole
     6					- oN(ip1,jm1,k)*quarter
     8					- uC(ip1,j,km1)*quarter
     7					+ oC(ip1,j,k)*half
     8					- uC(ip1,j,k)*quarter
     7					- oN(ip1,j,k)*quarter)
c
c
		tmp3_XoC =
     2			+ quarter * (- oE(i,j,km1)*half
     2					- oN(ip1,jm1,km1)*eighth
     2					+ oC(ip1,j,km1)*quarter
     2					- uC(ip1,j,km1)*half
     2					- oN(ip1,j,km1)*eighth)
c
     2			+ eighth*(-oE(i,jm1,kp1)*quarter
     2					-uC(ip1,jm1,k) * quarter
     2					+ oC(ip1,jm1,kp1)*eighth
     2					-oN(ip1,jm1,kp1)*quarter)

     2			+ eighth*(-oE(i,jp1,kp1)*quarter
     2					- oN(ip1,j,kp1) * quarter
     2					- uC(ip1,jp1,k) * quarter
     2					+oC(ip1,jp1,kp1)*eighth)
c
c
		tmp4_XoC =
     2			+ half*(-oE(im1,jm1,k)*quarter
     2					- uC(i,jm1,km1)*quarter
     2					+ oC(i,jm1,k) * half
     2					- uC(i,jm1,k) * quarter
     2					- oN(i,jm1,k) * whole
     2					- oE(i,jm1,k) * quarter)
c
     2			+quarter*(-oE(im1,jm1,km1)*eighth
     2					+ oC(i,jm1,km1) * quarter
     2					- uC(i,jm1,km1) * half
     2					- oN(i,jm1,km1) * half
     2					- oE(i,jm1,km1) *eighth)
c
     2			+ quarter*(-oE(i,jm1,k)*half
     2					- uC(ip1,jm1,km1)*eighth
     2					+ oC(ip1,jm1,k) * quarter
     2					- uC(ip1,jm1,k) * eighth
     2					- oN(ip1,jm1,k) * half)
c
c
		tmp5_XoC =
     2			+ quarter*(-oE(im1,jp1,km1)*eighth
     2					- oN(i,j,km1) * half
     2					+ oC(i,jp1,km1) * quarter
     2					- uC(i,jp1,km1) * half
     2					- oE(i,jp1,km1)*eighth)

     2			+ half*(-oE(im1,j,kp1)*quarter
     2					- oN(i,jm1,kp1) * quarter
     2					- uC(i,j,k) * whole
     2					+ oC(i,j,kp1) * half
     2					- oN(i,j,kp1) * quarter
     2					- oE(i,j,kp1) * quarter)

     2			+ whole * (- oE(im1,j,k) * half
     2					- oN(i,jm1,k) * half
     2					- uC(i,j,km1) * half
     2					+ oC(i,j,k) * whole
     2					- uC(i,j,k) * half
     2					- oN(i,j,k) * half
     2					- oE(i,j,k) * half)
c
c
		tmp6_XoC =
     2			+ half*(-oE(im1,j,km1)*quarter
     2					- oN(i,jm1,km1) * quarter
     2					+ oC(i,j,km1) * half
     2					- uC(i,j,km1) * whole
     2					- oN(i,j,km1) * quarter
     2					- oE(i,j,km1) * quarter)
c
     2   		+ quarter*(-oE(im1,jm1,kp1)*eighth
     2					- uC(i,jm1,k) * half
     2					+ oC(i,jm1,kp1)*quarter
     2					- oN(i,jm1,kp1)*half
     2					- oE(i,jm1,kp1)*eighth)
c
     2			+ eighth*(-oN(im1,j,kp1)*quarter
     2					- uC(im1,jp1,k) * quarter
     2					+ oC(im1,jp1,kp1)*eighth
     2					- oE(im1,jp1,kp1)*quarter)
c
c
		tmp7_XoC =
     2			+ quarter*(-oN(im1,jm1,km1)*eighth
     2					+ oC(im1,j,km1) * quarter
     2					- uC(im1,j,km1) * half
     2					- oN(im1,j,km1) * eighth
     2					- oE(im1,j,km1) * half)

     2			+eighth*(-uC(im1,jm1,k)*quarter
     2					+ oC(im1,jm1,kp1)*eighth
     2					- oN(im1,jm1,kp1)*quarter
     2					- oE(im1,jm1,kp1)*quarter)
c
     2			+quarter*(-uC(im1,jm1,km1)*eighth
     2					+ oC(im1,jm1,k) * quarter
     2					- uC(im1,jm1,k) * eighth
     2					- oN(im1,jm1,k) * half
     2					- oE(im1,jm1,k) * half)
c
c
		tmp8_XoC =
     2			+ quarter*(-oN(im1,j,k)*half
     2					- uC(im1,jp1,km1)*eighth
     2					+ oC(im1,jp1,k) * quarter
     2					- uC(im1,jp1,k) * eighth
     2					- oE(im1,jp1,k) * half)
c
     2			+ eighth*(-oN(im1,j,km1)*quarter
     2					+ oC(im1,jp1,km1)*eighth
     2					- uC(im1,jp1,km1)*quarter
     2					- oE(im1,jp1,km1)*quarter)

     2			+ half*(-oN(im1,jm1,k)*quarter
     2					- uC(im1,j,km1) * quarter
     2					+ oC(im1,j,k) * half
     2					- uC(im1,j,k) * quarter
     2					- oN(im1,j,k) * quarter
     2					- oE(im1,j,k) * whole)
c
c
		tmp9_XoC =
     2			+ quarter*(-oN(im1,jm1,kp1)*eighth
     2					- uC(im1,j,k) * half
     2					+ oC(im1,j,kp1) * quarter
     2					- oN(im1,j,kp1) * eighth
     2					- oE(im1,j,kp1) * half)
c
     2			+ quarter*(-oE(im1,jp1,kp1)*eighth
     2					- oN(i,j,kp1) * half
     2					- uC(i,jp1,k) * half
     2					+ oC(i,jp1,kp1) * quarter
     2					- oE(i,jp1,kp1) * eighth)
c
     2			+ half*(-oE(im1,jp1,k)*quarter
     2					- oN(i,j,k) * whole
     2					- uC(i,jp1,km1) * quarter
     2					+ oC(i,jp1,k) * half
     2					- uC(i,jp1,k) * quarter
     2					- oE(i,jp1,k) * quarter)
c
		XoC(ii,jj,kk) = tmp1_XoC + tmp2_XoC + tmp3_XoC +tmp4_XoC
     2		+ tmp5_XoC + tmp6_XoC + tmp7_XoC + tmp8_XoC + tmp9_XoC
c
c
c		XoE(ii,jj,kk)
c		=============================================================================
c
		tmp1_XoE =
     2		+ quarter * oE(i,jm1,km1) * eighth
     2			+ half * oE(i,jm1,k) * quarter
     2			+ quarter * oE(i,jm1,kp1)*eighth
     2			+ half * oE(i,j,km1) * quarter
     2			+ whole * oE(i,j,k) * half
     2			+ half * oE(i,j,kp1) * quarter
     2			+ quarter * oE(i,jp1,km1)*eighth
     2			+ half * oE(i,jp1,k) * quarter
     2			+ quarter * oE(i,jp1,kp1)*eighth
c
     2		- eighth * (  oC(ip1,jm1,km1) * eighth
     2				- uC(ip1,jm1,km1) * quarter
     2				- oN(ip1,jm1,km1) * quarter
     2				- oE(ip1,jm1,km1) * quarter)
c
c
      		tmp2_XoE =
     2		- quarter*(-uC(ip1,jm1,km1)*eighth
     2				+ oC(ip1,jm1,k) * quarter
     2				- uC(ip1,jm1,k) * eighth
     2				- oN(ip1,jm1,k) * half
     2				- oE(ip1,jm1,k) * half)
c
     2		- eighth*(- uC(ip1,jm1,k) * quarter
     2				+ oC(ip1,jm1,kp1) * eighth
     2				- oN(ip1,jm1,kp1) * quarter
     2				- oE(ip1,jm1,kp1) * quarter)
c
     2		- quarter * (- oN(ip1,jm1,km1) * eighth
     2				+ oC(ip1,j,km1) * quarter
     2				- uC(ip1,j,km1) * half
     2				- oN(ip1,j,km1) * eighth
     2				- oE(ip1,j,km1) * half)
c
c
		tmp3_XoE =
     2		- half * (- oN(ip1,jm1,k) * quarter
     2				- uC(ip1,j,km1) * quarter
     2				+ oC(ip1,j,k) * half
     2				- uC(ip1,j,k) * quarter
     2				- oN(ip1,j,k) * quarter
     2				- oE(ip1,j,k) * whole)
c
     2		- quarter * (- oN(ip1,jm1,kp1) * eighth
     2				- uC(ip1,j,k) * half
     2				+ oC(ip1,j,kp1) * quarter
     2				- oN(ip1,j,kp1) * eighth
     2				- oE(ip1,j,kp1) * half)
c
     2		- eighth * (- oN(ip1,j,km1) * quarter
     2				+ oC(ip1,jp1,km1) * eighth
     2				- uC(ip1,jp1,km1) * quarter
     2				- oE(ip1,jp1,km1) * quarter)
c
c
		tmp4_XoE =
     2		- quarter * (- oN(ip1,j,k) * half
     2				- uC(ip1,jp1,km1) * eighth
     2				+ oC(ip1,jp1,k) * quarter
     2				- uC(ip1,jp1,k) * eighth
     2				- oE(ip1,jp1,k) * half)
c
     2		- eighth * (- oN(ip1,j,kp1) * quarter
     2				- uC(ip1,jp1,k) * quarter
     2				+ oC(ip1,jp1,kp1) * eighth
     2				- oE(ip1,jp1,kp1) * quarter)
c
		XoE(ii,jj,kk) = tmp1_XoE+tmp2_XoE+tmp3_XoE+tmp4_XoE
		if(ii.eq.nxm1) XoE(ii,jj,kk) = 0.d0
c
c
c		XoN(ii,jj,kk)
c		=============================================================================
c
		tmp1_XoN =
     2		+ quarter * oN(im1,j,km1) * eighth
     2		+ half * oN(im1,j,k) * quarter
     2		+ quarter * oN(im1,j,kp1) * eighth
c
     2		- eighth*(oC(im1,jp1,km1)* eighth
     2				- uC(im1,jp1,km1) * quarter
     2				- oN(im1,jp1,km1) * quarter
     2				- oE(im1,jp1,km1) * quarter)

     2		- quarter * (- uC(im1,jp1,km1)*eighth
     2				+ oC(im1,jp1,k) * quarter
     2				- uC(im1,jp1,k) * eighth
     2				- oN(im1,jp1,k) * half
     2				- oE(im1,jp1,k) * half)
c
c
		tmp2_XoN =
     2		- eighth * (- uC(im1,jp1,k) * quarter
     2				+ oC(im1,jp1,kp1) * eighth
     2				- oN(im1,jp1,kp1) * quarter
     2				- oE(im1,jp1,kp1) * quarter)
c
     2		+ half * oN(i,j,km1) * quarter
     2		+ whole * oN(i,j,k) * half
     2		+ half * oN(i,j,kp1) * quarter
c
     2		- quarter * (- oE(im1,jp1,km1) * eighth
     2				+ oC(i,jp1,km1) * quarter
     2				- uC(i,jp1,km1) * half
     2				- oN(i,jp1,km1) * half
     2				- oE(i,jp1,km1) * eighth)
c
c
		tmp3_XoN =
     2		- half * (- oE(im1,jp1,k) * quarter
     2				- uC(i,jp1,km1) * quarter
     2				+ oC(i,jp1,k) * half
     2				- uC(i,jp1,k) * quarter
     2				- oN(i,jp1,k) * whole
     2				- oE(i,jp1,k) * quarter)
c
     2		- quarter * (- oE(im1,jp1,kp1) * eighth
     2				- uC(i,jp1,k) * half
     2				+ oC(i,jp1,kp1) * quarter
     2				- oN(i,jp1,kp1) * half
     2				- oE(i,jp1,kp1) * eighth)
c
     2		+ quarter * oN(ip1,j,km1) * eighth
     2		+ half * oN(ip1,j,k) * quarter
     2		+ quarter * oN(ip1,j,kp1) * eighth
c
c
		tmp4_XoN =
     2		- eighth * (- oE(i,jp1,km1) * quarter
     2				+ oC(ip1,jp1,km1) * eighth
     2				- uC(ip1,jp1,km1) * quarter
     2				- oN(ip1,jp1,km1) * quarter)
c
     2		- quarter * (- oE(i,jp1,k) * half
     2				- uC(ip1,jp1,km1) * eighth
     2				+ oC(ip1,jp1,k) * quarter
     2				- uC(ip1,jp1,k) * eighth
     2				- oN(ip1,jp1,k) * half)
c
     2		- eighth * (- oE(i,jp1,kp1) * quarter
     2				- uC(ip1,jp1,k) * quarter
     2				+ oC(ip1,jp1,kp1) * eighth
     2				- oN(ip1,jp1,kp1) * quarter)
c
		XoN(ii,jj,kk) = tmp1_XoN + tmp2_XoN + tmp3_XoN+tmp4_XoN
		if(jj.eq.nym1) XoN(ii,jj,kk) = 0.d0
c
c
c		XuC(ii,jj,kk)
c		=============================================================================
c
		tmp1_XuC =
     2		+ quarter * uC(im1,jm1,k) * eighth
c
     2		- eighth * (  oC(im1,jm1,kp1) * eighth
     2				- uC(im1,jm1,kp1) * quarter
     2				- oN(im1,jm1,kp1) * quarter
     2				- oE(im1,jm1,kp1) * quarter)
c
     2		+ half * uC(im1,j,k) * quarter
c
     2		- quarter * (- oN(im1,jm1,kp1) * eighth
     2				+ oC(im1,j,kp1) * quarter
     2				- uC(im1,j,kp1) * half
     2				- oN(im1,j,kp1) * eighth
     2				- oE(im1,j,kp1) * half)
c
     2		+ quarter * uC(im1,jp1,k) * eighth
c
c
		tmp2_XuC =
     2		- eighth * (- oN(im1,j,kp1) * quarter
     2				+ oC(im1,jp1,kp1) * eighth
     2				- uC(im1,jp1,kp1) * quarter
     2				- oE(im1,jp1,kp1) * quarter)
c
     2		+ half * uC(i,jm1,k) * quarter
c
     2		- quarter * (- oE(im1,jm1,kp1) * eighth
     2				+ oC(i,jm1,kp1) * quarter
     2				- uC(i,jm1,kp1) * half
     2				- oN(i,jm1,kp1) * half
     2				- oE(i,jm1,kp1) * eighth)
c
     2		+ whole * uC(i,j,k) * half
c
     2		- half * (- oE(im1,j,kp1) * quarter
     2				- oN(i,jm1,kp1) * quarter
     2				+ oC(i,j,kp1) * half
     2				- uC(i,j,kp1) * whole
     2				- oN(i,j,kp1) * quarter
     2				- oE(i,j,kp1) * quarter)
c
c
		tmp3_XuC =
     2		+ half * uC(i,jp1,k) * quarter
c
     2		- quarter * (- oE(im1,jp1,kp1) * eighth
     2				- oN(i,j,kp1) * half
     2				+ oC(i,jp1,kp1) * quarter
     2				- uC(i,jp1,kp1) * half
     2				- oE(i,jp1,kp1) * eighth)
c
     2		+ quarter * uC(ip1,jm1,k) * eighth
c
     2		- eighth * (- oE(i,jm1,kp1) * quarter
     2				+ oC(ip1,jm1,kp1) * eighth
     2				- uC(ip1,jm1,kp1) * quarter
     2				- oN(ip1,jm1,kp1) * quarter)
c
c
		tmp4_XuC =
     2		+ half * uC(ip1,j,k) * quarter
c
     2		- quarter * (- oE(i,j,kp1) * half
     2				- oN(ip1,jm1,kp1) * eighth
     2				+ oC(ip1,j,kp1) * quarter
     2				- uC(ip1,j,kp1) * half
     2				- oN(ip1,j,kp1) * eighth)
c
     2		+ quarter * uC(ip1,jp1,k) * eighth
c
     2		- eighth * (- oE(i,jp1,kp1) * quarter
     2				- oN(ip1,j,kp1) * quarter
     2				+ oC(ip1,jp1,kp1) * eighth
     2				- uC(ip1,jp1,kp1) * quarter)
c
		XuC(ii,jj,kk) = tmp1_XuC + tmp2_XuC + tmp3_XuC+tmp4_XuC
c
		if(kk.eq.nzm1) XuC(ii,jj,kk) = 0.d0
c
c
c		XoNE(ii,jj,kk)
c		=============================================================================
c
		XoNE(ii,jj,kk) =
     2		+ quarter * oE(i,jp1,km1) * eighth
     2		+ half * oE(i,jp1,k) * quarter
     2		+ quarter * oE(i,jp1,kp1) * eighth
     2		+ quarter * oN(ip1,j,km1) * eighth
     2		+ half * oN(ip1,j,k) * quarter
     2		+ quarter * oN(ip1,j,kp1) * eighth
c
     2		- eighth * (  oC(ip1,jp1,km1)*eighth
     2				- uC(ip1,jp1,km1) * quarter
     2				- oN(ip1,jp1,km1) * quarter
     2				- oE(ip1,jp1,km1) * quarter)
c
     2		- quarter * (- uC(ip1,jp1,km1)*eighth
     2				+ oC(ip1,jp1,k) * quarter
     2				- uC(ip1,jp1,k) * eighth
     2				- oN(ip1,jp1,k) * half
     2				- oE(ip1,jp1,k) * half)

     2		- eighth * (- uC(ip1,jp1,k) * quarter
     2				+ oC(ip1,jp1,kp1) * eighth
     2				- oN(ip1,jp1,kp1) * quarter
     2				- oE(ip1,jp1,kp1) * quarter)
c
		if(ii.eq.nxm1.or.jj.eq.nym1) XoNE(ii,jj,kk) = 0.d0
c
c
c		XoNW(ii,jj,kk)
c		=============================================================================
c
		XoNW(ii,jj,kk) =
     2		+ quarter * oN(im1,j,km1) * eighth
     2		+ half * oN(im1,j,k) * quarter
     2		+ quarter * oN(im1,j,kp1) * eighth
c
     2		- eighth * (- oE(im2,jp1,km1) *quarter
     2				+ oC(im1,jp1,km1) * eighth
     2				- uC(im1,jp1,km1) * quarter
     2				- oN(im1,jp1,km1) * quarter)
c
     2		- quarter * (- oE(im2,jp1,k) * half
     2				- uC(im1,jp1,km1) * eighth
     2				+ oC(im1,jp1,k) * quarter
     2				- uC(im1,jp1,k) * eighth
     2				- oN(im1,jp1,k) * half)
c
     2		- eighth * (- oE(im2,jp1,kp1) *quarter
     2				- uC(im1,jp1,k) * quarter
     2				+ oC(im1,jp1,kp1) * eighth
     2				- oN(im1,jp1,kp1) * quarter)
c
     2		+ quarter * oE(im1,jp1,km1) * eighth
     2		+ half * oE(im1,jp1,k) * quarter
     2		+ quarter * oE(im1,jp1,kp1) * eighth
c
		if(ii.eq.2) XoNW(ii,jj,kk) = 0
		if(jj.eq.nym1) XoNW(ii,jj,kk) = 0.d0
c
c
c		XuE(ii,jj,kk)
c		=============================================================================
c
		XuE(ii,jj,kk) =
     2		+ quarter * oE(i,jm1,kp1) * eighth
     2		+ half * oE(i,j,kp1) * quarter
     2		+ quarter * oE(i,jp1,kp1) * eighth
     2		+ quarter * uC(ip1,jm1,k) * eighth
c
     2		- eighth * (  oC(ip1,jm1,kp1)*eighth
     2				- uC(ip1,jm1,kp1) * quarter
     2				- oN(ip1,jm1,kp1) * quarter
     2				- oE(ip1,jm1,kp1) * quarter)

     2		+ half * uC(ip1,j,k) * quarter

     2		- quarter * (- oN(ip1,jm1,kp1)*eighth
     2				+ oC(ip1,j,kp1) * quarter
     2				- uC(ip1,j,kp1) * half
     2				- oN(ip1,j,kp1) * eighth
     2				- oE(ip1,j,kp1) * half)
c
     2		+ quarter * uC(ip1,jp1,k) * eighth
c
     2		- eighth * (- oN(ip1,j,kp1) * quarter
     2				+ oC(ip1,jp1,kp1) * eighth
     2				- uC(ip1,jp1,kp1) * quarter
     2				- oE(ip1,jp1,kp1) * quarter)
c
		if(ii.eq.nxm1.or.kk.eq.nzm1) XuE(ii,jj,kk) = 0.d0
c
c
c		XuW(ii,jj,kk)
c		=============================================================================
c
		XuW(ii,jj,kk) =
     2		+ quarter * uC(im1,jm1,k) * eighth
c
     2		- eighth * (- oE(im2,jm1,kp1)*quarter
     2				+ oC(im1,jm1,kp1) * eighth
     2				- uC(im1,jm1,kp1) * quarter
     2				- oN(im1,jm1,kp1) * quarter)
c
     2		+ half * uC(im1,j,k) * quarter
c
     2		- quarter * (- oE(im2,j,kp1) * half
     2				- oN(im1,jm1,kp1) * eighth
     2				+ oC(im1,j,kp1) * quarter
     2				- uC(im1,j,kp1) * half
     2				- oN(im1,j,kp1) * eighth)
c
     2		+ quarter * uC(im1,jp1,k) * eighth
c
     2		- eighth * (- oE(im2,jp1,kp1) *quarter
     2				- oN(im1,j,kp1) * quarter
     2				+ oC(im1,jp1,kp1) * eighth
     2				- uC(im1,jp1,kp1) * quarter)
c
     2		+ quarter * oE(im1,jm1,kp1) * eighth
     2		+ half * oE(im1,j,kp1) * quarter
     2		+ quarter * oE(im1,jp1,kp1) * eighth
c
		if(ii.eq.2) XuW(ii,jj,kk) = 0
		if(kk.eq.nzm1) XuW(ii,jj,kk) = 0.d0
c
c
c		XuN(ii,jj,kk)
c		=============================================================================
c
		XuN(ii,jj,kk) =
     2		+ quarter * oN(im1,j,kp1) * eighth
     2		+ quarter * uC(im1,jp1,k) * eighth
c
     2		- eighth * (  oC(im1,jp1,kp1)*eighth
     2				- uC(im1,jp1,kp1) * quarter
     2				- oN(im1,jp1,kp1) * quarter
     2				- oE(im1,jp1,kp1) * quarter)
c
     2		+ half * oN(i,j,kp1) * quarter
     2		+ half * uC(i,jp1,k) * quarter
c
     2		- quarter * (- oE(im1,jp1,kp1) *eighth
     2				+ oC(i,jp1,kp1) * quarter
     2				- uC(i,jp1,kp1) * half
     2				- oN(i,jp1,kp1) * half
     2				- oE(i,jp1,kp1) * eighth)
c
     2		+ quarter * oN(ip1,j,kp1) * eighth
     2		+ quarter * uC(ip1,jp1,k) * eighth
c
     2		- eighth * (- oE(i,jp1,kp1) * quarter
     2				+ oC(ip1,jp1,kp1) * eighth
     2				- uC(ip1,jp1,kp1) * quarter
     2				- oN(ip1,jp1,kp1) * quarter)
c
		if(jj.eq.nym1.or.kk.eq.nzm1) XuN(ii,jj,kk) = 0
c
c
c		XuS(ii,jj,kk)
c		=============================================================================
c
		XuS(ii,jj,kk) =
     2		+ quarter * uC(im1,jm1,k) * eighth
c
     2		- eighth * (- oN(im1,jm2,kp1)*quarter
     2				+ oC(im1,jm1,kp1) * eighth
     2				- uC(im1,jm1,kp1) * quarter
     2				- oE(im1,jm1,kp1) * quarter)
c
     2		+ quarter * oN(im1,jm1,kp1) * eighth
     2		+ half * uC(i,jm1,k) * quarter
c
     2		- quarter * (- oE(im1,jm1,kp1)*eighth
     2				- oN(i,jm2,kp1) * half
     2				+ oC(i,jm1,kp1) * quarter
     2				- uC(i,jm1,kp1) * half
     2				- oE(i,jm1,kp1) * eighth)
c
     2		+ half * oN(i,jm1,kp1) * quarter
     2		+ quarter * uC(ip1,jm1,k) * eighth
c
     2		- eighth * (- oE(i,jm1,kp1) * quarter
     2				- oN(ip1,jm2,kp1) * quarter
     2				+ oC(ip1,jm1,kp1) * eighth
     2				- uC(ip1,jm1,kp1) * quarter)
c
     2		+ quarter * oN(ip1,jm1,kp1) * eighth
c
		if(jj.eq.2) XuS(ii,jj,kk) = 0.d0
		if(kk.eq.nzm1) XuS(ii,jj,kk) = 0
c
c
c		XuNE(ii,jj,kk)
c		=============================================================================
c
		XuNE(ii,jj,kk) =
     2		+ quarter * oE(i,jp1,kp1) * eighth
     2		+ quarter * oN(ip1,j,kp1) * eighth
     2		+ quarter * uC(ip1,jp1,k) * eighth
c
     2		- eighth * (oC(ip1,jp1,kp1)*eighth
     2				- uC(ip1,jp1,kp1) * quarter
     2				- oN(ip1,jp1,kp1) * quarter
     2				- oE(ip1,jp1,kp1) * quarter)
c
		if(ii.eq.nxm1.or.jj.eq.nym1.or.kk.eq.nzm1) 
     1			XuNE(ii,jj,kk) = 0.d0
c
c
c		XuNW(ii,jj,kk)
c		=============================================================================
c
		XuNW(ii,jj,kk) =
     2		+ quarter * oN(im1,j,kp1) * eighth
     2		+ quarter * uC(im1,jp1,k) * eighth
c
     2		- eighth *(-oE(im2,jp1,kp1)*quarter
     2				+oC(im1,jp1,kp1) * eighth
     2				-uC(im1,jp1,kp1) * quarter
     2				-oN(im1,jp1,kp1) * quarter)
c
     2		+ quarter * oE(im1,jp1,kp1) * eighth
c
		if(ii.eq.2) XuNW(ii,jj,kk) = 0
		if(jj.eq.nym1.or.kk.eq.nzm1) XuNW(ii,jj,kk) = 0.d0
c
c
c		XuSE(ii,jj,kk)
c		=============================================================================
c
		XuSE(ii,jj,kk) =
     2		+ quarter * oE(i,jm1,kp1) * eighth
     2		+ quarter * uC(ip1,jm1,k) * eighth
c
     2		- eighth * (-oN(ip1,jm2,kp1)*quarter
     2				+ oC(ip1,jm1,kp1) * eighth
     2				- uC(ip1,jm1,kp1) * quarter
     2				- oE(ip1,jm1,kp1) * quarter)
c
     2		+ quarter * oN(ip1,jm1,kp1) * eighth
c
		if(jj.eq.2) XuSE(ii,jj,kk) = 0
		if(ii.eq.nxm1.or.kk.eq.nzm1) XuSE(ii,jj,kk) = 0.d0
c
c
c		XuSW(ii,jj,kk)
c		=============================================================================
c
		XuSW(ii,jj,kk) =
     2		+ quarter * uC(im1,jm1,k) * eighth
c
     2		- eighth*(- oE(im2,jm1,kp1)*quarter
     2				- oN(im1,jm2,kp1) * quarter
     2				+ oC(im1,jm1,kp1) * eighth
     2				- uC(im1,jm1,kp1)*quarter)
c
     2		+ quarter * oN(im1,jm1,kp1) * eighth
     2		+ quarter * oE(im1,jm1,kp1) * eighth
c
		if(ii.eq.2.or.jj.eq.2) XuSW(ii,jj,kk) = 0
		if(kk.eq.nzm1) XuSW(ii,jj,kk) = 0.d0
c
c
c		End 3 main loops
c
100	    continue
c
200	  continue
c
300	continue
c
	return
	end
c
c ==========================================================================================
c ==========================================================================================
c	Build_G27.f
c ==========================================================================================
c ==========================================================================================
c
c	Compute a 27-point Galerkin coarse grid matrix from a 27-point fine grid matrix
c
c
	subroutine build_G27(nxf,nyf,nzf,nx,ny,nz,
     4		oC,oE,oN,uC,oNE,oNW,uE,uW,uN,uS,uNE,uNW,uSE,uSW,
     5		XoC,XoE,XoN,XuC,XoNE,XoNW,XuE,XuW,XuN,XuS,XuNE,XuNW,
     6		XuSE,XuSW)
c
c ==========================================================================================
c	Declare all variables
c ==========================================================================================
c
c	Input/Output
c ==========================================================================================
c
	integer	nxf,nyf,nzf,nx,ny,nz
c
	real*8	oC(nxf,nyf,nzf),oE(nxf,nyf,nzf),oN(nxf,nyf,nzf)
	real*8	uC(nxf,nyf,nzf),uE(nxf,nyf,nzf),uN(nxf,nyf,nzf)
	real*8	oNE(nxf,nyf,nzf),oNW(nxf,nyf,nzf),uW(nxf,nyf,nzf)
	real*8	uNE(nxf,nyf,nzf),uNW(nxf,nyf,nzf),uSW(nxf,nyf,nzf)
	real*8	uSE(nxf,nyf,nzf),uS(nxf,nyf,nzf)
c
	real*8	XoC(nx,ny,nz),XoE(nx,ny,nz),XoN(nx,ny,nz)
	real*8	XuC(nx,ny,nz),XuE(nx,ny,nz),XuN(nx,ny,nz)
	real*8	XoNE(nx,ny,nz),XoNW(nx,ny,nz),XuW(nx,ny,nz)
	real*8	XuS(nx,ny,nz),XuNE(nx,ny,nz),XuNW(nx,ny,nz)
	real*8	XuSE(nx,ny,nz),XuSW(nx,ny,nz)
c
c	Local variable
c ==========================================================================================
c
	integer	i,j,k,ii,jj,kk
	integer	im1,ip1,im2,ip2,jm1,jp1,jm2,jp2,km1,kp1,km2,kp2
	integer	nxm1,nym1,nzm1
c
	real*8	whole,half,quarter,eighth
c
	real*8	tmp1_XoC,tmp2_XoC,tmp3_XoC,tmp4_XoC
	real*8	tmp5_XoC,tmp6_XoC,tmp7_XoC,tmp8_XoC
	real*8	tmp9_XoC,tmp10_XoC,tmp11_XoC,tmp12_XoC,tmp13_XoC
	real*8	tmp14_XoC,tmp15_XoC,tmp16_XoC,tmp17_XoC,tmp18_XoC
	real*8	tmp19_XoC,tmp20_XoC,tmp21_XoC,tmp22_XoC,tmp23_XoC
	real*8	tmp24_XoC,tmp25_XoC,tmp26_XoC,tmp27_XoC
	real*8	tmp1_XoE,tmp2_XoE,tmp3_XoE,tmp4_XoE
	real*8	tmp5_XoE,tmp6_XoE,tmp7_XoE,tmp8_XoE
	real*8	tmp9_XoE,tmp10_XoE,tmp11_XoE,tmp12_XoE
	real*8	tmp1_XoN,tmp2_XoN,tmp3_XoN,tmp4_XoN
	real*8	tmp5_XoN,tmp6_XoN,tmp7_XoN,tmp8_XoN
	real*8	tmp9_XoN,tmp10_XoN,tmp11_XoN,tmp12_XoN
	real*8	tmp1_XuC,tmp2_XuC,tmp3_XuC,tmp4_XuC
	real*8	tmp5_XuC,tmp6_XuC,tmp7_XuC,tmp8_XuC
	real*8	tmp9_XuC,tmp10_XuC,tmp11_XuC,tmp12_XuC
	real*8	tmp1_XoNE,tmp2_XoNE,tmp3_XoNE,tmp4_XoNE
	real*8	tmp5_XoNE,tmp6_XoNE
	real*8	tmp1_XoNW,tmp2_XoNW,tmp3_XoNW,tmp4_XoNW
	real*8	tmp5_XoNW,tmp6_XoNW
	real*8	tmp1_XuE,tmp2_XuE,tmp3_XuE,tmp4_XuE
	real*8	tmp5_XuE,tmp6_XuE
	real*8	tmp1_XuW,tmp2_XuW,tmp3_XuW,tmp4_XuW
	real*8	tmp5_XuW,tmp6_XuW
	real*8	tmp1_XuN,tmp2_XuN,tmp3_XuN,tmp4_XuN
	real*8	tmp5_XuN,tmp6_XuN
	real*8	tmp1_XuS,tmp2_XuS,tmp3_XuS,tmp4_XuS
	real*8	tmp5_XuS,tmp6_XuS
	real*8	tmp1_XuNE,tmp2_XuNE,tmp1_XuNW,tmp2_XuNW
	real*8	tmp1_XuSE,tmp2_XuSE,tmp1_XuSW,tmp2_XuSW
c
c ==========================================================================================
c	Build the operator
c ==========================================================================================
c
c	Setup
c
	whole = 1.d0
	half = whole/2.d0
	quarter = half/2.d0
	eighth = quarter/2.d0
c
	nxm1 = nx - 1
	nym1 = ny - 1
	nzm1 = nz - 1
c
	do 300 kk = 2,nzm1
	  k = 2*kk-1
c
	  do 200 jj =2,nym1
	    j = 2*jj-1
c
	    do 100 ii = 2,nxm1
		i = 2*ii - 1
c
c		Set indicies
c              =============================================================================
c
		im1 = i-1
		ip1 = i+1
		im2 = i-2
		ip2 = i+2
		jm1 = j-1
		jp1 = j+1
		jm2 = j-2
		jp2 = j+2
		km1 = k-1
		kp1 = k+1
		km2 = k-2
		kp2 = k+2
c
c		XoC(ii,jj,kk)
c		=============================================================================
c
		tmp1_XoC =
     2		+ half * (- uNE(im1,j,km1) * quarter
     2			- oNE(im1,j,k) * half
     2			- uSW(i,jp1,k) * quarter
     2			- uE(im1,jp1,km1) * eighth
     2			- oE(im1,jp1,k) * quarter
     2			- uW(i,jp1,k) * eighth
     2			- uN(i,j,km1) * half
     2			- oN(i,j,k) * whole
     2			- uS(i,jp1,k) * half
     2			- uC(i,jp1,km1) * quarter
     2			+ oC(i,jp1,k) * half
     2			- uC(i,jp1,k) * quarter
     2			- uNW(ip1,j,km1) * quarter
     2			- oNW(ip1,j,k) * half
     2			- uSE(i,jp1,k) * quarter
     2			- uW(ip1,jp1,km1) * eighth
     2			- oE(i,jp1,k) * quarter
     2			- uE(i,jp1,k) * eighth) 
c
c
		tmp2_XoC =
     2		+ quarter * (- oNE(im1,j,km1) * quarter
     2			- uSW(i,jp1,km1) * half
     2			- oE(im1,jp1,km1) * eighth
     2			- uW(i,jp1,km1) * quarter 
     2			- oN(i,j,km1) * half
     2			- uS(i,jp1,km1) * whole 
     2			+ oC(i,jp1,km1) * quarter
     2			- uC(i,jp1,km1) * half 
     2			- oNW(ip1,j,km1) * quarter
     2			- uSE(i,jp1,km1) * half 
     2			- oE(i,jp1,km1) * eighth
     2			- uE(i,jp1,km1) * quarter) 
c
c
		tmp3_XoC =
     2		+ half * (- oNE(im1,jm1,km1) * eighth 
     2			- uSW(i,j,km1) * quarter
     2			- oE(im1,j,km1) * quarter 
     2			- uW(i,j,km1) * half
     2			- oNW(i,j,km1) * eighth 
     2			- uNW(i,j,km1) * quarter
     2			- oN(i,jm1,km1) * quarter 
     2			- uS(i,j,km1) * half
     2			+ oC(i,j,km1) * half 
     2			- uC(i,j,km1) * whole
     2			- oN(i,j,km1) * quarter 
     2			- uN(i,j,km1) * half
     2			- oNW(ip1,jm1,km1) * eighth 
     2			- uSE(i,j,km1) * quarter
     2			- oE(i,j,km1) * quarter 
     2			- uE(i,j,km1) * half
     2			- oNE(i,j,km1) * eighth 
     2			- uNE(i,j,km1) * quarter) 
c
c
		tmp4_XoC =
     2		+ half * (- uNE(im1,jm1,k) * quarter
     2			- oNE(im1,jm1,kp1) * eighth 
     2			- uE(im1,j,k) * half
     2			- oE(im1,j,kp1) * quarter 
     2			- uSE(im1,jp1,k) * quarter
     2			- oNW(i,j,kp1) * eighth 
     2			- uN(i,jm1,k) * half
     2			- oN(i,jm1,kp1) * quarter 
     2			- uC(i,j,k) * whole
     2			+ oC(i,j,kp1) * half 
     2			- uS(i,jp1,k) * half
     2			- oN(i,j,kp1) * quarter 
     2			- uNW(ip1,jm1,k) * quarter
     2			- oNW(ip1,jm1,kp1) * eighth 
     2			- uW(ip1,j,k) * half
     2			- oE(i,j,kp1) * quarter 
     2			- uSW(ip1,jp1,k) * quarter
     2			- oNE(i,j,kp1) * eighth) 
c
c
		tmp5_XoC =
     2		+ whole * (- uW(ip1,j,km1) * quarter 
     2			- oE(im1,j,k) * half
     2			- uSE(im1,jp1,km1) * eighth
     2			- uNE(im1,jm1,km1) * eighth
     2			- uN(i,jm1,km1) * quarter 
     2			- oNE(im1,jm1,k) * quarter
     2			- uE(im1,j,km1) * quarter 
     2			- oNW(ip1,jm1,k) * quarter
     2			- uC(i,j,km1) * half 
     2			- uNW(ip1,jm1,km1) * eighth 
     2			- uSW(ip1,jp1,km1) * eighth
     2		- uS(i,jp1,km1)*quarter-oN(i,jm1,k)*half
     2		- uNE(i,j,k)*eighth-oNE(i,j,k)*quarter
     2		- uE(i,j,k)*quarter-uSE(i,j,k) * eighth
     2		- oN(i,j,k)*half-oE(i,j,k) * half
     2		- uS(i,j,k)*quarter+oC(i,j,k) * whole
     2		- uSW(i,j,k)*eighth-uN(i,j,k) * quarter
     2		- uC(i,j,k)*half-uW(i,j,k) * quarter
     2		- oNW(i,j,k)*quarter-uNW(i,j,k)*eighth) 
c
c
		tmp6_XoC =
     2		+ quarter * (- uE(im1,jm1,k) * quarter
     2			- oE(im1,jm1,kp1) * eighth 
     2			- uSE(im1,j,k) * half
     2			- oNW(i,jm1,kp1) * quarter 
     2			- uC(i,jm1,k) * half
     2			+ oC(i,jm1,kp1) * quarter 
     2			- uS(i,j,k) * whole
     2			- oN(i,jm1,kp1) * half 
     2			- uW(ip1,jm1,k) * quarter
     2			- oE(i,jm1,kp1) * eighth 
     2			- uSW(ip1,j,k) * half
     2			- oNE(i,jm1,kp1) * quarter) 
c
c
		tmp7_XoC =
     2		+ half * (- uE(im1,jm1,km1) * eighth 
     2			- oE(im1,jm1,k) * quarter
     2			- uW(i,jm1,k) * eighth 
     2			- uSE(im1,j,km1) * quarter
     2			- oNW(i,jm1,k) * half 
     2			- uNW(i,jm1,k) * quarter
     2			- uC(i,jm1,km1) * quarter 
     2			+ oC(i,jm1,k) * half
     2			- uC(i,jm1,k) * quarter 
     2			- uS(i,j,km1) * half
     2			- oN(i,jm1,k) * whole 
     2			- uN(i,jm1,k) * half
     2			- uW(ip1,jm1,km1) * eighth 
     2			- oE(i,jm1,k) * quarter
     2			- uE(i,jm1,k) * eighth 
     2			- uSW(ip1,j,km1) * quarter
     2			- oNE(i,jm1,k) * half 
     2			- uNE(i,jm1,k) * quarter) 
c
c
		tmp8_XoC =
     2		+ quarter * (- oE(im1,jm1,km1) * eighth
     2			- uW(i,jm1,km1) * quarter 
     2			- oNW(i,jm1,km1) * quarter
     2			- uNW(i,jm1,km1) * half 
     2			+ oC(i,jm1,km1) * quarter
     2			- uC(i,jm1,km1) * half 
     2			- oN(i,jm1,km1) * half
     2			- uN(i,jm1,km1) * whole 
     2			- oE(i,jm1,km1) * eighth
     2			- uE(i,jm1,km1) * quarter 
     2			- oNE(i,jm1,km1) * quarter
     2			- uNE(i,jm1,km1) * half) 
c
c
		tmp9_XoC =
     2		+ eighth * (- uN(im1,j,k) * half 
     2			 - oN(im1,j,kp1) * quarter
     2			 - uC(im1,jp1,k) * quarter
     2			 + oC(im1,jp1,kp1) * eighth 
     2			 - uNW(i,j,k) * whole
     2			 - oNW(i,j,kp1) * half 
     2			 - uW(i,jp1,k) * half
     2			 - oE(im1,jp1,kp1) * quarter) 
c
c
		tmp10_XoC =
     2		+ quarter * (- uN(im1,j,km1) * quarter 
     2			 - oN(im1,j,k) * half
     2			 - uS(im1,jp1,k) * quarter 
     2			 - uC(im1,jp1,km1) * eighth
     2			 + oC(im1,jp1,k) * quarter 
     2			 - uC(im1,jp1,k) * eighth
     2			 - uNW(i,j,km1) * half 
     2			 - oNW(i,j,k) * whole
     2			 - uSE(im1,jp1,k) * half 
     2			 - uW(i,jp1,km1) * quarter
     2			 - oE(im1,jp1,k) * half 
     2			 - uE(im1,jp1,k) * quarter) 
c
c
		tmp11_XoC =
     2		+ quarter * (- uN(im1,jm1,k) * quarter
     2			- oN(im1,jm1,kp1) * eighth 
     2			- uC(im1,j,k) * half
     2			+ oC(im1,j,kp1) * quarter 
     2			- uS(im1,jp1,k) * quarter
     2			- oN(im1,j,kp1) * eighth 
     2			- uNW(i,jm1,k) * half
     2			- oNW(i,jm1,kp1) * quarter 
     2			- uW(i,j,k) * whole
     2			- oE(im1,j,kp1) * half 
     2			- uSW(i,jp1,k) * half
     2			- oNE(im1,j,kp1) * quarter) 
c
c
		tmp12_XoC =
     2		+ eighth * (- oN(im1,j,km1) * quarter 
     2			 - uS(im1,jp1,km1) * half
     2			 + oC(im1,jp1,km1) * eighth
     2			 - uC(im1,jp1,km1) * quarter 
     2			 - oNW(i,j,km1) * half
     2			 - uSE(im1,jp1,km1) * whole
     2			 - oE(im1,jp1,km1) * quarter
     2			 - uE(im1,jp1,km1) * half) 
c
c
		tmp13_XoC =
     2		+ half * (- uN(im1,jm1,km1) * eighth 
     2			- oN(im1,jm1,k) * quarter
     2			- uS(im1,j,k) * eighth 
     2			- uC(im1,j,km1) * quarter
     2			+ oC(im1,j,k) * half 
     2			- uC(im1,j,k) * quarter
     2			- uS(im1,jp1,km1) * eighth 
     2			- oN(im1,j,k) * quarter
     2			- uN(im1,j,k) * eighth 
     2			- uNW(i,jm1,km1) * quarter
     2			- oNW(i,jm1,k) * half 
     2			- uSE(im1,j,k) * quarter
     2			- uW(i,j,km1) * half 
     2			- oE(im1,j,k) * whole
     2			- uE(im1,j,k) * half 
     2			- uSW(i,jp1,km1) * quarter
     2			- oNE(im1,j,k) * half 
     2			- uNE(im1,j,k) * quarter) 
c
c
		tmp14_XoC =
     2		+ eighth * (- uC(im1,jm1,k) * quarter
     2			 + oC(im1,jm1,kp1) * eighth 
     2			 - uS(im1,j,k) * half
     2			 - oN(im1,jm1,kp1) * quarter 
     2			 - uW(i,jm1,k) * half
     2			 - oE(im1,jm1,kp1) * quarter 
     2			 - uSW(i,j,k) * whole
     2			 - oNE(im1,jm1,kp1) * half) 
c
c
		tmp15_XoC =
     2		+ quarter * (- uC(im1,jm1,km1) * eighth 
     2			 + oC(im1,jm1,k) * quarter
     2			 - uC(im1,jm1,k) * eighth 
     2			 - uS(im1,j,km1) * quarter
     2			 - oN(im1,jm1,k) * half 
     2			 - uN(im1,jm1,k) * quarter
     2			 - uW(i,jm1,km1) * quarter 
     2			 - oE(im1,jm1,k) * half
     2			 - uE(im1,jm1,k) * quarter 
     2			 - uSW(i,j,km1) * half
     2			 - oNE(im1,jm1,k) * whole 
     2			 - uNE(im1,jm1,k) * half)
c
c
		tmp16_XoC =
     2		+ quarter * (- oN(im1,jm1,km1) * eighth
     2			- uS(im1,j,km1) * quarter 
     2			+ oC(im1,j,km1) * quarter
     2			- uC(im1,j,km1) * half 
     2			- oN(im1,j,km1) * eighth
     2			- uN(im1,j,km1) * quarter 
     2			- oNW(i,jm1,km1) * quarter
     2			- uSE(im1,j,km1) * half 
     2			- oE(im1,j,km1) * half
     2			- uE(im1,j,km1) * whole 
     2			- oNE(im1,j,km1) * quarter
     2			- uNE(im1,j,km1) * half) 
c
c
		tmp17_XoC =
     2		+ eighth * (- uNE(i,j,k) * whole 
     2			 - oNE(i,j,kp1) * half
     2			 - uE(i,jp1,k) * half 
     2			 - oE(i,jp1,kp1) * quarter
     2			 - uN(ip1,j,k) * half 
     2			 - oN(ip1,j,kp1) * quarter
     2			 - uC(ip1,jp1,k) * quarter
     2			 + oC(ip1,jp1,kp1) * eighth) 
c
c
		tmp18_XoC =
     2		+ quarter * (- uNE(i,jm1,k) * half 
     2			- oNE(i,jm1,kp1) * quarter
     2			- uE(i,j,k) * whole 
     2			- oE(i,j,kp1) * half
     2			- uSE(i,jp1,k) * half 
     2			- oNW(ip1,j,kp1) * quarter
     2			- uN(ip1,jm1,k) * quarter
     2			- oN(ip1,jm1,kp1) * eighth 
     2			- uC(ip1,j,k) * half
     2			+ oC(ip1,j,kp1) * quarter 
     2			- uS(ip1,jp1,k) * quarter
     2			- oN(ip1,j,kp1) * eighth) 
c
c
		tmp19_XoC =
     2		+ eighth * (- oNE(i,j,km1) * half 
     2			 - uSW(ip1,jp1,km1) * whole
     2			 - oE(i,jp1,km1) * quarter 
     2			 - uW(ip1,jp1,km1) * half
     2			 - oN(ip1,j,km1) * quarter 
     2			 - uS(ip1,jp1,km1) * half
     2			 + oC(ip1,jp1,km1) * eighth
     2			 - uC(ip1,jp1,km1) * quarter) 
c
c
		tmp20_XoC =
     2		+ quarter * (- uNE(i,j,km1) * half 
     2			 - oNE(i,j,k) * whole
     2			 - uSW(ip1,jp1,k) * half 
     2			 - uE(i,jp1,km1) * quarter
     2			 - oE(i,jp1,k) * half 
     2			 - uW(ip1,jp1,k) * quarter
     2			 - uN(ip1,j,km1) * quarter 
     2			 - oN(ip1,j,k) * half
     2			 - uS(ip1,jp1,k) * quarter 
     2			 - uC(ip1,jp1,km1) * eighth
     2			 + oC(ip1,jp1,k) * quarter 
     2			 - uC(ip1,jp1,k) * eighth)
c
c
		tmp21_XoC =
     2		+ quarter * (- uE(i,jm1,km1) * quarter
     2			 - oE(i,jm1,k) * half 
     2			 - uW(ip1,jm1,k) * quarter
     2			 - uSE(i,j,km1) * half 
     2			 - oNW(ip1,jm1,k) * whole
     2			 - uNW(ip1,jm1,k) * half
     2			 - uC(ip1,jm1,km1) * eighth
     2			 + oC(ip1,jm1,k) * quarter 
     2			 - uC(ip1,jm1,k) * eighth
     2			 - uS(ip1,j,km1) * quarter 
     2			 - oN(ip1,jm1,k) * half
     2			 - uN(ip1,jm1,k) * quarter) 
c
c
		tmp22_XoC =
     2		+ eighth * (- oE(i,jm1,km1) * quarter 
     2			 - uW(ip1,jm1,km1) * half
     2			 - oNW(ip1,jm1,km1) * half
     2			 - uNW(ip1,jm1,km1) * whole
     2			 + oC(ip1,jm1,km1) * eighth
     2			 - uC(ip1,jm1,km1) * quarter
     2			 - oN(ip1,jm1,km1) * quarter
     2			 - uN(ip1,jm1,km1) * half) 
c
c
		tmp23_XoC =
     2		+ eighth * (- uE(i,jm1,k) * half 
     2			 - oE(i,jm1,kp1) * quarter
     2			 - uSE(i,j,k) * whole 
     2			 - oNW(ip1,jm1,kp1) * half
     2			 - uC(ip1,jm1,k) * quarter
     2			 + oC(ip1,jm1,kp1) * eighth 
     2			 - uS(ip1,j,k) * half
     2			 - oN(ip1,jm1,kp1) * quarter) 
c
c
		tmp24_XoC =
     2		+ half * (- uNE(i,jm1,km1) * quarter 
     2			- oNE(i,jm1,k) * half
     2			- uSW(ip1,j,k) * quarter 
     2			- uE(i,j,km1) * half
     2			- oE(i,j,k) * whole 
     2			- uW(ip1,j,k) * half
     2			- uSE(i,jp1,km1) * quarter 
     2			- oNW(ip1,j,k) * half
     2			- uNW(ip1,j,k) * quarter 
     2			- uN(ip1,jm1,km1) * eighth
     2			- oN(ip1,jm1,k) * quarter 
     2			- uS(ip1,j,k) * eighth
     2			- uC(ip1,j,km1) * quarter 
     2			+ oC(ip1,j,k) * half
     2			- uC(ip1,j,k) * quarter 
     2			- uS(ip1,jp1,km1) * eighth
     2			- oN(ip1,j,k) * quarter 
     2			- uN(ip1,j,k) * eighth) 
c
c
		tmp25_XoC =
     2		+ quarter * (- oNE(i,jm1,km1) * quarter
     2			- uSW(ip1,j,km1) * half 
     2			- oE(i,j,km1) * half
     2			- uW(ip1,j,km1) * whole 
     2			- oNW(ip1,j,km1) * quarter
     2			- uNW(ip1,j,km1) * half
     2			- oN(ip1,jm1,km1) * eighth
     2			- uS(ip1,j,km1) * quarter 
     2			+ oC(ip1,j,km1) * quarter
     2			- uC(ip1,j,km1) * half 
     2			- oN(ip1,j,km1) * eighth
     2			- uN(ip1,j,km1) * quarter) 
c
c
		tmp26_XoC =
     2		+ quarter * (- uNE(im1,j,k) * half 
     2			- oNE(im1,j,kp1) * quarter
     2			- uE(im1,jp1,k) * quarter
     2			- oE(im1,jp1,kp1) * eighth 
     2			- uN(i,j,k) * whole
     2			- oN(i,j,kp1) * half 
     2			- uC(i,jp1,k) * half
     2			+ oC(i,jp1,kp1) * quarter 
     2			- uNW(ip1,j,k) * half
     2			- oNW(ip1,j,kp1) * quarter 
     2			- uW(ip1,jp1,k) * quarter
     2			- oE(i,jp1,kp1) * eighth) 
c
c
		tmp27_XoC =
     2		+ eighth * ( oC(im1,jm1,km1) * eighth
     2			- uC(im1,jm1,km1) * quarter
     2			- oN(im1,jm1,km1) * quarter
     2			- uN(im1,jm1,km1) * half
     2			- oE(im1,jm1,km1) * quarter
     2			- uE(im1,jm1,km1) * half
     2			- oNE(im1,jm1,km1) * half
     2			- uNE(im1,jm1,km1) * whole)
c
		XoC(ii,jj,kk) = tmp1_XoC + tmp2_XoC + tmp3_XoC +tmp4_XoC
     2		+ tmp5_XoC + tmp6_XoC + tmp7_XoC + tmp8_XoC + tmp9_XoC 
     3		+ tmp10_XoC + tmp11_XoC + tmp12_XoC+tmp13_XoC+tmp14_XoC 
     4		+ tmp15_XoC + tmp16_XoC + tmp17_XoC+tmp18_XoC+tmp19_XoC 
     5		+ tmp20_XoC + tmp21_XoC + tmp22_XoC+tmp23_XoC+tmp24_XoC 
     6		+ tmp25_XoC + tmp26_XoC + tmp27_XoC
c
c		XoE(ii,jj,kk)
c		=============================================================================
c
		tmp1_XoE =
     2		- quarter * (- oE(i,jm1,km1) * eighth
     2			- uE(i,jm1,km1) * quarter
     2			- oNE(i,jm1,km1) * quarter
     2			- uNE(i,jm1,km1) * half) 
c
     2		- half * (- uW(ip1,jm1,km1) * eighth
     2			- oE(i,jm1,k) * quarter 
     2			- uE(i,jm1,k) * eighth
     2			- uSW(ip1,j,km1) * quarter
     2			- oNE(i,jm1,k) * half 
     2			- uNE(i,jm1,k) * quarter) 
c
     2		- quarter * (- uW(ip1,jm1,k) * quarter
     2			- oE(i,jm1,kp1) * eighth
     2			- uSW(ip1,j,k) * half
     2			- oNE(i,jm1,kp1) * quarter) 
c
c
		tmp2_XoE =
     2		- half * (- oNW(ip1,jm1,km1) * eighth
     2			- uSE(i,j,km1) * quarter 
     2			- oE(i,j,km1) * quarter
     2			- uE(i,j,km1) * half 
     2			- oNE(i,j,km1) * eighth
     2			- uNE(i,j,km1) * quarter) 
c
     2		- whole * (- uNW(ip1,jm1,km1) * eighth
     2			- oNW(ip1,jm1,k) * quarter 
     2			- uSE(i,j,k) * eighth
     2			- uW(ip1,j,km1) * quarter 
     2			- oE(i,j,k) * half
     2			- uE(i,j,k) * quarter
     2			- uSW(ip1,jp1,km1) * eighth
     2			- oNE(i,j,k) * quarter 
     2			- uNE(i,j,k) * eighth) 
c
c
		tmp3_XoE =
     2		- half * (- uNW(ip1,jm1,k) * quarter
     2			- oNW(ip1,jm1,kp1) * eighth
     2			- uW(ip1,j,k) * half 
     2			- oE(i,j,kp1) * quarter
     2			- uSW(ip1,jp1,k) * quarter
     2			- oNE(i,j,kp1) * eighth) 
c
     2		- quarter * (- oNW(ip1,j,km1) * quarter
     2			- uSE(i,jp1,km1) * half
     2			- oE(i,jp1,km1) * eighth
     2			- uE(i,jp1,km1) * quarter) 
c
     2		- half * (- uNW(ip1,j,km1) * quarter 
     2			- oNW(ip1,j,k) * half
     2			- uSE(i,jp1,k) * quarter
     2			- uW(ip1,jp1,km1) * eighth
     2			- oE(i,jp1,k) * quarter 
     2			- uE(i,jp1,k) * eighth) 
c
c
		tmp4_XoE =
     2		- quarter * (- uNW(ip1,j,k) * half
     2			- oNW(ip1,j,kp1) * quarter
     2			- uW(ip1,jp1,k) * quarter
     2			- oE(i,jp1,kp1) * eighth) 
c
     2		- eighth * (  oC(ip1,jm1,km1) * eighth
     2			 - uC(ip1,jm1,km1) * quarter
     2			 - oN(ip1,jm1,km1) * quarter
     2			 - uN(ip1,jm1,km1) * half
     2			 - oE(ip1,jm1,km1) * quarter
     2			 - uE(ip1,jm1,km1) * half
     2			 - oNE(ip1,jm1,km1) * half
     2			 - uNE(ip1,jm1,km1) * whole) 
c
c
		tmp5_XoE =
     2		- quarter * (- uC(ip1,jm1,km1) * eighth
     2			 + oC(ip1,jm1,k) * quarter
     2			 - uC(ip1,jm1,k) * eighth
     2			 - uS(ip1,j,km1) * quarter
     2			 - oN(ip1,jm1,k) * half
     2			 - uN(ip1,jm1,k) * quarter
     2			 - uW(ip2,jm1,km1) * quarter
     2			 - oE(ip1,jm1,k) * half
     2			 - uE(ip1,jm1,k) * quarter
     2			 - uSW(ip2,j,km1) * half
     2			 - oNE(ip1,jm1,k) * whole
     2			 - uNE(ip1,jm1,k) * half) 
c
c
		tmp6_XoE =
     2		- eighth * (- uC(ip1,jm1,k) * quarter
     2			 + oC(ip1,jm1,kp1) * eighth
     2			 - uS(ip1,j,k) * half
     2			 - oN(ip1,jm1,kp1) * quarter
     2			 - uW(ip2,jm1,k) * half
     2			 - oE(ip1,jm1,kp1) * quarter
     2			 - uSW(ip2,j,k) * whole
     2			 - oNE(ip1,jm1,kp1) * half) 
c
c
		tmp7_XoE =
     2		- quarter * (- oN(ip1,jm1,km1) * eighth
     2			- uS(ip1,j,km1) * quarter
     2			+ oC(ip1,j,km1) * quarter
     2			- uC(ip1,j,km1) * half
     2			- oN(ip1,j,km1) * eighth
     2			- uN(ip1,j,km1) * quarter
     2			- oNW(ip2,jm1,km1) * quarter
     2			- uSE(ip1,j,km1) * half
     2			- oE(ip1,j,km1) * half
     2			- uE(ip1,j,km1) * whole
     2			- oNE(ip1,j,km1) * quarter
     2			- uNE(ip1,j,km1) * half) 
c
c
		tmp8_XoE =
     2		- half * (- uN(ip1,jm1,km1) * eighth
     2			- oN(ip1,jm1,k) * quarter
     2			- uS(ip1,j,k) * eighth 
     2			- uC(ip1,j,km1) * quarter
     2			+ oC(ip1,j,k) * half 
     2			- uC(ip1,j,k) * quarter
     2			- uS(ip1,jp1,km1) * eighth
     2			- oN(ip1,j,k) * quarter 
     2			- uN(ip1,j,k) * eighth
     2			- uNW(ip2,jm1,km1) * quarter
     2			- oNW(ip2,jm1,k) * half
     2			- uSE(ip1,j,k) * quarter 
     2			- uW(ip2,j,km1) * half
     2			- oE(ip1,j,k) * whole 
     2			- uE(ip1,j,k) * half
     2			- uSW(ip2,jp1,km1) * quarter
     2			- oNE(ip1,j,k) * half 
     2			- uNE(ip1,j,k) * quarter) 
c
c
		tmp9_XoE =
     2		- quarter * (- uN(ip1,jm1,k) * quarter
     2			- oN(ip1,jm1,kp1) * eighth
     2			- uC(ip1,j,k) * half 
     2			+ oC(ip1,j,kp1) * quarter
     2			- uS(ip1,jp1,k) * quarter
     2			- oN(ip1,j,kp1) * eighth
     2			- uNW(ip2,jm1,k) * half
     2			- oNW(ip2,jm1,kp1) * quarter
     2			- uW(ip2,j,k) * whole 
     2			- oE(ip1,j,kp1) * half
     2			- uSW(ip2,jp1,k) * half
     2			- oNE(ip1,j,kp1) * quarter) 
c
c
		tmp10_XoE =
     2		- eighth * (- oN(ip1,j,km1) * quarter
     2			 - uS(ip1,jp1,km1) * half
     2			 + oC(ip1,jp1,km1) * eighth
     2			 - uC(ip1,jp1,km1) * quarter
     2			 - oNW(ip2,j,km1) * half
     2			 - uSE(ip1,jp1,km1) * whole
     2			 - oE(ip1,jp1,km1) * quarter
     2			 - uE(ip1,jp1,km1) * half)
c
c
		tmp11_XoE =
     2		- quarter * (- uN(ip1,j,km1) * quarter 
     2			 - oN(ip1,j,k) * half
     2			 - uS(ip1,jp1,k) * quarter
     2			 - uC(ip1,jp1,km1) * eighth
     2			 + oC(ip1,jp1,k) * quarter
     2			 - uC(ip1,jp1,k) * eighth
     2			 - uNW(ip2,j,km1) * half
     2			 - oNW(ip2,j,k) * whole
     2			 - uSE(ip1,jp1,k) * half
     2			 - uW(ip2,jp1,km1) * quarter
     2			 - oE(ip1,jp1,k) * half
     2			 - uE(ip1,jp1,k) * quarter) 
c
c
		tmp12_XoE =
     2		- eighth * (- uN(ip1,j,k) * half 
     2			 - oN(ip1,j,kp1) * quarter
     2			 - uC(ip1,jp1,k) * quarter
     2			 + oC(ip1,jp1,kp1) * eighth
     2			 - uNW(ip2,j,k) * whole
     2			 - oNW(ip2,j,kp1) * half
     2			 - uW(ip2,jp1,k) * half
     2			 - oE(ip1,jp1,kp1) * quarter)
c
		XoE(ii,jj,kk) = tmp1_XoE + tmp2_XoE + tmp3_XoE +tmp4_XoE
     2		+ tmp5_XoE + tmp6_XoE + tmp7_XoE + tmp8_XoE + tmp9_XoE
     3		+ tmp10_XoE + tmp11_XoE + tmp12_XoE
c
		if(ii.eq.nxm1) XoE(ii,jj,kk) = 0.d0
c
c		XoN(ii,jj,kk)
c		=============================================================================
c
		tmp1_XoN =
     2		- quarter * (- oN(im1,j,km1) * eighth
     2			- uN(im1,j,km1) * quarter
     2			- oNE(im1,j,km1) * quarter
     2			- uNE(im1,j,km1) * half) 
c
     2		- half * (- uS(im1,jp1,km1) * eighth
     2			- oN(im1,j,k) * quarter 
     2			- uN(im1,j,k) * eighth
     2			- uSW(i,jp1,km1) * quarter
     2			- oNE(im1,j,k) * half 
     2			- uNE(im1,j,k) * quarter) 
c
     2		- quarter * (- uS(im1,jp1,k) * quarter
     2			- oN(im1,j,kp1) * eighth
     2			- uSW(i,jp1,k) * half
     2			- oNE(im1,j,kp1) * quarter) 
c
c
		tmp2_XoN =
     2		- eighth * (  oC(im1,jp1,km1) * eighth
     2			 - uC(im1,jp1,km1) * quarter
     2			 - oN(im1,jp1,km1) * quarter
     2			 - uN(im1,jp1,km1) * half
     2			 - oE(im1,jp1,km1) * quarter
     2			 - uE(im1,jp1,km1) * half
     2			 - oNE(im1,jp1,km1) * half
     2			 - uNE(im1,jp1,km1) * whole) 
c
c
		tmp3_XoN =
     2		- quarter * (- uC(im1,jp1,km1) * eighth
     2			 + oC(im1,jp1,k) * quarter
     2			 - uC(im1,jp1,k) * eighth
     2			 - uS(im1,jp2,km1) * quarter
     2			 - oN(im1,jp1,k) * half
     2			 - uN(im1,jp1,k) * quarter
     2			 - uW(i,jp1,km1) * quarter
     2			 - oE(im1,jp1,k) * half
     2			 - uE(im1,jp1,k) * quarter
     2			 - uSW(i,jp2,km1) * half
     2			 - oNE(im1,jp1,k) * whole
     2			 - uNE(im1,jp1,k) * half) 
c
c
		tmp4_XoN =
     2		- eighth * (- uC(im1,jp1,k) * quarter
     2			 + oC(im1,jp1,kp1) * eighth
     2			 - uS(im1,jp2,k) * half
     2			 - oN(im1,jp1,kp1) * quarter
     2			 - uW(i,jp1,k) * half
     2			 - oE(im1,jp1,kp1) * quarter
     2			 - uSW(i,jp2,k) * whole
     2			 - oNE(im1,jp1,kp1) * half) 
c
     2		- half * (- oNW(i,j,km1) * eighth 
     2			- uNW(i,j,km1) * quarter
     2			- oN(i,j,km1) * quarter 
     2			- uN(i,j,km1) * half
     2			- oNE(i,j,km1) * eighth 
     2			- uNE(i,j,km1) * quarter)
c
c
		tmp5_XoN =
     2		- whole * (- uSE(im1,jp1,km1) * eighth
     2			- oNW(i,j,k) * quarter 
     2			- uNW(i,j,k) * eighth
     2			- uS(i,jp1,km1) * quarter 
     2			- oN(i,j,k) * half
     2			- uN(i,j,k) * quarter
     2			- uSW(ip1,jp1,km1) * eighth
     2			- oNE(i,j,k) * quarter 
     2			- uNE(i,j,k) * eighth) 
c
     2		- half * (- uSE(im1,jp1,k) * quarter
     2			- oNW(i,j,kp1) * eighth 
     2			- uS(i,jp1,k) * half
     2			- oN(i,j,kp1) * quarter
     2			- uSW(ip1,jp1,k) * quarter
     2			- oNE(i,j,kp1) * eighth) 
c
c
		tmp6_XoN =
     2		- quarter * (- oE(im1,jp1,km1) * eighth
     2			- uW(i,jp1,km1) * quarter
     2			- oNW(i,jp1,km1) * quarter
     2			- uNW(i,jp1,km1) * half
     2			+ oC(i,jp1,km1) * quarter
     2			- uC(i,jp1,km1) * half
     2			- oN(i,jp1,km1) * half
     2			- uN(i,jp1,km1) * whole
     2			- oE(i,jp1,km1) * eighth
     2			- uE(i,jp1,km1) * quarter
     2			- oNE(i,jp1,km1) * quarter
     2			- uNE(i,jp1,km1) * half) 
c
c
		tmp7_XoN =
     2		- half * (- uE(im1,jp1,km1) * eighth
     2			- oE(im1,jp1,k) * quarter
     2			- uW(i,jp1,k) * eighth
     2			- uSE(im1,jp2,km1) * quarter
     2			- oNW(i,jp1,k) * half 
     2			- uNW(i,jp1,k) * quarter
     2			- uC(i,jp1,km1) * quarter 
     2			+ oC(i,jp1,k) * half
     2			- uC(i,jp1,k) * quarter 
     2			- uS(i,jp2,km1) * half
     2			- oN(i,jp1,k) * whole 
     2			- uN(i,jp1,k) * half
     2			- uW(ip1,jp1,km1) * eighth
     2			- oE(i,jp1,k) * quarter 
     2			- uE(i,jp1,k) * eighth
     2			- uSW(ip1,jp2,km1) * quarter
     2			- oNE(i,jp1,k) * half 
     2			- uNE(i,jp1,k) * quarter) 
c
c
		tmp8_XoN =
     2		- quarter * (- uE(im1,jp1,k) * quarter
     2			- oE(im1,jp1,kp1) * eighth
     2			- uSE(im1,jp2,k) * half
     2			- oNW(i,jp1,kp1) * quarter 
     2			- uC(i,jp1,k) * half
     2			+ oC(i,jp1,kp1) * quarter 
     2			- uS(i,jp2,k) * whole
     2			- oN(i,jp1,kp1) * half
     2			- uW(ip1,jp1,k) * quarter
     2			- oE(i,jp1,kp1) * eighth
     2			- uSW(ip1,jp2,k) * half
     2			- oNE(i,jp1,kp1) * quarter) 
c
     2		- quarter * (- oNW(ip1,j,km1) * quarter
     2			- uNW(ip1,j,km1) * half
     2			- oN(ip1,j,km1) * eighth
     2			- uN(ip1,j,km1) * quarter) 
c
c
		tmp9_XoN =
     2		- half * (- uSE(i,jp1,km1) * quarter 
     2			- oNW(ip1,j,k) * half
     2			- uNW(ip1,j,k) * quarter
     2			- uS(ip1,jp1,km1) * eighth
     2			- oN(ip1,j,k) * quarter 
     2			- uN(ip1,j,k) * eighth) 
c
     2		- quarter * (- uSE(i,jp1,k) * half
     2			- oNW(ip1,j,kp1) * quarter
     2			- uS(ip1,jp1,k) * quarter
     2			- oN(ip1,j,kp1) * eighth) 
c
c
		tmp10_XoN =
     2		- eighth * (- oE(i,jp1,km1) * quarter
     2			 - uW(ip1,jp1,km1) * half
     2			 - oNW(ip1,jp1,km1) * half
     2			 - uNW(ip1,jp1,km1) * whole
     2			 + oC(ip1,jp1,km1) * eighth
     2			 - uC(ip1,jp1,km1) * quarter
     2			 - oN(ip1,jp1,km1) * quarter
     2			 - uN(ip1,jp1,km1) * half)
c
c
		tmp11_XoN =
     2		- quarter * (- uE(i,jp1,km1) * quarter 
     2			  - oE(i,jp1,k) * half
     2			  - uW(ip1,jp1,k) * quarter
     2			  - uSE(i,jp2,km1) * half
     2			  - oNW(ip1,jp1,k) * whole
     2			  - uNW(ip1,jp1,k) * half
     2			  - uC(ip1,jp1,km1) * eighth
     2			  + oC(ip1,jp1,k) * quarter
     2			  - uC(ip1,jp1,k) * eighth
     2			  - uS(ip1,jp2,km1) * quarter
     2			  - oN(ip1,jp1,k) * half
     2			  - uN(ip1,jp1,k) * quarter) 
c
c
		tmp12_XoN =
     2		- eighth * (- uE(i,jp1,k) * half 
     2			 - oE(i,jp1,kp1) * quarter
     2			 - uSE(i,jp2,k) * whole
     2			 - oNW(ip1,jp1,kp1) * half
     2			 - uC(ip1,jp1,k) * quarter
     2			 + oC(ip1,jp1,kp1) * eighth
     2			 - uS(ip1,jp2,k) * half
     2			 - oN(ip1,jp1,kp1) * quarter)
c
		XoN(ii,jj,kk) = tmp1_XoN + tmp2_XoN + tmp3_XoN +tmp4_XoN
     2		+ tmp5_XoN + tmp6_XoN + tmp7_XoN + tmp8_XoN + tmp9_XoN
     3		+ tmp10_XoN + tmp11_XoN + tmp12_XoN
c
		if(jj.eq.nym1) XoN(ii,jj,kk) = 0.d0
c
c
c		XuC(ii,jj,kk)
c		=============================================================================
c
		tmp1_XuC =
     2		- quarter * (- uC(im1,jm1,k) * eighth
     2			 - uN(im1,jm1,k) * quarter
     2			 - uE(im1,jm1,k) * quarter
     2			 - uNE(im1,jm1,k) * half) 
c
     2		- eighth * (  oC(im1,jm1,kp1) * eighth
     2			 - uC(im1,jm1,kp1) * quarter
     2			 - oN(im1,jm1,kp1) * quarter
     2			 - uN(im1,jm1,kp1) * half
     2			 - oE(im1,jm1,kp1) * quarter
     2			 - uE(im1,jm1,kp1) * half
     2			 - oNE(im1,jm1,kp1) * half
     2			 - uNE(im1,jm1,kp1) * whole) 
c
c
		tmp2_XuC =
     2		- half * (- uS(im1,j,k) * eighth 
     2			- uC(im1,j,k) * quarter
     2			- uN(im1,j,k) * eighth 
     2			- uSE(im1,j,k) * quarter
     2			- uE(im1,j,k) * half 
     2			- uNE(im1,j,k) * quarter) 
c
c
		tmp3_XuC =
     2		- quarter * (- oN(im1,jm1,kp1) * eighth
     2			- uS(im1,j,kp1) * quarter
     2			+ oC(im1,j,kp1) * quarter
     2			- uC(im1,j,kp1) * half
     2			- oN(im1,j,kp1) * eighth
     2			- uN(im1,j,kp1) * quarter
     2			- oNW(i,jm1,kp1) * quarter
     2			- uSE(im1,j,kp1) * half
     2			- oE(im1,j,kp1) * half
     2			- uE(im1,j,kp1) * whole
     2			- oNE(im1,j,kp1) * quarter
     2			- uNE(im1,j,kp1) * half) 
c
     2		- quarter * (- uS(im1,jp1,k) * quarter
     2			 - uC(im1,jp1,k) * eighth
     2			 - uSE(im1,jp1,k) * half
     2			 - uE(im1,jp1,k) * quarter) 
c
c
		tmp4_XuC =
     2		- eighth * (- oN(im1,j,kp1) * quarter
     2			 - uS(im1,jp1,kp1) * half
     2			 + oC(im1,jp1,kp1) * eighth
     2			 - uC(im1,jp1,kp1) * quarter
     2			 - oNW(i,j,kp1) * half
     2			 - uSE(im1,jp1,kp1) * whole
     2			 - oE(im1,jp1,kp1) * quarter
     2			 - uE(im1,jp1,kp1) * half)
c
     2		- half * (- uW(i,jm1,k) * eighth 
     2			- uNW(i,jm1,k) * quarter
     2			- uC(i,jm1,k) * quarter 
     2			- uN(i,jm1,k) * half
     2			- uE(i,jm1,k) * eighth 
     2			- uNE(i,jm1,k) * quarter) 
c
c
		tmp5_XuC =
     2		- quarter * (- oE(im1,jm1,kp1) * eighth
     2			- uW(i,jm1,kp1) * quarter
     2			- oNW(i,jm1,kp1) * quarter
     2			- uNW(i,jm1,kp1) * half
     2			+ oC(i,jm1,kp1) * quarter
     2			- uC(i,jm1,kp1) * half
     2			- oN(i,jm1,kp1) * half
     2			- uN(i,jm1,kp1) * whole
     2			- oE(i,jm1,kp1) * eighth
     2			- uE(i,jm1,kp1) * quarter
     2			- oNE(i,jm1,kp1) * quarter
     2			- uNE(i,jm1,kp1) * half) 
c
c
		tmp6_XuC =
     2		- whole * (- uSW(i,j,k) * eighth 
     2			- uW(i,j,k) * quarter
     2			- uNW(i,j,k) * eighth 
     2			- uS(i,j,k) * quarter
     2			- uC(i,j,k) * half 
     2			- uN(i,j,k) * quarter
     2			- uSE(i,j,k) * eighth 
     2			- uE(i,j,k) * quarter
     2			- uNE(i,j,k) * eighth) 
c
c
		tmp7_XuC =
     2		- half * (- oNE(im1,jm1,kp1) * eighth
     2			- uSW(i,j,kp1) * quarter
     2			- oE(im1,j,kp1) * quarter 
     2			- uW(i,j,kp1) * half
     2			- oNW(i,j,kp1) * eighth 
     2			- uNW(i,j,kp1) * quarter
     2			- oN(i,jm1,kp1) * quarter 
     2			- uS(i,j,kp1) * half
     2			+ oC(i,j,kp1) * half 
     2			- uC(i,j,kp1) * whole
     2			- oN(i,j,kp1) * quarter 
     2			- uN(i,j,kp1) * half
     2			- oNW(ip1,jm1,kp1) * eighth
     2			- uSE(i,j,kp1) * quarter 
     2			- oE(i,j,kp1) * quarter
     2			- uE(i,j,kp1) * half 
     2			- oNE(i,j,kp1) * eighth
     2			- uNE(i,j,kp1) * quarter) 
c
c
		tmp8_XuC =
     2		- half * (- uSW(i,jp1,k) * quarter 
     2			- uW(i,jp1,k) * eighth
     2			- uS(i,jp1,k) * half 
     2			- uC(i,jp1,k) * quarter
     2			- uSE(i,jp1,k) * quarter 
     2			- uE(i,jp1,k) * eighth) 
c
c
		tmp9_XuC =
     2		- quarter * (- oNE(im1,j,kp1) * quarter
     2			- uSW(i,jp1,kp1) * half
     2			- oE(im1,jp1,kp1) * eighth
     2			- uW(i,jp1,kp1) * quarter 
     2			- oN(i,j,kp1) * half
     2			- uS(i,jp1,kp1) * whole
     2			+ oC(i,jp1,kp1) * quarter
     2			- uC(i,jp1,kp1) * half
     2			- oNW(ip1,j,kp1) * quarter
     2			- uSE(i,jp1,kp1) * half
     2			- oE(i,jp1,kp1) * eighth
     2			- uE(i,jp1,kp1) * quarter) 
c
     2		- quarter * (- uW(ip1,jm1,k) * quarter
     2			 - uNW(ip1,jm1,k) * half
     2			 - uC(ip1,jm1,k) * eighth
     2			 - uN(ip1,jm1,k) * quarter) 
c
c
		tmp10_XuC =
     2		- eighth * (- oE(i,jm1,kp1) * quarter
     2			 - uW(ip1,jm1,kp1) * half
     2			 - oNW(ip1,jm1,kp1) * half
     2			 - uNW(ip1,jm1,kp1) * whole
     2			 + oC(ip1,jm1,kp1) * eighth
     2			 - uC(ip1,jm1,kp1) * quarter
     2			 - oN(ip1,jm1,kp1) * quarter
     2			 - uN(ip1,jm1,kp1) * half)
c
     2		- half * (- uSW(ip1,j,k) * quarter 
     2			- uW(ip1,j,k) * half
     2			- uNW(ip1,j,k) * quarter 
     2			- uS(ip1,j,k) * eighth
     2			- uC(ip1,j,k) * quarter 
     2			- uN(ip1,j,k) * eighth) 
c
c
		tmp11_XuC =
     2		- quarter * (- oNE(i,jm1,kp1) * quarter
     2			- uSW(ip1,j,kp1) * half 
     2			- oE(i,j,kp1) * half
     2			- uW(ip1,j,kp1) * whole
     2			- oNW(ip1,j,kp1) * quarter
     2			- uNW(ip1,j,kp1) * half
     2			- oN(ip1,jm1,kp1) * eighth
     2			- uS(ip1,j,kp1) * quarter
     2			+ oC(ip1,j,kp1) * quarter
     2			- uC(ip1,j,kp1) * half
     2			- oN(ip1,j,kp1) * eighth
     2			- uN(ip1,j,kp1) * quarter) 
c
c
		tmp12_XuC =
     2		- quarter * (- uSW(ip1,jp1,k) * half
     2			 - uW(ip1,jp1,k) * quarter
     2			 - uS(ip1,jp1,k) * quarter
     2			 - uC(ip1,jp1,k) * eighth)
c
     2		- eighth * (- oNE(i,j,kp1) * half
     2			 - uSW(ip1,jp1,kp1) * whole
     2			 - oE(i,jp1,kp1) * quarter
     2			 - uW(ip1,jp1,kp1) * half
     2			 - oN(ip1,j,kp1) * quarter
     2			 - uS(ip1,jp1,kp1) * half
     2			 + oC(ip1,jp1,kp1) * eighth
     2			 - uC(ip1,jp1,kp1) * quarter)
c
		XuC(ii,jj,kk) = tmp1_XuC + tmp2_XuC + tmp3_XuC +tmp4_XuC
     2		+ tmp5_XuC + tmp6_XuC + tmp7_XuC + tmp8_XuC + tmp9_XuC
     3		+ tmp10_XuC + tmp11_XuC + tmp12_XuC
c
		if(kk.eq.nzm1) XuC(ii,jj,kk) = 0.d0
c
c
c		XoNE(ii,jj,kk)
c		=============================================================================
c
		tmp1_XoNE =
     2		- half * (- oNE(i,j,km1) * eighth
     2			- uNE(i,j,km1) * quarter) 
c
     2		- whole * (- uSW(ip1,jp1,km1)*eighth
     2			- oNE(i,j,k) * quarter 
     2			- uNE(i,j,k) * eighth)
c
     2		- half * (- uSW(ip1,jp1,k) * quarter
     2			- oNE(i,j,kp1) * eighth) 
c
     2		- quarter * (- oE(i,jp1,km1) * eighth
     2			- uE(i,jp1,km1) * quarter
     2			- oNE(i,jp1,km1) * quarter
     2			- uNE(i,jp1,km1) * half)
c
c
		tmp2_XoNE =
     2		- half * (- uW(ip1,jp1,km1)*eighth
     2			- oE(i,jp1,k) * quarter
     2			- uE(i,jp1,k) * eighth
     2			- uSW(ip1,jp2,km1) * quarter
     2			- oNE(i,jp1,k) * half
     2			- uNE(i,jp1,k) * quarter) 
c
     2		- quarter * (- uW(ip1,jp1,k) * quarter
     2			- oE(i,jp1,kp1) * eighth
     2			- uSW(ip1,jp2,k) * half
     2			- oNE(i,jp1,kp1) * quarter)
c
     2		- quarter * (- oN(ip1,j,km1) * eighth
     2			 - uN(ip1,j,km1) * quarter
     2			 - oNE(ip1,j,km1) * quarter
     2			 - uNE(ip1,j,km1) * half)
c
c
		tmp3_XoNE =
     2		- half * (- uS(ip1,jp1,km1)*eighth
     2			 - oN(ip1,j,k) * quarter
     2			 - uN(ip1,j,k) * eighth
     2			 - uSW(ip2,jp1,km1) * quarter
     2			 - oNE(ip1,j,k) * half
     2			 - uNE(ip1,j,k) * quarter) 
c
     2		- quarter * (- uS(ip1,jp1,k) * quarter
     2			- oN(ip1,j,kp1) * eighth
     2			- uSW(ip2,jp1,k) * half
     2			- oNE(ip1,j,kp1) * quarter)
c
c
		tmp4_XoNE =
     2		- eighth * (  oC(ip1,jp1,km1)*eighth
     2			  - uC(ip1,jp1,km1) * quarter
     2			  - oN(ip1,jp1,km1) * quarter
     2			  - uN(ip1,jp1,km1) * half
     2			  - oE(ip1,jp1,km1) * quarter
     2			  - uE(ip1,jp1,km1) * half
     2			  - oNE(ip1,jp1,km1) * half
     2			  - uNE(ip1,jp1,km1) * whole) 
c
c
		tmp5_XoNE =
     2		- quarter * (- uC(ip1,jp1,km1)*eighth
     2			 + oC(ip1,jp1,k) * quarter
     2			 - uC(ip1,jp1,k) * eighth
     2			 - uS(ip1,jp2,km1) * quarter
     2			 - oN(ip1,jp1,k) * half
     2			 - uN(ip1,jp1,k) * quarter
     2			 - uW(ip2,jp1,km1) * quarter
     2			 - oE(ip1,jp1,k) * half
     2			 - uE(ip1,jp1,k) * quarter
     2			 - uSW(ip2,jp2,km1) * half
     2			 - oNE(ip1,jp1,k) * whole
     2			 - uNE(ip1,jp1,k) * half) 
c
c
		tmp6_XoNE =
     2		- eighth * (- uC(ip1,jp1,k) * quarter
     2			 + oC(ip1,jp1,kp1) * eighth
     2			 - uS(ip1,jp2,k) * half
     2			 - oN(ip1,jp1,kp1) * quarter
     2			 - uW(ip2,jp1,k) * half
     2			 - oE(ip1,jp1,kp1) * quarter
     2			 - uSW(ip2,jp2,k) * whole
     2			 - oNE(ip1,jp1,kp1) * half)
c
		XoNE(ii,jj,kk) = tmp1_XoNE + tmp2_XoNE + tmp3_XoNE 
     2			+ tmp4_XoNE + tmp5_XoNE + tmp6_XoNE
c
		if(ii.eq.nxm1.or.jj.eq.nym1) XoNE(ii,jj,kk) = 0.d0
c
c
c		XoNW(ii,jj,kk)
c		=============================================================================
c
		tmp1_XoNW =
     2		- quarter * (- oNW(im1,j,km1) * quarter
     2			- uNW(im1,j,km1) * half
     2			- oN(im1,j,km1) * eighth
     2			- uN(im1,j,km1) * quarter) 
c
     2		- half * (- uSE(im2,jp1,km1)*quarter
     2			- oNW(im1,j,k) * half
     2			- uNW(im1,j,k) * quarter
     2			- uS(im1,jp1,km1) * eighth
     2			- oN(im1,j,k) * quarter
     2			- uN(im1,j,k) * eighth) 
c
c
		tmp2_XoNW =
     2		- quarter * (- uSE(im2,jp1,k) * half
     2			- oNW(im1,j,kp1) * quarter
     2			- uS(im1,jp1,k) * quarter
     2			- oN(im1,j,kp1) * eighth) 
c
     2		- eighth * (- oE(im2,jp1,km1)*quarter
     2			 - uW(im1,jp1,km1) * half
     2			 - oNW(im1,jp1,km1) * half
     2			 - uNW(im1,jp1,km1) * whole
     2			 + oC(im1,jp1,km1) * eighth
     2			 - uC(im1,jp1,km1) * quarter
     2			 - oN(im1,jp1,km1) * quarter
     2			 - uN(im1,jp1,km1) * half) 
c
c
		tmp3_XoNW =
     2		- quarter * (- uE(im2,jp1,km1)*quarter
     2			 - oE(im2,jp1,k) * half
     2			 - uW(im1,jp1,k) * quarter
     2			 - uSE(im2,jp2,km1) * half
     2			 - oNW(im1,jp1,k) * whole
     2			 - uNW(im1,jp1,k) * half
     2			 - uC(im1,jp1,km1) * eighth
     2			 + oC(im1,jp1,k) * quarter
     2			 - uC(im1,jp1,k) * eighth
     2			 - uS(im1,jp2,km1) * quarter
     2			 - oN(im1,jp1,k) * half
     2			 - uN(im1,jp1,k) * quarter) 
c
c
		tmp4_XoNW =
     2		- eighth * (- uE(im2,jp1,k) * half
     2			 - oE(im2,jp1,kp1) * quarter
     2			 - uSE(im2,jp2,k) * whole
     2			 - oNW(im1,jp1,kp1) * half
     2			 - uC(im1,jp1,k) * quarter
     2			 + oC(im1,jp1,kp1) * eighth
     2			 - uS(im1,jp2,k) * half
     2			 - oN(im1,jp1,kp1) * quarter) 
c
     2		- half * (- oNW(i,j,km1) * eighth
     2			- uNW(i,j,km1) * quarter) 
c
c
		tmp5_XoNW =
     2		- whole * (- uSE(im1,jp1,km1)*eighth
     2			- oNW(i,j,k) * quarter 
     2			- uNW(i,j,k) * eighth)
c
     2		- half * (- uSE(im1,jp1,k) * quarter
     2			- oNW(i,j,kp1) * eighth) 
c
     2		- quarter * (- oE(im1,jp1,km1)*eighth
     2			- uW(i,jp1,km1) * quarter
     2			- oNW(i,jp1,km1) * quarter
     2			- uNW(i,jp1,km1) * half) 
c
c
		tmp6_XoNW =
     2		- half * (- uE(im1,jp1,km1)*eighth
     2			- oE(im1,jp1,k) * quarter
     2			- uW(i,jp1,k) * eighth
     2			- uSE(im1,jp2,km1) * quarter
     2			- oNW(i,jp1,k) * half
     2			- uNW(i,jp1,k) * quarter) 
c
     2		- quarter * (- uE(im1,jp1,k) * quarter
     2			- oE(im1,jp1,kp1) * eighth
     2			- uSE(im1,jp2,k) * half
     2			- oNW(i,jp1,kp1) * quarter)
c
		XoNW(ii,jj,kk) = tmp1_XoNW+tmp2_XoNW+tmp3_XoNW+tmp4_XoNW
     2				+ tmp5_XoNW + tmp6_XoNW
c
		if(ii.eq.2) XoNW(ii,jj,kk) = 0
		if(jj.eq.nym1) XoNW(ii,jj,kk) = 0
c
c
c		XuE(ii,jj,kk)
c		=============================================================================
c
		tmp1_XuE =
     2		- half * (- uE(i,jm1,k) * eighth 
     2			- uNE(i,jm1,k) * quarter) 
c
     2		- quarter * (- oE(i,jm1,kp1) * eighth
     2			- uE(i,jm1,kp1) * quarter
     2			- oNE(i,jm1,kp1) * quarter
     2			- uNE(i,jm1,kp1) * half) 
c
     2		- whole * (- uSE(i,j,k) * eighth 
     2			- uE(i,j,k) * quarter
     2			- uNE(i,j,k) * eighth) 
c
c
		tmp2_XuE =
     2		- half * (- oNW(ip1,jm1,kp1)*eighth
     2			- uSE(i,j,kp1) * quarter
     2			- oE(i,j,kp1) * quarter 
     2			- uE(i,j,kp1) * half
     2			- oNE(i,j,kp1) * eighth
     2			- uNE(i,j,kp1) * quarter) 
c
     2		- half * (- uSE(i,jp1,k) * quarter 
     2			- uE(i,jp1,k) * eighth) 
c
     2		- quarter * (- oNW(ip1,j,kp1) * quarter
     2			- uSE(i,jp1,kp1) * half
     2			- oE(i,jp1,kp1) * eighth
     2			- uE(i,jp1,kp1) * quarter) 
c
c
		tmp3_XuE =
     2		- quarter * (- uC(ip1,jm1,k) * eighth
     2			 - uN(ip1,jm1,k) * quarter
     2			 - uE(ip1,jm1,k) * quarter
     2			 - uNE(ip1,jm1,k) * half) 
c
     2		- eighth * (  oC(ip1,jm1,kp1)*eighth
     2			 - uC(ip1,jm1,kp1) * quarter
     2			 - oN(ip1,jm1,kp1) * quarter
     2			 - uN(ip1,jm1,kp1) * half
     2			 - oE(ip1,jm1,kp1) * quarter
     2			 - uE(ip1,jm1,kp1) * half
     2			 - oNE(ip1,jm1,kp1) * half
     2			 - uNE(ip1,jm1,kp1) * whole) 
c
c
		tmp4_XuE =
     2		- half * (- uS(ip1,j,k) * eighth 
     2			- uC(ip1,j,k) * quarter
     2			- uN(ip1,j,k) * eighth
     2			- uSE(ip1,j,k) * quarter 
     2			- uE(ip1,j,k) * half
     2			- uNE(ip1,j,k) * quarter) 
c
c
		tmp5_XuE =
     2		- quarter * (- oN(ip1,jm1,kp1)*eighth
     2			- uS(ip1,j,kp1) * quarter
     2			+ oC(ip1,j,kp1) * quarter
     2			- uC(ip1,j,kp1) * half
     2			- oN(ip1,j,kp1) * eighth
     2			- uN(ip1,j,kp1) * quarter
     2			- oNW(ip2,jm1,kp1) * quarter
     2			- uSE(ip1,j,kp1) * half
     2			- oE(ip1,j,kp1) * half
     2			- uE(ip1,j,kp1) * whole
     2			- oNE(ip1,j,kp1) * quarter
     2			- uNE(ip1,j,kp1) * half) 
c
c
		tmp6_XuE =
     2		- quarter * (- uS(ip1,jp1,k) * quarter
     2			 - uC(ip1,jp1,k) * eighth
     2			 - uSE(ip1,jp1,k) * half
     2			 - uE(ip1,jp1,k) * quarter) 
c
     2		- eighth * (- oN(ip1,j,kp1) * quarter
     2			 - uS(ip1,jp1,kp1) * half
     2			 + oC(ip1,jp1,kp1) * eighth
     2			 - uC(ip1,jp1,kp1) * quarter
     2			 - oNW(ip2,j,kp1) * half
     2			 - uSE(ip1,jp1,kp1) * whole
     2			 - oE(ip1,jp1,kp1) * quarter
     2			 - uE(ip1,jp1,kp1) * half)
c
		XuE(ii,jj,kk) = tmp1_XuE + tmp2_XuE + tmp3_XuE +tmp4_XuE
     2				+ tmp5_XuE + tmp6_XuE
c
		if(ii.eq.nxm1.or.kk.eq.nzm1) XuE(ii,jj,kk) = 0.d0
c
c
c		XuW(ii,jj,kk)
c		=============================================================================
c
		tmp1_XuW =
     2		- quarter * (- uW(im1,jm1,k) * quarter
     2			 - uNW(im1,jm1,k) * half
     2			 - uC(im1,jm1,k) * eighth
     2			 - uN(im1,jm1,k) * quarter) 
c
     2		- eighth * (- oE(im2,jm1,kp1) *quarter
     2			 - uW(im1,jm1,kp1) * half
     2			 - oNW(im1,jm1,kp1) * half
     2			 - uNW(im1,jm1,kp1) * whole
     2			 + oC(im1,jm1,kp1) * eighth
     2			 - uC(im1,jm1,kp1) * quarter
     2			 - oN(im1,jm1,kp1) * quarter
     2			 - uN(im1,jm1,kp1) * half) 
c
c
		tmp2_XuW =
     2		- half * (- uSW(im1,j,k) * quarter 
     2			- uW(im1,j,k) * half
     2			- uNW(im1,j,k) * quarter
     2			- uS(im1,j,k) * eighth 
     2			- uC(im1,j,k) * quarter
     2			- uN(im1,j,k) * eighth) 
c
c
		tmp3_XuW =
     2		- quarter * (- oNE(im2,jm1,kp1)*quarter
     2			- uSW(im1,j,kp1) * half
     2			- oE(im2,j,kp1) * half
     2			- uW(im1,j,kp1) * whole
     2			- oNW(im1,j,kp1) * quarter
     2			- uNW(im1,j,kp1) * half
     2			- oN(im1,jm1,kp1) * eighth
     2			- uS(im1,j,kp1) * quarter
     2			+ oC(im1,j,kp1) * quarter
     2			- uC(im1,j,kp1) * half
     2			- oN(im1,j,kp1) * eighth
     2			- uN(im1,j,kp1) * quarter) 
c
c
		tmp4_XuW =
     2		- quarter * (- uSW(im1,jp1,k) * half
     2			  - uW(im1,jp1,k) * quarter
     2			  - uS(im1,jp1,k) * quarter
     2			  - uC(im1,jp1,k) * eighth) 
c
     2		- eighth * (- oNE(im2,j,kp1) * half
     2			  - uSW(im1,jp1,kp1) * whole
     2			  - oE(im2,jp1,kp1) * quarter
     2			  - uW(im1,jp1,kp1) * half
     2			  - oN(im1,j,kp1) * quarter
     2			  - uS(im1,jp1,kp1) * half
     2			  + oC(im1,jp1,kp1) * eighth
     2			  - uC(im1,jp1,kp1) * quarter) 
c
c
		tmp5_XuW =
     2		- half * (- uW(i,jm1,k) * eighth 
     2			- uNW(i,jm1,k) * quarter) 
c
     2		- quarter * (- oE(im1,jm1,kp1)*eighth
     2			- uW(i,jm1,kp1) * quarter
     2			- oNW(i,jm1,kp1) * quarter
     2			- uNW(i,jm1,kp1) * half) 
c
     2		- whole * (- uSW(i,j,k) * eighth 
     2			- uW(i,j,k) * quarter
     2			- uNW(i,j,k) * eighth) 
c
c
		tmp6_XuW =
     2		- half * (- oNE(im1,jm1,kp1)*eighth
     2			- uSW(i,j,kp1) * quarter
     2			- oE(im1,j,kp1) * quarter
     2			- uW(i,j,kp1) * half
     2			- oNW(i,j,kp1) * eighth
     2			- uNW(i,j,kp1) * quarter) 
c
     2		- half * (- uSW(i,jp1,k) * quarter 
     2			- uW(i,jp1,k) * eighth) 
c
     2		- quarter * (- oNE(im1,j,kp1) * quarter
     2			- uSW(i,jp1,kp1) * half
     2			- oE(im1,jp1,kp1) * eighth
     2			- uW(i,jp1,kp1) * quarter)
c
		XuW(ii,jj,kk) = tmp1_XuW + tmp2_XuW + tmp3_XuW+tmp4_XuW
     2				+ tmp5_XuW + tmp6_XuW
c
		if(ii.eq.2) XuW(ii,jj,kk) = 0
		if(kk.eq.nzm1) XuW(ii,jj,kk) = 0.d0
c
c
c		XuN(ii,jj,kk)
c		=============================================================================
c
		tmp1_XuN =
     2		- half * (- uN(im1,j,k) * eighth 
     2			- uNE(im1,j,k) * quarter) 
c
     2		- quarter * (- oN(im1,j,kp1) * eighth
     2			- uN(im1,j,kp1) * quarter
     2			- oNE(im1,j,kp1) * quarter
     2			- uNE(im1,j,kp1) * half) 
c
     2		- quarter * (- uC(im1,jp1,k) * eighth
     2			 - uN(im1,jp1,k) * quarter
     2			 - uE(im1,jp1,k) * quarter
     2			 - uNE(im1,jp1,k) * half) 
c
c
		tmp2_XuN =
     2		- eighth * (  oC(im1,jp1,kp1)*eighth
     2			 - uC(im1,jp1,kp1) * quarter
     2			 - oN(im1,jp1,kp1) * quarter
     2			 - uN(im1,jp1,kp1) * half
     2			 - oE(im1,jp1,kp1) * quarter
     2			 - uE(im1,jp1,kp1) * half
     2			 - oNE(im1,jp1,kp1) * half
     2			 - uNE(im1,jp1,kp1) * whole)
c
     2		- whole * (- uNW(i,j,k) * eighth 
     2			 - uN(i,j,k) * quarter
     2			 - uNE(i,j,k) * eighth) 
c
c
		tmp3_XuN =
     2		- half * (- oNW(i,j,kp1) * eighth
     2			- uNW(i,j,kp1) * quarter
     2			- oN(i,j,kp1) * quarter 
     2			- uN(i,j,kp1) * half
     2			- oNE(i,j,kp1) * eighth
     2			- uNE(i,j,kp1) * quarter) 
c
     2		- half * (- uW(i,jp1,k) * eighth 
     2			- uNW(i,jp1,k) * quarter
     2			- uC(i,jp1,k) * quarter 
     2			- uN(i,jp1,k) * half
     2			- uE(i,jp1,k) * eighth
     2			- uNE(i,jp1,k) * quarter) 
c
c
		tmp4_XuN =
     2		- quarter * (- oE(im1,jp1,kp1)*eighth
     2			- uW(i,jp1,kp1) * quarter
     2			- oNW(i,jp1,kp1) * quarter
     2			- uNW(i,jp1,kp1) * half
     2			+ oC(i,jp1,kp1) * quarter
     2			- uC(i,jp1,kp1) * half
     2			- oN(i,jp1,kp1) * half
     2			- uN(i,jp1,kp1) * whole
     2			- oE(i,jp1,kp1) * eighth
     2			- uE(i,jp1,kp1) * quarter
     2			- oNE(i,jp1,kp1) * quarter
     2			- uNE(i,jp1,kp1) * half) 
c
c
		tmp5_XuN =
     2		- half * (- uNW(ip1,j,k) * quarter 
     2			- uN(ip1,j,k) * eighth) 
c
     2		- quarter * (- oNW(ip1,j,kp1) * quarter
     2			- uNW(ip1,j,kp1) * half
     2			- oN(ip1,j,kp1) * eighth
     2			- uN(ip1,j,kp1) * quarter) 
c
     2		- quarter * (- uW(ip1,jp1,k) * quarter
     2			 - uNW(ip1,jp1,k) * half
     2			 - uC(ip1,jp1,k) * eighth
     2			 - uN(ip1,jp1,k) * quarter) 
c
c
		tmp6_XuN =
     2		- eighth * (- oE(i,jp1,kp1) * quarter
     2			 - uW(ip1,jp1,kp1) * half
     2			 - oNW(ip1,jp1,kp1) * half
     2			 - uNW(ip1,jp1,kp1) * whole
     2			 + oC(ip1,jp1,kp1) * eighth
     2			 - uC(ip1,jp1,kp1) * quarter
     2			 - oN(ip1,jp1,kp1) * quarter
     2			 - uN(ip1,jp1,kp1) * half)
c
		XuN(ii,jj,kk) = tmp1_XuN + tmp2_XuN + tmp3_XuN+tmp4_XuN
     2				+ tmp5_XuN + tmp6_XuN
c
		if(jj.eq.nym1.or.kk.eq.nzm1) XuN(ii,jj,kk) = 0.d0
c
c
c		XuS(ii,jj,kk)
c		=============================================================================
c
		tmp1_XuS =
     2		- quarter * (- uS(im1,jm1,k) * quarter
     2			 - uC(im1,jm1,k) * eighth
     2			 - uSE(im1,jm1,k) * half
     2			 - uE(im1,jm1,k) * quarter) 
c
     2		- eighth * (- oN(im1,jm2,kp1)*quarter
     2			 - uS(im1,jm1,kp1) * half
     2			 + oC(im1,jm1,kp1) * eighth
     2			 - uC(im1,jm1,kp1) * quarter
     2			 - oNW(i,jm2,kp1) * half
     2			 - uSE(im1,jm1,kp1) * whole
     2			 - oE(im1,jm1,kp1) * quarter
     2			 - uE(im1,jm1,kp1) * half)
c
c
		tmp2_XuS =
     2		- half * (- uS(im1,j,k) * eighth 
     2			- uSE(im1,j,k) * quarter) 
c
     2		- quarter * (- oN(im1,jm1,kp1)*eighth
     2			- uS(im1,j,kp1) * quarter
     2			- oNW(i,jm1,kp1) * quarter
     2			- uSE(im1,j,kp1) * half) 
c
     2		- half * (- uSW(i,jm1,k) * quarter 
     2			- uW(i,jm1,k) * eighth
     2			- uS(i,jm1,k) * half 
     2			- uC(i,jm1,k) * quarter
     2			- uSE(i,jm1,k) * quarter
     2			- uE(i,jm1,k) * eighth) 
c
c
		tmp3_XuS =
     2		- quarter * (- oNE(im1,jm2,kp1)*quarter
     2			- uSW(i,jm1,kp1) * half
     2			- oE(im1,jm1,kp1) * eighth
     2			- uW(i,jm1,kp1) * quarter
     2			- oN(i,jm2,kp1) * half
     2			- uS(i,jm1,kp1) * whole
     2			+ oC(i,jm1,kp1) * quarter
     2			- uC(i,jm1,kp1) * half
     2			- oNW(ip1,jm2,kp1) * quarter
     2			- uSE(i,jm1,kp1) * half
     2			- oE(i,jm1,kp1) * eighth
     2			- uE(i,jm1,kp1) * quarter) 
c
c
		tmp4_XuS =
     2		- whole * (- uSW(i,j,k) * eighth 
     2			- uS(i,j,k) * quarter
     2			- uSE(i,j,k) * eighth) 

     2		- half * (- oNE(im1,jm1,kp1)*eighth
     2			- uSW(i,j,kp1) * quarter
     2			- oN(i,jm1,kp1) * quarter
     2			- uS(i,j,kp1) * half
     2			- oNW(ip1,jm1,kp1) * eighth
     2			- uSE(i,j,kp1) * quarter) 
c
c
		tmp5_XuS =
     2		- quarter * (- uSW(ip1,jm1,k) * half
     2			 - uW(ip1,jm1,k) * quarter
     2			 - uS(ip1,jm1,k) * quarter
     2			 - uC(ip1,jm1,k) * eighth) 
c
     2		- eighth * (- oNE(i,jm2,kp1) * half
     2			 - uSW(ip1,jm1,kp1) * whole
     2			 - oE(i,jm1,kp1) * quarter
     2			 - uW(ip1,jm1,kp1) * half
     2			 - oN(ip1,jm2,kp1) * quarter
     2			 - uS(ip1,jm1,kp1) * half
     2			 + oC(ip1,jm1,kp1) * eighth
     2			 - uC(ip1,jm1,kp1) * quarter)
c
c
		tmp6_XuS =
     2		- half * (- uSW(ip1,j,k) * quarter 
     2			 - uS(ip1,j,k) * eighth) 
c
     2		- quarter * (- oNE(i,jm1,kp1) * quarter
     2			 - uSW(ip1,j,kp1) * half
     2			 - oN(ip1,jm1,kp1) * eighth
     2			 - uS(ip1,j,kp1) * quarter)
c
		XuS(ii,jj,kk) = tmp1_XuS + tmp2_XuS + tmp3_XuS+tmp4_XuS
     2				+ tmp5_XuS + tmp6_XuS
c
		if(jj.eq.2) XuS(ii,jj,kk) = 0.d0
		if(kk.eq.nzm1) XuS(ii,jj,kk) = 0.d0
c
c
c		XuNE(ii,jj,kk)
c		=============================================================================
c
		tmp1_XuNE =
     2		+ whole * uNE(i,j,k) * eighth
c
     2		- half * (- oNE(i,j,kp1) * eighth
     2			- uNE(i,j,kp1) * quarter)
c
     2		- half * (- uE(i,jp1,k) * eighth
     2			- uNE(i,jp1,k) * quarter)
c
     2		- quarter * (- oE(i,jp1,kp1)*eighth
     2			- uE(i,jp1,kp1) * quarter
     2			- oNE(i,jp1,kp1) * quarter
     2			- uNE(i,jp1,kp1) * half)
c
     2		- half * (- uN(ip1,j,k) * eighth
     2			- uNE(ip1,j,k) * quarter)
c
c
		tmp2_XuNE =
     2		- quarter * (- oN(ip1,j,kp1)*eighth
     2			- uN(ip1,j,kp1) * quarter
     2			- oNE(ip1,j,kp1) * quarter
     2			- uNE(ip1,j,kp1) * half)
c
     2		- quarter * (- uC(ip1,jp1,k)*eighth
     2			 - uN(ip1,jp1,k) * quarter
     2			 - uE(ip1,jp1,k) * quarter
     2			 - uNE(ip1,jp1,k) * half)
c
     2		- eighth*(  oC(ip1,jp1,kp1)*eighth
     2			 - uC(ip1,jp1,kp1) * quarter
     2			 - oN(ip1,jp1,kp1) * quarter
     2			 - uN(ip1,jp1,kp1) * half
     2			 - oE(ip1,jp1,kp1) * quarter
     2			 - uE(ip1,jp1,kp1) * half
     2			 - oNE(ip1,jp1,kp1) * half
     2			 - uNE(ip1,jp1,kp1) * whole)
c
		XuNE(ii,jj,kk) = tmp1_XuNE + tmp2_XuNE
c
		if(ii.eq.nxm1.or.jj.eq.nym1.or.kk.eq.nzm1) 
     1			XuNE(ii,jj,kk) = 0.d0
c
c
c		XuNW(ii,jj,kk)
c		=============================================================================
c
		tmp1_XuNW =
     2		- half * (- uNW(im1,j,k) * quarter
     2			- uN(im1,j,k) * eighth)
c
     2		- quarter * (- oNW(im1,j,kp1)*quarter
     2			- uNW(im1,j,kp1) * half
     2			- oN(im1,j,kp1) * eighth
     2			- uN(im1,j,kp1) * quarter)
c
     2		- quarter * (- uW(im1,jp1,k)*quarter
     2			 - uNW(im1,jp1,k) * half
     2			 - uC(im1,jp1,k) * eighth
     2			 - uN(im1,jp1,k) * quarter)
c
c
		tmp2_XuNW =
     2		- eighth*(- oE(im2,jp1,kp1)*quarter
     2			 - uW(im1,jp1,kp1) * half
     2			 - oNW(im1,jp1,kp1) * half
     2			 - uNW(im1,jp1,kp1) * whole
     2			 + oC(im1,jp1,kp1) * eighth
     2			 - uC(im1,jp1,kp1) * quarter
     2			 - oN(im1,jp1,kp1) * quarter
     2			 - uN(im1,jp1,kp1) * half)
c
     2		+ whole * uNW(i,j,k) * eighth
c
     2		- half * (- oNW(i,j,kp1) * eighth
     2			- uNW(i,j,kp1) * quarter)
c
     2		- half * (- uW(i,jp1,k) * eighth
     2			- uNW(i,jp1,k) * quarter)
c
     2		- quarter*(- oE(im1,jp1,kp1)*eighth
     2			- uW(i,jp1,kp1) * quarter
     2			- oNW(i,jp1,kp1) * quarter
     2			- uNW(i,jp1,kp1) * half)
c
		XuNW(ii,jj,kk) = tmp1_XuNW + tmp2_XuNW
c
		if(ii.eq.2) XuNW(ii,jj,kk) = 0.d0
		if(jj.eq.nym1.or.kk.eq.nzm1) XuNW(ii,jj,kk) = 0.d0
c
c
c		XuSE(ii,jj,kk)
c		=============================================================================
c
		tmp1_XuSE =
     2		- half * (- uSE(i,jm1,k) * quarter
     2			- uE(i,jm1,k) * eighth)
c
     2		- quarter*(- oNW(ip1,jm2,kp1)*quarter
     2			- uSE(i,jm1,kp1) * half 
     2			- oE(i,jm1,kp1) * eighth
     2			- uE(i,jm1,kp1) * quarter)
c
     2		+ whole * uSE(i,j,k) * eighth

     2		- half*(- oNW(ip1,jm1,kp1)*eighth
     2			- uSE(i,j,kp1) * quarter)
c
c
		tmp2_XuSE =
     2		- quarter*(- uS(ip1,jm1,k) * quarter
     2			 - uC(ip1,jm1,k) * eighth
     2			 - uSE(ip1,jm1,k) * half 
     2			 - uE(ip1,jm1,k) * quarter)
c
     2		- eighth*(- oN(ip1,jm2,kp1)*quarter
     2			 - uS(ip1,jm1,kp1) * half
     2			 + oC(ip1,jm1,kp1) * eighth
     2			 - uC(ip1,jm1,kp1) * quarter
     2			 - oNW(ip2,jm2,kp1) * half
     2			 - uSE(ip1,jm1,kp1) * whole
     2			 - oE(ip1,jm1,kp1) * quarter
     2			 - uE(ip1,jm1,kp1) * half)
c
     2		- half * (- uS(ip1,j,k) * eighth
     2			- uSE(ip1,j,k) * quarter)
c
     2		- quarter*(- oN(ip1,jm1,kp1)*eighth
     2			- uS(ip1,j,kp1) * quarter
     2			- oNW(ip2,jm1,kp1) * quarter
     2			- uSE(ip1,j,kp1) * half)
c
		XuSE(ii,jj,kk) = tmp1_XuSE + tmp2_XuSE
c
		if(jj.eq.2) XuSE(ii,jj,kk) = 0
		if(ii.eq.nxm1.or.kk.eq.nzm1) XuSE(ii,jj,kk) = 0.d0
c
c
c		XuSW(ii,jj,kk)
c		=============================================================================
c
		tmp1_XuSW =
     2		- quarter * (- uSW(im1,jm1,k)*half
     2			 - uW(im1,jm1,k) * quarter
     2			 - uS(im1,jm1,k) * quarter
     2			 - uC(im1,jm1,k) * eighth)
c
     2		- eighth*(- oNE(im2,jm2,kp1)*half
     2			 - uSW(im1,jm1,kp1) * whole
     2			 - oE(im2,jm1,kp1) * quarter
     2			 - uW(im1,jm1,kp1) * half
     2			 - oN(im1,jm2,kp1) * quarter
     2			 - uS(im1,jm1,kp1) * half
     2			 + oC(im1,jm1,kp1) * eighth
     2			 - uC(im1,jm1,kp1) * quarter)
c
     2		- half * (- uSW(im1,j,k) * quarter
     2			- uS(im1,j,k) * eighth)
c
c
		tmp2_XuSW =
     2		- quarter*(- oNE(im2,jm1,kp1)*quarter
     2			- uSW(im1,j,kp1) * half
     2			- oN(im1,jm1,kp1) * eighth
     2			- uS(im1,j,kp1) * quarter)
c
     2		- half * (- uSW(i,jm1,k) * quarter
     2			- uW(i,jm1,k) * eighth)
c
     2		- quarter*(- oNE(im1,jm2,kp1)*quarter
     2			- uSW(i,jm1,kp1) * half
     2			- oE(im1,jm1,kp1) * eighth
     2			- uW(i,jm1,kp1) * quarter)
c
     2		+ whole * uSW(i,j,k) * eighth
c
     2		- half*(- oNE(im1,jm1,kp1)*eighth
     2			- uSW(i,j,kp1) * quarter)
c
		XuSW(ii,jj,kk) = tmp1_XuSW + tmp2_XuSW
c
		if(ii.eq.2.or.jj.eq.2) XuSW(ii,jj,kk) = 0.d0
		if(kk.eq.nzm1) XuSW(ii,jj,kk) = 0.d0
c
c
c		End 3 main loops
c
100	    continue
c
200	  continue
c
300	continue
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c	Update_Galerkin.f
c ==========================================================================================
c ==========================================================================================
c
c	Form the Galerkin coarse grid system
c
c	(note that the coarse grid is always 27 diagonal (with only 14 stored, as
c	the matrix is symmetric)
c
	subroutine update_Galerkin(nxf,nyf,nzf,nxc,nyc,nzc,pcF,ipcF,
     1				acF,ac)
c
c
c ==========================================================================================
c	Declare all variables
c ==========================================================================================
c
c	Input/Output
c ==========================================================================================
c
	integer	nxf,nyf,nzf,nxc,nyc,nzc
	integer	numdia_local
c
	integer	ipcF(*)
c
	real*8	pcF(*),acF(*)
	real*8	ac(*)
c
c ==========================================================================================
c	Call Galerkin routine
c ==========================================================================================
c
	numdia_local = ipcF(11)
	call build_G(nxf,nyf,nzf,nxc,nyc,nzc,numdia_local,acF,ac)
c
	return
	end
