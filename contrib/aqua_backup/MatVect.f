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
c	MatVect.f
c
c	Version 1:      7/4/07
c
c	Author: Patrice Koehl (in collaboration with Marc Delarue)
c
c	This file contains a set of routines for grid operations
c
c ==========================================================================================
c ==========================================================================================
c	Matvec_band.f
c
c	This subroutines computes a matrix-vector multiplication, when the matrix
c	is stored in diagonal (banded) format
c ==========================================================================================
c ==========================================================================================
c
	subroutine matvec_band(nx,ny,nz,ipc,rpc,ac,cc,x,y)
c
c ==========================================================================================
c	Declaration 
c ==========================================================================================
c
	integer	nx,ny,nz,numdia
c
	integer	ipc(*)
c
	real*8	rpc(*),ac(nx*ny*nz,*),cc(*),x(*),y(*)
c
c ==========================================================================================
c	Call proper routine based on number of diagonals
c ==========================================================================================
c
	numdia = ipc(11)
c
	if(numdia.eq.7) then
		call matvec_band7(nx,ny,nz,ipc,rpc,ac(1,1),cc,
     1			ac(1,2),ac(1,3),ac(1,4),
     2			x,y)
	elseif(numdia.eq.27) then
		call matvec_band27(nx,ny,nz,ipc,rpc,ac(1,1),cc,
     1		ac(1,2),ac(1,3),ac(1,4),ac(1,5),ac(1,6),ac(1,7),ac(1,8),
     2		ac(1,9),ac(1,10),ac(1,11),ac(1,12),ac(1,13),ac(1,14),
     3		x,y)
	else
		write(6,*) 'Matrix-vector: bad stencil type...'
		stop
	endif
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c	Matvec_band7.f
c
c	This subroutines computes a matrix-vector multiplication, when the matrix
c	is stored with 7-diagonal (banded) format
c ==========================================================================================
c ==========================================================================================
c
	subroutine matvec_band7(nx,ny,nz,ipc,rpc,oC,cc,oE,oN,uC,x,y)
c
c ==========================================================================================
c	Declaration 
c ==========================================================================================
c
	integer	nx,ny,nz,i,j,k
c
	integer	ipc(*)
c
	real*8	rpc(*)
	real*8	oE(nx,ny,nz),oN(nx,ny,nz),uC(nx,ny,nz)
	real*8	cc(nx,ny,nz),oC(nx,ny,nz)
	real*8	x(nx,ny,nz),y(nx,ny,nz)
c
c ==========================================================================================
c	Perform matrix-vector multiplication
c ==========================================================================================
c
	do 300 k = 2,nz-1
		do 200 j = 2,ny-1
			do 100 i = 2,nx-1
c
				y(i,j,k) =
     1					- oN(i,j,k)	*x(i,j+1,k)
     2					- oN(i,j-1,k)	*x(i,j-1,k)
     3					- oE(i,j,k)	*x(i+1,j,k)
     4					- oE(i-1,j,k)	*x(i-1,j,k)
     5					- uC(i,j,k)	*x(i,j,k+1)
     6					- uC(i,j,k-1)	*x(i,j,k-1)
     7					+ (oC(i,j,k)+cc(i,j,k))*x(i,j,k)
c
100			continue
200		continue
300	continue
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c	Matvec_band27.f
c
c	This subroutines computes a matrix-vector multiplication, when the matrix
c	is stored with 27-diagonal (banded) format
c ==========================================================================================
c ==========================================================================================
c
	subroutine matvec_band27(nx,ny,nz,ipc,rpc,oC,cc,
     1	oE,oN,uC,oNE,oNW,uE,uW,uN,uS,uNE,uNW,uSE,uSW,x,y)
c
c ==========================================================================================
c	Declaration 
c ==========================================================================================
c
	integer	nx,ny,nz,i,j,k
c
	integer	ipc(*)
c
	real*8	tmp0,tmpU,tmpD
c
	real*8	rpc(*)
	real*8	oE(nx,ny,nz),oN(nx,ny,nz),uC(nx,ny,nz)
	real*8	oNE(nx,ny,nz),oNW(nx,ny,nz),uE(nx,ny,nz)
	real*8	uW(nx,ny,nz),uN(nx,ny,nz),uS(nx,ny,nz)
	real*8	uNE(nx,ny,nz),uNW(nx,ny,nz),uSE(nx,ny,nz)
	real*8	uSW(nx,ny,nz)
	real*8	cc(nx,ny,nz),oC(nx,ny,nz)
	real*8	x(nx,ny,nz),y(nx,ny,nz)
c
c ==========================================================================================
c	Perform matrix-vector multiplication
c ==========================================================================================
c
	do 300 k=2,nz-1
		do 200 j = 2,ny-1
			do 100 i = 2,nx-1
c
				tmp0 =
     1				- oN(i,j,k)	* x(i,j+1,k)
     2				- oN(i,j-1,k)	* x(i,j-1,k)
     3				- oE(i,j,k)	* x(i+1,j,k)
     4				- oE(i-1,j,k)	* x(i-1,j,k)
     3				- oNE(i,j,k)	* x(i+1,j+1,k)
     4				- oNE(i-1,j-1,k)* x(i-1,j-1,k)
     3				- oNW(i,j,k)	* x(i-1,j+1,k)
     4				- oNW(i+1,j-1,k)* x(i+1,j-1,k)
c
				tmpU =
     1				- uC(i,j,k)	* x(i,j,k+1)
     2				- uN(i,j,k)	* x(i,j+1,k+1)
     3				- uS(i,j,k)	* x(i,j-1,k+1)
     4				- uE(i,j,k)	* x(i+1,j,k+1)
     3				- uW(i,j,k)	* x(i-1,j,k+1)
     4				- uNE(i,j,k)	* x(i+1,j+1,k+1)
     3				- uNW(i,j,k)	* x(i-1,j+1,k+1)
     4				- uSE(i,j,k)	* x(i+1,j-1,k+1)
     4				- uSW(i,j,k)	* x(i-1,j-1,k+1)
c
				tmpD =
     1				- uC(i,j,k-1)	* x(i,j,k-1)
     2				- uS(i,j+1,k-1)	* x(i,j+1,k-1)
     3				- uN(i,j-1,k-1)	* x(i,j-1,k-1)
     4				- uW(i+1,j,k-1)	* x(i+1,j,k-1)
     3				- uE(i-1,j,k-1)	* x(i-1,j,k-1)
     4				- uSW(i+1,j+1,k-1)* x(i+1,j+1,k-1)
     3				- uSE(i-1,j+1,k-1)* x(i-1,j+1,k-1)
     4				- uNW(i+1,j-1,k-1)* x(i+1,j-1,k-1)
     4				- uNE(i-1,j-1,k-1)* x(i-1,j-1,k-1)
c
				y(i,j,k) = tmp0 + tmpU + tmpD
     1				+ (oC(i,j,k) + cc(i,j,k))*x(i,j,k)
c
100			continue
200		continue
300	continue
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c	Mresidual.f
c
c	This subroutines computes the residual for the linear equation
c ==========================================================================================
c ==========================================================================================
c
	subroutine mresidual(nx,ny,nz,ipc,rpc,ac,cc,fc,x,r)
c
c ==========================================================================================
c	Declaration 
c ==========================================================================================
c
	integer	nx,ny,nz,numdia
c
	integer	ipc(*)
c
	real*8	rpc(*),ac(nx*ny*nz,*),cc(*),fc(*),x(*),r(*)
c
c ==========================================================================================
c	Call proper routine based on number of diagonals
c ==========================================================================================
c
	numdia = ipc(11)
c
	if(numdia.eq.7) then
		call mresidual7(nx,ny,nz,ipc,rpc,ac(1,1),cc,fc,
     1			ac(1,2),ac(1,3),ac(1,4),
     2			x,r)
	elseif(numdia.eq.27) then
		call mresidual27(nx,ny,nz,ipc,rpc,ac(1,1),cc,fc,
     1		ac(1,2),ac(1,3),ac(1,4),ac(1,5),ac(1,6),ac(1,7),ac(1,8),
     2		ac(1,9),ac(1,10),ac(1,11),ac(1,12),ac(1,13),ac(1,14),
     3		x,r)
	else
		write(6,*) 'Mresidual: bad stencil type...'
		stop
	endif
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c	Mresidual7.f
c
c	This subroutines computes the residual when the matrix is stored with 
c	7-diagonal (banded) format
c ==========================================================================================
c ==========================================================================================
c
	subroutine mresidual7(nx,ny,nz,ipc,rpc,oC,cc,fc,oE,oN,uC,x,r)
c
c ==========================================================================================
c	Declaration 
c ==========================================================================================
c
	integer	nx,ny,nz,i,j,k
c
	integer	ipc(*)
c
	real*8	rpc(*)
	real*8	oE(nx,ny,nz),oN(nx,ny,nz),uC(nx,ny,nz)
	real*8	cc(nx,ny,nz),fc(nx,ny,nz),oC(nx,ny,nz)
	real*8	x(nx,ny,nz),r(nx,ny,nz)
c
c ==========================================================================================
c	Compute residual
c ==========================================================================================
c
	do 300 k = 2,nz-1
		do 200 j = 2,ny-1
			do 100 i = 2,nx-1
c
				r(i,j,k) = fc(i,j,k)
     1					+ oN(i,j,k)	*x(i,j+1,k)
     2					+ oN(i,j-1,k)	*x(i,j-1,k)
     3					+ oE(i,j,k)	*x(i+1,j,k)
     4					+ oE(i-1,j,k)	*x(i-1,j,k)
     5					+ uC(i,j,k)	*x(i,j,k+1)
     6					+ uC(i,j,k-1)	*x(i,j,k-1)
     7					- (oC(i,j,k)+cc(i,j,k))*x(i,j,k)
c
100			continue
200		continue
300	continue
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c	Mresidual27.f
c
c	This subroutines computes the residual of the linear equation when A
c	is stored with 27-diagonal (banded) format
c ==========================================================================================
c ==========================================================================================
c
	subroutine mresidual27(nx,ny,nz,ipc,rpc,oC,cc,fc,
     1	oE,oN,uC,oNE,oNW,uE,uW,uN,uS,uNE,uNW,uSE,uSW,x,r)
c
c ==========================================================================================
c	Declaration 
c ==========================================================================================
c
	integer	nx,ny,nz,i,j,k
c
	integer	ipc(*)
c
	real*8	tmp0,tmpU,tmpD
c
	real*8	rpc(*)
	real*8	oE(nx,ny,nz),oN(nx,ny,nz),uC(nx,ny,nz)
	real*8	oNE(nx,ny,nz),oNW(nx,ny,nz),uE(nx,ny,nz)
	real*8	uW(nx,ny,nz),uN(nx,ny,nz),uS(nx,ny,nz)
	real*8	uNE(nx,ny,nz),uNW(nx,ny,nz),uSE(nx,ny,nz)
	real*8	uSW(nx,ny,nz)
	real*8	cc(nx,ny,nz),oC(nx,ny,nz),fc(nx,ny,nz)
	real*8	x(nx,ny,nz),r(nx,ny,nz)
c
c ==========================================================================================
c	Compute the residual
c ==========================================================================================
c
	do 300 k=2,nz-1
		do 200 j = 2,ny-1
			do 100 i = 2,nx-1
c
				tmp0 =
     1				+ oN(i,j,k)	* x(i,j+1,k)
     2				+ oN(i,j-1,k)	* x(i,j-1,k)
     3				+ oE(i,j,k)	* x(i+1,j,k)
     4				+ oE(i-1,j,k)	* x(i-1,j,k)
     3				+ oNE(i,j,k)	* x(i+1,j+1,k)
     4				+ oNE(i-1,j-1,k)* x(i-1,j-1,k)
     3				+ oNW(i,j,k)	* x(i-1,j+1,k)
     4				+ oNW(i+1,j-1,k)* x(i+1,j-1,k)
c
				tmpU =
     1				+ uC(i,j,k)	* x(i,j,k+1)
     2				+ uN(i,j,k)	* x(i,j+1,k+1)
     3				+ uS(i,j,k)	* x(i,j-1,k+1)
     4				+ uE(i,j,k)	* x(i+1,j,k+1)
     3				+ uW(i,j,k)	* x(i-1,j,k+1)
     4				+ uNE(i,j,k)	* x(i+1,j+1,k+1)
     3				+ uNW(i,j,k)	* x(i-1,j+1,k+1)
     4				+ uSE(i,j,k)	* x(i+1,j-1,k+1)
     4				+ uSW(i,j,k)	* x(i-1,j-1,k+1)
c
				tmpD =
     1				+ uC(i,j,k-1)	* x(i,j,k-1)
     2				+ uS(i,j+1,k-1)	* x(i,j+1,k-1)
     3				+ uN(i,j-1,k-1)	* x(i,j-1,k-1)
     4				+ uW(i+1,j,k-1)	* x(i+1,j,k-1)
     3				+ uE(i-1,j,k-1)	* x(i-1,j,k-1)
     4				+ uSW(i+1,j+1,k-1)* x(i+1,j+1,k-1)
     3				+ uSE(i-1,j+1,k-1)* x(i-1,j+1,k-1)
     4				+ uNW(i+1,j-1,k-1)* x(i+1,j-1,k-1)
     4				+ uNE(i-1,j-1,k-1)* x(i-1,j-1,k-1)
c
				r(i,j,k) = fc(i,j,k)+tmp0 + tmpU + tmpD
     1				- (oC(i,j,k) + cc(i,j,k))*x(i,j,k)
c
100			continue
200		continue
300	continue
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c	Nmatvec_band.f
c
c	This subroutines computes a matrix-vector multiplication, when the matrix
c	is stored in diagonal (banded) format
c ==========================================================================================
c ==========================================================================================
c
	subroutine nmatvec_band(nx,ny,nz,ipc,rpc,ac,cc,x,y,w1)
c
c ==========================================================================================
c	Declaration 
c ==========================================================================================
c
	integer	nx,ny,nz,numdia
c
	integer	ipc(*)
c
	real*8	rpc(*),ac(nx*ny*nz,*),cc(*),x(*),y(*),w1(*)
c
c ==========================================================================================
c	Call proper routine based on number of diagonals
c ==========================================================================================
c
	numdia = ipc(11)
c
	if(numdia.eq.7) then
		call nmatvec_band7(nx,ny,nz,ipc,rpc,ac(1,1),cc,
     1			ac(1,2),ac(1,3),ac(1,4),
     2			x,y,w1)
	elseif(numdia.eq.27) then
		call nmatvec_band27(nx,ny,nz,ipc,rpc,ac(1,1),cc,
     1		ac(1,2),ac(1,3),ac(1,4),ac(1,5),ac(1,6),ac(1,7),ac(1,8),
     2		ac(1,9),ac(1,10),ac(1,11),ac(1,12),ac(1,13),ac(1,14),
     3		x,y,w1)
	else
		write(6,*) 'Matrix-vector: bad stencil type...'
		stop
	endif
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c	Nmatvec_band7.f
c
c	This subroutines computes a matrix-vector multiplication, when the matrix
c	is stored with 7-diagonal (banded) format
c ==========================================================================================
c ==========================================================================================
c
	subroutine nmatvec_band7(nx,ny,nz,ipc,rpc,oC,cc,oE,oN,uC,x,y,w1)
c
c ==========================================================================================
c	Declaration 
c ==========================================================================================
c
	integer	nx,ny,nz,i,j,k,ipkey
c
	integer	ipc(*)
c
	real*8	rpc(*)
	real*8	oE(nx,ny,nz),oN(nx,ny,nz),uC(nx,ny,nz)
	real*8	cc(nx,ny,nz),oC(nx,ny,nz)
	real*8	x(nx,ny,nz),y(nx,ny,nz),w1(nx,ny,nz)
c
c ==========================================================================================
c	Perform matrix-vector multiplication
c ==========================================================================================
c
	ipkey = ipc(10)
	call c_vec(cc,x,w1,nx,ny,nz,ipkey)
c
	do 300 k = 2,nz-1
		do 200 j = 2,ny-1
			do 100 i = 2,nx-1
c
				y(i,j,k) =
     1					- oN(i,j,k)	*x(i,j+1,k)
     2					- oN(i,j-1,k)	*x(i,j-1,k)
     3					- oE(i,j,k)	*x(i+1,j,k)
     4					- oE(i-1,j,k)	*x(i-1,j,k)
     5					- uC(i,j,k)	*x(i,j,k+1)
     6					- uC(i,j,k-1)	*x(i,j,k-1)
     7					+ oC(i,j,k)*x(i,j,k) + w1(i,j,k)
c
100			continue
200		continue
300	continue
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c	Nmatvec_band27.f
c
c	This subroutines computes a matrix-vector multiplication, when the matrix
c	is stored with 27-diagonal (banded) format
c ==========================================================================================
c ==========================================================================================
c
	subroutine nmatvec_band27(nx,ny,nz,ipc,rpc,oC,cc,
     1	oE,oN,uC,oNE,oNW,uE,uW,uN,uS,uNE,uNW,uSE,uSW,x,y,w1)
c
c ==========================================================================================
c	Declaration 
c ==========================================================================================
c
	integer	nx,ny,nz,i,j,k,ipkey
c
	integer	ipc(*)
c
	real*8	tmp0,tmpU,tmpD
c
	real*8	rpc(*)
	real*8	oE(nx,ny,nz),oN(nx,ny,nz),uC(nx,ny,nz)
	real*8	oNE(nx,ny,nz),oNW(nx,ny,nz),uE(nx,ny,nz)
	real*8	uW(nx,ny,nz),uN(nx,ny,nz),uS(nx,ny,nz)
	real*8	uNE(nx,ny,nz),uNW(nx,ny,nz),uSE(nx,ny,nz)
	real*8	uSW(nx,ny,nz)
	real*8	cc(nx,ny,nz),oC(nx,ny,nz)
	real*8	x(nx,ny,nz),y(nx,ny,nz),w1(nx,ny,nz)
c
c ==========================================================================================
c	Perform matrix-vector multiplication
c ==========================================================================================
c
	ipkey = ipc(10)
	call c_vec(cc,x,w1,nx,ny,nz,ipkey)
c
	do 300 k=2,nz-1
		do 200 j = 2,ny-1
			do 100 i = 2,nx-1
c
				tmp0 =
     1				- oN(i,j,k)	* x(i,j+1,k)
     2				- oN(i,j-1,k)	* x(i,j-1,k)
     3				- oE(i,j,k)	* x(i+1,j,k)
     4				- oE(i-1,j,k)	* x(i-1,j,k)
     3				- oNE(i,j,k)	* x(i+1,j+1,k)
     4				- oNE(i-1,j-1,k)* x(i-1,j-1,k)
     3				- oNW(i,j,k)	* x(i-1,j+1,k)
     4				- oNW(i+1,j-1,k)* x(i+1,j-1,k)
c
				tmpU =
     1				- uC(i,j,k)	* x(i,j,k+1)
     2				- uN(i,j,k)	* x(i,j+1,k+1)
     3				- uS(i,j,k)	* x(i,j-1,k+1)
     4				- uE(i,j,k)	* x(i+1,j,k+1)
     3				- uW(i,j,k)	* x(i-1,j,k+1)
     4				- uNE(i,j,k)	* x(i+1,j+1,k+1)
     3				- uNW(i,j,k)	* x(i-1,j+1,k+1)
     4				- uSE(i,j,k)	* x(i+1,j-1,k+1)
     4				- uSW(i,j,k)	* x(i-1,j-1,k+1)
c
				tmpD =
     1				- uC(i,j,k-1)	* x(i,j,k-1)
     2				- uS(i,j+1,k-1)	* x(i,j+1,k-1)
     3				- uN(i,j-1,k-1)	* x(i,j-1,k-1)
     4				- uW(i+1,j,k-1)	* x(i+1,j,k-1)
     3				- uE(i-1,j,k-1)	* x(i-1,j,k-1)
     4				- uSW(i+1,j+1,k-1)* x(i+1,j+1,k-1)
     3				- uSE(i-1,j+1,k-1)* x(i-1,j+1,k-1)
     4				- uNW(i+1,j-1,k-1)* x(i+1,j-1,k-1)
     4				- uNE(i-1,j-1,k-1)* x(i-1,j-1,k-1)
c
				y(i,j,k) = tmp0 + tmpU + tmpD
     1				+ oC(i,j,k)*x(i,j,k) + w1(i,j,k)
c
100			continue
200		continue
300	continue
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c	residual.f
c
c	This subroutines computes the residual for the non linear equation
c ==========================================================================================
c ==========================================================================================
c
	subroutine residual(nx,ny,nz,ipc,rpc,ac,cc,fc,x,r,w1)
c
c ==========================================================================================
c	Declaration 
c ==========================================================================================
c
	integer	nx,ny,nz,numdia
c
	integer	ipc(*)
c
	real*8	rpc(*),ac(nx*ny*nz,*),cc(*),fc(*),x(*),r(*),w1(*)
c
c ==========================================================================================
c	Call proper routine based on number of diagonals
c ==========================================================================================
c
	numdia = ipc(11)
c
	if(numdia.eq.7) then
		call residual7(nx,ny,nz,ipc,rpc,ac(1,1),cc,fc,
     1			ac(1,2),ac(1,3),ac(1,4),
     2			x,r,w1)
	elseif(numdia.eq.27) then
		call residual27(nx,ny,nz,ipc,rpc,ac(1,1),cc,fc,
     1		ac(1,2),ac(1,3),ac(1,4),ac(1,5),ac(1,6),ac(1,7),ac(1,8),
     2		ac(1,9),ac(1,10),ac(1,11),ac(1,12),ac(1,13),ac(1,14),
     3		x,r,w1)
	else
		write(6,*) 'Residual: bad stencil type...'
		stop
	endif
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c	residual7.f
c
c	This subroutines computes the residual when the matrix is stored with 
c	7-diagonal (banded) format
c ==========================================================================================
c ==========================================================================================
c
	subroutine residual7(nx,ny,nz,ipc,rpc,oC,cc,fc,oE,oN,uC,x,r,w1)
c
c ==========================================================================================
c	Declaration 
c ==========================================================================================
c
	integer	nx,ny,nz,i,j,k,ipkey
c
	integer	ipc(*)
c
	real*8	rpc(*)
	real*8	oE(nx,ny,nz),oN(nx,ny,nz),uC(nx,ny,nz)
	real*8	cc(nx,ny,nz),fc(nx,ny,nz),oC(nx,ny,nz)
	real*8	x(nx,ny,nz),r(nx,ny,nz),w1(nx,ny,nz)
c
c ==========================================================================================
c	Compute residual
c ==========================================================================================
c
	ipkey = ipc(10)
c
	call c_vec(cc,x,w1,nx,ny,nz,ipkey)
c
	do 300 k = 2,nz-1
		do 200 j = 2,ny-1
			do 100 i = 2,nx-1
c
				r(i,j,k) = fc(i,j,k)
     1					+ oN(i,j,k)	*x(i,j+1,k)
     2					+ oN(i,j-1,k)	*x(i,j-1,k)
     3					+ oE(i,j,k)	*x(i+1,j,k)
     4					+ oE(i-1,j,k)	*x(i-1,j,k)
     5					+ uC(i,j,k)	*x(i,j,k+1)
     6					+ uC(i,j,k-1)	*x(i,j,k-1)
     7					- oC(i,j,k)*x(i,j,k) -w1(i,j,k)
c
100			continue
200		continue
300	continue
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c	residual27.f
c
c	This subroutines computes the residual of the linear equation when A
c	is stored with 27-diagonal (banded) format
c ==========================================================================================
c ==========================================================================================
c
	subroutine residual27(nx,ny,nz,ipc,rpc,oC,cc,fc,
     1	oE,oN,uC,oNE,oNW,uE,uW,uN,uS,uNE,uNW,uSE,uSW,x,r,w1)
c
c ==========================================================================================
c	Declaration 
c ==========================================================================================
c
	integer	nx,ny,nz,i,j,k,ipkey
c
	integer	ipc(*)
c
	real*8	tmp0,tmpU,tmpD
c
	real*8	rpc(*)
	real*8	oE(nx,ny,nz),oN(nx,ny,nz),uC(nx,ny,nz)
	real*8	oNE(nx,ny,nz),oNW(nx,ny,nz),uE(nx,ny,nz)
	real*8	uW(nx,ny,nz),uN(nx,ny,nz),uS(nx,ny,nz)
	real*8	uNE(nx,ny,nz),uNW(nx,ny,nz),uSE(nx,ny,nz)
	real*8	uSW(nx,ny,nz)
	real*8	cc(nx,ny,nz),oC(nx,ny,nz),fc(nx,ny,nz)
	real*8	x(nx,ny,nz),r(nx,ny,nz),w1(nx,ny,nz)
c
c ==========================================================================================
c	Compute the residual
c ==========================================================================================
c
	ipkey = ipc(10)
c
	call c_vec(cc,x,w1,nx,ny,nz,ipkey)
c
	do 300 k=2,nz-1
		do 200 j = 2,ny-1
			do 100 i = 2,nx-1
c
				tmp0 =
     1				+ oN(i,j,k)	* x(i,j+1,k)
     2				+ oN(i,j-1,k)	* x(i,j-1,k)
     3				+ oE(i,j,k)	* x(i+1,j,k)
     4				+ oE(i-1,j,k)	* x(i-1,j,k)
     3				+ oNE(i,j,k)	* x(i+1,j+1,k)
     4				+ oNE(i-1,j-1,k)* x(i-1,j-1,k)
     3				+ oNW(i,j,k)	* x(i-1,j+1,k)
     4				+ oNW(i+1,j-1,k)* x(i+1,j-1,k)
c
				tmpU =
     1				+ uC(i,j,k)	* x(i,j,k+1)
     2				+ uN(i,j,k)	* x(i,j+1,k+1)
     3				+ uS(i,j,k)	* x(i,j-1,k+1)
     4				+ uE(i,j,k)	* x(i+1,j,k+1)
     3				+ uW(i,j,k)	* x(i-1,j,k+1)
     4				+ uNE(i,j,k)	* x(i+1,j+1,k+1)
     3				+ uNW(i,j,k)	* x(i-1,j+1,k+1)
     4				+ uSE(i,j,k)	* x(i+1,j-1,k+1)
     4				+ uSW(i,j,k)	* x(i-1,j-1,k+1)
c
				tmpD =
     1				+ uC(i,j,k-1)	* x(i,j,k-1)
     2				+ uS(i,j+1,k-1)	* x(i,j+1,k-1)
     3				+ uN(i,j-1,k-1)	* x(i,j-1,k-1)
     4				+ uW(i+1,j,k-1)	* x(i+1,j,k-1)
     3				+ uE(i-1,j,k-1)	* x(i-1,j,k-1)
     4				+ uSW(i+1,j+1,k-1)* x(i+1,j+1,k-1)
     3				+ uSE(i-1,j+1,k-1)* x(i-1,j+1,k-1)
     4				+ uNW(i+1,j-1,k-1)* x(i+1,j-1,k-1)
     4				+ uNE(i-1,j-1,k-1)* x(i-1,j-1,k-1)
c
				r(i,j,k) = fc(i,j,k)+tmp0 + tmpU + tmpD
     1				- (oC(i,j,k) + cc(i,j,k))*x(i,j,k)
c
100			continue
200		continue
300	continue
c
	return
	end
c
c ==========================================================================================
c ==========================================================================================
c	Interpolate.f
c
c	This subroutine applies the prolongation operator
c ==========================================================================================
c ==========================================================================================
c
	subroutine interpolate(nxc,nyc,nzc,nxf,nyf,nzf,xin,xout)
c
c ==========================================================================================
c	Declaration 
c ==========================================================================================
c
	integer	nxf,nyf,nzf,nxc,nyc,nzc
	integer	i,j,k,ii,jj,kk
c
	real*8	one,half,quarter,eighth
	real*8	xin(nxc,nyc,nzc),xout(nxf,nyf,nzf)
c
c ==========================================================================================
c	Setup 
c ==========================================================================================
c
c	Defines constants
c
	one = 1.d0
	half = one / 2.d0
	quarter = half/2.d0
	eighth = quarter/2.d0
c
c	Set boundary of inputs to zero
c
	call set_bound_zero(nxc,nyc,nzc,xin)
c
c ==========================================================================================
c	Apply prolongator
c ==========================================================================================
c
	do 300 k = 1,nzf-2,2
	   kk = (k-1)/2 + 1
c
	   do 200 j = 1,nyf-2,2
		jj = (j-1)/2 + 1
c
		do 100 i = 1,nxf-2,2
			ii = (i-1)/2 + 1
c
c			type 1: fine grid points common to a coarse grid point
c                       ====================================================================
c
			xout(i,j,k) = xin(ii,jj,kk)
c
c			type 2: fine grid points common to a coarse grid plane
c                       ====================================================================
c
c			fine grid points common only to y-z planes on coarse grid
c			(intermediate points between 2 grid points on x-row)
c
			xout(i+1,j,k) = half*(xin(ii,jj,kk)+
     1						xin(ii+1,jj,kk))
c
c			fine grid points common only to x-z planes on coarse grid
c			(intermediate points between 2 grid points on y-row)
c
			xout(i,j+1,k)=half*(xin(ii,jj,kk)+
     1						xin(ii,jj+1,kk))
c
c			fine grid points common only to x-y planes on coarse grid
c			(intermediate points between 2 grid points on z-row)
c
			xout(i,j,k+1) = half*(xin(ii,jj,kk)
     1				      + xin(ii,jj,kk+1))
c
c			type 3: fine grid points common to a coarse grid line
c                       ====================================================================
c
c			fine grid points common only to z planes on coarse grid
c			(intermediate points between 4 grid points on x-y plane)
c
			xout(i+1,j+1,k) = quarter*(xin(ii,jj,kk)
     1				 + xin(ii+1,jj,kk)
     1				 + xin(ii,jj+1,kk)
     1				 + xin(ii+1,jj+1,kk))
c
c			fine grid points common only to y planes on coarse grid
c			(intermediate points between 4 grid points on x-z plane)
c
			xout(i+1,j,k+1) = quarter*(xin(ii,jj,kk)
     1				 + xin(ii+1,jj,kk)
     1				 + xin(ii,jj,kk+1)
     1				 + xin(ii+1,jj,kk+1))
c
c			fine grid points common only to x planes on coarse grid
c			(intermediate points between 4 grid points on y-z plane)
c
			xout(i,j+1,k+1) = quarter*(xin(ii,jj,kk)
     1				 + xin(ii,jj+1,kk)
     1				 + xin(ii,jj,kk+1)
     1				 + xin(ii,jj+1,kk+1))
c
c			type 4: fine grid points common not common to coarse grid
c				points/line/planes
c                       ====================================================================
c
			xout(i+1,j+1,k+1) = eighth*(
     1			xin(ii,jj,kk)+ xin(ii+1,jj,kk)+ xin(ii,jj+1,kk)
     4			+ xin(ii+1,jj+1,kk)+ xin(ii,jj,kk+1)
     6			+ xin(ii+1,jj,kk+1)+ xin(ii,jj+1,kk+1)
     8			+ xin(ii+1,jj+1,kk+1))
c
100		continue
200	   continue
300	continue
c
c	Set boundary of outputs to zero
c
	call set_bound_zero(nxf,nyf,nzf,xout)
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c	Restrict.f
c
c	This subroutine applies the restriction operator
c
c	Thus is the adjoint of the prolongation operator:
c
c	R = 2**dim P^T
c
c	where dim=3 is the number of spatial dimensions, and P the prolongation
c	operator
c
c ==========================================================================================
c ==========================================================================================
c
	subroutine restrict(nxf,nyf,nzf,nxc,nyc,nzc,xin,xout)
c
c ==========================================================================================
c	Declaration 
c ==========================================================================================
c
	integer	idim
	parameter (idim=3)
c
	integer	nxf,nyf,nzf,nxc,nyc,nzc
	integer	i,j,k,ii,jj,kk
c
	real*8	xin(nxf,nyf,nzf),xout(nxc,nyc,nzc)
c
	real*8	dimfac
	real*8	tmpO,tmpU,tmpD
	real*8	one,half,quarter,eighth
c
c ==========================================================================================
c	Setup 
c ==========================================================================================
c
c	Defines constants
c
	one = 1.d0
	half = one / 2.d0
	quarter = half/2.d0
	eighth = quarter/2.d0
c
c	Make sure boundary conditions are 0 on fine grid
c
	call set_bound_zero(nxf,nyf,nzf,xin)
c
c	Define dimension factor
c
	dimfac = 2.d0**idim
c
c ==========================================================================================
c	Compute the restriction factor
c ==========================================================================================
c
	do 300 k = 2,nzc -1
	   kk = (k-1)*2 + 1
	   do 200 j = 2,nyc - 1
		jj = (j-1)*2 + 1
		do 100 i = 2,nxc-1
			ii = (i-1)*2 + 1
c
			tmpO =
     1			+ one	*xin(ii,jj,kk)
     2			+ half	*xin(ii,jj+1,kk)
     3			+ half	*xin(ii,jj-1,kk)
     4			+ half	*xin(ii+1,jj,kk)
     5			+ half	*xin(ii-1,jj,kk)
     6			+ quarter	*xin(ii+1,jj+1,kk)
     7			+ quarter	*xin(ii-1,jj+1,kk)
     8			+ quarter	*xin(ii+1,jj-1,kk)
     9			+ quarter	*xin(ii-1,jj-1,kk)
c
			tmpU =
     1			+ half	*xin(ii,jj,kk+1)
     2			+ quarter	*xin(ii,jj+1,kk+1)
     3			+ quarter	*xin(ii,jj-1,kk+1)
     4			+ quarter   	*xin(ii+1,jj,kk+1)
     5			+ quarter    	*xin(ii-1,jj,kk+1)
     6			+ eighth	*xin(ii+1,jj+1,kk+1)
     7			+ eighth     	*xin(ii-1,jj+1,kk+1)
     8			+ eighth     	*xin(ii+1,jj-1,kk+1)
     9			+ eighth     	*xin(ii-1,jj-1,kk+1)
c
			tmpD =
     1			+ half      	*xin(ii,jj,kk-1)
     2			+ quarter   	*xin(ii,jj+1,kk-1)
     3			+ quarter   	*xin(ii,jj-1,kk-1)
     4			+ quarter       *xin(ii+1,jj,kk-1)
     5			+ quarter       *xin(ii-1,jj,kk-1)
     6			+ eighth     	*xin(ii+1,jj+1,kk-1)
     7			+ eighth     	*xin(ii-1,jj+1,kk-1)
     8			+ eighth     	*xin(ii+1,jj-1,kk-1)
     9			+ eighth     	*xin(ii-1,jj-1,kk-1)
c
			xout(i,j,k) = tmpO + tmpU + tmpD
c
100		continue
200	    continue
300	continue
c
c	Make sure boundary conditions are 0 on coarse grid
c
	call set_bound_zero(nxc,nyc,nzc,xout)
c
	return
	end
