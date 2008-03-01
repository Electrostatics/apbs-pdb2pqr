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
c	BlasWrapper.f
c
c	Version 1:      7/4/07
c
c	Author: Patrice Koehl (in collaboration with Marc Delarue)
c
c	This file contains a set of BLAS-like routines for grid operation
c	as well as some BLAS routines.
c
c ==========================================================================================
c ==========================================================================================
c
c
c ==========================================================================================
c ==========================================================================================
c  	Vectsmult.f	
c
c	This subroutine replaces a vector x with x + a*y, where a is a constant
c ==========================================================================================
c ==========================================================================================
c
	subroutine vectsmult(nx,ny,nz,alpha,x,y)
c
c ==========================================================================================
c	Declaration
c ==========================================================================================
c
	integer	nx,ny,nz,i,j,k
c
	real*8	x(nx,ny,nz),y(nx,ny,nz),alpha
c
c ==========================================================================================
c	Do vector operation
c ==========================================================================================
c
	do 300 k = 2,nz-1
		do 200 j = 2,ny-1
			do 100 i = 2,nx-1
				y(i,j,k) = y(i,j,k) + alpha*x(i,j,k)
100			continue
200		continue
300	continue

	return
	end
c
c ==========================================================================================
c ==========================================================================================
c  	Vectcopy.f	
c
c	This subroutine copies a vector into another 
c ==========================================================================================
c ==========================================================================================
c
	subroutine vectcopy(nx,ny,nz,x,y)
c
c ==========================================================================================
c	Declaration
c ==========================================================================================
c
	integer	nx,ny,nz,i,j,k
c
	real*8	x(nx,ny,nz),y(nx,ny,nz)
c
c ==========================================================================================
c	copy vector
c ==========================================================================================
c
	do 300 k = 2,nz-1
		do 200 j = 2,ny-1
			do 100 i = 2,nx-1
				y(i,j,k) = x(i,j,k)
100			continue
200		continue
300	continue
c
	return
	end
c
c ==========================================================================================
c ==========================================================================================
c  	Vectcopy_small.f	
c
c	This subroutine copies a vector into another 
c ==========================================================================================
c ==========================================================================================
c
	subroutine vectcopy_small(nx,ny,nz,x,y)
c
c ==========================================================================================
c	Declaration
c ==========================================================================================
c
	integer	nx,ny,nz,i,j,k
c
	real*8	x(nx,ny,nz),y(nx-2,ny-2,nz-2)
c
c ==========================================================================================
c	copy vector
c ==========================================================================================
c
	do 300 k = 2,nz-1
		do 200 j = 2,ny-1
			do 100 i = 2,nx-1
				y(i-1,j-1,k-1) = x(i,j,k)
100			continue
200		continue
300	continue
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c  	Vectcopy_large.f	
c
c	This subroutine copies a vector into another 
c ==========================================================================================
c ==========================================================================================
c
	subroutine vectcopy_large(nx,ny,nz,x,y)
c
c ==========================================================================================
c	Declaration
c ==========================================================================================
c
	integer	nx,ny,nz,i,j,k
c
	real*8	x(nx-2,ny-2,nz-2),y(nx,ny,nz)
c
c ==========================================================================================
c	copy vector
c ==========================================================================================
c
	do 300 k = 2,nz-1
		do 200 j = 2,ny-1
			do 100 i = 2,nx-1
				y(i,j,k) = x(i-1,j-1,k-1)
100			continue
200		continue
300	continue
c
	return
	end
c
c ==========================================================================================
c ==========================================================================================
c  	xnorm1.f	
c
c	This subroutine computes the 1-norm of a vector on a grid (excluding boundaries)
c ==========================================================================================
c ==========================================================================================
c
	function xnorm1(nx,ny,nz,x)
c
c ==========================================================================================
c	Declaration
c ==========================================================================================
c
	integer	nx,ny,nz,i,j,k
c
	real*8	xnorm1
	real*8	x(nx,ny,nz)
c
c ==========================================================================================
c	Compute 1-norm
c ==========================================================================================
c
	xnorm1 = 0
	do 300 k = 2,nz-1
		do 200 j = 2,ny-1
			do 100 i = 2,nx-1
				xnorm1 = xnorm1 + dabs(x(i,j,k))
100			continue
200		continue
300	continue
c
	return
	end
c
c ==========================================================================================
c ==========================================================================================
c  	xnorm2.f	
c
c	This function computes the 2-norm of a vector on a grid (excluding boundaries)
c ==========================================================================================
c ==========================================================================================
c
	function xnorm2(nx,ny,nz,x)
c
c ==========================================================================================
c	Declaration
c ==========================================================================================
c
	integer	nx,ny,nz,i,j,k
c
	real*8	xnorm2
	real*8	x(nx,ny,nz)
c
c ==========================================================================================
c	Compute 2-norm
c ==========================================================================================
c
	xnorm2 = 0
	do 300 k = 2,nz-1
		do 200 j = 2,ny-1
			do 100 i = 2,nx-1
				xnorm2 = xnorm2 + x(i,j,k)*x(i,j,k)
100			continue
200		continue
300	continue
c
	xnorm2 = dsqrt(xnorm2)
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c  	dotvect.f	
c
c	This function computes the dot product of two grid vectors
c ==========================================================================================
c ==========================================================================================
c
	function dotvect(nx,ny,nz,x,y)
c
c ==========================================================================================
c	Declaration
c ==========================================================================================
c
	integer	nx,ny,nz,i,j,k
c
	real*8  dotvect
	real*8	x(nx,ny,nz),y(nx,ny,nz)
c
c ==========================================================================================
c	Compute 2-norm
c ==========================================================================================
c
	dotvect = 0
	do 300 k = 2,nz-1
		do 200 j = 2,ny-1
			do 100 i = 2,nx-1
				dotvect = dotvect + x(i,j,k)*y(i,j,k)
100			continue
200		continue
300	continue
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c  	Set_all_zeros.f	
c
c	This subroutine sets all values of a grid function to 0
c ==========================================================================================
c
	subroutine set_all_zeros(nx,ny,nz,x)
c
c ==========================================================================================
c	Declaration
c ==========================================================================================
c
	integer	nx,ny,nz,n,ipara,nproc
	integer	i,ii
	parameter	(nproc = 1)
c
	real*8	x(*)
c
c ==========================================================================================
c	Set all values to 0; allows for parallelisation
c ==========================================================================================
c
	n = nx*ny*nz
	ipara = n / nproc
c
	do 200 ii = 1,nproc
		do 100 i = 1+(ipara*(ii-1)),ipara*ii
			x(i) = 0.d0
100		continue
200	continue
c
	do 300 i = ipara*nproc + 1,n
		x(i) = 0.d0
300	continue
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c  	Scalevect.f	
c
c	This subroutine scales a vector with a constant
c ==========================================================================================
c ==========================================================================================
c
	subroutine scalevect(nx,ny,nz,fac,x)
c
c ==========================================================================================
c	Declaration
c ==========================================================================================
c
	integer	nx,ny,nz,i,j,k
c
	real*8	fac
	real*8	x(nx,ny,nz)
c
c ==========================================================================================
c	Scale vector
c ==========================================================================================
c
	do 300 k = 2,nz-1
		do 200 j = 2,ny-1
			do 100 i = 2,nx-1
				x(i,j,k) = fac*x(i,j,k)
100			continue
200		continue
300	continue
c
	return
	end
