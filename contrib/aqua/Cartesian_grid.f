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
c	Cartesian_Grid.f
c       (Discretization of structured grids for PDE solver, based on Michael Holst MG package)
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
c 	Maxlevel.f
c ==========================================================================================
c ==========================================================================================
c
c 	Find maximum number of levels in the grid hierarchy for a multigrid solver, for
c	a given fine grid.
c	We assume that coarsening goes through factor of 2.
c
	function maxlevaqua(nx,ny,nz)
c
c ==========================================================================================
c 	Declare all variables
c ==========================================================================================
c
	integer	nx,ny,nz,nxtry,nytry,nztry
	integer	idone,level,maxlevaqua,ival
c
c ==========================================================================================
c 	Find max level
c ==========================================================================================
c
	idone = 0
	level = 0
	ival = 1
c
100	continue
	nxtry = (nx-1)/ival + 1
	nytry = (ny-1)/ival + 1
	nztry = (nz-1)/ival + 1
	if( ( (nxtry-1)*ival .ne. (nx-1)) .or. (nxtry .le.2)) idone = 1
	if( ( (nytry-1)*ival .ne. (ny-1)) .or. (nytry .le.2)) idone = 1
	if( ( (nztry-1)*ival .ne. (nz-1)) .or. (nztry .le.2)) idone = 1
	if(idone.ne.1) then
		level = level + 1
		ival = 2*ival
		goto 100
	endif
	maxlevaqua = level - 1
c
	return
	end
c
c ==========================================================================================
c ==========================================================================================
c 	Multigrid_size.f
c ==========================================================================================
c ==========================================================================================
c
c	This subroutine computes the sizes of the real and integer work arrays for
c	the multigrid solver
c	
c	Input:
c		mgcoar	:	coarsening technique
c				0= standard discretization
c				1= averaged coefficient + standard discretization
c				2=algebraic Galerkin coarsening
c
c		mgdisc	:	discretization technique
c				0=box method
c				1=fem method
c
c		mgsolv	:	coarse grid solver:
c				0 = conjugate gradient
c				1 = symmetric banded Linpack solver
c
c		nx,ny,nz:	grid dimensions
c
c		nlevel:		number of multigrid levels
c
c	Output:
c		nxc,nyc,nzc:	grid size on coarsest level
c		nf		number of unknowns on the finest mesh (nf = nx*ny*nz)
c		nc		number of unknowns on the coarsest mesh
c		narr		storage for one vector on all meshes
c		narrc		storage for one vector in all meshes, except finest
c		n_rpc		size of array containing real info			
c		n_ipc		size of array containing integer info			
c		n_iz		size of pointer array
c		ireal_tot	minimum size required for work array rwork
c		iint_tot	minimum size required for work array iwork
c
c
c	The work arrays rwork and iwork are chopped into smaller pieces according to:
c
c	real*8	ac(STORE)		: operators (i.e. A matrix) on all levels
c	real*8	helm(narr),fc(narr)	: Helmholtz terms and rhs, all levels
c	real*8	rpc(100*(nlevel+1))	: real info; all levels
c	integer ipc(100*(nlevel+1))	: integer info; all levels
c	integer	iz(50,nlevel+1)		: pointers to delimit rwork and iwork for each level
c
c	STORE depends on the discretization and coarsening of grid solver:
c
c	STORE = 4*nf + 4*narrc + NBAND*nc 	(mgdisc=box, mgcoar=stan/harm)
c	STORE = 4*nf +14*narrc + NBAND*nc 	(mgdisc=box, mgcoar=gal)
c	STORE =14*nf +14*narrc + NBAND*nc 	(mgdisc=fem, mgcoar=stan/harm/gal)
c
c	and NBAND =
c		0			  	(mgsolv=iterative)
c		1+(nxc-2)*(nyc-2)	  	(mgsolv = 7-pt banded Linpack)
c		1+(nxc-2)*(nyc-2)+(nxc-2)+1 	(mgsolv = 27-pt banded Limpack)
c
	subroutine mgszaqua(mgcoar,mgdisc,mgsolv,nx,ny,nz,nlevel,
     1		nxc,nyc,nzc,nf,nc,narr,narrc,n_rpc,n_iz,n_ipc,
     2		ireal_tot,iint_tot)
c
c ==========================================================================================
c 	Declare all variables
c ==========================================================================================
c
c	Input variables
c ==========================================================================================
c
	integer	mgcoar,mgdisc,mgsolv
	integer	nx,ny,nz,nlevel
c
c	parameters
c ==========================================================================================
c
	integer	num_nf,num_narr,num_narrc
	parameter	(num_narr = 2)
	parameter	(num_nf   = 2)
	parameter	(num_narrc = 27)
c
c	Output variables
c ==========================================================================================
c
	integer	nxc,nyc,nzc,nf,nc,narr,narrc
	integer	n_rpc,n_ipc,n_iz,ireal_tot,iint_tot
c
c	Miscellaneous:
c ==========================================================================================
c
	integer	i
	integer	nc_band,num_band,n_band,nxf,nyf,nzf,level
	integer	num_nf_oper,num_narrc_oper
c
c ==========================================================================================
c 	Get array sizes      
c ==========================================================================================
c
	nf	= nx * ny * nz
	narr	= nf
	nxf	= nx
	nyf	= ny
	nzf	= nz
	nxc	= nx
	nyc	= ny
	nzc	= nz
c
c	Loop over all levels, where level 1 is the finest
c ==========================================================================================
c
	do 100 i = 2,nlevel
c
c		find new grid size
c
		call make_coarse(1,nxf,nyf,nzf,nxc,nyc,nzc)
c
c		New grid size
c
		nxf = nxc
		nyf = nyc
		nzf = nzc
c
c		Add unknowns on this level to total
c
		narr = narr + nxf*nyf*nzf
c
100	continue
c
	nc = nxc * nyc *nzc
	narrc = narr - nf
c
c	Flag, based on grid type: box (i.e. structured mesh), or fem (i.e. unstructured) mesh?
c ==========================================================================================
c
	if(mgdisc.eq.0) then
		num_nf_oper = 4
	elseif(mgdisc.eq.1) then
		num_nf_oper = 14
	else
		write(6,*) 
     1		'Problem with solver: invalid discretization choice!'
	endif
c
c	Flag, based on coarsening technique: Galerkin, or standard?
c ==========================================================================================
c
	if (( (mgcoar.eq.0).or.(mgcoar.eq.1)).and.(mgdisc.eq.0)) then
		num_narrc_oper = 4
	elseif(mgcoar.eq.2) then
		num_narrc_oper = 14
	else
		write(6,*) 
     1		'Problem with solver: invalid coarsening choice!'
	endif
c
c	Storage for system solver on coarser grid
c ==========================================================================================
c
	if(mgsolv.eq.0) then
		n_band = 0
	elseif (mgsolv.eq.1) then
		if (( (mgcoar.eq.0).or.(mgcoar.eq.1)).and.
     1		(mgdisc.eq.0)) then
			num_band = 1+(nxc-2)*(nyc-2)
		else
			num_band = 1+(nxc-2)*(nyc-2)+(nxc-2)+1
		endif
		nc_band = (nxc-2)*(nyc-2)*(nzc-2)
		n_band = nc_band *num_band
	else
		write(6,*) 
     1		'Problem with solver: invalid choice for system solver!'
	endif
c
c	Size of work and int arrays:
c ==========================================================================================
c
	n_rpc = 100*(nlevel+1)
c
c	Resulting total (real storage) 
c
	ireal_tot = num_narr*narr
     1			+(num_nf + num_nf_oper)*nf
     2			+ (num_narrc_oper)*narrc
     4			+ n_band + n_rpc
c
c	Resulting total (integer storage)
c
	n_iz = 50*(nlevel+1)
	n_ipc = 100*(nlevel+1)
c
	iint_tot = n_iz + n_ipc
c
	return
	end
c
c ==========================================================================================
c ==========================================================================================
c 	Make_coarse.f
c ==========================================================================================
c ==========================================================================================
c
c	This subroutine computes the number of grid points in the coarser grid, given
c	the number of grid points in a finer grid
c
	subroutine make_coarse(numlev,nx_fine,ny_fine,nz_fine,nx_coar,
     1		ny_coar,nz_coar)
c
c ==========================================================================================
c 	Declare all variables
c ==========================================================================================
c
c	Input variables
c ==========================================================================================
c
	integer	numlev,nx_fine,ny_fine,nz_fine
c
c	Output variables
c ==========================================================================================
c
	integer	nx_coar,ny_coar,nz_coar
c
c	Local variables
c ==========================================================================================
c
	integer	i
c
c ==========================================================================================
c 	Find size of coarser grid, knowing fine grid and number of levels
c ==========================================================================================
c
	nx_coar = nx_fine
	ny_coar = ny_fine
	nz_coar = nz_fine
c
	do 100 i = 1,numlev
		nx_coar = (nx_coar-1)/2 + 1
		ny_coar = (ny_coar-1)/2 + 1
		nz_coar = (nz_coar-1)/2 + 1
100	continue
c
	return
	end
c
c ==========================================================================================
c ==========================================================================================
c 	Make_fine.f
c ==========================================================================================
c ==========================================================================================
c
c	This subroutine computes the number of grid points in a fine grid, given
c	the number of grid points in a coarser grid
c
	subroutine make_fine(numlev,nx_coar,ny_coar,nz_coar,nx_fine,
     1		ny_fine,nz_fine)
c
c ==========================================================================================
c 	Declare all variables
c ==========================================================================================
c
c	Input variables
c ==========================================================================================
c
	integer	numlev,nx_fine,ny_fine,nz_fine
c
c	Output variables
c ==========================================================================================
c
	integer	nx_coar,ny_coar,nz_coar
c
c	Local variables
c ==========================================================================================
c
	integer	i
c
c ==========================================================================================
c 	Find size of finer grid, knowing coarse grid and number of levels
c ==========================================================================================
c
	nx_fine = nx_coar
	ny_fine = ny_coar
	nz_fine = nz_coar
c
	do 100 i = 1,numlev
		nx_fine = (nx_fine-1)*2 + 1
		ny_fine = (ny_fine-1)*2 + 1
		nz_fine = (nz_fine-1)*2 + 1
100	continue
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c 	Build_structure.f
c ==========================================================================================
c ==========================================================================================
c
c	This subroutine builds the nested multigrid data structure in the array
c	iz
c
c	Note: iz(*,) indexes into the grid function arrays for each level
c	i=(1,...,nlevel+1) as follows:
c
c	fun(i)	= 	fun(iz(1,i))
c	bndx(i)	=	bndx(iz(2,i))
c	bndy(i)	=	bndx(iz(3,i))
c	bndz(i)	=	bndx(iz(4,i))
c	ipc(i)	=	ipc(iz(5,i))
c	rpc(i)	=	rpc(iz(6,i))
c	oper(i) = 	oper(iz(7,i))
c	grdx(i)	=	brdx(iz(8,i))
c	grdy(i)	=	brdy(iz(9,i))
c	grdz(i)	=	brdz(iz(10,i))
c
	subroutine build_structure(nx,ny,nz,nlevel,iz)
c
c ==========================================================================================
c 	Declare all variables
c ==========================================================================================
c
c	Input variables
c ==========================================================================================
c
	integer	nx,ny,nz,nlevel,level,n
	integer	nxold,nyold,nzold,nxnew,nynew,nznew
c
	integer	iz(50,*)
c
c ==========================================================================================
c	Setup
c ==========================================================================================
c
	nxnew	= nx
	nynew	= ny
	nznew	= nz
	n	= nxnew*nynew*nznew
c
c	Start with level 1
c ==========================================================================================
c
	level	= 1
c
c	Mark beginning of everything at level 1
c ==========================================================================================
c
	iz(1,level)	= 1
	iz(2,level)	= 1
	iz(3,level)	= 1
	iz(4,level)	= 1
	iz(5,level)	= 1
	iz(6,level)	= 1
	iz(7,level)	= 1
	iz(8,level)	= 1
	iz(9,level)	= 1
	iz(10,level)	= 1
c
c	Now mark beginning of everything at level 2
c
	iz(1,level+1)	= iz(1,level)	+ n
	iz(2,level+1)	= iz(2,level)	+ 2*nynew*nznew
	iz(3,level+1)	= iz(3,level)	+ 2*nxnew*nznew
	iz(4,level+1)	= iz(4,level)	+ 2*nxnew*nynew
	iz(5,level+1)	= iz(5,level)	+ 100
	iz(6,level+1)	= iz(6,level)	+ 100
	iz(7,level+1)	= iz(7,level)	+ 4*n
	iz(8,level+1)	= iz(8,level)	+ nxnew
	iz(9,level+1)	= iz(9,level)	+ nynew
	iz(10,level+1)	= iz(10,level)	+ nznew
c
c	Now mark the beginning of everything at level (level+1)
c
	do 100 level = 2,nlevel
c
		nxold = nxnew
		nyold = nynew
		nzold = nznew
c
		call make_coarse(1,nxold,nyold,nzold,nxnew,nynew,nznew)
c
		n = nxnew*nynew*nznew
c
		iz(1,level+1)	= iz(1,level)	+ n
		iz(2,level+1)	= iz(2,level)	+ 2*nynew*nznew
		iz(3,level+1)	= iz(3,level)	+ 2*nxnew*nznew
		iz(4,level+1)	= iz(4,level)	+ 2*nxnew*nynew
		iz(5,level+1)	= iz(5,level)	+ 100
		iz(6,level+1)	= iz(6,level)	+ 100
		iz(7,level+1)	= iz(7,level)	+ 14*n
		iz(8,level+1)	= iz(8,level)	+ nxnew
		iz(9,level+1)	= iz(9,level)	+ nynew
		iz(10,level+1)	= iz(10,level)	+ nznew
c
100	continue
c
	return
	end
c
c ==========================================================================================
c ==========================================================================================
c 	Build_operator.f
c ==========================================================================================
c ==========================================================================================
c
c	This subroutine build the operators, the boundary arrays, and modify affine vectors.
c
c	if (ido=0)	do only fine level
c	if (ido=1)	do only coarse levels
c	if (ido=2)	do all levels
c	if (ido=3)	rebuild the second operator at the coarsest level
c
c	Note: the fine level must be built before any coarse levels
c
	subroutine build_operator(nx,ny,nz,nlevel,ipkey,iinfo,ido,iz,
     1			mgprol,mgcoar,mgsolv,mgdisc,
     2			ipc,rpc,ac,cc,fc,xpos,ypos,zpos,gxcf,gycf,
     3			gzcf,a1cf,a2cf,a3cf,ccf,fcf)
c
c ==========================================================================================
c 	Declare all variables
c ==========================================================================================
c
c	Input variables
c ==========================================================================================
c
	integer	nx,ny,nz,nlevel,iinfo,ido,ipkey
	integer	mgprol,mgcoar,mgsolv,mgdisc
c
	integer	ipc(*),iz(50,*)
c
	real*8	rpc(*),ac(*),cc(*),fc(*)
	real*8	a1cf(*),a2cf(*),a3cf(*),ccf(*),fcf(*)
	real*8	gxcf(*),gycf(*),gzcf(*),xpos(*),ypos(*),zpos(*)
c
c	Local variables
c ==========================================================================================
c
	integer	nxx,nyy,nzz,nxold,nyold,nzold,numdia,key,k
	integer	level
c
c ==========================================================================================
c	Setup
c ==========================================================================================
c
	nxx	= nx
	nyy	= ny
	nzz	= nz
c
c ==========================================================================================
c	Build operator a on the finest level
c ==========================================================================================
c
	if(ido.eq.0.or.ido.eq.2) then
c
		level = 1
c
		call build_A(nxx,nyy,nzz,ipkey,mgdisc,numdia,
     1		ipc(iz(5,level)),rpc(iz(6,level)),ac(iz(7,level)),
     2		cc(iz(1,level)),fc(iz(1,level)),xpos(iz(8,level)),
     3		ypos(iz(9,level)),zpos(iz(10,level)),gxcf(iz(2,level)),
     4		gycf(iz(3,level)),gzcf(iz(4,level)),a1cf(iz(1,level)),
     5		a2cf(iz(1,level)),a3cf(iz(1,level)),ccf(iz(1,level)),
     6		fcf(iz(1,level)))
c
		iz(7,level+1) = iz(7,level) + numdia * nxx * nyy * nzz
c
	endif
c
c ==========================================================================================
c	Build the (nlevel-1) level operators
c ==========================================================================================
c
	if((ido.eq.1).or.(ido.eq.2).or.(ido.eq.3)) then
c
	    do 100 level = 2, nlevel
c
		nxold = nxx
		nyold = nyy
		nzold = nzz
c
		call make_coarse(1,nxold,nyold,nzold,nxx,nyy,nzz)
c
		if(ido.ne.3) then
c
c ==========================================================================================
c		Differential operators (based on mgcoar:
c		if mgcoar = 0, use "standard" discretization
c		if mgcoar = 1, use harmonic discretization
c		if mgcoar = 2, use Galerkin discretization
c ==========================================================================================
c
		   if(mgcoar.eq.0) then
c
			call build_copy0(nxx,nyy,nzz,nxold,nyold,nzold,
     1			xpos(iz(8,level)),ypos(iz(9,level)),
     2			zpos(iz(10,level)),gxcf(iz(2,level)),
     3			gycf(iz(3,level)),gzcf(iz(4,level)),
     4			a1cf(iz(1,level)),a2cf(iz(1,level)),
     5			a3cf(iz(1,level)),ccf(iz(1,level)),
     6			fcf(iz(1,level)),
     7			xpos(iz(8,level-1)),ypos(iz(9,level-1)),
     8			zpos(iz(10,level-1)),gxcf(iz(2,level-1)),
     9			gycf(iz(3,level-1)),gzcf(iz(4,level-1)),
     S			a1cf(iz(1,level-1)),a2cf(iz(1,level-1)),
     S			a3cf(iz(1,level-1)),ccf(iz(1,level-1)),
     S			fcf(iz(1,level-1)))
c
			call build_A(nxx,nyy,nzz,ipkey,mgdisc,numdia,
     1			ipc(iz(5,level)),rpc(iz(6,level)),
     2			ac(iz(7,level)),cc(iz(1,level)),fc(iz(1,level)),
     3			xpos(iz(8,level)),ypos(iz(9,level)),
     4			zpos(iz(10,level)),gxcf(iz(2,level)),
     5			gycf(iz(3,level)),gzcf(iz(4,level)),
     6			a1cf(iz(1,level)),a2cf(iz(1,level)),
     7			a3cf(iz(1,level)),ccf(iz(1,level)),
     8			fcf(iz(1,level)))
c
		   elseif(mgcoar.eq.1) then
c
			call build_harm0(nxx,nyy,nzz,nxold,nyold,nzold,
     1			xpos(iz(8,level)),ypos(iz(9,level)),
     2			zpos(iz(10,level)),gxcf(iz(2,level)),
     3			gycf(iz(3,level)),gzcf(iz(4,level)),
     4			a1cf(iz(1,level)),a2cf(iz(1,level)),
     5			a3cf(iz(1,level)),ccf(iz(1,level)),
     6			fcf(iz(1,level)),
     7			xpos(iz(8,level-1)),ypos(iz(9,level-1)),
     8			zpos(iz(10,level-1)),gxcf(iz(2,level-1)),
     9			gycf(iz(3,level-1)),gzcf(iz(4,level-1)),
     S			a1cf(iz(1,level-1)),a2cf(iz(1,level-1)),
     S			a3cf(iz(1,level-1)),ccf(iz(1,level-1)),
     S			fcf(iz(1,level-1)))
c
			call build_A(nxx,nyy,nzz,ipkey,mgdisc,numdia,
     1			ipc(iz(5,level)),rpc(iz(6,level)),
     2			ac(iz(7,level)),cc(iz(1,level)),fc(iz(1,level)),
     3			xpos(iz(8,level)),ypos(iz(9,level)),
     4			zpos(iz(10,level)),gxcf(iz(2,level)),
     5			gycf(iz(3,level)),gzcf(iz(4,level)),
     6			a1cf(iz(1,level)),a2cf(iz(1,level)),
     7			a3cf(iz(1,level)),ccf(iz(1,level)),
     8			fcf(iz(1,level)))
c
		   elseif(mgcoar.eq.2) then
c
			call build_Galerkin(nxold,nyold,nzold,nxx,nyy,nzz,
     1			ipkey,numdia,
     2			ipc(iz(5,level-1)),rpc(iz(6,level-1)),
     3			ac(iz(7,level-1)),cc(iz(1,level-1)),
     4			fc(iz(1,level-1)),
     5			ipc(iz(5,level)),rpc(iz(6,level)),
     6			ac(iz(7,level)),cc(iz(1,level)),
     7			fc(iz(1,level)))
c
		   else
c
			write(6,*) 'Bad mgcoar value'
			stop
c
		   endif
c
		   iz(7,level+1) = iz(7,level) + numdia*nxx*nyy*nzz
c
		endif
c
100	    continue
c
c ==========================================================================================
c		Build a sparse format coarse grid operator
c ==========================================================================================
c
	    if(mgsolv.eq.1) then
c
		level = nlevel
		call build_band(key,nxx,nyy,nzz,
     1		ipc(iz(5,level)),rpc(iz(6,level)),ac(iz(7,level)),
     1		ipc(iz(5,level+1)),rpc(iz(6,level+1)),ac(iz(7,level+1)))
c
	   endif
c
	endif
c
	return
	end
c
c ==========================================================================================
c ==========================================================================================
c 	Build_A.f
c ==========================================================================================
c ==========================================================================================
c
c	This subroutine discretizes a 3D PDE on a tensor product three dimensional 
c	mesh
c
	subroutine build_A(nx,ny,nz,ipkey,mgdisc,numdia,ipc,rpc,ac,cc,
     1		fc,xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf)
c
c ==========================================================================================
c 	Declare all variables
c ==========================================================================================
c
c	Input variables
c ==========================================================================================
c
	integer	nx,ny,nz,ipkey,mgdisc,numdia
c
	integer	ipc(*)
c
	real*8	a1cf(*),a2cf(*),a3cf(*)
	real*8	ccf(*),fcf(*),ac(nx*ny*nz,*),fc(*),cc(*)
	real*8	rpc(*)
	real*8	gxcf(*),gycf(*),gzcf(*)
	real*8	xf(*),yf(*),zf(*)
c
c ==========================================================================================
c 	Call appropriate routine (mgdisc = 1: regular grid, mgdisc =2, fem)
c ==========================================================================================
c
	if(mgdisc.eq.0) then
c
		call build_A_box(nx,ny,nz,ipkey,numdia,ipc,rpc,ac(1,1),
     1		cc,fc,ac(1,2),ac(1,3),ac(1,4),
     2		xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf)
c
	else
c
		write(6,*) 'Invalid discretization in BuilA'
		stop
c
	endif
c
	return
	end
c
c ==========================================================================================
c ==========================================================================================
c 	Build_A_box.f
c ==========================================================================================
c ==========================================================================================
c
c	This subroutine discretizes a 3D PDE on a tensor product three dimensional 
c	mesh (regular mesh)
c
	subroutine build_A_box(nx,ny,nz,ipkey,numdia,ipc,rpc,oC,cc,fc,
     1		oE,oN,uC,xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,
     2		fcf)
c
c ==========================================================================================
c 	Declare all variables
c ==========================================================================================
c
c	Input variables
c ==========================================================================================
c
	integer	nx,ny,nz,ipkey,numdia
c
	integer	ipc(*)
c
	real*8	a1cf(nx,ny,nz),a2cf(nx,ny,nz),a3cf(nx,ny,nz)
	real*8	ccf(nx,ny,nz),fcf(nx,ny,nz)
	real*8	rpc(*),oE(nx,ny,nz),oN(nx,ny,nz),uC(nx,ny,nz)
	real*8	cc(nx,ny,nz),fc(nx,ny,nz),oC(nx,ny,nz)
	real*8	gxcf(ny,nz,*),gycf(nx,nz,*),gzcf(nx,ny,*)
	real*8	xf(*),yf(*),zf(*)
c
c	Local variables
c ==========================================================================================
c
	integer	i,j,k,ike,jke,kke
	integer	nxm1,nym1,nzm1
c
	real*8	hx,hy,hz
	real*8	hxm1,hym1,hzm1,coef_fc
	real*8	bc_cond_e,bc_cond_w,bc_cond_n,bc_cond_s
	real*8	bc_cond_u,bc_cond_d
	real*8	coef_oE,coef_oN,coef_uC
	real*8	coef_oEm1,coef_oNm1,coef_uCm1
	real*8	diag
c
c ==========================================================================================
c 	Setup
c ==========================================================================================
c
c
c	Save problem key
c
	ipc(10) = ipkey
c
c	Note number of nonzeros in this stencil
c
	ipc(11) = 7
	ipc(12) = 1
	numdia  = 4
c
c	Determine n and define number of mesh points
c
	nxm1 = nx -1
	nym1 = ny -1
	nzm1 = nz -1
c
c	Diagonal scale factor
c
	diag = 1.0d0
c
c ==========================================================================================
c 	Build the operator
c ==========================================================================================
c
	do 300 k = 2,nz-1
c
	   hzm1 = zf(k) - zf(k-1)
	   hz   = zf(k+1) - zf(k)
c
	   do 200 j = 2,ny-1
c
		hym1 = yf(j) - yf(j-1)
		hy   = yf(j+1) - yf(j)
c
		do 100 i = 2,nx-1
c
			hxm1 = xf(i) - xf(i-1)
			hx   = xf(i+1) - xf(i)
c
c			Define coefficients
c                       =====================
c
			coef_oE  = diag *(hym1+hy)*(hzm1+hz)/(4.d0*hx)
			coef_oEm1= diag *(hym1+hy)*(hzm1+hz)/(4.d0*hxm1)
			coef_oN  = diag *(hxm1+hx)*(hzm1+hz)/(4.d0*hy)
			coef_oNm1= diag *(hxm1+hx)*(hzm1+hz)/(4.d0*hym1)
			coef_uC  = diag *(hxm1+hx)*(hym1+hy)/(4.d0*hz)
			coef_uCm1= diag *(hxm1+hx)*(hym1+hy)/(4.d0*hzm1)
			coef_fc  = diag *(hxm1+hx)*(hym1+hy)*(hzm1+hz)
     1					/8.d0
c
			fc(i,j,k) = coef_fc * fcf(i,j,k)
			cc(i,j,k) = coef_fc * ccf(i,j,k)
c
			oC(i,j,k) =
     1				coef_oE     * a1cf(i,j,k)
     2                          + coef_oEm1 * a1cf(i-1,j,k)
     3				+ coef_oN   * a2cf(i,j,k)
     4                          + coef_oNm1 * a2cf(i,j-1,k)
     5                          + coef_uC   * a3cf(i,j,k)
     6                          + coef_uCm1 * a3cf(i,j,k-1)
c
c			East neighbour
c                       =====================
c
			ike = min0(1,iabs(i-nxm1))
c
			oE(i,j,k) = ike *coef_oE*a1cf(i,j,k)
			bc_cond_e = (1-ike)*coef_oE*a1cf(i,j,k)*
     1					gxcf(j,k,2)
			fc(i,j,k) = fc(i,j,k) + bc_cond_e
c
c			North neighbor
c                       =====================
c
			jke = min0(1,iabs(j-nym1))
			oN(i,j,k) = jke*coef_oN*a2cf(i,j,k)
			bc_cond_n = (1-jke)*coef_oN*a2cf(i,j,k)*
     1					gycf(i,k,2)
			fc(i,j,k) = fc(i,j,k) + bc_cond_n
c
c			Up neighbour
c                       =====================
c
			kke = min0(1,iabs(k-nzm1))
			uC(i,j,k) = kke * coef_uC*a3cf(i,j,k)
			bc_cond_u = (1-kke)*coef_uC*a3cf(i,j,k)*
     1					gzcf(i,j,2)
			fc(i,j,k) = fc(i,j,k) + bc_cond_u
c
c			West neighbour 
c                       =====================
c
			ike = min0(1,iabs(i-2))
			bc_cond_w = (1-ike)*coef_oEm1*a1cf(i-1,j,k)*
     1					gxcf(j,k,1)
			fc(i,j,k) = fc(i,j,k) + bc_cond_w
c
			jke = min0(1,iabs(j-2))
			bc_cond_s = (1-jke)*coef_oNm1*a2cf(i,j-1,k)*
     1					gycf(i,k,1)
			fc(i,j,k) = fc(i,j,k) + bc_cond_s
c
			kke = min0(1,iabs(k-2))
			bc_cond_d = (1-kke)*coef_uCm1*a3cf(i,j,k-1)*
     1					gzcf(i,j,1)
			fc(i,j,k) = fc(i,j,k) + bc_cond_d
c
100		continue
c
200	   continue
c
300	continue
c
	return
	end
c
c ==========================================================================================
c ==========================================================================================
c 	Build_copy0.f
c ==========================================================================================
c ==========================================================================================
c
c	This subroutine translates information from a fine grid to a coarser grid
c	It does a simple copy.
c
	subroutine build_copy0(nxc,nyc,nzc,nxf,nyf,nzf,
     1			xc,yc,zc,gxc,gyc,gzc,a1c,a2c,a3c,ccc,fcc,
     2			xf,yf,zf,gxf,gyf,gzf,a1f,a2f,a3f,ccf,fcf)
c
c	Note:
c	Variables ending in "c" refer to the coarse grid, while variables ending in "f"
c	refer to the fine grid
c
c ==========================================================================================
c 	Declare all variables
c ==========================================================================================
c
c	Input variables
c ==========================================================================================
c
	integer	nxc,nyc,nzc,nxf,nyf,nzf
c
c	Coarse grid:
c
	real*8	xc(nxc),yc(nyc),zc(nzc)
	real*8	gxc(nyc,nzc,*),gyc(nxc,nzc,*),gzc(nxc,nyc,*)
	real*8	a1c(nxc,nyc,nzc),a2c(nxc,nyc,nzc),a3c(nxc,nyc,nzc)
	real*8	ccc(nxc,nyc,nzc),fcc(nxc,nyc,nzc)
c
c	Fine grid:
c
	real*8	xf(nxf),yf(nyf),zf(nzf)
	real*8	gxf(nyf,nzf,*),gyf(nxf,nzf,*),gzf(nxf,nyf,*)
	real*8	a1f(nxf,nyf,nzf),a2f(nxf,nyf,nzf),a3f(nxf,nyf,nzf)
	real*8	ccf(nxf,nyf,nzf),fcf(nxf,nyf,nzf)
c
c	Local variables
c ==========================================================================================
c
	integer	iadd,jadd,kadd,i,j,k,ii,jj,kk
c
c ==========================================================================================
c 	Compute coarse grid pde coefficients
c ==========================================================================================
c
c	How far to step in
c
	iadd	= (nxf-1)/(nxc-1)
	jadd	= (nyf-1)/(nyc-1)
	kadd	= (nzf-1)/(nzc-1)
c
	if(iadd.ne.2.or.jadd.ne.2.or.kadd.ne.2) then
		write(6,*) 'Problem with grid dimension in build_copy0'
		stop
	endif
c
	do 300 k = 1,nzc
c
		kk = 2*k -1
		zc(k) = zf(kk)
c
		do 200 j = 1,nyc
c
			jj = 2*j - 1
			yc(j) = yf(jj)
c
			do 100 i = 1,nxc
c
				ii = 2*i - 1
				xc(i) = xf(ii)
c
c	Copy Helmholtz coefficient
c ==========================================================================================
c
				ccc(i,j,k) = ccf(ii,jj,kk)
c
c	Copy source function (fixed charges)
c ==========================================================================================
c
				fcc(i,j,k) = fcf(ii,jj,kk)
c
c	East/west Laplacian coefficients     
c ==========================================================================================
c
				a1c(i,j,k) = a1f(ii,jj,kk)
c
c	North/South Laplacian coefficients     
c ==========================================================================================
c
				a2c(i,j,k) = a2f(ii,jj,kk)
c
c	Up/Down Laplacian coefficients     
c ==========================================================================================
c
				a3c(i,j,k) = a3f(ii,jj,kk)
c
100			continue
c
200		continue
c
300	continue
c
c ==========================================================================================
c	Now copy boundary values
c ==========================================================================================
c
c	Planes i=1 and i=nx
c ==========================================================================================
c
	do 500 k = 1,nzc
c
		kk = 2*k-1
c
		do 400 j = 1,nyc
c
			jj = 2 * j - 1
c
			gxc(j,k,1) = gxf(jj,kk,1)
			gxc(j,k,2) = gxf(jj,kk,2)
c
400		continue
c
500	continue
c
c	Planes j=1 and j=ny
c ==========================================================================================
c
	do 700 k = 1,nzc
c
		kk = 2*k-1
c
		do 600 i = 1,nxc
c
			ii = 2 * i - 1
c
			gyc(i,k,1) = gyf(ii,kk,1)
			gyc(i,k,2) = gyf(ii,kk,2)
c
600		continue
c
700	continue
c
c	Planes k=1 and k=nz
c ==========================================================================================
c
	do 900 j = 1,nyc
c
		jj = 2*j-1
c
		do 800 i = 1,nxc
c
			ii = 2 * i - 1
c
			gzc(i,j,1) = gzf(ii,jj,1)
			gzc(i,j,2) = gzf(ii,jj,2)
c
800		continue
c
900	continue
c
c ==========================================================================================
c	Return to calling program
c ==========================================================================================
c
	return
	end
c
c ==========================================================================================
c ==========================================================================================
c 	Build_harm0.f
c ==========================================================================================
c ==========================================================================================
c
c	This subroutine translates information from a fine grid to a coarser grid
c	It does a copy plus harmonic averaging
c
	subroutine build_harm0(nxc,nyc,nzc,nxf,nyf,nzf,
     1			xc,yc,zc,gxc,gyc,gzc,a1c,a2c,a3c,ccc,fcc,
     2			xf,yf,zf,gxf,gyf,gzf,a1f,a2f,a3f,ccf,fcf)
c
c	Note:
c	Variables ending in "c" refer to the coarse grid, while variables ending in "f"
c	refer to the fine grid
c
c ==========================================================================================
c 	Declare all variables
c ==========================================================================================
c
c	Input variables
c ==========================================================================================
c
	integer	nxc,nyc,nzc,nxf,nyf,nzf
c
c	Coarse grid:
c
	real*8	xc(nxc),yc(nyc),zc(nzc)
	real*8	gxc(nyc,nzc,*),gyc(nxc,nzc,*),gzc(nxc,nyc,*)
	real*8	a1c(nxc,nyc,nzc),a2c(nxc,nyc,nzc),a3c(nxc,nyc,nzc)
	real*8	ccc(nxc,nyc,nzc),fcc(nxc,nyc,nzc)
c
c	Fine grid:
c
	real*8	xf(nxf),yf(nyf),zf(nzf)
	real*8	gxf(nyf,nzf,*),gyf(nxf,nzf,*),gzf(nxf,nyf,*)
	real*8	a1f(nxf,nyf,nzf),a2f(nxf,nyf,nzf),a3f(nxf,nyf,nzf)
	real*8	ccf(nxf,nyf,nzf),fcf(nxf,nyf,nzf)
c
c	Local variables
c ==========================================================================================
c
	integer	iadd,jadd,kadd,i,j,k,ii,jj,kk
c
c	Harmonic functions
c ==========================================================================================
c
	real*8	harmo2
	real*8	a,b
c
	harmo2(a,b)	= 2.d0 * a * b / (a+b)
c
c ==========================================================================================
c 	Compute coarse grid pde coefficients
c ==========================================================================================
c
c	How far to step in
c
	iadd	= (nxf-1)/(nxc-1)
	jadd	= (nyf-1)/(nyc-1)
	kadd	= (nzf-1)/(nzc-1)
c
	if(iadd.ne.2.or.jadd.ne.2.or.kadd.ne.2) then
		write(6,*) 'Problem with grid dimension in build_harm0'
		stop
	endif
c
	do 300 k = 1,nzc
c
		kk = 2*k -1
		zc(k) = zf(kk)
c
		do 200 j = 1,nyc
c
			jj = 2*j - 1
			yc(j) = yf(jj)
c
			do 100 i = 1,nxc
c
				ii = 2*i - 1
				xc(i) = xf(ii)
c
c	Copy Helmholtz coefficient
c ==========================================================================================
c
				ccc(i,j,k) = ccf(ii,jj,kk)
c
c	Copy source function (fixed charges)
c ==========================================================================================
c
				fcc(i,j,k) = fcf(ii,jj,kk)
c
c	East/west Laplacian coefficients     
c ==========================================================================================
c
				a1c(i,j,k) = (
     1				+0.5d0 *harmo2(a1f(ii,jj,kk),
     1				a1f(min0(nxf,ii+1),jj,kk))
     2				+0.125d0*harmo2(a1f(ii,jj,max0(1,kk-1)),
     2				a1f(min0(nxf,ii+1),jj,max0(1,kk-1)))
     3				+0.125d0*harmo2(a1f(ii,jj,min0(nzf,
     3				kk+1)),
     3				a1f(min0(nxf,ii+1),jj,min0(nzf,kk+1)))
     4				+0.125d0*harmo2(a1f(ii,max0(1,jj-1),kk),
     4				a1f(min0(nxf,ii+1),max0(1,jj-1),kk))
     5				+0.125d0*harmo2(a1f(ii,min0(nyf,jj+1),
     5				kk),
     5				a1f(min0(nxf,ii+1),min0(nyf,jj+1),kk)))
c
c	North/South Laplacian coefficients     
c ==========================================================================================
c
				a2c(i,j,k) = (
     1				+0.5d0 *harmo2(a2f(ii,jj,kk),
     1				a2f(ii,min0(nyf,jj+1),kk))
     2				+0.125d0*harmo2(a2f(ii,jj,max0(1,kk-1)),
     2				a2f(ii,min0(nyf,jj+1),max0(1,kk-1)))
     3				+0.125d0*harmo2(a2f(ii,jj,min0(nzf,
     3				kk+1)),
     3				a2f(ii,min0(nyf,jj+1),min0(nzf,kk+1)))
     4				+0.125d0*harmo2(a2f(max0(1,ii-1),jj,kk),
     4				a2f(max0(1,ii-1),min0(nyf,jj+1),kk))
     5				+0.125d0*harmo2(a2f(min0(nxf,ii+1),jj,
     5				kk),
     5				a2f(min0(nxf,ii+1),min0(nyf,jj+1),kk)))
c
c	Up/Down Laplacian coefficients     
c ==========================================================================================
c
				a3c(i,j,k) = (
     1				+0.5d0 *harmo2(a3f(ii,jj,kk),
     1				a3f(ii,jj,min0(nzf,kk+1)))
     2				+0.125d0*harmo2(a3f(ii,max0(1,jj-1),kk),
     2				a3f(ii,max0(1,jj-1),min0(nzf,kk+1)))
     3				+0.125d0*harmo2(a3f(ii,min0(nyf,jj+1),
     3				kk),
     3				a3f(ii,min0(nyf,jj+1),min0(nzf,kk+1)))
     4				+0.125d0*harmo2(a3f(max0(1,ii-1),jj,kk),
     4				a3f(max0(1,ii-1),jj,min0(nzf,kk+1)))
     5				+0.125d0*harmo2(a3f(min0(nxf,ii+1),jj,
     5				kk),
     5				a3f(min0(nxf,ii+1),jj,min0(nzf,kk+1))))
c
100			continue
c
200		continue
c
300	continue
c
c ==========================================================================================
c	Now copy boundary values
c ==========================================================================================
c
c	Planes i=1 and i=nx
c ==========================================================================================
c
	do 500 k = 1,nzc
c
		kk = 2*k-1
c
		do 400 j = 1,nyc
c
			jj = 2 * j - 1
c
			gxc(j,k,1) = gxf(jj,kk,1)
			gxc(j,k,2) = gxf(jj,kk,2)
c
400		continue
c
500	continue
c
c	Planes j=1 and j=ny
c ==========================================================================================
c
	do 700 k = 1,nzc
c
		kk = 2*k-1
c
		do 600 i = 1,nxc
c
			ii = 2 * i - 1
c
			gyc(i,k,1) = gyf(ii,kk,1)
			gyc(i,k,2) = gyf(ii,kk,2)
c
600		continue
c
700	continue
c
c	Planes k=1 and k=nz
c ==========================================================================================
c
	do 900 j = 1,nyc
c
		jj = 2*j-1
c
		do 800 i = 1,nxc
c
			ii = 2 * i - 1
c
			gzc(i,j,1) = gzf(ii,jj,1)
			gzc(i,j,2) = gzf(ii,jj,2)
c
800		continue
c
900	continue
c
c ==========================================================================================
c	Return to calling program
c ==========================================================================================
c
	return
	end
c
c ==========================================================================================
c ==========================================================================================
c 	Build_band.f
c ==========================================================================================
c ==========================================================================================
c
c	This subroutine builds and factor a banded matrix given a matrix in diagonal form.
c	It uses the Linpack program DPBFA
c
	subroutine build_band(key,nx,ny,nz,ipc,rpc,ac,ipcB,rpcB,acB)
c
c ==========================================================================================
c 	Declare all variables
c ==========================================================================================
c
c	Input variables
c ==========================================================================================
c
	integer	key,nx,ny,nz
	integer	ipc(*)
c
	real*8	rpc(*),ac(nx*ny*nz,*)
c
c	Output variables
c ==========================================================================================
c
	integer	ipcB(*)
c
	real*8	rpcB(*),acB(*)
c
c	Local variables
c ==========================================================================================
c
	integer	numdia,n,m,lda,info
c
c ==========================================================================================
c	Build banded matrix
c ==========================================================================================
c
	numdia = ipc(11)	! Number of diagonals in banded matrix
c
	if(numdia.eq.7) then
c
		n = (nx-2)*(ny-2)*(nz-2)
		m = (nx-2)*(ny-2)
c
		lda = m+1
c
		call build_band_7(nx,ny,nz,ipc,rpc,ac(1,1),ac(1,2),
     1		ac(1,3),ac(1,4),ipcB,rpcB,acB,n,m,lda)
c
	elseif(numdia.eq.27) then
c
		n = (nx-2)*(ny-2)*(nz-2)
		m = (nx-2)*(ny-2) + (nx-2) + 1
		lda = m+1
c
		call build_band_27(nx,ny,nz,ipc,rpc,ac(1,1),ac(1,2),
     1		ac(1,3),ac(1,4),ac(1,5),ac(1,6),ac(1,7),ac(1,8),
     2		ac(1,9),ac(1,10),ac(1,11),ac(1,12),ac(1,13),ac(1,14),
     3		ipcB,rpcB,acB,n,m,lda)
c
	else
c
		write(6,*) 'In build_band: invalid stencil type ...'
		stop
c
	endif
c
c ==========================================================================================
c	Factor banded matrix
c ==========================================================================================
c
	key = 0
	info = 0
c
	call dpbfa(acB,lda,n,m,info)
c
	ipcB(4) = 1
c
	if(info.ne.0) then
		write(6,*) 'Build_band: Problem in dpbfa :',info
		stop
	endif
c
	return
	end
c
c ==========================================================================================
c ==========================================================================================
c 	Build_band_7.f
c ==========================================================================================
c ==========================================================================================
c
c	This subroutine builds the operator in banded form given its 7-diagonals
c
	subroutine build_band_7(nx,ny,nz,ipc,rpc,oC,oE,oN,uC,ipcB,rpcB,
     1			acB,n,m,lda)
c
c ==========================================================================================
c 	Declare all variables
c ==========================================================================================
c
c	Input variables
c ==========================================================================================
c
	integer	nx,ny,nz,n,m,lda
c
	integer	ipc(*)
c
	real*8	rpc(*),oC(nx,ny,nz),oE(nx,ny,nz),oN(nx,ny,nz)
	real*8	uC(nx,ny,nz)
c
c	Output variables
c ==========================================================================================
c
	integer	ipcB(*)
c
	real*8	acB(lda,*),rpcB(*)
c
c	Local variables
c ==========================================================================================
c
	integer	i,j,k,ii,jj,kk
c
c ==========================================================================================
c 	Build the operator
c ==========================================================================================
c
	ipcB(1) = n
	ipcB(2) = m
	ipcB(3) = lda
	ipcB(4) = 0
c
	jj = 0
	do 300 k = 2,nz-1
		do 200 j = 2,ny-1
			do 100 i = 2,nx-1
c
				jj = jj + 1
c
c	Diagonal term
c ==========================================================================================
c
				ii = jj
				kk = ii - jj + m + 1
				acB(kk,jj) = oC(i,j,k)
c
c	East neighbour
c ==========================================================================================
c
				ii = jj - 1
				kk = ii - jj + m + 1
				acB(kk,jj) = -oE(i-1,j,k)
c
c	North neighbour
c ==========================================================================================
c
				ii = jj - (nx-2)
				kk = ii - jj + m + 1
				acB(kk,jj) = -oN(i,j-1,k)
c
c	Up Neighbour
c ==========================================================================================
c
				ii = jj - (nx-2)*(ny-2)
				kk = ii - jj + m + 1
				acB(kk,jj) = -uC(i,j,k-1)
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
c 	Build_band_27.f
c ==========================================================================================
c ==========================================================================================
c
c	This subroutine builds the operator in banded form given its 27-diagonals
c
	subroutine build_band_27(nx,ny,nz,ipc,rpc,oC,oE,oN,uC,
     1		oNE,oNW,uE,uW,uN,uS,uNE,uNW,uSE,uSW,ipcB,rpcB,
     1			acB,n,m,lda)
c
c ==========================================================================================
c 	Declare all variables
c ==========================================================================================
c
c	Input variables
c ==========================================================================================
c
	integer	nx,ny,nz,n,m,lda
c
	integer	ipc(*),ipcB(*)
c
	real*8	rpc(*),oC(nx,ny,nz),oE(nx,ny,nz),oN(nx,ny,nz)
	real*8	uC(nx,ny,nz),oNE(nx,ny,nz),oNW(nx,ny,nz),uE(nx,ny,nz)
	real*8	uW(nx,ny,nz),uN(nx,ny,nz),uS(nx,ny,nz),uNE(nx,ny,nz)
	real*8	uNW(nx,ny,nz),uSE(nx,ny,nz),uSW(nx,ny,nz)
c
c	Output variables
c ==========================================================================================
c
	real*8	acB(lda,*),rpcB(*)
c
c	Local variables
c ==========================================================================================
c
	integer	i,j,k,ii,jj,kk
c
c ==========================================================================================
c 	Build the operator
c ==========================================================================================
c
	ipcB(1) = n
	ipcB(2) = m
	ipcB(3) = lda
	ipcB(4) = 0
c
	jj = 0
	do 300 k = 2,nz-1
		do 200 j = 2,ny-1
			do 100 i = 2,nx-1
c
				jj = jj + 1
c
c	Diagonal term
c ==========================================================================================
c
				ii = jj
				kk = ii - jj + m + 1
				acB(kk,jj) = oC(i,j,k)
c
c	East neighbour
c ==========================================================================================
c
				ii = jj - 1
				kk = ii - jj + m + 1
				acB(kk,jj) = -oE(i-1,j,k)
c
c	North neighbour
c ==========================================================================================
c
				ii = jj - (nx-2)
				kk = ii - jj + m + 1
				acB(kk,jj) = -oN(i,j-1,k)
c
c	North East neighbour
c ==========================================================================================
c
				ii = jj - (nx-2) + 1
				kk = ii - jj + m + 1
				acB(kk,jj) = -oNE(i,j-1,k)
c
c	North West neighbour
c ==========================================================================================
c
				ii = jj - (nx-2) - 1
				kk = ii - jj + m + 1
				acB(kk,jj) = -oNW(i,j-1,k)
c
c	Up Neighbour
c ==========================================================================================
c
				ii = jj - (nx-2)*(ny-2)
				kk = ii - jj + m + 1
				acB(kk,jj) = -uC(i,j,k-1)
c
c	Up East Neighbour
c ==========================================================================================
c
				ii = jj - (nx-2)*(ny-2)+1
				kk = ii - jj + m + 1
				acB(kk,jj) = -uE(i,j,k-1)
c
c	Up West Neighbour
c ==========================================================================================
c
				ii = jj - (nx-2)*(ny-2)-1
				kk = ii - jj + m + 1
				acB(kk,jj) = -uW(i,j,k-1)
c
c	Up North Neighbour
c ==========================================================================================
c
				ii = jj - (nx-2)*(ny-2) + (nx-2)
				kk = ii - jj + m + 1
				acB(kk,jj) = -uN(i,j,k-1)
c
c	Up South Neighbour
c ==========================================================================================
c
				ii = jj - (nx-2)*(ny-2) - (nx-2)
				kk = ii - jj + m + 1
				acB(kk,jj) = -uS(i,j,k-1)
c
c	Up North East Neighbour
c ==========================================================================================
c
				ii = jj - (nx-2)*(ny-2) + (nx-2) + 1
				kk = ii - jj + m + 1
				acB(kk,jj) = -uNE(i,j,k-1)
c
c	Up North West Neighbour
c ==========================================================================================
c
				ii = jj - (nx-2)*(ny-2) + (nx-2) - 1
				kk = ii - jj + m + 1
				acB(kk,jj) = -uNW(i,j,k-1)
c
c	Up South East Neighbour
c ==========================================================================================
c
				ii = jj - (nx-2)*(ny-2) - (nx-2) + 1
				kk = ii - jj + m + 1
				acB(kk,jj) = -uSE(i,j,k-1)
c
c	Up South West Neighbour
c ==========================================================================================
c
				ii = jj - (nx-2)*(ny-2) - (nx-2) - 1
				kk = ii - jj + m + 1
				acB(kk,jj) = -uSW(i,j,k-1)
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
c 	Set_bound_zero.f
c ==========================================================================================
c ==========================================================================================
c
c	This subroutine sets the potential on the boundary of the grid to 0
c
	subroutine set_bound_zero(nx,ny,nz,u)
c
c ==========================================================================================
c 	Declare all variables
c ==========================================================================================
c
c	Input variables:
c ==========================================================================================
c
	integer	nx,ny,nz
c
	real*8	u(nx,ny,nz)
c
c ==========================================================================================
c	Local variables
c ==========================================================================================
c
	integer	i,j,k
c
c ==========================================================================================
c 	Set values on the different boundary planes
c ==========================================================================================
c
c	The (i=1) and (i=nx) planes
c ==========================================================================================
c
	do 200 k = 1,nz
		do 100 j = 1,ny
			u(1,j,k) = 0.d0
			u(nx,j,k) = 0.d0
100		continue
200	continue
c
c	The (j=1) and (j=ny) planes
c ==========================================================================================
c
	do 400 k = 1,nz
		do 300 i = 1,nx
			u(i,1,k) = 0.d0
			u(i,ny,k) = 0.d0
300		continue
400	continue
c
c	The (k=1) and (k=ny) planes
c ==========================================================================================
c
	do 600 j = 1,ny
		do 500 i = 1,nx
			u(i,j,1) = 0.d0
			u(i,j,nz) = 0.d0
500		continue
600	continue
c
	return
	end
c
c ==========================================================================================
c ==========================================================================================
c 	Set_bound.f
c ==========================================================================================
c ==========================================================================================
c
c	This subroutine sets the potential on the boundary to the given boundary values
c
	subroutine set_bound(nx,ny,nz,u,gx,gy,gz)
c
c ==========================================================================================
c 	Declare all variables
c ==========================================================================================
c
c	Input variables:
c ==========================================================================================
c
	integer	nx,ny,nz
c
	real*8	u(nx,ny,nz)
	real*8	gx(ny,nz,2),gy(nx,nz,2),gz(nx,ny,2)
c
c ==========================================================================================
c	Local variables
c ==========================================================================================
c
	integer	i,j,k
c
c ==========================================================================================
c 	Set values on the different boundary planes
c ==========================================================================================
c
c	The (i=1) and (i=nx) planes
c ==========================================================================================
c
	do 200 k = 1,nz
		do 100 j = 1,ny
			u(1,j,k) = gx(j,k,1)
			u(nx,j,k) = gx(j,k,2)
100		continue
200	continue
c
c	The (j=1) and (j=ny) planes
c ==========================================================================================
c
	do 400 k = 1,nz
		do 300 i = 1,nx
			u(i,1,k) = gy(i,k,1)
			u(i,ny,k) = gy(i,k,2)
300		continue
400	continue
c
c	The (k=1) and (k=ny) planes
c ==========================================================================================
c
	do 600 j = 1,ny
		do 500 i = 1,nx
			u(i,j,1) = gz(i,j,1)
			u(i,j,nz) = gz(i,j,2)
500		continue
600	continue
c
	return
	end
c
c ==========================================================================================
c ==========================================================================================
c 	Update_A.f
c ==========================================================================================
c ==========================================================================================
c
c	This subroutine updates the matrix operator (finite volume method)
c
	subroutine update_A(nx,ny,nz,mgdisc,ipc,rpc,ac,xf,yf,zf,
     1			a1cf,a2cf,a3cf)
c
c ==========================================================================================
c 	Declare all variables
c ==========================================================================================
c
c	Input variables
c ==========================================================================================
c
	integer	nx,ny,nz,mgdisc
c
	integer	ipc(*)
c
	real*8	a1cf(*),a2cf(*),a3cf(*)
	real*8	ac(nx*ny*nz,*)
	real*8	rpc(*)
	real*8	xf(*),yf(*),zf(*)
c
c ==========================================================================================
c 	Call appropriate routine (mgdisc = 1: regular grid, mgdisc =2, fem)
c ==========================================================================================
c
	if(mgdisc.eq.0) then
c
		call update_A_box(nx,ny,nz,ipc,rpc,ac(1,1),
     1		ac(1,2),ac(1,3),ac(1,4),xf,yf,zf,a1cf,a2cf,a3cf)
c
	else
c
		write(6,*) 'Invalid discretization in BuilA'
		stop
c
	endif
c
	return
	end
c
c ==========================================================================================
c ==========================================================================================
c 	Update_A_box.f
c ==========================================================================================
c ==========================================================================================
c
c	This subroutine updates the discretization of a 3D PDE on a tensor product three dimensional 
c	mesh (regular mesh)
c
	subroutine update_A_box(nx,ny,nz,ipc,rpc,oC,oE,oN,uC,xf,yf,zf,
     1				a1cf,a2cf,a3cf)
c
c ==========================================================================================
c 	Declare all variables
c ==========================================================================================
c
c	Input variables
c ==========================================================================================
c
	integer	nx,ny,nz
c
	integer	ipc(*)
c
	real*8	a1cf(nx,ny,nz),a2cf(nx,ny,nz),a3cf(nx,ny,nz)
	real*8	rpc(*),oE(nx,ny,nz),oN(nx,ny,nz),uC(nx,ny,nz)
	real*8	oC(nx,ny,nz)
	real*8	xf(*),yf(*),zf(*)
c
c	Local variables
c ==========================================================================================
c
	integer	i,j,k,ike,jke,kke
	integer	nxm1,nym1,nzm1
c
	real*8	hx,hy,hz
	real*8	hxm1,hym1,hzm1
	real*8	coef_oE,coef_oN,coef_uC
	real*8	coef_oEm1,coef_oNm1,coef_uCm1
	real*8	diag
c
c ==========================================================================================
c 	Setup
c ==========================================================================================
c
c	Determine n and define number of mesh points
c
	nxm1 = nx -1
	nym1 = ny -1
	nzm1 = nz -1
c
c	Diagonal scale factor
c
	diag = 1.0d0
c
c ==========================================================================================
c 	Update the operator
c ==========================================================================================
c
	do 300 k = 2,nz-1
c
	   hzm1 = zf(k) - zf(k-1)
	   hz   = zf(k+1) - zf(k)
c
	   do 200 j = 2,ny-1
c
		hym1 = yf(j) - yf(j-1)
		hy   = yf(j+1) - yf(j)
c
		do 100 i = 2,nx-1
c
			hxm1 = xf(i) - xf(i-1)
			hx   = xf(i+1) - xf(i)
c
c			Define coefficients
c                       =====================
c
			coef_oE  = diag *(hym1+hy)*(hzm1+hz)/(4.d0*hx)
			coef_oEm1= diag *(hym1+hy)*(hzm1+hz)/(4.d0*hxm1)
			coef_oN  = diag *(hxm1+hx)*(hzm1+hz)/(4.d0*hy)
			coef_oNm1= diag *(hxm1+hx)*(hzm1+hz)/(4.d0*hym1)
			coef_uC  = diag *(hxm1+hx)*(hym1+hy)/(4.d0*hz)
			coef_uCm1= diag *(hxm1+hx)*(hym1+hy)/(4.d0*hzm1)
c
c			oC(i,j,k) = oC(i,j,k) +
			oC(i,j,k) = 
     1				coef_oE     * a1cf(i,j,k)
     2                          + coef_oEm1 * a1cf(i-1,j,k)
     3				+ coef_oN   * a2cf(i,j,k)
     4                          + coef_oNm1 * a2cf(i,j-1,k)
     5                          + coef_uC   * a3cf(i,j,k)
     6                          + coef_uCm1 * a3cf(i,j,k-1)
c
c			East neighbour
c                       =====================
c
			ike = min0(1,iabs(i-nxm1))
c
c			oE(i,j,k) = oE(i,j,k) + ike *coef_oE*a1cf(i,j,k)
			oE(i,j,k) = ike *coef_oE*a1cf(i,j,k)
c
c			North neighbor
c                       =====================
c
			jke = min0(1,iabs(j-nym1))
c			oN(i,j,k) = oN(i,j,k) + jke*coef_oN*a2cf(i,j,k)
			oN(i,j,k) = jke*coef_oN*a2cf(i,j,k)
c
c			Up neighbour
c                       =====================
c
			kke = min0(1,iabs(k-nzm1))
c			uC(i,j,k) = uC(i,j,k) + kke*coef_uC*a3cf(i,j,k)
			uC(i,j,k) = kke*coef_uC*a3cf(i,j,k)
c
100		continue
c
200	   continue
c
300	continue
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c 	Update_operator.f
c ==========================================================================================
c ==========================================================================================
c
c	This subroutine build the operators, the boundary arrays, and modify affine vectors.
c
c	Note: the fine level must be built before any coarse levels
c
	subroutine update_operator(nx,ny,nz,nlevel,iz,mgprol,mgcoar,
     1			mgsolv,mgdisc,ipc,rpc,ac,xpos,ypos,zpos,
     3			a1cf,a2cf,a3cf)
c
c ==========================================================================================
c 	Declare all variables
c ==========================================================================================
c
c	Input variables
c ==========================================================================================
c
	integer	nx,ny,nz,nlevel
	integer	mgprol,mgcoar,mgsolv,mgdisc
c
	integer	ipc(*),iz(50,*)
c
	real*8	rpc(*),ac(*)
	real*8	a1cf(*),a2cf(*),a3cf(*)
	real*8	xpos(*),ypos(*),zpos(*)
c
c	Local variables
c ==========================================================================================
c
	integer	nxx,nyy,nzz,nxold,nyold,nzold,key,k
	integer	level
c
c ==========================================================================================
c	Setup
c ==========================================================================================
c
	nxx	= nx
	nyy	= ny
	nzz	= nz
c
c ==========================================================================================
c	Update operator a on the finest level
c ==========================================================================================
c
	level = 1
c
	call update_A(nxx,nyy,nzz,mgdisc,
     1	ipc(iz(5,level)),rpc(iz(6,level)),ac(iz(7,level)),
     2	xpos(iz(8,level)),ypos(iz(9,level)),zpos(iz(10,level)),
     4	a1cf(iz(1,level)),a2cf(iz(1,level)),a3cf(iz(1,level)))
c
c ==========================================================================================
c	Update the (nlevel-1) level operators
c ==========================================================================================
c
	do 100 level = 2, nlevel
c
		nxold = nxx
		nyold = nyy
		nzold = nzz
c
		call make_coarse(1,nxold,nyold,nzold,nxx,nyy,nzz)
c
c ==========================================================================================
c		Differential operators (based on mgcoar:
c		if mgcoar = 0, use "standard" discretization
c		if mgcoar = 1, use harmonic discretization
c		if mgcoar = 2, use Galerkin discretization
c ==========================================================================================
c
		if(mgcoar.eq.0) then
c
			call update_copy0(nxx,nyy,nzz,nxold,nyold,nzold,
     2			a1cf(iz(1,level)),a2cf(iz(1,level)),
     3			a3cf(iz(1,level)),
     5			a1cf(iz(1,level-1)),
     6			a2cf(iz(1,level-1)),a3cf(iz(1,level-1)))
c
			call update_A(nxx,nyy,nzz,mgdisc,
     1			ipc(iz(5,level)),rpc(iz(6,level)),
     2			ac(iz(7,level)),
     3			xpos(iz(8,level)),ypos(iz(9,level)),
     4			zpos(iz(10,level)),
     6			a1cf(iz(1,level)),a2cf(iz(1,level)),
     7			a3cf(iz(1,level)))
c
		elseif(mgcoar.eq.1) then
c
			call update_harm0(nxx,nyy,nzz,nxold,nyold,nzold,
     1			a1cf(iz(1,level)),a2cf(iz(1,level)),
     2			a3cf(iz(1,level)),
     3			a1cf(iz(1,level-1)),a2cf(iz(1,level-1)),
     4			a3cf(iz(1,level-1)))
c
			call update_A(nxx,nyy,nzz,mgdisc,
     1			ipc(iz(5,level)),rpc(iz(6,level)),
     2			ac(iz(7,level)),
     3			xpos(iz(8,level)),ypos(iz(9,level)),
     4			zpos(iz(10,level)),
     6			a1cf(iz(1,level)),a2cf(iz(1,level)),
     7			a3cf(iz(1,level)))
c
		elseif(mgcoar.eq.2) then
c
			call update_Galerkin(nxold,nyold,nzold,
     1			nxx,nyy,nzz,
     2			ipc(iz(5,level-1)),
     3			ac(iz(7,level-1)),ac(iz(7,level)))
c
		else
c
			write(6,*) 'Bad mgcoar value'
			stop
c
		endif
c
100	continue
c
c ==========================================================================================
c		Build a sparse format coarse grid operator
c ==========================================================================================
c
	if(mgsolv.eq.1) then
c
		level = nlevel
		call build_band(key,nxx,nyy,nzz,ipc(iz(5,level)),
     1		rpc(iz(6,level)),ac(iz(7,level)),
     2		ipc(iz(5,level+1)),rpc(iz(6,level+1)),ac(iz(7,level+1)))
c
	endif
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c 	Update_copy0.f
c ==========================================================================================
c ==========================================================================================
c
c	This subroutine translates information from a fine grid to a coarser grid
c	It does a simple copy.
c
	subroutine update_copy0(nxc,nyc,nzc,nxf,nyf,nzf,
     1			a1c,a2c,a3c,a1f,a2f,a3f)
c
c	Note:
c	Variables ending in "c" refer to the coarse grid, while variables ending in "f"
c	refer to the fine grid
c
c ==========================================================================================
c 	Declare all variables
c ==========================================================================================
c
c	Input variables
c ==========================================================================================
c
	integer	nxc,nyc,nzc,nxf,nyf,nzf
c
c	Coarse grid:
c
	real*8	a1c(nxc,nyc,nzc),a2c(nxc,nyc,nzc),a3c(nxc,nyc,nzc)
c
c	Fine grid:
c
	real*8	a1f(nxf,nyf,nzf),a2f(nxf,nyf,nzf),a3f(nxf,nyf,nzf)
c
c	Local variables
c ==========================================================================================
c
	integer	iadd,jadd,kadd,i,j,k,ii,jj,kk
c
c ==========================================================================================
c 	Compute coarse grid pde coefficients
c ==========================================================================================
c
c	How far to step in
c
	iadd	= (nxf-1)/(nxc-1)
	jadd	= (nyf-1)/(nyc-1)
	kadd	= (nzf-1)/(nzc-1)
c
	if(iadd.ne.2.or.jadd.ne.2.or.kadd.ne.2) then
		write(6,*) 'Problem with grid dimension in build_copy0'
		stop
	endif
c
	do 300 k = 1,nzc
c
		kk = 2*k -1
c
		do 200 j = 1,nyc
c
			jj = 2*j - 1
c
			do 100 i = 1,nxc
c
				ii = 2*i - 1
c
c	East/west Laplacian coefficients     
c ==========================================================================================
c
				a1c(i,j,k) = a1f(ii,jj,kk)
c
c	North/South Laplacian coefficients     
c ==========================================================================================
c
				a2c(i,j,k) = a2f(ii,jj,kk)
c
c	Up/Down Laplacian coefficients     
c ==========================================================================================
c
				a3c(i,j,k) = a3f(ii,jj,kk)
c
100			continue
c
200		continue
c
300	continue
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c 	Update_harm0.f
c ==========================================================================================
c ==========================================================================================
c
c	This subroutine translates information from a fine grid to a coarser grid
c	It does a copy plus harmonic averaging
c
	subroutine update_harm0(nxc,nyc,nzc,nxf,nyf,nzf,
     1			a1c,a2c,a3c,a1f,a2f,a3f)
c
c	Note:
c	Variables ending in "c" refer to the coarse grid, while variables ending in "f"
c	refer to the fine grid
c
c ==========================================================================================
c 	Declare all variables
c ==========================================================================================
c
c	Input variables
c ==========================================================================================
c
	integer	nxc,nyc,nzc,nxf,nyf,nzf
c
c	Coarse grid:
c
	real*8	a1c(nxc,nyc,nzc),a2c(nxc,nyc,nzc),a3c(nxc,nyc,nzc)
c
c	Fine grid:
c
	real*8	a1f(nxf,nyf,nzf),a2f(nxf,nyf,nzf),a3f(nxf,nyf,nzf)
c
c	Local variables
c ==========================================================================================
c
	integer	iadd,jadd,kadd,i,j,k,ii,jj,kk
c
c	Harmonic functions
c ==========================================================================================
c
	real*8	harmo2
	real*8	a,b
c
	harmo2(a,b)	= 2.d0 * a * b / (a+b)
c
c ==========================================================================================
c 	Compute coarse grid pde coefficients
c ==========================================================================================
c
c	How far to step in
c
	iadd	= (nxf-1)/(nxc-1)
	jadd	= (nyf-1)/(nyc-1)
	kadd	= (nzf-1)/(nzc-1)
c
	if(iadd.ne.2.or.jadd.ne.2.or.kadd.ne.2) then
		write(6,*) 'Problem with grid dimension in build_harm0'
		stop
	endif
c
	do 300 k = 1,nzc
c
		kk = 2*k -1
c
		do 200 j = 1,nyc
c
			jj = 2*j - 1
c
			do 100 i = 1,nxc
c
				ii = 2*i - 1
c
c	East/west Laplacian coefficients     
c ==========================================================================================
c
				a1c(i,j,k) = (
     1				+0.5d0 *harmo2(a1f(ii,jj,kk),
     1				a1f(min0(nxf,ii+1),jj,kk))
     2				+0.125d0*harmo2(a1f(ii,jj,max0(1,kk-1)),
     2				a1f(min0(nxf,ii+1),jj,max0(1,kk-1)))
     3				+0.125d0*harmo2(a1f(ii,jj,min0(nzf,
     3				kk+1)),
     3				a1f(min0(nxf,ii+1),jj,min0(nzf,kk+1)))
     4				+0.125d0*harmo2(a1f(ii,max0(1,jj-1),kk),
     4				a1f(min0(nxf,ii+1),max0(1,jj-1),kk))
     5				+0.125d0*harmo2(a1f(ii,min0(nyf,jj+1),
     5				kk),
     5				a1f(min0(nxf,ii+1),min0(nyf,jj+1),kk)))
c
c	North/South Laplacian coefficients     
c ==========================================================================================
c
				a2c(i,j,k) = (
     1				+0.5d0 *harmo2(a2f(ii,jj,kk),
     1				a2f(ii,min0(nyf,jj+1),kk))
     2				+0.125d0*harmo2(a2f(ii,jj,max0(1,kk-1)),
     2				a2f(ii,min0(nyf,jj+1),max0(1,kk-1)))
     3				+0.125d0*harmo2(a2f(ii,jj,min0(nzf,
     3				kk+1)),
     3				a2f(ii,min0(nyf,jj+1),min0(nzf,kk+1)))
     4				+0.125d0*harmo2(a2f(max0(1,ii-1),jj,kk),
     4				a2f(max0(1,ii-1),min0(nyf,jj+1),kk))
     5				+0.125d0*harmo2(a2f(min0(nxf,ii+1),jj,
     5				kk),
     5				a2f(min0(nxf,ii+1),min0(nyf,jj+1),kk)))
c
c	Up/Down Laplacian coefficients     
c ==========================================================================================
c
				a3c(i,j,k) = (
     1				+0.5d0 *harmo2(a3f(ii,jj,kk),
     1				a3f(ii,jj,min0(nzf,kk+1)))
     2				+0.125d0*harmo2(a3f(ii,max0(1,jj-1),kk),
     2				a3f(ii,max0(1,jj-1),min0(nzf,kk+1)))
     3				+0.125d0*harmo2(a3f(ii,min0(nyf,jj+1),
     3				kk),
     3				a3f(ii,min0(nyf,jj+1),min0(nzf,kk+1)))
     4				+0.125d0*harmo2(a3f(max0(1,ii-1),jj,kk),
     4				a3f(max0(1,ii-1),jj,min0(nzf,kk+1)))
     5				+0.125d0*harmo2(a3f(min0(nxf,ii+1),jj,
     5				kk),
     5				a3f(min0(nxf,ii+1),jj,min0(nzf,kk+1))))
c
100			continue
c
200		continue
c
300	continue
c
	return
	end
