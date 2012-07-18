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
c	LinMG.f
c       (Linear multigrid method for elliptic PDEs, based on Michael Holst MG package)
c
c	Version 1:      7/4/07
c
c	Author: Patrice Koehl (in collaboration with Marc Delarue)
c
c	This file contains a series of subroutines for solving modified Poisson-Boltzmann
c	equations. It is inspired form Michael Holst MG code, which is available at:
c	http://www.fetk.org/codes/pmg/index.html
c
c ==========================================================================================
c ==========================================================================================
c 	lin_mg_vcycle.f
c ==========================================================================================
c ==========================================================================================
c
c	Linear multilevel method.
c
c	Input:
c		fine and coarse grid discrete linear operators L_h, L_H
c		fine grid rhs : f_h
c		fine grid approximate solution: u_h
c
c	Output:
c		fine grid improved solution: u_h
c
c	The two grid algorithm is:
c		(1) pre-smoothing	    u1_h = smooth(L_h,f_h,u_h)
c		(2) restrict defect	    d_H = r(L_h(u1_h) - f_h)
c		(3) Solve for correction    c_H = L_H^{-1} d_H
c		(4) prolongate and correct  u2_h = u1_h - p(c_H)
c		(5) post-smoothing	    u_h = smooth(L_h,f_h,u2_h)
c
c	Notes:
c		(1) u1_h must be kept on each level until "c_H" is computed
c		(2) u_h on each level is stored in the "x" array
c		(3) d_H is f_h on the next coarser grid
c		(4) c_h is u_h on the next coarser grid
c		(5) d_H is stored in the r array (must be kept for post-smooth)
c		(6) f_h is stored in the fc array
c		(7) L_h on all levels is stored in the ac array
c		(8) signs maybe reversed; i.e., residual is used in place of the
c		    defect, ...
c
	subroutine lin_mg_vcycle(nx,ny,nz,x,iz,w0,w1,w2,w3,
     1		istop,itmax,iters,ierror,nlev,ilev,nlev_real,mgsolv,
     2		iok,iinfo,epsilon,errtol,omega,nu1,nu2,mgsmoo,
     3		ipc,rpc,ac,cc,fc)
c
c ==========================================================================================
c 	Declare all variables
c ==========================================================================================
c
c	Input and local variables:
c ==========================================================================================
c
	integer	iok,ilev,iinfo,nlev,level,lev
	integer	itmax,iters,ierror,istop,nu1,nu2,mgsmoo
	integer	itmax_s,iters_s,nuuu,ivariable_v,mgsmoo_s,iresid
	integer	nx,ny,nz,nxf,nyf,nzf,nxc,nyc,nzc
	integer	lpv,m,n,lda,mgsolv,nlev_real,iadjoint
c
	integer	ipc(*),iz(50,*)
c
	real*8	xtest
	real*8	omega,errtol,epsilon,errtol_s
	real*8	rsden,rsnrm,orsnrm,xnorm1,xnum,xden,xdamp,dotvect
	real*8	x(*),w0(*),w1(*),w2(*),w3(*)
	real*8	rpc(*),ac(*),cc(*),fc(*)
c
	real*4	start,tend
c
c ==========================================================================================
c 	Setup
c ==========================================================================================
c
c	Recover level information
c ==========================================================================================
c
	level = 1
	lev = (ilev-1) + level
c
c	Recover gridsizes
c ==========================================================================================
c
	nxf = nx
	nyf = ny
	nzf = nz
c
	call make_coarse(nlev-1,nxf,nyf,nzf,nxc,nyc,nzc)
c
	if(iok.eq.1) then
		call print_residual(iok,-1,0.d0,0.d0,0.d0)
	endif
c
c	Compute denominator for stopping criterion
c ==========================================================================================
c
	if(iok.ne.0) then
		if(istop.eq.0) then
			rsden = 1.d0
		elseif(istop.eq.1) then
			rsden = xnorm1(nxf,nyf,nzf,fc(iz(1,lev)))
		else
			write(6,*) 'nlin_mg: bad istop value'
			stop
		endif
		if(rsden.eq.0.d0) then
			rsden = 1.0d0
			write(6,*) 'Warning...in nlin_mg...rhs is 0'
		endif
		rsnrm = rsden
		orsnrm = rsnrm
		if(iok.eq.1) call print_residual(iok,0,rsnrm,rsden,
     1			orsnrm)
	endif
c
c ==========================================================================================
c 	Solve directly if nlev = 1
c ==========================================================================================
c
	if(nlev.eq.1) then
c
		if(mgsolv.eq.0) then
c
c			Use iterative method
c			=====================
c
			iresid = 0
			iadjoint = 0
			itmax_s = 100
			iters_s = 0
			errtol_s = epsilon
			mgsmoo_s = 4
			call set_all_zeros(nxf,nyf,nzf,x(iz(1,lev)))
			call smoothing(nxf,nyf,nzf,
     1			ipc(iz(5,lev)),rpc(iz(6,lev)),
     2			ac(iz(7,lev)),cc(iz(1,lev)),fc(iz(1,lev)),
     3			x(iz(1,lev)),w1,w2,w3,
     4			itmax_s,iters_s,errtol_s,omega,iresid,iadjoint,
     5			mgsmoo_s)
c
			if(iters_s .ge.itmax_s) then
				write(6,*) 
     1				'problem on the coarse grid: iters_s =',
     2				iters_s
			endif
c
		elseif(mgsolv.eq.1) then
c
c			Use direct method
c			=====================
c
			lpv = lev + 1
c
c			Setup for band format
c
			n = ipc((iz(5,lpv)-1)+1)
			m = ipc((iz(5,lpv)-1)+2)
			lda = ipc((iz(5,lpv)-1)+3)
c
			call vectcopy_small(nxf,nyf,nzf,fc(iz(1,lev)),
     1			w1)
			call dpbsl(ac(iz(7,lpv)),lda,n,m,w1)
			call vectcopy_large(nxf,nyf,nzf,w1,x(iz(1,lev)))
			call set_bound_zero(nxf,nyf,nzf,x(iz(1,lev)))
c
		else
c
			write(6,*) 'Invalid coarse solver in lin_mg'
			stop
c
		endif
c
c		Compute stopping criteria
c
		iters = 1
		if(iok.eq.0) then
			orsnrm = rsnrm
			if(istop.eq.1) then
				call mresidual(nxf,nyf,nzf,
     1				ipc(iz(5,lev)),rpc(iz(6,lev)),
     2				ac(iz(7,lev)),cc(iz(1,lev)),
     3				fc(iz(1,lev)),x(iz(1,lev)),w1)
				rsnrm = xnorm1(nxf,nyf,nzf,w1)
			else
				write(6,*) 'In lin_mg, bad istop...'
				stop
			endif
			if(iok.eq.1) call print_residual(iok,iters,
     1					rsnrm,rsden,orsnrm)
		endif
c
c		Return
c
		goto 500
c
	endif
c
c ==========================================================================================
c 	Begin mg iteration
c ==========================================================================================
c
c	Remark: nxf,nyf,nzf change during loop
c
	iters = 0
100	continue
c
c		fine level initialization
c		==========================
c
		level = 1
		lev = (ilev-1) + level
c
c		nu1 pre-smoothing on fine grid
c		==============================
c
		iresid   = 1
		iadjoint = 0
		iters_s  = 0
		errtol_s = 0.d0
		nuuu     = ivariable_v(nu1,lev)
c
		call smoothing(nxf,nyf,nzf,
     1		ipc(iz(5,lev)),rpc(iz(6,lev)),
     2		ac(iz(7,lev)),cc(iz(1,lev)),fc(iz(1,lev)),
     3		x(iz(1,lev)),w2,w3,w1,
     4		nuuu,iters_s,errtol_s,omega,iresid,iadjoint,
     5		mgsmoo)
c
		call vectcopy(nxf,nyf,nzf,w1,w0(iz(1,lev)))
c
c ==========================================================================================
c 	Begin cycling down to coarse grid
c ==========================================================================================
c
		do 200 level = 2,nlev
c
			lev = (ilev-1)+level
c
c			Find new grid size
c			===================
c
			call make_coarse(1,nxf,nyf,nzf,nxc,nyc,nzc)
c
c			Restrict residual on coarser grid
c			=====================================
c
			call restrict(nxf,nyf,nzf,nxc,nyc,nzc,
     1				w1,w0(iz(1,level)))
c
c			New grid size
c
			nxf = nxc
			nyf = nyc
			nzf = nzc
c
c			If not yet on coarsest level ....
c			=====================================
c
			if(level.ne.nlev) then
c
c				nu1 pre-smoothing, with computation of residual
c				===============================================
c
				call set_all_zeros(nxf,nyf,nzf,
     1					x(iz(1,lev)))
c
				iresid   = 1
				iadjoint = 0
				iters_s  = 0
				errtol_s = 0.d0
				nuuu     = ivariable_v(nu1,lev)
c
				call smoothing(nxf,nyf,nzf,
     1				ipc(iz(5,lev)),rpc(iz(6,lev)),
     2				ac(iz(7,lev)),cc(iz(1,lev)),
     3				w0(iz(1,lev)),x(iz(1,lev)),w2,w3,w1,
     4				nuuu,iters_s,errtol_s,omega,iresid,
     5				iadjoint,mgsmoo)
c
			endif
c
200		continue
c
c ==========================================================================================
c 	Solve on coarse grid
c ==========================================================================================
c
		level = nlev
		lev = (ilev-1)+level
c
		if(mgsolv.eq.0) then
c
c			Use iterative method
c			=====================
c
			iresid = 0
			iadjoint = 0
			itmax_s = 100
			iters_s = 0
			errtol_s = epsilon
			mgsmoo_s = 4
c
			call set_all_zeros(nxf,nyf,nzf,x(iz(1,lev)))
			call smoothing(nxf,nyf,nzf,
     1			ipc(iz(5,lev)),rpc(iz(6,lev)),
     2			ac(iz(7,lev)),cc(iz(1,lev)),w0(iz(1,lev)),
     3			x(iz(1,lev)),w1,w2,w3,
     4			itmax_s,iters_s,errtol_s,omega,iresid,iadjoint,
     5			mgsmoo_s)
c
			if(iters_s .ge.itmax_s) then
			write(6,*) 
     1			'problem on the coarse grid: iters_s =',
     2			iters_s
			endif
c
		elseif(mgsolv.eq.1) then
c
c		Use direct method
c		=====================
c
			lpv = lev + 1
c
c			Setup for band format
c
			n = ipc((iz(5,lpv)-1)+1)
			m = ipc((iz(5,lpv)-1)+2)
			lda = ipc((iz(5,lpv)-1)+3)
c
			call vectcopy_small(nxf,nyf,nzf,fc(iz(1,lev)),
     1			w1)
			call dpbsl(ac(iz(7,lpv)),lda,n,m,w1)
			call vectcopy_large(nxf,nyf,nzf,w1,x(iz(1,lev)))
			call set_bound_zero(nxf,nyf,nzf,x(iz(1,lev)))
c
		else
c
			write(6,*) 'Invalid coarse solver in lin_mg'
			stop
c
		endif
c
c ==========================================================================================
c 	Cycle back to fine grid
c ==========================================================================================
c
		do 300 level = nlev-1,1,-1
c
			lev = (ilev-1)+level
c
c			Find new grid size
c			=====================
c
			call make_fine(1,nxf,nyf,nzf,nxc,nyc,nzc)
c
c			Interpolate to next finer grid
c			================================
c
			call interpolate(nxf,nyf,nzf,nxc,nyc,nzc,
     1			x(iz(1,lev+1)),w1)
c
c			Compute damping parameter
c			=========================
c			(equivalent to standard linear cg steplength)
c
			call matvec_band(nxf,nyf,nzf,
     1			ipc(iz(5,lev+1)),rpc(iz(6,lev+1)),
     2			ac(iz(7,lev+1)),cc(iz(1,lev+1)),
     3			x(iz(1,lev+1)),w2)
c
			xnum = dotvect(nxf,nyf,nzf,x(iz(1,lev+1)),
     1				w0(iz(1,lev+1)))
			xden = dotvect(nxf,nyf,nzf,x(iz(1,lev+1)),w2)
			xdamp = xnum / xden
c
c			New grid size
c			=============
c
			nxf = nxc
			nyf = nyc
			nzf = nzc
c
c			coarse grid correction
c			======================
c
			call vectsmult(nxf,nyf,nzf,xdamp,w1,
     1			x(iz(1,lev)))
c
c			nu2 post-smooting for correction (no residual)
c			==============================================
c
			iresid   = 0
			iadjoint = 1
			iters_s  = 0
			errtol_s = 0.d0
			nuuu     = ivariable_v(nu2,lev)
c
			if(level.eq.1) then
				call smoothing(nxf,nyf,nzf,
     1				ipc(iz(5,lev)),rpc(iz(6,lev)),
     2				ac(iz(7,lev)),cc(iz(1,lev)),
     3				fc(iz(1,lev)),x(iz(1,lev)),w1,w2,w3,
     4				nuuu,iters_s,errtol_s,omega,iresid,
     5				iadjoint,mgsmoo)
			else
				call smoothing(nxf,nyf,nzf,
     1				ipc(iz(5,lev)),rpc(iz(6,lev)),
     2				ac(iz(7,lev)),cc(iz(1,lev)),
     3				w0(iz(1,lev)),x(iz(1,lev)),w1,w2,w3,
     4				nuuu,iters_s,errtol_s,omega,iresid,
     5				iadjoint,mgsmoo)
			endif
c
300		continue
c
c ==========================================================================================
c 	Iteration complete: do some i/o
c ==========================================================================================
c
		iters = iters + 1
c
		if(iok.ne.0) then
c
			orsnrm = rsnrm
			if(istop.eq.0) then
				call mresidual(nxf,nyf,nzf,
     1				ipc(iz(5,lev)),rpc(iz(6,lev)),
     2				ac(iz(7,lev)),cc(iz(1,lev)),
     3				fc(iz(1,lev)),x(iz(1,lev)),w1)
				rsnrm = xnorm1(nxf,nyf,nzf,w1)
			elseif(istop.eq.1) then
				call mresidual(nxf,nyf,nzf,
     1				ipc(iz(5,lev)),rpc(iz(6,lev)),
     2				ac(iz(7,lev)),cc(iz(1,lev)),
     3				fc(iz(1,lev)),x(iz(1,lev)),w1)
				rsnrm = xnorm1(nxf,nyf,nzf,w1)
			else
				write(6,*) 'Bad istop in lin_mg'
				stop
			endif
c
			if(iok.eq.1) call print_residual(iok,iters,
     1					rsnrm,rsden,orsnrm)
c
			if(rsnrm/rsden .le. errtol) goto 500
c
		endif
c
		if(iters.ge.itmax) goto 400
c
	goto 100
c
400	continue
	ierror = 1
c
500	continue
c
	return
	end
c
c ==========================================================================================
c ==========================================================================================
c 	ivariable_v.f
c ==========================================================================================
c ==========================================================================================
c
c	This function defines the number of smoothings for a particular level in
c	a variable v-cycle method
c
	function ivariable_v(nu,level)
c
	integer	ivariable_v,nu,level
c
	ivariable_v = nu
c
	return
	end
