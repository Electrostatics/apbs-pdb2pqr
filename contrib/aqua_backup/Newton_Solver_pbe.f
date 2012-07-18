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
c	Newton_Solver_pbe.f
c       (NEwton-based PDE solver for PB equations, based on Michael Holst MG package)
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
c 	Newton_Solver_pbe.f
c ==========================================================================================
c ==========================================================================================
c
c	Newton_Solver implements the inexact Newton Multilevel solver applied by Michael
c	Holst for solving the Poisson-Boltzmann equation, and currently implemented in
c	APBS.
c
	subroutine newdrivaqua(iparam,rparam,iwork,rwork,u,
     1				xpos,ypos,zpos,bndx,bndy,bndz,
     2				eps_X,eps_Y,eps_Z,helm,rhs)
c
c	Input:
c		iparam:		integer array: gives flags to Solver
c		rparam:		real array: constants for solver
c		iwork		integer work array
c		rwork		real work array
c
c	One dimensional:
c		xpos		positions of grid points along X axis
c		ypos		positions of grid points along Y axis
c		zpos		positions of grid points along Z axis
c
c	Two dimensional:
c		bndx		boundary conditions (X=1 and X=Nx planes)
c		bndy		boundary conditions (Y=1 and Y=Ny planes)
c		bndz		boundary conditions (Z=1 and Z=Nz planes)
c
c	Three dimensional:
c		eps_X		epsilon map, defined at positions (i+0.5,j,k)
c		eps_Y		epsilon map, defined at positions (i,j+0.5,k)
c		eps_Z		epsilon map, defined at positions (i,j,k+0.5)
c
c		helm		Helmholtz term in the equation
c		rhs		right hand side of the equation
c
c	Output:
c		u		the dimensionless electrostatic potential map
c
c	Dimensions of all arrays have been set in the calling program
c
c ==========================================================================================
c 	Declare all variables
c ==========================================================================================
c
c	Input variables:
c ==========================================================================================
c
	integer		iparam(*),iwork(*)
	real*8		rparam(*),rwork(*),u(*)
	real*8		xpos(*),ypos(*),zpos(*),bndx(*),bndy(*),bndz(*)
	real*8		eps_X(*),eps_Y(*),eps_Z(*),helm(*),rhs(*)
c
c	Local variables:
c ==========================================================================================
c
	integer	i,nx,ny,nz,niwk,nrwk,nlevel,ierror,maxlev,mxlevel
	integer	mgcoar,mgdisc,mgsolv
	integer k_iz,k_w1,k_w2
	integer k_ipc,k_rpc,k_ac,k_cc,k_fc
c
c	Variables from multigrid size
c ==========================================================================================
c
	integer	nxc,nyc,nzc,nf,nc,narr,narrc
	integer	n_rpc,n_ipc,n_iz,ireal_tot,iint_tot
c
c ==========================================================================================
c 	Decode input array iparam, and do some sanity check
c ==========================================================================================
c
	nrwk	= iparam(1)
	niwk	= iparam(2)
	nx	= iparam(3)
	ny	= iparam(4)
	nz	= iparam(5)
	nlevel	= iparam(6)
	mgcoar	= iparam(18)
	mgdisc	= iparam(19)
	mgsolv	= iparam(21)
c
	if(nlevel.le.0.or.nx.le.0.or.ny.le.0.or.nz.le.0) then
		write(6,*) 'nx,ny,nz and nlevel must be strictly positive!!'
		ierror = -1
		iparam(51) = ierror
		return
	endif
c
c	Check that the number of levels proposed is compatible with the grid size
c
	mxlevel = maxlev(nx,ny,nz)
	if(nlevel.gt.mxlevel) then
		write(6,*) 'The max. # of levels for your grid is: ',
     1			mxlevel
		ierror = -2
		iparam(51) = ierror
		return
	endif
c
c ==========================================================================================
c	Get grid size info for all grid levels: this will be used to partition the arrays
c	iwork and rwork so that they can store the whole hierarchy of grids
c ==========================================================================================
c
	call mgszaqua(mgcoar,mgdisc,mgsolv,nx,ny,nz,nlevel,
     1		nxc,nyc,nzc,nf,nc,narr,narrc,n_rpc,n_iz,n_ipc,
     2		ireal_tot,iint_tot)
c
c	Allocate space for two additional work vectors
c ==========================================================================================
c
	do 100 i = 1,ireal_tot
		rwork(i) = 0
100	continue
c
c	Check that work arrays (both real and integer) have enough space
c ==========================================================================================
c
	if(nrwk.lt.ireal_tot) then
		write(6,*) 'Real work space mut be at least: ',ireal_tot
		ierror = -3
		iparam(51) = ierror
		return
	endif
c
	if(niwk.lt.iint_tot) then
		write(6,*) 'Integer work space mut be at least: ',
     1			iint_tot
		ierror = -3
		iparam(51) = ierror
		return
	endif
c
c	Split the integer work array
c ==========================================================================================
c
	k_iz = 1
	k_ipc = k_iz + n_iz
c
c	Split the real work array
c ==========================================================================================
c
	k_rpc = 1
	k_cc  = k_rpc + n_rpc
	k_fc  = k_cc + narr
	k_w1  = k_fc + narr
	k_w2  = k_w1 + nf
	k_ac  = k_w2 + nf
c
c ==========================================================================================
c	Start the multigrid solver
c ==========================================================================================
c
	call newton_mg_pbe(iparam,rparam,nx,ny,nz,u,iwork(k_iz),
     1		rwork(k_w1),rwork(k_w2),iwork(k_ipc),rwork(k_rpc),
     2		rwork(k_ac),rwork(k_cc),rwork(k_fc),
     3		xpos,ypos,zpos,bndx,bndy,bndz,eps_X,eps_Y,eps_Z,helm,
     4		rhs)
c
c ==========================================================================================
c	The program should have found the solution...return to the
c	main program
c ==========================================================================================
c
	return
	end
c
c ==========================================================================================
c ==========================================================================================
c 	Newton_mg_pbe.f
c ==========================================================================================
c ==========================================================================================
c
c	This subroutine uses a Newton method, combined with a linear multigrid
c	solver, to solve modified Poisson-Boltzmann equations.
c
	subroutine newton_mg_pbe(iparam,rparam,nx,ny,nz,u,iz,w1,w2,
     1		ipc,rpc,ac,cc,fc,xpos,ypos,zpos,gxcf,gycf,gzcf,
     2		a1cf,a2cf,a3cf,ccf,fcf)
c
c ==========================================================================================
c 	Declare all variables
c ==========================================================================================
c
c	Input variables:
c ==========================================================================================
c
	integer		iparam(*),ipc(*),iz(*)
	real*8		rparam(*),rpc(*),ac(*),cc(*),fc(*)
	real*8		u(*),w1(*),w2(*)
	real*8		xpos(*),ypos(*),zpos(*),gxcf(*),gycf(*),gzcf(*)
	real*8		a1cf(*),a2cf(*),a3cf(*),ccf(*),fcf(*)
c
c	Local variables:
c ==========================================================================================
c
	integer		mgkey,nlevel,itmax,iok,iinfo,istop,ipkey,nu1,nu2
	integer		nx,ny,nz,ilev,ido,iters,ierror,nlevel_real
	integer		mgprol,mgcoar,mgsolv,mgdisc,mgsmoo
c
c	real*8		epsmach
	real*8		epsilon,errtol,omegal,omegan
c
c ==========================================================================================
c 	Decode the iparam and rparam arrays
c ==========================================================================================
c
	nlevel = iparam(6)	! Number of mesh level
	nu1    = iparam(7)	! number of pre-smoothings
	nu2    = iparam(8)	! number of post-smoothings
	mgkey  = iparam(9)	! multigrid method: 0=variable, v-cycle, 1=nest operation
	itmax  = iparam(10)	! maximum number of iterations
	istop  = iparam(11)	! stopping criteria
	iinfo  = iparam(12)	! runtime status messages :
				! 0=none, 1=some, 2=many, 3=more
	ipkey  = iparam(14)	! Defines type of equation considered
				! 0: Laplace, 1: linear PB, 2: non linear PB
				! 3: Size modified PB, 4: Generalized PB Langevin
	mgprol = iparam(17)	! type of prolongation (0=trilinear)
	mgcoar = iparam(18)	! coarse equations (0=stan, 1=harm, 2=Galerkin)
	mgdisc = iparam(19)	! discretization (0=fv, 1=fe)
	mgsmoo = iparam(20)	! type of smoothing (0=wjac, 1=gsrb, 2=sor, 3=rich, 4=cghs)
	mgsolv = iparam(21)	! coarse equation solver (0=cghs, 1=nanded Linpack)
c
	errtol = rparam(1)	! tolerance factor
	omegal = rparam(9)	! linear relaxation parameter
	omegan = rparam(10)	! non linear relaxation parameter
c
c ==========================================================================================
c 	Build the multigrid data structure in array iz
c ==========================================================================================
c
	call build_structure(nx,ny,nz,nlevel,iz)
c
c ==========================================================================================
c 	Build operator and right-hand-size on fine grid
c ==========================================================================================
c
	ido = 0
	call build_operator(nx,ny,nz,nlevel,ipkey,iinfo,ido,iz,mgprol,
     1	mgcoar,mgsolv,mgdisc,ipc,rpc,ac,cc,fc,xpos,ypos,zpos,
     2	gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf)
c
c ==========================================================================================
c 	Build operator and right-hand-size on all coarse grids
c ==========================================================================================
c
	ido = 1
	call build_operator(nx,ny,nz,nlevel,ipkey,iinfo,ido,iz,mgprol,
     1	mgcoar,mgsolv,mgdisc,ipc,rpc,ac,cc,fc,xpos,ypos,zpos,
     2	gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf)
c
c ==========================================================================================
c 	Prepare for Newton solver and solve
c ==========================================================================================
c
c	Determine machine epsilon (hard coded)
c ==========================================================================================
c
c	epsilon = epsmach(0)
	epsilon = 1.d-9
c
c	Impose zero Dirichlet boundary conditions
c ==========================================================================================
c
	call set_bound_zero(nx,ny,nz,u)
c
c	Call multigrid solver
c ==========================================================================================
c
	nlevel_real = nlevel
	iok = 1
	ilev = 1
c
	if(mgkey.eq.0) then
		call newton_vcycle_pbe(nx,ny,nz,u,iz,ccf,fcf,w1,w2,istop,
     1			itmax,iters,ierror,nlevel,ilev,nlevel_real,
     2			mgsolv,iok,iinfo,epsilon,errtol,omegan,nu1,nu2,
     3			mgsmoo,a1cf,a2cf,a3cf,ipc,rpc,ac,cc,fc)
	elseif(mgkey.eq.1) then
		call newton_nested_pbe(nx,ny,nz,u,iz,ccf,fcf,w1,w2,istop,
     1			itmax,iters,ierror,nlevel,ilev,nlevel_real,
     2			mgsolv,iok,iinfo,epsilon,errtol,omegan,nu1,nu2,
     3			mgsmoo,a1cf,a2cf,a3cf,ipc,rpc,ac,cc,fc)
	else
		write(6,*) 'Bad mgkey!'
	endif
c
c	Restore boundary conditions
c ==========================================================================================
c
	call set_bound(nx,ny,nz,u,gxcf,gycf,gzcf)
c
c ==========================================================================================
c 	Return and end
c ==========================================================================================
c
	return
	end
c
c ==========================================================================================
c ==========================================================================================
c 	newton_nested_pbe.f
c ==========================================================================================
c ==========================================================================================
c
c	Newton_nested uses nested iteration for an inexact Newton multilevel method
c
	subroutine newton_nested_pbe(nx,ny,nz,u,iz,w0,w1,w2,w3,istop,itmax,
     1			iters,ierror,nlevel,ilev,nlevel_real,mgsolv,
     2			iok,iinfo,epsilon,errtol,omega,nu1,nu2,mgsmoo,
     3			cprime,rhs,utemp,ipc,rpc,ac,cc,fc)
c
c ==========================================================================================
c 	Declare all variables
c ==========================================================================================
c
c	Input variables:
c ==========================================================================================
c
	integer	nx,ny,nz
	integer	istop,itmax,iters,ierror,level,nlevel,ilev,nlevel_real
	integer	mgsolv,mgsmoo,iok,iinfo,nu1,nu2
c
	integer	ipc(*),iz(50,*)
c
	real*8	epsilon,errtol,omega
	real*8	u(*),w0(*),w1(*),w2(*),w3(*)
	real*8	rpc(*),ac(*),cc(*),fc(*)
	real*8	cprime(*),rhs(*),utemp(*)
c
c	Local variables:
c ==========================================================================================
c
	integer	nxf,nyf,nzf,nxc,nyc,nzc
	integer	itmxd,nleveld,iterd,iokd,istpd
c
	real*8	errd
c
c	Recover grid size
c
	nxf = nx
	nyf = ny
	nzf = nz
c
	call make_coarse(nlevel-1,nxf,nyf,nzf,nxc,nyc,nzc)
c
c	Move up grids: interpolate solution, and use Newton_vcycle
c
	do 100 level=nlevel_real,ilev+1,-1
c
c		errd = errtol
		itmxd = nlevel_real
		nleveld = nlevel_real - level + 1
		iterd = 0
		iokd = iok
		istpd = istop
		call newton_vcycle_pbe(nxc,nyc,nzc,u,iz,w0,w1,w2,w3,istpd,
     1		itmxd,iterd,ierror,nleveld,level,nlevel_real,mgsolv,
     2		iokd,iinfo,epsilon,errd,omega,nu1,nu2,mgsmoo,
     3		cprime,rhs,utemp,ipc,rpc,ac,cc,fc)
c
		call make_fine(1,nxc,nyc,nzc,nxf,nyf,nzf)
c
		call interpolate(nxc,nyc,nzc,nxf,nyf,nzf,u(iz(1,level)),
     1		u(iz(1,level-1)))
c
		nxc = nxf
		nyc = nyf
		nzc = nzf
c
100	continue
c
	level = ilev
	call newton_vcycle_pbe(nxc,nyc,nzc,u,iz,w0,w1,w2,w3,istop,
     1	itmax,iters,ierror,nlevel,level,nlevel_real,mgsolv,
     2	iok,iinfo,epsilon,errtol,omega,nu1,nu2,mgsmoo,
     3		cprime,rhs,utemp,ipc,rpc,ac,cc,fc)
c
	return
	end
c
c ==========================================================================================
c ==========================================================================================
c 	newton_vcycle_pbe.f
c ==========================================================================================
c ==========================================================================================
c
c	Newton_vcycle is the inexact Newton multilevel solver
c
	subroutine newton_vcycle_pbe(nx,ny,nz,x,iz,w0,w1,w2,w3,istop,itmax,
     1			iters,ierror,nlev,ilev,nlevel_real,mgsolv,
     2			iok,iinfo,epsilon,errtol,omega,nu1,nu2,mgsmoo,
     3			cprime,rhs,xtemp,ipc,rpc,ac,cc,fc)
c
c ==========================================================================================
c 	Declare all variables
c ==========================================================================================
c
c	Input variables:
c ==========================================================================================
c
	integer	nx,ny,nz
	integer	istop,itmax,iters,ierror,nlev,ilev,nlevel_real
	integer	mgsolv,iok,iinfo,nu1,nu2,mgsmoo
c
	integer	ipc(*),iz(50,*)
c
	real*8	epsilon,errtol,omega
	real*8	x(*),w0(*),w1(*),w2(*),w3(*)
	real*8	rpc(*),ac(*),cc(*),fc(*)
	real*8	cprime(*),rhs(*),xtemp(*)
c
c	Local variables
c ==========================================================================================
c
	integer	level,lev,i
	integer	itmax_s, iters_s, ierror_s, iok_s, iinfo_s,istop_s
	integer	iter_d,itmax_d,mode,idamp,ipkey
c
	real*8	xnorm1
	real*8	bigc,ord
	real*8	rsden,rsnrm,orsnrm
	real*8	xnorm_old,xnorm_new,damp,xnorm_med
	real*8	errtol_s
c
	real*4	start,tend
c
c ==========================================================================================
c	Initialization
c ==========================================================================================
c
c	Recover level information
c ==========================================================================================
c
	level = 1
	lev  = (ilev-1) + level
c
c	Start clock and printout
c ==========================================================================================
c	
	if(iok.ne.0) then
		call print_ini
		call print_residual(iok,-1,0.0d0,0.0d0,0.0d0)
	endif
c
c	Compute denominator for stopping criteria
c ==========================================================================================
c	
c	We only implement the relative residual
c
	if(istop.eq.1) then
c
c		Compute initial residual with zero initial guess
c		This is analogous to simply taking the norm of the rhs of the
c		equation
c
		call set_all_zeros(nx,ny,nz,w1)
		call residual(nx,ny,nz,
     1			ipc(iz(5,lev)),rpc(iz(6,lev)),
     2			ac(iz(7,lev)),cc(iz(1,lev)),fc(iz(1,lev)),
     3			w1,w2,w3)
		rsden = xnorm1(nx,ny,nz,w2)
	else
		write(6,*) 'Bad istop value ....'
		stop
	endif
c
	if(rsden.eq.0.0d0) then
		rsden = 1.d0
		write(6,*) 'Warning: rhs is zero ....'
	endif
c
	rsnrm = rsden
	orsnrm = rsnrm
c
	if(iok.ne.0) then
		call print_residual(iok,0,rsnrm,rsden,orsnrm)
	endif
c
c ==========================================================================================
c	Newton iteration
c ==========================================================================================
c
c	Begin: compute residual with the initial guess
c ==========================================================================================
c
	call set_all_zeros(nx,ny,nz,x(iz(1,lev)))
	call residual(nx,ny,nz,
     1		ipc(iz(5,lev)),rpc(iz(6,lev)),
     2		ac(iz(7,lev)),cc(iz(1,lev)),fc(iz(1,lev)),
     3		x(iz(1,lev)),w0,w2)
	xnorm_old = xnorm1(nx,ny,nz,w0)
c
c	Loop over all levels
c ==========================================================================================
c
	idamp = 1
	iters = 0
c
100	continue
c
		iters = iters + 1
c
c		Compute current Jacobian of the system:
c
c		Since we consider the "standard" PBE, the Jacobian is simply
c		the sum of the operator AC and the derivative (diagonal) of
c		the non linear term
c
		ipkey = ipc(10)
c
		call get_jacobian(nx,ny,nz,nlevel_real,iz,ilev,ipkey,
     1		x,w0,cprime,rhs,cc)
c
c
c		Define tolerance, to get global-superlinear convergence
c               ========================================================
c
		bigc = 1.0d0
		ord  = 2.0d0
		errtol_s = dmin1((0.9d0*xnorm_old),
     1				 (bigc*xnorm_old**ord))
c
c		Solve the Newton equations using linear multigrid
c               ========================================================
c
		call set_all_zeros(nx,ny,nz,xtemp(iz(1,lev)))
c
c		Set parameters for multigrid
c		
		itmax_s  = 100
		istop_s  = 0
		iters_s  = 0
		ierror_s = 0
c
		iok_s    = 2
		iinfo_s  = iinfo
c
		call lin_mg_vcycle(nx,ny,nz,xtemp,iz,w0,w1,w2,w3,
     1		istop_s,itmax_s,iters_s,ierror_s,
     2		nlev,ilev,nlevel_real,mgsolv,
     3		iok_s,iinfo_s,epsilon,errtol_s,omega,nu1,nu2,mgsmoo,
     4		ipc,rpc,ac,cprime,rhs)
c
c		(note that now rhs and cprime are available as storage arrays)
c
c		Apply damping, if still enabled
c               ========================================================
c
		if(idamp.eq.1) then
c
c		    try the correction
c		    ==================
c
		    call vectcopy(nx,ny,nz,x(iz(1,lev)),w1)
c
		    damp = 1.d0
c
		    call vectsmult(nx,ny,nz,damp,xtemp(iz(1,lev)),w1)
c
		    call residual(nx,ny,nz,
     1			ipc(iz(5,lev)),rpc(iz(6,lev)),
     2			ac(iz(7,lev)),cc(iz(1,lev)),fc(iz(1,lev)),
     3			w1,w0,rhs(iz(1,lev)))
c
		    xnorm_new = xnorm1(nx,ny,nz,w0)
c
		    damp    = 1.0d0
		    iter_d  = 0
		    itmax_d = 6
		    mode    = 0
c
200		    continue
c
			if(iter_d.ge.itmax_d) goto 300
c
			if(mode.eq.0) then
				if(xnorm_new.lt.xnorm_old) mode = 1
			else
				if(xnorm_new.gt.xnorm_med) goto 300
			endif
c
c			Store old solution and residual
c               	===============================================
c
			call vectcopy(nx,ny,nz,w1,w2)
			call vectcopy(nx,ny,nz,w0,w3)
c
			xnorm_med = xnorm_new
c
c			New damped correction, residual, and its norm
c               	===============================================
c
			call vectcopy(nx,ny,nz,x(iz(1,lev)),w1)
			damp = damp / 2.0d0
			call vectsmult(nx,ny,nz,damp,xtemp(iz(1,lev)),
     1			w1)
c
			call residual(nx,ny,nz,
     1			ipc(iz(5,lev)),rpc(iz(6,lev)),
     2			ac(iz(7,lev)),cc(iz(1,lev)),fc(iz(1,lev)),
     3			w1,w0,rhs(iz(1,lev)))
c
			xnorm_new = xnorm1(nx,ny,nz,w0)
c
			iter_d = iter_d + 1
c
		    goto 200
c
300		    continue
c
c		Accept the previous Newton step
c               ========================================================
c
		    call vectcopy(nx,ny,nz,w2,x(iz(1,lev)))
		    call vectcopy(nx,ny,nz,w3,w0)
c
		    xnorm_new = xnorm_med
		    xnorm_old = xnorm_new
c
c		Determine whether or not to disable damping
c               ========================================================
c
		    if((iter_d-1).eq.0.and.ipkey.le.2) then
			idamp = 0
		    endif
c
		else
c
c		    Damping is disabled. Accept the Newton step
c                   ========================================================
c
		    damp = 1.d0
		    call vectsmult(nx,ny,nz,damp,xtemp(iz(1,lev)),
     1		    x(iz(1,lev)))
c
		    call residual(nx,ny,nz,
     1			ipc(iz(5,lev)),rpc(iz(6,lev)),
     2			ac(iz(7,lev)),cc(iz(1,lev)),fc(iz(1,lev)),
     3			x(iz(1,lev)),w0,rhs(iz(1,lev)))
c
		    xnorm_new = xnorm1(nx,ny,nz,w0)
		    xnorm_old = xnorm_new
c
		endif
c
c		Check stopping criteria
c
		if(iok.ne.0) then
			orsnrm = rsnrm
			if(istop.eq.1) then
				rsnrm = xnorm_new
			else
				write(6,*) 'Bad istop value ....'
			endif
			call print_residual(iok,iters,rsnrm,rsden,
     1			orsnrm)
			if((rsnrm/rsden).le.errtol) goto 400
		endif
c
		if(iters.ge.itmax) goto 400
c
	goto 100
c
400	continue
c
	return
	end
c
c ==========================================================================================
c ==========================================================================================
c 	Get_jacobian.f
c ==========================================================================================
c ==========================================================================================
c
c	Find the Jacobian corresponding to the non linear equation
c
	subroutine get_jacobian(nx,ny,nz,nlev_real,iz,lev,ipkey,x,r,
     1			cprime,rhs,cc)
c
c ==========================================================================================
c 	Declare all variables
c ==========================================================================================
c
c	Input and local variables:
c ==========================================================================================
c
	integer	nx,ny,nz,nlev_real,lev,level
	integer	nxx,nyy,nzz,nxold,nyold,nzold,ipkey
c
	integer	iz(50,*)
c
	real*8	x(*),r(*),cprime(*),rhs(*),cc(*)
c
c ==========================================================================================
c 	Build Jacobian
c ==========================================================================================
c
c	Setup
c ==========================================================================================
c
	nxx = nx
	nyy = ny
	nzz = nz
c
c	Form the rhs of the Newton system: just current residual
c
	call vectcopy(nx,ny,nz,r,rhs(iz(1,lev)))
c
c	Get Non linear part of the Jacobian operator
c
	call dc_vec(cc(iz(1,lev)),x(iz(1,lev)),cprime(iz(1,lev)),
     1		nx,ny,nz,ipkey)
c
c	Build the nlev-1 level operators
c
	do 100 level = lev+1, nlev_real
c
		nxold = nxx
		nyold = nyy
		nzold = nzz
c
		call make_coarse(1,nxold,nyold,nzold,nxx,nyy,nzz)
c
		call restrict(nxold,nyold,nzold,nxx,nyy,nzz,
     1		rhs(iz(1,level-1)),rhs(iz(1,level)))
c
		call restrict(nxold,nyold,nzold,nxx,nyy,nzz,
     1		cprime(iz(1,level-1)),cprime(iz(1,level)))
c
100	continue
c
	return
	end
