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
c	Smoothing.f
c
c	Version 1:      7/4/07
c
c	Author: Patrice Koehl (in collaboration with Marc Delarue)
c
c	This file contains a set of routines for linear smoothing
c
c ==========================================================================================
c ==========================================================================================
c
c
c ==========================================================================================
c ==========================================================================================
c	smoothing.f
c
c	This is the driver routine for the smoothing procedure. We only implement the fast
c	diagonal iterative method
c ==========================================================================================
c ==========================================================================================
c
	subroutine smoothing(nx,ny,nz,ipc,rpc,ac,cc,fc,x,w1,w2,r,itmax,
     1		iters,errtol,omega,iresid,iadjoint,meth)
c
c ==========================================================================================
c	Declaration
c ==========================================================================================
c
	integer	nx,ny,nz,itmax,iters,iresid,iadjoint,meth
c
	integer	ipc(*)
c
	real*8	omega,errtol
	real*8	rpc(*),ac(nx*ny*nz,*),cc(nx,ny,nz),fc(nx,ny,nz)
	real*8	x(nx,ny,nz),w1(nx,ny,nz),w2(nx,ny,nz)
	real*8	r(nx,ny,nz)
c
c ==========================================================================================
c	Call the appropriate smoothing routine
c ==========================================================================================
c
	if(meth.eq.1) then
		call GaussSeidel_RB(nx,ny,nz,ipc,rpc,ac,cc,fc,x,
     1		w1,w2,r,itmax,iters,errtol,omega,iresid,iadjoint)
	elseif(meth.eq.4) then
		call conjgrad_lin(nx,ny,nz,ipc,rpc,ac,cc,fc,x,
     1		w1,w2,r,itmax,iters,errtol,omega,iresid,iadjoint)
	else
		write(6,*) 'Bad smoothing method'
		stop
	endif
c
	return
	end
c
c ==========================================================================================
c ==========================================================================================
c	GaussSeidel_RB.f
c
c	This is the driver routine for the Gauss Seidel red-black smoothing procedure. 
c ==========================================================================================
c ==========================================================================================
c
	subroutine GaussSeidel_RB(nx,ny,nz,ipc,rpc,ac,cc,fc,x,w1,w2,r,
     1		itmax,iters,errtol,omega,iresid,iadjoint)
c
c ==========================================================================================
c	Declaration
c ==========================================================================================
c
	integer	nx,ny,nz,itmax,iters,iresid,iadjoint
	integer	numdia
c
	integer	ipc(*)
c
	real*8	omega,errtol
	real*8	rpc(*),ac(nx*ny*nz,*),cc(nx,ny,nz),fc(nx,ny,nz)
	real*8	x(nx,ny,nz),w1(nx,ny,nz),w2(nx,ny,nz)
	real*8	r(nx,ny,nz)
c
c ==========================================================================================
c	Call the proper routine based on Stencil size
c ==========================================================================================
c
	numdia = ipc(11)
c
	if(numdia.eq.7) then
		call GaussSeidel7(nx,ny,nz,ipc,rpc,ac(1,1),cc,fc,
     1		ac(1,2),ac(1,3),ac(1,4),
     2		x,w1,w2,r,
     3		itmax,iters,errtol,omega,iresid,iadjoint)
	elseif(numdia.eq.27) then
		call GaussSeidel27(nx,ny,nz,ipc,rpc,ac(1,1),cc,fc,
     1		ac(1,2),ac(1,3),ac(1,4),ac(1,5),ac(1,6),ac(1,7),ac(1,8),
     1		ac(1,9),ac(1,10),ac(1,11),ac(1,12),ac(1,13),ac(1,14),
     2		x,w1,w2,r,
     3		itmax,iters,errtol,omega,iresid,iadjoint)
	else
		write(6,*) 'In GaussSeidel_RB : bad stencil type...'
	endif
c
	return
	end
c
c ==========================================================================================
c ==========================================================================================
c	GaussSeidel7.f
c
c	7-diagonal Gauss Seidel routine
c ==========================================================================================
c ==========================================================================================
c
	subroutine GaussSeidel7(nx,ny,nz,ipc,rpc,
     1			oC,cc,fc,oE,oN,uC,
     2			x,w1,w2,r,itmax,iters,errtol,omega,iresid,
     3			iadjoint)
c
c ==========================================================================================
c	Declaration
c ==========================================================================================
c
	integer	nx,ny,nz,itmax,iters,iresid,iadjoint
	integer	i,j,k,i1,i2,j1,j2,k1,k2,istep
c
	integer	ipc(*)
c
	real*8	omega,errtol
c
	real*8	rpc(*)
	real*8	oE(nx,ny,nz),oN(nx,ny,nz),uC(nx,ny,nz)
	real*8	fc(nx,ny,nz),oC(nx,ny,nz),cc(nx,ny,nz)
	real*8	x(nx,ny,nz),w1(nx,ny,nz),w2(nx,ny,nz)
	real*8	r(nx,ny,nz)
c
c ==========================================================================================
c	Do Gauss Seidel iterations itmax times
c ==========================================================================================
c
	i1 = (1-iadjoint) * 2 + iadjoint	*(nx-1)
	i2 = iadjoint     * 2 + (1-iadjoint)	*(nx-1)
	j1 = (1-iadjoint) * 2 + iadjoint	*(ny-1)
	j2 = iadjoint     * 2 + (1-iadjoint)	*(ny-1)
	k1 = (1-iadjoint) * 2 + iadjoint	*(nz-1)
	k2 = iadjoint     * 2 + (1-iadjoint)	*(nz-1)
c
	istep = iadjoint*(-1) + (1-iadjoint)*(1)
c
	do 400 iters = 1,itmax
c
		do 300 k=k1,k2,istep
		   do 200 j = j1,j2,istep
			do 100 i = i1,i2,istep
c
				x(i,j,k) = (fc(i,j,k) + (
     1				+ oN(i,j,k)	* x(i,j+1,k)
     2				+ oN(i,j-1,k)	* x(i,j-1,k)
     3				+ oE(i,j,k)	* x(i+1,j,k)
     4				+ oE(i-1,j,k)	* x(i-1,j,k)
     5				+ uC(i,j,k)	* x(i,j,k+1)
     6				+ uC(i,j,k-1)	* x(i,j,k-1)
     7				)) /(oC(i,j,k) + cc(i,j,k))
c
100			continue
200		   continue
300		continue
c
400	continue
c
c	If needed, return the new residual
c
	if(iresid.eq.1) then
		call mresidual7(nx,ny,nz,ipc,rpc,oC,cc,fc,oE,oN,uC,x,r)
	endif
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c	GaussSeidel27.f
c
c	27-diagonal Gauss Seidel routine
c ==========================================================================================
c ==========================================================================================
c
	subroutine GaussSeidel27(nx,ny,nz,ipc,rpc,
     1			oC,cc,fc,oE,oN,uC,
     2			oNE,oNW,uE,uW,uN,uS,uNE,uNW,uSE,uSW,
     2			x,w1,w2,r,itmax,iters,errtol,omega,iresid,
     3			iadjoint)
c
c ==========================================================================================
c	Declaration
c ==========================================================================================
c
	integer	nx,ny,nz,itmax,iters,iresid,iadjoint
	integer	i,j,k,i1,i2,j1,j2,k1,k2,istep
c
	integer	ipc(*)
c
	real*8	omega,errtol
c
	real*8	tmp0,tmpU,tmpD
	real*8	rpc(*)
	real*8	oE(nx,ny,nz),oN(nx,ny,nz),uC(nx,ny,nz)
	real*8	oNE(nx,ny,nz),oNW(nx,ny,nz),uE(nx,ny,nz)
	real*8	uW(nx,ny,nz),uN(nx,ny,nz),uS(nx,ny,nz)
	real*8	uNE(nx,ny,nz),uNW(nx,ny,nz),uSE(nx,ny,nz)
	real*8	uSW(nx,ny,nz)
	real*8	fc(nx,ny,nz),oC(nx,ny,nz),cc(nx,ny,nz)
	real*8	x(nx,ny,nz),w1(nx,ny,nz),w2(nx,ny,nz)
	real*8	r(nx,ny,nz)
c
c ==========================================================================================
c	Do Gauss Seidel iterations itmax times
c ==========================================================================================
c
	i1 = (1-iadjoint) * 2 + iadjoint	*(nx-1)
	i2 = iadjoint     * 2 + (1-iadjoint)	*(nx-1)
	j1 = (1-iadjoint) * 2 + iadjoint	*(ny-1)
	j2 = iadjoint     * 2 + (1-iadjoint)	*(ny-1)
	k1 = (1-iadjoint) * 2 + iadjoint	*(nz-1)
	k2 = iadjoint     * 2 + (1-iadjoint)	*(nz-1)
c
	istep = iadjoint*(-1) + (1-iadjoint)*(1)
c
	do 400 iters = 1,itmax
c
		do 300 k=k1,k2,istep
		   do 200 j = j1,j2,istep
			do 100 i = i1,i2,istep
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
				x(i,j,k) = (fc(i,j,k)+tmp0+tmpU+tmpD)/
     1				(oC(i,j,k) + cc(i,j,k))
c
100			continue
200		    continue
300		continue
c
400	continue
c
c	If needed, return the new residual
c
	if(iresid.eq.1) then
		call mresidual27(nx,ny,nz,ipc,rpc,oC,cc,fc,
     1			oE,oN,uC,oNE,oNW,uE,uW,uN,uS,uNE,uNW,uSE,uSW,
     3			x,r)
	endif
c
	return
	end
c ==========================================================================================
c ==========================================================================================
c	conjgrad_lin.f
c
c	This routine solves the system of equation Ax=b using conjugate gradient
c ==========================================================================================
c ==========================================================================================
c
	subroutine conjgrad_lin(nx,ny,nz,ipc,rpc,ac,cc,fc,x,p,ap,r,
     1		itmax,iters,errtol,omega,iresid,iadjoint)
c
c ==========================================================================================
c	Declaration
c ==========================================================================================
c
	integer	nx,ny,nz,itmax,iters,iresid,iadjoint
	integer	i,j,k
c
	integer	ipc(*)
c
	real*8	omega,errtol,rsnrm,pAp,denom,xnorm2
	real*8	rhok1,rhok2,alpha,beta,dotvect
c
	real*8	rpc(*),ac(nx,ny,nz,*),fc(nx,ny,nz),cc(nx,ny,nz)
	real*8	x(nx,ny,nz),r(nx,ny,nz),p(nx,ny,nz),ap(nx,ny,nz)
c
c ==========================================================================================
c	Setup
c ==========================================================================================
c
	iters = 0
	if(iters.ge.itmax.and.iresid.eq.0) return
c
	call mresidual(nx,ny,nz,ipc,rpc,ac,cc,fc,x,r)
c
	denom = xnorm2(nx,ny,nz,r)
c
	if(denom.eq.0.d0) return
	if(iters.ge.itmax) return
c
c ==========================================================================================
c	ConjGrad loop
c ==========================================================================================
c
100	continue
c
c		Compute current stopping criteria
c
		rhok2 = dotvect(nx,ny,nz,r,r)
		rsnrm = dsqrt(rhok2)
c
		if(rsnrm/denom.le.errtol) goto 200
		if(iters.ge.itmax) goto 200
c
c		Form new direction vector
c
		if(iters.eq.0) then
			call vectcopy(nx,ny,nz,r,p)
		else
			beta = rhok2/rhok1
			call vectsmult(nx,ny,nz,(1.d0/beta),r,p)
			call scalevect(nx,ny,nz,beta,p)
		endif
c
c		Find alpha that minimizes error
c
		call matvec_band(nx,ny,nz,ipc,rpc,ac,cc,p,ap)
		pAp = dotvect(nx,ny,nz,p,ap)
		alpha = rhok2/pAp
c
c		Save rhok2
c
		rhok1 = rhok2
c
c		Update solution in direction p of length alpha
c
		call vectsmult(nx,ny,nz,alpha,p,x)
c
c		Update residual
c
		call vectsmult(nx,ny,nz,(-alpha),ap,r)
c
		iters = iters + 1
c
		goto 100
c
200	continue
c
	return
	end
