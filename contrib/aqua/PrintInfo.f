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
c	PrintInfo.f
c       (Information while running PDE solver, based on Michael Holst MG package)
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
c 	Print_ini.f
c ==========================================================================================
c ==========================================================================================
c
c	This subroutine prints out information at the beginning of refinement
c
	subroutine print_ini
c
c ==========================================================================================
c 	Declare all variables
c ==========================================================================================
c
c	Input variables:
c ==========================================================================================
c
	integer	l_log
c
c	Local variables
c ==========================================================================================
c
	integer	i
c
	character*80 header
	character*80 line1
	character*80 line2
c
c	Common blocks
c ==========================================================================================
c
c	Log file for the run
c
	common /log_file/	l_log
c
1	format(a80)
c
c ==========================================================================================
c	Build output
c ==========================================================================================
c
	do 100 i = 1,80
		header(i:i) = '='
100	continue
c
	line1(1:20)  = 'Iteration count     '
	line1(21:40) = 'Relative residual   '
	line1(41:60) = 'Contraction number  '
	line1(61:80) = 'Cumulative CPU time '
c
	line2(1:20)  = '                    '
	line2(21:40) = '  (1-norm)          '
	line2(41:60) = '                    '
	line2(61:80) = '  (seconds)         '
c
	write(6,1) header
	write(6,1) line1
	write(6,1) line2
	write(6,1) header
c
	write(l_log,1) header
	write(l_log,1) line1
	write(l_log,1) line2
	write(l_log,1) header
c
	return
	end
c
c
c ==========================================================================================
c ==========================================================================================
c 	Print_residual.f
c ==========================================================================================
c ==========================================================================================
c
c	This subroutine prints out information during iterations
c
	subroutine print_residual(iok,iters,resid_nrm,residc,old_resid)
c
c ==========================================================================================
c 	Declare all variables
c ==========================================================================================
c
c	Input variables:
c ==========================================================================================
c
	integer	iok,iters
c
	real*8	resid_nrm,residc,old_resid
c
c	Local variables
c ==========================================================================================
c
	real*8	relative_resid,contraction
	real*8	cpu
c
	real*4	start,tend
c
	save start,tend
c
c	Log file for the run
c ==========================================================================================
c
	integer	l_log
c
	common /log_file/ l_log
c
1	format(2x,i5,15x,1pe12.5,8x,1pe12.5,8x,1pe8.2)
2	format('%',i5,15x,1pe12.5,8x,1pe12.5,8x,1pe8.2)
c
c ==========================================================================================
c 	Build output
c ==========================================================================================
c
c	Initialize if number of iterations set to -99
c ==========================================================================================
c
	if(iters.eq.-99) then
c
		call cpu_time(start)
		cpu = 0.d0
c
	elseif(iters.eq.-1) then
c
		call cpu_time(tend)
		cpu = tend-start
		if(iok.eq.1) then
			write(6,1) -1,0.d0,0.d0,cpu
			write(l_log,1) -1,0.d0,0.d0,cpu
		endif
c
	else
c
		call cpu_time(tend)
		cpu = tend-start
c
		if(residc.eq.0.d0) then
			relative_resid = 1.e10
		else
			relative_resid = resid_nrm/residc
		endif
c
		if(old_resid.eq.0) then
			contraction = 1.e10
		else
			contraction = resid_nrm/old_resid
		endif
c
		if(iok.eq.1) then
			write(6,1) iters,relative_resid,contraction,cpu
			write(l_log,1) iters,relative_resid,contraction,
     1					cpu
c		elseif(iok.eq.2) then
c			write(6,2) iters,relative_resid,contraction,cpu
		endif
c
	endif
c
	return
	end
