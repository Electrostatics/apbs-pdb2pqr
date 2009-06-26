c* ///////////////////////////////////////////////////////////////////////////
c* @file    mgdrvd.f
c* @author  Michael Holst
c* @brief   Driver for the multigrid routines.
c* @version $Id: mgdrvd.f 868 2006-04-05 18:32:05Z todd_dolinsky $
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

      subroutine mgdriv(iparm,rparm,iwork,rwork,u,
     2   xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf)
c* *********************************************************************
c* purpose:
c*
c*    multilevel solver driver.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
c*
c*    *** input parameters ***
      integer          iparm(*),iwork(*)
      double precision rparm(*),rwork(*),u(*)
      double precision xf(*),yf(*),zf(*),gxcf(*),gycf(*),gzcf(*)
      double precision a1cf(*),a2cf(*),a3cf(*),ccf(*),fcf(*),tcf(*)
c*
c*    *** variables returned from mgsz ***
      integer          nxc,nyc,nzc,nf,nc,narr,narrc,n_rpc
      integer          n_iz,n_ipc,iretot,iintot
c*
c*    *** misc variables ***
      integer          nrwk,niwk,nx,ny,nz,nlev,ierror,maxlev,mxlv
      integer          mgcoar,mgdisc,mgsolv
      integer          k_iz
      integer          k_ipc,k_rpc,k_ac,k_cc,k_fc,k_pc
c*
c*    *** decode some parameters ***
      nrwk   = iparm(1)
      niwk   = iparm(2)
      nx     = iparm(3)
      ny     = iparm(4)
      nz     = iparm(5)
      nlev   = iparm(6)
c*
c*    *** some checks on input ***
      if ((nlev.le.0).or.(nx.le.0).or.(ny.le.0).or.(nz.le.0)) then
         call vnmprt(2,'% MGDRIV:  nx ny nz and nlev ', 29)
         call vnmprt(2,'% MGDRIV:  must be positive...', 30)
         ierror = -1
         iparm(51) = ierror 
         return
      endif
      mxlv = maxlev(nx,ny,nz)
      if (nlev.gt.mxlv) then
         call vnmprt(2,'% MGDRIV:  max levels for your ', 31)
         call vnmpri(2,'% MGDRIV:  grid size is: ',25 ,mxlv)
         ierror = -2
         iparm(51) = ierror 
         return
      endif
c*
c*    *** basic grid sizes, etc. ***
      mgcoar = iparm(18)
      mgdisc = iparm(19)
      mgsolv = iparm(21)
      call mgsz(mgcoar,mgdisc,mgsolv,nx,ny,nz,nlev,nxc,nyc,nzc,
     2   nf,nc,narr,narrc,n_rpc,n_iz,n_ipc,iretot,iintot)
c*
c*    *** some more checks on input ***
      if ((nrwk.lt.iretot) .or. (niwk.lt.iintot)) then
         call vnmpri(2,'% MGDRIV: real    work space must be: ',
     2      38, iretot)
         call vnmpri(2,'% MGDRIV: integer work space must be: ',
     2      38, iintot)
         ierror = -3
         iparm(51) = ierror 
         return
      endif
c*
c*    *** split up the integer work array ***
      k_iz  = 1
      k_ipc = k_iz + n_iz
c*
c*    *** split up the real work array ***
      k_rpc = 1
      k_cc  = k_rpc + n_rpc
      k_fc  = k_cc  + narr
      k_pc  = k_fc  + narr
      k_ac  = k_pc  + 27*narrc
c* ***k_ac_after =  4*nf +  4*narrc
c* ***k_ac_after =  4*nf + 14*narrc
c* ***k_ac_after = 14*nf + 14*narrc
c*
c*    *** call the multigrid driver ***
      call mgdriv2(iparm,rparm,nx,ny,nz,u,iwork(k_iz),
     2   iwork(k_ipc),rwork(k_rpc),
     3   rwork(k_pc),rwork(k_ac),rwork(k_cc),rwork(k_fc),
     4   xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf)
c*
c*    *** return and end ***
      return
      end
      subroutine mgdriv2(iparm,rparm,nx,ny,nz,u,iz,
     2   ipc,rpc,pc,ac,cc,fc,
     3   xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf)
c* *********************************************************************
c* purpose:
c*
c*    this routine uses a multigrid method to solve the following
c*    three-dimensional, 2nd order elliptic partial differential
c*    equation:
c*
c*         lu = f, u in omega
c*          u = g, u on boundary of omega
c*    where
c*
c*         omega = [xmin,xmax]x[ymin,ymax]x[zmin,zmax]
c*
c*    the multigrid code requires the operator in the form:
c*
c*         - \nabla \cdot (a \nabla u) + c(u) = f
c*
c*    with
c*
c*        a(x,y,z),f(x,y,z), scalar functions (possibly discontinuous) 
c*        on omega.  (discontinuities must be along fine grid lines).
c*        boundary function g(x,y,z) is smooth on boundary of omega.
c*
c*        the function c(u) is a possibly nonlinear function of the
c*        unknown u, and varies (possibly discontinuously) with the
c*        spatial position also.
c*
c* user inputs:
c*
c*    the user must provide the coefficients of the differential
c*    operator, some initial parameter settings in an integer and a
c*    real parameter array, and various work arrays.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
c*
c*    *** input parameters ***
      integer          iparm(*),ipc(*),iz(50,*)
      double precision rparm(*),rpc(*),pc(*),ac(*),cc(*),fc(*)
      double precision u(*)
      double precision xf(*),yf(*),zf(*),gxcf(*),gycf(*),gzcf(*)
      double precision a1cf(*),a2cf(*),a3cf(*),ccf(*),fcf(*),tcf(*)
c*
c*    *** misc variables ***
      integer          mgkey,nlev,itmax,iok,iinfo,istop,ipkey,nu1,nu2
      integer          nx,ny,nz,ilev,ido,iters,ierror,nlev_real,ibound
      integer          mgprol,mgcoar,mgsolv,mgdisc,mgsmoo,iperf,mode
      double precision epsiln,epsmac,errtol,omegal,omegan
      double precision bf,oh,tsetupf,tsetupc,tsolve
c*
c*    *** more misc variables ***
      integer          itmax_p,iters_p,iok_p,iinfo_p
      double precision errtol_p,rho_p
      double precision rho_min,rho_max,rho_min_mod,rho_max_mod
      integer          nxf,nyf,nzf,nxc,nyc,nzc,level,nlevd
c*
c*    *** decode the iparm array ***
      nlev   = iparm(6)
      nu1    = iparm(7)
      nu2    = iparm(8)
      mgkey  = iparm(9)
      itmax  = iparm(10)
      istop  = iparm(11)
      iinfo  = iparm(12)
      ipkey  = iparm(14)
      mode   = iparm(16)
      mgprol = iparm(17)
      mgcoar = iparm(18)
      mgdisc = iparm(19)
      mgsmoo = iparm(20)
      mgsolv = iparm(21)
      iperf  = iparm(22)
      errtol = rparm(1)
      omegal = rparm(9)
      omegan = rparm(10)
c*
c*    *** intitialize the iteration timer ***
      call prtstp(0,-99,0.0d0,0.0d0,0.0d0)
c*
c*    *** build the multigrid data structure in iz ***
      call buildstr (nx,ny,nz,nlev,iz)
c*
c*    *** start timer ***
      call vtstrt(30, 'MGDRIV2: fine problem setup', 17)
c*
c*    *** build op and rhs on fine grid ***
      ido = 0
      call buildops (nx,ny,nz,nlev,ipkey,iinfo,ido,iz,
     2   mgprol,mgcoar,mgsolv,mgdisc,
     3   ipc,rpc,pc,ac,cc,fc,
     4   xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf)
c*
c*    *** stop timer ***
      call vtstop(30, 'MGDRIV2: fine problem setup', 17)
c*
c*    *** start timer ***
      call vtstrt(30, 'MGDRIV2: coarse problem setup', 19)
c*
c*    *** build op and rhs on all coarse grids ***
      ido = 1
      call buildops (nx,ny,nz,nlev,ipkey,iinfo,ido,iz,
     2   mgprol,mgcoar,mgsolv,mgdisc,
     3   ipc,rpc,pc,ac,cc,fc,
     4   xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf)
c*
c*    *** stop timer ***
      call vtstop(30, 'MGDRIV2: coarse problem setup', 19)
c*
c*    *** determine machine epsilon ***
      epsiln = epsmac(0)
c*
c* ******************************************************************
c* *** analysis ***
c* *** note: we destroy the rhs function "fc" here in "mpower" ***
c* ******************************************************************
c*
c*    *** errtol and itmax ***
      itmax_p   = 1000
      iok_p     = 0
      nlev_real = nlev
      nlevd     = nlev_real
c*
c*    *** finest level initialization ***
      nxf  = nx
      nyf  = ny
      nzf  = nz
c*
c*    *** go down grids: compute max/min eigenvalues of all operators ***
      do 40 level = 1, nlev_real
         nlevd = nlev_real - level + 1
c*
c*       *** move down the grids ***
         if (level .ne. 1) then
c*
c*          *** find new grid size ***
            call mkcors(1,nxf,nyf,nzf,nxc,nyc,nzc)
c*
c*          *** new grid size ***
            nxf = nxc
            nyf = nyc
            nzf = nzc
         endif
c*       print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
c*   2         ,'%%%%%%%%%%'
c*       write(6,200)'% MGDRIV2: ANALYSIS==>',nxf,nyf,nzf
 200     format(a,' [',i3,',',i3,',',i3,'] ')
c*
c*       *** largest eigenvalue of the system matrix A ***
         if ((iperf .eq. 1) .or. (iperf .eq. 3)) then
            call vnmprt(2,'% MGDRIV2: power calculating rho(A)...',38)
            iters_p   = 0
            iinfo_p   = iinfo
            errtol_p  = 1.0e-4
            call power(nxf,nyf,nzf,iz,level,
     2         ipc,rpc,ac,cc,
     3         a1cf,a2cf,a3cf,ccf,
     4         rho_max,rho_max_mod,errtol_p,itmax_p,iters_p,iinfo_p)
            call vnmpri(2,'% MGDRIV2: power iters   = ',27,iters_p)
            call vnmprd(2,'% MGDRIV2: power eigmax  = ',27,rho_max)
            call vnmprd(2,'% MGDRIV2: power (MODEL) = ',27,rho_max_mod)
c*
c*          *** smallest eigenvalue of the system matrix A ***
            call vnmprt(2,'% MGDRIV2: ipower calculating ', 30)
            call vnmprt(2,'% MGDRIV2: lambda_min(A)...', 27)
            iters_p   = 0
            iinfo_p   = iinfo
            errtol_p  = 1.0e-4
            call azeros(nxf,nyf,nzf,u)
            call ipower(nxf,nyf,nzf,u,iz,
     2         a1cf,a2cf,a3cf,ccf,fcf,
     3         rho_min,rho_min_mod,errtol_p,itmax_p,iters_p,
     4         nlevd,level,nlev_real,mgsolv,
     5         iok_p,iinfo_p,epsiln,errtol,omegal,nu1,nu2,mgsmoo,
     6         ipc,rpc,pc,ac,cc,tcf)
            call vnmpri(2,'% MGDRIV2: ipower iters   = ',28,iters_p)
            call vnmprd(2,'% MGDRIV2: ipower eigmin  = ',28,rho_min)
            call vnmprd(2,'% MGDRIV2: ipower (MODEL) = ',28,rho_min_mod)
c*
c*          *** condition number estimate ***
            call vnmprd(2,'% MGDRIV2: condition number  = ',
     2         31, rho_max/rho_min)
            call vnmprd(2,'% MGDRIV2: condition (MODEL) = ',
     2         31, rho_max_mod/rho_min_mod)
         endif
c*
c*       *** spectral radius of the multigrid operator M ***
c*       *** NOTE: due to lack of vectors, we destroy "fc" in mpower... ***
         if ((iperf .eq. 2) .or. (iperf .eq. 3)) then
            call vnmprt(2,'% MGDRIV2: mpower calculating rho(M)...',39)
            iters_p   = 0
            iinfo_p   = iinfo
            errtol_p  = epsiln
            call azeros(nxf,nyf,nzf,u(iz(1,level)))
            call mpower(nxf,nyf,nzf,u,iz,
     2         a1cf,a2cf,a3cf,ccf,fcf,
     3         rho_p,errtol_p,itmax_p,iters_p,
     4         nlevd,level,nlev_real,mgsolv,
     5         iok_p,iinfo_p,epsiln,errtol,omegal,nu1,nu2,mgsmoo,
     6         ipc,rpc,pc,ac,cc,fc,tcf)
            call vnmpri(2,'% MGDRIV2: mpower iters  = ',27,iters_p)
            call vnmprd(2,'% MGDRIV2: mpower rho(M) = ',27,rho_p)
         endif
c*
c*       *** reinitialize the solution function ***
         call azeros(nxf,nyf,nzf,u(iz(1,level)))
c*
c*    *** next grid ***
 40   continue
c*    print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
c*
c*    *** reinitialize the solution function ***
      call azeros(nx,ny,nz,u)
c*
c* ******************************************************************
c* *** this overwrites the rhs array provided by pde specification
c* ****** compute an algebraically produced rhs for the given tcf ***
      if ((istop .eq. 4) .or. (istop .eq. 5) .or. (iperf .ne. 0)) then
         call vnmprt(2,'% MGDRIV2: generating algebraic ', 32)
         call vnmprt(2,'% MGDRIV2: RHS from your soln...', 32)
         call buildALG (nx,ny,nz,mode,nlev,iz,
     2      ipc,rpc,ac,cc,ccf,tcf,fc,fcf)
      endif
c* ******************************************************************
c*
c*    *** impose zero dirichlet boundary conditions (now in source fcn) ***
      call fboundPMG00(nx,ny,nz,u)
c*
c*    *** MATLAB ***
c*    print*,' mg = [ '
c*
c*    *** start timer ***
      call vtstrt(30, 'MGDRIV2: solve', 14)
c*
c*    *** call specified multigrid method ***
      if ((mode .eq. 0) .or. (mode .eq. 2)) then
c*       print*,'% MGDRIV2: linear mode...'
         nlev_real = nlev
         iok  = 1
         ilev = 1
         if (mgkey .eq. 0) then
            call mvcs(nx,ny,nz,u,iz,a1cf,a2cf,a3cf,ccf,
     2         istop,itmax,iters,ierror,nlev,ilev,nlev_real,mgsolv,
     3         iok,iinfo,epsiln,errtol,omegal,nu1,nu2,mgsmoo,
     4         ipc,rpc,pc,ac,cc,fc,tcf)
         elseif (mgkey .eq. 1) then
            call fmvcs(nx,ny,nz,u,iz,a1cf,a2cf,a3cf,ccf,
     2         istop,itmax,iters,ierror,nlev,ilev,nlev_real,mgsolv,
     3         iok,iinfo,epsiln,errtol,omegal,nu1,nu2,mgsmoo,
     4         ipc,rpc,pc,ac,cc,fc,tcf)
         else
            call vnmprt(2,'% MGDRIV2: bad mgkey given... ',30)
         endif
      endif
      if ((mode .eq. 1) .or. (mode .eq. 2)) then
c*       print*,'% MGDRIV2: nonlinear mode...'
         nlev_real = nlev
         iok  = 1
         ilev = 1
         if (mgkey .eq. 0) then
            call mvfas(nx,ny,nz,u,iz,a1cf,a2cf,a3cf,ccf,fcf,
     2         istop,itmax,iters,ierror,nlev,ilev,nlev_real,mgsolv,
     3         iok,iinfo,epsiln,errtol,omegan,nu1,nu2,mgsmoo,
     4         ipc,rpc,pc,ac,cc,fc,tcf)
         else if (mgkey .eq. 1) then
            call fmvfas(nx,ny,nz,u,iz,a1cf,a2cf,a3cf,ccf,fcf,
     2         istop,itmax,iters,ierror,nlev,ilev,nlev_real,mgsolv,
     3         iok,iinfo,epsiln,errtol,omegan,nu1,nu2,mgsmoo,
     4         ipc,rpc,pc,ac,cc,fc,tcf)
         else
            call vnmprt(2,'% MGDRIV2: bad mgkey given... ',30)
         endif
      endif
c*
c*    *** stop timer ***
      call vtstop(30, 'MGDRIV2: solve', 14)
c*
c*    *** MATLAB ***
c*    write(*,100) 'mg_sf',tsetupf,'mg_sc',tsetupc,
c*   2   'mg_st',(tsetupf+tsetupc),'mg_so',tsolve
 100  format(' ];',4(' ',a7,'=',1pe9.3,';'))
c*
c*    *** restore boundary conditions ***
      ibound = 1
      call fboundPMG(ibound,nx,ny,nz,u,gxcf,gycf,gzcf)
c*
c*    *** return and end ***
      return
      end
      subroutine mgsz(mgcoar,mgdisc,mgsolv,nx,ny,nz,nlev,nxc,nyc,nzc,
     2   nf,nc,narr,narrc,n_rpc,n_iz,n_ipc,iretot,iintot)
c* *********************************************************************
c* purpose:
c*
c*   this routine computes the required sizes of the real and integer
c*   work arrays for the multigrid code.  these two sizes are a
c*   (complicated) function of input parameters.
c*   
c*   the work arrays must have been declared in the calling program as:
c*
c*       double precision rwork(iretot)
c*       integer          iwork(iintot)
c*
c*   where:
c*
c*       iretot   = function_of(mgcoar,mgdisc,mgsolv,nx,ny,nz,nlev)
c*       iintot   = function_of(mgcoar,mgdisc,mgsolv,nx,ny,nz,nlev)
c*
c*       mgcoar   = coarsening technique:
c*                  0=standard discretization
c*                  1=averaged coefficient + standard discretization
c*                  2=algebraic galerkin coarsening
c*
c*       mgdisc   = discretization technique:
c*                  0=box method
c*                  1=fem method
c*
c*       mgsolv   = coarse grid solver:
c*                  0=conjugate gradients
c*                  1=symmetric banded linpack solver
c*
c*       nx,ny,nz = grid dimensions in each direction, 
c*                  including boundary points
c*
c*       nlev     = the number of multigrid levels desired for the 
c*                  method.
c*
c*   other parameters:
c*
c*       nf       = number of unknowns on the finest mesh
c*                = nx * ny * nz
c*
c*       nc       = number of unknowns on the coarsest mesh
c*
c*       narr     = storage for one vector on all the meshes
c*
c*       narrc    = storage for one vector on all the meshes but the finest
c*
c*   the work arrays rwork and iwork will be chopped into smaller 
c*   pieces according to:
c*
c*       double precision ac(STORE)         (system operators on all levels)
c*       double precision pc(27*narrc)      (prol. opers for coarse levels)
c*       double precision cc(narr),fc(narr) (helmholtz term, rhs -- all levels)
c*       double precision rpc(100*(nlev+1)) (real info for all levels)
c*       integer          ipc(100*(nlev+1)) (integer info for all levels)
c*       integer          iz(50,nlev+1),    (pointers into ac,pc,cc,fc,etc.)
c*
c*   where STORE depends on the discretization, coarsening, and coarse
c*   grid solver:
c*
c*       STORE =  4*nf +  4*narrc + NBAND*nc (mgdisc=box, mgcoar=stan/harm)
c*          or =  4*nf + 14*narrc + NBAND*nc (mgdisc=box, mgcoar=gal)
c*          or = 14*nf + 14*narrc + NBAND*nc (mgdisc=fem, mgcoar=stan/harm/gal)
c*
c*       NBAND = 0                           (mgsolv=iterative)
c*          or = 1+(nxc-2)*(nyc-2)           (mgsolv=7-pt banded linpack)
c*          or = 1+(nxc-2)*(nyc-2)+(nxc-2)+1 (mgsolv=27-pt banded linpack)
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
c*
c*    *** input parameters ***
      integer          mgcoar,mgdisc,mgsolv,nx,ny,nz,nlev
c*
c*    *** parameters: num of different types of arrays in mg code ***
      integer          num_nf,num_narr,num_narrc
      parameter        (num_narr  = 2)
      parameter        (num_nf    = 0)
      parameter        (num_narrc = 27)
c*
c*    *** variables returned from mgsz ***
      integer          nxc,nyc,nzc,nf,nc,narr,narrc,n_rpc
      integer          n_iz,n_ipc,iretot,iintot
c*
c*    *** misc variables ***
      integer          nc_band,num_band,n_band,nxf,nyf,nzf,level
      integer          num_nf_oper,num_narrc_oper
c*
c*    *** go down grids: compute max/min eigenvalues of all operators ***
	  nf   = nx * ny * nz
      narr = nf 
      nxf  = nx
      nyf  = ny
      nzf  = nz
      nxc  = nx
      nyc  = ny
      nzc  = nz
	  
      do 10 level = 2, nlev
c*
c*       *** find new grid size ***
         call mkcors(1,nxf,nyf,nzf,nxc,nyc,nzc)
c*
c*       *** new grid size ***
         nxf = nxc
         nyf = nyc
         nzf = nzc
c*
c*       *** add the unknowns on this level to the total ***
         narr = narr + (nxf * nyf * nzf)
 10   continue
      nc = nxc * nyc * nzc
      narrc = narr - nf
c*
c*    *** box or fem on fine grid? ***
      if (mgdisc .eq. 0) then
         num_nf_oper = 4
      elseif (mgdisc .eq. 1) then
         num_nf_oper = 14
      else
         call vnmpri(2,'% MGSZ: invalid mgdisc parameter... ',
     2      36, mgdisc)
      endif
c*
c*    *** galerkin or standard coarsening? ***
      if (((mgcoar .eq. 0) .or. (mgcoar .eq. 1))
     2      .and. (mgdisc .eq. 0)) then
         num_narrc_oper = 4
      elseif (mgcoar .eq. 2) then
         num_narrc_oper = 14
      else
         call vnmprt(2,'% MGSZ: invalid mgcoar parameter...',35)
      endif
c*
c*    *** symmetric banded linpack storage on coarse grid ***
      if (mgsolv .eq. 0) then
         n_band = 0
      elseif (mgsolv .eq. 1) then
         if (((mgcoar .eq. 0) .or. (mgcoar .eq. 1))
     2      .and. (mgdisc .eq. 0)) then
            num_band = 1 + (nxc-2)*(nyc-2)
         else
            num_band = 1 + (nxc-2)*(nyc-2) + (nxc-2) + 1
         endif
         nc_band = (nxc-2)*(nyc-2)*(nzc-2)
         n_band  = nc_band * num_band
      else
         call vnmprt(2,'% MGSZ: invalid mgsolv parameter...',35)
      endif
c*
c*    *** info work array required storage ***
      n_rpc  = 100*(nlev+1)
c*
c*    *** resulting total required real storage for method ***
      iretot = num_narr*narr 
     2       + (num_nf    + num_nf_oper)*nf 
     3       + (num_narrc + num_narrc_oper)*narrc 
     4       + n_band
     5       + n_rpc
c*
c*    *** the integer storage parameters ***
      n_iz  = 50*(nlev+1)
      n_ipc = 100*(nlev+1)
c*
c*    *** resulting total required integer storage for method ***
      iintot = n_iz + n_ipc
c*
c*    *** return and end ***
      return
      end
      subroutine mgsize(mgcoar,mgdisc,mgsolv,nx,ny,nz,nlev)
c* *********************************************************************
c* purpose:
c*
c*   this routine prints out the required sizes of the real and integer
c*   work arrays for the multigrid code.  these two sizes are a
c*   (complicated) function of the four input parameters:
c*
c*        nx,ny,nz ==> number of grid points in each grid direction.
c*        nlev     ==> number of multigrid levels to employ.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
c*
c*    *** input parameters ***
      integer          mgcoar,mgdisc,mgsolv,nx,ny,nz,nlev
c*
c*    *** variables returned from mgsz ***
      integer          nxc,nyc,nzc,nf,nc,narr,narrc,n_rpc
      integer          n_iz,n_ipc,iretot,iintot
c*
c*    *** basic grid sizes, etc. ***
      call mgsz(mgcoar,mgdisc,mgsolv,nx,ny,nz,nlev,nxc,nyc,nzc,
     2   nf,nc,narr,narrc,n_rpc,n_iz,n_ipc,iretot,iintot)
c*
c*    *** resulting total required for real/integer storage ***
c*    print*,'% MGSIZE: number of unknowns on finest level:    ',nf
c*    print*,'% MGSIZE: number of unknowns on coarsest level:  ',nc
c*    print*,'% MGSIZE: storage for a vector on all levels:    ',narr
c*    print*,'% MGSIZE: storage for a vector on coarse levels: ',narrc
c*    print*,'% MGSIZE: REQUIRED floating point array size:    ',iretot
c*    print*,'% MGSIZE: REQUIRED integer array size:           ',iintot
c*
c*    *** return and end ***
      return
      end
