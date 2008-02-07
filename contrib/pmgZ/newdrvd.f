c* ///////////////////////////////////////////////////////////////////////////
c* @file    newdrvd.f
c* @author  Michael Holst
c* @brief   Driver for the Newton solver.
c* @version $Id: newdrvd.f 868 2006-04-05 18:32:05Z todd_dolinsky $
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

      subroutine newdriv(iparm,rparm,iwork,rwork,u,
     2   xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf)
c* *********************************************************************
c* purpose:
c*
c*    driver for a screaming inexact-newton-multilevel solver.
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
      integer          k_iz,k_w1,k_w2
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
       call vnmprt(2, 'NEWDRIV:  nx,ny,nz,nlev must be positive...',
     1    48)
         ierror = -1
         iparm(51) = ierror 
         return
      endif
      mxlv = maxlev(nx,ny,nz)
      if (nlev.gt.mxlv) then
         call vnmpri(2, 'NEWDRIV:  max lev for your grid size is: ',
     1      12, mxlv)
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
c*    *** allocate space for two additional work vectors ***
      iretot = iretot + 2*nf
c*
c*    *** some more checks on input ***
      if ((nrwk.lt.iretot) .or. (niwk.lt.iintot)) then
         call vnmpri(2,'NEWDRIV: real    work space must be: ',
     1     45, iretot)
         call vnmpri(2,'NEWDRIV: integer work space must be: ',
     1     45, iintot)
         ierror = -3
         iparm(51) = ierror 
         return
      endif
c*
c*    *** split up the integer work array ***
      k_iz   = 1
      k_ipc  = k_iz   + n_iz
c*
c*    *** split up the real work array ***
      k_rpc  = 1
      k_cc   = k_rpc  + n_rpc
      k_fc   = k_cc   + narr
      k_w1   = k_fc   + narr
      k_w2   = k_w1   + nf
      k_pc   = k_w2   + nf
      k_ac   = k_pc   + 27*narrc
c* ***k_ac_after = 4*nf + 14*narrc
c*
c*    *** call the multigrid driver ***
      call newdriv2(iparm,rparm,nx,ny,nz,u,iwork(k_iz),
     3   rwork(k_w1),rwork(k_w2),
     5   iwork(k_ipc),rwork(k_rpc),
     6   rwork(k_pc),rwork(k_ac),rwork(k_cc),rwork(k_fc),
     8   xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf)
c*
c*    *** return and end ***
      return
      end
      subroutine newdriv2(iparm,rparm,nx,ny,nz,u,iz,
     2   w1,w2,
     3   ipc,rpc,pc,ac,cc,fc,
     4   xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf)
c* *********************************************************************
c* purpose:
c*
c*    this routine uses a newton's method, combined with a linear
c*    multigrid iteration, to solve the following three-dimensional, 
c*    2nd order elliptic partial differential equation:
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
      integer          iparm(*),ipc(*),iz(*)
      double precision rparm(*),rpc(*),pc(*),ac(*),cc(*),fc(*)
      double precision u(*),w1(*),w2(*)
      double precision xf(*),yf(*),zf(*),gxcf(*),gycf(*),gzcf(*)
      double precision a1cf(*),a2cf(*),a3cf(*),ccf(*),fcf(*),tcf(*)
c*
c*    *** misc variables ***
      integer          mgkey,nlev,itmax,iok,iinfo,istop,ipkey,nu1,nu2
      integer          nx,ny,nz,ilev,ido,iters,ierror,nlev_real,ibound
      integer          mgprol,mgcoar,mgsolv,mgdisc,mgsmoo,mode
      double precision epsiln,epsmac,errtol,omegal,omegan
      double precision bf,oh,tsetupf,tsetupc,tsolve
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
      mgprol = iparm(17)
      mgcoar = iparm(18)
      mgdisc = iparm(19)
      mgsmoo = iparm(20)
      mgsolv = iparm(21)
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
      call vtstrt(30, 'NEWDRIV2: fine problem setup', 28)
c*
c*    *** build op and rhs on fine grid ***
      ido = 0
      call buildops (nx,ny,nz,nlev,ipkey,iinfo,ido,iz,
     2   mgprol,mgcoar,mgsolv,mgdisc,
     3   ipc,rpc,pc,ac,cc,fc,
     4   xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf)
c*
c*    *** stop timer ***
      call vtstop(30, 'NEWDRIV2: fine problem setup', 28)
c*
c*    *** start timer ***
      call vtstrt(30, 'NEWDRIV2: coarse problem setup', 30)
c*
c*    *** build op and rhs on all coarse grids ***
      ido = 1
      call buildops (nx,ny,nz,nlev,ipkey,iinfo,ido,iz,
     2   mgprol,mgcoar,mgsolv,mgdisc,
     3   ipc,rpc,pc,ac,cc,fc,
     4   xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf)
c*
c*    *** stop timer ***
      call vtstop(30, 'NEWDRIV2: coarse problem setup', 30)
c*
c* ******************************************************************
c* *** this overwrites the rhs array provided by pde specification
c* ****** compute an algebraically produced rhs for the given tcf ***
      mode = 1
      if ((istop .eq. 4) .or. (istop .eq. 5)) then
         call buildALG (nx,ny,nz,mode,nlev,iz,
     2      ipc,rpc,ac,cc,ccf,tcf,fc,fcf)
      endif
c* ******************************************************************
c*
c*    *** determine machine epsilon ***
      epsiln = epsmac(0)
c*
c*    *** impose zero dirichlet boundary conditions (now in source fcn) ***
      call fboundPMG00(nx,ny,nz,u)
c*
c*    *** MATLAB ***
c*    print*,' new = [ '
c*
c*    *** start timer ***
      call vtstrt(30, 'NEWDRIV2: solve', 15)
c*
c*    *** call specified multigrid method ***
      nlev_real = nlev
      iok  = 1
      ilev = 1
      if (mgkey .eq. 0) then
         call newton(nx,ny,nz,u,iz,ccf,fcf,w1,w2,
     2      istop,itmax,iters,ierror,nlev,ilev,nlev_real,mgsolv,
     3      iok,iinfo,epsiln,errtol,omegan,nu1,nu2,mgsmoo,
     4      a1cf,a2cf,a3cf,
     5      ipc,rpc,pc,ac,cc,fc,tcf)
      else if (mgkey .eq. 1) then
         call fnewton(nx,ny,nz,u,iz,ccf,fcf,w1,w2,
     2      istop,itmax,iters,ierror,nlev,ilev,nlev_real,mgsolv,
     3      iok,iinfo,epsiln,errtol,omegan,nu1,nu2,mgsmoo,
     4      a1cf,a2cf,a3cf,
     5      ipc,rpc,pc,ac,cc,fc,tcf)
      else
         call vnmprt(2,'NEWDRIV2: bad mgkey given ', 26)
      endif
c*
c*    *** stop timer ***
      call vtstop(30, 'NEWDRIV2: solve', 15)
c*
c*    *** MATLAB ***
c*    write(*,100) 'new_sf',tsetupf,'new_sc',tsetupc,
c*   2   'new_st',(tsetupf+tsetupc),'new_so',tsolve
c* 100  format(' ];',4(' ',a7,'=',1pe9.3,';'))
c*
c*    *** restore boundary conditions ***
      ibound = 1
      call fboundPMG(ibound,nx,ny,nz,u,gxcf,gycf,gzcf)
c*
c*    *** return and end ***
      return
      end
