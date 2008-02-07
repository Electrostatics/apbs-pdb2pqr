c* ///////////////////////////////////////////////////////////////////////////
c* @file    mgcsd.f
c* @author  Michael Holst
c* @brief   The core linear (correction scheme) multigrid routines.
c* @version $Id: mgcsd.f 868 2006-04-05 18:32:05Z todd_dolinsky $
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

      subroutine fmvcs(nx,ny,nz,x,iz,w0,w1,w2,w3,
     2   istop,itmax,iters,ierror,nlev,ilev,nlev_real,mgsolv,
     3   iok,iinfo,epsiln,errtol,omega,nu1,nu2,mgsmoo,
     4   ipc,rpc,pc,ac,cc,fc,tru)
c* *********************************************************************
c* purpose:
c*
c*    nested iteration for a linear multilevel method.
c*
c*    algorithm:  linear multigrid iteration (cs)
c*
c*    this routine is the full multigrid front-end for a multigrid 
c*    v-cycle solver.  in other words, at repeatedly calls the v-cycle
c*    multigrid solver on successively finer and finer grids.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
c*
c*    *** other declarations ***
      integer          ipc(*),iz(50,*),iok,ilev,iinfo,nlev,itmax
      integer          iters,ierror,level,itmxd,nlevd,iterd,iokd,istop
      integer          nx,ny,nz,nxf,nyf,nzf,nxc,nyc,nzc,nlev_real,istpd
      integer          nu1,nu2,mgsmoo,iinfod,mgsolv
      double precision epsiln,errd,errtol,omega
      double precision x(*),w0(*),w1(*),w2(*),w3(*)
      double precision rpc(*),pc(*),ac(*),cc(*),fc(*),tru(*)
c*
c*    *** recover gridsizes ***
      nxf = nx
      nyf = ny
      nzf = nz
      call mkcors(nlev-1,nxf,nyf,nzf,nxc,nyc,nzc)
c*
c*    *** move up grids: interpolate solution to finer, do v cycle ***
      if (iinfo.ne.0) then
c*       write(6,100)'% FMVCS: starting:  ',nxf,nyf,nzf,nxc,nyc,nzc
c* 100     format(a,2(2x,' [',i3,',',i3,',',i3,'] '))
      endif
      do 10 level = nlev_real, ilev+1, -1
c*
c*       *** call mv cycle ***
         errd   = 1.0e-5
         itmxd  = 1
         nlevd  = nlev_real - level + 1
         iterd  = 0
         iokd   = 2
         iinfod = iinfo
         istpd  = istop
         call mvcs(nxc,nyc,nzc,x,iz,w0,w1,w2,w3,
     2      istpd,itmxd,iterd,ierror,nlevd,level,nlev_real,mgsolv,
     3      iokd,iinfod,epsiln,errtol,omega,nu1,nu2,mgsmoo,
     4      ipc,rpc,pc,ac,cc,fc,tru)
c*
c*       *** find new grid size ***
         call mkfine(1,nxc,nyc,nzc,nxf,nyf,nzf)
c*
c*       *** interpolate to next finer grid ***
         call interpPMG(nxc,nyc,nzc,nxf,nyf,nzf,
     2      x(iz(1,level)),x(iz(1,level-1)),pc(iz(11,level-1)))
c*
c*       *** new grid size ***
         nxc = nxf
         nyc = nyf
         nzc = nzf
 10   continue
c*
c*    *** call mv cycle ***
      level = ilev
      call mvcs(nxf,nyf,nzf,x,iz,w0,w1,w2,w3,
     2   istop,itmax,iters,ierror,nlev,level,nlev_real,mgsolv,
     3   iok,iinfo,epsiln,errtol,omega,nu1,nu2,mgsmoo,
     4   ipc,rpc,pc,ac,cc,fc,tru)
c*
c*    *** return and end ***
      return
      end
      subroutine mvcs(nx,ny,nz,x,iz,w0,w1,w2,w3,
     2   istop,itmax,iters,ierror,nlev,ilev,nlev_real,mgsolv,
     3   iok,iinfo,epsiln,errtol,omega,nu1,nu2,mgsmoo,
     4   ipc,rpc,pc,ac,cc,fc,tru)
c* *********************************************************************
c* purpose:
c*
c*    screaming linear multilevel method.
c*
c*    algorithm:  linear multigrid iteration (cs)
c*
c*    multigrid v-cycle solver.
c*
c*    input:  
c*       (1) fine and coarse grid discrete linear operators: L_h, L_H
c*       (2) fine grid source function: f_h
c*       (3) fine grid approximate solution: u_h
c*
c*    output:
c*       (1) fine grid improved solution: u_h
c*
c*    the two-grid algorithm is:
c*       (1) pre-smooth:               u1_h = smooth(L_h,f_h,u_h)
c*       (2) restrict defect:          d_H  = r(L_h(u1_h) - f_h)
c*       (3) solve for correction:     c_H  = L_H^{-1}(d_H)
c*       (4) prolongate and correct:   u2_h = u1_h - p(c_H)
c*       (5) post-smooth:              u_h  = smooth(L_h,f_h,u2_h)
c*
c*    (of course, c_H is determined with another two-grid algorithm)
c*
c*    implementation notes:
c*       (0) "u1_h" must be kept on each level until "c_H" is computed,
c*           and then both are used to compute "u2_h".
c*       (1) "u_h" (and then "u1_h") on all levels is stored in the "x" array.
c*       (2) "d_H" is identically "f_h" for f_h on the next coarser grid.
c*       (3) "c_h" is identically "u_h" for u_h on the next coarser grid.
c*       (4) "d_H" is stored in the "r" array (must be kept for post-smooth).
c*       (5) "f_h" is stored in the "fc" array.
c*       (6) "L_h" on all levels is stored in the "ac" array.
c*       (7) signs may be reveresed; i.e., residual is used in place
c*           of the defect in places, etc.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
c*
c*    *** other declarations ***
      integer          ipc(*),iz(50,*),iok,ilev,iinfo,nlev,level,lev
      integer          itmax,iters,ierror,istop,nu1,nu2,mgsmoo
      integer          itmax_s,iters_s,nuuu,ivariv,mgsmoo_s,iresid
      integer          nx,ny,nz,nxf,nyf,nzf,nxc,nyc,nzc
      integer          lpv,n,m,lda,mgsolv,nlev_real,iadjoint
      double precision omega,errtol,epsiln,errtol_s
      double precision rsden,rsnrm,orsnrm,xnrm1,xnrm2,xdot,xnum,xden
      double precision xdamp
      double precision x(*),w0(*),w1(*),w2(*),w3(*)
      double precision rpc(*),pc(*),ac(*),cc(*),fc(*),tru(*)
c*
c*    *** recover level information ***
      level = 1
      lev   = (ilev-1)+level
c*
c*    *** recover gridsizes ***
      nxf = nx
      nyf = ny
      nzf = nz
      call mkcors(nlev-1,nxf,nyf,nzf,nxc,nyc,nzc)
c*
c*    *** do some i/o if requested ***
      if (iinfo.ne.0) then
c*         write(6,100)'% MVCS: starting:   ',nxf,nyf,nzf,nxc,nyc,nzc
c* 100     format(a,2(2x,' [',i3,',',i3,',',i3,'] '))
      endif
c*
c*    *** initial wall clock ***
      if (iok.ne.0) then
         call prtini(istop)
         call prtstp(iok,-1,0.0d0,0.0d0,0.0d0)
      endif
c*
c*    **************************************************************
c*    *** note: if (iok.ne.0) then:  use a stopping test.        ***
c*    ***       else:  use just the itmax to stop iteration.     ***
c*    **************************************************************
c*    *** istop=0 most efficient (whatever it is)                ***
c*    *** istop=1 relative residual                              ***
c*    *** istop=2 rms difference of successive iterates          ***
c*    *** istop=3 relative true error (provided for testing)     ***
c*    **************************************************************
c*
c*    *** compute denominator for stopping criterion ***
      if (iok.ne.0) then
         if (istop .eq. 0) then
            rsden = 1.0d0
         elseif (istop .eq. 1) then
            rsden = xnrm1(nxf,nyf,nzf,fc(iz(1,lev)))
         elseif (istop .eq. 2) then
            rsden = dsqrt(dble(nxf*nyf*nzf))
         elseif (istop .eq. 3) then
            rsden = xnrm2(nxf,nyf,nzf,tru(iz(1,lev)))
         elseif (istop .eq. 4) then
            rsden = xnrm2(nxf,nyf,nzf,tru(iz(1,lev)))
         elseif (istop .eq. 5) then
            call matvec(nxf,nyf,nzf,
     2         ipc(iz(5,lev)),rpc(iz(6,lev)),
     3         ac(iz(7,lev)),cc(iz(1,lev)),
     4         tru(iz(1,lev)),w1)
            rsden = dsqrt(xdot(nxf,nyf,nzf,tru(iz(1,lev)),w1))
         else
            call vnmprt(2,'% MVCS: bad istop value... ',27)
         endif
         if (rsden.eq.0.0d0) then
            rsden = 1.0d0
            call vnmprt(2,'% MVCS: rhs is zero on finest level ',36)
         endif
         rsnrm = rsden
         orsnrm = rsnrm
         call prtstp (iok,0,rsnrm,rsden,orsnrm)
      endif
c*
c* *********************************************************************
c* *** solve directly if nlev = 1 
c* *********************************************************************
c*
c*    *** solve directly if on the coarse grid ***
      if (nlev .eq. 1) then
c*
c*       *** use iterative method? ***
         if (mgsolv .eq. 0) then
c*
c*          *** solve on coarsest grid with cghs, mgsmoo_s=4 (no residual) ***
            iresid = 0
            iadjoint = 0
            itmax_s  = 100
            iters_s  = 0
            errtol_s = epsiln
            mgsmoo_s = 4
            call azeros(nxf,nyf,nzf,x(iz(1,lev)))
            call smooth (nxf,nyf,nzf,
     2         ipc(iz(5,lev)),rpc(iz(6,lev)),
     3         ac(iz(7,lev)),cc(iz(1,lev)),fc(iz(1,lev)),
     4         x(iz(1,lev)),w1,w2,w3,
     5         itmax_s,iters_s,errtol_s,omega,
     6         iresid,iadjoint,mgsmoo_s)
c*    
c*          *** check for trouble on the coarse grid ***
            if (iters_s .ge. itmax_s) then
               call vnmpri(2,'% MVCS: > max iters on coarse grid: ',
     2        36, iters_s)
            endif
c*
c*       *** use direct method? ***
         elseif (mgsolv .eq. 1) then
c*
c*          *** setup lpv to access the factored/banded operator ***
            lpv = lev+1
c*
c*          *** setup for banded format ***
            n   = ipc((iz(5,lpv)-1)+1)
            m   = ipc((iz(5,lpv)-1)+2)
            lda = ipc((iz(5,lpv)-1)+3)
c*
c*          *** call dpbsl to solve ***
            call xcopy_small(nxf,nyf,nzf,fc(iz(1,lev)),w1)
            call dpbsl(ac(iz(7,lpv)),lda,n,m,w1)
            call xcopy_large(nxf,nyf,nzf,w1,x(iz(1,lev)))
            call fboundPMG00(nxf,nyf,nzf,x(iz(1,lev)))
         else
            call vnmprt(2,'% MVCS: invalid coarse solver requested...',
     2        42)
         endif
c*      
c*       *** compute the stopping test ***
         iters = 1 
         if (iok.ne.0) then
            orsnrm = rsnrm
            if (istop .eq. 0) then
               call mresid(nxf,nyf,nzf,
     2            ipc(iz(5,lev)),rpc(iz(6,lev)),
     3            ac(iz(7,lev)),cc(iz(1,lev)),fc(iz(1,lev)),
     4            x(iz(1,lev)),w1)
               rsnrm = xnrm1(nxf,nyf,nzf,w1)
            elseif (istop .eq. 1) then
               call mresid(nxf,nyf,nzf,
     2            ipc(iz(5,lev)),rpc(iz(6,lev)),
     3            ac(iz(7,lev)),cc(iz(1,lev)),fc(iz(1,lev)),
     4            x(iz(1,lev)),w1)
               rsnrm = xnrm1(nxf,nyf,nzf,w1)
            elseif (istop .eq. 2) then
               call xcopy(nxf,nyf,nzf,tru(iz(1,lev)),w1)
               call xaxpy(nxf,nyf,nzf,(-1.0d0),x(iz(1,lev)),w1)
               rsnrm = xnrm1(nxf,nyf,nzf,w1)
               call xcopy(nxf,nyf,nzf,x(iz(1,lev)),tru(iz(1,lev)))
            elseif (istop .eq. 3) then
               call xcopy(nxf,nyf,nzf,tru(iz(1,lev)),w1)
               call xaxpy(nxf,nyf,nzf,(-1.0d0),x(iz(1,lev)),w1)
               rsnrm = xnrm2(nxf,nyf,nzf,w1)
            elseif (istop .eq. 4) then
               call xcopy(nxf,nyf,nzf,tru(iz(1,lev)),w1)
               call xaxpy(nxf,nyf,nzf,(-1.0d0),x(iz(1,lev)),w1)
               rsnrm = xnrm2(nxf,nyf,nzf,w1)
            elseif (istop .eq. 5) then
               call xcopy(nxf,nyf,nzf,tru(iz(1,lev)),w1)
               call xaxpy(nxf,nyf,nzf,(-1.0d0),x(iz(1,lev)),w1)
               call matvec(nxf,nyf,nzf,
     2            ipc(iz(5,lev)),rpc(iz(6,lev)),
     3            ac(iz(7,lev)),cc(iz(1,lev)),
     4            w1,w2)
               rsnrm = dsqrt(xdot(nxf,nyf,nzf,w1,w2))
            else
               call vnmprt(2,'% MVCS: bad istop value... ',27)
            endif
            call prtstp (iok,iters,rsnrm,rsden,orsnrm)
         endif
c*
c*       *** return now ***
         goto 99
      endif
c*
c* *********************************************************************
c* *** begin mg iteration (note nxf,nyf,nzf changes during loop)
c* *********************************************************************
c*
c*    *** setup for the v-cycle looping ***
      iters = 0 
 30   continue
c*
c*       *** finest level initialization ***
         level = 1
         lev   = (ilev-1)+level
c*
c*       *** nu1 pre-smoothings on fine grid (with residual) ***
         iresid = 1
         iadjoint = 0
         iters_s  = 0
         errtol_s = 0.0d0
         nuuu = ivariv (nu1,lev)
         call smooth(nxf,nyf,nzf,
     2      ipc(iz(5,lev)),rpc(iz(6,lev)),
     3      ac(iz(7,lev)),cc(iz(1,lev)),fc(iz(1,lev)),
     4      x(iz(1,lev)),w2,w3,w1,
     5      nuuu,iters_s,errtol_s,omega,
     6      iresid,iadjoint,mgsmoo)
         call xcopy(nxf,nyf,nzf,w1,w0(iz(1,lev)))
c*
c* *********************************************************************
c* begin cycling down to coarse grid
c* *********************************************************************
c*
c*       *** go down grids: restrict resid to coarser and smooth ***
         do 40 level = 2, nlev
            lev = (ilev-1)+level
c*
c*          *** find new grid size ***
            call mkcors(1,nxf,nyf,nzf,nxc,nyc,nzc)
c*
c*          *** restrict residual to coarser grid ***
            call restrc(nxf,nyf,nzf,nxc,nyc,nzc,
     2         w1,w0(iz(1,lev)),pc(iz(11,lev-1)))
c*
c*          *** new grid size ***
            nxf = nxc
            nyf = nyc
            nzf = nzc
c*
c*          *** if not on coarsest level yet... ***
            if (level .ne. nlev) then
c*
c*             *** nu1 pre-smoothings on this level (with residual) ***
c*             *** (w1 has residual...) ***
               call azeros(nxf,nyf,nzf,x(iz(1,lev)))
               iresid = 1
               iadjoint = 0
               iters_s  = 0
               errtol_s = 0.0d0
               nuuu = ivariv (nu1,lev)
               call smooth(nxf,nyf,nzf,
     2            ipc(iz(5,lev)),rpc(iz(6,lev)),
     3            ac(iz(7,lev)),cc(iz(1,lev)),w0(iz(1,lev)),
     4            x(iz(1,lev)),w2,w3,w1,
     5            nuuu,iters_s,errtol_s,omega,
     6            iresid,iadjoint,mgsmoo)
            endif
c*
c*       *** end of cycling down to coarse grid loop ***
 40      continue
c*
c* *********************************************************************
c* begin coarse grid
c* *********************************************************************
c*
c*       *** coarsest level ***
         level = nlev
         lev = (ilev-1)+level
c*
c*       *** use iterative method? ***
         if (mgsolv .eq. 0) then
c*
c*          *** solve on coarsest grid with cghs, mgsmoo_s=4 (no residual) ***
            iresid = 0
            iadjoint = 0
            itmax_s  = 100
            iters_s  = 0
            errtol_s = epsiln
            mgsmoo_s = 4
            call azeros(nxf,nyf,nzf,x(iz(1,lev)))
            call smooth (nxf,nyf,nzf,
     2         ipc(iz(5,lev)),rpc(iz(6,lev)),
     3         ac(iz(7,lev)),cc(iz(1,lev)),w0(iz(1,lev)),
     4         x(iz(1,lev)),w1,w2,w3,
     5         itmax_s,iters_s,errtol_s,omega,
     6         iresid,iadjoint,mgsmoo_s)
c*    
c*          *** check for trouble on the coarse grid ***
            if (iters_s .ge. itmax_s) then
               call vnmpri(2,'% MVCS: iters on coarse grid: ',
     2            30, iters_s)
            endif
c*
c*       *** use direct method? ***
         elseif (mgsolv .eq. 1) then
c*
c*          *** setup lpv to access the factored/banded operator ***
            lpv = lev+1
c*
c*          *** setup for banded format ***
            n   = ipc((iz(5,lpv)-1)+1)
            m   = ipc((iz(5,lpv)-1)+2)
            lda = ipc((iz(5,lpv)-1)+3)
c*
c*          *** call dpbsl to solve ***
            call xcopy_small(nxf,nyf,nzf,w0(iz(1,lev)),w1)
            call dpbsl(ac(iz(7,lpv)),lda,n,m,w1)
            call xcopy_large(nxf,nyf,nzf,w1,x(iz(1,lev)))
            call fboundPMG00(nxf,nyf,nzf,x(iz(1,lev)))
         else
            call vnmprt(2,'% MVCS: invalid coarse solver requested...',
     2         42)
         endif
c*      
c* *********************************************************************
c* begin cycling back to fine grid
c* *********************************************************************
c*
c*       *** move up grids: interpolate resid to finer and smooth ***
         do 70 level = nlev-1, 1, -1
            lev   = (ilev-1)+level
c*
c*          *** find new grid size ***
            call mkfine(1,nxf,nyf,nzf,nxc,nyc,nzc)
c*
c*          *** interpolate to next finer grid ***
            call interpPMG(nxf,nyf,nzf,nxc,nyc,nzc,
     2         x(iz(1,lev+1)),w1,pc(iz(11,lev)))
c*
c*          *** compute the hackbusch/reusken damping parameter ***
c*          *** which is equivalent to the standard linear cg steplength ***
            call matvec(nxf,nyf,nzf,
     2         ipc(iz(5,lev+1)),rpc(iz(6,lev+1)),
     3         ac(iz(7,lev+1)),cc(iz(1,lev+1)),
     4         x(iz(1,lev+1)),w2)
            xnum = xdot(nxf,nyf,nzf,x(iz(1,lev+1)),w0(iz(1,lev+1)))
            xden = xdot(nxf,nyf,nzf,x(iz(1,lev+1)),w2)
            xdamp = xnum / xden
c*
c*          *** new grid size ***
            nxf = nxc
            nyf = nyc
            nzf = nzc
c*
c*          *** perform the coarse grid correction ***
CZZZZZ      xdamp = 1.0d0
            call xaxpy(nxf,nyf,nzf,xdamp,w1,x(iz(1,lev)))
c*
c*          *** nu2 post-smoothings for correction (no residual) ***
            iresid = 0
            iadjoint = 1
            iters_s  = 0
            errtol_s = 0.0d0
            nuuu = ivariv (nu2,lev)
            if (level .eq. 1) then
               call smooth(nxf,nyf,nzf,
     2            ipc(iz(5,lev)),rpc(iz(6,lev)),
     3            ac(iz(7,lev)),cc(iz(1,lev)),fc(iz(1,lev)),
     4            x(iz(1,lev)),w1,w2,w3,
     5            nuuu,iters_s,errtol_s,omega,
     6            iresid,iadjoint,mgsmoo)
            else
               call smooth(nxf,nyf,nzf,
     2            ipc(iz(5,lev)),rpc(iz(6,lev)),
     3            ac(iz(7,lev)),cc(iz(1,lev)),w0(iz(1,lev)),
     4            x(iz(1,lev)),w1,w2,w3,
     5            nuuu,iters_s,errtol_s,omega,
     6            iresid,iadjoint,mgsmoo)
            endif
 70      continue
c*
c* *********************************************************************
c* iteration complete: do some i/o
c* *********************************************************************
c*
c*       *** increment the iteration counter ***
         iters = iters + 1
c*
c*       *** compute/check the current stopping test ***
         if (iok.ne.0) then
            orsnrm = rsnrm
            if (istop .eq. 0) then
               call mresid(nxf,nyf,nzf,
     2            ipc(iz(5,lev)),rpc(iz(6,lev)),
     3            ac(iz(7,lev)),cc(iz(1,lev)),fc(iz(1,lev)),
     4            x(iz(1,lev)),w1)
               rsnrm = xnrm1(nxf,nyf,nzf,w1)
            elseif (istop .eq. 1) then
               call mresid(nxf,nyf,nzf,
     2            ipc(iz(5,lev)),rpc(iz(6,lev)),
     3            ac(iz(7,lev)),cc(iz(1,lev)),fc(iz(1,lev)),
     4            x(iz(1,lev)),w1)
               rsnrm = xnrm1(nxf,nyf,nzf,w1)
            elseif (istop .eq. 2) then
               call xcopy(nxf,nyf,nzf,tru(iz(1,lev)),w1)
               call xaxpy(nxf,nyf,nzf,(-1.0d0),x(iz(1,lev)),w1)
               rsnrm = xnrm1(nxf,nyf,nzf,w1)
               call xcopy(nxf,nyf,nzf,x(iz(1,lev)),tru(iz(1,lev)))
            elseif (istop .eq. 3) then
               call xcopy(nxf,nyf,nzf,tru(iz(1,lev)),w1)
               call xaxpy(nxf,nyf,nzf,(-1.0d0),x(iz(1,lev)),w1)
               rsnrm = xnrm2(nxf,nyf,nzf,w1)
            elseif (istop .eq. 4) then
               call xcopy(nxf,nyf,nzf,tru(iz(1,lev)),w1)
               call xaxpy(nxf,nyf,nzf,(-1.0d0),x(iz(1,lev)),w1)
               rsnrm = xnrm2(nxf,nyf,nzf,w1)
            elseif (istop .eq. 5) then
               call xcopy(nxf,nyf,nzf,tru(iz(1,lev)),w1)
               call xaxpy(nxf,nyf,nzf,(-1.0d0),x(iz(1,lev)),w1)
               call matvec(nxf,nyf,nzf,
     2            ipc(iz(5,lev)),rpc(iz(6,lev)),
     3            ac(iz(7,lev)),cc(iz(1,lev)),
     4            w1,w2)
               rsnrm = dsqrt(xdot(nxf,nyf,nzf,w1,w2))
            else
               call vnmprt(2,'% MVCS: bad istop value... ',27)
            endif
            call prtstp (iok,iters,rsnrm,rsden,orsnrm)
            if ((rsnrm/rsden) .le. errtol) goto 99
         endif
         if (iters .ge. itmax) goto 91
      goto 30
c*
c*    *** problems ***
 91   continue
      ierror = 1
c*
c*    *** return and end ***
 99   continue
      return
      end

