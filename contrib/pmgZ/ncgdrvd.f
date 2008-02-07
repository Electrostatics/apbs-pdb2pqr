c* ///////////////////////////////////////////////////////////////////////////
c* @file    ncgdrvd.f
c* @author  Michael Holst
c* @brief   Driver for the nonlinear CG methods.
c* @version $Id: ncgdrvd.f 868 2006-04-05 18:32:05Z todd_dolinsky $
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

      subroutine ncghsdriv(iparm,rparm,iwork,rwork,u,
     2   xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf)
c* *********************************************************************
c* purpose:
c*
c*    linear/nonlinear conjugate gradient driver (fletcher-reeves).
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
      integer          iretot,iintot
c*
c*    *** misc variables ***
      integer          nrwk,niwk,nx,ny,nz,nlev,ierror,n
      integer          k_iz,k_w0
      integer          k_ipc,k_rpc,k_ac,k_cc,k_fc
      integer          n_iz,n_ipc,n_rpc
c*
c*    *** decode some parameters ***
      nrwk    = iparm(1)
      niwk    = iparm(2)
      nx      = iparm(3)
      ny      = iparm(4)
      nz      = iparm(5)
      nlev    = iparm(6)
      n       = nx * ny * nz
      n_iz    = 10*(nlev+1)
      n_ipc   = 100*(nlev+1)
      n_rpc   = 100*(nlev+1)
c*
c*    *** compute required work array sizes ***
      iintot  = n_iz + n_ipc
      iretot  = n_rpc + (3*n) + (4*n)
c*
c*    *** some more checks on input ***
      if ((nrwk.lt.iretot) .or. (niwk.lt.iintot)) then
         call vnmpri(2,'% NCGHSDRIV: real    work space must be: ',
     2      41, iretot)
         call vnmpri(2,'% NCGHSDRIV: integer work space must be: ',
     2      41, iintot)
         ierror = -3
         iparm(51) = ierror 
         return
      endif
c*
c*    *** iwork offsets ***
      k_iz  = 1
      k_ipc = k_iz  + n_iz
c*
c*    *** rwork offsets ***
      k_rpc = 1
      k_cc  = k_rpc + n_rpc
      k_fc  = k_cc  + n
      k_w0  = k_fc  + n
      k_ac  = k_w0  + n
c* ***k_ac_after  = k_ac + 4*n
c*
c*    *** call the multigrid driver ***
      call ncghsdriv2(iparm,rparm,nx,ny,nz,u,iwork(k_iz),
     2   rwork(k_w0),
     4   iwork(k_ipc),rwork(k_rpc),
     5   rwork(k_ac),rwork(k_cc),rwork(k_fc),
     8   xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf)
c*
c*    *** return and end ***
      return
      end
      subroutine ncghsdriv2(iparm,rparm,nx,ny,nz,u,iz,
     2   w0,
     3   ipc,rpc,ac,cc,fc,
     4   xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf)
c* *********************************************************************
c* purpose:
c*
c*    linear/nonlinear conjugate gradient driver (fletcher-reeves).
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
c*
c*    *** input parameters ***
      integer          iparm(*),ipc(*),iz(*)
      double precision rparm(*),rpc(*)
      double precision u(*),w0(*),ac(*),cc(*),fc(*)
      double precision xf(*),yf(*),zf(*),gxcf(*),gycf(*),gzcf(*)
      double precision a1cf(*),a2cf(*),a3cf(*),ccf(*),fcf(*),tcf(*)
c*
c*    *** misc variables ***
      integer          mgkey,nlev,itmax,iok,iinfo,istop,ipkey,nu1,nu2
      integer          nx,ny,nz,ilev,ido,iters,ierror,ibound,mode
      integer          mgprol,mgcoar,mgsolv,mgdisc
      double precision epsiln,epsmac,errtol,omegal,omegan
      double precision pc_dumm
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
      mode   = iparm(16)
      mgprol = iparm(17)
      mgcoar = iparm(18)
      mgdisc = iparm(19)
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
      call vtstrt(30, 'NCGHSDRIV2: fine problem setup', 30)
c*
c*    *** build op and rhs on fine grid ***
      ido = 0
      call buildops (nx,ny,nz,nlev,ipkey,iinfo,ido,iz,
     2   mgprol,mgcoar,mgsolv,mgdisc,
     3   ipc,rpc,pc_dumm,ac,cc,fc,
     4   xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf)
c*
c*    *** stop timer ***
      call vtstrt(30, 'NCGHSDRIV2: fine problem setup', 30)
      tsetupc = 0.0d0
c*
c* ******************************************************************
c* *** this overwrites the rhs array provided by pde specification
c* ****** compute an algebraically produced rhs for the given tcf ***
      if ((istop .eq. 4) .or. (istop .eq. 5)) then
         if ((mode .eq. 1) .or. (mode .eq. 2)) then
            call nmatvec(nx,ny,nz,ipc,rpc,ac,cc,tcf,fc,w0)
         else
            call matvec(nx,ny,nz,ipc,rpc,ac,cc,tcf,fc)
         endif
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
      call vnmprt(2,' cg = [ ',8)
c*
c*    *** start timer ***
      call vtstrt(30, 'NCGHSDRIV2: solve', 17)
c*
c*    *** call specified multigrid method ***
      if ((mode .eq. 0) .or. (mode .eq. 2)) then
         call vnmprt(2,'% NCGHSDRIV2: linear mode...',28)
         iok  = 1
         ilev = 1
         call cghsgo(nx,ny,nz,u,w0,a1cf,a2cf,a3cf,ccf,fcf,
     2      istop,itmax,iters,ierror,
     3      iok,iinfo,epsiln,errtol,omegal,
     4      ipc,rpc,ac,cc,fc,tcf)
      endif
      if ((mode .eq. 1) .or. (mode .eq. 2)) then
         call vnmprt(2,'% NCGHSDRIV2: nonlinear mode...',31)
         iok  = 1
         ilev = 1
         call ncghsgo(nx,ny,nz,u,w0,a1cf,a2cf,a3cf,ccf,fcf,
     2      istop,itmax,iters,ierror,
     3      iok,iinfo,epsiln,errtol,omegan,
     4      ipc,rpc,ac,cc,fc,tcf)
      endif
c*
c*    *** stop timer ***
      call vtstop(30, 'NCGHSDRIV2: solve', 17)
c*
c*    *** MATLAB ***
      write(*,100) 'cg_sf',tsetupf,'cg_sc',tsetupc,
     2   'cg_st',(tsetupf+tsetupc),'cg_so',tsolve
 100  format(' ];',4(' ',a7,'=',1pe9.3,';'))
c*
c*    *** restore boundary conditions ***
      ibound = 1
      call fboundPMG(ibound,nx,ny,nz,u,gxcf,gycf,gzcf)
c*
c*    *** return and end ***
      return
      end
      subroutine ncghsgo(nx,ny,nz,x,r,p,ap,zk,zkp1,tmp,
     2   istop,itmax,iters,ierror,
     3   iok,iinfo,epsiln,errtol,omega,
     4   ipc,rpc,ac,cc,fc,tru)
c* *********************************************************************
c* purpose:
c*
c*    nonlinear conjugate gradients (fletcher-reeves).
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
c*
c*    *** other declarations ***
      integer          ipc(*),iok,iinfo
      integer          itmax,iters,ierror
      integer          istop
      integer          nx,ny,nz
      double precision omega,errtol,epsiln
      double precision rsden,rsnrm,orsnrm,xnrm1,xnrm2,xdot
      double precision x(*),r(*),p(*),ap(*),zk(*),zkp1(*),tmp(*)
      double precision rpc(*),ac(*),cc(*),fc(*),tru(*)
      double precision alpha,rhok2,beta,rhok1
c*
cmdir 0 0
c*
c*    *** do some i/o if requested ***
      if (iinfo.ne.0) then
         write(6,100)'% NCGHSGO: starting:',nx,ny,nz
 100     format(a,(2x,' [',i3,',',i3,',',i3,'] '))
      endif
c*
c*    *** initial wall clock ***
      call prtini(istop)
      call prtstp(iok,-1,0.0d0,0.0d0,0.0d0)
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
      if (istop .eq. 0) then
         rsden = 1.0d0
      elseif (istop .eq. 1) then
c*       *** compute initial residual with zero initial guess ***
c*       *** this is analogous to the linear case where one can ***
c*       *** simply take norm of rhs for a zero initial guess ***
         call azeros(nx,ny,nz,tmp)
         call nmresid(nx,ny,nz,ipc,rpc,ac,cc,fc,tmp,r,zk)
         rsden = xnrm1(nx,ny,nz,r)
      elseif (istop .eq. 2) then
         rsden = dsqrt(dble((nx-2)*(ny-2)*(nz-2)))
      elseif (istop .eq. 3) then
         rsden = xnrm2(nx,ny,nz,tru)
      elseif (istop .eq. 4) then
         rsden = xnrm2(nx,ny,nz,tru)
      elseif (istop .eq. 5) then
         call nmatvec(nx,ny,nz,ipc,rpc,ac,cc,tru,r,zk)
         rsden = dsqrt(xdot(nx,ny,nz,tru,r))
      else
         call vnmprt(2,'% NCGHSGO: bad istop value... ',30)
      endif
      if (rsden.eq.0.0d0) then
         rsden = 1.0d0
         call vnmprt(2,'% NCGHSGO: rhs is zero ',23)
      endif
      rsnrm = rsden
      orsnrm = rsnrm
      call prtstp (iok,0,rsnrm,rsden,orsnrm)
c*
c*    *** now compute residual with the initial guess ***
      call nmresid(nx,ny,nz,ipc,rpc,ac,cc,fc,x,r,zk)
c*
c*    *** setup for the looping ***
      iters = 0 
 30   continue
c*
c*       *** save iterate if stop test will use it on next iter ***
         if (istop .eq. 2) call xcopy(nx,ny,nz,x,tru)
c*
c*       *** form new direction vector from old one and residual ***
         rhok2 = xdot(nx,ny,nz,r,r)
         if (iters .eq. 0) then
            call xcopy(nx,ny,nz,r,p)
         else
            beta = rhok2 / rhok1
            call xaxpy(nx,ny,nz,((1.0d0)/beta),r,p)
            call xscal(nx,ny,nz,beta,p)
         endif
c*
c*       *** nonlinear case: do a line search ***
c*       *** (note: "ap,zk,zkp1" passed back from line search as desired) ***
         call xcopy(nx,ny,nz,r,tmp)
         call linesearch(nx,ny,nz,alpha,
     2      ipc,rpc,ac,cc,fc,p,x,tmp,ap,zk,zkp1)
c*
c*       *** save rhok2 for next iteration ***
         rhok1 = rhok2
c*
c*       *** update solution in direction p of length alpha ***
         call xaxpy(nx,ny,nz,alpha,p,x)
c*
c*       *** update residual ***
         call xaxpy(nx,ny,nz,(-alpha),ap,r)
         call xaxpy(nx,ny,nz,(1.0d0),zk,r)
         call xaxpy(nx,ny,nz,(-1.0d0),zkp1,r)
c*
c* ***** *** switch to descent if necessary ***
c* ***** itmax_d = 10
c* ***** iter_d = 0
c* ***** alpha = -1.0d0
c* ***** rsnrm_tmp = xnrm1(nx,ny,nz,r)
c* ***** 18 continue
c* *****    if ((rsnrm_tmp.lt.rsnrm).or.(iter_d.gt.itmax_d)) then
c* *****       print*,'% finished with descent: r_o,r_n',rsnrm,rsnrm_tmp
c* *****       if (iter_d .gt. 0) call xcopy(nx,ny,nz,tmp4,r)
c* *****       goto 19
c* *****    endif
c* *****    print*,'% trying a descent:      r_o,r_n',rsnrm,rsnrm_tmp
c* *****    call xcopy(nx,ny,nz,tmp2,tmp)
c* *****    call xaxpy(nx,ny,nz,alpha,tmp3,tmp)
c* *****    call nmresid(nx,ny,nz,ipc,rpc,ac,cc,fc,tmp,tmp4,ap)
c* *****    rsnrm_tmp = xnrm1(nx,ny,nz,tmp4)
c* *****    alpha = alpha / 2.0d0
c* *****    iter_d = iter_d + 1
c* ***** goto 18
c* ***** 19 continue
c*
c*       *** some bookkeeping ***
         iters = iters + 1
c*
c*       *** compute/check the current stopping test ***
         orsnrm = rsnrm
         if (istop .eq. 0) then
            rsnrm = xnrm1(nx,ny,nz,r)
         elseif (istop .eq. 1) then
            rsnrm = xnrm1(nx,ny,nz,r)
         elseif (istop .eq. 2) then
            call xcopy(nx,ny,nz,tru,tmp)
            call xaxpy(nx,ny,nz,(-1.0d0),x,tmp)
            rsnrm = xnrm1(nx,ny,nz,tmp)
         elseif (istop .eq. 3) then
            call xcopy(nx,ny,nz,tru,tmp)
            call xaxpy(nx,ny,nz,(-1.0d0),x,tmp)
            rsnrm = xnrm2(nx,ny,nz,tmp)
         elseif (istop .eq. 4) then
            call xcopy(nx,ny,nz,tru,tmp)
            call xaxpy(nx,ny,nz,(-1.0d0),x,tmp)
            rsnrm = xnrm2(nx,ny,nz,tmp)
         elseif (istop .eq. 5) then
            call xcopy(nx,ny,nz,tru,tmp)
            call xaxpy(nx,ny,nz,(-1.0d0),x,tmp)
            call nmatvec(nx,ny,nz,ipc,rpc,ac,cc,tmp,zk,zkp1)
            rsnrm = dsqrt(xdot(nx,ny,nz,tmp,zk))
         else
            call vnmprt(2,'% NCGHSGO: bad istop value... ',30)
         endif
         call prtstp (iok,iters,rsnrm,rsden,orsnrm)
         if ((rsnrm/rsden) .le. errtol) goto 99
         if (iters .ge. itmax) goto 99
      goto 30
c*
c*    *** return and end ***
 99   continue
      return
      end
      subroutine cghsgo(nx,ny,nz,x,r,p,ap,zk,zkp1,tmp,
     2   istop,itmax,iters,ierror,
     3   iok,iinfo,epsiln,errtol,omega,
     4   ipc,rpc,ac,cc,fc,tru)
c* *********************************************************************
c* purpose:
c*
c*    linear conjugate gradients (hestenes-steifel).
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
c*
c*    *** other declarations ***
      integer          ipc(*),iok,iinfo
      integer          itmax,iters,ierror
      integer          istop
      integer          nx,ny,nz
      double precision omega,errtol,epsiln
      double precision rsden,rsnrm,orsnrm,xnrm1,xnrm2,xdot
      double precision x(*),r(*),p(*),ap(*),zk(*),zkp1(*),tmp(*)
      double precision rpc(*),ac(*),cc(*),fc(*),tru(*)
      double precision alpha,rhok2,beta,rhok1,pAp
c*
cmdir 0 0
c*
c*    *** do some i/o if requested ***
      if (iinfo.ne.0) then
         write(6,100)'% CGHSGO: starting: ',nx,ny,nz
 100     format(a,(2x,' [',i3,',',i3,',',i3,'] '))
      endif
c*
c*    *** initial wall clock ***
      call prtini(istop)
      call prtstp(iok,-1,0.0d0,0.0d0,0.0d0)
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
      if (istop .eq. 0) then
         rsden = 1.0d0
      elseif (istop .eq. 1) then
         rsden = xnrm1(nx,ny,nz,fc)
      elseif (istop .eq. 2) then
         rsden = dsqrt(dble((nx-2)*(ny-2)*(nz-2)))
      elseif (istop .eq. 3) then
         rsden = xnrm2(nx,ny,nz,tru)
      elseif (istop .eq. 4) then
         rsden = xnrm2(nx,ny,nz,tru)
      elseif (istop .eq. 5) then
         call matvec(nx,ny,nz,ipc,rpc,ac,cc,tru,r)
         rsden = dsqrt(xdot(nx,ny,nz,tru,r))
      else
         call vnmprt(2,'% CGHSGO: bad istop value... ',29)
      endif
      if (rsden.eq.0.0d0) then
         rsden = 1.0d0
         call vnmprt(2,'% CGHSGO: rhs is zero ',22)
      endif
      rsnrm = rsden
      orsnrm = rsnrm
      call prtstp (iok,0,rsnrm,rsden,orsnrm)
c*
c*    *** now compute residual with the initial guess ***
      call mresid(nx,ny,nz,ipc,rpc,ac,cc,fc,x,r)
c*
c*    *** setup for the looping ***
      iters = 0 
 30   continue
c*
c*       *** save iterate if stop test will use it on next iter ***
         if (istop .eq. 2) call xcopy(nx,ny,nz,x,tru)
c*
c*       *** form new direction vector from old one and residual ***
         rhok2 = xdot(nx,ny,nz,r,r)
         if (iters .eq. 0) then
            call xcopy(nx,ny,nz,r,p)
         else
            beta = rhok2 / rhok1
            call xaxpy(nx,ny,nz,((1.0d0)/beta),r,p)
            call xscal(nx,ny,nz,beta,p)
         endif
c*
c*       *** linear case: alpha which minimizes energy norm of error ***
         call matvec(nx,ny,nz,ipc,rpc,ac,cc,p,ap)
         pAp = xdot(nx,ny,nz,p,ap)
         alpha = rhok2 / pAp
c*
c*       *** save rhok2 for next iteration ***
         rhok1 = rhok2
c*
c*       *** update solution in direction p of length alpha ***
         call xaxpy(nx,ny,nz,alpha,p,x)
c*
c*       *** update residual ***
         call xaxpy(nx,ny,nz,(-alpha),ap,r)
c*
c*       *** some bookkeeping ***
         iters = iters + 1
c*
c*       *** compute/check the current stopping test ***
         orsnrm = rsnrm
         if (istop .eq. 0) then
            rsnrm = xnrm1(nx,ny,nz,r)
         elseif (istop .eq. 1) then
            rsnrm = xnrm1(nx,ny,nz,r)
         elseif (istop .eq. 2) then
            call xcopy(nx,ny,nz,tru,tmp)
            call xaxpy(nx,ny,nz,(-1.0d0),x,tmp)
            rsnrm = xnrm1(nx,ny,nz,tmp)
         elseif (istop .eq. 3) then
            call xcopy(nx,ny,nz,tru,tmp)
            call xaxpy(nx,ny,nz,(-1.0d0),x,tmp)
            rsnrm = xnrm2(nx,ny,nz,tmp)
         elseif (istop .eq. 4) then
            call xcopy(nx,ny,nz,tru,tmp)
            call xaxpy(nx,ny,nz,(-1.0d0),x,tmp)
            rsnrm = xnrm2(nx,ny,nz,tmp)
         elseif (istop .eq. 5) then
            call xcopy(nx,ny,nz,tru,tmp)
            call xaxpy(nx,ny,nz,(-1.0d0),x,tmp)
            call matvec(nx,ny,nz,ipc,rpc,ac,cc,tmp,zk)
            rsnrm = dsqrt(xdot(nx,ny,nz,tmp,zk))
         else
            call vnmprt(2,'% CGHSGO: bad istop value... ',29)
         endif
         call prtstp (iok,iters,rsnrm,rsden,orsnrm)
         if ((rsnrm/rsden) .le. errtol) goto 99
         if (iters .ge. itmax) goto 99
      goto 30
c*
c*    *** return and end ***
 99   continue
      return
      end
