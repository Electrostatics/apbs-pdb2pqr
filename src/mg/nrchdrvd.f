c* ///////////////////////////////////////////////////////////////////////////
c* MG/XMG -- Multilevel nonlinear scalar elliptic PDE solver and X interface
c* Copyright (C) 1995  Michael Holst
c*
c* This program is free software; you can redistribute it and/or modify
c* it under the terms of the GNU General Public License as published by
c* the Free Software Foundation; either version 2 of the License, or
c* (at your option) any later version.
c*
c* This program is distributed in the hope that it will be useful,
c* but WITHOUT ANY WARRANTY; without even the implied warranty of
c* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c* GNU General Public License for more details.
c*
c* You should have received a copy of the GNU General Public License
c* along with this program; if not, write to the Free Software
c* Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
c*
c* MG/XMG was developed by:
c*
c*    Michael Holst                TELE:  (858) 534-4899
c*    Department of Mathematics    FAX:   (858) 534-5273
c*    UC San Diego, AP&M 5739      EMAIL: mholst@math.ucsd.edu
c*    La Jolla, CA 92093 USA       WEB:   http://www.scicomp.ucsd.edu/~mholst
c*
c* See the file "maind.f" for more information and pointers to papers.
c*
c* IMPORTANT: If you intend to use or modify this code, make sure you 
c* understand your responsibilities under the GNU license.
c* ///////////////////////////////////////////////////////////////////////////

      subroutine nrichdriv(iparm,rparm,iwork,rwork,u,
     2   xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf)
c* *********************************************************************
c* purpose:
c*
c*    driver for linear/nonlinear richardson's iteration.
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
         print*,'% NRICHDRIV: real    work space must be: ',iretot
         print*,'% NRICHDRIV: integer work space must be: ',iintot
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
      call nrichdriv2(iparm,rparm,nx,ny,nz,u,iwork(k_iz),
     2   rwork(k_w0),
     3   iwork(k_ipc),rwork(k_rpc),
     4   rwork(k_ac),rwork(k_cc),rwork(k_fc),
     5   xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf)
c*
c*    *** return and end ***
      return
      end
      subroutine nrichdriv2(iparm,rparm,nx,ny,nz,u,iz,
     2   w0,
     3   ipc,rpc,ac,cc,fc,
     4   xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf)
c* *********************************************************************
c* purpose:
c*
c*    driver for linear/nonlinear richardson's iteration.
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
      call vtstrt(30, 'NRICHDRIV2: fine problem setup', 30)
c*
c*    *** build op and rhs on fine grid ***
      ido = 0
      call buildops (nx,ny,nz,nlev,ipkey,iinfo,ido,iz,
     2   mgprol,mgcoar,mgsolv,mgdisc,
     3   ipc,rpc,pc_dumm,ac,cc,fc,
     4   xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf)
c*
c*    *** stop timer ***
      call vtstop(30, 'NRICHDRIV2: fine problem setup', 30)
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
      print*,' rich = [ '
c*
c*    *** start timer ***
      call vtstrt(30, 'NRICHDRIV2: solve', 17)
c*
c*    *** call specified multigrid method ***
      if ((mode .eq. 0) .or. (mode .eq. 2)) then
         print*,'% NRICHDRIV2: linear mode...'
         iok  = 1
         ilev = 1
         call richgo(nx,ny,nz,u,w0,a1cf,a2cf,
     2      istop,itmax,iters,ierror,
     3      iok,iinfo,epsiln,errtol,omegal,
     4      ipc,rpc,ac,cc,fc,tcf)
      endif
      if ((mode .eq. 1) .or. (mode .eq. 2)) then
         print*,'% NRICHDRIV2: nonlinear mode...'
         iok  = 1
         ilev = 1
         call nrichgo(nx,ny,nz,u,w0,a1cf,a2cf,
     2      istop,itmax,iters,ierror,
     3      iok,iinfo,epsiln,errtol,omegan,
     4      ipc,rpc,ac,cc,fc,tcf)
      endif
c*
c*    *** stop timer ***
      call vtstop(30, 'NRICHDRIV2: solve', 17)
c*
c*    *** MATLAB ***
      write(*,100) 'rich_sf',tsetupf,'rich_sc',tsetupc,
     2   'rich_st',(tsetupf+tsetupc),'rich_so',tsolve
 100  format(' ];',4(' ',a7,'=',1pe9.3,';'))
c*
c*    *** restore boundary conditions ***
      ibound = 1
      call fboundPMG(ibound,nx,ny,nz,u,gxcf,gycf,gzcf)
c*
c*    *** return and end ***
      return
      end
      subroutine nrichgo(nx,ny,nz,x,r,w1,w2,
     2   istop,itmax,iters,ierror,
     3   iok,iinfo,epsiln,errtol,omega,
     4   ipc,rpc,ac,cc,fc,tru)
c* *********************************************************************
c* purpose:
c*
c*    nonlinear richardson's iteration.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
c*
c*    *** other declarations ***
      integer          ipc(*),iok,iinfo
      integer          itmax,iters,ierror
      integer          iresid,iadjoint,istop,itmax_s,iters_s
      integer          nx,ny,nz
      double precision omega,errtol,epsiln,errtol_s,omega_s
      double precision rsden,rsnrm,orsnrm,xnrm1,xnrm2,xdot
      double precision x(*),r(*),w1(*),w2(*)
      double precision rpc(*),ac(*),cc(*),fc(*),tru(*)
c*
c*    *** do some i/o if requested ***
      if (iinfo.ne.0) then
         write(6,100)'% NRICHGO: starting:',nx,ny,nz
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
         call azeros(nx,ny,nz,w1)
         call nmresid(nx,ny,nz,ipc,rpc,ac,cc,fc,w1,r,w2)
         rsden = xnrm1(nx,ny,nz,r)
      elseif (istop .eq. 2) then
         rsden = dsqrt(dble(nx*ny*nz))
      elseif (istop .eq. 3) then
         rsden = xnrm2(nx,ny,nz,tru)
      elseif (istop .eq. 4) then
         rsden = xnrm2(nx,ny,nz,tru)
      elseif (istop .eq. 5) then
         call nmatvec(nx,ny,nz,ipc,rpc,ac,cc,tru,w1,w2)
         rsden = dsqrt(xdot(nx,ny,nz,tru,w1))
      else
         print*,'% NRICHGO: bad istop value... '
      endif
      if (rsden.eq.0.0d0) then
         rsden = 1.0d0
         print*,'% NRICHGO: rhs is zero '
      endif
      rsnrm = rsden
      orsnrm = rsnrm
      call prtstp (iok,0,rsnrm,rsden,orsnrm)
c*
c*    *** setup for the looping ***
      iters  = 0 
 30   continue
         iters  = iters  + 10
c*
c*       *** save iterate if stop test will use it on next iter ***
         if (istop .eq. 2) call xcopy(nx,ny,nz,x,tru)
c*
c*       *** do 10 iterations ***
         iresid = 1
         iadjoint = 0
         errtol_s = 0.0d0
         itmax_s  = 10
         omega_s  = omega
         call nrich(nx,ny,nz,ipc,rpc,ac,cc,fc,x,w1,w2,r,
     2      itmax_s,iters_s,errtol_s,omega_s,iresid,iadjoint)
c*
c*       *** compute/check the current stopping test ***
         orsnrm = rsnrm
         if (istop .eq. 0) then
            rsnrm = xnrm1(nx,ny,nz,r)
         elseif (istop .eq. 1) then
            rsnrm = xnrm1(nx,ny,nz,r)
         elseif (istop .eq. 2) then
            call xcopy(nx,ny,nz,tru,w1)
            call xaxpy(nx,ny,nz,(-1.0d0),x,w1)
            rsnrm = xnrm1(nx,ny,nz,w1)
         elseif (istop .eq. 3) then
            call xcopy(nx,ny,nz,tru,w1)
            call xaxpy(nx,ny,nz,(-1.0d0),x,w1)
            rsnrm = xnrm2(nx,ny,nz,w1)
         elseif (istop .eq. 4) then
            call xcopy(nx,ny,nz,tru,w1)
            call xaxpy(nx,ny,nz,(-1.0d0),x,w1)
            rsnrm = xnrm2(nx,ny,nz,w1)
         elseif (istop .eq. 5) then
            call xcopy(nx,ny,nz,tru,w1)
            call xaxpy(nx,ny,nz,(-1.0d0),x,w1)
            call nmatvec(nx,ny,nz,ipc,rpc,ac,cc,w1,w2,r)
            rsnrm = dsqrt(xdot(nx,ny,nz,w1,w2))
         else
            print*,'% NRICHGO: bad istop value... '
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
      subroutine richgo(nx,ny,nz,x,r,w1,w2,
     2   istop,itmax,iters,ierror,
     3   iok,iinfo,epsiln,errtol,omega,
     4   ipc,rpc,ac,cc,fc,tru)
c* *********************************************************************
c* purpose:
c*
c*    linear richardson's iteration.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
c*
c*    *** other declarations ***
      integer          ipc(*),iok,iinfo
      integer          itmax,iters,ierror
      integer          iresid,iadjoint,istop,itmax_s,iters_s
      integer          nx,ny,nz
      double precision omega,errtol,epsiln,errtol_s,omega_s
      double precision rsden,rsnrm,orsnrm,xnrm1,xnrm2,xdot
      double precision x(*),r(*),w1(*),w2(*)
      double precision rpc(*),ac(*),cc(*),fc(*),tru(*)
c*
c*    *** do some i/o if requested ***
      if (iinfo.ne.0) then
         write(6,100)'% RICHGO: starting: ',nx,ny,nz
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
         rsden = dsqrt(dble(nx*ny*nz))
      elseif (istop .eq. 3) then
         rsden = xnrm2(nx,ny,nz,tru)
      elseif (istop .eq. 4) then
         rsden = xnrm2(nx,ny,nz,tru)
      elseif (istop .eq. 5) then
         call matvec(nx,ny,nz,ipc,rpc,ac,cc,tru,w1)
         rsden = dsqrt(xdot(nx,ny,nz,tru,w1))
      else
         print*,'% RICHGO: bad istop value... '
      endif
      if (rsden.eq.0.0d0) then
         rsden = 1.0d0
         print*,'% RICHGO: rhs is zero '
      endif
      rsnrm = rsden
      orsnrm = rsnrm
      call prtstp (iok,0,rsnrm,rsden,orsnrm)
c*
c*    *** setup for the looping ***
      iters  = 0 
 30   continue
         iters  = iters  + 10
c*
c*       *** save iterate if stop test will use it on next iter ***
         if (istop .eq. 2) call xcopy(nx,ny,nz,x,tru)
c*
c*       *** do 10 iterations ***
         iresid = 1
         iadjoint = 0
         errtol_s = 0.0d0
         itmax_s  = 10
         omega_s  = omega
         call rich(nx,ny,nz,ipc,rpc,ac,cc,fc,x,w1,w2,r,
     2      itmax_s,iters_s,errtol_s,omega_s,iresid,iadjoint)
c*
c*       *** compute/check the current stopping test ***
         orsnrm = rsnrm
         if (istop .eq. 0) then
            rsnrm = xnrm1(nx,ny,nz,r)
         elseif (istop .eq. 1) then
            rsnrm = xnrm1(nx,ny,nz,r)
         elseif (istop .eq. 2) then
            call xcopy(nx,ny,nz,tru,w1)
            call xaxpy(nx,ny,nz,(-1.0d0),x,w1)
            rsnrm = xnrm1(nx,ny,nz,w1)
         elseif (istop .eq. 3) then
            call xcopy(nx,ny,nz,tru,w1)
            call xaxpy(nx,ny,nz,(-1.0d0),x,w1)
            rsnrm = xnrm2(nx,ny,nz,w1)
         elseif (istop .eq. 4) then
            call xcopy(nx,ny,nz,tru,w1)
            call xaxpy(nx,ny,nz,(-1.0d0),x,w1)
            rsnrm = xnrm2(nx,ny,nz,w1)
         elseif (istop .eq. 5) then
            call xcopy(nx,ny,nz,tru,w1)
            call xaxpy(nx,ny,nz,(-1.0d0),x,w1)
            call matvec(nx,ny,nz,ipc,rpc,ac,cc,w1,w2)
            rsnrm = dsqrt(xdot(nx,ny,nz,w1,w2))
         else
            print*,'% RICHGO: bad istop value... '
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
