c* ///////////////////////////////////////////////////////////////////////////
c* @file    powerd.f
c* @author  Michael Holst
c* @brief   Power and Inverse Power Methods for doing spectral diagnostics.
c* @version $Id: powerd.f 868 2006-04-05 18:32:05Z todd_dolinsky $
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

      subroutine power(nx,ny,nz,iz,ilev,ipc,rpc,ac,cc,
     2   w1,w2,w3,w4,eigmax,eigmax_model,tol,itmax,iters,iinfo)
c* *********************************************************************
c* purpose:
c*
c*    standard power method for maximum eigenvalue estimation
c*    of the system matrix A.
c*
c* notes:
c*
c*    to test, note that the 3d laplacean has min/max eigenvalues:
c*
c*       lambda_min = 6 - 2*dcos(pi/(nx-1))
c*                      - 2*dcos(pi/(ny-1))
c*                      - 2*dcos(pi/(nz-1))
c*
c*       lambda_max = 6 - 2*dcos((nx-2)*pi/(nx-1))
c*                      - 2*dcos((ny-2)*pi/(ny-1))
c*                      - 2*dcos((nz-2)*pi/(nz-1))
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),iz(50,*),ilev,lev,nx,ny,nz,itmax,iters
      integer          iinfo,level
      double precision rpc(*),ac(*),cc(*)
      double precision w1(*),w2(*),w3(*),w4(*)
      double precision eigmax,eigmax_model,denom,fac,tol,rho,oldrho
      double precision error,relerr,xnrm2,xdot
c*
c*    *** some parameters ***
      double precision pi
      pi               = 4.0d0 * datan(1.0d0)
c*
c*    *** recover level information ***
      level = 1
      lev   = (ilev-1)+level
c*
c*    *** seed vector: random to contain all components ***
      call axrand(nx,ny,nz,w1)
      call azeros(nx,ny,nz,w2)
      call azeros(nx,ny,nz,w3)
      call azeros(nx,ny,nz,w4)
c*
c*    *** compute raleigh quotient with the seed vector ***
      denom = xnrm2(nx,ny,nz,w1)
      fac = 1.0 / denom
      call xscal(nx,ny,nz,fac,w1)
      call matvec(nx,ny,nz,ipc(iz(5,lev)),rpc(iz(6,lev)),
     2   ac(iz(7,lev)),cc(iz(1,lev)),w1,w2)
      oldrho = xdot(nx,ny,nz,w1,w2)
c*
c*    *** i/o ***
      if (oldrho .eq. 0.0e0) then
         if (iinfo .gt. 3) then
            write(6,510) 'POWER: iter; estimate = ',iters,oldrho
 510        format('% ',a,i5,1pe11.4)
         endif
         rho = oldrho
         goto 99
      endif
c*
c*    *** main iteration ***
      iters = 0
 20   continue
         iters = iters + 1
c*
c*       *** apply the matrix A ***
         call matvec(nx,ny,nz,ipc(iz(5,lev)),rpc(iz(6,lev)),
     2      ac(iz(7,lev)),cc(iz(1,lev)),w1,w2)
         call xcopy(nx,ny,nz,w2,w1)
c*
c*       *** normalize the new vector ***
         denom = xnrm2(nx,ny,nz,w1)
         fac = 1.0 / denom
         call xscal(nx,ny,nz,fac,w1)
c*
c*       *** compute the new raleigh quotient ***
         call matvec(nx,ny,nz,ipc(iz(5,lev)),rpc(iz(6,lev)),
     2      ac(iz(7,lev)),cc(iz(1,lev)),w1,w2)
         rho = xdot(nx,ny,nz,w1,w2)
c*
c*       *** stopping test ***
c*       *** w2=A*x, w1=x, stop = 2-norm(A*x-lamda*x) ***
         call xcopy(nx,ny,nz,w1,w3)
         call xcopy(nx,ny,nz,w2,w4)
         call xscal(nx,ny,nz,rho,w3)
         call xaxpy(nx,ny,nz,(-1.0d0),w3,w4)
         error = xnrm2(nx,ny,nz,w4)
         relerr = dabs( rho - oldrho ) / dabs( rho )
c*
c*       *** i/o ***
         if (iinfo .gt. 3) then
            write(6,500)
     2         'POWER:  iter; error; relerr; estimate = ',
     3         iters,error,relerr,rho
 500           format('% ',a,i5,2(1pe9.2),1pe11.4)
         endif
         if ((relerr .lt. tol) .or. (iters .eq. itmax)) then
            goto 99
         endif
         oldrho = rho
      goto 20
c*
c*    *** return some stuff ***
 99   continue
      eigmax = rho
      fac = 2.0e0**(ilev-1)
      eigmax_model = fac * (6.0d0 
     2   - 2.0d0*dcos(dble(nx-2)*pi/dble(nx-1))
     2   - 2.0d0*dcos(dble(ny-2)*pi/dble(ny-1))
     2   - 2.0d0*dcos(dble(nz-2)*pi/dble(nz-1)) )
c*
c*    *** return and end ***
      return
      end
      subroutine ipower(nx,ny,nz,u,iz,w0,w1,w2,w3,w4,
     2      eigmin,eigmin_model,tol,itmax,iters,
     3      nlev,ilev,nlev_real,mgsolv,
     4      iok,iinfo,epsiln,errtol,omega,nu1,nu2,mgsmoo,
     5      ipc,rpc,pc,ac,cc,tru)
c* *********************************************************************
c* purpose:
c*
c*    standard inverse power method for minimum eigenvalue estimation
c*    of the system matrix A.
c*
c* notes:
c*
c*    to test, note that the 3d laplacean has min/max eigenvalues:
c*
c*       lambda_min = 6 - 2*dcos(pi/(nx-1))
c*                      - 2*dcos(pi/(ny-1))
c*                      - 2*dcos(pi/(nz-1))
c*
c*       lambda_max = 6 - 2*dcos((nx-2)*pi/(nx-1))
c*                      - 2*dcos((ny-2)*pi/(ny-1))
c*                      - 2*dcos((nz-2)*pi/(nz-1))
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
c*
c*    *** other declarations ***
      integer          ipc(*),iz(50,*),iok,ilev,iinfo,nx,ny,nz
      integer          nu1,nu2,mgsmoo,mgsolv
      integer          itmax,iters,level,lev,nlev,nlev_real
      double precision epsiln,errtol,omega
      double precision u(*),w0(*)
      double precision w1(*),w2(*),w3(*),w4(*)
      double precision rpc(*),pc(*),ac(*),cc(*),tru(*)
c*
c*    *** some more ***
      double precision eigmin,eigmin_model,denom,fac,tol,rho,oldrho
      double precision error,relerr,xnrm2,xdot,errtol_s
      integer          itmax_s,iters_s,ierror_s,iok_s,iinfo_s,istop_s
      integer          nu1_s,nu2_s,mgsmoo_s
c*
c*    *** some parameters ***
      double precision pi
      pi               = 4.0d0 * datan(1.0d0)
c*
c*    *** recover level information ***
      level = 1
      lev   = (ilev-1)+level
c*
c*    *** seed vector: random to contain all components ***
      call axrand(nx,ny,nz,w1)
      call azeros(nx,ny,nz,w2)
      call azeros(nx,ny,nz,w3)
      call azeros(nx,ny,nz,w4)
      call azeros(nx,ny,nz,w0(iz(1,lev)))
      call azeros(nx,ny,nz,u(iz(1,lev)))
c*
c*    *** compute raleigh quotient with the seed vector ***
      denom = xnrm2(nx,ny,nz,w1)
      fac = 1.0 / denom
      call xscal(nx,ny,nz,fac,w1)
      call matvec(nx,ny,nz,ipc(iz(5,lev)),rpc(iz(6,lev)),
     2   ac(iz(7,lev)),cc(iz(1,lev)),w1,w2)
      oldrho = xdot(nx,ny,nz,w1,w2)
c*
c*    *** i/o ***
      if (oldrho .eq. 0.0e0) then
         if (iinfo .gt. 3) then
            write(6,510) 'IPOWER: iter; estimate = ',iters,oldrho
 510        format('% ',a,i5,1pe11.4)
         endif
         rho = oldrho
         goto 99
      endif
c*
c*    *** main iteration ***
      iters = 0
 20   continue
         iters = iters + 1
c*
c*       *** apply the matrix A^{-1} (using MG solver) ***
         itmax_s = 100
         iters_s = 0
         ierror_s = 0
         iok_s = 0
         iinfo_s = 0
         istop_s = 0
         mgsmoo_s = 1
         nu1_s = 1
         nu2_s = 1
         errtol_s = epsiln
         call xcopy(nx,ny,nz,w1,w0(iz(1,lev)))
         call mvcs(nx,ny,nz,u,iz,w1,w2,w3,w4,
     2      istop_s,itmax_s,iters_s,ierror_s,
     3      nlev,ilev,nlev_real,mgsolv,
     4      iok_s,iinfo_s,epsiln,errtol_s,omega,nu1_s,nu2_s,mgsmoo_s,
     5      ipc,rpc,pc,ac,cc,w0,tru)
         call xcopy(nx,ny,nz,u(iz(1,lev)),w1)
c*
c*       *** normalize the new vector ***
         denom = xnrm2(nx,ny,nz,w1)
         fac = 1.0 / denom
         call xscal(nx,ny,nz,fac,w1)
c*
c*       *** compute the new raleigh quotient ***
         call matvec(nx,ny,nz,ipc(iz(5,lev)),rpc(iz(6,lev)),
     2      ac(iz(7,lev)),cc(iz(1,lev)),w1,w2)
         rho = xdot(nx,ny,nz,w1,w2)
c*
c*       *** stopping test ***
c*       *** w2=A*x, w1=x, stop = 2-norm(A*x-lamda*x) ***
         call xcopy(nx,ny,nz,w1,w3)
         call xcopy(nx,ny,nz,w2,w4)
         call xscal(nx,ny,nz,rho,w3)
         call xaxpy(nx,ny,nz,(-1.0d0),w3,w4)
         error = xnrm2(nx,ny,nz,w4)
         relerr = dabs( rho - oldrho ) / dabs( rho )
c*
c*       *** i/o ***
         if (iinfo .gt. 3) then
            write(6,500)
     2         'IPOWER: iter; error; relerr; estimate = ',
     3         iters,error,relerr,rho
 500           format('% ',a,i5,2(1pe9.2),1pe11.4)
         endif
         if ((relerr .lt. tol) .or. (iters .eq. itmax)) then
            goto 99
         endif
         oldrho = rho
      goto 20
c*
c*    *** return some stuff ***
 99   continue
      eigmin = rho
      fac = 2.0e0**(ilev-1)
      eigmin_model = fac * (6.0d0 
     2   - 2.0d0*dcos(pi/dble(nx-1))
     2   - 2.0d0*dcos(pi/dble(ny-1))
     2   - 2.0d0*dcos(pi/dble(nz-1)) )
c*
c*    *** return and end ***
      return
      end
      subroutine mpower(nx,ny,nz,u,iz,w0,w1,w2,w3,w4,
     2      eigmax,tol,itmax,iters,
     3      nlev,ilev,nlev_real,mgsolv,
     4      iok,iinfo,epsiln,errtol,omega,nu1,nu2,mgsmoo,
     5      ipc,rpc,pc,ac,cc,fc,tru)
c* *********************************************************************
c* purpose:
c*
c*    standard power method for maximimum eigenvalue estimation
c*    of the multigrid operator M.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
c*
c*    *** other declarations ***
      integer          ipc(*),iz(50,*),iok,ilev,iinfo,nx,ny,nz
      integer          nu1,nu2,mgsmoo,mgsolv
      integer          itmax,iters,level,lev,nlev,nlev_real
      double precision epsiln,errtol,omega
      double precision u(*),w0(*)
      double precision w1(*),w2(*),w3(*),w4(*)
      double precision rpc(*),pc(*),ac(*),cc(*),fc(*),tru(*)
c*
c*    *** some more ***
      double precision eigmax,denom,fac,tol,rho,oldrho,error
      double precision relerr,xnrm2,xdot
      integer          itmax_s,iters_s,ierror_s,iok_s,iinfo_s,istop_s
c*
c*    *** recover level information ***
      level = 1
      lev   = (ilev-1)+level
c*
c*    *** seed vector: random to contain all components ***
      call axrand(nx,ny,nz,w1)
      call azeros(nx,ny,nz,w2)
      call azeros(nx,ny,nz,w3)
      call azeros(nx,ny,nz,w4)
      call azeros(nx,ny,nz,u(iz(1,lev)))
c*
c*    *** NOTE: we destroy "fc" on this level due to lack of vectors... ***
      call azeros(nx,ny,nz,fc(iz(1,lev)))
c*
c*    *** normalize the seed vector ***
      denom = xnrm2(nx,ny,nz,w1)
      fac = 1.0 / denom
      call xscal(nx,ny,nz,fac,w1)
c*
c*    *** compute raleigh quotient with the seed vector ***
      call xcopy(nx,ny,nz,w1,u(iz(1,lev)))
      itmax_s = 1
      iters_s = 0
      ierror_s = 0
      iok_s = 0
      iinfo_s = 0
      istop_s = 1
      call mvcs(nx,ny,nz,u,iz,w0,w2,w3,w4,
     2   istop_s,itmax_s,iters_s,ierror_s,
     3   nlev,ilev,nlev_real,mgsolv,
     4   iok_s,iinfo_s,epsiln,errtol,omega,nu1,nu2,mgsmoo,
     5   ipc,rpc,pc,ac,cc,fc,tru)
      oldrho = xdot(nx,ny,nz,w1,u(iz(1,lev)))
c*
c*    *** i/o ***
      if (oldrho .eq. 0.0e0) then
         if (iinfo .gt. 3) then
            write(6,510) 'MPOWER: iter; estimate = ',iters,oldrho
 510        format('% ',a,i5,1pe11.4)
         endif
         rho = oldrho
         goto 99
      endif
c*
c*    *** main iteration ***
      iters = 0
 20   continue
         iters = iters + 1
c*
c*       *** apply the matrix M ***
         call xcopy(nx,ny,nz,w1,u(iz(1,lev)))
         itmax_s = 1
         iters_s = 0
         ierror_s = 0
         iok_s = 0
         iinfo_s = 0
         istop_s = 1
         call mvcs(nx,ny,nz,u,iz,w1,w2,w3,w4,
     2      istop_s,itmax_s,iters_s,ierror_s,
     3      nlev,ilev,nlev_real,mgsolv,
     4      iok_s,iinfo_s,epsiln,errtol,omega,nu1,nu2,mgsmoo,
     5      ipc,rpc,pc,ac,cc,fc,tru)
         call xcopy(nx,ny,nz,u(iz(1,lev)),w1)
c*
c*       *** normalize the new vector ***
         denom = xnrm2(nx,ny,nz,w1)
         fac = 1.0 / denom
         call xscal(nx,ny,nz,fac,w1)
c*
c*       *** compute the new raleigh quotient ***
         call xcopy(nx,ny,nz,w1,u(iz(1,lev)))
         itmax_s = 1
         iters_s = 0
         ierror_s = 0
         iok_s = 0
         iinfo_s = 0
         istop_s = 1
         call mvcs(nx,ny,nz,u,iz,w0,w2,w3,w4,
     2      istop_s,itmax_s,iters_s,ierror_s,
     3      nlev,ilev,nlev_real,mgsolv,
     4      iok_s,iinfo_s,epsiln,errtol,omega,nu1,nu2,mgsmoo,
     5      ipc,rpc,pc,ac,cc,fc,tru)
         call xcopy(nx,ny,nz,u(iz(1,lev)),w2)
         rho = xdot(nx,ny,nz,w1,w2)
c*
c*       *** stopping test ***
c*       *** w2=A*x, w1=x, stop = 2-norm(A*x-lamda*x) ***
         call xcopy(nx,ny,nz,w1,w3)
         call xcopy(nx,ny,nz,w2,w4)
         call xscal(nx,ny,nz,rho,w3)
         call xaxpy(nx,ny,nz,(-1.0d0),w3,w4)
         error = xnrm2(nx,ny,nz,w4)
         relerr = dabs( rho - oldrho ) / dabs( rho )
c*
c*       *** i/o ***
         if (iinfo .gt. 3) then
            write(6,500)
     2         'MPOWER: iter; error; relerr; estimate = ',
     3         iters,error,relerr,rho
 500           format('% ',a,i5,2(1pe9.2),1pe11.4)
         endif
         if ((relerr .lt. tol) .or. (iters .eq. itmax)) then
            goto 99
         endif
         oldrho = rho
      goto 20
c*
c*    *** return some stuff ***
 99   continue
      eigmax = rho
c*
c*    *** return and end ***
      return
      end

