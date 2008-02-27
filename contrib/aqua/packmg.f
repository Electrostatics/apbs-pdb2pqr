      subroutine packmg (iparm,rparm,
     2   nrwk,niwk,nx,ny,nz,nlev,nu1,nu2,mgkey,itmax,istop,ipcon,
     3   nonlin,mgsmoo,mgprol,mgcoar,mgsolv,mgdisc,iinfo,errtol,ipkey,
     4   omegal,omegan,irite,iperf)
c* *********************************************************************
c* purpose:
c*
c*    this routine reads in some initial values.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          iparm(*),nrwk,niwk,nx,ny,nz,nlev,nu1,nu2,mgkey
      integer          itmax,istop,ipcon,nonlin,iinfo,irite,ipkey,iperf
      integer          mgsmoo,mgprol,mgcoar,mgsolv,mgdisc
      double precision rparm(*),errtol,omegal,omegan
c*
c*    *** encode iparm parameters ***
      iparm(1)  = nrwk
      iparm(2)  = niwk
      iparm(3)  = nx
      iparm(4)  = ny
      iparm(5)  = nz
      iparm(6)  = nlev
      iparm(7)  = nu1
      iparm(8)  = nu2
      iparm(9)  = mgkey
      iparm(10) = itmax
      iparm(11) = istop
      iparm(12) = iinfo
      iparm(13) = irite
      iparm(14) = ipkey
      iparm(15) = ipcon
      iparm(16) = nonlin
      iparm(17) = mgprol
      iparm(18) = mgcoar
      iparm(19) = mgdisc
      iparm(20) = mgsmoo
      iparm(21) = mgsolv
      iparm(22) = iperf
c*
c*    *** encode rparm parameters ***
      rparm(1)  = errtol
      rparm(9)  = omegal
      rparm(10) = omegan
c*
c*    *** return and end ***
      return
      end
