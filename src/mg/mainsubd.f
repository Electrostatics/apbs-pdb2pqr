      subroutine readit (iparm,rparm,nx,ny,nz,nlev,nrwk,niwk,key,meth)
c* *********************************************************************
c* purpose:
c*
c*    this routine reads in some initial values.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          iparm(*),nrwk,niwk,nx,ny,nz,nlev,nu1,nu2,mgkey
      integer          istop,iinfo,key,meth,itmax,ipcon
      integer          ipkey,iperf
      integer          nonlin,mgprol,mgcoar,mgdisc,mgsmoo,mgsolv
      double precision rparm(*),errtol,omegal,omegan
c*
c*    *** parameters ***
      integer          iread,irite
      parameter        (iread=7, irite=8)
c*
c*    *** if not interactive mode then open i/o files ***
      open(unit=iread,  file='in',   status='unknown')
      open(unit=irite,  file='outt', status='unknown')
      rewind(iread)
      rewind(irite)
c*
c*    *** input the controling parameters ***
      read (iread,10)
      read (iread,10)
      read (iread,30) nx
      read (iread,30) ny
      read (iread,30) nz
      read (iread,20) errtol
      read (iread,30) itmax
      read (iread,30) istop
      read (iread,30) iinfo
      read (iread,30) ipkey
      read (iread,30) key
      read (iread,30) iperf
      read (iread,10)
      read (iread,10)
      read (iread,10)
      read (iread,30) meth
      read (iread,30) nonlin
      read (iread,30) mgkey
      read (iread,30) nlev
      read (iread,30) nu1
      read (iread,30) nu2
      read (iread,30) mgsmoo
      read (iread,30) mgprol
      read (iread,30) mgcoar
      read (iread,30) mgsolv
      read (iread,30) mgdisc
      read (iread,20) omegal
      read (iread,20) omegan
      read (iread,30) ipcon
c*
c*    *** pack iparm/rparm correctly for desired method ***
      call packmg (iparm,rparm,
     2   nrwk,niwk,nx,ny,nz,nlev,nu1,nu2,mgkey,itmax,istop,ipcon,
     3   nonlin,mgsmoo,mgprol,mgcoar,mgsolv,mgdisc,iinfo,errtol,ipkey,
     4   omegal,omegan,irite,iperf)
c*
c*    *** do a little output now ***
      write(*,40) '% READIT: done reading input file... '
c*
c*    *** format statements ***
 10   format()
 20   format(e10.1)
 30   format(i10)
 40   format (a)
c*
c*    *** return and end ***
      return
      end
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
      subroutine writit (iparm,rparm,nx,ny,nz,u,
     2   xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf,key)
c* *********************************************************************
c* purpose:
c*
c*    this routine prints out the solution
c*
c* author:  michael holst
c* *********************************************************************
      implicit     none
c*
c*    *** dimensions of parameters ***
      integer          iparm(*),nx,ny,nz,key,i,j,k
      double precision rparm(*),xmin,xmax,ymin,ymax,zmin,zmax,hx,hy,hz
      double precision xf(*),yf(*),zf(*)
      double precision u(nx,ny,nz),tcf(nx,ny,nz)
      double precision a1cf(nx,ny,nz),a2cf(nx,ny,nz),a3cf(nx,ny,nz)
      double precision ccf(nx,ny,nz),fcf(nx,ny,nz)
      double precision gxcf (ny,nz,2)
      double precision gycf (nx,nz,2)
      double precision gzcf (nx,ny,2)
      character*20     title
      integer          iounit
      parameter        (iounit=9)
c*
c*    *** get some stuff from iparm and rparm ***
      xmin   = rparm(3)
      xmax   = rparm(4)
      ymin   = rparm(5)
      ymax   = rparm(6)
      zmin   = rparm(7)
      zmax   = rparm(8)
c*
c*    *** compute a few things ***
      hx = (xmax-xmin) / dble(nx-1)
      hy = (ymax-ymin) / dble(ny-1)
      hz = (zmax-zmin) / dble(nz-1)
      title = 'elliptic multigrid'
c*
c*    *** print out the solution values in the ncsa format ***
      if (key.eq.1) then
c*
c*       *** open and rewind the file ***
         open(unit=iounit, file='outn', status='unknown') 
         rewind(iounit)
c*
c*       *** write it out ***
         write(iounit,605) title
         write(iounit,600) nx,ny,nz,0.0d0
         do 10 k = 1, nz
            do 11 j = 1, ny
               do 12 i = 1, nx
                  write (iounit,610) u(i,j,k)
 12            continue
 11         continue
 10      continue
c*
c*       *** close the file ***
         close(iounit)
c*
c*    *** print out the solution values in a normal format ***
      elseif (key.eq.2) then
c*
c*       *** open and rewind the file ***
         open(unit=iounit, file='outn', status='unknown') 
         rewind(iounit)
c*
c*       *** write it out ***
         write(iounit,605) title
         write(iounit,600) nx,ny,nz,0.0d0
         do 20 k = 1, nz
            do 21 j = 1, ny
               do 22 i = 1, nx
                  write (iounit,620) xf(i),yf(j),zf(k),u(i,j,k)
 22            continue
 21         continue
 20      continue
c*
c*       *** close the file ***
         close(iounit)
c*
c*    *** print out the solution values in a second normal format ***
      elseif (key.eq.3) then
c*
c*       *** open and rewind the file ***
         open(unit=iounit, file='outn', status='unknown') 
         rewind(iounit)
c*
c*       *** write it out ***
         write(iounit,605) title
         write(iounit,600) nx,ny,nz,0.0d0
         do 30 k = 1, nz
            do 31 j = 1, ny
               do 32 i = 1, nx
                  write (iounit,630) i,j,k,u(i,j,k)
 32            continue
 31         continue
 30      continue
c*
c*       *** close the file ***
         close(iounit)
c*
c*    *** print out the solution values to check against true solution ***
      elseif (key.eq.4) then
c*
c*       *** open and rewind the file ***
         open(unit=iounit, file='outn', status='unknown') 
         rewind(iounit)
c*
c*       *** write it out ***
         do 40 i = 33, 33
            do 41 j = 33, 33
               do 42 k = 28, 38
                  write (6,630) i,j,k,u(i,j,k),tcf(i,j,k)
 42            continue
 41         continue
 40      continue
         write(6,*) ' '
         do 50 i = 1, 1
            do 51 j = 1, 1
               do 52 k = 1, 3
                  write (6,630) i,j,k,u(i,j,k),tcf(i,j,k)
 52            continue
 51         continue
 50      continue
         write(6,*) ' '
         do 60 i = 2, 2
            do 61 j = 2, 2
               do 62 k = 2, 4
                  write (6,630) i,j,k,u(i,j,k),tcf(i,j,k)
 62            continue
 61         continue
 60      continue
c*
c*       *** close the file ***
         close(iounit)
c*
c*    *** write out some slices only ***
      elseif (key.eq.5) then
c*
c*       *** re-fill the coefficient arrays since we overwrote them ***
         call fillco (iparm,rparm,nx,ny,nz,
     2      xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf)
c*
c*       *** open and rewind the file ***
         open(unit=iounit, file='outnx.m', status='unknown') 
         rewind(iounit)
c*
c*       *** write it out ***
         write(iounit,650) title
         write(iounit,660) nx
         write(iounit,670) 'xi','yj','zk','x','a','c','f','t'
         write(iounit,605) 'Z = ['
         do 70 k = nz/2, nz/2
            do 71 j = ny/2, ny/2
               do 72 i = 1, nx
                  write (iounit,640) xf(i),yf(j),zf(k),u(i,j,k),
     2               (a1cf(i,j,k)+a1cf(i+1,j,k)
     3               +a2cf(i,j,k)+a2cf(i,j+1,k)
     4               +a3cf(i,j,k)+a3cf(i,j,k+1))/6.0d0,
     5               ccf(i,j,k),fcf(i,j,k),tcf(i,j,k)
 72            continue
 71         continue
 70      continue
         write(iounit,605) '];'
c*
c*       *** close the file ***
         close(iounit)
c*
c*       *** open and rewind the file ***
         open(unit=iounit, file='outny.m', status='unknown') 
         rewind(iounit)
c*
c*       *** write it out ***
         write(iounit,650) title
         write(iounit,661) ny
         write(iounit,670) 'xi','yj','zk','x','a','c','f','t'
         write(iounit,605) 'Z = ['
         do 80 k = nz/2, nz/2
            do 81 j = 1, ny
               do 82 i = nx/2, nx/2
                  write (iounit,640) xf(i),yf(j),zf(k),u(i,j,k),
     2               (a1cf(i,j,k)+a1cf(i+1,j,k)
     3               +a2cf(i,j,k)+a2cf(i,j+1,k)
     4               +a3cf(i,j,k)+a3cf(i,j,k+1))/6.0d0,
     5               ccf(i,j,k),fcf(i,j,k),tcf(i,j,k)
 82            continue
 81         continue
 80      continue
         write(iounit,605) '];'
c*
c*       *** close the file ***
         close(iounit)
c*
c*       *** open and rewind the file ***
         open(unit=iounit, file='outnz.m', status='unknown') 
         rewind(iounit)
c*
c*       *** write it out ***
         write(iounit,650) title
         write(iounit,662) nz
         write(iounit,670) 'xi','yj','zk','x','a','c','f','t'
         write(iounit,605) 'Z = ['
         do 90 k = 1, nz
            do 91 j = ny/2, ny/2
               do 92 i = nx/2, nx/2
                  write (iounit,640) xf(i),yf(j),zf(k),u(i,j,k),
     2               (a1cf(i,j,k)+a1cf(i+1,j,k)
     3               +a2cf(i,j,k)+a2cf(i,j+1,k)
     4               +a3cf(i,j,k)+a3cf(i,j,k+1))/6.0d0,
     5               ccf(i,j,k),fcf(i,j,k),tcf(i,j,k)
 92            continue
 91         continue
 90      continue
         write(iounit,605) '];'
c*
c*       *** close the file ***
         close(iounit)
c*
c*    *** write out in a fast binary format ***
      elseif (key.eq.6) then
c*
c*       *** open and rewind the file ***
         open(unit=iounit, file='outn', form='unformatted') 
         rewind(iounit)
c*
c*       *** write it out ***
         write(iounit) title
         write(iounit) nx,ny,nz,0.0d0
         write(iounit) u
c*
c*       *** close the file ***
         close(iounit)
c*
c*    *** print out the mesh for display/etc. ***
      elseif (key.eq.7) then
c*
c*       *** open and rewind the file ***
         open(unit=iounit, file='outn', status='unknown') 
         rewind(iounit)
c*
c*       *** write it out ***
         write(iounit,605) "the mesh"
         write(iounit,680) nx,ny,nz
         do 100 i = 1, nx
            write (iounit,690) xf(i)
 100     continue
         do 110 j = 1, ny
            write (iounit,690) yf(j)
 110     continue
         do 120 k = 1, nz
            write (iounit,690) zf(k)
 120     continue
c*
c*       *** close the file ***
         close(iounit)
c*
c*    *** that's all ***
      endif
c*
c*    *** format statements ***
 600  format (3(i10),2(1pe10.2))
 605  format (a)
 610  format (2(1pe10.2))
 620  format (5(1pe10.2))
 630  format (3(i10),2(1pe15.7))
 640  format (8(1pe10.2))
 650  format ('% ',a)
 660  format ('nx = ',i10,';')
 661  format ('ny = ',i10,';')
 662  format ('nz = ',i10,';')
 670  format ('% ',a8,7(a10))
 680  format (3(i10))
 690  format (1pe10.2)
c*
c*    *** return and end ***
      return
      end
