c* ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
c*
c* DIFFERENTIAL EQUATION: Poisson-Boltzmann
c*
c* BOUNDARY CONDITIONS:
c*
c*     East   Face (xmin):  Dirichlet, homogeneous
c*     West   Face (xmax):  Dirichlet, homogeneous
c*     North  Face (ymin):  Dirichlet, homogeneous
c*     South  Face (ymax):  Dirichlet, homogeneous
c*     Top    Face (zmin):  Dirichlet, homogeneous
c*     Bottom Face (zmax):  Dirichlet, homogeneous
c*
c* MESH:                  hx  = (xmax-xmin) / (nx-1)
c*                        hy  = (ymax-ymin) / (ny-1)
c*                        hz  = (zmax-zmin) / (nz-1)
c*                        xi = xmin + (i-1) * hx,  i=1,...,nx
c*                        yi = ymin + (j-1) * hy,  j=1,...,ny
c*                        zi = zmin + (k-1) * hk,  k=1,...,nz
c*
c* CHOSEN TRUE SOLUTION:  numerical solution provided by delphi tools
c*
c* author:  michael holst
c* ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
c*
c* *********************************************************************
c* notes:
c*
c*    to control overflow in the hyperbolic and exp functions, note
c*    that the following are the argument limits of the various 
c*    functions on various machines after which overflow occurs:
c*
c*    Convex C240, Sun 3/60, Sun SPARC, IBM RS/6000:
c*
c*       sinh, cosh, exp:     maximal argument (abs value) =  88.0d0
c*       dsinh, dcosh, dexp:  maximal argument (abs value) = 709.0d0
c*
c* author:  michael holst
c* *********************************************************************

      subroutine mypdefinit(tnion,tcharge,tkappa)
c* *********************************************************************
c* Purpose:
c*
c*    Set up the ionic species to be used in later calculations.  This 
c*    must be called before any other of the routines in this file.
c*  
c* Arguments:
c*
c*    tnion   = number of ionic species
c*    tcharge = charge in electrons
c*    tkappa  = prefactor for counterion Boltzmann distribution terms,
c*              basically:
c*                 -(zkappa2/bulkIonicStrength)*ionConc[i]
c*
c* author: Nathan Baker
c* *********************************************************************
      integer          nion, tnion
      integer          MAXION
      integer          i
      parameter        (MAXION = 50)
      double precision charge(MAXION), kappa(MAXION) 
      double precision tcharge(*), tkappa(*)
      common /MYPDEF/  charge, kappa
      common /MYPDEF/  nion

      nion = tnion
      if (nion.gt.MAXION) then
          call vnmprt(2, 'mypdef_init: Error: too many ion species', 43)
          call vnmprt(2, 'mypdef_init:  Ignoring the extra ones!', 38)
          nion = MAXION
      endif
      do 10 i = 1, nion
          charge(i) = tcharge(i)
          kappa(i) = tkappa(i)
10    continue

      return
      end

      subroutine mypdefclear()
c* *********************************************************************
c* Purpose:
c*
c*    Clears out arrays
c*  
c* author: Nathan Baker
c* *********************************************************************
      integer          nion
      integer          MAXION
      integer          i
      parameter        (MAXION = 50)
      double precision charge(MAXION), kappa(MAXION) 
      common /MYPDEF/  charge, kappa
      common /MYPDEF/  nion

      do 10 i = 1, nion
          charge(i) = 0.
          kappa(i) = 0.
10    continue
      nion = 0

      return
      end

 
      function c_scal(coef,u,ipkey)
c* *********************************************************************
c* purpose:
c*
c*    define the nonlinearity (scalar version)
c*
c* author:  Nathan Baker and Michael Holst
c* *********************************************************************
      implicit         none
      double precision coef,u,c_scal,am_zero,am_neg,am_pos,argument
      double precision poly,fact,coef2,c_scal2,u2
      integer          ipkey,ideg,iion
      integer          ichopped,ichopped_neg,ichopped_pos
      integer          MAXPOLY
      double precision ZSMALL,ZLARGE,SINH_MIN,SINH_MAX
      parameter        (MAXPOLY = 50)
      parameter        (ZSMALL   = 1.0e-20, ZLARGE = 1.0e20)
      parameter        (SINH_MIN = -85.0, SINH_MAX = 85.0)
c* Added by NAB to allow different ions with different charges
      integer          nion
      integer          MAXION
      parameter        (MAXION = 50)
      double precision charge(MAXION)
      double precision kappa(MAXION)
      common /MYPDEF/ charge
      common /MYPDEF/ kappa
      common /MYPDEF/ nion

      call vnmprt(2, 'MYPDEF: C_SCAL NOT SUPPORTED! USE NEWTON SOLVER',
     2     47)
      stop

cc*    print *, 'HELLO FROM C_SCAL'
c      c_scal2 = 0.0
c
cc*    Loop over all the ions
c      do 39 iion = 1, nion
cc*       Assemble the ion-specific coefficient
c         coef2 = -1.0 * coef * kappa(iion)
cc*       Assemble the ion-specific potential value
c         u2 = -1.0 * u * charge(iion)
cc*
cc*       *** check if full exp requested ***
c         if (ipkey .eq. 0) then
cc*
cc*         *** initialize chopped counter ***
c           ichopped = 0
cc*
cc*         *** am_zero is 0 if coef2 zero, and 1 if coef nonzero ***
c           am_zero = dmin1(ZSMALL,dabs(coef2)) * ZLARGE
cc*
cc*         *** am_neg is chopped u if u negative, 0 if u positive ***
c           am_neg = dmax1( dmin1(u2,0.0d0) , SINH_MIN)
cc*
cc*         *** am_neg is chopped u if u positive, 0 if u negative ***
c           am_pos = dmin1( dmax1(u2,0.0d0) , SINH_MAX)
cc*
cc*         *** finally determine the function value ***
c           argument = am_zero * ( am_neg + am_pos )
c           c_scal2 = c_scal2 + coef2 * exp ( argument )
cc*
cc*         *** count chopped values ***
c           ichopped_neg = idnint(aint(am_neg / SINH_MIN))
c           ichopped_pos = idnint(aint(am_pos / SINH_MAX))
c           ichopped = ichopped
c     2         + idnint(am_zero) * (ichopped_neg + ichopped_pos)
cc*
cc*      *** else if polynomial requested ***
c        elseif ((ipkey .gt. 1) .and. (mod(ipkey,2) .eq. 1) .and.
c     2          (ipkey .le. MAXPOLY)) then
c          call vnmprt(2, 'MYPDEF: POLYNOMIAL APPROXIMATION UNAVAILABLE',
c     2      44)
cc*      *** else return linear approximation ***
c        else
c           c_scal2 = c_scal2 + coef2 * u2
cc* *****   print*,'% C_SCAL: using linear approximation...'
c        endif
c 39   continue
c
c      c_scal = c_scal2
c
c*
c* ****** info ***
c* ***if (ichopped .gt. 0) then
c* ***   print*,'% C_SCAL: trapped hyperbolic sine overflow...'
c* ***endif
c*
c*    *** end it ***
      return
      end
      function dc_scal(coef,u,ipkey)
c* *********************************************************************
c* purpose:
c*
c*    define the derivative of the nonlinearity (scalar version)
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      double precision coef,u,dc_scal,am_zero,am_neg,am_pos,argument
      double precision poly,fact,coef2,u2,dc_scal2
      integer          ipkey,ideg, iion
      integer          ichopped,ichopped_neg,ichopped_pos
      integer          MAXPOLY
      double precision ZSMALL,ZLARGE,SINH_MIN,SINH_MAX
      parameter        (MAXPOLY = 50)
      parameter        (ZSMALL   = 1.0e-20, ZLARGE = 1.0e20)
      parameter        (SINH_MIN = -85.0, SINH_MAX = 85.0)
c* Added by NAB to allow different ions with different charges
      integer          nion
      integer          MAXION
      parameter        (MAXION = 50)
      double precision charge(MAXION)
      double precision kappa(MAXION)
      common /MYPDEF/ charge
      common /MYPDEF/ kappa
      common /MYPDEF/ nion

      call vnmprt(2, 'MYPDEF: DC_SCAL NOT SUPPORTED! USE NEWTON SOLVER',
     2    48)
      stop

      dc_scal2 = 0.0
cc*    print *, 'HELLO FROM DC_SCAL'
c
c      do 39 iion = 1, nion
c        coef2 = coef * kappa(iion) * charge(iion)
c        u2 = -1.0 * u * charge(iion)
cc*
cc*      *** check if full exp requested ***
c        if (ipkey .eq. 0) then
cc*
cc*         *** initialize chopped counter ***
c           ichopped = 0
cc*
cc*         *** am_zero is 0 if coef zero, and 1 if coef nonzero ***
c           am_zero = dmin1(ZSMALL,dabs(coef2)) * ZLARGE
cc*
cc*         *** am_neg is chopped u if u negative, 0 if u positive ***
c           am_neg = dmax1( dmin1(u2,0.0d0) , SINH_MIN)
cc*
cc*         *** am_neg is chopped u if u positive, 0 if u negative ***
c           am_pos = dmin1( dmax1(u2,0.0d0) , SINH_MAX)
cc*
cc*         *** finally determine the function value ***
c           argument = am_zero * ( am_neg + am_pos )
c           dc_scal2 = dc_scal2 + coef2 * exp ( argument )
cc*
cc*         *** count chopped values ***
c           ichopped_neg = idnint(aint(am_neg / SINH_MIN))
c           ichopped_pos = idnint(aint(am_pos / SINH_MAX))
c           ichopped = ichopped
c     2         + idnint(am_zero) * (ichopped_neg + ichopped_pos)
cc*
cc*      *** else if polynomial requested ***
c        elseif ((ipkey .gt. 1) .and. (mod(ipkey,2) .eq. 1) .and.
c     2          (ipkey .le. MAXPOLY)) then
c          call vnmprt(2, 'MYPDEF: POLYNOMIAL APPROXIMATION UNAVAILABLE',
c     2      44)
cc*
cc*      *** else return linear approximation ***
c        else
c           dc_scal2 = dc_scal2 + coef*kappa(iion)*charge(iion)
cc* ***** print*,'% DC_SCAL: using linear approximation...'
c        endif
c
c 39   continue
c
c      dc_scal = dc_scal2
c*
c* ****** info ***
c* ***if (ichopped .gt. 0) then
c* ***   print*,'% DC_SCAL: trapped hyperbolic cosine overflow...'
c* ***endif
c*
c*    *** end it ***
      return
      end
      subroutine c_vec(coef,uin,uout,nx,ny,nz,ipkey)
c* *********************************************************************
c* purpose:
c*
c*    define the nonlinearity (vector version)
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          nx,ny,nz,ipkey,ideg
      double precision uin(*),uout(*),coef(*),zlin,zexp
      double precision am_zero,am_neg,am_pos,argument,poly,fact
      integer          ichopped,ichopped_neg,ichopped_pos,iion
      integer          MAXPOLY
      double precision ZSMALL,ZLARGE,SINH_MIN,SINH_MAX
      parameter        (MAXPOLY = 50)
      parameter        (ZSMALL   = 1.0e-20, ZLARGE = 1.0e20)
      parameter        (SINH_MIN = -85.0, SINH_MAX = 85.0)
      integer          n,i,ii
c* Added by NAB to allow different ions with different charges
      integer          nion
      integer          MAXION
      parameter        (MAXION = 50)
      double precision charge(MAXION)
      double precision kappa(MAXION)
      common /MYPDEF/ charge
      common /MYPDEF/ kappa
      common /MYPDEF/ nion
 
c*    print *, 'HELLO FROM C_VEC'

c*
c*    *** find parallel loops (ipara), remainder (ivect) ***
      n     = nx * ny * nz

      do 38 i = 1, n
        uout(i) = 0.0d0
 38   continue

      do 39 iion = 1, nion
c*      Assemble the linear multiplier
c*        zlin * coef(i) * exp (zexp * u(i) )
        zlin = -1.0 * kappa(iion) * charge(iion)
c*      Assemble the exponent multiplier:  
c*        zlin * coef(i) * exp (zexp * u(i) )
        zexp = -1.0 * charge(iion)

c*
c*
c*      *** check if full exp requested ***
        if (ipkey .eq. 0) then
c*
c*         *** info ***
c* ***** print*,'% C_VEC: using full exponential'
c*
c*         *** initialize chopped counter ***
           ichopped = 0
c*
c*         *** do parallel loops ***
           do 10 i = 1, n
               uout(i) = uout(i) + zlin * coef(i) * 
     2             exp( zexp*uin(i) )
c*
 10        continue
c*
c*      *** else if polynomial requested ***
        elseif ((ipkey .gt. 1) .and. (mod(ipkey,2) .eq. 1) .and.
     2          (ipkey .le. MAXPOLY)) then
          call vnmprt(2, 'MYPDEF: POLYNOMIAL APPROXIMATION UNAVAILABLE',
     2      44)
c*
c*      *** else return linear approximation ***
        else
c*
c*       *** do parallel loops ***
         do 50 i = 1, n
c*           *** make the linear term ***
             uout(i) = uout(i) + zlin * coef(i) * zexp * uin(i)
 50      continue
c*
        endif
 39   continue
c*
c*    *** end it ***
      return
      end
      subroutine dc_vec(coef,uin,uout,nx,ny,nz,ipkey)
c* *********************************************************************
c* purpose:
c*
c*    define the derivative of the nonlinearity (vector version)
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          nx,ny,nz,ipkey,ideg,iion
      double precision uin(*),uout(*),coef(*),zlin,zexp
      double precision am_zero,am_neg,am_pos,argument,poly,fact
      integer          ichopped,ichopped_neg,ichopped_pos
      integer          MAXPOLY
      double precision ZSMALL,ZLARGE,SINH_MIN,SINH_MAX
      parameter        (MAXPOLY = 50)
      parameter        (ZSMALL   = 1.0e-20, ZLARGE = 1.0e20)
      parameter        (SINH_MIN = -85.0, SINH_MAX = 85.0)
      integer          n,i,ii
c* Added by NAB to allow different ions with different charges
      integer          nion
      integer          MAXION
      parameter        (MAXION = 50)
      double precision charge(MAXION)
      double precision kappa(MAXION)
      common /MYPDEF/ charge
      common /MYPDEF/ kappa
      common /MYPDEF/ nion

c*
c*    *** find parallel loops (ipara), remainder (ivect) ***
      n     = nx * ny * nz

c*    print *, 'HELLO FROM DC_VEC'

      do 38 i = 1, n
        uout(i) = 0.0
 38   continue

      do 39 iion = 1, nion
c*      Assemble the linear multiplier
c*        zlin * zexp * coef(i) * exp (zexp * u(i) )
        zlin = -1.0 * kappa(iion) * charge(iion)
c*      Assemble the exponent multiplier:  
c*        zlin * zexp * coef(i) * exp (zexp * u(i) )
        zexp = -1.0 * charge(iion)
c*
c*
c*      *** check if full exp requested ***
        if (ipkey .eq. 0) then
c*
c*         *** info ***
c* *****   print*,'% DC_VEC: using full exponential'
c*
c*         *** initialize chopped counter ***
           ichopped = 0
c*
c*         *** do parallel loops ***
           do 10 i = 1, n
c*
             uout(i) = uout(i)+zlin*zexp*coef(i)*exp(zexp*uin(i))
c*
 10        continue
c*
c*      *** else if polynomial requested ***
        elseif ((ipkey .gt. 1) .and. (mod(ipkey,2) .eq. 1) .and.
     2          (ipkey .le. MAXPOLY)) then

          call vnmprt(2, 'MYPDEF: POLYNOMIAL APPROXIMATION UNAVAILABLE',
     2      44)
*
c*      *** else return linear approximation ***
        else
c*
c*         *** do parallel loops ***
           do 50 i = 1, n
             uout(i) = uout(i) + zlin * coef(i) * zexp * uin(i)
 50        continue
c*      *** end if ***
        endif

 39   continue
c*
c*    *** end it ***
      return
      end
