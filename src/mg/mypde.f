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
c* author:  Nathan Baker and michael holst 
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
c* author:  Nathan Baker and michael holst
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
c*                 -(zkappa2/bulkIonicStrength)*ionConc[i]*ionQ[i]
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
c*    print *, charge
c*    print *, kappa

      return
      end




      function c_scal(coef,u,ipkey)
c* *********************************************************************
c* purpose:
c*
c*    define the nonlinearity (scalar version)
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      double precision coef,u,c_scal,am_zero,am_neg,am_pos,argument
      double precision poly,fact
      integer          ipkey,ideg
      integer          ichopped,ichopped_neg,ichopped_pos
      integer          MAXPOLY
      double precision ZSMALL,ZLARGE,SINH_MIN,SINH_MAX
      parameter        (MAXPOLY = 50)
      parameter        (ZSMALL   = 1.0e-20, ZLARGE = 1.0e20)
      parameter        (SINH_MIN = -85.0, SINH_MAX = 85.0)
c*
c*    *** check if full sinh requested ***
      if (ipkey .eq. 0) then
c*
c*       *** initialize chopped counter ***
         ichopped = 0
c*
c*       *** am_zero is 0 if coef zero, and 1 if coef nonzero ***
         am_zero = dmin1(ZSMALL,dabs(coef)) * ZLARGE
c*
c*       *** am_neg is chopped u if u negative, 0 if u positive ***
         am_neg = dmax1( dmin1(u,0.0d0) , SINH_MIN)
c*
c*       *** am_neg is chopped u if u positive, 0 if u negative ***
         am_pos = dmin1( dmax1(u,0.0d0) , SINH_MAX)
c*
c*       *** finally determine the function value ***
         argument = am_zero * ( am_neg + am_pos )
         c_scal = coef * sinh ( argument )
c*
c*       *** count chopped values ***
         ichopped_neg = idnint(aint(am_neg / SINH_MIN))
         ichopped_pos = idnint(aint(am_pos / SINH_MAX))
         ichopped = ichopped
     2       + idnint(am_zero) * (ichopped_neg + ichopped_pos)
c*
c*    *** else if polynomial requested ***
      elseif ((ipkey .gt. 1) .and. (mod(ipkey,2) .eq. 1) .and.
     2        (ipkey .le. MAXPOLY)) then
c*
c*       *** am_zero is 0 if coef zero, and 1 if coef nonzero ***
         am_zero = dmin1(ZSMALL,dabs(coef)) * ZLARGE
         argument = am_zero * u
c*
c*       *** sum the polynomial (large terms first for stability) ***
         poly = argument
         fact = 1.0d0
         do 40 ideg = 3,ipkey,2
            fact = fact * dble(ideg-1)
            fact = fact * ideg
            poly = poly + (1.0d0/fact) * argument**ideg
 40      continue
         c_scal = coef * poly
c*
c*    *** else return linear approximation ***
      else
         c_scal = coef * u
c* ***** print*,'% C_SCAL: using linear approximation...'
      endif
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
      double precision poly,fact
      integer          ipkey,ideg
      integer          ichopped,ichopped_neg,ichopped_pos
      integer          MAXPOLY
      double precision ZSMALL,ZLARGE,SINH_MIN,SINH_MAX
      parameter        (MAXPOLY = 50)
      parameter        (ZSMALL   = 1.0e-20, ZLARGE = 1.0e20)
      parameter        (SINH_MIN = -85.0, SINH_MAX = 85.0)
c*
c*    *** check if full sinh requested ***
      if (ipkey .eq. 0) then
c*
c*       *** initialize chopped counter ***
         ichopped = 0
c*
c*       *** am_zero is 0 if coef zero, and 1 if coef nonzero ***
         am_zero = dmin1(ZSMALL,dabs(coef)) * ZLARGE
c*
c*       *** am_neg is chopped u if u negative, 0 if u positive ***
         am_neg = dmax1( dmin1(u,0.0d0) , SINH_MIN)
c*
c*       *** am_neg is chopped u if u positive, 0 if u negative ***
         am_pos = dmin1( dmax1(u,0.0d0) , SINH_MAX)
c*
c*       *** finally determine the function value ***
         argument = am_zero * ( am_neg + am_pos )
         dc_scal = coef * cosh ( argument )
c*
c*       *** count chopped values ***
         ichopped_neg = idnint(aint(am_neg / SINH_MIN))
         ichopped_pos = idnint(aint(am_pos / SINH_MAX))
         ichopped = ichopped
     2       + idnint(am_zero) * (ichopped_neg + ichopped_pos)
c*
c*    *** else if polynomial requested ***
      elseif ((ipkey .gt. 1) .and. (mod(ipkey,2) .eq. 1) .and.
     2        (ipkey .le. MAXPOLY)) then
c*
c*       *** am_zero is 0 if coef zero, and 1 if coef nonzero ***
         am_zero = dmin1(ZSMALL,dabs(coef)) * ZLARGE
         argument = am_zero * u
c*
c*       *** sum the polynomial (large terms first for stability) ***
         poly = 1.0d0
         fact = 1.0d0
         do 40 ideg = 3,ipkey,2
            fact = fact * dble(ideg-1)
            fact = fact * ideg
            poly = poly + ideg*(1.0d0/fact)*argument**(ideg-1)
 40      continue
         dc_scal = coef * poly
c*
c*    *** else return linear approximation ***
      else
         dc_scal = coef 
c* ***** print*,'% DC_SCAL: using linear approximation...'
      endif
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
      double precision uin(*),uout(*),coef(*)
      double precision am_zero,am_neg,am_pos,argument,poly,fact
      integer          ichopped,ichopped_neg,ichopped_pos
      integer          MAXPOLY
      double precision ZSMALL,ZLARGE,SINH_MIN,SINH_MAX
      parameter        (MAXPOLY = 50)
      parameter        (ZSMALL   = 1.0e-20, ZLARGE = 1.0e20)
      parameter        (SINH_MIN = -85.0, SINH_MAX = 85.0)
      integer          n,i,ii,nproc,ipara,ivect
      parameter        (nproc=1)
c*
c*    *** find parallel loops (ipara), remainder (ivect) ***
      n     = nx * ny * nz
      ipara = n / nproc
      ivect = mod(n,nproc)
c*
c*
c*    *** check if full sinh requested ***
      if (ipkey .eq. 0) then
c*
c*       *** info ***
c* ***** print*,'% C_VEC: using full exponential'
c*
c*       *** initialize chopped counter ***
         ichopped = 0
c*
c*       *** do parallel loops ***
         do 10 ii = 1, nproc
            do 11 i = 1+(ipara*(ii-1)), ipara*ii
c*
c*             *** am_zero is 0 if coef zero, and 1 if coef nonzero ***
               am_zero = dmin1(ZSMALL,dabs(coef(i))) * ZLARGE
c*
c*             *** am_neg is chopped u if u negative, 0 if u positive ***
               am_neg = dmax1( dmin1(uin(i),0.0d0) , SINH_MIN)
c*
c*             *** am_neg is chopped u if u positive, 0 if u negative ***
               am_pos = dmin1( dmax1(uin(i),0.0d0) , SINH_MAX)
c*
c*             *** finally determine the function value ***
               argument = am_zero * ( am_neg + am_pos )
               uout(i) = coef(i) * sinh( argument )
c*
c*             *** count chopped values ***
               ichopped_neg = idnint(aint(am_neg / SINH_MIN))
               ichopped_pos = idnint(aint(am_pos / SINH_MAX))
               ichopped = ichopped
     2            + idnint(am_zero) * (ichopped_neg + ichopped_pos)
 11         continue
 10      continue
c*
c*       *** do vector loops ***
         do 20 i = ipara*nproc+1, n
c*          *** am_zero is 0 if coef zero, and 1 if coef nonzero ***
            am_zero = dmin1(ZSMALL,dabs(coef(i))) * ZLARGE
c*
c*          *** am_neg is chopped u if u negative, 0 if u positive ***
            am_neg = dmax1( dmin1(uin(i),0.0d0) , SINH_MIN)
c*
c*          *** am_neg is chopped u if u positive, 0 if u negative ***
            am_pos = dmin1( dmax1(uin(i),0.0d0) , SINH_MAX)
c*
c*          *** finally determine the function value ***
            argument = am_zero * ( am_neg + am_pos )
            uout(i) = coef(i) * sinh( argument )
c*
c*          *** count chopped values ***
            ichopped_neg = idnint(aint(am_neg / SINH_MIN))
            ichopped_pos = idnint(aint(am_pos / SINH_MAX))
            ichopped = ichopped
     2         + idnint(am_zero) * (ichopped_neg + ichopped_pos)
 20      continue
c*
c*       *** info ***
         if (ichopped .gt. 0) then
         print*,'% C_VEC: trapped sinh overflows:          ',ichopped
         endif
c*
c*    *** else if polynomial requested ***
      elseif ((ipkey .gt. 1) .and. (mod(ipkey,2) .eq. 1) .and.
     2        (ipkey .le. MAXPOLY)) then
c*
c*       *** create the polynomial term ***
c*
c*       *** info ***
c* ***** print*,'% C_VEC: using polynomial of degree:      ',ipkey
c*
c*       *** do parallel loops ***
         do 30 ii = 1, nproc
            do 31 i = 1+(ipara*(ii-1)), ipara*ii
c*
c*             *** am_zero is 0 if coef zero, and 1 if coef nonzero ***
               am_zero = dmin1(ZSMALL,dabs(coef(i))) * ZLARGE
               argument = am_zero * uin(i)
c*
c*             *** sum the polynomial (large terms first for stability) ***
               poly = argument
               fact = 1.0d0
               do 32 ideg = 3,ipkey,2
                  fact = fact * dble(ideg-1)
                  fact = fact * ideg
                  poly = poly + (1.0d0/fact) * argument**ideg
 32            continue
               uout(i) = coef(i) * poly
 31         continue
 30      continue
c*
c*       *** do vector loops ***
         do 40 i = ipara*nproc+1, n
c*
c*          *** am_zero is 0 if coef zero, and 1 if coef nonzero ***
            am_zero = dmin1(ZSMALL,dabs(coef(i))) * ZLARGE
            argument = am_zero * uin(i)
c*
c*          *** sum the polynomial (large terms first for stability) ***
            poly = argument
            fact = 1.0d0
            do 41 ideg = 3,ipkey,2
               fact = fact * dble(ideg-1)
               fact = fact * ideg
               poly = poly + (1.0d0/fact) * argument**ideg
 41         continue
            uout(i) = coef(i) * poly
 40      continue
c*
c*    *** else return linear approximation ***
      else
c*
c*       *** info ***
c* ***** print*,'% C_VEC: using linear for degree given:   ',ipkey
c*
c*       *** do parallel loops ***
         do 50 ii = 1, nproc
            do 51 i = 1+(ipara*(ii-1)), ipara*ii
c*
c*             *** am_zero is 0 if coef zero, and 1 if coef nonzero ***
               am_zero = dmin1(ZSMALL,dabs(coef(i))) * ZLARGE
               argument = am_zero * uin(i)
c*
c*             *** make the linear term ***
               uout(i) = coef(i) * argument
 51         continue
 50      continue
c*
c*       *** do vector loops ***
         do 60 i = ipara*nproc+1, n
c*
c*          *** am_zero is 0 if coef zero, and 1 if coef nonzero ***
            am_zero = dmin1(ZSMALL,dabs(coef(i))) * ZLARGE
            argument = am_zero * uin(i)
c*
c*          *** make the linear term ***
            uout(i) = coef(i) * argument
 60      continue
c*
c*    *** end if ***
      endif
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
      integer          nx,ny,nz,ipkey,ideg
      double precision uin(*),uout(*),coef(*)
      double precision am_zero,am_neg,am_pos,argument,poly,fact
      integer          ichopped,ichopped_neg,ichopped_pos
      integer          MAXPOLY
      double precision ZSMALL,ZLARGE,SINH_MIN,SINH_MAX
      parameter        (MAXPOLY = 50)
      parameter        (ZSMALL   = 1.0e-20, ZLARGE = 1.0e20)
      parameter        (SINH_MIN = -85.0, SINH_MAX = 85.0)
      integer          n,i,ii,nproc,ipara,ivect
      parameter        (nproc=1)
c*
c*    *** find parallel loops (ipara), remainder (ivect) ***
      n     = nx * ny * nz
      ipara = n / nproc
      ivect = mod(n,nproc)
c*
c*
c*    *** check if full sinh requested ***
      if (ipkey .eq. 0) then
c*
c*       *** info ***
c* ***** print*,'% DC_VEC: using full exponential'
c*
c*       *** initialize chopped counter ***
         ichopped = 0
c*
c*       *** do parallel loops ***
         do 10 ii = 1, nproc
            do 11 i = 1+(ipara*(ii-1)), ipara*ii
c*
c*             *** am_zero is 0 if coef zero, and 1 if coef nonzero ***
               am_zero = dmin1(ZSMALL,dabs(coef(i))) * ZLARGE
c*
c*             *** am_neg is chopped u if u negative, 0 if u positive ***
               am_neg = dmax1( dmin1(uin(i),0.0d0) , SINH_MIN)
c*
c*             *** am_neg is chopped u if u positive, 0 if u negative ***
               am_pos = dmin1( dmax1(uin(i),0.0d0) , SINH_MAX)
c*
c*             *** finally determine the function value ***
               argument = am_zero * ( am_neg + am_pos )
               uout(i) = coef(i) * cosh( argument )
c*
c*             *** count chopped values ***
               ichopped_neg = idnint(aint(am_neg / SINH_MIN))
               ichopped_pos = idnint(aint(am_pos / SINH_MAX))
               ichopped = ichopped
     2            + idnint(am_zero) * (ichopped_neg + ichopped_pos)
 11         continue
 10      continue
c*
c*       *** do vector loops ***
         do 20 i = ipara*nproc+1, n
c*          *** am_zero is 0 if coef zero, and 1 if coef nonzero ***
            am_zero = dmin1(ZSMALL,dabs(coef(i))) * ZLARGE
c*
c*          *** am_neg is chopped u if u negative, 0 if u positive ***
            am_neg = dmax1( dmin1(uin(i),0.0d0) , SINH_MIN)
c*
c*          *** am_neg is chopped u if u positive, 0 if u negative ***
            am_pos = dmin1( dmax1(uin(i),0.0d0) , SINH_MAX)
c*
c*          *** finally determine the function value ***
            argument = am_zero * ( am_neg + am_pos )
            uout(i) = coef(i) * cosh( argument )
c*
c*          *** count chopped values ***
            ichopped_neg = idnint(aint(am_neg / SINH_MIN))
            ichopped_pos = idnint(aint(am_pos / SINH_MAX))
            ichopped = ichopped
     2         + idnint(am_zero) * (ichopped_neg + ichopped_pos)
 20      continue
c*
c*       *** info ***
         if (ichopped .gt. 0) then
         print*,'% DC_VEC: trapped sinh overflows:         ',ichopped
         endif
c*
c*    *** else if polynomial requested ***
      elseif ((ipkey .gt. 1) .and. (mod(ipkey,2) .eq. 1) .and.
     2        (ipkey .le. MAXPOLY)) then
c*
c*       *** create the polynomial term ***
c*
c*       *** info ***
c* ***** print*,'% DC_VEC: using polynomial of degree:     ',ipkey
c*
c*       *** do parallel loops ***
         do 30 ii = 1, nproc
            do 31 i = 1+(ipara*(ii-1)), ipara*ii
c*
c*             *** am_zero is 0 if coef zero, and 1 if coef nonzero ***
               am_zero = dmin1(ZSMALL,dabs(coef(i))) * ZLARGE
               argument = am_zero * uin(i)
c*
c*             *** sum the polynomial (large terms first for stability) ***
               poly = 1.0d0
               fact = 1.0d0
               do 32 ideg = 3,ipkey,2
                  fact = fact * dble(ideg-1)
                  fact = fact * ideg
                  poly = poly + ideg*(1.0d0/fact) * argument**(ideg-1)
 32            continue
               uout(i) = coef(i) * poly
 31         continue
 30      continue
c*
c*       *** do vector loops ***
         do 40 i = ipara*nproc+1, n
c*
c*          *** am_zero is 0 if coef zero, and 1 if coef nonzero ***
            am_zero = dmin1(ZSMALL,dabs(coef(i))) * ZLARGE
            argument = am_zero * uin(i)
c*
c*          *** sum the polynomial (large terms first for stability) ***
            poly = 1.0d0
            fact = 1.0d0
            do 41 ideg = 3,ipkey,2
               fact = fact * dble(ideg-1)
               fact = fact * ideg
               poly = poly + ideg*(1.0d0/fact) * argument**(ideg-1)
 41         continue
            uout(i) = coef(i) * poly
 40      continue
c*
c*    *** else return linear approximation ***
      else
c*
c*       *** info ***
c* ***** print*,'% DC_VEC: using linear for degree given:  ',ipkey
c*
c*       *** do parallel loops ***
         do 50 ii = 1, nproc
            do 51 i = 1+(ipara*(ii-1)), ipara*ii
c*
c*             *** make the linear term ***
               uout(i) = coef(i)
 51         continue
 50      continue
c*
c*       *** do vector loops ***
         do 60 i = ipara*nproc+1, n
c*
c*          *** make the linear term ***
            uout(i) = coef(i)
 60      continue
c*
c*    *** end if ***
      endif
c*
c*    *** end it ***
      return
      end
