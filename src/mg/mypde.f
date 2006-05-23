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

      subroutine mypdefinit(tnion,tcharge,tsconc)
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
c*    tsconc  = prefactor for counterion Boltzmann distribution terms,
c*              basically a scaled concentration:
c*                 -(ion concentration/bulkIonicStrength)/2
c*
c* author: Nathan Baker
c* *********************************************************************
      integer          nion, tnion
      integer          MAXION
      integer          i
      parameter        (MAXION = 50)
      double precision charge(MAXION), sconc(MAXION) 
      double precision tcharge(*), tsconc(*)
      common /MYPDEF/  charge, sconc
      common /MYPDEF/  nion

      nion = tnion
      if (nion.gt.MAXION) then
          call vnmprt(2, 'mypdef_init: Error: too many ion species', 43)
          call vnmprt(2, 'mypdef_init:  Ignoring the extra ones!', 38)
          nion = MAXION
      endif
      do 10 i = 1, nion
          charge(i) = tcharge(i)
          sconc(i) = tsconc(i)
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
      double precision charge(MAXION), sconc(MAXION) 
      common /MYPDEF/  charge, sconc
      common /MYPDEF/  nion

      do 10 i = 1, nion
          charge(i) = 0.
          sconc(i) = 0.
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
      parameter        (ZSMALL   = 1.0d-20, ZLARGE = 1.0d20)
      parameter        (SINH_MIN = -85.0, SINH_MAX = 85.0)
c* Added by NAB to allow different ions with different charges
      integer          nion
      integer          MAXION
      parameter        (MAXION = 50)
      double precision charge(MAXION)
      double precision sconc(MAXION)
      common /MYPDEF/ charge
      common /MYPDEF/ sconc
      common /MYPDEF/ nion

      call vnmprt(2, 'MYPDEF: C_SCAL NOT SUPPORTED! USE NEWTON SOLVER',
     2     47)
      stop

c*    print *, 'HELLO FROM C_SCAL'
c*      c_scal2 = 0.0
c*
c*    Loop over all the ions
c*      do 39 iion = 1, nion
c*       Assemble the ion-specific coefficient
c*         coef2 = -1.0 * coef * sconc(iion)
c*       Assemble the ion-specific potential value
c*         u2 = -1.0 * u * charge(iion)
c*
c*       *** check if full exp requested ***
c*         if (ipkey .eq. 0) then
c*
c*         *** initialize chopped counter ***
c*           ichopped = 0
c*
c*         *** am_zero is 0 if coef2 zero, and 1 if coef nonzero ***
c*           am_zero = dmin1(ZSMALL,dabs(coef2)) * ZLARGE
c*
c*         *** am_neg is chopped u if u negative, 0 if u positive ***
c*           am_neg = dmax1( dmin1(u2,0.0d0) , SINH_MIN)
c*
c*         *** am_neg is chopped u if u positive, 0 if u negative ***
c*           am_pos = dmin1( dmax1(u2,0.0d0) , SINH_MAX)
c*
c*         *** finally determine the function value ***
c*           argument = am_zero * ( am_neg + am_pos )
c*           c_scal2 = c_scal2 + coef2 * exp ( argument )
c*
c*         *** count chopped values ***
c*           ichopped_neg = idnint(aint(am_neg / SINH_MIN))
c*           ichopped_pos = idnint(aint(am_pos / SINH_MAX))
c*           ichopped = ichopped
c*     2         + idnint(am_zero) * (ichopped_neg + ichopped_pos)
c*
c*      *** else if polynomial requested ***
c*        elseif ((ipkey .gt. 1) .and. (mod(ipkey,2) .eq. 1) .and.
c*     2          (ipkey .le. MAXPOLY)) then
c*          call vnmprt(2, 'MYPDEF: POLYNOMIAL APPROXIMATION UNAVAILABLE',
c*     2      44)
c*      *** else return linear approximation ***
c*        else
c*           c_scal2 = c_scal2 + coef2 * u2
c* *****   print*,'% C_SCAL: using linear approximation...'
c*        endif
c*  39   continue
c*
c*      c_scal = c_scal2
c*
c*
c* ****** info ***
c* ***if (ichopped .gt. 0) then
c* ***   print*,'% C_SCAL: trapped hyperbolic sine overflow...'
c* ***endif
c*
c*    *** end it ***
      c_scal = 0
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
      parameter        (ZSMALL   = 1.0d-20, ZLARGE = 1.0d20)
      parameter        (SINH_MIN = -85.0, SINH_MAX = 85.0)
c* Added by NAB to allow different ions with different charges
      integer          nion
      integer          MAXION
      parameter        (MAXION = 50)
      double precision charge(MAXION)
      double precision sconc(MAXION)
      common /MYPDEF/ charge
      common /MYPDEF/ sconc
      common /MYPDEF/ nion

      call vnmprt(2, 'MYPDEF: DC_SCAL NOT SUPPORTED! USE NEWTON SOLVER',
     2    48)
      stop

c*      dc_scal2 = 0.0
c*    print *, 'HELLO FROM DC_SCAL'
c*
c*      do 39 iion = 1, nion
c*        coef2 = coef * sconc(iion) * charge(iion)
c*        u2 = -1.0 * u * charge(iion)
c*
c*      *** check if full exp requested ***
c*        if (ipkey .eq. 0) then
c*
c*         *** initialize chopped counter ***
c*           ichopped = 0
c*
c*         *** am_zero is 0 if coef zero, and 1 if coef nonzero ***
c*           am_zero = dmin1(ZSMALL,dabs(coef2)) * ZLARGE
c*
c*         *** am_neg is chopped u if u negative, 0 if u positive ***
c*           am_neg = dmax1( dmin1(u2,0.0d0) , SINH_MIN)
c*
c*         *** am_neg is chopped u if u positive, 0 if u negative ***
c*           am_pos = dmin1( dmax1(u2,0.0d0) , SINH_MAX)
c*
c*         *** finally determine the function value ***
c*           argument = am_zero * ( am_neg + am_pos )
c*           dc_scal2 = dc_scal2 + coef2 * exp ( argument )
c*
c*         *** count chopped values ***
c*           ichopped_neg = idnint(aint(am_neg / SINH_MIN))
c*           ichopped_pos = idnint(aint(am_pos / SINH_MAX))
c*           ichopped = ichopped
c*     2         + idnint(am_zero) * (ichopped_neg + ichopped_pos)
c*
c*      *** else if polynomial requested ***
c*        elseif ((ipkey .gt. 1) .and. (mod(ipkey,2) .eq. 1) .and.
c*     2          (ipkey .le. MAXPOLY)) then
c*          call vnmprt(2, 'MYPDEF: POLYNOMIAL APPROXIMATION UNAVAILABLE',
c*     2      44)
c*
c*      *** else return linear approximation ***
c*        else
c*           dc_scal2 = dc_scal2 + coef*sconc(iion)*charge(iion)
c* ***** print*,'% DC_SCAL: using linear approximation...'
c*        endif
c*
c*  39   continue
c*
c*      dc_scal = dc_scal2
c*
c* ****** info ***
c* ***if (ichopped .gt. 0) then
c* ***   print*,'% DC_SCAL: trapped hyperbolic cosine overflow...'
c* ***endif
c*
c*    *** end it ***
      dc_scal = 0
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
      double precision uin(*),uout(*),coef(*),zcf2,zu2
      double precision am_zero,am_neg,am_pos,argument,poly,fact
      integer          ichopped,ichopped_neg,ichopped_pos,iion
      integer          MAXPOLY
      double precision ZSMALL,ZLARGE,SINH_MIN,SINH_MAX
      parameter        (MAXPOLY = 50)
      parameter        (ZSMALL   = 1.0d-20, ZLARGE = 1.0d20)
      parameter        (SINH_MIN = -85.0, SINH_MAX = 85.0)
      integer          n,i,ii,nproc,ipara,ivect
      parameter        (nproc=1)
c* Added by NAB to allow different ions with different charges
      integer          nion
      integer          MAXION
      parameter        (MAXION = 50)
      double precision charge(MAXION)
      double precision sconc(MAXION)
      common /MYPDEF/ charge
      common /MYPDEF/ sconc
      common /MYPDEF/ nion
 
c*    print *, 'HELLO FROM C_VEC'

c*
c*    *** find parallel loops (ipara), remainder (ivect) ***
      n     = nx * ny * nz
      ipara = n / nproc
      ivect = mod(n,nproc)

      do 38 i = 1, n
        uout(i) = 0
 38   continue

      do 39 iion = 1, nion
c*      Assemble the ion-specific coefficient
        zcf2 = -1.0 * sconc(iion) * charge(iion)
c*      Assemble the ion-specific potential value
        zu2 = -1.0 * charge(iion)

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
           do 10 ii = 1, nproc
              do 11 i = 1+(ipara*(ii-1)), ipara*ii
c*
c*               *** am_zero is 0 if coef zero, and 1 if coef nonzero**
                 am_zero = dmin1(ZSMALL,dabs(zcf2*coef(i)))*ZLARGE
c*
c*               *** am_neg is chopped u if u negative, 0 if u positive**
                 am_neg = dmax1( dmin1(zu2*uin(i),0.0d0) , SINH_MIN)
c*
c*               *** am_neg is chopped u if u positive, 0 if u negative ***
                 am_pos = dmin1( dmax1(zu2*uin(i),0.0d0) , SINH_MAX)
c*
c*               *** finally determine the function value ***
                 argument = am_zero * ( am_neg + am_pos )
                 uout(i) = uout(i) + zcf2 * coef(i) * exp( argument )
c*
c*               *** count chopped values ***
                 ichopped_neg = idnint(aint(am_neg / SINH_MIN))
                 ichopped_pos = idnint(aint(am_pos / SINH_MAX))
                 ichopped = ichopped
     2              + idnint(am_zero) * (ichopped_neg + ichopped_pos)
 11           continue
 10        continue
c*
c*         *** do vector loops ***
           do 20 i = ipara*nproc+1, n
c*            *** am_zero is 0 if coef zero, and 1 if coef nonzero ***
              am_zero = dmin1(ZSMALL,dabs(zcf2*coef(i))) * ZLARGE
c*
c*            *** am_neg is chopped u if u negative, 0 if u positive ***
              am_neg = dmax1( dmin1(zu2*uin(i),0.0d0) , SINH_MIN)
c*
c*            *** am_neg is chopped u if u positive, 0 if u negative ***
              am_pos = dmin1( dmax1(zu2*uin(i),0.0d0) , SINH_MAX)
c*
c*            *** finally determine the function value ***
              argument = am_zero * ( am_neg + am_pos )
              uout(i) = uout(i) + zcf2*coef(i)*exp( argument )
c*
c*            *** count chopped values ***
              ichopped_neg = idnint(aint(am_neg / SINH_MIN))
              ichopped_pos = idnint(aint(am_pos / SINH_MAX))
              ichopped = ichopped
     2           + idnint(am_zero) * (ichopped_neg + ichopped_pos)
 20        continue
c*
c*         *** info ***
           if (ichopped .gt. 0) then
           call vnmpri(2,'% C_VEC: trapped exp overflows:          ',
     2      41, ichopped)
           endif
c*
c*      *** else if polynomial requested ***
        elseif ((ipkey .gt. 1) .and. (mod(ipkey,2) .eq. 1) .and.
     2          (ipkey .le. MAXPOLY)) then
          call vnmprt(2, 'MYPDEF: POLYNOMIAL APPROXIMATION UNAVAILABLE',
     2      44)
          stop
c*
c*      *** else return linear approximation ***
        else
           call vnmprt(2, 'MYPDEF: LINEAR APPROXIMATION UNAVAILABLE',
     2      40)
           stop
c*
c*         *** info ***
c* *****   print*,'% C_VEC: using linear for degree given:   ',ipkey
c*
c*         *** do parallel loops ***
c*           do 50 ii = 1, nproc
c*              do 51 i = 1+(ipara*(ii-1)), ipara*ii
c*
c*               *** am_zero is 0 if coef zero, and 1 if coef nonzero ***
c*                 am_zero = dmin1(ZSMALL,dabs(zcf2*coef(i))) * ZLARGE
c*                 argument = am_zero * zu2 * uin(i)
c*
c*               *** make the linear term ***
c*                 uout(i) = uout(i) + zcf2 * coef(i) * argument
c* 51         continue
c* 50      continue
c*
c*       *** do vector loops ***
c*         do 60 i = ipara*nproc+1, n
c*
c*          *** am_zero is 0 if coef zero, and 1 if coef nonzero ***
c*            am_zero = dmin1(ZSMALL,dabs(coef(i))) * ZLARGE
c*            argument = am_zero * uin(i)
c*
c*          *** make the linear term ***
c*            uout(i) = uout(i) + zcf2 * coef(i) * argument
c* 60      continue
c*
c*      *** end if ***
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
      double precision uin(*),uout(*),coef(*),zcf2,zu2
      double precision am_zero,am_neg,am_pos,argument,poly,fact
      integer          ichopped,ichopped_neg,ichopped_pos
      integer          MAXPOLY
      double precision ZSMALL,ZLARGE,SINH_MIN,SINH_MAX
      parameter        (MAXPOLY = 50)
      parameter        (ZSMALL   = 1.0d-20, ZLARGE = 1.0d20)
      parameter        (SINH_MIN = -85.0, SINH_MAX = 85.0)
      integer          n,i,ii,nproc,ipara,ivect
      parameter        (nproc=1)
c* Added by NAB to allow different ions with different charges
      integer          nion
      integer          MAXION
      parameter        (MAXION = 50)
      double precision charge(MAXION)
      double precision sconc(MAXION)
      common /MYPDEF/ charge
      common /MYPDEF/ sconc
      common /MYPDEF/ nion

c*
c*    *** find parallel loops (ipara), remainder (ivect) ***
      n     = nx * ny * nz
      ipara = n / nproc
      ivect = mod(n,nproc)

c*    print *, 'HELLO FROM DC_VEC'

      do 38 i = 1, n
        uout(i) = 0.0
 38   continue

      do 39 iion = 1, nion
        zcf2 = sconc(iion) * charge(iion) * charge(iion)
        zu2 = -1.0 * charge(iion)
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
           do 10 ii = 1, nproc
              do 11 i = 1+(ipara*(ii-1)), ipara*ii
c*
c*               *** am_zero is 0 if coef zero, and 1 if coef nonzero ***
                 am_zero = dmin1(ZSMALL,dabs(zcf2*coef(i)))*ZLARGE
c*
c*               *** am_neg is chopped u if u negative, 0 if u positive ***
                 am_neg = dmax1(dmin1(zu2*uin(i),0.0d0),SINH_MIN)
c*
c*               *** am_neg is chopped u if u positive, 0 if u negative ***
                 am_pos = dmin1(dmax1(zu2*uin(i),0.0d0) , SINH_MAX)
c*
c*               *** finally determine the function value ***
                 argument = am_zero * ( am_neg + am_pos )
                 uout(i) = uout(i) + zcf2*coef(i)*exp( argument )
c*
c*               *** count chopped values ***
                 ichopped_neg = idnint(aint(am_neg / SINH_MIN))
                 ichopped_pos = idnint(aint(am_pos / SINH_MAX))
                 ichopped = ichopped
     2              + idnint(am_zero) * (ichopped_neg + ichopped_pos)
 11           continue
 10        continue
c*
c*         *** do vector loops ***
           do 20 i = ipara*nproc+1, n
c*            *** am_zero is 0 if coef zero, and 1 if coef nonzero ***
              am_zero = dmin1(ZSMALL,dabs(zcf2*coef(i))) * ZLARGE
c*
c*            *** am_neg is chopped u if u negative, 0 if u positive ***
              am_neg = dmax1( dmin1(zu2*uin(i),0.0d0) , SINH_MIN)
c*
c*            *** am_neg is chopped u if u positive, 0 if u negative ***
              am_pos = dmin1( dmax1(zu2*uin(i),0.0d0) , SINH_MAX)
c*
c*            *** finally determine the function value ***
              argument = am_zero * ( am_neg + am_pos )
              uout(i) = uout(i) + zcf2*coef(i) * exp( argument )
c*
c*            *** count chopped values ***
              ichopped_neg = idnint(aint(am_neg / SINH_MIN))
              ichopped_pos = idnint(aint(am_pos / SINH_MAX))
              ichopped = ichopped
     2           + idnint(am_zero) * (ichopped_neg + ichopped_pos)
 20        continue
c*
c*         *** info ***
           if (ichopped .gt. 0) then
           call vnmpri(2,'% DC_VEC: trapped exp overflows:         ',
     2       41, ichopped)
           endif
c*
c*      *** else if polynomial requested ***
        elseif ((ipkey .gt. 1) .and. (mod(ipkey,2) .eq. 1) .and.
     2          (ipkey .le. MAXPOLY)) then

          call vnmprt(2, 'MYPDEF: POLYNOMIAL APPROXIMATION UNAVAILABLE',
     2      44)
          stop
*
c*      *** else return linear approximation ***
        else
           call vnmprt(2, 'MYPDEF: LINEAR APPROXIMATION UNAVAILABLE',
     2      40)
           stop
c*
c*         *** info ***
c* *****   print*,'% DC_VEC: using linear for degree given:  ',ipkey
c*
c*         *** do parallel loops ***
c*           do 50 ii = 1, nproc
c*              do 51 i = 1+(ipara*(ii-1)), ipara*ii
c*
c*               *** make the linear term ***
c*                 uout(i) = uout(i)+coef(i)*sconc(iion)*charge(iion)
c* 51           continue
c* 50        continue
c*
c*         *** do vector loops ***
c*           do 60 i = ipara*nproc+1, n
c*
c*            *** make the linear term ***
c*              uout(i) = uout(i) + coef(i)*sconc(iion)*charge(iion)
c* 60        continue
c*
c*      *** end if ***
        endif

 39   continue
c*
c*    *** end it ***
      return
      end
