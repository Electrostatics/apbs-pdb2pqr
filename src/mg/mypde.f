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

      subroutine mypdefinitlpbe(tnion,tcharge,tsconc)
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

	  double precision v1, v2, v3, conc1, conc2, conc3, vol, relSize
      common /MYPDEF/  v1, v2, v3, conc1, conc2, conc3, vol, relSize

      common /MYPDEF/  charge, sconc
      common /MYPDEF/  nion
	  
	  nion = tnion
	  if (nion.gt.MAXION) then
		  call vnmprt(2, 'mypdef_init: Error: too many ion 
     .						species', 43)
		  call vnmprt(2, 'mypdef_init:  Ignoring the extra 
     .						ones', 38)
		  nion = MAXION
	  endif
	  do i = 1, nion
		  charge(i) = tcharge(i)
		  sconc(i) = tsconc(i)
	  end do
	  
      return
      end

      subroutine mypdefinitnpbe(tnion,tcharge,tsconc)
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
	  
	  double precision v1, v2, v3, conc1, conc2, conc3, vol, relSize
      common /MYPDEF/  v1, v2, v3, conc1, conc2, conc3, vol, relSize
	  
      common /MYPDEF/  charge, sconc
      common /MYPDEF/  nion
		  
	  nion = tnion
	  if (nion.gt.MAXION) then
		  call vnmprt(2, 'mypdef_init: Error: too many ion 
     .						species', 43)
		  call vnmprt(2, 'mypdef_init:  Ignoring the extra 
     .						ones', 38)
		  nion = MAXION
	  endif
	  do i = 1, nion
		  charge(i) = tcharge(i)
		  sconc(i) = tsconc(i)
	  end do

      return
      end

      subroutine mypdefinitsmpbe(tnion,tcharge,tsconc,smvolume,smsize)
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
	  
	  integer		   ipkey
	  double precision smvolume, smsize

	  double precision v1, v2, v3, conc1, conc2, conc3, vol, relSize
      common /MYPDEF/  v1, v2, v3, conc1, conc2, conc3, vol, relSize

      common /MYPDEF/  charge, sconc
      common /MYPDEF/  nion
		  
	  !nion = tnion
	  if (tnion.gt.3) then
		 call vnmprt(2, 'SMPBE: modified theory handles only', 36)
		 call vnmprt(2, 'SMPBE: three ion species. Ignoring ', 36)
		 call vnmprt(2, 'SMPBE: the rest of the ions!       ', 36)
		 call vnmprt(2, 'SMPBE: (mypde.f::mypdefinit)       ', 36)
	  endif

	  v1 = tcharge(1)
	  v2 = tcharge(2)
	  v3 = tcharge(3)
	  conc1 = tsconc(1)
	  conc2 = tsconc(2)
	  conc3 = tsconc(3)
	  
	  vol = smvolume
	  relSize = smsize

	  !call vnmprd(2, '!! v1        = ', 15, v1)
	  !call vnmprd(2, '!! v2        = ', 15, v2)
	  !call vnmprd(2, '!! v3        = ', 15, v3)
	  !call vnmprd(2, '!! conc1     = ', 15, conc1)
	  !call vnmprd(2, '!! conc2     = ', 15, conc2)
	  !call vnmprd(2, '!! conc3     = ', 15, conc3)
	  !call vnmprd(2, '!! SMPBE Vol = ', 15, smvolume)
	  !call vnmprd(2, '!! SMPBE Size= ', 15, smsize)
	  
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
	  
	  !double precision v1, v2, v3, conc1, conc2, conc3, vol, relSize
      !common /MYPDEF/  v1, v2, v3, conc1, conc2, conc3, vol, relSize
	  
      common /MYPDEF/  charge, sconc
      common /MYPDEF/  nion

	  ! THIS FUNCTION CALL IS DEPRECATED DUE TO A BUG IN THE COMMON BLOCK
	  ! USAGE. IT WILL BE MYPDEFINIT WILL BE REPLACED WITH EXPLICIT CALLS
	  ! FOR EACH TYPE OF PDE INITIALIZATION

      do i = 1, nion
          charge(i) = 0.
          sconc(i) = 0.
	  end do
      nion = 0
	  
	  !v1 = 0.0
	  !v2 = 0.0
	  !v3 = 0.0
	  !conc1 = 0.0
	  !conc2 = 0.0
	  !conc3 = 0.0
	  !vol = 0.0
	  !relSize = 0.0
	  
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
	  
	  if(ipkey .eq. -2) then
		! Do Nothing
		! This was added since we don't know what the code path is for
		! for non SMPBE calculations (in particular NPBE)
	  else
		  call vnmprt(2, 'MYPDEF: C_SCAL NOT SUPPORTED. 
     2     USE NEWTON SOLVER', 47)
		  stop

		  c_scal = 0
	  endif
	  
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

      if(ipkey .eq. -2) then
		! Do Nothing
		! This was added since we don't know what the code path is for
		! for non SMPBE calculations
	  else
		  call vnmprt(2, 'MYPDEF: DC_SCAL NOT SUPPORTED. 
     2     USE NEWTON SOLVER', 47)
		  stop

		  dc_scal = 0
	  endif

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
	  
	  !SMPBE Added
	  if(ipkey .eq. -2) then
		  call c_vecsmpbe(coef,uin,uout,nx,ny,nz,ipkey)
	  else
		  call c_vecpmg(coef,uin,uout,nx,ny,nz,ipkey)
	  end if
	  
      return
      end

      subroutine c_vecpmg(coef,uin,uout,nx,ny,nz,ipkey)
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
      integer          n,i

c* Added by NAB to allow different ions with different charges
      integer          nion
      integer          MAXION
      parameter        (MAXION = 50)
      double precision charge(MAXION)
      double precision sconc(MAXION)
	  
	  double precision v1, v2, v3, conc1, conc2, conc3, vol, relSize
      common /MYPDEF/  v1, v2, v3, conc1, conc2, conc3, vol, relSize
	  
      common /MYPDEF/ charge
      common /MYPDEF/ sconc
      common /MYPDEF/ nion

      !*** find parallel loops (ipara), remainder (ivect) !***
      n     = nx * ny * nz

      do i = 1, n
        uout(i) = 0
	  end do
	  
	  do iion = 1, nion
		!Assemble the ion-specific coefficient
		zcf2 = -1.0 * sconc(iion) * charge(iion)
		!Assemble the ion-specific potential value
		zu2 = -1.0 * charge(iion)

		!*** check if full exp requested !***
		if (ipkey .eq. 0) then
  
		   !*** initialize chopped counter !***
		   ichopped = 0
		   
!$OMP  PARALLEL DO 
!$OMP& DEFAULT(shared) 
!$OMP& PRIVATE(i,ichopped_neg,ichopped_pos,am_zero,am_neg,
!$OMP&          am_pos,argument)
!$OMP& REDUCTION (+:ichopped)  
		   !*** do parallel loops !***
		   do i = 1, n
  
			 !*** am_zero is 0 if coef zero, and 1 if coef nonzero**
			 am_zero = dmin1(ZSMALL,dabs(zcf2*coef(i)))*ZLARGE

			 !*** am_neg is chopped u if u negative, 0 if u positive**
			 am_neg = dmax1( dmin1(zu2*uin(i),0.0d0) , SINH_MIN)

			 !*** am_neg is chopped u if u positive, 0 if u negative !***
			 am_pos = dmin1( dmax1(zu2*uin(i),0.0d0) , SINH_MAX)

			 !*** finally determine the function value !***
			 argument = am_zero * ( am_neg + am_pos )
			 uout(i) = uout(i) + zcf2 * coef(i) * exp(argument)

			 !*** count chopped values !***
			 ichopped_neg = idnint(aint(am_neg / SINH_MIN))
			 ichopped_pos = idnint(aint(am_pos / SINH_MAX))
			 ichopped = ichopped + idnint(am_zero) *
     2                        (ichopped_neg + ichopped_pos)
		   end do
!$OMP END PARALLEL DO
		   
		   !*** info !***
		   if (ichopped .gt. 0) then
		   call vnmpri(2,'% C_VEC: trapped exp overflows:          ',
     2						41, ichopped)
		   endif
  
		!*** else if polynomial requested !***
		elseif ((ipkey .gt. 1) .and. (mod(ipkey,2) .eq. 1) .and.
     2          (ipkey .le. MAXPOLY)) then
		  call vnmprt(2, 'MYPDEF: POLYNOMIAL APPROXIMATION
     2                UNAVAILABLE',44)
		  stop
  
		!*** else return linear approximation !***
		else
		   call vnmprt(2, 'MYPDEF: LINEAR APPROXIMATION
     2                UNAVAILABLE',40)
		   stop
		endif
	  end do

      !*** end it !***
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
	  
	  ! SMPBE Added
 	  if(ipkey .eq. -2) then 
		  call dc_vecsmpbe(coef,uin,uout,nx,ny,nz,ipkey)
 	  else
 		  call dc_vecpmg(coef,uin,uout,nx,ny,nz,ipkey)
 	  end if
	  
      return
      end

      subroutine dc_vecpmg(coef,uin,uout,nx,ny,nz,ipkey)
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
      integer          n,i
c* Added by NAB to allow different ions with different charges
      integer          nion
      integer          MAXION
      parameter        (MAXION = 50)
      double precision charge(MAXION)
      double precision sconc(MAXION)

      double precision v1, v2, v3, conc1, conc2, conc3, vol, relSize
      common /MYPDEF/  v1, v2, v3, conc1, conc2, conc3, vol, relSize
	  
      common /MYPDEF/ charge
      common /MYPDEF/ sconc
      common /MYPDEF/ nion

c*
c*    *** find parallel loops (ipara), remainder (ivect) ***
      n     = nx * ny * nz

      do i = 1, n
        uout(i) = 0.0
	  end do
	   
	  do iion = 1, nion
		zcf2 = sconc(iion) * charge(iion) * charge(iion)
		zu2 = -1.0 * charge(iion)
  
		!*** check if full exp requested !***
		if (ipkey .eq. 0) then

		   !*** initialize chopped counter !***
		   ichopped = 0
  
!$OMP  PARALLEL DO 
!$OMP& DEFAULT(shared) 
!$OMP& PRIVATE(i,ichopped_neg,ichopped_pos,am_zero,am_neg,
!$OMP&          am_pos,argument)
!$OMP& REDUCTION (+:ichopped) 
		   !*** do parallel loops !***
		   do i = 1, n

			 !*** am_zero is 0 if coef zero, and 1 if coef nonzero !***
			 am_zero = dmin1(ZSMALL,dabs(zcf2*coef(i)))*ZLARGE

			 !*** am_neg is chopped u if u negative, 0 if u positive !***
			 am_neg = dmax1(dmin1(zu2*uin(i),0.0d0),SINH_MIN)

			 !*** am_neg is chopped u if u positive, 0 if u negative !***
			 am_pos = dmin1(dmax1(zu2*uin(i),0.0d0) , SINH_MAX)

			 !*** finally determine the function value !***
			 argument = am_zero * ( am_neg + am_pos )
			 uout(i) = uout(i) + zcf2*coef(i)*exp( argument )

			 !*** count chopped values !***
			 ichopped_neg = idnint(aint(am_neg / SINH_MIN))
			 ichopped_pos = idnint(aint(am_pos / SINH_MAX))
			 ichopped = ichopped + idnint(am_zero) * 
     2								(ichopped_neg + ichopped_pos)
		   end do
!$OMP END PARALLEL DO
  
		   !*** info !***
		   if (ichopped .gt. 0) then
		   call vnmpri(2,'% DC_VEC: trapped exp overflows:   ',
     2					41, ichopped)
		   endif
  
		!*** else if polynomial requested !***
		elseif ((ipkey .gt. 1) .and. (mod(ipkey,2) .eq. 1) .and.
     2          (ipkey .le. MAXPOLY)) then

		  call vnmprt(2, 'MYPDEF: POLYNOMIAL APPROXIMATION
     2				UNAVAILABLE',44)
		  stop
		  
		!*** else return linear approximation !***
		else
		   call vnmprt(2, 'MYPDEF: LINEAR APPROXIMATION
     2				UNAVAILABLE',40)
		   stop

		endif
	  end do

	  
      return
      end

      subroutine c_vecsmpbe(coef,uin,uout,nx,ny,nz,ipkey)
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

c* Added by DG SMPBE variables and common blocks
	  double precision Na, fracOccA, fracOccB, fracOccC, phi, ionStr
      double precision z1, z2, z3, ca, cb, cc, a, k
      double precision a1_neg, a1_pos, a2_neg, a2_pos
      double precision a3_neg, a3_pos, a1, a2, a3
      double precision f, g, gpark, alpha
      parameter        (Na = 6.022045000e-04)
	  
	  double precision v1, v2, v3, conc1, conc2, conc3, vol, relSize
      common /MYPDEF/  v1, v2, v3, conc1, conc2, conc3, vol, relSize
	  
      common /MYPDEF/ charge
      common /MYPDEF/ sconc
      common /MYPDEF/ nion

      !*** find parallel loops (ipara), remainder (ivect) !***
      n     = nx * ny * nz
      ipara = n / nproc
      ivect = mod(n,nproc)

      do i = 1, n
        uout(i) = 0
	  end do
	  
	  !*** initialize the chopped counter !***
	  ichopped = 0

	  
	  z1 = v1
	  z2 = v2
	  z3 = v3
	  ca = conc1
	  cb = conc2
	  cc = conc3
	  a = vol
	  k = relSize

	  if ((k-1).lt. ZSMALL) then
	  call vnmprt(2, '!!  C_VEC: k=1, using special routine', 37) 
	  endif

	  !derived quantities
	  fracOccA = Na*ca*a**3
	  fracOccB = Na*cb*a**3
	  fracOccC = Na*cc*a**3
	  phi = (fracOccA/k) + fracOccB + fracOccC
	  alpha = (fracOccA/k)/(1-phi)
	  ionStr = 0.5*(ca*z1**2 + cb*z2**2 + cc*z3**2)

	  do i = 1, n

		 am_zero = dmin1(ZSMALL,dabs(coef(i)))*ZLARGE

		 !compute the arguments for exp(-z*u) term
		 a1_neg = dmax1( dmin1(-1.0*z1*uin(i),0.0d0) , SINH_MIN)
		 a1_pos = dmin1( dmax1(-1.0*z1*uin(i),0.0d0) , SINH_MAX)

		 !compute the arguments for exp(-u) term
		 a2_neg = dmax1( dmin1(-1.0*z2*uin(i),0.0d0) , SINH_MIN)
		 a2_pos = dmin1( dmax1(-1.0*z2*uin(i),0.0d0) , SINH_MAX)

		 !compute the arguments for exp(u) term
		 a3_neg = dmax1( dmin1(-1.0*z3*uin(i),0.0d0) , SINH_MIN)
		 a3_pos = dmin1( dmax1(-1.0*z3*uin(i),0.0d0) , SINH_MAX)

		 a1 = am_zero * (a1_neg + a1_pos)
		 a2 = am_zero * (a2_neg + a2_pos)
		 a3 = am_zero * (a3_neg + a3_pos)

		 gpark = ((1+alpha*exp(a1))/(1+alpha))

		 if ((k-1).lt. ZSMALL) then
		 f = z1*ca*exp(a1) + z2*cb*exp(a2) + z3*cc*exp(a3)
		 g = 1-phi+fracOccA*exp(a1) + fracOccB*exp(a2)
     2        + fracOccC*exp(a3)
		 else
		 f = z1*ca*exp(a1)*gpark**(k-1)+z2*cb*exp(a2)+z3*cc*exp(a3)
		 g = (1-phi+(fracOccA/k))*gpark**k 
     2        + fracOccB*exp(a2)+fracOccC*exp(a3)
		 endif

		 uout(i) = -1.0*coef(i)*(0.5/ionStr)*(f/g)

		 !*** count chopped values !***
		 ichopped_neg = idnint(aint((a1_neg + a2_neg+a3_neg) /
     2							SINH_MIN))
		 ichopped_pos = idnint(aint((a2_pos + a2_pos+a3_pos) /
     2							SINH_MAX))
		 ichopped = ichopped
     2        + idnint(am_zero) * (ichopped_neg + ichopped_pos)

	  end do 

	  !*** info !***
	  if (ichopped .gt. 0) then
	  call vnmpri(2,'% C_VEC: trapped exp overflows:          ',
     2        41, ichopped)
	  endif
	  
      !*** end it !***
      return
      end
	  
      subroutine dc_vecsmpbe(coef,uin,uout,nx,ny,nz,ipkey)
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
	  
c* Added by DG SMPBE variables and common blocks
      double precision Na, fracOccA, fracOccB, fracOccC, phi, ionStr
      double precision z1, z2, z3, ca, cb, cc, a, k
      double precision a1_neg, a1_pos, a2_neg, a2_pos
      double precision a3_neg, a3_pos, a1, a2, a3
      double precision f, g, fprime, gprime, gpark, alpha
      parameter        (Na = 6.022045000e-04)

      double precision v1, v2, v3, conc1, conc2, conc3, vol, relSize
      common /MYPDEF/  v1, v2, v3, conc1, conc2, conc3, vol, relSize
	  
      common /MYPDEF/ charge
      common /MYPDEF/ sconc
      common /MYPDEF/ nion

c*
c*    *** find parallel loops (ipara), remainder (ivect) ***
      n     = nx * ny * nz
      ipara = n / nproc
      ivect = mod(n,nproc)

      do i = 1, n
        uout(i) = 0.0
	  end do
	  
c*    *** initialize the chopped counter ***
	  ichopped = 0
	  
	  z1 = v1
	  z2 = v2
	  z3 = v3
	  ca = conc1
	  cb = conc2
	  cc = conc3
	  a = vol
	  k = relSize

	  if ((k-1).lt. ZSMALL) then
	  call vnmprt(2, '!! DC_VEC: k=1, using special routine', 37) 
	  endif

c* derived quantities
	  fracOccA = Na*ca*a**3
	  fracOccB = Na*cb*a**3
	  fracOccC = Na*cc*a**3
	  phi = (fracOccA/k) + fracOccB + fracOccC
	  alpha = (fracOccA/k)/(1-phi)
	  ionStr = 0.5*(ca*z1**2 + cb*z2**2 + cc*z3**2)

	  do i = 1, n

		 am_zero = dmin1(ZSMALL,dabs(coef(i)))*ZLARGE

c*       compute the arguments for exp(-z*u) term
		 a1_neg = dmax1( dmin1(-1.0*z1*uin(i),0.0d0) , SINH_MIN)
		 a1_pos = dmin1( dmax1(-1.0*z1*uin(i),0.0d0) , SINH_MAX)

c*       compute the arguments for exp(-u) term
		 a2_neg = dmax1( dmin1(-1.0*z2*uin(i),0.0d0) , SINH_MIN)
		 a2_pos = dmin1( dmax1(-1.0*z2*uin(i),0.0d0) , SINH_MAX)

c*       compute the arguments for exp(u) term
		 a3_neg = dmax1( dmin1(-1.0*z3*uin(i),0.0d0) , SINH_MIN)
		 a3_pos = dmin1( dmax1(-1.0*z3*uin(i),0.0d0) , SINH_MAX)

		 a1 = am_zero * (a1_neg + a1_pos)
		 a2 = am_zero * (a2_neg + a2_pos)
		 a3 = am_zero * (a3_neg + a3_pos)

		 gpark = ((1+alpha*exp(a1))/(1+alpha))

		 if ((k-1).lt. ZSMALL) then
		 f = z1*ca*exp(a1) + z2*cb*exp(a2) + z3*cc*exp(a3)
		 g = 1-phi+fracOccA*exp(a1) + fracOccB*exp(a2)
     2        + fracOccC*exp(a3)
		 fprime = -z1**2*ca*exp(a1) - z2**2*cb*exp(a2)
     2        - z3**2*cc*exp(a3)
		 gprime = -z1*fracOccA*exp(a1) - z2*fracOccB*exp(a2)
     2        - z3*fracOccC*exp(a3)
		 else
		 f = z1*ca*exp(a1)*gpark**(k-1)+z2*cb*exp(a2)+z3*cc*exp(a3)
		 g = (1-phi+(fracOccA/k))*gpark**k 
     2        + fracOccB*exp(a2)+fracOccC*exp(a3)
		 fprime = -z1**2*ca*exp(a1)*gpark**(k-2)*(gpark 
     2        + (k-1)*(alpha/(1+alpha))*exp(a1)) - z2**2*cb*exp(a2)
     3        - z3**2*cc*exp(a3)
		 gprime = -k*z1*(alpha/(1+alpha))*exp(a1)
     2        *(1-phi+(fracOccA/k))*gpark**(k-1)-z2*fracOccB*exp(a2)
     3		  -z3*fracOccC*exp(a3)

		 endif

		 uout(i) = -1.0*coef(i)*(0.5/ionStr)*
     2				(fprime*g - gprime*f)/g**2

c*       *** count chopped values ***
		 ichopped_neg = idnint(aint((a1_neg + a2_neg+a3_neg) 
     2									/ SINH_MIN))
		 ichopped_pos = idnint(aint((a2_pos + a2_pos+a3_pos) 
     2									/ SINH_MAX))
		 ichopped = ichopped
     2        + idnint(am_zero) * (ichopped_neg + ichopped_pos)

	  end do

c*    *** info ***
	  if (ichopped .gt. 0) then
	  call vnmpri(2,'% DC_VEC: trapped exp overflows:          ',
     2        42, ichopped)
	  endif
	  
      return
      end
