!**********************************************
! coded by Zhan Chen 
! This main program is for the differential geometry in eulerian formulation
! Michigan state university 
!**********************************************

!************this module is for the global computational domain informatio
 MODULE COMDATA

   character(100) :: fname
   integer :: nx,ny,nz
   real*8 :: xleft, xright, yleft, yright, zleft, zright, DELTAX,DELTAY,DELTAZ, dcel
   real*8, dimension(:), allocatable :: xc, yc, zc
   real*8 :: pi

 END MODULE 

 !********this module is only used for surface generation parameter***. 
 MODULE LJPARAMETERS 

   REAL*8 TAUVAL                 ! the tau value to build a computational region
   REAL*8 PROB 
   REAL*8:: VDWDISPERSION,SIGMAS, RORO,CONMS,DENSITY ! VDWDISPERSION previously called REPULSIVE 
   INTEGER,PARAMETER:: IOSETAR=1, IOSETAA=1,IWCA=1	! D.G.T. changed value of IOSETAR to 1 on 3/8/2011 
   INTEGER :: FFMODEL			! added by D.G.T. on 8/1/2011
   REAL*8 :: EPSILONW
 END MODULE 
     
!*************************MAIN PROGRAM*************************
program mib_3d_chen
USE COMDATA
USE LJPARAMETERS

    Implicit Real*8(A-H, O-Z)
	INTEGER, PARAMETER:: MAXATOM=15000                            ! pre-set maximum number of atom
    REAL*8:: XYZR(4,MAXATOM),PQR(MAXATOM),POS(3)                 ! protein information
	REAL*8:: LJEPSILON(MAXATOM)		! line added by D.G.T. on 7/25/2011
	REAL*8, DIMENSION(:,:,:),ALLOCATABLE::PHI,PHIX,SURFU,EPS,PHIVOC     ! potential and surface level set value 
	REAL*8, DIMENSION(:,:,:),ALLOCATABLE:: phi_c,delphix_c,delphiy_c,delphiz_c! coulombic potential and its gradient; added by DGT 4/25/2013
	REAL*8, DIMENSION(:,:,:),ALLOCATABLE:: surfux,surfuy,surfuz ! gradient components of S; added by DGT 4/25/2013
	REAL*8, DIMENSION(:,:,:),ALLOCATABLE:: phi_rf	! reaction field potential; added by DGT 4/25/2013
	CHARACTER(100)::atomname,fnames, FILNAME(MAXATOM)  ! for the group of protein
	INTEGER :: NMOL                                              ! number of proteins in a group 
	REAL*8,DIMENSION(:),ALLOCATABLE:: SOLV(:),BG(:),BGUNI(:)               ! iteration solvation
	REAL*8, ALLOCATABLE:: corlocqt(:,:,:),charget(:,:)
	INTEGER,ALLOCATABLE:: loc_qt(:,:,:)  
	REAL*8, DIMENSION(MAXATOM)::ELEC,SAS,SV,SVINT,SUMPOT               ! targeted quatities
	REAL*8:: EXPV(100),EXPVALUE,RADIUS                                 ! experiment data 
	INTEGER:: NATOM, NCHR
	INTEGER:: LENFNAME,IMORD
	REAL*8 :: pqrvalue	! Added by D.G.T. on 1/21/2011
	REAL*8 :: epsilon	! Added by D.G.T. on 7/25/2011
	!****indicator and parameter
	INTEGER:: IOUTPUTP, IOUTPUTS  
	INTEGER:: IDACSL
	REAL*8 :: EPSILONS,EPSILONP,SOLENG(2),ALPHA
	REAL*8 :: ERR, MAXERR ,WEIT    
	INTEGER:: IPREC,IPATH,ISTEP, IADI 
	INTEGER:: IPBIN, IGFIN 
    REAL*8 :: FSURF, FPB,TSURF,TPB, TBEGIN,TEND
	         ! time span calculation for first iteration (fsurf and fpb) and total (tsurf, tpb)
	INTEGER:: ITERF,ITERT           ! iteration number for first iteration and total   
	INTEGER:: IDEFAULT      
	REAL*8 :: GAMA,GAMA_I,PRES_I,GAMA_STEP,PRES_STEP		! added by D.G.T. on 3/7/2011 (Monday)
							! GAMA_STEP added by D.G.T. on 3/16/2011
	REAL*8 :: EXTVALUE 	! added by D.G.T. on 8/11/2011
	REAL*8 :: PRES		! added by D.G.T. on 3/11/2011 (Friday)
	REAL*8 :: DIFFENERGY	! added by D.G.T. on 5/11/2011 (Wed)
	INTEGER :: npiter,indpres,ngiter, indgama	! added by D.G.T. on 3/8/2011	
	CHARACTER*2 :: filenumstr
	INTEGER :: IFLAG_REGULARIZED    !added by D.G.T. on 4/26/2013
	
	real:: start,finish	
	real total_time

!*********main body*************************************************************

	call cpu_time(start)

!	WRITE(*,*)'please enter the txt file name containing names of compounds'
!    READ(*,*)fnames  	
	fnames="17set"

!     FNAMES="complex"                ! file name 
    	read(12,*) NMOL		! no. of molecule files
	read(12,*) PRES_I
        read(12,*) GAMA_I
        read(12,*) npiter
        read(12,*) ngiter
	read(12,*) TAUVAL
	read(12,*) PROB
        read(12,*) FFMODEL
        read(12,*) SIGMAS
        read(12,*) EPSILONW
        read(12,*) VDWDISPERSION
        read(12,*) EXTVALUE
	read(12,*) IPREC	! flag to indicte the usage of preconditioner iprec=1(yes) iprec =0(no)
        read(12,*) istep
 	read(12,*) IADI	 ! 0: explicit scheme, 1: ADI scheme 
	read(12,*) ALPHA ! weight of previous solution to change the next solution in geometry flow.
	read(12,*) IPBIN ! enter start guess for PB 1: inherit 0
	read(12,*) TOL
	read(12,*) TOTTF ! total time
	read(12,*) DCEL  
	read(12,*) MAXSTEP ! maximum step of iteration
	read(12,*) EPSILONS
	read(12,*) EPSILONP
	read(12,*) RADEXP
	read(12,*) CREVALUE ! don't forget to change the name of the output files.
	read(12,*) idacsl
    read(12,*) density  ! added by DGT on April 30,2012
    read(12,*) IFLAG_REGULARIZED    ! added by DGT on 4/25/2013
!*********parameters mainly for algorithms********
	 
 	IOUTPUTP=0; IOUTPUTS=0    ! ioutputp: 1 for output potential, ioutputs: 1 for oupt isosurface value

	PRES=PRES_I

	OPEN(5,FILE="press_gama_rms.txt")     ! added by D.G.T. on 3/15/2011
	OPEN(6,FILE="free_energy.txt")		! added by D.G.T. on 6/2/2011
	do IMORD=1,NMOL
	write(filenumstr,'(I2)') IMORD+9
	OPEN(IMORD+30,FILE='output_mol'//trim(filenumstr)//'.txt')
	WRITE(IMORD+30,'(4X,A,4X,A,4X,A,4X,A,5X,A,2X,A,4X,A,4X,A,4X,A,4X,A)') "pressure","gamma","nonpolar","electro ","totl_energy", "area","volume","uattint","expt","error"
	enddo
	
	 WRITE(5,'(4X,A,4X,A,4X,A,4X,A)') "pressure","gamma","rmse","rmse_error"

    if(IFLAG_REGULARIZED.eq.0) then
    write(*,*) "Using original formulation of the generalized Poisson equation"
    else if (IFLAG_REGULARIZED.eq.1) then
    write(*,*) "Using the regularized formulation of the generalized Poisson equation"
    endif

	do indpres=1,npiter

		if((PRES.lt.0.001d0)) then
	PRES_STEP = 0.0001d0
		endif
		if((PRES.ge.0.001d0).and.(PRES.lt.0.01d0)) then
	PRES_STEP = 0.001d0
		endif
	if(PRES.ge.0.01d0) then
        PRES_STEP = 0.005d0
                endif
	
	if(indpres.gt.1) then
	PRES = PRES+PRES_STEP
	endif

	GAMA=GAMA_I
	do indgama = 1,ngiter	! added by D.G.T. on 3/8/2011: do loop to repeat calculations with different fitted
!				values of gama

	if ((GAMA.ge.0.00001d0).and.(GAMA.lt.0.0001d0)) then
        GAMA_STEP = 0.00001d0
        endif

	if ((GAMA.ge.0.0001d0).and.(GAMA.lt.0.001d0)) then
	GAMA_STEP = 0.0001d0
	endif
	
	if ((GAMA.ge.0.001d0).and.(GAMA.lt.0.01d0)) then
        GAMA_STEP = 0.001d0
        endif
	
	if ((GAMA.ge.0.01d0).and.(GAMA.le.0.055)) then
        GAMA_STEP = 0.005d0
        endif

	if (GAMA.gt.0.055) then
	GAMA_STEP = 0.005d0
	endif
	
	write(*,*) 'IFLAG_REGULARIZED = ', IFLAG_REGULARIZED
	write(*,*) 'EPSILONP = ', EPSILONP
	write(*,*) 'DCEL = ', DCEL	
	write(*,*) 'GAMA_STEP= ',GAMA_STEP

	if(indgama.gt.1) then
	GAMA=GAMA+GAMA_STEP
	endif
	
!	RORO= 0.03346d0/GAMA	! added by D.G.T. on 3/8/2011		
        RORO = density/GAMA     ! added by D.G.T. on 4/30/2012 (Monday)  
 	POTCOE = 1.d0/GAMA	! added by D.G.T. on 3/7/2011 (Monday)
	CONMS=PRES/GAMA		! added by D.G.T. on 3/11/2011 (Friday)

	write(*,*) 'GAMA = ', GAMA
        write(*,*) 'density = ', density !added by DGT on 4/30/2012
	write(*,*) 'RORO= ', RORO
	write(*,*) 'POTCOE= ',POTCOE
	write(*,*) 'PRES= ', PRES
	write(*,*) 'CONMS= ', CONMS
	

	NMOL=0   

! **************************************
!    CALL CPU_TIME(TBEGIN)
    IGFIN=1
    FSURF=0.d0; FPB=0.d0;TSURF=0.d0;TPB=0.d0 
    ITERF=0; ITERT=0
	
	deltax=dcel
	deltay=dcel
	deltaz=dcel
 	PI=DACOS(-1.D0)  
 
 ! **** input the name of file which contains a group of system***
    
	LENFNAME = LEN(FNAMES) 
    DO WHILE (FNAMES(LENFNAME:LENFNAME) .EQ. ' ')
        LENFNAME = LENFNAME - 1
    END DO

 !*****************read the file
	OPEN(1,file=FNAMES(1:lenfname)//".txt")            
    NMOL = 0
    DO 
        READ(1,*,IOSTAT = MEOF) ATOMNAME,EXPVALUE
        IF(MEOF .LT. 0) EXIT
            NMOL = NMOL+ 1
			FILNAME(NMOL)=ATOMNAME
            EXPV(NMOL)=EXPVALUE	  
    END DO         
    CLOSE(1)

!**************loop at the level of system begin here***********

   DO IMORD=1,NMOL    ! to obtain the electrostatic energy
       
       FNAME=FILNAME(IMORD)                          
!**********obtain the information of protein
	   LENFNAME = LEN(FNAME) 
       DO WHILE (FNAME(LENFNAME:LENFNAME) .EQ. ' ')
            LENFNAME = LENFNAME - 1
       END DO
 
		write(*,*) 'molecule coordinate file', FNAME 
!**********read atoms     `
       OPEN(2,file=FNAME(1:lenfname)//".xyzr")            
       NATM = 0
  	pqr = 0.d0 	! Added this line by D.G.T. on 1/21/2011
	
	SELECT CASE(FFMODEL)
	
	CASE(1)
       DO 	
           READ(2,*,IOSTAT = MEOF)POS(1), POS(2), POS(3), RADIUS, pqrvalue  ! added ", pqrvalue" by D.G.T. on 1/21/2011
           IF(MEOF .LT. 0) EXIT
           NATM = NATM + 1
           XYZR(1:3,NATM) = POS;  XYZR(4,NATM) = RADIUS
	   pqr(NATM) = pqrvalue	! this line added by D.G.T. on 1/21/2011
	       
		
       END DO         
       CLOSE(2)

	CASE(2)
       DO
         READ(2,*,IOSTAT = MEOF)POS(1), POS(2), POS(3), RADIUS, pqrvalue,epsilon  ! added ", pqrvalue" by D.G.T. on 1/21/2011
           IF(MEOF .LT. 0) EXIT
           NATM = NATM + 1
           XYZR(1:3,NATM) = POS;  XYZR(4,NATM) = RADIUS
           pqr(NATM) = pqrvalue ! this line added by D.G.T. on 1/21/2011
           LJEPSILON(NATM)=epsilon      ! this line added by D.G.T. on 7/24/2011
       END DO
       CLOSE(2)

	END SELECT

	   XYZR(4,:) = XYZR(4,:)*RADEXP	

       CALL DOMAININI(XYZR,NATM,EXTVALUE)          ! initialize the domain   
	   
	   nchr=natm                          ! Charge distribution. 
       allocate(corlocqt(nchr,8,3),charget(nchr,8),loc_qt(nchr,8,3))
	   corlocqt=0.d0;charget=0.d0;loc_qt=0
       do iatm=1,nchr
	       call chargedist(xyzr,pqr,nchr,charget,corlocqt,loc_qt,iatm)    
       enddo              

       ALLOCATE(XC(NX),YC(NY),ZC(NZ),STAT=ierr)
       DO IX=1,NX
         XC(IX) = XLEFT+(IX-1)*DCEL
       END DO

       DO IY=1,NY
         YC(IY) = YLEFT+(IY-1)*DCEL
       END DO

       DO IZ=1,NZ
          ZC(IZ) = ZLEFT+(IZ-1)*DCEL
       END DO

	   ALLOCATE(SURFU(NX,NY,NZ),PHI(NX,NY,NZ),PHIX(NX,NY,NZ),EPS(NX,NY,NZ),PHIVOC(NX,NY,NZ))
	   ALLOCATE(BG(NX*NY*NZ),BGUNI(NX*NY*NZ)) 
       
	ALLOCATE(phi_c(NX,NY,NZ)) ! added by DGT on 4/25/2013
	ALLOCATE(delphix_c(NX,NY,NZ),delphiy_c(NX,NY,NZ),delphiz_c(NX,NY,NZ))
	ALLOCATE(phi_rf(NX,NY,NZ))	! added by DGT on 4/25/2013
	ALLOCATE(surfux(NX,NY,NZ),surfuy(NX,NY,NZ),surfuz(NX,NY,NZ))    ! added by DGT on 4/25/2013

	SURFU=0.D0
       PHI=0.D0; PHIX=0.D0; PHIVOC=0.D0   

	phi_c = 0.d0	! added by DGT on 4/25/2013
	delphix_c = 0.d0
	delphiy_c = 0.d0
	delphiz_c = 0.d0
	phi_rf = 0.d0	  		! added by DGT " "
    surfux = 0.d0
    surfuy = 0.d0
    surfuz = 0.d0

    CALL computeCoulombPotential(phi_c,EPSILONP,natm,xyzr,pqr)  ! added by DGT on 4/25/2013
	
	! compute the components of the coulombic potential gradient
	! added by DGT on 4/25/2013

	 WEIT=1.D0/(2.D0*DCEL)
           do ix=2,nx-1
              do iy=2,ny-1
                 do iz=2,nz-1
                 delphix_c(ix,iy,iz) = (phi_c(ix+1,iy,iz)-phi_c(ix-1,iy,iz))*WEIT
                 delphiy_c(ix,iy,iz) = (phi_c(ix,iy+1,iz) - phi_c(ix,iy-1,iz))*WEIT
                 delphiz_c(ix,iy,iz) = (phi_c(ix,iy,iz+1) - phi_c(ix,iy,iz-1))*WEIT
                 enddo
              enddo
           enddo



       ILOOP=0
       ALLOCATE(SOLV(MAXSTEP+1))
	SOLV = 0.D0
       DIFFENERGY=100.D0
!	   PRINT*,'begin the self-consistency process'
  ! ***************iteration coupling surface generation and poisson solver**********
       DO WHILE(ILOOP<MAXSTEP.AND.DIFFENERGY>CREVALUE) 
           ILOOP=ILOOP+1             

	       IF(IADI==0)DELTAT=(DCEL)**2/4.5D0 
		   IF(IPATH==0)THEN
		     TOTT=MAX(TOTTF-(ILOOP-1)*1.D0,1.D0) 
		   ELSE 
		      IF(ILOOP==1)THEN
		          TOTT=TOTTF
		      ELSE 
                        IF(IADI==0)THEN
		               TOTT=ISTEP*(DCEL)**2/4.5D0 
                        ELSE 
                             TOTT=ISTEP*(DCEL)**2/4.5D0 
                        ENDIF
		      ENDIF   
		   ENDIF       

	write(*,*) 'ILOOP = ',ILOOP, ' TOTT = ',TOTT         
	        
	    ! surface generation	   	                     
!	       CALL CPU_TIME(TIME1) 
		   AREA=0.D0
		   VOLUME=0.D0
		   ATTINT=0.D0
!		   PRINT*,'Begin surface generation'


	       CALL YHSURFACE(XYZR,LJEPSILON,NATM,TOTT,DELTAT,PHIX,SURFU,ILOOP,AREA,VOLUME,ATTINT,ALPHA,IADI,IGFIN,IOUTPUTS)                   
		
		!EPS
		   EPS=0.D0

		   DO IX=1,NX
		      DO IY=1,NY
			     DO IZ=1,NZ

				    if(SURFU(IX,IY,IZ)>1000.D0) then    ! post-process filtering by JC 02/11/2011
					SURFU(IX,IY,IZ)=1000.D0         ! post-process filtering by JC 02/11/2011
				    elseif(SURFU(IX,IY,IZ)<0.D0) then   ! post-process filtering by JC 02/11/2011
					SURFU(IX,IY,IZ)=0.D0            ! post-process filtering by JC 02/11/2011
				    endif                            ! post-process filtering by JC 02/11/2011

				    EPS(IX,IY,IZ)=EPSILONP+(EPSILONS-EPSILONP)*((1000.D0-SURFU(IX,IY,IZ))/1000.D0)
!		if(IMORD==1) then
!		write(20,'(f10.5,f10.5,f10.5,f10.5)') XC(IX), YC(IY),ZC(IZ),EPS(IX,IY,IZ)
!		ENDIF
				 ENDDO 
			  ENDDO 
		   ENDDO

    ! compute the gradient components of S function
    ! added by DGT on 4/25/2013

     WEIT=1.D0/(2.D0*DCEL)
           do ix=2,nx-1
              do iy=2,ny-1
                 do iz=2,nz-1
                 surfux(ix,iy,iz) = (SURFU(ix+1,iy,iz) - SURFU(ix-1,iy,iz))*WEIT/1000.d0
                 surfuy(ix,iy,iz) = (SURFU(ix,iy+1,iz) - SURFU(ix,iy-1,iz))*WEIT/1000.d0
                 surfuz(ix,iy,iz) = (SURFU(ix,iy,iz+1) - SURFU(ix,iy,iz-1))*WEIT/1000.d0
                 enddo
              enddo
           enddo



		
	    ! PB
	     if(iloop==1)then
    	!	      call seteqb(bg,xyzr,pqr,natm,charget,corlocqt,epsilons)

		  CALL seteqb(bg,xyzr,pqr,natm,charget,corlocqt,epsilons,epsilonp,epsilons,     &
		       surfux,surfuy,surfuz,delphix_c,delphiy_c,delphiz_c,IFLAG_REGULARIZED)

	     endif 
  		
     100   if(idacsl==1)then      ! check pb solver for one ball
	           do i=1,nx
			        do j=1,ny
					    do k=1,nz
						   x=xvalue(i)
                           y=yvalue(j)
                           z=zvalue(k)
						   eps(i,j,k)=(x+y+z)
						   ijk=(i-1)*nz*ny+(j-1)*nz+k
						   if(i==1.or.i==nx.or.j==1.or.j==ny.or.k==1.or.k==nz)then
 								bg(ijk)=cos(x)*cos(y)*cos(z)
						   else 
						        bg(ijk)=(-sin(x)*cos(y)*cos(z)-cos(x)*sin(y)*cos(z)-cos(x)*cos(y)*sin(z))&
								        -3.d0*(x+y+z)*cos(x)*cos(y)*cos(z)
						   endif 
						enddo 
				    enddo 
			    enddo 
	       endif 
           
!		   PRINT*,'Begin solving the PBs'
            if(IFLAG_REGULARIZED.eq.0) then
		   CALL PBSOLVER(EPS,PHI,BG,NX,NY,NZ,DCEL,IPREC,TOL,IPBIN,ITER,TITER)

		   else if(IFLAG_REGULARIZED.eq.1) then
		    CALL PBSOLVER(EPS,PHI_RF,BG,NX,NY,NZ,DCEL,IPREC,TOL,IPBIN,ITER,TITER)
		    PHI = PHI_C + PHI_RF
		   endif

		   IF(ILOOP==1)THEN 
		      FPB=TITER
			  ITERF=ITER
		   ENDIF 
		   TPB=TPB+TITER
		   ITERT=ITERT+ITER
		   

		     
		   if(idacsl==1)then               ! for check the pb solver
		      maxerr=0.d0
			  do i=1,nx
			        do j=1,ny
					    do k=1,nz
						   x=xvalue(i)
                           y=yvalue(j)
                           z=zvalue(k)
						   err=abs(phi(i,j,k)-cos(x)*cos(y)*cos(z))
						   maxerr=max(maxerr,err)
						enddo 
					enddo 
			   enddo  
!			   print*,'the maximum err=',maxerr
			   goto 200 
		    
		   endif  
	      
          IF(ILOOP==1)THEN

           if(IFLAG_REGULARIZED.eq.0) then  ! added by DGT on 4/26/2013
             DO I=1,NX
			     DO J=1,NY
				    DO K=1,NZ
                                   IJK=(I-1)*NY*NZ+(J-1)*NZ+K
                       if(i<2.or.i>nx-1.or.j<2.or.j>ny-1.or.k<2.or.k>nz-1)then
			                
							BGUNI(IJK)=BG(IJK)*EPSILONS/EPSILONP
					   ELSE 
                            BGUNI(IJK)=BG(IJK)
					   ENDIF 
			        ENDDO 
				 ENDDO 
			 ENDDO 
	    ! call seteqb(bguni,xyzr,pqr,natm,charget,corlocqt,epsilonp)
	    ! vocuum potential 
		     DO IX=1,NX
		      DO IY=1,NY
			     DO IZ=1,NZ
	    			    EPS(IX,IY,IZ)=EPSILONP
				 ENDDO 
			  ENDDO 
		     ENDDO  
           

	   CALL PBSOLVER(EPS,PHIVOC,BGUNI,NX,NY,NZ,DCEL,IPREC,TOL,IPBIN,ITER,TITER)  ! DGT inserted this line on June 15, 2012


	! DGT changed argument BG to BGUNI in the above line on March 15, 2013

     !       CALL PBSOLVER(EPS,PHIVOC,BGUNI,NX,NY,NZ,DCEL) 
     !        CALL POISOLVE3D(BGUNI,PHIVOC,NX,NY,NZ,XLEFT,XRIGHT,YLEFT,YRIGHT,ZLEFT,ZRIGHT) ! DGT commented this line on June 15, 2012
 
	   else if(IFLAG_REGULARIZED.eq.1) then


               PHIVOC = PHI_C
       endif

	     ENDIF 	 
	    
        
		   WEIT=1.D0/(2.D0*DCEL)
		   DO IX=2,NX-1
		      DO IY=2,NY-1
			     DO IZ=2,NZ-1
				    PHIXX=(PHI(IX+1,IY,IZ)-PHI(IX-1,IY,IZ))*WEIT
					PHIXY=(PHI(IX,IY+1,IZ)-PHI(IX,IY-1,IZ))*WEIT
					PHIXZ=(PHI(IX,IY,IZ+1)-PHI(IX,IY,IZ-1))*WEIT
				    PHIX(IX,IY,IZ)=0.5D0*(EPSILONS-EPSILONP)*(PHIXX**2+PHIXY**2+PHIXZ**2)*POTCOE  !!?
				 ENDDO 
			  ENDDO 
		   ENDDO  

        !SOLVATION
		   SOLENG(1)=0.D0

           do i=1,nchr
		       do j=1,8
		           i1=loc_qt(i,j,1)
			       j1=loc_qt(i,j,2)
			       k1=loc_qt(i,j,3)
		           ijk=(i1-1)*ny*nz+(j1-1)*nz+k1
		           if(charget(i,j).ne.0.d0)then
		               soleng(1)=soleng(1)+0.5d0*PHI(i1,j1,k1)*charget(i,j)
			       endif
           !        print*,i1,j1,k1,charget(i,j),phi(i1,j1,k1)
	           enddo
	       enddo


	       soleng(2)=0.d0
           do i=1,nchr
		       do j=1,8
		           i1=loc_qt(i,j,1)
			       j1=loc_qt(i,j,2)
			       k1=loc_qt(i,j,3)
		           ijk=(i1-1)*ny*nz+(j1-1)*nz+k1
		           if(charget(i,j).ne.0.d0)then
		               soleng(2)=soleng(2)+0.5d0*phivoc(i1,j1,k1)*charget(i,j)
		           endif
				   
	           enddo
	       enddo


	    !    SOLV(ILOOP)=(soleng(1)-soleng(2))*332.0716d0  ! DGT commented this line on May 1, 2013
	    SOLV(ILOOP)=(soleng(1)-soleng(2))*332.06364d0       ! DGT changed 332.0716 to 332.06364 on May 1, 2013

	 	ELEC(IMORD) = SOLV(ILOOP)
		
		SOLV(ILOOP)=ELEC(IMORD)+GAMA*(AREA+VOLUME*CONMS+ATTINT*RORO)
		

	       IF(ILOOP>1)DIFFENERGY=DABS(SOLV(ILOOP)-SOLV(ILOOP-1))  !convergence criteria

 	
		ENDDO 

	! output the potential distribution

  !electrostatic energy, surface area, volume and attractive integral 	     
	    SAS(IMORD)=AREA
	    SV(IMORD)=VOLUME
	    SVINT(IMORD)=ATTINT  
	    SUMPOT(IMORD)=(SAS(IMORD)+SV(IMORD)*CONMS+SVINT(IMORD)*RORO) ! modified by JC 03/08/2011
	    call cpu_time(t2)
		deallocate(corlocqt,charget,loc_qt)
		deallocate(surfu,phi,phix,eps,phivoc)
		deallocate(phi_c,delphix_c,delphiy_c,delphiz_c,phi_rf)  ! DGT
		deallocate(surfux,surfuy,surfuz)    ! DGT
        DEALLOCATE(XC,YC,ZC)
		deallocate(bg,solv,bguni)

    ENDDO 

!******************* optimize the surface tension coefficience rho
   IF(idAcSl==0.AND.NMOL.GT.1)THEN   

!	OPEN(5,FILE="17setwca.dat")	! added by D.G.T. on 3/8/2011
   write(6,*) "PRESSURE = ", PRES
   write(6,*) "GAMA = ",GAMA

        CALL MAXRMS(SUMPOT,EXPV,ELEC,NMOL,GAMA,ngiter,PRES,SAS,SV,SVINT)	! D.G.T. added GAMA as an argument on 3/7/2011
!	CALL MAXRMS_PGAMA(SUMPOT,EXPV,ELEC,SAS,SV,SVINT,NMOL,GAMA,PRES,ngiter)	!added by D.G.T. on 3/11/2011
!	WRITE(5,*) "GAMA 	RMS"
	CALL WRITERMS_GAMA(SUMPOT,EXPV,ELEC,NMOL,GAMA,ngiter,PRES)	! added by D.G.T. on 3/14/2011
   ELSE
!       PRINT*,'total free energy of solvation is', SUMPOT(1)*0.0065D0+ELEC(1)
!	   PRINT*,'Nonpolar solvation energy', SUMPOT(1)*0.0065D0
	PRINT*,'total free energy of solvation is',SUMPOT(1)*GAMA+ELEC(1)
	   PRINT*,'Nonpolar solvation energy',SUMPOT(1)*GAMA
	   PRINT*,'Electrostatic solvaiton free energy',ELEC(1)
	   PRINT*,'reference system energy (kcal/mol)',soleng(2)*332.06364d0
	    PRINT*,'actual system energy (kcal/mol)',soleng(1)*332.06364d0

   ENDIF 

	enddo	! end of loop "do indgama =1, ngiter": added by D.G.T. on 3/8/2011
	enddo	! end of loop "do indpres = 1, npiter" added by D.G.T. on 3/18/2011


	CLOSE(5)	! added by D.G.T. on 3/8/2011
	CLOSE(6) 	! added by D.G.T. on 6/02/2011
	do IMORD=1,NMOL
	CLOSE(IMORD+30)
	enddo

!***************************************************
200   continue 

	call cpu_time(finish)
	print*, 'Time = ',finish-start, ' seconds.'

end program mib_3d_chen


SUBROUTINE WRITERMS_GAMA(SURF,SEXP,SELEC,NMOL,GAMA,ngiter,PRES)
!***********************************************************************************
! purpose To write gama vs. RMS values
! surf: nonpolar term
! sexp: the experiment value(total)
! selec: electrastatic value
! D.G.T. took part of the code from MAXRMS to create this subroutine on 3/14/2011
!***************************************************************************************

USE COMDATA

  IMPLICIT REAL*8(A-H,O-Z)
  INTEGER:: NMOL,ngiter
  REAL*8::SURF(50),SELEC(50),SEXP(50),SNUM(50),ERR(50)
  REAL*8::SDEN,SUMN
  REAL*8:: GAMA,RMS,PRES,AREA,VOLUME,ATTINT
  REAL*8:: RMS_STDV 
	RMS_STDV=0.D0
        RMS=0.D0
        DO I=1, NMOL          ! obtain the rms

       SNUM(I)=(GAMA)*SURF(I)+SELEC(I)
           ERR(I)=(SNUM(I)-SEXP(I))
!           PRINT*,GAMA*SURF(I),SELEC(I),SNUM(I),ERR(I)
       RMS=RMS+ERR(I)**2

        ENDDO
	
	RMS=DSQRT(RMS/NMOL)
	DO I=1,NMOL
	RMS_STDV = RMS_STDV+(DSQRT(ERR(I)**2)-RMS)**2	
	ENDDO

	RMS_STDV = DSQRT(RMS_STDV/NMOL)

        WRITE(5,'(4F12.7)') PRES,GAMA,RMS,RMS_STDV

END SUBROUTINE

!*******************************************************************

SUBROUTINE MAXRMS(SURF,SEXP,SELEC,NMOL,GAMA,ngiter,PRES,SAS,SV,SVINT) 
!***********************************************************************************
! purpose By knowing the surface area, electrastical value and experiment value,
!         we try to find a optimization value for surface tension coefficience
!         Gnp=Gamma*surf and Sexp is approximated by Gnp plus Selec
! surf: nonpolar term
! sexp: the experiment value(total) 
! selec: electrastatic value 
!***************************************************************************************      

!	D.G.T. included GAMA as an argument on 3/7/2011
!	D.G.T. created new integer ngiter on 3/8/2011

USE COMDATA 

  IMPLICIT REAL*8(A-H,O-Z)
  INTEGER:: NMOL,ngiter
  REAL*8::SURF(50),SELEC(50),SEXP(50),SNUM(50),ERR(50),SAS(50),SV(50),SVINT(50)
  REAL*8::SDEN,SUMN
  REAL*8:: GAMA,RMS,PRES
 
  
!	D.G.T. moved the RMS calculation before calculating the optimized gamma value

        RMS=0.D0
        DO I=1, NMOL          ! obtain the rms

       SNUM(I)=(GAMA)*SURF(I)+SELEC(I)
           ERR(I)=(SNUM(I)-SEXP(I))
       RMS=RMS+ERR(I)**2

        ENDDO
        RMS=DSQRT(RMS/NMOL)

       WRITE(6,'(4X,A,4X,A,5X,A,2X,A,4X,A,4X,A,4X,A,4X,A)') "nonpolar","electro ","totl_energy", "area","volume","uattint","expt","error"
        DO I=1,NMOL
       WRITE(6,'(9F13.7)') GAMA*SURF(I),SELEC(I),SNUM(I),SAS(I),SV(I),SVINT(I),SEXP(I),ERR(I)

	WRITE(I+30,'(10F13.7)') PRES,GAMA,GAMA*SURF(I),SELEC(I),SNUM(I),SAS(I),SV(I),SVINT(I),SEXP(I),ERR(I)	
        ENDDO
        WRITE(6,'(A,F12.7)')"the root mean square=",RMS

 END SUBROUTINE 

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
SUBROUTINE MAXRMS_PGAMA(SURF,SEXP,SELEC,SSAS,SSV,SSVINT,NMOL,GAMA,PRES,ngiter) 
USE COMDATA

  IMPLICIT NONE
  INTEGER:: NMOL,ngiter,I
  REAL*8::SURF(50),SELEC(50),SEXP(50),SNUM(50),ERR(50),X(50),Y(50),Z(50)
  REAL*8:: SSAS(50),SSV(50),SSVINT(50)
  REAL*8:: SUM_XX, SUM_YY, SUM_XY, SUM_XZ, SUM_YZ
  REAL*8::SDEN,SUMN
  REAL*8:: GAMA,RMS,PRES
	

        write(5,*) 'Number of iterations, ngiter: ', ngiter
        write(5,'(A,F12.4)') 'current gamma value, GAMA: ', GAMA
	write(5,'(A,F12.4)') 'current pressure value, PRES: ', PRES


 	 RMS=0.D0
        DO I=1, NMOL          ! obtain the rms

       SNUM(I)=(GAMA)*SURF(I)+SELEC(I)
           ERR(I)=(SNUM(I)-SEXP(I))
           PRINT*,GAMA*SURF(I),SELEC(I),SNUM(I),ERR(I)
       RMS=RMS+ERR(I)**2

        ENDDO
        RMS=DSQRT(RMS/NMOL)

	    WRITE(5,'(4X,A,5X,A,2X,A,4X,A)')"nonpolar","electro ","totl_energy","error"
        DO I=1,NMOL
       WRITE(5,'(5F12.4)')GAMA*SURF(I),SELEC(I),SNUM(I),ERR(I)
        ENDDO
        WRITE(5,'(A,F12.4)')"GAMA=",GAMA
	WRITE(5,'(A,F12.4)')"PRES = ", PRES
        WRITE(5,'(A,F12.4)')"the root mean square=",RMS
        WRITE(5,*)

	X = 0.d0
	Y=0.d0	
	Z=0.d0
	
	do I=1,NMOL
	
	X(I)=SSAS(I)
	Y(I)=SSV(I)
	Z(I) = SEXP(I)-SELEC(I)-0.03346d0*SSVINT(I)
	enddo
	

	SUM_XY=0.d0
	SUM_YY=0.d0
	SUM_XY=0.d0
	SUM_XZ=0.d0
	SUM_YZ=0.d0

	do I=1,NMOL
	SUM_XY =  SUM_XX+X(I)*X(I)
	SUM_YY = SUM_YY+Y(I)*Y(I)
	SUM_XY = SUM_XY+X(I)*Y(I)
	SUM_XZ = SUM_XZ +X(I)*Z(I)
	SUM_YZ = SUM_YZ +Y(I)*Z(I)
	enddo
	
	PRES = (SUM_XZ*SUM_XY-SUM_YZ*SUM_XX)/(SUM_XY*SUM_XY-SUM_YY*SUM_XX)
	GAMA = (SUM_XZ-PRES*SUM_XY)/SUM_XX


        PRINT*,"there are test case",nmol
        PRINT*," new gamma value= ",GAMA
	PRINT*,"new pressure value= ", PRES

	WRITE(5,'(A,F12.4)')"New gamma value, GAMA=",GAMA
        WRITE(5,'(A,F12.4)')"New pressure value, PRES= ", PRES

 END SUBROUTINE

