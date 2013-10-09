

!-----$-------$-------$-------$-------$-------$-------$-------$-------$
!      Molecular Surface Generation by pde geometry flow 
!-----$-------$-------$-------$-------$-------$-------$-------$-------$

 SUBROUTINE DOMAININI(XYZR,NATM,EXTVALUE)
 USE COMDATA 
 IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-N)
 REAL*8  :: XYZR(4,NATM), ATOM_X(NATM),ATOM_Y(NATM),ATOM_Z(NATM),ATOM_R(NATM)
 REAL*8  :: EXTVALUE,TMPMIN,TMPMAX
 REAL	:: LEFT,RIGHT     
	 
     DO I=1,NATM
	     
         ATOM_X(I) =XYZR(1,I)
         ATOM_Y(I) =XYZR(2,I)
         ATOM_Z(I) =XYZR(3,I)
         ATOM_R(I) =XYZR(4,I)

	  ENDDO        
!*****************grid domain setting**************************    
!      EXTVALUE=1.9D0	! D.G.T. commented on August 8, 2011

      TMPMIN=ATOM_X(1)-ATOM_R(1) - EXTVALUE   !X-GRID                       !!!!*
      TMPMAX=ATOM_X(1)+ATOM_R(1) + EXTVALUE
      DO I=2,NATM
         
            TMP=ATOM_X(I)+ATOM_R(I) + EXTVALUE
            IF (TMPMAX.LT.TMP) TMPMAX=TMP
            TMP=ATOM_X(I)-ATOM_R(I) - EXTVALUE
            IF (TMPMIN.GT.TMP) TMPMIN=TMP
        
      END DO
	LEFT = TMPMIN/DCEL
	RIGHT = TMPMAX/DCEL
      XLEFT  = FLOOR(LEFT)*DCEL- EXTVALUE
      XRIGHT = CEILING(RIGHT)*DCEL+ EXTVALUE


      TMPMIN=ATOM_Y(1)-ATOM_R(1) - EXTVALUE   !Y-GRID
      TMPMAX=ATOM_Y(1)+ATOM_R(1) + EXTVALUE

      DO I=2,NATM
         
            TMP=ATOM_Y(I)+ATOM_R(I) + EXTVALUE
            IF (TMPMAX.LT.TMP) TMPMAX=TMP
            TMP=ATOM_Y(I)-ATOM_R(I) - EXTVALUE
            IF (TMPMIN.GT.TMP) TMPMIN=TMP
        
      END DO
	LEFT = TMPMIN/DCEL
        RIGHT = TMPMAX/DCEL
      YLEFT  = FLOOR(LEFT)*DCEL- EXTVALUE
      YRIGHT = CEILING(RIGHT)*DCEL + EXTVALUE


      TMPMIN=ATOM_Z(1)-ATOM_R(1) - EXTVALUE    !Z-GRID
      TMPMAX=ATOM_Z(1)+ATOM_R(1) + EXTVALUE
      DO I=2,NATM
         
            TMP=ATOM_Z(I)+ATOM_R(I) + EXTVALUE
            IF (TMPMAX.LT.TMP) TMPMAX=TMP
            TMP=ATOM_Z(I)-ATOM_R(I) - EXTVALUE
            IF (TMPMIN.GT.TMP) TMPMIN=TMP 
         
      END DO
	LEFT = TMPMIN/DCEL
        RIGHT = TMPMAX/DCEL
      ZLEFT  = FLOOR(LEFT)*DCEL- EXTVALUE
      ZRIGHT = CEILING(RIGHT)*DCEL+ EXTVALUE


      NX = (XRIGHT - XLEFT)/DCEL + 1
      NY = (YRIGHT - YLEFT)/DCEL + 1 
      NZ = (ZRIGHT - ZLEFT)/DCEL + 1   

	  XRIGHT=XLEFT+(NX-1)*DCEL
	  YRIGHT=YLEFT+(NY-1)*DCEL
	  ZRIGHT=ZLEFT+(NZ-1)*DCEL

 RETURN 
 END 
           
 SUBROUTINE yhsurface(XYZR,LJEPSILON,NATM,TOTT,DT,PHITOTX,SURFU,ILOOP,AREA,VOLUME,ATTINT,ALPHA,IADI,IGFIN,IOUTPUTS)
 
!****************************************************
! INPUT: xyzr(1:3,:) coordinator, xyzr(4,:) radii 
!        natm: number of atoms in the system
!        tott: total time span for geometry flow
!        Dt:   time step size
!        phitotx: electrostatic potential energy
!        alpha: the weight factor of previous surface location
!        iadi: flag to indicate whether adi is used.
!        ioutputs: flag for whether output the s function
!	 LJEPSILON(NATM): epsilon parameters for different atom types to be used for LJ potential
! OUTPUT: surfu: inout  in: the initial guess. out: the value of level set function 
!         area,volume and attint
!**************************************************
      USE COMDATA       ! grid domain and grid size
      USE LJPARAMETERS  ! only for the lj potential parameters 

      IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-N)
      
      REAL*8  :: SURFU(NX,NY,NZ)
      REAL*8  :: TT,RATIO,RP,TOTT,DT
      REAL*8  :: AREA,VOLUME,ATTINT,SUX,SUY,SUZ,ALPHA
      INTEGER :: I,J,IX,IY,IZ,ILOOP
 
      REAL*8  :: XYZR(4,NATM), ATOM_X(NATM),ATOM_Y(NATM),ATOM_Z(NATM),ATOM_R(NATM)
      REAL*8  :: PHITOTX(NX,NY,NZ)
      REAL*8,  ALLOCATABLE :: SU(:,:,:)
      REAL*8,  ALLOCATABLE :: G(:,:,:)  
            
	  ! LJ potential term
	  REAL*8,  ALLOCATABLE :: POTR(:,:,:),POTA(:,:,:),SIGMA(:),SETA12(:),SETA6(:),EXPAN(:)
	  REAL*8,  ALLOCATABLE :: FINTEGR(:,:,:)
      INTEGER :: NIND(5),NT,IADI,IGFIN

     REAL*8 :: LJEPSILON(NATM),EPSILON(NATM)
     REAL*8 :: RCFACTOR			! entered first on Jan 26, 2012 by DGT 	  
      DO I=1,NATM
	     
         ATOM_X(I) =XYZR(1,I)
         ATOM_Y(I) =XYZR(2,I)
         ATOM_Z(I) =XYZR(3,I)
         ATOM_R(I) =XYZR(4,I)

	  ENDDO        
      
      DDX =DCEL               
      NX2 =NX        !2*NX-1    inherited from previous code resulting in introduction of nx2       
      NY2 =NY        !2*NY - 1
      NZ2 =NZ        !2*NZ - 1
     
 !    WRITE(*,*) 'Resolution for the molecular surface'
 !    WRITE(*,'(3(A10,I4))') 'NX2=',NX2, 'NY2=',NY2, 'NZ2=',NZ2
      LENFNAME = LEN(FNAME) 
      DO WHILE (FNAME(LENFNAME:LENFNAME) .EQ. ' ')
          LENFNAME = LENFNAME - 1
      END DO
      
  !************LJ parameter setting******************************************
       
	  DO I=1,NATM                            ! change the atom radi
	  	 ATOM_R(I)=ATOM_R(I)     !-DCEL/4.d0
	  ENDDO 
	  
      ALLOCATE(SU(NX2,NY2,NZ2),STAT=ierr)
      ALLOCATE(G(NX2,NY2,NZ2),STAT=ierr)
      CALL INITIAL(NX2,NY2,NZ2,XLEFT,YLEFT,ZLEFT,DDX,NATM, &
                   ATOM_X,ATOM_Y,ATOM_Z,ATOM_R,G,SU)
      
!  *********impose the previous solution
      IF(ILOOP>1.AND.IGFIN==1)THEN
	      SU=SURFU
	  ENDIF 
!*********************************************
      ALLOCATE(SIGMA(NATM),SETA12(NATM),SETA6(NATM),EXPAN(NATM))  
	
	SELECT CASE(FFMODEL)

	CASE(1)
!       for ZAP-9/AM1-BCCv

!	write(*,*) "surfconcz.f90/yhsurface: model used is ZAP-9/AM1-BCCv"
	RCFACTOR = 1.0
      DO I=1,NATM 
     	                            ! set the value for sigma and seta
         SIGMA(I)=ATOM_R(I)         !!*2.D0**(1.D0/6.D0)   
         EXPAN(I)=ATOM_R(I)+PROB
	     SIGMA(I)=SIGMA(I)+SIGMAS
!	two lines below are part of the original code sent by Zhan. D.G.T commented them on March 7, 2011

!	     SETA12(I)=IOSETAR*VDWDISPERSION/((SIGMA(I)/EXPAN(I))**12-2.D0*(SIGMA(I)/EXPAN(I))**6)/RORO
!	     SETA6(I)=IOSETAA*VDWDISPERSION/((SIGMA(I)/EXPAN(I))**12-2.D0*(SIGMA(I)/EXPAN(I))**6)/RORO                  

!	two lines below added by D.G.T. on March 7, 2011 (Monday)

	
	EPSILON(I) = 1.D0/((SIGMA(I)/EXPAN(I))**12-2.D0*(SIGMA(I)/EXPAN(I))**6)

	if(VDWDISPERSION==0) then
	EPSILON(I) = 0.D0
	endif

	SETA12(I)=IOSETAR*VDWDISPERSION*EPSILON(I)
        SETA6(I)=2.D0*IOSETAA*VDWDISPERSION*EPSILON(I)

      ENDDO 
      
	CASE(2)
!	write(*,*) "surfconcz.f90/yhsurface: model used is OPLS/AA"
	RCFACTOR = 2.D0**(1.D0/6.D0)
	    DO I=1,NATM
         SIGMA(I)=DSQRT(2.D0*ATOM_R(I)*2*SIGMAS)      
	EPSILON(I) = DSQRT(LJEPSILON(I)*EPSILONW)

	if(VDWDISPERSION==0) then
	EPSILON(I) = 0.D0
	endif

	SETA12(I) = 4.D0*EPSILON(I)
	SETA6(I) = 4.D0*EPSILON(I)
	ENDDO
	
	END SELECT

	

	  ALLOCATE(POTR(NX2,NY2,NZ2),POTA(NX2,NY2,NZ2),FINTEGR(NX,NY,NZ),STAT=ierr)
      POTR=0.D0;POTA=0.D0        ! potr: repulsive term   pota: attractive term 
      DO IX=2,NX2-1
         DO IY=2,NY2-1
	        DO IZ=2,NZ2-1

		        IF(G(IX,IY,IZ).NE.0.d0)THEN 

                   XI=XLEFT+(IX-1)*DDX
                   YI=YLEFT+(IY-1)*DDX
		           ZI=ZLEFT+(IZ-1)*DDX
		           PR=0.D0
		           PA=0.D0 
		           DO IATO=1,NATM   
				                              !!!***
			          DIS=SQRT((XI-ATOM_X(IATO))**2+(YI-ATOM_Y(IATO))**2+(ZI-ATOM_Z(IATO))**2)+PROB
                              IF(DIS==0)THEN
                                  RATIO=1
                              ELSE 
             			       RATIO=SIGMA(IATO)/DIS
                              ENDIF
				      IF(IWCA==1)THEN
					     IF((RATIO*RCFACTOR).GT.1)THEN
				            PR=PR+SETA12(IATO)*(RATIO)**12.D0+EPSILON(IATO)-SETA6(IATO)*(RATIO)**6.D0
		                    PA=PA-EPSILON(IATO)   
				         ELSE 
				            PA=PA-SETA6(IATO)*(RATIO)**6.D0+SETA12(IATO)*(RATIO)**12.D0    ! this has more advanced expression for lj (wca) 
                         ENDIF 
			           ELSE 
				            PR=PR+SETA12(IATO)*(RATIO)**12.D0
					        PA=PA-SETA6(IATO)*(RATIO)**6.D0
			           ENDIF 

			        ENDDO 
   			    	POTR(IX,IY,IZ)=PR
			        POTA(IX,IY,IZ)=PA
		         ENDIF
             ENDDO 
          ENDDO 
      ENDDO 

  ! ******************end of lj potential term************************

      CALL CPU_TIME(TIME11)     

      IF(IWCA==1)POTR=0.D0                    !!?
	  DO I=1,NX                   ! update phitotx to be total potential
	     DO J=1,NY
		     DO K=1,NZ
                       PHITOTX(I,J,K)=-CONMS-PHITOTX(I,J,K)+RORO*POTR(I,J,K)+RORO*POTA(I,J,K)
	            ENDDO
	     ENDDO 
	  ENDDO 
    
!	write(*,*) "yhsurface: IADI = ", IADI 
       
      IF(IADI==0.OR.ILOOP>1)THEN  ! explicit scheme
               DT=DDX**2/4.5D0
               NT   = ANINT(TOTT/DT)+1   
	        CALL UPWINDING(NX2,NY2,NZ2,DDX,DDX**2/4.5D0,NT,G,SU,PHITOTX) 
!	if (ILOOP ==2) then
!	write(*,*) "yhsurface: ILOOP = ", ILOOP
!	 write(*,*) "NX = ", NX
!       write(*,*) "NY = ", NY
!       write(*,*) "NZ = ", NZ
!               DO IX=1,NX
!                       DO IY=1,NY
!                               DO IZ=1,NZ
!               write(50,'(6F12.4)') IX, IY, IZ, G(IX,IY,IZ),SU(IX,IY,IZ), POTA(IX,IY,IZ)
!                               ENDDO
!                       ENDDO
!               ENDDO
!
!	endif

	  ELSE             ! adi scheme 
	
            NT= ANINT(TOTT/DT)+1   
	 
	        DO I=1,NT
		       CALL ADICOR(SU,G,DDX,DT,NZ2,NY2,NX2,PHITOTX) 
	        ENDDO 		

	  ENDIF 

	 
!************* output isosurface function value*******************************
!     IF(ILOOP==1.AND.IOUTPUTS==1)THEN
!         OPEN(1,FILE=FNAME(1:lenfname)//".grid")
!         WRITE(1,'(6F12.4,3I6)') XLEFT,XRIGHT,YLEFT,YRIGHT,ZLEFT,ZRIGHT,NX2,NY2,NZ2
!         CLOSE(1)
     
!         OPEN(3,FILE=FNAME(1:lenfname)//".dat")
!         DO IX=1,NX2
!             DO IY=1,NY2
!                DO IZ=1,NZ2
!                   WRITE(3,*) SU(IX,IY,IZ)
!                END DO
!             END DO
!         END DO
!         CLOSE(3)
!      ENDIF 

      IF(ILOOP>1)THEN

          DO I=1,NX
	          DO J=1,NY
		           DO K=1,NZ

                      SURFU(I,J,K)=SURFU(I,J,K)*ALPHA+SU(I,J,K)*(1.D0-ALPHA)
	                  SU(I,J,K)=SURFU(I,J,K)

	               ENDDO
		      ENDDO 
	       ENDDO 
	  ELSE 
	       SURFU=SU
	  ENDIF 
   
	  CALL CPU_TIME(TIME21)	  
!*****************************************************************	  
! calculate the volum enclosed inside the surface and the volum 
! integration upon outside the surface. 
!*****************************************************************
	  ! VOLUME      
	  VOLUME=0.D0
	  FINTEGR=SU
      CALL VOLUMINTEGRATION(FINTEGR,NX,NY,NZ,DCEL,VOLUME)

      !AREA 
	  AREA=0.D0   
	  WEITGHT=1.D0/(2.D0*DCEL)
	  FINTEGR=0.D0
	  DO IX=2,NX-1
	     DO IY=2,NY-1
		    DO IZ=2,NZ-1
			       
			   	    SUX=(SU(IX+1,IY,IZ)-SU(IX-1,IY,IZ))*WEITGHT
					SUY=(SU(IX,IY+1,IZ)-SU(IX,IY-1,IZ))*WEITGHT	
					SUZ=(SU(IX,IY,IZ+1)-SU(IX,IY,IZ-1))*WEITGHT
					FINTEGR(IX,IY,IZ)=SQRT(SUX**2+SUY**2+SUZ**2)
			  
			ENDDO 
		 ENDDO 
	  ENDDO   
      CALL VOLUMINTEGRATION(FINTEGR,NX,NY,NZ,DCEL,AREA)
	  
	  ! attractive integral      
	  DO IX=2,NX2-1       ! add the potential value outside sas 
         DO IY=2,NY2-1
	         DO IZ=2,NZ2-1

!			        IF(G(IX,IY,IZ).EQ.0.d0.AND.SU(IX,IY,IZ).EQ.1.0d3)THEN
!        	IF(G(IX,IY,IZ).EQ.0.d0.AND.SU(IX,IY,IZ).EQ.0)THEN               ! DGT changed 1.0d3 to 0 on Nov. 28,2011
                IF(G(IX,IY,IZ).NE.0.d0) THEN  ! DGT commented out the above line and add this line on April 30,2012
                      XI=XLEFT+(IX-1)*DDX
                      YI=YLEFT+(IY-1)*DDX
		              ZI=ZLEFT+(IZ-1)*DDX
		              PA=0.D0
					  PR=0.D0
		              DO IATO=1,NATM                              !!!***
			              DIS=SQRT((XI-ATOM_X(IATO))**2+(YI-ATOM_Y(IATO))**2+(ZI-ATOM_Z(IATO))**2)+PROB
		   IF (DIS==0)THEN
            RATIO = 1
            ELSE
                   RATIO=SIGMA(IATO)/DIS
		ENDIF

                   IF(IWCA==1)THEN
					          IF((RATIO*RCFACTOR).GT.1)THEN
	                PR=PR+SETA12(IATO)*(RATIO)**12.D0+EPSILON(IATO)-SETA6(IATO)*(RATIO)**6.D0
		                         PA=PA-EPSILON(IATO)   
				              ELSE 
				                 PA=PA-SETA6(IATO)*(RATIO)**6.D0+SETA12(IATO)*(RATIO)**12.D0    ! this has more advanced expression for lj (wca) 
                             ENDIF 
			              ELSE 
				             PR=PR+SETA12(IATO)*(RATIO)**12.D0
					         PA=PA-SETA6(IATO)*(RATIO)**6.D0
			              ENDIF 
					   ENDDO 
				   POTR(IX,IY,IZ)=PR		
			           POTA(IX,IY,IZ)=PA			

                 ENDIF 
              ENDDO 
		   ENDDO 
 	    ENDDO 
			           
		  
	   DO I=1,NX
	     DO J=1,NY
		    DO K=1,NZ
			   IF(IWCA==1)THEN
		           FINTEGR(I,J,K)=POTA(I,J,K)*(1000.D0-SU(I,J,K))    
               ELSE 
			       FINTEGR(I,J,K)=(POTA(I,J,K) +POTR(I,J,K))*(1000.D0-SU(I,J,K)) 
			   ENDIF 
			ENDDO 
		 ENDDO 
	   ENDDO 
      
	  ATTINT=0.D0
	  CALL VOLUMINTEGRATION(FINTEGR,NX,NY,NZ,DCEL,ATTINT)

!*******************************************************************************************	                                               
      deallocate(SU,STAT=ierr) 
	  deallocate(G,STAT=ierr)
	  deallocate(FINTEGR,EXPAN,STAT=ierr)
	  deallocate(potr,pota,sigma,seta12,seta6,STAT=ierr)
         
END SUBROUTINE YHSURFACE

!------------------------------------------------------------------------------
     SUBROUTINE VOLUMINTEGRATION(F,NX,NY,NZ,DCEL,VOLUM)
     IMPLICIT REAL*8 (A-H,O-Z)
     INTEGER:: NX, NY ,NZ,INDSUM
	 REAL*8 :: F(NX,NY,NZ)
	 REAL*8 :: DCEL,VOLUM

	 !************************************************
	 ! volum integration with ingetrand f and inside outside indicator io
	 !************************************************
	       
     VOLUM=0.D0
  	 F=F/1000.D0

	 DO IX=1,NX
	    DO IY=1,NY
		   DO IZ=1,NZ                
				 VOLUM=VOLUM+F(IX,IY,IZ)*DCEL**3
		   ENDDO 
		ENDDO 
	 ENDDO 
	 
     RETURN 
     END 

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!	initialization
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      SUBROUTINE INITIAL(NX,NY,NZ,XL,YL,ZL,DX,N_ATOM,ATOM_X, &
                         ATOM_Y,ATOM_Z,ATOM_R,G,PHI)
         
	  USE LJPARAMETERS  ! there is variable tauval

      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER :: N_ATOM    
      DIMENSION G(NX,NY,NZ),PHI(NX,NY,NZ),XI(NX),YI(NY),ZI(NZ)      
      DIMENSION ATOM_X(N_ATOM),ATOM_Y(N_ATOM),ATOM_Z(N_ATOM),ATOM_R(N_ATOM)
      DIMENSION FFF(NX*NY*NZ)    
      REAL*8  :: RP  
      REAL  :: x_ceil,y_ceil,z_ceil,x_floor,y_floor,z_floor

      ALPHA = 1.0D3	!scaling factor of initial value
      EPS   = 1.0D-7
      
      VMIN  = 0.0D0	!fixed Dirichlet BC
      PHI   = VMIN
      VMAX  = ALPHA	!fixed edge values

      DO IX=1,NX
        XI(IX)=XL+(IX-1.0D0)*DX
      END DO
      DO IY=1,NY
        YI(IY)=YL+(IY-1.0D0)*DX
      END DO
      DO IZ=1,NZ
        ZI(IZ)=ZL+(IZ-1.0D0)*DX
      END DO

!      G=0.D0   ! DGT commented this line on April 20, 2012
      G = 1.D0  ! DGT added this line on April 20,2012

!      TAU = TAUVAL                    ! rp 
!	
!	PHI = VMIN	! DGT entered on Nov 18, 2011
!
!      DO IA=1,N_ATOM
!         R = ATOM_R(IA) +  TAU !ATOM_R(IA)*BETA
!         R2=R**2
!	z_ceil = (ATOM_Z(IA)-ZL-R)/DX+1.0
!	z_floor = (ATOM_Z(IA)-ZL+R)/DX + 1.0

!!         DO IZ=CEILING((ATOM_Z(IA)-ZL-R)/DX+1),FLOOR((ATOM_Z(IA)-ZL+R)/DX+1)
!	  DO IZ = CEILING(z_ceil), FLOOR(z_floor)
!            DIST2=(ZL+(IZ-1)*DX-ATOM_Z(IA))**2
!            RXY2=ABS(R2-DIST2)
!            RXY=SQRT(RXY2)

!	y_ceil = (ATOM_Y(IA)-YL-RXY)/DX+1.0
!	y_floor = (ATOM_Y(IA)-YL+RXY)/DX+1.0

!!            DO IY=CEILING((ATOM_Y(IA)-YL-RXY)/DX+1),FLOOR((ATOM_Y(IA)-YL+RXY)/DX+1)
!	DO IY=CEILING(y_ceil),FLOOR(y_floor)
!               DIST2=(YL+(IY-1)*DX-ATOM_Y(IA))**2
!               RX=SQRT(ABS(RXY2-DIST2))
!        x_ceil = (ATOM_X(IA)-XL-RX)/DX+1.0
!	x_floor = (ATOM_X(IA)-XL+RX)/DX+1.0
    
!!	   DO IX=CEILING((ATOM_X(IA)-XL-RX)/DX+1),FLOOR((ATOM_X(IA)-XL+RX)/DX+1)
! 	   DO IX = CEILING(x_ceil), FLOOR(x_floor)
!	                 PHI(IX,IY,IZ) = VMAX
!                  G  (IX,IY,IZ) = 1.D0
!               END DO
!            END DO
!         END DO
!      END DO 


      DO IA=1,N_ATOM
         R=ATOM_R(IA)
         R2=R**2

	z_ceil = (ATOM_Z(IA)-ZL-R)/DX+1.0
        z_floor = (ATOM_Z(IA)-ZL+R)/DX + 1.0
!         DO IZ=CEILING((ATOM_Z(IA)-ZL-R)/DX+1),FLOOR((ATOM_Z(IA)-ZL+R)/DX+1)
	DO IZ = CEILING(z_ceil), FLOOR(z_floor)
            DIST2=(ZL+(IZ-1)*DX-ATOM_Z(IA))**2
            RXY2=ABS(R2-DIST2)
            RXY=SQRT(RXY2)        
	  y_ceil = (ATOM_Y(IA)-YL-RXY)/DX+1.0
          y_floor = (ATOM_Y(IA)-YL+RXY)/DX+1.0
 
!            DO IY=CEILING((ATOM_Y(IA)-YL-RXY)/DX+1),FLOOR((ATOM_Y(IA)-YL+RXY)/DX+1)

	DO IY=CEILING(y_ceil),FLOOR(y_floor)	
               DIST2=(YL+(IY-1)*DX-ATOM_Y(IA))**2
               RX=SQRT(ABS(RXY2-DIST2))
	
		x_ceil = (ATOM_X(IA)-XL-RX)/DX+1.0
        	x_floor = (ATOM_X(IA)-XL+RX)/DX+1.0

!               DO IX=CEILING((ATOM_X(IA)-XL-RX)/DX+1),FLOOR((ATOM_X(IA)-XL+RX)/DX+1)
		 DO IX = CEILING(x_ceil), FLOOR(x_floor)
                  G(IX,IY,IZ)=0.D0
                  PHI(IX,IY,IZ) = VMAX   ! DGT added this line on April 20, 2012 
               END DO
            END DO
         END DO
       END DO       

!	--- DGT commented out the following block of lines ---
!       PHI = VMIN

!       DO IA=1,N_ATOM
!         R = ATOM_R(IA) + tau !*BETA
!         R2=R**2
!	    z_ceil = (ATOM_Z(IA)-ZL-R)/DX+1.0
!        z_floor = (ATOM_Z(IA)-ZL+R)/DX + 1.0
!
!!         DO IZ=CEILING((ATOM_Z(IA)-ZL-R)/DX+1),FLOOR((ATOM_Z(IA)-ZL+R)/DX+1)
!		 DO IZ = CEILING(z_ceil), FLOOR(z_floor)
!            DIST2=(ZL+(IZ-1)*DX-ATOM_Z(IA))**2
!            RXY2=ABS(R2-DIST2)
!            RXY=SQRT(RXY2)
!          y_ceil = (ATOM_Y(IA)-YL-RXY)/DX+1.0
!          y_floor = (ATOM_Y(IA)-YL+RXY)/DX+1.0
!
!!            DO IY=CEILING((ATOM_Y(IA)-YL-RXY)/DX+1),FLOOR((ATOM_Y(IA)-YL+RXY)/DX+1)
!	  DO IY=CEILING(y_ceil), FLOOR(y_floor)	
!               DIST2=(YL+(IY-1)*DX-ATOM_Y(IA))**2
!               RX=SQRT(ABS(RXY2-DIST2))
!		 x_ceil = (ATOM_X(IA)-XL-RX)/DX+1.0
!                x_floor = (ATOM_X(IA)-XL+RX)/DX+1.0
!
!               DO IX=CEILING((ATOM_X(IA)-XL-RX)/DX+1),FLOOR((ATOM_X(IA)-XL+RX)/DX+1)
!		DO IX = CEILING(x_ceil), FLOOR(x_floor)
!                  PHI(IX,IY,IZ) = VMAX
!               END DO
!            END DO
!         END DO
!      END DO

!	--- end of commented block of lines---

!	write(*,*) "NX = ", NX
!	write(*,*) "NY = ", NY
!	write(*,*) "NZ = ", NZ		
!		DO IX=1,NX
!			DO IY=1,NY
!				DO IZ=1,NZ
!		write(50,'(5F12.4)') XL+(IX-1)*DX, YL+(IY-1)*DX, ZL+(IZ-1)*DX, G(IX,IY,IZ),PHI(IX,IY,IZ)
!				ENDDO
!			ENDDO
!		ENDDO 
!
!	STOP

      RETURN
      END 

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!	upwinding scheme, Neumann BC for PHI
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       SUBROUTINE UPWINDING(NX,NY,NZ,DX,DT,NT,G,PHI,PHITOTX)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Input: phi: the isosurface function 
!        potr: repulsive lj potential 
!        pota: attractive lj potential 
!        conms: pressure term 
! Output: renewed phi
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       IMPLICIT REAL*8 (A-H,O-Z)
       INTEGER :: NX,NY,NZ,NT
       DIMENSION G(NX,NY,NZ)
       DIMENSION PHI(NX,NY,NZ),PHINEW(NX,NY,NZ) 
       DIMENSION W1(-1:1),W2(-1:1),WXY(-1:1,-1:1),UP(0:1),DOWN(-1:0)
       REAL*8 ::DX,DT
       REAL*8 :: POTR(NX,NY,NZ),POTA(NX,NY,NZ)
	   REAL*8 :: PHIX,PHIY,PHIZ,DPHI,PHIXX,PHIYY,PHIZZ,PHIXY,PHIXZ,PHIYZ
	   REAL*8 :: PHITOTX(NX,NY,NZ)
	   REAL*8 :: PROT
 
       W1(-1)=-5.D-1/DX
       W1(0) =0.0D0
       W1(1) =5.D-1/DX

       W2(-1)=1.0D0/DX/DX
       W2(0) =-2.0D0/DX/DX
       W2(1) =1.0D0/DX/DX

       WXY(-1,-1)=2.5D-1/DX/DX
       WXY(-1,0)=0.0D0
       WXY(-1,1)=-2.5D-1/DX/DX
       WXY(0,-1)=0.0D0
       WXY(0,0)=0.0D0
       WXY(0,1)=0.0D0
       WXY(1,-1)=-2.5D-1/DX/DX
       WXY(1,0)=0.0D0
       WXY(1,1)=2.5D-1/DX/DX

       UP(0)=-1.0D0/DX
       UP(1)=1.0D0/DX

       DOWN(0)=1.0D0/DX
       DOWN(-1)=-1.0D0/DX

       PHINEW=PHI
DO IT=1,NT
	   DO IX=2,NX-1
	      DO IY=2,NY-1
		     DO IZ=2,NZ-1
			    IJK=(IX-1)*NZ*NY+(IY-1)*NZ+IZ
                IF(G(IX,IY,IZ).GT.2.D-2)THEN

				   PHIX=W1(-1)*PHI(IX-1,IY,IZ)+W1(1)*PHI(IX+1,IY,IZ)
                   PHIY=W1(-1)*PHI(IX,IY-1,IZ)+W1(1)*PHI(IX,IY+1,IZ)
				   PHIZ=W1(-1)*PHI(IX,IY,IZ-1)+W1(1)*PHI(IX,IY,IZ+1)
                   
				   PHIZZ=W2(-1)*PHI(IX,IY,IZ-1)+W2(0)*PHI(IX,IY,IZ)    &
                           +W2(1) *PHI(IX,IY,IZ+1)
				   PHIYY=W2(-1)*PHI(IX,IY-1,IZ)+W2(0)*PHI(IX,IY,IZ)    &
                           +W2(1) *PHI(IX,IY+1,IZ)
				   PHIXX=W2(-1)*PHI(IX-1,IY,IZ)+W2(0)*PHI(IX,IY,IZ)    &
                           +W2(1) *PHI(IX+1,IY,IZ)
						   
				   DPHI=(1.0D0+PHIX**2+PHIY**2)*PHIZZ+(1.0D0+PHIX**2+PHIZ**2)*PHIYY+  &
				        (1.0D0+PHIY**2+PHIZ**2)*PHIXX
				   
				   PHIXZ=WXY(-1,-1)*PHI(IX-1,IY,IZ-1)    &
                           +WXY(-1,1) *PHI(IX-1,IY,IZ+1)    &
                           +WXY(1,-1) *PHI(IX+1,IY,IZ-1)    &
                           +WXY(1,1)  *PHI(IX+1,IY,IZ+1)
				   PHIXY=WXY(-1,-1)*PHI(IX-1,IY-1,IZ)    &
                           +WXY(-1,1) *PHI(IX-1,IY+1,IZ)    &
                           +WXY(1,-1) *PHI(IX+1,IY-1,IZ)    &
                           +WXY(1,1)  *PHI(IX+1,IY+1,IZ)
				   PHIYZ=WXY(-1,-1)*PHI(IX,IY-1,IZ-1)    &
                           +WXY(-1,1) *PHI(IX,IY-1,IZ+1)    &
                           +WXY(1,-1) *PHI(IX,IY+1,IZ-1)    &
                           +WXY(1,1)  *PHI(IX,IY+1,IZ+1)
						   
				   DPHI=DPHI-2.D0*(PHIXY*PHIX*PHIY+PHIXZ*PHIX*PHIZ+PHIYZ*PHIY*PHIZ)		
				   GRAM=(1.D0+PHIX**2+PHIY**2+PHIZ**2)
                   DPHI=DPHI/GRAM+SQRT(GRAM)*PHITOTX(IX,IY,IZ)  ! add force field to the GF
				  
				                       !!!***			   
				   PHINEW(IX,IY,IZ)	=PHI(IX,IY,IZ)+DT*DPHI
				   if(PHINEW(IX,IY,IZ)>1000) PHINEW(IX,IY,IZ)=1000  ! modified by JC 02/18/2011
				   if(PHINEW(IX,IY,IZ)<0) PHINEW(IX,IY,IZ)=0	 ! modified by JC 02/18/2011
	   
				ENDIF 
			     
			 ENDDO 
		  ENDDO 
	   ENDDO   
	  
       DO IZ=2,NZ-1
         DO IY=2,NY-1
           DO IX=2,NX-1
             PHI(IX,IY,IZ)=PHINEW(IX,IY,IZ)
           END DO
         END DO
       END DO

  ENDDO 

  RETURN
  END 


