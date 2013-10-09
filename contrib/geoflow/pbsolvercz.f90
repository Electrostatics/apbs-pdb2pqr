 
 MODULE MOD_3D
 
 REAL*8, DIMENSION(:),ALLOCATABLE::SA
 INTEGER:: NMAX
 INTEGER,ALLOCATABLE:: IJA(:)

 END MODULE 


 SUBROUTINE PBSOLVER(EPS,PHI,BGF,NX,NY,NZ,DCEL,IPREC,TOL,IPBIN,ITER,TITER)

 USE MOD_3D   ! for poisson solver 

 IMPLICIT REAL*8(A-H,O-Z)
     
	 INTEGER:: NX,NY,NZ,IPREC
     REAL*8 :: BGF(NX*NY*NZ),EPS(NX,NY,NZ),PHI(NX,NY,NZ)
     REAL*8 :: EPSHALF(NX,NY,NZ,3),U(NX*NY*NZ)
	 INTEGER:: ITOL,ITMAX,ITER,ICOUNT,IP(6),IERR,IPBIN
     REAL*8 :: TOL,ERR,DCEL,TITER
	 REAL*8 :: WEIT(6)
	 real*8 :: maxerr,idac  ! for test 
     
	 idac=0
	 NMAX=NX*NY*NZ*7
	 ALLOCATE(IJA(NMAX),SA(NMAX))

	 DO I=1,NX-1  ! epsilon value setting
	     DO J=1,NY-1
		     DO K=1,NZ-1
			    
                EPSHALF(I,J,K,1)=(EPS(I+1,J,K)+EPS(I,J,K))/2.D0
                EPSHALF(I,J,K,2)=(EPS(I,J+1,K)+EPS(I,J,K))/2.D0
				EPSHALF(I,J,K,3)=(EPS(I,J,K+1)+EPS(I,J,K))/2.D0
               
			 ENDDO 
		 ENDDO 
	 ENDDO 

	 ICOUNT=NX*NY*NZ+1                              ! write sa ija
     IJA(1)=ICOUNT+1
     DO I=1,NX
        DO J=1,NY
            DO K=1,NZ
	           IJK=(I-1)*NZ*NY+(J-1)*NZ+K 
               IF(I==1.OR.I==NX.OR.J==1.OR.J==NY.OR.K==1.OR.K==NZ)THEN   ! boundary 
		           SA(IJK)=1.D0 
		       ELSE   
			       SA(IJK)=-(EPSHALF(I,J,K,1)+EPSHALF(I-1,J,K,1)+&
				        EPSHALF(I,J,K,2)+EPSHALF(I,J-1,K,2)+&
						EPSHALF(I,J,K,3)+EPSHALF(I,J,K-1,3))/DCEL**2
				   WEIT(1)=EPSHALF(I-1,J,K,1)
				   WEIT(2)=EPSHALF(I,J-1,K,2)
				   WEIT(3)=EPSHALF(I,J,K-1,3)
				   WEIT(4)=EPSHALF(I,J,K,3)
				   WEIT(5)=EPSHALF(I,J,K,2)
				   WEIT(6)=EPSHALF(I,J,K,1)

	               IP(1)=(IJK-NZ*NY)
	               IP(2)=(IJK-NZ)
	               IP(3)=(IJK-1)
	               IP(4)=(IJK+1)
	               IP(5)=IJK+NZ
	               IP(6)=IJK+NZ*NY
	               DO II=1,6
	                    ICOUNT=ICOUNT+1
		                IJA(ICOUNT)=IP(II)
		                SA(ICOUNT)=WEIT(II)/DCEL**2     !!?
	               ENDDO 
                ENDIF
                IJA(IJK+1)=ICOUNT+1
		     ENDDO 
		  ENDDO 
	  ENDDO  
      
      ! check the accuracy before the solver!!!!!!!!!!!!!!!!!!!!
	  if(idac==1)then
	  maxerr=0.d0
	  sum=0.d0
      do i=1,nx
	      do j=1,ny
		      do k=1,nz
			      
			      x=xvalue(i)
				  y=yvalue(j)
				  z=zvalue(k)
			      ijk=(i-1)*ny*nz+(j-1)*nz+k
				  jj=ija(ijk+1)-ija(ijk)
				  if(jj>0)then
				      sum=sa(ijk)*cos(x)*cos(y)*cos(z)
                      sum=sum+sa(ija(ijk))*cos(x-dcel)*cos(y)*cos(z)
					  sum=sum+sa(ija(ijk)+1)*cos(x)*cos(y-dcel)*cos(z)
					  sum=sum+sa(ija(ijk)+2)*cos(x)*cos(y)*cos(z-dcel)
					  sum=sum+sa(ija(ijk)+3)*cos(x)*cos(y)*cos(z+dcel)
					  sum=sum+sa(ija(ijk)+4)*cos(x)*cos(y+dcel)*cos(z)
					  sum=sum+sa(ija(ijk)+5)*cos(x+dcel)*cos(y)*cos(z)
                      err=abs(bgf(ijk)-sum)
!					  if(err>5.d0)then
!					     print*,i,j,k
!					  endif 
				      maxerr=max(err,maxerr)
				  else
				      sum=cos(x)*cos(y)*cos(z)
				      err=abs(bgf(ijk)-sum)
				      maxerr=max(err,maxerr)
				  endif
			  enddo 
		  enddo 
	  enddo 

!	  print*,'the precheck of the matrix',maxerr
      endif 

	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO I=1,NX
	     DO J=1,NY
		    DO K=1,NZ
			   IJK=(I-1)*NY*NZ+(J-1)*NZ+K
			   U(IJK)=PHI(I,J,K)
			ENDDO 
		 ENDDO 
	  ENDDO 
      IF(IPBIN==0)U=0.D0                ! set all initials zero
      
      CALL CPU_TIME(TIME1)
      ITOL=1
      ITMAX=1000
	  IF(IPREC==0)THEN
           CALL LINBCG(NX*NY*NZ,BGF,U,ITOL,TOL,ITMAX,ITER,ERR)
	  ELSEIF(IPREC==1)THEN 
	       CALL PRECBCG(NX*NY*NZ,ICOUNT,SA,IJA,BGF,U,ITOL,TOL,ITMAX,ITER, ERR, IERR)
	  ELSE 
	       PRINT*,"please input the option of linear system"
	  ENDIF 
      CALL CPU_TIME(TIME2) 
      TITER=TIME2-TIME1
!      PRINT*,ITER,ERR


      DO I=1,NX
	     DO J=1,NY
		    DO K=1,NZ
			   IJK=(I-1)*NY*NZ+(J-1)*NZ+K
			   PHI(I,J,K)=U(IJK)
			ENDDO 
		 ENDDO 
	  ENDDO 

	  DEALLOCATE(IJA,SA)
    
  RETURN 
  END SUBROUTINE 

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE linbcg(n,b,x,itol,tol,itmax,iter,err)
use mod_3d
Implicit Double Precision(A-H, O-Z)
PARAMETER (EPS=1.d-20)
!----------------------------------------------------------------------------------------
!USES atimes,asolve,snrm
!Solves A ¡¤ x = b for x(1:n), given b(1:n), by the iterative biconjugate gradient method.
!On input x(1:n) should be set to an initial guess of the solution (or all zeros); itol is
!1,2,3, or 4, specifying which convergence test is applied (see text); itmax is the maximum
!number of allowed iterations; and tol is the desired convergence tolerance. On output,
!x(1:n) is reset to the improved solution, iter is the number of iterations actually taken,
!and err is the estimated error. The matrix A is referenced only through the user-supplied
!routines atimes, which computes the product of either A or its transpose on a vector; and
!asolve, which solves .A ¡¤ x = b or .AT
!¡¤ x = b for some preconditioner matrix .A (possibly
!the trivial diagonal part of A).
!--------------------------------------------------------------------------------------------
real*8,allocatable :: p(:),pp(:),r(:),rr(:),z(:),zz(:)
dimension b(N),x(N)

allocate(p(NMAX),pp(NMAX),r(NMAX),rr(NMAX),z(NMAX),zz(NMAX))

iter=0								!Calculate initial residual.
call atimes(n,x,r,0)				!Input to atimes is x(1:n), output is r(1:n);
            						!the final 0 indicates that the matrix (not
									!its transpose) is to be used.
do j=1,n
	r(j)=b(j)-r(j)
	rr(j)=r(j)
enddo
!call atimes(n,r,rr,0)				!Uncomment this line to get the ¡°minimum
									!residual¡± variant of the algorithm. 
if (itol.eq.1) then
	bnrm=snrm(n,b,itol)
	call asolve(n,r,z,0)			!Input to asolve is r(1:n), output is z(1:n);
									!the final 0 indicates that the matrix .A
									!(not its transpose) is to be used.
elseif (itol.eq.2) then
	call asolve(n,b,z,0)
	bnrm=snrm(n,z,itol)
	call asolve(n,r,z,0)
elseif (itol.eq.3.or.itol.eq.4) then
	call asolve(n,b,z,0)
	bnrm=snrm(n,z,itol)
	call asolve(n,r,z,0)
	znrm=snrm(n,z,itol)
else
		pause 'illegal itol in linbcg'
endif

100 if (iter.le.itmax) then			!Main loop.
		iter=iter+1
		call asolve(n,rr,zz,1)		!Final 1 indicates use of transpose matrix .AT.
		bknum=0.d0
		do j=1,n					!Calculate coefficient bk and direction vectors p and pp. 

			bknum=bknum+z(j)*rr(j)
		enddo
		if (iter.eq.1) then
		do j=1,n
			p(j)=z(j)
			pp(j)=zz(j)
		enddo
	else
		bk=bknum/bkden
		do j=1,n
			p(j)=bk*p(j)+z(j)
			pp(j)=bk*pp(j)+zz(j)
		enddo
	endif
bkden=bknum							!Calculate coefficient ak, new iterate x, and
									!new residuals r and rr. 


call atimes(n,p,z,0)
akden=0.d0
do j=1,n
	akden=akden+z(j)*pp(j)
enddo
ak=bknum/akden

call atimes(n,pp,zz,1)
do j=1,n
	x(j)=x(j)+ak*p(j)
	r(j)=r(j)-ak*z(j)
	rr(j)=rr(j)-ak*zz(j)
enddo
call asolve(n,r,z,0)				!Solve .A¡¤z = r and check stopping criterion.

if (itol.eq.1) then
	err=snrm(n,r,itol)/bnrm
elseif(itol.eq.2)then
	err=snrm(n,z,itol)/bnrm
elseif(itol.eq.3.or.itol.eq.4)then
	zm1nrm=znrm
	znrm=snrm(n,z,itol)
	if(abs(zm1nrm-znrm).gt.EPS*znrm) then
		dxnrm=abs(ak)*snrm(n,p,itol)
		err=znrm/abs(zm1nrm-znrm)*dxnrm
	else
		err=znrm/bnrm				!Error may not be accurate, so loop again.
		goto 100
	endif
	xnrm=snrm(n,x,itol)
	if(err.le.0.5d0*xnrm) then
		err=err/xnrm
	else
		err=znrm/bnrm				!Error may not be accurate, so loop again.
	goto 100
	endif
endif

!write (*,*) 'iter=', iter, 'err=',err
if(err.gt.tol) goto 100
endif

deallocate(p,pp,r,rr,z,zz)
return
END
!
!--------------------------------------------------------------------------
!
!The routine linbcg uses this short utility for computing vector norms:
FUNCTION snrm(n,sx,itol)
use mod_3d
Implicit Double Precision(A-H, O-Z)
INTEGER n,itol,i,isamax
DOUBLE PRECISION sx(n),snrm
!Compute one of two norms for a vector sx(1:n), as signaled by itol. Used by linbcg.
if (itol.le.3)then
	snrm=0.d0
	do i=1,n						!Vector magnitude norm.
		snrm=snrm+sx(i)**2
	enddo
	!print*,snrm
	snrm=dsqrt(snrm)
else
	isamax=1
	do i=1,n						!Largest component norm.
		if(dabs(sx(i)).gt.dabs(sx(isamax))) isamax=i
	enddo
	snrm=dabs(sx(isamax))
endif
return
END
!
!-------------------------------------------------------------------------------------------------------------------
!
SUBROUTINE atimes(n,x,r,itrnsp)
use mod_3d
Implicit Double Precision(A-H, O-Z)
dimension x(n),r(n)
!---------------------------------------------------------------------------
! USES dsprsax,dsprstx DOUBLE PRECISION versions of sprsax and sprstx.
!--------------------------------------------------------------------------
if (itrnsp.eq.0) then
	call dsprsax(x,r,n)
else
	call dsprstx(x,r,n)
endif
return
END
!
!------------------------------------------------------------------------------------------------------------
!
SUBROUTINE asolve(n,b,x,itrnsp)
use mod_3d
Implicit Double Precision(A-H, O-Z)
dimension x(n),b(n)
do i=1,n
	x(i)=b(i)/sa(i)					!The matrix .A is the diagonal part of A, stored in
									!the first n elements of sa. Since the transpose
									!matrix has the same diagonal, the flag itrnsp is not used.
enddo
return
END
!
!------------------------------------------------------------------------------------
!
SUBROUTINE dsprsax(x,b,n)
use mod_3d
Implicit Double Precision(A-H, O-Z)
dimension b(n),x(n)
!------------------------------------------------------------------------------------------------
!Multiply a matrix in row-index sparse storage arrays sa and ija by a vector x(1:n), giving
!a vector b(1:n).
!-----------------------------------------------------------------------------------------------
! if (ija(1).ne.n+2)then
 !    print *,ija(1),n+2
!	 pause 'mismatched vector and matrix in sprsax'
! endif
do i=1,n
	b(i)=sa(i)*x(i)					!Start with diagonal term.
	do k=ija(i),ija(i+1)-1			!Loop over off-diagonal terms.
		b(i)=b(i)+sa(k)*x(ija(k))
	enddo
enddo
return
END
!
!---------------------------------------------------------------------------------------------------
!
SUBROUTINE dsprstx(x,b,n)
use mod_3d
Implicit Double Precision(A-H, O-Z)
dimension b(n),x(n)
!------------------------------------------------------------------------------------------------
!Multiply the transpose of a matrix in row-index sparse storage arrays sa and ija by a
!vector x(1:n), giving a vector b(1:n).
!------------------------------------------------------------------------------------------------

if (ija(1).ne.n+2) pause 'mismatched vector and matrix in sprstx'
	do i=1,n						!Start with diagonal terms.
		b(i)=sa(i)*x(i)
	enddo
	do i=1,n						!Loop over off-diagonal terms.
		do k=ija(i),ija(i+1)-1
			j=ija(k)
			b(j)=b(j)+sa(k)*x(i)
		enddo
	enddo
return
END
!
!--------------------------------------------------------------------------
!
SUBROUTINE dsprsin(a,n,np,thresh)
use mod_3d
Implicit Double Precision(A-H, O-Z)
dimension a(np,np)
!----------------------------------------------------------------------------------------
!Converts a square matrix a(1:n,1:n) with physical dimension np into row-indexed sparse
!storage mode. Only elements of a with magnitude ¡Ýthresh are retained. Output is in
!two linear arrays with physical dimension nmax (an input parameter): sa(1:) contains
!array values, indexed by ija(1:). The logical sizes of sa and ija on output are both
!ija(ija(1)-1)-1 (see text).
!-------------------------------------------------------------------------------------
do j=1,n								!Store diagonal elements.
	sa(j)=a(j,j)
enddo 
ija(1)=n+2								!Index to 1st row off-diagonal element, if any.
k=n+1
do i=1,n								!Loop over rows.
	do j=1,n							!Loop over columns.
		if(dabs(a(i,j)).ge.thresh)then
			if(i.ne.j)then				!Store off-diagonal elements and their columns.
				k=k+1
				if(k.gt.nmax) pause 'nmax too small in sprsin'
				sa(k)=a(i,j)
				ija(k)=j
			endif
		endif
	enddo
	ija(i+1)=k+1						!As each row is completed, store index to next.
enddo
return
END

SUBROUTINE PRECBCG(NDIM,NMAX,SA,IJA,BR,U,ITOL,TOL,ITMAX,ITER, ERR, IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  A routine to solver the linear system with preconditioner from BCG matrix information
!   Input:
!         Ndim : in    dimension of matrix=nx*ny*nz
!         Nmax : in    dimension of sa ija chart 
!         Sa, IJa: in   store the information of non-zero element 
!         Br     : in    vector of right hand side of AX=B 
!         U      : inout  in vector of initial guess
!                :        out  solution of linear system 
!         Itol   : in   how to determine the stop creteria 
!         Tol    : in   stop creteria 
!         Itmax  : in   maximum steps of iteration 
!         Iter   : out  number of iteration 
!         err    : out  err when iteration stops 
!         ierr   : out  =0 iteration goes well 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPLICIT NONE 
    INTEGER:: NDIM, NMAX,IJA(NMAX),ITOL,ITMAX, ITER,IERR  
    REAL*8 :: SA(NMAX), U(NDIM),BR(NDIM),TOL,ERR
    INTEGER:: NELT,ISYM,NSAVE,IUNIT,LENW,LENIW
    REAL*8 ,DIMENSION(:),ALLOCATABLE:: A
    REAL*8,ALLOCATABLE::RWORK(:)
    INTEGER, ALLOCATABLE:: IWORK(:),IA(:),JA(:)
    INTEGER:: ICOUNT, II, JJ, I
    
       NELT=IJA(NDIM+1)-2
	   ALLOCATE(A(NELT))
	   ALLOCATE(IA(NELT),JA(NELT))
       ISYM=0
       NSAVE=1
       IUNIT=0     ! no writing

       LENW=18*NDIM
       LENIW=10*NDIM
       ALLOCATE(RWORK(LENW))
       ALLOCATE(IWORK(LENIW))
       
	   ! transfer the ija sa matrix information to A IA(row number) JA(column)
       ICOUNT=1              
	   DO I=1,NDIM

	      A(ICOUNT)=SA(I)     ! input the diognal element
		  IA(ICOUNT)=I
		  JA(ICOUNT)=I
		  ICOUNT=ICOUNT+1

		  JJ=IJA(I+1)-IJA(I)  ! input non-diognal element
		  IF(JJ.GT.0)THEN
              II=1
			  DO WHILE(II.LT.JJ+1)
			       A(ICOUNT)=SA(IJA(I)+II-1)
				   IA(ICOUNT)=I
				   JA(ICOUNT)=IJA(IJA(I)+II-1)
				   ICOUNT=ICOUNT+1
				   II=II+1
			  ENDDO 
		   ENDIF
	    ENDDO
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		! linear algebra solver
		! detail of this linear solver can be seen in website 
        !  http://sdphca.ucsd.edu/slatec_top/source/dsluom.f
		CALL DSLUOM (NDIM, BR, U, NELT, IA, JA, A, ISYM, NSAVE, ITOL, &
         TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW)

        
        DEALLOCATE(A,IA,JA,RWORK,IWORK)
        RETURN 

       
END SUBROUTINE PRECBCG
