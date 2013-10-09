     
SUBROUTINE ADICOR(U,G,DELTAH,DELTAT,NX,NY,NZ,UX)
   !IMPLICIT NONE 
   IMPLICIT REAL*8(A-H,O-Z)
   
   REAL*8,INTENT(IN)::G(NZ,NY,NX)
   REAL*8::U(NZ,NY,NX),DELTAH,DELTAT,V(NX*NY*NZ),PHIX,PHIY,PHIZ,PHIXZ,PHIYZ,PHIXY
   REAL*8::A2(NZ,NY,NX),A3(NZ,NY,NX)
   REAL*8::A1XYZ(NX*NY*NZ),A2XYZ(NX*NY*NZ),CXYZ(NY*NX*NZ)
   REAL*8::UO(NZ,NY,NX),UOXYZ(NZ*NX*NY),A,B
   REAL*8::UX(NZ,NY,NX)

!	WRITE(*,*)"ENTER ADICOR"
	 NXYZ=NY*NX*NZ 
	 A=2.D0*DELTAH
	 B=(4.D0*DELTAH**2)
	 DO K=2,NZ-1   
       DO I=2,NY-1     !INITIALIZE THE COEFFICIENT OF A1,A2
	    DO J=2,NX-1
		   M=(K-1)*NX*NY+NX*(I-1)+J
		   KK=K+1
		   KK1=K-1
		   II=I+1
		   II1=I-1
		   JJ=J+1
		   JJ1=J-1
		   PHIX=(U(K,I,JJ)-U(K,I,JJ1))/A
           PHIY=(U(K,II,J)-U(K,II1,J))/A
		   PHIZ=(U(KK,I,J)-U(KK1,I,J))/A
		   PHIXY=(U(K,II,JJ)-U(K,II,JJ1)-U(K,II1,JJ)+U(K,II1,JJ1))/B
		   PHIXZ=(U(KK,I,JJ)-U(KK,I,JJ1)-U(KK1,I,JJ)+U(KK1,I,JJ1))/B
		   PHIYZ=(U(KK,II,J)-U(KK,II1,J)-U(KK1,II,J)+U(KK1,II1,J))/B
		   DE=(1+PHIX**2+PHIY**2+PHIZ**2)
		   A3(K,I,J)=G(K,I,J)*DELTAT*(1+PHIX**2+PHIY**2)/DE
		   A2(K,I,J)=G(K,I,J)*DELTAT*(1+PHIX**2+PHIZ**2)/(DE)
		   A1XYZ(M)=G(K,I,J)*DELTAT*(1+PHIY**2+PHIZ**2)/(DE)
		   CXYZ(M)=G(K,I,J)*(2.D0*(PHIX*PHIY*PHIXY+PHIX*PHIZ*PHIXZ+PHIY*PHIZ*PHIYZ)*DELTAT/DE&
                                 -DELTAT*UX(K,I,J)*SQRT(DE)) 
		END DO 
	   END DO
	  END DO  
	 
	 	  
	 DO K=1,NZ
	   DO I=1,NY       
	    DO J=1,NX
		  M=(K-1)*NX*NY+NX*(I-1)+J
		  UO(K,I,J)=U(K,I,J)
		  V(M)=U(K,I,J)
		  A2XYZ(M)=A2(K,I,J)
		  UOXYZ(M)=A3(K,I,J)
		  
		END DO 
	   END DO 
	 END DO 
	!WRITE(*,*)"ENTER ADI1"
       CALL ADI1(V,A1XYZ,A2XYZ,UOXYZ,CXYZ,NX,NY,NZ,NXYZ,DELTAH) !STEP 1 IN ADI

	 DO K=1,NZ
	   DO I=1,NY                   !CHANGE THE DIRECTION OF V TO Y
	     DO J=1,NX
		   U(K,I,J)=V((K-1)*NX*NY+(I-1)*NX+J)
		 END DO 
	    END DO 
	 END DO 
     
	 DO K=1,NZ                      
       DO J=1,NX
	      DO I=1,NY
		     M=(K-1)*NX*NY+(J-1)*NY+I
		     V(M)=U(K,I,J)
		     A2XYZ(M)=A2(K,I,J)
		     UOXYZ(M)=UO(K,I,J)
		  END DO 
	   END DO 
	  END DO 

	 CALL ADI2(V,A2XYZ,NY,NX,NZ,NXYZ,DELTAH,UOXYZ) !STEP 2 IN ADI-CREATE THE N^(k+1)

	 DO K=1,NZ
	    DO J=1,NX
	       DO I=1,NY
		      U(K,I,J)=V((K-1)*NX*NY+(J-1)*NY+I)
		   END DO 
		END DO 
	 END DO 

	 DO I=1,NY
	    DO J=1,NX
		   DO K=1,NZ
		      M=(I-1)*NX*NZ+(J-1)*NZ+K
		      V(M)=U(K,I,J)
			  A1XYZ(M)=A3(K,I,J)
			  UOXYZ(M)=UO(K,I,J)
  		   END DO
		END DO 
	  END DO 
	  
	  CALL ADI2(V,A1XYZ,NZ,NX,NY,NXYZ,DELTAH,UOXYZ)
	  
	  DO I=1,NY
	     DO J=1,NX
		    DO K=1,NZ
			   M=(I-1)*NZ*NX+(J-1)*NZ+K
			   U(K,I,J)=V(M)
			END DO 
		 END DO 
	  END DO 
	   
   RETURN 
  END SUBROUTINE ADICOR


  SUBROUTINE ADI1(V,A1XYZ,A2XYZ,A3XYZ,CXYZ,NX,NY,NZ,NXYZ,DELTAH)
  !IMPLICIT NONE 
  IMPLICIT REAL*8(A-H,O-Z)
    REAL*8,INTENT(IN)::A1XYZ(NXYZ),A2XYZ(NXYZ),A3XYZ(NXYZ),CXYZ(NXYZ)
  REAL*8:: R(NXYZ),ELEMENT,A(NXYZ),B(NXYZ),V(NXYZ),DELTAH
  
  ELEMENT=1/DELTAH**2
  !WRITE(*,*)"ENTER ADI1"  
                   !compute a,b,c
    B=1.D0
    A=0.D0
 

  DO K=2,NZ-1
    DO I=2,NY-1
       DO J=2,NX-1
	    M=(K-1)*NX*NY+(I-1)*NX+J
	    B(M)=1.0D0+ELEMENT*A1XYZ(M)
		A(M)=-ELEMENT*A1XYZ(M)/2.D0
	   END DO 
	 END DO 
  END DO 
  
  K=(NZ-1)*NX*NY
  DO I=1,NY
     DO J=1,NX
	    M=(I-1)*NX+J            !FOR THE TOP AND BOTTOM SURFACE POINTS
   	    R(M)=V(M)
	    M=M+K
	    R(M)=V(M)
      END DO 
  END DO 
  
  I=(NY-1)*NX
  DO K=1,NZ                      !FOR THE FRONT AND BEHIND SURFACE POINTS
     DO J=1,NX
	    M=(K-1)*NX*NY+J
		R(M)=V(M)
	    M=M+I
	    R(M)=V(M)
	 END DO 
  END DO 

  J=(NX-1)                       !FOR THE LEFT AND RIGHT SURFACE POINTS
  DO K=1,NZ
    DO  I=1,NY
	   M=(K-1)*NX*NY+(I-1)*NX+1
	   R(M)=V(M)
	   M=M+J
	   R(M)=V(M)
	 END DO 
   END DO 
  
  DO K=2,NZ-1 
    DO I=2,NY-1   
      DO J=2,NX-1
	    M=(K-1)*NX*NY+(I-1)*NX+J
	    R(M)=(1.D0-2.D0*ELEMENT*A2XYZ(M)-2.D0*ELEMENT*A3XYZ(M))*V(M)-ELEMENT*A1XYZ(M)*V(M)&
		      +ELEMENT*A2XYZ(M)*(V(M-NX)+V(M+NX))+ELEMENT*A3XYZ(M)*(V(M-NX*NY)+V(M+NX*NY))&
                    +0.5D0*ELEMENT*A1XYZ(M)*(V(M-1)+V(M+1))-CXYZ(M)
	   END DO 
	 END DO 
  END DO 

  CALL TRIDAG(A,B,A,R,V,NXYZ)  
  RETURN 
  END SUBROUTINE ADI1


  SUBROUTINE ADI2(V,AXYZ,NX,NY,NZ,NXYZ,DELTAH,UOXYZ)
  IMPLICIT REAL*8(A-H,O-Z)
  !IMPLICIT NONE 
  INTEGER::NX,NY,NZ,NXYZ,I,J,K,M
  REAL*8,INTENT(IN)::AXYZ(NXYZ),UOXYZ(NXYZ)
  REAL*8:: R(NXYZ),ELEMENT,A(NXYZ),B(NXYZ),V(NXYZ),DELTAH
  
  ELEMENT=1/DELTAH**2
  !WRITE(*,*)"ENTER ADI2"  
                    !compute a,b,c
       B=1.0D0
	   A=0.0D0
	   
    
    DO K=2,NZ-1
       DO I=2,NY-1
          DO J=2,NX-1
	         M=(K-1)*NX*NY+(I-1)*NX+J
	         B(M)=1.0D0+ELEMENT*AXYZ(M)
		     A(M)=-ELEMENT*AXYZ(M)/2.D0
		   END DO 
	 END DO 
  END DO 
  
  K=(NZ-1)*NX*NY
  DO I=1,NY
     DO J=1,NX
	    M=(I-1)*NX+J            !FOR THE TOP AND BOTTOM SURFACE POINTS
   	    R(M)=V(M)
	    M=M+K
	    R(M)=V(M)
      END DO 
  END DO 
  
  I=(NY-1)*NX
  DO K=1,NZ                      !FOR THE FRONT AND BEHIND SURFACE POINTS
     DO J=1,NX
	    M=(K-1)*NX*NY+J
		R(M)=V(M)
	    M=M+I
	    R(M)=V(M)
	 END DO 
  END DO 

  J=(NX-1)                       !FOR THE LEFT AND RIGHT SURFACE POINTS
  DO K=1,NZ
    DO  I=1,NY
	   M=(K-1)*NX*NY+(I-1)*NX+1
	   R(M)=V(M)
	   M=M+J
	   R(M)=V(M)
	 END DO 
   END DO 
  
  DO K=2,NZ-1 
    DO I=2,NY-1   
      DO J=2,NX-1
	    M=(K-1)*NX*NY+(I-1)*NX+J
	    R(M)=V(M)+ELEMENT*AXYZ(M)*UOXYZ(M)-ELEMENT*AXYZ(M)*(UOXYZ(M-1)+UOXYZ(M+1))/2.D0
	   END DO 
	 END DO 
  END DO 

  CALL TRIDAG(A,B,A,R,V,NXYZ)

  RETURN 
  END SUBROUTINE ADI2

   SUBROUTINE tridag(a,b,c,r,u,n)
  IMPLICIT REAL*8(A-H,O-Z)

  INTEGER:: n,NMAX
  REAL*8:: a(n),b(n),c(n),r(n),u(n)
  PARAMETER (NMAX=1000000)
!Solves f or a vector u(1:n) of length n the tridiagonal linear set given by equation (2.4.1).
!a(1:n), b(1:n), c(1:n), and r(1:n) are input vectors and are not modified.
!Parameter: NMAX is the maximum expected value of n.
  INTEGER::j
  REAL*8:: bet,gam(NMAX) !One vector of workspace, gam is needed.
  if(b(1).eq.0.)pause !'tridag: rewrite equations'
!If this happens then you should rewrite your equations as a set of order N ? 1, with u2
!trivially eliminated.
  bet=b(1)
  u(1)=r(1)/bet
  do j=2,n ! Decomposition and forward substitution.
     gam(j)=c(j-1)/bet
     bet=b(j)-a(j)*gam(j)
     if(bet.eq.0.)pause !'tridag failed' Algorithm fails; see below.
     u(j)=(r(j)-a(j)*u(j-1))/bet
  enddo 
  do  j=n-1,1,-1 !Backsubstitution.
        u(j)=u(j)-gam(j+1)*u(j+1)
	  
  enddo 
  return
END 


