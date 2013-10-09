subroutine seteqb(bg,xyzr,pqr,natm,charget,corlocqt,epsilonsp,  &
            epsilon_p,epsilon_s,surfux,surfuy,surfuz,delphix_c, &
            delphiy_c,delphiz_c,iflag)
 use comdata
 implicit double precision(a-h,o-z)
 real*8,dimension(nx*ny*nz):: bg
 real*8:: xyzr(4,natm),pqr(natm),corlocqt(natm,8,3),charget(natm,8)
 real*8 :: vbdn,fp,epsilonsp,epsilon_p,epsilon_s
 real*8 :: surfux(nx,ny,nz),surfuy(nx,ny,nz),surfuz(nx,ny,nz)
 real*8 :: delphix_c(nx,ny,nz),delphiy_c(nx,ny,nz),delphiz_c(nx,ny,nz)
 integer :: iflag

    ! -----------
    !iflag = 0: sets the B matrix of the non-regularized generalized Poission equation
    !           that is solved for computing the electrostatic potential.
    !
    !iflag = 1: sets the B matrix for the regularized Poisson equation that is
    !           solved for computing the reaction field potential.
    !
    !-----------
 bg=0.d0

 do i=1,nx
    do j=1,ny
	    do k=1,nz
		    ijk=(i-1)*ny*nz+(j-1)*nz+k
		    x=xvalue(i)
	        y=yvalue(j)
			z=zvalue(k)
		    if(i<2.or.i>nx-1.or.j<2.or.j>ny-1.or.k<2.or.k>nz-1)then
			      x=xvalue(i)
                  y=yvalue(j)
                  z=zvalue(k)
                  vbdn=0.d0

                  do i1=1,natm
                      x_q=xyzr(1,i1)
	                  y_q=xyzr(2,i1)
	                  z_q=xyzr(3,i1)
	                  q_q=pqr(i1)
                      rr=sqrt((x-x_q)**2+(y-y_q)**2+(z-z_q)**2)

                      ! IF statement added by DGT on 4/26/2013
                      if(iflag.eq.0) then
	                  vbdn=vbdn+q_q/(epsilonsp*rr)
	                  else if(iflag.eq.1) then

                        ! added by DGT
	                  vbdn = vbdn + (q_q/rr)*(1.d0/epsilon_s - 1.d0/epsilon_p)
	                  endif

                  enddo
			      bg(ijk)=vbdn
			  else
			      fp=0.d0
			      if(iflag.eq.0) then
                  do jj=1,natm
                      do ii=1,8
                         if((x-corlocqt(jj,ii,1))**2+(y-corlocqt(jj,ii,2))**2+(z-corlocqt(jj,ii,3))**2.le.1.d-13)then
	                          fp=fp-4.d0*pi*charget(jj,ii)/deltax/deltay/deltaz
	                     endif
	                  enddo
                   enddo
                   else if(iflag.eq.1) then
                   fp = -(epsilon_p - epsilon_s)*(surfux(i,j,k)*delphix_c(i,j,k)    &
                   + surfuy(i,j,k)*delphiy_c(i,j,k) + surfuz(i,j,k)*delphiz_c(i,j,k))
                   endif
				   bg(ijk)=fp
				
				endif

            enddo
	     enddo
      enddo



     return
     end

    

     function xvalue(i)

     use comdata
     implicit double precision(a-h,o-z)

     xvalue=(i-1)*deltax+xleft 

     return
     end
!
!----------------------------------------------------------------------
!
     function yvalue(j)

     use comdata
     implicit double precision(a-h,o-z)

     yvalue=(j-1)*deltay+yleft 

     return
     end
!
!----------------------------------------------------------------------
!
     function zvalue(k)

     use comdata
     implicit double precision(a-h,o-z)

     zvalue=(k-1)*deltaz+zleft 
 
     return
     end
!
!-------------------------------------------

!
!------------------------------------------------------------------------
!

SUBROUTINE computeCoulombPotential(phi_c,epsilonsp,natm,xyzr,pqr)
    ! -------
    ! written by DGT
    ! last updated by DGT on 4/25/2013
    !-------
use comdata
implicit real*8(A-H,O-Z),integer(I-N)

    real*8::phi_c(nx,ny,nz)
    integer:: ix,iy,iz,iatm
    real*8:: xyzr(4,natm)
    real*8 :: pqr(natm)

    real*8 :: sum_term;
    real*8 :: xd,yd,zd,dist;


    do ix = 1,nx
        do iy = 1,ny
            do iz = 1,nz

        sum_term = 0.d0
                do iatm = 1,natm
                ! calculate the distance between a grid point and an atom location
               xd = xc(ix) - xyzr(1,iatm)
               yd = yc(iy) - xyzr(2,iatm)
               zd = zc(iz) - xyzr(3,iatm)

           dist = sqrt(xd*xd + yd*yd + zd*zd)

           ! calculate the coulombic potential at the grid point due to natm charges
           ! if the grid point is located at the atom charge, then do not consider
           ! that charge in the calculation.
           if(dist.ge.1.d-13) then
           sum_term = sum_term + pqr(iatm)/(epsilonsp*dist)

            endif
                enddo
                phi_c(ix,iy,iz) = sum_term
            enddo
         enddo
      enddo


return
END

