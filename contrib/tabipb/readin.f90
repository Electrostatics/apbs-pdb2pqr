!Read protein structure from msms
subroutine readin(idens,ichrpos)
!use IFPORT 
!use mesh_procedures
use molecule
use comdata
implicit double precision(a-h,o-z)
real*8 pos(3),vector(3)
integer nind(5)
CHARACTER(100) :: FHEAD


!Obtain path
	pathname='test_proteins/'
        lenpath = len(pathname)
	do while (pathname(lenpath:lenpath) .eq. ' ')
		lenpath = lenpath - 1
	enddo  


!Obtain protein name
	!write(*,*) 'Please input the protein name:'
	!read(*,*) fname
	fname='1a63'
	!write(*,*) 'Please input MSMS triangulation density in # per astrong^2:'
	!read(*,*) den
	den='10'	
	lenfname = len(fname)
	do while (fname(lenfname:lenfname) .eq. ' ')
		lenfname = lenfname - 1
	enddo 


!Read atom coordinates and partial charges
	OPEN(1,FILE=pathname(1:lenpath)//FNAME(1:lenfname)//".xyzr")            
	NATM = 0
		DO 
			READ(1,*,IOSTAT = MEOF) xxx, yyy, zzz, rrr
			IF(MEOF .LT. 0) EXIT
			NATM = NATM + 1
		END DO         
	CLOSE(1)
	
	OPEN(1,FILE=pathname(1:lenpath)//FNAME(1:lenfname)//".pqr")            
	nchr = 0
		DO 
			READ(1,*,IOSTAT = MEOF) xxx, yyy, zzz, rrr
			IF(MEOF .LT. 0) EXIT
			nchr = nchr + 1
		END DO         
	CLOSE(1)
	
	if (nchr .ne. natm) print *,'atoms completely inside atoms are deleted'
	
    !nchr=natm
	allocate(atmpos(3,natm),atmrad(natm),atmchr(nchr),chrpos(3,nchr),STAT=ierr)
	IF (ierr .NE. 0) THEN
    		WRITE(6,*) 'Error allocating atmpos, atmrad, atmchr, chrpos!'
    		STOP
	END IF


	open(1,FILE=pathname(1:lenpath)//FNAME(1:lenfname)//".xyzr")
	do i=1,natm
		read(1,*) atmpos(:,i), atmrad(i)
	enddo
	close(1)
    
    ! For sphere only
    rds=atmrad(1)
	open(1,FILE=pathname(1:lenpath)//FNAME(1:lenfname)//".pqr")
	do i=1,nchr
		read(1,*) chrpos(:,i),atmchr(i)
	enddo
	close(1)
	
	!chrpos(1,1)=(ichrpos-1)*0.1d0                               ! For multiple runs
	

    rslt=system('msms -if '//pathname(1:lenpath)//fname(1:lenfname)//".xyzr"//' -prob 1.4 -de ' &
    //den(1:5)//' -of '//pathname(1:lenpath)//fname(1:lenfname))    
      ! read the surface points
      OPEN(2,FILE=pathname(1:lenpath)//FNAME(1:lenfname)//".vert")

      READ(2,*) FHEAD
      READ(2,*) FHEAD
      READ(2,*) NSPT, ppp, qqq, rrr

      ALLOCATE(SPTPOS(3,NSPT), SPTNRM(3,NSPT), NATMAFF(NSPT), NSFTYPE(NSPT), STAT= ierr)
      IF (ierr .NE. 0) THEN
          WRITE(6,*) 'Error allocating SPTPOS, SPTNRM, NATMAFF, NSFTYPE!'
      STOP
      END IF

	SPTPOS=0.D0; SPTNRM=0.D0; NATMAFF=0; NSFTYPE=0;	
      
      DO I=1,NSPT
         READ(2,*) POS(1:3), VECTOR(1:3), KK, NAFF, NAFFT 
		 
         
	 SPTPOS(:,I) = POS;   SPTNRM(:,I) = VECTOR
         NATMAFF(I)  = NAFF;  NSFTYPE(I)  = NAFFT
      END DO
	 
      CLOSE(2)

! read the surface triangulization

      OPEN(3,FILE=pathname(1:lenpath)//FNAME(1:lenfname)//".face")

      READ(3,*) FHEAD
      READ(3,*) FHEAD
      READ(3,*) NFACE, PPP, QQQ, RRR

      ALLOCATE(NVERT(3,NFACE), MFACE(NFACE), STAT=ierr)
      IF (ierr .NE. 0) THEN
          WRITE(6,*) 'Error allocating NVERT, MFACE'
      STOP
      END IF

      NVERT=0; MFACE=0
      
      DO I=1,NFACE 
         READ(3,*) NIND(1:5) 
         NVERT(1:3,I) = NIND(1:3);  MFACE(I) = NIND(4)
      END DO
      CLOSE(3)
      call surface_area(s_area) ! the post-MSMS code
      print *,'surface area=', real(s_area)


End
!------------------------------------------------------------------------
subroutine surface_area(s_area)
use molecule
implicit double precision(a-h,o-z)
integer iface(3),jface(3),nfacenew,ialert
real*8 face(3,3),s_area,face_old(3,3),xx(3),yy(3),cpu1,cpu2
real*8,allocatable:: nvert_copy(:,:)

        print *,'# of surfaces=',nface,' # of surface points=',nspt
        call cpu_time(cpu1)
        s_area=0.d0
        nfacenew=nface
        allocate(nvert_copy(3,nface))
        nvert_copy=nvert
        do i=1,nface
                iface=nvert(:,i)
                xx=0.d0
                ialert=0;
                do j=1,3
                    face(:,j)=sptpos(:,iface(j))
                    xx=xx+1/3.d0*(face(:,j))
                enddo
                aa=sqrt(dot_product(face(:,1)-face(:,2),face(:,1)-face(:,2)))
                bb=sqrt(dot_product(face(:,1)-face(:,3),face(:,1)-face(:,3)))
                cc=sqrt(dot_product(face(:,2)-face(:,3),face(:,2)-face(:,3)))
                area_local=triangle_area(aa,bb,cc)
                
                do ii=max(1,i-10),i-1
                    jface=nvert(:,ii)
                    yy=0.d0
                    do j=1,3
                        face_old(:,j)=sptpos(:,jface(j))
                        yy=yy+1/3.d0*(face_old(:,j))
                    enddo
                    dist_local=dot_product(xx-yy,xx-yy)
                    if (dist_local<1.d-5) then
                       ialert=1
                       print *,i,ii,'particles are too close',dist_local
                    endif
                enddo
                if (area_local < 1.d-5 .or. ialert==1) then
                    print *,i,j,'small area=', area_local
                    ichanged=nface-nfacenew
                    nvert_copy(:,(i-ichanged):(nface-1))=nvert_copy(:,(i-ichanged+1):nface)
                    nfacenew=nfacenew-1
                endif
                s_area=s_area+area_local
        enddo
        print *,nface-nfacenew,' ugly faces are deleted'
        
        nface=nfacenew
        deallocate(nvert,STAT=ierr)
        IF (ierr .NE. 0) THEN
            WRITE(6,*) 'Error deallocating nvert'
            STOP
        EndIF
       
        allocate(nvert(3,nface),STAT=ierr)
        IF (ierr .NE. 0) THEN
            WRITE(6,*) 'Error allocating nvert'
            STOP
        EndIF
        nvert=nvert_copy(:,1:nface)
        deallocate(nvert_copy, STAT=ierr)
        print *,3
        IF (ierr .NE. 0) THEN
            WRITE(6,*) 'Error deallocating nvert_copy'
            STOP
        call cpu_time(cpu2)
        EndIF
        print *, 'total MSMS post-processing time =',cpu2-cpu1
end



!------------------------------------------------------------------------
function triangle_area(aa,bb,cc)
implicit double precision(a-h,o-z)
	s=0.5d0*(aa+bb+cc)
	triangle_area=sqrt(s*(s-aa)*(s-bb)*(s-cc))
end 

!---------------------------------------------------------------------------------
! For a new vertex, compared with all the stored vertex
! ifind .ne. 0:	if the vertex is already stored
! ifind=0:		if the vertex is completely new
Subroutine mesh_find_vertex(indx_sptpos,ifind,sptpos,nspt,sptpos_new)
implicit none
real*8 sptpos(3,nspt),sptpos_new(3),diff(3)
integer indx_sptpos,ifind,nspt,isptpos

ifind=0

do isptpos=1,indx_sptpos
	diff=sptpos_new-sptpos(:,isptpos)
	if (sqrt(dot_product(diff,diff))<1.d-10) then
		ifind = isptpos
		return
	endif
enddo

End

!#############################################################################################


!The face file contains three header lines followed by one triangle per line. 
!The first header line provides a comment and the file name of the sphere set. 
!The second header line holds comments about the content of the third line. 
!The third header line provides the number of triangles, the number of spheres in the set, 
!the triangulation density and the probe sphere radius. 

!The first three numbers are (1 based) vertex indices. 

!The next field can be: 
!1 for a triangle in a toric reen trant face, 
!2 for a triangle in a spheric reentrant face and 
!3 for a triangle in a contact face. 

!The last # on the line is the (1 based) face number in the analytical description of the solvent excluded surface. 
!These values are written in the following format ``%6d %6d %6d %2d %6d''.

!The vertex file contains three header lines (similar to the header in the .face file) 
!followed by one vertex per line and provides the coordinates (x,y,z) and the normals (nx,ny,nz) 
!followed by the number of the face (in the analytical description of the solvent excluded surface) 
!to which the vertex belongs. The vertices of the analytical surface have a value 0 in that field 
!and the vertices lying on edges of this surface have negative values. 
!The next field holds the (1 based) index of the closest sphere. 
!The next field is 
!1 for vertices which belong to toric reentrant faces (including ver tices of the analytical surface), 
!2 for vertices inside reentrant faces and 
!3 for vertices inside contact faces. 
!These values are written in the following format ``%9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %7d %7d %2d''.
