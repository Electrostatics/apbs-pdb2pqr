! charge distribution 

subroutine chargedist(atmpos,chratm,natm,charget,corlocqt,loc_qt,iatm)
use comdata

implicit double precision(a-h,o-z)
dimension iloc(3),charge(8),loc_q(8,3),corlocq(8,3)
dimension charget(natm,8),corlocqt(natm,8,3),loc_qt(natm,8,3)
dimension atmpos(4,natm),chratm(natm)

charge=0.d0
iloc=1

x_q=atmpos(1,iatm)
y_q=atmpos(2,iatm)
z_q=atmpos(3,iatm)
q_q=chratm(iatm)

i_q=inverx(x_q) ; j_q=invery(y_q) ; k_q=inverz(z_q)
x=xvalue(i_q)   ; y=yvalue(j_q)   ; z=zvalue(k_q)
xd1=x_q-x       ; yd1=y_q-y       ; zd1=z_q-z

if(xd1==0.d0) iloc(1)=0
if(yd1==0.d0) iloc(2)=0
if(zd1==0.d0) iloc(3)=0

do i=0,1
	do j=0,1
		do k=0,1
			loc_q(k*4+j*2+i+1,:)=(/i_q+i,j_q+j,k_q+k/)
			corlocq(k*4+j*2+i+1,:)=(/xvalue(i_q+i),yvalue(j_q+j),zvalue(k_q+k)/)
		enddo
	enddo
enddo

if(sum(iloc(:))==0)then
    charge(1)=1.d0
elseif(sum(iloc(:))==3)then
	do i=0,1
		do j=0,1
			do k=0,1
				xd=i*deltax-xd1
				yd=j*deltay-yd1
				zd=k*deltaz-zd1
				charge(4*k+2*j+i+1)=1.d0/abs(xd*yd*zd)
			enddo
		enddo
	enddo
elseif(sum(iloc(:))==2)then
    if(iloc(1)==0)then
		do j=0,1
			do k=0,1
				yd=j*deltay-yd1
				zd=k*deltaz-zd1                   !  weihua found a bug.
				charge(j*2+4*k+1)=1.d0/abs(yd*zd)
			enddo
		enddo
	elseif(iloc(2)==0)then
	    do i=0,1
		    do k=0,1
			    xd=i*deltax-xd1
				zd=k*deltaz-zd1
				charge(4*k+i+1)=1.d0/abs(xd*zd)
			enddo
		enddo
	elseif(iloc(3)==0)then
	    do i=0,1
		    do j=0,1
			    xd=i*deltax-xd1
				yd=j*deltay-yd1
				charge(i+2*j+1)=1.d0/abs(xd*yd)
			enddo
		enddo
	endif
elseif(sum(iloc(:))==1)then
    if(iloc(1)==1)then
	    charge(1)=1.d0/xd1
		charge(2)=1.d0/(deltax-xd1)
	elseif(iloc(2)==1)then
	    charge(1)=1.d0/yd1
		charge(3)=1.d0/(deltay-yd1)
    elseif(iloc(3)==1)then
	    charge(1)=1.d0/zd1
		charge(5)=1.d0/(deltaz-zd1)
	endif
endif

asum=sum(charge(:))
charge=q_q*charge/asum
corlocqt(iatm,:,:)=corlocq(:,:)
charget(iatm,:)=charge(:)
loc_qt(iatm,:,:)=loc_q(:,:)

return
end

!******************************************************
function inverx(x)
use mod_3d
use comdata
implicit double precision(a-h,o-z)

inverx=int((x-xleft)/deltax)+1

return
end
!
!-------------------------------------------
!
function invery(y)
use mod_3d
use comdata
implicit double precision(a-h,o-z)

invery=int((y-yleft)/deltay)+1

return
end
!
!-------------------------------------------
!
function inverz(z)
use mod_3d
use comdata
implicit double precision(a-h,o-z)

inverz=int((z-zleft)/deltaz)+1

return
end