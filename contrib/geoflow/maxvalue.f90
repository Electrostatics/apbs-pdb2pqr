!
!
	program READMAX
	implicit none
	real*8 ::	surfu,max,min
	integer :: ix, iy, iz
	integer :: n, ind

	write(*,*) "enter number of values to read: "
	read(*,*) n
	min=0.d0
	max=0.d0
	do ind =1,n
		read(21,*) surfu
		if (surfu > max) then
			max = surfu
		endif
		
		if(surfu < min) then
			min = surfu
		endif
	enddo
	
!	write(*,*) ind
!	write(*,*) "max = ",max
!	write(*,*) "min = ",min

	end
	
