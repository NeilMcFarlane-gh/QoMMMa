SUBROUTINE det_sm_coord()
use nrtype ; use coordinates; use optimdata
implicit none

! Using the data in coordinates, determines the full geometry of the small system.

integer(i4b) :: i, j, k, ii, jj, img_num
real(sp) :: vec(3)

! Loop over all images
do img_num=1,nimg
	! Check if update needed
	if (update_geom(img_num)) then
		! First, the easy bit: copy relevant bits of x to xq.
		! Copy coordinates of image first
		x(:)=fullx(img_num,:)

		do i=1,nq
			 j=3*(i-1)+1
			 k=3*(qm(i)-1)+1
			 xq(j:j+2)=x(k:k+2)
		end do

		! Copies coordinats of QM atoms into full QM atoms coordinate array
		fullxq(img_num,:)=xq(:)

		! then the slightly more difficult bit.
		! Rem each link H atom is placed along the QM-MM axis, at a distance of
		!       rQM-MM * (rQM-H_0 / rQM-MM_0)
		!  where rQM-H_0 and rQM-MM_0 are the "ideal" QM-L and QM-MM bond lengths.
		if (nl .ne. 0) then
		  do i=1,nl
			 j=3*(links(i,1)-1)+1
			 k=3*(links(i,2)-1)+1
			 vec = lratio(i) * (x(j:j+2)-x(k:k+2))
			 j=3*(i-1)+1
			 xl(j:j+2)=x(k:k+2)+vec
		  end do
		  
		! Copies coordinates of link atoms into full link atoms coordinate array
		  fullxl(img_num,:)=xl(:)
		end if
	end if
! End of loop over all images.
end do
return

END SUBROUTINE det_sm_coord


