SUBROUTINE write_sm_coord()
use nrtype ; use coordinates ; use optimdata
implicit none

! This subroutine writes an xyz block suitable for inclusion in Jaguar(Gaussian/etc.) input.

integer(i4b) :: i, j, k, ii, img_num

! Loop over all images
do img_num=1,nimg
	
	! Check whether update necessary
	if (update_geom(img_num)) then
	! Copies coordinates of QM and link atoms into single-image array
	!xq(:)=fullxq(img_num,:)
	!xl(:)=fullxl(img_num,:)

	! Then Start to Write the file
	! and check that the number of the various atom types is correct

	open(unit=90,file=("nqmgeom"//trim(img_string(img_num))//".xyz"),status="replace")
	
	do i=1,nq
		 j = 3*(i-1)+1
		 write(unit=90,fmt='(A2,2X,3F12.6)') qlabel(i),xq(j),xq(j+1),xq(j+2)
	end do
	do i=1,nl
		 j=3*(i-1)+1
		 write(unit=90,fmt='(A2,2X,3F12.6)') llabel(i),xl(j),xl(j+1),xl(j+2)
	end do

	close(90)

	end if

! End of loop over all images
end do

END SUBROUTINE write_sm_coord


