SUBROUTINE write_hessgeom()
use nrtype ; use coordinates ; use optimdata
implicit none

! This subroutine writes an xyz block of the geometry included in the Hessian Optimization.

integer(i4b) :: i, j, k, img_num

! Loop over all images
do img_num=1,nimg
	if (update_geom(img_num)) then
	xopt(:)=fullxopt(img_num,:)

	open(unit=9,file=("nhessgeom"//trim(img_string(img_num))//".xyz"),status="replace")

	do i=1,nopt
		 j = 3*(i-1)+1
		 write(unit=9,fmt='(3F14.8)') (xopt(k),k=j,j+2)
	end do

	close(9)

	end if
	
! End of loop over images
end do

END SUBROUTINE write_hessgeom


