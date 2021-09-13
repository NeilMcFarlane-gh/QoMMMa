SUBROUTINE read_hessgeom()
use nrtype ; use coordinates ; use optimdata
implicit none

! This subroutine reads an xyz block as output by the Hessian optimization cycle.

integer(i4b) :: i, j, k, ii, img_num

! Loop over all images
do img_num=1,nimg
if (update_geom(img_num)) then
! Then Start to Write the file
! and check that the number of the various atom types is correct

open(unit=9,file=("nhessgeom"//trim(img_string(img_num))//".xyz"))

do i=1,nopt
     j = 3*(i-1)+1
     read(unit=9,fmt=*) xopt(j:j+2)
end do

close(9)

fullxopt(img_num,:)=xopt(:)

end if
! End of loop over all images
end do
END SUBROUTINE read_hessgeom


