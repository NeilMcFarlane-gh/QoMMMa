SUBROUTINE read_hessgeom_simple()
use nrtype ; use coordinates ; use optimdata
implicit none

! This subroutine reads an xyz block as output by the Hessian optimization cycle.

integer(i4b) :: i, j, k, ii, img_num

! Then Start to Write the file
! and check that the number of the various atom types is correct

open(unit=9,file="nhessgeom.xyz")

do i=1,nopt
     j = 3*(i-1)+1
     read(unit=9,fmt=*) xopt(j:j+2)
end do

close(9)

END SUBROUTINE read_hessgeom_simple


