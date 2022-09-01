SUBROUTINE write_resp()
use nrtype ; use coordinates ; use optimdata
implicit none

! This subroutine writes a set of coordinates (and weights) for computation of the
! electrostatic potential for each image. These points are at six positions for 
! each MM atom (except the link atoms), slightly displaced in both directions 
! along each of the x, y and z axes.

integer(i4b) :: i, j, k, ii, img_num
real(dp) :: vec0(3), vec(3)
real(dp), parameter :: dx = bohr * .5d-3, eps = 1.d-4

! Start loop over all images
do img_num=1,nimg
	! Update only if necessary
	if (update_geom(img_num)) then
		! Copies charges, coordinates for the image from all-image table into single-image table
		x(:)=fullx(img_num,:)

		! Then Start to Write the file
		! and check that the number of the various atom types is correct
		open(unit=9,file=("points"//trim(img_string(img_num))//".pts"),status="replace")
		open(unit=11,file=("points_chg"//trim(img_string(img_num))//".pts"),status="replace")
		open(unit=12,file=("points"//trim(img_string(img_num))//"gau.pts"),status="replace")

		atoms: do i=1,n
			 do j = 1, nq
				 if (qm(j) .eq. i) cycle atoms
			 end do
			 do j = 1, nl
				 if (links(j,1) .eq. i) cycle atoms
			 end do
			 if (inact(i)) cycle atoms
			 if (abs(chg(i)).lt.eps) cycle atoms
			 vec0 = x(3*(i-1)+1:3*i)
			 vec = vec0
			 write(11,*) chg(i)
			 write(12,'(3F20.6)') vec
			 do j = 1, 3
				  vec(j) = vec0(j) + dx
				  write(unit=9,fmt='(3F12.6,F5.1)') vec, 0.d0
				  vec(j) = vec0(j) - dx
				  write(unit=9,fmt='(3F12.6,F5.1)') vec, 0.d0
				  vec(j) = vec0(j)
			 end do
		end do atoms

		close(9)
		close(11)
		close(12)

	end if
	
! End of loop over all images
end do
END SUBROUTINE write_resp


