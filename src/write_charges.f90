SUBROUTINE write_charges()
use nrtype ; use coordinates ; use optimdata
implicit none

! This subroutine writes an block of charges suitable for inclusion in Jaguar(Gaussian/etc.) input.
! This subroutine is QM-program SPECIFIC as yet.

integer(i4b) :: i, j, k, ii, img_num
real(sp) :: eps=1d-4

!Loop over all images
do img_num=1,nimg
	if (update_geom(img_num)) then
		! Copies charges from all-image array into single-image array
		x(:)=fullx(img_num,:)
		
		! Start to Write the file
		! and check that the number of the various atom types is correct
		open(unit=9,file=("charges"//trim(img_string(img_num))//".xyz"),status="replace")
		open(unit=10,file=("charges"//trim(img_string(img_num))//"mol.xyz"),status="replace")

		write (unit=9,fmt='(A)') "&pointch"
		write (unit=10,fmt='(A)') "&pointch"
		atoms: do i=1,n
			 do j = 1, nq
				 if (qm(j) .eq. i) cycle atoms
			 end do
			 do j = 1, nl
				 if (links(j,1) .eq. i) cycle atoms
			 end do
			 if (abs(chg(i)).lt.eps) cycle atoms
			 ii = 3*(i-1)+1
			 write(unit=9,fmt='(F8.4,2X,3F12.6)') chg(i),(x(k),k=ii,ii+2)
			 if (inact(i)) then
			   write(unit=10,fmt='(F8.4,2X,3F12.6,I5)') chg(i),(x(k),k=ii,ii+2),0
			 else
			   write(unit=10,fmt='(F8.4,2X,3F12.6,I5)') chg(i),(x(k),k=ii,ii+2),1
			 end if
		end do atoms
		write (unit=9,fmt='(A)') "&"
		write (unit=10,fmt='(A)') "&"

		close(9)
		close(10)

	end if
! End of the loop over all images
end do

END SUBROUTINE write_charges


