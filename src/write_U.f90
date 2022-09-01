SUBROUTINE write_U
use nrtype ; use coordinates ; use optimdata
implicit none

! This subroutine write the U matrix to a file: 'Umat'.

integer(i4b) :: rstat, i, j, k, ii, img_num
integer(i4b) :: work(nprim,ndlc)
character(80) :: dummy

if (coordtype .eq. 1) then
	! Write the U matrix to the file 'Umat'.
	
	open(unit=6,file="Umat",status="replace")
	write(unit=6,fmt=*) nprim, ndlc

	do i=1,nprim
		 write(unit=6,fmt=*) ( Umat(i,k),k=1,ndlc )
	end do

	close(6)

end if

END SUBROUTINE write_U


