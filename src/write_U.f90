SUBROUTINE write_U
use nrtype ; use coordinates ; use optimdata
implicit none

! This subroutine write the U matrix to a file: 'Umat'.

integer(i4b) :: rstat, i, j, k, ii, img_num
integer(i4b) :: work(nprim,ndlc)
character(80) :: dummy

if (coordtype .eq. 1) then
	! Write the U matrix to the file 'Umat'.
	open(unit=60,file="Umat",status="replace")
	write(unit=60,fmt=*) nprim, ndlc
	do i=1,nprim
		if (ncon_prim .eq. 0) then
			write(unit=60,fmt=*) ( Umat(i,k),k=1,ndlc )
		else if (ncon_prim .gt. 0) then
			write(unit=60,fmt=*) ( Vmat(i,k),k=1,ndlc )
		end if
	end do
	close(60)

end if

END SUBROUTINE write_U


