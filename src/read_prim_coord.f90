SUBROUTINE read_prim_coord()
use nrtype ; use coordinates ; use optimdata
implicit none

! This subroutine reads the file "prim_list".
! Only used if primitive internal coordinates are used.

integer(i4b) :: rstat, i, j, k, ii, img_num
integer(i4b) :: work(nprim,4)
character(80) :: dummy

if (coordtype .eq. 1) then
	! Write the primitive internal coordinate definitions to "prim_list".
	
	open(unit=7,file="prim_list",status="replace")
	read(unit=7,fmt=*) nprim

	do i=1,nprim
		read(unit=7,fmt=*) ( prim_list(i,k),k=1,4 )
	end do

	close(7)

end if

END SUBROUTINE read_prim_coord


