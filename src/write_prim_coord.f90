SUBROUTINE write_prim_coord()
use nrtype ; use coordinates ; use optimdata
implicit none

! This subroutine writes two files: "prims" and "prim_list".
! Only used if primitive internal coordinates are used.

integer(i4b) :: rstat, i, j, k, ii, img_num
integer(i4b) :: work(nprim,4)
character(80) :: dummy

if (coordtype .eq. 1) then
	! Write the primitive coordinates themselves to the file "prims".
	! Also, write the primitive internal coordinate definitions to "prim_list".
	
	open(unit=6,file="prims",status="replace")
	open(unit=7,file="prim_list",status="replace")
	write(unit=6,fmt=*) nprim
	write(unit=7,fmt=*) nprim

	do i=1,nprim
		 write(unit=6,fmt='(F12.6)') prims(i)
		 write(unit=7,fmt='(4I6)') ( prim_list(i,k),k=1,4 )
	end do

	close(6)
	close(7)

end if

END SUBROUTINE write_prim_coord


