SUBROUTINE update_step_geometry()
use nrtype ; use coordinates ; use optimdata
implicit none
                                                                                                                                             
integer(i4b) :: i, j, k, img_num
character(80) :: dummy

! Updates the geometry arrays with the QM region geometry following the update from a primitive step.
if (disp_prim .gt. 0) then
	! This is done in a fairly crude, but most simple way. "geom1.xyz" is simply overwritten.
	open(unit=8,file=("geom1.xyz"),status="replace")
	write(unit=8,fmt='(I6,2X,A)') n, "Title"
	do i=1,n
		 j = 3*(i-1)+1
		 write(unit=8,fmt='(I6,2x,A3,3F12.6,X,I5,10I6)') i,label(i), (x(k),k=j,j+2), &
			   & attyp(i),(bonds(i,k),k=1,nbonds(i))
	     !print *, i, label(i), (x(k),k=j,j+2), attyp(i), (bonds(i,k),k=1,nbonds(i))
	end do
	close(8)

end if
                                                                                                                                           
END SUBROUTINE update_step_geometry



