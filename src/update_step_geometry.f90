SUBROUTINE update_step_geometry()
use nrtype ; use coordinates ; use optimdata
implicit none
                                                                                                                                             
integer(i4b) :: i, j, k, img_num
character(80) :: dummy

! Updates the geometry arrays with the QM region geometry following the update from a primitive step.
if (ncon_prim .gt. 0) then
	! First, copy over the new cartesian coordinates to the array containing the whole coordinates.
	do i = 1, nopt
		j = 3*i-2
		k = 3*opt(i)-2
		x(k:k+2) = newx(j:j+2)
	end do
	fullxopt(img_num,:)=xopt(:)
	fullx(img_num,:)=x(:)

	! This is done in a fairly crude, but most simple way. "geom_expl.xyz" is simply overwritten.
	open(unit=8,file=("geom_expl"//trim(img_string(img_num))//".xyz"),status="replace")
	write(unit=8,fmt='(I6,2X,A)') n, "Title"
	do i=1,n
		 j = 3*(i-1)+1
		 write(unit=8,fmt='(I6,2x,A3,3F12.6,X,I5,10I6)') i,label(i), (x(k),k=j,j+2), &
			   & attyp(i),(bonds(i,k),k=1,nbonds(i))
	end do
	close(8)
end if
                                                                                                                                             
END SUBROUTINE update_step_geometry



