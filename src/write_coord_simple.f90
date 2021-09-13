SUBROUTINE write_coord_simple()
use nrtype ; use coordinates ; use optimdata
implicit none

! This subroutine writes the "ngeom.xyz" tinker coordinates file
! It can only work after "alloc_coord" - so one needs to check that the arrays are allocated.

integer(i4b) :: rstat, i, j, k, ii, img_num
character(80) :: dummy

! First check whether the arrays are in fact already allocated...

IF (.not.(allocated(qm).and.allocated(bonds).and.allocated(attyp) &
   & .and.allocated(x).and.allocated(nbonds).and.allocated(label))) THEN
      write (*,*) "Arrays not yet allocated in write_coord. ERROR."
      stop
END IF

! Then Start to dump a new "FullGeom" file
! and check that the number of the various atom types is correct

open(unit=8,file="ngeom.xyz",status="replace")
write(unit=8,fmt='(I6,2X,A)') n, "Title"

do i=1,n
     j = 3*(i-1)+1
     write(unit=8,fmt='(I6,2x,A3,3F12.6,X,I5,10I6)') i,label(i), (x(k),k=j,j+2), &
           & attyp(i),(bonds(i,k),k=1,nbonds(i))
end do

close(8)

END SUBROUTINE write_coord_simple


