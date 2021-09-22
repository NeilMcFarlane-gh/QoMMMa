SUBROUTINE write_atomic()
use nrtype ; use coordinates
implicit none

! This subroutine writes an "&atomic" group suitable for inclusion in Jaguar.
! It is QM-code SPECIFIC!!

integer(i4b) :: i, j, k, ii

! Then Start to Write the file
! and check that the number of the various atom types is correct
! Unchanged for NEB - the same atomic section for all images

open(unit=9,file="atomic_prelim",status="replace")

do i=1,nq
     write(unit=9,fmt='(A2)') qlabel(i)
end do
do i=1,nl
     write(unit=9,fmt='(A2)') llabel(i)
end do

close(9)

END SUBROUTINE write_atomic


