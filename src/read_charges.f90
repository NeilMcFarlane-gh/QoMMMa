SUBROUTINE read_charges()
use nrtype ; use coordinates
implicit none

! This subroutine reads a file produced by "analyze" containing the default points
!  charges located on each atom (including the QM ones). This array is modified thereafter
!  by cancelling all charges on QM atoms (automatically) and changing charges on L and connected
!  MM atoms (by reading in a list of changes from InitFile).

character(80) :: dummy
integer(i4b) :: i, j, k, rstat
real(dp) :: tc

open(unit=8,file="DefaultCharges")

do i=1,n
     read(unit=8,fmt=*,IOSTAT=rstat) j,chg(i)
     if (rstat .ne. 0) then
         write (*,*) "error reading DefaultCharges. ERROR."
         close(8)
         stop
     end if
     if (j .ne. i) then
         write (*,*) "problem with atom numbering in DefaultCharges. ERROR."
         close(8)
         stop
     end if
end do
close(8)


END SUBROUTINE read_charges


