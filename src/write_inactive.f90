subroutine write_inactive()
use nrtype ; use coordinates
implicit none

! Writes lines to be included in the tinker "key" file
! which set the QM atoms as "inactive" and zeroes their charges.
! Also writes any modified charges for MM atoms.
! Only one inactive list created, as inactives in all images must be the same

integer(i4b) :: nql, i, j, k, ii(ninact)

open(unit=8,file="inactive_list",status="replace")

write (unit=8,fmt=*) "EXTRATERM NONE"

nql=nq/10

do i=1,nql
     k=10*(i-1)+1
     write(8,'(A,10I6)') "INACTIVE ",(qm(j),j=k,k+9)
end do
k=10*nql+1
if (k .le. nq) then
   write (8,'(A,9I6)') "INACTIVE ",(qm(j),j=k,nq)
end if
if (nl .gt. 0) then
   nql=nl/10
   do i=1,nql
       k=10*(i-1)+1
       write (8,'(A,10I6)') "INACTIVE ",(links(j,1),j=k,k+9)
   end do
   k=10*nql+1
   if (k.le.nl) then
       write (8,'(A,9I6)') "INACTIVE ",(links(j,1),j=k,nl)
   end if
end if

! Print list of inactive MM atoms.

j=0
do i=1,n
   if (inact(i)) then
       j=j+1
       ii(j)=i
   end if
end do
nql=ninact/10
do i=1,nql
     k=10*(i-1)+1
     write(8,'(A,10I6)') "INACTIVE ",(ii(j),j=k,k+9)
end do
k=10*nql+1
if (k .le. ninact) then
    write (8,'(A,9I6)') "INACTIVE ",(ii(j),j=k,ninact)
end if

! Then write a list of modified charges.

do i=1,n
   if (modchg(i)) then
       write (unit=8,fmt='(A,I8,F8.3)') "charge",-i, chg(i)
   end if
end do

write (unit=8,fmt=*) ""

close(8)

end subroutine write_inactive



