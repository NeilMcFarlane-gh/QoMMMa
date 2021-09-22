subroutine write_mulliken()
use nrtype ; use coordinates ; use optimdata
implicit none

! writes lines to be included in the tinker "key" file for the microiteration optimization.
! Sets the QM (and Link) atoms as "inactive", and gives the Mulliken Charges for the QM atoms.
! Also writes any modified charges for MM atoms.

integer(i4b) :: nql, i, j, k, kk, ii(ninact), img_num

! Loop over all images
do img_num=1,nimg

	! Check whether update necessary
	if (update_geom(img_num)) then

	mull(:)=fullmull(img_num,:)
	open(unit=8,file=("inactive_microiter"//trim(img_string(img_num))),status="replace")

	write (8,*) "EXTRATERM"
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
	kk=nopt-nq-nl
	if (kk.gt.0) then
	   nql=kk/10
	   do i=1,nql
		  k=10*(i-1)+1
		  write(8,'(A,10I6)') "INACTIVE ",(opt(nq+nl+j),j=k,k+9)
	   end do
	   k=10*nql+1
	   if (k .le. kk) then
		  write(8,'(A,10I6)') "INACTIVE ",(opt(nq+nl+j),j=k,kk)
	   end if
	end if

	! Print list of inactive MM atoms.
	if (ninact.gt.0) then
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
		   write (8,'(A,10I6)') "INACTIVE ",(ii(j),j=k,ninact)
	   end if
	end if

	do i=1,nq
	   write (unit=8,fmt='(A,I8,F12.5)') "charge",-qm(i), mull(i)
	end do

	outer: do i=1,n
	   do j = 1, nq
		   if(i.eq.qm(j)) cycle outer
	   end do
	   if (modchg(i)) then
		   write (unit=8,fmt='(A,I8,F8.3)') "charge",-i, chg(i)
	   end if
	end do outer

	write (unit=8,fmt=*) ""

	close(8)

	end if

! End of loop over images
end do

end subroutine write_mulliken



