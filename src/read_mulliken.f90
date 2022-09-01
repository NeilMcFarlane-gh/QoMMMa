SUBROUTINE read_mulliken()
use nrtype ; use coordinates ; use optimdata
implicit none

! This subroutine reads the "mulliken" file containing the Mulliken charges
!   on the QM atoms. Charges on H Link atoms are added in to the charges on the
!   corresponding QM atoms.

integer(i4b) :: i, j, k, img_num
real(dp) :: chgt

!Loop over all images
do img_num=1,nimg
	open(unit=8,file=("mulliken"//trim(img_string(img_num))),position="rewind")

	! Clear single-image table of charges
	mull=0.d0

	! read in the charges on QM atoms

	do i = 1, nq
		read(unit=8,fmt=*) mull(i)
	end do

	! Read in the charges on link atoms; add to corresponding QM charge.

	do i = 1, nl
		k=0
		inner: do j = 1, nq
			if (qm(j).eq.links(i,2)) then
				k = j
				exit inner
			end if
		end do inner
		if (k.eq.0) then
			write (*,*) "not found k in read_mulliken. Error."
			stop
		end if
		read(unit=8,fmt=*) chgt
		mull(k)=mull(k)+chgt
	end do

	close(8)

	! Copies Mulliken charges into all-image table

	fullmull(img_num,:)=mull(:)

! End of loop over all images
end do
return

END SUBROUTINE read_mulliken


