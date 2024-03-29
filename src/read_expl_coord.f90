SUBROUTINE read_expl_coord()
use nrtype ; use coordinates
implicit none

! This subroutine reads the "geom_expl.xyz" augmented tinker coordinates file
! It can only work after "alloc_coord" - so one needs to check that the arrays are allocated.

! Variable pointing on a specific image
integer(i4b) :: img_num
integer(i4b) :: rstat
integer(i4b) :: i, j, k, ii

! Loop generating sequential strings for file names

do img_num=1,nimg
	write (img_string(img_num),fmt='(i4)') img_num
	img_string(img_num)=adjustl(img_string(img_num))
end do


! Loop over all images
do img_num=1,nimg

	! First check whether the arrays are in fact already allocated...
	IF (.not.(allocated(bonds).and.allocated(attyp) &
	   & .and.allocated(x).and.allocated(nbonds).and.allocated(label))) THEN
		  write (*,*) "Arrays not yet allocated in read_expl_coord. ERROR."
		  stop
	END IF

	! Clear single-image coordinate table
	x=0.d0


	! Then open groups the "geom.xyz" Tinker file
	! and check that the number of atoms is correct
	open(unit=8,file=("geom_expl"//trim(img_string(img_num))//".xyz"), position="rewind")
	read(unit=8,fmt=*) i

	if (i.ne.n) THEN
		  write (*,*) "Mismatch in number of atoms in read_expl_coord. ERROR."
		  stop
	end if

	nbonds = 0
	bonds = 0
	do i=1,n
		 read(unit=8,fmt=*) nbonds(i)
		 j = 3*(i-1)+1
		 if (nbonds(i) .eq. 0) then
			  read(unit=8,fmt=*) ii,label(i),(x(k),k=j,j+2),attyp(i)
		 else
			  read(unit=8,fmt=*) ii,label(i),(x(k),k=j,j+2),attyp(i), &
				   & (bonds(i,k),k=1,nbonds(i))
		 end if
		 if (i.ne.ii) then
			   write (*,*) "Atom numbering wrong in fullgeom. ERROR."
			   stop
		 end if
	end do
	
	close(8)

	! Copies coorindates of the image into the main array of coordinates
	fullx(img_num,:)=x(:)

! End of loop over all images
end do

END SUBROUTINE read_expl_coord


