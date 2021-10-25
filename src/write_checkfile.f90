SUBROUTINE write_checkfile()
use nrtype ; use coordinates ; use optimdata
implicit none

! This subroutine writes the "CheckFile"
! It can only work after "alloc_coord" - so one needs to check that the arrays are allocated.

integer(i4b) :: i, j, k, ii, img_num
character(80) :: dummy

! First check whether the arrays are in fact already allocated...

IF (.not.(allocated(qm).and.allocated(bonds).and.allocated(attyp) &
   & .and.allocated(x).and.allocated(nbonds).and.allocated(label))) THEN
      write (*,*) "Arrays not yet allocated in write_coord. ERROR."
      stop
END IF

! Loop over all images
do img_num=1,nimg

	! First copy data for the image

	! Use old values for images that are not to be updated
	if (update_geom(img_num)) then
		xopt(:)=fullxopt(img_num,:)
		optg(:)=fulloptg(img_num,:)
		e=fulle(img_num)
		h(:,:)=fullh(img_num,:,:)
	else
		xopt(:)=fullox(img_num,:)
		optg(:)=fullog(img_num,:)
		!Will store a new energy, for convergence test to be zero
		!e=fulloe(img_num)
		e=fulle(img_num)
		h(:,:)=fulloh(img_num,:,:)
	end if

	if ((line_search).and.(nstep.ne.1)) then
		 xopt(:)=fullox(img_num,:)
		 optg(:)=fullog(img_num,:)
		 e=fulloe(img_num)
		 h(:,:)=fulloh(img_num,:,:)
		 if (nstep.eq.1) then
			 e=fulle(img_num)
		 end if
	end if

	if ((resetting).and.(nstep.eq.2)) then
		xopt(:)=fullxopt(img_num,:)
		optg(:)=fulloptg(img_num,:)
		e=fulle(img_num)
		h(:,:)=fullh(img_num,:,:)
	end if

	! Then Start to dump a new "FullGeom" file
	! and check that the number of the various atom types is correct

	open(unit=8,file=("nCheckFile"//trim(img_string(img_num))),status="replace")
	write(unit=8,fmt='(A)') "This file contains info needed for execution of qommma.x"
	write(unit=8,fmt='(A)') "First, the number of steps already taken:"
	write(unit=8,fmt='(I3)') nstep
	write(unit=8,fmt='(A)') "Then, the total, QM, Link, and Hessian-optimized number of atoms:"
	write(unit=8,fmt='(4I6)') n, nq, nl, nopt
	write(unit=8,fmt='(A)') "Do we want dispersion energy between QM atoms to be calculated?"
	write(unit=8,fmt='(I6)') disp
	if ((nebtype.eq.4).or.(nebtype.eq.3).or.(nebtype.eq.5)) then
		write(unit=8,fmt='(A)') "Then the number of images,type of calculation, spring force constant, scaling factor"
		write(unit=8,fmt='(2I6, F10.2, F10.5)') nimg, nebtype, kspring*hart_kcal, weight
	else
		write(unit=8,fmt='(A)') "Then the number of images,type of calculation, spring force constant"
		write(unit=8,fmt='(2I6, F10.2)') nimg, nebtype, kspring*hart_kcal
	end if
	write(unit=8,fmt='(A)') "Then the type of growing string method which is to be used"
	write(unit=8,fmt='(2I6, F10.2)') gsmtype
	write(unit=8,fmt='(A)') "Then the coordinate selection for optimisation"
	write(unit=8,fmt='(2I6, F10.2)') coordtype
	if (coordtype .eq. 1) then
		write(unit=8,fmt='(A)') "Then the number of primitive internal coordinates and their definition."
		write(unit=8,fmt='(2I6, F10.2)') nprim
		do i = 1, nprim
			write (unit=8,fmt='(2I6,3X)') (prim_list(i,j),j=1,2)
		end do
	end if
	write(unit=8,fmt='(A)') "Then the number of constraints and the type (1=Harmonic or 2=tanh)"
	write(unit=8,fmt='(2I6)') ncon, kcnstyp
	write(unit=8,fmt='(A)') "Then the list of which atoms are QM:"
	do i = 1, nq
		 write (unit=8,fmt='(I6,3X,A2)') qm(i),qlabel(i)
	end do
	write(unit=8,fmt='(A)') "Link atom details:"
	if (nl .gt. 0) then
		 do i = 1, nl
			  write (unit=8,fmt='(2I6,F12.6,2X,A2)') (links(i,j),j=1,2), lratio(i), llabel(i)
		 end do
	end if
	write(unit=8,fmt='(A)') "Then a list of all non-QM and non-link atoms included in Hessian optimization:"
	j = nopt - nq - nl
	if (j .gt. 0) then
		 do i = 1, j
			  write (unit=8,fmt='(I6)') opt(nq+nl+i)
		 end do
	end if

	write(unit=8,fmt='(A)') "Constraint Details. For each constraint, first the nature &
			&  of the constrained coordinate:"
	write(unit=8,fmt='(A)') "(1 = r(A-B), 2 = r(A-B) - r(C-D), 3 = r(A-B) + r(C-D),&
			& 4 = r(A-B) + r(C-D) - r(E-F))"
	write(unit=8,fmt='(A)') "Then the force constant in kcal/mol / Angstrom^2 and the ideal value."
	write(unit=8,fmt='(A)') "Then (in a new line) the Atoms to which the constraint&
			& should apply (A, B, (C, D, (E, F)))"
	if (ncon.gt.0) then
	   do i=1,ncon
			write(unit=8,fmt='(I6,F10.2,F10.4)') cnstyp(i),kcns(i),cnsidl(i)
			write(unit=8,fmt='(13I6)') (cnsat(i,j),j=1,ncnsat(i))
	   end do
	end if

	write(unit=8,fmt='(A)') "Atomic Charges (on all atoms - not used on QM atoms! - flag for non-standard)"
	do i = 1, n
		 ii = 0
		 if (modchg(i)) ii = 1
		 write (unit=8,fmt='(I6,F12.6,1X,I1)') i, chg(i), ii
	end do
	write(unit=8,fmt='(A)') "List of MM atoms to be held frozen (inactive). First, number of such atoms:"
	write (unit=8,fmt='(I8)') ninact
	write(unit=8,fmt='(A)') "Then the list of such atoms:"
	do i = 1, n
		 if (inact(i)) then
				write (unit=8,fmt='(I8)') i
		 end if
	end do

	write(unit=8,fmt='(A)') "Energy at the previous step:"
	write(unit=8,fmt='(F15.8)') e
	write(unit=8,fmt='(A)') "Geometry of the Hessian Part at the previous step:"
	write(unit=8,fmt='(3F12.6)') (xopt(i),i=1,noptx)
	write(unit=8,fmt='(A)') "Gradient at the previous step:"
	write(unit=8,fmt='(3F12.6)') (optg(i),i=1,noptx)
	write(unit=8,fmt='(A)') "Approximate Inverse Hessian at the previous step (lower half!):"
	if (coordtype .eq. 0) then
		do i = 1, noptx
			write (unit=8,fmt='(5F20.6)') (h(i,j),j=1,i)
		end do
	else if (coordtype .eq. 1) then
		do i = 1, nprim
			write (unit=8,fmt='(5F20.6)') (h_p(i,j),j=1,i)
		end do
	end if

	close(8)

! End of loop over all images
end do

END SUBROUTINE write_checkfile


