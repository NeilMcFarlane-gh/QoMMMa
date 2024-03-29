SUBROUTINE write_first_checkfile()
use nrtype ; use coordinates ; use optimdata
implicit none

! This subroutine writes the "CheckFile" for the first time.
! It can only work after "alloc_coord" - so one needs to check that the arrays are allocated.

integer(i4b) :: i, j, k, ii, img_num
character(80) :: dummy

! First check whether the arrays are in fact already allocated...

IF (.not.(allocated(qm).and.allocated(bonds).and.allocated(attyp) &
   & .and.allocated(x).and.allocated(nbonds).and.allocated(label))) THEN
      write (*,*) "Arrays not yet allocated in write_coord. ERROR."
      stop
END IF

! Start loop over all images
do img_num=1,nimg
	! Then Start to dump a new "FullGeom" file
	! and check that the number of the various atom types is correct
	open(unit=8,file=("CheckFile"//trim(img_string(img_num))),status="replace")
	write(unit=8,fmt='(A)') "This file contains info needed for execution of qmmm.x"
	write(unit=8,fmt='(A)') "First, the number of steps already taken:"
	write(unit=8,fmt='(I3)') 0
	write(unit=8,fmt='(A)') "Then, the total, QM, Link, and Hessian-optimized number of atoms:"
	write(unit=8,fmt='(5I6)') n, nq, nl, nopt
	write(unit=8,fmt='(A)') "Do we want dispersion energy between QM atoms to be calculated?"
	write(unit=8,fmt='(I6)') disp
	write(unit=8,fmt='(A)') "Then the number of images,type of calculation, spring force constant"
	write(unit=8,fmt='(2I6, F10.2)') nimg, nebtype, kspring
	write(unit=8,fmt='(A)') "Then the type of growing string method which is being used"
	write(unit=8,fmt='(2I6, F10.2)') gsmtype
	write(unit=8,fmt='(A)') "Then the current phase of the growing string method"
	write(unit=8,fmt='(2I6, F10.2)') gsmphase
	write(unit=8,fmt='(A)') "Then the coordinate selection for optimisation"
	write(unit=8,fmt='(2I6, F10.2)') coordtype
	write(unit=8,fmt='(A)') "Then the trust radius for optimisation"
	write(unit=8,fmt='(F10.2)') trust_radius
	write(unit=8,fmt='(A)') "If DLC are used, the type of primitive internal coordinates used to generate the DLC"
	write(unit=8,fmt='(2I6, F10.2)') primtype
	write(unit=8,fmt='(A)') "Then the number of cartesian constraints and the type (1=Harmonic or 2=tanh)"
	write(unit=8,fmt='(2I6)') ncon_cart, kcnstyp
	write(unit=8,fmt='(A)') "Then the number of primitive internal coordinate constraints"
	write(unit=8,fmt='(2I6)') ncon_prim
	if (coordtype .eq. 1) then
		write(unit=8,fmt='(A)') "Then the number of primitive internal coordinates and their definition."
		write(unit=8,fmt='(I6)') nprim
		do i = 1, nprim
			write (unit=8,fmt='(4I6,3X)') (prim_list(i,j),j=1,4)
		end do
		write(unit=8,fmt='(A)') "Then the number of delocalised internal coordinates and their transformation matrix."
		write(unit=8,fmt='(I6)') ndlc
		do i = 1, nprim
			if (ncon_prim .eq. 0) then
				write (unit=8,fmt=*) (Umat(i,j),j=1,ndlc)
			else if (ncon_prim .gt. 0) then
				write (unit=8,fmt=*) (Vmat(i,j),j=1,ndlc)
			end if
		end do
	end if
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
	
	if (coordtype .eq. 0) then
		write(unit=8,fmt='(A)') "Cartesian constraint Details. For each constraint, first the nature &
			&  of the constrained coordinate:"
		write(unit=8,fmt='(A)') "(1 = r(A-B), 2 = r(A-B) - r(C-D), 3 = r(A-B) + r(C-D),&
			& 4 = r(A-B) + r(C-D) - r(E-F))"
		write(unit=8,fmt='(A)') "Then the force constant in kcal/mol / Angstrom^2, and the ideal value."
		write(unit=8,fmt='(A)') "Then (in a new line) the atoms to which the constraint&
			& should apply (A, B, (C, D, (E, F)))"
		if (ncon_cart.gt.0) then
		   do i=1,ncon_cart
				write(unit=8,fmt='(I6,F10.2,F10.4)') cnstyp(i),kcns(i),cnsidl(i)
				write(unit=8,fmt='(13I6)') (cnsat(i,j),j=1,ncnsat(i))
		   end do
		end if
	else if (coordtype .eq. 1) then
		write(unit=8,fmt='(A)') "Primitive internal coordinate constraint Details. For each constraint, first the change, if necessary."
		write(unit=8,fmt='(A)') "Then (in a new line) the atoms to which the constraint should apply."
		if (ncon_prim.gt.0) then
		   do i=1,ncon_prim
				write(unit=8,fmt='(A,I6,A)') "The number of shaling coefficients for constraint number ", i, ' is...'
				write(unit=8,fmt='(I4)') cns_n_coeff_p(i)
				write(unit=8,fmt='(A)') "Now what are the details of the primitive internal coordinate linear combination constraint..."
				do j=1,cns_n_coeff_p(i)
					write(unit=8,fmt='(F10.2,4I6)') cnscoeff_p(i,j), (cnsat_p(i,j,k),k=1,4)
				end do
		   end do
		end if
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

	close(8)

! End of loop over all images
end do

END SUBROUTINE write_first_checkfile


