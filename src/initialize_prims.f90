SUBROUTINE initialize_prims()
use nrtype ; use coordinates ; use optimdata ; use primitive ; use math
implicit none

! In this subroutine, the primitive internal coordinate set is initialised and is used throughout optimisation.
! This is to avoid discontinuity in the DLC, as the primitive internal coordinates are defined once at the start.
! With each evaluation of DLC, the actual values of the primitives are re-calculated, but their definition remains the same.
! If DLC are not used, then primitive internal coordinates are not defined as they are not required.

integer(i4b) :: i, j, k
	
if (coordtype .eq. 1) then
	! If the primitive internal coordinates have already been defined, then we can add any additional primitives and move on.
	if (ALLOCATED(prim_list)) then
		call add_primitives()
		return
	end if
	
	! The cartesian coordinates of the DLC region are used to define the prims in define_prims_full.
	! First, add the QM atoms, and then add the link atom coordinates.
	xopt = 0.d0
	do i=1, nq
		j = (3 * (i-1)) + 1
		k = (3 * (opt(i)-1)) + 1
		xopt(j:j+2) = x(k:k+2)
	end do
	j = (3 * (nq-1)) + 1
	do i=1, nl
		k = (3 * (i-1)) + 1
		j = j + 3
		xopt(j:j+2) = xl(k:k+2)
	end do

	! Now, the primitive internal coordinates can be defined.
	if (primtype .eq. 0) then ! Total connectivity
		call define_prims_TC(nopt, opt, prim_list)
	else if (primtype .eq. 1) then ! Full definition
		call define_prims_full(nopt, opt, xopt, prim_list)	
	end if
	
	!do i=1, SIZE(prim_list,1)
	!	print *, prim_list(i,:)
	!end do
	!stop
	
	! Lastly, add any additional primitive internal coordinates which have not been found automatically.
	call add_primitives()
end if

end subroutine initialize_prims


subroutine add_primitives()
use nrtype ; use coordinates ; use optimdata ; use primitive ; use math
implicit none

! Subroutine which adds the additional primitives defined to prim_list.
! If the additional primitive has been added already by the generation algorithm or was read in via fortinput, then it will not be added.

integer(i4b) :: i, j, to_add, prim_temp1(4), prim_temp2(4), prim_list_save(nprim,4)
logical :: need_to_add(add_prims)

! Find all the primitives which need to be included.
to_add = 0
do i=1, add_prims
		prim_temp1 = prim_add_list(i,:)
		do j=1, nprim
			prim_temp2 = prim_list(j,:)
			if (MAXVAL(ABS(prim_temp1 - prim_temp2)) .eq. 0) then ! the primtive internal coordinate is already there so no need to add it.
				need_to_add = .False.
				exit
			else if (j .eq. nprim) then
				to_add = to_add + 1
				need_to_add(i) = .True. ! we need to add a new one to prim_list.
			end if
		end do
end do

! Deallocate and subsequently re-allocate prim_list with space allowed for the new primitives to be included.
prim_list_save(:,:) = prim_list(:,:)
deallocate(prim_list)
allocate(prim_list(nprim + to_add, 4))

! Now add the original and additional primitived to prim_list
prim_list(1:nprim,:) = prim_list_save(:,:) ! the original primitives
to_add = 0
do i=1, add_prims ! the additional primitives
	if (need_to_add(i) .eqv. .True.) then
		to_add = to_add + 1
		prim_list(nprim + to_add,:) = prim_add_list(i,:)
	end if
end do
nprim = nprim + to_add ! update the number of primitives

end subroutine add_primitives