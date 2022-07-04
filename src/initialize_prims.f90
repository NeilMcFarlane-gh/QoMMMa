SUBROUTINE initialize_prims()
use nrtype ; use coordinates ; use optimdata ; use primitive ; use math
implicit none

! In this subroutine, the primitive internal coordinate set is initialised and is used throughout optimisation.
! This is to avoid discontinuity in the DLC, as the primitive internal coordinates are defined once at the start.
! With each evaluation of DLC, the actual values of the primitives are re-calculated, but their definition remains the same.
! If DLC are not used, then primitive internal coordinates are not defined as they are not required.

integer(i4b) :: i, j, k
	
if (coordtype .eq. 1) then
	! If the primitive internal coordinates have already been defined, then we can skip this subroutine.
	if (ALLOCATED(prim_list)) then
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
	do i=1, nl
		j = (3 * nq) + i
		xopt(j:j+2) = xl(i:i+2)
	end do

	! Now, the primitive internal coordinates can be defined.
	if (primtype .eq. 0) then ! Total connectivity
		call define_prims_TC(nopt, opt, prim_list)
	else if (primtype .eq. 1) then ! Full definition
		call define_prims_full(nopt, opt, xopt, prim_list)	
	end if
end if

end subroutine initialize_prims