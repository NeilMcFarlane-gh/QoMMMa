SUBROUTINE evaluate_constraint_dlc
use nrtype ; use coordinates ; use optimdata ; use primitive ; use dlc_constraint
implicit none
integer(i4b) :: i, j, k

! This subroutine generates the constraint matrix used for dlc constraints, cdat.

! To calculate the primitive internal coordinates, the QM region coordinates need to be obtained first.
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

! Now calculate the primitive internal coordinates.
call calc_prims((nq + nl), nprim, prims, opt, xopt, prim_list)

! Now, generate the constrained coordinate array.
call gen_cons(ncon_prim, nprim, cdat, prim_list, prims)

! To avoid problems down the line, deallocate prims.
deallocate(prims)

END SUBROUTINE evaluate_constraint_dlc