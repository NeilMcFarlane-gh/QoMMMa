SUBROUTINE initialize_DLC()
use nrtype ; use coordinates ; use optimdata ; use primitive ; use delocalised ; use math
implicit none

integer(i4b) :: i, j, k

if (coordtype .eq. 1) then
	if (ncon_prim.gt.0) then
		call evaluate_constraint_dlc()
	end if
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
	call refresh_DLC(nopt, ndlc, xopt, cdat)
	
end if

end subroutine initialize_DLC
