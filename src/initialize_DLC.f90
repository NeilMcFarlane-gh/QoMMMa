SUBROUTINE initialize_DLC()
use nrtype ; use coordinates ; use optimdata ; use primitive ; use delocalised ; use math
implicit none

if (ncon_prim.gt.0) then
	call evaluate_constraint_dlc()
end if

call refresh_DLC(nopt, ndlc, xopt, cdat)

end subroutine initialize_DLC
