SUBROUTINE initialize_DLC()
use nrtype ; use coordinates ; use optimdata ; use primitive ; use delocalised ; use math
implicit none

call refresh_DLC(nopt, ndlc, xopt, cdat)
call write_U()

end subroutine initialize_DLC