PROGRAM setup_first
use nrtype
implicit none

! This program is run at the beginning of the QM/MM calculation.
! It reads the Tinker and QM/MM input file (geom.xyz and InitFile)
! It initializes output, and helps prepare the first input files.

call alloc_coord()
call alloc_optim()
call read_expl_coord()
call read_charges()
call read_initfile()
call initialize_update_geom()
call det_sm_coord()
call write_sm_coord()
call write_atomic()
call write_charges()
call write_resp()
call write_inactive()
call write_first_report()
call write_first_checkfile()

END PROGRAM setup_first
