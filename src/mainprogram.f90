PROGRAM qommmafort6
use nrtype ; use coordinates ; use optimdata
implicit none

! This program is run at each cycle of the QM/MM calculation.
! It reads the Tinker and QM/MM scratch files (geom.xyz and ProgFile)
! It reads the QM gradient file.
! It reads the QM resp output and converts this into the relevant gradients.
! It reads the QM Mulliken charges, and computes the difference between the
!    correct QM/MM force on the MM atoms, and that derived from Mulliken charges.
! It reads the Tinker gradient file.
! It writes out overall gradients.
! It writes the correction to the charge-charge QM/MM gradient.
! It writes a report.
! It can carry out constrained optimization using a penalty function.
! Up to three constraints may be used.
! The penalty function can be either harmonic or based on a tanh gradient
! This is version 8.

call alloc_coord()
call read_converg()
call alloc_optim()
call read_expl_coord()
call read_checkfile()
call read_tinkergrad()
call read_ab_initio()
call read_qmlatgrad()
call read_mulliken()
call evaluate_overall_grad()
call write_gradcorrection()
call evaluate_spring_force()
call update_opt_geometry()
call write_report()
call write_checkfile()
call update_full_geometry()
call det_sm_coord()
call write_sm_coord()
call write_hessgeom()
call write_mulliken()

END PROGRAM qommmafort6
