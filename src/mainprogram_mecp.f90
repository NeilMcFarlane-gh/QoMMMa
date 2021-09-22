PROGRAM qommmafort6
use nrtype
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
call alloc_optim_mecp()
call read_expl_coord()
call read_checkfile_mecp()
call read_tinkergrad()
call read_ab_initio_mecp()
call evaluate_overall_grad_mecp()
call Effective_Gradient_mecp()
call update_opt_geometry_mecp()
call read_qmlatgrad_mecp()
call read_mulliken_mecp()
call write_gradcorrection()
call evaluate_spring_force()
call write_report_mecp()
call write_checkfile_mecp()
call update_full_geometry()
call det_sm_coord()
call write_sm_coord()
call write_hessgeom()
call write_mulliken()

END PROGRAM qommmafort6
