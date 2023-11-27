PROGRAM reparam_string
use nrtype ; use coordinates ; use optimdata
implicit none

! This program is regularly called in both the double- and single-ended growing string method.
! In both growth and optimisation for DE, and only optimisation for SE.
! For a given node, this program generates new cartesian coordinates from a change in primitive internal coordinates.
! Details of what string reparameterisation actually is are found in the Python source code files.
! All this program does is read in the geometry, checkfile, and desired change in primitive internal coordinates.
! The output of this program is just new geometry files (cartesian and primitive) for the subsequent stages of string evolution.

call alloc_coord()
call alloc_optim()
call read_expl_coord()
call read_initfile()
call initialize_step()
call update_step_geometry()

END PROGRAM reparam_string
