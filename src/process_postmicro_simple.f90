PROGRAM process_postmicro_simple
use nrtype
implicit none

! This program is run at each cycle of the QM/MM calculation, after the Tinker microiterations.
! It reads: * the Tinker geometry created by the microiterations.
!           * the new geometry from the Hessian optimization region.
! Then writes: * the collated geometry to "geom.xyz"
!              * the array of charges for inclusion in the next QM job.
!              * the array of points for calculating the electrostatic potential.

call alloc_coord()
call alloc_optim()
call read_expl_coord_simple()
call short_read_checkfile_simple()
call read_hessgeom_simple()
call collategeoms_simple()
call write_coord_simple()

END PROGRAM process_postmicro_simple
