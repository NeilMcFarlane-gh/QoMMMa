MODULE delocalised
use nrtype ; use coordinates ; use optimdata ; use math ; use primitive
implicit none

contains


	subroutine gen_Gmat(atom_num, n_prims, Bmat_p, Gmat)
	! Here, the G matrix used in the generatation of the DLC subspace is generated.
	!
	! ARGUMENTS:    atom_num    : Integer which represents the total number of atoms to be delocalised.
	!               n_prims     : Integer which represents the total number of primitive internal coordinates.
    !               Bmat_p      : 2D array which contains the primitive Wilson B matrix.	
    !               Gmat        : 2D array which contains the G matrix used to generate the DLC subspace.	
	
	implicit none
	integer(i4b) :: atom_num, n_prims
	real(sp) :: Bmat_p((3 * atom_num), n_prims)
	real(sp), allocatable :: Gmat(:,:)
	
	! The G matrix is simply the Wilson B matrix multiplied by its transpose.
	! By definition, it has dimensions of n_prims x n_prims, so it is allocated accordingly.
	if (.not. ALLOCATED(Gmat)) allocate(Gmat(n_prims, n_prims))
	Gmat(:,:) = 0.0
	Gmat = MATMUL(TRANSPOSE(Bmat_p), Bmat_p)

	end subroutine gen_Gmat
	
	
	subroutine diag_Gmat(atom_num, n_prims, Gmat, Umat, Rmat)
	! Here, the G matrix used in the generatation of the final DLC is diagaonalised, and the resulting non-redundant (U matrix) and redundant (R matrix - not used) subspaces are separated.
	!
	! ARGUMENTS:    atom_num    : Integer which represents the total number of atoms to be delocalised.
	!               n_prims     : Integer which represents the total number of primitive internal coordinates.
    !               Gmat        : 2D array which contains the G matrix used to generated the DLC subspace.
	!               Umat        : 2D array containing the eigenvectors with eigenvalues greater than zero. 
	!                             This is the DLC transformation vector set.
	!               Rmat        : 2D array containing the eigenvectors with eigenvalues equal to zero (with numerical precision considerations).
	!                             This can be used to transform to a redundant internal coordinate set, but this is not presently used.
	
	implicit none
	integer(i4b) :: atom_num, n_prims, j, i, U_vector_counter, R_vector_counter
	real(sp) :: Gmat(n_prims, n_prims), eigenvals(n_prims), eigenvecs(n_prims, n_prims), temp_eval
	real(sp), allocatable :: Umat(:,:), Rmat(:,:)
	
	! By definition, the U matrix should contain 3N-6 (where N is the total number of atoms) eigenvectors.
	! The R matrix, should contain the remaining eigenvectors, which must be at least 6.
	! These matrices are allocated appropriately.
	if (.not. ALLOCATED(Umat)) allocate(Umat(n_prims, ((3 * atom_num) - 6)))
	if (.not. ALLOCATED(Rmat)) allocate(Rmat(n_prims, (n_prims - ((3 * atom_num) - 6))))
	Umat(:,:) = 0.0
	Rmat(:,:) = 0.0

	! Extracting eigenvalues and eigenvectors...
	eigenvals = EVALS(Gmat, SIZE(Gmat,1))
	eigenvecs = EVECS(Gmat, SIZE(Gmat,1))

	U_vector_counter = 1
	R_vector_counter = 1
	do j=1, SIZE(eigenvals)
	    temp_eval = eigenvals(j)
		if (ABS(temp_eval) .lt. 1E-2) then
		    Rmat(:,R_vector_counter) = eigenvecs(:,j)
			R_vector_counter = R_vector_counter + 1
		else
		    Umat(:,U_vector_counter) = eigenvecs(:,j)
			U_vector_counter = U_vector_counter + 1
		end if
	end do
	
	end subroutine diag_Gmat
	
	
	subroutine gen_Bmat_DLC(atom_num, n_prims, Bmat_p, Umat, Bmat_dlc)
	! Here, the Wilson B matrix in primitive internal coordinate subspace is updated to DLC subspace.
	!
	! ARGUMENTS:    atom_num   : Integer which represents the total number of atoms to be delocalised.
	!               n_prims    : Integer which represents the total number of primitive internal coordinates.
	!               Bmat_p     : 2D array which contains the primitive Wilson B matrix.
	!				Umat       : 2D array containing the eigenvectors with eigenvalues greater than zero. 
	!                            This is the DLC transformation vector set. 	
	!               Bmat_dlc   : 2D array which contains the DLC Wilson B matrix.
	
	implicit none
	integer(i4b) :: atom_num, n_prims
	real(sp) :: Bmat_p((3 * atom_num), n_prims), Umat(n_prims, ((3 * atom_num) - 6))
	real(sp), allocatable :: Bmat_dlc(:,:)
	
	! The calculation of the B matrix in DLC subspace is a straightforward multiplication.
	if (.not. ALLOCATED(Bmat_dlc)) allocate(Bmat_dlc(((3 * atom_num) - 6), (3 * atom_num)))
	Bmat_dlc(:,:) = 0.0
	Bmat_dlc = MATMUL(TRANSPOSE(Umat), TRANSPOSE(Bmat_p))
	
	end subroutine gen_Bmat_DLC
	
	
	subroutine gen_DLC(atom_num, n_prims, Umat, prims, dlc)
	! Here, the DLC are actually generated from linear combinations of primitive internal coordinates and the U matrix.
	!
	! ARGUMENTS:  	atom_num : Integer which represents the total number of atoms to be delocalised.
	!               n_prims  : Integer which represents the total number of primitive internal coordinates.		
	!               Umat     : 2D array containing the eigenvectors with eigenvalues greater than zero. 
	!                          This is the DLC transformation vector set. 
	!               prims    : 1D array containing all the primitive internal coordinates associated with this input.				
	!				dlc      : 1D array containing the delocalised internal coordinate set.
	
	implicit none
	integer(i4b) :: n_prims, atom_num, i, j
	real(sp) :: Umat(n_prims, ((3 * atom_num) - 6)), prims(n_prims)
	real(sp), allocatable :: dlc(:)
	
	! The number of DLC is equal to the number of eigenvectors in the U matrix.
	if (.not. ALLOCATED(dlc)) allocate(dlc(SIZE(Umat,2)))
	dlc(:) = 0.0
	do i=1, SIZE(Umat,2)
		do j=1, SIZE(prims,1)
			dlc(i) = dlc(i) + (Umat(j,i) * prims(j))
		end do
	end do
	
	end subroutine gen_DLC
	
	
	subroutine refresh_DLC(atom_num, coords)
	! Here, the DLC are refreshed or generated for the first time.
	! All arrays used in the generation are deallocated so that they can be allocated appropriately again.
	!
	! ARGUMENTS:    atom_num    : Integer which represents the total number of atoms to be delocalised.
	!               coords      : 1D array containing all the cartesian coordinates of the system.
	
	implicit none
	integer(i4b) :: atom_num
	real(sp) :: coords(atom_num * 3)
	
	! First, (if necessary) deallocate all arrays used in DLC generation.
	if (ALLOCATED(prims) .or. ALLOCATED(Bmat_p) &
	& .or. ALLOCATED(Bmat_dlc) .or. ALLOCATED(Gmat) .or. ALLOCATED(Umat) &
	& .or. ALLOCATED(Rmat) .or. ALLOCATED(to_generate)) then
			deallocate(prim_list, prims, Bmat_p, Bmat_dlc, Gmat, Umat, Rmat, to_generate)
	end if

	! Now, all the DLC subroutines can be called to reallocate the arrays.
	call gen_prims(atom_num, to_generate, coords, prims, prim_list)
	call gen_Bmat_prims(atom_num, nprim, coords, prim_list, Bmat_p)
	call gen_Gmat(atom_num, nprim, Bmat_p, Gmat)
	call diag_Gmat(atom_num, nprim, Gmat, Umat, Rmat)
	call gen_DLC(atom_num, nprim, Umat, prims, dlc)
	call gen_Bmat_DLC(atom_num, nprim, Bmat_p, Umat, Bmat_dlc)

	end subroutine refresh_DLC
	
	subroutine gen_grad_cart_to_DLC(atom_num, n_prims, Bmat_dlc, g, g_dlc)
	! Here, the gradient array in cartesian subspace is updated to DLC subspace.
	!
	! ARGUMENTS:    atom_num : Integer which represents the total number of atoms to be delocalised.
	!               n_prims  : Integer which represents the total number of primitive internal coordinates.
	!               Bmat_dlc : 2D array which contains the DLC Wilson B matrix.
	!				g        : 1D array containing the gradient in cartesian subspace.
	!               g_dlc    : 1D array containing the gradient in DLC subspace.
	
	implicit none
	integer(i4b) :: atom_num, n_prims
	real(sp) :: Bmat_dlc((3 * atom_num), ((3 * atom_num) - 6)), g(atom_num * 3)
	real(sp), allocatable :: g_dlc(:)
	real(sp) :: BT_Ginv(((3 * atom_num) - 6), (3 * atom_num))

	if (.not. ALLOCATED(g_dlc)) allocate(g_dlc((3 * atom_num) - 6))

	! Using single value decomposition, the Moore-Penrose inverse is constructed and used to convert the gradient array.
	BT_Ginv = MATMUL(TRANSPOSE(Bmat_dlc), SVD_INVERSE(MATMUL(Bmat_dlc, TRANSPOSE(Bmat_dlc)), SIZE(Bmat_dlc, 1), SIZE(Bmat_dlc, 1)))
	g_dlc = MATMUL(g, TRANSPOSE(BT_Ginv))
	
	end subroutine gen_grad_cart_to_DLC
	
	
	subroutine gen_hess_prim_to_DLC(atom_num, n_prims, Umat, hess, hess_dlc)
	! Here, the cartesian hessian matrix is updated to DLC subspace.
	!
	! ARGUMENTS:    atom_num : Integer which represents the total number of atoms to be delocalised.
	!               n_prims  : Integer which represents the total number of primitive internal coordinates.
	!				Umat     : 2D array containing the eigenvectors with eigenvalues greater than zero. 
	!                          This is the DLC transformation vector set. 	
	!               hess     : 2D array which contains the hessian matrix in primitive subspace.
	!               hess_dlc : 2D array which contains the hessian matrix in DLC subspace.
	
	implicit none
	integer(i4b) :: atom_num, n_prims
	real(sp) :: Umat(n_prims, ((3 * atom_num) - 6)), hess(n_prims, n_prims)
	real(sp), allocatable :: hess_dlc(:,:)

	! First, the hessian matrix is DLC subspace should be allocated.
	if (.not. ALLOCATED(hess_dlc)) allocate(hess_dlc(((3 * atom_num) - 6), ((3 * atom_num) - 6)))

	! Now, the hessian can be calculated in DLC subspace by a simple multiplication procedure.
	hess_dlc = MATMUL(TRANSPOSE(Umat), MATMUL(hess, Umat))

	end subroutine gen_hess_prim_to_DLC
	
	
	subroutine DLC_to_cart(atom_num, n_prims, dS, dlc, x_1, x_2, Bmat_dlc)
	! Here, the DLC are converted to cartesian coordinates via an iterative procedure.
	! This is used at the end of an optimisation cycle to restore the cartesian coordinates for the next gradient evaluation.
	!
	! ARGUMENTS:    atom_num : Integer which represents the total number of atoms to be delocalised.
	!               n_prims  : Integer which represents the total number of primitive internal coordinates.	
	!           	dS       : 1D array containing the change in DLC.
	!				S_1      : 1D array containing the delocalised internal coordinate set of the starting point.  
    !               x_1      : 2D array containing all the cartesian coordinates of the starting point.	
	!				x_2      : 2D array containing all the cartesian coordinates after conversion.
	!				Bmat_dlc : 2D array containing the Wilson B matrix used to convert between cartesian and DLC.
	
	implicit none
	integer(i4b) :: atom_num, n_prims, i, iter_counter
	real(sp) :: dS((3 * atom_num) - 6), dlc((3 * atom_num) - 6), old_dlc((3 * atom_num) - 6), x_1(3 * atom_num), x_2(3 * atom_num)
	real(sp), allocatable :: S_n(:), Bmat_dlsc_n(:,:)
	real(sp) :: BT_Ginv(((3 * atom_num) - 6), (3 * atom_num)), Bmat_dlc((3 * atom_num), ((3 * atom_num) - 6))
	real(sp) :: dS_actual((3 * atom_num) - 6), check((3 * atom_num) - 6), temp_check, step_scale, temp_dS
	real(sp) :: init_dS((3 * atom_num) - 6), target_dlc((3 * atom_num) - 6), dx(3 * atom_num)
	real(sp) :: xyz_rms_1, xyz_rms_2
	logical :: convergence
	
    ! Since cartesians are rectilinear and internal coordinates are curvilinear, a simple transformation cannot be used.
    ! Instead, an iterative transformation procedure must be used.
    ! The expression B(transpose) * G(inverse) is initialised as it is used to convert between coordinate systems.
	BT_Ginv = MATMUL(TRANSPOSE(Bmat_dlc), SVD_INVERSE(MATMUL(Bmat_dlc, TRANSPOSE(Bmat_dlc)), SIZE(Bmat_dlc, 1), SIZE(Bmat_dlc, 1)))

	! Stashing some values for convergence criteria.
	step_scale = 1.0
	convergence = .FALSE.
	xyz_rms_1 = 0
	xyz_rms_2 = 0
	init_dS(:) = dS(:)
	target_dlc(:) = dlc(:) + init_dS(:)

100	do while (convergence .eqv. .FALSE.) 
		! First, the DLC from the previous iteration is saved.
		old_dlc(:) = dlc(:)
	
		! The change in cartesian coordinates associated with the change in DLC is calculated.
		dx(:) = 0.0
		dx = MATMUL(TRANSPOSE(BT_Ginv), dS)
		
		! The root-mean-square change is used as a convergence criteria, so it is evaluated.
		xyz_rms_2 = RMSD_CALC(dx, x_1, SIZE(x_1))
		
		! The new cartesian geometry is evaluated, and the new DLC geometry is obtained.
		x_2 = x_1 + dx
		dlc(:) = 0.0
		call refresh_DLC(atom_num, x_2)
		
		! The Moore-Penrose inverse is constructed for the next iteration.
		!BT_Ginv(:,:) = 0.0
		!BT_Ginv = MATMUL(TRANSPOSE(Bmat_dlc), SVD_INVERSE(MATMUL(Bmat_dlc, TRANSPOSE(Bmat_dlc)), SIZE(Bmat_dlc, 1), SIZE(Bmat_dlc, 1)))

		! The change in DLC for the next iteration is evaluated.
		dS(:) = 0.0
		dS = (target_dlc - dlc) * step_scale
		do i=1, SIZE(dS)
			if (ABS(dS(i)) > ABS(init_dS(i))) then
				step_scale = step_scale * 0.75
				go to 100
			end if
		end do
		
		! Now, the three exit conditions should be checked...
		! The first ending condition for this transformation is when the root-mean-square change in cartesians is less than 10^-6.
        ! The second ending condition for this transformation is when the difference in root-mean-square change in cartesians between iteration i and i+1 is less than 10^-12.
        ! The third ending condition for this transformation is when the difference between the target DIC and the calculated DLC is less than 10^-6.
		check = target_dlc - dlc
		if (ABS(xyz_rms_2) < 1E-06) then
			convergence = .TRUE.
		else if (ABS(xyz_rms_2 - xyz_rms_1) < 1E-12) then
			convergence = .TRUE.
		end if
		do i=1, SIZE(check)
			temp_check = check(i)
			if (ABS(temp_check) < 1E-12) then
				convergence = .TRUE.
			end if
		end do
		
		! Values which are calculated from the iterative procedure are updated to be ith property for the next iteration.
		x_1(:) = x_2(:)
		xyz_rms_1 = xyz_rms_2

		! In some cases, cartesians cannot be solved, so an exit condition must exist for this case.
		iter_counter = iter_counter + 1
		if (iter_counter == 1000) then
			print *, "COULDN'T SOLVE DLC"
			exit
		end if
	end do

	end subroutine DLC_to_cart
	
	
END MODULE delocalised