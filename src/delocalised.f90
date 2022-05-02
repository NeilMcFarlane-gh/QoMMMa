MODULE delocalised
use nrtype ; use coordinates ; use optimdata ; use math ; use primitive
implicit none

contains


	subroutine gen_Gmat(atom_num, n_dlc, n_prims, Bmat_p, Gmat)
	! Here, the G matrix used in the generatation of the DLC subspace is generated.
	!
	! ARGUMENTS:    atom_num    : Integer which represents the total number of atoms to be delocalised.
	!               n_dlc       : Integer which represents the number of delocalised internal coordinates (by definition, 3N-6).
	!               n_prims     : Integer which represents the total number of primitive internal coordinates.
    !               Bmat_p      : 2D array which contains the primitive Wilson B matrix.	
    !               Gmat        : 2D array which contains the G matrix used to generate the DLC subspace.	
	
	implicit none
	integer(i4b) :: atom_num, n_prims, n_dlc
	real(sp) :: Bmat_p((3 * atom_num), n_prims), deter
	real(sp), allocatable :: Gmat(:,:)
	
	! The G matrix is simply the Wilson B matrix multiplied by its transpose.
	! By definition, it has dimensions of n_prims x n_prims, so it is allocated accordingly.
	if (.not. ALLOCATED(Gmat)) allocate(Gmat(n_prims, n_prims))
	Gmat(:,:) = 0.0
	Gmat = MATMUL(TRANSPOSE(Bmat_p), Bmat_p)
	
	! By definition, the G matrix must have a determinant of zero. 
	! This criteria is checked to mitigate errors down the line.
	deter = DETERMINANT(Gmat, n_prims)
	if (deter .gt. 1E-6) then
		write (*,*) "Error; the determinant of the G matrix is not zero - there must be some problem."
        stop
	end if

	end subroutine gen_Gmat
	
	
	subroutine diag_Gmat(atom_num, n_dlc, n_prims, Gmat, Umat, Rmat)
	! Here, the G matrix used in the generatation of the final DLC is diagaonalised, and the resulting non-redundant (U matrix) and redundant (R matrix - not used) subspaces are separated.
	!
	! ARGUMENTS:    atom_num    : Integer which represents the total number of atoms to be delocalised.
	!               n_dlc       : Integer which represents the number of delocalised internal coordinates (by definition, 3N-6).
	!               n_prims     : Integer which represents the total number of primitive internal coordinates.
    !               Gmat        : 2D array which contains the G matrix used to generated the DLC subspace.
	!               Umat        : 2D array containing the eigenvectors with eigenvalues greater than zero. 
	!                             This is the DLC transformation vector set.
	!               Rmat        : 2D array containing the eigenvectors with eigenvalues equal to zero (with numerical precision considerations).
	!                             This can be used to transform to a redundant internal coordinate set, but this is not presently used.
	
	implicit none
	integer(i4b) :: atom_num, n_prims, n_dlc, j, U_vector_counter, R_vector_counter
	real(sp) :: Gmat(n_prims, n_prims), eigenvals(n_prims), eigenvecs(n_prims, n_prims), temp_eval, temp
	real(sp), allocatable :: Umat(:,:), Rmat(:,:)
	
	! By definition, the U matrix should contain 3N-6 (where N is the total number of atoms) eigenvectors.
	! The R matrix, should contain the remaining eigenvectors, which must be at least 6.
	! These matrices are allocated appropriately.
	if (.not. ALLOCATED(Umat)) allocate(Umat(n_dlc, n_prims))
	if (.not. ALLOCATED(Rmat)) allocate(Rmat((n_prims - n_dlc), n_prims))
	Umat(:,:) = 0.0
	Rmat(:,:) = 0.0

	! Extracting eigenvalues and eigenvectors...
	eigenvals = EVALS(Gmat, n_prims)
	eigenvecs = EVECS(Gmat, n_prims)
	U_vector_counter = 1
	R_vector_counter = 1
	do j=1, n_prims
	    temp_eval = eigenvals(j)
		if (ABS(temp_eval) .lt. 1E-10) then
		    Rmat(R_vector_counter,:) = eigenvecs(:,j)
			R_vector_counter = R_vector_counter + 1
		else
		    Umat(U_vector_counter,:) = eigenvecs(:,j)
			U_vector_counter = U_vector_counter + 1
		end if
	end do

	end subroutine diag_Gmat
	
	
	subroutine gen_Bmat_DLC(atom_num, n_dlc, n_prims, Bmat_p, Umat, Bmat_dlc)
	! Here, the Wilson B matrix in primitive internal coordinate subspace is updated to DLC subspace.
	!
	! ARGUMENTS:    atom_num   : Integer which represents the total number of atoms to be delocalised.
	!               n_dlc       : Integer which represents the number of delocalised internal coordinates (by definition, 3N-6).
	!               n_prims    : Integer which represents the total number of primitive internal coordinates.
	!               Bmat_p     : 2D array which contains the primitive Wilson B matrix.
	!				Umat       : 2D array containing the eigenvectors with eigenvalues greater than zero. 
	!                            This is the DLC transformation vector set. 	
	!               Bmat_dlc   : 2D array which contains the DLC Wilson B matrix.
	
	implicit none
	integer(i4b) :: atom_num, n_prims, n_dlc
	real(sp) :: Bmat_p((3 * atom_num), n_prims), Umat(n_prims, n_dlc)
	real(sp), allocatable :: Bmat_dlc(:,:)

	! The Wilson B matrix in DLC subspace is allocated. By definition, its dimensions are (3N-6) x (3N), where N is the number of atoms.
	if (.not. ALLOCATED(Bmat_dlc)) allocate(Bmat_dlc((3 * atom_num), n_dlc))
	Bmat_dlc(:,:) = 0.0
	
	! The calculation of the B matrix in DLC subspace is a straightforward multiplication.	
	Bmat_dlc = MATMUL(Bmat_p, Umat)

	end subroutine gen_Bmat_DLC
	
	
	subroutine gen_DLC(atom_num, n_dlc, n_prims, Umat, prims, dlc)
	! Here, the DLC are actually generated from linear combinations of primitive internal coordinates and the U matrix.
	!
	! ARGUMENTS:  	atom_num : Integer which represents the total number of atoms to be delocalised.
	!               n_dlc    : Integer which represents the number of delocalised internal coordinates (by definition, 3N-6).
	!               n_prims  : Integer which represents the total number of primitive internal coordinates.		
	!               Umat     : 2D array containing the eigenvectors with eigenvalues greater than zero. 
	!                          This is the DLC transformation vector set. 
	!               prims    : 1D array containing all the primitive internal coordinates associated with this input.				
	!				dlc      : 1D array containing the delocalised internal coordinate set.
	
	implicit none
	integer(i4b) :: n_prims, atom_num, n_dlc, i, j
	real(sp) :: Umat(n_prims, n_dlc), prims(n_prims)
	real(sp), allocatable :: dlc(:)
	
	! The DLC matrix is allocated.
	if (.not. ALLOCATED(dlc)) allocate(dlc(n_dlc))
	dlc(:) = 0.0

	! The DLC are solved simply by linear combinations of the U matrix with the primitive interal coordinates.
	do i=1, n_dlc
		do j=1, n_prims
			dlc(i) = dlc(i) + (Umat(j,i) * prims(j))
		end do
	end do

	end subroutine gen_DLC
	
	
	subroutine refresh_DLC(atom_num, n_dlc, coords)
	! Here, the DLC are refreshed or generated for the first time.
	! All arrays used in the generation are deallocated so that they can be allocated appropriately again.
	!
	! ARGUMENTS:    atom_num    : Integer which represents the total number of atoms to be delocalised.
	!               n_dlc       : Integer which represents the number of delocalised internal coordinates (by definition, 3N-6).
	!               coords      : 1D array containing all the cartesian coordinates of the system.
	
	implicit none
	integer(i4b) :: atom_num, n_dlc
	real(sp) :: coords(atom_num * 3)
	
	! First, (if necessary) deallocate all arrays used in DLC generation.
	if (ALLOCATED(prims) .or. ALLOCATED(Bmat_p) &
	& .or. ALLOCATED(Bmat_dlc) .or. ALLOCATED(Gmat) .or. ALLOCATED(Umat) &
	& .or. ALLOCATED(Rmat)) then
			deallocate(prims, Bmat_p, Bmat_dlc, Gmat, Umat, Rmat)
	end if

	! Now, all the DLC subroutines can be called to reallocate the arrays.
	call calc_prims(atom_num, nprim, prims, coords, prim_list)
	call gen_Bmat_prims(atom_num, nprim, coords, prim_list, Bmat_p)
	call gen_Gmat(atom_num, n_dlc, nprim, Bmat_p, Gmat)
	call diag_Gmat(atom_num, n_dlc, nprim, Gmat, Umat, Rmat)
	call gen_DLC(atom_num, n_dlc, nprim, Umat, prims, dlc)
	call gen_Bmat_DLC(atom_num, n_dlc, nprim, Bmat_p, Umat, Bmat_dlc)

	end subroutine refresh_DLC
	
	
	subroutine gen_grad_cart_to_DLC(atom_num, n_dlc, n_prims, Bmat_dlc, g, g_dlc)
	! Here, the gradient array in cartesian subspace is updated to DLC subspace.
	!
	! ARGUMENTS:    atom_num : Integer which represents the total number of atoms to be delocalised.
	!               n_dlc    : Integer which represents the number of delocalised internal coordinates (by definition, 3N-6).
	!               n_prims  : Integer which represents the total number of primitive internal coordinates.
	!               Bmat_dlc : 2D array which contains the DLC Wilson B matrix.
	!				g        : 1D array containing the gradient in cartesian subspace.
	!               g_dlc    : 1D array containing the gradient in DLC subspace.
	
	implicit none
	integer(i4b) :: atom_num, n_prims, n_dlc
	real(sp) :: Bmat_dlc((3 * atom_num), n_dlc), g(atom_num * 3)
	real(sp), allocatable :: g_dlc(:)
	real(sp) :: BT_Ginv(n_dlc, (3 * atom_num)), Gmat(n_dlc, n_dlc)
	
	! The primitive gradient array is allocated. By definition, its dimensions are the number of delocalised interal coordinates.
	if (.not. ALLOCATED(g_dlc)) allocate(g_dlc(n_dlc))
	g_dlc(:) = 0.0
	
	! Using single value decomposition, the Moore-Penrose inverse is constructed and used to convert the gradient array.
	Gmat = MATMUL(TRANSPOSE(Bmat_dlc), Bmat_dlc)
	BT_Ginv = MATMUL(SVD_INVERSE(Gmat, n_dlc, n_dlc), TRANSPOSE(Bmat_dlc))
	g_dlc = MATMUL(g, TRANSPOSE(BT_Ginv))

	end subroutine gen_grad_cart_to_DLC
	
	
	subroutine gen_hess_prim_to_DLC(atom_num, n_dlc, n_prims, Umat, hess, hess_dlc)
	! Here, the cartesian hessian matrix is updated to DLC subspace.
	!
	! ARGUMENTS:    atom_num : Integer which represents the total number of atoms to be delocalised.
	!               n_dlc    : Integer which represents the number of delocalised internal coordinates (by definition, 3N-6).
	!               n_prims  : Integer which represents the total number of primitive internal coordinates.
	!				Umat     : 2D array containing the eigenvectors with eigenvalues greater than zero. 
	!                          This is the DLC transformation vector set. 	
	!               hess     : 2D array which contains the hessian matrix in primitive subspace.
	!               hess_dlc : 2D array which contains the hessian matrix in DLC subspace.
	
	implicit none
	integer(i4b) :: atom_num, n_prims, n_dlc
	real(sp) :: Umat(n_prims, n_dlc), hess(n_prims, n_prims)
	real(sp), allocatable :: hess_dlc(:,:)

	! First, the hessian matrix is DLC subspace should be allocated.
	if (.not. ALLOCATED(hess_dlc)) allocate(hess_dlc(n_dlc, n_dlc))
	hess_dlc(:,:) = 0.0

	! Now, the hessian can be calculated in DLC subspace by a simple multiplication procedure.
	hess_dlc = MATMUL(TRANSPOSE(Umat), MATMUL(hess, Umat))

	end subroutine gen_hess_prim_to_DLC
	
	
	subroutine DLC_to_cart(atom_num, n_dlc, n_prims, dS, dlc, x_1, x_2, Bmat_dlc)
	! Here, the DLC are converted to cartesian coordinates via an iterative procedure.
	! This is used at the end of an optimisation cycle to restore the cartesian coordinates for the next gradient evaluation.
	! A note on notation: changes in dlc are often called 'dS' here as 'S' is another notation widely used for DLC.
	!
	! ARGUMENTS:    atom_num : Integer which represents the total number of atoms to be delocalised.
	!               n_dlc    : Integer which represents the number of delocalised internal coordinates (by definition, 3N-6).
	!               n_prims  : Integer which represents the total number of primitive internal coordinates.
	!				dS       : 1D array containing the change in DLC for which cartesians are to be solved for.
	!				dlc      : 1D array containing the delocalised internal coordinate set.
	!               x_1      : 1D array containing all the cartesian coordinates of the system prior to the change in DLC.
	!               x_1      : 1D array containing all the cartesian coordinates of the system following the change in DLC.
	!               Bmat_dlc : 2D array which contains the DLC Wilson B matrix.
	
	implicit none
	integer(i4b) :: k, atom_num, n_prims, n_dlc
	integer(i4b), parameter :: resolution = 500
	real(sp) :: dS_norm_save, dS_norm, dS_norm_init
	real(sp) :: x_1(3 * atom_num), x_2(3 * atom_num), dx(3 * atom_num)
	real(sp) :: dx_step(3 * atom_num), dx_temp(3 * atom_num), temp_x(3 * atom_num), dx_save(3 * atom_num)
	real(sp) :: BT_Ginv(n_dlc, (3 * atom_num)), Bmat_dlc((3 * atom_num), n_dlc), Gmat(n_dlc, n_dlc)
	real(sp) :: init_dlc(n_dlc), init_dS(n_dlc), target_dlc(n_dlc), dS(n_dlc)
	real(sp) :: dS_temp(n_dlc), dlc(n_dlc), dS_save(n_dlc)

    ! The expression B(transpose) * G(inverse) is initialised as it is used to convert between coordinate systems.
	Gmat = MATMUL(TRANSPOSE(Bmat_dlc), Bmat_dlc)
	BT_Ginv = MATMUL(SVD_INVERSE(Gmat, n_dlc, n_dlc), TRANSPOSE(Bmat_dlc))
	dx = MATMUL(TRANSPOSE(BT_Ginv), dS)

	! Stashing the initial and target values.
	init_dS(:) = dS(:)
	init_dlc(:) = dlc(:)
	target_dlc(:) = dlc(:) + init_dS(:)

	! In the original implementation by Baker, the transformation procedure is iterative using the equation above.
	! However, this procedure is rather unstable when it comes to large changes in DLC, or near 180 dihedral angles.
	! In the research group of P. Ayers, they use a procedure of calculating a number cartesian coordinates and comparing to DLC.
	! The cartesian set closest to the target DLC is then taken as the new coordinates.
	! In this implementation, the algorithm is inspired by this method where a 'resolution' is defined.
	! This resolution defines how many small steps in cartesians are made to find the DLC: higher resolution implies more accurate DLC.
	! While this method is not as accurate or fast as the original implementation, it is much more stable.
	dx_step = dx / resolution
	do k=1, resolution
		if (k == 1) then
			dx_temp = dx_step
			temp_x = dx_temp + x_1
			call refresh_DLC(atom_num, n_dlc, temp_x)
			dS_temp = target_dlc - dlc
			dS_norm_save = NORM2(dS_temp)
			dx_save = dx_temp
			dS_save = dS_temp
		else
			dx_temp = dx_temp + dx_step
			temp_x = dx_temp + x_1
			call refresh_DLC(atom_num, n_dlc, temp_x)
			dS_temp = target_dlc - dlc
			dS_norm = NORM2(dS_temp)
			if (dS_norm .lt. dS_norm_save) then
				dx_save = dx_temp
				dS_save = dS_temp
				dS_norm_save = dS_norm
			end if
		end if
	end do
		
	! The new DLC set (which should be close to init_dlc) and the new cartesian coordinates are saved.
	x_2 = x_1 + dx_save
	
	end subroutine DLC_to_cart
	
	
END MODULE delocalised