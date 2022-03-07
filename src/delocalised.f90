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
	if (.not. ALLOCATED(Umat)) allocate(Umat(((3 * atom_num) - 6), n_prims))
	if (.not. ALLOCATED(Rmat)) allocate(Rmat((n_prims - ((3 * atom_num) - 6)), n_prims))
	Umat(:,:) = 0.0
	Rmat(:,:) = 0.0

	! Extracting eigenvalues and eigenvectors...
	eigenvals = EVALS(Gmat, SIZE(Gmat,1))
	eigenvecs = EVECS(Gmat, SIZE(Gmat,1))
	U_vector_counter = 1
	R_vector_counter = 1
	do j=1, SIZE(eigenvals)
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
	
	! The Wilson B matrix in DLC subspace is allocated. By definition, its dimensions are (3N-6) x (3N), where N is the number of atoms.
	if (.not. ALLOCATED(Bmat_dlc)) allocate(Bmat_dlc(((3 * atom_num) - 6), (3 * atom_num)))
	Bmat_dlc(:,:) = 0.0
	
	! The calculation of the B matrix in DLC subspace is a straightforward multiplication.	
	Bmat_dlc = MATMUL(Bmat_p, Umat)
	
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
	
	! The DLC matrix is allocated. The number of DLC is equal to the number of eigenvectors in the U matrix.
	if (.not. ALLOCATED(dlc)) allocate(dlc(SIZE(Umat,2)))
	dlc(:) = 0.0

	! The DLC are solved simply by linear combinations of the U matrix with the primitive interal coordinates.
	do i=1, SIZE(Umat,2)
		do j=1, SIZE(prims,1)
			dlc(i) = dlc(i) + (Umat(j,i) * prims(j))
		end do
	end do

	end subroutine gen_DLC
	
	
	subroutine refresh_DLC(atom_num, coords, init)
	! Here, the DLC are refreshed or generated for the first time.
	! All arrays used in the generation are deallocated so that they can be allocated appropriately again.
	!
	! ARGUMENTS:    atom_num    : Integer which represents the total number of atoms to be delocalised.
	!               coords      : 1D array containing all the cartesian coordinates of the system.
	
	implicit none
	integer(i4b) :: atom_num
	real(sp) :: coords(atom_num * 3)
	logical :: init

	! First, (if necessary) deallocate all arrays used in DLC generation.
	if (ALLOCATED(prims) .or. ALLOCATED(Bmat_p) &
	& .or. ALLOCATED(Bmat_dlc) .or. ALLOCATED(Gmat) .or. ALLOCATED(Umat) &
	& .or. ALLOCATED(Rmat)) then
			deallocate(prims, Bmat_p, Bmat_dlc, Gmat, Umat, Rmat)
	end if

	! Now, all the DLC subroutines can be called to reallocate the arrays.
	call calc_prims(atom_num, nprim, prims, coords, prim_list)
	call gen_Bmat_prims(atom_num, nprim, coords, prim_list, Bmat_p)
	call gen_Gmat(atom_num, nprim, Bmat_p, Gmat)
	call diag_Gmat(atom_num, nprim, Gmat, Umat, Rmat)
	call gen_DLC(atom_num, nprim, Umat, prims, dlc)
	call gen_Bmat_DLC(atom_num, nprim, Bmat_p, Umat, Bmat_dlc)
	!if (init .eqv. .True.) then
		!print *, "calculated", dlc
	!end if

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

	! The primitive gradient array is allocated. By definition, its dimensions are the number of delocalised interal coordinates.
	if (.not. ALLOCATED(g_dlc)) allocate(g_dlc((3 * atom_num) - 6))
	g_dlc(:) = 0.0

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
	hess_dlc(:,:) = 0.0

	! Now, the hessian can be calculated in DLC subspace by a simple multiplication procedure.
	hess_dlc = MATMUL(TRANSPOSE(Umat), MATMUL(hess, Umat))

	end subroutine gen_hess_prim_to_DLC
	
	
	subroutine DLC_to_cart(atom_num, n_prims, dS, dlc, x_1, x_2, Bmat_dlc)
	! Here, the DLC are converted to cartesian coordinates via an iterative procedure.
	! This is used at the end of an optimisation cycle to restore the cartesian coordinates for the next gradient evaluation.
	!
	! ARGUMENTS:    atom_num : Integer which represents the total number of atoms to be delocalised.
	!               n_prims  : Integer which represents the total number of primitive internal coordinates.
	!				dS       : 1D array containing the change in DLC for which cartesians are to be solved for.
	!				dlc      : 1D array containing the delocalised internal coordinate set.
	!               x_1      : 1D array containing all the cartesian coordinates of the system prior to the change in DLC.
	!               x_1      : 1D array containing all the cartesian coordinates of the system following the change in DLC.
	!               Bmat_dlc : 2D array which contains the DLC Wilson B matrix.
	
	implicit none
	integer(i4b) :: k, atom_num, n_prims
	integer(i4b), parameter :: resolution = 10000
	real(sp) :: dS_norm_save, dS_norm
	real(sp) :: dS((3 * atom_num) - 6), dlc((3 * atom_num) - 6), x_1(3 * atom_num), x_2(3 * atom_num)
	real(sp) :: BT_Ginv(((3 * atom_num) - 6), (3 * atom_num)), Bmat_dlc((3 * atom_num), ((3 * atom_num) - 6))
	real(sp) :: init_dlc((3 * atom_num) - 6), init_dS((3 * atom_num) - 6), target_dlc((3 * atom_num) - 6)
	real(sp) :: dx(3 * atom_num), dx_step(3 * atom_num), dx_temp(3 * atom_num), dS_temp((3 * atom_num) - 6)
	real(sp) :: temp_x(3 * atom_num), dx_save(3 * atom_num), dS_save((3 * atom_num) - 6)
	logical :: init = .False.

    ! The expression B(transpose) * G(inverse) is initialised as it is used to convert between coordinate systems.
	BT_Ginv = MATMUL(TRANSPOSE(Bmat_dlc), SVD_INVERSE(MATMUL(Bmat_dlc, TRANSPOSE(Bmat_dlc)), SIZE(Bmat_dlc, 1), SIZE(Bmat_dlc, 1)))

	! Stashing the initial and target values.
	init_dS(:) = dS(:)
	init_dlc(:) = dlc(:)
	target_dlc(:) = dlc(:) + init_dS(:)

	! The change in cartesian coordinates associated with the change in DLC is calculated.
	! This is the initial drving force for the algorithm below.
	dx(:) = 0.0
	dx = MATMUL(TRANSPOSE(BT_Ginv), dS)
	
	! In the original implementation by Baker, the transformation procedure is iterative using the equation above.
	! However, this procedure is rather unstable when it comes to large changes in DLC, or near 180 dihedral angles.
	! In the research group of P. Ayers, they use a procedure of calculating a number cartesian coordinates and comparing to DLC.
	! The cartesian set closest to the target DLC is then taken as the new coordinates.
	! In this implementation, the algorithm is inspired by this method where a 'resolution' is defined.
	! This resolution defines how many small steps in cartesians are made to find the DLC - higher resolution = more accurate DLC.
	! While this method is not as accurate or fast as the original implementation, it is much more stable.
	dx_step = dx / resolution
	do k=1, resolution
		if (k == 1) then
			dx_temp = dx_step
			temp_x = dx_temp + x_1
			call refresh_DLC(atom_num, temp_x, init)
			dS_temp = target_dlc - dlc
			dS_norm_save = NORM2(dS_temp)
			dx_save = dx_temp
			dS_save = dS_temp
		else
			dx_temp = dx_temp + dx_step
			temp_x = dx_temp + x_1
			call refresh_DLC(atom_num, temp_x, init)
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
	print *, target_dlc
	call refresh_DLC(atom_num, x_2, init)
	print *, dlc
	
	end subroutine DLC_to_cart
	
	
END MODULE delocalised