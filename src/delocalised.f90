MODULE delocalised
use nrtype ; use coordinates ; use optimdata ; use math ; use primitive ; use dlc_constraint
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
	real(dp) :: Bmat_p(n_prims, (3 * atom_num)), deter
	real(dp), allocatable :: Gmat(:,:)
	
	! The G matrix is simply the Wilson B matrix multiplied by its transpose.
	! By definition, it has dimensions of n_prims x n_prims, so it is allocated accordingly.
	if (.not. ALLOCATED(Gmat)) allocate(Gmat(n_prims, n_prims))
	Gmat(:,:) = 0.0
	Gmat = MATMUL(Bmat_p, TRANSPOSE(Bmat_p))
	
	! By definition, the G matrix must have a determinant of zero. 
	! This criteria is checked to mitigate errors down the line.
	deter = DETERMINANT(Gmat, n_prims)
	if (deter .gt. 1E-6) then
		write (*,*) "WARNING; the determinant of the G matrix is not zero - either there is a problem, &
		& or there is exactly 3N-6 primitive coordinates."
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
	real(dp) :: Gmat(n_prims, n_prims), eigenvals(n_prims), eigenvecs(n_prims, n_prims), temp_eval, temp
	real(dp), allocatable :: Umat(:,:), Rmat(:,:)
	integer(i4b) :: IPIV(n_prims), INFO
	real(dp) :: WORK(n_prims*64), evecs_save(n_prims, n_prims), Gmat_test(n_prims, n_prims), evals_diag(n_prims, n_prims)
	real(dp) :: eigenvals_temp(n_prims), eigenvecs_temp(n_prims, n_prims)
	
	! By definition, the U matrix should contain 3N-6 (where N is the total number of atoms) eigenvectors.
	! The R matrix, should contain the remaining eigenvectors, which must be at least 6.
	! These matrices are allocated appropriately.
	if (.not. ALLOCATED(Umat)) allocate(Umat(n_prims, n_dlc))
	if (.not. ALLOCATED(Rmat)) allocate(Rmat(n_prims, (n_prims - n_dlc)))
	Umat(:,:) = 0.0
	Rmat(:,:) = 0.0

	! Extracting eigenvalues and eigenvectors...
	eigenvals = 0.0
	eigenvecs = 0.0
	eigenvals = EVALS(Gmat, n_prims)
	eigenvecs = EVECS(Gmat, n_prims)

	! Perform quick check to verify the consistency of the eigendecomposition.
	! We should be able to regain the G matrix by: Q * A * Q^(-1), where Q is the eigenvectors, and A is a diagonal matrix of eigenvalues.
	! The eigenvalues and eigenvectors are stored in increasing order.
	! This is not what we want, so these are flipped.
	eigenvecs_temp = eigenvecs(:,n_prims:1:-1)
	eigenvals_temp = eigenvals(n_prims:1:-1)

	! Save the eigenvectors and then make a diagonal matrix of the eigenvalues.
	evecs_save = 0.0
	evecs_save = eigenvecs_temp
	evals_diag = 0.0
	do j=1, SIZE(eigenvals_temp,1)
		evals_diag(j,j) = eigenvals_temp(j)
	end do
	
	! Diagonalise the eigenvectors.
	call DGETRF(n_prims, n_prims, eigenvecs_temp, n_prims, IPIV, INFO)
	call DGETRI(n_prims, eigenvecs_temp, n_prims, IPIV, WORK, 64*n_prims, INFO)
	
	! Now perform the check.
	Gmat_test = Gmat - MATMUL(MATMUL(evecs_save, evals_diag), eigenvecs_temp)
	if (MAXVAL(Gmat_test) .gt. 1E-08) then
		print *, 'Error; the eigendecomposition of the G matrix has gone wrong...'
		stop
	end if	
	
	! Lastly, create arrays for U and R.
	U_vector_counter = 1
	R_vector_counter = 1
	do j=1, n_prims
	    temp_eval = eigenvals(j)
		if (ABS(temp_eval) .lt. 1E-10) then
		    Rmat(:,R_vector_counter) = eigenvecs(:,j)
			R_vector_counter = R_vector_counter + 1
		else
		    Umat(:,U_vector_counter) = eigenvecs(:,j)
			U_vector_counter = U_vector_counter + 1
		end if
	end do

	end subroutine diag_Gmat
	
	
	subroutine gen_Bmat_DLC(atom_num, n_dlc, n_prims, Bmat_p, U_V_mat, Bmat_dlc)
	! Here, the Wilson B matrix in primitive internal coordinate subspace is updated to DLC subspace.
	!
	! ARGUMENTS:    atom_num   : Integer which represents the total number of atoms to be delocalised.
	!               n_dlc      : Integer which represents the number of delocalised internal coordinates (by definition, 3N-6).
	!               n_prims    : Integer which represents the total number of primitive internal coordinates.
	!               Bmat_p     : 2D array which contains the primitive Wilson B matrix.
	!				U_V_mat    : 2D array containing the eigenvectors with eigenvalues greater than zero. 
	!                            This is the DLC transformation vector set. 	
	!               Bmat_dlc   : 2D array which contains the DLC Wilson B matrix.
	
	implicit none
	integer(i4b) :: atom_num, n_prims, n_dlc
	real(dp) :: Bmat_p(n_prims, (3 * atom_num)), U_V_mat(n_prims, n_dlc)
	real(dp), allocatable :: Bmat_dlc(:,:)

	! The Wilson B matrix in DLC subspace is allocated. By definition, its dimensions are (3N-6) x (3N), where N is the number of atoms.
	if (.not. ALLOCATED(Bmat_dlc)) allocate(Bmat_dlc(n_dlc, (3 * atom_num)))
	Bmat_dlc(:,:) = 0.0

	! The calculation of the B matrix in DLC subspace is a straightforward multiplication.	
	Bmat_dlc = MATMUL(TRANSPOSE(U_V_mat), Bmat_p)
	
	end subroutine gen_Bmat_DLC
	
	
	subroutine gen_DLC(atom_num, n_dlc, n_prims, U_V_mat, prims, dlc)
	! Here, the DLC are actually generated from linear combinations of primitive internal coordinates and the U matrix.
	!
	! ARGUMENTS:  	atom_num : Integer which represents the total number of atoms to be delocalised.
	!               n_dlc    : Integer which represents the number of delocalised internal coordinates (by definition, 3N-6).
	!               n_prims  : Integer which represents the total number of primitive internal coordinates.		
	!               U_V_mat  : 2D array containing the eigenvectors with eigenvalues greater than zero. 
	!                          This is the DLC transformation vector set. 
	!               prims    : 1D array containing all the primitive internal coordinates associated with this input.				
	!				dlc      : 1D array containing the delocalised internal coordinate set.
	
	implicit none
	integer(i4b) :: n_prims, atom_num, n_dlc, i, j
	real(dp) :: U_V_mat(n_prims, n_dlc), prims(n_prims)
	real(dp), allocatable :: dlc(:)
	
	! The DLC matrix is allocated.
	if (.not. ALLOCATED(dlc)) allocate(dlc(n_dlc))
	dlc(:) = 0.0

	! The DLC are solved simply by linear combinations of the U matrix with the primitive interal coordinates.
	do i=1, n_dlc
		do j=1, n_prims
			dlc(i) = dlc(i) + (U_V_mat(j,i) * prims(j))
		end do
	end do

	end subroutine gen_DLC
	
	
	subroutine maintain_DLC(atom_num, n_dlc, coords)
	! Here, the DLC are maintained but not refreshed (i.e., the G matrix is not diagonalised!).
	! All arrays used in the generation are deallocated so that they can be allocated appropriately again.
	!
	! ARGUMENTS:    atom_num    : Integer which represents the total number of atoms to be delocalised.
	!               n_dlc       : Integer which represents the number of delocalised internal coordinates (by definition, 3N-6).
	!               coords      : 1D array containing all the cartesian coordinates of the system.
	
	implicit none
	integer(i4b) :: atom_num, n_dlc, i
	real(dp) :: coords(atom_num * 3)
	logical :: is_cons
	
	! First, (if necessary) deallocate all arrays used for DLC.
	if (ALLOCATED(prims)) deallocate(prims)
	if (ALLOCATED(Bmat_p)) deallocate(Bmat_p)
	if (ALLOCATED(Bmat_dlc)) deallocate(Bmat_dlc)
	if (ALLOCATED(dlc)) deallocate(dlc)
	
	! Now, check if there is any constraints...
	is_cons = .False.
	if (ncon_prim .gt. 0) then
		is_cons = .True.
	end if

	! Now, call DLC subroutines to re-populate the primitive arrays based on the Cartesian coordinates.
	call calc_prims(atom_num, nprim, prims, opt, coords, prim_list)
	call gen_Bmat_prims(atom_num, nprim, opt, coords, prim_list, Bmat_p)

	! The DLC are simply generated from linear combinations, and the Wilson B matrix converted to DLC subspace.
	if (is_cons .eqv. .False.) then
		! Generating delocalised internal coordinates...
		call gen_DLC(atom_num, n_dlc, nprim, Umat, prims, dlc)
		call gen_Bmat_DLC(atom_num, n_dlc, nprim, Bmat_p, Umat, Bmat_dlc)
	else if (is_cons .eqv. .True.) then
		! Generating delocalised internal coordinates...
		call gen_DLC(atom_num, n_dlc, nprim, Vmat, prims, dlc)
		call gen_Bmat_DLC(atom_num, n_dlc, nprim, Bmat_p, Vmat, Bmat_dlc)	
	end if

	end subroutine maintain_DLC
	
	subroutine refresh_DLC(atom_num, n_dlc, coords, cdat)
	! Here, the DLC are refreshed or generated for the first time.
	! All arrays used in the generation are deallocated so that they can be allocated appropriately again.
	!
	! ARGUMENTS:    atom_num    : Integer which represents the total number of atoms to be delocalised.
	!               n_dlc       : Integer which represents the number of delocalised internal coordinates (by definition, 3N-6).
	!               coords      : 1D array containing all the cartesian coordinates of the system.
	
	implicit none
	integer(i4b) :: atom_num, n_dlc, i
	real(dp) :: coords(atom_num * 3), cdat(ncon_prim, nprim)
	logical :: is_cons
	
	! First, (if necessary) deallocate all arrays used for DLC.
	if (ALLOCATED(prims)) deallocate(prims)
	if (ALLOCATED(Bmat_p)) deallocate(Bmat_p)
	if (ALLOCATED(Bmat_dlc)) deallocate(Bmat_dlc)
	if (ALLOCATED(dlc)) deallocate(dlc)
	
	! Now, (if required) deallocate the G and U + R matrices as they will be repopulated when refreshing the DLC.
	if (ALLOCATED(Gmat)) deallocate(Gmat)
	if (ALLOCATED(Umat)) deallocate(Umat)
	if (ALLOCATED(Rmat)) deallocate(Rmat)

	! Now, check if there is any constraints and deallocate the V matrix accordingly...
	is_cons = .False.
	if (ncon_prim .gt. 0) then
		is_cons = .True.
		if (ALLOCATED(Vmat)) then 
			deallocate(Vmat)
		end if
	end if

	! Now, all the DLC subroutines can be called to reallocate the arrays.
	call calc_prims(atom_num, nprim, prims, opt, coords, prim_list)
	call gen_Bmat_prims(atom_num, nprim, opt, coords, prim_list, Bmat_p)
	call gen_Gmat(atom_num, n_dlc, nprim, Bmat_p, Gmat)
	call diag_Gmat(atom_num, n_dlc, nprim, Gmat, Umat, Rmat)

	! If there is not any constraints, then it is relatively simple as the DLC are simply generated from linear combinations, and the Wilson B matrix converted to DLC subspace.
	if (is_cons .eqv. .False.) then
	
		! Generating delocalised internal coordinates...
		call gen_DLC(atom_num, n_dlc, nprim, Umat, prims, dlc)
		call gen_Bmat_DLC(atom_num, n_dlc, nprim, Bmat_p, Umat, Bmat_dlc)

	! If there are constraints, then the constraint must first be projected and then added to the active coordinate set.
	! Then, the set is Gram-Schmidt orthogonalised to regain the 3N-6 set of coordinates which contains ncon_prim constraint vectors.
	! These vectors can then be fixed in the optimisation and thus the constraint achieved.
	! The constraint vector(s) are moved to the back of the new active working matrix, Vmat, and the unprojected vectors are restored.
	! Lastly, the DLC can be generated, and the Wilson B matrix updated to DLC subspace.
	else if (is_cons .eqv. .True.) then
	
		! Generating constraints...
		call proj_cons(atom_num, n_dlc, ncon_prim, nprim, cdat, cdat_unproj, Umat)
		call gen_Vmat(atom_num, n_dlc, ncon_prim, nprim, cdat, Umat, Vmat)
		call ortho_mat(atom_num, n_dlc, ncon_prim, nprim, Vmat)
		call move_cons(atom_num, n_dlc, ncon_prim, nprim, Vmat)
		call unproj_cons(atom_num, n_dlc, ncon_prim, nprim, Vmat, cdat_unproj)

		! Generating delocalised internal coordinates...
		call gen_DLC(atom_num, n_dlc, nprim, Vmat, prims, dlc)
		call gen_Bmat_DLC(atom_num, n_dlc, nprim, Bmat_p, Vmat, Bmat_dlc)

	end if

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
	real(dp) :: Bmat_dlc(n_dlc, (3 * atom_num)), g(atom_num * 3)
	real(dp), allocatable :: g_dlc(:)
	real(dp) :: BT_inv(n_dlc, (3 * atom_num)), Gmat_dlc(n_dlc, n_dlc)
	
	! The primitive gradient array is allocated. By definition, its dimensions are the number of delocalised interal coordinates.
	if (.not. ALLOCATED(g_dlc)) allocate(g_dlc(n_dlc))
	g_dlc(:) = 0.0
	
	! Using single value decomposition, the Moore-Penrose inverse is constructed and used to convert the gradient array.
	Gmat_dlc = MATMUL(Bmat_dlc, TRANSPOSE(Bmat_dlc))
	BT_inv = MATMUL(SVD_INVERSE(Gmat_dlc, n_dlc, n_dlc), Bmat_dlc)
	g_dlc = MATMUL(BT_inv, g)

	end subroutine gen_grad_cart_to_DLC
	
	
	subroutine gen_hess_prim_to_DLC(atom_num, n_dlc, n_prims, U_V_mat, hess, hess_dlc)
	! Here, the cartesian hessian matrix is updated to DLC subspace.
	!
	! ARGUMENTS:    atom_num : Integer which represents the total number of atoms to be delocalised.
	!               n_dlc    : Integer which represents the number of delocalised internal coordinates (by definition, 3N-6).
	!               n_prims  : Integer which represents the total number of primitive internal coordinates.
	!				U_V_mat  : 2D array containing the eigenvectors with eigenvalues greater than zero. 
	!                          This is the DLC transformation vector set. 	
	!               hess     : 2D array which contains the hessian matrix in primitive subspace.
	!               hess_dlc : 2D array which contains the hessian matrix in DLC subspace.
	
	implicit none
	integer(i4b) :: atom_num, n_prims, n_dlc
	real(dp) :: U_V_mat(n_prims, n_dlc), hess(n_prims, n_prims)
	real(dp), allocatable :: hess_dlc(:,:)

	! First, the hessian matrix is DLC subspace should be allocated.
	if (.not. ALLOCATED(hess_dlc)) allocate(hess_dlc(n_dlc, n_dlc))
	hess_dlc(:,:) = 0.0

	! Now, the hessian can be calculated in DLC subspace by a simple multiplication procedure.
	hess_dlc = MATMUL(TRANSPOSE(U_V_mat), MATMUL(hess, U_V_mat))

	end subroutine gen_hess_prim_to_DLC
	
	
	subroutine DLC_to_cart_old(atom_num, n_dlc, n_prims, dS, dlc, x_1, x_2, Bmat_dlc)
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
	real(dp) :: dS_norm_save, dS_norm, dS_norm_init
	real(dp) :: x_1(3 * atom_num), x_2(3 * atom_num), dx(3 * atom_num)
	real(dp) :: dx_step(3 * atom_num), dx_temp(3 * atom_num), temp_x(3 * atom_num), dx_save(3 * atom_num)
	real(dp) :: BT_inv(n_dlc, (3 * atom_num)), Bmat_dlc(n_dlc, (3 * atom_num)), Gmat_dlc(n_dlc, n_dlc)
	real(dp) :: init_dlc(n_dlc), init_dS(n_dlc), target_dlc(n_dlc), dS(n_dlc)
	real(dp) :: dS_temp(n_dlc), dlc(n_dlc), dS_save(n_dlc)
	
    ! The expression B(transpose) * G(inverse) is initialised as it is used to convert between coordinate systems.
	Gmat_dlc = MATMUL(Bmat_dlc, TRANSPOSE(Bmat_dlc))
	BT_inv = MATMUL(SVD_INVERSE(Gmat_dlc, n_dlc, n_dlc), Bmat_dlc)
	
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
			call maintain_DLC(atom_num, n_dlc, temp_x)
			dS_temp = target_dlc - dlc
			dS_norm_save = NORM2(dS_temp)
			dx_save = dx_temp
			dS_save = dS_temp
		else
			dx_temp = dx_temp + dx_step
			temp_x = dx_temp + x_1
			call maintain_DLC(atom_num, n_dlc, temp_x)
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
	
	end subroutine DLC_to_cart_old
	
	subroutine DLC_to_cart(atom_num, n_dlc, n_prims, dq, q, x_1, x_2, INFO)
	! Here, the primitive internal coordinates are converted to cartesian coordinates using an iterative procedure.
	! This is used to generate the cartesian coordinates following a linear interpolation in primitive internal coordinates.
	!
	! ARGUMENTS:    atom_num : Integer which represents the total number of atoms to be delocalised.
	!               n_prims  : Integer which represents the total number of primitive internal coordinates.	
	!           	dq       : 1D array containing the change in primitive internal coordinates.
	!				q_1      : 1D array containing the primitive internal coordinate set of the starting point.  
    !               x_1      : 2D array containing all the cartesian coordinates of the starting point.	
	!				x_2      : 2D array containing all the cartesian coordinates after conversion.
	!				Bmat_p   : 2D array containing the Wilson B matrix used to convert between cartesian and primitive internal coordinates.
	!               prim_list   : 2D array containing the details of each primitive internal coordinate in the form ([1,2],[2,3], etc..).
	
	implicit none
	integer(i4b) :: atom_num, n_dlc, n_prims, i, k, iter_counter, counter
	real(dp) :: dq(n_dlc), q(n_dlc), q_2(n_dlc), x_1(3 * atom_num), x_2(3 * atom_num), init_x_1(3 * atom_num), init_init_dq(n_dlc)
	real(dp) :: BT_inv(n_dlc, (3 * atom_num)), temp_vec1(n_dlc), temp_vec2(n_dlc), temp_real, saved_best_x(3 * atom_num)
	real(dp) :: dq_actual(n_dlc), check(n_dlc), temp_check, xyz_rms_1, xyz_rms_2, old_dlc(n_dlc), temp_norm, scale_by
	real(dp) :: init_dq(n_dlc), target_q(n_dlc), dx(3 * atom_num), Gmat_dlc(n_dlc,n_dlc), saved_best_dlc(n_dlc)
	real(dp), allocatable :: temp_q(:)
	logical :: convergence, is_cons
	integer(i4b) :: IPIV(n_dlc), INFO
	real(dp) :: WORK(n_dlc*64)
	
    ! Since cartesians are rectilinear and internal coordinates are curvilinear, a simple transformation cannot be used.
    ! Instead, an iterative transformation procedure must be used.
    ! The expression G^(-1)* B is initialised as it is used to convert between coordinate systems.
	call maintain_DLC(atom_num, n_dlc, x_1)
	Gmat_dlc = MATMUL(Bmat_dlc, TRANSPOSE(Bmat_dlc))
	BT_inv = MATMUL(SVD_INVERSE(Gmat_dlc, n_dlc, n_dlc), Bmat_dlc)

	! Stashing some values for convergence criteria.
	convergence = .FALSE.
	xyz_rms_1 = 0
	xyz_rms_2 = 0
	INFO = 0
	init_x_1(:) = x_1(:)
	init_dq(:) = dq(:)
	init_init_dq(:) = init_dq(:)
	old_dlc(:) = dlc(:)
	target_q(:) = dlc(:) + init_dq(:)
	iter_counter = 0
	scale_by = 1.0

	! Checking if there are any constraints...
	is_cons = .False.
	if (ncon_prim .gt. 0) then
		is_cons = .True.
	end if
	
	do while (convergence .eqv. .FALSE.) 
		! The change in cartesian coordinates associated with the change in primitive internal coordinates is calculated.
		dx(:) = 0.0
		dx = MATMUL(TRANSPOSE(BT_inv), dq)

		! The root-mean-square change is used as a convergence criteria, so it is evaluated.
		xyz_rms_2 = RMSD_CALC(dx, x_1, (atom_num * 3))

		! The new cartesian geometry is evaluated, the new primitive internal coordinate set, Wilson B matrix, and DLC are obtained.
		x_2 = x_1 + dx
		call maintain_DLC(atom_num, n_dlc, x_2)
		
		! The transformation can become unstable with large changes in DLC.
		! In these cases, the first evaluation is probably the best guess as it will be in the right direction (at least).
		if (iter_counter .eq. 0) then
			saved_best_dlc(:) = dlc(:)
			saved_best_x(:) = x_2(:)
		end if

		! The Moore-Penrose inverse is constructed for the next iteration.
		BT_inv(:,:) = 0.0
		Gmat_dlc(:,:) = 0.0
		Gmat_dlc = MATMUL(Bmat_dlc, TRANSPOSE(Bmat_dlc))
		BT_inv = MATMUL(SVD_INVERSE(Gmat_dlc, n_dlc, n_dlc), Bmat_dlc)
		
        ! The change in primitive internals for the next iteration is evaluated.
		dq(:) = 0.0
		dq = target_q - dlc

		! Now, the three exit conditions should be checked...
		! The first ending condition for this transformation is when the root-mean-square change in cartesians is less than 10^-6.
        ! The second ending condition for this transformation is when the difference in root-mean-square change in cartesians between iteration i and i+1 is less than 10^-12.
        ! The third ending condition for this transformation is when the difference between the target DLC and the calculated DLC is less than 10^-6.
		check = target_q - dlc
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
		old_dlc(:) = dlc(:)
		xyz_rms_1 = xyz_rms_2

		! In some cases, cartesians cannot be solved, so an exit condition must exist for this case.
		iter_counter = iter_counter + 1
		if (iter_counter == 10000) then
			print *, "Error; could not solve cartesians from the change in delocalised internal coordinates. &
			& The change in delocalised internal coordinates was probably too large. &
			& Saving the most first evaluation as it's probably the best guess."
			INFO = 1
			dlc = saved_best_dlc
			x_1(:) = saved_best_x(:)
			exit
		end if
	end do

	! The new cartesian coordinates are saved.
	x_2(:) = x_1(:)

	end subroutine DLC_to_cart
	
	
END MODULE delocalised
