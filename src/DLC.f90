MODULE DLC
use nrtype ; use coordinates ; use optimdata ; use math
implicit none

contains


	subroutine gen_Gmat(atom_num, n_prims, Bmat_p, Gmat)
	! Here, the G matrix used in the generatation of the DLC subspace is generated.
	!
	! ARGUMENTS:    atom_num    : Integer which represents the total number of atoms in the system.
	!               n_prims     : Integer which represents the total number of primitive internal coordinates.
    !               Bmat_p      : 2D array which contains the primitive Wilson B matrix.	
    !               Gmat        : 2D array which contains the G matrix used to generate the DLC subspace.	
	
	implicit none
	integer(i4b) :: atom_num, n_prims
	real(sp) :: Bmat_p((3 * atom_num), n_prims)
	real(sp), allocatable :: Gmat(:,:)

	! The G matrix is simply the Wilson B matrix multiplied by its transpose.
	! By definition, it has dimensions of n_prims x n_prims, so it is allocated accordingly.
	allocate(Gmat(n_prims, n_prims))
	Gmat = MATMUL(TRANSPOSE(Bmat_p), Bmat_p)

	end subroutine gen_Gmat
	
	
	subroutine diag_Gmat(atom_num, n_prims, Gmat, Umat, Rmat)
	! Here, the G matrix used in the generatation of the final DLC is diagaonalised, and the resulting non-redundant (U matrix) and redundant (R matrix - not used) subspaces are separated.
	!
	! ARGUMENTS:    atom_num    : Integer which represents the total number of atoms in the system.
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
	allocate(Umat(n_prims, ((3 * atom_num) - 6)))
	allocate(Rmat(n_prims, (n_prims - ((3 * atom_num) - 6))))

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
	
	
	subroutine gen_Bmat_DLC(atom_num, n_prims, Bmat_p, Umat, Bmat_S)
	! Here, the Wilson B matrix in primitive internal coordinate subspace is updated to DLC subspace.
	!
	! ARGUMENTS:    atom_num : Integer which represents the total number of atoms in the system.
	!               n_prims  : Integer which represents the total number of primitive internal coordinates.
	!               Bmat_p   : 2D array which contains the primitive Wilson B matrix.
	!				Umat     : 2D array containing the eigenvectors with eigenvalues greater than zero. 
	!                          This is the DLC transformation vector set. 	
	!               Bmat_S   : 2D array which contains the DLC Wilson B matrix.
	
	implicit none
	integer(i4b) :: atom_num, n_prims
	real(sp) :: Bmat_p((3 * atom_num), n_prims), Umat(n_prims, ((3 * atom_num) - 6))
	real(sp), allocatable :: Bmat_S(:,:)
	
	! The calculation of the B matrix in DLC subspace is a straightforward multiplication.
	allocate(Bmat_S((3 * atom_num), ((3 * atom_num) - 6)))
	Bmat_S = MATMUL(Umat, TRANSPOSE(Bmat_p))
	
	end subroutine gen_Bmat_DLC
	
	
	subroutine gen_DLC(Umat, prims, n_prims, atom_num, S)
	! Here, the DLC are actually generated from linear combinations of primitive internal coordinates and the U matrix.
	!
	! ARGUMENTS:  	Umat     : 2D array containing the eigenvectors with eigenvalues greater than zero. 
	!                          This is the DLC transformation vector set. 
	!               prims    : 1D array containing all the primitive internal coordinates associated with this input.
	!               n_prims  : Integer which represents the total number of primitive internal coordinates.	
	! 				atom_num : Integer which represents the total number of atoms in the system.
	!				S        : 1D array containing the delocalised internal coordinate set.
	
	implicit none
	integer(i4b) :: n_prims, atom_num, i, j
	real(sp) :: Umat(n_prims, ((3 * atom_num) - 6)), prims(n_prims)
	real(sp), allocatable :: S(:)
	
	! The number of DLC is equal to the number of eigenvectors in the U matrix.
	allocate(S(SIZE(Umat,2)))
	do i=1, SIZE(Umat,2)
		do j=1, SIZE(prims,1)
			S(i) = S(i) + (Umat(j,i) * prims(j))
		end do
	end do
	
	end subroutine gen_DLC
	
	
	subroutine gen_grad_cart_to_DLC
	! Here, the gradient array in cartesian subspace is updated to DLC subspace.
	
	end subroutine gen_grad_cart_to_DLC
	
	
	subroutine gen_hess_cart_to_DLC
	! Here, the cartesian hessian matrix is updated to DLC subspace.
	
	end subroutine gen_hess_cart_to_DLC
	
	
	subroutine gen_hess_prim_to_DLC
	! Here, the primitive internal coordinate hessian matrix is updated to DLC subspace.
	
	end subroutine gen_hess_prim_to_DLC
	
	
	subroutine DLC_to_cart(atom_num, n_prims, dS, S_1, x_1, Bmat_S)
	! Here, the DLC are converted to cartesian coordinates via an iterative procedure.
	! This is used at the end of an optimisation cycle to restore the cartesian coordinates for the next gradient evaluation.
	!
	! ARGUMENTS:    atom_num : Integer which represents the total number of atoms in the system.
	!               n_prims  : Integer which represents the total number of primitive internal coordinates.	
	!           	dS       : 1D array containing the change in DLC.
	!				S_1      : 1D array containing the delocalised internal coordinate set of the starting point.  
    !               x_1      : 2D array containing all the cartesian coordinates of the starting point.	
	!				Bmat_S   : 2D array containing the Wilson B matrix used to convert between cartesian and DLC.
	
	implicit none
	integer(i4b) :: atom_num, n_prims, to_generate(21), i
	real(sp) :: dS((3 * atom_num) - 6), S_1((3 * atom_num) - 6), x_1(3, atom_num) 
	real(sp) :: Bmat_S((3 * atom_num), n_prims), BT_Ginv(n_prims, ((3 * atom_num) - 6))
	real(sp) :: dS_n((3 * atom_num) - 6), S_2((3 * atom_num) - 6), x_2(3, atom_num)
	real(sp) :: init_dS((3 * atom_num) - 6), target_S((3 * atom_num) - 6), dx(3, atom_num)
	real(sp) :: xyz_rms_1, xyz_rms_2, iter_counter
	logical :: convergence
	
    ! Since cartesians are rectilinear and internal coordinates are curvilinear, a simple transformation cannot be used.
    ! Instead, an iterative transformation procedure must be used.
    ! The expression B(transpose) * G(inverse) is initialised as it is used to convert between coordinate systems.
	BT_Ginv = MATMUL(Bmat_S, SVD_INVERSE(MATMUL(TRANSPOSE(Bmat_S), Bmat_S), SIZE(Bmat_S, 1), SIZE(Bmat_S, 2)))
	
	! Stashing some values for convergence criteria.
	convergence = .FALSE.
	xyz_rms_1 = 0
	xyz_rms_2 = 0
	init_dS(:) = dS(:)
	target_S(:) = S_1(:) + init_dS(:)
	
	do while (convergence .eqv. .FALSE.)
		! The change in cartesian coordinates associated with the change in DLC is calculated.
		dx = MATMUL(BT_Ginv, x_1)
		
		! The root-mean-square change is used as a convergence criteria, so it is evaluated.
		!xyz_rms_2 = RMSD_CALC(dx, x_1, SIZE(x_1, 1))
		
		! The new cartesian geometry is evaluated, and the new DLC geometry is obtained.
		x_2 = x_1 + dx
		to_generate = (/(i, i=1,21, 1)/) ! List of atoms to delocalise - for testing purposes...
		!call gen_prims(SIZE(x_2,2), to_generate, SIZE(to_generate), x, prims, prim_list)
		!call gen_Bmat_prims(SIZE(to_generate), x, prim_list, SIZE(prim_list, 1), Bmat_p)
		!call gen_Gmat(SIZE(to_generate), SIZE(prim_list,1), Bmat_p, Gmat)
		!call diag_Gmat(SIZE(to_generate), SIZE(prim_list,1), Gmat, Umat, Rmat)
		!call gen_DLC(Umat, prims, SIZE(prim_list,1), SIZE(to_generate), S)
		!call gen_Bmat_DLC(SIZE(to_generate), SIZE(prim_list,1), Bmat_p, Umat, Bmat_S)
	
	end do
	
	end subroutine DLC_to_cart
	
	
END MODULE DLC