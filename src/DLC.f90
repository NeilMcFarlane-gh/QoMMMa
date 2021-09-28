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
    !               Gmat        : 2D array which contains the G matrix used to generated the DLC subspace.	
	
	implicit none
	integer(i4b) :: atom_num, n_prims
	real(sp) :: Bmat_p(n_prims, (3 * atom_num)), det
	real(sp), allocatable :: Gmat(:,:)
	
	! The G matrix is allocated. By definition, its dimensions are (number of prims) x (number of prims).	
	allocate(Gmat(n_prims, n_prims))
	
	! The G matrix is calculated simply as the B matrix multiplied by its transpose.
	Gmat = MATMUL(TRANSPOSE(Bmat_p), Bmat_p)

	! By definition, the G matrix is singular and so has a determinant of zero.
	! This criteria is checked to minimise errors down the line.
	det = DETERMINANT(Gmat, SIZE(Gmat,1))

	if (det .gt. 0.000) then
	    print *, "WARNING: The G Matrix is not singlular. This can either mean an entirely non-redundant primitive set, or some error."
	end if
	
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
	print *, "total U and R", ((3 * atom_num) - 6), (n_prims - ((3 * atom_num) - 6))

	! Extracting eigenvalues and eigenvectors...
	eigenvals = EVALS(Gmat, SIZE(Gmat,1))
	eigenvecs = EVECS(Gmat, SIZE(Gmat,1))
	
	U_vector_counter = 1
	R_vector_counter = 1
	do j=1, SIZE(eigenvals)
	    temp_eval = eigenvals(j)
		if (ABS(temp_eval) .lt. 1E-01) then
			print *, "redundant...", ABS(temp_eval), R_vector_counter
		    Rmat(:,R_vector_counter) = eigenvecs(:,j)
			R_vector_counter = R_vector_counter + 1
		else
			print *, "non-redundant...", ABS(temp_eval), U_vector_counter
		    Umat(:,U_vector_counter) = eigenvecs(:,j)
			U_vector_counter = U_vector_counter + 1
		end if
	end do
	print *, "done!"
	!print *, "Umat...", Umat
	!print *, "Rmat...", Rmat	
	
	end subroutine diag_Gmat
	
	
	subroutine gen_Bmat_DLC
	! Here, the Wilson B matrix in primitive internal coordinate subspace is updated to DLC subspace.
	
	end subroutine gen_Bmat_DLC
	
	
	subroutine gen_DLC
	! Here, the DLC are actually generated from linear combinations of primitive internal coordinates and the U matrix.

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
	
	
	subroutine DLC_to_cart
	! Here, the DLC are converted to cartesian coordinates using an iterative procedure.
	
	end subroutine DLC_to_cart
	
	
END MODULE DLC