MODULE primitive
use nrtype ; use coordinates ; use optimdata
implicit none

contains


	subroutine gen_prims
	! Here, the primitive internal coordinates are generated for a given cartesian coordinate set using a total connectivity scheme with a distance cutoff.
	
	end subroutine gen_prims
	
	
	subroutine gen_Bmat_prims
	! Here, the Wilson B matrix for the conversion between cartesian and primitive internal coordinates is created.
	
	end subroutine gen_Bmat_prims

	
	subroutine gen_grad_cart_to_prims
	! Here, the gradient array in cartesian subspace is updated to primitive internal coordinate subspace.
	
	end subroutine gen_grad_cart_to_prims
	
	
	subroutine gen_hess_cart_to_prims
	! Here, the cartesian hessian matrix is updated to primitive internal coordinate subspace.
	
	end subroutine gen_hess_cart_to_prims
	
	
	subroutine prims_to_cart
	! Here, the primitive internal coordinates are converted to cartesian coordinates using an iterative procedure.
	
	end subroutine prims_to_cart
	

END MODULE primitive
