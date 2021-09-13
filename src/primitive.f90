MODULE primitive
use nrtype ; use coordinates ; use optimdata
implicit none

! Maybe have some declarations here...

contains


	subroutine gen_prims
	! Here, the primitive internal coordinates are generated for a given cartesian coordinate set using a total connectivity scheme with a distance cutoff.
	! TO-DO : Investigate Tinker's program files to grep covalent radii.
	
	end subroutine gen_prims
	
	
	subroutine gen_Bmat_prims
	! Here, the Wilson B matrix for the conversion between cartesian and primitive internal coordinates is created.
	
	end subroutine gen_Bmat_prims

	
	subroutine gen_grad_cartprims
	! Here, the gradient array in cartesian subspace is updated to primitive internal coordinate subspace.
	
	end subroutine gen_grad_cartprims
	
	
	subroutine gen_hess_cartprims
	! Here, the cartesian hessian matrix is updated to primitive internal coordinate subspace.
	
	end subroutine gen_hess_cartprims
	
	
	subroutine prims_to_cart
	! Here, the primitive internal coordinates are converted to cartesian coordinates using an iterative procedure.
	
	end subroutine prims_to_cart
	

END MODULE primitive
