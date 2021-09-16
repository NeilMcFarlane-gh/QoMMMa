MODULE DLC
use nrtype ; use coordinates ; use optimdata
implicit none

contains


	subroutine gen_Gmat
	! Here, the G matrix used in the generatation of the final DLC is generated.
	
	end subroutine gen_Gmat
	
	
	subroutine diag_Gmat
	! Here, the G matrix used in the generatation of the final DLC is diagaonalised, and the resulting non-redundant (U matrix) and redundant (R matrix - not used) subspaces are separated.
	
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