MODULE primitive
use nrtype ; use coordinates ; use optimdata ; use math
implicit none

contains


	subroutine gen_prims(bonds, atom_num, to_generate, prim_num, x, prims, prim_list)
	! Here, the primitive internal coordinates are generated for a given cartesian coordinate set using a total connectivity scheme with a distance cutoff.
	!
	! ARGUMENTS:	bonds       : 2D array containing all the bonds found in the Tinker xyz file.
	!                             The atoms are in numerical order on a rowwise basis.
	!               atom_num    : Integer which represents the total number of atoms in the system.
	!               to_generate : 1D array containing the list of atom indices which primitive internal coordinates are to be generated for.
	!               prim_num    : Integer which represents the total number of atoms which primitive internal coordinates are to be generated for.
	!               x           : 2D array containing all the cartesian coordinates of the system.
	!               prims       : 1D array containing all the primitive internal coordinates associated with this input.
	!               prim_list   : 2D array containing the details of each primitive internal coordinate in the form ([1,2],[2,3], etc..).
	
	implicit none
	integer(i4b) :: atom_num, prim_num, i, j, k, n, bonds(atom_num,5), to_generate(prim_num), prims_temp(prim_num,5)
	integer(i4b), allocatable :: prim_list(:,:)
	real(sp) :: x(3,atom_num), coords_1(3), coords_2(3), r
	real(sp), allocatable :: prims(:)

	! The array prims_temp is used to store the details of which prims are to be generated.
	do j=1, SIZE(to_generate)
	    i = to_generate(j)
		prims_temp(j,:) = bonds(i,:)
	end do
	
	! Before defining the primitive internal coordinates, the array prim_list must be allocated.
	n = 0
	do j=1, SIZE(prims_temp, 1)
	    do i=1, SIZE(prims_temp, 2)
		    if (prims_temp(j,i) /= 0) then
			    n = n + 1
			end if
		end do
	end do
	allocate(prim_list(n, 2))
	
	! The array prim_list is populated with the details of primitive internal coordinates.
	! Columns 1 and 2 correspond to the atom indices of the primitive internal coordinates.
	k = 1
	do j=1, SIZE(prims_temp, 1)
	    do i=1, SIZE(prims_temp, 2)
		    if (prims_temp(j,i) /= 0) then
			    prim_list(k,1) = to_generate(j)
				prim_list(k,2) = prims_temp(j,i)
				k = k + 1
			end if
		end do
	end do

	! Now, the primitive internal coordinate array, prims, can be populated.
	! The items in the array prims are in the same order as that of prim_list.
	allocate(prims(SIZE(prim_list, 1)))
	do j=1, SIZE(prim_list, 1)
		coords_1(:) = x(:, prim_list(j,1))
		coords_2(:) = x(:, prim_list(j,2))
		r = ATOM_DISTANCE(coords_1, coords_2)
		prims(j) = r
	end do
	
	end subroutine gen_prims
	
	
	subroutine gen_Bmat_prims(atom_num, x, prim_list, n_prims, Bmat_p)
	! Here, the Wilson B matrix for the conversion between cartesian and primitive internal coordinates is created.
	!
	! ARGUMENTS:    atom_num    : Integer which represents the total number of atoms in the system.
	!               x           : 2D array containing all the cartesian coordinates of the system.
	!               n_prims     : Integer which represents the total number of primitive internal coordinates.
	!               prim_list   : 2D array containing the details of each primitive internal coordinate in the form ([1,2],[2,3], etc..).
    !               Bmat_p      : 2D array which contains the primitive Wilson B matrix.	
	
	implicit none
	integer(i4b) :: i, j, k, atom_num, n_prims, prim_list(n_prims, 2)
	real(sp) :: x(3, atom_num), grad(2,3), coords_1(3), coords_2(3)
	real(sp), allocatable :: Bmat_p(:,:)
	
	! The Wilson B matrix is allocated. By definition, its dimensions are (number of prims) x (3N), where N is the number of atoms.
	allocate(Bmat_p(n_prims, (3 * atom_num)))
	
	! The Wilson B matrix is populated with the relevant second derivative terms.
	do j=1, SIZE(prim_list, 1)
		coords_1(:) = x(:, prim_list(j,1))
		coords_2(:) = x(:, prim_list(j,2))	    
	    grad = ATOM_DISTANCE_GRAD(coords_1, coords_2)
		do i=1, 2
		    k = (3 * prim_list(j,i))
			print *, k
			Bmat_p(j,k:) = grad(i,:)
		end do
	end do
	
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
