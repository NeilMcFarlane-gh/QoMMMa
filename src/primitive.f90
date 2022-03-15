MODULE primitive
use nrtype ; use coordinates ; use optimdata ; use math
implicit none

contains


	subroutine define_prims_TC(atom_num, to_generate, coords, prim_list)
	! Here, the primitive internal coordinates are defined for a given cartesian coordinate set using a total connectivity scheme with a distance cutoff.
	!
	! ARGUMENTS:    atom_num    : Integer which represents the total number of atoms to generate primitive internal coordinates for.
	!               to_generate : 1D array containing the list of atom indices which primitive internal coordinates are to be generated for.
	!               coords      : 1D array containing all the cartesian coordinates of the system.
	!               prim_list   : 2D array containing the details of each primitive internal coordinate in the form ([1,2],[2,3], etc..).

	implicit none
	integer(i4b) :: i, j, atom_num, prims_counter, alloc_counter, coord_counter_1, coord_counter_2, prim_type
	integer(i4b), allocatable :: prim_list(:,:), prim_list_temp(:,:), to_generate(:)
	real(sp) :: coords_1(3), coords_2(3), coords(atom_num * 3), r, temp_prim
	real(sp), allocatable :: prims_temp(:)
	real(sp), parameter :: cut_off = 25 !Angstroms

	! The list of primitive internal coordinates is rather simple to define when every atom is connected to every atom.
	prim_list_temp = COMBINATIONS_2(to_generate, SIZE(to_generate))
	
	! Next, using prim_list_temp, all the primitive coordinates can be generated, and their size can be checked.
	! Each of these distances are subsquently checked with a distance criteria so that the system makes chemical sense.
	! If they are greater than the predefined cut-off, then they are set to zero to allow for subsequent removal.
	if (.not. ALLOCATED(prims_temp)) allocate(prims_temp(SIZE(prim_list_temp, 1)))
	do j=1, SIZE(prim_list_temp, 1)
		coord_counter_1 = (prim_list_temp(j,1) - 1) * 3
		coord_counter_1 = coord_counter_1 + 1
		coord_counter_2 = (prim_list_temp(j,2) - 1) * 3
		coord_counter_2 = coord_counter_2 + 1
		coords_1(:) = coords(coord_counter_1:(coord_counter_1 + 2))
		coords_2(:) = coords(coord_counter_2:(coord_counter_2 + 2))
		r = ATOM_DISTANCE(coords_1, coords_2)
		if (r .lt. cut_off) then
			prims_temp(j) = r
		else
			prims_temp(j) = 0.0
		end if
	end do
	
	! The output arrays prims and prim_list must be allocated so that they can be populated.
	alloc_counter = 0
	do j=1, SIZE(prims_temp)
		temp_prim = prims_temp(j)
		if (temp_prim .gt. 0.0) then
			alloc_counter = alloc_counter + 1
		end if
	end do
	if (.not. ALLOCATED(prim_list)) allocate(prim_list(alloc_counter, 2))
	
	! The global integer, nprim, is also initialised here.
	nprim = alloc_counter
	
	! Lastly, the output arrays prims and prim_list can be populated.
	prim_list(:,:) = 0
	prims_counter = 1
	do j=1, SIZE(prims_temp)
		temp_prim = prims_temp(j)
		if (temp_prim .gt. 0.0) then
			prim_list(prims_counter,:) = prim_list_temp(j,:)
			prims_counter = prims_counter + 1
		end if
	end do
	print *, prim_list
	end subroutine define_prims_TC
	
	
	subroutine define_prims_full(atom_num, to_generate, coords, prim_list)
	! Here, the primitive internal coordinates are defined for a given cartesian coordinate set using a bond lengths, angles and torsions.
	!
	! ARGUMENTS:    atom_num    : Integer which represents the total number of atoms to generate primitive internal coordinates for.
	!               to_generate : 1D array containing the list of atom indices which primitive internal coordinates are to be generated for.
	!               coords      : 1D array containing all the cartesian coordinates of the system.
	!               prim_list   : 2D array containing the details of each primitive internal coordinate in the form ([1,2],[2,3], etc..).

	implicit none
	integer(i4b) :: i, j, atom_num, prims_counter, alloc_counter, coord_counter_1, coord_counter_2, prim_type
	integer(i4b), allocatable :: prim_list(:,:), prim_list_temp(:,:), to_generate(:)
	real(sp) :: coords_1(3), coords_2(3), coords(atom_num * 3), r, temp_prim
	real(sp), allocatable :: prims_temp(:)
	real(sp), parameter :: cut_off = 3.0 !Angstroms

	! The list of primitive internal coordinates is rather simple to define when every atom is connected to every atom.
	prim_list_temp = COMBINATIONS_2(to_generate, SIZE(to_generate))
	
	! Next, using prim_list_temp, all the primitive coordinates can be generated, and their size can be checked.
	! Each of these distances are subsquently checked with a distance criteria so that the system makes chemical sense.
	! If they are greater than the predefined cut-off, then they are set to zero to allow for subsequent removal.
	if (.not. ALLOCATED(prims_temp)) allocate(prims_temp(SIZE(prim_list_temp, 1)))
	do j=1, SIZE(prim_list_temp, 1)
		coord_counter_1 = (prim_list_temp(j,1) - 1) * 3
		coord_counter_1 = coord_counter_1 + 1
		coord_counter_2 = (prim_list_temp(j,2) - 1) * 3
		coord_counter_2 = coord_counter_2 + 1
		coords_1(:) = coords(coord_counter_1:(coord_counter_1 + 2))
		coords_2(:) = coords(coord_counter_2:(coord_counter_2 + 2))
		r = ATOM_DISTANCE(coords_1, coords_2)
		if (r .lt. cut_off) then
			prims_temp(j) = r
		else
			prims_temp(j) = 0.0
		end if
	end do
	
	! The output arrays prims and prim_list must be allocated so that they can be populated.
	alloc_counter = 0
	do j=1, SIZE(prims_temp)
		temp_prim = prims_temp(j)
		if (temp_prim .gt. 0.0) then
			alloc_counter = alloc_counter + 1
		end if
	end do
	if (.not. ALLOCATED(prim_list)) allocate(prim_list(alloc_counter, 2))
	
	! The global integer, nprim, is also initialised here.
	nprim = alloc_counter
	
	! Lastly, the output arrays prims and prim_list can be populated.
	prim_list(:,:) = 0
	prims_counter = 1
	do j=1, SIZE(prims_temp)
		temp_prim = prims_temp(j)
		if (temp_prim .gt. 0.0) then
			prim_list(prims_counter,:) = prim_list_temp(j,:)
			prims_counter = prims_counter + 1
		end if
	end do
	
	end subroutine define_prims_full
	
	subroutine calc_prims(atom_num, n_prims, prims, coords, prim_list) 
	! Here, the primitive internal coordinates are calculated from a given primitive coordinate set.
	!
	! ARGUMENTS:    atom_num    : Integer which represents the total number of atoms to generate primitive internal coordinates for.
	!               n_prims     : Integer which represents the total number of primitive internal coordinates.
	!               prims       : 1D array containing all the primitive internal coordinates associated with this input.
	!               coords      : 1D array containing all the cartesian coordinates of the system.
	!               prim_list   : 2D array containing the details of each primitive internal coordinate in the form ([1,2],[2,3], etc..).
	
	implicit none
	integer(i4b) :: i, j, atom_num, prims_counter, alloc_counter, n_prims
	integer(i4b) :: coord_counter_1, coord_counter_2, coord_counter_3, coord_counter_4
	integer(i4b), allocatable :: prim_list(:,:)
	real(sp) :: coords_1(3), coords_2(3), coords_3(3), coords_4(3), coords(atom_num * 3)
	real(sp) :: r, theta, phi, temp_prim
	real(sp), allocatable :: prims(:)

	! Using the earlier defined prim_list, each of the primitive internal coordinates can be calculated.
	if (.not. ALLOCATED(prims)) allocate(prims(n_prims))
	do j=1, n_prims
		if (SIZE(prim_list,j) .eq. 2) then
			! Stashing coordinates.
			coord_counter_1 = (prim_list(j,1) - 1) * 3
			coord_counter_1 = coord_counter_1 + 1
			coord_counter_2 = (prim_list(j,2) - 1) * 3
			coord_counter_2 = coord_counter_2 + 1
			coords_1(:) = coords(coord_counter_1:(coord_counter_1 + 2))
			coords_2(:) = coords(coord_counter_2:(coord_counter_2 + 2))
			
			! Calculating the primitive and adding to prims.
			r = ATOM_DISTANCE(coords_1, coords_2)
			prims(j) = r
		else if (SIZE(prim_list,j) .eq. 3) then
			! Stashing coordinates.
			coord_counter_1 = (prim_list(j,1) - 1) * 3
			coord_counter_1 = coord_counter_1 + 1
			coord_counter_2 = (prim_list(j,2) - 1) * 3
			coord_counter_2 = coord_counter_2 + 1
			coord_counter_3 = (prim_list(j,3) - 1) * 3
			coord_counter_3 = coord_counter_3 + 1
			coords_1(:) = coords(coord_counter_1:(coord_counter_1 + 2))
			coords_2(:) = coords(coord_counter_2:(coord_counter_2 + 2))
			coords_3(:) = coords(coord_counter_3:(coord_counter_3 + 2))
			
			! Calculating the primitive and adding to prims.
			theta = ATOM_ANGLE(coords_1, coords_2, coords_3)
			prims(j) = theta
		else if (SIZE(prim_list,j) .eq. 4) then
			! Stashing coordinates.
			coord_counter_1 = (prim_list(j,1) - 1) * 3
			coord_counter_1 = coord_counter_1 + 1
			coord_counter_2 = (prim_list(j,2) - 1) * 3
			coord_counter_2 = coord_counter_2 + 1
			coord_counter_3 = (prim_list(j,3) - 1) * 3
			coord_counter_3 = coord_counter_3 + 1
			coord_counter_4 = (prim_list(j,4) - 1) * 3
			coord_counter_4 = coord_counter_4 + 1
			coords_1(:) = coords(coord_counter_1:(coord_counter_1 + 2))
			coords_2(:) = coords(coord_counter_2:(coord_counter_2 + 2))
			coords_3(:) = coords(coord_counter_3:(coord_counter_3 + 2))
			coords_4(:) = coords(coord_counter_4:(coord_counter_4 + 2))
			
			! Calculating the primitive and adding to prims.
			phi = ATOM_DIHEDRAL(coords_1, coords_2, coords_3, coords_4)
			prims(j) = phi
		end if
	end do
	
	end subroutine calc_prims
	
	
	subroutine gen_Bmat_prims(atom_num, n_prims, coords, prim_list, Bmat_p)
	! Here, the Wilson B matrix for the conversion between cartesian and primitive internal coordinates is created.
	!
	! ARGUMENTS:    atom_num    : Integer which represents the total number of atoms to generate primitive internal coordinates for.
	!               n_prims     : Integer which represents the total number of primitive internal coordinates.
	!               coords      : 1D array containing all the cartesian coordinates of the system.
	!               prim_list   : 2D array containing the details of each primitive internal coordinate in the form ([1,2],[2,3], etc..).
    !               Bmat_p      : 2D array which contains the primitive Wilson B matrix.	
	
	implicit none
	integer(i4b) :: i, j, k, number_reset, atom_num, n_prims, prim_list(n_prims, 2), coord_counter_1, coord_counter_2
	real(sp) :: coords(atom_num * 3), grad(3,2), coords_1(3), coords_2(3)
	real(sp), allocatable :: Bmat_p(:,:)
	
	! The Wilson B matrix is allocated. By definition, its dimensions are (number of prims) x (3N), where N is the number of atoms.
	if (.not. ALLOCATED(Bmat_p)) allocate(Bmat_p((3 * atom_num), n_prims))
	Bmat_p(:,:) = 0.0
	
	! The Wilson B matrix is populated with the relevant second derivative terms.
	number_reset = prim_list(1,1)
	do j=1, SIZE(prim_list, 1)
		coord_counter_1 = (prim_list(j,1) - 1) * 3
		coord_counter_1 = coord_counter_1 + 1
		coord_counter_2 = (prim_list(j,2) - 1) * 3
		coord_counter_2 = coord_counter_2 + 1
		coords_1(:) = coords(coord_counter_1:(coord_counter_1 + 2))
		coords_2(:) = coords(coord_counter_2:(coord_counter_2 + 2))    
	    grad = ATOM_DISTANCE_GRAD(coords_1, coords_2)
		do i=1, 2
			if ((prim_list(j,i) - number_reset) == 0) then
			    k = 1
			else
				k = (3 * (prim_list(j,i) - number_reset)) + 1
			end if
			Bmat_p(k,j) = grad(1,i)
			Bmat_p(k+1,j) = grad(2,i)
			Bmat_p(k+2,j) = grad(3,i)
		end do
	end do

	end subroutine gen_Bmat_prims


	subroutine update_bfgs_p(atom_num, n_prims, hess_1, hess_2, g2_p, g1_p, prims_2, prims_1) 
	! Here, the primtive hessian matrix is updated using the BFGS method.
	! In the terminology in this subroutine, state 2 refers to the newly calculated values, and state 1 to the ones from the previous iteration.
	!
	! ARGUMENTS:    atom_num    : Integer which represents the total number of atoms to generate primitive internal coordinates for.
	!               n_prims     : Integer which represents the total number of primitive internal coordinates.
	!               hess_1      : 2D array which contains the hessian matrix in primitive subspace for state 1.
	!               hess_2      : 2D array which contains the hessian matrix in primitive subspace for state 2.
	!               g2_p        : 1D array containing the gradient in primitive subspace for state 2.
    !               g1_p        : 1D array containing the gradient in primitive subspace for state 1.
	!               prims_2     : 1D array containing all the primitive internal coordinates associated with state 2.
	!               prims_1     : 1D array containing all the primitive internal coordinates associated with state 1.
	
	implicit none
	integer(i4b) :: atom_num, n_prims, i, j
	real(sp) :: hess_1(n_prims, n_prims), hess_2(n_prims, n_prims), g2_p(n_prims), g1_p(n_prims)
	real(sp) :: dg(n_prims), dq(n_prims), dhess(n_prims, n_prims), prims_2(n_prims), prims_1(n_prims), dgdq, dqHdq

	! The changes associated with the gradient and coordinates from the current and previous iteration are calculated.
	dg = g2_p - g1_p
	dq = prims_2 - prims_1

	! To ensure that values in the BFGS update do not become unreasonably large, some scaling can be required.
	! First, initialise the values that may be scaled.
	dgdq = INNER_PRODUCT(dg, dq, n_prims)
	dqHdq = INNER_PRODUCT(MATMUL(dq, hess_1), dq, n_prims)
	if (dgdq < 0.001) dgdq = 0.001
	if (dqHdq < 0.001) dqHdq = 0.001

	! The change in the hessian can now be calculated.
	dhess = ((OUTER_PRODUCT(dg, dg, n_prims, n_prims)) / dgdq) &
	& - ((OUTER_PRODUCT(MATMUL(hess_1, dq), MATMUL(dq, hess_1), n_prims, n_prims))  / dqHdq)

	do i=1, n_prims
		do j=1, n_prims
			if (ABS(dhess(i,j)) .gt. 0.5) then
!				print *, "adjusting hessian in position...", i, j
				dhess(i,j) = h_p(i,j) * 0.0
			end if
		end do
	end do

	! Lastly, the newly updated hessian is obtained.
	hess_2 = hess_1 + dhess

	prims_1(:) = 0.0
	prims_2(:) = 0.0
	
	end subroutine update_bfgs_p
		
		
	subroutine gen_grad_cart_to_prim(atom_num, n_prims, Bmat_p, g, g_p)
	! Here, the gradient array in cartesian subspace is updated to primitive subspace.
	!
	! ARGUMENTS:    atom_num : Integer which represents the total number of atoms to be delocalised.
	!               n_prims  : Integer which represents the total number of primitive internal coordinates.
	!               Bmat_p   : 2D array which contains the primitive Wilson B matrix.
	!				g        : 1D array containing the gradient in cartesian subspace.
	!               g_p      : 1D array containing the gradient in primitive subspace.
	
	implicit none
	integer(i4b) :: atom_num, n_prims
	real(sp) :: Bmat_p((3 * atom_num), n_prims), g(atom_num * 3)
	real(sp), allocatable :: g_p(:)
	real(sp) :: BT_Ginv(n_prims, (3 * atom_num))

	! The primitive gradient array is allocated. By definition, its dimensions are the number of prims.
	if (.not. ALLOCATED(g_p)) allocate(g_p(n_prims))
	g_p(:) = 0.0

	! Using single value decomposition, the Moore-Penrose inverse is constructed and used to convert the gradient array.
	BT_Ginv = MATMUL(TRANSPOSE(Bmat_p), SVD_INVERSE(MATMUL(Bmat_p, TRANSPOSE(Bmat_p)), SIZE(Bmat_p, 1), SIZE(Bmat_p, 1)))
	g_p = MATMUL(g, TRANSPOSE(BT_Ginv))
	
	end subroutine gen_grad_cart_to_prim
	
	
	subroutine gen_hess_cart_to_prims
	! Here, the cartesian hessian matrix is updated to primitive internal coordinate subspace.
	
	end subroutine gen_hess_cart_to_prims
	
	
	subroutine prims_to_cart
	! Here, the primitive internal coordinates are converted to cartesian coordinates using an iterative procedure.
	
	end subroutine prims_to_cart
	

END MODULE primitive
