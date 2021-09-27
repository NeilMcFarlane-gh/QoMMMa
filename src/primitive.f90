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
	
	
	subroutine gen_prims2(bonds, atom_num, to_generate, prim_num, x, prims, prim_list)
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
	integer(i4b) :: atom_num, prim_num, i, j, k, l, m, n, prims_tot, prims_counter, alloc_counter, temp
	integer(i4b) :: bonds(atom_num,5), to_generate(prim_num)
	integer(i4b), allocatable :: prim_list(:,:), prims_temp(:,:)
	real(sp) :: x(3,atom_num), coords_1(3), coords_2(3), r, temp_prim, temp_prim2
	real(sp), allocatable :: prims(:), prims_all(:), prims_all2(:)
	real(sp), parameter :: cut_off = 2.5 !Angstroms
	
	! Before defining the primitive internal coordinates, the array prims_temp must be allocated.
	! For connecting every atom to every atom, allocating this is relatively straightforward.
	temp = SIZE(to_generate) * (SIZE(to_generate) - 1)
	allocate(prims_temp(temp, 2))
	
	! The array prims_temp is used to store the details of which prims are to be generated.
	prims_tot = 1
	do j=1, SIZE(to_generate)
	    i = to_generate(j)
		do l=1, SIZE(to_generate)
			m = to_generate(l)
			if (i .ne. m) then
				prims_temp(prims_tot,1) = i
				prims_temp(prims_tot,2) = m
				prims_tot = prims_tot + 1
			end if
		end do
	end do
	
	! Now, the primitive internal coordinate array, prims_all, can be populated with every primitive interal coordinate.
	! Then, any of these primitive internal coordinates which are greater than the predefined cut-off are set to 0.0 in the array prims_all.
	! The array prims_all2 is also created to store all the primitives for later use.
	allocate(prims_all(SIZE(prims_temp, 1)))
	allocate(prims_all2(SIZE(prims_temp, 1)))
	do j=1, SIZE(prims_temp, 1)
	    if ((prims_temp(j,1) .ne. 0) .and. (prims_temp(j,2) .ne. 0)) then
		    coords_1(:) = x(:, prims_temp(j,1))
		    coords_2(:) = x(:, prims_temp(j,2))
		    r = ATOM_DISTANCE(coords_1, coords_2)
		    if (r .lt. cut_off) then
			    prims_all(j) = r
				prims_all2(j) = r
		    else
			    prims_all(j) = 0.0
				prims_all2(j) = 0.0
		    end if
		end if
	end do

	! In the code above, primitives are often added in duplicate. i.e., (1,5) is the same as (5,1) - these duplicates are removed.
	do j=1, SIZE(prims_all)
		temp_prim = prims_all(j)
		do i=1, SIZE(prims_all2)
			temp_prim2 = prims_all2(i)
			if  (temp_prim == temp_prim2) then
				prims_all(j) = 0.0
				exit
			end if
		end do
	end do
	    
	! The output arrays prims and prim_list must be allocated so that they can be populated.
	alloc_counter = 0
	do j=1, SIZE(prims_all)
	    temp_prim = prims_all(j)
		if (temp_prim .gt. 0.0) then
		    alloc_counter = alloc_counter + 1
		end if
	end do
	allocate(prims(alloc_counter))
	allocate(prim_list(alloc_counter, 2))
	
	! Lastly, the output arrays prims and prim_list can be populated with appropriate details.
	prims_counter = 1
	do j=1, SIZE(prims_all)
	    temp_prim = prims_all(j)
		print *, temp_prim
		if (temp_prim .gt. 0.0) then
		    prims(prims_counter) = prims_all(j)
			prim_list(prims_counter,:) = prims_temp(j,:)
			prims_counter = prims_counter + 1
		end if
	end do
	print *, prim_list
	end subroutine gen_prims2
	
	
	subroutine gen_Bmat_prims(atom_num, x, prim_list, n_prims, Bmat_p)
	! Here, the Wilson B matrix for the conversion between cartesian and primitive internal coordinates is created.
	!
	! ARGUMENTS:    atom_num    : Integer which represents the total number of atoms in the system.
	!               x           : 2D array containing all the cartesian coordinates of the system.
	!               n_prims     : Integer which represents the total number of primitive internal coordinates.
	!               prim_list   : 2D array containing the details of each primitive internal coordinate in the form ([1,2],[2,3], etc..).
    !               Bmat_p      : 2D array which contains the primitive Wilson B matrix.	
	
	implicit none
	integer(i4b) :: i, j, k, number_reset, atom_num, n_prims, prim_list(n_prims, 2)
	real(sp) :: x(3, atom_num), grad(2,3), coords_1(3), coords_2(3)
	real(sp), allocatable :: Bmat_p(:,:)
	
	! The Wilson B matrix is allocated. By definition, its dimensions are (number of prims) x (3N), where N is the number of atoms.
	allocate(Bmat_p(n_prims, (3 * atom_num)))
	
	! The Wilson B matrix is populated with the relevant second derivative terms.
	number_reset = prim_list(1,1)
	do j=1, SIZE(prim_list, 1)
		coords_1(:) = x(:, prim_list(j,1))
		coords_2(:) = x(:, prim_list(j,2))	    
	    grad = ATOM_DISTANCE_GRAD(coords_1, coords_2)
		do i=1, 2
			if ((prim_list(j,i) - number_reset) == 0) then
			    k = 1
			else
				k = (3 * (prim_list(j,i) - number_reset))
			end if
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
