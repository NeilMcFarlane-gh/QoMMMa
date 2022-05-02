MODULE primitive
use nrtype ; use coordinates ; use optimdata ; use math
implicit none

contains


	subroutine define_prims_TC(atom_num, to_generate, prim_list)
	! Here, the primitive internal coordinates are defined for a given cartesian coordinate set using a total connectivity scheme.
	! This method of primitive internal coordinate definition works well in complex ring systems, as typical primitive internal coordinates struggle with these systems.
	!
	! ARGUMENTS:    atom_num    : Integer which represents the total number of atoms to generate primitive internal coordinates for.
	!               to_generate : 1D array containing the list of atom indices which primitive internal coordinates are to be generated for.
	!               prim_list   : 2D array containing the details of each primitive internal coordinate in the form ([1,2],[2,3], etc..).

	implicit none
	integer(i4b) :: atom_num, nprim, to_generate(atom_num)
	integer(i4b), allocatable :: prim_list(:,:)

	! The list of primitive internal coordinates is rather simple to define when every atom is connected to every atom.
	prim_list = COMBINATIONS_2(to_generate, atom_num)
	
	! The global integer, nprim, is also initialised here.
	nprim = SIZE(prim_list, 1)

	end subroutine define_prims_TC
	
	subroutine define_prims_full(atom_num, to_generate, coords, prim_list) 
	! Here, the primitive internal coordinates are defined for a given cartesian coordinate set using bond lengths, angles and dihedral torsions.
	!
	! ARGUMENTS:    atom_num      : Integer which represents the total number of atoms to generate primitive internal coordinates for.
	!               to_generate   : 1D array containing the list of atom indices which primitive internal coordinates are to be generated for.
	!               coords        : 1D array containing all the cartesian coordinates of the system.
	!               prim_list     : 2D array containing the details of each primitive internal coordinate in the form ([1,2],[2,3,4],[2,3,4,5] etc..).

	implicit none
	integer(i4b) :: i, j, k, m, atom_num, to_generate(atom_num), neigh_pos
	integer(i4b) :: coord_counter_1, coord_counter_2, opt_id, neighbours(atom_num,maxbond), total_prims
	integer(i4b) :: actual_bond(2), actual_bond2(2), actual_ang(3), actual_ang2(3), possible_di(4), possible_di2(4)
	integer(i4b) :: atom_1, atom_2, prim_counter, ang_counter, di_atom, temp_prim(4), temp_prim_1(4), temp_prim_2(4)
	integer(i4b), allocatable :: TC_list(:,:), prim_list(:,:), prim_list_ordered(:,:), angle_atoms(:)
	integer(i4b), allocatable :: prim_list_temp(:,:), combo_list(:,:)
	real(sp) :: coords_1(3), coords_2(3), coords(atom_num * 3), r, dist
	real(sp), allocatable :: distances(:)
	logical :: is_dupe = .FALSE.
	logical :: is_allocated = .FALSE.

	! This algorithm starts very similarly to define_prims_TC, where every atom is connected to every atom.
	! Then, the distances between every atom can be calculated so that the bonding can be established.
	TC_list = COMBINATIONS_2_DUPE(to_generate, atom_num)
	allocate(distances(SIZE(TC_list,1)))
	do i=1, SIZE(TC_list,1)
		! Stashing coordinates.
		coord_counter_1 = (TC_list(i,1) - 1) * 3
		coord_counter_1 = coord_counter_1 + 1
		coord_counter_2 = (TC_list(i,2) - 1) * 3
		coord_counter_2 = coord_counter_2 + 1
		coords_1(:) = coords(coord_counter_1:(coord_counter_1 + 2))
		coords_2(:) = coords(coord_counter_2:(coord_counter_2 + 2))

		! Calculating the distance and adding to distances.
		r = ATOM_DISTANCE(coords_1, coords_2)
		distances(i) = r
	end do

	! Now, an array containing the neighbours for each atom can be initialised.
	! Here, a predefined cut-off determines whether two atoms are bonded or not.
	! TO-DO : This should be based on the VdW radius so that larger atoms (e.g., Fe) can have neighbouring atoms.
	neighbours = 0
	do i=1, SIZE(TC_list,1)
		dist = distances(i)
		if (dist .lt. cut_off) then
			opt_id = TC_list(i,1)
			do j=1, maxbond
				neigh_pos = neighbours(opt_id,j)
				if (neigh_pos .eq. 0) then
					neighbours(opt_id,j) = TC_list(i,2)
					exit
				end if
			end do
		end if
	end do

	! Due to the fact that the number of primitive internal coordinates is unknown at initial execution, the algorithm must be performed twice.
	! Indeed, this is not computationally efficient, but this avoids problems with memory allocation.
	do m=1, 2
		! The primitive counter is reset.
		prim_counter = 1
		
		! First, define all bonds. This is simple as it has already been defined in the array neighbours.
		do i=1, atom_num
			atom_1 = i
			do j=1, maxbond
				atom_2 = neighbours(i,j)
				if (atom_2 .ne. 0) then
					if (is_allocated .eqv. .FALSE.) then
						prim_counter = prim_counter + 1
					else
						actual_bond(1) = atom_1
						actual_bond(2) = atom_2
						prim_list_temp(prim_counter,1:2) = actual_bond
						actual_bond2 = ASCEND_ORDER(actual_bond,2)
						prim_list_ordered(prim_counter,1:2) = actual_bond2
						prim_counter = prim_counter + 1
					end if
				end if
			end do
		end do
		
		! Now, the angles and dihedral torsions can can be defined.
		! Angles are relatively easily to define as they are essentially given in the array neighbours.
		! Dihedrals are not directly defined in neighbours but angles are used as for a dihedral to exist, three atoms must be bonded already.
		do i=1, atom_num
			atom_1 = i
			
			! Firstly, place all the atoms that atom i is bonded to into the array angle_atoms.
			ang_counter = 0
			do j=1, maxbond
				if (neighbours(i,j) .ne. 0) then
					ang_counter = ang_counter + 1
				end if
			end do
			if (ALLOCATED(angle_atoms)) deallocate(angle_atoms)
			allocate(angle_atoms(ang_counter))
			ang_counter = 0
			do j=1, maxbond
				if (neighbours(i,j) .ne. 0) then
					ang_counter = ang_counter + 1
					angle_atoms(ang_counter) = neighbours(i,j)
				end if
			end do
			
			! If there is only one or zero atom bonded to atom i, then no angles are possible so the loop is exited.
			if (SIZE(angle_atoms) .le. 1) then
				exit
			end if
			
			! Rather like defining a total connectivity scheme, we can use COMBINATIONS_2 of the atoms which atom_i is bonded to and generate the possible angles.
			combo_list = COMBINATIONS_2(angle_atoms, SIZE(angle_atoms))
			
			! Lastly, we can generate each angle definition by iterating through combo_list.
			! For the calculation of the angle, atom i is placed in the middle of the definition.
			do j=1, SIZE(combo_list,1)
				actual_ang(1) = combo_list(j,1)
				actual_ang(2) = atom_1
				actual_ang(3) = combo_list(j,2)
				if (is_allocated .eqv. .FALSE.) then
					prim_counter = prim_counter + 1
					
					! Now checking if any dihedral torsions exist for this angle.
					do k=1, maxbond
						di_atom = neighbours(actual_ang(1),k)
						if ((di_atom .ne. 0) .AND. (di_atom .ne. actual_ang(2)) .AND. (di_atom .ne. actual_ang(3))) then
							prim_counter = prim_counter + 1
						di_atom = neighbours(actual_ang(3),k)
						else if ((di_atom .ne. 0) .AND. (di_atom .ne. actual_ang(1)) .AND. (di_atom .ne. actual_ang(2))) then
							prim_counter = prim_counter + 1
						end if
					end do
				else
					prim_list_temp(prim_counter,1:3) = actual_ang
					actual_ang2 = ASCEND_ORDER(actual_ang, 3)
					prim_list_ordered(prim_counter,1:3) = actual_ang2
					prim_counter = prim_counter + 1
					
					! Now checking if any dihedral torsions exist for this angle.
					do k=1, maxbond
						di_atom = neighbours(actual_ang(1),k)
						if ((di_atom .ne. 0) .AND. (di_atom .ne. actual_ang(2)) .AND. (di_atom .ne. actual_ang(3))) then
							possible_di(1) = di_atom
							possible_di(2) = actual_ang(1)
							possible_di(3) = actual_ang(2)
							possible_di(4) = actual_ang(3)
							prim_list_temp(prim_counter,:) = possible_di
							possible_di2 = ASCEND_ORDER(possible_di, 4)
							prim_list_ordered(prim_counter,:) = possible_di2
							prim_counter = prim_counter + 1
						end if
					end do
					do k=1, maxbond
						di_atom = neighbours(actual_ang(3),k)
						if ((di_atom .ne. 0) .AND. (di_atom .ne. actual_ang(2)) .AND. (di_atom .ne. actual_ang(1))) then
							possible_di(1) = actual_ang(1)
							possible_di(2) = actual_ang(2)
							possible_di(3) = actual_ang(3)
							possible_di(4) = di_atom
							prim_list_temp(prim_counter,:) = possible_di
							possible_di2 = ASCEND_ORDER(possible_di, 4)
							prim_list_ordered(prim_counter,:) = possible_di2
							prim_counter = prim_counter + 1
						end if
					end do
				end if
			end do
		end do
		
		! Now, the primitive interal coordinate arrays can be allocated and the second iteration of the cycle performed.
		if (is_allocated .eqv. .FALSE.) then
			allocate(prim_list_temp(prim_counter-1,4))
			allocate(prim_list_ordered(prim_counter-1,4))
			prim_list_temp = 0
			prim_list_ordered = 0
			is_allocated = .TRUE.
		end if
	end do

	! The algorithm is not perfect, so can produce duplicates. These are first located.
	total_prims = SIZE(prim_list_temp,1)
	do i=1, SIZE(prim_list_temp,1)
		temp_prim_1(:) = prim_list_ordered(i,:)
		do j=1, SIZE(prim_list_temp,1)
			temp_prim_2(:) = prim_list_ordered(j,:)
			if ((i .ne. j) .AND. (SUM(temp_prim_1) .ne. 0) .AND. (SUM(temp_prim_2) .ne. 0)) then
				temp_prim = temp_prim_1 - temp_prim_2
				do k=1, 4
					if (ABS(temp_prim(k)) .gt. 0) then
						is_dupe = .FALSE.
						exit
					else if (k == 4) then
						is_dupe = .TRUE.
					end if
				end do
				if (is_dupe .eqv. .TRUE.) then
					prim_list_ordered(j,:) = 0
					prim_list_temp(j,:) = 0
					total_prims = total_prims - 1
				end if
			end if
		end do
	end do

	! Lastly, the final set of primitive internal coordinates is placed in an array.
	allocate(prim_list(total_prims,4))
	prim_list = 0
	prim_counter = 1
	do i=1, SIZE(prim_list_temp,1)
		if (prim_list_temp(i,1) .ne. 0) then
			prim_list(prim_counter,:) = prim_list_temp(i,:)
			prim_counter = prim_counter + 1
		end if
	end do
	
	! The global integer, nprim, is also initialised here.
	nprim = SIZE(prim_list,1)

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
	integer(i4b) :: j, k, atom_num, n_prims, zero_count, temp_int
	integer(i4b) :: coord_counter_1, coord_counter_2, coord_counter_3, coord_counter_4
	integer(i4b) :: prim_list(n_prims, 4), work_prim(4)
	real(sp) :: coords_1(3), coords_2(3), coords_3(3), coords_4(3), coords(atom_num * 3)
	real(sp) :: r, theta, phi
	real(sp), allocatable :: prims(:)

	! First, allocate (if necessary) and zero the primitive internal coordinates.
	if (.not. ALLOCATED(prims)) allocate(prims(n_prims))
	prims(:) = 0.0

	! Using the earlier defined prim_list, each of the primitive internal coordinates can be calculated.
	do j=1, n_prims
		! First, decide which kind of primitive internal coordinate we are dealing with.
		zero_count = 0
		work_prim = prim_list(j,:)
		do k=1, 4
			temp_int = work_prim(k)
			if (temp_int .eq. 0) then
				zero_count = zero_count + 1
			end if
		end do
		if (zero_count .eq. 2) then
			! Stashing coordinates.
			coord_counter_1 = (prim_list(j,1) - 1) * 3
			coord_counter_1 = coord_counter_1 + 1
			coord_counter_2 = (prim_list(j,2) - 1) * 3
			coord_counter_2 = coord_counter_2 + 1
			coords_1(:) = coords(coord_counter_1:(coord_counter_1 + 2))
			coords_2(:) = coords(coord_counter_2:(coord_counter_2 + 2))
			
			! Calculating the primitive and adding to prims.
			r = 0.0
			r = ATOM_DISTANCE(coords_1, coords_2)
			prims(j) = r
		else if (zero_count .eq. 1) then
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
			theta = 0.0
			theta = ATOM_ANGLE(coords_1, coords_2, coords_3)
			prims(j) = theta
		else if (zero_count .eq. 0) then
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
			phi = 0.0
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
	integer(i4b) :: i, j, k, number_reset, atom_num, n_prims, zero_count, temp_int
	integer(i4b) :: coord_counter_1, coord_counter_2, coord_counter_3, coord_counter_4
	integer(i4b) :: prim_list(n_prims, 4), work_prim(4)
	real(sp) :: grad_r(3,2), grad_theta(3,3), grad_phi(3,4)
	real(sp) :: coords_1(3), coords_2(3), coords_3(3), coords_4(3), coords(atom_num * 3)
	real(sp), allocatable :: Bmat_p(:,:)
	
	! The Wilson B matrix is allocated. By definition, its dimensions are (number of prims) x (3N), where N is the number of atoms.
	if (.not. ALLOCATED(Bmat_p)) allocate(Bmat_p((3 * atom_num), n_prims))
	Bmat_p(:,:) = 0.0
	
	! The Wilson B matrix is populated with the relevant second derivative terms.
	number_reset = 1
	do j=1, n_prims
		! First, decide which kind of primitive internal coordinate we are dealing with.
		zero_count = 0
		work_prim = prim_list(j,:)
		do k=1, 4
			temp_int = work_prim(k)
			if (temp_int .eq. 0) then
				zero_count = zero_count + 1
			end if
		end do
		if (zero_count .eq. 2) then
			! Stashing coordinates.
			coord_counter_1 = (prim_list(j,1) - 1) * 3
			coord_counter_1 = coord_counter_1 + 1
			coord_counter_2 = (prim_list(j,2) - 1) * 3
			coord_counter_2 = coord_counter_2 + 1
			coords_1(:) = coords(coord_counter_1:(coord_counter_1 + 2))
			coords_2(:) = coords(coord_counter_2:(coord_counter_2 + 2))  

			! Calculating the analytical first derivatives and adding them to the Wilson B matrix.
			grad_r = ATOM_DISTANCE_GRAD(coords_1, coords_2)
			do i=1, 2
				if ((prim_list(j,i) - number_reset) == 0) then
					k = 1
				else
					k = (3 * (prim_list(j,i) - number_reset)) + 1
				end if
				Bmat_p(k,j) = grad_r(1,i)
				Bmat_p(k+1,j) = grad_r(2,i)
				Bmat_p(k+2,j) = grad_r(3,i)
			end do
		else if (zero_count .eq. 1) then
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

			! Calculating the analytical first derivatives and adding them to the Wilson B matrix.
			grad_theta = ATOM_ANGLE_GRAD(coords_1, coords_2, coords_3)
			do i=1, 3
				if ((prim_list(j,i) - number_reset) == 0) then
					k = 1
				else
					k = (3 * (prim_list(j,i) - number_reset)) + 1
				end if
				Bmat_p(k,j) = grad_theta(1,i)
				Bmat_p(k+1,j) = grad_theta(2,i)
				Bmat_p(k+2,j) = grad_theta(3,i)
			end do
		else if (zero_count .eq. 0) then
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

			! Calculating the analytical first derivatives and adding them to the Wilson B matrix.
			grad_phi = ATOM_DIHEDRAL_GRAD(coords_1, coords_2, coords_3, coords_4)
			do i=1, 4
				if ((prim_list(j,i) - number_reset) == 0) then
					k = 1
				else
					k = (3 * (prim_list(j,i) - number_reset)) + 1
				end if
				Bmat_p(k,j) = grad_phi(1,i)
				Bmat_p(k+1,j) = grad_phi(2,i)
				Bmat_p(k+2,j) = grad_phi(3,i)
			end do
		end if
	end do

	end subroutine gen_Bmat_prims


	subroutine update_bfgs_p(atom_num, n_prims, hess, g2_p, g1_p, prims_2, prims_1) 
	! Here, the primtive hessian matrix is updated using the BFGS method.
	! In the terminology in this subroutine, state 2 refers to the newly calculated values, and state 1 to the ones from the previous iteration.
	!
	! ARGUMENTS:    atom_num    : Integer which represents the total number of atoms to generate primitive internal coordinates for.
	!               n_prims     : Integer which represents the total number of primitive internal coordinates.
	!               hess        : 2D array which contains the hessian matrix in primitive subspace.
	!               g2_p        : 1D array containing the gradient in primitive subspace for state 2.
    !               g1_p        : 1D array containing the gradient in primitive subspace for state 1.
	!               prims_2     : 1D array containing all the primitive internal coordinates associated with state 2.
	!               prims_1     : 1D array containing all the primitive internal coordinates associated with state 1.
	
	implicit none
	integer(i4b) :: atom_num, n_prims, i, j
	real(sp) :: hess(n_prims, n_prims), g2_p(n_prims), g1_p(n_prims)
	real(sp) :: dg(n_prims), dq(n_prims), dhess(n_prims, n_prims), prims_2(n_prims), prims_1(n_prims), dgdq, dqHdq

	! The changes associated with the gradient and coordinates from the current and previous iteration are calculated.
	dg = g2_p - g1_p
	dq = prims_2 - prims_1

	! To ensure that values in the BFGS update do not become unreasonably large, some scaling can be required.
	! First, initialise the values that may be scaled.
	dgdq = INNER_PRODUCT(dg, dq, n_prims)
	dqHdq = INNER_PRODUCT(MATMUL(dq, hess), dq, n_prims)
	if (dgdq < 0.001) dgdq = 0.001
	if (dqHdq < 0.001) dqHdq = 0.001

	! The change in the hessian can now be calculated.
	dhess = ((OUTER_PRODUCT(dg, dg, n_prims, n_prims)) / dgdq) &
			- ((OUTER_PRODUCT(MATMUL(hess, dq), MATMUL(dq, hess), n_prims, n_prims))  / dqHdq)
	
    ! In some cases, the change in Hessian becomes unreasonably small. 
    ! When this is the case, the magnitude of the value is increased, but it's direction (positive or negative) remains the same.
	! In other cases, the change in Hessian is unreasonably large.
	! When this is the case, the change in the given Hessian element is simply set to zero.
	do i=1, n_prims
		do j=1, n_prims
			if (ABS(dhess(i,j)) .lt. 0.05) then
				dhess(i,j) = dhess(i,j) * 2
			else if (ABS(dhess(i,j)) .gt. 2) then
				dhess(i,j) = 0.0
			end if
		end do
	end do

	! Lastly, the newly updated hessian is obtained.
	hess = hess + dhess
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
	real(sp) :: BT_Ginv(n_prims, (3 * atom_num)), Gmat(n_prims, n_prims)

	! The primitive gradient array is allocated. By definition, its dimensions are the number of prims.
	if (.not. ALLOCATED(g_p)) allocate(g_p(n_prims))
	g_p(:) = 0.0
	
	! Using single value decomposition, the Moore-Penrose inverse is constructed and used to convert the gradient array.
	Gmat = MATMUL(TRANSPOSE(Bmat_p), Bmat_p)
	BT_Ginv = MATMUL(SVD_INVERSE(Gmat, n_prims, n_prims), TRANSPOSE(Bmat_p))
	g_p = MATMUL(g, TRANSPOSE(BT_Ginv))

	end subroutine gen_grad_cart_to_prim
	

	subroutine prims_to_cart(atom_num, n_prims, dq, q, x_1, x_2, Bmat_p, prim_list)
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
	integer(i4b) :: atom_num, n_prims, i, iter_counter, prim_list(n_prims, 4)
	real(sp) :: dq(n_prims), q(n_prims), q_2(n_prims), x_1(3 * atom_num), x_2(3 * atom_num)
	real(sp) :: BT_Ginv(n_prims, (3 * atom_num)), Bmat_p((3 * atom_num), n_prims)
	real(sp) :: dq_actual(n_prims), check(n_prims), temp_check, xyz_rms_1, xyz_rms_2
	real(sp) :: init_dq(n_prims), target_q(n_prims), dx(3 * atom_num), Gmat(n_prims, n_prims)
	real(sp), allocatable :: temp_q(:)
	logical :: convergence
	
    ! Since cartesians are rectilinear and internal coordinates are curvilinear, a simple transformation cannot be used.
    ! Instead, an iterative transformation procedure must be used.
    ! The expression B(transpose) * G(inverse) is initialised as it is used to convert between coordinate systems.
	Gmat = MATMUL(TRANSPOSE(Bmat_p), Bmat_p)
	BT_Ginv = MATMUL(SVD_INVERSE(Gmat, n_prims, n_prims), TRANSPOSE(Bmat_p))

	! Stashing some values for convergence criteria.
	convergence = .FALSE.
	xyz_rms_1 = 0
	xyz_rms_2 = 0
	init_dq(:) = dq(:)
	target_q(:) = q(:) + init_dq(:)

	do while (convergence .eqv. .FALSE.) 
		! The change in cartesian coordinates associated with the change in primitive internal coordinates is calculated.
		dx(:) = 0.0
		dx = MATMUL(TRANSPOSE(BT_Ginv), dq)
		
		! The root-mean-square change is used as a convergence criteria, so it is evaluated.
		xyz_rms_2 = RMSD_CALC(dx, x_1, (atom_num * 3))
		
		! The new cartesian geometry is evaluated, and the new primitive internal coordinate set and Wilson B matrix are obtained.
		x_2 = x_1 + dx
		q(:) = 0.0
		Bmat_p(:,:) = 0.0
		if (ALLOCATED(temp_q)) then 
			deallocate(temp_q)
		end if
		call calc_prims(atom_num, n_prims, temp_q, x_2, prim_list) 
		!call gen_Bmat_prims(atom_num, n_prims, x_2, prim_list, Bmat_p)
		q(:) = temp_q(:)
		
		! The Moore-Penrose inverse is constructed for the next iteration.
		!BT_Ginv(:,:) = 0.0
		!Gmat(:,:) = 0.0
		!Gmat = MATMUL(TRANSPOSE(Bmat_p), Bmat_p)
		!BT_Ginv = MATMUL(SVD_INVERSE(Gmat, n_prims, n_prims), TRANSPOSE(Bmat_p))

        ! The change in primitive internals for the next iteration is evaluated, and any which do not change in the original change in primitive internals is set to zero.
		! This mitigates any risk of primitive internal coordinates changing which should remain constant, thus making the interpolation linear.
		dq(:) = 0.0
		dq = (target_q - q)
		do i=1, SIZE(dq)
			if (init_dq(i) .eq. 0.0) then
				dq(i) = 0.0
			end if
		end do
		
		! Now, the three exit conditions should be checked...
		! The first ending condition for this transformation is when the root-mean-square change in cartesians is less than 10^-6.
        ! The second ending condition for this transformation is when the difference in root-mean-square change in cartesians between iteration i and i+1 is less than 10^-12.
        ! The third ending condition for this transformation is when the difference between the target DIC and the calculated DLC is less than 10^-6.
		check = target_q - q
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
		xyz_rms_1 = xyz_rms_2

		! In some cases, cartesians cannot be solved, so an exit condition must exist for this case.
		iter_counter = iter_counter + 1
		if (iter_counter == 1000) then
			print *, "Error; could not solve cartesians from the change in primitive internal coordinates."
			exit
		end if
	end do

	end subroutine prims_to_cart
	

END MODULE primitive
