MODULE math
use nrtype ; use coordinates ; use optimdata
implicit none

! Maybe have some declarations here...

contains
	
	
	function atom_distance(coords_i, coords_j) result(r)
	! Here, two sets of cartesian coordinates are taken, and the distance between the points is calculated.
	!
	! ARGUMENTS:	coords_i : 1D array containing the x,y,z coordinates of atom i.
	!				coords_j : 1D array containing the x,y,z coordinates of atom j.
	
	implicit none
	integer(i4b) :: i
	real(sp), intent(in) :: coords_i(3), coords_j(3)
	real(sp) :: r
	
	! Calculating the distance between the two coordinates.
	r = 0.0
	do i = 1, 3
		r = r + (coords_i(i) - coords_j(i))**2
	end do
	r = SQRT(r)
	
	end function atom_distance
	
	
	function atom_distance_grad(coords_i, coords_j) result(grad)
	! Here, two sets of cartesian coordinates are taken, and the first derivatives for the distance between the points are calculated.
	! This is used to fill the Wilson B matrix in generation of a primitive internal coordinate space.
	!
	! ARGUMENTS:	coords_i : 1D array containing the x,y,z coordinates of atom i.
	!				coords_j : 1D array containing the x,y,z coordinates of atom j.
	
	implicit none
	integer(i4b) :: i
	real(sp) :: unit_vec(3), vec(3)
	real(sp), intent(in) :: coords_i(3), coords_j(3)
	real(sp) :: grad(2,3)
	
	! Calculating vector between both coordinates.
	do i = 1, 3
		vec(i) = coords_i(i) - coords_j(i)
	end do
	
	! Calculating analytical for first derivatives.
	unit_vec = unit_vector(vec, SIZE(vec))
	grad(1,:) = unit_vec
	grad(2,:) = -1 * unit_vec
	
	end function atom_distance_grad
	
	
	function vector_project(vec1, vec2, cols) result(proj_vec)
	! Here, the first vector taken is projected onto the second vector taken.
	! In this case, the vectors must be of the same dimensions.
	!
	! ARGUMENTS:	vec1 : 1D array containing the vector which is to be projected.
	!				vec2 : 1D array containing the vector which vec1 is projected upon.
	!               cols : integer which represents the number of elements in both vectors (which must of necessity be the same).
		
	implicit none
	integer(i4b), intent(in) :: cols
	real(sp) :: proj_vec(cols), unit_vec(cols)
	real(sp), intent(in) :: vec1(cols), vec2(cols)
	
	! The projection space is normalised.
    unit_vec = unit_vector(vec2, SIZE(vec2))
    
    ! The vector is projected by the usual formula.
    proj_vec = DOT_PRODUCT(vec1, unit_vec) * unit_vec

	end function vector_project
	
	
	function unit_vector(vec, cols) result(unit_vec)
	! Here, a vector is taken and is transformed to a unit vector.
	!
	! ARGUMENTS:	vec  : 1D array containing the vector which is made to a unit vector.
	!               cols : integer which represents the number of elements in the vector.	
	
	implicit none
	integer(i4b), intent(in) :: cols
	real(sp), intent(in) :: vec(cols)
	real(sp) :: unit_vec(cols)

	! The unit vector is obtained by the usual formula.
	unit_vec = vec / NORM2(vec)

	end function unit_vector
	
	
	!function Gram_Schmidt(vecs) result(res)
	! Here, an array of vectors is orthonormalised by the Gram Schmidt methodology.
	! The first vector in the array is taken as the first, and thus it will not change, and the last vector will drop out.
	
	! TO-DO : Will write later when it comes to constrained optimisation process.
	
	!end function Gram_Schmidt
	
	
	function SVD_inverse(A, cols, rows) result(inv_arr)
	! Here, a matrix is taken as input and the generalised inverse is calculated using single value decomposition.
	!
	! ARGUMENTS:	A    : 2D array containing the array which is to be inverted by single value decomposition.
	!               cols : integer which represents the number of columns.
    !				rows : integer which represents the number of rows.
	
	implicit none
	integer(i4b), intent(in) :: cols, rows
	integer(i4b) :: LDA, LDU, LWORK, LDVT, INFO
	character :: JOBU, JOBVT
	real(sp), intent(in) :: A(cols, rows)
	real(sp) :: U(cols, cols), V(rows, rows), S(rows), inv_arr(rows, cols)
	real(sp), allocatable :: WORK(:)
	
	JOBU = 'A'
	JOBVT = 'A'
	LDA = cols
	LDU = cols
	LDVT = rows
	
	LWORK = MAX(1, (3 * MIN(cols, rows) + MAX(cols, rows)), (5 * min(cols, rows)))
	
	allocate(WORK(LWORK))
	
	call DGESVD(JOBU, JOBVT, rows, cols, A, LDA, S, U, LDU, V, LDVT, WORK, LWORK, INFO)
	
	end function SVD_inverse
	
	
	function is_orthog(vecs, length, height) result(orthogonality)
	! Here, an array of vectors is taken as input and their orthogonality is verified.
	! This is useful as a check after an orthogonalisation procedure.
	!
	! ARGUMENTS:	vecs : 2D array containing the set of vectors whose orthogonality is checked.
	!               length : integer which represents the length of the vectors.
    !				height : integer which represents the number of vectors in the set.
	
	implicit none
	integer(i4b) :: i
	integer(i4b), intent(in) :: length, height
	logical :: orthogonality
	real(sp), intent(in) :: vecs(length, height)
	real(sp) :: temp_dot
	
	! Takes the dot product of each vector pair, and if the result is close to zero, then the two vectors are orthogonal.
	do i = 1, (height - 1)
		temp_dot = 0.0
		
		! The dot product always has some numerical precision remainder, so a margin of 1E-5 is used.
		! Therefore, if the dot product between the two vectors is greater than this value, then the set is not orthogonal.
		temp_dot = DOT_PRODUCT(vecs(:,i), vecs(:, i+1))
		if (ABS(temp_dot) > (1E-05)) then
			orthogonality = .FALSE.
			RETURN
		end if
	end do
	
	! The set of vectors is found to be orthogonal if the exit condition in the above loop is not satisfied.
	orthogonality = .TRUE.
	
	end function is_orthog
	
	
END MODULE math