MODULE math
use nrtype ; use coordinates ; use optimdata
implicit none

! In this module, a collection of maths functions are found.
! These functions are useful in general, but specifically designed for the generation of delocalised internal coordinates.

contains
	
	
	function ATOM_DISTANCE(coords_i, coords_j) result(r)
	! Here, two sets of cartesian coordinates are taken, and the distance between the points is calculated.
	! This is used to obtain an internal coordinate set.
	!
	! ARGUMENTS:	coords_i : 1D array containing the x,y,z coordinates of atom i.
	!				coords_j : 1D array containing the x,y,z coordinates of atom j.
	
	implicit none
	integer(i4b) :: i
	real(sp), intent(in) :: coords_i(3), coords_j(3)
	real(sp) :: r, v1(3)
	
	! Calculating vector between both coordinates.
	do i = 1, 3
		v1(i) = coords_i(i) - coords_j(i)
	end do
	
	! Calculating the distance between the two coordinates.
	r = 0.0
	r = NORM2(v1)

	end function ATOM_DISTANCE
	
	
	function ATOM_DISTANCE_GRAD(coords_i, coords_j) result(grad_r)
	! Here, two sets of cartesian coordinates are taken, and the first derivatives for the distance between the points are calculated.
	! This is used to fill the Wilson B matrix in generation of a primitive internal coordinate space.
	!
	! ARGUMENTS:	coords_i : 1D array containing the x,y,z coordinates of atom i.
	!				coords_j : 1D array containing the x,y,z coordinates of atom j.
	
	implicit none
	integer(i4b) :: i
	real(sp) :: unit_vec(3), v1(3)
	real(sp), intent(in) :: coords_i(3), coords_j(3)
	real(sp) :: grad_r(3,2)
	
	! Calculating vector between both coordinates.
	do i = 1, 3
		v1(i) = coords_i(i) - coords_j(i)
	end do
	
	! Calculating analytical first derivatives.
	unit_vec = unit_vector(v1, SIZE(v1))
	grad_r = 0.0
	grad_r(:,1) = unit_vec
	grad_r(:,2) = -1 * unit_vec

	end function ATOM_DISTANCE_GRAD
	
	
	function ATOM_ANGLE(coords_i, coords_j, coords_k) result(theta)
	! Here, three sets of cartesian coordinates are taken, and the angle between the points is calculated.
	! This is used to obtain an internal coordinate set.
	!
	! ARGUMENTS:	coords_i : 1D array containing the x,y,z coordinates of atom i.
	!				coords_j : 1D array containing the x,y,z coordinates of atom j.	
	!				coords_k : 1D array containing the x,y,z coordinates of atom k.	
	
	implicit none
	integer(i4b) :: i
	real(sp), intent(in) :: coords_i(3), coords_j(3), coords_k(3)
	real(sp) :: theta, v1(3), v2(3), cos_val
	
	! First, evaluate the vectors between coordinates i and j, and k and j.
	do i = 1, 3
		v1(i) = coords_i(i) - coords_j(i)
		v2(i) = coords_k(i) - coords_j(i)
	end do
	
	! Now, the angle theta can be calculated by the usual scalar product cosine relationship.
	! It is possible that the value within the cosine evaluation is either <-1 or >1.
	! This is unlikely, but in the event that is the case then theta is set manually.
	cos_val = INNER_PRODUCT(v1,v2,SIZE(v1)) / (NORM2(v1) * NORM2(v2))
	theta = 0.0
	if (cos_val .le. -1) then
		theta = PI
	else if (cos_val .ge. 1) then
		theta = 0.0
	else
		theta = ACOS(cos_val)
	end if
	
	end function ATOM_ANGLE
	
	
	function ATOM_ANGLE_GRAD(coords_i, coords_j, coords_k) result(grad_theta)
	! Here, three sets of cartesian coordinates are taken, and the first derivatives for the angle between the points is calculated.
	! This is used to obtain an internal coordinate set.
	!
	! ARGUMENTS:	coords_i : 1D array containing the x,y,z coordinates of atom i.
	!				coords_j : 1D array containing the x,y,z coordinates of atom j.	
	!				coords_k : 1D array containing the x,y,z coordinates of atom k.	
	
	implicit none
	integer(i4b) :: i
	real(sp), intent(in) :: coords_i(3), coords_j(3), coords_k(3)
	real(sp) :: grad_theta(3,3), theta, v1(3), v2(3), cos_val	
	
	! First, evaluate the vectors between coordinates i and j, and k and j.
	v1 = 0.0
	v2 = 0.0
	do i = 1, 3
		v1(i) = coords_i(i) - coords_j(i)
		v2(i) = coords_k(i) - coords_j(i)
	end do
	
	! Now, the angle theta can be calculated by the usual scalar product cosine relationship.
	! It is possible that the value within the cosine evaluation is either <-1 or >1.
	! This is unlikely, but in the event that is the case then theta is set manually.
	cos_val = 0.0
	cos_val = INNER_PRODUCT(v1,v2,SIZE(v1)) / (NORM2(v1) * NORM2(v2))
	theta = 0.0
	if (cos_val .le. -1) then
		theta = PI
	else if (cos_val .ge. 1) then
		theta = 0.0
	else
		theta = ACOS(cos_val)
	end if

	! If the angle is very close to pi, then pi is used to calculate the gradient so that there are no problems with linear molecules.
	! Otherwise, the typical analytical first derivatives are used.
	! These terms are shamelessly ripped from PyBerny...
	grad_theta = 0.0
	if (ABS(theta) .gt. (PI - 1E-6)) then
		grad_theta(:,1) = (PI - theta) / (2 * NORM2(v1) ** 2) * v1
		grad_theta(:,2) = (1 / NORM2(v1) - 1 / NORM2(v2)) * (PI - theta) / (2 * NORM2(v1)) * v1
		grad_theta(:,3) = (PI - theta) / (2 * NORM2(v2) ** 2) * v2
	else
		grad_theta(:,1) = 1 / TAN(theta) * v1 / NORM2(v1) ** 2 - v2 / (NORM2(v1) * NORM2(v2) * SIN(theta))
	    grad_theta(:,2) = (v1 + v2) / (NORM2(v1) * NORM2(v2) * SIN(theta)) &
                           - 1 / TAN(theta) * (v1 / norm2(v1) ** 2 + v2 / norm2(v2) ** 2)
	    grad_theta(:,3) = 1 / TAN(theta) * v2 / NORM2(v2) ** 2 &
 		                   - v1 / (NORM2(v1) * NORM2(v2) * SIN(theta))
	end if
	
	end function ATOM_ANGLE_GRAD
	
	
	function ATOM_DIHEDRAL(coords_i, coords_j, coords_k, coords_l) result(phi)
	! Here, three sets of cartesian coordinates are taken, and the dihedral torsion between the points is calculated.
	! This is used to obtain an internal coordinate set.
	!
	! ARGUMENTS:	coords_i : 1D array containing the x,y,z coordinates of atom i.
	!				coords_j : 1D array containing the x,y,z coordinates of atom j.	
	!				coords_k : 1D array containing the x,y,z coordinates of atom k.	
	!				coords_l : 1D array containing the x,y,z coordinates of atom l.	
	
	implicit none
	integer(i4b) :: i
	real(sp), intent(in) :: coords_i(3), coords_j(3), coords_k(3), coords_l(3)
	real(sp) :: phi, v1(3), v2(3), w(3), ew(3), a1(3), a2(3)
	real(sp) :: det_array(3,3), sgn, deter, dot_prod

	! First, evaluate the vectors between coordinates i and j, l and k, and k and j.
	v1 = 0.0
	v2 = 0.0
	w = 0.0
	do i = 1, 3
		v1(i) = coords_i(i) - coords_j(i)
		v2(i) = coords_l(i) - coords_k(i)
		w(i)  = coords_k(i) - coords_j(i) 
	end do
	
	! The calculation of the dihedral involves the normal vectors to the planes containing v1 and w, and v2 and w, respectively.
	ew = 0.0
	a1 = 0.0
	a2 = 0.0
    ew = w / NORM2(w)
    a1 = v1 - INNER_PRODUCT(v1,ew,SIZE(v1)) * ew
    a2 = v2 - INNER_PRODUCT(v2,ew,SIZE(v1)) * ew	
	
	! Now, the determinant of the normal vectors is calculated as its sign is used to evaluate the sign of the dihedral angle.
	det_array = 0.0
	det_array(1,:) = v2
	det_array(2,:) = v1
	det_array(3,:) = w
	deter = 0.0
	deter = DETERMINANT(det_array, SIZE(det_array,1))
	if (deter .le. 0.0) then 
		sgn = 1
	else if (deter .gt. 0.0) then 
		sgn = -1
	end if
	
	! To ensure that the dihedral angle stays in the range -PI < 0 < PI, some conditions are enforced on its value.
	dot_prod = 0.0
	dot_prod = INNER_PRODUCT(a1,a2,SIZE(a1)) / (NORM2(a1) * NORM2(a2))
	if (dot_prod .lt. -1.0) then
		dot_prod = -1.0
	else if (dot_prod .gt. 1.0) then
		dot_prod = 1.0
	end if
	
	! Now, the dihedral angle can be calculated.	
	phi = 0.0
	phi = ACOS(dot_prod) * sgn
	
	end function ATOM_DIHEDRAL
	
	
	function ATOM_DIHEDRAL_GRAD(coords_i, coords_j, coords_k, coords_l) result(grad_phi)
	! Here, four sets of cartesian coordinates are taken, and the first derivatives for the dihedral torsion between the points is calculated.
	! This is used to obtain an internal coordinate set.
	!
	! ARGUMENTS:	coords_i : 1D array containing the x,y,z coordinates of atom i.
	!				coords_j : 1D array containing the x,y,z coordinates of atom j.	
	!				coords_k : 1D array containing the x,y,z coordinates of atom k.	
	!				coords_l : 1D array containing the x,y,z coordinates of atom l.	
	
	implicit none
	integer(i4b) :: i
	real(sp), intent(in) :: coords_i(3), coords_j(3), coords_k(3), coords_l(3)
	real(sp) :: phi, v1(3), v2(3), w(3), ew(3), a1(3), a2(3), g(3), A, B, grad_phi(3,4)
	real(sp) :: det_array(3,3), sgn, deter, dot_prod
	
	! First, evaluate the vectors between coordinates i and j, l and k, and k and j.
	v1 = 0.0
	v2 = 0.0
	w = 0.0
	do i = 1, 3
		v1(i) = coords_i(i) - coords_j(i)
		v2(i) = coords_l(i) - coords_k(i)
		w(i)  = coords_k(i) - coords_j(i) 
	end do
	
	! The calculation of the dihedral involves the normal vectors to the planes containing v1 and w, and v2 and w, respectively.
	ew = 0.0
	a1 = 0.0
	a2 = 0.0
    ew = w / NORM2(w)
    a1 = v1 - INNER_PRODUCT(v1,ew,SIZE(v1)) * ew
    a2 = v2 - INNER_PRODUCT(v2,ew,SIZE(v1)) * ew	
	
	! Now, the determinant of the normal vectors is calculated as its sign is used to evaluate the sign of the dihedral angle.
	det_array = 0.0
	det_array(1,:) = v2
	det_array(2,:) = v1
	det_array(3,:) = w
	deter = 0.0
	deter = DETERMINANT(det_array, SIZE(det_array,1))
	if (deter .le. 0.0) then 
		sgn = 1
	else if (deter .gt. 0.0) then 
		sgn = -1
	end if
	
	! To ensure that the dihedral angle stays in the range -PI < 0 < PI, some conditions are enforced on its value.
	dot_prod = 0.0
	dot_prod = INNER_PRODUCT(a1,a2,SIZE(a1)) / (NORM2(a1) * NORM2(a2))
	if (dot_prod .lt. -1.0) then
		dot_prod = -1.0
	else if (dot_prod .gt. 1.0) then
		dot_prod = 1.0
	end if
	
	! Now, the dihedral angle can be calculated.	
	phi = 0.0
	phi = ACOS(dot_prod) * sgn

	! If the angle is very close to zero or pi, then alternative terms are used to calculate the gradient so that there are no problems with linear molecules.
	! Otherwise, the typical analytical first derivatives are used.
	! These terms are shamelessly ripped from PyBerny...
    g = CROSS_PRODUCT(w,a1)
    g = g / NORM2(g)
    A = INNER_PRODUCT(v1,ew,SIZE(v1)) / NORM2(w)
    B = INNER_PRODUCT(v2,ew,SIZE(v1)) / NORM2(w)
	if (ABS(phi) .gt. (PI - 1E-6)) then
        grad_phi(:,1) = g / (NORM2(g) * NORM2(a1))
        grad_phi(:,2) = -((1 - A) / NORM2(a1) - B / NORM2(a2)) * g
        grad_phi(:,3) = -((1 + B) / NORM2(a2) + A / NORM2(a1)) * g
        grad_phi(:,4) = g / (NORM2(g) * NORM2(a2))
    else if (ABS(phi) .lt. 1E-6) then
        grad_phi(:,1) = g / (NORM2(g) * NORM2(a1))
        grad_phi(:,2) = -((1 - A) / NORM2(a1) + B / NORM2(a2)) * g
        grad_phi(:,3) = ((1 + B) / NORM2(a2) - A / NORM2(a1)) * g
        grad_phi(:,4) = -g / (NORM2(g) * NORM2(a2))
    else
        grad_phi(:,1) = 1 / TAN(phi) * a1 / NORM2(a1) ** 2 - a2 / (NORM2(a1) * NORM2(a2) * SIN(phi))
        grad_phi(:,2) = ((1 - A) * a2 - B * a1) / (NORM2(a1) * NORM2(a2) * SIN(phi)) &
						- 1 / TAN(phi) * ((1 - A) * a1 / NORM2(a1) ** 2 - B * a2 / NORM2(a2) ** 2)
        grad_phi(:,3) = ((1 + B) * a1 + A * a2) / (NORM2(a1) * NORM2(a2) * SIN(phi)) &
						- 1 / TAN(phi) * ((1 + B) * a2 / NORM2(a2) ** 2 + A * a1 / NORM2(a1) ** 2)	
        grad_phi(:,4) = 1 / TAN(phi) * a2 / NORM2(a2) ** 2 - a1 / (NORM2(a1) * NORM2(a2) * SIN(phi))
	end if
	
	end function ATOM_DIHEDRAL_GRAD	
	
	
	function VECTOR_PROJECT(vec1, vec2, length) result(proj_vec)
	! Here, the first vector taken is projected onto the second vector taken.
	! In this case, the vectors must of necessity be of the same dimensions.
	!
	! ARGUMENTS:	vec1   : 1D array containing the vector which is to be projected.
	!				vec2   : 1D array containing the vector which vec1 is projected upon.
	!               length : integer which represents the number of elements in both vectors.
		
	implicit none
	integer(i4b), intent(in) :: length
	real(sp), intent(in) :: vec1(length), vec2(length)
	real(sp) :: proj_vec(length), unit_vec(length)
	
	! The projection space is normalised.
    unit_vec = UNIT_VECTOR(vec2, SIZE(vec2))
    
    ! The vector is projected by the usual formula.
    proj_vec = DOT_PRODUCT(vec1, unit_vec) * unit_vec

	end function VECTOR_PROJECT
	
	
	function UNIT_VECTOR(vec, length) result(unit_vec)
	! Here, a vector is taken and is transformed to a unit vector.
	!
	! ARGUMENTS:	vec    : 1D array containing the vector which is made to a unit vector.
	!               length : integer which represents the number of elements in the vector.	
	
	implicit none
	integer(i4b), intent(in) :: length
	real(sp), intent(in) :: vec(length)
	real(sp) :: unit_vec(length)

	! The unit vector is calculated by the usual formula.
	unit_vec = vec / NORM2(vec)

	end function UNIT_VECTOR
	
	
	function INNER_PRODUCT(vec1, vec2, length) result(inner)
	! Here, the inner product of two 1D vectors is calculated.
	!
	! ARGUMENTS:	vec1   : 1D array containing vector 1.
	!				vec2   : 1D array containing vector 2.
	!               length : integer which represents the number of elements in both vectors.
	
	implicit none
	integer(i4b), intent(in) :: length
	integer(i4b) :: i
	real(sp) :: temp, inner
	real(sp), intent(in) :: vec1(length), vec2(length)
	
	! The inner product is calculated simply as:
	! [a1, a2, a3] * [b1, b2, b3] = SUM[a1*b1, a2*b2, a3*b3]
	inner = 0.0
	temp = 0.0
	do i=1, length
		temp = vec1(i) * vec2(i)
		inner = inner + temp
	end do
	
	end function INNER_PRODUCT
	
	
	function OUTER_PRODUCT(vec1, vec2, n1, n2) result(outer)
	! Here, the outer product of two 1D vectors is calculated.
	!
	! ARGUMENTS:	vec1   : 1D array containing vector 1.
	!				vec2   : 1D array containing vector 2.
	!               n1     : integer which represents the number of elements in vector 1.
	!               n2     : integer which represents the number of elements in vector 2.	
	
	implicit none
	integer(i4b), intent(in) :: n1, n2
	integer(i4b) :: i, j
	real(sp) :: outer(n1,n2)
	real(sp), intent(in) :: vec1(n1), vec2(n2)
	
	! The outer product is simply calculated by the usual formula.
	outer = 0.0
	do i=1, n1
		do j=1, n2
			outer(i,j) = vec1(i) * vec2(j)
		end do
	end do
	outer = TRANSPOSE(outer)
	
	end function OUTER_PRODUCT
	
	
	function CROSS_PRODUCT(vec1, vec2) result(cross)
	! Here, the cross product of two 1D vectors of length 3 is calculated.
	!
	! ARGUMENTS:	vec1   : 1D array containing vector 1.
	!				vec2   : 1D array containing vector 2.
	
	implicit none
	real(sp) :: cross(3)
	real(sp), intent(in) :: vec1(3), vec2(3)
	
	! The cross product is simply calculated by the usual forumla.
	cross = 0.0
	cross(1) = (vec1(2) * vec2(3)) - (vec1(3) * vec2(2))
	cross(2) = (vec1(3) * vec2(1)) - (vec1(1) * vec2(3))
	cross(3) = (vec1(1) * vec2(2)) - (vec1(2) * vec2(1))
	
	end function CROSS_PRODUCT
		
	
	function SVD_INVERSE(A, rows, cols) result(A_inv)
	! Here, a square matrix is taken as input and the generalised inverse is calculated by using singlular value decomposition (SVD).
	!
	! ARGUMENTS:	A    : 2D array containing the array which is to be inverted by SVD.
	!               rows : integer which represents the number of rows.
    !				cols : integer which represents the number of columns.
	
	implicit none
	integer(i4b), intent(in) :: rows, cols
	integer(i4b) :: LDA, LDU, LWORK, LDVT, INFO, i, j
	real(sp), intent(in) :: A(rows, cols)
	real(sp), allocatable :: WORK(:)
	real(sp) :: U(rows, rows), VT(cols, cols), S(cols), S_inv(rows, cols), A_inv(rows, cols)

	! Initialising the leading dimensions for the SVD arrays.
	LDA = rows
	LDU = rows
	LDVT = cols
	
	! Zero the arrays used in this function.
	S = 0.0
	VT = 0.0
	U = 0.0
	S_inv = 0.0
	A_inv = 0.0

	! A working array is created for use in calculation of SVD.
	LWORK = MAX(1, (3 * MIN(rows, cols) + MAX(rows, cols)), (5 * min(rows, cols)))
	allocate(WORK(LWORK))
	
	! SVD defined in LAPACK is performed to obtain matrices S, U and VT.
	call DGESVD('A', 'A', rows, cols, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO)

	! To calculate the inverse of matrix A, the vector S must first be inverted.
	do i=1, SIZE(S, 1)
		if (INT(S(i)) .ne. 0) then
			S(i) = 1 / S(i)
		end if
	end do

	! Next, the vector S must be converted to a diagonal matrix.
	do j=1, SIZE(S_inv, 1)
		S_inv(j,j) = S(j)
	end do
	
	! Lastly, the inverse is calculated by the usual formula.
	A_inv = MATMUL(MATMUL(TRANSPOSE(VT), S_inv), TRANSPOSE(U))

	deallocate(WORK)
	
	end function SVD_INVERSE

	
	function IS_ORTHOG(vecs, length, n) result(orthogonality)
	! Here, an array of vectors is taken, and their orthogonality is verified.
	! This is useful as a check after an orthogonalisation procedure.
	!
	! ARGUMENTS:	vecs   : 2D array containing the set of vectors whose orthogonality is checked.
	!               length : integer which represents the length of the vectors.
    !				n      : integer which represents the number of vectors in the set.

	implicit none
	integer(i4b) :: i
	integer(i4b), intent(in) :: length, n
	real(sp), intent(in) :: vecs(length, n)
	real(sp) :: temp_dot
	logical :: orthogonality
	
	! Takes the dot product of each vector pair, and if the result is close to zero, then the two vectors are orthogonal.
	do i = 1, (n - 1)
		temp_dot = 0.0
		
		! The dot product always has some numerical precision remainder, so a margin of 1E-5 is used.
		! Therefore, if the dot product between the two vectors is greater than this value, then the set is not orthogonal.
		temp_dot = DOT_PRODUCT(vecs(:,i), vecs(:, i+1))
		if (ABS(temp_dot) > (1E-05)) then
			orthogonality = .FALSE.
			return
		end if
	end do
	
	! The set of vectors is found to be orthogonal if the exit condition in the above loop is not satisfied.
	orthogonality = .TRUE.
	
	end function IS_ORTHOG
	
	
	function DETERMINANT(matrix, n) result(det)
	! Here, the determinant of a square matrix is calculated. 
	!
	! ARGUMENTS:	matrix : 2D array containing the matrix which the determinant of will be calculated.
	!               n      : integer which represents the number of rows/columns (it doesn't matter which as the matrix is square).
	
	implicit none
	integer(i4b), intent(in) :: n
	integer(i4b) :: i, j, k, l
    real(sp), intent(in) :: matrix(n,n)
    real(sp) :: m, temp, det, work_mat(n,n)
    logical :: DetExists
	
	! Copying the input matrix to the working matrix...
	work_mat(:,:) = matrix(:,:)
	
	! Initialising some values.
	DetExists = .TRUE.
    l = 1
	
    ! The matrix is converted to upper diagonal form.
    do k=1, (n - 1)
        if (INT(work_mat(k,k)) == 0) then
            DetExists = .FALSE.
            do i=k+1, n
                if (INT(work_mat(i,k)) .ne. 0) then
                    do j=1, n
                        temp = work_mat(i,j)
                        work_mat(i,j) = work_mat(k,j)
                        work_mat(k,j) = temp
                    end do
                    DetExists = .TRUE.
                    l = -l
                    exit
                end if
            end do
            if (DetExists .EQV. .FALSE.) then
                det = 0
                return
            end if
        end if
        do j=k+1, n
            m = work_mat(j,k) / work_mat(k,k)
            do i=k+1, n
                work_mat(j,i) = work_mat(j,i) - m*work_mat(k,i)
            end do
        end do
    end do
   
    ! The determinant is calculated by finding the product of diagonal elements.
    det = l
    do i = 1, n
        det = det * work_mat(i,i)
    end do
   
	end function DETERMINANT
	
	
	function EVALS(matrix, n) result(eigenvals)
	! Here, the eigenvalues from a square, symmetric, real matrix are calculated.
	! The eigenvectors associated with said eigenvalues are not given as a result, but can be obtained with the function EVEVS.
	!
	! ARGUMENTS:	matrix : 2D array containing the matrix which the eigenvalues of will be calculated.
	!               n      : integer which represents the number of rows/columns (it doesn't matter which as the matrix is square).
    	
	implicit none
	integer(i4b), intent(in) :: n
	integer(i4b) :: LDA, LWORK, INFO
	real(sp), intent(in) :: matrix(n, n)
	real(sp), allocatable :: WORK(:)
	real(sp) :: eigenvals(n), work_mat(n,n)
	character :: JOBZ, UPLO

	! Assigning working matrix...
	work_mat(:,:) = 0.0
	work_mat(:,:) = matrix(:,:)

	! Initialising some values....
	JOBZ = 'N'
	UPLO = 'U'
	LDA = n
	INFO = 0
	
	! Allocating working array...
	LWORK = MAX(1, (3 * n) - 1)
	allocate(WORK(LWORK))
	WORK(:) = 0.0
	
	! Obtaining eigenvalues...
	eigenvals(:) = 0.0
	call DSYEV(JOBZ, UPLO, n, work_mat, LDA, eigenvals, WORK, LWORK, INFO)
	
	deallocate(WORK)

	end function EVALS
	
	
	function EVECS(matrix, n) result(eigenvecs)
	! Here, the eigenvectors of a square, symmetric, real matrix are calculated.
	! The eigenvalues associated with said eigenvectors are not given as a result, but can be obtained with the function EVALS.
	!
	! ARGUMENTS:	matrix : 2D array containing the matrix which the eigenvalues of will be calculated.
	!               n      : integer which represents the number of rows/columns (it doesn't matter which as the matrix is square).
	   	
	implicit none
	integer(i4b), intent(in) :: n
	integer(i4b) :: LDA, LWORK, INFO
	real(sp), intent(in) :: matrix(n, n)
	real(sp), allocatable :: WORK(:)
	real(sp) :: eigenvecs(n, n), eigenvals(n), work_mat(n,n)
	character :: JOBZ, UPLO

	! Assigning working matrix...
	work_mat(:,:) = 0.0
	work_mat(:,:) = matrix(:,:)

	! Initialising some values....
	JOBZ = 'V'
	UPLO = 'U'
	LDA = n
	INFO = 0
	
	! Allocating working array....
	LWORK = MAX(1, (3 * n) - 1)
	allocate(WORK(LWORK))
	WORK(:) = 0.0

	! Obtaining eigenvectors...
	eigenvals(:) = 0.0
	call DSYEV(JOBZ, UPLO, n, work_mat, LDA, eigenvals, WORK, LWORK, INFO)
	eigenvecs(:,:) = 0.0
	eigenvecs(:,:) = work_mat(:,:)

	deallocate(WORK)
	
	end function EVECS

	
	function COMBINATIONS_2(integers, n) result(combos_2)
	! Here, a series of integers is taken, and outputs all possible 2-number combinations of these integers is generated.
	! Critically, the algorithm does not produce duplicates, so is useful in defining a total connectivity scheme.
	!
	! ARGUMENTS:      integers : 1D array of integers for which the possible combinations of are to be generated.
	!                 n        : integer which represents the length of the array integers.
	
	implicit none
	integer(i4b), intent(in) :: n, integers(n)
	integer(i4b), allocatable :: combos_2(:,:)
	integer(i4b) :: i, j, k, l, combo_alloc, combos_counter, integers_save(n), integers_work(n)
	
	! Saving the integers for future use...
	integers_save(:) = integers(:)
	integers_work(:) = integers(:)
	
	! Firstly, the output array, combos_2, must be allocated.
	combo_alloc = 0
	j = SIZE(integers_work) - 1
	do i=1, SIZE(integers_work)
		combo_alloc = combo_alloc + j
		j = j - 1
	end do
	allocate(combos_2(combo_alloc, 4))
	
	! The 2-integer combinations are generated.
	combos_counter = 1
	do i=1, SIZE(integers_save)
		k = integers_save(i)
		do j=1, SIZE(integers_work)
			l = integers_work(j)
			
			! If the integer is 0 (which is set below) or the integers are the same, then the cycle is skipped.
			! i.e., only [1,2] and not [2,1] is saved.
			if ((l == 0) .or. (l == k)) then
			    cycle
			end if
			
			! Populating the output array...
			combos_2(combos_counter,1) = k
			combos_2(combos_counter,2) = l
			combos_2(combos_counter,3) = 0
			combos_2(combos_counter,4) = 0
			combos_counter = combos_counter + 1
		end do
		integers_work(i) = 0
	end do

	end function COMBINATIONS_2
	
	
	function COMBINATIONS_2_DUPE(integers, n) result(combos_2)
	! Here, a series of integers is taken, and outputs all possible 2-number combinations of these integers is generated.
	! Unlike COMBINATIONS_2, the algorithm does produce duplicates, which is useful in generation defining a full primitive internal coordinate scheme.
	!
	! ARGUMENTS:      integers : 1D array of integers for which the possible combinations of are to be generated.
	!                 n        : integer which represents the length of the array integers.
	
	implicit none
	integer(i4b), intent(in) :: n, integers(n)
	integer(i4b), allocatable :: combos_2(:,:)
	integer(i4b) :: i, j, k, l, combo_alloc, combos_counter
	
	! Firstly, the output array, combos_2, must be allocated.
	combo_alloc = 0
	do i=1, SIZE(integers)
		combo_alloc = combo_alloc + (1 * (n-1))
	end do
	allocate(combos_2(combo_alloc, 2))
	
	! Now, generating all possible combinations, including reverse duplicates.
	! i.e., both [1,2] and [2,1] are saved.
	combos_counter = 1
	do i=1, n
		k = integers(i)
		do j=1, n
			l = integers(j)
			if (k .ne. l) then
				combos_2(combos_counter,1) = k
				combos_2(combos_counter,2) = l
				combos_counter = combos_counter + 1
			end if
		end do
	end do		
		
	end function COMBINATIONS_2_DUPE
	
	
	function ASCEND_ORDER(integers, n) result(integers_order)
	! Here, a series of integers is taken, and a 1D array of the same size is returned except it is in numerical order.
	!
	! ARGUMENTS:      integers : 1D array of integers which are to be put in ascending order.
	!                 n        : integer which represents the length of the array integers.
	
	implicit none
	integer(i4b), intent(in) :: n, integers(n)
	integer(i4b) :: i, j, temp_int, integers_order(n)

	! The numbers are made into ascending order through a simple do loop construct.
	integers_order(:) = integers(:)
	do i=1, n-1
		do j=i+1, n
			if (integers_order(i) .gt. integers_order(j)) then
				temp_int = integers_order(i)
				integers_order(i) = integers_order(j)
				integers_order(j) = temp_int
			end if
		end do
	end do

	
	end function ASCEND_ORDER
	
	
	function RMSD_CALC(dx, x, n) result(RMSD)
	! Here, a series of values is taken (typically coordinates, gradients, or forces) and the root-mean-square-deviation is calculated.
	!
	! ARGUMENTS:     dx : 1D array containing the change in x.
	!				 x  : 1D array containing the quantity for the change to be calculated against.
	!                n  : Length of both 1D arrays which must, by definition, be the same.
	
	implicit none
	integer(i4b), intent(in) :: n
	real(sp), intent(in) :: dx(n), x(n)
	real(sp) :: x_n(n), RMSD
	
	! The new x is calculated.
	! Indeed, this expression simply cancels out in the subsequent equations, but it is included for clarity.
	x_n(:) = x(:) + dx(:)
	
	! Now, the RMSD can be calculated.
	RMSD = 0.0
	RMSD = SQRT((SUM((x_n - x)**2)) / n)
	
	end function RMSD_CALC
	
	
END MODULE math