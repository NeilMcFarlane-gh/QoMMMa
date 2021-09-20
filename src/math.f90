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
	real(sp) :: r
	
	! Calculating the distance between the two coordinates.
	r = 0.0
	do i = 1, 3
		r = r + (coords_i(i) - coords_j(i))**2
	end do
	r = SQRT(r)
	
	end function ATOM_DISTANCE
	
	
	function ATOM_DISTANCE_GRAD(coords_i, coords_j) result(grad)
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
	
	! Calculating analytical first derivatives.
	unit_vec = unit_vector(vec, SIZE(vec))
	grad(1,:) = unit_vec
	grad(2,:) = -1 * unit_vec
	
	end function ATOM_DISTANCE_GRAD
	
	
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
    unit_vec = unit_vector(vec2, SIZE(vec2))
    
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
	
	
	!function GRAM_SCHMIDT(vecs) result(res)
	! Here, an array of vectors is orthonormalised by the Gram Schmidt methodology.
	! The first vector in the array is taken as the first, and thus it will not change, and the last vector will drop out.
	
	! TO-DO : Will write later when it comes to programming the constrained optimisation process.
	
	!end function GRAM_SCHMIDT
	
	
	function SVD_INVERSE(A, rows, cols) result(A_inv)
	! Here, a matrix is taken as input and the generalised inverse is calculated by using singlular value decomposition (SVD).
	!
	! ARGUMENTS:	A    : 2D array containing the array which is to be inverted by SVD.
	!               rows : integer which represents the number of rows.
    !				cols : integer which represents the number of columns.
	
	implicit none
	integer(i4b), intent(in) :: rows, cols
	integer(i4b) :: LDA, LDU, LWORK, LDVT, INFO, i, j
	real(sp), intent(in) :: A(rows, cols)
	real(sp), allocatable :: WORK(:), A_inv(:,:)
	real(sp) :: U(rows, rows), VT(cols, cols), S(cols), S_inv(rows, cols)

	! Initialising the leading dimensions for the SVD arrays.
	LDA = rows
	LDU = rows
	LDVT = cols
	
	! A working array is created for use in calculation of SVD.
	LWORK = MAX(1, (3 * MIN(rows, cols) + MAX(rows, cols)), (5 * min(rows, cols)))
	allocate(WORK(LWORK))
	
	! SVD defined in LAPACK is performed to obtain matrices S, U and VT.
	call DGESVD('A', 'A', rows, cols, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO)
	
	! To calculate the inverse of matrix A, the vector S must first be inverted.
	do i=1, SIZE(S, 1)
		if (S(i) .ne. 0) then
			S(i) = 1 / S(i)
		end if
	end do

	! Next, the vector S must be converted to a diagonal matrix.
	S_inv(:,:) = 0.0
	do j=1, SIZE(S_inv, 1)
		S_inv(j,j) = S(j)
	end do
	
	! Lastly, the inverse is calculated by the usual formula.
	! NOTE: This only works for square matrices, which is likely all that will be necessary.
	!       If non-square matrices are to be used, the code will have to be updated.
	A_inv = MATMUL(MATMUL(TRANSPOSE(VT), S_inv), TRANSPOSE(U))
	
	end function SVD_INVERSE
	
	
	function IS_ORTHOG(vecs, length, n) result(orthogonality)
	! Here, an array of vectors is taken, and their orthogonality is verified.
	! This is useful as a check after an orthogonalisation procedure.
	!
	! ARGUMENTS:	vecs   : 2D array containing the set of vectors whose orthogonality is checked.
	!               length : integer which represents the length of the vectors.
    !				n      : integer which represents the number of vectors in the set.
	
	! NOTE: The integers length and n may be the wrong way around - needs to be checked when it comes to utilisation...
	
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
    real(sp) :: m, temp, det
    logical :: DetExists
	
	! Initialising some values.
	DetExists = .TRUE.
    l = 1
	
    ! The matrix is converted to upper diagonal form.
    do k=1, (n - 1)
        if (matrix(k,k) == 0) then
            DetExists = .FALSE.
            do i=k+1, n
                if (matrix(i,k) .ne. 0) then
                    do j=1, n
                        temp = matrix(i,j)
                        matrix(i,j) = matrix(k,j)
                        matrix(k,j) = temp
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
            m = matrix(j,k) / matrix(k,k)
            do i=k+1, n
                matrix(j,i) = matrix(j,i) - m*matrix(k,i)
            end do
        end do
    end do
   
    ! The determinant is calculated by finding the product of diagonal elements.
    det = l
    do i = 1, n
        det = det * matrix(i,i)
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
	real(sp) :: eigenvals(n)
	character :: JOBZ, UPLO
	
	! Initialising some values....
	JOBZ = 'N'
	UPLO = 'U'
	LDA = n
	
	! Allocating working array...
	LWORK = MAX(1, (3 * n) -1)
	allocate(WORK(LWORK))

	! Obtaining eigenvalues...
	call DSYEV(JOBZ, UPLO, n, matrix, LDA, eigenvals, WORK, LWORK, INFO)

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
	real(sp) :: eigenvecs(n, n), eigenvals(n)
	character :: JOBZ, UPLO

	
	! Initialising some values....
	JOBZ = 'V'
	UPLO = 'U'
	LDA = n
	
	! Allocating working array....
	LWORK = MAX(1, (3 * n) -1)
	allocate(WORK(LWORK))

	! Obtaining eigenvectors...
	call DSYEV(JOBZ, UPLO, n, matrix, LDA, eigenvals, WORK, LWORK, INFO)
	eigenvecs = matrix
	
	end function EVECS
	
	
END MODULE math