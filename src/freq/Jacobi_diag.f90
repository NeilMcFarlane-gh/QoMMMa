SUBROUTINE Diag(N,M1,M2,Vec)
 implicit none

! This subroutine diagonalizes the N * N input matrix, M1.
!	Required on Input:
!		N, the size of the (square) matrix
!		M1, the matrix
!	Output:
!		M1 is unchanged
!		M2 is the matrix of eigenvectors
!		Vec is a vector with the eigenvalues
!			These are sorted from smallest to largest

INTEGER, INTENT(IN) :: N
DOUBLE PRECISION, DIMENSION(N,N), INTENT(IN) :: M1
DOUBLE PRECISION, DIMENSION(N,N), INTENT(OUT) :: M2
DOUBLE PRECISION, DIMENSION(N), INTENT(OUT) :: Vec

INTEGER :: nrot
DOUBLE PRECISION :: M3(N,N), M4(N,N), eigs(N)

M3 = M1

CALL Jacobi(M3,N,eigs,M4,nrot)
CALL sort_eigenvalues(eigs,M4,N)

Vec = eigs
M2 = M4

END SUBROUTINE Diag



SUBROUTINE sort_eigenvalues(E,M,N)
implicit none

!	This subroutine sorts the eigenvalues & eigenvectors produced by a diagonalisation
!	On Input : N is the dimension, E the vector of eigenvalues, M the matrix with the
!		eigenvectors in the columns M(:,i)

INTEGER, INTENT(IN) :: N
DOUBLE PRECISION, INTENT(INOUT), DIMENSION(N) :: E
DOUBLE PRECISION, INTENT(INOUT), DIMENSION(N,N) :: M

INTEGER :: i, j, k
DOUBLE PRECISION :: p

DO i = 1, N - 1
	k = i
	p = E(i)
	DO j = i+1, n
		if (E(j) .ge. p) THEN
			k = j
			p = E(j)
		END IF
	END DO
	if (k .ne. i) THEN
		E(k) = E(i)
		E(i) = p
		DO J = 1, n
			p = M(j,i)
			M(j,i) = M(j,k)
			M(j,k) = p
		END DO
	end if
END DO

END SUBROUTINE sort_eigenvalues

