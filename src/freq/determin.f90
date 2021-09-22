SUBROUTINE determin(N,A,DD)

INTEGER, INTENT(IN) :: N
DOUBLE PRECISION, INTENT(IN) :: A(N,N)
REAL, INTENT(OUT) :: DD
integer ID, INDX(N)
call LUDCMP(A,N,INDX,ID) 
DD=DFLOAT(ID)
do J=1, N
  DD=DD*A(J,J)
end do
write(6,*) 'determinant is ',DD
END SUBROUTINE Determin


Subroutine LUDCMP(A,N,INDX,D)
!*******************************************************
!*    LU decomposition routines used by test_lu.f90    *
!*                                                     *
!*                 F90 version by J-P Moreau, Paris    *
!* --------------------------------------------------- *
!* Reference:                                          *
!*                                                     *
!* "Numerical Recipes by W.H. Press, B. P. Flannery,   *
!*  S.A. Teukolsky and W.T. Vetterling, Cambridge      *
!*  University Press, 1986".                           *
!*                                                     * 
!*******************************************************
!  ***************************************************************
!  * Given an N x N matrix A, this routine replaces it by the LU *
!  * decomposition of a rowwise permutation of itself. A and N   *
!  * are input. INDX is an output vector which records the row   *
!  * permutation effected by the partial pivoting; D is output   *
!  * as -1 or 1, depending on whether the number of row inter-   *
!  * changes was even or odd, respectively. 			 *
!  ***************************************************************


 PARAMETER(NMAX=100,TINY=1.0D-20)
 REAL*8  AMAX,DUM, SUM, A(N,N),VV(NMAX)
 INTEGER D, INDX(N)
 D=1; 
 DO I=1,N
   AMAX=0.d0
   DO J=1,N
     IF (DABS(A(I,J)).GT.AMAX) AMAX=DABS(A(I,J))
   END DO ! j loop
   IF(AMAX.LT.TINY) THEN
     write(6,*) 'Error, Singular matrix'
     stop
   END IF
   VV(I) = 1.d0 / AMAX
 END DO ! i loop
 DO J=1,N
   DO I=1,J-1
     SUM = A(I,J)
     DO K=1,I-1
       SUM = SUM - A(I,K)*A(K,J) 
     END DO ! k loop
     A(I,J) = SUM
   END DO ! i loop
   AMAX = 0.d0
   DO I=J,N
     SUM = A(I,J)
     DO K=1,J-1
       SUM = SUM - A(I,K)*A(K,J) 
     END DO ! k loop
     A(I,J) = SUM
     DUM = VV(I)*DABS(SUM)
     IF(DUM.GE.AMAX) THEN
       IMAX = I
       AMAX = DUM
     END IF
   END DO ! i loop  
   IF(J.NE.IMAX) THEN
     DO K=1,N
       DUM = A(IMAX,K)
       A(IMAX,K) = A(J,K)
       A(J,K) = DUM
     END DO ! k loop
     D = -D
     VV(IMAX) = VV(J)
   END IF
   INDX(J) = IMAX
   IF(DABS(A(J,J)) < TINY) A(J,J) = TINY
   IF(J.NE.N) THEN
     DUM = 1.d0 / A(J,J)
     DO I=J+1,N
       A(I,J) = A(I,J)*DUM
     END DO ! i loop
   END IF 
 END DO ! j loop
 RETURN
 END subroutine LUDCMP

