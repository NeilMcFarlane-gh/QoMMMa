!------------------------------------------------------------------------------
! subroutine gen_cons
!
! Generates a constraints matrix in the primitive internal basis from the
! specification in icons(iseq,*)
!
! Arguments:
! cdat(np,nc): constraints matrix (out)
! icons(6,nc): constraints specification (in)
! nc:          number of constraints (in)
! np:          dimension of the primitive space
!------------------------------------------------------------------------------

  SUBROUTINE gen_cons(cdat,icons,nc,np)

! args
    INTEGER nc, np, icons(6,nc)
    REAL (rk), DIMENSION (np,nc) :: cdat

! local vars
    INTEGER ip, ic
! begin
    DO ic = 1, nc
      DO ip = 1, np
        cdat(ip,ic) = 0.0D0
      END DO
      cdat(icons(6,ic),ic) = 1.0D0
! jk: take care of bond difference
! the positive contribution is stored in icons(6,ic)
! the negative contribution is stored in icons(4,ic)      
      IF (icons(5,ic)==5) THEN
        cdat(icons(4,ic),ic) = -dsqrt(2.D0)
        cdat(icons(6,ic),ic) = dsqrt(2.D0)
      END IF
    END DO

! end gen_cons
  END SUBROUTINE gen_cons

!------------------------------------------------------------------------------
! subroutine proj_cons
!
! Projects the vectors in the matrix c into the space spanned by utmat
!
! Arguments:
! cdat(np,nc): unprojected constraints matrix (in) / projected matrix (out)
! utdat(m,np): m-dimensional space of vectors of dimension np (in)
! work(m):     scratch array, used for dp_{ic,j} (in)
! m:           dimension of the space spanned by utmat (nonredundant) (in)
! nc:          number of constraints (in)
! np:          dimension of the space in which utmat is represented (in)
!------------------------------------------------------------------------------

  SUBROUTINE proj_cons(cdat,utdat,work,m,nc,np)

! args
    INTEGER nc, np, m
    REAL (rk), DIMENSION (m) :: work
    REAL (rk), DIMENSION (np,nc) :: cdat
    REAL (rk), DIMENSION (m,np) :: utdat

! local vars
    INTEGER ic, ip, j

! begin, dp_{ic,j} = <C_ic|U_j>
    DO ic = 1, nc
      DO j = 1, m
        work(j) = 0.0D0
        DO ip = 1, np
          work(j) = work(j) + utdat(j,ip)*cdat(ip,ic)
        END DO
      END DO

! C_ic = 0
      DO ip = 1, np
        cdat(ip,ic) = 0.0D0
      END DO

! C_ic = sum_j dp_{ic,j}*U_j
      DO j = 1, m
        DO ip = 1, np
          cdat(ip,ic) = cdat(ip,ic) + work(j)*utdat(j,ip)
        END DO
      END DO
    END DO ! ic = 1,nc

! end proj_cons
  END SUBROUTINE proj_cons

!------------------------------------------------------------------------------
! subroutine merge_cons
!
! Composes a new matrix V out of [C,Ut]
!
! Arguments:
! v_mat:  matrix (np,m+nc) of [cmat,utmat transposed] (out)
! c_mat:  matrix (np,nc) of projected constraints (in)
! ut_mat: matrix (m,np) of U transposed (in)
! work:   work array(np) (in)
! m:      dimension of the space spanned by utmat (nonredundant) (in)
! nc:     number of constraints (in)
! np:     dimension of the space in which utmat is represented (in)
!------------------------------------------------------------------------------

  SUBROUTINE merge_cons(v_mat,c_mat,ut_mat,work,m,nc,np)

! args
    INTEGER m, nc, np
    TYPE (matrix) :: v_mat, c_mat, ut_mat
    REAL (rk), DIMENSION (np) :: work

! local vars
    INTEGER i, idum, j

! begin
    DO i = 1, nc
      idum = matrix_get_column(c_mat,size(work),work,i)
      idum = matrix_set_column(v_mat,size(work),work,i)
    END DO
    j = nc
    DO i = 1, m
      j = j + 1
      idum = matrix_get_row(ut_mat,size(work),work,i)
      idum = matrix_set_column(v_mat,size(work),work,j)
    END DO

! end merge_cons
  END SUBROUTINE merge_cons

!------------------------------------------------------------------------------
! subroutine ortho_mat
!
! Applies Schmidt orthogonalisation to the columns of vmat
! Taken from mankyopt: orthog.F
!
! Arguments:
! v_dat(np,nc+m): [cmat,utmat transposed] (in), orthogonalised (out)
! work(np):       work array (in)
! m, nc, np:      see above (in)
!------------------------------------------------------------------------------

  SUBROUTINE ortho_mat(v_dat,work,m,nc,np)

! args
    INTEGER m, nc, np
    REAL (rk), DIMENSION (np,nc+m) :: v_dat
    REAL (rk), DIMENSION (np) :: work

! local vars
    INTEGER i, j, k, nelem, nvec
    REAL (rk) dnorm, scapro, tol

! data
    DATA tol/1.0D-10/

! begin, orthogonalise vectors i = 2,nvec
    nvec = m + nc
    nelem = np
    DO i = 1, nvec
      DO j = 1, nelem
        work(j) = 0.0D0
      END DO

! make vector i orthogonal to vectors  k = 1,i-1
      DO k = 1, i - 1
        scapro = 0.0D0
        DO j = 1, nelem
          scapro = scapro + v_dat(j,i)*v_dat(j,k)
        END DO
        DO j = 1, nelem
          work(j) = work(j) + v_dat(j,k)*scapro
        END DO
      END DO

! subtract the collinear vector to make vector i orthogonal to k = 1,i-1
      DO j = 1, nelem
        v_dat(j,i) = v_dat(j,i) - work(j)
      END DO

! normalise vector i
      dnorm = 0.0D0
      DO j = 1, nelem
        dnorm = dnorm + v_dat(j,i)*v_dat(j,i)
      END DO
      IF (abs(dnorm)<tol) THEN
        dnorm = 0.0D0
      ELSE
        dnorm = 1.0D0/dsqrt(dnorm)
      END IF
      DO j = 1, nelem
        v_dat(j,i) = v_dat(j,i)*dnorm
      END DO
    END DO

! report resulting matrix
!   do i = 1,nvec
!      work(i) = 0.0D0
!      do j = 1,nelem
!         work(i) = work(i) + v_dat(j,i)*v_dat(j,i)
!      end do
!   end do

! end ortho_mat
  END SUBROUTINE ortho_mat

!------------------------------------------------------------------------------
! subroutine move_cons
!
! Moves the constraints behind the active space in the V matrix, transposes
! V and stores it in Ut
!
! Arguments: see merge_cons
!------------------------------------------------------------------------------

  SUBROUTINE move_cons(v_mat,ut_mat,work,m,nc,np)

! args
    INTEGER m, nc, np
    TYPE (matrix) :: v_mat, ut_mat
    REAL (rk), DIMENSION (np) :: work

! local vars
    INTEGER is, it, i
    REAL (rk) dnorm, tol

! data
    DATA tol/1.0D-10/

! begin
    DO is = 1, nc
      i = matrix_get_column(v_mat,size(work),work,is)
      i = matrix_set_row(ut_mat,size(work),work,m-nc+is)
    END DO
    it = 0
    DO is = nc + 1, nc + m
      i = matrix_get_column(v_mat,size(work),work,is)
      dnorm = 0.0D0
      DO i = 1, np
        dnorm = dnorm + work(i)*work(i)
      END DO
      IF (dnorm>tol) THEN
!         write(stdout,'("JK dnorm ",g15.7," tol ",g15.7)') dnorm,tol
        it = it + 1
        IF (it>m-nc) THEN
          WRITE (stdout,'(A,I5)') 'Too many active vectors, required: ', &
            m - nc
!            write(stdout,'("JK dnorm ",g15.7," tol ",g15.7)') dnorm,tol
          CALL hdlc_errflag('Constraints error','abort')
        ELSE
          i = matrix_set_row(ut_mat,size(work),work,it)
        END IF
      END IF
    END DO

! end move_cons
  END SUBROUTINE move_cons

!------------------------------------------------------------------------------
! subroutine unproj_cons
!
! Replaces projected constraints by the unprojected ones in the matrix
! U transposed
!
! Arguments: see merge_cons
!------------------------------------------------------------------------------

  SUBROUTINE unproj_cons(ut_mat,c_mat,work,m,nc,np)

! args
    INTEGER m, nc, np
    TYPE (matrix) :: ut_mat, c_mat
    REAL (rk), DIMENSION (np) :: work

! local vars
    INTEGER is, i

! begin
    DO is = 1, nc
      i = matrix_get_column(c_mat,size(work),work,is)
      i = matrix_set_row(ut_mat,size(work),work,m-nc+is)
    END DO

! end unproj_cons
  END SUBROUTINE unproj_cons

END MODULE dlfhdlc_constraint