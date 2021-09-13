PROGRAM coker_watts_for_ground_state
implicit none

	! This program reads a hessian file and cart coords (in Angstrom) from a file
	! called 'hessian_input'.
	! It calculates vibrational frequencies.

INTEGER :: ns,n,i, j,k,jj,ll,l
character(len=2), allocatable :: at(:)
DOUBLE PRECISION, ALLOCATABLE :: Q(:),aw(:),h(:,:),mwhess(:,:),eigvx(:,:),eigs(:),grd(:),mwgrd(:)
DOUBLE PRECISION, ALLOCATABLE :: freqau(:),freqwn(:),eigv(:,:),pmhes(:,:)
CHARACTER(LEN=10), ALLOCATABLE :: imag(:)
DOUBLE PRECISION, PARAMETER :: au_to_wn = 219474.0

CHARACTER(LEN=2), allocatable :: lab(:)


call read_dimensions(ns)

n=3*ns
allocate(q(n),h(n,n),aw(n),at(ns),mwhess(n,n),eigvx(n,n),eigs(n),grd(n),mwgrd(n),pmhes(n,n))
allocate(freqau(n),freqwn(n),imag(n),lab(ns),eigv(n,n))

CALL ReadInput_prj(ns,n,q,h,aw,lab,grd)

write (*,*) "Input Read OK ..."

DO I = 1,n
        mwgrd(I)=grd(I)/sqrt(aw(I))
	DO J=1,I
		mwhess(I,J) = H(I,J) / SQRT(aw(I) * aw(J))
		mwhess(J,I) = mwhess(I,J)
	END DO
END DO
call write_hessian(n,mwhess)

CALL projection(n,mwgrd,mwhess,pmhes)
! Project the direction of the gradient out of the Hessian.
! Proj Hess = Pt x mwHess x P ; where P is the outer product of the mass-weighted gradient

CALL Diag(n,pmhes,EigV,eigs)
write (*,*) "Hessian diagonalized successfully..."
DO I = 1, N
	eigvx(:,I) = eigv(:,I) / aw
	eigvx(:,I) = eigvx(:,I) / SQRT(SUM(eigvx(:,I)**2))
        freqau(I) = SQRT(ABS(eigs(I)))
        IF (eigs(I) .lt. 0.d0) THEN
                imag(I) = " IMAGINARY"
        ELSE
                imag(I) = ""
        END IF
	freqwn(I) = freqau(I) * au_to_wn
END DO

OPEN(UNIT=8,FILE="freqs_results",ACTION="write")

WRITE (Unit=8,FMT=*) "RESULT of the FRequency calcs"
Write(Unit=8,FMT=*) ""
write (8,*) ""
write (8,*) "The frequencies: (au then wn)"
DO J=1,n
        write (UNIT=8,FMT='(F20.10,F20.3,A,A)') freqau(J), freqwn(J),"   ",imag(J)
END DO
write (8,*) ""
write (8,*) ""
write (8,*) "The Zero-point energy:"
write (8,*) ""
write (8,'(F20.10)') .5d0 * SUM(freqau(1:N-6))

write (8,*) "The Eigenvectors(as Nvectors of N numbers):"
write (8,*) ""
write (8,*) "A. In cartesians (Xi)"
write (8,*) ""
DO I = 1, N
	DO J = 1, N
		write (8,'(F20.10)') EigVX(J,I)
	END DO
END DO
write (8,*) ""
write (8,*) ""
write (8,*) "B. In mass-weighted cartesians (qi)"
write (8,*) ""
DO I = 1, N
        DO J = 1, N
                write (8,'(F20.10)') EigV(J,I)
        END DO
END DO

CLOSE(8)
OPEN(UNIT=8,FILE="freqs.xyz")
DO I = 1, N
        write (8,*) ns
	write (8,'(A,I4,A,F10.3,A)') "Mode number",I," Frequency:",freqwn(I),imag(I)
	DO J = 1, N/3
		JJ = 3 * (J - 1) + 1
		write (8,'(A2,3F12.8,3F8.4)') lab(j), (q(LL),LL=JJ,JJ+2), (EigV(LL,I),LL=JJ,JJ+2)
	END DO
END DO

CLOSE(8)

write (*,*) "Everything is now over !!"

END PROGRAM coker_watts_for_ground_state




