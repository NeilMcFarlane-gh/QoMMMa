SUBROUTINE write_hessian(nx,h)
implicit none

INTEGER, INTENT(IN) :: nx
DOUBLE PRECISION, INTENT(IN) :: h(nx,nx)

integer :: I, ii, J, jj, kk, K, nxt

OPEN(UNIT=8,FILE="test_write_hessian")

do i=1, nx
	j=i/5
	do k=1,j
		write (8,'(I3,5E14.6)') i, (h(i,ii),ii=5*k-4,5*k)
	end do
	if (i.ne.5*j) then
		write (8,'(I3,5E14.6)') i, (h(i,ii),ii=5*j+1,i)
	end if
end do

CLOSE(8)


END SUBROUTINE write_hessian


