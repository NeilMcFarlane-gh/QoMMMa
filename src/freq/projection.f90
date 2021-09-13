SUBROUTINE projection(nx,mgrd,mhess,pmhess)
implicit none

INTEGER, INTENT(IN) :: nx
DOUBLE PRECISION, INTENT(IN) :: mgrd(nx),mhess(nx,nx)
DOUBLE PRECISION, INTENT(OUT) :: pmhess(nx,nx)
integer :: i, J, k
DOUBLE PRECISION :: drc(nx),dr_drt(nx,nx),unim(nx,nx),pmat(nx,nx)

! projection with respect to gradient
drc=mgrd/sqrt(sum((mgrd)**2))
do i=1,nx
  do j=1,nx
    dr_drt(i,j)=drc(i)*drc(j)
  enddo
enddo

unim=0.d0
do i=1,nx
  unim(i,i)=1.0d0
enddo 

do i=1,nx
  do j=1,nx
    pmat(i,j)=unim(i,j)-dr_drt(i,j)
  enddo
enddo

pmhess=matmul(pmat,matmul(mhess,transpose(pmat)))

End SUBROUTINE projection  


