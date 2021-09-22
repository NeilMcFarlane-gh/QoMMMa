SUBROUTINE collategeoms_simple()
use nrtype ; use coordinates ; use optimdata
implicit none

! This subroutine folds in two distinct geometries into a single new one.

integer(i4b) :: i, j, k, ii, jj, hh

outer: do i=1,n
     ! If this atom is part of the Hessian optimization group, copy xopt to x.
     do ii=1,nopt
         if (opt(ii).eq.i) then
              j=3*(i-1)+1
              jj=3*(ii-1)+1
              x(j:j+2)=xopt(jj:jj+2)
              cycle outer
         end if
     end do
     ! Otherwise, do nothing.
end do outer


END SUBROUTINE collategeoms_simple


