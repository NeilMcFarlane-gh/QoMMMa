SUBROUTINE collategeoms()
use nrtype ; use coordinates ; use optimdata
implicit none

! This subroutine folds in two distinct geometries into a single new one.

integer(i4b) :: i, j, k, ii, jj, hh, img_num

! Loop over all images
do img_num=1, nimg
if (update_geom(img_num)) then

xopt(:)=fullxopt(img_num,:)
x(:)=fullx(img_num,:)
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

fullx(img_num,:)=x(:)

end if

! End of loop over all images
end do

END SUBROUTINE collategeoms


