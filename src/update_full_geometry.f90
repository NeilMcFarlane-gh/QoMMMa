SUBROUTINE update_full_geometry()
use nrtype ; use coordinates ; use optimdata
implicit none
                                                                                                                                             
integer(i4b) :: i, j, k, img_num
character(80) :: dummy
                              
! Loop over all images
do img_num=1,nimg
    if (update_geom(img_num)) then
        newx(:)=fullnewx(img_num,:)
        x(:)=fullx(img_num,:)
        xopt(:)=fullxopt(img_num,:)
    ! for the remaining parts of the program, update "xopt" to "newx"
    ! Also update the relevant coordinates within "x"
        xopt=newx
        
        do i = 1, nopt
            j = 3*i-2
            k = 3*opt(i)-2
            x(k:k+2) = xopt(j:j+2)
        end do
        fullxopt(img_num,:)=xopt(:)
        fullx(img_num,:)=x(:)
    end if
end do
return
                                                                                                                                             
END SUBROUTINE update_full_geometry



