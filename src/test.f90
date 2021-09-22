SUBROUTINE update_conj_opt()
use nrtype ; use coordinates ; use optimdata
implicit none


! Specially Adapted Polak-Ribiere routine from Numerical Recipes
! Version for any number of images
! TO DO: checinkg for small displacement and freezeing of images
! Is it correct to use data of all images (including first and last
! for calculating search direction etc.?

integer(i4b) :: I,j, img_num
real(sp), allocatable :: DelG(:,:), ChgeX(:,:), DelX(:,:) 
real(sp), allocatable :: dire(:,:), norm_dire(:,:)
real(sp) ::  maxchgx, mod_dire, mod_delx, weight
real(sp),parameter ::  eps = 1.d-6  ! eps = 1.d-5 ! eps = 3.d-8 !
real(sp) :: gg, dgg, ll, dll, gam, check_restart, hh, ratio, init_step_size, reversing_step_size, resetting_step_size, new_step_size
logical :: line_search, reversing, resetting, if_climbing
allocate(delx(nimg,noptx),ChgeX(nimg,noptx),DelX(nimg,noptx),DelG(nimg,noptx))
allocate(dire(nimg,noptx),norm_dire(nimg,noptx))

update_geom=.true.
line_search=.true.  ! This value determines which geometry/gradient will be stored in Checkfile
reversing = .false.
resetting = .false.

if ((nebtype.eq.4).or.(nebtype.eq.3).or.(nebtype.eq.5)) then ! if minimization with conjugate gradient

    if ((nebtype.eq.3).or.(nebtype.eq.5)) then
        ! First image fixed!
        fulloe(1)=fulle(1)
        update_geom(1)=.false.
        ! Last image fixed!
        fulloe(nimg) = fulle(nimg)
        update_geom(nimg) = .false.
    end if

    IF (Nstep.eq.0) THEN ! first step at all
        init_step_size = 0.7
        ChgeX = -init_step_size*fulloptg(:,:) ! Initial, arbitrary step size, negative as we want to go downhill 
        ! If doing neb, fix the edges
        if ((nebtype.eq.3).or.(nebtype.eq.5)) then
            ! First image fixed!
            ChgeX(1,:) = 0
            ! Last image fixed!
            ChgeX(nimg,:) = 0
        end if
        weight = 1/init_step_size
        call check_step_length(chgex(:,:),weight) ! Do check for the first step - had problems with CM minimizations
        fullnewx(:,:) = fullxopt(:,:) + ChgeX
    else
      ! Find line search direction leading to this point
      ! as delx for first and last image will be always zero, they will be effectively fixed
        DelX = fullxopt(:,:) - fullox(:,:)
      ! Accordingly, dire for first and last will be zero
        dire = weight*delx  !Scale search direction vector with appropriate weight
        DelG = fulloptg(:,:) - fullog(:,:)
      !calculate some vector lengths
        mod_dire = 0.
        mod_delx = 0.
        do img_num=1,nimg     ! these should be null contributions from first and last image in neb
            mod_dire = mod_dire + sum(dire(img_num,:)**2)
            mod_delx = mod_delx + sum(delx(img_num,:)**2)
        end do 
        mod_dire = sqrt(mod_dire)
        mod_delx = sqrt(mod_delx)
        !norm_dire = dire/sqrt(sum(dire**2)) ! Normalized search direction
        norm_dire = dire/mod_dire ! Normalized search direction
      ! Changing dot product into sum of products
      ! Skipping first and last image
        if (nimg.ne.1) then
            hh = (sum(fulloptg(2:nimg-1,:)*norm_dire(2:nimg-1,:)) - sum(fullog(2:nimg-1,:)*norm_dire(2:nimg-1,:)))/mod_delx
            !hh = (sum(fulloptg(:,:)*norm_dire) - sum(fullog(:,:)*norm_dire))/sqrt(sum(delx**2))
        else
            hh = (sum(fulloptg(:,:)*norm_dire(:,:)) - sum(fullog(:,:)*norm_dire(:,:)))/mod_delx
        end if
        ! Scale hessian element if large
        !if (hh.gt.1) then
        !    hh = 0.1*hh
        !    write (*,*) " Scaling hh"
        !end if

      !Check the change in gradient along search direction
        !ratio=dot_product(fulloptg(1,:),norm_dire) / dot_product(fullog(1,:),norm_dire)
        ! Changing into sum again

      ! For all images, as norm_dire for first and last one will be zero
        ratio=sum(fulloptg(:,:)*norm_dire) / sum(fullog(:,:)*norm_dire) ! norm_dire or dire - wouldn't make a difference.
                ! In case you are doing SD, make it easier to accept
                !if ((ratio.gt. 0.33).or.(ratio.lt.-0.33)) then
        if_climbing=.false.
        do img_num=1,nimg
             if_climbing=if_climbing .or. climbing(img_num)
        end do
        ! This conditional statement was rewritten to avoid problems with g95 compiler
        !if (.not. sum(climbing) .and. (hh.lt.0.)) then ! Not good, check not necessary in CI-NEB
        if (.not. if_climbing .and. (hh.lt.0.)) then ! Not good, check not necessary in CI-NEB
            dire = -fullog(:,:) ! resetting search direction - steepest descent
            ! ChgeX for first and last image will be zero, as og for them is zero
            resetting_step_size = 0.3
            ChgeX = resetting_step_size  * dire ! For debugging purposes fix the guess Hessian for new direction
            weight = 1/resetting_step_size
            ! Check for small changes for each image?
             do img_num=2,nimg           
                 if (maxval(abs(ChgeX(img_num,:))).lt.(0.25*toldxmax_org)) then           
                     ChgeX(img_num,:) = 0.0d0           
                     update_geom(img_num)=.false.           
                     h = oh           
                 end if           
             end do           
            fullnewx(:,:) = fullox(:,:) + ChgeX(:,:)  !using old geometry
            !if (nstep.eq.1) then            ! if it happens to start with, use current geometry
            !    fullnewx(:,:) = fullxopt(:,:) + ChgeX(:,:)
            !end if
            resetting = .true.
            !call check_step_length(chgex,weight)
                !else if (ratio.gt. 0.33) then
        else if (ratio.gt. 0.6) then
                !if ((ratio.gt. 0.6).or.(ratio.lt.-0.33)) then
                !chgex = 0.7 * dire
            weight = hh
            ChgeX = 1/hh * dire ! for the first and last image this is equal to zero 
            !call check_step_length(chgex(:,:),weight)
            ! this is just for testing purposes!
            ! it should be fullox in this format!
            do img_num=2,nimg
                if (maxval(abs(ChgeX(img_num,:))).lt.(0.25*toldxmax_org)) then
                    ChgeX(img_num,:) = 0.0d0
                    update_geom(img_num)=.false.
                    h = oh
                end if
            end do
            !fullnewx(:,:) = fullxopt(:,:) + ChgeX(:,:)
            fullnewx(:,:) = fullox(:,:) + ChgeX(:,:)
        ! if ratio negative, we've overshoot
        else if (ratio.lt.-0.6) then
        !else if (ratio.lt.-0.33) then
            !if (ratio.lt.-5.) then  ! if horribly overshoot, take tiny step
            !    reversing_step_size = 0.1
            !else
                reversing_step_size = 0.5
            !end if
            ChgeX = reversing_step_size*DelX
            weight = weight/reversing_step_size
            !call check_step_length(chgex(:,:),weight)
            reversing = .true.
            fullnewx(:,:) = fullox(:,:) + ChgeX(:,:)
        else
            ! if gradient decreased sufficiently, construct new conjugate gradient
            ! do not include end-point images in here
            if (nimg.eq.1) then
                gg = sum(fullog(:,:)*fullog(:,:))
                dgg = sum((fulloptg(:,:) - fullog(:,:))* fulloptg(:,:)) ! Polak-Ribiere gradient, CHECK SIGNS
                ll = sum(fulloptg(:,:)*fulloptg(:,:))
                dll = sum(fullog(:,:)* fulloptg(:,:))
            else
                gg = sum(fullog(2:nimg-1,:)*fullog(2:nimg-1,:))
                dgg = sum((fulloptg(2:nimg-1,:) - fullog(2:nimg-1,:))* fulloptg(2:nimg-1,:)) ! Polak-Ribiere gradient, CHECK SIGNS
                ll = sum(fulloptg(2:nimg-1,:)*fulloptg(2:nimg-1,:))
                dll = sum(fullog(2:nimg-1,:)* fulloptg(2:nimg-1,:))
            end if
            if (gg == 0.0) return ! It should never happen 
            gam = dgg/gg
            check_restart=dll/ll
            ! this was initially 0.5, change to something larger, and reset to steepest descent
            ! changed, because it was (for ci-neb) resetting it almost all the time
            if (gam .gt. 1.0) gam=1.0  ! reduce contribution of old direction
            !if (gam .gt. 0.5) gam=0.  ! make steepest descent step in strange regions
            if ((check_restart .gt. 0.5) .or. (check_restart .lt. -0.5)) gam=0.0 
            !if (gam .gt. 0.5) gam=0.5 ! decrease gamma in strange regions
                !if (gam .gt. 0.5) gam=0. ! make steepest descent step in strange regions
            ! Added 10/10/05 JZ
            if (gam .lt. 0.) gam=0.
            ! Flag that you are doing SD somehow
            if (nimg.ne.1) then
                dire(2:nimg-1,:)=-fulloptg(2:nimg-1,:)+gam*dire(2:nimg-1,:)
            else
                dire = -fulloptg(:,:)+gam*dire 
            end if
            if (nimg.ne.1) then
                dire(1,:) = 0
                dire(nimg,:) = 0
            end if
            mod_dire = 0.
            do img_num=1,nimg
                   mod_dire = mod_dire + sum(dire(img_num,:)**2)
            end do
            mod_dire = sqrt(mod_dire)
            ! I don't know why I calculate norm_dire here ?? Probably to be able to calculated new hh?
            norm_dire = dire/mod_dire
            !ChgeX = 0.7  * dire ! For debugging purposes fix the guess Hessian for new direction
            new_step_size = 0.3
            !new_step_size = 0.2/hh
            !if (new_step_size .lt. 0.3
            !new_step_size = 1/hh
            ChgeX = new_step_size*dire
            weight = 1/new_step_size
            !call check_step_length(chgex(:,:),weight)
            line_search=.false. ! End of line search, new direction of the search
            ! Don't apply this check to climbing image; increased again from 0.1 to 0.25
            do img_num=2,nimg
                if ((maxval(abs(ChgeX(img_num,:))).lt.(0.25*toldxmax_org)).and.(.not.climbing(img_num))) then
                    ChgeX(img_num,:) = 0.0d0
                    update_geom(img_num)=.false.
                    h = oh
                end if
            end do
            fullnewx(:,:) = fullxopt(:,:) + ChgeX(:,:)
        end if
    end if
end if


! evaluate convergence tests
do img_num=1,nimg
    convs = " NO"
    converged=.false.
    i=0
    conv(1)=fulle(img_num)-fulloe(img_num)
    conv(2)=maxval(abs(ChgeX(img_num,:)))
    conv(3)=sqrt(sum(chgex(img_num,:)**2)/real(noptx,sp))
    conv(4)=maxval(abs(fulloptg(img_num,:)))
    conv(5)=sqrt(sum(fulloptg(img_num,:)**2)/real(noptx,sp))

    ! Tighter convergence test for climbing image
    if (climbing(img_num)) then
        if(abs(conv(1)) .lt. tolde_org) then
            convs(1)="YES" ; i=i+1
        end if
        if (conv(2) .lt. toldxmax_org) then
            convs(2)="YES" ; i=i+1
        end if
        if (conv(3) .lt. toldxrms_org) then
            convs(3)="YES" ; i=i+1
        end if
        if (conv(4) .lt. tolgmax_org) then
            convs(4)="YES" ; i=i+1
        end if
        if (conv(5) .lt. tolgrms_org) then
            convs(5)="YES" ; i=i+1
        end if
    else
        if (abs(conv(1)) .lt. tolde) then
             convs(1)="YES" ; i=i+1
        end if
        if (conv(2) .lt. toldxmax) then
             convs(2)="YES" ; i=i+1
        end if
        if (conv(3) .LT. toldxrms) then
             convs(3)="YES" ; i=i+1
        end if
        if (conv(4) .lt. tolgmax) then
             convs(4)="YES" ; i=i+1
        end if
        if (conv(5) .lt. tolgrms) then
             convs(5)="YES" ; i=i+1
        end if
        h = 0.0  ! No hessian matrix used in conjugate gradient minimization
    end if

    if (i .eq. 5) converged = .true.

    fullconverged(img_num)=converged
    fullconv(img_num,:)=conv(:)
    fullconvs(img_num,:)=convs(:)

    open(unit=8,file=("add_to_update"//trim(img_string(img_num))),status="replace")
    write (unit=8,fmt='(A,5F13.6)') " HH, gamma, weight, ratio, check_restart: ", hh, gam, weight, ratio, check_restart
    close(8)

end do

! increment nstep
nstep=nstep + 1

! If doing line search, keep reference geometry and gradient - redundant?(already present in write_checkfile ?)
!if ((line_search).and.(nstep.ne.1)) then
!    fullxopt(:,:)=fullox(:,:)
!    fulloptg(:,:)=fullog(:,:)
!    ! This line is not really necessary, but it would require altering files for reading and writing Checkfile
!    fullh(:,:,:)=fulloh(:,:,:)
!end if


return
END SUBROUTINE update_conj_opt

SUBROUTINE check_step_length(change, new_weight)
use nrtype ; use coordinates ; use optimdata
implicit none

real(sp) :: Change(nimg,noptx)
real(sp) :: new_weight
real(sp) :: stpl, lstep, stpmax

! I can't allocate it
!allocate(Change(nimg,noptx))

stpmax = STPMX * REAL(noptx,sp)

    stpl = SQRT(SUM(change**2))
    IF (stpl .gt. stpmax) THEN
        change = change / stpl * stpmax
        new_weight = new_weight * stpl/stpmax
        write (*,*) "Changing (1) Step Length"
    END IF

    lstep=maxval(abs(change))
    IF (lstep .gt. STPMX) THEN
        change = change / lstep * STPMX
        new_weight = new_weight*lstep/(stpmx)
        write (*,*)"Changing (2) Step Length"
    END IF

END SUBROUTINE check_step_length
