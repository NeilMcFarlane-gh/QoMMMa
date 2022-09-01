subroUTINE update_opt_geometry_mecp()
use nrtype ; use coordinates ; use optimdata
implicit none


! Specially Adapted BFGS routine from Numerical Recipes

integer(i4b) :: I,j, img_num
real(dp) :: DelG(noptx), HDelG(noptx), ChgeX(noptx), DelX(noptx), w(noptx)
real(dp) :: fac, fad, fae, sumdg, sumdx, stpl, lstep, stpmax, maxchgx
real(dp),parameter ::  eps = 1.d-6  ! eps = 1.d-5 ! eps = 3.d-8 !
stpmax = STPMX * REAL(noptx,sp)

line_search=.false.
if ((nebtype.eq.4).or.(nebtype.eq.3).or.(nebtype.eq.5)) then
    call update_conj_opt()
else

    update_geom=.true.
    do img_num=1,nimg

    ! Copy data first
    oea=fulloea(img_num)
    oeb=fulloeb(img_num)
    ea=fullea(img_num)
    eb=fulleb(img_num)
    optg(:)=fulloptg(img_num,:)
    og(:)=fullog(img_num,:)
    xopt(:)=fullxopt(img_num,:)
    ox(:)=fullox(img_num,:)
    oh(:,:)=fulloh(img_num,:,:)

    if (nebtype.eq.0) then
        IF (Nstep.eq.0) THEN
            ChgeX = -0.7d0*optg
            h = oh
        else
            DelG = optg - og
            DelX = xopt - ox	! xi in NR
            HDelG = MATMUL(oh,DelG)
            fac = 1.d0 / SUM(DelG * DelX)	! NR has a check to see this is not too small
            fae = SUM(DelG * HDelG)
            fad = 1.d0 / fae
            sumdg = SUM(DelG**2)
            sumdx = SUM(DelX**2)
            w = fac * DelX - fad * HDelG	! see Schlegel p. 473; called dg in NR
            DO I = 1, noptx
                h(i,i:noptx) = oh(i,i:noptx) + fac * DelX(I) * DelX(I:noptx) - &
                & fad * HDelG(I) * HDelG(I:Noptx) + fae * w(I) * w(I:noptx)
                h(i:noptx,i) = H(I,I:noptx)
            END DO
            ChgeX = - MATMUL(h,optg)
        end if
    else
        if (img_num.eq.1) then
            WRITE (*,*)  'First image fixed'
            ChgeX = 0.0d0
            h = oh
            fulloea(img_num)=ea
            fulloeb(img_num)=eb
            update_geom(img_num)=.false.
        else if (img_num.eq.nimg) then
            WRITE (*,*)  'Last image fixed'
            ChgeX = 0.0d0
            h = oh
            fulloea(img_num)=ea
            fulloeb(img_num)=eb
            update_geom(img_num)=.false.
        else
            if (Nstep.eq.0) then
                ChgeX = -0.7d0*optg
                h = oh
            else
                DelG = optg - og
                DelX = xopt - ox	! xi in NR
                HDelG = MATMUL(oh,DelG)
                fac = SUM(DelG * DelX)	
                sumdg = SUM(DelG**2)
                sumdx = SUM(DelX**2)
    ! Test fac - skip update if not sufficiently positive - from NR
    !            if (fac**2.gt.(eps * sumdg * sumdx)) then
    !            if (abs(fac).gt.eps) then
                    fac = 1.d0 / fac
                    fae = SUM(DelG * HDelG)
                    fad = 1.d0 / fae
                    w = fac * DelX - fad * HDelG	! see Schlegel p. 473; called dg in NR
                    DO I = 1, noptx
                        h(i,i:noptx) = oh(i,i:noptx) + fac * DelX(I) * DelX(I:noptx) - &
                        & fad * HDelG(I) * HDelG(I:Noptx) + fae * w(I) * w(I:noptx)
                        h(i:noptx,i) = H(I,I:noptx)
                    END DO
                    ChgeX = - MATMUL(h,optg)
                    maxchgx=maxval(abs(ChgeX))
                    if (maxval(abs(ChgeX)).lt.(0.25*toldxmax_org)) then
                        ChgeX = 0.0d0
                        update_geom(img_num)=.false.
                        h = oh
                    end if
    !           else
    !           
    !                ChgeX = 0.0d0
    !                update_geom(img_num)=.false.
    !            end if
            end if
        end if 
    end if

    stpl = SQRT(SUM(ChgeX**2))
    IF (stpl .gt. stpmax) THEN
        ChgeX = ChgeX / stpl * stpmax
        write (*,*) "Changing (1) Step Length"
    END IF

    lstep=maxval(abs(chgex))
    IF (lstep .gt. STPMX) THEN
        ChgeX = ChgeX / lstep * STPMX
        write (*,*)"Changing (2) Step Length"
    END IF
    newx = xopt + ChgeX

    ! evaluate convergence tests

    convs = " NO"
    converged=.false.
    i=0
!   conv(1)=e-oe
    conv(1)=ea-eb
    conv(2)=maxval(abs(ChgeX))
    conv(3)=sqrt(sum(chgex**2)/real(noptx,sp))
    conv(4)=maxval(abs(optg))
    conv(5)=sqrt(sum(optg**2)/real(noptx,sp))

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
        if (conv(5) .lt. tolgrms) then
            convs(5)="YES" ; i=i+1
        end if
    else
        if (abs(conv(1)) .lt. tolde) then
             convs(1)="YES" ; i=i+1
        end if
        if (conv(2) .lt. toldxmax) then
             convs(2)="YES" ; i=i+1
        end if
        if (conv(3) .lt. toldxrms) then
             convs(3)="YES" ; i=i+1
        end if
        if (conv(4) .lt. tolgmax) then
             convs(4)="YES" ; i=i+1
        end if
        if (conv(5) .lt. tolgrms) then
             convs(5)="YES" ; i=i+1
        end if
    end if

    if (i .eq. 5) converged = .true.

    fullnewx(img_num,:)=newx(:)
    fullh(img_num,:,:)=h(:,:)
    fullconverged(img_num)=converged
    fullconv(img_num,:)=conv(:)
    fullconvs(img_num,:)=convs(:)

    open(unit=8,file=("add_to_update"//trim(img_string(img_num))),status="replace")
    !write (unit=8,fmt='(A,4F9.6,3F12.6)') " DelX,Delg, DelX*Delg, maxChgX, fac, fad, fae: ", sum(DelX), sum(DelG), sum(DelX*DelG), maxval(abs(ChgeX)), fac, fad, fae
    write (unit=8,fmt='(A,3F9.6, F13.6)') " DelX, Delg, DelX*Delg, maxChgX: ", sum(DelX), sum(DelG), sum(DelX*DelG), maxchgx
    close(8)

!End of loop over images
end do
! increment nstep
nstep=nstep + 1

end if
return
END SUBROUTINE update_opt_geometry_mecp
