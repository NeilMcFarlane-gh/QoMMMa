SUBROUTINE update_opt_geometry()
use nrtype ; use coordinates ; use optimdata
implicit none


! Specially Adapted BFGS routine from Numerical Recipes

integer(i4b) :: I,j, img_num
real(dp) :: DelG(noptx), HDelG(noptx), ChgeX(noptx), DelX(noptx), w(noptx)
real(dp) :: fac, fad, fae, sumdg, sumdx, stpl, lstep, stpmax
real(dp),parameter ::  eps = 3.d-8 !eps = 1.d-5 
stpmax = STPMX * REAL(noptx,sp)

update_geom=.true.
do img_num=1,nimg

! Copy data first
oe=fulloe(img_num)
e=fulle(img_num)
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
            & fad * HDelG(I) * HDelG(I:Nx) + fae * w(I) * w(I:noptx)
            h(i:noptx,i) = H(I,I:noptx)
        END DO
        ChgeX = - MATMUL(h,optg)
    end if
else
    if (img_num.eq.1) then
        WRITE (*,*)  'First image fixed'
        ChgeX = 0.0d0
        h = oh
        fulloe(img_num)=e
        update_geom(img_num)=.false.
    else if (img_num.eq.nimg) then
        WRITE (*,*)  'Last image fixed'
        ChgeX = 0.0d0
        h = oh
        fulloe(img_num)=e
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
            if ((fac**2).gt.(eps * sumdg * sumdx)) then
            ! if (abs(fac).gt.eps) then
                fac = 1.d0 / fac
                fae = SUM(DelG * HDelG)
                fad = 1.d0 / fae
                w = fac * DelX - fad * HDelG	! see Schlegel p. 473; called dg in NR
                DO I = 1, noptx
                    h(i,i:noptx) = oh(i,i:noptx) + fac * DelX(I) * DelX(I:noptx) - &
                    & fad * HDelG(I) * HDelG(I:Nx) + fae * w(I) * w(I:noptx)
                    h(i:noptx,i) = H(I,I:noptx)
                END DO
                ChgeX = - MATMUL(h,optg)
            else
                ChgeX = 0.0d0
                update_geom(img_num)=.false.
            end if
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
conv(1)=e-oe
if (abs(conv(1)) .lt. tolde) then
     convs(1)="YES" ; i=i+1
end if
conv(2)=maxval(abs(ChgeX))
if (conv(2) .lt. toldxmax) then
     convs(2)="YES" ; i=i+1
end if
conv(3)=sqrt(sum(chgex**2)/real(noptx,sp))
if (conv(3) .lt. toldxrms) then
     convs(3)="YES" ; i=i+1
end if
conv(4)=maxval(abs(optg))
if (conv(4) .lt. tolgmax) then
     convs(4)="YES" ; i=i+1
end if
conv(5)=sqrt(sum(optg**2)/real(noptx,sp))
if (conv(5) .lt. tolgrms) then
     convs(5)="YES" ; i=i+1
end if

if (i .eq. 5) converged = .true.

fullnewx(img_num,:)=newx(:)
fullh(img_num,:,:)=h(:,:)
fullconverged(img_num)=converged
fullconv(img_num,:)=conv(:)
fullconvs(img_num,:)=convs(:)

open(unit=8,file=("add_to_update"//trim(img_string(img_num))),status="replace")
write (unit=8,fmt='(A,4F15.6)') " Delg*DelX, fac, fad, fae: ", sum(DelG*DelX), fac, fad, fae
close(8)

!End of loop over images
end do
! increment nstep
nstep=nstep + 1

return
END SUBROUTINE update_opt_geometry
