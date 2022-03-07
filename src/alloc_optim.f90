SUBROUTINE alloc_optim()
use nrtype ; use optimdata ; use coordinates
implicit none

character(80) :: dummy

! First check whether the arrays are in fact already allocated...

IF (allocated(qg).or.allocated(tg)) THEN
      write (*,*) "Arrays already allocated in alloc_optim. ERROR."
      stop
END IF

allocate(qg(nx),tg(nx),g(nx),optg(noptx),og(noptx),h(noptx,noptx), &
      & oh(noptx,noptx),ox(noptx),newx(noptx),mull(nq),x_copy(noptx))
allocate(fulloe(nimg),fullox(nimg,noptx),fullog(nimg,noptx), &
      & fulloh(nimg,noptx,noptx), fullte(nimg),fulltg(nimg,nx))
allocate(fullqg(nimg,nx),fullqe(nimg),fullmull(nimg,nq),fulle(nimg), &
      & fulltotcnsen(nimg),fulloptg(nimg,noptx),norm_per_force(nimg))
allocate(fullconv(nimg,5),fullconverged(nimg),fullh(nimg,noptx,noptx), &
      & fullnewx(nimg,noptx),fullconvs(nimg,6),update_geom(nimg),climbing(nimg))
allocate(dispgrad(nx))
allocate(fullqend(nimg))

if (ncon.gt.0) then
      allocate(cnstyp(ncon),cnsat(ncon,maxcnsat),kcns(ncon),cnsval(ncon), &
         & ncnsat(ncon),cnsidl(ncon),cnsen(ncon),cnsg(ncon), & 
         & fullcnsen(nimg,ncon),fullcnsg(nimg,ncon),fullcnsval(nimg,ncon))
end if

! Assign values for convergence tests 
! TO-DO : Need to check what convergence criteria are required for gsm.
if (nimg.eq.1) then
   tolde=tolde_org
   tolgmax=tolgmax_org 
   tolgrms=tolgrms_org 
   toldxmax=toldxmax_org 
   toldxrms=toldxrms_org
else if (coordtype .eq. 1) then
   tolde=tolde_org*2
   tolgmax=tolgmax_org*2
   tolgrms=tolgrms_org*2
   toldxmax=toldxmax_org*2
   toldxrms=toldxrms_org*2
else
   tolde=tolde_org*2
   tolgmax=tolgmax_org*2
   tolgrms=tolgrms_org*2
   toldxmax=toldxmax_org*2
   toldxrms=toldxrms_org*2
end if

return

END SUBROUTINE alloc_optim


