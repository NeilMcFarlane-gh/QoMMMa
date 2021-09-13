SUBROUTINE alloc_coord()
use nrtype ; use coordinates
implicit none

character(80) :: dummy
integer(i4b), parameter :: nconsmax = 10 ! not a good idea to have more constraints than this
integer(i4b) :: disp

! First check whether the arrays are in fact already allocated...

IF (allocated(qm).or.allocated(bonds).or.allocated(attyp) &
   & .or.allocated(x).or.allocated(xq).or.allocated(xl).or.allocated(lratio) &
   & .or.allocated(label).or.allocated(nbonds).or.allocated(links).or.allocated(xopt) &
   & .or.allocated(opt).or.allocated(llabel).or.allocated(chg).or.allocated(modchg)) THEN
      write (*,*) "Arrays already allocated in alloc_coord. ERROR."
      stop
END IF

open(unit=8,file="fortinput",position="rewind")
read(8,*) dummy
read(8,*) n, nq, nl, nopt
read(8,*) dummy
read(8,*) disp
read(8,*) dummy
read(8,*) ncon, kcnstyp
if ((ncon .lt. 0).or.(ncon.gt.nconsmax)) then
      write (*,*) "Too few/many constraints. ERROR"
      stop
else if ((kcnstyp.ne.1).and.(kcnstyp.ne.2)) then
      write (*,*) "Invalid constraint type. ERROR"
      stop
end if
read(8,*) dummy
read(8,*) nimg, nebtype, kspring
read(8,*) dummy
read(8,*) gsmtype
read(8,*) dummy
read(8,*) coordtype
close(8)

! then allocate

nx=3*n
nqx=3*nq
nlx=3*nl
noptx=3*nopt

allocate(qm(nq),nbonds(n),chg(n),bonds(n,maxbond),attyp(n),modchg(n),inact(n),opt(nopt)) 
allocate(fullx(nimg,nx),fullxq(nimg,nqx),fullxl(nimg,nlx),fullxopt(nimg,noptx))

! If DLC are used (either pure or hybrid), then the appropriate arrays are allocated.
! The length of a DLC array is (by definition) 3N - 6.
IF ((coordtype==1).or.(coordtype==2)) THEN
	allocate(fullx_dlc(nimg,nx-6),fullxq_dlc(nimg,nqx-6),x_dlc(nx-6),xq_dlc(nqx-6))
END IF

allocate(x(nx),xq(nqx),xl(nlx),xopt(noptx),lratio(nl),llabel(nl),label(n),links(nl,2)) 
allocate(qlabel(nq))
allocate(img_string(nimg))

END SUBROUTINE alloc_coord


