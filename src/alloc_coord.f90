SUBROUTINE alloc_coord()
use nrtype ; use coordinates
implicit none

character(80) :: dummy
integer(i4b), parameter :: nconsmax_cart = 10 ! not a good idea to have more constraints than this
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
read(8,*) ncon_cart, kcnstyp 
if ((ncon_cart .lt. 0).or.(ncon_cart.gt.nconsmax_cart)) then
      write (*,*) "Too few/many constraints. ERROR"
      stop
else if ((kcnstyp.ne.1).and.(kcnstyp.ne.2)) then
      write (*,*) "Invalid constraint type. ERROR"
      stop
end if
read(8,*) dummy
read(8,*) ncon_prim
if (ncon_prim .lt. 0) then
	write (*,*) "Too few constraints. ERROR"
    stop
end if
read(8,*) dummy
read(8,*) nimg, nebtype, kspring
read(8,*) dummy
read(8,*) gsmtype
read(8,*) dummy
read(8,*) gsmphase
read(8,*) dummy
read(8,*) coordtype
read(8,*) dummy
read(8,*) primtype
close(8)

! then allocate

nx=3*n
nqx=3*nq
nlx=3*nl
noptx=3*nopt

! The number of atoms delocalised internal coordinates (if used) is initialised.
! By definition, this is 3N-6 (with N representing the number of atoms in the QM region - including link atoms!!)
if (coordtype .eq. 1) then
	ndlc = (nopt * 3) - 6
end if

allocate(qm(nq),nbonds(n),chg(n),bonds(n,maxbond),bonds_xopt(nopt,maxbond),attyp(n),modchg(n),inact(n),opt(nopt)) 
allocate(fullx(nimg,nx),fullxq(nimg,nqx),fullxl(nimg,nlx),fullxopt(nimg,noptx))
allocate(x(nx),xq(nqx),xl(nlx),xopt(noptx),lratio(nl),llabel(nl),label(n),links(nl,2)) 
allocate(qlabel(nq))
allocate(img_string(nimg))

END SUBROUTINE alloc_coord


