SUBROUTINE read_ab_initio_mecp()
use nrtype ; use coordinates ; use optimdata
implicit none

! This subroutine reads the "ab_initio" file containing the total energy
!   and the forces (in Hartree/bohr). These are converted to gradients in Hartree/Å.

integer(i4b) :: i, j, k, ii, img_num
character(20) :: dummy
real(sp) :: vec(3)

do img_num=1,nimg

open(unit=8,file=("ab_initio"//trim(img_string(img_num))),position="rewind")

! Start by setting qg to 0.d0 to avoid surprises!

qga=0.d0
qgb=0.d0
read(unit=8,fmt=*) dummy
read(unit=8,fmt=*) qea
read(unit=8,fmt=*) dummy
do i = 1, nq
    j = 3*(qm(i)-1)+1
    read(unit=8,fmt=*) ii, (vec(k),k=1,3)
    qga(j:j+2) = vec / bohr
end do

! Read in the QM gradient on the link atoms.
! These are "spread" between the corresponding QM and MM atoms.
! Given that rL = lratio * rMM + (1 - lratio) * rQM, the gradient is spread in the same way

do i = 1, nl
    read(unit=8,fmt=*) ii, (vec(k),k=1,3)
    j = 3*(links(i,1)-1)+1  ! the MM atom.
    k = 3*(links(i,2)-1)+1  ! the QM atom.
    qga(j:j+2) = lratio(i) * (vec)/bohr
    qga(k:k+2) = qga(k:k+2) + (1.d0-lratio(i)) * (vec)/bohr
end do
! for second state
read(unit=8,fmt=*) dummy
read(unit=8,fmt=*) qeb
read(unit=8,fmt=*) dummy
do i = 1, nq
    j = 3*(qm(i)-1)+1
    read(unit=8,fmt=*) ii, (vec(k),k=1,3)
    qgb(j:j+2) = vec / bohr
end do

! Read in the QM gradient on the link atoms.
! These are "spread" between the corresponding QM and MM atoms.
! Given that rL = lratio * rMM + (1 - lratio) * rQM, the gradient is spread in the same way

do i = 1, nl
    read(unit=8,fmt=*) ii, (vec(k),k=1,3)
    j = 3*(links(i,1)-1)+1  ! the MM atom.
    k = 3*(links(i,2)-1)+1  ! the QM atom.
    qgb(j:j+2) = lratio(i) * (vec)/bohr
    qgb(k:k+2) = qgb(k:k+2) + (1.d0-lratio(i)) * (vec)/bohr
end do

close(8)

fullqga(img_num,:)=qga(:)
fullqgb(img_num,:)=qgb(:)
fullqea(img_num)=qea
fullqeb(img_num)=qeb

! End of loop over images
end do
return

END SUBROUTINE read_ab_initio_mecp


