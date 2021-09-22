SUBROUTINE read_qmlatgrad_mecp()
use nrtype ; use coordinates ; use optimdata
implicit none

! This subroutine reads the "qmlatgrad.out" file containing the latice gradient
! for both state and takes the average
! These are converted to Hartree/Å, and assigned as qg()

integer(i4b) :: i, j, k, ii, jj, img_num, jfst,jend
character(20) :: dummy
real(sp) :: a, esp1, esp2, vec1(3), vec2(3), grad(3)
real(sp), parameter :: eps=1.d-4

do img_num=1,nimg

! Copies charges and gradients of the image to single-image tables
qg(:)=fullqg(img_num,:)

open(unit=8,file=("qmlatgrad"//trim(img_string(img_num))//".out"),position="rewind")

! reading for first state
read(unit=8,fmt=*) dummy
read(unit=8,fmt=*) i

atoms: do i = 1, n
    jfst=0
    do j = 1, nq
        if (qm(j) .eq. i) cycle atoms
! no resp contribution on qm atoms!
    end do
    do j = 1, nl
        if (links(j,1) .eq. i) cycle atoms
! no resp contribution on link atoms as charge = 0!
    end do
    if (inact(i)) cycle atoms
! no resp contribution on frozen MM atoms as we do not need the gradient.
    if (abs(chg(i)).lt.eps) cycle atoms
! no resp contribution to gradient on MM atoms with no charge.
    jfst=3*(i-1)+1
    read(unit=8,fmt=*) (qga(j),j=jfst,jfst+2)
    do k=jfst,jfst+2
      qga(k)=qga(k)/bohr
    enddo 
end do atoms

! for second state
read(unit=8,fmt=*) dummy
read(unit=8,fmt=*) i
atomsb: do i = 1, n
    jfst=0
    do j = 1, nq
        if (qm(j) .eq. i) cycle atomsb
! no resp contribution on qm atoms!
    end do
    do j = 1, nl
        if (links(j,1) .eq. i) cycle atomsb
! no resp contribution on link atoms as charge = 0!
    end do
    if (inact(i)) cycle atomsb
! no resp contribution on frozen MM atoms as we do not need the gradient.
    if (abs(chg(i)).lt.eps) cycle atomsb
! no resp contribution to gradient on MM atoms with no charge.
    jfst=3*(i-1)+1
    read(unit=8,fmt=*) (qgb(j),j=jfst,jfst+2)
    do k=jfst,jfst+2
      qgb(k)=qgb(k)/bohr
    enddo 
end do atomsb
close(8)

! taking average and assigning to qg
atomst: do i = 1, n
    jfst=0
    do j = 1, nq
        if (qm(j) .eq. i) cycle atomst
! no resp contribution on qm atoms!
    end do
    do j = 1, nl
        if (links(j,1) .eq. i) cycle atomst
! no resp contribution on link atoms as charge = 0!
    end do
    if (inact(i)) cycle atomst
! no resp contribution on frozen MM atoms as we do not need the gradient.
    if (abs(chg(i)).lt.eps) cycle atomst
! no resp contribution to gradient on MM atoms with no charge.
    jfst=3*(i-1)+1
    do k=jfst,jfst+2
      qg(k)=(qga(k)+qgb(k))/2
    enddo 
end do atomst

! Copies gradients back into all-image table
fullqg(img_num,:)=qg(:)

! End of the loop over all images
end do
return

END SUBROUTINE read_qmlatgrad_mecp


