SUBROUTINE read_mulliken_mecp()
use nrtype ; use coordinates ; use optimdata
implicit none

! This subroutine reads the "mulliken" file containing the Mulliken charges
! for both states for MECP
!   on the QM atoms. Charges on H Link atoms are added in to the charges on the
!   corresponding QM atoms.

integer(i4b) :: i, j, k, img_num
real(dp) :: chgt
character(50) :: dummy
DOUBLE PRECISION :: mulla(nq), mullb(nq)

do img_num=1,nimg
open(unit=8,file=("mulliken"//trim(img_string(img_num))),position="rewind")

! Clear single-image table of charges
mull=0.d0
mulla=0.d0
mullb=0.d0

! read in the charges on QM atoms
! for first state
do i = 1, nq
    read(unit=8,fmt=*) mulla(i)
end do

! Read in the charges on link atoms; add to corresponding QM charge.

do i = 1, nl
    k=0
    inner: do j = 1, nq
        if (qm(j).eq.links(i,2)) then
            k = j
            exit inner
        end if
    end do inner
    if (k.eq.0) then
        write (*,*) "not found k in read_mulliken. Error."
        stop
    end if
    read(unit=8,fmt=*) chgt
    mull(k)=mulla(k)+chgt
end do
read(unit=8,fmt=*)dummy
read(unit=8,fmt=*)dummy

! for second state
do i = 1, nq
    read(unit=8,fmt=*) mullb(i)
end do

! Read in the charges on link atoms; add to corresponding QM charge.

do i = 1, nl
    k=0
    innerb: do j = 1, nq
        if (qm(j).eq.links(i,2)) then
            k = j
            exit innerb
        end if
    end do innerb
    if (k.eq.0) then
        write (*,*) "not found k in read_mulliken. Error."
        stop
    end if
    read(unit=8,fmt=*) chgt
    mull(k)=mullb(k)+chgt
end do
close(8)
do i = 1, nq
    mull(i)=(mulla(i)+mullb(i))/2
end do
! Copies Mulliken charges into all-image table

fullmull(img_num,:)=mull(:)

! End of loop over all images
end do
return

END SUBROUTINE read_mulliken_mecp


