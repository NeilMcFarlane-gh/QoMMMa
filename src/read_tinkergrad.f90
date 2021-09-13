SUBROUTINE read_tinkergrad()
use nrtype ; use coordinates ; use optimdata
implicit none

! This subroutine reads the "mm_grad" with gradients in kcal/mol/Å
!   they are then converted to Hartree/Å.

integer(i4b) :: i, j, k, ii, img_num
character(3) :: tmplab
character(10) :: dummy

!Loop over all images
do img_num=1,nimg

	open(unit=8,file=("mm_grad"//trim(img_string(img_num))),position="rewind")

	tg=0.d0
	read(unit=8,fmt=*) dummy, te
	do i = 1, n
		j = 3*(i-1)+1
		read(unit=8,fmt='(I5,1X,A3,2X,3F18.10)') ii, tmplab, (tg(k),k=j,j+2)
	end do

	close(8)

	te = te / hart_kcal
	tg = tg / hart_kcal

	fullte(img_num)=te
	fulltg(img_num,:)=tg(:)

! End of loop over all images
end do
return

END SUBROUTINE read_tinkergrad


