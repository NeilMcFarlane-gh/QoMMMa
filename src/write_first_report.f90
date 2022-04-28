subroutine write_first_report()
use nrtype ; use coordinates; use optimdata
implicit none

integer(i4b) :: i,j,ccl(4), img_num

! Start loop over all images
do img_num=1,nimg

	open(unit=8,file=("report_header"//trim(img_string(img_num))),status="replace")

	ccl = 0
	write (unit=8,fmt='(A)') ""
	write (unit=8,fmt='(A)') "       ****************************"
	write (unit=8,fmt='(A)') "                  QoMMMa"
	write (unit=8,fmt='(A)') "                  v.8.07"
	write (unit=8,fmt='(A)') "           J.N. Harvey, KU Leuven"
	write (unit=8,fmt='(A)') "       Modified version 16 Aug 2021"
	write (unit=8,fmt='(A)') "            Please see file" 
	write (unit=8,fmt='(A)') "     changes_made_to_qommma_aug2021.txt"
	write (unit=8,fmt='(A)') "            for more details"
	write (unit=8,fmt='(A)') "       ****************************"
	write (unit=8,fmt='(A)') ""
	write (unit=8,fmt='(A)') "     Using a modified version of Tinker"
	write (unit=8,fmt='(A)') "          Version 6.1.01"
	write (unit=8,fmt='(A)') "   Copyright 1990-2012 Jay William Ponder"
	write (unit=8,fmt='(A)') "      see http://dasher.wustl.edu/tinker/"
	write (unit=8,fmt='(A)') ""
	write (unit=8,fmt='(A,I5,A,I3,A)') "The system is made up of ",n," atoms of which ",nq," are QM."
	write (unit=8,fmt='(A,I4,A)') "There are ",nl," link atoms connecting QM and MM parts."
	write (unit=8,fmt='(A,I4,A)') "A total of ",nopt," atoms will be optimized using the QoMMMa Hessian"
	write (unit=8,fmt='(A,I4,A)') "This is equivalent to ",nopt-nq-nl," atoms on top of the QM & Link ones"

	if (nimg.gt.1) then
	   write (unit=8,fmt='(A,1X,I4,A)') "You have ",nimg," images of your system."
	   write (unit=8,fmt='(A,1X,I4)') "This is a report file for image", img_num
	end if

	if (ncon_cart.gt.0) then
	   write (unit=8,fmt='(A,1X,I4,A)') "You have put",ncon_cart," cartesian constraints on the geometry."
	end if

	if (disp.eq.1) then
	   write (unit=8,fmt='(A)') "You have chosen to calculate an empirical dispersion correction for the QM atoms" 
	   write (unit=8,fmt='(A)') "Currently only available for B3LYP and BP86 (both have scaling factor s6 = 1.05)"
	   write (unit=8,fmt='(A)') "Uses the equations and constants from Grimme, J. Comp. Chem. 2006 (DFT-D v2)"
	end if

	write (unit=8,fmt='(A)') ""

	close(8)

! End of loop over all images
end do

! Now write a header for master_report
if (nimg.gt.1) then
	open(unit=8,file="master_report",status="replace")
	write (unit=8,fmt='(A)') " This file contains summary data of a nudged elastic band simulation"
	write (unit=8,fmt='(A)') " After each step of QoMMMa, for each image following data is printed: "
	write (unit=8,fmt='(A)') " distance from previous image, cosine of angle formed by an image with"
	write (unit=8,fmt='(A)') " its two neighbouring images, next this angle given in degrees, "
	write (unit=8,fmt='(A)') " cosine of an angle formed by current and previous image with direction "
	write (unit=8,fmt='(A)') " of the whole pathway, then total energy, and convergence flags for "
	write (unit=8,fmt='(A)') " energy, maximum and average displacement, maximum and average gradient"
	write (unit=8,fmt='(A)') ""
	write (unit=8,fmt='(A)') ""

	close(8)
end if

END subroutine write_first_report
