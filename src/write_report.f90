SUBROUTINE write_report()
use nrtype ; use coordinates ; use optimdata
implicit none

! This subroutine initialises output.

integer(i4b) :: rstat, i, j, k, ii, img_num, conv_counter
real(sp) :: dist, angle, arccos, tot_angle, qe_no_disp
logical :: any
conv_counter=0


! Loop over all images
do img_num=1,nimg

	! Copy data of the image
	e=fulle(img_num)
	oe=fulloe(img_num)
	qe=fullqe(img_num)
	te=fullte(img_num)
	conv(:)=fullconv(img_num,:)
	convs(:)=fullconvs(img_num,:)
	converged=fullconverged(img_num)
	qe_no_disp=fullqend(img_num)

	if (ncon.gt.0) then
	totcnsen=fulltotcnsen(img_num)
	cnsval(:)=fullcnsval(img_num,:)
	cnsen(:)=fullcnsen(img_num,:)
	cnsg(:)=fullcnsg(img_num,:)
	end if

	open(unit=8,file=("add_to_report"//trim(img_string(img_num))),status="replace")
	print *, 'addtoreport created'
	write (unit=8,fmt='(A,I5)') "Situation at step:",nstep
	if ((nebtype.eq.4).or.(nebtype.eq.3).or.(nebtype.eq.5)) then
		if (line_search) then
			write (unit=8,fmt='(A)')  "Carrying on line search"
		end if
		if (reversing) then
			write (unit=8, fmt='(A)') " Reversing search direction"
		end if
		if (resetting) then 
			write (unit=8, fmt='(A)') " Resetting search direction to steepest descent"
		end if
	end if

	write (unit=8,fmt='(A,1X,2F13.6)') "Present, previous energy:",e,oe
	write (unit=8,fmt='(A25,1X,2F13.6)') "QM, MM, energies:",qe,te
	if (disp.eq.1) then
	   write (unit=8,fmt='(A,1X,F13.6)') "QM energy without dispersion:",qe_no_disp
	end if
	if (ncon.gt.0) then
	   write (unit=8,fmt='(A,1X,2F13.6)') "Constr. energy, Constraint-free energy:",totcnsen,e-totcnsen
	   write (unit=8,fmt='(A)') "Constr. no.;  Value  ; Ideal val. ; Energy contrib. ; Gradient"
	   do i=1,ncon
		   write (unit=8,fmt='(2X,I4,6X,F8.4,3X,F8.4,7X,F10.6,4X,F10.6)') i,cnsval(i),cnsidl(i), &
			& cnsen(i),cnsg(i)
	   end do
	end if

	if (climbing(img_num)) then
		write (unit=8,fmt=100) "Change in energy:",conv(1),tolde_org,convs(1)
		write (unit=8,fmt=100) "Maximum change of X:",conv(2),toldxmax_org,convs(2)
		write (unit=8,fmt=100) "RMS change of X:",conv(3),toldxrms_org,convs(3)
		write (unit=8,fmt=100) "Maximum gradient element:",conv(4),tolgmax_org,convs(4)
		write (unit=8,fmt=100) "RMS gradient element:",conv(5),tolgrms_org,convs(5)
	else
		write (unit=8,fmt=100) "Change in energy:",conv(1),tolde,convs(1)
		write (unit=8,fmt=100) "Maximum change of X:",conv(2),toldxmax,convs(2)
		write (unit=8,fmt=100) "RMS change of X:",conv(3),toldxrms,convs(3)
		write (unit=8,fmt=100) "Maximum gradient element:",conv(4),tolgmax,convs(4)
		write (unit=8,fmt=100) "RMS gradient element:",conv(5),tolgrms,convs(5)
	end if
	! If doing NEB, print RMS of perpendicular part - measure of convergence of elastic band
	if ((nebtype.eq.3).or.(nebtype.eq.5).or.(nebtype.eq.1).or.(nebtype.eq.2)) then
		write (unit=8,fmt=100) "RMS perpendicular component:",norm_per_force(img_num),tolper,convs(6)
	end if
	write (unit=8,fmt=*) ""

	100 FORMAT (A39,1X,F13.6,"    (",F8.6,") ",A3)


	if ((nimg.eq.1).and.(converged)) then
			write (unit=8,fmt='(A)') "The optimization has completed."
			write (unit=8,fmt='(A)') "Congratulations!!"
			stop
	end if

	! Don't check convergence for edge points - this will allow NEB to start at an arbitrary point on the path
	if (nimg.ne.1) then
		if (img_num.eq.1) then
			conv_counter=conv_counter+1
		else if (img_num.eq.nimg) then
			conv_counter=conv_counter+1
		else if (converged) then
			conv_counter=conv_counter+1
			write (unit=8,fmt='(A)') "The optimization converged"
		end if
	end if

	close(8)

	! Create update marker for bash script
	if (update_geom(img_num)) then
		open(unit=10,file="update_geom"//trim(img_string(img_num)),status="replace")
		close(10)
	end if

! End of loop over all images
end do

if ((nebtype.eq.1).or.(nebtype.eq.2).or.(nebtype.eq.3).or.(nebtype.eq.5)) then
    open(unit=9,file="add_to_master_report",status="replace")
    write(unit=9,fmt='(A,I5)') "QoMMMa step number ", nstep
    write(unit=9,fmt='(A20,A14,A2,A14,A11,A7,1X,A13,5A7,A14)') " ","Distance from", " ","Cos. of angle", &
          "Angle   ", "Total", "Present   ", "Energy", "Max ", "RMS ", "Max ", "RMS ", "Perpendicular"
    write(unit=9,fmt='(A20,A14,A2,A14,A11,A7,1X,A13,5A7,A14)') "","previous image", " ", "between images", &
           "in degrees", "angle", "energy   ", "change", "dX ", "dX ", "grad ", "grad ", "component"
    do img_num=1,nimg
        if (img_num.eq.1) then
            dist=0
            angle=0
            arccos=0
            tot_angle=0
        else
            dist=sqrt(sum((fullxopt(img_num,:)-fullxopt(img_num-1,:))**2))
            if (img_num.eq.nimg) then
                angle=0
                arccos=0
                !tot_angle=0
            else
! lines were splitted in two lines to compile in g95
                !tot_angle=dot_product(((fullxopt(1,:)-fullxopt(nimg,:))/sqrt(sum((fullxopt(1,:)-fullxopt(nimg,:))**2)))&
                 !  ,((fullxopt(img_num-1,:)-fullxopt(img_num,:))/sqrt(sum((fullxopt(img_num-1,:)-fullxopt(img_num,:))**2))))
                angle=dot_product(((fullxopt(img_num-1,:)-fullxopt(img_num,:))/sqrt(sum((fullxopt(img_num-1,:)-&
                  fullxopt(img_num,:))**2))),((fullxopt(img_num,:)-fullxopt(img_num+1,:))&
                        /sqrt(sum((fullxopt(img_num,:)-fullxopt(img_num+1,:))**2))))
                arccos=acos(angle)
            end if
            ! calculate total angle for last image as well
             tot_angle=dot_product(((fullxopt(1,:)-fullxopt(nimg,:))/sqrt(sum((fullxopt(1,:)-fullxopt(nimg,:))**2))),&
                ((fullxopt(img_num-1,:)-fullxopt(img_num,:))/sqrt(sum((fullxopt(img_num-1,:)-fullxopt(img_num,:))**2))))
        end if
        convs(:)=fullconvs(img_num,:)
        write (unit=9,fmt='(A,I4,A,F14.6, F14.3, F12.3, F9.3, F13.6,5A7, F13.6)') "Image number ", img_num, ": ",dist, angle, &
              arccos*180/pi, tot_angle, fulle(img_num),convs(1), convs(2), & 
        & convs(3), convs(4), convs(5), norm_per_force(img_num)
    end do
    write (unit=9,fmt=*) ""
    if (conv_counter.eq.nimg) then
        write (unit=9,fmt='(A)') "The optimization of all images has completed."
        write (unit=9,fmt='(A)') "Congratulations!!"
        stop
    end if
!   if (.not. sum(update_geom(:))) then
!   this is modified as like below to compile in G95
    any=.false.
    do i=1,nimg
       any=any.or.update_geom(i)
    enddo
    if(.not.any) then
        write (unit=9,fmt='(A)') "Geomtries of all images fixed - termination of Nudged Elastic Band calculation. "
        stop
    end if
end if

close(9)

END SUBROUTINE write_report

