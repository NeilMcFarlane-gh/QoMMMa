SUBROUTINE short_read_checkfile()
use nrtype ; use coordinates ; use optimdata
implicit none

! This subroutine reads the "CheckFile" briefly, so as to get the atomic charges.
! It can only work after "alloc_coord" - so one needs to check that the arrays are allocated.

character(3) :: tmplab
character(80) :: dummy
integer(i4b) :: rstat
integer(i4b) :: i, j, k, ii, jj, img_num, q
real(sp) :: bb

! Loop over all images
do img_num=1,nimg
	
	! First check whether the arrays are in fact already allocated...

	open(unit=8,file=("CheckFile"//trim(img_string(img_num))))

	read(unit=8,fmt=*) dummy
	read(unit=8,fmt=*) dummy
	read(unit=8,fmt=*) i
	read(unit=8,fmt=*) dummy
	read(unit=8,fmt=*) i, j, k, ii

	! Dispersion correction for QM atoms

	read(unit=8,fmt=*) dummy
	read(unit=8,fmt=*) disp

	read(unit=8,fmt=*) dummy
	! This line actually contains data about images, but processing is not needed
	! Image data is extracted in alloc_coord
	if (nstep.eq.0) then
		read(unit=8,fmt=*) dummy
	else if ((nebtype.eq.4).or.(nebtype.eq.3).or.(nebtype.eq.5)) then
		read(unit=8,fmt=*) bb, bb, bb, weight
	else
		read(unit=8,fmt=*) dummy
	end if
	
	! GSM
	read(unit=8,fmt=*) dummy
	read(unit=8,fmt=*) gsmtype

	! Coordinate selection
	read(unit=8,fmt=*) dummy
	read(unit=8,fmt=*) coordtype
	
	! Primitive internal coordinate definitions.
	if (coordtype .eq. 1) then
		read(unit=8,fmt=*) dummy
		read(unit=8,fmt=*) nprim
		do q = 1, nprim
			read (unit=8,fmt=*) (prim_list(q,j),j=1,2)
		end do
	end if

	! Now number & type of constraints
	read(unit=8,fmt=*) dummy
	read(unit=8,fmt=*) ncon, kcnstyp

	! Now read through list of QM atoms.
	read(unit=8,fmt=*) dummy
	do i=1,nq
	!     read(unit=8,fmt=*) opt(i)
		 read(unit=8,fmt=*,IOSTAT=rstat) j, qlabel(i)
		 qm(i) = j
		 opt(i) = j
	end do

	! Now read list of Link atoms and their details.

	read(unit=8,fmt=*) dummy
	! read through links: MM atom; QM atom; rQL/rQM ideal ratio; label.
	if (nl.gt.0) then
	   do i = 1, nl
	!        read(unit=8,fmt=*) opt(nq+i)
			read(unit=8,fmt=*) (links(i,j),j=1,2), lratio(i), llabel(i)
			k = 1
			do j = 1,nq
				if (qm(j).eq.links(i,2)) then
					k = 0
					exit
				end if
			end do
			if (k.ne.0) then
				write (*,*) "Link atom not connected to a QM atom. ERROR."
				close(8)
				stop
			end if
			k = 0
			do j = 1,nq
				if (qm(j).eq.links(i,1)) then
					k = 1
					exit
				end if
			end do
			if (k.ne.0) then
				write (*,*) "Link atom is a QM atom. ERROR."
				close(8)
				stop
			end if
			opt(nq+i)=links(i,1)
	   end do
	end if

	! Now read list of atoms included in Hessian Optimization which are neither QM nor link (if any).
	read(unit=8,fmt=*) dummy
	ii = nopt-nq-nl
	if (ii.gt.0) then
		 do i = 1, ii
			  read (unit=8,fmt=*) opt(nq+nl+i)
		  end do
	end if

	! Now read in details of any constraint to apply to the geometry

	read(unit=8,fmt=*) dummy
	read(unit=8,fmt=*) dummy
	read(unit=8,fmt=*) dummy
	read(unit=8,fmt=*) dummy

	if (ncon.gt.0) then
	! otherwise move on....
	cnsat=0
	do i=1,ncon
			read(unit=8,fmt=*) cnstyp(i),kcns(i),cnsidl(i)
			if ((cnstyp(i).lt.1).or.(cnstyp(i).gt.4)) then
					write (*,*) "Impossible Bond Constraint - Should be between 1-4. ERROR."
					close(8)
					stop
			end if
			if (kcns(i).lt.0.) then
					write (*,*) "Constraint constant kcns negative. ERROR."
					close(8)
					stop
			end if
			if (cnstyp(i).eq.1) ncnsat(i)=2
			if ((cnstyp(i).eq.2).or.cnstyp(i).eq.3) ncnsat(i)=4
			if (cnstyp(i).eq.4) ncnsat(i)=6
			read(unit=8,fmt=*) (cnsat(i,j),j=1,ncnsat(i))
	end do
	! Check that all constrained atoms are QM, Link or HessOpt
	! This is not necessary in this version of file.
	! call check_constrained_atoms()
	end if

	! Read list of atomic charges.
	read(unit=8,fmt=*) dummy
	modchg=.false.
	do i = 1,n
		 read (unit=8,fmt=*) j,chg(i), ii
		 if (ii .eq. 1) modchg(i) = .true.
	end do

	! Read in the list of MM atoms which are frozen.
	! Not allowed for QM or MM atoms, of course!
	read(unit=8,fmt=*) dummy
	read(unit=8,fmt=*) ninact
	read(unit=8,fmt=*) dummy
	inact = .false.
	do i = 1, ninact
		 read (unit=8,fmt=*) jj
		 do j = 1, nq
			 if (qm(j) .eq. jj) then
				   write (*,*) "Error, Trying to inactivate a QM atom :",qm(j)
				   close(8)
				   stop
			  end if
		 end do
		 do j = 1, nl
			 if (links(j,1) .eq. jj) then
				   write (*,*) "Error, Trying to inactivate a link atom :",links(j,1)
				   close(8)
				   stop
			  end if
		 end do
		 inact(jj) = .true.
	end do

	! that's all that's needed!

	close(8)

	inquire (file="update_geom"//trim(img_string(img_num)), exist=update_geom(img_num))

! End of loop over all images
end do
END SUBROUTINE short_read_checkfile


