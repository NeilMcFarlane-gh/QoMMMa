SUBROUTINE read_checkfile()
use nrtype ; use coordinates ; use optimdata
implicit none

! This subroutine reads the "CheckFile" with details of the QMMM job
! It can only work after "alloc_coord" - so one needs to check that the arrays are allocated.

character(3) :: tmplab
character(80) :: dummy
integer(i4b) :: rstat
integer(i4b) :: i, j, k, ii, jj, img_num
real(sp) :: bb

! Start loop over all images
do img_num=1,nimg

	! First check whether the arrays are in fact already allocated...

	open(unit=8,file=("CheckFile"//trim(img_string(img_num))))

	read(unit=8,fmt=*) dummy
	read(unit=8,fmt=*) dummy
	read(unit=8,fmt=*) nstep
	read(unit=8,fmt=*) dummy
	read(unit=8,fmt=*) i, j, k, ii

	if ((i.ne.n).or.(j.ne.nq).or.(k.ne.nl).or.(ii.ne.nopt)) THEN
		  write (*,*) "Mismatch in number of atoms in read_checkfile. ERROR."
		  close(8)
		  stop
	end if

	! Dispersion correction

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

	! GSM type

	read(unit=8,fmt=*) dummy
	read(unit=8,fmt=*) gsmtype
	
	! Coordinate selection for optimisation

	read(unit=8,fmt=*) dummy
	read(unit=8,fmt=*) coordtype

	! Number and overall type of restraints

	read(unit=8,fmt=*) dummy
	read(unit=8,fmt=*) ncon,kcnstyp

	! Now read list of QM atoms.
	read(unit=8,fmt=*) dummy
	do i=1,nq
		 read(unit=8,fmt=*,IOSTAT=rstat) j, qlabel(i)
		 if (rstat .ne. 0) then
			 write (*,*) "error reading checkfile. ERROR."
			 close(8)
			 stop
		 end if
		 if (j .gt. n) then
			 write (*,*) "impossible number for qm atom. ERROR."
			 close(8)
			 stop
		 end if
		 qm(i) = j
		 opt(i) = j
	end do

	! Now read list of Link atoms and their details.

	read(unit=8,fmt=*) dummy
	! read links: MM atom; QM atom; rQL/rQM ideal ratio; label.
	if (nl.gt.0) then
	   do i = 1, nl
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
			  read (unit=8,fmt=*) jj
			  k = 0
			  do j = 1, nq
				  if (jj.eq.qm(j)) then
						k = 1
						exit
				  end if
			  end do
			  if (k.ne.0) then
				  write (*,*) "QM atoms are automatically included in Hessian Optimization."
				  write (*,*) "You need not specify them in the list of atoms. ERROR."
				  stop
			  end if
			  do j = 1, nl
				  if (jj.eq.links(j,1)) then
						k = 1
						exit
				  end if
			  end do
			  if (k.ne.0) then
				  write (*,*) "Link atoms are automatically included in Hessian Optimization."
				  write (*,*) "You need not specify them in the list of atoms. ERROR."
				  stop
			  end if
			  opt(nq+nl+i)=jj
		  end do
	end if

	! Now read in details of constraints to apply to the geometry
																												 
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
	call check_constrained_atoms()
	end if
																															   
	! Read list of atomic charges.
	read(unit=8,fmt=*) dummy
	modchg=.false.
	! Clear single-image table of charges first
	chg(:)=0.d0
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

	! copy over the coordinates of the "opt" atoms.
	! Copy coordinates first and clear the xopt, ox, og  and oh tables first

	oe=0.d0
	xopt=0.d0
	ox = 0.d0
	og = 0.d0
	oh = 0.d0
	x = 0.d0
	x(:)=fullx(img_num,:)
	do i = 1, nopt
		j=3*(i-1)+1
		k=3*(opt(i)-1)+1
		xopt(j:j+2)=x(k:k+2)
	end do

	if (nstep .eq. 0) then
		 do i = 1, noptx
			 oh(i,i) = 0.7d0
		 end do
		 ox = xopt
	else
		 ! start reading: previous energy; previous geometry; previous gradient; previous inverse Hessian.
		 ! Only one half of the inverse Hessian is saved & read. Only "opt" atoms are included.
		 ! Note that the current geometry has already been read in from "geom_expl.xyz".
		 read (unit=8,fmt=*) dummy
		 read (unit=8,fmt='(F15.8)') oe
		 read (unit=8,fmt=*) dummy
		 read (unit=8,fmt='(3F12.6)') (ox(i),i=1,noptx)
		 read (unit=8,fmt=*) dummy
		 read (unit=8,fmt='(3F12.6)') (og(i),i=1,noptx)
	!     if (nebtype.eq.4) then
	!        read (unit=8,fmt=*) dummy
	!        read (unit=8,fmt='(3F12.6)') (dire(i),i=1,noptx)
	!     end if
		 read (unit=8,fmt=*) dummy
		 do i = 1, noptx
			  read (unit=8,fmt='(5F20.6)') (oh(i,j),j=1,i)
		 end do
		 do i = 1, noptx-1
			  oh(i,i+1:noptx) = oh(i+1:noptx,i)
		 end do 
	end if
	close(8)

	! Copy data from single-image arrays into all-image array
	fullxopt(img_num,:)=xopt(:)
	fulloe(img_num)=oe
	fullox(img_num,:)=ox(:)
	fullog(img_num,:)=og(:)
	fulloh(img_num,:,:)=oh(:,:)
	
! End of loop over all images
end do
END SUBROUTINE read_checkfile


