SUBROUTINE read_initfile()
use nrtype ; use coordinates ; use optimdata
implicit none

! This subroutine reads the "fortinput" to get details needed for QM/MM job
! It can only work after "alloc_coord" - so one needs to check that the arrays are allocated.

character(3) :: tmplab
character(80) :: dummy
integer(i4b) :: rstat, prima, gsma, gsmb, coorda
integer(i4b) :: q, i, j, k, ii, jj, kk, iii, nmodif, neba, nebb, nebc
real(dp) :: rql, rqm, tmpchg

! First check whether the arrays are in fact already allocated...
IF (.not.(allocated(lratio).and.allocated(links))) THEN
      write (*,*) "Arrays not yet allocated in read_initfile. ERROR."
      stop
END IF

open(unit=8,file="fortinput",position="rewind")
read(unit=8,fmt=*) dummy
read(unit=8,fmt=*) i, j, k, ii

if ((i.ne.n).or.(j.ne.nq).or.(k.ne.nl).or.(ii.ne.nopt)) THEN
      write (*,*) "Mismatch in number of atoms in input file. ERROR. Check your input or fortinput file."
      close(8)
      stop
end if

! check to see if dispersion is to be calculated
read(unit=8,fmt=*) dummy
read(unit=8,fmt=*) disp

read(unit=8,fmt=*) dummy
read(unit=8,fmt=*) i, j

read(unit=8,fmt=*) dummy
read(unit=8,fmt=*) i

! NEB
read(unit=8,fmt=*) dummy
read(unit=8,fmt=*) neba, nebb, nebc

! GSM
read(unit=8,fmt=*) dummy
read(unit=8,fmt=*) gsma

read(unit=8,fmt=*) dummy
read(unit=8,fmt=*) gsmb

! Coordinate selection
read(unit=8,fmt=*) dummy
read(unit=8,fmt=*) coorda

! If DLC are used, the primitive internal coordinate type
read(unit=8,fmt=*) dummy
read(unit=8,fmt=*) prima

! Now, if they have been pre-defined, read in the primitive internal coordinate definitions.
read(unit=8,fmt=*) dummy
read(unit=8,fmt=*) nprim
if (nprim .gt. 0) then
	allocate(prim_list(nprim,4))
	do q = 1, nprim
		read (unit=8,fmt=*) (prim_list(q,j),j=1,4)
	end do	
end if	

! now reads the QM atoms
read(unit=8,fmt=*) dummy
do i=1,nq
     read(unit=8,fmt=*,IOSTAT=rstat) j, qlabel(i)
     if (rstat .ne. 0) then
         write (*,*) "Error; Check your input or fortinput file in jobfiles directory."
         close(8)
         stop
     end if
     if (j .gt. n) then
         write (*,*) "impossible number for qm atom in input file. Check your input or fortinput file", j
         close(8)
         stop
     end if
     qm(i) = j
     opt(i) = j
end do

! Now read the list of Link atoms
! read links: MM atom; QM atom; rQM ideal, rQL ideal; label.
if (nl.gt.0) then
! link_details will be used only in frequency calculations in qomfreq.py
 open(unit=9,file="link_details",position="rewind")
 read(unit=8,fmt=*) dummy
 do i = 1, nl
     read(unit=8,fmt=*) (links(i,j),j=1,2), rqm, rql, llabel(i)
     lratio(i) = rql/rqm
     k = 1
     do j = 1,nq
         if (qm(j).eq.links(i,2)) then
             k = 0
             write(9,*) j,lratio(i)
             exit
         end if
     end do
     if (k.ne.0) then
         write (*,*) "Link atom not connected to a QM atom. Check your input or fortinput file."
         close(8)
         close(9)
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
         write (*,*) "Link atom is a QM atom. Check your input or fortinput file ",qm(j)
         close(8)
         close(9)
         stop
     end if
     opt(nq+i)=links(i,1)
 end do
 close(9)
end if

! Now read in the list (if any!) of atoms to be included in the QM optimization.
jj=nopt-nq-nl
if (jj.gt.0) then
     read(unit=8,fmt=*) dummy
     do i=1,jj
          read(unit=8,fmt=*) kk
          k = 0
          do ii = 1,nq
              if (qm(ii).eq.kk) then
                  k = 1
                  exit
              end if
          end do
          if (k.ne.0) then
              write (*,*) "You have specified a HessOpt atom which is a QM atom. Check your input or fortinput file",kk
              close(8)
              stop
          end if
          k = 0
          do ii = 1,nl
              if (links(ii,1).eq.kk) then
                  k = 1
                  exit
              end if
          end do
          if (k.ne.0) then
              write (*,*) "You have specified a HessOpt atom which is a Link atom. Check your input or fortinput file", kk
              close(8)
              stop
          end if
          opt(nq+nl+i)=kk
      end do
end if

! Now read in details of constraints (cartesian or primitive) to apply to the geometry
if (ncon_cart.gt.0) then
    read(unit=8,fmt=*) dummy
	do i=1,ncon_cart
		read(unit=8,fmt=*) cnstyp(i),kcns(i),cnsidl(i)
		if ((cnstyp(i).lt.1).or.(cnstyp(i).gt.4)) then
			write (*,*) "Impossible Bond Constraint - Should be between 1-4. Check your input."
			close(8)
			stop
		end if
		if (kcns(i).lt.0.) then
			write (*,*) "Constraint constant kcns negative. Check your input."
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
else if (ncon_prim .gt. 0) then
	read(unit=8,fmt=*) dummy
	do i=1,ncon_prim
		read(unit=8,fmt=*) cnsdq_p(i)
		read(unit=8,fmt=*) (cnsat_p(i,j),j=1,4)
	end do
	! Check that all constrained atoms are QM, Link or HessOpt
	call check_constrained_atoms()
end if

! read in list of atoms on which the charge needs to be changed.
! not allowed to do this for a QM atom!! 
modchg=.false.
read(unit=8,fmt=*) dummy
read(unit=8,fmt=*) nmodif
if (nmodif.gt.0) then
 do i = 1,nmodif
     read (unit=8,fmt=*) jj,tmpchg
     do j = 1, nq
          if (qm(j) .eq. jj) then
               write (*,*) "Trying to set charge on QM atom. Check your input or fortinput file", qm(j)
               close(8)
               stop
          end if
     end do
     chg(jj)=tmpchg
     modchg(jj) = .true.
 end do
end if

! Read in the list of MM atoms which are frozen.
! Not allowed for QM or MM atoms, of course!
read(unit=8,fmt=*) dummy
read(unit=8,fmt=*) ninact
inact = .false.
if (ninact.gt.0) then
 do i = 1, ninact
     read (unit=8,fmt=*) jj
     do j = 1, nq
         if (qm(j) .eq. jj) then
               write (*,*) "ERROR. Check your input or fortinput file. Trying to inactivate a QM atom :", qm(j)
               close(8)
               stop
          end if
     end do
     do j = 1, nl
         if (links(j,1) .eq. jj) then
              write (*,*) "ERROR. Check your input or fortinput file.  Trying to inactivate a link atom :",links(j,1)
              close(8)
              stop
          end if
     end do
     inact(jj) = .true.
 end do
end if
close(8)

! set qm charges to zero.
do i = 1,nq
     chg(qm(i)) = 0.
     modchg(qm(i)) = .true.
end do

return

END SUBROUTINE read_initfile
