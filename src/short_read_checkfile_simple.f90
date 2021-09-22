SUBROUTINE short_read_checkfile_simple()
use nrtype ; use coordinates ; use optimdata
implicit none

! This subroutine reads the "CheckFile" briefly, so as to get the atomic charges.
! It can only work after "alloc_coord" - so one needs to check that the arrays are allocated.

character(3) :: tmplab
character(80) :: dummy
integer(i4b) :: rstat
integer(i4b) :: i, j, k, ii, jj, img_num

! First check whether the arrays are in fact already allocated...

open(unit=8,file=("CheckFile"))

read(unit=8,fmt=*) dummy
read(unit=8,fmt=*) dummy
read(unit=8,fmt=*) i
read(unit=8,fmt=*) dummy
read(unit=8,fmt=*) i, j, k, ii

read(unit=8,fmt=*) dummy
! This line actually contains data about images, but processing is not needed
read(unit=8,fmt=*) dummy

! Now number & type of constraints
read(unit=8,fmt=*) dummy
read(unit=8,fmt=*) ncon, kcnstyp

! Now read through list of QM atoms.
read(unit=8,fmt=*) dummy
do i=1,nq
     read(unit=8,fmt=*) opt(i)
end do

! Now read list of Link atoms and their details.

read(unit=8,fmt=*) dummy
! read through links: MM atom; QM atom; rQL/rQM ideal ratio; label.
if (nl.gt.0) then
   do i = 1, nl
        read(unit=8,fmt=*) opt(nq+i)
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
               write (*,*) "Trying to inactivate a QM atom. ERROR."
               close(8)
               stop
          end if
     end do
     do j = 1, nl
         if (links(j,1) .eq. jj) then
               write (*,*) "Trying to inactivate a link atom. ERROR."
               close(8)
               stop
          end if
     end do
     inact(jj) = .true.
end do

! that's all that's needed!

close(8)

END SUBROUTINE short_read_checkfile_simple


