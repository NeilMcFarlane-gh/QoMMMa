SUBROUTINE write_gradcorrection()
use nrtype ; use coordinates ; use optimdata
implicit none

! This subroutine computes then writes the correction to the qm/mm gradient
!   on the mm atoms, in kcal/mol/Å.

integer(i4b) :: i, j, k, kk, img_num
real(dp) :: vec(3), mulgrad(3), gcorr(3), rmsgcorr, maxgcorr, qxq, lvec
! real(dp) :: mg(nx),cg(nx)

!Loop over all images
do img_num=1,nimg

	! Copies coordinates into single-image array
	x(:)=fullx(img_num,:)
	qg(:)=fullqg(img_num,:)
	! and Mulliken charges
	mull(:)=fullmull(img_num,:)

	open(unit=8,file=("gradcorrection"//trim(img_string(img_num))),status="replace")

	write (unit=8,fmt='(I6)') n

	! If ALL the atoms are in the Hessian Optimization, there is no gradient correction term.
	if (n.eq.nopt) return

	! Otherwise, ...
	rmsgcorr=0.d0
	maxgcorr=0.d0
	outer: do i=1,n
		 ! First check if the atom is included in the Hessian Optimization.
		 ! If yes, the correction gradient is irrelevant and should be set to 0.
		 kk=3*(i-1)+1
		 do j = 1, nopt
			 if (i.eq.opt(j)) then
				  write (unit=8,fmt='(6F16.8)') (x(k),k=kk,kk+2),0.d0,0.d0,0.d0
				  cycle outer
			 end if
		 end do
		 ! Then check if the atom is one of the MM inactive ones - no need to have any gradient
		 !  for these.
		 if (inact(i)) then
			  write (unit=8,fmt='(6F16.8)') (x(k),k=kk,kk+2),0.d0,0.d0,0.d0
			  cycle outer
		 end if
		 ! Otherwise, compute the predicted gradient as a sum over the QM atoms.
		 mulgrad = 0.d0
		 do j = 1, nq
			 vec = (x(3*i-2:3*i)-x(3*qm(j)-2:3*qm(j)))/bohr
			 lvec = sqrt(sum(vec**2))
			 qxq = chg(i)*mull(j)
			 mulgrad = mulgrad - qxq * vec / lvec**3
		 end do
		 
		 ! and convert from Hartree/Bohr to Hartree/Angstrom
		 mulgrad = mulgrad / bohr
		 gcorr = qg(i*3-2:i*3)-mulgrad
		 
		 ! now convert the resulting gradient correction to kcal/mol/Angstrom used by Tinker.
		 gcorr = gcorr * hart_kcal
		 lvec = sum(gcorr**2)
		 if (lvec .gt. maxgcorr) maxgcorr=lvec
		 rmsgcorr = rmsgcorr + lvec
		 write (unit=8,fmt='(6F16.8)') (x(k),k=kk,kk+2),(gcorr(k),k=1,3)
	end do outer

	! report ave. and max magnitude of the correction term:
	rmsgcorr=sqrt(rmsgcorr/(n-nopt))
	maxgcorr=sqrt(maxgcorr)
	write (unit=8,fmt='(A)') "RMS and Max. value of gradient correction term:"
	write (unit=8,fmt='(2F12.6)') rmsgcorr, maxgcorr

	write (8,*) ""

	close(8)

! End of loop over all images
end do
return

END SUBROUTINE write_gradcorrection


