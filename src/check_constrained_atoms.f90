SUBROUTINE check_constrained_atoms()
use nrtype ; use coordinates ; use optimdata
implicit none

! Checks that all constrained atoms are in the Hess-optimized ones.

integer(i4b) :: i, j, k, cnsok

cnsok=1

! Checked only for one image, as all images have the same constraints
do i = 1, ncon
	do j = 1, ncnsat(i)
		cnsok=0
		do k=1,nopt
			if (opt(k).eq.cnsat(i,j)) then
				cnsok=1
			end if
		end do
		if (cnsok.eq.0) then
			write (*,*) "You have tried to constrain a non-Hessopt atom."
			write (*,*) "This is not allowed. ERROR."
			stop
		end if
	end do
end do
return
END SUBROUTINE check_constrained_atoms

