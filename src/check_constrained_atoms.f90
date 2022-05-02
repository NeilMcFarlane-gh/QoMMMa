SUBROUTINE check_constrained_atoms()
use nrtype ; use coordinates ; use optimdata
implicit none

! Checks that all constrained atoms are in the Hess-optimized ones.

integer(i4b) :: i, j, k, cnsok

cnsok=1
! TO-DO : Make it work for prims
! Checked only for one image, as all images have the same constraints

if (coordtype .eq. 0) then
	do i = 1, ncon_cart
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
else if (coordtype .eq. 1) then
	!do i = 1, ncon_prim
		!do j = 1, 4
			!cnsok=0
			!do k=1,nopt
				!print *, cnsat_dlc(i,j), cnsok
			    !if (cnsat_dlc(i,j) .ne. 0) then
					!if (opt(k).eq.cnsat_dlc(i,j)) then
						!cnsok=1
					!end if
				!end if
			!end do
			!if (cnsok.eq.0) then
				!write (*,*) "You have tried to constrain a non-Hessopt atom."
				!write (*,*) "This is not allowed. ERROR."
				!stop
			!end if
		!end do
	!end do
end if

return
END SUBROUTINE check_constrained_atoms

