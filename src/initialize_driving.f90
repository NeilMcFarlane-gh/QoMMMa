SUBROUTINE initialize_driving
use nrtype ; use coordinates ; use optimdata ; use primitive ; use math
implicit none

! In this subroutine, the driving coordinates which were provided in qommma.in are converted to a format consistent with the primitive definitions.
! This means that the atom indices which were used are reset to start at 1 and are placed in numerical order (as is the case in initialize_prims.).
! If a primitive driving coordinate is provided which is not automatically generated by initialize_prims, then it is added to prim_list.
! This can happen in cases where a bond formation is being studied, as this will not be generated by initialize_prims.

integer(i4b) :: i, j, k, l, m, ii
integer(i4b) :: drive_temp, qm_temp, driving_work(ndriv,4), prim_list_work(nprim,4), to_add
integer(i4b) :: prim_temp1(4), prim_temp2(4)
real(sp) :: drive_coord(3), qm_coord(3)

if (gsmtype .eq. 2) then
	! The first part is relatively simple. The coordinates in xopt are compared to the indices of the QM atoms.
	! After comparison, the numbers can simply be reset to their relative indice starting from 1.
	
	! First, get the QM region coordinates.
	do l=1, nq
		j = (3 * (l-1)) + 1
		k = (3 * (qm(l)-1)) + 1
		xq(j:j+2) = x(k:k+2)
	end do

	! Now, a series of loops and if-statements compares coordinates to generate the reordered driving coordinates.
	driving_work = driving_coords
	do i=1, SIZE(driving_work,1)
		do j=1, 4
			drive_temp = driving_work(i,j)
			if (drive_temp .ne. 0) then
				k = (3 * (drive_temp-1)) + 1
				drive_coord = x(k:k+2)
				do l=1, nq
					qm_temp = qm(l)
					if (drive_temp .eq. qm_temp) then
						do m=1, nq
							ii = (3 * (m-1)) + 1
							qm_coord = xq(ii:ii+2)
							if (MAXVAL(drive_coord - qm_coord) .eq. 0) then
								driving_coords(i,j) = l
							end if
						end do
					end if
				end do
			else 
				cycle
			end if
		end do
	end do	
	
	! If there are any primitive internal coordinates which were not automatically generated, then these are added to prim_list.
	! First the new version of prim_list must be allocated.
	prim_list_work = prim_list
	to_add = 0
	do i=1, SIZE(driving_coords,1)
		prim_temp1 = driving_coords(i,:)
		do j=1, SIZE(prim_list_work,1)
			prim_temp2 = prim_list_work(j,:)
			if (MAXVAL(prim_temp1 - prim_temp2) .eq. 0) then ! the primtive internal coordinate has already been generated so no need to add it.
				exit
			else if (j .eq. SIZE(prim_list_work,1)) then
				to_add = to_add + 1 ! we need to add a new one to prim_list.
			end if
		end do
	end do
	if (to_add .gt. 0) then
		nprim = nprim + to_add
		deallocate(prim_list)
		allocate(prim_list(nprim,4))
	end if
	
	! Now, the new primitives can be added to the array using the same method as above.
	to_add = 0
	do i=1, SIZE(prim_list_work,1)
		prim_list(i,:) = prim_list_work(i,:)
	end do
	do i=1, SIZE(driving_coords,1)
		prim_temp1 = driving_coords(i,:)
		do j=1, SIZE(prim_list_work,1)
			prim_temp2 = prim_list_work(j,:)
			if (MAXVAL(prim_temp1 - prim_temp2) .eq. 0) then
				exit
			else if (j .eq. SIZE(prim_list_work,1)) then
				to_add = to_add + 1
				prim_list((SIZE(prim_list_work,1) + to_add),:) = prim_temp1
			end if
		end do
	end do
	
end if	

end subroutine initialize_driving