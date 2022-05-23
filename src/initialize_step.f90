SUBROUTINE initialize_step
use nrtype ; use coordinates ; use optimdata ; use primitive ; use math
implicit none

! In this subroutine, the step(s) associated with the constrained primitive internal coordinates given in qommma.in is/are taken.
! The useful thing is that within the Fortran code, both the double- and single-ended variants of GSM are the same.
! All the generation of primitive internal coordinate constraints is handled within the Python code.
! Overall, this leans up the Fortran code, and allows for easier manipulation of GSM-specific code.
!
! Practically speaking, the primitive internal coordinate constraints which were provided in qommma.in are converted to a format consistent with the primitive definitions.
! This means that the atom indices which were used are reset to start at 1 and are placed in numerical order (as is the case in initialize_prims).
! If a primitive driving coordinate (for SE-GSM) is provided which is not automatically generated by initialize_prims, then it is added to prim_list.
! This can happen in cases where a bond formation is being studied, as this will not be generated by initialize_prims.


integer(i4b) :: i, j, k, l, m, ii, in_pos, old_nprim
integer(i4b) :: cons_temp, qm_temp, cons_work(ncon_prim,4), prim_list_work(nprim,4), to_add
integer(i4b) :: prim_temp1(4), prim_temp2(4)
real(sp) :: temp_coord(3), qm_coord(3), dq
real(sp), allocatable :: complete_dq(:)

if (ncon_prim .gt. 0) then
	! The first part is relatively simple. The coordinates in xopt are compared to the indices of the QM atoms.
	! After comparison, the numbers can simply be reset to their relative indice starting from 1.
	
	! First, get the QM region coordinates.
	xq = 0.0
	do l=1, nq
		j = (3 * (l-1)) + 1
		k = (3 * (qm(l)-1)) + 1
		xq(j:j+2) = x(k:k+2)
	end do
	
	! Now, get the optimised region coordinates (includes link atoms).
	xopt = 0.d0
	do i=1, nq
		j = (3 * (i-1)) + 1
		k = (3 * (opt(i)-1)) + 1
		xopt(j:j+2) = x(k:k+2)
	end do
	do i=1, nl
		j = (3 * nq) + i
		xopt(j:j+2) = xl(i:i+2)
	end do

	! Now, a series of loops and if-statements compares coordinates to generate the reordered coordinate steps.
	cons_work(:,:) = cnsat_p(:,:)
	do i=1, ncon_prim
		do j=1, 4
			cons_temp = cons_work(i,j)
			if (cons_temp .ne. 0) then
				k = (3 * (cons_temp-1)) + 1
				temp_coord = x(k:k+2)
				do l=1, nq
					qm_temp = qm(l)
					if (cons_temp .eq. qm_temp) then
						do m=1, nq
							ii = (3 * (m-1)) + 1
							qm_coord = xq(ii:ii+2)
							if (MAXVAL(temp_coord - qm_coord) .eq. 0) then
								cons_work(i,j) = l
							end if
						end do
					end if
				end do
			else 
				cycle
			end if
		end do
	end do
	
	! Updating the constrained atom array with the correct integers.
	cnsat_p(:,:) = cons_work(:,:)
	
	! If there are any primitive internal coordinates which were not automatically generated, then these are added to prim_list.
	! First the new version of prim_list must be allocated.
	prim_list_work = prim_list
	old_nprim = nprim
	to_add = 0
	do i=1, ncon_prim
		prim_temp1 = cons_work(i,:)
		do j=1, nprim
			prim_temp2 = prim_list_work(j,:)
			if (MAXVAL(prim_temp1 - prim_temp2) .eq. 0) then ! the primtive internal coordinate has already been generated so no need to add it.
				exit
			else if (j .eq. nprim) then
				to_add = to_add + 1 ! we need to add a new one to prim_list.
			end if
		end do
	end do
	if (to_add .gt. 0) then
		nprim = nprim + to_add
		deallocate(prim_list)
		allocate(prim_list(nprim,4))
	end if
	
	! Now, the new primitives can be added to the array prim_list using the same method as above.
	to_add = 0
	do i=1, old_nprim
		prim_list(i,:) = prim_list_work(i,:)
	end do
	do i=1, ncon_prim
		prim_temp1 = cons_work(i,:)
		do j=1, old_nprim
			prim_temp2 = prim_list_work(j,:)
			if (MAXVAL(prim_temp1 - prim_temp2) .eq. 0) then
				! At this point, the position which the constrained coordinate is at is recorded so that the step can be recorded.
				in_pos = j
				cnspos_p(i) = in_pos
				exit
			else if (j .eq. old_nprim) then
				to_add = to_add + 1
				prim_list((old_nprim + to_add),:) = prim_temp1
		
				! At this point, the position which the constrained coordinate is at is recorded so that the step can be recorded.
				in_pos = old_nprim + to_add
				cnspos_p(i) = in_pos
			end if
		end do
	end do
	
	! A new array which is the complete list of primitive coordinates which change, including all of those which are zero, is created.
	allocate(complete_dq(nprim))
	complete_dq(:) = 0.0
	do i=1, ncon_prim
		j = cnspos_p(i)
		complete_dq(j) = cnsdq_p(i)
	end do
	
	! Lastly, the primitive internal coordinates have been integrated into prim_list, so the step(s) can be taken.
	! If the constrained coordinate is to be kept fixed at its current value, then the value of dq is zero and nothing changes.
	call calc_prims(nopt, nprim, prims, xopt, prim_list)
	call prims_to_cart(nopt, nprim, complete_dq, prims, xopt, newx, Bmat_p, prim_list)

	! Remember to copy over the new cartesian coordinates to the cartesian region!
	do l=1, nq
		j = (3 * (l-1)) + 1
		k = (3 * (qm(l)-1)) + 1
		x(k:k+2) = newx(j:j+2)
	end do
	do i=1,nq
		j=3*(i-1)+1
		k=3*(qm(i)-1)+1
		xq(j:j+2)=x(k:k+2)
	end do
	if (nl .ne. 0) then
	    do i=1,nl
		    j=3*(i-1)+1
		    k=(nq*3) + (3*(i-1)+1)
		    xl(j:j+2)=newx(k:k+2)
	    end do
	end if 

end if	

end subroutine initialize_step