subroutine update_opt_geometry()
use primitive ; use delocalised ; use math ; use nrtype ; use coordinates ; use optimdata 
implicit none


! Specially Adapted BFGS routine from Numerical Recipes

integer(i4b) :: I,j, img_num, k
real(dp) :: DelG(noptx), HDelG(noptx), ChgeX(noptx), ChgeS(ndlc), DelX(noptx), w(noptx)
real(dp) :: fac, fad, fae, sumdg, sumdx, stpl, lstep, stpmax, stpmx, maxchgx, temp_x(noptx), temp_dlc
real(dp),parameter ::  eps = 1.d-6  ! eps = 1.d-5 ! eps = 3.d-8 !
real(dp),parameter :: minstep = 0.1
logical :: is_cons

line_search=.false.
if ((nebtype.eq.4).or.(nebtype.eq.3).or.(nebtype.eq.5)) then
    call update_conj_opt()
else
    update_geom=.true.


if (coordtype .eq. 0) then	
	do img_num=1,nimg
		! Initialise the maximum step.
		stpmax = stpmax_cart
		stpmx = stpmax_cart * REAL(noptx,sp)
		
		! Copy data first
		oe=fulloe(img_num)
		e=fulle(img_num)
		optg(:)=fulloptg(img_num,:)
		og(:)=fullog(img_num,:)
		xopt(:)=fullxopt(img_num,:)
		ox(:)=fullox(img_num,:)
		oh(:,:)=fulloh(img_num,:,:)

		if (nebtype.eq.0) then
			IF (Nstep.eq.0) THEN
				ChgeX = -0.7d0*optg
				h = oh
			else
				DelG = optg - og
				DelX = xopt - ox	! xi in NR
				HDelG = MATMUL(oh,DelG)
				fac = 1.d0 / SUM(DelG * DelX)	! NR has a check to see this is not too small
				fae = SUM(DelG * HDelG)
				fad = 1.d0 / fae
				sumdg = SUM(DelG**2)
				sumdx = SUM(DelX**2)
				w = fac * DelX - fad * HDelG	! see Schlegel p. 473; called dg in NR
				DO I = 1, noptx
					h(i,i:noptx) = oh(i,i:noptx) + fac * DelX(I) * DelX(I:noptx) - &
					& fad * HDelG(I) * HDelG(I:Noptx) + fae * w(I) * w(I:noptx)
					h(i:noptx,i) = H(I,I:noptx)
				END DO
				ChgeX = - MATMUL(h,optg)
			end if
		else
			if (img_num.eq.1) then
				WRITE (*,*)  'First image fixed'
				ChgeX = 0.0d0
				h = oh
				fulloe(img_num)=e
				update_geom(img_num)=.false.
			else if (img_num.eq.nimg) then
				WRITE (*,*)  'Last image fixed'
				ChgeX = 0.0d0
				h = oh
				fulloe(img_num)=e
				update_geom(img_num)=.false.
			else
				if (Nstep.eq.0) then
					ChgeX = -0.7d0*optg
					h = oh
				else
					DelG = optg - og
					DelX = xopt - ox	! xi in NR
					HDelG = MATMUL(oh,DelG)
					fac = SUM(DelG * DelX)	
					sumdg = SUM(DelG**2)
					sumdx = SUM(DelX**2)
		! Test fac - skip update if not sufficiently positive - from NR
		!            if (fac**2.gt.(eps * sumdg * sumdx)) then
		!            if (abs(fac).gt.eps) then
						fac = 1.d0 / fac
						fae = SUM(DelG * HDelG)
						fad = 1.d0 / fae
						w = fac * DelX - fad * HDelG	! see Schlegel p. 473; called dg in NR
						DO I = 1, noptx
							h(i,i:noptx) = oh(i,i:noptx) + fac * DelX(I) * DelX(I:noptx) - &
							& fad * HDelG(I) * HDelG(I:Noptx) + fae * w(I) * w(I:noptx)
							h(i:noptx,i) = H(I,I:noptx)
						END DO
						ChgeX = - MATMUL(h,optg)
						maxchgx=maxval(abs(ChgeX))
						if (maxval(abs(ChgeX)).lt.(0.25*toldxmax_org)) then
							ChgeX = 0.0d0
							update_geom(img_num)=.false.
							h = oh
						end if
		!           else
		!           
		!                ChgeX = 0.0d0
		!                update_geom(img_num)=.false.
		!            end if
				end if
			end if 
		end if

		stpl = SQRT(SUM(ChgeX**2))
		IF (stpl .gt. stpmax) THEN
			ChgeX = ChgeX / stpl * stpmax
			write (*,*) "Changing (1) Step Length"
		END IF

		lstep=maxval(abs(chgex))
		IF (lstep .gt. STPMX) THEN
			ChgeX = ChgeX / lstep * STPMX
			write (*,*)"Changing (2) Step Length"
		END IF
		newx = xopt + ChgeX

		! evaluate convergence tests

		convs = " NO"
		converged=.false.
		i=0
		conv(1)=e-oe
		conv(2)=maxval(abs(ChgeX))
		conv(3)=sqrt(sum(chgex**2)/real(noptx,sp))
		conv(4)=maxval(abs(optg))
		conv(5)=sqrt(sum(optg**2)/real(noptx,sp))

		! Tighter convergence test for climbing image
		if (climbing(img_num)) then
			if(abs(conv(1)) .lt. tolde_org) then
				convs(1)="YES" ; i=i+1
			end if
			if (conv(2) .lt. toldxmax_org) then
				convs(2)="YES" ; i=i+1
			end if
			if (conv(3) .lt. toldxrms_org) then
				convs(3)="YES" ; i=i+1
			end if
			if (conv(4) .lt. tolgmax_org) then
				convs(4)="YES" ; i=i+1
			end if
			if (conv(5) .lt. tolgrms) then
				convs(5)="YES" ; i=i+1
			end if
		else
			if (abs(conv(1)) .lt. tolde) then
				 convs(1)="YES" ; i=i+1
			end if
			if (conv(2) .lt. toldxmax) then
				 convs(2)="YES" ; i=i+1
			end if
			if (conv(3) .lt. toldxrms) then
				 convs(3)="YES" ; i=i+1
			end if
			if (conv(4) .lt. tolgmax) then
				 convs(4)="YES" ; i=i+1
			end if
			if (conv(5) .lt. tolgrms) then
				 convs(5)="YES" ; i=i+1
			end if
		end if

		if (i .eq. 5) converged = .true.

		fullnewx(img_num,:)=newx(:)
		fullh(img_num,:,:)=h(:,:)
		fullconverged(img_num)=converged
		fullconv(img_num,:)=conv(:)
		fullconvs(img_num,:)=convs(:)

		open(unit=8,file=("add_to_update"//trim(img_string(img_num))),status="replace")
		!write (unit=8,fmt='(A,4F9.6,3F12.6)') " DelX,Delg, DelX*Delg, maxChgX, fac, fad, fae: ", sum(DelX), sum(DelG), sum(DelX*DelG), maxval(abs(ChgeX)), fac, fad, fae
		write (unit=8,fmt='(A,3F9.6, F13.6)') " DelX, Delg, DelX*Delg, maxChgX: ", sum(DelX), sum(DelG), sum(DelX*DelG), maxchgx
		close(8)

	!End of loop over images
	end do
	
else if (coordtype .eq. 1) then
	do img_num=1,nimg
		! Initialise the maximum step.
		stpmax = 0.05
		stpmx = 0.05 * REAL(noptx,sp)
		
		! Copy data first
		oe=fulloe(img_num)
		e=fulle(img_num)
		optg(:)=fulloptg(img_num,:)
		og(:)=fullog(img_num,:)
		xopt(:)=fullxopt(img_num,:)
		ox(:)=fullox(img_num,:)
		oh(:,:)=fulloh(img_num,:,:)
		x_copy = xopt

		! The primitives, DLC and the B matrices are maintained...
		call maintain_DLC(nopt, ndlc, xopt, cdat)

		! A selection of values are preserved as they are used in the BFGS update.
		if (.not. ALLOCATED(old_prims)) allocate(old_prims(nprim))
		old_prims(:) = prims(:)
		if (.not. ALLOCATED(old_Bmat_dlc)) allocate(old_Bmat_dlc(ndlc, (3 * nopt)))
		old_Bmat_dlc(:,:) = Bmat_dlc(:,:)
		if (.not. ALLOCATED(old_Bmat_p)) allocate(old_Bmat_p(nprim, (nopt * 3)))
		old_Bmat_p(:,:) = Bmat_p(:,:)

		! Now, the BFGS algorithm can be used to generate the change in DLC from the calculated gradient.
		! Firstly, the gradients from both the current and previous step must be updated to DLC subspace.
		! In addition, the gradients are updated to primitive subspace as these are used in the updating of the hessian matrix.
		call gen_grad_cart_to_DLC(nopt, ndlc, nprim, Bmat_dlc, optg, optg_dlc)
		call gen_grad_cart_to_DLC(nopt, ndlc, nprim, old_Bmat_dlc, og, og_dlc)
		call gen_grad_cart_to_prim(nopt, nprim, Bmat_p, optg, optg_p)
		call gen_grad_cart_to_prim(nopt, nprim, old_Bmat_p, og, og_p)

		! To avoid the repeated matrix diagonalisation that would be necessary to continually update the hessian matrix in DLC subspace, we update the primitive hessian.
		! However, on the first optimisation step, the initial matrix is generated simply as a weighted unit matrix.
		! Subsequent iterations of the optimisation will use the BFGS scheme to update this primitive hessian matrix.
		if (Nstep .eq. 0) then
			allocate(h_p(nprim, nprim))
			allocate(oh_p(nprim, nprim))
			oh_p(:,:) = 0.0
			h_p(:,:) = 0.0
			do i=1, nprim
				do j=1, nprim
					if (i == j) then
						oh_p(i,j) = 1
						h_p(i,j) = 1
					end if
				end do
			end do
		else
			allocate(h_p(nprim, nprim))
			h_p(:,:) = oh_p(:,:)
		end if

		! Checking if there are any constraints...
		is_cons = .False.
		if (allocated(cdat) .eqv. .True.) then
			is_cons = .True.
		end if

		if (Nstep .le. 100) then
			!###################
			! Steepest descent.#
			!###################
			
			! Scaled by 0.7 to avoid particularly large steps.
			ChgeS = optg_dlc * 0.7 * (-1)

			! Now, if there are any constraints, given elements of the DLC should not change.
			if (is_cons .eqv. .True.) then
				ChgeS(ndlc - (ncon_prim - 1):ndlc) = 0.0
			end if
			
			! The change is DLC is scaled using the maximum step length.
			stpl = RMSD_CALC(ChgeS, dlc, ndlc)
			IF (stpl .gt. stpmax) THEN
				ChgeS = ChgeS / stpl * stpmax
				write (*,*) "Changing (1) Step Length"
			END IF
			lstep = maxval(abs(ChgeS))
			IF (lstep .gt. STPMX) THEN
				ChgeS = ChgeS / lstep * STPMX
				write (*,*)"Changing (2) Step Length"
			END IF

			!The new DLC and, more importantly, cartesian coordinates can now be evaluated.
			temp_x(:) = xopt(:)
			call DLC_to_cart(nopt, ndlc, nprim, ChgeS, dlc, xopt, newx)
			call calc_prims(nopt, nprim, prims, opt, newx, prim_list)
			ChgeX = newx(:) - temp_x(:)
		else
			!###############
			! Quasi newton.#
			!###############
			
			! The hessian matrix is updated to DLC subspace.
			call gen_hess_prim_to_DLC(nopt, ndlc, nprim, Umat, h_p, h_dlc)

			! The change in DLC can now be evaluated using a quasi-Newton methodology.
			ChgeS = MATMUL(TRANSPOSE(h_dlc), optg_dlc) * (-1)
			
		    ! Now, if there are any constraints, given elements of the DLC should not change.
			if (is_cons .eqv. .True.) then
				ChgeS(ndlc - (ncon_prim - 1):ndlc) = 0.0
			end if
		
			! The change is DLC is scaled using the maximum step length.
			stpl = RMSD_CALC(ChgeS, dlc, ndlc)
			IF (stpl .gt. stpmax) THEN
				ChgeS = ChgeS / stpl * stpmax
				write (*,*) "Changing (1) Step Length"
			END IF
			lstep = maxval(abs(ChgeS))
			IF (lstep .gt. STPMX) THEN
				ChgeS = ChgeS / lstep * STPMX
				write (*,*)"Changing (2) Step Length"
			END IF
			
			! The new DLC and, more importantly, cartesian coordinates can now be evaluated.
			temp_x(:) = xopt(:)
			call DLC_to_cart(nopt, ndlc, nprim, ChgeS, dlc, xopt, newx)
			ChgeX = newx(:) - temp_x(:)
			
			! Lastly, using the BFGS method, the primitive hessian for the next optimisation cycle is calculated.
			! To ensure continuity, a new set of primitive internal coordinates is calculated for newx.
			call calc_prims(nopt, nprim, prims, opt, newx, prim_list)
			call update_bfgs_p(nopt, nprim, h_p, optg_p, og_p, prims, old_prims)
		end if
		
		! evaluate convergence tests
		convs = " NO"
		converged=.false.
		i=0
		conv(1)=e-oe
		conv(2)=maxval(abs(ChgeX))
		conv(3)=sqrt(sum(chgex**2)/real(noptx,sp))
		conv(4)=maxval(abs(optg_dlc(1:(ndlc - ncon_prim))))
		conv(5)=sqrt(sum(optg_dlc(1:(ndlc - ncon_prim))**2)/real((ndlc - ncon_prim),sp))

		! Tighter convergence test for climbing image
		if (climbing(img_num)) then
			if(abs(conv(1)) .lt. tolde_org) then
				convs(1)="YES" ; i=i+1
			end if
			if (conv(2) .lt. toldxmax_org) then
				convs(2)="YES" ; i=i+1
			end if
			if (conv(3) .lt. toldxrms_org) then
				convs(3)="YES" ; i=i+1
			end if
			if (conv(4) .lt. tolgmax_org) then
				convs(4)="YES" ; i=i+1
			end if
			if (conv(5) .lt. tolgrms) then
				convs(5)="YES" ; i=i+1
			end if
		! When in the growth phase of the string, the convergence criteria are loosened by a factor of 10.
		else if (gsmphase .eq. 1) then
			if (abs(conv(1)) .lt. (tolde * 10)) then
				 convs(1)="YES" ; i=i+1
			end if
			if (conv(2) .lt. (toldxmax * 10)) then
				 convs(2)="YES" ; i=i+1
			end if
			if (conv(3) .lt. (toldxrms * 10)) then
				 convs(3)="YES" ; i=i+1
			end if
			if (conv(4) .lt. (tolgmax * 10)) then
				 convs(4)="YES" ; i=i+1
			end if
			if (conv(5) .lt. (tolgrms * 10)) then
				 convs(5)="YES" ; i=i+1
			end if
		else
			if (abs(conv(1)) .lt. tolde) then
				 convs(1)="YES" ; i=i+1
			end if
			if (conv(2) .lt. toldxmax) then
				 convs(2)="YES" ; i=i+1
			end if
			if (conv(3) .lt. toldxrms) then
				 convs(3)="YES" ; i=i+1
			end if
			if (conv(4) .lt. tolgmax) then
				 convs(4)="YES" ; i=i+1
			end if
			if (conv(5) .lt. tolgrms) then
				 convs(5)="YES" ; i=i+1
			end if
		end if

		if (i .eq. 5) converged = .true.

		fullnewx(img_num,:)=newx(:)
		fullh(img_num,:,:)=h(:,:)
		fullconverged(img_num)=converged
		fullconv(img_num,:)=conv(:)
		fullconvs(img_num,:)=convs(:)

		open(unit=8,file=("add_to_update"//trim(img_string(img_num))),status="replace")
		write (unit=8,fmt='(A,3F9.6, F13.6)') " DelX, Delg, DelX*Delg, maxChgX: ", sum(DelX), sum(DelG), sum(DelX*DelG), maxchgx
		close(8)

	!End of loop over images
	end do
end if	

! increment nstep
nstep=nstep + 1

end if
return
END SUBROUTINE update_opt_geometry
