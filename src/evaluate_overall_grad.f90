SUBROUTINE evaluate_overall_grad()
use nrtype ; use coordinates ; use optimdata
implicit none

! This subroutine evaluates the total gradient, including constraint terms, then writes it,
!    the qm and the mm gradients in Hartree/Å.

integer(i4b) :: i, j, k, img_num
real(sp) :: vec(3)

!Loop over all images
do img_num=1,nimg
	! evaluate the total gradient as a sum of the MM (tg) and QM (qg) ones.
	! First zero out the gradient
	g=0.d0

	! Copy energies and gradients from all-image arrays
	tg(:)=fulltg(img_num,:)
	qg(:)=fullqg(img_num,:)
	te=fullte(img_num)
	qe=fullqe(img_num)
	x(:)=fullx(img_num,:)
	! And now calculated gradient and energy
	g=tg+qg
	e=te+qe

	! Now assign the gradient on the Hessian-optimized atoms.

	! Clear optg table first
	optg(:)=0.d0
	 
	do i=1, nopt
		j=3*opt(i)-2
		k=3*i-2
		optg(k:k+2)=g(j:j+2)
	end do

	open(unit=8,file=("gradients"//trim(img_string(img_num))),status="replace")

	write (unit=8,fmt='(A)') "Total Gradient (no constraints):"
	do i=1,n
		 j=3*(i-1)+1
		 write (unit=8,fmt='(A3,1X,3F21.15)') label(i), (g(k),k=j,j+2)
	end do
	write (unit=8,fmt='(A)') "MM Gradient (no constraints):"
	do i=1,n
		 j=3*(i-1)+1
		 write (unit=8,fmt='(A3,1X,3F21.15)') label(i), (tg(k),k=j,j+2)
	end do
	write (unit=8,fmt='(A)') "QM Gradient (no constraints):"
	do i=1,n
		 j=3*(i-1)+1
		 write (unit=8,fmt='(A3,1X,3F21.15)') label(i), (qg(k),k=j,j+2)
	end do

	vec(1)=sum(g(1:nx:3))
	vec(2)=sum(g(2:nx:3))
	vec(3)=sum(g(3:nx:3))
	write (unit=8,fmt='(A)') "Sum of gradients:"
	write (unit=8,fmt='(3F12.6)') vec

	vec(1)=sum(tg(1:nx:3))
	vec(2)=sum(tg(2:nx:3))
	vec(3)=sum(tg(3:nx:3))
	write (unit=8,fmt='(A)') "Sum of MM gradients:"
	write (unit=8,fmt='(3F12.6)') vec

	vec(1)=sum(qg(1:nx:3))
	vec(2)=sum(qg(2:nx:3))
	vec(3)=sum(qg(3:nx:3))
	write (unit=8,fmt='(A)') "Sum of QM gradients:"
	write (unit=8,fmt='(3F12.6)') vec

	close(8)

	! Now, if appropriate, evaluate the constraint effect on energy & gradients.
	if (ncon_cart.gt.0) then
		totcnsen=0.d0
		do i=1,ncon_cart
			call evaluate_constraint(i)
			fullcnsval(img_num,i)=cnsval(i)
			fullcnsg(img_num,i)=cnsg(i)
			fullcnsen(img_num,i)=cnsen(i)
		end do
			totcnsen=sum(cnsen(1:ncon_cart))
			fulltotcnsen(img_num)=totcnsen
			e = e + totcnsen
	else if (ncon_prim.gt.0) then
		! TO-DO : Implement prim/DLC constraint here?
		!placeholder
	end if


	fulloptg(img_num,:)=optg(:)
	fulle(img_num)=e

! End of loop over all images
end do    

return

END SUBROUTINE evaluate_overall_grad


