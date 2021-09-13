SUBROUTINE evaluate_overall_grad_mecp()
use nrtype ; use coordinates ; use optimdata
implicit none

! This subroutine evaluates the total gradient, including constraint terms, then writes it,
!    the qm and the mm gradients in Hartree/Å.
! evaluate the total gradient as a sum of the MM (tg) and QM (qg) ones.

integer(i4b) :: i, j, k, img_num
real(sp) :: vec(3)

do img_num=1,nimg
! First zero out the gradient
ga=0.d0
gb=0.d0

! Copy energies and gradients from all-image arrays
tg(:)=fulltg(img_num,:)
qga(:)=fullqga(img_num,:)
qgb(:)=fullqgb(img_num,:)
te=fullte(img_num)
qea=fullqea(img_num)
qeb=fullqeb(img_num)
x(:)=fullx(img_num,:)
! And now calculated gradient and energy

ga=tg+qga
gb=tg+qgb
ea=te+qea
eb=te+qeb

! Now assign the gradient on the Hessian-optimized atoms.

! Clear optg table first
optga(:)=0.d0
optgb(:)=0.d0
 
do i=1, nopt
    j=3*opt(i)-2
    k=3*i-2
    optga(k:k+2)=ga(j:j+2)
    optgb(k:k+2)=gb(j:j+2)
end do

open(unit=8,file=("gradients"//trim(img_string(img_num))),status="replace")

write (unit=8,fmt='(A)') "Total Gradient (no constraints) of state A:"
do i=1,n
     j=3*(i-1)+1
     write (unit=8,fmt='(A3,1X,3F21.15)') label(i), (ga(k),k=j,j+2)
end do
write (unit=8,fmt='(A)') "Total Gradient (no constraints) of state B:"
do i=1,n
     j=3*(i-1)+1
     write (unit=8,fmt='(A3,1X,3F21.15)') label(i), (gb(k),k=j,j+2)
end do
write (unit=8,fmt='(A)') "MM Gradient (no constraints):"
do i=1,n
     j=3*(i-1)+1
     write (unit=8,fmt='(A3,1X,3F21.15)') label(i), (tg(k),k=j,j+2)
end do
write (unit=8,fmt='(A)') "QM Gradient (no constraints) of state A:"
do i=1,n
     j=3*(i-1)+1
     write (unit=8,fmt='(A3,1X,3F21.15)') label(i), (qga(k),k=j,j+2)
end do
write (unit=8,fmt='(A)') "QM Gradient (no constraints) of state B:"
do i=1,n
     j=3*(i-1)+1
     write (unit=8,fmt='(A3,1X,3F21.15)') label(i), (qgb(k),k=j,j+2)
end do

vec(1)=sum(ga(1:nx:3))
vec(2)=sum(ga(2:nx:3))
vec(3)=sum(ga(3:nx:3))
write (unit=8,fmt='(A)') "Sum of gradients of state A:"
write (unit=8,fmt='(3F12.6)') vec

vec(1)=sum(gb(1:nx:3))
vec(2)=sum(gb(2:nx:3))
vec(3)=sum(gb(3:nx:3))
write (unit=8,fmt='(A)') "Sum of gradients of state B:"
write (unit=8,fmt='(3F12.6)') vec

vec(1)=sum(tg(1:nx:3))
vec(2)=sum(tg(2:nx:3))
vec(3)=sum(tg(3:nx:3))
write (unit=8,fmt='(A)') "Sum of MM gradients:"
write (unit=8,fmt='(3F12.6)') vec

vec(1)=sum(qga(1:nx:3))
vec(2)=sum(qga(2:nx:3))
vec(3)=sum(qga(3:nx:3))
write (unit=8,fmt='(A)') "Sum of QM gradients of state A:"
write (unit=8,fmt='(3F12.6)') vec

vec(1)=sum(qgb(1:nx:3))
vec(2)=sum(qgb(2:nx:3))
vec(3)=sum(qgb(3:nx:3))
write (unit=8,fmt='(A)') "Sum of QM gradients of state B:"
write (unit=8,fmt='(3F12.6)') vec

close(8)

! Now, if appropriate, evaluate the constraint effect on energy & gradients.

totcnsen=0.d0
if (ncon.ne.0) then

do i=1,ncon
    call evaluate_constraint(i)
    fullcnsval(img_num,i)=cnsval(i)
    fullcnsg(img_num,i)=cnsg(i)
    fullcnsen(img_num,i)=cnsen(i)
end do
totcnsen=sum(cnsen(1:ncon))

fulltotcnsen(img_num)=totcnsen
ea = ea + totcnsen
eb = eb + totcnsen

end if

fulloptga(img_num,:)=optga(:)
fulloptgb(img_num,:)=optgb(:)
fullea(img_num)=ea
fulleb(img_num)=eb

! End of loop over all images
end do    

return

END SUBROUTINE evaluate_overall_grad_mecp


