SUBROUTINE evaluate_constraint(c)
use nrtype ; use coordinates ; use optimdata
implicit none

! This subroutine evaluates the constraint term for constraint c,
!    including the contribution to the gradient
!    This latter term is then added to the overall gradient.

integer(i4b), parameter :: maxcns=maxcnsat/2
integer(i4b), intent(in) :: c
integer(i4b) :: i, j, k, kk, l(maxcnsat)
real(sp) :: vecs(maxcns,3), cg(maxcns,3), lvec(maxcns), kcnsau

! Reminder of the constraint types:
! 1 = r(A-B)
! 2 = r(A-B) - r(C-D)
! 3 = r(A-B) + r(C-D)
! 4 = r(A-B) + r(C-D) - r(E-F))

kcnsau=kcns(c)/Hart_kcal
! How many atoms are affected by this constraint?
kk=ncnsat(c)/2

! Evaluate the relevant vectors and their lengths
do i=1,kk
    ! integer position in the FULL geometry array of the two atoms
    l(2*i-1) = 3*cnsat(c,2*i-1)-2
    l(2*i) = 3*cnsat(c,2*i)-2
    vecs(i,:) = x(l(2*i):l(2*i)+2) - x(l(2*i-1):l(2*i-1)+2)
    lvec(i)=sqrt(sum(vecs(i,:)**2))
end do

! Then evaluate the constraint value.

if (cnstyp(c).eq.1) then
    cnsval(c)=lvec(1)
else if (cnstyp(c).eq.2) then
    cnsval(c)=lvec(1)-lvec(2)
else if (cnstyp(c).eq.3) then
    cnsval(c)=lvec(1)+lvec(2)
else if (cnstyp(c).eq.4) then
    cnsval(c)=lvec(1)+lvec(2)-lvec(3)
end if

! Now compute the corresponding energy contribution, and scalar gradient contribution

call calc_constraint(cnsen(c),cnsg(c),kcnstyp,cnsval(c)-cnsidl(c),kcnsau)

! And now work out the corresponding gradient contribution!

cg=0._sp
if (cnstyp(c).eq.1) then
    cg(1,:) = vecs(1,:)/lvec(1) * cnsg(c)
else if (cnstyp(c).eq.2) then
    cg(1,:) = vecs(1,:)/lvec(1) * cnsg(c)
    cg(2,:) = -vecs(2,:)/lvec(2) *cnsg(c)
else if (cnstyp(c).eq.3) then
    cg(1,:) = vecs(1,:)/lvec(1) * cnsg(c)
    cg(2,:) = vecs(2,:)/lvec(2) *cnsg(c)
else if (cnstyp(c).eq.4) then
    cg(1,:) = vecs(1,:)/lvec(1) * cnsg(c)
    cg(2,:) = vecs(2,:)/lvec(2) *cnsg(c)
    cg(3,:) = -vecs(3,:)/lvec(3) *cnsg(c)
end if

! and spread it out on the atoms in the array of gradients for Hess-optimized atoms

! write (*,*) "nopt is",nopt,"and the opt are",(opt(j),j=1,nopt)
! write (*,*) "Before"
! do i=1,nopt
! write (*,'(A3,3F13.6)') label(opt(i)),(optg(i*3-2:i*3))
! end do
do i=1,kk
    do j=1, nopt
        if (opt(j).eq.cnsat(c,2*i-1)) then
             k=3*j-2
             optg(k:k+2)=optg(k:k+2)-cg(i,:)
! write (*,*) "Constraint on ",label(opt(j))," number ",opt(j)," opt. number ",j
        else if (opt(j).eq.cnsat(c,2*i)) then
             k=3*j-2
             optg(k:k+2)=optg(k:k+2)+cg(i,:)
! write (*,*) "Constraint on ",label(opt(j))," number ",opt(j)," opt. number ",j
        end if
    end do
end do

! write (*,*) "After"
! do i=1,nopt
! write (*,'(A3,3F13.6)') label(opt(i)),(optg(i*3-2:i*3))
! end do
! rem opt(i) stores the number, in the full geometry, of the i-th optimized atom

return

END SUBROUTINE evaluate_constraint


SUBROUTINE calc_constraint(en,grad,typ,x,k)
use nrtype
implicit none

! This subroutine evaluates the energy (in au) and (scalar) gradient
! for a displacement x, a force constant k (in au/Ang/Ang), for either a
! harmonic constraint (kcnstyp=1), or a tanh-based one (kcnstyp=2)

! The tanh function is chosen in the following way:
!  The gradient for small x is k * x
!  The gradient plateaus at 0.2 au/angstrom. (tanhmax option)

real(sp), intent(out) :: en, grad
integer(i4b),intent(in) :: typ
real(sp), intent(in) :: x, k
real(sp), parameter :: tanhmax=.2_sp

if (typ.eq.1) then
! easy! Harmonic
     en=.5_sp*k*x**2
     grad=k*x
else if (typ.eq.2) then
     grad=tanhmax * tanh(k/tanhmax*x)
     en= tanhmax**2/k * log(cosh(k/tanhmax*x))
else
     write (*,*) "wrong k_constraint type. ERROR."
     stop
end if

return
END SUBROUTINE calc_constraint
