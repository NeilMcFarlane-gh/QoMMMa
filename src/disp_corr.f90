SUBROUTINE disp_corr()
use coordinates; use nrtype; use optimdata  ! optimdata is where you store things like edistot and dispgrad
implicit none

! This subroutine calculates the dispersion energy correction to the QM energy
! It also calculates dispersion gradient corrections
! Equations and constants taken from Grimme S, J. Comp. Chem., 2006, 27, 1787

real(dp), dimension(nq) :: rzero, csix     ! Each is an array containing the constants for each atom
real(dp), parameter :: ssix = 1.05d0        ! This parameter is specific to B3LYP and BP86
real(dp), parameter :: hart_j_fac = 2625499.621d0      ! Required to get the answer in hartrees
integer(i4b) :: i, j, ii, jj, k, l
real(dp), dimension(nq,nq) :: dist, fdmp, edis, gfac
real(dp), dimension(3) :: vecs
real(dp) :: qc(3*nq)
! real(dp) :: edistot

! -----------------------------------------------------------
! assign values of rzero and csix to arrays rzero() and csix() - more elements can be added
! ----------------------------------------------------------

do i=1,nq
  if (qlabel(i).eq."H ") then
    csix(i) = 0.14d0
    rzero(i) = 1.001d0
  else if (qlabel(i).eq."He") then
    csix(i) = 0.08d0
    rzero(i) = 1.012d0
  else if (qlabel(i).eq."Li") then
    csix(i) = 1.61d0
    rzero(i) = 0.825d0
  else if (qlabel(i).eq."Be") then
    csix(i) = 1.61d0
    rzero(i) = 1.408d0
  else if (qlabel(i).eq."B ") then
    csix(i) = 3.13d0
    rzero(i) = 1.485d0
  else if (qlabel(i).eq."C ") then
    csix(i) = 1.75d0
    rzero(i) = 1.452d0
  else if (qlabel(i).eq."N ") then
    csix(i) = 1.23d0
    rzero(i) = 1.397d0
  else if (qlabel(i).eq."O ") then
    csix(i) = 0.70d0
    rzero(i) = 1.342d0
  else if (qlabel(i).eq."F ") then
    csix(i) = 0.75d0
    rzero(i) = 1.287d0
  else if (qlabel(i).eq."Ne") then
    csix(i) = 0.63d0
    rzero(i) = 1.243d0
  else if (qlabel(i).eq."Na") then
    csix(i) = 5.71d0
    rzero(i) = 1.144d0
  else if (qlabel(i).eq."Mg") then
    csix(i) = 5.71d0
    rzero(i) = 1.364d0
  else if (qlabel(i).eq."Al") then
    csix(i) = 10.79d0
    rzero(i) = 1.639d0
  else if (qlabel(i).eq."Si") then
    csix(i) = 9.23d0
    rzero(i) = 1.716d0
  else if (qlabel(i).eq."P ") then
    csix(i) = 7.84d0
    rzero(i) = 1.705d0
  else if (qlabel(i).eq."S ") then
    csix(i) = 5.57d0
    rzero(i) = 1.683d0
  else if (qlabel(i).eq."Cl") then
    csix(i) = 5.07d0
    rzero(i) = 1.639d0
  else if (qlabel(i).eq."Ar") then
    csix(i) = 4.61d0
    rzero(i) = 1.595d0
  else if (qlabel(i).eq."K ") then
    csix(i) = 10.80d0
    rzero(i) = 1.485d0
  else if (qlabel(i).eq."Ca") then
    csix(i) = 10.80d0
    rzero(i) = 1.474d0
  else if (qlabel(i).eq."Sc") then
    csix(i) = 10.80d0
    rzero(i) = 1.562d0
  else if (qlabel(i).eq."Ti") then
    csix(i) = 10.80d0
    rzero(i) = 1.562d0
  else if (qlabel(i).eq."V ") then
    csix(i) = 10.80d0
    rzero(i) = 1.562d0
  else if (qlabel(i).eq."Cr") then
    csix(i) = 10.80d0
    rzero(i) = 1.562d0
  else if (qlabel(i).eq."Mn") then
    csix(i) = 10.80d0
    rzero(i) = 1.562d0
  else if (qlabel(i).eq."Fe") then
    csix(i) = 10.80d0
    rzero(i) = 1.562d0
  else if (qlabel(i).eq."Co") then
    csix(i) = 10.80d0
    rzero(i) = 1.562d0
  else if (qlabel(i).eq."Ni") then
    csix(i) = 10.80d0
    rzero(i) = 1.562d0
  else if (qlabel(i).eq."Cu") then
    csix(i) = 10.80d0
    rzero(i) = 1.562d0
  else if (qlabel(i).eq."Zn") then
    csix(i) = 10.80d0
    rzero(i) = 1.562d0
  else if (qlabel(i).eq."Ga") then
    csix(i) = 16.99d0
    rzero(i) = 1.650d0
  else if (qlabel(i).eq."Ge") then
    csix(i) = 17.10d0
    rzero(i) = 1.727d0
  else if (qlabel(i).eq."As") then
    csix(i) = 16.37d0
    rzero(i) = 1.760d0
  else if (qlabel(i).eq."Se") then
    csix(i) = 12.64d0
    rzero(i) = 1.771d0
  else if (qlabel(i).eq."Br") then
    csix(i) = 12.47d0
    rzero(i) = 1.749d0
  else if (qlabel(i).eq."Kr") then
    csix(i) = 12.01d0
    rzero(i) = 1.727d0
  else if (qlabel(i).eq."Rb") then
    csix(i) = 24.67d0
    rzero(i) = 1.628d0
  else if (qlabel(i).eq."Sr") then
    csix(i) = 24.67d0
    rzero(i) = 1.606d0
  else if (qlabel(i).eq."Y ") then
    csix(i) = 24.67d0
    rzero(i) = 1.639d0
  else if (qlabel(i).eq."Zr") then
    csix(i) = 24.67d0
    rzero(i) = 1.639d0
  else if (qlabel(i).eq."Nb") then
    csix(i) = 24.67d0
    rzero(i) = 1.639d0
  else if (qlabel(i).eq."Mo") then
    csix(i) = 24.67d0
    rzero(i) = 1.639d0
  else if (qlabel(i).eq."Tc") then
    csix(i) = 24.67d0
    rzero(i) = 1.639d0
  else if (qlabel(i).eq."Ru") then
    csix(i) = 24.67d0
    rzero(i) = 1.639d0
  else if (qlabel(i).eq."Rh") then
    csix(i) = 24.67d0
    rzero(i) = 1.639d0
  else if (qlabel(i).eq."Pd") then
    csix(i) = 24.67d0
    rzero(i) = 1.639d0
  else if (qlabel(i).eq."Ag") then
    csix(i) = 24.67d0
    rzero(i) = 1.639d0
  else if (qlabel(i).eq."Cd") then
    csix(i) = 24.67d0
    rzero(i) = 1.639d0
  else if (qlabel(i).eq."In") then
    csix(i) = 37.32d0
    rzero(i) = 1.672d0
  else if (qlabel(i).eq."Sn") then
    csix(i) = 38.71d0
    rzero(i) = 1.804d0
  else if (qlabel(i).eq."Sb") then
    csix(i) = 38.44d0
    rzero(i) = 1.881d0
  else if (qlabel(i).eq."Te") then
    csix(i) = 31.74d0
    rzero(i) = 1.892d0
  else if (qlabel(i).eq."I ") then
    csix(i) = 31.50d0
    rzero(i) = 1.892d0
  else if (qlabel(i).eq."Xe") then
    csix(i) = 29.99d0
    rzero(i) = 1.881d0
  else
    write(*,*) "Dispersion parameters not available for atom type", qlabel(i)
    STOP
  end if
! write (*,'(A,I5,A,F12.6)') "In disp_corr, I have just set csix for atom",i,"to ",csix(i)
end do

! Convert the csix parameters to Hartrees Angstrom**6
csix = csix * 1.d6 / hart_j_fac

! write (*,*) "QM Coordinates:"

do i=1,nq
   j=3*(i-1)+1
   k=3*(qm(i)-1)+1
   qc(j:j+2) = x(k:k+2) 
!   write (*,'(A3,3F12.6)') qlabel(i),qc(j:j+2)
end do

! write (*,'(A,I5)') "Number of quantum atoms, nq= ", nq

! ------------------------------
! calculate everything in one go 
! ------------------------------

dispgrad=0.
do i=1,nq-1
  do j=i+1,nq
    ii = 3*(i-1) + 1
    jj = 3*(j-1) + 1
    k = 3*(qm(i)-1) + 1
    l = 3*(qm(j)-1) + 1
    vecs = qc(ii:ii+2) - qc(jj:jj+2)
!    write (*,'(A,2I5,A,3F12.6)') "Treating coords starting at:" ,ii, jj,"that is:", vecs
    dist(i,j) = sqrt(sum((vecs)**2))
!    write (*,'(A,I5,A,I5,A,F12.6)') "distance between atom i= ",i," and j= ",j," is: ",dist(i,j)
    fdmp(i,j) = 1 / (1 + exp(-20.d0 * (dist(i,j) / (rzero(i) + rzero(j)) - 1.d0)))
!    write (*,'(A,F14.10)') "damping function", fdmp(i,j)
    edis(i,j) = (-1.d0 * ssix * (sqrt(csix(i) * csix(j))) * fdmp(i,j) / (dist(i,j)**6))
!    write (*,'(A,F12.4)') "dispersion correction in: ",edis(i,j)
    gfac(i,j) = -6.d0 * edis(i,j) / dist(i,j) + &
        & (edis(i,j) * fdmp(i,j) * 20.d0 * exp(-20.d0 * (dist(i,j) / (rzero(i) + rzero(j)) - 1.d0))) & 
        & / (rzero(i) + rzero(j))
!    write (*,'(A,F18.12)') "gradient factor", gfac(i,j)
    dispgrad(k:k+2) = dispgrad(k:k+2) + gfac(i,j)*vecs/dist(i,j)   
!    write (*,*) "dispersion gradient on atom", i, "due to atom", j
!    write (*,*) dispgrad(ii:ii+2)
    dispgrad(l:l+2) = dispgrad(l:l+2) - gfac(i,j)*vecs/dist(i,j)
!    write (*,*) "dispersion gradient on atom", j, "due to atom", i
!    write (*,*) dispgrad(jj:jj+2)
  end do
end do

edistot = SUM(edis)

write (*,'(A,F18.12,A)') "Dispersion energy is: ", edistot, " Hartree"
! write (*,'(A)') "Dispersion gradient (Hartrees/Ang):"
! write (*,*) dispgrad

END SUBROUTINE disp_corr
