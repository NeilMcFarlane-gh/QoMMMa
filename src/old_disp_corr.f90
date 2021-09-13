SUBROUTINE disp_corr
use coordinates; use nrtype
implicit none

! This subroutine calculates the dispersion energy correction to the QM energy
! It also calculates dispersion gradient corrections
! Equations and constants taken from Grimme S, J. Comp. Chem., 2006, 27, 1787
! The dispersion energy will be SUBTRACTED from the total QM energy

real dimension(nq) :: rzero, csix ! Each is an array containing the constants for each atom
real parameter :: ssix = 1.05    ! This parameter is specific to B3LYP
real parameter :: d = 20
integer(i4b) :: a, b
real dimension(nq,nq) :: dist, fdmp, edis
real dimension(3) :: vecs
real :: edistot

! -----------------------------------------------------------
! assign values of rzero and csix to arrays rzero() and csix() - more elements can be added
! -----------------------------------------------------------

do a=1,nq
  if (qlabel(a).eq."H ") then
    csix(a) = 0.14
    rzero(a) = 1.001
  else if (qlabel(a).eq."C ") then
    csix(a) = 1.75
    rzero(a) = 1.452
  else if (qlabel(a).eq."N ") then
    csix(a) = 1.23
    rzero(a) = 1.397
  else if (qlabel(a).eq."O ") then
    csix(a) = 0.70
    rzero(a) = 1.342
  else if (qlabel(a).eq."S ") then
    csix(a) = 5.57
    rzero(a) = 1.683
  else if (qlabel(a).eq."Fe") then
    csix(a) = 10.80
    rzero(a) = 1.562
  else
    write(*,*) "Dispersion parameters not available for atom type", qlabel(a)
  end if
end do

! ------------------------------------------------------------------
! create a distance matrix dist - diagonal elements will be zero a=b
! ------------------------------------------------------------------

do a=1,nq
  do b=1,nq
  if (a.ne.b) then
    j = 3*(a-1) + 1
    k = 3*(b-1) + 1
    vecs = xq(j:j+2) - xq(k:k+2)
    dist(a,b) = sqrt((sum(vecs))**2)
  else
    dist(a,b) = 0
  end if
  end do
end do

! ----------------------------------------------
! calculate a matrix of the damping forces, fdmp
! ----------------------------------------------

do a=1,nq                                
  do b=1,nq                              
  if (a.ne.b) then                       
    fdmp(a,b) = 1/1+exp(-d*(dist(a,b)/(rzero(a)+rzero(b))))
  else                                   
    fdmp(a,b) = 0                        
  end if
  end do
end do

! --------------------------------------------
! calculate the dispersion energy matrix, edis
! --------------------------------------------

do a=1,nq
  do b=1,nq
  if (a.ne.b) then
   edis(a,b) = s6*sqrt(csix(a)*csix(b))*fdmp(a,b)/(dist(a,b)**6)
  else
   edis(a,b)
  end if
  end do
end do

! -------------------------------------------------------------------------------
! calculate the dispersion energy, edistot - divide by 2 to avoid double-counting
! -------------------------------------------------------------------------------

edistot = SUM(edis)/2

! -------------------------------------------------------------------
! calculate a matrix of the gradients - 3 dimensions (x,y,z) per atom
! -------------------------------------------------------------------






END SUBROUTINE disp_corr
