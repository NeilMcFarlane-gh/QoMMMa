SUBROUTINE Effective_Gradient_mecp()
use nrtype ; use coordinates ; use optimdata
implicit none
       
! Computes the parallel and perpendicular compenents of the Effective Gradient
! As well as the effective gradient itself.
       
integer i,img_num

double precision :: ParG(nx),PerpG(nx),npg,pp
double precision facPP, facP
parameter (facPP=140.d0,facP=1.d0)

!  These factors are only really important for the first step
!  The "difference gradient" is typically ca. 0.075 Hartree/Bohr.
!      i.e. 0.14 Hartree/Angstrom.
!  Assuming this is constant, this gives a Hessian term for the func (Ea-Eb)**2
!     of ca. 0.01 Hartree**2 / Angstrom**2  (0.14**2 / 2)
!  The Hessian term along "normal" coordinates is empirically about 1.4 Hartree / Angstrom**2
!  Using facPP ~ 140 /Hartree means that along this coordinate, too, the Hessian is about right.

optg=0.0d0
do img_num=1,nimg
   optga(:)=fulloptga(img_num,:)
   optgb(:)=fulloptgb(img_num,:)
   ea=fullea(img_num)
   eb=fulleb(img_num)
   npg = 0.d0
   pp = 0.d0
   do i = 1, nq
      PerpG(i) = optga(i) - optgb(i)
      npg = npg + PerpG(i)**2
      pp = pp + optga(i) * PerpG(i)
   end do
   npg = sqrt(npg)
   pp = pp / npg
   do i = 1, nq
      ParG(i) = optga(i) - PerpG(i) / npg * pp
      optg(i) = (ea - eb) * facPP * PerpG(i) + facP * ParG(i)
   end do

! effective gradient is copied as optg in all-image table

   fulloptg(img_num,:)=optg(:)

! End of loop over images
end do

return
END SUBROUTINE Effective_Gradient_mecp
       
