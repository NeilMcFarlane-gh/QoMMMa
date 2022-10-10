MODULE optimdata
use nrtype
implicit none

integer(i4b) :: nstep, disp
logical :: converged, line_search, reversing, resetting

logical, allocatable :: fullconverged(:), update_geom(:), climbing(:)
character(len=3) :: convs(6)
character(len=3), allocatable :: fullconvs(:,:)
integer(i4b), allocatable :: cnstyp(:), cnsat(:,:), ncnsat(:), cnsat_p(:,:,:), cns_n_coeff_p(:)
integer(i4b), allocatable :: dq_at_p(:,:), dq_pos_p(:)
real(dp) :: te, qe, e, oe, conv(5), totcnsen, weight
real(dp) :: tolde, tolgmax, tolgrms, toldxmax, toldxrms
real(dp), allocatable :: fulle(:), fulloe(:), fullte(:), fullqe(:), fulltotcnsen(:), &
          fullconv(:,:), fullcnsen(:,:), fullcnsg(:,:), fullcnsval(:,:), cnscoeff_p(:,:)

real(dp), allocatable :: tg(:), qg(:), g(:), optg(:), og(:), h(:,:), oh(:,:), &
     &  ox(:), newx(:), mull(:), kcns(:), cnsidl(:), cnsval(:), cnsen(:), cnsg(:), &
	 &  cdat(:,:), cdat_unproj(:,:), dq_p(:)
real(dp), allocatable :: fullog(:,:), fulloh(:,:,:), fullox(:,:), fulltg(:,:), &
     &  fullqg(:,:), fullmull(:,:), fulloptg(:,:), fullh(:,:,:), fullnewx(:,:), norm_per_force(:)
	 
! for DLC
real(dp), allocatable :: optg_dlc(:), og_dlc(:), h_dlc(:,:), oh_dlc(:,:), new_dlc(:), fullnew_dlc(:,:)
real(dp), allocatable :: optg_p(:), og_p(:)
real(dp), allocatable :: fulloptg_dlc(:,:), fullog_dlc(:,:), fullh_dlc(:,:,:), fulloh_dlc(:,:,:)
real(dp), allocatable :: h_p(:,:), oh_p(:,:), fullh_p(:,:,:), fulloh_p(:,:,:)

! for mecp
real(dp),allocatable :: qga(:),qgb(:),fullea(:),fulleb(:)
real(dp),allocatable :: fullqga(:,:),fullqgb(:,:),fullqea(:),fullqeb(:)     
real(dp),allocatable :: fulloea(:),fulloeb(:),ga(:),gb(:),optga(:),optgb(:)
real(dp),allocatable :: fulloptga(:,:),fulloptgb(:,:)
real(dp) :: qea,qeb,ea,eb,oea,oeb

real(dp), parameter :: hart_kcal = 627.5095d0
real(dp), parameter :: stpmax_dlc = 0.1
real(dp), parameter :: stpmax_cart = 0.1
integer(i4b), parameter :: maxcnsat_cart = 10 ! maximum constraint of 10 atoms
integer(i4b), parameter :: maxcnsat_dlc = 10  ! maximum constraint of 10 primitive internal coordindates
											  ! In principle, we could constrain any number, but keep it simple.

! These parts apply to the dispersion calculation in disp_corr()

real(dp), allocatable :: dispgrad(:)
real(dp) :: edistot
real(dp), allocatable :: fullqend(:)
END MODULE optimdata


