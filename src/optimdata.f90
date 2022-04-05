MODULE optimdata
use nrtype
implicit none

integer(i4b) :: nstep, disp
logical :: converged, line_search, reversing, resetting

logical, allocatable :: fullconverged(:), update_geom(:), climbing(:)
character(len=3) :: convs(6)
character(len=3), allocatable :: fullconvs(:,:)
integer(i4b), allocatable :: cnstyp(:), cnsat(:,:), ncnsat(:),cnstyp_dlc(:), cnsat_dlc(:,:), ncnsat_dlc(:)
real(sp) :: te, qe, e, oe, conv(5), totcnsen, weight
real(sp) :: tolde, tolgmax, tolgrms, toldxmax, toldxrms
real(sp), allocatable :: fulle(:), fulloe(:), fullte(:), fullqe(:), fulltotcnsen(:), &
          fullconv(:,:), fullcnsen(:,:), fullcnsg(:,:), fullcnsval(:,:), &
          & fullcnsen_dlc(:,:), fullcnsg_dlc(:,:), fullcnsval_dlc(:,:)

real(sp), allocatable :: tg(:), qg(:), g(:), optg(:), og(:), h(:,:), oh(:,:), &
     &  ox(:), newx(:), mull(:), kcns(:), cnsidl(:), cnsval(:), cnsen(:), cnsg(:), &
	 & cnsidl_dlc(:), cnsval_dlc(:), cnsen_dlc(:), cnsg_dlc(:)
real(sp), allocatable :: fullog(:,:), fulloh(:,:,:), fullox(:,:), fulltg(:,:), &
     &  fullqg(:,:), fullmull(:,:), fulloptg(:,:), fullh(:,:,:), fullnewx(:,:), norm_per_force(:)
	 
! for DLC
real(sp), allocatable :: optg_dlc(:), og_dlc(:), h_dlc(:,:), oh_dlc(:,:), new_dlc(:), fullnew_dlc(:,:)
real(sp), allocatable :: optg_p(:), og_p(:)
real(sp), allocatable :: fulloptg_dlc(:,:), fullog_dlc(:,:), fullh_dlc(:,:,:), fulloh_dlc(:,:,:)
real(sp), allocatable :: h_p(:,:), oh_p(:,:), fullh_p(:,:,:), fulloh_p(:,:,:)

! for mecp
real(sp),allocatable :: qga(:),qgb(:),fullea(:),fulleb(:)
real(sp),allocatable :: fullqga(:,:),fullqgb(:,:),fullqea(:),fullqeb(:)     
real(sp),allocatable :: fulloea(:),fulloeb(:),ga(:),gb(:),optga(:),optgb(:)
real(sp),allocatable :: fulloptga(:,:),fulloptgb(:,:)
real(sp) :: qea,qeb,ea,eb,oea,oeb

real(sp), parameter :: hart_kcal = 627.5095d0
real(sp), parameter :: stpmx = 0.5
integer(i4b), parameter :: maxcnsat_cart = 10
integer(i4b), parameter :: maxcnsat_dlc = 4
! The most complicated constraint will apply to less than 10 atoms
! All the following convergence criteria are read from user input through read_converg
!real(sp), parameter :: tolde_org = 1.d-4
!real(sp), parameter :: tolgmax_org = 3.d-3
!real(sp), parameter :: tolgrms_org = 2.d-3
!real(sp), parameter :: toldxmax_org = 7.d-3
!real(sp), parameter :: toldxrms_org = 4.d-3
!real(sp), parameter :: tolper = 2.d-3

! These parts apply to the dispersion calculation in disp_corr()

real(sp), allocatable :: dispgrad(:)
real(sp) :: edistot
real(sp), allocatable :: fullqend(:)
END MODULE optimdata


