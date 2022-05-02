MODULE optimdata
use nrtype
implicit none

integer(i4b) :: nstep, disp
logical :: converged, line_search, reversing, resetting

logical, allocatable :: fullconverged(:), update_geom(:), climbing(:)
character(len=3) :: convs(6)
character(len=3), allocatable :: fullconvs(:,:)
integer(i4b), allocatable :: cnstyp(:), cnsat(:,:), ncnsat(:), cnsat_p(:,:)
real(sp) :: te, qe, e, oe, conv(5), totcnsen, weight
real(sp) :: tolde, tolgmax, tolgrms, toldxmax, toldxrms
real(sp), allocatable :: fulle(:), fulloe(:), fullte(:), fullqe(:), fulltotcnsen(:), &
          fullconv(:,:), fullcnsen(:,:), fullcnsg(:,:), fullcnsval(:,:)

real(sp), allocatable :: tg(:), qg(:), g(:), optg(:), og(:), h(:,:), oh(:,:), &
     &  ox(:), newx(:), mull(:), kcns(:), cnsidl(:), cnsval(:), cnsen(:), cnsg(:), &
	 & cnsdq_p(:), cnspos_p(:)
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
real(sp), parameter :: stpmax_dlc = 0.5
real(sp), parameter :: stpmax_cart = 0.1
integer(i4b), parameter :: maxcnsat_cart = 10
integer(i4b), parameter :: maxcnsat_dlc = 4

! These parts apply to the dispersion calculation in disp_corr()

real(sp), allocatable :: dispgrad(:)
real(sp) :: edistot
real(sp), allocatable :: fullqend(:)
END MODULE optimdata


