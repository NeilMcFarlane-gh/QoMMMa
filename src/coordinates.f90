MODULE coordinates
use nrtype
implicit none

integer(i4b), parameter :: maxbond = 8
integer(i4b) :: n, nx, nq, nqx, nl, nlx, nopt, noptx, ninact, ncon, kcnstyp, nimg, nebtype, gsmtype, coordtype
real(sp) :: tolde_org, tolgmax_org, tolgrms_org, toldxmax_org, toldxrms_org, tolper

real(sp) :: kspring

logical, allocatable :: modchg(:), inact(:)

real(sp), parameter :: bohr=.529177d0
integer(i4b), allocatable :: nbonds(:),bonds(:,:),attyp(:),links(:,:), &
    &  qm(:),opt(:)

! tables of coordinate of various groups of atoms, of all images 
real(sp), allocatable :: fullx(:,:), fullxq(:,:), fullxl(:,:), fullxopt(:,:), fullx_dlc(:,:), fullxq_dlc(:,:)

! tables of coordinates of various groups of atoms, of a particular image
real(sp), allocatable :: x(:), xq(:), x_dlc(:), xq_dlc(:), xl(:), xopt(:), lratio(:), chg(:)

character(4), allocatable :: img_string(:)
character(3), allocatable :: label(:)
character(2), allocatable :: llabel(:), qlabel(:)

END MODULE coordinates


