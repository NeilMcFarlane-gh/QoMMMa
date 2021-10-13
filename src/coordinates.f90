MODULE coordinates
use nrtype
implicit none

integer(i4b), parameter :: maxbond = 8
integer(i4b) :: n, nx, nq, nqx, nl, nlx, nopt, noptx, ndlc, nprim, ninact, ncon, kcnstyp, nimg, nebtype, gsmtype, coordtype
real(sp) :: tolde_org, tolgmax_org, tolgrms_org, toldxmax_org, toldxrms_org, tolper

real(sp) :: kspring

logical, allocatable :: modchg(:), inact(:)

real(sp), parameter :: bohr=.529177d0
integer(i4b), allocatable :: nbonds(:),bonds(:,:),attyp(:),links(:,:), &
    &  qm(:),opt(:)

! tables of coordinate of various groups of atoms, of all images 
real(sp), allocatable :: fullx(:,:), fullxq(:,:), fullxl(:,:), fullxopt(:,:), full_dlc(:,:)

! tables of coordinates of various groups of atoms, of a particular image
real(sp), allocatable :: x(:), xq(:), xl(:), xopt(:), dlc(:), lratio(:), chg(:)

! matrices used in the generation of DLC.
integer(i4b), allocatable :: full_prim_list(:,:,:), prim_list(:,:) ! Primitive coordinate indice arrays.
integer(i4b), allocatable :: to_generate(:) ! Used dynamically to assign temporary atom indices.
real(sp), allocatable :: full_prims(:,:), prims(:) ! Primitive internal coordinate arrays.
real(sp), allocatable :: full_Bmat_p(:,:,:), Bmat_p(:,:) ! Primitive Wilson B matrix arrays.
real(sp), allocatable :: full_Bmat_dlc(:,:,:), Bmat_dlc(:,:) ! DLC Wilson B matrix arrays.
real(sp), allocatable :: full_Gmat(:,:,:), Gmat(:,:) ! G matrix arrays.
real(sp), allocatable :: full_Umat(:,:,:), Umat(:,:) ! U matrix arrays.
real(sp), allocatable :: full_Rmat(:,:,:), Rmat(:,:) ! R matrix arrays.

character(4), allocatable :: img_string(:)
character(3), allocatable :: label(:)
character(2), allocatable :: llabel(:), qlabel(:)

END MODULE coordinates


