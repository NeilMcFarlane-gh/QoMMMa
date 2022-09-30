MODULE coordinates
use nrtype
implicit none

integer(i4b), parameter :: maxbond = 8
integer(i4b) :: n, nx, nq, nqx, nl, nlx, nopt, noptx, ndlc, nprim, ninact, kcnstyp, nimg, nebtype 
integer(i4b) :: gsmtype, gsmphase, coordtype, primtype, ncon_prim, ncon_cart, disp_prim, add_prims
real(dp) :: tolde_org, tolgmax_org, tolgrms_org, toldxmax_org, toldxrms_org, tolper

real(dp) :: kspring

logical, allocatable :: modchg(:), inact(:)

real(dp), parameter :: bohr=.529177d0
real(dp), parameter :: cut_off = 1.5 ! Angstroms
integer(i4b), allocatable :: nbonds(:),bonds(:,:), bonds_xopt(:,:),attyp(:),links(:,:), &
    &  qm(:),opt(:)

! tables of coordinate of various groups of atoms, of all images 
real(dp), allocatable :: fullx(:,:), fullxq(:,:), fullxl(:,:), fullxopt(:,:)

! tables of coordinates of various groups of atoms, of a particular image
real(dp), allocatable :: x(:), xq(:), xl(:), xopt(:), dlc(:), lratio(:), chg(:), x_copy(:)

! matrices used in the generation and constraining of DLC.
! these matrices are frequently deallocated and reallocated.
integer(i4b), allocatable :: prim_list(:,:), prim_add_list(:,:) ! Primitive coordinate indice array.
integer(i4b), allocatable :: to_generate(:) ! Used dynamically to assign temporary atom indices.
real(dp), allocatable :: prims(:), old_prims(:), prims_save(:) ! Primitive internal coordinate array.
real(dp), allocatable :: Bmat_p(:,:), old_Bmat_p(:,:), Bmat_p_save(:,:) ! Primitive Wilson B matrix array.
real(dp), allocatable :: Bmat_dlc(:,:), old_Bmat_dlc(:,:) ! DLC Wilson B matrix array.
real(dp), allocatable :: Gmat(:,:) ! G matrix array.
real(dp), allocatable :: Umat(:,:) ! U matrix array.
real(dp), allocatable :: Vmat(:,:) ! V matrix array.
real(dp), allocatable :: Rmat(:,:) ! R matrix array.

! for NEB
character(4), allocatable :: img_string(:)
character(3), allocatable :: label(:)
character(2), allocatable :: llabel(:), qlabel(:)

END MODULE coordinates


