MODULE coordinates
use nrtype
implicit none

integer(i4b), parameter :: maxbond = 8
integer(i4b) :: n, nx, nq, nqx, nl, nlx, nopt, noptx, ndlc, nprim, ninact, ncon, kcnstyp, nimg, nebtype 
integer(i4b) :: gsmtype, coordtype, primtype
real(sp) :: tolde_org, tolgmax_org, tolgrms_org, toldxmax_org, toldxrms_org, tolper

real(sp) :: kspring

logical, allocatable :: modchg(:), inact(:)

real(sp), parameter :: bohr=.529177d0
real(sp), parameter :: cut_off = 1.5 ! Angstroms
integer(i4b), allocatable :: nbonds(:),bonds(:,:), bonds_xopt(:,:),attyp(:),links(:,:), &
    &  qm(:),opt(:)

! tables of coordinate of various groups of atoms, of all images 
real(sp), allocatable :: fullx(:,:), fullxq(:,:), fullxl(:,:), fullxopt(:,:), full_dlc(:,:)

! tables of coordinates of various groups of atoms, of a particular image
real(sp), allocatable :: x(:), xq(:), xl(:), xopt(:), dlc(:), lratio(:), chg(:), x_copy(:)

! matrices used in the generation of DLC.
! these matrices are frequently deallocated and reallocated.
integer(i4b), allocatable :: prim_list(:,:) ! Primitive coordinate indice array.
integer(i4b), allocatable :: to_generate(:) ! Used dynamically to assign temporary atom indices.
real(sp), allocatable :: prims(:), old_prims(:), d_prims(:) ! Primitive internal coordinate array.
real(sp), allocatable :: Bmat_p(:,:), old_Bmat_p(:,:) ! Primitive Wilson B matrix array.
real(sp), allocatable :: Bmat_dlc(:,:), old_Bmat_dlc(:,:) ! DLC Wilson B matrix array.
real(sp), allocatable :: Gmat(:,:) ! G matrix array.
real(sp), allocatable :: Umat(:,:) ! U matrix array.
real(sp), allocatable :: Rmat(:,:) ! R matrix array.

character(4), allocatable :: img_string(:)
character(3), allocatable :: label(:)
character(2), allocatable :: llabel(:), qlabel(:)

END MODULE coordinates


