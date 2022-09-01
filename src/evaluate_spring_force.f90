SUBROUTINE evaluate_spring_force()
use nrtype ; use coordinates ; use optimdata
implicit none

integer(i4b) img_num, max_img
real(dp), allocatable :: tangent(:,:), par_spring(:), per_force(:), par_force(:)
real(dp) :: delta_min, delta_max, max_en

allocate(tangent(nimg,noptx))
allocate(per_force(noptx))
allocate(par_force(noptx))
allocate(par_spring(noptx))

! This subroutine evaluates parallel part of spring force acting between images
! and perpendicular part of the true force acting on images
! The total force is equal to parallel spring part minus perpendicular true part
! And gradient is negative total force

kspring=kspring/Hart_kcal

! Clear logical table pointing on climbing image
climbing(:)=.false.

! If only one geometry, don't evaluate this contribution
if (nimg.eq.1) then
    return
end if

! If Climbing Image NEB specified, find he highest-energy one
if ((nebtype.eq.2).or.(nebtype.eq.5)) then
    max_en=fulle(1)
    max_img=0
    do img_num=2,nimg-1
        if (fulle(img_num).gt. max_en) then
            max_en=fulle(img_num)
            max_img=img_num
        end if
    end do
    ! The highest-energy image marked for climbing
    climbing(max_img)=.true.
    WRITE (*,*)  'Image', max_img, ' is climbing'
end if

tangent(:,:)=0.d0
par_spring(:)=0.d0
do img_num=1,nimg
! Calculate tangent first
    ! For first no calculation
    if (img_num.eq.1) cycle
    ! For last one likewise
    if (img_num.eq.nimg) cycle
    ! Expression for tangent_plus (going uphill)
    if ((fulle(img_num-1).lt.fulle(img_num)).and.(fulle(img_num).lt.fulle(img_num+1))) then
        tangent(img_num,:)=fullxopt(img_num+1,:)-fullxopt(img_num,:)
    ! Expression for tangent_minus (going downhill)
    else if((fulle(img_num-1).gt.fulle(img_num)).and.(fulle(img_num).gt.fulle(img_num+1))) then
        tangent(img_num,:)=fullxopt(img_num,:)-fullxopt(img_num-1,:)
    ! Extremum
    else if(fulle(img_num-1).lt.fulle(img_num+1)) then
        call evaluate_delta_v(img_num,delta_max,delta_min)
        tangent(img_num,:)=(fullxopt(img_num+1,:)-fullxopt(img_num,:))*delta_max + &
                (fullxopt(img_num,:)-fullxopt(img_num-1,:))*delta_min
    ! Extremum
    else if(fulle(img_num-1).gt.fulle(img_num+1) )then
        call evaluate_delta_v(img_num,delta_max,delta_min)
        tangent(img_num,:)=(fullxopt(img_num+1,:)-fullxopt(img_num,:))*delta_min + &
             (fullxopt(img_num,:)-fullxopt(img_num-1,:))*delta_max
    end if

! Normalize tangent
tangent(img_num,:)=tangent(img_num,:)/sqrt(sum((tangent(img_num,:)**2)))

if (climbing(img_num)) then
    
    ! Real force acting on an image is a sum of perpendicular and parallel parts 
    ! true_force = parallel_true_force + perpendicular__true_force
    
    ! In climbing image NEB, spring force does not act on highest energy image
    ! This is the full force due to the potential with the component along the elastic band inverted
    ! so total_force = perpendicular_true force + (-parallel_true_force)
    !    total_force = true_force - parallel_true_force - parallel_true_force
    !    total_force = true_force - 2*parallel_true_force
    
    ! The image moves up the potential energy surface along the elastic band
    ! and down the potential surface perpendicular to the band
    ! Should calculate perpendicular part before gradient is adjusted
    par_force(:) = dot_product(fulloptg(img_num,:), tangent(img_num,:))*tangent(img_num,:)
    per_force(:) = fulloptg(img_num,:) - par_force(:)
    norm_per_force(img_num) = sqrt(sum(per_force(:)**2/real(noptx,sp)))  ! RMS of perpendicular part of the true force
    fulloptg(img_num,:) = per_force(:) - par_force(:)
else
    ! Evaluate parallel part of spring force
    par_spring(:)=kspring*(sqrt(sum((fullxopt(img_num+1,:)-fullxopt(img_num,:))**2)) - &
          sqrt(sum((fullxopt(img_num,:)-fullxopt(img_num-1,:))**2)))*tangent(img_num,:)

    ! Evaluate perpendicular part of true force

    ! In nudged elastic band, total force acting on image is a sum of perpendicular part of true force vector 
    ! and parallel part of spring force vector: 
    ! total_force = parallel_spring_force + perpendicular_true_force

    ! per_force variable is perpendicular_part_of_true_force
    per_force(:) = -fulloptg(img_num,:) + (dot_product(fulloptg(img_num,:), tangent(img_num,:)))*tangent(img_num,:)
    norm_per_force(img_num) = sqrt(sum(per_force(:)**2/real(noptx,sp)))  ! RMS of perpendicular part of the true force
    ! Force = -gradient, so total_gradient= -perpendicular_true_force - parallel spring_force
    fulloptg(img_num,:)=-per_force(:) - par_spring(:)
end if
end do
END SUBROUTINE evaluate_spring_force

SUBROUTINE evaluate_delta_v(c, dmax, dmin)
use nrtype ; use coordinates ; use optimdata
implicit none

integer(i4b)  c

real(dp) :: dmax, dmin
if (abs(fulle(c+1)-fulle(c)).gt.abs(fulle(c-1)-fulle(c))) then
    dmax=abs(fulle(c+1)-fulle(c))
else
    dmax=abs(fulle(c)-fulle(c-1))
end if

if (abs(fulle(c+1)-fulle(c)).lt.abs(fulle(c-1)-fulle(c))) then
    dmin=abs(fulle(c+1)-fulle(c))
else
    dmin=abs(fulle(c)-fulle(c-1))
end if

END SUBROUTINE evaluate_delta_v

