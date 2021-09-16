program math_test
use math ; use nrtype

real(sp) :: distance
real(sp), allocatable :: grad(:,:), unit_vec(:), proj_vec(:), inv_arr(:,:)
real(sp) :: vectors(3,3), coords_1(3), coords_2(3), arr(2,2), det
logical :: orthogonality


! Random coordinate set...
coords_1(1) = 2
coords_2(1) = 4
coords_1(2) = 6
coords_2(2) = 1
coords_1(3) = 5
coords_2(3) = 2


! Set of orthogonal vectors...
!Vector 1
vectors(1,1) = 1
vectors(2,1) = 0
vectors(3,1) = -1
!Vector 2
vectors(1,2) = 1
vectors(2,2) = SQRT(2.0)
vectors(3,2) = 1
!Vector 2
vectors(1,3) = 1
vectors(2,3) = SQRT(2.0) * (-1)
vectors(3,3) = 1


! Array for testing SVD...
arr(:,:) = 0
arr(1,1) = 1
arr(1,2) = 2
arr(2,1) = 3
arr(2,2) = 4


det = DETERMINANT(arr, SIZE(arr, 1))

print *, det

inv_arr = SVD_inverse(arr, SIZE(arr, 1), SIZE(arr, 2))

!print *, "(pseudo)inverse: ", inv_arr

orthogonality = is_orthog(vectors, SIZE(vectors, 1), SIZE(vectors, 2))

!print *, "orthogonal? ", orthogonality

distance = atom_distance(coords_1, coords_2)

!print *, "distance: ", distance

grad = atom_distance_grad(coords_1, coords_2)

!print *, "grad: ", grad

unit_vec = unit_vector(coords_1, SIZE(coords_1))

!print *, "unit vector: ", unit_vec

proj_vec = vector_project(coords_1, coords_2, SIZE(coords_1))

!print *, "projected vector: ", proj_vec

end program math_test
      
