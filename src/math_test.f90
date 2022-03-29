program math_test
use math ; use nrtype

real(sp) :: distance, theta, phi
real(sp), allocatable :: unit_vec(:), proj_vec(:), inv_arr(:,:), output_mat(:,:)
integer(i4b) :: i
integer(i4b), allocatable :: combinations(:,:)
real(sp) :: vectors(3,3), eigens(3,3), coords_1(3), coords_2(3), arr(2,2), det, eigenvals(3), eigenvecs(3,3)
real(sp) :: coords_3(3), grad_r(3,2), grad_theta(3,3), grad_phi(3,4), array(3,3), coords_4(3)
logical :: orthogonality

! Random coordinate set...
coords_1(1) = 0
coords_2(1) = 0
coords_3(1) = 0.819
coords_4(1) = -0.819

coords_1(2) = 0.7375
coords_2(2) = -0.7375
coords_3(2) = 0.817
coords_4(2) = -0.817

coords_1(3) = -0.0528
coords_2(3) = -0.0528
coords_3(3) = 0.422
coords_4(3) = 0.422


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


! Array for testing eigenvalues...
eigens(1,1) = 3
eigens(2,1) = 2
eigens(3,1) = 4
eigens(1,2) = 2
eigens(2,2) = 0
eigens(3,2) = 2
eigens(1,3) = 4
eigens(2,3) = 2
eigens(3,3) = 3


! Array for testing SVD...
arr(:,:) = 0
arr(1,1) = 1
arr(1,2) = 2
arr(2,1) = 3
arr(2,2) = 4


! Array for testing determinant...
array(1,1) = 2
array(1,2) = -3
array(1,3) = 1
array(2,1) = 2
array(2,2) = 0
array(2,3) = -1
array(3,1) = 1
array(3,2) = 4
array(3,3) = 5


eigenvecs = EVECS(eigens, SIZE(eigens, 1))

!print *, "eigenvectors: ", eigenvecs

eigenvals = EVALS(eigens, SIZE(eigens, 1))

!print *, "eigenvalues: ", eigenvals

det = DETERMINANT(array, SIZE(array, 1))

!print *, "determinant: ", det

inv_arr = SVD_inverse(arr, SIZE(arr, 1), SIZE(arr, 2))

!print *, "(pseudo)inverse: ", inv_arr

orthogonality = is_orthog(vectors, SIZE(vectors, 1), SIZE(vectors, 2))

!print *, "orthogonal? ", orthogonality

distance = atom_distance(coords_1, coords_2)

!print *, "distance: ", distance

grad_r = atom_distance_grad(coords_1, coords_2)

!print *, "grad: ", grad_r

theta = atom_angle(coords_2, coords_1, coords_3)

!print *, "angle: ", theta

grad_theta = atom_angle_grad(coords_1, coords_2, coords_3)

!print *, "grad: ", grad_theta

phi = atom_dihedral(coords_3, coords_1, coords_2, coords_4)

print *, "dihedral: ", phi * 57

grad_phi = atom_dihedral_grad(coords_4, coords_1, coords_2, coords_3)

!print *, "grad: ", grad_phi

unit_vec = unit_vector(coords_1, SIZE(coords_1))

!print *, "unit vector: ", unit_vec

proj_vec = vector_project(coords_1, coords_2, SIZE(coords_1))

!print *, "projected vector: ", proj_vec

end program math_test
      
