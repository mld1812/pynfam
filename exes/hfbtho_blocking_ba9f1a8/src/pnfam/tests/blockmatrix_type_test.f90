!------------------------------------------------------------------------------
! blockmatrix_type_test.f90
!
! Test code for the module blockmatrix_type.
!
! M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
!------------------------------------------------------------------------------
program blockmatrix_test
   use blockmatrix_type
   implicit none
   type(blockmatrix) :: A, B, C, M, U
   real, parameter :: tolerance = 0.00001
   
   ! Define the grid for the test matrices
   nb = 3
   allocate(db(nb))
   db = [3, 2, 1]
   
   ! Matrix A:
   !  1  2  2  0  0  0
   !  1  0  3  0  0  0
   ! -1  2  1  0  0  0
   !  0  0  0  5  1  0
   !  0  0  0  6  2  0
   !  0  0  0  0  0 20
   call allocate_blockmatrix(A, 14)
   A%ir2c = [1, 2, 3] ; A%ic2r = [1, 2, 3]
   A%ir2m = [1, 10, 14]; A%ic2m = [1, 10, 14]
   A%elem = [1,1,-1, 2,0,2, 2,3,1, 5,6, 1,2, 20]
   
   ! Allocate a matrix M large enough to hold any of the multiplication results
   call allocate_blockmatrix(M, 14)
   
   ! Matrix B:
   !  0  0  0  0  0  6
   !  0  0  0  0  0  2
   !  0  0  0  0  0 -1
   !  3  2  1  0  0  0
   ! -1 -2  0  0  0  0
   !  0  0  0  0  0  0
   call allocate_blockmatrix(B, 9)
   B%ir2c = [3, 1, 0] ; B%ic2r = [2, 0, 1]
   B%ir2m = [1, 4, 0]; B%ic2m = [4, 0, 1]
   B%elem = [6,2,-1, 3,-1, 2,-2, 1,0]
   
   ! Matrix C:
   !  2  3  1  0  0  0
   !  1 -1  0  0  0  0
   !  0  0  2  0  0  0
   !  0  0  0  1  0  0
   !  0  0  0 -1  1  0
   !  0  0  0  0  0  3
   call allocate_blockmatrix(C, 14)
   C%ir2c = [1, 2, 3] ; C%ic2r = [1, 2, 3]
   C%ir2m = [1, 10, 14]; C%ic2m = [1, 10, 14]
   C%elem = [2,1,0, 3,-1,0, 1,0,2, 1,-1, 0,1, 3]
   
   ! Identity matrix U
   call allocate_blockmatrix(U, 14)
   U%ir2c = [1, 2, 3] ; U%ic2r = [1, 2, 3]
   U%ir2m = [1, 10, 14]; U%ic2m = [1, 10, 14]
   U%elem = [1,0,0, 0,1,0, 0,0,1, 1,0, 0,1, 1]
   
   ! Test 1: no transposes, compare result to a hand-computed one
   call triprod('N', A, 'N', B, 'N', C, 1d0, 0d0, M)
   if (all(M%elem - [24,9,-9, 36,40, 34,40, 24,28, 0,0,0,0,0] < tolerance) &
      .and. equal_structure(B,M)) then
      write(*,*) "Test 1 passed"
   else
      write(*,*) "Test 1 FAILED"
   end if
   
   ! Test 2: transpose B only
   call triprod('N', U, 'T', B, 'N', U, 1d0, 0d0, M)
   if (all(M%elem - [3,2,1, -1,-2,0, 6,2,-1, 0,0,0,0,0] < tolerance) &
      .and. .not.equal_structure(B,M)) then
      write(*,*) "Test 2 passed"
   else
      write(*,*) "Test 2 FAILED"
   end if
   
   ! Test 3: A^T B C^T, compare result to a hand-computed one
   call triprod('T', A, 'N', B, 'T', C, 1d0, 0d0, M)
   if (all(M%elem - [27,30,51, 17,-3, 11,3, 10,2, 0,0,0,0,0] < tolerance) &
      .and. equal_structure(B,M)) then
      write(*,*) "Test 3 passed"
   else
      write(*,*) "Test 3 FAILED"
      write(*,*) "Equal structure ok?", equal_structure(B,M)
      write(*,*) "Result:", M%elem(1:9)
   end if
   
end program
