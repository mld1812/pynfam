!------------------------------------------------------------------------------
! blockmatrix_type.f90
!
! A custom type to hold the FAM block matrices.
!
! M.T. Mustonen and T. Shafer, UNC Chapel Hill 2013-15
!------------------------------------------------------------------------------
module blockmatrix_type
   implicit none
   integer, parameter, private :: dp = kind(1d0)
   
   ! The grid for the block structure
   ! This is the same for every block matrix we have
   integer :: nb
   integer, allocatable :: db(:)
   
   ! The storage order is column-major within blocks, blocks stored from top to
   ! bottom (increasing row index).
   type blockmatrix
      integer, allocatable :: ic2r(:), ir2c(:)
      integer, allocatable :: ic2m(:), ir2m(:)
      real(dp), allocatable :: elem(:)
   end type blockmatrix
   
contains
   
   !---------------------------------------------------------------------------
   ! Allocates a block matrix a to be able to hold n elements and zeroes its
   ! contents.
   !---------------------------------------------------------------------------
   subroutine allocate_blockmatrix(a, n)
      implicit none
      type(blockmatrix), intent(inout) :: a
      integer, intent(in) :: n
      
      allocate(a%ic2r(nb), a%ir2c(nb), a%ir2m(nb), a%ic2m(nb))
      allocate(a%elem(n))
      
      a%ic2r = 0 ; a%ir2c = 0 ; a%ir2m = 0 ; a%ic2m = 0
      a%elem = 0
      
   end subroutine
   
   
   subroutine deallocate_blockmatrix(a)
      implicit none
      type(blockmatrix), intent(inout) :: a
      
      if (allocated(a%ic2r)) deallocate(a%ic2r, a%ir2c, a%ir2m, a%ic2m, a%elem)
      
   end subroutine
   
   
   !---------------------------------------------------------------------------
   ! Copies the structure of A to B
   !---------------------------------------------------------------------------
   subroutine copy_block_structure(a, b)
      implicit none
      type(blockmatrix), intent(in) :: a
      type(blockmatrix) :: b
      
      b%ic2r = a%ic2r ; b%ir2c = a%ir2c ; b%ic2m = a%ic2m ; b%ir2m = a%ir2m
      
   end subroutine
   
   
   !---------------------------------------------------------------------------
   ! BLAS-accelerated product of three block matrices
   ! IMPORTANT: ABC is expected to have been allocated beforehand and A and C
   ! must have square blocks (B can have rectangular blocks)
   !---------------------------------------------------------------------------
   subroutine triprod(transa, a, transb, b, transc, c, alpha, beta, abc)
      implicit none
      type(blockmatrix), intent(in) :: a, b, c
      character, intent(in) :: transa, transb, transc
      real(dp), intent(in) :: alpha, beta
      type(blockmatrix), intent(inout) :: abc
      
      integer :: i, j, k, l, ipt, ipa, ipb, ipc, nda, ndc, ldb
      real(dp) :: aux(size(abc%elem))
      
      ! Clear the index arrays of the result matrix
      abc%ir2c(:) = 0 ; abc%ic2r(:) = 0
      abc%ir2m(:) = 0 ; abc%ic2m(:) = 0
      
      ! Iterate in the storage order of the result matrix ABC
      ipt = 1
      do i = 1, nb  ! the block row index of ABC
         
         ! Determine the block indices in ABC_[i][j] = A_[i][k] B_[k][l] C_[l][j]
         ! and the starting index of the block in the storage of each matrix
         if (transa == 'n' .or. transa == 'N') then
            k = a%ir2c(i)
            ipa = a%ir2m(i)
         else
            k = a%ic2r(i)
            ipa = a%ic2m(i)
         end if
         if (k == 0) cycle
         
         if (transb == 'n' .or. transb == 'N') then
            l = b%ir2c(k)
            ipb = b%ir2m(k)
         else
            l = b%ic2r(k)
            ipb = b%ic2m(k)
         end if
         if (l == 0) cycle
         
         if (transc == 'n' .or. transc == 'N') then
            j = c%ir2c(l)
            ipc = c%ir2m(l)
         else
            j = c%ic2r(l)
            ipc = c%ic2m(l)
         end if
         if (j == 0) cycle
         
         ! Store the indices for ABC
         abc%ir2c(i) = j ; abc%ic2r(j) = i
         abc%ir2m(i) = ipt ; abc%ic2m(j) = ipt
         
         ! Determine the dimensions of the blocks
         nda = db(i) ; ndc = db(j)
         if (transb == 'n' .or. transb == 'N') then
            ldb = nda
         else
            ldb = ndc
         end if
         
         ! The block of A is nda x nda, the block of B nda x ndc, so AB is nda x ndc
         call dgemm(transa, transb, nda,ndc,nda, alpha,a%elem(ipa),nda, b%elem(ipb), &
            ldb, 0d0, aux(ipt),nda)
         ! The block of AB is nda x ndc, the block of C ndc x ndc, so ABC is nda x ndc
         call dgemm('n', transc, nda,ndc,ndc, 1d0,aux(ipt),nda, c%elem(ipc),ndc, beta, &
            abc%elem(ipt),nda)
         
         ipt = ipt + nda*ndc
         
      end do
      
   end subroutine triprod
   
   !---------------------------------------------------------------------------
   ! Debugging routine to check if two block matrices A and B have the same
   ! block structure.
   !---------------------------------------------------------------------------
   logical function equal_structure(a, b)
      implicit none
      type(blockmatrix), intent(in) :: a, b
      
      ! Actually the following all-statements should be redundant, but if any
      ! of them fails it's a sign of inconsistency in the code.
      equal_structure = all(a%ir2c == b%ir2c) .and. all(a%ic2r == b%ic2r) .and. &
         all(a%ir2m == b%ir2m) .and. all(a%ic2m == b%ic2m)
      
   end function
   
   
   !---------------------------------------------------------------------------
   ! Elemental multiplication of two complex block matrices split to real and
   ! imaginary arrays: (reab,imab) = (rea,ima)*(reb,imb)
   !---------------------------------------------------------------------------
   subroutine complex_multiply(rea, ima, reb, imb, reab, imab)
      implicit none
      type(blockmatrix), intent(in) :: rea, ima, reb, imb
      type(blockmatrix), intent(inout) :: reab, imab
      
      reab%elem(:) = rea%elem*reb%elem - ima%elem*imb%elem
      imab%elem(:) = ima%elem*reb%elem + imb%elem*rea%elem
      
   end subroutine complex_multiply
   
end module blockmatrix_type
