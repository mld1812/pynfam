!------------------------------------------------------------------------------
!> This module contains a derived type to hold the FAM block matrices, and
!> routines for manipulating them. Importantly, this type assumes for every
!> block row, there is only one non-zero block column.
!>
!> @authors M.T. Mustonen and T. Shafer, UNC Chapel Hill 2013-15
!> @authors E.M. Ney, UNC Chapel Hill, 2018-
!------------------------------------------------------------------------------
module type_blockmatrix
   implicit none

   integer, parameter, private :: dp = kind(1.0d0)
   logical, private :: init=.false.

   integer              :: nb         !> Number of blocks (only one block per [row,column])
   integer              :: dqp        !> Number of rows/cols in entire matrix
   integer, allocatable :: db(:)      !> Number of basis states per block
   integer, allocatable :: isstart(:) !> Index of first state (i.e. first col or row) in block i

   ! The storage order is column-major within blocks.
   ! blocks stored from top to bottom (increasing row index).
   type blockmatrix
      integer, allocatable :: ic2r(:), ir2c(:)
      integer, allocatable :: ic2m(:), ir2m(:)
      real(dp), allocatable :: elem(:)
   end type blockmatrix

contains

   !---------------------------------------------------------------------------
   ! I rescind my previous thoughts on nb/db here. I'm going to treat them as
   ! "class variables", since they are the same for all blockmatrices.
   ! However, a method to initialize the class still seems like a good idea.
   ! I also include isstart here, as they are just functions of nb/db.
   !---------------------------------------------------------------------------
   subroutine init_blockmatrix_type(nb_in, db_in, ierr)
      implicit none
      integer, intent(in)  :: nb_in
      integer, intent(in)  :: db_in(:)
      integer, intent(out) :: ierr
      integer :: i

      ierr = 0
      if(.not. init) then
         nb = nb_in
         db = db_in
         If(Allocated(isstart)) Deallocate(isstart)
         allocate(isstart(nb_in))
         isstart(1) = 1
         do i=1,nb_in-1
            isstart(i+1) = sum(db_in(1:i))+1
         end do
         dqp = sum(db)
         init = .true.
      else
         if (nb/=nb_in .or. .not. all(db==db_in)) then
            ierr = 1
         end if
      end if

   end subroutine init_blockmatrix_type

   !---------------------------------------------------------------------------
   ! Allocates a block matrix a to be able to hold n elements and zeroes its
   ! contents.
   !---------------------------------------------------------------------------
   subroutine allocate_blockmatrix(a, n)
      implicit none
      type(blockmatrix), intent(inout) :: a
      integer, intent(in) :: n

      if(.not.init) return
      
      call deallocate_blockmatrix(a)
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
   ! Sets the block matrix to the identity
   !---------------------------------------------------------------------------
   subroutine set_identity(a)
      implicit none
      type(blockmatrix), intent(inout) :: a
      integer :: i, im, ib, nd, ir, bij

      ! Block diagonal structure
      a%ir2c = [ (i, i=1, nb) ]; a%ic2r = a%ir2c
      a%ir2m(1) = 1
      do i=2,nb
        a%ir2m(i) = a%ir2m(i-1) + db(i-1)**2
      end do
      a%ic2m = a%ir2m

      ! Fill diagonal elements of each block with 1s
      a%elem = 0
      im = 0
      do ib=1,nb
         nd = db(ib)
         do ir=1,nd
            bij = ir + (ir-1)*nd
            a%elem(im + bij) = 1.0_dp
         end do
         im = im + nd*nd
      end do

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

      if(.not.init) return
      
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

   !---------------------------------------------------------------------------
   ! Subroutine to reorder basis states of a block matrix type, thus allowing
   ! reordering of matrices which are not block diagonal.
   !---------------------------------------------------------------------------
   subroutine reorder_blockmatrix_basis(U, new_order, rc)
      type(blockmatrix), intent(inout) :: U
      integer, intent(in) :: new_order(dqp)
      character(len=1), intent(in) :: rc
      integer :: ibr, ibc, ir, ic, im, dr, nor, noc
      real(dp), allocatable, dimension(:) :: U_new

      if (allocated(U_new)) deallocate(U_new)
      allocate(U_new(size(U%elem)))

      im = 0
      do ibr=1, nb
         ibc = U%ir2c(ibr)
         if (ibc == 0) cycle
         do ic=1, db(ibc)
            do ir=1, db(ibr)
               dr  = db(ibr)
               nor = new_order(isstart(ibr)-1+ir) - (isstart(ibr)-1)
               noc = new_order(isstart(ibc)-1+ic) - (isstart(ibc)-1)
               im  = U%ir2m(ibr)-1
               select case(rc)
                  case('r')
                     U_new(im + ir + (ic-1)*dr) = U%elem(im + nor + (ic-1)*dr)
                  case('c')
                     U_new(im + ir + (ic-1)*dr) = U%elem(im + ir + (noc-1)*dr)
                  case('a')
                     U_new(im + ir + (ic-1)*dr) = U%elem(im + nor + (noc-1)*dr)
                  case default
                     write(*,*) "Error: value error for parameter rc"
                     stop
               end select
            end do
         end do
      end do

      U%elem = U_new
      deallocate(U_new)

   end subroutine reorder_blockmatrix_basis

end module type_blockmatrix
