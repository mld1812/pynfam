!------------------------------------------------------------------------------
!> This module contains a derived type to store matrix elements of the
!> two-body current mean fields, separated by direct, exchange, and LEC
!>
!> @authors E.M. Ney, UNC Chapel Hill, 2018-
!------------------------------------------------------------------------------
module type_gamdel_2bc
   implicit none
   integer, parameter, private :: dp = kind(1.0d0)

   type gamdel_2bc
      real(dp), dimension(:), allocatable :: c3d, c3e, c4d, c4e, cpd, cpe
   end type gamdel_2bc

contains

   !----------------------------------------------------------------------------
   ! Allocate pieces
   !----------------------------------------------------------------------------
   subroutine allocate_gamdel(field, n)
      implicit none
      type(gamdel_2bc), intent(inout) :: field
      integer, intent(in) :: n

      call deallocate_gamdel(field)
      allocate(field%c3d(n)) ; allocate(field%c3e(n))
      allocate(field%c4d(n)) ; allocate(field%c4e(n))
      allocate(field%cpd(n)) ; allocate(field%cpe(n))

   end subroutine

   !----------------------------------------------------------------------------
   ! Deallocate pieces
   !----------------------------------------------------------------------------
   subroutine deallocate_gamdel(field)
      implicit none
      type(gamdel_2bc), intent(inout) :: field

      if (allocated(field%c3d)) deallocate(field%c3d)
      if (allocated(field%c4d)) deallocate(field%c4d)
      if (allocated(field%cpd)) deallocate(field%cpd)
      if (allocated(field%c3e)) deallocate(field%c3e)
      if (allocated(field%c4e)) deallocate(field%c4e)
      if (allocated(field%cpe)) deallocate(field%cpe)

   end subroutine

   !----------------------------------------------------------------------------
   ! Zero matrix elements
   !----------------------------------------------------------------------------
   subroutine zero_gamdel(field)
      implicit none
      type(gamdel_2bc), intent(inout) :: field

      field%c3d = 0 
      field%c3e = 0 
      field%c4d = 0 
      field%c4e = 0 
      field%cpd = 0 
      field%cpe = 0 

   end subroutine

   !----------------------------------------------------------------------------
   ! Fill matrix elements
   !----------------------------------------------------------------------------
   subroutine fill_gamdel(field, i, mtxels)
      implicit none
      type(gamdel_2bc), intent(inout) :: field
      integer, intent(in) :: i
      real(dp), intent(in) :: mtxels(6)

      field%c3d(i) = mtxels(1)
      field%c3e(i) = mtxels(2)
      field%c4d(i) = mtxels(3)
      field%c4e(i) = mtxels(4)
      field%cpd(i) = mtxels(5)
      field%cpe(i) = mtxels(6)

   end subroutine

   !----------------------------------------------------------------------------
   ! Combine into full mean field
   !----------------------------------------------------------------------------
   subroutine calc_totgamdel(field, tot, tde)
      implicit none
      type(gamdel_2bc), intent(in) :: field
      real(dp), dimension(:), intent(inout) :: tot
      character(1), intent(in), optional :: tde

      integer :: itde

      ! Handle optional input for (total, dir, exc)
      itde = 0
      if (present(tde)) then
         if (tde=='d') then
            itde = 1
         else if (tde=='e') then
            itde = -1
         else if (tde/='t') then
            stop 
         end if
      end if

      if (size(field%c3d) /= size(tot)) stop

      ! Combine direct, exchange, or both
      tot = 0
      if (itde >= 0) then
         tot = tot + field%c3d + field%c4d + field%cpd 
      end if
      if (itde <= 0) then
         tot = tot + field%c3e + field%c4e + field%cpe
      end if

   end subroutine

   !----------------------------------------------------------------------------
   ! Strip LECs
   !----------------------------------------------------------------------------
   subroutine strip_gamdel_lecs(field, c3, c4)
      implicit none
      type(gamdel_2bc), intent(inout) :: field
      real(dp), intent(in) :: c3, c4

      field%c3d = field%c3d/c3
      field%c3e = field%c3e/c3
      field%c4d = field%c4d/c4
      field%c4e = field%c4e/c4

   end subroutine

   !----------------------------------------------------------------------------
   ! Add in LECs
   !----------------------------------------------------------------------------
   subroutine mult_gamdel_lecs(field, c3, c4)
      implicit none
      type(gamdel_2bc), intent(inout) :: field
      real(dp), intent(in) :: c3, c4

      field%c3d = field%c3d*c3
      field%c3e = field%c3e*c3
      field%c4d = field%c4d*c4
      field%c4e = field%c4e*c4

   end subroutine

   !----------------------------------------------------------------------------
   ! Zero matrix elements
   !----------------------------------------------------------------------------
   subroutine reorder_gamdel_basis(field, op_mat, new_order, rc)
      use type_blockmatrix
      implicit none
      type(gamdel_2bc), intent(inout) :: field
      type(blockmatrix), intent(in) :: op_mat
      integer, intent(in) :: new_order(dqp)
      character(len=1), intent(in) :: rc
      type(blockmatrix) :: mat

      call copy_block_structure(op_mat, mat)

      mat%elem = field%c3d
      call reorder_blockmatrix_basis(mat, new_order, rc)
      field%c3d = mat%elem

      mat%elem = field%c3e
      call reorder_blockmatrix_basis(mat, new_order, rc)
      field%c3e = mat%elem

      mat%elem = field%c4d
      call reorder_blockmatrix_basis(mat, new_order, rc)
      field%c4d = mat%elem

      mat%elem = field%c4e
      call reorder_blockmatrix_basis(mat, new_order, rc)
      field%c4e = mat%elem

      mat%elem = field%cpd
      call reorder_blockmatrix_basis(mat, new_order, rc)
      field%cpd = mat%elem

      mat%elem = field%cpe
      call reorder_blockmatrix_basis(mat, new_order, rc)
      field%cpe = mat%elem

   end subroutine

end module type_gamdel_2bc
