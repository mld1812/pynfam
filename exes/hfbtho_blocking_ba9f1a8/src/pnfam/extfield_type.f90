!------------------------------------------------------------------------------
! extfield_type.f90
!
! A custom type to store a pnFAM external field.
!
! M.T. Mustonen and T. Shafer, UNC Chapel Hill 2013-15
!------------------------------------------------------------------------------
module external_field_type
   use blockmatrix_type
   implicit none
   integer, parameter, private :: dp = kind(1.0d0)
   
   type external_field
      character(len=80) :: label
      integer :: k
      integer :: rank
      logical :: beta_minus
      logical :: parity_even
      type(blockmatrix) :: mat
   end type external_field

end module external_field_type
