!------------------------------------------------------------------------------
!> This module contains a derived type to store properties of a charge-changing
!> pnFAM external field.
!>
!> @authors M.T. Mustonen and T. Shafer, UNC Chapel Hill 2013-15
!> @authors E.M. Ney, UNC Chapel Hill, 2018-
!------------------------------------------------------------------------------
module type_extfield
   use type_blockmatrix, only : blockmatrix
   use type_gamdel_2bc, only : gamdel_2bc
   implicit none
   integer, parameter, private :: dp = kind(1.0d0)

   type external_field
      character(len=80) :: label
      integer :: k
      integer :: rank
      logical :: beta_minus
      logical :: parity_even
      integer :: use_2bc(6)
      type(blockmatrix) :: mat
      type(blockmatrix) :: mat12
      type(gamdel_2bc) :: gam
      type(gamdel_2bc) :: del
   end type external_field

end module type_extfield
