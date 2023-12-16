!------------------------------------------------------------------------------
! hfb_solution_type.f90
!
! A custom type to store an HFB solution.  Unfortunately the data is scattered
! around.  This type only holds the data which is not publicly available from
! hfbtho_basis.f90.
!
! M.T. Mustonen and T. Shafer, UNC Chapel Hill 2015
!------------------------------------------------------------------------------
module hfb_solution_type
   use blockmatrix_type
   implicit none
   integer, parameter, private :: dp = kind(1.0d0)

   type hfb_solution
      real(dp), allocatable :: en(:), ep(:)
      type(blockmatrix) :: up, un, vp, vn
   end type

end module hfb_solution_type
