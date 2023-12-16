!------------------------------------------------------------------------------
! betadecay.f90
!
! This file only wraps the primary routine of the betadecay module into a
! Fortran program (so that the module itself can be called from other codes
! as well).
! 
! M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
!------------------------------------------------------------------------------
program betadecay_main
   use betadecay, only : betadecay_main_program
   implicit none
   
   call betadecay_main_program
   
end program
