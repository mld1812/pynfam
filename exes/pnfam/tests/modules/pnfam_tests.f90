!------------------------------------------------------------------------------
! pnfam_tests.f90
!
! M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
!------------------------------------------------------------------------------
module pnfam_tests

   implicit none
   
   integer, parameter :: tsp = kind(1.0)
   integer, parameter :: tdp = kind(1.0d0)
   
   ! A little bling to make it easier to tspot any errors (using ANSI colors)
   character(len=*), parameter :: escape_char = achar(27)
   character(len=*), parameter :: cc_reset    = escape_char // "[0m"
   character(len=*), parameter :: fg_green    = escape_char // "[0;32m"
   character(len=*), parameter :: fg_red      = escape_char // "[0;31m"
   
   ! Generic testing interface
   interface test_equal
      module procedure test_complex_equal_double, test_complex_equal_single, &
         test_real_equal_double, test_real_equal_single
   end interface
   
contains
   
   !----------------------------------------------------------------------------
   ! From Mika, complex_quadrature_test.f90 
   !----------------------------------------------------------------------------
   subroutine test_complex_equal_double(a, b, tol, testname, passed)
      implicit none
      
      complex(tdp),      intent(in)  :: a, b
      real(tdp),         intent(in)  :: tol
      character(len=*),  intent(in)  :: testname
      logical, optional, intent(out) :: passed
      
      logical :: did_pass = .false.
      
      write(*,*) "Test ........ : ", testname
      write(*,*) "Result ...... : ", a
      write(*,*) "Expected .... : ", b
      write(*,*) "Tolerance ... : ", tol
      
      if (abs(a-b) < tol) then
         write(*,*) fg_green, "Test passed."
         did_pass = .true.
      else
         write(*,*) fg_red, "*** TEST FAILED! ***"
         did_pass = .false.
      end if
      write(*,*) cc_reset
      
      if (present(passed)) passed = did_pass
      
   end subroutine

   !----------------------------------------------------------------------------
   ! From Mika, complex_quadrature_test.f90 
   !----------------------------------------------------------------------------
   subroutine test_complex_equal_single(a, b, tol, testname, passed)
      implicit none
      
      complex(tsp),      intent(in)  :: a, b
      real(tsp),         intent(in)  :: tol
      character(len=*),  intent(in)  :: testname
      logical, optional, intent(out) :: passed
      
      logical :: did_pass = .false.
      
      write(*,*) "Test ........ : ", testname
      write(*,*) "Result ...... : ", a
      write(*,*) "Expected .... : ", b
      write(*,*) "Tolerance ... : ", tol
      
      if (abs(a-b) < tol) then
         write(*,*) fg_green, "Test passed."
         did_pass = .true.
      else
         write(*,*) fg_red, "*** TEST FAILED! ***"
         did_pass = .false.
      end if
      write(*,*) cc_reset
      
      if (present(passed)) passed = did_pass
      
   end subroutine

   !----------------------------------------------------------------------------
   ! From Mika, complex_quadrature_test.f90 
   !----------------------------------------------------------------------------
   subroutine test_real_equal_double(a, b, tol, testname, passed)
      implicit none
      
      real(tdp),         intent(in)  :: a, b
      real(tdp),         intent(in)  :: tol
      character(len=*),  intent(in)  :: testname
      logical, optional, intent(out) :: passed
      
      logical :: did_pass = .false.
      
      write(*,*) "Test ........ : ", testname
      write(*,*) "Result ...... : ", a
      write(*,*) "Expected .... : ", b
      write(*,*) "Tolerance ... : ", tol
      
      if (abs(a-b) < tol) then
         write(*,*) fg_green, "Test passed."
         did_pass = .true.
      else
         write(*,*) fg_red, "*** TEST FAILED! ***"
         did_pass = .false.
      end if
      write(*,*) cc_reset
      
      if (present(passed)) passed = did_pass
      
   end subroutine
   
   !----------------------------------------------------------------------------
   ! From Mika, complex_quadrature_test.f90 
   !----------------------------------------------------------------------------
   subroutine test_real_equal_single(a, b, tol, testname, passed)
      implicit none
      
      real(tsp),         intent(in)  :: a, b
      real(tsp),         intent(in)  :: tol
      character(len=*),  intent(in)  :: testname
      logical, optional, intent(out) :: passed
      
      logical :: did_pass = .false.
      
      write(*,*) "Test ........ : ", testname
      write(*,*) "Result ...... : ", a
      write(*,*) "Expected .... : ", b
      write(*,*) "Tolerance ... : ", tol
      
      if (abs(a-b) < tol) then
         write(*,*) fg_green, "Test passed."
         did_pass = .true.
      else
         write(*,*) fg_red, "*** TEST FAILED! ***"
         did_pass = .false.
      end if
      write(*,*) cc_reset
      
      if (present(passed)) passed = did_pass
      
   end subroutine
   
end module pnfam_tests
