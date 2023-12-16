!------------------------------------------------------------------------------
! complex_quadrature_test.f90
!
! Simple analytically verifiable test cases for the complex contour integration
! routines.
!
! M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
!------------------------------------------------------------------------------
program test_complex_quadrature
   
   use complex_quadrature
   
   implicit none
   
   integer, parameter :: dp = kind(1d0)
   integer, parameter :: npts = 141
   real(dp), parameter :: pi = 3.14159265359d0
   complex(dp) :: z(npts), f(npts), dz(npts)
   real(dp) :: t(npts)
   real(dp) :: dt
   integer :: i
   
   integer, parameter :: cgpts = 5
   real(dp) :: w(cgpts), x(cgpts)
   complex(dp) :: f2(cgpts)
   
   write(*,*) ""
   
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ! Contour: line on the real axis
   dt = 1d0/(npts-1)
   t(:) = (/ (i*dt, i=0, npts-1) /)
   z(:) = t(:)
   dz(:) = 1  ! derivative of z = t w.r.t. t
   
   ! Integral over x is x^2/2
   f(:) = z(:)
   call test_complex_equal(cquad_trapezoidal(f, dz, t), (.5d0,0), 1d-5, "x^2")
   
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ! Contour: the unit circle
   dt = 2*pi/(npts-1)
   t(:) = (/ (i*dt, i=0, npts-1) /)
   z(:) = cos(t(:)) + (0,1)*sin(t(:))
   dz(:) = -sin(t(:)) + (0,1)*cos(t(:))
   
   ! Integral over exp(2/z) is 4 pi i
   f(:) = exp(2d0/z(:))
   call test_complex_equal(cquad_trapezoidal(f, dz, t), 4*pi*(0,1), 1d-5, "exp(2/z)")
   
   ! Integral over exp(z)/z^5 is pi i/12
   f(:) = exp(z(:))/(z(:)**5)
   call test_complex_equal(cquad_trapezoidal(f, dz, t), pi*(0,1)/12d0, 1d-5, "exp(z)/z^5")
   
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ! Test GL quadrature
   call gl_weights_abscissae(cgpts, -1d0, 1d0, w, x)
   f2(:) = x(:)
   call test_complex_equal(sum(w(:)*f2(:)), (0d0,0d0), 1d-6, "GL x")
   
   call gl_weights_abscissae(cgpts, 0d0, pi, w, x)
   f2(:) = sin(x(:))
   call test_complex_equal(sum(w(:)*f2(:)), (2d0,0d0), 1d-6, "GL sin x")
   
contains
   
   subroutine test_complex_equal(a, b, tol, testname)
      implicit none
      complex(dp), intent(in) :: a, b
      real(dp), intent(in) :: tol
      character(len=*), intent(in) :: testname
      
      ! A little bling to make it easier to spot any errors (using ANSI colors)
      character, parameter :: escape_char = achar(27)
      character(len=*), parameter :: cc_reset = escape_char // "[0m"
      character(len=*), parameter :: fg_green = escape_char // "[0;32m"
      character(len=*), parameter :: fg_red = escape_char // "[0;31m"
      
      write(*,*) "Test: ", testname
      write(*,*) "Result:   ", a
      write(*,*) "Expected: ", b
      if (abs(a-b) < tol) then
         write(*,*) fg_green, "Test passed."
      else
         write(*,*) fg_red, "*** TEST FAILED! ***"
      end if
      write(*,*) cc_reset
      
   end subroutine
   
end program
