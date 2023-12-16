!------------------------------------------------------------------------------
! complex_quadrature.f90
!
! Some simple quadrature routines for contour integrals of complex functions.
!
! M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
!------------------------------------------------------------------------------
module complex_quadrature
   implicit none
   integer, parameter, private :: dp = kind(1d0)
   
contains
   
   !----------------------------------------------------------------------------
   ! Simple trapezoidal integration
   ! The array f contains the values of the integrand in points z along the
   ! contour; the array dz contains the derivative dz/dt in those points, and
   ! t is the contour parameter in those points
   !----------------------------------------------------------------------------
   complex(dp) function cquad_trapezoidal(f, dz, t)
      implicit none
      complex(dp), intent(in) :: f(:), dz(:)
      real(dp), intent(in) :: t(:)
      real(dp) :: re, im
      complex(dp) :: aux(size(f))
      
      aux(:) = dz(:)*f(:)
      re = quad_trapezoidal(real(aux,dp), t)
      im = quad_trapezoidal(aimag(aux), t)
      cquad_trapezoidal = re + im*(0,1)
      
   end function cquad_trapezoidal
   
   ! The real version
   real(dp) function quad_trapezoidal(f, x)
      implicit none
      real(dp), intent(in) :: f(:), x(:)
      
      integer :: i
      
      if (size(f) /= size(x)) stop "ERROR: f and x have mismatching sizes in quad_trapezoidal!"
      
      quad_trapezoidal = 0
      do i=1,size(f)-1
         quad_trapezoidal = quad_trapezoidal + (f(i+1) + f(i))/2 * (x(i+1) - x(i))
      end do
      
   end function quad_trapezoidal
   
   !----------------------------------------------------------------------------
   ! cquad_simpson
   ! Integration of complex function f(t) wrt real parameter t;
   ! e.g. f(t) = dz/dt*f(z(t))
   !
   ! Applies Simpson's 3/8 rule as far as possible, then falls back to the usual
   ! Simpson's rule or trapezoidal rule as necessary. Assumes equidistant
   ! spacing of the parameterization points t.
   !----------------------------------------------------------------------------
   function cquad_simpson(f, dz, t)
      implicit none
      real(dp),    dimension(:), intent(in) :: t
      complex(dp), dimension(:), intent(in) :: f, dz
      complex(dp) :: cquad_simpson
      
      complex(dp), dimension(size(f)) :: aux
      integer  :: npts, nint, ii, i1, i2, i3, i4
      real(dp) :: h

      cquad_simpson = 0
      
      if ((size(t) /= size(f)) .or. (size(t) /= size(dz))) then
         write(*,'(/A)') 'ERROR in cquad_simpson: dimensions of f, t, and dz are unequal.'
         stop
      end if
      
      npts = size(t)
      nint = npts-1
      h = t(2)-t(1)
      aux = dz(:)*f(:)
      
      ! Warn if there are a bad number of intervals
      if (mod(nint,3) /= 0) then
         write(*,'(A)') 'WARNING in fun.cquad_simpson: nint%3 /= 0.&
            & Errors will be larger than expected.'
      end if
      
      ! Integrate over all available sets of 3 intervals
      do ii=1, nint/3
         i1 = (ii-1)*3 + 1
         i2 = (ii-1)*3 + 2
         i3 = (ii-1)*3 + 3
         i4 = (ii-1)*3 + 4
         cquad_simpson = cquad_simpson + aux(i1) + 3*aux(i2) + 3*aux(i3) + aux(i4)
      end do
      cquad_simpson = 3.0_dp/8.0_dp*h*cquad_simpson
      
      ! What happens next depends on the number of remaining points. This
      ! quantity is npts-i4, and since i4 is the left-hand anchor of the next
      ! interval, npts-i4 is *also* the number of remaining intervals.
      select case (npts-i4)
         case (0)
            continue
         ! Trapezoid rule
         case (1)
            i1 = i4
            i2 = i4 + 1
            cquad_simpson = cquad_simpson + 0.5_dp*h*(aux(i1) + aux(i2))
         ! Simpson's rule
         case (2)
            i1 = i4
            i2 = i4 + 1
            i3 = i4 + 2
            cquad_simpson = cquad_simpson + h/3.0_dp*(aux(i1) + 4*aux(i2) + aux(i3))
         ! We messed up the counting
         case default
            write(*,'(/A)') 'WARNING in cquad_simpson: nint_remain > 2.'
      end select
   end function cquad_simpson
   
   !----------------------------------------------------------------------------
   ! Integrate a real function defined on the arrays X and Y.
   ! Use Simpson's 3/8 rule for all possible sets of 3 intervals, then handle
   ! the boundaries with either Simpson's rule (2 intervals) or the trapezoidal
   ! rule (1 interval) if necessary.
   !----------------------------------------------------------------------------
   function quad_simpson(x, y)
      implicit none
      real(dp), dimension(:), intent(in) :: x, y
      real(dp) :: quad_simpson
      
      real(dp), parameter :: one_three   = 1.0_dp/3.0_dp
      real(dp), parameter :: three_eight = 3.0_dp/8.0_dp
      
      integer :: npts, nintv, ii, i1, i2, i3, i4
      real(dp) :: dx, y1, y2, y3, y4
      
      if (size(x) /= size(y)) then
         write(*,'(/,A,1I0,A,1I0,A)') 'ERROR in fun.quad_simpson: dim(X) != dim(Y).&
            & Dim(X) = ', size(x), ', dim(Y) = ', size(y), '.'
         stop
      end if
      
      quad_simpson = 0
      npts  = size(x)
      nintv = npts - 1
      
      ! Warn if there are a bad number of intervals
      if (mod(nintv,3) /= 0) then
         write(*,'(A)') 'WARNING in fun.quad_simpson: nintv%3 /= 0.&
            & Errors will be larger than expected.'
      end if
      
      ! Check for constant dx
      do ii=3, npts
         if (abs((x(ii)-x(ii-1))**2 - (x(ii-1)-x(ii-2))**2) > 1d-10) then
            write(*,'(A)') 'WARNING in fun.quad_simpson: dx values vary by more&
               & than 1e-10.'
         end if
      end do
      
      ! Average dx
      dx = (x(npts)-x(1))/(1.0_dp*nintv)
      
      ! Simpson's 3/8 rule for the body of the work
      do ii=1, nintv/3
         i1 = (ii-1)*3 + 1;  y1 = y(i1)
         i2 = (ii-1)*3 + 2;  y2 = y(i2)
         i3 = (ii-1)*3 + 3;  y3 = y(i3)
         i4 = (ii-1)*3 + 4;  y4 = y(i4)
         
         quad_simpson = quad_simpson + three_eight*dx*(y1 + y4 + 3*(y2 + y3))
      end do
      
      ! Leftovers are handled with either Simpson's rule (2 intervals/3 points)
      ! or the trapezoidal rule (1 interval/2 points)
      select case (nintv-3*(ii-1))
         ! Simpson's rule
         case (2)
            i1 = i4;    y1 = y(i1)
            i2 = i4+1;  y2 = y(i2)
            i3 = i4+2;  y3 = y(i3)
            
            dx = (x(i3)-x(i1))*0.5_dp
            quad_simpson = quad_simpson + dx*one_three*(y1 + 4*y2 + y3)
         
         ! Trapezoidal rule
         case (1)
            i1 = i4;    y1 = y(i1)
            i2 = i4+1;  y2 = y(i2)
            
            quad_simpson = quad_simpson + 0.5_dp*(x(i2)-x(i1))*(y1+y2)
         
         ! Nothing to do
         case (0)
            continue
         
         ! Out of bounds
         case default
            write(*,'(A,1X,1I0,A)') 'ERROR in sub.integrate_array: N_leftover not in [0,2].', &
               'N_leftover = ', nintv-3*(ii-1), '.'
            stop
      end select
   end function quad_simpson
   
   !---------------------------------------------------------------------------
   ! Gauss-Legendre weights w and abscissae x (n-point rule) for the interval
   ! [a, b]
   !---------------------------------------------------------------------------
   subroutine gl_weights_abscissae(n, a, b, w, x)
      implicit none
      integer, intent(in) :: n
      real(dp), intent(out) :: w(n), x(n)
      real(dp), intent(in) :: a, b
      
      real(dp), parameter :: pi = 3.14159265359d0
      real(dp), parameter :: tol = 1d-10
      
      integer :: i
      real(dp) :: xi, diff, f(0:n), df(0:n), fm1(0:n), fm2(0:n)
      
      ! Construct the nth Legendre polynomial and its derivative using recursion
      ! rules
      f(:) = 0
      fm1(:) = 0
      fm1(0) = 1  ! P_0(x) = 1
      f(1) = 1  ! P_1(x) = x
      do i=2,n
         fm2(:) = fm1(:)
         fm1(:) = f(:)
         f(0) = 0
         f(1:n) = (2*i-1)*fm1(0:n-1)
         f(:) = f(:) - (i-1)*fm2(:)
         f(:) = f(:)/i
      end do
      df(:) = 0
      do i=1,n
         df(i-1) = i*f(i)
      end do
      
      ! Find the abscissae and weights
      do i=1,n
         ! Initial guess
         xi = cos(pi*(i-0.25d0)/(n+0.5d0))
         ! Newton-Rhapsod for finding the root
         do
            diff = polyeval(f, xi)/polyeval(df, xi)
            xi = xi - diff
            if (abs(diff) < tol) then
               ! Compute the weigth
               w(i) = 2/((1-xi*xi)*polyeval(df, xi)**2)
               x(i) = xi
               exit
            end if
         end do
      end do      
      
      ! Transform from the interval [-1,1] to [a,b]
      w(:) = w(:)*(b-a)/2
      x(:) = x(:)*(b-a)/2 + (a+b)/2
      
   end subroutine gl_weights_abscissae
   
   real(dp) function polyeval(f, x)
      implicit none
      real(dp), intent(in) :: f(:), x
      
      integer :: i
      
      polyeval = f(size(f))
      do i=size(f)-1,1,-1
         polyeval = x*polyeval + f(i)
      end do
      
   end function polyeval
   
end module complex_quadrature
