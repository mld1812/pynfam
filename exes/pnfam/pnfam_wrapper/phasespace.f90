!------------------------------------------------------------------------------
! phasespace.f90
!
! This module contains the analytic continuations of the phase-space integrals
! for computing the beta decay rates (or half-lives).
!
! M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
!------------------------------------------------------------------------------

module phasespace
   use rational_interp
   implicit none
   integer, parameter, private :: dp = kind(1d0)
   
contains

   !============================================================================
   !                        POLYNOMIAL FIT ROUTINES
   !============================================================================
   !----------------------------------------------------------------------------
   ! A wrapper for the C routine which performs a chi square polynomial fit
   ! and returns the polynomial functions
   !----------------------------------------------------------------------------
   subroutine poly_fit(nr_fit, nr_data, x, y, c)
      use iso_c_binding
      implicit none
      interface
         integer function fit_poly(nr_fit, nr_data_points, xs, ys, cs) &
            bind(C, name="fit_poly")
            use iso_c_binding
            implicit none
            integer(c_int), intent(in), value :: nr_fit
            integer(c_int), intent(in), value :: nr_data_points
            real(c_double), intent(in) :: xs(nr_data_points)
            real(c_double), intent(in) :: ys(nr_data_points)
            real(c_double), intent(out) :: cs(nr_fit)
         end function
      end interface
      integer, intent(in) :: nr_fit, nr_data
      real(dp), intent(in) :: x(nr_data), y(nr_data)
      real(dp), intent(out) :: c(nr_fit)
      
      integer :: ierr
      real(c_double) :: xs(nr_data), ys(nr_data), cs(nr_fit)
      
      xs = x ; ys = y
      ierr = fit_poly(nr_fit, nr_data, xs, ys, cs)
      c = cs
   end subroutine

   ! Evaluating the analytic continuation of the above fit
   complex(dp) function polyfit_eval(order, cs, z)
      implicit none
      integer, intent(in) :: order
      real(dp), intent(in) :: cs(order)
      complex(dp), intent(in) :: z
      
      integer :: i
      
      polyfit_eval = cs(order)*(z-1)
      do i=order-1,1,-1
         polyfit_eval = (z-1)*(polyfit_eval + cs(i))
      end do
   end function polyfit_eval
   
   !----------------------------------------------------------------------------
   ! Construct a polynomial fit to a phase-space integral 
   !----------------------------------------------------------------------------
   function get_polyfit_coeffs(nz, na, order, npts, w0max, nfi, si) result(cs)
      implicit none

      integer,  intent(in) :: nz, na, order, npts, nfi
      real(dp), intent(in) :: w0max, si
      real(dp), dimension(order) :: cs
      
      real(dp), parameter :: w0min = 1.0_dp
      integer  :: i
      real(dp) :: dw0, w0(npts), ps(npts)
      
      cs(:) = 0
      
      ! Compute w0-values and phase-space integrals
      dw0 = (w0max-w0min)/(npts-1)
      w0 = 0;  ps = 0
      
      do i=1, npts
         w0(i) = w0min + (i-1)*dw0
         ps(i) = real_psi(nz, na, w0(i), nfi, si)
      end do
      
      call poly_fit(order, npts, w0, ps, cs)
   end function get_polyfit_coeffs


   !============================================================================
   !                  RATIONAL FUNCTION INTERPOLATION ROUTINES
   !============================================================================
   !----------------------------------------------------------------------------
   !> Construct a rational function interpolation to phase-space integral points
   !> calculated along the real axis.
   !> The grid for w0 are chebychev roots, which improves the rational function
   !> interpolation as long as there are no poles near the center of our
   !> interpolation interval.
   !----------------------------------------------------------------------------
   function get_ratint_coeffs(nz, na, npts, w0max, gn, nleg) result(cs)
      use complex_quadrature, only: gaucheby
      implicit none

      integer,  intent(in) :: nz, na, npts, gn, nleg
      real(dp), intent(in) :: w0max
      real(dp), dimension(npts) :: cs

      real(dp), parameter :: w0min = 1.0_dp
      real(dp), dimension(npts) :: w0_grid, ps
      real(dp) :: ps_check, dy_dummy, weight_dummy(npts), dw0
      integer  :: i

      cs(:) = 0

      ! Construct the grid to compute interpolation points on
      call gaucheby(npts, w0min, w0max, weight_dummy, w0_grid)

      ! Compute phase space integral on the grid
      ps = 0
      do i=1,npts
         ps(i) = real_psi_leg(nz, na, w0_grid(i), gn, nleg)
      end do

      ! Get the rational function thiele coefficients
      call thiele_coeffs(npts, w0_grid, ps, cs)

      ! Check the function actually passes through the interpolation points
      do i=1,npts
         call ratint_eval_real(npts, w0_grid, cs, w0_grid(i), ps_check, dy_dummy)
         if ((ps(i) - ps_check)/(ps(i)+1e-25) > 1e-5) then
             write(*,'(A,1es10.3,A)') "WARNING: Rational function interpolation does not pass &
                & through interpolation points precisely (relative difference is ", &
                (ps(i) - ps_check)/(ps(i)+1e-25), ")."
         end if
      end do
   end function get_ratint_coeffs

   !----------------------------------------------------------------------------
   !> Evaluate a rational function interpolation of psi's along the real axis at
   !> with a COMPLEX argument. This function is a wrapper so we can pass w0max
   !> in place of the whole grid of 'x' values.
   !>
   !> NB!: make sure the 'x' here matches the 'x' in 'get_ratint_coeffs'!!!
   !----------------------------------------------------------------------------
   function ratint_eval_psi(npts, w0max, c, z) result(y)
      use complex_quadrature, only: gaucheby
      implicit none
      integer, intent(in) :: npts
      real(dp), intent(in) :: c(npts), w0max
      complex(dp), intent(in) :: z

      complex(dp) :: y, dy_dummy
      real(dp), parameter :: w0min = 1.0_dp
      real(dp) :: x(npts), weight_dummy(npts)
      integer :: i

      ! the grid (should be the same as used for get_ratint_coeffs)
      call gaucheby(npts, w0min, w0max, weight_dummy, x)
      ! Evaluate for a complex argument
      call ratint_eval_cmplx(npts, x, c, z, y, dy_dummy)
   end function ratint_eval_psi
   !----------------------------------------------------------------------------
   !> Evaluate a rational function interpolation of psi's along the real axis at
   !> with a REAL argument. This function is a wrapper so we can pass w0max
   !> in place of the whole grid of 'x' values.
   !>
   !> NB!: make sure the 'x' here matches the 'x' in 'get_ratint_coeffs'!!!
   !----------------------------------------------------------------------------
   function ratint_eval_psi_real(npts, w0max, c, z) result(y)
      use complex_quadrature, only: gaucheby
      implicit none
      integer, intent(in) :: npts
      real(dp), intent(in) :: c(npts), w0max
      real(dp), intent(in) :: z

      real(dp) :: y, dy_dummy
      real(dp), parameter :: w0min = 1.0_dp
      real(dp) :: x(npts), weight_dummy(npts)
      integer :: i

      ! the grid (should be the same as used for get_ratint_coeffs)
      call gaucheby(npts, w0min, w0max, weight_dummy, x)
      ! Evaluate for a complex argument
      call ratint_eval_real(npts, x, c, z, y, dy_dummy)
   end function ratint_eval_psi_real


   !============================================================================
   !                        PHASE SPACE INTEGRAL ROUTINES
   !============================================================================
   !----------------------------------------------------------------------------
   !> The Primakoff-Rosen approximation for the allowed phase-space integral.
   !> This function is no longer used.
   !----------------------------------------------------------------------------
   complex(dp) elemental function pr_psi(Z, E)
      implicit none
      complex(dp), intent(in) :: E
      integer, intent(in) :: Z
      real(dp), parameter :: alpha = 0.00729735256d0  ! the fine-structure constant
      real(dp), parameter :: pi = 3.14159265359d0
      
      pr_psi = (E**5 - 10*E*E + 15*E - 6)/30
      pr_psi = pr_psi*(2*pi*alpha*Z)/(1-exp(-2*pi*alpha*Z))
   end function pr_psi

   !----------------------------------------------------------------------------
   !> Compute a phase-space integral along the real axis from w=1,w0. For the
   !> contour integration, w0 becomes Re(EQRPA) and we compute the psi on a grid
   !> of w0's.
   !>
   !> NOTE: this does not check for NaN values the way a previous routine did. As
   !> a result, we may want to implement an interation cap or NaN test.
   !>
   !> nz    - Z (daughter)
   !> na    - A
   !> w0    - maximum electron energy scaled by mec2.
   !> n     - index of the shape; these correspond to J. Suhonen, Nucl. Phys. A
   !>         563, 205 (1993).
   !> si    - integration tolerance; 1d-5 should get 3-4 significant figures.
   !> debug - pass .true. to see crude statistics.
   !>
   !>----------------------------------------------------------------------------
   !>
   !> Compute the phase space integral from 1 to w0
   !>
   !> Method: Simpsons 3/8 with mesh refinement until convergence criterion
   !----------------------------------------------------------------------------
   function real_psi(nz, na, w0, n, si, debug)
      use complex_quadrature, only: quad_simpson
      use fermi, only: L0_us, F0_us, p2lambda2_us, g1 => gamma_1
   
      implicit none
   
      integer,  parameter :: dp = kind(1.0d0)
      real(dp), parameter :: wmin = 1.0_dp

      integer,  intent(in) :: nz, na, n
      real(dp), intent(in) :: w0, si
      logical,  intent(in), optional :: debug
   
      real(dp) :: real_psi
      logical  :: debug_on = .false.
   
      logical  :: encountered_nan
      integer  :: npts, is, i
      real(dp) :: dw, result, prev_result, si_achieved
      real(dp) :: wx, px
      real(dp), allocatable :: w(:), f(:)

   
      if (present(debug)) debug_on = debug
   
      real_psi = 0
      
      ! Initial value for having found a NaN point
      encountered_nan = .false.
   
      ! Initial value for npts can be some hand-wavy function of w0
      npts = 100
      if (mod(npts,3) == 0) npts = npts+1
      if (mod(npts,3) == 2) npts = npts+2
   
      ! Initial values
      ! is = loop index, prev_result = storage for last iteration's result,
      ! si_achieved = integration delta
      is = 0;  prev_result = 0;  si_achieved = 0
      do
         is = is+1
      
         ! Increment the loop by some number of points each iteration;
         ! should be a multiple of 3
         npts = npts + 30
            
         allocate(w(npts), f(npts))
         dw = (w0-wmin)/(npts-1)
      
         ! Compute the integral on the mesh and integrate
         do i=1, npts
            w(i) = wmin + (i-1)*dw
         
            wx = w(i)
            px = sqrt(wx**2-1)
         
            select case (n)
               case (1)
                  f(i) = px*(w0-wx)**2*L0_us(nz)*F0_us(nz,na,wx)*g1(nz)
               case (2)
                  f(i) = px*(w0-wx)**2*L0_us(nz)*F0_us(nz,na,wx)*wx
               case (3)
                  f(i) = px*(w0-wx)**2*L0_us(nz)*F0_us(nz,na,wx)*wx**2
               case (4)
                  f(i) = px*(w0-wx)**2*L0_us(nz)*F0_us(nz,na,wx)*wx**3
               case (5)
                  f(i) = px*(w0-wx)**2*L0_us(nz)*F0_us(nz,na,wx)*wx*(w0-wx)**2
               case (6)
                  f(i) = px*(w0-wx)**2*L0_us(nz)*F0_us(nz,na,wx)*wx*p2lambda2_us(nz,na,wx)
            end select
         end do
         
         ! Fun NaN check from StackOverflow
         ! http://stackoverflow.com/a/17391586/656740, acc. 12/12/2013
         ! If we see any NaN results, log them
         if (any(f /= f)) then
            ! If we haven't seen one before, print out that we found one
            if (.not.encountered_nan) then
               encountered_nan = .true.
               write(*,'(a)') '! Found NaN value in fun.real_psi'
            end if
            ! If there is a NaN, cycle the computation (don't record this point)
            deallocate(w,f)
            cycle
         end if
      
         result = quad_simpson(x=w, y=f)
         deallocate(w,f)
      
         ! SI computation (for iterations > 1)
         ! If we are converging to zero, SI is abs of the difference;
         ! If non-zero, SI is the relative difference
         if (is > 1) then
            if (prev_result == 0d0) then
               si_achieved = abs(result-prev_result)
            else
               si_achieved = abs((result-prev_result)/prev_result)
            end if
         end if
      
         ! If we achieves the tolerance requested, exit;
         ! otherwise sotre the result and re-cycle
         if (is > 1 .and. si_achieved <= si) then
            if (debug_on) then
               write(*,*) is, npts, si_achieved, prev_result, result
            end if
            exit
         else
            prev_result = result
            cycle
         end if
      end do

      real_psi = result
      return
   end function real_psi

   !----------------------------------------------------------------------------
   !> Compute the phase space integral from 1 to w0
   !>
   !> Method: gauss-legendre quadrature
   !----------------------------------------------------------------------------
   function real_psi_leg(nz, na, w0, n, npts, debug)
      use complex_quadrature, only: gauleg
      use fermi, only: L0_us, F0_us, p2lambda2_us, g1 => gamma_1
   
      implicit none

      integer,  parameter :: dp = kind(1.0d0)
      real(dp), parameter :: wmin = 1.0_dp

      integer,  intent(in) :: nz, na, n, npts
      real(dp), intent(in) :: w0
      logical,  intent(in), optional :: debug

      real(dp) :: real_psi_leg
      logical  :: debug_on = .false.

      integer  :: is, i
      real(dp) :: dw, result, prev_result
      real(dp) :: wx, px
      real(dp), allocatable :: w(:), f(:)

      real(dp), allocatable :: weights(:)

      allocate(w(npts), f(npts), weights(npts))
      w=0.0_dp; f=0.0_dp; weights=0.0_dp
      result = 0.0_dp

      call gauleg(npts, wmin, w0, weights, w)

      do i=1, npts
         wx = w(i)
         px = sqrt(wx**2.0_dp - 1.0_dp)
         select case (n)
            case (1)
               f(i) = px*(w0-wx)**2*L0_us(nz)*F0_us(nz,na,wx)*g1(nz)
            case (2)
               f(i) = px*(w0-wx)**2*L0_us(nz)*F0_us(nz,na,wx)*wx
            case (3)
               f(i) = px*(w0-wx)**2*L0_us(nz)*F0_us(nz,na,wx)*wx**2
            case (4)
               f(i) = px*(w0-wx)**2*L0_us(nz)*F0_us(nz,na,wx)*wx**3
            case (5)
               f(i) = px*(w0-wx)**2*L0_us(nz)*F0_us(nz,na,wx)*wx*(w0-wx)**2
            case (6)
               f(i) = px*(w0-wx)**2*L0_us(nz)*F0_us(nz,na,wx)*wx*p2lambda2_us(nz,na,wx)
         end select

         ! NaN is not equal to anything, even itself
         if (f(i) /= f(i)) then
            write(*,'(a)') 'WARNING: Found NaN value in real_psi_leg. Setting to zero.'
            f(i)=0.0
         end if

         result = result + weights(i)*f(i)
      end do

      real_psi_leg = result
      return
   end function real_psi_leg

end module phasespace
