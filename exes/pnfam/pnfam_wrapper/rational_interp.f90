!------------------------------------------------------------------------------
! rational_interp.f90
!
! This module contains subroutines to perform a rational function interpolation
! on a data set and evaluate that function at real or complex values.
!
! E.M. Ney UNC Chapel Hill, 2018
!------------------------------------------------------------------------------
module rational_interp
   implicit none
   integer, parameter, private :: dp = kind(1d0)

   ! A genereic interface for real or complex arguments to ratint_eval
   interface ratint_eval
      module procedure ratint_eval_real, ratint_eval_cmplx
   end interface ratint_eval

contains
   !----------------------------------------------------------------------------
   ! Compute the continuous fraction Thiele coefficients corresponding to a
   ! rational function passing through the supplied points (x,y). The order of the
   ! rational function numerator and denominator depends on 'npts'.
   !    npts.....: odd           even
   !    O(Num)...: ('npts'-1)/2  ('npts'/2)
   !    O(Den)...: (numerator)   (numerator-1)
   !
   ! - Inputs:  Arrays of length 'npts' for known 'x' and 'y' values of the function
   !            to be interpolated later.
   ! - Outputs: Thiele coefficients 'c' (there are npts coefficients)
   !----------------------------------------------------------------------------
   subroutine thiele_coeffs(npts, x, y, c)
      integer, intent(in) :: npts
      real(dp), intent(in) :: x(npts), y(npts)
      real(dp), intent(out) :: c(npts)
      real(dp) :: aux, den
      integer :: i,j,k
      do i=1, npts
         aux = y(i)
         j = 0
         do k=2, i
            j = j+1
            den = aux - c(j)
            if(den.eq.0.0) then
                write(*,'(A)') "Fatal: Division by zero attempted in Thiele continued fraction."
                stop
            end if
            aux = (x(i) - x(j))/den
         end do
         c(i) = aux
      end do
   end subroutine thiele_coeffs

   !----------------------------------------------------------------------------
   ! A wrapper for ratint_cmplx for to evaluate at real values
   !----------------------------------------------------------------------------
   subroutine ratint_eval_real(npts, x, c, z, y_int, dy_int)
      ! IO variables
      integer,  intent(in) :: npts
      real(dp), intent(in) :: x(npts), c(npts)
      real(dp), intent(in) :: z
      real(dp), intent(out):: y_int, dy_int
      ! Local variables
      complex(dp) :: zz
      complex(dp) :: yy_int, dyy_int

      ! Compute via complex routine
      zz = cmplx(z, 0.0, kind=dp)
      call ratint_eval_cmplx(npts, x, c, zz, yy_int, dyy_int)

      ! Check something weird didn't happend and our answer is real
      if (abs(aimag(yy_int)) > epsilon(1.0d0)) then
         write(*,'(/,"Warning in ratint_real: Im[result] != 0. Value = ",1es9.1,".")') aimag(yy_int)
      end if

      ! Save the real part
      y_int  = real(yy_int, kind=dp)
      dy_int = real(dyy_int, kind=dp)
   end subroutine ratint_eval_real

   !----------------------------------------------------------------------------
   ! Evaluate a rational function at a complex abscissa from its thiele coeffs
   ! by calculating the thiele continuous fraction.
   !
   ! - Inputs:  Arrays of length 'npts' for the thiele coefficients 'c' and the
   !            abscissae 'x' used to compute them, and the desired abscissa 'z'
   ! - Outputs: The interpolated value 'y_int' and an estimated error 'dy_int'
   !----------------------------------------------------------------------------
   subroutine ratint_eval_cmplx(npts, x, c, z, y_int, dy_int)
      ! IO variables
      integer,     intent(in) :: npts
      real(dp),    intent(in) :: x(npts), c(npts)
      complex(dp), intent(in) :: z
      complex(dp), intent(out):: y_int, dy_int
      ! Local variables
      real(dp)    :: TINY=1.0e-25, BIG=HUGE(0.0d0)
      complex(dp) :: s(npts), q(npts), aux
      complex(dp) :: y(npts), dy(npts-1)
      integer     :: k

      ! Initialize
      s(1) = c(1); q(1) = 1.0_dp; y(1) = s(1)
      s(2) = c(1)*c(2) + (z-x(1)); q(2) = c(2)
      if (q(2).eq.0.0) then
          write(*,*) "Fatal: Divide by zero attempted in ratint_eval."
          stop
      end if
      y(2) = s(2)/q(2)
      ! Calculate
      do k=3,npts
         aux = z - x(k-1)
         s(k) = c(k)*s(k-1) + aux*s(k-2)
         q(k) = c(k)*q(k-1) + aux*q(k-2)
         if (q(k).eq.0.0) stop
         y(k) = s(k)/q(k)
         if ( (abs(s(k)).gt.BIG .and. abs(q(k)).gt.BIG) .or. (abs(s(k)).lt.TINY.and.abs(q(k)).lt.TINY) ) then
            s(k-1) = s(k-1)/q(k)
            q(k-1) = q(k-1)/q(k)
            s(k) = y(k)
            q(k) = 1
         end if
      end do
      ! Get the "error" for each step
      do k=2,npts
         dy(k-1) = y(k) - y(k-1)
      end do
      ! Save the result
      y_int = y(npts)
      dy_int = dy(npts-1)
   end subroutine ratint_eval_cmplx

   !----------------------------------------------------------------------------
   ! Rational Function Interpolation via Bulirsch and Stoer algorithm, which
   ! uses diagonal rational functions. Source: Numerical Recipes.
   !
   ! - Inputs:  Arrays of length 'n' for known abscissae 'xa' and ordinates 'ya',
   !            and the desired abscissa 'x'
   ! - Outputs: interpolated value 'y' and error estimate 'dy'
   !----------------------------------------------------------------------------
   subroutine ratint(xa, ya, n, x, y, dy)
      integer, parameter  :: NMAX=10
      real(dp), parameter :: TINY=1.e-25
      integer   :: n, i, m, ns
      real(dp)  :: xa(n), ya(n)
      real(dp)  :: w, c(NMAX), d(NMAX), hh, h
      real(dp) :: dy, x ,y,t, dd

      ! Find the closest xa to desired x
      ns=1
      hh = abs(x-xa(1))
      do i=1,n
         h= abs(x-xa(i))
         if (h.eq.0) then
            y=ya(i)
            dy= 0.0_dp
            return
         else if (h.lt.hh) then
            ns=i
            hh=h
         end if
         c(i)=ya(i)
         ! TINY is needed to prevent a rare zero over zero condition
         d(i)=ya(i)+TINY
      end do
      y=ya(ns)
      ns=ns-1
      do m=1, n-1
         do i=1, n-m
            w=c(i+1)-d(i)
            h=xa(i+m)-x
            ! h is never zero (otherise x is in xa and we called return above)
            t=(xa(i)-x)*d(i)/h
            dd=t-c(i+1)
            ! This condition indicates the interpolating function as a pole at x
            if (dd.eq.0.0) then
               print *, 'failure in ratint'
               return
            end if
            dd=w/dd
            d(i)=c(i+1)*dd
            c(i)=t*dd
         end do
         if (2*ns .lt. n-m) then
            dy=c(ns+1)
         else
            dy=d(ns)
            ns=ns-1
         end if
         y=y+dy
      end do
      return
   end subroutine ratint

end module
