!------------------------------------------------------------------------------
! fermi.f90
!
! This module contains the Coulomb functions necessary for
! beta-decay rate computations.
!
! References for the future:
! [1] H. Behrens and W. Buhring, Electron Radial Wave Functions and Nuclear
!     Beta-decay (International Series of Monographs on Physics), Oxford
!     University Press, USA, 1982).
! [2] C. Lanczos, SIAP: Series B, Numerical Analysis 1, pp. 86 (1964).
!
! M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
!------------------------------------------------------------------------------
module fermi
   use constants, only : iu, pi, alpha, hbar_mec, r0
   
   implicit none
   public
   
   ! Kinds
   integer, parameter :: sp = kind(1.0)
   integer, parameter :: dp = kind(1.0d0)
   integer, parameter :: sc = kind(cmplx(1.0,1.0))
   integer, parameter :: dc = kind(cmplx(1.0,1.0,kind=dp))
      
   ! Screening parameter
   real(dp), parameter :: v0_tilde = 1.43_dp
   
   ! Generic interfaces
   interface L0_2_us
      module procedure L0_2_us_real, L0_2_us_cmplx
   end interface L0_2_us
   interface F0_us
      module procedure F0_us_real, F0_us_cmplx
   end interface F0_us
   interface w_tilde
      module procedure w_tilde_real, w_tilde_cmplx
   end interface w_tilde
   interface F0_sc
      module procedure F0_sc_real, F0_sc_cmplx
   end interface F0_sc
   interface F1_us
      module procedure F1_us_real, F1_us_cmplx
   end interface F1_us
   interface p2lambda2_us
      module procedure p2lambda2_us_real, p2lambda2_us_cmplx
   end interface p2lambda2_us
   interface p2lambda2_sc
      module procedure p2lambda2_sc_real, p2lambda2_sc_cmplx
   end interface p2lambda2_sc
   
   contains
   
   !----------------------------------------------------------------------------
   ! lanczos_gamma
   ! Lanczos approximation to the gamma function (\gamma = 5) (ref. [2]) in the
   ! domain Re(z) >= 1.
   !
   ! Parameter z is a double-precision complex number.
   !----------------------------------------------------------------------------
   function lanczos_gamma(z)
      implicit none
      complex(dc), intent(in) :: z
      complex(dc) :: lanczos_gamma
      complex(dc) :: zz, factor, series
      
      real(dp), dimension(7), parameter :: A5 = (/1.000000000178_dp, 76.180091729406_dp, &
                            -86.505320327112_dp, 24.014098222230_dp, -1.231739516140_dp, &
                              0.001208580030_dp, -0.000005363820_dp/)

      lanczos_gamma = 0
      zz = z - 1.0_dp
      
      ! Lanczos Gamma approximation is only defined for Re[z] >= 1
      if (real(z, kind=dp) < 1.0_dp) then
         write(*,'(/A)') 'ERROR in lanczos_gamma: Re[z] < 1.'
         stop
      end if
      
      factor = zz + 5.0_dp + 0.5_dp
      series = A5(1) + A5(2)/(zz+1.0_dp) + A5(3)/(zz+2.0_dp) + A5(4)/(zz+3.0_dp) &
                     + A5(5)/(zz+4.0_dp) + A5(6)/(zz+5.0_dp) + A5(7)/(zz+6.0_dp)   
   
      lanczos_gamma = factor**(zz+0.5_dp)*exp(-factor)*sqrt(2.0_dp*pi)*series
   end function lanczos_gamma
   
   !----------------------------------------------------------------------------
   ! lanczos_gamma_reflection
   ! Extend Gamma(z) to the left complex half-plane with the reflection theorem.
   ! See ref. [2], Eq. (4). This is done by testing the domain of z and either
   ! applying lanczos_gamma(z) directly or by using the reflection formula.
   !
   ! Parameter z is a double-precision complex number.
   !----------------------------------------------------------------------------
   function lanczos_gamma_reflection(z)
      implicit none
      complex(dc), intent(in) :: z
      complex(dc) :: lanczos_gamma_reflection
      
      if (real(z, kind=dp) >= 1.0_dp) then
         lanczos_gamma_reflection = lanczos_gamma(z)
      else
         lanczos_gamma_reflection = pi*(1-z)/sin(pi*(1-z))/lanczos_gamma(2-z)
      end if
   end function lanczos_gamma_reflection
   
   !----------------------------------------------------------------------------
   ! V0_shift
   ! Calculate the (positive) energy shift due to charge screening.
   !
   ! z is the charge number of the DAUGHTER nucleus.
   !----------------------------------------------------------------------------
   function V0_shift(z)
      implicit none
      integer, intent(in) :: z
      real(dp) :: V0_shift
      
      ! Throw and error for a positron calculation
      if (z < 0) then
         write(*,'(/A)') 'ERROR in V0_shift: Z < 0 unsupported.'
         stop
      end if
      
      V0_shift = v0_tilde*alpha**2*(z-1)**(4.0_dp/3.0_dp)
   end function V0_shift
   
   !----------------------------------------------------------------------------
   ! gamma_1
   ! Appears in the definition of the Fermi functions in e.g. [1].
   ! Defined in [1, eq. (4.8b)].
   !
   ! Here z is the charge number of the DAUGHTER nucleus.
   !----------------------------------------------------------------------------
   function gamma_1(z)
      implicit none
      integer, intent(in) :: z
      real(dp) :: gamma_1
      
      gamma_1 = sqrt(1.0_dp - (alpha*z)**2)
   end function gamma_1
   
   !----------------------------------------------------------------------------
   ! gamma_2
   ! Appears in the definition of the Fermi functions in e.g. [1, eq. (4.23)].
   ! Defined in the footnote of [1, pg. 105].
   !
   ! Here z is the charge number of the DAUGHTER nucleus.
   !----------------------------------------------------------------------------
   function gamma_2(z)
      implicit none
      integer, intent(in) :: z
      real(dp) :: gamma_2
      
      gamma_2 = sqrt(4.0_dp - (alpha*z)**2)
   end function gamma_2
   
   !----------------------------------------------------------------------------
   ! L0_us
   ! Simplest approximation for L0 without charge-screening [1, eq. (4.28)].
   !
   ! Parameter z is the charge number of the DAUGHTER nucleus.
   !----------------------------------------------------------------------------
   function L0_us(z)
      implicit none
      integer, intent(in) :: z
      real(dp) :: L0_us
      
      L0_us = 0.5_dp*(1.0_dp + gamma_1(z))
   end function L0_us
   
   !----------------------------------------------------------------------------
   ! L0_2_us_cmplx
   ! Approximation to L0 for |Z| <~ 15 in [1, eq. (4.126)], set up for double
   ! complex values of w.
   !
   ! z is the charge number of the DAUGHTER nucleus
   ! a is the mass number
   ! w is the electron total energy scaled by its mass-energy, mec2.
   !----------------------------------------------------------------------------
   function L0_2_us_cmplx(z, a, w)
      implicit none
      integer,  intent(in) :: z, a
      complex(dc), intent(in) :: w
      complex(dc) :: L0_2_us_cmplx
      real(dp)    :: R
      
      R = (r0/hbar_mec)*a**(1.0_dp/3.0_dp)
      L0_2_us_cmplx = 1 - alpha*z*w*R + (13.0_dp/60.0_dp)*(alpha*z)**2 - 0.5_dp*(alpha*z*R/w)
   end function L0_2_us_cmplx
   
   !----------------------------------------------------------------------------
   ! L0_2_us_real
   ! Approximation to L0 for |Z| <~ 15 in [1, eq. (4.126)], set up for real 
   ! values of w.
   !
   ! This routine calls L0_2_us_cmplx.
   !----------------------------------------------------------------------------
   function L0_2_us_real(z, a, w)
      implicit none
      integer,  intent(in) :: z, a
      real(dp), intent(in) :: w
      real(dp) :: L0_2_us_real
      complex(dc) :: aux
      
      aux = L0_2_us_cmplx(z=z, a=a, w=cmplx(w, 0, kind=dp))
      
      ! Test the complex part strictly
      if (abs(aimag(aux)) > epsilon(1.0d0)) then
         write(*,'(/,"Warning in L0_2_us_real: Im[result] != 0. Value = ",1es9.1,".")') aimag(aux)
      end if
      
      L0_2_us_real = real(aux, kind=dp)
   end function L0_2_us_real  
   
   !----------------------------------------------------------------------------
   ! F0_us_cmplx
   ! Fermi function F0 without charge-screening; see [1, eq. (4.8a)].
   ! NOTE: this does not include L0, so multiply that separately.
   ! NOTE: this is an updated version to accommodate complex energies w.
   !
   ! z is the charge number of the DAUGHTER nucleus.
   ! a is the mass number of the isotopes involved.
   ! w is the electron total energy scaled by its mass-energy, mec2.
   !
   ! NOTE: it is likely that double-precision is not enough as Re[y] grows
   ! large. This occurs for p small (w --> 1.0+). Double-precision can
   ! accommodate as far down as about w = 1.000 000 6 in my testing. With the
   ! old MQRPA code, I needed quad precision at times to avoid blowup. Maybe not
   ! necessary with FAM and its points at round numbers, but something to watch.
   !----------------------------------------------------------------------------
   function F0_us_cmplx(z, a, w)
      implicit none
      integer,     intent(in) :: z, a
      complex(dc), intent(in) :: w
      complex(dc) :: F0_us_cmplx
      complex(dc) :: p, y, expy, factor, gamma_arg
      real(dp)    :: R, g1, num, den
      
      F0_us_cmplx = 0
      
      ! This is not programmed for positron decay
      if (z < 0) then
         write(*,'(a)') 'ERROR in F0_us_cmplx: z < 0 not allowed.'
         stop
      end if
      
      ! Warn if Re(w) < 1
      if (real(w, kind=dp) < 1.0_dp) then
         write(*,'(/A)') 'WARNING: Re(w) < 1 in F0_us_cmplx'
      end if
      
      ! Special case: when p --> 0, the shape should go to 0 so we can safely return 0
      if (abs(w-1.0_dp) < epsilon(1.0_dp)) then
         F0_us_cmplx = 0
         return
      end if
      
      R = r0/hbar_mec*a**(1.0_dp/3.0_dp)
      p = sqrt(w**2 - 1)
      y = alpha*z*w/p
      g1 = gamma_1(z)
      
      factor = 4*(2*p*R)**(2*g1-2)
      expy = exp(pi*y)
      
      ! y is now complex, not necessarily real
      gamma_arg = g1 + iu*y
      
      ! Use the Fortran intrinsic for the denominator, since 2*g1 + 1 is real.
      ! Now using the all-z Lanczos approximation.
      num = abs(lanczos_gamma_reflection(gamma_arg))**2
      den = gamma(2*g1+1)**2
      
      F0_us_cmplx = factor*expy*num/den
   end function F0_us_cmplx
   
   !----------------------------------------------------------------------------
   ! F0_us_real
   ! Fermi function F0 without charge-screening; see [1, eq. (4.8a)].
   ! NOTE: this does not include L0, so multiply that separately.
   !
   ! This routine calls F0_us_cmplx.
   !----------------------------------------------------------------------------
   function F0_us_real(z, a, w)
      implicit none
      integer,  intent(in) :: z, a
      real(dp), intent(in) :: w
      real(dp)    :: F0_us_real
      complex(dc) :: f0_cmplx
      
      F0_us_real = 0
      
      ! This is not programmed for positron decay
      if (z < 0) then
         write(*,'(a)') 'ERROR in F0_us_real: z < 0 not allowed.'
         stop
      end if
      
      ! Warn if Re(w) < 1
      if (w < 1.0_dp) then
         write(*,'(/A)') 'WARNING: w < 1 in F0_us_real'
      end if
      
      ! Compute the complex version and assign the real part as output
      f0_cmplx = F0_us_cmplx(z, a, cmplx(w, 0, kind=dp))
      F0_us_real = real(f0_cmplx, kind=dp)
      
      ! Make sure we received a real result as a cross-check
      ! This is a very strict tolerance.
      if (abs(aimag(f0_cmplx)) > epsilon(1d0)) then
         write(*,'(/,"WARNING in F0_us_real: Im[result] != 0. Value = ",1es9.1,".")') aimag(f0_cmplx)
      end if
   end function F0_us_real
   
   !----------------------------------------------------------------------------
   ! w_tilde_cmplx
   ! Calculate the shifted electron energy due to charge screening. This version
   ! works for complex values of the energy.
   !
   ! z is the charge number of the DAUGHTER nucleus.
   ! w is the unshifted electron total energy scaled by its mass-energy, mec2.
   !----------------------------------------------------------------------------
   function w_tilde_cmplx(z, w)
      integer,     intent(in) :: z
      complex(dc), intent(in) :: w
      complex(dc) :: w_tilde_cmplx
      
      w_tilde_cmplx = w - V0_shift(z)
   end function w_tilde_cmplx

   !----------------------------------------------------------------------------
   ! w_tilde_real
   ! Calculate the shifted electron energy due to charge screening.
   !
   ! z is the charge number of the DAUGHTER nucleus.
   ! w is the unshifted electron total energy scaled by its mass-energy, mec2.
   !----------------------------------------------------------------------------
   function w_tilde_real(z, w)
      integer,  intent(in) :: z
      real(dp), intent(in) :: w
      real(dp)  :: w_tilde_real
      
      w_tilde_real = w - V0_shift(z)
   end function w_tilde_real

   !----------------------------------------------------------------------------
   ! F0_sc_cmplx
   ! Fermi function F0 with charge-screening; see [1, eq. (4.166)]. The screened
   ! version is, up to a multiplying factor, the unscreened version evaluated at
   ! a slightly different energy. NOTE: this does not include L0, so multiply
   ! that separately. This routine also multiplies the factors (p~/p) and
   ! (W~/W), so those are NOT required for L0 as well.
   ! NOTE: this routine is updated to handle complex energy w.
   !
   ! z is the charge number of the DAUGHTER nucleus.
   ! a is the mass number of the isotopes involved.
   ! w is the electron total energy scaled by its mass-energy, mec2.
   !
   ! NOTE: this routine inherits the blowup issues described for F0_us.
   !----------------------------------------------------------------------------
   function F0_sc_cmplx(z, a, w)
      implicit none
      integer,     intent(in) :: z, a
      complex(dc), intent(in) :: w
      complex(dc) :: F0_sc_cmplx, wt, p, pt
      
      F0_sc_cmplx = 0
      
      ! This is not programmed for positron decay
      if (z < 0) then
         write(*,'(a)') 'ERROR in F0_sc_cmplx: z < 0 not allowed.'
         stop
      end if
      
      ! Warn if Re(w) < 1
      if (real(w, kind=dp) < 1.0_dp) then
         write(*,'(/A)') 'WARNING: Re(w) < 1 in F0_sc_cmplx'
      end if      

      ! Shifted energy
      wt = w_tilde(z, w)
            
      ! Energy conservation
      if (real(wt, kind=dp) < 1.0_dp) then
         F0_sc_cmplx = 0
         return
      end if
      
      p  = sqrt(w**2 - 1)
      pt = sqrt(wt**2 - 1)
      
      ! F0_sc_cmplx = factor * F0_us(W-V0)
      F0_sc_cmplx = (pt/p)*(wt/w)*F0_us(z, a, wt)
   end function F0_sc_cmplx
   
   !----------------------------------------------------------------------------
   ! F0_sc_real
   ! Fermi function F0 with charge-screening; see [1, eq. (4.166)]. The screened
   ! version is, up to a multiplying factor, the unscreened version evaluated at
   ! a slightly different energy. NOTE: this does not include L0, so multiply
   ! that separately. This routine also multiplies the factors (p~/p) and
   ! (W~/W), so those are NOT required for L0 as well.
   !
   ! This routine calls F0_sc_cmplx.
   !----------------------------------------------------------------------------
   function F0_sc_real(z, a, w)
      implicit none
      integer,  intent(in) :: z, a
      real(dp), intent(in) :: w
      real(dp)    :: F0_sc_real
      complex(dc) :: f0_cmplx
      
      F0_sc_real = 0
      
      ! This is not programmed for positron decay
      if (z < 0) then
         write(*,'(a)') 'ERROR in F0_sc_real: z < 0 not allowed.'
         stop
      end if
      
      ! Warn if w < 1
      if (w < 1.0_dp) then
         write(*,'(/A)') 'WARNING: Re(w) < 1 in F0_sc_real'
      end if
      
      ! Call the complex version
      f0_cmplx = F0_sc_cmplx(z, a, cmplx(w, 0, kind=dp))
      F0_sc_real = real(f0_cmplx, kind=dp)
      
      ! Make sure we received a real result as a cross-check
      ! This is a very strict tolerance.
      if (abs(aimag(f0_cmplx)) > epsilon(1d0)) then
         write(*,'(/,"WARNING in F0_sc_real: Im[result] != 0. Value = ",1es9.1,".")') aimag(f0_cmplx)
      end if
   end function F0_sc_real
   
   !----------------------------------------------------------------------------
   ! F1_us_cmplx
   ! Generalized Fermi function F1 without charge-screening [1, eq. (4.23)].
   ! NOTE: this does not include L0 or other Coulomb functions, so multiply
   ! those separately. 
   ! NOTE: this routine is updated for complex energy w.
   !
   ! z is the charge number of the DAUGHTER nucleus.
   ! a is the mass number of the isotopes involved.
   ! w is the electron total energy scaled by its mass-energy, mec2.
   !
   ! NOTE: the same issues noted for F0_us_cmplx may occur in F1_us_cmplx.
   !----------------------------------------------------------------------------
   function F1_us_cmplx(z, a, w)
      implicit none
      integer,     intent(in) :: z, a
      complex(dc), intent(in) :: w
      complex(dc) :: F1_us_cmplx, p, y, expy, factor, gamma_arg
      real(dp) :: R, g2, num, den
      
      F1_us_cmplx = 0
      
      ! This is not programmed for positron decay
      if (z < 0) then
         write(*,'(a)') 'ERROR in F1_us_cmplx: z < 0 not allowed.'
         stop
      end if
      
      ! Warn if Re(w) < 1
      if (real(w, kind=dp) < 1.0_dp) then
         write(*,'(/A)') 'WARNING: Re(w) < 1 in F1_us_cmplx'
      end if
      
      ! Special case: when p --> 0, the shape should go to 0 so we can safely return 0
      if (abs(w-1.0_dp) < epsilon(1.0_dp)) then
         F1_us_cmplx = 0
         return
      end if
      
      R = r0/hbar_mec*a**(1.0_dp/3.0_dp)
      p = sqrt(w**2 - 1)
      y = alpha*z*w/p
      g2 = gamma_2(z)
      
      factor = 576*(2*p*R)**(2*g2-4)
      expy = exp(pi*y)
      
      ! y is now complex
      gamma_arg = g2 + iu*y
      
      ! Use the Fortran intrinsic for the denominator, since 2*g2 + 1 is real
      num = abs(lanczos_gamma_reflection(gamma_arg))**2
      den = gamma(2*g2+1)**2
   
      F1_us_cmplx = factor*expy*num/den
   end function F1_us_cmplx
   
   !----------------------------------------------------------------------------
   ! F1_us_real
   ! Generalized Fermi function F1 without charge-screening [1, eq. (4.23)].
   ! NOTE: this does not include L0 or other Coulomb functions, so multiply
   ! those separately. 
   ! 
   ! This routine calls F1_us_cmplx.
   !----------------------------------------------------------------------------
   function F1_us_real(z, a, w)
      implicit none
      integer,  intent(in) :: z, a
      real(dp), intent(in) :: w
      real(dp)    :: F1_us_real
      complex(dc) :: f1_cmplx
      
      F1_us_real = 0
      
      ! This is not programmed for positron decay
      if (z < 0) then
         write(*,'(a)') 'ERROR in F1_us_real: z < 0 not allowed.'
         stop
      end if
      
      ! Warn if w < 1
      if (w < 1.0_dp) then
         write(*,'(/A)') 'WARNING: w < 1 in F1_us_real'
      end if
      
      ! Handle with the complex routine
      f1_cmplx = F1_us_cmplx(z, a, cmplx(w, 0, kind=dp))
      F1_us_real = real(f1_cmplx, kind=dp)
      
      ! Make sure we received a real result as a cross-check
      ! This is a very strict tolerance.
      if (abs(aimag(f1_cmplx)) > epsilon(1d0)) then
         write(*,'(/,"WARNING in F1_us_real: Im[result] != 0. Value = ",1es9.1,".")') aimag(f1_cmplx)
      end if
   end function F1_us_real
   
   !----------------------------------------------------------------------------
   ! p2lambda2_us_cmplx
   ! First approximation to the Coulomb function p^2*lambda2 defined without
   ! charge-screening in [1, eq. (4.30)].
   ! NOTE: now updated for complex energy w.
   !
   ! z is the charge number of the DAUGHTER nucleus.
   ! a is the mass number of the isotopes involved.
   ! w is the electron total energy scaled by its mass-energy, mec2.
   !----------------------------------------------------------------------------
   function p2lambda2_us_cmplx(z, a, w)
      implicit none
      integer,     intent(in) :: z, a
      complex(dc), intent(in) :: w
      complex(dc) :: p2lambda2_us_cmplx
      complex(dc) :: f0, f1, p
      real(dp)    :: g1, g2
      
      p2lambda2_us_cmplx = 0
      
      ! This is not programmed for positron decay
      if (z < 0) then
         write(*,'(a)') 'ERROR in p2lambda2_us_cmplx: z < 0 not allowed.'
         stop
      end if
      
      ! Warn if Re(w) < 1
      if (real(w, kind=dp) < 1.0_dp) then
         write(*,'(/A)') 'WARNING: Re(w) < 1 in p2lambda2_us_cmplx'
      end if
      
      ! Special case: when p --> 0, the shape should go to 0 so we can safely return 0
      if (abs(w-1.0_dp) < epsilon(1.0_dp)) then
         p2lambda2_us_cmplx = 0
         return
      end if
      
      f0 = F0_us(z, a, w)
      f1 = F1_us(z, a, w)
      g1 = gamma_1(z)
      g2 = gamma_2(z)
      p  = sqrt(w**2 - 1)

      p2lambda2_us_cmplx = p**2*(f1/f0)*(2+g2)/(2+2*g1)
   end function p2lambda2_us_cmplx

   !----------------------------------------------------------------------------
   ! p2lambda2_us_real
   ! First approximation to the Coulomb function p^2*lambda2 defined without
   ! charge-screening in [1, eq. (4.30)].
   !
   ! This routine callse p2lambda2_us_cmplx
   !----------------------------------------------------------------------------
   function p2lambda2_us_real(z, a, w)
      implicit none
      integer,  intent(in) :: z, a
      real(dp), intent(in) :: w
      real(dp) :: p2lambda2_us_real
      complex(dc) :: p2l2_cmplx
      
      p2lambda2_us_real = 0
      
      ! This is not programmed for positron decay
      if (z < 0) then
         write(*,'(a)') 'ERROR in p2lambda2_us_real: z < 0 not allowed.'
         stop
      end if
      
      ! Warn if w < 1
      if (w < 1.0_dp) then
         write(*,'(/A)') 'WARNING: w < 1 in p2lambda2_us_real'
      end if
      
      ! Handle with the complex routine
      p2l2_cmplx = p2lambda2_us_cmplx(z, a, cmplx(w, 0, kind=dp))
      p2lambda2_us_real = real(p2l2_cmplx, kind=dp)
      
      ! Make sure we received a real result as a cross-check
      ! This is a very strict tolerance.
      if (abs(aimag(p2l2_cmplx)) > epsilon(1d0)) then
         write(*,'(/,"WARNING in p2lambda2_us_real: Im[result] != 0. Value = ",1es9.1,".")') aimag(p2l2_cmplx)
      end if
   end function p2lambda2_us_real

   !----------------------------------------------------------------------------
   ! p2lambda2_sc_cmplx
   ! First approximation to the Coulomb function p^2*lambda2 WITH charge
   ! screening. See [1, bottom of page 149]. There is some question of whether
   ! this is actually justified or not.
   ! NOTE: this routine updated for complex energy w.
   !
   ! For charge screening, make the change p^2*lambda2(W) -> (p~)^2*lambda(W~).
   !
   ! z is the charge number of the DAUGHTER nucleus.
   ! a is the mass number of the isotopes involved.
   ! w is the electron total energy scaled by its mass-energy, mec2.
   !----------------------------------------------------------------------------
   function p2lambda2_sc_cmplx(z, a, w)
      implicit none
      integer,     intent(in) :: z, a
      complex(dc), intent(in) :: w
      complex(dc) :: p2lambda2_sc_cmplx, wt
      
      p2lambda2_sc_cmplx = 0
      
      ! This is not programmed for positron decay
      if (z < 0) then
         write(*,'(a)') 'ERROR in p2lambda2_sc_cmplx: z < 0 not allowed.'
         stop
      end if
      
      ! Warn if Re(w) < 1
      if (real(w, kind=dp) < 1.0_dp) then
         write(*,'(/A)') 'WARNING: Re(w) < 1 in p2lambda2_sc_cmplx'
      end if
      
      ! Shifted energy
      wt = w_tilde(z, w)
      
      ! Energy conservation
      if (real(wt, kind=dp) < 1.0_dp) then
         p2lambda2_sc_cmplx = 0
         return
      end if
      
      ! Simply the unscreened version evaluated at W = W~.
      p2lambda2_sc_cmplx = p2lambda2_us_cmplx(z, a, wt)
   end function p2lambda2_sc_cmplx
   
   !----------------------------------------------------------------------------
   ! p2lambda2_sc_real
   ! First approximation to the Coulomb function p^2*lambda2 WITH charge
   ! screening. See [1, bottom of page 149]. There is some question of whether
   ! this is actually justified or not.
   !
   ! For charge screening, make the change p^2*lambda2(W) -> (p~)^2*lambda(W~).
   !
   ! z is the charge number of the DAUGHTER nucleus.
   ! a is the mass number of the isotopes involved.
   ! w is the electron total energy scaled by its mass-energy, mec2.
   !----------------------------------------------------------------------------
   function p2lambda2_sc_real(z, a, w)
      implicit none
      integer,  intent(in) :: z, a
      real(dp), intent(in) :: w
      real(dp) :: p2lambda2_sc_real, wt
      
      p2lambda2_sc_real = 0
      
      ! This is not programmed for positron decay
      if (z < 0) then
         write(*,'(a)') 'ERROR in p2lambda2_sc_real: z < 0 not allowed.'
         stop
      end if
      
      ! Warn if w < 1
      if (w < 1.0_dp) then
         write(*,'(/A)') 'WARNING: w < 1 in p2lambda2_sc_real'
      end if
      
      ! Shifted energy
      wt = w_tilde(z, w)
      
      ! Energy conservation
      if (wt < 1.0_dp) then
         p2lambda2_sc_real = 0
         return
      end if
      
      ! Simply the unscreened version evaluated at W = W~.
      p2lambda2_sc_real = p2lambda2_us_real(z, a, wt)
   end function p2lambda2_sc_real

end module fermi
