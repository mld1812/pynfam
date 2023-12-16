!------------------------------------------------------------------------------
! fermi_test.f90
!
! Test the Fermi module.
!
! M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
!------------------------------------------------------------------------------
program test_fermi

   use pnfam_tests
   use fermi

   implicit none
   
   real(tdp) :: lb_data(11,5,2)
   integer   :: npassed, nfailed
   
   npassed = 0
   nfailed = 0
   
   
   call gamma_test(npassed, nfailed)      ! Test Gamma function
   
   call init_landolt_bornstein
   call l0_us_test(npassed, nfailed)      ! Test unscreened Coulomb function L0
   call f0_us_test(npassed, nfailed)      ! Test unscreeneed Fermi function
   call lambda2_us_test(npassed, nfailed) ! Test p^2\lambda_2


   ! Summary
   write(*,*) 
   write(*,'(1x,a,/)') repeat('=',30)
   write(*,'(1x,"NO. PASSED: ",1i0)') npassed
   write(*,'(1x,"NO. FAILED: ",1i0)') nfailed
   
   if (nfailed == 0) then
      write(*,'(/,1x,3a)') fg_green, 'ALL TESTS PASSED!', cc_reset
   else
      write(*,'(/,1x,a,1i0,1x,3a)') fg_red, nfailed, &
         trim(merge('TEST ', 'TESTS', nfailed==1)), ' FAILED!!!', cc_reset
   end if
   
contains
   
   !----------------------------------------------------------------------------
   ! Test the gamma function against Mathematica 
   !----------------------------------------------------------------------------
   subroutine gamma_test(npassed, nfailed)
      implicit none
      
      integer, intent(inout) :: npassed, nfailed
      
      real(tdp), parameter :: tol = 1.0d-7
      
      real(tdp)    :: d1, input
      complex(tdp) :: dc1, dc2, inputc
      logical      :: pass
      integer      :: np0, nf0
      
      np0 = npassed
      nf0 = nfailed
      
      write(*,'(a)') 'TESTING LANCZOS APPROXIMATION TO GAMMA'
      write(*,'(a)') repeat('-',60)
      
      ! Gamma of a positive real number, double precision (random from Python)
      input = 9.820015072335234_tdp
      d1    = 242380.0_tdp
      dc1   = lanczos_gamma_reflection(cmplx(input, 0, kind=tdp))
      
      ! Dial tolerance down to about 5 digits
      call test_equal(real(dc1,kind=tdp), d1, tol*100*d1, 'Real part, GAMMA(x>0)', pass)
      if (pass) then
         npassed = npassed+1
      else
         nfailed = nfailed+1
      end if
      
      call test_equal(aimag(dc1), 0.0_dp, tol, 'Imag part, GAMMA(x>0)', pass)
      if (pass) then
         npassed = npassed+1
      else
         nfailed = nfailed+1
      end if
      
      ! Gamma of negative real number, double precision
      input = -0.13125734951754509_tdp
      d1    = -8.343848114466135270_tdp
      dc1   = lanczos_gamma_reflection(cmplx(input, 0, kind=tdp))
      
      call test_equal(real(dc1,kind=tdp), d1, tol, 'Real part, GAMMA(x<0)', pass)
      if (pass) then
         npassed = npassed+1
      else
         nfailed = nfailed+1
      end if
      call test_equal(aimag(dc1), 0.0_dp, tol, 'Imag part, GAMMA(x<0)', pass)
      if (pass) then
         npassed = npassed+1
      else
         nfailed = nfailed+1
      end if

      
      ! Gamma of a ++ complex number
      inputc = cmplx(6.723314931578964_tdp, 5.422529485263395_tdp, kind=tdp)
      dc1    = cmplx(-25.0703170963108514_tdp, -44.9521892231538128_tdp, kind=tdp)
      dc2    = lanczos_gamma_reflection(inputc)
      call test_equal(real(dc1,kind=tdp), real(dc2,kind=tdp), tol, 'Real part, GAMMA([+] + I[+])', pass)
      if (pass) then
         npassed = npassed+1
      else
         nfailed = nfailed+1
      end if
      call test_equal(aimag(dc1), aimag(dc2), tol, 'Imag part, GAMMA([+] + I[+])', pass)
      if (pass) then
         npassed = npassed+1
      else
         nfailed = nfailed+1
      end if

      ! Gamma of a +- complex number
      inputc = cmplx(6.723314931578964_tdp, -5.422529485263395_tdp, kind=tdp)
      dc1    = cmplx(-25.0703170963108514_tdp, 44.9521892231538128_tdp, kind=tdp)
      dc2    = lanczos_gamma_reflection(inputc)
      call test_equal(real(dc1,kind=tdp), real(dc2,kind=tdp), tol, 'Real part, GAMMA([+] + I[-])', pass)
      if (pass) then
         npassed = npassed+1
      else
         nfailed = nfailed+1
      end if
      call test_equal(aimag(dc1), aimag(dc2), tol, 'Imag part, GAMMA([+] + I[-])', pass)
      if (pass) then
         npassed = npassed+1
      else
         nfailed = nfailed+1
      end if

      ! Gamma of a -+ complex number
      inputc = cmplx(-6.723314931578964_tdp, 5.422529485263395_tdp, kind=tdp)
      dc1    = cmplx(2.8628249457700249d-10, 4.868176166114194d-10, kind=tdp)
      dc2    = lanczos_gamma_reflection(inputc)
      call test_equal(real(dc1,kind=tdp), real(dc2,kind=tdp), tol, 'Real part, GAMMA([-] + I[+])', pass)
      if (pass) then
         npassed = npassed+1
      else
         nfailed = nfailed+1
      end if
      call test_equal(aimag(dc1), aimag(dc2), tol, 'Imag part, GAMMA([-] + I[+])', pass)
      if (pass) then
         npassed = npassed+1
      else
         nfailed = nfailed+1
      end if

      ! Gamma of a -- complex number
      inputc = cmplx(-6.723314931578964_tdp, -5.422529485263395_tdp, kind=tdp)
      dc1    = cmplx(2.8628249457700249d-10, -4.868176166114194d-10, kind=tdp)
      dc2    = lanczos_gamma_reflection(inputc)
      call test_equal(real(dc1,kind=tdp), real(dc2,kind=tdp), tol, 'Real part, GAMMA([-] + I[-])', pass)
      if (pass) then
         npassed = npassed+1
      else
         nfailed = nfailed+1
      end if
      call test_equal(aimag(dc1), aimag(dc2), tol, 'Imag part, GAMMA([-] + I[-])', pass)
      if (pass) then
         npassed = npassed+1
      else
         nfailed = nfailed+1
      end if
      
      ! Summary
      if ((nfailed-nf0) == 0) then
         write(*,'(1x,3a)') fg_green, 'ALL TESTS PASSED!', cc_reset
      else
         write(*,'(1x,a,1i0,1x,3a)') fg_red, nfailed-nf0, &
            trim(merge('TEST ', 'TESTS', nfailed-nf0==1)), ' FAILED!!!', cc_reset
      end if

   end subroutine gamma_test
   
   !----------------------------------------------------------------------------
   ! Test F0 (unscreened) against Landolt Bornstein
   !----------------------------------------------------------------------------
   subroutine f0_us_test(npassed, nfailed)
      
      implicit none
      
      integer, intent(inout) :: npassed, nfailed
      
      integer   :: itest, ipt, z, a, np0, nf0
      real(tdp) :: f0_lb, f0_pnfam, p_m, percent
      
      np0 = npassed
      nf0 = nfailed
      
      write(*,'(/)') 
      write(*,'(1x,a)') 'TESTING FERMI FUNCTION F0 W/O SCREENING'
      
      do itest=1,2
         
         z = merge(32, 84, itest==1)
         a = merge(72,210, itest==1)
         
         write(*,'(1x,3(a,1i0))') 'Test ', itest, ': Z=', z, ', N=', a-z
         write(*,'(1x,a)') repeat('-',60)
         write(*,'(1x,1a5,2a12,1a11)') 'p/m', 'L0[L-B]', 'L0[PNFAM]', 'delta'
         
         do ipt=1, size(lb_data, 1)
            
            p_m      = lb_data(ipt,1,itest)
            f0_lb    = lb_data(ipt,4,itest)
            f0_pnfam = F0_us(z, a, sqrt(p_m**2+1))
            percent  = (f0_pnfam-f0_lb)/f0_lb*100
            
            ! CUTOFF = 0.1%
            if (abs(percent) < 0.1_dp) then
               write(*,'(1f6.1,2f12.6,1f9.2," %",3a)') p_m, f0_lb, f0_pnfam, percent, fg_green, '    OK', cc_reset
               npassed = npassed+1
            else
               write(*,'(1f6.1,2f12.6,1f9.2," %",3a)') p_m, f0_lb, f0_pnfam, percent, fg_red, '    INCORRECT', cc_reset
               nfailed = nfailed+1
            end if
         end do
         write(*,*) 
      end do
      
      ! Summary
      if ((nfailed-nf0) == 0) then
         write(*,'(1x,3a)') fg_green, 'ALL TESTS PASSED!', cc_reset
      else
         write(*,'(1x,a,1i0,1x,3a)') fg_red, nfailed-nf0, &
            trim(merge('TEST ', 'TESTS', nfailed-nf0==1)), ' FAILED!!!', cc_reset
      end if
      
   end subroutine f0_us_test

   !----------------------------------------------------------------------------
   ! Test L0 (unscreened) against Landolt Bornstein
   !----------------------------------------------------------------------------
   subroutine l0_us_test(npassed, nfailed)
      
      implicit none
      
      integer, intent(inout) :: npassed, nfailed
      
      integer   :: itest, ipt, z, a, np0, nf0
      real(tdp) :: l0_lb, l0_pnfam, p_m, percent
      
      np0 = npassed
      nf0 = nfailed
      
      write(*,'(/)') 
      write(*,'(1x,a)') 'TESTING COULOMB FUNCTION L0 W/O SCREENING'
      
      do itest=1,2
         
         z = merge(32, 84, itest==1)
         a = merge(72,210, itest==1)
         
         write(*,'(1x,3(a,1i0))') 'Test ', itest, ': Z=', z, ', N=', a-z
         write(*,'(1x,a)') repeat('-',60)
         write(*,'(1x,1a5,2a12,1a11)') 'p/m', 'F0[L-B]', 'F0[PNFAM]', 'delta'
         
         do ipt=1, size(lb_data, 1)
            p_m      = lb_data(ipt,1,itest)           ! p/m_e
            l0_lb    = lb_data(ipt,3,itest)           ! L0[LB]
            l0_pnfam = l0_2_us(z, a, sqrt(p_m**2+1))  ! L0[PNFAM]
            percent  = (l0_pnfam-l0_lb)/l0_lb*100
            
            if (abs(percent) < 25.0_dp) then
               write(*,'(1f6.1,2f12.6,1f9.2," %",3a)') p_m, l0_lb, l0_pnfam, percent, fg_green, '    OK', cc_reset
               npassed = npassed+1
            else
               write(*,'(1f6.1,2f12.6,1f9.2," %",3a)') p_m, l0_lb, l0_pnfam, percent, fg_red, '    INCORRECT', cc_reset
               nfailed = nfailed+1
            end if
         end do
         write(*,*) 
      end do
      
      ! Summary
      if ((nfailed-nf0) == 0) then
         write(*,'(1x,3a)') fg_green, 'ALL TESTS PASSED!', cc_reset
      else
         write(*,'(1x,a,1i0,1x,3a)') fg_red, nfailed-nf0, &
            trim(merge('TEST ', 'TESTS', nfailed-nf0==1)), ' FAILED!!!', cc_reset
      end if
      
   end subroutine l0_us_test
   
   !----------------------------------------------------------------------------
   ! Test \lambda_2 (unscreened) against Landolt Bornstein
   !----------------------------------------------------------------------------
   subroutine lambda2_us_test(npassed, nfailed)
      
      implicit none
      
      integer, intent(inout) :: npassed, nfailed
      
      integer   :: itest, ipt, z, a, np0, nf0
      real(tdp) :: lb, pnfam, p_m, percent
      
      np0 = npassed
      nf0 = nfailed
      
      write(*,'(/)') 
      write(*,'(1x,a)') 'TESTING COULOMB FUNCTION LAMBDA2 W/O SCREENING'
      
      do itest=1,2
         
         z = merge(32, 84, itest==1)
         a = merge(72,210, itest==1)
         
         write(*,'(1x,3(a,1i0))') 'Test ', itest, ': Z=', z, ', N=', a-z
         write(*,'(1x,a)') repeat('-',60)
         write(*,'(1x,1a5,2a12,1a11)') 'p/m', 'Lam2[L-B]', 'Lam2[PNFAM]', 'delta'
         
         do ipt=1, size(lb_data, 1)
            
            p_m    = lb_data(ipt,1,itest)
            lb     = lb_data(ipt,5,itest)
            pnfam  = p2lambda2_us(z, a, sqrt(p_m**2+1))/p_m**2
            percent = (pnfam-lb)/lb*100
            
            ! 35% CUTOFF
            if (abs(percent) < 35.0_dp) then
               write(*,'(1f6.1,2f12.6,1f9.2," %",3a)') p_m, lb, pnfam, percent, fg_green, '    OK', cc_reset
               npassed = npassed+1
            else
               write(*,'(1f6.1,2f12.6,1f9.2," %",3a)') p_m, lb, pnfam, percent, fg_red, '    INCORRECT', cc_reset
               nfailed = nfailed+1
            end if
         end do
         write(*,*) 
      end do
      
      ! Summary
      if ((nfailed-nf0) == 0) then
         write(*,'(1x,3a)') fg_green, 'ALL TESTS PASSED!', cc_reset
      else
         write(*,'(1x,a,1i0,1x,3a)') fg_red, nfailed-nf0, &
            trim(merge('TEST ', 'TESTS', nfailed-nf0==1)), ' FAILED!!!', cc_reset
      end if
      
   end subroutine lambda2_us_test
   
   !----------------------------------------------------------------------------
   ! Load LB test data
   ! 
   ! This test data is sampled from the Landolt-Bornstein database:
   ! Numerical Tables for Beta-Decay and Electron Capture (Springer-Verlag, 1969)
   ! H. Behrens, J. Janecke
   ! H. Schopper, ed.
   ! http://dx.doi.org/10.1007/b19939
   !----------------------------------------------------------------------------
   subroutine init_landolt_bornstein
      implicit none
      
      lb_data(:,:,:) = 0
      
      ! Point, (p/m, F0L0, L0, F0, Lambda2), Nucleus
      ! Set 1: Z=32, N=40, A=72
      lb_data( 1,:,1) = [ 0.1_tdp, 2.1749D+01, 1.0064_tdp, 2.1611D+01, 5.5523_tdp]
      lb_data( 2,:,1) = [ 0.5_tdp, 4.9124D+00, 1.0062_tdp, 4.8821D+00, 1.0796_tdp]
      lb_data( 3,:,1) = [ 1.0_tdp, 3.3163D+00, 1.0055_tdp, 3.2982D+00, 0.9524_tdp]
      lb_data( 4,:,1) = [ 2.0_tdp, 2.7375D+00, 1.0032_tdp, 2.7288D+00, 0.9344_tdp]
      lb_data( 5,:,1) = [ 3.0_tdp, 2.5788D+00, 1.0004_tdp, 2.5778D+00, 0.9396_tdp]
      lb_data( 6,:,1) = [ 4.0_tdp, 2.4983D+00, 0.9973_tdp, 2.5051D+00, 0.9464_tdp]
      lb_data( 7,:,1) = [ 5.0_tdp, 2.4451D+00, 0.9942_tdp, 2.4594D+00, 0.9531_tdp]
      lb_data( 8,:,1) = [10.0_tdp, 2.2962D+00, 0.9781_tdp, 2.3476D+00, 0.9801_tdp]
      lb_data( 9,:,1) = [20.0_tdp, 2.1333D+00, 0.9462_tdp, 2.2546D+00, 1.0199_tdp]
      lb_data(10,:,1) = [25.0_tdp, 2.0723D+00, 0.9307_tdp, 2.2266D+00, 1.0367_tdp]
      lb_data(11,:,1) = [30.0_tdp, 2.0181D+00, 0.9157_tdp, 2.2039D+00, 1.0525_tdp]
      
      ! Set 2: Z=84, N=126, A=210
      lb_data( 1,:,2) = [ 0.1_tdp, 3.8588D+02, 1.0111_tdp, 3.8164D+02, 15.6186_tdp]
      lb_data( 2,:,2) = [ 0.5_tdp, 8.1241D+01, 1.0095_tdp, 8.0476D+01,  1.1168_tdp]
      lb_data( 3,:,2) = [ 1.0_tdp, 4.5608D+01, 1.0049_tdp, 4.5386D+01,  0.6928_tdp]
      lb_data( 4,:,2) = [ 2.0_tdp, 2.8911D+01, 0.9902_tdp, 2.9197D+01,  0.6440_tdp]
      lb_data( 5,:,2) = [ 3.0_tdp, 2.3064D+01, 0.9729_tdp, 2.3706D+01,  0.6798_tdp]
      lb_data( 6,:,2) = [ 4.0_tdp, 1.9777D+01, 0.9551_tdp, 2.0707D+01,  0.7222_tdp]
      lb_data( 7,:,2) = [ 5.0_tdp, 1.7554D+01, 0.9374_tdp, 1.8726D+01,  0.7637_tdp]
      lb_data( 8,:,2) = [10.0_tdp, 1.1881D+01, 0.8569_tdp, 1.3865D+01,  0.9429_tdp]
      lb_data( 9,:,2) = [20.0_tdp, 7.5877D+00, 0.7338_tdp, 1.0340D+01,  1.2198_tdp]
      lb_data(10,:,2) = [25.0_tdp, 6.4652D+00, 0.6868_tdp, 9.4135D+00,  1.3340_tdp]
      lb_data(11,:,2) = [30.0_tdp, 5.6414D+00, 0.6471_tdp, 8.7180D+00,  1.4352_tdp]
      
   end subroutine init_landolt_bornstein
   
end program test_fermi
