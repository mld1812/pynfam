!------------------------------------------------------------------------------
! phasespace_test.f90
!
! Test the phasespace module.
!
! M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
!------------------------------------------------------------------------------
program test_phasespace

   use pnfam_tests
   use phasespace

   implicit none
   
   real(tdp), parameter :: wmax = 27.766732_tdp
   
   integer :: npassed, nfailed
   
   npassed = 0
   nfailed = 0
   
   call real_psi_test
   call polyfit_test
   
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
   ! Compare our real_psi routine against Mathematica for
   ! Z_final = 57, A = 174 
   !----------------------------------------------------------------------------
   subroutine real_psi_test
      implicit none

      real(tdp), parameter :: mathematica_values(6) = [2.209608199614178d5, &
         2.577925271039772d6, 3.470883516733893d7, 5.398097056657456d8,     &
         5.998661241258061d8, 5.178321508233396d8]

      integer :: np0, nf0, i
      real(tdp) :: result, delta
      
      np0 = npassed
      nf0 = nfailed
      
      write(*,'(1x,a,/,1x,a)') 'TEST OF REAL_PSI FOR Z_final=57, A=174, Wmax=27.766732, SI=1.0D-5', &
         repeat('-', 60)
      write(*,'(1x,1a5,3a20)') 'n', 'pnfam', 'mathematica', 'frac. change'

      do i=1, 6
         result = real_psi(57, 174, wmax, i, 1d-5)
         delta  = (result-mathematica_values(i))/mathematica_values(i)
         
         ! CUTOFF = 5E-4
         if (abs(delta) < 5d-4) then
            write(*,'(1i6,3es20.10,3a)') i, result, mathematica_values(i), delta, fg_green, '    OK', cc_reset
            npassed = npassed+1
         else
            write(*,'(1i6,3es20.10,3a)') i, result, mathematica_values(i), delta, fg_red, '    FAILED', cc_reset
            nfailed = nfailed+1
         end if
      end do
      write(*,*) 
      
      ! Summary
      if ((nfailed-nf0) == 0) then
         write(*,'(1x,3a)') fg_green, 'ALL TESTS PASSED!', cc_reset
      else
         write(*,'(1x,a,1i0,1x,3a)') fg_red, nfailed-nf0, &
            trim(merge('TEST ', 'TESTS', nfailed-nf0==1)), ' FAILED!!!', cc_reset
      end if

   end subroutine real_psi_test
   
   !----------------------------------------------------------------------------
   ! Test the polyfit routine
   !----------------------------------------------------------------------------
   subroutine polyfit_test
      
      use constants, only : iu, pi
      use complex_quadrature, only : cquad_simpson
      
      implicit none
      
      integer,   parameter :: order = 10
      integer,   parameter :: npts  = 200
      integer,   parameter :: ipts  = 55
      
      ! Energy values to evaluate the phase space integral
      real(tdp), parameter :: w(3) = [3.8843147869239845_tdp, 23.78515610154606_tdp, 8.366549403020539_tdp]
      
      ! Values of the total integral (sum of PSI(n; w)) from Mathematica
      real(tdp), parameter :: mathematica_values(6) = [1.793991047896786d5, &
         1.703821841693513d6, 1.874713489994236d7, 2.399812415210511d8,     &
         2.661499325243474d8, 2.256819886214595d8]
      
      integer      :: np0, nf0, i, ip, ip2
      real(tdp)    :: coeffs(order), int1, theta, t(ipts), result, delta, im, delta2
      complex(tdp) :: f(ipts), z, psi, int2, dz(ipts)
      
      np0 = npassed
      nf0 = nfailed
      
      write(*,'(/)') 
      write(*,'(1x,a)') 'TEST OF CONTOUR INTEGRATION FOR Z_final=57, A=174, Wmax=27.766732, SI=1.0D-5'
      write(*,'(1x,a)') 'Unit strength at w = [3.884, 8.367, 23.785]'
      write(*,'(1x,a)') 'Using production values (ORDER=10, NPTS=200, COUNTOUR POINTS=55)'
      write(*,'(1x,a)') repeat('-', 60)
      write(*,'(1x,1a5,5a20)') 'n', 'contour', 'real_psi', 'mathematica', 'f.c.(real_psi)', 'f.c.(mathematica)'
      
      do i=1, 6
         
         ! Function is unit strength at 3 values of w
         ! wmax = 1 + (eqrpamax/me), w = eqrpa/me
         int1 = 0
         do ip=1, 3
            int1 = int1 + real_psi(57, 174, wmax-w(ip), i, 1d-5)
         end do
         
         ! Approximate with the polyfit
         coeffs = get_polyfit_coeffs(57, 174, order, npts, wmax, i, 1d-5)
         
         do ip=1, ipts
            theta  = 2*pi/(ipts-1)*(ip-1)
            t(ip)  = theta
            
            z      = (wmax+1.0_tdp)/2.0_tdp + (wmax-1.0_tdp)/2.0_tdp*exp(iu*theta)
            dz(ip) = iu*(wmax-1.0_tdp)/2.0_tdp*exp(iu*theta)
            
            ! Evaluate at 1 + (EQRPAmax-E)/me
            psi = polyfit_eval(order, coeffs, wmax-z)
            
            ! Sum of three poles
            f(ip) = 0
            do ip2=1, 3
               f(ip) = f(ip) + psi/(z-w(ip2))
            end do
         end do
         
         ! Evaluation
         int2   = cquad_simpson(f, dz, t)
         im     = real(int2, kind=tdp)
         result = aimag(int2)/(2*pi)
         delta  = (result-int1)/int1
         delta2 = (result-mathematica_values(i))/mathematica_values(i)
         
         ! CUTOFF = 5E-3
         if (all([abs(delta),abs(delta2)] < 5d-3)) then
            write(*,'(1i6,5es20.10,3a)') i, result, int1, mathematica_values(i), delta, delta2, fg_green, '    OK', cc_reset
            npassed = npassed+1
         else
            write(*,'(1i6,5es20.10,3a)') i, result, int1, mathematica_values(i), delta, delta2, fg_red, '    FAILED', cc_reset
            nfailed = nfailed+1
         end if
      end do
      write(*,*) 
      
      ! Summary
      if ((nfailed-nf0) == 0) then
         write(*,'(1x,3a)') fg_green, 'ALL TESTS PASSED!', cc_reset
      else
         write(*,'(1x,a,1i0,1x,3a)') fg_red, nfailed-nf0, &
            trim(merge('TEST ', 'TESTS', nfailed-nf0==1)), ' FAILED!!!', cc_reset
      end if
      
   end subroutine polyfit_test
   
end program test_phasespace