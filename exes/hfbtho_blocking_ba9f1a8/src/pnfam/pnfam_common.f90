!------------------------------------------------------------------------------
! pnfam_common.f90
!
! Parts of the main program common to the serial and the parallel version.
! Mostly I/O.  Here to reduce code duplication.
!
! M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
!------------------------------------------------------------------------------
module pnfam_common
!$ use omp_lib, only : omp_get_max_threads, omp_get_wtime
   use iso_fortran_env,     only : error_unit, output_unit
   use constants,           only : kappa, element_symbols, pi, mec2, iu, ln2_kappa, ln2, &
                                   translate_uppercase, version
   use blockmatrix_type,    only : blockmatrix
   use external_field_type, only : external_field
   use broyden_mixer,       only : broyden_history_size, fam_si => si
   use hfbtho_basis,        only : ft_active, ft_temp, ft_fn, ft_fp, hfb_npr, pwi_qp_p, pwi_qp_n, &
                                   get_HFBTHO_solution, get_hfbtho_betadecay_properties
   use interaction,         only : init_skyrme, interaction_name, print_interaction_params
   use extfield,            only : init_external_field, print_crossterms, setup_crossterms
   use pnfam_solver,        only : moment_of_inertia
   use complex_quadrature,  only : cquad_simpson, cquad_trapezoidal
   use phasespace,          only : get_polyfit_coeffs, polyfit_eval
   use logger
   implicit none
   
   integer, parameter :: dp = kind(1d0)
   
   ! File number to avoid littering the code with hard-coded constants too much
   integer, parameter :: log_file_number = 32
   
   ! Namelist parameters not defined elsewhere ---
   ! GENERAL
   character(len=12)   :: fam_mode
   character(len=200) :: hfb_input_filename, fam_output_filename, log_filename

   ! EXT_FIELD
   character(len=80) :: operator_name
   integer :: k
   logical :: compute_crossterms

   ! SOLVER
   integer  :: max_iter
   real(dp) :: convergence_epsilon
   
   ! STR_PARAMETERS
   real(dp) :: energy_step

   ! CONTOUR_PARAMETERS
   integer  :: poly_order, poly_fit_points
   real(dp) :: q_value, rot_correction, contour_anchor, poly_fit_si
   logical  :: include_endpoint, rotate_contour
   
   ! FINDMAX_PARAMETERS
   integer  :: max_points
   real(dp) :: search_step, delta, energy_tolerance
   
   ! Shared parameters
   integer  :: nr_points
   real(dp) :: energy_start, half_width
   
   ! Namelist definitions ---
   namelist /general/ fam_mode, hfb_input_filename, fam_output_filename, log_filename
   namelist /ext_field/ operator_name, k, compute_crossterms
   namelist /str_parameters/ energy_start, energy_step, nr_points, half_width
   namelist /contour_parameters/ nr_points, q_value, rot_correction, include_endpoint, &
      rotate_contour, contour_anchor, poly_order, poly_fit_points, poly_fit_si
   namelist /findmax_parameters/ energy_start, half_width, delta, energy_tolerance, &
      search_step, max_points
   namelist /solver/ max_iter, broyden_history_size, convergence_epsilon
   
   ! Program variables ---
   ! HFB solutions
   integer  :: z_hfb, a_hfb, k_hfb, par_hfb
   real(dp) :: q_hfb, egs_hfb
   real(dp), allocatable :: Ep(:), En(:)
   type(blockmatrix) :: Up, Vp, Un, Vn
   
   ! External field and crossterms
   integer :: ixterms, nxterms
   type(external_field) :: f
   type(external_field), dimension(:), allocatable :: g
   
   ! Local program variables
   character(len=80) :: parameter_filename
   integer :: i, iter, fnml, fno, ierr, date_values(8), a, z_final, istr
   complex(dp) :: psi, cint
   character(len=160) :: st
   real(dp) :: t_start, t_finish, wt_start, wt_now, moi, omega
   real(dp),    dimension(:), allocatable :: str
   complex(dp), dimension(:), allocatable :: cstr
   
   ! Contour mode
   integer     :: nr_compute
   real(dp)    :: ctr_bound_left, ctr_bound_right, r0, r
   complex(dp) :: ci1, ci2, btot_simp, btot_trap, rate_simp, rate_trap, psi_simp, psi_trap
   real(dp),    allocatable :: ctr_t(:), poly_coeffs(:)
   complex(dp), allocatable :: ctr_z(:), ctr_dzdt(:), ctr_strength(:,:), ctr_psi(:)
   
   ! Findmax mode
   logical :: interval_found
   real(dp) :: deriv, y1, y2
   integer :: dir, firstdir
   
   ! OpenMP
   integer :: nthreads
   logical :: using_openmp
   
contains
   
   !---------------------------------------------------------------------------
   ! (Mostly) common initializations from the input files
   !---------------------------------------------------------------------------
   subroutine init_common ! except the interaction parameters
      implicit none
      
      wt_start = get_timer()
      
      
      ! Namelist initialization ---
      ! GENERAL
      fam_mode = ''; fam_output_filename = ''; hfb_input_filename = ''; log_filename = ''
      ! EXT_FIELD
       compute_crossterms = .false.; k = 0; operator_name = ''
      ! SOLVER
      broyden_history_size = 50; convergence_epsilon = 1.0d-7; max_iter = 200
      ! STR_PARAMETERS
      energy_start = 0; energy_step = 0.2_dp; half_width = 0.5_dp; nr_points = 101
      ! CONTOUR_PARAMETERS
      contour_anchor = -1.0_dp; include_endpoint = .false.; nr_points = 55
      poly_fit_points = 200; poly_fit_si = 1.0d-5; poly_order = 6
      q_value = -1.0_dp; rotate_contour = .false.; rot_correction = -1.0_dp
      ! FINDMAX_PARAMETERS
      delta = 0.002_dp; energy_start = 20.0_dp; energy_tolerance = 0.005_dp
      half_width = 5.0_dp; max_points = 20; search_step = 2.0_dp      

      fnml = 10
      open(fnml,file=parameter_filename)
      read(fnml,nml=general)
      read(fnml,nml=ext_field)
      call translate_uppercase(fam_mode) ! Let fam_mode be case-insensitive
      select case (fam_mode)
         case("STR")
            read(fnml,nml=str_parameters)
         case("CONTOUR")
            read(fnml,nml=contour_parameters)
         case("FINDMAX","FINDMAXDEF","FINDSDR","FINDSDRDEF")
            read(fnml,nml=findmax_parameters)
         case default
            write(error_unit,'(3a)') "ERROR: Unimplemented mode: ", trim(fam_mode), "!"
            call abort
      end select
      read(fnml,nml=solver)
      
      ! If a log file was specified, use that in place of stdout
      ! In the serial code, this is not very useful, as you could just pipe the
      ! output to a file, but in the parallel code we want to separate each
      ! nucleus/operator/K
      call openlog(log_filename)
      
      ! Read in the HFBTHO solution: qp energies, U and V matrices, equivalent s.p.
      ! energies, and the s.p. cutoff
      call get_HFBTHO_solution(fn=hfb_input_filename, Ep=Ep, Vp=Vp, Up=Up, En=En, Vn=Vn, Un=Un)
   
      ! Compute the external field in the single-particle basis
      call init_external_field(label=operator_name, k=k, op=f)
   
      ! If there are cross-terms, compute them here
      nxterms = 0
      if (compute_crossterms) then
         call setup_crossterms(op=f, crossterms=g, n=nxterms)
      end if
      ! The strength array is of dimension nxterms + 1 (nxterms + strength function)
      allocate(str(1:nxterms+1),cstr(1:nxterms+1))
      str(:) = 0; cstr(:) = 0
   
      ! Initialize the interaction (this is where the 'interaction' namelist is read)
      call init_skyrme(fh=fnml, vp=vp%elem, vn=vn%elem, up=up%elem, un=un%elem)
      close(fnml)
      
   end subroutine
   
   !----------------------------------------------------------------------------
   ! Read an input file and extract only the parameter fam_mode from the
   ! namelist GENERAL 
   !----------------------------------------------------------------------------
   function get_fam_mode(nmlfile) result(mode)
      implicit none
      character(len=*), intent(in) :: nmlfile
      character(len=len(nmlfile))  :: mode
      
      integer, parameter :: myfh = 45
      integer :: myierr
      ! Hacky to duplicate this definition, but I don't want to set the globals
      character(len=80)   :: fam_mode
      character(len=200) :: hfb_input_filename, fam_output_filename, log_filename
      namelist /general/ fam_mode, hfb_input_filename, fam_output_filename, log_filename
      
      mode = ""; fam_mode = ""
      open(myfh, file=trim(nmlfile), status='old', iostat=myierr)
      if (myierr /= 0) then
         write(error_unit,'(a)') 'Error setting expected_mode (opening file)'
         call abort
      end if
      read(myfh,nml=general,iostat=myierr)
      if (myierr /= 0) then
         write(error_unit,'(a)') 'Error setting expected_mode (reading NML)'
         call abort
      end if
      close(myfh)
      mode = adjustl(fam_mode)
      
   end function get_fam_mode
   
   
   !---------------------------------------------------------------------------
   ! Log that may go to standard output
   !---------------------------------------------------------------------------
   subroutine log_header
      implicit none
      character(len=160) :: st ! temporary strings for formatting
      
      ! Log the input parameters to standard output
      call writelog("")
      call writelog("Charge-changing Finite Amplitude Method (pnFAM)")
      
      ! Write the number of OpenMP threads if applicable
      nthreads = 0
      !$ nthreads = omp_get_max_threads()
      if (nthreads > 0) then
         write(st,'(3A)') "Version: ", version, " compiled with OpenMP" ; call writelog(st)
         write(st,'(A,1I0)') "OpenMP threads: ", nthreads ; call writelog(st)
      else
         write(st,'(2A)') "Version: ", version ; call writelog(st)
      end if
      call writelog("")
      call date_and_time(values=date_values)
      write(st,'("Run date: ",1i0,"/",1i0,"/",1i0,1x,1i0,":",1i2.2,":",1i2.2)') date_values(2), &
         date_values(3), date_values(1), date_values(5:7) ; call writelog(st)
      call writelog("")
      write(st,'(2A)') "Mode: ", fam_mode ; call writelog(st)
      call writelog("")
      write(st,'(3A)') "Input parameter file name: '", trim(parameter_filename), "'" ; call writelog(st)
      write(st,'(3A)') "HFB input data file name:  '", trim(hfb_input_filename), "'" ; call writelog(st)
      write(st,'(3A)') "FAM output data file name: '", trim(fam_output_filename), "'" ; call writelog(st)
      
      call writelog("")
      write(st,'(3A,I0)') "Operator: ", trim(operator_name), " for K = ", k ; call writelog(st)
   
      select case(fam_mode)
      
         case("STR")
            write(st,'(A,F6.4,A)') "Half-width gamma:      ", half_width, " MeV" ; call writelog(st)
            write(st,'(A,2(F0.4,A),I0,A)') "Energy range:          ", energy_start, " ... ", &
               energy_start+energy_step*(nr_points-1), " MeV (", nr_points, " points)" ; call writelog(st)
               
         case("CONTOUR")
            write(st,'(A,I0)')  "No. contour points:  ", nr_points ; call writelog(st)
         
            call writelog("")
            write(st,'(A,I0)')    "Maximum number of iterations:   ", max_iter ; call writelog(st)
            write(st,'(A,I0)')    "Broyden history size:           ", broyden_history_size ; call writelog(st)
            write(st,'(A,ES7.1)') "Convergence limit:              ", convergence_epsilon ; call writelog(st)
            write(st,'(A,1I0)')   "Phase-space polynomial order:   ", poly_order ; call writelog(st)
            write(st,'(A,1I0)')   "Polynomial fit points:          ", poly_fit_points ; call writelog(st)
            write(st,'(2A)')      "Include high-omega endpoint:    ", merge("Yes", "No ", include_endpoint) ; call writelog(st)
            write(st,'(2A)')      "Rotate contour to avoid poles:  ", merge("Yes", "No ", rotate_contour) ; call writelog(st)
            write(st,'(2A)')      "Manually adjust contour anchor: ", merge("Yes", "No ", contour_anchor >= 0) ; call writelog(st)
            if (contour_anchor >= 0) then
               write(st,'(A,1F0.4,A)') "Anchor point: ", contour_anchor, " MeV" ; call writelog(st)
            end if
            
         case("FINDMAX")
            write(st,'(A,F8.5,A)') "Initial guess for the peak:           ", energy_start, " MeV" ; call writelog(st)
            write(st,'(A,F8.5,A)') "Initial search step size:             ", search_step, " MeV" ; call writelog(st)
            write(st,'(A,I0)')      "Maximum number of derivatives:        ", max_points ; call writelog(st)
            write(st,'(A,F8.5,A)') "Half-width (gamma):                   ", half_width, " MeV" ; call writelog(st)
            write(st,'(A,ES7.1,A)') "Delta for estimating derivatives:     ", delta, " MeV" ; call writelog(st)
            write(st,'(A,ES7.1,A)') "Requested accuracy for peak location: ", energy_tolerance, " MeV" ; call writelog(st)
   
      end select
   
      call writelog("")
      write(st,'(A,I0)') "Maximum number of iterations: ", max_iter ; call writelog(st)
      write(st,'(A,I0)') "Broyden history size:         ", broyden_history_size ; call writelog(st)
      write(st,'(A,ES7.1)') "Convergence limit:            ", convergence_epsilon ; call writelog(st)
      call writelog("")
      write(st,'(2A)') 'Compute cross-terms:  ', merge('Yes', 'No ', compute_crossterms) ; call writelog(st)
      if (nxterms > 0) call print_crossterms(g%label)
      call writelog("")
      write(st,'(2A)') 'Finite-temperature active: ', merge('Yes', 'No ', ft_active) ; call writelog(st)
      
      ! Finite temperature
      if (ft_active) then
         write(st,'(A,1F8.4,A)')         'Temperature:      ', ft_temp, ' MeV' ; call writelog(st)
         write(st,'(A,1F8.4,A,1F0.4,A)') 'Approx. max. f_p: ', maxval(ft_fp), &
            ' (E_p = ', Ep(pwi_qp_p(maxloc(ft_fp))), ' MeV)' ; call writelog(st)
         write(st,'(A,1F8.4,A,1F0.4,A)') 'Approx. max. f_n: ', maxval(ft_fn), &
            ' (E_n = ', En(pwi_qp_n(maxloc(ft_fn))), ' MeV)' ; call writelog(st)
      end if
      
      call writelog("")
      call print_interaction_params
      call writelog("")
      
   end subroutine
   
   
   subroutine log_str_point
      implicit none
      character(len=160) :: st
      
      ! Print the result (strength function only!) to STDOUT w/statistics
      write(st,'(1F11.6,6X,1F16.6,6X,1I10,4X,1F8.2,4X,1F12.2)') omega, str(1), &
         iter, t_finish-t_start, (t_finish-t_start)/iter ; call writelog(st)

      ! At least let us know the convergence if there was a failure
      if (iter == -1) then
         write(st,'(A,1ES12.6)') 'FAM was interrupted at iteration limit. si = ', fam_si ; call writelog(st)
      end if
      
   end subroutine
   
   
   subroutine log_findmax_point
      implicit none
      character(len=160) :: st
      
      write(st,'(3(1F11.6,4X),I0)') omega, deriv, (y1+y2)/2, iter ; call writelog(st)
      
   end subroutine
   
   
   subroutine log_contour_point(z, b)
      implicit none
      complex(dp), intent(in) :: z, b
      character(len=160) :: st
      
      ! Print the result (strength function only!) to STDOUT w/statistics
      write(st,'(1F11.6,1X,1F11.6,1X,1F16.6,1X,1F16.6,1X,1I9,4X,1F12.2,5X,1F12.2)') &
         real(z,kind=dp), aimag(z), real(b,kind=dp),  aimag(b), &
         iter, t_finish-t_start, wt_now-wt_start
      call writelog(st)
      
      ! At least let us know the convergence if there was a failure
      if (iter == -1) then
         write(st,'(A,1ES12.6)') 'FAM was interrupted at iteration limit. si = ', fam_si
         call writelog(st)
      end if
      
   end subroutine
   
   
   !----------------------------------------------------------------------------
   ! Until we implement beta+ in the phase-space integrals, manually change the
   ! output for beta+ computations (e.g. need beta+ strength for a sum rule)
   !----------------------------------------------------------------------------
   subroutine log_contour_footer
      implicit none
      character(len=200) :: st

      call writelog(repeat('=',101))
      write(st,'(33x,a,11x,a,12x,a,8x,a)') "Simpson's 3/8","Trapezoidal","Delta","Percent"
      call writelog(st)
      call writelog(repeat('=',101))

      call contour_results
      
      ! Only available with beta-
      if (f%beta_minus) then
         write(st,'(2a)') 'Im[Phase space] ........ :', contour_compare(aimag(psi_simp),aimag(psi_trap))
         call writelog(st)
         write(st,'(2a)') 'Re[Contour integral] ... :', &
            contour_compare(-2/ln2_kappa*real(rate_simp,kind=dp),-2/ln2_kappa*real(rate_trap,kind=dp))
         call writelog(st)
      else
         write(st,'(a)') 'Im[Phase space] ........ :         N/A; Beta+'; call writelog(st)
         write(st,'(a)') 'Re[Contour integral] ... :         N/A; Beta+'; call writelog(st)
      end if
      
      call writelog(repeat('=',101))
      write(st,'(2a)') 'Transition strength .... :', &
         contour_compare(aimag(btot_simp),aimag(btot_trap)) ; call writelog(st)
      
      ! Only available with beta-
      if (f%beta_minus) then
         write(st,'(2a)') 'Total rate .........(/s) :', &
            contour_compare(aimag(rate_simp),aimag(rate_trap))
         call writelog(st)
         write(st,'(2a)') 'Half-life ...........(s) :', &
            contour_compare(ln2/aimag(rate_simp),ln2/aimag(rate_trap))
         call writelog(st)         
      else
         write(st,'(2a)') 'Total rate .........(/s) :         N/A; Beta+'; call writelog(st)
         write(st,'(2a)') 'Half-life ...........(s) :         N/A; Beta+'; call writelog(st)
      end if
      
      call writelog(repeat('=',101))
      
      call contour_tests
      
      wt_now = get_timer()
      call writelog("")
      call writelog("fin.")
      write(st,'(a,1f0.3,a)') 'wall-time was ', wt_now-wt_start, ' seconds.'
      call writelog(st)
      
   end subroutine
   
   
   subroutine log_findmax_footer
      implicit none
      character(len=160) :: st
      
      call writelog(repeat('-',50))
      call writelog("")
      
      if (search_step < energy_tolerance) then
         write(st,*) "Local energy maximum: ", omega ; call writelog(st)
      else
         write(st,*) "Energy maximum not located within ", max_points, " iterations." ; call writelog(st)
         call writelog("Please increase the number of iterations or improve the initial guess.")
      end if
      
   end subroutine
   
   
   subroutine log_str_footer
      implicit none
      
      ! Intentionally empty
      
   end subroutine
   
   
   !---------------------------------------------------------------------------
   ! Contour computation initializations:  Q-value, contour points and such
   !---------------------------------------------------------------------------
   subroutine init_contour_mode
      implicit none
      character(len=160) :: st
      logical :: lqnml, lanchornml, lrotnml
      
      ! HFB quantities 
      !-------------------------------------------------------------------------
      call writelog("")
      call get_HFBTHO_betadecay_properties(Ep=Ep, En=En, Vp=Vp, Vn=Vn, Q=q_hfb, &
         E_gs=egs_hfb, K_gs=k_hfb, P_gs=par_hfb, lpr=.true.)
      
      z_hfb = hfb_npr(2); a_hfb = hfb_npr(3); a = a_hfb
      if (f%beta_minus) then
         z_final = z_hfb + 1
      else
         z_final = z_hfb  - 1
      end if
      
      call writelog("")
      call writelog('CALCULATION PARAMETERS:')
      call writelog(repeat('-',30))
      write(st,'(2(A,I0,A))') "Mother nucleus:   ", A, element_symbols(z_hfb) ; call writelog(st)
      write(st,'(2(A,I0,A))') "Daughter nucleus: ", A, element_symbols(Z_final) ; call writelog(st)
      
      ! Contour energy interval
      !-------------------------------------------------------------------------
      lqnml = (q_value >= 0)
      lanchornml = (contour_anchor >= 0)
      ! Right bound --- q_value is set in the namelist
      ctr_bound_right = 0
      if (lqnml) then
         ctr_bound_right = egs_hfb + q_value         
      else
         ctr_bound_right = egs_hfb + q_hfb
      end if
      ! Left bound --- contour_anchor is set in the namelist
      ctr_bound_left = 0
      if (lanchornml) then
         ctr_bound_left = contour_anchor
      else
         ! Finite temperature ---
         if (ft_active) ctr_bound_left = 0.5_dp*egs_hfb
      end if
      
      call writelog("")
      if (lqnml) then
         write(st,'("Effective Q-value: ",1F9.6," MeV  (set in namelist CONTOUR_PARAMETERS)")') &
            ctr_bound_right-egs_hfb
         call writelog(st)
      else
         write(st,'("Effective Q-value: ",1F9.6," MeV  (computed from HFB solution)")') &
            ctr_bound_right-egs_hfb
         call writelog(st)
      end if
      if (lanchornml) then
         write(st,'("Contour L bound:   ",1F9.6," MeV  (set in namelist CONTOUR_PARAMETERS)")') &
            ctr_bound_left
         call writelog(st)
      else
         if (ft_active) then
            write(st,'("Contour L bound:   ",1F9.6," MeV  (set to 1/2 of HFB g.s. energy)")') &
               ctr_bound_left
            call writelog(st)
         else
            write(st,'("Contour L bound:   ",1F9.6," MeV")') ctr_bound_left
            call writelog(st)
         end if
      end if
      write(st,'("Contour R bound:   ",1F9.6," MeV")') ctr_bound_right
      call writelog(st)
      
      if (ctr_bound_left > egs_hfb) then
         call writelog('WARNING: EQRPA_min > E(g.s).')
      end if
      if (ft_active .and. egs_hfb < 1.0_dp) then
         call writelog('WARNING: E(g.s.) < 1 MeV. Treat the (0,0) pole with care.')
      end if
      call writelog("")
      
      ! Integration contour setup
      call contour_setup

      ! Phase space integral
      allocate(poly_coeffs(poly_order))
      poly_coeffs(:) = 0
      
      ! Beta+ is not considered yet --- only set up the polynomial fit if in beta- mode
      if (f%beta_minus) then
         poly_coeffs = get_polyfit_coeffs(nz=z_final, na=a, order=poly_order, nfi=2, &
            npts=poly_fit_points, w0max=1+(ctr_bound_right-ctr_bound_left)/mec2, si=poly_fit_si)
      else
         call writelog('WARNING: beta+ computations are not available in CONTOUR mode. &
            & As a result, all rates will be zero.')
      end if
      
      ! Rotational energy correction --- rot_correction is set in the namelist
      lrotnml = (rot_correction >= 0)
      if (.not.lrotnml) then
         moi = moment_of_inertia(ep, up, vp) + moment_of_inertia(en, un, vn)
         rot_correction = 0.5d0/moi*f%rank*(f%rank + 1)         
      end if
      
      if (.not.lrotnml) then
         write(st,'("Moment of inertia: ",1F9.6," MeV^-1")') moi
         call writelog(st)
      else
         write(st,'("Moment of inertia: N/A")')
         call writelog(st)
      end if
      
      if (.not.lrotnml) then
         write(st,'("Rotational energy correction: ",1F9.6," MeV  (computed automatically)")') &
            rot_correction
         call writelog(st)
      else
         write(st,'("Rotational energy correction: ",1F9.6," MeV  (set manually)")') &
            rot_correction
         call writelog(st)
      end if
      
      if (rot_correction > 1.0_dp) then
         call writelog("*** THIS SEEMS HIGH, PLEASE CHECK ***")
      end if
      
      ! Header for STDOUT
      call writelog("")
      write(st,'(A,I0)') "Number of non-trivial density matrix elements: ", size(f%mat%elem) ; call writelog(st)
      call writelog("")
      call writelog("       Energy [MeV]                  Strength [MeV^-1]  &
                      &       Iter.        Time (s)        Wtime (s)")
      call writelog("     real        imag             real             imag")
      call writelog(repeat('-',101))
         
   end subroutine
   
   
   subroutine init_str_mode
      implicit none
      character(len=160) :: st
      
      write(st,'(A,I0)') "Number of non-trivial density matrix elements: ", size(f%mat%elem) ; call writelog(st)
      call writelog("")
      call writelog("Energy [MeV]     Strength [MeV^-1]     Iterations     Time (s)     T/Iter. (s)")
      call writelog("------------------------------------------------------------------------------")
      
   end subroutine
   
   
   subroutine init_findmax_mode
      implicit none
      character(len=160) :: st
      
      write(st,'(A,I0)') "Number of non-trivial density matrix elements: ", size(f%mat%elem)
      call writelog(st)
      call writelog("")
      call writelog("    Energy      derivative       strength  ~iter")
      call writelog(repeat('-',50))
      
   end subroutine
   
   
   !---------------------------------------------------------------------------
   ! ASCII output file
   !---------------------------------------------------------------------------
   subroutine write_output_file_header
      implicit none
      
      fno = 11
      open(fno,file=fam_output_filename)
      write(fno,'(2A)')             "# pnFAM code version: ", version
      write(fno,'(2A)')             "# Residual interaction: ", interaction_name
      write(fno,'(3A,I0)')          "# Operator: ", trim(operator_name), " with K = ", k
      
      select case(fam_mode)
      case("STR")
         write(fno,'(A,F0.4)')         "# Gamma (half-width): ", half_width
         write(fno,'(A)')              "#"
         write(fno,'(A)',advance='no') "# Energy [MeV]          Strength [MeV^-1]"
         ! Add extra headers for cross-terms
         do istr=2, size(str)
            write(fno,'(4x,a23)',advance='no') trim(g(istr-1)%label)//" [MeV^-1]"
         end do
         write(fno,'(A)') ""
      
      case("CONTOUR")
         write(fno,'(A)') "#"
         write(fno,'(6A)',advance='no') "#", txtr("Parameter", 29), txtr("Re(EQRPA)", 30), &
            txtr("Im(EQRPA)", 30), txtr("Re(Strength)", 34), txtr("Im(Strength)", 34)
         do ixterms=2, size(str)
           write(fno,'(2A)',advance='no') txtr("Re("//trim(g(ixterms-1)%label)//")", 34), &
              txtr("Im("//trim(g(ixterms-1)%label)//")", 34)
         end do
         write(fno,'(2A)') txtr("Re(f(E))", 34), txtr("Im(f(E))", 34)
         
      end select
   end subroutine
   
   
   subroutine write_output_str_point
      implicit none
      ! Write the full result to the output file, including cross-terms
      write(fno,'(1F14.6)',advance='no') omega
      do istr=1, size(str)
         write(fno,'(4X,1ES23.16)',advance='no') str(istr)
      end do
      write(fno,'(A)') ""
      flush(fno)
   end subroutine
   
   
   subroutine write_output_contour_point(ip)
      implicit none
      integer, intent(in) :: ip
      integer :: i
      
      ! Write the sampled points in the output file for diagnostic purposes
      write(fno,'(3f30.19)',advance='no') ctr_t(ip), ctr_z(ip)
      
      ! Write the strengths
      do i=1, size(cstr)
         write(fno,'(2ES34.19)',advance='no') ctr_strength(ip,i)
      end do
      
      ! Write the phase-space integral value here
      write(fno,'(2ES34.19)') ctr_psi(ip)
      
   end subroutine
   
   
   !----------------------------------------------------------------------------
   ! Until we implement beta+ in the phase-space integrals, manually change the
   ! output for beta+ computations (e.g. need beta+ strength for a sum rule)
   !----------------------------------------------------------------------------
   subroutine write_output_contour_footer
      implicit none
      
      write(fno,*) 
      write(fno,'("# NR_POINTS: ",1i0)') nr_points
      write(fno,'("# NR_XTERMS: ",1i0)') nxterms
      write(fno,'("# ",a)') repeat('-',98)
      write(fno,'("# ",a,27x,a,16x,a)') 'Quantity', "Simpson's Rule", "Trapezoidal Rule"
      write(fno,'("# ",a)') repeat('-',98)
      write(fno,'("# Transition strength ... : ",2ES30.15)') aimag(btot_simp),     aimag(btot_trap)
      
      if (f%beta_minus) then
         write(fno,'("# Decay rate ..... (s^-1) : ",2ES30.15)') aimag(rate_simp),     aimag(rate_trap)
         write(fno,'("# Half-life ..........(s) : ",2ES30.15)') ln2/aimag(rate_simp), ln2/aimag(rate_trap)
      else
         write(fno,'("# Decay rate ..... (s^-1) :          N/A; Beta+")')
         write(fno,'("# Half-life ..........(s) :          N/A; Beta+")')
      end if
      
      close(fno) ! close the output file
      
   end subroutine
   
   
   subroutine write_output_findmax_result
      implicit none
      
      write(fno,'(A/A)') '#', '# Energy'
      if (search_step < energy_tolerance) then
         write(fno,'(1F0.16)') omega
      else
         write(fno,'(1I8)') -1
         write(fno,'(A,1I0,A)') "# Maximum not located within ", max_points, &
            merge(" iterations.", " iteration. ", max_points > 1)
      end if
      close(fno) ! close the output file
      
   end subroutine
   
   
   !----------------------------------------------------------------------------
   ! Set up a circular contour
   !
   ! By default, the contour is anchored with the first point at
   ! z = (ctr_bound_right, 0). If rotate_contour is .TRUE., then the contour is
   ! spun so that the first point is dt/2 radians above the real line.
   !
   ! Practically, this means:
   ! nr_points even: one point on the real axis; the rightmost point if
   !                 rotate_contour = F, the leftmost if rotate_contour = T
   ! nr_points odd:  either two or zero points on the real axis
   !----------------------------------------------------------------------------
   subroutine contour_setup
      implicit none
      integer  :: icp, nv
      real(dp) :: t0
      
      allocate(ctr_z(nr_points), ctr_t(nr_points), ctr_dzdt(nr_points))
      allocate(ctr_psi(nr_points), ctr_strength(nr_points,size(cstr)))
      ctr_z = 0; ctr_t = 0; ctr_dzdt = 0; ctr_psi = 0; ctr_strength = 0
      
      r0 = (ctr_bound_right+ctr_bound_left)/2
      r  = ctr_bound_right-r0
      nv = nr_points-1
      t0 = 0
      
      ! Spin the contour away from the real line
      if (rotate_contour) t0 = pi/nv
      
      do icp=1, nr_points
         ctr_t(icp) = t0 + (icp-1)*2*pi/nv
      end do
      
      ctr_z(:) = r0 + r*exp(iu*ctr_t(:))
      ctr_z(nr_points) = ctr_z(1)
      ctr_dzdt(:) = iu*(ctr_z(:) - r0)

      ! The number of points to compute depends on the orientation, as do the
      ! symmetry properties of the contour. Integer division gets this right.
      if (rotate_contour) then
         nr_compute = nr_points/2
      else
         nr_compute = (nr_points+1)/2
      end if
      
   end subroutine contour_setup
   
   !----------------------------------------------------------------------------
   ! Use symmetry to fill in the contour, strength function, and phase space  
   !----------------------------------------------------------------------------
   subroutine contour_finalize_computation
      implicit none
      integer     :: ip, iflip
      complex(dp) :: ener, psint
      
      ctr_strength(nr_points,:) = ctr_strength(1,:)
      do ip=nr_compute+1, nr_points-1
         if (rotate_contour) then
            iflip = nr_points-ip
         else
            iflip = nr_points-ip+1
         end if
         ctr_strength(ip,:) = conjg(ctr_strength(iflip,:))
      end do
      
      do ip=1, nr_points
         ener  = 1.0_dp+(ctr_bound_right-ctr_z(ip))/mec2
         psint = polyfit_eval(order=poly_order, cs=poly_coeffs(:), z=ener)
         ctr_psi(ip) = psint
         
         call write_output_contour_point(ip)
      end do
      
   end subroutine contour_finalize_computation
   
   !----------------------------------------------------------------------------
   ! Calculate various numerical results 
   !----------------------------------------------------------------------------
   subroutine contour_results
      implicit none
      
      real(dp), parameter :: half = 0.5_dp
      
      ! Total strength in the contour
      btot_trap = cquad_trapezoidal(f=ctr_strength(:,1), dz=ctr_dzdt, t=ctr_t)
      btot_simp = cquad_simpson(f=ctr_strength(:,1), dz=ctr_dzdt, t=ctr_t)
      ! Phase space only
      psi_trap = cquad_trapezoidal(f=ctr_psi, dz=ctr_dzdt, t=ctr_t)
      psi_simp = cquad_simpson(f=ctr_psi, dz=ctr_dzdt, t=ctr_t)
      ! Total rate
      rate_trap = cquad_trapezoidal(f=ctr_strength(:,1)*ctr_psi, dz=ctr_dzdt, t=ctr_t)
      rate_simp = cquad_simpson(f=ctr_strength(:,1)*ctr_psi, dz=ctr_dzdt, t=ctr_t)
      
      ! Constants
      ! Now, B = Im(Btot) and Rate = Im(Rate)
      btot_trap = -half*btot_trap
      btot_simp = -half*btot_simp
      psi_trap  = -half*psi_trap
      psi_simp  = -half*psi_simp
      rate_trap = -half*ln2_kappa*rate_trap
      rate_simp = -half*ln2_kappa*rate_simp
      
   end subroutine contour_results
   
   !----------------------------------------------------------------------------
   ! Compare two real quantities 
   !----------------------------------------------------------------------------
   function contour_compare(a, b, fmt_exp)
      implicit none
      real(dp), intent(in) :: a, b
      logical, optional, intent(in) :: fmt_exp
      character(len=100) :: fmt
      character(len=160) :: contour_compare
      
      fmt = '(1f19.9,1f23.9,1es20.4,1f11.2,"%")'
      if (present(fmt_exp)) then
         if (fmt_exp) fmt = '(2es23.9,1es16.4,1f11.2,"%")'
      end if
      write(contour_compare,trim(fmt)) a, b, b-a, 100*(b-a)/a

   end function contour_compare
   
   !----------------------------------------------------------------------------
   ! Check various numerical properties of the integrals 
   !----------------------------------------------------------------------------
   subroutine contour_tests
      implicit none
      real(dp), parameter :: tol = 1.0d-4
      real(dp) :: simp, trap
      logical  :: lsimp, ltrap
      
      ! Strength function --- 
      ! Residues are all real, so the real part of B_tot should vanish!
      trap = real(btot_trap,kind=dp); ltrap = (abs(trap) > tol)
      simp = real(btot_simp,kind=dp); lsimp = (abs(simp) > tol)
      if (lsimp .or. ltrap) then
         call writelog("")
         call writelog(txtc(repeat('*',95),101))
         if (lsimp .and. (.not.ltrap)) then
            call writelog(txtc("NOTE: Re[B] is non-zero with Simpson's&
               & rule but seems OK with Trapezoidal.",101))
         else if (ltrap .and. (.not.lsimp)) then
            call writelog(txtc("NOTE: Re[B] is non-zero with Trapezoidal&
               & rule but seems OK with Simpson's rule.",101))
         else
            call writelog(txtc('WARNING: Re[B] is non-zero with both&
               & quadratures. CHECK VERY CAREFULLY.',101))
         end if
         call writelog(txtc(repeat('*',95),101))
      end if

      ! Phase space ---
      ! This should be a closed integral of an entire function and evaluate to 0
      ! Real part
      trap = real(psi_trap,kind=dp); ltrap = (abs(trap) > tol)
      simp = real(psi_simp,kind=dp); lsimp = (abs(simp) > tol)
      if (lsimp .or. ltrap) then
         call writelog(txtc(repeat('*',95),101))
         if (lsimp .and. (.not.ltrap)) then
            call writelog(txtc("NOTE: Re[Phase space] is non-zero with Simpson's&
               & rule but seems OK with Trapezoidal.",101))
         else if (ltrap .and. (.not.lsimp)) then
            call writelog(txtc("NOTE: Re[Phase space] is non-zero with Trapezoidal&
               & rule but seems OK with Simpson's rule.",101))
         else
            call writelog(txtc('WARNING: Re[Phase space] is non-zero with both&
               & quadratures. CHECK VERY CAREFULLY.',101))
         end if
         call writelog(txtc(repeat('*',95),101))
      end if
      ! Imaginary part
      trap = aimag(psi_trap); ltrap = (abs(trap) > tol)
      simp = aimag(psi_simp); lsimp = (abs(simp) > tol)
      if (lsimp .or. ltrap) then
         call writelog(txtc(repeat('*',95),101))
         if (lsimp .and. (.not.ltrap)) then
            call writelog(txtc("NOTE: Im[Phase space] is non-zero with Simpson's&
               & rule but seems OK with Trapezoidal.",101))
         else if (ltrap .and. (.not.lsimp)) then
            call writelog(txtc("NOTE: Im[Phase space] is non-zero with Trapezoidal&
               & rule but seems OK with Simpson's rule.",101))
         else
            call writelog(txtc('WARNING: Im[Phase space] is non-zero with both&
               & quadratures. CHECK VERY CAREFULLY.',101))
         end if
         call writelog(txtc(repeat('*',95),101))
      end if

      
      ! Real part of contour ---
      trap = -2/ln2_kappa*real(rate_trap,kind=dp); ltrap = (abs(trap) > tol)
      simp = -2/ln2_kappa*real(rate_simp,kind=dp); lsimp = (abs(simp) > tol)
      if (lsimp .or. ltrap) then
         call writelog(txtc(repeat('*',95),101))
         if (lsimp .and. (.not.ltrap)) then
            call writelog(txtc("NOTE: Re[Contour integral] is non-zero with Simpson's&
               & rule but seems OK with Trapezoidal.",101))
         else if (ltrap .and. (.not.lsimp)) then
            call writelog(txtc("NOTE: Re[Contour integral] is non-zero with Trapezoidal&
               & rule but seems OK with Simpson's rule.",101))
         else
            call writelog(txtc('WARNING: Re[Contour integral] is non-zero with both&
               & quadratures. CHECK VERY CAREFULLY.',101))
         end if
         call writelog(txtc(repeat('*',95),101))
      end if
      
      ! Sign of decay rate
      trap = aimag(rate_trap); ltrap = (trap < 0)
      simp = aimag(rate_simp); lsimp = (simp < 0)
      if (lsimp .or. ltrap) then
         call writelog(txtc(repeat('*',95),101))
         if (lsimp .and. (.not.ltrap)) then
            call writelog(txtc("NOTE: Decay rate is negative with Simpson's&
               & rule but seems OK with Trapezoidal.",101))
         else if (ltrap .and. (.not.lsimp)) then
            call writelog(txtc("NOTE: Decay rate is negative with Trapezoidal&
               & rule but seems OK with Simpson's rule.",101))
         else
            call writelog(txtc('WARNING: Decay rate is negative with both&
               & quadratures. CHECK VERY CAREFULLY.',101))
         end if
         call writelog(txtc(repeat('*',95),101))
      end if

   end subroutine contour_tests
   
   !---------------------------------------------------------------------------
   ! Binary output file for the further processing of contour calculation
   !---------------------------------------------------------------------------
   subroutine contour_binary_output
      implicit none
      integer, parameter :: myfh = 44
      
      ! Store the output data to a special binary file
      !
      ! FORMAT:
      ! integer             VERSION
      ! integer             Z, A (parent)
      ! character(80)       OPERATOR
      ! integer             K, NR_POINTS, NXTERMS
      ! character(80) array XTERM_LABELS(:)        *** or EMPTY if no crossterms ***
      ! double        array TVALS(:)
      ! double        array DZ(:)
      ! dcomplex      array EQRPA(:)
      ! dcomplex      array SVALS(:,:)
      open(myfh, file=trim(fam_output_filename)//'.ctr', form='unformatted', &
         status='unknown', access='sequential', iostat=ierr)
      write(myfh) 2
      write(myfh) z_final-1, a
      write(myfh) f%label
      write(myfh) f%k, nr_points, nxterms
      if (nxterms > 0) then
         write(myfh) g(:)%label
      else
         write(myfh) 
      end if
      write(myfh) ctr_t(:)
      write(myfh) ctr_dzdt(:)
      write(myfh) ctr_z(:)
      write(myfh) ctr_strength(:,:)
      close(myfh)
         
   end subroutine

   
   !---------------------------------------------------------------------------
   ! An auxiliary formatting routine: padding with spaces on the left
   !---------------------------------------------------------------------------
   function txtr(str, n)
      implicit none
      integer,          intent(in) :: n
      character(len=*), intent(in) :: str
      character(len=n) :: txtr
      txtr = repeat(" ",n-len_trim(str))//trim(str)
   end function txtr
   
   !----------------------------------------------------------------------------
   ! Center a string in a field of a fixed width
   !----------------------------------------------------------------------------
   function txtc(text, width)
      implicit none
      character(len=*), intent(in) :: text
      integer,          intent(in) :: width
      character(len=width) :: txtc
      integer :: lpad
      
      ! Extra padding on left for odd number
      lpad = (width-len_trim(text))/2 + mod(width-len_trim(text),2)
      txtc = repeat(' ', lpad)//adjustl(text)
   end function txtc
   
   !----------------------------------------------------------------------------
   ! Get the wallclock time with or without OpenMP 
   !----------------------------------------------------------------------------
   function get_timer() result(t)
      implicit none
      integer  :: it, ir
      real(dp) :: t
      if (using_openmp .eqv. .true.) then
         !$ t = omp_get_wtime()
      else
         call system_clock(count=it, count_rate=ir)
         t = 1.0_dp*it/ir
      end if
   end function get_timer
      
end module pnfam_common

