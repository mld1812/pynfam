!----------------------------------------------------------------------------
! Temporary place to store high level routines as I refactor the code
!----------------------------------------------------------------------------
   ! The moment of inertia of the HFB solution in the cranking approximation
   ! Call this separately for protons and neutrons, add the results up
   function moment_of_inertia(E, U, V) result(moi)
      use external_field_type
      use pnfam_extfield
      use hfb_solution
      implicit none
      real(dp), intent(in) :: E(:)
      type(blockmatrix), intent(in) :: U, V  ! qp energies and amplitudes
      real(dp) :: moi
      
      type(blockmatrix) :: M1, M2
      type(external_field) :: Jplus, Jminus
      integer :: ipt, ibx1, ibx2, i1, i2, nd1, nd2, ix1, ix2
      
      ! Finite temperature / odd A
      if (ft_active) then
         call writelog('Warning: moment_of_inertia has not been updated for finite-temperature.')
      end if
      if (hfb_blo_active) then
         call writelog('Warning: moment_of_inertia has not been updated for odd-mass nuclei.')
      end if
      call writelog('Warning: moment_of_inertia has not been checked for using the PWI.')
      
      ! Because the blockmatrix type requires there to be only one nonzero block
      ! in each row and column, we need to separate the matrix into two:
      ! M1 = U^dagger J_+ V^* - V^dagger J_- U^*
      ! M2 = U^dagger J_- V^* - V^dagger J_+ U^*
      
      ! Construct raising and lowering operators of the total angular momentum
      Jplus%k = 1 ; Jplus%parity_even = .true. ; Jplus%rank = 1
      call init_fam_mapping(Jplus)
      ipt = 1
      do ibx1 = 1,nb
         ibx2 = Jplus%mat%ir2c(ibx1)
         if (ibx2 == 0) cycle
         nd1 = db(ibx1) ; nd2 = db(ibx2)
         if (nd1 == 0 .or. nd2 == 0) cycle
         
         do i2=1, nd2
            ! Calculate the true index of the state |2>
            ix2 = i2 + isstart(ibx2) - 1

            do i1=1, nd1
               ! Calculate the true index of the state <1|
               ix1 = i1 + isstart(ibx1) - 1
               
               !if (nr(ix1)==nr(ix2).and.nz(ix1)==nz(ix2)) then
                  if (nl(ix1)==nl(ix2) .and. ns(ix1)==ns(ix2)+2) then
                     Jplus%mat%elem(ipt) = dot_product(wf(:,ix1), wf(:,ix2))  ! ok
                  end if
                  if (nl(ix1)==nl(ix2)+1 .and. ns(ix1)==ns(ix2)) then
                     Jplus%mat%elem(ipt) = dot_product(-wf(:,ix1)/y(:), wfdz(:,ix2)) &
                        + dot_product(wf(:,ix1)*z(:), wfdr(:,ix2)) &
                        - nl(ix2)*dot_product(wf(:,ix1)*z(:)*y(:), wf(:,ix2))
                  end if
               !end if
               
               ipt = ipt + 1
            end do
         end do
      end do
      
      Jminus%k = -1 ; Jminus%parity_even = .true. ; Jminus%rank = 1
      call init_fam_mapping(Jminus)
      ipt = 1
      do ibx1 = 1,nb
         ibx2 = Jminus%mat%ir2c(ibx1)
         if (ibx2 == 0) cycle
         nd1 = db(ibx1) ; nd2 = db(ibx2)
         if (nd1 == 0 .or. nd2 == 0) cycle
         
         do i2=1, nd2
            ! Calculate the true index of the state |2>
            ix2 = i2 + isstart(ibx2) - 1

            do i1=1, nd1
               ! Calculate the true index of the state <1|
               ix1 = i1 + isstart(ibx1) - 1
               
               !if (nr(ix1)==nr(ix2).and.nz(ix1)==nz(ix2)) then
                  if (nl(ix1)==nl(ix2) .and. ns(ix1)==ns(ix2)-2) then
                     Jminus%mat%elem(ipt) = dot_product(wf(:,ix1)/y(:), wfdz(:,ix2)) &
                        - dot_product(z(:)*wf(:,ix1), wfdr(:,ix2)) &
                        - nl(ix2)*dot_product(wf(:,ix1)*z(:)*y(:), wf(:,ix2))
                  end if
                  if (nl(ix1)==nl(ix2)-1 .and. ns(ix1)==ns(ix2)) then
                     Jminus%mat%elem(ipt) = dot_product(wf(:,ix1), wf(:,ix2))
                  end if
               !end if

               ipt = ipt + 1
            end do
         end do
      end do
      
      ! Prep the matrices M1 and M2
      call allocate_blockmatrix(M1,size(Jplus%mat%elem))
      call allocate_blockmatrix(M2,size(Jplus%mat%elem))
      
      ! Transform them into the qp basis
      call triprod('T', U, 'N', Jplus%mat, 'N', V, 1d0, 0d0, M1)
      call triprod('T', V, 'N', Jminus%mat, 'N', U, -1d0, 1d0, M1)
      call triprod('T', U, 'N', Jminus%mat, 'N', V, 1d0, 0d0, M2)
      call triprod('T', V, 'N', Jplus%mat, 'N', U, -1d0, 1d0, M2)
      
      ! Sum over all qp pairs, block by block
      moi = 0
      ipt = 1
      do ibx1 = 1,nb
         ibx2 = M1%ir2c(ibx1)
         if (ibx2 == 0) cycle
         nd1 = db(ibx1) ; nd2 = db(ibx2)
         if (nd1 == 0 .or. nd2 == 0) cycle
         
         do i2=1, nd2
            ix2 = i2 + isstart(ibx2) - 1
            do i1=1, nd1
               ix1 = i1 + isstart(ibx1) - 1
               moi = moi + M1%elem(ipt)**2/(E(ix1) + E(ix2))
               ipt = ipt + 1
            end do
         end do
      end do
      ipt = 1
      do ibx1 = 1,nb
         ibx2 = M2%ir2c(ibx1)
         if (ibx2 == 0) cycle
         nd1 = db(ibx1) ; nd2 = db(ibx2)
         if (nd1 == 0 .or. nd2 == 0) cycle
         
         do i2=1, nd2
            ix2 = i2 + isstart(ibx2) - 1
            do i1=1, nd1
               ix1 = i1 + isstart(ibx1) - 1
               moi = moi + M2%elem(ipt)**2/(E(ix1) + E(ix2))
               ipt = ipt + 1
            end do
         end do
      end do
      
      moi = moi/4
      
   end function moment_of_inertia
   !============================================================================
   ! PNFAM HIGH LEVEL MODE - PNFAM COMMON VARIABLES
   !============================================================================
   ! contour mode
   integer  :: ratint_psi_gl_points, ratint_points
   real(dp) :: q_value, rot_correction, contour_anchor
   logical  :: include_endpoint, rotate_contour, use_gauleg_contour
   ! findmax mode
   integer  :: max_points
   real(dp) :: search_step, delta, energy_tolerance

   namelist /contour_parameters/ nr_points, q_value, rot_correction, include_endpoint, &
      rotate_contour, contour_anchor, use_gauleg_contour, ratint_psi_gl_points, ratint_points
   namelist /findmax_parameters/ energy_start, half_width, delta, energy_tolerance, &
      search_step, max_points

   ! Contour mode
   integer     :: nr_compute
   real(dp)    :: ctr_bound_left, ctr_bound_right, r0, r
   complex(dp) :: ci1, ci2, btot_simp, btot_trap, btot_gleg, rate_simp, rate_trap, rate_gleg, psi_simp, psi_trap, psi_gleg
   real(dp),    allocatable :: ctr_t(:), ratint_coeffs(:), ctr_gauleg_weights(:)
   complex(dp), allocatable :: ctr_z(:), ctr_dzdt(:), ctr_strength(:,:), ctr_psi_rat(:)
   ! Findmax mode
   logical :: interval_found
   real(dp) :: deriv, y1, y2
   integer :: dir, firstdir

   ! CONTOUR_PARAMETERS
   contour_anchor = 0.0_dp; include_endpoint = .false.; nr_points = 55
   use_gauleg_contour=.true.; ratint_psi_gl_points = 75; ratint_points = 20
   q_value = -1.0_dp; rotate_contour = .false.; rot_correction = -1.0_dp
   ! FINDMAX_PARAMETERS
   delta = 0.002_dp; energy_start = 20.0_dp; energy_tolerance = 0.005_dp
   half_width = 5.0_dp; max_points = 20; search_step = 2.0_dp
   ! STR_PARAMETERS
   energy_start = 0; energy_step = 0.2_dp; half_width = 0.5_dp; nr_points = 101

   !---------------------------------------------------------------------------
   ! ASCII output file
   !---------------------------------------------------------------------------
   subroutine write_output_file_header
      implicit none
      integer :: ixterms

      case("CONTOUR")
         write(fout,'(A)') "#"
         write(fout,'(6A)',advance='no') "#", txtr("Parameter", 29), txtr("Re(EQRPA)", 30), &
            txtr("Im(EQRPA)", 30), txtr("Re(Strength)", 34), txtr("Im(Strength)", 34)
         do ixterms=2, size(cstr)
           write(fout,'(2A)',advance='no') txtr("Re("//trim(g(ixterms-1)%label)//")", 34), &
              txtr("Im("//trim(g(ixterms-1)%label)//")", 34)
         end do
         write(fout,'(2A)', advance='no') txtr("Re(f(E)_rat)", 34), txtr("Im(f(E)_rat)", 34)
         write(fout,'(3A)') txtr("Re(dz/dt)", 30), txtr("Im(dz/dt)", 30), txtr("GL_Weights",30)
   end subroutine

   !============================================================================
   ! PNFAM SERIAL - CONTOUR AND FINDMAX CASES
   !============================================================================

   select case(fam_mode)
      case("CONTOUR")
         
         ! Initialize the contour and fit the phase-space polynomial
         call init_contour_mode(ierr)
         if (ierr /= 0) then
            write(error_unit,*) "Q-value is negative. Aborting beta decay calculation."
            stop
         end if
         
         ! Compute as many points as symmetrization requires, allowing for
         ! skipping the point on the right bound real axis if it has no phase space.
         ctr_strength(:,:) = 0
         do i=1, nr_compute
            if ((abs(aimag(ctr_z(i))) < epsilon(1.0)) &
                .and. abs((ctr_bound_right - real(ctr_z(i)))) < epsilon(1.0) &
                .and. (.not.rotate_contour) .and. (.not.include_endpoint)) then
               write(*,'(A,1F0.5,",",1F0.5,A)') 'N.b.: skipping point (', real(ctr_z(i)), &
                  aimag(ctr_z(i)), ') which should have no phase space'
               iter = 0; t_start = 0; t_finish = 0; str(:) = 0;  cstr(:) = 0
            else
               call cpu_time(t_start)
               call pnfam_solve(real(ctr_z(i),kind=dp), aimag(ctr_z(i)), Ep, En,            &
                  Vp, Up, Vn, Un, max_iter, convergence_epsilon, quench_residual_int, &
                  f%beta_minus, str(:), iter, cstr(:))
               call cpu_time(t_finish)
            end if
            wt_now = get_timer()
            
            ctr_strength(i,:) = cstr(:)
            call log_contour_point(ctr_z(i), cstr(1))
         end do
         
         call contour_finalize_computation
         call log_contour_footer
         call write_output_contour_footer
         call contour_binary_output
   
      case("FINDMAX")
         
         call init_findmax_mode
         
         ! The main loop: try to find a maximum
         omega = energy_start
         do i=1, max_points
      
            ! Solve the pnFAM equations for two points and approximate the derivative
            call cpu_time(t_start)
            call pnfam_solve(omega+delta, half_width, Ep, En, vp, up, vn, un, max_iter, &
               convergence_epsilon, quench_residual_int, f%beta_minus, str(:), iter)
            if (iter == -1) then
               write(*,'(A,1ES12.6)') 'FAM was interrupted at iteration limit. si = ', fam_si
               stop ! No point in trying to estimate the derivative, if no convergence
            end if
            y2 = str(1)
            
            call pnfam_solve(omega-delta, half_width, Ep, En, vp, up, vn, un, max_iter, &
               convergence_epsilon, quench_residual_int, f%beta_minus, str(:), iter)
            if (iter == -1) then
               write(*,'(A,1ES12.6)') 'FAM was interrupted at iteration limit. si = ', fam_si
               stop
            end if
            call cpu_time(t_finish)
            
            y1 = str(1)
            deriv = (y2 - y1)/(2*delta)
            
            call log_findmax_point
      
            if (deriv > 0) then
               dir = 1
            else
               dir = -1
            end if
            
            ! Check if we found an interval where the sign of the derivative flips
            if (i == 1) then
               firstdir = dir
               interval_found = .false.
            else
               if (dir /= firstdir) interval_found = .true.
            end if
            
            ! When we know the interval, we keep bisecting it...
            if (interval_found) search_step = search_step/2
            omega = omega + dir*search_step
      
            ! ...until we know the location of the maximum to the requested tolerance
            if (search_step < energy_tolerance) exit
      
         end do
   
         call write_output_findmax_result
         call log_findmax_footer
         call closelog
         
   end select

   !============================================================================
   ! THIS ROUTINE IS ONLY USED IN PARALLEL MODE... CAN BE DEPRECATED, JUST HAVE
   ! RANK 0 READ NAMELIST...
   !============================================================================
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
      character(len=80)  :: fam_mode
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
   !============================================================================
   ! THESE ROUTINES ARE USED FOR HIGH-LEVEL MODES
   !============================================================================
   subroutine log_header
      implicit none
      real(dp) :: estop=0

      write(st,'(1x,A)') "Mode-specific pnFAM details"; call writelog(st)
      write(st,'(1x,A)') repeat("-",w); call writelog(st)

      select case(fam_mode)

         case("STR")
            estop = energy_start+energy_step*(nr_points-1)
            write(st,'(1x,A,I0)') &
               "Number fam points:               ", nr_points;          call writelog(st)
            write(st,'(1x,A,F8.4," MeV")') &
               "Half-width (gamma):              ", half_width;         call writelog(st)
            write(st,'(1x,A,F0.4," MeV")') &
               "Energy step:        ", energy_step;                     call writelog(st)
            write(st,'(1x,A,F0.4," ... ",F0.4," MeV")') &
               "Energy range:                    ", energy_start,estop; call writelog(st)
            write(st,'(1x,A,I0)') &
               "Block matrix dimension:          ", dqp;                call writelog(st)
            write(st,'(1x,A,I0)') &
               "Number matrix blocks:            ", nb;                 call writelog(st)
            write(st,'(1x,A,I0)') &
               "Non-trivial HFB matrix elements: ", dmat;               call writelog(st)
            write(st,'(1x,A,I0)') &
               "Non-trivial FAM matrix elements: ", size(f%mat%elem);   call writelog(st)


         case("CONTOUR")
            write(st,'(1x,A,I0)') &
               "Number fam points:             ", nr_points;                             call writelog(st)
            write(st,'(1x,2A)') &
               "Use Gauss-Legendre contour:    ", merge("Yes","No ",use_gauleg_contour); call writelog(st)
            write(st,'(1x,A,1I0)') &
               "Rational-fct interp. Points:   ", ratint_points;                         call writelog(st)
            write(st,'(1x,A,1I0)') &
               "PSI Quadrature Points          ", ratint_psi_gl_points;                  call writelog(st)
            write(st,'(1x,2A)') &
               "Include high-omega endpoint:   ", merge("Yes", "No ", include_endpoint); call writelog(st)
            write(st,'(1x,2A)') &
               "Rotate contour to avoid poles: ", merge("Yes", "No ", rotate_contour);   call writelog(st)
            write(st,'(1x,A,1F7.3," MeV")') &
               "Contour left anchor point:     ", contour_anchor;                        call writelog(st)

         case("FINDMAX")
            write(st,'(1x,A,I0)') &
               "Maximum number fam points:            ", max_points ;       call writelog(st)
            write(st,'(1x,A,F8.5,A," MeV")') &
               "Half-width (gamma):                   ", half_width ;       call writelog(st)
            write(st,'(1x,A,F8.5,A," MeV")') &
               "Initial guess for the peak:           ", energy_start ;     call writelog(st)
            write(st,'(1x,A,F8.5,A," MeV")') &
               "Initial search step size:             ", search_step ;      call writelog(st)
            write(st,'(1x,A,ES7.1,A," MeV")') &
               "Delta for estimating derivatives:     ", delta ;            call writelog(st)
            write(st,'(1x,A,ES7.1,A," MeV")') &
               "Requested accuracy for peak location: ", energy_tolerance ; call writelog(st)

      end select

   end subroutine log_header

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

      call writelog(repeat('=',140))
      write(st,'(26x, a20,a20,a20, a15,a12,a15,a12)') "Simpson's 3/8","Trapezoidal","Gauss-Legendre", &
          "(GL-S)", "(GL-S)/S", "(GL-T)", "(GL-T)/T"

      call writelog(st)
      call writelog(repeat('=',140))

      call contour_results
      
      ! Only available with beta-
      if (f%beta_minus) then
         write(st,'(2a)') 'Im[Phase space] ........ :', contour_compare(aimag(psi_simp),aimag(psi_trap), &
             aimag(psi_gleg))
         call writelog(st)
         write(st,'(2a)') 'Re[Contour integral] ... :', &
            contour_compare(-2/ln2_kappa*real(rate_simp,kind=dp),-2/ln2_kappa*real(rate_trap,kind=dp), &
            -2/ln2_kappa*real(rate_gleg,kind=dp))
         call writelog(st)
      else
         write(st,'(a)') 'Im[Phase space] ........ :         N/A; Beta+'; call writelog(st)
         write(st,'(a)') 'Re[Contour integral] ... :         N/A; Beta+'; call writelog(st)
      end if
      
      call writelog(repeat('=',140))
      write(st,'(2a)') 'Transition strength .... :', &
         contour_compare(aimag(btot_simp),aimag(btot_trap), aimag(btot_gleg))
         call writelog(st)
      
      ! Only available with beta-
      if (f%beta_minus) then
         write(st,'(2a)') 'Total rate .........(/s) :', &
            contour_compare(aimag(rate_simp),aimag(rate_trap), aimag(rate_gleg))
         call writelog(st)
         write(st,'(2a)') 'Half-life ...........(s) :', &
            contour_compare(ln2/aimag(rate_simp),ln2/aimag(rate_trap), ln2/aimag(rate_gleg))
         call writelog(st)         
      else
         write(st,'(2a)') 'Total rate .........(/s) :         N/A; Beta+'; call writelog(st)
         write(st,'(2a)') 'Half-life ...........(s) :         N/A; Beta+'; call writelog(st)
      end if
      
      call writelog(repeat('=',140))
      
      call contour_tests
      
      wt_now = get_timer()
      call writelog("")
      call writelog("fin.")
      write(st,'(a,1f0.3,a)') 'wall-time was ', wt_now-wt_start, ' seconds.'
      call writelog(st)
      
   end subroutine
   
   !---------------------------------------------------------------------------
   ! Contour computation initializations:  Q-value, contour points and such
   !---------------------------------------------------------------------------
   subroutine init_contour_mode(qerr)
      implicit none
      character(len=160) :: st
      logical :: lqnml, lrotnml
      integer, intent(out) :: qerr

      qerr = 0
      
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
      ! Right bound --- q_value is set in the namelist
      ctr_bound_right = 0
      if (lqnml) then
         ctr_bound_right = egs_hfb + q_value         
      else
         if (q_hfb <= 0) then
            write(st,'(" FATAL ERROR: Q-value computed from HFB solution is negative. Value = ",1F12.6," MeV.")') &
            q_hfb
            call writelog(st)
            qerr = 1
            return
         end if
         ctr_bound_right = egs_hfb + q_hfb
      end if
      ! Left bound --- contour_anchor is set in the namelist
      ctr_bound_left = contour_anchor

      ! Adjust the bounds if the q.p. energies have been shifted
      if (abs(energy_shift_neut) > 1d-10 .or. abs(energy_shift_prot) > 1d-10) then
         write(st,'(1A)') '*** NB: QP energies are shifted, so the contour will shift to compensate ***'
         call writelog("")
         call writelog(st)
      end if
      ctr_bound_left  = ctr_bound_left  + energy_shift_neut + energy_shift_prot
      ctr_bound_right = ctr_bound_right + energy_shift_neut + energy_shift_prot

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
      write(st,'("Contour L bound:   ",1F9.6," MeV  (set in namelist CONTOUR_PARAMETERS)")') &
         ctr_bound_left
      call writelog(st)
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
      allocate(ratint_coeffs(ratint_points))
      ratint_coeffs(:) = 0.0
      
      ! Beta+ is not considered yet --- only set up the polynomial fit if in beta- mode
      if (f%beta_minus) then
         ! Calc the psi corresponding to rantin_x, then get the interp coeffs
         ratint_coeffs = get_ratint_coeffs(nz=z_final, na=a, npts=ratint_points, &
             w0max=1+(ctr_bound_right-ctr_bound_left)/mec2, gn=2, nleg=ratint_psi_gl_points)
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
   subroutine write_output_contour_point(ip)
      implicit none
      integer, intent(in) :: ip
      integer :: i
      
      ! Write the sampled points in the output file for diagnostic purposes
      write(fno,'(3f30.19)',advance='no') ctr_t(ip), &
         ctr_z(ip)
      
      ! Write the strengths
      do i=1, size(cstr)
         write(fno,'(2ES34.19)',advance='no') ctr_strength(ip,i)
      end do
      
      ! Write the phase-space integral value here
      write(fno,'(2ES34.19)', advance='no') ctr_psi_rat(ip)

      ! Write the derivative of the contour (so we have everything we need
      ! in the outputfile to reconstruct the integrands by hand)
      write(fno,'(2f30.19)',advance='no') ctr_dzdt(ip)

      ! Write the gauss-legendre weights if active
      if (use_gauleg_contour) then
         write(fno, '(1ES30.19)') ctr_gauleg_weights(ip)
      else
         write(fno,*)
      end if
      
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
      write(fno,'("# ",a)') repeat('-',128)
      write(fno,'("# ",a,9x,a,16x,a,16x,a)') 'Quantity (RFI Phase Space)', "Simpson's Rule", "Trapezoidal Rule", "Gauss-Legendre"
      write(fno,'("# ",a)') repeat('-',128)
      write(fno,'("# Transition strength ... : ",3ES30.15)') aimag(btot_simp), aimag(btot_trap), aimag(btot_gleg)
      
      if (f%beta_minus) then
         write(fno,'("# Decay rate ..... (s^-1) : ",3ES30.15)') aimag(rate_simp),     aimag(rate_trap), aimag(rate_gleg)
         write(fno,'("# Half-life ..........(s) : ",3ES30.15)') ln2/aimag(rate_simp), ln2/aimag(rate_trap), ln2/aimag(rate_gleg)
      else
         write(fno,'("# Decay rate ..... (s^-1) :          N/A; Beta+")')
         write(fno,'("# Half-life ..........(s) :          N/A; Beta+")')
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
      use complex_quadrature, only: gauleg
      implicit none
      integer  :: icp, nv
      real(dp) :: t0
      
      ! Gauleg mode cannot rotate contour (risk of ruining reflection symmetry)
      ! so ignore namelist input and force false.
      if (use_gauleg_contour) then
         if (rotate_contour) call writelog("WARNING: Rotate contour forced to false for gauleg contour.")
         rotate_contour = .false.
      end if

      allocate(ctr_z(nr_points), ctr_t(nr_points), ctr_dzdt(nr_points))
      allocate(ctr_psi_rat(nr_points), ctr_strength(nr_points,size(cstr)))
      allocate(ctr_gauleg_weights(nr_points))

      ctr_z = 0; ctr_t = 0; ctr_dzdt = 0; ctr_strength = 0
      ctr_psi_rat = 0
      
      r0 = (ctr_bound_right+ctr_bound_left)/2
      r  = ctr_bound_right-r0

      
      ! Construct the contour grid
      if (.not. use_gauleg_contour) then
         ! Spin the contour away from the real line if desired
         nv = nr_points-1
         t0 = 0
         if (rotate_contour) t0 = pi/nv

         ! Evenly spaced theta grid
         do icp=1, nr_points
           ctr_t(icp) = t0 + (icp-1)*2*pi/nv
         end do

         ! Circular contour and its derivative
         ctr_z(:)              = r0 + r*exp(iu*ctr_t(:))
         ctr_z(nr_points)      = ctr_z(1)
         ctr_dzdt(:) = iu*(ctr_z(:) - r0)

         ctr_gauleg_weights(:) = 0 ! This is unused in this mode
      else
         ! Force gauleg to start at origin and go ccw so we don't get points very close 
         ! to real axis at right bound which take forever to converge.
         t0 = pi

         ! legendre theta grid. NB: weights only depend on difference between end points
         ! so the shift only affects ctr_t, which is solely used to get desired ctr_z.
         call gauleg(nr_points, t0, t0+2*pi, ctr_gauleg_weights, ctr_t)

         ! Circular contour and its derivative
         ctr_z(:) = r0 + r*exp(iu*ctr_t(:))
         ctr_dzdt(:) = iu*(ctr_z(:) - r0)
      end if

      ! The number of points to compute depends on the orientation, as do the
      ! symmetry properties of the contour. Integer division gets this right.
      if (rotate_contour) then
         nr_compute = nr_points/2
      else
         nr_compute = (nr_points+1)/2
      end if
   end subroutine contour_setup
   
   !----------------------------------------------------------------------------
   !> Use symmetry to fill in the strength function via 
   !> \f[ S(\omega) = S*(\omega*) \f]. 
   !> Also calculate the phase space along the contour to calculate the 
   !> contribution to the rate from the operator that was run.
   !----------------------------------------------------------------------------
   subroutine contour_finalize_computation
      implicit none
      integer     :: ip, iflip
      complex(dp) :: ener, psint_rat
      
      ! NB: Take special care of the endpoints! If ctr_t(1)=ctr_t(nr_points)+2pi then 
      ! we are okay to force equal strengths since they are the same contour point. 
      ! With legendre ctr_t starting at 0.0 going ccw, this is not the case. The 1st 
      ! point is slightly below the real axis and the last is slightly above, so they 
      ! should be conjugates.
      if (ctr_z(1) == ctr_z(nr_points)) then
          ctr_strength(nr_points,:) = ctr_strength(1,:)
      else
          ctr_strength(nr_points,:) = conjg(ctr_strength(1,:))
      end if

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
         psint_rat  = ratint_eval_psi(npts=ratint_points, &
             w0max=1+(ctr_bound_right-ctr_bound_left)/mec2, &
             c=ratint_coeffs(:), z=ener)
         ctr_psi_rat(ip) = psint_rat
         
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
      if (use_gauleg_contour) then
         btot_gleg = cquad_gauss(f=ctr_strength(:,1), dz=ctr_dzdt(:), &
             wts=ctr_gauleg_weights(:))
      else
         btot_gleg = 0
      end if

      ! Phase space only
      psi_trap = cquad_trapezoidal(f=ctr_psi_rat, dz=ctr_dzdt, t=ctr_t)
      psi_simp = cquad_simpson(f=ctr_psi_rat, dz=ctr_dzdt, t=ctr_t)
      if (use_gauleg_contour) then
         psi_gleg = cquad_gauss(f=ctr_psi_rat, dz=ctr_dzdt(:), wts=ctr_gauleg_weights(:))
      else
         btot_gleg = 0
      end if
      ! Total rate
      rate_trap = cquad_trapezoidal(f=ctr_strength(:,1)*ctr_psi_rat, dz=ctr_dzdt, t=ctr_t)
      rate_simp = cquad_simpson(f=ctr_strength(:,1)*ctr_psi_rat, dz=ctr_dzdt, t=ctr_t)
      if (use_gauleg_contour) then
         rate_gleg = cquad_gauss(f=ctr_strength(:,1)*ctr_psi_rat, dz=ctr_dzdt(:), &
             wts=ctr_gauleg_weights(:))
      else
         btot_gleg = 0
      end if

      ! Constants
      ! Now, B = Im(Btot) and Rate = Im(Rate)
      btot_trap = -half*btot_trap
      btot_simp = -half*btot_simp
      btot_gleg = -half*btot_gleg
      psi_trap  = -half*psi_trap
      psi_simp  = -half*psi_simp
      psi_gleg  = -half*psi_gleg
      rate_trap = -half*ln2_kappa*rate_trap
      rate_simp = -half*ln2_kappa*rate_simp
      rate_gleg = -half*ln2_kappa*rate_gleg
   end subroutine contour_results
   
   !----------------------------------------------------------------------------
   ! Compare two real quantities 
   !----------------------------------------------------------------------------
   function contour_compare(a, b, c, fmt_exp)
      implicit none
      real(dp), intent(in) :: a, b, c
      logical, optional, intent(in) :: fmt_exp
      character(len=100) :: fmt
      character(len=160) :: contour_compare
      
      fmt = '(1es20.10,1es20.10,1es20.10, 1es15.4,1f11.2"%",1es15.4,1f11.2"%")'
      write(contour_compare,trim(fmt)) a, b, c, c-a, 100*(c-a)/a, c-b, 100*(c-b)/b

   end function contour_compare
   
   !----------------------------------------------------------------------------
   ! Check various numerical properties of the integrals 
   !----------------------------------------------------------------------------
   subroutine contour_tests
      implicit none
      real(dp), parameter :: tol = 1.0d-4
      real(dp) :: simp, trap, gleg
      logical  :: lsimp, ltrap, lgleg
      
      ! Strength function --- 
      ! Residues are all real, so the real part of B_tot should vanish!
      trap = real(btot_trap,kind=dp); ltrap = (abs(trap) > tol)
      simp = real(btot_simp,kind=dp); lsimp = (abs(simp) > tol)
      gleg = real(btot_gleg,kind=dp); lgleg = (abs(gleg) > tol)
      if (lsimp .or. ltrap .or. lgleg) then
         call writelog("")
         call writelog(txtc(repeat('*',95),101))
         if (lgleg) then
            call writelog(txtc("NOTE: Re[B] is non-zero with Gauss-Legendre.",101))
         end if
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
      gleg = real(psi_gleg,kind=dp); lgleg = (abs(gleg) > tol)

      if (lsimp .or. ltrap .or. lgleg) then
         call writelog(txtc(repeat('*',95),101))
         if (lgleg) then
            call writelog(txtc("NOTE: Re[Phase space] is non-zero with Gauss-Legendre.",101))
         end if
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
      gleg = aimag(psi_gleg); lgleg = (abs(gleg) > tol)
      if (lsimp .or. ltrap .or. lgleg) then
         call writelog(txtc(repeat('*',95),101))
         if (lgleg) then
            call writelog(txtc("NOTE: Im[Phase space] is non-zero with Gauss-Legendre.",101))
         end if
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
      gleg = -2/ln2_kappa*real(rate_gleg,kind=dp); lgleg = (abs(gleg) > tol)
      if (lsimp .or. ltrap .or. lgleg) then
         call writelog(txtc(repeat('*',95),101))
         if (lgleg) then
            call writelog(txtc("NOTE: Re[Contour integral] is non-zero with Gauss-Legendre.",101))
         end if
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
      gleg = aimag(rate_gleg); lgleg = (gleg < 0)
      if (lsimp .or. ltrap .or. lgleg) then
         call writelog(txtc(repeat('*',95),101))
         if (lgleg) then
            call writelog(txtc("NOTE: Decay rate is negative with Gauss-Legendre.",101))
         end if
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
      write(myfh) 3
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
      write(myfh) merge(1,0, use_gauleg_contour)
      write(myfh) ctr_gauleg_weights(:)
      close(myfh)
         
   end subroutine

   !============================================================================
   ! THESE ROUTINES ARE USED HIGH LEVEL MODES (FINDMAX)
   !============================================================================
   subroutine log_findmax_point
      implicit none
      character(len=160) :: st
      
      write(st,'(3(1F11.6,4X),I0)') omega, deriv, (y1+y2)/2, iter ; call writelog(st)
      
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

   subroutine init_findmax_mode

      implicit none
      character(len=160) :: st
      
      write(st,'(A,I0)') "Number of non-trivial density matrix elements: ", size(f%mat%elem)
      call writelog(st)
      call writelog("")
      call writelog("    Energy      derivative       strength  ~iter")
      call writelog(repeat('-',50))
      
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

   !============================================================================
   ! THESE ROUTINES ARE USED FOR BETADECAY
   !============================================================================
   !----------------------------------------------------------------------------
   ! Compute beta decay quantities using an already-loaded HFB solution.
   ! Input:  Ep(:), En(:)
   ! Output: Q, E_gs, K_gs, P_gs (Q-value & g.s. energy, K-projection, parity)
   !----------------------------------------------------------------------------
   subroutine get_hfbtho_betadecay_properties(Ep, En, Vp, Vn, Q, E_gs, K_gs, P_gs, lpr)
      use logger
      use constants, only : dmassnH, IT_NEUTRON, IT_PROTON
      implicit none

      type(blockmatrix), intent(in)  :: Vp, Vn
      real(dp),          intent(in)  :: Ep(:), En(:)
      logical, optional, intent(in)  :: lpr
      real(dp),          intent(out) :: Q, E_gs
      integer,           intent(out) :: K_gs, P_gs

      real(dp), parameter :: occ_tolerance = 5.0d-09

      integer  :: iqp, iprot, ineut
      real(dp) :: eqrpa_max, Op(size(Ep)), On(size(En))
      logical  :: mask_p(size(Ep)), mask_n(size(En))

      ! PWI
      integer :: ik
      logical :: kactive

      character(len=160) :: st

      ! Check for an HFB solution
      if (abs(hfb_lambda(1)) < epsilon(1.0_dp) .and. abs(hfb_lambda(2)) < epsilon(1.0_dp)) then
         call writelog('ERROR in sub.get_HFBTHO_betadecay_properties: no HFB solution found.&
            & Both lambda values are zero.')
         stop
      end if

      ! Q and E_gs via J. Engel et al., Phys. Rev. C 60, 014302 (1999)
      ! Note: When pairing vanishes, lambda is arbitrary. It is contained in
      ! all QP energies and used to define EQRPA_max. In principle, we should
      ! use a different value (e.g. the last bound/hole s.p. energy), and shift
      ! particle/hole energies in the appropriate direction. However, this adds
      ! unneccessary complication.
      !   1. Which new lambda to choose is ambiguous, particularly for odd nuclei
      !   2. Changing lambda would shift Egs and EQRPA_max by the same amount, so the
      !      relevant energy interval (for the contour) and Q is unchanged.
      !   3. The QRPA phonon spectrum (i.e. QRPA eigenvalues) shift by the same amount
      !      as E_gs and EQRPA_max, so the relative locations of peaks do not change.
      !      This also offsets the strength from being centered at zero, which is less ideal.
      !   4. The beta-decay integrals (at least for beta-minus... haven't
      !      checked beta-plus) are unchanged with a change of lambdas.
      ! As a result, there's no need to make any correction.
      eqrpa_max = dmassnH + hfb_lambda(1) - hfb_lambda(2)

      !-------------------------------------------------------------------------
      ! Compute ground-state energy
      !-------------------------------------------------------------------------
      mask_p(:) = .true.
      mask_n(:) = .true.

      ! Compute q.p. occupations (norm of lower component, at least) to mask the q.p. energies
      call get_hfbtho_occupations(Vp, IT_PROTON, Op)
      call get_hfbtho_occupations(Vn, IT_NEUTRON, On)

      ! Exclude full proton levels, empty neutron levels, and above PWI
      do iqp=1, size(Op)
         ! Neutron PWI
         kactive = .true.
         do ik=1, sum(pwi_dim_n)
            if (pwi_qp_n(ik) == iqp) then
               kactive = .true.
               exit
            end if
         end do

         if ((.not.kactive) .or. (abs(On(iqp)) < occ_tolerance)) then
            mask_n(iqp) = .false.
         end if

         ! Proton PWI
         kactive = .true.
         do ik=1, sum(pwi_dim_p)
            if (pwi_qp_p(ik) == iqp) then
               kactive = .true.
               exit
            end if
         end do

         if ((.not.kactive) .or. (abs(1.0_dp-Op(iqp)) < occ_tolerance)) then
            mask_p(iqp) = .false.
         end if
      end do

      ! Indices of lowest-energy proton and neutron levels which satisfy
      ! occupation requirements
      iprot = minloc(Ep, 1, mask_p)
      ineut = minloc(En, 1, mask_n)

      ! Ground-state energy is the lowest HFB 2-QP energy
      E_gs = Ep(iprot) + En(ineut)

      ! Q-value is now simply EQRPA_max - E_gs
      Q = eqrpa_max - E_gs

      ! K and parity
      K_gs = nl(iprot) + nl(ineut) + (ns(iprot) + ns(ineut))/2
      P_gs = (-1)**(npar(iprot)+npar(ineut))

      ! If lpr is True, print HFB properties
      if (present(lpr)) then
         if (lpr) then
            write(st,'(a)') 'HFB PARAMETERS:';                                    call writelog(st)
            write(st,'(a)') repeat('-',30);                                       call writelog(st)
            write(st,'(a,3x,1i0)')    'Z (parent) ............ :', hfb_npr(2);    call writelog(st)
            write(st,'(a,3x,1i0)')    'N (parent) ............ :', hfb_npr(1);    call writelog(st)
            write(st,'(a,3x,1i0)')    'A ..................... :', hfb_npr(3);    call writelog(st)
            write(st,'(a,1x,1f10.6)') 'Lambda prot. (ALA) (MeV):', hfb_lambda(2); call writelog(st)
            write(st,'(a,1x,1f10.6)') 'Lambda neut. (ALA) (MeV):', hfb_lambda(1); call writelog(st)
            write(st,'(a,1x,1f10.6)') 'Delta prot. ...... (MeV):', hfb_delta(2);  call writelog(st)
            write(st,'(a,1x,1f10.6)') 'Delta neut. ...... (MeV):', hfb_delta(1);  call writelog(st)
            write(st,'(a,1x,1f10.6)') 'MIN(Ep) .......... (MeV):', Ep(iprot);     call writelog(st)
            write(st,'(a,1x,1f10.6)') 'Occupation (vp**2) .... :', Op(iprot);     call writelog(st)
            write(st,'(a,1x,1f10.6)') 'MIN(En) .......... (MeV):', En(ineut);     call writelog(st)
            write(st,'(a,1x,1f10.6)') 'Occupation (vn**2) .... :', On(ineut);     call writelog(st)
            write(st,'(a,1x,1f10.6,2x,"[",1i0,a,"]")') 'Ground state ..... (MeV):', E_gs, &
               K_gs, merge('+', '-', P_gs > 1); call writelog(st)
            write(st,'(a,1x,1f10.6)') 'Q-value .......... (MeV):', Q;         call writelog(st)
            write(st,'(a,1x,1f10.6)') 'Max EQRPA ........ (MeV):', eqrpa_max; call writelog(st)
         end if
      end if

   end subroutine get_hfbtho_betadecay_properties



   !----------------------------------------------------------------------------
   ! Return the lowest-energy q.p. from the list *energies* while only
   ! considering q.p.'s below the PWI.  Optionally request only partially-
   ! occupied or partially-unoccupied levels.  Variable *it* determines whether
   ! to mask using the neutron (IT_NEUTRON) or proton (IT_PROTON) PWI arrays.
   !----------------------------------------------------------------------------
   function lowest_energy_qp(energies, v, it, mask_k_positive, mask_k_negative, &
      mask_2k, mask_parity, mask_occ, mask_unocc) result(iqp)

      use constants, only : IT_PROTON, IT_NEUTRON
      use logger

      implicit none

      real(dp), parameter :: occ_cut = 5d-9

      real(dp),          intent(in) :: energies(:)
      type(blockmatrix), intent(in) :: v
      integer,           intent(in) :: it
      integer, optional, intent(in) :: mask_2k, mask_parity
      logical, optional, intent(in) :: mask_k_positive, mask_k_negative, mask_occ, mask_unocc

      integer  :: iqp, i, n_degen
      integer, pointer :: pwi_start(:), pwi_dim(:), pwi_qp(:)
      real(dp) :: occ(size(v%elem))
      logical  :: active, mask(size(energies))
      character(len=200) :: s

      if (it /= IT_NEUTRON .and. it /= IT_PROTON) then
         write(s,'("ERROR in f.lowest_energy_qp: unrecognized isospin --- ",1i0)') it
         call writelog(s)
         stop
      else
         if (it == IT_PROTON) then
            pwi_start => pwi_start_p
            pwi_dim   => pwi_dim_p
            pwi_qp    => pwi_qp_p
         else if (it == IT_NEUTRON) then
            pwi_start => pwi_start_n
            pwi_dim   => pwi_dim_n
            pwi_qp    => pwi_qp_n
         end if
      end if

      call get_hfbtho_occupations(v,it,occ)

      mask(:) = .true.
      do i=1, size(energies)
         ! PWI
         active = ((it == IT_NEUTRON .and. pwi_active_n(i)) .or. &
                   (it == IT_PROTON  .and. pwi_active_p(i)))
         if (.not.active) then
            mask(i) = .false.
         end if
         ! Occupation
         if (present(mask_occ)) then
            if (mask_occ .and. (abs(occ(i)) < occ_cut)) mask(i) = .false.
         end if
         if (present(mask_unocc)) then
            if (mask_unocc .and. (abs(1.0_dp-occ(i)) < occ_cut)) mask(i) = .false.
         end if
         ! K > 0 option
         if (present(mask_k_positive)) then
            if (mask_k_positive .and. (2*nl(i)+ns(i) < 0)) mask(i) = .false.
         end if
         ! K < 0 option
         if (present(mask_k_negative)) then
            if (mask_k_negative .and. (2*nl(i)+ns(i) > 0)) mask(i) = .false.
         end if
         ! K option
         if (present(mask_2k)) then
            if ((2*nl(i)+ns(i)) /= mask_2k) mask(i) = .false.
         end if
         ! Parity option
         if (present(mask_parity)) then
            if (npar(i) /= mask_parity) mask(i) = .false.
         end if
      end do

      iqp = minloc(energies, 1, mask)

      ! Degeneracy should be a warning
      n_degen = 0
      do i=1, size(energies)
         if (.not.mask(i)) cycle
         if (abs(energies(i)-energies(iqp)) < 1d-6) n_degen = n_degen+1
      end do
      if (n_degen > 1) then
         write(s,'("WARNING in f.lowest_energy_qp: found ",1i0," degenerate not-fully-occupied qp levels")') n_degen
         call writelog(s)
         call writelog('')
      end if

   end function lowest_energy_qp

   !============================================================================
   ! THESE ROUTINES ARE NOT REALLY NEEDED, SINCE WE CAN GET THEM FROM HFBTHO LIBRARY
   !
   ! HOWEVER, THEY COULD BE USEFUL IN RECONSTRUCTING QUANTITIES FROM U,V
   ! FOR OTHER SOLVERS - NEED TO ELIMINATE DEPENDENCE ON PWI THOUGH.
   !============================================================================
   !----------------------------------------------------------------------------
   ! Compute q.p. occupations using blockmatrix V.
   ! This routine requires get_HFBTHO_solution to have been called previously
   ! to set the blockmatrix dimension variables db and nb.
   !----------------------------------------------------------------------------
   subroutine get_hfbtho_occupations(v,it, occ)
      use constants, only : IT_PROTON, IT_NEUTRON

      implicit none

      type(blockmatrix), intent(in) :: v
      integer, intent(in) :: it
      real(dp), intent(out) :: occ(sum(db(:)))

      integer :: ib, iqp, istart, iend, nme, nqp, iqp_all
      integer :: blo_qp, blo_qpr

      occ(:) = 0
      nme = 0
      nqp = 0

      ! Blocked QP indices
      blo_qp = 0; blo_qpr = 0
      blo_qp = hfb_blo_qp(it)
      if (blo_qp /= 0) blo_qpr = qp_t_reverse(blo_qp)

      do ib=1, nb
         ! Recall V are V_{p\pi}, i.e. q.p.s are column vectors.
         ! Norm (occupation/2) should be <column|column> = sum(col(:)**2).
         do iqp=1, db(ib)
            istart = nme + (iqp-1)*db(ib) + 1
            iend   = nme + iqp*db(ib)
            iqp_all= nqp + iqp
            ! Because input is the blockmatrix V, there is no factor of 2 for
            ! the time-reversed partner. It is treated explicitly.
            occ(iqp_all) = sum(v%elem(istart:iend)**2)
            ! EFA occupancies are V^2 - 1/2(V^2 - U^2) = 1/2 (since U^2 + V^2 = 1)
            if (hfb_blo_active .and. blo_qp /= 0) then
               if (iqp_all == blo_qp .or. iqp_all == blo_qpr) then
                  occ(iqp_all) = 0.5_dp
               end if
            end if
         end do

         nme = nme + db(ib)**2
         nqp = nqp + db(ib)
      end do

   end subroutine get_hfbtho_occupations

   !---------------------------------------------------------------------------
   ! Compute the normal density of the mean field in COORDINATE space.
   ! Extended to finite temperature and the PWI cut from HFBTHO.
   ! it = 1 neutron density
   !    = 2 proton density
   !    = 3 total density
   !
   ! If blocking is active, the density will be the equal-filling s.p. density.
   !---------------------------------------------------------------------------
   subroutine hfb_density(it, vp, vn, up, un, rho)
      use constants, only : IT_NEUTRON, IT_PROTON, IT_ISOSCALAR
      use blockmatrix_type
      implicit none
      integer,  intent(in)  :: it
      real(dp), intent(in)  :: vp(:), vn(:), up(:), un(:)
      real(dp), intent(out) :: rho(nghl)

      real(dp) :: rho_ab, f
      integer  :: ib, nd, ir, ic, i1, i2, k, iuv, iqp, blo_n, blo_p, blo_nr, blo_pr

      ! Check of isospin
      if (it /= IT_NEUTRON .and. it /= IT_PROTON .and. it /= IT_ISOSCALAR) then
         write(*,'(A)') 'ERROR in hfb_density(): Choose from neutron,&
            & proton, or total density.'
         stop
      end if
      
      ! No work has been done to handle T>0 and odd nuclei simultaneously
      if (hfb_blo_active .and. ft_active) then
         write(*,'(a)') 'ERROR: this code cannot handle T>0 and odd-Z/odd-N simultaneously.'
         stop
      end if
      
      ! Blocking q.p. indices
      blo_nr = 0;  blo_pr = 0
      
      blo_n = hfb_blo_qp(IT_NEUTRON)
      blo_p = hfb_blo_qp(IT_PROTON)
      
      if (blo_n /= 0) blo_nr = qp_t_reverse(blo_n)
      if (blo_p /= 0) blo_pr = qp_t_reverse(blo_p)
      

      rho(:) = 0

      do ib=1, nb
         nd = db(ib)

         do ic=1, nd
            do ir=1, nd
               ! sp w.f. indices
               i1 = isstart(ib)-1 + ir
               i2 = isstart(ib)-1 + ic

               !----------------------------------------------------------------
               ! Neutrons
               !----------------------------------------------------------------
               if (it == IT_NEUTRON .or. it == IT_ISOSCALAR) then
                  rho_ab = 0
                  do k=pwi_start_n(ib), pwi_start_n(ib)+pwi_dim_n(ib)-1
                     iuv = pwi_uv_n(k)
                     iqp = pwi_qp_n(k)

                     ! Finite temperature and odd nucleus
                     f = 0
                     if (ft_active) then
                        f = ft_fn(k)   ! Access f variables through pwi indices, not qp indices!
                     else if (hfb_blo_active .and. blo_n /= 0) then
                        if (iqp == blo_n .or. iqp == blo_nr) then
                           f = 0.5_dp  ! Equal filling
                        end if
                     end if

                     rho_ab = rho_ab + vn(iuv+ir)*vn(iuv+ic)*(1.0_dp-f) + un(iuv+ir)*un(iuv+ic)*f
                  end do

                  if (ns(i1) == ns(i2) .and. nl(i1) == nl(i2) .and. npar(i1) == npar(i2)) then
                     rho(:) = rho(:) + rho_ab*wf(:,i1)*wf(:,i2)
                  end if
               end if
               !----------------------------------------------------------------
               ! Protons
               !----------------------------------------------------------------
               if (it == IT_PROTON .or. it == IT_ISOSCALAR) then
                  rho_ab = 0
                  do k=pwi_start_p(ib), pwi_start_p(ib)+pwi_dim_p(ib)-1
                     iuv = pwi_uv_p(k)
                     iqp = pwi_qp_p(k)

                     f = 0
                     if (ft_active) then
                        f = ft_fp(k)   ! Access f variables through pwi indices, not qp indices!
                     else if (hfb_blo_active .and. blo_p /= 0) then
                        if (iqp == blo_p .or. iqp == blo_pr) then
                           f = 0.5_dp  ! Equal filling
                        end if
                     end if

                     rho_ab = rho_ab + vp(iuv+ir)*vp(iuv+ic)*(1.0_dp-f) + up(iuv+ir)*up(iuv+ic)*f
                  end do

                  if (ns(i1) == ns(i2) .and. nl(i1) == nl(i2) .and. npar(i1) == npar(i2)) then
                     rho(:) = rho(:) + rho_ab*wf(:,i1)*wf(:,i2)
                  end if
               end if
            end do
         end do
      end do

      rho(:) = rho(:)*wdcori(:)  ! remove integration weights
   end subroutine hfb_density
