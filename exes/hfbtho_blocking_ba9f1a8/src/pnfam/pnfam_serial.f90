!------------------------------------------------------------------------------
! pnfam_serial.f90
!
! Charge-changing FAM code to compute pnQRPA strength functions or contours
! using the HFBTHO solutions.  The serial/OpenMP version.
!
! M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
!------------------------------------------------------------------------------
program pnfam
   use pnfam_common
   use interaction,  only : skip_residual_interaction
   use pnfam_solver, only : init_pnfam_solver, pnfam_solve
   implicit none

   !---------------------------------------------------------------------------
   ! Most variable and namelist declarations are in pnfam_common to reduce
   ! code duplication between the serial and the parallel version.
   !---------------------------------------------------------------------------
   
   ! Associate logging procedure pointers to the serial versions
   openlog => openlog_stdout
   writelog => writelog_stdout
   closelog => closelog_stdout
   abort => abort_nompi
   
   ! Read the input parameter files and initialize the interaction and the
   ! external field
   parameter_filename = "default.in"
   if (command_argument_count() == 1) then
      call get_command_argument(1, parameter_filename, status=ierr)
      if (ierr /= 0) then
         write(error_unit,*) "Error reading the command-line argument (too long?)"
         stop
      end if
   else
      if (command_argument_count() > 1) then
         write(error_unit,*) "More than one command-line argument! Aborting."
         stop
      end if
   end if
   
   ! Namelists, HFB solution, external field & cross-terms, and interaction
   call init_common

   ! Write the output header with all the input parameters
   call log_header
   
   ! Initialize the nonlinear solver
   call init_pnfam_solver(f=f%mat, g=g, bminus=f%beta_minus, vp=Vp, up=Up, vn=Vn, un=Un)
   
   ! Start the output file by logging the most important input parameters as comments
   call write_output_file_header
   
   
   select case(fam_mode)
      case("STR")
         
         call init_str_mode
   
         ! The main loop: iterate over a range of energies
         omega = energy_start
         do i=1, nr_points
      
            ! Solve the pnFAM equations for this energy
            call cpu_time(t_start)
            call pnfam_solve(omega, half_width, Ep, En, vp, up, vn, un, max_iter, &
               convergence_epsilon, skip_residual_interaction, f%beta_minus, str(:), iter)
            call cpu_time(t_finish)
            
            ! Log and save the result
            call log_str_point
            call write_output_str_point
            
            ! Next energy value
            omega = omega + energy_step
      
         end do
         
         call log_str_footer
   
      case("CONTOUR")
         
         ! Initialize the contour and fit the phase-space polynomial
         call init_contour_mode
         
         ! Compute as many points as symmetrization requires, allowing for
         ! skipping the first point if it has no phase space.
         ctr_strength(:,:) = 0
         do i=1, nr_compute
            if ((i == 1) .and. (aimag(ctr_z(i)) < epsilon(1.0)) &
            .and. (.not.rotate_contour) .and. (.not.include_endpoint)) then
               write(*,'(A,1F0.5,",",1F0.5,A)') 'N.b.: skipping point (', real(ctr_z(i)), &
                  aimag(ctr_z(i)), ') which should have no phase space'
               iter = 0; t_start = 0; t_finish = 0; str(:) = 0;  cstr(:) = 0
            else
               call cpu_time(t_start)
               call pnfam_solve(real(ctr_z(i),kind=dp), aimag(ctr_z(i)), Ep, En,            &
                  Vp, Up, Vn, Un, max_iter, convergence_epsilon, skip_residual_interaction, &
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
               convergence_epsilon, skip_residual_interaction, f%beta_minus, str(:), iter)
            if (iter == -1) then
               write(*,'(A,1ES12.6)') 'FAM was interrupted at iteration limit. si = ', fam_si
               stop ! No point in trying to estimate the derivative, if no convergence
            end if
            y2 = str(1)
            
            call pnfam_solve(omega-delta, half_width, Ep, En, vp, up, vn, un, max_iter, &
               convergence_epsilon, skip_residual_interaction, f%beta_minus, str(:), iter)
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
   
end program pnfam
