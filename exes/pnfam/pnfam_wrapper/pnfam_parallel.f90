!------------------------------------------------------------------------------
! pnfam_parallel.f90
!
! Charge-changing FAM code to compute pnQRPA strength functions or contours
! using the HFBTHO solutions.  The MPI/hybrid version.  Uses ADLB for load
! balancing.
!
! M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
!------------------------------------------------------------------------------
program pnfam
   use pnfam_common!, psi_serial => psi
   !use phasespace,   only : polyfit_eval, ratint_eval_psi
   use pnfam_solver, only : init_pnfam_solver, pnfam_solve
   !use constants,    only : ln2_kappa
   !use complex_quadrature
   use logger
   use logger_mpi
   implicit none
   !include 'mpif.h' now included in logger_mpi
#if(USE_ADLB==1)
   include 'adlbf.h'
#endif
   !---------------------------------------------------------------------------
   ! Most variable and namelist declarations are in pnfam_common to reduce
   ! code duplication between the serial and the parallel version.
   !---------------------------------------------------------------------------
   
   ! Namelist PARALLEL
   integer :: nr_adlb_servers, nr_parameter_files
   namelist /parallel/ nr_parameter_files, nr_adlb_servers
   
   ! Variables for parallellization
#if(USE_ADLB==1)
   integer :: work_handle(ADLB_HANDLE_SIZE)
#endif
   integer :: aprint_flag, am_server, am_debug_server, use_debug_server, &
      work_len, work_prio, answer_rank,   &
      work_type, req_types(16)
   integer :: rank, nrprocs, dqp, dmat, nuv, istatus(MPI_STATUS_SIZE)
   real(dp), parameter :: server_mem_max = 1d6
   integer, allocatable :: task_types(:)
   character(len=80) :: expected_mode
   character(len=80), allocatable :: parameter_filenames(:)
   integer :: comm_group, comm_all_workers, rank_in_group, &
      my_group, mytask, nr_procs_group
   real(dp) :: aux
      
   ! Mode STR
   real(dp), allocatable :: strengths(:,:), str_table(:,:), omegas(:)
   
   ! Mode CONTOUR
   type contour_task
      integer :: ix         ! index of the point
      complex(dp) :: omega  ! point on the complex plane
   end type contour_task
   type(contour_task) :: myctask, newtask
   complex(dp), allocatable :: ctr_strength_dist(:,:)
   
   ! Mode FINDMAX
   logical :: continue_loop
   
   ! Mode FINDSDRDEF
   integer :: irec, kfactor
   
   ! Initialize the MPI
   call MPI_INIT(ierr)
   if (ierr /= MPI_SUCCESS) stop "Fatal error: MPI initialization failed."
   call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, nrprocs, ierr)
   
   ! Set up the MPI log file procedure pointers
   openlog => openlog_mpi
   writelog => writelog_mpi
   closelog => closelog_mpi
   abort => abort_mpi
   
   ! Input parameters
   if (rank == 0) then
      if (command_argument_count() == 1) then
#if(USE_ADLB==1)
         nr_adlb_servers    = 1
#else
         nr_adlb_servers    = 0
#endif
         nr_parameter_files = 1
         allocate(parameter_filenames(1))
         call get_command_argument(1, parameter_filenames(1), status=ierr)
         if (ierr /= 0) then
            write(error_unit,*) "Error reading the command-line argument (too long?)"
            call MPI_ABORT(MPI_COMM_WORLD, -1, ierr)
         end if         
      else
         if (command_argument_count() > 1) then
            write(error_unit,*) "More than one command-line argument! Aborting."
            call MPI_ABORT(MPI_COMM_WORLD, -1, ierr)
         end if
         nr_adlb_servers    = -1
         nr_parameter_files = -1
         read(*,nml=parallel)
#if(USE_ADLB==0)
         if (nr_adlb_servers .ne. 0) then
            write(error_unit,*) "pnfam_mpi.x not compiled with adlb. nr_adlb_servers must be 0."
            call MPI_ABORT(MPI_COMM_WORLD, -1, ierr)
         end if
#endif
         if (nr_adlb_servers < 0 .or. nr_parameter_files < 1) then
            if (nr_adlb_servers < 0)    write(error_unit,*) "nr_adlb_servers cannot be < 0"
            if (nr_parameter_files < 1) write(error_unit,*) "nr_parameter_files cannot be < 1"
            call MPI_ABORT(MPI_COMM_WORLD, -1, ierr)
         end if

         allocate(parameter_filenames(nr_parameter_files))
         do i=1,nr_parameter_files
            read(*,*) parameter_filenames(i)
         end do
      end if
      expected_mode = get_fam_mode(parameter_filenames(1))
      call translate_uppercase(expected_mode)
   
      ! The resonance energy bisecting modes do not benefit from ADLB; automatically
      ! override nr_adlb_servers
      if (expected_mode == 'FINDMAX' .or. expected_mode == 'FINDMAXDEF' &
          .or. expected_mode == 'FINDSDR' .or. expected_mode == 'FINDSDRDEF') nr_adlb_servers = 0
      
      ! Check that we have enough MPI threads in our disposal
      ! The work distribution logic below likely fails otherwise
      if (nrprocs < (nr_adlb_servers + nr_parameter_files)) then
         write(error_unit,*) nrprocs, nr_adlb_servers, nr_parameter_files
         write(error_unit,*) "Too few MPI processes reserved; Need at least one per&
            & parameter file in addition to ADLB servers."
         call MPI_ABORT(MPI_COMM_WORLD, -1, ierr)
      end if
   end if
   
   ! Sync parameter file names (because we don't yet know who ends up reading what,
   ! everyone needs to get the memo) and nr_adlb_servers
   call MPI_BCAST(expected_mode, len(expected_mode), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(nr_parameter_files, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(nr_adlb_servers, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   if (nr_parameter_files > 1) then
      if (rank /= 0) allocate(parameter_filenames(nr_parameter_files))
      call MPI_BCAST(parameter_filenames, len(parameter_filenames(1))*nr_parameter_files, &
         MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
   end if
   
   ! If we want to use ADLB, initialize it
   ! If not, the MPI communicator for all workers is MPI_COMM_WORLD
#if(USE_ADLB==1)
   if (nr_adlb_servers > 0) then
      use_debug_server = 0  ! 0 = false in C
      aprint_flag = 0
      allocate(task_types(nr_parameter_files))
      do i=1,nr_parameter_files  ! Weird: using implied do here would segfault in ADLB_init.  Bug?
         task_types(i) = i
      end do
      
      call ADLB_init(nr_adlb_servers, use_debug_server, aprint_flag, nr_parameter_files, &
         task_types(:), am_server, am_debug_server, comm_all_workers, ierr)
      if (ierr /= ADLB_SUCCESS) call MPI_ABORT(MPI_COMM_WORLD, -1, ierr)
   endif
#endif
   if (nr_adlb_servers == 0) then
      call MPI_COMM_DUP(MPI_COMM_WORLD, comm_all_workers, ierr)
      am_server = 0
   end if
   
   ! Segregate the ADLB server(s) off the worker processes
#if(USE_ADLB==1)
   if (am_server /= 0) then
      call ADLB_server(server_mem_max , 0d0, ierr)
   end if
#endif
   if (am_server == 0) then
      ! If we have several input files, split further into worker groups, as
      ! evenly as possible (assigning color on a round-robin basis)
      if (nr_parameter_files > 1) then
         call MPI_COMM_RANK(comm_all_workers, rank_in_group, ierr)
         my_group = mod(rank_in_group,nr_parameter_files) + 1
         call MPI_COMM_SPLIT(comm_all_workers, my_group, 1, comm_group, ierr)
      else
         call MPI_COMM_DUP(comm_all_workers, comm_group, ierr)
         my_group = 1   ! We're #1! Just because we're the only group :(
      end if
      req_types(:) = my_group
      
      ! Get our rank amongst our group
      call MPI_COMM_RANK(comm_group, rank_in_group, ierr)
      
      ! If a log file was specified, use that in place of stdout
      ! In the serial code, this is not very useful, as you could just pipe the
      ! output to a file, but in the parallel code we want to separate each
      ! nucleus/operator/K
      call set_mpi_file_communicator(comm_group, rank_in_group)
      
      ! If our rank in the group is 0, we'll be responsible for the writing the header
      if (rank_in_group == 0) then
         
         ! Read the input parameter files and initialize the interaction and the
         ! external field
         parameter_filename = parameter_filenames(my_group)

         call init_common ! calls openlog...

         if (fam_mode /= expected_mode) then
            write(st,'(4A)') "Parameter file: ", trim(parameter_filename), &
               " is not for the expected mode ", trim(expected_mode)
            call writelog(st)
            call closelog
            call MPI_ABORT(MPI_COMM_WORLD, -1, ierr)
         end if

         ! Write the output header with all the input parameters
         call log_header
   
         ! Initialize the nonlinear solver
         call init_pnfam_solver(f=f%mat, g=g, bminus=f%beta_minus, vp=Vp, up=Up, vn=Vn, un=Un)
   
         ! Start the output file by logging the most important input parameters as comments
         call write_output_file_header
         
      else
         call openlog("")  ! only rank_in_group 0 knows the name
      end if
      
      ! Broadcast the data read from input files to all workers in the group
      call sync_initial_data(comm_group, rank_in_group)
      
      if (fam_mode == 'FINDMAXDEF') then
         K = rank_in_group
         if (rank_in_group > 0) then
            ! In FINDMAXDEF mode, we haven't initialized the external field beyond rank_in_group 0
            call init_external_field(label=operator_name, k=k, op=f)
            call init_pnfam_solver(f=f%mat, g=g, bminus=f%beta_minus, vp=Vp, up=Up, vn=Vn, un=Un)
         end if
      end if
      
      if (fam_mode == 'FINDSDR') then
         if (rank_in_group > 0) then
            ! In FINDSDR mode, we haven't initialized the external field beyond rank_in_group 0
            if (rank_in_group == 1) operator_name = 'RS1-'
            if (rank_in_group == 2) operator_name = 'RS2-'
            call init_external_field(label=operator_name, k=0, op=f)
            call init_pnfam_solver(f=f%mat, g=g, bminus=f%beta_minus, vp=Vp, up=Up, vn=Vn, un=Un)
         end if
      end if
      
      if (fam_mode == 'FINDSDRDEF') then
         if (rank_in_group > 0) then
            ! In FINDSDR mode, we haven't initialized the external field beyond rank_in_group 0
            if (rank_in_group == 1 .or. rank_in_group == 2) operator_name = 'RS1-'
            if (rank_in_group > 2) operator_name = 'RS2-'
            if (rank_in_group == 1 .or. rank_in_group == 4) k = 1
            if (rank_in_group == 5) k = 2
            call init_external_field(label=operator_name, k=k, op=f)
            call init_pnfam_solver(f=f%mat, g=g, bminus=f%beta_minus, vp=Vp, up=Up, vn=Vn, un=Un)
         end if
      end if
      
      if (fam_mode == 'FINDMAX' .or. fam_mode == 'FINDMAXDEF' .or. fam_mode == 'FINDSDR' &
         .or. fam_mode == 'FINDSDRDEF') then
         call MPI_BCAST(delta, 1, MPI_DOUBLE_PRECISION, 0, comm_group, ierr)
      end if
      
      ! A non-elegant safeguard:
      ! In the SDR modes, the operator in the input file is directly used to initialize
      ! rank #0; if it's not RS0-, our results would be completely wrong. Abort now
      ! if this is about to happen.
      if (fam_mode == 'FINDSDR' .or. fam_mode == 'FINDSDRDEF') then
         if (rank_in_group == 0) then
            if (operator_name /= 'RS0-') then
               call writelog("ERROR: SDR modes require RS0- in the input file!")
               call closelog
               call MPI_ABORT(MPI_COMM_WORLD, -1, ierr)
            end if
         end if
      end if
      
      ! At this point every worker should be in sync with their group, the ADLB
      ! servers (if any) are running, and we can fill up the work pool and
      ! then get computing
      select case(fam_mode)
         case("STR")
            
            ! Every worker in the group should initialize these the same way
            allocate(omegas(nr_points), strengths(nr_points,nxterms+1))
            strengths = 0;  omegas = 0
            omega = energy_start
            do i=1, nr_points
               omegas(i) = omega
               omega = omega + energy_step
            end do
            
            ! Only one worker fills the pool -- if ADLB is active
#if(USE_ADLB==1)
            if (nr_adlb_servers > 0 .and. rank_in_group == 0) then
               do i=1, nr_points
                  call ADLB_put(i, storage_size(i)/8, -1, -1, my_group, 5, ierr)
                  if (ierr /= ADLB_SUCCESS) call ADLB_ABORT(1, ierr)
               end do
            end if
#endif
      
            ! This barrier makes sure that the others wait for the pool to fill up
            call MPI_BARRIER(comm_group, ierr)
            
            if (rank_in_group == 0) call init_str_mode
            
            ! If ADLB disabled, prepare to assign points by round robin
            if (nr_adlb_servers == 0) then
               call MPI_COMM_SIZE(comm_group, nr_procs_group, ierr)
               i = rank_in_group + 1
            end if
            
            ! The main loop: request and complete tasks until there is none left
            do
#if(USE_ADLB==1)
               if (nr_adlb_servers > 0) then
                  ! If ADLB is active, get the point from the shared pool
                  call ADLB_RESERVE(req_types, work_type, work_prio, work_handle, work_len, answer_rank, ierr)
                  if (ierr /= ADLB_SUCCESS) exit
                  call ADLB_GET_RESERVED(mytask, work_handle, ierr)
                  if (ierr /= ADLB_SUCCESS) exit
                  
                  ! Assertion -- this should never happen, but if it did, our
                  ! results would be nonsense, so let's be careful and just abort
                  if (work_type /= my_group) then
                     write(st,*) "ERROR: RECEIVED WRONG TYPE OF WORK"
                     call writelog(st)
                     call closelog
                     call MPI_ABORT(MPI_COMM_WORLD, -1, ierr)
                  end if
               end if
#endif
               if (nr_adlb_servers == 0) then
                  ! If ADLB is disabled, choose the point in round robin
                  if (i > nr_points) exit
                  mytask = i
                  i = i + nr_procs_group
               end if
               
               call cpu_time(t_start)
         
               omega = omegas(mytask)
               call pnfam_solve(omega, half_width, Ep, En, vp, up, vn, un, max_iter, &
                  convergence_epsilon, quench_residual_int, f%beta_minus, str(:), iter)
               strengths(mytask,:) = str(:)
         
               call cpu_time(t_finish)
         
               ! Log the result
               call log_str_point
         
               ! At least let us know the convergence if there was a failure
               if (iter == -1) then
                  write(st,'(A,1ES12.6)') 'FAM was interrupted at iteration limit. si = ', fam_si
                  call writelog(st)
               end if
         
            end do
            
            ! Collect the results and write them out
            allocate(str_table(nr_points,nxterms+1))
            call MPI_REDUCE(strengths, str_table, nr_points*(nxterms+1), MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm_group, ierr)
            if (rank_in_group == 0) then
               do i=1,nr_points
                  ! Write the full result to the output file, including cross-terms
                  write(fno,'(1F14.6)',advance='no') omegas(i)
                  do istr=1, size(str)
                     write(fno,'(4X,1ES23.16)',advance='no') str_table(i,istr)
                  end do
                  write(fno,'(A)') ""
                  flush(fno)
               end do
               close(fno)        ! close the output file
            end if
            
            if (rank_in_group == 0) call log_str_footer
   
         case("CONTOUR")
         
            if (rank_in_group == 0) then
               ! Contour setup, but only on rank 0 (exit if Q < 0)
               call init_contour_mode(ierr)
               if (ierr /= 0) then
                  write(error_unit,*) "Q-value is negative. Aborting beta decay calculation."
                  call MPI_ABORT(MPI_COMM_WORLD, -1, ierr)
               end if
            
               ! Submit the computation to ADLB if active
#if(USE_ADLB==1)
               if (nr_adlb_servers > 0) then
                  do i=1, nr_compute
                     newtask%ix    = i
                     newtask%omega = ctr_z(i)
                     work_prio     = nr_compute+1-i
                     call ADLB_PUT(newtask, storage_size(newtask)/8, -1, -1, my_group, work_prio, ierr)
                     if (ierr /= ADLB_SUCCESS) call ADLB_ABORT(1, ierr)
                  end do
               end if
#endif
            end if
            ! If ADLB off, share initializations from init_contour_mode
            if (nr_adlb_servers == 0) then
               call MPI_BCAST(nr_compute, 1, MPI_INTEGER, 0, comm_group, ierr)
               call MPI_BCAST(ctr_bound_right, 1, MPI_DOUBLE_PRECISION, 0, comm_group, ierr)
               if (.not.allocated(ctr_z)) allocate(ctr_z(nr_points))
               call MPI_BCAST(ctr_z, nr_points, MPI_COMPLEX, 0, comm_group, ierr)
               call MPI_COMM_SIZE(comm_group, nr_procs_group, ierr)
               i = rank_in_group + 1
            end if
   
            call MPI_BARRIER(comm_group, ierr)
            allocate(ctr_strength_dist(nr_points, size(cstr)))
            if (.not.allocated(ctr_strength)) allocate(ctr_strength(nr_points, size(cstr)))
            
            ctr_strength = 0
            ctr_strength_dist = 0
            req_types(:) = my_group
            do
#if(USE_ADLB==1)
               if (nr_adlb_servers > 0) then
                  call ADLB_RESERVE(req_types, work_type, work_prio, work_handle, work_len, answer_rank, ierr)
                  if (ierr /= ADLB_SUCCESS) exit
                  call ADLB_GET_RESERVED(myctask, work_handle, ierr)
                  if (ierr /= ADLB_SUCCESS) exit
                  ! Assertion -- this should never happen, but if it did, our
                  ! results would be nonsense, so let's be careful and just abort
                  if (work_type /= my_group) then
                     write(st,*) "ERROR: RECEIVED WRONG TYPE OF WORK"
                     call writelog(st)
                     call closelog
                     call MPI_ABORT(MPI_COMM_WORLD, -1, ierr)
                  end if
               end if
#endif
               if (nr_adlb_servers == 0) then
                  ! If ADLB is disabled, choose the point in round robin
                  if (i > nr_compute) exit
                  myctask%ix = i
                  myctask%omega = ctr_z(i) ! this is nr_points long, but we only use the first nr_compute elements
                  i = i + nr_procs_group
               end if
            
               iter = 0; t_start = 0; t_finish = 0; str = 0; cstr = 0
               
               ! Compute as many points as symmetrization requires, allowing for
               ! skipping the point on the right bound real axis if it has no phase space.
               if ((abs(aimag(myctask%omega)) < epsilon(1.0)) &
                   .and. abs((ctr_bound_right - real(myctask%omega))) < epsilon(1.0) &
                   .and. (.not.rotate_contour) .and. (.not.include_endpoint)) then
                  write(st,'(A,1F0.5,",",1F0.5,A)') 'N.b.: skipping point (', real(myctask%omega), &
                     aimag(myctask%omega), ') which should have no phase space'
                  call writelog(st)
               else
                  call cpu_time(t_start)
                  call pnfam_solve(real(myctask%omega,kind=dp), aimag(myctask%omega), Ep, En,  &
                     Vp, Up, Vn, Un, max_iter, convergence_epsilon, quench_residual_int, &
                     f%beta_minus, str(:), iter, cstr(:))
                  call cpu_time(t_finish)
               end if
               wt_now = get_timer()
            
               ctr_strength_dist(myctask%ix,:) = cstr(:)
               call log_contour_point(myctask%omega, cstr(1))
            end do
         
            ! Collect to rank 0
            call MPI_REDUCE(ctr_strength_dist, ctr_strength, nr_points*size(cstr), &
               MPI_DOUBLE_COMPLEX, MPI_SUM, 0, comm_group, ierr)
         
            ! Finalize
            if (rank_in_group == 0) then            
               call contour_finalize_computation
               call log_contour_footer
               call write_output_contour_footer
               call contour_binary_output
            end if
   
         case("FINDMAX")
            
            ! This mode does not benefit from ADLB
            ! Furthermore, using more than two MPI process / group is a waste
            call MPI_COMM_SIZE(comm_group, nr_procs_group, ierr)
            if (nr_procs_group > 2) then
               call writelog("STOP: Wastefully using more than two MPI processes per file in FINDMAX.")
               call writelog("Aborting since this is likely unintentional.")
               call closelog
               call MPI_ABORT(MPI_COMM_WORLD, -1, ierr)
            end if
            
            if (rank_in_group == 0) call init_findmax_mode
         
            ! The main loop: try to find a maximum
            continue_loop = .false.
            omega = energy_start
            i = 1
            do
   
               ! Solve the pnFAM equations for two points and approximate the derivative
               call cpu_time(t_start)
               
               ! If we have two processes in the group, split the effort
               if (nr_procs_group == 2) then
                  if (rank_in_group == 0) then
                     call MPI_SEND(omega, 1, MPI_DOUBLE_PRECISION, 1, 1, comm_group, ierr)
                  else
                     call MPI_RECV(omega, 1, MPI_DOUBLE_PRECISION, 0, 1, comm_group, istatus, ierr)
                  end if
               end if
               
               ! Rank 0 in group always computes the first point for the derivative approx.
               if (rank_in_group == 0) then
                  call pnfam_solve(omega+delta, half_width, Ep, En, vp, up, vn, un, max_iter, &
                     convergence_epsilon, quench_residual_int, f%beta_minus, str(:), iter)
                  if (iter == -1) then
                     write(*,'(A,1ES12.6)') 'FAM was interrupted at iteration limit. si = ', fam_si
                     call MPI_ABORT(MPI_COMM_WORLD, -1, ierr) ! No point in trying to estimate the derivative, if no convergence
                  end if
                  y2 = str(1)
               end if
               
               ! The second point is computed by the other process, if present, otherwise by the first
               if (nr_procs_group == 1 .or. rank_in_group == 1) then
                  call pnfam_solve(omega-delta, half_width, Ep, En, vp, up, vn, un, max_iter, &
                     convergence_epsilon, quench_residual_int, f%beta_minus, str(:), iter)
                  if (iter == -1) then
                     write(*,'(A,1ES12.6)') 'FAM was interrupted at iteration limit. si = ', fam_si
                     call MPI_ABORT(MPI_COMM_WORLD, -1, ierr)
                  end if
                  y1 = str(1)
               end if
               
               ! If we're the second process (in the group), send the result back
               if (rank_in_group == 1) then
                  call MPI_SEND(y1, 1, MPI_DOUBLE_PRECISION, 0, 2, comm_group, ierr)
               else if (nr_procs_group == 2) then
                  call MPI_RECV(y1, 1, MPI_DOUBLE_PRECISION, 1, 2, comm_group, istatus, ierr)
               end if
               
               call cpu_time(t_finish)
               
               if (rank_in_group == 0) then
               
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
                  if (search_step < energy_tolerance .or. i > max_points) continue_loop = .true.
                  i = i + 1
                  
                  ! Let the other process (if it exists) know if we should exit the loop
                  if (nr_procs_group == 2) then
                     call MPI_SEND(continue_loop, 1, MPI_LOGICAL, 1, 3, comm_group, ierr)
                  end if
               else
                  call MPI_RECV(continue_loop, 1, MPI_LOGICAL, 0, 3, comm_group, istatus, ierr)
               end if
               
               if (continue_loop) exit
               
            end do
   
            if (rank_in_group == 0) then
               call write_output_findmax_result
               call log_findmax_footer
            end if
         
         case("FINDMAXDEF")
         
            ! This mode does not benefit from ADLB
            ! Furthermore, we require exactly two processes / group, because
            ! otherwise things would get _really_ complicated with two different K
            call MPI_COMM_SIZE(comm_group, nr_procs_group, ierr)
            if (nr_procs_group /= 2) then
               call writelog("STOP: FINDMAXDEF requires exactly two processes per group.")
               call writelog("Aborting.")
               call closelog
               call MPI_ABORT(MPI_COMM_WORLD, -1, ierr)
            end if
            
            if (rank_in_group == 0) call init_findmax_mode
            ! The main loop: try to find a maximum
            continue_loop = .false.
            omega = energy_start
            i = 1
            do
   
               ! Solve the pnFAM equations for two points and approximate the derivative
               call cpu_time(t_start)
               
               ! Split the effort: rank K in group computes projection K
               ! The rank 0 in group is in charge
               if (rank_in_group == 0) then
                  call MPI_SEND(omega, 1, MPI_DOUBLE_PRECISION, 1, 1, comm_group, ierr)
               else
                  call MPI_RECV(omega, 1, MPI_DOUBLE_PRECISION, 0, 1, comm_group, istatus, ierr)
               end if
               
               ! Computes the contribution from your K to the derivative approx. points
               
               call pnfam_solve(omega+delta, half_width, Ep, En, vp, up, vn, un, max_iter, &
                  convergence_epsilon, quench_residual_int, f%beta_minus, str(:), iter)
               if (iter == -1) then
                  write(*,'(A,1ES12.6)') 'FAM was interrupted at iteration limit. si = ', fam_si
                  call MPI_ABORT(MPI_COMM_WORLD, -1, ierr) ! No point in trying to estimate the derivative, if no convergence
               end if
               y2 = str(1)
            
               call pnfam_solve(omega-delta, half_width, Ep, En, vp, up, vn, un, max_iter, &
                  convergence_epsilon, quench_residual_int, f%beta_minus, str(:), iter)
               if (iter == -1) then
                  write(*,'(A,1ES12.6)') 'FAM was interrupted at iteration limit. si = ', fam_si
                  call MPI_ABORT(MPI_COMM_WORLD, -1, ierr)
               end if
               y1 = str(1)
               
               ! If we're the second process (in the group), send the result back
               if (rank_in_group == 1) then
                  call MPI_SEND(y1, 1, MPI_DOUBLE_PRECISION, 0, 2, comm_group, ierr)
                  call MPI_SEND(y2, 1, MPI_DOUBLE_PRECISION, 0, 4, comm_group, ierr)
               else
                  call MPI_RECV(aux, 1, MPI_DOUBLE_PRECISION, 1, 2, comm_group, istatus, ierr)
                  y1 = y1 + 2*aux
                  call MPI_RECV(aux, 1, MPI_DOUBLE_PRECISION, 1, 4, comm_group, istatus, ierr)
                  y2 = y2 + 2*aux
               end if
               
               call cpu_time(t_finish)
               
               if (rank_in_group == 0) then
               
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
                  if (search_step < energy_tolerance .or. i > max_points) continue_loop = .true.
                  i = i + 1
                  
                  ! Let the other process (if it exists) know if we should exit the loop
                  if (nr_procs_group == 2) then
                     call MPI_SEND(continue_loop, 1, MPI_LOGICAL, 1, 3, comm_group, ierr)
                  end if
               else
                  call MPI_RECV(continue_loop, 1, MPI_LOGICAL, 0, 3, comm_group, istatus, ierr)
               end if
               
               if (continue_loop) exit
               
            end do
   
            if (rank_in_group == 0) then
               call write_output_findmax_result
               call log_findmax_footer
            end if
            
         case("FINDSDR")
         
            ! This mode does not benefit from ADLB
            ! Furthermore, we require exactly three processes / group
            ! #0 computes RS0, #1 computes RS1, #2 computes RS2
            call MPI_COMM_SIZE(comm_group, nr_procs_group, ierr)
            if (nr_procs_group /= 3) then
               call writelog("STOP: FINDSDR requires exactly three processes per group.")
               call writelog("Aborting.")
               call closelog
               call MPI_ABORT(MPI_COMM_WORLD, -1, ierr)
            end if
            
            if (rank_in_group == 0) call init_findmax_mode
            ! The main loop: try to find a maximum
            continue_loop = .false.
            omega = energy_start
            i = 1
            do
   
               ! Solve the pnFAM equations for two points and approximate the derivative
               call cpu_time(t_start)
               
               ! Split the effort: rank J in group computes the operator [RS]_J
               ! The rank 0 in group is in charge
               if (rank_in_group == 0) then
                  call MPI_SEND(omega, 1, MPI_DOUBLE_PRECISION, 1, 1, comm_group, ierr)
                  call MPI_SEND(omega, 1, MPI_DOUBLE_PRECISION, 2, 1, comm_group, ierr)
               else
                  call MPI_RECV(omega, 1, MPI_DOUBLE_PRECISION, 0, 1, comm_group, istatus, ierr)
               end if
               
               ! Compute the contribution from your operator to the derivative approx. points
               
               call pnfam_solve(omega+delta, half_width, Ep, En, vp, up, vn, un, max_iter, &
                  convergence_epsilon, quench_residual_int, f%beta_minus, str(:), iter)
               if (iter == -1) then
                  write(*,'(A,1ES12.6)') 'FAM was interrupted at iteration limit. si = ', fam_si
                  call MPI_ABORT(MPI_COMM_WORLD, -1, ierr) ! No point in trying to estimate the derivative, if no convergence
               end if
               y2 = str(1)
            
               call pnfam_solve(omega-delta, half_width, Ep, En, vp, up, vn, un, max_iter, &
                  convergence_epsilon, quench_residual_int, f%beta_minus, str(:), iter)
               if (iter == -1) then
                  write(*,'(A,1ES12.6)') 'FAM was interrupted at iteration limit. si = ', fam_si
                  call MPI_ABORT(MPI_COMM_WORLD, -1, ierr)
               end if
               y1 = str(1)
               
               ! If we're not the first process (in the group), send the result back
               if (rank_in_group > 0) then
                  call MPI_SEND(y1, 1, MPI_DOUBLE_PRECISION, 0, 2, comm_group, ierr)
                  call MPI_SEND(y2, 1, MPI_DOUBLE_PRECISION, 0, 4, comm_group, ierr)
               else
                  ! Collect the results; factors 3 and 5 come from summing over K
                  call MPI_RECV(aux, 1, MPI_DOUBLE_PRECISION, 1, 2, comm_group, istatus, ierr)
                  y1 = y1 + 3*aux
                  call MPI_RECV(aux, 1, MPI_DOUBLE_PRECISION, 1, 4, comm_group, istatus, ierr)
                  y2 = y2 + 3*aux
                  call MPI_RECV(aux, 1, MPI_DOUBLE_PRECISION, 2, 2, comm_group, istatus, ierr)
                  y1 = y1 + 5*aux
                  call MPI_RECV(aux, 1, MPI_DOUBLE_PRECISION, 2, 4, comm_group, istatus, ierr)
                  y2 = y2 + 5*aux
               end if
               
               call cpu_time(t_finish)
               
               if (rank_in_group == 0) then
               
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
                  if (search_step < energy_tolerance .or. i > max_points) continue_loop = .true.
                  i = i + 1
                  
                  ! Let the other processes know if we should exit the loop
                  call MPI_SEND(continue_loop, 1, MPI_LOGICAL, 1, 3, comm_group, ierr)
                  call MPI_SEND(continue_loop, 1, MPI_LOGICAL, 2, 3, comm_group, ierr)
               else
                  call MPI_RECV(continue_loop, 1, MPI_LOGICAL, 0, 3, comm_group, istatus, ierr)
               end if
               
               if (continue_loop) exit
               
            end do
   
            if (rank_in_group == 0) then
               call write_output_findmax_result
               call log_findmax_footer
            end if
         
         case("FINDSDRDEF")
         
            ! This mode does not benefit from ADLB
            ! Furthermore, we require exactly six processes / group
            ! #0 computes RS0, #1-2 compute RS1, #3-5 compute RS2
            call MPI_COMM_SIZE(comm_group, nr_procs_group, ierr)
            if (nr_procs_group /= 6) then
               call writelog("STOP: FINDSDRDEF requires exactly six processes per group.")
               call writelog("Aborting.")
               call closelog
               call MPI_ABORT(MPI_COMM_WORLD, -1, ierr)
            end if
            
            if (rank_in_group == 0) call init_findmax_mode
            ! The main loop: try to find a maximum
            continue_loop = .false.
            omega = energy_start
            i = 1
            do
   
               ! Solve the pnFAM equations for two points and approximate the derivative
               call cpu_time(t_start)
               
               ! Split the effort: each rank computes one of the K's for one of
               ! the operators [RS]_J
               ! The rank 0 in group is in charge
               if (rank_in_group == 0) then
                  do irec=1,5
                     call MPI_SEND(omega, 1, MPI_DOUBLE_PRECISION, irec, 1, comm_group, ierr)
                  end do
               else
                  call MPI_RECV(omega, 1, MPI_DOUBLE_PRECISION, 0, 1, comm_group, istatus, ierr)
               end if
               
               ! Compute the contribution from your operator to the derivative approx. points
               
               call pnfam_solve(omega+delta, half_width, Ep, En, vp, up, vn, un, max_iter, &
                  convergence_epsilon, quench_residual_int, f%beta_minus, str(:), iter)
               if (iter == -1) then
                  write(*,'(A,1ES12.6)') 'FAM was interrupted at iteration limit. si = ', fam_si
                  call MPI_ABORT(MPI_COMM_WORLD, -1, ierr) ! No point in trying to estimate the derivative, if no convergence
               end if
               y2 = str(1)
            
               call pnfam_solve(omega-delta, half_width, Ep, En, vp, up, vn, un, max_iter, &
                  convergence_epsilon, quench_residual_int, f%beta_minus, str(:), iter)
               if (iter == -1) then
                  write(*,'(A,1ES12.6)') 'FAM was interrupted at iteration limit. si = ', fam_si
                  call MPI_ABORT(MPI_COMM_WORLD, -1, ierr)
               end if
               y1 = str(1)
               
               ! If we're not the first process (in the group), send the result back
               if (rank_in_group > 0) then
                  call MPI_SEND(y1, 1, MPI_DOUBLE_PRECISION, 0, 2, comm_group, ierr)
                  call MPI_SEND(y2, 1, MPI_DOUBLE_PRECISION, 0, 4, comm_group, ierr)
               else
                  ! Collect the results; note the factor of two for K>0 (ranks 2, 4, and 5)
                  do irec=1,5
                     kfactor = 1
                     if (rank_in_group == 2 .or. rank_in_group == 4 .or. rank_in_group == 5) kfactor = 2
                     call MPI_RECV(aux, 1, MPI_DOUBLE_PRECISION, irec, 2, comm_group, istatus, ierr)
                     y1 = y1 + kfactor*aux
                     call MPI_RECV(aux, 1, MPI_DOUBLE_PRECISION, irec, 4, comm_group, istatus, ierr)
                     y2 = y2 + kfactor*aux
                  end do
               end if
               
               call cpu_time(t_finish)
               
               if (rank_in_group == 0) then
               
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
                  if (search_step < energy_tolerance .or. i > max_points) continue_loop = .true.
                  i = i + 1
                  
                  ! Let the other processes know if we should exit the loop
                  do irec=1,5
                     call MPI_SEND(continue_loop, 1, MPI_LOGICAL, irec, 3, comm_group, ierr)
                  end do
               else
                  call MPI_RECV(continue_loop, 1, MPI_LOGICAL, 0, 3, comm_group, istatus, ierr)
               end if
               
               if (continue_loop) exit
               
            end do
   
            if (rank_in_group == 0) then
               call write_output_findmax_result
               call log_findmax_footer
            end if
      end select
      
      call closelog
      
   end if
#if(USE_ADLB==1)
   if (nr_adlb_servers > 0) call ADLB_FINALIZE(ierr)
#endif
   call MPI_FINALIZE(ierr)
   
contains
   
   ! The following subroutine relies on the visibility of the main program
   ! variables.  It's isolated here just to make the main program more readable,
   ! since there are so many MPI calls to be made.
   subroutine sync_initial_data(group_comm, rank)
      use hfb_solution
      use pnfam_interaction
      implicit none
      integer :: ixt
      integer, intent(in) :: group_comm, rank
      
      ! Master sends data to other processes, which now need to allocate all needed
      ! arrays and data structures
      ! --- start syncronization ---
      call MPI_BCAST(fam_mode, len(fam_mode), MPI_CHARACTER, 0, group_comm, ierr)
      call MPI_BCAST(nr_points, 1, MPI_INTEGER, 0, group_comm, ierr)
      call MPI_BCAST(energy_start, 1, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      call MPI_BCAST(energy_step, 1, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      call MPI_BCAST(max_iter, 1, MPI_INTEGER, 0, group_comm, ierr)
      call MPI_BCAST(convergence_epsilon, 1, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      call MPI_BCAST(half_width, 1, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      call MPI_BCAST(quench_residual_int, 1, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      call MPI_BCAST(rot_correction, 1, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      call MPI_BCAST(include_endpoint, 1, MPI_LOGICAL, 0, group_comm, ierr)
      call MPI_BCAST(rotate_contour, 1, MPI_LOGICAL, 0, group_comm, ierr)
      call MPI_BCAST(energy_shift_prot, 1, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      call MPI_BCAST(energy_shift_neut, 1, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      
      ! Timer
      call MPI_BCAST(wt_start, 1, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
   
      call MPI_BCAST(nb, 1, MPI_INTEGER, 0, group_comm, ierr)
      if (rank > 0) allocate(db(nb))
      call MPI_BCAST(db, nb, MPI_INTEGER, 0, group_comm, ierr)
      dqp = sum(db); dmat = sum(db*db)
   
      call MPI_BCAST(nghl, 1, MPI_INTEGER, 0, group_comm, ierr)
      if (rank > 0) then
         call init_blockmatrix_type(nb,db)
         allocate(Ep(dqp), En(dqp))
         call allocate_blockmatrix(Up, dmat)
         call allocate_blockmatrix(Vp, dmat)
         call allocate_blockmatrix(Un, dmat)
         call allocate_blockmatrix(Vn, dmat)
         allocate(nr(dqp), nz(dqp), nl(dqp), ns(dqp), npar(dqp))
         allocate(wdcori(nghl), y(nghl), z(nghl))
         allocate(wf(nghl,dqp), wfdr(nghl,dqp), wfdz(nghl,dqp), wfd2(nghl,dqp))
         allocate(imstart(nb),isstart(nb))
#ifdef USE_HBLAS
#if USE_HBLAS==1
      allocate(num_spin_up(nb))
      allocate(wfdp(nghl,dqp), wfd2_all(nghl,dqp))
#endif
#endif
      end if
      call MPI_BCAST(Ep, dqp, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      call MPI_BCAST(En, dqp, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      call MPI_BCAST(Up%ir2c, nb, MPI_INTEGER, 0, group_comm, ierr)
      call MPI_BCAST(Up%ic2r, nb, MPI_INTEGER, 0, group_comm, ierr)
      call MPI_BCAST(Up%ir2m, nb, MPI_INTEGER, 0, group_comm, ierr)
      call MPI_BCAST(Up%ic2m, nb, MPI_INTEGER, 0, group_comm, ierr)
      call MPI_BCAST(Up%elem, dmat, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      call MPI_BCAST(Un%ir2c, nb, MPI_INTEGER, 0, group_comm, ierr)
      call MPI_BCAST(Un%ic2r, nb, MPI_INTEGER, 0, group_comm, ierr)
      call MPI_BCAST(Un%ir2m, nb, MPI_INTEGER, 0, group_comm, ierr)
      call MPI_BCAST(Un%ic2m, nb, MPI_INTEGER, 0, group_comm, ierr)
      call MPI_BCAST(Un%elem, dmat, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      call MPI_BCAST(Vp%ir2c, nb, MPI_INTEGER, 0, group_comm, ierr)
      call MPI_BCAST(Vp%ic2r, nb, MPI_INTEGER, 0, group_comm, ierr)
      call MPI_BCAST(Vp%ir2m, nb, MPI_INTEGER, 0, group_comm, ierr)
      call MPI_BCAST(Vp%ic2m, nb, MPI_INTEGER, 0, group_comm, ierr)
      call MPI_BCAST(Vp%elem, dmat, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      call MPI_BCAST(Vn%ir2c, nb, MPI_INTEGER, 0, group_comm, ierr)
      call MPI_BCAST(Vn%ic2r, nb, MPI_INTEGER, 0, group_comm, ierr)
      call MPI_BCAST(Vn%ir2m, nb, MPI_INTEGER, 0, group_comm, ierr)
      call MPI_BCAST(Vn%ic2m, nb, MPI_INTEGER, 0, group_comm, ierr)
      call MPI_BCAST(Vn%elem, dmat, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      call MPI_BCAST(nr, dqp, MPI_INTEGER, 0, group_comm, ierr)
      call MPI_BCAST(nz, dqp, MPI_INTEGER, 0, group_comm, ierr)
      call MPI_BCAST(nl, dqp, MPI_INTEGER, 0, group_comm, ierr)
      call MPI_BCAST(ns, dqp, MPI_INTEGER, 0, group_comm, ierr)
      call MPI_BCAST(npar, dqp, MPI_INTEGER, 0, group_comm, ierr)
      call MPI_BCAST(wdcori, nghl, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      call MPI_BCAST(y, nghl, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      call MPI_BCAST(z, nghl, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      call MPI_BCAST(wf, nghl*dqp, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      call MPI_BCAST(wfdr, nghl*dqp, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      call MPI_BCAST(wfdz, nghl*dqp, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      call MPI_BCAST(wfd2, nghl*dqp, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      call MPI_BCAST(imstart, nb, MPI_INTEGER, 0, group_comm, ierr)
      call MPI_BCAST(isstart, nb, MPI_INTEGER, 0, group_comm, ierr)
#ifdef USE_HBLAS
#if USE_HBLAS==1
      call MPI_BCAST(num_spin_up, nb, MPI_INTEGER, 0, group_comm, ierr)
      call MPI_BCAST(wfdp, nghl*dqp, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      call MPI_BCAST(wfd2_all, nghl*dqp, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
#endif
#endif
      
      ! Pairing window
      if (rank > 0) then
         allocate(pwi_start_n(nb), pwi_start_p(nb), pwi_dim_n(nb), pwi_dim_p(nb))
         allocate(pwi_qp_p(dqp), pwi_uv_p(dqp), pwi_qp_n(dqp), pwi_uv_n(dqp))
         allocate(pwi_active_n(dqp), pwi_active_p(dqp))
      end if
      call MPI_BCAST(pwi_start_n,   nb, MPI_INTEGER, 0, group_comm, ierr)
      call MPI_BCAST(pwi_start_p,   nb, MPI_INTEGER, 0, group_comm, ierr)
      call MPI_BCAST(pwi_dim_n,     nb, MPI_INTEGER, 0, group_comm, ierr)
      call MPI_BCAST(pwi_dim_p,     nb, MPI_INTEGER, 0, group_comm, ierr)
      call MPI_BCAST(pwi_qp_n,     dqp, MPI_INTEGER, 0, group_comm, ierr)
      call MPI_BCAST(pwi_qp_p,     dqp, MPI_INTEGER, 0, group_comm, ierr)
      call MPI_BCAST(pwi_uv_n,     dqp, MPI_INTEGER, 0, group_comm, ierr)
      call MPI_BCAST(pwi_uv_p,     dqp, MPI_INTEGER, 0, group_comm, ierr)
      call MPI_BCAST(pwi_active_n, dqp, MPI_LOGICAL, 0, group_comm, ierr)
      call MPI_BCAST(pwi_active_p, dqp, MPI_LOGICAL, 0, group_comm, ierr)
      
      ! Finite temperature
      call MPI_BCAST(ft_active,    1, MPI_LOGICAL,          0, group_comm, ierr)
      if (ft_active) then
         call MPI_BCAST(ft_temp,      1, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
         call MPI_BCAST(ft_entropy,   3, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
         if (rank > 0) allocate(ft_fp(dqp),ft_fn(dqp))
         call MPI_BCAST(ft_fp,      dqp, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
         call MPI_BCAST(ft_fn,      dqp, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      end if
      
      ! Odd nucleus
      call MPI_BCAST(hfb_blo_active, 1, MPI_LOGICAL, 0, comm_group, ierr)
      call MPI_BCAST(hfb_blo_qp,     2, MPI_INTEGER, 0, comm_group, ierr)
   
      ! External field; not broadcasted in the FINDMAXDEF mode, because there
      ! different workers have a different K
      if (expected_mode /= 'FINDMAXDEF' .and. expected_mode /= 'FINDSDR' &
         .and. expected_mode /= 'FINDSDRDEF') then
         call MPI_BCAST(f%label,       80, MPI_CHARACTER, 0, group_comm, ierr)
         call MPI_BCAST(f%k,            1, MPI_INTEGER,   0, group_comm, ierr)
         call MPI_BCAST(f%beta_minus,   1, MPI_LOGICAL,   0, group_comm, ierr)
         call MPI_BCAST(f%parity_even,  1, MPI_LOGICAL,   0, group_comm, ierr)
         if (rank == 0) nuv = size(f%mat%elem)
         call MPI_BCAST(nuv, 1, MPI_INTEGER, 0, group_comm, ierr)
         if (rank > 0) call allocate_blockmatrix(f%mat, nuv)
         call MPI_BCAST(f%mat%ir2c, nb, MPI_INTEGER, 0, group_comm, ierr)
         call MPI_BCAST(f%mat%ic2r, nb, MPI_INTEGER, 0, group_comm, ierr)
         call MPI_BCAST(f%mat%ir2m, nb, MPI_INTEGER, 0, group_comm, ierr)
         call MPI_BCAST(f%mat%ic2m, nb, MPI_INTEGER, 0, group_comm, ierr)
         call MPI_BCAST(f%mat%elem, nuv, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      else
         call MPI_BCAST(operator_name, len(operator_name), MPI_CHARACTER, 0, group_comm, ierr)
      end if

      ! Cross-terms fields
      call MPI_BCAST(nxterms, 1, MPI_INTEGER, 0, group_comm, ierr)
      if ((rank > 0) .and. nxterms > 0) allocate(g(1:nxterms))
   
      do ixt=1, nxterms
         call MPI_BCAST(g(ixt)%label,       80, MPI_CHARACTER, 0, group_comm, ierr)
         call MPI_BCAST(g(ixt)%k,            1, MPI_INTEGER,   0, group_comm, ierr)
         call MPI_BCAST(g(ixt)%beta_minus,   1, MPI_LOGICAL,   0, group_comm, ierr)
         call MPI_BCAST(g(ixt)%parity_even,  1, MPI_LOGICAL,   0, group_comm, ierr)
         if (rank > 0) call allocate_blockmatrix(g(ixt)%mat, nuv)
         call MPI_BCAST(g(ixt)%mat%ir2c, nb, MPI_INTEGER, 0, group_comm, ierr)
         call MPI_BCAST(g(ixt)%mat%ic2r, nb, MPI_INTEGER, 0, group_comm, ierr)
         call MPI_BCAST(g(ixt)%mat%ir2m, nb, MPI_INTEGER, 0, group_comm, ierr)
         call MPI_BCAST(g(ixt)%mat%ic2m, nb, MPI_INTEGER, 0, group_comm, ierr)
         call MPI_BCAST(g(ixt)%mat%elem, nuv, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)      
      end do
   
      if (rank > 0) allocate(str(1:nxterms+1))
      if (rank > 0) allocate(cstr(1:nxterms+1))
      if ((expected_mode /= 'FINDMAXDEF' .and. expected_mode /= 'FINDSDR' &
         .and. expected_mode /= 'FINDSDRDEF') .and. (rank > 0)) &
         call init_pnfam_solver(f=f%mat, g=g, bminus=f%beta_minus, vp=Vp, up=Up, vn=Vn, un=Un)

      if (rank > 0) allocate(crho(nghl), cs(nghl), cpair(nghl), cspair(nghl))
      call MPI_BCAST(crho, nghl, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      call MPI_BCAST(cs, nghl, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      call MPI_BCAST(cpair, nghl, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      call MPI_BCAST(cspair, nghl, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      call MPI_BCAST(cdrho, 1, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      call MPI_BCAST(ctau, 1, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      call MPI_BCAST(ctj1, 1, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      call MPI_BCAST(ctj2, 1, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      call MPI_BCAST(crdj, 1, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      call MPI_BCAST(cds, 1, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      call MPI_BCAST(ct, 1, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      call MPI_BCAST(cj, 1, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      call MPI_BCAST(cgs, 1, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      call MPI_BCAST(ctj0, 1, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      call MPI_BCAST(cf, 1, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
      call MPI_BCAST(csdj, 1, MPI_DOUBLE_PRECISION, 0, group_comm, ierr)
   
      call MPI_BCAST(broyden_history_size, 1, MPI_INTEGER, 0, group_comm, ierr)
   
   end subroutine sync_initial_data
   
end program pnfam
