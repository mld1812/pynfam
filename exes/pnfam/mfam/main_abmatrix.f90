!------------------------------------------------------------------------------
! pnfam_matrix/main_abmatrix.f90
!
! T. Shafer, UNC Chapel Hill, 2015
!------------------------------------------------------------------------------
program pnfam_matrix_ab
!$ use omp_lib,             only : omp_get_max_threads, omp_get_wtime
   use logger
   use blockmatrix_type,    only : blockmatrix, allocate_blockmatrix, deallocate_blockmatrix
   use hfbtho_basis,        only : get_hfbtho_solution
   use hfb_solution_type,   only : hfb_solution
   use interaction,         only : init_skyrme, skip_residual_interaction
   use external_field_type, only : external_field
   use extfield,            only : deallocate_extfield
   use config_type,         only : matrix_config, read_namelist_from_file
   use pnfam_matrix,        only : init_sp_operator, init_qp_operator, qp_reverse,        &
                                   input_file_from_cli, make_two_qp_basis, ab_file_count, &
                                   mpi_open_ab_files, get_qp_indices, mpi_close_ab_files, &
                                   log_time, residual_interaction, time_reverse_basis,    &
                                   write_main_header

   use mpi

   implicit none
   
   integer, parameter :: si = kind(1)
   integer, parameter :: di = selected_int_kind(11)
   integer, parameter :: dp = kind(1.0d0)
   
   ! MPI and OpenMP
   !include 'mpif.h'
   integer :: myrank, nprocs, nthreads
   integer :: ierr
   
   ! Configuration
   type(matrix_config) :: config
   type(hfb_solution)  :: hfb_soln
   
   ! AB matrix storage
   integer :: recl, nfiles
   
   ! Two-quasiparticle basis
   integer, allocatable :: twoqp_20(:,:), twoqp_02(:,:)
   integer :: num_ab
   
   ! Reference operators
   type(external_field) :: f
   type(blockmatrix)    :: F20, F02
   integer, allocatable :: f20_lookup(:,:), f02_lookup(:,:)
   
   ! Calculation
   integer :: icalc, istart, iend, ic, ir, iqpp, iqpn, ik
   integer, allocatable :: fh(:)
   real(dp) :: time_start, time_end, time_total_start, time_total_end
   real(dp), allocatable :: H20(:), matcol(:)
   
   ! Misc
   character(len=210) :: infile, ms
   
   !----------------------------------------------------------------------------
   !----------------------------------------------------------------------------

   ! Init
   call mpi_init(ierr)
   call mpi_comm_rank(mpi_comm_world, myrank, ierr)
   call mpi_comm_size(mpi_comm_world, nprocs, ierr)
   
   ! Logging
   openlog  => openlog_stdout
   writelog => writelog_stdout
   closelog => closelog_stdout
   call openlog('')
   
   ! Set up the calculation on rank 0
   if (myrank == 0) then
      
      ! Check for correct integer precision
      if (di == -1) then
         write(*,'(1x,a)') 'ERROR! this machine does not support integers of&
            & selected_int_kind(11)!'
         call mpi_abort(mpi_comm_world, -1, ierr)
      end if
      
      ! Configuration and setup
      call input_file_from_cli(infile, default='matrix.in')
      
      open(1, file=trim(infile), status='old', access='sequential', iostat=ierr)
      if (ierr /= 0) stop 'Error opening Namelist file'
      
      call read_namelist_from_file(fh=1, config=config)
      
      call get_hfbtho_solution(fn=trim(config%file_hfb), &
         ep=hfb_soln%ep, vp=hfb_soln%vp, up=hfb_soln%up, &
         en=hfb_soln%en, vn=hfb_soln%vn, un=hfb_soln%un)
      
      call init_skyrme(fh=1, vp=hfb_soln%vp%elem, vn=hfb_soln%vn%elem, &
         up=hfb_soln%up%elem, un=hfb_soln%un%elem)
      
      close(1)
      
      ! Header, etc.
      call write_main_header(config, infile=infile, nmpi=nprocs, nomp=get_num_omp_threads(), &
         label='QRPA Matrix Construction')

      ! Make the *FORWARD* two-quasiparticle basis
      call make_two_qp_basis(class="20", k=config%basis_k, parity=config%basis_parity, &
         ecut=config%basis_cutoff, hfb=hfb_soln, basis_out=twoqp_20)
      
      ! IMPORTANT: the reverse (02) basis is always simply the time-reversed (20) basis!
      call time_reverse_basis(in=twoqp_20, out=twoqp_02)
      
      nthreads = get_num_omp_threads()
      num_ab   = size(twoqp_20, 2)
      nfiles   = ab_file_count(num_ab)
      
      inquire(iolength=recl) 1.0_dp
      
      write(*,'(a)')     'AB MATRIX PARAMETERS'
      write(*,'(a)')     repeat('-', 33)
      write(*,'(a,1i0)') 'Two-QP dimension ........ : ', num_ab
      write(*,'(a,1i0)') 'No. Calculations ........ : ', 2*num_ab
      write(*,'(a,1i0)') 'Record length (DP) ...... : ', recl
      write(*,'(a,1i0)') 'No. output files ........ : ', nfiles
      
      ! Notice re: odd orbital exclusion
      if (config%odd_block_n(1) /= 0 .or. config%odd_block_p(1) /= 0) then
         write(*,'(/,a)') "Notice: all orbitals are always included in calculation of&
            & A and B despite settings block_n and block_p."
      end if
   end if
   
   
   ! Distribute everything
   if (myrank == 0) then
      write(*,*) 
      call log_time('Distributing calculation details')
   end if
   call distribute_calculation_details
 
   if (myrank == 0) call log_time('Distributing PNFAM arrays')
   call distribute_pnfam_quantities
   call mpi_barrier(mpi_comm_world, ierr)
   
   
   ! Initialize reference operator and reverse-lookup arrays for F(20) and F(02)
   ! These are needed to calculate A and B
   call init_sp_operator(k=config%basis_k, parity=config%basis_parity, f_out=f)
   call init_qp_operator(class='20', f_sp=f, hfb=hfb_soln, f_qp_out=F20)
   call init_qp_operator(class='02', f_sp=f, hfb=hfb_soln, f_qp_out=F02)
   call get_qp_indices(F20, f20_lookup)
   call get_qp_indices(F02, f02_lookup)
 
   
   ! Each rank gets its own start and end indices
   call get_calculation_indices(dim=2*num_ab, ncpu=nprocs, is=istart, ie=iend)
   
   allocate(matcol(num_ab))
   matcol = 0
   
   ! Each rank opens the file for direct access
   !call mpi_open_ab_files(num_ab, config%file_basename, mpi_comm_world, &
   !   mpi_mode_wronly+mpi_mode_create+mpi_mode_excl, fh)
   
   call mpi_barrier(mpi_comm_world, ierr)
   
   if (myrank == 0) call log_time('Beginning the calculation')
   
   
   !----------------------------------------------------------------------------
   call get_timer(time_total_start)
   
   ! A matrix if odd, B matrix if even
   do icalc=istart, iend
      
      ! A matrix
      if (mod(icalc,2) == 1) then   
         ic = 1+(icalc-1)/2
         matcol(:) = 0
         
         call get_timer(time_start)
         
         ! Q.P. enegies on the diagonal
         iqpn = twoqp_20(1, ic)
         iqpp = twoqp_20(2, ic)
         matcol(ic) = hfb_soln%en(iqpn) + hfb_soln%ep(iqpp)
         
         ! Interaction
         if (.not.skip_residual_interaction) then
            ! Calculate the derivative
            ik = f20_lookup(iqpn, iqpp)
            if (ik == 0) then
               write(*,'(a)') 'ERROR: ic was unexpectedly zero entering residual_interaction'
               call mpi_abort(mpi_comm_world, -1, ierr)
            end if
            
            call residual_interaction(class='X', hfb=hfb_soln, ik=ik, delta_r_qp=F20, &
               f_sp=f%mat, f20_qp=F20, f02_qp=F02, h20_vector=H20)
            
            ! Assign the matrix elements of H20 to rows of the A matrix
            if (size(H20) == 0) then
               write(*,'(1x,a,1i0,a)') 'ERROR: H20(A) was zero unexpectedly on rank ', myrank, '!'
               call mpi_abort(mpi_comm_world, -1, ierr)
            else
               do ir=1, num_ab
                  iqpn = twoqp_20(1, ir)
                  iqpp = twoqp_20(2, ir)
                  ik = f20_lookup(iqpn, iqpp)
                  
                  if (ik == 0) then
                     write(*,'(1x,a,1i0,a)') 'WARNING: ik(A) was zero unexpectedly on rank ', myrank, '!'
                     cycle
                  end if
                  matcol(ir) = matcol(ir) + H20(ik)
               end do
            end if
         end if
         
         ! Store
         call mpi_store_ab_lower(matcol(ic:num_ab), ic, num_ab, 'A', fh(:))
         call get_timer(time_end)
         
         if (ic == 1) then
            write(ms,'("A matrix vector time: ",1f0.3," s")') time_end-time_start
            call log_time(trim(ms))
            write(ms,'("Estimated calculation time: ",1f0.3," minutes")') (time_end-time_start)*2*num_ab/nprocs/60.0_dp
            call log_time(trim(ms))
         end if
   
   
      ! B matrix
      else if (mod(icalc,2) == 0) then
         ic = icalc/2
         matcol(:) = 0
         
         call get_timer(time_start)         
         
         ! Interaction
         if (.not.skip_residual_interaction) then
            ! Calculate the derivative
            iqpn = twoqp_02(1, ic)
            iqpp = twoqp_02(2, ic)
            ik = f02_lookup(iqpn, iqpp)
            if (ik == 0) then
               write(*,'(a)') 'ERROR: ic was unexpectedly zero entering residual_interaction'
               call mpi_abort(mpi_comm_world, -1, ierr)
            end if
            
            call residual_interaction(class='Y', hfb=hfb_soln, ik=ik, delta_r_qp=F02, &
               f_sp=f%mat, f20_qp=F20, f02_qp=F02, h20_vector=H20)
            
            ! Assign the matrix elements of H20 to rows of the A matrix
            if (size(H20) == 0) then
               write(*,'(1x,a,1i0,a)') 'ERROR: H20(A) was zero unexpectedly on rank ', myrank, '!'
               call mpi_abort(mpi_comm_world, -1, ierr)
            else
               do ir=1, num_ab
                  iqpn = twoqp_20(1, ir)
                  iqpp = twoqp_20(2, ir)
                  ik = f20_lookup(iqpn, iqpp)
                  
                  if (ik == 0) then
                     write(*,'(1x,a,1i0,a)') 'WARNING: ik(A) was zero unexpectedly on rank ', myrank, '!'
                     cycle
                  end if
                  matcol(ir) = matcol(ir) + H20(ik)
               end do
            end if
         end if
         
         ! Store
         call mpi_store_ab_lower(matcol(ic:num_ab), ic, num_ab, 'B', fh(:))
         call get_timer(time_end)
         
         if (ic == 1) then
            write(ms,'("B matrix vector time: ",1f0.3," s")') time_end-time_start
            call log_time(trim(ms))
         end if
      end if
   end do
   
   call get_timer(time_total_end)
   
   write(ms,'("Rank ",1i5," completed after ",1f9.3," seconds.")') myrank, time_total_end-time_total_start
   call log_time(trim(ms))
   !----------------------------------------------------------------------------
   
   
   call mpi_close_ab_files(fh)
   call mpi_finalize(ierr)
 
contains
   
   !----------------------------------------------------------------------------
   ! MPI_BCAST everything from Rank 0
   !----------------------------------------------------------------------------
   subroutine distribute_calculation_details
      implicit none
      
      ! Configuration
      call mpi_bcast(config%file_hfb,              200, mpi_character,        0, mpi_comm_world, ierr)
      call mpi_bcast(config%file_basename,         200, mpi_character,        0, mpi_comm_world, ierr)
      call mpi_bcast(config%file_ab,               200, mpi_character,        0, mpi_comm_world, ierr)
      call mpi_bcast(config%file_xy,               200, mpi_character,        0, mpi_comm_world, ierr)
      call mpi_bcast(config%file_me_one,           200, mpi_character,        0, mpi_comm_world, ierr)
      call mpi_bcast(config%file_me_two,           200, mpi_character,        0, mpi_comm_world, ierr)
      call mpi_bcast(config%file_vv,               200, mpi_character,        0, mpi_comm_world, ierr)
      call mpi_bcast(config%basis_k,               1,   mpi_integer,          0, mpi_comm_world, ierr)
      call mpi_bcast(config%basis_parity,          1,   mpi_integer,          0, mpi_comm_world, ierr)
      call mpi_bcast(config%basis_cutoff,          1,   mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(config%diag_ev_limit,         1,   mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(config%diag_blocking,         1,   mpi_integer,          0, mpi_comm_world, ierr)
      call mpi_bcast(config%diag_nprow,            1,   mpi_integer,          0, mpi_comm_world, ierr)
      call mpi_bcast(config%me_check_norm,         1,   mpi_logical,          0, mpi_comm_world, ierr)
      call mpi_bcast(config%me_check_complete,     1,   mpi_logical,          0, mpi_comm_world, ierr)
      call mpi_bcast(config%me_disable_correction, 1,   mpi_logical,          0, mpi_comm_world, ierr)
      call mpi_bcast(config%odd_block_n,           2,   mpi_integer,          0, mpi_comm_world, ierr)
      call mpi_bcast(config%odd_block_p,           2,   mpi_integer,          0, mpi_comm_world, ierr)
      call mpi_bcast(config%odd_even_basename,     200, mpi_character,        0, mpi_comm_world, ierr)
      call mpi_bcast(config%odd_even_file_ab,      200, mpi_character,        0, mpi_comm_world, ierr)
      call mpi_bcast(config%odd_even_file_xy,      200, mpi_character,        0, mpi_comm_world, ierr)
      
      call mpi_bcast(nthreads, 1, mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast(num_ab,   1, mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast(recl,     1, mpi_integer, 0, mpi_comm_world, ierr)
      
      if (myrank /= 0) allocate(twoqp_20(2, num_ab), twoqp_02(2, num_ab))
      call mpi_bcast(twoqp_20, 2*num_ab, mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast(twoqp_02, 2*num_ab, mpi_integer, 0, mpi_comm_world, ierr)
      
   end subroutine distribute_calculation_details
   
   !----------------------------------------------------------------------------
   ! MPI_BCAST everything for pnFAM from Rank 0 
   !----------------------------------------------------------------------------
   subroutine distribute_pnfam_quantities
      use blockmatrix_type
      use hfbtho_basis
      use interaction
      
      implicit none
      
      integer :: nqp, nme
      
      ! HFBTHO
      !-------------------------------------------------------------------------
      ! Blocks
      call mpi_bcast(nb,  1, mpi_integer, 0, mpi_comm_world, ierr)
      if (myrank /= 0) allocate(db(nb))
      call mpi_bcast(db, nb, mpi_integer, 0, mpi_comm_world, ierr)
      
      nqp = sum(db)
      nme = sum(db**2)
      
      ! Energies
      if (myrank /= 0) allocate(hfb_soln%ep(nqp), hfb_soln%en(nqp))
      call mpi_bcast(hfb_soln%ep, nqp, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(hfb_soln%en, nqp, mpi_double_precision, 0, mpi_comm_world, ierr)
      ! Matrices
      if (myrank /= 0) then
         call allocate_blockmatrix(hfb_soln%up, nme)
         call allocate_blockmatrix(hfb_soln%un, nme)
         call allocate_blockmatrix(hfb_soln%vp, nme)
         call allocate_blockmatrix(hfb_soln%vn, nme)
      end if
      call mpi_bcast(hfb_soln%up%ir2c,  nb,          mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast(hfb_soln%up%ic2r,  nb,          mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast(hfb_soln%up%ir2m,  nb,          mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast(hfb_soln%up%ic2m,  nb,          mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast(hfb_soln%up%elem, nme, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(hfb_soln%un%ir2c,  nb,          mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast(hfb_soln%un%ic2r,  nb,          mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast(hfb_soln%un%ir2m,  nb,          mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast(hfb_soln%un%ic2m,  nb,          mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast(hfb_soln%un%elem, nme, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(hfb_soln%vp%ir2c,  nb,          mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast(hfb_soln%vp%ic2r,  nb,          mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast(hfb_soln%vp%ir2m,  nb,          mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast(hfb_soln%vp%ic2m,  nb,          mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast(hfb_soln%vp%elem, nme, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(hfb_soln%vn%ir2c,  nb,          mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast(hfb_soln%vn%ic2r,  nb,          mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast(hfb_soln%vn%ir2m,  nb,          mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast(hfb_soln%vn%ic2m,  nb,          mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast(hfb_soln%vn%elem, nme, mpi_double_precision, 0, mpi_comm_world, ierr)
      ! Arrays
      if (myrank /= 0) then
         allocate(nr(nqp), nz(nqp), nl(nqp), ns(nqp), npar(nqp))
         allocate(imstart(nb), isstart(nb))
      end if
      call mpi_bcast(nr,      nqp, mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast(nz,      nqp, mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast(nl,      nqp, mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast(ns,      nqp, mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast(npar,    nqp, mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast(imstart,  nb, mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast(isstart,  nb, mpi_integer, 0, mpi_comm_world, ierr)
      ! Physics
      call mpi_bcast(hfb_npr,    3,          mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast(hfb_lambda, 2, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(hfb_elast,  2, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(hfb_delta,  2, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(hfb_beta2,  1, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(hfb_energy, 1, mpi_double_precision, 0, mpi_comm_world, ierr)
      ! Pairing window
      call mpi_bcast(hfb_pairing_window, 1, mpi_double_precision, 0, mpi_comm_world, ierr)
      if (myrank /= 0) then
         allocate(pwi_start_n(nb), pwi_start_p(nb), pwi_dim_n(nb), pwi_dim_p(nb))
         allocate(pwi_qp_p(nqp), pwi_uv_p(nqp), pwi_qp_n(nqp), pwi_uv_n(nqp))
         allocate(pwi_active_p(nqp), pwi_active_n(nqp))
      end if
      call mpi_bcast(pwi_start_n,   nb, mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast(pwi_start_p,   nb, mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast(pwi_dim_n,     nb, mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast(pwi_dim_p,     nb, mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast(pwi_qp_n,     nqp, mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast(pwi_qp_p,     nqp, mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast(pwi_uv_n,     nqp, mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast(pwi_uv_p,     nqp, mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast(pwi_active_n, nqp, mpi_logical, 0, mpi_comm_world, ierr)
      call mpi_bcast(pwi_active_p, nqp, mpi_logical, 0, mpi_comm_world, ierr)
      ! Wave functions
      call mpi_bcast(nghl, 1, mpi_integer, 0, mpi_comm_world, ierr)
      if (myrank /= 0) then
         allocate(wf(nghl,nqp), wfdr(nghl,nqp), wfdz(nghl,nqp), wfd2(nghl,nqp))
         allocate(y(nghl), z(nghl), wdcori(nghl))
      end if
      call mpi_bcast(wf,     nghl*nqp, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(wfdr,   nghl*nqp, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(wfdz,   nghl*nqp, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(wfd2,   nghl*nqp, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(y,          nghl, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(z,          nghl, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(wdcori,     nghl, mpi_double_precision, 0, mpi_comm_world, ierr)
      ! Pairing
      call mpi_bcast(hfb_cpair,      2, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(hfb_alpha_pair, 2, mpi_double_precision, 0, mpi_comm_world, ierr)
      ! Coupling constants
      call mpi_bcast(hfb_cr0,         1, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(hfb_crr,         1, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(hfb_cdrho,       1, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(hfb_ctau,        1, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(hfb_ctj,         1, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(hfb_crdj,        1, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(hfb_use_j2terms, 1, mpi_logical,          0, mpi_comm_world, ierr)
      ! Finite-temperature
      call mpi_bcast(ft_active, 1, mpi_logical, 0, mpi_comm_world, ierr)
      if (ft_active) then
         call mpi_bcast(ft_temp,    1, mpi_double_precision, 0, mpi_comm_world, ierr)
         call mpi_bcast(ft_entropy, 3, mpi_double_precision, 0, mpi_comm_world, ierr)
         if (myrank > 0) allocate(ft_fp(nqp), ft_fn(nqp))
         call mpi_bcast(ft_fp,    nqp, mpi_double_precision, 0, mpi_comm_world, ierr)
         call mpi_bcast(ft_fn,    nqp, mpi_double_precision, 0, mpi_comm_world, ierr)
      end if
      !-------------------------------------------------------------------------
 
      ! Interaction
      !-------------------------------------------------------------------------
      if (myrank /= 0) allocate(crho(nghl), cs(nghl), cpair(nghl), cspair(nghl))
      call mpi_bcast(crho,   nghl, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(cs,     nghl, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(cpair,  nghl, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(cspair, nghl, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(cdrho,     1, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(ctau,      1, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(ctj0,      1, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(ctj1,      1, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(ctj2,      1, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(crdj,      1, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(cds,       1, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(ct,        1, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(cj,        1, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(cgs,       1, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(cf,        1, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(csdj,      1, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(interaction_name,         80, mpi_character, 0, mpi_comm_world, ierr)
      call mpi_bcast(skip_residual_interaction, 1,   mpi_logical, 0, mpi_comm_world, ierr)
      !-------------------------------------------------------------------------
 
   end subroutine distribute_pnfam_quantities
   
   !----------------------------------------------------------------------------
   ! How many OpenMP threads do we have? 
   !----------------------------------------------------------------------------
   function get_num_omp_threads() result(n)
      implicit none
      integer :: n
      n = 0
      !$ n = omp_get_max_threads()
   end function get_num_omp_threads
   
   !----------------------------------------------------------------------------
   ! Get the wallclock time with or without OpenMP 
   !----------------------------------------------------------------------------
   subroutine get_timer(t)
      implicit none
      real(dp), intent(out) :: t
      integer :: it, ir
      if (get_num_omp_threads() == 0) then
         call system_clock(count=it, count_rate=ir)
         t = 1.0_dp*it/ir
      else
         !$ t = omp_get_wtime()
      end if
   end subroutine get_timer
 
   !----------------------------------------------------------------------------
   ! Get the integers for the calculation 
   ! MPI things are global (myrank, mpi_comm_world)
   !----------------------------------------------------------------------------
   subroutine get_calculation_indices(dim, ncpu, is, ie)
      implicit none
      
      integer, intent(in)  :: dim, ncpu
      integer, intent(out) :: is, ie      
      integer :: chunk_size, num_leftover, my_count
      
      chunk_size = dim/ncpu
      num_leftover = mod(dim, ncpu)
      
      is = chunk_size*myrank + 1
      if ((myrank.ge.1) .and. (myrank.le.num_leftover)) then
         is = is + myrank
      else if ((myrank.ge.1) .and. (myrank.gt.num_leftover)) then
         is = is + num_leftover
      end if
      
      my_count = chunk_size
      if (myrank < num_leftover) then
         my_count = my_count + 1
      end if
      ie = is - 1 + my_count
 
   end subroutine get_calculation_indices
   
   !----------------------------------------------------------------------------
   ! Use MPI-IO to store a column vector to one or more files
   ! using the file logic from my note.
   !----------------------------------------------------------------------------
   subroutine mpi_store_ab_lower(vector, istart, lda, class, fh_list)
      use constants,    only : translate_uppercase
      use pnfam_matrix, only : io_max_rec

      use mpi

      implicit none
      
      real(dp),         intent(in) :: vector(:)
      integer,          intent(in) :: istart, lda
      integer,          intent(in) :: fh_list(:)
      character(len=1), intent(in) :: class
      
      !include 'mpif.h'
      
      integer          :: if1, if2, irec, i, istatus(mpi_status_size)
      integer(di)      :: offset
      character(len=1) :: c
      
      c = class
      call translate_uppercase(c)
      
      ! The B matrix is the second lower-triangular matrix,
      ! occuring after num_ab*(1+num_ab)/2 records.
      if (c == 'A') then
         offset = 0
      else if (c == 'B') then
         offset = 1_di*lda*(1_di+lda)/2
      else
         offset = 0  ! Satisfy -Wmaybe-uninitialized
         write(*,'(1x,2a)') 'ERROR: unknown class in store_ab_lower_parallel --- ', c
         call mpi_abort(mpi_comm_world, -1, ierr)
      end if
      
      ! If everything is on the same file, write in one go
      if1  = ceiling((istart + 0.5_dp*(2_di*lda-istart)*(istart-1) + offset)/io_max_rec, kind=si)
      if2  = ceiling((lda    + 0.5_dp*(2_di*lda-istart)*(istart-1) + offset)/io_max_rec, kind=si)
      irec = int(istart + (2_di*lda-istart)*(istart-1)/2 + offset - (if1-1)*io_max_rec, kind=si)
      
      if (if1 == if2) then
         call mpi_file_write_at(fh_list(if1), int(irec-1, kind=mpi_offset_kind), &
            vector, size(vector), mpi_double_precision, istatus, ierr)
         
         if (ierr /= 0) then
            write(*,'(1x,a,1i0,2a,2(a,1i0))') &
               'ERROR in MPI_FILE_WRITE_AT (AT-ONCE): ierr = ', ierr, &
               ', class = ', c, ', ifile = ', if1, ', irec = ', irec
         end if
      ! If the vector straddles two files, write entry-by-entry to get it right
      else
         do i=1, size(vector)
            if1  = ceiling((istart+i-1 + 0.5_dp*(2_di*lda-istart)*(istart-1) + offset)/io_max_rec, kind=si)
            irec = int(istart+i-1 + (2_di*lda-istart)*(istart-1)/2 + offset - (if1-1)*io_max_rec, kind=si)
            
            call mpi_file_write_at(fh_list(if1), int(irec-1, kind=mpi_offset_kind), &
               vector(i), 1, mpi_double_precision, istatus, ierr)
            
               if (ierr /= 0) then
                  write(*,'(1x,a,1i0,2a,2(a,1i0))') &
                     'ERROR in MPI_FILE_WRITE_AT (ONE-BY-ONE): ierr = ', ierr, &
                     ', class = ', c, ', ifile = ', if1, ', irec = ', irec
               end if
         end do
      end if
 
   end subroutine mpi_store_ab_lower

end program pnfam_matrix_ab
