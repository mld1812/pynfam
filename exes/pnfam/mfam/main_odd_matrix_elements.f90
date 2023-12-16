!------------------------------------------------------------------------------
! pnfam_matrix/main_odd_matrix_elements.f90
!
! T. Shafer, UNC Chapel Hill, 2015
!------------------------------------------------------------------------------
program pnfam_matrix_transitions
   !$ use omp_lib,          only : omp_get_max_threads, omp_get_wtime
   use logger
   use constants,           only : IT_PROTON, IT_NEUTRON
   use hfbtho_basis,        only : get_hfbtho_solution
   use interaction,         only : init_skyrme
   use blockmatrix_type,    only : blockmatrix
   use external_field_type, only : external_field
   use hfb_solution_type,   only : hfb_solution
   use config_type,         only : matrix_config, read_namelist_from_file
   use pnfam_matrix,        only : input_file_from_cli, write_main_header, find_blocked_level,  &
                                   make_two_qp_basis, init_sp_operator, init_qp_operator,       &
                                   get_qp_indices, find_1to1_transitions, residual_interaction, &
                                   log_time, display_qp_level
   implicit none

   integer, parameter :: dp = kind(1.0d0)

   ! Configuration
   type(matrix_config) :: config
   type(hfb_solution)  :: hfb_soln

   ! Two-quasiparticle basis and blocking
   integer :: bo_index, bo_isospin, num_20
   integer, allocatable :: twoqp_20(:,:)

   ! Reference operators
   type(external_field) :: ref_f_sp
   type(blockmatrix)    :: ref_f20, ref_f02
   integer, allocatable :: f20_lookup(:,:)

   ! 1-to-1 settings
   type(blockmatrix)    :: ref_delta_r
   character(len=1)     :: class_1to1
   integer, allocatable :: twoqp_1to1(:,:), lookup_1to1(:,:)
   integer :: num_1to1, iqpn, iqpp, ik

   ! Calculations
   integer :: icalc, recl, irec
   real(dp), allocatable :: h20(:)


   ! Misc
   integer :: ierr
   character(len=200) :: infile, ls

   !----------------------------------------------------------------------------
   !----------------------------------------------------------------------------

   ! Logging
   openlog  => openlog_stdout
   writelog => writelog_stdout
   closelog => closelog_stdout
   call openlog('')

   ! Configuration and setup
   call input_file_from_cli(infile, default='matrix.in')

   open(1, file=trim(infile), status='old', iostat=ierr)
   if (ierr /= 0) stop 'Error opening Namelist file'

   call read_namelist_from_file(fh=1, config=config)

   call get_hfbtho_solution(fn=trim(config%file_hfb), &
      ep=hfb_soln%ep, vp=hfb_soln%vp, up=hfb_soln%up, &
      en=hfb_soln%en, vn=hfb_soln%vn, un=hfb_soln%un)

   call init_skyrme(fh=1,                       &
      vp=hfb_soln%vp%elem, vn=hfb_soln%vn%elem, &
      up=hfb_soln%up%elem, un=hfb_soln%un%elem)

   close(1)

   ! Header, etc.
   call write_main_header(config, infile=infile, nmpi=0, nomp=get_num_omp_threads(), &
      label='QRPA Odd-Nucleus Driving-Force Matrix Elements')

   ! Identify the "blocked" level to be used for the 1-to-1 transitions
   if (config%odd_block_n(1) == 0 .and. config%odd_block_p(1) == 0) then
      write(*,'(a)') 'ERROR: no odd orbital was specified in the namelist.'
      stop
   else if (config%odd_block_n(1) /= 0 .and. config%odd_block_p(1) /= 0) then
      write(*,'(a)') "ERROR: the program can only handle a single blocked level."
      stop
   else
      if (config%odd_block_n(1) /= 0) then
         bo_isospin = IT_NEUTRON
         bo_index   = find_blocked_level(bk=config%odd_block_n(1), &
            bparity=config%odd_block_n(2), bt=IT_NEUTRON, hfb=hfb_soln)
      else if (config%odd_block_p(1) /= 0) then
         bo_isospin = IT_PROTON
         bo_index   = find_blocked_level(bk=config%odd_block_p(1), &
            bparity=config%odd_block_p(2), bt=IT_PROTON,  hfb=hfb_soln)
      end if
   end if

   write(*,'(a,1x,a)') 'ODD NUCLEON/BLOCKING'
   write(*,'(a)') repeat('-', 33)
   call display_qp_level(iqp=bo_index, it=bo_isospin, hfb=hfb_soln)
   write(*,*)

   ! Two-QP basis for H(20)
   call make_two_qp_basis(class='20', k=config%basis_k, parity=config%basis_parity, &
      ecut=config%basis_cutoff, hfb=hfb_soln, basis_out=twoqp_20)

   num_20 = size(twoqp_20, 2)

   ! Reference operators
   call init_sp_operator(k=config%basis_k, parity=config%basis_parity, f_out=ref_f_sp)
   call init_qp_operator(class='20', f_sp=ref_f_sp, hfb=hfb_soln, f_qp_out=ref_f20)
   call init_qp_operator(class='02', f_sp=ref_f_sp, hfb=hfb_soln, f_qp_out=ref_f02)

   ! Lookup arrays
   call get_qp_indices(ref_f20, f20_lookup)

   write(*,'(a)')     'CALCULATION DETAILS'
   write(*,'(a)')     repeat('-',33)
   write(*,'(a,1i0)') 'Two-QP basis size (AB) ... : ', num_20

   inquire(iolength=recl) 1.0_dp
   write(*,'(a,1i0)') 'I/O record length ........ : ', recl

   ! 1-to-1 transitions
   if (bo_isospin == IT_NEUTRON) then
      class_1to1 = 'P'
      call find_1to1_transitions(class='11', ib=bo_index, tb=bo_isospin, k=config%basis_k, &
         parity=config%basis_parity, ecut=config%basis_cutoff, hfb=hfb_soln, basis_out=twoqp_1to1)
      call init_qp_operator(class='11', f_sp=ref_f_sp, hfb=hfb_soln, f_qp_out=ref_delta_r)
   else if (bo_isospin == IT_PROTON) then
      class_1to1 = 'Q'
      call find_1to1_transitions(class='11~', ib=bo_index, tb=bo_isospin, k=config%basis_k, &
         parity=config%basis_parity, ecut=config%basis_cutoff, hfb=hfb_soln, basis_out=twoqp_1to1)
      call init_qp_operator(class='11~', f_sp=ref_f_sp, hfb=hfb_soln, f_qp_out=ref_delta_r)
   else
      stop 'Unrecognized bo_isospin'
   end if

   num_1to1 = size(twoqp_1to1, 2)
   call get_qp_indices(ref_delta_r, lookup_1to1)

   write(*,'(2a)')    '1-to-1 matrix class ...... : ', class_1to1
   write(*,'(a,1i0)') 'No. 1-to-1 transitions ... : ', num_1to1
   write(*,'(a)') ''

   !----------------------------------------------------------------------------
   !----------------------------------------------------------------------------

   write(ls,'(3a)') 'Calculating matrix elements and writing to file "', trim(config%file_vv), '"'
   call log_time(ls)

   ! Main calculations
   open(1, file=trim(config%file_vv), recl=recl, status='unknown', &
      form='unformatted', access='direct', iostat=ierr)

   irec = 0

   do icalc=1, num_1to1
      iqpn = twoqp_1to1(1, icalc)
      iqpp = twoqp_1to1(2, icalc)
      ik   = lookup_1to1(iqpn, iqpp)

      write(ls,'("Transition ",1i0,": Ep = ",1f0.1," MeV, En = ",1f0.1," MeV")') &
         icalc, hfb_soln%ep(iqpp), hfb_soln%en(iqpn)
      call log_time(ls)

      call residual_interaction(class=class_1to1, hfb=hfb_soln, ik=ik, &
         delta_r_qp=ref_delta_r, f_sp=ref_f_sp%mat, f20_qp=ref_f20, f02_qp=ref_f02, &
         h20_vector=h20)

      do ik=1, num_20
         iqpn = twoqp_20(1, ik)
         iqpp = twoqp_20(2, ik)
         irec = irec+1
         write(1, rec=irec) h20(f20_lookup(iqpn,iqpp))
      end do
   end do

   close(1)

contains

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

end program pnfam_matrix_transitions
