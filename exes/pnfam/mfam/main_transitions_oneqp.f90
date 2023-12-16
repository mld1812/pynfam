!------------------------------------------------------------------------------
! pnfam_matrix/main_transitions_oneqp.f90
!
! Single-quasiparticle matrix elements, eventually corrected via the QRPA
! for the polarization due to the odd nucleon.
!
! T. Shafer, UNC Chapel Hill, 2015
!------------------------------------------------------------------------------
program pnfam_matrix_transitions
   use logger
   use constants,         only : IT_PROTON, IT_NEUTRON
   use blockmatrix_type,  only : blockmatrix
   use hfbtho_basis,      only : get_hfbtho_solution, lowest_energy_qp
   use hfb_solution_type, only : hfb_solution
   use config_type,       only : matrix_config, read_namelist_from_file
   use pnfam_matrix,      only : qp_reverse, input_file_from_cli, time_reversal_phases,           &
                                 get_qp_indices, log_time, find_blocked_level, make_two_qp_basis, &
                                 setup_all_operators, find_1to1_transitions, write_main_header,   &
                                 remove_levels_from_basis, read_xye, display_qp_level,            &
                                 display_1to1_energy
   implicit none

   integer, parameter :: dp = kind(1.0d0)

   ! Configuration
   type(matrix_config) :: config
   type(hfb_solution)  :: hfb_soln

   ! Two-quasiparticle basis and blocking
   integer, allocatable :: twoqp_odd(:,:), twoqp_ab(:,:), tphase_ab(:)
   integer :: bo_index, bo_isospin, bo_orig_index, num_ab

   ! Operators
   integer :: n_operators
   integer, allocatable :: f11_lookup(:,:), f11t_lookup(:,:), f20_lookup(:,:), f02_lookup(:,:)
   type(blockmatrix), allocatable :: F11(:), F11t(:), F20(:), F02(:)

   ! Transitions
   integer  :: iop, nme, iqpn, iqpp, im, ixy, nxy, ik, iqp1, iqp2, ix, iy, tp
   real(dp) :: factor1, factor2, factor3, factor4
   real(dp), allocatable :: energy(:), me(:,:), eqrpa(:), x(:,:), y(:,:), odd_me(:,:)

   ! Misc
   integer :: i, ierr
   character(len=210) :: infile, ls

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

   close(1)

   ! Header, etc.
   call write_main_header(config, infile=infile, nmpi=0, nomp=0, &
      label='QRPA One-Quasiparticle Matrix Elements (Odd Nucleus)')

   ! Identify the "blocked" level to be used for the 1-to-1 transitions
   if (config%odd_block_n(1) == 0 .and. config%odd_block_p(1) == 0) then
      write(*,'(a)') 'ERROR: no odd orbital was specified in the namelist.'
      stop
   else if (config%odd_block_n(1) /= 0 .and. config%odd_block_p(1) /= 0) then
      write(*,'(a)') "ERROR: the program can only handle a single blocked level."
      stop
   else
      if (config%odd_block_n(1) /= 0) then
         bo_isospin    = IT_NEUTRON
         bo_orig_index = config%odd_block_n(1)
         bo_index      = find_blocked_level(bk=config%odd_block_n(1), &
            bparity=config%odd_block_n(2), bt=IT_NEUTRON, hfb=hfb_soln)
      else if (config%odd_block_p(1) /= 0) then
         bo_isospin    = IT_PROTON
         bo_orig_index = config%odd_block_p(1)
         bo_index      = find_blocked_level(bk=config%odd_block_p(1), &
            bparity=config%odd_block_p(2), bt=IT_PROTON,  hfb=hfb_soln)
      end if
   end if

   write(*,'(a,1x,a)') 'ODD NUCLEON/BLOCKING'
   write(*,'(a)') repeat('-', 33)
   call display_qp_level(iqp=bo_index, it=bo_isospin, hfb=hfb_soln)
   call display_1to1_energy(iqp=bo_index, it=bo_isospin, hfb=hfb_soln, sign=sign(1, bo_orig_index))
   write(*,*)

   ! Two-QP basis used in the diagonalization
   call make_two_qp_basis(class='20', k=config%basis_k, parity=config%basis_parity, &
      ecut=config%basis_cutoff, hfb=hfb_soln, basis_out=twoqp_ab)
   call time_reversal_phases(twoqp_ab, tphase_ab)
   num_ab = size(twoqp_ab, 2)

   write(ls,'("NUM_AB = ",1i0)') num_ab
   call log_time(ls)


   ! Make sure we have been given the even-even XY amplitudes
   if (len_trim(config%odd_even_basename) == 0) then
      write(*,'(a)') "Error: no even-even XY file was specified; required for the correction"
      stop
   end if

   ! Operators
   call log_time('Building operators')

   call setup_all_operators(class='11',  k=config%basis_k, parity=config%basis_parity,&
      hfb=hfb_soln, fqp_out=F11)
   call setup_all_operators(class='11~', k=config%basis_k, parity=config%basis_parity,&
      hfb=hfb_soln, fqp_out=F11t)
   call setup_all_operators(class='20',  k=config%basis_k, parity=config%basis_parity,&
      hfb=hfb_soln, fqp_out=F20)
   call setup_all_operators(class='02',  k=config%basis_k, parity=config%basis_parity,&
      hfb=hfb_soln, fqp_out=F02)

   n_operators = size(F11)

   if (bo_isospin == IT_NEUTRON) then
      call find_1to1_transitions(class='11',  k=config%basis_k, parity=config%basis_parity, &
         ecut=config%basis_cutoff, hfb=hfb_soln, ib=bo_index, tb=bo_isospin, basis_out=twoqp_odd)
   else if (bo_isospin == IT_PROTON) then
      call find_1to1_transitions(class='11~', k=config%basis_k, parity=config%basis_parity, &
         ecut=config%basis_cutoff, hfb=hfb_soln, ib=bo_index, tb=bo_isospin, basis_out=twoqp_odd)
   end if

   nme = size(twoqp_odd, 2)

   write(ls,'("NME = ",1i0)') nme
   call log_time(ls)

   ! Read in X and Y (from the EVEN nucleus... no blocked levels)
   if (.not.config%me_disable_correction) then
      call log_time('Loading X and Y vectors')
      call read_xye(file=config%odd_even_file_xy, n_basis=num_ab, x_out=x, y_out=y, e_out=eqrpa)

      nxy = size(x, 2)

      write(ls,'("NXY = ",1i0)') nxy
      call log_time(ls)

      if (size(x,1) /= num_ab) then
         write(*,'(/,2(a,1i0),a)') "ERROR: XY QP dimension does not match NUM_AB (", &
            size(x,1), " != ", num_ab ,")"
         stop
      end if

      ! Read in the matrix elements of the odd quasiparticle
      call log_time('Loading odd quasiparticle matrix elements')
      call load_odd_me(file=config%file_vv, n_basis=num_ab, me_out=odd_me)

      write(ls,'("NVV = ",1i0)') size(odd_me, 2)
      call log_time(ls)

      if (size(odd_me, 2) < 1) then
         write(*,'(/,a)') 'ERROR: could not load polarization matrix elements.'
         stop
      end if
   end if


   ! Calculate matrix elements
   if (.not.config%me_disable_correction) then
      call log_time('Calculating matrix elements with polarization correction')
   else
      call log_time('Calculating matrix elements at 0th order (HFB matrix elements)')
   end if

   allocate(energy(nme), me(nme, n_operators))
   energy = 0;  me = 0

   ! Matrix elements
   do iop=1, n_operators
      call get_qp_indices(F11(iop),  f11_lookup)
      call get_qp_indices(F11t(iop), f11t_lookup)
      call get_qp_indices(F20(iop),  f20_lookup)
      call get_qp_indices(F02(iop),  f02_lookup)

      do im=1, nme
         iqpn = twoqp_odd(1,im)
         iqpp = twoqp_odd(2,im)

         ! <p| F |n> with beta-
         if (bo_isospin == IT_NEUTRON .and. iop <= n_operators/2) then
            me(im, iop) = F11(iop)%elem(f11_lookup(iqpn,iqpp))
         ! <p| F |n> with beta+
         else if (bo_isospin == IT_NEUTRON .and. iop  > n_operators/2) then
            me(im, iop) = F11(iop)%elem(f11_lookup(iqpn,iqpp))
         ! <n| F |p> with beta+
         else if (bo_isospin == IT_PROTON  .and. iop  > n_operators/2) then
            me(im, iop) = F11t(iop)%elem(f11t_lookup(iqpn,iqpp))
         ! <n| F |p> with beta-
         else if (bo_isospin == IT_PROTON  .and. iop <= n_operators/2) then
            me(im, iop) = F11t(iop)%elem(f11t_lookup(iqpn,iqpp))
         end if

         ! Energies
         if (iop == 1) then
            if (bo_isospin == IT_PROTON) then
               energy(im) = hfb_soln%en(iqpn)-hfb_soln%ep(iqpp)
            else if (bo_isospin == IT_NEUTRON) then
               energy(im) = hfb_soln%ep(iqpp)-hfb_soln%en(iqpn)
            end if
         end if

         ! Correction
         if (.not.config%me_disable_correction) then
            !$omp  parallel do default(shared) reduction(+:me) &
            !$omp& private(ixy,factor1,factor2,factor3,factor4,ik,iqp1,iqp2,ix,iy,tp)
            do ixy=1, nxy

               ! "Forward" correction
               factor1 = 0;  factor2 = 0

               do ik=1, num_ab
                  iqp1 = twoqp_ab(1, ik)
                  iqp2 = twoqp_ab(2, ik)
                  ix   = f20_lookup(iqp1, iqp2)
                  iy   = f02_lookup(qp_reverse(iqp1), qp_reverse(iqp2))
                  tp   = tphase_ab(ik)

                  factor1 = factor1 + x(ik,ixy)*f20(iop)%elem(ix) + tp*y(ik,ixy)*f02(iop)%elem(iy)
                  factor2 = factor2 + x(ik,ixy)*odd_me(ik,im)     + y(ik,ixy)*odd_me(ik,im)
               end do

               ! "Reverse" correction
               factor3 = 0;  factor4 = 0

               do ik=1, num_ab
                  iqp1 = twoqp_ab(1, ik)
                  iqp2 = twoqp_ab(2, ik)

                  if (config%basis_k == 0) then
                     ix   = f20_lookup(qp_reverse(iqp1), qp_reverse(iqp2))
                     iy   = f02_lookup(iqp1, iqp2)
                     tp   = tphase_ab(ik)

                     factor3 = factor3 + tp*y(ik,ixy)*f20(iop)%elem(ix) + x(ik,ixy)*f02(iop)%elem(iy)
                     factor4 = factor4 +    y(ik,ixy)*odd_me(ik,im)     + x(ik,ixy)*odd_me(ik,im)
                  else
                     ix   = f20_lookup(iqp1, iqp2)
                     iy   = f02_lookup(qp_reverse(iqp1), qp_reverse(iqp2))
                     tp   = tphase_ab(ik)

                     factor3 = factor3 + y(ik,ixy)*f20(iop)%elem(ix) + tp*x(ik,ixy)*f02(iop)%elem(iy)
                     factor4 = factor4 + y(ik,ixy)*odd_me(ik,im)     +    x(ik,ixy)*odd_me(ik,im)
                  end if
               end do

               me(im,iop) = me(im,iop) + factor1*factor2/(energy(im)-eqrpa(ixy)) &
                                       - factor3*factor4/(energy(im)+eqrpa(ixy))
            end do
            !$omp end parallel do
         end if
      end do
   end do


   ! Storage in plain text
   call log_time('Storing results')

   open(1, file=trim(config%file_me_one), status='replace', iostat=ierr)
   if (ierr /= 0) then
      write(*,'(3a)') 'ERROR opening file ', trim(config%file_me_one), ' for writing.'
      stop
   end if

   do i=1, nme
      write(1,'(1f22.16)',advance='no') energy(i)
      do iop=1, n_operators
         write(1,'(1es28.16e3)',advance='no') me(i,iop)
      end do
      write(1,'("")')
   end do
   close(1)

contains

   !----------------------------------------------------------------------------
   ! Load the odd-nucleon matrix elements
   !----------------------------------------------------------------------------
   subroutine load_odd_me(file, n_basis, me_out)
      implicit none

      character(len=*),      intent(in)  :: file
      integer,               intent(in)  :: n_basis
      real(dp), allocatable, intent(out) :: me_out(:,:)

      integer  :: recl, ierr, i, ime, nme, irec


      if (allocated(me_out)) deallocate(me_out)

      inquire(iolength=recl) 1.0_dp

      open(11, file=trim(file), recl=recl, status='old', &
         form='unformatted', access='direct', iostat=ierr)

      ! On an error, return a zero-length array
      if (ierr /= 0) then
         allocate(me_out(0,0))
         return
      end if

      ! Count records
      i = 0
      do
         read(11, rec=i+1, iostat=ierr)

         if (ierr /= 0) then
            exit
         else
            i = i+1
         end if
      end do

      ! The number of records should be a multiple of n_basis
      if (mod(i, n_basis) /= 0) then
         write(*,'(a,1i0)') 'ERROR in load_odd_me: n_records != [int] x n_basis. MOD = ', mod(i, n_basis)
         stop
      end if

      ! The number of transitions is just (# records)/(# basis)
      nme = i/n_basis
      allocate(me_out(n_basis, nme))
      me_out(:,:) = 0

      ! Store the vector
      irec = 0
      do ime=1, nme
         do i=1, n_basis
            irec = irec+1
            read(11, rec=irec) me_out(i,ime)
         end do
      end do

      close(11)

   end subroutine load_odd_me


end program pnfam_matrix_transitions
