!------------------------------------------------------------------------------
! pnfam_matrix/main_transitions_twoqp.f90
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
   use pnfam_matrix,      only : input_file_from_cli, make_two_qp_basis, log_time, &
                                 time_reversal_phases, qp_reverse, get_qp_indices, &
                                 find_blocked_level, remove_levels_from_basis,     &
                                 setup_all_operators, write_main_header, read_xye, &
                                 display_qp_level, display_1to1_energy
   implicit none

   integer, parameter :: dp = kind(1.0d0)

   ! Configuration
   type(matrix_config) :: config
   type(hfb_solution)  :: hfb_soln

   ! Two-quasiparticle basis (forward-basis only!)
   integer, allocatable :: twoqp(:,:), tphase(:), tmp_twoqp(:,:)
   integer :: num_diag, bo_index, bo_isospin, bo_orig_index

   ! X and Y amplitudes
   integer  :: i, ierr, recl_xy, nxy
   real(dp), allocatable :: eqrpa(:), x(:,:), y(:,:)

   ! Operators
   integer :: num_operators
   integer, allocatable :: f20_lookup(:,:), f02_lookup(:,:)
   type(blockmatrix), allocatable :: F20(:), F02(:)

   ! Matrix elements
   integer  :: iop, iev, ixy, iqp1, iqp2, iqpr1, iqpr2, ix, iy, tp
   real(dp), allocatable :: me(:,:), sr(:)

   ! Misc
   character(len=210) :: infile

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
      label='QRPA Two-Quasiparticle Matrix Elements')

   ! Basis
   call make_two_qp_basis(class='20', k=config%basis_k, parity=config%basis_parity, &
      ecut=config%basis_cutoff, hfb=hfb_soln, basis_out=twoqp)

   ! Option to remove a level from the QRPA for the odd nucleus
   if (config%odd_block_n(1) /= 0 .or. config%odd_block_p(1) /= 0) then
      if (config%odd_block_n(1) /= 0 .and. config%odd_block_p(1) /= 0) then
         write(*,'(a)') "ERROR: the program can only handle a single blocked level."
         stop
      end if

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

      call remove_levels_from_basis(in=twoqp, out=tmp_twoqp, it=bo_isospin, &
         levels=[bo_index, qp_reverse(bo_index)])

      ! Copy the reduced basis to twoqp(:,:)
      deallocate(twoqp)
      allocate(twoqp(2, size(tmp_twoqp, 2)))
      twoqp(:,:) = tmp_twoqp(:,:)
      deallocate(tmp_twoqp)

      write(*,'(a,1x,a)') 'ODD NUCLEON/BLOCKING'
      write(*,'(a)') repeat('-', 33)
      call display_qp_level(iqp=bo_index, it=bo_isospin, hfb=hfb_soln)
      call display_qp_level(iqp=qp_reverse(bo_index), it=bo_isospin, hfb=hfb_soln)
      call display_1to1_energy(iqp=bo_index, it=bo_isospin, hfb=hfb_soln, sign=sign(1, bo_orig_index))
   end if

   num_diag = size(twoqp,2)
   write(*,'(a,1i0)') 'Two-QP dimension ....... : ', num_diag
   call time_reversal_phases(basis=twoqp, phases=tphase)
   write(*,*)


   ! X and Y
   !----------------------------------------------------------------------------
   inquire(iolength=recl_xy) [(1.0_dp, i=1,1+2*num_diag)]

   call log_time('Reading in X and Y')
   call read_xye(file=config%file_xy, n_basis=num_diag, x_out=x, y_out=y, e_out=eqrpa)
   nxy = size(x, 2)


   ! Test X and Y
   if (config%me_check_norm) then
      call log_time('Checking orthonormality')
      call check_xy_orthonormal
   end if

   if (config%me_check_complete) then
      call log_time('Checking completeness')
      call check_xy_complete
   end if


   ! Operators in the q.p. basis
   call setup_all_operators(class='20', k=config%basis_k, parity=config%basis_parity, &
      hfb=hfb_soln, fqp_out=F20)
   call setup_all_operators(class='02', k=config%basis_k, parity=config%basis_parity, &
      hfb=hfb_soln, fqp_out=F02)

   num_operators = size(F20)

   call log_time('Calculating matrix elements')


   ! Transition matrix elements ...
   ! These look the same for both K = 0 and K != 0
   allocate(me(nxy, num_operators))
   me(:,:) = 0

   ! For sum rules
   allocate(sr(num_operators))
   sr = 0

   do iop=1, num_operators
      call get_qp_indices(F20(iop), f20_lookup)
      call get_qp_indices(F02(iop), f02_lookup)

      !$omp  parallel do default(shared) reduction(+:me) reduction(+:sr) &
      !$omp& private(iev,ixy,iqp1,iqp2,iqpr1,iqpr2,ix,iy,tp)
      do iev=1, nxy
         do ixy=1, num_diag
            iqp1  = twoqp(1,ixy)
            iqp2  = twoqp(2,ixy)
            iqpr1 = qp_reverse(iqp1)
            iqpr2 = qp_reverse(iqp2)

            ix = f20_lookup(iqp1, iqp2)
            iy = f02_lookup(iqpr1,iqpr2)

            tp = tphase(ixy)
            me(iev,iop) = me(iev,iop) + x(ixy,iev)*F20(iop)%elem(ix) + tp*y(ixy,iev)*F02(iop)%elem(iy)
         end do
         sr(iop) = sr(iop) + me(iev,iop)**2
      end do
      !omp end parallel do
   end do


   ! Storage in plain text
   call log_time('Storing the solution')

   open(1, file=trim(config%file_me_two), status='replace', iostat=ierr)

   if (ierr /= 0) then
      write(*,'(3a)') 'ERROR opening file ', trim(config%file_me_two), ' for writing.'
      stop
   end if

   do iev=1, nxy
      write(1,'(1f22.16)',advance='no') eqrpa(iev)

      do iop=1, num_operators
         write(1,'(1es28.16e3)',advance='no') me(iev,iop)
      end do

      write(1,'("")')
   end do

   close(1)


   ! Sum rules
   call log_time('Calculating sum rules')
   write(*,*)
   write(*,'(a)') "OPERATOR         Beta-         Beta+         Sum Rule      Expected      Delta         Percent"
   write(*,'(a)') repeat('-', 100)

   select case (config%basis_parity)

      ! Allowed operators
      case (1)
         select case (config%basis_k)
            case (0)
               call sumrule(k=config%basis_k, op="F",  label="Fermi .... :", hfb=hfb_soln, srm=sr(1), srp=sr(3))
               call sumrule(k=config%basis_k, op="GT", label="GT ....... :", hfb=hfb_soln, srm=sr(2), srp=sr(4))
            case (1, -1)
               call sumrule(k=config%basis_k, op="GT", label="GT ....... :", hfb=hfb_soln, srm=sr(1), srp=sr(2))
            case default
               continue
         end select

      ! Forbidden
      case (-1)
         select case (config%basis_k)
            case (0)
               call sumrule(k=config%basis_k, op="RS0", label="[RS]0 ... :", hfb=hfb_soln, srm=sr(1), srp=sr(7))
               call sumrule(k=config%basis_k, op="RS1", label="[RS]1 ... :", hfb=hfb_soln, srm=sr(3), srp=sr(9))
               call sumrule(k=config%basis_k, op="RS2", label="[RS]2 ... :", hfb=hfb_soln, srm=sr(6), srp=sr(12))
               call sumrule(k=config%basis_k, op="R",   label="R ....... :", hfb=hfb_soln, srm=sr(4), srp=sr(10))
               call sumrule(k=config%basis_k, op="P",   label="P ....... :", hfb=hfb_soln, srm=sr(5), srp=sr(11))
               call sumrule(k=config%basis_k, op="PS0", label="[PS]0 ... :", hfb=hfb_soln, srm=sr(2), srp=sr(8))
            case (1, -1)
               call sumrule(k=config%basis_k, op="RS1", label="[RS]1 ... :", hfb=hfb_soln, srm=sr(1), srp=sr(5))
               call sumrule(k=config%basis_k, op="RS2", label="[RS]2 ... :", hfb=hfb_soln, srm=sr(4), srp=sr(8))
               call sumrule(k=config%basis_k, op="R",   label="R ....... :", hfb=hfb_soln, srm=sr(2), srp=sr(6))
               call sumrule(k=config%basis_k, op="P",   label="P ....... :", hfb=hfb_soln, srm=sr(3), srp=sr(7))
            case (2, -2)
               call sumrule(k=config%basis_k, op="RS2", label="[RS]2 ... :", hfb=hfb_soln, srm=sr(1), srp=sr(2))
            case default
               continue
         end select

      case default
         continue
   end select

contains

   !----------------------------------------------------------------------------
   ! Calculate e.g. the Ikeda sum rule
   !----------------------------------------------------------------------------
   subroutine sumrule(k, op, label, hfb, srm, srp)
      use constants,         only : translate_uppercase, IT_NEUTRON, IT_PROTON, IT_ISOSCALAR
      use hfbtho_basis,      only : nghl, wdcori, y, z
      use hfb_solution_type, only : hfb_solution
      implicit none

      integer,            intent(in) :: k
      character(len=*),   intent(in) :: op, label
      type(hfb_solution), intent(in) :: hfb
      real(dp),           intent(in) :: srm, srp

      character(len=len(op)) :: op_
      character(len=200) :: fmt
      real(dp) :: expected, rp(nghl), rn(nghl), r(nghl), w(nghl)

      op_ = op
      call translate_uppercase(op_)

      r(:) = 1.0_dp/y(:)
      w(:) = 1.0_dp/wdcori(:)

      ! Calculate HFB expectation value
      select case (op_)
         case ("F")
            call rz_density(hfb=hfb, rho=rn, it=IT_NEUTRON)
            call rz_density(hfb=hfb, rho=rp, it=IT_PROTON)
            expected = sum(w*(rn-rp))

         case ("GT")
            call rz_density(hfb=hfb, rho=rn, it=IT_NEUTRON)
            call rz_density(hfb=hfb, rho=rp, it=IT_PROTON)
            expected = sum(w*(rn-rp))

         case ("R")
            call rz_density(hfb=hfb, rho=rn, it=IT_NEUTRON)
            call rz_density(hfb=hfb, rho=rp, it=IT_PROTON)
            select case (k)
               case (0)
                  expected = sum(w*(rn-rp) * z**2)
               case (1, -1)
                  expected = 0.5_dp*sum(w*(rn-rp) * r**2)
               end select

         case ("P")
            select case (k)
               case (0)
                  call rz_density(hfb=hfb, rho=rn, it=IT_NEUTRON, dz=.true.)
                  call rz_density(hfb=hfb, rho=rp, it=IT_PROTON,  dz=.true.)
                  expected = sum(w*(rn-rp))
               case (1, -1)
                  call rz_density(hfb=hfb, rho=rn, it=IT_NEUTRON, dr=.true.)
                  call rz_density(hfb=hfb, rho=rp, it=IT_PROTON,  dr=.true.)
                  expected = 0.5_dp*sum(w*(rn-rp))
            end select

         case ("RS0")
            call rz_density(hfb=hfb, rho=rn, it=IT_NEUTRON)
            call rz_density(hfb=hfb, rho=rp, it=IT_PROTON)
            expected = sum(w*(rn-rp) * (z**2 + r**2))

         case ("PS0")
            call rz_density(hfb=hfb, rho=rn, it=IT_NEUTRON, dz=.true., dr=.true.)
            call rz_density(hfb=hfb, rho=rp, it=IT_PROTON,  dz=.true., dr=.true.)
            expected = sum(w*(rn-rp))

         case ("RS1")
            call rz_density(hfb=hfb, rho=rn, it=IT_NEUTRON)
            call rz_density(hfb=hfb, rho=rp, it=IT_PROTON)
            select case (k)
               case (1, -1)
                  expected = 3.0_dp/4.0_dp*sum(w*(rn-rp) * (2*z**2 + r**2))
               case (0)
                  expected = 3.0_dp/2.0_dp*sum(w*(rn-rp) * r**2)
            end select

         case ("RS2")
            call rz_density(hfb=hfb, rho=rn, it=IT_NEUTRON)
            call rz_density(hfb=hfb, rho=rp, it=IT_PROTON)
            select case (k)
               case (0)
                  expected = 1.0_dp/2.0_dp*sum(w*(rn-rp) * (4*z**2 + r**2))
               case (1, -1)
                  expected = 3.0_dp/4.0_dp*sum(w*(rn-rp) * (2*z**2 + r**2))
               case (2, -2)
                  expected = 3.0_dp/2.0_dp*sum(w*(rn-rp) * r**2)
            end select
      end select

      fmt = '(a,5f14.6,1f10.2,"%")'
      write(*,fmt) trim(label), srm, srp, srm-srp, expected, (srm-srp)-expected, 100*((srm-srp)-expected)/expected

   end subroutine sumrule

   !----------------------------------------------------------------------------
   ! Calculate the HFB coordinate-space density, allowing for
   ! selecting states of a certain spin or isospin.
   !
   ! Operates by calculating VV', then selecting only states with matching spin.
   ! See matrix code appendix.
   !----------------------------------------------------------------------------
   subroutine rz_density(hfb, rho, it, dz, dr)
      use constants, only : IT_PROTON, IT_NEUTRON
      use blockmatrix_type,  only : nb, db, blockmatrix, triprod, allocate_blockmatrix, copy_block_structure
      use hfbtho_basis,      only : nghl, isstart, imstart, wf, wdcori, ns, wfdz, y, nl, wfdr
      use hfb_solution_type, only : hfb_solution
      implicit none

      type(hfb_solution), intent(in)  :: hfb
      real(dp),           intent(out) :: rho(nghl)
      integer,            intent(in)  :: it
      logical, optional,  intent(in)  :: dz, dr

      integer  :: ib, ir, ic, i1, i2, ik
      type(blockmatrix) :: rhoab, one

      ! Define the unit HFB matrix
      call allocate_blockmatrix(one, size(hfb%un%elem))
      call copy_block_structure(hfb%un, one)
      one%elem(:) = 0
      do ib=1, nb
         do ir=1, db(ib)
            one%elem(imstart(ib)-1 + ir + (ir-1)*db(ib)) = 1.0_dp
         end do
      end do

      ! Calculate the neutron or proton density in configuration space
      call allocate_blockmatrix(rhoab, size(hfb%un%elem))
      if (it == IT_NEUTRON) then
         call triprod('n', hfb%vn, 'n', one, 't', hfb%vn, 1d0, 0d0, rhoab)
      else if (it == IT_PROTON) then
         call triprod('n', hfb%vp, 'n', one, 't', hfb%vp, 1d0, 0d0, rhoab)
      end if

      ! Make the coordinate-space density
      rho(:) = 0
      !$omp  parallel do default(shared) private(ib,ic,ir,i1,i2,ik) reduction(+:rho)
      do ib=1, nb
         do ic=1, db(ib)
            i2 = isstart(ib)-1 + ic
            do ir=1, db(ib)
               i1 = isstart(ib)-1 + ir
               ik = imstart(ib)-1 + ir + (ic-1)*db(ib)
               ! Spin selection
               if (ns(i1) == ns(i2)) then
                  if ((.not.present(dz)) .and. (.not.present(dr))) then
                     rho(:) = rho(:) + rhoab%elem(ik)*wf(:,i1)*wf(:,i2)
                  else
                     if (present(dz) .and. (dz)) then
                        rho(:) = rho(:) + rhoab%elem(ik)*wfdz(:,i1)*wfdz(:,i2)
                     end if
                     if (present(dr) .and. (dr)) then
                        rho(:) = rho(:) + rhoab%elem(ik)*(wfdr(:,i1)*wfdr(:,i2) &
                           - wfdr(:,i1)*wf(:,i2)*y(:) + wf(:,i1)*wf(:,i2)*nl(i1)*nl(i2)*y(:)**2)
                     end if
                  end if
               end if
            end do
         end do
      end do
      !$omp end parallel do

      rho(:) = rho(:)*wdcori(:)

   end subroutine rz_density

   !----------------------------------------------------------------------------
   ! Check that X and Y are orthonormal
   !
   ! x and y are globals
   !----------------------------------------------------------------------------
   subroutine check_xy_orthonormal
      implicit none

      real(dp), parameter :: tol = 5.0d-13

      integer  :: i1, i2, nxy, nev
      real(dp) :: check

      nxy = size(x,1)
      nev = size(x,2)

      do i1=1, nev
         do i2=i1, nev
            check = dot_product(x(:,i1),x(:,i2))-dot_product(y(:,i1),y(:,i2))

            if ((i1 == i2) .and. (abs(check-1.0_dp) > tol)) then
               write(*,'(a,2i8,1es12.1)') 'Warning: X,Y may not be orthonormal!', i1, i2, check
            else if ((i1 /= i2) .and. (abs(check) > tol)) then
               write(*,'(a,2i8,1es12.1)') 'Warning: X,Y may not be orthonormal!', i1, i2, check
            end if
         end do
      end do

   end subroutine check_xy_orthonormal

   !----------------------------------------------------------------------------
   ! Check that X and Y are complete
   !
   ! x and y are globals
   !----------------------------------------------------------------------------
   subroutine check_xy_complete
      implicit none

      real(dp), parameter :: tol = 5.0d-13

      integer  :: i1, i2, nxy, nev
      real(dp) :: check

      nxy = size(x,1)
      nev = size(x,2)

      do i1=1, nxy
         do i2=i1, nxy
            check = dot_product(x(i1,:),x(i2,:))-dot_product(y(i1,:),y(i2,:))

            if ((i1 == i2) .and. (abs(check-1.0_dp) > tol)) then
               write(*,'(a,2i8,1es12.1)') 'Warning: X,Y may not be complete!', i1, i2, check
            else if ((i1 /= i2) .and. (abs(check) > tol)) then
               write(*,'(a,2i8,1es12.1)') 'Warning: X,Y may not be complete!', i1, i2, check
            end if
         end do
      end do

   end subroutine check_xy_complete

end program pnfam_matrix_transitions
