!------------------------------------------------------------------------------
! pnfam_matrix/mod_common.f90
!
! T. Shafer, UNC Chapel Hill, 2015
!------------------------------------------------------------------------------
module pnfam_matrix

   implicit none
   public

   integer, parameter, private :: di = selected_int_kind(11)
   integer, parameter, private :: dp = kind(1.0d0)

   ! Biggest IO record number
   integer(di), parameter :: io_max_rec = 2**30-1

   ! Starting file unit for AB matrices (io_ab_unit, io_ab_unit+1, ...)
   integer, parameter :: io_ab_unit = 60

contains

   !----------------------------------------------------------------------------
   ! Take a QP index and change it from -K to +K or vice-versa
   !----------------------------------------------------------------------------
   function qp_reverse(iqp)
      use blockmatrix_type, only : db
      integer, intent(in) :: iqp
      integer :: qp_reverse

      if (iqp > sum(db)/2) then
         qp_reverse = iqp - sum(db)/2
      else
         qp_reverse = iqp + sum(db)/2
      end if
   end function qp_reverse

   !----------------------------------------------------------------------------
   ! Get QP indices of an arbitrary FAM matrix
   ! Returns array(dim(n),dim(p)) so array(iqpn,iqpp) = i_me.
   ! Or, in the more general case, returns array(icol,irow) = i_me.
   !----------------------------------------------------------------------------
   subroutine get_qp_indices(mat, array)
      use hfbtho_basis,     only : isstart
      use blockmatrix_type, only : blockmatrix, nb, db

      implicit none

      type(blockmatrix),    intent(in)  :: mat
      integer, allocatable, intent(out) :: array(:,:)
      integer :: ibr, ibc, ic, ir

      if (allocated(array)) deallocate(array)
      allocate(array(sum(db), sum(db)))
      array(:,:) = 0

      ! Block row
      do ibr=1, nb
         ! Block column
         ibc = mat%ir2c(ibr)

         if (ibc == 0) cycle

         ! Row-major within blocks
         do ic=1, db(ibc)
            do ir=1, db(ibr)
               array(isstart(ibc)+ic-1,isstart(ibr)+ir-1) = mat%ir2m(ibr) + (ic-1)*db(ibr) + ir - 1
            end do
         end do
      end do

   end subroutine get_qp_indices

   !----------------------------------------------------------------------------
   ! Read input file name from the command line, falling back to 'default'
   !----------------------------------------------------------------------------
   subroutine input_file_from_cli(file, default)
      implicit none

      character(len=*), intent(out) :: file
      character(len=*), intent(in)  :: default
      integer :: nargs

      nargs = command_argument_count()
      if (nargs == 0) then
         file = adjustl(default)
      else
         if (nargs > 1) then
            stop 'This program only allows one input file per run.'
         else
            call get_command_argument(1, file)
         end if
      end if

   end subroutine input_file_from_cli

   !----------------------------------------------------------------------------
   ! Determine the action of T on each 2-qp state |p n>.
   ! E.g. T|p n> = |p- n->, but T|p n-> = -|p- n>, etc. due to the double
   ! action of T-reversal.
   !
   ! 'phases' is of dimension 2*num_ab, including the K<0 basis.
   !
   ! This depends on the HFB solution already being loaded via hfbtho_basis!
   !----------------------------------------------------------------------------
   subroutine time_reversal_phases(basis, phases)
      use hfbtho_basis, only : nl, ns
      implicit none

      integer,              intent(in)  :: basis(:,:)
      integer, allocatable, intent(out) :: phases(:)

      integer :: i, iqpp, iqpn, iqrp, iqrn, dim

      dim = size(basis, 2)

      if (allocated(phases)) deallocate(phases)
      allocate(phases(2*dim))

      phases(:) = 1

      do i=1, dim
         iqpn = basis(1,i)
         iqpp = basis(2,i)
         iqrn = qp_reverse(basis(1,i))
         iqrp = qp_reverse(basis(2,i))

         ! This works since the quasiparticles preserve K (though not L or S)
         phases(i)     = merge(-1, 1, (2*nl(iqpn) + ns(iqpn)) < 0) * merge(-1, 1, (2*nl(iqpp) + ns(iqpp)) < 0)
         phases(i+dim) = merge(-1, 1, (2*nl(iqrn) + ns(iqrn)) < 0) * merge(-1, 1, (2*nl(iqrp) + ns(iqrp)) < 0)
      end do

   end subroutine time_reversal_phases

   !----------------------------------------------------------------------------
   ! Write a log message to STDOUT with the time
   !----------------------------------------------------------------------------
   subroutine log_time(msg)
      implicit none

      character(len=*), intent(in) :: msg
      integer :: values(8)

      call date_and_time(values=values)
      write(*,'("[",1i2.2,":",1i2.2,":",1i2.2,"] ",a)') values(5), values(6), values(7), trim(msg)

   end subroutine log_time

   !----------------------------------------------------------------------------
   ! Determine how many files to open and write to, based upone the size of the
   ! two-quasiparticle basis. Fortran seems to cap the number of records
   ! at 2^31-1 (32-bit integer).
   !----------------------------------------------------------------------------
   function ab_file_count(basis_size) result(n)
      implicit none

      integer, intent(in) :: basis_size
      real(dp) :: real_basis_size, current_record_number
      integer  :: n

      real_basis_size = 1.0_dp*basis_size
      current_record_number = real_basis_size*(1+real_basis_size)
      n = ceiling(current_record_number/io_max_rec)

   end function ab_file_count

   !----------------------------------------------------------------------------
   ! Open N AB matrix files with name filename_base.i.ab
   !----------------------------------------------------------------------------
   subroutine open_ab_files(basis_size, recl, filename_base, status, unit)
      implicit none

      integer,           intent(in) :: basis_size, recl
      character(len=*),  intent(in) :: filename_base, status
      integer, optional, intent(in) :: unit

      integer :: ifile, ierr, nfiles, my_unit
      character(len=10) :: file_index
      character(len=len(filename_base)+10) :: filename

      nfiles = ab_file_count(basis_size)

      ! Added optional unit for the cutdown code which needs
      ! two sets of files to be open at the same time.
      if (present(unit)) then
         my_unit = unit
      else
         my_unit = io_ab_unit
      end if

      do ifile=1, nfiles
         ! filename_base.i.ab
         write(file_index, '(".",1i0)') ifile
         filename = trim(filename_base)//trim(file_index)//'.ab'

         open(my_unit+ifile-1, file=trim(filename), status=trim(status), form='unformatted', &
            access='direct', recl=recl, iostat=ierr)

         if (ierr /= 0) then
            write(*,'("Error opening AB matrix file ",1i0," (nfiles = ",1i0,").")') ifile, nfiles
            stop
         end if
      end do

   end subroutine open_ab_files

   !----------------------------------------------------------------------------
   ! Open N AB matrix files with name filename_base.i.ab using MPI IO
   !----------------------------------------------------------------------------
   subroutine mpi_open_ab_files(basis_size, filename_base, comm, mode, fh_list)

      use mpi

      implicit none

      integer,              intent(in)  :: basis_size, comm, mode
      character(len=*),     intent(in)  :: filename_base
      integer, allocatable, intent(out) :: fh_list(:)

      integer :: ifile, ierr, nfiles
      character(len=10) :: file_index
      character(len=len(filename_base)+10) :: filename

      !include 'mpif.h'

      nfiles = ab_file_count(basis_size)

      if (allocated(fh_list)) deallocate(fh_list)
      allocate(fh_list(nfiles))
      fh_list = 0

      do ifile=1, nfiles
         ! filename_base.i.ab
         write(file_index, '(".",1i0)') ifile
         filename = trim(filename_base)//trim(file_index)//'.ab'

         call mpi_file_open(comm, trim(filename), mode, mpi_info_null, fh_list(ifile), ierr)

         if (ierr /= 0) then
            write(*,*) 'IERR(MPI_FILE_OPEN) =', ierr
            call mpi_abort(comm, -1, ierr)
         end if

         call mpi_file_set_view(fh_list(ifile), int(0, kind=mpi_offset_kind), &
            mpi_double_precision, mpi_double_precision, 'native', mpi_info_null, ierr)

         if (ierr /= 0) then
            write(*,*) 'IERR(MPI_FILE_SET_VIEW) =', ierr
            call mpi_abort(comm, -1, ierr)
         end if
      end do

   end subroutine mpi_open_ab_files

   !----------------------------------------------------------------------------
   ! Close the AB matrix files
   !----------------------------------------------------------------------------
   subroutine close_ab_files(basis_size, unit)
      implicit none

      integer,           intent(in) :: basis_size
      integer, optional, intent(in) :: unit

      integer :: ifile, nfiles, my_unit

      ! Added optional unit for the cutdown code which needs
      ! two sets of files to be open at the same time.
      if (present(unit)) then
         my_unit = unit
      else
         my_unit = io_ab_unit
      end if

      nfiles = ab_file_count(basis_size)
      do ifile=1, nfiles
         close(my_unit+ifile-1)
      end do

   end subroutine close_ab_files

   !----------------------------------------------------------------------------
   ! Close the AB matrix files opened with MPI
   !----------------------------------------------------------------------------
   subroutine mpi_close_ab_files(fh_list)

      use mpi

      implicit none

      integer, allocatable, intent(inout) :: fh_list(:)
      integer :: ifile, ierr

      !include 'mpif.h'


      do ifile=1, size(fh_list)
         call mpi_file_close(fh_list(ifile), ierr)
      end do

      deallocate(fh_list)

   end subroutine mpi_close_ab_files

   !----------------------------------------------------------------------------
   ! Calculate matrix elements of a single-particle operator.
   ! Essentially a wrapper for extfield.init_external_field.
   !
   ! Unless an operator is specified via op, the routine assumes a generic
   ! beta-minus operator with the correct K and parity.
   !----------------------------------------------------------------------------
   subroutine init_sp_operator(k, parity, op, f_out)
      use constants,           only : translate_uppercase
      use external_field_type, only : external_field
      use extfield,            only : init_external_field, init_fam_mapping
      implicit none

      integer,                    intent(in)  :: k, parity
      character(len=*), optional, intent(in)  :: op
      type(external_field),       intent(out) :: f_out

      character(len=80) :: op_

      if (parity == 1) then
         op_ = 'GT-'
      else if (parity == -1) then
         op_ = 'RS2-'
      else
         stop 'Error with parity'
      end if

      ! If we specified an operator, don't simply use an empty operator
      if (present(op)) then
         op_ = adjustl(op)
         call translate_uppercase(op_)
         call init_external_field(label=op_, k=k, op=f_out)
      else
         f_out%label = ""
         f_out%k = k
         f_out%rank = 0
         f_out%beta_minus = .true.

         if (parity == 1) then
            f_out%parity_even = .true.
         else if (parity == -1) then
            f_out%parity_even = .false.
         else
            stop 'Incorrect parity'
         end if

         call init_fam_mapping(f_out)
         f_out%mat%elem(:) = 0
      end if

   end subroutine init_sp_operator

   !----------------------------------------------------------------------------
   ! Convert a single-particle matrix f_sp to a quasiparticle one f_qp_out,
   ! selecting which secotr of an HFB matrix is desired (20, 02, 11, 11~).
   !
   ! Always returns a p-n matrix (i.e. we transpose a beta+ matrix)
   !----------------------------------------------------------------------------
   subroutine init_qp_operator(class, f_sp, hfb, f_qp_out)
      use constants,           only : translate_uppercase
      use hfb_solution_type,   only : hfb_solution
      use external_field_type, only : external_field
      use blockmatrix_type,    only : blockmatrix, allocate_blockmatrix, triprod
      implicit none

      character(len=*),     intent(in)  :: class
      type(external_field), intent(in)  :: f_sp
      type(hfb_solution),   intent(in)  :: hfb
      type(blockmatrix),    intent(out) :: f_qp_out

      integer :: dim
      character(len=len_trim(class)) :: class_

      class_ = trim(class)
      call translate_uppercase(class_)

      dim = size(f_sp%mat%elem)
      call allocate_blockmatrix(f_qp_out, dim)
      f_qp_out%elem = 0

      if (f_sp%beta_minus) then
         if (class == '20') then
            call triprod('t', hfb%up, 'n', f_sp%mat, 'n', hfb%vn,  1d0, 0d0, f_qp_out)
         else if (class == '02') then
            call triprod('t', hfb%vp, 'n', f_sp%mat, 'n', hfb%un, -1d0, 0d0, f_qp_out)
         else if (class == '11') then
            call triprod('t', hfb%up, 'n', f_sp%mat, 'n', hfb%un,  1d0, 0d0, f_qp_out)
         else if (class == '11T' .or. class == '11~') then
            call triprod('t', hfb%vp, 'n', f_sp%mat, 'n', hfb%vn, -1d0, 0d0, f_qp_out)
         else
            stop 'Unknown class'
         end if
      else
         if (class == '20') then
            call triprod('t', hfb%vp, 't', f_sp%mat, 'n', hfb%un, -1d0, 0d0, f_qp_out)
         else if (class == '02') then
            call triprod('t', hfb%up, 't', f_sp%mat, 'n', hfb%vn,  1d0, 0d0, f_qp_out)
         else if (class == '11') then
            call triprod('t', hfb%vp, 't', f_sp%mat, 'n', hfb%vn, -1d0, 0d0, f_qp_out)
         else if (class == '11T' .or. class == '11~') then
            call triprod('t', hfb%up, 't', f_sp%mat, 'n', hfb%un,  1d0, 0d0, f_qp_out)
         else
            stop 'Unknown class'
         end if
      end if

   end subroutine init_qp_operator

   !----------------------------------------------------------------------------
   ! Generate a set of two-quasiparticle states (n,p) which correspond to a
   ! given sector of an HFB matrix (20, 02, 11, 11~) and respect a given K,
   ! parity, and cutoff energy.
   !
   ! Return a list basis_out(1:2, 1:dim)
   ! basis(1,i) = neutron iqp
   ! basis(2,i) = proton iqp
   !----------------------------------------------------------------------------
   subroutine make_two_qp_basis(class, k, parity, ecut, hfb, basis_out)
      use constants,           only : translate_uppercase
      use external_field_type, only : external_field
      use blockmatrix_type,    only : blockmatrix, deallocate_blockmatrix
      use hfb_solution_type,   only : hfb_solution
      use hfbtho_basis,        only : pwi_active_p, pwi_active_n
      implicit none

      character(len=*),     intent(in)  :: class
      integer,              intent(in)  :: k, parity
      real(dp),             intent(in)  :: ecut
      type(hfb_solution),   intent(in)  :: hfb
      integer, allocatable, intent(out) :: basis_out(:,:)

      integer :: iqpn, iqpp, n2qp
      integer, allocatable :: basis_tmp(:,:)

      type(external_field) :: f_sp
      type(blockmatrix)    :: f_qp
      character(len=3)     :: class_

      class_ = adjustl(class)
      call translate_uppercase(class_)

      ! Reference operator --- using default operator for basis
      call init_sp_operator(k=k, parity=parity, f_out=f_sp)

      ! Q.P. operator
      call init_qp_operator(class=trim(class_), f_sp=f_sp, hfb=hfb, f_qp_out=f_qp)

      ! Basis with cutoffs
      call get_qp_indices(f_qp, basis_tmp)

      ! Count possible two-qp states below PWI ("active") and below cutoff
      n2qp = 0
      do iqpn=1, size(basis_tmp, 1)
         if ((.not.pwi_active_n(iqpn)) .or. hfb%en(iqpn) > ecut) cycle
         do iqpp=1, size(basis_tmp, 2)
            if ((.not.pwi_active_p(iqpp)) .or. hfb%ep(iqpp) > ecut) cycle
            if (basis_tmp(iqpn, iqpp) /= 0) n2qp = n2qp+1
         end do
      end do

      if (allocated(basis_out)) deallocate(basis_out)
      allocate(basis_out(2, n2qp))
      basis_out(:,:) = 0

      ! Build the final array
      n2qp = 0
      do iqpn=1, size(basis_tmp, 1)
         if ((.not.pwi_active_n(iqpn)) .or. hfb%en(iqpn) > ecut) cycle
         do iqpp=1, size(basis_tmp, 2)
            if ((.not.pwi_active_p(iqpp)) .or. hfb%ep(iqpp) > ecut) cycle
            if (basis_tmp(iqpn, iqpp) /= 0) then
               n2qp = n2qp+1
               basis_out(1, n2qp) = iqpn
               basis_out(2, n2qp) = iqpp
            end if
         end do
      end do

      deallocate(basis_tmp)
      call deallocate_blockmatrix(f_sp%mat)
      call deallocate_blockmatrix(f_qp)

   end subroutine make_two_qp_basis

   !----------------------------------------------------------------------------
   ! Take in a (pn) basis and return the time-reversed basis
   !----------------------------------------------------------------------------
   subroutine time_reverse_basis(in, out)
      implicit none
      integer,              intent(in)  :: in(:,:)
      integer, allocatable, intent(out) :: out(:,:)

      integer :: i

      if (allocated(out)) deallocate(out)
      allocate(out(2, size(in, 2)))
      out(:,:) = 0

      do i=1, size(in, 2)
         out(1,i) = qp_reverse(in(1,i))
         out(2,i) = qp_reverse(in(2,i))
      end do

   end subroutine time_reverse_basis

   !----------------------------------------------------------------------------
   ! Calculate the residual interaction between two-quasiparticle states.
   ! Input:
   !    class      - which quadrant of the density is modified ('20', '02', '11', '11~')
   !    hfb        - HFB solution (U, V, E)
   !    ik         - matrix element to perturb
   !    delta_r_qp - blockmatrix holding the q.p. structure for the perturbed density
   !    f_sp       - blockmatrix with the s.p. operator structure
   !    f20_qp     - blockmatrix with the (20) q.p. operator structure
   !    f02_qp     - blockmatrix with the (02) q.p. operator structure
   !
   ! Output:
   !    h20_vector - vector which holds the output H(20)
   !----------------------------------------------------------------------------
   subroutine residual_interaction(class, hfb, ik, delta_r_qp, f_sp, f20_qp, f02_qp, h20_vector)
      use constants,         only : translate_uppercase
      use blockmatrix_type,  only : blockmatrix, allocate_blockmatrix, deallocate_blockmatrix, &
                                    copy_block_structure, triprod
      use hfb_solution_type, only : hfb_solution
      use density_set_type,  only : density_set, allocate_density, deallocate_density, zero_density
      use hamiltonian,       only : density, meanfield, pairingfield
      implicit none

      character(len=1),      intent(in)  :: class
      type(hfb_solution),    intent(in)  :: hfb
      integer,               intent(in)  :: ik
      type(blockmatrix),     intent(in)  :: delta_r_qp, f_sp, f20_qp, f02_qp
      real(dp), allocatable, intent(out) :: h20_vector(:)

      integer           :: i
      character(len=1)  :: class_
      type(density_set) :: rhox
      type(blockmatrix) :: reh_np, imh_np, reh_pn, imh_pn, rerho_pn, imrho_pn,  &
         rerho_np, imrho_np, redp, imdp, redm, imdm, rekp, imkp, rekm, imkm,    &
         delta_r_re, h20re


      ! Common allocations
      call allocate_blockmatrix(reh_np,   size(f20_qp%elem))
      call allocate_blockmatrix(imh_np,   size(f20_qp%elem))
      call allocate_blockmatrix(reh_pn,   size(f20_qp%elem))
      call allocate_blockmatrix(imh_pn,   size(f20_qp%elem))
      call allocate_blockmatrix(rerho_pn, size(f20_qp%elem))
      call allocate_blockmatrix(imrho_pn, size(f20_qp%elem))
      call allocate_blockmatrix(rerho_np, size(f20_qp%elem))
      call allocate_blockmatrix(imrho_np, size(f20_qp%elem))
      call allocate_blockmatrix(redp,     size(f20_qp%elem))
      call allocate_blockmatrix(imdp,     size(f20_qp%elem))
      call allocate_blockmatrix(redm,     size(f20_qp%elem))
      call allocate_blockmatrix(imdm,     size(f20_qp%elem))
      call allocate_blockmatrix(rekp,     size(f20_qp%elem))
      call allocate_blockmatrix(imkp,     size(f20_qp%elem))
      call allocate_blockmatrix(rekm,     size(f20_qp%elem))
      call allocate_blockmatrix(imkm,     size(f20_qp%elem))
      call allocate_blockmatrix(h20re,    size(f20_qp%elem))
      call allocate_density(rhox)

      ! Common configuration
      call copy_block_structure(f_sp,   reh_pn)
      call copy_block_structure(f_sp,   imh_pn)
      call copy_block_structure(f_sp,   reh_np)
      call copy_block_structure(f_sp,   imh_np)
      call copy_block_structure(f20_qp, redp)
      call copy_block_structure(f20_qp, imdp)
      call copy_block_structure(f02_qp, redm)
      call copy_block_structure(f02_qp, imdm)
      call copy_block_structure(f20_qp, h20re)

      ! Init
      reh_pn%elem   = 0;  imh_pn%elem   = 0;  reh_np%elem   = 0;  imh_np%elem   = 0
      redp%elem     = 0;  imdp%elem     = 0;  redm%elem     = 0;  imdm%elem     = 0
      rerho_pn%elem = 0;  imrho_pn%elem = 0;  rerho_np%elem = 0;  imrho_np%elem = 0
      rekp%elem     = 0;  imkp%elem     = 0;  rekm%elem     = 0;  imkm%elem     = 0
      h20re%elem    = 0
      call zero_density(rhox)

      ! Construct the perturbed density based on the \delta R we are using
      call allocate_blockmatrix(delta_r_re, size(f20_qp%elem))
      call copy_block_structure(delta_r_qp, delta_r_re)
      delta_r_re%elem(:)  = 0
      delta_r_re%elem(ik) = 1.0_dp

      class_ = class
      call translate_uppercase(class_)

      select case (class_)
         case ('X')
            ! rho_np      = -un x(-)'  vp'
            ! rho_pn      =  up x(-)   vn'
            ! kappa(+)_np = -un x(-)'  up'
            ! kappa(-)_np =  vn x(-)*' vp'
            call triprod('n', hfb%un, 't', delta_r_re, 't', hfb%vp, -1d0, 0d0, rerho_np)
            call triprod('n', hfb%up, 'n', delta_r_re, 't', hfb%vn,  1d0, 0d0, rerho_pn)
            call triprod('n', hfb%un, 't', delta_r_re, 't', hfb%up, -1d0, 0d0, rekp)
            call triprod('n', hfb%vn, 't', delta_r_re, 't', hfb%vp,  1d0, 0d0, rekm)

         case ('Y')
            ! rho_np =       vn y(-)'  up'
            ! rho_pn =      -vp y(-)   un'
            ! kappa(+)_np =  vn y(-)'  vp'
            ! kappa(-)_np = -un y(-)*' up'
            call triprod('n', hfb%vn, 't', delta_r_re, 't', hfb%up,  1d0, 0d0, rerho_np)
            call triprod('n', hfb%vp, 'n', delta_r_re, 't', hfb%un, -1d0, 0d0, rerho_pn)
            call triprod('n', hfb%vn, 't', delta_r_re, 't', hfb%vp,  1d0, 0d0, rekp)
            call triprod('n', hfb%un, 't', delta_r_re, 't', hfb%up, -1d0, 0d0, rekm)

         case ('P')
            ! rho_pn      =  up p(-)   un'
            ! rho_np      = -vn p(-)'  vp'
            ! kappa(+)_np = -vn p(-)'  up'
            ! kappa(-)_np =  un p(-)'* vp'
            call triprod('n', hfb%up, 'n', delta_r_re, 't', hfb%un,  1d0, 0d0, rerho_pn)
            call triprod('n', hfb%vn, 't', delta_r_re, 't', hfb%vp, -1d0, 0d0, rerho_np)
            call triprod('n', hfb%vn, 't', delta_r_re, 't', hfb%up, -1d0, 0d0, rekp)
            call triprod('n', hfb%un, 't', delta_r_re, 't', hfb%vp,  1d0, 0d0, rekm)

         case ('Q')
            ! rho_pn      = -vp q(-)   vn'
            ! rho_np      =  un q(-)'  up'
            ! kappa(+)_np =  un q(-)'  vp'
            ! kappa(-)_np = -vn q(-)'* up'
            call triprod('n', hfb%vp, 'n', delta_r_re, 't', hfb%vn, -1d0, 0d0, rerho_pn)
            call triprod('n', hfb%un, 't', delta_r_re, 't', hfb%up,  1d0, 0d0, rerho_np)
            call triprod('n', hfb%un, 't', delta_r_re, 't', hfb%vp,  1d0, 0d0, rekp)
            call triprod('n', hfb%vn, 't', delta_r_re, 't', hfb%up, -1d0, 0d0, rekm)

         case default
            write(*,'("ERROR: unrecognized class in mod_common.residual_interaction --- ",a)') class_
            stop
      end select

      call deallocate_blockmatrix(delta_r_re)

      ! Density matrices and mean field
      call density(rerho_pn, imrho_pn, rekp, imkp, rhox)
      call meanfield(rhox, reh_pn, imh_pn)
      call pairingfield(rhox, redp, imdp)

      call density(rerho_np, imrho_np, rekm, imkm, rhox)
      call meanfield(rhox, reh_np, imh_np)
      call pairingfield(rhox, redm, imdm)

      ! Test assumption of real matrix elements
      do i=1, size(imh_pn%elem)
         if (abs(imh_pn%elem(i)) > 1d-16 .or. &
             abs(imh_np%elem(i)) > 1d-16 .or. &
             abs(imdp%elem(i))   > 1d-16 .or. &
             abs(imdm%elem(i))   > 1d-16) then
            write(*,'(a)') 'ERROR: h_sp was not strictly real in mod_common.residual_interaction'
            stop
         end if
      end do

      ! Put stuff together
      call triprod('t', hfb%up, 'n', reh_pn, 'n', hfb%vn,  1d0, 0d0, h20re)
      call triprod('t', hfb%vp, 't', reh_np, 'n', hfb%un, -1d0, 1d0, h20re)
      call triprod('t', hfb%up, 'n', redp,   'n', hfb%un,  1d0, 1d0, h20re)
      call triprod('t', hfb%vp, 'n', redm,   'n', hfb%vn, -1d0, 1d0, h20re)

      if (allocated(h20_vector)) deallocate(h20_vector)
      allocate(h20_vector(size(f_sp%elem)))
      h20_vector(:) = H20re%elem(:)

      ! Common deallocations
      call deallocate_blockmatrix(reh_np)
      call deallocate_blockmatrix(imh_np)
      call deallocate_blockmatrix(reh_pn)
      call deallocate_blockmatrix(imh_pn)
      call deallocate_blockmatrix(rerho_pn)
      call deallocate_blockmatrix(imrho_pn)
      call deallocate_blockmatrix(rerho_np)
      call deallocate_blockmatrix(imrho_np)
      call deallocate_blockmatrix(redp)
      call deallocate_blockmatrix(imdp)
      call deallocate_blockmatrix(redm)
      call deallocate_blockmatrix(imdm)
      call deallocate_blockmatrix(rekp)
      call deallocate_blockmatrix(imkp)
      call deallocate_blockmatrix(rekm)
      call deallocate_blockmatrix(imkm)
      call deallocate_blockmatrix(h20re)
      call deallocate_density(rhox)

   end subroutine residual_interaction

   !----------------------------------------------------------------------------
   ! Determine the quasiparticle level to block using the HFB solution
   ! Always block K > 0.
   !----------------------------------------------------------------------------
   function find_blocked_level(bk, bparity, bt, hfb) result(iqp_out)
      use constants,         only : IT_PROTON, IT_NEUTRON, hfbtho_parity
      use hfb_solution_type, only : hfb_solution
      use hfbtho_basis,      only : lowest_energy_qp
      implicit none

      integer,            intent(in) :: bk, bparity, bt
      type(hfb_solution), intent(in) :: hfb

      integer :: iqp_out
      logical :: block_hole, block_particle

      iqp_out = -1

      ! Partial safeguard against using HFBTHO parity integers unexpectedly
      if (bparity /= 1 .and. bparity /= -1) then
         write(*,'(a,1i0)') 'ERROR in mod_common.find_blocked_level: unknown parity integer --- ', bparity
         stop
      end if

      ! Find the lowest-energy quasiparticles which match specified K and parity.
      ! The lowest-lying one should be correct, but check with HFBTHO.
      ! If K < 0, block the (quasi-)hole state with K = |K|.
      ! If K > 0, block the (quasi-)particle state with K = |K|.
      block_hole     = .false.
      block_particle = .false.
      if (bk > 0) then
         block_particle = .true.
      else
         block_hole = .true.
      end if

      if (bt == IT_NEUTRON) then
         iqp_out = lowest_energy_qp(energies=hfb%en, v=hfb%vn, it=bt, mask_2k=abs(bk), &
            mask_parity=hfbtho_parity(bparity), mask_unocc=block_particle, mask_occ=block_hole)
      else if (bt == IT_PROTON) then
         iqp_out = lowest_energy_qp(energies=hfb%ep, v=hfb%vp, it=bt, mask_2k=abs(bk), &
            mask_parity=hfbtho_parity(bparity), mask_unocc=block_particle, mask_occ=block_hole)
      end if

   end function find_blocked_level

   !----------------------------------------------------------------------------
   ! Remove all references of a set of levels from a two-qp basis
   ! and return a copy with the new data.
   !----------------------------------------------------------------------------
   subroutine remove_levels_from_basis(in, levels, it, out)
      use constants, only : IT_PROTON, IT_NEUTRON
      implicit none

      integer,              intent(in)  :: in(:,:), levels(:), it
      integer, allocatable, intent(out) :: out(:,:)

      logical :: found
      integer :: i, j, nout
      integer, allocatable :: tmp(:,:)

      ! For each basis state in in(:,:), add it to tmp(:,:)
      ! if it should not be excluded
      allocate(tmp(2, size(in, 2)))
      nout = 0;  tmp(:,:) = 0

      do i=1, size(in, 2)
         found = .false.

         do j=1, size(levels)
            if ((it == IT_NEUTRON .and. in(1,i) == levels(j)) .or. &
                (it == IT_PROTON  .and. in(2,i) == levels(j))) then
               found = .true.
               exit
            end if
         end do

         if (.not.found) then
            nout = nout+1
            tmp(1, nout) = in(1, i)
            tmp(2, nout) = in(2, i)
         end if
      end do

      if (allocated(out)) deallocate(out)
      allocate(out(2, nout))
      out(:,:) = tmp(:, 1:nout)
      deallocate(tmp)

   end subroutine remove_levels_from_basis

   !----------------------------------------------------------------------------
   ! Essentially the reverse of remove_level_...
   !
   ! Given an index and isospin, only return basis states which
   ! involve the level
   !----------------------------------------------------------------------------
   subroutine find_1to1_transitions(class, ib, tb, k, parity, ecut, hfb, basis_out)
      use constants, only : IT_PROTON, IT_NEUTRON
      use hfb_solution_type, only : hfb_solution
      implicit none

      character(len=*),     intent(in)  :: class
      integer,              intent(in)  :: ib, tb, k, parity
      real(dp),             intent(in)  :: ecut
      type(hfb_solution),   intent(in)  :: hfb
      integer, allocatable, intent(out) :: basis_out(:,:)

      integer :: i, nme, kt, kti
      integer, allocatable :: basis_tmp(:,:)


      call make_two_qp_basis(class=class, k=k, parity=parity, ecut=ecut, &
         hfb=hfb, basis_out=basis_tmp)

      ! Index in the basis to fix (kti is justthe "other" index)
      if (tb == IT_NEUTRON) then
         kt  = 1
         kti = 2
      else
         kt  = 2
         kti = 1
      end if

      ! Count states
      nme = 0
      do i=1, size(basis_tmp, 2)
         if (basis_tmp(kt, i) == ib) nme = nme+1
      end do

      allocate(basis_out(2, nme))
      basis_out = 0

      ! Save states
      nme = 0
      do i=1, size(basis_tmp, 2)
         if (basis_tmp(kt, i) == ib) then
            nme = nme+1
            basis_out(kt,  nme) = ib
            basis_out(kti, nme) = basis_tmp(kti, i)
         end if
      end do

      deallocate(basis_tmp)

   end subroutine find_1to1_transitions

   !----------------------------------------------------------------------------
   ! Generate all the operators we have for a specific K and parity.
   ! Return a list F(op 1(-), op 2(-), ... op n(-), op 1(+), ....).
   !----------------------------------------------------------------------------
   subroutine setup_all_operators(class, k, parity, hfb, fqp_out)
      use blockmatrix_type,    only : blockmatrix
      use external_field_type, only : external_field
      use hfb_solution_type,   only : hfb_solution
      implicit none

      character(len=*),               intent(in)  :: class
      integer,                        intent(in)  :: k, parity
      type(hfb_solution),             intent(in)  :: hfb
      type(blockmatrix), allocatable, intent(out) :: fqp_out(:)

      type(external_field) :: f

      ! Check
      if (k < 0) write(*,'(a)') 'Notice: K < 0 has not been rigorously tested.'

      if (parity == 1) then
         select case (k)
            case (0)
               allocate(fqp_out(4))
               call init_sp_operator(k=k, parity=1, op='f-',  f_out=f)
               call init_qp_operator(class=class, f_sp=f, hfb=hfb, f_qp_out=fqp_out(1))
               call init_sp_operator(k=k, parity=1, op='gt-', f_out=f)
               call init_qp_operator(class=class, f_sp=f, hfb=hfb, f_qp_out=fqp_out(2))
               call init_sp_operator(k=k, parity=1, op='f+',  f_out=f)
               call init_qp_operator(class=class, f_sp=f, hfb=hfb, f_qp_out=fqp_out(3))
               call init_sp_operator(k=k, parity=1, op='gt+', f_out=f)
               call init_qp_operator(class=class, f_sp=f, hfb=hfb, f_qp_out=fqp_out(4))

            case (1, -1)
               allocate(fqp_out(2))
               call init_sp_operator(k=k, parity=1, op='gt-', f_out=f)
               call init_qp_operator(class=class, f_sp=f, hfb=hfb, f_qp_out=fqp_out(1))
               call init_sp_operator(k=k, parity=1, op='gt+', f_out=f)
               call init_qp_operator(class=class, f_sp=f, hfb=hfb, f_qp_out=fqp_out(2))

         end select
      else if (parity == -1) then
         select case (k)
            case (0)
               allocate(fqp_out(12))
               call init_sp_operator(k=k, parity=-1, op='rs0-', f_out=f)
               call init_qp_operator(class=class, f_sp=f, hfb=hfb, f_qp_out=fqp_out(1))
               call init_sp_operator(k=k, parity=-1, op='ps0-', f_out=f)
               call init_qp_operator(class=class, f_sp=f, hfb=hfb, f_qp_out=fqp_out(2))
               call init_sp_operator(k=k, parity=-1, op='rs1-', f_out=f)
               call init_qp_operator(class=class, f_sp=f, hfb=hfb, f_qp_out=fqp_out(3))
               call init_sp_operator(k=k, parity=-1, op='r-',   f_out=f)
               call init_qp_operator(class=class, f_sp=f, hfb=hfb, f_qp_out=fqp_out(4))
               call init_sp_operator(k=k, parity=-1, op='p-',   f_out=f)
               call init_qp_operator(class=class, f_sp=f, hfb=hfb, f_qp_out=fqp_out(5))
               call init_sp_operator(k=k, parity=-1, op='rs2-', f_out=f)
               call init_qp_operator(class=class, f_sp=f, hfb=hfb, f_qp_out=fqp_out(6))
               call init_sp_operator(k=k, parity=-1, op='rs0+', f_out=f)
               call init_qp_operator(class=class, f_sp=f, hfb=hfb, f_qp_out=fqp_out(7))
               call init_sp_operator(k=k, parity=-1, op='ps0+', f_out=f)
               call init_qp_operator(class=class, f_sp=f, hfb=hfb, f_qp_out=fqp_out(8))
               call init_sp_operator(k=k, parity=-1, op='rs1+', f_out=f)
               call init_qp_operator(class=class, f_sp=f, hfb=hfb, f_qp_out=fqp_out(9))
               call init_sp_operator(k=k, parity=-1, op='r+',   f_out=f)
               call init_qp_operator(class=class, f_sp=f, hfb=hfb, f_qp_out=fqp_out(10))
               call init_sp_operator(k=k, parity=-1, op='p+',   f_out=f)
               call init_qp_operator(class=class, f_sp=f, hfb=hfb, f_qp_out=fqp_out(11))
               call init_sp_operator(k=k, parity=-1, op='rs2+', f_out=f)
               call init_qp_operator(class=class, f_sp=f, hfb=hfb, f_qp_out=fqp_out(12))

            case (1, -1)
               allocate(fqp_out(8))
               call init_sp_operator(k=k, parity=-1, op='rs1-', f_out=f)
               call init_qp_operator(class=class, f_sp=f, hfb=hfb, f_qp_out=fqp_out(1))
               call init_sp_operator(k=k, parity=-1, op='r-',   f_out=f)
               call init_qp_operator(class=class, f_sp=f, hfb=hfb, f_qp_out=fqp_out(2))
               call init_sp_operator(k=k, parity=-1, op='p-',   f_out=f)
               call init_qp_operator(class=class, f_sp=f, hfb=hfb, f_qp_out=fqp_out(3))
               call init_sp_operator(k=k, parity=-1, op='rs2-', f_out=f)
               call init_qp_operator(class=class, f_sp=f, hfb=hfb, f_qp_out=fqp_out(4))
               call init_sp_operator(k=k, parity=-1, op='rs1+', f_out=f)
               call init_qp_operator(class=class, f_sp=f, hfb=hfb, f_qp_out=fqp_out(5))
               call init_sp_operator(k=k, parity=-1, op='r+',   f_out=f)
               call init_qp_operator(class=class, f_sp=f, hfb=hfb, f_qp_out=fqp_out(6))
               call init_sp_operator(k=k, parity=-1, op='p+',   f_out=f)
               call init_qp_operator(class=class, f_sp=f, hfb=hfb, f_qp_out=fqp_out(7))
               call init_sp_operator(k=k, parity=-1, op='rs2+', f_out=f)
               call init_qp_operator(class=class, f_sp=f, hfb=hfb, f_qp_out=fqp_out(8))

            case (2, -2)
               allocate(fqp_out(2))
               call init_sp_operator(k=k, parity=-1, op='rs2-', f_out=f)
               call init_qp_operator(class=class, f_sp=f, hfb=hfb, f_qp_out=fqp_out(1))
               call init_sp_operator(k=k, parity=-1, op='rs2+', f_out=f)
               call init_qp_operator(class=class, f_sp=f, hfb=hfb, f_qp_out=fqp_out(2))
         end select
      end if

   end subroutine setup_all_operators

   !----------------------------------------------------------------------------
   ! Read X, Y, and QRPA energies into allocatable vectors
   !----------------------------------------------------------------------------
   subroutine read_xye(file, n_basis, x_out, y_out, e_out)
      implicit none

      character(len=*),      intent(in)  :: file
      integer,               intent(in)  :: n_basis
      real(dp), allocatable, intent(out) :: x_out(:,:), y_out(:,:), e_out(:)

      integer :: recl, ierr, i, nxy


      if (allocated(x_out)) deallocate(x_out)
      if (allocated(y_out)) deallocate(y_out)
      if (allocated(e_out)) deallocate(e_out)

      ! Storage is [ENERGY, X(1:N), Y(1:N)]
      inquire(iolength=recl) [(1.0_dp, i=1, 1 + 2*n_basis)]

      open(21, file=trim(file), status='old', access='direct', form='unformatted', &
         recl=recl, iostat=ierr)

      if (ierr /= 0) then
         write(*,'(2a)') "ERROR in mod_common.read_xye: could not open file ", trim(file)
         stop
      end if

      ! Count records
      i = 0;  nxy = 0

      do
         i = i+1
         read(21, rec=i, iostat=ierr)

         if (ierr /= 0) then
            exit
         else
            nxy = nxy+1
         end if
      end do

      ! Read in X and Y
      allocate(x_out(n_basis, nxy), y_out(n_basis, nxy), e_out(nxy))
      x_out(:,:) = 0;  y_out(:,:) = 0;  e_out(:) = 0

      do i=1, nxy
         read(21,rec=i) e_out(i), x_out(:,i), y_out(:,i)

         if (abs(e_out(i)) < 1d-12) then
            write(*,'(1x,"ERROR: an eigenvalue is zero! irec = ",1i0,", EQRPA(i) = ",1f0.6)')  i, e_out(i)
            stop
         end if
      end do

      close(21)

   end subroutine read_xye

   !----------------------------------------------------------------------------
   ! Header for the output files
   !----------------------------------------------------------------------------
   subroutine write_main_header(config, infile, nmpi, nomp, label)
      use config_type,  only : matrix_config
      implicit none

      type(matrix_config), intent(in) :: config
      integer,             intent(in) :: nmpi, nomp
      character(len=*),    intent(in) :: infile, label

      integer, parameter :: txtw = 100

      integer :: date_values(8)
      character(len=80)  :: date

      write(*,'(a)') repeat('*', txtw)
      write(*,'(a)') txtc(width=txtw,text='')
      write(*,'(a)') txtc(width=txtw,text='Quasiparticle Random-Phase Approximation for Even and Odd Nuclei')
      write(*,'(a)') txtc(width=txtw,text='Applying the Charge-Changing Finite Amplitude Method (PnFAM)')
      write(*,'(a)') txtc(width=txtw,text='')
      write(*,'(a)') txtc(width=txtw,text='Thomas Shafer and Mika T. Mustonen')
      write(*,'(a)') txtc(width=txtw,text='The University of North Carolina at Chapel Hill')
      write(*,'(a)') txtc(width=txtw,text='')

      call date_and_time(values=date_values)
      write(date,'(1i0,"/",1i0,"/",1i0,1x,1i0,":",1i2.2)') date_values(2), date_values(3), &
         date_values(1), date_values(5), date_values(6)

      write(*,'(a)') txtc(width=txtw,text='Run Date: '//trim(date))
      write(*,'(a)') txtc(width=txtw,text='')
      write(*,'(a/)') repeat('*', txtw)

      write(*,'(a)') 'USING NAMELIST FILE:'
      write(*,'(a)') repeat('-', txtw/3)
      write(*,'(a)')  trim(infile)
      write(*,'(a)')  ''

      write(*,'(a)') 'Namelist &FILENAMES'
      write(*,'(a)') repeat('-', txtw/3)
      write(*,'(a," : ",a)') 'HFB  file name ...', trim(config%file_hfb)
      write(*,'(a," : ",a)') 'Base file name ...', trim(config%file_basename)
      write(*,'(a," : ",a)') 'AB   file name ...', trim(config%file_ab)
      write(*,'(a," : ",a)') 'XY   file name ...', trim(config%file_xy)
      write(*,'(a," : ",a)') '1QP  file name ...', trim(config%file_me_one)
      write(*,'(a," : ",a)') '2QP  file name ...', trim(config%file_me_two)
      write(*,'(a," : ",a)') 'VV   file name ...', trim(config%file_vv)
      write(*,'(a)') ''

      write(*,'(a)') 'Namelist &BASIS_OPTIONS'
      write(*,'(a)') repeat('-', txtw/3)
      write(*,'(a," : ",1i2)')          'K .........', config%basis_k
      write(*,'(a," : ",1i2)')          'Parity ....', config%basis_parity
      write(*,'(a," : ",1f5.1," MeV")') 'Ecutoff ...', config%basis_cutoff
      write(*,'(a)') ''

      write(*,'(a)') 'Namelist &DIAGONALIZATION_OPTIONS'
      write(*,'(a)') repeat('-', txtw/3)
      write(*,'(a," : ",1f5.1)') 'Max. eigenvalue ...', config%diag_ev_limit
      write(*,'(a," : ",1i0)')   'Blocking factor ...', config%diag_blocking
      write(*,'(a," : ",1i0)')   'Proc. per row .....', config%diag_nprow
      write(*,'(a)') ''

      write(*,'(a)') 'Namelist &MATRIX_ELEMENT_OPTIONS'
      write(*,'(a)') repeat('-', txtw/3)
      write(*,'(a," : ",a)') 'Check XY orthonormality ...', trim(merge('Yes', 'No ', config%me_check_norm))
      write(*,'(a," : ",a)') 'Check XY completeness .....', trim(merge('Yes', 'No ', config%me_check_complete))
      write(*,'(a," : ",a)') 'Disable 1-QP correction ...', trim(merge('Yes', 'No ', config%me_disable_correction))
      write(*,'(a)') ''

      write(*,'(a)') 'Namelist &ODD_NUCLEUS_OPTIONS'
      write(*,'(a)') repeat('-', txtw/3)
      write(*,'(a," : ",1i2,1i4)')      'Blocking (N) .........', config%odd_block_n
      write(*,'(a," : ",1i2,1i4)')      'Blocking (P) .........', config%odd_block_p
      write(*,'(a," : ")',advance='no') 'Even file basename ...'
      if (len_trim(config%odd_even_basename) == 0) then
         write(*,'(a)') '<None>'
         write(*,'(a," : ",a)')         'Even AB filenames ....', '<None>'
         write(*,'(a," : ",a)')         'Even XY filename .....', '<None>'
      else
         write(*,'(a)') trim(config%odd_even_basename)
         write(*,'(a," : ",a)')         'Even AB filenames ....', trim(config%odd_even_file_ab)
         write(*,'(a," : ",a)')         'Even XY filename .....', trim(config%odd_even_file_xy)
      end if
      write(*,'(a)') ''

      write(*,'(a)') repeat('*', txtw)
      write(*,'(a)') txtc('BEGINNING PROGRAM: '//trim(label), txtw)
      write(*,'(a)') repeat('*', txtw)
      write(*,'(a)') ''

      write(*,'(a)') 'PARALLELIZATION'
      write(*,'(a)') repeat('-', txtw/3)
      write(*,'(a," : ",1i0)') 'No. OpenMP threads ...', nomp
      write(*,'(a," : ",1i0)') 'No. MPI processes ....', nmpi
      write(*,'(a)') ''


   end subroutine write_main_header

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
   ! Display some information abotu a quasiparticle level
   !----------------------------------------------------------------------------
   subroutine display_qp_level(iqp, it, hfb)
      use constants,         only : strpar, IT_PROTON, IT_NEUTRON
      use blockmatrix_type,  only : nb, db
      use hfb_solution_type, only : hfb_solution
      use hfbtho_basis,      only : nl, ns, npar, get_hfbtho_occupations
      implicit none

      integer, intent(in) :: iqp, it
      type(hfb_solution), intent(in) :: hfb

      integer  :: is, ib
      real(dp) :: occ(size(hfb%ep))

      if (it /= IT_NEUTRON .and. it /= IT_PROTON) then
         write(*,'(a,1i0)') 'Error: unrecognized isospin -- ', it
         stop
      end if

      ! Find which block and state in the block for comparing with HFBTHO
      is = 0
      do ib=1, nb
         if ((is+db(ib)) > iqp) then
            exit
         else
            is = is+db(ib)
         end if
      end do
      is = iqp-sum(db(1:ib-1))

      ! Lower norm
      if (it == IT_NEUTRON) then
         call get_hfbtho_occupations(hfb%vn, occ)
      else if (it == IT_PROTON) then
         call get_hfbtho_occupations(hfb%vp, occ)
      end if

      write(*,'(2a)') 'Isospin ................ : ', merge('P', 'N', it == IT_PROTON)
      write(*,'(a,1i0," [block=",1i0,", state=",1i0,"]")') 'Index .................. : ', iqp, ib, is
      write(*,'(a,1i2,"/2",a)') 'Quantum numbers ........ :', 2*nl(iqp)+ns(iqp), strpar(npar(iqp))
      write(*,'(a,1f7.4," MeV")') 'Energy ................. :', merge(hfb%ep(iqp), hfb%en(iqp), it == IT_PROTON)
      write(*,'(a,1f9.6)') 'Lower norm ............. :', occ(iqp)

   end subroutine display_qp_level

   !----------------------------------------------------------------------------
   ! Calculate the 1-to-1 g.s. energy
   !----------------------------------------------------------------------------
   subroutine display_1to1_energy(iqp, it, hfb, sign)
      use constants,         only : IT_NEUTRON, IT_PROTON
      use hfbtho_basis,      only : get_hfbtho_occupations, lowest_energy_qp
      use hfb_solution_type, only : hfb_solution
      implicit none

      integer, intent(in) :: iqp, it, sign
      type(hfb_solution), intent(in) :: hfb

      integer  :: iqp2
      real(dp) :: eqp1, eqp2, occ(size(hfb%ep))
      logical  :: cut_occ, cut_unocc

      cut_occ   = .false.
      cut_unocc = .false.

      if (sign > 0) then
         cut_unocc = .true.
      else
         cut_occ = .true.
      end if

      if (it == IT_PROTON) then
         eqp1 = hfb%ep(iqp)
         iqp2 = lowest_energy_qp(hfb%en, hfb%vn, IT_NEUTRON, mask_k_positive=.true., &
            mask_occ=cut_occ, mask_unocc=cut_unocc)
         eqp2 = hfb%en(iqp2)
         call get_hfbtho_occupations(hfb%vn, occ)
      else if (it == IT_NEUTRON) then
         eqp1 = hfb%en(iqp)
         iqp2 = lowest_energy_qp(hfb%ep, hfb%vp, IT_PROTON, mask_k_positive=.true., &
            mask_occ=cut_occ, mask_unocc=cut_unocc)
         eqp2 = hfb%ep(iqp2)
         call get_hfbtho_occupations(hfb%vp, occ)
      else
         write(*,'(a)') 'ERROR in mod_common.display_1to1_energy: unknown isospin integer'
         stop
      end if

      write(*,'(a,1f7.4," MeV")') "Min. 1-to-1 energy ..... :", eqp2-eqp1
      write(*,'(a,1f9.6)')        "Final lower norm ....... :", occ(iqp2)

   end subroutine display_1to1_energy

end module pnfam_matrix
