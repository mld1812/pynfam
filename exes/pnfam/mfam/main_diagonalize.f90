!------------------------------------------------------------------------------
! pnfam_matrix/main_diagonalize.f90
!
! T. Shafer, UNC Chapel Hill, 2015
!------------------------------------------------------------------------------
program pnfam_matrix_diagonalize
   use logger
   use constants,         only : IT_PROTON, IT_NEUTRON
   use blockmatrix_type,  only : blockmatrix
   use hfbtho_basis,      only : get_hfbtho_solution
   use hfb_solution_type, only : hfb_solution
   use config_type,       only : matrix_config, read_namelist_from_file
   use pnfam_matrix,      only : ab_file_count, io_max_rec, qp_reverse, input_file_from_cli, &
                                 time_reversal_phases, log_time, mpi_close_ab_files,         &
                                 mpi_open_ab_files, make_two_qp_basis, find_blocked_level,   &
                                 remove_levels_from_basis, write_main_header, display_qp_level
   use mpi

   implicit none

   integer, parameter :: si = kind(1)
   integer, parameter :: di = selected_int_kind(11)
   integer, parameter :: dp = kind(1.0d0)

   ! Configuration
   type(matrix_config) :: config
   type(hfb_solution)  :: hfb_soln

   ! Storage
   integer(di) :: offset
   integer :: recl, nfiles, iab, iproc, ic, ir, ifile, irec
   integer, allocatable  :: fhandles(:)
   real(dp) :: temporary
   real(dp), allocatable :: amat(:), bmat(:), minusmat(:), plusmat(:)

   ! Two-quasiparticle *FORWARD* basis and orbital blocking
   integer, allocatable :: twoqp_ab(:,:), tphase_ab(:), twoqp_diag(:,:)
   integer :: num_ab, num_diag
   integer :: bo_index = -1, bo_isospin = -1

   ! Parallelization
   !include 'mpif.h'
   integer :: ierr, myrank, nprocs
   integer :: ictxt, nprow, npcol, mb, nb, rsrc, csrc, myrow, mycol
   integer :: istatus(mpi_status_size)

   ! Scalapack
   integer :: locc, locr, lld, info, local_nrows, local_ncols
   integer :: descab(9)
   integer, external :: numroc
   integer, allocatable :: ir_map_ab(:), ic_map_ab(:), ir_map_diag(:), ic_map_diag(:)

   ! Diagonalization
   real(dp), allocatable :: dmat(:), dvec(:), smat(:), cmat(:), wc(:), vc(:)
   real(dp), allocatable :: tm1(:), tm2(:), tm3(:), idmat(:)

   ! Final storage
   integer  :: ir_global, ic_global, col_comm, col_rank, row_comm, fh
   real(dp) :: norm
   real(dp), allocatable :: xy(:)

   ! Misc
   character(len=210) :: infile

   !----------------------------------------------------------------------------
   !----------------------------------------------------------------------------

   ! Logging
   openlog  => openlog_stdout
   writelog => writelog_stdout
   closelog => closelog_stdout
   call openlog('')

   ! Init parallel
   call mpi_init(ierr)
   call mpi_comm_rank(mpi_comm_world, myrank, ierr)
   call mpi_comm_size(mpi_comm_world, nprocs, ierr)

   ! Reading of namelists and setup on rank 0
   if (myrank == 0) then

      ! Check for correct integer precision
      if (di == -1) then
         write(*,'(1x,a)') 'ERROR! this machine does not support integers of&
            & selected_int_kind(11)!'
         call mpi_abort(mpi_comm_world, -1, ierr)
      end if

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
      call write_main_header(config, infile=infile, nmpi=nprocs, nomp=0, &
         label='QRPA Matrix Diagonalization')

      ! Basis for the A and B matrices
      call make_two_qp_basis(class='20', k=config%basis_k, parity=config%basis_parity, &
         ecut=config%basis_cutoff, hfb=hfb_soln, basis_out=twoqp_ab)

      call time_reversal_phases(basis=twoqp_ab, phases=tphase_ab)

      num_ab   = size(twoqp_ab, 2)
      num_diag = num_ab ! Size of the matrix to diagonalize ... could be overwritten below

      inquire(iolength=recl) 1.0_dp
      nfiles = ab_file_count(num_ab)

      write(*,'(a)')     'DIAGONALIZATION PARAMETERS'
      write(*,'(a)')     repeat('-', 33)
      write(*,'(a,1i0)') 'Two-QP dimension ........ : ', num_ab
      write(*,'(a,1i0)') 'Record length (DP) ...... : ', recl
      write(*,'(a,1i0)') 'No. AB output files ..... : ', nfiles
      write(*,*)

      ! Option to remove a level from the QRPA for the odd nucleus
      if (config%odd_block_n(1) /= 0 .or. config%odd_block_p(1) /= 0) then
         if (config%odd_block_n(1) /= 0 .and. config%odd_block_p(1) /= 0) then
            write(*,'(a)') "ERROR: the program can only handle a single blocked level."
            stop
         else
            ! Find the level to block
            if (config%odd_block_n(1) /= 0) then
               bo_isospin = IT_NEUTRON
               bo_index   = find_blocked_level(bk=config%odd_block_n(1), &
                  bparity=config%odd_block_n(2), bt=IT_NEUTRON, hfb=hfb_soln)
            else if (config%odd_block_p(1) /= 0) then
               bo_isospin = IT_PROTON
               bo_index   = find_blocked_level(bk=config%odd_block_p(1), &
                  bparity=config%odd_block_p(2), bt=IT_PROTON, hfb=hfb_soln)
            end if
         end if

         ! Update the basis
         if (bo_index < 1) then
            allocate(twoqp_diag(2, size(twoqp_ab, 2)))
            twoqp_diag(:,:) = twoqp_ab(:,:)
         else
            call remove_levels_from_basis(in=twoqp_ab, levels=[bo_index, qp_reverse(bo_index)], &
               it=bo_isospin, out=twoqp_diag)
            num_diag = size(twoqp_diag, 2)
         end if

         write(*,'(a,1x,a)') 'ODD NUCLEON/BLOCKING'
         write(*,'(a)') repeat('-', 33)
         call display_qp_level(iqp=bo_index, it=bo_isospin, hfb=hfb_soln)
         call display_qp_level(iqp=qp_reverse(bo_index), it=bo_isospin, hfb=hfb_soln)
         write(*,'(a,1i0)') 'New two-qp dimension ... : ', num_diag
         write(*,*)

         ! Now we will use the even-even AB matrices !!!
         if (len_trim(config%odd_even_basename) < 1) then
            write(*,'(a)') "ERROR: no even-solution basename was specified to load the A and B matrices."
            stop
         else
            write(*,'(3a,/)') "NB: Using even-even A and B matrices '", trim(config%odd_even_file_ab), &
               "' since blocking is activated."
         end if
      end if

      ! Check processor count
      if (nprocs /= config%diag_nprow**2) then
         write(*,'(/,"ERROR: nprocs /= nprow**2. nprocs = ", 1i0, ", nprow = ", 1i0)') nprocs, config%diag_nprow
         call mpi_abort(mpi_comm_world, -1, ierr)
      end if
   end if

   ! MPI distribution
   call distribute_calculation_details
   call mpi_barrier(mpi_comm_world, ierr)

   ! Scalapack init
   nprow = config%diag_nprow
   npcol = config%diag_nprow
   mb    = config%diag_blocking
   nb    = config%diag_blocking
   rsrc  = 0   ! Pin the first block on the (0,0) process
   csrc  = 0   ! Pin the first block on the (0,0) process

   call sl_init(ictxt, nprow, npcol)
   call blacs_gridinit(ictxt, 'R', nprow, npcol)
   call blacs_gridinfo(ictxt, nprow, npcol, myrow, mycol)


   ! Block-cyclic distribution of the A+-B matrices
   !----------------------------------------------------------------------------
   if (myrank == 0) call log_time('Reading in the QRPA matrix')

   if (bo_index /= -1) then
      call mpi_open_ab_files(num_ab, config%odd_even_basename, mpi_comm_world, mpi_mode_rdonly, fhandles)
   else
      call mpi_open_ab_files(num_ab, config%file_basename, mpi_comm_world, mpi_mode_rdonly, fhandles)
   end if

   ! Size of the matrix block on this process
   locr = numroc(num_diag, mb, myrow, rsrc, nprow)
   locc = numroc(num_diag, nb, mycol, csrc, npcol)
   lld  = locr

   ! Scalapack descriptor --- all matrices are distributed identically
   call descinit(descab, num_diag, num_diag, mb, nb, rsrc, csrc, ictxt, lld, info)

   if (info == -9) then
      lld = descab(9)
      write(*,'(3(a,1i0),a,/)') 'Notice: DESCINIT re-sized a local matrix on rank ', &
         myrank, ' so that LLD_A = ', descab(9), ' (was ', lld, ')'
   else if (info /= 0) then
      write(*,'(a,1i0)') 'ERROR: DESCINIT failed with INFO = ', info
      call mpi_abort(mpi_comm_world, -1, ierr)
   end if

   ! 2D block-cyclic distribution: all square matrices are distributed
   ! according to the same scheme. Global row and column indices are contained
   ! in ir_map_ab(ir) and ic_map_ab(ic) for the A,B matrices. Global row and
   ! column indices for the matrix to be diagonalized (<= A,B) are in i?_map_diag.
   call block_cyclic(local_nrows, local_ncols, ir_map_ab, ic_map_ab, ir_map_diag, ic_map_diag)

   if (local_nrows /= lld)  write(*,'(a)') 'WARNING: local_nrows /= lld  for (pr,pc) =', myrow, mycol
   if (local_ncols /= locc) write(*,'(a)') 'WARNING: local_ncols /= locc for (pr,pc) =', myrow, mycol

   ! Initially we need only A-B and A+B
   allocate(minusmat(locc*lld), plusmat(locc*lld), amat(locc*lld), bmat(locc*lld))
   minusmat(:) = 0; plusmat(:) = 0; amat(:) = 0; bmat(:) = 0

   ! Lower-triangular storage: iab=0 => A, iab=1 => B
   do iab=0, 1
      offset = iab*num_ab*(1_di+num_ab)/2_di
      iproc = 0
      do ic=1,locc
         do ir=1,locr
            iproc = iproc+1
            if (ir_map_ab(ir) < ic_map_ab(ic)) then
               cycle
            else
               ifile = ceiling((ir_map_ab(ir) + 0.5_dp*(2_di*num_ab-ic_map_ab(ic))*(ic_map_ab(ic)-1) + offset)/io_max_rec, kind=si)
               irec  = int(ir_map_ab(ir) + (2_di*num_ab-ic_map_ab(ic))*(ic_map_ab(ic)-1)/2 + offset - (ifile-1)*io_max_rec, kind=si)

               call mpi_file_read_at(fhandles(ifile), int(irec-1, kind=mpi_offset_kind), &
                  temporary, 1, mpi_double_precision, istatus, ierr)

               if (ierr /= 0) then
                  write(*,*) 'IERR(MPI_FILE_READ_AT) =', ierr, myrank
                  call mpi_abort(mpi_comm_world, -1, ierr)
               end if

               if (iab == 0) then
                  amat(iproc) = temporary
               else
                  bmat(iproc) = tphase_ab(ic_map_ab(ic))*temporary
               end if
            end if
         end do
      end do
   end do

   call mpi_barrier(mpi_comm_world, ierr)
   call mpi_close_ab_files(fhandles)

   minusmat(:) = amat(:)-bmat(:)
   plusmat(:)  = amat(:)+bmat(:)

   deallocate(amat, bmat)


   ! Diagonalization
   !----------------------------------------------------------------------------

   ! Step 1: diagonalize (A-B)
   ! smat  - orthogonal matrix of eigenvectors of (A-B) (distributed)
   ! dmat  - diagonal matrix with sqrt(<eigenvalue>) along the diagonal (distributed)
   ! idmat - diagonal matrix with 1/sqrt(<eigenvalue>) along the diagonal (distributed)
   ! dvec  - array(num_diag) containing <eigenvalues of (A-B)>
   if (myrank == 0) call log_time('Diagonalizing (A-B) matrix')

   allocate(smat(locc*lld), dmat(locc*lld), idmat(locc*lld))
   smat = 0; dmat = 0; idmat = 0

   ! Unfortunately, I have to bail out to ORFAC=1d-5 for large space :/
   call diagonalize_symmetric(minusmat, dvec, smat, orfac=1.0d-5)

   do ic=1, locc
      do ir=1, lld
         iproc = ir + (ic-1)*lld
         if (ic_map_diag(ic) == ir_map_diag(ir)) then
            dmat(iproc)  = sqrt(dvec(ir_map_diag(ir)))
            idmat(iproc) = 1.0_dp/sqrt(dvec(ir_map_diag(ir)))
         end if
      end do
   end do

   deallocate(dvec)  ! Allocated by diagonalize_symmetric

   ! Step 2: construct C = D**1/2 S' (A+B) S D**1/2
   ! In practice: dmat smat' (A+B) smat dmat
   if (myrank == 0) call log_time('Building C matrix')

   allocate(cmat(locc*lld), tm1(locc*lld), tm2(locc*lld), tm3(locc*lld))

   ! D**1/2 S' --- GEMM for S transpose
   if (myrank == 0) call log_time('... D**1/2 S')
   call pdgemm('N', 'T', num_diag, num_diag, num_diag, 1.0_dp, dmat, 1, 1, descab, &
      smat, 1, 1, descab, 0.0_dp, tm1, 1, 1, descab)

   ! (D**1/2 S') (A+B)
   if (myrank == 0) call log_time('... D**1/2 S (A+B)')
   call pdsymm('R', 'L', num_diag, num_diag, 1.0_dp, plusmat, 1, 1, descab, &
      tm1, 1, 1, descab, 0.0_dp, tm2, 1, 1, descab)

   ! (D**1/2 S' (A+B)) S
   if (myrank == 0) call log_time('... D**1/2 S (A+B) S')
   call pdgemm('N', 'N', num_diag, num_diag, num_diag, 1.0_dp, tm2, 1, 1, descab, &
      smat, 1, 1, descab, 0.0_dp, tm1, 1, 1, descab)

   ! (D**1/2 S' (A+B) S) D**1/2
   if (myrank == 0) call log_time('... D**1/2 S (A+B) S D**1/2')
   call pdsymm('R', 'L', num_diag, num_diag, 1.0_dp, dmat, 1, 1, descab, &
      tm1, 1, 1, descab, 0.0_dp, cmat, 1, 1, descab)

   ! Step 3: diagonalize C
   ! wc - eigenvalues of C (initially eqrpa**2)
   ! vc - eigenvectors of C (distributed)
   if (myrank == 0) call log_time('Diagonalizing C matrix')

   allocate(vc(locc*lld))

   vc = 0
   call diagonalize_symmetric(cmat, wc, vc, evmax=config%diag_ev_limit, orfac=1.0d-5)
   wc = sqrt(wc)

   ! Step 4: compute the quantities which form X and Y
   ! tm1 - S D**1/2 V
   ! tm2 - S D**(-1/2) V
   if (myrank == 0) call log_time('Untangling eigenvectors')

   ! D**1/2 V
   if (myrank == 0) call log_time('... D**1/2 V')
   call pdgemm('N', 'N', num_diag, num_diag, num_diag, 1.0_dp, dmat, 1, 1, descab, &
      vc, 1, 1, descab, 0.0_dp, tm3, 1, 1, descab)

   ! S (D**1/2 V)
   if (myrank == 0) call log_time('... S D**1/2 V')
   call pdgemm('N', 'N', num_diag, num_diag, num_diag, 1.0_dp, smat, 1, 1, descab, &
      tm3, 1, 1, descab, 0.0_dp, tm1, 1, 1, descab)

   ! S (D**(-1/2) V)
   if (myrank == 0) call log_time('... D**(-1/2) V')
   call pdgemm('N', 'N', num_diag, num_diag, num_diag, 1.0_dp, idmat, 1, 1, descab, &
      vc, 1, 1, descab, 0.0_dp, tm3, 1, 1, descab)

   ! S (D**(-1/2) V)
   if (myrank == 0) call log_time('... S D**(-1/2) V')
   call pdgemm('N', 'N', num_diag, num_diag, num_diag, 1.0_dp, smat, 1, 1, descab, &
      tm3, 1, 1, descab, 0.0_dp, tm2, 1, 1, descab)

   deallocate(smat, dmat, idmat, cmat, vc, tm3)

   ! Step 5: collect eigenvectors to form X and Y. Each process column
   ! will collect the sets of X and Y it holds as a way to reduce
   ! the overall memory cost.
   if (myrank == 0) call log_time('Storing the result')

   allocate(xy(2*num_diag))

   ! New communicator col_comm per-process-column
   call mpi_comm_split(mpi_comm_world, mycol, myrow, col_comm, ierr)
   call mpi_comm_rank(col_comm, col_rank, ierr)

   ! New communicator across top row for writing to XY file
   if (myrow == 0) then
      call mpi_comm_split(mpi_comm_world, 0, mycol, row_comm, ierr)
   else
      call mpi_comm_split(mpi_comm_world, mpi_undefined, mycol, row_comm, ierr)
   end if
   call mpi_barrier(mpi_comm_world, ierr)

   if (myrow == 0) then
      call mpi_file_open(row_comm, trim(config%file_xy), mpi_mode_wronly+mpi_mode_create, &
         mpi_info_null, fh, ierr)

      call mpi_file_set_view(fh, int(0, kind=mpi_offset_kind), &
         mpi_double_precision, mpi_double_precision, 'native', mpi_info_null, ierr)
   end if

   ! For each local column: (1) calculate my parts of X and Y,
   ! (2) reduce X and Y to row 0, and (3) store the result
   do ic=1, locc
      ic_global = ic_map_diag(ic)

      if (ic_global <= size(wc)) then
         xy(:) = 0

         do ir=1, lld
            ir_global = ir_map_diag(ir)
            iproc = ir + (ic-1)*lld

            xy(ir_global)          = 0.5_dp*(tm1(iproc)/sqrt(wc(ic_global)) + tm2(iproc)*sqrt(wc(ic_global)))
            xy(ir_global+num_diag) = 0.5_dp*(tm1(iproc)/sqrt(wc(ic_global)) - tm2(iproc)*sqrt(wc(ic_global)))
         end do

         call mpi_barrier(col_comm, ierr)
         if (col_rank == 0) then
            call mpi_reduce(mpi_in_place, xy, 2*num_diag, mpi_double_precision, mpi_sum, 0, col_comm, ierr)
         else
            call mpi_reduce(xy, xy, 2*num_diag, mpi_double_precision, mpi_sum, 0, col_comm, ierr)
         end if
         call mpi_barrier(col_comm, ierr)

         ! Re-normalize and collect
         if (myrow == 0) then
            norm = sum(xy(1:num_diag)**2 - xy(1+num_diag:2*num_diag)**2)

            if (abs(norm-1.0_dp) > 1.0d-10) then
              write(*,'(1x,"Strange init. normalization for ixy = ",1i0,": ",1es13.6)') ic_global, norm
            end if

            call mpi_file_write_at(fh, int((ic_global-1)*(1+2*num_diag),kind=mpi_offset_kind),&
               wc(ic_global),1,mpi_double_precision,istatus,ierr)

            call mpi_file_write_at(fh, int((ic_global-1)*(1+2*num_diag)+1,kind=mpi_offset_kind),&
               xy,2*num_diag,mpi_double_precision,istatus,ierr)
         end if
      end if
   end do

   if (myrow == 0) call mpi_file_close(fh,ierr)

   ! Finalize
   call blacs_gridexit(ictxt)
   call blacs_exit(1)
   call mpi_finalize(ierr)

contains

   !----------------------------------------------------------------------------
   ! MPI_BCAST everything from Rank 0
   !----------------------------------------------------------------------------
   subroutine distribute_calculation_details
      use blockmatrix_type, only : nb, db

      implicit none

      ! blockmatrix_type distribution
      call mpi_bcast(nb, 1,  mpi_integer, 0, mpi_comm_world, ierr)
      if (myrank /= 0) allocate(db(nb))
      call mpi_bcast(db, nb, mpi_integer, 0, mpi_comm_world, ierr)

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

      ! Other local
      call mpi_bcast(num_ab,   1, mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast(num_diag, 1, mpi_integer, 0, mpi_comm_world, ierr)

      if (myrank /= 0) allocate(twoqp_ab(2,num_ab), tphase_ab(2*num_ab))
      call mpi_bcast(twoqp_ab,  2*num_ab, mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast(tphase_ab, 2*num_ab, mpi_integer, 0, mpi_comm_world, ierr)

      ! QP removal
      call mpi_bcast(bo_index, 1, mpi_integer, 0, mpi_comm_world, ierr)
      if (bo_index /= -1) then
         if (myrank /= 0) allocate(twoqp_diag(2,num_diag))
         call mpi_bcast(twoqp_diag, 2*num_diag, mpi_integer, 0, mpi_comm_world, ierr)
         call mpi_bcast(bo_isospin, 1,          mpi_integer, 0, mpi_comm_world, ierr)
      end if

   end subroutine distribute_calculation_details

   !----------------------------------------------------------------------------
   ! Calculate the block-cyclic logic for this process
   !
   ! This ASSUMES the following globals have been set:
   ! mycol    - starting from 0, my process column
   ! myrow    - starting from 0, my process row
   ! npcol    - # of processes per column
   ! nprow    - # of processes per row
   ! nb       - blocking factor for columns
   ! mb       - blocking factor for rows
   ! num_ab   - dimension of A,B matrices
   ! num_diag - dimension of (possibly to-be-cutdown) QRPA matrix
   !
   ! It also assumes blocks are distributed in the normal way (row-major). The
   ! subroutine returns a few things:
   !     nrows_local: # rows on this process (check against lld or locr)
   !     ncols_local: # cols on this process (check against locc)
   !     *_ab:        for local irow or icol, what is irow or icol in the A,B matrices?
   !     *_diag:      for local irow or icol, what is irow or icol in the matrices to diagonalize?
   !
   ! If we don't exclude an odd orbital then *_diag = *_ab since the matrices are
   ! identical. If we DO exclude an orbital from the QRPA, *_ab are required to
   ! map the A,B matrix elements to our smaller matrix and *_diag are required
   ! to access elements of these smaller matrices (previously, *_ab handled both cases).
   !----------------------------------------------------------------------------
   subroutine block_cyclic(nrows_local, ncols_local, ir_ab, ic_ab, ir_diag, ic_diag)
      implicit none

      integer,              intent(out) :: nrows_local, ncols_local
      integer, allocatable, intent(out) :: ir_ab(:), ic_ab(:), ir_diag(:), ic_diag(:)

      integer :: nb_col_total, nb_row_total, ir_local, ic_local
      integer :: ib, db, ifound, iab, iqpp, iqpn
      logical :: skip


      ! Number of process block matrices which will be distributed on the grid
      nb_col_total = ceiling(1.0_dp*num_diag/nb)
      nb_row_total = ceiling(1.0_dp*num_diag/mb)

      ! How many columns and rows this rank's process block will handle.
      ! We get some number of full (nb * mb) blocks, plus maybe leftovers.
      nrows_local = 0
      ncols_local = 0

      ! E.g. blocks 0, NPCOL, 2*NPCOL, ...
      do ib=mycol, nb_col_total-1, npcol
         ! Full block
         if (ib /= (nb_col_total-1)) then
            ncols_local = ncols_local + nb
         ! Leftovers
         else
            ! Edge case where the "leftover" block could be a full block
            ! Remember num_diag <= num_ab is the dimension of the total matrix
            ! to be diagonalized
            if (mod(num_diag, nb) == 0) then
               ncols_local = ncols_local + nb
            else
               ncols_local = ncols_local + mod(num_diag, nb)
            end if
         end if
      end do

      do ib=myrow, nb_row_total-1, nprow
         if (ib /= (nb_row_total-1)) then
            nrows_local = nrows_local + mb
         else
            if (mod(num_diag, mb) == 0) then
               nrows_local = nrows_local + mb
            else
               nrows_local = nrows_local + mod(num_diag, mb)
            end if
         end if
      end do

      allocate(ir_ab(nrows_local), ic_ab(ncols_local), ir_diag(nrows_local), ic_diag(ncols_local))

      ir_ab = 0; ir_diag = 0; ir_local = 0
      ic_ab = 0; ic_diag = 0; ic_local = 0

      ! Mapping:
      ! Map a column or row on this process to a) the total A matrix and
      ! b) the reduced A matrix derived from removing a q.p. level or two.
      ! This is done by looping over the submatrix blocks
      ! Columns
      do ib=mycol, nb_col_total-1, npcol
         ! Number of columns in this block
         if (ib < (nb_col_total-1)) then
            db = nb
         else
            if (mod(num_diag, nb) == 0) then
               db = nb
            else
               db = mod(num_diag, nb)
            end if
         end if

         ! Mapping is easy if there are no excluded orbitals
         if (bo_index == -1) then
            do ifound=1, db
               ic_local = ic_local + 1
               ic_ab(ic_local)   = ib*nb + ifound
               ic_diag(ic_local) = ib*nb + ifound
            end do
         ! Harder with blocking
         else
            ! Now we have to loop over all possible 2QP states and exclude ones we don't want.
            ! This is slow but easier to keep track of
            ifound = 0
            do iab=1, num_ab
               iqpn = twoqp_ab(1, iab)
               iqpp = twoqp_ab(2, iab)

               skip = ((bo_isospin == IT_PROTON  .and. any([bo_index, qp_reverse(bo_index)] == iqpp)) &
                  .or. (bo_isospin == IT_NEUTRON .and. any([bo_index, qp_reverse(bo_index)] == iqpn)))

               ! If the level is excluded, go to the next one
               if (skip) then
                  cycle
               else
                  ifound = ifound + 1
                  ! Haven't made it to this block yet
                  if (ifound <= ib*nb) then
                     cycle
                  ! Past this block and ready for the next one
                  else if (ifound > (ib*nb + db)) then
                     exit
                  ! Assign the level
                  else
                     ic_local = ic_local + 1
                     ic_ab(ic_local)   = iab
                     ic_diag(ic_local) = ifound
                  end if
               end if
            end do
         end if
      end do

      ! Rows
      do ib=myrow, nb_row_total-1, nprow
         ! Number of rowumns in this block
         if (ib < (nb_row_total-1)) then
            db = mb
         else
            if (mod(num_diag, mb) == 0) then
               db = mb
            else
               db = mod(num_diag, mb)
            end if
         end if

         if (bo_index == -1) then
            do ifound=1, db
               ir_local = ir_local + 1
               ir_ab(ir_local)   = ib*mb + ifound
               ir_diag(ir_local) = ib*mb + ifound
            end do
         else
            ifound = 0
            do iab=1, num_ab
               iqpn = twoqp_ab(1, iab)
               iqpp = twoqp_ab(2, iab)

               skip = ((bo_isospin == IT_PROTON  .and. any([bo_index, qp_reverse(bo_index)] == iqpp)) &
                  .or. (bo_isospin == IT_NEUTRON .and. any([bo_index, qp_reverse(bo_index)] == iqpn)))

               if (skip) then
                  cycle
               else
                  ifound = ifound + 1
                  if (ifound <= ib*mb) then
                     cycle
                  else if (ifound > (ib*mb + db)) then
                     exit
                  else
                     ir_local = ir_local + 1
                     ir_ab(ir_local)   = iab
                     ir_diag(ir_local) = ifound
                  end if
               end if
            end do
         end if
      end do

   end subroutine block_cyclic

   !----------------------------------------------------------------------------
   ! Diagonalize a symmetric matrix using Scalapack
   !
   ! This ASSUMES the following globals:
   ! num_diag       - dimension of the matrix (2-qp basis size)
   ! descab         - descriptor for all NUM_DIAG x NUM_DIAG matrices
   ! npcol          - # of processes per column
   ! nprow          - # of processes per row
   ! ictxt          - BLACS/Scalapack context
   ! mpi_comm_world - MPI communicator
   !----------------------------------------------------------------------------
   subroutine diagonalize_symmetric(mat, w, v, evmax, orfac)
      implicit none

      real(dp),              intent(inout) :: mat(:)
      real(dp), allocatable, intent(out)   :: w(:)
      real(dp),              intent(out)   :: v(:)
      real(dp), optional,    intent(in)    :: evmax, orfac

      integer :: neigenvalues, neigenvectors, info, lwork, liwork
      integer, allocatable :: iwork(:), ifail(:), icluster(:)
      real(dp) :: abstol, vu_, orfac_
      real(dp), external :: pdlamch
      real(dp), allocatable :: work(:), gap(:), eigenvalues(:)
      character(len=1) :: evmode
      character(len=100) :: ostr

      ! Scalapack parameters
      abstol = 2*pdlamch(ictxt, 'S')
      orfac_ = 1.0d-4
      vu_    = 0.0_dp
      evmode = 'A'

      ! Workspace
      lwork  = int(5_di*get_lwork(), kind=si)
      liwork = get_liwork()

      ! Set mode
      if (present(evmax)) then
         vu_ = evmax**2
         evmode = 'V'
      end if

      ! Set orfac
      if (present(orfac)) then
         orfac_ = orfac
      end if

      if (myrank == 0) then
         write(ostr,'("... eigenvector mode = ",1a)') evmode
         call log_time(ostr)

         if (evmode == 'V') then
         write(ostr,'("... SQRT[max. eigenvalue] = ",1f0.4)') sqrt(vu_)
            call log_time(ostr)
         end if

         write(ostr,'("... orfac = ",1es7.1e2)') orfac_
         call log_time(ostr)

         write(ostr,'("... lwork = ",1i0)') lwork
         call log_time(ostr)

         write(ostr,'("... liwork = ",1i0)') liwork
         call log_time(ostr)
      end if

      allocate(eigenvalues(num_diag), ifail(num_diag), icluster(2*nprow*npcol), gap(nprow*npcol))
      allocate(work(lwork), iwork(liwork))

      eigenvalues = 0; ifail = 0; icluster = 0; gap = 0; work = 0; iwork = 0

      ! Diagonalize
      if (myrank == 0) call log_time('... diagonalizing')

      call pdsyevx('V', evmode, 'L', num_diag, mat, 1, 1, descab, 0d0, vu_, 0, 0, abstol, &
         neigenvalues, neigenvectors, eigenvalues, orfac_, v, 1, 1, descab, work, lwork, &
         iwork, liwork, ifail, icluster, gap, info)

      call blacs_barrier(ictxt, 'A')

      if (info /= 0) call handle_diagonalization_error(info, ifail, icluster, gap)

      if (neigenvalues /= neigenvectors) then
         write(*,'(1x,2(a,1i0),a)') 'Warning: neigenvalues != neigenvectors (neigenvalues = ', &
            neigenvalues, ', neigenvectors = ', neigenvectors, '.'
      end if

      if (myrank == 0) then
         write(ostr,'("... found ",1i0, " eigenvalues")') neigenvalues
         call log_time(ostr)
      end if

      if (allocated(w)) deallocate(w)
      allocate(w(neigenvalues))

      w = 0
      w(1:neigenvalues) = eigenvalues(1:neigenvalues)

      deallocate(eigenvalues, ifail, icluster, gap, work, iwork)

   end subroutine diagonalize_symmetric

   !----------------------------------------------------------------------------
   ! Compute LWORK as described in the Scalapack reference.
   !
   ! Assumes nprow = npcol
   ! Assumes the globals:
   !  - num_diag
   !  - nb
   !  - nprow
   !----------------------------------------------------------------------------
   function get_lwork() result(lwork)
      implicit none

      integer :: lwork
      integer :: nn, np0

      nn  = max(num_diag, nb, 2)
      np0 = numroc(nn, nb, 0, 0, nprow)

      lwork = 5*num_diag + max(5*nn, np0**2 + 2*nb**2) + nn*ceiling(1.0_dp*num_diag/nprow**2)

   end function get_lwork

   !----------------------------------------------------------------------------
   ! Compute LIWORK as described in the Scalapack reference.
   !
   ! Assumes nprow = npcol
   ! Assumes the globals:
   !  - num_diag
   !  - nprow
   !----------------------------------------------------------------------------
   function get_liwork() result(liwork)
      implicit none
      integer :: liwork
      liwork = 6*max(num_diag, 1+nprow**2, 4)

   end function get_liwork

   !----------------------------------------------------------------------------
   ! Error out and write the correct output files for Scalapack errors
   !----------------------------------------------------------------------------
   subroutine handle_diagonalization_error(info, ifail, icluster, gap)
      implicit none

      integer,  intent(in) :: info, ifail(:), icluster(:)
      real(dp), intent(in) :: gap(:)

      ! Global error check on rank 0
      if (myrank == 0 .and. info /= 0) then
         write(*,'(a,1i0,a)') 'ERROR: INFO = ', info, ' in PDSYGVX. INFO mods are:'

         if (info > 0) then
            write(*,'(5i4)') mod(info,2), mod(info/2,2), mod(info/4,2), mod(info/8,2), mod(info/16,2)
         end if

         ! If lack of convergence was an issue
         if (mod(info,2) /= 0) then
            open(81, file='ifail.err', status='unknown')
            write(81,*) ifail(:)
            close(81)
         end if

         ! If re-orthogonalization was an issue
         if (mod(info/2,2) /= 0) then
            open(80, file='icluster.err', status='unknown')
            open(82, file='gap.err', status='unknown')

            write(80,*) icluster(:)
            write(82,*) gap(:)

            close(80)
            close(82)
         end if

         call mpi_abort(mpi_comm_world, -1, ierr)
      end if

   end subroutine handle_diagonalization_error

end program pnfam_matrix_diagonalize
