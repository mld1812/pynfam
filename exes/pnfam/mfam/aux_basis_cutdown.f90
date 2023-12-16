!------------------------------------------------------------------------------
! pnfam_matrix/aux_basis_cutdown.f90
!
! T. Shafer, UNC Chapel Hill, 2015
!------------------------------------------------------------------------------
program pnfam_matrix_cutdown
   use logger
   use blockmatrix_type,  only : blockmatrix
   use hfbtho_basis,      only : get_hfbtho_solution
   use hfb_solution_type, only : hfb_solution
   use config_type,       only : matrix_config, read_namelist_from_file
   use pnfam_matrix,      only : io_max_rec, input_file_from_cli, open_ab_files,   &
                                 close_ab_files, ab_file_count, make_two_qp_basis, &
                                 write_main_header, log_time
   implicit none

   integer, parameter :: si = kind(1)
   integer, parameter :: di = selected_int_kind(11)
   integer, parameter :: dp = kind(1.0d0)

   ! Configuration
   type(matrix_config) :: config
   type(hfb_solution)  :: hfb_soln

   ! Storage
   integer :: recl, nfiles, ierr
   character(len=210) :: infile = '', filename_ab_new = ''

   ! New parameter
   character(len=10) :: char_new_ecutoff
   real(dp) :: new_ecutoff

   ! Two-quasiparticle basis (forward-basis only!)
   integer, allocatable :: twoqp_orig(:,:), twoqp_new(:,:)
   integer :: num_ab_orig, num_ab_new

   ! Matrices
   integer     :: ic1, ic2, inc1, ipc1, inc2, ipc2, imult, if_orig, if_new, frec_new, frec_orig
   integer(di) :: irec_orig, irec_new
   real(dp)    :: tmpval
   logical     :: found

   !----------------------------------------------------------------------------
   !----------------------------------------------------------------------------

   ! Logging
   openlog  => openlog_stdout
   writelog => writelog_stdout
   closelog => closelog_stdout
   call openlog('')

   ! Check for correct integer precision
   if (di == -1) then
      write(*,'(1x,a)') 'ERROR! this machine does not support integers of&
         & selected_int_kind(11)!'
      stop
   end if

   ! Configuration and setup
   call input_file_from_cli(infile, default='matrix.in')

   open(1, file=trim(infile), status='old', access='sequential', iostat=ierr)
   if (ierr /= 0) stop 'Error opening Namelist file'

   call read_namelist_from_file(fh=1, config=config)

   call get_hfbtho_solution(fn=trim(config%file_hfb), &
      ep=hfb_soln%ep, vp=hfb_soln%vp, up=hfb_soln%up, &
      en=hfb_soln%en, vn=hfb_soln%vn, un=hfb_soln%un)

   close(1)

   ! Read in new cutoff
   write(*,'(a,1f0.1)') 'Original Ecutoff = ', config%basis_cutoff

   write(*,'(a)') 'New cutoff?'
   read(*,*) new_ecutoff

   write(char_new_ecutoff,'(1f10.1)') new_ecutoff
   filename_ab_new = trim(config%file_basename)//'-cutdown-'//trim(adjustl(char_new_ecutoff))//'.n.ab'

   ! Header
   write(*,*)
   call write_main_header(config, infile=infile, nmpi=0, nomp=0, &
      label='QRPA Basis Truncation')

   write(*,'(a)') 'TRUNCATION PARAMETERS'
   write(*,'(a)') repeat('-', 33)
   write(*,'(a,1f5.1)') 'Old Ecutoff ............... :', config%basis_cutoff
   write(*,'(a,1f5.1)') 'New Ecutoff ............... :', new_ecutoff
   write(*,'(a,a)')     'Old matrix storage ........ : ', trim(config%file_ab)
   write(*,'(a,a)')     'New matrix storage ........ : ', trim(filename_ab_new)
   write(*,*)

   ! Make original 2QP basis
   call make_two_qp_basis(class='20', k=config%basis_k, parity=config%basis_parity, &
      ecut=config%basis_cutoff, hfb=hfb_soln, basis_out=twoqp_orig)
   call make_two_qp_basis(class='20', k=config%basis_k, parity=config%basis_parity, &
      ecut=new_ecutoff, hfb=hfb_soln, basis_out=twoqp_new)

   num_ab_orig = size(twoqp_orig,2)
   num_ab_new  = size(twoqp_new, 2)

   inquire(iolength=recl) 1.0_dp

   write(*,'(a,1i0)') 'Two-QP dimension (orig) ... : ', size(twoqp_orig,2)
   write(*,'(a,1i0)') 'Two-QP dimension (new) .... : ', size(twoqp_new, 2)
   write(*,'(a,1i0)') 'Record length (DP) ........ : ', recl

   nfiles = ab_file_count(num_ab_orig)
   write(*,'(a,1i0)')       'No. output files (orig) ... : ', nfiles
   nfiles = ab_file_count(num_ab_new)
   write(*,'(a,1i0)')       'No. output files (new) .... : ', nfiles

   write(*,*)
   call log_time('Masking basis states')

   ! Mask out the basis states we don't want anymore
   do ic1=1, num_ab_orig
      inc1 = twoqp_orig(1,ic1)
      ipc1 = twoqp_orig(2,ic1)
      found = .false.
      do ic2=1, num_ab_new
         inc2 = twoqp_new(1,ic2)
         ipc2 = twoqp_new(2,ic2)
         if (inc2 == inc1 .and. ipc2 == ipc1) then
            found = .true.
            exit
         end if
      end do
      if (.not. found) then
         twoqp_orig(1,ic1) = -1
         twoqp_orig(2,ic1) = -1
      end if
   end do

   call log_time('Reading and storing')

   ! Reading and storage
   call open_ab_files(size(twoqp_orig,2), recl, config%file_basename, 'old', unit=60)
   call open_ab_files(size(twoqp_new,2),  recl, &
      trim(config%file_basename)//'-cutdown-'//trim(adjustl(char_new_ecutoff)), &
      'replace', unit=70)

   irec_orig = 0
   irec_new  = 0

   ! This allows iteration over both A (1->num_ab**2) and B(num_ab**2+1 -> 2*num_ab**2)
   do imult=0, 1
      do ic1=1, num_ab_orig
         inc1 = twoqp_orig(1,ic1)
         ipc1 = twoqp_orig(2,ic1)

         if (inc1 == -1 .or. ipc1 == -1) then
            irec_orig = irec_orig+(num_ab_orig-ic1+1)
            cycle
         end if

         do ic2=ic1, num_ab_orig
            inc2 = twoqp_orig(1,ic2)
            ipc2 = twoqp_orig(2,ic2)

            irec_orig = irec_orig+1

            if (inc2 == -1 .or. ipc2 == -1) then
               cycle
            else
               irec_new = irec_new+1

               if_orig = ceiling(1.0_dp*irec_orig/io_max_rec)
               if_new  = ceiling(1.0_dp*irec_new/io_max_rec)

               frec_orig = int(irec_orig - (if_orig-1)*io_max_rec, kind=si)
               frec_new  = int(irec_new  - (if_new -1)*io_max_rec, kind=si)

               read (60-1+if_orig,rec=frec_orig) tmpval
               write(70-1+if_new, rec=frec_new)  tmpval
            end if
         end do
      end do
   end do

   call close_ab_files(size(twoqp_orig,2), unit=60)
   call close_ab_files(size(twoqp_new,2),  unit=70)

end program pnfam_matrix_cutdown
