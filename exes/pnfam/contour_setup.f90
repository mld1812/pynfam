!------------------------------------------------------------------------------
!> This module contains the variables and routines used by the contour program.
!>
!> @authors M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
!> @authors E.M. Ney, UNC Chapel Hill, 2018-
!------------------------------------------------------------------------------
module contour_setup
   use pnfam_setup
   use pnfam_logger
   implicit none
   integer, parameter, private :: dp = kind(1d0)

   real(dp) :: energy_start, energy_step, half_width
   integer :: nr_points, operator_active(5)
   character(len=2) :: operator_groups(5)
   character(len=200) :: fam_input_filename
   character(len=12) :: fam_mode

   namelist /ctr_general/ fam_mode, fam_input_filename
   namelist /ctr_extfield/ operator_groups, operator_active
   namelist /str_parameters/ energy_start, energy_step, nr_points, half_width

   integer :: fout=11
   character(len=200) :: ctr_output_txtfile, ctr_parameter_filename

   real(dp), allocatable :: real_eqrpa_ctr(:), imag_eqrpa_ctr(:)
   integer, allocatable :: kval_list(:)
   character(len=10), allocatable :: ops_list(:)

   type ctr_task
      ! Inputs
      real(dp) :: real_eqrpa, imag_eqrpa
      integer  :: i, op_k
      character(len=10) :: op_name
      ! Outputs
      real(dp), allocatable :: re_str(:), im_str(:)
      integer :: conv
   end type ctr_task

   type(ctr_task), allocatable :: tasks(:)

contains
   subroutine setup_contour
      implicit none
      integer :: i

      ! Setup the output file procedures
      openlog  => openlog_stdout
      writelog => writelog_stdout
      closelog => closelog_stdout
      abort    => abort_nompi

      ! Open the log for error reporting
      write(ctr_output_txtfile,'(a,"_errors.dat")') trim(fam_output_filename)
      call openlog(output_txtfile, print_stdout)

      call setup_energy_contour
      call setup_operators
      call setup_tasks

      ! Close the log and delete if no errors were written to it
      call closelog('delete')

   end subroutine


   !---------------------------------------------------------------------------
   subroutine init_ctr_namelist
      implicit none

      energy_start = 0
      energy_step = 0.2_dp
      half_width = 0.5_dp
      nr_points = 5
      operator_groups = (/ '0+', '1+', '0-', '1-', '2-' /)
      operator_active = (/ 0, 0, 0, 0, 0 /)

   end subroutine

   subroutine read_ctr_namelist(fn,use_cli)
      implicit none
      character(len=*), intent(in) :: fn
      logical, intent(in) :: use_cli
      integer, parameter :: fnml = 12

      ctr_parameter_filename = fn
      if (use_cli) then
         if (command_argument_count() == 1) then
            call get_command_argument(1, parameter_filename, status=ierr)
            if (ierr /= 0) then
               write(*,'(1x,A)') "Error reading the command-line argument (too long?)"
               call abort
            end if
         else
            if (command_argument_count() > 1) then
               write(*,'(1x,A)') "More than one command-line argument! Aborting."
               call abort
            end if
         end if
      end if

      open(fnml,file=ctr_parameter_filename,iostat=ierr)

      read(fnml,nml=ctr_general,iostat=ierr)
      read(fnml,nml=ctr_extfield,iostat=ierr)
      read(fnml,nml=str_parameters,iostat=ierr)

      close(fnml)

      if (ierr /= 0) then
         write(*,'(1x,a,i0)') 'ERROR: problem reading contour namelist:', ierr
         call abort
      end if

   end subroutine


   !---------------------------------------------------------------------------
   subroutine setup_energy_contour
      implicit none

      integer :: i

      allocate(real_eqrpa_ctr(nr_points))
      allocate(imag_eqrpa_ctr(nr_points))
      real_eqrpa_ctr = 0
      imag_eqrpa_ctr = 0
      do i=1, nr_points
         real_eqrpa_ctr(i) = energy_start + (i-1)*energy_step
         imag_eqrpa_ctr(i) = half_width
      end do
   end subroutine

   subroutine setup_operators
      implicit none
      integer :: nops, i, j
      integer, parameter :: nops_per_group(5) = (/ 1,2,2,6,3 /)

      ! Count number of operators to compute
      nops = 0
      do i = 1, size(operator_groups)
         if (operator_active(i)/=0) nops = nops + nops_per_group(i)
      end do

      if (nops==0) then
         ! If none active, compute the operator specified in pnfam_NAMELIST.dat
         nops = 1
         allocate(ops_list(nops), kval_list(nops))
         ops_list(1) = operator_name; kval_list(1) = operator_k
      else
         ! Jpi   N_ops  Operator(K)
         ! '0+'  1      F(1)
         ! '1+'  2      GT(0),  GT(1)
         ! '0-'  2      RS0(0), PS0(0)
         ! '1-'  6      R(0),   R(1),   P(0),   P(1),   RS1(0), RS1(1)
         ! '2-'  3      RS2(0), RS2(1), RS2(2)
         allocate(ops_list(nops), kval_list(nops))
         j=0
         do i = 1, size(operator_groups)
            if (operator_active(i)/=0) then
               select case(i)
               case(1)
                  ops_list(j+1) = 'F'  ; kval_list(j+1) = sign(0,operator_k)
               case(2)
                  ops_list(j+1) = 'GT' ; kval_list(j+1) = sign(0,operator_k)
                  ops_list(j+2) = 'GT' ; kval_list(j+2) = sign(1,operator_k)
               case(3)
                  ops_list(j+1) = 'RS0'; kval_list(j+1) = sign(0,operator_k)
                  ops_list(j+2) = 'PS0'; kval_list(j+2) = sign(0,operator_k)
               case(4)
                  ops_list(j+1) = 'R'  ; kval_list(j+1) = sign(0,operator_k)
                  ops_list(j+2) = 'R'  ; kval_list(j+2) = sign(1,operator_k)
                  ops_list(j+3) = 'P'  ; kval_list(j+3) = sign(0,operator_k)
                  ops_list(j+4) = 'P'  ; kval_list(j+4) = sign(1,operator_k)
                  ops_list(j+5) = 'RS1'; kval_list(j+5) = sign(0,operator_k)
                  ops_list(j+6) = 'RS1'; kval_list(j+6) = sign(1,operator_k)
               case(5)
                  ops_list(j+1) = 'RS2'; kval_list(j+1) = sign(0,operator_k)
                  ops_list(j+2) = 'RS2'; kval_list(j+2) = sign(1,operator_k)
                  ops_list(j+3) = 'RS2'; kval_list(j+3) = sign(2,operator_k)
               end select
               j = j + nops_per_group(i)
            end if
         end do
      end if

   end subroutine


   !---------------------------------------------------------------------------
   subroutine write_ctr_output_header
      use pnfam_constants, only : txtr
      implicit none
      integer :: ixterms
      character(len=160) :: st

      open(fout,file=ctr_output_txtfile)

      write(fout,'(2A)')             "# pnFAM code version: ", version
      write(fout,'(2A)')             "# pnFAM code commit:  ", commit
      write(fout,'(2A)')             "# Contour input file name: ", trim(ctr_parameter_filename)
      write(fout,'(2A)')             "# Residual interaction: ", trim(interaction_name)
      write(fout,'(3A,I0)')          "# Operator: ", trim(operator_name), " with K = ", operator_k
      write(fout,'(A,F0.4)')         "# Gamma (half-width): ", half_width
      write(fout,'(A)')              "#"

      write(fout,'(5A)',advance='no') "#", txtr("Conv", 5), txtr("Re(EQRPA)", 29), &
         txtr("Re(Strength)", 34), txtr("Im(Strength)", 34)
      do ixterms=2, size(re_str)
        write(fout,'(2A)',advance='no') txtr("Re("//trim(g(ixterms-1)%label)//")", 34), &
           txtr("Im("//trim(g(ixterms-1)%label)//")", 34)
      end do
      write(fout,'(A)') ""
   end subroutine

   subroutine write_ctr_output_point(task)
      implicit none
      type(ctr_task), intent(in) :: task
      integer :: istr
      write(fout,'(1i5)',advance='no') task%conv
      write(fout,'(1F30.19)',advance='no') task%real_eqrpa
      do istr=1, size(task%re_str)
         write(fout,'(2ES34.19)',advance='no') task%re_str(istr), task%im_str(istr)
      end do
      write(fout,'(A)') ""
      flush(fout)
   end subroutine

   subroutine write_ctr_output_footer
      implicit none
      close(fout)
   end subroutine


   !---------------------------------------------------------------------------
   subroutine setup_tasks
      implicit none
      integer :: i, j, npts, nops, ipt

      npts = size(real_eqrpa_ctr)
      nops = size(ops_list)

      allocate(tasks(npts*nops))

      ipt = 0
      do i=1,nops
         do j=1,npts
            ipt = ipt + 1
            tasks(ipt)%real_eqrpa = real_eqrpa_ctr(j)
            tasks(ipt)%imag_eqrpa = imag_eqrpa_ctr(j)
            tasks(ipt)%op_name    = ops_list(i)
            tasks(ipt)%op_k       = kval_list(i)
            tasks(ipt)%i          = j-1 ! start at zero
         end do
      end do
   end subroutine

   subroutine apply_task_values(task)
      implicit none
      type(ctr_task) :: task

      real_eqrpa    = task%real_eqrpa
      imag_eqrpa    = task%imag_eqrpa
      operator_name = task%op_name
      operator_k    = task%op_k

      if (fam_output_filename /= "") then
         write(fam_output_filename,'(2a,"K",i0,"_",i6.6)') &
             trim(operator_name),beta_type,operator_k,task%i
      end if

   end subroutine



end module
