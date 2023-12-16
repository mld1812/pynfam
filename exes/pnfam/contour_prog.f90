!------------------------------------------------------------------------------
!> This program wraps the primary pnfam program to solve the fam equations for
!> a set of complex energy points (an energy contour). Different contour shapes
!> are implemented through the fam_mode input.
!>
!> @authors M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
!> @authors E.M. Ney, UNC Chapel Hill, 2018-
!------------------------------------------------------------------------------
program contour_main
   use pnfam_solver, only : pnfam_solve
   use contour_setup

   implicit none
   integer, parameter :: dp = kind(1d0)
   character(len=80), parameter :: file_contour  = "pnfam_CONTOUR.dat"

   integer :: i, j

   call init_ctr_namelist
   call read_ctr_namelist(file_contour, .true.)

   call init_pnfam_namelist
   call read_pnfam_namelist(fam_input_filename, .false.)

   ! Setup the pnFAM tasks
   call setup_contour

   ! Run the calculations
   do i=1, size(tasks)
      call apply_task_values(tasks(i))
      call pnfam_solve
      tasks(i)%re_str = re_str
      tasks(i)%im_str = im_str
      tasks(i)%conv = max(iter_conv/abs(iter_conv),0)
   end do

   ! Record summary file
   do i=1, size(ops_list)
      operator_name = ops_list(i); operator_k = kval_list(i)
      write(ctr_output_txtfile,'(2a,"K",i0,".out")') &
             trim(operator_name),beta_type,operator_k
      call write_ctr_output_header
      do j=1, size(tasks)
         if (tasks(j)%op_name==operator_name .and. &
             tasks(j)%op_k==operator_k) then
            call write_ctr_output_point(tasks(j))
         end if
      end do
      call write_ctr_output_footer
   end do

end program contour_main
