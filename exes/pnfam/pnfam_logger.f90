!------------------------------------------------------------------------------
!> This module contains the interfaces for logging routines use in pnFAM.
!> These are implemented by the subroutines below for serial calls to pnfam.
!>
!> @authors M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
!> @authors E.M. Ney, UNC Chapel Hill, 2018-
!------------------------------------------------------------------------------
module pnfam_logger
   use iso_fortran_env
   implicit none
   
   ! By default, we call the standard output routines
   procedure(openlog_interface),  pointer :: openlog
   procedure(writelog_interface), pointer :: writelog
   procedure(closelog_interface), pointer :: closelog
   procedure(abort_interface),    pointer :: abort

   character(len=200), private :: log_filename
   integer, parameter, private :: log_no = 32
   logical, private :: print_stdout = .false.
   
   abstract interface
      
      subroutine openlog_interface(file_name,stdout)
         implicit none
         character(len=*), intent(in) :: file_name
         logical, intent(in) :: stdout
      end subroutine
      
      subroutine writelog_interface(st,no_comment)
         implicit none
         character(len=*), intent(in) :: st
         logical, optional, intent(in) :: no_comment
      end subroutine
      
      subroutine closelog_interface(fstatus)
         implicit none
         character(len=*), optional, intent(in) :: fstatus
      end subroutine
      
      subroutine abort_interface(st)
         implicit none
         character(len=*), optional, intent(in) :: st
      end subroutine
      
   end interface
   
contains
   
   subroutine openlog_stdout(filename,stdout)
      implicit none
      character(len=*), intent(in) :: filename
      logical, intent(in) :: stdout
      
      log_filename = trim(filename)
      print_stdout = stdout
      
      if (log_filename /= "") then
         open(unit=log_no, file=log_filename, status='unknown')
      end if
      
   end subroutine
   
   subroutine writelog_stdout(st,no_comment)
      use iso_fortran_env, only : output_unit
      implicit none
      character(len=*), intent(in) :: st
      logical, optional, intent(in) :: no_comment
      character(len=10) :: fmt
      
      ! By default, everything is a comment
      fmt = '("#",A)'
      if (present(no_comment)) fmt = '(1x,A)'

      if (log_filename /= "") then
         write(log_no,fmt) trim(st)
      end if

      if (print_stdout) then
         write(*,fmt) trim(st)
         flush(output_unit)
      end if
      
   end subroutine
   
   subroutine closelog_stdout(fstatus)
      implicit none
      character(len=*), optional, intent(in) :: fstatus
      character(len=6) :: file_status
      
      file_status='keep'
      if (present(fstatus)) file_status=trim(fstatus)
      
      if (log_filename /= "") close(unit=log_no, status=trim(file_status))
      
   end subroutine
   
   subroutine abort_nompi(st)
      implicit none
      character(len=*), optional, intent(in) :: st
      
      if (present(st)) call writelog(st)
      call closelog
      stop
   end subroutine
   
end module pnfam_logger
