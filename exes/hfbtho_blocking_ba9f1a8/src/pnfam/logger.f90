!------------------------------------------------------------------------------
! logger.f90
!
! Logging routines for pnFAM.  These should be called via the pointers defined,
! so that they can be overridden with more sophisticated routines with MPI I/O
! in the case of the parallel code.
!
! M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
!------------------------------------------------------------------------------
module logger
   use iso_fortran_env
   implicit none
   
   ! By default, we call the standard output routines
   procedure(openlog_interface), pointer :: openlog
   procedure(writelog_interface), pointer :: writelog
   procedure(closelog_interface), pointer :: closelog
   procedure(abort_interface), pointer :: abort
   
   abstract interface
      
      subroutine openlog_interface(file_name)
         implicit none
         character(len=*), intent(in) :: file_name
      end subroutine
      
      subroutine writelog_interface(st)
         implicit none
         character(len=*), intent(in) :: st
      end subroutine
      
      subroutine closelog_interface
         implicit none
      end subroutine
      
      subroutine abort_interface
         implicit none
      end subroutine
      
   end interface
   
contains
   
   subroutine openlog_stdout(file_name)
      implicit none
      character(len=*), intent(in) :: file_name
      
      ! No need to open stdout
      ! The filename is just a dummy variable
      
   end subroutine
   
   subroutine writelog_stdout(st)
      use iso_fortran_env, only : output_unit
      implicit none
      character(len=*), intent(in) :: st
      
      write(*,'(A)') trim(st)
      flush(output_unit)
      
   end subroutine
   
   subroutine closelog_stdout
      implicit none
      
      ! No need to close stdout
      
   end subroutine
   
   subroutine abort_nompi
      implicit none
      stop -1
   end subroutine
   
end module logger
