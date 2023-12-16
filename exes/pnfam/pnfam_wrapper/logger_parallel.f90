!------------------------------------------------------------------------------
! logger_parallel.f90
!
! Parallel logging routines, using MPI I/O.
!
! M.T. Mustonen, UNC Chapel Hill, 2013-15
!------------------------------------------------------------------------------
module logger_mpi
   use iso_fortran_env, only : output_unit, error_unit
   use mpi
   implicit none
   
   integer, private :: file_handle
   integer, private :: comm, rank
   logical, private :: comm_set = .false.
   
contains
   
   ! This is separated to keep the parameter list of openlog_mpi compatible
   ! with the abstract interface openlog_interface
   subroutine set_mpi_file_communicator(c,r)
      implicit none
      integer, intent(in) :: c,r
      
      comm = c
      rank = r
      
   end subroutine
   
   
   subroutine openlog_mpi(filename)
      use logger
      implicit none
      character(len=*), intent(in) :: filename
      
      integer :: ierr, ierr1, resultlength=120
      character(len=120) :: fn
      character(len=mpi_max_error_string) :: msg
      
      ! Only rank (in group) 0 knows the filename, so it must share
      fn = filename
      call MPI_BCAST(fn, 120, MPI_CHARACTER, 0, comm, ierr)
      
      if (fn /= ' ') then
         call MPI_FILE_OPEN(comm, fn, MPI_MODE_WRONLY + MPI_MODE_CREATE + &
            MPI_MODE_SEQUENTIAL, MPI_INFO_NULL, file_handle, ierr)
         if (ierr /= 0) then
            call MPI_ERROR_STRING(ierr, msg, resultlength, ierr1)
            write(error_unit,*) msg
            write(error_unit,*) "Opening log file without sequential mode"
            call MPI_FILE_OPEN(comm, fn, MPI_MODE_WRONLY + MPI_MODE_CREATE, &
            MPI_INFO_NULL, file_handle, ierr)
            if (ierr /= 0) then
               write(error_unit,*) "ERROR: Unable to open the log file for writing!"
               call MPI_ABORT(MPI_COMM_WORLD, -1, ierr)
            end if
         end if
      else
         ! If filename was empty, default write and close commands back to stdout logger
         writelog => writelog_stdout
         closelog => closelog_stdout
         call openlog_stdout(filename)
      end if
      
   end subroutine
   
   
   subroutine closelog_mpi
      implicit none
      integer :: ierr
      
      call MPI_FILE_CLOSE(file_handle, ierr)
      
   end subroutine
   
   
   subroutine writelog_mpi(line)
      implicit none
      character(len=*), intent(in) :: line
      
      integer :: istatus(MPI_STATUS_SIZE), ierr
      character(len=200) :: ln
      
      write(ln,'(2A)') trim(line), char(10)
      
      call MPI_FILE_WRITE_SHARED(file_handle, trim(ln), len_trim(ln), MPI_CHARACTER, &
         istatus, ierr)
      if (ierr /= 0) then
         write(error_unit,*) "ERROR writing to the MPI file!"
         call MPI_ABORT(MPI_COMM_WORLD, -1, ierr)
      end if
      
   end subroutine
   
   
   subroutine abort_mpi
      implicit none
      integer :: ierr
      
      call MPI_ABORT(MPI_COMM_WORLD, -1, ierr)
      
   end subroutine
   
end module logger_mpi
