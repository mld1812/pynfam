!------------------------------------------------------------------------------
!> Charge-changing FAM code to compute pnQRPA strength functions at a single
!> complex energy point, using a solution from HFBTHO.
!>
!> @authors M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
!> @authors E.M. Ney, UNC Chapel Hill, 2018-
!------------------------------------------------------------------------------
program pnfam_prog
   use pnfam_solver
   implicit none
   character(len=80), parameter :: file_namelist = "pnfam_NAMELIST.dat"

   call init_pnfam_namelist
   call read_pnfam_namelist(file_namelist, .true.)

   call pnfam_solve

end program pnfam_prog
