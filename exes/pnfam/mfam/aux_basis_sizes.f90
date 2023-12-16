!------------------------------------------------------------------------------
! pnfam_matrix/aux_basis_sizes.f90
!
! T. Shafer, UNC Chapel Hill, 2015
!------------------------------------------------------------------------------
program pnfam_matrix_cutoffs

   use logger
   use blockmatrix_type,  only : blockmatrix
   use hfbtho_basis,      only : get_hfbtho_solution
   use hfb_solution_type, only : hfb_solution
   use pnfam_matrix,      only : make_two_qp_basis, input_file_from_cli
   implicit none

   integer, parameter :: dp = kind(1.0d0)

   integer,  parameter :: k_configs(5)  = [0, 1, 0, 1, 2]
   integer,  parameter :: p_configs(5)  = [1, 1,-1,-1,-1]
   real(dp), parameter :: e_configs(13) = [5.0d0, 10.0d0, 15.0d0, 20.0d0, 25.0d0, 30.0d0, &
                                          35.0d0, 40.0d0, 45.0d0, 50.0d0, 55.0d0, 60.0d0, 240.0d0]

   ! HFB solution
   character(len=200) :: hfbfile
   type(hfb_solution) :: hfb_soln

   integer  :: iconfig, icut
   real(dp) :: ecutoff
   integer, allocatable :: twoqp(:,:)

   !----------------------------------------------------------------------------
   !----------------------------------------------------------------------------

   openlog  => openlog_stdout
   writelog => writelog_stdout
   closelog => closelog_stdout
   call openlog('')

   call input_file_from_cli(hfbfile, 'hfb.solution')
   write(*,'("Using HFB solution file: ",1a)') trim(hfbfile)

   call get_hfbtho_solution(fn=trim(hfbfile),         &
      ep=hfb_soln%ep, vp=hfb_soln%vp, up=hfb_soln%up, &
      en=hfb_soln%en, vn=hfb_soln%vn, un=hfb_soln%un)

   do iconfig=1, 5
      write(*,'(/"Config #",1i0,": ",1i0,1a)') iconfig, k_configs(iconfig), merge('+', '-', p_configs(iconfig)>0)

      do icut=1, size(e_configs)
         ecutoff = e_configs(icut)

         call make_two_qp_basis(class='20', k=k_configs(iconfig), parity=p_configs(iconfig), &
            ecut=ecutoff, hfb=hfb_soln, basis_out=twoqp)

         write(*,'("Ecutoff = ",1f5.1," MeV,    NUM_AB = ",1i0)') ecutoff, size(twoqp,2)
      end do
   end do

end program pnfam_matrix_cutoffs
