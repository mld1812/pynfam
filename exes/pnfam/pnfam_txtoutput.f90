!------------------------------------------------------------------------------
!> This module contains routines which write the text output file for a pnfam
!> solution.
!>
!> @authors M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
!> @authors E.M. Ney, UNC Chapel Hill, 2018-
!------------------------------------------------------------------------------
module pnfam_txtoutput
   use pnfam_setup
   use pnfam_constants
   private

   integer, parameter, private :: dp = kind(1d0)

   integer, parameter :: w=51    ! header width
   integer, parameter :: fout=11 ! file handle

   integer :: istr

   public :: log_header, log_result, txtr, txtc, txtl, print_main_header

contains
   !---------------------------------------------------------------------------
   ! Textfile output header
   !---------------------------------------------------------------------------
   subroutine log_header
      implicit none

      call print_solver_details

      call print_interaction_params

      !call print_hfb_details

      call print_statfam_details

   end subroutine

   !----------------------------------------------------------------------------
   ! Primary FAM code header
   !----------------------------------------------------------------------------
   subroutine print_main_header
      use pnfam_constants, only : element_symbols
      implicit none
      integer :: zf, nf, l
      character(3) :: sn

      ! Nucleus
      l = len_trim(operator_name)
      zf=0; nf=0
      if (.not. any(hfb_npr == 0)) then
      if (beta_type == '-') then
         zf = hfb_npr(2) + 1
         nf = hfb_npr(1) - 1
      else if (beta_type == '+') then
         zf = hfb_npr(2) - 1
         nf = hfb_npr(1) + 1
      end if
      end if

      call writelog("")
      write(st,'(1x,A)') repeat("=",w); call writelog(st)

      write(st,'(A)')  "pnFAM";                                   call writelog(txtc(st,w))
      write(st,'(A)')  "Charge-changing Finite Amplitude Method"; call writelog(txtc(st,w))
      write(st,'(2A)') "Version: ", version;                      call writelog(txtc(st,w))
      if (nthreads > 0) then
         write(st,'(A,1I0,A)') "Compiled with OpenMP (threads=",nthreads,")"; call writelog(txtc(st,w))
      end if
      call writelog("")
      write(st,'(2A)') "Commit: ", trim(commit);                 call writelog(txtc(st,w))
      call date_and_time(values=date_values)
      write(st,'("Run date: ",1i0,"/",1i0,"/",1i0,1x,1i0,":",1i2.2,":",1i2.2)') date_values(2), &
         date_values(3), date_values(1), date_values(5:7); call writelog(txtc(st,w))

      write(st,'(1x,A)') repeat("=",w); call writelog(st)

      write(st,'(2A,A1," for K = ",I0)') "Operator: ", trim(operator_name), beta_type, operator_k
      call writelog(txtc(st,w))
      sn = " + "; if(imag_eqrpa<0) sn = " - "
      write(st,'(A,"(",F8.4,A,F8.4,"i)")') "Energy: ", real_eqrpa, sn, abs(imag_eqrpa)
      call writelog(txtc(st,w))
      write(st,'(A,I0,A,A,I0,A,I0,A)') "Parent nucleus:   ", hfb_npr(3), element_symbols(hfb_npr(2)), &
          "(N=",hfb_npr(1), ", Z=", hfb_npr(2), ")"; call writelog(txtc(st,w))
      write(st,'(A,I0,A,A,I0,A,I0,A)') "Daughter nucleus: ", hfb_npr(3), element_symbols(zf), &
          "(N=", nf, ", Z=", zf, ")"; call writelog(txtc(st,w))

      write(st,'(1x,A)') repeat("=",w); call writelog(st)
      call writelog("")

   end subroutine print_main_header

   !----------------------------------------------------------------------------
   ! FAM solver details
   !----------------------------------------------------------------------------
   subroutine print_solver_details
      implicit none
      real(dp) :: dp_tol=1e-10
      character(len=100) :: stx="None"

      write(st,'(1x,A)') repeat("-",w); call writelog(st)
      write(st,'(1x,A)') "pnFAM solver details"; call writelog(txtc(st,w))
      write(st,'(1x,A)') repeat("-",w); call writelog(st)

      !write(st,'(3A)') "HFB input data file name:  '", trim(hfb_input_filename), "'"; call writelog(txtp(st,p))
      write(st,'(1x,A,A)') &
         "FAM input parameter file name:   ", "'"//trim(parameter_filename)//"'";  call writelog(st)
      write(st,'(1x,A,A)') &
         "FAM output file name base:       ", "'"//trim(fam_output_filename)//"'"; call writelog(st)

      call writelog("")
      write(st,'(1x,2A,A1,2x,"K=",1x,I0)') &
         "Operator:                        ", trim(operator_name), beta_type, operator_k ; call writelog(st)
      write(st,'(1x,1A,I0,1x,"(",I0,")")') &
         "Two-body current mode (use_p):   ", two_body_current_mode, merge(1,0,two_body_current_usep) ; call writelog(st)
      if (two_body_current_mode /= 0) then
         write(st,'(1x,A,sp,3(F8.4,1x))') &
         "  LECs (c3, c4, cd):             ", two_body_current_lecs ; call writelog(st)
      end if
      write(st,'(1x,2A)') &
         "Compute cross-terms:             ", merge('Yes', 'No ', compute_crossterms) ; call writelog(st)
      if (compute_crossterms) then
         if (nxterms > 0) stx = txt_clist(g%label)
         write(st,'(1x,2A)') &
         "  Cross-terms computed:          ", stx ; call writelog(st)
       end if
      write(st,'(1x,A,F8.4," MeV")') &
         "Re(EQRPA):                       ", real_eqrpa;       call writelog(st)
      write(st,'(1x,A,F8.4," MeV")') &
         "Im(EQRPA):                       ", imag_eqrpa;         call writelog(st)

      call writelog("")
      write(st,'(1x,A,I0)') &
         "Number shells:                   ", nshells;            call writelog(st)
      write(st,'(1x,A,I0)') &
         "Basis size:                      ", sum(db);            call writelog(st)
      write(st,'(1x,A,I0)') &
         "Number matrix blocks:            ", nb;                 call writelog(st)
      write(st,'(1x,A,I0)') &
         "Non-trivial HFB matrix elements: ", sum(db*db);         call writelog(st)
      write(st,'(1x,A,I0)') &
         "Non-trivial FAM matrix elements: ", nxy;                call writelog(st)

      call writelog("")
      write(st,'(1x,A,I0)') &
         "Maximum iterations:              ", max_iter;              call writelog(st)
      write(st,'(1x,A,I0)') &
         "Broyden history size:            ", broyden_history_size;  call writelog(st)
      write(st,'(1x,A,f4.2)') &
         "Broyden mixing factor:           ", qrpa_alphamix;         call writelog(st)
      write(st,'(1x,A,ES8.1)') &
         "Convergence limit:               ", convergence_epsilon;   call writelog(st)
      if (abs(quench_residual_int-1.0_dp)>dp_tol) then
         write(st,'(1x,1A,1F6.3," (x dH)")') &
         "Interaction Quenching:           ", quench_residual_int;   call writelog(st)
      end if
      if (abs(energy_shift_prot)>dp_tol .or. abs(energy_shift_neut)>dp_tol) then
         write(st,'(1x,1A,2F6.3," MeV")') &
         "QP energy shift:                 ", energy_shift_prot, energy_shift_neut; call writelog(st)
      end if
      call writelog("")

   end subroutine

   !----------------------------------------------------------------------------
   ! Statistical FAM details
   !----------------------------------------------------------------------------
   subroutine print_statfam_details
      use pnfam_constants,    only : IT_NEUTRON, IT_PROTON
      implicit none

      write(st,'(1x,A)') repeat("-",w); call writelog(st)
      write(st,'(1x,A)') "Statistical pnFAM details"; call writelog(txtc(st,w))
      write(st,'(1x,A)') repeat("-",w); call writelog(st)

      write(st,'(1x,2A)') &
         "Finite-temperature active: ", trim(merge("Yes", "No ", ft_active)) ; call writelog(st)

      if (ft_active) then
         write(st,'(3x,A,1F8.4,A, " MeV")') &
         "Temperature .............: ", ft_temp ; call writelog(st)
         write(st,'(3x,A,1F8.4,1x,"(Ep=",1F0.4,"MeV)")') &
         "Approx. max. f_p ........: ", maxval(qp_fp), Ep(maxloc(qp_fp)) ; call writelog(st)
         write(st,'(3x,A,1F8.4,1x,"(En=",1F0.4,"MeV)")') &
         "Approx. max. f_n ........: ", maxval(qp_fn), En(maxloc(qp_fn)) ; call writelog(st)
      end if

      write(st,'(1x,2A)') &
         "Odd-nucleus EFA active ..: ", trim(merge("Yes", "No ", hfb_blo_active)) ; call writelog(st)

      if (hfb_blo_active) then
         if (hfb_blo_qp(IT_NEUTRON) /= 0) then
            write(st,'(3x,"Neutron blocking ..... : On")')
            call writelog(st)
            write(st,'(5x,"Blocked level ...... : ", a)') &
                hfb_blo_label(IT_NEUTRON)
            call writelog(st)
            write(st,'(5x,"QP energy .......... : ", 1f6.4, " MeV")') En(hfb_blo_qp(IT_NEUTRON))
            call writelog(st)
            write(st,'(5x,"QP index ........... : ", 1i0, " (block = ", 1i0, ", state = ", 1i0, ")")') &
               hfb_blo_qp(IT_NEUTRON), hfb_blo_ib(IT_NEUTRON), hfb_blo_is(IT_NEUTRON)
            call writelog(st)
         else
            write(st,'(3x,"Neutron blocking ... : Off")')
            call writelog(st)
         end if
         if (hfb_blo_qp(IT_PROTON) /= 0) then
            write(st,'(3x,"Proton blocking ...... : On")')
            call writelog(st)
            write(st,'(5x,"Blocked level ...... : ", a)') &
               hfb_blo_label(IT_PROTON)
            call writelog(st)
            write(st,'(5x,"QP energy .......... : ", 1f6.4, " MeV")') Ep(hfb_blo_qp(IT_PROTON))
            call writelog(st)
            write(st,'(5x,"QP index ........... : ", 1i0, " (block = ", 1i0, ", state = ", 1i0, ")")') &
               hfb_blo_qp(IT_PROTON), hfb_blo_ib(IT_PROTON), hfb_blo_is(IT_PROTON)
            call writelog(st)
         else
            write(st,'(3x,"Proton blocking .... : Off")')
            call writelog(st)
         end if
      end if
      call writelog("")

   end subroutine print_statfam_details

   !----------------------------------------------------------------------------
   ! Write the strength with all the digits in the output footer
   !----------------------------------------------------------------------------
   subroutine log_result
      implicit none
      integer :: ixterms
      character(len=1024) :: sto

      ! Note: Maximum xterm label is 7 chars (3,"x",3)
      call writelog("")
      write(sto,'(1x,A)') "Result (Energy [MeV], strength and cross-terms [MeV^-1]):"
         call writelog(sto)
      write(sto,'(1x,A)') repeat("-",78)
         call writelog(sto)
      write(sto,'(1x,10x,A,A)') txtr("Real",34), txtr("Imag",34)
         call writelog(sto,.true.)
      write(sto,'(1x,A,2ES34.19)') txtl("Energy",10), real_eqrpa, imag_eqrpa
         call writelog(sto,.true.)
      write(sto,'(1x,A,2ES34.19)') txtl("Strength",10), re_str(1), im_str(1)
         call writelog(sto,.true.)
      do ixterms=2, size(re_str)
         write(sto,'(1x,A,2ES34.19)') txtl(trim(g(ixterms-1)%label),10), re_str(ixterms), im_str(ixterms)
         call writelog(sto,.true.)
      end do
      call writelog("",.true.)

   end subroutine

end module pnfam_txtoutput
