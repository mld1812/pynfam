!------------------------------------------------------------------------------
!> This module contains all the variables/interfaces used by the pnfam solver and
!> the setup routines to populate/implement them.
!>
!> @authors M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
!> @authors E.M. Ney, UNC Chapel Hill, 2018-
!------------------------------------------------------------------------------
module pnfam_setup
!$ use omp_lib, only : omp_get_max_threads
   use pnfam_logger
   use type_extfield
   use type_blockmatrix
   use pnfam_constants, only : translate_uppercase, version, commit, nthreads
   use pnfam_broyden, only : broyden_history_size, qrpa_alphamix, si,&
      qrpa_bbroyden, nbroy
   use pnfam_interaction, only : print_interaction_params,&
      interaction_name, skip_residual_interaction,&
      require_gauge_invariance, require_self_consistency,&
      force_j2_terms, vpair0, vpair1, vpair_t0, vpair_t1,&
      override_cs0, override_csr, override_cds, override_ct,&
      override_cf, override_cgs, override_cj, override_csdj

   implicit none

   integer, parameter, private :: dp = kind(1d0)

   !===========================================================================
   !                        PNFAM SOLVER VARIABLES
   !        (populate these using external routines for solver to run)
   !===========================================================================
   ! setup_hfb_solution()
   !---------------------------------------------------------------------------
   integer               :: hfb_npr(3)            !> N, Z, A
   real(dp), allocatable :: Ep(:), En(:)          !> QP energy arrays
   type(blockmatrix)     :: Up, Vp, Un, Vn        !> QP wavefunction matrices
   ! (optional, statistical HFB)
   logical               :: ft_active=.false.     !> Is finite temp active
   logical               :: hfb_blo_active=.false.!> Is blocking active in N and/or P
   real(dp), allocatable :: qp_fp(:), qp_fn(:)    !> QP occupations (FT or EFA)
   real(dp)              :: ft_temp               !> Temperature in MeV
   integer               :: hfb_blo_qp(2)         !> Blocked QP indentity: overall index
   integer               :: hfb_blo_ib(2), hfb_blo_is(2) !> Blocked QP indentity: block and state within block
   character(13)         :: hfb_blo_label(2)      !> Blocked QP asymptotic quantum numbers
   integer               :: nshells

   ! setup_extfield()
   !---------------------------------------------------------------------------
   integer :: nxterms
   type(external_field) :: f
   type(external_field), dimension(:), allocatable :: g

   ! setup_hamiltonian()
   !---------------------------------------------------------------------------
   ! - Just provide the interface: input is densities, output is hamiltonian.
   ! - Implement the interface in the setup routine below.
   procedure(calc_dHsp_interface), pointer :: calc_hamiltonian
   abstract interface
      subroutine calc_dHsp_interface(rerho_pn, imrho_pn, rekp, imkp,&
                           rerho_np, imrho_np, rekm, imkm,&
                           reh_pn, imh_pn, redp, imdp,&
                           reh_np, imh_np, redm, imdm)
      use type_blockmatrix, only : blockmatrix
      implicit none
      ! Perturbed generalized density (sp basis)
      type(blockmatrix), intent(inout) :: rerho_pn, imrho_pn, rekp, imkp,&
                                          rerho_np, imrho_np, rekm, imkm
      ! Perturbed hamiltonian (sp basis)
      type(blockmatrix), intent(inout) :: reh_pn, imh_pn, redp, imdp,&
                                          reh_np, imh_np, redm, imdm
      end subroutine calc_dHsp_interface
   end interface

   !===========================================================================
   !                        PNFAM INTERNAL VARIABLES
   !===========================================================================
   ! namelist parameters
   !---------------------------------------------------------------------------
   ! general
   logical            :: print_stdout
   character(len=200) :: fam_output_filename
   integer            :: use_fam_storage
   real(dp)           :: real_eqrpa, imag_eqrpa
   ! external field
   character(len=1)  :: beta_type
   character(len=80) :: operator_name
   integer           :: operator_k
   logical           :: compute_crossterms
   integer           :: two_body_current_mode
   logical           :: two_body_current_usep
   real(dp)          :: two_body_current_lecs(3)
   ! solver
   integer  :: max_iter
   real(dp) :: convergence_epsilon
   ! allow use in setup routines without reading namelist by initizlizing
   real(dp) :: energy_shift_prot=0, energy_shift_neut=0
   real(dp) :: quench_residual_int=1

   namelist /general/ fam_output_filename, print_stdout, use_fam_storage,&
            real_eqrpa, imag_eqrpa
   namelist /ext_field/ beta_type, operator_name, operator_k, compute_crossterms,&
            two_body_current_mode, two_body_current_usep, two_body_current_lecs
   namelist /solver/ max_iter, broyden_history_size, convergence_epsilon,&
            energy_shift_prot, energy_shift_neut, quench_residual_int
   namelist /interaction/ interaction_name, require_self_consistency,   &
            require_gauge_invariance, force_j2_terms,                   &
            vpair0, vpair1, vpair_t0, vpair_t1, override_cs0,           &
            override_csr, override_cds,  override_ct, override_cf,      &
            override_cgs, override_cj, override_csdj

   ! setup_pnfam/ifam variables
   !---------------------------------------------------------------------------
   integer :: ierr, iter_conv, nxy
   real(dp), dimension(:), allocatable :: re_str, im_str

   ! pnfam_output variables
   !---------------------------------------------------------------------------
   character(len=160) :: st
   character(len=80) :: parameter_filename ! input nml
   character(len=200) :: output_txtfile="",output_binfile="",extfield_binfile=""
   integer :: date_values(8)

contains

   !===========================================================================
   ! Primary pnfam setup routines
   !
   ! These routines should not need to be changed, they just call the separate,
   ! setup routines below, and perform some pnfam-specific initializations, like
   ! handling the namelist.
   !===========================================================================

   ! Call all the setup routines
   !----------------------------------------------
   subroutine setup_pnfam
      implicit none

      !$ nthreads = omp_get_max_threads()

      ! Setup the output file procedures
      openlog  => openlog_stdout
      writelog => writelog_stdout
      closelog => closelog_stdout
      abort    => abort_nompi

      ! Open the log for error reporting
      if (fam_output_filename /= "") then
         write(output_binfile,'(a,".bin")') trim(fam_output_filename)
         write(output_txtfile,'(a,".dat")') trim(fam_output_filename)
         write(extfield_binfile,'(a,".tbc")') trim(fam_output_filename)
      else
         print_stdout = .true.
         use_fam_storage = 0
      end if
      call openlog(output_txtfile, print_stdout)

      ! Prepare everything needed for the solver
      call setup_hfb_solution(ierr)
      if(ierr /= 0) then
         call abort(" ERROR reading hfbtho solution.")
      end if
      call setup_hamiltonian(ierr)
      if(ierr /= 0) then
         call abort(" ERROR setting up hamiltonian routine.")
      end if

      ! Check no EFA and FT
      if (hfb_blo_active .and. ft_active) then
         call abort(" ERROR: this code cannot handle T>0 and odd-Z/odd-N simultaneously.")
      end if

      ! Close the log and delete if nothing was written to it
      call closelog('delete')

   end subroutine

   ! Default values for namelist parameters
   !----------------------------------------------
   subroutine init_pnfam_namelist
      implicit none

      ! GENERAL
      fam_output_filename = ''; print_stdout = .true.; use_fam_storage = 0
      real_eqrpa = 0; imag_eqrpa = 0.5_dp
      ! EXT_FIELD
      beta_type = '-'; operator_k = 0; operator_name = 'F'
      compute_crossterms = .false.
      two_body_current_mode = 0
      two_body_current_usep = .false.
      two_body_current_lecs(1) = -3.4
      two_body_current_lecs(2) = 5.4
      two_body_current_lecs(3) = 0
      ! SOLVER
      broyden_history_size = 50; convergence_epsilon = 1.0d-7; max_iter = 200
      energy_shift_prot = 0.0_dp; energy_shift_neut = 0.0_dp
      quench_residual_int = 1.0_dp
      ! INTERACTION
      interaction_name = 'NONE'; force_j2_terms = .false.
      require_gauge_invariance = .false.; require_self_consistency = .true.
      vpair0 = -huge(1.0d0); vpair_t0 = -huge(1.0d0)
      vpair1 = -huge(1.0d0); vpair_t1 = -huge(1.0d0)
      override_cs0 = huge(1.0d0)
      override_csr = huge(1.0d0)
      override_cds = huge(1.0d0)
      override_ct  = huge(1.0d0)
      override_cf  = huge(1.0d0)
      override_cgs = huge(1.0d0)
      override_cj  = huge(1.0d0)
      override_csdj= huge(1.0d0)

   end subroutine init_pnfam_namelist

   ! Read the namelist
   !----------------------------------------------
   subroutine read_pnfam_namelist(fn,use_cli)
      implicit none
      character(len=*), intent(in) :: fn
      logical, intent(in) :: use_cli
      integer, parameter :: fnml = 10

      parameter_filename = fn
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

      open(fnml,file=parameter_filename,iostat=ierr)

      read(fnml,nml=general,iostat=ierr)
      read(fnml,nml=ext_field,iostat=ierr)
      read(fnml,nml=interaction,iostat=ierr)
      read(fnml,nml=solver,iostat=ierr)

      close(fnml)

      if (ierr /= 0) then
         write(*,'(1x,a,1x,i0)') 'ERROR: problem reading pnfam namelist:', ierr
         call abort
      end if

   end subroutine read_pnfam_namelist

   !===========================================================================
   ! Series of small setup functions to populate variables needed by the solver
   !
   ! We separate the setup functions to provide a modicum of flexibility for
   ! the FAM to work with other ground state solvers. These functions can be
   ! replaced, as long as we populate the public variables in this module the
   ! FAM will have what it needs to run the iterative solver. Since HFBTHO is
   ! integerated we simply pull the necessary variables from the hfbtho module here.
   !
   ! It's a little complicated, because the interaction and the external field
   ! are so engrained into hfbtho and the pnfam namelist.
   !===========================================================================

   ! Get the U,V,E and create blockmatrices
   !----------------------------------------------
   subroutine setup_hfb_solution(ierr)
      use pnfam_constants, only : IT_NEUTRON, IT_PROTON
      use hfb_solution, only : get_raw_hfb_solution
      use hfb_solution, only : hnpr=>hfb_npr, hnb, hdb, n00
      use hfb_solution, only : Epe=>Ep,Ene=>En
      use hfb_solution, only : Vpe=>Vp,Vne=>Vn,Upe=>Up,Une=>Un
      ! Statistical HFB
      use hfb_solution, only : qp_t_reverse
      use hfb_solution, only : hfb_ft_active=>ft_active, hfb_ft_temp=>ft_temp,&
                               hfb_ft_fp=>ft_fp, hfb_ft_fn=>ft_fn,&
                               tb, numax, blok1k2d, keyblo,&
                               bloblo, blo123, pwi_qp_n, pwi_qp_p

      implicit none
      integer :: i, dmat, dqp, it
      integer, intent(out) :: ierr
      ierr = 0

      ! Load the hfb solution from hfbtho program
      call get_raw_hfb_solution(ierr)
      if (ierr /= 0) return
      hfb_npr = hnpr
      nshells = n00

      ! Setup block matrix structure
      call init_blockmatrix_type(hnb, hdb, ierr)
      if (ierr /= 0) return
      dmat = sum(hdb*hdb)
      dqp = sum(hdb)

      ! Setup quasiparticle solution
      Ep = Epe + energy_shift_prot
      En = Ene + energy_shift_neut

      call allocate_blockmatrix(Up, dmat); Up%elem = Upe
      call allocate_blockmatrix(Vp, dmat); Vp%elem = Vpe
      call allocate_blockmatrix(Un, dmat); Un%elem = Une
      call allocate_blockmatrix(Vn, dmat); Vn%elem = Vne

      ! Block structure: U is block diagonal (square blocks)
      Un%ir2c = [ (i, i=1, hnb) ]; Un%ic2r = Un%ir2c
      Un%ir2m(1) = 1
      do i=2,hnb
        Un%ir2m(i) = Un%ir2m(i-1) + hdb(i-1)**2
      end do
      Un%ic2m = Un%ir2m
      call copy_block_structure(Un, Up)

      ! Block structure: V connects K to -K (square blocks)
      Vn%ir2c(1:hnb/2) = [ (i+hnb/2, i=1, hnb/2) ]
      Vn%ir2c(hnb/2+1:hnb) = [ (i, i=1, hnb/2) ]
      Vn%ic2r = Vn%ir2c
      Vn%ir2m = Un%ir2m
      Vn%ic2m = Vn%ir2m(Vn%ic2r(:))
      call copy_block_structure(Vn, Vp)

      ! Finite temperature quantities
      ft_active = .false.
      ft_temp   = 0.0_dp
      if (hfb_ft_active) then
         ft_active = .true.
         ft_temp = hfb_ft_temp
         qp_fp = hfb_ft_fp
         qp_fn = hfb_ft_fn
      end if

      ! EFA-Blocking Quantities
      hfb_blo_active = .false.
      hfb_blo_label(:) = ""
      hfb_blo_qp(:) = 0; hfb_blo_ib(:) = 0;  hfb_blo_is(:) = 0
      do it=1,2
         if (keyblo(it) /= 0) then
            hfb_blo_active = .true.
            hfb_blo_label(it) = tb(numax(blok1k2d(it),it))
            if (it==IT_NEUTRON) then
               hfb_blo_qp(it) = pwi_qp_n(blok1k2d(it))
            else if (it==IT_PROTON) then
               hfb_blo_qp(it) = pwi_qp_p(blok1k2d(it))
            end if
            hfb_blo_ib(it) = bloblo(keyblo(it),it)
            hfb_blo_is(it) = blo123(keyblo(it),it)
         end if
      end do
      if (hfb_blo_active) then
         if (allocated(qp_fn)) deallocate(qp_fn, qp_fp)
         allocate(qp_fn(dqp), qp_fp(dqp))
         qp_fn = 0; qp_fp = 0
      end if
      if (hfb_blo_qp(IT_NEUTRON) /= 0) then
         qp_fn(hfb_blo_qp(IT_NEUTRON)) = 0.5_dp
         qp_fn(qp_t_reverse(hfb_blo_qp(IT_NEUTRON))) = 0.5_dp
      end if
      if (hfb_blo_qp(IT_PROTON) /= 0) then
         qp_fp(hfb_blo_qp(IT_PROTON)) = 0.5_dp
         qp_fp(qp_t_reverse(hfb_blo_qp(IT_PROTON))) = 0.5_dp
      end if

   end subroutine setup_hfb_solution



   ! Initialize the interaction couplings
   ! (depends on namelist, hfb couplings, and coord space density)
   !----------------------------------------------
   subroutine setup_hamiltonian(ierr)
      use pnfam_interaction, only : init_interaction
      use pnfam_hamiltonian, only : calc_dHsp=>calc_hamiltonian
      implicit none
      integer, intent(out) :: ierr
      ierr = 0

      call init_interaction
      if (skip_residual_interaction) quench_residual_int = 0
      calc_hamiltonian => calc_dHsp

   end subroutine setup_hamiltonian

end module pnfam_setup
