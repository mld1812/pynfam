!------------------------------------------------------------------------------
! betadecay.f90
!
! This module takes in strength functions and cross-terms on a complex contour
! and integrates them with the lepton phase space to obtain a beta-decay rate.
! 
! M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
!------------------------------------------------------------------------------
module betadecay
   use iso_fortran_env, only : output_unit
   use constants,       only : gv, hbar_mec, r0, version
   use logger

   implicit none
   
   integer,  parameter :: dp = kind(1.0d0)
   
   ! Read-in/HFB-based beta-decay quantities
   integer  :: zinum, zfnum, anum
   real(dp) :: qvalue, gsenergy, eqrpamax, ga_input
   
   ! Calculated beta-decay quantities
   real(dp) :: lambda, r
   
   ! Contour integral source files
   character(len=250) :: srcfile_k0_f, srcfile_k0_gt, srcfile_k1_gt, srcfile_k0_rs0, &
      srcfile_k0_ps0, srcfile_k0_rs2, srcfile_k1_rs2, srcfile_k2_rs2, srcfile_k0_r, &
      srcfile_k1_r, srcfile_k0_p, srcfile_k1_p, srcfile_k0_rs1, srcfile_k1_rs1
   
   ! Phase space parameterization
   ! Polynomial constants are coeffs(gn=1:6,order=1:x) where x=ratint_pts or psi_order
   character(len=10) :: quadrature
   logical  :: psi_use_polyfit
   integer  :: psi_npts, psi_ncoeffs, rat_order_n, rat_order_d
   integer  :: psi_gl_npts
   real(dp) :: psi_si
   real(dp), allocatable :: coeffs(:,:)
   
   ! Computation input and output files
   ! Defaults are beta.in (IN), beta.out (OUT)
   character(len=250) :: infile, outfile
   
   ! Contour storage
   ! The contour data has a real parameterization, complex eigenvalues, and
   ! an array complex results (depending on cross-terms).
   type contour_data
      integer :: k
      character(len=80) :: label
      character(len=80), dimension(:),   allocatable :: op
      complex(dp),       dimension(:),   allocatable :: dz
      complex(dp),       dimension(:),   allocatable :: eqrpa
      complex(dp),       dimension(:,:), allocatable :: str
      real(dp),          dimension(:),   allocatable :: t
      real(dp),          dimension(:),   allocatable :: gl_wts
   end type contour_data
   
   ! Actual storage of contour data
   type(contour_data), target :: data_k0_f, data_k0_gt, data_k1_gt,           &
      data_k0_rs0, data_k0_ps0, data_k0_r, data_k1_r, data_k0_p, data_k1_p,   &
      data_k0_rs1, data_k1_rs1, data_k0_rs2, data_k1_rs2, data_k2_rs2

   ! Decay rates
   real(dp), target :: rate_k0_f, rate_k0_gt, rate_k1_gt, rate_k0_ff, & 
      rate_k1_ff, rate_k2_ff, rate_j0_ff, rate_j1_ff, rate_j2_ff,     &
      rate_j0_k0_ff, rate_j1_k0_ff, rate_j1_k1_ff, rate_j2_k0_ff,     &
      rate_j2_k1_ff, rate_j2_k2_ff


contains
   
   subroutine betadecay_main_program
      implicit none
   
      ! Initialization
      openlog => openlog_stdout; writelog => writelog_stdout; closelog => closelog_stdout
      call openlog("")
      call show_header
      call read_namelists_and_init

      ! Calculate some basic properties from input
      ! N.b. lambda = gA/gV when gA > 0, cf. Behrens & Buhring.  The essential point is
      ! that lambda should be positive.
      lambda = ga_input/gv
   
      if (lambda < 0) then
         write(*,'(a,1f0.2,a)') 'ERROR: lambda = gA/gV should be a positive parameter&
            & (lambda = ', lambda, ')'
         stop
      end if
   
      ! Nuclear radius = r0 A^(1/3) [fm] * (me c^2) [MeV] / (hbar c) [MeV fm]
      r = r0*anum**(1.0_dp/3.0_dp) / hbar_mec
   
      ! Show calculation details
      write(*,*)
      call show_parameters
   
      ! Read in the contour data
      write(*,'(/,a,/,a)') 'INPUT FILES:', repeat('-',30)
      call load_contour_data(file=srcfile_k0_f,   storage=data_k0_f,   label='K=0 Fermi')
      call load_contour_data(file=srcfile_k0_gt,  storage=data_k0_gt,  label='K=0 Gamow-Teller')
      call load_contour_data(file=srcfile_k1_gt,  storage=data_k1_gt,  label='K=1 Gamow-Teller')
      write(*,*) 
      call load_contour_data(file=srcfile_k0_rs0, storage=data_k0_rs0, label='K=0 [RS]0')
      call load_contour_data(file=srcfile_k0_ps0, storage=data_k0_ps0, label='K=0 [PS]0')
      write(*,*) 
      call load_contour_data(file=srcfile_k0_rs1, storage=data_k0_rs1, label='K=0 [RS]1')
      call load_contour_data(file=srcfile_k1_rs1, storage=data_k1_rs1, label='K=1 [RS]1')
      call load_contour_data(file=srcfile_k0_r,   storage=data_k0_r,   label='K=0 R')
      call load_contour_data(file=srcfile_k1_r,   storage=data_k1_r,   label='K=1 R')
      call load_contour_data(file=srcfile_k0_p,   storage=data_k0_p,   label='K=0 P')
      call load_contour_data(file=srcfile_k1_p,   storage=data_k1_p,   label='K=1 P')
      write(*,*) 
      call load_contour_data(file=srcfile_k0_rs2, storage=data_k0_rs2, label='K=0 [RS]2')
      call load_contour_data(file=srcfile_k1_rs2, storage=data_k1_rs2, label='K=1 [RS]2')
      call load_contour_data(file=srcfile_k2_rs2, storage=data_k2_rs2, label='K=2 [RS]2')   
   
      ! Compute the allowed rate
      ! NOTE -- these rates already include the appropriate factors of 2 to account
      ! for being matrix elements in the instrinsic nuclear frame
      rate_k0_f = 0;  rate_k0_gt = 0;  rate_k1_gt = 0
      write(*,'(2/,a)') 'ALLOWED DECAY RATE:'
      call compute_allowed_rate
   
      ! Compute the forbidden rate
      ! NOTE -- these rates already include the appropriate factors of 2 to account
      ! for being matrix elements in the instrinsic nuclear frame
      rate_k0_ff = 0;     rate_k1_ff = 0;     rate_k2_ff = 0
      rate_j0_ff = 0;     rate_j1_ff = 0;     rate_j2_ff = 0   
      rate_j0_k0_ff = 0;  rate_j1_k0_ff = 0;  rate_j1_k1_ff = 0
      rate_j2_k0_ff = 0;  rate_j2_k1_ff = 0;  rate_j2_k2_ff = 0
      write(*,'(2/,a)') 'FORBIDDEN DECAY RATE:'
      call compute_forbidden_rate
   
      ! Write the results
      call show_results(store=.true.)
   
   end subroutine betadecay_main_program
   
   !----------------------------------------------------------------------------
   ! Initial reading of parameters from a set of namelists.
   !
   ! Assigns values to the globals:
   !  zinum 
   !  zfnum 
   !  anum
   !  qvalue
   !  gsenergy
   !  eqrpamax
   !  ga_input
   !  infile
   !  outfile
   !  srcfile_...    
   !
   ! Also initializes and allocates:
   !  cpoly
   !----------------------------------------------------------------------------
   subroutine read_namelists_and_init
      use constants, only : gA, translate_uppercase
      use blockmatrix_type
      use pnfam_common, only : setup_hfb_solution, get_hfbtho_betadecay_properties,&
                               Ep,En,Vp,Vn,Up,Un, hfb_npr
      implicit none
      
      integer  :: ierr
      character(len=250) :: namelist_file
      
      ! Namelist "beta_parameters"
      integer  :: z_init, a_init
      real(dp) :: q_value, gs_energy, ga_effective
      character(len=250) :: hfbtho_solution_file, output_file
      
      ! Namelist "phase_space"
      character(len=10) :: quadrature_type
      logical  :: use_polyfit
      real(dp) :: polynomial_fit_si
      integer  :: polynomial_order, polynomial_fit_points
      integer  :: ratint_psi_gl_points, ratint_points
      
      ! Namelist "input_files"
      character(len=250) :: k0_f, k0_gt, k1_gt, k0_rs0, k0_ps0, k0_rs1, k1_rs1, &
         k0_r, k1_r, k0_p, k1_p, k0_rs2, k1_rs2, k2_rs2
   
      namelist /beta_parameters/ z_init, a_init, q_value, gs_energy, &
         ga_effective, hfbtho_solution_file, output_file
      namelist /phase_space/ quadrature_type, use_polyfit, polynomial_order, &
         polynomial_fit_points, polynomial_fit_si, ratint_psi_gl_points, ratint_points
      namelist /input_files/ k0_f, k0_gt, k1_gt, k0_rs0, k0_ps0, k0_rs1, k1_rs1, &
         k0_r, k1_r, k0_p, k1_p, k0_rs2, k1_rs2, k2_rs2
      
      ! Separate storage in case we read from the HFB file
      integer  :: k_hfb, pi_hfb
      real(dp) :: q_hfb, egs_hfb
      
      
      ! Default values
      z_init = -1;  a_init = -1;  q_value = -1;  gs_energy = -1;  ga_effective = gA
      hfbtho_solution_file = "";  output_file = "beta.out"
      
      quadrature_type=""; use_polyfit = .false.; polynomial_order = 5
      polynomial_fit_points = 200;  polynomial_fit_si = 1d-5
      ratint_psi_gl_points = 75 ; ratint_points = 20
      
      ! Setting the filenames to empty strings allows the use of len_trim later
      k0_f   = "";  k0_gt  = "";  k1_gt = "";  k0_rs0 = "";  k0_ps0 = "";  k0_rs1 = ""
      k1_rs1 = "";  k0_r   = "";  k1_r  = "";  k0_p   = "";  k1_p   = "";  k0_rs2 = ""
      k1_rs2 = "";  k2_rs2 = ""
      

      ! Read the input file --- default name is 'beta.in'
      if (command_argument_count() < 1) then
         namelist_file = 'beta.in'
      else
         call get_command_argument(1, namelist_file)
      end if
      
      open(1, file=trim(namelist_file), status='old', iostat=ierr)
      
      if (ierr /= 0) then
         write(*,'(a,1i0,3a)') 'ERROR in sub.read_parameters: ierr = ', ierr, &
            ' while opening the input file "', trim(namelist_file), '".'
         stop
      end if
      
      read(1, nml=beta_parameters)
      read(1, nml=phase_space)
      read(1, nml=input_files)

      close(1)
      
      
      ! If an HFB solution is given, re-set the parameters
      if (len_trim(hfbtho_solution_file) > 0) then
         write(*,'(/,3a)') 'Reading input parameters from HFB solution "', &
            trim(hfbtho_solution_file), '".'
         
         ! Reset the important quantities
         q_hfb = -1;  egs_hfb = -1
         
         ! Now read and store the solution in order to not duplicate code
         call setup_hfb_solution(ierr)
         if (ierr /= 0) then
            call writelog(" Error reading hfb solution.")
            call closelog
            call abort
         end if
         
         ! Beta decay
         write(*,*)
         call get_hfbtho_betadecay_properties(ep=ep, en=en, vp=vp, vn=vn, &
            q=q_hfb, e_gs=egs_hfb, k_gs=k_hfb, p_gs=pi_hfb, lpr=.true.)
         
         ! If the namelist values are less than 0, re-set them
         if (z_init < 0) then
            z_init = hfb_npr(2)
         else
            write(*,'(a)') 'WARNING in sub.read_parameters: did not re-set Z_init from HFB'
         end if
         
         if (a_init < 0) then
            a_init = hfb_npr(3)
         else
            write(*,'(a)') 'WARNING in sub.read_parameters: did not re-set A_init from HFB'
         end if
         if (q_value < 0) then
            q_value = q_hfb
         else
            write(*,'(a)') 'WARNING in sub.read_parameters: did not re-set Q_value from HFB'
         end if
         if (gs_energy < 0) then
            gs_energy = egs_hfb
         else
            write(*,'(a)') 'WARNING in sub.read_parameters: did not re-set GS_energy from HFB'
         end if
      end if
      
      
      ! Basic tests of input parameters for consistency
      ierr = 0
      if (z_init < 0) then
         write(*,'(a,1i0,a)') 'ERROR in sub.read_parameters: z_init < 1 (z_init = ', z_init,').'
         ierr = ierr + 1
      end if
      if (a_init < 1) then
         write(*,'(a,1i0,a)') 'ERROR in sub.read_parameters: a_init < 1 (a_init = ', a_init,').'
         ierr = ierr + 1
      end if
      if (a_init < z_init) then
         write(*,'(a,1i0,a,1i0,a)') 'ERROR in sub.read_parameters: a_init < z_init (z_init = ', &
            z_init, ', a_init = ', a_init,').'
            ierr = ierr + 1
      end if
      if (q_value < 0.0_dp) then
         write(*,'(a,1i0,a)') 'ERROR in sub.read_parameters: Q < 0.'
         ierr = ierr + 1
      end if
      if (gs_energy < 0.0_dp) then
         write(*,'(a,1i0,a)') 'ERROR in sub.read_parameters: E_gs < 0.'
         ierr = ierr + 1
      end if
      
      if (ierr /= 0) stop
      
      
      ! INIT
      !-------------------------------------------------------------------------
      ! Assign the computations parameters to global variables
      
      ! GLOBAL beta decay constants
      zinum    = z_init
      zfnum    = z_init + 1
      anum     = a_init
      qvalue   = q_value
      gsenergy = gs_energy
      eqrpamax = gs_energy+q_value
      ga_input = ga_effective
      
      ! GLOBAL file names
      infile  = trim(namelist_file)
      outfile = trim(output_file)

      ! GLOBAL phase space parameters
      psi_use_polyfit = use_polyfit
      call translate_uppercase(quadrature_type)
      quadrature = quadrature_type
      if (psi_use_polyfit) then
         psi_ncoeffs = polynomial_order
         psi_npts    = polynomial_fit_points
         psi_si      = polynomial_fit_si    ! Tolerance on psi iterative simpsons
         ! These are unused in poly-fit mode
         psi_gl_npts = 0
         rat_order_n = 0
         rat_order_d = 0
      else
         psi_ncoeffs = ratint_points        ! number Thiele coeffs = number points
         psi_npts    = ratint_points
         psi_gl_npts = ratint_psi_gl_points ! Number gauleg pts for psi quadrature
         ! For printout only:
         ! Numerator and denominator polynomial order from number of rat-int points
         rat_order_n = ratint_points/2
         if (rat_order_n*2 /= ratint_points) then
            rat_order_d = rat_order_n
         else
            rat_order_d = rat_order_n - 1
         end if
         ! This is unused in rat-int mode
         psi_si      = 0.0
      end if
      
      ! GLOBAL - Allocate the polynomial fitting storage
      allocate(coeffs(1:6,1:psi_ncoeffs))
      coeffs(:,:) = 0
      
      ! GLOBAL contour integral file names
      srcfile_k0_f   = adjustl(k0_f)
      srcfile_k0_gt  = adjustl(k0_gt)
      srcfile_k1_gt  = adjustl(k1_gt)
      srcfile_k0_rs0 = adjustl(k0_rs0)
      srcfile_k0_ps0 = adjustl(k0_ps0)
      srcfile_k0_rs1 = adjustl(k0_rs1)
      srcfile_k1_rs1 = adjustl(k1_rs1)
      srcfile_k0_r   = adjustl(k0_r)
      srcfile_k1_r   = adjustl(k1_r)
      srcfile_k0_p   = adjustl(k0_p)
      srcfile_k1_p   = adjustl(k1_p)
      srcfile_k0_rs2 = adjustl(k0_rs2)
      srcfile_k1_rs2 = adjustl(k1_rs2)
      srcfile_k2_rs2 = adjustl(k2_rs2)
   end subroutine read_namelists_and_init
   
   !----------------------------------------------------------------------------
   ! Display parameters of the calculation
   !
   ! Uses the globals:
   !  zinum, zfnum, anum, gsenergy, qvalue, eqrpamax
   !  psi_ncoeffs, psi_npts, psi_si, rat_order_n, rat_order_d, psi_gl_npts
   !  ga_input, lambda, r
   !----------------------------------------------------------------------------
   subroutine show_parameters
      use constants, only : alpha, mec2, mn, esy => element_symbols
      
      implicit none
      write(*,'(a,/,a)') 'OVERVIEW PARAMETERS:', repeat('-', 30)
      write(*,'(a,1x,3a)')    'Configuration file ........ :', '"', trim(infile), '"'
      write(*,'(a,1x,3a)')    'Output storage file ....... :', '"', trim(outfile), '"'
      write(*,'(a,1x,1i3,1x,a,"(Z = ",1i3,", N = ",1i3,")")') &
                              'Initial nucleus ........... :', anum, esy(zinum), zinum, anum-zinum
      write(*,'(a,1x,1i3,1x,a,"(Z = ",1i3,", N = ",1i3,")")') &
                              'Final nucleus ............. :', anum, esy(zfnum), zfnum, anum-zfnum

      write(*,*) 

      select case (quadrature)
         case("GAUSS")
            write(*,'(a,1x,a)')   'Contour int. method ....... :', 'Gauss-Legendre'
         case("SIMPSON")
            write(*,'(a,1x,a)')   'Contour int. method ....... :', 'Simpsons 3/8'
         case("TRAP")
            write(*,'(a,1x,a)')   'Contour int. method ....... :', 'Trapezoidal'
         case default
            write(*,'(a,1x,a)')   'Contour int. method ....... :', 'ERROR'
      end select
      if (psi_use_polyfit) then
          write(*,'(a,1x,a)')     'Phase space method ........ :', 'Polynomial-fit'
          write(*,'(a,1x,1i0)')   'Poly-fit order ............ :', psi_ncoeffs
          write(*,'(a,1x,1i0)')   'Number fit points ......... :', psi_npts
          write(*,'(a,1x,1es7.1)')'PSI tolerance ............. :', psi_si
      else
          write(*,'(a,1x,a)')        'Phase space method ........ :', 'Rational-fct-interp'
          write(*,'(a,1x,1i0,a,1i0)')'Rat-fct order (num,den) ... :', rat_order_n,', ', rat_order_d
          write(*,'(a,1x,1i0)')      'Number fit points ......... :', psi_npts
          write(*,'(a,1x,1i0)')      'PSI quadrature points ..... :', psi_gl_npts
      end if

      write(*,*)

      write(*,'(a,/,a)') 'CALCULATION PARAMETERS:', repeat('-', 30)
      write(*,'(a,1x,1f9.6)') 'Ground-state EQRPA ... (MeV):', gsenergy
      write(*,'(a,1x,1f9.6)') 'Q-value .............. (MeV):', qvalue
      write(*,'(a,1x,1f9.6)') 'Max EQRPA ............ (MeV):', eqrpamax

      write(*,*) 
      
      write(*,'(a,1x,1f4.2)')  'Vector coupling (gV) ...... :', gv
      write(*,'(a,1x,1f4.2)')  'Axial coupling (|gA|) ..... :', ga_input
      write(*,'(a,1x,1f4.2)')  'Lambda .................... :', lambda
      
      write(*,*) 
      
      write(*,'(a,1x,1f10.6)') 'R .................... (nat):', r
      write(*,'(a,1x,1f10.6)') 'Xi ................... (nat):', alpha*zfnum/2/r
      write(*,'(a,1x,1f10.6)') 'W0 max ............... (nat):', 1+qvalue/mec2
      write(*,'(a,1x,1f0.1)')  'Nucleon mass ......... (MeV):', mn
      
      write(*,*) 
      
      write(*,'(a,1f0.1," (", 1f8.6,")")') 'alpha*Z .............. (nat): 1/', 1/(alpha*zfnum), alpha*zfnum
      write(*,'(a,1f0.1," (", 1f8.6,")")') 'W0*R max ............. (nat): 1/', 1/((1+qvalue/mec2)*r), (1+qvalue/mec2)*r
      
      write(*,*) 
   end subroutine show_parameters
   
   !----------------------------------------------------------------------------
   ! Initialization header
   ! Author info, version, and run date.
   !
   ! Uses the globals:
   !  version
   !----------------------------------------------------------------------------
   subroutine show_header
      implicit none
      integer :: date_values(8)
      character(len=80) :: date
      
      write(*,'(a)') repeat('*', 80)
      write(*,'(a)') txtc(width=80,text='')
      write(*,'(a)') txtc(width=80,text='Allowed & First-Forbidden Beta-Decay Half-Lives Applying')
      write(*,'(a)') txtc(width=80,text='the Charge-Changing Finite Amplitude Method (pnFAM)')
      write(*,'(a)') txtc(width=80,text='')
      write(*,'(a)') txtc(width=80,text='M.T. Mustonen and T. Shafer')
      write(*,'(a)') txtc(width=80,text='The University of North Carolina at Chapel Hill')
      write(*,'(a)') txtc(width=80,text='')
      write(*,'(a)') txtc(width=80,text='Version: '//trim(version))
      write(*,'(a)') txtc(width=80,text='')
      
      call date_and_time(values=date_values)
      write(date,'(1i0,"/",1i0,"/",1i0,1x,1i0,":",1i2.2)') date_values(2), date_values(3), &
         date_values(1), date_values(5), date_values(6)
      
      write(*,'(a)') txtc(width=80,text='Run Date: '//trim(date))
      write(*,'(a)') txtc(width=80,text='')
      write(*,'(a)') repeat('*', 80)
   end subroutine show_header
   
   !----------------------------------------------------------------------------
   ! Center a string in a field of a fixed width
   !----------------------------------------------------------------------------
   function txtc(text, width)
      implicit none
      character(len=*), intent(in) :: text
      integer,          intent(in) :: width
      character(len=width) :: txtc
      integer :: lpad
      
      ! Extra padding on left for odd number
      lpad = (width-len_trim(text))/2 + mod(width-len_trim(text),2)
      txtc = repeat(' ', lpad)//adjustl(text)
   end function txtc
   
   !----------------------------------------------------------------------------
   ! Load a contour file
   !
   ! Uses globals:
   !  zinum
   !  anum
   !----------------------------------------------------------------------------
   subroutine load_contour_data(file, storage, label)
      implicit none
      
      character(len=*),   intent(in)  :: file
      character(len=*),   intent(in)  :: label
      type(contour_data), intent(out) :: storage
      
      character(len=20) :: ctmp, fam_op_name
      integer  :: ierr, fam_version, fam_k, fam_npoints, fam_nxterms, nstrengths, fam_z, fam_a
      real(dp) :: rtmp, dz_check
      logical  :: stop_import = .false.
      integer  :: gauss_avail, max_ind

      ! First, write the file under consideration to STDOUT
      ctmp = label//':'
      write(*,'(a20)', advance='no') ctmp
      
      if (len_trim(file) == 0) then
         write(*,'(a)') '--'
      else
         write(*,'(a)') trim(file)
      end if
      
      ! If no filename was passed, exit here
      if (len_trim(file) == 0) then
         return
      ! Otherwise, try and load the file
      else
         
         !----------------------------------------------------------------------
         ! Allocate the storage and read the data 
         !----------------------------------------------------------------------
         
         open(2, file=trim(file), status='old', form='unformatted', &
            access='sequential', iostat=ierr)
         
         ! Error reading => stop the calculation
         if (ierr /= 0) then
            write(*,'(3a)') 'ERROR in sub.load_contour_data: could not open the file "', &
               trim(file), '" for reading.'
            stop
         end if
         
         ! Format from pnfam_contour.f90:
         ! Version 1:
         !     integer             VERSION
         !     character(80)       OPERATOR
         !     integer             K, NR_POINTS, NXTERMS
         !     character(80) array XTERM_LABELS(:)        *** or EMPTY if no crossterms ***
         !     double        array TVALS(:)
         !     dcomplex      array DZ(:)
         !     dcomplex      array EQRPA(:)
         !     dcomplex      array SVALS(:,:)
         !
         ! Version 2:
         !     integer             VERSION
         !     integer             Z (parent), A
         !     character(80)       OPERATOR
         !     integer             K, NR_POINTS, NXTERMS
         !     character(80) array XTERM_LABELS(:)        *** or EMPTY if no crossterms ***
         !     double        array TVALS(:)
         !     dcomplex      array DZ(:)
         !     dcomplex      array EQRPA(:)
         !     dcomplex      array SVALS(:,:)         

         ! Version 3:
         !     version 2 + :
         !     integer             USE_GAULEG_CONTOUR (0 = no, 1= yes)
         !     double        array CTR_GAULEG_WEIGHTS(:)
         
         ! Version --- warn if the versions are different
         read(2) fam_version
         
         if (fam_version /= 1 .and. fam_version /= 2 .and. fam_version /= 3) then
            write(*,'(2(a,1i0),a)') 'WARNING in sub.load_contour_data: contour data format version&
               & is unexpected. Version = ', fam_version, ', expected 1, 2 or 3.'
         end if
         
         ! Version 2 and 3 only -- read Z and A and check them
         if (fam_version >= 2) then
            read(2) fam_z, fam_a
            
            ! zinum and anum are globals
            if (fam_z /= zinum) then
               stop_import = .true.
               write(*,'(/a)') 'ERROR in sub.load_contour_data: Z(parent) from&
                  & QRPA /= Z(parent) for this calculation.'
               write(*,'(2(a,1i0),a)') 'Z(QRPA) = ', fam_z, ', Z(beta) = ', zinum, '.'
            end if
            if (fam_a /= anum) then
               stop_import = .true.
               write(*,'(/a)') 'ERROR in sub.load_contour_data: A from&
                  & QRPA /= A for this calculation.'
               write(*,'(2(a,1i0),a)') 'A(QRPA) = ', fam_a, ', A(beta) = ', anum, '.'
               stop
            end if
            
            ! If there was an error, halt
            if (stop_import) stop
         end if
         
         ! Ignore the operator label for now
         read(2) fam_op_name
         
         ! Projection, number of points, and number of crossterms
         read(2) fam_k, fam_npoints, fam_nxterms
         
         ! Allocate permanent storage
         ! nstrengths = strength function istelf + crossterms
         nstrengths = fam_nxterms + 1
         
         allocate(storage%op(nstrengths), storage%t(fam_npoints), storage%dz(fam_npoints), &
            storage%eqrpa(fam_npoints), storage%str(fam_npoints, nstrengths), &
            storage%gl_wts(fam_npoints))
         
         ! BEGIN ASSIGNMENTS
         storage%label = trim(label)
         storage%k = fam_k
         
         storage%op(:) = ""
         storage%op(1) = trim(fam_op_name)
         
         storage%t(:) = 0;  storage%dz(:) = 0;  storage%eqrpa(:) = 0;  storage%str(:,:) = 0
         storage%gl_wts(:) = 0
         
         ! If there are crossterms, their labels live here. Otherwise, it's an empty line.
         if (fam_nxterms > 0) then
            read(2) (storage%op(ierr), ierr=2,nstrengths)
         else
            read(2)
         end if
         
         ! t-parameterization
         read(2) storage%t(:)
         
         ! dz/dt
         read(2) storage%dz(:)
         
         ! EQRPA values
         read(2) storage%eqrpa(:)
         
         ! Strengths and cross-terms
         read(2) storage%str(:,:)

         if (fam_version >= 3) then
            ! Gauss quadrature weights (if they're there...)
            read(2) gauss_avail

            if (quadrature == 'GAUSS') then
                if (gauss_avail == 0) then
                   write(*,'(a)') "ERROR: Cannot use GAUSS quadrature. Binary does not contain &
                       &points on a Gauss-Legendre contour."
                   stop
                else
                   read(2) storage%gl_wts(:)
                end if
            end if
         else
            if (quadrature == 'GAUSS') then
               write(*,'(a)') "ERROR: Cannot use GAUSS quadrature. Binary version incompatible."
               stop
            end if
         end if

         close(2)
         

         !----------------------------------------------------------------------
         ! Check that the data is valid
         !----------------------------------------------------------------------
         
         ! Check against our own values for EQRPA_max (at the 1 eV level)
         max_ind = maxloc(real(storage%eqrpa),1)

         if ( abs(aimag(storage%eqrpa(max_ind))) < epsilon(1.0) ) then
            rtmp = real(storage%eqrpa(max_ind),kind=dp) - eqrpamax
            if (abs(rtmp) > 1d-6) then
               write(*,'(a,1f0.6,a)') 'WARNING in sub.load_contour_data: EQRPA_max of the source&
                  & file does not match our EQRPA_max within 1 eV. Deviation = ', rtmp, ' MeV.'
            end if
         else
             write(*,'(a,1f0.6,a)') 'WARNING in sub.load_contour_data: max(EQRPA) does not lie on real axis.&
                 & Skipping check vs eqrpamax_hfb.'
         end if

         ! Check that dz/dt is sane
         dz_check = dz_from_data(storage)
         if (dz_check > 1d-7) then
            write(*,'(a,1es8.1,a)') 'WARNING in sub.load_contour_data: dz/dt did not&
               & pass consistency check, delta = ', dz_check, '.'
         end if
      end if
   end subroutine load_contour_data
   
   !----------------------------------------------------------------------------
   ! Compute the allowed decay rate.
   !----------------------------------------------------------------------------
   subroutine compute_allowed_rate
      use constants,  only : pi, ln2, ln2_kappa, mec2
      use phasespace, only : get_polyfit_coeffs, polyfit_eval, get_ratint_coeffs, ratint_eval_psi
      implicit none
      
      ! Pointers
      integer :: ip
      real(dp), pointer :: pt_rate
      type(contour_data), pointer :: pt_data
      
      ! Various constants
      integer     :: i, dim
      real(dp)    :: theta_k, const, psi_residue
      complex(dp) :: wmax, int
      complex(dp), dimension(:), allocatable :: psi

      ! Integration parameter
      real(dp), dimension(:), allocatable :: int_param
      
      
      ! Exit if no work to do
      if (.not.allocated(data_k0_f%t) .and. .not.allocated(data_k0_gt%t) .and. &
          .not.allocated(data_k1_gt%t)) then
         write(*,'(A)') 'No allowed-decay input data.'
         return
      end if
      
      
      ! Compute the coefficients for the allowed phase space approximation.
      ! Stored as the n=2 integral.
      ! Since all omega are > 0, the upper bound of electron energy should be
      ! m_e + EQRPA_max => 1 + EQRPAmax/m_e for the plynomial fit.
      write(*,'(a)') '* Computing phase space approximation with n = 2.'
      if (psi_use_polyfit) then
         coeffs(2,:) = get_polyfit_coeffs(zfnum, anum, psi_ncoeffs, psi_npts, 1+eqrpamax/mec2, 2, psi_si)
      else
         coeffs(2,:)  = get_ratint_coeffs(zfnum, anum, psi_npts, 1+eqrpamax/mec2, 2, psi_gl_npts)
      end if
      
      
      ! Header
      write(*,'(a)') repeat('-',80)
      write(*,'(1x,a,12x,a,4x,a,5x,a,12x,a,9x,a)') 'Contrib.', 'Theta_K', 'Dim', &
         'Int.Str.', 'Rate', 'HL (s)'
      write(*,'(a)') repeat('-',80)
      
      
      ! Set pointers to the data for each case and compute the rate
      ! ip=1: Fermi, K=0
      ! ip=2: GT, K=0
      ! ip=3: GT, K=1      
      nullify(pt_rate, pt_data)
      
      ! theta_k: 1 (K=0) or \sqrt{2} (K>0)
      ! const:   multiplying factor of the strength
      ! pt_data: pointer to this channel's data
      ! pt_rate: pointer to this channel's decay rate (polyfit + simps/trap)
      do ip=1, 3
         select case (ip)
            case (1)
               if (.not.allocated(data_k0_f%t)) cycle
               theta_k = 1
               const = theta_k**2
               pt_data => data_k0_f
               pt_rate => rate_k0_f
               call check_input(storage=data_k0_f, k=0, op='F', coord=pt_data)
            case (2)
               if (.not.allocated(data_k0_gt%t)) cycle
               theta_k = 1
               const = theta_k**2*lambda**2
               pt_data => data_k0_gt
               pt_rate => rate_k0_gt
               call check_input(storage=data_k0_gt, k=0, op='GT', coord=pt_data)
            case (3)
               if (.not.allocated(data_k1_gt%t)) cycle
               theta_k = sqrt(2.0_dp)
               const = theta_k**2*lambda**2
               pt_data => data_k1_gt
               pt_rate => rate_k1_gt
               call check_input(storage=data_k1_gt, k=1, op='GT', coord=pt_data)
            case default
               write(*,'(a)') 'ERROR in sub.compute_allowed_rate: invalid pointer number'
               stop
         end select
         
         ! Allocate space for the phase-space integrals
         dim = size(pt_data%t)
         allocate(psi(dim))
         psi(:) = 0

         ! Integration parameter is gl weights or theta depending on quadrature method
         allocate(int_param(dim))
         int_param(:) = 0
         if (quadrature == 'GAUSS') then
            int_param = pt_data%gl_wts(:)
         else
            int_param = pt_data%t(:)
         end if
         
         ! Compute the phase-space integrals
         ! Important: W0  =  1 + (Q-Ex)/m_e
         !                =  1 + (EQRPAmax - EQRPA)/m_e
         ! ---------------------------------------------------------------------
         do i=1, dim
            wmax = 1 + (eqrpamax-pt_data%eqrpa(i))/mec2
            if (psi_use_polyfit) then
               psi(i) = polyfit_eval(psi_ncoeffs, coeffs(2,:), wmax)
            else
               psi(i) = ratint_eval_psi(psi_npts, 1+eqrpamax/mec2, coeffs(2,:), wmax)
            end if
         end do  
         
         ! Integrate the contour
         ! Transition strength is const*str(:,1)*psi(:)
         int = cquad_wrapper(const*pt_data%str(:,1)*psi(:), pt_data%dz(:), int_param)
         pt_rate = -ln2_kappa*0.5d0*aimag(int)

         write(*,'(1x,a20,1f7.5,3x,1i4,1es13.2,4x,1es12.4,3x,1es12.5)') pt_data%label, theta_k, &
            dim, -0.5_dp*aimag(int), pt_rate, ln2/pt_rate
         
         ! Check the phase space integral
         psi_residue = -aimag(cquad_wrapper(psi(:), pt_data%dz(:), int_param))/(2*pi)

         if (abs(psi_residue) > 1d-8) then
            write(*,'(a,1es9.2,a)') 'WARNING in sub.compute_allowed_rate: phase space residue&
               & above 1e-8. Actual value: ', psi_residue, '.'
         end if
         
         deallocate(psi, int_param)
         nullify(pt_rate, pt_data)
      end do
   end subroutine compute_allowed_rate
   
   !----------------------------------------------------------------------------
   ! Compute the forbidden rates 
   !----------------------------------------------------------------------------
   subroutine compute_forbidden_rate
      use constants,  only : ln2, alpha, hbar_mec, mec2, ln2_kappa, mn
      use phasespace, only : get_polyfit_coeffs, polyfit_eval, get_ratint_coeffs, ratint_eval_psi
      use fermi,      only : gamma_1
      
      
      implicit none
      
      integer  :: k, dim, ii
      logical  :: work_to_do(0:2)
      real(dp) :: drate0, drate1, drate2
      complex(dp) :: wmax, int0, int1, int2, xp, xm
      
      ! Strength functions
      complex(dp) :: brs0, bps0, brs1, br, bp, brs2
      complex(dp) :: mrs0ps0, mrrs1, mprs1, mpr
      
      ! Fermi integrals
      complex(dp) :: fi1, fi2, fi3, fi4, fi5, fi6
      
      ! Muliplying matrix-element factors
      real(dp) :: theta_k
      real(dp) :: cps0, crs0, cp, cr, crs1, crs2
      
      ! Shape factors per-J (denoted CJ)
      complex(dp), dimension(:), allocatable :: C0, C1, C2
      
      ! Pointers
      real(dp) , pointer :: pt_rate
      type(contour_data), pointer :: pt_data

      ! Integration parameter
      real(dp), dimension(:), allocatable :: int_param
      
      
      ! Nb - these are .AND. so the calculation is skipped if any strength file
      ! is missing. This is to avoid weird stuff with missing operators, but it
      ! could be improved by spearating everything by J since we do not have
      ! cross-terms between operators with different J.
      work_to_do(0) = allocated(data_k0_rs0%t) .and. allocated(data_k0_ps0%t)  &
                .and. allocated(data_k0_rs1%t) .and. allocated(data_k0_r%t)    &
                .and. allocated(data_k0_p%t)   .and. allocated(data_k0_rs2%t)
      work_to_do(1) = allocated(data_k1_rs1%t) .and. allocated(data_k1_r%t)    &
                .and. allocated(data_k1_p%t)   .and. allocated(data_k1_rs2%t)
      work_to_do(2) = allocated(data_k2_rs2%t)
      
      ! Exit if there's nothing to do
      if (.not.any(work_to_do(:))) then
         write(*,'(A)') 'No forbidden-decay input data.'
         return
      end if
      
      
      ! Compute the coefficients for the forbidden phase space approximations,
      ! skipping n=2 if we already computed it (i.e. leading coefficient is not 0).
      ! Again, a safe domain for the approximation is W0 \in [1,1 + EQRPAmax/m_e].
      do ii=1, 6
         if (coeffs(ii,1) /= 0) cycle
         write(*,'(a,1i0,a)') '* Computing phase space approximation with n = ', ii, '.'
         if (psi_use_polyfit) then
            coeffs(ii,:) = get_polyfit_coeffs(zfnum, anum, psi_ncoeffs, psi_npts, 1+eqrpamax/mec2, ii, psi_si)
         else
            coeffs(ii,:)  = get_ratint_coeffs(zfnum, anum, psi_npts, 1+eqrpamax/mec2, ii, psi_gl_npts)
         end if
      end do
      
      
      ! Clear the pointers
      nullify(pt_rate, pt_data)
      
      
      ! Header
      write(*,'(A)') repeat('-',80)
      write(*,'(1x,a,3x,a,3x,a,5x,a,7x,a,7x,a,8x,a,16x,a)') 'K', 'Theta_K', 'Dim', 'J=0', &
         'J=1', 'J=2', 'Rate', 'HL (s)'
      write(*,'(A)') repeat('-',80)
      
      
      ! Main loop for computations, iterating over K projection
      do k=0, 2
         
         ! Skip if there isn't work for this K
         if (.not.work_to_do(k)) cycle
         
         ! Setup values based on K
         select case (k)
            case (0)
               pt_rate => rate_k0_ff
               theta_k = 1
               
               ! The pointer pt_data is used for computing dim and wmax,
               ! as well as for the t-parameterization of the contour integrals.
               if (allocated(data_k0_rs2%t)) then
                  pt_data => data_k0_rs2
               else if (allocated(data_k0_rs0%t)) then
                  pt_data => data_k0_rs0
               else if (allocated(data_k0_ps0%t)) then
                  pt_data => data_k0_ps0
               else if (allocated(data_k0_rs1%t)) then
                  pt_data => data_k0_rs1
               else if (allocated(data_k0_r%t)) then
                  pt_data => data_k0_r
               else if (allocated(data_k0_p%t)) then
                  pt_data => data_k0_p
               else
                  stop 'Inconsistency --- no allocated K=0 data.'
               end if
               
               ! Test to make sure we haven't linked the wrong strength
               call check_input(storage=data_k0_rs0, k=0, op='RS0', coord=pt_data)
               call check_input(storage=data_k0_ps0, k=0, op='PS0', coord=pt_data)
               call check_input(storage=data_k0_rs1, k=0, op='RS1', coord=pt_data)
               call check_input(storage=data_k0_r,   k=0, op='R',   coord=pt_data)
               call check_input(storage=data_k0_p,   k=0, op='P',   coord=pt_data)
               call check_input(storage=data_k0_rs2, k=0, op='RS2', coord=pt_data)
               
            case (1)
               pt_rate => rate_k1_ff
               theta_k = sqrt(2.0_dp)          
      
               ! The pointer pt_data is used for computing dim, dz, and wmax
               ! as well as for the t-parameterization of the contour integrals.
               if (allocated(data_k1_rs2%t)) then
                  pt_data => data_k1_rs2
               else if (allocated(data_k1_rs1%t)) then
                  pt_data => data_k1_rs1
               else if (allocated(data_k1_r%t)) then
                  pt_data => data_k1_r
               else if (allocated(data_k1_p%t)) then
                  pt_data => data_k1_p
               else
                  stop 'Inconsistency --- no allocated K=1 data.'
               end if
               
               ! Test to make sure we haven't linked the wrong strength
               call check_input(storage=data_k1_rs1, k=1, op='RS1', coord=pt_data)
               call check_input(storage=data_k1_r,   k=1, op='R',   coord=pt_data)
               call check_input(storage=data_k1_p,   k=1, op='P',   coord=pt_data)
               call check_input(storage=data_k1_rs2, k=1, op='RS2', coord=pt_data)
               
            case (2)
               pt_rate => rate_k2_ff
               theta_k = sqrt(2.0_dp)
               
               ! The pointer pt_data is used for computing dim, dz, and wmax
               ! as well as for the t-parameterization of the contour integrals.
               if (allocated(data_k2_rs2%t)) then
                  pt_data => data_k2_rs2
               else
                  stop 'Inconsistency --- no allocated K=2 data.'
               end if
               
               ! Test to make sure we haven't linked the wrong strength
               call check_input(storage=data_k2_rs2, k=2, op='RS2', coord=pt_data)
               
         end select
         
         
         dim = size(pt_data%t)
                  
         ! Multiplying factors
         ! See notes for where these come from, but in particular I've removed
         ! the 1/2 for the momentum operators since (pp + pn)/2M -> <p>/M.
         crs0 =  1.0_dp  / hbar_mec
         crs1 = -theta_k / hbar_mec
         crs2 =  theta_k / hbar_mec
         cr   =  theta_k * sqrt(3.0) / hbar_mec
         cps0 = -1.0_dp  * hbar_mec  * mec2 / Mn
         cp   = -theta_k * hbar_mec  * mec2 / Mn
         
         
         ! Construct the contour integrand, splitting off J=0, J=1, and J=2 parts
         allocate(C0(dim), C1(dim), C2(dim))
         C0(:) = 0;  C1(:) = 0;  C2(:) = 0

         ! Integration parameter is gl weights or t depending on quadrature method
         allocate(int_param(dim))
         int_param(:) = 0
         if (quadrature == 'GAUSS') then
            int_param = pt_data%gl_wts(:)
         else
            int_param = pt_data%t(:)
         end if

         do ii=1, dim
            ! Maximum decay energy and related quantities
            ! Important: W0  =  1 + (Q-Ex)/m_e
            !                =  1 + (EQRPAmax - EQRPA)/m_e
            ! ---------------------------------------------------------------------
            wmax = 1 + (eqrpamax-pt_data%eqrpa(ii))/mec2
            
            xp = wmax/3.0_dp + 0.5_dp*alpha*zfnum/r
            xm = wmax/3.0_dp - 0.5_dp*alpha*zfnum/r
            
            brs0 = 0;  bps0 = 0;  brs1 = 0;  br = 0;  bp = 0;  brs2 = 0
            mrs0ps0 = 0;  mrrs1 = 0;  mprs1 = 0;  mpr = 0
            
            ! Strength functions and cross-terms
            ! Note on dimensions ---
            ! 1) Matrix elements of r should have units of [fm] so they have to
            !    be divided by the reduced electron Compton wavelength (hbar c)/(me c^2)
            ! 2) Matrix elements of p should have units of [fm^-1] so they
            !    should be divided by the nucleon mass [MeV/c^2] and multiplied
            !    by (hbar c) [MeV fm].
            select case (k)
               case (0)
                  ! Strength functions
                  if (allocated(data_k0_rs0%t)) brs0 = crs0**2 * data_k0_rs0%str(ii,1)
                  if (allocated(data_k0_ps0%t)) bps0 = cps0**2 * data_k0_ps0%str(ii,1)
                  if (allocated(data_k0_rs1%t)) brs1 = crs1**2 * data_k0_rs1%str(ii,1)
                  if (allocated(data_k0_r%t))   br   = cr**2   * data_k0_r%str(ii,1)
                  if (allocated(data_k0_p%t))   bp   = cp**2   * data_k0_p%str(ii,1)
                  if (allocated(data_k0_rs2%t)) brs2 = crs2**2 * data_k0_rs2%str(ii,1)
                  
                  ! Cross-terms
                  call assert_equal_strength(data_k0_rs0%str(ii,3), data_k0_ps0%str(ii,2), data_k0_rs0%op(3))
                  call assert_equal_strength(data_k0_r%str(ii,3),   data_k0_rs1%str(ii,2), data_k0_r%op(3))
                  call assert_equal_strength(data_k0_p%str(ii,3),   data_k0_rs1%str(ii,4), data_k0_p%op(3))
                  call assert_equal_strength(data_k0_r%str(ii,4),   data_k0_p%str(ii,2),   data_k0_r%op(4))

                  ! [RS]0 x [PS]0
                  if ((trim(data_k0_rs0%op(3)) /= 'RS0xPS0')  .or. &
                      (trim(data_k0_ps0%op(2)) /= 'PS0xRS0')) then
                     stop 'Forbidden cross-term label is inconsistent.'
                  end if
                  
                  if (allocated(data_k0_rs0%t)) then
                     mrs0ps0 = crs0*cps0 * data_k0_rs0%str(ii,3)
                  else if (allocated(data_k0_ps0%t)) then
                     mrs0ps0 = crs0*cps0 * data_k0_ps0%str(ii,2)
                  end if
                  
                  ! R x [RS]1
                  if ((trim(data_k0_r%op(3))   /= 'RxRS1') .or. &
                      (trim(data_k0_rs1%op(2)) /= 'RS1xR')) then
                     stop 'Forbidden cross-term label is inconsistent.'
                  end if
                  
                  if (allocated(data_k0_r%t)) then
                     mrrs1 = cr*crs1 * data_k0_r%str(ii,3)
                  else if (allocated(data_k0_rs1%t)) then
                     mrrs1 = cr*crs1 * data_k0_rs1%str(ii,2)
                  end if

                  ! P x [RS]1
                  if ((trim(data_k0_p%op(3))   /= 'PxRS1') .or. &
                      (trim(data_k0_rs1%op(4)) /= 'RS1xP')) then
                     stop 'Forbidden cross-term label is inconsistent.'
                  end if
                  
                  if (allocated(data_k0_p%t)) then
                     mprs1 = cp*crs1 * data_k0_p%str(ii,3)
                  else if (allocated(data_k0_rs1%t)) then
                     mprs1 = cp*crs1 * data_k0_rs1%str(ii,4)
                  end if
                  
                  ! P x R
                  if ((trim(data_k0_r%op(4)) /= 'RxP') .or. &
                      (trim(data_k0_p%op(2)) /= 'PxR')) then
                     stop 'Forbidden cross-term label is inconsistent.'
                  end if
                  
                  if (allocated(data_k0_r%t)) then
                     mpr = cp*cr * data_k0_r%str(ii,4)
                  else if (allocated(data_k0_p%t)) then
                     mpr = cp*cr * data_k0_p%str(ii,2)
                  end if
                  
               case (1)
                  ! Strength functions
                  if (allocated(data_k1_rs1%t)) brs1 = crs1**2 * data_k1_rs1%str(ii,1)
                  if (allocated(data_k1_r%t))   br = cr**2 * data_k1_r%str(ii,1)
                  if (allocated(data_k1_p%t))   bp = cp**2 * data_k1_p%str(ii,1)
                  if (allocated(data_k1_rs2%t)) brs2 = crs2**2 * data_k1_rs2%str(ii,1)
                  
                  ! Cross-terms
                  call assert_equal_strength(data_k1_r%str(ii,3), data_k1_rs1%str(ii,2), data_k1_r%op(3))
                  call assert_equal_strength(data_k1_p%str(ii,3), data_k1_rs1%str(ii,4), data_k1_p%op(3))
                  call assert_equal_strength(data_k1_r%str(ii,4), data_k1_p%str(ii,2),   data_k1_r%op(4))

                  ! R x [RS]1
                  if ((trim(data_k1_r%op(3))   /= 'RxRS1') .or. &
                      (trim(data_k1_rs1%op(2)) /= 'RS1xR')) then
                     stop 'Forbidden cross-term label is inconsistent.'
                  end if
                  
                  if (allocated(data_k1_r%t)) then
                     mrrs1 = cr*crs1 * data_k1_r%str(ii,3)
                  else if (allocated(data_k1_rs1%t)) then
                     mrrs1 = cr*crs1 * data_k1_rs1%str(ii,2)
                  end if
                  
                  ! P x [RS]1
                  if ((trim(data_k1_p%op(3))   /= 'PxRS1') .or. &
                      (trim(data_k1_rs1%op(4)) /= 'RS1xP')) then
                     stop 'Forbidden cross-term label is inconsistent.'
                  end if
                  
                  if (allocated(data_k1_p%t)) then
                     mprs1 = cp*crs1 * data_k1_p%str(ii,3)
                  else if (allocated(data_k1_rs1%t)) then
                     mprs1 = cp*crs1 * data_k1_rs1%str(ii,4)
                  end if
                  
                  ! P x R
                  if ((trim(data_k1_r%op(4)) /= 'RxP') .or. &
                      (trim(data_k1_p%op(2)) /= 'PxR')) then
                     stop 'Forbidden cross-term label is inconsistent.'
                  end if
                  
                  if (allocated(data_k1_r%t)) then
                     mpr = cp*cr * data_k1_r%str(ii,4)
                  else if (allocated(data_k1_p%t)) then
                     mpr = cp*cr * data_k1_p%str(ii,2)
                  end if
               case (2)  
                  if (allocated(data_k2_rs2%t)) brs2 = crs2**2 * data_k2_rs2%str(ii,1)
            end select
            
            ! Fermi integrals
            if (psi_use_polyfit) then
               fi1 = polyfit_eval(psi_ncoeffs, coeffs(1,:), wmax)
               fi2 = polyfit_eval(psi_ncoeffs, coeffs(2,:), wmax)
               fi3 = polyfit_eval(psi_ncoeffs, coeffs(3,:), wmax)
               fi4 = polyfit_eval(psi_ncoeffs, coeffs(4,:), wmax)
               fi5 = polyfit_eval(psi_ncoeffs, coeffs(5,:), wmax)
               fi6 = polyfit_eval(psi_ncoeffs, coeffs(6,:), wmax)
            else
               fi1 = ratint_eval_psi(psi_npts, 1+eqrpamax/mec2, coeffs(1,:), wmax)
               fi2 = ratint_eval_psi(psi_npts, 1+eqrpamax/mec2, coeffs(2,:), wmax)
               fi3 = ratint_eval_psi(psi_npts, 1+eqrpamax/mec2, coeffs(3,:), wmax)
               fi4 = ratint_eval_psi(psi_npts, 1+eqrpamax/mec2, coeffs(4,:), wmax)
               fi5 = ratint_eval_psi(psi_npts, 1+eqrpamax/mec2, coeffs(5,:), wmax)
               fi6 = ratint_eval_psi(psi_npts, 1+eqrpamax/mec2, coeffs(6,:), wmax)
            end if

            ! Complete shape factor
            C0(ii) = 0;  C1(ii) = 0;  C2(ii) = 0
            ! J = 0
            C0(ii) = C0(ii) + lambda**2*(bps0 + brs0*(1/9.0_dp + xp**2) + mrs0ps0*2*xp)*fi2
            C0(ii) = C0(ii) - 2/3.0_dp*lambda**2*(brs0*xp + mrs0ps0)*fi1
            ! J = 1
            C1(ii) = C1(ii) - 2/9.0_dp*(xp*br - 2*lambda**2*xm*brs1    &
                            + lambda*sqrt(2.0_dp)*(xp-xm)*mrrs1        &
                            - sqrt(3.0_dp)*mpr - lambda*sqrt(6.0_dp)*mprs1)*fi1
            C1(ii) = C1(ii) + (bp + 1/3.0_dp*xp**2*br + 2/3.0_dp*xm**2*lambda**2*brs1)*fi2
            C1(ii) = C1(ii) + sqrt(2.0_dp/3.0_dp)*(2*lambda*xm*mprs1 - sqrt(2.0_dp)*xp*mpr &
                            - 2/sqrt(3.0_dp)*lambda*xm*xp*mrrs1)*fi2
            C1(ii) = C1(ii) - 8/27.0_dp*(lambda**2*brs1 + lambda/sqrt(2.0_dp)*mrrs1)*gamma_1(zfnum)*fi2
            C1(ii) = C1(ii) + 1/27.0_dp*(br + 2*lambda**2*brs1 + 2*sqrt(2.0_dp)*lambda*mrrs1)*fi2
            C1(ii) = C1(ii) + (4*sqrt(2.0_dp)/9.0_dp*lambda*xp*mrrs1 - 8/9.0_dp*lambda**2*xm*brs1 &
                            -  4*sqrt(2.0_dp)/3.0_dp/sqrt(3.0_dp)*lambda*mprs1)*fi3
            C1(ii) = C1(ii) + 8/27.0_dp*lambda**2*brs1*fi4
            C1(ii) = C1(ii) + 1/27.0_dp*(2*br + brs1*lambda**2 + 2*sqrt(2.0_dp)*mrrs1*lambda)*fi5
            C1(ii) = C1(ii) + 1/27.0_dp*(2*br + brs1*lambda**2 - 2*sqrt(2.0_dp)*mrrs1*lambda)*fi6
            ! J = 2
            C2(ii) = C2(ii) + 1/9.0_dp*lambda**2*brs2*(fi5+fi6)
         end do
         
         ! Integrate
         int0 = cquad_wrapper(C0, pt_data%dz, int_param)
         int1 = cquad_wrapper(C1, pt_data%dz, int_param)
         int2 = cquad_wrapper(C2, pt_data%dz, int_param)
         
         ! These hold the J-specific decay rates. Recall that the factor of
         ! 1/pi is already removed in the pnFAM solver
         drate0  = -ln2_kappa*0.5_dp*aimag(int0)
         drate1  = -ln2_kappa*0.5_dp*aimag(int1)
         drate2  = -ln2_kappa*0.5_dp*aimag(int2)
         
         ! The total rate for this value of K is the sum over J
         pt_rate = drate0+drate1+drate2
         
         ! J-specific total rates
         rate_j0_ff = rate_j0_ff + drate0
         rate_j1_ff = rate_j1_ff + drate1
         rate_j2_ff = rate_j2_ff + drate2
         
         ! J- and K-specific total rates
         select case (k)
            case (0)
               rate_j0_k0_ff = drate0
               rate_j1_k0_ff = drate1
               rate_j2_k0_ff = drate2
            case (1)
               rate_j1_k1_ff = drate1
               rate_j2_k1_ff = drate2
            case (2)
               rate_j2_k2_ff = drate2
         end select
         
         ! Status update
         write(*,'(1x,1i1,3x,1f7.5,2x,1i4,3f10.5,4x,1es10.4,3x,1f14.5)') k, theta_k, dim, &
            drate0, drate1, drate2, pt_rate, ln2/pt_rate
         
         ! Check for the correct J-contributions and positive rates
         if (drate0 < 0 .or. drate1 < 0 .or. drate2 < 0) then
            write(*,'(5x,a)') '>>> ERROR: a computed rate is less than zero!'
         end if
         if (k > 0 .and. drate0 > 0) then
            write(*,'(5x,a)') '>>> WARNING: incorrect J=0 contribution!'
         end if
         if (k > 1 .and. drate1 > 0) then
            write(*,'(5x,a)') '>>> WARNING: incorrect J=1 contribution!'
         end if
         
         deallocate(C0, C1, C2, int_param)
         nullify(pt_rate, pt_data)
      end do
      
   end subroutine compute_forbidden_rate
   
   !----------------------------------------------------------------------------
   ! Display the computed quantities 
   !----------------------------------------------------------------------------
   subroutine show_results(store)
      use constants, only : ln2, alpha, mec2, mn
      implicit none
      
      logical, intent(in) :: store
      
      ! Quantities
      real(dp) :: rate0p, rate1p, rate0m, rate1m, rate2m, hl0p, hl1p, hl0m, hl1m, hl2m
      real(dp) :: rate_allowed, rate_forbidden, rate_total, per0p, per1p, per0m, per1m, per2m
      real(dp) :: hl_allowed, hl_forbidden, hl_total, per_allowed, per_forbidden
      real(dp) :: ratej0p, ratej1p, ratej0m, ratej1m, ratej2m
      
      ! Formats
      character(len=*), parameter :: &
         fmt1 = '("Rate ",a," = ",1f9.6," s^(-1)  [ ",1f6.2,"% ]    t1/2 ",a," = ",1f0.6," s")'
      character(len=*), parameter :: &
         fmt2 = '(a," = ",1f0.6," s  [ ",1f0.2,"% ]")'
      character(len=*), parameter :: &
         fmt3 = '(a," = ",1f0.6," s")'
      character(len=*), parameter :: &
         fmt4 = '(a," = ",1f0.6," s  [ ",sp,1f0.2,"% ]")'
      character(len=*), parameter :: &
         fmt5 = '("Rate ",a," = ",1f9.6," s^(-1)  [ ",1f6.2,"% ]    t1/2 ",a," = ",1f0.6," s")'      
      
      ! For storing files
      integer :: ierr, date_values(8)
      
      
      ! Broken down by K
      ! ------------------------------------------------------------------------
      rate0p = rate_k0_f + rate_k0_gt
      rate1p = rate_k1_gt
      rate0m = rate_k0_ff
      rate1m = rate_k1_ff
      rate2m = rate_k2_ff

      !! For debugging only --- avoid dividing by zero
      !if (rate0p == 0) rate0p = epsilon(1.0d0)
      !if (rate1p == 0) rate1p = epsilon(1.0d0)
      !if (rate0m == 0) rate0m = epsilon(1.0d0)
      !if (rate1m == 0) rate1m = epsilon(1.0d0)
      !if (rate2m == 0) rate2m = epsilon(1.0d0)
      
      hl0p = ln2/rate0p
      hl1p = ln2/rate1p
      hl0m = ln2/rate0m
      hl1m = ln2/rate1m
      hl2m = ln2/rate2m
      
      rate_allowed   = rate0p + rate1p
      rate_forbidden = rate0m + rate1m + rate2m
      rate_total = rate_allowed + rate_forbidden
      
      hl_allowed   = ln2/rate_allowed
      hl_forbidden = ln2/rate_forbidden
      hl_total = ln2/rate_total
      
      per0p = 100*rate0p/rate_total
      per1p = 100*rate1p/rate_total
      per0m = 100*rate0m/rate_total
      per1m = 100*rate1m/rate_total
      per2m = 100*rate2m/rate_total
      
      per_allowed   = 100*rate_allowed/rate_total
      per_forbidden = 100*rate_forbidden/rate_total
      ! ------------------------------------------------------------------------
      
      ! Broken down by J
      ! ------------------------------------------------------------------------
      ratej0p = rate_k0_f
      ratej1p = rate_k0_gt + rate_k1_gt
      ratej0m = rate_j0_ff
      ratej1m = rate_j1_ff
      ratej2m = rate_j2_ff
      
      !! For debugging only --- avoid dividing by zero
      !if (ratej0p == 0) ratej0p = epsilon(1.0d0)
      !if (ratej1p == 0) ratej1p = epsilon(1.0d0)
      !if (ratej0m == 0) ratej0m = epsilon(1.0d0)
      !if (ratej1m == 0) ratej1m = epsilon(1.0d0)
      !if (ratej2m == 0) ratej2m = epsilon(1.0d0)
      ! ------------------------------------------------------------------------
      
   
      write(*,'(/)')
      write(*,'(a,/,a)') 'SUMMARY OF RESULTS:', repeat('-',30)
      write(*,'(a)') 'ALLOWED by K'
      write(*,fmt1) '0+', rate0p, per0p, '0+', hl0p
      write(*,fmt1) '1+', rate1p, per1p, '1+', hl1p
      
      write(*,'(a)') ''
      write(*,'(a)') 'FORBIDDEN by K'
      write(*,fmt1) '0-', rate0m, per0m, '0-', hl0m
      write(*,fmt1) '1-', rate1m, per1m, '1-', hl1m
      write(*,fmt1) '2-', rate2m, per2m, '2-', hl2m
      
      write(*,'(a)') ''
      write(*,'(a)') 'CONTRIBUTIONS by J'
      write(*,fmt5) '0+', ratej0p, 100*ratej0p/rate_total, '0+', ln2/ratej0p
      write(*,fmt5) '1+', ratej1p, 100*ratej1p/rate_total, '1+', ln2/ratej1p
      write(*,fmt5) '0-', ratej0m, 100*ratej0m/rate_total, '0-', ln2/ratej0m
      write(*,fmt5) '1-', ratej1m, 100*ratej1m/rate_total, '1-', ln2/ratej1m
      write(*,fmt5) '2-', ratej2m, 100*ratej2m/rate_total, '2-', ln2/ratej2m
      
      write(*,'(/)')
      write(*,'(a,/,a)') 'PARTIAL HALF-LIVES:', repeat('-',30)
      write(*,fmt2) 'Allowed  ', hl_allowed, per_allowed
      write(*,fmt2) 'Forbidden', hl_forbidden, per_forbidden
      
      write(*,'(/)')      
      write(*,'(a,/,a)') 'TOTAL HALF-LIFE:', repeat('-',30)
      write(*,'("Allowed + FF = ", 1f0.8," s")') hl_total
      write(*,*)
      
      
      ! If we are to store the data, do so
      if (store .eqv. .true.) then
         open(10, file=trim(outfile), status='unknown', iostat=ierr)
      
         ! Header information
         write(10,'(a)')  '# Nuclear Beta Decay Rates and Half-Lives'
         write(10,'(2a)') '# pnFAM code version: ', version
      
         call date_and_time(values=date_values)
         write(10,'(a,1i0,"/",1i0,"/",1i0,1x,1i0,":",1i2.2)') '# Run date: ', date_values(2), &
            date_values(3), date_values(1), date_values(5), date_values(6)
      
         write(10,'(a)') '# '
      
         select case (quadrature)
            case("GAUSS")
               write(10,'(a,1x,a)') '# Contour integration method ..... :', 'Gauss-Legendre'
            case("SIMPSON")
               write(10,'(a,1x,a)') '# Contour integration method ..... :', 'Simpsons 3/8'
            case("TRAP")
               write(10,'(a,1x,a)') '# Contour integration method ..... :', 'Trapezoidal'
            case default
               write(10,'(a,1x,a)') '# Contour integration method ..... :', 'ERROR'
         end select
         if (psi_use_polyfit) then
             write(10,'(a,1x,a)')   '# Phase space method ............. :', 'Polynomial-fit'
             write(10,'(a,1i10)',  advance='no')'#  poly_order  =', psi_ncoeffs
             write(10,'(a,1i10)',  advance='no')',     fit_pts =', psi_npts
             write(10,'(a,1es10.1)',advance='no')',    psi_tol =', psi_si
         else
             write(10,'(a,1x,a)')  '# Phase space method ............. :', 'Rational-fct-interp'
             write(10,'(a,1i10)',advance='no')  '#  rfi_order_n =',       rat_order_n
             write(10,'(a,1i10)',advance='no')  ', rfi_order_d =', rat_order_d
             write(10,'(a,1i10)',advance='no')  ', interp_pts =',  psi_npts
             write(10,'(a,1i10)',advance='no')  ', psi_quad_pts=',  psi_gl_npts
         end if

         write(10,'(a)') 
         write(10,'(a,  i10)',advance='no')     '#  Z_(ini)     =', zinum
         write(10,'(a,  i10)',advance='no')     ',           A =', anum
         write(10,'(a,  i10)',advance='no')     ',    Z_(fin) =', zfnum
         write(10,'(a)') 
         write(10,'(a,f10.5)',advance='no')     '#  GS_energy   =', gsenergy
         write(10,'(a,f10.5)',advance='no')     ',     Q-value =', qvalue
         write(10,'(a,f10.5)',advance='no')     ',   EQRPAmax =', eqrpamax
         write(10,'(a,f10.5)',advance='no')     ',     lambda  =', lambda
         write(10,'(a)') 
         write(10,'(a,f10.5)',advance='no')     '#  g_A         =', ga_input
         write(10,'(a,f10.5)',advance='no')     ',         g_V =', gv
         write(10,'(a,f10.5)',advance='no')     ',    Nucl. M =', mn
         write(10,'(a,f10.5)',advance='no')     ',     alpha*Z =', alpha*zfnum
         write(10,'(a)')
         write(10,'(a,f10.5)',advance='no')     '#  R           =', r
         write(10,'(a,f10.5)',advance='no')     ',          Xi =', alpha*zfnum/2/r
         write(10,'(a,f10.5)',advance='no')     ',     W0 max =', 1+qvalue/mec2
         write(10,'(a,f10.5)',advance='no')     ',        W0*R =', (1+qvalue/mec2)*r
         write(10,'(a)') 
         write(10,'(a)') '# '
         write(10,'(a)') '#            [RATE] (s^-1)           [HALF-LIFE] (s)'
         
         !! For debugging only --- avoid dividing by zero
         !if (rate_total     == 0) rate_total     = epsilon(1.0d0)
         !if (rate_allowed   == 0) rate_allowed   = epsilon(1.0d0)
         !if (rate_forbidden == 0) rate_forbidden = epsilon(1.0d0)
         !if (rate_k0_f      == 0) rate_k0_f      = epsilon(1.0d0)
         !if (rate_k0_gt     == 0) rate_k0_gt     = epsilon(1.0d0)
         !if (rate_k1_gt     == 0) rate_k1_gt     = epsilon(1.0d0)
         !if (rate_k0_ff     == 0) rate_k0_ff     = epsilon(1.0d0)
         !if (rate_k1_ff     == 0) rate_k1_ff     = epsilon(1.0d0)
         !if (rate_k2_ff     == 0) rate_k2_ff     = epsilon(1.0d0)
         !if (rate_j0_ff     == 0) rate_j0_ff     = epsilon(1.0d0)
         !if (rate_j1_ff     == 0) rate_j1_ff     = epsilon(1.0d0)
         !if (rate_j2_ff     == 0) rate_j2_ff     = epsilon(1.0d0)         
         
         ! Rates
         write(10,'(2es26.16,4x,1a)') rate_total,     ln2/rate_total,     '# Total'
         write(10,'(2es26.16,4x,1a)') rate_allowed,   ln2/rate_allowed,   '# Allowed'
         write(10,'(2es26.16,4x,1a)') rate_forbidden, ln2/rate_forbidden, '# Forbidden'
         write(10,'(2es26.16,4x,1a)') rate_k0_f,      ln2/rate_k0_f,      '# Fermi'
         write(10,'(2es26.16,4x,1a)') rate_k0_gt,     ln2/rate_k0_gt,     '# Gamow-Teller, K=0'
         write(10,'(2es26.16,4x,1a)') rate_k1_gt,     ln2/rate_k1_gt,     '# Gamow-Teller, K=1'
         write(10,'(2es26.16,4x,1a)') rate_k0_ff,     ln2/rate_k0_ff,     '# Forbidden, K=0'
         write(10,'(2es26.16,4x,1a)') rate_k1_ff,     ln2/rate_k1_ff,     '# Forbidden, K=1'
         write(10,'(2es26.16,4x,1a)') rate_k2_ff,     ln2/rate_k2_ff,     '# Forbidden, K=2'
         write(10,'(2es26.16,4x,1a)') rate_j0_ff,     ln2/rate_j0_ff,     '# Forbidden, J=0'         
         write(10,'(2es26.16,4x,1a)') rate_j1_ff,     ln2/rate_j1_ff,     '# Forbidden, J=1'         
         write(10,'(2es26.16,4x,1a)') rate_j2_ff,     ln2/rate_j2_ff,     '# Forbidden, J=2'         
         write(10,'(2es26.16,4x,1a)') rate_j0_k0_ff,  ln2/rate_j0_k0_ff,  '# Forbidden, (J,K)=(0,0)'
         write(10,'(2es26.16,4x,1a)') rate_j1_k0_ff,  ln2/rate_j1_k0_ff,  '# Forbidden, (J,K)=(1,0)'
         write(10,'(2es26.16,4x,1a)') rate_j1_k1_ff,  ln2/rate_j1_k1_ff,  '# Forbidden, (J,K)=(1,1)'
         write(10,'(2es26.16,4x,1a)') rate_j2_k0_ff,  ln2/rate_j2_k0_ff,  '# Forbidden, (J,K)=(2,0)'
         write(10,'(2es26.16,4x,1a)') rate_j2_k1_ff,  ln2/rate_j2_k1_ff,  '# Forbidden, (J,K)=(2,1)'
         write(10,'(2es26.16,4x,1a)') rate_j2_k2_ff,  ln2/rate_j2_k2_ff,  '# Forbidden, (J,K)=(2,2)'
         
         close(10)
      end if
   end subroutine show_results
   
   !----------------------------------------------------------------------------
   ! Derive dz/dt from a strength array as a test of consistency
   ! NB: only for circular contours with z = z0 + r exp(it)
   !----------------------------------------------------------------------------
   function dz_from_data(s) result(delta)
      use constants, only : iu
      implicit none
      
      type(contour_data), intent(in) :: s
      real(dp) :: delta
      
      real(dp)    :: rr
      complex(dp) :: z0
      
      delta = 0
      
      ! Compute r and z0 by sampling two points using the ansatz z = z0 + r exp(it)
      !         z1 - r*exp(i*t1) = z2 - r*exp(i*t2) --> solve for r

      rr = real((s%eqrpa(1)-s%eqrpa(2))/(exp(iu*s%t(1))-exp(iu*s%t(2))), kind=dp)
      z0 = s%eqrpa(1)-rr*exp(iu*s%t(1))
      
      ! Max deviation between expected dz/dt and stored dz/dt
      delta = maxval(abs(s%dz(:)-iu*rr*exp(iu*s%t(:))))

   end function dz_from_data
   
   !----------------------------------------------------------------------------
   ! Check a data storage type against our expectations.
   ! This should catch most of the errors we could make in linking the wrong
   ! file in the namelist.
   !----------------------------------------------------------------------------
   subroutine check_input(storage, k, op, coord)
      implicit none
      
      integer,            intent(in) :: k
      character(len=*),   intent(in) :: op
      type(contour_data), intent(in) :: storage, coord
         
      integer :: ierr = 0
      
      ! Only perform the tests if the storage is being used
      if (.not.allocated(storage%t)) then
         return
      else
         ! Check K
         if (storage%k /= k) then
            ierr = ierr + 1
            write(*,'(2(a,1i0),a)') 'ERROR in sub.check_input: Ks unequal. K_file = ', &
               storage%k, ', K_expect = ', k, '.'
            write(*,'(3a)') 'Operator = "', trim(storage%op(1)), '".'
         end if
         
         ! Check the operator code
         if (trim(storage%op(1)) /= trim(op)) then
            ierr = ierr + 1
            write(*,'(5a)') 'ERROR in sub.check_input: OPs unequal. OP_file = "', &
               trim(storage%op(1)), '", OP_expect = "', trim(op), '".'
         end if
         
         ! Check the matrix dimension
         if (size(storage%t) /= size(coord%t)) then
            ierr = ierr + 1
            write(*,'(2(a,1i0),a)') 'ERROR in sub.check_input: storage dimensions unequal.&
               & DIM_file = ', size(storage%t), ', DIM_expect = ', size(coord%t), '.'            
         end if
         
         ! Check the t-parameterization
         if (any(abs((storage%t(:)-coord%t(:))/coord%t(:)) > 1d-10) .eqv. .true.) then
            ierr = ierr + 1
            write(*,'(2(a,1i0),a)') 'ERROR in sub.check_input: t-parameterization unequal.&
               & Relative difference greater than 1e-10.'            
         end if
         
         ! Check the dz/dt arrays
         if (any(abs((storage%dz(:)-coord%dz(:))/coord%dz(:)) > 1d-10) .eqv. .true.) then
            ierr = ierr + 1
            write(*,'(2(a,1i0),a)') 'ERROR in sub.check_input: dz/dt unequal.&
               & Relative difference greater than 1e-10.'            
         end if
         
         ! Check the eqrpa arrays
         if (any(abs((storage%eqrpa(:)-coord%eqrpa(:))/coord%eqrpa(:)) > 1d-10) .eqv. .true.) then
            ierr = ierr + 1
            write(*,'(2(a,1i0),a)') 'ERROR in sub.check_input: eqrpa vectors unequal.&
               & Relative difference greater than 1e-10.'            
         end if
         
         ! If there were any errors, stop execution
         if (ierr > 0) stop
      end if
   end subroutine check_input
   
   !----------------------------------------------------------------------------
   ! A wrapper to use the correct integrator depending on the global variable
   ! 'quadrature'.
   !----------------------------------------------------------------------------
   function cquad_wrapper(f, dz, int_param) result(value)
      use complex_quadrature, only : cquad_simpson, cquad_trapezoidal, cquad_gauss
      implicit none
   
      real(dp),    dimension(:), intent(in) :: int_param
      complex(dp), dimension(:), intent(in) :: f, dz
      complex(dp) :: value
      value = 0

      select case (quadrature)
         case("GAUSS")
            value = cquad_gauss(f=f, dz=dz, wts=int_param)
         case("SIMPSON")
            value = cquad_simpson(f=f, dz=dz, t=int_param)
         case("TRAP")
            value = cquad_trapezoidal(f=f, dz=dz, t=int_param)
         case default
            write(error_unit,'(3a)') "ERROR: Unimplemented quadrature type: ", &
                trim(quadrature), "!"
            stop
      end select
      
   end function cquad_wrapper

   !----------------------------------------------------------------------------
   ! Check the real and imaginary parts of the cross-terms separately with the
   ! goal of verifying that e.g. [RS]0 x [PS]0 == [PS]0 x [RS]0. I'm trying to
   ! check that there are ~4 digits of agreement, i.e.
   ! (magnitude) x (rel. difference) < 0.1%.
   !----------------------------------------------------------------------------
   subroutine assert_equal_strength(str_a, str_b, label)
      implicit none
      
      real(dp), parameter :: tolerance_fail = 1d-3
      real(dp), parameter :: tolerance_warn = 5d-4
      
      character(len=*), parameter :: fmt = &
         '(A,1X,A,". |REL DIFF| =",1ES9.2,", |ABS DIFF| =",1ES9.2,", |DELTA| =",1ES9.2,".")'
      
      complex(dp),      intent(in) :: str_a, str_b
      character(len=*), intent(in) :: label
      
      logical  :: failed = .false.
      real(dp) :: rel_diff, abs_diff, delta, mean
      
      ! Real
      mean     = 0.5_dp*abs(real(str_a,kind=dp)+real(str_b,kind=dp))
      rel_diff = abs((real(str_b,kind=dp)-real(str_a,kind=dp))/mean)
      abs_diff = abs(real(str_b,kind=dp)-real(str_a,kind=dp))
      delta    = rel_diff*mean
      
      if (delta > tolerance_warn) then
         if (delta > tolerance_fail) then
            failed = .true.
            write(*,fmt) 'ERROR: Re[X] are inconsistent for', trim(label), rel_diff, abs_diff, delta
         else
            write(*,fmt) 'WARNING: Re[X] may be inconsistent for', trim(label), rel_diff, abs_diff, delta
         end if
      end if
      
      ! Imaginary
      mean     = 0.5_dp*abs(aimag(str_a)+aimag(str_b))
      rel_diff = abs((aimag(str_b)-aimag(str_a))/mean)
      abs_diff = abs(aimag(str_b)-aimag(str_a))
      delta    = rel_diff*mean
      
      if (delta > tolerance_warn) then
         if (delta > tolerance_fail) then
            failed = .true.
            write(*,fmt) 'ERROR: Im[X] are inconsistent for', trim(label), rel_diff, abs_diff, delta
         else
            write(*,fmt) 'WARNING: Im[X] may be inconsistent for', trim(label), rel_diff, abs_diff, delta
         end if
      end if
      
      if (failed) stop
      
   end subroutine assert_equal_strength

end module betadecay
