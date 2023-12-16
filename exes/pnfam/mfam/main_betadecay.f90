!------------------------------------------------------------------------------
! pnfam_matrix/main_betadecay.f90
!
! T. Shafer, UNC Chapel Hill, 2015
!------------------------------------------------------------------------------
program pnfam_matrix_betadecay
   use logger
   use hfbtho_basis,      only : get_hfbtho_solution
   use hfb_solution_type, only : hfb_solution
   implicit none

   integer, parameter :: dp = kind(1.0d0)

   type matrix_element_set
      logical :: is_loaded
      integer :: k
      real(dp), allocatable :: energy(:), me(:,:)
   end type matrix_element_set

   type beta_config
      real(dp) :: ga_effective, phasespace_si
      logical  :: disable_axial_symmetry, disable_fermi
      character(len=250) :: output_file, hfb_solution_file
      character(len=250) :: file_00p, file_p1p, file_m1p, file_00m, file_p1m, file_m1m, file_p2m, file_m2m
   end type beta_config

   type(beta_config)        :: config
   type(hfb_solution)       :: hfb
   type(matrix_element_set) :: me_k00p, me_kp1p, me_k00m, me_kp1m, me_kp2m

   logical  :: sw_calc_fermi, sw_calc_gt_k0, sw_calc_gt_kp1, sw_calc_ff_k0,   &
      sw_calc_ff_kp1, sw_calc_ff_kp2

   real(dp) :: rate_k0_fermi, rate_k0_gt, rate_k1_gt, rate_k0_j0_ff, rate_k0_j1_ff, &
      rate_k0_j2_ff, rate_k1_j1_ff, rate_k1_j2_ff, rate_k2_j2_ff

   ! Global beta-decay constants
   integer  :: nnum, znum, anum, zfnum
   real(dp) :: eqrpa_max, r, lambda

   !----------------------------------------------------------------------------
   !----------------------------------------------------------------------------

   ! Logging
   openlog  => openlog_stdout
   writelog => writelog_stdout
   closelog => closelog_stdout
   call openlog('')

   ! Input namelist
   call read_namelist_and_init(config)

   ! Read in HFB data (Lambda values)
   call get_hfbtho_solution(fn=trim(config%hfb_solution_file), ep=hfb%ep, vp=hfb%vp, &
      up=hfb%up, en=hfb%en, vn=hfb%vn, un=hfb%un)

   ! Calculate beta-decay constants
   call calculate_constants(config=config, emax=eqrpa_max, radius=r, &
      xlam=lambda, n=nnum, z=znum, a=anum)
   zfnum = znum+1

   write(*,*)
   write(*,'(a,1f10.6," MeV")') 'EQRPA max ......... : ', eqrpa_max
   write(*,'(a,1f10.6)')        'Nuclear radius .... : ', r
   write(*,'(a,1f10.6)')        'Lambda (-gA/gV) ... : ', lambda

   ! Read in strength data
   ! Allowed
   me_k00p%k = 0
   if (len_trim(config%file_00p) > 0) then
      call load_matrix_element_data(file=trim(config%file_00p), ncols=5, me_set=me_k00p)
   end if
   me_kp1p%k = 1
   if (len_trim(config%file_p1p) > 0) then
      call load_matrix_element_data(file=trim(config%file_p1p), ncols=3, me_set=me_kp1p)
   end if

   ! Forbidden
   me_k00m%k = 0
   if (len_trim(config%file_00m) > 0) then
      call load_matrix_element_data(file=trim(config%file_00m), ncols=7, me_set=me_k00m)
   end if
   me_kp1m%k = 1
   if (len_trim(config%file_p1m) > 0) then
      call load_matrix_element_data(file=trim(config%file_p1m), ncols=5, me_set=me_kp1m)
   end if
   me_kp2m%k = 2
   if (len_trim(config%file_p2m) > 0) then
      call load_matrix_element_data(file=trim(config%file_p2m), ncols=3, me_set=me_kp2m)
   end if

   ! Main switches
   sw_calc_fermi   = ((.not.config%disable_fermi) .and. me_k00p%is_loaded)
   sw_calc_gt_k0   = me_k00p%is_loaded
   sw_calc_gt_kp1  = me_kp1p%is_loaded
   sw_calc_ff_k0   = me_k00m%is_loaded
   sw_calc_ff_kp1  = me_kp1m%is_loaded
   sw_calc_ff_kp2  = me_kp2m%is_loaded

   write(*,*)

   !----------------------------------------------------------------------------
   ! Actual beta-decay calculations
   !----------------------------------------------------------------------------
   ! Allowed
   rate_k0_fermi = 0;  rate_k0_gt = 0;  rate_k1_gt = 0

   if (sw_calc_fermi) then
      call allowed_rate(me_set=me_k00p, const=1.0_dp, col=1, emax=eqrpa_max, rate=rate_k0_fermi)
      write(*,'("0+ RATE, K=0  =",1f14.6)') rate_k0_fermi
   end if
   if (sw_calc_gt_k0) then
      call allowed_rate(me_set=me_k00p, const=1.0_dp, col=2, emax=eqrpa_max, rate=rate_k0_gt)
      write(*,'("1+ RATE, K=0  =",1f14.6)') rate_k0_gt
   end if
   if (sw_calc_gt_kp1) then
      call allowed_rate(me_set=me_kp1p, const=sqrt(2.0_dp), col=1, emax=eqrpa_max, rate=rate_k1_gt)
      write(*,'("1+ RATE, K=1  =",1f14.6)') rate_k1_gt
   end if

   ! Forbidden ...
   ! These implicitly use the globals we set before (r, lambda, eqrpa_max, etc.)
   ! even though they aren't declared.
   rate_k0_j0_ff = 0;  rate_k0_j1_ff = 0;  rate_k0_j2_ff = 0
   rate_k1_j1_ff = 0;  rate_k1_j2_ff = 0;  rate_k2_j2_ff = 0

   if (sw_calc_ff_k0) then
      call forbidden_rate(me_set=me_k00m, const=1.0_dp, &
         rate_j0=rate_k0_j0_ff, rate_j1=rate_k0_j1_ff, rate_j2=rate_k0_j2_ff)
      write(*,'("FF RATE, K=0  =",3f14.6)') rate_k0_j0_ff, rate_k0_j1_ff, rate_k0_j2_ff
   end if
   if (sw_calc_ff_kp1) then
      call forbidden_rate(me_set=me_kp1m, const=sqrt(2.0_dp), &
         rate_j1=rate_k1_j1_ff, rate_j2=rate_k1_j2_ff)
      write(*,'("FF RATE, K=0  =",3f14.6)') 0d0, rate_k1_j1_ff, rate_k1_j2_ff
   end if
   if (sw_calc_ff_kp2) then
      call forbidden_rate(me_set=me_kp2m, const=sqrt(2.0_dp), &
         rate_j2=rate_k2_j2_ff)
      write(*,'("FF RATE, K=2  =",3f14.6)') 0d0, 0d0, rate_k2_j2_ff
   end if

   ! Print result
   write(*,*)
   call show_results(store=.true.)

contains

   !----------------------------------------------------------------------------
   ! Initialize the calculation
   !----------------------------------------------------------------------------
   subroutine read_namelist_and_init(config)
      use pnfam_matrix, only : input_file_from_cli
      implicit none

      type(beta_config), intent(out) :: config

      real(dp) :: ga_effective
      logical  :: disable_axial_symmetry, disable_fermi
      character(len=250) :: output_file, hfb_solution_file
      namelist /parameters/ ga_effective, disable_axial_symmetry, disable_fermi, &
         output_file, hfb_solution_file

      real(dp) :: phasespace_si
      namelist /debug/ phasespace_si

      character(len=250) :: file_00p, file_p1p, file_m1p, file_00m, file_p1m, file_m1m, file_p2m, file_m2m
      namelist /input_files/ file_00p, file_p1p, file_m1p, file_00m, file_p1m, file_m1m, file_p2m, file_m2m

      character(len=250) :: infile
      integer :: ierr

      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------

      call input_file_from_cli(infile, default='matrix-beta.in')

      open(1, file=trim(infile), status='old', iostat=ierr)
      if (ierr /= 0) then
         write(*,'(2a)') "ERROR: could not open the file ", trim(infile)
         stop
      end if

      ga_effective = 1.0
      disable_axial_symmetry = .false.
      disable_fermi = .false.
      output_file = ""
      hfb_solution_file = ""
      read(1, nml=parameters)
      config%ga_effective           = ga_effective
      config%disable_axial_symmetry = disable_axial_symmetry
      config%disable_fermi          = disable_fermi
      config%output_file            = trim(output_file)
      config%hfb_solution_file      = trim(hfb_solution_file)


      phasespace_si = 1.0d-5
      read(1, nml=debug)
      config%phasespace_si = phasespace_si

      file_00p = "";  file_p1p = "";  file_m1p = "";  file_00m = "";  file_p1m = ""
      file_m1m = "";  file_p2m = "";  file_m2m = "";
      read(1, nml=input_files)
      config%file_00p = trim(file_00p)
      config%file_p1p = trim(file_p1p)
      config%file_m1p = trim(file_m1p)
      config%file_00m = trim(file_00m)
      config%file_p1m = trim(file_p1m)
      config%file_m1m = trim(file_m1m)
      config%file_p2m = trim(file_p2m)
      config%file_m2m = trim(file_m2m)
      close(1)

      ! Display
      write(*,'(a)') 'MATRIX-FAM BETA-DECAY CODE'
      write(*,*)
      write(*,'(a)') "INPUT PARAMETERS:"
      write(*,'(a)') repeat('-', 33)
      write(*,'(a,1f4.1)')  "gA ....................... :", config%ga_effective
      write(*,'(a,a)')      "Disable axial symmetry ... : ", trim(merge('Yes', 'No ', config%disable_axial_symmetry))
      write(*,'(a,a)')      "Disable Fermi decay ...... : ", trim(merge('Yes', 'No ', config%disable_fermi))
      write(*,'(a,a)')      "HFB solution file ........ : ", trim(config%hfb_solution_file)
      write(*,'(a,a)')      "Output file .............. : ", trim(config%output_file)
      write(*,*)
      write(*,'(a,1es7.1)') "Phase space SI ........... : ", config%phasespace_si
      write(*,*)
      write(*,'(a,a)')      "File  0+ ................. : ", trim(config%file_00p)
      write(*,'(a,a)')      "File  1+ ................. : ", trim(config%file_p1p)
      write(*,'(a,a)')      "File -1+ ................. : ", trim(config%file_m1p)
      write(*,'(a,a)')      "File  0- ................. : ", trim(config%file_00m)
      write(*,'(a,a)')      "File  1- ................. : ", trim(config%file_p1m)
      write(*,'(a,a)')      "File -1- ................. : ", trim(config%file_m1m)
      write(*,'(a,a)')      "File  2- ................. : ", trim(config%file_p2m)
      write(*,'(a,a)')      "File -2- ................. : ", trim(config%file_m2m)

   end subroutine read_namelist_and_init

   !----------------------------------------------------------------------------
   ! Calculate a few beta-decay constants for later
   !----------------------------------------------------------------------------
   subroutine calculate_constants(config, emax, radius, xlam, n, z, a)
      use constants,    only : dmassnH, r0, hbar_mec
      use hfbtho_basis, only : hfb_npr, hfb_lambda
      implicit none

      type(beta_config), intent(in)  :: config
      integer,           intent(out) :: n, z, a
      real(dp),          intent(out) :: emax, radius, xlam

      ! Particle number
      n = hfb_npr(1)
      z = hfb_npr(2)
      a = hfb_npr(3)

      ! Max. EQRPA
      emax = hfb_lambda(1)-hfb_lambda(2)+dmassnH

      ! Nuclear radius (natural units)
      radius = r0*hfb_npr(3)**(1.0_dp/3.0_dp)/hbar_mec

      ! Lambda = -gA/gV must be >= 0!
      ! Note, though, we enter gA as a positive number
      xlam = config%ga_effective

   end subroutine calculate_constants

   !----------------------------------------------------------------------------
   ! Read a strength file, expected to have a certain number of columns
   ! including EQRPA.
   !----------------------------------------------------------------------------
   subroutine load_matrix_element_data(file, ncols, me_set)
      implicit none

      character(len=*),         intent(in)    :: file
      integer,                  intent(in)    :: ncols
      type(matrix_element_set), intent(inout) :: me_set

      integer :: ierr, nrows, i, j

      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------

      me_set%is_loaded = .false.

      open(2, file=trim(file), status='old', form='formatted', iostat=ierr)
      if (ierr /= 0) then
         write(*,'(3a)') 'WARNING: could not open ', file, ' for reading.'
         return
      end if

      if (allocated(me_set%energy)) deallocate(me_set%energy)
      if (allocated(me_set%me))     deallocate(me_set%me)

      ! Count lines
      nrows = 0
      do
         read(2, *, iostat=ierr)
         if (ierr /= 0) then
            exit
         else
            nrows = nrows+1
         end if
      end do

      ! Column 0 is EQRPA, so me_set has N-1 columns
      allocate(me_set%energy(nrows), me_set%me(nrows, ncols-1))
      me_set%energy(:) = 0;  me_set%me(:,:) = 0

      rewind(2)
      do i=1, nrows
         read(2,*) me_set%energy(i), (me_set%me(i,j), j=1, ncols-1)
      end do
      close(2)

      me_set%is_loaded = .true.

   end subroutine load_matrix_element_data

   !----------------------------------------------------------------------------
   ! Calculate an allowed-shape beta-decay rate without charge screening via
   ! rate = ln(2)/kappa * B * f0(E).
   !----------------------------------------------------------------------------
   subroutine allowed_rate(me_set, const, col, emax, rate)
      use constants,    only : mec2, ln2_kappa, kappa, ln2
      use phasespace,   only : real_psi
      use hfbtho_basis, only : hfb_npr
      implicit none

      type(matrix_element_set), intent(in) :: me_set
      integer,  intent(in)  :: col
      real(dp), intent(in)  :: const, emax
      real(dp), intent(out) :: rate

      integer :: i
      real(dp) :: me, b, eqrpa, f, r, logft

      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------

      rate = 0
      if (.not.me_set%is_loaded) return

      !write(*,*)
      !write(*,'("Using column ", 1i0, " with K constant ", 1f0.4)') col, const**2
      !write(*,'("Maximum EQRPA = ",1f0.4," MeV")') emax
      !write(*,'(a)') repeat('-',80)
      !write(*,'(a)') "       I      EQRPA         B(GT-)        RATE          F0              LOGFT"
      !write(*,'(a)') repeat('-',80)

      do i=1, size(me_set%energy)
         eqrpa = me_set%energy(i)
         if (eqrpa > emax) exit

         ! Transition strength is ME**2 * CONST**2
         me = me_set%me(i, col)
         b  = (const*me)**2

         ! Integral of the usual Fermi function
         f = real_psi(nz=hfb_npr(2)+1, na=hfb_npr(3), w0=1+(emax-eqrpa)/mec2, n=2, si=config%phasespace_si)

         ! Single decay rate, total decay rate, and log(ft)
         r = ln2_kappa*b*f
         rate  = rate+r
         logft = log10(kappa/b)

         ! Log the result
         ! if (logft <= 10 .and. r > 1d-3) then
         !    write(*,'(1i8,3f14.6,1es16.4,1f14.6)') i, eqrpa, b, r, f, log10(kappa/b)
         ! end if
      end do

   end subroutine allowed_rate

   !----------------------------------------------------------------------------
   ! Forbidden rate
   ! Uses globals implicitly ... eqrpa_max, particle number, lambda, r, etc.
   !----------------------------------------------------------------------------
   subroutine forbidden_rate(me_set, const, rate_j0, rate_j1, rate_j2)
      use constants,    only : hbar_mec, mec2, Mn, alpha, ln2_kappa
      use phasespace,   only : real_psi
      use fermi,        only : gamma_1
      implicit none

      type(matrix_element_set), intent(in)  :: me_set
      real(dp),                 intent(in)  :: const
      real(dp), optional,       intent(out) :: rate_j0, rate_j1, rate_j2

      integer  :: k, i
      real(dp) :: crs0, crs1, crs2, cr, cp, cps0
      real(dp) :: tmp_rate_j0, tmp_rate_j1, tmp_rate_j2
      real(dp) :: eqrpa, wmax, xm, xp
      real(dp) :: brs0, brs1, brs2, bps0, br, bp, mrs0ps0, mrrs1, mprs1, mpr
      real(dp) :: fi1, fi2, fi3, fi4, fi5, fi6
      real(dp) :: c0, c1, c2

      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------

      ! K value is determined from the matrix elements we have passed
      k = abs(me_set%k)

      ! Constants of integration
      crs0 =  1.0_dp / hbar_mec
      crs1 = -const  / hbar_mec
      crs2 =  const  / hbar_mec
      cr   =  const  * sqrt(3.0) / hbar_mec
      cps0 = -1.0_dp * hbar_mec  * mec2 / Mn
      cp   = -const  * hbar_mec  * mec2 / Mn

      ! Individual J rates
      tmp_rate_j0 = 0;  tmp_rate_j1 = 0;  tmp_rate_j2 = 0

      ! Loop over energetically-accessible states
      do i=1, size(me_set%energy)
         eqrpa = me_set%energy(i)
         if (eqrpa > eqrpa_max) exit

         ! Beta-decay quantities
         wmax = 1.0_dp + (eqrpa_max-eqrpa)/mec2
         xp = wmax/3.0_dp + 0.5_dp*alpha*zfnum/r
         xm = wmax/3.0_dp - 0.5_dp*alpha*zfnum/r

         ! Strength functions and cross-terms
         brs0 = 0;  bps0    = 0;  brs1  = 0;  br    = 0;  bp  = 0
         brs2 = 0;  mrs0ps0 = 0;  mrrs1 = 0;  mprs1 = 0;  mpr = 0
         select case (k)
            case (0)
               brs0    = crs0**2     * me_set%me(i,1)**2
               bps0    = cps0**2     * me_set%me(i,2)**2
               brs1    = crs1**2     * me_set%me(i,3)**2
               br      = cr**2       * me_set%me(i,4)**2
               bp      = cp**2       * me_set%me(i,5)**2
               brs2    = crs2**2     * me_set%me(i,6)**2
               mrs0ps0 = crs0 * cps0 * me_set%me(i,1)*me_set%me(i,2)
               mrrs1   = cr   * crs1 * me_set%me(i,4)*me_set%me(i,3)
               mpr     = cr   * cp   * me_set%me(i,5)*me_set%me(i,4)
               mprs1   = cp   * crs1 * me_set%me(i,5)*me_set%me(i,3)
            case (1, -1)
               brs1  = crs1**2   * me_set%me(i,1)**2
               br    = cr**2     * me_set%me(i,2)**2
               bp    = cp**2     * me_set%me(i,3)**2
               brs2  = crs2**2   * me_set%me(i,4)**2
               mrrs1 = cr * crs1 * me_set%me(i,2)*me_set%me(i,1)
               mpr   = cr * cp   * me_set%me(i,2)*me_set%me(i,3)
               mprs1 = cp * crs1 * me_set%me(i,3)*me_set%me(i,1)
            case (2, -2)
               brs2 = crs2**2 * me_set%me(i,1)**2
         end select

         ! Phase space integration
         fi1 = real_psi(nz=zfnum, na=anum, w0=wmax, n=1, si=config%phasespace_si)
         fi2 = real_psi(nz=zfnum, na=anum, w0=wmax, n=2, si=config%phasespace_si)
         fi3 = real_psi(nz=zfnum, na=anum, w0=wmax, n=3, si=config%phasespace_si)
         fi4 = real_psi(nz=zfnum, na=anum, w0=wmax, n=4, si=config%phasespace_si)
         fi5 = real_psi(nz=zfnum, na=anum, w0=wmax, n=5, si=config%phasespace_si)
         fi6 = real_psi(nz=zfnum, na=anum, w0=wmax, n=6, si=config%phasespace_si)

         ! Complete shape factor
         c0 = 0;  c1 = 0;  c2 = 0

         ! J = 0
         c0 = c0 + lambda**2*(bps0 + brs0*(1/9.0_dp + xp**2) + mrs0ps0*2*xp)*fi2
         c0 = c0 - 2/3.0_dp*lambda**2*(brs0*xp + mrs0ps0)*fi1
         ! J = 1
         c1 = c1 - 2/9.0_dp*(xp*br - 2*lambda**2*xm*brs1    &
                 + lambda*sqrt(2.0_dp)*(xp-xm)*mrrs1        &
                 - sqrt(3.0_dp)*mpr - lambda*sqrt(6.0_dp)*mprs1)*fi1
         c1 = c1 + (bp + 1/3.0_dp*xp**2*br + 2/3.0_dp*xm**2*lambda**2*brs1)*fi2
         c1 = c1 + sqrt(2.0_dp/3.0_dp)*(2*lambda*xm*mprs1 - sqrt(2.0_dp)*xp*mpr &
                 - 2/sqrt(3.0_dp)*lambda*xm*xp*mrrs1)*fi2
         c1 = c1 - 8/27.0_dp*(lambda**2*brs1 + lambda/sqrt(2.0_dp)*mrrs1)*gamma_1(zfnum)*fi2
         c1 = c1 + 1/27.0_dp*(br + 2*lambda**2*brs1 + 2*sqrt(2.0_dp)*lambda*mrrs1)*fi2
         c1 = c1 + (4*sqrt(2.0_dp)/9.0_dp*lambda*xp*mrrs1 - 8/9.0_dp*lambda**2*xm*brs1 &
                 -  4*sqrt(2.0_dp)/3.0_dp/sqrt(3.0_dp)*lambda*mprs1)*fi3
         c1 = c1 + 8/27.0_dp*lambda**2*brs1*fi4
         c1 = c1 + 1/27.0_dp*(2*br + brs1*lambda**2 + 2*sqrt(2.0_dp)*mrrs1*lambda)*fi5
         c1 = c1 + 1/27.0_dp*(2*br + brs1*lambda**2 - 2*sqrt(2.0_dp)*mrrs1*lambda)*fi6
         ! J = 2
         c2 = c2 + 1/9.0_dp*lambda**2*brs2*(fi5+fi6)

         ! Storage
         tmp_rate_j0 = tmp_rate_j0+c0
         tmp_rate_j1 = tmp_rate_j1+c1
         tmp_rate_j2 = tmp_rate_j2+c2
      end do

      ! Convert to decay rates and store
      tmp_rate_j0 = ln2_kappa*tmp_rate_j0
      tmp_rate_j1 = ln2_kappa*tmp_rate_j1
      tmp_rate_j2 = ln2_kappa*tmp_rate_j2

      if (.not.present(rate_j0) .and. abs(tmp_rate_j0) > 1e-16) then
         write(*,'("Warning: J=0 not requested but non-zero [K=", 1i0, "]: ", 1es9.1)') k, tmp_rate_j0
      end if
      if (.not.present(rate_j1) .and. abs(tmp_rate_j1) > 1e-16) then
         write(*,'("Warning: J=1 not requested but non-zero [K=", 1i0, "]: ", 1es9.1)') k, tmp_rate_j0
      end if
      if (.not.present(rate_j2) .and. abs(tmp_rate_j2) > 1e-16) then
         write(*,'("Warning: J=2 not requested but non-zero [K=", 1i0, "]: ", 1es9.1)') k, tmp_rate_j0
      end if

      if (present(rate_j0)) rate_j0 = tmp_rate_j0
      if (present(rate_j1)) rate_j1 = tmp_rate_j1
      if (present(rate_j2)) rate_j2 = tmp_rate_j2

   end subroutine forbidden_rate

   !----------------------------------------------------------------------------
   ! Copied from the i-FAM betadecay file.
   ! Uses globals quite heavily. :/
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
      real(dp) :: rate_k0_ff, rate_k1_ff, rate_k2_ff, rate_j0_ff, rate_j1_ff, rate_j2_ff

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
      rate0p = rate_k0_fermi + rate_k0_gt
      rate1p = rate_k1_gt
      rate0m = rate_k0_j0_ff+rate_k0_j1_ff+rate_k0_j2_ff
      rate1m = rate_k1_j1_ff+rate_k1_j2_ff
      rate2m = rate_k2_j2_ff

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
      ratej0p = rate_k0_fermi
      ratej1p = rate_k0_gt + rate_k1_gt
      ratej0m = rate_k0_j0_ff
      ratej1m = rate_k0_j1_ff+rate_k1_j1_ff
      ratej2m = rate_k0_j2_ff+rate_k1_j2_ff+rate_k2_j2_ff
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
         rate_k0_ff = rate_k0_j0_ff+rate_k0_j1_ff+rate_k0_j2_ff
         rate_k1_ff = rate_k1_j1_ff+rate_k1_j2_ff
         rate_k2_ff = rate_k2_j2_ff
         rate_j0_ff = rate_k0_j0_ff
         rate_j1_ff = rate_k0_j1_ff+rate_k1_j1_ff
         rate_j2_ff = rate_k0_j2_ff+rate_k1_j2_ff+rate_k2_j2_ff

         open(10, file=trim(config%output_file), status='unknown', iostat=ierr)

         ! Header information
         write(10,'(a)')  '# Nuclear Beta Decay Rates and Half-Lives'
         write(10,'(a)')  '# MATRIX-FAM CODE'

         call date_and_time(values=date_values)
         write(10,'(a,1i0,"/",1i0,"/",1i0,1x,1i0,":",1i2.2)') '# Run date: ', date_values(2), &
            date_values(3), date_values(1), date_values(5), date_values(6)

         !!write(10,'(a)') '# '
         !!write(10,'(a,  i10)',advance='no') '# Z (ini) =', zinum
         !!write(10,'(a,  i10)',advance='no') ',     A       =', anum
         !!write(10,'(a,  i10)',advance='no') ',     Z (fin) =', zfnum
         !!write(10,'(a)')
         !!write(10,'(a,f10.5)',advance='no') '# GS en.  =', gsenergy
         !!write(10,'(a,f10.5)',advance='no') ',     Q-value =', qvalue
         !!write(10,'(a,f10.5)',advance='no') ',     EQRPAmax=', eqrpamax
         !!write(10,'(a,f10.5)',advance='no') ',     lambda  =', lambda
         !!write(10,'(a)')
         !!write(10,'(a,f10.5)',advance='no') '# g_A     =', ga_input
         !!write(10,'(a,f10.5)',advance='no') ',     g_V     =', gv
         !!write(10,'(a,f10.5)',advance='no') ',     Nucl. M =', mn
         !!write(10,'(a,f10.5)',advance='no') ',     alpha*Z =', alpha*zfnum
         !!write(10,'(a)')
         !!write(10,'(a,f10.5)',advance='no') '# R       =', r
         !!write(10,'(a,f10.5)',advance='no') ',     Xi      =', alpha*zfnum/2/r
         !!write(10,'(a,f10.5)',advance='no') ',     W0 max  =', 1+qvalue/mec2
         !!write(10,'(a,f10.5)',advance='no') ',     W0*R    =', (1+qvalue/mec2)*r
         !!write(10,'(a)')
         !!write(10,'(a)') '# '
         write(10,'(a)') '# FILE FORMAT:  [RATE] (s^-1)    [HALF-LIFE] (s)'

         ! Rates
         write(10,'(2es26.16,4x,1a)') rate_total,     ln2/rate_total,     '# Total'
         write(10,'(2es26.16,4x,1a)') rate_allowed,   ln2/rate_allowed,   '# Allowed'
         write(10,'(2es26.16,4x,1a)') rate_forbidden, ln2/rate_forbidden, '# Forbidden'
         write(10,'(2es26.16,4x,1a)') rate_k0_fermi,  ln2/rate_k0_fermi,  '# Fermi'
         write(10,'(2es26.16,4x,1a)') rate_k0_gt,     ln2/rate_k0_gt,     '# Gamow-Teller, K=0'
         write(10,'(2es26.16,4x,1a)') rate_k1_gt,     ln2/rate_k1_gt,     '# Gamow-Teller, K=1'
         write(10,'(2es26.16,4x,1a)') rate_k0_ff,     ln2/rate_k0_ff,     '# Forbidden, K=0'
         write(10,'(2es26.16,4x,1a)') rate_k1_ff,     ln2/rate_k1_ff,     '# Forbidden, K=1'
         write(10,'(2es26.16,4x,1a)') rate_k2_ff,     ln2/rate_k2_ff,     '# Forbidden, K=2'
         write(10,'(2es26.16,4x,1a)') rate_j0_ff,     ln2/rate_j0_ff,     '# Forbidden, J=0'
         write(10,'(2es26.16,4x,1a)') rate_j1_ff,     ln2/rate_j1_ff,     '# Forbidden, J=1'
         write(10,'(2es26.16,4x,1a)') rate_j2_ff,     ln2/rate_j2_ff,     '# Forbidden, J=2'
         write(10,'(2es26.16,4x,1a)') rate_k0_j0_ff,  ln2/rate_k0_j0_ff,  '# Forbidden, (J,K)=(0,0)'
         write(10,'(2es26.16,4x,1a)') rate_k0_j1_ff,  ln2/rate_k0_j1_ff,  '# Forbidden, (J,K)=(1,0)'
         write(10,'(2es26.16,4x,1a)') rate_k1_j1_ff,  ln2/rate_k1_j1_ff,  '# Forbidden, (J,K)=(1,1)'
         write(10,'(2es26.16,4x,1a)') rate_k0_j2_ff,  ln2/rate_k0_j2_ff,  '# Forbidden, (J,K)=(2,0)'
         write(10,'(2es26.16,4x,1a)') rate_k1_j2_ff,  ln2/rate_k1_j2_ff,  '# Forbidden, (J,K)=(2,1)'
         write(10,'(2es26.16,4x,1a)') rate_k2_j2_ff,  ln2/rate_k2_j2_ff,  '# Forbidden, (J,K)=(2,2)'

         close(10)
      end if
   end subroutine show_results

end program pnfam_matrix_betadecay
