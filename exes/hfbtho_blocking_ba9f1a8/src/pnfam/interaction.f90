!------------------------------------------------------------------------------
! interaction.f90
!
! Contains the coupling constants and the subroutines for computing them.
!
! M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
!------------------------------------------------------------------------------
module interaction
   use logger
   implicit none
   private

   integer,  parameter :: dp = kind(1.0d0)
   real(dp), parameter :: sc_tolerance = 5.0d-12
   
   ! Public coupling constants
   real(dp), allocatable, public :: crho(:), cs(:)
   real(dp), allocatable, public :: cpair(:), cspair(:)
   real(dp), public  :: cdrho, ctau, ctj0, ctj1, ctj2, crdj, cds, ct, cj, cgs, cf, csdj

   ! Private constants which make up Crho, Cs, and Cpair
   real(dp) :: cr0, crr, cs0, csr, sigma_r, sigma_s
   real(dp) :: cpair0, cpairr, sigma_pair, cspair0, cspairr
   
   ! Effective mass is a function of the Skyrme parameters
   real(dp) :: effmass
   
   ! This is the mean density that enters the pairing calculation in HFBTHO
   real(dp), parameter :: rho_c = 0.160_dp

   ! Public interaction parameters
   character(len=80), public :: interaction_name
   logical,           public :: skip_residual_interaction

   ! Private interaction parameters
   ! Namelist options:
   ! 1. force_j2_terms allows J**2 terms even if the Skyrme interaction
   !    wouldn't usually include them in the functional. E.g.
   !    force_j2_terms = .true. will make SkM* keep the J**2 terms.
   ! 2. require_gauge_invariance requires the functional satisfy the 5 gauge
   !    invariance properties. These should be satisfied automatically by the
   !    Skyrme interaction, so it really tests externally-loaded parameter sets.
   !    See E. Perlinska et al., Phys. Rev. C 69, 014316 (2004) for the gauge
   !    invariance.
   ! 3. require_self_consistency REQUIRES the FAM coupling constants to be
   !    within some tolerance of the HFBTHO coupling constants from the parent
   !    ground state calculation --- otherwise, the program halts. We are also
   !    aware of the HFBTHO flag "use_j2terms" and check against it. Note:
   !    ctj1(FAM) is (and I think should be) equal to ctj(HFBTHO)/2. See E.
   !    Perlinska et al., Phys. Rev. C 69, 014316 (2004).
   ! 4. vpair0, vpair1 are the original PNFAM pairing strengths (cf. HFBTHO). We
   !    don't allow altering the form of the pairing functional (e.g.
   !    surface/mixed/volume pairing) from the default namelist for reasons of
   !    self-consistency.  NB! It seems vpair0 is off by a factor of 2 --
   !    hence the introduction (6) of vpair_t0 and vpair_t1 to allow old
   !    calculations to still function.
   ! 5. the override_* values are used to override built-in time-odd constants
   !    which come with Skyrme-derived EDFs by comparing their values to their
   !    default values [huge(1d0)].
   ! 6. new pairing strengths vpair_t0 and vpair_t1 which should correctly
   !    mirror the HFBTHO pairing interaction. Still some question about
   !    factors of 2 related to kappa, but these seem OK.
   logical  :: force_j2_terms
   logical  :: require_gauge_invariance
   logical  :: require_self_consistency
   real(dp) :: vpair0, vpair1, vpair_t0, vpair_t1
   real(dp) :: override_cs0, override_csr, override_cds, override_ct, override_cf, &
      override_cgs, override_cj, override_csdj

   ! External interaction
   logical :: use_external_interaction
   character(len=200) :: external_interaction_file

   ! This is used by the subroutine print_interaction_params
   logical :: skyrme_uses_j2_terms

   ! Publicly avaliable subroutines
   public :: init_skyrme, print_interaction_params

contains

   !---------------------------------------------------------------------------
   ! Use the input namelist to define the residual interaction.
   !---------------------------------------------------------------------------
   subroutine init_skyrme(fh, vp, vn, up, un)
      use constants,    only : translate_uppercase
      use hfbtho_basis, only : nghl, hfb_cpair, hfb_alpha_pair, hfb_density
      implicit none

      integer,  intent(in) :: fh
      real(dp), intent(in), dimension(:) :: vp, vn, up, un

      integer :: ierr
      real(dp), dimension(:), allocatable :: rho
      
      character(len=80)  :: tmpname
      character(len=160) :: st

      ! Interaction namelist definition
      namelist /interaction/ interaction_name, force_j2_terms,             &
               require_gauge_invariance, require_self_consistency,         &
               vpair0, vpair1, vpair_t0, vpair_t1, override_cs0,           &
               override_csr, override_cds,  override_ct, override_cf,      &
               override_cgs, override_cj, override_csdj

      ! Allocate & zero the coupling constants
      call init_coupling_constants
      
      ! Namelist defaults are set here
      interaction_name = 'NONE'
      force_j2_terms = .false.
      require_gauge_invariance = .false.
      require_self_consistency = .true.
      ! Pairing strengths
      vpair0   = -huge(1.0d0)
      vpair1   = -huge(1.0d0)
      vpair_t0 = -huge(1.0d0)
      vpair_t1 = -huge(1.0d0)
      ! Override couplings (e.g. to customize a Skyrme interaction's GTR)
      override_cs0 = huge(1.0d0)
      override_csr = huge(1.0d0)
      override_cds = huge(1.0d0)
      override_ct  = huge(1.0d0)
      override_cf  = huge(1.0d0)
      override_cgs = huge(1.0d0)
      override_cj = huge(1.0d0)
      override_csdj = huge(1.0d0)

      rewind(fh)
      read(fh, nml=interaction, iostat=ierr)

      if (ierr /=0) then
         write(st,'(a,i0)') 'ERROR: could not read namelist INTERACTION: ', ierr
         call writelog(st)
         call closelog
         stop
      end if

      interaction_name = adjustl(interaction_name)
      tmpname = interaction_name
      call translate_uppercase(tmpname)
      
      use_external_interaction = .false.
      
      !-------------------------------------------------------------------------
      ! If the interaction is "NONE", set the flag to skip the residual
      ! interaction and return. Allow some variation of case.
      ! tmpname = uppercase(interaction_name)
      !-------------------------------------------------------------------------
      if (trim(tmpname) == 'NONE') then
         skip_residual_interaction = .true.
         return
      !-------------------------------------------------------------------------
      ! If there is an interaction, load_ph_interaction() handles the ph part.
      ! The pp-part is set below as well. 
      !-------------------------------------------------------------------------
      else
         skip_residual_interaction = .false.
         call load_ph_interaction(skyrme_uses_j2_terms, use_external_interaction)

         !----------------------------------------------------------------------
         ! The following happens only if the user chooses a built-in functional.
         ! If they have used custom_interaction.dat, none of this applies; we
         ! assume they have set the J^2 terms/gauge invariance/pairing
         ! appropriately there.
         !----------------------------------------------------------------------
         if (use_external_interaction .eqv. .false.) then
            ! If we've used an included Skyrme functional, erase the J**2 terms
            ! if required.
            if ((skyrme_uses_j2_terms .eqv. .false.) .and. (force_j2_terms .eqv. .false.)) then
               ctj0 = 0;  ctj1 = 0;  ctj2 = 0
            end if
            
            ! The coupling constants CT and CF are linked to the J**2
            ! couplings via gauge-invariance.  The CJ are fixed by HFBTHO.
            if (require_gauge_invariance .eqv. .true.) then
               ! CF is uniquely fixed by the difference between CJ(1) and CJ(2)
               cf = 2*ctj1-ctj2
               ! CT is then fixed 
               ct = -2*ctj1 + 0.5_dp*cf
            end if
            
            ! Set the pairing via the namelist parameters. Relying on absurd
            ! default values for the pairing strengths & floating-point
            ! comparisons. Should be OK since we define vpair* = -huge(1d0).
            
            ! T=0 pairing
            ! NB! Now using the (correct I think) vpair_t0 and warning if the
            ! user passes the old pairing strengths (but the code will still run).
            if (abs(vpair_t0+huge(1.0d0)) <= epsilon(1d0)) then
               ! Check for legacy pairing
               if (abs(vpair0+huge(1.0d0)) <= epsilon(1d0)) then
                  ! If not set, default to 0 isoscalar pairing
                  vpair_t0 = 0
               else
                  ! Translate
                  vpair_t0 = 2*vpair0
                  ! Warn the user
                  write(st,'(a)') 'Warning: input parameter vpair0 is wrong by a factor of 2&
                     & and deprecated --- use vpair_t0 instead'
                  call writelog(st)
               end if
            else
               ! In this case, both vpair0 AND vpair_t0 are set, which is an error
               if (abs(vpair0+huge(1.0d0)) > epsilon(1d0)) then
                  write(st,'(a,i0)') 'ERROR: cannot set both vpair0 and vpair_t0!'
                  call writelog(st)
                  call closelog
                  stop
               end if
            end if
            
            ! T=1 pairing
            if (abs(vpair_t1+huge(1.0d0)) <= epsilon(1d0)) then
               ! Check for legacy pairing
               if (abs(vpair1+huge(1.0d0)) <= epsilon(1d0)) then
                  ! If not set, default to the HFBTHO average value
                  vpair_t1 = 0.5_dp*sum(hfb_cpair)
               else
                  ! Translate
                  vpair_t1 = vpair1
                  ! Warn the user
                  write(st,'(a)') 'Warning: input parameter vpair1 is deprecated ---&
                      & use vpair_t1 instead'
                  call writelog(st)
               end if
            else
               ! In this case, both vpair0 AND vpair_t0 are set, which is an error
               if (abs(vpair1+huge(1.0d0)) > epsilon(1d0)) then                  
                  write(st,'(a,i0)') 'ERROR: cannot set both vpair1 and vpair_t1!'
                  call writelog(st)
                  call closelog
                  stop
               end if
            end if
            
            ! Assemble the pairing coupling constants using the same functional form as HFBTHO.
            ! Sigma pair should really always be 1 from HFBTHO, but allow modifications
            if (abs(sigma_pair) < 1d0*epsilon(1.0)) then
               sigma_pair = 1.0_dp
            end if
            
            ! cpair0  - isovector pairing strength
            ! cspair0 - (negative) isoscalar pairing strength
            cpair0  =    vpair_t1/8.0_dp ! cpair0  =  v1
            cspair0 =    vpair_t0/8.0_dp ! cspair0 =  v0
            
            ! The density-dependence is controlled by HFBTHO
            cpairr  =  cpair0  * (-1*hfb_alpha_pair(1)/rho_c**sigma_pair)
            cspairr =  cspair0 * (-1*hfb_alpha_pair(1)/rho_c**sigma_pair)
            
            ! If the custom couplings are set, update them here
            if (abs(override_cs0-huge(1d0)) > epsilon(1d0)) then
               cs0 = override_cs0
            end if
            if (abs(override_csr-huge(1d0)) > epsilon(1d0)) then
               csr = override_csr
            end if
            if (abs(override_cds-huge(1d0)) > epsilon(1d0)) then
               cds = override_cds
            end if
            if (abs(override_ct-huge(1d0)) > epsilon(1d0)) then
               ct = override_ct
            end if
            if (abs(override_cf-huge(1d0)) > epsilon(1d0)) then
               cf = override_cf
            end if
            if (abs(override_cgs-huge(1d0)) > epsilon(1d0)) then
               cgs = override_cgs
            end if
            if (abs(override_cj-huge(1d0)) > epsilon(1d0)) then
               cj = override_cj
            end if
            if (abs(override_csdj-huge(1d0)) > epsilon(1d0)) then
               csdj = override_csdj
            end if
         end if
         
         ! Tests to run
         if (require_self_consistency .eqv. .true.) then
            call assert_self_consistency(tolerance=sc_tolerance)
         end if

         if ((require_gauge_invariance .eqv. .true.)  .and. &
             (gauge_invariance_is_ok() .eqv. .false.)) then
            call print_interaction_params
            call writelog("")
            call writelog('ERROR: gauge invariance was requested, but is not fulfilled.')
            call closelog
            stop
         end if

         ! Allocate and calculate the isoscalar density (it=3 => rho_n + rho_p)
         if (allocated(rho) .eqv. .false.) allocate(rho(nghl))
         call hfb_density(it=3, vp=vp, vn=vn, up=up, un=un, rho=rho)

         ! Put together the density-dependent coupling constants
         ! ph channel
         crho(:) = cr0 + crr*rho(:)**sigma_r
         cs(:)   = cs0 + csr*rho(:)**sigma_s
         
         ! pp channel:
         ! cpair[0r]  - isovector pairing strength
         ! cspair[0r] - (negative) isoscalar pairing strength
         cpair(:)  = cpair0  + cpairr*rho(:)**sigma_pair
         cspair(:) = cspair0 + cspairr*rho(:)**sigma_pair
         return
      end if
   end subroutine init_skyrme


   !---------------------------------------------------------------------------
   ! Given the global "interaction_name", load t's and x's. Then convert them
   ! to isovector coupling constants.
   !---------------------------------------------------------------------------
   subroutine load_ph_interaction(skyrme_uses_j2_terms, ext_interaction)
      use constants, only : translate_uppercase
      implicit none
      
      logical, intent(out) :: skyrme_uses_j2_terms, ext_interaction
      
      ! Use the most general tensor terms (i.e. CJ0, CJ1, CJ2 not proportional).
      ! In general, HFBTHO does *not* support this
      logical :: skyrme_uses_extended_tensor_terms

      ! Skyrme parameters
      real(dp) :: t0, t1, t2, t3, x0, x1, x2, x3, b4, b4p, w, alpha, tto, tte
      
      ! Auxiliary variables
      character(len=80)  :: interaction_name_ucase
      character(len=160) :: st

      ! Namelist for reading the couplings from external file. This has ZERO
      ! error checking, and should maybe be cleaned up later.
      integer :: ierr
      namelist /custom_interaction/ cr0, crr, sigma_r, cs0, csr, sigma_s,  &
               cpair0, cpairr, cspair0, cspairr, sigma_pair, cdrho, ctau,  &
               ctj0, ctj1, ctj2, crdj, cds, ct, cj, csdj, cf, cgs

      ! By default, use J**2 terms
      skyrme_uses_j2_terms = .true.
      
      ! By default, use the simplified J^2 terms (single coupling CJ)
      skyrme_uses_extended_tensor_terms = .false.
      
      ! Conversion
      interaction_name_ucase = interaction_name
      call translate_uppercase(interaction_name_ucase)

      ! Set up the t's and x's
      ! If loading from an external file, read in the coupling constants
      if (interaction_name_ucase(1:5) == 'FILE:') then

         ext_interaction = .true.
         interaction_name(1:5) = 'FILE:'
         external_interaction_file = adjustl(interaction_name(6:))

         ! Read the contents
         open(unit=45, file=trim(external_interaction_file), status='old', iostat=ierr)

         if (ierr /= 0) then
            write(st,'(3a,i0)') 'ERROR: could not open interaction file "',   &
               trim(external_interaction_file), '".  Error number ', ierr
            call writelog(st)
            call closelog
            stop
         end if

         read(45, nml=custom_interaction, iostat=ierr)

         if (ierr /= 0) then
            write(st,'(3a,i0)') 'ERROR: could not read interaction file "',   &
               trim(external_interaction_file), '".  Error number ', ierr
            call writelog(st)
            call closelog
            stop
         end if

         close(45)
         return
      
      ! If not loading from an external file, try an included Skyrme functional.
      else
         ! No problem with just uppercasing the whole thing
         interaction_name = interaction_name_ucase
         
         ! Listing of included functionals
         select case (trim(interaction_name))
            ! SIII
            ! M. Beiner et al., Nuclear Physics A 238, 29 (1975).
            ! NOTE: This reference doesn't seem to include alpha or x3. As a
            ! result, I fall back to the HFBTHOv2 parameters.
            case('SIII', 'S3')
               interaction_name = 'SIII'
               t0 = -1128.75_dp
               t1 =  395.0_dp
               t2 = -95.0_dp
               t3 =  14000.0_dp
               x0 =  0.45_dp
               x1 =  0.0_dp
               x2 =  0.0_dp
               x3 =  1.0_dp
               w   = 120.0
               b4  = w/2.0_dp
               b4p = w/2.0_dp
               alpha = 1.0_dp
               tto = 0.0_dp
               tte = 0.0_dp
               skyrme_uses_j2_terms = .false.

            ! SGII
            ! N. V. Giai and H. Sagawa, Physics Letters B 106, 379 (1981).
            case('SGII', 'SG2')
               interaction_name = 'SGII'
               t0 = -2645.0_dp
               t1 =  340.0_dp
               t2 = -41.9_dp
               t3 =  15595.0_dp
               x0 =  0.09_dp
               x1 = -0.0588_dp
               x2 =  1.425_dp
               x3 =  0.06044_dp
               w   = 105.0_dp
               b4  = w/2.0_dp
               b4p = w/2.0_dp
               alpha = 1.0_dp/6.0_dp
               tto = 0.0_dp
               tte = 0.0_dp
               skyrme_uses_j2_terms = .false.

            ! SkM*
            ! J. Bartel et al., Nuclear Physics A 386, 79 (1982) via
            ! H. Krivine et al., Nuclear Physics A 336, 155 (1980)
            case('SKM*')
               interaction_name = 'SkM*'
               t0 = -2645.0_dp
               t1 =  410.0_dp
               t2 = -135.0_dp
               t3 =  15595.0_dp
               x0 =  0.09_dp
               x1 =  0.0_dp
               x2 =  0.0_dp
               x3 =  0.0_dp
               w   = 130.0_dp
               b4  = w/2.0_dp
               b4p = w/2.0_dp
               alpha = 1.0_dp/6.0_dp
               tto = 0.0_dp
               tte = 0.0_dp
               skyrme_uses_j2_terms = .false.

            ! SkO
            ! P.-G. Reinhard et al., Phys. Rev. C 60, 014316 (1999).
            case('SKO')
               interaction_name = 'SkO'
               t0 = -2103.653_dp
               t1 =  303.352_dp
               t2 =  791.674_dp
               t3 =  13553.252_dp
               x0 = -0.210701_dp
               x1 = -2.810752_dp
               x2 = -1.461595_dp
               x3 = -0.429881_dp
               w   =  huge(1.0d0)
               b4  =  176.578_dp
               b4p = -198.7490_dp
               alpha = 1.0_dp/4.0_dp
               tto = 0.0_dp
               tte = 0.0_dp
               skyrme_uses_j2_terms = .false.

            ! SkO'
            ! P.-G. Reinhard et al., Phys. Rev. C 60, 014316 (1999).
            case('SKOP', "SKO'")
               interaction_name = "SkO'"
               t0 = -2099.419_dp
               t1 =  301.531_dp
               t2 =  154.781_dp
               t3 =  13526.464_dp
               x0 = -0.029503_dp
               x1 = -1.325732_dp
               x2 = -2.323439_dp
               x3 = -0.147404_dp
               w   =  huge(1.0d0)
               b4  =  143.895_dp
               b4p = -82.8888_dp
               alpha = 1.0_dp/4.0_dp
               tto = 0.0_dp
               tte = 0.0_dp
               skyrme_uses_j2_terms = .true.

            ! SLy4
            ! E. Chabanat et al., Nuclear Physics A 635, 231 (1998).
            ! NOTE: HFBTHOv2 + others (Reinhard99) have slightly different
            ! values from the reference above. These values seems to come from
            ! K. Bennaceur and J. Dobaczewski, Computer Physics Communications
            ! 168, 96 (2005). I use these to be consistent with HFBTHO.
            case('SLY4')
               interaction_name = 'SLy4'
               t0 = -2488.913_dp
               t1 =  486.818_dp
               t2 = -546.395_dp
               t3 =  13777.0_dp
               x0 =  0.834_dp
               x1 = -0.344_dp
               x2 = -1.0_dp
               x3 =  1.354_dp
               w   = 123.0_dp
               b4  = w/2.0_dp
               b4p = w/2.0_dp
               alpha = 1.0_dp/6.0_dp
               tto = 0.0_dp
               tte = 0.0_dp
               skyrme_uses_j2_terms = .false.

            ! SLy4
            ! E. Chabanat et al., Nuclear Physics A 635, 231 (1998).
            ! N.b. These are the published values!
            case ('SLY4PUB')
               interaction_name = 'SLy4 (published)'
               t0 = -2488.91_dp
               t1 =  486.820_dp
               t2 = -546.390_dp
               t3 =  13777.0_dp
               x0 =  0.834_dp
               x1 = -0.344_dp
               x2 = -1.000_dp
               x3 =  1.354_dp
               w   = 123.0_dp
               b4  = w/2.0_dp
               b4p = w/2.0_dp
               alpha = 1.0_dp/6.0_dp
               tto = 0.0_dp
               tte = 0.0_dp
               skyrme_uses_j2_terms = .false.

            ! SLy5
            ! E. Chabanat et al., Nuclear Physics A 635, 231 (1998).
            ! NOTE: HFBTHOv2 + others (Reinhard99) have slightly different
            ! values from the reference above. These values seems to come from
            ! K. Bennaceur and J. Dobaczewski, Computer Physics Communications
            ! 168, 96 (2005). I use these to be consistent with HFBTHO.
            case ('SLY5')
               interaction_name = 'SLy5 (HFBTHO)'
               t0  = -2483.45_dp
               t1  =   484.23_dp
               t2  =  -556.69_dp
               t3  = 13757.00_dp
               x0  =  0.776_dp
               x1  = -0.317_dp
               x2  = -1.000_dp
               x3  =  1.263_dp
               w   =  125.0_dp
               b4  = w/2.0_dp
               b4p = w/2.0_dp
               alpha = 1.0_dp/6.0_dp
               tto = 0.0_dp
               tte = 0.0_dp
               skyrme_uses_j2_terms = .true.

            ! SLy5
            ! E. Chabanat et al., Nuclear Physics A 635, 231 (1998).
            ! N.b. These are the published values!
            case ('SLY5PUB')
               interaction_name = 'SLy5 (published)'
               t0  = -2484.88_dp
               t1  =   483.13_dp
               t2  =  -549.40_dp
               t3  = 13763.00_dp
               x0  =  0.778_dp
               x1  = -0.328_dp
               x2  = -1.000_dp
               x3  =  1.267_dp
               w   =  126.0_dp
               b4  = w/2.0_dp
               b4p = w/2.0_dp
               alpha = 1.0_dp/6.0_dp
               tto = 0.0_dp
               tte = 0.0_dp
               skyrme_uses_j2_terms = .true.

            ! T43
            ! from Lesinski et al., PRC 76 (2007) 014312 supplementary material
            case('T43')
               interaction_name = 'T43'
               t0 = -2490.275_dp
               t1 =  494.608_dp
               t2 = -255.534_dp
               t3 =  13847.12_dp
               x0 =  0.698702_dp
               x1 = -0.781655_dp
               x2 = -0.646302_dp
               x3 =  1.135795_dp
               w   = 153.103_dp
               b4  = w/2.0_dp
               b4p = w/2.0_dp
               alpha = 1.0_dp/6.0_dp
               tto = -49.160_dp
               tte = 196.868_dp
               skyrme_uses_j2_terms = .true.
            
            ! SV-min
            ! from Klüpfel et al., PRC 79 (2009) 034310
            ! and personal comm. from P.-G. Reinhard.
            ! Nb - B4 and b4p have odd definitions in the paper; check carefully
            case('SVMIN', 'SV-MIN')
               interaction_name = 'SV-min'
               t0  = -2112.248_dp 
               t1  =   295.781_dp 
               t2  =   142.268_dp 
               t3  = 13988.567_dp 
               x0  =  0.243886_dp
               x1  = -1.434926_dp
               x2  = -2.625899_dp
               x3  =  0.258070_dp
               w   =   111.291_dp
               b4  =  w/2.0_dp
               b4p =  45.93615_dp
               tto = 0
               tte = 0
               alpha = 0.255368_dp
               skyrme_uses_j2_terms = .false.
            
            ! UNEDF0
            ! From personal comm. of supplemental material for arXiv:1312.1746.
            case('UNEDF0', 'UNE0')
               interaction_name = 'UNEDF0'
               
               t0  = -1883.68781034247695000000_dp
               t1  =  277.50021223893111300000_dp
               t2  =  608.43090559153438300000_dp
               t3  =  13901.94834463432560000000_dp
               x0  =  0.00974374651612197606_dp
               x1  = -1.77784394560870229000_dp
               x2  = -1.67699034528797575000_dp
               x3  = -0.38079041463310026000_dp
               w   =  huge(1.0d0)
               b4  =  125.16100000000000000000_dp
               b4p = -91.26040000000000000000_dp
               alpha = 0.32195598958826443500_dp
               
               ! The tensor terms are used in UNEDF0 to set CJ=0. We don't need
               ! them here and they affect other terms, so remove them.
               tte = 0 !-184.46155200385393000000_dp
               tto = 0 !-118.27541333333331900000_dp
               skyrme_uses_j2_terms = .false.

            ! UNEDF1
            ! From personal comm. of supplemental material for arXiv:1312.1746.
            case('UNEDF1', 'UNE1')
               interaction_name = 'UNEDF1'
               
               t0  = -2078.3280232556408000000_dp
               t1  =  239.4008120415220450000_dp
               t2  =  1575.1195418975798900000_dp
               t3  =  14263.6462470775950000000_dp
               x0  =  0.0537569206858470316_dp
               x1  = -5.0772323818769367100_dp
               x2  = -1.3665056139429860900_dp
               x3  = -0.1624911680899563390_dp
               w   =  huge(1.0d0)
               b4  =  38.3680720616682009000_dp
               b4p =  71.3165222295833985000_dp
               alpha = 0.2700180115027076000_dp

               ! The tensor terms are used in UNEDF1 to set CJ=0. We don't need
               ! them here and they affect other terms, so remove them.
               tte = 0 !-470.3621981635254770000_dp
               tto = 0 !-203.2184521923139470000_dp
               skyrme_uses_j2_terms = .false.
               
            ! UNEDF2
            ! From personal comm. of supplemental material for arXiv:1312.1746.
            case('UNEDF2', 'UNE2')
               interaction_name = 'UNEDF2'
               
               t0  = -1735.4568519085034900000_dp
               t1  =  262.8144562353139690000_dp
               t2  =  1183.3126324330401100000_dp
               t3  =  12293.2432941125371000000_dp
               x0  =  0.1722472276241011620_dp
               x1  = -3.7669268504694102300_dp
               x2  = -1.3834981267788097900_dp
               x3  =  0.0528642727436126059_dp
               w   =  huge(1.0d0)
               b4  =  25.6586677306483040000_dp
               b4p =  77.3003893702710059000_dp
               alpha = 0.3514551325554836070_dp
               tte = -240.1405982155587540000_dp
               tto = -266.9306623187136440000_dp
               skyrme_uses_j2_terms = .true.
            
            ! UNEDF1-HFB from Schunck et al., J. Phys. G: Nucl. Part. Phys. 42 (2015) 034024
            case('UNEDF1HFB', 'UNE1HFB', 'UNEDF1-HFB', 'UNE1-HFB')   
               interaction_name = 'UNEDF1-HFB'

               t0  =  -1664.12235570454641_dp
               t1  =    255.225795498529067_dp
               t2  =   1521.27542948382984_dp
               t3  =  12064.2262021683473_dp
               x0  =      0.146274492359324770_dp
               x1  =     -4.82248318434292855_dp
               x2  =     -1.35211048250546950_dp
               x3  =      0.00431763101908833047_dp
               w   =  huge(1.0d0)
               b4  =    22.0338399999999979_dp
               b4p =   103.825096000000002_dp
               alpha =   0.377442191373274338_dp
               ! The time-odd terms are not zero, but 1) they don't affect anything
               ! (they are used to put CJ=0 which we force to be true) and 2) I'm choosing a
               ! "simple" t-odd functional and I don't want this term here.
               tte = 0 !-455.380419597060154_dp
               tto = 0 !-202.170492799999977_dp
               ! This must be .FALSE. if tte = tto = 0 above
               skyrme_uses_j2_terms = .false.
              
            case default
               write(st,'(3a)') 'ERROR: interaction "', trim(interaction_name), '" not found.'
               call writelog(st)
               call closelog
               stop
         end select
         
         ! Convert the t's and x's to coupling constants
         ! E. Perlinska et al., Phys. Rev. C 69, 014316 (2004)
         cr0 = -1.0_dp/8.0_dp  * t0 * (2.0_dp*x0 + 1.0_dp)
         crr = -1.0_dp/48.0_dp * t3 * (2.0_dp*x3 + 1.0_dp)
         cs0 = -1.0_dp/8.0_dp  * t0
         csr = -1.0_dp/48.0_dp * t3
         sigma_r = alpha;  sigma_s = alpha

         cdrho =  3.0_dp/32.0_dp *  t1*(x1 + 0.5_dp) + 1.0_dp/32.0_dp * t2*(x2 + 0.5_dp)
         ctau  = -1.0_dp/8.0_dp  *  t1*(x1 + 0.5_dp) + 1.0_dp/8.0_dp  * t2*(x2 + 0.5_dp)
         ctj0  =  1.0_dp/48.0_dp * (t1 - t2 + 10.0_dp*tte - 10.0_dp*tto)
         ctj1  =  1.0_dp/32.0_dp * (t1 - t2 - 5.0_dp*tte + 5.0_dp*tto)
         ctj2  =  1.0_dp/16.0_dp * (t1 - t2 + tte - tto)
         crdj  = -0.5_dp*b4p

         cds  =  1.0_dp/64.0_dp * (3.0_dp*t1 + t2 - 6.0_dp*tte - 2.0_dp*tto)
         ct   = -1.0_dp/16.0_dp * (t1 - t2 - 2.0_dp*tte + 2.0_dp*tto)
         cj   =  1.0_dp/8.0_dp  *  t1*(x1 + 0.5_dp) - 1.0_dp/8.0_dp*t2*(x2 + 0.5_dp)
         cgs  = -3.0_dp/32.0_dp * (3.0_dp*tte + tto)
         cf   = -3.0_dp/8.0_dp  * (tte - tto)
         csdj = -0.5_dp*b4p
         
         ! HFBTHO cannot compute individual CJ0, CJ1, and CJ2. In the presence
         ! of Skyrme's tensor interaction, our CJi are not self-consistent
         ! unless we simplify things.
         if (skyrme_uses_extended_tensor_terms .eqv. .false.) then
            ! Re-set CJ0 and CJ2 to the right multiples of CJ1
            ctj0 = 2/3.0_dp*ctj1
            ctj2 = 2*ctj1
         end if
         
         ! m*/m is a simple function of t's and x's
         ! N.b. - using hb0 = 20.735530 until that variable is passed in HFB solution
         effmass = (1+1/16.0_dp/20.73553_dp*rho_c*(3*t1+t2*(5+4*x2)))**(-1)
         return
      end if
   end subroutine load_ph_interaction


   !---------------------------------------------------------------------------
   ! Check the list of 5 coupling constant combinations from Perlinska, et al.
   !---------------------------------------------------------------------------
   function gauge_invariance_is_ok()
      implicit none
      logical :: gauge_invariance_is_ok
      real(dp) :: gauge_tolerance = 1.0d-10
      
      if ((abs(cj+ctau)        > gauge_tolerance)  .or.   &
          (abs(csdj-crdj)      > gauge_tolerance)  .or.   &
          (abs(3*ctj0+ct+2*cf) > gauge_tolerance)  .or.   &
          (abs(4*ctj1+2*ct-cf) > gauge_tolerance)  .or.   &
          (abs(2*ctj2+2*ct+cf) > gauge_tolerance)) then
         gauge_invariance_is_ok = .false.
      else
         gauge_invariance_is_ok = .true.
      end if
   end function gauge_invariance_is_ok


   !---------------------------------------------------------------------------
   ! Allocate and zero out the coupling constants.
   !---------------------------------------------------------------------------
   subroutine init_coupling_constants
      use hfbtho_basis, only: nghl

      implicit none

      if (allocated(crho)) deallocate(crho, cs, cpair, cspair)
      allocate(crho(nghl), cs(nghl), cpair(nghl), cspair(nghl))

      crho(:) = 0;  cr0 = 0;  crr = 0;  sigma_r = 0
      cs(:)   = 0;  cs0 = 0;  csr = 0;  sigma_s = 0
      cpair(:)  = 0;  cpair0  = 0;  cpairr  = 0;  sigma_pair = 0
      cspair(:) = 0;  cspair0 = 0;  cspairr = 0
      cdrho = 0;  ctau = 0;  ctj0 = 0;  ctj1 = 0;  ctj2 = 0;  crdj = 0
      cds   = 0;  ct   = 0;  cj   = 0;  csdj = 0;  cgs  = 0;  cf   = 0
   end subroutine init_coupling_constants


   !---------------------------------------------------------------------------
   ! Check against the HFBTHO solution for self-consistency.
   ! If the namelist parameter 'require_self_consistency' is .true., the
   ! program will halt if the coupling constants do not match.
   !---------------------------------------------------------------------------
   subroutine assert_self_consistency(tolerance)
      use hfbtho_basis, only: hfb_cr0, hfb_crr, hfb_cdrho, hfb_ctau, hfb_ctj,       &
                              hfb_crdj, hfb_use_j2terms, hfb_cpair, hfb_alpha_pair

      implicit none

      real(dp), intent(in) :: tolerance

      character(len=51), parameter :: fmt = "(a,2x,1f12.6,1x,a,1f12.6,1x,'[DELTA = ',1es9.2,']')"
      logical :: consistent = .true.
      character(len=160) :: st

      ! Crho[0]
      if (abs(hfb_cr0-cr0) > tolerance) then
         consistent = .false.
         write(st,fmt) 'CONSISTENCY WARNING: Crho[0] /= Crho[0](HFB).     ',   &
            cr0, 'vs.', hfb_cr0, (cr0-hfb_cr0) ; call writelog(st)
      end if

      ! Crho[r]
      if (abs(hfb_crr-crr) > tolerance) then
         consistent = .false.
         write(st,fmt) 'CONSISTENCY WARNING: Crho[r] /= Crho[r](HFB).     ',   &
            crr, 'vs.', hfb_crr, (crr-hfb_crr) ; call writelog(st)
      end if

      ! Cdrho
      if (abs(hfb_cdrho-cdrho) > tolerance) then
         consistent = .false.
         write(st,fmt) 'CONSISTENCY WARNING: Cdrho /= Cdrho(HFB).         ',   &
            cdrho, 'vs.', hfb_cdrho, (cdrho-hfb_cdrho) ; call writelog(st)
      end if

      ! Ctau
      if (abs(hfb_ctau-ctau) > tolerance) then
         consistent = .false.
         write(st,fmt) 'CONSISTENCY WARNING: Ctau /= Ctau(HFB).           ',   &
            ctau, 'vs.', hfb_ctau, (ctau-hfb_ctau) ; call writelog(st)
      end if

      ! CJ (J^2 terms)
      if ((hfb_use_j2terms .eqv. .false.)) then
         if (abs(ctj1) > tolerance) then
            consistent = .false.
            write(st,fmt) 'CONSISTENCY WARNING: 2*CJ1 /= CJ(HFB).            ',   &
               2*ctj1, 'vs.', 0.0d0, 2*ctj1 ; call writelog(st)
         end if
      else
         if (abs(hfb_ctj-2*ctj1) > tolerance) then
            consistent = .false.
            write(st,fmt) 'CONSISTENCY WARNING: 2*CJ1 /= CJ(HFB).            ',      &
               2*ctj1, 'vs.', hfb_ctj, (2*ctj1-hfb_ctj) ; call writelog(st)
         end if
         if (abs(hfb_ctj-ctj2) > tolerance) then
            consistent = .false.
            write(st,fmt) 'CONSISTENCY WARNING: CJ2 /= CJ(HFB).              ',      &
               ctj2, 'vs.', hfb_ctj, (ctj2-hfb_ctj) ; call writelog(st)   
         end if
      end if

      ! CrdJ
      if (abs(hfb_crdj-crdj) > tolerance) then
         consistent = .false.
         write(st,fmt) 'CONSISTENCY WARNING: CrdJ /= CrdJ(HFB).           ',      &
            crdj, 'vs.', hfb_crdj, (crdj-hfb_crdj) ; call writelog(st)
      end if

      ! Isovector pairing: not-density-dependent part
      if (abs(0.5_dp*sum(hfb_cpair)-8*cpair0) > tolerance) then
         consistent = .false.
         write(st,fmt) 'CONSISTENCY WARNING: 8*Crho~ /= AVG(CpV0(HFB)).   ',      &
            8*cpair0, 'vs.', 0.5_dp*sum(hfb_cpair), (8*cpair0-0.5_dp*sum(hfb_cpair))
         call writelog(st)
      end if
      
      ! Isovector pairing: surface/mixed/volume pairing
      if (abs(cpair0) > 0.0_dp) then
         if (abs(hfb_alpha_pair(1) + cpairr*rho_c**sigma_pair/cpair0) > tolerance) then
            consistent = .false.
            write(st,fmt) 'CONSISTENCY WARNING: pp_alpha /= pp_alpha(HFB).   ',   &
               -cpairr*rho_c**sigma_pair/cpair0, 'vs.', hfb_alpha_pair(1),       &
               (-cpairr*rho_c**sigma_pair/cpair0 - hfb_alpha_pair(1))
            call writelog(st)
         end if
      end if
      
      ! Isoscalar pairing: surface/mixed/volume pairing
      if (abs(cspair0) > 0.0_dp) then
         if (abs(hfb_alpha_pair(1) + cspairr*rho_c**sigma_pair/cspair0) > tolerance) then
            consistent = .false.
            write(st,fmt) 'CONSISTENCY WARNING: pp_alpha /= pp_alpha(HFB).   ',   &
               -cspairr*rho_c**sigma_pair/cspair0, 'vs.', hfb_alpha_pair(1),     &
               (-cspairr*rho_c**sigma_pair/cspair0 - hfb_alpha_pair(1))
            call writelog(st)
         end if
      end if

      if (require_self_consistency .and. (.not. consistent)) then
         write(st,'(a)') 'ERROR: this is not a self-consistent calculation.'
         call writelog(st)
         call closelog
         stop
      end if
   end subroutine assert_self_consistency


   !---------------------------------------------------------------------------
   ! Primary display of the coupling constants for STDOUT.
   !---------------------------------------------------------------------------
   subroutine print_interaction_params
      implicit none
      character(len=160) :: st
      
      ! constants and auxiliary variables for converting to natural units
      ! See PRC 82 (2010) 011304(R)
      ! Note: HFBTHO uses lambda = 700, hbarc = 197.33, causing slight mismatch in values
      real(dp), parameter :: fpi = 93  ! pion decay constant (MeV)
      real(dp), parameter :: lambda = 687  ! breakdown scale of the chiral effective theory (MeV)
      real(dp), parameter :: hbarc = 197.33  ! MeV fm
      real(dp) :: s1, s2s, s2r, s3  ! scaling factors for coupling constants with different dims
      
      ! Header
      write(st,'(a,1x,a)') 'pnFAM Interaction:', trim(interaction_name) ; call writelog(st)
      call writelog(repeat('=', 60))

      ! If no interaction, say only that
      if (skip_residual_interaction .eqv. .true.) then
         call writelog('[!] The residual interaction is OFF')

      ! If there is a residual interaction, print the coupling constants
      else

         ! If loading from file, we don't care about flags
         if (use_external_interaction .eqv. .true.) then
            call writelog('[!] Parameters loaded from external file.')
         else
            ! J^2 status
            if (force_j2_terms .eqv. .true.) then
               if (skyrme_uses_j2_terms .eqv. .true.) then
                  call writelog('[!] J^2 terms are INCLUDED (default).')
               else
                  call writelog('[!] J^2 terms are INCLUDED (normally EXCLUDED).')
               end if
            end if
            
            ! Self-consistency
            if (require_self_consistency .eqv. .true.) then
               call writelog('[!] Self-consistency is enforced')
            else
               call writelog('[!] Self-consistency is RELAXED!')
            end if

            ! Gauge invariance status
            if (require_gauge_invariance .eqv. .true.) then
               call writelog('[!] Gauge invariance is enforced')
            else
               call writelog('[!] Gauge invariance is RELAXED!')
            end if
            
            ! Custom coupling-constant overrides
            if (abs(override_cs0-huge(1d0)) > epsilon(1d0)) then
               call writelog('[!] Coupling constant Cs0 has been overridden manually')
            end if
            if (abs(override_csr-huge(1d0)) > epsilon(1d0)) then
               call writelog('[!] Coupling constant Csr has been overridden manually')
            end if
            if (abs(override_cds-huge(1d0)) > epsilon(1d0)) then
               call writelog('[!] Coupling constant Cds has been overridden manually')
            end if
            if (abs(override_ct-huge(1d0)) > epsilon(1d0)) then
               call writelog('[!] Coupling constant CT has been overridden manually')
            end if
            if (abs(override_cf-huge(1d0)) > epsilon(1d0)) then
               call writelog('[!] Coupling constant CF has been overridden manually')
            end if
            if (abs(override_cgs-huge(1d0)) > epsilon(1d0)) then
               call writelog('[!] Coupling constant Cgs has been overridden manually')
            end if
            if (abs(override_cj-huge(1d0)) > epsilon(1d0)) then
               call writelog('[!] Coupling constant Cj has been overridden manually')
            end if
            if (abs(override_csdj-huge(1d0)) > epsilon(1d0)) then
               call writelog('[!] Coupling constant Csdj has been overridden manually')
            end if
         end if

         ! Self-consistency
         call assert_self_consistency(tolerance=sc_tolerance)
         
         ! Special update for SV-min and UNEDF1-HFB
         if (trim(interaction_name) == 'SV-min' .or. trim(interaction_name) == 'UNEDF1-HFB') then
            call writelog("")
            call writelog(repeat('*',82))
            call writelog('N.b. This interaction should have all t-odd couplings zero&
               & up to gauge invariance.')
            call writelog(repeat('*',82))
         end if
         
         call writelog("")
         call writelog('Coupling Constants')
         call writelog(repeat('-',60))
         write(st,'("Crho[0] = ",1f15.9,";   Crho[r] = ",1f15.9,";   sigma = ",1f12.9)')   &
            cr0, crr, sigma_r ; call writelog(st)
         write(st,'("Cs[0]   = ",1f15.9,";   Cs[r]   = ",1f15.9,";   sigma = ",1f12.9)')   &
            cs0, csr, sigma_s ; call writelog(st)

         !write(lg,*)
         write(st,'("Cdr     = ",1f15.9,";   Cds     = ",1f15.9)') cdrho, cds ; call writelog(st)
         write(st,'("Ctau    = ",1f15.9,";   CT      = ",1f15.9)') ctau, ct   ; call writelog(st)
         write(st,'("CtJ0    = ",1f15.9,";   Cj      = ",1f15.9)') ctj0, cj   ; call writelog(st)

         !write(lg,*)
         write(st,'("CtJ1    = ",1f15.9,";   Cgs     = ",1f15.9)') ctj1, cgs  ; call writelog(st)
         write(st,'("CtJ2    = ",1f15.9,";   CF      = ",1f15.9)') ctj2, cf   ; call writelog(st)
         write(st,'("CrdJ    = ",1f15.9,";   Csdj    = ",1f15.9)') crdj, csdj ; call writelog(st)

         call writelog("")
         write(st,'("Cpr[0]  = ",1f15.9,";   Cpr[r]  = ",1f15.9,";   sigma = ",1f12.9)')   &
            cpair0, cpairr, sigma_pair ; call writelog(st)
         write(st,'("Cps[0]  = ",1f15.9,";   Cps[r]  = ",1f15.9)') cspair0, cspairr ; call writelog(st)
         
         ! Scaling factors to translate coupling constants to natural units
         s1 = fpi**2 * hbarc**(-3)
         s2r = fpi**(2*(sigma_r + 1)) * lambda**sigma_r * hbarc**(-3*(1+sigma_r))
         s2s = fpi**(2*(sigma_s + 1)) * lambda**sigma_s * hbarc**(-3*(1+sigma_s))
         s3 = fpi**2 * lambda**2 * hbarc**(-5)
         
         call writelog("")
         write(st,'(A,F0.1,A)') 'Coupling Constants in natural units (using Lambda = ', lambda, ' MeV)'
         call writelog(st)
         call writelog(repeat('-',60))
         write(st,'("Crho[0] = ",1f15.9,";   Crho[r] = ",1f15.9,";   sigma = ",1f12.9)')   &
            cr0*s1, crr*s2r, sigma_r ; call writelog(st)
         write(st,'("Cs[0]   = ",1f15.9,";   Cs[r]   = ",1f15.9,";   sigma = ",1f12.9)')   &
            cs0*s1, csr*s2s, sigma_s ; call writelog(st)

         !write(lg,*)
         write(st,'("Cdr     = ",1f15.9,";   Cds     = ",1f15.9)') cdrho*s3, cds*s3 ; call writelog(st)
         write(st,'("Ctau    = ",1f15.9,";   CT      = ",1f15.9)') ctau*s3, ct*s3   ; call writelog(st)
         write(st,'("CtJ0    = ",1f15.9,";   Cj      = ",1f15.9)') ctj0*s3, cj*s3   ; call writelog(st)

         !write(lg,*)
         write(st,'("CtJ1    = ",1f15.9,";   Cgs     = ",1f15.9)') ctj1*s3, cgs*s3  ; call writelog(st)
         write(st,'("CtJ2    = ",1f15.9,";   CF      = ",1f15.9)') ctj2*s3, cf*s3   ; call writelog(st)
         write(st,'("CrdJ    = ",1f15.9,";   Csdj    = ",1f15.9)') crdj*s3, csdj*s3 ; call writelog(st)

      end if
      
      ! Landau parameters of the calculation
      call writelog("")
      call print_landau_params
   end subroutine print_interaction_params
   
   !----------------------------------------------------------------------------
   ! Compute approximate Landau parameters in the spin-isospin channel 
   !----------------------------------------------------------------------------
   subroutine print_landau_params
      use constants, only : pi
      implicit none
      character(len=160) :: st
      ! Unused, so commented out to avoid warnings -TS
      ! real(dp), parameter :: beta = (3*pi**2/2.0_dp)**(2/3.0_dp)
      real(dp) :: kf, n0, g0p, g1p, h0p, aux
      
      call writelog("Approximate Landau Parameters")
      call writelog(repeat('-',60))
      
      if (skip_residual_interaction .eqv. .true.) then
         call writelog('* No Landau parameters -- residual interaction is off.')
      else if (interaction_name(1:5) == 'FILE:') then
         call writelog('* No Landau parameters -- interaction loaded from file.')
      else
         ! Fermi energy from saturation density (not computed! taken as 0.16 fm**-3)
         kf = (3*pi**2/2.0_dp*rho_c)**(1/3.0_dp)
      
         ! Normalization using hb0 = 20.73553 MeV fm**2,
         n0 = (20.73553_dp*pi**2/effmass/kf)**(-1)
      
         ! Parameters
         g0p =    n0*(2*(cs0+csr*rho_c**sigma_s) + 2*ct*kf*kf + 2*cf*kf*kf/3)
         g1p = -2*n0*ct*kf*kf - 2*cf*kf*kf/3
         h0p = n0*2*kf*kf*cf/3
      
         write(st,'("rho_00  = ",1f15.9,";   k_F     = ",1f15.9)') rho_c, kf ; call writelog(st)
         write(st,'("m*/m    = ",1f15.9,";   N0^(-1) = ",1f15.9)') effmass, 1/n0 ; call writelog(st)
         write(st,'("g0''     = ",1f15.9,";   g1''     = ",1f15.9)') g0p, g1p ; call writelog(st)
         write(st,'("h0''     = ",1f15.9)') h0p ; call writelog(st)
         
         ! Warn if the Landau parameters imply unstable infinite nuclear matter
         ! Stability conditions derived in NPA 312 (1979) 10,
         ! also listed in PRC 81 (2010) 044302
         aux = 1 + g1p/3 - 10*h0p/3
         if (aux <= 0) then
            call writelog("[!] Stability condition 1 + 1/3*g1' - 10/3*h0' > 0 not fulfilled!")
            write(st,'(1f15.9," <= 0")') aux; call writelog(st)
         end if
         aux = 1 + g1p/3 + 5*h0p/3
         if (aux <= 0) then
            call writelog("[!] Stability condition 1 + 1/3*g1' + 5/3*h0' > 0 not fulfilled!")
            write(st,'(1f15.9," <= 0")') aux; call writelog(st)
         end if
         aux = 1 + g1p/3 - h0p/3
         if (aux <= 0) then
            call writelog("[!] Stability condition 1 + 1/3*g1' - 1/3*h0' > 0 not fulfilled!")
            write(st,'(1f15.9," <= 0")') aux; call writelog(st)
         end if
         aux = 1 + g0p/2 + sqrt(g0p*g0p + 8*h0p*h0p)/2
         if (aux <= 0) then
            call writelog("[!] Stability condition 1 + 1/2*g0' + sqrt(g0'^2 + 8*h0'^2) > 0 not fulfilled!")
            write(st,'(1f15.9," <= 0")') aux; call writelog(st)
         end if
         aux = 1 + g0p/2 - sqrt(g0p*g0p + 8*h0p*h0p)/2
         if (aux <= 0) then
            call writelog("[!] Stability condition 1 + 1/2*g0' - sqrt(g0'^2 + 8*h0'^2) > 0 not fulfilled!")
            write(st,'(1f15.9," <= 0")') aux; call writelog(st)
         end if
         
      end if
   end subroutine print_landau_params

end module interaction
