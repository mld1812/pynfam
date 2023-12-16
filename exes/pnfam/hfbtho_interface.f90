!------------------------------------------------------------------------------
!> An interface for the hfbtho library to initialize the solution and basis.
!> Requires the hfbtho_NAMELIST.dat and corresponding hfbtho_output.wel file.
!>
!> This module is a pared-down version of hfbtho_main that runs hfbtho with
!> zero iterations, no constraints, etc. This recreates all basis quantities, reads
!> in the fields from the wel file, and uses them to construct and diagonlize the
!> hfb matrix to get the solution (U's, V's, E's) needed by pnfam. With my updates
!> to hfbtho_solver (6/23/2020) this should result in the exact solution used to
!> populate thoout.dat.
!>
!> @authors E.M. Ney, UNC Chapel Hill, 2018-
!------------------------------------------------------------------------------
Module hfbtho_interface

  Implicit None

Contains

  Subroutine HFBTHO_program(iexit)
    Use HFBTHO_utilities
    Use HFBTHO
#if(READ_FUNCTIONAL==1)
    Use HFBTHO_read_functional
#endif
    Use HFBTHO_solver
    Use HFBTHO_storage
    Use HFBTHO_io, only : Z_r
  
    Implicit None
    Integer(ipr), intent(out) :: iexit

    Character(1024), parameter :: filename_hfbtho='hfbtho_NAMELIST.dat'
    Character(1024), parameter :: filename_unedf='UNEDF_NAMELIST.dat'
  
    Integer(ipr) :: iblocase(2),nkblocase(2,5)
    Integer(ipr) :: i,it,icount,jcount,l,noForce,icalc
    Integer(ipr) :: j, lambda
  
    !-------------------------------------------------------------------
    ! Read input and namelist data for the requested nucleus
    !-------------------------------------------------------------------
    Call initialize_HFBTHO_NAMELIST
    Call read_HFBTHO_NAMELIST(filename_hfbtho)
#if(READ_FUNCTIONAL==1)
    Call read_HFBTHO_Functional
#endif
  
    !-------------------------------------------------------------------
    ! Override inputs for PNFAM to just reconstruct solution
    !-------------------------------------------------------------------
    lout=6; lfile=6
    write_hel = .false. ! don't write result
    !do_print = 0        ! don't print
    number_of_shells       = -abs(number_of_shells) ! don't print
    lambda_active          = (/ 0, 0, 0, 0, 0, 0, 0, 0 /)
    expectation_values     = (/ 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr /)
    number_iterations      = 0
    restart_file           = -1
    localization_functions = .false.
    collective_inertia     = .false.
    fission_fragments      = .false.
  
    n00_INI                   = number_of_shells       ! number of shells
    b0_INI                    = oscillator_length      ! oscillator length
    q_INI                     = basis_deformation      ! deformation beta_2 of the basis
    npr_INI(1)                = neutron_number         ! N
    npr_INI(2)                = proton_number          ! Z
    kindhfb_INI               = type_of_calculation    ! 1: HFB, -1: HFB+LN
    !
    b2_0                      = beta2_deformation      ! beta2 parameter of the initial WS solution
    b3_0                      = beta3_deformation      ! beta3 parameter of the initial WS solution
    b4_0                      = beta4_deformation      ! beta4 parameter of the initial WS solution
    !
    MAX_ITER_INI              = number_iterations      ! max number of iterations
    epsi_INI                  = accuracy               ! convergence of iterations
    inin_INI                  = restart_file           ! restart from file
    !
    skyrme_INI                = TRIM(functional)       ! functional
    Add_Pairing_INI           = add_initial_pairing    ! add pairing starting from file
    icou_INI                  = type_of_coulomb        ! coul: no-(0), dir.only-(1), plus exchange-(2)
    use_3N_couplings          = include_3N_force       ! Include 3N force on certain DME functionals
    !
    set_pairing               = user_pairing           ! pairing is defined by user if .True.
    V0n_INI                   = vpair_n                ! pairing strength for neutrons
    V0p_INI                   = vpair_p                ! pairing strength for protons
    pwi_INI                   = pairing_cutoff         ! pairing q.p. cutoff
    cpv1_INI                  = pairing_feature        ! Type of pairing: volume, surface, mixed
    !
    nkblocase(1,:)            = neutron_blocking       ! config. of neutron blocked state
    nkblocase(2,:)            = proton_blocking        ! config. of proton blocked state
    !
    iLST_INI                  = switch_to_THO          ! 0:HO, -1:HO->THO, 1:THO
    keypjn_INI                = gauge_points           ! PNP: number of gauge points
    keypjp_INI                = gauge_points           ! PNP: number of gauge points
    iproj_INI                 = projection_is_on       ! projecting on different nucleus
    npr1pj_INI                = delta_N                ! its neutron number
    npr2pj_INI                = delta_Z                ! its proton number
    !
    switch_on_temperature     = set_temperature        ! switches on temperature mode
    temper                    = temperature            ! value of the temperature
    !
    collective_inertia        = collective_inertia     ! calculate collective mass and zero-point energy
    fission_fragments         = fission_fragments      ! calculate fission fragment characteristics
    pairing_regularization    = pairing_regularization ! activates the regularization of the pairing force
    automatic_basis           = automatic_basis        ! computes localization functions
    localization_functions    = localization_functions ! computes localization functions
    !
    neck_constraints          = set_neck_constrain     ! activate the constraint on the neck
    neckRequested             = neck_value             ! set the requested value for the neck
    !
    ngh_INI                   = number_Gauss           ! number of Gauss-Hermite points for z-direction
    ngl_INI                   = number_Laguerre        ! number of Gauss-Laguerre points for rho-direction
    nleg_INI                  = number_Legendre        ! number of Gauss-Legendre points for Coulomb
    basis_HFODD_INI           = compatibility_HFODD    ! flag to enforce same basis as HFODD
    nstate_INI                = number_states          ! total number of states in basis
    Parity_INI                = force_parity           ! reflection symmetry
    IDEBUG_INI                = print_time             ! debug
    DO_FITT_INI               = .False.                ! calculates quantities for reg.optimization
    Print_HFBTHO_Namelist_INI = .False.                ! Print Namelist
  
    !-------------------------------------------------------------------
    ! Perform some checks and initialziations
    !-------------------------------------------------------------------
    ! Checking consistency of *_INI variables
    Call check_consistency
    If(ierror_flag/=0) Then
       Write(lout,'(1x,a22)') 'ERRORS IN Main_Program'
       Do i=1,ierror_flag
          Write(lout,'(1x,a11,i2,2x,a)') 'error_flag=',i,ierror_info(i)
       End Do
       Write(lout,'(1x,"Terminating very early...")')
       Return
    End If
    ! Check if there is at least one constraint
    icount=0
    Do l=1,lambdaMax
       If(lambda_active(l)>0) icount=icount+1
    End Do
    ! If there is at least one constraint, check if any breaks parity
    If(icount>0) Then
       jcount=0
       Do l=1,lambdaMax,2
          If(lambda_active(l)>0) jcount=jcount+1
       End Do
       If(jcount>0) Parity_INI=.False.
    Else
       collective_inertia = .False.
    End If
    If(fission_fragments) Parity_INI=.False.
    ! For fission fragment properties: nearest integer (just for printing)
    If(fragment_properties) Then
       tz_fragments(1) = real_N; tz_fragments(2) = real_Z
       npr_INI(1) = Nint(tz_fragments(1))
       npr_INI(2) = Nint(tz_fragments(2))
    End If
    ! Enforces no-temperature mode if T<1.e-10
    If(set_temperature .And. Abs(temper)<=1.e-10_pr) switch_on_temperature=.False.
  
    !-------------------------------------------------------------------
    ! Setup functional
    !-------------------------------------------------------------------
    Call read_UNEDF_NAMELIST(skyrme_INI,noForce)
    If(skyrme_INI=='SeaLL1' .Or. skyrme_INI=='SLY4mod') Then
       pairing_regularization = .True.; pwi_INI=100.0
    End If
    ! If functional is used, projection automatically switched off
    If(noForce==0) iproj_INI=0
    Call set_all_gaussians(icou_INI)
    ! Read parameters of the energy functionals from a file
#if(READ_FUNCTIONAL==1)
    Call replace_functional
#endif
  
    !-------------------------------------------------------------------
    ! Setup blocking
    !-------------------------------------------------------------------
    iblocase=0; bloqpdif=zero
    Do it=1,2
       ! Blocking all qp within an energy window
       If(nkblocase(it,1)/=0.And.nkblocase(it,2)==0) Then
          If(it==1) Then
             iblocase(1)=iblocase(1)+1
             If(iblocase(1)>blomax(1)) iblocase(1)=1
          Else
             If(iblocase(1)<=1) iblocase(2)=iblocase(2)+1
          End If
          nkblo_INI(it,1)=Sign(iblocase(it),nkblocase(it,1))
          nkblo_INI(it,2)=0
       Else
          ! Blocking a single qp
          nkblo_INI(it,:)=nkblocase(it,:)
       Endif
    End Do
  
    !--------------------------------------------------------------------
    ! Run the solver to populate variables
    !--------------------------------------------------------------------
    iexit=0; Z_r = -1
    Call HFBTHO_DFT_SOLVER
    ! Hack to check if we successfully read from wel...
    ! (otherwise we just have a wood saxon solution)
    If(Z_r == -1) iexit = 1
  
    !--------------------------------------------------------------------
    ! Display error messages in case of problems
    !--------------------------------------------------------------------
    If(lout<lfile) Then
       If (ierror_flag/=0) Then
          Write(lout,*)
          Write(lout,'(a)') ' ERRORS IN HFBTHO_SOLVER'
          Do i=1,ierror_flag
             Write(lout,'(a,i2,2x,a)') ' error_flag=',i,ierror_info(i)
          End Do
          Write(lout,*)
       Else
          Write(lout,*)
          Write(lout,'(a)') ' HFBTHO_SOLVER ended without errors'
          Write(lout,*)
       End If
  
       Close(lfile)! close thoout.dat
       Close(lout) ! close output
    End If
  
  End Subroutine HFBTHO_program

End Module hfbtho_interface
