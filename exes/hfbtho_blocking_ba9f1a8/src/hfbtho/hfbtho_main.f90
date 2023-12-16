!***********************************************************************
!
!    Copyright (c) 2020, Lawrence Livermore National Security, LLC.
!                        Produced at the Lawrence Livermore National
!                        Laboratory.
!                        Lead developer: Nicolas Schunck, schunck1@llnl.gov
!    HFODD
!    -----
!      LLNL-CODE-710577 All rights reserved.
!      LLNL-CODE-470611 All rights reserved.
!
!      Copyright 2017, N. Schunck, J. Dobaczewski, W. Satula, P. Baczyk,
!                      J. Dudek, Y. Gao, M. Konieczka, K. Sato, Y. Shi,
!                      X.B. Wang, and T.R. Werner
!      Copyright 2012, N. Schunck, J. Dobaczewski, J. McDonnell,
!                      W. Satula, J.A. Sheikh, A. Staszczak,
!                      M. Stoitsov, P. Toivanen
!      Copyright 2009, J. Dobaczewski, W. Satula, B.G. Carlsson, J. Engel,
!                      P. Olbratowski, P. Powalowski, M. Sadziak,
!                      J. Sarich, N. Schunck, A. Staszczak, M. Stoitsov,
!                      M. Zalewski, H. Zdunczuk
!      Copyright 2004, 2005, J. Dobaczewski, P. Olbratowski
!      Copyright 1997, 2000, J. Dobaczewski, J. Dudek
!
!    HFBTHO
!    -----
!      LLNL-CODE-728299 All rights reserved.
!      LLNL-CODE-573953 All rights reserved.
!
!      Copyright 2017, R. Navarro Perez, N. Schunck, R. Lasseri, C. Zhang,
!                      J. Sarich
!      Copyright 2012, M.V. Stoitsov, N. Schunck, M. Kortelainen, H.A. Nam,
!                      N. Michel, J. Sarich, S. Wild
!      Copyright 2005, M.V. Stoitsov, J. Dobaczewski, W. Nazarewicz, P.Ring
!
!
!    This file is part of DFTNESS.
!
!    DFTNESS is free software: you can redistribute it and/or modify it
!    under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    DFTNESS is distributed in the hope that it will be useful, but
!    WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with DFTNESS. If not, see <http://www.gnu.org/licenses/>.
!
!    OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC
!    LICENSE
!
!    Our Preamble Notice
!
!      A. This notice is required to be provided under our contract
!         with the U.S. Department of Energy (DOE). This work was
!         produced at the Lawrence Livermore National Laboratory under
!         Contract No. DE-AC52-07NA27344 with the DOE.
!      B. Neither the United States Government nor Lawrence Livermore
!         National Security, LLC nor any of their employees, makes any
!         warranty, express or implied, or assumes any liability or
!         responsibility for the accuracy, completeness, or usefulness
!         of any information, apparatus, product, or process disclosed,
!         or represents that its use would not infringe privately-owned
!         rights.
!      C. Also, reference herein to any specific commercial products,
!         process, or services by trade name, trademark, manufacturer
!         or otherwise does not necessarily constitute or imply its
!         endorsement, recommendation, or favoring by the United States
!         Government or Lawrence Livermore National Security, LLC. The
!         views and opinions of authors expressed herein do not
!         necessarily state or reflect those of the United States
!         Government or Lawrence Livermore National Security, LLC, and
!         shall not be used for advertising or product endorsement
!         purposes.
!
!    The precise terms and conditions for copying, distribution and
!    modification are contained in the file COPYING.
!
!***********************************************************************

Program hfbthoprog
  Use HFBTHO
  Use HFBTHO_utilities
  Character(1024) :: filename_hfbtho
  Character(1024) :: filename_unedf

  filename_hfbtho = 'hfbtho_NAMELIST.dat'
  filename_unedf  = 'UNEDF_NAMELIST.dat'

  Call Main_Program(filename_hfbtho, filename_unedf)
End Program hfbthoprog
  !=======================================================================
  !
  !=======================================================================
  Subroutine Main_Program(filename_hfbtho, filename_unedf)
    Use HFBTHO_utilities
    Use HFBTHO
#if(USE_QRPA==1)
    Use HFBTHO_storage
#endif
#if(USE_MPI==2)
    Use HFBTHO_mpi_communication
#endif
#if(DO_MASSTABLE==1 || DRIP_LINES==1 || DO_PES==1)
    Use HFBTHO_large_scale
#endif
#if(READ_FUNCTIONAL==1)
    Use HFBTHO_read_functional
#endif
    Use HFBTHO_localization
    Use HFBTHO_solver
    Use HFBTHO_projections
    Implicit None
#if(USE_MPI==1)
    Include 'mpif.h'
#endif
    Character(*), Intent(IN) :: filename_hfbtho
    Character(*), Intent(IN) :: filename_unedf

    Integer(ipr) :: iblocase(2),nkblocase(2,5)
    Integer(ipr) :: i,it,icount,jcount,l,noForce,icalc
    Integer(ipr) :: j, lambda
    ! Initialize MPI environment and possibly create subcommunicators
#if(USE_MPI==2)
    call MPI_INIT(ierr_mpi)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, ierr_mpi)
    call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_taskid, ierr_mpi)
#if(DRIP_LINES==1)
    call Create_MPI_Teams(number_deformations)
#endif
#if(DO_PES==1)
    Call Create_MPI_Teams(HFB_cores)
#endif
#else
    mpi_size = 1
    mpi_taskid = 0
#endif
    !
    Write(ID_string,'(i6.6)') mpi_taskid
    task_error = 0
    !-------------------------------------------------------------------
    ! Read input and namelist data for the requested nucleus
    !-------------------------------------------------------------------
    Call initialize_HFBTHO_NAMELIST
    !-------------------------------------------------------------------
    !  Process 0 reads general input data. In serial mode, that is all
    !  the data needed.
    !-------------------------------------------------------------------
    if(mpi_taskid == 0) then
       Call read_HFBTHO_NAMELIST(filename_hfbtho)
#if(DO_MASSTABLE==1)
       call read_HFBTHO_MassTable
#endif
#if(DO_PES==1)
       Call read_HFBTHO_PES
#endif
#if(DRIP_LINES==1)
       call read_HFBTHO_StableLine
#endif
#if(READ_FUNCTIONAL==1)
       call read_HFBTHO_Functional
#endif
    End If
    !-------------------------------------------------------------------
    ! Broadcast of process-independent input data in parallel mode
    !-------------------------------------------------------------------
#if(USE_MPI==2)
    if(mpi_taskid == 0) then
       call Construct_Vectors
       allocate(task_error_gthr(0:mpi_size-1))
    else
       allocate(task_error_gthr(0:0))
    End If
    call broadcast_vectors
    if(nrows > 0) then
       if(mpi_taskid > 0) call Deconstruct_Vectors
    End If
#if(READ_FUNCTIONAL==1)
    call broadcast_functional
#endif
#endif
    ! Overwrite basis characteristics if so requested
#if(USE_MPI==0)
    If(lambda_active(2) > 0 .And. automatic_basis) &
       Call adjust_basis(expectation_values(2),.False.,proton_number,neutron_number)
#endif
    !-------------------------------------------------------------------
    ! Allocation of vectors for mass tables
    !-------------------------------------------------------------------
#if(DO_MASSTABLE==1)
    if(nrows > 0) then
#if(USE_MPI==2)
       call allocate_out_vectors
#else
       call allocate_mass_table
#endif
    End If
#endif
    !-------------------------------------------------------------------
    ! Mass table calculation
    !-------------------------------------------------------------------
#if(DO_MASSTABLE==1)
    pi = 4*atan(1.0_pr)
    icalc = 0
    !team leader creates a file for bookkeeping
    if(mpi_taskid == 0) then
       open(127,file='TableLog.dat')
    End If
    !loop over elements of the mass table
    do iRow = 0,nRows
       Z_chain = Z_masstable(iRow)
       N_chain = N_masstable(iRow)
       A_chain = Z_chain+N_chain
       beta_deformation = beta_masstable(iRow)
       Q20 = Q20_masstable(iRow)!beta_deformation*sqrt(5/pi)*(A_chain)**(5/3._pr)/100._pr
       Write(row_string,'("_",i6.6)') irow
       if(irow == 0.and.nrows > 0) cycle
       if(mpi_taskid == 0) then
          Write(127,'(a7,2i5,2f15.8)') row_string,Z_chain,N_chain,Q20,beta_deformation
       End If
       !only do the calculations that correspond to your task id
       if(mod(irow,mpi_size) /= mpi_taskid) cycle
       proton_number  = Z_chain
       neutron_number = N_chain
       expectation_values(2) = Q20
       basis_deformation = beta_deformation
#endif
    !-------------------------------------------------------------------
    ! Potential energy surface calculation
    !-------------------------------------------------------------------
#if(DO_PES==1)
    icalc = 0
    ! Team leader creates a file for bookkeeping
    If(mpi_taskid == 0) Then
       Open(127,file='TableLog.dat')
    End If
    ! Loop over points in the potential energy surface
    Do iRow = 1,npoints
       Z_chain = Z_PES(iRow)
       N_chain = N_PES(iRow)
       A_chain = Z_chain+N_chain
       Write(row_string,'("_",i6.6)') iRow
       If(mpi_taskid == 0) Then
          Write(127,'(a7,2i5,9f10.3)') row_string,Z_chain,N_chain,(Q_PES(iRow,j),j=1,ndefs)
       End If
       ! Only do the calculation that correspond to your team
       If(Mod(iRow,number_teams) /= team_color) Cycle
       proton_number  = Z_chain
       neutron_number = N_chain
       ! Default basis deformation and WS deformation based on input file
       If(bet2_PES(iRow) > -8.0) Then
          basis_deformation = bet2_PES(iRow)
          beta2_deformation = bet2_PES(iRow)
       End If
       If(bet3_PES(iRow) > -8.0) beta3_deformation = bet3_PES(iRow)
       If(bet4_PES(iRow) > -8.0) beta4_deformation = bet4_PES(iRow)
       Do j=1,ndefs
          lambda = lambda_PES(j)
          ! More advanced fit based on value of Q2 only
          If(lambda == 2 .And. automatic_basis) Call adjust_basis(Q_PES(iRow,j),.False.,proton_number,neutron_number)
          expectation_values(lambda) = Q_PES(iRow,j)
          lambda_active(lambda) = 1
       End Do
#endif
    !-------------------------------------------------------------------
    ! Dripline calculation
    !-------------------------------------------------------------------
#if(DRIP_LINES==1)
    beta_step = 1.0_pr/real(number_deformations-1,kind=pr)
    pi = 4*atan(1.0_pr)
    Write(team_string,'(1i3.3)') team_color
    !team leader allocates array to recieve energies and creates a file for bookkeeping
    if(team_rank == 0) then
       allocate(energy_chain_gthr(0:team_size-1))
       open(127,file='TeamTable'//team_string//'.dat')
    End If
    calc_counter = 0
    !loop over the nuclei in the "stable line"
    do iRow = 1,nRows
       !only do the chains that correspond to your team
       if(mod(irow,number_teams) /= team_color) cycle
       Minimum_Energy_Prev = 100
       Z_chain = Z_stable_line(iRow)
       N_chain = N_stable_line(iRow)
       !move along the isotopic (or isotonic) chain until the drip line is reached
       do
          A_chain = Z_chain + N_chain
          beta_deformation = -0.5_pr - beta_step
          !loop over the different basis deformations of each nucleus
          do i_deformation = 1,number_deformations
             beta_deformation = beta_deformation + beta_step
             calc_counter = calc_counter + 1
             Q20 = beta_deformation*sqrt(5.0_pr/pi)*(A_chain)**(5.0_pr/3.0_pr)/100._pr
             Write(row_string,'("_",a3,"_",i6.6)') team_string,calc_counter
             !team leader writes type of calculation for bookkeeping
             if(team_rank == 0) then
                if(direction_sl(irow) == 1) then
                   direction_str = ' isotopic'
                else
                   direction_str = ' isotonic'
                End If
                Write(127,'(a11,2i5,2f15.8,a9)') row_string,z_chain,n_chain,Q20,beta_deformation,direction_str
             End If
             if(number_deformations == 1) beta_deformation = 0._pr
             !only calculate what corresponds to each process
             if(mod(i_deformation,team_size) /= team_rank) cycle
             proton_number  = Z_chain
             neutron_number = N_chain
             expectation_values(2) =  Q20
             basis_deformation = beta_deformation
#endif
       !-------------------------------------------------------------------
       ! Regular HFBTHO execution, whether or not mass table or dripline
       ! mode is activated
       !-------------------------------------------------------------------
       !memo: Namelist /HFBTHO_GENERAL/ number_of_shells,oscillator_length, basis_deformation, &
       !                                proton_number,neutron_number,type_of_calculation
       !      Namelist /HFBTHO_INITIAL/ beta2_deformation, beta3_deformation, beta4_deformation
       !      Namelist /HFBTHO_ITERATIONS/ number_iterations, accuracy, restart_file
       !      Namelist /HFBTHO_FUNCTIONAL/ functional, add_initial_pairing, type_of_coulomb, include_3N_force
       !      Namelist /HFBTHO_PAIRING/ user_pairing, vpair_n, vpair_p, pairing_cutoff, pairing_feature
       !      Namelist /HFBTHO_CONSTRAINTS/ lambda_values, lambda_active, expectation_values
       !      Namelist /HFBTHO_BLOCKING/ proton_blocking(1:5), neutron_blocking(1:5)
       !      Namelist /HFBTHO_PROJECTION/ switch_to_THO,projection_is_on,gauge_points,delta_Z,delta_N
       !      Namelist /HFBTHO_TEMPERATURE/ set_temperature, temperature
       !      Namelist /HFBTHO_FEATURES/ collective_inertia, fission_fragments, pairing_regularization, &
       !                                 automatic_basis, localization_functions
       !      Namelist /HFBTHO_TDDFT/ filter, fragment_properties, real_Z, real_N
       !      Namelist /HFBTHO_NECK/ set_neck_constrain, neck_value
       !      Namelist /HFBTHO_DEBUG/ number_Gauss, number_Laguerre, number_Legendre, &
       !                              compatibility_HFODD, number_states, force_parity, print_time
       !
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
       !
       ! Checking consistency of *_INI variables
       Call check_consistency
       If(ierror_flag /= 0) Then
          Write(lout,'(a33,i4)') ' ERRORS IN Main_Program, process',mpi_taskid
          Write(lout,'(2(a3,i4),a18,f6.1)')' Z=',npr_INI(2),' N=',npr_INI(1), ' basis deformation', q_INI
          Do i=1,ierror_flag
             Write(lout,'(a8,i4,a12,i2,2x,a)') ' process',mpi_taskid,', error_flag=',i,ierror_info(i)
          End Do
          Write(lout,'("Terminating very early...")')
#if((DO_PES==1 || DO_MASSTABLE==1 || DRIP_LINES==1))
          Cycle
#else
          Return
#endif
       End If
       ! Check if there is at least one constraint
       icount=0
       Do l=1,lambdaMax
          If(lambda_active(l) > 0) icount=icount+1
       End Do
       ! If there is at least one constraint, check if any breaks parity
       If(icount > 0) Then
          jcount=0
          Do l=1,lambdaMax,2
             If(lambda_active(l) > 0) jcount=jcount+1
          End Do
          If(jcount > 0) Parity_INI=.False.
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
       If(set_temperature .And. Abs(temper) <= 1.e-10_pr) switch_on_temperature=.False.
       !
       Call read_UNEDF_NAMELIST(skyrme_INI,noForce)
       If(skyrme_INI == 'SeaLL1' .Or. skyrme_INI == 'SLY4mod') Then
          pairing_regularization = .True.; pwi_INI=100.0
       End If
       ! If functional is used, projection automatically switched off
       If(noForce == 0) iproj_INI=0
       call set_all_gaussians(icou_INI)
       ! Read parameters of the energy functionals from a file
#if(READ_FUNCTIONAL==1)
       call replace_functional
#endif
       !---------------------------------------------------------------------------
       !                                BLOCKING
       !  Blocking candidates are selected from HFB results in the even-even parent nucleus
       !  (with Z+/- 1 or/and N +/- 1 particles).
       !   - If nkblocase(it,2) = 0, we block all blocking candidates within an
       !     energy window of pwiblo = 1 MeV (set in preparer()) around the Fermi level of
       !     the even-even parent nucleus
       !   - If nkblocase(it,2)  /=  0, we block the single candidate identified by the quantum
       !     numbers contained in the array nkblocase
       !---------------------------------------------------------------------------
       iblocase=0; bloqpdif=zero ! blomax will be charged from the previous solution
       Do it=1,2
          ! Blocking all qp within an energy window
          If(nkblocase(it,1) /= 0 .And. nkblocase(it,2) == 0) Then
             If(it == 1) Then
                iblocase(1)=iblocase(1)+1
                If(iblocase(1) > blomax(1)) iblocase(1)=1
             Else
                If(iblocase(1) <= 1) iblocase(2)=iblocase(2)+1
             End If
             nkblo_INI(it,1)=Sign(iblocase(it),nkblocase(it,1))
             nkblo_INI(it,2)=0
          Else
             ! Blocking a single qp
             nkblo_INI(it,:)=nkblocase(it,:)
          End If
       End Do
       !--------------------------------------------------------------------
       ! Run the solver in all cases EVEN/ODDS, FITS/NO-FITS
       !--------------------------------------------------------------------
       Call HFBTHO_DFT_SOLVER
       !--------------------------------------------------------------------
       ! Restore broken symmetries of HFB states
       !--------------------------------------------------------------------
       If(PNP_is_on >= 1 .Or. AMP_is_on >= 1) Call HFBTHO_restore()
       !--------------------------------------------------------------------
       ! Calculate localization functions
       !--------------------------------------------------------------------
       If(localization_functions) Then
          Call localization()
       End If
       !--------------------------------------------------------------------
       ! Display error messages in case of problems
       !--------------------------------------------------------------------
#if(USE_MPI==1 || USE_MPI==2)
       If(do_print == 1) Then
          If(ierror_flag == 0) Then
             Write(lout,*)
             Write(lout,'(a)') ' HFBTHO_SOLVER ended without errors'
             Write(lout,*)
          Else
             task_error = 1
             Write(lout,*)
             Write(lout,'(a33,i4)') ' ERRORS IN HFBTHO_SOLVER, process',mpi_taskid
             Do i=1,ierror_flag
                Write(lout,'(a8,i4,a12,i2,2x,a)') ' process',mpi_taskid,', error_flag=',i,ierror_info(i)
             End Do
             Write(lout,*)
          End If
          ierror_flag = 0
       End If
#else
       If(ierror_flag /= 0) Then
          task_error = 1
          Write(lout,*)
          Write(lout,'(a)') ' ERRORS IN HFBTHO_SOLVER'
          Do i=1,ierror_flag
             Write(lout,'(a,i2,2x,a)') ' error_flag=',i,ierror_info(i)
          End Do
          Write(lout,*)
          ierror_flag = 0
       End If
#endif
    !-------------------------------------------------------------------
    ! Mass table calculation
    !-------------------------------------------------------------------
#if(DO_MASSTABLE==1)
       if(nrows > 0) then
          Write(*,'("task ",a6," finished row ",a6)') ID_string, row_string(2:7)
#if(USE_MPI==2)
          call fill_out_vectors(icalc)
#else
          call fill_mass_table(icalc)
#endif
          icalc = icalc + 1
       End If
       Close(lfile) ! close the output
    End Do
    if(nrows > 0) then
#if(USE_MPI==2)
       call gather_results
#endif
       call print_mass_table
    End If
    if(team_rank == 0) close(127)
#endif
    !-------------------------------------------------------------------
    ! Potential energy surface calculation
    !-------------------------------------------------------------------
#if(DO_PES==1)
       If(npoints > 0) Then
          Write(*,'("task ",a6," finished row ",a6)') ID_string, row_string(2:7)
          icalc = icalc + 1
       End If
       If(do_print == 1) Close(lfile) ! close the output
    End Do ! End of loop over points in the PES
    If(do_print == 1) Close(127)
#endif
    !-------------------------------------------------------------------
    ! Dripline calculation
    !-------------------------------------------------------------------
#if(DRIP_LINES==1)
             !Calculations without Lipkin-Nogami
             if(kindhfb_INI > 0) then
                Energy_chain = ehfb
             !Calculations with Lipkin-Nogami
             else
                Energy_chain = etot
             End If
             Close(lfile) ! close the output
          ! end loop over deformations
          End Do
          call find_minimum_energy
          separation_2N = Minimum_Energy_Prev - Minimum_Energy
          Minimum_Energy_Prev = Minimum_Energy
          if(separation_2N < 0._pr) then
             ! Drip line has been reached
             exit
          else
             ! Drip line has not been reached
             if(direction_sl(irow) == 1) then
                if(N_chain >= 310) exit
                N_chain = N_chain + 2
             else
                if(Z_chain >= 120) exit
                Z_chain = Z_chain + 2
             End If
          End If
       ! end isotopic (or isotonic) chain
       End Do
       !team leader announces isotopic (or isotonic) chain finished
       if(team_rank == 0) then
          if(direction_sl(irow) == 1) then
             Write(*,'(a4,i3,a27,i4,a6,i4)') 'team',team_color,' finished isotopic chain Z=',Z_chain,' at N=',N_chain
          else
             Write(*,'(a4,i3,a27,i4,a6,i4)') 'team',team_color,' finished isotonic chain N=',N_chain,' at Z=',Z_chain
          End If
       End If
       Close(lfile) ! close the output
    !end loop over nuclei inside the "stable line
    End Do
    !team leader closes bookkeeping file
    if(mpi_taskid == 0) close(127)
#endif

#if(USE_MPI==2)
    call mpi_gather(task_error,1,mpi_integer,task_error_gthr,1,mpi_integer,0,MPI_COMM_WORLD,ierr_mpi)
    if(mpi_taskid == 0) then
       if(sum(task_error_gthr) == 0) then
          Write(*,*)
          Write(*,'(a)') ' Parallel execution ended without errors'
          Write(*,*)
       End If
    End If
#endif
#if(DO_PES==0 && DO_MASSTABLE==0 && DRIP_LINES==0)
    If(lout < lfile) Close(lfile) ! close the output
    If(do_print == 1) Close(lout) ! close the output
#endif

#if(USE_MPI==1 || USE_MPI==2)
    ! Wait here until all processes are done
    Call mpi_barrier(MPI_COMM_WORLD, ierr_mpi)
    ! Close MPI environment
    Call mpi_finalize(ierr_mpi)
#endif
    !
  End Subroutine Main_Program
