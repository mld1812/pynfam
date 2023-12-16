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

! ==================================================================== !
!                                                                      !
!                  LARGE-SCALE CALCULATION PACKAGE                     !
!                                                                      !
! ==================================================================== !
!----------------------------------------------------------------------
!>  This module contains various utility functions to perform large
!>  scale calculations with HFBTHO, be it full nuclear chart from
!>  dripline to dripline or a subset of nuclei with different
!>  deformations. Calculations of this type require USE_MPI=2 (the
!>  value USE_MPI=1 is reserved for optimization with the ANL code).
!>  Refer to hfbtho_main.f90 for further information.
!>
!>  @author
!>  Rodrigo Navarro Perez
!----------------------------------------------------------------------
!  Subroutines: - read_HFBTHO_MassTable
!               - allocate_mass_table
!               - fill_mass_table
!               - print_mass_table
!----------------------------------------------------------------------
Module HFBTHO_large_scale
  Use HFBTHO_utilities
  Use HFBTHO
  Implicit None
Contains

#if(DO_MASSTABLE==1)
  !=======================================================================
  !>  Opens and reads the 'hfbtho_MASSTABLE.dat' file that indicates
  !>  the nuclei (with different deformations) to be calculated with
  !>  HFBTHO
  !>    - The first integer read from the file gives the number of
  !>      calculations to be made (called nrows along the source code)
  !>    - Each row contains the number of protons (Z), number of neutrons
  !>      (N), requested value for the quadrupole moment (Q20) and basis
  !>      deformation (beta)
  !>    - If nrows is 0 (or negative) the rest of the file will be
  !>      ignored and a regular HFBTHO calculation will be made with the
  !>      parameters given in the 'hfbtho_NAMELIST.dat' file.
  !>    - If nrows is greater than zero the first nrows calculations will
  !>      be made with the parameters given in each row. Other input
  !>      parameters will be take from the 'hfbtho_NAMELIST.dat' file.
  !>  nrows should not be larger than the actual number of rows on the
  !>  file.
  !=======================================================================
  subroutine read_HFBTHO_MassTable
    implicit none
    integer(ipr) :: lmasstable=15,i
    open(lmasstable,file='hfbtho_MASSTABLE.dat')
    read(lmasstable,*)
    read(lmasstable,*) nRows
    read(lmasstable,*)
    allocate(   Z_masstable(0:nRows))
    allocate(   N_masstable(0:nRows))
    allocate( Q20_masstable(0:nRows))
    allocate(beta_masstable(0:nRows))
       Z_masstable(0) =  proton_number
       N_masstable(0) = neutron_number
     Q20_masstable(0) = expectation_values(2)
    beta_masstable(0) = basis_deformation
    do i=1,nRows
       read(lmasstable,*) Z_masstable(i),N_masstable(i),Q20_masstable(i)&
            ,beta_masstable(i)
    enddo
    close(lmasstable)
  end subroutine read_HFBTHO_MassTable
  !=======================================================================
  !>  Allocates the arrays necessary to construct the output mass table
  !=======================================================================
  subroutine allocate_mass_table
    implicit none
    allocate(ierrors_out(1:nrows))
    allocate(Z_out(1:nrows))
    allocate(N_out(1:nrows))
    allocate(Q20_in(1:nrows))
    allocate(beta_in(1:nrows))
    allocate(E_HFB_out(1:nrows))
    allocate(Q20Z_out(1:nrows))
    allocate(Q20N_out(1:nrows))
    allocate(Q20T_out(1:nrows))
  end subroutine allocate_mass_table
  !=======================================================================
  !>  The results from the HFBTHO calculation are stored in arrays.
  !>  This subroutine is used in the case that the mass table calculated
  !>  without MPI parallelism
  !=======================================================================
  subroutine fill_mass_table(icalc)
    implicit none
    integer, intent(in) :: icalc !< Calculations performed so far
    ierrors_out(icalc+1) = ierror_flag
    Z_out(icalc+1) = Z_masstable(iRow)
    N_out(icalc+1) = N_masstable(iRow)
    Q20_in(icalc+1) = Q20_masstable(iRow)
    beta_in(icalc+1) = beta_masstable(iRow)
    if(kindhfb_INI.gt.0) then
       E_HFB_out(icalc+1) = ehfb
    else
       E_HFB_out(icalc+1) = etot
    endif
    Q20Z_out(icalc+1) = qmoment(2,1)
    Q20N_out(icalc+1) = qmoment(2,2)
    Q20T_out(icalc+1) = qmoment(2,3)
  end subroutine fill_mass_table
  !=======================================================================
  !>  The output mass table is written to the 'MassTableOut.dat' file
  !>  and on screen
  !=======================================================================
  subroutine print_mass_table
    implicit none
    integer :: i
    if(mpi_taskid.eq.0) then
       open(100,file='MassTableOut.dat')
       do i = 1,nrows
          write(  *,'(4i5,6f15.8)') i,ierrors_out(i),Z_out(i),N_out(i),Q20_in(i),beta_in(i),E_HFB_out(i),Q20Z_out(i),Q20N_out(i),Q20T_out(i)
          write(100,'(4i5,6f15.8)') i,ierrors_out(i),Z_out(i),N_out(i),Q20_in(i),beta_in(i),E_HFB_out(i),Q20Z_out(i),Q20N_out(i),Q20T_out(i)
       enddo
       close(100)
    endif
  end subroutine print_mass_table
#endif
#if(DRIP_LINES==1)
  !=======================================================================
  !>  This subroutine reads the location of the valley of stability as
  !>  well as the direction (proton or neutron) in which dripline
  !>  calculations should proceed.
  !=======================================================================
  subroutine read_HFBTHO_StableLine
    implicit none
    integer(ipr) :: ldripline=15,i
    open(ldripline,file='hfbtho_STABLELINE.dat')
    read(ldripline,*)
    read(ldripline,*) nRows
    read(ldripline,*)
    allocate(Z_stable_line(0:nRows))
    allocate(N_stable_line(0:nRows))
    allocate(direction_sl(0:nRows))
    Z_stable_line(0) =  proton_number
    N_stable_line(0) = neutron_number
    do i=1,nRows
       read(ldripline,*) Z_stable_line(i),N_stable_line(i),direction_sl(i)
    enddo
    close(ldripline)
  end subroutine read_HFBTHO_StableLine
#endif
#if(DO_PES==1)
  !=======================================================================
  !>  This subroutine defines a potential energy surface (PES) by reading
  !>  the file hfbtho_PES.dat. This routine is only relevant when the
  !>  DO_PES flag is activated
  !=======================================================================
  Subroutine read_HFBTHO_PES
    Implicit None
    Integer(ipr) :: unit_PES=15,i,j,lambda
    !
    Open(unit_PES,file='hfbtho_PES.dat')
    Read(unit_PES,*)
    ! Read the total number of points in the PES and the number of constraints
    Read(unit_PES,*) npoints,ndefs
    nRows = npoints; n_real_masstable_in = 3+ndefs
    Read(unit_PES,*)
    ! Read the multipolarity lambda of each constraint
    Allocate(lambda_PES(1:ndefs))
    Read(unit_PES,*) (lambda_PES(i),i=1,ndefs)
    Read(unit_PES,*)
    ! Allocate arrays for Z, N, collective variables and initial deformations. The default
    ! value at index 0 corresponds to the initialization of the Namelist
    Allocate(Z_PES(0:npoints))
    Allocate(N_PES(0:npoints))
    Allocate(Q_PES(0:npoints,1:ndefs))
    Allocate(bet2_PES(0:npoints))
    Allocate(bet3_PES(0:npoints))
    Allocate(bet4_PES(0:npoints))
    Z_PES(0) = proton_number
    N_PES(0) = neutron_number
    Do j=1,ndefs
       lambda = lambda_PES(j)
       Q_PES(0,j) = expectation_values(lambda)
    End Do
    bet2_PES(0) = beta2_deformation
    bet3_PES(0) = beta3_deformation
    bet4_PES(0) = beta4_deformation
    Do i=1,npoints
       Read(unit_PES,*) Z_PES(i),N_PES(i),(Q_PES(i,j),j=1,ndefs),bet2_PES(i),bet3_PES(i),bet4_PES(i)
    End Do
    !
    Close(unit_PES)
  End Subroutine read_HFBTHO_PES
#endif

end module HFBTHO_large_scale


#if(USE_MPI==2)
!----------------------------------------------------------------------
!>  This module contain utility functions and subroutines to
!>  exchange information among MPI processes involved in a mass table
!>  or dripline calculation.
!>
!> @author
!> Rodrigo Navarro Perez
!----------------------------------------------------------------------
!  Subroutines: - Create_MPI_teams
!               - Construct_Vectors
!               - broadcast_vectors
!               - allocate_mpi_vectors
!               - allocate_out_vectors
!               - Deconstruct_Vectors
!               - gather_results
!               - fill_out_vectors(icalc)
!----------------------------------------------------------------------
Module HFBTHO_mpi_communication
  Use HFBTHO_utilities
  Use HFBTHO
  Implicit None
  Include 'mpif.h'
  Integer, PRIVATE, SAVE :: debug = 1
Contains

  !=======================================================================
  !>  Creates new communicators for teams of MPI processes to calculate
  !>  each nucleus with different deformations. The processes are
  !>  distributed as evenly as possible and each team will have at
  !>  most as many members as the number of deformations that will be used
  !>  (having more members would be inefficient)
  !>
  !>  A group and communicator are created with one member of each team
  !>  that will act as team leader and communicate with a 'super-leader'
  !>  to distribute the nucleus to be calculated by each group
  !=======================================================================
  subroutine Create_MPI_Teams(n_cores)
    use mpi
    implicit none
    integer, intent(in) :: n_cores
    !Determine number of teams (each team must have at most n_cores members)
#if(DO_PES==1)
    number_teams = mpi_size/n_cores
#else
    number_teams = (mpi_size-1)/n_cores+1
#endif
    team_color = mod(mpi_taskid,number_teams)
    !Create Team communicators by splitting the world communicator
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,team_color,mpi_taskid,COMM_team,ierr_mpi)
    !Get size team and rank within the team
    call MPI_COMM_SIZE(COMM_team, team_size, ierr_mpi)
    call MPI_COMM_RANK(COMM_team, team_rank, ierr_mpi)
    If(debug > 0) Then
       If(mpi_taskid == 0) Then
          Write(6,'("mpi_size = ",i6," number of teams",i6)') mpi_size,number_teams
       End If
       Write(6,'("Process",i6," has rank ",i6," in team of color",i6)') mpi_taskid,team_rank,team_color
    End If
  end subroutine Create_MPI_Teams
  !=======================================================================
  !>  After the master process has read the 'hfbtho_NAMELIST.dat' and
  !>  'hfbtho_MASSTABLE.dat' files, it constructs (allocates and fills out)
  !>  a set of vectors to be broadcast to the other processes.
  !>  - vector_int_mpi contains all the integers to broadcasted
  !>  - vector_real_mpi contains all the reals to broadcasted
  !>  - vector_logicals_mpi contains all the logicals to broadcasted
  !>
  !>  Vector_sizes contains 4 integers, the size of the three previous
  !>  arrays and the number of rows on the mass table.
  !=======================================================================
  subroutine Construct_Vectors
    implicit none
    integer :: i,j,ntot,ncons
#if(DO_PES==1)
    allocate(vector_sizes(1:5))
    ncons=ndefs
#else
    allocate(vector_sizes(1:4))
    ncons=0
#endif
    vector_sizes(1) = 34 + 2*lambdamax +  n_int_masstable_in*nRows + ncons! integers
    vector_sizes(2) = 12 +   lambdamax + n_real_masstable_in*nRows ! reals
    vector_sizes(3) = 12 ! boolean
    vector_sizes(4) =  nRows
#if(DO_PES==1)
    vector_sizes(5) =  ndefs
#endif
    allocate( vector_int_mpi(1:vector_sizes(1)))
    allocate(vector_real_mpi(1:vector_sizes(2)))
    allocate( vector_log_mpi(1:vector_sizes(3)))
    ! Inputs of type 'integer'
    ntot = 1
    vector_int_mpi(ntot) = number_of_shells
    ntot = ntot + 1
    vector_int_mpi(ntot) = proton_number
    ntot = ntot + 1
    vector_int_mpi(ntot) = neutron_number
    ntot = ntot + 1
    vector_int_mpi(ntot) = type_of_calculation
    ntot = ntot + 1
    vector_int_mpi(ntot) = number_iterations
    ntot = ntot + 1
    vector_int_mpi(ntot) = type_of_coulomb
    ntot = ntot + 1
    vector_int_mpi(ntot) = restart_file
    ntot = ntot + 1
    vector_int_mpi(ntot) = projection_is_on
    ntot = ntot + 1
    vector_int_mpi(ntot) = gauge_points
    ntot = ntot + 1
    vector_int_mpi(ntot) = delta_Z
    ntot = ntot + 1
    vector_int_mpi(ntot) = delta_N
    ntot = ntot + 1
    vector_int_mpi(ntot) = switch_to_THO
    ntot = ntot + 1
    vector_int_mpi(ntot) = number_Gauss
    ntot = ntot + 1
    vector_int_mpi(ntot) = number_Laguerre
    ntot = ntot + 1
    vector_int_mpi(ntot) = number_Legendre
    ntot = ntot + 1
    vector_int_mpi(ntot) = number_states
    ntot = ntot + 1
    vector_int_mpi(ntot) = print_time
    do i = 1,5
       vector_int_mpi(ntot+i)   =  proton_blocking(i)
       vector_int_mpi(ntot+5+i) = neutron_blocking(i)
    enddo
    ntot = ntot + 2*5 ! should be 27 now
    do i = 1,lambdamax
       vector_int_mpi(ntot+i) = lambda_values(i)
       vector_int_mpi(ntot+lambdamax+i) = lambda_active(i)
    enddo
    ntot = ntot + 2*lambdamax ! should be 45 by now
    ! Projection
    ntot = ntot + 1
    vector_int_mpi(ntot) = PNP_is_on
    ntot = ntot + 1
    vector_int_mpi(ntot) = number_of_gauge_points
    ntot = ntot + 1
    vector_int_mpi(ntot) = delta_neutrons
    ntot = ntot + 1
    vector_int_mpi(ntot) = delta_protons
    ntot = ntot + 1
    vector_int_mpi(ntot) = AMP_is_on
    ntot = ntot + 1
    vector_int_mpi(ntot) = number_of_rotational_angles
    ntot = ntot + 1
    vector_int_mpi(ntot) = maximal_angular_momentum
#if(DO_MASSTABLE==1)
    do i = 1,nRows
       vector_int_mpi(ntot + n_int_masstable_in*(i-1)+1) = Z_masstable(i)
       vector_int_mpi(ntot + n_int_masstable_in*(i-1)+2) = N_masstable(i)
    enddo
#endif
#if(DO_PES==1)
    Do i = 1,npoints
       vector_int_mpi(ntot + n_int_masstable_in*(i-1)+1) = Z_PES(i)
       vector_int_mpi(ntot + n_int_masstable_in*(i-1)+2) = N_PES(i)
    End Do
    Do j = 1, ndefs
       vector_int_mpi(ntot + n_int_masstable_in*npoints + j) = lambda_PES(j)
    End Do
#endif
#if(DRIP_LINES==1)
    do i = 1,nRows
       vector_int_mpi(ntot + n_int_masstable_in*(i-1)+1) = Z_stable_line(i)
       vector_int_mpi(ntot + n_int_masstable_in*(i-1)+2) = N_stable_line(i)
       vector_int_mpi(ntot + n_int_masstable_in*(i-1)+3) = direction_sl(i)
    enddo
#endif
    ! Inputs of type 'real'
    ntot = 1
    vector_real_mpi(ntot) = oscillator_length
    ntot = ntot + 1
    vector_real_mpi(ntot) = basis_deformation
    ntot = ntot + 1
    vector_real_mpi(ntot) = beta2_deformation
    ntot = ntot + 1
    vector_real_mpi(ntot) = beta3_deformation
    ntot = ntot + 1
    vector_real_mpi(ntot) = beta4_deformation
    ntot = ntot + 1
    vector_real_mpi(ntot) = accuracy
    ntot = ntot + 1
    vector_real_mpi(ntot) = temperature
    ntot = ntot + 1
    vector_real_mpi(ntot) = vpair_n
    ntot = ntot + 1
    vector_real_mpi(ntot) = vpair_p
    ntot = ntot + 1
    vector_real_mpi(ntot) = pairing_cutoff
    ntot = ntot + 1
    vector_real_mpi(ntot) = pairing_feature
    ntot = ntot + 1
    vector_real_mpi(ntot) = neck_value
    do i = 1,lambdamax
       vector_real_mpi(ntot+i) = expectation_values(i)
    enddo
    ntot = ntot + lambdamax
#if(DO_MASSTABLE==1)
    do i = 1,nRows
       vector_real_mpi(ntot+n_real_masstable_in*(i-1)+1) = Q20_masstable(i)
       vector_real_mpi(ntot+n_real_masstable_in*(i-1)+2) = beta_masstable(i)
    enddo
#endif
#if(DO_PES==1)
    Do j=1,ndefs
       Do i=1,npoints
          vector_real_mpi(ntot + (j-1)*npoints + i) = Q_PES(i,j)
       End Do
    End Do
    Do i=1,npoints
       vector_real_mpi(ntot + ndefs   *npoints + i) = bet2_PES(i)
       vector_real_mpi(ntot +(ndefs+1)*npoints + i) = bet3_PES(i)
       vector_real_mpi(ntot +(ndefs+2)*npoints + i) = bet4_PES(i)
    End Do
#endif
    ! Inputs of type 'logical'
    ntot = 1
    vector_log_mpi(ntot) = add_initial_pairing
    ntot = ntot + 1
    vector_log_mpi(ntot) = set_temperature
    ntot = ntot + 1
    vector_log_mpi(ntot) = collective_inertia
    ntot = ntot + 1
    vector_log_mpi(ntot) = fission_fragments
    ntot = ntot + 1
    vector_log_mpi(ntot) = pairing_regularization
    ntot = ntot + 1
    vector_log_mpi(ntot) = localization_functions
    ntot = ntot + 1
    vector_log_mpi(ntot) = set_neck_constrain
    ntot = ntot + 1
    vector_log_mpi(ntot) = compatibility_HFODD
    ntot = ntot + 1
    vector_log_mpi(ntot) = force_parity
    ntot = ntot + 1
    vector_log_mpi(ntot) = user_pairing
    ntot = ntot + 1
    vector_log_mpi(ntot) = automatic_basis
    ntot = ntot + 1
    vector_log_mpi(ntot) = include_3N_force
  end subroutine Construct_Vectors
  !=======================================================================
  !>  The master process broadcasts first the vector_sizes arrays, the
  !>  other processes allocate the other arrays (int, reals and logicals)
  !>  then the master process broadcasts integers, reals, logicals and
  !>  a single string
  !=======================================================================
  subroutine broadcast_vectors
    use mpi
    implicit none
#if(DO_PES==1)
    if(mpi_taskid > 0) allocate(vector_sizes(1:5))
    call mpi_bcast(vector_sizes,5,mpi_integer,0,mpi_comm_world,ierr_mpi)
    if(mpi_taskid > 0) then
       nRows = vector_sizes(4); npoints=nRows; ndefs = vector_sizes(5)
       call allocate_mpi_vectors
    endif
#else
    if(mpi_taskid > 0) allocate(vector_sizes(1:4))
    call mpi_bcast(vector_sizes,4,mpi_integer,0,mpi_comm_world,ierr_mpi)
    if(mpi_taskid > 0) then
       nRows = vector_sizes(4)
       call allocate_mpi_vectors
    endif
#endif
    call mpi_bcast( vector_int_mpi,vector_sizes(1),mpi_integer,0,mpi_comm_world,ierr_mpi)
    call mpi_bcast(vector_real_mpi,vector_sizes(2),mpi_double_precision,0,mpi_comm_world,ierr_mpi)
    call mpi_bcast( vector_log_mpi,vector_sizes(3),mpi_logical,0,mpi_comm_world,ierr_mpi)
    call mpi_bcast(functional,30,mpi_character,0,mpi_comm_world,ierr_mpi)
  end subroutine broadcast_vectors
  !=======================================================================
  !>  Non-master processes allocate the arrays that will be received
  !>  from the master's broadcast.
  !=======================================================================
  subroutine allocate_mpi_vectors
    implicit none
    allocate( vector_int_mpi(1:vector_sizes(1)))
    allocate(vector_real_mpi(1:vector_sizes(2)))
    allocate( vector_log_mpi(1:vector_sizes(3)))
#if(DO_MASSTABLE==1)
    allocate(   Z_masstable(0:nRows))
    allocate(   N_masstable(0:nRows))
    allocate( Q20_masstable(0:nRows))
    allocate(beta_masstable(0:nRows))
#endif
#if(DO_PES==1)
    Allocate(Z_PES(0:npoints))
    Allocate(N_PES(0:npoints))
    Allocate(lambda_PES(1:ndefs))
    Allocate(Q_PES(0:npoints,1:ndefs))
    Allocate(bet2_PES(0:npoints))
    Allocate(bet3_PES(0:npoints))
    Allocate(bet4_PES(0:npoints))
#endif
#if(DRIP_LINES==1)
    allocate(Z_stable_line(0:nRows))
    allocate(N_stable_line(0:nRows))
    allocate(direction_sl(0:nRows))
#endif
  end subroutine allocate_mpi_vectors
  !=======================================================================
  !>  After receiving the master's broadcast, the non-master processes
  !>  'deconstruct' the arrays by setting the parameters to be used
  !>  in the HFBTHO calculations
  !=======================================================================
  subroutine Deconstruct_Vectors
    implicit none
    integer :: i,j,ntot
    ! Inputs of type 'integer'
    ntot = 1
    number_of_shells    = vector_int_mpi(ntot)
    ntot = ntot + 1
    proton_number       = vector_int_mpi(ntot)
    ntot = ntot + 1
    neutron_number      = vector_int_mpi(ntot)
    ntot = ntot + 1
    type_of_calculation = vector_int_mpi(ntot)
    ntot = ntot + 1
    number_iterations   = vector_int_mpi(ntot)
    ntot = ntot + 1
    type_of_coulomb     = vector_int_mpi(ntot)
    ntot = ntot + 1
    restart_file        = vector_int_mpi(ntot)
    ntot = ntot + 1
    projection_is_on    = vector_int_mpi(ntot)
    ntot = ntot + 1
    gauge_points        = vector_int_mpi(ntot)
    ntot = ntot + 1
    delta_Z             = vector_int_mpi(ntot)
    ntot = ntot + 1
    delta_N             = vector_int_mpi(ntot)
    ntot = ntot + 1
    switch_to_THO       = vector_int_mpi(ntot)
    ntot = ntot + 1
    number_Gauss        = vector_int_mpi(ntot)
    ntot = ntot + 1
    number_Laguerre     = vector_int_mpi(ntot)
    ntot = ntot + 1
    number_Legendre     = vector_int_mpi(ntot)
    ntot = ntot + 1
    number_states       = vector_int_mpi(ntot)
    ntot = ntot + 1
    print_time          = vector_int_mpi(ntot)
    do i = 1,5
       proton_blocking(i)  = vector_int_mpi(ntot+i)
       neutron_blocking(i) = vector_int_mpi(ntot+5+i)
    enddo
    ntot = ntot + 2*5 ! should be 27 now
    do i = 1,lambdamax
       lambda_values(i) = vector_int_mpi(ntot+i)
       lambda_active(i) = vector_int_mpi(ntot+lambdamax+i)
    enddo
    ntot = ntot + 2*lambdamax ! should be 27 now
    ! Projection
    ntot = ntot + 1
    PNP_is_on                   = vector_int_mpi(ntot)
    ntot = ntot + 1
    number_of_gauge_points      = vector_int_mpi(ntot)
    ntot = ntot + 1
    delta_neutrons              = vector_int_mpi(ntot)
    ntot = ntot + 1
    delta_protons               = vector_int_mpi(ntot)
    ntot = ntot + 1
    AMP_is_on                   = vector_int_mpi(ntot)
    ntot = ntot + 1
    number_of_rotational_angles = vector_int_mpi(ntot)
    ntot = ntot + 1
    maximal_angular_momentum    = vector_int_mpi(ntot)
#if(DO_MASSTABLE==1)
    do i = 1,nRows
       Z_masstable(i) = vector_int_mpi(ntot + n_int_masstable_in*(i-1)+1)
       N_masstable(i) = vector_int_mpi(ntot + n_int_masstable_in*(i-1)+2)
    enddo
#endif
#if(DO_PES==1)
    Do i = 1,npoints
       Z_PES(i) = vector_int_mpi(ntot + n_int_masstable_in*(i-1)+1)
       N_PES(i) = vector_int_mpi(ntot + n_int_masstable_in*(i-1)+2)
    End Do
    Do j = 1, ndefs
       lambda_PES(j) = vector_int_mpi(ntot + n_int_masstable_in*npoints + j)
    End Do
#endif
#if(DRIP_LINES==1)
    do i = 1,nRows
       Z_stable_line(i) = vector_int_mpi(ntot + n_int_masstable_in*(i-1)+1)
       N_stable_line(i) = vector_int_mpi(ntot + n_int_masstable_in*(i-1)+2)
       direction_sl(i)  = vector_int_mpi(ntot + n_int_masstable_in*(i-1)+3)
    enddo
#endif
    ! Inputs of type 'real'
    ntot = 1
    oscillator_length = vector_real_mpi(ntot)
    ntot = ntot + 1
    basis_deformation = vector_real_mpi(ntot)
    ntot = ntot + 1
    beta2_deformation = vector_real_mpi(ntot)
    ntot = ntot + 1
    beta3_deformation = vector_real_mpi(ntot)
    ntot = ntot + 1
    beta4_deformation = vector_real_mpi(ntot)
    ntot = ntot + 1
    accuracy          = vector_real_mpi(ntot)
    ntot = ntot + 1
    temperature       = vector_real_mpi(ntot)
    ntot = ntot + 1
    vpair_n           = vector_real_mpi(ntot)
    ntot = ntot + 1
    vpair_p           = vector_real_mpi(ntot)
    ntot = ntot + 1
    pairing_cutoff    = vector_real_mpi(ntot)
    ntot = ntot + 1
    pairing_feature   = vector_real_mpi(ntot)
    ntot = ntot + 1
    neck_value        = vector_real_mpi(ntot)
    do i = 1,lambdamax
       expectation_values(i) =  vector_real_mpi(ntot+i)
    enddo
    ntot = ntot + lambdamax
#if(DO_MASSTABLE==1)
    do i = 1,nRows
        Q20_masstable(i) = vector_real_mpi(ntot+n_real_masstable_in*(i-1)+1)
       beta_masstable(i) = vector_real_mpi(ntot+n_real_masstable_in*(i-1)+2)
    enddo
#endif
#if(DO_PES==1)
    Do j=1,ndefs
       Do i=1,npoints
          Q_PES(i,j) = vector_real_mpi(ntot + (j-1)*npoints + i)
       End Do
    End Do
    Do i=1,npoints
        bet2_PES(i) = vector_real_mpi(ntot + ndefs   *npoints + i)
        bet3_PES(i) = vector_real_mpi(ntot +(ndefs+1)*npoints + i)
        bet4_PES(i) = vector_real_mpi(ntot +(ndefs+2)*npoints + i)
    End Do
#endif
    ! Inputs of type 'logical'
    ntot = 1
    add_initial_pairing    = vector_log_mpi(ntot)
    ntot = ntot + 1
    set_temperature        = vector_log_mpi(ntot)
    ntot = ntot + 1
    collective_inertia     = vector_log_mpi(ntot)
    ntot = ntot + 1
    fission_fragments      = vector_log_mpi(ntot)
    ntot = ntot + 1
    pairing_regularization = vector_log_mpi(ntot)
    ntot = ntot + 1
    localization_functions = vector_log_mpi(ntot)
    ntot = ntot + 1
    set_neck_constrain     = vector_log_mpi(ntot)
    ntot = ntot + 1
    compatibility_HFODD    = vector_log_mpi(ntot)
    ntot = ntot + 1
    force_parity           = vector_log_mpi(ntot)
    ntot = ntot + 1
    user_pairing           = vector_log_mpi(ntot)
    ntot = ntot + 1
    automatic_basis        = vector_log_mpi(ntot)
    ntot = ntot + 1
    include_3N_force       = vector_log_mpi(ntot)
  end subroutine Deconstruct_Vectors
#if(DRIP_LINES==1)
  !=======================================================================
  !>  Get and broadcast the value of the minimum energy for a given nucleus
  !=======================================================================
  subroutine find_minimum_energy()
    use mpi
    implicit none
    !team leader gathers energies and looks for the minimum
    call mpi_gather(Energy_Chain,1,mpi_double_precision,Energy_Chain_gthr,1,mpi_double_precision,0,COMM_team,ierr_mpi)
    if(team_rank.eq.0) then
       Minimum_Energy = minval(Energy_Chain_gthr,1)
    endif
    !team leader broadcasts the minimum
    call mpi_bcast(Minimum_Energy,1,mpi_double_precision,0,COMM_team,ierr_mpi)
  end subroutine find_minimum_energy
#endif
#if(DO_MASSTABLE==1)
  !=======================================================================
  !>  Every process allocates the vectors that will be gathered by the master.
  !=======================================================================
  subroutine allocate_out_vectors
    implicit none
    deallocate(vector_sizes,vector_int_mpi,vector_real_mpi,vector_log_mpi)
    nrows_task = nrows/mpi_size
    if(mpi_taskid.gt.0.and.mpi_taskid.le.mod(nrows,mpi_size)) then
       nrows_task = nrows_task + 1
    endif
    allocate(vector_sizes(1:2))
    vector_sizes(1) = n_int_masstable_out*nrows_task
    vector_sizes(2) = n_real_masstable_out*nrows_task
    allocate( vector_int_mpi(1:vector_sizes(1)))
    allocate(vector_real_mpi(1:vector_sizes(2)))
  end subroutine allocate_out_vectors
  !=======================================================================
  !>  After each HFBTHO calculation, each process saves the results in
  !>  vectors that will be gathered later by the master.
  !=======================================================================
  subroutine fill_out_vectors(icalc)
    implicit none
    integer, intent(in) :: icalc
       vector_int_mpi(n_int_masstable_out*icalc+1) = ierror_flag
       vector_int_mpi(n_int_masstable_out*icalc+2) = Z_masstable(iRow)
       vector_int_mpi(n_int_masstable_out*icalc+3) = N_masstable(iRow)
       vector_real_mpi(n_real_masstable_out*icalc+1) = Q20_masstable(iRow)
       vector_real_mpi(n_real_masstable_out*icalc+2) = beta_masstable(iRow)
       if(kindhfb_INI.gt.0) then
          vector_real_mpi(n_real_masstable_out*icalc+3) = ehfb
       else
          vector_real_mpi(n_real_masstable_out*icalc+3) = etot
       endif
       vector_real_mpi(n_real_masstable_out*icalc+4) = qmoment(2,1)
       vector_real_mpi(n_real_masstable_out*icalc+5) = qmoment(2,2)
       vector_real_mpi(n_real_masstable_out*icalc+6) = qmoment(2,3)
  end subroutine fill_out_vectors
  !=======================================================================
  !>  The master process gathers the result from all calculations
  !>  and fills out the resulting mass table including number of protons,
  !>  number of neutrons, input quadrupole moment, basis deformation, HFB
  !>  energy, output quadrupole moment from protons, neutrons and total
  !=======================================================================
  subroutine gather_results
    use mpi
    implicit none
    integer :: i,j
    if(mpi_taskid.eq.0) then
       allocate(vector_sizes_gthr(1:2*mpi_size))
       allocate(vector_sizes_int_gthr(0:mpi_size-1))
       allocate(vector_sizes_real_gthr(0:mpi_size-1))
       allocate(vector_disp_int_gthr(0:mpi_size-1))
       allocate(vector_disp_real_gthr(0:mpi_size-1))
    endif
    call mpi_gather(vector_sizes,2,mpi_integer,vector_sizes_gthr,2,&
         mpi_integer,0,mpi_comm_world,ierr_mpi)

    if(mpi_taskid.eq.0) then
       disp_int = 0
       disp_real = 0
       do i = 0,mpi_size-1
          vector_sizes_int_gthr(i)  = vector_sizes_gthr(2*i+1)
          vector_sizes_real_gthr(i) = vector_sizes_gthr(2*i+2)
          if(i.eq.0) cycle
          vector_disp_int_gthr(i)  = disp_int
          vector_disp_real_gthr(i) = disp_real
          disp_int = disp_int + vector_sizes_gthr(2*i+1)
          disp_real = disp_real + vector_sizes_gthr(2*i+2)
       enddo
       vector_disp_int_gthr(0)  = disp_int
       vector_disp_real_gthr(0) = disp_real
       allocate( vector_int_gthr(1:sum( vector_sizes_int_gthr)))
       allocate(vector_real_gthr(1:sum(vector_sizes_real_gthr)))
    endif

    call mpi_gatherv(vector_int_mpi,vector_sizes(1),mpi_integer,&
         vector_int_gthr,vector_sizes_int_gthr,vector_disp_int_gthr,&
         mpi_integer,0,mpi_comm_world,ierr_mpi)

    call mpi_gatherv(vector_real_mpi,vector_sizes(2),&
         mpi_double_precision,vector_real_gthr,vector_sizes_real_gthr,&
         vector_disp_real_gthr,mpi_double_precision,0,mpi_comm_world,&
         ierr_mpi)

    if(mpi_taskid.eq.0) then
       allocate(ierrors_out(1:nrows))
       allocate(Z_out(1:nrows))
       allocate(N_out(1:nrows))
       allocate(Q20_in(1:nrows))
       allocate(beta_in(1:nrows))
       allocate(E_HFB_out(1:nrows))
       allocate(Q20Z_out(1:nrows))
       allocate(Q20N_out(1:nrows))
       allocate(Q20T_out(1:nrows))
       do i = 1,mpi_size
          do j = 1,vector_sizes_int_gthr(mod(i,mpi_size))/n_int_masstable_out
             ierrors_out(mpi_size*(j-1)+i) = vector_int_gthr(vector_disp_int_gthr(mod(i,mpi_size))+n_int_masstable_out*(j-1)+1)
                   Z_out(mpi_size*(j-1)+i) = vector_int_gthr(vector_disp_int_gthr(mod(i,mpi_size))+n_int_masstable_out*(j-1)+2)
                   N_out(mpi_size*(j-1)+i) = vector_int_gthr(vector_disp_int_gthr(mod(i,mpi_size))+n_int_masstable_out*(j-1)+3)
          enddo
          do j = 1,vector_sizes_real_gthr(mod(i,mpi_size))/n_real_masstable_out
                Q20_in(mpi_size*(j-1)+i) = vector_real_gthr(vector_disp_real_gthr(mod(i,mpi_size))+n_real_masstable_out*(j-1)+1)
               beta_in(mpi_size*(j-1)+i) = vector_real_gthr(vector_disp_real_gthr(mod(i,mpi_size))+n_real_masstable_out*(j-1)+2)
             E_HFB_out(mpi_size*(j-1)+i) = vector_real_gthr(vector_disp_real_gthr(mod(i,mpi_size))+n_real_masstable_out*(j-1)+3)
              Q20Z_out(mpi_size*(j-1)+i) = vector_real_gthr(vector_disp_real_gthr(mod(i,mpi_size))+n_real_masstable_out*(j-1)+4)
              Q20N_out(mpi_size*(j-1)+i) = vector_real_gthr(vector_disp_real_gthr(mod(i,mpi_size))+n_real_masstable_out*(j-1)+5)
              Q20T_out(mpi_size*(j-1)+i) = vector_real_gthr(vector_disp_real_gthr(mod(i,mpi_size))+n_real_masstable_out*(j-1)+6)
          enddo
       enddo
    endif

  end subroutine gather_results
#endif
#if(DO_PES==1)
  !=======================================================================
  !>  After the master process has read the binary file, it constructs (allocates
  !>  and fills out) a set of vectors containing all the data to be broadcast to
  !>  the other processes.
  !>  - vector_log_mpi contains all the booleans to broadcast
  !>  - vector_int_mpi contains all the integers to broadcast
  !>  - vector_real_mpi  contains all the reals to broadcast
  !>  Sizes of each vector must also be broadcast since other MPI processes don't
  !>  know about it before
  !=======================================================================
  Subroutine pack()
    Implicit None
    Integer :: i, k, ihl, n_bool, n_int, n_real

    ! Expected sizes are used for allocation
    If(Trim(functional) == 'READ') Then
       n_bool = 6
       n_int  = 4 + 6*(bloall+1)
       n_real = 68 + lambdaMax + nghl*(lambdaMax+1) + 3*2*nghl + 24*nghl + 2*(bloall+1) + 2*nqp
    Else
       n_bool = 0
       n_int  = 4 + 6*(bloall+1)
       n_real = 34 + lambdaMax + nghl*(lambdaMax+1) + 3*2*nghl + 24*nghl + 2*(bloall+1) + 2*nqp
    End If
    If(Trim(functional) == 'READ') Then
       If(Allocated(vector_log_mpi)) Deallocate(vector_log_mpi)
       Allocate(vector_log_mpi(1:n_bool))
    End If
    If(Allocated(vector_int_mpi)) Deallocate(vector_int_mpi,vector_real_mpi)
    Allocate(vector_int_mpi(1:n_int),vector_real_mpi(1:n_real))
    If(Allocated(vector_sizes)) Deallocate(vector_sizes)
    Allocate(vector_sizes(1:3))
    vector_sizes(1) = n_bool
    vector_sizes(2) = n_int
    vector_sizes(3) = n_real

    ! BOOLEANS
    i = 0
    ! 'SkyFunct'
    If(Trim(functional) == 'READ') Then
       i = i+1
       vector_log_mpi(i) = use_INM
       i = i+1
       vector_log_mpi(i) = use_cm_cor
       i = i+1
       vector_log_mpi(i) = use_j2terms
       i = i+1
       vector_log_mpi(i) = force_is_dme
       i = i+1
       vector_log_mpi(i) = finite_range
       i = i+1
       vector_log_mpi(i) = hb0_charge_dependent
    End If
    ! Total number of booleans (converted to integer for broadcast)
    If(debug > 0) Write(6,'("n_bool = ",i6," expected = ",i6)') i,n_bool

    ! INTEGERS
    i = 0
    ! 'Blocking'
    Do k=0,bloall
       i = i+1
       vector_int_mpi(i) = bloblo(k,1)  ! (bloall+1, 2)
       i = i+1
       vector_int_mpi(i) = bloblo(k,2)
       i = i+1
       vector_int_mpi(i) = blo123(k,1)  ! (bloall+1, 2)
       i = i+1
       vector_int_mpi(i) = blo123(k,2)
       i = i+1
       vector_int_mpi(i) = blok1k2(k,1) ! (bloall+1, 2)
       i = i+1
       vector_int_mpi(i) = blok1k2(k,2)
    End Do
    i = i+1
    vector_int_mpi(i) = blomax(1) ! (2, 1)
    i = i+1
    vector_int_mpi(i) = blomax(2)
    ! 'Temperat'
    i = i+1
    vector_int_mpi(i) = nqp !  (1, 1)
    i = i+1
    vector_int_mpi(i) = nuv !  (1, 1)
    ! Total number of integers
    If(debug > 0) Write(6,'("n_int = ",i6," expected = ",i6)') i, n_int

    ! REALS
    i = 0
    ! 'SkyFunct'
    If(Trim(functional) == 'READ') Then
       i = i+1
       vector_real_mpi(i) = pwi       ! (1, 1)
       i = i+1
       vector_real_mpi(i) = E_NM      ! (1, 1)
       i = i+1
       vector_real_mpi(i) = K_NM      ! (1, 1)
       i = i+1
       vector_real_mpi(i) = SMASS_NM  ! (1, 1)
       i = i+1
       vector_real_mpi(i) = RHO_NM    ! (1, 1)
       i = i+1
       vector_real_mpi(i) = ASS_NM    ! (1, 1)
       i = i+1
       vector_real_mpi(i) = LASS_NM   ! (1, 1)
       i = i+1
       vector_real_mpi(i) = VMASS_NM  ! (1, 1)
       i = i+1
       vector_real_mpi(i) = P_NM      ! (1, 1)
       i = i+1
       vector_real_mpi(i) = KA_NM     ! (1, 1)
       i = i+1
       vector_real_mpi(i) = Crho(0)   ! (2, 1)
       i = i+1
       vector_real_mpi(i) = Crho(1)
       i = i+1
       vector_real_mpi(i) = Cdrho(0)  ! (2, 1)
       i = i+1
       vector_real_mpi(i) = Cdrho(1)
       i = i+1
       vector_real_mpi(i) = Ctau(0)   ! (2, 1)
       i = i+1
       vector_real_mpi(i) = Ctau(1)
       i = i+1
       vector_real_mpi(i) = CrDr(0)   ! (2, 1)
       i = i+1
       vector_real_mpi(i) = CrDr(1)
       i = i+1
       vector_real_mpi(i) = CrdJ(0)   ! (2, 1)
       i = i+1
       vector_real_mpi(i) = CrdJ(1)
       i = i+1
       vector_real_mpi(i) = CJ(0)     ! (2, 1)
       i = i+1
       vector_real_mpi(i) = CJ(1)
       i = i+1
       vector_real_mpi(i) = CpV0(0)   ! (2, 1)
       i = i+1
       vector_real_mpi(i) = CpV0(1)
       i = i+1
       vector_real_mpi(i) = CpV1(0)   ! (2, 1)
       i = i+1
       vector_real_mpi(i) = CpV1(1)
       i = i+1
       vector_real_mpi(i) = sigma     ! (1, 1)
       i = i+1
       vector_real_mpi(i) = hbzero    ! (1, 1)
       i = i+1
       vector_real_mpi(i) = hb0       ! (1, 1)
       i = i+1
       vector_real_mpi(i) = hb0n      ! (1, 1)
       i = i+1
       vector_real_mpi(i) = hb0p      ! (1, 1)
       ! 'HO-Basis'
       i = i+1
       vector_real_mpi(i) = b0 ! (1, 1)
       i = i+1
       vector_real_mpi(i) = bz ! (1, 1)
       i = i+1
       vector_real_mpi(i) = bp ! (1, 1)
    End If
    ! 'Various.'
    i = i+1
    vector_real_mpi(i) = si            ! (1, 1)
    i = i+1
    vector_real_mpi(i) = etot          ! (1, 1)
    i = i+1
    vector_real_mpi(i) = rms(1)        ! (3, 1)
    i = i+1
    vector_real_mpi(i) = rms(2)
    i = i+1
    vector_real_mpi(i) = rms(3)
    i = i+1
    vector_real_mpi(i) = bet           ! (1, 1)
    i = i+1
    vector_real_mpi(i) = pwi           ! (1, 1)
    i = i+1
    vector_real_mpi(i) = xmix          ! (1, 1)
    i = i+1
    vector_real_mpi(i) = del(1)        ! (2, 1)
    i = i+1
    vector_real_mpi(i) = del(2)
    i = i+1
    vector_real_mpi(i) = ept(1)        ! (3, 1)
    i = i+1
    vector_real_mpi(i) = ept(2)
    i = i+1
    vector_real_mpi(i) = ept(3)
    i = i+1
    vector_real_mpi(i) = ala(1)        ! (2, 1)
    i = i+1
    vector_real_mpi(i) = ala(2)
    i = i+1
    vector_real_mpi(i) = ala2(1)       ! (2, 1)
    i = i+1
    vector_real_mpi(i) = ala2(2)
    i = i+1
    vector_real_mpi(i) = alast(1)      ! (2, 1)
    i = i+1
    vector_real_mpi(i) = alast(2)
    i = i+1
    vector_real_mpi(i) = tz(1)         ! (2, 1)
    i = i+1
    vector_real_mpi(i) = tz(2)
    i = i+1
    vector_real_mpi(i) = varmas        ! (1, 1)
    i = i+1
    vector_real_mpi(i) = varmasNZ(1)   ! (2, 1)
    i = i+1
    vector_real_mpi(i) = varmasNZ(2)
    i = i+1
    vector_real_mpi(i) = pjmassNZ(1)   ! (2, 1)
    i = i+1
    vector_real_mpi(i) = pjmassNZ(2)
    i = i+1
    vector_real_mpi(i) = ass(1)        ! (2, 1)
    i = i+1
    vector_real_mpi(i) = ass(2)
    i = i+1
    vector_real_mpi(i) = skass         ! (1, 1)
    ! 'Constrai'
    Do k=1,lambdaMax
       i = i+1
       vector_real_mpi(i) = multLag(k)       ! (lambdaMax, 1)
    End Do
    i = i+1
    vector_real_mpi(i) = neckLag             ! (1, 1)
    Do k=0,lambdaMax
       Do ihl=1,nghl
          i = i+1
          vector_real_mpi(i) = qfield(ihl,k)   ! (nghl, lambdaMax+1)
       End Do
    End Do
    ! 'Densits.
    Do k=1,nghl
       i = i+1
       vector_real_mpi(i) = ro(k,1)    ! (nghl, 2)
       i = i+1
       vector_real_mpi(i) = ro(k,2)
       i = i+1
       vector_real_mpi(i) = aka(k,1)   ! (nghl, 2)
       i = i+1
       vector_real_mpi(i) = aka(k,2)
    End Do
    ! 'FieldsN.'
    Do k=1,nghl
       i = i+1
       vector_real_mpi(i) = vn(k)      ! (nghl, 1)
       i = i+1
       vector_real_mpi(i) = vhbn(k)    ! (nghl, 1)
       i = i+1
       vector_real_mpi(i) = vrn(k)     ! (nghl, 1)
       i = i+1
       vector_real_mpi(i) = vzn(k)     ! (nghl, 1)
       i = i+1
       vector_real_mpi(i) = vdn(k)     ! (nghl, 1)
       i = i+1
       vector_real_mpi(i) = vsn(k)     ! (nghl, 1)
       i = i+1
       vector_real_mpi(i) = vSFIZn(k)  ! (nghl, 1)
       i = i+1
       vector_real_mpi(i) = vSZFIn(k)  ! (nghl, 1)
       i = i+1
       vector_real_mpi(i) = vSFIRn(k)  ! (nghl, 1)
       i = i+1
       vector_real_mpi(i) = vSRFIn(k)  ! (nghl, 1)
       i = i+1
       vector_real_mpi(i) = dvn(k)     ! (nghl, 1)
    End Do
    ! 'FieldsP.'
    Do k=1,nghl
       i = i+1
       vector_real_mpi(i) = vp(k)      ! (nghl, 1)
       i = i+1
       vector_real_mpi(i) = vhbp(k)    ! (nghl, 1)
       i = i+1
       vector_real_mpi(i) = vrp(k)     ! (nghl, 1)
       i = i+1
       vector_real_mpi(i) = vzp(k)     ! (nghl, 1)
       i = i+1
       vector_real_mpi(i) = vdp(k)     ! (nghl, 1)
       i = i+1
       vector_real_mpi(i) = vsp(k)     ! (nghl, 1)
       i = i+1
       vector_real_mpi(i) = vSFIZp(k)  ! (nghl, 1)
       i = i+1
       vector_real_mpi(i) = vSZFIp(k)  ! (nghl, 1)
       i = i+1
       vector_real_mpi(i) = vSFIRp(k)  ! (nghl, 1)
       i = i+1
       vector_real_mpi(i) = vSRFIp(k)  ! (nghl, 1)
       i = i+1
       vector_real_mpi(i) = dvp(k)     ! (nghl, 1)
    End Do
    ! 'Blocking'
    Do k=0,bloall
       i = i+1
       vector_real_mpi(i) = bloqpdif(k,1) ! (bloall+1, 2)
       i = i+1
       vector_real_mpi(i) = bloqpdif(k,2)
    End Do
    ! 'Temperat'
    i = i+1
    vector_real_mpi(i) = temper      ! (1, 1)
    i = i+1
    vector_real_mpi(i) = entropy(1)  ! (3, 1)
    i = i+1
    vector_real_mpi(i) = entropy(2)
    i = i+1
    vector_real_mpi(i) = entropy(3)
    Do k=1,nqp
       i = i+1
       vector_real_mpi(i) = fp_T(k)  ! (nqp, 1)
       i = i+1
       vector_real_mpi(i) = fn_T(k)  ! (nqp, 1)
    End Do
    ! 'Regular.'
    Do k=1,nghl
       i = i+1
       vector_real_mpi(i) = MEFFn(k)       ! (nghl, 1)
       i = i+1
       vector_real_mpi(i) = MEFFp(k)       ! (nghl, 1)
       i = i+1
       vector_real_mpi(i) = geff_inv(k,1)  ! (nghl, 2)
       i = i+1
       vector_real_mpi(i) = geff_inv(k,2)
    End Do
    ! Total number of real numbers
    If(debug > 0) Write(6,'("n_real = ",i6," expected = ",i6)') i, n_real

  End Subroutine pack
  !=======================================================================
  !>  After the master process has read the binary file, it constructs (allocates
  !>  and fills out) a set of vectors containing all the data to be broadcast to
  !>  the other processes.
  !>  - vector_log_mpi contains all the booleans to broadcast
  !>  - vector_int_mpi contains all the integers to broadcast
  !>  - vector_real_mpi  contains all the reals to broadcast
  !>  Sizes of each vector must also be broadcast since other MPI processes don't
  !>  know about it before
  !=======================================================================
  Subroutine unpack()
    Implicit None
    Integer :: i, k,ihl,n_bool,n_int,n_real

    ! BOOLEANS
    i = 0
    ! 'SkyFunct'
    If(Trim(functional) == 'READ') Then
       i = i+1
       use_INM              = vector_log_mpi(i)
       i = i+1
       use_cm_cor           = vector_log_mpi(i)
       i = i+1
       use_j2terms          = vector_log_mpi(i)
       i = i+1
       force_is_dme         = vector_log_mpi(i)
       i = i+1
       finite_range         = vector_log_mpi(i)
       i = i+1
       hb0_charge_dependent = vector_log_mpi(i)
    End If

    ! INTEGERS
    i = 0
    ! 'Blocking'
    Do k=0,bloall
       i = i+1
       bloblo(k,1)  = vector_int_mpi(i) ! (bloall+1, 2)
       i = i+1
       bloblo(k,2)  = vector_int_mpi(i)
       i = i+1
       blo123(k,1)  = vector_int_mpi(i) ! (bloall+1, 2)
       i = i+1
       blo123(k,2)  = vector_int_mpi(i)
       i = i+1
       blok1k2(k,1) = vector_int_mpi(i) ! (bloall+1, 2)
       i = i+1
       blok1k2(k,2) = vector_int_mpi(i)
    End Do
    i = i+1
    blomax(1) = vector_int_mpi(i) ! (2, 1)
    i = i+1
    blomax(2) = vector_int_mpi(i)
    ! 'Temperat'
    i = i+1
    nqp = vector_int_mpi(i) !  (1, 1)
    i = i+1
    nuv = vector_int_mpi(i) !  (1, 1)

    ! REALS
    i = 0
    ! 'SkyFunct'
    If(Trim(functional) == 'READ') Then
       i = i+1
       pwi      = vector_real_mpi(i) ! (1, 1)
       i = i+1
       E_NM     = vector_real_mpi(i) ! (1, 1)
       i = i+1
       K_NM     = vector_real_mpi(i) ! (1, 1)
       i = i+1
       SMASS_NM = vector_real_mpi(i) ! (1, 1)
       i = i+1
       RHO_NM   = vector_real_mpi(i) ! (1, 1)
       i = i+1
       ASS_NM   = vector_real_mpi(i) ! (1, 1)
       i = i+1
       LASS_NM  = vector_real_mpi(i) ! (1, 1)
       i = i+1
       VMASS_NM = vector_real_mpi(i) ! (1, 1)
       i = i+1
       P_NM     = vector_real_mpi(i) ! (1, 1)
       i = i+1
       KA_NM    = vector_real_mpi(i) ! (1, 1)
       i = i+1
       Crho(0)  = vector_real_mpi(i) ! (2, 1)
       i = i+1
       Crho(1)  = vector_real_mpi(i)
       i = i+1
       Cdrho(0) = vector_real_mpi(i) ! (2, 1)
       i = i+1
       Cdrho(1) = vector_real_mpi(i)
       i = i+1
       Ctau(0)  = vector_real_mpi(i) ! (2, 1)
       i = i+1
       Ctau(1)  = vector_real_mpi(i)
       i = i+1
       CrDr(0)  = vector_real_mpi(i) ! (2, 1)
       i = i+1
       CrDr(1)  = vector_real_mpi(i)
       i = i+1
       CrdJ(0)  = vector_real_mpi(i) ! (2, 1)
       i = i+1
       CrdJ(1)  = vector_real_mpi(i)
       i = i+1
       CJ(0)    = vector_real_mpi(i) ! (2, 1)
       i = i+1
       CJ(1)    = vector_real_mpi(i)
       i = i+1
       CpV0(0)  = vector_real_mpi(i) ! (2, 1)
       i = i+1
       CpV0(1)  = vector_real_mpi(i)
       i = i+1
       CpV1(0)  = vector_real_mpi(i) ! (2, 1)
       i = i+1
       CpV1(1)  = vector_real_mpi(i)
       i = i+1
       sigma     = vector_real_mpi(i) ! (1, 1)
       i = i+1
       hbzero    = vector_real_mpi(i) ! (1, 1)
       i = i+1
       hb0       = vector_real_mpi(i) ! (1, 1)
       i = i+1
       hb0n      = vector_real_mpi(i) ! (1, 1)
       i = i+1
       hb0p      = vector_real_mpi(i) ! (1, 1)
       ! 'HO-Basis'
       i = i+1
       b0 = vector_real_mpi(i)! (1, 1)
       i = i+1
       bz = vector_real_mpi(i)! (1, 1)
       i = i+1
       bp = vector_real_mpi(i)! (1, 1)
    End If
    ! 'Various.'
    i = i+1
    si          = vector_real_mpi(i) ! (1, 1)
    siold=si
    i = i+1
    etot        = vector_real_mpi(i) ! (1, 1)
    i = i+1
    rms(1)      = vector_real_mpi(i) ! (3, 1)
    i = i+1
    rms(2)      = vector_real_mpi(i)
    i = i+1
    rms(3)      = vector_real_mpi(i)
    i = i+1
    bet         = vector_real_mpi(i) ! (1, 1)
    i = i+1
    xmix        = vector_real_mpi(i) ! (1, 1)
    i = i+1
    pwi         = vector_real_mpi(i) ! (1, 1)
    i = i+1
    del(1)      = vector_real_mpi(i) ! (2, 1)
    i = i+1
    del(2)      = vector_real_mpi(i)
    i = i+1
    ept(1)      = vector_real_mpi(i) ! (3, 1)
    i = i+1
    ept(2)      = vector_real_mpi(i)
    i = i+1
    ept(3)      = vector_real_mpi(i)
    i = i+1
    ala(1)      = vector_real_mpi(i) ! (2, 1)
    i = i+1
    ala(2)      = vector_real_mpi(i)
    i = i+1
    ala2(1)     = vector_real_mpi(i) ! (2, 1)
    i = i+1
    ala2(2)     = vector_real_mpi(i)
    i = i+1
    alast(1)    = vector_real_mpi(i) ! (2, 1)
    i = i+1
    alast(2)    = vector_real_mpi(i)
    i = i+1
    tz(1)       = vector_real_mpi(i) ! (2, 1)
    i = i+1
    tz(2)       = vector_real_mpi(i)
    i = i+1
    varmas      = vector_real_mpi(i) ! (1, 1)
    i = i+1
    varmasNZ(1) = vector_real_mpi(i) ! (2, 1)
    i = i+1
    varmasNZ(2) = vector_real_mpi(i)
    i = i+1
    pjmassNZ(1) = vector_real_mpi(i) ! (2, 1)
    i = i+1
    pjmassNZ(2) = vector_real_mpi(i)
    i = i+1
    ass(1)      = vector_real_mpi(i) ! (2, 1)
    i = i+1
    ass(2)      = vector_real_mpi(i)
    i = i+1
    skass       = vector_real_mpi(i) ! (1, 1)
    ! 'Constrai'
    Do k=1,lambdaMax
       i = i+1
       multLag(k) = vector_real_mpi(i)      ! (lambdaMax, 1)
    End Do
    i = i+1
    neckLag = vector_real_mpi(i)            ! (1, 1)
    Do k=0,lambdaMax
       Do ihl=1,nghl
          i = i+1
          qfield(ihl,k) = vector_real_mpi(i)  ! (nghl, lambdaMax+1)
       End Do
    End Do
    ! 'Densits.
    Do k=1,nghl
       i = i+1
       ro(k,1)  = vector_real_mpi(i)  ! (nghl, 2)
       i = i+1
       ro(k,2)  = vector_real_mpi(i)
       i = i+1
       aka(k,1) = vector_real_mpi(i)  ! (nghl, 2)
       i = i+1
       aka(k,2) = vector_real_mpi(i)
    End Do
    ! 'FieldsN.'
    Do k=1,nghl
       i = i+1
       vn(k)     = vector_real_mpi(i) ! (nghl, 1)
       i = i+1
       vhbn(k)   = vector_real_mpi(i) ! (nghl, 1)
       i = i+1
       vrn(k)    = vector_real_mpi(i) ! (nghl, 1)
       i = i+1
       vzn(k)    = vector_real_mpi(i) ! (nghl, 1)
       i = i+1
       vdn(k)    = vector_real_mpi(i) ! (nghl, 1)
       i = i+1
       vsn(k)    = vector_real_mpi(i) ! (nghl, 1)
       i = i+1
       vSFIZn(k) = vector_real_mpi(i) ! (nghl, 1)
       i = i+1
       vSZFIn(k) = vector_real_mpi(i) ! (nghl, 1)
       i = i+1
       vSFIRn(k) = vector_real_mpi(i) ! (nghl, 1)
       i = i+1
       vSRFIn(k) = vector_real_mpi(i) ! (nghl, 1)
       i = i+1
       dvn(k)    = vector_real_mpi(i) ! (nghl, 1)
    End Do
    ! 'FieldsP.'
    Do k=1,nghl
       i = i+1
       vp(k)     = vector_real_mpi(i) ! (nghl, 1)
       i = i+1
       vhbp(k)   = vector_real_mpi(i) ! (nghl, 1)
       i = i+1
       vrp(k)    = vector_real_mpi(i) ! (nghl, 1)
       i = i+1
       vzp(k)    = vector_real_mpi(i) ! (nghl, 1)
       i = i+1
       vdp(k)    = vector_real_mpi(i) ! (nghl, 1)
       i = i+1
       vsp(k)    = vector_real_mpi(i) ! (nghl, 1)
       i = i+1
       vSFIZp(k) = vector_real_mpi(i) ! (nghl, 1)
       i = i+1
       vSZFIp(k) = vector_real_mpi(i) ! (nghl, 1)
       i = i+1
       vSFIRp(k) = vector_real_mpi(i) ! (nghl, 1)
       i = i+1
       vSRFIp(k) = vector_real_mpi(i) ! (nghl, 1)
       i = i+1
       dvp(k)    = vector_real_mpi(i) ! (nghl, 1)
    End Do
    ! 'Blocking'
    Do k=0,bloall
       i = i+1
       bloqpdif(k,1) = vector_real_mpi(i) ! (bloall+1, 2)
       i = i+1
       bloqpdif(k,2) = vector_real_mpi(i) ! (bloall+1, 2)
    End Do
    ! 'Temperat'
    i = i+1
    temper     = vector_real_mpi(i) ! (1, 1)
    i = i+1
    entropy(1) = vector_real_mpi(i) ! (3, 1)
    i = i+1
    entropy(2) = vector_real_mpi(i)
    i = i+1
    entropy(3) = vector_real_mpi(i)
    Do k=1,nqp
       i = i+1
       fp_T(k) = vector_real_mpi(i) ! (nqp, 1)
       i = i+1
       fn_T(k) = vector_real_mpi(i) ! (nqp, 1)
    End Do
    ! 'Regular.'
    Do k=1,nghl
       i = i+1
       MEFFn(k)      = vector_real_mpi(i) ! (nghl, 1)
       i = i+1
       MEFFp(k)      = vector_real_mpi(i) ! (nghl, 1)
       i = i+1
       geff_inv(k,1) = vector_real_mpi(i) ! (nghl, 2)
       i = i+1
       geff_inv(k,2) = vector_real_mpi(i)
    End Do

  End Subroutine unpack
  !=======================================================================
  !>  Within each team, the master process broadcasts first the vector_sizes
  !>  arrays so that the other processes in the team can do all allocations.
  !>  Then the master process broadcasts integers, reals, logicals and
  !>  a single string
  !=======================================================================
  Subroutine broadcast_binary_to_team(n_sizes)
    Use mpi
    Implicit None
    Integer, Intent(In) :: n_sizes
    Integer :: n_bool, n_int, n_real
    !
    If(do_print == 0) Then
       If(Allocated(vector_sizes)) Deallocate(vector_sizes)
       Allocate(vector_sizes(1:n_sizes))
    End If
    Call mpi_bcast(vector_sizes, n_sizes, mpi_integer, 0, COMM_team, ierr_mpi)
    If(do_print == 0) then
       n_bool = vector_sizes(1)
       n_int  = vector_sizes(2)
       n_real = vector_sizes(3)
       If(Trim(functional) == 'READ') Then
          If(Allocated(vector_log_mpi)) Deallocate(vector_log_mpi)
          Allocate(vector_log_mpi(1:n_bool))
       End If
       If(Allocated(vector_int_mpi)) Deallocate(vector_int_mpi,vector_real_mpi)
       Allocate(vector_int_mpi(1:n_int),vector_real_mpi(1:n_real))
    End If
    If(Trim(functional) == 'READ') Then
       Call mpi_bcast(vector_log_mpi, vector_sizes(1), mpi_logical, 0, COMM_team, ierr_mpi)
    End If
    Call mpi_bcast(vector_int_mpi, vector_sizes(2), mpi_integer, 0, COMM_team, ierr_mpi)
    Call mpi_bcast(vector_real_mpi,vector_sizes(3), mpi_double_precision, 0, COMM_team, ierr_mpi)
    If(debug > 1) Then
       If(do_print == 1) Then
          Write(6,'("bcast_binary: mpi_taskid = ",i4," team_rank = ",i4)') mpi_taskid,team_rank
          Write(6,'("bcast_binary: n_bool = ",i6," n_int = ",i6," n_real = ",i6)') n_bool,n_int,n_real
        End If
    End If
    !
  End Subroutine broadcast_binary_to_team
#endif

end module HFBTHO_mpi_communication
#endif
