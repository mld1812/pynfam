!***********************************************************************
!
!    Copyright (c) 2016, Lawrence Livermore National Security, LLC.
!                        Produced at the Lawrence Livermore National
!                        Laboratory.
!                        Written by Nicolas Schunck, schunck1@llnl.gov
!
!    LLNL-CODE-728299 All rights reserved.
!    LLNL-CODE-573953 All rights reserved.
!
!    Copyright 2017, R. Navarro Perez, N. Schunck, R. Lasseri, C. Zhang,
!                    J. Sarich
!    Copyright 2012, M.V. Stoitsov, N. Schunck, M. Kortelainen, H.A. Nam,
!                    N. Michel, J. Sarich, S. Wild
!    Copyright 2005, M.V. Stoitsov, J. Dobaczewski, W. Nazarewicz, P.Ring
!
!    This file is part of HFBTHO.
!
!    HFBTHO is free software: you can redistribute it and/or modify it
!    under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    HFBTHO is distributed in the hope that it will be useful, but
!    WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with HFBTHO. If not, see <http://www.gnu.org/licenses/>.
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
!              RESTORATION OF BROKEN SYMMETRIES PACKAGE                !
!                                                                      !
! ==================================================================== !

!----------------------------------------------------------------------
!> This module restores broken symmetries (particle number, angular
!> momentum) of HFB states by employing projection techniques.
!----------------------------------------------------------------------
Module HFBTHO_projections

  Use math
  Use HFBTHO_utilities
  Use HFBTHO_SOLVER

  Implicit None

Contains
  !==============================================================
  ! The main subroutine for performing restoration of symmetries.
  !==============================================================
  Subroutine HFBTHO_restore()
    Implicit None
    Integer(ipr) ibet
    !
    Call initialize_projections()
    Call initialize_wavefunctions()
    !
#if(USE_MPI==2)
    Do ibet=1,maxbet
       If( (ibet-1)/(maxbet/team_size) /= team_rank) Cycle
       Call initialize_angle()
       Call calculate_rotation_matrix(ibet)
       Call calculate_overlaps(ibet)
       Call calculate_densities(ibet)
       Call calculate_energies(ibet)
       Call finalize_angle()
    End Do
    !
    ! Team leader gathers rotated overlaps and energies
    call mpi_barrier(mpi_comm_world,ierr_mpi)
    call mpi_gather(all_overlaps,2*beta_size*maxphi,mpi_double_complex,&
         all_overlaps_gthr,2*beta_size*maxphi,mpi_double_complex,0,COMM_team,ierr_mpi)
    call mpi_gather(all_energies,17*betaphi_size,mpi_double_complex,&
         all_energies_gthr,17*betaphi_size,mpi_double_complex,0,COMM_team,ierr_mpi)
#else
    Do ibet=1,maxbet
       Call initialize_angle()
       Call calculate_rotation_matrix(ibet)
       Call calculate_overlaps(ibet)
       Call calculate_densities(ibet)
       Call calculate_energies(ibet)
       Call finalize_angle()
    End Do
    !
    all_overlaps_gthr=all_overlaps
    all_energies_gthr=all_energies
#endif
    !
    If(team_rank.Eq.0) Then
       Call project()
       Call print_project()
    End If
    !
    Call finalize_projections()
    !
  End Subroutine HFBTHO_restore
  !============================================================
  ! Subroutine that initializes symmetry restoration procedure.
  !============================================================
  Subroutine initialize_projections()
    Implicit None
    Integer(ipr) ibet,it,kcount,ib,nd,kaib,kdib,ki,kf,k1,iphicyl,ih,il,ihil,ihil_iphicyl,iosc1,iosc2,iosc_pair,iit,iphi,iit_iphi,nn_neut,nn_prot
    Real(pr) xmin,xmax,s
    Real(pr),Allocatable:: betamesh(:)
    Integer(ipr) tstart_ibet,tfinish_ibet
    Character(50) :: projfname
    ! This flag is set to False by default. It should be set to True only when the large-scale calculations are carried out, and the number of initially
    ! considered angles needs to be guessed in order to eventually keep and calculate about N_{beta}=30 angles. Therefore, this option overwrites the
    ! 'number_of_rotational_angles' value from the input. The coded 'fit' corresponds to the case of 30 oscillator shells and 1100 basis states cut-off.
    Logical automatic_angle
    Real(pr) nrot
    automatic_angle=.False.
    !
    ! Currently we consider only projections for even-even nuclei
    If(blocking_mode(1).ne.0 .or. blocking_mode(2).ne.0) Then
       Stop 'STOP: Blocking not allowed in projections.'
    End If
    !
    ! Initialize timing variables
    Call system_clock(tstart,clock_rate)
    Call system_clock(tstart_ibet,clock_rate)
    time_ibet=zero
    !
#if(DO_PES==0)
    Write(projfname, '("projections",a4)') '.out'
#elif(DO_PES==1)
    Write(projfname, '("projections",a7,a4)') row_string,'.out'
#endif
    Open(lproj,file=projfname,status='unknown')
    !
  If(automatic_angle) Then
      If(beta0.Lt.0.15d0) Then
         nrot=30.d0
      Else If(beta0.Ge.0.15d0 .And. beta0.Lt.0.673d0) Then
         nrot=12.12d0+102.44d0*beta0
      Else If(beta0.Ge.0.673d0 .And. beta0.Lt.0.867d0) Then
         nrot=1.48d0+112.88d0*beta0
      Else If(beta0.Ge.0.867d0 .And. beta0.Lt.0.929d0) Then
         nrot=100.d0
      Else If(beta0.Ge.0.929d0 .And. beta0.Lt.1.091d0) Then
         nrot=102.d0
      Else If(beta0.Ge.1.091) Then
         nrot=108.d0
      End If
      number_of_rotational_angles=nint(nrot)
      If(Mod(number_of_rotational_angles,2).Eq.1) number_of_rotational_angles=number_of_rotational_angles-1
      If(number_of_rotational_angles.Lt.30) number_of_rotational_angles=30
    End If
    !
    maxphi = number_of_gauge_points
    maxbet = number_of_rotational_angles
    maxj   = maximal_angular_momentum
    maxN   = delta_neutrons
    maxP   = delta_protons
    !
    Write(lproj,'(a)')     '              ===================================================='
    Write(lproj,'(a)')     '              Restoration of broken symmetries module, HFBTHO code'
    Write(lproj,'(a)')     '              ===================================================='
    Write(lproj,'(a44,2x,i2)') '  Number of oscillator shells:              ',n00
    Write(lproj,'(a41,2x,f5.3)') '  Basis deformation beta0:               ',beta0
    Write(lproj,'(a)') '  ---'
    If(AMP_is_on .Eq. 1) Then
       Write(lproj, '(a)') '  Angular momentum projection is on.'
       If(force_parity) Then
          If(Mod(maxbet,2).Eq.0) Then
             maxbet=maxbet/2
          Else
             maxbet=(maxbet+1)/2
          End If
       End If
       Write(lproj,'(a44,x,i3)')   '  Number of rotational angles given:         ', number_of_rotational_angles
       Write(lproj,'(a44,x,i3)')   '  Number of rotational angles considered:    ', maxbet
       Write(lproj,'(a44,x,i3)')   '  Number of rotational angles calculated:    ', maxbet
       Write(lproj,'(a44,x,i3)')   '  Maximal angular momentum considered:       ', maxj
    Else
       Write(lproj, '(a)') '  Angular momentum projection is off.'
       maxbet=1
    End If
    Write(lproj,'(a)') '  ---'
    If(PNP_is_on .Gt. 0) Then
       Write(lproj, '(a)') '  Particle number projection is on.'
       If(Mod(maxphi,2).Eq.0) maxphi=maxphi+1
       maxphi_eff=(maxphi+1)/2
       If(PNP_is_on .Eq. 1) Then
          Write(lproj, '(a)') '  Mixed prescription for non-integer powers of density.'
       Else If(PNP_is_on .eq. 2) Then
          Write(lproj, '(a)') '  Projected prescription in PNP for non-integer powers of density.'
       End If
       Write(lproj,'(a44,x,i3)')   '  Number of gauge angles given:             ', number_of_gauge_points
       Write(lproj,'(a44,x,i3)')   '  Number of gauge angles considered:        ', maxphi
       Write(lproj,'(a44,x,i3)')   '  Number of gauge angles calculated:        ', maxphi_eff
    Else
       Write(lproj, '(a)') '  Particle number projection is off.'
       maxphi=1; maxphi_eff=1
    End If
    Write(lproj,'(a)') '  ---'
    Write(lproj,'(a44,2x,l2)')   '  Parity enforced:                          ', force_parity
    Write(lproj,'(a42,2x,a)')      '  Energy functional:                      ',skyrme
    !-----------------------------------------------------------
    ! Counting the number of qp states for neutrons and protons.
    !-----------------------------------------------------------
    Do it=1,2
       kcount=0
       Do ib=1,nb
          nd=id(ib); kaib=ka(ib,it); kdib=kd(ib,it); ki=kaib+1; kf=kaib+kdib
          Do k1=ki,kf
             kcount=kcount+1
          End Do
       End Do
       kdim(it)=kcount
    End Do
    kdim(3)=nt; kdim(4)=nt
    ! ---------------------------------
    ! Setting up AMP mesh in beta angle
    ! ---------------------------------
    If(maxbet.Eq.1) Then
       Allocate(betabs(1),betaweight(1))
       betabs(1)=zero; betaweight(1)=one
    Else
       If(force_parity) Then
          integration_prefactor=two ! integration from 0 to pi/2
          xmin=0.d0; xmax=1.d0
          Allocate(betabs(maxbet),betamesh(maxbet),betaweight(maxbet))
          Call gauleg(xmin,xmax,betamesh,betaweight,maxbet)
          Do ibet=1,maxbet
             betabs(ibet)=Acos(betamesh(ibet))
          End Do
          Deallocate(betamesh)
       Else
          integration_prefactor=one ! integration from 0 to pi
          xmin=-1.d0; xmax=1.d0
          Allocate(betabs(maxbet),betamesh(maxbet),betaweight(maxbet))
          Call gauleg(xmin,xmax,betamesh,betaweight,maxbet)
          Do ibet=1,maxbet
             betabs(ibet)=Acos(betamesh(ibet))
          End Do
          Deallocate(betamesh)
       End If
    End If
    ! --------------------------------
    ! Setting up PNP mesh in phi angle
    ! --------------------------------
    If(Mod(maxN,2).Ne.0) maxN=maxN+1 ! Number of neutrons is even
    If(Mod(maxP,2).Ne.0) maxP=maxP+1 ! Number of protons is even
    Allocate(phiabs(maxphi),ephi(maxphi),ephic(maxphi),ephicN(maxphi,-maxN/2:maxN/2),ephicP(maxphi,-maxP/2:maxP/2))
    Do iphi=1,maxphi
       phiabs(iphi)=(iphi-1)*pi/maxphi
       ephi(iphi)=Cmplx(Cos(phiabs(iphi)),Sin(phiabs(iphi)))
       ephic(iphi)=Cmplx(Cos(phiabs(iphi)),-Sin(phiabs(iphi)))
       Do nn_neut=-maxN/2,maxN/2
          ephicN(iphi,nn_neut)=Cmplx(Cos(phiabs(iphi)*(tz(1)+two*nn_neut)),-Sin(phiabs(iphi)*(tz(1)+two*nn_neut)))
       End Do
       Do nn_prot=-maxP/2,maxP/2
          ephicP(iphi,nn_prot)=Cmplx(Cos(phiabs(iphi)*(tz(2)+two*nn_prot)),-Sin(phiabs(iphi)*(tz(2)+two*nn_prot)))
       End Do
    End Do
    !------------------------------------------------------------------
    ! Allocating memory for different global arrays, initializing them.
    !------------------------------------------------------------------
    Allocate(rotated_overlap(maxbet,maxphi,2),detR(maxbet))
    Allocate(projected_overlap(0:maxj),projected_ekinN(0:maxj),projected_ekinP(0:maxj),projected_ecodi(0:maxj),projected_ecoex(0:maxj),&
             projected_EVOL_rho_tau(0:maxj),projected_EVOL_rho_rho(0:maxj),projected_ESURF_rho_drho(0:maxj),projected_ETENS(0:maxj),&
             projected_ESO_rho_nablaj(0:maxj),projected_eptN(0:maxj),projected_eptP(0:maxj),projected_xn1(0:maxj),projected_xn2(0:maxj),&
             projected_rms1(0:maxj),projected_rms2(0:maxj),projected_delN(0:maxj),projected_delP(0:maxj),&
             projected_NP(0:maxj,-maxN/2:maxN/2,-maxP/2:maxP/2),projected_NP_norm(0:maxj))
    Allocate(nz_sim(nt),nr_sim(nt),nl_sim(nt))
    Allocate(phicyl(ngphi))
    Allocate(ihil_convert(nghl*ngphi),iphicyl_convert(nghl*ngphi),ihil_iphicyl_convert(nghl,ngphi))
    Allocate(beta_active(maxbet))
    Allocate(iosc1_pair(nt*(nt+1)/2),iosc2_pair(nt*(nt+1)/2))
    Allocate(itiphi_pair1(4*maxphi_eff),itiphi_pair2(4*maxphi_eff))
    !
    rotated_overlap=czero;detR=czero
    projected_overlap=zero;projected_ekinN=zero;projected_ekinP=zero;projected_ecodi=zero;projected_ecoex=zero;projected_EVOL_rho_tau=zero;
    projected_EVOL_rho_rho=zero;projected_ESURF_rho_drho=zero;projected_ETENS=zero;projected_ESO_rho_nablaj=zero;projected_eptN=zero;
    projected_eptP=zero;projected_xn1=zero;projected_xn2=zero;projected_rms1=zero;projected_rms2=zero;
    projected_delN=czero;projected_delP=czero;projected_NP=zero;projected_NP_norm=zero
    beta_active=1
    ! ----------------------------
    ! Setting up the simplex basis
    ! ----------------------------
    Call simplex_basis(.False.)
    ! ----------------------------------------
    ! Setting up mesh in cylindrical angle phi
    ! ----------------------------------------
    If(ngphi.Eq.1) Then
       phicyl(1)=zero
    Else
       Do iphicyl=1,ngphi
          phicyl(iphicyl)=(iphicyl-1)*two*pi/(dble(ngphi)-one)
       End Do
    End If
    phicyl_integration_step=two*pi/(dble(ngphi)-one)
    ! ----------------------------------------------------------------
    ! Setting up conversion arrays for (z,perp) x (phicyl) coordinates
    ! ----------------------------------------------------------------
    ihil_iphicyl=0
    Do iphicyl=1,ngphi
       Do ihil=1,nghl
          ihil_iphicyl=ihil_iphicyl+1
          !
          ihil_convert(ihil_iphicyl)=ihil
          iphicyl_convert(ihil_iphicyl)=iphicyl
          !
          ihil_iphicyl_convert(ihil,iphicyl)=ihil_iphicyl
       End Do
    End Do
    ! ------------------------------------
    ! Determining relevant angular momenta
    ! ------------------------------------
    If(force_parity) Then
       jjstep=2
    Else
       jjstep=1
       If(Abs(qmoment(3,3)).Lt.1.d-6) jjstep=2 ! No octupole deformation
    End If
    ! --------------------------------------
    ! Determining pairs of oscillator states
    ! --------------------------------------
    iosc_pair=0
    Do iosc1=1,nt
       Do iosc2=iosc1,nt
          iosc_pair=iosc_pair+1
          iosc1_pair(iosc_pair)=iosc1;iosc2_pair(iosc_pair)=iosc2
       End Do
    End Do
    nt_pair=iosc_pair
    ! --------------------------
    ! Determining iit-iphi pairs
    ! --------------------------
    iit_iphi=0
    Do iit=1,4
       Do iphi=1,maxphi_eff
          iit_iphi=iit_iphi+1
          itiphi_pair1(iit_iphi)=iit
          itiphi_pair2(iit_iphi)=iphi
       End Do
    End Do
    ! ----------------------------------
    ! Setting up xl_ihil, xh_ihil meshes
    ! ----------------------------------
    If(Allocated(xl_ihil)) Deallocate(xl_ihil,xh_ihil)
    Allocate(xl_ihil(nghl),xh_ihil(nghl))
    Do ih=1,ngh
       Do il=1,ngl
          ihil=ih+(il-1)*ngh
          xl_ihil(ihil)=xl(il); xh_ihil(ihil)=xh(ih)
       End Do
    End Do
    ! --------------------
    ! MPI and broadcasting
    ! --------------------
#if(USE_MPI==2)
    beta_size=maxbet/team_size          ! # of beta angles calculated by each single process
    betaphi_size=beta_size*maxphi**2    ! # of beta and phi angles calculated by each single process
#else
    beta_size=maxbet
    betaphi_size=maxbet*maxphi**2
#endif
    !
    If(team_rank.Eq.0) Then
       Do it=1,2
          s=two*sum(ro_normalization(:,it))
          piu(it)=tz(it)/s
       End Do
    End If
    !
#if(USE_MPI==2)
    Call mpi_bcast(piu,2,mpi_double_precision,0,COMM_team,ierr_mpi) 
    Call mpi_bcast(RVqpN,nuv,mpi_double_precision,0,COMM_team,ierr_mpi)
    Call mpi_bcast(RVqpP,nuv,mpi_double_precision,0,COMM_team,ierr_mpi)
    Call mpi_bcast(RUqpN,nuv,mpi_double_precision,0,COMM_team,ierr_mpi)
    Call mpi_bcast(RUqpP,nuv,mpi_double_precision,0,COMM_team,ierr_mpi)
    Call mpi_bcast(KpwiN,nqp,mpi_integer,0,COMM_team,ierr_mpi)
    Call mpi_bcast(KpwiP,nqp,mpi_integer,0,COMM_team,ierr_mpi)
    Call mpi_bcast(ia,nbx,mpi_integer,0,COMM_team,ierr_mpi)
    Call mpi_bcast(id,nbx,mpi_integer,0,COMM_team,ierr_mpi)
    Call mpi_bcast(ka,2*nbx,mpi_integer,0,COMM_team,ierr_mpi)
    Call mpi_bcast(kd,2*nbx,mpi_integer,0,COMM_team,ierr_mpi)
    Call mpi_bcast(kdim,4,mpi_integer,0,COMM_team,ierr_mpi)
#endif
    !
    Allocate(all_overlaps(2*beta_size*maxphi),all_energies(17*betaphi_size))
    all_overlaps=czero; all_energies=czero
    If(team_rank.Eq.0) Then
       Allocate(all_overlaps_gthr(2*beta_size*maxphi*team_size)); all_overlaps_gthr=czero
       Allocate(all_energies_gthr(17*betaphi_size*team_size)); all_energies_gthr=czero
    End If
    ! ------
    ! Timing
    ! ------
    Call system_clock(tfinish_ibet,clock_rate)
    time_ibet(1)=time_ibet(1)+Real(tfinish_ibet-tstart_ibet)/Real(clock_rate)
    !
  End Subroutine initialize_projections
  !=================================================================================
  ! Subroutine that initializes non-rotated U and V matrices in the simplex-y basis.
  !=================================================================================
  Subroutine initialize_wavefunctions()
    Implicit None
    Integer(ipr) it,iit,kcount,kcount_full,i_uv,ib,nd,i0,kaib,kdib,ki,kf,k1,nn,nn1,iosc
    Integer(ipr), Pointer :: KpwiPo(:)
    Real(pr), Pointer     :: VqpPo(:),UqpPo(:)
    Integer(ipr) tstart_ibet,tfinish_ibet
    ! ------
    ! Timing
    ! ------
    Call system_clock(tstart_ibet,clock_rate)
    !
    Allocate(VmatrixN1(1:nt,1:kdim(1)),UmatrixN1(1:nt,1:kdim(1)),VmatrixN2(1:nt,1:kdim(3)),UmatrixN2(1:nt,1:kdim(3)),&
             VmatrixP1(1:nt,1:kdim(2)),UmatrixP1(1:nt,1:kdim(2)),VmatrixP2(1:nt,1:kdim(4)),UmatrixP2(1:nt,1:kdim(4)))
    VmatrixN1=czero; UmatrixN1=czero; VmatrixN2=czero; UmatrixN2=czero
    VmatrixP1=czero; UmatrixP1=czero; VmatrixP2=czero; UmatrixP2=czero
    !
    Do it=1,2
       If(it.Eq.1) Then ! Neutrons
          kcount=0; kcount_full=0; i_uv=0
          Do ib=1,nb ! Loop over blocks
             nd=id(ib); i0=ia(ib); kaib=ka(ib,it); kdib=kd(ib,it); ki=kaib+1; kf=kaib+kdib
             ! Bogoliubov matrices with a cut-off
             Do k1=ki,kf ! Loop over quasiparticle states
                kcount=kcount+1
                Do nn=1,nd ! Loop over basis states within the block ib
                   nn1=KpwiN(k1)+nn; iosc=i0+nn
                   If(nl_sim(iosc).Ge.0) Then ! Simplex-y componets of U and V matrices
                      UmatrixN1(iosc,kcount)=Cmplx(RUqpN(nn1),zero); VmatrixN1(iosc,kcount)=-Cmplx(RVqpN(nn1),zero)
                   Else
                      UmatrixN1(iosc,kcount)=Cmplx(zero,RUqpN(nn1)); VmatrixN1(iosc,kcount)=-Cmplx(zero,RVqpN(nn1))
                   End If
                End Do ! nn
             End Do ! k1
             ! Full Bogoliubov matrices (for overlap calculation)
             Do k1=1,nd
                kcount_full=kcount_full+1
                Do nn=1,nd
                   i_uv=i_uv+1; iosc=i0+nn
                   If(nl_sim(iosc).Ge.0) Then ! Simplex-y componets of U and V matrices
                      UmatrixN2(iosc,kcount_full)=Cmplx(RUqpN(i_uv),zero); VmatrixN2(iosc,kcount_full)=-Cmplx(RVqpN(i_uv),zero)
                   Else
                      UmatrixN2(iosc,kcount_full)=Cmplx(zero,RUqpN(i_uv)); VmatrixN2(iosc,kcount_full)=-Cmplx(zero,RVqpN(i_uv))
                   End If
                End Do
             End Do
          End Do ! ib
          !
       Else ! Protons
          kcount=0; kcount_full=0; i_uv=0
          Do ib=1,nb ! Loop over blocks
             nd=id(ib); i0=ia(ib); kaib=ka(ib,it); kdib=kd(ib,it); ki=kaib+1; kf=kaib+kdib
             ! Bogoliubov matrices with a cut-off
             Do k1=ki,kf ! Loop over quasiparticle states
                kcount=kcount+1
                Do nn=1,nd ! Loop over basis states within the block ib
                   nn1=KpwiP(k1)+nn; iosc=i0+nn
                   If(nl_sim(iosc).Ge.0) Then ! Simplex-y componets of U and V matrices
                      UmatrixP1(iosc,kcount)=Cmplx(RUqpP(nn1),zero); VmatrixP1(iosc,kcount)=-Cmplx(RVqpP(nn1),zero)
                   Else
                      UmatrixP1(iosc,kcount)=Cmplx(zero,RUqpP(nn1)); VmatrixP1(iosc,kcount)=-Cmplx(zero,RVqpP(nn1))
                   End If
                End Do ! nn
             End Do ! k1
             ! Full Bogoliubov matrices (for overlap calculation)
             Do k1=1,nd
                kcount_full=kcount_full+1
                Do nn=1,nd
                   i_uv=i_uv+1; iosc=i0+nn
                   If(nl_sim(iosc).Ge.0) Then ! Simplex-y componets of U and V matrices
                      UmatrixP2(iosc,kcount_full)=Cmplx(RUqpP(i_uv),zero); VmatrixP2(iosc,kcount_full)=-Cmplx(RVqpP(i_uv),zero)
                   Else
                      UmatrixP2(iosc,kcount_full)=Cmplx(zero,RUqpP(i_uv)); VmatrixP2(iosc,kcount_full)=-Cmplx(zero,RVqpP(i_uv))
                   End If
                End Do
             End Do
          End Do ! ib
          !
       End If
    End Do ! it
    ! ------
    ! Timing
    ! ------
    Call system_clock(tfinish_ibet,clock_rate)
    time_ibet(1)=time_ibet(1)+Real(tfinish_ibet-tstart_ibet)/Real(clock_rate)
    !
  End Subroutine initialize_wavefunctions
  !===============================================================
  ! Subroutine that initializes quantities specific to each angle.
  !===============================================================
  Subroutine initialize_angle()
    Implicit None
    Integer(ipr) iosc,it,iphi,isimplex
    !-----------------------------
    ! Allocating rotation matrices
    !-----------------------------
    Allocate(rotation_matrix(nt,nt),inverse_rotation_matrix(nt,nt))
    Allocate(rotated_density(2,maxphi_eff,2),rotated_kappa(2,maxphi_eff,2))
    Do it=1,2
       Do iphi=1,maxphi_eff
          Do isimplex=1,2
             Allocate(rotated_density(it,iphi,isimplex)%arr(1:nt,1:nt)); rotated_density(it,iphi,isimplex)%arr(1:nt,1:nt)=czero
             Allocate(rotated_kappa(it,iphi,isimplex)%arr(1:nt,1:nt)); rotated_kappa(it,iphi,isimplex)%arr(1:nt,1:nt)=czero
          End Do
       End Do
    End Do
    !
    Allocate(rotated_ro(maxphi_eff,nghl*ngphi,2),rotated_tau(maxphi_eff,nghl*ngphi,2),rotated_dj(maxphi_eff,nghl*ngphi,2),rotated_dro(maxphi_eff,nghl*ngphi,2),rotated_sfiz(maxphi_eff,nghl*ngphi,2),&
             rotated_sfir(maxphi_eff,nghl*ngphi,2),rotated_srfi(maxphi_eff,nghl*ngphi,2),rotated_szfi(maxphi_eff,nghl*ngphi,2),rotated_szz(maxphi_eff,nghl*ngphi,2),rotated_srz(maxphi_eff,nghl*ngphi,2),&
             rotated_srr(maxphi_eff,nghl*ngphi,2),rotated_szr(maxphi_eff,nghl*ngphi,2),rotated_sfifi(maxphi_eff,nghl*ngphi,2),rotated_aka(maxphi_eff,nghl*ngphi,2),ro_projected(nghl*ngphi,2))
    Allocate(iosc1_contributing(nt*nt,maxphi_eff),iosc2_contributing(nt*nt,maxphi_eff),nt_contributing(maxphi_eff))
    Allocate(cou_projected(maxphi_eff,nghl))
    rotation_matrix=czero;inverse_rotation_matrix=czero
    rotated_ro=czero;rotated_tau=czero;rotated_dj=czero;rotated_dro=czero;rotated_sfiz=czero;rotated_sfir=czero;rotated_srfi=czero;
    rotated_szfi=czero;rotated_szz=czero;rotated_srz=czero;rotated_srr=czero;rotated_szr=czero;rotated_sfifi=czero;rotated_aka=czero;ro_projected=czero
    iosc1_contributing=0;iosc2_contributing=0;nt_contributing=0
    !
    Do iosc=1,nt
       inverse_rotation_matrix(iosc,iosc)=cone
    End Do
    !
  End Subroutine initialize_angle
  !=======================================================================================
  ! Subroutine that calculates the rotation matrix exp(-i*beta*Jy) in the simplex-y basis.
  !=======================================================================================
  Subroutine calculate_rotation_matrix(ibet)
    Use math
    Implicit None
    Integer(ipr) ibet,iosc1,iosc2,iosc_pair,nza1,nra1,nla1,nza2,nra2,nla2,ifl
    Integer(ipr) tstart_ibet,tfinish_ibet
    Real(pr)  beta,rotel1,rotel2,rotel3,rotel4
    Real(pr), Parameter:: tol_one=one-1.d-8, tol_zero=1.d-8
    Complex(pr) det
    Integer(ipr) info
    ! ------
    ! Timing
    ! ------
    Call system_clock(tstart_ibet,clock_rate)
    !--------------------------
    ! Defining rotational angle
    !--------------------------
    beta=betabs(ibet)
    !-----------------------------
    ! Loop over oscillator states
    !-----------------------------
!$OMP  PARALLEL DO       &
!$OMP& DEFAULT(NONE)     &
!$OMP& SCHEDULE(DYNAMIC) &
!$OMP& SHARED(nt_pair,beta,bp,bz,nz_sim,nr_sim,nl_sim,iv,rotation_matrix,iosc1_pair,iosc2_pair) &
!$OMP& PRIVATE(iosc_pair,iosc1,nza1,nra1,nla1,iosc2,nza2,nra2,nla2,rotel1,rotel2,rotel3,rotel4)
    Do iosc_pair=1,nt_pair
       !
       iosc1=iosc1_pair(iosc_pair); nza1=nz_sim(iosc1); nra1=nr_sim(iosc1); nla1=nl_sim(iosc1)
       iosc2=iosc2_pair(iosc_pair); nza2=nz_sim(iosc2); nra2=nr_sim(iosc2); nla2=nl_sim(iosc2)
       !
       rotel1=zero; rotel2=zero; rotel3=zero; rotel4=zero
       call calculate_ry(beta,bp,bz,bp,bz,nza1,nra1, nla1,nza2,nra2,nla2,rotel1)
       call calculate_ry(beta,bp,bz,bp,bz,nza1,nra1, nla1,nza2,nra2,-nla2,rotel3)
       !
       rotation_matrix(iosc1,iosc2)=Cmplx(Cos(half*beta)*rotel1,Sin(half*beta)*rotel3)
       rotation_matrix(iosc2,iosc1)=iv(nla1+nla2)*Cmplx(Cos(half*beta)*rotel1,Sin(half*beta)*rotel3) ! R(b,a)=iv(la+lb)*R(a,b)
    End Do ! iosc_pair
!$OMP END PARALLEL DO
    !
    !---------------------------------------------------------------
    ! Tests
    ! ------
    ! Some of the properties that the rotation matrix should satisfy
    ! when the basis is closed under rotation (spherical basis)
    !  1) det(R)=1
    !  2) R(alpha)*R(beta)=R(alpha+beta)
    !  3) (R^{-1})^T = R^* (Robledo, 1994)
    !---------------------------------------------------------------
    inverse_rotation_matrix=rotation_matrix
    !
    Call calculate_inverse(nt,inverse_rotation_matrix,det,info)
    If(info.Eq.0 .and. Abs(det).Gt.1.d-300) Then
       detR(ibet)=det
       beta_active(ibet)=1
    Else
       detR(ibet)=Cmplx(0.,0.)
       beta_active(ibet)=0
    End If
    ! Test if R*R^-1 = identity
!    Call Zgemm('N','N',nt,nt,nt,cone,rotation_matrix,nt,inverse_rotation_matrix,nt,czero,identity,nt)
!    Do iosc1=1,nt
!    Do iosc2=1,nt
!       If(iosc1.Ne.iosc2 .And. Zabs(identity(iosc1,iosc2)).Gt.tol_zero) Then
!          Write(lproj,'(2x,a29,i3,x,i3,a10,f12.8,a1,f13.9)') 'Error in R x R^{-1} element (',iosc1,iosc2,') at angle',betabs(ibet),':',Zabs(identity(iosc1,iosc2))
!       End If
!       If(iosc1.Eq.iosc2 .And. Zabs(identity(iosc1,iosc2)).Lt.tol_one) Then
!          Write(lproj,'(2x,a29,i3,x,i3,a10,f12.8,a1,f13.9)') 'Error in R x R^{-1} element (',iosc1,iosc2,') at angle',betabs(ibet),':',Zabs(identity(iosc1,iosc2))
!       End If
!    End Do
!    End Do
    !
    ! ------
    ! Timing
    ! ------
    Call system_clock(tfinish_ibet,clock_rate)
    time_ibet(2)=time_ibet(2)+Real(tfinish_ibet-tstart_ibet)/Real(clock_rate)
    !
  End Subroutine calculate_rotation_matrix
  !=======================================================================================
  ! Subroutine that calculates overlaps between rotated states and rotated density matrix.
  !=======================================================================================
  Subroutine calculate_overlaps(ibet)
    Implicit None
    Integer(ipr) ibet,it,ifl,k1,iosc1,iosc2,iosc_contributing,iit,iphi,iit_iphi,isimplex
    Integer(ipr) tstart_ibet,tfinish_ibet
    Real(pr), Parameter:: eps=1.d-15
    Complex(pr), Allocatable :: OverlapMatrixTmp(:,:),AuxiliaryMatrix(:,:),AuxiliaryDensity(:,:),aUmatrix(:,:),aVmatrix(:,:)
    Type(ptr_to_cmplx2darray), Allocatable :: rotated_Vmatrix(:,:),rotated_Umatrix(:,:),InverseOverlapMatrix(:,:,:),Vmatrix(:),Umatrix(:)
    Complex(pr) prefac,detA(2)
    !
    If(beta_active(ibet).Eq.0) Return
    ! ------
    ! Timing
    ! ------
    Call system_clock(tstart_ibet,clock_rate)
    ! ---------------------------------------- !
    !                AMP only                  !
    ! ---------------------------------------- !
    If(AMP_is_on.Eq.1 .And. PNP_is_on.Eq.0) Then
       !--------------------------------
       ! Allocating local overlap arrays
       !--------------------------------
       Allocate(rotated_Vmatrix(4,maxphi_eff),rotated_Umatrix(4,maxphi_eff),InverseOverlapMatrix(4,maxphi_eff,1),Vmatrix(4),Umatrix(4))
       Do iit=1,4
          Allocate(Vmatrix(iit)%arr(1:nt,1:kdim(iit)),Umatrix(iit)%arr(1:nt,1:kdim(iit)))
          Do iphi=1,maxphi_eff
             Allocate(rotated_Vmatrix(iit,iphi)%arr(1:nt,1:kdim(iit))); Allocate(rotated_Umatrix(iit,iphi)%arr(1:nt,1:kdim(iit)))
             isimplex=1
             Allocate(InverseOverlapMatrix(iit,iphi,isimplex)%arr(1:kdim(iit),1:kdim(iit)))
             InverseOverlapMatrix(iit,iphi,isimplex)%arr(1:kdim(iit),1:kdim(iit))=czero
             Do k1=1,kdim(iit)
                InverseOverlapMatrix(iit,iphi,isimplex)%arr(k1,k1)=cone
             End Do
          End Do
       End Do
       !
       Vmatrix(1)%arr(1:nt,1:kdim(1))=VmatrixN1(1:nt,1:kdim(1))
       Vmatrix(2)%arr(1:nt,1:kdim(2))=VmatrixP1(1:nt,1:kdim(2))
       Vmatrix(3)%arr(1:nt,1:kdim(3))=VmatrixN2(1:nt,1:kdim(3))
       Vmatrix(4)%arr(1:nt,1:kdim(4))=VmatrixP2(1:nt,1:kdim(4))
       Umatrix(1)%arr(1:nt,1:kdim(1))=UmatrixN1(1:nt,1:kdim(1))
       Umatrix(2)%arr(1:nt,1:kdim(2))=UmatrixP1(1:nt,1:kdim(2))
       Umatrix(3)%arr(1:nt,1:kdim(3))=UmatrixN2(1:nt,1:kdim(3))
       Umatrix(4)%arr(1:nt,1:kdim(4))=UmatrixP2(1:nt,1:kdim(4))
       !
       !---------------------
       ! Calculating overlaps
       !---------------------
!$OMP  PARALLEL DO       &
!$OMP& DEFAULT(NONE)     &
!$OMP& SCHEDULE(DYNAMIC) &
!$OMP& SHARED(itiphi_pair1,itiphi_pair2,nt,kdim,ibet,rotation_matrix,inverse_rotation_matrix,Vmatrix,Umatrix,rotated_Vmatrix,rotated_Umatrix,&
!$OMP&        InverseOverlapMatrix,detR,phiabs,ephi,ephic,lproj,maxphi_eff,maxphi,all_overlaps,team_rank,team_color,beta_size) &
!$OMP& PRIVATE(iit,iphi,ifl,detA,OverlapMatrixTmp,aUmatrix,aVmatrix,isimplex,prefac)
       Do iit_iphi=1,4*maxphi_eff
          iit=itiphi_pair1(iit_iphi); iphi=itiphi_pair2(iit_iphi)
          Allocate(aUmatrix(kdim(iit),kdim(iit)),aVmatrix(kdim(iit),kdim(iit))); aUmatrix=czero; aVmatrix=czero
          ! Rotating V matrix
          Call Zgemm('N','N',nt,kdim(iit),nt,cone,rotation_matrix,nt,Vmatrix(iit)%arr(1:nt,1:kdim(iit)),nt,czero,rotated_Vmatrix(iit,iphi)%arr(1:nt,1:kdim(iit)),nt)                
          Call Zgemm('T','N',kdim(iit),kdim(iit),nt,cone,Vmatrix(iit)%arr(1:nt,1:kdim(iit)),nt,conjg(rotated_Vmatrix(iit,iphi)%arr(1:nt,1:kdim(iit))),nt,czero,aVmatrix(1:kdim(iit),1:kdim(iit)),kdim(iit))
          ! Rotating U matrix
          Call Zgemm('T','N',nt,kdim(iit),nt,cone,inverse_rotation_matrix,nt,conjg(Umatrix(iit)%arr(1:nt,1:kdim(iit))),nt,czero,rotated_Umatrix(iit,iphi)%arr(1:nt,1:kdim(iit)),nt)    
          Call Zgemm('T','N',kdim(iit),kdim(iit),nt,cone,Umatrix(iit)%arr(1:nt,1:kdim(iit)),nt,rotated_Umatrix(iit,iphi)%arr(1:nt,1:kdim(iit)),nt,cone,aUmatrix(1:kdim(iit),1:kdim(iit)),kdim(iit))
          !
          Allocate(OverlapMatrixTmp(kdim(iit),kdim(iit))); OverlapMatrixTmp=czero
          isimplex=1
          OverlapMatrixTmp(1:kdim(iit),1:kdim(iit))=ephic(iphi)*aUmatrix(1:kdim(iit),1:kdim(iit))+ephi(iphi)*aVmatrix(1:kdim(iit),1:kdim(iit))
          ! Calculating overlap          
          Call calculate_inverse(kdim(iit),OverlapMatrixTmp,detA(isimplex),ifl)
          InverseOverlapMatrix(iit,iphi,isimplex)%arr(1:kdim(iit),1:kdim(iit))=OverlapMatrixTmp(1:kdim(iit),1:kdim(iit))
          Deallocate(OverlapMatrixTmp)
          If(iit.Eq. 3 .Or. iit .Eq. 4) Then
             prefac=Cmplx(Cos(phiabs(iphi)*kdim(iit)),Sin(phiabs(iphi)*kdim(iit)))
             all_overlaps((iit-3)*maxphi*beta_size+(iphi-1)*beta_size+ibet-team_rank*beta_size)=detA(1)*detR(ibet)*prefac ! saving overlap
          End If
          Deallocate(aUmatrix,aVmatrix)
       End Do
!$OMP  End PARALLEL DO
       !-----------------------------------------------------
       ! Calculating rotated densities in configuration space
       !-----------------------------------------------------
!$OMP  PARALLEL DO       &
!$OMP& DEFAULT(NONE)     &
!$OMP& SCHEDULE(DYNAMIC) &
!$OMP& SHARED(nt,kdim,rotated_Vmatrix,InverseOverlapMatrix,Vmatrix,rotated_density,Umatrix,rotated_kappa,itiphi_pair1,itiphi_pair2,ephi,maxphi_eff,maxphi) &
!$OMP& PRIVATE(it,iit_iphi,iphi,AuxiliaryMatrix,AuxiliaryDensity)
       Do iit_iphi=1,2*maxphi_eff
          it=itiphi_pair1(iit_iphi); iphi=itiphi_pair2(iit_iphi)
          Allocate(AuxiliaryMatrix(nt,kdim(it)),AuxiliaryDensity(nt,nt))
          ! rho density
          ! simplex ++
          Call Zgemm('N','N',nt,kdim(it),kdim(it),cone,rotated_Vmatrix(it,iphi)%arr(1:nt,1:kdim(it)),nt,conjg(InverseOverlapMatrix(it,iphi,1)%arr(1:kdim(it),1:kdim(it))),kdim(it),czero,AuxiliaryMatrix,nt)
          Call Zgemm('N','C',nt,nt,kdim(it),cone,AuxiliaryMatrix,nt,Vmatrix(it)%arr(1:nt,1:kdim(it)),nt,czero,AuxiliaryDensity,nt)
          rotated_density(it,iphi,1)%arr(1:nt,1:nt)=ephi(iphi)*AuxiliaryDensity(1:nt,1:nt)
          ! simplex --
          rotated_density(it,iphi,2)%arr(1:nt,1:nt)=conjg(rotated_density(it,iphi,1)%arr(1:nt,1:nt))
          ! kappa tensor
          ! simplex +-
          Call Zgemm('N','N',nt,kdim(it),kdim(it),cone,rotated_Vmatrix(it,iphi)%arr(1:nt,1:kdim(it)),nt,conjg(InverseOverlapMatrix(it,iphi,1)%arr(1:kdim(it),1:kdim(it))),kdim(it),czero,AuxiliaryMatrix,nt)
          Call Zgemm('N','C',nt,nt,kdim(it),-cone,AuxiliaryMatrix,nt,Umatrix(it)%arr(1:nt,1:kdim(it)),nt,czero,AuxiliaryDensity,nt)
          rotated_kappa(it,iphi,1)%arr(1:nt,1:nt)=ephi(iphi)*AuxiliaryDensity(1:nt,1:nt)
          ! simplex -+
          rotated_kappa(it,iphi,2)%arr(1:nt,1:nt)=-conjg(rotated_kappa(it,iphi,1)%arr(1:nt,1:nt))
          Deallocate(AuxiliaryMatrix,AuxiliaryDensity)
       End Do
!$OMP  End PARALLEL DO
       !-----------------------------------------------------
       ! Determining non-zero density and kappa contributions
       !-----------------------------------------------------
       Do iphi=1,maxphi_eff
          iosc_contributing=0
          Do iosc2=1,nt
          Do iosc1=1,nt
             If(Abs(rotated_density(1,iphi,1)%arr(iosc1,iosc2)).Gt.eps .Or. Abs(rotated_density(2,iphi,1)%arr(iosc1,iosc2)).Gt.eps .Or. &
                Abs(rotated_kappa(1,iphi,1)%arr(iosc1,iosc2)).Gt.eps  .Or. Abs(rotated_kappa(2,iphi,1)%arr(iosc1,iosc2)).Gt.eps) Then
                iosc_contributing=iosc_contributing+1
                iosc1_contributing(iosc_contributing,iphi)=iosc1
                iosc2_contributing(iosc_contributing,iphi)=iosc2
             End If
          End Do
          End Do
          nt_contributing(iphi)=iosc_contributing
       End Do
       !
       Deallocate(rotated_Vmatrix,rotated_Umatrix,InverseOverlapMatrix,Umatrix,Vmatrix)
       !
    ! ---------------------------------------- !
    !             PNP or PNP&AMP               !
    ! ---------------------------------------- !
    Else
       !--------------------------------
       ! Allocating local overlap arrays
       !--------------------------------
       Allocate(rotated_Vmatrix(4,maxphi_eff),rotated_Umatrix(4,maxphi_eff),InverseOverlapMatrix(4,maxphi_eff,2),Vmatrix(4),Umatrix(4))
       Do iit=1,4
          Allocate(Vmatrix(iit)%arr(1:nt,1:kdim(iit)),Umatrix(iit)%arr(1:nt,1:kdim(iit)))
          Do iphi=1,maxphi_eff
             Allocate(rotated_Vmatrix(iit,iphi)%arr(1:nt,1:kdim(iit))); Allocate(rotated_Umatrix(iit,iphi)%arr(1:nt,1:kdim(iit)))
             Do isimplex=1,2
                Allocate(InverseOverlapMatrix(iit,iphi,isimplex)%arr(1:kdim(iit),1:kdim(iit)))
                InverseOverlapMatrix(iit,iphi,isimplex)%arr(1:kdim(iit),1:kdim(iit))=czero
                Do k1=1,kdim(iit)
                   InverseOverlapMatrix(iit,iphi,isimplex)%arr(k1,k1)=cone
                End Do
             End Do
          End Do
       End Do
       Vmatrix(1)%arr(1:nt,1:kdim(1))=VmatrixN1(1:nt,1:kdim(1))
       Vmatrix(2)%arr(1:nt,1:kdim(2))=VmatrixP1(1:nt,1:kdim(2))
       Vmatrix(3)%arr(1:nt,1:kdim(3))=VmatrixN2(1:nt,1:kdim(3))
       Vmatrix(4)%arr(1:nt,1:kdim(4))=VmatrixP2(1:nt,1:kdim(4))
       Umatrix(1)%arr(1:nt,1:kdim(1))=UmatrixN1(1:nt,1:kdim(1))
       Umatrix(2)%arr(1:nt,1:kdim(2))=UmatrixP1(1:nt,1:kdim(2))
       Umatrix(3)%arr(1:nt,1:kdim(3))=UmatrixN2(1:nt,1:kdim(3))
       Umatrix(4)%arr(1:nt,1:kdim(4))=UmatrixP2(1:nt,1:kdim(4))
       !---------------------
       ! Calculating overlaps
       !---------------------
!$OMP  PARALLEL DO       &
!$OMP& DEFAULT(NONE)     &
!$OMP& SCHEDULE(DYNAMIC) &
!$OMP& SHARED(itiphi_pair1,itiphi_pair2,nt,kdim,ibet,rotation_matrix,inverse_rotation_matrix,Vmatrix,Umatrix,rotated_Vmatrix,rotated_Umatrix,&
!$OMP&        InverseOverlapMatrix,detR,phiabs,ephi,ephic,lproj,maxphi_eff,maxphi,all_overlaps,beta_size,team_rank) &
!$OMP& PRIVATE(iit,iphi,ifl,detA,OverlapMatrixTmp,aUmatrix,aVmatrix,isimplex,prefac)
       Do iit_iphi=1,4*maxphi_eff
          iit=itiphi_pair1(iit_iphi); iphi=itiphi_pair2(iit_iphi)
          Allocate(aUmatrix(kdim(iit),kdim(iit)),aVmatrix(kdim(iit),kdim(iit))); aUmatrix=czero; aVmatrix=czero
          ! Rotating V matrix
          Call Zgemm('N','N',nt,kdim(iit),nt,cone,rotation_matrix,nt,Vmatrix(iit)%arr(1:nt,1:kdim(iit)),nt,czero,rotated_Vmatrix(iit,iphi)%arr(1:nt,1:kdim(iit)),nt)               
          Call Zgemm('T','N',kdim(iit),kdim(iit),nt,cone,Vmatrix(iit)%arr(1:nt,1:kdim(iit)),nt,conjg(rotated_Vmatrix(iit,iphi)%arr(1:nt,1:kdim(iit))),nt,czero,aVmatrix(1:kdim(iit),1:kdim(iit)),kdim(iit))
          ! Rotating U matrix
          Call Zgemm('T','N',nt,kdim(iit),nt,cone,inverse_rotation_matrix,nt,conjg(Umatrix(iit)%arr(1:nt,1:kdim(iit))),nt,czero,rotated_Umatrix(iit,iphi)%arr(1:nt,1:kdim(iit)),nt) 
          Call Zgemm('T','N',kdim(iit),kdim(iit),nt,cone,Umatrix(iit)%arr(1:nt,1:kdim(iit)),nt,rotated_Umatrix(iit,iphi)%arr(1:nt,1:kdim(iit)),nt,cone,aUmatrix(1:kdim(iit),1:kdim(iit)),kdim(iit))
          ! Loop over two simplex blocks
          Allocate(OverlapMatrixTmp(kdim(iit),kdim(iit)))
          Do isimplex=1,2
             OverlapMatrixTmp=czero
             If(isimplex .Eq. 1) Then
                OverlapMatrixTmp(1:kdim(iit),1:kdim(iit))=ephic(iphi)*aUmatrix(1:kdim(iit),1:kdim(iit))+ephi(iphi)*aVmatrix(1:kdim(iit),1:kdim(iit))
             Else
                OverlapMatrixTmp(1:kdim(iit),1:kdim(iit))=ephic(iphi)*conjg(aUmatrix(1:kdim(iit),1:kdim(iit)))+ephi(iphi)*conjg(aVmatrix(1:kdim(iit),1:kdim(iit)))
             End If
             ! Calculating overlap          
             Call calculate_inverse(kdim(iit),OverlapMatrixTmp,detA(isimplex),ifl)
             InverseOverlapMatrix(iit,iphi,isimplex)%arr(1:kdim(iit),1:kdim(iit))=OverlapMatrixTmp(1:kdim(iit),1:kdim(iit))   
          End Do
          Deallocate(OverlapMatrixTmp)
          If(iit.Eq. 3 .Or. iit .Eq. 4) Then
             prefac=Cmplx(Cos(phiabs(iphi)*kdim(iit)),Sin(phiabs(iphi)*kdim(iit)))
             all_overlaps((iit-3)*maxphi*beta_size+(iphi-1)*beta_size+ibet-team_rank*beta_size)=detA(1)*detR(ibet)*prefac ! saving overlap
             If(iphi.Gt.1) Then
                all_overlaps((iit-3)*maxphi*beta_size+(maxphi-iphi+2-1)*beta_size+ibet-team_rank*beta_size)=conjg(detA(1)*detR(ibet)*prefac)
             End If
          End If
          Deallocate(aUmatrix,aVmatrix)
       End Do
!$OMP End PARALLEL DO
       !-----------------------------------------------------
       ! Calculating rotated densities in configuration space
       !-----------------------------------------------------
!$OMP  PARALLEL DO       &
!$OMP& DEFAULT(NONE)     &
!$OMP& SCHEDULE(DYNAMIC) &
!$OMP& SHARED(nt,kdim,rotated_Vmatrix,InverseOverlapMatrix,Vmatrix,rotated_density,Umatrix,rotated_kappa,itiphi_pair1,itiphi_pair2,ephi,maxphi_eff,maxphi) &
!$OMP& PRIVATE(it,iit_iphi,iphi,AuxiliaryMatrix,AuxiliaryDensity)
       Do iit_iphi=1,2*maxphi_eff
          it=itiphi_pair1(iit_iphi); iphi=itiphi_pair2(iit_iphi)
          Allocate(AuxiliaryMatrix(nt,kdim(it)),AuxiliaryDensity(nt,nt))
          ! rho density
          ! simplex ++
          Call Zgemm('N','N',nt,kdim(it),kdim(it),cone,rotated_Vmatrix(it,iphi)%arr(1:nt,1:kdim(it)),nt,InverseOverlapMatrix(it,iphi,2)%arr(1:kdim(it),1:kdim(it)),kdim(it),czero,AuxiliaryMatrix,nt)
          Call Zgemm('N','C',nt,nt,kdim(it),cone,AuxiliaryMatrix,nt,Vmatrix(it)%arr(1:nt,1:kdim(it)),nt,czero,AuxiliaryDensity,nt)
          rotated_density(it,iphi,1)%arr(1:nt,1:nt)=ephi(iphi)*AuxiliaryDensity(1:nt,1:nt)
          ! simplex --
          Call Zgemm('N','N',nt,kdim(it),kdim(it),cone,conjg(rotated_Vmatrix(it,iphi)%arr(1:nt,1:kdim(it))),nt,InverseOverlapMatrix(it,iphi,1)%arr(1:kdim(it),1:kdim(it)),kdim(it),czero,AuxiliaryMatrix,nt)
          Call Zgemm('N','T',nt,nt,kdim(it),cone,AuxiliaryMatrix,nt,Vmatrix(it)%arr(1:nt,1:kdim(it)),nt,czero,AuxiliaryDensity,nt)
          rotated_density(it,iphi,2)%arr(1:nt,1:nt)=ephi(iphi)*AuxiliaryDensity(1:nt,1:nt)
          ! kappa tensor
          ! simplex +-
          Call Zgemm('N','N',nt,kdim(it),kdim(it),cone,rotated_Vmatrix(it,iphi)%arr(1:nt,1:kdim(it)),nt,InverseOverlapMatrix(it,iphi,2)%arr(1:kdim(it),1:kdim(it)),kdim(it),czero,AuxiliaryMatrix,nt)
          Call Zgemm('N','C',nt,nt,kdim(it),-cone,AuxiliaryMatrix,nt,Umatrix(it)%arr(1:nt,1:kdim(it)),nt,czero,AuxiliaryDensity,nt)
          rotated_kappa(it,iphi,1)%arr(1:nt,1:nt)=ephi(iphi)*AuxiliaryDensity(1:nt,1:nt)
          ! simplex -+
          Call Zgemm('N','N',nt,kdim(it),kdim(it),cone,conjg(rotated_Vmatrix(it,iphi)%arr(1:nt,1:kdim(it))),nt,InverseOverlapMatrix(it,iphi,1)%arr(1:kdim(it),1:kdim(it)),kdim(it),czero,AuxiliaryMatrix,nt)
          Call Zgemm('N','T',nt,nt,kdim(it),cone,AuxiliaryMatrix,nt,Umatrix(it)%arr(1:nt,1:kdim(it)),nt,czero,AuxiliaryDensity,nt)
          rotated_kappa(it,iphi,2)%arr(1:nt,1:nt)=ephi(iphi)*AuxiliaryDensity(1:nt,1:nt)
          Deallocate(AuxiliaryMatrix,AuxiliaryDensity)
       End Do
!$OMP End PARALLEL DO
       !-----------------------------------------------------
       ! Determining non-zero density and kappa contributions
       !-----------------------------------------------------
       Do iphi=1,maxphi_eff
          iosc_contributing=0
          Do iosc2=1,nt
          Do iosc1=1,nt
             If(Abs(rotated_density(1,iphi,1)%arr(iosc1,iosc2)).Gt.eps .Or. Abs(rotated_density(1,iphi,2)%arr(iosc1,iosc2)).Gt.eps .Or. &
                Abs(rotated_density(2,iphi,1)%arr(iosc1,iosc2)).Gt.eps .Or. Abs(rotated_density(2,iphi,2)%arr(iosc1,iosc2)).Gt.eps .Or. &
                Abs(rotated_kappa(1,iphi,1)%arr(iosc1,iosc2)).Gt.eps .Or. Abs(rotated_kappa(1,iphi,2)%arr(iosc1,iosc2)).Gt.eps .Or. &
                Abs(rotated_kappa(2,iphi,1)%arr(iosc1,iosc2)).Gt.eps .Or. Abs(rotated_kappa(2,iphi,2)%arr(iosc1,iosc2)).Gt.eps) Then
                iosc_contributing=iosc_contributing+1
                iosc1_contributing(iosc_contributing,iphi)=iosc1
                iosc2_contributing(iosc_contributing,iphi)=iosc2
             End If
          End Do
          End Do
          nt_contributing(iphi)=iosc_contributing
       End Do
       !
       Deallocate(rotated_Vmatrix,rotated_Umatrix,InverseOverlapMatrix)
    End If
    !
    !-------
    ! Timing
    !-------
    Call system_clock(tfinish_ibet,clock_rate)
    time_ibet(3)=time_ibet(3)+Real(tfinish_ibet-tstart_ibet)/Real(clock_rate)
    ! ----------------------------
    ! ZGEMM subroutine explanation
    ! ----------------------------
    ! >  SUBROUTINE ZGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
    ! >  performs : C := alpha*op( A )*op( B ) + beta*C
    ! >  arguments:
    ! >      - TRANSA and TRANSB define operations: 'N', 'T', 'C'
    ! >      - M is integer defining number of rows of matrix op(A) and of matrix C
    ! >      - N is integer defining number of colums of matrix op(B) and of matrix C
    ! >      - K is integer defining number of columns of matrix op(A) and number of rows of matrix op(B)
    ! >      - LDA and LDB leading dimensions (number of rows in definition) of matrices A and B
    ! >      - LDC specifies the first dimension of C as declared in  the  calling  (sub)  program.
    ! >      - alpha and beta are complex numbers
    ! >      - 'T' transposes the matrix, 'C' acts as a dagger (transpose and complex conjugate).
  End Subroutine calculate_overlaps
  !==============================================
  ! Subroutine that calculates rotated densities.
  !==============================================
  Subroutine calculate_densities(ibet)
   Implicit None
   Integer(ipr) ibet,ihil_iphicyl,ihil,iphicyl,iosc,iosc1,iosc2,it,nla1,nla2,iphi,iphin,iphip,ihil1,ihil2
   Integer(ipr) tstart_ibet,tfinish_ibet
   Complex(pr) ro_sumN,ro_sumP,tau_sumN,tau_sumP,dj_sumN,dj_sumP,dro_sumN,dro_sumP,sfir_sumN,sfir_sumP,s,sfiz_sumN,sfiz_sumP,srfi_sumN,srfi_sumP,szfi_sumN,&
            szfi_sumP,szz_sumN,szz_sumP,srz_sumN,srz_sumP,srr_sumN,srr_sumP,szr_sumN,szr_sumP,sfifi_sumN,sfifi_sumP,&
            aka_sumN,aka_sumP,fac_ro,fac_tau,fac_dj1,fac_dj2,fac_sfiz,fac_sfir,fac_srfi,fac_szfi,fac_szz,fac_srz,fac_srr,fac_szr,&
            fac_sfifi,fac_dro,fac_aka,ro_projN,ro_projP,ro_normN,ro_normP,ro_norm
   Real(pr) y,y2,phy
   Real(pr) q_h0l0a_h0l0b,q_h1l0a_h1l0b,q_h0l1a_h0l1b,q_h0l0a_h0l1b,q_h0l1a_h0l0b,q_h0l0a_h1l0b,q_h1l0a_h0l0b,q_h1l0a_h0l1b,q_h0l1a_h1l0b,&
            rhoN_real,rhoP_real,rhoN_imag,rhoP_imag,akaN_real,akaP_real
   Real(pr) bpi,bpi2,bzi,bzi2
   Complex(pr) densNpl,densPpl,densNmi,densPmi,kappNmi,kappPmi
   Complex(pr) rotated_ro_N,rotated_ro_P
   Complex(pr), Allocatable :: vc_pr(:,:)
   !
   If(beta_active(ibet).Eq.0) Return
   !-------
   ! Timing
   !-------
   Call system_clock(tstart_ibet,clock_rate)
   !
   bpi=one/bp; bpi2=bpi**2; bzi=one/bz; bzi2=bzi**2
   !
!$OMP  PARALLEL DO       &
!$OMP& DEFAULT(NONE)     &
!$OMP& SCHEDULE(DYNAMIC) &
!$OMP& SHARED(bpi2,bzi2,nghl,ihil_convert,iphicyl_convert,phicyl,y_opt,nt,nl_sim,rotated_density,rotated_kappa,&
!$OMP&       qhla_opt,fi1z_opt,fi1r_opt,rotated_ro,rotated_tau,rotated_dj,rotated_dro,rotated_sfiz,rotated_sfir,&
!$OMP&       rotated_srfi,rotated_szfi,rotated_szz,rotated_srz,rotated_srr,rotated_szr,rotated_sfifi,rotated_aka,xl_ihil,xh_ihil,nr,nz,&
!$OMP&       nt_contributing,iosc1_contributing,iosc2_contributing,maxphi,maxphi_eff,ibet,lproj,ro_projected,ephicN,ephicP,rotated_overlap,team_color) &
!$OMP& PRIVATE(ihil_iphicyl,ihil,iphicyl,phy,y,y2,ro_sumN,ro_sumP,tau_sumN,tau_sumP,dj_sumN,dj_sumP,dro_sumN,dro_sumP,sfir_sumN,sfir_sumP,&
!$OMP&         sfiz_sumN,sfiz_sumP,srfi_sumN,srfi_sumP,szfi_sumN,szfi_sumP,szz_sumN,szz_sumP,srz_sumN,srz_sumP,srr_sumN,srr_sumP,szr_sumN,&
!$OMP&         szr_sumP,sfifi_sumN,sfifi_sumP,aka_sumN,aka_sumP,iosc,iosc1,iosc2,nla1,nla2,rhoN_real,rhoP_real,rhoN_imag,rhoP_imag,akaN_real,&
!$OMP&         akaP_real,q_h0l0a_h0l0b,q_h1l0a_h1l0b,q_h0l1a_h0l1b,q_h0l0a_h0l1b,q_h0l1a_h0l0b,q_h0l0a_h1l0b,q_h1l0a_h0l0b,q_h1l0a_h0l1b,&
!$OMP&         q_h0l1a_h1l0b,fac_ro,fac_tau,fac_dj1,fac_dj2,fac_dro,fac_sfiz,fac_sfir,fac_srfi,fac_szfi,fac_szz,fac_srz,fac_srr,fac_szr,&
!$OMP&         fac_sfifi,fac_aka,densNpl,densPpl,densNmi,densPmi,kappNmi,kappPmi,iphi,ro_projN,ro_projP,ro_normN,ro_normP)
   Do ihil_iphicyl=1,nghl*ngphi
      ihil=ihil_convert(ihil_iphicyl); iphicyl=iphicyl_convert(ihil_iphicyl)
      !
      phy=phicyl(iphicyl)
      y=y_opt(ihil);y2=y*y
      !
      Do iphi=1,maxphi_eff
         ro_sumN=czero;ro_sumP=czero;tau_sumN=czero;tau_sumP=czero
         sfiz_sumN=czero;sfiz_sumP=czero;sfir_sumN=czero;sfir_sumP=czero;srfi_sumN=czero;srfi_sumP=czero;szfi_sumN=czero;szfi_sumP=czero
         szz_sumN=czero;szz_sumP=czero;srz_sumN=czero;srz_sumP=czero;srr_sumN=czero;srr_sumP=czero;szr_sumN=czero;szr_sumP=czero
         sfifi_sumN=czero;sfifi_sumP=czero;dro_sumN=czero;dro_sumP=czero;dj_sumN=czero;dj_sumP=czero;aka_sumN=czero;aka_sumP=czero
         Do iosc=1,nt_contributing(iphi)
            iosc1=iosc1_contributing(iosc,iphi); iosc2=iosc2_contributing(iosc,iphi)
            nla1=nl_sim(iosc1); nla2=nl_sim(iosc2)
            !
            densNpl=rotated_density(1,iphi,1)%arr(iosc1,iosc2)+rotated_density(1,iphi,2)%arr(iosc1,iosc2)
            densNmi=rotated_density(1,iphi,1)%arr(iosc1,iosc2)-rotated_density(1,iphi,2)%arr(iosc1,iosc2)
            densPpl=rotated_density(2,iphi,1)%arr(iosc1,iosc2)+rotated_density(2,iphi,2)%arr(iosc1,iosc2)
            densPmi=rotated_density(2,iphi,1)%arr(iosc1,iosc2)-rotated_density(2,iphi,2)%arr(iosc1,iosc2)
            kappNmi=rotated_kappa(1,iphi,1)%arr(iosc1,iosc2)-rotated_kappa(1,iphi,2)%arr(iosc1,iosc2)
            kappPmi=rotated_kappa(2,iphi,1)%arr(iosc1,iosc2)-rotated_kappa(2,iphi,2)%arr(iosc1,iosc2)
            !
            q_h0l0a_h0l0b=qhla_opt(iosc1,ihil)*qhla_opt(iosc2,ihil)
            q_h1l0a_h1l0b=fi1z_opt(iosc1,ihil)*fi1z_opt(iosc2,ihil); q_h0l1a_h0l1b=fi1r_opt(iosc1,ihil)*fi1r_opt(iosc2,ihil)
            q_h0l0a_h0l1b=qhla_opt(iosc1,ihil)*fi1r_opt(iosc2,ihil); q_h0l1a_h0l0b=fi1r_opt(iosc1,ihil)*qhla_opt(iosc2,ihil)
            q_h0l0a_h1l0b=qhla_opt(iosc1,ihil)*fi1z_opt(iosc2,ihil); q_h1l0a_h0l0b=fi1z_opt(iosc1,ihil)*qhla_opt(iosc2,ihil)
            q_h1l0a_h0l1b=fi1z_opt(iosc1,ihil)*fi1r_opt(iosc2,ihil); q_h0l1a_h1l0b=fi1r_opt(iosc1,ihil)*fi1z_opt(iosc2,ihil)
            !
            ! ro density
            fac_ro=q_h0l0a_h0l0b*Cos((nla2-nla1)*phy)
            ro_sumN=ro_sumN+densNpl*fac_ro; ro_sumP=ro_sumP+densPpl*fac_ro
            ! tau density
            fac_tau=(q_h1l0a_h1l0b+q_h0l1a_h0l1b+nla1*nla2*y2*q_h0l0a_h0l0b)*Cos((nla2-nla1)*phy)
            tau_sumN=tau_sumN+densNpl*fac_tau; tau_sumP=tau_sumP+densPpl*fac_tau
            ! dj density
            fac_dj1=half*y*(nla1*q_h0l0a_h0l1b+nla2*q_h0l1a_h0l0b)*Cos((nla2-nla1)*phy)
            fac_dj2=half*(y*(nla2*q_h1l0a_h0l0b-nla1*q_h0l0a_h1l0b)-q_h1l0a_h0l1b+q_h0l1a_h1l0b)*Cos((nla1+nla2+1)*phy)*cimag
            !
            dj_sumN=dj_sumN+densNpl*fac_dj1+densNmi*fac_dj2; dj_sumP=dj_sumP+densPpl*fac_dj1+densPmi*fac_dj2
            ! dro density
            fac_dro=two*(xl_ihil(ihil)*bpi2+(xh_ihil(ihil))**2*bzi2-two*bpi2*(nr(iosc1)+nr(iosc2)+0.5d0*(abs(nla1)+abs(nla2))+1)&
                    -one*bzi2*(nz(iosc1)+nz(iosc2)+1)+nla1*nla2*y2)*q_h0l0a_h0l0b+two*(q_h0l1a_h0l1b+q_h1l0a_h1l0b)
            fac_dro=fac_dro*half*Cos((nla2-nla1)*phy)
            dro_sumN=dro_sumN+densNpl*fac_dro; dro_sumP=dro_sumP+densPpl*fac_dro
            ! tensor components
            ! sfiz
            fac_sfiz=half*q_h0l0a_h0l0b*y*(nla1+nla2)*Cos((nla2-nla1)*phy)
            sfiz_sumN=sfiz_sumN+densNpl*fac_sfiz; sfiz_sumP=sfiz_sumP+densPpl*fac_sfiz
            ! sfir
            fac_sfir=half*cimag*q_h0l0a_h0l0b*y*(nla1-nla2)*Cos((nla2+nla1+1)*phy)
            sfir_sumN=sfir_sumN+densNmi*fac_sfir; sfir_sumP=sfir_sumP+densPmi*fac_sfir
            ! srfi
            fac_srfi=half*cimag*(q_h0l1a_h0l0b-q_h0l0a_h0l1b)*Cos((nla1+nla2+1)*phy)
            srfi_sumN=srfi_sumN+densNmi*fac_srfi; srfi_sumP=srfi_sumP+densPmi*fac_srfi
            ! szfi
            fac_szfi=half*cimag*(q_h1l0a_h0l0b-q_h0l0a_h1l0b)*Cos((nla1+nla2+1)*phy)
            szfi_sumN=szfi_sumN+densNmi*fac_szfi; szfi_sumP=szfi_sumP+densPmi*fac_szfi
            ! szz
            fac_szz=half*(q_h0l0a_h1l0b-q_h1l0a_h0l0b)*Sin((nla2-nla1)*phy)
            szz_sumN=szz_sumN+densNpl*fac_szz; szz_sumP=szz_sumP+densPpl*fac_szz
            ! srz
            fac_srz=half*(q_h0l0a_h0l1b-q_h0l1a_h0l0b)*Sin((nla2-nla1)*phy)
            srz_sumN=srz_sumN+densNpl*fac_srz; srz_sumP=srz_sumP+densPpl*fac_srz
            ! srr
            fac_srr=half*cimag*(q_h0l1a_h0l0b-q_h0l0a_h0l1b)*Sin((nla1+nla2+1)*phy)
            srr_sumN=srr_sumN+densNmi*fac_srr; srr_sumP=srr_sumP+densPmi*fac_srr
            ! szr
            fac_szr=half*cimag*(q_h1l0a_h0l0b-q_h0l0a_h1l0b)*Sin((nla1+nla2+1)*phy)
            szr_sumN=szr_sumN+densNmi*fac_szr; szr_sumP=szr_sumP+densPmi*fac_szr
            ! sfifi
            fac_sfifi=half*cimag*q_h0l0a_h0l0b*y*(nla1-nla2)*Sin((nla1+nla2+1)*phy)
            sfifi_sumN=sfifi_sumN+densNmi*fac_sfifi; sfifi_sumP=sfifi_sumP+densPmi*fac_sfifi
            ! aka density
            fac_aka=half*q_h0l0a_h0l0b*Cos((nla2-nla1)*phy)
            aka_sumN=aka_sumN+kappNmi*fac_aka; aka_sumP=aka_sumP+kappPmi*fac_aka
            !
         End Do ! iosc
         !
         rotated_ro(iphi,ihil_iphicyl,1)=ro_sumN;       rotated_ro(iphi,ihil_iphicyl,2)=ro_sumP
         rotated_tau(iphi,ihil_iphicyl,1)=tau_sumN;     rotated_tau(iphi,ihil_iphicyl,2)=tau_sumP
         rotated_dj(iphi,ihil_iphicyl,1)=dj_sumN;       rotated_dj(iphi,ihil_iphicyl,2)=dj_sumP
         rotated_dro(iphi,ihil_iphicyl,1)=dro_sumN;     rotated_dro(iphi,ihil_iphicyl,2)=dro_sumP
         rotated_sfiz(iphi,ihil_iphicyl,1)=sfiz_sumN;   rotated_sfiz(iphi,ihil_iphicyl,2)=sfiz_sumP
         rotated_sfir(iphi,ihil_iphicyl,1)=sfir_sumN;   rotated_sfir(iphi,ihil_iphicyl,2)=sfir_sumP
         rotated_srfi(iphi,ihil_iphicyl,1)=srfi_sumN;   rotated_srfi(iphi,ihil_iphicyl,2)=srfi_sumP
         rotated_szfi(iphi,ihil_iphicyl,1)=szfi_sumN;   rotated_szfi(iphi,ihil_iphicyl,2)=szfi_sumP
         rotated_szz(iphi,ihil_iphicyl,1)=szz_sumN;     rotated_szz(iphi,ihil_iphicyl,2)=szz_sumP
         rotated_srz(iphi,ihil_iphicyl,1)=srz_sumN;     rotated_srz(iphi,ihil_iphicyl,2)=srz_sumP
         rotated_srr(iphi,ihil_iphicyl,1)=srr_sumN;     rotated_srr(iphi,ihil_iphicyl,2)=srr_sumP
         rotated_szr(iphi,ihil_iphicyl,1)=szr_sumN;     rotated_szr(iphi,ihil_iphicyl,2)=szr_sumP
         rotated_sfifi(iphi,ihil_iphicyl,1)=sfifi_sumN; rotated_sfifi(iphi,ihil_iphicyl,2)=sfifi_sumP
         rotated_aka(iphi,ihil_iphicyl,1)=aka_sumN;     rotated_aka(iphi,ihil_iphicyl,2)=aka_sumP
        !
      End Do ! iphi
      !
   End Do ! ihil_iphicyl
!$OMP End PARALLEL DO
   !
!$OMP  PARALLEL DO       &
!$OMP& DEFAULT(NONE)     &
!$OMP& SCHEDULE(DYNAMIC) &
!$OMP& SHARED(rotated_ro,nghl,maxphi,maxphi_eff,ibet,ephicN,ephicP,rotated_overlap,ro_projected) &
!$OMP& PRIVATE(ihil_iphicyl,ro_projN,ro_projP,ro_norm,iphin,iphip,rotated_ro_N,rotated_ro_P)
   Do ihil_iphicyl=1,nghl*ngphi
      ro_projN=czero;ro_projP=czero;ro_norm=czero
      Do iphin=1,maxphi
      Do iphip=1,maxphi
         If(iphin .Le. maxphi_eff) Then
            rotated_ro_N=rotated_ro(iphin,ihil_iphicyl,1)
         Else
            rotated_ro_N=conjg(rotated_ro(maxphi-iphin+2,ihil_iphicyl,1))
         End If
         !
         If(iphip .Le. maxphi_eff) Then
            rotated_ro_P=rotated_ro(iphip,ihil_iphicyl,2)
         Else
            rotated_ro_P=conjg(rotated_ro(maxphi-iphip+2,ihil_iphicyl,2))
         End If
         ro_projN=ro_projN+ephicN(iphin,0)*ephicP(iphip,0)*rotated_overlap(ibet,iphin,1)*rotated_overlap(ibet,iphip,2)*rotated_ro_N
         ro_projP=ro_projP+ephicN(iphin,0)*ephicP(iphip,0)*rotated_overlap(ibet,iphin,1)*rotated_overlap(ibet,iphip,2)*rotated_ro_P
         ro_norm=ro_norm+ephicN(iphin,0)*ephicP(iphip,0)*rotated_overlap(ibet,iphin,1)*rotated_overlap(ibet,iphip,2)
      End Do
      End Do
      ro_projected(ihil_iphicyl,1)=ro_projN/ro_norm; ro_projected(ihil_iphicyl,2)=ro_projP/ro_norm
   End Do
!$OMP End PARALLEL DO
   !
   Do it=1,2
     Do ihil_iphicyl=1,nghl*ngphi
        Do iphi=1,maxphi_eff
           rotated_ro(iphi,ihil_iphicyl,it)=rotated_ro(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
           rotated_tau(iphi,ihil_iphicyl,it)=rotated_tau(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
           rotated_sfiz(iphi,ihil_iphicyl,it)=rotated_sfiz(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
           rotated_sfir(iphi,ihil_iphicyl,it)=rotated_sfir(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
           rotated_srfi(iphi,ihil_iphicyl,it)=rotated_srfi(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
           rotated_szfi(iphi,ihil_iphicyl,it)=rotated_szfi(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
           rotated_szz(iphi,ihil_iphicyl,it)=rotated_szz(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
           rotated_srz(iphi,ihil_iphicyl,it)=rotated_srz(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
           rotated_srr(iphi,ihil_iphicyl,it)=rotated_srr(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
           rotated_szr(iphi,ihil_iphicyl,it)=rotated_szr(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
           rotated_sfifi(iphi,ihil_iphicyl,it)=rotated_sfifi(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
           rotated_aka(iphi,ihil_iphicyl,it)=rotated_aka(iphi,ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
           rotated_dj(iphi,ihil_iphicyl,it)=rotated_dj(iphi,ihil_iphicyl,it)*two*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
           rotated_dro(iphi,ihil_iphicyl,it)=rotated_dro(iphi,ihil_iphicyl,it)*two*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
        End Do
        ro_projected(ihil_iphicyl,it)=ro_projected(ihil_iphicyl,it)*piu(it)*wdcori(ihil_convert(ihil_iphicyl))
     End Do
   End Do
   !
   ! Calculation of the Coulomb field 
   Allocate(vc_pr(nghl,nghl))
   Do ihil1=1,nghl
   Do ihil2=1,nghl
      vc_pr(ihil1,ihil2)=Cmplx(vc(ihil1,ihil2),zero)
   End Do
   End Do
   !
   Do iphi=1,maxphi_eff
      cou_projected(iphi,:)=czero
      Do iphicyl=1,ngphi
         Call Zgemm('N','N',nghl,1,nghl,cone/dble(ngphi),vc_pr,nghl,rotated_ro(iphi,(iphicyl-1)*nghl+1:iphicyl*nghl,2),nghl,cone,&
                    cou_projected(iphi,:),nghl)
      End Do
   End Do
   Deallocate(vc_pr)
   ! ------
   ! Timing
   ! ------
   Call system_clock(tfinish_ibet,clock_rate)
   time_ibet(4)=time_ibet(4)+Real(tfinish_ibet-tstart_ibet)/Real(clock_rate)
   !
  End Subroutine calculate_densities
  !=============================================
  ! Subroutine that calculates rotated energies.
  !=============================================
  Subroutine calculate_energies(ibet)
   Implicit None
   Integer(ipr) ibet,ihil,iphicyl,ihil_iphicyl,iphin,iphip,t,iphin_iphip,stride_size,loc
   Integer(ipr) tstart_ibet,tfinish_ibet
   Complex(pr) rn,rp,tnt,tpt,dn,dp,djn,djp,sfizn,sfizp,sfirn,sfirp,szfin,szfip,srfin,srfip,szzn,szzp,srzn,srzp,&
               srrn,srrp,szrn,szrp,sfifin,sfifip,akn,akp,akn2,akp2,rprojn,rprojp,adn,adp
   Complex(pr) RHO_0,RHO_1,TAU_0,TAU_1,DRHO_0,DRHO_1,DJ_0,DJ_1,SFIZ_0,SFIZ_1,SFIR_0,SFIR_1,SZFI_0,SZFI_1,SRFI_0,SRFI_1,SZZ_0,SZZ_1,&
               SRZ_0,SRZ_1,SRR_0,SRR_1,SZR_0,SZR_1,SFIFI_0,SFIFI_1,J2_0,J2_1,RHOPROJ_0,RHOPROJ_1
   Complex(pr) ekinN_phicyl(ngphi),ekinP_phicyl(ngphi),ecodi_phicyl(ngphi),ecoex_phicyl(ngphi),row1_phicyl(ngphi),row2_phicyl(ngphi),&
            EVOL_rho_tau_phicyl(ngphi),EVOL_rho_rho_phicyl(ngphi),ESURF_rho_drho_phicyl(ngphi),ETENS_phicyl(ngphi),&
            ESO_rho_nablaj_phicyl(ngphi),eptN_phicyl(ngphi),eptP_phicyl(ngphi),delN_phicyl(ngphi),delP_phicyl(ngphi)
   Real(pr) EVOL_rho_tau_real,EVOL_rho_tau_imag,EVOL_rho_rho_real,EVOL_rho_rho_imag,ESURF_rho_drho_real,ESURF_rho_drho_imag,ekinN_real,ekinN_imag,&
            ekinP_real,ekinP_imag,ESO_rho_nablaj_real,ESO_rho_nablaj_imag,ecodi_real,ecodi_imag,ecoex_real,ecoex_imag,ETENS_real,ETENS_imag,&
	    row1_real,row1_imag,row2_real,row2_imag,eptN_real,eptN_imag,eptP_real,eptP_imag,delN_real,delN_imag,delP_real,delP_imag
   Complex(pr) ekinN,ekinP,ecodi,ecoex,EVOL_rho_tau,EVOL_rho_rho,ESURF_rho_drho,ETENS,ESO_rho_nablaj,eptN,eptP,xn1,rms1,q21,xn2,rms2,q22,delN,delP
   Real(pr) z,zz,rrr,p2,twopii,whl
   Complex(pr) Urhorho_pr(0:1),Urhotau_pr(0:1),UrhoDrho_pr(0:1),UJJ_pr(0:1),UrhonablaJ_pr(0:1)
   Complex(pr) rho_coupling,rho_exchange,cou_projected_element
   !
   If(beta_active(ibet).Eq.0) Return
   ! ------
   ! Timing
   ! ------
   Call system_clock(tstart_ibet,clock_rate)
   !
   twopii=one/(two*pi)
   iphin_iphip=0
   Do iphin=1,maxphi
   Do iphip=1,maxphi
      iphin_iphip=iphin_iphip+1
      ekinN=czero;ekinP=czero;ecodi=czero;ecoex=czero;EVOL_rho_tau=czero;EVOL_rho_rho=czero;ESURF_rho_drho=czero;ETENS=czero;ESO_rho_nablaj=czero
      eptN=czero;eptP=czero;xn1=czero;rms1=czero;q21=czero;xn2=czero;rms2=czero;q22=czero;delN=czero;delP=czero
      Do ihil=1,nghl
         whl=wdcor(ihil);
         z=fh(ihil); zz=z*z; rrr=zz+fl(ihil)**2; p2=p32*zz-half*rrr;
         !
         Do iphicyl=1,ngphi
            ihil_iphicyl=ihil_iphicyl_convert(ihil,iphicyl)
            ! Projected density
            rprojn=ro_projected(ihil_iphicyl,1); rprojp=ro_projected(ihil_iphicyl,2)
            RHOPROJ_0=rprojn+rprojp; RHOPROJ_1=rprojn-rprojp
            !------------------------
            ! rms and deformations
            !------------------------
            z=fh(ihil); zz=z*z; rrr=zz+fl(ihil)**2
            p2=p32*zz-half*rrr;
            !------------------
            ! np-representation
            !------------------
            If(iphin .Le. maxphi_eff) Then
               rn=rotated_ro(iphin,ihil_iphicyl,1)
               tnt=rotated_tau(iphin,ihil_iphicyl,1)
               dn=rotated_dro(iphin,ihil_iphicyl,1)
               djn=rotated_dj(iphin,ihil_iphicyl,1)
               sfizn=rotated_sfiz(iphin,ihil_iphicyl,1)
               sfirn=rotated_sfir(iphin,ihil_iphicyl,1)
               szfin=rotated_szfi(iphin,ihil_iphicyl,1)
               srfin=rotated_srfi(iphin,ihil_iphicyl,1)
               szzn=rotated_szz(iphin,ihil_iphicyl,1)
               srzn=rotated_srz(iphin,ihil_iphicyl,1)
               srrn=rotated_srr(iphin,ihil_iphicyl,1)
               szrn=rotated_szr(iphin,ihil_iphicyl,1)
               sfifin=rotated_sfifi(iphin,ihil_iphicyl,1)
               akn=rotated_aka(iphin,ihil_iphicyl,1); akn2=akn*conjg(akn); adn=akn*rn
            Else
               rn=conjg(rotated_ro(maxphi-iphin+2,ihil_iphicyl,1))
               tnt=conjg(rotated_tau(maxphi-iphin+2,ihil_iphicyl,1))
               dn=conjg(rotated_dro(maxphi-iphin+2,ihil_iphicyl,1))
               djn=conjg(rotated_dj(maxphi-iphin+2,ihil_iphicyl,1))
               sfizn=conjg(rotated_sfiz(maxphi-iphin+2,ihil_iphicyl,1))
               sfirn=conjg(rotated_sfir(maxphi-iphin+2,ihil_iphicyl,1))
               szfin=conjg(rotated_szfi(maxphi-iphin+2,ihil_iphicyl,1))
               srfin=conjg(rotated_srfi(maxphi-iphin+2,ihil_iphicyl,1))
               szzn=conjg(rotated_szz(maxphi-iphin+2,ihil_iphicyl,1))
               srzn=conjg(rotated_srz(maxphi-iphin+2,ihil_iphicyl,1))
               srrn=conjg(rotated_srr(maxphi-iphin+2,ihil_iphicyl,1))
               szrn=conjg(rotated_szr(maxphi-iphin+2,ihil_iphicyl,1))
               sfifin=conjg(rotated_sfifi(maxphi-iphin+2,ihil_iphicyl,1))
               akn=conjg(rotated_aka(maxphi-iphin+2,ihil_iphicyl,1)); akn2=akn*conjg(akn); adn=akn*rn
            End If
            !
            If(iphip .Le. maxphi_eff) Then
               rp=rotated_ro(iphip,ihil_iphicyl,2)
               tpt=rotated_tau(iphip,ihil_iphicyl,2)
               dp=rotated_dro(iphip,ihil_iphicyl,2)
               djp=rotated_dj(iphip,ihil_iphicyl,2)
               sfizp=rotated_sfiz(iphip,ihil_iphicyl,2)
               sfirp=rotated_sfir(iphip,ihil_iphicyl,2)
               szfip=rotated_szfi(iphip,ihil_iphicyl,2)
               srfip=rotated_srfi(iphip,ihil_iphicyl,2)
               szzp=rotated_szz(iphip,ihil_iphicyl,2)
               srzp=rotated_srz(iphip,ihil_iphicyl,2)
               srrp=rotated_srr(iphip,ihil_iphicyl,2)
               szrp=rotated_szr(iphip,ihil_iphicyl,2)
               sfifip=rotated_sfifi(iphip,ihil_iphicyl,2)
               akp=rotated_aka(iphip,ihil_iphicyl,2); akp2=akp*conjg(akp); adp=akp*rp
               cou_projected_element=cou_projected(iphip,ihil)
            Else
               rp=conjg(rotated_ro(maxphi-iphip+2,ihil_iphicyl,2))
               tpt=conjg(rotated_tau(maxphi-iphip+2,ihil_iphicyl,2))
               dp=conjg(rotated_dro(maxphi-iphip+2,ihil_iphicyl,2))
               djp=conjg(rotated_dj(maxphi-iphip+2,ihil_iphicyl,2))
               sfizp=conjg(rotated_sfiz(maxphi-iphip+2,ihil_iphicyl,2))
               sfirp=conjg(rotated_sfir(maxphi-iphip+2,ihil_iphicyl,2))
               szfip=conjg(rotated_szfi(maxphi-iphip+2,ihil_iphicyl,2))
               srfip=conjg(rotated_srfi(maxphi-iphip+2,ihil_iphicyl,2))
               szzp=conjg(rotated_szz(maxphi-iphip+2,ihil_iphicyl,2))
               srzp=conjg(rotated_srz(maxphi-iphip+2,ihil_iphicyl,2))
               srrp=conjg(rotated_srr(maxphi-iphip+2,ihil_iphicyl,2))
               szrp=conjg(rotated_szr(maxphi-iphip+2,ihil_iphicyl,2))
               sfifip=conjg(rotated_sfifi(maxphi-iphip+2,ihil_iphicyl,2))
               akp=conjg(rotated_aka(maxphi-iphip+2,ihil_iphicyl,2)); akp2=akp*conjg(akp); adp=akp*rp
               cou_projected_element=conjg(cou_projected(maxphi-iphip+2,ihil))
            End If
            !-----------------
            ! t-representation
            !-----------------
            RHO_0=rn+rp; RHO_1=rn-rp
            TAU_0=tnt+tpt; TAU_1=tnt-tpt
            DRHO_0=dn+dp; DRHO_1=dn-dp
            DJ_0=djn+djp; DJ_1=djn-djp
            SFIZ_0=sfizn+sfizp; SFIZ_1=sfizn-sfizp
            SFIR_0=sfirn+sfirp; SFIR_1=sfirn-sfirp
            SZFI_0=szfin+szfip; SZFI_1=szfin-szfip;
            SRFI_0=srfin+srfip; SRFI_1=srfin-srfip;
            SZZ_0=szzn+szzp; SZZ_1=szzn-szzp
            SRZ_0=srzn+srzp; SRZ_1=srzn-srzp
            SRR_0=srrn+srrp; SRR_1=srrn-srrp
            SZR_0=szrn+szrp; SZR_1=szrn-szrp
            SFIFI_0=sfifin+sfifip; SFIFI_1=sfifin-sfifip
            J2_0=SFIZ_0**2+SFIR_0**2+SZFI_0**2+SRFI_0**2+SZZ_0**2+SRZ_0**2+SRR_0**2+SZR_0**2+SFIFI_0**2
            J2_1=SFIZ_1**2+SFIR_1**2+SZFI_1**2+SRFI_1**2+SZZ_1**2+SRZ_1**2+SRR_1**2+SZR_1**2+SFIFI_1**2
            !
            If(PNP_is_on .Eq. 0 .Or. PNP_is_on .Eq. 1) Then
               ! Mixed prescription
               rho_coupling=RHO_0
               rho_exchange=rp
            Else If(PNP_is_on .Eq. 2) Then
               ! Projected prescription
               rho_coupling=RHOPROJ_0
               rho_exchange=rprojp
            End If
            !--------------------
            ! Coupling parameters
            !--------------------
            Urhorho_pr=czero;Urhotau_pr=czero;UrhoDrho_pr=czero;UJJ_pr=czero;UrhonablaJ_pr=czero
            Do t=0,1
               Urhorho_pr(t)=Cmplx(Crho(t),zero)+Cdrho(t)*rho_coupling**sigma
               Urhotau_pr(t)=Cmplx(Ctau(t),zero)
               UrhoDrho_pr(t)=Cmplx(Crdr(t),zero)
               If(use_j2terms) UJJ_pr(t)=Cmplx(CJ(t),zero)
               UrhonablaJ_pr(t)=Cmplx(Crdj(t),zero)
            End Do
            !------------------------------------------------------
            ! Energies for fixed (z,rperp) point and varying phicyl
            !------------------------------------------------------
            ! Volume energy
            EVOL_rho_tau_phicyl(iphicyl)=Urhotau_pr(0)*RHO_0*TAU_0+Urhotau_pr(1)*RHO_1*TAU_1           ! volume rho tau
            EVOL_rho_rho_phicyl(iphicyl)=Urhorho_pr(0)*RHO_0**2+Urhorho_pr(1)*RHO_1**2                 ! volume density dependent
            ! Surface energy
            ESURF_rho_drho_phicyl(iphicyl)=UrhoDrho_pr(0)*RHO_0*DRHO_0+UrhoDrho_pr(1)*RHO_1*DRHO_1     ! surface rho delta rho
            ! Kinetic energy
            ekinN_phicyl(iphicyl)=hb0n*(TAU_0+TAU_1)*half*facECM                                      ! kinetic neutrons
            ekinP_phicyl(iphicyl)=hb0p*(TAU_0-TAU_1)*half*facECM                                      ! kinetic protons
            ! Spin orbit
            ESO_rho_nablaj_phicyl(iphicyl)=UrhonablaJ_pr(0)*RHO_0*DJ_0+UrhonablaJ_pr(1)*RHO_1*DJ_1    ! spin-orbit rho Nabla . J
            ! Coulomb energy
            ecodi_phicyl(iphicyl)=czero; ecoex_phicyl(iphicyl)=czero
            If (icou.Ge.1) ecodi_phicyl(iphicyl)=half*cou_projected_element*rp                        ! Coulomb direct
            If (icou.Eq.2.Or.icou.Eq.-4) ecoex_phicyl(iphicyl)=CExPar*cex*rho_exchange**p43           ! Coulomb exchange, Slater approximation 
            ! Tensor energy
            ETENS_phicyl(iphicyl)=UJJ_pr(0)*J2_0+UJJ_pr(1)*J2_1
            ! Pairing energy
            If(pairing_regularization .and. All(geff_inv .Ne. 0.0_pr)) Then
!               eptN_phicyl(iphicyl)=akn2
!               eptP_phicyl(iphicyl)=akp2
            Else
               eptN_phicyl(iphicyl)=CpV0(0)*(one-(RHO_0/rho_c)*CpV1(0))*akn2
               eptP_phicyl(iphicyl)=CpV0(1)*(one-(RHO_0/rho_c)*CpV1(1))*akp2
               !
               delN_phicyl(iphicyl)=CpV0(0)*(one-(RHO_0/rho_c)*CpV1(0))*adn
               delP_phicyl(iphicyl)=CpV0(1)*(one-(RHO_0/rho_c)*CpV1(1))*adp
            End If
            ! Rms and deformations
            row1_phicyl(iphicyl)=rn; row2_phicyl(iphicyl)=rp
         End Do ! iphicyl
         !
         ! Volume energy
         Call integrate_simpson(Dble(EVOL_rho_tau_phicyl(:)),ngphi,phicyl_integration_step,EVOL_rho_tau_real)
         Call integrate_simpson(Aimag(EVOL_rho_tau_phicyl(:)),ngphi,phicyl_integration_step,EVOL_rho_tau_imag)
         EVOL_rho_tau=EVOL_rho_tau+Cmplx(EVOL_rho_tau_real,EVOL_rho_tau_imag)*whl
         Call integrate_simpson(Dble(EVOL_rho_rho_phicyl(:)),ngphi,phicyl_integration_step,EVOL_rho_rho_real)
         Call integrate_simpson(Aimag(EVOL_rho_rho_phicyl(:)),ngphi,phicyl_integration_step,EVOL_rho_rho_imag)
         EVOL_rho_rho=EVOL_rho_rho+Cmplx(EVOL_rho_rho_real,EVOL_rho_rho_imag)*whl
         ! Surface energy
         Call integrate_simpson(Dble(ESURF_rho_drho_phicyl(:)),ngphi,phicyl_integration_step,ESURF_rho_drho_real)
         Call integrate_simpson(Aimag(ESURF_rho_drho_phicyl(:)),ngphi,phicyl_integration_step,ESURF_rho_drho_imag)
         ESURF_rho_drho=ESURF_rho_drho+Cmplx(ESURF_rho_drho_real,ESURF_rho_drho_imag)*whl
         ! Kinetic energy
         Call integrate_simpson(Dble(ekinN_phicyl(:)),ngphi,phicyl_integration_step,ekinN_real)
         Call integrate_simpson(Aimag(ekinN_phicyl(:)),ngphi,phicyl_integration_step,ekinN_imag)
         Call integrate_simpson(Dble(ekinP_phicyl(:)),ngphi,phicyl_integration_step,ekinP_real)
         Call integrate_simpson(Aimag(ekinP_phicyl(:)),ngphi,phicyl_integration_step,ekinP_imag)
         ekinN=ekinN+Cmplx(ekinN_real,ekinN_imag)*whl; ekinP=ekinP+Cmplx(ekinP_real,ekinP_imag)*whl
         ! Spin orbit
         Call integrate_simpson(Dble(ESO_rho_nablaj_phicyl(:)),ngphi,phicyl_integration_step,ESO_rho_nablaj_real)
         Call integrate_simpson(Aimag(ESO_rho_nablaj_phicyl(:)),ngphi,phicyl_integration_step,ESO_rho_nablaj_imag)
         ESO_rho_nablaj=ESO_rho_nablaj+Cmplx(ESO_rho_nablaj_real,ESO_rho_nablaj_imag)*whl
         ! Coulomb energy
         Call integrate_simpson(Dble(ecodi_phicyl(:)),ngphi,phicyl_integration_step,ecodi_real)
         Call integrate_simpson(Aimag(ecodi_phicyl(:)),ngphi,phicyl_integration_step,ecodi_imag)
         Call integrate_simpson(Dble(ecoex_phicyl(:)),ngphi,phicyl_integration_step,ecoex_real)
         Call integrate_simpson(Aimag(ecoex_phicyl(:)),ngphi,phicyl_integration_step,ecoex_imag)
         ecodi=ecodi+Cmplx(ecodi_real,ecodi_imag)*whl; ecoex=ecoex-Cmplx(ecoex_real,ecoex_imag)*whl
         ! Tensor energy
         Call integrate_simpson(Dble(ETENS_phicyl(:)),ngphi,phicyl_integration_step,ETENS_real)
         Call integrate_simpson(Aimag(ETENS_phicyl(:)),ngphi,phicyl_integration_step,ETENS_imag)
         ETENS=ETENS+Cmplx(ETENS_real,ETENS_imag)*whl
         ! Pairing energy
         Call integrate_simpson(Dble(eptN_phicyl(:)),ngphi,phicyl_integration_step,eptN_real)
         Call integrate_simpson(Aimag(eptN_phicyl(:)),ngphi,phicyl_integration_step,eptN_imag)
         Call integrate_simpson(Dble(eptP_phicyl(:)),ngphi,phicyl_integration_step,eptP_real)
         Call integrate_simpson(Aimag(eptP_phicyl(:)),ngphi,phicyl_integration_step,eptP_imag)
         !
         Call integrate_simpson(Dble(delN_phicyl(:)),ngphi,phicyl_integration_step,delN_real)
         Call integrate_simpson(Aimag(delN_phicyl(:)),ngphi,phicyl_integration_step,delN_imag)
         Call integrate_simpson(Dble(delP_phicyl(:)),ngphi,phicyl_integration_step,delP_real)
         Call integrate_simpson(Aimag(delP_phicyl(:)),ngphi,phicyl_integration_step,delP_imag)
         If(pairing_regularization .and. All(geff_inv .Ne. 0.0_pr)) Then
!          eptN=eptN+Cmplx(eptN_real,eptN_imag)*whl/geff_inv(ihil,1)    ! Asumming that geff_inv depend on (z,rperp) but not phi
!          eptP=eptP+Cmplx(eptP_real,eptP_imag)*whl/geff_inv(ihil,2)    ! Asumming that geff_inv depend on (z,rperp) but not phi
         Else
          eptN=eptN+Cmplx(eptN_real,eptN_imag)*whl
          eptP=eptP+Cmplx(eptP_real,eptP_imag)*whl
          !
          delN=delN+Cmplx(delN_real,delN_imag)*whl
          delP=delP+Cmplx(delP_real,delP_imag)*whl
         End If
         ! Rms and deformation
         Call integrate_simpson(Dble(row1_phicyl(:)),ngphi,phicyl_integration_step,row1_real)
         Call integrate_simpson(Aimag(row1_phicyl(:)),ngphi,phicyl_integration_step,row1_imag)
         Call integrate_simpson(Dble(row2_phicyl(:)),ngphi,phicyl_integration_step,row2_real)
         Call integrate_simpson(Aimag(row2_phicyl(:)),ngphi,phicyl_integration_step,row2_imag)
         xn1=xn1+Cmplx(row1_real,row1_imag)*whl; rms1=rms1+Cmplx(row1_real,row1_imag)*rrr*whl; q21=q21+Cmplx(row1_real,row1_imag)*p2*whl
         xn2=xn2+Cmplx(row2_real,row2_imag)*whl; rms2=rms2+Cmplx(row2_real,row2_imag)*rrr*whl; q22=q22+Cmplx(row2_real,row2_imag)*p2*whl
         !
      End Do ! ihil
      !
      stride_size=17
      loc=(ibet-team_rank*beta_size-1)*maxphi*maxphi*stride_size+(iphin_iphip-1)*stride_size
      all_energies(loc+1)=ekinN*twopii
      all_energies(loc+2)=ekinP*twopii
      all_energies(loc+3)=ecodi*twopii
      all_energies(loc+4)=ecoex*twopii
      all_energies(loc+5)=EVOL_rho_tau*twopii
      all_energies(loc+6)=EVOL_rho_rho*twopii
      all_energies(loc+7)=ESURF_rho_drho*twopii
      all_energies(loc+8)=ETENS*twopii
      all_energies(loc+9)=ESO_rho_nablaj*twopii
      all_energies(loc+10)=eptN*twopii
      all_energies(loc+11)=eptP*twopii
      all_energies(loc+12)=xn1*twopii
      all_energies(loc+13)=xn2*twopii
      all_energies(loc+14)=Sqrt(rms1/xn1)
      all_energies(loc+15)=Sqrt(rms2/xn2)
      all_energies(loc+16)=delN/xn1
      all_energies(loc+17)=delP/xn2
      !
   End Do ! iphip
   End Do ! iphin
   ! ------
   ! Timing
   ! ------
   Call system_clock(tfinish_ibet,clock_rate)
   time_ibet(5)=time_ibet(5)+Real(tfinish_ibet-tstart_ibet)/Real(clock_rate)
   !
  End Subroutine calculate_energies
  !===============================================================
  ! Subroutine that deallocates quantities specific to each angle.
  !===============================================================
  Subroutine finalize_angle()
    Implicit None
    !
    Deallocate(rotation_matrix,inverse_rotation_matrix,rotated_density,rotated_kappa)
    Deallocate(rotated_ro,rotated_tau,rotated_dj,rotated_dro,rotated_sfiz,rotated_sfir,rotated_srfi,rotated_szfi,rotated_szz,&
               rotated_srz,rotated_srr,rotated_szr,rotated_sfifi,rotated_aka,ro_projected)
    Deallocate(iosc1_contributing,iosc2_contributing,nt_contributing)
    Deallocate(cou_projected)
    !
  End Subroutine finalize_angle
  !==============================================================
  ! Subroutine that calculates the symmetry-projected quantities.
  !==============================================================
  Subroutine project()
    Implicit None
    Integer(ipr) ibet,iphi,jj,iphin,iphip,jjmax,ineu,ipro,nn_min,nn_max,nn_neut,nn_prot,icount,itask,i_beta,it
    Real(pr) beta,facj,betafac,num_neu,num_pro
    Complex(pr) fac_overlap,cj,ekinN,ekinP,ecodi,ecoex,EVOL_rho_tau,EVOL_rho_rho,ESURF_rho_drho,ETENS,ESO_rho_nablaj,&
                eptN,eptP,xn1,xn2,rms1,rms2,delN,delP
    Complex(pr) cj_phi,ekinN_phi,ekinP_phi,ecodi_phi,ecoex_phi,EVOL_rho_tau_phi,EVOL_rho_rho_phi,ESURF_rho_drho_phi,ETENS_phi,ESO_rho_nablaj_phi,&
                eptN_phi,eptP_phi,xn1_phi,xn2_phi,rms1_phi,rms2_phi,delN_phi,delP_phi,cj_jj
    Complex(pr), Allocatable:: cj_NP(:,:)
    Integer(ipr) tstart_ibet,tfinish_ibet
    Complex(pr), Allocatable :: ekinN_rotated(:,:,:),ekinP_rotated(:,:,:),ecodi_rotated(:,:,:),ecoex_rotated(:,:,:),EVOL_rho_tau_rotated(:,:,:),&
                                EVOL_rho_rho_rotated(:,:,:),ESURF_rho_drho_rotated(:,:,:),ETENS_rotated(:,:,:),ESO_rho_nablaj_rotated(:,:,:),eptN_rotated(:,:,:),&
                                eptP_rotated(:,:,:),xn1_rotated(:,:,:),xn2_rotated(:,:,:),rms1_rotated(:,:,:),rms2_rotated(:,:,:),delN_rotated(:,:,:),delP_rotated(:,:,:)
    !
    ! ------
    ! Timing
    ! ------
    Call system_clock(tstart_ibet,clock_rate)
    !
    If(AMP_is_on .Eq. 1) Then
       jjmax=maxj
    Else
       jjmax=0
    End If
    ! ----------------------------------------
    ! Unwrapping rotated overlaps and energies
    ! ---------------------------------------
    Allocate(ekinN_rotated(maxbet,maxphi,maxphi),ekinP_rotated(maxbet,maxphi,maxphi),ecodi_rotated(maxbet,maxphi,maxphi),&
             ecoex_rotated(maxbet,maxphi,maxphi),EVOL_rho_tau_rotated(maxbet,maxphi,maxphi),EVOL_rho_rho_rotated(maxbet,maxphi,maxphi),&
             ESURF_rho_drho_rotated(maxbet,maxphi,maxphi),ETENS_rotated(maxbet,maxphi,maxphi),ESO_rho_nablaj_rotated(maxbet,maxphi,maxphi),&
             eptN_rotated(maxbet,maxphi,maxphi),eptP_rotated(maxbet,maxphi,maxphi),xn1_rotated(maxbet,maxphi,maxphi),xn2_rotated(maxbet,maxphi,maxphi),&
             rms1_rotated(maxbet,maxphi,maxphi),rms2_rotated(maxbet,maxphi,maxphi),delN_rotated(maxbet,maxphi,maxphi),delP_rotated(maxbet,maxphi,maxphi))
    !
    icount=1
    Do itask=1,maxbet/beta_size
       Do it=1,2
          Do iphi=1,maxphi
             Do i_beta=1,beta_size
                ibet=(itask-1)*beta_size+i_beta
                rotated_overlap(ibet,iphi,it)=all_overlaps_gthr(icount); icount=icount+1
             End Do
          End Do
       End Do
    End Do
    !
    icount=1
    Do itask=1,maxbet/beta_size
       Do i_beta=1,beta_size
          ibet=(itask-1)*beta_size+i_beta
          Do iphin=1,maxphi
          Do iphip=1,maxphi
             ekinN_rotated(ibet,iphin,iphip)=all_energies_gthr(icount); icount=icount+1
             ekinP_rotated(ibet,iphin,iphip)=all_energies_gthr(icount); icount=icount+1
             ecodi_rotated(ibet,iphin,iphip)=all_energies_gthr(icount); icount=icount+1
             ecoex_rotated(ibet,iphin,iphip)=all_energies_gthr(icount); icount=icount+1
             EVOL_rho_tau_rotated(ibet,iphin,iphip)=all_energies_gthr(icount); icount=icount+1
             EVOL_rho_rho_rotated(ibet,iphin,iphip)=all_energies_gthr(icount); icount=icount+1
             ESURF_rho_drho_rotated(ibet,iphin,iphip)=all_energies_gthr(icount); icount=icount+1
             ETENS_rotated(ibet,iphin,iphip)=all_energies_gthr(icount); icount=icount+1
             ESO_rho_nablaj_rotated(ibet,iphin,iphip)=all_energies_gthr(icount); icount=icount+1
             eptN_rotated(ibet,iphin,iphip)=all_energies_gthr(icount); icount=icount+1
             eptP_rotated(ibet,iphin,iphip)=all_energies_gthr(icount); icount=icount+1
             xn1_rotated(ibet,iphin,iphip)=all_energies_gthr(icount); icount=icount+1
             xn2_rotated(ibet,iphin,iphip)=all_energies_gthr(icount); icount=icount+1
             rms1_rotated(ibet,iphin,iphip)=all_energies_gthr(icount); icount=icount+1
             rms2_rotated(ibet,iphin,iphip)=all_energies_gthr(icount); icount=icount+1
             delN_rotated(ibet,iphin,iphip)=all_energies_gthr(icount); icount=icount+1
             delP_rotated(ibet,iphin,iphip)=all_energies_gthr(icount); icount=icount+1
          End Do
          End Do
       End Do
    End Do ! itask
    !---------------------------------------------------------------------------
    ! Projecting on desired number of particles and on different angular momenta
    !---------------------------------------------------------------------------
    Do jj=0,jjmax,jjstep
       facj=half*(two*jj+one)*integration_prefactor
       cj=czero;ekinN=czero;ekinP=czero;ecodi=czero;ecoex=czero;EVOL_rho_tau=czero;EVOL_rho_rho=czero;ESURF_rho_drho=czero;ETENS=czero
       ESO_rho_nablaj=czero;eptN=czero;eptP=czero;xn1=czero;xn2=czero;rms1=czero;rms2=czero;delN=czero;delP=czero
       ! Angular momentum projection
       Do ibet=1,maxbet
          beta=betabs(ibet);betafac=betaweight(ibet)*wigner(jj,0,0,beta)
          ! Particle number projection
          cj_phi=czero;ekinN_phi=czero;ekinP_phi=czero;ecodi_phi=czero;ecoex_phi=czero;EVOL_rho_tau_phi=czero;EVOL_rho_rho_phi=czero;ESURF_rho_drho_phi=czero
          ETENS_phi=czero;ESO_rho_nablaj_phi=czero;eptN_phi=czero;eptP_phi=czero;xn1_phi=czero;xn2_phi=czero;rms1_phi=czero;rms2_phi=czero;delN_phi=czero;delP_phi=czero
          Do iphin=1,maxphi
          Do iphip=1,maxphi
             fac_overlap=rotated_overlap(ibet,iphin,1)*rotated_overlap(ibet,iphip,2)*ephicN(iphin,0)*ephicP(iphip,0)
             cj_phi=cj_phi+fac_overlap
             ekinN_phi=ekinN_phi+fac_overlap*ekinN_rotated(ibet,iphin,iphip)
             ekinP_phi=ekinP_phi+fac_overlap*ekinP_rotated(ibet,iphin,iphip)
             ecodi_phi=ecodi_phi+fac_overlap*ecodi_rotated(ibet,iphin,iphip)
             ecoex_phi=ecoex_phi+fac_overlap*ecoex_rotated(ibet,iphin,iphip)
             EVOL_rho_tau_phi=EVOL_rho_tau_phi+fac_overlap*EVOL_rho_tau_rotated(ibet,iphin,iphip)
             EVOL_rho_rho_phi=EVOL_rho_rho_phi+fac_overlap*EVOL_rho_rho_rotated(ibet,iphin,iphip)
             ESURF_rho_drho_phi=ESURF_rho_drho_phi+fac_overlap*ESURF_rho_drho_rotated(ibet,iphin,iphip)
             ETENS_phi=ETENS_phi+fac_overlap*ETENS_rotated(ibet,iphin,iphip)
             ESO_rho_nablaj_phi=ESO_rho_nablaj_phi+fac_overlap*ESO_rho_nablaj_rotated(ibet,iphin,iphip)
             eptN_phi=eptN_phi+fac_overlap*eptN_rotated(ibet,iphin,iphip)
             eptP_phi=eptP_phi+fac_overlap*eptP_rotated(ibet,iphin,iphip)
             xn1_phi=xn1_phi+fac_overlap*xn1_rotated(ibet,iphin,iphip)
             xn2_phi=xn2_phi+fac_overlap*xn2_rotated(ibet,iphin,iphip)
             rms1_phi=rms1_phi+fac_overlap*rms1_rotated(ibet,iphin,iphip)
             rms2_phi=rms2_phi+fac_overlap*rms2_rotated(ibet,iphin,iphip)
             delN_phi=delN_phi+fac_overlap*delN_rotated(ibet,iphin,iphip)
             delP_phi=delP_phi+fac_overlap*delP_rotated(ibet,iphin,iphip)
             !
          End Do ! iphip
          End Do ! iphin
          cj_phi=cj_phi/maxphi**2; ekinN_phi=ekinN_phi/maxphi**2; ekinP_phi=ekinP_phi/maxphi**2; ecodi_phi=ecodi_phi/maxphi**2; ecoex_phi=ecoex_phi/maxphi**2
          EVOL_rho_tau_phi=EVOL_rho_tau_phi/maxphi**2; EVOL_rho_rho_phi=EVOL_rho_rho_phi/maxphi**2; ESURF_rho_drho_phi=ESURF_rho_drho_phi/maxphi**2
          ETENS_phi=ETENS_phi/maxphi**2; ESO_rho_nablaj_phi=ESO_rho_nablaj_phi/maxphi**2; eptN_phi=eptN_phi/maxphi**2; eptP_phi=eptP_phi/maxphi**2
          xn1_phi=xn1_phi/maxphi**2; xn2_phi=xn2_phi/maxphi**2; rms1_phi=rms1_phi/maxphi**2; rms2_phi=rms2_phi/maxphi**2
          delN_phi=delN_phi/maxphi**2; delP_phi=delP_phi/maxphi**2
          !
          If(AMP_is_on .Eq. 1) Then
             cj=cj+cj_phi*betafac
             ekinN=ekinN+betafac*ekinN_phi; ekinP=ekinP+betafac*ekinP_phi
             ecodi=ecodi+betafac*ecodi_phi; ecoex=ecoex+betafac*ecoex_phi
             EVOL_rho_tau=EVOL_rho_tau+betafac*EVOL_rho_tau_phi
             EVOL_rho_rho=EVOL_rho_rho+betafac*EVOL_rho_rho_phi
             ESURF_rho_drho=ESURF_rho_drho+betafac*ESURF_rho_drho_phi
             ETENS=ETENS+betafac*ETENS_phi
             ESO_rho_nablaj=ESO_rho_nablaj+betafac*ESO_rho_nablaj_phi
             eptN=eptN+betafac*eptN_phi
             eptP=eptP+betafac*eptP_phi
             xn1=xn1+betafac*xn1_phi; xn2=xn2+betafac*xn2_phi
             rms1=rms1+betafac*rms1_phi; rms2=rms2+betafac*rms2_phi
             delN=delN+betafac*delN_phi; delP=delP+betafac*delP_phi
          Else
             projected_overlap(jj)=Dble(cj_phi)
             projected_ekinN(jj)=Dble(ekinN_phi)/projected_overlap(jj);projected_ekinP(jj)=Dble(ekinP_phi)/projected_overlap(jj)
             projected_ecodi(jj)=Dble(ecodi_phi)/projected_overlap(jj);projected_ecoex(jj)=Dble(ecoex_phi)/projected_overlap(jj)
             projected_EVOL_rho_tau(jj)=Dble(EVOL_rho_tau_phi)/projected_overlap(jj)
             projected_EVOL_rho_rho(jj)=Dble(EVOL_rho_rho_phi)/projected_overlap(jj)
             projected_ESURF_rho_drho(jj)=Dble(ESURF_rho_drho_phi)/projected_overlap(jj)
             projected_ETENS(jj)=Dble(ETENS_phi)/projected_overlap(jj)
             projected_ESO_rho_nablaj(jj)=Dble(ESO_rho_nablaj_phi)/projected_overlap(jj)
             projected_eptN(jj)=Dble(eptN_phi)/projected_overlap(jj);projected_eptP(jj)=Dble(eptP_phi)/projected_overlap(jj)
             projected_xn1(jj)=Dble(xn1_phi)/projected_overlap(jj);projected_xn2(jj)=Dble(xn2_phi)/projected_overlap(jj)
             projected_rms1(jj)=Dble(rms1_phi)/projected_overlap(jj);projected_rms2(jj)=Dble(rms2_phi)/projected_overlap(jj)
             projected_delN(jj)=Dble(delN_phi)/projected_overlap(jj);projected_delP(jj)=Dble(delP_phi)/projected_overlap(jj)
          End If
          !
       End Do ! ibet
       !
       If(AMP_is_on .Eq. 1) Then
          projected_overlap(jj)=Dble(cj)*facj
          projected_ekinN(jj)=Dble(ekinN)*facj/projected_overlap(jj);projected_ekinP(jj)=Dble(ekinP)*facj/projected_overlap(jj)
          projected_ecodi(jj)=Dble(ecodi)*facj/projected_overlap(jj);projected_ecoex(jj)=Dble(ecoex)*facj/projected_overlap(jj)
          projected_EVOL_rho_tau(jj)=Dble(EVOL_rho_tau)*facj/projected_overlap(jj)
          projected_EVOL_rho_rho(jj)=Dble(EVOL_rho_rho)*facj/projected_overlap(jj)
          projected_ESURF_rho_drho(jj)=Dble(ESURF_rho_drho)*facj/projected_overlap(jj)
          projected_ETENS(jj)=Dble(ETENS)*facj/projected_overlap(jj)
          projected_ESO_rho_nablaj(jj)=Dble(ESO_rho_nablaj)*facj/projected_overlap(jj)
          projected_eptN(jj)=Dble(eptN)*facj/projected_overlap(jj);projected_eptP(jj)=Dble(eptP)*facj/projected_overlap(jj)
          projected_xn1(jj)=Dble(xn1)*facj/projected_overlap(jj);projected_xn2(jj)=Dble(xn2)*facj/projected_overlap(jj)
          projected_rms1(jj)=Dble(rms1)*facj/projected_overlap(jj);projected_rms2(jj)=Dble(rms2)*facj/projected_overlap(jj)
          projected_delN(jj)=Dble(delN)*facj/projected_overlap(jj);projected_delP(jj)=Dble(delP)*facj/projected_overlap(jj)
       End If
       !
    End Do ! jj
    !
    !---------------------------------------------
    ! Projecting on different numbers of particles
    !---------------------------------------------
    If(PNP_is_on .Gt. 0) Then
    Allocate(cj_NP(-maxN/2:maxN/2,-maxP/2:maxP/2))
       Do jj=0,jjmax,jjstep
          facj=half*(two*jj+one)*integration_prefactor
          cj_NP=czero
          Do ibet=1,maxbet
             beta=betabs(ibet);betafac=betaweight(ibet)*wigner(jj,0,0,beta)
             Do nn_neut=-maxN/2,maxN/2
             Do nn_prot=-maxP/2,maxP/2
                cj_phi=czero
                Do iphin=1,maxphi
                Do iphip=1,maxphi
                   fac_overlap=rotated_overlap(ibet,iphin,1)*rotated_overlap(ibet,iphip,2)*ephicN(iphin,nn_neut)*ephicP(iphip,nn_prot)
                   cj_phi=cj_phi+fac_overlap
                End Do
                End Do
                cj_phi=cj_phi/maxphi**2
                If(AMP_is_on .Eq. 1) Then
                   cj_NP(nn_neut,nn_prot)=cj_NP(nn_neut,nn_prot)+cj_phi*betafac
                Else
                   projected_NP(jj,nn_neut,nn_prot)=Dble(cj_phi)
                End If
             End Do ! nn_neut
             End Do ! nn_prot
          End Do ! ibet
          If(AMP_is_on .Eq. 1) Then
             Do nn_neut=-maxN/2,maxN/2
             Do nn_prot=-maxP/2,maxP/2
                projected_NP(jj,nn_neut,nn_prot)=Dble(cj_NP(nn_neut,nn_prot))*facj
             End Do
             End Do
          End If
       End Do ! jj
       !
       Do jj=0,jjmax,jjstep
          cj_jj=czero
          Do nn_neut=-maxN/2,maxN/2
          Do nn_prot=-maxP/2,maxP/2
             cj_jj=cj_jj+projected_NP(jj,nn_neut,nn_prot)
          End Do
          End Do
          projected_NP_norm(jj)=Dble(cj_jj)
       End Do
       !
    Deallocate(cj_NP)
    End If
    !
    Deallocate(ekinN_rotated,ekinP_rotated,ecodi_rotated,ecoex_rotated,EVOL_rho_tau_rotated,EVOL_rho_rho_rotated,ESURF_rho_drho_rotated,ETENS_rotated,&
               ESO_rho_nablaj_rotated,eptN_rotated,eptP_rotated,xn1_rotated,xn2_rotated,rms1_rotated,rms2_rotated,delN_rotated,delP_rotated)

   ! ------
   ! Timing
   ! ------
   Call system_clock(tfinish_ibet,clock_rate)
   time_ibet(6)=time_ibet(6)+Real(tfinish_ibet-tstart_ibet)/Real(clock_rate)
    !
  End Subroutine project
  !==================================================================
  ! Subroutine that prints to file the symmetry-projected quantities.
  !==================================================================
  Subroutine print_project()
    Implicit None
    Integer(ipr) ibet,jj,beta_active_sum,jjmax,iphin,iphip,iphi,nn_neut,nn_prot
    Real(pr) cjtot,cnztot,ekt(3),evol,esurf,eso,ecoul,etens,eproj,ept(3),phi_sum,cjnorm,cjnormtot,cnznorm,cnznormtot
    Real(pr) xn(3),rms(3),del(2),r212,r222,RpRMSsq,RnRMSsq,DarwinFoldy,rc
    Logical large_printout
    Integer(ipr) lenproj
    Character(50) enfname
    !
    ! Angular momentum projection only
    If(AMP_is_on .Eq. 1 .and. PNP_is_on .Eq. 0) Then
       Write(lproj,'(a)')     '  ============================================================================'
       Write(lproj,'(a)')     '  =================       Mesh in rotational angle beta      ================='
       Write(lproj,'(a)')     '  ============================================================================'
       Write(lproj,'(3x,a,4x,a,2x,a,6x,a,12x,a)') 'ibet','angle','active','overlap neutrons','overlap protons'
       Write(lproj,'(3x,a,4x,a,2x,a,5x,a,5x,a)') '---','-------','---','-----------------------','-----------------------'
       beta_active_sum=0
       Do ibet=1,maxbet
          Write(lproj,'(2x,i3,3x,f9.6,2x,i2,3x,f12.9,x,f12.9,"i",2x,f12.9,x,f12.9,"i")') ibet,betabs(ibet),beta_active(ibet),rotated_overlap(ibet,1,1),&
                                                                  rotated_overlap(ibet,1,2)
          beta_active_sum=beta_active_sum+beta_active(ibet)
       End Do
       Write(lproj,'(19x,a)') '---'
       Write(lproj,'(18x,i3)') beta_active_sum
       Write(lproj,'(a)')     '  ============================================================================'
       Write(lproj,'(a)')     '  =============     Angular momentum content of an HFB state     ============='
       Write(lproj,'(a)')     '  ============================================================================'
       Write(lproj,'(22x,a,10x,a,7x,a)') 'J','|cJ|^2','partial sum'
       Write(lproj,'(21x,a,8x,a,7x,a)') '---','--------','----------'
       !
       cjtot=zero
       Do jj=0,maxj,jjstep
          cjtot=cjtot+abs(projected_overlap(jj))
          Write(lproj,'(20x,i3,7x,f10.7,7x,f10.7)')  jj,projected_overlap(jj),cjtot
       End Do
       Write(lproj,'(47x,a)') ' ---------'
       Write(lproj,'(47x,f10.7)') cjtot
    End If
    ! Particle number projection only
    If(AMP_is_on .Eq. 0 .and. PNP_is_on .Gt. 0) Then
       Write(lproj,'(a)')     '  ============================================================================'
       Write(lproj,'(a)')     '  ===================       Mesh in gauge angle phi        ==================='
       Write(lproj,'(a)')     '  ============================================================================'
       Write(lproj,'(3x,a,4x,a,10x,a,12x,a)') 'iphi','angle','overlap neutrons','overlap protons'
       Write(lproj,'(3x,a,3x,a,6x,a,5x,a)') '---','--------','-----------------------','-----------------------'
       Do iphi=1,maxphi
          Write(lproj,'(2x,i3,3x,f9.6,3x,f12.9,x,f12.9,"i",2x,f12.9,x,f12.9,"i")') iphi,phiabs(iphi),rotated_overlap(1,iphi,1),&
                                                                                                      rotated_overlap(1,iphi,2)
       End Do
       Write(lproj,'(a)')     '  ============================================================================'
       Write(lproj,'(a)')     '  ===========     Particle-number decomposition of an HFB state    ==========='
       Write(lproj,'(a)')     '  ============================================================================'
       Write(lproj,'(20x,a,7x,a,7x,a,5x,a)') 'N', 'Z', '|cNZ|^2', 'partial sum'
       Write(lproj,'(19x,a,5x,a,6x,a,5x,a)') '---', '---', '---------', '---------'
       cnztot=zero
       Do nn_neut=-maxN/2,maxN/2
       Do nn_prot=-maxP/2,maxP/2
          cnztot=cnztot+projected_NP(0,nn_neut,nn_prot)
          Write(lproj,'(16x,f6.2,2x,f6.2,3x,f12.9,2x,f12.9)') tz(1)+nn_neut*2, tz(2)+nn_prot*2, projected_NP(0,nn_neut,nn_prot),cnztot
       End Do
       End Do
       Write(lproj,'(50x,a)') '---------'
       Write(lproj,'(47x,f12.9)') cnztot
    End If
    ! Angular momentum and Particle number projection
    If(AMP_is_on .Eq. 1 .and. PNP_is_on .Gt. 0) Then
       Write(lproj,'(a)')     '  ============================================================================'
       Write(lproj,'(a)')     '  ==========   Mesh in rotational angle beta and gauge angle phi    =========='
       Write(lproj,'(a)')     '  ============================================================================'
       Write(lproj,'(3x,a,4x,a,3x,a,2x,a,4x,a,10xa,13x,a)') 'ibet','beta','active','iphi','phi','overlap neutrons','overlap protons'
       Write(lproj,'(3x,a,4x,a,2x,a,4x,a,3x,a,5x,a,5x,a)') '---','-------','---','---','-------','-----------------------','-----------------------'
       beta_active_sum=0
       Do ibet=1,maxbet
       beta_active_sum=beta_active_sum+beta_active(ibet)
       Do iphi=1,maxphi
          Write(lproj,'(2x,i3,3x,f9.6,2x,i2,4x,i3,3x,f9.6,2x,f12.9,x,f12.9,"i",2x,f12.9,x,f12.9,"i")') ibet,betabs(ibet),beta_active(ibet),&
                iphi,phiabs(iphi),rotated_overlap(ibet,iphi,1),rotated_overlap(ibet,iphi,2)
       End Do
       End Do
       Write(lproj,'(19x,a)') '---'
       Write(lproj,'(18x,i3)') beta_active_sum
       Write(lproj,'(a)')     '  ============================================================================'
       Write(lproj,'(a)')     '  =============    Angular momentum content of PNP-HFB state     ============='
       Write(lproj,'(a)')     '  ============================================================================'
       Write(lproj,'(7x,a,10x,a,7x,a,6x,a,4x,a)') 'J','|cJ|^2','partial sum','|cJ|^2 (norm)','partial sum'
       Write(lproj,'(6x,a,8x,a,7x,a,7x,a,7x,a)') '---','--------','----------','----------','----------'
       cjnorm=zero
       Do jj=0,maxj,jjstep
          cjnorm=cjnorm+abs(projected_overlap(jj))
       End Do
       cjtot=zero;cjnormtot=zero
       Do jj=0,maxj,jjstep
          cjtot=cjtot+abs(projected_overlap(jj))
          cjnormtot=cjnormtot+abs(projected_overlap(jj))/cjnorm
          Write(lproj,'(5x,i3,7x,f10.7,7x,f10.7,7x,f10.7,7x,f10.7)')  jj,projected_overlap(jj),cjtot,projected_overlap(jj)/cjnorm,cjnormtot
       End Do
       Write(lproj,'(32x,a,24x,a)') ' ---------',' ---------'
       Write(lproj,'(32x,f10.7,24x,f10.7)') cjtot,cjnormtot
       Write(lproj,'(a)')     '  ============================================================================'
       Write(lproj,'(a)')     '  =======     Particle-number decomposition of AMP-projected states    ======='
       Write(lproj,'(a)')     '  ============================================================================'
       jjmax=maxj
       If(maxj.gt.40) jjmax=40
       Do jj=0,jjmax,jjstep
          Write(lproj, '(32x,a10)') '----------'
          Write(lproj,'(30x,a7,x,i3,a2)') '| J = ',jj,' |'
          Write(lproj, '(32x,a10)') '----------'
          Write(lproj,'(6x,a)') '--------------------------------------------------------------------'
          Write(lproj,'(7x,a,7x,a,7x,a,5x,a,2x,a,x,a)') 'N', 'Z', '|cNZ|^2', 'partial sum','|cNZ|^2 (norm)','partial sum'
          Write(lproj,'(6x,a,5x,a,6x,a,5x,a,5x,a,5x,a)') '---', '---', '---------', '---------','---------', '---------'
          cnznorm=zero
          Do nn_neut=-maxN/2,maxN/2
          Do nn_prot=-maxP/2,maxP/2
             cnznorm=cnznorm+projected_NP(jj,nn_neut,nn_prot)
          End Do
          End Do
          cnztot=zero;cnznormtot=zero
          Do nn_neut=-maxN/2,maxN/2
          Do nn_prot=-maxP/2,maxP/2
             cnztot=cnztot+projected_NP(jj,nn_neut,nn_prot)
             cnznormtot=cnznormtot+projected_NP(jj,nn_neut,nn_prot)/cnznorm
             Write(lproj,'(3x,f6.2,2x,f6.2,3x,f12.9,2x,f12.9,2x,f12.9,2x,f12.9)') tz(1)+nn_neut*2, tz(2)+nn_prot*2,&
                   projected_NP(jj,nn_neut,nn_prot),cnztot,projected_NP(jj,nn_neut,nn_prot)/cnznorm,cnznormtot
          End Do
          End Do
          Write(lproj,'(37x,a,19x,a)') '---------','---------'
          Write(lproj,'(34x,f12.9,16x,f12.9)') cnztot,cnznormtot
       End Do
       !
    End If
    !
    Write(lproj,'(a)')     '  ============================================================================'
    Write(lproj,'(a)')     '  ====================     Symmetry-projected energies    ===================='
    Write(lproj,'(a)')     '  ============================================================================'
    If(AMP_is_on .Eq. 1) Then
       jjmax=maxj
       If(maxj.Gt.60) jjmax=60
    Else
       jjmax=0
    End If
    Do jj=0,jjmax,jjstep
       ! Energies
       ekt(1)=projected_ekinN(jj); ekt(2)=projected_ekinP(jj); ekt(3)=ekt(1)+ekt(2)
       evol=projected_EVOL_rho_tau(jj)+projected_EVOL_rho_rho(jj); esurf=projected_ESURF_rho_drho(jj)
       eso=projected_ESO_rho_nablaj(jj); ecoul=projected_ecodi(jj)+projected_ecoex(jj); etens=projected_ETENS(jj)
       ept(1)=projected_eptN(jj); ept(2)=projected_eptP(jj); ept(3)=ept(1)+ept(2)
       eproj=ekt(3)+evol+esurf+eso+ecoul+etens+ept(3)
       ! Radii
       xn(1)=projected_xn1(jj); xn(2)=projected_xn2(jj); xn(3)=xn(1)+xn(2)
       rms(1)=projected_rms1(jj); rms(2)=projected_rms2(jj)
       r212=rms(1)**2; r222=rms(2)**2
       rms(3)=Sqrt((xn(1)*r212+xn(2)*r222)/amas)
       del(1)=projected_delN(jj); del(2)=projected_delP(jj)
       ! Charge radius, from Adv. Nucl. Phys. 8, 219 (1975)
       RpRMSsq=0.769_pr
       RnRMSsq=-0.1161_pr   ! J. Phys. G 33, 1 (2006)
       DarwinFoldy=0.033_pr ! Phys. Rev. A 56, 4579 (1997)
       rc=Sqrt(r222+RpRMSsq+(xn(1)/xn(2))*RnRMSsq+DarwinFoldy)
       !
       If(AMP_is_on .Eq. 1) Then
          Write(lproj, '(12x,a10)') '----------'
          Write(lproj,'(10x,a7,x,i3,a2)') '| J = ',jj,' |'
          Write(lproj, '(12x,a10)') '----------'
       End If
       Write(lproj,'(a,3f15.6)') '  rms-radius ..........',rms
       Write(lproj,'(a,15x,2f15.6)') '  charge-radius, r0 ...',rc,r00
       Write(lproj,'(a,2f15.6)')    '  delta(n,p) ..........',del(1),del(2)
       Write(lproj,'(a,3f15.6)')    '  kinetic energy ......',ekt(1),ekt(2),ekt(3)
       Write(lproj,'(a,3f15.6)')    '  pairing energy ......',ept
       Write(lproj,'(a,30x,f15.6)') '  volume energy .......',evol
       Write(lproj,'(a,30x,f15.6)') '        rho_tau .......',projected_EVOL_rho_tau(jj)
       Write(lproj,'(a,30x,f15.6)') '        rho_rho .......',projected_EVOL_rho_rho(jj)
       Write(lproj,'(a,30x,f15.6)') '  surface energy ......',esurf
       Write(lproj,'(a,30x,f15.6)') '   rho_DELTA_rho ......',projected_ESURF_rho_drho(jj)
       Write(lproj,'(a,30x,f15.6)') '  spin-orbit energy ...',eso
       Write(lproj,'(a,30x,f15.6)') '        rho_NABLA_J ...',projected_ESO_rho_nablaj(jj)
       Write(lproj,'(a,30x,f15.6)') '  coulomb energy ......',ecoul
       Write(lproj,'(a,30x,f15.6)') '          direct ......',projected_ecodi(jj)
       Write(lproj,'(a,30x,f15.6)') '          exchange ....',projected_ecoex(jj)
       Write(lproj,'(a,30x,f15.6)') '  tensor energy .......',etens
       Write(lproj,'(a,30x,f15.6)') '  Projected Energy ....',eproj
       !
    End Do
    ! ------
    ! Timing
    ! ------
    Call system_clock(tfinish,clock_rate)
    Write(lproj,'(a)')     '  ============================================================================'
    Write(lproj,'(a)')     '  =======================      System clock times      ======================='
    Write(lproj,'(a)')     '  ============================================================================'
    Write(lproj,'(43x,a,7x,a,8x,a)') 'seconds','minutes','hours'
    Write(lproj,'(43x,a,7x,a,7x,a)') '-------','-------','-------'
    Write(lproj,'(a37,x,f12.2,2x,f12.2,2x,f12.2)') '  Total clock time:                  ', Real(tfinish-tstart)/Real(clock_rate),Real(tfinish-tstart)/(Real(clock_rate)*60_pr),&
									          Real(tfinish-tstart)/(Real(clock_rate)*60_pr*60_pr)
    Write(lproj,'(43x,a)') '-----------------------------------'
    Write(lproj,'(a37,x,f12.2,2x,f12.2,2x,f12.2)') '  Time in initialization routines:   ',time_ibet(1),time_ibet(1)/60_pr,time_ibet(1)/(60_pr*60_pr)
    Write(lproj,'(a37,x,f12.2,2x,f12.2,2x,f12.2)') '  Time in calculate_rotation_matrix: ',time_ibet(2),time_ibet(2)/60_pr,time_ibet(2)/(60_pr*60_pr)
    Write(lproj,'(a37,x,f12.2,2x,f12.2,2x,f12.2)') '  Time in calculate_overlaps:        ',time_ibet(3),time_ibet(3)/60_pr,time_ibet(3)/(60_pr*60_pr)
    Write(lproj,'(a37,x,f12.2,2x,f12.2,2x,f12.2)') '  Time in calculate_densities:       ',time_ibet(4),time_ibet(4)/60_pr,time_ibet(4)/(60_pr*60_pr)
    Write(lproj,'(a37,x,f12.2,2x,f12.2,2x,f12.2)') '  Time in calculate_energies:        ',time_ibet(5),time_ibet(5)/60_pr,time_ibet(5)/(60_pr*60_pr)
    Write(lproj,'(a37,x,f12.2,2x,f12.2,2x,f12.2)') '  Time in project:                   ',time_ibet(6),time_ibet(6)/60_pr,time_ibet(6)/(60_pr*60_pr)
    !
    ! --------------------
    ! Large scale printout
    ! --------------------
    large_printout=.False.
    If(large_printout) Then
#if(DO_PES==0)
       Write(enfname, '("energies",a4)') '.out'
#elif(DO_PES==1)
       Write(enfname, '("energies",a7,a4)') row_string,'.out'
#endif
       lenproj=9
       Open(lenproj,file=enfname,status='unknown')
       If(AMP_is_on .Eq. 1) Then
          jjmax=maxj
          If(maxj.Gt.40) jjmax=40
       Else
          jjmax=0
       End If
       Write(lenproj,'(f7.2,f7.2)') qmoment(2,3),qmoment(3,3)
       Do jj=0,jjmax,jjstep
          ekt(1)=projected_ekinN(jj); ekt(2)=projected_ekinP(jj); ekt(3)=ekt(1)+ekt(2)
          evol=projected_EVOL_rho_tau(jj)+projected_EVOL_rho_rho(jj); esurf=projected_ESURF_rho_drho(jj)
          eso=projected_ESO_rho_nablaj(jj); ecoul=projected_ecodi(jj)+projected_ecoex(jj); etens=projected_ETENS(jj)
          ept(1)=projected_eptN(jj); ept(2)=projected_eptP(jj); ept(3)=ept(1)+ept(2)
          eproj=ekt(3)+evol+esurf+eso+ecoul+etens+ept(3)
          Write(lenproj,'(f12.6,x)',Advance='No') eproj
       End Do
       Close(lenproj)
    End If
    !
  End Subroutine print_project
  !==========================================================
  ! Subroutine that finalizes symmetry restoration procedure.
  !==========================================================
  Subroutine finalize_projections()
    Implicit None
    !
    Close(lproj)
    Deallocate(betabs,betaweight)
    Deallocate(phiabs,ephi,ephic,ephicN,ephicP)
    Deallocate(rotated_overlap,detR)
    Deallocate(projected_overlap,projected_ekinN,projected_ekinP,projected_ecodi,projected_ecoex,projected_EVOL_rho_tau,&
               projected_EVOL_rho_rho,projected_ESURF_rho_drho,projected_ETENS,projected_ESO_rho_nablaj,projected_eptN,&
               projected_eptP,projected_xn1,projected_xn2,projected_rms1,projected_rms2,projected_delN,projected_delP,&
               projected_NP,projected_NP_norm)
    Deallocate(nz_sim,nr_sim,nl_sim)
    Deallocate(phicyl)
    Deallocate(ihil_convert,iphicyl_convert,ihil_iphicyl_convert)
    Deallocate(beta_active)
    Deallocate(iosc1_pair,iosc2_pair)
    Deallocate(itiphi_pair1,itiphi_pair2)
    Deallocate(ro_normalization,xl_ihil,xh_ihil)
    Deallocate(VmatrixN1,UmatrixN1,VmatrixN2,UmatrixN2,VmatrixP1,UmatrixP1,VmatrixP2,UmatrixP2)
    Deallocate(all_overlaps,all_energies)
    If(team_rank.Eq.0) Deallocate(all_overlaps_gthr,all_energies_gthr)
    !
  End Subroutine finalize_projections
  !===========================================================================
  ! Subroutine that calculates basis of eigenstates of the simplex-y operator.
  !===========================================================================
  Subroutine simplex_basis(lpr)
    Use math
    Implicit None
    !
    Logical lpr
    Integer(ipr) iosc
    !
    Do iosc=1,nt
       nz_sim(iosc)=nz(iosc)
       nr_sim(iosc)=nr(iosc)
       If(ns(iosc).Eq.1) Then
          nl_sim(iosc)=nl(iosc)
       Else
          nl_sim(iosc)=-nl(iosc)
       End If
    End Do
    !
    If(lpr) Then
       Do iosc=1,nt
          If(iosc.Eq.1) Then
          Write(lproj,'(a)')     '  ================================================================'
          Write(lproj,'(a)')     '  ======= Correspondence between HO and simplex-y bases  ========='
          Write(lproj,'(a)')     '  ================================================================'
          Write(lproj,'(2x,a,9x,a8,9x,a,9x,a15,5x,a)')   '|', 'HO basis', '|','Simplex-y basis','|'
          Write(lproj,'(2x,a)')   '----------------------------------------------------------'
          Write(lproj,'(4x,a,4x,a2,3x,a2,3x,a2,3x,a2,12x,a,3x,a2,3x,a2,3x,a2)') 'N','nz', &
                       'nr','nl','ns','N','nz','nr','nl'
          End If
          Write(lproj,'(i5,3x,i3,2x,i3,2x,i3,2x,i3,8x,i5,2x,i3,2x,i3,2x,i3)') iosc,nz(iosc),nr(iosc),nl(iosc), &
                       ns(iosc),iosc,nz_sim(iosc),nr_sim(iosc),nl_sim(iosc)
       End Do
    End If
    !
  End Subroutine simplex_basis
  !========================================================================================
  ! Subroutine that calculates the rotation matrix element between two spin-less HO states,
  ! based on Nazmiditinov et al. (1996).
  !========================================================================================
  Subroutine calculate_ry(bb,bp1,bz1,bp2,bz2,nz1,nr1,nl1,nz2,nr2,nl2,rotel)
    Use math
    Implicit None
    !
    Real(pr):: bb,bp1,bz1,bp2,bz2,rotel
    Integer(ipr):: nz1,nr1,nl1,nz2,nr2,nl2
    Real(pr):: eta,mub,muc,omega,deltap,deltam,M(1:3,1:3),Minv(1:3,1:3),detm
    Real(pr):: fmat11,fmat13,fmat22,fmat33,ffmat11,ffmat22,ffmat33,ffmat12,ffmat13,ffmat23,gmat11,gmat13,gmat22,gmat33, &
              ggmat11,ggmat22,ggmat33,ggmat12,ggmat13,ggmat23,kmat11,kmat13,kmat22,kmat31,kmat33,kkmat11,kkmat12,kkmat13, &
              kkmat21,kkmat22,kkmat23,kkmat31,kkmat32,kkmat33
    Integer(ipr):: ifl,np1,nm1,nn1,np2,nm2,nn2,ip1,mm1,mm2,mm3,mm4,mm5,mm6,lim1,lim2,lim3,ii1,jj1,kk1,ll1
    Real(pr):: p1,prefac,bfac,bsum1,bsum2,bsum3,bsum4,cfac,csum1,csum2,csum3,csum4,afac,asum1,asum2,asum3,asum4,asum5,asum6
    !
    rotel=zero
    nn1=2*nr1+Abs(nl1)+nz1; nn2=2*nr2+Abs(nl2)+nz2
    If(nn1.Ne.nn2) Return
    !------------------------------
    ! Defining F, G, and K matrices
    !------------------------------
    eta=(bp2/bp1)**2; mub=(bp1/bz1)**2; muc=(bp2/bz2)**2; omega=mub*eta+muc/eta+(one+mub*muc)*(sin(bb))**2+(mub+muc)*(cos(bb))**2
    deltap=(one-mub*muc)*(sin(bb))**2+(mub-muc)*(cos(bb))**2; deltam=(one-mub*muc)*(sin(bb))**2-(mub-muc)*(cos(bb))**2
    !
    fmat11=one/omega*(mub*eta+deltam-muc/eta); fmat13=-one/omega*sqrt(mub)*(muc-one)*sin(two*bb)
    fmat22=(eta-one)/(eta+one); fmat33=one/omega*(mub*eta-deltam-muc/eta)
    !
    ffmat11=half*(fmat11-fmat22); ffmat22=half*(fmat11-fmat22); ffmat33=fmat33
    ffmat12=-half*(fmat11+fmat22); ffmat13=-sqrt(half)*fmat13; ffmat23=sqrt(half)*fmat13
    !
    gmat11=-one/omega*(mub*eta-deltap-muc/eta); gmat13=one/omega*sqrt(muc)*(mub-one)*sin(two*bb)
    gmat22=-fmat22; gmat33=-one/omega*(mub*eta+deltap-muc/eta)
    !
    ggmat11=half*(gmat11-gmat22); ggmat22=half*(gmat11-gmat22); ggmat33=gmat33
    ggmat12=-half*(gmat11+gmat22); ggmat13=-sqrt(half)*gmat13; ggmat23=sqrt(half)*gmat13
    !
    kmat11=two/omega*(mub*sqrt(eta)+muc/sqrt(eta))*cos(bb); kmat13=-two/omega*sqrt(mub)*(muc/sqrt(eta)+sqrt(eta))*sin(bb)
    kmat22=two*sqrt(eta)/(eta+one); kmat31=two/omega*sqrt(muc)*(one/sqrt(eta)+mub*sqrt(eta))*sin(bb)
    kmat33=two/omega*sqrt(mub*muc)*(sqrt(eta)+one/sqrt(eta))*cos(bb)
    !
    kkmat11=half*(kmat11+kmat22); kkmat12=half*(-kmat11+kmat22); kkmat13=-sqrt(half)*kmat13
    kkmat21=half*(-kmat11+kmat22); kkmat22=half*(kmat11+kmat22); kkmat23=sqrt(half)*kmat13
    kkmat31=-sqrt(half)*kmat31; kkmat32=sqrt(half)*kmat31; kkmat33=kmat33
    !
    !----------------------------------------
    ! Calculating determinant of the M matrix
    !----------------------------------------
    M=zero; Minv=zero;
    Minv(1,1)=one; Minv(2,2)=one; Minv(3,3)=one
    M(1,1)=half*(((bp1/bp2)*cos(bb))**2+((bp1/bz2)*sin(bb))**2+one)
    M(1,3)=half*(-(bz1/bp2)*(bp1/bp2)+(bz1/bz2)*(bp1/bz2))*cos(bb)*sin(bb)
    M(2,2)=half*((bp1/bp2)**2+one)
    M(3,1)=half*(-(bz1/bp2)*(bp1/bp2)+(bz1/bz2)*(bp1/bz2))*cos(bb)*sin(bb)
    M(3,3)=half*(((bz1/bz2)*cos(bb))**2+((bz1/bp2)*sin(bb))**2+one)
    !
    Call lingd(3,3,3,3,M,Minv,detm,ifl)
    !-----------------------------------------------------
    ! Establishing quantum numbers np1,nm1,np2,nm2,nn1,nn2
    !-----------------------------------------------------
    If(nl1.Ge.0) Then
       np1=nr1+nl1; nm1=nr1
    Else
       np1=nr1; nm1=nr1+Abs(nl1)
    End If
    !
    If(nl2.Ge.0) Then
       np2=nr2+nl2; nm2=nr2
    Else
       np2=nr2; nm2=nr2+Abs(nl2)
    End If
    !
    !----------------------
    ! Calculating prefactor
    !----------------------
    ip1=iv(np1+np2+nr1+nr2)
    p1=(bp1/bp2)*Sqrt(bz1/bz2)/Sqrt(detm*two**(nn1+nn2))*Sqrt(fak(np1)*fak(nm1)*fak(nz1)*fak(np2)*fak(nm2)*fak(nz2))
    prefac=ip1*p1
    !
    Do mm1=0,np1  ! sum over n1
    Do mm2=0,nm1  ! sum over n2
    Do mm3=0,nz1  ! sum over n3
       bfac=zero
       If (Mod(np1+nm1+nz1-mm1-mm2-mm3,2).Eq.0) Then ! N-N_tilde is even
          lim1=(np1+nm1+nz1-mm1-mm2-mm3)/2-np1+mm1
          lim2=(np1+nm1+nz1-mm1-mm2-mm3)/2-nm1+mm2
          lim3=(np1+nm1+nz1-mm1-mm2-mm3)/2-nz1+mm3
          Do ii1=0,(lim2+lim3)/2
          Do jj1=0,(lim1+lim3)/2
          Do kk1=0,(lim1+lim2)/2
             If(lim3-ii1-jj1+kk1.Ge.0 .And. lim2-ii1-kk1+jj1.Ge.0 .And. lim1-jj1-kk1+ii1.Ge.0) Then
                bsum1=(ffmat11**ii1)*(ffmat22**jj1)*(ffmat33**kk1)*fi(ii1)*fi(jj1)*fi(kk1)*two**(-ii1-jj1-kk1) 
                bsum2=ffmat12**(lim3-ii1-jj1+kk1)*fi(lim3-ii1-jj1+kk1)
                bsum3=ffmat13**(lim2-ii1-kk1+jj1)*fi(lim2-ii1-kk1+jj1)
                bsum4=ffmat23**(lim1-jj1-kk1+ii1)*fi(lim1-jj1-kk1+ii1)
                bfac=bfac+bsum1*bsum2*bsum3*bsum4
             End If
          End Do ! kk1
          End Do ! jj1
          End Do ! ii1
          bfac=bfac*two**((np1+nm1+nz1-mm1-mm2-mm3)/2)
       End If ! Mod(np1+nm1+nz1-mm1-mm2-mm3,2).Eq.0
       !
       Do mm4=0,np2  ! sum over n1'
       Do mm5=0,nm2  ! sum over n2'
       Do mm6=0,nz2 ! sum over n3'
          cfac=zero;afac=zero
          If(Mod(np2+nm2+nz2-mm4-mm5-mm6,2).Eq.0) Then
          If(mm1+mm2+mm3.Eq.mm4+mm5+mm6) Then ! Condition on A, that N=N'. Otherwise, A component will be zero, and the entire term as well.
             lim1=(np2+nm2+nz2-mm4-mm5-mm6)/2-np2+mm4
             lim2=(np2+nm2+nz2-mm4-mm5-mm6)/2-nm2+mm5
             lim3=(np2+nm2+nz2-mm4-mm5-mm6)/2-nz2+mm6
             Do ii1=0,(lim2+lim3)/2
             Do jj1=0,(lim1+lim3)/2
             Do kk1=0,(lim1+lim2)/2
                If(lim3-ii1-jj1+kk1.Ge.0 .And. lim2-ii1-kk1+jj1.Ge.0 .And. lim1-jj1-kk1+ii1.Ge.0) Then
                   csum1=(ggmat11**ii1)*(ggmat22**jj1)*(ggmat33**kk1)*fi(ii1)*fi(jj1)*fi(kk1)*two**(-ii1-jj1-kk1)
                   csum2=ggmat12**(lim3-ii1-jj1+kk1)*fi(lim3-ii1-jj1+kk1)
                   csum3=ggmat13**(lim2-ii1-kk1+jj1)*fi(lim2-ii1-kk1+jj1)
                   csum4=ggmat23**(lim1-jj1-kk1+ii1)*fi(lim1-jj1-kk1+ii1)
                   cfac=cfac+csum1*csum2*csum3*csum4
                End If
             End Do
             End Do
             End Do
             cfac=cfac*two**((np2+nm2+nz2-mm4-mm5-mm6)/2)
             !
             Do ii1=0,mm1
             Do jj1=0,mm2
             Do kk1=0,mm1
             Do ll1=0,mm2
                If(mm4-ii1-jj1.Ge.0 .And. mm5-kk1-ll1.Ge.0 .And. mm1-ii1-kk1.Ge.0 .And. mm2-jj1-ll1.Ge.0 .And. mm6-mm1-mm2+ii1+jj1+kk1+ll1.Ge.0) Then
                   asum1=(kkmat11**ii1)*(kkmat12**jj1)*(kkmat21**kk1)*(kkmat22**ll1)*fi(ii1)*fi(jj1)*fi(kk1)*fi(ll1)
                   asum2=kkmat13**(mm4-ii1-jj1)*fi(mm4-ii1-jj1)
                   asum3=kkmat23**(mm5-kk1-ll1)*fi(mm5-kk1-ll1)
                   asum4=kkmat31**(mm1-ii1-kk1)*fi(mm1-ii1-kk1)
                   asum5=kkmat32**(mm2-jj1-ll1)*fi(mm2-jj1-ll1)
                   asum6=kkmat33**(mm6-mm1-mm2+ii1+jj1+kk1+ll1)*fi(mm6-mm1-mm2+ii1+jj1+kk1+ll1)
                   afac=afac+asum1*asum2*asum3*asum4*asum5*asum6
                End If
             End Do ! ll1
             End Do ! kk1
             End Do ! jj1
             End Do ! ii1
             afac=afac*two**(mm1+mm2+mm3)
          End If ! N=N'
          End If ! Mod(np2+nm2+nz2-mm4-mm5-mm6,2).Eq.0
          !
          rotel=rotel+bfac*afac*cfac
       End Do ! mm6
       End Do ! mm5
       End Do ! mm4
       !
    End Do ! mm3
    End Do ! mm2
    End Do ! mm1
    !
    rotel=prefac*rotel
  !
  End Subroutine calculate_ry
  !======================================================================================
  ! Function thats calculates small Wigner function for an arbitrary value of angle beta
  !======================================================================================
  Function wigner(j,mm1,mm2,beta) 
    Implicit None
    Integer(ipr) j,mm1,mm2,m,imax,i,i1,i2,i3
    Real(pr) wigner,beta,fac1,fac2,sum1
    !
    m=max(mm1,mm2)
    imax=j+mm2
    !
    sum1=zero
    Do i=0,imax
       i1=j+mm2-i
       i2=mm1-mm2+i
       i3=j-mm1-i
       If(i1.Ge.0 .And. i2.Ge.0 .And. i3.Ge.0) Then
          fac1=fi(i)*wf(j+mm1)*fi(i1)*wf(j-mm1); fac2=fi(i2)*wf(j+mm2)*fi(i3)*wf(j-mm2)
          sum1=sum1+iv(i1)*fac1*fac2*(Cos(half*beta))**(2*j+mm2-mm1-2*i)*(Sin(half*beta))**(mm1-mm2+2*i)
       End If
    End Do
    wigner=sum1*iv(j)
    !
  End Function wigner
  !=============================================================================
  ! Subroutine that performs Simpson's rule integration of a function
  ! defined by an array of equidistant values.
  ! Parameters are:
  ! > f --- array of values of the function
  ! > n --- number of points
  ! > h --- uniform spacing between values
  ! > r --- estimate of the integral that is returned to caller
  !=============================================================================
  Subroutine integrate_simpson(f,n,h,r)
    Implicit None
    Integer(ipr) n,npanels,nbegin,nend,nn
    Real(pr) h,r,x,f(n)
    Real(pr), Parameter:: c3d8=0.3750_pr
    !
    ! No integration for n=1
    If(n.Eq.1) Then
       r=f(1)
       Return
    End If
    !
    npanels=n-1
    nbegin=1
    r=zero
    If(Mod(npanels,2).Ne.0) Then
    ! Number of panels is odd. Use 3/8 rule on the first three panels.
       r=h*c3d8*(f(1)+3_pr*(f(2)+f(3))+f(4))
       If(n.Eq.3) Return
       nbegin=4
    Endif
    ! Apply 1/3 rule for the remaining panels.
    r=r+h*third*(f(nbegin)+4_pr*f(nbegin+1)+f(n))
    nbegin=nbegin+2
    If(nbegin.Eq.n) Return
    x=zero
    nend=n-2
    Do nn=nbegin,nend,2
       x=x+f(nn)+two*f(nn+1)
    End Do
    r=r+h*two*third*x
    !
  End Subroutine integrate_simpson
  !=============================================================================
  ! Subroutine that performs trapezoidal rule integration of a function
  ! defined by an array of values that are not necessarily equidistant.
  ! Parameters are:
  ! > f(n) --- array of values of the function
  ! > x(n) --- array of points in which function is evaluated
  ! > n --- number of points
  ! > r --- estimate of the integral that is returned to caller
  !=============================================================================
  Subroutine integrate_trapez(x,f,n,r)
    Implicit None
    Integer(ipr) n,k
    Real(pr) r,x(n),f(n),xstep
    !
    r=zero
    Do k=1,N-1
       xstep=x(k+1)-x(k)
       r=r+(f(k)+f(k+1))*half*xstep
    End Do
    !
  End Subroutine integrate_trapez
  !=============================================================================
  ! >
  ! >
  !=============================================================================
  !
End Module HFBTHO_projections
