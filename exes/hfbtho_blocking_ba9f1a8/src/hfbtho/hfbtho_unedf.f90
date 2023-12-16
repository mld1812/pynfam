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
!                             UNEDF PACKAGE                            !
!                                                                      !
! ==================================================================== !

!----------------------------------------------------------------------
!> This module provides the main \theCode DFT solver. It includes
!> routines for the calculation and diagonalization of the HFB matrix;
!> the definition of densities on the quadrature mesh and in configuration
!> space; the self-consistent loop the calculation of expectation values of
!> observables. It also includes a package of a few routines to perform
!> particle number projection.
!>
!>  @author
!>    Markus Kortelainen, Mario Stoitsov, Nicolas Schunck
!----------------------------------------------------------------------
Module UNEDF
  Use HFBTHO_utilities
  Implicit None
  !
  Character(16), Private :: Version='17'
  !
  ! Version History
  !--------------------------------------------------------------------------------------
  ! ver#17:(Mario)   use_TMR_pairing=0/1 standard/TMR pairing added
  !                  to Namelist. Using:
  !                  CpV0(0)=G,    CpV0(1)=a
  !                  CpV1(0)=vfacn,CpV1(1)=vfacp
  ! ver#16:(Mario)   #ifndef hide_dme preprocessing directive included
  ! ver#15:(Markus)  Added parameter CExPar, used in Coul. excange term.
  !                  Also, all the channels included in direct Hartree
  ! ver#14:(Markus)  Added function Vexternal for the external field,
  !                  and use_j2terms to switch off tensor terms.
  !                  Direct Hartree set to zero.
  ! Ver#13:(Mario)   Added ac2,ac3,acoord
  ! ver#12:(Mario)   hartree term temprorary dropped. rDr NNN terms taken
  !                  with a factor of 1/2
  ! ver#11:(Mario)   Gaussian approximation to the Hartree term added,
  ! [3/10/2010]      hatree_rc removed. NB! Function HartreeV is an
  !                  elemental function with possible array arguments
  ! ver#10: (Markus) Added e2charg (e^2 for Coulomb) to the public variables
  ! ver#9: (Mario)   Hartree 'CHrho' calculated in INM with rc='hatree_rc'
  ! [2/2/2010]       is subtracted from Crho(0)at DMEorder >= 0.
  !                  CHrho added to the public list, 'hatree_rc' added
  !                  to interaction parameters and the namelist.
  !                  In the case DMEorder=-1 (standard Skyrme)
  !                  both, 'CHrho' and 'hatree_rc', do not play.
  !                  New function HartreeV(u) defines Hatree energy as
  !                  E(Hartree)=(1/2)*Int[rho_0(r)*V(|r-r'|)*rho_0(r')]
  !                  HartreeV(u) is zero for u=<'hatree_rc'
  ! ver#8: (Markus)  Hartree DME terms dropped out.
  ! ver#7: (Markus)  Added switch to turn off the 3N terms.
  !        (Mario)   Added Abs to density and gradient dependent LDA
  !                  Public :: DMEorder,DMElda,use_DME3N_terms
  ! ver#6: (Mario)   Skyrme transformation added.
  ! ver#5: (Mario)   Print_Namelist=T/F added to the namelist
  ! ver#4: (Markus)  Added natural units to the module. Used only for printing.
  ! ver#3: (Mario)   Uamplitudes(0:3,0:7) in normal order
  !
  ! t for Uamplitudes(t,*)
  ! 0 -> 0,0
  ! 1 -> 1,1
  ! 2 -> 0,1
  ! 3 -> 1,0
  ! n for Uamplitudes(*,n)
  ! 0 -> U
  ! 1 -> dU/dRHO_0
  ! 2 -> dU/dRHO_1
  ! 3 -> d2U/(dRHO_0*dRHO_0)
  ! 4 -> d2U/(dRHO_1*dRHO_1)
  ! 5 -> d2U/(dRHO_0*dRHO_1)
  ! 6 -> dU/d(TAU_0)
  ! 7 -> dU/d(Delta RHO_0)
  !
  ! TESTED MATTHEMETICA<=>BIRUC & SCOTT; MATTHEMETICA<=>Module UNEDF (energy amplitudes only)
  ! ver#2: (Mario) Pairing included
  !  - set_functional_parameters(fname,lpr)
  !  - pairing incorporated into CpV0(0:1),CpV1(0:1)
  !    as public variables also serving two public amplitudes
  !     Urhorhopr(0:1,0)=CpV0(0:1)+CpV1(0:1)*rho(0)
  !     Urhorhopr(0:1,1)=CpV1(0:1)
  !     so, they can be used with appropriate values by the DME solver
  !  -need improvement later,
  !      currently HFBTHO uses CpV0(0:1), CpV0(0:1)  as before
  !      just substituting V0,V1 in pn-representation
  !      CpV0*(1-CpV1/0.16*rho_0)and this defines
  !      the default values in the module CpV0=V0,CpV1=1/2)
  !  -NAMELIST and input/output modified. RESERVED NAMES ARE:
  !      -namelist forbiden:
  !          'UNRDF'  - best UNEDF
  !          'SKYRME' - best SKYRME
  !      -namelist inforced but not for C-parameters (use_INM=F)
  !       or NM-parameters (use_INM=T) defined by the solver
  !          'FITS'
  !      -namelist inforced (one can overwrite all):
  !          'ANY OTHER NAME'
  !       i.e., the solver defines C-/NM- only using 'FITS'
  ! ver#1: (Mario) Complete rewrite consistent with HFBTHO
  !  -CB-LDA added
  !  -INM added
  !  -HFBTHO BENCHMARK: LN, ZR(110) prolate solution with SLY4,
  !   mixed pairing and tensor terms. Agreement with previouse
  !   implemetation to the last significant digit in the cases:
  !      - Standard Skyrme
  !      - LO+LDA
  !      - LO+CB-LDA
  !      - (NrNr=0,rDj=0), (rDr=0,jDr=0), 0.5(NrNr=-rDr,jDr=-rDj)
  !   -use_j2terms removed, i.e., in the SKYRME case CJ=0 removes all
  !    tensor terms, while in DME tensor terms are always present
  ! ver#0: (Marcus) Basic coding from scratch
  !   -DME(u) consistent with Mathematica numbers
  !   -including small 'u' approximation
  !--------------------------------------------------------------------------------------
  !
  ! === PUBLIC VARIABLES ===
  !
  ! Use pointers to prevent conflicts with UNEDF public variabes
  ! Example: Use UNEDF, pr=>my_pr, ipr=>my_ipr, Crho=>my_Crho ...
  !
  !--------------------------------------------------------------------------------------
  !
  ! === PUBLIC VARIABLES ===
  !
  Logical, Public :: use_charge_density, use_cm_cor,use_DME3N_terms,   &
                     use_j2terms,use_full_cm_cor,use_INM,use_Namelist, &
                     Print_Namelist,finite_range,hb0_charge_dependent, &
                     force_is_DME,use_3N_couplings,override_3N_couplings, &
                     coulomb_gaussian
  Integer(ipr), Public :: DMEorder,DMElda,use_TMR_pairing
  Real(pr), Public, Dimension(0:3,0:7) :: Urhorho,Urhotau,UrhoDrho,Unablarho  ! ph DME amplitudes
  Real(pr), Public, Dimension(0:3,0:7) :: UJnablarho,UrhonablaJ,UJJ,UJabJba
  Real(pr), Public, Dimension(0:3,0:7) :: Urhorhopr                           ! pp amplitudes
  Real(pr), Public, Dimension(0:1) :: UEnonstdr,UFnonstdr,URnonstdr           ! Other amplitudes
  Real(pr), Public :: hbzero,sigma,e2charg,CExPar                             ! hbr^2/2m, DD sigma, e^2 charge, coul.exch.
  Real(pr), Public :: hbzeron,hbzerop                                         ! hbr^2/2m_n,hbr^2/2m_p
  Real(pr), Public, Dimension(0:1) :: Crho,Cdrho,Ctau,CrDr,CrdJ,CJ,CpV0,CpV1  ! basic coupling constants
  Real(pr), Public, Dimension(0:1) :: Cnrho,CJdr
  Real(pr), Public :: E_NM,K_NM,SMASS_NM,RHO_NM,ASS_NM,LASS_NM,VMASS_NM,P_NM,KA_NM
  Real(pr), Public :: CHrho                                                   ! Crho(0) from the Hartree term in NM
  Real(pr), Public :: mpi,gA,fpi,c1,c3,c4,cd,ce,LambdaX
  Real(pr), PUBLIC :: t0s,t0a,drs,dra,ts,ta,t3alp,t3al0,t3alm,t324,alp,alm,wla0, &
                      wla1,TA7,TA8,TB7,TB8,tv1,tv2,tv3,tv4,tv5,tv6,ts1,ts2,t4o3
  Real(pr), PUBLIC :: t0_pub,t1_pub,t2_pub,t3_pub,x0_pub,x1_pub,x2_pub,x3_pub,&
                      b4_pub,b4p_pub,te_pub,to_pub
  !Gogny parameters
  real(pr), Public, allocatable, dimension(:) :: mu_g, W_g, B_g, H_g, M_g, mu_g_all
  integer(ipr) :: n_g=0, n_g_all=0
  !Coulomb by gaussians parameters
  integer(ipr), Public, parameter :: n_g_coul = 5
  real(pr),  dimension(1:n_g_coul) :: &
       mu_g_coul = [53.56502433354500,10.60650459879994,4.35698887528154,2.08951260631933,0.79888844310443],&
       V_g_coul  = [ 0.05331491213328, 0.11653242893617,0.19644220323999,0.41149438105314,1.91615663177677]
  !Bulgac functionals
  Logical, Public :: is_NEDF
  Real(pr), Public, Dimension(0:2) :: a_NEDF,b_NEDF,c_NEDF
  Real(pr), Public :: eta_NEDF , W0_NEDF !
  !
  ! DME variables
  integer(ipr) :: ndme
  real(pr), allocatable, dimension(:,:) :: grhorho,grhotau,grhodelrho,&
       gJabJab,gJabJba,h_rho_rho_rho,h_rho_rho_tau,h_rho_rho_delrho,&
       h_rho_nablarho,h_rho_J_nablarho,h_rho_Jab_Jab,h_rho_Jab_Jba

  real(pr) :: alfa_dme
  logical :: g_couplings_allocated = .false.
  logical :: h_couplings_allocated = .false.
  ! === PRIVATE VARIABLES ===
  !
  Logical, Private :: coulomb_regulator = .False.
  Real(pr), Private, Dimension(0:1) :: nuCrho,nuCdrho,nuCtau,nuCrDr  ! basic coupling constants in natural units
  Real(pr), Private, Dimension(0:1) :: nuCrdJ,nuCJ,nuCpV0,nuCpV1     !
  Real(pr), Private :: t0,t1,t2,t3,x0,x1,x2,x3,b4,b4p,te,to
  Real(pr), Private :: nuLambda,nufpi                                ! parameters associated to natural units
  Integer(ipr), Private :: i_cut                                     ! dmeorder: -1=Standard Skyrme, 0=LO, 1=NLO, 2=N2LO
  Real(pr), Private :: Pi,eps                                        ! dmelda: 0=Kf-LDA, 1=CB-LDA
  Real(pr), Private :: kfconst,CK                                    ! (3Pi^2/2)^(1/3)
  Real(pr), Parameter, Private :: mevfm=197.30_pr;
  Real(pr), Private :: rho(0:1),tau(0:1),nrho2(0:1),lrho(0:1)
  Real(pr), Private :: mpi2,fpi2,fpi4,gA2,gA4,gA6,CHartree
  Real(pr), Private :: arhorho,brhorho,arhodrho,brhodrho,arhotau,brhotau,ajj,bjj,adrdr,bdrdr
  Real(pr), Private :: darhorho,dbrhorho,darhodrho,dbrhodrho,darhotau,dbrhotau,dajj,dbjj,dadrdr,dbdrdr
  Real(pr), Private :: ddarhodrho,ddbrhodrho,ddarhotau,ddbrhotau,ddarhorho,ddbrhorho
  Real(pr), Private :: hrho0rho0,hrho1rho1,hdr0dr0,hdr1dr1,hrho0Drho0,hrho1Drho0, &
                       hrho1Drho1,hrho0tau0,hrho1tau0,hrho1tau1,hJ0dr0,hrho0DJ0,hJ1dr1,hrho1DJ1, &
                       hJ0dr1,hrho1DJ0,hJ1dr0,hJ0J0,hJ0J1,hJ1J1
  Real(pr), Private :: dhrho0rho0,dhrho1rho1,dhdr0dr0,dhdr1dr1,dhrho0Drho0, &
                       dhrho1Drho0,dhrho1Drho1,dhrho0tau0,dhrho1tau0,dhrho1tau1,dhJ0dr0,dhrho0DJ0, &
                       dhJ1dr1,dhrho1DJ1,dhJ0dr1,dhrho1DJ0,dhJ1dr0,dhJ0J0,dhJ0J1,dhJ1J1
  Real(pr), Private :: ddhrho0rho0,ddhrho1rho1,ddhrho0Drho0,ddhrho1Drho0, &
                       ddhrho1Drho1,ddhrho0tau0,ddhrho1tau0,ddhrho1tau1
  Real(pr), Private, Dimension(3,3,33) :: ctr0r0,ctr1r1,ctdr0dr0,ctdr1dr1, & ! coefficients for 3N part
                                          ctr0Dr0,ctr1Dr0,ctr1Dr1,ctr0t0,ctr1t0,ctr1t1,ctJ0dr0,ctr0dJ0,ctJ1dr1, &
                                          ctr1dJ1,ctJ0dr1,ctr1dJ0,ctJ1dr0,ctJ0J0,ctJ0J1,ctJ1J1
  Real(pr), Private :: u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12
  Real(pr), Private :: ual,lual,atu,asqu,asqu4
  Real(pr), Private :: ac2,ac3,acoord
  Parameter(acoord=0.50_pr,ac2=4.0_pr*(acoord**2-acoord+0.50_pr),ac3=2.0_pr*(acoord**2-acoord+0.50_pr))
  Character (30) :: FunctionalName
  !
  Real(pr), Private :: A1_1,A1_2,A1_3,A1_4,A1_5,b1_1,b1_2,b1_3,b1_4,b1_5
  Real(pr), Private :: A3_1,A3_2,A3_3,A3_4,A3_5,b3_1,b3_2,b3_3,b3_4,b3_5
  Real(pr), Private :: h0mpi6,h0mpi6c1,h0mpi6c3,h0mpi6NM,h0mpi6c1NM,h0mpi6c3NM
  !
  Namelist /UNEDF_NAMELIST/ FunctionalName, DMEorder, DMElda, use_INM, hbzero, use_TMR_pairing,   &
                            Crho, Cdrho, Ctau, CrDr, CrdJ, CJ, sigma, CpV0, CpV1, e2charg,        &
                            E_NM, K_NM, SMASS_NM, RHO_NM, ASS_NM, LASS_NM, VMASS_NM,              &
                            mpi, gA, fpi, c1, c3, c4, cd, ce, LambdaX,                            &
                            use_cm_cor, use_charge_density, use_DME3N_terms, use_j2terms, CExPar, &
                            Print_Namelist
 Contains
  !
  !=======================================================================
  !
  !=======================================================================
  Subroutine calculate_U_parameters(rho0_in,rho1_in,tau0_in,tau1_in,laprho0,laprho1,nablarho0s,nablarho1s)
    Implicit None
    Real(pr), Intent(in) :: rho0_in,rho1_in,tau0_in,tau1_in
    Real(pr), Intent(in), Optional :: &
         nablarho0s,nablarho1s,laprho0,laprho1
    Integer(ipr) :: t,i,j,k,l
    Real(pr) :: u,du,ddu,dtu,dlu
    Real(pr) :: ph,aux,daux,ddaux
    Real(pr) :: y,dy,ddy,marc,dmarc,ddmarc,mlog,dmlog,ddmlog
    Real(pr) :: ucut,ucut3n
    !
    ucut=0.1_pr; ucut3n=0.6_pr
    !
    rho(0)=rho0_in; rho(1)=rho1_in;
    tau(0)=tau0_in; tau(1)=tau1_in;
    !
    lrho=0.0_pr; nrho2=0.0_pr;
    If (Present(laprho0)) lrho(0)=laprho0
    If (Present(laprho1)) lrho(1)=laprho1
    If (Present(nablarho0s)) nrho2(0)=nablarho0s
    If (Present(nablarho1s)) nrho2(1)=nablarho1s
    !
    arhorho=0.0_pr; darhorho=0.0_pr; ddarhorho=0.0_pr
    brhorho=0.0_pr; dbrhorho=0.0_pr; ddbrhorho=0.0_pr
    arhodrho=0.0_pr; darhodrho=0.0_pr; ddarhodrho=0.0_pr
    brhodrho=0.0_pr; dbrhodrho=0.0_pr; ddbrhodrho=0.0_pr
    arhotau=0.0_pr; darhotau=0.0_pr; ddarhotau=0.0_pr
    brhotau=0.0_pr; dbrhotau=0.0_pr; ddbrhotau=0.0_pr
    adrdr=0.0_pr; dadrdr=0.0_pr
    bdrdr=0.0_pr; dbdrdr=0.0_pr
    ajj=0.0_pr; dajj=0.0_pr
    bjj=0.0_pr; dbjj=0.0_pr
    !
    hrho0rho0=0.0_pr; hrho1rho1=0.0_pr; hdr0dr0=0.0_pr; hdr1dr1=0.0_pr
    hrho0Drho0=0.0_pr; hrho1Drho0=0.0_pr; hrho1Drho1=0.0_pr
    hrho0tau0=0.0_pr; hrho1tau0=0.0_pr; hrho1tau1=0.0_pr
    hJ0dr0=0.0_pr; hrho0DJ0=0.0_pr; hJ1dr1=0.0_pr; hrho1DJ1=0.0_pr
    hJ0dr1=0.0_pr; hrho1DJ0=0.0_pr; hJ1dr0=0.0_pr
    hJ0J0=0.0_pr; hJ0J1=0.0_pr; hJ1J1=0.0_pr
    dhrho0rho0=0.0_pr; dhrho1rho1=0.0_pr; dhdr0dr0=0.0_pr; dhdr1dr1=0.0_pr
    dhrho0Drho0=0.0_pr; dhrho1Drho0=0.0_pr; dhrho1Drho1=0.0_pr
    dhrho0tau0=0.0_pr; dhrho1tau0=0.0_pr; dhrho1tau1=0.0_pr
    dhJ0dr0=0.0_pr; dhrho0DJ0=0.0_pr; dhJ1dr1=0.0_pr; dhrho1DJ1=0.0_pr
    dhJ0dr1=0.0_pr; dhrho1DJ0=0.0_pr; dhJ1dr0=0.0_pr
    dhJ0J0=0.0_pr; dhJ0J1=0.0_pr; dhJ1J1=0.0_pr
    ddhrho0rho0=0.0_pr; ddhrho1rho1=0.0_pr
    ddhrho0Drho0=0.0_pr; ddhrho1Drho0=0.0_pr; ddhrho1Drho1=0.0_pr
    ddhrho0tau0=0.0_pr; ddhrho1tau0=0.0_pr; ddhrho1tau1=0.0_pr
    !
    u=0.0_pr; du=0.0_pr; ddu=0.0_pr; dtu=0.0_pr; dlu=0.0_pr
    !
    ! U and partial derivatives with respect of rho_0 in Thomas Fermi approximation
    !
    If (dmeorder.Ge.0) Then
       If (dmelda.Eq.0) Then
          ! density dependent LDA
          u=(kfconst/mpi)*rho(0)**(1.0_pr/3.0_pr)
          du=(1.0_pr/3.0_pr)*u/(rho(0)+eps)                 ! u'(RHO_0)
          ddu=-(2.0_pr/9.0_pr)*u/(rho(0)**2+eps)            ! u''(RHO_0)
          dtu=0.0_pr                                        ! u'(TAU_0)
          dlu=0.0_pr                                        ! u'(DeltaRHO_0)
       Else
          ! density and gradient dependent LDA
          u=Sqrt(Abs((5.0_pr/3.0_pr)*(tau(0)-0.250_pr*lrho(0))/(rho(0)+eps)))/mpi
          du=-0.50_pr*u/(rho(0)+eps)                        ! u'(RHO_0)
          ddu=0.750_pr*u/(rho(0)**2+eps)                    ! u''(RHO_0)
          dtu=0.50_pr*u/(Abs(tau(0)-0.250_pr*lrho(0))+eps)  ! u'(TAU_0)
          dlu=-0.250_pr*dtu                                 ! u'(DeltaRHO_0)
       End If
    End If
    !
    ! Partial optimization
    u2=u*u; u3=u2*u; u4=u3*u; u5=u4*u; u6=u5*u; u7=u6*u; u8=u7*u; u9=u8*u; u10=u9*u; u11=u10*u; u12=u11*u;
    ual=1.0_pr+4.0_pr*u2; lual=Log(ual); atu=Atan(2.0_pr*u); asqu=Sqrt(1.0_pr+u2); asqu4=Sqrt(4.0_pr+u2)
    !
    ! A and B functions and their partial derivatives with respect of u
    !
    !  LO, 2N terms
    !
    If (dmeorder.Ge.0) Then
       If (u.Gt.ucut) Then
          arhorho=(gA2*(4.0_pr*u2*(21.0_pr-498.0_pr*u2-64.0_pr*u4+16.0_pr*u6)+48.0_pr*u3 &
               *(35.0_pr+4.0_pr*u2)*atu+3.0_pr*(-7.0_pr+16.0_pr*u2*(-8.0_pr+9.0_pr*u2))*lual)) &
               /(1024.0_pr*fpi2*u8)
          darhorho=(gA2*(4.0_pr*u2*(-21.0_pr+279.0_pr*u2+16.0_pr*u4)-6.0_pr*u3*(175.0_pr &
               +12.0_pr*u2)*atu+3.0_pr*(7.0_pr+96.0_pr*u2-72.0_pr*u4)*lual))/(128.0_pr*fpi2*u9)
          ddarhorho=(-3.0_pr*gA2*(4.0_pr*u2*(-63.0_pr+504.0_pr*u2+16.0_pr*u4)-12.0_pr*u3*(175.0_pr &
               +8.0_pr*u2)*atu+(63.0_pr+672.0_pr*u2-360.0_pr*u4)*lual))/(128.0_pr*fpi2*u10)
          brhorho=2.0_pr*arhorho
          dbrhorho=2.0_pr*darhorho
          ddbrhorho=2.0_pr*ddarhorho
          !
          arhodrho=(35.0_pr*gA2*(-4.0_pr*u2*(-3.0_pr+72.0_pr*u2+4.0_pr*u4-60.0_pr*u*atu) &
               +3.0_pr*(-1.0_pr-18.0_pr*u2+24.0_pr*u4)*lual))/(12288.0_pr*fpi2*mpi2*u10)
          darhodrho=(35.0_pr*gA2*(4.0_pr*u2*(-15.0_pr+234.0_pr*u2+8.0_pr*u4 &
               -210.0_pr*u*atu)+3.0_pr*(5.0_pr+72.0_pr*u2-72.0_pr*u4)*lual)) &
               /(6144.0_pr*fpi2*mpi2*u11)
          ddarhodrho=(35.0_pr*gA2*(4.0_pr*u2*(165.0_pr-1746.0_pr*u2-40.0_pr*u4+1680.0_pr*u &
               *atu)+3.0_pr*(-55.0_pr+72.0_pr*u2*(-9.0_pr+7.0_pr*u2))*lual))/(6144.0_pr*fpi2*mpi2*u12)
          !
          brhodrho=2.0_pr*arhodrho
          dbrhodrho=2.0_pr*darhodrho
          ddbrhodrho=2.0_pr*ddarhodrho
          !
          arhotau=(35.0_pr*gA2*(4.0_pr*u2*(-3.0_pr+72.0_pr*u2+4.0_pr*u4)-240.0_pr*u3 &
               *atu+(3.0_pr+54.0_pr*u2-72.0_pr*u4)*lual))/(3072.0_pr*fpi2*mpi2*u10)
          darhotau=(35.0_pr*gA2*(4.0_pr*u2*(15.0_pr-234.0_pr*u2-8.0_pr*u4+210.0_pr*u &
               *atu)+3.0_pr*(-5.0_pr+72.0_pr*u2*(-1.0_pr+u2))*lual))/(1536.0_pr*fpi2*mpi2*u11)
          ddarhotau=(35.0_pr*gA2*(4.0_pr*u2*(-165.0_pr+1746.0_pr*u2+40.0_pr*u4-1680.0_pr*u*atu) &
               +3.0_pr*(55.0_pr+72.0_pr*u2*(9.0_pr-7.0_pr*u2))*lual))/(1536.0_pr*fpi2*mpi2*u12)
          brhotau=2.0_pr*arhotau
          dbrhotau=2.0_pr*darhotau
          ddbrhotau=2.0_pr*ddarhotau
          !
          ajj=(3.0_pr*gA2*(-4.0_pr*u2+8.0_pr*u4+lual))/(512.0_pr*fpi2*mpi2*u6)
          dajj=(-3.0_pr*gA2*(4.0_pr*u2*(-3.0_pr-6.0_pr*u2+8.0_pr*u4)+3.0_pr*(ual)*lual)) &
               /(256.0_pr*fpi2*mpi2*u7*(ual))
          bjj=2.0_pr*ajj
          dbjj=2.0_pr*dajj
       Else
          Arhorho=(gA2*u4*(121.0_pr-448.0_pr*u2))/(1155.0_pr*fpi2)
          dArhorho=(4.0_pr*gA2*u3*(121.0_pr-672.0_pr*u2))/(1155.0_pr*fpi2)
          ddArhorho=(4.0_pr*gA2*u2*(1573.0_pr+560.0_pr*u2*(-26.0_pr+153.0_pr*u2)))/(5005.0_pr*fpi2)
          Brhorho=2.0_pr*Arhorho
          dBrhorho=2.0_pr*dArhorho
          ddBrhorho=2.0_pr*ddArhorho
          !
          Arhotau=(gA2*(1287.0_pr-4004.0_pr*u2+11232.0_pr*u4-31680.0_pr*u6))/(10296.0_pr*fpi2*mpi2)
          dArhotau=(gA2*u*(-1001.0_pr+432.0_pr*u2*(13.0_pr-55.0_pr*u2)))/(1287.0_pr*fpi2*mpi2)
          ddArhotau=(gA2*(-1001.0_pr+16.0_pr*u2*(1053.0_pr-7425.0_pr*u2+40040.0_pr*u4)))/(1287.0_pr*fpi2*mpi2)
          Brhotau=2.0_pr*Arhotau
          dBrhotau=2.0_pr*dArhotau
          ddBrhotau=2.0_pr*ddArhotau
          !
          ArhoDrho=-Arhotau/4.0_pr
          dArhoDrho=-dArhotau/4.0_pr
          ddArhoDrho=-ddArhotau/4.0_pr
          BrhoDrho=2.0_pr*ArhoDrho
          dBrhoDrho=2.0_pr*dArhoDrho
          ddBrhoDrho=2.0_pr*ddArhoDrho
          !
          AJJ=(gA2*(5.0_pr-15.0_pr*u2+48.0_pr*u4-160.0_pr*u6))/(40.0_pr*fpi2*mpi2)
          dAJJ=(-3.0_pr*gA2*u*(5.0_pr-32.0_pr*u2+160.0_pr*u4))/(20.0_pr*fpi2*mpi2)
          Bjj=2.0_pr*Ajj;
          dBjj=2.0_pr*dAjj
       End If
    End If
    !
    ! NLO, 2N terms
    !
    If (dmeorder.Ge.1) Then
       If (u.Gt.ucut) Then
          arhorho=arhorho+(mpi2*(4.0_pr*u2*(11025.0_pr*(-1.0_pr + 10.0_pr*gA2 + 127.0_pr*gA4) &
               - 525.0_pr*(-839.0_pr + 3014.0_pr*gA2 + 50489.0_pr*gA4)*u2 - 4200.0_pr*(-17.0_pr &
               - 6.0_pr*gA2 + 495.0_pr*gA4)*u4 - 140.0_pr*(-133.0_pr - 718.0_pr*gA2 &
               + 107.0_pr*gA4)*u6 + 1536.0_pr*(1.0_pr + 10.0_pr*gA2 + 13.0_pr*gA4)*u8) &
               - 525.0_pr*Log(1.0_pr + 2.0_pr*u*(u + asqu))*(4.0_pr*u*asqu*(21.0_pr*(-1.0_pr &
               + 10.0_pr*gA2 + 127.0_pr*gA4) - 2.0_pr*(-167.0_pr + 518.0_pr*gA2 + 9305.0_pr*gA4)*u2 &
               + 8.0_pr*(7.0_pr + 10.0_pr*gA2 - 153.0_pr*gA4)*u4 + 16.0_pr*(1.0_pr + 6.0_pr*gA2 &
               + gA4)*u6) + 3.0_pr*(7.0_pr - 70.0_pr*gA2 - 889.0_pr*gA4 - 64.0_pr*(-1.0_pr + 6.0_pr*gA2 &
               + 83.0_pr*gA4)*u2 + 48.0_pr*(-1.0_pr + 2.0_pr*gA2 + 47.0_pr*gA4)*u4)*Log(1.0_pr &
               + 2.0_pr*u*(u + asqu)))))/(1.72032e7_pr*fpi4*Pi**2*u8)

          darhorho=darhorho+(mpi2*(4.0_pr*u2*(-3675.0_pr*(-1.0_pr+10.0_pr*gA2+127.0_pr*gA4) &
               +175.0_pr*(-503.0_pr+1766.0_pr*gA2+29897.0_pr*gA4)*u2+2800.0_pr*(-3.0_pr-2.0_pr*gA2 &
               +81.0_pr*gA4)*u4-700.0_pr*(1.0_pr+6.0_pr*gA2+gA4)*u6+128.0_pr*(1.0_pr+10.0_pr*gA2 &
               +13.0_pr*gA4)*u8)*(asqu+2.0_pr*u*(1.0_pr+u*(u+asqu))) &
               -700.0_pr*u*(1.0_pr+u2)*(21.0_pr-2*u2*(95.0_pr+8.0_pr*u2)+gA2*(-210.0_pr+556.0_pr*u2 &
               -32.0_pr*u4)+gA4*(-2667.0_pr+10306.0_pr*u2+304.0_pr*u4))*(1.0_pr+2.0_pr*u*(u+Sqrt(1.0_pr &
               +u2)))*Log(1.0_pr+2.0_pr*u*(u+asqu))+525.0_pr*(7.0_pr-70.0_pr*gA2 &
               -889.0_pr*gA4-48.0_pr*(-1.0_pr+6.0_pr*gA2+83.0_pr*gA4)*u2+24.0_pr*(-1.0_pr+2.0_pr*gA2 &
               +47.0_pr*gA4)*u4)*(asqu+2.0_pr*u*(1.0_pr+u*(u+asqu)))*Log(1.0_pr &
               +2.0_pr*u*(u+asqu))**2))/(716800.0_pr*fpi4*Pi**2*u9*asqu*(1.0_pr &
               +2.0_pr*u*(u+asqu)))

          ddarhorho=ddarhorho+(mpi2*(4.0_pr*u2*(33075.0_pr*(-1.0_pr+10.0_pr*gA2+127.0_pr*gA4) &
               -525.0_pr*(-965.0_pr+3314.0_pr*gA2+56699.0_pr*gA4)*u2-2800.0_pr*(-11.0_pr-10.0_pr*gA2 &
               +281.0_pr*gA4)*u4+700.0_pr*(1.0_pr+6.0_pr*gA2+gA4)*u6+128.0_pr*(1.0_pr+10.0_pr*gA2 &
               +13.0_pr*gA4)*u8)*(asqu+4.0_pr*(u+u3)*(1.0_pr+2.0_pr*u*(u &
               +asqu)))+525.0_pr*Log(1.0_pr+2.0_pr*u*(u+asqu))*(4.0_pr*u*(1.0_pr &
               +u2)*(63.0_pr-346.0_pr*u2-16.0_pr*u4+gA2*(-630.0_pr+964.0_pr*u2-32.0_pr*u4) &
               +gA4*(-8001.0_pr+18406.0_pr*u2+304.0_pr*u4))*(1.0_pr+4.0_pr*u*(asqu &
               +2.0_pr*u*(1.0_pr+u*(u+asqu))))-3.0_pr*(-21.0_pr*(-1.0_pr+10.0_pr*gA2 &
               +127.0_pr*gA4)-112.0_pr*(-1.0_pr+6.0_pr*gA2+83.0_pr*gA4)*u2+40.0_pr*(-1.0_pr &
               +2.0_pr*gA2+47.0_pr*gA4)*u4)*(asqu+4.0_pr*(u+u3)*(1.0_pr+2.0_pr*u*(u &
               +asqu)))*Log(1.0_pr+2.0_pr*u*(u+asqu)))))/(716800.0_pr*fpi4*Pi**2*u10*asqu &
               *(1.0_pr+2.0_pr*u*(u+asqu))**2)

          brhorho=brhorho+(mpi2*(4.0_pr*u2*(11025.0_pr*(-1.0_pr + 10.0_pr*gA2 + 43.0_pr*gA4) &
               - 525.0_pr*(-839.0_pr + 3014.0_pr*gA2 + 12269.0_pr*gA4)*u2 - 4200.0_pr*(-17.0_pr &
               - 6.0_pr*gA2 + 123.0_pr*gA4)*u4 - 140.0_pr*(-133.0_pr - 718.0_pr*gA2 &
               + 1223.0_pr*gA4)*u6 - 1536.0_pr*(-1.0_pr - 10.0_pr*gA2 + 23.0_pr*gA4)*u8) &
               + 525.0_pr*Log(1.0_pr + 2.0_pr*u*(u + asqu))*(4.0_pr*u*asqu*(21.0_pr &
               - 2.0_pr*u2*(167.0_pr + 28.0_pr*u2 + 8.0_pr*u4) - 2.0_pr*gA2*(105.0_pr &
               - 518.0_pr*u2 + 40.0_pr*u4 + 48.0_pr*u6) + gA4*(-903.0_pr + 4378.0_pr*u2 &
               + 360.0_pr*u4 + 176.0_pr*u6)) - 3.0_pr*(7.0_pr - 70.0_pr*gA2 - 301.0_pr*gA4 &
               - 64.0_pr*(-1.0_pr + 6.0_pr*gA2 + 23.0_pr*gA4)*u2 + 48.0_pr*(-1.0_pr + 2.0_pr*gA2 &
               + 11.0_pr*gA4)*u4)*Log(1.0_pr + 2.0_pr*u*(u + asqu)))))/(8.6016e6_pr*fpi4*Pi**2*u8)

          dbrhorho=dbrhorho+(mpi2*(-4.0_pr*u2*(3675.0_pr*(-1.0_pr+10.0_pr*gA2+43.0_pr*gA4) &
               -175.0_pr*(-503.0_pr+1766.0_pr*gA2+7229.0_pr*gA4)*u2-2800.0_pr*(-3.0_pr-2.0_pr*gA2 &
               +21.0_pr*gA4)*u4-700.0_pr*(-1.0_pr-6.0_pr*gA2+11.0_pr*gA4)*u6+128.0_pr*(-1.0_pr &
               -10.0_pr*gA2+23.0_pr*gA4)*u8)*(asqu+2.0_pr*u*(1.0_pr+u*(u+Sqrt(1.0_pr &
               +u2))))-700.0_pr*u*(1.0_pr+u2)*(21.0_pr-2.0_pr*u2*(95.0_pr+8.0_pr*u2)+gA2*(-210.0_pr &
               +556.0_pr*u2-32.0_pr*u4)+gA4*(-903.0_pr+2410.0_pr*u2+112.0_pr*u4))*(1.0_pr+2.0_pr*u*(u &
               +asqu))*Log(1.0_pr+2.0_pr*u*(u+asqu))+525.0_pr*(7.0_pr-70.0_pr*gA2 &
               -301.0_pr*gA4-48.0_pr*(-1.0_pr+6.0_pr*gA2+23.0_pr*gA4)*u2+24.0_pr*(-1.0_pr+2.0_pr*gA2 &
               +11.0_pr*gA4)*u4)*(asqu+2.0_pr*u*(1.0_pr+u*(u+asqu)))*Log(1.0_pr &
               +2.0_pr*u*(u+asqu))**2))/(358400.0_pr*fpi4*Pi**2*u9*asqu*(1.0_pr &
               +2.0_pr*u*(u+asqu)))

          ddbrhorho=ddbrhorho+(mpi2*(-4.0_pr*u2*(-33075.0_pr*(-1.0_pr+10.0_pr*gA2+43.0_pr*gA4) &
               +525.0_pr*(-965.0_pr+3314.0_pr*gA2+13655.0_pr*gA4)*u2+2800.0_pr*(-11.0_pr-10.0_pr*gA2 &
               +77.0_pr*gA4)*u4+700.0_pr*(-1.0_pr-6.0_pr*gA2+11.0_pr*gA4)*u6+128.0_pr*(-1.0_pr &
               -10.0_pr*gA2+23.0_pr*gA4)*u8)*(asqu+4.0_pr*(u+u3)*(1.0_pr+2.0_pr*u*(u &
               +asqu)))+525.0_pr*Log(1.0_pr+2.0_pr*u*(u+asqu))*(4.0_pr*u*(1.0_pr &
               +u2)*(-63.0_pr*(-1.0_pr+10.0_pr*gA2+43.0_pr*gA4)+2.0_pr*(-173.0_pr+482.0_pr*gA2 &
               +2159.0_pr*gA4)*u2+16.0_pr*(-1.0_pr-2.0_pr*gA2+7.0_pr*gA4)*u4)*(1.0_pr &
               +4.0_pr*u*(asqu+2.0_pr*u*(1.0_pr+u*(u+asqu))))-3.0_pr*( &
               -21.0_pr*(-1.0_pr+10.0_pr*gA2+43.0_pr*gA4)-112.0_pr*(-1.0_pr+6.0_pr*gA2 &
               +23.0_pr*gA4)*u2+40.0_pr*(-1.0_pr+2.0_pr*gA2+11.0_pr*gA4)*u4)*(asqu &
               +4.0_pr*(u+u3)*(1.0_pr+2.0_pr*u*(u+asqu)))*Log(1.0_pr+2.0_pr*u*(u &
               +asqu)))))/(358400.0_pr*fpi4*Pi**2*u10*asqu*(1.0_pr &
               +2.0_pr*u*(u+asqu))**2)

          arhodrho=arhodrho -(-44100.0_pr*(-1.0_pr+10.0_pr*gA2+127.0_pr*gA4)*u2+4200.0_pr*(-421.0_pr &
               +1516.0_pr*gA2+25369.0_pr*gA4)*u4+420.0_pr*(-283.0_pr-30.0_pr*gA2+8641.0_pr*gA4)*u6 &
               +1680.0_pr*(-19.0_pr-106.0_pr*gA2+5.0_pr*gA4)*u8+464.0_pr*(1.0_pr+10.0_pr*gA2 &
               +13.0_pr*gA4)*u10-420.0_pr*u*asqu*(105.0_pr-1685.0_pr*u2-134.0_pr*u4 &
               -8.0_pr*u6+16.0_pr*u8+10.0_pr*gA2*(-105.0_pr+527.0_pr*u2-22.0_pr*u4-8.0_pr*u6 &
               +16.0_pr*u8)+gA4*(-13335.0_pr+94295.0_pr*u2+2738.0_pr*u4-104.0_pr*u6+208.0_pr*u8)) &
               *Log(1.0_pr+2.0_pr*u*(u+asqu))+11025.0_pr*(1.0_pr-10.0_pr*gA2-127.0_pr*gA4 &
               -9.0_pr*(-1.0_pr+6.0_pr*gA2+83.0_pr*gA4)*u2+8.0_pr*(-1.0_pr+2.0_pr*gA2 &
               +47.0_pr*gA4)*u4)*Log(1.0_pr+2.0_pr*u*(u+asqu))**2-128.0_pr*u3*(33600.0_pr &
               *(-1.0_pr+4.0_pr*gA2+16.0_pr*gA4)*u+560.0_pr*(-29.0_pr-10.0_pr*gA2+212.0_pr*gA4)*u3 &
               +420.0_pr*(-3.0_pr-16.0_pr*gA2+27.0_pr*gA4)*u5+33.0_pr*(-1.0_pr-10.0_pr*gA2 &
               +23.0_pr*gA4)*u7+16800.0_pr*(-1.0_pr+4.0_pr*gA2+16.0_pr*gA4)*asqu4*Log(2.0_pr) &
               +1680.0_pr*(-4.0_pr-5.0_pr*gA2+22.0_pr*gA4)*u2*asqu4*Log(2.0_pr)+210.0_pr &
               *(-1.0_pr-10.0_pr*gA2+23.0_pr*gA4)*u4*asqu4*Log(2.0_pr)+105.0_pr*(1.0_pr &
               +10.0_pr*gA2-23.0_pr*gA4)*u6*asqu4*Log(2.0_pr)+105.0_pr*asqu4 &
               *(-((-10.0_pr+u2)*(4.0_pr+u2)**2)-10.0_pr*gA2*(64.0_pr-8.0_pr*u2-2.0_pr*u4+u6) &
               +gA4*(-2560.0_pr-352.0_pr*u2-46.0_pr*u4+23.0_pr*u6))*Log(2.0_pr+u*(u+Sqrt(4.0_pr &
               +u2)))))/(4.128768e7_pr*fpi4*Pi**2*u10)

          darhodrho=darhodrho+(-420.0_pr*u*(1.0_pr+u2)*asqu4*(15.0_pr-158.0_pr*u2 &
               -8.0_pr*u4-2.0_pr*gA2*(75.0_pr-238.0_pr*u2+8.0_pr*u4)+gA4*(-1905.0_pr+8690.0_pr*u2 &
               +152.0_pr*u4))*Log(1.0_pr+2.0_pr*u*(u+asqu))+315.0_pr*asqu &
               *asqu4*(5.0_pr-50.0_pr*gA2-635.0_pr*gA4-36.0_pr*(-1.0_pr+6.0_pr*gA2 &
               +83.0_pr*gA4)*u2+24.0_pr*(-1.0_pr+2.0_pr*gA2+47.0_pr*gA4)*u4)*Log(1.0_pr+2.0_pr*u*(u &
               +asqu))**2+4.0_pr*u2*asqu*(1575.0_pr*asqu4+gA4 &
               *(-200025.0_pr*asqu4+u*(u*(845775.0_pr*asqu4+4.0_pr*u*(u*(-47705.0_pr &
               *asqu4+3.0_pr*u*(u*asqu4*(-1253.0_pr+236.0_pr*u2) &
               -7840.0_pr*Log(2.0_pr)))-268800.0_pr*Log(2.0_pr)))-3440640.0_pr*Log(2.0_pr))) &
               +u*(64575.0_pr*u*asqu4+33740.0_pr*u3*asqu4+1092.0_pr*u5*asqu4 &
               -48.0_pr*u7*asqu4+13440.0_pr*(4.0_pr+u2)**2*Log(2.0_pr))-2.0_pr*gA2*(7875.0_pr &
               *asqu4+u*(138915.0_pr*u*asqu4-8540.0_pr*u3*asqu4 &
               -3276.0_pr*u5*asqu4+240.0_pr*u7*asqu4-13440.0_pr*(-8.0_pr+u2) &
               *(4.0_pr+u2)*Log(2.0_pr)))+13440.0_pr*u*(-(4.0_pr+u2)**2+gA2*(64.0_pr+8.0_pr*u2 &
               -2.0_pr*u4)+gA4*(256.0_pr+80.0_pr*u2+7.0_pr*u4))*Log(2.0_pr+u*(u+Sqrt(4.0_pr &
               +u2)))))/(589824.0_pr*fpi4*Pi**2*u11*asqu*asqu4)

          ddarhodrho=ddarhodrho+(420.0_pr*u*asqu*(4.0_pr+u2)**1.5_pr*(-165.0_pr*(-1.0_pr &
               +10.0_pr*gA2+127.0_pr*gA4)+2.0_pr*(-593.0_pr+1730.0_pr*gA2+32183.0_pr*gA4)*u2 &
               +40.0_pr*(-1.0_pr-2.0_pr*gA2+19.0_pr*gA4)*u4)*Log(1.0_pr+2.0_pr*u*(u+asqu)) &
               -315.0_pr*(4.0_pr+u2)**1.5_pr*(-55.0_pr*(-1.0_pr+10.0_pr*gA2+127.0_pr*gA4)-324.0_pr &
               *(-1.0_pr+6.0_pr*gA2+83.0_pr*gA4)*u2+168.0_pr*(-1.0_pr+2.0_pr*gA2+47.0_pr*gA4)*u4) &
               *Log(1.0_pr+2.0_pr*u*(u+asqu))**2-4.0_pr*u2*(-((4.0_pr+u2)*(-17325.0_pr*asqu4 &
               +u*(u*(-526365.0_pr*asqu4+4.0_pr*u*(u*(-48475.0_pr*asqu4+3.0_pr*u &
               *(u*asqu4*(-273.0_pr+4.0_pr*u2)-5600.0_pr*Log(2.0_pr)))-174720.0_pr*Log(2.0_pr))) &
               -1720320.0_pr*Log(2.0_pr))))+gA4*(-8801100.0_pr*asqu4+21899745.0_pr*u2*Sqrt(4.0_pr &
               +u2)+1906205.0_pr*u4*asqu4-1290772.0_pr*u6*asqu4-33780.0_pr*u8 &
               *asqu4+2832.0_pr*u10*asqu4-13440.0_pr*u*(8192.0_pr+4224.0_pr*u2 &
               +672.0_pr*u4+35.0_pr*u6)*Log(2.0_pr))-2.0_pr*gA2*(4.0_pr+u2)*(86625.0_pr*asqu4 &
               +u*(3440640.0_pr*Log(2.0_pr)+u*(1137465.0_pr*asqu4+4.0_pr*u*(u*(-16975.0_pr*Sqrt(4.0_pr &
               +u2)+3.0_pr*u*(u*asqu4*(-819.0_pr+20.0_pr*u2)-5600.0_pr*Log(2.0_pr))) &
               +107520.0_pr*Log(2.0_pr)))))+13440.0_pr*u*(512.0_pr*(-1.0_pr+4.0_pr*gA2+16.0_pr*gA4) &
               +48.0_pr*(-7.0_pr+16.0_pr*gA2+88.0_pr*gA4)*u2+24.0_pr*(-3.0_pr+gA2+28.0_pr*gA4)*u4 &
               +5.0_pr*(-1.0_pr-2.0_pr*gA2+7.0_pr*gA4)*u6)*Log(2.0_pr+u*(u+asqu4)))) &
               /(589824.0_pr*fpi4*Pi**2*u12*(4.0_pr+u2)**1.5_pr)

          brhodrho=brhodrho -(-44100.0_pr*(-1.0_pr+10.0_pr*gA2+43.0_pr*gA4)*u2+4200.0_pr*(-421.0_pr &
               +1516.0_pr*gA2+6169.0_pr*gA4)*u4+420.0_pr*(-283.0_pr-30.0_pr*gA2+2029.0_pr*gA4)*u6 &
               +1680.0_pr*(-19.0_pr-106.0_pr*gA2+185.0_pr*gA4)*u8-464.0_pr*(-1.0_pr-10.0_pr*gA2 &
               +23.0_pr*gA4)*u10+420.0_pr*u*asqu*(-105.0_pr+1685.0_pr*u2+134.0_pr*u4 &
               +8.0_pr*u6-16.0_pr*u8+10.0_pr*gA2*(105.0_pr-527.0_pr*u2+22.0_pr*u4+8.0_pr*u6 &
               -16.0_pr*u8)+gA4*(4515.0_pr-22235.0_pr*u2-842.0_pr*u4-184.0_pr*u6+368.0_pr*u8)) &
               *Log(1.0_pr+2.0_pr*u*(u+asqu))+11025.0_pr*(1.0_pr+9.0_pr*u2-8.0_pr*u4 &
               +2.0_pr*gA2*(-5.0_pr-27.0_pr*u2+8.0_pr*u4)+gA4*(-43.0_pr-207.0_pr*u2+88.0_pr*u4)) &
               *Log(1.0_pr+2.0_pr*u*(u+asqu))**2+64.0_pr*u3*(33600.0_pr*(-1.0_pr+4.0_pr*gA2 &
               +16.0_pr*gA4)*u+560.0_pr*(-29.0_pr-10.0_pr*gA2+212.0_pr*gA4)*u3+420.0_pr*(-3.0_pr &
               -16.0_pr*gA2+27.0_pr*gA4)*u5+33.0_pr*(-1.0_pr-10.0_pr*gA2+23.0_pr*gA4)*u7 &
               +16800.0_pr*(-1.0_pr+4.0_pr*gA2+16.0_pr*gA4)*asqu4*Log(2.0_pr)+1680.0_pr*(-4.0_pr &
               -5.0_pr*gA2+22.0_pr*gA4)*u2*asqu4*Log(2.0_pr)+210.0_pr*(-1.0_pr-10.0_pr*gA2 &
               +23.0_pr*gA4)*u4*asqu4*Log(2.0_pr)+105.0_pr*(1.0_pr+10.0_pr*gA2 &
               -23.0_pr*gA4)*u6*asqu4*Log(2.0_pr)+105.0_pr*asqu4*(-((-10.0_pr &
               +u2)*(4.0_pr+u2)**2)-10.0_pr*gA2*(64.0_pr-8*u2-2*u4+u6)+gA4*(-2560.0_pr &
               -352.0_pr*u2-46.0_pr*u4+23.0_pr*u6))*Log(2.0_pr+u*(u+asqu4)))) &
               /(2.064384e7_pr*fpi4*Pi**2*u10)

          dbrhodrho=dbrhodrho+(-420.0_pr*u*(1.0_pr+u2)*asqu4*(15.0_pr-158.0_pr*u2 &
               -8.0_pr*u4-2.0_pr*gA2*(75.0_pr-238.0_pr*u2+8.0_pr*u4)+gA4*(-645.0_pr+2042.0_pr*u2 &
               +56.0_pr*u4))*Log(1.0_pr+2.0_pr*u*(u+asqu))+315.0_pr*asqu &
               *asqu4*(5.0_pr-50.0_pr*gA2-215.0_pr*gA4-36.0_pr*(-1.0_pr+6.0_pr*gA2 &
               +23.0_pr*gA4)*u2+24.0_pr*(-1.0_pr+2.0_pr*gA2+11.0_pr*gA4)*u4)*Log(1.0_pr+2.0_pr*u*(u &
               +asqu))**2-4.0_pr*u2*asqu*(-1575.0_pr*asqu4 &
               +gA4*(67725.0_pr*asqu4+u*(u*(-1482075.0_pr*asqu4+4.0_pr*u*(u &
               *(-35035.0_pr*asqu4+3.0_pr*u*(u*asqu4*(-847.0_pr+184.0_pr*u2) &
               -3920.0_pr*Log(2.0_pr)))-134400.0_pr*Log(2.0_pr)))-1720320.0_pr*Log(2.0_pr))) &
               +u*(96705.0_pr*u*asqu4+20020.0_pr*u3*asqu4+924.0_pr*u5*asqu4 &
               -96.0_pr*u7*asqu4+6720.0_pr*(4.0_pr+u2)**2*Log(2.0_pr))+2.0_pr*gA2*(7875.0_pr &
               *asqu4+u*(-183645.0_pr*u*asqu4+4900.0_pr*u3*asqu4 &
               +2772.0_pr*u5*asqu4-480.0_pr*u7*asqu4+6720.0_pr*(-8.0_pr+u2) &
               *(4.0_pr+u2)*Log(2.0_pr)))+6720.0_pr*u*(-(4.0_pr+u2)**2+gA2*(64.0_pr+8.0_pr*u2 &
               -2.0_pr*u4)+gA4*(256.0_pr+80.0_pr*u2+7.0_pr*u4))*Log(2.0_pr+u*(u+asqu4)))) &
               /(294912.0_pr*fpi4*Pi**2*u11*asqu*asqu4)

          ddbrhodrho=ddbrhodrho+(420.0_pr*u*asqu*(4.0_pr+u2)**1.5_pr*(-165.0_pr*(-1.0_pr &
               +10.0_pr*gA2+43.0_pr*gA4)+2.0_pr*(-593.0_pr+1730.0_pr*gA2+7571.0_pr*gA4)*u2 &
               +40.0_pr*(-1.0_pr-2.0_pr*gA2+7.0_pr*gA4)*u4)*Log(1.0_pr+2.0_pr*u*(u+asqu)) &
               -315.0_pr*(4.0_pr+u2)**1.5_pr*(-55.0_pr*(-1.0_pr+10.0_pr*gA2+43.0_pr*gA4)-324.0_pr*(-1.0_pr &
               +6.0_pr*gA2+23.0_pr*gA4)*u2+168.0_pr*(-1.0_pr+2.0_pr*gA2+11.0_pr*gA4)*u4)*Log(1.0_pr &
               +2.0_pr*u*(u+asqu))**2+4.0_pr*u2*(-((4.0_pr+u2)*(17325.0_pr*asqu4 &
               +u*(u*(-763875.0_pr*asqu4+4.0_pr*u*(u*(-28805.0_pr*asqu4+3.0_pr*u*(u &
               *asqu4*(-231.0_pr+8*u2)-2800.0_pr*Log(2.0_pr)))-87360.0_pr*Log(2.0_pr))) &
               -860160.0_pr*Log(2.0_pr))))+ gA4*(2979900.0_pr*asqu4-45909045.0_pr*u2*Sqrt(4.0_pr &
               +u2)-14728385.0_pr*u4*asqu4-928508.0_pr*u6*asqu4-21660.0_pr*u8 &
               *asqu4+2208.0_pr*u10*asqu4-6720.0_pr*u*(8192.0_pr+4224.0_pr*u2 &
               +672.0_pr*u4+35.0_pr*u6)*Log(2.0_pr))-2.0_pr*gA2*(4.0_pr+u2)*(-86625.0_pr*asqu4 &
               +u*(1720320.0_pr*Log(2.0_pr)+u*(1443015.0_pr*asqu4+4.0_pr*u*(u*(-9905.0_pr*Sqrt(4.0_pr &
               +u2)+3.0_pr*u*(u*asqu4*(-693.0_pr+40.0_pr*u2)-2800.0_pr*Log(2.0_pr))) &
               +53760.0_pr*Log(2.0_pr)))))+6720.0_pr*u*(512.0_pr*(-1.0_pr+4.0_pr*gA2+16.0_pr*gA4) &
               +48.0_pr*(-7.0_pr+16.0_pr*gA2+88.0_pr*gA4)*u2+24.0_pr*(-3.0_pr+gA2+28.0_pr*gA4)*u4 &
               +5.0_pr*(-1.0_pr-2.0_pr*gA2+7.0_pr*gA4)*u6)*Log(2.0_pr+u*(u+asqu4)))) &
               /(294912.0_pr*fpi4*Pi**2*u12*(4.0_pr+u2)**1.5_pr)

          adrdr=adrdr +0.0_pr
          dadrdr=dadrdr +0.0_pr
          bdrdr=bdrdr +0.0_pr
          dbdrdr=dbdrdr +0.0_pr

          arhotau=arhotau+(-44100.0_pr*(-1.0_pr+10.0_pr*gA2+127.0_pr*gA4)*u2+4200.0_pr*(-421.0_pr &
               +1516.0_pr*gA2+25369.0_pr*gA4)*u4+420.0_pr*(-283.0_pr-30.0_pr*gA2+8641.0_pr*gA4)*u6 &
               +1680.0_pr*(-19.0_pr-106.0_pr*gA2+5.0_pr*gA4)*u8+464.0_pr*(1.0_pr+10.0_pr*gA2 &
               +13.0_pr*gA4)*u10-420.0_pr*u*asqu*(105.0_pr-1685.0_pr*u2-134.0_pr*u4 &
               -8.0_pr*u6+16.0_pr*u8+10.0_pr*gA2*(-105.0_pr+527.0_pr*u2-22.0_pr*u4-8.0_pr*u6 &
               +16.0_pr*u8)+gA4*(-13335.0_pr+94295.0_pr*u2+2738.0_pr*u4-104.0_pr*u6+208.0_pr*u8)) &
               *Log(1.0_pr+2.0_pr*u*(u+asqu))+11025.0_pr*(1.0_pr-10.0_pr*gA2-127.0_pr*gA4 &
               -9.0_pr*(-1.0_pr+6.0_pr*gA2+83.0_pr*gA4)*u2+8.0_pr*(-1.0_pr+2.0_pr*gA2 &
               +47.0_pr*gA4)*u4)*Log(1.0_pr+2.0_pr*u*(u+asqu))**2)/(1.032192e7_pr*fpi4*Pi**2*u10)

          darhotau=darhotau+(-4.0_pr*u2*(-525.0_pr*(-1.0_pr+10.0_pr*gA2+127.0_pr*gA4)+35.0_pr &
               *(-409.0_pr+1450.0_pr*gA2+24439.0_pr*gA4)*u2+140.0_pr*(-5.0_pr-2.0_pr*gA2 &
               +143.0_pr*gA4)*u4-84.0_pr*(1.0_pr+6.0_pr*gA2+gA4)*u6+16.0_pr*(1.0_pr+10.0_pr*gA2 &
               +13.0_pr*gA4)*u8)*(asqu+2.0_pr*u*(1.0_pr+u*(u+asqu)))+35.0_pr &
               *Log(1.0_pr+2.0_pr*u*(u+asqu))*(4.0_pr*u*(1.0_pr+u2)*(15.0_pr-158.0_pr*u2 &
               -8.0_pr*u4-2.0_pr*gA2*(75.0_pr-238.0_pr*u2+8.0_pr*u4)+gA4*(-1905.0_pr+8690.0_pr*u2 &
               +152.0_pr*u4))*(1.0_pr+2.0_pr*u*(u+asqu))-3.0_pr*(5.0_pr-50.0_pr*gA2 &
               -635.0_pr*gA4-36.0_pr*(-1.0_pr+6.0_pr*gA2+83.0_pr*gA4)*u2+24.0_pr*(-1.0_pr+2.0_pr*gA2 &
               +47.0_pr*gA4)*u4)*(asqu+2.0_pr*u*(1.0_pr+u*(u+asqu)))*Log(1.0_pr &
               +2.0_pr*u*(u+asqu))))/(49152.0_pr*fpi4*Pi**2*u11*asqu*(1.0_pr &
               +2.0_pr*u*(u+asqu)))

          ddarhotau=ddarhotau+(4.0_pr*u2*(-5775.0_pr*(-1.0_pr+10.0_pr*gA2+127.0_pr*gA4)+35.0_pr &
               *(-3179.0_pr+11102.0_pr*gA2+188453.0_pr*gA4)*u2+140.0_pr*(-29.0_pr-18.0_pr*gA2 &
               +791.0_pr*gA4)*u4-252.0_pr*(1.0_pr+6.0_pr*gA2+gA4)*u6+16.0_pr*(1.0_pr+10.0_pr*gA2 &
               +13.0_pr*gA4)*u8)*(asqu+4.0_pr*(u+u3)*(1.0_pr+2.0_pr*u*(u+Sqrt(1.0_pr &
               +u2))))-35.0_pr*Log(1.0_pr+2.0_pr*u*(u+asqu))*(4.0_pr*u*(1.0_pr+u2)*(-165.0_pr &
               *(-1.0_pr+10.0_pr*gA2+127.0_pr*gA4)+2.0_pr*(-593.0_pr+1730.0_pr*gA2+32183.0_pr*gA4)*u2 &
               +40.0_pr*(-1.0_pr-2.0_pr*gA2+19.0_pr*gA4)*u4)*(1.0_pr+4.0_pr*u*(asqu &
               +2.0_pr*u*(1.0_pr+u*(u+asqu))))-3.0_pr*(-55.0_pr*(-1.0_pr+10.0_pr*gA2 &
               +127.0_pr*gA4)-324.0_pr*(-1.0_pr+6.0_pr*gA2+83.0_pr*gA4)*u2+168.0_pr*(-1.0_pr &
               +2.0_pr*gA2+47.0_pr*gA4)*u4)*(asqu+4.0_pr*(u+u3)*(1.0_pr+2.0_pr*u*(u &
               +asqu)))*Log(1.0_pr+2.0_pr*u*(u+asqu))))/(49152.0_pr*fpi4*Pi**2*u12 &
               *asqu*(1.0_pr+2.0_pr*u*(u+asqu))**2)

          brhotau=brhotau+(-44100.0_pr*(-1.0_pr+10.0_pr*gA2+43.0_pr*gA4)*u2+4200.0_pr*(-421.0_pr &
               +1516.0_pr*gA2+6169.0_pr*gA4)*u4+420.0_pr*(-283.0_pr-30.0_pr*gA2+2029.0_pr*gA4)*u6 &
               +1680.0_pr*(-19.0_pr-106.0_pr*gA2+185.0_pr*gA4)*u8-464.0_pr*(-1.0_pr-10.0_pr*gA2 &
               +23.0_pr*gA4)*u10+105.0_pr*Log(1.0_pr+2.0_pr*u*(u+asqu))*(4.0_pr*u*Sqrt(1.0_pr &
               +u2)*(-105.0_pr+1685.0_pr*u2+134.0_pr*u4+8.0_pr*u6-16.0_pr*u8+10.0_pr*gA2*(105.0_pr &
               -527.0_pr*u2+22.0_pr*u4+8.0_pr*u6-16.0_pr*u8)+gA4*(4515.0_pr-22235.0_pr*u2 &
               -842.0_pr*u4-184.0_pr*u6+368.0_pr*u8))+105.0_pr*(1.0_pr+9.0_pr*u2-8.0_pr*u4 &
               +2.0_pr*gA2*(-5.0_pr-27.0_pr*u2+8.0_pr*u4)+gA4*(-43.0_pr-207.0_pr*u2+88.0_pr*u4)) &
               *Log(1.0_pr+2.0_pr*u*(u+asqu))))/(5.16096e6*fpi4*Pi**2*u10)

          dbrhotau=dbrhotau+(4.0_pr*u2*(525.0_pr*(-1.0_pr+10.0_pr*gA2+43.0_pr*gA4)-35.0_pr*(-409.0_pr &
               +1450.0_pr*gA2+5923.0_pr*gA4)*u2+140.0_pr*(5.0_pr+2.0_pr*gA2-35.0_pr*gA4)*u4 &
               +84.0_pr*(1.0_pr+6.0_pr*gA2-11.0_pr*gA4)*u6+16.0_pr*(-1.0_pr-10.0_pr*gA2 &
               +23.0_pr*gA4)*u8)*(asqu+2.0_pr*u*(1.0_pr+u*(u+asqu)))+35.0_pr &
               *Log(1.0_pr+2.0_pr*u*(u+asqu))*(4.0_pr*u*(1.0_pr+u2)*(15.0_pr-158.0_pr*u2 &
               -8.0_pr*u4-2.0_pr*gA2*(75.0_pr-238.0_pr*u2+8.0_pr*u4)+gA4*(-645.0_pr+2042.0_pr*u2 &
               +56.0_pr*u4))*(1.0_pr+2.0_pr*u*(u+asqu))-3.0_pr*(5.0_pr-50.0_pr*gA2 &
               -215.0_pr*gA4-36.0_pr*(-1.0_pr+6.0_pr*gA2+23.0_pr*gA4)*u2+24.0_pr*(-1.0_pr+2.0_pr*gA2 &
               +11.0_pr*gA4)*u4)*(asqu+2.0_pr*u*(1.0_pr+u*(u+asqu)))*Log(1.0_pr &
               +2.0_pr*u*(u+asqu))))/(24576.0_pr*fpi4*Pi**2*u11*asqu*(1.0_pr &
               +2.0_pr*u*(u+asqu)))

          ddbrhotau=ddbrhotau+(-4.0_pr*u2*(5775.0_pr*(-1.0_pr+10.0_pr*gA2+43.0_pr*gA4)-35.0_pr &
               *(-3179.0_pr+11102.0_pr*gA2+45545.0_pr*gA4)*u2+140.0_pr*(29.0_pr+18.0_pr*gA2 &
               -203.0_pr*gA4)*u4-252.0_pr*(-1.0_pr-6.0_pr*gA2+11.0_pr*gA4)*u6+16.0_pr*(-1.0_pr &
               -10.0_pr*gA2+23.0_pr*gA4)*u8)*(asqu+4.0_pr*(u+u3)*(1.0_pr+2.0_pr*u*(u &
               +asqu)))-35.0_pr*Log(1.0_pr+2.0_pr*u*(u+asqu))*(4.0_pr*u*(1.0_pr &
               +u2)*(-165.0_pr*(-1.0_pr+10.0_pr*gA2+43.0_pr*gA4)+2.0_pr*(-593.0_pr+1730.0_pr*gA2 &
               +7571.0_pr*gA4)*u2+40.0_pr*(-1.0_pr-2.0_pr*gA2+7.0_pr*gA4)*u4)*(1.0_pr+4.0_pr*u &
               *(asqu+2.0_pr*u*(1.0_pr+u*(u+asqu))))-3.0_pr*(-55.0_pr*(-1.0_pr &
               +10.0_pr*gA2+43.0_pr*gA4)-324.0_pr*(-1.0_pr+6.0_pr*gA2+23.0_pr*gA4)*u2 &
               +168.0_pr*(-1.0_pr+2.0_pr*gA2+11.0_pr*gA4)*u4)*(asqu+4.0_pr*(u+u3)*(1.0_pr &
               +2.0_pr*u*(u+asqu)))*Log(1.0_pr+2.0_pr*u*(u+asqu)))) &
               /(24576.0_pr*fpi4*Pi**2*u12*asqu*(1.0_pr+2.0_pr*u*(u+asqu))**2)

          ajj=ajj+(4.0_pr*u2*(-9.0_pr+24.0_pr*u2+u4+gA4*(18.0_pr-363.0_pr*u2-113.0_pr*u4) &
               +2.0_pr*gA2*(9.0_pr+66.0_pr*u2+5.0_pr*u4))-12.0_pr*u*asqu*(-3.0_pr+u2 &
               -2.0_pr*u4+2.0_pr*gA2*(3.0_pr+5.0_pr*u2-10.0_pr*u4)+gA4*(6.0_pr-59.0_pr*u2 &
               +10.0_pr*u4))*Log(1.0_pr+2.0_pr*u*(u+asqu))+9.0_pr*(-1.0_pr+2.0_pr*gA2 &
               +2.0_pr*gA4+3.0_pr*(-1.0_pr-2.0_pr*gA2+gA4)*u2)*Log(1.0_pr+2.0_pr*u*(u &
               +asqu))**2)/(73728.0_pr*fpi4*Pi**2*u6)

          dajj=dajj+(-4.0_pr*u2*(-9.0_pr+9.0_pr*u2-2.0_pr*u4+2.0_pr*gA2*(9.0_pr+27.0_pr*u2 &
               -10.0_pr*u4)+2.0_pr*gA4*(9.0_pr-90.0_pr*u2+5.0_pr*u4))*(asqu+2.0_pr*u &
               *(1.0_pr+u*(u+asqu)))-9.0_pr*Log(1.0_pr+2.0_pr*u*(u+asqu)) &
               *(4.0_pr*u*(1.0_pr+u2-2.0_pr*gA2*(1.0_pr+u2)+gA4*(-2.0_pr+8.0_pr*u2+6.0_pr*u4)) &
               *(1.0_pr+2.0_pr*u*(u+asqu))+(-1.0_pr-2.0_pr*u2+2.0_pr*gA2*(1.0_pr+gA2 &
               +(-2.0_pr+gA2)*u2))*(asqu+2.0_pr*u*(1.0_pr+u*(u+asqu))) &
               *Log(1.0_pr+2.0_pr*u*(u+asqu))))/(12288.0_pr*fpi4*Pi**2*u7*asqu &
               *(1.0_pr+2.0_pr*u*(u+asqu)))

          bjj=bjj+(4.0_pr*u2*(-9.0_pr+24.0_pr*u2+u4+gA4*(99.0_pr-228.0_pr*u2-23.0_pr*u4) &
               +2.0_pr*gA2*(9.0_pr+66.0_pr*u2+5.0_pr*u4))-12.0_pr*u*asqu*(-3.0_pr+u2 &
               -2.0_pr*u4+2.0_pr*gA2*(3.0_pr+5.0_pr*u2-10.0_pr*u4)+gA4*(33.0_pr-23.0_pr*u2 &
               +46.0_pr*u4))*Log(1.0_pr+2.0_pr*u*(u+asqu))+9.0_pr*(-1.0_pr+2.0_pr*gA2 &
               +11.0_pr*gA4+3.0_pr*(-1.0_pr-2.0_pr*gA2+7.0_pr*gA4)*u2)*Log(1.0_pr+2.0_pr*u*(u &
               +asqu))**2)/(36864.0_pr*fpi4*Pi**2*u6)

          dbjj=dbjj+(-4.0_pr*u2*(-9.0_pr+9.0_pr*u2-2.0_pr*u4+2.0_pr*gA2*(9.0_pr+27.0_pr*u2 &
               -10.0_pr*u4)+gA4*(99.0_pr-99.0_pr*u2+46.0_pr*u4))*(asqu+2.0_pr*u*(1.0_pr &
               +u*(u+asqu)))+36.0_pr*u*(-1.0_pr-u2+2.0_pr*gA2*(1.0_pr+u2)+gA4*(11.0_pr &
               +7.0_pr*u2))*(1.0_pr+2.0_pr*u*(u+asqu))*Log(1.0_pr+2.0_pr*u*(u+Sqrt(1.0_pr &
               +u2)))-9.0_pr*(-1.0_pr+2.0_pr*gA2+11.0_pr*gA4+2.0_pr*(-1.0_pr-2.0_pr*gA2 &
               +7.0_pr*gA4)*u2)*(asqu+2.0_pr*u*(1.0_pr+u*(u+asqu)))*Log(1.0_pr &
               +2.0_pr*u*(u+asqu))**2)/(6144.0_pr*fpi4*Pi**2*u7*asqu*(1.0_pr &
               +2.0_pr*u*(u+asqu)))
       Else    !! u smaller than ucut
          arhorho=arhorho +((-1.0_pr - 4.0_pr*gA2 + 8.0_pr*gA**4)*mpi2)/(384.0_pr*fpi4*Pi**2) &
               + (11.0_pr*(1.0_pr + 14.0_pr*gA2 + 17.0_pr*gA**4)*mpi2*u4)/(201600.0_pr*fpi4*Pi**2) &
               - ((1.0_pr + 18.0_pr*gA2 + 13.0_pr*gA**4)*mpi2*u6)/(69300.0_pr*fpi4*Pi**2)
          darhorho=darhorho+(mpi2*u3*(121.0_pr*(1.0_pr+14.0_pr*gA2+17.0_pr*gA4)-48.0_pr*(1.0_pr &
               +18.0_pr*gA2+13.0_pr*gA4)*u2))/(554400.0_pr*fpi4*Pi**2)
          ddarhorho=ddarhorho+(mpi2*u2*(1573.0_pr*(1.0_pr+14.0_pr*gA2+17.0_pr*gA4)-1040.0_pr*(1.0_pr &
               +18.0_pr*gA2+13.0_pr*gA4)*u2+680.0_pr*(1.0_pr+22.0_pr*gA2+gA4)*u4)) &
               /(2.4024e6_pr*fpi4*Pi**2)

          brhorho=brhorho +((-1.0_pr - 4.0_pr*gA2 + 8.0_pr*gA**4)*mpi2)/(192.0_pr*fpi4*Pi**2) &
               - (11.0_pr*(-1.0_pr - 14.0_pr*gA2 + 43.0_pr*gA**4)*mpi2*u4)/(100800.0_pr*fpi4*Pi**2) &
               + ((-1.0_pr - 18.0_pr*gA2 + 71.0_pr*gA**4)*mpi2*u6)/(34650.0_pr*fpi4*Pi**2)
          dbrhorho=dbrhorho+(mpi2*u3*(121.0_pr*(1.0_pr+14.0_pr*gA2-43.0_pr*gA4)+48.0_pr*(-1.0_pr &
               -18.0_pr*gA2+71.0_pr*gA4)*u2))/(277200.0_pr*fpi4*Pi**2)
          ddbrhorho=ddbrhorho+(mpi2*u2*(-1573.0_pr*(-1.0_pr-14.0_pr*gA2+43.0_pr*gA4)+1040.0_pr &
               *(-1.0_pr-18.0_pr*gA2+71.0_pr*gA4)*u2-680.0_pr*(-1.0_pr-22.0_pr*gA2 &
               +107.0_pr*gA4)*u4))/(1.2012e6_pr*fpi4*Pi**2)

          arhodrho=arhodrho+(120120.0_pr*(-2.0_pr-17.0_pr*gA2+88.0_pr*gA4)+8008.0_pr*(1.0_pr &
               +14.0_pr*gA2+167.0_pr*gA4)*u2-39.0_pr*(109.0_pr+1962.0_pr*gA2+4357.0_pr*gA4)*u4 &
               +540.0_pr*(3.0_pr+66.0_pr*gA2+31.0_pr*gA4)*u6)/(5.5351296e8_pr*fpi4*Pi**2)
          darhodrho=darhodrho+(4004.0_pr*(1.0_pr+14.0_pr*gA2+167.0_pr*gA4)*u-39.0_pr*(109.0_pr &
               +1962.0_pr*gA2+4357.0_pr*gA4)*u3+810.0_pr*(3.0_pr+66.0_pr*gA2+31.0_pr*gA4)*u5) &
               /(1.3837824e8_pr*fpi4*Pi**2)
          ddarhodrho=ddarhodrho+(12012.0_pr*(1.0_pr+14.0_pr*gA2+167.0_pr*gA4)-351.0_pr*(109.0_pr &
               +1962.0_pr*gA2+4357.0_pr*gA4)*u2+12150.0_pr*(3.0_pr+66.0_pr*gA2+31.0_pr*gA4)*u4 &
               +35.0_pr*(-811.0_pr-21086.0_pr*gA2+12637.0_pr*gA4)*u6)/(4.1513472e8_pr*fpi4*Pi**2)

          brhodrho=brhodrho -(480480.0_pr*(-2.0_pr-17.0_pr*gA2+34.0_pr*gA4)+76076.0_pr*(-1.0_pr &
               -14.0_pr*gA2+43.0_pr*gA4)*u2-12597.0_pr*(-1.0_pr-18.0_pr*gA2+71.0_pr*gA4)*u4 &
               +3660.0_pr*(-1.0_pr-22.0_pr*gA2+107.0_pr*gA4)*u6)/(5.5351296e8_pr*fpi4*Pi**2)
          dbrhodrho=dbrhodrho+(-38038.0_pr*(-1.0_pr-14.0_pr*gA2+43.0_pr*gA4)*u+12597.0_pr*(-1.0_pr &
               -18.0_pr*gA2+71.0_pr*gA4)*u3-5490.0_pr*(-1.0_pr-22.0_pr*gA2+107.0_pr*gA4)*u5) &
               /(1.3837824e8_pr*fpi4*Pi**2)
          ddbrhodrho=ddbrhodrho+(-114114.0_pr*(-1.0_pr-14.0_pr*gA2+43.0_pr*gA4)+113373.0_pr*(-1.0_pr &
               -18.0_pr*gA2+71.0_pr*gA4)*u2-82350.0_pr*(-1.0_pr-22.0_pr*gA2+107.0_pr*gA4)*u4 &
               +58975.0_pr*(-1.0_pr-26.0_pr*gA2+151.0_pr*gA4)*u6)/(4.1513472e8_pr*fpi4*Pi**2)

          adrdr=adrdr +0.0_pr
          dadrdr=dadrdr +0.0_pr
          bdrdr=bdrdr +0.0_pr
          dbdrdr=dbdrdr +0.0_pr

          arhotau=arhotau -(30030.0_pr*(2.0_pr+17.0_pr*gA2+20.0_pr*gA4)+7007.0_pr*(1.0_pr+14.0_pr*gA2 &
               +17.0_pr*gA4)*u2-1404.0_pr*(1.0_pr+18.0_pr*gA2+13.0_pr*gA4)*u4+440.0_pr*(1.0_pr &
               +22.0_pr*gA2+gA4)*u6)/(3.459456e7_pr*fpi4*Pi**2)
          darhotau=darhotau+(-7007.0_pr*(1.0_pr+14.0_pr*gA2+17.0_pr*gA4)*u+2808.0_pr*(1.0_pr &
               +18.0_pr*gA2+13.0_pr*gA4)*u3-1320.0_pr*(1.0_pr+22.0_pr*gA2+gA4)*u5) &
               /(1.729728e7_pr*fpi4*Pi**2)
          ddarhotau=ddarhotau -(21021.0_pr*(1.0_pr+14.0_pr*gA2+17.0_pr*gA4)-25272.0_pr*(1.0_pr &
               +18.0_pr*gA2+13.0_pr*gA4)*u2+19800.0_pr*(1.0_pr+22.0_pr*gA2+gA4)*u4+14560.0_pr &
               *(-1.0_pr-26.0_pr*gA2+19.0_pr*gA4)*u6)/(5.189184e7_pr*fpi4*Pi**2)

          brhotau=brhotau+(30030.0_pr*(-2.0_pr-17.0_pr*gA2+34.0_pr*gA4)+7007.0_pr*(-1.0_pr &
               -14.0_pr*gA2+43.0_pr*gA4)*u2-1404.0_pr*(-1.0_pr-18.0_pr*gA2+71.0_pr*gA4)*u4 &
               +440.0_pr*(-1.0_pr-22.0_pr*gA2+107.0_pr*gA4)*u6)/(1.729728e7_pr*fpi4*Pi**2)
          dbrhotau=dbrhotau+(7007.0_pr*(-1.0_pr-14.0_pr*gA2+43.0_pr*gA4)*u-2808.0_pr*(-1.0_pr &
               -18.0_pr*gA2+71.0_pr*gA4)*u3+1320.0_pr*(-1.0_pr-22.0_pr*gA2+107.0_pr*gA4)*u5) &
               /(8.64864e6_pr*fpi4*Pi**2)
          ddbrhotau=ddbrhotau+(21021.0_pr*(-1.0_pr-14.0_pr*gA2+43.0_pr*gA4)-25272.0_pr*(-1.0_pr &
               -18.0_pr*gA2+71.0_pr*gA4)*u2+19800.0_pr*(-1.0_pr-22.0_pr*gA2+107.0_pr*gA4)*u4 &
               -14560.0_pr*(-1.0_pr-26.0_pr*gA2+151.0_pr*gA4)*u6)/(2.594592e7_pr*fpi4*Pi**2)

          ajj=ajj+(-1050.0_pr*(-2.0_pr-17.0_pr*gA2+7.0_pr*gA4)+315.0_pr*(1.0_pr+14.0_pr*gA2 &
               -28.0_pr*gA4)*u2+81.0_pr*(-1.0_pr-18.0_pr*gA2+57.0_pr*gA4)*u4+16.0_pr*(2.0_pr &
               +44.0_pr*gA2-187.0_pr*gA4)*u6)/(2.4192e6_pr*fpi4*Pi**2)
          dajj=dajj+(105.0_pr*(1.0_pr+14.0_pr*gA2-28.0_pr*gA4)*u+54.0_pr*(-1.0_pr-18.0_pr*gA2 &
               +57.0_pr*gA4)*u3+16.0_pr*(2.0_pr+44.0_pr*gA2-187.0_pr*gA4)*u5)/(403200.0_pr*fpi4*Pi**2)

          bjj=bjj+(-1050.0_pr*(-2.0_pr-17.0_pr*gA2+34.0_pr*gA4)-315.0_pr*(-1.0_pr-14.0_pr*gA2 &
               +43.0_pr*gA4)*u2+81.0_pr*(-1.0_pr-18.0_pr*gA2+71.0_pr*gA4)*u4+32.0_pr*(1.0_pr &
               +22.0_pr*gA2-107.0_pr*gA4)*u6)/(1.2096e6_pr*fpi4*Pi**2)
          dbjj=dbjj+(105.0_pr*(1.0_pr+14.0_pr*gA2-43.0_pr*gA4)*u+54.0_pr*(-1.0_pr-18.0_pr*gA2 &
               +71.0_pr*gA4)*u3+32.0_pr*(1.0_pr+22.0_pr*gA2-107.0_pr*gA4)*u5)/(201600.0_pr*fpi4*Pi**2)

       End If     !! if (u.gt.ucut...
    End If    !! if (dmeorder.ge.1....
    !
    ! N2LO, 2N terms
    !
    If (dmeorder.Ge.2) Then
       If (u.Gt.ucut) Then
          arhorho=arhorho +(gA2*mpi**3*(u2*(210.0_pr*(-129.0_pr*c3 + 56.0_pr*c4) + 6.0_pr*(38607.0_pr*c3 &
               - 27325.0_pr*c4)*u2 + 7.0_pr*(3165.0_pr*c3 - 2014.0_pr*c4)*u4 + 2714.0_pr*(3.0_pr*c3 &
               + c4)*u6 + 570.0_pr*(3.0_pr*c3 + c4)*u8 - 66.0_pr*c1*(1050.0_pr - 6366.0_pr*u2 + 299.0_pr*u4 &
               + 237.0_pr*u6)) + 3.0_pr*u3*(-242550.0_pr*c1 + 6237.0_pr*(2.0_pr*c1 - c3)*u4 &
               + 572.0_pr*(3.0_pr*c1 - 3.0_pr*c3 - c4)*u6 + 190.0_pr*(3.0_pr*c3 + c4)*u8 - 1155.0_pr*(3.0_pr*c3 &
               - 2.0_pr*c4)*(35.0_pr + 2.0_pr*u2))*Atan(u) - 6.0_pr*(-35.0_pr*(330.0_pr*c1 + 129.0_pr*c3 &
               - 56.0_pr*c4) - 352.0_pr*(162.0_pr*c1 + 69.0_pr*c3 - 40.0_pr*c4)*u2 + 891.0_pr*(14.0_pr*c1 &
               + 11.0_pr*c3 - 8.0_pr*c4)*u4)*Log(1.0_pr + u2)))/(110880.0_pr*fpi4*Pi*u8)

          darhorho=darhorho+(gA2*mpi**3*(1680.0_pr*(330.0_pr*c1+129.0_pr*c3-56.0_pr*c4)*u2 &
               -3.0_pr*(620862.0_pr*c1+351039.0_pr*c3-250970.0_pr*c4)*u4+(71742.0_pr*c1-56163.0_pr*c3 &
               +30482.0_pr*c4)*u6+18.0_pr*(286.0_pr*c1-127.0_pr*(3.0_pr*c3+c4))*u8+1710.0_pr*(3.0_pr*c3 &
               +c4)*u10+3.0_pr*u3*(1212750.0_pr*c1-6237.0_pr*(2.0_pr*c1-c3)*u4+572.0_pr*(3.0_pr*c1 &
               -3.0_pr*c3-c4)*u6+570.0_pr*(3.0_pr*c3+c4)*u8+1155.0_pr*(3.0_pr*c3-2.0_pr*c4)*(175.0_pr &
               +6.0_pr*u2))*Atan(u)+24.0_pr*(-70.0_pr*(330.0_pr*c1+129.0_pr*c3-56.0_pr*c4)-528.0_pr*(162.0_pr*c1 &
               +69.0_pr*c3-40.0_pr*c4)*u2+891.0_pr*(14.0_pr*c1+11.0_pr*c3-8.0_pr*c4)*u4)*Log(1.0_pr+u2))) &
               /(110880.0_pr*fpi4*Pi*u9)

          ddarhorho=ddarhorho+(gA2*mpi**3*(-2520.0_pr*(330.0_pr*c1+129.0_pr*c3-56.0_pr*c4)*u2 &
               +21.0_pr*(79002.0_pr*c1+45789.0_pr*c3-32950.0_pr*c4)*u4-9.0_pr*(4774.0_pr*c1-3657.0_pr*c3 &
               +1630.0_pr*c4)*u6-190.0_pr*(3.0_pr*c3+c4)*u8+570.0_pr*(3.0_pr*c3+c4)*u10 &
               +3.0_pr*(u3*(-1212750.0_pr*c1-2079.0_pr*(-2.0_pr*c1+c3)*u4+190.0_pr*(3.0_pr*c3+c4)*u8 &
               -1155.0_pr*(3.0_pr*c3-2.0_pr*c4)*(175.0_pr+4.0_pr*u2))*Atan(u)+4.0_pr*(210.0_pr*(330.0_pr*c1 &
               +129.0_pr*c3-56.0_pr*c4)+1232.0_pr*(162.0_pr*c1+69.0_pr*c3-40.0_pr*c4)*u2-1485.0_pr*(14.0_pr*c1 &
               +11.0_pr*c3-8.0_pr*c4)*u4)*Log(1.0_pr+u2))))/(18480.0_pr*fpi4*Pi*u10)

          brhorho=brhorho +(c4*gA2*mpi**3*(u2*(5880.0_pr - 81975.0_pr*u2 - 7049.0_pr*u4 + 1357.0_pr*u6 &
               + 285.0_pr*u8) + 3.0_pr*u3*(40425.0_pr + 2310.0_pr*u2 - 286.0_pr*u6 + 95.0_pr*u8)*Atan(u) &
               + 24.0_pr*(-245.0_pr - 1760.0_pr*u2 + 891.0_pr*u4)*Log(1.0_pr + u2)))/(27720.0_pr*fpi4*Pi*u8)

          dbrhorho=dbrhorho+(c4*gA2*mpi**3*(u2*(-47040.0_pr+376455.0_pr*u2+15241.0_pr*u4 &
               -1143.0_pr*u6+855.0_pr*u8)+3.0_pr*u3*(-202125.0_pr-6930.0_pr*u2-286.0_pr*u6 &
               +285.0_pr*u8)*Atan(u)-96.0_pr*(-490.0_pr-2640.0_pr*u2+891.0_pr*u4)*Log(1.0_pr+u2))) &
               /(27720.0_pr*fpi4*Pi*u9)

          ddbrhorho=ddbrhorho+(c4*gA2*mpi**3*(u2*(14112.0_pr-69195.0_pr*u2-1467.0_pr*u4-19.0_pr*u6 &
               +57.0_pr*u8)+3.0_pr*u3*(40425.0_pr+924.0_pr*u2+19.0_pr*u8)*Atan(u)+48.0_pr*(-294.0_pr &
               -1232.0_pr*u2+297.0_pr*u4)*Log(1.0_pr+u2)))/(924.0_pr*fpi4*Pi*u10)

          arhodrho=arhodrho+(gA2*mpi*(10395.0_pr*u3*(4.0_pr+u2)*(8.0_pr*c1*(40.0_pr-10.0_pr*u2+u4) &
               +c3*(160.0_pr+8.0_pr*u2+4.0_pr*u4-3.0_pr*u6))*Atan(0.5_pr*u)-2.0_pr*(480.0_pr*(330.0_pr*c1 &
               +129.0_pr*c3-56.0_pr*c4)*u2+48.0_pr*(49038.0_pr*c1+23511.0_pr*c3+7870.0_pr*c4)*u4 &
               +32.0_pr*(-7953.0_pr*c1+10545.0_pr*c3+479.0_pr*c4)*u6-12.0_pr*(5874.0_pr*c1-2769.0_pr*c3 &
               +232.0_pr*c4)*u8+27.0_pr*(147.0_pr*c3+16.0_pr*c4)*u10-48.0_pr*u3*(-5775.0_pr*(6.0_pr*c1 &
               +3.0_pr*c3-2.0_pr*c4)+297.0_pr*(2.0_pr*c1-c3)*u4+88.0_pr*(-3.0_pr*c1+3.0_pr*c3+c4)*u6 &
               +90.0_pr*(3.0_pr*c3+c4)*u8)*Atan(u)+48.0_pr*(-1290.0_pr*c3+560.0_pr*c4+33.0_pr*(-100.0_pr*c1 &
               -3.0_pr*(162.0_pr*c1+69.0_pr*c3-40.0_pr*c4)*u2+9.0_pr*(14.0_pr*c1+11.0_pr*c3 &
               -8.0_pr*c4)*u4))*Log(1.0_pr+u2))))/(1.216512e6_pr*fpi4*Pi*u10)

          darhodrho=darhodrho+(gA2*mpi*(-10395.0_pr*u3*(8960.0_pr*c1+72.0_pr*(-2.0_pr*c1+c3)*u4 &
               +8.0_pr*(c1-c3)*u6+3.0_pr*c3*u8+320.0_pr*c3*(14.0_pr+3.0_pr*u2))*Atan(0.5_pr*u) &
               +2.0_pr*(4800.0_pr*(330.0_pr*c1+129.0_pr*c3-56.0_pr*c4)*u2+48.0_pr*(354354.0_pr*c1+169473.0_pr*c3 &
               +51970.0_pr*c4)*u4+8.0_pr*(-226050.0_pr*c1+177369.0_pr*c3+7676.0_pr*c4)*u6 &
               -12.0_pr*(5874.0_pr*c1-8979.0_pr*c3+472.0_pr*c4)*u8+135.0_pr*(-135.0_pr*c3+32.0_pr*c4)*u10 &
               +48.0_pr*u3*(40425.0_pr*(6.0_pr*c1+3.0_pr*c3-2.0_pr*c4)+891.0_pr*(-2.0_pr*c1+c3)*u4 &
               +88.0_pr*(3.0_pr*c1-3.0_pr*c3-c4)*u6+90.0_pr*(3.0_pr*c3+c4)*u8)*Atan(u) &
               +96.0_pr*(-50.0_pr*(330.0_pr*c1+129.0_pr*c3-56.0_pr*c4)-396.0_pr*(162.0_pr*c1+69.0_pr*c3 &
               -40.0_pr*c4)*u2+891.0_pr*(14.0_pr*c1+11.0_pr*c3-8.0_pr*c4)*u4)*Log(1.0_pr+u2)))) &
               /(1.216512e6_pr*fpi4*Pi*u11)

          ddarhodrho=ddarhodrho+(gA2*mpi*(-2400.0_pr*(330.0_pr*c1+129.0_pr*c3-56.0_pr*c4)*u2 &
               -24.0_pr*(273126.0_pr*c1+129687.0_pr*c3+34280.0_pr*c4)*u4-4.0_pr*(229530.0_pr*c1+264813.0_pr*c3 &
               +56912.0_pr*c4)*u6+16.0_pr*(8052.0_pr*c1-6513.0_pr*c3-188.0_pr*c4)*u8+6.0_pr*(267.0_pr*(c1 &
               -c3)+16.0_pr*c4)*u10+3.0_pr*(4.0_pr+u2)*(315.0_pr*u3*(2240.0_pr*(2.0_pr*c1+c3) &
               +360.0_pr*c3*u2+18.0_pr*(-2.0_pr*c1+c3)*u4+(c1-c3)*u6)*Atan(0.5_pr*u) &
               -4.0_pr*(2.0_pr*u3*(7350.0_pr*(6.0_pr*c1+3.0_pr*c3-2.0_pr*c4)+81.0_pr*(-2.0_pr*c1+c3)*u4 &
               -4.0_pr*(-3.0_pr*c1+3.0_pr*c3+c4)*u6)*Atan(u)+(-50.0_pr*(330.0_pr*c1+129.0_pr*c3-56.0_pr*c4) &
               -324.0_pr*(162.0_pr*c1+69.0_pr*c3-40.0_pr*c4)*u2+567.0_pr*(14.0_pr*c1+11.0_pr*c3 &
               -8.0_pr*c4)*u4)*Log(1.0_pr+u2)))))/(6912.0_pr*fpi4*Pi*u12*(4.0_pr+u2))

          brhodrho=brhodrho+(gA2*mpi*(-99.0_pr*u3*(33600.0_pr*(2.0_pr*c1+c3)*u+560.0_pr*(-10.0_pr*c1 &
               +13.0_pr*c3)*u3+840.0_pr*(-2.0_pr*c1+c3)*u5+54.0_pr*c3*u7-105.0_pr*(4.0_pr+u2) &
               *(8.0_pr*c1*(40.0_pr-10.0_pr*u2+u4)+c3*(160.0_pr+8.0_pr*u2+4.0_pr*u4-3.0_pr*u6)) &
               *Atan(0.5_pr*u))+64.0_pr*c4*(u2*(1680.0_pr-23610.0_pr*u2-958.0_pr*u4+174.0_pr*u6 &
               -27.0_pr*u8)+6.0_pr*u3*(5775.0_pr+44.0_pr*u6+45.0_pr*u8)*Atan(u)+24.0_pr*(-70.0_pr &
               +99.0_pr*u2*(-5.0_pr+3.0_pr*u2))*Log(1.0_pr+u2))))/(1.216512e6_pr*fpi4*Pi*u10)

          dbrhodrho=dbrhodrho+(gA2*mpi*(2.0_pr*u2*(11642400.0_pr*(2.0_pr*c1+c3)*u2+64.0_pr*c4*( &
               -8400.0_pr+77955.0_pr*u2+1919.0_pr*u4-177.0_pr*u6+135.0_pr*u8)-3465.0_pr*u4 &
               *(8.0_pr*c1*(70+3*u2)+c3*(-440.0_pr+9.0_pr*u2*(-4.0_pr+u2))))-3.0_pr*(3465.0_pr*u3 &
               *(8960.0_pr*c1+72.0_pr*(-2.0_pr*c1+c3)*u4+8.0_pr*(c1-c3)*u6+3.0_pr*c3*u8+320.0_pr*c3 &
               *(14.0_pr+3.0_pr*u2))*Atan(0.5_pr*u)-128.0_pr*c4*(u3*(-40425.0_pr-44.0_pr*u6+45.0_pr*u8) &
               *Atan(u)+8.0_pr*(350.0_pr+99.0_pr*u2*(20.0_pr-9.0_pr*u2))*Log(1.0_pr+u2))))) &
               /(1.216512e6_pr*fpi4*Pi*u11)

          ddbrhodrho=ddbrhodrho+(gA2*mpi*(2.0_pr*u2*(32.0_pr*c4*(4.0_pr+u2)*(1050.0_pr-6690.0_pr*u2 &
               -106.0_pr*u4+3.0_pr*u6)+315.0_pr*u2*(-6720.0_pr*(2.0_pr*c1+c3)-40.0_pr*(56.0_pr*c1 &
               +55.0_pr*c3)*u2+2.0_pr*(110.0_pr*c1-89.0_pr*c3)*u4+3.0_pr*(c1-c3)*u6))+3.0_pr*(4.0_pr &
               +u2)*(315.0_pr*u3*(2240.0_pr*(2.0_pr*c1+c3)+360.0_pr*c3*u2+18.0_pr*(-2.0_pr*c1+c3)*u4+(c1 &
               -c3)*u6)*Atan(0.5_pr*u)+64.0_pr*c4*(u3*(3675.0_pr+u6)*Atan(u)+(-350.0_pr+81.0_pr*u2 &
               *(-20.0_pr+7.0_pr*u2))*Log(1.0_pr+u2)))))/(6912.0_pr*fpi4*Pi*u12*(4.0_pr+u2))

          adrdr=adrdr+0.0_pr
          dadrdr=dadrdr+0.0_pr
          bdrdr=bdrdr+0.0_pr
          dbdrdr=dbdrdr+0.0_pr

          arhotau=arhotau+(gA2*mpi*(60.0_pr*(330.0_pr*c1+129.0_pr*c3-56.0_pr*c4)*u2-6.0_pr*(20262.0_pr*c1 &
               +11139.0_pr*c3-7870.0_pr*c4)*u4+(2838.0_pr*c1-2865.0_pr*c3+1916.0_pr*c4)*u6 &
               +12.0_pr*(132.0_pr*c1-29.0_pr*(3.0_pr*c3+c4))*u8+54.0_pr*(3.0_pr*c3+c4)*u10 &
               -6.0_pr*u3*(-5775.0_pr*(6.0_pr*c1+3.0_pr*c3-2.0_pr*c4)+297.0_pr*(2.0_pr*c1-c3)*u4 &
               +88.0_pr*(-3.0_pr*c1+3.0_pr*c3+c4)*u6+90.0_pr*(3.0_pr*c3+c4)*u8)*Atan(u) &
               +6.0_pr*(-1290.0_pr*c3+560.0_pr*c4+33.0_pr*(-100.0_pr*c1-3.0_pr*(162.0_pr*c1+69.0_pr*c3 &
               -40.0_pr*c4)*u2+9.0_pr*(14.0_pr*c1+11.0_pr*c3-8.0_pr*c4)*u4))*Log(1.0_pr+u2))) &
               /(19008.0_pr*fpi4*Pi*u10)

          darhotau=darhotau -(gA2*mpi*(300.0_pr*(330.0_pr*c1+129.0_pr*c3-56.0_pr*c4)*u2 &
               -3.0_pr*(130746.0_pr*c1+73077.0_pr*c3-51970.0_pr*c4)*u4+(8250.0_pr*c1-6603.0_pr*c3 &
               +3838.0_pr*c4)*u6+6.0_pr*(132.0_pr*c1-59.0_pr*(3.0_pr*c3+c4))*u8+270.0_pr*(3.0_pr*c3 &
               +c4)*u10+3.0_pr*u3*(40425.0_pr*(6.0_pr*c1+3.0_pr*c3-2.0_pr*c4)+891.0_pr*(-2.0_pr*c1 &
               +c3)*u4+88.0_pr*(3.0_pr*c1-3.0_pr*c3-c4)*u6+90.0_pr*(3.0_pr*c3+c4)*u8)*Atan(u) &
               +6.0_pr*(-50.0_pr*(330.0_pr*c1+129.0_pr*c3-56.0_pr*c4)-396.0_pr*(162.0_pr*c1+69.0_pr*c3 &
               -40.0_pr*c4)*u2+891.0_pr*(14.0_pr*c1+11.0_pr*c3-8.0_pr*c4)*u4)*Log(1.0_pr+u2))) &
               /(9504.0_pr*fpi4*Pi*u11)

          ddarhotau=ddarhotau+(gA2*mpi*(150.0_pr*(330.0_pr*c1+129.0_pr*c3-56.0_pr*c4)*u2 &
               -3.0_pr*(43962.0_pr*c1+24969.0_pr*c3-17840.0_pr*c4)*u4+(2154.0_pr*c1-1695.0_pr*c3 &
               +848.0_pr*c4)*u6-24.0_pr*(-3.0_pr*c1+3.0_pr*c3+c4)*u8+6.0_pr*u3*(7350.0_pr*(6.0_pr*c1 &
               +3.0_pr*c3-2.0_pr*c4)+81.0_pr*(-2.0_pr*c1+c3)*u4-4.0_pr*(-3.0_pr*c1+3.0_pr*c3+c4)*u6) &
               *Atan(u)+3.0_pr*(-50.0_pr*(330.0_pr*c1+129.0_pr*c3-56.0_pr*c4)-324.0_pr*(162.0_pr*c1+69.0_pr*c3 &
               -40.0_pr*c4)*u2+567.0_pr*(14.0_pr*c1+11.0_pr*c3-8.0_pr*c4)*u4)*Log(1.0_pr+u2))) &
               /(432.0_pr*fpi4*Pi*u12)

          brhotau=brhotau+(c4*gA2*mpi*(u2*(-1680.0_pr+23610.0_pr*u2+958.0_pr*u4-174.0_pr*u6 &
               +27.0_pr*u8)-6.0_pr*u3*(5775.0_pr+44.0_pr*u6+45.0_pr*u8)*Atan(u)+24.0_pr*(70.0_pr &
               +99.0_pr*u2*(5.0_pr-3.0_pr*u2))*Log(1.0_pr+u2)))/(4752.0_pr*fpi4*Pi*u10)

          dbrhotau=dbrhotau+(c4*gA2*mpi*(u2*(8400.0_pr-77955.0_pr*u2-1919.0_pr*u4+177.0_pr*u6 &
               -135.0_pr*u8)+3.0_pr*u3*(40425.0_pr+44.0_pr*u6-45.0_pr*u8)*Atan(u)+24.0_pr*(-350.0_pr &
               +99.0_pr*u2*(-20.0_pr+9.0_pr*u2))*Log(1.0_pr+u2)))/(2376.0_pr*fpi4*Pi*u11)

          ddbrhotau=ddbrhotau+(c4*gA2*mpi*(-1050.0_pr*u2+6690.0_pr*u4+106.0_pr*u6-3.0_pr*u8 &
               -3.0_pr*u3*(3675.0_pr+u6)*Atan(u)+3.0_pr*(350.0_pr+81.0_pr*u2*(20.0_pr-7.0_pr*u2)) &
               *Log(1.0_pr+u2)))/(27.0_pr*fpi4*Pi*u12)

          ajj=ajj+(gA2*mpi*(u2*(-84.0_pr*(c1+2.0_pr*c1*u2)+6.0_pr*c3*(-11.0_pr+16.0_pr*u2+u4) &
               -c4*(12.0_pr+11.0_pr*u2*(2.0_pr+u2)))+2.0_pr*u3*(210.0_pr*c1-105.0_pr*c3 &
               +28.0_pr*(-3.0_pr*c1+3.0_pr*c3+c4)*u2+12.0_pr*(9.0_pr*c3+c4)*u4)*Atan(u)+(84.0_pr*c1 &
               +66.0_pr*c3+12.0_pr*c4+7.0_pr*(-30.0_pr*c1+21.0_pr*c3+4.0_pr*c4)*u2)*Log(1.0_pr+u2))) &
               /(4480.0_pr*fpi4*Pi*u6)

          dajj=dajj+(gA2*mpi*(u3*(315.0_pr*(-2.0_pr*c1+c3)+28.0_pr*(3.0_pr*c1-3.0_pr*c3-c4)*u2 &
               +12.0_pr*(9.0_pr*c3+c4)*u4)*Atan(u)+2.0_pr*(9.0_pr*(14.0_pr*c1+11.0_pr*c3+2.0_pr*c4)*u2 &
               +(42.0_pr*c1-60.0_pr*c3+19.0_pr*c4)*u4+6.0_pr*(9.0_pr*c3+c4)*u6+(-9.0_pr*(14.0_pr*c1 &
               +11.0_pr*c3+2.0_pr*c4)+7.0_pr*(30.0_pr*c1-21.0_pr*c3-4.0_pr*c4)*u2)*Log(1.0_pr+u2)))) &
               /(2240.0_pr*fpi4*Pi*u7)

          bjj=bjj+(c4*gA2*mpi*(-(u2*(12.0_pr+11.0_pr*u2*(2.0_pr+u2)))+8.0_pr*u5*(7.0_pr &
               +3.0_pr*u2)*Atan(u)+4.0_pr*(3.0_pr+7.0_pr*u2)*Log(1.0_pr+u2)))/(2240.0_pr*fpi4*Pi*u6)

          dbjj=dbjj+(c4*gA2*mpi*(u2*(18.0_pr+19.0_pr*u2+6.0_pr*u4)+2.0_pr*u5*(-7.0_pr &
               +3.0_pr*u2)*Atan(u)-2.0_pr*(9.0_pr+14.0_pr*u2)*Log(1.0_pr+u2)))/(560.0_pr*fpi4*Pi*u7)

       Else  !! if (u.gt.ucut...
          arhorho=arhorho +(3.0_pr*(2.0_pr*c1 - c3)*gA2*mpi**3)/(64.0_pr*fpi4*Pi) + (11.0_pr*(42.0_pr*c1 &
               + 129.0_pr*c3 + 40.0_pr*c4)*gA2*mpi**3*u4)/(100800.0_pr*fpi4*Pi) - ((162.0_pr*c1 &
               + 213.0_pr*c3 + 56.0_pr*c4)*gA2*mpi**3*u6)/(69300.0_pr*fpi4*Pi)
          darhorho=darhorho+(gA2*mpi**3*u3*(121.0_pr*(42.0_pr*c1+129.0_pr*c3+40.0_pr*c4) &
               -24.0_pr*(162.0_pr*c1+213.0_pr*c3+56.0_pr*c4)*u2))/(277200.0_pr*fpi4*Pi)
          ddarhorho=ddarhorho+(gA2*mpi**3*u2*(1573.0_pr*(42.0_pr*c1+129.0_pr*c3+40.0_pr*c4) &
               -520.0_pr*(162.0_pr*c1+213.0_pr*c3+56.0_pr*c4)*u2+765.0_pr*(110.0_pr*c1+107.0_pr*c3 &
               +24.0_pr*c4)*u4))/(1.2012e6_pr*fpi4*Pi)

          brhorho=brhorho +(11.0_pr*c4*gA2*mpi**3*u4)/(1260.0_pr*fpi4*Pi) &
               - (4.0_pr*c4*gA2*mpi**3*u6)/(2475.0_pr*fpi4*Pi)
          dbrhorho=dbrhorho+(c4*gA2*mpi**3*u3*(605.0_pr-168.0_pr*u2))/(17325.0_pr*fpi4*Pi)
          ddbrhorho=ddbrhorho+(c4*gA2*mpi**3*u2*(1573.0_pr-728.0_pr*u2+459.0_pr*u4)) &
               /(15015.0_pr*fpi4*Pi)

          arhodrho=arhodrho+(gA2*mpi*(1081080.0_pr*(10.0_pr*c1-11.0_pr*c3+4.0_pr*c4)+16016.0_pr &
               *(42.0_pr*c1+129.0_pr*c3+140.0_pr*c4)*u2-117.0_pr*(5886.0_pr*c1+7739.0_pr*c3+2688.0_pr*c4)*u4 &
               +135.0_pr*(2970.0_pr*c1+2889.0_pr*c3+704.0_pr*c4)*u6))/(5.5351296e8_pr*fpi4*Pi)
          darhodrho=darhodrho+(gA2*mpi*u*(16016.0_pr*(42.0_pr*c1+129.0_pr*c3+140.0_pr*c4)-234.0_pr &
               *(5886.0_pr*c1+7739.0_pr*c3+2688.0_pr*c4)*u2+405.0_pr*(2970.0_pr*c1+2889.0_pr*c3 &
               +704.0_pr*c4)*u4))/(2.7675648e8_pr*fpi4*Pi)
          ddarhodrho=ddarhodrho+(gA2*mpi*(192192.0_pr*(42.0_pr*c1+129.0_pr*c3+140.0_pr*c4)-8424.0_pr &
               *(5886.0_pr*c1+7739.0_pr*c3+2688.0_pr*c4)*u2+24300.0_pr*(2970.0_pr*c1+2889.0_pr*c3 &
               +704.0_pr*c4)*u4-175.0_pr*(442806.0_pr*c1+367383.0_pr*c3+73216.0_pr*c4)*u6)) &
               /(3.32107776e9_pr*fpi4*Pi)

          brhodrho=brhodrho+(gA2*mpi*(2162160.0_pr*(10.0_pr*c1-11.0_pr*c3+4*c4)-40040.0_pr*(42.0_pr*c1 &
               +129.0_pr*c3-112.0_pr*c4)*u2+819.0_pr*(270.0_pr*c1+355.0_pr*c3-768.0_pr*c4)*u4 &
               -45.0_pr*(770.0_pr*c1+749.0_pr*c3-4224.0_pr*c4)*u6))/(5.5351296e8_pr*fpi4*Pi)
          dbrhodrho=dbrhodrho+(gA2*mpi*u*(-40040.0_pr*(42.0_pr*c1+129.0_pr*c3-112.0_pr*c4)+1638.0_pr &
               *(270.0_pr*c1+355.0_pr*c3-768.0_pr*c4)*u2-135.0_pr*(770.0_pr*c1+749.0_pr*c3 &
               -4224.0_pr*c4)*u4))/(2.7675648e8_pr*fpi4*Pi)
          ddbrhodrho=ddbrhodrho+(gA2*mpi*(-480480.0_pr*(42.0_pr*c1+129.0_pr*c3-112.0_pr*c4)+58968.0_pr &
               *(270.0_pr*c1+355.0_pr*c3-768.0_pr*c4)*u2-8100.0_pr*(770.0_pr*c1+749.0_pr*c3 &
               -4224.0_pr*c4)*u4+175.0_pr*(11466.0_pr*c1+9513.0_pr*c3-146432.0_pr*c4)*u6)) &
               /(3.32107776e9_pr*fpi4*Pi)

          adrdr=adrdr+0.0_pr
          dadrdr=dadrdr+0.0_pr
          bdrdr=bdrdr+0.0_pr
          dbdrdr=dbdrdr+0.0_pr

          arhotau=arhotau+(gA2*mpi*(135135.0_pr*(10.0_pr*c1-11.0_pr*c3-4.0_pr*c4)-7007.0_pr*(42.0_pr*c1 &
               +129.0_pr*c3+40.0_pr*c4)*u2+702.0_pr*(162.0_pr*c1+213.0_pr*c3+56.0_pr*c4)*u4-495.0_pr &
               *(110.0_pr*c1+107.0_pr*c3+24.0_pr*c4)*u6))/(1.729728e7_pr*fpi4*Pi)
          darhotau=darhotau+(gA2*mpi*u*(-7007.0_pr*(42.0_pr*c1+129.0_pr*c3+40.0_pr*c4)+1404.0_pr &
               *(162.0_pr*c1+213.0_pr*c3+56.0_pr*c4)*u2-1485.0_pr*(110.0_pr*c1+107.0_pr*c3+24.0_pr*c4)*u4)) &
               /(8.64864e6_pr*fpi4*Pi)
          ddarhotau=ddarhotau+(gA2*mpi*(-21021.0_pr*(42.0_pr*c1+129.0_pr*c3+40.0_pr*c4)+12636.0_pr &
               *(162.0_pr*c1+213.0_pr*c3+56.0_pr*c4)*u2-22275.0_pr*(110.0_pr*c1+107.0_pr*c3+24.0_pr*c4)*u4 &
               +4550.0_pr*(546.0_pr*c1+453.0_pr*c3+88.0_pr*c4)*u6))/(2.594592e7_pr*fpi4*Pi)

          brhotau=brhotau -(c4*gA2*mpi*(135135.0_pr+70070.0_pr*u2-9828.0_pr*u4+2970.0_pr*u6)) &
               /(2.16216e6_pr*fpi4*Pi)
          dbrhotau=dbrhotau -(c4*gA2*mpi*u*(35035.0_pr-9828.0_pr*u2+4455.0_pr*u4))/(540540.0_pr*fpi4*Pi)
          ddbrhotau=ddbrhotau+(c4*gA2*mpi*(-105105.0_pr+88452.0_pr*u2-66825.0_pr*u4+50050.0_pr*u6)) &
               /(1.62162e6_pr*fpi4*Pi)

          ajj=ajj+(gA2*mpi*(-1050.0_pr*(10.0_pr*c1-11.0_pr*c3-2.0_pr*c4)+70.0_pr*(42.0_pr*c1 &
               +129.0_pr*c3+10.0_pr*c4)*u2-3.0_pr*(486.0_pr*c1+639.0_pr*c3+28.0_pr*c4)*u4+8.0_pr &
               *(110.0_pr*c1+107.0_pr*c3+3.0_pr*c4)*u6))/(268800.0_pr*fpi4*Pi)
          dajj=dajj+(gA2*mpi*u*(35.0_pr*(42.0_pr*c1+129.0_pr*c3+10.0_pr*c4)-3.0_pr*(486.0_pr*c1 &
               +639.0_pr*c3+28.0_pr*c4)*u2+12.0_pr*(110.0_pr*c1+107.0_pr*c3+3.0_pr*c4)*u4)) &
               /(67200.0_pr*fpi4*Pi)

          bjj=bjj+(c4*gA2*mpi*(525.0_pr+175.0_pr*u2-21.0_pr*u4+6.0_pr*u6))/(33600.0_pr*fpi4*Pi)
          dbjj=dbjj+(c4*gA2*mpi*u*(175.0_pr-42.0_pr*u2+18.0_pr*u4))/(16800.0_pr*fpi4*Pi)
       End If  !! if (u.gt.ucut...
       !
       ! N2LO, 3N terms
       If (use_DME3N_terms) Then
          If (u.Gt.ucut3n) Then
             y=1.0_pr/u
             dy=-y*y
             ddy=2.0_pr*y*y*y
             marc=atu
             dmarc=2.0_pr/(1.0_pr+ 4.0_pr*u2)
             ddmarc=-16.0_pr*u/((1.0_pr+ 4.0_pr*u2)**2)
             mlog=lual
             dmlog=8.0_pr*u/(1.0_pr+ 4.0_pr*u2)
             ddmlog=(8.0_pr-32.0_pr*u2)/((1.0_pr+ 4.0_pr*u2)**2)
             Do l=1,33
                Do j=1,3
                   Do k=1,3
                      aux=((y**(l-1))*(atu**(j-1))*(lual**(k-1)))*mevfm
                      daux=(Real(-1+j,pr)*marc**(-2+j)*mlog**(-1+k)*y**(-1+l)*dmarc &
                           +Real(-1+k,pr)*marc**(-1+j)*mlog**(-2+k)*y**(-1+l)*dmlog &
                           +Real(-1+l,pr)*marc**(-1+j)*mlog**(-1+k)*y**(-2+l)*dy)*mevfm
                      ddaux=((2.0_pr*Real(-1+l,pr)*y**(-2+l)*(Real(-1+j,pr) &
                           *marc**(-2+j)*mlog**(-1+k)*dmarc+Real(-1+k,pr)*marc**(-1+j)*mlog**(-2+k) &
                           *dmlog)*dy+y**(-1+l)*(2.0_pr*Real((-1+j)*(-1+k),pr)*marc**(-2+j)*mlog**(-2+k) &
                           *dmarc*dmlog+mlog**(-1+k)*(Real((-2+j)*(-1+j),pr)*marc**(-3+j)*dmarc**2 &
                           +Real(-1+j,pr)*marc**(-2+j)*ddmarc)+marc**(-1+j)*(Real((-2+k)*(-1+k),pr) &
                           *mlog**(-3+k)*dmlog**2+Real(-1+k,pr)*mlog**(-2+k)*ddmlog))+marc**(-1+j) &
                           *mlog**(-1+k)*(Real((-2+l)*(-1+l),pr)*y**(-3+l)*dy**2+Real(-1+l,pr) &
                           *y**(-2+l)*ddy)))*mevfm
                      !
                      hrho0rho0 =hrho0rho0 +ctr0r0(k,j,l)*aux
                      hrho1rho1 =hrho1rho1 +ctr1r1(k,j,l)*aux
                      hdr0dr0   =hdr0dr0+ctdr0dr0(k,j,l)*aux
                      hdr1dr1   =hdr1dr1+ctdr1dr1(k,j,l)*aux
                      hrho0Drho0=hrho0Drho0+ctr0Dr0(k,j,l)*aux
                      hrho1Drho0=hrho1Drho0+ctr1Dr0(k,j,l)*aux
                      hrho1Drho1=hrho1Drho1+ctr1Dr1(k,j,l)*aux
                      hrho0tau0 =hrho0tau0+ctr0t0(k,j,l)*aux
                      hrho1tau0 =hrho1tau0+ctr1t0(k,j,l)*aux
                      hrho1tau1 =hrho1tau1+ctr1t1(k,j,l)*aux
                      hJ0dr0    =hJ0dr0+ctJ0dr0(k,j,l)*aux
                      hrho0DJ0  =hrho0DJ0+ctr0dJ0(k,j,l)*aux
                      hJ1dr1    =hJ1dr1+ctJ1dr1(k,j,l)*aux
                      hrho1DJ1  =hrho1DJ1+ctr1dJ1(k,j,l)*aux
                      hJ0dr1    =hJ0dr1+ctJ0dr1(k,j,l)*aux
                      hrho1DJ0  =hrho1DJ0+ctr1dJ0(k,j,l)*aux
                      hJ1dr0    =hJ1dr0+ctJ1dr0(k,j,l)*aux
                      hJ0J0     =hJ0J0+ctJ0J0(k,j,l)*aux
                      hJ0J1     =hJ0J1+ctJ0J1(k,j,l)*aux
                      hJ1J1     =hJ1J1+ctJ1J1(k,j,l)*aux
                      !
                      dhrho0rho0 =dhrho0rho0 +ctr0r0(k,j,l)*daux
                      dhrho1rho1 =dhrho1rho1 +ctr1r1(k,j,l)*daux
                      dhdr0dr0   =dhdr0dr0+ctdr0dr0(k,j,l)*daux
                      dhdr1dr1   =dhdr1dr1+ctdr1dr1(k,j,l)*daux
                      dhrho0Drho0=dhrho0Drho0+ctr0Dr0(k,j,l)*daux
                      dhrho1Drho0=dhrho1Drho0+ctr1Dr0(k,j,l)*daux
                      dhrho1Drho1=dhrho1Drho1+ctr1Dr1(k,j,l)*daux
                      dhrho0tau0 =dhrho0tau0+ctr0t0(k,j,l)*daux
                      dhrho1tau0 =dhrho1tau0+ctr1t0(k,j,l)*daux
                      dhrho1tau1 =dhrho1tau1+ctr1t1(k,j,l)*daux
                      dhJ0dr0    =dhJ0dr0+ctJ0dr0(k,j,l)*daux
                      dhrho0DJ0  =dhrho0DJ0+ctr0dJ0(k,j,l)*daux
                      dhJ1dr1    =dhJ1dr1+ctJ1dr1(k,j,l)*daux
                      dhrho1DJ1  =dhrho1DJ1+ctr1dJ1(k,j,l)*daux
                      dhJ0dr1    =dhJ0dr1+ctJ0dr1(k,j,l)*daux
                      dhrho1DJ0  =dhrho1DJ0+ctr1dJ0(k,j,l)*daux
                      dhJ1dr0    =dhJ1dr0+ctJ1dr0(k,j,l)*daux
                      dhJ0J0     =dhJ0J0+ctJ0J0(k,j,l)*daux
                      dhJ0J1     =dhJ0J1+ctJ0J1(k,j,l)*daux
                      dhJ1J1     =dhJ1J1+ctJ1J1(k,j,l)*daux
                      !
                      ddhrho0rho0 =ddhrho0rho0 +ctr0r0(k,j,l)*ddaux
                      ddhrho1rho1 =ddhrho1rho1 +ctr1r1(k,j,l)*ddaux
                      ddhrho0Drho0=ddhrho0Drho0+ctr0Dr0(k,j,l)*ddaux
                      ddhrho1Drho0=ddhrho1Drho0+ctr1Dr0(k,j,l)*ddaux
                      ddhrho1Drho1=ddhrho1Drho1+ctr1Dr1(k,j,l)*ddaux
                      ddhrho0tau0 =ddhrho0tau0+ctr0t0(k,j,l)*ddaux
                      ddhrho1tau0 =ddhrho1tau0+ctr1t0(k,j,l)*ddaux
                      ddhrho1tau1 =ddhrho1tau1+ctr1t1(k,j,l)*ddaux
                   End Do
                End Do
             End Do  !! end of l,j,k loops
          Else  !! if (u.gt.ucut3n..
             hrho0rho0=mevfm*(2.0555896e-8_pr*ce*(-15.058202_pr + u)*(13.851257_pr + u)*(218.02762_pr &
                  - 15.661674_pr*u + u2)*(200.58157_pr + 13.247785_pr*u + u2) + gA*(-0.042893802_pr*cd*(-1.5948724_pr &
                  + u)*(-0.01657426_pr + u)*(1.8884833_pr - 2.1431293_pr*u + u2)*(0.029641853_pr - 0.28437579_pr*u &
                  + u2) + gA*LambdaX*(6.0619551_pr*c3*(-1.1935357_pr + u)*(-0.30214838_pr + u)*(0.9281019_pr &
                  - 1.5964577_pr*u + u2)*(0.15916821_pr - 0.63445995_pr*u + u2) + 3.9928539_pr*c1*(-1.2431741_pr &
                  + u)*(-0.19155888_pr + u)*(1.1456958_pr - 1.8138841_pr*u + u2)*(0.19109424_pr - 0.36771486_pr*u &
                  + u2) - 0.06402411_pr*c4*(-1.3982905_pr + u)*(0.1480126_pr + u)*(1.6779335_pr - 2.2114154_pr*u &
                  + u2)*(0.041603529_pr - 0.35726143_pr*u + u2))))/(fpi4*LambdaX)

             hrho1rho1=mevfm*(2.0555896e-8_pr*ce*(-15.058202_pr + u)*(13.851257_pr + u)*(218.02762_pr &
                  - 15.661674_pr*u + u2)*(200.58157_pr + 13.247785_pr*u + u2) + gA*(-0.042893802_pr*cd*(-1.5948724_pr &
                  + u)*(-0.01657426_pr + u)*(1.8884833_pr - 2.1431293_pr*u + u2)*(0.029641853_pr - 0.28437579_pr*u &
                  + u2) + gA*LambdaX*(6.0619551_pr*c3*(-1.1935357_pr + u)*(-0.30214838_pr + u)*(0.9281019_pr &
                  - 1.5964577_pr*u + u2)*(0.15916821_pr - 0.63445995_pr*u + u2) + 3.9928539_pr*c1*(-1.2431741_pr &
                  + u)*(-0.19155888_pr + u)*(1.1456958_pr - 1.8138841_pr*u + u2)*(0.19109424_pr - 0.36771486_pr*u &
                  + u2) - 0.06402411_pr*c4*(-1.3982905_pr + u)*(0.1480126_pr + u)*(1.6779335_pr - 2.2114154_pr*u &
                  + u2)*(0.041603529_pr - 0.35726143_pr*u + u2))))/(fpi4*LambdaX)

             hdr0dr0=0.0_pr
             hdr1dr1=0.0_pr

             hrho0Drho0=mevfm*(gA*(0.12666864_pr*cd*(-1.2705742_pr + u)*(0.21652154_pr + u)*(1.4119605_pr &
                  - 2.1693922_pr*u + u2)*(0.89941341_pr - 1.1919922_pr*u + u2) + gA*LambdaX*(29.501949_pr*c1*( &
                  -1.0215543_pr + u)*(-0.15813736_pr + u)*(0.84479502_pr - 1.6998479_pr*u + u2)*(0.34169853_pr &
                  - 0.81578924_pr*u + u2) - 7.1108952_pr*c3*(1.1461026_pr - 2.0823268_pr*u + u2)*(0.62212787_pr &
                  - 1.1642433_pr*u + u2)*(0.074190685_pr - 0.4633705_pr*u + u2) + 0.18892679_pr*c4*(1.7714304_pr &
                  - 2.5914415_pr*u + u2)*(1.2438813_pr - 1.6683493_pr*u + u2)*(0.004127797_pr - 0.089415137_pr*u &
                  + u2))))/(fpi4*LambdaX*mpi2)

             hrho1Drho0=mevfm*(gA*(-0.042222876_pr*cd*(-1.2705742_pr + u)*(0.21652155_pr + u)*(1.4119605_pr &
                  - 2.1693922_pr*u + u2)*(0.89941343_pr - 1.1919922_pr*u + u2) + gA*LambdaX*(-0.064338673_pr*c1*( &
                  -1.9302829_pr + u)*(0.24441857_pr + u)*(2.1104361_pr - 2.7720023_pr*u + u2)*(1.3417224_pr &
                  - 1.5814849_pr*u + u2) + 0.049739131_pr*c3*(1.8967607_pr - 2.6742154_pr*u + u2)*(1.316844_pr &
                  - 1.6740809_pr*u + u2)*(0.011061538_pr - 0.13317506_pr*u + u2) - 0.073720916_pr*c4*(1.6783675_pr &
                  - 2.5203342_pr*u + u2)*(1.1903356_pr - 1.6163511_pr*u + u2)*(0.0074741816_pr - 0.10943007_pr*u &
                  + u2))))/(fpi4*LambdaX*mpi2)

             hrho1Drho1=mevfm*(gA*(-0.084445761_pr*cd*(-1.2705742_pr + u)*(0.21652154_pr + u)*(1.4119605_pr &
                  - 2.1693922_pr*u + u2)*(0.89941341_pr - 1.1919922_pr*u + u2) + gA*LambdaX*(-7.0472043_pr*c1*( &
                  -1.0402373_pr + u)*(-0.15610243_pr + u)*(0.87447681_pr - 1.7182557_pr*u + u2)*(0.33211898_pr &
                  - 0.74186157_pr*u + u2) + 8.0775475_pr*c3*(0.95424896_pr - 1.9113384_pr*u + u2)*(0.50394896_pr &
                  - 1.1611245_pr*u + u2)*(0.11279874_pr - 0.61804716_pr*u + u2) - 0.1259512_pr*c4*(1.7714304_pr &
                  - 2.5914415_pr*u + u2)*(1.2438813_pr - 1.6683493_pr*u + u2)*(0.004127797_pr - 0.089415137_pr*u &
                  + u2))))/(fpi4*LambdaX*mpi2)

             hrho0tau0=mevfm*(ga*(-0.25333726_pr*cd*(-1.2705742_pr+u)*(0.21652155_pr+u)*(1.4119605_pr &
                  +(-2.1693922_pr+u)*u)*(0.89941343_pr+(-1.1919922_pr+u)*u)+ga*LambdaX*(-4.0039062_pr*c1 &
                  *(-1.2352256_pr+u)*(0.064279963_pr+u)*(1.2516998_pr+(-1.9866561_pr+u)*u)*(0.62481483_pr &
                  +(-0.6277641_pr+u)*u)-0.43884277_pr*c4*(1.682325_pr+(-2.523347_pr+u)*u)*(1.1930406_pr &
                  +(-1.6187124_pr+u)*u)*(0.0072756559_pr+(-0.10828832_pr+u)*u)-1.3046875_pr*c3*(-1.3097655_pr &
                  +u)*(-0.12575226_pr+u)*(1.337577_pr+(-1.9544714_pr+u)*u)*(0.15931999_pr+u*(0.30346213_pr+u))))) &
                  /(fpi4*LambdaX*mpi2)

             hrho1tau0=mevfm*(ga*(0.084445755_pr*cd*(-1.2705742_pr+u)*(0.21652155_pr+u)*(1.4119605_pr &
                  +(-2.1693922_pr+u)*u)*(0.89941343_pr+(-1.1919922_pr+u)*u)+ga*LambdaX*(0.12695312_pr*c1 &
                  *(-1.9629389_pr+u)*(0.24499103_pr+u)*(2.1016515_pr+(-2.7692262_pr+u)*u)*(1.3407016_pr &
                  +(-1.5830183_pr+u)*u)-0.10644531_pr*c3*(1.8516394_pr+(-2.6412277_pr+u)*u)*(1.2914728_pr &
                  +(-1.6504171_pr+u)*u)*(0.012431285_pr+(-0.14092406_pr+u)*u)+0.14648438_pr*c4*(1.6816422_pr &
                  +(-2.5228331_pr+u)*u)*(1.1925378_pr+(-1.6182745_pr+u)*u)*(0.007307093_pr+(-0.10847576_pr &
                  +u)*u))))/(fpi4*LambdaX*mpi2)

             hrho1tau1=mevfm*(ga*(0.16889151_pr*cd*(-1.2705742_pr+u)*(0.21652155_pr+u)*(1.4119605_pr &
                  +(-2.1693922_pr+u)*u)*(0.89941343_pr+(-1.1919922_pr+u)*u)+ga*LambdaX*(1.0390625_pr*c1 &
                  *(-1.2278564_pr+u)*(0.017559313_pr+u)*(1.2287066_pr+(-1.9487414_pr+u)*u)*(0.66412923_pr &
                  +(-0.34133743_pr+u)*u)-0.19921875_pr*c3*(1.86082_pr+(-2.5829231_pr+u)*u)*(1.9047944_pr &
                  +(-1.441557_pr+u)*u)*(0.025980569_pr+(-0.2323826_pr+u)*u)+0.29296875_pr*c4*(1.6816422_pr &
                  +(-2.5228331_pr+u)*u)*(1.1925378_pr+(-1.6182745_pr+u)*u)*(0.007307093_pr+(-0.10847576_pr &
                  +u)*u))))/(fpi4*LambdaX*mpi2)

             hJ0dr0=0.0_pr

             hrho0DJ0=mevfm*(gA2*(0.87488139_pr*c1*(-1.5075509_pr+u)*(0.26375366_pr+u)*(1.6807042_pr &
                  +(-2.4330383_pr+u)*u)*(0.96407512_pr+(-1.3041604_pr+u)*u)+0.10945568_pr*c4*(1.6868327_pr &
                  +(-2.5235789_pr+u)*u)*(1.16087_pr+(-1.5798009_pr+u)*u)*(0.0058674843_pr+(-0.099090543_pr+u)*u) &
                  -0.65769729_pr*c3*(1.7722509_pr+(-2.5823827_pr+u)*u)*(1.1921659_pr+(-1.56672_pr+u)*u) &
                  *(0.0031379608_pr+(-0.079284979_pr+u)*u)))/(fpi4*mpi2)

             hJ1dr1=0.0_pr

             hrho1DJ1=mevfm*(gA2*(-0.43565926_pr*c1*(-1.2316322_pr+u)*(0.20486088_pr+u)*(1.3126963_pr &
                  +(-2.0776323_pr+u)*u)*(0.80436474_pr+(-1.0724984_pr+u)*u)-0.07297045_pr*c4*(1.6868327_pr &
                  +(-2.5235789_pr+u)*u)*(1.1608701_pr+(-1.5798009_pr+u)*u)*(0.0058674831_pr+(-0.099090535_pr &
                  +u)*u)+0.1942317_pr*c3*(1.7488022_pr+(-2.5635214_pr+u)*u)*(1.1781067_pr+(-1.5454099_pr+u)*u) &
                  *(0.0037519521_pr+(-0.082717466_pr+u)*u)))/(fpi4*mpi2)

             hJ0dr1=0.0_pr

             hrho1DJ0=mevfm*(gA2*(-0.49880616_pr*c1*(-1.2658085_pr+u)*(0.20574562_pr+u)*(1.4174942_pr &
                  +(-2.1861549_pr+u)*u)*(0.9423041_pr+(-1.257487_pr+u)*u)-0.036485225_pr*c4*(1.6868327_pr &
                  +(-2.5235789_pr+u)*u)*(1.1608701_pr+(-1.5798009_pr+u)*u)*(0.0058674831_pr+(-0.099090536_pr &
                  +u)*u)+0.10945568_pr*c3*(1.6868327_pr+(-2.5235789_pr+u)*u)*(1.1608701_pr+(-1.5798009_pr+u)*u) &
                  *(0.0058674829_pr+(-0.099090534_pr+u)*u)))/(fpi4*mpi2)

             hJ1dr0=0.0_pr

             hJ0J0=mevfm*(ga*(-0.1536163_pr*cd*(-1.3967122_pr+u)*(0.25293929_pr+u)*(1.6675202_pr+(-2.3021441_pr &
                  +u)*u)*(0.98775885_pr+(-1.1088032_pr+u)*u)+ga*LambdaX*(0.37890625_pr*c1*(-1.9915809_pr+u) &
                  *(-1.0930869_pr+u)*(-0.66138682_pr+u)*(0.62037829_pr+u)*(0.88993989_pr+(-1.6302515_pr+u)*u) &
                  -0.82470703_pr*c3*(-0.1151617_pr+u)*(0.19519195_pr+u)*(1.3633551_pr+(-2.2799728_pr+u)*u) &
                  *(0.92255094_pr+(-1.4901995_pr+u)*u)-0.31591797_pr*c4*(1.6836216_pr+(-2.5209696_pr+u)*u) &
                  *(1.1565142_pr+(-1.5751323_pr+u)*u)*(0.0053910501_pr+(-0.095764354_pr+u)*u)))) &
                  /(fpi4*LambdaX*mpi2)

             hJ0J1=mevfm*(ga*(0.10241087_pr*cd*(-1.3967122_pr+u)*(0.25293929_pr+u)*(1.6675202_pr+(-2.3021441_pr &
                  +u)*u)*(0.98775885_pr+(-1.1088032_pr+u)*u)+ga*LambdaX*(0.19000244_pr*c3*(-0.078125417_pr+u) &
                  *(0.032069683_pr+u)*(1.5685634_pr+(-2.4365485_pr+u)*u)*(1.0625904_pr+(-1.537248_pr+u)*u) &
                  -0.69067383_pr*c1*(-1.1098034_pr+u)*(0.18815441_pr+u)*(1.1047942_pr+(-1.9744345_pr+u)*u) &
                  *(0.7614037_pr+(-1.2677199_pr+u)*u)+0.21069336_pr*c4*(1.6833299_pr+(-2.5207351_pr+u)*u) &
                  *(1.1564102_pr+(-1.5750045_pr+u)*u)*(0.0053966658_pr+(-0.095801508_pr+u)*u)))) &
                  /(fpi4*LambdaX*mpi2)

             hJ1J1=mevfm*(ga*(0.051205435_pr*cd*(-1.3967122_pr+u)*(0.25293929_pr+u)*(1.6675202_pr+(-2.3021441_pr &
                  +u)*u)*(0.98775885_pr+(-1.1088032_pr+u)*u)+ga*LambdaX*(0.19219971_pr*c3*(-0.12213139_pr+u) &
                  *(0.2951145_pr+u)*(2.0521898_pr+(-2.8194623_pr+u)*u)*(1.2084447_pr+(-1.7226538_pr+u)*u) &
                  -0.54418945_pr*c1*(-1.7111589_pr+u)*(0.30206023_pr+u)*(1.421294_pr+(-2.2931824_pr+u)*u) &
                  *(0.85528206_pr+(-1.3180868_pr+u)*u)+0.10534668_pr*c4*(1.6833299_pr+(-2.5207351_pr+u)*u) &
                  *(1.1564102_pr+(-1.5750045_pr+u)*u)*(0.0053966658_pr+(-0.095801508_pr+u)*u)))) &
                  /(fpi4*LambdaX*mpi2)

             dhrho0rho0=mevfm*(-3.0696805e-6_pr*ce*(-0.63337387_pr + u)*(-0.6_pr + u)*(0.37421487_pr &
                  - 1.2218226_pr*u + u2)*(0.33405954_pr - 1.1537321_pr*u + u2) + gA*(0.27912917_pr*cd*(&
                  1.7611632_pr - 2.5821217_pr*u + u2)*(1.2266141_pr - 1.6507188_pr*u + u2)*(0.029143369_pr &
                  - 0.28917996_pr*u + u2) + gA*LambdaX*(-1491.7558_pr*c1*(0.67170782_pr - 1.62489_pr*u &
                  + u2)*(0.4525586_pr - 1.2543926_pr*u + u2)*(0.17617549_pr - 0.80139711_pr*u + u2) &
                  - 536.43232_pr*c3*(0.76316673_pr - 1.7218477_pr*u + u2)*(0.43361641_pr - 1.1835313_pr*u &
                  + u2)*(0.16367811_pr - 0.78054496_pr*u + u2) + 0.30876998_pr*c4*(-1.5784873_pr + u)*( &
                  -1.0937068_pr + u)*(1.340452_pr - 2.04474_pr*u + u2)*(0.019873806_pr - 0.25332432_pr*u &
                  + u2))))/(fpi4*LambdaX)

             dhrho1rho1=mevfm*(-3.0696805e-6_pr*ce*(-0.63337387_pr + u)*(-0.6_pr + u)*(0.37421487_pr &
                  - 1.2218226_pr*u + u2)*(0.33405954_pr - 1.1537321_pr*u + u2) + gA*(0.27912917_pr*cd*(&
                  1.7611632_pr - 2.5821217_pr*u + u2)*(1.2266141_pr - 1.6507188_pr*u + u2)*(0.029143369_pr &
                  - 0.28917996_pr*u + u2) + gA*LambdaX*(-1491.7558_pr*c1*(0.67170782_pr - 1.62489_pr*u &
                  + u2)*(0.4525586_pr - 1.2543926_pr*u + u2)*(0.17617549_pr - 0.80139711_pr*u + u2) &
                  - 536.43232_pr*c3*(0.76316673_pr - 1.7218477_pr*u + u2)*(0.43361641_pr - 1.1835313_pr*u &
                  + u2)*(0.16367811_pr - 0.78054496_pr*u + u2) + 0.30876998_pr*c4*(-1.5784873_pr + u)*( &
                  -1.0937068_pr + u)*(1.340452_pr - 2.04474_pr*u + u2)*(0.019873806_pr - 0.25332432_pr*u &
                  + u2))))/(fpi4*LambdaX)

             dhdr0dr0=0.0_pr
             dhdr1dr1=0.0_pr

             dhrho0Drho0=mevfm*(gA*(-0.52109887_pr*cd*(-1.3082043_pr + u)*(-0.071686644_pr + u)*(1.4677685_pr &
                  - 2.2766738_pr*u + u2)*(0.93721725_pr - 1.4019146_pr*u + u2) + gA*LambdaX*(-1.0505013_pr*c4*( &
                  -0.6756734_pr + u)*(-0.093758568_pr + u)*(1.59727_pr - 2.4899591_pr*u + u2)*(1.0589336_pr &
                  - 1.6598738_pr*u + u2) + 410.5668_pr*c3*(-0.82531999_pr + u)*(-0.33506221_pr + u)*(0.60822284_pr &
                  - 1.5136238_pr*u + u2)*(0.31121003_pr - 0.98997717_pr*u + u2) - 1232.118_pr*c1*(-0.8439295_pr &
                  + u)*(-0.39116927_pr + u)*(0.58509589_pr - 1.4743171_pr*u + u2)*(0.28266076_pr - 0.97625946_pr*u &
                  + u2))))/(fpi4*LambdaX*mpi2)

             dhrho1Drho0=mevfm*(gA*(0.17369909_pr*cd*(-1.308206_pr + u)*(-0.071686588_pr + u)*(1.4677702_pr &
                  - 2.2766756_pr*u + u2)*(0.9372174_pr - 1.4019153_pr*u + u2) + gA*LambdaX*(0.0083633968_pr*c3*( &
                  -1.1273967_pr + u)*(-0.75552081_pr + u)*(-0.07564609_pr + u)*(36.068809_pr + u)*(1.033597_pr &
                  - 1.7618699_pr*u + u2) + 0.36997868_pr*c4*(-0.67567294_pr + u)*(-0.091577817_pr + u)*(1.4806242_pr &
                  - 2.3819186_pr*u + u2)*(1.075619_pr - 1.644617_pr*u + u2) - 0.68655766_pr*c1*(-0.063343379_pr &
                  + u)*(0.97369203_pr + u)*(1.3870134_pr - 2.3075684_pr*u + u2)*(1.0431426_pr - 1.6344504_pr*u &
                  + u2))))/(fpi4*LambdaX*mpi2)

             dhrho1Drho1=mevfm*(gA*(0.34739925_pr*cd*(-1.3082043_pr + u)*(-0.071686644_pr + u)*(1.4677685_pr &
                  - 2.2766738_pr*u + u2)*(0.93721725_pr - 1.4019146_pr*u + u2) + gA*LambdaX*(0.7003342_pr*c4*( &
                  -0.6756734_pr + u)*(-0.093758568_pr + u)*(1.59727_pr - 2.4899591_pr*u + u2)*(1.0589336_pr &
                  - 1.6598738_pr*u + u2) - 352.88094_pr*c3*(-0.80527266_pr + u)*(-0.38439748_pr + u)*(0.56781414_pr &
                  - 1.468113_pr*u + u2)*(0.30457752_pr - 1.0189296_pr*u + u2) + 388.57352_pr*c1*(-0.83954784_pr &
                  + u)*(-0.4308896_pr + u)*(0.57229209_pr - 1.4557589_pr*u + u2)*(0.25241536_pr - 0.93218873_pr*u &
                  + u2))))/(fpi4*LambdaX*mpi2)

             dhrho0tau0=mevfm*(ga*(1.0421949_pr*cd*(-1.3082058_pr+u)*(-0.071686595_pr+u)*(1.46777_pr &
                  +(-2.2766754_pr+u)*u)*(0.93721738_pr+(-1.4019152_pr+u)*u)+ga*LambdaX*(2.1171875_pr*c4 &
                  *(-0.6756739_pr+u)*(-0.090082162_pr+u)*(1.5215283_pr+(-2.4175155_pr+u)*u)*(1.0854982_pr &
                  +(-1.6625771_pr+u)*u)+224.0_pr*c1*(-0.94326159_pr+u)*(-0.29152041_pr+u)*(0.71677264_pr &
                  +(-1.5877363_pr+u)*u)*(0.25328282_pr+(-0.77748171_pr+u)*u)+48.0_pr*c3*(0.88506991_pr &
                  +(-1.8677526_pr+u)*u)*(0.70765726_pr+(-1.4865692_pr+u)*u)*(0.079022134_pr+(-0.37067813_pr &
                  +u)*u))))/(fpi4*LambdaX*mpi2)

             dhrho1tau0=mevfm*(ga*(-0.34739831_pr*cd*(-1.3082058_pr+u)*(-0.071686595_pr+u)*(1.46777_pr &
                  +(-2.2766754_pr+u)*u)*(0.93721738_pr+(-1.4019152_pr+u)*u)+ga*LambdaX*(0.3046875_pr*c3 &
                  *(-0.75536572_pr+u)*(-0.092756738_pr+u)*(2.1347117_pr+(-2.8953677_pr+u)*u)*(1.2457416_pr &
                  +(-1.8629201_pr+u)*u)-0.6953125_pr*c4*(-0.67567303_pr+u)*(-0.089866675_pr+u)*(1.5384947_pr &
                  +(-2.4324091_pr+u)*u)*(1.0870855_pr+(-1.6682029_pr+u)*u)+1.65625_pr*c1*(-0.058760606_pr+u) &
                  *(0.81179488_pr+u)*(1.3491153_pr+(-2.2757768_pr+u)*u)*(1.0158947_pr+(-1.6126348_pr+u)*u)))) &
                  /(fpi4*LambdaX*mpi2)

             dhrho1tau1=mevfm*(ga*(-0.69479662_pr*cd*(-1.3082058_pr+u)*(-0.071686595_pr+u)*(1.46777_pr &
                  +(-2.2766754_pr+u)*u)*(0.93721738_pr+(-1.4019152_pr+u)*u)+ga*LambdaX*(-1.390625_pr*c4 &
                  *(-0.67567303_pr+u)*(-0.089866675_pr+u)*(1.5384947_pr+(-2.4324091_pr+u)*u)*(1.0870855_pr &
                  +(-1.6682029_pr+u)*u)-12.0_pr*c3*(-0.2307474_pr+u)*(0.76898229_pr+(-1.7297015_pr+u)*u) &
                  *(0.4881705_pr+(-1.0916344_pr+u)*u)-192.0_pr*c1*(-0.89170462_pr+u)*(-0.40511316_pr+u) &
                  *(0.6397114_pr+(-1.5139956_pr+u)*u)*(0.21220628_pr+(-0.78918661_pr+u)*u)))) &
                  /(fpi4*LambdaX*mpi2)

             dhJ0dr0=0.0_pr

             dhrho0DJ0=mevfm*(gA2*(-0.59709561_pr*c4*(-0.69515033_pr+u)*(-0.089020409_pr+u)*(1.4150964_pr &
                  +(-2.3274531_pr+u)*u)*(1.0076893_pr+(-1.5882571_pr+u)*u)+4.125_pr*c3*(-0.74030576_pr+u) &
                  *(-0.084829197_pr+u)*(1.3684104_pr+(-2.2855723_pr+u)*u)*(0.95092906_pr+(-1.51202_pr+u)*u) &
                  -20.0_pr*c1*(-1.0733108_pr+u)*(-0.097980712_pr+u)*(0.99914461_pr+(-1.8457061_pr+u)*u) &
                  *(0.66726597_pr+(-1.0330025_pr+u)*u)))/(fpi4*mpi2)

             dhJ1dr1=0.0_pr

             dhrho1DJ1=mevfm*(gA2*(0.39806361_pr*c4*(-0.69515033_pr+u)*(-0.089020402_pr+u)*(1.4150967_pr &
                  +(-2.3274533_pr+u)*u)*(1.0076893_pr+(-1.5882572_pr+u)*u)- 1.125_pr*c3*(-0.74849274_pr+u) &
                  *(-0.080617806_pr+u)*(1.3642178_pr+(-2.2788157_pr+u)*u)*(0.96130355_pr+(-1.5129071_pr+u)*u) &
                  -5.0_pr*c1*(1.0009125_pr+(-1.9641805_pr+u)*u)*(0.6986326_pr+(-1.3789259_pr+u)*u)*(0.041493177_pr &
                  +u*(0.16810641_pr+u))))/(fpi4*mpi2)

             dhJ0dr1=0.0_pr

             dhrho1DJ0=mevfm*(gA2*(0.19903181_pr*c4*(-0.69515033_pr+u)*(-0.089020402_pr+u)*(1.4150967_pr &
                  +(-2.3274533_pr+u)*u)*(1.0076893_pr+(-1.5882572_pr+u)*u)-0.59709555_pr*c3*(-0.69515033_pr+u) &
                  *(-0.089020406_pr+u)*(1.4150964_pr+(-2.3274531_pr+u)*u)*(1.0076893_pr+(-1.5882571_pr+u)*u) &
                  +1.4395216_pr*c1*(-1.6910704_pr+u)*(-0.073466262_pr+u)*(1.5355865_pr+(-2.3975331_pr+u)*u) &
                  *(0.99872048_pr+(-1.5169799_pr+u)*u)))/(fpi4*mpi2)

             dhJ1dr0=0.0_pr

             dhJ0J0=mevfm*(ga*(0.45817098_pr*cd*(-1.6176653_pr+u)*(-0.066866022_pr+u)*(1.7318817_pr &
                  +(-2.5016073_pr+u)*u)*(1.0139312_pr+(-1.425551_pr+u)*u)+ga*LambdaX*(1.7695313_pr*c4 &
                  *(-0.69501265_pr+u)*(-0.088068636_pr+u)*(1.3968124_pr+(-2.3120281_pr+u)*u)*(0.99464383_pr &
                  +(-1.5741511_pr+u)*u)+28.0_pr*c3*(-0.52322774_pr+u)*(-0.14878162_pr+u)*(0.83449184_pr &
                  +(-1.7952477_pr+u)*u)*(0.5636391_pr+(-1.2398858_pr+u)*u)-48.0_pr*c1*(-0.86033423_pr+u) &
                  *(-0.098046912_pr+u)*(0.71197564_pr+(-1.6379394_pr+u)*u)*(0.50452901_pr+(-1.1286795_pr &
                  +u)*u))))/(fpi4*LambdaX*mpi2)

             dhJ0J1=mevfm*(ga*(-0.30544732_pr*cd*(-1.6176653_pr+u)*(-0.066866022_pr+u)*(1.7318817_pr &
                  +(-2.5016073_pr+u)*u)*(1.0139312_pr+(-1.425551_pr+u)*u)+ga*LambdaX*(-1.1796875_pr*c4 &
                  *(-0.69501214_pr+u)*(-0.088062144_pr+u)*(1.3949349_pr+(-2.3100216_pr+u)*u)*(0.99576969_pr &
                  +(-1.5747849_pr+u)*u)-1.5332031_pr*c3*(-0.63533554_pr+u)*(-0.083219793_pr+u)*(1.2047515_pr &
                  +(-2.1485114_pr+u)*u)*(0.85574473_pr+(-1.4759269_pr+u)*u)+4.515625_pr*c1*(-1.0974449_pr+u) &
                  *(-0.076646725_pr+u)*(1.0942329_pr+(-1.9814273_pr+u)*u)*(0.7814348_pr+(-1.3562458_pr &
                  +u)*u))))/(fpi4*LambdaX*mpi2)

             dhJ1J1=mevfm*(ga*(-0.15272366_pr*cd*(-1.6176653_pr+u)*(-0.066866022_pr+u)*(1.7318817_pr &
                  +(-2.5016073_pr+u)*u)*(1.0139312_pr+(-1.425551_pr+u)*u)+ga*LambdaX*(-0.58984375_pr*c4 &
                  *(-0.69501214_pr+u)*(-0.088062144_pr+u)*(1.3949349_pr+(-2.3100216_pr+u)*u)*(0.99576969_pr &
                  +(-1.5747849_pr+u)*u)-4.0_pr*c3*(-0.59708253_pr+u)*(-0.031303824_pr+u)*(1.029716_pr &
                  +(-1.979428_pr+u)*u)*(0.81581869_pr+(-1.3671856_pr+u)*u)-16.0_pr*c1*(0.93008679_pr &
                  +(-1.8984052_pr+u)*u)*(0.65743897_pr+(-1.3677632_pr+u)*u)*(0.048712336_pr+u*(0.22866833_pr &
                  +u)))))/(fpi4*LambdaX*mpi2)

             ddhrho0rho0=mevfm*(3.5082063e-6_pr*ce*(-4.2415409_pr + u)*(-0.63600549_pr + u)*(0.36430595_pr &
                  - 1.2064733_pr*u + u2)*(0.36384534_pr - 1.2034803_pr*u + u2) + gA*(-0.10526469_pr*cd*(-15.717486_pr &
                  + u)*(-1.1063208_pr + u)*(-0.79219406_pr + u)*(-0.15821661_pr + u)*(0.99686003_pr - 1.7359229_pr*u &
                  + u2) + gA*LambdaX*(-23.140197_pr*c4*(-0.54937402_pr + u)*(-0.25617028_pr + u)*(0.87042847_pr &
                  - 1.8342813_pr*u + u2)*(0.55052287_pr - 1.1703549_pr*u + u2) + 34417.155_pr*c3*(-0.75220388_pr &
                  + u)*(-0.45847411_pr + u)*(0.49790159_pr - 1.3890715_pr*u + u2)*(0.31420079_pr - 1.0830145_pr*u &
                  + u2) + 150844.82_pr*c1*(0.5097462_pr - 1.4229463_pr*u + u2)*(0.38532363_pr - 1.2176864_pr*u &
                  + u2)*(0.27148791_pr - 1.0367952_pr*u + u2))))/(fpi4*LambdaX)

             ddhrho1rho1=mevfm*(3.5082063e-6_pr*ce*(-4.2415409_pr + u)*(-0.63600549_pr + u)*(0.36430595_pr &
                  - 1.2064733_pr*u + u2)*(0.36384534_pr - 1.2034803_pr*u + u2) + gA*(-0.10526469_pr*cd*(-15.717486_pr &
                  + u)*(-1.1063208_pr + u)*(-0.79219406_pr + u)*(-0.15821661_pr + u)*(0.99686003_pr - 1.7359229_pr*u &
                  + u2) + gA*LambdaX*(-23.140197_pr*c4*(-0.54937402_pr + u)*(-0.25617028_pr + u)*(0.87042847_pr &
                  - 1.8342813_pr*u + u2)*(0.55052287_pr - 1.1703549_pr*u + u2) + 34417.155_pr*c3*(-0.75220388_pr &
                  + u)*(-0.45847411_pr + u)*(0.49790159_pr - 1.3890715_pr*u + u2)*(0.31420079_pr - 1.0830145_pr*u &
                  + u2) + 150844.82_pr*c1*(0.5097462_pr - 1.4229463_pr*u + u2)*(0.38532363_pr - 1.2176864_pr*u &
                  + u2)*(0.27148791_pr - 1.0367952_pr*u + u2))))/(fpi4*LambdaX)

             ddhrho0Drho0=mevfm*(gA*(-3.1120954_pr*cd*(-0.30605132_pr + u)*(1.2680842_pr + u)*(1.0861033_pr &
                  - 2.0506502_pr*u + u2)*(0.79950278_pr - 1.5067227_pr*u + u2) + gA*LambdaX*(-24514.431_pr*c3*( &
                  -0.76031522_pr + u)*(-0.46096893_pr + u)*(0.49710047_pr - 1.3844361_pr*u + u2)*(0.29795986_pr &
                  - 1.055072_pr*u + u2) + 22229.743_pr*c1*(-0.74615955_pr + u)*(-0.61188783_pr + u)*(0.46154959_pr &
                  - 1.3315221_pr*u + u2)*(0.25673397_pr - 0.99266282_pr*u + u2) + 110.78144_pr*c4*(-0.8714267_pr &
                  + u)*(-0.31340402_pr + u)*(0.6979966_pr - 1.6133393_pr*u + u2)*(0.37679963_pr - 0.97187797_pr*u &
                  + u2))))/(fpi4*LambdaX*mpi2)

             ddhrho1Drho0=mevfm*(gA*(1.0374387_pr*cd*(-0.30605117_pr + u)*(1.2679976_pr + u)*(1.0860961_pr &
                  - 2.0506436_pr*u + u2)*(0.79949758_pr - 1.5067185_pr*u + u2) + gA*LambdaX*(7.4928249_pr*c1*( &
                  -0.37847001_pr + u)*(0.20959238_pr + u)*(1.3329648_pr - 2.2776684_pr*u + u2)*(0.96918292_pr &
                  - 1.6804594_pr*u + u2) + 0.50449085_pr*c4*(-0.25756953_pr + u)*(4.642167_pr + u)*(1.0022598_pr &
                  - 1.987984_pr*u + u2)*(0.89301348_pr - 1.6569412_pr*u + u2) - 8.2338865_pr*c3*(-0.24176677_pr &
                  + u)*(-0.021931235_pr + u)*(0.86956532_pr - 1.8470984_pr*u + u2)*(0.70044277_pr - 1.4791249_pr*u &
                  + u2))))/(fpi4*LambdaX*mpi2)

             ddhrho1Drho1=mevfm*(gA*(2.0747303_pr*cd*(-0.30605132_pr + u)*(1.2680842_pr + u)*(1.0861033_pr &
                  - 2.0506502_pr*u + u2)*(0.79950278_pr - 1.5067227_pr*u + u2) + gA*LambdaX*(15284.793_pr*c3*( &
                  -0.74604143_pr + u)*(-0.48901435_pr + u)*(0.48038865_pr - 1.3647474_pr*u + u2)*(0.30523896_pr &
                  - 1.0777475_pr*u + u2) - 73.854292_pr*c4*(-0.8714267_pr + u)*(-0.31340402_pr + u)*(0.6979966_pr &
                  - 1.6133393_pr*u + u2)*(0.37679963_pr - 0.97187797_pr*u + u2) - 14523.513_pr*c1*(0.51615934_pr &
                  - 1.4339319_pr*u + u2)*(0.42639291_pr - 1.2743751_pr*u + u2)*(0.23481515_pr - 0.95144731_pr*u &
                  + u2))))/(fpi4*LambdaX*mpi2)

             ddhrho0tau0=mevfm*(ga*(6.2246138_pr*cd*(-0.30605117_pr+u)*(1.2680013_pr+u)*(1.0860964_pr &
                  +(-2.0506438_pr+u)*u)*(0.79949779_pr+(-1.5067187_pr+u)*u)+ga*LambdaX*(7.25_pr*c4*(-0.25490416_pr &
                  +u)*(1.9314619_pr+u)*(0.97416501_pr+(-1.9598149_pr+u)*u)*(0.85402263_pr+(-1.6270877_pr+u)*u) &
                  -4096.0_pr*c1*(-0.82877841_pr+u)*(-0.61205527_pr+u)*(0.55294169_pr+(-1.4322212_pr+u)*u) &
                  *(0.21133499_pr+(-0.85194512_pr+u)*u)-1024.0_pr*c3*(-1.5110662_pr+u)*(-0.81113397_pr+u) &
                  *(0.48220899_pr+(-1.3083223_pr+u)*u)*(0.13559828_pr+(-0.71947752_pr+u)*u))))/(fpi4*LambdaX*mpi2)

             ddhrho1tau0=mevfm*(ga*(-2.0748713_pr*cd*(-0.30605117_pr+u)*(1.2680013_pr+u)*(1.0860964_pr &
                  +(-2.0506438_pr+u)*u)*(0.79949779_pr+(-1.5067187_pr+u)*u)+ga*LambdaX*(-2.1875_pr*c4 &
                  *(-0.25526539_pr+u)*(2.1267165_pr+u)*(0.98082879_pr+(-1.9667092_pr+u)*u)*(0.86027632_pr &
                  +(-1.6333133_pr+u)*u)+4.0_pr*c3*(-0.27578203_pr+u)*(0.73913647_pr+u)*(0.98670704_pr &
                  +(-1.9658775_pr+u)*u)*(0.85430872_pr+(-1.6131019_pr+u)*u)-29.0_pr*c1*(-0.37708308_pr+u) &
                  *(0.063860656_pr+u)*(1.1337991_pr+(-2.0962155_pr+u)*u)*(0.84554078_pr+(-1.5461655_pr+u)*u)))) &
                  /(fpi4*LambdaX*mpi2)

             ddhrho1tau1=mevfm*(ga*(-4.1497425_pr*cd*(-0.30605117_pr+u)*(1.2680013_pr+u)*(1.0860964_pr &
                  +(-2.0506438_pr+u)*u)*(0.79949779_pr+(-1.5067187_pr+u)*u)+ga*LambdaX*(-4.375_pr*c4 &
                  *(-0.25526539_pr+u)*(2.1267165_pr+u)*(0.98082879_pr+(-1.9667092_pr+u)*u)*(0.86027632_pr &
                  +(-1.6333133_pr+u)*u)-9728.0_pr*c1*(-0.72811509_pr+u)*(-0.40683586_pr+u)*(0.51101775_pr &
                  +(-1.4100026_pr+u)*u)*(0.32463579_pr+(-1.0747833_pr+u)*u)+1536.0_pr*c3*(0.61725539_pr &
                  +(-1.5571261_pr+u)*u)*(0.39877632_pr+(-1.1610973_pr+u)*u)*(0.15498418_pr+(-0.77760994_pr &
                  +u)*u))))/(fpi4*LambdaX*mpi2)

          End If !! if (u.gt.ucut3n..
       End If !! if(use_dme3N_terms
       !
    End If !! if (dmeorder.ge.2....
    !
    Urhorho=0.0_pr   ; Urhotau=0.0_pr
    UrhoDrho=0.0_pr  ; Unablarho=0.0_pr
    UJnablarho=0.0_pr; UrhonablaJ=0.0_pr
    Urhorhopr=0.0_pr ; UJJ=0.0_pr
    UJabJba=0.0_pr
    UEnonstdr=0.0_pr ; UFnonstdr=0.0_pr ; URnonstdr=0.0_pr
    !
    ! Notations for Uamplitudes(0:3,0:7)
    ! t for Uamplitudes(t,*)
    ! 0 -> 0,0
    ! 1 -> 1,1
    ! 2 -> 0,1
    ! 3 -> 1,0
    ! n for Uamplitudes(*,n)
    ! 0 -> U
    ! 1 -> dU/dRHO_0
    ! 2 -> dU/dRHO_1
    ! 3 -> d2U/(dRHO_0*dRHO_0)
    ! 4 -> d2U/(dRHO_1*dRHO_1)
    ! 5 -> d2U/(dRHO_0*dRHO_1)
    ! 6 -> dU/d(TAU_0)
    ! 7 -> dU/d(Delta RHO_0)
    !
    !! 2N terms
    Do t=0,1
       ph=1.0_pr
       If(t.Eq.1) ph=-1.0_pr
       Urhorho(t,0)=Crho(t)+Cdrho(t)*rho(0)**sigma &
            +0.50_pr*(arhorho+ph*brhorho)*mevfm
       Urhotau(t,0)=Ctau(t)+0.50_pr*(arhotau+ph*brhotau)*mevfm        !! These two determine the
       UrhoDrho(t,0)=Crdr(t)+ac2*0.50_pr*(arhoDrho+ph*brhoDrho)*mevfm !! effective mass (when recoupling to p-n)??
       UJJ(t,0)=CJ(t)+0.50_pr*(ajj+ph*bjj)*mevfm
       Unablarho(t,0)=Cnrho(t)+0.50_pr*(adrdr+ph*bdrdr)*mevfm
       UrhonablaJ(t,0)=Crdj(t)
       UJnablarho(t,0)=Cjdr(t)

       Urhorho(t,1)=sigma*Cdrho(t)*(rho(0)**sigma)/(rho(0)+eps) &
            +0.50_pr*(darhorho+ph*dbrhorho)*du*mevfm
       Urhotau(t,1)=0.50_pr*(darhotau+ph*dbrhotau)*du*mevfm
       UrhoDrho(t,1)=ac2*0.50_pr*(darhoDrho+ph*dbrhoDrho)*du*mevfm
       UJJ(t,1)=0.50_pr*(dajj+ph*dbjj)*du*mevfm
       Unablarho(t,1)=0.50_pr*(dadrdr+ph*dbdrdr)*du*mevfm

       Urhorho(t,6)=0.50_pr*(darhorho+ph*dbrhorho)*dtu*mevfm
       Urhotau(t,6)=0.50_pr*(darhotau+ph*dbrhotau)*dtu*mevfm
       UrhoDrho(t,6)=ac2*0.50_pr*(darhoDrho+ph*dbrhoDrho)*dtu*mevfm
       UJJ(t,6)=0.50_pr*(dajj+ph*dbjj)*dtu*mevfm
       Unablarho(t,6)=0.50_pr*(dadrdr+ph*dbdrdr)*dtu*mevfm

       Urhorho(t,7)=0.50_pr*(darhorho+ph*dbrhorho)*dlu*mevfm
       Urhotau(t,7)=0.50_pr*(darhotau+ph*dbrhotau)*dlu*mevfm
       UrhoDrho(t,7)=ac2*0.50_pr*(darhoDrho+ph*dbrhoDrho)*dlu*mevfm
       UJJ(t,7)=0.50_pr*(dajj+ph*dbjj)*dlu*mevfm
       Unablarho(t,7)=0.50_pr*(dadrdr+ph*dbdrdr)*dlu*mevfm

       Urhorho(t,3)=sigma*(sigma-1.0_pr)*Cdrho(t)*(rho(0)**sigma)/(rho(0)**2+eps) &
            +0.50_pr*(darhorho+ph*dbrhorho)*ddu*mevfm &
            +0.50_pr*(ddarhorho+ph*ddbrhorho)*du*du*mevfm
       Urhotau(t,3)=0.50_pr*(darhotau+ph*dbrhotau)*ddu*mevfm &
            +0.50_pr*(ddarhotau+ph*ddbrhotau)*du*du*mevfm
       UrhoDrho(t,3)=ac2*0.50_pr*(darhoDrho+ph*dbrhoDrho)*ddu*mevfm &
            +ac2*0.50_pr*(ddarhoDrho+ph*ddbrhoDrho)*du*du*mevfm

    End Do
    Urhorhopr(0,0) = CpV0(0)*(1.0_pr-CpV1(0)*rho(0)/0.16_pr)          &
         +CpV0(1)*(1.0_pr-CpV1(1)*rho(0)/0.16_pr)
    Urhorhopr(1,0) = CpV0(0)*(1.0_pr-CpV1(0)*rho(0)/0.16_pr)          &
         +CpV0(1)*(1.0_pr-CpV1(1)*rho(0)/0.16_pr)
    Urhorhopr(2,0) = (CpV0(0)*(1.0_pr-CpV1(0)*rho(0)/0.16_pr)         &
         -CpV0(1)*(1.0_pr-CpV1(1)*rho(0)/0.16_pr))*2.0_pr
    Urhorhopr(0,1) = (-CpV0(0)*CpV1(0)-CpV0(1)*CpV1(1))/0.16_pr
    Urhorhopr(1,1) = (-CpV0(0)*CpV1(0)-CpV0(1)*CpV1(1))/0.16_pr
    Urhorhopr(2,1) = 2.0_pr*(-CpV0(0)*CpV1(0)+CpV0(1)*CpV1(1))/0.16_pr
    Urhorhopr=Urhorhopr/16.0_pr
    !! 3N terms
    If (use_DME3N_terms) Then

       Urhorho(0,0)=Urhorho(0,0)+hrho0rho0*rho(0)
       Urhorho(1,0)=Urhorho(1,0)+hrho1rho1*rho(0)
       Unablarho(0,0)=Unablarho(0,0)+hdr0dr0*rho(0)
       Unablarho(1,0)=Unablarho(1,0)+hdr1dr1*rho(0)
       UrhoDrho(0,0)=UrhoDrho(0,0)+hrho0Drho0*rho(0)*ac3
       UrhoDrho(1,0)=UrhoDrho(1,0)+hrho1Drho1*rho(0)*ac3
       UrhoDrho(3,0)=UrhoDrho(3,0)+hrho1Drho0*rho(1)*ac3
       Urhotau(0,0)=Urhotau(0,0)+hrho0tau0*rho(0)
       Urhotau(1,0)=Urhotau(1,0)+hrho1tau1*rho(0)
       Urhotau(3,0)=Urhotau(3,0)+hrho1tau0*rho(1)
       UJnablarho(0,0)=UJnablarho(0,0)+hJ0dr0*rho(0)
       UJnablarho(1,0)=UJnablarho(1,0)+hJ1dr1*rho(0)
       UJnablarho(2,0)=UJnablarho(2,0)+hJ0dr1*rho(1)
       UJnablarho(3,0)=UJnablarho(3,0)+hJ1dr0*rho(1)
       UrhonablaJ(0,0)=UrhonablaJ(0,0)+hrho0dJ0*rho(0)
       UrhonablaJ(1,0)=UrhonablaJ(1,0)+hrho1dJ1*rho(0)
       UrhonablaJ(3,0)=UrhonablaJ(3,0)+hrho1dJ0*rho(1)
       UJJ(0,0)=UJJ(0,0)+hJ0J0*rho(0)
       UJJ(1,0)=UJJ(1,0)+hJ1J1*rho(0)
       UJJ(2,0)=UJJ(2,0)+hJ0J1*rho(1)

       Urhorho(0,1)=Urhorho(0,1)+dhrho0rho0*du*rho(0) +hrho0rho0
       Urhorho(1,1)=Urhorho(1,1)+dhrho1rho1*du*rho(0) +hrho1rho1
       Unablarho(0,1)=Unablarho(0,1)+dhdr0dr0*du*rho(0) +hdr0dr0
       Unablarho(1,1)=Unablarho(1,1)+dhdr1dr1*du*rho(0) +hdr1dr1
       UrhoDrho(0,1)=UrhoDrho(0,1)+(dhrho0Drho0*du*rho(0) +hrho0Drho0)*ac3
       UrhoDrho(1,1)=UrhoDrho(1,1)+(dhrho1Drho1*du*rho(0) +hrho1Drho1)*ac3
       UrhoDrho(3,1)=UrhoDrho(3,1)+dhrho1Drho0*du*rho(1)*ac3
       Urhotau(0,1)=Urhotau(0,1)+dhrho0tau0*du*rho(0) +hrho0tau0
       Urhotau(1,1)=Urhotau(1,1)+dhrho1tau1*du*rho(0) +hrho1tau1
       Urhotau(3,1)=Urhotau(3,1)+dhrho1tau0*du*rho(1)
       UJnablarho(0,1)=UJnablarho(0,1)+dhJ0dr0*du*rho(0) +hJ0dr0
       UJnablarho(1,1)=UJnablarho(1,1)+dhJ1dr1*du*rho(0) +hJ1dr1
       UJnablarho(2,1)=UJnablarho(2,1)+dhJ0dr1*du*rho(1)
       UJnablarho(3,1)=UJnablarho(3,1)+dhJ1dr0*du*rho(1)
       UrhonablaJ(0,1)=UrhonablaJ(0,1)+dhrho0dJ0*du*rho(0) +hrho0dJ0
       UrhonablaJ(1,1)=UrhonablaJ(1,1)+dhrho1dJ1*du*rho(0) +hrho1dJ1
       UrhonablaJ(3,1)=UrhonablaJ(3,1)+dhrho1dJ0*du*rho(1)
       UJJ(0,1)=UJJ(0,1)+dhJ0J0*du*rho(0) +hJ0J0
       UJJ(1,1)=UJJ(1,1)+dhJ1J1*du*rho(0) +hJ1J1
       UJJ(2,1)=UJJ(2,1)+dhJ0J1*du*rho(1)

       Urhorho(0,6)=Urhorho(0,6)+dhrho0rho0*dtu*rho(0)
       Urhorho(1,6)=Urhorho(1,6)+dhrho1rho1*dtu*rho(0)
       Unablarho(0,6)=Unablarho(0,6)+dhdr0dr0*dtu*rho(0)
       Unablarho(1,6)=Unablarho(1,6)+dhdr1dr1*dtu*rho(0)
       UrhoDrho(0,6)=UrhoDrho(0,6)+dhrho0Drho0*dtu*rho(0)*ac3
       UrhoDrho(1,6)=UrhoDrho(1,6)+dhrho1Drho1*dtu*rho(0)*ac3
       UrhoDrho(3,6)=UrhoDrho(3,6)+dhrho1Drho0*dtu*rho(1)*ac3
       Urhotau(0,6)=Urhotau(0,6)+dhrho0tau0*dtu*rho(0)
       Urhotau(1,6)=Urhotau(1,6)+dhrho1tau1*dtu*rho(0)
       Urhotau(3,6)=Urhotau(3,6)+dhrho1tau0*dtu*rho(1)
       UJnablarho(0,6)=UJnablarho(0,6)+dhJ0dr0*dtu*rho(0)
       UJnablarho(1,6)=UJnablarho(1,6)+dhJ1dr1*dtu*rho(0)
       UJnablarho(2,6)=UJnablarho(2,6)+dhJ0dr1*dtu*rho(1)
       UJnablarho(3,6)=UJnablarho(3,6)+dhJ1dr0*dtu*rho(1)
       UrhonablaJ(0,6)=UrhonablaJ(0,6)+dhrho0dJ0*dtu*rho(0)
       UrhonablaJ(1,6)=UrhonablaJ(1,6)+dhrho1dJ1*dtu*rho(0)
       UrhonablaJ(3,6)=UrhonablaJ(3,6)+dhrho1dJ0*dtu*rho(1)
       UJJ(0,6)=UJJ(0,6)+dhJ0J0*dtu*rho(0)
       UJJ(1,6)=UJJ(1,6)+dhJ1J1*dtu*rho(0)
       UJJ(2,6)=UJJ(2,6)+dhJ0J1*dtu*rho(1)

       Urhorho(0,7)=Urhorho(0,7)+dhrho0rho0*dlu*rho(0)
       Urhorho(1,7)=Urhorho(1,7)+dhrho1rho1*dlu*rho(0)
       Unablarho(0,7)=Unablarho(0,7)+dhdr0dr0*dlu*rho(0)
       Unablarho(1,7)=Unablarho(1,7)+dhdr1dr1*dlu*rho(0)
       UrhoDrho(0,7)=UrhoDrho(0,7)+dhrho0Drho0*dlu*rho(0)*ac3
       UrhoDrho(1,7)=UrhoDrho(1,7)+dhrho1Drho1*dlu*rho(0)*ac3
       UrhoDrho(3,7)=UrhoDrho(3,7)+dhrho1Drho0*dlu*rho(1)*ac3
       Urhotau(0,7)=Urhotau(0,7)+dhrho0tau0*dlu*rho(0)
       Urhotau(1,7)=Urhotau(1,7)+dhrho1tau1*dlu*rho(0)
       Urhotau(3,7)=Urhotau(3,7)+dhrho1tau0*dlu*rho(1)
       UJnablarho(0,7)=UJnablarho(0,7)+dhJ0dr0*dlu*rho(0)
       UJnablarho(1,7)=UJnablarho(1,7)+dhJ1dr1*dlu*rho(0)
       UJnablarho(2,7)=UJnablarho(2,7)+dhJ0dr1*dlu*rho(1)
       UJnablarho(3,7)=UJnablarho(3,7)+dhJ1dr0*dlu*rho(1)
       UrhonablaJ(0,7)=UrhonablaJ(0,7)+dhrho0dJ0*dlu*rho(0)
       UrhonablaJ(1,7)=UrhonablaJ(1,7)+dhrho1dJ1*dlu*rho(0)
       UrhonablaJ(3,7)=UrhonablaJ(3,7)+dhrho1dJ0*dlu*rho(1)
       UJJ(0,7)=UJJ(0,7)+dhJ0J0*dlu*rho(0)
       UJJ(1,7)=UJJ(1,7)+dhJ1J1*dlu*rho(0)
       UJJ(2,7)=UJJ(2,7)+dhJ0J1*dlu*rho(1)

       Urhorho(0,3)=Urhorho(0,3)+2.0_pr*dhrho0rho0*du &
                   +ddhrho0rho0*du*du*rho(0) +dhrho0rho0*ddu*rho(0)
       Urhorho(1,3)=Urhorho(1,3)+2.0_pr*dhrho1rho1*du &
                   +ddhrho1rho1*du*du*rho(0) +dhrho1rho1*ddu*rho(0)
       UrhoDrho(0,3)=UrhoDrho(0,3)+(2.0_pr*dhrho0Drho0*du &
                    +ddhrho0Drho0*du*du*rho(0) +dhrho0Drho0*ddu*rho(0))*ac3
       UrhoDrho(1,3)=UrhoDrho(1,3)+(2.0_pr*dhrho1Drho1*du &
                    +ddhrho1Drho1*du*du*rho(0)+dhrho1Drho1*ddu*rho(0))*ac3
       Urhotau(0,3)=Urhotau(0,3)+2.0_pr*dhrho0tau0*du &
                   +ddhrho0tau0*du*du*rho(0) +dhrho0tau0*ddu*rho(0)
       Urhotau(1,3)=Urhotau(1,3)+2.0_pr*dhrho1tau1*du &
                   +ddhrho1tau1*du*du*rho(0) +dhrho1tau1*ddu*rho(0)
       UrhoDrho(3,3)=UrhoDrho(3,3) +(ddhrho1Drho0*du*du*rho(1) &
                    +dhrho1Drho0*ddu*rho(1))*ac3
       Urhotau(3,3)=Urhotau(3,3) +ddhrho1tau0*du*du*rho(1) &
                   +dhrho1tau0*ddu*rho(1)

       UrhoDrho(3,2)=UrhoDrho(3,2)+hrho1Drho0*ac3
       Urhotau(3,2)=Urhotau(3,2)+hrho1tau0
       UJnablarho(3,2)=UJnablarho(3,2)+hJ1dr0
       UJnablarho(2,2)=UJnablarho(2,2)+hJ0dr1
       UrhonablaJ(3,2)=UrhonablaJ(3,2)+hrho1dJ0
       UJJ(2,2)=UJJ(2,2)+hJ0J1

       UrhoDrho(3,5)=UrhoDrho(3,5) +dhrho1Drho0*du*ac3
       Urhotau(3,5)=Urhotau(3,5) +dhrho1tau0*du

    End If   !! if(use3Nterm..

    UEnonstdr=0.0_pr; UFnonstdr=0.0_pr; URnonstdr=0.0_pr
    if(force_is_dme) then
       call Add_gcoupling_to_U(grhorho,rho(0),Urhorho)
       call Add_gcoupling_to_U(grhotau,rho(0),Urhotau)
       call Add_gcoupling_to_U(grhodelrho,rho(0),UrhoDrho)
       call Add_gcoupling_to_U(gJabJab,rho(0),UJJ)
       call Add_gcoupling_to_U(gJabJba,rho(0),UJabJba)
       if(use_3N_couplings) then
          call Add_hcoupling_to_U(h_rho_rho_rho,rho,Urhorho)
          call Add_hcoupling_to_U(h_rho_rho_tau,rho,Urhotau)
          call Add_hcoupling_to_U(h_rho_rho_delrho,rho,UrhoDrho)
          call Add_hcoupling_to_U(h_rho_nablarho,rho,Unablarho)
          call Add_hcoupling_to_U(h_rho_J_nablarho,rho,UJnablarho)
          call Add_hcoupling_to_U(h_rho_Jab_Jab,rho,UJJ)
          call Add_hcoupling_to_U(h_rho_Jab_Jba,rho,UJabJba)
       endif
       call Calculate_gh_for_INM(rho(0))
    endif
    !
    if(is_nedf) then
       do  t = 0,1
          Urhorho(t,0) =   a_NEDF(t)/(rho(0)**(1/three)+eps) + b_NEDF(t) + c_NEDF(t)*rho(0)**( 1/three)
          Urhorho(t,1) = (-a_NEDF(t)/(rho(0)**(4/three)+eps)             + c_NEDF(t)/(rho(0)**(2/three)+eps))/three
       enddo
       Urhorho(2,0) =  ( a_NEDF(2)/(rho(0)**(10/three)+eps) + &
                         b_NEDF(2)/(rho(0)**(3)+eps) + &
                         c_NEDF(2)/(rho(0)**( 8/three)+eps))*rho(1)**3
       Urhorho(2,1) =  (-10*a_NEDF(2)/(rho(0)**(13/three)+eps) - &
                          9*b_NEDF(2)/(rho(0)**(4)+eps) - &
                          8*c_NEDF(2)/(rho(0)**(11/three)+eps))*rho(1)**3/three
       Urhorho(2,2) =3*( a_NEDF(2)/(rho(0)**(10/three)+eps) + &
                         b_NEDF(2)/(rho(0)**(3)+eps) + &
                         c_NEDF(2)/(rho(0)**( 8/three)+eps))*rho(1)**2
    endif
    !
    If (.Not.use_j2terms) Then
       UJJ= zero
       UJabJba = zero
    End If
    !
  End Subroutine calculate_U_parameters

  !=====================================================================
  !> Adds he density dependent 2N couplings \f$ g(\rho_0) \f$ to the
  !> contact couplings stored in the U_parameter input array.
  !>
  !> The density dependent coupling is parametrized as
  !> \f[
  !>    g(\rho_0) = g_0 + \sum_{i=1}^{n_{\rm dme}} a_i
  !>    \left[\tan^{-1}(b_i \rho_0^{c_i})\right]^i.
  !> \f]
  !> The parameters \f$g_0\f$, \f$a_i\f$, \f$b_i\f$, \f$c_i\f$ are
  !> stored in the input array g_coupling.
  !>
  !> The first and second derivates of \f$g(\rho_0)\f$ with respect
  !> of \f$\rho_0\f$ are also stored in U_parameter
  !=====================================================================
  subroutine Add_gcoupling_to_U(g_coupling,rho0,U_parameter)
    implicit none
    real(pr), dimension(0:1,1:ndme), intent(in) :: g_coupling !< Array with the parameters defining the 2N coupling
    real(pr), intent(in) :: rho0 !< Isoscalar density at which the coupling is evaluated
    real(pr), dimension(0:3,0:7), intent(inout) :: U_parameter !< Array containing the functional couplings and its derivatives
    integer :: it,idme,n
    real(pr) :: f1,f2,f3
    do it = 0,1
       U_parameter(it,0) = U_parameter(it,0) + g_coupling(it,1)*alfa_dme
       do idme = 2,ndme,3
          f1 = g_coupling(it,idme)
          f2 = g_coupling(it,idme+1)
          f3 = g_coupling(it,idme+2)
          n = (idme+1)/3
          U_parameter(it,0) = U_parameter(it,0) + f1*atan(f2*rho0**f3)**n*alfa_dme
          U_parameter(it,1) = U_parameter(it,1) + &
               (f1*f2*f3*n*(rho0)**(f3-1)*&
               atan(f2*rho0**f3)**(n-1))/(1 + f2**2*rho0**(2*f3))*alfa_dme
          U_parameter(it,3) = U_parameter(it,3) + &
               (f1*f2*f3*n*rho0**(f3-2)*atan(f2*rho0**f3)**(n-2)*&
               (f2*f3*(n-1)*rho0**f3 + (f3-1-f2**2*(f3+1)*rho0**(2*f3))*atan(f2*rho0**f3)))/&
               (1 + f2**2*rho0**(2*f3))**2*alfa_dme
       enddo
    enddo
  end subroutine Add_gcoupling_to_U

  !=====================================================================
  !> Adds he density dependent 3N couplings \f$ \rho_{j} h(\rho_0)\f$,
  !> with \f$j=0,1\f$, to the couplings stored in the U_parameter input
  !> array.
  !>
  !> The density dependent coupling is parametrized as
  !> \f[
  !>    h(\rho_0) = h_0 + \sum_{i=1}^{n_{\rm dme}} a_i
  !>    \left[\tan^{-1}(b_i \rho_0^{c_i})\right]^i.
  !> \f]
  !> The parameters \f$h_0\f$, \f$a_i\f$, \f$b_i\f$, \f$c_i\f$ are
  !> stored in the input array h_coupling.
  !>
  !> The first and second derivates of \f$\rho_{j}h(\rho_0)\f$ with
  !> respect of \f$rho_0\f$ and \f$\rho_1\f$ are also stored in
  !> U_parameter
  !=====================================================================
  subroutine Add_hcoupling_to_U(h_coupling,rho_in,U_parameter)
    implicit none
    real(pr), dimension(0:3,1:ndme), intent(in) :: h_coupling !<  Array with the parameters defining the 2N coupling
    real(pr), dimension(0:1), intent(in) :: rho_in !< Array with the isoscalar and isovector densities
    real(pr), dimension(0:3,0:7), intent(inout) :: U_parameter !< Array containing the functional couplings and its derivatives
    real(pr), dimension(0:3,0:7) :: h_parameter
    real(pr) :: rho0, rho1
    integer :: it,idme,n
    real(pr) :: f1,f2,f3
    rho0 = rho_in(0)
    rho1 = rho_in(1)
    h_parameter = 0._pr
    do it = 0,3
       h_parameter(it,0) = h_coupling(it,1)
       do idme = 2,ndme,3
          f1 = h_coupling(it,idme)
          f2 = h_coupling(it,idme+1)
          f3 = h_coupling(it,idme+2)
          n = (idme+1)/3
          h_parameter(it,0) = h_parameter(it,0) + f1*atan(f2*rho0**f3)**n
          h_parameter(it,1) = h_parameter(it,1) + &
               (f1*f2*f3*n*(rho0)**(f3-1)*&
               atan(f2*rho0**f3)**(n-1))/(1 + f2**2*rho0**(2*f3))
          h_parameter(it,3) = h_parameter(it,3) + &
               (f1*f2*f3*n*rho0**(f3-2)*atan(f2*rho0**f3)**(n-2)*&
               (f2*f3*(n-1)*rho0**f3 + (f3-1-f2**2*(f3+1)*rho0**(2*f3))*atan(f2*rho0**f3)))/&
               (1 + f2**2*rho0**(2*f3))**2
       enddo
    enddo
    do it = 0,1
       h_parameter(it,3) = 2*h_parameter(it,1) + rho0*h_parameter(it,3)
       h_parameter(it,1) = h_parameter(it,0) + rho0*h_parameter(it,1)
       h_parameter(it,0) = rho0*h_parameter(it,0)
    enddo
    do it = 2,3
       h_parameter(it,5) = h_parameter(it,1)
       h_parameter(it,2) = h_parameter(it,0)
       h_parameter(it,3) = rho1*h_parameter(it,3)
       h_parameter(it,1) = rho1*h_parameter(it,1)
       h_parameter(it,0) = rho1*h_parameter(it,0)
    enddo
    U_parameter = U_parameter + h_parameter*alfa_dme
  end subroutine Add_hcoupling_to_U

  !=====================================================================
  !> Calculates the density dependent couplings \f$g(\rho_0)\f$ and
  !> \f$h(\rho_0)\f$ along with the first and second derivatives that
  !> are used in the calculate_C_from_NM() and calculate_NM_properties()
  !> subroutines.
  !>
  !> While the couplings are a function of \f$\rho_0\f$, the derivates
  !> are calculated with respect of \f$ u = (k^F/m_\pi)\rho_0^{1/3}\f$
  !> since those are the derivatives explicitely used in the mentioned
  !> subroutines.
  !>
  !> The specific couplings calculated here are \f$g^{\rho_0\rho_0}\f$,
  !> \f$g^{\rho_1\rho_1}\f$, \f$g^{\rho_0\tau_0}\f$,
  !> \f$g^{\rho_1\tau_1}\f$, \f$h^{\rho_0\rho_0}\f$,
  !> \f$h^{\rho_1\rho_1}\f$, \f$h^{\rho_0\tau_0}\f$,
  !> \f$h^{\rho_1\tau_1}\f$, \f$h^{\rho_1\tau_0}\f$.
  !=====================================================================
  Subroutine Calculate_gh_for_INM(rho0)
    implicit none
    real(pr), intent(in) :: rho0 !< Isoscalar density at which the couplings are evaluated
    real(pr) :: g_rho0rho0, dg_rho0rho0, ddg_rho0rho0
    real(pr) :: g_rho1rho1, dg_rho1rho1, ddg_rho1rho1
    real(pr) :: g_rho0tau0, dg_rho0tau0, ddg_rho0tau0
    real(pr) :: g_rho1tau1, dg_rho1tau1, ddg_rho1tau1
    real(pr) :: h_rho0rho0, dh_rho0rho0, ddh_rho0rho0
    real(pr) :: h_rho1rho1, dh_rho1rho1, ddh_rho1rho1
    real(pr) :: h_rho0tau0, dh_rho0tau0, ddh_rho0tau0
    real(pr) :: h_rho1tau1, dh_rho1tau1, ddh_rho1tau1
    real(pr) :: h_rho1tau0, dh_rho1tau0, ddh_rho1tau0
    real(pr) :: aRho0Rho0,daRho0Rho0,ddaRho0Rho0,aRho1Rho1,daRho1Rho1,ddaRho1Rho1
    real(pr) :: aRho0Tau0,daRho0Tau0,ddaRho0Tau0,aRho1Tau1,daRho1Tau1,ddaRho1Tau1
    real(pr) :: u, drho_du, ddrho_duu, CkF
    call calculate_ghcoupling(grhorho,1,0,rho0,g_rho0rho0,dg_rho0rho0,ddg_rho0rho0)
    call calculate_ghcoupling(grhorho,1,1,rho0,g_rho1rho1,dg_rho1rho1,ddg_rho1rho1)
    call calculate_ghcoupling(grhotau,1,0,rho0,g_rho0tau0,dg_rho0tau0,ddg_rho0tau0)
    call calculate_ghcoupling(grhotau,1,1,rho0,g_rho1tau1,dg_rho1tau1,ddg_rho1tau1)
    if(use_3N_couplings) then
       call calculate_ghcoupling(h_rho_rho_rho,3,0,rho0,h_rho0rho0,dh_rho0rho0,ddh_rho0rho0)
       call calculate_ghcoupling(h_rho_rho_rho,3,1,rho0,h_rho1rho1,dh_rho1rho1,ddh_rho1rho1)
       call calculate_ghcoupling(h_rho_rho_tau,3,0,rho0,h_rho0tau0,dh_rho0tau0,ddh_rho0tau0)
       call calculate_ghcoupling(h_rho_rho_tau,3,1,rho0,h_rho1tau1,dh_rho1tau1,ddh_rho1tau1)
       call calculate_ghcoupling(h_rho_rho_tau,3,3,rho0,h_rho1tau0,dh_rho1tau0,ddh_rho1tau0)
    else
       h_rho0rho0 = zero; dh_rho0rho0 = zero; ddh_rho0rho0 = zero
       h_rho1rho1 = zero; dh_rho1rho1 = zero; ddh_rho1rho1 = zero
       h_rho0tau0 = zero; dh_rho0tau0 = zero; ddh_rho0tau0 = zero
       h_rho1tau1 = zero; dh_rho1tau1 = zero; ddh_rho1tau1 = zero
       h_rho1tau0 = zero; dh_rho1tau0 = zero; ddh_rho1tau0 = zero
    endif
    CkF = kfconst/mpi
    u = CkF*rho0**(1._pr/3._pr) ! => rho0(u) = (u/CkF)**3
    drho_du = 3*u**2/(CkF**3)   ! rho0'(u)
    ddrho_duu = 6*u/(CkF**3)    ! rho0''(u)
    aRho0Rho0 = g_rho0rho0
    aRho1Rho1 = g_rho1rho1
    aRho0Tau0 = g_rho0tau0
    aRho1Tau1 = g_rho1tau1
    hRho0Rho0 = h_rho0rho0
    hRho1Rho1 = h_rho1rho1
    hRho0Tau0 = h_rho0tau0
    hRho1Tau1 = h_rho1tau1
    hRho1Tau0 = h_rho1tau0
    daRho0Rho0 = dg_rho0rho0*drho_du
    daRho1Rho1 = dg_rho1rho1*drho_du
    daRho0Tau0 = dg_rho0tau0*drho_du
    daRho1Tau1 = dg_rho1tau1*drho_du
    dhRho0Rho0 = dh_rho0rho0*drho_du
    dhRho1Rho1 = dh_rho1rho1*drho_du
    dhRho0Tau0 = dh_rho0tau0*drho_du
    dhRho1Tau1 = dh_rho1tau1*drho_du
    dhRho1Tau0 = dh_rho1tau0*drho_du
    ddaRho0Rho0 = ddg_rho0rho0*drho_du**2 + dg_rho0rho0*ddrho_duu
    ddaRho1Rho1 = ddg_rho1rho1*drho_du**2 + dg_rho1rho1*ddrho_duu
    ddaRho0Tau0 = ddg_rho0tau0*drho_du**2 + dg_rho0tau0*ddrho_duu
    ddaRho1Tau1 = ddg_rho1tau1*drho_du**2 + dg_rho1tau1*ddrho_duu
    ddhRho0Rho0 = ddh_rho0rho0*drho_du**2 + dh_rho0rho0*ddrho_duu
    ddhRho1Rho1 = ddh_rho1rho1*drho_du**2 + dh_rho1rho1*ddrho_duu
    ddhRho0Tau0 = ddh_rho0tau0*drho_du**2 + dh_rho0tau0*ddrho_duu
    ddhRho1Tau1 = ddh_rho1tau1*drho_du**2 + dh_rho1tau1*ddrho_duu
    ddhRho1Tau0 = ddh_rho1tau0*drho_du**2 + dh_rho1tau0*ddrho_duu

    aRhoRho = (aRho0Rho0+aRho1Rho1)/mevfm
    bRhoRho = (aRho0Rho0-aRho1Rho1)/mevfm
    aRhoTau = (aRho0Tau0+aRho1Tau1)/mevfm
    bRhoTau = (aRho0Tau0-aRho1Tau1)/mevfm
    daRhoRho = (daRho0Rho0+daRho1Rho1)/mevfm
    dbRhoRho = (daRho0Rho0-daRho1Rho1)/mevfm
    daRhoTau = (daRho0Tau0+daRho1Tau1)/mevfm
    dbRhoTau = (daRho0Tau0-daRho1Tau1)/mevfm
    ddaRhoRho = (ddaRho0Rho0+ddaRho1Rho1)/mevfm
    ddbRhoRho = (ddaRho0Rho0-ddaRho1Rho1)/mevfm
    ddaRhoTau = (ddaRho0Tau0+ddaRho1Tau1)/mevfm
    ddbRhoTau = (ddaRho0Tau0-ddaRho1Tau1)/mevfm
  End subroutine Calculate_gh_for_INM

  !=====================================================================
  !> Given the set of parameters that define the 2N or 3N density
  !> dependent coupling, calculates the corresponding coupling at a
  !> given density \f$\rho_0\f$ along with the first and second
  !> derivatives with respect of \f$\rho_0\f$.
  !>
  !> The size of the gh_par array is given by np and depends on the
  !> coupling comming from a 2N (np = 1) or 3N (np = 3) force.
  !> The index it indicates which coupling is calculated. For example:
  !>  - it = 0, \f$h^{\rho_0\tau_0}\f$
  !>  - it = 1, \f$h^{\rho_1\tau_1}\f$
  !>  - it = 2, \f$h^{\rho_0\tau_1}\f$
  !>  - it = 3, \f$h^{\rho_1\tau_0}\f$
  !=====================================================================
  Subroutine calculate_ghcoupling(gh_par,np,it,rho0,gh,dgh,ddgh)
    implicit none
    real(pr), dimension(0:np,1:ndme), intent(in) :: gh_par !< Array with the parameters defining the 2N or 3N coupling.
    integer, intent(in) :: np !< Size of the gh_par array.
    integer, intent(in) :: it !< Index indicating which coupling to calculate.
    real(pr), intent(in) :: rho0 !< Isoscalar density at which the coupling is evaluated.
    real(pr), intent(out) :: gh !< 2N or 3N coupling.
    real(pr), intent(out) :: dgh !< Derivative of gh with respect of rho0.
    real(pr), intent(out) :: ddgh !< Second derivative of gh with respect of rho0.
    integer :: idme,n
    real(pr) :: f1,f2,f3
    gh = zero; dgh = zero; ddgh = zero
    gh = gh_par(it,1)
    do idme = 2,ndme,3
       f1 = gh_par(it,idme)
       f2 = gh_par(it,idme+1)
       f3 = gh_par(it,idme+2)
       n = (idme+1)/3
       gh = gh + f1*atan(f2*rho0**f3)**n
       dgh = dgh + (f1*f2*f3*n*(rho0)**(f3-1)*&
            atan(f2*rho0**f3)**(n-1))/(1 + f2**2*rho0**(2*f3))
       ddgh = ddgh + (f1*f2*f3*n*rho0**(f3-2)*atan(f2*rho0**f3)**(n-2)*&
            (f2*f3*(n-1)*rho0**f3 + (f3-1-f2**2*(f3+1)*rho0**(2*f3))*atan(f2*rho0**f3)))/&
            (1 + f2**2*rho0**(2*f3))**2
    enddo
    gh = gh*alfa_dme;
    dgh = dgh*alfa_dme;
    ddgh = ddgh*alfa_dme;
  End Subroutine calculate_ghcoupling
  !=======================================================================
  !
  !=======================================================================
  Subroutine load_tables
    Implicit None
    !
    !  Urho0rho0
    !
    ctr0r0=0.0_pr
    ctr0r0(1,1,1) = (-0.964246762695464_pr*c3*gA2)/fpi4 - &
         (0.0263462308309038_pr*c4*gA2)/fpi4 - (0.1875_pr*ce)/(fpi4*LambdaX) + &
         (0.046875_pr*cd*ga)/(fpi4*LambdaX)
    ctr0r0(1,1,3) = (-6.84195361981036_pr*c1*gA2)/fpi4 + &
         (5.01410456224811_pr*c3*gA2)/fpi4 + (0.0732776804846939_pr*c4*gA2)/fpi4 &
         - (0.06328125_pr*cd*ga)/(fpi4*LambdaX)
    ctr0r0(2,1,3) = (0.147698198400372_pr*c3*gA2)/fpi4 + &
         (0.0703526785714286_pr*c4*gA2)/fpi4 - &
         (0.02109375_pr*cd*ga)/(fpi4*LambdaX)
    ctr0r0(1,2,4) = (-0.1875_pr*c3*gA2)/fpi4 - (0.375_pr*c4*gA2)/fpi4 + &
         (0.140625_pr*cd*ga)/(fpi4*LambdaX)
    ctr0r0(1,1,5) = (-23.1434396482781_pr*c1*gA2)/fpi4 + &
         (12.2465870768616_pr*c3*gA2)/fpi4 - (0.18795405731824_pr*c4*gA2)/fpi4 + &
         (0.017578125_pr*cd*ga)/(fpi4*LambdaX)
    ctr0r0(2,1,5) = (35.0618960855328_pr*c1*gA2)/fpi4 - &
         (16.8781682906031_pr*c3*gA2)/fpi4 + (0.0543163105867347_pr*c4*gA2)/fpi4 &
         - (0.06328125_pr*cd*ga)/(fpi4*LambdaX)
    ctr0r0(1,2,6) = (-325.455281192602_pr*c1*gA2)/fpi4 + &
         (83.7334611805583_pr*c3*gA2)/fpi4 + (1.3574300625_pr*c4*gA2)/fpi4
    ctr0r0(2,2,6) = (32.9911875098339_pr*c3*gA2)/fpi4 + &
         (0.16875_pr*c4*gA2)/fpi4
    ctr0r0(1,1,7) = (-5761.43588948783_pr*c1*gA2)/fpi4 - &
         (869.283572934387_pr*c3*gA2)/fpi4 - (13.4068326170281_pr*c4*gA2)/fpi4
    ctr0r0(2,1,7) = (2192.10779645869_pr*c1*gA2)/fpi4 + &
         (263.842952909344_pr*c3*gA2)/fpi4 + (3.99388560267857_pr*c4*gA2)/fpi4 - &
         (0.00439453125_pr*cd*ga)/(fpi4*LambdaX)
    ctr0r0(3,1,7) = (-198.153843790565_pr*c1*gA2)/fpi4 - &
         (58.742085950407_pr*c3*gA2)/fpi4 - (0.4904296875_pr*c4*gA2)/fpi4
    ctr0r0(1,3,7) = (-526.213527565944_pr*c3*gA2)/fpi4 - &
         (2.84068526785714_pr*c4*gA2)/fpi4
    ctr0r0(1,2,8) = (11194.2459766741_pr*c1*gA2)/fpi4 - &
         (901.889210563564_pr*c3*gA2)/fpi4 + (1.90516445471939_pr*c4*gA2)/fpi4
    ctr0r0(2,2,8) = (-751.043571428571_pr*c1*gA2)/fpi4 + &
         (1336.68785503856_pr*c3*gA2)/fpi4 + (4.15185267857143_pr*c4*gA2)/fpi4
    ctr0r0(1,1,9) = (-19577.5085802329_pr*c1*gA2)/fpi4 - &
         (20061.3662796935_pr*c3*gA2)/fpi4 - (402.203416845285_pr*c4*gA2)/fpi4
    ctr0r0(2,1,9) = (6927.60941249884_pr*c1*gA2)/fpi4 + &
         (10204.912096572_pr*c3*gA2)/fpi4 + (157.45919188058_pr*c4*gA2)/fpi4
    ctr0r0(3,1,9) = (-1688.44574624066_pr*c1*gA2)/fpi4 - &
         (1848.73711933536_pr*c3*gA2)/fpi4 - (17.251846749442_pr*c4*gA2)/fpi4
    ctr0r0(1,3,9) = (-14639.9109308036_pr*c1*gA2)/fpi4 - &
         (2073.4495546875_pr*c3*gA2)/fpi4 + (12.406921875_pr*c4*gA2)/fpi4
    ctr0r0(1,2,10) = (-6051.26397950415_pr*c1*gA2)/fpi4 + &
         (25318.5711035127_pr*c3*gA2)/fpi4 + (1111.77481741071_pr*c4*gA2)/fpi4
    ctr0r0(2,2,10) = (20263.462438058_pr*c1*gA2)/fpi4 - &
         (645.573171294818_pr*c3*gA2)/fpi4 - (220.10902734375_pr*c4*gA2)/fpi4
    ctr0r0(1,1,11) = (-117612.589511906_pr*c1*gA2)/fpi4 - &
         (22681.6420225835_pr*c3*gA2)/fpi4 + (1494.47871529855_pr*c4*gA2)/fpi4
    ctr0r0(2,1,11) = (78131.8227221231_pr*c1*gA2)/fpi4 + &
         (11306.3195229362_pr*c3*gA2)/fpi4 - (1044.42429449456_pr*c4*gA2)/fpi4
    ctr0r0(3,1,11) = (-19244.1407277278_pr*c1*gA2)/fpi4 - &
         (4031.58809797695_pr*c3*gA2)/fpi4 + (143.952209472656_pr*c4*gA2)/fpi4
    ctr0r0(1,3,11) = (-12782.365790625_pr*c1*gA2)/fpi4 - &
         (19898.5801884375_pr*c3*gA2)/fpi4 - (763.061715_pr*c4*gA2)/fpi4
    ctr0r0(1,2,12) = (112532.892207631_pr*c1*gA2)/fpi4 + &
         (2673.93101274094_pr*c3*gA2)/fpi4 - (3075.39951940848_pr*c4*gA2)/fpi4
    ctr0r0(2,2,12) = (-4174.98147321429_pr*c1*gA2)/fpi4 + &
         (16959.475747232_pr*c3*gA2)/fpi4 + (1249.48156138393_pr*c4*gA2)/fpi4
    ctr0r0(1,1,13) = (-62220.7516704385_pr*c1*gA2)/fpi4 - &
         (37645.4468013683_pr*c3*gA2)/fpi4 - (1458.88122707611_pr*c4*gA2)/fpi4
    ctr0r0(2,1,13) = (49392.5657264478_pr*c1*gA2)/fpi4 + &
         (35938.6848928339_pr*c3*gA2)/fpi4 + (1614.83538151428_pr*c4*gA2)/fpi4
    ctr0r0(3,1,13) = (-21590.8338635122_pr*c1*gA2)/fpi4 - &
         (12971.4689072793_pr*c3*gA2)/fpi4 - (402.76586151123_pr*c4*gA2)/fpi4
    ctr0r0(1,3,13) = (-66293.2611893973_pr*c1*gA2)/fpi4 - &
         (8811.70601861049_pr*c3*gA2)/fpi4 + (1385.87636601563_pr*c4*gA2)/fpi4
    ctr0r0(1,2,14) = (-8634.9144059954_pr*c1*gA2)/fpi4 + &
         (13584.0523217795_pr*c3*gA2)/fpi4 + (1941.14579905622_pr*c4*gA2)/fpi4
    ctr0r0(2,2,14) = (46092.3086405378_pr*c1*gA2)/fpi4 + &
         (6421.96926354563_pr*c3*gA2)/fpi4 - (1169.36919344657_pr*c4*gA2)/fpi4
    ctr0r0(1,1,15) = (-54823.1734859566_pr*c1*gA2)/fpi4 - &
         (21478.1701524643_pr*c3*gA2)/fpi4 + (30.6758968031093_pr*c4*gA2)/fpi4
    ctr0r0(2,1,15) = (79796.8451847046_pr*c1*gA2)/fpi4 + &
         (29550.3190922257_pr*c3*gA2)/fpi4 - (265.149560200718_pr*c4*gA2)/fpi4
    ctr0r0(3,1,15) = (-34549.4910141965_pr*c1*gA2)/fpi4 - &
         (12487.0361763301_pr*c3*gA2)/fpi4 + (129.866342965262_pr*c4*gA2)/fpi4
    ctr0r0(1,3,15) = (-6841.585265625_pr*c1*gA2)/fpi4 - &
         (8394.50421527972_pr*c3*gA2)/fpi4 - (637.487162214007_pr*c4*gA2)/fpi4
    ctr0r0(1,2,16) = (-966.144667675781_pr*c1*gA2)/fpi4 - &
         (1635.80437776282_pr*c3*gA2)/fpi4 - (20.2837476928711_pr*c4*gA2)/fpi4
    ctr0r0(2,2,16) = (10019.4643315377_pr*c1*gA2)/fpi4 + &
         (5806.54125301064_pr*c3*gA2)/fpi4 + (162.342586296387_pr*c4*gA2)/fpi4
    ctr0r0(1,1,17) = (-27328.2660786696_pr*c1*gA2)/fpi4 - &
         (11536.8965413_pr*c3*gA2)/fpi4 - (0.153939957765852_pr*c4*gA2)/fpi4
    ctr0r0(2,1,17) = (51084.8335772702_pr*c1*gA2)/fpi4 + &
         (21317.5006264206_pr*c3*gA2)/fpi4 - (5.04213459735325_pr*c4*gA2)/fpi4
    ctr0r0(3,1,17) = (-23176.714608858_pr*c1*gA2)/fpi4 - &
         (9591.77377013768_pr*c3*gA2)/fpi4 - (7.01374572372437_pr*c4*gA2)/fpi4
    ctr0r0(1,3,17) = (-839.926762520926_pr*c1*gA2)/fpi4
    ctr0r0(1,2,18) = (-71.6191333138602_pr*c1*gA2)/fpi4 - &
         (29.5562121168625_pr*c3*gA2)/fpi4
    ctr0r0(2,2,18) = (1117.2724960968_pr*c1*gA2)/fpi4 + &
         (423.729200499137_pr*c3*gA2)/fpi4 + (5.07093692321777_pr*c4*gA2)/fpi4
    ctr0r0(1,1,19) = (-9253.22855094582_pr*c1*gA2)/fpi4 - &
         (4005.49423137719_pr*c3*gA2)/fpi4
    ctr0r0(2,1,19) = (21870.1711895189_pr*c1*gA2)/fpi4 + &
         (9310.62118848288_pr*c3*gA2)/fpi4 + (0.0769699788829259_pr*c4*gA2)/fpi4
    ctr0r0(3,1,19) = (-11369.3247545611_pr*c1*gA2)/fpi4 - &
         (4767.97117814042_pr*c3*gA2)/fpi4 - (0.656709900856018_pr*c4*gA2)/fpi4
    ctr0r0(2,2,20) = (17.9047833284651_pr*c1*gA2)/fpi4 + &
         (7.38905302921563_pr*c3*gA2)/fpi4
    ctr0r0(1,1,21) = (-2139.34177995088_pr*c1*gA2)/fpi4 - &
         (951.965547431256_pr*c3*gA2)/fpi4
    ctr0r0(2,1,21) = (6585.92097742191_pr*c1*gA2)/fpi4 + &
         (2872.41393631472_pr*c3*gA2)/fpi4
    ctr0r0(3,1,21) = (-4206.97567815619_pr*c1*gA2)/fpi4 - &
         (1804.69962149452_pr*c3*gA2)/fpi4 - &
         (0.00962124736036573_pr*c4*gA2)/fpi4
    ctr0r0(1,1,23) = (-320.600347402137_pr*c1*gA2)/fpi4 - &
         (147.828203709359_pr*c3*gA2)/fpi4
    ctr0r0(2,1,23) = (1373.03667928303_pr*c1*gA2)/fpi4 + &
         (615.472966679238_pr*c3*gA2)/fpi4
    ctr0r0(3,1,23) = (-1139.80667854591_pr*c1*gA2)/fpi4 - &
         (500.614173373695_pr*c3*gA2)/fpi4
    ctr0r0(1,1,25) = (-28.0192319249075_pr*c1*gA2)/fpi4 - &
         (13.7137479594838_pr*c3*gA2)/fpi4
    ctr0r0(2,1,25) = (187.596940698279_pr*c1*gA2)/fpi4 + &
         (87.2234077809484_pr*c3*gA2)/fpi4
    ctr0r0(3,1,25) = (-216.193884091206_pr*c1*gA2)/fpi4 - &
         (97.5972978305762_pr*c3*gA2)/fpi4
    ctr0r0(1,1,27) = (-1.08369739154545_pr*c1*gA2)/fpi4 - &
         (0.619027285041159_pr*c3*gA2)/fpi4
    ctr0r0(2,1,27) = (15.0933133539992_pr*c1*gA2)/fpi4 + &
         (7.47177985304339_pr*c3*gA2)/fpi4
    ctr0r0(3,1,27) = (-27.1326378098226_pr*c1*gA2)/fpi4 - &
         (12.7192853287926_pr*c3*gA2)/fpi4
    ctr0r0(1,1,29) = (-0.00618211760948692_pr*c3*gA2)/fpi4
    ctr0r0(2,1,29) = (0.541848695772724_pr*c1*gA2)/fpi4 + &
         (0.315695760130066_pr*c3*gA2)/fpi4
    ctr0r0(3,1,29) = (-2.02212634319308_pr*c1*gA2)/fpi4 - &
         (1.01238124519548_pr*c3*gA2)/fpi4
    ctr0r0(2,1,31) = (0.00309105880474346_pr*c3*gA2)/fpi4
    ctr0r0(3,1,31) = (-0.0677310869715905_pr*c1*gA2)/fpi4 - &
         (0.0402347347174442_pr*c3*gA2)/fpi4
    ctr0r0(3,1,33) = (-0.000386382350592933_pr*c3*gA2)/fpi4
    !
    !  Urho1rho1
    !
    ctr1r1=0.0_pr
    ctr1r1(1,1,1) = (-0.964246762695464_pr*c3*gA2)/fpi4 - &
         (0.0263462308309038_pr*c4*gA2)/fpi4 - (0.1875_pr*ce)/(fpi4*LambdaX) + &
         (0.046875_pr*cd*ga)/(fpi4*LambdaX)
    ctr1r1(1,1,3) = (-6.84195361981036_pr*c1*gA2)/fpi4 + &
         (5.01410456224811_pr*c3*gA2)/fpi4 + (0.0732776804846939_pr*c4*gA2)/fpi4 &
         - (0.06328125_pr*cd*ga)/(fpi4*LambdaX)
    ctr1r1(2,1,3) = (0.147698198400372_pr*c3*gA2)/fpi4 + &
         (0.0703526785714286_pr*c4*gA2)/fpi4 - &
         (0.02109375_pr*cd*ga)/(fpi4*LambdaX)
    ctr1r1(1,2,4) = (-0.1875_pr*c3*gA2)/fpi4 - (0.375_pr*c4*gA2)/fpi4 + &
         (0.140625_pr*cd*ga)/(fpi4*LambdaX)
    ctr1r1(1,1,5) = (-23.1434396482781_pr*c1*gA2)/fpi4 + &
         (12.2465870768616_pr*c3*gA2)/fpi4 - (0.18795405731824_pr*c4*gA2)/fpi4 + &
         (0.017578125_pr*cd*ga)/(fpi4*LambdaX)
    ctr1r1(2,1,5) = (35.0618960855328_pr*c1*gA2)/fpi4 - &
         (16.8781682906031_pr*c3*gA2)/fpi4 + (0.0543163105867347_pr*c4*gA2)/fpi4 &
         - (0.06328125_pr*cd*ga)/(fpi4*LambdaX)
    ctr1r1(1,2,6) = (-325.455281192602_pr*c1*gA2)/fpi4 + &
         (83.7334611805583_pr*c3*gA2)/fpi4 + (1.3574300625_pr*c4*gA2)/fpi4
    ctr1r1(2,2,6) = (32.9911875098339_pr*c3*gA2)/fpi4 + &
         (0.16875_pr*c4*gA2)/fpi4
    ctr1r1(1,1,7) = (-5761.43588948783_pr*c1*gA2)/fpi4 - &
         (869.283572934387_pr*c3*gA2)/fpi4 - (13.4068326170281_pr*c4*gA2)/fpi4
    ctr1r1(2,1,7) = (2192.10779645869_pr*c1*gA2)/fpi4 + &
         (263.842952909344_pr*c3*gA2)/fpi4 + (3.99388560267857_pr*c4*gA2)/fpi4 - &
         (0.00439453125_pr*cd*ga)/(fpi4*LambdaX)
    ctr1r1(3,1,7) = (-198.153843790565_pr*c1*gA2)/fpi4 - &
         (58.742085950407_pr*c3*gA2)/fpi4 - (0.4904296875_pr*c4*gA2)/fpi4
    ctr1r1(1,3,7) = (-526.213527565944_pr*c3*gA2)/fpi4 - &
         (2.84068526785714_pr*c4*gA2)/fpi4
    ctr1r1(1,2,8) = (11194.2459766741_pr*c1*gA2)/fpi4 - &
         (901.889210563564_pr*c3*gA2)/fpi4 + (1.90516445471939_pr*c4*gA2)/fpi4
    ctr1r1(2,2,8) = (-751.043571428571_pr*c1*gA2)/fpi4 + &
         (1336.68785503856_pr*c3*gA2)/fpi4 + (4.15185267857143_pr*c4*gA2)/fpi4
    ctr1r1(1,1,9) = (-19577.5085802329_pr*c1*gA2)/fpi4 - &
         (20061.3662796935_pr*c3*gA2)/fpi4 - (402.203416845285_pr*c4*gA2)/fpi4
    ctr1r1(2,1,9) = (6927.60941249884_pr*c1*gA2)/fpi4 + &
         (10204.912096572_pr*c3*gA2)/fpi4 + (157.45919188058_pr*c4*gA2)/fpi4
    ctr1r1(3,1,9) = (-1688.44574624066_pr*c1*gA2)/fpi4 - &
         (1848.73711933536_pr*c3*gA2)/fpi4 - (17.251846749442_pr*c4*gA2)/fpi4
    ctr1r1(1,3,9) = (-14639.9109308036_pr*c1*gA2)/fpi4 - &
         (2073.4495546875_pr*c3*gA2)/fpi4 + (12.406921875_pr*c4*gA2)/fpi4
    ctr1r1(1,2,10) = (-6051.26397950415_pr*c1*gA2)/fpi4 + &
         (25318.5711035127_pr*c3*gA2)/fpi4 + (1111.77481741071_pr*c4*gA2)/fpi4
    ctr1r1(2,2,10) = (20263.462438058_pr*c1*gA2)/fpi4 - &
         (645.573171294818_pr*c3*gA2)/fpi4 - (220.10902734375_pr*c4*gA2)/fpi4
    ctr1r1(1,1,11) = (-117612.589511906_pr*c1*gA2)/fpi4 - &
         (22681.6420225835_pr*c3*gA2)/fpi4 + (1494.47871529855_pr*c4*gA2)/fpi4
    ctr1r1(2,1,11) = (78131.8227221231_pr*c1*gA2)/fpi4 + &
         (11306.3195229362_pr*c3*gA2)/fpi4 - (1044.42429449456_pr*c4*gA2)/fpi4
    ctr1r1(3,1,11) = (-19244.1407277278_pr*c1*gA2)/fpi4 - &
         (4031.58809797695_pr*c3*gA2)/fpi4 + (143.952209472656_pr*c4*gA2)/fpi4
    ctr1r1(1,3,11) = (-12782.365790625_pr*c1*gA2)/fpi4 - &
         (19898.5801884375_pr*c3*gA2)/fpi4 - (763.061715_pr*c4*gA2)/fpi4
    ctr1r1(1,2,12) = (112532.892207631_pr*c1*gA2)/fpi4 + &
         (2673.93101274094_pr*c3*gA2)/fpi4 - (3075.39951940848_pr*c4*gA2)/fpi4
    ctr1r1(2,2,12) = (-4174.98147321429_pr*c1*gA2)/fpi4 + &
         (16959.475747232_pr*c3*gA2)/fpi4 + (1249.48156138393_pr*c4*gA2)/fpi4
    ctr1r1(1,1,13) = (-62220.7516704385_pr*c1*gA2)/fpi4 - &
         (37645.4468013683_pr*c3*gA2)/fpi4 - (1458.88122707611_pr*c4*gA2)/fpi4
    ctr1r1(2,1,13) = (49392.5657264478_pr*c1*gA2)/fpi4 + &
         (35938.6848928339_pr*c3*gA2)/fpi4 + (1614.83538151428_pr*c4*gA2)/fpi4
    ctr1r1(3,1,13) = (-21590.8338635122_pr*c1*gA2)/fpi4 - &
         (12971.4689072793_pr*c3*gA2)/fpi4 - (402.76586151123_pr*c4*gA2)/fpi4
    ctr1r1(1,3,13) = (-66293.2611893973_pr*c1*gA2)/fpi4 - &
         (8811.70601861049_pr*c3*gA2)/fpi4 + (1385.87636601563_pr*c4*gA2)/fpi4
    ctr1r1(1,2,14) = (-8634.9144059954_pr*c1*gA2)/fpi4 + &
         (13584.0523217795_pr*c3*gA2)/fpi4 + (1941.14579905622_pr*c4*gA2)/fpi4
    ctr1r1(2,2,14) = (46092.3086405378_pr*c1*gA2)/fpi4 + &
         (6421.96926354563_pr*c3*gA2)/fpi4 - (1169.36919344657_pr*c4*gA2)/fpi4
    ctr1r1(1,1,15) = (-54823.1734859566_pr*c1*gA2)/fpi4 - &
         (21478.1701524643_pr*c3*gA2)/fpi4 + (30.6758968031093_pr*c4*gA2)/fpi4
    ctr1r1(2,1,15) = (79796.8451847046_pr*c1*gA2)/fpi4 + &
         (29550.3190922257_pr*c3*gA2)/fpi4 - (265.149560200718_pr*c4*gA2)/fpi4
    ctr1r1(3,1,15) = (-34549.4910141965_pr*c1*gA2)/fpi4 - &
         (12487.0361763301_pr*c3*gA2)/fpi4 + (129.866342965262_pr*c4*gA2)/fpi4
    ctr1r1(1,3,15) = (-6841.585265625_pr*c1*gA2)/fpi4 - &
         (8394.50421527972_pr*c3*gA2)/fpi4 - (637.487162214007_pr*c4*gA2)/fpi4
    ctr1r1(1,2,16) = (-966.144667675781_pr*c1*gA2)/fpi4 - &
         (1635.80437776282_pr*c3*gA2)/fpi4 - (20.2837476928711_pr*c4*gA2)/fpi4
    ctr1r1(2,2,16) = (10019.4643315377_pr*c1*gA2)/fpi4 + &
         (5806.54125301064_pr*c3*gA2)/fpi4 + (162.342586296387_pr*c4*gA2)/fpi4
    ctr1r1(1,1,17) = (-27328.2660786696_pr*c1*gA2)/fpi4 - &
         (11536.8965413_pr*c3*gA2)/fpi4 - (0.153939957765852_pr*c4*gA2)/fpi4
    ctr1r1(2,1,17) = (51084.8335772702_pr*c1*gA2)/fpi4 + &
         (21317.5006264206_pr*c3*gA2)/fpi4 - (5.04213459735325_pr*c4*gA2)/fpi4
    ctr1r1(3,1,17) = (-23176.714608858_pr*c1*gA2)/fpi4 - &
         (9591.77377013768_pr*c3*gA2)/fpi4 - (7.01374572372437_pr*c4*gA2)/fpi4
    ctr1r1(1,3,17) = (-839.926762520926_pr*c1*gA2)/fpi4
    ctr1r1(1,2,18) = (-71.6191333138602_pr*c1*gA2)/fpi4 - &
         (29.5562121168625_pr*c3*gA2)/fpi4
    ctr1r1(2,2,18) = (1117.2724960968_pr*c1*gA2)/fpi4 + &
         (423.729200499137_pr*c3*gA2)/fpi4 + (5.07093692321777_pr*c4*gA2)/fpi4
    ctr1r1(1,1,19) = (-9253.22855094582_pr*c1*gA2)/fpi4 - &
         (4005.49423137719_pr*c3*gA2)/fpi4
    ctr1r1(2,1,19) = (21870.1711895189_pr*c1*gA2)/fpi4 + &
         (9310.62118848288_pr*c3*gA2)/fpi4 + (0.0769699788829259_pr*c4*gA2)/fpi4
    ctr1r1(3,1,19) = (-11369.3247545611_pr*c1*gA2)/fpi4 - &
         (4767.97117814042_pr*c3*gA2)/fpi4 - (0.656709900856018_pr*c4*gA2)/fpi4
    ctr1r1(2,2,20) = (17.9047833284651_pr*c1*gA2)/fpi4 + &
         (7.38905302921563_pr*c3*gA2)/fpi4
    ctr1r1(1,1,21) = (-2139.34177995088_pr*c1*gA2)/fpi4 - &
         (951.965547431256_pr*c3*gA2)/fpi4
    ctr1r1(2,1,21) = (6585.92097742191_pr*c1*gA2)/fpi4 + &
         (2872.41393631472_pr*c3*gA2)/fpi4
    ctr1r1(3,1,21) = (-4206.97567815619_pr*c1*gA2)/fpi4 - &
         (1804.69962149452_pr*c3*gA2)/fpi4 - &
         (0.00962124736036573_pr*c4*gA2)/fpi4
    ctr1r1(1,1,23) = (-320.600347402137_pr*c1*gA2)/fpi4 - &
         (147.828203709359_pr*c3*gA2)/fpi4
    ctr1r1(2,1,23) = (1373.03667928303_pr*c1*gA2)/fpi4 + &
         (615.472966679238_pr*c3*gA2)/fpi4
    ctr1r1(3,1,23) = (-1139.80667854591_pr*c1*gA2)/fpi4 - &
         (500.614173373695_pr*c3*gA2)/fpi4
    ctr1r1(1,1,25) = (-28.0192319249075_pr*c1*gA2)/fpi4 - &
         (13.7137479594838_pr*c3*gA2)/fpi4
    ctr1r1(2,1,25) = (187.596940698279_pr*c1*gA2)/fpi4 + &
         (87.2234077809484_pr*c3*gA2)/fpi4
    ctr1r1(3,1,25) = (-216.193884091206_pr*c1*gA2)/fpi4 - &
         (97.5972978305762_pr*c3*gA2)/fpi4
    ctr1r1(1,1,27) = (-1.08369739154545_pr*c1*gA2)/fpi4 - &
         (0.619027285041159_pr*c3*gA2)/fpi4
    ctr1r1(2,1,27) = (15.0933133539992_pr*c1*gA2)/fpi4 + &
         (7.47177985304339_pr*c3*gA2)/fpi4
    ctr1r1(3,1,27) = (-27.1326378098226_pr*c1*gA2)/fpi4 - &
         (12.7192853287926_pr*c3*gA2)/fpi4
    ctr1r1(1,1,29) = (-0.00618211760948692_pr*c3*gA2)/fpi4
    ctr1r1(2,1,29) = (0.541848695772724_pr*c1*gA2)/fpi4 + &
         (0.315695760130066_pr*c3*gA2)/fpi4
    ctr1r1(3,1,29) = (-2.02212634319308_pr*c1*gA2)/fpi4 - &
         (1.01238124519548_pr*c3*gA2)/fpi4
    ctr1r1(2,1,31) = (0.00309105880474346_pr*c3*gA2)/fpi4
    ctr1r1(3,1,31) = (-0.0677310869715905_pr*c1*gA2)/fpi4 - &
         (0.0402347347174442_pr*c3*gA2)/fpi4
    ctr1r1(3,1,33) = (-0.000386382350592933_pr*c3*gA2)/fpi4
    !
    !  UDrho0Drho0
    !
    ctdr0dr0=0.0_pr
    !
    !  UDrho1Drho1
    !
    ctdr1dr1=0.0_pr
    !
    !  Urho0DDrho0
    !
    ctr0Dr0=0.0_pr
    ctr0Dr0(1,1,3) = (-0.682809650873706_pr*c3*gA2)/(fpi4*mpi2) + &
         (0.0239679050048591_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,1,5) = (-5.25041224998774_pr*c1*gA2)/(fpi4*mpi2) + &
         (3.77245681259736_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.166505020727041_pr*c4*gA2)/(fpi4*mpi2) + &
         (0.03515625_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctr0Dr0(2,1,5) = (0.12308183200031_pr*c3*gA2)/(fpi4*mpi2) + &
         (0.0586272321428571_pr*c4*gA2)/(fpi4*mpi2) - &
         (0.017578125_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctr0Dr0(1,1,7) = (-26.0336741183036_pr*c1*gA2)/(fpi4*mpi2) + &
         (12.9213033479161_pr*c3*gA2)/(fpi4*mpi2) + &
         (0.0724766691246811_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,1,7) = (29.218246737944_pr*c1*gA2)/(fpi4*mpi2) - &
         (13.9991785771689_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.0572336176658163_pr*c4*gA2)/(fpi4*mpi2) - &
         (0.0087890625_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctr0Dr0(1,2,8) = (-252.855350558036_pr*c1*gA2)/(fpi4*mpi2) + &
         (56.1626339963032_pr*c3*gA2)/(fpi4*mpi2) + &
         (0.49978546875_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,2,8) = (27.4926562581949_pr*c3*gA2)/(fpi4*mpi2) + &
         (0.140625_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,1,9) = (-4737.41799564353_pr*c1*gA2)/(fpi4*mpi2) - &
         (732.936270957556_pr*c3*gA2)/(fpi4*mpi2) - &
         (10.8826855699936_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,1,9) = (1812.39349332203_pr*c1*gA2)/(fpi4*mpi2) + &
         (230.621448631588_pr*c3*gA2)/(fpi4*mpi2) + &
         (3.48044224330357_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(3,1,9) = (-165.128203158804_pr*c1*gA2)/(fpi4*mpi2) - &
         (48.9517382920058_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.40869140625_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,3,9) = (-412.495287388393_pr*c3*gA2)/(fpi4*mpi2) - &
         (1.89848772321429_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,2,10) = (8987.88523861607_pr*c1*gA2)/(fpi4*mpi2) - &
         (631.430929909134_pr*c3*gA2)/(fpi4*mpi2) + &
         (1.32764820631378_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,2,10) = (-625.869642857143_pr*c1*gA2)/(fpi4*mpi2) + &
         (1045.17490521998_pr*c3*gA2)/(fpi4*mpi2) + &
         (3.10831473214286_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,1,11) = (-17022.0529432446_pr*c1*gA2)/(fpi4*mpi2) - &
         (15425.6442847704_pr*c3*gA2)/(fpi4*mpi2) - &
         (325.368422118792_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,1,11) = (6144.54246637738_pr*c1*gA2)/(fpi4*mpi2) + &
         (7927.36933047731_pr*c3*gA2)/(fpi4*mpi2) + &
         (127.556172677176_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(3,1,11) = (-1407.03812186722_pr*c1*gA2)/(fpi4*mpi2) - &
         (1455.27334812581_pr*c3*gA2)/(fpi4*mpi2) - &
         (13.9398324148996_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,3,11) = (-11514.015703125_pr*c1*gA2)/(fpi4*mpi2) - &
         (1727.87462890625_pr*c3*gA2)/(fpi4*mpi2) + &
         (10.3391015625_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,2,12) = (-694.446624441964_pr*c1*gA2)/(fpi4*mpi2) + &
         (18247.3109917788_pr*c3*gA2)/(fpi4*mpi2) + &
         (899.399912946429_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,2,12) = (15092.8155175781_pr*c1*gA2)/(fpi4*mpi2) - &
         (116.351373250285_pr*c3*gA2)/(fpi4*mpi2) - &
         (178.470087890625_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,1,13) = (-84816.595261174_pr*c1*gA2)/(fpi4*mpi2) - &
         (21005.6538813155_pr*c3*gA2)/(fpi4*mpi2) + &
         (1207.97363186907_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,1,13) = (56932.8731897877_pr*c1*gA2)/(fpi4*mpi2) + &
         (11075.3299701944_pr*c3*gA2)/(fpi4*mpi2) - &
         (844.741128183943_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(3,1,13) = (-14559.3117129491_pr*c1*gA2)/(fpi4*mpi2) - &
         (3504.49103169193_pr*c3*gA2)/(fpi4*mpi2) + &
         (116.596984863281_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,3,13) = (-12860.985515625_pr*c1*gA2)/(fpi4*mpi2) - &
         (14398.998496875_pr*c3*gA2)/(fpi4*mpi2) - &
         (616.763278125_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,2,14) = (71603.420010498_pr*c1*gA2)/(fpi4*mpi2) + &
         (8297.07769839308_pr*c3*gA2)/(fpi4*mpi2) - &
         (2483.1196515904_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,2,14) = (2438.09865025112_pr*c1*gA2)/(fpi4*mpi2) + &
         (10960.7457584461_pr*c3*gA2)/(fpi4*mpi2) + &
         (1009.19850167411_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,1,15) = (-55860.3376653923_pr*c1*gA2)/(fpi4*mpi2) - &
         (25835.4946732797_pr*c3*gA2)/(fpi4*mpi2) - &
         (1174.87146717911_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,1,15) = (47592.1515045526_pr*c1*gA2)/(fpi4*mpi2) + &
         (24786.4294824231_pr*c3*gA2)/(fpi4*mpi2) + &
         (1301.50458021698_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(3,1,15) = (-19018.8059498995_pr*c1*gA2)/(fpi4*mpi2) - &
         (9429.01449770894_pr*c3*gA2)/(fpi4*mpi2) - &
         (324.792194366455_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,3,15) = (-44074.7843378906_pr*c1*gA2)/(fpi4*mpi2) - &
         (9754.29660498047_pr*c3*gA2)/(fpi4*mpi2) + &
         (1118.35204980469_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,2,16) = (1466.81527313232_pr*c1*gA2)/(fpi4*mpi2) + &
         (4924.95415195852_pr*c3*gA2)/(fpi4*mpi2) + &
         (1563.11953388149_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,2,16) = (29864.464335022_pr*c1*gA2)/(fpi4*mpi2) + &
         (7829.64002841679_pr*c3*gA2)/(fpi4*mpi2) - &
         (942.380516836984_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,1,17) = (-42578.4126891706_pr*c1*gA2)/(fpi4*mpi2) - &
         (17202.0887780327_pr*c3*gA2)/(fpi4*mpi2) + &
         (24.7012874970572_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,1,17) = (61948.0344894364_pr*c1*gA2)/(fpi4*mpi2) + &
         (24461.7278379897_pr*c3*gA2)/(fpi4*mpi2) - &
         (213.463917761666_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(3,1,17) = (-26551.291288457_pr*c1*gA2)/(fpi4*mpi2) - &
         (10248.0190226489_pr*c3*gA2)/(fpi4*mpi2) + &
         (104.612325701032_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,3,17) = (-7839.31645019531_pr*c1*gA2)/(fpi4*mpi2) - &
         (4409.58359893799_pr*c3*gA2)/(fpi4*mpi2) - &
         (513.293134852818_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,2,18) = (-2406.80499142456_pr*c1*gA2)/(fpi4*mpi2) - &
         (982.897425834421_pr*c3*gA2)/(fpi4*mpi2) - &
         (16.3321086730957_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,2,18) = (8676.01512762451_pr*c1*gA2)/(fpi4*mpi2) + &
         (3669.79377386557_pr*c3*gA2)/(fpi4*mpi2) + &
         (130.679305718994_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,1,19) = (-21929.7613258823_pr*c1*gA2)/(fpi4*mpi2) - &
         (9003.06173155003_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.123949706677028_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,1,19) = (41084.26511068_pr*c1*gA2)/(fpi4*mpi2) + &
         (16712.9487662382_pr*c3*gA2)/(fpi4*mpi2) - &
         (4.06063970530374_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(3,1,19) = (-18477.8744897258_pr*c1*gA2)/(fpi4*mpi2) - &
         (7453.07953912234_pr*c3*gA2)/(fpi4*mpi2) - &
         (5.64185216903687_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,2,20) = (601.70124785614_pr*c1*gA2)/(fpi4*mpi2) + &
         (245.724356458605_pr*c3*gA2)/(fpi4*mpi2) + &
         (4.08302716827393_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,1,21) = (-7429.54966559941_pr*c1*gA2)/(fpi4*mpi2) - &
         (3087.22338867281_pr*c3*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,1,21) = (17522.7958314409_pr*c1*gA2)/(fpi4*mpi2) + &
         (7223.97599005467_pr*c3*gA2)/(fpi4*mpi2) + &
         (0.0619748533385141_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(3,1,21) = (-9051.54263718816_pr*c1*gA2)/(fpi4*mpi2) - &
         (3700.96455501676_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.528670542240143_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,1,23) = (-1722.74624112667_pr*c1*gA2)/(fpi4*mpi2) - &
         (722.002478964633_pr*c3*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,1,23) = (5292.28471362693_pr*c1*gA2)/(fpi4*mpi2) + &
         (2204.46546742604_pr*c3*gA2)/(fpi4*mpi2)
    ctr0Dr0(3,1,23) = (-3370.28377443923_pr*c1*gA2)/(fpi4*mpi2) - &
         (1394.125888934_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.00774685666731426_pr*c4*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,1,25) = (-258.704980303594_pr*c1*gA2)/(fpi4*mpi2) - &
         (108.950617217712_pr*c3*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,1,25) = (1106.14939761973_pr*c1*gA2)/(fpi4*mpi2) + &
         (464.078201799952_pr*c3*gA2)/(fpi4*mpi2)
    ctr0Dr0(3,1,25) = (-916.533683372315_pr*c1*gA2)/(fpi4*mpi2) - &
         (382.506711741246_pr*c3*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,1,27) = (-22.6474359092818_pr*c1*gA2)/(fpi4*mpi2) - &
         (9.54811915326822_pr*c3*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,1,27) = (151.415132381586_pr*c1*gA2)/(fpi4*mpi2) + &
         (63.7775488277401_pr*c3*gA2)/(fpi4*mpi2)
    ctr0Dr0(3,1,27) = (-174.235171472089_pr*c1*gA2)/(fpi4*mpi2) - &
         (73.1584858358234_pr*c3*gA2)/(fpi4*mpi2)
    ctr0Dr0(1,1,29) = (-0.877190519238951_pr*c1*gA2)/(fpi4*mpi2) - &
         (0.368818401576571_pr*c3*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,1,29) = (12.2009084738798_pr*c1*gA2)/(fpi4*mpi2) + &
         (5.14287797821068_pr*c3*gA2)/(fpi4*mpi2)
    ctr0Dr0(3,1,29) = (-21.9040194562317_pr*c1*gA2)/(fpi4*mpi2) - &
         (9.22717823122213_pr*c3*gA2)/(fpi4*mpi2)
    ctr0Dr0(2,1,31) = (0.438595259619476_pr*c1*gA2)/(fpi4*mpi2) + &
         (0.184409200788286_pr*c3*gA2)/(fpi4*mpi2)
    ctr0Dr0(3,1,31) = (-1.63476237413985_pr*c1*gA2)/(fpi4*mpi2) - &
         (0.688962047473406_pr*c3*gA2)/(fpi4*mpi2)
    ctr0Dr0(3,1,33) = (-0.0548244074524344_pr*c1*gA2)/(fpi4*mpi2) - &
         (0.0230511500985357_pr*c3*gA2)/(fpi4*mpi2)
    !
    !  Urho1DDrho0
    !
    ctr1Dr0=0.0_pr
    ctr1Dr0(1,1,3) = (0.0185237085762877_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.00798930166828636_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(1,1,5) = (0.0956155117984694_pr*c1*gA2)/(fpi4*mpi2) - &
         (0.0234414046556122_pr*c3*gA2)/(fpi4*mpi2) + &
         (0.0555016735756803_pr*c4*gA2)/(fpi4*mpi2) - &
         (0.01171875_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctr1Dr0(2,1,5) = (-0.00389508928571429_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.0195424107142857_pr*c4*gA2)/(fpi4*mpi2) + &
         (0.005859375_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctr1Dr0(1,1,7) = (-0.382716029575893_pr*c1*gA2)/(fpi4*mpi2) + &
         (0.159635060586735_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.024158889708227_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(2,1,7) = (-0.0256136798469388_pr*c1*gA2)/(fpi4*mpi2) + &
         (0.00456883769132653_pr*c3*gA2)/(fpi4*mpi2) + &
         (0.0190778725552721_pr*c4*gA2)/(fpi4*mpi2) + &
         (0.0029296875_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctr1Dr0(1,2,8) = (0.706878683035714_pr*c1*gA2)/(fpi4*mpi2) - &
         (0.423425189732143_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.16659515625_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(2,2,8) = (-0.0234375_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.046875_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(1,1,9) = (11.4067307059152_pr*c1*gA2)/(fpi4*mpi2) + &
         (0.609889416055485_pr*c3*gA2)/(fpi4*mpi2) + &
         (3.62756185666454_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(2,1,9) = (-3.54638253348214_pr*c1*gA2)/(fpi4*mpi2) - &
         (0.0995982142857143_pr*c3*gA2)/(fpi4*mpi2) - &
         (1.16014741443452_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(3,1,9) = (0.25048828125_pr*c1*gA2)/(fpi4*mpi2) + &
         (0.04833984375_pr*c3*gA2)/(fpi4*mpi2) + &
         (0.13623046875_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(1,3,9) = (0.949243861607143_pr*c3*gA2)/(fpi4*mpi2) + &
         (0.632829241071429_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(1,2,10) = (-22.8885294642857_pr*c1*gA2)/(fpi4*mpi2) + &
         (18.9939639150191_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.442549402104592_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(2,2,10) = (2.53108258928571_pr*c1*gA2)/(fpi4*mpi2) - &
         (4.85298549107143_pr*c3*gA2)/(fpi4*mpi2) - &
         (1.03610491071429_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(1,1,11) = (220.437076035854_pr*c1*gA2)/(fpi4*mpi2) + &
         (298.906979832589_pr*c3*gA2)/(fpi4*mpi2) + &
         (108.456140706264_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(2,1,11) = (-68.8102811802455_pr*c1*gA2)/(fpi4*mpi2) - &
         (121.491437360491_pr*c3*gA2)/(fpi4*mpi2) - &
         (42.5187242257254_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(3,1,11) = (6.13759068080357_pr*c1*gA2)/(fpi4*mpi2) + &
         (12.8441554478237_pr*c3*gA2)/(fpi4*mpi2) + &
         (4.64661080496652_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(1,3,11) = (11.318203125_pr*c1*gA2)/(fpi4*mpi2) - &
         (23.38048828125_pr*c3*gA2)/(fpi4*mpi2) - &
         (3.4463671875_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(1,2,12) = (-732.624496372768_pr*c1*gA2)/(fpi4*mpi2) - &
         (842.4480234375_pr*c3*gA2)/(fpi4*mpi2) - &
         (299.799970982143_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(2,2,12) = (117.953701171875_pr*c1*gA2)/(fpi4*mpi2) + &
         (175.760551757812_pr*c3*gA2)/(fpi4*mpi2) + &
         (59.490029296875_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(1,1,13) = (-1470.77295530744_pr*c1*gA2)/(fpi4*mpi2) - &
         (1177.66729674159_pr*c3*gA2)/(fpi4*mpi2) - &
         (402.65787728969_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(2,1,13) = (906.809157364328_pr*c1*gA2)/(fpi4*mpi2) + &
         (809.938923974958_pr*c3*gA2)/(fpi4*mpi2) + &
         (281.580376061314_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(3,1,13) = (-110.580871582031_pr*c1*gA2)/(fpi4*mpi2) - &
         (112.623596191406_pr*c3*gA2)/(fpi4*mpi2) - &
         (38.8656616210938_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(1,3,13) = (585.845296875_pr*c1*gA2)/(fpi4*mpi2) + &
         (588.08008125_pr*c3*gA2)/(fpi4*mpi2) + &
         (205.587759375_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(1,2,14) = (3251.26034692383_pr*c1*gA2)/(fpi4*mpi2) + &
         (2508.40495063477_pr*c3*gA2)/(fpi4*mpi2) + &
         (827.706550530134_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(2,2,14) = (-1184.58378557478_pr*c1*gA2)/(fpi4*mpi2) - &
         (995.956718052455_pr*c3*gA2)/(fpi4*mpi2) - &
         (336.399500558036_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(1,1,15) = (1877.99099063873_pr*c1*gA2)/(fpi4*mpi2) + &
         (1378.17489546073_pr*c3*gA2)/(fpi4*mpi2) + &
         (391.623822393036_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(2,1,15) = (-1930.62984776132_pr*c1*gA2)/(fpi4*mpi2) - &
         (1432.9768855957_pr*c3*gA2)/(fpi4*mpi2) - &
         (433.834860072327_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(3,1,15) = (440.209643554687_pr*c1*gA2)/(fpi4*mpi2) + &
         (341.656391143799_pr*c3*gA2)/(fpi4*mpi2) + &
         (108.264064788818_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(1,3,15) = (-1549.60388085937_pr*c1*gA2)/(fpi4*mpi2) - &
         (1159.75179052734_pr*c3*gA2)/(fpi4*mpi2) - &
         (372.784016601562_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(1,2,16) = (-2582.26988983154_pr*c1*gA2)/(fpi4*mpi2) - &
         (1875.22398799962_pr*c3*gA2)/(fpi4*mpi2) - &
         (521.039844627162_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(2,2,16) = (1459.93622772217_pr*c1*gA2)/(fpi4*mpi2) + &
         (1058.4680922808_pr*c3*gA2)/(fpi4*mpi2) + &
         (314.126838945661_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(1,1,17) = (-56.9118636474609_pr*c1*gA2)/(fpi4*mpi2) - &
         (35.9791338062831_pr*c3*gA2)/(fpi4*mpi2) - &
         (8.23376249901908_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(2,1,17) = (396.575549217224_pr*c1*gA2)/(fpi4*mpi2) + &
         (276.498193133763_pr*c3*gA2)/(fpi4*mpi2) + &
         (71.1546392538888_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(3,1,17) = (-182.383999263218_pr*c1*gA2)/(fpi4*mpi2) - &
         (125.990329873221_pr*c3*gA2)/(fpi4*mpi2) - &
         (34.8707752336775_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(1,3,17) = (871.035161132812_pr*c1*gA2)/(fpi4*mpi2) + &
         (629.940514133998_pr*c3*gA2)/(fpi4*mpi2) + &
         (171.097711617606_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(1,2,18) = (38.1077882995605_pr*c1*gA2)/(fpi4*mpi2) + &
         (24.4976645050049_pr*c3*gA2)/(fpi4*mpi2) + &
         (5.44403622436523_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(2,2,18) = (-244.521582824707_pr*c1*gA2)/(fpi4*mpi2) - &
         (173.383349386597_pr*c3*gA2)/(fpi4*mpi2) - &
         (43.559768572998_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(1,1,19) = (0.340248109817505_pr*c1*gA2)/(fpi4*mpi2) + &
         (0.240601839556013_pr*c3*gA2)/(fpi4*mpi2) + &
         (0.0413165688923427_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(2,1,19) = (9.06178956413269_pr*c1*gA2)/(fpi4*mpi2) + &
         (5.50013281108311_pr*c3*gA2)/(fpi4*mpi2) + &
         (1.35354656843458_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(3,1,19) = (10.8524751663208_pr*c1*gA2)/(fpi4*mpi2) + &
         (8.40921586990356_pr*c3*gA2)/(fpi4*mpi2) + &
         (1.88061738967896_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(2,2,20) = (-9.52694707489014_pr*c1*gA2)/(fpi4*mpi2) - &
         (6.12441612625122_pr*c3*gA2)/(fpi4*mpi2) - &
         (1.36100905609131_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(2,1,21) = (-0.170124054908752_pr*c1*gA2)/(fpi4*mpi2) - &
         (0.120300919778006_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.0206582844461714_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(3,1,21) = (1.29154408693314_pr*c1*gA2)/(fpi4*mpi2) + &
         (0.873662660121918_pr*c3*gA2)/(fpi4*mpi2) + &
         (0.176223514080048_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr0(3,1,23) = (0.0212655068635941_pr*c1*gA2)/(fpi4*mpi2) + &
         (0.0150376149722508_pr*c3*gA2)/(fpi4*mpi2) + &
         (0.00258228555577142_pr*c4*gA2)/(fpi4*mpi2)
    !
    !  Urho1DDrho1
    !
    ctr1Dr1=0.0_pr
    ctr1Dr1(1,1,3) = (0.225058111718187_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.0159786033365727_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,1,5) = (1.76926855346105_pr*c1*gA2)/(fpi4*mpi2) - &
         (1.25129754669487_pr*c3*gA2)/(fpi4*mpi2) + &
         (0.111003347151361_pr*c4*gA2)/(fpi4*mpi2) - &
         (0.0234375_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctr1Dr1(2,1,5) = (-0.0449223666191509_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.0390848214285714_pr*c4*gA2)/(fpi4*mpi2) + &
         (0.01171875_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctr1Dr1(1,1,7) = (8.62173784319196_pr*c1*gA2)/(fpi4*mpi2) - &
         (4.30128936533428_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.0483177794164541_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,1,7) = (-9.7579980091616_pr*c1*gA2)/(fpi4*mpi2) + &
         (4.67405890125441_pr*c3*gA2)/(fpi4*mpi2) + &
         (0.0381557451105442_pr*c4*gA2)/(fpi4*mpi2) + &
         (0.005859375_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctr1Dr1(1,2,8) = (84.4669955357143_pr*c1*gA2)/(fpi4*mpi2) - &
         (18.660361268857_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.3331903125_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,2,8) = (-9.18765625273164_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.09375_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,1,9) = (1573.18444149334_pr*c1*gA2)/(fpi4*mpi2) + &
         (246.450668847167_pr*c3*gA2)/(fpi4*mpi2) + &
         (7.25512371332908_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,1,9) = (-602.738093849159_pr*c1*gA2)/(fpi4*mpi2) - &
         (77.5124918057675_pr*c3*gA2)/(fpi4*mpi2) - &
         (2.32029482886905_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(3,1,9) = (54.959238292518_pr*c1*gA2)/(fpi4*mpi2) + &
         (16.3655859410853_pr*c3*gA2)/(fpi4*mpi2) + &
         (0.2724609375_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,3,9) = (137.182014508929_pr*c3*gA2)/(fpi4*mpi2) + &
         (1.26565848214286_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,2,10) = (-2950.89558816964_pr*c1*gA2)/(fpi4*mpi2) + &
         (192.771436581372_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.885098804209184_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,2,10) = (202.294921875_pr*c1*gA2)/(fpi4*mpi2) - &
         (345.017723242968_pr*c3*gA2)/(fpi4*mpi2) - &
         (2.07220982142857_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,1,11) = (6073.45819058092_pr*c1*gA2)/(fpi4*mpi2) + &
         (4996.23922290246_pr*c3*gA2)/(fpi4*mpi2) + &
         (216.912281412528_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,1,11) = (-2212.39051736854_pr*c1*gA2)/(fpi4*mpi2) - &
         (2578.84096274467_pr*c3*gA2)/(fpi4*mpi2) - &
         (85.0374484514509_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(3,1,11) = (484.396391719877_pr*c1*gA2)/(fpi4*mpi2) + &
         (478.331459408286_pr*c3*gA2)/(fpi4*mpi2) + &
         (9.29322160993304_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,3,11) = (3776.4796875_pr*c1*gA2)/(fpi4*mpi2) + &
         (595.892330729167_pr*c3*gA2)/(fpi4*mpi2) - &
         (6.892734375_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,2,12) = (-1047.22431947545_pr*c1*gA2)/(fpi4*mpi2) - &
         (5677.84872382209_pr*c3*gA2)/(fpi4*mpi2) - &
         (599.599941964286_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,2,12) = (-4730.6293359375_pr*c1*gA2)/(fpi4*mpi2) - &
         (56.7874345025092_pr*c3*gA2)/(fpi4*mpi2) + &
         (118.98005859375_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,1,13) = (26158.9112133749_pr*c1*gA2)/(fpi4*mpi2) + &
         (7519.86242027358_pr*c3*gA2)/(fpi4*mpi2) - &
         (805.315754579381_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,1,13) = (-17547.6432543683_pr*c1*gA2)/(fpi4*mpi2) - &
         (4066.12243819318_pr*c3*gA2)/(fpi4*mpi2) + &
         (563.160752122628_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(3,1,13) = (4641.24227101558_pr*c1*gA2)/(fpi4*mpi2) + &
         (1225.12479783611_pr*c3*gA2)/(fpi4*mpi2) - &
         (77.7313232421875_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,3,13) = (5176.02796875_pr*c1*gA2)/(fpi4*mpi2) + &
         (4535.422284375_pr*c3*gA2)/(fpi4*mpi2) + &
         (411.17551875_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,2,14) = (-19673.8226904297_pr*c1*gA2)/(fpi4*mpi2) - &
         (3790.5753035822_pr*c3*gA2)/(fpi4*mpi2) + &
         (1655.41310106027_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,2,14) = (-2415.02747628348_pr*c1*gA2)/(fpi4*mpi2) - &
         (3234.1744717979_pr*c3*gA2)/(fpi4*mpi2) - &
         (672.799001116071_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,1,15) = (20527.2675777519_pr*c1*gA2)/(fpi4*mpi2) + &
         (8144.9046506711_pr*c3*gA2)/(fpi4*mpi2) + &
         (783.247644786072_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,1,15) = (-17967.065175939_pr*c1*gA2)/(fpi4*mpi2) - &
         (7747.02659701377_pr*c3*gA2)/(fpi4*mpi2) - &
         (867.669720144653_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(3,1,15) = (6843.01542324125_pr*c1*gA2)/(fpi4*mpi2) + &
         (3014.94458858161_pr*c3*gA2)/(fpi4*mpi2) + &
         (216.528129577637_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,3,15) = (12800.9049609375_pr*c1*gA2)/(fpi4*mpi2) + &
         (3707.04542871094_pr*c3*gA2)/(fpi4*mpi2) - &
         (745.568033203125_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,2,16) = (-3092.52618530273_pr*c1*gA2)/(fpi4*mpi2) - &
         (1011.19649671026_pr*c3*gA2)/(fpi4*mpi2) - &
         (1042.07968925432_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,2,16) = (-8378.82125244141_pr*c1*gA2)/(fpi4*mpi2) - &
         (2987.23517760232_pr*c3*gA2)/(fpi4*mpi2) + &
         (628.253677891323_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,1,17) = (14135.5592743402_pr*c1*gA2)/(fpi4*mpi2) + &
         (5746.10194483582_pr*c3*gA2)/(fpi4*mpi2) - &
         (16.4675249980382_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,1,17) = (-20256.3592693034_pr*c1*gA2)/(fpi4*mpi2) - &
         (8245.07712593439_pr*c3*gA2)/(fpi4*mpi2) + &
         (142.309278507778_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(3,1,17) = (8660.61555253689_pr*c1*gA2)/(fpi4*mpi2) + &
         (3459.42204756108_pr*c3*gA2)/(fpi4*mpi2) - &
         (69.7415504673549_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,3,17) = (3484.14064453125_pr*c1*gA2)/(fpi4*mpi2) + &
         (1259.881028268_pr*c3*gA2)/(fpi4*mpi2) + &
         (342.195423235212_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,2,18) = (840.376118774414_pr*c1*gA2)/(fpi4*mpi2) + &
         (319.466587109805_pr*c3*gA2)/(fpi4*mpi2) + &
         (10.8880724487305_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,2,18) = (-3131.19715759277_pr*c1*gA2)/(fpi4*mpi2) - &
         (1166.81519764533_pr*c3*gA2)/(fpi4*mpi2) - &
         (87.1195371459961_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,1,19) = (7310.26069007057_pr*c1*gA2)/(fpi4*mpi2) + &
         (3000.94037657016_pr*c3*gA2)/(fpi4*mpi2) + &
         (0.0826331377846854_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,1,19) = (-13685.5267014613_pr*c1*gA2)/(fpi4*mpi2) - &
         (5572.85595346119_pr*c3*gA2)/(fpi4*mpi2) + &
         (2.70709313686916_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(3,1,19) = (6169.2187169198_pr*c1*gA2)/(fpi4*mpi2) + &
         (2481.7781746784_pr*c3*gA2)/(fpi4*mpi2) + &
         (3.76123477935791_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,2,20) = (-210.094029693604_pr*c1*gA2)/(fpi4*mpi2) - &
         (79.8666467774513_pr*c3*gA2)/(fpi4*mpi2) - &
         (2.72201811218262_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,1,21) = (2476.5165551998_pr*c1*gA2)/(fpi4*mpi2) + &
         (1029.07446289094_pr*c3*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,1,21) = (-5841.10206786854_pr*c1*gA2)/(fpi4*mpi2) - &
         (2407.9518963783_pr*c3*gA2)/(fpi4*mpi2) - &
         (0.0413165688923427_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(3,1,21) = (3018.45160491616_pr*c1*gA2)/(fpi4*mpi2) + &
         (1233.36858750781_pr*c3*gA2)/(fpi4*mpi2) + &
         (0.352447028160095_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,1,23) = (574.248747042225_pr*c1*gA2)/(fpi4*mpi2) + &
         (240.667492988211_pr*c3*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,1,23) = (-1764.09490454231_pr*c1*gA2)/(fpi4*mpi2) - &
         (734.821822475346_pr*c3*gA2)/(fpi4*mpi2)
    ctr1Dr1(3,1,23) = (1123.44919031994_pr*c1*gA2)/(fpi4*mpi2) + &
         (464.703617106343_pr*c3*gA2)/(fpi4*mpi2) + &
         (0.00516457111154284_pr*c4*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,1,25) = (86.2349934345314_pr*c1*gA2)/(fpi4*mpi2) + &
         (36.3168724059041_pr*c3*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,1,25) = (-368.716465873243_pr*c1*gA2)/(fpi4*mpi2) - &
         (154.692733933317_pr*c3*gA2)/(fpi4*mpi2)
    ctr1Dr1(3,1,25) = (305.511227790772_pr*c1*gA2)/(fpi4*mpi2) + &
         (127.502237247082_pr*c3*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,1,27) = (7.54914530309393_pr*c1*gA2)/(fpi4*mpi2) + &
         (3.18270638442274_pr*c3*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,1,27) = (-50.4717107938621_pr*c1*gA2)/(fpi4*mpi2) - &
         (21.25918294258_pr*c3*gA2)/(fpi4*mpi2)
    ctr1Dr1(3,1,27) = (58.0783904906964_pr*c1*gA2)/(fpi4*mpi2) + &
         (24.3861619452745_pr*c3*gA2)/(fpi4*mpi2)
    ctr1Dr1(1,1,29) = (0.292396839746317_pr*c1*gA2)/(fpi4*mpi2) + &
         (0.12293946719219_pr*c3*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,1,29) = (-4.06696949129328_pr*c1*gA2)/(fpi4*mpi2) - &
         (1.71429265940356_pr*c3*gA2)/(fpi4*mpi2)
    ctr1Dr1(3,1,29) = (7.30133981874389_pr*c1*gA2)/(fpi4*mpi2) + &
         (3.07572607707404_pr*c3*gA2)/(fpi4*mpi2)
    ctr1Dr1(2,1,31) = (-0.146198419873159_pr*c1*gA2)/(fpi4*mpi2) - &
         (0.0614697335960952_pr*c3*gA2)/(fpi4*mpi2)
    ctr1Dr1(3,1,31) = (0.54492079137995_pr*c1*gA2)/(fpi4*mpi2) + &
         (0.229654015824469_pr*c3*gA2)/(fpi4*mpi2)
    ctr1Dr1(3,1,33) = (0.0182748024841448_pr*c1*gA2)/(fpi4*mpi2) + &
         (0.00768371669951191_pr*c3*gA2)/(fpi4*mpi2)
    !
    !  Urho0tau0
    !
    ctr0t0=0.0_pr
    ctr0t0(1,1,3)=(1.36561930174741_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0479358100097182_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,1,5)=(10.5008244999755_pr*c1*gA2)/(fpi4*mpi2)- &
         (7.54491362519471_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.333010041454082_pr*c4*gA2)/(fpi4*mpi2)- &
         (0.0703125_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctr0t0(2,1,5)=(-0.24616366400062_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.117254464285714_pr*c4*gA2)/(fpi4*mpi2)+ &
         (0.03515625_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctr0t0(1,1,7)=(52.0673482366071_pr*c1*gA2)/(fpi4*mpi2)- &
         (25.8426066958322_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.144953338249362_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(2,1,7)=(-58.436493475888_pr*c1*gA2)/(fpi4*mpi2)+ &
         (27.9983571543377_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.114467235331633_pr*c4*gA2)/(fpi4*mpi2)+ &
         (0.017578125_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctr0t0(1,2,8)=(505.710701116071_pr*c1*gA2)/(fpi4*mpi2)- &
         (112.325267992606_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.9995709375_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(2,2,8)=(-54.9853125163899_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.28125_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,1,9)=(9474.83599128706_pr*c1*gA2)/(fpi4*mpi2)+ &
         (1465.87254191511_pr*c3*gA2)/(fpi4*mpi2)+ &
         (21.7653711399872_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(2,1,9)=(-3624.78698664406_pr*c1*gA2)/(fpi4*mpi2)- &
         (461.242897263176_pr*c3*gA2)/(fpi4*mpi2)- &
         (6.96088448660714_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(3,1,9)=(330.256406317608_pr*c1*gA2)/(fpi4*mpi2)+ &
         (97.9034765840117_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.8173828125_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,3,9)=(824.990574776786_pr*c3*gA2)/(fpi4*mpi2)+ &
         (3.79697544642857_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,2,10)=(-17975.7704772321_pr*c1*gA2)/(fpi4*mpi2)+ &
         (1262.86185981827_pr*c3*gA2)/(fpi4*mpi2)- &
         (2.65529641262755_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(2,2,10)=(1251.73928571429_pr*c1*gA2)/(fpi4*mpi2)- &
         (2090.34981043995_pr*c3*gA2)/(fpi4*mpi2)- &
         (6.21662946428571_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,1,11)=(34044.1058864892_pr*c1*gA2)/(fpi4*mpi2)+ &
         (30851.2885695409_pr*c3*gA2)/(fpi4*mpi2)+ &
         (650.736844237584_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(2,1,11)=(-12289.0849327548_pr*c1*gA2)/(fpi4*mpi2)- &
         (15854.7386609546_pr*c3*gA2)/(fpi4*mpi2)- &
         (255.112345354353_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(3,1,11)=(2814.07624373444_pr*c1*gA2)/(fpi4*mpi2)+ &
         (2910.54669625161_pr*c3*gA2)/(fpi4*mpi2)+ &
         (27.8796648297991_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,3,11)=(23028.03140625_pr*c1*gA2)/(fpi4*mpi2)+ &
         (3455.7492578125_pr*c3*gA2)/(fpi4*mpi2)- &
         (20.678203125_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,2,12)=(1388.89324888393_pr*c1*gA2)/(fpi4*mpi2)- &
         (36494.6219835576_pr*c3*gA2)/(fpi4*mpi2)- &
         (1798.79982589286_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(2,2,12)=(-30185.6310351562_pr*c1*gA2)/(fpi4*mpi2)+ &
         (232.70274650057_pr*c3*gA2)/(fpi4*mpi2)+ &
         (356.94017578125_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,1,13)=(169633.190522348_pr*c1*gA2)/(fpi4*mpi2)+ &
         (42011.3077626309_pr*c3*gA2)/(fpi4*mpi2)- &
         (2415.94726373814_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(2,1,13)=(-113865.746379575_pr*c1*gA2)/(fpi4*mpi2)- &
         (22150.6599403888_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1689.48225636789_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(3,1,13)=(29118.6234258982_pr*c1*gA2)/(fpi4*mpi2)+ &
         (7008.98206338386_pr*c3*gA2)/(fpi4*mpi2)- &
         (233.193969726563_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,3,13)=(25721.97103125_pr*c1*gA2)/(fpi4*mpi2)+ &
         (28797.99699375_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1233.52655625_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,2,14)=(-143206.840020996_pr*c1*gA2)/(fpi4*mpi2)- &
         (16594.1553967862_pr*c3*gA2)/(fpi4*mpi2)+ &
         (4966.2393031808_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(2,2,14)=(-4876.19730050223_pr*c1*gA2)/(fpi4*mpi2)- &
         (21921.4915168923_pr*c3*gA2)/(fpi4*mpi2)- &
         (2018.39700334821_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,1,15)=(111720.675330785_pr*c1*gA2)/(fpi4*mpi2)+ &
         (51670.9893465594_pr*c3*gA2)/(fpi4*mpi2)+ &
         (2349.74293435822_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(2,1,15)=(-95184.3030091052_pr*c1*gA2)/(fpi4*mpi2)- &
         (49572.8589648463_pr*c3*gA2)/(fpi4*mpi2)- &
         (2603.00916043396_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(3,1,15)=(38037.6118997991_pr*c1*gA2)/(fpi4*mpi2)+ &
         (18858.0289954179_pr*c3*gA2)/(fpi4*mpi2)+ &
         (649.58438873291_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,3,15)=(88149.5686757812_pr*c1*gA2)/(fpi4*mpi2)+ &
         (19508.5932099609_pr*c3*gA2)/(fpi4*mpi2)- &
         (2236.70409960937_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,2,16)=(-2933.63054626465_pr*c1*gA2)/(fpi4*mpi2)- &
         (9849.90830391705_pr*c3*gA2)/(fpi4*mpi2)- &
         (3126.23906776297_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(2,2,16)=(-59728.9286700439_pr*c1*gA2)/(fpi4*mpi2)- &
         (15659.2800568336_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1884.76103367397_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,1,17)=(85156.8253783411_pr*c1*gA2)/(fpi4*mpi2)+ &
         (34404.1775560654_pr*c3*gA2)/(fpi4*mpi2)- &
         (49.4025749941145_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(2,1,17)=(-123896.068978873_pr*c1*gA2)/(fpi4*mpi2)- &
         (48923.4556759794_pr*c3*gA2)/(fpi4*mpi2)+ &
         (426.927835523333_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(3,1,17)=(53102.5825769139_pr*c1*gA2)/(fpi4*mpi2)+ &
         (20496.0380452978_pr*c3*gA2)/(fpi4*mpi2)- &
         (209.224651402065_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,3,17)=(15678.6329003906_pr*c1*gA2)/(fpi4*mpi2)+ &
         (8819.16719787598_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1026.58626970564_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,2,18)=(4813.60998284912_pr*c1*gA2)/(fpi4*mpi2)+ &
         (1965.79485166884_pr*c3*gA2)/(fpi4*mpi2)+ &
         (32.6642173461914_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(2,2,18)=(-17352.030255249_pr*c1*gA2)/(fpi4*mpi2)- &
         (7339.58754773114_pr*c3*gA2)/(fpi4*mpi2)- &
         (261.358611437988_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,1,19)=(43859.5226517645_pr*c1*gA2)/(fpi4*mpi2)+ &
         (18006.1234631001_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.247899413354056_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(2,1,19)=(-82168.5302213599_pr*c1*gA2)/(fpi4*mpi2)- &
         (33425.8975324765_pr*c3*gA2)/(fpi4*mpi2)+ &
         (8.12127941060747_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(3,1,19)=(36955.7489794515_pr*c1*gA2)/(fpi4*mpi2)+ &
         (14906.1590782447_pr*c3*gA2)/(fpi4*mpi2)+ &
         (11.2837043380737_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(2,2,20)=(-1203.40249571228_pr*c1*gA2)/(fpi4*mpi2)- &
         (491.448712917211_pr*c3*gA2)/(fpi4*mpi2)- &
         (8.16605433654785_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,1,21)=(14859.0993311988_pr*c1*gA2)/(fpi4*mpi2)+ &
         (6174.44677734561_pr*c3*gA2)/(fpi4*mpi2)
    ctr0t0(2,1,21)=(-35045.5916628818_pr*c1*gA2)/(fpi4*mpi2)- &
         (14447.9519801093_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.123949706677028_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(3,1,21)=(18103.0852743763_pr*c1*gA2)/(fpi4*mpi2)+ &
         (7401.92911003352_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1.05734108448029_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,1,23)=(3445.49248225335_pr*c1*gA2)/(fpi4*mpi2)+ &
         (1444.00495792927_pr*c3*gA2)/(fpi4*mpi2)
    ctr0t0(2,1,23)=(-10584.5694272539_pr*c1*gA2)/(fpi4*mpi2)- &
         (4408.93093485208_pr*c3*gA2)/(fpi4*mpi2)
    ctr0t0(3,1,23)=(6740.56754887846_pr*c1*gA2)/(fpi4*mpi2)+ &
         (2788.251777868_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0154937133346285_pr*c4*gA2)/(fpi4*mpi2)
    ctr0t0(1,1,25)=(517.409960607188_pr*c1*gA2)/(fpi4*mpi2)+ &
         (217.901234435425_pr*c3*gA2)/(fpi4*mpi2)
    ctr0t0(2,1,25)=(-2212.29879523946_pr*c1*gA2)/(fpi4*mpi2)- &
         (928.156403599904_pr*c3*gA2)/(fpi4*mpi2)
    ctr0t0(3,1,25)=(1833.06736674463_pr*c1*gA2)/(fpi4*mpi2)+ &
         (765.013423482493_pr*c3*gA2)/(fpi4*mpi2)
    ctr0t0(1,1,27)=(45.2948718185636_pr*c1*gA2)/(fpi4*mpi2)+ &
         (19.0962383065364_pr*c3*gA2)/(fpi4*mpi2)
    ctr0t0(2,1,27)=(-302.830264763173_pr*c1*gA2)/(fpi4*mpi2)- &
         (127.55509765548_pr*c3*gA2)/(fpi4*mpi2)
    ctr0t0(3,1,27)=(348.470342944178_pr*c1*gA2)/(fpi4*mpi2)+ &
         (146.316971671647_pr*c3*gA2)/(fpi4*mpi2)
    ctr0t0(1,1,29)=(1.7543810384779_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.737636803153143_pr*c3*gA2)/(fpi4*mpi2)
    ctr0t0(2,1,29)=(-24.4018169477597_pr*c1*gA2)/(fpi4*mpi2)- &
         (10.2857559564214_pr*c3*gA2)/(fpi4*mpi2)
    ctr0t0(3,1,29)=(43.8080389124633_pr*c1*gA2)/(fpi4*mpi2)+ &
         (18.4543564624443_pr*c3*gA2)/(fpi4*mpi2)
    ctr0t0(2,1,31)=(-0.877190519238951_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.368818401576571_pr*c3*gA2)/(fpi4*mpi2)
    ctr0t0(3,1,31)=(3.2695247482797_pr*c1*gA2)/(fpi4*mpi2)+ &
         (1.37792409494681_pr*c3*gA2)/(fpi4*mpi2)
    ctr0t0(3,1,33)=(0.109648814904869_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.0461023001970714_pr*c3*gA2)/(fpi4*mpi2)
    !
    !  Urho1tau0
    !
    ctr1t0=0.0_pr
    ctr1t0(1,1,3)=(-0.0370474171525753_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0159786033365727_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(1,1,5)=(-0.191231023596939_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.0468828093112245_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.111003347151361_pr*c4*gA2)/(fpi4*mpi2)+ &
         (0.0234375_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctr1t0(2,1,5)=(0.00779017857142857_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0390848214285714_pr*c4*gA2)/(fpi4*mpi2)- &
         (0.01171875_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctr1t0(1,1,7)=(0.765432059151786_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.319270121173469_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0483177794164541_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(2,1,7)=(0.0512273596938776_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.00913767538265306_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0381557451105442_pr*c4*gA2)/(fpi4*mpi2)- &
         (0.005859375_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctr1t0(1,2,8)=(-1.41375736607143_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.846850379464286_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.3331903125_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(2,2,8)=(0.046875_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.09375_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(1,1,9)=(-22.8134614118304_pr*c1*gA2)/(fpi4*mpi2)- &
         (1.21977883211097_pr*c3*gA2)/(fpi4*mpi2)- &
         (7.25512371332908_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(2,1,9)=(7.09276506696429_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.199196428571429_pr*c3*gA2)/(fpi4*mpi2)+ &
         (2.32029482886905_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(3,1,9)=(-0.5009765625_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0966796875_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.2724609375_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(1,3,9)=(-1.89848772321429_pr*c3*gA2)/(fpi4*mpi2)- &
         (1.26565848214286_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(1,2,10)=(45.7770589285714_pr*c1*gA2)/(fpi4*mpi2)- &
         (37.9879278300383_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.885098804209184_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(2,2,10)=(-5.06216517857143_pr*c1*gA2)/(fpi4*mpi2)+ &
         (9.70597098214286_pr*c3*gA2)/(fpi4*mpi2)+ &
         (2.07220982142857_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(1,1,11)=(-440.874152071708_pr*c1*gA2)/(fpi4*mpi2)- &
         (597.813959665179_pr*c3*gA2)/(fpi4*mpi2)- &
         (216.912281412528_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(2,1,11)=(137.620562360491_pr*c1*gA2)/(fpi4*mpi2)+ &
         (242.982874720982_pr*c3*gA2)/(fpi4*mpi2)+ &
         (85.0374484514509_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(3,1,11)=(-12.2751813616071_pr*c1*gA2)/(fpi4*mpi2)- &
         (25.6883108956473_pr*c3*gA2)/(fpi4*mpi2)- &
         (9.29322160993304_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(1,3,11)=(-22.63640625_pr*c1*gA2)/(fpi4*mpi2)+ &
         (46.7609765625_pr*c3*gA2)/(fpi4*mpi2)+ &
         (6.892734375_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(1,2,12)=(1465.24899274554_pr*c1*gA2)/(fpi4*mpi2)+ &
         (1684.896046875_pr*c3*gA2)/(fpi4*mpi2)+ &
         (599.599941964286_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(2,2,12)=(-235.90740234375_pr*c1*gA2)/(fpi4*mpi2)- &
         (351.521103515625_pr*c3*gA2)/(fpi4*mpi2)- &
         (118.98005859375_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(1,1,13)=(2941.54591061489_pr*c1*gA2)/(fpi4*mpi2)+ &
         (2355.33459348319_pr*c3*gA2)/(fpi4*mpi2)+ &
         (805.315754579381_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(2,1,13)=(-1813.61831472866_pr*c1*gA2)/(fpi4*mpi2)- &
         (1619.87784794992_pr*c3*gA2)/(fpi4*mpi2)- &
         (563.160752122628_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(3,1,13)=(221.161743164063_pr*c1*gA2)/(fpi4*mpi2)+ &
         (225.247192382813_pr*c3*gA2)/(fpi4*mpi2)+ &
         (77.7313232421875_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(1,3,13)=(-1171.69059375_pr*c1*gA2)/(fpi4*mpi2)- &
         (1176.1601625_pr*c3*gA2)/(fpi4*mpi2)- &
         (411.17551875_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(1,2,14)=(-6502.52069384766_pr*c1*gA2)/(fpi4*mpi2)- &
         (5016.80990126953_pr*c3*gA2)/(fpi4*mpi2)- &
         (1655.41310106027_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(2,2,14)=(2369.16757114955_pr*c1*gA2)/(fpi4*mpi2)+ &
         (1991.91343610491_pr*c3*gA2)/(fpi4*mpi2)+ &
         (672.799001116071_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(1,1,15)=(-3755.98198127747_pr*c1*gA2)/(fpi4*mpi2)- &
         (2756.34979092146_pr*c3*gA2)/(fpi4*mpi2)- &
         (783.247644786072_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(2,1,15)=(3861.25969552264_pr*c1*gA2)/(fpi4*mpi2)+ &
         (2865.95377119141_pr*c3*gA2)/(fpi4*mpi2)+ &
         (867.669720144653_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(3,1,15)=(-880.419287109375_pr*c1*gA2)/(fpi4*mpi2)- &
         (683.312782287598_pr*c3*gA2)/(fpi4*mpi2)- &
         (216.528129577637_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(1,3,15)=(3099.20776171875_pr*c1*gA2)/(fpi4*mpi2)+ &
         (2319.50358105469_pr*c3*gA2)/(fpi4*mpi2)+ &
         (745.568033203125_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(1,2,16)=(5164.53977966309_pr*c1*gA2)/(fpi4*mpi2)+ &
         (3750.44797599923_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1042.07968925432_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(2,2,16)=(-2919.87245544434_pr*c1*gA2)/(fpi4*mpi2)- &
         (2116.93618456159_pr*c3*gA2)/(fpi4*mpi2)- &
         (628.253677891323_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(1,1,17)=(113.823727294922_pr*c1*gA2)/(fpi4*mpi2)+ &
         (71.9582676125663_pr*c3*gA2)/(fpi4*mpi2)+ &
         (16.4675249980382_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(2,1,17)=(-793.151098434448_pr*c1*gA2)/(fpi4*mpi2)- &
         (552.996386267526_pr*c3*gA2)/(fpi4*mpi2)- &
         (142.309278507778_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(3,1,17)=(364.767998526437_pr*c1*gA2)/(fpi4*mpi2)+ &
         (251.980659746443_pr*c3*gA2)/(fpi4*mpi2)+ &
         (69.7415504673549_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(1,3,17)=(-1742.07032226563_pr*c1*gA2)/(fpi4*mpi2)- &
         (1259.881028268_pr*c3*gA2)/(fpi4*mpi2)- &
         (342.195423235212_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(1,2,18)=(-76.2155765991211_pr*c1*gA2)/(fpi4*mpi2)- &
         (48.9953290100098_pr*c3*gA2)/(fpi4*mpi2)- &
         (10.8880724487305_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(2,2,18)=(489.043165649414_pr*c1*gA2)/(fpi4*mpi2)+ &
         (346.766698773193_pr*c3*gA2)/(fpi4*mpi2)+ &
         (87.1195371459961_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(1,1,19)=(-0.68049621963501_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.481203679112026_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0826331377846854_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(2,1,19)=(-18.1235791282654_pr*c1*gA2)/(fpi4*mpi2)- &
         (11.0002656221662_pr*c3*gA2)/(fpi4*mpi2)- &
         (2.70709313686916_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(3,1,19)=(-21.7049503326416_pr*c1*gA2)/(fpi4*mpi2)- &
         (16.8184317398071_pr*c3*gA2)/(fpi4*mpi2)- &
         (3.76123477935791_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(2,2,20)=(19.0538941497803_pr*c1*gA2)/(fpi4*mpi2)+ &
         (12.2488322525024_pr*c3*gA2)/(fpi4*mpi2)+ &
         (2.72201811218262_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(2,1,21)=(0.340248109817505_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.240601839556013_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0413165688923427_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(3,1,21)=(-2.58308817386627_pr*c1*gA2)/(fpi4*mpi2)- &
         (1.74732532024384_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.352447028160095_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t0(3,1,23)=(-0.0425310137271881_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0300752299445016_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00516457111154284_pr*c4*gA2)/(fpi4*mpi2)
    !
    !  Urho1tau1
    !
    ctr1t1=0.0_pr
    ctr1t1(1,1,3)=(-0.450116223436374_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0319572066731454_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,1,5)=(-3.5385371069221_pr*c1*gA2)/(fpi4*mpi2)+ &
         (2.50259509338973_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.222006694302721_pr*c4*gA2)/(fpi4*mpi2)+ &
         (0.046875_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctr1t1(2,1,5)=(0.0898447332383019_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0781696428571429_pr*c4*gA2)/(fpi4*mpi2)- &
         (0.0234375_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctr1t1(1,1,7)=(-17.2434756863839_pr*c1*gA2)/(fpi4*mpi2)+ &
         (8.60257873066855_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0966355588329082_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(2,1,7)=(19.5159960183232_pr*c1*gA2)/(fpi4*mpi2)- &
         (9.34811780250883_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0763114902210884_pr*c4*gA2)/(fpi4*mpi2)- &
         (0.01171875_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctr1t1(1,2,8)=(-168.933991071429_pr*c1*gA2)/(fpi4*mpi2)+ &
         (37.320722537714_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.666380625_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(2,2,8)=(18.3753125054633_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.1875_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,1,9)=(-3146.36888298668_pr*c1*gA2)/(fpi4*mpi2)- &
         (492.901337694334_pr*c3*gA2)/(fpi4*mpi2)- &
         (14.5102474266582_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(2,1,9)=(1205.47618769832_pr*c1*gA2)/(fpi4*mpi2)+ &
         (155.024983611535_pr*c3*gA2)/(fpi4*mpi2)+ &
         (4.6405896577381_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(3,1,9)=(-109.918476585036_pr*c1*gA2)/(fpi4*mpi2)- &
         (32.7311718821706_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.544921875_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,3,9)=(-274.364029017857_pr*c3*gA2)/(fpi4*mpi2)- &
         (2.53131696428571_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,2,10)=(5901.79117633929_pr*c1*gA2)/(fpi4*mpi2)- &
         (385.542873162743_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1.77019760841837_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(2,2,10)=(-404.58984375_pr*c1*gA2)/(fpi4*mpi2)+ &
         (690.035446485936_pr*c3*gA2)/(fpi4*mpi2)+ &
         (4.14441964285714_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,1,11)=(-12146.9163811618_pr*c1*gA2)/(fpi4*mpi2)- &
         (9992.47844580492_pr*c3*gA2)/(fpi4*mpi2)- &
         (433.824562825056_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(2,1,11)=(4424.78103473708_pr*c1*gA2)/(fpi4*mpi2)+ &
         (5157.68192548934_pr*c3*gA2)/(fpi4*mpi2)+ &
         (170.074896902902_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(3,1,11)=(-968.792783439754_pr*c1*gA2)/(fpi4*mpi2)- &
         (956.662918816572_pr*c3*gA2)/(fpi4*mpi2)- &
         (18.5864432198661_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,3,11)=(-7552.959375_pr*c1*gA2)/(fpi4*mpi2)- &
         (1191.78466145833_pr*c3*gA2)/(fpi4*mpi2)+ &
         (13.78546875_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,2,12)=(2094.44863895089_pr*c1*gA2)/(fpi4*mpi2)+ &
         (11355.6974476442_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1199.19988392857_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(2,2,12)=(9461.258671875_pr*c1*gA2)/(fpi4*mpi2)+ &
         (113.574869005018_pr*c3*gA2)/(fpi4*mpi2)- &
         (237.9601171875_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,1,13)=(-52317.8224267498_pr*c1*gA2)/(fpi4*mpi2)- &
         (15039.7248405472_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1610.63150915876_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(2,1,13)=(35095.2865087366_pr*c1*gA2)/(fpi4*mpi2)+ &
         (8132.24487638635_pr*c3*gA2)/(fpi4*mpi2)- &
         (1126.32150424526_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(3,1,13)=(-9282.48454203116_pr*c1*gA2)/(fpi4*mpi2)- &
         (2450.24959567222_pr*c3*gA2)/(fpi4*mpi2)+ &
         (155.462646484375_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,3,13)=(-10352.0559375_pr*c1*gA2)/(fpi4*mpi2)- &
         (9070.84456875_pr*c3*gA2)/(fpi4*mpi2)- &
         (822.3510375_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,2,14)=(39347.6453808594_pr*c1*gA2)/(fpi4*mpi2)+ &
         (7581.1506071644_pr*c3*gA2)/(fpi4*mpi2)- &
         (3310.82620212054_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(2,2,14)=(4830.05495256696_pr*c1*gA2)/(fpi4*mpi2)+ &
         (6468.34894359579_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1345.59800223214_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,1,15)=(-41054.5351555039_pr*c1*gA2)/(fpi4*mpi2)- &
         (16289.8093013422_pr*c3*gA2)/(fpi4*mpi2)- &
         (1566.49528957214_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(2,1,15)=(35934.130351878_pr*c1*gA2)/(fpi4*mpi2)+ &
         (15494.0531940275_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1735.33944028931_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(3,1,15)=(-13686.0308464825_pr*c1*gA2)/(fpi4*mpi2)- &
         (6029.88917716322_pr*c3*gA2)/(fpi4*mpi2)- &
         (433.056259155273_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,3,15)=(-25601.809921875_pr*c1*gA2)/(fpi4*mpi2)- &
         (7414.09085742188_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1491.13606640625_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,2,16)=(6185.05237060547_pr*c1*gA2)/(fpi4*mpi2)+ &
         (2022.39299342052_pr*c3*gA2)/(fpi4*mpi2)+ &
         (2084.15937850865_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(2,2,16)=(16757.6425048828_pr*c1*gA2)/(fpi4*mpi2)+ &
         (5974.47035520464_pr*c3*gA2)/(fpi4*mpi2)- &
         (1256.50735578265_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,1,17)=(-28271.1185486804_pr*c1*gA2)/(fpi4*mpi2)- &
         (11492.2038896716_pr*c3*gA2)/(fpi4*mpi2)+ &
         (32.9350499960763_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(2,1,17)=(40512.7185386068_pr*c1*gA2)/(fpi4*mpi2)+ &
         (16490.1542518688_pr*c3*gA2)/(fpi4*mpi2)- &
         (284.618557015555_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(3,1,17)=(-17321.2311050738_pr*c1*gA2)/(fpi4*mpi2)- &
         (6918.84409512217_pr*c3*gA2)/(fpi4*mpi2)+ &
         (139.48310093471_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,3,17)=(-6968.2812890625_pr*c1*gA2)/(fpi4*mpi2)- &
         (2519.76205653599_pr*c3*gA2)/(fpi4*mpi2)- &
         (684.390846470424_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,2,18)=(-1680.75223754883_pr*c1*gA2)/(fpi4*mpi2)- &
         (638.933174219611_pr*c3*gA2)/(fpi4*mpi2)- &
         (21.7761448974609_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(2,2,18)=(6262.39431518555_pr*c1*gA2)/(fpi4*mpi2)+ &
         (2333.63039529067_pr*c3*gA2)/(fpi4*mpi2)+ &
         (174.239074291992_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,1,19)=(-14620.5213801411_pr*c1*gA2)/(fpi4*mpi2)- &
         (6001.88075314032_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.165266275569371_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(2,1,19)=(27371.0534029225_pr*c1*gA2)/(fpi4*mpi2)+ &
         (11145.7119069224_pr*c3*gA2)/(fpi4*mpi2)- &
         (5.41418627373832_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(3,1,19)=(-12338.4374338396_pr*c1*gA2)/(fpi4*mpi2)- &
         (4963.55634935681_pr*c3*gA2)/(fpi4*mpi2)- &
         (7.52246955871582_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(2,2,20)=(420.188059387207_pr*c1*gA2)/(fpi4*mpi2)+ &
         (159.733293554903_pr*c3*gA2)/(fpi4*mpi2)+ &
         (5.44403622436523_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,1,21)=(-4953.03311039961_pr*c1*gA2)/(fpi4*mpi2)- &
         (2058.14892578187_pr*c3*gA2)/(fpi4*mpi2)
    ctr1t1(2,1,21)=(11682.2041357371_pr*c1*gA2)/(fpi4*mpi2)+ &
         (4815.90379275659_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0826331377846854_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(3,1,21)=(-6036.90320983232_pr*c1*gA2)/(fpi4*mpi2)- &
         (2466.73717501561_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.70489405632019_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,1,23)=(-1148.49749408445_pr*c1*gA2)/(fpi4*mpi2)- &
         (481.334985976422_pr*c3*gA2)/(fpi4*mpi2)
    ctr1t1(2,1,23)=(3528.18980908462_pr*c1*gA2)/(fpi4*mpi2)+ &
         (1469.64364495069_pr*c3*gA2)/(fpi4*mpi2)
    ctr1t1(3,1,23)=(-2246.89838063988_pr*c1*gA2)/(fpi4*mpi2)- &
         (929.407234212685_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0103291422230857_pr*c4*gA2)/(fpi4*mpi2)
    ctr1t1(1,1,25)=(-172.469986869063_pr*c1*gA2)/(fpi4*mpi2)- &
         (72.6337448118083_pr*c3*gA2)/(fpi4*mpi2)
    ctr1t1(2,1,25)=(737.432931746486_pr*c1*gA2)/(fpi4*mpi2)+ &
         (309.385467866635_pr*c3*gA2)/(fpi4*mpi2)
    ctr1t1(3,1,25)=(-611.022455581544_pr*c1*gA2)/(fpi4*mpi2)- &
         (255.004474494164_pr*c3*gA2)/(fpi4*mpi2)
    ctr1t1(1,1,27)=(-15.0982906061879_pr*c1*gA2)/(fpi4*mpi2)- &
         (6.36541276884548_pr*c3*gA2)/(fpi4*mpi2)
    ctr1t1(2,1,27)=(100.943421587724_pr*c1*gA2)/(fpi4*mpi2)+ &
         (42.51836588516_pr*c3*gA2)/(fpi4*mpi2)
    ctr1t1(3,1,27)=(-116.156780981393_pr*c1*gA2)/(fpi4*mpi2)- &
         (48.7723238905489_pr*c3*gA2)/(fpi4*mpi2)
    ctr1t1(1,1,29)=(-0.584793679492634_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.245878934384381_pr*c3*gA2)/(fpi4*mpi2)
    ctr1t1(2,1,29)=(8.13393898258656_pr*c1*gA2)/(fpi4*mpi2)+ &
         (3.42858531880712_pr*c3*gA2)/(fpi4*mpi2)
    ctr1t1(3,1,29)=(-14.6026796374878_pr*c1*gA2)/(fpi4*mpi2)- &
         (6.15145215414809_pr*c3*gA2)/(fpi4*mpi2)
    ctr1t1(2,1,31)=(0.292396839746317_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.12293946719219_pr*c3*gA2)/(fpi4*mpi2)
    ctr1t1(3,1,31)=(-1.0898415827599_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.459308031648938_pr*c3*gA2)/(fpi4*mpi2)
    ctr1t1(3,1,33)=(-0.0365496049682896_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0153674333990238_pr*c3*gA2)/(fpi4*mpi2)
    !
    !  UJ0Drho0
    !
    ctJ0dr0=0.0_pr
    !
    !  Urho0DJ0
    !
    ctr0DJ0=0.0_pr
    ctr0DJ0(1,1,3)=(-0.244419766459837_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0234375000151846_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,1,5)=(-1.03158597015212_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.780470479470306_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0351562492028061_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,1,7)=(-21.0022652575419_pr*c1*gA2)/(fpi4*mpi2)+ &
         (2.27019651550668_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00439427564971301_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(2,1,7)=(6.64452901902785_pr*c1*gA2)/(fpi4*mpi2)- &
         (1.38867403439195_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0175780851403061_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,2,8)=(-0.000291666666666667_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.0000478125_pr*c3*gA2)/(fpi4*mpi2)- &
         (4.01785714285714e-7_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,1,9)=(-307.281664803321_pr*c1*gA2)/(fpi4*mpi2)- &
         (17.4248922644585_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0131892671795281_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(2,1,9)=(167.072976616587_pr*c1*gA2)/(fpi4*mpi2)+ &
         (4.91318494493787_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00878976004464286_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(3,1,9)=(-20.6458593811462_pr*c1*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,2,10)=(0.07489375_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.000520145089285714_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0000100207270408163_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(2,2,10)=(-0.009375_pr*c1*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,1,11)=(-1417.24568001336_pr*c1*gA2)/(fpi4*mpi2)- &
         (148.086566637738_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00426528076171875_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(2,1,11)=(1006.30895614752_pr*c1*gA2)/(fpi4*mpi2)+ &
         (93.3959309491168_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0102003609793527_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(3,1,11)=(-165.352324283285_pr*c1*gA2)/(fpi4*mpi2)- &
         (12.9178818130112_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00332205636160714_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,3,11)=(-0.0875_pr*c1*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,2,12)=(-0.663125_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.09263671875_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00322935267857143_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(2,2,12)=(0.22640625_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.01634765625_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00052734375_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,1,13)=(-3324.31750521588_pr*c1*gA2)/(fpi4*mpi2)- &
         (384.460747217616_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00614578334263393_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(2,1,13)=(2946.4643595763_pr*c1*gA2)/(fpi4*mpi2)+ &
         (332.165470161666_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00216459803989955_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(3,1,13)=(-602.212375765244_pr*c1*gA2)/(fpi4*mpi2)- &
         (64.5897786600606_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.001190185546875_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,3,13)=(0.37125_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0631546875_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0026578125_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,2,14)=(1.72600911458333_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.135301904296875_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0143205845424107_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(2,2,14)=(-0.62890625_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.0793212890625_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00530482700892857_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,1,15)=(-4612.5195280616_pr*c1*gA2)/(fpi4*mpi2)- &
         (513.154765943116_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00963807367597307_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(2,1,15)=(5043.02746283907_pr*c1*gA2)/(fpi4*mpi2)+ &
         (573.956470083127_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00939304275512695_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(3,1,15)=(-1259.59323184652_pr*c1*gA2)/(fpi4*mpi2)- &
         (143.759310262187_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00231170654296875_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,3,15)=(-0.71465625_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.0487265625_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0069609375_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,2,16)=(-0.13190234375_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.084516064453125_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0133438685825893_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(2,2,16)=(0.345556640625_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0366448974609375_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0070806884765625_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,1,17)=(-4081.8063003269_pr*c1*gA2)/(fpi4*mpi2)- &
         (405.584681811802_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.000288566207885742_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(2,1,17)=(5511.17933250622_pr*c1*gA2)/(fpi4*mpi2)+ &
         (582.019749761264_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00207645034790039_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(3,1,17)=(-1676.42040558742_pr*c1*gA2)/(fpi4*mpi2)- &
         (184.111758668806_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.000897445678710938_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,3,17)=(-0.031904296875_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00455775669642857_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,2,18)=(-0.0106455078125_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.00226593017578125_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00019940185546875_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(2,2,18)=(0.03829833984375_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.0119082458496094_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00132149047851563_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,1,19)=(-2389.93945013977_pr*c1*gA2)/(fpi4*mpi2)- &
         (200.758956740934_pr*c3*gA2)/(fpi4*mpi2)+ &
         (2.18095779418945e-6_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(2,1,19)=(4026.89040316843_pr*c1*gA2)/(fpi4*mpi2)+ &
         (373.930398299718_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0000424012184143066_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(3,1,19)=(-1499.34103770493_pr*c1*gA2)/(fpi4*mpi2)- &
         (149.680535490279_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0000683426856994629_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(2,2,20)=(0.002661376953125_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.000566482543945312_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0000498504638671875_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,1,21)=(-937.343180466548_pr*c1*gA2)/(fpi4*mpi2)- &
         (62.9319292613857_pr*c3*gA2)/(fpi4*mpi2)
    ctr0DJ0(2,1,21)=(2012.17746809803_pr*c1*gA2)/(fpi4*mpi2)+ &
         (156.764429796522_pr*c3*gA2)/(fpi4*mpi2)- &
         (1.09047889709473e-6_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(3,1,21)=(-928.899230533649_pr*c1*gA2)/(fpi4*mpi2)- &
         (80.7120141836349_pr*c3*gA2)/(fpi4*mpi2)+ &
         (7.43508338928223e-6_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,1,23)=(-243.60536380909_pr*c1*gA2)/(fpi4*mpi2)- &
         (12.1417483739015_pr*c3*gA2)/(fpi4*mpi2)
    ctr0DJ0(2,1,23)=(690.008776700697_pr*c1*gA2)/(fpi4*mpi2)+ &
         (42.8121424298444_pr*c3*gA2)/(fpi4*mpi2)
    ctr0DJ0(3,1,23)=(-403.794812461118_pr*c1*gA2)/(fpi4*mpi2)- &
         (29.2881400946876_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1.36309862136841e-7_pr*c4*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,1,25)=(-40.2494438142034_pr*c1*gA2)/(fpi4*mpi2)- &
         (1.31603730589479_pr*c3*gA2)/(fpi4*mpi2)
    ctr0DJ0(2,1,25)=(159.713027711956_pr*c1*gA2)/(fpi4*mpi2)+ &
         (7.34601767825556_pr*c3*gA2)/(fpi4*mpi2)
    ctr0DJ0(3,1,25)=(-122.828777785615_pr*c1*gA2)/(fpi4*mpi2)- &
         (7.0783524478032_pr*c3*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,1,27)=(-3.82761492020625_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0613407218849266_pr*c3*gA2)/(fpi4*mpi2)
    ctr0DJ0(2,1,27)=(23.8460141906355_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.719359374832321_pr*c3*gA2)/(fpi4*mpi2)
    ctr0DJ0(3,1,27)=(-25.6066641016362_pr*c1*gA2)/(fpi4*mpi2)- &
         (1.09298032666628_pr*c3*gA2)/(fpi4*mpi2)
    ctr0DJ0(1,1,29)=(-0.159483955008594_pr*c1*gA2)/(fpi4*mpi2)
    ctr0DJ0(2,1,29)=(2.07329141511172_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.0306703609424633_pr*c3*gA2)/(fpi4*mpi2)
    ctr0DJ0(3,1,29)=(-3.48578429802332_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0975875120896559_pr*c3*gA2)/(fpi4*mpi2)
    ctr0DJ0(2,1,31)=(0.0797419775042969_pr*c1*gA2)/(fpi4*mpi2)
    ctr0DJ0(3,1,31)=(-0.279096921265039_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.00383379511780791_pr*c3*gA2)/(fpi4*mpi2)
    ctr0DJ0(3,1,33)=(-0.00996774718803712_pr*c1*gA2)/(fpi4*mpi2)
    !
    !  UJ1Drho1
    !
    ctJ1dr1=0.0_pr
    !
    !  Urho1DJ1
    !
    ctr1dJ1=0.0_pr
    ctr1dJ1(1,1,3)=(0.0736607554815508_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0156250000101231_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,1,5)=(0.308705740050707_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.248438076755833_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0234374994685374_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,1,7)=(6.96559883584731_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.755267413285654_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00292951709980867_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(2,1,7)=(-2.17968675634262_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.45703198308388_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0117187234268707_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,2,8)=(0.0000972222222222222_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0000158035714285714_pr*c3*gA2)/(fpi4*mpi2)+ &
         (2.67857142857143e-7_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,1,9)=(102.418432538607_pr*c1*gA2)/(fpi4*mpi2)+ &
         (5.80390099909298_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00879284478635204_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(2,1,9)=(-55.6646250180288_pr*c1*gA2)/(fpi4*mpi2)- &
         (1.63479839496441_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0058598400297619_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(3,1,9)=(6.87316406454873_pr*c1*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,2,10)=(-0.0249645833333333_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.000170041454081633_pr*c3*gA2)/(fpi4*mpi2)+ &
         (6.68048469387755e-6_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(2,2,10)=(0.003125_pr*c1*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,1,11)=(472.415226671119_pr*c1*gA2)/(fpi4*mpi2)+ &
         (49.3607671189922_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0028435205078125_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(2,1,11)=(-335.431924184589_pr*c1*gA2)/(fpi4*mpi2)- &
         (31.1285768627125_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00680024065290179_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(3,1,11)=(55.1130468965117_pr*c1*gA2)/(fpi4*mpi2)+ &
         (4.30485325221653_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00221470424107143_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,3,11)=(0.0291666666666667_pr*c1*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,2,12)=(0.221041666666667_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0298024553571429_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00215290178571429_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(2,2,12)=(-0.07546875_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.0052734375_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0003515625_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,1,13)=(1108.10583507196_pr*c1*gA2)/(fpi4*mpi2)+ &
         (128.155631000319_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00409718889508929_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(2,1,13)=(-982.154786525432_pr*c1*gA2)/(fpi4*mpi2)- &
         (110.722544919902_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00144306535993304_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(3,1,13)=(200.736909272009_pr*c1*gA2)/(fpi4*mpi2)+ &
         (21.5295294915046_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00079345703125_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,3,13)=(-0.12375_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.020165625_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.001771875_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,2,14)=(-0.575336371527778_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.0403271065848214_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00954705636160714_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(2,2,14)=(0.209635416666667_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0246721540178571_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00353655133928571_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,1,15)=(1537.50650935387_pr*c1*gA2)/(fpi4*mpi2)+ &
         (171.04837595648_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00642538245064872_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(2,1,15)=(-1681.00915427969_pr*c1*gA2)/(fpi4*mpi2)- &
         (191.31569234679_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00626202850341797_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(3,1,15)=(419.864410615508_pr*c1*gA2)/(fpi4*mpi2)+ &
         (47.9189995185479_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0015411376953125_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,3,15)=(0.23821875_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.013921875_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.004640625_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,2,16)=(0.0439674479166667_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0237240652901786_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00889591238839286_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(2,2,16)=(-0.115185546875_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.009854736328125_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.004720458984375_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,1,17)=(1360.60210010897_pr*c1*gA2)/(fpi4*mpi2)+ &
         (135.194990126003_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.000192377471923828_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(2,1,17)=(-1837.05977750207_pr*c1*gA2)/(fpi4*mpi2)- &
         (194.007275403871_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00138430023193359_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(3,1,17)=(558.806801862472_pr*c1*gA2)/(fpi4*mpi2)+ &
         (61.370885371495_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.000598297119140625_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,3,17)=(0.00911551339285714_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00303850446428571_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,2,18)=(0.00354850260416667_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.0006888427734375_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0001329345703125_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(2,2,18)=(-0.01276611328125_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.00352891845703125_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00088099365234375_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,1,19)=(796.646483379922_pr*c1*gA2)/(fpi4*mpi2)+ &
         (66.919651519992_pr*c3*gA2)/(fpi4*mpi2)- &
         (1.45397186279297e-6_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(2,1,19)=(-1342.29680105614_pr*c1*gA2)/(fpi4*mpi2)- &
         (124.643480233646_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0000282674789428711_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(3,1,19)=(499.780345901644_pr*c1*gA2)/(fpi4*mpi2)+ &
         (49.8934890491979_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0000455617904663086_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(2,2,20)=(-0.000887125651041667_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.000172210693359375_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.000033233642578125_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,1,21)=(312.447726822183_pr*c1*gA2)/(fpi4*mpi2)+ &
         (20.9773097537952_pr*c3*gA2)/(fpi4*mpi2)
    ctr1dJ1(2,1,21)=(-670.725822699343_pr*c1*gA2)/(fpi4*mpi2)- &
         (52.2548095686812_pr*c3*gA2)/(fpi4*mpi2)+ &
         (7.26985931396484e-7_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(3,1,21)=(309.63307684455_pr*c1*gA2)/(fpi4*mpi2)+ &
         (26.9040022495172_pr*c3*gA2)/(fpi4*mpi2)- &
         (4.95672225952148e-6_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,1,23)=(81.2017879363635_pr*c1*gA2)/(fpi4*mpi2)+ &
         (4.04724945796715_pr*c3*gA2)/(fpi4*mpi2)
    ctr1dJ1(2,1,23)=(-230.002925566899_pr*c1*gA2)/(fpi4*mpi2)- &
         (14.2707141432815_pr*c3*gA2)/(fpi4*mpi2)
    ctr1dJ1(3,1,23)=(134.598270820373_pr*c1*gA2)/(fpi4*mpi2)+ &
         (9.76271331945924_pr*c3*gA2)/(fpi4*mpi2)- &
         (9.08732414245605e-8_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,1,25)=(13.4164812714011_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.43867910196493_pr*c3*gA2)/(fpi4*mpi2)
    ctr1dJ1(2,1,25)=(-53.2376759039853_pr*c1*gA2)/(fpi4*mpi2)- &
         (2.44867255941852_pr*c3*gA2)/(fpi4*mpi2)
    ctr1dJ1(3,1,25)=(40.9429259285383_pr*c1*gA2)/(fpi4*mpi2)+ &
         (2.3594508159344_pr*c3*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,1,27)=(1.27587164006875_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.0204469072949755_pr*c3*gA2)/(fpi4*mpi2)
    ctr1dJ1(2,1,27)=(-7.94867139687851_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.23978645827744_pr*c3*gA2)/(fpi4*mpi2)
    ctr1dJ1(3,1,27)=(8.53555470054539_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.364326775555427_pr*c3*gA2)/(fpi4*mpi2)
    ctr1dJ1(1,1,29)=(0.053161318336198_pr*c1*gA2)/(fpi4*mpi2)
    ctr1dJ1(2,1,29)=(-0.691097138370574_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0102234536474878_pr*c3*gA2)/(fpi4*mpi2)
    ctr1dJ1(3,1,29)=(1.16192809934111_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.032529170696552_pr*c3*gA2)/(fpi4*mpi2)
    ctr1dJ1(2,1,31)=(-0.026580659168099_pr*c1*gA2)/(fpi4*mpi2)
    ctr1dJ1(3,1,31)=(0.0930323070883464_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.00127793170593597_pr*c3*gA2)/(fpi4*mpi2)
    ctr1dJ1(3,1,33)=(0.00332258239601237_pr*c1*gA2)/(fpi4*mpi2)
    !
    !  UJ0Drho1
    !
    ctJ0dr1=0.0_pr
    !
    !  Urho1DJ0
    !
    ctr1dJ0=0.0_pr
    ctr1dJ0(1,1,3)=(0.0234375000151846_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00781250000506155_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(1,1,5)=(0.10546875_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0351562492028061_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0117187497342687_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(1,1,7)=(0.10546875_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.00439427564971301_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00146475854990434_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(2,1,7)=(-0.10546875_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.0175780851403061_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00585936171343537_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(1,2,8)=(-4.01785714285714e-7_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1.33928571428571e-7_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(1,1,9)=(0.0263671875_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.0131892671795281_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00439642239317602_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(2,1,9)=(-0.0791015625_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.00878976004464286_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00292992001488095_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(3,1,9)=(0.0263671875_pr*c1*gA2)/(fpi4*mpi2)
    ctr1dJ0(1,2,10)=(-0.0000100207270408163_pr*c3*gA2)/(fpi4*mpi2)+ &
         (3.34024234693878e-6_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(1,1,11)=(0.00426528076171875_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00142176025390625_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(2,1,11)=(-0.01318359375_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0102003609793527_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00340012032645089_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(3,1,11)=(0.01318359375_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.00332205636160714_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00110735212053571_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(1,2,12)=(-0.00322935267857143_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00107645089285714_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(2,2,12)=(0.00052734375_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00017578125_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(1,1,13)=(-0.00614578334263393_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00204859444754464_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(2,1,13)=(0.00216459803989955_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.000721532679966518_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(3,1,13)=(0.00164794921875_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.001190185546875_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.000396728515625_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(1,3,13)=(0.0026578125_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0008859375_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(1,2,14)=(0.0143205845424107_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00477352818080357_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(2,2,14)=(-0.00530482700892857_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00176827566964286_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(1,1,15)=(0.00963807367597307_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00321269122532436_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(2,1,15)=(-0.00939304275512695_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00313101425170898_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(3,1,15)=(0.00231170654296875_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00077056884765625_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(1,3,15)=(-0.0069609375_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0023203125_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(1,2,16)=(-0.0133438685825893_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00444795619419643_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(2,2,16)=(0.0070806884765625_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0023602294921875_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(1,1,17)=(-0.000288566207885742_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0000961887359619141_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(2,1,17)=(0.00207645034790039_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.000692150115966797_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(3,1,17)=(-0.000897445678710938_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.000299148559570312_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(1,3,17)=(0.00455775669642857_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00151925223214286_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(1,2,18)=(0.00019940185546875_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00006646728515625_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(2,2,18)=(-0.00132149047851563_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.000440496826171875_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(1,1,19)=(2.18095779418945e-6_pr*c3*gA2)/(fpi4*mpi2)- &
         (7.26985931396484e-7_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(2,1,19)=(0.0000424012184143066_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0000141337394714355_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(3,1,19)=(0.0000683426856994629_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0000227808952331543_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(2,2,20)=(-0.0000498504638671875_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0000166168212890625_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(2,1,21)=(-1.09047889709473e-6_pr*c3*gA2)/(fpi4*mpi2)+ &
         (3.63492965698242e-7_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(3,1,21)=(7.43508338928223e-6_pr*c3*gA2)/(fpi4*mpi2)- &
         (2.47836112976074e-6_pr*c4*gA2)/(fpi4*mpi2)
    ctr1dJ0(3,1,23)=(1.36309862136841e-7_pr*c3*gA2)/(fpi4*mpi2)- &
         (4.54366207122803e-8_pr*c4*gA2)/(fpi4*mpi2)
    !
    !  UJ1Drho0
    !
    ctJ1dr0=0.0_pr
    !
    !  UJ0J0
    !
    ctJ0J0=0.0_pr
    ctJ0J0(1,1,3)=(-0.0388151948706056_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0689519291560374_pr*c4*gA2)/(fpi4*mpi2)+ &
         (0.03515625_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctJ0J0(1,1,5)=(-0.148030393439646_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.125249890279611_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.124165184550383_pr*c4*gA2)/(fpi4*mpi2)- &
         (0.017578125_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctJ0J0(1,2,6)=(0.266832652280755_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0685045280612245_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(1,1,7)=(-50.9400757481906_pr*c1*gA2)/(fpi4*mpi2)+ &
         (1.49639517552504_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00598016955317283_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(2,1,7)=(9.4436896594477_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.387304024062572_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00211338089923469_pr*c4*gA2)/(fpi4*mpi2)+ &
         (0.00439453125_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctJ0J0(1,2,8)=(60.4392047247024_pr*c1*gA2)/(fpi4*mpi2)+ &
         (21.2791752723318_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0653900223214286_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(1,1,9)=(-6519.80035123157_pr*c1*gA2)/(fpi4*mpi2)+ &
         (316.045390495029_pr*c3*gA2)/(fpi4*mpi2)+ &
         (2.16832575424506_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(2,1,9)=(2094.21472214372_pr*c1*gA2)/(fpi4*mpi2)- &
         (79.8751918140941_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.347175990513393_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(3,1,9)=(-174.923789099377_pr*c1*gA2)/(fpi4*mpi2)
    ctJ0J0(1,3,9)=(-136.232770647321_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.632829241071429_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(1,2,10)=(22280.2381832589_pr*c1*gA2)/(fpi4*mpi2)- &
         (3725.75432559614_pr*c3*gA2)/(fpi4*mpi2)- &
         (40.7846229073661_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(2,2,10)=(-3482.03035714286_pr*c1*gA2)/(fpi4*mpi2)+ &
         (857.700259443584_pr*c3*gA2)/(fpi4*mpi2)+ &
         (7.79614955357143_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(1,1,11)=(52284.4157057857_pr*c1*gA2)/(fpi4*mpi2)- &
         (22649.0853151035_pr*c3*gA2)/(fpi4*mpi2)- &
         (389.622252285779_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(2,1,11)=(-29244.8802258893_pr*c1*gA2)/(fpi4*mpi2)+ &
         (11240.6895935096_pr*c3*gA2)/(fpi4*mpi2)+ &
         (170.39136171439_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(3,1,11)=(3224.34088479694_pr*c1*gA2)/(fpi4*mpi2)- &
         (1377.64430659016_pr*c3*gA2)/(fpi4*mpi2)- &
         (18.1685616629464_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(1,3,11)=(-20240.849921875_pr*c1*gA2)/(fpi4*mpi2)+ &
         (3674.53847265625_pr*c3*gA2)/(fpi4*mpi2)+ &
         (49.9097109375_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(1,2,12)=(-142538.810602121_pr*c1*gA2)/(fpi4*mpi2)+ &
         (51959.7177006278_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1056.48083035714_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(2,2,12)=(47556.1071386719_pr*c1*gA2)/(fpi4*mpi2)- &
         (13973.2343116668_pr*c3*gA2)/(fpi4*mpi2)- &
         (251.396513671875_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(1,1,13)=(-141614.511333235_pr*c1*gA2)/(fpi4*mpi2)+ &
         (35973.0682685409_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1263.90306398228_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(2,1,13)=(120423.751289578_pr*c1*gA2)/(fpi4*mpi2)- &
         (33366.8468838635_pr*c3*gA2)/(fpi4*mpi2)- &
         (938.242939340534_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(3,1,13)=(-23741.2058036637_pr*c1*gA2)/(fpi4*mpi2)+ &
         (6070.38539035745_pr*c3*gA2)/(fpi4*mpi2)+ &
         (145.681594848633_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(1,3,13)=(76403.207390625_pr*c1*gA2)/(fpi4*mpi2)- &
         (27942.3675271205_pr*c3*gA2)/(fpi4*mpi2)- &
         (671.186349441964_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(1,2,14)=(204907.907447184_pr*c1*gA2)/(fpi4*mpi2)- &
         (56211.484659242_pr*c3*gA2)/(fpi4*mpi2)- &
         (2422.99460657436_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(2,2,14)=(-95398.9531424386_pr*c1*gA2)/(fpi4*mpi2)+ &
         (30129.2960691746_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1018.34175868443_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(1,1,15)=(9978.05620114632_pr*c1*gA2)/(fpi4*mpi2)- &
         (6474.12892378917_pr*c3*gA2)/(fpi4*mpi2)- &
         (968.323984416389_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(2,1,15)=(-36487.3307862002_pr*c1*gA2)/(fpi4*mpi2)+ &
         (13937.926974747_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1128.44835742705_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(3,1,15)=(13562.8790229223_pr*c1*gA2)/(fpi4*mpi2)- &
         (5370.52492057074_pr*c3*gA2)/(fpi4*mpi2)- &
         (292.012623596191_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(1,3,15)=(-76061.0427714844_pr*c1*gA2)/(fpi4*mpi2)+ &
         (21292.048190918_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1052.55900878906_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(1,2,16)=(-26755.8098873901_pr*c1*gA2)/(fpi4*mpi2)+ &
         (11051.0158225074_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1284.3833242048_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(2,2,16)=(33104.1130819702_pr*c1*gA2)/(fpi4*mpi2)- &
         (10790.492615363_pr*c3*gA2)/(fpi4*mpi2)- &
         (813.688233476911_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(1,1,17)=(-13525.0298722901_pr*c1*gA2)/(fpi4*mpi2)+ &
         (3920.12706464951_pr*c3*gA2)/(fpi4*mpi2)+ &
         (21.2402742906025_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(2,1,17)=(18744.3022469218_pr*c1*gA2)/(fpi4*mpi2)- &
         (5272.85370859068_pr*c3*gA2)/(fpi4*mpi2)- &
         (174.676497182301_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(3,1,17)=(-7131.99279047117_pr*c1*gA2)/(fpi4*mpi2)+ &
         (1930.52177440425_pr*c3*gA2)/(fpi4*mpi2)+ &
         (88.4860045787266_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(1,3,17)=(7839.31645019531_pr*c1*gA2)/(fpi4*mpi2)- &
         (2869.73837755476_pr*c3*gA2)/(fpi4*mpi2)- &
         (419.970217895508_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(1,2,18)=(355.059586151123_pr*c1*gA2)/(fpi4*mpi2)+ &
         (183.100210884983_pr*c3*gA2)/(fpi4*mpi2)- &
         (13.999078478132_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(2,2,18)=(-1327.89377142334_pr*c1*gA2)/(fpi4*mpi2)+ &
         (15.4343164854196_pr*c3*gA2)/(fpi4*mpi2)+ &
         (105.873926083374_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(1,1,19)=(-8867.21170169755_pr*c1*gA2)/(fpi4*mpi2)+ &
         (2936.70582294006_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.109368735306604_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(2,1,19)=(13819.8747331675_pr*c1*gA2)/(fpi4*mpi2)- &
         (4438.03614592496_pr*c3*gA2)/(fpi4*mpi2)- &
         (3.51122917092868_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(3,1,19)=(-4725.00077792637_pr*c1*gA2)/(fpi4*mpi2)+ &
         (1532.04177968309_pr*c3*gA2)/(fpi4*mpi2)- &
         (4.33029295563698_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(2,2,20)=(-88.7648965377808_pr*c1*gA2)/(fpi4*mpi2)- &
         (45.7750527212457_pr*c3*gA2)/(fpi4*mpi2)+ &
         (3.49976961953299_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(1,1,21)=(-3950.27983470918_pr*c1*gA2)/(fpi4*mpi2)+ &
         (1347.42215988918_pr*c3*gA2)/(fpi4*mpi2)
    ctJ0J0(2,1,21)=(7835.22431087985_pr*c1*gA2)/(fpi4*mpi2)- &
         (2624.64064001171_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0546843676533018_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(3,1,21)=(-3337.91716486204_pr*c1*gA2)/(fpi4*mpi2)+ &
         (1111.14574092555_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.449709850430489_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(1,1,23)=(-1142.83677022112_pr*c1*gA2)/(fpi4*mpi2)+ &
         (401.209483659988_pr*c3*gA2)/(fpi4*mpi2)
    ctJ0J0(2,1,23)=(3005.12399409001_pr*c1*gA2)/(fpi4*mpi2)- &
         (1034.46527332476_pr*c3*gA2)/(fpi4*mpi2)
    ctJ0J0(3,1,23)=(-1635.837501786_pr*c1*gA2)/(fpi4*mpi2)+ &
         (553.405150406837_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00683554595666272_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J0(1,1,25)=(-207.173615193637_pr*c1*gA2)/(fpi4*mpi2)+ &
         (74.5943777265557_pr*c3*gA2)/(fpi4*mpi2)
    ctJ0J0(2,1,25)=(765.620232650329_pr*c1*gA2)/(fpi4*mpi2)- &
         (270.43438082633_pr*c3*gA2)/(fpi4*mpi2)
    ctJ0J0(3,1,25)=(-549.802354657625_pr*c1*gA2)/(fpi4*mpi2)+ &
         (190.708530684211_pr*c3*gA2)/(fpi4*mpi2)
    ctJ0J0(1,1,27)=(-21.371516474291_pr*c1*gA2)/(fpi4*mpi2)+ &
         (7.86480912750357_pr*c3*gA2)/(fpi4*mpi2)
    ctJ0J0(2,1,27)=(124.32036907328_pr*c1*gA2)/(fpi4*mpi2)- &
         (44.92276431339_pr*c3*gA2)/(fpi4*mpi2)
    ctJ0J0(3,1,27)=(-125.00166164342_pr*c1*gA2)/(fpi4*mpi2)+ &
         (44.3795879210135_pr*c3*gA2)/(fpi4*mpi2)
    ctJ0J0(1,1,29)=(-0.956932496743248_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.358850516087132_pr*c3*gA2)/(fpi4*mpi2)
    ctJ0J0(2,1,29)=(11.6426907338887_pr*c1*gA2)/(fpi4*mpi2)- &
         (4.29125507983892_pr*c3*gA2)/(fpi4*mpi2)
    ctJ0J0(3,1,29)=(-18.3709744429036_pr*c1*gA2)/(fpi4*mpi2)+ &
         (6.65825509945955_pr*c3*gA2)/(fpi4*mpi2)
    ctJ0J0(2,1,31)=(0.478466248371624_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.179425258043566_pr*c3*gA2)/(fpi4*mpi2)
    ctJ0J0(3,1,31)=(-1.574952903829_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.581263199490756_pr*c3*gA2)/(fpi4*mpi2)
    ctJ0J0(3,1,33)=(-0.059808281046453_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.0224281572554458_pr*c3*gA2)/(fpi4*mpi2)
    !
    !  UJ0J1
    !
    ctJ0J1=0.0_pr
    ctJ0J1(1,1,3)=(0.0285278380152656_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0459679527706916_pr*c4*gA2)/(fpi4*mpi2)- &
         (0.0234375_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctJ0J1(1,1,5)=(-0.00310602359693878_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0516361625744048_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0827767897002551_pr*c4*gA2)/(fpi4*mpi2)+ &
         (0.01171875_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctJ0J1(1,2,6)=(0.0475446853741497_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0456696853741497_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(1,1,7)=(0.147463309151786_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0358081265611713_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00398677970211522_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(2,1,7)=(0.00200860969387755_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0111190874787415_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0014089205994898_pr*c4*gA2)/(fpi4*mpi2)- &
         (0.0029296875_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctJ0J1(1,2,8)=(-0.363757366071429_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.177937700892857_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0435933482142857_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(1,1,9)=(11.9185698381696_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.18312424087213_pr*c3*gA2)/(fpi4*mpi2)- &
         (1.44555050283004_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(2,1,9)=(-2.81250837053571_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0791064453125_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.231450660342262_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(3,1,9)=(0.17578125_pr*c1*gA2)/(fpi4*mpi2)
    ctJ0J1(1,3,9)=(-0.210943080357143_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.421886160714286_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(1,2,10)=(-90.1323160714286_pr*c1*gA2)/(fpi4*mpi2)+ &
         (9.21932595663265_pr*c3*gA2)/(fpi4*mpi2)+ &
         (27.1897486049107_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(2,2,10)=(12.6565848214286_pr*c1*gA2)/(fpi4*mpi2)- &
         (1.095703125_pr*c3*gA2)/(fpi4*mpi2)- &
         (5.19743303571429_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(1,1,11)=(-798.881085665458_pr*c1*gA2)/(fpi4*mpi2)+ &
         (160.074257447634_pr*c3*gA2)/(fpi4*mpi2)+ &
         (259.748168190519_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(2,1,11)=(328.414995954241_pr*c1*gA2)/(fpi4*mpi2)- &
         (61.5935256713867_pr*c3*gA2)/(fpi4*mpi2)- &
         (113.594241142927_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(3,1,11)=(-30.7629743303571_pr*c1*gA2)/(fpi4*mpi2)+ &
         (5.57625906808036_pr*c3*gA2)/(fpi4*mpi2)+ &
         (12.1123744419643_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(1,3,11)=(123.05109375_pr*c1*gA2)/(fpi4*mpi2)- &
         (15.2580234375_pr*c3*gA2)/(fpi4*mpi2)- &
         (33.273140625_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(1,2,12)=(2557.41305524554_pr*c1*gA2)/(fpi4*mpi2)- &
         (494.193056919643_pr*c3*gA2)/(fpi4*mpi2)- &
         (704.320553571429_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(2,2,12)=(-600.61833984375_pr*c1*gA2)/(fpi4*mpi2)+ &
         (105.762021484375_pr*c3*gA2)/(fpi4*mpi2)+ &
         (167.59767578125_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(1,1,13)=(4226.57441403285_pr*c1*gA2)/(fpi4*mpi2)- &
         (759.654990867397_pr*c3*gA2)/(fpi4*mpi2)- &
         (842.602042654855_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(2,1,13)=(-2859.96228445522_pr*c1*gA2)/(fpi4*mpi2)+ &
         (519.521237971814_pr*c3*gA2)/(fpi4*mpi2)+ &
         (625.495292893689_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(3,1,13)=(423.723815917969_pr*c1*gA2)/(fpi4*mpi2)- &
         (74.1513366699219_pr*c3*gA2)/(fpi4*mpi2)- &
         (97.1210632324219_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(1,3,13)=(-1778.06559375_pr*c1*gA2)/(fpi4*mpi2)+ &
         (349.290122544643_pr*c3*gA2)/(fpi4*mpi2)+ &
         (447.457566294643_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(1,2,14)=(-8387.96795947266_pr*c1*gA2)/(fpi4*mpi2)+ &
         (1557.82783909738_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1615.32973771624_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(2,2,14)=(3204.65585239955_pr*c1*gA2)/(fpi4*mpi2)- &
         (610.780761369978_pr*c3*gA2)/(fpi4*mpi2)- &
         (678.894505789621_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(1,1,15)=(-3814.31004524231_pr*c1*gA2)/(fpi4*mpi2)+ &
         (782.002794581876_pr*c3*gA2)/(fpi4*mpi2)+ &
         (645.54932294426_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(2,1,15)=(4206.02934884295_pr*c1*gA2)/(fpi4*mpi2)- &
         (834.729188615308_pr*c3*gA2)/(fpi4*mpi2)- &
         (752.298904951368_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(3,1,15)=(-1006.82687988281_pr*c1*gA2)/(fpi4*mpi2)+ &
         (200.402084350586_pr*c3*gA2)/(fpi4*mpi2)+ &
         (194.675082397461_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(1,3,15)=(3781.37963671875_pr*c1*gA2)/(fpi4*mpi2)- &
         (709.003098632812_pr*c3*gA2)/(fpi4*mpi2)- &
         (701.706005859375_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(1,2,16)=(5207.17552185059_pr*c1*gA2)/(fpi4*mpi2)- &
         (1065.41645752128_pr*c3*gA2)/(fpi4*mpi2)- &
         (856.255549469866_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(2,2,16)=(-3152.00038513184_pr*c1*gA2)/(fpi4*mpi2)+ &
         (618.68313205392_pr*c3*gA2)/(fpi4*mpi2)+ &
         (542.458822317941_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(1,1,17)=(114.489910766602_pr*c1*gA2)/(fpi4*mpi2)- &
         (21.7326762047359_pr*c3*gA2)/(fpi4*mpi2)- &
         (14.1601828604017_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(2,1,17)=(-785.971121017456_pr*c1*gA2)/(fpi4*mpi2)+ &
         (158.614487102509_pr*c3*gA2)/(fpi4*mpi2)+ &
         (116.450998121534_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(3,1,17)=(379.629753897531_pr*c1*gA2)/(fpi4*mpi2)- &
         (73.6791959381104_pr*c3*gA2)/(fpi4*mpi2)- &
         (58.9906697191511_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(1,3,17)=(-1742.07032226563_pr*c1*gA2)/(fpi4*mpi2)+ &
         (357.745064784459_pr*c3*gA2)/(fpi4*mpi2)+ &
         (279.980145263672_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(1,2,18)=(-76.2155765991211_pr*c1*gA2)/(fpi4*mpi2)+ &
         (14.7764228733608_pr*c3*gA2)/(fpi4*mpi2)+ &
         (9.33271898542132_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(2,2,18)=(478.384230102539_pr*c1*gA2)/(fpi4*mpi2)- &
         (98.7791618408203_pr*c3*gA2)/(fpi4*mpi2)- &
         (70.582617388916_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(1,1,19)=(-0.68049621963501_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.150680578790392_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0729124902044024_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(2,1,19)=(-18.4566708641052_pr*c1*gA2)/(fpi4*mpi2)+ &
         (3.32744608689717_pr*c3*gA2)/(fpi4*mpi2)+ &
         (2.34081944728579_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(3,1,19)=(-19.8544406890869_pr*c1*gA2)/(fpi4*mpi2)+ &
         (4.63924562215805_pr*c3*gA2)/(fpi4*mpi2)+ &
         (2.88686197042465_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(2,2,20)=(19.0538941497803_pr*c1*gA2)/(fpi4*mpi2)- &
         (3.69410571834019_pr*c3*gA2)/(fpi4*mpi2)- &
         (2.33317974635533_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(2,1,21)=(0.340248109817505_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0753402893951961_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0364562451022012_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(3,1,21)=(-2.54145170688629_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.526430741071701_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.299806566953659_pr*c4*gA2)/(fpi4*mpi2)
    ctJ0J1(3,1,23)=(-0.0425310137271881_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.00941753617439951_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00455703063777515_pr*c4*gA2)/(fpi4*mpi2)
    !
    !  UJ1J1
    !
    ctJ1J1=0.0_pr
    ctJ1J1(1,1,3)=(0.0496322121163276_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0229839763853458_pr*c4*gA2)/(fpi4*mpi2)- &
         (0.01171875_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctJ1J1(1,1,5)=(0.196140339479882_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.0113169982018643_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0413883948501276_pr*c4*gA2)/(fpi4*mpi2)+ &
         (0.005859375_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctJ1J1(1,2,6)=(-0.0889442174269183_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0228348426870748_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(1,1,7)=(16.7237752493969_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.34790459899251_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.00199338985105761_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(2,1,7)=(-3.22524030314923_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.137722860274259_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.000704460299744898_pr*c4*gA2)/(fpi4*mpi2)- &
         (0.00146484375_pr*cd*ga)/(fpi4*mpi2*LambdaX)
    ctJ1J1(1,2,8)=(-19.6214015749008_pr*c1*gA2)/(fpi4*mpi2)- &
         (7.57700061161061_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0217966741071429_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(1,1,9)=(2190.64598296261_pr*c1*gA2)/(fpi4*mpi2)- &
         (106.868359765483_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.722775251415019_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(2,1,9)=(-703.063761547905_pr*c1*gA2)/(fpi4*mpi2)+ &
         (27.158281478954_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.115725330171131_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(3,1,9)=(58.6594921997924_pr*c1*gA2)/(fpi4*mpi2)
    ctJ1J1(1,3,9)=(46.67658203125_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.210943080357143_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(1,2,10)=(-7494.70074858631_pr*c1*gA2)/(fpi4*mpi2)+ &
         (1278.61760582159_pr*c3*gA2)/(fpi4*mpi2)+ &
         (13.5948743024554_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(2,2,10)=(1169.53616071429_pr*c1*gA2)/(fpi4*mpi2)- &
         (294.126983802623_pr*c3*gA2)/(fpi4*mpi2)- &
         (2.59871651785714_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(1,1,11)=(-17607.1420353921_pr*c1*gA2)/(fpi4*mpi2)+ &
         (7994.24713374194_pr*c3*gA2)/(fpi4*mpi2)+ &
         (129.87408409526_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(2,1,11)=(9843.68403362975_pr*c1*gA2)/(fpi4*mpi2)- &
         (3932.01024951878_pr*c3*gA2)/(fpi4*mpi2)- &
         (56.7971205714634_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(3,1,11)=(-1084.01759961981_pr*c1*gA2)/(fpi4*mpi2)+ &
         (478.820795649102_pr*c3*gA2)/(fpi4*mpi2)+ &
         (6.05618722098214_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(1,3,11)=(6819.79372395833_pr*c1*gA2)/(fpi4*mpi2)- &
         (1268.16076692708_pr*c3*gA2)/(fpi4*mpi2)- &
         (16.6365703125_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(1,2,12)=(48059.0188986235_pr*c1*gA2)/(fpi4*mpi2)- &
         (18566.944349986_pr*c3*gA2)/(fpi4*mpi2)- &
         (352.160276785714_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(2,2,12)=(-16034.3911816406_pr*c1*gA2)/(fpi4*mpi2)+ &
         (4929.07689946187_pr*c3*gA2)/(fpi4*mpi2)+ &
         (83.798837890625_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(1,1,13)=(47847.3513627872_pr*c1*gA2)/(fpi4*mpi2)- &
         (13686.6719432792_pr*c3*gA2)/(fpi4*mpi2)- &
         (421.301021327427_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(2,1,13)=(-40664.4224147226_pr*c1*gA2)/(fpi4*mpi2)+ &
         (12306.5684431232_pr*c3*gA2)/(fpi4*mpi2)+ &
         (312.747646446845_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(3,1,13)=(8015.01712823947_pr*c1*gA2)/(fpi4*mpi2)- &
         (2193.04572012566_pr*c3*gA2)/(fpi4*mpi2)- &
         (48.5605316162109_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(1,3,13)=(-25770.923296875_pr*c1*gA2)/(fpi4*mpi2)+ &
         (10166.4482434152_pr*c3*gA2)/(fpi4*mpi2)+ &
         (223.728783147321_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(1,2,14)=(-69245.3594485406_pr*c1*gA2)/(fpi4*mpi2)+ &
         (22270.458788223_pr*c3*gA2)/(fpi4*mpi2)+ &
         (807.664868858119_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(2,2,14)=(32217.3951881045_pr*c1*gA2)/(fpi4*mpi2)- &
         (11458.4663920128_pr*c3*gA2)/(fpi4*mpi2)- &
         (339.44725289481_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(1,1,15)=(-3355.18276569786_pr*c1*gA2)/(fpi4*mpi2)+ &
         (4003.15120252837_pr*c3*gA2)/(fpi4*mpi2)+ &
         (322.77466147213_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(2,1,15)=(12334.8284220602_pr*c1*gA2)/(fpi4*mpi2)- &
         (6594.07536966716_pr*c3*gA2)/(fpi4*mpi2)- &
         (376.149452475684_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(3,1,15)=(-4584.16347069416_pr*c1*gA2)/(fpi4*mpi2)+ &
         (2259.89314979311_pr*c3*gA2)/(fpi4*mpi2)+ &
         (97.3375411987305_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(1,3,15)=(25694.7668613281_pr*c1*gA2)/(fpi4*mpi2)- &
         (8712.71905517578_pr*c3*gA2)/(fpi4*mpi2)- &
         (350.853002929688_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(1,2,16)=(8939.92116689046_pr*c1*gA2)/(fpi4*mpi2)- &
         (6189.35971202371_pr*c3*gA2)/(fpi4*mpi2)- &
         (428.127774734933_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(2,2,16)=(-11150.7683255005_pr*c1*gA2)/(fpi4*mpi2)+ &
         (5032.6588526575_pr*c3*gA2)/(fpi4*mpi2)+ &
         (271.22941115897_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(1,1,17)=(4508.67638249921_pr*c1*gA2)/(fpi4*mpi2)- &
         (1354.76069989184_pr*c3*gA2)/(fpi4*mpi2)- &
         (7.08009143020085_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(2,1,17)=(-6244.51076026543_pr*c1*gA2)/(fpi4*mpi2)+ &
         (2125.28532690204_pr*c3*gA2)/(fpi4*mpi2)+ &
         (58.225499060767_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(3,1,17)=(2384.7618078426_pr*c1*gA2)/(fpi4*mpi2)- &
         (812.91389298321_pr*c3*gA2)/(fpi4*mpi2)- &
         (29.4953348595755_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(1,3,17)=(-2613.10548339844_pr*c1*gA2)/(fpi4*mpi2)+ &
         (1796.50318320138_pr*c3*gA2)/(fpi4*mpi2)+ &
         (139.990072631836_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(1,2,18)=(-118.353195383708_pr*c1*gA2)/(fpi4*mpi2)- &
         (28.3697180204175_pr*c3*gA2)/(fpi4*mpi2)+ &
         (4.66635949271066_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(2,2,18)=(437.301789367676_pr*c1*gA2)/(fpi4*mpi2)- &
         (234.978396185244_pr*c3*gA2)/(fpi4*mpi2)- &
         (35.291308694458_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(1,1,19)=(2955.73723389918_pr*c1*gA2)/(fpi4*mpi2)- &
         (978.581137073306_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.0364562451022012_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(2,1,19)=(-4606.79145692376_pr*c1*gA2)/(fpi4*mpi2)+ &
         (1486.71857443532_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1.17040972364289_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(3,1,19)=(1575.92551413057_pr*c1*gA2)/(fpi4*mpi2)- &
         (499.689660100295_pr*c3*gA2)/(fpi4*mpi2)+ &
         (1.44343098521233_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(2,2,20)=(29.5882988459269_pr*c1*gA2)/(fpi4*mpi2)+ &
         (7.09242950510437_pr*c3*gA2)/(fpi4*mpi2)- &
         (1.16658987317766_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(1,1,21)=(1316.75994490306_pr*c1*gA2)/(fpi4*mpi2)- &
         (449.140719963059_pr*c3*gA2)/(fpi4*mpi2)
    ctJ1J1(2,1,21)=(-2611.74143695995_pr*c1*gA2)/(fpi4*mpi2)+ &
         (874.719811383881_pr*c3*gA2)/(fpi4*mpi2)- &
         (0.0182281225511006_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(3,1,21)=(1112.6598731875_pr*c1*gA2)/(fpi4*mpi2)- &
         (369.221981860559_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.14990328347683_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(1,1,23)=(380.945590073708_pr*c1*gA2)/(fpi4*mpi2)- &
         (133.736494553329_pr*c3*gA2)/(fpi4*mpi2)
    ctJ1J1(2,1,23)=(-1001.70799803_pr*c1*gA2)/(fpi4*mpi2)+ &
         (344.82175777492_pr*c3*gA2)/(fpi4*mpi2)
    ctJ1J1(3,1,23)=(545.279167261999_pr*c1*gA2)/(fpi4*mpi2)- &
         (184.448333224776_pr*c3*gA2)/(fpi4*mpi2)+ &
         (0.00227851531888757_pr*c4*gA2)/(fpi4*mpi2)
    ctJ1J1(1,1,25)=(69.0578717312122_pr*c1*gA2)/(fpi4*mpi2)- &
         (24.8647925755186_pr*c3*gA2)/(fpi4*mpi2)
    ctJ1J1(2,1,25)=(-255.206744216776_pr*c1*gA2)/(fpi4*mpi2)+ &
         (90.1447936087768_pr*c3*gA2)/(fpi4*mpi2)
    ctJ1J1(3,1,25)=(183.267451552542_pr*c1*gA2)/(fpi4*mpi2)- &
         (63.5695102280704_pr*c3*gA2)/(fpi4*mpi2)
    ctJ1J1(1,1,27)=(7.12383882476366_pr*c1*gA2)/(fpi4*mpi2)- &
         (2.62160304250119_pr*c3*gA2)/(fpi4*mpi2)
    ctJ1J1(2,1,27)=(-41.4401230244268_pr*c1*gA2)/(fpi4*mpi2)+ &
         (14.97425477113_pr*c3*gA2)/(fpi4*mpi2)
    ctJ1J1(3,1,27)=(41.6672205478068_pr*c1*gA2)/(fpi4*mpi2)- &
         (14.7931959736712_pr*c3*gA2)/(fpi4*mpi2)
    ctJ1J1(1,1,29)=(0.318977498914416_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.119616838695711_pr*c3*gA2)/(fpi4*mpi2)
    ctJ1J1(2,1,29)=(-3.88089691129625_pr*c1*gA2)/(fpi4*mpi2)+ &
         (1.43041835994631_pr*c3*gA2)/(fpi4*mpi2)
    ctJ1J1(3,1,29)=(6.12365814763455_pr*c1*gA2)/(fpi4*mpi2)- &
         (2.21941836648652_pr*c3*gA2)/(fpi4*mpi2)
    ctJ1J1(2,1,31)=(-0.159488749457208_pr*c1*gA2)/(fpi4*mpi2)+ &
         (0.0598084193478553_pr*c3*gA2)/(fpi4*mpi2)
    ctJ1J1(3,1,31)=(0.524984301276333_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.193754399830252_pr*c3*gA2)/(fpi4*mpi2)
    ctJ1J1(3,1,33)=(0.019936093682151_pr*c1*gA2)/(fpi4*mpi2)- &
         (0.00747605241848192_pr*c3*gA2)/(fpi4*mpi2)
  End Subroutine load_tables
  !=======================================================================
  !> Preset default values for Skyrme-like functionals
  !=======================================================================
  Subroutine default_UNEDF_NAMELIST()
    Implicit None
    ! Parameters
    eps     = Spacing(one)
    Pi      = four*Atan(one)
    kfconst = (1.5_pr*PI**2)**(one/three)
    CK      = three/five*kfconst**2
    ! Some default values for all cases
    Print_Namelist=.False.
    FunctionalName="Bla-Bla"
    ! Kind of the functional
    use_INM                = .False.
    use_DME3N_terms        = .False.
    use_charge_density     = .False.
    use_cm_cor             = .False.
    use_full_cm_cor        = .False.
    use_j2terms            = .False.
    hb0_charge_dependent   = .False.
    finite_range           = .False.
    force_is_DME           = .False.
    override_3N_couplings  = .False.
    alfa_dme               =  zero
    is_nedf                = .False.
    coulomb_gaussian       = .False.
    use_TMR_pairing        =  0
    DMEorder               = -1
    DMElda                 =  0
    ! Coupling constants: ph channel (UNEDF1 values)
    Crho(0)=   -779.373008720865300_pr; Crho(1)  =   287.7221315832867960_pr
    CDrho(0)=   891.477890442349690_pr; CDrho(1) =  -200.5877743178848790_pr
    Ctau(0)=  -0.989915057807676746_pr; Ctau(1)  =   -33.6320970701835549_pr
    CrDr(0)=   -45.1351310222373030_pr; CrDr(1)  =  -145.3821679080569990_pr
    CrdJ(0)=   -74.0263331764599002_pr; CrdJ(1)  =   -35.6582611147916992_pr
    CJ(0)=      0.00000000000000000_pr; CJ(1)    =     0.0000000000000000_pr
    ! Coupling constants: pp channel (UNEDF1 values)
    CpV0(0)=   -186.065399575123990_pr; CpV0(1)  =  -206.5795938901059970_pr
    CpV1(0)=   0.500000000000000000_pr; CpV1(1)  =     0.5000000000000000_pr
    ! Various
    Cnrho    =    0.000000000000000_pr; CJdr     =     0.0000000000000000_pr
    sigma    =    0.2700180115027076_pr;
    hbzero   =   20.7355300000000007_pr;
    e2charg  =    1.4399784085965135_pr ; CExPar = 1.0_pr
    ! DME
    mpi=   138.03_pr/197.3_pr; fpi = 92.4_pr/197.3_pr; gA = 1.29_pr
    c1 =    -0.81_pr/1000.0_pr*197.3_pr
    c3 =    -3.40_pr/1000.0_pr*197.3_pr
    c4 =     3.40_pr/1000.0_pr*197.3_pr
    cd = -2062.00_pr/1000.0_pr
    ce =  -625.00_pr/1000.0_pr
    ! Natural units
    LambdaX  = 700.0_pr/197.3_pr; nuLambda = 700.0_pr; nufpi = 93.0_pr
    ! Nuclear matter (UNEDF1 values)
    E_NM     = -15.80000000000000000_pr; RHO_NM   =  0.158706769332587039_pr
    K_NM     = 220.00000000000002800_pr; SMASS_NM =  0.992423332283363990_pr
    ASS_NM   =  28.98678905777210700_pr; LASS_NM  = 40.004790480413632300_pr
    VMASS_NM =   1.24983857423226996_pr
    ! NEDF
    a_NEDF   = zero
    b_NEDF   = zero
    c_NEDF   = zero
    eta_NEDF = zero
    W0_NEDF  = zero
  End Subroutine default_UNEDF_NAMELIST
  !===============================================================================================
  !
  !===============================================================================================
  Subroutine read_UNEDF_NAMELIST(fname,noForces,filename_nml)
    Use HFBTHO_utilities, Only: lout
    !--------------------------------------------------------------------------------
    ! RESERVED NAMES ARE:
    !  -namelist forbiden:
    !          'UNEDF'  - best UNEDF
    !          'SKYRME' - best SKYRME
    !  -namelist inforced but not for C-parameters (use_INM=F)
    !   or NM-parameters (use_INM=T) defined by the solver
    !          'FITS'
    !  -namelist inforced (one can overwrite all):
    !          'ANY OTHER NAME'
    ! i.e., the DME solver defines C-/NM- only using 'FITS'
    !--------------------------------------------------------------------------------
    Implicit None
    Character (30), Intent(inout)          :: fname
    Integer(ipr),   Intent(out)            :: noForces
    Character (*),  Intent(in),   Optional :: filename_nml
    Character (30) :: inforcedname
    Logical        :: regularization
    Integer(ipr)   :: ios,lnamelist=16
    ! TODO: Make this length a program-wide constant
    Character(256) :: filename

    If (present(filename_nml)) then
        filename = filename_nml
    Else
        filename = 'UNEDF_NAMELIST.dat'
    End If
    !
    ! parameters
    eps     = Spacing(1.0_pr)
    Pi      = 4.0_pr*Atan(1.0_pr)
    kfconst =(1.50_pr*Pi**2)**(1.0_pr/3.0_pr)    ! (3Pi^2/2)^(1/3)
    CK      = 3.0_pr/5.0_pr*kfconst**2
    !
    use_Namelist=.True.
    Do
       !---------------------------------------------------------------------
       ! Some default values for all cases
       !---------------------------------------------------------------------
       Call default_UNEDF_NAMELIST()
       !---------------------------------------------------------------------
       ! Select the functional: start with interaction
       !---------------------------------------------------------------------
       noForces=0 ! No forces to start with
       Call skforce(fname,noForces)
       !
       If (noForces.Eq.1) Then
           inforcedname='FORCE'
           use_Namelist=.False.
       Else
           FUNCTIONAL: Select Case (Trim(fname))
           Case ('READ')
              inforcedname='READ'
              use_Namelist=.False.
              ! Set to the same
              use_INM  = .True.
              ! Surface coefficients
              CrDr(0)  =  -45.599763_pr
              CrDr(1)  = -143.935074_pr
              CpV0(0)  = -234.380010_pr
              CpV0(1)  = -260.437001_pr
              CrdJ(0)  =  -73.946388_pr
              CrdJ(1)  =  -51.912548_pr
              CJ(0)    =    0.000000_pr
              CJ(1)    =    0.000000_pr
              CExPar   =    1.000000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.156247_pr
              E_NM     =  -15.800000_pr
              P_NM     =    0.000000_pr
              K_NM     =  244.668204_pr
              ASS_NM   =   28.668204_pr
              LASS_NM  =   40.109081_pr
              SMASS_NM =    1.067970_pr
              VMASS_NM =    1.249838574232270_pr
              if(use_3N_couplings) then
                 use_3N_couplings = .False.
                 override_3N_couplings = .True.
              endif
           Case ('RAND')
              inforcedname='RAND'
              use_Namelist=.False.
              ! Set to the same
              use_INM  = .True.
              ! Surface coefficients
              CrDr(0)  =  -45.599763_pr
              CrDr(1)  = -143.935074_pr
              CpV0(0)  = -234.380010_pr
              CpV0(1)  = -260.437001_pr
              CrdJ(0)  =  -73.946388_pr
              CrdJ(1)  =  -51.912548_pr
              CJ(0)    =    0.000000_pr
              CJ(1)    =    0.000000_pr
              CExPar   =    1.000000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.156247_pr
              E_NM     =  -15.800000_pr
              P_NM     =    0.000000_pr
              K_NM     =  244.668204_pr
              ASS_NM   =   28.668204_pr
              LASS_NM  =   40.109081_pr
              SMASS_NM =    1.067970_pr
              VMASS_NM =    1.249838574232270_pr
              if(use_3N_couplings) then
                 use_3N_couplings = .False.
                 override_3N_couplings = .True.
              endif
           Case ('FITS')
              inforcedname='FITS'
              use_Namelist=.True.
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              P_NM     =    0.000000000000000_pr
              VMASS_NM =    1.249838000000000_pr
              !
              force_is_dme   = .False.
              finite_range   = .False.
              use_cm_cor     = .False.
              use_j2terms    = .True.
              CExPar         = 1.0000000000000000_pr
              use_INM        = .True.
              VMASS_NM       =    0.9916634575083850_pr
              RHO_NM         =    0.1545314206232500_pr
              E_NM           =  -14.9428151261574005_pr
              K_NM           =  254.6493849869999906_pr
              ASS_NM         =   28.2555644286320984_pr
              LASS_NM        =   64.9345887215435056_pr
              SMASS_NM       =    0.8026870082347940_pr
              CrDr(0)        =  -40.1069828003495985_pr
              CrDr(1)        =  -35.8476436441168005_pr
              CpV0(0)        = -233.9624836153689955_pr
              CpV0(1)        = -284.1118846887910081_pr
              CrdJ(0)        =  -74.0943831060370002_pr
              CrdJ(1)        =   10.3610837684245993_pr
              CJ(0)          =  -81.3789370214482943_pr
              CJ(1)          =  -66.9270694667199990_pr
            Case ('UNE0')
              inforcedname='UNE0'
              use_Namelist=.False.
              ! kind of the functional
              use_INM    = .True.
              use_cm_cor = .True.
              ! Surface coefficients
              CrDr(0)  =  -55.260600000000000_pr
              CrDr(1)  =  -55.622600000000000_pr
              CpV0(0)  = -170.374000000000000_pr
              CpV0(1)  = -199.202000000000000_pr
              CrdJ(0)  =  -79.530800000000000_pr
              CrdJ(1)  =   45.630200000000000_pr
              CJ(0)    =    0.000000000000000_pr
              CJ(1)    =    0.000000000000000_pr
              CExPar   =    1.000000000000000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.160526000000000_pr
              E_NM     =  -16.055900000000000_pr
              P_NM     =    0.000000000000000_pr
              K_NM     =  230.000000000000000_pr
              ASS_NM   =   30.542900000000000_pr
              LASS_NM  =   45.080400000000000_pr
              SMASS_NM =    0.900000000000000_pr
              VMASS_NM =    1.249838000000000_pr
              if(use_3N_couplings) then
                 use_3N_couplings = .False.
                 override_3N_couplings = .True.
              endif
           Case ('UNE1')
              inforcedname='UNE1'
              use_Namelist=.False.
              ! kind of the functional
              use_INM  = .True.
              ! Surface coefficients
              CrDr(0)  =  -45.135131022237300_pr
              CrDr(1)  = -145.382167908057000_pr
              CpV0(0)  = -186.065399575124000_pr
              CpV0(1)  = -206.579593890106000_pr
              CrdJ(0)  =  -74.026333176459900_pr
              CrdJ(1)  =  -35.658261114791700_pr
              CJ(0)    =    0.000000000000000_pr
              CJ(1)    =    0.000000000000000_pr
              CExPar   =    1.000000000000000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.158706769332587_pr
              E_NM     =  -15.800000000000000_pr
              P_NM     =    0.000000000000000_pr
              K_NM     =  220.000000000000000_pr
              ASS_NM   =   28.986789057772100_pr
              LASS_NM  =   40.004790480413600_pr
              SMASS_NM =    0.992423332283364_pr
              VMASS_NM =    1.249838574232270_pr
              if(use_3N_couplings) then
                 use_3N_couplings = .False.
                 override_3N_couplings = .True.
              endif
           Case ('UNE2')
              inforcedname='UNE2'
              use_Namelist=.False.
              ! kind of the functional
              use_INM     = .True.
              use_j2terms = .True.
              ! Surface coefficients
              CrDr(0)  =  -46.831409147060600_pr
              CrDr(1)  = -113.163790795259000_pr
              CpV0(0)  = -208.889001962571000_pr
              CpV0(1)  = -230.329984038628000_pr
              CrdJ(0)  =  -64.308862415783800_pr
              CrdJ(1)  =  -38.650194685135500_pr
              CJ(0)    =  -54.433363597372100_pr
              CJ(1)    =  -65.903031044593800_pr
              CExPar   =    1.000000000000000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.156310622197074_pr
              E_NM     =  -15.800000000000000_pr
              P_NM     =    0.000000000000000_pr
              K_NM     =  239.929568022437000_pr
              ASS_NM   =   29.131006470773700_pr
              LASS_NM  =   40.000000000000000_pr
              SMASS_NM =    1.073763804147980_pr
              VMASS_NM =    1.249838574232270_pr
              if(use_3N_couplings) then
                 use_3N_couplings = .False.
                 override_3N_couplings = .True.
              endif
           Case ('HFB1')
              ! J. Phys. G: Nucl. Part. Phys. 42 (2015) 034024
              inforcedname='HFB1'
              use_Namelist=.False.
              ! kind of the functional
              use_INM  = .True.
              ! Surface coefficients
              CrDr(0)  =  -45.599763_pr
              CrDr(1)  = -143.935074_pr
              CpV0(0)  = -234.380010_pr
              CpV0(1)  = -260.437001_pr
              CrdJ(0)  =  -73.946388_pr
              CrdJ(1)  =  -51.912548_pr
              CJ(0)    =    0.000000_pr
              CJ(1)    =    0.000000_pr
              CExPar   =    1.000000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.156247_pr
              E_NM     =  -15.800000_pr
              P_NM     =    0.000000_pr
              K_NM     =  244.668204_pr
              ASS_NM   =   28.668204_pr
              LASS_NM  =   40.109081_pr
              SMASS_NM =    1.067970_pr
              VMASS_NM =    1.249838574232270_pr
              if(use_3N_couplings) then
                 use_3N_couplings = .False.
                 override_3N_couplings = .True.
              endif
           Case ('N0LO')
              inforcedname='N0LO'
              use_Namelist=.False.
              ! kind of the functional
              use_INM         = .True.
              use_j2terms     = .False.
              use_DME3N_terms = .False.
              DMEorder        = 0
              ! Surface coefficients
              CrDr(0)  =  -67.437_pr
              CrDr(1)  =   21.551_pr
              CpV0(0)  = -241.203_pr
              CpV0(1)  = -252.818_pr
              CrdJ(0)  =  -95.451_pr
              CrdJ(1)  =  -65.906_pr
              CJ(0)    =    0.000_pr
              CJ(1)    =    0.000_pr
              CExPar   =    1.000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.1595_pr
              E_NM     =  -15.9700_pr
              P_NM     =    0.0000_pr
              K_NM     =  229.9000_pr
              ASS_NM   =   32.0000_pr
              LASS_NM  =   45.9600_pr
              SMASS_NM =    1.4400_pr
              VMASS_NM =    1.2500_pr
              if(use_3N_couplings) then
                 use_3N_couplings = .False.
                 override_3N_couplings = .True.
              endif
           Case ('N1LO')
              inforcedname='N1LO'
              use_Namelist=.False.
              ! kind of the functional
              use_INM         = .True.
              use_j2terms     = .False.
              use_DME3N_terms = .False.
              DMEorder        = 1
              ! Surface coefficients
              CrDr(0)  =  -63.996_pr
              CrDr(1)  =    -9.276_pr
              CpV0(0)  = -241.484_pr
              CpV0(1)  = -252.222
              CrdJ(0)  =  -95.463_pr
              CrdJ(1)  =  -60.800_pr
              CJ(0)    =    0.000_pr
              CJ(1)    =    0.000_pr
              CExPar   =    1.000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.1595_pr
              E_NM     =  -15.9700_pr
              P_NM     =    0.0000_pr
              K_NM     =  229.9000_pr
              ASS_NM   =   32.0000_pr
              LASS_NM  =   45.9600_pr
              SMASS_NM =    1.4400_pr
              VMASS_NM =    1.2500_pr
              if(use_3N_couplings) then
                 use_3N_couplings = .False.
                 override_3N_couplings = .True.
              endif
           Case ('N2LO')
              inforcedname='N2LO'
              use_Namelist=.False.
              ! kind of the functional
              use_INM         = .True.
              use_j2terms     = .False.
              use_DME3N_terms = .True.
              DMEorder        = 2
              ! Surface coefficients
              CrDr(0)  = -197.132_pr
              CrDr(1)  =  -12.503_pr
              CpV0(0)  = -272.164_pr
              CpV0(1)  = -193.188_pr
              CrdJ(0)  = -193.188_pr
              CrdJ(1)  =   37.790_pr
              CJ(0)    =    0.000_pr
              CJ(1)    =    0.000_pr
              CExPar   =    1.000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.1595_pr
              E_NM     =  -15.9700_pr
              P_NM     =    0.0000_pr
              K_NM     =  229.9000_pr
              ASS_NM   =   32.0000_pr
              LASS_NM  =   45.9600_pr
              SMASS_NM =    1.4400_pr
              VMASS_NM =    1.2500_pr
              if(use_3N_couplings) then
                 use_3N_couplings = .False.
                 override_3N_couplings = .True.
              endif
           !---------------------------------------------------------------------
           !NEDF Nuclear energy density functional
           !---------------------------------------------------------------------
           Case ('NEDF')
              inforcedname='NEDF'
              use_Namelist=.False.
              ! kind of functional
              is_nedf = .True.
              hb0_charge_dependent = .True.
              a_NEDF=[   0._pr         ,31.4703220539_pr, -64.0583220539_pr]
              b_NEDF=[-674.797522134_pr,59.2294137863_pr, 500.1281083477_pr]
              c_NEDF=[ 804.984434780_pr, 0._pr          ,-695.8644347800_pr]
              eta_NEDF = 90.300256558980_pr
              W0_NEDF = 76.214174675_pr
              ! Coupling constants: ph channel
              Crho  = zero
              CDrho = zero
              Ctau  = zero
              CrDr  = zero
              CrdJ  = zero
              CJ    = zero
              Cnrho = eta_NEDF/two
              Cjdr(0)  = W0_NEDF; Cjdr(1)  = zero
              ! Coupling constants: pp channel
              CpV0  = -200._pr
              CpV1  = zero
              ! Various
              hbzero  = 20.735527840105213_pr
              hbzeron = 20.721246548757218_pr
              hbzerop = 20.749809131453205_pr
              e2charg = 1.4399644_pr
           Case ('SeaLL1')
              inforcedname='SeaLL1'
              use_Namelist=.False.
              ! kind of functional
              is_nedf = .True.
              hb0_charge_dependent = .True.
              a_NEDF=[   0._pr         ,  64.2474102072_pr, -96.8354102072_pr]
              b_NEDF=[-684.524043779_pr, 119.8621469590_pr, 449.221896820_pr]
              c_NEDF=[ 827.262878410_pr,-256.4927039210_pr,-461.650174489_pr]
              eta_NEDF = 81.3917529003_pr
              W0_NEDF  = 73.5210618422_pr
              ! Coupling constants: ph channel
              Crho  = zero
              CDrho = zero
              Ctau  = zero
              CrDr  = zero
              CrdJ  = zero
              CJ    = zero
              Cnrho = eta_NEDF/two
              Cjdr(0)  = W0_NEDF; Cjdr(1)  = zero
              ! Coupling constants: pp channel
              CpV0  = -200._pr
              CpV1  = zero
              ! Various
              hbzero  = 20.735527840105213_pr
              hbzeron = 20.721246548757218_pr
              hbzerop = 20.749809131453205_pr
              e2charg = 1.4399644_pr
           Case ('NEDF1')
              inforcedname='NEDF1'
              use_Namelist=.False.
              ! kind of functional
              is_nedf = .True.
              hb0_charge_dependent = .True.
              a_NEDF=[   0._pr         ,30.1396291921_pr, -62.7276291921_pr]
              b_NEDF=[-672.784776452_pr,56.7249539338_pr, 500.6198225180_pr]
              c_NEDF=[ 802.203802846_pr, 0._pr          ,-693.0838028460_pr]
              eta_NEDF = 105.629528244_pr
              W0_NEDF = 76.214174675_pr
              ! Coupling constants: ph channel
              Crho  = zero
              CDrho = zero
              Ctau  = zero
              CrDr  = zero
              CrdJ  = zero
              CJ    = zero
              Cnrho = eta_NEDF/two
              Cjdr(0)  = W0_NEDF; Cjdr(1)  = zero
              ! Coupling constants: pp channel
              CpV0  = -200._pr
              CpV1  = zero
              ! Various
              hbzero  = 20.735527840105213_pr
              hbzeron = 20.721246548757218_pr
              hbzerop = 20.749809131453205_pr
              e2charg = 1.4399644_pr
           Case ('NEDF2')
              inforcedname='NEDF2'
              use_Namelist=.False.
              ! kind of functional
              is_nedf = .True.
              hb0_charge_dependent = .True.
              a_NEDF=[   0._pr         ,30.3723779104_pr, -62.9603779104_pr]
              b_NEDF=[-672.212514042_pr,57.1630037930_pr, 499.6095102490_pr]
              c_NEDF=[ 801.413215568_pr, 0._pr          ,-692.2932155680_pr]
              eta_NEDF = 96.1410827861_pr
              W0_NEDF = 76.214174675_pr
              ! Coupling constants: ph channel
              Crho  = zero
              CDrho = zero
              Ctau  = zero
              CrDr  = zero
              CrdJ  = zero
              CJ    = zero
              Cnrho = eta_NEDF/two
              Cjdr(0)  = W0_NEDF; Cjdr(1)  = zero
              ! Coupling constants: pp channel
              CpV0  = -200._pr
              CpV1  = zero
              ! Various
              hbzero  = 20.735527840105213_pr
              hbzeron = 20.721246548757218_pr
              hbzerop = 20.749809131453205_pr
              e2charg = 1.4399644_pr
           Case ('NEDF3')
              inforcedname='NEDF3'
              use_Namelist=.False.
              ! kind of functional
              is_nedf = .True.
              hb0_charge_dependent = .True.
              a_NEDF=[   0._pr         ,30.3259301650_pr, -62.913930165_pr]
              b_NEDF=[-671.889100224_pr,57.0755857892_pr, 499.373514435_pr]
              c_NEDF=[ 800.966415558_pr, 0._pr          ,-691.846415558_pr]
              eta_NEDF = 91.896020403_pr
              W0_NEDF = 76.214174675_pr
              ! Coupling constants: ph channel
              Crho  = zero
              CDrho = zero
              Ctau  = zero
              CrDr  = zero
              CrdJ  = zero
              CJ    = zero
              Cnrho = eta_NEDF/two
              Cjdr(0)  = W0_NEDF; Cjdr(1)  = zero
              ! Coupling constants: pp channel
              CpV0  = -200._pr
              CpV1  = zero
              ! Various
              hbzero  = 20.735527840105213_pr
              hbzeron = 20.721246548757218_pr
              hbzerop = 20.749809131453205_pr
              e2charg = 1.4399644_pr
           Case ('NEDF4')
              inforcedname='NEDF4'
              use_Namelist=.False.
              ! kind of functional
              is_nedf = .True.
              hb0_charge_dependent = .True.
              a_NEDF=[   0._pr         ,31.5030085475_pr, -64.0910085475_pr]
              b_NEDF=[-672.625018435_pr,59.2909321226_pr, 497.8940863130_pr]
              c_NEDF=[ 801.983095259_pr, 0._pr          ,-692.8630952590_pr]
              eta_NEDF = 87.2919112779_pr
              W0_NEDF = 76.214174675_pr
              ! Coupling constants: ph channel
              Crho  = zero
              CDrho = zero
              Ctau  = zero
              CrDr  = zero
              CrdJ  = zero
              CJ    = zero
              Cnrho = eta_NEDF/two
              Cjdr(0)  = W0_NEDF; Cjdr(1)  = zero
              ! Coupling constants: pp channel
              CpV0  = -200._pr
              CpV1  = zero
              ! Various
              hbzero  = 20.735527840105213_pr
              hbzeron = 20.721246548757218_pr
              hbzerop = 20.749809131453205_pr
              e2charg = 1.4399644_pr
           Case ('NEDF5')
              inforcedname='NEDF5'
              use_Namelist=.False.
              ! kind of functional
              is_nedf = .True.
              hb0_charge_dependent = .True.
              a_NEDF=[   0.0_pr        , 43.8215204792_pr,  -76.4095204792_pr]
              b_NEDF=[-669.239299724_pr, 80.7199308819_pr,  473.0793688420_pr]
              c_NEDF=[ 797.305684868_pr,-92.7644849555_pr, -595.4211999120_pr]
              eta_NEDF = 82.95333556289553_pr
              W0_NEDF = 76.214174675_pr
              ! Coupling constants: ph channel
              Crho  = zero
              CDrho = zero
              Ctau  = zero
              CrDr  = zero
              CrdJ  = zero
              CJ    = zero
              Cnrho = eta_NEDF/two
              Cjdr(0)  = W0_NEDF; Cjdr(1)  = zero
              ! Coupling constants: pp channel
              CpV0  = -200._pr
              CpV1  = zero
              ! Various
              hbzero  = 20.735527840105213_pr
              hbzeron = 20.721246548757218_pr
              hbzerop = 20.749809131453205_pr
              e2charg = 1.4399644_pr
           !---------------------------------------------------------------------
           !DME To Finite range chiral potentials.
           !---------------------------------------------------------------------
           Case ('DME_LO')
              inforcedname='DME_LO'
              use_Namelist=.False.
              ! kind of functional
              use_INM      = .True.
              use_j2terms  = .True.
              force_is_dme = .true.
              alfa_dme     = 1._pr
              ! Surface coefficients
              CrDr(0)  =  -36.184280066031600_pr
              CrDr(1)  =  -70.270337907872700_pr
              CpV0(0)  = -194.765957487293000_pr
              CpV0(1)  = -227.129830227349000_pr
              CrdJ(0)  =  -62.015368161964000_pr
              CrdJ(1)  =  -81.261507317600400_pr
              CJ(0)    = -101.894527100110000_pr
              CJ(1)    =   34.453836064809500_pr
              CExPar   =    1.000000000000000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.154430976350843_pr
              E_NM     =  -15.802473677471800_pr
              P_NM     =    0.000000000000000_pr
              K_NM     =  258.653617357751000_pr
              ASS_NM   =   30.057808846781400_pr
              LASS_NM  =   41.957729056715900_pr
              SMASS_NM =    0.976297468233889_pr
              VMASS_NM =    1.249838574232270_pr
              call set_dme_couplings(fname)
              if(use_3N_couplings) then
                 use_3N_couplings = .False.
                 override_3N_couplings = .True.
              endif
           Case ('DME_NLO')
              inforcedname='DME_NLO'
              use_Namelist=.False.
              ! kind of functional
              use_INM      = .True.
              use_j2terms  = .True.
              finite_range = .True.
              force_is_dme = .True.
              alfa_dme     = 1._pr
              ! Surface coefficients
              CrDr(0)  =  -52.644745897289800_pr
              CrDr(1)  =  -61.145395220329500_pr
              CpV0(0)  = -163.965944945033000_pr
              CpV0(1)  = -189.080575493328000_pr
              CrdJ(0)  =  -62.193662887152700_pr
              CrdJ(1)  = -104.361601668738000_pr
              CJ(0)    =  -84.884220406100000_pr
              CJ(1)    =   31.074843204318200_pr
              CExPar   =    1.000000000000000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.154994136163034_pr
              E_NM     =  -15.800000000000000_pr
              P_NM     =    0.000000000000000_pr
              K_NM     =  254.956393939499000_pr
              ASS_NM   =   30.520106114644500_pr
              LASS_NM  =   42.994668332820700_pr
              SMASS_NM =    0.900000000000000_pr
              VMASS_NM =    1.249838574232270_pr
              ! W, B, H, M parameters of the finite-range two-body
              Call gogny_force(fname)
              call set_dme_couplings(fname)
              if(use_3N_couplings) then
                 use_3N_couplings = .False.
                 override_3N_couplings = .True.
              endif
           Case ('DME_N2LO')
              inforcedname='DME_N2LO'
              use_Namelist=.False.
              ! kind of functional
              use_INM          = .True.
              use_j2terms      = .True.
              finite_range     = .True.
              force_is_dme     = .True.
              ! use_3N_couplings = .True.
              alfa_dme         = 1._pr
              ! Surface coefficients
              CrDr(0)  =    6.989194393828900_pr
              CrDr(1)  =  -65.905177885039700_pr
              CpV0(0)  = -164.224186429086000_pr
              CpV0(1)  = -191.103480033817000_pr
              CrdJ(0)  =  -67.167603212928400_pr
              CrdJ(1)  =  -71.254770672344700_pr
              CJ(0)    = -100.486880561974000_pr
              CJ(1)    =  -41.598127580615400_pr
              CExPar   =    1.000000000000000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.158351677657629_pr
              E_NM     =  -15.835256707635300_pr
              P_NM     =    0.000000000000000_pr
              K_NM     =  221.241255926659000_pr
              ASS_NM   =   29.755365769845100_pr
              LASS_NM  =   40.000000000000000_pr
              SMASS_NM =    0.904818485586906_pr
              VMASS_NM =    1.249838574232270_pr
              ! W, B, H, M parameters of the finite-range two-body
              Call gogny_force(fname)
              call set_dme_couplings(fname)
              if(use_3N_couplings) then
                 ! Surface coefficients
                 CrDr(0)  =  -29.129313562207500_pr
                 CrDr(1)  =  -51.641409431444500_pr
                 CpV0(0)  = -164.145096920900000_pr
                 CpV0(1)  = -190.066440123089000_pr
                 CrdJ(0)  =  -72.837039627045200_pr
                 CrdJ(1)  =  -66.790139645005400_pr
                 CJ(0)    = -104.900812174163000_pr
                 CJ(1)    =  -11.068293023710000_pr
                 ! Associated INM parameters
                 RHO_NM   =    0.153029494019512_pr
                 E_NM     =  -15.813341547417200_pr
                 K_NM     =  250.355901971974000_pr
                 ASS_NM   =   29.264035166590900_pr
                 LASS_NM  =   40.250000000000000_pr
                 SMASS_NM =    0.900000000000000_pr
                 call set_3n_couplings(fname)
              endif
           Case ('DME_NLOD')
              inforcedname='DME_NLOD'
              use_Namelist=.False.
              ! kind of functional
              use_INM      = .True.
              use_j2terms  = .True.
              finite_range = .True.
              force_is_dme = .True.
              alfa_dme     = 1._pr
              ! Surface coefficients
              CrDr(0)  =  -13.786950951288400_pr
              CrDr(1)  =  -68.901631461819500_pr
              CpV0(0)  = -163.390139040488000_pr
              CpV0(1)  = -188.702817302331000_pr
              CrdJ(0)  =  -94.918899286392000_pr
              CrdJ(1)  =  -40.334637716943500_pr
              CJ(0)    =  -18.714467834162700_pr
              CJ(1)    =   39.632013138712000_pr
              CExPar   =    1.000000000000000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.158566222958099_pr
              E_NM     =  -15.861668539861900_pr
              P_NM     =    0.000000000000000_pr
              K_NM     =  223.030410024355000_pr
              ASS_NM   =   30.504161781633500_pr
              LASS_NM  =   44.207693483573800_pr
              SMASS_NM =    0.900000000000000_pr
              VMASS_NM =    1.249838000000000_pr
              ! W, B, H, M parameters of the finite-range two-body
              Call gogny_force(fname)
              call set_dme_couplings(fname)
              if(use_3N_couplings) then
                 ! Surface coefficients
                 CrDr(0)  =  -33.231209431539500_pr
                 CrDr(1)  =  -48.579735677238000_pr
                 CpV0(0)  = -164.720240152312000_pr
                 CpV0(1)  = -189.696969418118000_pr
                 CrdJ(0)  =  -75.885499961593400_pr
                 CrdJ(1)  =  -46.261173248881100_pr
                 CJ(0)    =  -86.849555488362700_pr
                 CJ(1)    =   14.931749554909900_pr
                 ! Associated INM parameters
                 RHO_NM   =    0.153928912056897_pr
                 E_NM     =  -15.813542795736700_pr
                 K_NM     =  250.013710221248000_pr
                 ASS_NM   =   29.692138213755300_pr
                 LASS_NM  =   40.000000000000000_pr
                 SMASS_NM =    0.900000000000000_pr
                 call set_3n_couplings(fname)
              endif
           Case ('DME_N2LOD')
              inforcedname='DME_N2LOD'
              use_Namelist=.False.
              ! kind of functional
              use_INM      = .True.
              use_j2terms  = .True.
              finite_range = .True.
              force_is_dme = .True.
              alfa_dme     = 1._pr
              ! Surface coefficients
              CrDr(0)  =    1.043960766688410_pr
              CrDr(1)  =  -69.111057901842400_pr
              CpV0(0)  = -163.191576998704000_pr
              CpV0(1)  = -189.026411466385000_pr
              CrdJ(0)  =  -64.064577174344700_pr
              CrdJ(1)  =  -37.698931086984300_pr
              CJ(0)    = -112.115747266403000_pr
              CJ(1)    =  -10.316426355326800_pr
              CExPar   =    1.000000000000000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.155590485949568_pr
              E_NM     =  -15.838536372452200_pr
              P_NM     =    0.000000000000000_pr
              K_NM     =  248.205818148798000_pr
              ASS_NM   =   29.698295354209200_pr
              LASS_NM  =   41.941194986544700_pr
              SMASS_NM =    0.900000000000000_pr
              VMASS_NM =    1.249838000000000_pr
              ! W, B, H, M parameters of the finite-range two-body
              Call gogny_force(fname)
              call set_dme_couplings(fname)
              if(use_3N_couplings) then
                 ! Surface coefficients
                 CrDr(0)  =  -35.368179901851600_pr
                 CrDr(1)  =   -4.157311061417660_pr
                 CpV0(0)  = -163.577306416090000_pr
                 CpV0(1)  = -190.270201375991000_pr
                 CrdJ(0)  =  -69.672413603996300_pr
                 CrdJ(1)  =  -64.976688763511900_pr
                 CJ(0)    = -115.910795064510000_pr
                 CJ(1)    =  -20.766273634196300_pr
                 ! Associated INM parameters
                 RHO_NM   =    0.152724608962683_pr
                 E_NM     =  -15.841688141661500_pr
                 K_NM     =  259.242325078061000_pr
                 ASS_NM   =   30.404007524673000_pr
                 LASS_NM  =   40.000000000000000_pr
                 SMASS_NM =    0.900038000058470_pr
                 call set_3n_couplings(fname)
              endif
           Case ('REG_LO')
              inforcedname='REG_LO'
              use_Namelist=.False.
              ! kind of functional
              use_INM      = .True.
              use_j2terms  = .True.
              force_is_dme = .true.
              alfa_dme     = 1._pr
              ! Surface coefficients
              CrDr(0)  =  -43.685848651174600_pr
              CrDr(1)  =  -85.136292247205000_pr
              CpV0(0)  = -203.294067105786000_pr
              CpV0(1)  = -231.625817447383000_pr
              CrdJ(0)  =  -62.384142876249500_pr
              CrdJ(1)  =  -78.149705107064500_pr
              CJ(0)    =  -67.806949472637900_pr
              CJ(1)    =   -2.479721269022790_pr
              CExPar   =    1.000000000000000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.155680484003159_pr
              E_NM     =  -15.800000000000000_pr
              P_NM     =    0.000000000000000_pr
              K_NM     =  246.542921116394000_pr
              ASS_NM   =   29.748490283602300_pr
              LASS_NM  =   40.000000000000000_pr
              SMASS_NM =    1.048977829021120_pr
              VMASS_NM =    1.249838574232270_pr
              call set_dme_couplings(fname)
              if(use_3N_couplings) then
                 use_3N_couplings = .False.
                 override_3N_couplings = .True.
              endif
           Case ('REG_NLO')
              inforcedname='REG_NLO'
              use_Namelist=.False.
              ! kind of functional
              use_INM      = .True.
              use_j2terms  = .True.
              finite_range = .True.
              force_is_dme = .True.
              alfa_dme     = 1._pr
              ! Surface coefficients
              CrDr(0)  =  -46.124117664456100_pr
              CrDr(1)  =  -76.114695060911800_pr
              CpV0(0)  = -172.431908169702000_pr
              CpV0(1)  = -195.222250131617000_pr
              CrdJ(0)  =  -68.331224308185300_pr
              CrdJ(1)  =  -89.319396586598400_pr
              CJ(0)    =  -33.797779281639400_pr
              CJ(1)    =   15.787970445575300_pr
              CExPar   =    1.000000000000000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.156195283941728_pr
              E_NM     =  -15.800000000000000_pr
              P_NM     =    0.000000000000000_pr
              K_NM     =  244.834911975064000_pr
              ASS_NM   =   29.731176752598100_pr
              LASS_NM  =   40.000000000000000_pr
              SMASS_NM =    0.990376134955352_pr
              VMASS_NM =    1.249838574232270_pr
              ! W, B, H, M parameters of the finite-range two-body
              Call gogny_force(fname)
              call set_dme_couplings(fname)
              if(use_3N_couplings) then
                 use_3N_couplings = .False.
                 override_3N_couplings = .True.
              endif
           Case ('REG_N2LO')
              inforcedname='REG_N2LO'
              use_Namelist=.False.
              ! kind of functional
              use_INM          = .True.
              use_j2terms      = .True.
              finite_range     = .True.
              force_is_dme     = .True.
              ! use_3N_couplings = .True.
              alfa_dme         = 1._pr
              ! Surface coefficients
              CrDr(0)  =  -37.431098840766600_pr
              CrDr(1)  =  -99.824966690444500_pr
              CpV0(0)  = -179.103001387633000_pr
              CpV0(1)  = -199.707935811866000_pr
              CrdJ(0)  =  -70.204086093431000_pr
              CrdJ(1)  =  -44.399705107064500_pr
              CJ(0)    =  -32.806949472637900_pr
              CJ(1)    =   -5.849671858453230_pr
              CExPar   =    1.000000000000000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.155249617997318_pr
              E_NM     =  -15.809456396888300_pr
              P_NM     =    0.000000000000000_pr
              K_NM     =  259.456169802079000_pr
              ASS_NM   =   28.803483548161700_pr
              LASS_NM  =   40.000000000000000_pr
              SMASS_NM =    1.046159468832010_pr
              VMASS_NM =    1.249838574232270_pr
              ! W, B, H, M parameters of the finite-range two-body
              Call gogny_force(fname)
              call set_dme_couplings(fname)
              if(use_3N_couplings) then
                 ! Surface coefficients
                 CrDr(0)  =  -40.114769373261700_pr
                 CrDr(1)  =  -73.002980697872600_pr
                 CpV0(0)  = -167.161104104339000_pr
                 CpV0(1)  = -191.810525198722000_pr
                 CrdJ(0)  =  -69.616816229645500_pr
                 CrdJ(1)  =  -38.304733250755400_pr
                 CJ(0)    =  -31.546085430728900_pr
                 CJ(1)    =   10.669432990558100_pr
                 ! Associated INM parameters
                 RHO_NM   =    0.154958309238237_pr
                 E_NM     =  -15.800000000000000_pr
                 K_NM     =  260.000000000000000_pr
                 ASS_NM   =   29.443425158154700_pr
                 LASS_NM  =   47.521039614767500_pr
                 SMASS_NM =    0.941084935108090_pr
                 call set_3n_couplings(fname)
              endif
           Case ('REG_NLOD')
              inforcedname='REG_NLOD'
              use_Namelist=.False.
              ! kind of functional
              use_INM      = .True.
              use_j2terms  = .True.
              finite_range = .True.
              force_is_dme = .True.
              alfa_dme     = 1._pr
              ! Surface coefficients
              CrDr(0)  =  -40.669080992703000_pr
              CrDr(1)  =  -86.175846888767000_pr
              CpV0(0)  = -175.303469310476000_pr
              CpV0(1)  = -197.648077406864000_pr
              CrdJ(0)  =  -72.419327948040800_pr
              CrdJ(1)  =  -66.848580025574400_pr
              CJ(0)    =  -23.706800159766100_pr
              CJ(1)    =    9.676548018459500_pr
              CExPar   =    1.000000000000000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.155635914985647_pr
              E_NM     =  -15.805150888847100_pr
              P_NM     =    0.000000000000000_pr
              K_NM     =  248.907039879481000_pr
              ASS_NM   =   29.340303195075200_pr
              LASS_NM  =   40.000000000000000_pr
              SMASS_NM =    1.023011075927570_pr
              VMASS_NM =    1.249838000000000_pr
              ! W, B, H, M parameters of the finite-range two-body
              Call gogny_force(fname)
              call set_dme_couplings(fname)
              if(use_3N_couplings) then
                 ! Surface coefficients
                 CrDr(0)  =  -40.169775516345400_pr
                 CrDr(1)  =  -63.392885785677800_pr
                 CpV0(0)  = -161.884301248489000_pr
                 CpV0(1)  = -186.906640221911000_pr
                 CrdJ(0)  =  -73.569208038118200_pr
                 CrdJ(1)  =  -61.606526569648100_pr
                 CJ(0)    =  -11.596288618941400_pr
                 CJ(1)    =    2.824960117986440_pr
                 ! Associated INM parameters
                 RHO_NM   =    0.155113516190737_pr
                 E_NM     =  -15.806936450605300_pr
                 K_NM     =  252.768695937858000_pr
                 ASS_NM   =   29.850582455590600_pr
                 LASS_NM  =   41.935497043472100_pr
                 SMASS_NM =    0.901074681502770_pr
                 call set_3n_couplings(fname)
              endif
           Case ('REG_N2LOD')
              inforcedname='REG_N2LOD'
              use_Namelist=.False.
              ! kind of functional
              use_INM      = .True.
              use_j2terms  = .True.
              finite_range = .True.
              force_is_dme = .True.
              alfa_dme     = 1._pr
              ! Surface coefficients
              CrDr(0)  =  -34.093328268374800_pr
              CrDr(1)  =  -73.080644860918500_pr
              CpV0(0)  = -162.206218600365000_pr
              CpV0(1)  = -187.472113009341000_pr
              CrdJ(0)  =  -66.681658875513100_pr
              CrdJ(1)  =  -49.706084478615000_pr
              CJ(0)    =  -34.848339320248600_pr
              CJ(1)    =  -12.501086595857100_pr
              CExPar   =    1.000000000000000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.154766929878851_pr
              E_NM     =  -15.800000000000000_pr
              P_NM     =    0.000000000000000_pr
              K_NM     =  260.000000000000000_pr
              ASS_NM   =   29.286557588572900_pr
              LASS_NM  =   40.000000000000000_pr
              SMASS_NM =    0.900000000000000_pr
              VMASS_NM =    1.249838000000000_pr
              ! W, B, H, M parameters of the finite-range two-body
              Call gogny_force(fname)
              call set_dme_couplings(fname)
              if(use_3N_couplings) then
                 ! Surface coefficients
                 CrDr(0)  =  -40.366672057367400_pr
                 CrDr(1)  =  -61.190787119369200_pr
                 CpV0(0)  = -161.756975562431000_pr
                 CpV0(1)  = -186.799578594996000_pr
                 CrdJ(0)  =  -61.137702817074100_pr
                 CrdJ(1)  =  -48.233026271880900_pr
                 CJ(0)    =  -52.719389559626500_pr
                 CJ(1)    =    9.736563542572240_pr
                 ! Associated INM parameters
                 RHO_NM   =    0.155186604293239_pr
                 E_NM     =  -15.800000000000000_pr
                 K_NM     =  254.124497711116000_pr
                 ASS_NM   =   29.490918546324300_pr
                 LASS_NM  =   40.000000000000000_pr
                 SMASS_NM =    0.900000000000000_pr
                 call set_3n_couplings(fname)
              endif
           Case default
              inforcedname=fname
              use_Namelist=.True.
           End Select FUNCTIONAL
       End If
       !---------------------------------------------------------------------
       ! Exit loop condition
       !---------------------------------------------------------------------
       If(.Not.use_Namelist) Exit
       !---------------------------------------------------------------------
       ! Read namelists
       !---------------------------------------------------------------------
       Open(UNIT=lnamelist,file=filename,status='OLD',delim='APOSTROPHE')
       Read(UNIT=lnamelist,NML=UNEDF_NAMELIST,iostat=ios)
       If(ios.Ne.0) Then
          ! WRong entry within UNEDF_NAMELIST.DAT file
          Write(*,'(1X,/,A)') 'ATTENTION: WRONG INPUT!'
          Write(*,*) 'THE INPUT DATA WITH LABEL FUNCTIONALNAME=''',Trim(INFORCEDNAME),''''
          Write(*,*) 'INSIDE THE UNEDF_NAMELIST.DAT FILE IS WRONG.'
          Write(*,*) 'PLEASE CORRECT AND TRY AGAIN!'
          Stop 'PROGRAM STOP IN read_UNEDF_NAMELIST'
       End If
       Close(lnamelist)
       If(Trim(FunctionalName).Eq.Trim(inforcedname)) Exit
    End Do
    !---------------------------------------------------------------------
    ! See what the namelists modified
    !---------------------------------------------------------------------
    INFORCED_FUNCTIONAL: Select Case (Trim(inforcedname))
    Case ("FORCE")
       FunctionalName='FORCE'
    Case ("UNE0")
       FunctionalName='UNE0'
    Case ("UNE1")
       FunctionalName='UNE1'
    Case ("UNE2")
       FunctionalName='UNE2'
    Case ("HFB1")
       FunctionalName='HFB1'
    Case ("N0LO")
       FunctionalName='N0LO'
    Case ("N1LO")
       FunctionalName='N1LO'
    Case ("N2LO")
       FunctionalName='N2LO'
    Case ("NEDF")
       FunctionalName='NEDF'
    Case ("SeaLL1")
       FunctionalName='SeaLL1'
    Case ("NEDF1")
       FunctionalName='NEDF1'
    Case ("NEDF2")
       FunctionalName='NEDF2'
    Case ("NEDF3")
       FunctionalName='NEDF3'
    Case ("NEDF4")
       FunctionalName='NEDF4'
    Case ("NEDF5")
       FunctionalName='NEDF5'
    Case ("READ")
       FunctionalName='READ'
    Case ("RAND")
       FunctionalName='RAND'
    Case ("FITS")
       FunctionalName='FITS'
    Case ("DME_LO")
       FunctionalName='DME_LO'
    Case ("DME_NLO")
       FunctionalName='DME_NLO'
    Case ("DME_N2LO")
       FunctionalName='DME_N2LO'
    Case ("DME_NLOD")
       FunctionalName='DME_NLOD'
    Case ("DME_N2LOD")
       FunctionalName='DME_N2LOD'
    Case ("REG_LO")
       FunctionalName='REG_LO'
    Case ("REG_NLO")
       FunctionalName='REG_NLO'
    Case ("REG_N2LO")
       FunctionalName='REG_N2LO'
    Case ("REG_NLOD")
       FunctionalName='REG_NLOD'
    Case ("REG_N2LOD")
       FunctionalName='REG_N2LOD'
    Case default
       ! Missing entry within hfbtho_NAMELIST.dat file
       If(Trim(FunctionalName).Ne.Trim(inforcedname)) Then
          Write(*,'(1X,/,A)') 'ATTENTION: MISSING INPUT!'
          Write(*,*) 'THE INPUT DATA WITH LABEL FUNCTIONALNAME=''',Trim(INFORCEDNAME),''''
          Write(*,*) 'IS MISSING INSIDE THE UNEDF_NAMELIST.DAT FILE.'
          Write(*,*) 'PLEASE CORECT AND TRY AGAIN!'
          Stop 'PROGRAM STOP IN SET_FUNCTIONAL_PARAMETERS'
       End If
    End Select INFORCED_FUNCTIONAL
    !
  End Subroutine read_UNEDF_NAMELIST
  !=======================================================================
  !> Set up parameters of the Gogny force
  !=======================================================================
  Subroutine gogny_force(fname)
    Implicit None
    Character (30), Intent(inout) :: fname
    !
    INTERACTION: Select Case (Trim(fname))
    !---------------------------------------------------------------------
    ! PRC 21, 1568 (1980)
    !---------------------------------------------------------------------
    Case ('D1')
       n_g = 2
       if(.not.allocated(mu_g)) &
            allocate(mu_g(1:n_g),W_g(1:n_g),B_g(1:n_g),H_g(1:n_g),M_g(1:n_g))
       mu_g= [ 0.7_pr, 1.2_pr ]
       W_g = [ -402.40_pr, -21.30_pr ]
       B_g = [ -100.00_pr, -11.77_pr ]
       H_g = [ -496.20_pr,  37.27_pr ]
       M_g = [  -23.56_pr, -68.81_pr ]
    !---------------------------------------------------------------------
    ! CPC 63, 365 (1991)
    !---------------------------------------------------------------------
    Case ('D1S')
       n_g = 2
       if(.not.allocated(mu_g)) &
            allocate(mu_g(1:n_g),W_g(1:n_g),B_g(1:n_g),H_g(1:n_g),M_g(1:n_g))
       mu_g= [ 0.7_pr, 1.2_pr ]
       W_g = [ -1720.30_pr,  103.64_pr ]
       B_g = [  1300.00_pr, -163.48_pr ]
       H_g = [ -1813.53_pr,  162.81_pr ]
       M_g = [  1397.60_pr, -223.93_pr ]
       !W_g = [ -1720.30_pr,  103.639_pr ]
       !B_g = [  1300.00_pr, -163.483_pr ]
       !H_g = [ -1813.53_pr,  162.812_pr ]
       !M_g = [  1397.60_pr, -223.934_pr ]
    !---------------------------------------------------------------------
    ! PRC 21, 1568 (1980)
    !---------------------------------------------------------------------
    Case ('D1p')
       n_g = 2
       if(.not.allocated(mu_g)) &
            allocate(mu_g(1:n_g),W_g(1:n_g),B_g(1:n_g),H_g(1:n_g),M_g(1:n_g))
       mu_g= [ 0.7_pr, 1.2_pr ]
       W_g = [ -402.40_pr, -21.30_pr ]
       B_g = [ -100.00_pr, -11.77_pr ]
       H_g = [ -496.20_pr,  37.27_pr ]
       M_g = [  -23.56_pr, -68.81_pr ]
    !---------------------------------------------------------------------
    !
    !---------------------------------------------------------------------
    Case ('D1N')
       n_g = 2
       if(.not.allocated(mu_g)) &
            allocate(mu_g(1:n_g),W_g(1:n_g),B_g(1:n_g),H_g(1:n_g),M_g(1:n_g))
       mu_g= [ 0.8_pr, 1.2_pr ]
       W_g = [ -2047.61_pr,  293.02_pr ]
       B_g = [  1700.00_pr, -300.78_pr ]
       H_g = [ -2414.93_pr,  414.59_pr ]
       M_g = [  1519.35_pr, -316.84_pr ]
    !---------------------------------------------------------------------
    !
    !---------------------------------------------------------------------
    Case ('DME_NLO')
       n_g = 5
       if(.not.allocated(mu_g)) &
            allocate(mu_g(1:n_g),W_g(1:n_g),B_g(1:n_g),H_g(1:n_g),M_g(1:n_g))
       ! Regulator: (1-exp(-(r/1.0)**2))**6 !arXiv:nucl-th/0703087 couplings
       mu_g= [ 1.74125407_pr,  1.08327461_pr,    0.59032535_pr,    0.56379030_pr, 0.54333625_pr]
       W_g = [ 0.54005090_pr,  6.91550331_pr, -344.00093544_pr,  659.22825793_pr, 0._pr]*alfa_dme
       B_g = zero
       H_g = [-1.08010179_pr,-13.83100661_pr,  688.00187089_pr, -1318.45651585_pr,0._pr]*alfa_dme
       M_g = zero
       W_g(n_g) = -sum(W_g(1:n_g-1))
       H_g(n_g) = -sum(H_g(1:n_g-1))
    !---------------------------------------------------------------------
    !
    !---------------------------------------------------------------------
    Case ('DME_N2LO')
       n_g = 5
       if(.not.allocated(mu_g)) &
            allocate(mu_g(1:n_g),W_g(1:n_g),B_g(1:n_g),H_g(1:n_g),M_g(1:n_g))
       ! Regulator: (1-exp(-(r/1.0)**2))**6 !arXiv:nucl-th/0703087 couplings
       mu_g= [  1.38358765_pr,   0.88231221_pr,   0.68955965_pr,     0.46579342_pr, 0.46000139_pr]
       W_g = [-11.20121649_pr,-142.22940980_pr, 274.17013591_pr, -1552.39921847_pr, 0._pr]*alfa_dme
       B_g = zero
       H_g = [ -5.67965601_pr, -15.73992290_pr,  47.48299561_pr,  -417.94912176_pr, 0._pr]*alfa_dme
       M_g = zero
       W_g(n_g) = -sum(W_g(1:n_g-1))
       H_g(n_g) = -sum(H_g(1:n_g-1))
    !---------------------------------------------------------------------
    !
    !---------------------------------------------------------------------
    Case ('DME_NLOD')
       n_g = 5
       if(.not.allocated(mu_g)) &
            allocate(mu_g(1:n_g),W_g(1:n_g),B_g(1:n_g),H_g(1:n_g),M_g(1:n_g))
       ! Regulator: (1-exp(-(r/1.0)**2))**6 !arXiv:nucl-th/0703087 couplings
       mu_g= [ 1.40150660_pr,   0.91011501_pr,   0.66261233_pr,     0.46234689_pr, 0.45624925_pr]
       W_g = [-7.25024547_pr, -97.37816779_pr, 219.61898151_pr, -1370.65240796_pr, 0._pr]*alfa_dme
       B_g = zero
       H_g = [-6.69300176_pr, -39.24095540_pr, 103.59850575_pr,  -740.89522639_pr, 0._pr]*alfa_dme
       M_g = zero
       W_g(n_g) = -sum(W_g(1:n_g-1))
       H_g(n_g) = -sum(H_g(1:n_g-1))
    !---------------------------------------------------------------------
    !
    !---------------------------------------------------------------------
    Case ('DME_N2LOD')
       n_g = 5
       if(.not.allocated(mu_g)) &
            allocate(mu_g(1:n_g),W_g(1:n_g),B_g(1:n_g),H_g(1:n_g),M_g(1:n_g))
       ! Regulator: (1-exp(-(r/1.0)**2))**6 !arXiv:nucl-th/0703087 couplings
       mu_g= [ 1.41924549_pr,    0.93739092_pr,   0.62201921_pr,     0.47289016_pr, 0.46154189_pr]
       W_g = [-6.44209382_pr, -105.47508640_pr, 340.70046970_pr, -1421.01017122_pr, 0._pr]*alfa_dme
       B_g = zero
       H_g = [-5.01983486_pr,  -43.51414535_pr, 151.96320622_pr,  -652.79409940_pr, 0._pr]*alfa_dme
       M_g = zero
       W_g(n_g) = -sum(W_g(1:n_g-1))
       H_g(n_g) = -sum(H_g(1:n_g-1))
    !---------------------------------------------------------------------
    !
    !---------------------------------------------------------------------
    Case ('REG_NLO')
       n_g = 5
       if(.not.allocated(mu_g)) &
            allocate(mu_g(1:n_g),W_g(1:n_g),B_g(1:n_g),H_g(1:n_g),M_g(1:n_g))
       ! Regulator: (1-exp(-(r/2.0)**2))**6 !arXiv:nucl-th/0703087 couplings
       mu_g= [ 2.37939092_pr, 1.65287253_pr, 1.27520880_pr, 1.05227325_pr,0.94589848_pr]
       W_g = [ 0.03922136_pr, 0.31333037_pr,-1.10781372_pr, 1.33022277_pr,0._pr]*alfa_dme
       B_g = zero
       H_g = [-0.07844272_pr,-0.62666073_pr, 2.21562744_pr,-2.66044555_pr,0._pr]*alfa_dme
       M_g = zero
       W_g(n_g) = -sum(W_g(1:n_g-1))
       H_g(n_g) = -sum(H_g(1:n_g-1))
    !---------------------------------------------------------------------
    !
    !---------------------------------------------------------------------
    Case ('REG_N2LO')
       n_g = 5
       if(.not.allocated(mu_g)) &
            allocate(mu_g(1:n_g),W_g(1:n_g),B_g(1:n_g),H_g(1:n_g),M_g(1:n_g))
       ! Regulator: (1-exp(-(r/2.0)**2))**6 !arXiv:nucl-th/0703087 couplings
       mu_g= [ 2.07116798_pr, 1.44082837_pr, 1.17971841_pr, 1.07792473_pr,0.80689959_pr]
       W_g = [-0.32944606_pr,-2.45084067_pr,10.53421272_pr,-8.28227070_pr,0._pr]*alfa_dme
       B_g = zero
       H_g = [-0.26899125_pr,-0.09173571_pr, 2.51928834_pr,-2.39784450_pr,0._pr]*alfa_dme
       M_g = zero
       W_g(n_g) = -sum(W_g(1:n_g-1))
       H_g(n_g) = -sum(H_g(1:n_g-1))
    !---------------------------------------------------------------------
    !
    !---------------------------------------------------------------------
    Case ('REG_NLOD')
       n_g = 5
       if(.not.allocated(mu_g)) &
            allocate(mu_g(1:n_g),W_g(1:n_g),B_g(1:n_g),H_g(1:n_g),M_g(1:n_g))
       ! Regulator: (1-exp(-(r/2.0)**2))**6 !arXiv:nucl-th/0703087 couplings
       mu_g= [ 2.04448163_pr, 1.42618471_pr, 1.16558729_pr, 1.08201918_pr,0.75977094_pr]
       W_g = [-0.22903475_pr,-2.08771132_pr, 9.57867662_pr,-7.52300104_pr,0._pr]*alfa_dme
       B_g = zero
       H_g = [-0.34062315_pr,-0.40443314_pr, 4.40965332_pr,-3.85817607_pr,0._pr]*alfa_dme
       M_g = zero
       W_g(n_g) = -sum(W_g(1:n_g-1))
       H_g(n_g) = -sum(H_g(1:n_g-1))
    !---------------------------------------------------------------------
    !
    !---------------------------------------------------------------------
    Case ('REG_N2LOD')
       n_g = 5
       if(.not.allocated(mu_g)) &
            allocate(mu_g(1:n_g),W_g(1:n_g),B_g(1:n_g),H_g(1:n_g),M_g(1:n_g))
       ! Regulator: (1-exp(-(r/2.0)**2))**6 !arXiv:nucl-th/0703087 couplings
       mu_g= [ 1.99519747_pr, 1.46863759_pr, 1.11530962_pr, 1.03566073_pr,0.72565489_pr]
       W_g = [-0.17566300_pr,-1.74758603_pr, 9.72532149_pr,-8.12099199_pr,0._pr]*alfa_dme
       B_g = zero
       H_g = [-0.30265081_pr,-0.24487933_pr, 3.68915445_pr,-3.30302305_pr,0._pr]*alfa_dme
       M_g = zero
       W_g(n_g) = -sum(W_g(1:n_g-1))
       H_g(n_g) = -sum(H_g(1:n_g-1))
    !---------------------------------------------------------------------
    ! Default
    !---------------------------------------------------------------------
    Case default
        !Write(6,'("No Gogny interaction defined in routine skforce()")')
    End Select INTERACTION
    !
    Return
  End Subroutine gogny_force

  !=====================================================================
  !> Set up parameters of 2N density dependent couplings
  !=====================================================================
  Subroutine set_dme_couplings(fname)
    Implicit None
    Character (30), Intent(inout) :: fname
    !
    INTERACTION: Select Case (Trim(fname))
    Case ('DME_LO')
       if(.not.g_couplings_allocated) then
          ndme = 10
          allocate(grhorho(0:1,1:ndme))
          allocate(grhotau(0:1,1:ndme))
          allocate(grhodelrho(0:1,1:ndme))
          allocate(gJabJab(0:1,1:ndme))
          allocate(gJabJba(0:1,1:ndme))
          g_couplings_allocated = .true.
       endif
       ! Regulator: (1-exp(-(r/1.0)**2))**6
       grhorho(0,1:ndme) = [-100.155039_pr, -28.438527_pr,&
            1.138265_pr, 0.768756_pr, -15.821948_pr, 11.854066_pr,&
            0.517911_pr, 52.338074_pr, 7.632898_pr, 0.393321_pr]
       grhorho(1,1:ndme) = [33.385013_pr, 9.479528_pr, 1.138263_pr,&
            0.768755_pr, 5.273979_pr, 11.854066_pr, 0.517911_pr,&
            -17.446026_pr, 7.632896_pr, 0.393321_pr]
       grhodelrho(0,1:ndme) = [-139.617047_pr, 144.233763_pr,&
            0.031183_pr, 0.341619_pr, 148.719020_pr, 7.594671_pr,&
            0.341634_pr, -62.300609_pr, 6.033415_pr, 0.338381_pr]
       grhodelrho(1,1:ndme) = [46.539016_pr, -48.197890_pr,&
            0.031106_pr, 0.341618_pr, -49.573796_pr, 7.594647_pr,&
            0.341635_pr, 20.767403_pr, 6.033455_pr, 0.338381_pr]
       grhotau(0,1:ndme) = [558.468188_pr, -579.814641_pr,&
            0.031031_pr, 0.341617_pr, -594.905278_pr, 7.594597_pr,&
            0.341636_pr, 249.222134_pr, 6.033540_pr, 0.338381_pr]
       grhotau(1,1:ndme) = [-186.156062_pr, 192.311590_pr,&
            0.031183_pr, 0.341619_pr, 198.291975_pr, 7.594671_pr,&
            0.341634_pr, -83.067442_pr, 6.033414_pr, 0.338381_pr]
       gJabJab(0,1:ndme) = [527.802147_pr, 37.346104_pr, 1.750514_pr,&
            0.577058_pr, -442.802549_pr, 7.402071_pr, 0.335753_pr,&
            129.216649_pr, 5.501146_pr, 0.344248_pr]
       gJabJab(1,1:ndme) = [-175.934049_pr, -12.451375_pr,&
            1.750427_pr, 0.577043_pr, 147.593359_pr, 7.402158_pr,&
            0.335752_pr, -43.065677_pr, 5.501176_pr, 0.344257_pr]
       gJabJba(0,1:ndme) = [-512.489732_pr, 3.076003_pr, 1.989115_pr,&
            1.895073_pr, 437.871655_pr, 7.580670_pr, 0.339210_pr,&
            -175.167608_pr, 3.863425_pr, 0.264389_pr]
       gJabJba(1,1:ndme) = [170.829911_pr, -1.025331_pr, 1.989111_pr,&
            1.895080_pr, -145.957209_pr, 7.580670_pr, 0.339210_pr,&
            58.389178_pr, 3.863427_pr, 0.264389_pr]
    Case ('DME_NLO')
       if(.not.g_couplings_allocated) then
          ndme = 10
          allocate(grhorho(0:1,1:ndme))
          allocate(grhotau(0:1,1:ndme))
          allocate(grhodelrho(0:1,1:ndme))
          allocate(gJabJab(0:1,1:ndme))
          allocate(gJabJba(0:1,1:ndme))
          g_couplings_allocated = .true.
       endif
       ! Regulator: (1-exp(-(r/1.0)**2))**6
       grhorho(0,1:ndme) = [-95.034273_pr, -35.620677_pr, 1.132264_pr,&
            0.686283_pr, 3.821740_pr, 31.842150_pr, 0.587428_pr,&
            40.304134_pr, 6.324965_pr, 0.384860_pr]
       grhorho(1,1:ndme) = [18.169396_pr, 16.055814_pr, 1.662214_pr,&
            0.717249_pr, -9.906771_pr, 6.663177_pr, 0.483036_pr,&
            -4.938342_pr, 11.793629_pr, 0.428712_pr]
       grhodelrho(0,1:ndme) = [-137.887947_pr, 23.350183_pr,&
            0.224458_pr, 0.357310_pr, 144.454268_pr, 7.633521_pr,&
            0.341008_pr, -60.360303_pr, 5.948286_pr, 0.338928_pr]
       grhodelrho(1,1:ndme) = [42.212513_pr, -10.806162_pr,&
            0.793588_pr, 0.432688_pr, -50.517888_pr, 8.084105_pr,&
            0.347838_pr, 26.490516_pr, 5.414234_pr, 0.312814_pr]
       grhotau(0,1:ndme) = [551.551688_pr, -148.298527_pr,&
            0.135207_pr, 0.355601_pr, -575.556468_pr, 7.634596_pr,&
            0.340849_pr, 239.556876_pr, 5.944447_pr, 0.339278_pr]
       grhotau(1,1:ndme) = [-168.850052_pr, 43.224935_pr, 0.793593_pr,&
            0.432689_pr, 202.072185_pr, 8.084118_pr, 0.347838_pr,&
            -105.962726_pr, 5.414234_pr, 0.312813_pr]
       gJabJab(0,1:ndme) = [547.082063_pr, -51.726688_pr, 0.057647_pr,&
            1.734899_pr, -538.288457_pr, 8.198387_pr, 0.345969_pr,&
            223.493357_pr, 5.191101_pr, 0.278005_pr]
       gJabJab(1,1:ndme) = [-180.870520_pr, 34.958255_pr, 0.040984_pr,&
            1.150272_pr, 179.710564_pr, 8.197049_pr, 0.346024_pr,&
            -75.586920_pr, 5.172580_pr, 0.277883_pr]
       gJabJba(0,1:ndme) = [-508.541109_pr, 4.591227_pr, 1.473392_pr,&
            1.645378_pr, 437.641722_pr, 7.576360_pr, 0.339277_pr,&
            -176.982238_pr, 3.863401_pr, 0.265233_pr]
       gJabJba(1,1:ndme) = [174.778902_pr, -0.330045_pr, 6.258210_pr,&
            2.943817_pr, -146.294655_pr, 7.598849_pr, 0.338717_pr,&
            57.348953_pr, 3.802699_pr, 0.258690_pr]
    Case ('DME_N2LO')
       if(.not.g_couplings_allocated) then
          ndme = 10
          allocate(grhorho(0:1,1:ndme))
          allocate(grhotau(0:1,1:ndme))
          allocate(grhodelrho(0:1,1:ndme))
          allocate(gJabJab(0:1,1:ndme))
          allocate(gJabJba(0:1,1:ndme))
          g_couplings_allocated = .true.
       endif
       ! Regulator: (1-exp(-(r/1.0)**2))**6
       grhorho(0,1:ndme) = [-81.307748_pr, -33.753980_pr, 1.306145_pr,&
            0.709767_pr, 16.186426_pr, 6.122874_pr, 0.512274_pr,&
            25.821490_pr, 9.137382_pr, 0.405681_pr]
       grhorho(1,1:ndme) = [72.019958_pr, -178.783260_pr, 1.151572_pr,&
            0.453792_pr, 175.092909_pr, 1.125708_pr, 0.229926_pr,&
            -11.188956_pr, 7.439195_pr, 0.385899_pr]
       grhodelrho(0,1:ndme) = [-134.316488_pr, 18.853075_pr,&
            0.502924_pr, 0.387684_pr, 144.393799_pr, 7.724463_pr,&
            0.342261_pr, -63.715222_pr, 5.787681_pr, 0.333528_pr]
       grhodelrho(1,1:ndme) = [54.725625_pr, -15.227928_pr,&
            2.300536_pr, 0.465151_pr, -140.409598_pr, 5.419387_pr,&
            0.360903_pr, 82.280273_pr, 6.443836_pr, 0.361026_pr]
       grhotau(0,1:ndme) = [537.265949_pr, -75.418371_pr, 0.502605_pr,&
            0.387674_pr, -577.498849_pr, 7.724523_pr, 0.342256_pr,&
            254.800947_pr, 5.787467_pr, 0.333537_pr]
       grhotau(1,1:ndme) = [-218.902498_pr, 60.911809_pr, 2.300531_pr,&
            0.465151_pr, 561.639989_pr, 5.419381_pr, 0.360903_pr,&
            -329.122175_pr, 6.443831_pr, 0.361026_pr]
       gJabJab(0,1:ndme) = [591.174485_pr, -6.828089_pr, 7.113982_pr,&
            0.843972_pr, -515.116816_pr, 8.093103_pr, 0.343036_pr,&
            196.468150_pr, 4.841304_pr, 0.270649_pr]
       gJabJab(1,1:ndme) = [-140.759938_pr, 39.171924_pr, 0.721912_pr,&
            0.651303_pr, 163.161638_pr, 7.874211_pr, 0.341125_pr,&
            -99.737083_pr, 3.717561_pr, 0.263900_pr]
       gJabJba(0,1:ndme) = [-499.601652_pr, 9.168367_pr, 1.041510_pr,&
            1.327319_pr, 435.546444_pr, 7.567443_pr, 0.339130_pr,&
            -182.058328_pr, 3.782662_pr, 0.264702_pr]
       gJabJba(1,1:ndme) = [171.798686_pr, -0.772782_pr, 2.736644_pr,&
            2.154478_pr, -145.626933_pr, 7.582073_pr, 0.338881_pr,&
            58.224499_pr, 3.791282_pr, 0.261566_pr]
    Case ('DME_NLOD')
       if(.not.g_couplings_allocated) then
          ndme = 10
          allocate(grhorho(0:1,1:ndme))
          allocate(grhotau(0:1,1:ndme))
          allocate(grhodelrho(0:1,1:ndme))
          allocate(gJabJab(0:1,1:ndme))
          allocate(gJabJba(0:1,1:ndme))
          g_couplings_allocated = .true.
       endif
       ! Regulator: (1-exp(-(r/1.0)**2))**6
       grhorho(0,1:ndme) = [-60.729225_pr, -47.026444_pr, 1.462277_pr,&
            0.733509_pr, 22.435449_pr, -6.327016_pr, 0.496394_pr,&
            21.160364_pr, 9.994769_pr, 0.413067_pr]
       grhorho(1,1:ndme) = [54.639601_pr, -76.745147_pr, 1.250050_pr,&
            0.458945_pr, 71.882433_pr, 1.258592_pr, 0.236504_pr,&
            -12.797102_pr, 6.969074_pr, 0.378133_pr]
       grhodelrho(0,1:ndme) = [-130.131592_pr, 31.390112_pr,&
            0.698767_pr, 0.411737_pr, 162.671192_pr, 8.250544_pr,&
            0.351491_pr, -84.611660_pr, 5.711832_pr, 0.310235_pr]
       grhodelrho(1,1:ndme) = [50.594162_pr, -5.762380_pr,&
            2.341080_pr, 0.473718_pr, -47.101672_pr, 7.776034_pr,&
            0.342374_pr, 19.775230_pr, 5.861644_pr, 0.335990_pr]
       grhotau(0,1:ndme) = [520.526369_pr, -125.560923_pr,&
            0.698765_pr, 0.411737_pr, -650.685773_pr, 8.250552_pr,&
            0.351491_pr, 338.447473_pr, 5.711837_pr, 0.310234_pr]
       grhotau(1,1:ndme) = [-202.376635_pr, 22.980123_pr, 2.349403_pr,&
            0.474296_pr, 187.577963_pr, 7.780659_pr, 0.342245_pr,&
            -78.518832_pr, 5.851443_pr, 0.336190_pr]
       gJabJab(0,1:ndme) = [595.017460_pr, -9.689236_pr, 5.673100_pr,&
            0.780652_pr, -517.813409_pr, 8.120002_pr, 0.343439_pr,&
            197.623810_pr, 4.921126_pr, 0.271641_pr]
       gJabJab(1,1:ndme) = [-146.711250_pr, 30.310738_pr, 0.716817_pr,&
            0.698707_pr, 164.756392_pr, 7.876200_pr, 0.341192_pr,&
            -94.620827_pr, 3.821893_pr, 0.264514_pr]
       gJabJba(0,1:ndme) = [-501.602430_pr, 7.980880_pr, 1.102319_pr,&
            1.384138_pr, 436.087270_pr, 7.569374_pr, 0.339174_pr,&
            -180.839087_pr, 3.804112_pr, 0.264915_pr]
       gJabJba(1,1:ndme) = [174.274117_pr, -0.379184_pr, 5.143002_pr,&
            2.756397_pr, -146.273289_pr, 7.596776_pr, 0.338880_pr,&
            57.325528_pr, 3.830719_pr, 0.260263_pr]
    Case ('DME_N2LOD')
       if(.not.g_couplings_allocated) then
          ndme = 10
          allocate(grhorho(0:1,1:ndme))
          allocate(grhotau(0:1,1:ndme))
          allocate(grhodelrho(0:1,1:ndme))
          allocate(gJabJab(0:1,1:ndme))
          allocate(gJabJba(0:1,1:ndme))
          g_couplings_allocated = .true.
       endif
       ! Regulator: (1-exp(-(r/1.0)**2))**6
       grhorho(0,1:ndme) = [-69.362363_pr, -3.716361_pr, 42.489137_pr,&
            0.815098_pr, -31.065394_pr, 1.821175_pr, 0.600350_pr,&
            38.257488_pr, 7.895239_pr, 0.395152_pr]
       grhorho(1,1:ndme) = [64.436286_pr, -157.267798_pr, 1.062161_pr,&
            0.459791_pr, 164.571284_pr, 1.054032_pr, 0.233574_pr,&
            -11.891157_pr, 7.275178_pr, 0.382997_pr]
       grhodelrho(0,1:ndme) = [-133.326346_pr, 56.534352_pr,&
            0.183237_pr, 0.377193_pr, 151.580052_pr, 7.841617_pr,&
            0.345227_pr, -70.503644_pr, 5.738554_pr, 0.324708_pr]
       grhodelrho(1,1:ndme) = [52.273597_pr, -14.020060_pr,&
            2.154094_pr, 0.458527_pr, -134.117215_pr, 5.481912_pr,&
            0.359922_pr, 78.282624_pr, 6.503408_pr, 0.361173_pr]
       grhotau(0,1:ndme) = [533.304686_pr, -122.473434_pr,&
            0.644681_pr, 0.400138_pr, -677.584364_pr, 8.483326_pr,&
            0.354348_pr, 349.940698_pr, 5.925302_pr, 0.307174_pr]
       grhotau(1,1:ndme) = [-209.094388_pr, 56.080420_pr, 2.154086_pr,&
            0.458527_pr, 536.472037_pr, 5.481899_pr, 0.359922_pr,&
            -313.132633_pr, 6.503397_pr, 0.361173_pr]
       gJabJab(0,1:ndme) = [600.702478_pr, 63.514580_pr, 1.005265_pr,&
            0.594125_pr, -438.289059_pr, 7.367382_pr, 0.334442_pr,&
            90.449147_pr, 6.434219_pr, 0.368371_pr]
       gJabJab(1,1:ndme) = [-142.697983_pr, 25.854127_pr, 0.811625_pr,&
            0.809808_pr, 163.126515_pr, 7.857003_pr, 0.340547_pr,&
            -96.364085_pr, 3.672051_pr, 0.261887_pr]
       gJabJba(0,1:ndme) = [-496.115206_pr, 11.117231_pr, 0.974858_pr,&
            1.255730_pr, 434.667538_pr, 7.564703_pr, 0.339059_pr,&
            -184.189064_pr, 3.747382_pr, 0.264327_pr]
       gJabJba(1,1:ndme) = [173.446826_pr, -0.487609_pr, 3.938371_pr,&
            2.501895_pr, -146.115172_pr, 7.591693_pr, 0.338957_pr,&
            57.522085_pr, 3.835298_pr, 0.261364_pr]
    Case ('REG_LO')
       if(.not.g_couplings_allocated) then
          ndme = 10
          allocate(grhorho(0:1,1:ndme))
          allocate(grhotau(0:1,1:ndme))
          allocate(grhodelrho(0:1,1:ndme))
          allocate(gJabJab(0:1,1:ndme))
          allocate(gJabJba(0:1,1:ndme))
          g_couplings_allocated = .true.
       endif
       ! Regulator: (1-exp(-(r/2.0)**2))**6
       grhorho(0,1:ndme) = [-54.763722_pr, -32.456823_pr, 4.257666_pr,&
            0.741074_pr, 105.635227_pr, 5.740131_pr, 0.383692_pr,&
            -44.012364_pr, 4.506705_pr, 0.261245_pr]
       grhorho(1,1:ndme) = [18.254574_pr, 10.818930_pr, 4.257668_pr,&
            0.741074_pr, -35.211737_pr, 5.740132_pr, 0.383692_pr,&
            14.670790_pr, 4.506705_pr, 0.261245_pr]
       grhodelrho(0,1:ndme) = [-118.301284_pr, 43.833034_pr,&
            1.512348_pr, 0.504812_pr, 115.686395_pr, 8.203794_pr,&
            0.340797_pr, -65.421946_pr, 4.911474_pr, 0.337662_pr]
       grhodelrho(1,1:ndme) = [39.433761_pr, -14.610981_pr,&
            1.512349_pr, 0.504812_pr, -38.562090_pr, 8.203796_pr,&
            0.340797_pr, 21.807265_pr, 4.911477_pr, 0.337662_pr]
       grhotau(0,1:ndme) = [473.205136_pr, -175.331613_pr,&
            1.512348_pr, 0.504812_pr, -462.745444_pr, 8.203793_pr,&
            0.340797_pr, 261.687351_pr, 4.911479_pr, 0.337662_pr]
       grhotau(1,1:ndme) = [-157.735045_pr, 58.444258_pr, 1.512347_pr,&
            0.504812_pr, 154.248894_pr, 8.203790_pr, 0.340797_pr,&
            -87.229658_pr, 4.911469_pr, 0.337661_pr]
       gJabJab(0,1:ndme) = [412.501669_pr, 38.995081_pr, 6.863103_pr,&
            0.730785_pr, -393.891190_pr, 7.555354_pr, 0.334689_pr,&
            124.356461_pr, 6.117904_pr, 0.377254_pr]
       gJabJab(1,1:ndme) = [-137.500556_pr, -12.999342_pr,&
            6.862950_pr, 0.730774_pr, 131.298832_pr, 7.555315_pr,&
            0.334689_pr, -41.452923_pr, 6.117889_pr, 0.377251_pr]
       gJabJba(0,1:ndme) = [-382.156208_pr, -41.565137_pr,&
            6.349674_pr, 0.721412_pr, 351.439623_pr, 7.420112_pr,&
            0.334476_pr, -104.333993_pr, 5.701018_pr, 0.377973_pr]
       gJabJba(1,1:ndme) = [127.385403_pr, 13.855224_pr, 6.349649_pr,&
            0.721410_pr, -117.146736_pr, 7.420108_pr, 0.334476_pr,&
            34.778056_pr, 5.701013_pr, 0.377973_pr]
    Case ('REG_NLO')
       if(.not.g_couplings_allocated) then
          ndme = 10
          allocate(grhorho(0:1,1:ndme))
          allocate(grhotau(0:1,1:ndme))
          allocate(grhodelrho(0:1,1:ndme))
          allocate(gJabJab(0:1,1:ndme))
          allocate(gJabJba(0:1,1:ndme))
          g_couplings_allocated = .true.
       endif
       ! Regulator: (1-exp(-(r/2.0)**2))**6
       grhorho(0,1:ndme) = [-54.131350_pr, -32.397872_pr, 4.296646_pr,&
            0.741389_pr, 104.141639_pr, 5.782615_pr, 0.384517_pr,&
            -43.112995_pr, 4.533136_pr, 0.261736_pr]
       grhorho(1,1:ndme) = [16.905426_pr, 11.048674_pr, 4.530510_pr,&
            0.747943_pr, -32.454503_pr, 5.987334_pr, 0.388793_pr,&
            12.873152_pr, 4.671826_pr, 0.264169_pr]
       grhodelrho(0,1:ndme) = [-117.753612_pr, 43.414626_pr,&
            1.525689_pr, 0.505053_pr, 115.109221_pr, 8.206171_pr,&
            0.340696_pr, -64.873880_pr, 4.947158_pr, 0.339169_pr]
       grhodelrho(1,1:ndme) = [38.422549_pr, -14.500775_pr,&
            1.586918_pr, 0.507663_pr, -37.425341_pr, 8.238574_pr,&
            0.340367_pr, 21.120368_pr, 5.052262_pr, 0.345004_pr]
       grhotau(0,1:ndme) = [471.014449_pr, -173.657969_pr,&
            1.525689_pr, 0.505053_pr, -460.436675_pr, 8.206170_pr,&
            0.340696_pr, 259.495034_pr, 4.947163_pr, 0.339169_pr]
       grhotau(1,1:ndme) = [-153.690196_pr, 58.003239_pr, 1.586917_pr,&
            0.507663_pr, 149.701708_pr, 8.238568_pr, 0.340367_pr,&
            -84.481811_pr, 5.052258_pr, 0.345004_pr]
       gJabJab(0,1:ndme) = [417.490394_pr, 41.164740_pr, 6.572937_pr,&
            0.720794_pr, -397.414189_pr, 7.543554_pr, 0.334732_pr,&
            124.580282_pr, 6.048107_pr, 0.375051_pr]
       gJabJab(1,1:ndme) = [-138.747244_pr, -13.530687_pr,&
            6.647836_pr, 0.723334_pr, 132.183489_pr, 7.546356_pr,&
            0.334722_pr, -41.516436_pr, 6.065462_pr, 0.375604_pr]
       gJabJba(0,1:ndme) = [-381.381484_pr, -41.165695_pr,&
            6.397553_pr, 0.723372_pr, 350.876000_pr, 7.422613_pr,&
            0.334471_pr, -104.307941_pr, 5.714211_pr, 0.378441_pr]
       gJabJba(1,1:ndme) = [128.160127_pr, 14.251434_pr, 6.211398_pr,&
            0.715758_pr, -117.699860_pr, 7.412904_pr, 0.334491_pr,&
            34.798144_pr, 5.662085_pr, 0.376601_pr]
    Case ('REG_N2LO')
       if(.not.g_couplings_allocated) then
          ndme = 10
          allocate(grhorho(0:1,1:ndme))
          allocate(grhotau(0:1,1:ndme))
          allocate(grhodelrho(0:1,1:ndme))
          allocate(gJabJab(0:1,1:ndme))
          allocate(gJabJba(0:1,1:ndme))
          g_couplings_allocated = .true.
       endif
       ! Regulator: (1-exp(-(r/2.0)**2))**6
       grhorho(0,1:ndme) = [-53.089688_pr, -32.641240_pr, 4.358885_pr,&
            0.742628_pr, 102.078739_pr, 5.842177_pr, 0.385747_pr,&
            -41.736661_pr, 4.574103_pr, 0.262458_pr]
       grhorho(1,1:ndme) = [20.073555_pr, 9.614783_pr, 4.020561_pr,&
            0.736588_pr, -38.074133_pr, 5.528728_pr, 0.378882_pr,&
            16.984576_pr, 4.329241_pr, 0.258309_pr]
       grhodelrho(0,1:ndme) = [-116.982517_pr, 43.299818_pr,&
            1.544140_pr, 0.505757_pr, 114.249352_pr, 8.214121_pr,&
            0.340582_pr, -64.337892_pr, 4.983001_pr, 0.341033_pr]
       grhodelrho(1,1:ndme) = [40.536967_pr, -13.975573_pr,&
            1.443725_pr, 0.499906_pr, -39.831537_pr, 8.160440_pr,&
            0.341188_pr, 22.146581_pr, 4.849738_pr, 0.331535_pr]
       grhotau(0,1:ndme) = [467.930070_pr, -173.198740_pr,&
            1.544140_pr, 0.505757_pr, -456.997263_pr, 8.214121_pr,&
            0.340581_pr, 257.351129_pr, 4.983005_pr, 0.341033_pr]
       grhotau(1,1:ndme) = [-162.147867_pr, 55.902582_pr, 1.443724_pr,&
            0.499907_pr, 159.326636_pr, 8.160436_pr, 0.341188_pr,&
            -88.586857_pr, 4.849732_pr, 0.331534_pr]
       gJabJab(0,1:ndme) = [425.134265_pr, 43.759921_pr, 6.156364_pr,&
            0.705696_pr, -402.659635_pr, 7.525135_pr, 0.334805_pr,&
            125.199240_pr, 5.945098_pr, 0.371458_pr]
       gJabJab(1,1:ndme) = [-131.812469_pr, 59.589223_pr, 1.491836_pr,&
            0.520415_pr, 119.224429_pr, 8.086509_pr, 0.339777_pr,&
            -72.455927_pr, 4.376595_pr, 0.341709_pr]
       gJabJba(0,1:ndme) = [-380.161048_pr, -40.617695_pr,&
            6.474960_pr, 0.726762_pr, 349.930262_pr, 7.427341_pr,&
            0.334460_pr, -104.180687_pr, 5.735868_pr, 0.379291_pr]
       gJabJba(1,1:ndme) = [127.753316_pr, 14.072570_pr, 6.283631_pr,&
            0.718923_pr, -117.396331_pr, 7.417283_pr, 0.334481_pr,&
            34.762245_pr, 5.682984_pr, 0.377415_pr]
    Case ('REG_NLOD')
       if(.not.g_couplings_allocated) then
          ndme = 10
          allocate(grhorho(0:1,1:ndme))
          allocate(grhotau(0:1,1:ndme))
          allocate(grhodelrho(0:1,1:ndme))
          allocate(gJabJab(0:1,1:ndme))
          allocate(gJabJba(0:1,1:ndme))
          g_couplings_allocated = .true.
       endif
       ! Regulator: (1-exp(-(r/2.0)**2))**6
       grhorho(0,1:ndme) = [-52.224958_pr, -33.346077_pr, 4.387677_pr,&
            0.743522_pr, 100.830014_pr, 5.874594_pr, 0.386511_pr,&
            -40.652826_pr, 4.606491_pr, 0.262927_pr]
       grhorho(1,1:ndme) = [18.993949_pr, 10.041829_pr, 4.206461_pr,&
            0.739530_pr, -36.110048_pr, 5.680770_pr, 0.382158_pr,&
            15.576752_pr, 4.435400_pr, 0.260262_pr]
       grhodelrho(0,1:ndme) = [-116.496527_pr, 43.850990_pr,&
            1.551795_pr, 0.506521_pr, 113.782635_pr, 8.222762_pr,&
            0.340570_pr, -64.408435_pr, 4.983880_pr, 0.341586_pr]
       grhodelrho(1,1:ndme) = [39.800019_pr, -13.973493_pr,&
            1.498647_pr, 0.502560_pr, -38.901221_pr, 8.181468_pr,&
            0.340775_pr, 21.615156_pr, 4.935325_pr, 0.337086_pr]
       grhotau(0,1:ndme) = [465.986107_pr, -175.403431_pr,&
            1.551795_pr, 0.506521_pr, -455.130393_pr, 8.222762_pr,&
            0.340570_pr, 257.633304_pr, 4.983884_pr, 0.341586_pr]
       grhotau(1,1:ndme) = [-159.200076_pr, 55.894172_pr, 1.498646_pr,&
            0.502560_pr, 155.605229_pr, 8.181464_pr, 0.340776_pr,&
            -86.461000_pr, 4.935320_pr, 0.337086_pr]
       gJabJab(0,1:ndme) = [424.650825_pr, 43.273594_pr, 6.182681_pr,&
            0.706271_pr, -402.413363_pr, 7.525045_pr, 0.334806_pr,&
            125.381234_pr, 5.951354_pr, 0.371550_pr]
       gJabJab(1,1:ndme) = [-133.238227_pr, 60.202612_pr, 1.468844_pr,&
            0.520204_pr, 120.189772_pr, 8.076596_pr, 0.339832_pr,&
            -73.280613_pr, 4.311912_pr, 0.339165_pr]
       gJabJba(0,1:ndme) = [-380.416965_pr, -40.730963_pr,&
            6.458032_pr, 0.726001_pr, 350.135623_pr, 7.426270_pr,&
            0.334462_pr, -104.213844_pr, 5.731220_pr, 0.379102_pr]
       gJabJba(1,1:ndme) = [128.024630_pr, 14.173716_pr, 6.236862_pr,&
            0.716778_pr, -117.604087_pr, 7.414094_pr, 0.334489_pr,&
            34.799360_pr, 5.668848_pr, 0.376834_pr]
    Case ('REG_N2LOD')
       if(.not.g_couplings_allocated) then
          ndme = 10
          allocate(grhorho(0:1,1:ndme))
          allocate(grhotau(0:1,1:ndme))
          allocate(grhodelrho(0:1,1:ndme))
          allocate(gJabJab(0:1,1:ndme))
          allocate(gJabJba(0:1,1:ndme))
          g_couplings_allocated = .true.
       endif
       ! Regulator: (1-exp(-(r/2.0)**2))**6
       grhorho(0,1:ndme) = [-53.452298_pr, -33.453735_pr, 4.291127_pr,&
            0.741554_pr, 103.568127_pr, 5.789835_pr, 0.384815_pr,&
            -42.271469_pr, 4.556213_pr, 0.261962_pr]
       grhorho(1,1:ndme) = [19.199099_pr, 9.582627_pr, 4.222692_pr,&
            0.738698_pr, -36.145988_pr, 5.689039_pr, 0.382123_pr,&
            15.812670_pr, 4.419759_pr, 0.260201_pr]
       grhodelrho(0,1:ndme) = [-117.508059_pr, 44.413146_pr,&
            1.523658_pr, 0.505775_pr, 115.020487_pr, 8.213545_pr,&
            0.340770_pr, -65.431325_pr, 4.921112_pr, 0.338579_pr]
       grhodelrho(1,1:ndme) = [39.828933_pr, -13.439654_pr,&
            1.507258_pr, 0.501551_pr, -38.847471_pr, 8.171707_pr,&
            0.340612_pr, 21.264083_pr, 4.990061_pr, 0.338706_pr]
       grhotau(0,1:ndme) = [470.032235_pr, -177.652058_pr,&
            1.523658_pr, 0.505775_pr, -460.081814_pr, 8.213545_pr,&
            0.340770_pr, 261.724868_pr, 4.921116_pr, 0.338579_pr]
       grhotau(1,1:ndme) = [-159.315733_pr, 53.758755_pr, 1.507256_pr,&
            0.501552_pr, 155.390206_pr, 8.171703_pr, 0.340612_pr,&
            -85.056650_pr, 4.990058_pr, 0.338705_pr]
       gJabJab(0,1:ndme) = [423.255884_pr, 42.133965_pr, 6.244602_pr,&
            0.707424_pr, -401.713049_pr, 7.525065_pr, 0.334806_pr,&
            125.787744_pr, 5.969373_pr, 0.371845_pr]
       gJabJab(1,1:ndme) = [-133.598848_pr, 62.064676_pr, 1.458995_pr,&
            0.520958_pr, 120.413940_pr, 8.080629_pr, 0.339912_pr,&
            -74.452135_pr, 4.259669_pr, 0.337773_pr]
       gJabJba(0,1:ndme) = [-379.710883_pr, -40.424934_pr,&
            6.502859_pr, 0.727984_pr, 349.581468_pr, 7.429107_pr,&
            0.334456_pr, -104.129653_pr, 5.743890_pr, 0.379606_pr]
       gJabJba(1,1:ndme) = [127.895798_pr, 14.112635_pr, 6.260246_pr,&
            0.717779_pr, -117.509108_pr, 7.415419_pr, 0.334486_pr,&
            34.791080_pr, 5.675426_pr, 0.377083_pr]
    Case default
        !Write(6,'("No Gogny interaction defined in routine skforce()")')
    End Select INTERACTION
    !
    Return
  End Subroutine set_dme_couplings

  !=====================================================================
  !> Set up parameters of 3N density dependent couplings
  !=====================================================================
  Subroutine set_3N_couplings(fname)
    Implicit None
    Character (30), Intent(inout) :: fname
    !
    INTERACTION: Select Case (Trim(fname))
    Case ('DME_N2LO')
       if(.not.h_couplings_allocated) then
          ndme = 10
          allocate(h_rho_rho_rho(0:3,1:ndme))
          allocate(h_rho_rho_tau(0:3,1:ndme))
          allocate(h_rho_rho_delrho(0:3,1:ndme))
          allocate(h_rho_nablarho(0:3,1:ndme))
          allocate(h_rho_J_nablarho(0:3,1:ndme))
          allocate(h_rho_Jab_Jab(0:3,1:ndme))
          allocate(h_rho_Jab_Jba(0:3,1:ndme))
          h_couplings_allocated = .true.
          h_rho_rho_rho = zero
          h_rho_rho_tau = zero
          h_rho_rho_delrho = zero
          h_rho_nablarho = zero
          h_rho_J_nablarho = zero
          h_rho_Jab_Jab = zero
          h_rho_Jab_Jba = zero
       endif
       ! Regulator: (1-exp(-(r/1.0)**2))**6
       h_rho_rho_rho(0,1:ndme) = [-191.466291_pr, -3.947512_pr,&
            1104.152403_pr, 1.016715_pr, -207.609004_pr, 2.779447_pr,&
            0.632152_pr, 192.150193_pr, 7.873038_pr,&
            0.447590_pr] !rho0_rho0_rho0
       h_rho_rho_rho(1,1:ndme) = [191.466323_pr, 4.848435_pr,&
            291.781485_pr, 0.824588_pr, -137.787569_pr, 7.789136_pr,&
            0.493806_pr, 31.925473_pr, 3.742476_pr,&
            0.674936_pr] !rho0_rho1_rho1

       h_rho_rho_tau(0,1:ndme) = [1430.524224_pr, -504.226873_pr,&
            1.112788_pr, 0.612691_pr, 1708.935293_pr, 3.335513_pr,&
            0.357872_pr,-1165.393475_pr, 11.553869_pr,&
            0.399463_pr] !rho0_rho0_tau0
       h_rho_rho_tau(1,1:ndme) = [-677.093081_pr, 486.414982_pr,&
            1.053993_pr, 0.512865_pr, 1341.666344_pr, 9.641634_pr,&
            0.396463_pr, -981.108094_pr, 5.315148_pr,&
            0.282057_pr] !rho0_rho1_tau1
       h_rho_rho_tau(3,1:ndme) = [-476.833580_pr, 53.888013_pr,&
            1.056395_pr, 0.953378_pr, 723.637555_pr, 6.608471_pr,&
            0.355932_pr, -431.820474_pr, 4.176922_pr,&
            0.261339_pr] !rho1_rho1_tau0

       h_rho_rho_delrho(0,1:ndme) = [255.323143_pr, 2335.780279_pr,&
            1.062769_pr, 0.485762_pr, 3525.434810_pr, 9.076207_pr,&
            0.391973_pr,-3579.933780_pr, 5.280914_pr,&
            0.315419_pr] !rho0_rho0_delrho0
       h_rho_rho_delrho(1,1:ndme) = [276.441260_pr,-1042.600506_pr,&
            1.137639_pr, 0.465403_pr,-5108.905879_pr, 4.495298_pr,&
            0.377863_pr, 3686.491849_pr, 5.439784_pr,&
            0.370533_pr] !rho0_rho1_delrho1
       h_rho_rho_delrho(3,1:ndme) = [-85.142612_pr, -146.099199_pr,&
            1.968557_pr, 0.542293_pr, -64.450297_pr, 10.242367_pr,&
            0.347482_pr, 124.369477_pr, 8.882779_pr,&
            0.429860_pr] !rho1_rho1_delrho0

       h_rho_nablarho(0,1:ndme) = [2009.771482_pr, 561.291449_pr,&
            9.708717_pr, 0.663387_pr,-3982.884779_pr, 5.265705_pr,&
            0.417405_pr, 1905.968525_pr, 5.138535_pr,&
            0.317006_pr] !rho0_nablarho_0
       h_rho_nablarho(1,1:ndme) = [-669.923827_pr, -187.097140_pr,&
            9.708717_pr, 0.663387_pr, 1327.628207_pr, 5.265705_pr,&
            0.417405_pr, -635.322809_pr, 5.138535_pr,&
            0.317006_pr] !rho0_nablarho_1

       h_rho_J_nablarho(0,1:ndme) = [-622.293378_pr, -749.001174_pr,&
            1.857057_pr, 0.511437_pr, -323.426397_pr, 12.399832_pr,&
            0.352760_pr, 696.093796_pr, 8.424590_pr,&
            0.407726_pr] !rho0_J0_nablarho0
       h_rho_J_nablarho(1,1:ndme) = [207.431126_pr, 249.667061_pr,&
            1.857057_pr, 0.511437_pr, 107.808798_pr, 12.399832_pr,&
            0.352760_pr, -232.031265_pr, 8.424590_pr,&
            0.407726_pr] !rho0_J1_nablarho1
       h_rho_J_nablarho(2,1:ndme) = [-69.153447_pr, 18.447175_pr,&
            0.963997_pr, 0.408548_pr, 232.792946_pr, 4.811111_pr,&
            0.358285_pr, -138.461081_pr, 5.771041_pr,&
            0.359235_pr] !rho1_J0_nablarho1
       h_rho_J_nablarho(3,1:ndme) = [207.446176_pr, 241.319879_pr,&
            2.094487_pr, 0.467694_pr, -18.526094_pr, 5.301030_pr,&
            0.235338_pr, -145.681899_pr, 8.186773_pr,&
            0.382250_pr] !rho1_J1_nablarho0

       h_rho_Jab_Jab(0,1:ndme) = [-164.147136_pr, -389.181766_pr,&
            1.246970_pr, 0.673621_pr, -623.159226_pr, 15.137484_pr,&
            0.478502_pr, 753.788559_pr, 4.820061_pr,&
            0.242603_pr] !rho0_Jab0_Jab0
       h_rho_Jab_Jab(1,1:ndme) = [-221.833524_pr, -302.154245_pr,&
            2.907415_pr, 0.533511_pr, -85.138375_pr, 15.749395_pr,&
            0.353931_pr, 238.150934_pr, 12.684553_pr,&
            0.433143_pr] !rho0_Jab1_Jab1
       h_rho_Jab_Jab(2,1:ndme) = [109.492615_pr, 27.843741_pr,&
            1.525099_pr, 0.821725_pr, -118.559159_pr, 6.620652_pr,&
            0.428829_pr, 33.138910_pr, 9.020888_pr,&
            0.352926_pr] !rho1_Jab0_Jab1

       h_rho_Jab_Jba(0,1:ndme) = [754.160514_pr, 324.620705_pr,&
            1.332036_pr, 0.592914_pr,-1792.904724_pr, 5.720526_pr,&
            0.337042_pr, 762.958743_pr, 10.427306_pr,&
            0.392348_pr] !rho0_Jab0_Jba0
       h_rho_Jab_Jba(1,1:ndme) = [25.197525_pr, 206.262715_pr,&
            2.961757_pr, 0.523516_pr, 713.618970_pr, 6.235678_pr,&
            0.368977_pr, -550.174851_pr, 8.831554_pr,&
            0.378343_pr] !rho0_Jab1_Jab1
       h_rho_Jab_Jba(2,1:ndme) = [-502.800185_pr, 293.968840_pr,&
            2.696277_pr, 0.524940_pr, -411.145254_pr, 2.632525_pr,&
            0.410116_pr, 303.806843_pr, 4.395629_pr,&
            0.248732_pr] !rho1_Jab0_Jab1
    Case ('DME_NLOD')
       if(.not.h_couplings_allocated) then
          ndme = 10
          allocate(h_rho_rho_rho(0:3,1:ndme))
          allocate(h_rho_rho_tau(0:3,1:ndme))
          allocate(h_rho_rho_delrho(0:3,1:ndme))
          allocate(h_rho_nablarho(0:3,1:ndme))
          allocate(h_rho_J_nablarho(0:3,1:ndme))
          allocate(h_rho_Jab_Jab(0:3,1:ndme))
          allocate(h_rho_Jab_Jba(0:3,1:ndme))
          h_couplings_allocated = .true.
          h_rho_rho_rho = zero
          h_rho_rho_tau = zero
          h_rho_rho_delrho = zero
          h_rho_nablarho = zero
          h_rho_J_nablarho = zero
          h_rho_Jab_Jab = zero
          h_rho_Jab_Jba = zero
       endif
       ! Regulator: (1-exp(-(r/1.0)**2))**6
       h_rho_rho_rho(0,1:ndme) = [-194.604086_pr, -278.025048_pr,&
            16.162057_pr, 0.508052_pr, 824.854507_pr, 2.396056_pr,&
            0.255962_pr, -265.336723_pr, 2.724172_pr,&
            0.484238_pr] !rho0_rho0_rho0
       h_rho_rho_rho(1,1:ndme) = [194.608138_pr, 3.381308_pr,&
            862.008907_pr, 0.996858_pr, -541.735544_pr, 6.475590_pr,&
            0.620901_pr, 296.242114_pr, 7.499135_pr,&
            0.595847_pr] !rho0_rho1_rho1

       h_rho_rho_tau(0,1:ndme) = [727.382945_pr, 781.724205_pr,&
            3.261927_pr, 0.539204_pr, 2906.910365_pr, 5.937161_pr,&
            0.363585_pr,-2358.803330_pr, 9.663158_pr,&
            0.385059_pr] !rho0_rho0_tau0
       h_rho_rho_tau(1,1:ndme) = [-484.937933_pr, 296.087577_pr,&
            1.106016_pr, 0.510073_pr, -771.557200_pr, 3.259216_pr,&
            0.327849_pr, 455.089208_pr, 11.550233_pr,&
            0.401192_pr] !rho0_rho1_tau1
       h_rho_rho_tau(3,1:ndme) = [-242.458343_pr, -72.343068_pr,&
            42.043819_pr, 0.656860_pr, -67.353535_pr, 3.104972_pr,&
            0.509565_pr, 136.033968_pr, 14.373256_pr,&
            0.431670_pr] !rho1_rho1_tau0

       h_rho_rho_delrho(0,1:ndme) = [-24.500072_pr, 488.450027_pr,&
            32.534772_pr, 0.612765_pr,-2812.972057_pr, 4.009923_pr,&
            0.396005_pr, 1736.809657_pr, 4.660834_pr,&
            0.294118_pr] !rho0_rho0_delrho0
       h_rho_rho_delrho(1,1:ndme) = [242.395045_pr,-1613.786112_pr,&
            10.468058_pr, 0.449001_pr, 1535.231432_pr, 3.008052_pr,&
            0.219266_pr, -263.715175_pr, 1.743027_pr,&
            0.385994_pr] !rho0_rho1_delrho1
       h_rho_rho_delrho(3,1:ndme) = [8.160888_pr, -140.350181_pr,&
            1.704407_pr, 0.504991_pr, -379.514455_pr, 5.217533_pr,&
            0.413307_pr, 302.008716_pr, 6.525130_pr,&
            0.406468_pr] !rho1_rho1_delrho0

       h_rho_nablarho(0,1:ndme) = [1017.721445_pr, 509.507138_pr,&
            9.572072_pr, 0.664230_pr,-3195.334160_pr, 5.393699_pr,&
            0.422769_pr, 1633.181957_pr, 5.964534_pr,&
            0.335873_pr] !rho0_nablarho_0
       h_rho_nablarho(1,1:ndme) = [-339.240482_pr, -169.835738_pr,&
            9.572070_pr, 0.664230_pr, 1065.111531_pr, 5.393699_pr,&
            0.422769_pr, -544.394072_pr, 5.964534_pr,&
            0.335873_pr] !rho0_nablarho_1

       h_rho_J_nablarho(0,1:ndme) = [0.032009_pr, -540.285458_pr,&
            1.843541_pr, 0.509241_pr, -590.076762_pr, 10.158896_pr,&
            0.343519_pr, 611.327049_pr, 8.010173_pr,&
            0.401350_pr] !rho0_J0_nablarho0
       h_rho_J_nablarho(1,1:ndme) = [-0.010670_pr, 180.095170_pr,&
            1.843541_pr, 0.509241_pr, 196.692244_pr, 10.158897_pr,&
            0.343519_pr, -203.775684_pr, 8.010173_pr,&
            0.401350_pr] !rho0_J1_nablarho1
       h_rho_J_nablarho(2,1:ndme) = [0.000583_pr, -14.748693_pr,&
            2.101918_pr, 0.530066_pr, -13.074286_pr, 10.048895_pr,&
            0.343077_pr, 14.609144_pr, 8.542098_pr,&
            0.410781_pr] !rho1_J0_nablarho1
       h_rho_J_nablarho(3,1:ndme) = [-0.004082_pr, 103.240849_pr,&
            2.101918_pr, 0.530066_pr, 91.520005_pr, 10.048895_pr,&
            0.343077_pr, -102.264009_pr, 8.542098_pr,&
            0.410781_pr] !rho1_J1_nablarho0

       h_rho_Jab_Jab(0,1:ndme) = [-145.511580_pr, 527.447301_pr,&
            3.617475_pr, 0.510767_pr, 2459.461878_pr, 6.495625_pr,&
            0.371798_pr,-1756.482627_pr, 9.161623_pr,&
            0.378175_pr] !rho0_Jab0_Jab0
       h_rho_Jab_Jab(1,1:ndme) = [48.502687_pr, -230.040943_pr,&
            2.746221_pr, 0.528503_pr, -277.509382_pr, 11.269698_pr,&
            0.340991_pr, 258.871063_pr, 11.254291_pr,&
            0.420947_pr] !rho0_Jab1_Jab1
       h_rho_Jab_Jab(2,1:ndme) = [96.973246_pr, -430.867908_pr,&
            32.126610_pr, 0.694087_pr, 1309.518458_pr, 10.980102_pr,&
            0.555141_pr, -681.383292_pr, 14.179154_pr,&
            0.572591_pr] !rho1_Jab0_Jab1

       h_rho_Jab_Jba(0,1:ndme) = [581.852059_pr,-2155.783384_pr,&
            38.659968_pr, 0.690948_pr, 6057.720933_pr, 13.586243_pr,&
            0.561051_pr,-3121.831404_pr, 17.540421_pr,&
            0.578630_pr] !rho0_Jab0_Jba0
       h_rho_Jab_Jba(1,1:ndme) = [-193.976807_pr, 185.016602_pr,&
            2.667999_pr, 0.504703_pr, 734.585053_pr, 7.371384_pr,&
            0.349844_pr, -494.189807_pr, 9.182438_pr,&
            0.377984_pr] !rho0_Jab1_Jab1
       h_rho_Jab_Jba(2,1:ndme) = [-387.931096_pr, -65.663719_pr,&
            1.337612_pr, 0.533581_pr, 325.724412_pr, 7.870413_pr,&
            0.333829_pr, -74.721335_pr, 8.831724_pr,&
            0.385628_pr] !rho1_Jab0_Jab1
    Case ('DME_N2LOD')
       if(.not.h_couplings_allocated) then
          ndme = 10
          allocate(h_rho_rho_rho(0:3,1:ndme))
          allocate(h_rho_rho_tau(0:3,1:ndme))
          allocate(h_rho_rho_delrho(0:3,1:ndme))
          allocate(h_rho_nablarho(0:3,1:ndme))
          allocate(h_rho_J_nablarho(0:3,1:ndme))
          allocate(h_rho_Jab_Jab(0:3,1:ndme))
          allocate(h_rho_Jab_Jba(0:3,1:ndme))
          h_couplings_allocated = .true.
          h_rho_rho_rho = zero
          h_rho_rho_tau = zero
          h_rho_rho_delrho = zero
          h_rho_nablarho = zero
          h_rho_J_nablarho = zero
          h_rho_Jab_Jab = zero
          h_rho_Jab_Jba = zero
       endif
       ! Regulator: (1-exp(-(r/1.0)**2))**6
       h_rho_rho_rho(0,1:ndme) = [-176.477347_pr, -3.443058_pr,&
            1110.088151_pr, 1.015263_pr, -195.965861_pr, 2.783642_pr,&
            0.631110_pr, 180.161628_pr, 7.830985_pr,&
            0.445219_pr] !rho0_rho0_rho0
       h_rho_rho_rho(1,1:ndme) = [176.477099_pr, 4.372092_pr,&
            287.487907_pr, 0.818505_pr, -127.613630_pr, 7.839144_pr,&
            0.492529_pr, 29.806141_pr, 3.778158_pr,&
            0.675998_pr] !rho0_rho1_rho1

       h_rho_rho_tau(0,1:ndme) = [1372.560396_pr, -167.717044_pr,&
            1.306008_pr, 0.843996_pr, 1224.584683_pr, 3.751918_pr,&
            0.389866_pr,-1014.385207_pr, 11.722851_pr,&
            0.400092_pr] !rho0_rho0_tau0
       h_rho_rho_tau(1,1:ndme) = [-638.445069_pr, 457.174806_pr,&
            1.056120_pr, 0.509365_pr, 1262.555443_pr, 9.602309_pr,&
            0.394956_pr, -919.545593_pr, 5.353821_pr,&
            0.283449_pr] !rho0_rho1_tau1
       h_rho_rho_tau(3,1:ndme) = [-457.515489_pr, 54.928393_pr,&
            1.018488_pr, 0.920695_pr, 673.167961_pr, 6.571835_pr,&
            0.353287_pr, -401.280845_pr, 4.078226_pr,&
            0.259898_pr] !rho1_rho1_tau0

       h_rho_rho_delrho(0,1:ndme) = [257.875308_pr, 2190.910923_pr,&
            1.064936_pr, 0.485090_pr, 3326.140570_pr, 9.025694_pr,&
            0.391031_pr,-3370.812730_pr, 5.312109_pr,&
            0.316429_pr] !rho0_rho0_delrho0
       h_rho_rho_delrho(1,1:ndme) = [255.516378_pr, -943.238889_pr,&
            1.137930_pr, 0.464167_pr,-4681.481264_pr, 4.535639_pr,&
            0.376255_pr, 3368.930976_pr, 5.493675_pr,&
            0.370194_pr] !rho0_rho1_delrho1
       h_rho_rho_delrho(3,1:ndme) = [-85.991287_pr, -140.105444_pr,&
            1.973925_pr, 0.540903_pr, -58.178034_pr, 10.533944_pr,&
            0.348160_pr, 118.083103_pr, 8.962276_pr,&
            0.430090_pr] !rho1_rho1_delrho0

       h_rho_nablarho(0,1:ndme) = [1923.297085_pr, 515.940408_pr,&
            9.766996_pr, 0.663875_pr,-3703.390061_pr, 5.259877_pr,&
            0.416940_pr, 1763.117691_pr, 5.079985_pr,&
            0.315391_pr] !rho0_nablarho_0
       h_rho_nablarho(1,1:ndme) = [-641.099028_pr, -171.980136_pr,&
            9.766996_pr, 0.663875_pr, 1234.463354_pr, 5.259877_pr,&
            0.416940_pr, -587.705897_pr, 5.079985_pr,&
            0.315391_pr] !rho0_nablarho_1

       h_rho_J_nablarho(0,1:ndme) = [-622.296003_pr, -703.094529_pr,&
            1.857766_pr, 0.511570_pr, -275.737733_pr, 12.729047_pr,&
            0.354174_pr, 645.736226_pr, 8.453211_pr,&
            0.408005_pr] !rho0_J0_nablarho0
       h_rho_J_nablarho(1,1:ndme) = [207.432001_pr, 234.364843_pr,&
            1.857766_pr, 0.511570_pr, 91.912578_pr, 12.729047_pr,&
            0.354174_pr, -215.245409_pr, 8.453211_pr,&
            0.408005_pr] !rho0_J1_nablarho1
       h_rho_J_nablarho(2,1:ndme) = [-69.153457_pr, 16.355997_pr,&
            0.984893_pr, 0.403363_pr, 222.053513_pr, 4.933411_pr,&
            0.357082_pr, -130.520127_pr, 5.921573_pr,&
            0.360097_pr] !rho1_J0_nablarho1
       h_rho_J_nablarho(3,1:ndme) = [207.440487_pr, 276.340339_pr,&
            2.027786_pr, 0.446018_pr, -67.233336_pr, 3.290201_pr,&
            0.231004_pr, -133.339385_pr, 8.155965_pr,&
            0.380408_pr] !rho1_J1_nablarho0

       h_rho_Jab_Jab(0,1:ndme) = [-148.232069_pr, -233.244302_pr,&
            1.353780_pr, 0.712015_pr, 1621.741926_pr, 4.955686_pr,&
            0.342365_pr, -837.652008_pr, 10.392114_pr,&
            0.390496_pr] !rho0_Jab0_Jab0
       h_rho_Jab_Jab(1,1:ndme) = [-227.164733_pr, -283.783817_pr,&
            2.915483_pr, 0.533505_pr, -64.830056_pr, 16.850891_pr,&
            0.357301_pr, 219.006236_pr, 12.796697_pr,&
            0.433678_pr] !rho0_Jab1_Jab1
       h_rho_Jab_Jab(2,1:ndme) = [98.833338_pr, 5.660629_pr,&
            80.019361_pr, 0.662630_pr, 20.800791_pr, 1.919970_pr,&
            0.612323_pr, -41.761515_pr, 9.382391_pr,&
            0.461789_pr] !rho1_Jab0_Jab1

       h_rho_Jab_Jba(0,1:ndme) = [701.193159_pr, 309.839169_pr,&
            1.332891_pr, 0.588648_pr,-1630.404932_pr, 5.695844_pr,&
            0.336866_pr, 683.326582_pr, 10.461696_pr,&
            0.392584_pr] !rho0_Jab0_Jba0
       h_rho_Jab_Jba(1,1:ndme) = [42.855463_pr, 185.692310_pr,&
            2.994533_pr, 0.528115_pr, 615.861218_pr, 6.187809_pr,&
            0.370486_pr, -483.864806_pr, 8.862534_pr,&
            0.378740_pr] !rho0_Jab1_Jab1
       h_rho_Jab_Jba(2,1:ndme) = [-467.482142_pr, 301.199365_pr,&
            2.712915_pr, 0.539177_pr, -419.908951_pr, 2.672278_pr,&
            0.417835_pr, 302.750666_pr, 4.113724_pr,&
            0.246609_pr] !rho1_Jab0_Jab1
    Case ('REG_N2LO')
       if(.not.h_couplings_allocated) then
          ndme = 10
          allocate(h_rho_rho_rho(0:3,1:ndme))
          allocate(h_rho_rho_tau(0:3,1:ndme))
          allocate(h_rho_rho_delrho(0:3,1:ndme))
          allocate(h_rho_nablarho(0:3,1:ndme))
          allocate(h_rho_J_nablarho(0:3,1:ndme))
          allocate(h_rho_Jab_Jab(0:3,1:ndme))
          allocate(h_rho_Jab_Jba(0:3,1:ndme))
          h_couplings_allocated = .true.
          h_rho_rho_rho = zero
          h_rho_rho_tau = zero
          h_rho_rho_delrho = zero
          h_rho_nablarho = zero
          h_rho_J_nablarho = zero
          h_rho_Jab_Jab = zero
          h_rho_Jab_Jba = zero
       endif
       ! Regulator: (1-exp(-(r/2.0)**2))**6
       h_rho_rho_rho(0,1:ndme) = [-53.174211_pr, -29.118344_pr,&
            2.784291_pr, 0.754214_pr, -34.906549_pr, 15.081911_pr,&
            0.814075_pr, 49.413198_pr, 19.652126_pr,&
            0.537669_pr] !rho0_rho0_rho0
       h_rho_rho_rho(1,1:ndme) = [53.177322_pr, -14.610570_pr,&
            9.781653_pr, 1.031210_pr, 57.314238_pr, 6.867877_pr,&
            0.601993_pr, -45.310072_pr, 12.196177_pr,&
            0.468010_pr] !rho0_rho1_rho1

       h_rho_rho_tau(0,1:ndme) = [670.438267_pr, 256.974743_pr,&
            19.402823_pr, 0.820759_pr, 836.407794_pr, 6.552856_pr,&
            0.409727_pr, -810.713757_pr, 14.497902_pr,&
            0.412362_pr] !rho0_rho0_tau0
       h_rho_rho_tau(1,1:ndme) = [-333.977890_pr, 156.579524_pr,&
            70.318147_pr, 0.681091_pr, 965.959380_pr, 7.285391_pr,&
            0.518284_pr, -591.091334_pr, 10.415134_pr,&
            0.532411_pr] !rho0_rho1_tau1
       h_rho_rho_tau(3,1:ndme) = [-228.585004_pr, 115.928500_pr,&
            2.835112_pr, 0.488852_pr, 322.330410_pr, 6.246973_pr,&
            0.373409_pr, -190.323813_pr, 7.071707_pr,&
            0.435951_pr] !rho1_rho1_tau0

       h_rho_rho_delrho(0,1:ndme) = [84.198052_pr, 658.433803_pr,&
            34.523142_pr, 0.583642_pr, -684.115181_pr, 8.947343_pr,&
            0.539837_pr, 148.660942_pr, 5.251022_pr,&
            0.490609_pr] !rho0_rho0_delrho0
       h_rho_rho_delrho(1,1:ndme) = [140.687390_pr, -1924.733792_pr,&
            15.563650_pr, 0.674046_pr, 2549.961431_pr, 16.320613_pr,&
            0.629805_pr,-879.769126_pr, 18.135101_pr,&
            0.607554_pr] !rho0_rho1_delrho1
       h_rho_rho_delrho(3,1:ndme) = [-31.851630_pr, -81.887814_pr,&
            7.217597_pr, 0.726171_pr, -191.727195_pr, 4.276095_pr,&
            0.335013_pr, 159.075041_pr, 8.906028_pr,&
            0.409362_pr] !rho1_rho1_delrho0

       h_rho_nablarho(0,1:ndme) = [919.803583_pr, 189.327785_pr,&
            3.493615_pr, 0.628479_pr, -1270.096117_pr, 7.399942_pr,&
            0.445714_pr, 524.136377_pr, 6.459491_pr,&
            0.325399_pr] !rho0_nablarho_0
       h_rho_nablarho(1,1:ndme) = [-306.601194_pr, -63.109262_pr,&
            3.493615_pr, 0.628479_pr, 423.365372_pr, 7.399942_pr,&
            0.445714_pr, -174.712126_pr, 6.459491_pr,&
            0.325399_pr] !rho0_nablarho_1

       h_rho_J_nablarho(0,1:ndme) = [-276.920794_pr, -77.313721_pr,&
            143.703220_pr, 0.617370_pr, 808.938061_pr, 33.768895_pr,&
            0.763596_pr,-412.424694_pr, 40.844845_pr,&
            0.754131_pr] !rho0_J0_nablarho0
       h_rho_J_nablarho(1,1:ndme) = [92.306931_pr, 25.771240_pr,&
            143.703220_pr, 0.617370_pr, -269.646020_pr, 33.768895_pr,&
            0.763596_pr, 137.474898_pr, 40.844845_pr,&
            0.754131_pr] !rho0_J1_nablarho1
       h_rho_J_nablarho(2,1:ndme) = [-30.777602_pr, 24.100241_pr,&
            13.042504_pr, 0.589283_pr, 161.864462_pr, 4.233415_pr,&
            0.368074_pr, -104.670140_pr, 6.290028_pr,&
            0.378272_pr] !rho1_J0_nablarho1
       h_rho_J_nablarho(3,1:ndme) = [92.323815_pr, -49.185836_pr,&
            50.484276_pr, 0.802471_pr, -12.397222_pr, 67.322874_pr,&
            0.648095_pr, 3.903974_pr, 23.444232_pr,&
            0.916689_pr] !rho1_J1_nablarho0

       h_rho_Jab_Jab(0,1:ndme) = [-344.045436_pr, 133.614976_pr,&
            122.250197_pr, 0.609219_pr, 69.903775_pr, 42.310449_pr,&
            0.854580_pr, -9.453227_pr, 15.870551_pr,&
            0.975860_pr] !rho0_Jab0_Jab0
       h_rho_Jab_Jab(1,1:ndme) = [-8.460857_pr, -399.533141_pr,&
            8.447334_pr, 0.581015_pr, 542.325689_pr, 9.657477_pr,&
            0.418069_pr, -191.870243_pr, 7.927430_pr,&
            0.282199_pr] !rho0_Jab1_Jab1
       h_rho_Jab_Jab(2,1:ndme) = [229.240303_pr, 126.982692_pr,&
            0.891682_pr, 0.629554_pr, 116.771560_pr, 5.379794_pr,&
            0.579739_pr, -250.238302_pr, 3.569201_pr,&
            0.230681_pr] !rho1_Jab0_Jab1

       h_rho_Jab_Jba(0,1:ndme) = [766.413575_pr, -1371.090224_pr,&
            1.312862_pr, 0.484497_pr, -5365.450258_pr, 6.315596_pr,&
            0.400326_pr, 3949.866339_pr, 6.424856_pr,&
            0.351808_pr] !rho0_Jab0_Jba0
       h_rho_Jab_Jba(1,1:ndme) = [-132.430725_pr, 70.253022_pr,&
            98.487233_pr, 0.604421_pr, -3.654890_pr, 10.021272_pr,&
            1.016388_pr, 8.163778_pr, 58.627106_pr,&
            0.934008_pr] !rho0_Jab1_Jab1
       h_rho_Jab_Jba(2,1:ndme) = [-510.888861_pr, 45.305977_pr,&
            2.162530_pr, 0.743794_pr, 386.350207_pr, 7.294524_pr,&
            0.332712_pr, -126.453730_pr, 5.329835_pr,&
            0.472339_pr] !rho1_Jab0_Jab1
    Case ('REG_NLOD')
       if(.not.h_couplings_allocated) then
          ndme = 10
          allocate(h_rho_rho_rho(0:3,1:ndme))
          allocate(h_rho_rho_tau(0:3,1:ndme))
          allocate(h_rho_rho_delrho(0:3,1:ndme))
          allocate(h_rho_nablarho(0:3,1:ndme))
          allocate(h_rho_J_nablarho(0:3,1:ndme))
          allocate(h_rho_Jab_Jab(0:3,1:ndme))
          allocate(h_rho_Jab_Jba(0:3,1:ndme))
          h_couplings_allocated = .true.
          h_rho_rho_rho = zero
          h_rho_rho_tau = zero
          h_rho_rho_delrho = zero
          h_rho_nablarho = zero
          h_rho_J_nablarho = zero
          h_rho_Jab_Jab = zero
          h_rho_Jab_Jba = zero
       endif
       ! Regulator: (1-exp(-(r/2.0)**2))**6
       h_rho_rho_rho(0,1:ndme) = [-76.051241_pr, -51.477332_pr,&
            4.039349_pr, 0.738297_pr, -4.214395_pr, -69.329914_pr,&
            1.310747_pr, 44.055041_pr, 24.320593_pr,&
            0.582315_pr]!rho0_rho0_rho0
       h_rho_rho_rho(1,1:ndme) = [76.061427_pr, 20.021634_pr,&
            15.312310_pr, 0.949240_pr, 20.417583_pr, 4.165537_pr,&
            0.698388_pr, -42.366849_pr, 14.134631_pr,&
            0.501484_pr]!rho0_rho1_rho1

       h_rho_rho_tau(0,1:ndme) = [623.591420_pr, -450.255625_pr,&
            101.802547_pr, 0.862006_pr, 25.519051_pr, 52.481499_pr,&
            1.130609_pr, 5.148062_pr, 14.333303_pr,&
            0.979128_pr]!rho0_rho0_tau0
       h_rho_rho_tau(1,1:ndme) = [-361.540235_pr, 303.235520_pr,&
            55.955961_pr, 0.830184_pr, 334.076909_pr, 21.161269_pr,&
            0.574379_pr, -241.625632_pr, 33.942466_pr,&
            0.630450_pr]!rho0_rho1_tau1
       h_rho_rho_tau(3,1:ndme) = [-180.807293_pr, 48.703996_pr,&
            3.666395_pr, 0.518113_pr, 101.521635_pr, 18.871314_pr,&
            0.505725_pr, -36.802534_pr, 9.756513_pr,&
            0.589149_pr]!rho1_rho1_tau0

       h_rho_rho_delrho(0,1:ndme) = [-113.981748_pr, 1109.517919_pr,&
            25.348295_pr, 0.680399_pr, 114.379836_pr, 2.983957_pr,&
            0.505308_pr, -494.472936_pr, 20.144559_pr,&
            0.552155_pr]!rho0_rho0_delrho0
       h_rho_rho_delrho(1,1:ndme) = [180.819514_pr, -507.488039_pr,&
            13.713453_pr, 0.616593_pr,-1159.785535_pr, 4.776191_pr,&
            0.415566_pr, 892.106502_pr, 8.375177_pr,&
            0.440582_pr]!rho0_rho1_delrho1
       h_rho_rho_delrho(3,1:ndme) = [38.005057_pr, -30.265617_pr,&
            7.572168_pr, 0.769078_pr, -87.927441_pr, 6.236297_pr,&
            0.334683_pr, 56.228087_pr, 10.523051_pr,&
            0.457443_pr]!rho1_rho1_delrho0

       h_rho_nablarho(0,1:ndme) = [471.522189_pr, 131.235174_pr,&
            3.270266_pr, 0.610652_pr, -898.823330_pr, 7.453586_pr,&
            0.453551_pr, 411.491341_pr, 7.902356_pr,&
            0.351050_pr]!rho0_nablarho_0
       h_rho_nablarho(1,1:ndme) = [-157.174063_pr, -43.745061_pr,&
            3.270266_pr, 0.610652_pr, 299.607785_pr, 7.453586_pr,&
            0.453551_pr, -137.163784_pr, 7.902356_pr,&
            0.351050_pr]!rho0_nablarho_1

       h_rho_J_nablarho(0,1:ndme) = [-0.008332_pr, -301.427710_pr,&
            5.586949_pr, 0.571752_pr,-1586.921772_pr, 5.110533_pr,&
            0.353036_pr, 1119.909164_pr, 8.398249_pr,&
            0.385408_pr]!rho0_J0_nablarho0
       h_rho_J_nablarho(1,1:ndme) = [0.002777_pr, 100.475903_pr,&
            5.586949_pr, 0.571752_pr, 528.973923_pr, 5.110533_pr,&
            0.353036_pr, -373.303054_pr, 8.398249_pr,&
            0.385408_pr]!rho0_J1_nablarho1
       h_rho_J_nablarho(2,1:ndme) = [0.001428_pr, -20.061149_pr,&
            5.909583_pr, 0.561183_pr, 26.494498_pr, 8.134837_pr,&
            0.428964_pr, -9.295572_pr, 7.303396_pr,&
            0.281121_pr]!rho1_J0_nablarho1
       h_rho_J_nablarho(3,1:ndme) = [-0.009994_pr, 140.427051_pr,&
            5.909595_pr, 0.561183_pr, -185.457213_pr, 8.134935_pr,&
            0.428966_pr, 65.066491_pr, 7.303476_pr,&
            0.281122_pr]!rho1_J1_nablarho0

       h_rho_Jab_Jab(0,1:ndme) = [-108.495946_pr, 2158.139340_pr,&
            16.710475_pr, 0.674808_pr,-3642.354056_pr, 18.230941_pr,&
            0.594344_pr, 1473.160215_pr, 21.534934_pr,&
            0.573123_pr]!rho0_Jab0_Jab0
       h_rho_Jab_Jab(1,1:ndme) = [36.230657_pr, -75.696790_pr,&
            145.394878_pr, 0.618494_pr, -11.253850_pr, 23.750561_pr,&
            0.922365_pr, 28.585014_pr, 95.003343_pr,&
            0.758060_pr]!rho0_Jab1_Jab1
       h_rho_Jab_Jab(2,1:ndme) = [72.330081_pr, -201.266693_pr,&
            1.240318_pr, 0.488448_pr, -870.417433_pr, 5.683143_pr,&
            0.399533_pr, 644.700150_pr, 6.148983_pr,&
            0.353980_pr]!rho1_Jab0_Jab1

       h_rho_Jab_Jba(0,1:ndme) = [433.934674_pr,-2857.331880_pr,&
            16.568987_pr, 0.671883_pr, 4056.080112_pr, 18.045239_pr,&
            0.601839_pr,-1537.619323_pr, 20.945897_pr,&
            0.576644_pr]!rho0_Jab0_Jba0
       h_rho_Jab_Jba(1,1:ndme) = [-144.669537_pr, 298.347664_pr,&
            1.143682_pr, 0.450997_pr, 2069.016697_pr, 5.616825_pr,&
            0.386810_pr,-1434.125651_pr, 6.726089_pr,&
            0.357540_pr]!rho0_Jab1_Jab1
       h_rho_Jab_Jba(2,1:ndme) = [-289.267425_pr, -26.729480_pr,&
            6.913764_pr, 0.827537_pr, 254.178299_pr, 8.757916_pr,&
            0.338676_pr, -78.206575_pr, 5.641968_pr,&
            0.320652_pr]!rho1_Jab0_Jab1
    Case ('REG_N2LOD')
       if(.not.h_couplings_allocated) then
          ndme = 10
          allocate(h_rho_rho_rho(0:3,1:ndme))
          allocate(h_rho_rho_tau(0:3,1:ndme))
          allocate(h_rho_rho_delrho(0:3,1:ndme))
          allocate(h_rho_nablarho(0:3,1:ndme))
          allocate(h_rho_J_nablarho(0:3,1:ndme))
          allocate(h_rho_Jab_Jab(0:3,1:ndme))
          allocate(h_rho_Jab_Jba(0:3,1:ndme))
          h_couplings_allocated = .true.
          h_rho_rho_rho = zero
          h_rho_rho_tau = zero
          h_rho_rho_delrho = zero
          h_rho_nablarho = zero
          h_rho_J_nablarho = zero
          h_rho_Jab_Jab = zero
          h_rho_Jab_Jba = zero
       endif
       ! Regulator: (1-exp(-(r/2.0)**2))**6
       h_rho_rho_rho(0,1:ndme) = [-47.769852_pr, -25.448519_pr,&
            2.708041_pr, 0.760389_pr, -35.433746_pr, 14.384527_pr,&
            0.797819_pr, 46.854503_pr, 19.001265_pr,&
            0.532089_pr] !rho0_rho0_rho0
       h_rho_rho_rho(1,1:ndme) = [47.772324_pr, -14.530087_pr,&
            9.464791_pr, 1.018103_pr, 55.431043_pr, 6.880737_pr,&
            0.598810_pr, -42.711642_pr, 11.944693_pr,&
            0.464166_pr] !rho0_rho1_rho1

       h_rho_rho_tau(0,1:ndme) = [676.940564_pr, -540.956703_pr,&
            68.675017_pr, 0.755874_pr, 69.508101_pr, 27.838452_pr,&
            0.876080_pr, 0.116984_pr, 22.941453_pr,&
            1.662470_pr] !rho0_rho0_tau0
       h_rho_rho_tau(1,1:ndme) = [-307.884191_pr, 144.248726_pr,&
            66.158384_pr, 0.669495_pr, 938.006113_pr, 7.078193_pr,&
            0.511240_pr, -575.239082_pr, 10.119625_pr,&
            0.524465_pr] !rho0_rho1_tau1
       h_rho_rho_tau(3,1:ndme) = [-215.526513_pr, 108.537896_pr,&
            2.861427_pr, 0.493919_pr, 313.219957_pr, 6.025323_pr,&
            0.367174_pr, -184.742339_pr, 6.955290_pr,&
            0.431736_pr] !rho1_rho1_tau0

       h_rho_rho_delrho(0,1:ndme) = [92.538030_pr, 616.755699_pr,&
            34.665007_pr, 0.582633_pr, -661.256520_pr, 8.842686_pr,&
            0.537571_pr, 148.746449_pr, 5.334915_pr,&
            0.490511_pr] !rho0_rho0_delrho0
       h_rho_rho_delrho(1,1:ndme) = [126.897074_pr,-1795.132051_pr,&
            15.864080_pr, 0.674019_pr, 2377.752571_pr, 16.458504_pr,&
            0.628249_pr, -819.115385_pr, 18.199139_pr,&
            0.605380_pr] !rho0_rho1_delrho1
       h_rho_rho_delrho(3,1:ndme) = [-34.295640_pr, -78.383423_pr,&
            7.232659_pr, 0.725697_pr, -180.229940_pr, 4.240251_pr,&
            0.335211_pr, 151.268724_pr, 8.915010_pr,&
            0.408491_pr] !rho1_rho1_delrho0

       h_rho_nablarho(0,1:ndme) = [879.759213_pr, -132.044406_pr,&
            53.754550_pr, 0.768163_pr, 954.965911_pr, 4.297249_pr,&
            0.429899_pr, -775.084796_pr, 8.223972_pr,&
            0.461349_pr] !rho0_nablarho_0
       h_rho_nablarho(1,1:ndme) = [-293.253071_pr, 44.014802_pr,&
            53.754550_pr, 0.768163_pr, -318.321972_pr, 4.297249_pr,&
            0.429899_pr, 258.361600_pr, 8.223972_pr,&
            0.461349_pr] !rho0_nablarho_1

       h_rho_J_nablarho(0,1:ndme) = [-277.005913_pr, -301.718637_pr,&
            6.223003_pr, 0.602245_pr,-1195.896886_pr, 4.660756_pr,&
            0.357836_pr, 945.680337_pr, 8.578861_pr,&
            0.387590_pr] !rho0_J0_nablarho0
       h_rho_J_nablarho(1,1:ndme) = [92.335304_pr, 100.572880_pr,&
            6.223003_pr, 0.602245_pr, 398.632296_pr, 4.660756_pr,&
            0.357836_pr, -315.226780_pr, 8.578861_pr,&
            0.387590_pr] !rho0_J1_nablarho1
       h_rho_J_nablarho(2,1:ndme) = [-30.777649_pr, 23.669237_pr,&
            13.210381_pr, 0.587749_pr, 161.638384_pr, 4.301179_pr,&
            0.368066_pr, -104.335801_pr, 6.385055_pr,&
            0.378479_pr] !rho1_J0_nablarho1
       h_rho_J_nablarho(3,1:ndme) = [92.335600_pr, -58.807254_pr,&
            3.124296_pr, 0.551469_pr, -52.707671_pr, -17.505747_pr,&
            0.503602_pr, 33.542765_pr, 8.354767_pr,&
            0.560914_pr] !rho1_J1_nablarho0

       h_rho_Jab_Jab(0,1:ndme) = [-328.781010_pr, 123.188729_pr,&
            119.658444_pr, 0.608246_pr, 76.162169_pr, 31.277196_pr,&
            0.786223_pr, -13.101046_pr, 14.503848_pr,&
            0.899731_pr] !rho0_Jab0_Jab0
       h_rho_Jab_Jab(1,1:ndme) = [-13.569985_pr, 164.096356_pr,&
            -10.852037_pr, 0.676647_pr, -339.926493_pr, -5.879727_pr,&
            0.337758_pr, -281.082584_pr, -11.945846_pr,&
            0.404328_pr] !rho0_Jab1_Jab1
       h_rho_Jab_Jab(2,1:ndme) = [219.066350_pr, 175.193120_pr,&
            4.958994_pr, 0.638368_pr, -716.320065_pr, 4.510572_pr,&
            0.351510_pr, 339.011243_pr, 4.658399_pr,&
            0.317814_pr] !rho1_Jab0_Jab1

       h_rho_Jab_Jba(0,1:ndme) = [723.864139_pr,-3516.992037_pr,&
            14.950205_pr, 0.675326_pr, 4747.108948_pr, 15.232994_pr,&
            0.593125_pr,-1787.585780_pr, 17.255714_pr,&
            0.561010_pr] !rho0_Jab0_Jba0
       h_rho_Jab_Jba(1,1:ndme) = [-118.262777_pr, 60.735628_pr,&
            100.189141_pr, 0.604730_pr, -4.140489_pr, 9.702234_pr,&
            0.978339_pr, 8.663766_pr, 49.686878_pr,&
            0.888894_pr] !rho0_Jab1_Jab1
       h_rho_Jab_Jba(2,1:ndme) = [-482.575230_pr, -400.931456_pr,&
            4.093531_pr, 0.551541_pr, 930.591563_pr, 4.799363_pr,&
            0.293095_pr, -369.876340_pr, 2.990499_pr,&
            0.188919_pr] !rho1_Jab0_Jab1
    case default
    end Select INTERACTION
  end Subroutine set_3N_couplings

  !=======================================================================
  !> Set up set with all gaussian ranges (gogny and/or coulomb)
  !=======================================================================
  Subroutine set_all_gaussians(icoul)
    implicit none
    integer(ipr), intent(in) :: icoul
    if(allocated(mu_g_all)) deallocate(mu_g_all)
    if(icoul.lt.0) then
       coulomb_gaussian = .true.
       n_g_all = n_g + n_g_coul
       allocate(mu_g_all(1:n_g_all))
       if(finite_range) then
          mu_g_all(1:n_g) = mu_g(1:n_g)
       endif
       mu_g_all(n_g+1:n_g_all) = mu_g_coul(1:n_g_coul)
    elseif(finite_range) then
       n_g_all = n_g
       allocate(mu_g_all(1:n_g_all))
       mu_g_all = mu_g
    endif
  End Subroutine set_all_gaussians

  !=======================================================================
  !> Set up Pairing & Skyrme force parameters and their combinations
  !=======================================================================
  Subroutine skforce(fname,noForces)
    Implicit None
    Integer(ipr) :: noForces
    Real(pr) :: A,wls,TA7,TA8
    Real(pr) :: zero,one,two,three,four,five,six,seven,eight,nine
    Real(pr) :: half,pp16,pp24
    Character (30), Intent(inout) :: fname
    !
    zero = 0.0_pr; one = 1.0_pr; two = 2.0_pr; three = 3.0_pr; four = 4.0_pr
    five = 5.0_pr; six = 6.0_pr; seven = 7.0_pr; eight = 8.0_pr; nine = 9.0_pr
    half = 0.5_pr; pp16 = 16.0_pr; pp24 = 24.0_pr
    !
    ! Default for all forces if not modified
    hbzero = 1.0d0/0.04823_pr ! DMSHB0=1/hbzero
    sigma = one
    t0 = zero; x0 = zero
    t1 = zero; x1 = zero
    t2 = zero; x2 = one
    t3 = zero; x3 = one
    wls= zero; b4 = wls/two; b4p=wls/two
    te = zero; to = zero
    CExPar=1.0_pr
    !
    noForces=0 ! No forces at all
    !
    INTERACTION: Select Case (Trim(fname))
    !---------------------------------------------------------------------
    ! SIII, Beiner et al., Nucl. Phys. A 238, 29 (1975)
    !---------------------------------------------------------------------
    Case ('SIII')
        ! ph-Force
        noForces=1
        use_cm_cor = .True.
        hbzero = 20.73533_pr
        t0 = -.1128750d+04; x0 = +0.4500000_pr
        t1 = +.3950000d+03; x1 = +0.0000000_pr
        t2 = -.9500000d+02; x2 = +0.0000000_pr
        t3 = +.1400000d+05; x3 = +1.0000000_pr
        wls= +.1200000d+03; sigma = one
        b4=wls/two; b4p=wls/two
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -265.2500_pr, -340.0625_pr /)
    !---------------------------------------------------------------------
    ! SKM*, Bartel et al., Nucl. Phys. A 386, 79 (1982)
    !---------------------------------------------------------------------
    Case ('SKM*')
        ! ph-Force
        noForces=1
        use_cm_cor = .True.
        hbzero = 20.73_pr
        t0 = -.2645000d+04; x0 = +.0900000_pr
        t1 = +.4100000d+03; x1 = +.0000000_pr
        t2 = -.1350000d+03; x2 = +.0000000_pr
        t3 = +.1559500d+05; x3 = +.0000000_pr
        wls= +.1300000d+03; sigma = one/six
        b4=wls/two; b4p=wls/two
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -265.2500_pr, -340.0625_pr /)
    !---------------------------------------------------------------------
    ! SKM* without center of mass correction
    !---------------------------------------------------------------------
    Case ('SKM*mod')
        ! ph-Force
        noForces=1
        use_cm_cor = .False.
        hbzero = 20.73_pr
        t0 = -.2645000d+04; x0 = +.0900000_pr
        t1 = +.4100000d+03; x1 = +.0000000_pr
        t2 = -.1350000d+03; x2 = +.0000000_pr
        t3 = +.1559500d+05; x3 = +.0000000_pr
        wls= +.1300000d+03; sigma = one/six
        b4=wls/two; b4p=wls/two
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -265.2500_pr, -340.0625_pr /)
    !---------------------------------------------------------------------
    ! SKP, Dobaczewski et al., Nucl. Phys. A 422, 103 (1984)
    !---------------------------------------------------------------------
    Case ('SKP')
        ! ph-Force
        noForces=1
        use_cm_cor  = .True.
        use_j2terms = .True.
        hbzero = 20.730_pr
        t0 =-0.2931696d+04; x0 = 0.2921515_pr
        t1 = 0.3206182d+03; x1 = 0.6531765_pr
        t2 =-0.3374091d+03; x2 =-0.5373230_pr
        t3 = 0.1870896d+05; x3 = 0.1810269_pr
        wls= 0.1000000d+03; sigma=one/six
        b4=wls/two; b4p=wls/two
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -265.2500_pr, -340.0625_pr /)
    !---------------------------------------------------------------------
    ! SLY4, Chabanat et al. Nucl. Phys. A 635, 231 (1998) (unrounded)
    !---------------------------------------------------------------------
    Case ('SLY4')
        ! ph-Force
        noForces=1
        use_cm_cor = .True.
        hbzero = 20.735530_pr
        t0 =-0.2488913d+04; x0 = 0.8340000_pr
        t1 = 0.4868180d+03; x1 =-0.3440000_pr
        t2 =-0.5463950d+03; x2 =-1.0000000_pr
        t3 = 0.1377700d+05; x3 = 1.3540000_pr
        wls= 0.1230000d+03; sigma=one/six
        b4 = wls/two; b4p=wls/two
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -325.2500_pr, -340.0625_pr /) ! HFB
    !---------------------------------------------------------------------
    ! SLY4, Chabanat et al. Nucl. Phys. A 635, 231 (1998) (unrounded)
    ! Version withouth center of mass correction used in TDHFB calculations
    ! with regularized pairing; see PRL 116, 122504 (2016)
    !---------------------------------------------------------------------
    Case ('SLY4mod')
        ! ph-Force
        noForces=1
        use_cm_cor = .False.
        hbzero = 20.735530_pr
        t0 =-0.2488913d+04; x0 = 0.8340000_pr
        t1 = 0.4868180d+03; x1 =-0.3440000_pr
        t2 =-0.5463950d+03; x2 =-1.0000000_pr
        t3 = 0.1377700d+05; x3 = 1.3540000_pr
        wls= 0.1230000d+03; sigma=one/six
        b4 = wls/two; b4p=wls/two
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -370.000_pr, -370.000_pr /) ! HFB
    !---------------------------------------------------------------------
    ! SLY5, Chabanat et al. Nucl. Phys. A 635, 231 (1998) (unrounded)
    !---------------------------------------------------------------------
    Case ('SLY5')
        ! ph-Force
        noForces=1
        use_cm_cor  = .True.
        use_j2terms = .True.
        hbzero = 20.73553_pr
        t0 =-0.2483450d+04; x0 = 0.7760000_pr
        t1 = 0.4842300d+03; x1 =-0.3170000_pr
        t2 =-0.5566900d+03; x2 =-1.0000000_pr
        t3 = 0.1375700d+05; x3 = 1.2630000_pr
        wls= 0.1250000d+03; sigma=one/six
        b4 = wls/two; b4p=wls/two
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -291.5000_pr, -297.7402_pr /) ! HFB
    !---------------------------------------------------------------------
    ! SLY6, Chabanat et al. Nucl. Phys. A 635, 231 (1998) (unrounded)
    !---------------------------------------------------------------------
    Case ('SLY6')
        ! ph-Force
        noForces=1
        use_cm_cor      = .True.
        use_full_cm_cor = .True.
        hbzero = 20.73553_pr
        t0 =-0.2479500d+04; x0 = 0.8250000_pr
        t1 = 0.4621800d+03; x1 =-0.4650000_pr
        t2 =-0.4486100d+03; x2 =-1.0000000_pr
        t3 = 0.1367300d+05; x3 = 1.3550000_pr
        wls= 0.1220000d+03; sigma=one/six
        b4=wls/two; b4p=wls/two
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -291.5000_pr, -297.7402_pr /) ! HFB
    !---------------------------------------------------------------------
    ! SLY7, Chabanat et al. Nucl. Phys. A 635, 231 (1998) (unrounded)
    !---------------------------------------------------------------------
    Case ('SLY7')
        ! ph-Force
        noForces=1
        use_cm_cor      = .True.
        use_j2terms     = .True.
        use_full_cm_cor = .True.
        hbzero = 20.73553_pr
        t0 =-0.2480800d+04; x0 = 0.8480000_pr
        t1 = 0.4612900d+03; x1 =-0.4920000_pr
        t2 =-0.4339300d+03; x2 =-1.0000000_pr
        t3 = 0.1366900d+05; x3 = 1.3930000_pr
        wls= 0.1250000d+03; sigma=one/six
        b4=wls/two; b4p=wls/two
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -291.5000_pr, -297.7402_pr /) ! HFB
    !---------------------------------------------------------------------
    ! SKI3, P.G.-Reinhard et al. Nucl. Phys. A584, 467  (1995)
    !---------------------------------------------------------------------
    Case ('SKI3')
        ! ph-Force
        noForces=1
        use_cm_cor      = .True.
        use_full_cm_cor = .True.
        hbzero = 20.7525d0
        t0 =-0.176288d+04; x0 = 0.30830_pr
        t1 = 0.561608d+03; x1 =-1.17220_pr
        t2 =-0.227090d+03; x2 =-1.09070_pr
        t3 = 0.810620d+04; x3 = 1.29260_pr
        sigma=one/four
        b4 = 94.254_pr; b4p=zero
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -357.2324_pr, -388.5625_pr /)
    !---------------------------------------------------------------------
    ! SKO, P.-G. Reinhard et al. Phys. Rev. C 60, 014316 (1999)
    !---------------------------------------------------------------------
    Case ('SKO')
        ! ph-Force
        noForces=1
        use_cm_cor      = .True.
        use_full_cm_cor = .True.
        hbzero = 20.735530_pr
        t0 =-0.21036530d+04; x0 = -0.2107010_pr
        t1 = 0.30335200d+03; x1 = -2.8107520_pr
        t2 = 0.79167400d+03; x2 = -1.4615950_pr
        t3 = 0.13553252d+05; x3 = -0.4298810_pr
        wls= 0.12300000d+03; sigma=one/four
        b4 = 0.17657800d+03; b4p=-0.1987490d+03
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -259.0391_pr, -274.8433_pr /)
    !---------------------------------------------------------------------
    ! SKO', P.-G. Reinhard et al. Phys. Rev. C 60, 014316 (1999)
    !---------------------------------------------------------------------
    Case ('SKOP')
        ! ph-Force
        noForces=1
        use_cm_cor      = .True.
        use_full_cm_cor = .True.
        use_j2terms     = .True.
        hbzero = 20.735530_pr
        t0 =-0.20994190d+04; x0 = -0.029503_pr
        t1 = 0.30153100d+03; x1 = -1.325732_pr
        t2 = 0.15478100d+03; x2 = -2.323439_pr
        t3 = 0.13526464d+05; x3 = -0.147404_pr
        wls= 0.12300000d+03; sigma=one/four
        b4 = 0.14389500d+03; b4p=-0.828888d+02
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -253.753488522_pr, -274.679033802_pr /) ! See NB above
    !---------------------------------------------------------------------
    ! SKX, A.Brown Phys.Rev. C 58, 220 (1998)
    !---------------------------------------------------------------------
    Case ('SKX')
        ! ph-Force
        noForces=1
        use_cm_cor  = .True.
        use_j2terms = .True.
        hbzero = 20.73_pr
        t0 = -1445.300_pr; x0 = 0.340_pr
        t1 =   246.900_pr; x1 = 0.580_pr
        t2 =  -131.800_pr; x2 = 0.127_pr
        t3 = 12103.900_pr; x3 = 0.030_pr
        sigma=one/two
        b4 = 0.0743d+03; b4p=zero
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -259.0391_pr, -274.8433_pr /)
    !---------------------------------------------------------------------
    ! Gogny D1, Decharge et al. PRC 21, 1568 (1980)
    !---------------------------------------------------------------------
    Case ('D1')
        ! ph-Force
        noForces=1
        use_cm_cor   = .True.
        use_j2terms  = .False.
        finite_range = .True.
        hbzero = 20.73667552957479_pr
        t3 = 6.0_pr*1350.00_pr; x3 = one;
        wls= 115.000_pr;
        sigma=one/three
        b4=wls/two; b4p=wls/two
        ! No delta-pairing here. The pairing will be in the finite range part
        CpV1= zero
        CpV0= zero
        ! W, B, H, M parameters of the finite-range two-body
        Call gogny_force(fname)
    !---------------------------------------------------------------------
    ! Gogny D1S, Berger et al. CPC 63, 365 (1991)
    !---------------------------------------------------------------------
    Case ('D1S')
        ! ph-Force
        noForces=1
        use_cm_cor   = .True.
        use_j2terms  = .False.
        finite_range = .True.
        hbzero = 20.73667552957479_pr
        t3 = 6.0_pr*1390.600_pr; x3 = one;
        wls= 130.000_pr;
        sigma=one/three
        b4=wls/two; b4p=wls/two
        ! No delta-pairing here. The pairing will be in the finite range part
        CpV1= zero
        CpV0= zero
        ! W, B, H, M parameters of the finite-range two-body
        Call gogny_force(fname)
    !---------------------------------------------------------------------
    ! Gogny D1', Decharge et al. PRC 21, 1568 (1980)
    !---------------------------------------------------------------------
    Case ('D1p')
        ! ph-Force
        noForces=1
        use_cm_cor   = .True.
        use_j2terms  = .False.
        finite_range = .True.
        hbzero = 20.73667552957479_pr
        t3 = 6.0_pr*1350.00_pr; x3 = one;
        wls= 130.000_pr;
        sigma=one/three
        b4=wls/two; b4p=wls/two
        ! No delta-pairing here. The pairing will be in the finite range part
        CpV1= zero
        CpV0= zero
        ! W, B, H, M parameters of the finite-range two-body
        Call gogny_force(fname)
    !---------------------------------------------------------------------
    !
    !---------------------------------------------------------------------
    Case ('D1N')
        ! ph-Force
        noForces=1
        use_cm_cor   = .True.
        use_j2terms  = .False.
        finite_range = .True.
        hbzero = 20.73667552957479_pr
        t3 = 6.0_pr*1609.50_pr; x3 = one;
        wls= 115.000_pr;
        sigma=one/three
        b4=wls/two; b4p=wls/two
        ! No delta-pairing here. The pairing will be in the finite range part
        CpV1= zero
        CpV0= zero
        ! W, B, H, M parameters of the finite-range two-body
        Call gogny_force(fname)
    !---------------------------------------------------------------------
    !T0X0 A very simple skyrme functional
    !---------------------------------------------------------------------
    Case ('T0X0')
        noForces=1
        use_cm_cor  = .True.
        t0 =  0.0_pr; x0 = 0.0_pr;
        ! t0 =  -1128.75_pr; x0 = 0.45_pr;
        hbzero = 20.73667622931579_pr
        ! pp-Forces
        CpV1= zero
        CpV0= zero !No pairing here.
    !---------------------------------------------------------------------
    ! Default
    !---------------------------------------------------------------------
    Case default
        ! Write(6,'("No Skyrme interaction defined in routine skforce()")')
    End Select INTERACTION
    !
    If (noForces.Eq.1) Then
        ! obtain coupling constants
        Call C_from_t()
        ! Frequent combinations entering the energy
        tv1   =  t0*(one+half*x0)*half;    tv2 = t0*(x0+half)*half
        tv3   =  t3*(one+half*x3)/12.0_pr; tv4 = t3*(x3+half)/12.0_pr
        tv5   = (t1*(one+half*x1)+t2*(one+half*x2))/four
        tv6   = (t2*(half+x2)-t1*(half+x1))/four
        ts1   = (t2*(one+half*x2)-three*t1*(one+half*x1))/pp16
        ts2   = (t1*(half+x1)*three+t2*(half+x2))/pp16
        t4o3  =  four/three; t324 = t3/pp24
        ! Frequent combinations entering the potential
        t0s   =  t0*(one-x0)*half; t0a = t0*(one+x0*half)
        drs   = (t2*(one+x2)-t1*(one-x1))*three/pp16
        dra   = (t2*(one+half*x2)-three*t1*(one+half*x1))/eight
        ts    = (t1*(one-x1) + three*t2*(one+x2))/eight
        ta    = (t1*(one+half*x1) + t2*(one+half*x2))/four
        t3alp = t3*(two+sigma)*(two+x3)/pp24
        t3al0 = t3*(x3+half)/six; t3alm = t3*sigma*(one+two*x3)/pp24
        alp   = one + sigma; alm = sigma - one
        wla0  = CrdJ(0)+CrdJ(1); wla1  = CrdJ(0)-CrdJ(1);
        TA7   = zero; TA8 = zero
        If(use_j2terms) Then
           TA7=(T1*(ONE-X1)-T2*(ONE+X2))/eight + five*to/four
           TA8=-(T1*X1+T2*X2)/four             + five*(te+to)/four
        End If
        TB7 = TA7; TB8 = TA8*half
        if(use_3N_couplings) then
           use_3N_couplings = .False.
           override_3N_couplings = .True.
        endif
    End If
    !
    Return
  End Subroutine skforce
  !=======================================================================
  !> Define functional parameters
  !=======================================================================
  Subroutine set_functional_parameters(fname,lpr)
    Implicit None
    Logical, Intent(in) :: lpr
    Character (30), Intent(inout) :: fname
    Logical :: regularization
    Integer(ipr), Parameter :: lin=15
    !
    ! parameters
    FunctionalName=fname
    eps=Spacing(1.0_pr)
    Pi=4.0_pr*Atan(1.0_pr)
    kfconst=(1.50_pr*Pi**2)**(1.0_pr/3.0_pr)    ! (3Pi^2/2)^(1/3)
    CK=3.0_pr/5.0_pr*kfconst**2
    nuLambda=700.0_pr ; nufpi = 93.0_pr
    !
    Call Make_Parameter_Free_Useful_Combinations()
    !
    ! exact Hartree CHrho from INM
    CHrho=0.0_pr; !!!!If (dmeorder.eq.3) Call CHrho_from_NM()
    !
    If(use_INM) Then
       Call calculate_C_from_NM(E_NM,K_NM,SMASS_NM,RHO_NM,ASS_NM,LASS_NM,VMASS_NM)
    Else
       Crho(0)=Crho(0)+CHrho
    End If
    Call calculate_NM_properties()
    !
    Crho(0)=Crho(0)-CHrho
    !
    Call calculate_natural_units()
    !
    ! Print output
    !If(lpr) Then
    !   Call print_functional_parameters()
    !End If
    !
  End Subroutine set_functional_parameters
  !=======================================================================
  !
  !=======================================================================
  Subroutine print_functional_parameters(fname)
    Use HFBTHO_utilities, Only: lout,lfile
    Implicit None
    Character (30), Optional, Intent(in) :: fname
    Integer(ipr) :: iw
    !
    Do iw=lout,lfile
       Write(iw,'("  ---------------------------------------")')
       Write(iw,'("           UNEDF Module Version: ",a)') Trim(Version)
       Write(iw,'("         M.Kortelainen & M.Stoitsov ")')
       Write(iw,'("  ---------------------------------------")')
       !
       Write(iw,*)
       If(Present(fname)) Then
          Write(iw,'(2x,a," functional")') Trim(fname)
       Else
          Write(iw,'(2x,a," functional")') Trim(FunctionalName)
       End If
       Write(iw,'("  ----------------------------------------")')
       If(.not.is_nedf) Then
          Write(iw,'("  Crho(0)= ",g26.18,"; Crho(1)= ",g26.18)') Crho
          Write(iw,'("  CDrho(0)=",g26.18,"; CDrho(1)=",g26.18)') CDrho
          Write(iw,'("  Ctau(0)= ",g26.18,"; Ctau(1)= ",g26.18)') Ctau
          Write(iw,'("  CrDr(0)= ",g26.18,"; CrDr(1)= ",g26.18)') Crdr
          Write(iw,'("  CrdJ(0)= ",g26.18,"; CrdJ(1)= ",g26.18)') CrdJ
          Write(iw,'("  CJ(0)=   ",g26.18,"; CJ(1)=   ",g26.18)') CJ
          Write(iw,'("  CpV0(0)= ",g26.18,"; CpV0(1)= ",g26.18)') CpV0
          Write(iw,'("  CpV1(0)= ",g26.18,"; CpV1(1)= ",g26.18)') CpV1
          Write(iw,'("  sigma=   ",g26.18,"; hbzero=  ",g26.18)') sigma,hbzero
          Write(iw,'("  functional has DME couplings: ",L1)') force_is_dme
          if(override_3N_couplings) then
             Write(iw,'("  functional ",a," has no 3N force, ignoring namelist input")') trim(FunctionalName)
          else
             Write(iw,'("  use 3-Nucleon DME couplings: ",L1)') use_3n_couplings
          endif
       Else
          Write(iw,'("  a0=      ",g26.18,"; a1=      ",g26.18,"; a2=      ",g26.18)') a_NEDF
          Write(iw,'("  b0=      ",g26.18,"; b1=      ",g26.18,"; b2=      ",g26.18)') b_NEDF
          Write(iw,'("  c0=      ",g26.18,"; c1=      ",g26.18,"; c2=      ",g26.18)') c_NEDF
          Write(iw,'("  eta=     ",g26.18,"; W0=      ",g26.18)') eta_NEDF,W0_NEDF
          Write(iw,'("  CpV0(0)= ",g26.18,"; CpV0(1)= ",g26.18)') CpV0
          Write(iw,'("  CpV1(0)= ",g26.18,"; CpV1(1)= ",g26.18)') CpV1
          Write(iw,'("  hbzeron= ",g26.18,"; hbzerop= ",g26.18)') hbzeron,hbzerop
       End If
       ! Finite-range force (Gogny force)
       If(finite_range) Then
          Write(iw,*)
          Write(iw,'("  Finite-range potential")')
          Write(iw,'("  ----------------------------------------")')
          Write(iw,'("  mu=",10f26.18)') mu_g
          Write(iw,'("  W=",10f26.18)') W_g
          Write(iw,'("  B=",10f26.18)') B_g
          Write(iw,'("  H=",10f26.18)') H_g
          Write(iw,'("  M=",10f26.18)') M_g
       End If
       ! Common features
       Write(iw,'("  e^2 chrg=",g26.18,"; CExPar=  ",g26.18)') e2charg,CExPar
       Write(iw,'("  c.m. correction: ",L1,", chr. density in direct Coul: ",L1)') use_cm_cor,use_charge_density
       Write(iw,'("  use tensor terms: ",L1)') use_j2terms
       ! Natural units
       If(.not.is_nedf.And..Not.finite_range) Then
          Write(iw,*)
          Write(iw,'("  Coupling constants in natural units")')
          Write(iw,'("  ----------------------------------------")')
          Write(iw,'("  Crho_nu(0)= ",g26.18,"; Crho_nu(1)= ",g26.18)') nuCrho
          Write(iw,'("  CDrho_nu(0)=",g26.18,"; CDrho_nu(1)=",g26.18)') nuCDrho
          Write(iw,'("  Ctau_nu(0)= ",g26.18,"; Ctau_nu(1)= ",g26.18)') nuCtau
          Write(iw,'("  CrDr_nu(0)= ",g26.18,"; CrDr_nu(1)= ",g26.18)') nuCrdr
          Write(iw,'("  CrdJ_nu(0)= ",g26.18,"; CrdJ_nu(1)= ",g26.18)') nuCrdJ
          Write(iw,'("  CJ_nu(0)=   ",g26.18,"; CJ_nu(1)=   ",g26.18)') nuCJ
          Write(iw,'("  CpV0_nu(0)= ",g26.18,"; CpV0_nu(1)= ",g26.18)') nuCpV0
          Write(iw,'("  CpV1_nu(0)= ",g26.18,"; CpV1_nu(1)= ",g26.18)') nuCpV1
          Write(iw,'("  fpi_nu=     ",g26.18,"; Lambda_nu=  ",g26.18)') nufpi,nuLambda
       End If
       ! DME
       If(dmeorder.Ge.0) Then
          Write(iw,*)
          Write(iw,'("  DME parameters")')
          Write(iw,'("  ----------------------------------------")')
          Write(iw,'("       gA=",f12.6," mpi [1/fm]=",f12.6," fpi [1/fm]=",f12.6)') gA,mpi,fpi
          Write(iw,'("  c1 [fm]=",f12.6,"    c3 [fm]=",f12.6,"    c4 [fm]=",f12.6)') c1,c3,c4
          Write(iw,'("       cd=",f12.6,"         ce=",f12.6," LamX[1/fm]=",f12.6)') cd,ce,LambdaX
          Write(iw,'("  ->CHrho=",f12.6)') CHrho
          If(dmeorder.Ge.2) Write(iw,'("  use 3N terms: ",L1)') use_DME3N_terms
       End If
       Write(iw,*)
       Write(iw,'("  Nuclear matter properties")')
       Write(iw,'("  ----------------------------------------")')
       Write(iw,'("  E_NM=    ",g26.18,"; K_NM=     ",g26.18)') E_NM,K_NM
       Write(iw,'("  P_NM=    ",g26.18,"; RHO_NM=   ",g26.18)') P_NM,RHO_NM
       Write(iw,'("  ASS_NM=  ",g26.18,"; LASS_NM=  ",g26.18)') ASS_NM,LASS_NM
       Write(iw,'("  SMASS_NM=",g26.18,"; VMASS_NM= ",g26.18)') SMASS_NM,VMASS_NM
       !
       Call t_from_C()
       !
       ! (t,x) parametrization of the Skyrme functional
       If(.not.is_nedf) Then
          Write(iw,*)
          Write(iw,'("  Associated (t,x)-coupling constants")')
          Write(iw,'("  ----------------------------------------")')
          If(.Not.finite_range.or.force_is_dme) Then
             Write(iw,'("  t0=    ",g26.18,"; x0=     ",g26.18)') t0,x0
             Write(iw,'("  t1=    ",g26.18,"; x1=     ",g26.18)') t1,x1
             Write(iw,'("  t2=    ",g26.18,"; x2=     ",g26.18)') t2,x2
         End If
          Write(iw,'("  t3=    ",g26.18,"; x3=     ",g26.18)') t3,x3
          Write(iw,'("  b4=    ",g26.18,"; b4p=    ",g26.18)') b4,b4p
          Write(iw,'("  te=    ",g26.18,"; to=     ",g26.18)') te,to
          Write(iw,'("  sigma= ",g26.18,"; hbzero= ",g26.18)') sigma,hbzero
       End If
       !
       If(Print_Namelist) Then
          Write(iw,*)
          SELECTED_FUNCTIONAL: Select Case (Trim(FunctionalName))
          Case ("UNEDF","SKYRME")
                Write(iw,'("NAMELIST CONTENT (cannot be modified for functional names UNEDF,SKYRME)")')
                Write(iw,'("-----------------------------------------------------------------------")')
          Case ("FITS")
                Write(iw,'("NAMELIST CONTENT (advanced usage: modify all but not C-, NM-, and more...)")')
                Write(iw,'("--------------------------------------------------------------------------")')
          Case default
                Write(iw,'("NAMELIST CONTENT (copy/past to UNEDF_NAMELIST.DAT and modify)")')
                Write(iw,'("-------------------------------------------------------------")')
          End Select SELECTED_FUNCTIONAL
          Write(*,'(" !NB: FUNCTIONALNAME should be always in quotations")')
          Write(*,UNEDF_NAMELIST)
       End If
    End Do
  End Subroutine print_functional_parameters
  !=======================================================================
  !> Calculates coupling constants in natural units
  !=======================================================================
  Subroutine calculate_natural_units
    Implicit None
    nuCrho = Crho*(nufpi**2)/(mevfm**3)
    nuCdrho = Cdrho*(nufpi**2)*((nuLambda*nufpi*nufpi)**sigma)/(mevfm**(3.0_pr*(1.0_pr+sigma)))
    nuCtau = Ctau*((nufpi*nuLambda)**2)/(mevfm**5)
    nuCrDr = CrDr*((nufpi*nuLambda)**2)/(mevfm**5)
    nuCrdJ = CrdJ*((nufpi*nuLambda)**2)/(mevfm**5)
    nuCJ = CJ*((nufpi*nuLambda)**2)/(mevfm**5)
    nuCpV0 = CpV0*(nufpi**2)/(mevfm**3)
    nuCpV1 = CpV1*(nufpi**4)*nuLambda/(mevfm**6)
  End Subroutine calculate_natural_units
  !=======================================================================
  !> Calculates volume C-constants (and sigma) form NM properties
  !>
  !> Input: E,K,SMASS,RHO,ASS,LASS,VMASS,sigma_NM(optional)
  !>
  !> Output: Crho(0),Crho(1),Cdrho(0),Cdrho(1),Ctau(0),Ctau(0),sigma(optional)
  !>
  !> Options:
  !>  - When sigma_NM exists then 'sigma'=sigma_NM
  !>  - When sigma_NM does not exist then 'sigma' is defined from NM
  !=======================================================================
  Subroutine calculate_C_from_NM(E,K,SMASS,RHO,ASS,LASS,VMASS,sigma_NM)
    Implicit None
    Real(pr), Intent(in) :: E,K,SMASS,RHO,ASS,LASS,VMASS
    Real(pr), Intent(in), Optional :: sigma_NM
    Real(pr) :: aRho0Rho0,daRho0Rho0,ddaRho0Rho0,aRho1Rho1,daRho1Rho1,ddaRho1Rho1
    Real(pr) :: aRho0Tau0,daRho0Tau0,ddaRho0Tau0,aRho1Tau1,daRho1Tau1,ddaRho1Tau1
    Real(pr) :: u,tauc,rho2
    Real(pr) :: E_fr,P_fr,K_fr,SMASS_fr,ASS_fr,LASS_fr,KA_fr,VMASS_fr
    Real(pr),Parameter :: c13=1.0_pr/3.0_pr,c23=2.0_pr/3.0_pr
    !
    tauc=CK*RHO**c23; u=(kfconst/mpi)*RHO**c13; rho2=rho**2
    !
    Call calculate_U_parameters(RHO,RHO,tauc*RHO,tauc*RHO,0.0_pr,0.0_pr)
    if(finite_range) then
       call calculate_finite_range_NM(rho,E_fr,P_fr,K_fr,SMASS_fr,ASS_fr,LASS_fr,KA_fr,VMASS_fr)
    else
       E_fr=0_pr; P_fr=0_pr; K_fr=0_pr; SMASS_fr=0_pr
       ASS_fr=0_pr; LASS_fr=0_pr; KA_fr=0._pr; VMASS_fr=0_pr
    end if
    !
    aRho0Rho0=0.50_pr*(aRhoRho+bRhoRho)*mevfm
    aRho1Rho1=0.50_pr*(aRhoRho-bRhoRho)*mevfm
    aRho0Tau0=0.50_pr*(aRhoTau+bRhoTau)*mevfm
    aRho1Tau1=0.50_pr*(aRhoTau-bRhoTau)*mevfm
    daRho0Rho0=0.50_pr*(daRhoRho+dbRhoRho)*mevfm
    daRho1Rho1=0.50_pr*(daRhoRho-dbRhoRho)*mevfm
    daRho0Tau0=0.50_pr*(daRhoTau+dbRhoTau)*mevfm
    daRho1Tau1=0.50_pr*(daRhoTau-dbRhoTau)*mevfm
    ddaRho0Rho0=0.50_pr*(ddaRhoRho+ddbRhoRho)*mevfm
    ddaRho1Rho1=0.50_pr*(ddaRhoRho-ddbRhoRho)*mevfm
    ddaRho0Tau0=0.50_pr*(ddaRhoTau+ddbRhoTau)*mevfm
    ddaRho1Tau1=0.50_pr*(ddaRhoTau-ddbRhoTau)*mevfm
    !
    ! set/calculate sigma
    If (Present(sigma_NM)) Then
        sigma=sigma_NM
    Else
        sigma=((1.0_pr/3.0_pr)*(-(K-K_fr)+tauc*hbzero*(-3.0_pr+4.0_pr*(SMASS-SMASS_fr))-9.0_pr*(E-E_fr+P_fr/rho)+9.0_pr*RHO2*hRho0Rho0 &
             +21.0_pr*tauc*RHO2*hRho0Tau0+u*RHO*(daRho0Rho0+5.0_pr*tauc*daRho0Tau0 &
             +7.0_pr*RHO*dhRho0Rho0+11.0_pr*tauc*RHO*dhRho0Tau0+u*ddaRho0Rho0 &
             +u*tauc*ddaRho0Tau0+u*RHO*ddhRho0Rho0+u*tauc*RHO*ddhRho0Tau0))) &
             /(tauc*hbzero*(-3.0_pr+2.0_pr*(SMASS-SMASS_fr))+3.0_pr*(E-E_fr+P_fr/rho)+3.0_pr*RHO2*hRho0Rho0 &
             +3.0_pr*tauc*RHO2*hRho0Tau0+u*RHO*(daRho0Rho0+tauc*daRho0Tau0 &
             + RHO*dhRho0Rho0+tauc*RHO*dhRho0Tau0))
    End If
    !
    Crho(0)=(c13*(tauc*hbzero*(-3.0_pr+(2.0_pr-3.0_pr*sigma)*(SMASS-SMASS_fr)) &
        +3.0_pr*(1.0_pr+sigma)*(E-E_fr)+3.0_pr*P_fr/rho-3.0_pr*sigma*RHO*aRho0Rho0 &
        +3.0_pr*(1.0_pr-sigma)*RHO2*hRho0Rho0+3.0_pr*tauc*RHO2*hRho0Tau0 &
        +u*RHO*(daRho0Rho0+tauc*daRho0Tau0+RHO*dhRho0Rho0 &
        +tauc*RHO*dhRho0Tau0)))/(sigma*RHO)
    Cdrho(0)=(c13*RHO**(-1.0_pr-sigma)*(tauc*hbzero*(3.0_pr-2.0_pr*(SMASS-SMASS_fr))&
        -3.0_pr*(E-E_fr+P_fr/rho)-3.0_pr*RHO**2*hRho0Rho0-3.0_pr*tauc*RHO2*hRho0Tau0&
        -u*RHO*(daRho0Rho0+tauc*daRho0Tau0+RHO*dhRho0Rho0 &
        +tauc*RHO*dhRho0Tau0)))/sigma
    Ctau(0)=(hbzero*(SMASS-SMASS_fr-1.0_pr)-RHO*(aRho0Tau0+RHO*hRho0Tau0))/RHO
    !
    Crho(1)=(27.0_pr*(ASS-ASS_fr)*(1.0_pr+sigma)-9.0_pr*(LASS-LASS_fr) &
        +5.0_pr*tauc*hbzero*(5.0_pr-6.0_pr*(VMASS-VMASS_fr)+3.0_pr*sigma*(-4.0_pr+3.0_pr*(VMASS-VMASS_fr))) &
        +20.0_pr*tauc*(2.0_pr-3.0_pr*sigma)*RHO*aRho0Tau0 &
        +RHO*(-27.0_pr*sigma*aRho1Rho1+5.0_pr*tauc*(11.0_pr-12.0_pr*sigma)*RHO*hRho0Tau0 &
        -27.0_pr*(-1.0_pr+sigma)*RHO*hRho1Rho1+9.0_pr*tauc*(5.0_pr-3.0_pr*sigma)*RHO*hRho1Tau0 &
        +45.0_pr*tauc*RHO*hRho1Tau1+40.0_pr*tauc*Ctau(0)-60.0_pr*tauc*sigma*Ctau(0) &
        +5.0_pr*u*tauc*daRho0Tau0+9.0_pr*u*daRho1Rho1+15.0_pr*u*tauc*daRho1Tau1 &
        +5.0_pr*u*tauc*RHO*dhRho0Tau0+9.0_pr*u*RHO*dhRho1Rho1+9.0_pr*u*tauc*RHO*dhRho1Tau0 &
        +15.0_pr*u*tauc*RHO*dhRho1Tau1))/(27.0_pr*sigma*RHO)
    Cdrho(1)=-(RHO**(-1.0_pr-sigma)*(27.0_pr*(ASS-ASS_fr)-9.0_pr*(LASS-LASS_fr) &
        +5.0_pr*tauc*hbzero*(5.0_pr-6.0_pr*(VMASS-VMASS_fr))+40.0_pr*tauc*RHO*aRho0Tau0 &
        +55.0_pr*tauc*RHO2*hRho0Tau0+27.0_pr*RHO**2*hRho1Rho1+45.0_pr*tauc*RHO2*hRho1Tau0 &
        +45.0_pr*tauc*RHO2*hRho1Tau1+40.0_pr*tauc*RHO*Ctau(0) +5.0_pr*u*tauc*RHO*daRho0Tau0 &
        +9.0_pr*u*RHO*daRho1Rho1+15.0_pr*u*tauc*RHO*daRho1Tau1 &
        +5.0_pr*u*tauc*RHO2*dhRho0Tau0+9.0_pr*u*RHO2*dhRho1Rho1 &
        +9.0_pr*u*tauc*RHO2*dhRho1Tau0 +15.0_pr*u*tauc*RHO2*dhRho1Tau1))/(27.0_pr*sigma)
    Ctau(1)=(hbzero-hbzero*(VMASS-VMASS_fr)+RHO*(aRho0Tau0-aRho1Tau1+RHO*hRho0Tau0-RHO*hRho1Tau1+Ctau(0)))/RHO
    !
  End Subroutine calculate_C_from_NM

  !======================================================================
  !> Calculates the contribution to the Nuclear Matter properties comming
  !> from the finite range part of the functional
  !>
  !> The formulas for E,P,K,ASS,LASS and KASS are from
  !> \cite sellahewa2014 . The formulas for SMASS and VMASS where derived
  !> from the energy density as linear combinations of the derivatives
  !> with respect of \f$\tau_n\f$ and \f$\tau_p\f$.
  !>
  !> The formulas are as follows
  !> \f{align*}{
  !>   E^{\rm NM}(\rho_0) & =
  !>     \frac{1}{2}\sum_{i=1,n_g}\left\{A_0^i\rho_0+B_0^ig_0(\mu_i k_F)
  !>     \right\}, \\
  !>   P^{\rm NM}(\rho_0) &=
  !>     \sum_{i=1,n_g} \left\{\frac{1}{2} A_0^i\rho_0^2 +
  !>     B_0^ip_0(\mu_ik_F)\rho_0\right\}, \\
  !>   K^{\rm NM}(\rho_0) &= -3\sum_{i=1,n_g} B_0^ik_0(\mu_ik_F), \\
  !>   M_s(\rho_0) &= \frac{2m_N}{\hbar^2}\frac{1}{2}\sum_{i=1,n_g}
  !>     B_0^im_0(\mu_i k_F)\frac{3}{5k_F^2}, \\
  !>   a^{\rm NM}_{\rm sym}(\rho_0) &= \sum_{i=1,n_g}\left\{\frac{1}{2}
  !>     A_1^i\rho_0+\frac{1}{6}\left[B_{nn}^i s_1(\mu_i k_F) + B_{np}^i
  !>     s_2(\mu_i k_F) \right] \right\}, \\
  !>   L^{\rm NM}_{\rm sym}(\rho_0) &= \sum_{i=1,n_g}\left\{\frac{3}{2}
  !>     A_1^i\rho_0+\frac{1}{6}\left[B_{nn}^i l_1(\mu_i k_F) + B_{np}^i
  !>     l_2(\mu_i k_F) \right] \right\}, \\
  !>   \Delta K^{\rm NM}(\rho_0) &= -\frac{1}{6}\sum_{i=1,n_g}\left\{
  !>     B_{nn}^i k_1(\mu_i k_F) + B_{np}^ik_2(\mu_i k_F)\right\}, \\
  !>   M_v(\rho_0) &= \frac{2m_N}{\hbar^2}\frac{1}{2}\sum_{i=1,n_g}
  !>    \left\{ B_0^im_0(\mu_i k_F) + B_{nn}^im_1(\mu_i k_F) +
  !>    B_{np}^im_2(\mu_i k_F) \right\}\frac{3}{5k_F^2},
  !> \f}
  !> where
  !> \f{align*}{
  !>  k_F &= \left(\frac{3\pi^2}{2} \rho_0 \right)^{1/3}, \\
  !>  A_0^i &= \frac{\pi^{3/2}\mu_i^3}{4}(4W_i+2B_i-2H_i-M_i), \\
  !>  A_1^i &= \frac{\pi^{3/2}\mu_i^3}{4}(-2H_i-M_i), \\
  !>  B_{nn}^i &= -\frac{1}{\sqrt{\pi}}(W_i+2B_i-H_i-2M_i), \\
  !>  B_{np}^i &= -\frac{1}{\sqrt{\pi}}(-H_i-2M_i), \\
  !>  B_0^i &= B_{nn}^i + B_{np}^i, \\
  !>  g_0(q) &= \frac{1}{q^3}\left[2-3q^2-\left(2-q^2\right)e^{-q^2}
  !>     \right] + \sqrt{\pi} {\rm erf}(q), \\
  !>  p_0(q) &= \frac{1}{2q^3}\left[-2+ q^2+\left(2+q^2\right)e^{-q^2}
  !>     \right], \\
  !>  k_0(q) &= \frac{1}{q^3}\left[-6+ 2q^2+\left(6+4q^2+q^4\right)
  !>     e^{-q^2} \right], \\
  !>  m_0(q) &= \frac{1}{2q^3}\left[-6+ 3q^2+\left(6+3q^2\right)e^{-q^2}
  !>     \right], \\
  !>  s_1(q) &=\frac{1}{q}\left[1-\left(1+q^2\right)e^{-q^2}\right], \\
  !>  s_2(q) &=\frac{1}{q}\left[1-q^2-e^{-q^2}\right], \\
  !>  l_1(q) &=\frac{1}{q}\left[-1+\left(1+q^2+2q^4\right)e^{-q^2}
  !>     \right], \\
  !>  l_2(q) &=\frac{1}{q}\left[-1-q^2+\left(1+2q^2\right)e^{-q^2}
  !>     \right], \\
  !>  k_1(q) &=\frac{1}{2q}\left[-2+\left(2+2q^2+q^4+2q^6\right)e^{-q^2}
  !>     \right], \\
  !>  k_2(q) &=\frac{1}{2q}\left[-2-q^2+\left(2+3q^2+2q^4\right)e^{-q^2}
  !>     \right], \\
  !>  m_1(q) &=\frac{1}{2q^3}\left[-10+3q^2+\left(10+7q^2+2q^4\right)
  !>     e^{-q^2} \right], \\
  !>  m_2(q) &=\frac{1}{q^3}\left[-2+q^4+\left(2+2q^2\right)e^{-q^2}
  !>     \right].
  !> \f}
  !======================================================================
  Subroutine calculate_finite_range_NM(rho,E,P,K,SMASS,ASS,LASS,KASS,VMASS)
    implicit none
    real(pr), intent(in) :: rho
    real(pr), intent(out) :: E,P,K,SMASS,ASS,LASS,KASS,VMASS
    integer :: ig
    real(pr) :: A0,A1,Bnn,Bnp,B0,W,B,H,M,mui,kF,q,tauc
    real(pr) :: g0,p0,k0,m0,s1,s2,l1,l2,k1,k2,m1,m2
    kF = (1.5_pr*pi**2*rho)**(1._pr/3._pr)
    tauc = 3*kF**2/5._pr
    E=0._pr; P=0._pr; K=0._pr; SMASS=0._pr
    ASS=0._pr; LASS=0._pr; KASS=0._pr; VMASS=0._pr
    do ig = 1,n_g
       mui = mu_g(ig)
       W = W_g(ig)
       B = B_g(ig)
       H = H_g(ig)
       M = M_g(ig)
       A0 =  0.25_pr*(sqrt(pi)*mui)**3*(4*W+2*B-2*H-M)
       A1 =  0.25_pr*(sqrt(pi)*mui)**3*(-2*H-M)
       if(force_is_dme) then
          Bnn = 0._pr
          Bnp = 0._pr
       else
          Bnn = -(W+2*B-H-2*M)/sqrt(pi)
          Bnp =  (H+2*M)/sqrt(pi)
       endif
       B0 = Bnn + Bnp
       q = mui*kF
       g0 = 2/q**3 - 3/q - (2/q**3 - 1/q)*exp(-q**2) + sqrt(pi)*erf(q)
       p0 = -1/q**3 + 1/(2*q) + (1/q**3 + 1/(2*q))*exp(-q**2)
       k0 = -6/q**3 + 2/q + (6/q**3 + 4/q + q)*exp(-q**2)
       m0 = 3*(-2+q**2+exp(-q**2)*(2+q**2))/(2*q**3)
       s1 = 1/q - (1/q + q)*exp(-q**2)
       s2 = 1/q - q -1/q*exp(-q**2)
       l1 = -1/q + (1/q + q + 2*q**3)*exp(-q**2)
       l2 = -1/q - q + (1/q + 2*q)*exp(-q**2)
       k1 = -1/q + (1/q + q + q**3/2._pr + q**5)*exp(-q**2)
       k2 = -1/q - q/2._pr + (1/q + 3/2._pr*q + q**3)*exp(-q**2)
       m1 = (-10+3*q**2+exp(-q**2)*(10+7*q**2+2*q**4))/(2*q**3)
       m2 = (-2+q**4+2*exp(-q**2)*(1+q**2))/(q**3)
       E = E + 0.5_pr*(A0*rho + B0*g0)
       P = P + 0.5_pr*A0*rho**2 + B0*p0*rho
       K = K - 3*B0*k0
       SMASS = SMASS + B0*m0/(2*tauc*hbzero)
       ASS = ASS + 0.5_pr*A1*rho + (Bnn*s1+Bnp*s2)/6._pr
       LASS = LASS + 1.5_pr*A1*rho + (Bnn*l1+Bnp*l2)/6._pr
       KASS = KASS - 2*(Bnn*k1+Bnp*k2)/3._pr
       VMASS = VMASS + (B0*m0+Bnn*m1+Bnp*m2)/(2*tauc*hbzero)
    enddo
  End Subroutine calculate_finite_range_NM

  !======================================================================
  !> Calculates the finite range contribution to the Nuclear Matter
  !> presure.
  !>
  !> See \cite sellahewa2014 for details
  !> @result \f$ P^{\rm NM}(\rho_0) = \sum_{i=1,n_g} \left\{\frac{1}{2}
  !>          A_0^i\rho_0^2 + B_0^ip_0(\mu_ik_F)\rho_0\right\}, \f$
  !======================================================================
  function P_SNM_FR(rho) result(P)
    implicit none
    real(pr), intent(in) :: rho
    real(pr) :: P
    integer :: ig
    real(pr) :: A0,B0,W,B,H,M,mui,kF,q,p0
    kF = (1.5_pr*pi**2*rho)**(1._pr/3._pr)
    P=0._pr
    do ig = 1,n_g
       mui = mu_g(ig)
       W = W_g(ig)
       B = B_g(ig)
       H = H_g(ig)
       M = M_g(ig)
       A0 =  0.25_pr*(sqrt(pi)*mui)**3*(4*W+2*B-2*H-M)
       if(force_is_dme) then
          B0 = 0._pr
       else
          B0 = -(W+2*B-2*H-4*M)/sqrt(pi)
       endif
       q = mui*kF
       p0 = -1/q**3 + 1/(2*q) + (1/q**3 + 1/(2*q))*exp(-q**2)
       P = P + 0.5_pr*A0*rho**2 + B0*p0*rho
    enddo
  end function P_SNM_FR
  !=======================================================================
  !> Calculates INM properties
  !=======================================================================
  Subroutine calculate_NM_properties()
    Implicit None
    Real(pr) :: aRho0Rho0,daRho0Rho0,ddaRho0Rho0,aRho1Rho1,daRho1Rho1,ddaRho1Rho1
    Real(pr) :: aRho0Tau0,daRho0Tau0,ddaRho0Tau0,aRho1Tau1,daRho1Tau1,ddaRho1Tau1
    Real(pr) :: u,tauc,rho_NM2
    Real(pr), Parameter :: c13=1.0_pr/3.0_pr,c23=2.0_pr/3.0_pr
    real(pr) :: E_fr,P_fr,K_fr,SMASS_fr,ASS_fr,LASS_fr,KA_fr,VMASS_fr
    !
    RHO_NM=find_NM_RHOC()
    !
    if(finite_range) then
       call calculate_finite_range_NM(rho_nm,E_fr,P_fr,K_fr,SMASS_fr,ASS_fr,LASS_fr,KA_fr,VMASS_fr)
    else
       E_fr=0_pr; P_fr=0_pr; K_fr=0_pr; SMASS_fr=0_pr
       ASS_fr=0_pr; LASS_fr=0_pr; KA_fr=0._pr; VMASS_fr=0_pr
    end if
    !
    aRho0Rho0=0.50_pr*(aRhoRho+bRhoRho)*mevfm
    aRho1Rho1=0.50_pr*(aRhoRho-bRhoRho)*mevfm
    aRho0Tau0=0.50_pr*(aRhoTau+bRhoTau)*mevfm
    aRho1Tau1=0.50_pr*(aRhoTau-bRhoTau)*mevfm
    daRho0Rho0=0.50_pr*(daRhoRho+dbRhoRho)*mevfm
    daRho1Rho1=0.50_pr*(daRhoRho-dbRhoRho)*mevfm
    daRho0Tau0=0.50_pr*(daRhoTau+dbRhoTau)*mevfm
    daRho1Tau1=0.50_pr*(daRhoTau-dbRhoTau)*mevfm
    ddaRho0Rho0=0.50_pr*(ddaRhoRho+ddbRhoRho)*mevfm
    ddaRho1Rho1=0.50_pr*(ddaRhoRho-ddbRhoRho)*mevfm
    ddaRho0Tau0=0.50_pr*(ddaRhoTau+ddbRhoTau)*mevfm
    ddaRho1Tau1=0.50_pr*(ddaRhoTau-ddbRhoTau)*mevfm
    tauc=CK*RHO_NM**c23; u=(kfconst/mpi)*RHO_NM**c13; rho_NM2=rho_NM**2
    !
    ! Symmetric Nuclear Matter
    E_NM=tauc*hbzero+RHO_NM*(aRho0Rho0+RHO_NM*hRho0Rho0+Crho(0)+RHO_NM**sigma*Cdrho(0)) &
      +tauc*RHO_NM*(aRho0Tau0+RHO_NM*hRho0Tau0+Ctau(0)) + E_fr
    P_NM=c13*RHO_NM**2*((2.0_pr*tauc*hbzero)/RHO_NM+3.0_pr*aRho0Rho0+5.0_pr*tauc*aRho0Tau0 &
      +6.0_pr*RHO_NM*hRho0Rho0+8.0_pr*tauc*RHO_NM*hRho0Tau0+3.0_pr*Crho(0) &
      +3.0_pr*(1+sigma)*RHO_NM**sigma*Cdrho(0)+5.0_pr*tauc*Ctau(0)+u*daRho0Rho0 &
      +u*tauc*daRho0Tau0+u*RHO_NM*dhRho0Rho0+u*tauc*RHO_NM*dhRho0Tau0) + P_fr
    SMASS_NM=1.0_pr+(RHO_NM*(aRho0Tau0+RHO_NM*hRho0Tau0+Ctau(0)))/hbzero + SMASS_fr
    K_NM=9.0_pr*sigma*(1+sigma)*RHO_NM**(1+sigma)*Cdrho(0) &
      +(-2.0_pr*tauc*hbzero+10.0_pr*tauc*RHO_NM*aRho0Tau0+18.0_pr*RHO_NM2*hRho0Rho0 &
      +40.0_pr*tauc*RHO_NM**2*hRho0Tau0+4.0_pr*u*RHO_NM*daRho0Rho0 &
      +RHO_NM*(10.0_pr*tauc*Ctau(0)+u*(8.0_pr*tauc*daRho0Tau0+u*ddaRho0Rho0 &
      +(10.0_pr*RHO_NM*dhRho0Rho0+14.0_pr*tauc*RHO_NM*dhRho0Tau0 &
      +(u*tauc*ddaRho0Tau0+u*RHO_NM*ddhRho0Rho0+u*tauc*RHO_NM*ddhRho0Tau0))))) + k_fr
    !
    ! Asymmetric Nuclear Matter
    ASS_NM=RHO_NM2*hRho1Rho1+RHO_NM*(aRho1Rho1+Crho(1)+RHO_NM**sigma*Cdrho(1)) &
       +(tauc*(5.0_pr*hbzero+RHO_NM*(5.0_pr*aRho0Tau0+15.0_pr*aRho1Tau1+5.0_pr*RHO_NM*hRho0Tau0 &
       +9.0_pr*RHO_NM*hRho1Tau0+5.0_pr*(3.0_pr*RHO_NM*hRho1Tau1+Ctau(0)+3.0_pr*Ctau(1)))))/9.0_pr + ASS_fr
    VMASS_NM=(hbzero+RHO_NM*(aRho0Tau0-aRho1Tau1+RHO_NM*hRho0Tau0-RHO_NM*hRho1Tau1+Ctau(0)-Ctau(1)))/hbzero + VMASS_fr
    LASS_NM=6.0_pr*RHO_NM2*hRho1Rho1+3.0_pr*RHO_NM*(aRho1Rho1+Crho(1)+(1.0_pr+sigma)*RHO_NM**sigma*Cdrho(1)) &
       +u*RHO_NM*daRho1Rho1 +u*RHO_NM2*dhRho1Rho1 &
       +(tauc*(10.0_pr*hbzero+8.0_pr*RHO_NM2*(5.0_pr*hRho0Tau0+9.0_pr*hRho1Tau0+15.0_pr*hRho1Tau1) &
       +25.0_pr*RHO_NM*(aRho0Tau0+3.0_pr*aRho1Tau1+Ctau(0)+3*Ctau(1)) &
       +5.0_pr*u*RHO_NM*(daRho0Tau0+3.0_pr*daRho1Tau1) &
       +u*RHO_NM2*(5.0_pr*dhRho0Tau0+9.0_pr*dhRho1Tau0+15.0_pr*dhRho1Tau1)))/9.0_pr + LASS_fr
    KA_NM=18.0_pr*RHO_NM2*hRho1Rho1+9.0_pr*sigma*(1.0_pr+sigma)*RHO_NM**(1.0_pr+sigma)*Cdrho(1) &
       +4.0_pr*u*RHO_NM*daRho1Rho1 +10.0_pr*u*RHO_NM2*dhRho1Rho1 &
       + u**2*RHO_NM*ddaRho1Rho1+u**2*RHO_NM2*ddhRho1Rho1 &
       +(tauc*(-10.0_pr*hbzero+40.0_pr*RHO_NM2*(5.0_pr*hRho0Tau0+9.0_pr*hRho1Tau0+15.0_pr*hRho1Tau1) &
       +50.0_pr*RHO_NM*(aRho0Tau0+3.0_pr*aRho1Tau1+Ctau(0)+3*Ctau(1)) &
       +40.0_pr*u*RHO_NM*(daRho0Tau0+3.0_pr*daRho1Tau1) &
       +14.0_pr*u*RHO_NM2*(5.0_pr*dhRho0Tau0+9.0_pr*dhRho1Tau0 &
       +15.0_pr*dhRho1Tau1)+5.0_pr*u**2*RHO_NM*(ddaRho0Tau0 &
       +3.0_pr*ddaRho1Tau1)+u**2*RHO_NM2*(5.0_pr*ddhRho0Tau0+9*ddhRho1Tau0+15*ddhRho1Tau1)))/9._pr + KA_fr
    !
  End Subroutine calculate_NM_properties
  !=======================================================================
  !> Find the INM saturation density RHO_NM using the Secant Method
  !=======================================================================
  Real(pr) Function find_NM_RHOC()
    Implicit None
    Integer(pr) :: iter
    Real(pr) :: aRho0Rho0,daRho0Rho0,ddaRho0Rho0,aRho1Rho1,daRho1Rho1,ddaRho1Rho1
    Real(pr) :: aRho0Tau0,daRho0Tau0,ddaRho0Tau0,aRho1Tau1,daRho1Tau1,ddaRho1Tau1
    Real(pr) :: kfconstmpi,u,tauc
    Real(pr) :: rhom0,rhom,rhom2,w,w0,step,energy,P_FR
    Real(pr),Parameter :: c13=1.0_pr/3.0_pr,c23=2.0_pr/3.0_pr
    !
    kfconstmpi=kfconst/mpi; step=-0.0010_pr; iter=0
    ! initial point
    rhom=0.170_pr; tauc=CK*rhom**c23; u=kfconstmpi*rhom**c13; rhom2=rhom**2
    !
    Call calculate_U_parameters(rhom,rhom,tauc*rhom,tauc*rhom,0.0_pr,0.0_pr)
    !
    if(finite_range) then
       P_FR = P_SNM_FR(rhom)
    else
       P_FR = 0
    endif
    aRho0Rho0=0.50_pr*(aRhoRho+bRhoRho)*mevfm; daRho0Rho0=0.50_pr*(daRhoRho+dbRhoRho)*mevfm
    aRho0Tau0=0.50_pr*(aRhoTau+bRhoTau)*mevfm; daRho0Tau0=0.50_pr*(daRhoTau+dbRhoTau)*mevfm
    w0=c13*rhom2*((2.0_pr*tauc*hbzero)/rhom+3.0_pr*aRho0Rho0+5.0_pr*tauc*aRho0Tau0 &
      +6.0_pr*rhom*hRho0Rho0+8.0_pr*tauc*rhom*hRho0Tau0+3.0_pr*Crho(0) &
      +3.0_pr*(1.0_pr+sigma)*rhom**sigma*Cdrho(0)+5.0_pr*tauc*Ctau(0)+u*daRho0Rho0 &
      +u*tauc*daRho0Tau0+u*rhom*dhRho0Rho0+u*tauc*rhom*dhRho0Tau0) + P_FR
    rhom0=rhom; rhom=rhom+step
    !
    ! secant method
    Do While(Abs(step).Ge.eps*100.0_pr)
       iter=iter+1
       tauc=CK*rhom**c23; u=kfconstmpi*rhom**c13; rhom2=rhom**2
       !
       Call calculate_U_parameters(rhom,rhom,tauc*rhom,tauc*rhom,0.0_pr,0.0_pr)
       !
       if(finite_range) then
          P_FR = P_SNM_FR(rhom)
       else
          P_FR = 0
       endif
       aRho0Rho0=0.50_pr*(aRhoRho+bRhoRho)*mevfm; daRho0Rho0=0.50_pr*(daRhoRho+dbRhoRho)*mevfm
       aRho0Tau0=0.50_pr*(aRhoTau+bRhoTau)*mevfm; daRho0Tau0=0.50_pr*(daRhoTau+dbRhoTau)*mevfm
       w=c13*rhom2*((2.0_pr*tauc*hbzero)/rhom+3.0_pr*aRho0Rho0+5.0_pr*tauc*aRho0Tau0 &
         +6.0_pr*rhom*hRho0Rho0+8.0_pr*tauc*rhom*hRho0Tau0+3.0_pr*Crho(0) &
         +3.0_pr*(1.0_pr+sigma)*rhom**sigma*Cdrho(0)+5.0_pr*tauc*Ctau(0)+u*daRho0Rho0 &
         +u*tauc*daRho0Tau0+u*rhom*dhRho0Rho0+u*tauc*rhom*dhRho0Tau0) + P_FR
       step=-w*(rhom-rhom0)/(w-w0)
       rhom0=rhom; w0=w; rhom=rhom+step
       If(iter.Gt.100) Stop 'STOP(In find_NM_RHOC)'
       !energy=tauc*hbzero+rhom*(aRho0Rho0+rhom*hRho0Rho0+Crho(0)+rhom**sigma*Cdrho(0)) &
       ! +tauc*rhom*(aRho0Tau0+rhom*hRho0Tau0+Ctau(0))
       !Write(6,'(a,15(1pg12.4))') ' rhom0,rhom,step,e,w=',rhom0,rhom,step,energy,w
    End Do
    find_NM_RHOC=rhom
  End Function find_NM_RHOC
  !=======================================================================
  !
  !=======================================================================
  Subroutine C_from_t()
    !--------------------------------------------------------------------------------
    ! C- from (t,x)-
    !--------------------------------------------------------------------------------
    Implicit None
    Crho(0)  =   3.0_pr/8.0_pr  * t0
    Cdrho(0) =  (1.0_pr/16.0_pr)* t3
    Crho(1)  = -(1.0_pr/4.0_pr) * t0*(0.50_pr+x0)
    Cdrho(1) = -(1.0_pr/24.0_pr)* t3*(0.50_pr+x3)
    Ctau(0)  =  (3.0_pr/16.0_pr)* t1+(1.0_pr/4.0_pr)*t2*(5.0_pr/4.0_pr+x2)
    Ctau(1)  = -(1.0_pr/8.0_pr) * t1*(0.5+x1)+(1.0_pr/8.0_pr)*t2*(0.50_pr+x2)
    CrDr(0)  =  (1.0_pr/16.0_pr)* t2*(5.0_pr/4.0_pr+x2)-(9.0_pr/64.0_pr)*t1
    CrDr(1)  =  (3.0_pr/32.0_pr)* t1*(0.5+x1)+(1.0_pr/32.0_pr)*t2*(0.50_pr+x2)
    CJ(0)    = -(1.0_pr/16.0_pr)*(t1*(2.0_pr*x1-1.0_pr)+t2*(2.0_pr*x2+1)-5*te-15*to)
    CJ(1)    = -(1.0_pr/16.0_pr)*(t2 -t1 + 5.0_pr*te -5.0_pr*to )
    CrdJ(0)  = -b4-(0.50_pr)*b4p
    CrdJ(1)  = -0.50_pr*b4p
  End Subroutine C_from_t
  !=======================================================================
  !
  !=======================================================================
  Subroutine t_from_C()
    !--------------------------------------------------------------------------------
    ! (t,x)- from C-
    !--------------------------------------------------------------------------------
    Implicit None
    t0  =  (8.0_pr/3)*Crho(0)
    t1  =  4.0_pr/3.0_pr*(Ctau(0)-4.0_pr*CrDr(0))
    t2  =  4.0_pr/3.0_pr*(3.0_pr*Ctau(0)-6.0_pr*Ctau(1)+4.0_pr*CrDr(0)-8.0_pr*CrDr(1))
    t3  =  16.0_pr*Cdrho(0)
    x0  = -0.50_pr*(3.0_pr*Crho(1)/Crho(0)+1.0_pr)
    x1  =  2.0_pr*(-Ctau(0)-3.0_pr*Ctau(1)+4.0_pr*CrDr(0)+12.0_pr*CrDr(1))/t1/3.0_pr
    x2  = -2.0_pr*(3.0_pr*Ctau(0)-15.0_pr*Ctau(1)+4.0_pr*CrDr(0)-20.0_pr*CrDr(1))/t2/3.0_pr
    x3  = -0.50_pr*(3.0_pr*Cdrho(1)/Cdrho(0)+1.0_pr)
    b4  =  CrdJ(1)-CrdJ(0)
    b4p = -2.0_pr*CrdJ(1)
    te  = (4.0_pr/15.0_pr)*(3.0_pr*CJ(0)-9.0_pr*CJ(1)-4.0_pr*CrDr(0)+12.0_pr*CrDr(1)-2.0_pr*Ctau(0)+6.0_pr*Ctau(1))
    to  = (4.0_pr/15.0_pr)*(3.0_pr*CJ(0)+3.0_pr*CJ(1)+4.0_pr*CrDr(0)+4.0_pr*CrDr(1))
    if(finite_range.and..not.force_is_dme) then
       x0 = zero
       x1 = zero
       x2 = zero
    endif
    ! Set public variables
    t0_pub = t0; t1_pub  = t1;  t2_pub = t2; t3_pub = t3
    x0_pub = x0; x1_pub  = x1;  x2_pub = x2; x3_pub = x3
    b4_pub = b4; b4p_pub = b4p; te_pub = te; to_pub = to
  End Subroutine t_from_C
  !=======================================================================
  !
  !=======================================================================
  Subroutine Make_Parameter_Free_Useful_Combinations()
    !--------------------------------------------------------------------------------
    ! Make Useful combinations
    !--------------------------------------------------------------------------------
    Implicit None
    !
    If(dmeorder.Ge.0) Then
       !
       mpi2=mpi**2
       gA2=gA**2; gA4=gA2**2; gA6=gA2**3;
       fpi2=fpi**2; fpi4=fpi2**2;
       CHartree =mevfm*(3.0_pr*gA2)/(32.0_pr*fpi4*Pi**2)
       h0mpi6=197.30_pr*(mpi**6)*(3.0_pr*gA*gA)/(32.0_pr*fpi**4*Pi**2)
       h0mpi6c1=h0mpi6*c1;         h0mpi6c3=h0mpi6*c3
       !
       h0mpi6NM=197.30_pr*(3.0_pr*(mpi**3)*gA2)/(64.0_pr*fpi**4*Sqrt(Pi))
       h0mpi6c1NM=h0mpi6NM*c1;     h0mpi6c3NM=h0mpi6NM*c3
       !
       A3_1=42.7132145164590_pr;   A3_2=0.670441422115440_pr; A3_3=0.0525713896514650_pr;
       A3_4=0.0012545731701320_pr; A3_5=5.81008627207380_pr
       b3_1=3.0809379008590_pr;    b3_2=0.905186811964580_pr; b3_3=0.474514509597610_pr;
       b3_4=0.228138177966090_pr;  b3_5=1.66931540698090_pr;
       !
       A1_1=2.5000830618386_pr;    A1_2=0.619542286897850_pr; A1_3=0.169682589033730_pr;
       A1_4=0.0276112113725470_pr; A1_5=0.00108164458809540_pr
       b1_1=1.75854210706510_pr;   b1_2=0.88882524524657_pr;  b1_3=0.46377235143756_pr;
       b1_4=0.247665887704790_pr;  b1_5=0.132222413002680_pr
       !
       Call load_tables()
       !
    End If
    !
  End Subroutine Make_Parameter_Free_Useful_Combinations
  !=======================================================================
  !
  !=======================================================================
  Elemental Function Vexternal(t,x,y,z)
    !
    Implicit None
    Integer(ipr), Intent(in) :: t  !! isospin index: 0=isoscalar, 1=isovector
    Real(pr), Intent(in) :: x,y,z  !! position in cartesian basis
    Real(pr) :: Vexternal
    !
    Vexternal = 0.0_pr
    !
  End Function Vexternal
  !=======================================================================
  !
  !=======================================================================
End Module UNEDF

