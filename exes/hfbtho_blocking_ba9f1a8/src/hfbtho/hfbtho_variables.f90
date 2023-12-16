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
!                         GLOBAL VARIABLES PACKAGE                     !
!                                                                      !
! ==================================================================== !

!----------------------------------------------------------------------
!> This module defines most of the global variables used throughout
!> the program. In particular, it contains most global arrays
!> pertaining to either configuration or coordinate space.
!----------------------------------------------------------------------
Module HFBTHO
  Use HFBTHO_VERSION
  Use UNEDF
  !
  ! Input for HFBiterations
  Integer(ipr)  :: n00_INI,iLST_INI,inin_INI,icou_INI
  Integer(ipr)  :: npr_INI(3),kindhfb_INI
  Integer(ipr)  :: keyblo1_INI,keyblo2_INI,IDEBUG_INI
  Integer(ipr)  :: ngh_INI,ngl_INI,nleg_INI,nstate_INI
  Real(pr)      :: b0_INI,bz_INI=0.0_pr,bp_INI=0.0_pr,q_INI
  Character(30) :: skyrme_INI
  Real(pr)      :: pwi_INI,V0n_INI,V0p_INI,cpv1_INI,epsi_INI
  Logical       :: basis_HFODD_INI,Add_Pairing_INI,Print_HFBTHO_Namelist_INI,DO_FITT_INI
  ! Output for regression optimization
  Real(pr)  :: efit_0
  Real(pr), Dimension(0:1) :: efit_rhorho,efit_rhorhoD,efit_rhotau,efit_rhoDrho
  Real(pr), Dimension(0:1) :: efit_rhonablaJ,efit_JJ,efitV0,dfitV0,efV_0
  ! Serial output (1:on/0:off)
  Integer(ipr) :: debug_solver !< debug_solver>=1: times of important subroutines are printed
                               !< debug_solver>=2: tests of HFB matrix in hfbdiag() and breakdown of energy contributions
                               !< debug_solver>=3: HO basis quantum numbers
                               !< debug_solver>=4: HO basis wave functions on Gauss quadrature mesh
  Logical :: DO_FITT
  ! For loop over used particle types. For normal nuclei min=1, max=2. For n droplets min=max=1.
  Integer(ipr) :: itmin,itmax,irestart
  ! Frequent constants
  Real(pr) :: PI,ffdef3,ffdef4,ffdef5,ffdef6,ffdef7
  ! Single constants
  Real(pr) :: bet,beta0,q,bp,bpp,bz,hom,hb0,b0,etot,coex,cex,ty20,vin,rin,ain,  &
              qin,pwi,si,siold,epsi,xmix,xmix0,xmax,alst,clst,sklst,alphi,amas, &
              skass,varmas,v0ws,akv,hqc,amu,amn,r0,r00,r02,r04,decay,rmm3,amm3, &
              bmm3,cmm3,chargee2,EBASECUT,rho_c,cdef,b1_0,b2_0,b3_0,b4_0,wct_gogny,hb0n,hb0p
  Integer(ipr) :: lin,lwin,lwou,lplo,lwel,lres,icstr,icou,ncut,iLST1,iLST, &
                  mini,maxi,iiter,inin,nzm,nrm,icacou,iqrpa,icacoupj,icahartree,nlm,  &
                  nb,nt,n00,itass,kindhfb,iappend,iError_in_HO,iError_in_THO,    &
                  ierest,esu,nstate,iter_fail
  Integer(ipr), Parameter :: iter_fail_max = 400 ! 40 x n, with n the number of HFB iterations
  Integer(ipr), Parameter :: n00max=50
  ! Results
  Integer(ipr), Parameter :: ieresu=50,ieresl=20,ieresj=50,ieresbl=6
  Integer(ipr), Parameter :: ieres=ieresu+ieresl+ieresj+ieresbl
  Real(pr) :: eres(ieres)
  Character(13) :: ereslbl(2)
  Character(2) :: nucname
  Real(pr)    :: eresu(ieresu),eresl(ieresl),eresbl(ieresbl),eresj(ieresj)
  Character(13) :: hlabels(ieres+3)
  ! Common small arrays
  Real(pr) :: alast(2),ala(2),ala1(2),tz(2),ass(2),drhoi(2),del(2),vso(2),r0v(2),av(2),   &
              rso(2),aso(2),Sumnz(2),Dispersion(2),v2min(2),v2minv(2),rms(3),ept(3),q2(3),&
              frdel(2),frept(3),Dnfactor(3),varmasNZ(2),pjmassnz(2),&
              codel(2), coept(3)
  Integer(ipr) :: npr(3),inz(2),ldel(2),nk(2),itbl(2),kbl(2),tpar(2),ipbl(2),nbl(2),ibbl(2), &
                  klmax(2),inner(2),iasswrong(3),lcc
  ! Lipkin-Nogami
  Real(pr) :: ala2(2),etr(3),ssln(3,2),Geff(2)
  ! Blocking
  !   - bloblo ...: index ib of the Omega block where the bloking is made
  !   - blo123 ...: index k0 of blocked state within the Omega block
  !   - blok1k2 ..: absolute position of blocked state within all Omega blocks
  !   - bloqpdif .: difference in energy between current qp energy and lowest one in the block
  Real(pr) :: pwiblo=2.0_pr, eqpmin(2)=0.0_pr
  Logical :: odd_noBlock = .False.
  Logical, Dimension(2) :: blocking_never_done
  Integer(ipr), Dimension(2) :: blocking_mode
  Integer(ipr) :: bloall; Parameter(bloall=200)
  Integer(ipr), Dimension(0:bloall,2) :: bloblo,blo123=0,blok1k2=0
  Real(pr), Dimension(0:bloall,2) :: bloqpdif
  Integer(ipr) :: iparenti(2),keyblo(3),nkblo_INI(2,5),nkblo(2,5)=0
  Integer(ipr) :: blocross(2),blomax(2),blo123d(2),blok1k2d(2),blocanon(2)
  ! Logical and character variables
  Character(1) :: tq,tp(2),tl(0:20),tis(2)
  Character(30) :: skyrme
  Character(8) :: tit(2)
  Character(7) :: protn(2)
  Data  protn/'neutron','proton '/
  ! Allocatable arrays
  Character(13), Allocatable  :: tb(:)
  Character(25), Allocatable  :: txb(:)
  Real(pr), Allocatable, Target :: rk(:,:),ak(:,:),hh0(:,:),de0(:,:),ddc(:,:,:),ddc1(:,:,:),&
                                   qh(:,:),qh1(:,:),ql(:,:,:),ql1(:,:,:),ek(:,:),dk(:,:),   &
                                   vk(:,:),vk1(:,:),uk(:,:),hfb1(:,:),vkmax(:,:)
  Real(pr), Allocatable :: fdsx(:),fdsy(:),fdsy1(:),fdsy2(:),fdsy3(:),fspb0(:),fspc0(:),fspd0(:),  &
                           fspb1(:),fspc1(:),fspd1(:),fspb2(:),fspc2(:),fspd2(:),fspb3(:),fspc3(:),&
                           fspd3(:),fak(:),fi(:),sq(:),sqi(:),wf(:),wfi(:),rkass(:,:)
  Integer(ipr), Allocatable :: id(:),ia(:),ikb(:),ipb(:),nz(:),nr(:),nl(:),ns(:),npar(:), &
                               ka(:,:),kd(:,:),numax(:,:),iv(:), lcanon(:,:)
  Integer(ipr), Allocatable :: nrr(:),nll(:),nss(:),noo(:),nzzx(:),nzz(:,:)
  Real(pr), Allocatable  :: AN(:),ANK(:),PFIU(:),PFID(:)
  Real(pr), Allocatable  :: FIU(:),FID(:),FIUR(:),FIDR(:)
  Real(pr), Allocatable  :: FIUD2N(:),FIDD2N(:),FIUZ(:),FIDZ(:)
  ! Constraints
  Integer(ipr), Parameter :: lambdaMax=8,ncons_max=20
  Integer(ipr) :: numberCons
  Integer(ipr), Dimension(1:ncons_max) :: multLambda
  Real(pr), Dimension(0:lambdaMax) :: q_units,multLag,multRequested
  Real(pr), Dimension(0:lambdaMax,1:3) :: qmoment
  Real(pr), Allocatable :: gaussian_neck(:)
  Real(pr) :: neckLag,neckValue,neckRequested,mixing_neck=0.5_pr
  ! Temperature
  Logical :: switch_on_temperature
  Real(pr) :: temper
  Real(pr), Dimension(3) :: entropy
  Real(pr), Allocatable, Target :: fn_T(:),fp_T(:)
  ! optimization arrays
  Real(pr), Allocatable  :: QHLA_opt(:,:),FI1R_opt(:,:),FI1Z_opt(:,:),FI2D_opt(:,:),y_opt(:)
  ! Arrays depending on mesh points
  Integer(ipr)  :: ngh,ngl,nleg,nghl,nbx,ntx,nzx,nrx,nlx,ndx,ndx2,ndxs,nqx,nnx,nrlx,nttx,nox
  Integer(ipr)  :: nhfbqx,nb2x,nhfbx,nkx,nzrlx,iqqmax
  Real(pr), Allocatable :: COEF00(:,:,:),HERFAC(:)
  Real(pr), Allocatable :: xh(:),wh(:),xl(:),sxl(:),wl(:),xleg(:),wleg(:),vc(:,:)
  Real(pr), Allocatable :: vhbn(:),vn(:),vrn(:),vzn(:),vdn(:),vsn(:),dvn(:)
  Real(pr), Allocatable :: vhbp(:),vp(:),vrp(:),vzp(:),vdp(:),vsp(:),dvp(:)
  Real(pr), Allocatable :: vSZFIn(:),vSFIZn(:),vSRFIn(:),vSFIRn(:)
  Real(pr), Allocatable :: vSZFIp(:),vSFIZp(:),vSRFIp(:),vSFIRp(:)
  Real(pr), Allocatable :: qfield(:,:)
  Real(pr), Allocatable, Target :: aka(:,:),ro(:,:),tau(:,:),dro(:,:),dj(:,:), &
                                   SZFI(:,:),SFIZ(:,:),SRFI(:,:),SFIR(:,:),NABLAR(:,:),NABLAZ(:,:)
  Real(pr), Allocatable  :: fl(:),fli(:),fh(:),fd(:),fp1(:),fp2(:),fp3(:),fp4(:),fp5(:),fp6(:),  &
                            fs1(:),fs2(:),fs3(:),fs4(:),fs5(:),fs6(:),wdcor(:),wdcori(:),cou(:)
  Real(pr), Allocatable  :: vDHartree(:,:),vhart00(:,:),vhart01(:,:),vhart11(:,:)
  ! PAV Projection
  Logical :: statistical_PNP = .False.
  Integer(ipr) :: keypjn,keypjp,ilpjn,ilpjp,ilpj,ilpjmax,ilpj2,ilnqx,ilnghl
  Integer(ipr) :: iproj,npr1pj,npr2pj
  Real(pr) :: rehfbcan,ehfb,retotpj,depnp
  Complex(pr) :: onei=(0.0_pr,1.0_pr)
  Complex(pr), Allocatable, Target :: phypj(:),sinphy(:),exp1iphy(:),exp1iphym(:),exp2iphy(:), &
                                      exp2iphym(:),coupj(:,:),ropj(:,:,:),taupj(:,:,:),dropj(:,:,:), &
                                      djpj(:,:,:),akapj(:,:,:),pjk(:,:),SZFIpj(:,:,:),SFIZpj(:,:,:), &
                                      SRFIpj(:,:,:),SFIRpj(:,:,:),epj(:,:),ddepj(:,:,:),cpj(:,:,:), &
                                      ypj(:,:,:),rpj(:,:,:)
  Real(pr), Allocatable, Target :: xg_norm(:,:),xg_phase(:,:)
  Real(pr) :: polem(2),polep(2)
  ! CMC
  Integer(ipr) :: ICMinput
  Real(pr) :: ECMHFB(3),ECMPAV(3)
  ! CRC
  Integer(ipr) :: ICRinput
  Real(pr) :: DEROT(3),SQUJ(3),CRAN(3),ERIGHFB(3)
  ! hfbdiagonal
  Real(pr), Allocatable :: erhfb(:),drhfb(:),erhfb1(:),drhfb1(:)
  Real(pr), Allocatable :: hfb(:,:),zhfb(:),evvk(:),hfbcan(:,:),evvkcan(:)
  ! Jason: def derived types
  Type :: ptr_to_2darray
     Real(pr),    Dimension(:,:),Allocatable :: arr
  End Type ptr_to_2darray
  Type :: ptr_to_array
     Real(pr),    Dimension(:),Allocatable :: arr
  End Type ptr_to_array
  Type :: ptr_to_iarray
     Integer(ipr), Dimension(:),Allocatable :: arr
  End Type ptr_to_iarray
  ! Jason: use derived types
  Type(ptr_to_2darray), Allocatable :: allhfb(:)
  Type(ptr_to_array),   Allocatable :: allevvk(:),allalwork(:)
  Type(ptr_to_iarray),  Allocatable :: alllwork(:),allISUPPZ(:)
  Integer(ipr),        Allocatable :: allibro(:),allIALWORK(:),allILWORK(:)
  Integer(ipr) :: oldnb
  Real(pr) :: cutoff_tol=1.d-6
  ! Broyden
  Character(1) ::  bbroyden
  Integer(ipr) :: nbroyden=7
  Real(pr) :: alphamix=0.70_pr
  Integer(ipr) :: nhhdim,nhhdim2,nhhdim3,nhhdim4,ialwork,ilwork
  Real(pr), Allocatable, Target  :: brout(:),brin(:)
  Real(pr), Allocatable :: alwork(:)
  Integer(ipr), Allocatable :: lwork(:)
  ! cm
  Real(pr)  :: facECM=1.0_pr
  ! fission fragments properties
  Real(pr) :: tz_fragments(3)
  ! new keys
  Logical :: set_pairing,basis_HFODD,Parity,Parity_INI
  Logical :: neck_constraints = .False.
  Logical :: Print_Screen=.True.
  Logical :: write_hel=.True.
  Logical :: full_HFB_matrix=.False.
  Logical :: Add_Pairing,Print_HFBTHO_Namelist
  Integer(ipr) :: MAX_ITER_INI,keypjn_INI,keypjp_INI,iproj_INI,npr1pj_INI,npr2pj_INI
  ! Eqp U,V
  Integer(ipr) :: nuv,nqp
  Real(pr), Allocatable, Target :: RVqpN(:),RVqpP(:),RUqpN(:),RUqpP(:),REqpN(:),REqpP(:)
  Integer(ipr), Allocatable, Target :: KpwiN(:),KpwiP(:),KqpN(:),KqpP(:)
  ! Pairing regularization
  Real(pr), Allocatable :: MEFFn(:),MEFFp(:)
  Real(pr), Allocatable :: geff_inv(:,:)
  ! error indicator
  Integer(ipr)  :: ierror_flag=0
  Character(60) :: ierror_info(0:11)
  ! --------------------------------
  ! Petar, Restoration of symmetries
  ! --------------------------------
  ! Number constants
  Real(pr), Parameter:: third=one/3_pr
  Complex(pr), Parameter:: czero=Dcmplx(zero,zero),cone=Dcmplx(one,zero),cimag=Dcmplx(zero,one)
  ! Symmetry restoration meshes
  Integer(ipr):: maxbet,maxphi,maxphi_eff,maxj,jjstep,maxN,maxP
  Real(pr) :: phicyl_integration_step,integration_prefactor
  Real(pr), Allocatable :: betabs(:),betaweight(:),phiabs(:)
  Complex(pr), Allocatable:: rotated_overlap(:,:,:),rotation_matrix(:,:),inverse_rotation_matrix(:,:)
  Complex(pr), Allocatable:: ephi(:),ephic(:),ephicN(:,:),ephicP(:,:)
  Complex*32, Allocatable:: detR(:)
  Integer(ipr), Allocatable:: beta_active(:)
  ! Cylindrical coordinate phi
  Real(pr), Allocatable :: phicyl(:)
  Integer(ipr), Parameter:: NGPHI=20
  ! (rperp, z) x (iphicyl) mesh
  Integer(ipr), Allocatable:: ihil_convert(:),iphicyl_convert(:),ihil_iphicyl_convert(:,:)
  ! V and U Matrices
  Integer(ipr) kdim(4)
  Type :: ptr_to_cmplx2darray
     Complex(pr), Dimension(:,:), Allocatable :: arr
  End Type ptr_to_cmplx2darray
!  Type(ptr_to_cmplx2darray), Allocatable :: Vmatrix(:),Umatrix(:)
  Complex(pr), Allocatable :: VmatrixN1(:,:),VmatrixN2(:,:),UmatrixN1(:,:),UmatrixN2(:,:),&
                              VmatrixP1(:,:),VmatrixP2(:,:),UmatrixP1(:,:),UmatrixP2(:,:)
  ! Rotated densities
  Complex(pr), Allocatable :: rotated_ro(:,:,:),rotated_tau(:,:,:),rotated_dj(:,:,:),rotated_dro(:,:,:),rotated_sfiz(:,:,:),&
                              rotated_sfir(:,:,:),rotated_srfi(:,:,:),rotated_szfi(:,:,:),rotated_szz(:,:,:),rotated_srz(:,:,:),&
                              rotated_srr(:,:,:),rotated_szr(:,:,:),rotated_sfifi(:,:,:),rotated_aka(:,:,:),ro_projected(:,:)
  Real(pr), Allocatable :: ro_normalization(:,:)
  Real(pr) piu(2)
  Type(ptr_to_cmplx2darray), Allocatable ::  rotated_density(:,:,:),rotated_kappa(:,:,:)
  Integer(ipr):: nt_pair
  Integer(ipr), Allocatable:: iosc1_contributing(:,:),iosc2_contributing(:,:),nt_contributing(:),iosc1_pair(:),iosc2_pair(:),itiphi_pair1(:),itiphi_pair2(:)
  ! Rotated and projected energies
!  Complex(pr), Allocatable :: ekinN_rotated(:,:,:),ekinP_rotated(:,:,:),ecodi_rotated(:,:,:),ecoex_rotated(:,:,:),EVOL_rho_tau_rotated(:,:,:),EVOL_rho_rho_rotated(:,:,:),&
!                              ESURF_rho_drho_rotated(:,:,:),ETENS_rotated(:,:,:),ESO_rho_nablaj_rotated(:,:,:),eptN_rotated(:,:,:),eptP_rotated(:,:,:),&
!                              xn1_rotated(:,:,:),xn2_rotated(:,:,:),rms1_rotated(:,:,:),rms2_rotated(:,:,:),delN_rotated(:,:,:),delP_rotated(:,:,:)
  Complex(pr), Allocatable :: all_overlaps(:),all_overlaps_gthr(:),all_energies(:),all_energies_gthr(:)
  Real(pr), Allocatable:: projected_overlap(:),projected_ekinN(:),projected_ekinP(:),projected_ecodi(:),projected_ecoex(:),&
                          projected_EVOL_rho_tau(:),projected_EVOL_rho_rho(:),projected_ESURF_rho_drho(:),projected_ETENS(:),&
                          projected_ESO_rho_nablaj(:),projected_eptN(:),projected_eptP(:), projected_xn1(:),projected_xn2(:),projected_rms1(:),projected_rms2(:),&
                          projected_delN(:),projected_delP(:),projected_NP(:,:,:),projected_NP_norm(:)
  ! Simplex basis
  Integer(ipr), Allocatable :: nz_sim(:),nr_sim(:),nl_sim(:)
  ! Rperp nad z in fm, with (ih,il) indices
  Real(pr), Allocatable:: xl_ihil(:),xh_ihil(:)
  ! Timing variables
  Integer(ipr) tstart,tfinish,clock_rate
  Real(pr) time_total,time_ibet(6)
  ! Coulomb
  Complex(pr), Allocatable::  cou_projected(:,:)
  ! MPI
  Integer(ipr) beta_size,betaphi_size
  !---------------------------------------------------------------------
  ! Namelists
  !---------------------------------------------------------------------
  Logical      :: add_initial_pairing, set_temperature, compatibility_HFODD, force_parity, &
                  user_pairing, collective_inertia, automatic_basis, fission_fragments, &
                  pairing_regularization, localization_functions, set_neck_constrain, &
                  filter, fragment_properties, include_3N_force
  Integer(ipr) :: number_of_shells, proton_number, neutron_number, type_of_calculation, &
                  number_iterations, type_of_coulomb, restart_file, projection_is_on,   &
                  gauge_points, delta_Z, delta_N, switch_to_THO, number_Gauss,          &
                  number_Laguerre, number_Legendre, number_states, print_time,          &
                  PNP_is_on, number_of_gauge_points, delta_neutrons, delta_protons,     &
                  AMP_is_on, number_of_rotational_angles, maximal_angular_momentum
  Integer(ipr) :: proton_blocking(1:5), neutron_blocking(1:5), lambda_values(1:lambdaMax), &
                                                               lambda_active(1:lambdaMax)
  Real(pr)     :: oscillator_length, basis_deformation, beta2_deformation, beta3_deformation, &
                  beta4_deformation, accuracy, temperature, vpair_n, vpair_p, pairing_cutoff, &
                  pairing_feature, neck_value, real_N, real_Z
  Real(pr)     :: expectation_values(1:lambdaMax)
  Character(Len=30) :: functional
  Namelist /HFBTHO_GENERAL/ number_of_shells, oscillator_length, basis_deformation, &
                            proton_number, neutron_number, type_of_calculation
  Namelist /HFBTHO_INITIAL/ beta2_deformation, beta3_deformation, beta4_deformation
  Namelist /HFBTHO_ITERATIONS/ number_iterations, accuracy, restart_file
  Namelist /HFBTHO_FUNCTIONAL/ functional, add_initial_pairing, type_of_coulomb, include_3N_force
  Namelist /HFBTHO_PAIRING/ user_pairing, vpair_n, vpair_p, pairing_cutoff, pairing_feature
  Namelist /HFBTHO_CONSTRAINTS/ lambda_values, lambda_active, expectation_values
  Namelist /HFBTHO_BLOCKING/ proton_blocking, neutron_blocking
  Namelist /HFBTHO_PROJECTION/ switch_to_THO,projection_is_on,gauge_points,delta_Z,delta_N
  Namelist /HFBTHO_TEMPERATURE/ set_temperature, temperature
  Namelist /HFBTHO_FEATURES/ collective_inertia, fission_fragments, pairing_regularization, &
                             automatic_basis, localization_functions
  Namelist /HFBTHO_TDDFT/ filter, fragment_properties, real_Z, real_N
  Namelist /HFBTHO_NECK/ set_neck_constrain, neck_value
  Namelist /HFBTHO_DEBUG/ number_Gauss, number_Laguerre, number_Legendre, &
                          compatibility_HFODD, number_states, force_parity, &
                          write_hel, print_time
  Namelist /HFBTHO_RESTORATION/ PNP_is_on, number_of_gauge_points, delta_neutrons, delta_protons, &
                                AMP_is_on, number_of_rotational_angles, maximal_angular_momentum
  !---------------------------------------------------------------------
  ! Mass tables
  !---------------------------------------------------------------------
#if(DO_MASSTABLE==1)
  Integer(ipr) :: ierr_mpi, mpi_taskid, mpi_size, nRows, iRow, nRows_task, disp_int, disp_real, task_error
  Integer(ipr) :: n_int_masstable_in = 2, n_real_masstable_in = 2
  Integer(ipr) :: n_int_masstable_out = 3, n_real_masstable_out = 6
  Integer(ipr) :: Z_chain,N_chain,A_chain,Calc_counter
  Integer(ipr), Allocatable :: Z_masstable(:),N_masstable(:)
  Integer(ipr), Allocatable :: Z_out(:),N_out(:),ierrors_out(:)
  Real(pr) :: beta_deformation,Q20
  Real(pr), Allocatable :: Q20_in(:),E_HFB_out(:),Q20Z_out(:),Q20N_out(:),Q20T_out(:),beta_in(:)
  Real(pr), Allocatable :: Q20_masstable(:),beta_masstable(:)
  Character(Len=6) :: ID_string
  Character(Len=7) :: row_string
  !---------------------------------------------------------------------
  ! Potential energy surfaces (includes sub-communicators dedicated to HFB calculations)
  !---------------------------------------------------------------------
#elif(DO_PES==1)
  Integer(ipr) :: ierr_mpi, mpi_taskid, mpi_size, nRows, iRow, nRows_task, disp_int, disp_real, task_error
  Integer(ipr) :: n_int_masstable_in = 2, n_real_masstable_in = 3 ! minimum number: it is going to be set in read_hfbtho_PES
  Integer(ipr) :: n_int_masstable_out = 3, n_real_masstable_out = 6
  Integer(ipr) :: Z_chain,N_chain,A_chain,Calc_counter
  Integer(ipr) :: npoints,ndefs
  Integer(ipr) :: HFB_cores = 1
  Integer(ipr), Allocatable :: Z_PES(:),N_PES(:),lambda_PES(:)
  Real(pr) :: beta_deformation
  Real(pr), Allocatable :: bet2_PES(:), bet3_PES(:), bet4_PES(:), Q_PES(:,:)
  Character(Len=6) :: ID_string
  Character(Len=7) :: row_string
  !---------------------------------------------------------------------
  ! Driplines (includes sub-communicators dedicated to HFB calculations)
  !---------------------------------------------------------------------
#elif(DRIP_LINES==1)
  Integer(ipr) :: ierr_mpi, mpi_taskid, mpi_size, nRows, iRow, nRows_task, disp_int, disp_real, task_error
  Integer(ipr) :: n_int_masstable_in = 3, n_real_masstable_in = 0
  Integer(ipr) :: n_int_masstable_out = 3, n_real_masstable_out = 6
  Integer(ipr) :: number_deformations=11,i_deformation
  Integer(ipr) :: Z_chain,N_chain,A_chain,Calc_counter
  Integer(ipr), Allocatable :: Z_stable_line(:), N_stable_line(:), direction_sl(:)
  Integer(ipr), Allocatable :: Z_out(:),N_out(:),ierrors_out(:)
  Real(pr) :: beta_deformation,Q20
  Real(pr) :: beta_step,separation_2N,Minimum_Energy,Minimum_Energy_Prev,Energy_chain
  Real(pr), Allocatable :: Q20_in(:),E_HFB_out(:),Q20Z_out(:),Q20N_out(:),Q20T_out(:),beta_in(:)
  Real(pr), Allocatable :: Energy_chain_gthr(:)
  Character(Len=6) :: ID_string
  Character(Len=3) :: team_string
  Character(Len=3) :: counter_string
  Character(Len=9) :: direction_str
  Character(Len=11) :: row_string
#else
  !---------------------------------------------------------------------
  ! Others
  !---------------------------------------------------------------------
  Integer(ipr) :: ierr_mpi, mpi_taskid, mpi_size, nRows, iRow, nRows_task, disp_int, disp_real, task_error
  Integer(ipr) :: n_int_masstable_in = 0,n_real_masstable_in = 0
  Character(Len=6) :: ID_string
#endif
  !---------------------------------------------------------------------
  ! MPI
  !---------------------------------------------------------------------
  Integer(ipr) :: team_color=0,number_teams=1,team_size=1,team_rank=0
  Integer(ipr) :: do_print=1
#if(USE_MPI==2)
  Integer(ipr), Allocatable :: task_error_gthr(:)
  Integer(ipr), Allocatable :: vector_sizes(:),vector_int_mpi(:)
  Integer(ipr), Allocatable :: vector_sizes_gthr(:),vector_sizes_int_gthr(:),vector_sizes_real_gthr(:)
  Integer(ipr), Allocatable :: vector_int_gthr(:),vector_disp_int_gthr(:),vector_disp_real_gthr(:)
  Real(pr), Allocatable :: vector_real_mpi(:),vector_real_gthr(:)
  Logical, Allocatable :: vector_log_mpi(:)
  ! MPI Subcommunicators
  Integer(ipr) :: COMM_team, COMM_leaders
#endif
  !---------------------------------------------------------------------
  !Read_Functional variables
  !---------------------------------------------------------------------
  Integer(ipr) :: n_func_param
  Real(pr), allocatable :: functional_vector(:)
End Module HFBTHO
