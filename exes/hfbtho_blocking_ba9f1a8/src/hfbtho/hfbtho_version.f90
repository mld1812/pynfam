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
!                         VERSION HISTORY                              !
!                                                                      !
! ==================================================================== !

!----------------------------------------------------------------------
!> This module defines the version of \theCode and contains a
!> summary changelog of past versions (up to 201)
!----------------------------------------------------------------------
Module HFBTHO_VERSION
  Implicit None
  Character(6) :: Version='3.00'
  !--------------------------------------------------------------------------------------
  ! Version History
  !--------------------------------------------------------------------------------------
  ! ver#201: Reinstated the DME module
  ! ver#200d: fixed bug in gfv; improved legibility and accuracy of coulom and coulom1
  ! ver#200c: added LLNL release number
  ! ver#200b: added module linear_algebra, analyzing THO, formatted output, fixed bug
  !           in calculation of entropy
  ! ver#200a: Restored option to compute all blocking configurations within given energy
  !           window, removed spurious output printing
  ! ver#199: added a few input options and a compatibility mode with HFODD. Release
  !          candidate before publication, last bug to fix is OpenMP in hfbdiag()
  ! ver#142: removed module pairing, module UNEDF, spurious preprocessing options
  !          for publication purposes; used the Lahey compiler to identify a few bugs,
  !          added routine check_consistency(), module HFBTHO_gauss, improved file
  !          handling in inout() and fixed bug in multipole moments
  ! ver#141: Reinstated the THO module
  ! ver#140: added a namelist for debugging purposes, added OpenMP in coulom()
  ! ver#139a: fixed bug in readjustment of constraint at finite T, cleaned up output
  !           system, fixed bug in calculation of Coulomb in parity breaking mode
  ! ver#139: added temperature
  ! ver#138: added new module HFBTHO_utilities to improve portability
  ! ver#137: new system of inputs based on several namelists contained in one unique
  !          file called hfbtho_NAMELIST.dat, multiple constraints up and running
  ! ver#136: ANL OpenMP optimizations included, and code clean-up
  ! ver#135: Tested with all previous versions of HFBTHO down to 101, full compatibility
  !          achieved
  ! ver#134: Begining of work toward publication
  ! ver#133: Single-file HFBTHO
  ! ver#130: tho.dat mdifications due to blocking, error indicator introduced
  ! ver#129: Even-even tested and equivalent with ptho101spt15sp.f90 used in ANL fit
  ! ver#128: EQP,U,V and their dimentions NUV,NEQ required for qrpa incoporated
  !          permanently in HFBTHO substituting old arrays eqp and uv
  ! ver#127: For easy development the module split in different F90 files which are
  !          invoked using INCLUDE statements (remove that when developement is over)
  !          Optimized qrpa_DENSIT_PLUS and qrpa_GAMDEL to be twice faster
  ! ver#126: Cleaning, optimizing, and isolating THO stuff
  !                parity good:   Time per iteration: 3.841 seconds
  !                parity broken: Time per iteration: 9.933 seconds
  !          HFB+THO tested and works in both parity regimes. iserial removed and
  !          substituted with Print_Screen, i.e, record results only when Nsh>0
  ! ver#125: Implemented and tested reflection symmetry as option. If parity is broken,
  !          computer time per iteration is almost 5 times bigger:
  !                parity good:   Time per iteration:  4.5276 seconds
  !                parity broken: Time per iteration: 19.4646 seconds
  !                difference in total energy:  0.001  keV
  ! ver#124: Rewrite to prepare for breaking reflection symmetry
  !          Preprocessing directives included: #ifndef hide_qrpa, hide_tho, hide_dme
  !          For preprocessing one needs -Dhide_qrpa -Dhide_tho -Dhide_dme
  !          If no preprocessing or -Uhide_qrpa -Uhide_tho -Uhide_dme then
  !          all modules are included
  ! ver#123: Playground for QRPA calculations [16/12/2010]:
  !          The changes are:
  !          - all variables in Module HFBTHO are Public
  !          - include_qrpa=0 is added to Module HFBTHO with asssociated
  !            declarations eventualy used by qrpa
  !          - Subroutine ByNucleus moved to PTHO_PROGRAM where is its place
  !            and it should be done long ago. Call Do_QRPA() is used only there.
  !          - So if the program is compiled with -Uhide_qrpa one can use Do_QRPA()
  !            to do qrpa calculations.
  ! ver#123: Fixed crash after iterations limit.
  !          Tested against anl version hfbtho101spt15.f90 - itterations go differently
  !          but the final results are identical.
  ! ver#122: (MK) Added CExPar for coulomb exchange. Parameter read from UNEDF module
  ! ver#121: (MK) Added possibility to use zero particle number for droplet calculations
  ! ver#120: (MK) Added external field, and all channels to direct Hartree. e^2 for Coulomb
  !          now read from the UNEDF module. Direct Hartree now always calculated based on
  !          module function regardless of the value of DMEorder parameter
  ! ver#117: Direct Hartree added when DME_order>-1
  ! ver#115: (MK) added use_cm_cor variable to hb0 calculation and (nabla rho)^2 terms to
  !          the calculate_U_parameters function calls
  ! ver#114: Name list, new tho.dat file, proton/neutron fields, confirms all results of
  !          recent published version after ANL optimization ptho101b_last_tested.f90,
  !          public/Public variables
  ! ver#113: Cleaning
  ! ver#112: No parameter functions
  ! ver#111: Main program detached from the file as PTHO_MAIN_PROGRAM.f90 which will not be
  !          versioned. ptho becomes jus a HFBTHO module.  Pairing constants V0(2),V1(2)
  !          replaced by CpV0(0:1), CpV1(0:1) coming as public from defined in UNEDF module
  !          Removed dalf and ippforce form the pairing. For compatibility, ippforce
  !          stays in the input file by now but the kind of pairing is given by CpV1 only
  !          Dropped corrections 'ecmcpavpj', 'erotcorrection' which should be added later.
  !          For compatibility inputfile stays the same. Added IDEUB.
  !          THO part in 'densit' (not densitpj), 'gamdel' commented HO/THO for speed
  ! ver#109: All public variables, expectpj works with a jump: not clear how UNEDF can work
  !          with complex numbers, just skip this part by now bu write results data
  !          New thodefh(iw1)
  ! ver#108: Removed all programs not used in ver#107
  !          expect contains a key DO_FITT:
  !            =0 calculare energy, delta, def  & rms only
  !            =1 the same+all integrals for the regression optimization
  !          V0,V1 pairing constant separated for neutrons and protons: v0(2),v1(2)
  !          HFBTHO collected in MODULE HFBTHO
  !          KOP3 removed
  ! ver#107: Towards UNEDF: complete rewrite based on Marcus to include N2LO
  !          LN for ZR110 at prolate solution with SLY4, mixed pairing and tensor terms:
  !          -agreement with previouse Skyrme implemetation to the last significant digit
  !          -agreement with previouse LO+LDA implemetation to the last significant digit
  !          -agreement with previouse LO+CB  implemetation to the last significant digit
  ! ver#106: Towards UNEDF: the standard functional rewritten in terms of UNEDF U-amplitudes
  !          The assumption U=U(tau_0,Delta rho_0,rho_0,rho_1) becomes possible after
  !          adding Nabla rho_ij terms  (STANDARD FUNCTIONAL ONLY)
  ! ver#105: Towards UNEDF: the standard functional rewritten in terns of UNEDF U-amplitudes
  !          The assumption U=U(rho_0,rho_1) becomes possible after adding Delta rho_ij terms
  ! ver#104: Broyden improved with linear search at negative curvature
  !          Implemented Agumented Lagrangian method for constraint calculations
  !          Manual blocking included and tested, key: manualBlocking
  ! ver#103: From this version on-no more support for VAP (VAP completely removed)
  !          The whole program in terms of C-parameters (including tenzor terms)
  ! ver#102: The whole program in terms of C-parameters (without tenzor terms)
  ! ver#101: Optimization in terms of nuclear matter: 'FITS' regime
  ! ver#100: Toward isovector pairing following Sagawa and Yamagami
  ! ver# 99: Subroutine HFBiterations. The isotopic line in tho.dat removed.
  !          Subroutines byNucleus, byConstraint, FitPairing, HFBTHO_HFODD isolated
  !          at the end and could be ported if necessary. skyrme='FITS' assumes the skyrme
  !          parameters as explicitely given. -N00 supresses completely the output and only
  !          hodef.dat and thodef.dat are charged (if iserial=0 even these files are supressed)
  !          HFBTHO_HFODD updated (think further about a constraint in Q2 terms)
  ! ver# 98: INOUT modified and added interface subroutine HFBTHO_HFODD
  !          Pairing fitted with MIX/(LN-NOLN) for SLY4,SKP,SKM* forces
  ! ver# 97: Corrected blocking candidates criteria
  ! ver# 96: Extended so term W0,W1 and SKLY4T forces
  ! ver# 95: Pairing regularization removed, linear HFBDIAG mixing for Lambda
  !          when blocking, DSYEV replaced by the faster DSYEVD, hfbdiag caculates
  !          canonical basis only at the last hfbdiag iteration, expect optimized
  !          new subroutines HFBiterations, FitPairing, byConstraint
  !          Work around a bug in LAPAK related to DSYEVD
  ! ver# 94: Misprints
  ! ver# 93: Removed hh and de matrices and related manipulations
  !          Broyden_min now escape maximum and inflex points
  ! ver# 92: Bug in blocking while no pairing cleaned
  ! ver# 91: BLAST  & LAPACK diagonalization
  ! ver# 90: If applying Broyden method to matrix elements then
  !          at 20 shells the total number of matrix elements
  !          is 2x2x65307=261228 or about 2.1 Mb and if one keeps
  !          8 iterations it will be  about 17 Mb-not too much
  !          This is the only way one can mix Lipkin-Nogami
  !
  !          If one uses the potentials at 30x30 grid points
  !          the numbers are 8x2x900=14400 or 115 Kb and if one keeps
  !          8 iterations it will be about 1 Mb-much better
  !          but Lipkin-Nogami is out of this scheme (?!)
  !          If one uses densities at 30x30 grid points they
  !          are 9x2x900-almost like the potential case.
  !          (sent to George)
  ! ver# 89: Reduced printout (no anymore lprinter)
  !          LN del+ala2 printed during the iterations
  !          Strength in the initial constraint calculations reduced to 0.3, requested
  !          deformation+/-0.3, untill si<1.1
  !          If too slow convergence (1000 iterations) and Lambda>0 interrupt iterations
  !          Odd nucleus right away from the even-even (even) solution
  !          When even solution missing/corrupt (even at inin<0) calculate it first
  !          and then odd one
  ! ver# 88: Synchronization for the parallel run
  !          FileLabel subroutine added.
  !          Modified inin control
  !            inin<0: Always start from a file if it exists, not corrupted and correct
  !                    otherwise inin=iabs(inin) and start from scratch
  !            inin>0: Always start from scratch
  !          SCRATCH calculations start with initial 20 constrained iteration if constrain
  !          is not requested (icstr=0). When icstr#0 standard constraint calculations.
  !          BY CHAIN calculations temporary removed due to blocking complications
  ! ver# 87: Approximate Blocking keeping time-reversal symmetry to PNP PAV
  ! ver# 86: Approximate Blocking keeping time-reversal symmetry and tested agings HFODD
  ! ver# 85: LN in canonical basis. Benchmark to HFODD
  ! ver# 84: Testing HFODD LN again HFB-HO
  ! ver# 83: Cleaning, SKM* mixed volume LN pairing fitted
  ! ver# 82: As 81 but prepared for jaguar
  ! ver# 81: Pairing regularization/renormalization. PAV done with unprojected v_k
  !          V0(Nsh=20,pwi=50) fitted for SLY4,SKP,Renormalized,Regularized,Mixed,Volume
  !          Removed delta and gamdel0 completely
  ! ver# 80: Accuracy for large number of shells increased by the number of gauss points
  !          Gaussian points now calculated
  !          Initial guest now from deformed Wood-Saxon
  !          Initial run now starts with requested shell number n00
  ! ver# 79: Cranking rotational correction implemented:
  !          Printed to screen, thoout.dat, hodef.dat and thodef.dat
  !          but not added to the energy
  ! ver# 78: Full CM correction implemented in HFB  & HFB(PAV)
  !          Printed to screen, thoout.dat, hodef.dat and thodef.dat
  !          but not added to the energy
  !          NB : ilpjnp(2) removed
  ! ver# 77: Automated Blocking:
  !          First is calculated N,Z without blocking, remembered *.hel *.tel files and
  !           determined blocking candidates according to pwablo criteria.
  !           if N(Z) is odd we have neutron(proton) blocking candidates.
  !           if both, N and Z, are odd we have both, neutron and proton, block.candidates
  !          Then we block state after state among the blocking candidates and calculate
  !           starting from the recorded unblocked (N,Z) solution
  !           if only N ( or Z) is odd then all neutron (or proton) blocking candidates
  !           are calculated
  !           if both N and Z are odd all pairs of proton and neutron blocking candidates
  !           are calculated                                               (PAV  &LN unclear)
  ! ver# 76: Manual Blocking for a particular state in a particular block
  !          overlap criteria used to avoid the level crossing problem    (PAV  &LN unclear)
  ! ver# 75: Manual Blocking for the minimal qpe within a given block.
  ! ver# 74: Bulgac procedure .. not done
  ! ver# 73: Cleanup, introducing the cpc notations, beyond unit circle removed (MARK 1)
  ! ver# 72: PNP: still valid version for integration over the unit circle
  ! ver# 71: PNP: towards beyond unit circle integration
  ! ver# 70: PNP: detached neutron from proton projection
  ! ver# 69: PNP: quadrupole constraint fixed to converge
  ! ver# 68: equivalent to var#67
  ! ver# 67: PNP: corrections to the tensor term and initial dumping factor
  ! ver# 66: PNP: V0 fitted to PLN energies at Sn126 for SKP mixed and volume at Nsh=20,HO
  ! ver# 65: deformed HO basis implemented and tested
  ! ver# 64: fixed byChain to go not 2 beyond the forced break
  ! ver# 63: iasswrong(3) fixed for correct multiprocessor run
  ! ver# 61/62: pthotop for Cheetah added at the end
  ! ver# 60: Thodef.dat header line fixed (added U:). Fixed P/N in ByChain calculations
  !          The 'Stop' is removed from  mishmatch conditions with alternative to use old one.
  !          Consistent pairing for SLY4 and SKP forces
  ! ver# 59: LST modified to accept negative aa-values. SLY4 And SKP with consistent
  !          (high densities regime) pairing for all cases. Old asymptotic prescription
  !          is used in the case of Mishmatch asymptotic. Temporary,still new SLY4
  !          pairing constants are commented.   ! ver# 58: Proton line in byChain
  !          calculations goes vertically. In case of wrong
  !          asymptotic parameter 'kindhfb' is recorded as 'kindhfb+100' in thodef.dat
  !          file where the results for this nucleus are substituted with HO results
  !          Only Nsh=20 pairing constants are already fitted to the higher density
  !          asymptotic prescription which is already enforced. (Temporary, still old
  !          SLY4 pairing constants are in the code).
  ! ver# 57: Partially refitted pairing constants according to the new asymptotic
  !          prescription. NB! old constants are still for the SLY4 force and not all of
  !          the cases with SKP are fitted. Old ass. regime is temporary enforced.
  ! ver# 56: Code optimization and checks
  ! ver# 55: Tensor term J.J implemented and tested
  ! ver# 54: Back to *.hel *.tel files; Including new hodef.dat like thodef.dat file
  ! ver# 53: LST is choosing the higher density in the asymptotic region
  !--------------------------------------------------------------------------------------
End Module HFBTHO_VERSION
