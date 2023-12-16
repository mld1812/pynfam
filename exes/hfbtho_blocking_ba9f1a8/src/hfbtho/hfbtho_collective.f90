!!***********************************************************************
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
!                      COLLECTIVE INERTIA PACKAGE                      !
!                                                                      !
! ==================================================================== !

!-------------------------------------------------------------------
!> This module calculates the collective inertia mass tensor at the
!> perturbative cranking approximation. Both the ATDHFB and GCM
!> tensors are computed. The GCM formula foor the zero-point energy
!> correction is used with either the ATDHFB or GCM mass tensor.
!>
!> @author
!> Nicolas Schunck
!-------------------------------------------------------------------
!  Subroutines: - calculate_collective_mass
!               - define_collective_space
!               - energy_moments(it)
!               - qp_basis_F20(nd,n_qp,Umatr,Vmatr,FOPQP)
!               - qp_basis_F11
!               - S2VCPY(A,N,NA,V,M,MV)
!               - V2SCPY(V,M,MV,A,N,NA)
!               - INVMUL(A,B,XA,N,NX)
!               - INVERT(A,Y,N,NP)
!               - ludcmp(a,n,np,indx,d)
!               - lubksb(a,n,np,indx,b)
!  Functions: - TRACEM(A,N,NX)
!----------------------------------------------------------------------!
Module HFBTHO_collective

  Use HFBTHO_utilities
  Use HFBTHO
  Use HFBTHO_multipole_moments

  Implicit None

  Integer(ipr), Parameter :: NDEFOR=lambdaMax, NDCOMP=(NDEFOR*(NDEFOR+1))/2, NDEMOM=3

  Integer(ipr), PRIVATE, SAVE :: debug_inertia = 0 !< Internal flag for debugging
  Real(pr), PRIVATE, SAVE :: COVMET(NDEFOR,NDEFOR,1:3) !< GCM metric; third index: 1 (neutrons), 2 (protons), 3 (total)
  Real(pr), PRIVATE, SAVE :: EMOMET(NDEFOR,NDEFOR,NDEMOM) !< Energy moments of order k=1,.., NDEMOM

  Integer(ipr), PUBLIC, SAVE :: NOCOMP
  Integer(ipr), PUBLIC, SAVE :: NOMULT
  Integer(ipr), PUBLIC, SAVE :: IIPAIR
  Integer(ipr), PUBLIC, SAVE :: INDMAS(NDCOMP,2)

  Real(pr), PUBLIC, SAVE :: SK1(NDCOMP,1:3) !< Energy moments of order 1; second index: 1 (neutrons), 2 (protons), 3 (total)
  Real(pr), PUBLIC, SAVE :: SK2(NDCOMP,1:3) !< Energy moments of order 2; second index: 1 (neutrons), 2 (protons), 3 (total)
  Real(pr), PUBLIC, SAVE :: SK3(NDCOMP,1:3) !< Energy moments of order 3; second index: 1 (neutrons), 2 (protons), 3 (total)

  Real(pr), PUBLIC, SAVE :: ATDMAS(NDCOMP,1:3) !< ATDHFB mass tensor: second index 1 for neutrons, 2 for protons, 3 for total
  Real(pr), PUBLIC, SAVE :: GCMMAS(NDCOMP,1:3) !< GCM mass tensor: second index 1 for neutrons, 2 for protons, 3 for total
  Real(pr), PUBLIC, SAVE :: E0_ATD(1:3) !< ATDHFB zero-point energy correction; first index: 1 (neutrons), 2 (protons), 3 (total)
  Real(pr), PUBLIC, SAVE :: E0_GCM(1:3) !< GCM zero-point energy correction; first index:  1 (neutrons), 2 (protons), 3 (total)

Contains

  !=======================================================================
  !> This routine computes the collective inertia mass tensor and the
  !> zero-point energy vibrational corrections for the current set of
  !> constraints. It is a copy-paste of the same routine in HFODD, only
  !> adapted to the particular case of HFBTHO.
  !=======================================================================
  Subroutine calculate_collective_mass()
    Implicit None
    !
    Integer(ipr) :: it
    Real(pr), Dimension(NDEFOR,NDEFOR) :: AUX1,AUX2,AUX3
    !
    Call define_collective_space()
    !
    Do it=itmin,itmax
       !
       ! Calculate energy moments (EMOMET)
       Call energy_moments(it)
       !
       ! Saving all the moments M1, M2 and M3 for the current isospin
       AUX1=EMOMET(:,:,1)
       Call S2VCPY(AUX1,NOMULT,NDEFOR,SK1(:,it),NOCOMP,NDCOMP)
       AUX2=EMOMET(:,:,2)
       Call S2VCPY(AUX2,NOMULT,NDEFOR,SK2(:,it),NOCOMP,NDCOMP)
       AUX3=EMOMET(:,:,3)
       Call S2VCPY(AUX3,NOMULT,NDEFOR,SK3(:,it),NOCOMP,NDCOMP)
       !
       !-------------------------------------------------------------------!
       !> The GCM metric is
       !> \f[
       !>   \mathbf{G} = \frac{1}{2} [\mathbf{M}^{(1)}]^{-1} \mathbf{M}^{(2)} [\mathbf{M}^{(1)}]^{-1}.
       !>   \label{eq:G_local}
       !> \f]
       !-------------------------------------------------------------------!
       AUX1=EMOMET(:,:,1)
       AUX2=EMOMET(:,:,2)
       Call INVMUL(AUX1,AUX2,AUX3,NOMULT,NDEFOR)
       AUX1=0.50_pr*AUX1
       !
       ! Backup of the metric gamma for current isospin
       COVMET(:,:,it)=AUX1
       !
       !-------------------------------------------------------------------!
       !          I.  ATDHF (Cranking approximation only)                  !
       !-------------------------------------------------------------------!
       !
       AUX1=EMOMET(:,:,1)
       AUX3=EMOMET(:,:,3)
       !
       !> The ATDHFB mass is given by
       !> \f[
       !>   M_{\mathrm{ATD}} = \hbar^{2}[\mathbf{M}^{(1)}]^{-1} \mathbf{M}^{(3)} [\mathbf{M}^{(1)}]^{-1}.
       !> \f]
       Call INVMUL(AUX1,AUX3,AUX2,NOMULT,NDEFOR)
       Call S2VCPY(AUX1,NOMULT,NDEFOR,ATDMAS(:,it),NOCOMP,NDCOMP)
       !
       ! Make B = M^-1
       Call V2SCPY(ATDMAS(:,it),NOCOMP,NDCOMP,AUX1,NOMULT,NDEFOR)
       Call INVERT(AUX1,AUX3,NOMULT,NDEFOR)
       !
       !> The zero-point energy (ZPE) correction for the ATDHFB mass is
       !> \f[
       !>   \epsilon_{\mathrm{ATD}} = -\frac{\hbar^2}{2}\mathrm{Tr}\left(\mathbf{B}_{\mathrm{ATD}}\mathbf{G}\right)
       !> \f]
       !> where \f$\mathbf{B} = \mathbf{M}^{-1} \f$
       AUX2=COVMET(:,:,it)
       AUX3=MATMUL(AUX2,AUX1)
       E0_ATD(it)=0.50_pr*TRACEM(AUX3,NOMULT,NDEFOR)
       !
       !-------------------------------------------------------------------!
       !          II.  GCM                                                 !
       !-------------------------------------------------------------------!
       !
       !> The GCM inertia tensor is given by
       !> \f[
       !>   \mathbf{B}_{\mathrm{GCM}} = \frac{1}{4}\mathbf{G}^{-1} [\mathbf{M}^{(1)}]^{-1} \mathbf{G}^{-1}.
       !> \f]
       !> and the GCM collective mass is simply the inverse of this matrix,
       !> \f$ \mathbf{M}_{\mathrm{GCM}} = \mathbf{B}^{-1}_{\mathrm{GCM}} \f$
       AUX1=EMOMET(:,:,1)
       Call INVERT(AUX1,AUX3,NOMULT,NDEFOR)
       AUX2=COVMET(:,:,it)
       Call INVERT(AUX2,AUX3,NOMULT,NDEFOR)
       AUX2=MATMUL(AUX2,MATMUL(AUX1,AUX2))
       AUX2=0.25_pr*AUX2
       !
       ! Save M_GCM = B_GCM^(-1)
       Call INVERT(AUX2,AUX3,NOMULT,NDEFOR)
       Call S2VCPY(AUX2,NOMULT,NDEFOR,GCMMAS(:,it),NOCOMP,NDCOMP)
       !
       ! The ZPE is hbar^2/2 Tr(gamma * B_GCM)
       Call V2SCPY(GCMMAS(:,it),NOCOMP,NDCOMP,AUX1,NOMULT,NDEFOR)
       Call INVERT(AUX1,AUX3,NOMULT,NDEFOR)
       AUX2=COVMET(:,:,it)
       AUX3=MATMUL(AUX2,AUX1)
       E0_GCM(it) = 0.50_pr*TRACEM(AUX3,NOMULT,NDEFOR)
       !
    End Do
    !
    !-------------------------------------------------------------------!
    !          III.  CALCULATIONS FOR NEUTRONS+PROTONS                  !
    !-------------------------------------------------------------------!
    !
    ! The total moments are the sum of the moments for each isospin
    SK1(:,3)=SK1(:,1)+SK1(:,2) ! M1_total
    SK2(:,3)=SK2(:,1)+SK2(:,2) ! M2_total
    Call V2SCPY(SK1(:,3),NOCOMP,NDCOMP,AUX1,NOMULT,NDEFOR)
    Call V2SCPY(SK2(:,3),NOCOMP,NDCOMP,AUX2,NOMULT,NDEFOR)
    ! Total metric
    Call INVMUL(AUX1,AUX2,AUX3,NOMULT,NDEFOR)
    COVMET(:,:,3)=0.50_pr*AUX1
    !
    ! ATDHF
    SK1(:,3)=SK1(:,1)+SK1(:,2) ! M1_total
    SK3(:,3)=SK3(:,1)+SK3(:,2) ! M3_total
    Call V2SCPY(SK1(:,3),NOCOMP,NDCOMP,AUX1,NOMULT,NDEFOR)
    Call V2SCPY(SK3(:,3),NOCOMP,NDCOMP,AUX3,NOMULT,NDEFOR)
    !
    ! Total ATDHFB mass from the sum of the moments
    Call INVMUL(AUX1,AUX3,AUX2,NOMULT,NDEFOR)
    Call S2VCPY(AUX1,NOMULT,NDEFOR,ATDMAS(:,3),NOCOMP,NDCOMP)
    ! Make B = M^-1 and get trace of 1/2 * (gamma * B)
    Call V2SCPY(ATDMAS(:,3),NOCOMP,NDCOMP,AUX1,NOMULT,NDEFOR)
    Call INVERT(AUX1,AUX3,NOMULT,NDEFOR)
    ! Get trace of 1/2 * (gamma * B)
    AUX2=COVMET(:,:,3)
    AUX3=MATMUL(AUX2,AUX1)
    E0_ATD(3)=0.50_pr*TRACEM(AUX3,NOMULT,NDEFOR)
    !
    ! GCM
    Call V2SCPY(SK1(:,3),NOCOMP,NDCOMP,AUX1,NOMULT,NDEFOR)
    Call INVERT(AUX1,AUX3,NOMULT,NDEFOR)
    AUX2=COVMET(:,:,3)
    Call INVERT(AUX2,AUX3,NOMULT,NDEFOR)
    AUX2=MATMUL(AUX2,MATMUL(AUX1,AUX2))
    AUX2=0.25_pr*AUX2
    Call INVERT(AUX2,AUX3,NOMULT,NDEFOR)
    Call S2VCPY(AUX2,NOMULT,NDEFOR,GCMMAS(:,3),NOCOMP,NDCOMP)
    ! Make B = M^-1 and get trace of 1/2 * (gamma * B)
    Call V2SCPY(GCMMAS(:,3),NOCOMP,NDCOMP,AUX1,NOMULT,NDEFOR)
    Call INVERT(AUX1,AUX3,NOMULT,NDEFOR)
    AUX2=COVMET(:,:,3)
    AUX3=MATMUL(AUX2,AUX1)
    E0_GCM(3)=0.50_pr*TRACEM(AUX3,NOMULT,NDEFOR)
    !
  End Subroutine calculate_collective_mass
  !=======================================================================
  !> This subroutine defines the collective space, that is, the number
  !> of active multipole moments used as constraints (Q10 does not count)
  !> and the total number of independent matrix elements of the
  !> collective mass tensor.
  !=======================================================================
  Subroutine define_collective_space()
    Implicit None
    !
    Integer(ipr) :: icons,jcons,NIMULT,lambda,lambdb

    ! Number of LAMBDA>1 constraints
    NIMULT=0
    ! Number of different components with LAMBDA>1
    NOCOMP=0
    !
    Do icons=1,numberCons
       lambda=multLambda(icons)
       If(lambda.Ne.1) Then
          NIMULT=NIMULT+1
          Do jcons=1,icons
             lambdb=multLambda(jcons)
             If(lambdb.Ne.1) Then
                NOCOMP=NOCOMP+1
                INDMAS(NOCOMP,1)=lambda
                INDMAS(NOCOMP,2)=lambdb
             End If
          End Do
       End If
    End Do
    ! Number of multipoles
    NOMULT=NIMULT
    !
  End Subroutine define_collective_space
  !=======================================================================
  !> Energy moments \f$ \mathbf{M}^{(K)} \f$ are defined by
  !> \f[
  !>    M_{ab}^{(K)} = \mathfrak{Re} \sum_{\mu\nu}
  !>      \frac{\langle \mu\nu | \hat{Q}_{a} | 0\rangle
  !>            \langle 0 | \hat{Q}_{b} | \mu\nu\rangle}{(E_{\mu}+E_{\nu})^{K}}.
  !> \f]
  !> In principle the summation over indices \f$\mu, \nu \f$ are
  !> unrestricted. In practice, we run only over half of the qp
  !> states, those with positive simplex. An extra factor 2 is thus
  !> added in routine INVEMO.
  !=======================================================================
  Subroutine energy_moments(it)
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Integer(ipr), Intent(in) :: it !< isospin index (1=neutrons, 2=protons)
    !
    Logical :: debug_robledo, debug_BIII
    Integer(ipr) :: ib,i_uv,nd,nd2,n_qp,k,n1,n2,kk
    Integer(ipr) :: i_coll,j_coll,N_COMP,lambda,lambdb,lambda_test,lambdb_test,mu,nu
    Real(pr) :: FACTOR,Factor_q2,Factor_q3,ELPROD,SUMMEQ,qvala,qvalb,dd1,dd2
    Real(pr), Pointer     :: EqpPo(:),VqpPo(:),UqpPo(:)
    Integer(ipr), Pointer :: KpwiPo(:),KqpPo(:)
    Real(pr), Allocatable :: Eqp(:),rho(:,:)
    Real(pr), Allocatable :: FOPQPA(:,:),FOPQPB(:,:)
    Real(pr), Allocatable :: Umatr(:,:),Vmatr(:,:)
    Real(pr), Allocatable :: multMatElems(:)
    Real(pr), Dimension(NDEMOM) :: aux
    !
    If(it.Eq.1) Then
       EqpPo=>REqpN; VqpPo=>RVqpN; UqpPo=>RUqpN; KpwiPo=>KpwiN; KqpPo=>KqpN
    Else
       EqpPo=>REqpP; VqpPo=>RVqpP; UqpPo=>RUqpP; KpwiPo=>KpwiP; KqpPo=>KqpP
    End If
    !
    EMOMET = zero
    !
    ! Loop over the collective variables
    N_COMP=0
    Do i_coll=1,NOMULT
       Do j_coll=1,i_coll
          !
          N_COMP=N_COMP+1
          lambda=INDMAS(N_COMP,1)
          lambdb=INDMAS(N_COMP,2)
          !
          qvala=zero; qvalb=zero; dd1=zero; dd2=zero
          !
          i_uv=0
          Do ib = 1,nb
             ! Get characteristics of the current block: size of HO basis and number of qp
             nd=id(ib); nd2=nd*nd; n_qp=kd(ib,it)
             If(n_qp.Gt.0) Then
                !
                Allocate(FOPQPA(n_qp,n_qp),FOPQPB(n_qp,n_qp))
                Allocate(multMatElems(nd2))
                ! Get quasiparticle energies and U and V for the current block
                Allocate(Eqp(n_qp)); Eqp=zero
                Allocate(Umatr(nd,n_qp)); Umatr=zero
                Allocate(Vmatr(nd,n_qp)); Vmatr=zero
                Do k=1,n_qp
                   kk=KqpPo(ka(ib,it)+k); Eqp(k)=EqpPo(kk)
                End Do
                Do k=1,n_qp
                   Do n2=1,nd
                      Vmatr(n2,k)=VqpPo(KpwiPo(ka(ib,it)+k)+n2)
                      Umatr(n2,k)=UqpPo(KpwiPo(ka(ib,it)+k)+n2)
                   End Do
                End Do
                !
                ! In configuration space, test that Tr(rho) = N (contained in dd1 and dd2), and that
                ! Tr(Q*rho) = Q (in qvala and qvalb)
                If(debug_inertia.Ge.2) Then
                   Allocate(rho(nd,nd))
                   Call dgemm('n','t',nd,nd,n_qp,one,Vmatr,nd,Vmatr,nd,zero,rho,nd)
                   lambda_test=2
                   Call moments_expectation(lambda_test,it,ib,qvala,rho,dd2,multMatElems)
                   If(ib.Eq.nb) Write(6,'("ib=",i2," lambda=",i2," qvala=",f14.7," dd2=",f14.7)') ib,lambda,qvala,dd2
                   lambdb_test=3
                   Call moments_expectation(lambdb_test,it,ib,qvalb,rho,dd1,multMatElems)
                   If(ib.Eq.nb) Write(6,'("ib=",i2," lambdb=",i2," qvalb=",f14.7," dd1=",f14.7)') ib,lambdb,qvalb,dd1
                End If
                !
                ! Calculate matrix 'f' of constraint operator in qp basis for variables lambda and lambdb
                FOPQPA=zero; multMatElems=zero
                Call moments_computeField(lambda,ib,multMatElems)
                Call qp_basis_F20(nd,n_qp,Umatr,Vmatr,FOPQPA,multMatElems)
                FOPQPB=zero; multMatElems=zero
                Call moments_computeField(lambdb,ib,multMatElems)
                Call qp_basis_F20(nd,n_qp,Umatr,Vmatr,FOPQPB,multMatElems)
                !
                aux=zero
                Do mu=1,n_qp
                   Do nu=1,n_qp
                      SUMMEQ=Eqp(mu)+Eqp(nu)
                      ELPROD=FOPQPA(mu,nu)*FOPQPB(mu,nu)
                      Do k=1,NDEMOM
                         aux(k) = aux(k) + ELPROD/SUMMEQ**k
                      End Do
                   End Do ! End nu
                End Do ! End mu
                !
                Deallocate(Eqp)
                Deallocate(FOPQPA,FOPQPB)
                Deallocate(multMatElems)
                Deallocate(Umatr,Vmatr)
                If(debug_inertia.Ge.2) Deallocate(rho)
                !
                ! Save inverse energy moments
                FACTOR=two
                Do K=1,NDEMOM
                   EMOMET(i_coll,j_coll,K)=EMOMET(i_coll,j_coll,K)+FACTOR*aux(K)
                End Do
                !
             End If ! n_qp > 0
             !
          End Do ! end of loop over ib
          !
       End Do ! end of loop over j_coll
    End Do ! end of loop over i_coll
    Do K=1,NDEMOM
       Do i_coll=1,NOMULT
          Do j_coll=1,i_coll
             EMOMET(j_coll,i_coll,K)=EMOMET(i_coll,j_coll,K)
          End Do
       End Do
    End Do
    !
    If(debug_inertia.Ge.1) Then
       Factor_q2 = one; Factor_q3 = one
       debug_robledo = .False.
       debug_BIII    = .False.
       If(debug_robledo) Then
          Write(6,'(/,/,"Comparing with Robledo results")')
          Factor_q2 = 0.5_pr * 100.0_pr
          Factor_q3 = sqrt(16.0_pr*Atan(1.0_pr)/7.0_pr) * 1000.0_pr
       End If
       If(debug_BIII) Then
          Write(6,'(/,/,"Comparing with BIII results")')
          Factor_q2 = 100.0_pr
          Factor_q3 = sqrt(16.0_pr*Atan(1.0_pr)/7.0_pr) * 1000.0_pr
       End If
       Write(6,'(/,/,"Assume 2D calculation with Q2 and Q3 and multiply moments by appropriate factors it=",i2)') it
       write(6,'("Factor_q2=",f10.5," Factor_q3=",f10.5)') Factor_q2,Factor_q3
       N_COMP=0
       Do i_coll=1,NOMULT
          Do j_coll=1,i_coll
             N_COMP=N_COMP+1
             lambda=INDMAS(N_COMP,1)
             lambdb=INDMAS(N_COMP,2)
                                             Factor = one
             If(lambda.Eq.2.And.lambdb.Eq.2) Factor = Factor_q2**2
             If(lambda.Eq.3.And.lambdb.Eq.2) Factor = Factor_q2*Factor_q3
             If(lambda.Eq.3.And.lambdb.Eq.3) Factor = Factor_q3**2
             Write(6,'("  lambda=",i2," lambdb=",i2," M^1=",e14.4)') &
                              lambda,lambdb,EMOMET(i_coll,j_coll,1)*Factor
             Write(6,'(21x," M^2=",e14.4)') EMOMET(i_coll,j_coll,2)*Factor
             Write(6,'(21x," M^3=",e14.4)') EMOMET(i_coll,j_coll,3)*Factor
          End Do
       End Do
    End If
    !
  End Subroutine energy_moments
  !=======================================================================
  !>  This routine calculates matrix elements
  !>        \f$ \langle 0| \beta_i \beta_j \hat{F} |0 \rangle \f$
  !>  where |0> is the HFB vacuum (possibly under constraints), \f$ \beta_i \f$
  !>  and \f$ \beta_j \f$ are qp-annihilation operators and F is a one-body
  !>  hermitian operator:
  !>  \f[
  !>     \hat{F} = \sum_{kl} F_{kl} a^{\dagger}_{k}a_{l}
  !>  \f]
  !>
  !>  with \f$ a^{\dagger}_{k} \f$ and \f$ a_{l} \f$ single-particle
  !>  creation and anihilation operators, respectively.
  !=======================================================================
  Subroutine qp_basis_F20(nd,n_qp,Umatr,Vmatr,FOPQP,multMatElems)
    Use HFBTHO_utilities
    Implicit None
    Integer(ipr), Intent(in) :: nd !< size of the current block in the HO basis
    Integer(ipr), Intent(in) :: n_qp !< number of quasiparticles in the current block (n_qp \f$ \leq \f$ nd)
    Real(pr), Dimension(nd,n_qp), Intent(In) :: Umatr !< Matrix U of the Bogoliubov transformation
    Real(pr), Dimension(nd,n_qp), Intent(In) :: Vmatr !< Matrix V of the Bogoliubov transformation
    Real(pr), Dimension(n_qp,n_qp), Intent(inout) :: FOPQP !< Matrix of the operator F in the q.p. basis
    Real(pr), Allocatable, Intent(INOUT) :: multMatElems(:)
    !
    Integer(ipr) :: j,n1,n2
    Real(pr) :: hla
    Real(pr), Dimension(n_qp,nd) :: doubln
    Real(pr), Dimension(nd,nd) :: dblmul
    !
    ! matrix 'f' of the constraint in HO basis (size nd x nd)
    dblmul=zero
    j=0
    Do n1=1,nd
       Do n2=1,n1
          j=j+1;hla=multMatElems(j)
          dblmul(n1,n2)=hla;dblmul(n2,n1)=hla
       End Do
    End Do
    !
    ! temporary matrix
    doubln=zero
    !
    ! second term: v^{+} f^{*} u^{*} = v^{T} f u
    Call dgemm('t','n',n_qp,nd,nd,one,Vmatr,nd,dblmul,nd,zero,doubln,n_qp)
    Call dgemm('n','n',n_qp,n_qp,nd,one,doubln,n_qp,Umatr,nd,zero,FOPQP,n_qp)
    !
    ! first term:  u^{+} f v^{*} = u^{T} f v
    Call dgemm('t','n',n_qp,nd,nd,one,Umatr,nd,dblmul,nd,zero,doubln,n_qp)
    Call dgemm('n','n',n_qp,n_qp,nd,one,doubln,n_qp,Vmatr,nd,one,FOPQP,n_qp)
    !
  End Subroutine qp_basis_F20
  !=======================================================================
  !> This routine calculates matrix elements \f$ <0| \beta_i \beta_j F |0> \f$
  !> where \f$ |0> \f$ is the HFB vacuum (possibly under constraints), \f$ \beta_i \f$
  !> and \f$ \beta_j \f$ are qp-annihilation operators and F is a one-body
  !> hermitian operator:
  !> \f[
  !>        \hat{F} = \sum_{kl} F_{kl} a^{+}_{k}a_{l}
  !> \f]
  !> with \f$ a^{+}_{k}\f$ and \f$ a_{l} \f$ single-particle creation and anihilation
  !> operators, respectively.
  !>
  !=======================================================================
  Subroutine qp_basis_F11()
    Implicit None
    !
    !
  End Subroutine qp_basis_F11
  !=======================================================================
  !> This subroutine just prints the components of the collective
  !> inertia mass tensor and of the GCM metric.
  !=======================================================================
  Subroutine print_collective(iw)
    Implicit None
    !
    Integer(ipr), Intent(in) :: iw !< Fortran unit number where the data is printed
    Integer(ipr) :: i_coll,j_coll,N_COMP,lambda,lambdb,IUNITS
    Character(Len=1) :: FUNITS
    Character(Len=2) :: DIVBY2
    Character(Len=15) :: inertia_units
    Real(pr), Dimension(NDEFOR,NDEFOR) :: AUX1
    Real(pr), Dimension(NDCOMP) :: GMETRI
    !
    N_COMP=0
    Write(iw,'("  Inertia tensor   ATDHFB        GCM         Units")')
    Do i_coll=1,NOMULT
       Do j_coll=1,i_coll
          N_COMP=N_COMP+1
          lambda=INDMAS(N_COMP,1)
          lambdb=INDMAS(N_COMP,2)
          IUNITS=lambda+lambdb
          If(Mod(IUNITS,2).Eq.0) Then
             DIVBY2='  '
             IUNITS=IUNITS/2
          Else
             DIVBY2='/2'
          End If
          Write(FUNITS,'(i1)') IUNITS
          inertia_units='   '//'MeV-1.b^-'//FUNITS//DIVBY2
          Write(iw,'(6x,"M_",2i1,2x,2e14.4,a15)') lambda,lambdb,ATDMAS(N_COMP,3),GCMMAS(N_COMP,3),inertia_units
       End Do
    End Do
    !
    ! Copy the metric into the format used for printing
    AUX1 = COVMET(:,:,3)
    Call S2VCPY(AUX1,NOMULT,NDEFOR,GMETRI,NOCOMP,NDCOMP)
    N_COMP=0
    Write(iw,'("  GCM metric       Value        Units")')
    Do i_coll=1,NOMULT
       Do j_coll=1,i_coll
          N_COMP=N_COMP+1
          lambda=INDMAS(N_COMP,1)
          lambdb=INDMAS(N_COMP,2)
          IUNITS=lambda+lambdb
          If(Mod(IUNITS,2).Eq.0) Then
             DIVBY2='  '
             IUNITS=IUNITS/2
          Else
             DIVBY2='/2'
          End If
          Write(FUNITS,'(i1)') IUNITS
          inertia_units='     b^-'//FUNITS//DIVBY2//'    '
          Write(iw,'(6x,"G_",2i1,2x,1e14.4,a15)') lambda,lambdb,GMETRI(N_COMP),inertia_units
       End Do
    End Do
    Write(iw,*)
    !
  End Subroutine print_collective
  !=======================================================================
  !>  This routine copies a square symmetric matrix A into a vector V
  !=======================================================================
  Subroutine S2VCPY(A,N,NA,V,M,MV)
    Implicit None
    Integer(ipr), Intent(In) :: N  !< Actual dimensions
    Integer(ipr), Intent(In) :: NA !< Physical dimensions of A
    Integer(ipr), Intent(In) :: M  !< Number of elements in B
    Integer(ipr), Intent(In) :: MV !< Physical dimension of V
    Real(pr), Dimension(NA,NA), Intent(In) :: A !< Matrix to be copied into B
    Real(pr), Dimension(MV), Intent(Out) :: V !< Resulting, lower triangular part of A
    Integer(ipr) :: i,j,k
    If(N.Gt.NA.Or.M.Gt.MV) Stop &
        'S2VCPY: N .GT. NA; Physical dimensions of A or V exceeded'
    K=0
    Do I=1,N
       Do J=1,I
          K=K+1
          V(K)=A(I,J)
       End Do
    End Do
    Return
  End Subroutine S2VCPY
  !=======================================================================
  !>  This routine copies the vector v(:) into a square matrix a(:,:)
  !=======================================================================
  Subroutine V2SCPY(V,M,MV,A,N,NA)
    Implicit None
    Integer(ipr), Intent(In) :: N  !< Actual dimensions
    Integer(ipr), Intent(In) :: NA !< Physical dimensions of A
    Integer(ipr), Intent(In) :: M  !< Number of elements in B
    Integer(ipr), Intent(In) :: MV !< Physical dimension of V
    Real(pr), Dimension(MV), Intent(In) :: V !< Vector to be copied
    Real(pr), Dimension(NA,NA), Intent(Out) :: A !< Resulting symmetric matrix
    Integer(ipr) :: i,j,k
    If(N.Gt.NA.Or.M.Gt.MV) Stop &
        'V2SCPY: N .GT. NA; Physical dimensions of A or V exceeded'
    K=0
    Do I=1,N
       Do J=1,I
          K=K+1
          A(I,J)=V(K)
          A(J,I)=V(K)
       End Do
    End Do
    Return
  End Subroutine V2SCPY
  !=======================================================================
  !>  This routine does the following special matrix multiplication
  !>   \f$ A := A^{-1} B A^{-1} \f$
  !=======================================================================
  Subroutine INVMUL(A,B,XA,N,NX)
    Implicit None
    Integer(ipr), Intent(In) :: N  !< Actual dimension of A
    Integer(ipr), Intent(In) :: NX !< Physical dimensions of A
    Real(pr), Dimension(NX,NX), Intent(In) :: B !< Input matrix B
    Real(pr), Dimension(NX,NX), Intent(Inout) :: A  !< Upon entry, input matrix A; overwritten upon exit
    Real(pr), Dimension(NX,NX), Intent(Inout) :: XA !< Auxiliary array of the same dimension as A
    If(N.Gt.NX) STOP &
        'INVMUL: N .GT. NX; Physical dimensions of A or B exceeded'
    Call INVERT(A,XA,N,NX)
    A=MATMUL(A,MATMUL(B,A))
    Return
  End Subroutine INVMUL
  !=======================================================================
  !>  This function calculates the trace of a matrix
  !=======================================================================
  Real(pr) Function TRACEM(A,N,NX)
    Implicit None
    Integer(ipr), Intent(In) :: N  !< Actual dimension of A
    Integer(ipr), Intent(In) :: NX !< Physical dimensions of A
    Real(pr), Dimension(NX,NX), Intent(In) :: A !< Input matrix A
    Integer(ipr) :: i
    If(N.Gt.NX) Stop 'TRACEM: Dimension N > NX!'
    TRACEM=0.0_pr
    Do I=1,N
       TRACEM=TRACEM+A(I,I)
    End Do
    Return
  End Function TRACEM
  !=======================================================================
  !>  This subroutine calculates inverse of a matrix. A is the matrix
  !>  to be inverted and is destroyed in the routine. Upon exit, matrix
  !>  A contains the inverse \f$ A^{-1} \f$.
  !>
  !>  Method ....: CHOLESKY LU decomposition + LU back substitution
  !>  Routines ..: LUDCMP, LUBKSB  See: Num. Recipes., p 40)
  !=======================================================================
  Subroutine INVERT(A,Y,N,NP)
    Implicit None
    Integer(ipr), Parameter :: NMAX=500
    Integer(ipr), Intent(In) :: N  !< Actual dimension of A
    Integer(ipr), Intent(In) :: NP !< Physical dimensions of A
    Real(pr), Dimension(NP,NP), Intent(Inout) :: A !< Upon entry, input matrix A; inverse \f$ A^{-1} \f$ upon exit
    Real(pr), Dimension(NP,NP), Intent(Inout) :: Y !< Auxiliary array
    !
    Integer(ipr) :: i,j
    Integer(ipr), Dimension(NMAX) :: INDX
    Real(pr) :: D
    If(N.Gt.NP.Or.N.Gt.NMAX) &
       Stop 'INVERT: N TOO LARGE. INCREASE NP OR NMAX IN INVERT.'
    ! Set up identity matrix
    Y(:,:)=0.0_pr
    Do I=1,N
       Y(I,I)=1.0_pr
    End Do
    !
    ! Decompose the matrix just once
    Call LUDCMP(A,N,NP,INDX,D)
    !
    Do J=1,N
       ! Find inverse by columns
       Call LUBKSB(A,N,NP,INDX,Y(1,J))
       ! Note that fortran stores two-dimensional matrices by column,
       ! so Y(1,J) is the address of the J-th column of Y.
    End Do
    !
    ! Rewrite Y to A
    A=Y
    !
    Return
  End Subroutine INVERT
  !=======================================================================
  !
  !=======================================================================
  Subroutine ludcmp(a,n,np,indx,d)
    Implicit None
    Integer(ipr), Parameter :: NMAX=500
    Real(pr), Parameter :: TINY=1.0e-20
    Integer(ipr), Intent(In) :: n,np
    Integer(ipr), Dimension(n), Intent(Inout) :: indx
    Real(pr), Intent(Inout) :: d
    Real(pr), Dimension(np,np), Intent(Inout) :: a
    !
    Integer(ipr) :: i,imax,j,k
    Real(pr) :: aamax,sum,dum
    Real(pr), Dimension(NMAX) :: vv
    !
    d=1.0_pr
    imax=1
    Do i=1,n
       aamax=0.0_pr
       Do j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
       End Do ! j
       If (aamax.eq.0.) Write(*,*) 'singular matrix in ludcmp'
       vv(i)=1.0_pr/aamax
    End do
    Do j=1,n
       Do i=1,j-1
          sum=a(i,j)
          Do k=1,i-1
             sum=sum-a(i,k)*a(k,j)
          End Do ! k
          a(i,j)=sum
       End Do ! i
       aamax=0.0_pr
       Do i=j,n
          sum=a(i,j)
          Do k=1,j-1
             sum=sum-a(i,k)*a(k,j)
          End Do ! k
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) Then
             imax=i
             aamax=dum
          End If
       End Do ! i
       If (j.ne.imax) Then
          do k=1,n
             dum=a(imax,k)
             a(imax,k)=a(j,k)
             a(j,k)=dum
          End Do
          d=-d
          vv(imax)=vv(j)
       End If
       indx(j)=imax
       If(a(j,j).eq.0.) a(j,j)=TINY
       If(j.ne.n) Then
          dum=1.0_pr/a(j,j)
          do i=j+1,n
             a(i,j)=a(i,j)*dum
          End Do
       End If
    End do ! j
    !
    Return
  End Subroutine ludcmp
  !=======================================================================
  !
  !=======================================================================
  Subroutine lubksb(a,n,np,indx,b)
    Integer(ipr), Intent(In) :: n,np
    Integer(ipr), Dimension(n), Intent(Inout) :: indx
    Real(pr), Dimension(np,np), Intent(Inout) :: a
    Real(pr), Dimension(n), Intent(Inout) :: b
    !
    Integer(ipr) :: i,ii,j,ll
    Real(pr) :: sum
    !
    ii=0
    Do i=1,n
       ll=indx(i)
       sum=b(ll)
       b(ll)=b(i)
       If(ii.Ne.0) Then
          Do j=ii,i-1
             sum=sum-a(i,j)*b(j)
          End Do
       Else If(sum.Ne.0.) Then
          ii=i
       End If
       b(i)=sum
    End Do
    Do i=n,1,-1
       sum=b(i)
       Do j=i+1,n
          sum=sum-a(i,j)*b(j)
       End Do
       b(i)=sum/a(i,i)
    End Do
    !
    Return
  End Subroutine lubksb
  !=======================================================================
  !
  !=======================================================================
End Module HFBTHO_collective
