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
!                  FISSION FRAGMENT PROPERTIES PACKAGE                 !
!                                                                      !
! ==================================================================== !

!----------------------------------------------------------------------!
!>  This module contains a number of routines and functions dealing
!>  with the properties of the fission fragments. It contains the
!>  determination of the position of the neck between the two
!>  prefragments, the calculation of the matrix elements and expecta
!>  tion value of the Gaussian neck operator, as well as the expecta
!>  tion value of the multipole moments in the intrinsic frame of each
!>  prefragment. For l = m = 0, the latter can be used to define the
!>  number of protons and total number of particles in each prefragment.
!>
!>  <b>Inputs:</b>
!>    - ro(:,it) ..: density for isospin it on the integration mesh
!>
!>  <b>Outputs:</b>
!>    - Z_NECK, Q_NECK ..: position and size of the neck (in QNFIND,
!>                         dimensionless)
!>    - CENLEF, CENRIG ..: position of the center of mass for the left
!>                         and right fragments (in center_of_mass,
!>                         dimensionless)
!>    - QLMLEF, QLMRIG ..: multipole moments in the left and right
!>                         fragments
!>
!> @author
!> Nicolas Schunck
!----------------------------------------------------------------------
!  Subroutines: - neck_computeField(ib)
!               - QNFIND
!               - DEFMAS(NZMAXX,ZPOINT,SFACTO)
!               - center_of_mass(Z_POSI,CENLEF,CENRIG)
!               - QLMFRA(Z_POSI,LAMACT,QLMLEF,QLMRIG,CENLEF,CENRIG,I_TYPE)
!               - DEFSPH(LAMACT,COSTHE,SPHHAR)
!               - print_moments(iw)
!               - SIMP38(FUNCTI,NINTEG,X_STEP,RESULT)
!               - SPLINE(X,Y,N,YP1,YPN,Y2,IERROR)
!  Functions: - DERIVE(Z_POSI)
!             - QMOM_K(Z_POSI,KORDER)
!             - PROINT(ZPOINT)
!             - DEFLEG(LAMACT,MIUACT,ZVALUE)
!             - ZBRENT(FONCTI,XLOWER,XUPPER,TOLERA,IFOUND)
!             - SPLINT(XENTRY,YENTRY,Y2A,NPOINT,XARGUM)
!----------------------------------------------------------------------!
Module hfbtho_fission_fragments

  Use HFBTHO
  Use HFBTHO_utilities
  Use HFBTHO_gauss

  Implicit None

  Logical, PUBLIC, SAVE :: qp_occupation = .False.

  Integer(ipr), PRIVATE, SAVE :: NPOINT=500
  Integer(ipr), PRIVATE, SAVE :: debug_fission=0
  Real(pr), PRIVATE :: ZNMINI=0.0_pr
  Real(pr), PRIVATE :: ZNMAXI=0.0_pr

  ! Constraints on the neck
  Real(pr), PUBLIC, SAVE :: AN_VAL=1.0_pr !< Range in fermis in the Gaussian neck operator
  Real(pr), PUBLIC, SAVE :: Q_NECK=0.0_pr !< Value of the Gaussian neck operator
  Real(pr), PUBLIC, SAVE :: Z_NECK=0.0_pr !< Position of the neck (dimensionless)
  Real(pr), PUBLIC, SAVE :: CENLEF !< Position of the center of mass for the left fragment (dimensionless)
  Real(pr), PUBLIC, SAVE :: CENRIG !< Position of the center of mass for the right fragment (dimensionless)
  ! Fission fragments
  Real(pr), Allocatable, PUBLIC, SAVE :: QLMTOT(:,:) !< Arrays containing \f$ \langle \hat{Q}_{\ell} \rangle \f$ for
                                                     !< isoscalar pseudo-densities. First index is \f$ \ell \f$, second index is 0
                                                     !, for the left fragment, 1 for the right fragment
  Real(pr), Allocatable, PUBLIC, SAVE :: QLMPRO(:,:) !< Same as QLMTOT(:,:) but expectation values are computed with charge pseudo-densities

Contains
  !=======================================================================
  !> Calculates the matrix of the Gaussian neck operator in the current
  !> block. The result is stored in the global variable gaussian_neck(:)
  !=======================================================================
  Subroutine neck_computeField(ib)
    Implicit None
    Integer(ipr), Intent(in) :: ib !< Index of the current block
    Integer(ipr) :: i,ih,il,nd,nd2,ihli,ihil,im,n1,n2
    Integer(ipr) :: ja,jb,nsa,nsb,nsab,ssu,ssd
    Real(pr)    :: qhla,vh,fiun1,fiun2,fidn1,fidn2,vnhl
    Real(pr), Allocatable :: Vmom(:)
    Real(pr), Dimension(0:8) :: Qval
    Real(pr) :: z,rrr
    Integer(ipr) :: ndxmax
    Parameter(ndxmax=(n00max+2)*(n00max+2)/4)
    Real(pr) :: OMPFIU(ndxmax),OMPFID(ndxmax)
    !
    Allocate(Vmom(1:nghl))
    !
    Qval=zero
    !
    ! Compute value of neck operator on integration mesh
    Do ihli = 1,nghl
       z=fh(ihli)
       Vmom(ihli)=Exp(-((z-Z_NECK*bz)/AN_VAL)**2)
    End Do !ihli
    !
    ! Form matrix of the Gaussian neck operator in HO basis
    nd=id(ib); nd2=nd*nd; im=ia(ib)
    ! sum over gauss integration points
    Do ihil=1,nghl
       vnhl=Vmom(ihil)
       ! scan over basis states
       Do n1=1,nd
          ja=n1+im; nsa=NS(ja); ssu=Max(nsa,0); ssd=Max(-nsa,0)
          QHLA=QHLA_opt(ja,ihil)
          OMPFIU(N1)=QHLA*ssu
          OMPFID(N1)=QHLA*ssd
       End Do
       i=0
       Do n1=1,nd
          ja=n1+im; nsa=NS(ja)
          fiun1=OMPFIU(N1); fidn1=OMPFID(N1)
          Do n2=1,n1
             i=i+1; nsb=NS(n2+im); nsab=nsa+nsb; vh=0.0_pr
             If (nsab.Ne.0) Then
                 If (nsb.Gt.0) Then
                     fiun2 = OMPFIU(N2)
                     vh    = fiun1*fiun2
                 Else
                     fidn2 = OMPFID(N2)
                     vh    = fidn1*fidn2
                 End If
                 gaussian_neck(i)=gaussian_neck(i)+vh*vnhl
             End If
          End Do !n2
       End Do !n1
    End Do !ihil
       !
    Deallocate(Vmom)
    !
    Return
  End Subroutine neck_computeField
  !====================================================================
  !>  The subroutine finds the position of the neck, where
  !>  \f$ \langle \hat{Q}_{N}\rangle \f$ is the lowest.  At the minimum,
  !>  all partial derivatives of \f$ \langle \hat{Q}_{N}\rangle \f$ with
  !>  respect to \f$ x_{N}, y_{N}, z_{N} \f$ (coordinates of the neck)
  !>  should be 0. We find these 0 by using a variant of the Newton
  !>  method.
  !====================================================================
  Subroutine QNFIND()
    Logical :: symmetrical,interpolate,enforced_symmetry,any_odd_multipole
    Integer(ipr) :: i,l,IERROR,IFOUND,N_SCAN,KORDER,found,il,ih,ihli
    Real(pr) :: QN_VAL,QBEFOR,ZBEFOR,DQNPRE,DQNCUR,XLOWER,XUPPER
    Real(pr) :: dq,Z_INIT_back,Qtotal
    Real(pr) :: Z_MINI,Z_MAXI,Z_INIT,Q_INIT,Z_POSI,TOLERA
    !
    ! Initialization of the neck coordinates from previously found solutions.
    ! At the very first iteration, (X_NECK,Y_NECK,Z_NECK)=(0,0,0)
    Z_INIT=Z_NECK; Q_INIT=0.0_pr
    !
    symmetrical = .True.
    !
    ! Test if symmetry flags enforce symmetric nucleus
       enforced_symmetry = .True.
    If(.Not.Parity) Then
       enforced_symmetry = .False.
    End If
    ! Test if all multipole moments are actually zero
    any_odd_multipole = .False.
    Do l=1,lambdaMax,2
       Qtotal = qmoment(l,1) + qmoment(l,2)
       If(Abs(Qtotal).Gt.1.e-4) Then
          any_odd_multipole =.True.
       End If
    End Do
    ! System is symmetric either if symmetry flags enforce it, or if self-consistency makes it so
    If(.Not.enforced_symmetry) Then
       symmetrical = .Not.any_odd_multipole
    Else
       symmetrical = enforced_symmetry
    End If
    !
    KORDER=0
    !
    ! Initial scanning of the values of QN: we need this to get a somewhat
    ! reliable initial estimate of the position of the neck (only if
    ! reflection symmetry is broken). We search around the last known
    ! position of the neck.
    If(.Not.symmetrical) Then
       !
       DQNPRE=0.1_pr; N_SCAN=60
       Z_MINI=-5.0_pr; Z_MAXI=+5.0_pr
       dq = (Z_MAXI-Z_MINI)/Real(N_SCAN-1,Kind=pr)
       ZBEFOR=Z_MINI; QBEFOR=QMOM_K(ZBEFOR-dq,KORDER)
       !
       found = 0
       Do i=1,N_SCAN
          Z_POSI=Z_MINI+dq*Real(i-1,Kind=pr)
          QN_VAL=QMOM_K(Z_POSI,KORDER)
          DQNCUR=DERIVE(Z_POSI)
          If(debug_fission.Ge.2) Then
             Write(6,'("z=",f20.14," QN=",f20.14," dQN(analytic)=",f20.14," dQN(numerical)=",f20.14)') &
                        Z_POSI,QN_VAL,DQNCUR,(QN_VAL-QBEFOR)/dq
          End If
          ! Found the minimum
          If(DQNCUR.Gt.0.0_pr.And.DQNPRE.Lt.0.0_pr) Then
             Z_INIT=0.5_pr*(Z_POSI+ZBEFOR)
             ZNMINI=Z_INIT-dq; ZNMAXI=Z_INIT+dq
             Q_INIT=QN_VAL
             found =found + 1
          End If
          ZBEFOR=Z_POSI
          QBEFOR=QN_VAL
          DQNPRE=DQNCUR
       End Do
       !
       If(debug_fission.Ge.1) Then
          Write(6,'("NUMITE=",i4," Z_NECK=",f20.14)') iiter,Z_NECK
          Write(6,'("ZNMINI=",f20.14," ZNMAXI=",f20.14)') ZNMINI,ZNMAXI
          Write(6,'("found=",i4," Z_POSI=",f20.14," QN=",f20.14)') found,Z_INIT,Q_INIT
       End If
       !
       ! Search for the z-coordinate of the neck (most likely to be different from 0)
       If(found.Gt.0) Then
          Z_INIT_back = Z_INIT ! backup of rough estimate
          XLOWER=ZNMINI; XUPPER=ZNMAXI; TOLERA=1.D-5; IFOUND=0
          Z_INIT=ZBRENT(DERIVE,XLOWER,XUPPER,TOLERA,IFOUND)
          If(debug_fission.Ge.1) Then
             Write(6,'("found (Newton)=",i4," Z_INIT=",f20.14)') ifound,Z_INIT
          End If
          ! If ZBRENT gives a value too different from the rough estimate, there is a problem...
          If(Abs(Z_INIT-Z_INIT_back).Gt.0.5_pr*Abs(dq)) Then
             If(debug_fission.Ge.1) Then
                Write(6,'("Initial guess: ",f20.14," ZBRENT found: ", F20.14)') Z_INIT_back, Z_INIT
             End If
             Z_INIT = Z_INIT_back
          End If
          If(IFOUND.EQ.0) Z_INIT = Z_INIT_back
       Else
          ! If no candidate for the neck position was found, we default to 0
          Z_INIT=0.0_pr
       End If
       !
       If(iiter.Ge.2.And.neck_constraints) Then
          Z_NECK=0.25*Z_NECK + (1.0-0.25)*Z_INIT
       Else
          Z_NECK=Z_INIT
       End If
       Q_NECK=QMOM_K(Z_NECK,KORDER)
       If(debug_fission.Ge.1) Then
          Write(6,'("Final value zN: ",f20.14," qN: ", 2F20.14)') Z_NECK, Q_NECK
       End If
       !
    Else
       ! Reflection symmetry is conserved, the neck is at z=0, we compute QN
       Z_NECK=0.0_pr
       Q_NECK=QMOM_K(Z_NECK,KORDER)
    End If
    !
    neckValue = Q_NECK
    !
  End Subroutine QNFIND
  !====================================================================
  !>  This function computes the first derivative of the Gaussian
  !>  neck expectation value with respect to the z-coordinate of the
  !>  neck.
  !====================================================================
  Real(pr) Function DERIVE(Z_POSI)
    Real(pr), INTENT(IN) :: Z_POSI !< Current position of the neck (dimensionless)
    !
    Integer(ipr) :: KORDER
    Real(pr) :: DERIV1,PREFAC
    !
    DERIV1=0.0_pr
    PREFAC=2.0_pr*bz*bz/AN_VAL**2
    KORDER=1; DERIV1=QMOM_K(Z_POSI,KORDER)
    KORDER=0; DERIV1=DERIV1-Z_POSI*QMOM_K(Z_POSI,KORDER)
    DERIVE=PREFAC*DERIV1

  End Function DERIVE
  !====================================================================
  !>  This function computes moments \f$ I_k \f$ of order k of the
  !>  Gaussian neck operator
  !====================================================================
  Real(pr) Function QMOM_K(Z_POSI,KORDER)
    Integer(ipr), INTENT(IN) :: KORDER
    Real(pr), INTENT(IN) :: Z_POSI
    !
    Integer(ipr) :: il,ih,ihli
    Real(pr) :: QN_VAL,ZNECKZ,COORDO,PREFAC
    !
    QN_VAL=0.0_pr
    Do il=1,ngl
       Do ih=1,ngh
          ihli = ih + (il-1)*ngh
          ZNECKZ=Exp(-((xh(ih)-Z_POSI)*bz/AN_VAL)**2)
          COORDO=xh(ih)
          If(KORDER.Eq.0) Then
             PREFAC=1.0_pr
          Else
             PREFAC=COORDO**KORDER
          End If
          QN_VAL=QN_VAL+wdcor(ihli)*(ro(ihli,1)+ro(ihli,2))*ZNECKZ*PREFAC
       End Do
    End Do
    !
    QMOM_K=QN_VAL
    !
  End Function QMOM_K
  !====================================================================
  !>  This function computes moments \f$ I_k \f$ of order k of the
  !>  Gaussian neck operator
  !====================================================================
  Real(pr) Function test_density()
    Integer(ipr) :: il,ih,ihli
    Real(pr) :: QN_VAL
    !
    QN_VAL=0.0_pr
    Do il=1,ngl
       Do ih=1,ngh
          ihli = ih + (il-1)*ngh
          QN_VAL=QN_VAL+wdcor(ihli)*(ro(ihli,1)+ro(ihli,2))
       End Do
    End Do
    !
    test_density=QN_VAL
    !
  End Function test_density
  !=======================================================================
  !
  !=======================================================================
  Subroutine neck_expectation(it,ib,qval,rho,dd)
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Integer(ipr), Intent(In) :: it,ib
    Real(pr), Allocatable, Intent(In) :: rho(:,:)
    Real(pr), Intent(Inout) :: qval,dd
    !
    Integer(ipr) :: nd,nd2,nhfb,j,n1,n2,i
    Real(pr) :: hla
    Real(pr), Allocatable :: dblmul(:,:),qblock(:,:)
    !
    nd=id(ib); nd2=nd*nd; nhfb=nd+nd;
    Do n1=1,nd
       dd=dd+rho(n1,n1)
    End Do

    Allocate(dblmul(nd,nd));dblmul=zero
    gaussian_neck=zero
    Call neck_computeField(ib)
    j=0
    Do n1=1,nd
       Do n2=1,n1
          j=j+1;hla=gaussian_neck(j)
          dblmul(n1,n2)=hla;dblmul(n2,n1)=hla
       End Do
    End Do
    Allocate(qblock(nd,nd))
    Call dgemm('n','n',nd,nd,nd,one,rho,nd,dblmul,nd,zero,qblock,nd)
    Do n1=1,nd
       qval=qval+qblock(n1,n1)
    End Do

    Return
  End Subroutine neck_expectation
  !====================================================================
  !
  !=======================================================================
  Subroutine neck()
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Integer(ipr) :: it,ib,nd,nd2,n_qp,k,n1,n2,kk
    Real(pr) :: qvala,dd2
    Real(pr), Allocatable :: rho(:,:), Umatr(:,:),Vmatr(:,:)
    !
    qvala=zero; dd2=zero
    !
    it = 1
    !
    Do ib = 1,nb
       ! Get characteristics of the current block: size of HO basis and number of qp
       nd=id(ib); nd2=nd*nd; n_qp=kd(ib,it)
       If(n_qp.Gt.0) Then
          Allocate(gaussian_neck(1:nd2))
          Allocate(Umatr(nd,n_qp)); Umatr=zero
          Allocate(Vmatr(nd,n_qp)); Vmatr=zero
          Do k=1,n_qp
             Do n2=1,nd
                Vmatr(n2,k)=RVqpN(KpwiN(ka(ib,it)+k)+n2)
                Umatr(n2,k)=RUqpN(KpwiN(ka(ib,it)+k)+n2)
             End Do
          End Do
          ! In configuration space, test that Tr(rho) = N (contained in dd1 and dd2), and that
          ! Tr(Q*rho) = Q (in qvala and qvalb)
          Allocate(rho(nd,nd)); rho=zero
          Call dgemm('n','t',nd,nd,n_qp,one,Vmatr,nd,Vmatr,nd,zero,rho,nd)
          Call neck_expectation(it,ib,qvala,rho,dd2)
          Deallocate(rho)
          !
          Deallocate(Umatr,Vmatr)
          Deallocate(gaussian_neck)
       End If
    End Do
    !
    it = 2
    !
    Do ib = 1,nb
       ! Get characteristics of the current block: size of HO basis and number of qp
       nd=id(ib); nd2=nd*nd; n_qp=kd(ib,it)
       If(n_qp.Gt.0) Then
          Allocate(gaussian_neck(1:nd2))
          Allocate(Umatr(nd,n_qp)); Umatr=zero
          Allocate(Vmatr(nd,n_qp)); Vmatr=zero
          Do k=1,n_qp
             Do n2=1,nd
                Vmatr(n2,k)=RVqpP(KpwiP(ka(ib,it)+k)+n2)
                Umatr(n2,k)=RUqpP(KpwiP(ka(ib,it)+k)+n2)
             End Do
          End Do
          ! In configuration space, test that Tr(rho) = N (contained in dd1 and dd2), and that
          ! Tr(Q*rho) = Q (in qvala and qvalb)
          Allocate(rho(nd,nd)); rho=zero
          Call dgemm('n','t',nd,nd,n_qp,one,Vmatr,nd,Vmatr,nd,zero,rho,nd)
          Call neck_expectation(it,ib,qvala,rho,dd2)
          Deallocate(rho)
          !
          Deallocate(Umatr,Vmatr)
          Deallocate(gaussian_neck)
       End If
    End Do
    !
    Write(6,'("qvala=",f14.7," dd2=",f14.7)') 2*qvala,2*dd2
    !
    Return
  End Subroutine neck
  !====================================================================
  !
  !=======================================================================
  Subroutine wave_localization(ZPOINT)
    Use HFBTHO_utilities
    Use HFBTHO
    Use HFBTHO_gauss
    Implicit None
    Real(pr), INTENT(IN) :: ZPOINT !< Position of the neck \f$ \xi \f$ (dimensionless)
    Integer(ipr) :: it,ib,nd,n_qp,k,n1,n2,kk,im
    Integer(ipr) :: nr2,nl2,ns2,nz2,nr1,nl1,ns1,nz1
    Real(pr) :: cN
    Real(pr), Allocatable :: SFACTO(:,:) !< Array containing \f$ S_{m_{z}n_{z}} \f$
    !
    Call DEVHER(nzx)

    If(.Not.Allocated(SFACTO)) Allocate(SFACTO(0:nzx,0:nzx))
    Call DEFMAS(nzx,ZPOINT,SFACTO)
    !
    it = 1
    !
    Do ib = 1,nb
       ! Get characteristics of the current block: size of HO basis, number of qp and running index im
       nd=id(ib); n_qp=kd(ib,it); im=ia(ib)
       If(n_qp.Gt.0) Then
          Do k=1,n_qp
             cN = zero
             Do n2=1,nd
                nr2=nr(im+n2); nl2=nl(im+n2); ns2=ns(im+n2); nz2=nz(im+n2)
                Do n1=1,nd
                   nr1=nr(im+n1); nl1=nl(im+n1); ns1=ns(im+n1); nz1=nz(im+n1)
                   If(nr1.Eq.nr2.And.nl1.Eq.nl2.And.ns1.Eq.ns2) Then
                      cN = cN + RVqpN(KpwiN(ka(ib,it)+k)+n1)*RVqpN(KpwiN(ka(ib,it)+k)+n2)*SFACTO(nz1,nz2)
                   End If
                End Do
             End Do
             If(uk(ka(ib,it)+k,it).Gt.1.0e-3_pr) Then
                Write(6,'("ib=",i3," k=",i5," Eqp=",f20.14," v2l=",e24.12," v2=",f17.14," diff=",e24.12)') &
                           ib,k,REqpN(KqpN(ka(ib,it)+k)),cN,uk(ka(ib,it)+k,it),uk(ka(ib,it)+k,it)-cN
             End If
          End Do
       End If
    End Do
    !
    Return
  End Subroutine wave_localization
  !====================================================================
  !>  This subroutine computes:
  !>   \f[
  !>     S_{m_{z}n_{z}} = \sum_{k}^{n_z+m_z} C_{n_{z}m_{z}}^{k}(00)
  !>              \int_{z_{N}}^{+\infty} d\xi H_{k}(\xi) e{-\xi^2}
  !>   \f]
  !>  This is needed to evaluate the  expectation value of the
  !>  density in the fission fragments. Inputs are NZMAXX, the
  !>  maximum value for either nz or mz, and ZPOINT, the lower
  !>  bound of the generalized integral
  !====================================================================
  Subroutine DEFMAS(NZMAXX,ZPOINT,SFACTO)
    Integer(ipr), INTENT(IN) :: NZMAXX !< Maximum value for either \f$ n_{z} \f$ or \f$ m_{z} \f$
    Real(pr), INTENT(IN) :: ZPOINT !< Position of the neck \f$ \xi \f$ (dimensionless)
    Real(pr), Allocatable, INTENT(INOUT) :: SFACTO(:,:) !< Array containing \f$ S_{m_{z}n_{z}} \f$

    Integer(ipr) :: NORDER,NZ,MZ,K,I
    Real(pr) :: PIARGU,VALINT,XPOINT,SUMVAL,F_INTE
    Real(pr), Allocatable :: PHERMI(:),DHERMI(:),ZHERMI(:)
    !
    PIARGU=four*Atan(one)
    !
    ! Computing the value of the probability distribution function
    VALINT=PROINT(ZPOINT)
    !
    ! Compute value of Hermite polynomials up to order Max(Nz+Mz) at z=0
    XPOINT=zero; NORDER=2*NZMAXX
    Allocate(PHERMI(1:NORDER+1)); Allocate(DHERMI(1:NORDER+1))
    Call D_HERM(XPOINT,NORDER,PHERMI,DHERMI,NORDER+1)
    !
    ! Compute value of Hermite polynomials up to order Max(Nz+Mz) at z=Abs(ZPOINT)
    XPOINT=Abs(ZPOINT); NORDER=2*NZMAXX
    Allocate(ZHERMI(1:NORDER+1))
    Call D_HERM(XPOINT,NORDER,ZHERMI,DHERMI,NORDER+1)
    !
    ! Loop over bra and ket z quantum number
    Do NZ=0,NZMAXX
       Do MZ=0,NZMAXX
          ! Initial value of the sum for K = 0
          SUMVAL=COEF00(0,MZ,NZ) * half*Sqrt(Sqrt(PIARGU))*(one - VALINT)
          If(ZPOINT.Lt.0.0_pr) SUMVAL=COEF00(0,MZ,NZ)*half*Sqrt(Sqrt(PIARGU))*(one + VALINT)
          Do K=1,NZ+MZ
             ! Hermite polynomial of order k is stored in PHERMI(k+1), but this shift does
             ! *not* apply to its normalization coefficient, which is stored in HERFAC(k)
             ! as expected. Below, we need Hermite polynomials of order k-1 and normalization
             ! coefficients of order k
             If(Mod(K,2).Eq.1) Then
                F_INTE = ZHERMI(K)*Exp(-ZPOINT*ZPOINT)/HERFAC(K)
             Else
                F_INTE =(ZHERMI(K)*Exp(-ZPOINT*ZPOINT) - PHERMI(K)) / HERFAC(K)
                If(ZPOINT.Lt.0.0_pr) F_INTE = - F_INTE
             End If
             SUMVAL=SUMVAL+COEF00(K,MZ,NZ) * F_INTE
          End Do
          SFACTO(MZ,NZ)=SUMVAL
       End Do
    End Do
    !
  End Subroutine DEFMAS
  !====================================================================
  !> This function computes position of the center of mass of the left
  !> and right fragments. The resulting quantity is in fermis.
  !====================================================================
  Subroutine center_of_mass(Z_POSI,CENLEF,CENRIG)
    Real(pr), INTENT(IN) :: Z_POSI !< Position of the neck (dimensionless)
    Real(pr), INTENT(INOUT) :: CENLEF !< Center of mass of the left fragment (dimensionless)
    Real(pr), INTENT(INOUT) :: CENRIG !< Center of mass of the right fragment (dimensionless)

    Integer(ipr) :: IX,IY,IZ,IERROR,i,ih,il,ihli
    Real(pr) :: DERIV1,DERIVN,XBEGIN,XFINIS,XARGUM,X_STEP,RHOTMP,ZPOTMP,dh
    Real(pr) :: Z_ACTU,ZPOINT,DENLOC,W_HERM,RHONOR,Z_NORM
    Real(pr), Allocatable :: FUNCTI(:),GUNCTI(:),XVALUE(:),ZINTER(:),AUXSTO(:)
    Real(pr), Allocatable :: DINTER(:,:)

    Z_ACTU=Z_POSI

    ! Integrations over the transverse coordinates are done using the dimensionless units
    ! xi_x and xi_y. In the longitudinal direction, the density is defined as a function
    ! of the dimensionless variable xi_z = bz * z. The operator O(z) = z must thus be
    ! defined on the grid z = xi_z/bz with integration limits given below.
    Allocate(FUNCTI(1:ngh),XVALUE(1:ngh))
    Allocate(AUXSTO(1:ngh),ZINTER(1:NPOINT))
    Allocate(DINTER(1:ngl,1:NPOINT))
    XBEGIN=xh(1)*bz + 1.D-14
    XFINIS=Z_ACTU*bz
    ! Boundary test: make sure there is a valid integration interval
    If(Z_ACTU.LT.XBEGIN) Then
       Write(6,'("ATTENTION! In center_of_mass(): changing neck position from ",f20.14," to ",f20.14)') &
                  Z_ACTU,XBEGIN + 1.D-12
       Z_ACTU = XBEGIN + 1.D-12
    End If
    !
    Do il=1,ngl
       Do ih=1,ngh
          ihli = ih + (il-1)*ngh
          XVALUE(ih)=xh(ih)*bz
          FUNCTI(ih)=ro(ihli,1)+ro(ihli,2)
       End Do
       ! Approximation of derivatives at the boundaries
       dh=XVALUE(2)-XVALUE(1) ! >0
       DERIV1=+0.5_pr*(-3.0_pr*FUNCTI(1)  +4.0_pr*FUNCTI(2)    -FUNCTI(3))    /dh
       DERIVN=-0.5_pr*(-3.0_pr*FUNCTI(ngh)+4.0_pr*FUNCTI(ngh-1)-FUNCTI(ngh-2))/dh
       ! Define spline coefficients
       IERROR=0; Call SPLINE(XVALUE,FUNCTI,ngh,DERIV1,DERIVN,AUXSTO,IERROR)
       ! Define interpolated mesh and value of the function of that mesh
       If(IERROR.Eq.0) Then
          Do i=1,NPOINT
             XARGUM=XBEGIN+(XFINIS-XBEGIN)*Real(i-1,Kind=pr)/Real(NPOINT-1,Kind=pr)
             ZINTER(i)=XARGUM
             DINTER(il,i)=SPLINT(XVALUE,FUNCTI,AUXSTO,ngh,XARGUM)
          End Do
       Else
          Write(6,'("In center_of_mass() - Error in interpolating the function (right fragment)!")')
       End If
    End Do
    X_STEP=(ZINTER(NPOINT)-ZINTER(1))/Real(NPOINT-1,Kind=pr)
    Deallocate(FUNCTI)
    !
    ! Computing the integrals \int \rho(z) and \int z\rho(z)
    !          -\infty < x < +\infty
    !          -\infty < y < +\infty
    !          -\infty < z < +z_{N}
    ! Integrations over x and y are performed "exactly" by Gauss-Hermite
    ! quadratures, integration over z is performed numerically using the
    ! Simpson 3/8 rule
    Allocate(FUNCTI(1:NPOINT),GUNCTI(1:NPOINT))
    Do i=1,NPOINT
       FUNCTI(i)=0.0_pr; RHOTMP=0.0_pr
       GUNCTI(i)=0.0_pr; ZPOTMP=0.0_pr
       Do il=1,ngl
          ZPOINT=ZINTER(i)
          DENLOC=DINTER(il,i)
          W_HERM=pi*wl(il)*bp*bp
          RHOTMP=RHOTMP+W_HERM*DENLOC
          ZPOTMP=ZPOTMP+W_HERM*DENLOC*ZPOINT
       End Do
       FUNCTI(i)=RHOTMP; GUNCTI(i)=ZPOTMP
    End Do
    ! Integrating over z, from -infty to zN
    RHONOR=0.0_pr; Call SIMP38(FUNCTI,NPOINT,X_STEP,RHONOR)
    Z_NORM=0.0_pr; Call SIMP38(GUNCTI,NPOINT,X_STEP,Z_NORM)
    CENLEF=Z_NORM/RHONOR
    Deallocate(FUNCTI,GUNCTI,XVALUE,DINTER,AUXSTO,ZINTER)
    !
    ! Repeating this whole procedure for the right fragment
    Allocate(FUNCTI(1:ngh),XVALUE(1:ngh))
    Allocate(AUXSTO(1:ngh),ZINTER(1:NPOINT))
    Allocate(DINTER(1:ngl,1:NPOINT))
    XBEGIN=Z_ACTU*bz
    XFINIS=xh(ngh)*bz - 1.D-14
    If(Z_ACTU.Gt.XFINIS) Then
       Write(6,'("ATTENTION! In center_of_mass(): changing neck position from ",f20.14," to ",f20.14)') &
                  Z_ACTU,XFINIS - 1.D-12
       Z_ACTU = XFINIS - 1.D-12
    End If
    Do il=1,ngl
       Do ih=1,ngh
          ihli = ih + (il-1)*ngh
          XVALUE(ih)=xh(ih)*bz
          FUNCTI(ih)=ro(ihli,1)+ro(ihli,2)
       End Do
       ! Approximation of derivatives at the boundaries
       dh=XVALUE(2)-XVALUE(1) ! >0
       DERIV1=+0.5_pr*(-3.0_pr*FUNCTI(1)  +4.0_pr*FUNCTI(2)    -FUNCTI(3))       /dh
       DERIVN=-0.5_pr*(-3.0_pr*FUNCTI(ngh)+4.0_pr*FUNCTI(ngh-1)-FUNCTI(ngh-2))/dh
       ! Define spline coefficients
       IERROR=0; Call SPLINE(XVALUE,FUNCTI,ngh,DERIV1,DERIVN,AUXSTO,IERROR)
       ! Define interpolated mesh and value of the function of that mesh
       If(IERROR.Eq.0) Then
          Do i=1,NPOINT
             XARGUM=XBEGIN+(XFINIS-XBEGIN)*Real(i-1,Kind=pr)/Real(NPOINT-1,Kind=pr)
             ZINTER(i)=XARGUM
             DINTER(il,i)=SPLINT(XVALUE,FUNCTI,AUXSTO,ngh,XARGUM)
          End Do
       Else
          Write(6,'("In center_of_mass() - Error in interpolating the function (left fragment)!")')
       End If
    End Do
    X_STEP=(ZINTER(NPOINT)-ZINTER(1))/Real(NPOINT-1,Kind=pr)
    Deallocate(FUNCTI)
    !
    Allocate(FUNCTI(1:NPOINT),GUNCTI(1:NPOINT))
    Do i=1,NPOINT
       FUNCTI(i)=0.0_pr; RHOTMP=0.0_pr
       GUNCTI(i)=0.0_pr; ZPOTMP=0.0_pr
       Do il=1,ngl
          ZPOINT=ZINTER(i)
          DENLOC=DINTER(il,i)
          W_HERM=pi*wl(il)*bp*bp
          RHOTMP=RHOTMP+W_HERM*DENLOC
          ZPOTMP=ZPOTMP+W_HERM*DENLOC*ZPOINT
       End Do
       FUNCTI(i)=RHOTMP; GUNCTI(i)=ZPOTMP
    End Do
    RHONOR=0.0_pr; Call SIMP38(FUNCTI,NPOINT,X_STEP,RHONOR)
    Z_NORM=0.0_pr; Call SIMP38(GUNCTI,NPOINT,X_STEP,Z_NORM)
    CENRIG=Z_NORM/RHONOR
    !
    Deallocate(FUNCTI,GUNCTI,AUXSTO,ZINTER,DINTER,XVALUE)
    !
  End Subroutine center_of_mass
  !====================================================================
  !> Function computes the expectation value of the multipole moment
  !> operators \f$ Q_{\ell m} \f$ in the left and right fragments,
  !> assuming the position of the neck is \f$ z_{N} \f$
  !====================================================================
  Subroutine QLMFRA(Z_POSI,LAMACT,QLMLEF,QLMRIG,CENLEF,CENRIG,I_TYPE)
    Integer(ipr), INTENT(IN) :: LAMACT !< \f$ \ell \f$ value
    Integer(ipr), INTENT(IN) :: I_TYPE !< m value
    Real(pr), INTENT(IN) :: Z_POSI !< Position of the neck (dimensionless)
    Real(pr), INTENT(IN) :: CENLEF !< Center of mass of the left fragment (dimensionless)
    Real(pr), INTENT(IN) :: CENRIG !< Center of mass of the right fragment (dimensionless)
    Real(pr), INTENT(INOUT) :: QLMLEF !< \f$ \langle \hat{Q}_{\ell m} \rangle \f$ for the left fragment (dimensionless)
    Real(pr), INTENT(INOUT) :: QLMRIG !< \f$ \langle \hat{Q}_{\ell m} \rangle \f$ for the right fragment (dimensionless)

    Integer(ipr) :: IX,IY,IZ,IERROR,i,ih,il,ihli
    Real(pr) :: DERIV1,DERIVN,XBEGIN,XFINIS,XARGUM,X_STEP,pi,dh
    Real(pr) :: Z_ACTU,ZPOINT,DENLOC,W_HERM,QLMSUM,QLMVAL,RESULT
    Real(pr) :: RADIUS,COSTHE,epsilon,one
    Real(pr), Dimension(1:2) :: SPHHAR
    Real(pr), Allocatable :: FUNCTI(:),XVALUE(:),ZINTER(:),AUXSTO(:)
    Real(pr), Allocatable :: DINTER(:,:)

    pi = 4.0_pr*Atan(1.0_pr)
    Z_ACTU=Z_POSI

    ! Define boundaries for the integration over z, see comments in subroutine center_of_mass()
    Allocate(FUNCTI(1:ngh),XVALUE(1:ngh))
    Allocate(AUXSTO(1:ngh),ZINTER(1:NPOINT))
    Allocate(DINTER(1:ngl,1:NPOINT))
    XBEGIN=xh(1)*bz + 1.D-14
    XFINIS=Z_ACTU*bz
    ! Boundary test: make sure there is a valid integration interval
    If(Z_ACTU.LT.XBEGIN) Then
       Write(6,'("ATTENTION! In QLMFRA(): changing neck position from ",f20.14," to ",f20.14)') &
                  Z_ACTU,XBEGIN + 1.D-12
       Z_ACTU = XBEGIN + 1.D-12
    End If
    !
    Do il=1,ngl
       ! See comments in subroutine center_of_mass()
       ! Total density
       If(I_TYPE.Eq.1) Then
          Do ih=1,ngh
             ihli = ih + (il-1)*ngh
             XVALUE(ih)=xh(ih)*bz
             FUNCTI(ih)=ro(ihli,1)+ro(ihli,2)
          End Do
       End If
       ! Proton (=charge in HFBTHO) density
       If(I_TYPE.Eq.2) Then
          Do ih=1,ngh
             ihli = ih + (il-1)*ngh
             XVALUE(ih)=xh(ih)*bz
             FUNCTI(ih)=ro(ihli,2)
          End Do
       End If
       ! Approximation of derivatives at the boundaries
       dh=XVALUE(2)-XVALUE(1) ! >0
       DERIV1=+0.5_pr*(-3.0_pr*FUNCTI(1)  +4.0_pr*FUNCTI(2)    -FUNCTI(3))    /dh
       DERIVN=-0.5_pr*(-3.0_pr*FUNCTI(ngh)+4.0_pr*FUNCTI(ngh-1)-FUNCTI(ngh-2))/dh
       ! Define spline coefficients
       IERROR=0; Call SPLINE(XVALUE,FUNCTI,ngh,DERIV1,DERIVN,AUXSTO,IERROR)
       ! Define interpolated mesh and value of the function of that mesh
       If(IERROR.Eq.0) Then
          Do i=1,NPOINT
             XARGUM=XBEGIN+(XFINIS-XBEGIN)*Real(i-1,Kind=pr)/Real(NPOINT-1,Kind=pr)
             ZINTER(i)=XARGUM
             DINTER(il,i)=SPLINT(XVALUE,FUNCTI,AUXSTO,ngh,XARGUM)
          End Do
       Else
          Write(6,'("In QLMFRA() - Error in interpolating the function (right fragment)!")')
       End If
    End Do
    X_STEP=(ZINTER(NPOINT)-ZINTER(1))/Real(NPOINT-1,Kind=pr)
    Deallocate(FUNCTI)
    !
    ! Computing the expectation values of multipole moment for the left fragment
    Allocate(FUNCTI(1:NPOINT))
    epsilon=1.D-14; one=1.0_pr
    Do i=1,NPOINT
       FUNCTI(i)=0.0_pr; QLMSUM=0.0_pr
       Do il=1,ngl
          ! Coordinates (x,y,z) for the multipole moments must be in fermis
          ZPOINT=ZINTER(i)-CENLEF ! shift with respect to the c.o.m. of left fragment
          RADIUS=Sqrt(xl(il)*bp**2+ZPOINT**2)
          ! Angle theta
          If(RADIUS.Le.epsilon) Then
             COSTHE=0.0_pr
          Else
             COSTHE=ZPOINT/RADIUS
          End If
          If(Abs(COSTHE-one).Le.epsilon) Then
             COSTHE=one-2.0_pr*epsilon
          End If
          If(Abs(COSTHE+one).Le.epsilon) Then
             COSTHE=2.0_pr*epsilon-one
          End If
          Call DEFSPH(LAMACT,COSTHE,SPHHAR)
          If(LAMACT.Eq.0) Then
             QLMVAL=1.0_pr
          Else
             If(RADIUS.Ge.epsilon) Then
                QLMVAL=RADIUS**(LAMACT)*SPHHAR(1)
             Else
                QLMVAL=0.0_pr
             End If
          End If
          DENLOC=DINTER(il,i)
          W_HERM=pi*wl(il)*bp*bp
          QLMSUM=QLMSUM+W_HERM*DENLOC*QLMVAL
       End Do
       FUNCTI(i)=QLMSUM
    End Do
    ! Integrating over z, from -infty to zN
    Call SIMP38(FUNCTI,NPOINT,X_STEP,RESULT)
    QLMLEF=RESULT
    Deallocate(FUNCTI,XVALUE,DINTER,AUXSTO,ZINTER)
    !
    ! Repeating the procedure for the right fragment
    Allocate(FUNCTI(1:ngh),XVALUE(1:ngh))
    Allocate(AUXSTO(1:ngh),ZINTER(1:NPOINT))
    Allocate(DINTER(1:ngl,1:NPOINT))
    XBEGIN=Z_ACTU*bz
    XFINIS=xh(ngh)*bz - 1.D-14
    ! Boundary test: make sure there is a valid integration interval
    If(Z_ACTU.Gt.XFINIS) Then
       Write(6,'("ATTENTION! In QLMFRA(): changing neck position from ",f20.14," to ",f20.14)') &
                  Z_ACTU,XFINIS - 1.D-12
       Z_ACTU = XFINIS - 1.D-12
    End If
    Do il=1,ngl
       If(I_TYPE.Eq.1) Then
          Do ih=1,ngh
             ihli = ih + (il-1)*ngh
             XVALUE(ih)=xh(ih)*bz
             FUNCTI(ih)=ro(ihli,1)+ro(ihli,2)
          End Do
       End If
       ! Proton (=charge in HFBTHO) density
       If(I_TYPE.Eq.2) Then
          Do ih=1,ngh
             ihli = ih + (il-1)*ngh
             XVALUE(ih)=xh(ih)*bz
             FUNCTI(ih)=ro(ihli,2)
          End Do
       End If
       ! Approximation of derivatives at the boundaries
       dh=XVALUE(2)-XVALUE(1) ! >0
       DERIV1=+0.5_pr*(-3.0_pr*FUNCTI(1)  +4.0_pr*FUNCTI(2)    -FUNCTI(3))    /dh
       DERIVN=-0.5_pr*(-3.0_pr*FUNCTI(ngh)+4.0_pr*FUNCTI(ngh-1)-FUNCTI(ngh-2))/dh
       ! Define spline coefficients
       IERROR=0; Call SPLINE(XVALUE,FUNCTI,ngh,DERIV1,DERIVN,AUXSTO,IERROR)
       ! Define interpolated mesh and value of the function of that mesh
       If(IERROR.Eq.0) Then
          Do i=1,NPOINT
             XARGUM=XBEGIN+(XFINIS-XBEGIN)*Real(i-1,Kind=pr)/Real(NPOINT-1,Kind=pr)
             ZINTER(i)=XARGUM
             DINTER(il,i)=SPLINT(XVALUE,FUNCTI,AUXSTO,ngh,XARGUM)
          End Do
       Else
          Write(6,'("In QLMFRA() - Error in interpolating the function (left fragment)!")')
       End If
    End Do
    X_STEP=(ZINTER(NPOINT)-ZINTER(1))/Real(NPOINT-1,Kind=pr)
    Deallocate(FUNCTI)
    !
    Allocate(FUNCTI(1:NPOINT))
    epsilon=1.D-14; one=1.0_pr
    Do i=1,NPOINT
       FUNCTI(i)=0.0_pr; QLMSUM=0.0_pr
       Do il=1,ngl
          ! Coordinates (x,y,z) for the multipole moments must be in fermis
          ZPOINT=ZINTER(i)-CENRIG ! shift with respect to the c.o.m. of left fragment
          RADIUS=Sqrt(xl(il)*bp**2+ZPOINT**2)
          ! Angle theta
          If(RADIUS.Le.epsilon) Then
             COSTHE=0.0_pr
          Else
             COSTHE=ZPOINT/RADIUS
          End If
          If(Abs(COSTHE-one).Le.epsilon) Then
             COSTHE=one-2.0_pr*epsilon
          End If
          If(Abs(COSTHE+one).Le.epsilon) Then
             COSTHE=2.0_pr*epsilon-one
          End If
          Call DEFSPH(LAMACT,COSTHE,SPHHAR)
          If(LAMACT.Eq.0) Then
             QLMVAL=1.0_pr
          Else
             If(RADIUS.Ge.epsilon) Then
                QLMVAL=RADIUS**(LAMACT)*SPHHAR(1)
             Else
                QLMVAL=0.0_pr
             End If
          End If
          DENLOC=DINTER(il,i)
          W_HERM=pi*wl(il)*bp*bp
          QLMSUM=QLMSUM+W_HERM*DENLOC*QLMVAL
       End Do
       FUNCTI(i)=QLMSUM
    End Do
    Call SIMP38(FUNCTI,NPOINT,X_STEP,RESULT)
    QLMRIG=RESULT
    !
    Deallocate(FUNCTI,XVALUE,DINTER,AUXSTO,ZINTER)
    !
  End Subroutine QLMFRA
  !====================================================================
  !> The routine computes the value of the multipole moment operators
  !> \f$ Q_{\ell 0}(\theta, \phi) = P_{\ell}(\cos\theta) \f$
  !====================================================================
  Subroutine DEFSPH(LAMACT,COSTHE,SPHHAR)
    Integer(ipr), INTENT(IN) :: LAMACT !< \f$ \ell \f$ value
    Real(pr), INTENT(IN) :: COSTHE !< \f$ \cos\theta \f$ value
    Real(pr), Dimension(1:2), INTENT(INOUT) :: SPHHAR !< Value \f$ P_{\ell}(\cos\theta) \f$

    Integer(ipr) :: i
    Real(pr) :: ZLEGPO,PIARGU,FACMUL
    Real(pr), Allocatable :: FACTOR(:)

    ! Computing P_{l,m}(cos(theta))
    ZLEGPO=DEFLEG(LAMACT,0,COSTHE)

    ! Defining the factorials
    Allocate(FACTOR(0:LAMACT))
    FACTOR(0)=1.0_pr
    Do i=1,LAMACT
       FACTOR(I)=FACTOR(i-1)*Real(i,Kind=pr)
    End Do

    ! Computing the spherical harmonics
    PIARGU=4.0_pr*Atan(1.0_pr)

    FACMUL = Sqrt(0.25_pr*Real(2*LAMACT+1,Kind=pr)/PIARGU)     &
           * q_units(LAMACT)

    SPHHAR(1) = FACMUL*ZLEGPO
    SPHHAR(2) = 0.0_pr

    Deallocate(FACTOR)

  End Subroutine DEFSPH
  !====================================================================
  !> The routine prints the characteristics of the fission fragments
  !====================================================================
  Subroutine print_moments(iw)
    Integer(ipr), Intent(in) :: iw !< Fortran unit where data is printed

    Write(iw,'(2X,"Fission fragments characteristics")')
    Write(iw,'(2X,"=================================")')
    Write(iw,'(2X,"Gaussian neck .......",2f12.4)') neckValue,Z_NECK
    Write(iw,'(2X,"Centers of mass .....",2f12.4)') CENLEF,CENRIG
    Write(iw,'(2X,"     Observable         Left fragment  Right fragment")')
    Write(iw,'(2X,"Charge Z ............",f12.4,f15.4)') QLMPRO(0,0),QLMPRO(0,1)
    Write(iw,'(2X,"Mass A ..............",f12.4,f15.4)') QLMTOT(0,0),QLMTOT(0,1)
    Write(iw,'(2X,"q10 [b^1/2] .........",2f15.6)') QLMTOT(1,0),QLMTOT(1,1)
    Write(iw,'(2X,"q20 [b] .............",2f15.6)') QLMTOT(2,0),QLMTOT(2,1)
    Write(iw,'(2X,"q30 [b^3/2] .........",2f15.6)') QLMTOT(3,0),QLMTOT(3,1)
    Write(iw,'(2X,"q40 [b^2] ...........",2f15.6)') QLMTOT(4,0),QLMTOT(4,1)
    Write(iw,'(2X,"q50 [b^5/2] .........",2f15.6)') QLMTOT(5,0),QLMTOT(5,1)
    Write(iw,'(2X,"q60 [b^3] ...........",2f15.6)') QLMTOT(6,0),QLMTOT(6,1)
    Write(iw,'(2X,"q70 [b^7/2] .........",2f15.6)') QLMTOT(7,0),QLMTOT(7,1)
    Write(iw,'(2X,"q80 [b^4] ...........",2f15.6)') QLMTOT(8,0),QLMTOT(8,1)

  End Subroutine print_moments
  !====================================================================
  !>  This function computes the probability integral at point z > 0
  !>   \f[
  !>      \Phi(z) = \frac{2}{\sqrt{\pi}} \int_{0}^{z} e^{-t^2}dt
  !>   \f]
  !>  by using the Simpson 3/8 integration rule
  !====================================================================
  Real(pr) Function PROINT(ZPOINT)
    Real(pr), INTENT(IN) :: ZPOINT !< Value of the z argument
    !
    Integer(ipr) :: i, Nx
    Real(pr) :: dx,x,PIARGU,VALRES
    Real(pr), Allocatable :: FUNCTI(:)
    !
    Nx = 100000
    Allocate(FUNCTI(1:Nx))
    !
    dx = Abs(ZPOINT/Real(Nx-1,Kind=pr))
    Do i=1,Nx
       x=dx*Real(i-1,Kind=pr); FUNCTI(i)=Exp(-x**2)
    End Do
    !
    VALRES=zero; Call SIMP38(FUNCTI,Nx,dx,VALRES)
    !
    Deallocate(FUNCTI)
    !
    PIARGU=four*Atan(one)
    PROINT=two*VALRES/Sqrt(PIARGU)
    !
  End Function PROINT
  !====================================================================
  !>  This routine performs the composite Simpson's rule integration
  !>  of a function F defined by a table of N equispaced values.
  !>  See: KOONIN, Computational physics, P.9
  !====================================================================
  Subroutine SIMP38(FUNCTI,NINTEG,X_STEP,RESULT)
    Integer(ipr), INTENT(IN) :: NINTEG !< Number of points N
    Real(pr), INTENT(IN) :: X_STEP !< Value of the uniform spacing h, \f$ h = x_{k+1} - x_{k} \f$
    Real(pr), Dimension(1:NINTEG), INTENT(IN) :: FUNCTI !< Vector of size 1:N containing the function values
    Real(pr), INTENT(INOUT) :: RESULT !< Result of the integral
    !
    Integer(ipr) :: NBEGIN,NPANEL,N_HALF,NENDIN,i
    Real(pr) :: VALSUM
    !
    ! Check to see if number of panels is even. Number of panels is n - 1.
    NBEGIN=1; NPANEL=NINTEG-NBEGIN; N_HALF=NPANEL/2; RESULT=0.0_pr
    !
    ! Number of panels is odd.  Use Simpson's 3/8 rule on first three
    ! panels, 1/3 rule on rest of them.
    If((NPANEL-2*N_HALF).Ne.0) Then
       RESULT = 3.0_pr*X_STEP*(FUNCTI(NBEGIN)               &
              + 3.0_pr*(FUNCTI(NBEGIN+1)+FUNCTI(NBEGIN+2))  &
                      + FUNCTI(NBEGIN+3))/8.0_pr
       If((NINTEG-NBEGIN).Eq.3) Return
       NBEGIN=NBEGIN+3
    End If
    !
    ! Apply 1/3 rule - add in first, second, last values
    RESULT = RESULT + X_STEP*(FUNCTI(NBEGIN)        &
                     + 4.0_pr*FUNCTI(NBEGIN+1)      &
                            + FUNCTI(NINTEG))/3.0_pr
    NBEGIN = NBEGIN+2
    !
    If(NBEGIN.Eq.NINTEG) Then
       Return
    Else
       VALSUM=0.0_pr; NENDIN=NINTEG - 1
       Do i=NBEGIN,NENDIN,2
          VALSUM=VALSUM+FUNCTI(i)+2.0_pr*FUNCTI(i+1)
       End Do
       RESULT=RESULT+2.0_pr*X_STEP*VALSUM/3.0_pr
       Return
    End If
    !
  End Subroutine SIMP38
  !====================================================================
  !>  This subroutine computes the value of the associated Legendre
  !>  polynomial \f$ P_{\ell m}(x) \f$
  !====================================================================
  Real(pr) Function DEFLEG(LAMACT,MIUACT,ZVALUE)
    Integer(ipr), INTENT(IN) :: LAMACT !< \f$ \ell \f$ value
    Integer(ipr), INTENT(IN) :: MIUACT !< m value
    Real(pr), INTENT(IN) :: ZVALUE !< x value
    !
    Integer(ipr) :: L,M
    Real(pr) :: ZVALPP,ZVAL_P,ZPOLYN,TMPVAL,epsilon,one
    !
    ! Argument z of polynomial P_{l,m}(z) must be lower than 1
    epsilon=1.D-14; one=1.0_pr
    If(Abs(Abs(ZVALUE)-one).Le.epsilon) Then
       Write(6,'("ZVALUE=",f20.16," one=",f20.16," epsilon=",f20.16)') ZVALUE,one,epsilon
       Stop 'Error in DEFLEG - ARGUMENT |z| = 1'
    End If
    ! Initialization
    ZVALPP = one; ZVAL_P = ZVALUE
    ! P{0,0}(z) = 1
    If (LAMACT.Eq.0) Then
        DEFLEG=one; Return
    End If
    ! P{1,0}(z) = 1
    If (LAMACT.Eq.1.And.MIUACT.Eq.0) Then
        DEFLEG=ZVALUE; Return
    End If
    ! P{1,1}(z) = 1
    If (LAMACT.Eq.1.And.MIUACT.Eq.1) Then
        DEFLEG=Sqrt(one-ZVALUE**2); Return
    End If
    !
    ! Obtaining P_{l,0}(z) by recurrence
    If (LAMACT.Gt.1) Then
        !
        Do L = 1,LAMACT-1
           ZPOLYN = (Real(2*L+1,Kind=pr)*ZVALUE*ZVAL_P - Real(L,Kind=pr)*ZVALPP) /Real(L+1,Kind=pr)
           TMPVAL=ZVAL_P; ZVAL_P=ZPOLYN; ZVALPP=TMPVAL
        End Do
        !
        If (MIUACT.Eq.0) Then
            DEFLEG=ZVAL_P
        Else
            ! Obtaining P_{l,1} from P_{l,0}(z) and P_{l-1,0}(z)
            !   - ZVAL_P contains P_{l,0}(z)
            !   - ZVALPP contains P_{l-1,0}(z)
            ZPOLYN =( Real(LAMACT,Kind=pr)*ZVALUE*ZVAL_P                 &
                    - Real(LAMACT,Kind=pr)*ZVALPP)/Sqrt(one - ZVALUE**2)
            TMPVAL=ZVAL_P; ZVAL_P=ZPOLYN; ZVALPP=TMPVAL
            ! Obtaining P_{l,m+1} from P_{l,m}(z) and P_{l,m-1}(z)
            !   - ZVAL_P contains P_{l,1}(z)
            !   - ZVALPP contains P_{l,0}(z)
            Do M = 1,MIUACT-1
               ZPOLYN =-2.0_pr*Real(M,Kind=pr)*ZVALUE*ZVAL_P/Sqrt(one - ZVALUE**2) &
                             - Real((LAMACT+M)*(LAMACT-M+1),Kind=pr)*ZVALPP
               TMPVAL=ZVAL_P; ZVAL_P=ZPOLYN; ZVALPP=TMPVAL
            End Do
            DEFLEG=(-1)**MIUACT * ZVAL_P
        End If
        !
    End If
    !
  End Function DEFLEG
  !====================================================================
  !>  Implementation of the Brent method to find the root of a function
  !====================================================================
  Real(pr) Function ZBRENT(FONCTI,XLOWER,XUPPER,TOLERA,IFOUND)
    Real(pr), INTENT(IN) :: XLOWER !< Lower bound for the root
    Real(pr), INTENT(IN) :: XUPPER !< Upper bound for the root
    Real(pr), INTENT(IN) :: TOLERA !< Numerical precision for the Brent method
    Integer(ipr), INTENT(INOUT) :: IFOUND !< Integer equal to 1 if minimum found
    !
    Integer(ipr), Parameter :: ITRMAX=100
    Integer(ipr) :: ITERAT
    Real(pr), Parameter :: EPSCPU=3.D-8
    Real(pr) :: A,B,C,D,E,FA,FB,FC,P,Q,R,S,EPSILO,XM
    !
    Interface
       Real(pr) Function FONCTI(X)
         Use HFBTHO
         Use HFBTHO_utilities
         Implicit None
         Real(pr), INTENT(IN) :: X
       End Function FONCTI
    End Interface
    !
    IFOUND=1; A=XLOWER; B=XUPPER; FA=FONCTI(A); FB=FONCTI(B)
    If((FA.Gt.0.0_pr.And.FB.Gt.0.0_pr).Or.(FA.Lt.0.0_pr.And.FB.Lt.0.0_pr)) Then
       IFOUND=0
       If(debug_fission.Ge.1) Then
          Write(6,'("A = ",F20.14," FA = ",f20.14)') A,FA
          Write(6,'("B = ",F20.14," FB = ",f20.14)') B,FB
          Write(6,'("ROOT MUST BE BRACKETED FOR ZBRENT")')
       End if
       Return
    End If
    !
    C=B; FC=FB
    !
    Do ITERAT=1,ITRMAX
       !
       If((FB.Gt.0.0_pr.And.FC.Gt.0.0_pr).Or.(FB.Lt.0.0_pr.And.FC.Lt.0.0_pr)) Then
         C=A; FC=FA; D=B-A; E=D
       End If
       !
       If(Abs(FC).LT.Abs(FB)) Then
         A=B; B=C; C=A; FA=FB; FB=FC; FC=FA
       End If
       !
       EPSILO=2.0_pr*EPSCPU*Abs(B) + 0.5_pr*TOLERA; XM=0.5_pr*(C - B)
       !
       If(Abs(XM).Le.EPSILO .Or. Abs(FB).Le.EPSCPU) Then
          ZBRENT=B; Return
       End If
       !
       If(Abs(E).Ge.EPSILO .And. Abs(FA).Ge.Abs(FB)) Then
           !
           S = FB/FA
           !
           If(Abs(A-C).Le.EPSCPU) Then
              P=2.0_pr*XM*S
              Q=1.0_pr - S
           Else
              ! Attempt inverse quadratic interpolation
              Q=FA/FC; R=FB/FC
              P=S*(2.0_pr*XM*Q*(Q - R) - (B - A)*(R - 1.0_pr))
              Q=(Q - 1.0_pr)*(R - 1.0_pr)*(S - 1.0_pr)
           End If
           !
           ! Check whether in bounds
           If (P.Gt.0.0_pr) Q = -Q
           !
           P=Abs(P)
           !
           ! Test quality of interpolation. If too bad, switch to bisection method
           If(2.0_pr*P .Lt. Min(3.0_pr*XM*Q-Abs(EPSILO*Q),Abs(E*Q))) Then
              E=D; D=P/Q
           Else
              D=XM; E=D
           End If
           !
       Else ! Bounds decreasing too slowly, use bisection
          !
          D=XM; E=D
          !
       End If
       !
       ! Move last best guess to A
       A=B; FA=FB
       !
       ! Evaluate new trial root
       If(Abs(D) .Gt. EPSILO) Then
          B=B+D
       Else
          B=B+Sign(EPSILO,XM)
       End If
       !
       FB=FONCTI(B)
       !
    End Do
    !
    IFOUND=2
    ZBRENT=B
    !
  End Function ZBRENT
  !====================================================================
  !>  Function returning the cubic spline interpolation at point x
  !====================================================================
  Real(pr) Function SPLINT(XENTRY,YENTRY,Y2A,NPOINT,XARGUM)
    Integer(ipr), INTENT(IN) :: NPOINT !< Number of data points N
    Real(pr), INTENT(IN) :: XARGUM !< Point x at which the function f is interpolated
    Real(pr), Dimension(1:NPOINT), INTENT(IN) :: XENTRY !< Vector of size 1:N of data alues \f$ x_{k} \f$
    Real(pr), Dimension(1:NPOINT), INTENT(IN) :: YENTRY !< Vector of size 1:N of data values \f$ y_{k} = f(x_{k}) \f$
    Real(pr), Dimension(1:NPOINT), INTENT(IN) :: Y2A !< Vector of size 1:N containing the second derivatives of f at point x

    Integer(ipr) :: KLOWER, KHIGHR, K
    Real(pr) :: HLNGTH,A,B

    KLOWER=1; KHIGHR=NPOINT
    Do While((KHIGHR-KLOWER).Gt.1)
       K=(KHIGHR+KLOWER)/2
       If(XENTRY(K).Gt.XARGUM) Then
          KHIGHR=K
       Else
          KLOWER=K
       End If
    End Do

    HLNGTH=XENTRY(KHIGHR)-XENTRY(KLOWER)

    If(HLNGTH.EQ.0) Stop 'SPLINT01'

    A=(XENTRY(KHIGHR)-XARGUM)/HLNGTH
    B=(XARGUM-XENTRY(KLOWER))/HLNGTH

    SPLINT=A*YENTRY(KLOWER)+B*YENTRY(KHIGHR) &
                      +((A**3-A)*Y2A(KLOWER) &
                      + (B**3-B)*Y2A(KHIGHR))*(HLNGTH**2)/6.0_pr

  End Function SPLINT
  !====================================================================
  !>  Subroutine calculating the vector of second derivatives needed
  !>  for the cubic spline interpolation ofthe function f
  !====================================================================
  Subroutine SPLINE(X,Y,N,YP1,YPN,Y2,IERROR)
    Integer(ipr), INTENT(IN) :: N !< Number of data points N
    Real(pr), INTENT(IN) :: YP1 !< Value of the first derivative of f at the first point \f$ x_{1} \f$
    Real(pr), INTENT(IN) :: YPN !< Value of the first derivative of f at the last point \f$ x_{N} \f$
    Real(pr), Dimension(1:N), INTENT(IN) :: X !< Vector of size 1:N of data alues \f$ x_{k} \f$
    Real(pr), Dimension(1:N), INTENT(IN) :: Y !< Vector of size 1:N of data values \f$ y_{k} = f(x_{k}) \f$
    Integer(ipr), INTENT(INOUT) :: IERROR !< Error flag is non-zero if something wrong happened
    Real(pr), Dimension(1:N), INTENT(OUT) :: Y2 !< Vector of size 1:N of second derivatives of the function f

    Integer(ipr), Parameter :: NDSPLN=200
    Integer(ipr) :: I,K
    Real(pr) :: SIG,P,QN,UN
    Real(pr), Dimension(1:NDSPLN) :: U

                              IERROR=0
    If(N.Lt.4.Or.N.Gt.NDSPLN) IERROR=1

    If(YP1.Gt.0.99D+30) Then
       Y2(1)=0.0_pr
       U(1)=0.0_pr
    Else
       Y2(1)=-0.5_pr
       U(1)=(3.0_pr/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
    End If

    Do I=2,N-1
       SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
       P=SIG*Y2(I-1)+2.0_pr
       Y2(I)=(SIG-1.0_pr)/P
       U(I)=(6.0_pr*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1)) /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
    End Do

    If(YPN.Gt.0.99D+30) Then
       QN=0.0_pr
       UN=0.0_pr
    Else
       QN=0.5_pr
       UN=(3.0_pr/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
    End If

    Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1)

    Do K=N-1,1,-1
       Y2(K)=Y2(K)*Y2(K+1)+U(K)
    End Do

  End Subroutine SPLINE
  !=======================================================================
  !
  !=======================================================================
End Module hfbtho_fission_fragments

