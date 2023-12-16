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
!                      MULITPOLE MOMENTS PACKAGE                       !
!                                                                      !
! ==================================================================== !

!-------------------------------------------------------------------
!> This module defines and computes the matrix elements and expectation
!> values of the axial multipole moments \f$ \hat{Q}_{\lambda 0} \f$
!> defined by
!>   \f[
!>        \hat{Q}_{\lambda 0}(r,\theta,\varphi)
!>        = \mathcal{N}_{\lambda}r^{\lambda} Y_{\lambda 0} (\theta,\varphi).
!>   \f]
!> Assuming axial symmetry, the multipole moments can be written
!>   \f[
!>        \hat{Q}_{\lambda 0}(r,\theta)
!>        = \mathcal{N}_{\lambda} \sqrt{\frac{2\lambda+1}{4\pi}}
!>                              r^{\lambda} P_{\lambda} (\cos\theta).
!>   \f]
!> In cylindrical coordinates, \f$ r^2 = \rho^{2} + z^2 \f$ and \f$ r\cos\theta = z \f$.
!> The multipole moments can be written as a function of \f$ \rho \f$ and \f$ z \f$ only.
!> The normalization constants \f$ \mathcal{N}_{\lambda} \f$ are the same as in the code
!> HFODD: all moments are expressed in powers of barns and the coefficients are thus
!>
!>   \f{align}{
!>        \mathcal{N}_{0} &  =\sqrt{4\pi} \\
!>        \mathcal{N}_{1} &  =\sqrt{\frac{4\pi}{3}}\frac{1}{10} \\
!>        \mathcal{N}_{2} &  =\sqrt{\frac{16\pi}{5}}\frac{1}{100} \\
!>        \mathcal{N}_{\lambda} &  =\frac{1}{10^{\lambda}},\ \ \lambda>2 \\
!>   \f}
!>
!> @author
!> Nicolas Schunck
!-------------------------------------------------------------------
!  Subroutines: - moments_setUnits
!               - moments_computeValue
!               - moments_valueMesh(z,rrr,Qval)
!               - moments_computeField(lambda,ib,debug)
!               - moments_expectation(lambda,it,ib,qval,rho,dd)
!----------------------------------------------------------------------!
Module HFBTHO_multipole_moments

  Use HFBTHO_utilities
  Use HFBTHO

  Implicit None

Contains
  !=======================================================================
  !> Defines standard units for multipole moments
  !=======================================================================
  Subroutine moments_setUnits()
    Implicit None
    Integer(ipr) :: lambda
    Real(pr) :: sqr4pi
    q_units = one

    sqr4pi=Sqrt(pp16*Atan(one))

    q_units(0)=+sqr4pi
    q_units(1)=+sqr4pi/Sqrt(three)
    q_units(2)=+sqr4pi/Sqrt(five)*two

    Do lambda=0,lambdaMax
       q_units(lambda)=q_units(lambda) / ten**lambda
    End Do

    Return
  End Subroutine moments_setUnits
  !=======================================================================
  !> Expectation value of multipole moments in coordinate space
  !=======================================================================
  Subroutine moments_computeValue()
    Implicit None
    Integer(ipr) :: lambda,ihli
    Real(pr), Dimension(0:lambdaMax) :: Qval
    Real(pr) :: sqr4pi,z,z2,z3,z4,z5,z6,z7,z8,rrr,rrr4,rrr6
    Real(pr) :: rown,rowp,whl,rn,rp
    sqr4pi=one/Sqrt(pp16*Atan(one))
    !
    qmoment=zero; Qval=zero
    !
    Do ihli=1,nghl
       !
       whl=wdcor(ihli)
       rn=ro(ihli,1); rp=ro(ihli,2)
       rown=whl*rn; rowp=whl*rp;
       z=fh(ihli); rrr=fl(ihli)**2
       !
       Call moments_valueMesh(z,rrr,Qval)
       !
       Do lambda=0,lambdaMax
          qmoment(lambda,1)=qmoment(lambda,1)+rown*Qval(lambda)
          qmoment(lambda,2)=qmoment(lambda,2)+rowp*Qval(lambda)
       End Do
       !
    End Do
    !
    Do lambda=0,lambdaMax
       qmoment(lambda,3)=qmoment(lambda,1)+qmoment(lambda,2)
    End Do

    Return
  End Subroutine moments_computeValue
  !=======================================================================
  !> Value of multipole moments at point \f$ (\rho,z)\f$ of the quadrature
  !> grid
  !=======================================================================
  Subroutine moments_valueMesh(z,rrr,Qval)
    Implicit None
    Integer(ipr) :: lambda
    Real(pr), Dimension(0:lambdaMax) :: Qval
    Real(pr) :: sqr4pi,z,z2,z3,z4,z5,z6,z7,z8,rrr,rrr4,rrr6,fq,betac,rc,ac
    !
    sqr4pi=one/Sqrt(pp16*Atan(one))
    !
    z2=z*z; z3=z2*z; z4=z3*z; z5=z4*z; z6=z5*z; z7=z6*z; z8=z7*z
    rrr4=rrr*rrr; rrr6=rrr4*rrr
    !
    Qval(0) =               sqr4pi
    Qval(1) = Sqrt(three)  *sqr4pi          * z
    Qval(2) = Sqrt(five)   *sqr4pi*half     * (two*z2-        rrr)
    Qval(3) = Sqrt(seven)  *sqr4pi*half     * (two*z3-three*z*rrr)
    Qval(4) = Sqrt(nine)   *sqr4pi*p18      * (eight*z4-24.0_pr*z2*rrr    + three *rrr4)
    Qval(5) = Sqrt(11.0_pr)*sqr4pi*p18      * (eight*z5-   pp40*z3*rrr    +pp15*z *rrr4)
    Qval(6) = Sqrt(13.0_pr)*sqr4pi/pp16     * (pp16*z6-120.0_pr*z4*rrr+ 90.0_pr*z2*rrr4-five     *rrr6)
    Qval(7) = Sqrt(15.0_pr)*sqr4pi/pp16     * (pp16*z7-168.0_pr*z5*rrr+210.0_pr*z3*rrr4-35.0_pr*z*rrr6)
    Qval(8) = Sqrt(17.0_pr)*sqr4pi/128.0_pr * (128.0_pr*z8-1792.0_pr*z6*rrr +3360.0_pr*z4*rrr4 &
                                                          -1120.0_pr*z2*rrr6+  35.0_pr*rrr4*rrr4)
    !
    If(Parity) Then
       Qval(1)=zero; Qval(3)=zero;Qval(5)=zero; Qval(7)=zero
    End If
    !
    fq = one
    If(filter) Then
       betac=1.5_pr; rc=ten; ac=one
       fq = one/(one + exp((Sqrt(z2/betac/betac + rrr) - rc)/ac))
    End If
    !
    Do lambda=0,lambdaMax
       Qval(lambda)=Qval(lambda)*q_units(lambda)*fq
    End Do
    !
    Return
  End Subroutine moments_valueMesh
  !=======================================================================
  !> Multipole moment in coordinate space (scalar fields at the Gauss
  !> quadrature points)
  !=======================================================================
  Subroutine moments_computeField(lambda,ib,multMatElems,debug)
    Implicit None
    Integer(ipr), Intent(in) :: lambda,ib
    Real(pr), Allocatable, Intent(INOUT) :: multMatElems(:)
    Logical, Optional, Intent(in) :: debug
    Integer(ipr) :: i,ih,il,nd,nd2,ihli,ihil,im,n1,n2
    Integer(ipr) :: ja,jb,nsa,nsb,nsab,ssu,ssd
    Real(pr) :: qhla,vh,fiun1,fiun2,fidn1,fidn2,vnhl
    Real(pr), Dimension(1:nghl) :: Vmom
    Real(pr), Dimension(0:lambdaMax) :: Qval
    Real(pr) :: z,rrr
    Integer(ipr) :: ndxmax
    Parameter(ndxmax=(n00max+2)*(n00max+2)/4)
    Real(pr) :: OMPFIU(ndxmax),OMPFID(ndxmax)
    !
    Qval=zero
    !
    ! Compute moment lambda on integration mesh
    Do ihli=1,nghl
       z=fh(ihli); rrr=fl(ihli)**2
       Call moments_valueMesh(z,rrr,Qval)
       Vmom(ihli)=Qval(lambda)
    End Do !ihli
    !
    ! Form matrix of the multipole constraint lambda in HO basis
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
          do n2=1,n1
             i=i+1; nsb=NS(n2+im); nsab=nsa+nsb; vh=0.0_pr
             If (nsab.Ne.0) Then
                 If (nsb.Gt.0) Then
                     fiun2 = OMPFIU(N2)
                     vh    = fiun1*fiun2
                 Else
                     fidn2 = OMPFID(N2)
                     vh    = fidn1*fidn2
                 End If
                 If (.Not.present(debug)) Then
                     multMatElems(i)=multMatElems(i)+vh*vnhl
                 End If
             End If
             If (present(debug)) Then
                 If(n2.eq.n1) multMatElems(i)=multMatElems(i)+vh
             End If
          End Do !n2
       End Do !n1
    End Do !ihil
    !
    Return
  End Subroutine moments_computeField
  !=======================================================================
  !> Expectation value of multipole moments in configuration space
  !=======================================================================
  Subroutine moments_expectation(lambda,it,ib,qval,rho,dd,multMatElems)
    Implicit None
    Integer(ipr), Intent(In) :: lambda,it,ib
    Real(pr), Allocatable, Intent(In) :: rho(:,:)
    Real(pr), Intent(Inout) :: qval,dd
    Real(pr), Allocatable, Intent(INOUT) :: multMatElems(:)
    !
    Integer(ipr) :: nd,nd2,nhfb,i0,m,j,n1,n2,i
    Real(pr) :: hla
    Real(pr), Allocatable :: dblmul(:,:),qblock(:,:)
    !
    nd=id(ib); nd2=nd*nd; nhfb=nd+nd; i0=ia(ib); m=ib+(it-1)*nbx
    multMatElems=zero
    Call moments_computeField(lambda,ib,multMatElems)
    ! Trace of the density should give particle number
    Do n1=1,nd
       dd=dd+rho(n1,n1)
    End Do
    ! Expectation value of multipole moment in configuration space
    Allocate(dblmul(nd,nd));dblmul=zero
    j=0
    Do n1=1,nd
       Do n2=1,n1
          j=j+1;hla=multMatElems(j)
          dblmul(n1,n2)=hla;dblmul(n2,n1)=hla
       End Do
    End Do
    Allocate(qblock(nd,nd))
    Call dgemm('n','n',nd,nd,nd,one,rho,nd,dblmul,nd,zero,qblock,nd)
    Do n1=1,nd
       qval=qval+qblock(n1,n1)
    End Do
    !
    Return
  End Subroutine moments_expectation
  !=======================================================================
  !
  !=======================================================================
End Module HFBTHO_multipole_moments
