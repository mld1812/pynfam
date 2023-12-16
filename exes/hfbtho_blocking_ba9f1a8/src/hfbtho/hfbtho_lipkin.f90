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
!                       LIPKIN-NOGAMI PACKAGE                          !
!                                                                      !
! ==================================================================== !

!-------------------------------------------------------------------
!> This module contains routines to compute densities corrected for the
!> Lipkin-Nogami prescription, as well as the expectation values of
!> several observables.
!----------------------------------------------------------------------
!  Subroutines: - tracesln
!               - tracesln_qp
!               - densitln
!----------------------------------------------------------------------
Module HFBTHO_Lipkin

  Use HFBTHO_utilities
  Use HFBTHO

  Implicit None

Contains
  !=======================================================================
  !> Computes the Lipkin-Nogami parameter \f$ \lambda_2 \f$ for a monopole
  !> pairing force of the type
  !>    \f[
  !>         \bar{v}_{ijkl} = -G \delta_{j\bar{i}}\delta_{l\bar{k}}
  !>                           \text{sign}(j)\text{sign}(l),
  !>    \f]
  !> The calculation is performed in the canonical basis
  !>    \f[
  !>        \lambda_{2} = \frac{G}{4}\frac{\displaystyle
  !>          \sum_{i>0} u_{i}v_{i}^{3}\sum_{i>0} u_{i}^{3}v_{i}
  !>          - \sum_{i>0} (u_{i}v_{i})^{4} }{ \displaystyle
  !>          \left( \sum_{i>0} u_{i}^{2}v_{i}^{2} \right)^{2}
  !>          - \sum_{i>0} (u_{i}v_{i})^{4} }.
  !>    \f]
  !=======================================================================
  Subroutine tracesln
    Implicit None
    Integer(ipr) :: iw,it,ib,k1,k2,kkk,k
    Real(pr) :: AAV,SNtor,SDtor
    Real(pr) :: S_U1V1,S_U1V3,S_U2V2,S_U3V1,S_U4V4
    Real(pr) :: U_ACTU,U_ACTU2,U_ACTU3,U_ACTU4
    Real(pr) :: V_ACTU,V_ACTU2,V_ACTU3,V_ACTU4
    !
    etr=zero
    Do it=itmin,itmax
       S_U1V1=ZERO; S_U1V3=ZERO; S_U2V2=ZERO; S_U3V1=ZERO; S_U4V4=ZERO
       Do ib=1,nb
          k1=ka(ib,it)+1; k2=ka(ib,it)+kd(ib,it)
          If(k1.Le.k2) Then
             kkk=lcanon(ib-1,it)
             Do k=1,id(ib)
                kkk=kkk+1; aav=vk(kkk,it)               ! v^2
                U_ACTU=Sqrt(AAV);       U_ACTU2=U_ACTU*U_ACTU
                U_ACTU3=U_ACTU2*U_ACTU; U_ACTU4=U_ACTU2*U_ACTU2
                V_ACTU=Sqrt(ONE-AAV);   V_ACTU2=V_ACTU*V_ACTU
                V_ACTU3=V_ACTU2*V_ACTU; V_ACTU4=V_ACTU2*V_ACTU2
                S_U1V1=S_U1V1+U_ACTU  * V_ACTU
                S_U1V3=S_U1V3+U_ACTU  * V_ACTU3
                S_U2V2=S_U2V2+U_ACTU2 * V_ACTU2     !Tr r (1-r)
                S_U3V1=S_U3V1+U_ACTU3 * V_ACTU
                S_U4V4=S_U4V4+U_ACTU4 * V_ACTU4     !Tr (1-r)^2 r^2
             End Do
          End If
       End Do !ib
       SNtor=8.0_pr*(S_U3V1*S_U1V3-S_U4V4)
       SDtor=32.0_pr*(S_U2V2*S_U2V2-S_U4V4)
       Geff(it)=del(it)**2/ept(it)
       ala2(it)=-Geff(it)*(SNtor/SDtor)
       If(ala2(it).Ge.10.0_pr) ala2(it)=4.0_pr ! ala2 goes to hell
       etr(it)=-four*ala2(it)*S_U2V2           ! total energy for curren particle number
    End Do !it
    etr(3)=etr(1)+etr(2) ! total energy
    Do iw=lout,lfile
       Write(iw,'(26x,a,2(1x,f7.3),a,3(1x,f9.3),a,2(1x,f7.3))')  &
            '  f: ala2(n,p)=',ala2,' #eln(n,p,t)=',etr,' #del+ala2=',del+ala2
    End Do
  End Subroutine tracesln
  !=======================================================================
  !> Computes the Lipkin-Nogami parameter \f$ \lambda_2 \f$ for a monopole
  !> pairing force of the type
  !>    \f[
  !>         \bar{v}_{ijkl} = -G \delta_{j\bar{i}}\delta_{l\bar{k}}
  !>                           \text{sign}(j)\text{sign}(l),
  !>    \f]
  !> The calculation is performed in the qp basis
  !>    \f[
  !>       \lambda_{2} = \frac{G}{4}
  !>         \frac{\text{Tr}^{<} \kappa^{*}\rho\;\text{Tr}^{<} (1-\rho)\kappa
  !>         -\sum_{ij} [ \rho(1-\rho) ]_{\bar{i}\bar{j}} [\rho(1-\rho)]_{ij}
  !>         }{ \displaystyle
  !>         \left( \text{Tr}^{>} \rho(1 - \rho) \right)^{2}
  !>         - \text{Tr}^{>} \rho^{2}(1 - \rho)^{2}.}
  !>    \f]
  !=======================================================================
  Subroutine tracesln_qp
    Implicit None
    Integer(ipr) :: iw,nd,ib,i1,i2,n2,ibitnb,it,i1n2nd,i2n2nd,i1i2nd
    Real(pr) :: frit,frit2,ftit
    Real(pr) :: etr2(2),trk(2),trk1(2),SNtor(2),SDtor(2),Sum(2)
    !
    ! initialization
    etr=zero; etr2=zero; trk=zero; trk1=zero
    ! loop over the blocks
    Do ib=1,nb
       nd=id(ib)
       ! Traces for neutrons and protons
       Do i2=1,nd ! index alpha
          Do i1=i2,nd ! index beta >= alpha
             sum=zero
             Do n2=1,nd
                i1n2nd=Max(i1,n2)+(Min(i1,n2)-1)*nd
                i2n2nd=Max(i2,n2)+(Min(i2,n2)-1)*nd
                Do it=itmin,itmax
                   ibitnb=ib+(it-1)*nbx
                   Sum(it)=Sum(it)+rk(i1n2nd,ibitnb)*rk(i2n2nd,ibitnb)*p14
                End Do !it
             End Do !n2
             i1i2nd=i1+(i2-1)*nd
             Do it=itmin,itmax
                ibitnb=ib+(it-1)*nbx
                frit=rk(i1i2nd,ibitnb)*half
                ftit=ak(i1i2nd,ibitnb)
                frit2=Sum(it)
                If(i1.Eq.i2) Then
                   etr(it)=etr(it)+frit-frit**2                  ! Tr r (1-r)
                   etr2(it)=etr2(it)+(one-two*frit+frit2)*frit2  ! Tr (1-r)^2 r^2
                   trk(it)=trk(it)+frit*ftit                     ! Tr r k
                   trk1(it)=trk1(it)+ftit -frit*ftit             ! Tr k (1-r)
                Else
                   etr(it)=etr(it) -two*frit**2                  ! Tr r (1-r)
                   etr2(it)=etr2(it)+two*(-two*frit+frit2)*frit2 ! Tr (1-r)^2 r^2
                   trk(it)=trk(it)+two*frit*ftit                 ! Tr r k
                   trk1(it)=trk1(it)-two*ftit*frit               ! Tr k (1-r)
                End If
             End Do !it
          End Do !i1
       End Do !i2
    End Do !ib
    ! total traces
    Do it=itmin,itmax
       SNtor(it)=8.0_pr*(trk1(it)*trk(it)-etr2(it))
       SDtor(it)=32.0_pr*(etr(it)**2     -etr2(it))
       Geff(it)=del(it)**2/ept(it)
       ala2(it)=-( SNtor(it)/SDtor(it) )*Geff(it)
       If(ala2(it).Ge.10.0_pr) ala2(it)=4.0_pr  ! in case ala2 goes to hell
       etr(it)=-four*ala2(it)*etr(it)           ! to total energy
    End Do
    etr(3)=etr(1)+etr(2)         !to total energy
    Do iw=lout,lfile
       Write(iw,'(26x,a,2(1x,f7.3),a,3(1x,f9.3),a,2(1x,f7.3))')  &
            '  #LN: ala2(n,p)=',ala2,' #eln(n,p,t)=',etr,' #del+ala2=',del+ala2
    End Do
  End Subroutine tracesln_qp
  !=======================================================================
  !> Calculates the densities in r-space at gauss-mesh points
  !> corrected due to Lipkin-Nogami
  !=======================================================================
  Subroutine densitln
    Implicit None
    Integer(ipr) :: ih,il,ib,nd,i0,i01,i02,n1,n2,nza,nzb,nra,nrb,nla,  &
         nlb,nsa,nsb,it,ml,ihli,k,k0(2),k00(2),k1,k2
    Real(pr) :: fr(2),vvs,vvc,ssln1,ssln2,ssln3,vks
    Real(pr) :: qla,qlb,qlab,qha,qhb,qhlab,qhab,sro
    !
    k0=0; k00=0
    ro=zero
    ! loop over the blocks
    Do ib=1,nb
       k00=k0
       nd=id(ib); i0=ia(ib)
       Do n2=1,nd
          i02=i0+n2; nzb=nz(i02); nrb=nr(i02);
          nlb=nl(i02); nsb=ns(i02)
          Do n1=1,n2
             i01=i0+n1; nza=nz(i01); nra=nr(i01)
             nla=nl(i01); nsa=ns(i01)
             k0=k00
             Do it=itmin,itmax
                k1=ka(ib,it)+1
                k2=ka(ib,it)+kd(ib,it)
                fr(it)=zero
                If(k1.Le.k2) Then
                   Do k=1,nd
                      k0(it)=k0(it)+1
                      ssln1=ssln(1,it)
                      ssln2=ssln(2,it)
                      ssln3=ssln(3,it)
                      vks=vk(k0(it),it)
                      vvc=vks
                      vvs=Abs(one-vks)
                      If(vvs.Ge.1.0d-40) Then
                         vvs=two*Sqrt(vks*vvs)   !2vu
                         vvc=vks+vvs**2*p14*ssln1*((two*vks-one)*ssln1-ssln2)/ssln3
                      End If
                      fr(it)=fr(it)+two*ddc(n2,k0(it),it)*ddc(n1,k0(it),it)*vvc
                   End Do
                   If (n1.Ne.n2) Then
                      fr(it)=two*fr(it)
                   End If
                End If
             End Do
             !---diagonal in spin
             If (nsa.Eq.nsb) Then
                ml=nla
                Do il=1,ngl
                   qla=ql (nra,ml,il);    qlb=ql (nrb,ml,il)
                   qlab=qla*qlb
                   Do ih=1,ngh
                      ihli=ih+(il-1)*ngh
                      qha=qh (nza,ih); qhb=qh (nzb,ih)
                      qhab=qha*qhb
                      qhlab=qhab*qlab; sro=qhlab
                      ro(ihli,:)=ro(ihli,:)+fr(:)*sro
                   End Do   !ih
                End Do  !il
             End If
          End Do !n2
       End Do !n1
    End Do !ib
    ! set the THO weights
    Do ihli=1,nghl
       ro(ihli,:)=ro(ihli,:)*wdcori(ihli)
    End Do
  End Subroutine densitln
  !=======================================================================
  !
  !=======================================================================
End Module HFBTHO_Lipkin