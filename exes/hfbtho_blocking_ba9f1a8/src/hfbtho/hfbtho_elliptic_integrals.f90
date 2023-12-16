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
!                      ELLIPTIC INTEGRAL PACKAGE                       !
!                                                                      !
! ==================================================================== !

!-------------------------------------------------------------------
!> This module provides a routine to compute the complete elliptic
!> integral of the second kind.
!>
!> Reference:
!>    Fukushima, T., "Fast Computation of Complete Elliptic Integrals
!>    and Jacobian Elliptic Functions," Celest. Mech. Dyn. Astron. 105,
!>    305-328 (2009b)
!-------------------------------------------------------------------
!  Subroutines: - CompleteEllipticFunction_2nd
!               - auxiliary(x,Emp,Kmp)
!               - elliptic_small_m(x)
!  Functions: - nome(x)
!----------------------------------------------------------------------!

Module EllipticIntegral

  Use HFBTHO_utilities

  Implicit None

Contains

  Real(pr) Function CompleteEllipticFunction_2nd(x)
    Use HFBTHO_utilities
    Implicit None
    Real(pr), INTENT(IN) :: x

    Real(pr) :: Emp,Kmp,Em,Km,qp,pi,x_eff

    pi = four*Atan(one)

    If(x.Lt.zero.Or.x.Gt.one) Stop 'Error in CompleteEllipticFunction_2nd'

    If(x.Lt.0.9_pr) Then
       Em = elliptic_small_m(x)
    Else
       Call auxiliary(x,Emp,Kmp)
       x_eff = one-x; qp = nome(x_eff)
       Km = -Log(qp)*Kmp/pi
       Em = Km + (half*pi - Emp*Km)/Kmp
    End If

    CompleteEllipticFunction_2nd = Em

  End Function CompleteEllipticFunction_2nd
  !=======================================================================
  !
  !=======================================================================
  Subroutine auxiliary(x,Emp,Kmp)
    Use HFBTHO_utilities
    Implicit None
    Real(pr), INTENT(IN) :: x
    Real(pr), INTENT(INOUT) :: Emp,Kmp

    Integer(ipr) :: JE, JK
    Parameter (JE=16, JK=16)
    Integer(ipr) :: i
    Real(pr) :: x0, x_eff
    Real(pr), Dimension(0:JE) :: Coeff_E
    Real(pr), Dimension(0:JK) :: Coeff_K

    x0 = 0.05_pr; x_eff = one - x

    Coeff_E(0) = +1.550973351780472328_pr
    Coeff_E(1) = -0.400301020103198524_pr
    Coeff_E(2) = -0.078498619442941939_pr
    Coeff_E(3) = -0.034318853117591992_pr
    Coeff_E(4) = -0.019718043317365499_pr
    Coeff_E(5) = -0.013059507731993309_pr
    Coeff_E(6) = -0.009442372874146547_pr
    Coeff_E(7) = -0.007246728512402157_pr
    Coeff_E(8) = -0.005807424012956090_pr
    Coeff_E(9) = -0.004809187786009338_pr
    Coeff_E(10)=  0.000000000000000000_pr
    Coeff_E(11)=  0.000000000000000000_pr
    Coeff_E(12)=  0.000000000000000000_pr
    Coeff_E(13)=  0.000000000000000000_pr
    Coeff_E(14)=  0.000000000000000000_pr
    Coeff_E(15)=  0.000000000000000000_pr
    Coeff_E(16)=  0.000000000000000000_pr

    Emp = 0.0_pr
    Do i=0,JE
       Emp = Emp + Coeff_E(i)*(x_eff - x0)**i
    End Do

    Coeff_K(0) = 1.591003453790792180_pr
    Coeff_K(1) = 0.416000743991786912_pr
    Coeff_K(2) = 0.245791514264103415_pr
    Coeff_K(3) = 0.179481482914906162_pr
    Coeff_K(4) = 0.144556057087555150_pr
    Coeff_K(5) = 0.123200993312427711_pr
    Coeff_K(6) = 0.108938811574293531_pr
    Coeff_K(7) = 0.098853409871592910_pr
    Coeff_K(8) = 0.091439629201749751_pr
    Coeff_K(9) = 0.085842591595413900_pr
    Coeff_K(10)= 0.081541118718303215_pr
    Coeff_K(11)= 0.000000000000000000_pr
    Coeff_K(12)= 0.000000000000000000_pr
    Coeff_K(13)= 0.000000000000000000_pr
    Coeff_K(14)= 0.000000000000000000_pr
    Coeff_K(15)= 0.000000000000000000_pr
    Coeff_K(16)= 0.000000000000000000_pr

    Kmp = 0.0_pr
    Do i=0,JK
       Kmp = Kmp + Coeff_K(i)*(x_eff - x0)**i
    End Do

  End Subroutine auxiliary
  !=======================================================================
  !
  !=======================================================================
  Real(pr) Function elliptic_small_m(x)
    Use HFBTHO_utilities
    Implicit None
    Real(pr), INTENT(IN) :: x

    Integer(ipr) :: JE
    Parameter (JE=16)
    Integer(ipr) :: i
    Real(pr) :: x0, Em
    Real(pr), Dimension(0:JE) :: Coeff_E

    If(x.Lt.0.1_pr) Then
       x0 = 0.05_pr
       Coeff_E(0) = +1.550973351780472328_pr
       Coeff_E(1) = -0.400301020103198524_pr
       Coeff_E(2) = -0.078498619442941939_pr
       Coeff_E(3) = -0.034318853117591992_pr
       Coeff_E(4) = -0.019718043317365499_pr
       Coeff_E(5) = -0.013059507731993309_pr
       Coeff_E(6) = -0.009442372874146547_pr
       Coeff_E(7) = -0.007246728512402157_pr
       Coeff_E(8) = -0.005807424012956090_pr
       Coeff_E(9) = -0.004809187786009338_pr
       Coeff_E(10)=  0.000000000000000000_pr
       Coeff_E(11)=  0.000000000000000000_pr
       Coeff_E(12)=  0.000000000000000000_pr
       Coeff_E(13)=  0.000000000000000000_pr
       Coeff_E(14)=  0.000000000000000000_pr
       Coeff_E(15)=  0.000000000000000000_pr
       Coeff_E(16)=  0.000000000000000000_pr
    End If

    If(x.Ge.0.1_pr.And.x.Lt.0.2_pr) Then
       x0 = 0.15_pr
       Coeff_E(0) = +1.510121832092819728_pr
       Coeff_E(1) = -0.417116333905867549_pr
       Coeff_E(2) = -0.090123820404774569_pr
       Coeff_E(3) = -0.043729944019084312_pr
       Coeff_E(4) = -0.027965493064761785_pr
       Coeff_E(5) = -0.020644781177568105_pr
       Coeff_E(6) = -0.016650786739707238_pr
       Coeff_E(7) = -0.014261960828842520_pr
       Coeff_E(8) = -0.012759847429264803_pr
       Coeff_E(9) = -0.011799303775587354_pr
       Coeff_E(10)= -0.011197445703074968_pr
       Coeff_E(11)=  0.000000000000000000_pr
       Coeff_E(12)=  0.000000000000000000_pr
       Coeff_E(13)=  0.000000000000000000_pr
       Coeff_E(14)=  0.000000000000000000_pr
       Coeff_E(15)=  0.000000000000000000_pr
       Coeff_E(16)=  0.000000000000000000_pr
    End If

    If(x.Ge.0.2_pr.And.x.Lt.0.3_pr) Then
       x0 = 0.25_pr
       Coeff_E(0) = +1.467462209339427155_pr
       Coeff_E(1) = -0.436576290946337775_pr
       Coeff_E(2) = -0.105155557666942554_pr
       Coeff_E(3) = -0.057371843593241730_pr
       Coeff_E(4) = -0.041391627727340220_pr
       Coeff_E(5) = -0.034527728505280841_pr
       Coeff_E(6) = -0.031495443512532783_pr
       Coeff_E(7) = -0.030527000890325277_pr
       Coeff_E(8) = -0.030916984019238900_pr
       Coeff_E(9) = -0.032371395314758122_pr
       Coeff_E(10)= -0.034789960386404158_pr
       Coeff_E(11)=  0.000000000000000000_pr
       Coeff_E(12)=  0.000000000000000000_pr
       Coeff_E(13)=  0.000000000000000000_pr
       Coeff_E(14)=  0.000000000000000000_pr
       Coeff_E(15)=  0.000000000000000000_pr
       Coeff_E(16)=  0.000000000000000000_pr
    End If

    If(x.Ge.0.3_pr.And.x.Lt.0.4_pr) Then
       x0 = 0.35_pr
       Coeff_E(0) = +1.422691133490879171_pr
       Coeff_E(1) = -0.459513519621048674_pr
       Coeff_E(2) = -0.125250539822061878_pr
       Coeff_E(3) = -0.078138545094409477_pr
       Coeff_E(4) = -0.064714278472050002_pr
       Coeff_E(5) = -0.062084339131730311_pr
       Coeff_E(6) = -0.065197032815572477_pr
       Coeff_E(7) = -0.072793895362578779_pr
       Coeff_E(8) = -0.084959075171781003_pr
       Coeff_E(9) = -0.102539850131045997_pr
       Coeff_E(10)= -0.127053585157696036_pr
       Coeff_E(11)= -0.160791120691274606_pr
       Coeff_E(12)=  0.000000000000000000_pr
       Coeff_E(13)=  0.000000000000000000_pr
       Coeff_E(14)=  0.000000000000000000_pr
       Coeff_E(15)=  0.000000000000000000_pr
       Coeff_E(16)=  0.000000000000000000_pr
    End If

    If(x.Ge.0.4_pr.And.x.Lt.0.5_pr) Then
       x0 = 0.45_pr
       Coeff_E(0) = +1.375401971871116291_pr
       Coeff_E(1) = -0.487202183273184837_pr
       Coeff_E(2) = -0.153311701348540228_pr
       Coeff_E(3) = -0.111849444917027833_pr
       Coeff_E(4) = -0.108840952523135768_pr
       Coeff_E(5) = -0.122954223120269076_pr
       Coeff_E(6) = -0.152217163962035047_pr
       Coeff_E(7) = -0.200495323642697339_pr
       Coeff_E(8) = -0.276174333067751758_pr
       Coeff_E(9) = -0.393513114304375851_pr
       Coeff_E(10)= -0.575754406027879147_pr
       Coeff_E(11)= -0.860523235727239756_pr
       Coeff_E(12)= -1.308833205758540162_pr
       Coeff_E(13)=  0.000000000000000000_pr
       Coeff_E(14)=  0.000000000000000000_pr
       Coeff_E(15)=  0.000000000000000000_pr
       Coeff_E(16)=  0.000000000000000000_pr
    End If

    If(x.Ge.0.5_pr.And.x.Lt.0.6_pr) Then
       x0 = 0.55_pr
       Coeff_E(0) = +1.325024497958230082_pr
       Coeff_E(1) = -0.521727647557566767_pr
       Coeff_E(2) = -0.194906430482126213_pr
       Coeff_E(3) = -0.171623726822011264_pr
       Coeff_E(4) = -0.202754652926419141_pr
       Coeff_E(5) = -0.278798953118534762_pr
       Coeff_E(6) = -0.420698457281005762_pr
       Coeff_E(7) = -0.675948400853106021_pr
       Coeff_E(8) = -1.136343121839229244_pr
       Coeff_E(9) = -1.976721143954398261_pr
       Coeff_E(10)= -3.531696773095722506_pr
       Coeff_E(11)= -6.446753640156048150_pr
       Coeff_E(12)= -11.97703130208884026_pr
       Coeff_E(13)=  0.000000000000000000_pr
       Coeff_E(14)=  0.000000000000000000_pr
       Coeff_E(15)=  0.000000000000000000_pr
       Coeff_E(16)=  0.000000000000000000_pr
    End If

    If(x.Ge.0.6_pr.And.x.Lt.0.7_pr) Then
       x0 = 0.65_pr
       Coeff_E(0) = +1.270707479650149744_pr
       Coeff_E(1) = -0.566839168287866583_pr
       Coeff_E(2) = -0.262160793432492598_pr
       Coeff_E(3) = -0.292244173533077419_pr
       Coeff_E(4) = -0.440397840850423189_pr
       Coeff_E(5) = -0.774947641381397458_pr
       Coeff_E(6) = -1.498870837987561088_pr
       Coeff_E(7) = -3.089708310445186667_pr
       Coeff_E(8) = -6.667595903381001064_pr
       Coeff_E(9) = -14.89436036517319078_pr
       Coeff_E(10)= -34.18120574251449024_pr
       Coeff_E(11)= -80.15895841905397306_pr
       Coeff_E(12)= -191.3489480762984920_pr
       Coeff_E(13)= -463.5938853480342030_pr
       Coeff_E(14)= -1137.380822169360061_pr
       Coeff_E(15)=  0.000000000000000000_pr
       Coeff_E(16)=  0.000000000000000000_pr
    End If

    If(x.Ge.0.7_pr.And.x.Lt.0.8_pr) Then
       x0 = 0.75_pr
       Coeff_E(0) = +1.211056027568459525_pr
       Coeff_E(1) = -0.630306413287455807_pr
       Coeff_E(2) = -0.387166409520669145_pr
       Coeff_E(3) = -0.592278235311934603_pr
       Coeff_E(4) = -1.237555584513049844_pr
       Coeff_E(5) = -3.032056661745247199_pr
       Coeff_E(6) = -8.181688221573590762_pr
       Coeff_E(7) = -23.55507217389693250_pr
       Coeff_E(8) = -71.04099935893064956_pr
       Coeff_E(9) = -221.8796853192349888_pr
       Coeff_E(10)= -712.1364793277635425_pr
       Coeff_E(11)= -2336.125331440396407_pr
       Coeff_E(12)= -7801.945954775964673_pr
       Coeff_E(13)= -26448.19586059191933_pr
       Coeff_E(14)= -90799.48341621365251_pr
       Coeff_E(15)= -315126.0406449163424_pr
       Coeff_E(16)= -1104011.344311591159_pr
    End If

    If(x.Ge.0.8_pr.And.x.Lt.0.85_pr) Then
       x0 = 0.825_pr
       Coeff_E(0) = +1.161307152196282836_pr
       Coeff_E(1) = -0.701100284555289548_pr
       Coeff_E(2) = -0.580551474465437362_pr
       Coeff_E(3) = -1.243693061077786614_pr
       Coeff_E(4) = -3.679383613496634879_pr
       Coeff_E(5) = -12.81590924337895775_pr
       Coeff_E(6) = -49.25672530759985272_pr
       Coeff_E(7) = -202.1818735434090269_pr
       Coeff_E(8) = -869.8602699308701437_pr
       Coeff_E(9) = -3877.005847313289571_pr
       Coeff_E(10)= -17761.70710170939814_pr
       Coeff_E(11)= -83182.69029154232061_pr
       Coeff_E(12)= -396650.4505013548170_pr
       Coeff_E(13)= -1920033.413682634405_pr
       Coeff_E(14)=  0.000000000000000000_pr
       Coeff_E(15)=  0.000000000000000000_pr
       Coeff_E(16)=  0.000000000000000000_pr
    End If

    If(x.Ge.0.85_pr.And.x.Lt.0.9_pr) Then
       x0 = 0.875_pr
       Coeff_E(0) = +1.124617325119752213_pr
       Coeff_E(1) = -0.770845056360909542_pr
       Coeff_E(2) = -0.844794053644911362_pr
       Coeff_E(3) = -2.490097309450394453_pr
       Coeff_E(4) = -10.23971741154384360_pr
       Coeff_E(5) = -49.74900546551479866_pr
       Coeff_E(6) = -267.0986675195705196_pr
       Coeff_E(7) = -1532.665883825229947_pr
       Coeff_E(8) = -9222.313478526091951_pr
       Coeff_E(9) = -57502.51612140314030_pr
       Coeff_E(10)= -368596.1167416106063_pr
       Coeff_E(11)= -2415611.088701091428_pr
       Coeff_E(12)= -16120097.81581656797_pr
       Coeff_E(13)= -109209938.5203089915_pr
       Coeff_E(14)= -749380758.1942496220_pr
       Coeff_E(15)= -5198725846.725541393_pr
       Coeff_E(16)= -36409256888.12139973_pr
    End If

    Em = 0.0_pr
    Do i=0,JE
       Em = Em + Coeff_E(i)*(x - x0)**i
    End Do

    elliptic_small_m = Em

  End Function elliptic_small_m
  !=======================================================================
  !
  !=======================================================================
  Real(pr) Function nome(x)
    Use HFBTHO_utilities
    Implicit None
    Real(pr), INTENT(IN) :: x
    Integer(ipr) :: Jq
    Parameter (Jq = 14)
    Integer(ipr) :: i
    Real(pr) :: qp,epsilon
    Real(pr), Dimension(1:Jq) :: Coeff_q

    epsilon = 1.D-14

    Coeff_q(1)  = 1.0_pr/16.0_pr
    Coeff_q(2)  = 1.0_pr/32.0_pr
    Coeff_q(3)  = 21.0_pr/1024.0_pr
    Coeff_q(4)  = 31.0_pr/2048.0_pr
    Coeff_q(5)  = 6257.0_pr/524288.0_pr
    Coeff_q(6)  = 10293.0_pr/1048576.0_pr
    Coeff_q(7)  = 279025.0_pr/33554432.0_pr
    Coeff_q(8)  = 483127.0_pr/67108864.0_pr
    Coeff_q(9)  = 435506703.0_pr/68719476736.0_pr
    Coeff_q(10) = 776957575.0_pr/137438953472.0_pr
    Coeff_q(11) = 22417045555.0_pr/4398046511104.0_pr
    Coeff_q(12) = 40784671953.0_pr/8796093022208.0_pr
    Coeff_q(13) = 9569130097211.0_pr/2251799813685248.0_pr
    Coeff_q(14) = 17652604545791.0_pr/4503599627370496.0_pr

    qp = 0.0_pr
    If(x.Gt.epsilon) Then
       Do i=1,Jq
          qp = qp + Coeff_q(i) * (x**i)
       End Do
    Else
       qp = epsilon
    End If

    nome = qp

  End Function nome
  !=======================================================================
  !
  !=======================================================================
End Module EllipticIntegral
