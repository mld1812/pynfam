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
!                      BESSEL FUNCTIONS PACKAGE                        !
!                                                                      !
! ==================================================================== !

!-------------------------------------------------------------------
!> This module provides a set of routines to compute the modified
!> Bessel functions of the first kind, with or without exponential
!> scaling.
!-------------------------------------------------------------------
!  Subroutines: - calci0 ( arg, result, jint )
!               - calci1(arg,result,jint)
!  Functions: - besei0(x)
!             - besei1(x)
!----------------------------------------------------------------------!

Module bessik

  Use HFBTHO_utilities

  Implicit None

Contains

  !=======================================================================
  !>  BESEI0 evaluates the exponentially scaled Bessel \f$ I_0(X) \f$ function.
  !>
  !>  Discussion:
  !>    This routine computes approximate values for the modified Bessel
  !>    function of the first kind of order zero multiplied by
  !>   \f$ \exp(-|x|) \f$.
  !>
  !>  Licensing:
  !>    This code is distributed under the GNU LGPL license.
  !>
  !>  Modified:
  !>    03 April 2007
  !>
  !>  Author:
  !>    Original FORTRAN77 version by William Cody.
  !>    FORTRAN90 version by John Burkardt.
  !>
  !>  Parameters:
  !>    Input, real ( kind = 8 ) X, the argument of the function.
  !>    Output, real ( kind = 8 ) BESEI0, the value of the function.
  !=======================================================================
  Function besei0(x)
    Implicit None

    Real(Kind=pr) :: besei0
    Integer(Kind=ipr) :: jint
    Real(Kind=pr) :: result,x

    jint = 2
    Call calci0(x,result,jint)
    besei0 = result

    Return
  End Function besei0
  !=======================================================================
  !> BESEI1 evaluates the exponentially scaled Bessel \f$ I_1(x) \f$
  !> function.
  !>
  !>  Discussion:
  !>    This routine computes approximate values for the modified Bessel
  !>    function of the first kind of order one multiplied by
  !>    \f$ \exp(-|x|) \f$.
  !>
  !>  Licensing:
  !>    This code is distributed under the GNU LGPL license.
  !>
  !>  Modified:
  !>    03 April 2007
  !>
  !>  Author:
  !>    Original FORTRAN77 version by William Cody.
  !>    FORTRAN90 version by John Burkardt.
  !>
  !>  Parameters:
  !>    Input, real ( kind = 8 ) X, the argument of the function.
  !>    Output, real ( kind = 8 ) BESEI1, the value of the function.
  !=======================================================================
  Function besei1(x)
    Implicit none

    Real(Kind=pr) :: besei1
    Integer(Kind=ipr) :: jint
    Real(Kind=pr) :: result,x

    jint = 2
    Call calci1 ( x, result, jint )
    besei1 = result

    Return
  End Function besei1
  !=======================================================================
  !> CALCI0 computes various I0 Bessel functions.
  !>
  !>  Discussion:
  !>    This routine computes modified Bessel functions of the first kind
  !>    and order zero, \f$ I_0(x) \f$ and \f$ \exp(-|x|)I_0(X) \f$, for
  !>    real arguments x. The main computation evaluates slightly modified
  !>    forms of minimax approximations generated by Blair and Edwards,
  !>    Chalk River (Atomic Energy of Canada Limited) Report AECL-4928,
  !>    October, 1974.
  !>
  !>  Licensing:
  !>    This code is distributed under the GNU LGPL license.
  !>
  !>  Modified:
  !>    03 April 2007
  !>
  !>  Author:
  !>    Original FORTRAN77 version by William Cody, Laura Stoltz.
  !>    FORTRAN90 version by John Burkardt.
  !>
  !>  Parameters:
  !>    - Input, real ( kind = 8 ) ARG, the argument.  If JINT = 1, then
  !>      the argument must be less than XMAX.
  !>    - Input, integer ( kind = 4 ) JINT, chooses the function to be computed.
  !>       1. \f$ I_0(x)\f$
  !>       2. \f$ \exp(-x) I_0(x) \f$
  !>    - Output, real ( kind = 8 ) RESULT, the value of the function,
  !>      which depends on the input value of JINT:
  !>       1. RESULT = \f$ I_0(x)\f$
  !>       2. RESULT = \f$ \exp(-x) I_0(x) \f$
  !=======================================================================
  subroutine calci0 ( arg, result, jint )
    Implicit None

    Real(Kind=pr) :: a,arg,b,exp40,forty
    Integer(Kind=ipr) :: i,jint
    Real(Kind=pr) :: one5,p(15),pp(8),q(5),qq(7),result,rec15
    Real(Kind=pr) :: sump,sumq,two25,x,xinf,xmax,xsmall,xx
    !  Mathematical constants
    Data one5 /15.0_pr/
    Data exp40 /2.353852668370199854d17/
    Data forty /40.0_pr/
    Data rec15 /6.6666666666666666666d-2/
    Data two25 /225.0_pr/
    !  Machine-dependent constants
    Data xsmall /5.55d-17/
    Data xinf /1.79d308/
    Data xmax /713.986d0/
    !  Coefficients for XSMALL <= ABS(ARG) < 15.0
    Data  p/-5.2487866627945699800d-18,-1.5982226675653184646d-14, &
            -2.6843448573468483278d-11,-3.0517226450451067446d-08, &
            -2.5172644670688975051d-05,-1.5453977791786851041d-02, &
            -7.0935347449210549190d+00,-2.4125195876041896775d+03, &
            -5.9545626019847898221d+05,-1.0313066708737980747d+08, &
            -1.1912746104985237192d+10,-8.4925101247114157499d+11, &
            -3.2940087627407749166d+13,-5.5050369673018427753d+14, &
            -2.2335582639474375249d+15/
    Data  q/-3.7277560179962773046d+03, 6.5158506418655165707d+06, &
            -6.5626560740833869295d+09, 3.7604188704092954661d+12, &
            -9.7087946179594019126d+14/
    !  Coefficients for 15.0 <= ABS(ARG)
    Data pp/-3.9843750000000000000d-01, 2.9205384596336793945d+00, &
            -2.4708469169133954315d+00, 4.7914889422856814203d-01, &
            -3.7384991926068969150d-03,-2.6801520353328635310d-03, &
             9.9168777670983678974d-05,-2.1877128189032726730d-06/
    Data qq/-3.1446690275135491500d+01, 8.5539563258012929600d+01, &
            -6.0228002066743340583d+01, 1.3982595353892851542d+01, &
            -1.1151759188741312645d+00, 3.2547697594819615062d-02, &
            -5.5194330231005480228d-04/

    x = Abs(arg)

    If(x <xsmall) Then

       result = one
    !
    !  XSMALL <= ABS(ARG) < 15.0.
    !
    Else If (x<one5) Then

       xx = x * x
       sump = p(1)
       Do i = 2, 15
          sump = sump * xx + p(i)
       End Do
       xx = xx - two25

       sumq = (((( xx + q(1) ) * xx + q(2) ) * xx + q(3) ) * xx + q(4) ) * xx + q(5)

       result = sump / sumq

       If(jint == 2 ) Then
          result = result * Exp(-x)
       End If

    Else If (one5<=x) Then

       If(jint==1 .and. xmax<x) Then
          result = xinf
       Else
          !
          !  15.0 <= ABS(ARG).
          !
          xx = one / x - rec15

          sump = (((((( pp(1) * xx + pp(2) ) * xx + pp(3) ) * xx + pp(4) ) &
               * xx + pp(5) ) * xx + pp(6) ) * xx + pp(7) ) * xx + pp(8)

          sumq = (((((( xx + qq(1) ) * xx + qq(2) ) * xx + qq(3) ) * xx + qq(4) ) &
                      * xx + qq(5) ) * xx + qq(6) ) * xx + qq(7)

          result = sump / sumq

          If(jint==2) Then
             result = (result - pp(1)) / Sqrt(x)
          Else
             !
             !  Calculation reformulated to avoid premature overflow.
             !
             If(x.le.(xmax-one5)) Then
                a = Exp(x)
                b = one
             Else
                a = Exp(x-forty)
                b = exp40
             End If

             result = ( (result*a - pp(1)*a) / Sqrt(x) ) * b

          End If

       End If

    End If

    Return
  End Subroutine calci0
  !=======================================================================
  !>  CALCI1 computes various I1 Bessel functions.
  !>
  !>  Discussion:
  !>    This routine computes modified Bessel functioons of the first kind
  !>    and order one, \f$ I_1(x) \f$ and \f$ \exp(-|x|) I_1(x) \f$, for real
  !>    arguments x. The main computation evaluates slightly modified forms
  !>    of minimax approximations generated by Blair and Edwards, Chalk
  !>    River (Atomic Energy of Canada Limited) Report AECL-4928,
  !>    October, 1974.
  !>
  !>  Licensing:
  !>    This code is distributed under the GNU LGPL license.
  !>
  !>  Modified:
  !>    03 April 2007
  !>
  !>  Author:
  !>    Original FORTRAN77 version by William Cody, Laura Stoltz.
  !>    FORTRAN90 version by John Burkardt.
  !>
  !>  Parameters:
  !>    - Input, real ( kind = 8 ) ARG, the argument.  If JINT = 1, then
  !>      the argument must be less than XMAX.
  !>    - Input, integer ( kind = 4 ) JINT, chooses the function to be computed.
  !>       1. \f$ I_1(x)\f$
  !>       2. \f$ \exp(-x) I_1(x) \f$
  !>    - Output, real ( kind = 8 ) RESULT, the value of the function,
  !>      which depends on the input value of JINT:
  !>       1. RESULT = \f$ I_1(x)\f$
  !>       2. RESULT = \f$ \exp(-x) I_1(x) \f$
  !=======================================================================
  Subroutine calci1(arg,result,jint)
    Implicit None

    Integer(Kind=ipr) :: j,jint
    Real(Kind=pr) :: a,arg,b,exp40,forty
    Real(Kind=pr) :: one5,p(15),pbar,pp(8),q(5),qq(6),rec15,result,sump
    Real(Kind=pr) :: sumq,two25,x,xinf,xmax,xsmall,xx
    !  Mathematical constants
    Data one5 /15.0d0/
    Data exp40 /2.353852668370199854d17/
    Data forty /40.0d0/
    Data rec15 /6.6666666666666666666d-2/
    Data two25 /225.0d0/
    !  Machine-dependent constants
    Data xsmall /5.55d-17/
    Data xinf /1.79d308/
    Data xmax /713.987d0/
    !  Coefficients for XSMALL <= ABS(ARG) < 15.0
    Data p/-1.9705291802535139930d-19,-6.5245515583151902910d-16, &
           -1.1928788903603238754d-12,-1.4831904935994647675d-09, &
           -1.3466829827635152875d-06,-9.1746443287817501309d-04, &
           -4.7207090827310162436d-01,-1.8225946631657315931d+02, &
           -5.1894091982308017540d+04,-1.0588550724769347106d+07, &
           -1.4828267606612366099d+09,-1.3357437682275493024d+11, &
           -6.9876779648010090070d+12,-1.7732037840791591320d+14, &
           -1.4577180278143463643d+15/
    Data q/-4.0076864679904189921d+03, 7.4810580356655069138d+06, &
           -8.0059518998619764991d+09, 4.8544714258273622913d+12, &
           -1.3218168307321442305d+15/
    !  Coefficients for 15.0 <= ABS(ARG)
    Data pp/-6.0437159056137600000d-02, 4.5748122901933459000d-01, &
            -4.2843766903304806403d-01, 9.7356000150886612134d-02, &
            -3.2457723974465568321d-03,-3.6395264712121795296d-04, &
             1.6258661867440836395d-05,-3.6347578404608223492d-07/
    Data qq/-3.8806586721556593450d+00, 3.2593714889036996297d+00, &
            -8.5017476463217924408d-01, 7.4212010813186530069d-02, &
            -2.2835624489492512649d-03, 3.7510433111922824643d-05/
    Data pbar/3.98437500d-01/

    x = Abs(arg)
    !
    !  Return for ABS(ARG) < XSMALL.
    !
    If(x<xsmall) Then

       result = half * x
    !
    !  XSMALL <= ABS(ARG) < 15.0.
    !
    Else If (x<one5) Then

       xx = x * x
       sump = p(1)
       Do j = 2, 15
          sump = sump * xx + p(j)
       End Do
       xx = xx - two25

       sumq = (((( xx + q(1) ) * xx + q(2) ) * xx + q(3) ) * xx + q(4) ) * xx + q(5)

       result = ( sump / sumq ) * x

       If(jint==2) Then
          result = result *Exp(-x)
       End If

    Else If (jint==1 .and. xmax<x) Then

       result = xinf

    Else
       !
       !  15.0 <= ABS(ARG).
       !
       xx = one / x - rec15

       sump = (((((( pp(1) * xx + pp(2) ) * xx + pp(3) ) * xx + pp(4) ) * xx + pp(5) ) &
                           * xx + pp(6) ) * xx + pp(7) ) * xx + pp(8)

       sumq = ((((( xx + qq(1) ) * xx + qq(2) ) * xx + qq(3) ) &
                  * xx + qq(4) ) * xx + qq(5) ) * xx + qq(6)

       result = sump / sumq

       If(jint/=1) Then
          result = (result + pbar) / Sqrt(x)
       Else
          !
          !  Calculation reformulated to avoid premature overflow.
          !
          If((xmax-one5)<x) Then
             a = Exp(x-forty)
             b = exp40
          Else
             a = Exp(x)
             b = one
          End If

          result = ( (result*a + pbar*a) / Sqrt(x) ) * b

       End If
    End If

    If(arg<zero) Then
       result = -result
    End If

    Return
  End Subroutine calci1
  !=======================================================================
  !
  !=======================================================================
End Module bessik
