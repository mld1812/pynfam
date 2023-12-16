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
!                      LOCALIZATION PACKAGE                            !
!                                                                      !
! ==================================================================== !

!-------------------------------------------------------------------
!> This module provides a set of routines to compute the localization
!> functions
!>
!> @author
!> Chunli Zhang
!-------------------------------------------------------------------
!  Subroutines: - localization
!----------------------------------------------------------------------!

Module HFBTHO_localization

  Use HFBTHO_utilities
  Use HFBTHO

  Implicit None

 Contains
  !=======================================================================
  !
  !=======================================================================
  Subroutine localization
    Real(pr), Parameter :: p35=three/five
    Real(pr), Parameter :: pi=four*Atan(one)
    !
    Integer(ipr) :: irz
    Real(pr), Allocatable :: tautf(:,:),loca(:,:)
    !
    If(Allocated(tautf)) Deallocate(tautf,loca)
    Allocate(tautf(nghl,2),loca(nghl,2))
    !
    Do irz=1,nghl
      tautf(irz,1)=p35*(six*pi**2)**p23 * ro(irz,1)**p53
      tautf(irz,2)=p35*(six*pi**2)**p23 * ro(irz,2)**p53

      ! All densities need to be divided by 2, which is given by p12**p23.
      loca(irz,1)=one  &
        /(one+((ro(irz,1)*tau(irz,1)-p14*(NABLAR(irz,1)**2+NABLAZ(irz,1)**2))  &
        /(ro(irz,1)*tautf(irz,1)*p12**p23))**2)

      loca(irz,2)=one  &
        /(one+((ro(irz,2)*tau(irz,2)-p14*(NABLAR(irz,2)**2+NABLAZ(irz,2)**2))  &
        /(ro(irz,2)*tautf(irz,2)*p12**p23))**2)
    End Do
    !
    Open(unit=300,file='localization.dat')
    !
    Write(300,*) 'r  z  rhoN/2  rhoP/2  localizarionN  localizationP '
    !
    Do irz=1,nghl
      ! Loca(:,1): neutron NLF for fixed spin; loca(:,2): proton NLF for fixed
      ! spin
      Write(300,"(6(1x,F16.8))") fl(irz),fh(irz),ro(irz,1)/two,ro(irz,2)/two,loca(irz,1),loca(irz,2)
    End Do
    !
    Close(300)
    !
    Deallocate(tautf,loca)
    !
    Return
  End Subroutine localization

End Module HFBTHO_localization
!==================================================================================================================================
!#END HFBTHO_localization
!==================================================================================================================================
