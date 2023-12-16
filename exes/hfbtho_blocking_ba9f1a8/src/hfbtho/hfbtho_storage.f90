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
!                        QRPA-pnFAM INTERFACE PACKAGE                  !
!                                                                      !
! ==================================================================== !

!----------------------------------------------------------------------
!> This module provides an interface between HFBTHO and the QRPA-pnFAM
!> code developed by the Chapel Hill group.
!>
!>  @author
!>    Mikka Mustonen, Thomas Shafer
!----------------------------------------------------------------------
!  Subroutines: - save_HFBTHO_solution(plaintext)
!----------------------------------------------------------------------
Module HFBTHO_storage
   implicit none

   ! Version number to check against
   integer, parameter :: FAM_VERSION = 10

   contains
   !=======================================================================
   ! Saves the relevant HFBTHO solution information to file for re-use by FAM.
   !=======================================================================
   subroutine save_HFBTHO_solution(plaintext)
      use UNEDF,  only: cr0 => crho, crr => cdrho, ctau, cdrho => crdr, ctj => cj, crdj, use_j2terms
      use HFBTHO, only: nb, id, REqpP, RUqpP, RVqpP, REqpN, RUqpN, RVqpN, nr, nz, nl, ns, npar,    &
                        nghl, wdcori, y_opt, fh, qhla_opt, fi1r_opt, fi1z_opt, fi2d_opt, npr, ala, &
                        alast, del, cpv0, cpv1, bet, ehfb, pwi, fn_T, fp_T, temper, entropy,       &
                        switch_on_temperature, nbx, ka, kd, KqpP, KpwiP, KqpN, KpwiN,              &
                        keyblo, blok1k2d, blo123d, bloblo, blok1k2, blo123, ntx

      implicit none

      logical, intent(in), optional :: plaintext
      logical :: save_as_text
      integer :: ifh, ierr, ib, is
      integer :: i, iexit
      integer, dimension(nghl, ntx) :: qhla_opt_t, fi1r_opt_t, fi1z_opt_t, fi2d_opt_t

      ! Option to save in plain text
      save_as_text = .false.
      if (present(plaintext)) then
         if (plaintext .eqv. .true.) then
            save_as_text = .true.
         end if
      end if

      ! Check and open file
      ifh = 77
      iexit = check_and_open(ifh,save_as_text)
      If(iexit.NE.0) Then
         write(*,'(a,i0)') ' Error saving HFB solution: could not open file to write. Error', ierr
         stop
      End If

      ! ------------------------------------------------------------------------
      ! Store the HFB solution details to 'solution.hfb(.txt if plain text)'
      ! ------------------------------------------------------------------------
      if (.not. save_as_text) then
         ! The solution version number
         write(ifh) FAM_VERSION

         ! Basic HFB quantities
         ! N, Z, A where applicable
         write(ifh) npr(:)             ! Particle number
         write(ifh) ala(:)             ! Lambdas ala (quasiparticles measured w.r.t. these, I think)
         write(ifh) alast(:)           ! Lambdas alast (last-bound s.p. energy)
         write(ifh) del(:)             ! Pairing gaps
         write(ifh) pwi                ! Pairing window
         write(ifh) CpV0(:), CpV1(:)   ! Pairing strengths
         write(ifh) bet                ! Total deformation
         write(ifh) ehfb               ! Binding energy

         ! The HFB quasiparticle energies and amplitudes
         write(ifh) nbx
         write(ifh) nb
         write(ifh) id(:)
         write(ifh) REqpP(:)
         write(ifh) RVqpP(:)
         write(ifh) RUqpP(:)
         write(ifh) REqpN(:)
         write(ifh) RVqpN(:)
         write(ifh) RUqpN(:)

         ! Pairing window active q.p. levels
         write(ifh) ka(:,:)
         write(ifh) kd(:,:)
         write(ifh) KqpP(:)
         write(ifh) KpwiP(:)
         write(ifh) KqpN(:)
         write(ifh) KpwiN(:)

         ! Basis state quantum numbers
         write(ifh) nr(:)
         write(ifh) nz(:)
         write(ifh) nl(:)
         write(ifh) ns(:)
         write(ifh) npar(:)

         ! Wave functions and integration data
         write(ifh) nghl                      ! number of integration points
         write(ifh) wdcori(:)                 ! inverse of integration weights
         write(ifh) y_opt(:)                  ! 1/rho in 'fm^(-1)'
         write(ifh) fh(:)                     ! z in 'fm'
         write(ifh) transpose(qhla_opt(:,:))  ! 1st index: state, 2nd index: coordinate (before transposing)
         write(ifh) transpose(fi1r_opt(:,:))  ! 1st index: state, 2nd index: coordinate (before transposing)
         write(ifh) transpose(fi1z_opt(:,:))  ! 1st index: state, 2nd index: coordinate (before transposing)
         write(ifh) transpose(fi2d_opt(:,:))  ! 1st index: state, 2nd index: coordinate (before transposing)

         ! Time-even isovector coupling constants from HFB mean field
         write(ifh) cr0(1), crr(1), cdrho(1), ctau(1), ctj(1), crdj(1)
         ! Store .TRUE. => 1, .FALSE. => 0
         write(ifh) merge(1, 0, use_j2terms)

         ! Temperature-dependence of the HFB solution
         ! Store .TRUE. => 1, .FALSE. => 0
         write(ifh) merge(1, 0, switch_on_temperature)
         write(ifh) temper
         write(ifh) entropy
         write(ifh) fp_T(:)
         write(ifh) fn_T(:)

         ! Neutron Blocking
         if (keyblo(1) == 0) then
            write(ifh) 0, 0, 0, 0, 0.0d0
         else
            ! Write down the block and state
            do ib=0, nb
               if (sum(id(1:ib)) >= KqpN(blok1k2d(1))) exit
            end do
            is = KqpN(blok1k2d(1)) - sum(id(1:ib-1))
            write(ifh) keyblo(1), KqpN(blok1k2d(1)), ib, is, REqpN(KqpN(blok1k2d(1)))
         end if
         ! Proton Blocking
         if (keyblo(2) == 0) then
            write(ifh) 0, 0, 0, 0, 0.0d0
         else
            do ib=0, nb
               if (sum(id(1:ib)) >= KqpP(blok1k2d(2))) exit
            end do
            is = KqpP(blok1k2d(2)) - sum(id(1:ib-1))
            write(ifh) keyblo(2), KqpP(blok1k2d(2)), ib, is, REqpP(KqpP(blok1k2d(2)))
         end if

      else

         ! The solution version number
         write(ifh,*) 'FAM_VERSION'
         write(ifh,*) FAM_VERSION
         write(ifh,*)

         ! Basic HFB quantities
         ! N, Z, A where applicable
         write(ifh,*) 'npr(:)'
         write(ifh,*) npr(:)     ! Particle number
         write(ifh,*) 'ala(:)'
         write(ifh,*) ala(:)     ! Lambdas ala (quasiparticles measured w.r.t. these, I think)
         write(ifh,*) 'alast(:)'
         write(ifh,*) alast(:)   ! Lambdas alast (last-bound s.p. energy)
         write(ifh,*) 'del(:)'
         write(ifh,*) del(:)     ! Pairing gaps
         write(ifh,*) 'pwi'
         write(ifh,*) pwi        ! Pairing window
         write(ifh,*) 'CpV0(:)'
         write(ifh,*) CpV0(:)    ! Pairing strength
         write(ifh,*) 'CpV1(:)'
         write(ifh,*) CpV1(:)    ! Pairing strength
         write(ifh,*) 'bet'
         write(ifh,*) bet        ! Total deformation
         write(ifh,*) 'ehfb'
         write(ifh,*) ehfb       ! Binding energy
         write(ifh,*)

         ! The HFB quasiparticle energies and amplitudes
         write(ifh,*) 'nbx'
         write(ifh,*) nbx
         write(ifh,*) 'nb'
         write(ifh,*) nb
         write(ifh,*) 'id(:)'
         write(ifh,*) id(:)
         write(ifh,'(2x,6(A,2x))') 'REqpP(:)', 'RVqpP(:)', 'RUqpP(:)', 'REqpN(:)', 'RVqpN(:)', 'RUqpN(:)'
         do i=0,size(REqpP)
            write(ifh,*) REqpP(i), RVqpP(i), RUqpP(i), REqpN(i), RVqpN(i), RUqpN(i)
         end do
         write(ifh,*)

         ! Pairing window active q.p. levels
         write(ifh,*) 'ka(:,:)'
         do i=1,2
            write(ifh,*) ka(:,i)
         end do
         write(ifh,*) 'kd(:,:)'
         do i=1,2
            write(ifh,*) kd(:,i)
         end do
         write(ifh,'(2x,4(A,2x))')  'KqpP(:)', 'KpwiP(:)', 'KqpN(:)', 'KpwiN(:)'
         do i=1,size(KqpP)
             write(ifh,*) KqpP(i), KpwiP(i), KqpN(i), KpwiN(i)
         end do
         write(ifh,*)

         ! Basis state quantum numbers
         write(ifh,'(2x,5(A,2x))') 'nr(:)', 'nz(:)', 'nl(:)', 'ns(:)', 'npar(:)'
         do i=1,size(nr)
             write(ifh,*) nr(i), nz(i), nl(i), ns(i), npar(i)
         end do
         write(ifh,*)

         ! Wave functions and integration data
         write(ifh,*) 'nghl'
         write(ifh,*) nghl  ! number of integration points
         write(ifh,'(2x,3(A,2x))') 'wdcori(:)', 'y_opt(:)', 'fh(:)'
         do i=1,nghl
            write(ifh,*) wdcori(i), y_opt(i), fh(i)  ! (inverse of integration weights), (1/rho in 'fm^(-1)'), (z in 'fm')
         end do
         write(ifh,*) 'trans(qhla_opt(:,:))'
         qhla_opt_t = transpose(qhla_opt(:,:))  ! 1st index: state, 2nd index: coordinate (before transposing)
         do i=1,nghl
            write(ifh,*) qhla_opt_t(i,:)
         end do
         write(ifh,*) 'trans(fi1r_opt(:,:))'
         fi1r_opt_t = transpose(fi1r_opt(:,:))  ! 1st index: state, 2nd index: coordinate (before transposing)
         do i=1,nghl
            write(ifh,*) fi1r_opt_t(i,:)
         end do
         write(ifh,*) 'trans(fi1z_opt(:,:))'
         fi1z_opt_t = transpose(fi1z_opt(:,:))  ! 1st index: state, 2nd index: coordinate (before transposing)
         do i=1,nghl
            write(ifh,*) fi1z_opt_t(i,:)
         end do
         write(ifh,*) 'trans(fi2d_opt(:,:))'
         fi2d_opt_t = transpose(fi2d_opt(:,:))  ! 1st index: state, 2nd index: coordinate (before transposing)
         do i=1,nghl
            write(ifh,*) fi2d_opt_t(i,:)
         end do
         write(ifh,*)

         ! Time-even isovector coupling constants from HFB mean field
         write(ifh,'(2x,6(A,2x))') 'cr0(1)', 'crr(1)', 'cdrho(1)', 'ctau(1)', 'ctj(1)', 'crdj(1)'
         write(ifh,*) cr0(1), crr(1), cdrho(1), ctau(1), ctj(1), crdj(1)
         ! Store .TRUE. => 1, .FALSE. => 0
         write(ifh,*) 'Use_J^2_terms?'
         write(ifh,*) merge(1, 0, use_j2terms)
         write(ifh,*)

         ! Temperature-dependence of the HFB solution
         ! Store .TRUE. => 1, .FALSE. => 0
         write(ifh,*) 'temp_on?'
         write(ifh,*) merge(1, 0, switch_on_temperature)
         write(ifh,*) 'temper'
         write(ifh,*) temper
         write(ifh,*) 'entropy'
         write(ifh,*) entropy
         write(ifh,'(2x,2(A,2x))') 'fp_T(:)', 'fn_T(:)'
         do i=1,size(fp_T)
            write(ifh,*) fp_T(i), fn_T(i)
         end do
         write(ifh,*)

         ! Neutron Blocking
         write(ifh,*) 'Neutron_Blocking'
         write(ifh,'(2x,5(A,2x))') 'keyblo(1)', 'KqpN', 'ib', 'is', 'REqpN'
         if (keyblo(1) == 0) then
            write(ifh,*) 0, 0, 0, 0, 0.0d0
         else
            ! Write down the block and state
            do ib=0, nb
               if (sum(id(1:ib)) >= KqpN(blok1k2d(1))) exit
            end do
            is = KqpN(blok1k2d(1)) - sum(id(1:ib-1))
            write(ifh, *) keyblo(1), KqpN(blok1k2d(1)), ib, is, REqpN(KqpN(blok1k2d(1)))
         end if
         ! Proton Blocking
         write(ifh,*) 'Proton_Blocking'
         write(ifh,'(2x,5(A,2x))') 'keyblo(2)', 'KqpP', 'ib', 'is', 'REqpP'
         if (keyblo(2) == 0) then
            write(ifh, *) 0, 0, 0, 0, 0.0d0
         else
            do ib=0, nb
               if (sum(id(1:ib)) >= KqpP(blok1k2d(2))) exit
            end do
            is = KqpP(blok1k2d(2)) - sum(id(1:ib-1))
            write(ifh, *) keyblo(2), KqpP(blok1k2d(2)), ib, is, REqpP(KqpP(blok1k2d(2)))
         end if
      end if

      close(ifh)
      write(*,'(a)') ' Storage completed.'
   end subroutine save_HFBTHO_solution

   !=======================================================================
   ! This function checks status of solution.hfb file (mirroring 'check_file'
   ! in hfbtho_io). If it exists and can be opened, it is opened and the
   ! function returns exit status 0. If it does not exit, a new file is opened
   ! the function returns exit status 0. Otherwise exit status is 1.
   !=======================================================================
   Integer Function check_and_open(fileunit,save_as_text)
      Implicit None
      Logical, intent(in) :: save_as_text
      Integer, intent(in) :: fileunit
      Logical :: file_exists,file_opened
      Integer :: ierr, iexit

      iexit = 0
      If (save_as_text) Then
         file_exists=.False.; inquire(file='solution.hfb.txt', exist=file_exists); ierr=0
         write(*,'(/a)') ' ### STORING HFB SOLUTION (PLAINTEXT format)'
         write(*,'(a)')  ' Filename: "solution.hfb.txt"'
         If(file_exists) Then
            file_opened=.False.; inquire(unit=fileunit, opened=file_opened)
            If(file_opened) Then
               Close(fileunit)
            End If
            open(unit=fileunit, file='solution.hfb.txt', status='old', iostat=ierr)
         Else
            open(unit=fileunit, file='solution.hfb.txt', status='new', iostat=ierr)
         End If
      Else
         file_exists=.False.; inquire(file='solution.hfb', exist=file_exists); ierr=0
         write(*,'(/a)') ' ### STORING HFB SOLUTION (BINARY format)'
         write(*,'(a)')  ' Filename: "solution.hfb"'
         If(file_exists) Then
            file_opened=.False.; inquire(unit=fileunit, opened=file_opened)
            If(file_opened) Then
               Close(fileunit)
            End If
            open(unit=fileunit, file='solution.hfb', status='old', form='unformatted', iostat=ierr)
         Else
            open(unit=fileunit, file='solution.hfb', status='new', form='unformatted', iostat=ierr)
         End If
      End If
      If(ierr.NE.0) Then
         iexit = 1
      End If
      check_and_open = iexit
   End Function check_and_open
end module HFBTHO_storage
