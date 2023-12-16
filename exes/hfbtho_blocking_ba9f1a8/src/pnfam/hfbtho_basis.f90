!------------------------------------------------------------------------------
! hfbtho_basis.f90
!
! This module contains the arrays holding the needed information about the
! [transformed] harmonic oscillator basis of the HFBTHO solution, as well as
! related subroutines.
!
! M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
!------------------------------------------------------------------------------
module hfbtho_basis
   use blockmatrix_type
   implicit none
   public

   integer, parameter, private :: dp = kind(1d0)

   integer :: nghl  ! number of 2D quadrature points in coordinate space

   ! The quantum numbers of the THO basis states:
   !   nr, nz = number of nodes in r, z direction (where r = distance from intrinsic z axis)
   !   nl = orbital angular momentum projection on the intrinsic axis
   !   ns = spin projection on the intrinsic axis *2 (so that it's an integer)
   !   npar = parity (1 for even? 2 for odd? CHECK THIS)
   integer, allocatable :: nr(:), nz(:), nl(:), ns(:), npar(:)

   ! The basis wave functions, their partial derivatives in the r and z direction,
   ! and their Laplacians, with quadrature weights mixed in to allow fast integration
   real(dp), allocatable :: wf(:,:), wfdr(:,:), wfdz(:,:), wfd2(:,:)

   ! The 1/rho and z in each quadrature point, and the inverse of quadrature weights times Jacobian
   real(dp), allocatable :: y(:), z(:), wdcori(:)

   ! The variables containing the block structure of the HFB solution
   integer :: hfb_nbx = -1         ! the number of blocks used to allocate ID(:) by HFBTHO
   !integer :: nb                  ! the number of blocks
   !integer, allocatable :: db(:)  ! the number of states in each block
   ! The following are useful in the general-K case, when we might not traverse
   ! the blocks in the order: imstart(i) tells the first index of the block i in
   ! the HFB matrices, isstart(i) the index of the first state of the block i
   integer, allocatable :: imstart(:), isstart(:)

   ! General HFB quantities
   ! In all cases, 1 => neutron, 2 => proton, 3 => total
   integer  :: hfb_npr(3) = -1
   real(dp) :: hfb_lambda(2) = 0, hfb_elast(2) = 0, hfb_delta(2) = 0
   real(dp) :: hfb_beta2 = 0, hfb_energy = 0

   ! HFB pairing window properties
   real(dp) :: hfb_pairing_window = 0
   integer, allocatable, target :: pwi_start_n(:), pwi_start_p(:), pwi_dim_n(:), pwi_dim_p(:), &
      pwi_qp_n(:), pwi_qp_p(:), pwi_uv_n(:), pwi_uv_p(:)
   logical, allocatable, target :: pwi_active_n(:), pwi_active_p(:)

   ! HFB T=1 pairing
   real(dp) :: hfb_cpair(2) = 0, hfb_alpha_pair(2) = 0

   ! Coupling constants and properties of the pn interaction
   real(dp) :: hfb_cr0 = 0, hfb_crr = 0, hfb_cdrho = 0, hfb_ctau = 0, hfb_ctj = 0, hfb_crdj = 0
   logical  :: hfb_use_j2terms = .false.

   ! Temperature-dependent quantities from finite-T HFB
   ! ft_active:    is T>0 active in HFB calculation?
   ! ft_temp:      temperature in the HFB (MeV)
   ! ft_entropy:   entropy of the calculation (neutron, proton, total)
   ! ft_fp, ft_fn: q.p. occupations as a function of temperature
   ! IMPORTANT: coming from HFBTHO, ft_fp and ft_fn are NOT ordered in terms of
   ! quasiparticle indices (cf. Eqp for a counterexample). The ft_fq's are
   ! ordered in terms of pwi_qp_q's --- the qp levels below-cutoff!
   logical  :: ft_active = .false.
   real(dp) :: ft_temp = 0
   real(dp) :: ft_entropy(3) = 0
   real(dp), allocatable, target :: ft_fp(:), ft_fn(:)

contains

   !---------------------------------------------------------------------------
   ! Load the HFBTHO solution and return the most important arrays
   !---------------------------------------------------------------------------
   subroutine get_HFBTHO_solution(fn, Ep, Vp, Up, En, Vn, Un)
      implicit none

      ! Minimum version of the binary solution
      integer, parameter :: minimum_version = 9

      character(len=*), intent(in) :: fn   ! input file name
      real(dp), allocatable, intent(out) :: Ep(:), En(:)   ! quasiparticle energies for protons and neutrons
      type(blockmatrix), intent(out) :: Vp, Up, Vn, Un   ! U and V matrices

      integer  :: i, ierr, iqp, iblock, ic, im, ik
      integer  :: f           ! file number
      integer  :: version     ! version number
      integer  :: dqp, dmat   ! dimension of the basis, number of non-trivial matrix elements in U, V, etc.
      integer  :: itmp        ! temporary storage for logical variables
      integer  :: ncut        ! cutoff for q.p. states below pwi
      logical  :: kactive     ! active q.p. level

      f = 13
      open(unit=f, file=fn, status='old', form='unformatted', iostat=ierr)

      if (ierr /= 0) then
         write(*,'(/3A)') 'ERROR in get_HFBTHO_solution: could not open file "', trim(fn), '".'
         stop
      end if

      ! Check version number of the solution file
      read(f) version
      if (version < minimum_version) then
         write(*,'(2(a,1i0))') 'Error in get_HFBTHO_solution(): expected version >= ', &
               minimum_version, ', got ', version
         write(*,'(a,a,a)') 'Was attempting to read in the file "', trim(fn), '"'
         stop
      end if

      !-------------------------------------------------------------------------
      ! Basic results
      !-------------------------------------------------------------------------
      read(f) hfb_npr(:)                           ! Particle number
      read(f) hfb_lambda(:)                        ! Lambdas ala (QEPs measured w.r.t. these)
      read(f) hfb_elast(:)                         ! Lambdas alast (last-bound s.p. energy)
      read(f) hfb_delta(:)                         ! Pairing gaps
      read(f) hfb_pairing_window                   ! Pairing window
      read(f) hfb_cpair(:), hfb_alpha_pair(:)      ! Pairing strengths
      read(f) hfb_beta2                            ! Total deformation
      read(f) hfb_energy                           ! Binding energy

      !-------------------------------------------------------------------------
      ! Matrices
      !-------------------------------------------------------------------------
      ! V. 9 introduces a check of nbx, which HFBTHO uses to allocate these arrays
      if (version >= 9) then
         read(f) hfb_nbx
         read(f) nb
         if (nb /= hfb_nbx) stop 'NB /= NBX. Block arrays may be misshapen!'
      else
         read(f) nb
      end if

      nb = 2*nb
      allocate(db(nb))
      read(f) db(1:nb/2)
      db(nb/2+1:) = db(1:nb/2)
      dqp = sum(db)      ! number of basis states (or, equivalently, qp states)
      dmat = sum(db*db)  ! number of nontrivial matrix elements in U, V
      allocate(Ep(dqp), En(dqp))
      call allocate_blockmatrix(Up, dmat)
      call allocate_blockmatrix(Vp, dmat)
      call allocate_blockmatrix(Un, dmat)
      call allocate_blockmatrix(Vn, dmat)

      read(f) Ep(1:dqp/2)        ; Ep(dqp/2+1:)       =  Ep(1:dqp/2)
      read(f) Vp%elem(dmat/2+1:) ; Vp%elem(1:dmat/2)  = -Vp%elem(dmat/2+1:)
      read(f) Up%elem(1:dmat/2)  ; Up%elem(dmat/2+1:) =  Up%elem(1:dmat/2)
      read(f) En(1:dqp/2)        ; En(dqp/2+1:)       =  En(1:dqp/2)
      read(f) Vn%elem(dmat/2+1:) ; Vn%elem(1:dmat/2)  = -Vn%elem(dmat/2+1:)
      read(f) Un%elem(1:dmat/2)  ; Un%elem(dmat/2+1:) =  Un%elem(1:dmat/2)

      !-------------------------------------------------------------------------
      ! V. 9 introduces the HFBTHO PWI cutoffs into all channels
      !-------------------------------------------------------------------------
      ! Documentation ---
      ! HFBTHO counts q.p. levels for each q=(n,p) below the pairing window
      ! (pwi).  These matrices connect those indices (kl in HFBTHO) to the
      ! "true" q.p. indices.
      ! ---
      ! pwi_start_q: ka(ib,it)+1 in HFBTHO
      !     Index of 1st q.p. below pwi in this block (ka is last q.p. in the
      !     previous block)
      ! pwi_dim_q: kd(ib,it) in HFBTHO
      !     Number of q.p. levels below pwi in this block. Typical iteration is
      !     then "do k=pwi_start_q(ib),pwi_start_q(ib)+pwi_dim_q(ib)-1"
      ! pwi_qp_q: KqpQ(kl) in HFBTHO
      !     Convert a "below-pwi" index (kl in HFBTHO) to a q.p. index
      ! pwi_uv_q: KpwiQ(kl) in HFBTHO
      !     For a given "below-pwi" index, how many U/V matrix elements have
      !     already been encountered? E.g. in calculation of particle density,
      !     V(ima) = V(pwi_uv_q(k) + ia)
      !     V(imb) = V(pwi_uv_q(k) + ib)
      ! pwi_active_q:
      !     Is a given q.p. index below the PWI?
      if (version >= 9) then
         allocate(pwi_start_n(nb), pwi_start_p(nb), pwi_dim_n(nb), pwi_dim_p(nb))
         allocate(pwi_qp_p(dqp), pwi_uv_p(dqp), pwi_qp_n(dqp), pwi_uv_n(dqp))

         pwi_start_n = 0;  pwi_dim_n = 0;  pwi_qp_n = 0;  pwi_uv_n = 0
         pwi_start_p = 0;  pwi_dim_p = 0;  pwi_qp_p = 0;  pwi_uv_p = 0

         ! From HFBTHO
         read(f) pwi_start_n(1:nb/2), pwi_start_p(1:nb/2)
         read(f) pwi_dim_n(1:nb/2),   pwi_dim_p(1:nb/2)
         read(f) pwi_qp_p(1:dqp/2)
         read(f) pwi_uv_p(1:dqp/2)
         read(f) pwi_qp_n(1:dqp/2)
         read(f) pwi_uv_n(1:dqp/2)

         ! NB! change from HFBTHO to make this match with the vectors imstart
         ! and isstart --- pwi_start_q contains the first q.p. in this block,
         ! not the last q.p. in the previous block
         pwi_start_n(1:nb/2) = pwi_start_n(1:nb/2) + 1
         pwi_start_p(1:nb/2) = pwi_start_p(1:nb/2) + 1

         ! Time-reversed blocks ---
         ! The counting of q.p.s continues as usual for the time-reversed blocks
         ! but the vectors which are functions of the number of q.p. levels
         ! below the pwi are more complicated. They aren't multiples of dqp, but
         ! rather have some length ncut <= dqp/2.
         ncut = sum(pwi_dim_n(1:nb/2))
         pwi_start_n(nb/2+1:nb)  = pwi_start_n(1:nb/2) + ncut
         pwi_dim_n(nb/2+1:nb)    = pwi_dim_n(1:nb/2)
         pwi_qp_n(ncut+1:2*ncut) = pwi_qp_n(1:ncut) + dqp/2
         pwi_uv_n(ncut+1:2*ncut) = pwi_uv_n(1:ncut) + dmat/2

         ncut = sum(pwi_dim_p(1:nb/2))
         pwi_start_p(nb/2+1:nb)  = pwi_start_p(1:nb/2) + ncut
         pwi_dim_p(nb/2+1:nb)    = pwi_dim_p(1:nb/2)
         pwi_qp_p(ncut+1:2*ncut) = pwi_qp_p(1:ncut) + dqp/2
         pwi_uv_p(ncut+1:2*ncut) = pwi_uv_p(1:ncut) + dmat/2

         ! Active arrays ---
         ! Which q.p. states are "active" (below the PWI)
         allocate(pwi_active_p(dqp), pwi_active_n(dqp))
         pwi_active_n = .false.
         pwi_active_p = .false.

         do ik=1, sum(pwi_dim_n)
            if (pwi_qp_n(ik) == 0) stop
            pwi_active_n(pwi_qp_n(ik)) = .true.
         end do
         do ik=1, sum(pwi_dim_p)
            if (pwi_qp_p(ik) == 0) stop
            pwi_active_p(pwi_qp_p(ik)) = .true.
         end do

         ! NB! Erase columns of U and V which weren't part of the minimization
         ! (that is, rho /= rho of U,V when e_k > pwi). This should give us the
         ! cutoff everywhere for free, thanks to Mika's matrix multiplication
         ! module.
         iqp = 0;  im = 0
         do iblock=1, nb
            do ic=1, db(iblock)
               iqp = iqp + 1

               ! Neutrons
               kactive = .false.
               do ik=pwi_start_n(iblock),pwi_start_n(iblock)+pwi_dim_n(iblock)-1
                  if (pwi_qp_n(ik) == iqp) then
                     kactive = .true.
                     exit
                  end if
               end do

               if (.not.kactive) then
                  un%elem(im+1:im+db(iblock)) = 0
                  vn%elem(im+1:im+db(iblock)) = 0
               else
                  if (im /= pwi_uv_n(ik)) stop 'im /= pwi_uv_n(ik).'
               end if

               ! Protons
               kactive = .false.
               do ik=pwi_start_p(iblock),pwi_start_p(iblock)+pwi_dim_p(iblock)-1
                  if (pwi_qp_p(ik) == iqp) then
                     kactive = .true.
                     exit
                  end if
               end do

               if (.not.kactive) then
                  up%elem(im+1:im+db(iblock)) = 0
                  vp%elem(im+1:im+db(iblock)) = 0
               else
                  if (im /= pwi_uv_p(ik)) stop 'im /= pwi_uv_n(ik).'
               end if

               im = im + db(iblock)
            end do
         end do
      end if

      ! U matrices are diagonal
      Un%ir2c = [ (i, i=1, nb) ]; Un%ic2r = Un%ir2c
      Un%ir2m(1) = 1
      do i=2,nb
        Un%ir2m(i) = Un%ir2m(i-1) + db(i-1)**2
      end do
      Un%ic2m = Un%ir2m
      call copy_block_structure(Un, Up)
      ! V matrices not quite
      Vn%ir2c(1:nb/2) = [ (i+nb/2, i=1, nb/2) ]
      Vn%ir2c(nb/2+1:nb) = [ (i, i=1, nb/2) ]
      Vn%ic2r = Vn%ir2c
      Vn%ir2m = Un%ir2m
      Vn%ic2m = Vn%ir2m(Vn%ic2r(:))
      call copy_block_structure(Vn, Vp)

      allocate(nr(dqp), nz(dqp), nl(dqp), ns(dqp), npar(dqp))
      read(f) nr(1:dqp/2)   ; nr(dqp/2+1:) = nr(1:dqp/2)
      read(f) nz(1:dqp/2)   ; nz(dqp/2+1:) = nz(1:dqp/2)
      read(f) nl(1:dqp/2)   ; nl(dqp/2+1:) = -nl(1:dqp/2)
      read(f) ns(1:dqp/2)   ; ns(dqp/2+1:) = -ns(1:dqp/2)
      read(f) npar(1:dqp/2) ; npar(dqp/2+1:) = npar(1:dqp/2)

      !-------------------------------------------------------------------------
      ! Wave functions
      !-------------------------------------------------------------------------
      read(f) nghl
      allocate(wdcori(nghl), y(nghl), z(nghl))
      allocate(wf(nghl,dqp), wfdr(nghl,dqp), wfdz(nghl,dqp), wfd2(nghl,dqp))
      read(f) wdcori(:)
      read(f) y(:)
      read(f) z(:)
      read(f) wf(:,1:dqp/2)    ! the transpose of the qhla_opt in HFBTHO
      read(f) wfdr(:,1:dqp/2)  ! the transpose of the fi1r_opt in HFBTHO
      read(f) wfdz(:,1:dqp/2)  ! the transpose of the fi1z_opt in HFBTHO
      read(f) wfd2(:,1:dqp/2)  ! the transpose of the fi2d_opt in HFBTHO
      do i=1,dqp/2
         if (ns(i) < 0) then
            wf(:,dqp/2+i)   = -wf(:,i)
            wfdr(:,dqp/2+i) = -wfdr(:,i)
            wfdz(:,dqp/2+i) = -wfdz(:,i)
            wfd2(:,dqp/2+i) = -wfd2(:,i)
         else
            wf(:,dqp/2+i)   = wf(:,i)
            wfdr(:,dqp/2+i) = wfdr(:,i)
            wfdz(:,dqp/2+i) = wfdz(:,i)
            wfd2(:,dqp/2+i) = wfd2(:,i)
         end if
      end do

      !-------------------------------------------------------------------------
      ! HFB EDF coupling constants
      !-------------------------------------------------------------------------
      read(f) hfb_cr0, hfb_crr, hfb_cdrho, hfb_ctau, hfb_ctj, hfb_crdj

      ! Beginning in ver. 8, logicals are stored as integers using:
      ! 1 => .TRUE., 0 => .FALSE.  Technically this means I should be able to
      ! simply read them as logical, but the following is safer.
      if (version >= 8) then
         itmp = -1
         read(f) itmp
         hfb_use_j2terms = HFBTHO_map_int_to_logical(itmp)
      else
         read(f) hfb_use_j2terms
      end if

      !-------------------------------------------------------------------------
      ! HFB finite-temperature quantities
      ! Introduced in v. 7 but completely wrong until v. 9
      !-------------------------------------------------------------------------
      ! Documentation ---
      ! Despite appearances and initialization in HFBTHO, the ft_fq's are not
      ! indexed by q.p. indices. Rather, they are indexed by levels below the
      ! pwi (kl in HFBTHO, pwi_start_q & pwi_dim_q in pnFAM). See e.g. DENSIT()
      ! in HFBTHO.
      ! ---
      ft_active = .false.
      if (version >= 9) then
         ! Beginning in ver. 8, logicals are stored as integers using:
         ! 1 => .TRUE., 0 => .FALSE.  Technically this means I should be able to
         ! simply read them as logical, but the following is safer.
         itmp = -1;  read(f) itmp
         ft_active = HFBTHO_map_int_to_logical(itmp)

         if (ft_active) then
            allocate(ft_fp(dqp), ft_fn(dqp))
            ft_fp(:) = 0;  ft_fn(:) = 0

            read(f) ft_temp
            read(f) ft_entropy(:)
            read(f) ft_fp(1:dqp/2)
            read(f) ft_fn(1:dqp/2)

            ncut = sum(pwi_dim_p(1:nb/2)); ft_fp(ncut+1:2*ncut) = ft_fp(1:ncut)
            ncut = sum(pwi_dim_n(1:nb/2)); ft_fn(ncut+1:2*ncut) = ft_fn(1:ncut)
         end if
      end if

      ! Initialize the auxiliary arrays imstart and isstart
      allocate(imstart(nb),isstart(nb))
      imstart(1) = 1 ; isstart(1) = 1
      do i=1,nb-1
         isstart(i+1) = sum(db(1:i))+1
         imstart(i+1) = sum(db(1:i)**2)+1
      end do

      close(f)

   end subroutine


   !----------------------------------------------------------------------------
   ! Compute beta decay quantities using an already-loaded HFB solution.
   ! Input:  Ep(:), En(:)
   ! Output: Q, E_gs, K_gs, P_gs (Q-value & g.s. energy, K-projection, parity)
   !----------------------------------------------------------------------------
   subroutine get_HFBTHO_betadecay_properties(Ep, En, Vp, Vn, Q, E_gs, K_gs, P_gs, lpr)
      use logger
      use constants, only : dmassnH
      implicit none

      type(blockmatrix), intent(in)  :: Vp, Vn
      real(dp),          intent(in)  :: Ep(:), En(:)
      logical, optional, intent(in)  :: lpr
      real(dp),          intent(out) :: Q, E_gs
      integer,           intent(out) :: K_gs, P_gs

      real(dp), parameter :: occ_tolerance = 5.0d-09

      integer  :: iqp, iprot, ineut
      real(dp) :: eqrpa_max, Op(size(Ep)), On(size(En))
      logical  :: mask_p(size(Ep)), mask_n(size(En))

      ! PWI
      integer :: ik
      logical :: kactive

      character(len=160) :: st

      ! Check for an HFB solution
      if (abs(hfb_lambda(1)) < epsilon(1.0_dp) .and. abs(hfb_lambda(2)) < epsilon(1.0_dp)) then
         call writelog('ERROR in sub.get_HFBTHO_betadecay_properties: no HFB solution found.&
            & Both lambda values are zero.')
         stop
      end if

      ! Q and E_gs via J. Engel et al., Phys. Rev. C 60, 014302 (1999)
      ! N.b.: yes, the EQRPA_max changes in a magic nucleus where we use the
      ! correct s.p. energies (rather than the arbtrary q.p. energies).
      ! However:
      !    1. The lambdas cancel in Q so there is no change w/a redefinition
      !    2. EQRPA_max and E_gs change in the same way w/a change of lambda
      !       so the interval satys the same
      !    3. The beta-decay integrals (at least for beta-minus... haven't
      !       checked beta-plus) are unchanged with a change of lambdas.
      !    4. The QRPA phonon spectrum shifts in the same way as E_gs and
      !       EQRPA_max, so the relative location of peaks does not change.
      ! As a result, there's no need to make a correction. See the note
      ! "Beta-decay Energetics when Pairing Gaps Collapse".
      eqrpa_max = dmassnH + hfb_lambda(1) - hfb_lambda(2)

      !-------------------------------------------------------------------------
      ! Compute ground-state energy
      !-------------------------------------------------------------------------
      mask_p(:) = .true.
      mask_n(:) = .true.

      ! Compute q.p. occupations (norm of lower component, at least) to mask the q.p. energies
      call get_HFBTHO_occupations(Vp, Op)
      call get_HFBTHO_occupations(Vn, On)

      ! Exclude full proton levels, empty neutron levels, and above PWI
      do iqp=1, size(Op)
         ! Neutron PWI
         kactive = .true.
         do ik=1, sum(pwi_dim_n)
            if (pwi_qp_n(ik) == iqp) then
               kactive = .true.
               exit
            end if
         end do

         if ((.not.kactive) .or. (abs(On(iqp)) < occ_tolerance)) then
            mask_n(iqp) = .false.
         end if

         ! Proton PWI
         kactive = .true.
         do ik=1, sum(pwi_dim_p)
            if (pwi_qp_p(ik) == iqp) then
               kactive = .true.
               exit
            end if
         end do

         if ((.not.kactive) .or. (abs(1.0_dp-Op(iqp)) < occ_tolerance)) then
            mask_p(iqp) = .false.
         end if
      end do

      ! Indices of lowest-energy proton and neutron levels which satisfy
      ! occupation requirements
      iprot = minloc(Ep, 1, mask_p)
      ineut = minloc(En, 1, mask_n)

      ! Ground-state energy is the lowest HFB 2-QP energy
      E_gs = Ep(iprot) + En(ineut)

      ! Q-value is now simply EQRPA_max - E_gs
      Q = eqrpa_max - E_gs

      ! K and parity
      K_gs = nl(iprot) + nl(ineut) + (ns(iprot) + ns(ineut))/2
      P_gs = (-1)**(npar(iprot)+npar(ineut))

      ! If lpr is True, print HFB properties
      if (present(lpr)) then
         if (lpr) then
            write(st,'(a)') 'HFB PARAMETERS:';                                    call writelog(st)
            write(st,'(a)') repeat('-',30);                                       call writelog(st)
            write(st,'(a,3x,1i0)')    'Z (parent) ............ :', hfb_npr(2);    call writelog(st)
            write(st,'(a,3x,1i0)')    'N (parent) ............ :', hfb_npr(1);    call writelog(st)
            write(st,'(a,3x,1i0)')    'A ..................... :', hfb_npr(3);    call writelog(st)
            write(st,'(a,1x,1f10.6)') 'Lambda prot. (ALA) (MeV):', hfb_lambda(2); call writelog(st)
            write(st,'(a,1x,1f10.6)') 'Lambda neut. (ALA) (MeV):', hfb_lambda(1); call writelog(st)
            write(st,'(a,1x,1f10.6)') 'Delta prot. ...... (MeV):', hfb_delta(2);  call writelog(st)
            write(st,'(a,1x,1f10.6)') 'Delta neut. ...... (MeV):', hfb_delta(1);  call writelog(st)
            write(st,'(a,1x,1f10.6)') 'MIN(Ep) .......... (MeV):', Ep(iprot);     call writelog(st)
            write(st,'(a,1x,1f10.6)') 'Occupation (vp**2) .... :', Op(iprot);     call writelog(st)
            write(st,'(a,1x,1f10.6)') 'MIN(En) .......... (MeV):', En(ineut);     call writelog(st)
            write(st,'(a,1x,1f10.6)') 'Occupation (vn**2) .... :', On(ineut);     call writelog(st)
            write(st,'(a,1x,1f10.6,2x,"[",1i0,a,"]")') 'Ground state ..... (MeV):', E_gs, &
               K_gs, merge('+', '-', P_gs > 1); call writelog(st)
            write(st,'(a,1x,1f10.6)') 'Q-value .......... (MeV):', Q;         call writelog(st)
            write(st,'(a,1x,1f10.6)') 'Max EQRPA ........ (MeV):', eqrpa_max; call writelog(st)
         end if
      end if

   end subroutine get_HFBTHO_betadecay_properties

   !----------------------------------------------------------------------------
   ! Compute q.p. occupations using blockmatrix V.
   ! This routine requires get_HFBTHO_solution to have been called previously
   ! to set the blockmatrix dimension variables db and nb.
   !----------------------------------------------------------------------------
   subroutine get_HFBTHO_occupations(v, occ)
      use blockmatrix_type, only : hfb_db => db, hfb_nb => nb

      implicit none

      type(blockmatrix), intent(in) :: v
      real(dp), intent(out) :: occ(sum(hfb_db(:)))

      integer :: ib, iqp, istart, iend, nme, nqp

      occ(:) = 0
      nme = 0
      nqp = 0

      do ib=1, hfb_nb
         ! Recall V are V_{p\pi}, i.e. q.p.s are column vectors.
         ! Norm (occupation/2) should be <column|column> = sum(col(:)**2).
         do iqp=1, hfb_db(ib)
            istart = nme + (iqp-1)*hfb_db(ib) + 1
            iend   = nme + iqp*hfb_db(ib)
            ! Because input is the blockmatrix V, there is no factor of 2 for
            ! the time-reversed partner. It is treated explicitly.
            occ(nqp+iqp) = sum(v%elem(istart:iend)**2)
         end do

         nme = nme + hfb_db(ib)**2
         nqp = nqp + hfb_db(ib)
      end do

   end subroutine get_HFBTHO_occupations


   !----------------------------------------------------------------------------
   ! Return the lowest-energy q.p. from the list *energies* while only
   ! considering q.p.'s below the PWI.  Optionally request only partially-
   ! occupied or partially-unoccupied levels.  Variable *it* determines whether
   ! to mask using the neutron (IT_NEUTRON) or proton (IT_PROTON) PWI arrays.
   !----------------------------------------------------------------------------
   function lowest_energy_qp(energies, v, it, mask_k_positive, mask_k_negative, &
      mask_2k, mask_parity, mask_occ, mask_unocc) result(iqp)

      use constants, only : IT_PROTON, IT_NEUTRON
      use logger

      implicit none

      real(dp), parameter :: occ_cut = 5d-9

      real(dp),          intent(in) :: energies(:)
      type(blockmatrix), intent(in) :: v
      integer,           intent(in) :: it
      integer, optional, intent(in) :: mask_2k, mask_parity
      logical, optional, intent(in) :: mask_k_positive, mask_k_negative, mask_occ, mask_unocc

      integer  :: iqp, i, n_degen
      integer, pointer :: pwi_start(:), pwi_dim(:), pwi_qp(:)
      real(dp) :: occ(size(v%elem))
      logical  :: active, mask(size(energies))
      character(len=200) :: s

      if (it /= IT_NEUTRON .and. it /= IT_PROTON) then
         write(s,'("ERROR in f.lowest_energy_qp: unrecognized isospin --- ",1i0)') it
         call writelog(s)
         stop
      else
         if (it == IT_PROTON) then
            pwi_start => pwi_start_p
            pwi_dim   => pwi_dim_p
            pwi_qp    => pwi_qp_p
         else if (it == IT_NEUTRON) then
            pwi_start => pwi_start_n
            pwi_dim   => pwi_dim_n
            pwi_qp    => pwi_qp_n
         end if
      end if

      call get_HFBTHO_occupations(v, occ)

      mask(:) = .true.
      do i=1, size(energies)
         ! PWI
         active = ((it == IT_NEUTRON .and. pwi_active_n(i)) .or. &
                   (it == IT_PROTON  .and. pwi_active_p(i)))
         if (.not.active) then
            mask(i) = .false.
         end if
         ! Occupation
         if (present(mask_occ)) then
            if (mask_occ .and. (abs(occ(i)) < occ_cut)) mask(i) = .false.
         end if
         if (present(mask_unocc)) then
            if (mask_unocc .and. (abs(1.0_dp-occ(i)) < occ_cut)) mask(i) = .false.
         end if
         ! K > 0 option
         if (present(mask_k_positive)) then
            if (mask_k_positive .and. (2*nl(i)+ns(i) < 0)) mask(i) = .false.
         end if
         ! K < 0 option
         if (present(mask_k_negative)) then
            if (mask_k_negative .and. (2*nl(i)+ns(i) > 0)) mask(i) = .false.
         end if
         ! K option
         if (present(mask_2k)) then
            if ((2*nl(i)+ns(i)) /= mask_2k) mask(i) = .false.
         end if
         ! Parity option
         if (present(mask_parity)) then
            if (npar(i) /= mask_parity) mask(i) = .false.
         end if
      end do

      iqp = minloc(energies, 1, mask)

      ! Degeneracy should be a warning
      n_degen = 0
      do i=1, size(energies)
         if (.not.mask(i)) cycle
         if (abs(energies(i)-energies(iqp)) < 1d-6) n_degen = n_degen+1
      end do
      if (n_degen > 1) then
         write(s,'("WARNING in f.lowest_energy_qp: found ",1i0," degenerate not-fully-occupied qp levels")') n_degen
         call writelog(s)
         call writelog('')
      end if

   end function lowest_energy_qp


   !---------------------------------------------------------------------------
   ! Compute the normal density of the mean field in coordinate space.
   ! Extended to finite temperature and the PWI cut from HFBTHO.
   ! it = 1 neutron density
   !    = 2 proton density
   !    = 3 total density
   !---------------------------------------------------------------------------
   subroutine hfb_density(it, vp, vn, up, un, rho)
      use constants, only : IT_NEUTRON, IT_PROTON, IT_ISOSCALAR
      use blockmatrix_type
      implicit none
      integer,  intent(in)  :: it
      real(dp), intent(in)  :: vp(:), vn(:), up(:), un(:)
      real(dp), intent(out) :: rho(nghl)

      real(dp) :: rho_ab, f
      integer :: ib, nd, ir, ic, i1, i2, k, iuv


      ! Check of isospin
      if (it /= IT_NEUTRON .and. it /= IT_PROTON .and. it /= IT_ISOSCALAR) then
         write(*,'(A)') 'ERROR in hfb_density(): Choose from neutron,&
            & proton, or total density.'
         stop
      end if

      rho(:) = 0

      do ib=1, nb
         nd = db(ib)

         do ic=1, nd
            do ir=1, nd
               ! sp w.f. indices
               i1 = isstart(ib)-1 + ir
               i2 = isstart(ib)-1 + ic

               !----------------------------------------------------------------
               ! Neutrons
               !----------------------------------------------------------------
               if (it == IT_NEUTRON .or. it == IT_ISOSCALAR) then
                  rho_ab = 0
                  do k=pwi_start_n(ib), pwi_start_n(ib)+pwi_dim_n(ib)-1
                     iuv = pwi_uv_n(k)

                     ! Finite temperature
                     if (ft_active) then
                        ! Access f variables through pwi indices, not qp indices!
                        f = ft_fn(k)
                     else
                        f = 0
                     end if

                     rho_ab = rho_ab + vn(iuv+ir)*vn(iuv+ic)*(1-f) + un(iuv+ir)*un(iuv+ic)*f
                  end do

                  if (ns(i1) == ns(i2) .and. nl(i1) == nl(i2) .and. npar(i1) == npar(i2)) then
                     rho(:) = rho(:) + rho_ab*wf(:,i1)*wf(:,i2)
                  end if
               end if
               !----------------------------------------------------------------
               ! Protons
               !----------------------------------------------------------------
               if (it == IT_PROTON .or. it == IT_ISOSCALAR) then
                  rho_ab = 0
                  do k=pwi_start_p(ib), pwi_start_p(ib)+pwi_dim_p(ib)-1
                     iuv = pwi_uv_p(k)

                     ! Finite temperature
                     if (ft_active) then
                        ! Access f variables through pwi indices, not qp indices!
                        f = ft_fp(k)
                     else
                        f = 0
                     end if

                     rho_ab = rho_ab + vp(iuv+ir)*vp(iuv+ic)*(1-f) + up(iuv+ir)*up(iuv+ic)*f
                  end do

                  if (ns(i1) == ns(i2) .and. nl(i1) == nl(i2) .and. npar(i1) == npar(i2)) then
                     rho(:) = rho(:) + rho_ab*wf(:,i1)*wf(:,i2)
                  end if
               end if
            end do
         end do
      end do

      rho(:) = rho(:)*wdcori(:)  ! remove integration weights
   end subroutine hfb_density


   !---------------------------------------------------------------------------
   ! Two-level toy model data (only used for debugging purposes, to have
   ! easily tractable mock data)
   !---------------------------------------------------------------------------
   subroutine get_toymodel(ep, vp, up, en, vn, un)
      implicit none
      real(dp), allocatable, intent(out) :: Ep(:), En(:)   ! quasiparticle energies for protons and neutrons
      type(blockmatrix) :: Vp, Up, Vn, Un ! U and V matrices

      integer :: i, dqp

      nb = 2
      allocate(db(2))
      db(:) = 1
      allocate(ep(2), en(2))
      call allocate_blockmatrix(up,2)
      call allocate_blockmatrix(un,2)
      call allocate_blockmatrix(vp,2)
      call allocate_blockmatrix(vn,2)
      ep = 2 ; en = 3
      vp%elem = [.5, -.5] ; up%elem = sqrt(1-vp%elem**2)
      vn%elem = [1, -1] ; un%elem = sqrt(1-vn%elem**2)
      up%ic2r = [1, 2] ; up%ir2c = [1,2] ; up%ir2m = [1,2] ; up%ic2m = [1,2]
      un%ic2r = [1, 2] ; un%ir2c = [1,2] ; un%ir2m = [1,2] ; un%ic2m = [1,2]
      vp%ic2r = [2, 1] ; vp%ir2c = [2,1] ; vp%ir2m = [2,1] ; vp%ic2m = [2,1]
      vn%ic2r = [2, 1] ; vn%ir2c = [2,1] ; vn%ir2m = [2,1] ; vn%ic2m = [2,1]

      allocate(nr(2), nz(2), nl(2), ns(2), npar(2))
      nr = 0 ; nz = 0; nl = 0; ns = [1,-1] ; npar = 1
      nghl = 1
      allocate(wdcori(nghl), y(nghl), wf(nghl,2), wfdr(nghl,2), wfdz(nghl,2), wfd2(nghl,2))
      wdcori = 1 ; y = 1 ; wf(1,1) = 1 ; wf(1,2) = 1 ; wfdr = 0 ; wfdz = 0 ; wfd2 = 0

      dqp = 2
      do i=1,dqp/2
         if (ns(i) < 0) then
            wf(:,dqp/2+i)   = -wf(:,i)
            wfdr(:,dqp/2+i) = -wfdr(:,i)
            wfdz(:,dqp/2+i) = -wfdz(:,i)
            wfd2(:,dqp/2+i) = -wfd2(:,i)
         else
            wf(:,dqp/2+i)   = wf(:,i)
            wfdr(:,dqp/2+i) = wfdr(:,i)
            wfdz(:,dqp/2+i) = wfdz(:,i)
            wfd2(:,dqp/2+i) = wfd2(:,i)
         end if
      end do

      allocate(isstart(2), imstart(2))
      isstart = [1,2]
      imstart = [1,2]

      allocate(pwi_start_n(nb), pwi_start_p(nb), pwi_dim_n(nb), pwi_dim_p(nb))
      allocate(pwi_qp_p(dqp), pwi_uv_p(dqp), pwi_qp_n(dqp), pwi_uv_n(dqp))

      pwi_start_n = 0;  pwi_dim_n = 0;  pwi_qp_n = 0;  pwi_uv_n = 0
      pwi_start_p = 0;  pwi_dim_p = 0;  pwi_qp_p = 0;  pwi_uv_p = 0

      pwi_start_n = [1,2]
      pwi_start_p = [1,2]
      pwi_dim_n   = [1,1]
      pwi_dim_p   = [1,1]
      pwi_qp_p    = [1,2]
      pwi_uv_p    = [1,2]
      pwi_qp_n    = [1,2]
      pwi_uv_n    = [1,2]

   end subroutine

   !----------------------------------------------------------------------------
   ! In hfbtho_storage.f90, logicals are stored using MERGE() via
   ! MERGE(1, 0, <logical>) as a way to allow compatibility between ifort and
   ! gfortran.  Here I simply back-convert them and warn if I see something
   ! unexpected. Default behavior if I encounter a negative integer is to
   ! take its absolute value and warn, thus converting iforts .TRUE. = -1 to
   ! whatever the local machine's default is.
   !----------------------------------------------------------------------------
   function HFBTHO_map_int_to_logical(hfb_int) result(hfb_bool)
      implicit none
      integer, intent(in) :: hfb_int
      logical :: hfb_bool

      ! The result is true iff |hfb_int| > 0
      hfb_bool = .false.
      if (abs(hfb_int) > 0) hfb_bool = .true.

      if (hfb_int < 0) then
         write(*,'(A)') 'Warning in HFBTHO_map_int_to_logical:&
            & integer value < 0 mapped to .TRUE.'
      end if
   end function HFBTHO_map_int_to_logical

end module hfbtho_basis
