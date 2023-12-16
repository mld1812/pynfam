!------------------------------------------------------------------------------
!> This modules bridges hfbtho and pnfam variables. It's main responsiblity is
!> renaming variables, and expanding the basis to include explicit time-reversed
!> states. It also implements a pairing cutoff but setting all QP quantities
!> above the cutoff to zero. This prevents them from contributing to the FAM result.
!>
!> @authors M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
!> @authors E.M. Ney, UNC Chapel Hill, 2018-
!> @authors T. Li, Michigan State Univ, 2020-
!------------------------------------------------------------------------------

module hfb_solution
   !use hfbtho_utilities, only : ipr, dp => pr

   !---------------------------------------------------------------------------
   ! Variables reused from hfbtho (some renamed)
   !---------------------------------------------------------------------------
   ! Coordinate space quantities
   use hfbtho, only : bz !> Harmonic oscillator length in z direction
   use hfbtho, only : bp !> Harmonic oscillator legnth in r direction
   use hfbtho, only : ngh !> Number of 1D coord. space quadrature points (Hermite=z-direction)
   use hfbtho, only : ngl !> Number of 1D coord. space quadrature points (Laguerre=r-direction)
   use hfbtho, only : nghl!> Total 2D coord. space quadrature points (ngh*ngl)
   use hfbtho, only : y=>y_opt !> 1/r (dimensionful) on the 2D quadrature grid
   use hfbtho, only : z=>fh !> z (dimensionful) on the 2D quadrature grid
   use hfbtho, only : xl !> Dimensionless r on 1D quadrature grid
   use hfbtho, only : xh !> Dimensionless z on 1D quadrature grid
   use hfbtho, only : wl !> Gauss-Laguerre weights on 1D r quadrature grid
   use hfbtho, only : wh !> Gauss-Hermite weights on 1D z quadrature grid
   use hfbtho, only : wdcori !> Inverse quadrature weights times Jacobian (for getting bare wavefunctions from wf)

   ! Configuration Space Quantities
   use hfbtho, only : n00 !> Number of HO shells in the basis
   use hfbtho, only : nzx !> Maximal nz quantum number
   use hfbtho, only : nrx !> Maximal nr quantum number
   use hfbtho, only : nlx !> Maximal nl quantum number
   use hfbtho, only : nox !> Maximal K (2Omega = 2K -1)
   use hfbtho, only : ntx !> Total number of basis states with K>0
   use hfbtho, only : nttx !> Total number of states with (k,r,l,s) < ntx = number with (k,r,l,s,pi,z)
   use hfbtho, only : nrlx !> Maximal quantum number between 2nr and nl
   use hfbtho, only : ndx !> (Maximal) Total number of basis states per block
   use hfbtho, only : nbx !> Number of HFB blocks (pnfams nb/2)
   use hfbtho, only : hnr => nr !> Number nodes in r direction (r=distance from intrinsic z axis)
   use hfbtho, only : hnz => nz!> Number nodes in z direction
   use hfbtho, only : hnl => nl!> Orbital angular momentum projection on intrinsic z axis
   use hfbtho, only : hns => ns!> Spin projection on intrinsic z axis *2 (so that it's an integer)
   use hfbtho, only : hnpar => npar !> Pairty (1 for even, 2 for odd)

   ! General HFB results
   use hfbtho, only : hfb_npr => npr !> Number of particles (N, Z, A)
   use hfbtho, only : hfb_lambda => ala !> N and P Fermi energy lagrange multiplier (arbitrary if pairing collapses)
   use hfbtho, only : hfb_elast => alast !> N and P Last bound SP energy (not used, though differs from lambda if pairing collapses)
   use hfbtho, only : hfb_delta => del !> N and P pairing gaps
   use hfbtho, only : hfb_beta2 => bet !> Total deformation beta2
   use hfbtho, only : hfb_energy => ehfb !> HFB solution binding energy
   use hfbtho, only : hfb_pairing_window => pwi !> Pairing energy cutoff (NB: HFBTHO smears this cutoff with a Fermi-fct)
   use hfbtho, only : rmat => rk !> Density matrix in HO basis
   use hfbtho, only : kmat => ak !> Pairing density matrix in HO basis
   use hfbtho, only : TAU !> kinetic density
   use hfbtho, only : d2rho => DRO !> Laplacian of density
   use hfbtho, only : drrho => NABLAR !> r derivative of density
   use hfbtho, only : dzrho => NABLAZ !> z derivative of density

   ! Functional (declared in UNEDF but used by HFBTHO)
   use hfbtho, only : hfb_cpair => CpV0 !> HFB T=1 pairing coupling
   use hfbtho, only : hfb_alpha_pair => CpV1 !> HFB T=1 pairing coupling density dependence exponent
   use hfbtho, only : hfb_use_j2terms => use_j2terms !> j^2 terms active in the functional?
   use hfbtho, only : hbzero !> \hbar^2/(2m_nucleon)
   use hfbtho, only : rho_nm !> nuclear matter density

   ! Statistical FAM (finite-temp or EFA)
   use hfbtho, only : ft_active => switch_on_temperature !> Is T>0 and active in HFB calculation
   use hfbtho, only : ft_temp => temper !> Temperature in MeV
   use hfbtho, only : ft_entropy => entropy ! (N,P,Total) entropy of solution
   use hfbtho, only : keyblo !> Index within the blocking candidates of blocked QP (0 if no blocking)
   use hfbtho, only : bloblo, blo123, blok1k2, blo123d, blok1k2d, tb, numax

   implicit none

   integer, parameter, private :: dp = kind(1.0d0)
   integer, parameter, private :: ipr = kind(1)

   !---------------------------------------------------------------------------
   ! New variables (or altered hfbtho variables, mostly doubling due to explicit TR states)
   !---------------------------------------------------------------------------
   ! Coordinate space quantities
   real(dp), allocatable, target :: wf(:,:) !> Dimensionless basis wave functions stored on 2D dimensionless quadrature grid
                                            !> with sqrt(quadrature weights) mixed in for fast integration
   real(dp), allocatable, target :: wfdr(:,:) !> Partial derivative of wf in r direction [1/L]
   real(dp), allocatable, target :: wfdz(:,:) !> Partial derivative of wf in z direction [1/L]
   real(dp), allocatable, target :: wfd2(:,:) !> Laplacian of wf (only r and z derivatives) [1/L^2]
   real(dp), allocatable, target :: wfdp(:,:) !> Derivative of phi, without imaginary unit "i" at the beginning
   real(dp), allocatable, target :: wfd2_all(:,:) !> Full Laplacian of wavefunction, with phi-derivative added

   ! Configuration space quantities
   integer(ipr), allocatable :: nr(:) !> Number nodes in r direction (r=distance from intrinsic z axis)
   integer(ipr), allocatable :: nz(:) !> Number nodes in z direction
   integer(ipr), allocatable :: nl(:) !> Orbital angular momentum projection on intrinsic z axis
   integer(ipr), allocatable :: ns(:) !> Spin projection on intrinsic z axis *2 (so that it's an integer)
   integer(ipr), allocatable :: npar(:) !> Pairty (1 for even, 2 for odd)

   ! Block matrix structure
   integer(ipr) :: hdqp !> dimension of the basis
   integer(ipr) :: hdmat !> number of non-trivial matrix elements in U, V, etc.
   integer(ipr) :: hnb !> the number of blocks (declared in blockmatrix_type)
   integer(ipr), allocatable :: hdb(:) !> the number of states (declared in blockmatrix_type)
   integer(ipr), allocatable :: ib_itx(:) !> the number of states (declared in blockmatrix_type)

   integer(ipr), allocatable, target :: pwi_start_n(:), pwi_start_p(:) !> Index of 1st active QP in block i, 
      !> not counting inactive states. In HFBTHO: ka(ib,it)+1.
   integer(ipr), allocatable, target :: pwi_dim_n(:), pwi_dim_p(:) !> Number of active states in block i
      !> Typical iteration is "do k=pwi_start_q(ib),pwi_start_q(ib)+pwi_dim_q(ib)-1". In HFBTHO: kd(ib,it).
   integer(ipr), allocatable, target :: pwi_qp_n(:), pwi_qp_p(:) !> Convert index in list of active QPs 
      !> to index among all QPs. In HFBTHO: KqpQ(kl).
   integer(ipr), allocatable, target :: pwi_uv_n(:), pwi_uv_p(:) !> Vectorized matrix index for active QP index i
      !> E.g. V(ima) = V(pwi_uv_q(k)+ia); V(imb) = V(pwi_uv_q(k) + ib). In HFBTHO: KqpQ(kl).
   logical(ipr), allocatable, target :: pwi_active_n(:), pwi_active_p(:) !> Boolean indicating if QP is below pwi.
   integer(ipr), allocatable :: num_spin_up(:) !> Number of basis states with spin up in each block
   integer, allocatable :: new_order(:)

   ! HFB Result and Functional
   real(dp) :: hfb_cr0 = 0, hfb_crr = 0, hfb_cdrho = 0 !> Just isovector T=1 couplings
   real(dp) :: hfb_ctau = 0, hfb_ctj = 0, hfb_crdj = 0 !> Just isovector T=1 couplings
   real(dp), allocatable :: hfb_occ_n(:), hfb_occ_p(:)
   real(dp), allocatable :: Ep(:), En(:) !> Quasiparticle energies
   real(dp), allocatable :: Vp(:), Up(:), Vn(:), Un(:) !> U and V matrices

   ! Statistical FAM (finite-temp or EFA)
   real(dp), allocatable, target :: ft_fp(:), ft_fn(:) !> FT QP occupations 1/(1+exp(E/T))

contains

   !---------------------------------------------------------------------------
   ! Load the HFBTHO solution
   !---------------------------------------------------------------------------
   subroutine get_raw_hfb_solution(ierr)
      use hfbtho_interface, only : hfbtho_program
      use pnfam_constants, only : IT_NEUTRON, IT_PROTON

      use hfbtho, only : fn_T, fp_T
      use hfbtho, only : REqpN, RVqpN, RUqpN
      use hfbtho, only : REqpP, RVqpP, RUqpP
      use hfbtho, only : uk
      use hfbtho, only : ka, kd, KqpN, KqpP, KpwiN, KpwiP
      use hfbtho, only : qhla_opt, fi1r_opt, fi1z_opt, fi2d_opt
      use hfbtho, only : nb, id
      use hfbtho, only : crho, cdrho, ctau, crdr, cj, crdj

      implicit none

      integer(ipr), intent(out) :: ierr
      integer(ipr) :: i, iqp, ib, ic, im, ik, it, ncut, ipt
      integer(ipr) :: blo_qp, blo_blo, blo_level, blo_blo_chk, blo_level_chk
      real(dp) :: blo_e, blo_e_chk, fT
      ierr = 0

      ! Reconstruct hfbtho solution
      !-------------------------------------------------------------------------
      call hfbtho_program(ierr)
      if(ierr /= 0) return

      ! Double the matrix sizes for explicit TR states
      hnb = nb*2
      hdqp = sum(id)*2
      hdmat = sum(id*id)*2

      ! Double block structure for explicit TR states
      !-------------------------------------------------------------------------
      if(allocated(hdb)) deallocate(hdb)
      allocate(hdb(hnb))
      hdb(1:hnb/2) = id(:); hdb(hnb/2+1:) = id(:)
      ! NB: DO NOT DEALLOCATE id, it is used in hfbtho_solver.base and pnfam_spatial_mtxels.prep_gaussian

      ! Given a state index, return the block
      If (allocated(ib_itx)) Deallocate(ib_itx)
      allocate(ib_itx(hdqp))
      ipt = 0
      do ib=1,hnb
         do i=1,hdb(ib)
            ipt = ipt + 1
            ib_itx(ipt) = ib
         end do
      end do

      ! Double quantum numbers and wavefunctions
      !-------------------------------------------------------------------------
      If(Allocated(nr)) Deallocate(nr,nz,nl,ns,npar)
      allocate(nr(hdqp), nz(hdqp), nl(hdqp), ns(hdqp), npar(hdqp))
      If(Allocated(wf)) Deallocate(wf,wfdr,wfdz,wfd2)
      allocate(wf(nghl,hdqp), wfdr(nghl,hdqp), wfdz(nghl,hdqp), wfd2(nghl,hdqp))
      nr(1:hdqp/2) = hnr(:)     ; nr(hdqp/2+1:) = nr(1:hdqp/2)
      nz(1:hdqp/2) = hnz(:)     ; nz(hdqp/2+1:) = nz(1:hdqp/2)
      nl(1:hdqp/2) = hnl(:)     ; nl(hdqp/2+1:) = -nl(1:hdqp/2)
      ns(1:hdqp/2) = hns(:)     ; ns(hdqp/2+1:) = -ns(1:hdqp/2)
      npar(1:hdqp/2) = hnpar(:) ; npar(hdqp/2+1:) = npar(1:hdqp/2)
      ! NB: DO NOT DEALLOCATE quantum numbers, they are used in hfbtho_solver.base and pnfam_spatial_mtxels.prep_gaussian

      wf(:,1:hdqp/2) = transpose(qhla_opt(:,:))
      wfdr(:,1:hdqp/2) = transpose(fi1r_opt(:,:))
      wfdz(:,1:hdqp/2) = transpose(fi1z_opt(:,:))
      wfd2(:,1:hdqp/2) = transpose(fi2d_opt(:,:))
      do i=1,hdqp/2
         if (ns(i) < 0) then
            wf(:,hdqp/2+i)   = -wf(:,i)
            wfdr(:,hdqp/2+i) = -wfdr(:,i)
            wfdz(:,hdqp/2+i) = -wfdz(:,i)
            wfd2(:,hdqp/2+i) = -wfd2(:,i)
         else
            wf(:,hdqp/2+i)   = wf(:,i)
            wfdr(:,hdqp/2+i) = wfdr(:,i)
            wfdz(:,hdqp/2+i) = wfdz(:,i)
            wfd2(:,hdqp/2+i) = wfd2(:,i)
         end if
      end do
      deallocate(qhla_opt, fi1r_opt, fi1z_opt, fi2d_opt)
      allocate(qhla_opt(1,1), fi1r_opt(1,1), fi1z_opt(1,1), fi2d_opt(1,1))



      ! Double quasiparticles
      !-------------------------------------------------------------------------
      If(Allocated(Vn)) Deallocate(Ep,Vp,Up, En,Vn,Un)
      allocate(Ep(hdqp), Vp(hdmat), Up(hdmat), En(hdqp), Vn(hdmat), Un(hdmat))
      Ep(1:hdqp/2)   = REqpP(:) ; Ep(hdqp/2+1:)  =  Ep(1:hdqp/2)
      Vp(hdmat/2+1:) = RVqpP(:) ; Vp(1:hdmat/2)  = -Vp(hdmat/2+1:)
      Up(1:hdmat/2)  = RUqpP(:) ; Up(hdmat/2+1:) =  Up(1:hdmat/2)
      En(1:hdqp/2)   = REqpN(:) ; En(hdqp/2+1:)  =  En(1:hdqp/2)
      Vn(hdmat/2+1:) = RVqpN(:) ; Vn(1:hdmat/2)  = -Vn(hdmat/2+1:)
      Un(1:hdmat/2)  = RUqpN(:) ; Un(hdmat/2+1:) =  Un(1:hdmat/2)
      deallocate(REqpP, RVqpP, RUqpP, REqpN, RVqpN, RUqpN)
      allocate(REqpP(1), RVqpP(1), RUqpP(1), REqpN(1), RVqpN(1), RUqpN(1))

      ! Occupations
      If(Allocated(hfb_occ_n)) Deallocate(hfb_occ_n, hfb_occ_p)
      allocate(hfb_occ_n(hdqp), hfb_occ_p(hdqp))
      hfb_occ_n(1:hdqp/2) = uk(:,IT_NEUTRON) ; hfb_occ_n(hdqp/2+1:) = hfb_occ_n(1:hdqp/2)
      hfb_occ_p(1:hdqp/2) = uk(:,IT_PROTON)  ; hfb_occ_p(hdqp/2+1:) = hfb_occ_p(1:hdqp/2)
      deallocate(uk)
      allocate(uk(1,1))

      ! Pairing Window
      !-------------------------------------------------------------------------
      If(Allocated(pwi_start_n)) Deallocate(pwi_start_n, pwi_start_p,&
          pwi_dim_n, pwi_dim_p, pwi_qp_p, pwi_uv_p, pwi_qp_n, pwi_uv_n)
      allocate(pwi_start_n(hnb), pwi_start_p(hnb), pwi_dim_n(hnb), pwi_dim_p(hnb),&
         pwi_qp_p(hdqp), pwi_uv_p(hdqp), pwi_qp_n(hdqp), pwi_uv_n(hdqp))

      pwi_start_n = 0;  pwi_dim_n = 0;  pwi_qp_n = 0;  pwi_uv_n = 0
      pwi_start_p = 0;  pwi_dim_p = 0;  pwi_qp_p = 0;  pwi_uv_p = 0

      pwi_start_n(1:hnb/2) = ka(:,IT_NEUTRON)
      pwi_start_p(1:hnb/2) = ka(:,IT_PROTON)
      pwi_dim_n(1:hnb/2) = kd(:,IT_NEUTRON)
      pwi_dim_p(1:hnb/2) = kd(:,IT_PROTON)
      pwi_qp_p(1:hdqp/2) = KqpP(:)
      pwi_uv_p(1:hdqp/2) = KpwiP(:)
      pwi_qp_n(1:hdqp/2) = KqpN(:)
      pwi_uv_n(1:hdqp/2) = KpwiN(:)
      deallocate(ka,kd,KqpP,KqpN,KpwiP,KpwiN)
      allocate(ka(1,1),kd(1,1),KqpP(1),KqpN(1),KpwiP(1),KpwiN(1))

      ! NB! change from HFBTHO to make this match with the vectors imstart
      ! and isstart --- pwi_start_q contains the first q.p. in this block,
      ! not the last q.p. in the previous block
      pwi_start_n(1:hnb/2) = pwi_start_n(1:hnb/2) + 1
      pwi_start_p(1:hnb/2) = pwi_start_p(1:hnb/2) + 1

      ! Time-reversed blocks ---
      ! The counting of q.p.s continues as usual for the time-reversed blocks
      ! but the vectors which are functions of the number of q.p. levels
      ! below the pwi are more complicated. They aren't multiples of hdqp, but
      ! rather have some length ncut <= hdqp/2.
      ncut = sum(pwi_dim_n(1:hnb/2))
      pwi_start_n(hnb/2+1:hnb)  = pwi_start_n(1:hnb/2) + ncut
      pwi_dim_n(hnb/2+1:hnb)    = pwi_dim_n(1:hnb/2)
      pwi_qp_n(ncut+1:2*ncut) = pwi_qp_n(1:ncut) + hdqp/2
      pwi_uv_n(ncut+1:2*ncut) = pwi_uv_n(1:ncut) + hdmat/2

      ncut = sum(pwi_dim_p(1:hnb/2))
      pwi_start_p(hnb/2+1:hnb)  = pwi_start_p(1:hnb/2) + ncut
      pwi_dim_p(hnb/2+1:hnb)    = pwi_dim_p(1:hnb/2)
      pwi_qp_p(ncut+1:2*ncut) = pwi_qp_p(1:ncut) + hdqp/2
      pwi_uv_p(ncut+1:2*ncut) = pwi_uv_p(1:ncut) + hdmat/2

      ! Active arrays ---
      ! Which q.p. states are "active" (below the PWI)
      If(Allocated(pwi_active_p)) Deallocate(pwi_active_p, pwi_active_n)
      allocate(pwi_active_p(hdqp), pwi_active_n(hdqp))
      pwi_active_n = .false.
      pwi_active_p = .false.
      do ik=1, sum(pwi_dim_n)
         if (pwi_qp_n(ik) == 0) then
            ierr=2; return
         end if
         pwi_active_n(pwi_qp_n(ik)) = .true.
      end do
      do ik=1, sum(pwi_dim_p)
         if (pwi_qp_p(ik) == 0) then
            ierr=2; return
         end if
         pwi_active_p(pwi_qp_p(ik)) = .true.
      end do

#ifdef USE_HBLAS
#if USE_HBLAS==1
      ! Set up basis for BLAS accelerated hamiltonian calculation
      !-------------------------------------------------------------------------
      ! Sort: put states with ns==1 at the beginning and those with ns==-1 at the end in each block
      ! The i-th element of array "new_order" gives the index of wf that should be put at position i
      If(Allocated(num_spin_up)) Deallocate(num_spin_up)
      allocate(num_spin_up(hnb))
      If(Allocated(new_order)) Deallocate(new_order)
      allocate(new_order(hdqp))
      new_order = 0; num_spin_up = 0
      im = 0
      do ib=1, hnb
         iqp = 0
         do ic=im+1, im+hdb(ib)
            if (ns(ic)>0) then
               iqp = iqp + 1
               new_order(im+iqp) = ic
            end if
         end do
         num_spin_up(ib) = iqp
         do ic=im+1, im+hdb(ib)
            if (ns(ic)<0) then
               iqp = iqp + 1
               new_order(im+iqp) = ic
            end if
         end do
         im = im + hdb(ib)
      end do
      nr(:) = nr(new_order); nz(:) = nz(new_order); nl(:) = nl(new_order)
      ns(:) = ns(new_order); npar(:) = npar(new_order)
      wf(:,:) = wf(:,new_order)
      wfdr(:,:) = wfdr(:,new_order)
      wfdz(:,:) = wfdz(:,new_order)
      wfd2(:,:) = wfd2(:,new_order)
      call reorder_hfbmatrix_row(Un, new_order)
      call reorder_hfbmatrix_row(Vn, new_order)
      call reorder_hfbmatrix_row(Up, new_order)
      call reorder_hfbmatrix_row(Vp, new_order)

      ! Pre-compute some quantities for module hamiltonian(blas)
      If(Allocated(wfdp)) Deallocate(wfdp, wfd2_all)
      allocate(wfdp(nghl,hdqp), wfd2_all(nghl,hdqp))
      do i=1, hdqp
         wfdp(:,i) = y(:)*nl(i)*wf(:,i)
         wfd2_all(:,i) = wfd2(:,i) - y(:)*nl(i)*wfdp(:,i)
      end do
#endif
#endif

      ! Functional - just isovector T=1 couplings
      !-------------------------------------------------------------------------
      hfb_cr0   = crho(1)
      hfb_crr   = cdrho(1)
      hfb_cdrho = crdr(1)
      hfb_ctau  = ctau(1)
      hfb_ctj   = cj(1)
      hfb_crdj  = crdj(1)

      ! Finite-Temperature
      !-------------------------------------------------------------------------
      ! EVAN: (10/27/2020)
      ! - In hfbdiag() ft_fq's are indexed exactly the same as Eqp, though they use
      !   1/2*(1-tanh(1/2*E/T)) to evaluate which gives different precision than using
      !   1/(1+exp(E/T)) at about 10^-16.
      ! - In ALambda() ft_fq's are OVERWRITTEN with below pwi indices
      ! - In DENSIT() the ft_fq's are accessed with below pwi indicies to construct
      !   the density using only below pwi occupations.
      ! To avoid confusion, I will recreate occupations here using same indexing as
      ! Eqp (all QPs not just below pwi).
      !-------------------------------------------------------------------------
      if (ft_active) then
         If(Allocated(ft_fp)) Deallocate(ft_fp, ft_fn)
         allocate(ft_fp(hdqp), ft_fn(hdqp))
         ft_fp(:) = 0;  ft_fn(:) = 0
         ! Reconstruct occupations to have same indexing as Eqp
         ! Note Eqp are already double for explicit TR states here
         do i = 1, hdqp
            fT = 0.5_dp*(1.0_dp-Tanh(0.5_dp*En(i)/ft_temp))
            ft_fn(i) = fT
            fT = 0.5_dp*(1.0_dp-Tanh(0.5_dp*Ep(i)/ft_temp))
            ft_fp(i) = fT
         end do
         deallocate(fp_T, fn_T)
         allocate(fp_T(1), fn_T(1))
      end if

      ! Set columns in U and V for QPs above PWI to zero
      Call erase_above_pwi(ierr)

   end subroutine

   !----------------------------------------------------------------------------
   ! Erase columns of U and V which weren't part of the minimization
   ! (that is, rho /= rho of U,V when e_k > pwi). This should give us the
   ! cutoff everywhere for free with block matrix multiplcation.
   !----------------------------------------------------------------------------
   subroutine erase_above_pwi(ierr)
      implicit none

      integer, intent(out) :: ierr
      integer(ipr) :: iqp, im, ic, ik, iblock
      logical :: kactive
      ierr = 0

      iqp = 0;  im = 0
      do iblock=1, hnb
         do ic=1, hdb(iblock)
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
               En(iqp) = 0
               if(ft_active) ft_fn(iqp) = 0
               Un(im+1:im+hdb(iblock)) = 0
               Vn(im+1:im+hdb(iblock)) = 0
            else
               if (im /= pwi_uv_n(ik)) then
                  ierr = 3; return !'im /= pwi_uv_n(ik).'
               end if
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
               Ep(iqp) = 0
               if(ft_active) ft_fp(iqp) = 0
               Up(im+1:im+hdb(iblock)) = 0
               Vp(im+1:im+hdb(iblock)) = 0
            else
               if (im /= pwi_uv_p(ik)) then
                  ierr=3; return !'im /= pwi_uv_n(ik).'
               end if
            end if

            im = im + hdb(iblock)
         end do
      end do

   end subroutine

   !----------------------------------------------------------------------------
   ! Subroutine to reorder basis states of HFBTHO matrices  according to ns.
   ! Assumes block diagonal.
   ! Input: U(block matrix), new_order
   ! Output: U(block matrix)
   !----------------------------------------------------------------------------
   subroutine reorder_hfbmatrix_row(U, new_order)
      implicit none

      real(dp), intent(inout) :: U(:)
      integer, intent(in) :: new_order(:)
      integer :: iblock, ir, ic, im

      ir = 0; im = 0
      do iblock=1, hnb
         do ic=1, hdb(iblock)
            U(im+1:im+hdb(iblock)) = U(im-ir+new_order(ir+1:ir+hdb(iblock)))
            im = im + hdb(iblock)
         end do
         ir = ir + hdb(iblock)
      end do
   end subroutine reorder_hfbmatrix_row

   !---------------------------------------------------------------------------
   ! Return the normal density in COORDINATE space. N(it=1), P(it=2), or N+P(it=3).
   !---------------------------------------------------------------------------
   subroutine hfb_density_coord(it, rho)
      use hfbtho, only : ro
      use pnfam_constants, only : IT_NEUTRON, IT_PROTON, IT_ISOSCALAR
      implicit none

      integer, intent(in) :: it
      real(dp), intent(out) :: rho(nghl)

      if (it /= IT_NEUTRON .and. it /= IT_PROTON .and. it /= IT_ISOSCALAR) then
         write(*,'(A)') 'ERROR in hfb_density_coord(): Choose from neutron,&
            & proton, or total density.'
         stop
      end if

      rho = 0
      if (it==IT_NEUTRON .or. it==IT_ISOSCALAR) then
         rho = rho + ro(:,IT_NEUTRON)
      end if
      if (it==IT_PROTON .or. it==IT_ISOSCALAR) then
         rho = rho + ro(:,IT_PROTON)
      end if

   end subroutine

   !----------------------------------------------------------------------------
   ! Convert a K>0 q.p. to K<0 or vice versa.
   !----------------------------------------------------------------------------
   function qp_t_reverse(iqp) result(iqp_r)
      implicit none
      
      integer, intent(in) :: iqp
      integer :: nqp, iqp_r
      
      nqp   = sum(hdb)
      iqp_r = 0
      
      if (iqp <= nqp/2) then
         iqp_r = iqp + nqp/2
      else
         iqp_r = iqp - nqp/2
      end if
      
   end function qp_t_reverse

   !----------------------------------------------------------------------------
   ! For a given IQP, find which block and state in the block it corresponds
   ! to for comparing with HFBTHO.  Ported from 'display_qp_level' in
   ! pnfam_matrix/mod_common.
   !----------------------------------------------------------------------------
   subroutine iqp_to_block_and_state(iqp, iblock, istate)
      implicit none
      integer, intent(in)  :: iqp
      integer, intent(out) :: iblock, istate

      istate = 0
      do iblock=1, hnb
         if ((istate+hdb(iblock)) > iqp) then
            exit
         else
            istate = istate+hdb(iblock)
         end if
      end do
      istate = iqp-sum(hdb(1:iblock-1))
      
   end subroutine iqp_to_block_and_state

   !---------------------------------------------------------------------------
   ! Two-level toy model data (only used for debugging purposes, to have
   ! easily tractable mock data)
   !---------------------------------------------------------------------------
   subroutine get_raw_toymodel
      implicit none

      integer :: i, hdqpt

      ! Matrix sizes
      hnb = 2
      if(allocated(hdb)) deallocate(hdb)
      allocate(hdb(hnb))
      hdb(:) = 1
      hdqp = sum(hdb)
      hdmat = sum(hdb*hdb)

      ! Quasiparticle solution
      If(Allocated(Vn)) Deallocate(Ep,Vp,Up, En,Vn,Un)
      allocate(Ep(hdqp), Vp(hdmat), Up(hdmat), En(hdqp), Vn(hdmat), Un(hdmat))
      ep = 2 ; en = 3
      vp = [.5, -.5] ; up = sqrt(1-vp**2)
      vn = [1, -1]   ; un = sqrt(1-vn**2)
      !up%ic2r = [1, 2] ; up%ir2c = [1,2] ; up%ir2m = [1,2] ; up%ic2m = [1,2]
      !un%ic2r = [1, 2] ; un%ir2c = [1,2] ; un%ir2m = [1,2] ; un%ic2m = [1,2]
      !vp%ic2r = [2, 1] ; vp%ir2c = [2,1] ; vp%ir2m = [2,1] ; vp%ic2m = [2,1]
      !vn%ic2r = [2, 1] ; vn%ir2c = [2,1] ; vn%ir2m = [2,1] ; vn%ic2m = [2,1]

      ! Quantum numbers
      If(Allocated(nr)) Deallocate(nr,nz,nl,ns,npar)
      allocate(nr(hdqp), nz(hdqp), nl(hdqp), ns(hdqp), npar(hdqp))
      nr = 0 ; nz = 0; nl = 0; ns = [1,-1] ; npar = 1

      ! Custom quadrature grid
      nghl = 1
      if(Allocated(wdcori)) Deallocate(wdcori, y)
      allocate(wdcori(nghl), y(nghl))
      wdcori = 1 ; y = 1

      ! Wavefunctions
      If(Allocated(wf)) Deallocate(wf,wfdr,wfdz,wfd2)
      allocate(wf(nghl,hdqp), wfdr(nghl,hdqp), wfdz(nghl,hdqp), wfd2(nghl,hdqp))
      wf(1,1) = 1 ; wf(1,2) = 1 ; wfdr = 0 ; wfdz = 0 ; wfd2 = 0

      ! Explicit TR states?
      do i=1,hdqp/2
         if (ns(i) < 0) then
            wf(:,hdqp/2+i)   = -wf(:,i)
            wfdr(:,hdqp/2+i) = -wfdr(:,i)
            wfdz(:,hdqp/2+i) = -wfdz(:,i)
            wfd2(:,hdqp/2+i) = -wfd2(:,i)
         else
            wf(:,hdqp/2+i)   = wf(:,i)
            wfdr(:,hdqp/2+i) = wfdr(:,i)
            wfdz(:,hdqp/2+i) = wfdz(:,i)
            wfd2(:,hdqp/2+i) = wfd2(:,i)
         end if
      end do

      If(Allocated(pwi_start_n)) Deallocate(pwi_start_n, pwi_start_p,&
          pwi_dim_n, pwi_dim_p, pwi_qp_p, pwi_uv_p, pwi_qp_n, pwi_uv_n)
      allocate(pwi_start_n(hnb), pwi_start_p(hnb), pwi_dim_n(hnb), pwi_dim_p(hnb),&
         pwi_qp_p(hdqp), pwi_uv_p(hdqp), pwi_qp_n(hdqp), pwi_uv_n(hdqp))

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

   subroutine get_wf_2nd_derivs(i, wfd2z, wfd2r, wfdrz)

      implicit none
      integer, intent(in) :: i !> HO wf index
      real(dp), intent(out) :: wfd2z(nghl) !> 2nd z derivative of wf [1/L^2], d/dz(d/dz)
      real(dp), intent(out) :: wfd2r(nghl) !> 2nd r derivative of wf [1/L^2], d/dr(d/dr)
      real(dp), intent(out) :: wfdrz(nghl) !> 2nd r derivative of wf [1/L^2], d/dr(d/dz)

      ! Compute some 2nd derivatives for 2BC
      wfd2z(:) = ((z(:)/bz)**2 - (2*nz(i) + 1))/(bz*bz)*wf(:,i)
      wfd2r(:) = wfd2(:,i) - wfd2z(:) - y(:)*wfdr(:,i)
      wfdrz(:) = wfdz(:,i)*wfdr(:,i)/wf(:,i)

   end subroutine

end module hfb_solution
