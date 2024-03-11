!------------------------------------------------------------------------------
!> This module implements the interface defined in pnfam_setup, which takes
!> perturbed densities in single particle configuration space and computes the
!> perturbed hamiltonian in single particle configuration space, in the
!> proton-neutron sector.
!>
!> It first converts single particle densities to coordinate-space, then simultaneously
!> computes the perturbed fields in coordinate space (as functional derivatives of the
!> energy density functional) and converts back to single particle configuration space.
!>
!> * This version of the module optimizes run-time by:
!>   Calculating mean fields and using BLAS operations. See README_hblas.md for details.
!>
!> @authors M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
!> @authors E.M. Ney, UNC Chapel Hill, 2018-
!> @authors T. Li, Michigan State Univ, 2020-
!------------------------------------------------------------------------------
module pnfam_hamiltonian
   use pnfam_constants, only : pi
   use pnfam_logger
   use hfb_solution
   use type_blockmatrix
   use type_densityset
   implicit none
   private

   integer, parameter :: dp = kind(1.0d0)
   integer            :: maxval_db = -1 ! Maximum value in db

   ! Coordinate space local densities, custom type for organization
   type(density_set) :: rhox

   type field
      real(dp), allocatable, dimension(:) :: re, im
   end type field

   type array_2d
      real(dp), allocatable, dimension(:,:) :: re, im
   end type array_2d

   public :: calc_hamiltonian

contains

   !---------------------------------------------------------------------------
   ! Wrapper routine to calculate the perturbed hamiltonian given perturbed densities,
   ! both in s.p. configuration space, in two steps:
   ! 1. Transform drho and dkappa from s.p. config space to coordinate space
   ! 2. Compute coordinate space dh and dDelta (functional derivatives of EDF)
   !    and transform back to s.p. config space
   !---------------------------------------------------------------------------
   subroutine calc_hamiltonian(rerho_pn, imrho_pn, rekp, imkp,&
                               rerho_np, imrho_np, rekm, imkm,&
                               reh_pn, imh_pn, redp, imdp,&
                               reh_np, imh_np, redm, imdm)
      implicit none
      ! Perturbed generalized density (sp basis)
      type(blockmatrix), intent(inout) :: rerho_pn, imrho_pn, rekp, imkp,&
                                          rerho_np, imrho_np, rekm, imkm
      ! Perturbed hamiltonian (sp basis)
      type(blockmatrix), intent(inout) :: reh_pn, imh_pn, redp, imdp,&
                                          reh_np, imh_np, redm, imdm
      if (maxval_db<0) maxval_db = Maxval(db)
      if (.not.rhox%alloc) call allocate_density(rhox, nghl)

      ! PN and + Fields
      call density(rerho_pn, imrho_pn, rekp, imkp, rhox)
      call meanfield(rhox, reh_pn, imh_pn) ! h_pn
      call pairingfield(rhox, redp, imdp)  ! Delta+ from rho-tilde+(r)

      ! NP and - Fields
      call density(rerho_np, imrho_np, rekm, imkm, rhox)
      call meanfield(rhox, reh_np, imh_np) ! h_np
      call pairingfield(rhox, redm, imdm)  ! Delta- from rho-tilde-(r)
   end subroutine calc_hamiltonian


   !---------------------------------------------------------------------------
   ! Some helper routines for custom types
   !---------------------------------------------------------------------------
   subroutine allocate_field(a, n)
      type(field) :: a
      integer, intent(in) :: n
      allocate(a%re(n), a%im(n))
   end subroutine allocate_field

   subroutine deallocate_field(a)
      type(field) :: a
      deallocate(a%re, a%im)
   end subroutine deallocate_field

   subroutine add_field(a, symbol, add)
      ! a(:) = a(:) + symbol*add(:)
      type(field), intent(inout) :: a
      integer, intent(in) :: symbol
      type(field), intent(in) :: add
      a%re(1:nghl) = a%re(1:nghl) + symbol*add%re(1:nghl)
      a%im(1:nghl) = a%im(1:nghl) + symbol*add%im(1:nghl)
   end subroutine add_field

   subroutine add_mul_1d(a, b, beta, c)
      ! c(:) = a(:)*b(:) + beta*c(:)
      real(dp), dimension(:), intent(in), contiguous :: a, b
      real(dp), intent(in) :: beta
      real(dp), dimension(:), intent(inout), contiguous :: c
      c(1:nghl) = a(1:nghl)*b(1:nghl) + beta*c(1:nghl)
   end subroutine add_mul_1d

   subroutine allocate_2darray(a, nr, nc)
      type(array_2d) :: a
      integer, intent(in) :: nr, nc
      allocate(a%re(nr, nc), a%im(nr, nc))
   end subroutine allocate_2darray

   subroutine deallocate_2darray(a)
      type(array_2d) :: a
      deallocate(a%re, a%im)
   end subroutine deallocate_2darray


   !---------------------------------------------------------------------------
   ! Computation of the local densities in coordinate space
   !---------------------------------------------------------------------------
   subroutine density(rerho, imrho, rek, imk, ds)

      implicit none

      type(blockmatrix), intent(in) :: rerho, imrho, rek, imk
      type(density_set), intent(inout) :: ds

      integer  :: ia, ib, ix, iy, iz
      real(dp), pointer, dimension(:) :: wf_b, dr_wf_b, dz_wf_b, dp_wf_b

      ! Result of Sum_a wf_a (or its derivatives) * rho_ab
      ! index gives the spin (0 is not used)
      type(array_2d), dimension(-1:1), save :: wfa_rhoab, dr_wfa_rhoab, dp_wfa_rhoab, dz_wfa_rhoab
      type(field), save :: wf_wf_prod, aux

      ! The following are needed to allow OpenMP reduction, which does not currently
      ! work for derived types.  Hacky :(
      real(dp), dimension(:), allocatable, save :: dsrerho, dsimrho, dsretau, dsimtau,   &
         dsretjrr, dsimtjrr, dsretjpr, dsimtjpr, dsretjzr, dsimtjzr, dsretjrp, dsimtjrp, &
         dsretjpp, dsimtjpp, dsretjzp, dsimtjzp, dsretjrz, dsimtjrz, dsretjpz, dsimtjpz, &
         dsretjzz, dsimtjzz, dsresr, dsimsr, dsresp, dsimsp, dsresz, dsimsz, dsretr,     &
         dsimtr, dsretp, dsimtp, dsretz, dsimtz, dsrejr, dsimjr, dsrejp, dsimjp,         &
         dsrejz, dsimjz, dsrefr, dsimfr, dsrefp, dsimfp, dsrefz, dsimfz, dsregs,         &
         dsimgs, dsrerb, dsimrb, dsresbr, dsimsbr, dsresbp, dsimsbp, dsresbz, dsimsbz
      
      ! Allocate auxiliary arrays
      ! Get around a problem with OpenMP on OS X by making these allocatable
      ! http://stackoverflow.com/a/13884056/656740
      if (.not.allocated(dsrerho)) then
         allocate(dsrerho(nghl), dsimrho(nghl), dsretau(nghl), dsimtau(nghl), dsretjrr(nghl),   &
            dsimtjrr(nghl), dsretjpr(nghl), dsimtjpr(nghl), dsretjzr(nghl), dsimtjzr(nghl),     &
            dsretjrp(nghl), dsimtjrp(nghl), dsretjpp(nghl), dsimtjpp(nghl), dsretjzp(nghl),     &
            dsimtjzp(nghl), dsretjrz(nghl), dsimtjrz(nghl), dsretjpz(nghl), dsimtjpz(nghl),     &
            dsretjzz(nghl), dsimtjzz(nghl), dsresr(nghl), dsimsr(nghl), dsresp(nghl),           &
            dsimsp(nghl), dsresz(nghl), dsimsz(nghl), dsretr(nghl), dsimtr(nghl), dsretp(nghl), &
            dsimtp(nghl), dsretz(nghl), dsimtz(nghl), dsrejr(nghl), dsimjr(nghl), dsrejp(nghl), &
            dsimjp(nghl), dsrejz(nghl), dsimjz(nghl), dsrefr(nghl), dsimfr(nghl), dsrefp(nghl), &
            dsimfp(nghl), dsrefz(nghl), dsimfz(nghl), dsregs(nghl), dsimgs(nghl), dsrerb(nghl), &
            dsimrb(nghl), dsresbr(nghl), dsimsbr(nghl), dsresbp(nghl), dsimsbp(nghl),           &
            dsresbz(nghl), dsimsbz(nghl))
         call allocate_field(wf_wf_prod,nghl)
         call allocate_field(aux,nghl)
         do iz = -1, 1, 2
            call allocate_2darray(wfa_rhoab(iz),nghl,maxval_db)
            call allocate_2darray(dr_wfa_rhoab(iz),nghl,maxval_db)
            call allocate_2darray(dp_wfa_rhoab(iz),nghl,maxval_db)
            call allocate_2darray(dz_wfa_rhoab(iz),nghl,maxval_db)
         end do
      end if

      ! Initialize the densities to zero; zeroing wf_wf_prod, aux is unnecessary.
      ! Zeroing wfa_rhoab, dr_wfa_rhoab, dp_wfa_rhoab and dz_wfa_rhoab is also unnecessary,
      ! which is ensured by conditions marked by (*) below in this subroutine.
      dsrerho  = 0 ; dsimrho  = 0 ; dsretau  = 0 ; dsimtau  = 0 ; dsretjrr = 0 ; dsimtjrr = 0
      dsretjpr = 0 ; dsimtjpr = 0 ; dsretjzr = 0 ; dsimtjzr = 0 ; dsretjrp = 0 ; dsimtjrp = 0
      dsretjpp = 0 ; dsimtjpp = 0
      dsretjzp = 0 ; dsimtjzp = 0 ; dsretjrz = 0 ; dsimtjrz = 0 ; dsretjpz = 0 ; dsimtjpz = 0
      dsretjzz = 0 ; dsimtjzz = 0
      dsresr   = 0 ; dsimsr   = 0 ; dsresp   = 0 ; dsimsp   = 0 ; dsresz   = 0 ; dsimsz   = 0
      dsretr   = 0 ; dsimtr   = 0 ; dsretp   = 0 ; dsimtp   = 0
      dsretz   = 0 ; dsimtz   = 0 ; dsrejr   = 0 ; dsimjr   = 0 ; dsrejp   = 0 ; dsimjp   = 0
      dsrejz   = 0 ; dsimjz   = 0 ; dsrefr   = 0 ; dsimfr   = 0
      dsrefp   = 0 ; dsimfp   = 0 ; dsrefz   = 0 ; dsimfz   = 0 ; dsregs   = 0 ; dsimgs   = 0
      dsrerb   = 0 ; dsimrb   = 0 ; dsresbr  = 0 ; dsimsbr  = 0 ; dsresbp  = 0 ; dsimsbp  = 0
      dsresbz  = 0 ; dsimsbz  = 0

      ! Loop over blocks
      do ix = 1, nb
         iy = rerho%ir2c(ix) ! column index of current block
         if (iy==0) cycle

         ! Calculate Sum_a wf_a (or its derivatives) * rho_ab
         ! Use OpenMP-enabled BLAS library to achieve best performance
         ia = isstart(ix)

         ! wf_a has spin up
         if (num_spin_up(ix)>0) then
            call dgemm('N','N',nghl,db(iy),num_spin_up(ix),1.0_dp,wf(1,ia),nghl,&
                       rerho%elem(rerho%ir2m(ix)),db(ix),0.0_dp,wfa_rhoab(1)%re(1,1),nghl)
            call dgemm('N','N',nghl,db(iy),num_spin_up(ix),1.0_dp,wf(1,ia),nghl,&
                       imrho%elem(imrho%ir2m(ix)),db(ix),0.0_dp,wfa_rhoab(1)%im(1,1),nghl)
            call dgemm('N','N',nghl,db(iy),num_spin_up(ix),1.0_dp,wfdr(1,ia),nghl,&
                       rerho%elem(rerho%ir2m(ix)),db(ix),0.0_dp,dr_wfa_rhoab(1)%re(1,1),nghl)
            call dgemm('N','N',nghl,db(iy),num_spin_up(ix),1.0_dp,wfdr(1,ia),nghl,&
                       imrho%elem(imrho%ir2m(ix)),db(ix),0.0_dp,dr_wfa_rhoab(1)%im(1,1),nghl)
            call dgemm('N','N',nghl,db(iy),num_spin_up(ix),1.0_dp,wfdp(1,ia),nghl,&
                       rerho%elem(rerho%ir2m(ix)),db(ix),0.0_dp,dp_wfa_rhoab(1)%re(1,1),nghl)
            call dgemm('N','N',nghl,db(iy),num_spin_up(ix),1.0_dp,wfdp(1,ia),nghl,&
                       imrho%elem(imrho%ir2m(ix)),db(ix),0.0_dp,dp_wfa_rhoab(1)%im(1,1),nghl)
            call dgemm('N','N',nghl,db(iy),num_spin_up(ix),1.0_dp,wfdz(1,ia),nghl,&
                       rerho%elem(rerho%ir2m(ix)),db(ix),0.0_dp,dz_wfa_rhoab(1)%re(1,1),nghl)
            call dgemm('N','N',nghl,db(iy),num_spin_up(ix),1.0_dp,wfdz(1,ia),nghl,&
                       imrho%elem(imrho%ir2m(ix)),db(ix),0.0_dp,dz_wfa_rhoab(1)%im(1,1),nghl)
         end if

         ! wf_a has spin down
         if (num_spin_dw(ix)>0) then
            ia = ia + num_spin_up(ix)
            call dgemm('N','N',nghl,db(iy),num_spin_dw(ix),1.0_dp,wf(1,ia),nghl,&
                       rerho%elem(rerho%ir2m(ix)+num_spin_up(ix)),db(ix),0.0_dp,wfa_rhoab(-1)%re(1,1),nghl)
            call dgemm('N','N',nghl,db(iy),num_spin_dw(ix),1.0_dp,wf(1,ia),nghl,&
                       imrho%elem(imrho%ir2m(ix)+num_spin_up(ix)),db(ix),0.0_dp,wfa_rhoab(-1)%im(1,1),nghl)
            call dgemm('N','N',nghl,db(iy),num_spin_dw(ix),1.0_dp,wfdr(1,ia),nghl,&
                       rerho%elem(rerho%ir2m(ix)+num_spin_up(ix)),db(ix),0.0_dp,dr_wfa_rhoab(-1)%re(1,1),nghl)
            call dgemm('N','N',nghl,db(iy),num_spin_dw(ix),1.0_dp,wfdr(1,ia),nghl,&
                       imrho%elem(imrho%ir2m(ix)+num_spin_up(ix)),db(ix),0.0_dp,dr_wfa_rhoab(-1)%im(1,1),nghl)
            call dgemm('N','N',nghl,db(iy),num_spin_dw(ix),1.0_dp,wfdp(1,ia),nghl,&
                       rerho%elem(rerho%ir2m(ix)+num_spin_up(ix)),db(ix),0.0_dp,dp_wfa_rhoab(-1)%re(1,1),nghl)
            call dgemm('N','N',nghl,db(iy),num_spin_dw(ix),1.0_dp,wfdp(1,ia),nghl,&
                       imrho%elem(imrho%ir2m(ix)+num_spin_up(ix)),db(ix),0.0_dp,dp_wfa_rhoab(-1)%im(1,1),nghl)
            call dgemm('N','N',nghl,db(iy),num_spin_dw(ix),1.0_dp,wfdz(1,ia),nghl,&
                       rerho%elem(rerho%ir2m(ix)+num_spin_up(ix)),db(ix),0.0_dp,dz_wfa_rhoab(-1)%re(1,1),nghl)
            call dgemm('N','N',nghl,db(iy),num_spin_dw(ix),1.0_dp,wfdz(1,ia),nghl,&
                       imrho%elem(imrho%ir2m(ix)+num_spin_up(ix)),db(ix),0.0_dp,dz_wfa_rhoab(-1)%im(1,1),nghl)
         end if

         ! Calculate local densities in coordinate space
         !!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(iy,isstart,db,ns,y,wf,wfdr,wfdp,wfdz) &
         !!$OMP& SHARED(wfa_rhoab,dr_wfa_rhoab,dp_wfa_rhoab,dz_wfa_rhoab) &
         !!$OMP& REDUCTION(+:dsrerho) REDUCTION(+:dsimrho) REDUCTION(+:dsretau) REDUCTION(+:dsimtau) &
         !!$OMP& REDUCTION(+:dsretjrr) REDUCTION(+:dsimtjrr) REDUCTION(+:dsretjpr) REDUCTION(+:dsimtjpr) &
         !!$OMP& REDUCTION(+:dsretjzr) REDUCTION(+:dsimtjzr) REDUCTION(+:dsretjrp) REDUCTION(+:dsimtjrp) &
         !!$OMP& REDUCTION(+:dsretjpp) REDUCTION(+:dsimtjpp) REDUCTION(+:dsretjzp) REDUCTION(+:dsimtjzp) &
         !!$OMP& REDUCTION(+:dsretjrz) REDUCTION(+:dsimtjrz) REDUCTION(+:dsretjpz) REDUCTION(+:dsimtjpz) &
         !!$OMP& REDUCTION(+:dsretjzz) REDUCTION(+:dsimtjzz) REDUCTION(+:dsresr) REDUCTION(+:dsimsr) &
         !!$OMP& REDUCTION(+:dsresp) REDUCTION(+:dsimsp) REDUCTION(+:dsresz) REDUCTION(+:dsimsz) &
         !!$OMP& REDUCTION(+:dsretr) REDUCTION(+:dsimtr) REDUCTION(+:dsretp) REDUCTION(+:dsimtp) &
         !!$OMP& REDUCTION(+:dsretz) REDUCTION(+:dsimtz) REDUCTION(+:dsrejr) REDUCTION(+:dsimjr) &
         !!$OMP& REDUCTION(+:dsrejp) REDUCTION(+:dsimjp) REDUCTION(+:dsrejz) REDUCTION(+:dsimjz) &
         !!$OMP& REDUCTION(+:dsrefr) REDUCTION(+:dsimfr) REDUCTION(+:dsrefp) REDUCTION(+:dsimfp) &
         !!$OMP& REDUCTION(+:dsrefz) REDUCTION(+:dsimfz) REDUCTION(+:dsregs) REDUCTION(+:dsimgs)
         do iz = 1, db(iy)

            ! Column vector |b>
            ib = isstart(iy)+iz-1
            wf_b => wf(:, ib)
            dr_wf_b => wfdr(:, ib)
            dp_wf_b => wfdp(:, ib)
            dz_wf_b => wfdz(:, ib)

            ! Diagonal in spin (|a> = |+>,  |b> = |+>; or |a> = |->,  |b> = |->)
            ! Condition below ensures no calc when no |a> with the required spin exists
            if (((ns(ib) == +1) .and. (num_spin_up(ix)>0)) .or. &
                ((ns(ib) == -1) .and. (num_spin_dw(ix)>0))) then ! (*)

               ! Auxiliary array
               wf_wf_prod%re(:) = wfa_rhoab(ns(ib))%re(:,iz)*wf_b(:)
               wf_wf_prod%im(:) = wfa_rhoab(ns(ib))%im(:,iz)*wf_b(:)

               ! rho: scalar density
               dsrerho(:) = dsrerho(:) + wf_wf_prod%re(:)
               dsimrho(:) = dsimrho(:) + wf_wf_prod%im(:)

               ! s_z: z-component of vector density
               dsresz(:) = dsresz(:) + ns(ib)*wf_wf_prod%re(:)
               dsimsz(:) = dsimsz(:) + ns(ib)*wf_wf_prod%im(:)

               ! j_phi: phi-component of current density
               aux%re(:) = (wfa_rhoab(ns(ib))%re(:,iz)*dp_wf_b(:)+dp_wfa_rhoab(ns(ib))%re(:,iz)*wf_b(:))/2
               aux%im(:) = (wfa_rhoab(ns(ib))%im(:,iz)*dp_wf_b(:)+dp_wfa_rhoab(ns(ib))%im(:,iz)*wf_b(:))/2
               dsrejp(:) = dsrejp(:) + aux%re(:)
               dsimjp(:) = dsimjp(:) + aux%im(:)

               ! J_{phi,z}: (phi,z)-component of tensor density [= J_{2,3}/rho]
               !aux%re(:) = (wfa_rhoab(ns(ib))%re(:,iz)*dp_wf_b(:)+dp_wfa_rhoab(ns(ib))%re(:,iz)*wf_b(:))/2
               !aux%im(:) = (wfa_rhoab(ns(ib))%im(:,iz)*dp_wf_b(:)+dp_wfa_rhoab(ns(ib))%im(:,iz)*wf_b(:))/2
               dsretjpz(:) = dsretjpz(:) + ns(ib)*aux%re(:)
               dsimtjpz(:) = dsimtjpz(:) + ns(ib)*aux%im(:)

               ! tau: kinetic density
               aux%re(:) = dr_wfa_rhoab(ns(ib))%re(:,iz)*dr_wf_b(:) + dp_wfa_rhoab(ns(ib))%re(:,iz)*dp_wf_b(:) &
                        + dz_wfa_rhoab(ns(ib))%re(:,iz)*dz_wf_b(:)
               aux%im(:) = dr_wfa_rhoab(ns(ib))%im(:,iz)*dr_wf_b(:) + dp_wfa_rhoab(ns(ib))%im(:,iz)*dp_wf_b(:) &
                        + dz_wfa_rhoab(ns(ib))%im(:,iz)*dz_wf_b(:)
               dsretau(:) = dsretau(:) + aux%re(:)
               dsimtau(:) = dsimtau(:) + aux%im(:)

               ! T_z: z-component of spin-kinetic density
               ! aux%re(:) = dr_wfa_rhoab(ns(ib))%re(:,iz)*dr_wf_b(:) + dp_wfa_rhoab(ns(ib))%re(:,iz)*dp_wf_b(:) &
               !           + dz_wfa_rhoab(ns(ib))%re(:,iz)*dz_wf_b(:)
               ! aux%im(:) = dr_wfa_rhoab(ns(ib))%im(:,iz)*dr_wf_b(:) + dp_wfa_rhoab(ns(ib))%im(:,iz)*dp_wf_b(:) &
               !           + dz_wfa_rhoab(ns(ib))%im(:,iz)*dz_wf_b(:)
               dsretz(:) = dsretz(:) + ns(ib)*aux%re(:)
               dsimtz(:) = dsimtz(:) + ns(ib)*aux%im(:)

               ! j_rho: rho-component of current density
               aux%re(:) = (wfa_rhoab(ns(ib))%re(:,iz)*dr_wf_b(:) - dr_wfa_rhoab(ns(ib))%re(:,iz)*wf_b(:))/2
               aux%im(:) = (wfa_rhoab(ns(ib))%im(:,iz)*dr_wf_b(:) - dr_wfa_rhoab(ns(ib))%im(:,iz)*wf_b(:))/2
               dsrejr(:) = dsrejr(:) - aux%im(:)
               dsimjr(:) = dsimjr(:) + aux%re(:)

               ! J_{rho,z}: (rho,z)-component of tensor density
               ! aux%re(:) = (wfa_rhoab(ns(ib))%re(:,iz)*dr_wf_b(:) - dr_wfa_rhoab(ns(ib))%re(:,iz)*wf_b(:))/2
               ! aux%im(:) = (wfa_rhoab(ns(ib))%im(:,iz)*dr_wf_b(:) - dr_wfa_rhoab(ns(ib))%im(:,iz)*wf_b(:))/2
               dsretjrz(:) = dsretjrz(:) - ns(ib)*aux%im(:)
               dsimtjrz(:) = dsimtjrz(:) + ns(ib)*aux%re(:)

               ! j_z: z-component of current density
               aux%re(:) = (wfa_rhoab(ns(ib))%re(:,iz)*dz_wf_b(:) - dz_wfa_rhoab(ns(ib))%re(:,iz)*wf_b(:))/2
               aux%im(:) = (wfa_rhoab(ns(ib))%im(:,iz)*dz_wf_b(:) - dz_wfa_rhoab(ns(ib))%im(:,iz)*wf_b(:))/2
               dsrejz(:) = dsrejz(:) - aux%im(:)
               dsimjz(:) = dsimjz(:) + aux%re(:)

               ! J_{z,z}: (z,z)-component of tensor density
               ! aux%re(:) = (wfa_rhoab(ns(ib))%re(:,iz)*dz_wf_b(:) - dz_wfa_rhoab(ns(ib))%re(:,iz)*wf_b(:))/2
               ! aux%im(:) = (wfa_rhoab(ns(ib))%im(:,iz)*dz_wf_b(:) - dz_wfa_rhoab(ns(ib))%im(:,iz)*wf_b(:))/2
               dsretjzz(:) = dsretjzz(:) - ns(ib)*aux%im(:)
               dsimtjzz(:) = dsimtjzz(:) + ns(ib)*aux%re(:)

               ! Del.s: z-component of divergence of s
               aux%re(:) = wfa_rhoab(ns(ib))%re(:,iz)*dz_wf_b(:) + dz_wfa_rhoab(ns(ib))%re(:,iz)*wf_b(:)
               aux%im(:) = wfa_rhoab(ns(ib))%im(:,iz)*dz_wf_b(:) + dz_wfa_rhoab(ns(ib))%im(:,iz)*wf_b(:)
               dsregs(:) = dsregs(:) + ns(ib)*aux%re(:)
               dsimgs(:) = dsimgs(:) + ns(ib)*aux%im(:)

               ! F_rho: rho-component of tensor-kinetic density
               aux%re(:) = dr_wfa_rhoab(ns(ib))%re(:,iz)*dz_wf_b(:) + dz_wfa_rhoab(ns(ib))%re(:,iz)*dr_wf_b(:)
               aux%im(:) = dr_wfa_rhoab(ns(ib))%im(:,iz)*dz_wf_b(:) + dz_wfa_rhoab(ns(ib))%im(:,iz)*dr_wf_b(:)
               dsrefr(:) = dsrefr(:) + ns(ib)*aux%re(:)/2
               dsimfr(:) = dsimfr(:) + ns(ib)*aux%im(:)/2

               ! F_phi: phi-component of tensor-kinetic density
               aux%re(:) = dp_wfa_rhoab(ns(ib))%re(:,iz)*dz_wf_b(:) - dz_wfa_rhoab(ns(ib))%re(:,iz)*dp_wf_b(:)
               aux%im(:) = dp_wfa_rhoab(ns(ib))%im(:,iz)*dz_wf_b(:) - dz_wfa_rhoab(ns(ib))%im(:,iz)*dp_wf_b(:)
               dsrefp(:) = dsrefp(:) - ns(ib)*aux%im(:)/2
               dsimfp(:) = dsimfp(:) + ns(ib)*aux%re(:)/2

               ! F_z: z-component of tensor-kinetic density
               aux%re(:) = dz_wfa_rhoab(ns(ib))%re(:,iz)*dz_wf_b(:)
               aux%im(:) = dz_wfa_rhoab(ns(ib))%im(:,iz)*dz_wf_b(:)
               dsrefz(:) = dsrefz(:) + ns(ib)*aux%re(:)
               dsimfz(:) = dsimfz(:) + ns(ib)*aux%im(:)
            
            end if

            ! Not-diagonal in spin

            ! |a> = |+>,  |b> = |->
            ! Condition below ensures no calc when no |a> with the required spin exists
            if ((ns(ib) == -1) .and. (num_spin_up(ix)>0)) then ! (*)

               ! Auxiliary array
               wf_wf_prod%re(:) = wfa_rhoab(-ns(ib))%re(:,iz)*wf_b(:)
               wf_wf_prod%im(:) = wfa_rhoab(-ns(ib))%im(:,iz)*wf_b(:)

               ! s_rho: rho-component of vector density
               dsresr(:) = dsresr(:) + wf_wf_prod%re(:)
               dsimsr(:) = dsimsr(:) + wf_wf_prod%im(:)

               ! s_phi: phi-component of vector density
               dsresp(:) = dsresp(:) - wf_wf_prod%im(:)
               dsimsp(:) = dsimsp(:) + wf_wf_prod%re(:)

               ! J_{phi,rho}: (phi,rho)-component of tensor density [= J_{2,1}/rho]
               aux%re(:) = (dp_wfa_rhoab(-ns(ib))%re(:,iz)*wf_b(:) + wfa_rhoab(-ns(ib))%re(:,iz)*dp_wf_b(:))/2
               aux%im(:) = (dp_wfa_rhoab(-ns(ib))%im(:,iz)*wf_b(:) + wfa_rhoab(-ns(ib))%im(:,iz)*dp_wf_b(:))/2
               dsretjpr(:) = dsretjpr(:) + aux%re(:)
               dsimtjpr(:) = dsimtjpr(:) + aux%im(:)

               ! J_{phi,phi}: (phi,phi)-component of tensor density [= J_{2,2}/rho^2]
               ! aux%re(:) = (dp_wfa_rhoab(-ns(ib))%re(:,iz)*wf_b(:) + wfa_rhoab(-ns(ib))%re(:,iz)*dp_wf_b(:))/2
               ! aux%im(:) = (dp_wfa_rhoab(-ns(ib))%im(:,iz)*wf_b(:) + wfa_rhoab(-ns(ib))%im(:,iz)*dp_wf_b(:))/2
               dsretjpp(:) = dsretjpp(:) - aux%im(:)
               dsimtjpp(:) = dsimtjpp(:) + aux%re(:)

               ! Del.s: phi-component of divergence of s
               aux%re(:) = wfa_rhoab(-ns(ib))%re(:,iz)*dp_wf_b(:) - dp_wfa_rhoab(-ns(ib))%re(:,iz)*wf_b(:) !&
                        !  - y(:)*wf_wf_prod%re(:)
               aux%im(:) = wfa_rhoab(-ns(ib))%im(:,iz)*dp_wf_b(:) - dp_wfa_rhoab(-ns(ib))%im(:,iz)*wf_b(:) !&
                        !  - y(:)*wf_wf_prod%im(:)
               dsregs(:) = dsregs(:) + aux%re(:)
               dsimgs(:) = dsimgs(:) + aux%im(:)

               ! T_rho: rho-component of spin-kinetic density
               aux%re(:) = dr_wfa_rhoab(-ns(ib))%re(:,iz)*dr_wf_b(:) + dp_wfa_rhoab(-ns(ib))%re(:,iz)*dp_wf_b(:) &
                        + dz_wfa_rhoab(-ns(ib))%re(:,iz)*dz_wf_b(:)
               aux%im(:) = dr_wfa_rhoab(-ns(ib))%im(:,iz)*dr_wf_b(:) + dp_wfa_rhoab(-ns(ib))%im(:,iz)*dp_wf_b(:) &
                        + dz_wfa_rhoab(-ns(ib))%im(:,iz)*dz_wf_b(:)
               dsretr(:) = dsretr(:) + aux%re(:)
               dsimtr(:) = dsimtr(:) + aux%im(:)

               ! T_phi: phi-component of spin-kinetic density
               ! aux%re(:) = dr_wfa_rhoab(-ns(ib))%re(:,iz)*dr_wf_b(:) + dp_wfa_rhoab(-ns(ib))%re(:,iz)*dp_wf_b(:) &
               !           + dz_wfa_rhoab(-ns(ib))%re(:,iz)*dz_wf_b(:)
               ! aux%im(:) = dr_wfa_rhoab(-ns(ib))%im(:,iz)*dr_wf_b(:) + dp_wfa_rhoab(-ns(ib))%im(:,iz)*dp_wf_b(:) &
               !           + dz_wfa_rhoab(-ns(ib))%im(:,iz)*dz_wf_b(:)
               dsretp(:) = dsretp(:) - aux%im(:)
               dsimtp(:) = dsimtp(:) + aux%re(:)

               ! J_{rho,rho}: (rho,rho)-component of tensor density
               aux%re(:) = (wfa_rhoab(-ns(ib))%re(:,iz)*dr_wf_b(:) - dr_wfa_rhoab(-ns(ib))%re(:,iz)*wf_b(:))/2
               aux%im(:) = (wfa_rhoab(-ns(ib))%im(:,iz)*dr_wf_b(:) - dr_wfa_rhoab(-ns(ib))%im(:,iz)*wf_b(:))/2
               dsretjrr(:) = dsretjrr(:) - aux%im(:)
               dsimtjrr(:) = dsimtjrr(:) + aux%re(:)

               ! J_{rho,phi}: (rho,phi)-component of tensor density [= J_{1,2}/rho]
               ! aux%re(:) = (wfa_rhoab(-ns(ib))%re(:,iz)*dr_wf_b(:) - dr_wfa_rhoab(-ns(ib))%re(:,iz)*wf_b(:))/2
               ! aux%im(:) = (wfa_rhoab(-ns(ib))%im(:,iz)*dr_wf_b(:) - dr_wfa_rhoab(-ns(ib))%im(:,iz)*wf_b(:))/2
               dsretjrp(:) = dsretjrp(:) - aux%re(:)
               dsimtjrp(:) = dsimtjrp(:) - aux%im(:)

               ! Del.s: rho-component of divergence of s
               aux%re(:) = wfa_rhoab(-ns(ib))%re(:,iz)*dr_wf_b(:) + dr_wfa_rhoab(-ns(ib))%re(:,iz)*wf_b(:) !&
                        !  + y(:)*wf_wf_prod%re(:)
               aux%im(:) = wfa_rhoab(-ns(ib))%im(:,iz)*dr_wf_b(:) + dr_wfa_rhoab(-ns(ib))%im(:,iz)*wf_b(:) !&
                        !  + y(:)*wf_wf_prod%im(:)
               dsregs(:) = dsregs(:) + aux%re(:)
               dsimgs(:) = dsimgs(:) + aux%im(:)

               ! J_{z,rho}: (z,rho)-component of tensor density
               aux%re(:) = (wfa_rhoab(-ns(ib))%re(:,iz)*dz_wf_b(:) - dz_wfa_rhoab(-ns(ib))%re(:,iz)*wf_b(:))/2
               aux%im(:) = (wfa_rhoab(-ns(ib))%im(:,iz)*dz_wf_b(:) - dz_wfa_rhoab(-ns(ib))%im(:,iz)*wf_b(:))/2
               dsretjzr(:) = dsretjzr(:) - aux%im(:)
               dsimtjzr(:) = dsimtjzr(:) + aux%re(:)

               ! J_{z,phi}: (z,phi)-component of tensor density [= J_{3,2}/rho]
               ! aux%re(:) = (wfa_rhoab(-ns(ib))%re(:,iz)*dz_wf_b(:) - dz_wfa_rhoab(-ns(ib))%re(:,iz)*wf_b(:))/2
               ! aux%im(:) = (wfa_rhoab(-ns(ib))%im(:,iz)*dz_wf_b(:) - dz_wfa_rhoab(-ns(ib))%im(:,iz)*wf_b(:))/2
               dsretjzp(:) = dsretjzp(:) - aux%re(:)
               dsimtjzp(:) = dsimtjzp(:) - aux%im(:)

               ! F_rho: rho-component of tensor-kinetic density
               aux%re(:) = dr_wfa_rhoab(-ns(ib))%re(:,iz)*(dr_wf_b(:) + 0.5_dp*dp_wf_b(:)) &
                        - 0.5_dp*dp_wfa_rhoab(-ns(ib))%re(:,iz)*dr_wf_b(:)
               aux%im(:) = dr_wfa_rhoab(-ns(ib))%im(:,iz)*(dr_wf_b(:) + 0.5_dp*dp_wf_b(:)) &
                        - 0.5_dp*dp_wfa_rhoab(-ns(ib))%im(:,iz)*dr_wf_b(:)
               dsrefr(:) = dsrefr(:) + aux%re(:)
               dsimfr(:) = dsimfr(:) + aux%im(:)

               ! F_phi: phi-component of tensor-kinetic density
               aux%re(:) = dp_wfa_rhoab(-ns(ib))%re(:,iz)*dp_wf_b(:) &
                        + 0.5_dp*(dp_wfa_rhoab(-ns(ib))%re(:,iz)*dr_wf_b(:) - dr_wfa_rhoab(-ns(ib))%re(:,iz)*dp_wf_b(:))
               aux%im(:) = dp_wfa_rhoab(-ns(ib))%im(:,iz)*dp_wf_b(:) &
                        + 0.5_dp*(dp_wfa_rhoab(-ns(ib))%im(:,iz)*dr_wf_b(:) - dr_wfa_rhoab(-ns(ib))%im(:,iz)*dp_wf_b(:))
               dsrefp(:) = dsrefp(:) - aux%im(:)
               dsimfp(:) = dsimfp(:) + aux%re(:)

               ! F_z: z-component of tensor-kinetic density
               aux%re(:) = dz_wfa_rhoab(-ns(ib))%re(:,iz)*(dr_wf_b(:) + dp_wf_b(:)) &
                        + dr_wfa_rhoab(-ns(ib))%re(:,iz)*dz_wf_b(:) - dp_wfa_rhoab(-ns(ib))%re(:,iz)*dz_wf_b(:)
               aux%im(:) = dz_wfa_rhoab(-ns(ib))%im(:,iz)*(dr_wf_b(:) + dp_wf_b(:)) &
                        + dr_wfa_rhoab(-ns(ib))%im(:,iz)*dz_wf_b(:) - dp_wfa_rhoab(-ns(ib))%im(:,iz)*dz_wf_b(:)
               dsrefz(:) = dsrefz(:) + aux%re(:)/2
               dsimfz(:) = dsimfz(:) + aux%im(:)/2

            ! |a> = |->,  |b> = |+>
            ! Condition below ensures no calc when no |a> with the required spin exists
            else if ((ns(ib) == +1) .and. (num_spin_dw(ix)>0)) then ! (*)

               ! Auxiliary array
               wf_wf_prod%re(:) = wfa_rhoab(-ns(ib))%re(:,iz)*wf_b(:)
               wf_wf_prod%im(:) = wfa_rhoab(-ns(ib))%im(:,iz)*wf_b(:)

               ! s_rho: rho-component of vector density
               dsresr(:) = dsresr(:) + wf_wf_prod%re(:)
               dsimsr(:) = dsimsr(:) + wf_wf_prod%im(:)

               ! s_phi: phi-component of vector density
               dsresp(:) = dsresp(:) + wf_wf_prod%im(:)
               dsimsp(:) = dsimsp(:) - wf_wf_prod%re(:)

               ! J_{phi,rho}: (phi,rho)-component of tensor density [= J_{2,1}/rho]
               aux%re(:) = (dp_wfa_rhoab(-ns(ib))%re(:,iz)*wf_b(:) + wfa_rhoab(-ns(ib))%re(:,iz)*dp_wf_b(:))/2
               aux%im(:) = (dp_wfa_rhoab(-ns(ib))%im(:,iz)*wf_b(:) + wfa_rhoab(-ns(ib))%im(:,iz)*dp_wf_b(:))/2
               dsretjpr(:) = dsretjpr(:) + aux%re(:)
               dsimtjpr(:) = dsimtjpr(:) + aux%im(:)

               ! J_{phi,phi}: (phi,phi)-component of tensor density [= J_{2,2}/rho]
               ! aux%re(:) = (dp_wfa_rhoab(-ns(ib))%re(:,iz)*wf_b(:) + wfa_rhoab(-ns(ib))%re(:,iz)*dp_wf_b(:))/2
               ! aux%im(:) = (dp_wfa_rhoab(-ns(ib))%im(:,iz)*wf_b(:) + wfa_rhoab(-ns(ib))%im(:,iz)*dp_wf_b(:))/2
               dsretjpp(:) = dsretjpp(:) + aux%im(:)
               dsimtjpp(:) = dsimtjpp(:) - aux%re(:)

               ! Del.s: phi-component of divergence of s
               aux%re(:) = wfa_rhoab(-ns(ib))%re(:,iz)*dp_wf_b(:) - dp_wfa_rhoab(-ns(ib))%re(:,iz)*wf_b(:) !&
                        !  + y(:)*wf_wf_prod%re(:)
               aux%im(:) = wfa_rhoab(-ns(ib))%im(:,iz)*dp_wf_b(:) - dp_wfa_rhoab(-ns(ib))%im(:,iz)*wf_b(:) !&
                        !  + y(:)*wf_wf_prod%im(:)
               dsregs(:) = dsregs(:) - aux%re(:)
               dsimgs(:) = dsimgs(:) - aux%im(:)

               ! T_rho: rho-component of spin-kinetic density
               aux%re(:) = dr_wfa_rhoab(-ns(ib))%re(:,iz)*dr_wf_b(:) + dp_wfa_rhoab(-ns(ib))%re(:,iz)*dp_wf_b(:) &
                        + dz_wfa_rhoab(-ns(ib))%re(:,iz)*dz_wf_b(:)
               aux%im(:) = dr_wfa_rhoab(-ns(ib))%im(:,iz)*dr_wf_b(:) + dp_wfa_rhoab(-ns(ib))%im(:,iz)*dp_wf_b(:) &
                        + dz_wfa_rhoab(-ns(ib))%im(:,iz)*dz_wf_b(:)
               dsretr(:) = dsretr(:) + aux%re(:)
               dsimtr(:) = dsimtr(:) + aux%im(:)

               ! T_phi: phi-component of spin-kinetic density
               ! aux%re(:) = dr_wfa_rhoab(-ns(ib))%re(:,iz)*dr_wf_b(:) + dp_wfa_rhoab(-ns(ib))%re(:,iz)*dp_wf_b(:) &
               !           + dz_wfa_rhoab(-ns(ib))%re(:,iz)*dz_wf_b(:)
               ! aux%im(:) = dr_wfa_rhoab(-ns(ib))%im(:,iz)*dr_wf_b(:) + dp_wfa_rhoab(-ns(ib))%im(:,iz)*dp_wf_b(:) &
               !           + dz_wfa_rhoab(-ns(ib))%im(:,iz)*dz_wf_b(:)
               dsretp(:) = dsretp(:) + aux%im(:)
               dsimtp(:) = dsimtp(:) - aux%re(:)

               ! J_{rho,rho}: (rho,rho)-component of tensor density
               aux%re(:) = (wfa_rhoab(-ns(ib))%re(:,iz)*dr_wf_b(:) - dr_wfa_rhoab(-ns(ib))%re(:,iz)*wf_b(:))/2
               aux%im(:) = (wfa_rhoab(-ns(ib))%im(:,iz)*dr_wf_b(:) - dr_wfa_rhoab(-ns(ib))%im(:,iz)*wf_b(:))/2
               dsretjrr(:) = dsretjrr(:) - aux%im(:)
               dsimtjrr(:) = dsimtjrr(:) + aux%re(:)

               ! J_{rho,phi}: (rho,phi)-component of tensor density [= J_{1,2}/rho]
               ! aux%re(:) = (wfa_rhoab(-ns(ib))%re(:,iz)*dr_wf_b(:) - dr_wfa_rhoab(-ns(ib))%re(:,iz)*wf_b(:))/2
               ! aux%im(:) = (wfa_rhoab(-ns(ib))%im(:,iz)*dr_wf_b(:) - dr_wfa_rhoab(-ns(ib))%im(:,iz)*wf_b(:))/2
               dsretjrp(:) = dsretjrp(:) + aux%re(:)
               dsimtjrp(:) = dsimtjrp(:) + aux%im(:)

               ! Del.s: rho-component of divergence of s
               aux%re(:) = wfa_rhoab(-ns(ib))%re(:,iz)*dr_wf_b(:) + dr_wfa_rhoab(-ns(ib))%re(:,iz)*wf_b(:) !&
                        !  + y(:)*wf_wf_prod%re(:)
               aux%im(:) = wfa_rhoab(-ns(ib))%im(:,iz)*dr_wf_b(:) + dr_wfa_rhoab(-ns(ib))%im(:,iz)*wf_b(:) !&
                        !  + y(:)*wf_wf_prod%im(:)
               dsregs(:) = dsregs(:) + aux%re(:)
               dsimgs(:) = dsimgs(:) + aux%im(:)

               ! J_{z,rho}: (z,rho)-component of tensor density
               aux%re(:) = (wfa_rhoab(-ns(ib))%re(:,iz)*dz_wf_b(:) - dz_wfa_rhoab(-ns(ib))%re(:,iz)*wf_b(:))/2
               aux%im(:) = (wfa_rhoab(-ns(ib))%im(:,iz)*dz_wf_b(:) - dz_wfa_rhoab(-ns(ib))%im(:,iz)*wf_b(:))/2
               dsretjzr(:) = dsretjzr(:) - aux%im(:)
               dsimtjzr(:) = dsimtjzr(:) + aux%re(:)

               ! J_{z,phi}: (z,phi)-component of tensor density [= J_{3,2}/rho]
               ! aux%re(:) = (wfa_rhoab(-ns(ib))%re(:,iz)*dz_wf_b(:) - dz_wfa_rhoab(-ns(ib))%re(:,iz)*wf_b(:))/2
               ! aux%im(:) = (wfa_rhoab(-ns(ib))%im(:,iz)*dz_wf_b(:) - dz_wfa_rhoab(-ns(ib))%im(:,iz)*wf_b(:))/2
               dsretjzp(:) = dsretjzp(:) + aux%re(:)
               dsimtjzp(:) = dsimtjzp(:) + aux%im(:)

               ! F_rho: rho-component of tensor-kinetic density
               aux%re(:) = dr_wfa_rhoab(-ns(ib))%re(:,iz)*(dr_wf_b(:) - 0.5_dp*dp_wf_b(:)) &
                        + 0.5_dp*dp_wfa_rhoab(-ns(ib))%re(:,iz)*dr_wf_b(:)
               aux%im(:) = dr_wfa_rhoab(-ns(ib))%im(:,iz)*(dr_wf_b(:) - 0.5_dp*dp_wf_b(:)) &
                        + 0.5_dp*dp_wfa_rhoab(-ns(ib))%im(:,iz)*dr_wf_b(:)
               dsrefr(:) = dsrefr(:) + aux%re(:)
               dsimfr(:) = dsimfr(:) + aux%im(:)

               ! F_phi: phi-component of tensor-kinetic density
               aux%re(:) = -dp_wfa_rhoab(-ns(ib))%re(:,iz)*dp_wf_b(:) &
                        + 0.5_dp*(dp_wfa_rhoab(-ns(ib))%re(:,iz)*dr_wf_b(:) - dr_wfa_rhoab(-ns(ib))%re(:,iz)*dp_wf_b(:))
               aux%im(:) = -dp_wfa_rhoab(-ns(ib))%im(:,iz)*dp_wf_b(:) &
                        + 0.5_dp*(dp_wfa_rhoab(-ns(ib))%im(:,iz)*dr_wf_b(:) - dr_wfa_rhoab(-ns(ib))%im(:,iz)*dp_wf_b(:))
               dsrefp(:) = dsrefp(:) - aux%im(:)
               dsimfp(:) = dsimfp(:) + aux%re(:)

               ! F_z: z-component of tensor-kinetic density
               aux%re(:) = dz_wfa_rhoab(-ns(ib))%re(:,iz)*(dr_wf_b(:) - dp_wf_b(:)) &
                        + dr_wfa_rhoab(-ns(ib))%re(:,iz)*dz_wf_b(:) + dp_wfa_rhoab(-ns(ib))%re(:,iz)*dz_wf_b(:)
               aux%im(:) = dz_wfa_rhoab(-ns(ib))%im(:,iz)*(dr_wf_b(:) - dp_wf_b(:)) &
                        + dr_wfa_rhoab(-ns(ib))%im(:,iz)*dz_wf_b(:) + dp_wfa_rhoab(-ns(ib))%im(:,iz)*dz_wf_b(:)
               dsrefz(:) = dsrefz(:) + aux%re(:)/2
               dsimfz(:) = dsimfz(:) + aux%im(:)/2

            ! else
               ! call abort(' Inconsistent spin in helpers.density')
            end if

         end do ! iz, sp state in one block
         !!$OMP END PARALLEL DO

      end do ! ix, block index

      ! The same for pairing densities (separate loop because kappas have
      ! a different block structure than rho)
      
      !-------------------------------------------------------------------------
      ! Pairing densities in the breve representation (Perlinska et al).
      ! rho-breve (rb): isovector pairing
      ! s-breve (sb):   isoscalar pairing
      !
      ! \breve\rho(r s t, r' s' t') = (-2s')(-2t') \sum_{ab} \kappa_{ab}
      !                                   x \phi_a(r s t) \phi_b(r' -s' -t')
      !-------------------------------------------------------------------------

      ! Loop over blocks
      do ix = 1, nb
         iy = rek%ir2c(ix)
         if (iy==0) cycle

         ! Calculate Sum_a wf_a * kappa_ab
         ! Use OpenMP-enabled BLAS library to achieve best performance
         ia = isstart(ix)

         ! wf_a has spin up
         if (num_spin_up(ix)>0) then
            call dgemm('N','N',nghl,db(iy),num_spin_up(ix),1.0_dp,wf(1,ia),nghl,&
                       rek%elem(rek%ir2m(ix)),db(ix),0.0_dp,wfa_rhoab(1)%re(1,1),nghl)
            call dgemm('N','N',nghl,db(iy),num_spin_up(ix),1.0_dp,wf(1,ia),nghl,&
                       imk%elem(imk%ir2m(ix)),db(ix),0.0_dp,wfa_rhoab(1)%im(1,1),nghl)
         end if

         ! wf_a has spin down
         if (num_spin_dw(ix)>0) then
            ia = ia + num_spin_up(ix)
            call dgemm('N','N',nghl,db(iy),num_spin_dw(ix),1.0_dp,wf(1,ia),nghl,&
                       rek%elem(rek%ir2m(ix)+num_spin_up(ix)),db(ix),0.0_dp,wfa_rhoab(-1)%re(1,1),nghl)
            call dgemm('N','N',nghl,db(iy),num_spin_dw(ix),1.0_dp,wf(1,ia),nghl,&
                       imk%elem(imk%ir2m(ix)+num_spin_up(ix)),db(ix),0.0_dp,wfa_rhoab(-1)%im(1,1),nghl)
         end if

         ! Calculate local pairing densities in coordinate space
         !!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(iy,isstart,db,ns,wf,wfa_rhoab)    &
         !!$OMP& REDUCTION(+:dsrerb)  REDUCTION(+:dsimrb)  REDUCTION(+:dsresbr) REDUCTION(+:dsimsbr) &
         !!$OMP& REDUCTION(+:dsresbp) REDUCTION(+:dsimsbp) REDUCTION(+:dsresbz) REDUCTION(+:dsimsbz)
         do iz = 1, db(iy)

            ! Column vector |b>
            ib = isstart(iy)+iz-1
            wf_b => wf(:, ib)

            ! Diagonal in spin
            ! ns(ia) == -ns(ib), note the sign flip w.r.t. normal densities
            ! |a> = |+>,  |b> = |->; or |a> = |->,  |b> = |+>
            ! Condition below ensures no calc when no |a> with the required spin exists
            if (((ns(ib) == -1) .and. (num_spin_up(ix)>0)) .or. &
                ((ns(ib) == +1) .and. (num_spin_dw(ix)>0))) then ! (*)
               
               ! Auxiliary array
               wf_wf_prod%re = 2*wfa_rhoab(-ns(ib))%re(:,iz)*wf_b
               wf_wf_prod%im = 2*wfa_rhoab(-ns(ib))%im(:,iz)*wf_b

               ! rho-breve_10
               dsrerb = dsrerb - ns(ib)*wf_wf_prod%re
               dsimrb = dsimrb - ns(ib)*wf_wf_prod%im
                  
               ! s-breve_0,z
               dsresbz = dsresbz + wf_wf_prod%re
               dsimsbz = dsimsbz + wf_wf_prod%im

            end if

            ! Not-diagonal in spin

            ! |a> = |->,  |b> = |->, note the sign flip w.r.t. normal densities
            ! Condition below ensures no calc when no |a> with the required spin exists
            if ((ns(ib) == -1) .and. (num_spin_dw(ix)>0)) then ! (*)

               ! Auxiliary array 
               wf_wf_prod%re = 2*wfa_rhoab(ns(ib))%re(:,iz)*wf_b
               wf_wf_prod%im = 2*wfa_rhoab(ns(ib))%im(:,iz)*wf_b

               ! s-breve_0,rho
               dsresbr = dsresbr + wf_wf_prod%re
               dsimsbr = dsimsbr + wf_wf_prod%im
                  
               ! s-breve_0,phi
               dsresbp = dsresbp + wf_wf_prod%im
               dsimsbp = dsimsbp - wf_wf_prod%re

            ! |a> = |+>, |b> = |+>, note the sign flip w.r.t. normal densities
            ! Condition below ensures no calc when no |a> with the required spin exists
            else if ((ns(ib) == +1) .and. (num_spin_up(ix)>0)) then ! (*)

               ! Auxiliary array
               wf_wf_prod%re = 2*wfa_rhoab(ns(ib))%re(:,iz)*wf_b
               wf_wf_prod%im = 2*wfa_rhoab(ns(ib))%im(:,iz)*wf_b

               ! s-breve_0,rho
               dsresbr = dsresbr - wf_wf_prod%re
               dsimsbr = dsimsbr - wf_wf_prod%im
                  
               ! s-breve_0,phi
               dsresbp = dsresbp + wf_wf_prod%im
               dsimsbp = dsimsbp - wf_wf_prod%re

            ! else
            !    call abort(' Inconsistent spin in helpers.density')
            end if
            
         end do ! iz, sp state in one block
         !!$OMP END PARALLEL DO
      end do ! ix, block index

      ! Strip off the integral weights
      ds%rerho(:) = dsrerho(:)*wdcori(:);    ds%imrho(:) = dsimrho(:)*wdcori(:)
      ds%retau(:) = dsretau(:)*wdcori(:);    ds%imtau(:) = dsimtau(:)*wdcori(:)
      ds%retjrr(:) = dsretjrr(:)*wdcori(:);  ds%imtjrr(:) = dsimtjrr(:)*wdcori(:)
      ds%retjpr(:) = dsretjpr(:)*wdcori(:);  ds%imtjpr(:) = dsimtjpr(:)*wdcori(:)
      ds%retjzr(:) = dsretjzr(:)*wdcori(:);  ds%imtjzr(:) = dsimtjzr(:)*wdcori(:)
      ds%retjrp(:) = dsretjrp(:)*wdcori(:);  ds%imtjrp(:) = dsimtjrp(:)*wdcori(:)
      ds%retjpp(:) = dsretjpp(:)*wdcori(:);  ds%imtjpp(:) = dsimtjpp(:)*wdcori(:)
      ds%retjzp(:) = dsretjzp(:)*wdcori(:);  ds%imtjzp(:) = dsimtjzp(:)*wdcori(:)
      ds%retjrz(:) = dsretjrz(:)*wdcori(:);  ds%imtjrz(:) = dsimtjrz(:)*wdcori(:)
      ds%retjpz(:) = dsretjpz(:)*wdcori(:);  ds%imtjpz(:) = dsimtjpz(:)*wdcori(:)
      ds%retjzz(:) = dsretjzz(:)*wdcori(:);  ds%imtjzz(:) = dsimtjzz(:)*wdcori(:)

      ds%resr(:) = dsresr(:)*wdcori(:);  ds%imsr(:) = dsimsr(:)*wdcori(:)
      ds%resp(:) = dsresp(:)*wdcori(:);  ds%imsp(:) = dsimsp(:)*wdcori(:)
      ds%resz(:) = dsresz(:)*wdcori(:);  ds%imsz(:) = dsimsz(:)*wdcori(:)
      ds%retr(:) = dsretr(:)*wdcori(:);  ds%imtr(:) = dsimtr(:)*wdcori(:)
      ds%retp(:) = dsretp(:)*wdcori(:);  ds%imtp(:) = dsimtp(:)*wdcori(:)
      ds%retz(:) = dsretz(:)*wdcori(:);  ds%imtz(:) = dsimtz(:)*wdcori(:)
      ds%rejr(:) = dsrejr(:)*wdcori(:);  ds%imjr(:) = dsimjr(:)*wdcori(:)
      ds%rejp(:) = dsrejp(:)*wdcori(:);  ds%imjp(:) = dsimjp(:)*wdcori(:)
      ds%rejz(:) = dsrejz(:)*wdcori(:);  ds%imjz(:) = dsimjz(:)*wdcori(:)
      ds%refr(:) = dsrefr(:)*wdcori(:);  ds%imfr(:) = dsimfr(:)*wdcori(:)
      ds%refp(:) = dsrefp(:)*wdcori(:);  ds%imfp(:) = dsimfp(:)*wdcori(:)
      ds%refz(:) = dsrefz(:)*wdcori(:);  ds%imfz(:) = dsimfz(:)*wdcori(:)
      ds%regs(:) = dsregs(:)*wdcori(:);  ds%imgs(:) = dsimgs(:)*wdcori(:)

      ! Pairing
      ds%rerb(:)  = dsrerb(:)*wdcori(:);   ds%imrb(:)  = dsimrb(:)*wdcori(:)
      ds%resbr(:) = dsresbr(:)*wdcori(:);  ds%imsbr(:) = dsimsbr(:)*wdcori(:)
      ds%resbp(:) = dsresbp(:)*wdcori(:);  ds%imsbp(:) = dsimsbp(:)*wdcori(:)
      ds%resbz(:) = dsresbz(:)*wdcori(:);  ds%imsbz(:) = dsimsbz(:)*wdcori(:)

   end subroutine density


   !---------------------------------------------------------------------------
   ! Computation of the Hamiltonian matrix from local densities
   !---------------------------------------------------------------------------
   subroutine meanfield(ds, reh, imh)
      use pnfam_interaction, only: crho, cs, ctau, cj, ct, cdrho, cds, crdj, csdj, ctj0, ctj1, ctj2, cf, cgs

      implicit none

      type(blockmatrix), intent(inout) :: reh, imh
      type(density_set), intent(in) :: ds

      integer  :: ia, ib, ix, iy, iz, ixx, iyy, izz
      real(dp) :: coeff
      
      type pointer_single
         real(dp), pointer :: p
      end type pointer_single

      type pointer_1d
         real(dp), pointer, dimension(:) :: p
      end type pointer_1d

      ! Arrays for mean fields
      ! first index: 0 wf_a, 1 dr_wf_a, 2 dp_wf_a = nl(ia)*y*wf_a, 3 dz_wf_a, 4 d2_wf_a_all
      ! second index: 0 wf_b, 1 dr_wf_b, 2 dp_wf_b = nl(ib)*y*wf_b, 3 dz_wf_b, 4 d2_wf_b_all
      ! third index: ns(ia). fourth index: ns(ib). (0 is not used in last two indices)
      type(field), save :: mf(0:4,0:4,-1:1,-1:1)
      type(array_2d), save :: hpsi(0:4,-1:1) ! Result of mf applied on column vector |b>
      type(pointer_single) :: wf_a_start(0:4)
      type(pointer_1d) :: wf_b_all(0:4) !, wf_a_all(0:4)
      type(field), save :: aux_fd
      integer, parameter :: max_ind(0:4) = (/ 4, 3, 3, 3, 0 /) ! max index for each row / col of mf / hpsi

      ! Optimized tensor force calculation
      type(field), save :: ctj0_tjrr_tjpp_tjzz, ctj1_tjzr_minus_tjrz, ctj1_tjpz_minus_tjzp, ctj2_tjrz_tjzr, ctj2_tjpz_tjzp

      ! reh%elem = 0; imh%elem = 0 ! Zeroing h is unnecessary 

      ! Allocate mean fields
      if (.not. allocated(aux_fd%re)) then
         call allocate_field(aux_fd, nghl)
         call allocate_field(ctj0_tjrr_tjpp_tjzz, nghl)
         call allocate_field(ctj1_tjzr_minus_tjrz, nghl)
         call allocate_field(ctj1_tjpz_minus_tjzp, nghl)
         call allocate_field(ctj2_tjrz_tjzr, nghl)
         call allocate_field(ctj2_tjpz_tjzp, nghl)
         do iy = -1, 1, 2 
            do ix = -1, 1, 2
               do ib = 0, 4
                  do ia = 0, max_ind(ib)
                     call allocate_field(mf(ia,ib,ix,iy),nghl)
                  end do
               end do
            end do
         end do
         do ixx = -1, 1, 2
            do iyy = 0, 4
               call allocate_2darray(hpsi(iyy,ixx),nghl,maxval_db)
            end do
         end do
      end if

      ! Optimized tensor force calculation
      ctj0_tjrr_tjpp_tjzz%re = ctj0*(ds%retjrr(:)+ds%retjpp(:)+ds%retjzz(:))
      ctj0_tjrr_tjpp_tjzz%im = ctj0*(ds%imtjrr(:)+ds%imtjpp(:)+ds%imtjzz(:))
      ctj1_tjzr_minus_tjrz%re = ctj1*(ds%retjzr(:)-ds%retjrz(:))
      ctj1_tjzr_minus_tjrz%im = ctj1*(ds%imtjzr(:)-ds%imtjrz(:))
      ctj1_tjpz_minus_tjzp%re = ctj1*(ds%retjpz(:)-ds%retjzp(:))
      ctj1_tjpz_minus_tjzp%im = ctj1*(ds%imtjpz(:)-ds%imtjzp(:))
      ctj2_tjrz_tjzr%re = ctj2*(ds%retjrz(:)+ds%retjzr(:))
      ctj2_tjrz_tjzr%im = ctj2*(ds%imtjrz(:)+ds%imtjzr(:))
      ctj2_tjpz_tjzp%re = ctj2*(ds%retjpz(:)+ds%retjzp(:))
      ctj2_tjpz_tjzp%im = ctj2*(ds%imtjpz(:)+ds%imtjzp(:))

      ! Calculate mean field: wf_a, wf_b, ns(ia)==ns(ib)
      mf(0,0,1,1)%re = 2.0_dp*crho(:)*ds%rerho(:) + ctau*ds%retau(:)
      mf(0,0,1,1)%im = 2.0_dp*crho(:)*ds%imrho(:) + ctau*ds%imtau(:)
      mf(0,0,-1,-1) = mf(0,0,1,1)
      aux_fd%re = 2.0_dp*cs(:)*ds%resz(:) + ct*ds%retz(:) + cf*ds%refz(:)
      aux_fd%im = 2.0_dp*cs(:)*ds%imsz(:) + ct*ds%imtz(:) + cf*ds%imfz(:)
      call add_field(mf(0,0,1,1), 1, aux_fd)
      call add_field(mf(0,0,-1,-1), -1, aux_fd)

      ! Calculate mean field: wf_a, dr_wf_b, ns(ia)==ns(ib); dr_wf_a, wf_b, ns(ia) == ns(ib)
      mf(0,1,1,1)%re = -crdj*(ds%retjpz(:)-ds%retjzp(:))
      mf(0,1,1,1)%im = -crdj*(ds%imtjpz(:)-ds%imtjzp(:))
      mf(0,1,-1,-1) = mf(0,1,1,1)
      aux_fd%re = -csdj*ds%rejp(:)
      aux_fd%im = -csdj*ds%imjp(:)
      call add_field(mf(0,1,1,1), 1, aux_fd)
      call add_field(mf(0,1,-1,-1), -1, aux_fd)
      mf(1,0,1,1) = mf(0,1,1,1)
      mf(1,0,-1,-1) = mf(0,1,-1,-1)
      ! aux_fd%re = -ctj1*(ds%imtjzr(:)-ds%imtjrz(:))+0.5_dp*ctj2*(ds%imtjrz(:)+ds%imtjzr(:))
      ! aux_fd%im = ctj1*(ds%retjzr(:)-ds%retjrz(:))-0.5_dp*ctj2*(ds%retjrz(:)+ds%retjzr(:))
      aux_fd%re = -ctj1_tjzr_minus_tjrz%im(:)+0.5_dp*ctj2_tjrz_tjzr%im(:)
      aux_fd%im = ctj1_tjzr_minus_tjrz%re(:)-0.5_dp*ctj2_tjrz_tjzr%re(:)
      call add_field(mf(0,1,1,1), 1, aux_fd)
      call add_field(mf(0,1,-1,-1), -1, aux_fd)
      call add_field(mf(1,0,1,1), -1, aux_fd)
      call add_field(mf(1,0,-1,-1), 1, aux_fd)
      aux_fd%re = cj*ds%imjr(:)
      aux_fd%im = -cj*ds%rejr(:)
      call add_field(mf(0,1,1,1), 1, aux_fd)
      call add_field(mf(0,1,-1,-1), 1, aux_fd)
      call add_field(mf(1,0,1,1), -1, aux_fd)
      call add_field(mf(1,0,-1,-1), -1, aux_fd)

      ! Calculate mean field: wf_a, y*nl(ib)*wf_b, ns(ia)==ns(ib); y*nl(ia)*wf_a, wf_b, ns(ia) == ns(ib)
      mf(0,2,1,1)%re = cj*ds%rejp(:)
      mf(0,2,1,1)%im = cj*ds%imjp(:)
      mf(0,2,-1,-1) = mf(0,2,1,1)
      ! aux_fd%re = ctj1*(ds%retjpz(:)-ds%retjzp(:))+0.5_dp*ctj2*(ds%retjpz(:)+ds%retjzp(:))
      ! aux_fd%im = ctj1*(ds%imtjpz(:)-ds%imtjzp(:))+0.5_dp*ctj2*(ds%imtjpz(:)+ds%imtjzp(:))
      aux_fd%re = ctj1_tjpz_minus_tjzp%re(:)+0.5_dp*ctj2_tjpz_tjzp%re(:)
      aux_fd%im = ctj1_tjpz_minus_tjzp%im(:)+0.5_dp*ctj2_tjpz_tjzp%im(:)
      call add_field(mf(0,2,1,1), 1, aux_fd)
      call add_field(mf(0,2,-1,-1), -1, aux_fd)
      mf(2,0,1,1) = mf(0,2,1,1)
      mf(2,0,-1,-1) = mf(0,2,-1,-1)
      aux_fd%re = -csdj*ds%imjr(:)
      aux_fd%im = csdj*ds%rejr(:)
      call add_field(mf(0,2,1,1), 1, aux_fd)
      call add_field(mf(0,2,-1,-1), -1, aux_fd)
      call add_field(mf(2,0,1,1), -1, aux_fd)
      call add_field(mf(2,0,-1,-1), 1, aux_fd)
      aux_fd%re = crdj*(ds%imtjzr(:)-ds%imtjrz(:))
      aux_fd%im = -crdj*(ds%retjzr(:)-ds%retjrz(:))
      call add_field(mf(0,2,1,1), 1, aux_fd)
      call add_field(mf(0,2,-1,-1), 1, aux_fd)
      call add_field(mf(2,0,1,1), -1, aux_fd)
      call add_field(mf(2,0,-1,-1), -1, aux_fd)

      ! Calculate mean field: wf_a, dz_wf_b, ns(ia)==ns(ib); dz_wf_a, wf_b, ns(ia) == ns(ib)
      mf(0,3,1,1)%re = -crdj*(ds%retjrp(:)-ds%retjpr(:))
      mf(0,3,1,1)%im = -crdj*(ds%imtjrp(:)-ds%imtjpr(:))
      mf(0,3,-1,-1) = mf(0,3,1,1)
      aux_fd%re = 2.0_dp*cgs*ds%regs(:)
      aux_fd%im = 2.0_dp*cgs*ds%imgs(:)
      call add_field(mf(0,3,1,1), 1, aux_fd)
      call add_field(mf(0,3,-1,-1), -1, aux_fd)
      mf(3,0,1,1) = mf(0,3,1,1)
      mf(3,0,-1,-1) = mf(0,3,-1,-1)
      aux_fd%re = ctj0_tjrr_tjpp_tjzz%im(:) & 
                ! ctj0*(ds%imtjrr(:)+ds%imtjpp(:)+ds%imtjzz(:)) &
                 -ctj2*(ds%imtjrr(:)+ds%imtjpp(:)-2.0_dp*ds%imtjzz(:))/3.0_dp
      aux_fd%im = -ctj0_tjrr_tjpp_tjzz%re(:) &
                ! -ctj0*(ds%retjrr(:)+ds%retjpp(:)+ds%retjzz(:)) &
                  +ctj2*(ds%retjrr(:)+ds%retjpp(:)-2.0_dp*ds%retjzz(:))/3.0_dp
      call add_field(mf(0,3,1,1), 1, aux_fd)
      call add_field(mf(0,3,-1,-1), -1, aux_fd)
      call add_field(mf(3,0,1,1), -1, aux_fd)
      call add_field(mf(3,0,-1,-1), 1, aux_fd)
      aux_fd%re = cj*ds%imjz(:)
      aux_fd%im = -cj*ds%rejz(:)
      call add_field(mf(0,3,1,1), 1, aux_fd)
      call add_field(mf(0,3,-1,-1), 1, aux_fd)
      call add_field(mf(3,0,1,1), -1, aux_fd)
      call add_field(mf(3,0,-1,-1), -1, aux_fd)

      ! Calculate mean field: wf_a, d2_wf_b_all, ns(ia)==ns(ib); d2_wf_a_all, wf_b, ns(ia) == ns(ib)
      mf(0,4,1,1)%re = 2.0_dp*cdrho*ds%rerho(:)
      mf(0,4,1,1)%im = 2.0_dp*cdrho*ds%imrho(:)
      mf(0,4,-1,-1) = mf(0,4,1,1)
      aux_fd%re = 2.0_dp*cds*ds%resz(:)
      aux_fd%im = 2.0_dp*cds*ds%imsz(:)
      call add_field(mf(0,4,1,1), 1, aux_fd)
      call add_field(mf(0,4,-1,-1), -1, aux_fd)
      mf(4,0,1,1) = mf(0,4,1,1)
      mf(4,0,-1,-1) = mf(0,4,-1,-1)

      ! Calculate mean field: dr_wf_a, y*nl(ib)*wf_b, ns(ia)==ns(ib); nl(ia)*y*wf_a, dr_wf_b, ns(ia) == ns(ib)
      mf(1,2,1,1)%re = csdj*ds%resz(:)
      mf(1,2,1,1)%im = csdj*ds%imsz(:)
      mf(1,2,-1,-1) = mf(1,2,1,1)
      aux_fd%re = crdj*ds%rerho(:)
      aux_fd%im = crdj*ds%imrho(:)
      call add_field(mf(1,2,1,1), 1, aux_fd)
      call add_field(mf(1,2,-1,-1), -1, aux_fd)
      mf(2,1,1,1) = mf(1,2,1,1)
      mf(2,1,-1,-1) = mf(1,2,-1,-1)

      ! Calculate mean field: dr_wf_a, dz_wf_b, ns(ia)==ns(ib); dz_wf_a, dr_wf_b, ns(ia) == ns(ib)
      mf(1,3,1,1)%re = 0.5_dp*cf*ds%resr(:)
      mf(1,3,1,1)%im = 0.5_dp*cf*ds%imsr(:)
      mf(1,3,-1,-1)%re = -mf(1,3,1,1)%re
      mf(1,3,-1,-1)%im = -mf(1,3,1,1)%im
      mf(3,1,1,1) = mf(1,3,1,1)
      mf(3,1,-1,-1) = mf(1,3,-1,-1)
      aux_fd%re = -csdj*ds%imsp(:)
      aux_fd%im = csdj*ds%resp(:)
      call add_field(mf(1,3,1,1), 1, aux_fd)
      call add_field(mf(1,3,-1,-1), 1, aux_fd)
      call add_field(mf(3,1,1,1), -1, aux_fd)
      call add_field(mf(3,1,-1,-1), -1, aux_fd)

      ! Calculate mean field: y*nl(ia)*wf_a, dz_wf_b, ns(ia)==ns(ib); dz_wf_a, y*nl(ib)*wf_b, ns(ia) == ns(ib)
      mf(2,3,1,1)%re = -csdj*ds%resr(:)
      mf(2,3,1,1)%im = -csdj*ds%imsr(:)
      mf(2,3,-1,-1) = mf(2,3,1,1)
      aux_fd%re = 0.5_dp*cf*ds%imsp(:)
      aux_fd%im = -0.5_dp*cf*ds%resp(:)
      call add_field(mf(2,3,1,1), 1, aux_fd)
      call add_field(mf(2,3,-1,-1), -1, aux_fd)
      mf(3,2,1,1) = mf(2,3,-1,-1)
      mf(3,2,-1,-1) = mf(2,3,1,1)

      ! Calculate mean field: dr_wf_a, dr_wf_b, ns(ia)==ns(ib); nl(ia)*y*wf_a, nl(ib)*y*wf_b, ns(ia)==ns(ib); 
      ! and dz_wf_a, dz_wf_b, ns(ia)==ns(ib)
      mf(1,1,1,1)%re = (4.0_dp*cdrho+ctau)*ds%rerho(:)
      mf(1,1,1,1)%im = (4.0_dp*cdrho+ctau)*ds%imrho(:)
      mf(1,1,-1,-1) = mf(1,1,1,1)
      aux_fd%re = (4.0_dp*cds+ct)*ds%resz(:)
      aux_fd%im = (4.0_dp*cds+ct)*ds%imsz(:)
      call add_field(mf(1,1,1,1), 1, aux_fd)
      call add_field(mf(1,1,-1,-1), -1, aux_fd)
      mf(2,2,1,1) = mf(1,1,1,1)
      mf(2,2,-1,-1) = mf(1,1,-1,-1)
      mf(3,3,1,1) = mf(1,1,1,1)
      mf(3,3,-1,-1) = mf(1,1,-1,-1)
      aux_fd%re = cf*ds%resz(:)
      aux_fd%im = cf*ds%imsz(:)
      call add_field(mf(3,3,1,1), 1, aux_fd)
      call add_field(mf(3,3,-1,-1), -1, aux_fd)

      ! Calculate mean field: wf_a, wf_b, ns(ia)/=ns(ib)
      mf(0,0,1,-1)%re = 2.0_dp*cs(:)*ds%resr(:)+ct*ds%retr(:)+cf*ds%refr(:)
      mf(0,0,1,-1)%im = 2.0_dp*cs(:)*ds%imsr(:)+ct*ds%imtr(:)+cf*ds%imfr(:)
      mf(0,0,-1,1) = mf(0,0,1,-1)
      aux_fd%re = 2.0_dp*cs(:)*ds%imsp(:)+ct*ds%imtp(:)+cf*ds%imfp(:)
      aux_fd%im = -2.0_dp*cs(:)*ds%resp(:)-ct*ds%retp(:)-cf*ds%refp(:)
      call add_field(mf(0,0,1,-1), 1, aux_fd)
      call add_field(mf(0,0,-1,1), -1, aux_fd)

      ! Calculate mean field: wf_a, dr_wf_b, ns(ia)/=ns(ib); dr_wf_a, wf_b, ns(ia)/=ns(ib)
      ! also part of: wf_a, nl(ib)*y*wf_b, ns(ia)/=ns(ib); nl(ia)*y*wf_a, wf_b, ns(ia)/=ns(ib)
      mf(0,1,1,-1)%re = 2.0_dp*cgs*ds%regs(:)
      mf(0,1,1,-1)%im = 2.0_dp*cgs*ds%imgs(:)
      mf(0,1,-1,1) = mf(0,1,1,-1)
      aux_fd%re = csdj*ds%imjz(:)
      aux_fd%im = -csdj*ds%rejz(:)
      call add_field(mf(0,1,1,-1), 1, aux_fd)
      call add_field(mf(0,1,-1,1), -1, aux_fd)
      mf(1,0,1,-1) = mf(0,1,1,-1)
      mf(1,0,-1,1) = mf(0,1,-1,1)
      mf(0,2,1,-1) = mf(0,1,1,-1)
      mf(2,0,-1,1) = mf(1,0,-1,1)
      mf(2,0,1,-1)%re = -mf(0,2,1,-1)%re; mf(2,0,1,-1)%im = -mf(0,2,1,-1)%im
      mf(0,2,-1,1)%re = -mf(2,0,-1,1)%re; mf(0,2,-1,1)%im = -mf(2,0,-1,1)%im
      aux_fd%re = -ctj1*(ds%retjrp(:)-ds%retjpr(:))
      aux_fd%im = -ctj1*(ds%imtjrp(:)-ds%imtjpr(:))
      call add_field(mf(0,1,1,-1), 1, aux_fd)
      call add_field(mf(0,1,-1,1), -1, aux_fd)
      call add_field(mf(1,0,1,-1), -1, aux_fd)
      call add_field(mf(1,0,-1,1), 1, aux_fd)
      call add_field(mf(0,2,1,-1), 1, aux_fd)
      call add_field(mf(0,2,-1,1), 1, aux_fd)
      call add_field(mf(2,0,1,-1), 1, aux_fd)
      call add_field(mf(2,0,-1,1), 1, aux_fd)
      ! aux_fd%re = ctj0*(ds%imtrjj(:)+ds%imtrpp(:)+ds%imtrzz(:))
      ! aux_fd%im = -ctj0*(ds%retrjj(:)+ds%retrpp(:)+ds%retrzz(:))
      aux_fd%re = ctj0_tjrr_tjpp_tjzz%im(:)
      aux_fd%im = -ctj0_tjrr_tjpp_tjzz%re(:)
      call add_field(mf(0,1,1,-1), 1, aux_fd)
      call add_field(mf(0,1,-1,1), 1, aux_fd)
      call add_field(mf(1,0,1,-1), -1, aux_fd)
      call add_field(mf(1,0,-1,1), -1, aux_fd)
      call add_field(mf(0,2,1,-1), 1, aux_fd)
      call add_field(mf(0,2,-1,1), -1, aux_fd)
      call add_field(mf(2,0,1,-1), 1, aux_fd)
      call add_field(mf(2,0,-1,1), -1, aux_fd)
      aux_fd%re = -0.5_dp*ctj2*(ds%retjrp(:)+ds%retjpr(:))
      aux_fd%im = -0.5_dp*ctj2*(ds%imtjrp(:)+ds%imtjpr(:))
      call add_field(mf(0,1,1,-1), 1, aux_fd)
      call add_field(mf(0,1,-1,1), -1, aux_fd)
      call add_field(mf(1,0,1,-1), -1, aux_fd)
      call add_field(mf(1,0,-1,1), 1, aux_fd)
      call add_field(mf(0,2,1,-1), -1, aux_fd)
      call add_field(mf(0,2,-1,1), -1, aux_fd)
      call add_field(mf(2,0,1,-1), -1, aux_fd)
      call add_field(mf(2,0,-1,1), -1, aux_fd)
      aux_fd%re = -ctj2*(-2.0_dp*ds%imtjrr(:)+ds%imtjpp(:)+ds%imtjzz(:))/3.0_dp
      aux_fd%im = ctj2*(-2.0_dp*ds%retjrr(:)+ds%retjpp(:)+ds%retjzz(:))/3.0_dp
      call add_field(mf(0,1,1,-1), 1, aux_fd)
      call add_field(mf(0,1,-1,1), 1, aux_fd)
      call add_field(mf(1,0,1,-1), -1, aux_fd)
      call add_field(mf(1,0,-1,1), -1, aux_fd)

      ! Remaining part of: wf_a, nl(ib)*y*wf_b, ns(ia)/=ns(ib); nl(ia)*y*wf_a, wf_b, ns(ia)/=ns(ib)
      aux_fd%re = -ctj2*(ds%imtjrr(:)-2.0_dp*ds%imtjpp(:)+ds%imtjzz(:))/3.0_dp
      aux_fd%im = ctj2*(ds%retjrr(:)-2.0_dp*ds%retjpp(:)+ds%retjzz(:))/3.0_dp
      call add_field(mf(0,2,1,-1), 1, aux_fd)
      call add_field(mf(0,2,-1,1), -1, aux_fd)
      call add_field(mf(2,0,1,-1), 1, aux_fd)
      call add_field(mf(2,0,-1,1), -1, aux_fd)

      ! Calculate mean field: wf_a, dz_wf_b, ns(ia)/=ns(ib); dz_wf_a, wf_b, ns(ia)/=ns(ib)
      mf(0,3,1,-1)%re = csdj*ds%rejp(:)
      mf(0,3,1,-1)%im = csdj*ds%imjp(:)
      mf(0,3,-1,1) = mf(0,3,1,-1)
      aux_fd%re = -csdj*ds%imjr(:)
      aux_fd%im = csdj*ds%rejr(:)
      call add_field(mf(0,3,1,-1), 1, aux_fd)
      call add_field(mf(0,3,-1,1), -1, aux_fd)
      mf(3,0,1,-1) = mf(0,3,1,-1)
      mf(3,0,-1,1) = mf(0,3,-1,1)
      ! aux_fd%re = ctj1*(ds%imtjzr(:)-ds%imtjrz(:))+0.5_dp*ctj2*(ds%imtjrz(:)+ds%imtjzr(:))
      ! aux_fd%im = -ctj1*(ds%retjzr(:)-ds%retjrz(:))-0.5_dp*ctj2*(ds%retjrz(:)+ds%retjzr(:))
      aux_fd%re = ctj1_tjzr_minus_tjrz%im(:)+0.5_dp*ctj2_tjrz_tjzr%im(:)
      aux_fd%im = -ctj1_tjzr_minus_tjrz%re(:)-0.5_dp*ctj2_tjrz_tjzr%re(:)
      call add_field(mf(0,3,1,-1), 1, aux_fd)
      call add_field(mf(0,3,-1,1), 1, aux_fd)
      call add_field(mf(3,0,1,-1), -1, aux_fd)
      call add_field(mf(3,0,-1,1), -1, aux_fd)
      ! aux_fd%re = ctj1*(ds%retjpz(:)-ds%retjzp(:))-0.5_dp*ctj2*(ds%retjpz(:)+ds%retjzp(:))
      ! aux_fd%im = ctj1*(ds%imtjpz(:)-ds%imtjzp(:))-0.5_dp*ctj2*(ds%imtjpz(:)+ds%imtjzp(:))
      aux_fd%re = ctj1_tjpz_minus_tjzp%re(:)-0.5_dp*ctj2_tjpz_tjzp%re(:)
      aux_fd%im = ctj1_tjpz_minus_tjzp%im(:)-0.5_dp*ctj2_tjpz_tjzp%im(:)
      call add_field(mf(0,3,1,-1), 1, aux_fd)
      call add_field(mf(0,3,-1,1), -1, aux_fd)
      call add_field(mf(3,0,1,-1), -1, aux_fd)
      call add_field(mf(3,0,-1,1), 1, aux_fd)

      ! Calculate mean field: wf_a, d2_wf_b_all, ns(ia)/=ns(ib); d2_wf_a_all, wf_b, ns(ia)/=ns(ib)
      mf(0,4,1,-1)%re = 2.0_dp*cds*ds%resr
      mf(0,4,1,-1)%im = 2.0_dp*cds*ds%imsr
      mf(0,4,-1,1) = mf(0,4,1,-1)
      aux_fd%re = 2.0_dp*cds*ds%imsp
      aux_fd%im = -2.0_dp*cds*ds%resp
      call add_field(mf(0,4,1,-1), 1, aux_fd)
      call add_field(mf(0,4,-1,1), -1, aux_fd)
      mf(4,0,1,-1) = mf(0,4,1,-1)
      mf(4,0,-1,1) = mf(0,4,-1,1)

      ! Calculate mean field: dr_wf_a, nl(ib)*y*wf_b, ns(ia)/=ns(ib); nl(ia)*y*wf_a, dr_wf_b, ns(ia)/=ns(ib)
      mf(1,2,1,-1)%re = -0.5_dp*cf*ds%imsp(:)
      mf(1,2,1,-1)%im = 0.5_dp*cf*ds%resp(:)
      mf(1,2,-1,1) = mf(1,2,1,-1)
      aux_fd%re = 0.5_dp*cf*ds%resr(:)
      aux_fd%im = 0.5_dp*cf*ds%imsr(:)
      call add_field(mf(1,2,1,-1), 1, aux_fd)
      call add_field(mf(1,2,-1,1), -1, aux_fd)
      mf(2,1,1,-1)%re = -mf(1,2,1,-1)%re; mf(2,1,1,-1)%im = -mf(1,2,1,-1)%im
      mf(2,1,-1,1)%re = -mf(1,2,-1,1)%re; mf(2,1,-1,1)%im = -mf(1,2,-1,1)%im

      ! Calculate mean field: dr_wf_a, dz_wf_a, ns(ia)/=ns(ib); dz_wf_a, dr_wf_b, ns(ia)/=ns(ib)
      ! also nl(ia)*y*wf_a, dz_wf_b, ns(ia)/=ns(ib); dz_wf_a, nl(ib)*y*wf_b, ns(ia)/=ns(ib)
      mf(1,3,1,-1)%re = 0.5_dp*cf*ds%resz(:)
      mf(1,3,1,-1)%im = 0.5_dp*cf*ds%imsz(:)
      mf(1,3,-1,1) = mf(1,3,1,-1)
      aux_fd%re = crdj*ds%rerho(:)
      aux_fd%im = crdj*ds%imrho(:)
      call add_field(mf(1,3,1,-1), 1, aux_fd)
      call add_field(mf(1,3,-1,1), -1, aux_fd)
      mf(3,1,1,-1) = mf(1,3,-1,1)
      mf(3,1,-1,1) = mf(1,3,1,-1)
      mf(3,2,1,-1) = mf(3,1,1,-1)
      mf(2,3,-1,1) = mf(1,3,-1,1)
      mf(2,3,1,-1)%re = -mf(1,3,1,-1)%re; mf(2,3,1,-1)%im = -mf(1,3,1,-1)%im
      mf(3,2,-1,1)%re = -mf(3,1,-1,1)%re; mf(3,2,-1,1)%im = -mf(3,1,-1,1)%im

      ! Calculate mean field: dr_wf_a, dr_wf_b, ns(ia)/=ns(ib); nl(ia)*y*wf_a, nl(ib)*y*wf_b, ns(ia)/=ns(ib); 
      ! and dz_wf_a, dz_wf_b, ns(ia)/=ns(ib)
      mf(3,3,1,-1)%re = (ct+4.0_dp*cds)*ds%resr(:)
      mf(3,3,1,-1)%im = (ct+4.0_dp*cds)*ds%imsr(:)
      mf(3,3,-1,1) = mf(3,3,1,-1)
      aux_fd%re = (ct+4.0_dp*cds)*ds%imsp(:)
      aux_fd%im = -(ct+4.0_dp*cds)*ds%resp(:)
      call add_field(mf(3,3,1,-1), 1, aux_fd)
      call add_field(mf(3,3,-1,1), -1, aux_fd)
      mf(1,1,1,-1) = mf(3,3,1,-1)
      mf(1,1,-1,1) = mf(3,3,-1,1)
      aux_fd%re = cf*ds%resr(:)
      aux_fd%im = cf*ds%imsr(:)
      call add_field(mf(1,1,1,-1), 1, aux_fd)
      call add_field(mf(1,1,-1,1), 1, aux_fd)
      mf(2,2,1,-1) = mf(3,3,1,-1)
      mf(2,2,-1,1) = mf(3,3,-1,1)
      aux_fd%re = cf*ds%imsp(:)
      aux_fd%im = -cf*ds%resp(:)
      call add_field(mf(2,2,1,-1), 1, aux_fd)
      call add_field(mf(2,2,-1,1), -1, aux_fd)


      ! Loop over blocks
      ! Calculate e.g. h_ab = Sum_{sa,sb} Int{ Conjg[wf_a(x,sa)] [ f(x;sa,sb) ] wf_b(x,sb) }

      do ix = 1, nb
         iy = reh%ir2c(ix) ! column index of current block
         if (iy==0) cycle

         ! Calculate [ f(x;sa,sb) ] * wf_b(x,sb)
         !!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ib,ixx,iyy,izz,wf_b_all,coeff)
         do iz = 1, db(iy)
            ! Column vector |b>
            ib = isstart(iy)+iz-1
            wf_b_all(0)%p => wf(:, ib)
            wf_b_all(1)%p => wfdr(:, ib)
            wf_b_all(2)%p => wfdp(:, ib) ! nl(ib)*y(:)*wf(:,ib)
            wf_b_all(3)%p => wfdz(:, ib)
            wf_b_all(4)%p => wfd2_all(:, ib)
            
            do ixx = -1, 1, 2
               do iyy = 0, 4
                  do izz = 0, max_ind(iyy)
                     coeff = merge(1.0_dp,0.0_dp,izz>0) ! Equivalent to "(izz>0) ? 1.0 : 0.0;" in C
                     ! hpsi(iyy,ixx)%re(:,iz) = hpsi(iyy,ixx)%re(:,iz) + mf(iyy,izz,ixx,ns(ib))%re(:) * wf_b_all(izz)%p(:)
                     call add_mul_1d(mf(iyy,izz,ixx,ns(ib))%re(:), wf_b_all(izz)%p(:), &
                                     coeff, hpsi(iyy,ixx)%re(:,iz))
                     ! hpsi(iyy,ixx)%im(:,iz) = hpsi(iyy,ixx)%im(:,iz) + mf(iyy,izz,ixx,ns(ib))%im(:) * wf_b_all(izz)%p(:)
                     call add_mul_1d(mf(iyy,izz,ixx,ns(ib))%im(:), wf_b_all(izz)%p(:), &
                                     coeff, hpsi(iyy,ixx)%im(:,iz))
                  end do
               end do
            end do
         end do ! iz, sp state in one block
         !!$OMP END PARALLEL DO
         
         ! Perform matrix-matrix multiplication to obtain reh and imh
         ! Use OpenMP-enabled BLAS library to achieve best performance
         ia = isstart(ix)

         ! wf_a has spin up
         if (num_spin_up(ix)>0) then
            wf_a_start(0)%p => wf(1, ia)
            wf_a_start(1)%p => wfdr(1, ia)
            wf_a_start(2)%p => wfdp(1, ia)
            wf_a_start(3)%p => wfdz(1, ia)
            wf_a_start(4)%p => wfd2_all(1, ia)
            do iz = 0, 4
               coeff = merge(1.0_dp,0.0_dp,iz>0) ! Equivalent to "(iz>0) ? 1.0 : 0.0;" in C
               call dgemm('T','N',num_spin_up(ix),db(iy),nghl,2.0_dp,wf_a_start(iz)%p,nghl,&
                          hpsi(iz,1)%re(1,1),nghl,coeff,reh%elem(reh%ir2m(ix)),db(ix))
               call dgemm('T','N',num_spin_up(ix),db(iy),nghl,2.0_dp,wf_a_start(iz)%p,nghl,&
                          hpsi(iz,1)%im(1,1),nghl,coeff,imh%elem(imh%ir2m(ix)),db(ix))
            end do
         end if

         ! wf_a has spin down
         if (num_spin_dw(ix)>0) then
            ia = ia + num_spin_up(ix)
            wf_a_start(0)%p => wf(1, ia)
            wf_a_start(1)%p => wfdr(1, ia)
            wf_a_start(2)%p => wfdp(1, ia)
            wf_a_start(3)%p => wfdz(1, ia)
            wf_a_start(4)%p => wfd2_all(1, ia)
            do iz = 0, 4
               coeff = merge(1.0_dp,0.0_dp,iz>0) ! Equivalent to "(iz>0) ? 1.0 : 0.0;" in C
               call dgemm('T','N',num_spin_dw(ix),db(iy),nghl,2.0_dp,wf_a_start(iz)%p,nghl,&
                          hpsi(iz,-1)%re(1,1),nghl,coeff,reh%elem(reh%ir2m(ix)+num_spin_up(ix)),db(ix))
               call dgemm('T','N',num_spin_dw(ix),db(iy),nghl,2.0_dp,wf_a_start(iz)%p,nghl,&
                          hpsi(iz,-1)%im(1,1),nghl,coeff,imh%elem(imh%ir2m(ix)+num_spin_up(ix)),db(ix))
            end do
         end if
      end do ! ix, block index
      
   end subroutine meanfield

   !---------------------------------------------------------------------------
   ! Compute the pairing potential Delta from the pairing densities
   ! rho-breve(r) and s-breve(r)
   !---------------------------------------------------------------------------
   subroutine pairingfield(ds, red, imd)
      use pnfam_interaction, only : cpair, cspair

      implicit none
      type(density_set), intent(in) :: ds
      type(blockmatrix), intent(inout) :: red, imd

      integer :: ia, ib, ix, iy, iz

      type(field), save :: aux_pp, aux_pm, aux_mp, aux_mm

      ! Results of pairing field applied on wave functions
      ! index gives the spin of wf_a
      type(array_2d), save :: dpsi(-1:1)

      if (.not. allocated(aux_pp%re)) then
         call allocate_field(aux_pp, nghl)
         call allocate_field(aux_pm, nghl)
         call allocate_field(aux_mp, nghl)
         call allocate_field(aux_mm, nghl)
         do iz = -1, 1, 2
            call allocate_2darray(dpsi(iz), nghl, maxval_db)
         end do
      end if

      ! red%elem = 0 ; imd%elem = 0 ! Zeroing Delta is unnecessary 

      ! Prepare pairing field
      aux_pp%re = cspair*(ds%resbr + ds%imsbp) ! |a> = |+>, |b> = |+>
      aux_pp%im = cspair*(ds%imsbr - ds%resbp) ! |a> = |+>, |b> = |+>
      aux_mp%re = -cpair*ds%rerb - cspair*ds%resbz ! |a> = |->, |b> = |+>
      aux_mp%im = -cpair*ds%imrb - cspair*ds%imsbz ! |a> = |->, |b> = |+>
      aux_pm%re = cpair*ds%rerb - cspair*ds%resbz ! |a> = |+>, |b> = |->
      aux_pm%im = cpair*ds%imrb - cspair*ds%imsbz ! |a> = |+>, |b> = |->
      aux_mm%re = cspair*(-ds%resbr + ds%imsbp) ! |a> = |->, |b> = |->
      aux_mm%im = cspair*(-ds%imsbr - ds%resbp) ! |a> = |->, |b> = |->

      ! Loop over blocks
      do ix = 1, nb
         iy = red%ir2c(ix) ! column index of current block
         if (iy==0) cycle

         ! Apply pairing field onto wave functions
         !!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ib)
         do iz = 1, db(iy)
            ! Coulomb vector |b>
            ib = isstart(iy)+iz-1

            if (ns(ib)==+1) then
               ! |a> = |+>, |b> = |+>
               dpsi(1)%re(:, iz) = aux_pp%re(:)*wf(:, ib)
               dpsi(1)%im(:, iz) = aux_pp%im(:)*wf(:, ib)

               ! |a> = |->, |b> = |+>
               dpsi(-1)%re(:, iz) = aux_mp%re(:)*wf(:, ib)
               dpsi(-1)%im(:, iz) = aux_mp%im(:)*wf(:, ib)

            else if (ns(ib)==-1) then
               ! |a> = |+>, |b> = |->
               dpsi(1)%re(:, iz) = aux_pm%re(:)*wf(:, ib)
               dpsi(1)%im(:, iz) = aux_pm%im(:)*wf(:, ib)

               ! |a> = |->, |b> = |->
               dpsi(-1)%re(:, iz) = aux_mm%re(:)*wf(:, ib)
               dpsi(-1)%im(:, iz) = aux_mm%im(:)*wf(:, ib)

            ! else
            !    call abort(' Inconsistent spin in helpers.pairingfield')
            end if
         end do

         ! Perform matrix-matrix multiplication to obtain red and imd
         ! Use OpenMP-enabled BLAS library to achieve best performance
         ia = isstart(ix)

         ! wf_a has spin up
         if (num_spin_up(ix)>0) then
            call dgemm('T','N',num_spin_up(ix),db(iy),nghl,2.0_dp,wf(1,ia),nghl,&
                       dpsi(1)%re(1,1),nghl,0.0_dp,red%elem(red%ir2m(ix)),db(ix))
            call dgemm('T','N',num_spin_up(ix),db(iy),nghl,2.0_dp,wf(1,ia),nghl,&
                       dpsi(1)%im(1,1),nghl,0.0_dp,imd%elem(imd%ir2m(ix)),db(ix))
         end if

         ! wf_a has spin down
         if (num_spin_dw(ix)>0) then
            ia = ia+num_spin_up(ix)
            call dgemm('T','N',num_spin_dw(ix),db(iy),nghl,2.0_dp,wf(1,ia),nghl,&
                       dpsi(-1)%re(1,1),nghl,0.0_dp,red%elem(red%ir2m(ix)+num_spin_up(ix)),db(ix))
            call dgemm('T','N',num_spin_dw(ix),db(iy),nghl,2.0_dp,wf(1,ia),nghl,&
                       dpsi(-1)%im(1,1),nghl,0.0_dp,imd%elem(imd%ir2m(ix)+num_spin_up(ix)),db(ix))
         end if
      end do ! ix, block index

   end subroutine pairingfield

end module pnfam_hamiltonian
