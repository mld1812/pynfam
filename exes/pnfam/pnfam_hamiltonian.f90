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
!>   Grouping similar array operations, performing only one loop.
!>
!> @authors M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
!> @authors E.M. Ney, UNC Chapel Hill, 2018-
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

   ! Coordinate space local densities, custom type for organization
   type(density_set) :: rhox

   ! The following arrays allow collapsing the do loops over blocks and states
   ! within blocks into one;  This is necessary for an effective OpenMP
   ! parallelization
   integer, allocatable :: ia_tab(:), ib_tab(:), ia_tabp(:), ib_tabp(:)
   logical, save :: module_initialized = .false.

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
   ! Initialization of ia_tab and ib_tab; Need to call only once
   !---------------------------------------------------------------------------
   subroutine init_hamiltonian_module(rerho, rek)
      implicit none
      type(blockmatrix), intent(in) :: rerho, rek
      integer :: n, ibx1, ibx2, nd1, nd2, i1a, i1b, i2a, i2b, ia, ib, im
      
      if (.not.allocated(ia_tab)) then
         n = size(rerho%elem)
         allocate(ia_tab(n), ib_tab(n))
         n = size(rek%elem)
         allocate(ia_tabp(n), ib_tabp(n))
      end if
      
      im = 0
      do ibx1 = 1,nb
         ibx2 = rerho%ir2c(ibx1)
         if (ibx2 == 0) cycle
         
         nd1 = db(ibx1) ; nd2 = db(ibx2)
         if (nd1 == 0 .or. nd2 == 0) cycle
         i1a = isstart(ibx1) ; i1b = i1a + nd1 - 1
         i2a = isstart(ibx2) ; i2b = i2a + nd2 - 1
         
         do ib = i2a, i2b
            do ia = i1a, i1b
               im = im + 1
               ia_tab(im) = ia
               ib_tab(im) = ib
            end do
         end do
         
      end do
      
      ! Same for the pairing tensor
      im = 0
      do ibx1 = 1, nb

         ! Determine the range of particle indices in this block
         ibx2 = rek%ir2c(ibx1)
         if (ibx2 == 0) cycle
         
         nd1 = db(ibx1) ; nd2 = db(ibx2)
         if (nd1 == 0 .or. nd2 == 0) cycle
         i1a = isstart(ibx1) ; i1b = i1a + nd1 - 1
         i2a = isstart(ibx2) ; i2b = i2a + nd2 - 1

         ! Column vector |b>
         do ib = i2a, i2b
            ! Row vector <a|
            do ia = i1a, i1b
               im = im + 1
               ia_tabp(im) = ia
               ib_tabp(im) = ib
            end do
         end do
      end do
      
      ! Inelegantly, we need to re-call this routine every time, because the
      ! pn and np sectors have a different matrix structure
      !module_initialized = .true.
      
   end subroutine

   !---------------------------------------------------------------------------
   ! Computation of the local densities in coordinate space
   !---------------------------------------------------------------------------
   subroutine density(rerho, imrho, rek, imk, ds)

      implicit none

      type(blockmatrix), intent(inout) :: rerho, imrho, rek, imk
      type(density_set), intent(inout) :: ds

      integer  :: ia, ib, im, xla, xlb, nqp
      real(dp) :: rerho_ab, imrho_ab
      real(dp), dimension(nghl) :: wf_a, wf_b, dr_wf_a, dr_wf_b, dz_wf_a, dz_wf_b, wf_wf, aux
      
      ! The following are needed to allow OpenMP reduction, which does not currently
      ! work for derived types.  Hacky :(
      real(dp), dimension(:), allocatable :: dsrerho, dsimrho, dsretau, dsimtau,         &
         dsretjrr, dsimtjrr, dsretjpr, dsimtjpr, dsretjzr, dsimtjzr, dsretjrp, dsimtjrp, &
         dsretjpp, dsimtjpp, dsretjzp, dsimtjzp, dsretjrz, dsimtjrz, dsretjpz, dsimtjpz, &
         dsretjzz, dsimtjzz, dsresr, dsimsr, dsresp, dsimsp, dsresz, dsimsz, dsretr,     &
         dsimtr, dsretp, dsimtp, dsretz, dsimtz, dsrejr, dsimjr, dsrejp, dsimjp,         &
         dsrejz, dsimjz, dsrefr, dsimfr, dsrefp, dsimfp, dsrefz, dsimfz, dsregs,         &
         dsimgs, dsrerb, dsimrb, dsresbr, dsimsbr, dsresbp, dsimsbp, dsresbz, dsimsbz
      
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
      end if

      ! Initialize the densities to zero
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
      if (.not.module_initialized) call init_hamiltonian_module(rerho,rek)

      ! Loop over the blocks.
      ! Calculate e.g. rho(x) = Sum_{ab} rho_ab wf_a(x) Conjg[wf_b(x)]
      
      !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(nl,ns,y,wf,wfdr,wfdz,rerho,imrho,ia_tab,ib_tab) &
      !$OMP& REDUCTION(+:dsrerho) REDUCTION(+:dsimrho) REDUCTION(+:dsretau) REDUCTION(+:dsimtau) &
      !$OMP& REDUCTION(+:dsretjrr) REDUCTION(+:dsimtjrr) REDUCTION(+:dsretjpr) REDUCTION(+:dsimtjpr) &
      !$OMP& REDUCTION(+:dsretjzr) REDUCTION(+:dsimtjzr) REDUCTION(+:dsretjrp) REDUCTION(+:dsimtjrp) &
      !$OMP& REDUCTION(+:dsretjpp) REDUCTION(+:dsimtjpp) REDUCTION(+:dsretjzp) REDUCTION(+:dsimtjzp) &
      !$OMP& REDUCTION(+:dsretjrz) REDUCTION(+:dsimtjrz) REDUCTION(+:dsretjpz) REDUCTION(+:dsimtjpz) &
      !$OMP& REDUCTION(+:dsretjzz) REDUCTION(+:dsimtjzz) REDUCTION(+:dsresr) REDUCTION(+:dsimsr) &
      !$OMP& REDUCTION(+:dsresp) REDUCTION(+:dsimsp) REDUCTION(+:dsresz) REDUCTION(+:dsimsz) &
      !$OMP& REDUCTION(+:dsretr) REDUCTION(+:dsimtr) REDUCTION(+:dsretp) REDUCTION(+:dsimtp) &
      !$OMP& REDUCTION(+:dsretz) REDUCTION(+:dsimtz) REDUCTION(+:dsrejr) REDUCTION(+:dsimjr) &
      !$OMP& REDUCTION(+:dsrejp) REDUCTION(+:dsimjp) REDUCTION(+:dsrejz) REDUCTION(+:dsimjz) &
      !$OMP& REDUCTION(+:dsrefr) REDUCTION(+:dsimfr) REDUCTION(+:dsrefp) REDUCTION(+:dsimfp) &
      !$OMP& REDUCTION(+:dsrefz) REDUCTION(+:dsimfz) REDUCTION(+:dsregs) REDUCTION(+:dsimgs)
      do im = 1, size(rerho%elem)
         ia = ia_tab(im)
         ib = ib_tab(im)
         
            ! Column vector |b>
            wf_b = wf(:, ib)
            dr_wf_b = wfdr(:, ib)
            dz_wf_b = wfdz(:, ib)

            ! Angular quantum number Lambda = Omega - s_z
            xlb = nl(ib)

               ! Row vector <a|
               wf_a = wf(:, ia)
               dr_wf_a = wfdr(:, ia)
               dz_wf_a = wfdz(:, ia)

               ! Angular quantum number Lambda = Omega - s_z
               xla = nl(ia)

               ! Pick the matrix element rho_ab
               rerho_ab = rerho%elem(im)
               imrho_ab = imrho%elem(im)
               
               ! Auxiliary array
               wf_wf(:) = wf_a(:)*wf_b(:)

               ! Diagonal in spin
               if (ns(ia) == ns(ib)) then
               
                  ! rho: scalar density
                  dsrerho(:) = dsrerho(:) + rerho_ab*wf_wf(:)
                  dsimrho(:) = dsimrho(:) + imrho_ab*wf_wf(:)

                  ! s_z: z-component of vector density
                  dsresz(:) = dsresz(:) + (ns(ia)*rerho_ab)*wf_wf(:)
                  dsimsz(:) = dsimsz(:) + (ns(ia)*imrho_ab)*wf_wf(:)
                  
                  ! j_phi: phi-component of current density
                  aux(:) = y(:)*wf_wf(:)
                  dsrejp(:) = dsrejp(:) + ((xla+xlb)*rerho_ab/2)*aux(:)
                  dsimjp(:) = dsimjp(:) + ((xla+xlb)*imrho_ab/2)*aux(:)
                  
                  ! J_{phi,z}: (phi,z)-component of tensor density [= J_{2,3}/rho]
                  !aux(:) = y(:)*wf_wf(:)
                  dsretjpz(:) = dsretjpz(:) + ((xla+xlb)*ns(ia)*rerho_ab/2)*aux(:)
                  dsimtjpz(:) = dsimtjpz(:) + ((xla+xlb)*ns(ia)*imrho_ab/2)*aux(:)

                  
                  ! tau: kinetic density
                  !aux(:) = dr_wf_a(:)*dr_wf_b(:) + dz_wf_a(:)*dz_wf_b(:)  &
                  !       + xla*xlb*y(:)**2*wf_wf(:)
                  aux(:) = dr_wf_a(:)*dr_wf_b(:) + dz_wf_a(:)*dz_wf_b(:)  &
                         + (xla*xlb)*y(:)*aux(:)
                  dsretau(:) = dsretau(:) + rerho_ab*aux(:)
                  dsimtau(:) = dsimtau(:) + imrho_ab*aux(:)
                  
                  ! T_z: z-component of spin-kinetic density
                  !aux(:) = dr_wf_a(:)*dr_wf_b(:) + dz_wf_a(:)*dz_wf_b(:)  &
                  !       + (xla*xlb)*y(:)**2*wf_wf(:)
                  dsretz(:) = dsretz(:) + (ns(ia)*rerho_ab)*aux(:)
                  dsimtz(:) = dsimtz(:) + (ns(ia)*imrho_ab)*aux(:)
                  
                  
                  ! j_rho: rho-component of current density
                  aux(:) = wf_a(:)*dr_wf_b(:) - dr_wf_a(:)*wf_b(:)
                  dsrejr(:) = dsrejr(:) - imrho_ab/2*aux(:)
                  dsimjr(:) = dsimjr(:) + rerho_ab/2*aux(:)
                  
                  ! J_{rho,z}: (rho,z)-component of tensor density
                  !aux(:) = wf_a(:)*dr_wf_b(:) - dr_wf_a(:)*wf_b(:)
                  dsretjrz(:) = dsretjrz(:) - (ns(ia)*imrho_ab/2)*aux(:)
                  dsimtjrz(:) = dsimtjrz(:) + (ns(ia)*rerho_ab/2)*aux(:)
                  
                  
                  ! j_z: z-component of current density
                  aux(:) = wf_a(:)*dz_wf_b(:) - dz_wf_a(:)*wf_b(:)
                  dsrejz(:) = dsrejz(:) - (imrho_ab/2)*aux(:)
                  dsimjz(:) = dsimjz(:) + (rerho_ab/2)*aux(:)

                  ! J_{z,z}: (z,z)-component of tensor density
                  !aux(:) = wf_a(:)*dz_wf_b(:) - dz_wf_a(:)*wf_b(:)
                  dsretjzz(:) = dsretjzz(:) - (ns(ia)*imrho_ab/2)*aux(:)
                  dsimtjzz(:) = dsimtjzz(:) + (ns(ia)*rerho_ab/2)*aux(:)
                  
                  ! Del.s: z-component of divergence of s
                  aux(:) = dz_wf_a(:)*wf_b(:) + wf_a(:)*dz_wf_b(:)
                  dsregs(:) = dsregs(:) + (ns(ia)*rerho_ab)*aux(:)
                  dsimgs(:) = dsimgs(:) + (ns(ia)*imrho_ab)*aux(:)


                  ! F_rho: rho-component of tensor-kinetic density
                  aux(:) = dr_wf_a(:)*dz_wf_b(:) + dz_wf_a(:)*dr_wf_b(:)
                  dsrefr(:) = dsrefr(:) + (ns(ia)*rerho_ab/2)*aux(:)
                  dsimfr(:) = dsimfr(:) + (ns(ia)*imrho_ab/2)*aux(:)

                  ! F_phi: phi-component of tensor-kinetic density
                  aux(:) = y(:)*(xla*wf_a(:)*dz_wf_b(:) - xlb*wf_b(:)*dz_wf_a(:))
                  dsrefp(:) = dsrefp(:) - (ns(ia)*imrho_ab/2)*aux(:)
                  dsimfp(:) = dsimfp(:) + (ns(ia)*rerho_ab/2)*aux(:)

                  ! F_z: z-component of tensor-kinetic density
                  aux(:) = dz_wf_a(:)*dz_wf_b(:)
                  dsrefz(:) = dsrefz(:) + (ns(ia)*rerho_ab)*aux(:)
                  dsimfz(:) = dsimfz(:) + (ns(ia)*imrho_ab)*aux(:)


               ! Not-diagonal in spin
               else

                  ! |a> = |+>,  |b> = |->
                  if (ns(ia) == 1 .and. ns(ib) == -1) then

                     ! s_rho: rho-component of vector density
                     dsresr(:) = dsresr(:) + rerho_ab*wf_wf(:)
                     dsimsr(:) = dsimsr(:) + imrho_ab*wf_wf(:)

                     ! s_phi: phi-component of vector density
                     dsresp(:) = dsresp(:) - imrho_ab*wf_wf(:)
                     dsimsp(:) = dsimsp(:) + rerho_ab*wf_wf(:)
                     
                     ! J_{phi,rho}: (phi,rho)-component of tensor density [= J_{2,1}/rho]
                     aux(:) = y(:)*wf_wf(:)
                     dsretjpr(:) = dsretjpr(:) + ((xla+xlb)*rerho_ab/2)*aux(:)
                     dsimtjpr(:) = dsimtjpr(:) + ((xla+xlb)*imrho_ab/2)*aux(:)

                     ! J_{phi,phi}: (phi,phi)-component of tensor density [= J_{2,2}/rho^2]
                     !aux(:) = y(:)*wf_wf(:)
                     dsretjpp(:) = dsretjpp(:) - (0.5_dp*imrho_ab*(xla+xlb))*aux(:)
                     dsimtjpp(:) = dsimtjpp(:) + (0.5_dp*rerho_ab*(xla+xlb))*aux(:)
                     
                     ! Del.s: phi-component of divergence of s
                     !aux(:) = y(:)*wf_wf(:)
                     dsregs(:) = dsregs(:) + (rerho_ab*(xlb-xla-1))*aux(:)
                     dsimgs(:) = dsimgs(:) + (imrho_ab*(xlb-xla-1))*aux(:)
                     
                     
                     ! T_rho: rho-component of spin-kinetic density
                     !aux(:) = dr_wf_a(:)*dr_wf_b(:) + dz_wf_a(:)*dz_wf_b(:)  &
                     !         + xla*xlb*y(:)*y(:)*wf_wf(:)
                     aux(:) = dr_wf_a(:)*dr_wf_b(:) + dz_wf_a(:)*dz_wf_b(:)  &
                              + xla*xlb*y(:)*aux(:)
                     dsretr(:) = dsretr(:) + rerho_ab*aux(:)
                     dsimtr(:) = dsimtr(:) + imrho_ab*aux(:)

                     ! T_phi: phi-component of spin-kinetic density
                     !aux(:) = dr_wf_a(:)*dr_wf_b(:) + dz_wf_a(:)*dz_wf_b(:) &
                     !         + xla*xlb*y(:)*y(:)*wf_wf(:)
                     dsretp(:) = dsretp(:) - imrho_ab*aux(:)
                     dsimtp(:) = dsimtp(:) + rerho_ab*aux(:)


                     ! J_{rho,rho}: (rho,rho)-component of tensor density
                     aux(:) = wf_a(:)*dr_wf_b(:) - dr_wf_a(:)*wf_b(:)
                     dsretjrr(:) = dsretjrr(:) - imrho_ab/2*aux(:)
                     dsimtjrr(:) = dsimtjrr(:) + rerho_ab/2*aux(:)

                     ! J_{rho,phi}: (rho,phi)-component of tensor density [= J_{1,2}/rho]
                     !aux(:) = wf_a(:)*dr_wf_b(:) - dr_wf_a(:)*wf_b(:)
                     dsretjrp(:) = dsretjrp(:) - rerho_ab/2*aux(:)
                     dsimtjrp(:) = dsimtjrp(:) - imrho_ab/2*aux(:)
                     
                     
                     ! Del.s: rho-component of divergence of s
                     aux(:) = dr_wf_a(:)*wf_b(:) + wf_a(:)*dr_wf_b(:) + y(:)*wf_wf(:)
                     dsregs(:) = dsregs(:) + rerho_ab*aux(:)
                     dsimgs(:) = dsimgs(:) + imrho_ab*aux(:)
                     
                     
                     ! J_{z,rho}: (z,rho)-component of tensor density
                     aux(:) = wf_a(:)*dz_wf_b(:) - dz_wf_a(:)*wf_b(:)
                     dsretjzr(:) = dsretjzr(:) - imrho_ab/2*aux(:)
                     dsimtjzr(:) = dsimtjzr(:) + rerho_ab/2*aux(:)

                     ! J_{z,phi}: (z,phi)-component of tensor density [= J_{3,2}/rho]
                     !aux(:) = wf_a(:)*dz_wf_b(:) - dz_wf_a(:)*wf_b(:)
                     dsretjzp(:) = dsretjzp(:) - (0.5_dp*rerho_ab)*aux(:)
                     dsimtjzp(:) = dsimtjzp(:) - (0.5_dp*imrho_ab)*aux(:)


                     ! F_rho: rho-component of tensor-kinetic density
                     aux(:) = dr_wf_a(:)*dr_wf_b(:) + 0.5_dp*y(:)*(xlb*dr_wf_a(:)*wf_b(:)  &
                                             - xla*wf_a(:)*dr_wf_b(:))
                     dsrefr(:) = dsrefr(:) + rerho_ab*aux(:)
                     dsimfr(:) = dsimfr(:) + imrho_ab*aux(:)

                     ! F_phi: phi-component of tensor-kinetic density
                     aux(:) = 0.5_dp*y(:)*(xla*wf_a(:)*dr_wf_b(:) - xlb*dr_wf_a(:)*wf_b(:))  &
                                             + xla*xlb*y(:)**2*wf_wf(:)
                     dsrefp(:) = dsrefp(:) - imrho_ab*aux(:)
                     dsimfp(:) = dsimfp(:) + rerho_ab*aux(:)

                     ! F_z: z-component of tensor-kinetic density
                     aux(:) = dz_wf_a(:)*dr_wf_b(:) + dr_wf_a(:)*dz_wf_b(:)  &
                            + y(:)*(xlb*dz_wf_a(:)*wf_b(:) - xla*wf_a(:)*dz_wf_b(:))
                     dsrefz(:) = dsrefz(:) + (0.5_dp*rerho_ab)*aux(:)
                     dsimfz(:) = dsimfz(:) + (0.5_dp*imrho_ab)*aux(:)

                  ! |a> = |->,  |b> = |+>
                  else if (ns(ia) == -1 .and. ns(ib) == 1) then

                     ! s_rho: rho-component of vector density
                     dsresr(:) = dsresr(:) + rerho_ab*wf_wf(:)
                     dsimsr(:) = dsimsr(:) + imrho_ab*wf_wf(:)

                     ! s_phi: phi-component of vector density
                     dsresp(:) = dsresp(:) + imrho_ab*wf_wf(:)
                     dsimsp(:) = dsimsp(:) - rerho_ab*wf_wf(:)

                     ! J_{phi,rho}: (phi,rho)-component of tensor density [= J_{2,1}/rho]
                     aux(:) = y(:)*wf_wf(:)
                     dsretjpr(:) = dsretjpr(:) + (0.5_dp*rerho_ab*(xla+xlb))*aux(:)
                     dsimtjpr(:) = dsimtjpr(:) + (0.5_dp*imrho_ab*(xla+xlb))*aux(:)

                     ! J_{phi,phi}: (phi,phi)-component of tensor density [= J_{2,2}/rho]
                     !aux(:) = y(:)*wf_wf(:)
                     dsretjpp(:) = dsretjpp(:) + (0.5_dp*imrho_ab*(xla+xlb))*aux(:)
                     dsimtjpp(:) = dsimtjpp(:) - (0.5_dp*rerho_ab*(xla+xlb))*aux(:)
                     
                     ! Del.s: phi-component of divergence of s
                     !aux(:) = y(:)*wf_wf(:)
                     dsregs(:) = dsregs(:) - (rerho_ab*(xlb-xla+1))*aux(:)
                     dsimgs(:) = dsimgs(:) - (imrho_ab*(xlb-xla+1))*aux(:)
                     
                     
                     ! T_rho: rho-component of spin-kinetic density
                     !aux(:) = dr_wf_a(:)*dr_wf_b(:) + dz_wf_a(:)*dz_wf_b(:)   &
                     !         + xla*xlb*y(:)**2*wf_wf(:)
                     aux(:) = dr_wf_a(:)*dr_wf_b(:) + dz_wf_a(:)*dz_wf_b(:)   &
                              + xla*xlb*y(:)*aux(:)
                     dsretr(:) = dsretr(:) + rerho_ab*aux(:)
                     dsimtr(:) = dsimtr(:) + imrho_ab*aux(:)

                     ! T_phi: phi-component of spin-kinetic density
                     !aux(:) = dr_wf_a(:)*dr_wf_b(:) + dz_wf_a(:)*dz_wf_b(:) &
                     !         + xla*xlb*y(:)**2*wf_wf(:)
                     dsretp(:) = dsretp(:) + imrho_ab*aux(:)
                     dsimtp(:) = dsimtp(:) - rerho_ab*aux(:)
                     
                     
                     ! J_{rho,rho}: (rho,rho)-component of tensor density
                     aux(:) = wf_a(:)*dr_wf_b(:) - dr_wf_a(:)*wf_b(:)
                     dsretjrr(:) = dsretjrr(:) - (0.5_dp*imrho_ab)*aux(:)
                     dsimtjrr(:) = dsimtjrr(:) + (0.5_dp*rerho_ab)*aux(:)

                     ! J_{rho,phi}: (rho,phi)-component of tensor density [= J_{1,2}/rho]
                     !aux(:) = wf_a(:)*dr_wf_b(:) - dr_wf_a(:)*wf_b(:)
                     dsretjrp(:) = dsretjrp(:) + (0.5_dp*rerho_ab)*aux(:)
                     dsimtjrp(:) = dsimtjrp(:) + (0.5_dp*imrho_ab)*aux(:)


                     ! Del.s: rho-component of divergence of s
                     aux(:) = dr_wf_a(:)*wf_b(:) + wf_a(:)*dr_wf_b(:) + y(:)*wf_wf(:)
                     dsregs(:) = dsregs(:) + rerho_ab*aux(:)
                     dsimgs(:) = dsimgs(:) + imrho_ab*aux(:)
                     
                  
                     ! J_{z,rho}: (z,rho)-component of tensor density
                     aux(:) = wf_a(:)*dz_wf_b(:) - dz_wf_a(:)*wf_b(:)
                     dsretjzr(:) = dsretjzr(:) - (0.5_dp*imrho_ab)*aux(:)
                     dsimtjzr(:) = dsimtjzr(:) + (0.5_dp*rerho_ab)*aux(:)

                     ! J_{z,phi}: (z,phi)-component of tensor density [= J_{3,2}/rho]
                     !aux(:) = wf_a(:)*dz_wf_b(:) - dz_wf_a(:)*wf_b(:)
                     dsretjzp(:) = dsretjzp(:) + (0.5_dp*rerho_ab)*aux(:)
                     dsimtjzp(:) = dsimtjzp(:) + (0.5_dp*imrho_ab)*aux(:)


                     ! F_rho: rho-component of tensor-kinetic density
                     aux(:) = dr_wf_a(:)*dr_wf_b(:) - 0.5_dp*y(:)*(xlb*dr_wf_a(:)*wf_b(:)  &
                                             - xla*wf_a(:)*dr_wf_b(:))
                     dsrefr(:) = dsrefr(:) + rerho_ab*aux(:)
                     dsimfr(:) = dsimfr(:) + imrho_ab*aux(:)

                     ! F_phi: phi-component of tensor-kinetic density
                     aux(:) = 0.5_dp*y(:)*(xla*wf_a(:)*dr_wf_b(:) - xlb*dr_wf_a(:)*wf_b(:)) &
                                             - y(:)*y(:)*xla*xlb*wf_wf(:)
                     dsrefp(:) = dsrefp(:) - imrho_ab*aux(:)
                     dsimfp(:) = dsimfp(:) + rerho_ab*aux(:)

                     ! F_z: z-component of tensor-kinetic density
                     aux(:) = dz_wf_a(:)*dr_wf_b(:) + dr_wf_a(:)*dz_wf_b(:)   &
                              - y(:)*(xlb*dz_wf_a(:)*wf_b(:) - xla*wf_a(:)*dz_wf_b(:))
                     dsrefz(:) = dsrefz(:) + (0.5_dp*rerho_ab)*aux(:)
                     dsimfz(:) = dsimfz(:) + (0.5_dp*imrho_ab)*aux(:)

                  else
                     call abort(' Inconsistent spin in helpers.density')
                  end if
               end if

            !end do   ! ia
         !end do      ! ib
      end do
      !$OMP END PARALLEL DO
      
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
      nqp = sum(db)
      
      !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(nl,ns,y,wf,wfdr,wfdz,rek,imk,ia_tabp,ib_tabp)    &
      !$OMP& REDUCTION(+:dsrerb)  REDUCTION(+:dsimrb)  REDUCTION(+:dsresbr) REDUCTION(+:dsimsbr) &
      !$OMP& REDUCTION(+:dsresbp) REDUCTION(+:dsimsbp) REDUCTION(+:dsresbz) REDUCTION(+:dsimsbz)
      do im = 1, size(rek%elem)
         ia = ia_tabp(im)
         ib = ib_tabp(im)

         ! Column vector |b>
         wf_b = wf(:, ib)

         ! Row vector <a|
         wf_a = wf(:, ia)
         
         wf_wf(:) = wf_a(:)*wf_b(:)
         rerho_ab = rek%elem(im)
         imrho_ab = imk%elem(im)

         ! Diagonal in spin
         if (ns(ia) == -ns(ib)) then  ! note the sign flip w.r.t. normal densities
            ! rho-breve_10
            dsrerb = dsrerb + 2*rerho_ab*ns(ia)*wf_wf
            dsimrb = dsimrb + 2*imrho_ab*ns(ia)*wf_wf
            
            ! s-breve_0,z
            dsresbz = dsresbz + 2*rerho_ab*wf_wf
            dsimsbz = dsimsbz + 2*imrho_ab*wf_wf
         ! Not-diagonal in spin
         else
            ! |a> = |+>,  |b> = |->
            if (ns(ia) == -1 .and. ns(ib) == -1) then  ! note the sign flip w.r.t. normal densities
               ! s-breve_0,rho
               dsresbr = dsresbr + 2*rerho_ab*wf_wf
               dsimsbr = dsimsbr + 2*imrho_ab*wf_wf
               
               ! s-breve_0,phi
               dsresbp = dsresbp + 2*imrho_ab*wf_wf
               dsimsbp = dsimsbp - 2*rerho_ab*wf_wf
            else
               ! s-breve_0,rho
               dsresbr = dsresbr - 2*rerho_ab*wf_wf
               dsimsbr = dsimsbr - 2*imrho_ab*wf_wf
               
               ! s-breve_0,phi
               dsresbp = dsresbp + 2*imrho_ab*wf_wf
               dsimsbp = dsimsbp - 2*rerho_ab*wf_wf
            end if
         end if
         
      end do
      !$OMP END PARALLEL DO

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

      integer  :: ia, ib, xla, xlb, ix
      real(dp), dimension(nghl) :: wf_a, wf_b, dr_wf_a, dr_wf_b, dz_wf_a, dz_wf_b, d2_wf_a, d2_wf_b

      ! Temporary arrays for optimized tensor integrations
      real(dp), dimension(nghl) :: rethrr, imthrp, rethrz, imthpr, rethpp
      real(dp), dimension(nghl) :: imthpz, rethzr, imthzp, rethzz
      
      real(dp), dimension(nghl) :: aux, wf_wf
      real(dp) :: re, im

      ! Zero out the mean field
      reh%elem = 0; imh%elem = 0

      ! Loop over the blocks.
      ! Calculate e.g. h_ab = Sum_{sa,sb} Int{ Conjg[wf_a(x,sa)] [ f(x;sa,sb) ] wf_b(x,sb) }
      
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ia,ib,wf_a,dr_wf_a,dz_wf_a,d2_wf_a) &
      !$OMP& PRIVATE(wf_b,dr_wf_b,dz_wf_b,d2_wf_b,xla,xlb,re,im,aux,wf_wf) &
      !$OMP& PRIVATE(rethrr,imthrp,rethrz,imthpr,rethpp,imthpz,rethzr,imthzp,rethzz)
      do ix = 1, size(reh%elem)
         ia = ia_tab(ix)
         ib = ib_tab(ix)
         
            wf_b = wf(:, ib)
            dr_wf_b = wfdr(:, ib)
            dz_wf_b = wfdz(:, ib)
            d2_wf_b = wfd2(:, ib)

            ! Angular quantum number Lambda = Omega - s_z
            xlb = nl(ib)

            ! Row vector <a|
               wf_a = wf(:, ia)
               dr_wf_a = wfdr(:, ia)
               dz_wf_a = wfdz(:, ia)
               d2_wf_a = wfd2(:, ia)
               
               wf_wf(:) = wf_a(:)*wf_b(:)

               ! Angular quantum number Lambda = Omega - s_z
               xla = nl(ia)
               
               ! To avoid repeated array lookups that the compiler might not know
               ! to optimize away, use temporary variables re and im:
               re = 0 ; im = 0

               ! Diagonal in spin
               if (ns(ia) == ns(ib)) then

                  ! rho-rho
                  aux(:) = crho(:)*wf_wf(:)
                  re = re + 4.0_dp*sum(ds%rerho(:)*aux(:))
                  im = im + 4.0_dp*sum(ds%imrho(:)*aux(:))
                  
                  ! z-component of s-s
                  aux(:) = cs(:)*wf_wf(:)
                  re = re + 4.0_dp*ns(ia)*sum(ds%resz(:)*aux(:))
                  im = im + 4.0_dp*ns(ia)*sum(ds%imsz(:)*aux(:))

                  ! z-component of s-T
                  re = re + 2.0_dp*ns(ia)*ct*sum(ds%retz(:)*wf_wf(:))
                  im = im + 2.0_dp*ns(ia)*ct*sum(ds%imtz(:)*wf_wf(:))
                  
                  ! z-component of s-F, part 1
                  re = re + 2.0_dp*cf*ns(ia)*sum(ds%refz(:)*wf_wf(:))
                  im = im + 2.0_dp*cf*ns(ia)*sum(ds%imfz(:)*wf_wf(:))
                  
                  
                  ! phi-component of j-j
                  aux(:) = y(:)*wf_wf(:)
                  re = re + 2.0_dp*cj*(xla+xlb)*sum(ds%rejp(:)*aux(:))
                  im = im + 2.0_dp*cj*(xla+xlb)*sum(ds%imjp(:)*aux(:))
                  
                  ! z-component of s-curl-j, part 1b
                  !aux(:) = y(:)*wf_wf(:)
                  re = re - 2.0_dp*csdj*ns(ia)*(xlb-xla)*sum(ds%imjr(:)*aux(:))
                  im = im + 2.0_dp*csdj*ns(ia)*(xlb-xla)*sum(ds%rejr(:)*aux(:))
                  
                  ! rho-divergence-J, part 1c
                  !aux(:) = y(:)*wf_wf(:)
                  re = re + 2.0_dp*crdj*(xlb-xla)*sum((ds%imtjzr(:)-ds%imtjrz(:))*aux(:))
                  im = im - 2.0_dp*crdj*(xlb-xla)*sum((ds%retjzr(:)-ds%retjrz(:))*aux(:))

                  ! Optimized tensor mean field (2)
                  !imthpz(:) = -ns(ia)*(xla+xlb)*y(:)*wf_wf(:)             ! phi,z; Imag
                  imthpz(:) = (-ns(ia)*(xla+xlb))*aux(:) 
                  
                  ! rho-component of j-j
                  aux(:) = dr_wf_a(:)*wf_b(:) - wf_a(:)*dr_wf_b(:)
                  re = re - 2.0_dp*cj*sum(ds%imjr(:)*aux(:))
                  im = im + 2.0_dp*cj*sum(ds%rejr(:)*aux(:))
                  
                  ! Optimized tensor mean field (1)
                  !rethrz(:) =  ns(ia)*(dr_wf_a(:)*wf_b(:) - wf_a(:)*dr_wf_b(:))  ! rho,z; Real
                  rethrz(:) =  ns(ia)*aux(:)


                  ! rho-divergence-J, part 1a
                  aux(:) = wf_a(:)*dr_wf_b(:) + dr_wf_a(:)*wf_b(:)
                  re = re - 2.0_dp*crdj*sum((ds%retjpz(:)-ds%retjzp(:))*aux(:))
                  im = im - 2.0_dp*crdj*sum((ds%imtjpz(:)-ds%imtjzp(:))*aux(:))
                  
                  ! z-component of s-curl-j, part 1a
                  !aux(:) = dr_wf_a(:)*wf_b(:) + wf_a(:)*dr_wf_b(:)
                  re = re - 2.0_dp*csdj*ns(ia)*sum(ds%rejp(:)*aux(:))
                  im = im - 2.0_dp*csdj*ns(ia)*sum(ds%imjp(:)*aux(:))
                  
                  
                  ! z-component of j-j
                  aux(:) = dz_wf_a(:)*wf_b(:) - wf_a(:)*dz_wf_b(:)
                  re = re - 2.0_dp*cj*sum(ds%imjz(:)*aux(:))
                  im = im + 2.0_dp*cj*sum(ds%rejz(:)*aux(:))
                  
                  ! Optimized tensor mean field (3)
                  !rethzz(:) =  ns(ia)*(dz_wf_a(:)*wf_b(:) - wf_a(:)*dz_wf_b(:))  ! z,z;   Real
                  rethzz(:) =  ns(ia)*aux(:)


                  ! rho-divergence-J, part 1b
                  aux(:) = wf_a(:)*dz_wf_b(:) + dz_wf_a(:)*wf_b(:)
                  re = re - 2.0_dp*crdj*sum((ds%retjrp(:)-ds%retjpr(:))*aux(:))
                  im = im - 2.0_dp*crdj*sum((ds%imtjrp(:)-ds%imtjpr(:))*aux(:))

                  ! z-component of (Del.s)^2
                  !aux(:) = dz_wf_a(:)*wf_b(:) + wf_a(:)*dz_wf_b(:)
                  re = re + 4.0_dp*cgs*ns(ia)*sum(ds%regs(:)*aux(:))
                  im = im + 4.0_dp*cgs*ns(ia)*sum(ds%imgs(:)*aux(:))


                  ! z-component of rho-divergence-J, part 2
                  aux(:) = y(:)*(xla*wf_a(:)*dr_wf_b(:) + xlb*dr_wf_a(:)*wf_b(:))
                  re = re + 2.0_dp*crdj*ns(ia)*sum(ds%rerho(:)*aux(:))
                  im = im + 2.0_dp*crdj*ns(ia)*sum(ds%imrho(:)*aux(:))

                  ! s-curl-j, part 2a
                  !aux(:) = y(:)*(xla*wf_a(:)*dr_wf_b(:) + xlb*dr_wf_a(:)*wf_b(:))
                  re = re + 2.0_dp*csdj*sum(ds%resz(:)*aux(:))
                  im = im + 2.0_dp*csdj*sum(ds%imsz(:)*aux(:))
                  
                  
                  ! s-curl-j, part 2b
                  aux(:) = y(:)*(xlb*dz_wf_a(:)*wf_b(:) + xla*wf_a(:)*dz_wf_b(:))
                  re = re - 2.0_dp*csdj*sum(ds%resr(:)*aux(:))
                  im = im - 2.0_dp*csdj*sum(ds%imsr(:)*aux(:))

                  ! z-component of s-F, part 2c
                  aux(:) = y(:)*(xlb*dz_wf_a(:)*wf_b(:) - xla*dz_wf_b(:)*wf_a(:))
                  re = re - cf*ns(ia)*sum(ds%imsp(:)*aux(:))
                  im = im + cf*ns(ia)*sum(ds%resp(:)*aux(:))
                  
                  
                  ! rho-laplacian-rho
                  aux(:) = d2_wf_a(:)*wf_b(:) + wf_a(:)*d2_wf_b(:) &
                         + 2.0_dp*(dr_wf_a(:)*dr_wf_b(:) + dz_wf_a(:)*dz_wf_b(:)) &
                         - (xla-xlb)**2*y(:)*y(:)*wf_wf(:)
                  re = re + 4.0_dp*cdrho*sum(ds%rerho(:)*aux(:))
                  im = im + 4.0_dp*cdrho*sum(ds%imrho(:)*aux(:))

                  ! z-component of s-laplacian-s
                  !aux(:) = d2_wf_a(:)*wf_b(:) + wf_a(:)*d2_wf_b(:) &
                  !       + 2.0_dp*(dr_wf_a(:)*dr_wf_b(:) + dz_wf_a(:)*dz_wf_b(:)) &
                  !       - (xla-xlb)**2*y(:)*y(:)*wf_wf(:)
                  re = re + 4.0_dp*cds*ns(ia)*sum(ds%resz(:)*aux(:))
                  im = im + 4.0_dp*cds*ns(ia)*sum(ds%imsz(:)*aux(:))
                  
                  ! z-component of s-F, part 2b
                  aux(:) = dz_wf_a(:)*dz_wf_b(:)
                  re = re + 2.0_dp*cf*ns(ia)*sum(ds%resz(:)*aux(:))
                  im = im + 2.0_dp*cf*ns(ia)*sum(ds%imsz(:)*aux(:))
                  
                  ! rho-tau
                  !aux(:) = dr_wf_a(:)*dr_wf_b(:) + dz_wf_a(:)*dz_wf_b(:) + xla*xlb*y(:)*y(:)*wf_wf(:)
                  aux(:) = aux(:) + dr_wf_a(:)*dr_wf_b(:) + xla*xlb*y(:)*y(:)*wf_wf(:)
                  re = re + 2.0_dp*ctau*(sum(ds%retau(:)*wf_wf(:)) + sum(ds%rerho(:)*aux(:)))
                  im = im + 2.0_dp*ctau*(sum(ds%imtau(:)*wf_wf(:)) + sum(ds%imrho(:)*aux(:)))

                  !aux(:) = dr_wf_a(:)*dr_wf_b(:) + dz_wf_a(:)*dz_wf_b(:) + xla*xlb*y(:)*y(:)*wf_wf(:)
                  re = re + 2.0_dp*ns(ia)*ct*sum(ds%resz(:)*aux(:))
                  im = im + 2.0_dp*ns(ia)*ct*sum(ds%imsz(:)*aux(:))


                  ! s-curl-j, part 2c
                  aux(:) = dr_wf_a(:)*dz_wf_b(:) - dz_wf_a(:)*dr_wf_b(:)
                  re = re - 2.0_dp*csdj*sum(ds%imsp*aux(:))
                  im = im + 2.0_dp*csdj*sum(ds%resp*aux(:))


                  ! z-component of s-F, part 2a
                  aux(:) = dz_wf_a(:)*dr_wf_b(:) + dr_wf_a(:)*dz_wf_b(:)
                  re = re + cf*ns(ia)*sum(ds%resr(:)*aux(:))
                  im = im + cf*ns(ia)*sum(ds%imsr(:)*aux(:))

                  
                  ! Optimized tensor mean field (4)
                  rethrr = 0;  imthrp = 0;  imthpr = 0;  rethpp = 0;  rethzr = 0;  imthzp = 0

               ! Not-diagonal in spin
               else

                  ! |a> = |+>,  |b> = |->
                  if (ns(ia) == 1 .and. ns(ib) == -1) then

                     ! rho-component of s-s
                     aux(:) = cs(:)*wf_wf(:)
                     re = re + 4.0_dp*sum(ds%resr(:)*aux(:))
                     im = im + 4.0_dp*sum(ds%imsr(:)*aux(:))
                     ! phi-component of s-s
                     !aux(:) = cs(:)*wf_wf(:)
                     re = re + 4.0_dp*sum(ds%imsp(:)*aux(:))
                     im = im - 4.0_dp*sum(ds%resp(:)*aux(:))

                     ! rho-component of s-T (1)
                     re = re + 2.0_dp*ct*sum(ds%retr(:)*wf_wf(:))
                     im = im + 2.0_dp*ct*sum(ds%imtr(:)*wf_wf(:))
                     ! phi-component of s-T (1)
                     re = re + 2.0_dp*ct*sum(ds%imtp(:)*wf_wf(:))
                     im = im - 2.0_dp*ct*sum(ds%retp(:)*wf_wf(:))
                      ! rho-component of s-F, part 1
                     re = re + 2.0_dp*cf*sum(ds%refr(:)*wf_wf(:))
                     im = im + 2.0_dp*cf*sum(ds%imfr(:)*wf_wf(:))
                     ! phi-component of s-F, part 1
                     re = re + 2.0_dp*cf*sum(ds%imfp(:)*wf_wf(:))
                     im = im - 2.0_dp*cf*sum(ds%refp(:)*wf_wf(:))
                     
                     ! rho-component of s-curl-j (2)
                     aux(:) = y(:)*wf_wf(:)
                     re = re + 2.0_dp*csdj*(xlb-xla)*sum(ds%imjz(:)*aux(:))
                     im = im - 2.0_dp*csdj*(xlb-xla)*sum(ds%rejz(:)*aux(:))
                     ! Optimized tensor mean field (2)
                     !imthpr(:) = -(xla+xlb)*y(:)*wf_wf(:)              ! phi,rho; Imag
                     !rethpp(:) = -(xla+xlb)*y(:)*wf_wf(:)              ! phi,phi; Real
                     imthpr(:) = -(xla+xlb)*aux(:)
                     rethpp(:) = imthpr(:)
                     ! phi-component of s-F, part 2b
                     !aux(:) = y(:)**2 * wf_wf(:)
                     aux(:) = y(:)*aux(:)
                     re = re + cf*2.0_dp*xla*xlb*sum(ds%imsp(:)*aux(:))
                     im = im - cf*2.0_dp*xla*xlb*sum(ds%resp(:)*aux(:))
                     

                     ! rho-component of s-curl-j (1)
                     aux(:) = dz_wf_a(:)*wf_b(:) + wf_a(:)*dz_wf_b(:)
                     re = re + 2.0_dp*csdj*sum(ds%rejp(:)*aux(:))
                     im = im + 2.0_dp*csdj*sum(ds%imjp(:)*aux(:))
                     ! phi-component of s-curl-j (2)
                     !aux(:) = dz_wf_a(:)*wf_b(:) + wf_a(:)*dz_wf_b(:)
                     re = re - 2.0_dp*csdj*sum(ds%imjr(:)*aux(:))
                     im = im + 2.0_dp*csdj*sum(ds%rejr(:)*aux(:))
                     
                     ! Optimized tensor mean field (3)
                     rethzr(:) =  dz_wf_a(:)*wf_b(:) - wf_a(:)*dz_wf_b(:)     ! z,rho;   Real
                     !imthzp(:) = -(dz_wf_a(:)*wf_b(:) - wf_a(:)*dz_wf_b(:))   ! z,phi;   Imag
                     imthzp(:) = -rethzr(:)
                     
                     ! phi-component of s-curl-j (1)
                     aux(:) = dr_wf_a(:)*wf_b(:) + wf_a(:)*dr_wf_b(:)
                     re = re + 2.0_dp*csdj*sum(ds%imjz(:)*aux(:))
                     im = im - 2.0_dp*csdj*sum(ds%rejz(:)*aux(:))
                     ! Optimized tensor mean field (1)
                     rethrr(:) =  dr_wf_a(:)*wf_b(:) - wf_a(:)*dr_wf_b(:)     ! rho,rho; Real
                     imthrp(:) = -(dr_wf_a(:)*wf_b(:) - wf_a(:)*dr_wf_b(:))   ! rho,phi; Imag
                     
                     ! rho- and phi-components of (Del.s)^2
                     !aux(:) = dr_wf_a(:)*wf_b(:) + wf_a(:)*dr_wf_b(:) - (xla-xlb)*y(:)*wf_wf(:)
                     aux(:) = aux(:) - (xla-xlb)*y(:)*wf_wf(:)
                     re = re + 4.0_dp*cgs*sum(ds%regs(:)*aux(:))
                     im = im + 4.0_dp*cgs*sum(ds%imgs(:)*aux(:))
                     
                     
                     ! rho-component of s-F, part 2c
                     aux(:) = y(:)*(xlb*dr_wf_a(:)*wf_b(:) - xla*wf_a(:)*dr_wf_b(:))
                     re = re - cf*sum(ds%imsp(:)*aux(:))
                     im = im + cf*sum(ds%resp(:)*aux(:))
                     ! phi-component of s-F, part 2a
                     !aux(:) = y(:)*(xla*wf_a(:)*dr_wf_b(:) - xlb*wf_b(:)*dr_wf_a(:))
                     !re = re - cf*sum(ds%resr(:)*aux(:))
                     !im = im - cf*sum(ds%imsr(:)*aux(:))
                     re = re + cf*sum(ds%resr(:)*aux(:))
                     im = im + cf*sum(ds%imsr(:)*aux(:))
                     
                     ! phi-component of s-F, part 2c
                     aux(:) = y(:)*(xla*wf_a(:)*dz_wf_b(:) - xlb*wf_b(:)*dz_wf_a(:))
                     re = re - cf*sum(ds%resz(:)*aux(:))
                     im = im - cf*sum(ds%imsz(:)*aux(:))
                     ! rho-component of rho-divergence-J, part 2
                     aux(:) = y(:)*(xla*wf_a(:)*dz_wf_b(:) + xlb*dz_wf_a(:)*wf_b(:))
                     re = re - 2.0_dp*crdj*sum(ds%rerho(:)*aux(:))
                     im = im - 2.0_dp*crdj*sum(ds%imrho(:)*aux(:))
                     

                     ! rho-component of s-F, part 2b
                     aux(:) = dr_wf_a(:)*dz_wf_b(:) + dz_wf_a(:)*dr_wf_b(:)
                     re = re + cf*sum(ds%resz(:)*aux(:))
                     im = im + cf*sum(ds%imsz(:)*aux(:))
                     
                     ! phi-component of rho-divergence-J, part 2
                     aux(:) = dr_wf_a(:)*dz_wf_b(:) - dz_wf_a(:)*dr_wf_b(:)
                     re = re + 2.0_dp*crdj*sum(ds%rerho(:)*aux(:))
                     im = im + 2.0_dp*crdj*sum(ds%imrho(:)*aux(:))
                     
                     ! rho-component of s-F, part 2a
                     aux(:) = dr_wf_a(:)*dr_wf_b(:)
                     re = re + cf*2.0_dp*sum(ds%resr(:)*aux(:))
                     im = im + cf*2.0_dp*sum(ds%imsr(:)*aux(:))
                     
                     ! rho-component of s-T (2)
                     !aux(:) = dr_wf_a(:)*dr_wf_b(:) + dz_wf_a(:)*dz_wf_b(:) + xla*xlb*y(:)*y(:)*wf_wf(:)
                     aux(:) = aux(:) + dz_wf_a(:)*dz_wf_b(:) + xla*xlb*y(:)*y(:)*wf_wf(:)
                     re = re + 2.0_dp*ct*sum(ds%resr(:)*aux(:))
                     im = im + 2.0_dp*ct*sum(ds%imsr(:)*aux(:))
                     ! phi-component of s-T (2)
                     !aux(:) = dr_wf_a(:)*dr_wf_b(:) + dz_wf_a(:)*dz_wf_b(:) + xla*xlb*y(:)*y(:)*wf_wf(:)
                     re = re + 2.0_dp*ct*sum(ds%imsp(:)*aux(:))
                     im = im - 2.0_dp*ct*sum(ds%resp(:)*aux(:))


                     ! rho-component of s-laplacian-s
                     aux(:) = d2_wf_a(:)*wf_b(:) + wf_a(:)*d2_wf_b(:) + 2.0_dp*(dr_wf_a(:)*dr_wf_b(:)  &
                                       + dz_wf_a(:)*dz_wf_b(:)) - (xla-xlb)**2*y(:)*y(:)*wf_wf(:)
                     re = re + 4.0_dp*cds*sum(ds%resr(:)*aux(:))
                     im = im + 4.0_dp*cds*sum(ds%imsr(:)*aux(:))
                     ! phi-component of s-laplacian-s
                     !aux(:) = d2_wf_a(:)*wf_b(:) + wf_a(:)*d2_wf_b(:) + 2.0_dp*(dr_wf_a(:)*dr_wf_b(:)  &
                     !                  + dz_wf_a(:)*dz_wf_b(:)) - (xla-xlb)**2*y(:)*y(:)*wf_wf(:)
                     re = re + 4.0_dp*cds*sum(ds%imsp(:)*aux(:))
                     im = im - 4.0_dp*cds*sum(ds%resp(:)*aux(:))

                     
                     ! Optimized tensor mean field (4)
                     rethrz = 0;  imthpz = 0;  rethzz = 0


                  ! |a> = |->,  |b> = |+>
                  else if (ns(ia) == -1 .and. ns(ib) == 1) then

                     ! rho-component of s-s
                     aux(:) = cs(:)*wf_wf(:)
                     re = re + 4.0_dp*sum(ds%resr(:)*aux(:))
                     im = im + 4.0_dp*sum(ds%imsr(:)*aux(:))
                     ! phi-component of s-s
                     !aux(:) = cs(:)*wf_wf(:)
                     re = re - 4.0_dp*sum(ds%imsp(:)*aux(:))
                     im = im + 4.0_dp*sum(ds%resp(:)*aux(:))

                     ! rho-component of s-T (1)
                     re = re + 2.0_dp*ct*sum(ds%retr(:)*wf_wf(:))
                     im = im + 2.0_dp*ct*sum(ds%imtr(:)*wf_wf(:))
                     ! phi-component of s-T (1)
                     re = re - 2.0_dp*ct*sum(ds%imtp(:)*wf_wf(:))
                     im = im + 2.0_dp*ct*sum(ds%retp(:)*wf_wf(:))
                     ! rho-component of s-F, part 1
                     re = re + 2.0_dp*cf*sum(ds%refr(:)*wf_wf(:))
                     im = im + 2.0_dp*cf*sum(ds%imfr(:)*wf_wf(:))
                     ! phi-component of s-F, part 1
                     re = re - 2.0_dp*cf*sum(ds%imfp(:)*wf_wf(:))
                     im = im + 2.0_dp*cf*sum(ds%refp(:)*wf_wf(:))
                     
                     ! rho-component of s-curl-j (2)
                     aux(:) = y(:)*wf_wf(:)
                     re = re + 2.0_dp*csdj*(xlb-xla)*sum(ds%imjz(:)*aux(:))
                     im = im - 2.0_dp*csdj*(xlb-xla)*sum(ds%rejz(:)*aux(:))
                     
                     ! Optimized tensor mean field (2)
                     !imthpr(:) = -(xla+xlb)*y(:)*wf_wf(:)           ! phi,rho; Imag
                     !rethpp(:) =  (xla+xlb)*y(:)*wf_wf(:)           ! phi,phi; Real
                     rethpp(:) =  (xla+xlb)*aux(:)           ! phi,phi; Real
                     imthpr(:) = -rethpp(:)           ! phi,rho; Imag
                     
                     ! phi-component of s-F, part 2b
                     !aux(:) = y(:)*y(:)*wf_wf(:)
                     aux(:) = y(:)*aux(:)
                     re = re - 2.0_dp*cf*xla*xlb*sum(ds%imsp(:)*aux(:))
                     im = im + 2.0_dp*cf*xla*xlb*sum(ds%resp(:)*aux(:))
                     
                     
                     ! rho-component of s-T (2)
                     aux(:) = dr_wf_a(:)*dr_wf_b(:) + dz_wf_a(:)*dz_wf_b(:) + xla*xlb*y(:)*y(:)*wf_wf(:)
                     re = re + 2.0_dp*ct*sum(ds%resr(:)*aux(:))
                     im = im + 2.0_dp*ct*sum(ds%imsr(:)*aux(:))
                     ! phi-component of s-T (2)
                     !aux(:) = dr_wf_a(:)*dr_wf_b(:) + dz_wf_a(:)*dz_wf_b(:) + xla*xlb*y(:)*y(:)*wf_wf(:)
                     re = re - 2.0_dp*ct*sum(ds%imsp(:)*aux(:))
                     im = im + 2.0_dp*ct*sum(ds%resp(:)*aux(:))

                     ! rho-component of s-laplacian-s
                     aux(:) = d2_wf_a(:)*wf_b(:) + wf_a(:)*d2_wf_b(:) + 2.0_dp*(dr_wf_a(:)*dr_wf_b(:)  &
                            + dz_wf_a(:)*dz_wf_b(:)) - (xla-xlb)**2*y(:)*y(:)*wf_wf(:)
                     re = re + 4.0_dp*cds*sum(ds%resr(:)*aux(:))
                     im = im + 4.0_dp*cds*sum(ds%imsr(:)*aux(:))
                     ! phi-component of s-laplacian-s
                     !aux(:) = d2_wf_a(:)*wf_b(:) + wf_a(:)*d2_wf_b(:) + 2.0_dp*(dr_wf_a(:)*dr_wf_b(:)  &
                     !       + dz_wf_a(:)*dz_wf_b(:)) - (xla-xlb)**2*y(:)*y(:)*wf_wf(:)
                     re = re - 4.0_dp*cds*sum(ds%imsp(:)*aux(:))
                     im = im + 4.0_dp*cds*sum(ds%resp(:)*aux(:))
                     
                     ! phi-component of rho-divergence-J, part 2
                     aux(:) = dr_wf_a(:)*dz_wf_b(:) - dz_wf_a(:)*dr_wf_b(:)
                     re = re - 2.0_dp*crdj*sum(ds%rerho(:)*aux(:))
                     im = im - 2.0_dp*crdj*sum(ds%imrho(:)*aux(:))
                     
                     ! rho-component of s-F, part 2a
                     aux(:) = dr_wf_a(:)*dr_wf_b(:)
                     re = re + cf*2.0_dp*sum(ds%resr(:)*aux(:))
                     im = im + cf*2.0_dp*sum(ds%imsr(:)*aux(:))
                     
                     ! rho-component of s-F, part 2b
                     aux(:) = dr_wf_a(:)*dz_wf_b(:) + dz_wf_a(:)*dr_wf_b(:)
                     re = re + cf*sum(ds%resz(:)*aux(:))
                     im = im + cf*sum(ds%imsz(:)*aux(:))
                     
                     
                     
                     ! rho-component of rho-divergence-J, part 2
                     aux(:) = y(:)*(xla*wf_a(:)*dz_wf_b(:) + xlb*dz_wf_a(:)*wf_b(:))
                     re = re - 2.0_dp*crdj*sum(ds%rerho(:)*aux(:))
                     im = im - 2.0_dp*crdj*sum(ds%imrho(:)*aux(:))
                     
                     ! phi-component of s-F, part 2c
                     aux(:) = y(:)*(xla*wf_a(:)*dz_wf_b(:) - xlb*wf_b(:)*dz_wf_a(:))
                     re = re + cf*sum(ds%resz(:)*aux(:))
                     im = im + cf*sum(ds%imsz(:)*aux(:))
                     
                     ! phi-component of s-curl-j (1)
                     aux(:) = dr_wf_a(:)*wf_b(:) + wf_a(:)*dr_wf_b(:)
                     re = re - 2.0_dp*csdj*sum(ds%imjz(:)*aux(:))
                     im = im + 2.0_dp*csdj*sum(ds%rejz(:)*aux(:))
                     
                     ! rho- and phi-components of (Del.s)^2
                     !aux(:) = dr_wf_a(:)*wf_b(:) + wf_a(:)*dr_wf_b(:) + y(:)*(xla-xlb)*wf_wf(:)
                     aux(:) = aux(:) + (xla-xlb)*(y(:)*wf_wf(:))
                     re = re + 4.0_dp*cgs*sum(ds%regs(:)*aux(:))
                     im = im + 4.0_dp*cgs*sum(ds%imgs(:)*aux(:))
                     
                     
                     ! rho-component of s-curl-j (1)
                     aux(:) = dz_wf_a(:)*wf_b(:) + wf_a(:)*dz_wf_b(:)
                     re = re + 2.0_dp*csdj*sum(ds%rejp(:)*aux(:))
                     im = im + 2.0_dp*csdj*sum(ds%imjp(:)*aux(:))
                     ! phi-component of s-curl-j (2)
                     !aux(:) = dz_wf_a(:)*wf_b(:) + wf_a(:)*dz_wf_b(:)
                     re = re + 2.0_dp*csdj*sum(ds%imjr(:)*aux(:))
                     im = im - 2.0_dp*csdj*sum(ds%rejr(:)*aux(:))

                     ! Optimized tensor mean field (1)
                     rethrr(:) =  dr_wf_a(:)*wf_b(:) - wf_a(:)*dr_wf_b(:)  ! rho,rho; Real
                     !imthrp(:) =  dr_wf_a(:)*wf_b(:) - wf_a(:)*dr_wf_b(:)  ! rho,phi; Imag
                     imthrp(:) = rethrr(:)
                     
                     ! Optimized tensor mean field (3)
                     rethzr(:) =  dz_wf_a(:)*wf_b(:) - wf_a(:)*dz_wf_b(:)  ! z,rho;   Real
                     !imthzp(:) =  dz_wf_a(:)*wf_b(:) - wf_a(:)*dz_wf_b(:)  ! z,phi;   Imag
                     imthzp(:) = rethzr(:)
                     
                     ! rho-component of s-F, part 2c
                     aux(:) = y(:)*(xlb*dr_wf_a(:)*wf_b(:) - xla*wf_a(:)*dr_wf_b(:))
                     re = re - cf*sum(ds%imsp(:)*aux(:))
                     im = im + cf*sum(ds%resp(:)*aux(:))
                     ! phi-component of s-F, part 2a
                     !aux(:) = y(:)*(xla*wf_a(:)*dr_wf_b(:) - xlb*wf_b(:)*dr_wf_a(:))
                     !re = re + cf*sum(ds%resr(:)*aux(:))
                     !im = im + cf*sum(ds%imsr(:)*aux(:))
                     re = re - cf*sum(ds%resr(:)*aux(:))
                     im = im - cf*sum(ds%imsr(:)*aux(:))
                     
                     ! Optimized tensor mean field (4)
                     rethrz = 0;  imthpz = 0;  rethzz = 0

                  else
                     call abort(' Inconsistent spin in helpers.meanfield')
                  end if
               end if


               ! Tensor terms in the mean field
               ! Pseudoscalar J**2 term
               aux(:) = rethrr(:) + rethpp(:) + rethzz(:)
               re = re - 2.0_dp*ctj0*sum((ds%imtjrr(:)+ds%imtjpp(:)+ds%imtjzz(:))*aux(:))
               im = im + 2.0_dp*ctj0*sum((ds%retjrr(:)+ds%retjpp(:)+ds%retjzz(:))*aux(:))

               ! Vector J**2 term
               re = re - 2.0_dp*ctj1*(sum((ds%retjpz(:)-ds%retjzp(:))*(imthpz(:)-imthzp(:)))   &
                                 + sum((ds%retjrp(:)-ds%retjpr(:))*(imthrp(:)-imthpr(:)))   &
                                 + sum((ds%imtjzr(:)-ds%imtjrz(:))*(rethzr(:)-rethrz(:))))
               im = im - 2.0_dp*ctj1*(sum((ds%imtjpz(:)-ds%imtjzp(:))*(imthpz(:)-imthzp(:)))   &
                                 + sum((ds%imtjrp(:)-ds%imtjpr(:))*(imthrp(:)-imthpr(:)))   &
                                 - sum((ds%retjzr(:)-ds%retjrz(:))*(rethzr(:)-rethrz(:))))

               ! Pseudotensor J**2 term
               ! rho,rho-component
               aux(:) = 2.0_dp*rethrr(:) - rethpp(:) - rethzz(:)
               re = re - 2.0_dp/9.0_dp*ctj2*sum((2.0_dp*ds%imtjrr(:)-ds%imtjpp(:)-ds%imtjzz(:)) * aux(:))
               im = im + 2.0_dp/9.0_dp*ctj2*sum((2.0_dp*ds%retjrr(:)-ds%retjpp(:)-ds%retjzz(:)) * aux(:))

               ! phi,phi-component
               aux(:) = 2.0_dp*rethpp(:) - rethzz(:) - rethrr(:)
               re = re - 2.0_dp/9.0_dp*ctj2*sum((2.0_dp*ds%imtjpp(:)-ds%imtjzz(:)-ds%imtjrr(:)) * aux(:))
               im = im + 2.0_dp/9.0_dp*ctj2*sum((2.0_dp*ds%retjpp(:)-ds%retjzz(:)-ds%retjrr(:)) * aux(:))

               ! z,z-component
               aux(:) = 2.0_dp*rethzz(:) - rethrr(:) - rethpp(:)
               re = re - 2.0_dp/9.0_dp*ctj2*sum((2.0_dp*ds%imtjzz(:)-ds%imtjrr(:)-ds%imtjpp(:)) * aux(:))
               im = im + 2.0_dp/9.0_dp*ctj2*sum((2.0_dp*ds%retjzz(:)-ds%retjrr(:)-ds%retjpp(:)) * aux(:))

               ! (rho,z + z,rho)-component
               aux(:) = rethrz(:) + rethzr(:)
               re = re - ctj2*sum((ds%imtjrz(:)+ds%imtjzr(:))*aux(:))
               im = im + ctj2*sum((ds%retjrz(:)+ds%retjzr(:))*aux(:))

               ! (rho,phi + phi,rho)-component
               aux(:) = imthrp(:) + imthpr(:)
               re = re - ctj2*sum((ds%retjrp(:)+ds%retjpr(:))*aux(:))
               im = im - ctj2*sum((ds%imtjrp(:)+ds%imtjpr(:))*aux(:))

               ! (phi,z + z,phi)-component
               aux(:) = imthpz(:)+imthzp(:)
               re = re - ctj2*sum((ds%retjpz(:)+ds%retjzp(:))*aux(:))
               im = im - ctj2*sum((ds%imtjpz(:)+ds%imtjzp(:))*aux(:))
               
               reh%elem(ix) = re
               imh%elem(ix) = im

      end do
      !$OMP END PARALLEL DO
      
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

      integer :: im, ia, ib, nqp
      real(dp) :: wf_a(nghl), wf_b(nghl), aux(nghl)

      red%elem = 0 ; imd%elem = 0
      nqp = sum(db)

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ia,ib,wf_a,wf_b,aux)
      do im = 1, size(red%elem)
         ia = ia_tabp(im)
         ib = ib_tabp(im)
         
         ! Column vector |b>
         wf_b = wf(:, ib)

         ! Row vector <a|
         wf_a = wf(:, ia)

         ! Diagonal in spin
         if (ns(ia) == -ns(ib)) then
            ! rho-breve_10
            aux = 2*cpair*wf_a*wf_b*ns(ia)
            red%elem(im) = red%elem(im) + sum(aux*ds%rerb)
            imd%elem(im) = imd%elem(im) + sum(aux*ds%imrb)
         
            ! s-breve_{0,z}
            aux = 2*cspair*wf_a*wf_b
            red%elem(im) = red%elem(im) - sum(aux*ds%resbz)
            imd%elem(im) = imd%elem(im) - sum(aux*ds%imsbz) 
         else
            if (ns(ia) == -1 .and. ns(ib) == -1) then
               ! s-breve_{0,rho}
               aux = 2*cspair*wf_a*wf_b
               red%elem(im) = red%elem(im) - sum(aux*ds%resbr)
               imd%elem(im) = imd%elem(im) - sum(aux*ds%imsbr)
      
               ! s-breve_{0,phi}
               red%elem(im) = red%elem(im) + sum(aux*ds%imsbp)
               imd%elem(im) = imd%elem(im) - sum(aux*ds%resbp)
            else
               ! s-breve_{0,rho}
               aux = 2*cspair*wf_a*wf_b
               red%elem(im) = red%elem(im) + sum(aux*ds%resbr)
               imd%elem(im) = imd%elem(im) + sum(aux*ds%imsbr)
         
               ! s-breve_{0,phi}
               red%elem(im) = red%elem(im) + sum(aux*ds%imsbp)
               imd%elem(im) = imd%elem(im) - sum(aux*ds%resbp)
            end if
         end if
      
      end do
      !$OMP END PARALLEL DO

   end subroutine pairingfield

end module pnfam_hamiltonian
