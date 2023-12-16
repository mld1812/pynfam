!------------------------------------------------------------------------------
!> This module contains a derived type to store spatial pieces of the 2 body
!> current matrix elements and methods to calculate the z and r components of them.
!> It also contains the methods to combine these matrix elements in a way that
!> is consistent with the contraction over spin and isospin for the mean field.
!>
!> @authors E.M. Ney, UNC Chapel Hill, 2018-
!------------------------------------------------------------------------------
module type_extfield_2bc
   use pnfam_setup, only : two_body_current_lecs, two_body_current_usep, st
   use pnfam_spatial_mtxels
   use pnfam_constants, only : Fpi, gA, Mpi, Mn, hbarc
   use pnfam_logger
   implicit none
   integer, parameter, private :: dp = kind(1.0d0)
   logical, private :: init=.false.

   type extfield_2bc
      real(dp) :: rkp, rkz, rkm
      real(dp) :: rpz, rpm, rzm
      real(dp) :: r0p, r0z, r0m
      real(dp) :: rp0, rz0, rm0
      real(dp) :: rpz2, rzm2
      integer :: k
   end type extfield_2bc

   ! Shortcuts for double precision numbers
   real(dp), parameter :: q=0.25_dp, h=0.5_dp, t=2.0_dp, cr2=sqrt(2.0_dp)
   real(dp), parameter :: c2r2=2.0_dp*cr2, c4r2=4.0_dp*cr2

   ! Low energy constants
   logical :: use_p
   real(dp) :: c3, c4, c4p14, d1, d2, cd, caux
   ! Combinations of low energy constants
   real(dp) :: cY0p, cY0p_c2r2, cY0p_2, cY0p_4
   real(dp) :: cdel, cdel_c2r2, cdel_2
   real(dp) :: caux_d1_c4r2, caux_d1_8
   real(dp) :: caux_c3_c4r2, caux_c3_8
   real(dp) :: caux_cr2, caux_2

   type(extfield_2bc), Public, allocatable :: MEJz(:,:)
   type(extfield_2bc), Public, allocatable :: MEJr(:,:)
   real(dp) :: NumVz, NumVr

contains

   subroutine init_extfield_2bc_type(mtxels)
      implicit none
      logical, intent(in) :: mtxels

         ! PRC 90 031301(R) 2018: (c3,c4)[Ref] = (-3.2,5.4)[37], (-4.78,3.96)[38], (-3.4,3.4)[39] in 1/GeV!
         c3    = two_body_current_lecs(1)
         c4    = two_body_current_lecs(2)
         c4p14 = c4 + 0.25_dp ! c4 + 1/4
         ! PRC 90 031301(R) 2018: cd varied from -2.5 to 2.5
         cd    = two_body_current_lecs(3)*(-0.25_dp) ! Extra factor of -1/4 in contact term (Eq.5.8, Annals of Phys. 378, (2017) 317-395)
         d1    = cd*0.50_dp ! arbitrarily split so cd = d1 + 2d2 (doesn't matter unless debug>1)
         d2    = cd*0.25_dp
         use_p = two_body_current_usep

         ! Some auxiliary constants (Convert all MeV to 1/fm)
         caux = (hbarc*hbarc*hbarc)/(2.0_dp*Mn*Fpi*Fpi) ! factor out gA
         cY0p = caux*(+4.0_dp*c4p14) ! J_rsx/caux = 4i*[ -(c4 + 1/4)i*J*Y0 ] --> +4cp14
         cdel = caux*(-8.0_dp*d2)    ! J_rsx/caux = 4i*[ +2d2*i*sx*delta ] --> -8d2

         ! Further auxialliary constants to avoid operations under loops
         cY0p_c2r2    = cY0p*c2r2
         cY0p_2       = cY0p*2.0_dp
         cY0p_4       = cY0p*4.0_dp
         cdel_c2r2    = cdel*c2r2
         cdel_2       = cdel*2.0_dp
         caux_d1_c4r2 = caux*d1*c4r2
         caux_d1_8    = caux*d1*8.0_dp
         caux_c3_c4r2 = caux*c3*c4r2
         caux_c3_8    = caux*c3*8.0_dp
         caux_cr2     = caux*cr2
         caux_2       = caux*2.0_dp

      if (.not. init .and. mtxels) then
         call prep_gaussian
         init = .true.
      end if

   end subroutine

   subroutine zero_extfield_2bc(J)
      implicit none
      type(extfield_2bc), intent(inout) :: J

      real(dp), parameter :: zero=0.0_dp

      J%rkp=zero ; J%rkz=zero ; J%rkm=zero
      J%rpz=zero ; J%rpm=zero ; J%rzm=zero
      J%r0p=zero ; J%r0z=zero ; J%r0m=zero
      J%rp0=zero ; J%rz0=zero ; J%rm0=zero
      J%rpz2=zero; J%rzm2=zero

   end subroutine

   subroutine multiply_extfield_2bc(J1, J2)
      implicit none
      type(extfield_2bc), intent(inout) :: J1
      type(extfield_2bc), intent(in) :: J2

      J1%rkp = J1%rkp*J2%rkp
      J1%rpz = J1%rpz*J2%rpz
      J1%r0p = J1%r0p*J2%r0p
      J1%rp0 = J1%rp0*J2%rp0
      J1%rkz = J1%rkz*J2%rkz
      J1%rpm = J1%rpm*J2%rpm
      J1%r0z = J1%r0z*J2%r0z
      J1%rz0 = J1%rz0*J2%rz0
      J1%rkm = J1%rkm*J2%rkm
      J1%rzm = J1%rzm*J2%rzm
      J1%r0m = J1%r0m*J2%r0m
      J1%rm0 = J1%rm0*J2%rm0

      J1%rpz2 = J1%rpz2*J2%rpz2
      J1%rzm2 = J1%rzm2*J2%rzm2

   end subroutine

   subroutine add_extfield_2bc(J1, J2)
      implicit none
      type(extfield_2bc), intent(inout) :: J1
      type(extfield_2bc), intent(in) :: J2

      J1%rkp = J1%rkp + J2%rkp
      J1%rpz = J1%rpz + J2%rpz
      J1%r0p = J1%r0p + J2%r0p
      J1%rp0 = J1%rp0 + J2%rp0
      J1%rkz = J1%rkz + J2%rkz
      J1%rpm = J1%rpm + J2%rpm
      J1%r0z = J1%r0z + J2%r0z
      J1%rz0 = J1%rz0 + J2%rz0
      J1%rkm = J1%rkm + J2%rkm
      J1%rzm = J1%rzm + J2%rzm
      J1%r0m = J1%r0m + J2%r0m
      J1%rm0 = J1%rm0 + J2%rm0

      J1%rpz2 = J1%rpz2 + J2%rpz2
      J1%rzm2 = J1%rzm2 + J2%rzm2

   end subroutine

   function scalar_mult_extfield_2bc(J1, s) result(J2)
      implicit none
      type(extfield_2bc), intent(in) :: J1
      real(dp), intent(in) :: s

      type(extfield_2bc) :: J2

      call zero_extfield_2bc(J2)
      J2%k = J1%k

      J2%rkp = J1%rkp*s
      J2%rpz = J1%rpz*s
      J2%r0p = J1%r0p*s
      J2%rp0 = J1%rp0*s
      J2%rkz = J1%rkz*s
      J2%rpm = J1%rpm*s
      J2%r0z = J1%r0z*s
      J2%rz0 = J1%rz0*s
      J2%rkm = J1%rkm*s
      J2%rzm = J1%rzm*s
      J2%r0m = J1%r0m*s
      J2%rm0 = J1%rm0*s

      J2%rpz2 = J1%rpz2*s
      J2%rzm2 = J1%rzm2*s

   end function

   subroutine stop_not_init
      implicit none
      call abort(" Error, type_extfield_2bc routine called before initialization")
   end subroutine

   !=========================================================================================
   ! Populate the spatial matrix elements needed for a given spherical component
   ! of the current. Calculates the Z component of the matrix elements.
   subroutine calc_Jz_opt(J,mu,ni,nj,nk,nl)

      implicit none
      type(extfield_2bc), intent(inout) :: J
      integer, intent(in) :: mu, ni,nj,nk,nl
      real(dp) :: me_p1dp, me_p1d0, me_p1dm
      real(dp) :: me_p2dp, me_p2d0, me_p2dm
      real(dp) :: me_dpp, me_dp0, me_dpm, me_d00, me_d0m, me_dmm, me_del
      character(1240) :: st

      if(.not.init) call stop_not_init

      call zero_extfield_2bc(J)

      ! P1 and P2 terms (treated same as ddY0)
      ! (Only difference is the component of P=K. Do it here to save repeat code
      !  and to easily turn off P terms)
      if (use_p) then
         me_p1dp = MatrixElement_ddGz(ni,nj,nk,nl,mu, J%k,+1, 1)
         me_p1d0 = MatrixElement_ddGz(ni,nj,nk,nl,mu, J%k, 0, 1)
         me_p1dm = MatrixElement_ddGz(ni,nj,nk,nl,mu, J%k,-1, 1)
         me_p2dp = MatrixElement_ddGz(ni,nj,nk,nl,mu, J%k,+1, 2)
         me_p2d0 = MatrixElement_ddGz(ni,nj,nk,nl,mu, J%k, 0, 2)
         me_p2dm = MatrixElement_ddGz(ni,nj,nk,nl,mu, J%k,-1, 2)
         ! p1-terms
         J%r0p = +caux_cr2*me_p1dm
         J%r0z = +caux_2  *me_p1d0
         J%r0m = -caux_cr2*me_p1dp
         ! p2-terms
         J%rp0 = +caux_cr2*me_p2dm
         J%rz0 = +caux_2  *me_p2d0
         J%rm0 = -caux_cr2*me_p2dp
      end if

      ! Matrix elements depend on the gaussian, so we must loop through terms
      ! and re-calculate each matrix element for each gaussian.
      select case (J%k)
      case default
          call abort(" ERROR, invalid K for extfield 2bc")
      !----------------------------------------------------------------------
      case(+1)

            ! 1st ddY0 terms
            me_dpm  = MatrixElement_ddGz(ni,nj,nk,nl,mu,  +1,-1, 0)
            me_d00  = MatrixElement_ddGz(ni,nj,nk,nl,mu,   0, 0, 0)
            me_dp0  = MatrixElement_ddGz(ni,nj,nk,nl,mu,  +1, 0, 0)
            me_dpp  = me_dpm ! MatrixElement_ddGz(ni,nj,nk,nl,mu,  +1,+1, 0)

            ! 1-spin terms
            J%rkp = -caux_c3_c4r2*me_dpm
            J%rkz = -caux_c3_8   *me_dp0
            J%rkm = +caux_c3_c4r2*me_dpp

            ! 2-spin upper triangle
            J%rpz = +cY0p_c2r2*me_d00
            J%rpm = -cY0p_2   *me_dp0
            J%rzm = -cY0p_c2r2*me_dpp

            ! 2nd ddY0 terms (only in 2-spin terms)
            J%rpz2 = -cY0p_c2r2*me_dpm

      !----------------------------------------------------------------------
      case(0)

            ! 1st ddY0 terms
            me_dpm  = MatrixElement_ddGz(ni,nj,nk,nl,mu,  +1,-1, 0)
            me_d00  = MatrixElement_ddGz(ni,nj,nk,nl,mu,   0, 0, 0)
            me_d0m  = MatrixElement_ddGz(ni,nj,nk,nl,mu,   0,-1, 0)
            me_dp0  = me_d0m ! MatrixElement_ddGz(ni,nj,nk,nl,mu,  +1, 0, 0)

            ! 1-spin terms
            J%rkp = -caux_c3_c4r2*me_d0m
            J%rkz = -caux_c3_8   *me_d00
            J%rkm = +caux_c3_c4r2*me_dp0

            ! 2-spin upper triangle
            J%rpz = +cY0p_c2r2*me_d0m
            J%rpm = -cY0p_4   *me_dpm
            J%rzm = -cY0p_c2r2*me_dp0

            ! no 2nd ddY0 terms for 0th component

      !----------------------------------------------------------------------
      case(-1)

            ! 1st ddY0 terms
            me_dpm  = MatrixElement_ddGz(ni,nj,nk,nl,mu,  +1,-1, 0)
            me_d00  = MatrixElement_ddGz(ni,nj,nk,nl,mu,   0, 0, 0)
            me_d0m  = MatrixElement_ddGz(ni,nj,nk,nl,mu,   0,-1, 0)
            me_dmm  = me_dpm ! MatrixElement_ddGz(ni,nj,nk,nl,mu,  -1,-1, 0)

            ! 1-spin terms
            J%rkp = -caux_c3_c4r2*me_dmm
            J%rkz = -caux_c3_8   *me_d0m
            J%rkm = +caux_c3_c4r2*me_dpm

            ! 2-spin upper triangle
            J%rpz = +cY0p_c2r2*me_dmm
            J%rpm = -cY0p_2   *me_d0m
            J%rzm = -cY0p_c2r2*me_d00

            ! 2nd ddY0 terms (only in 2-spin terms)
            J%rzm2 = +cY0p_c2r2*me_dpm

      end select

   end subroutine calc_Jz_opt

   ! Populate the spatial matrix elements needed for a given spherical component
   ! of the current. Calculates the radial component of the matrix elements.
   ! NOTE: prefactors are absent here!!!, Jz holds all prefactors, Jr holds just
   ! matrix elements such that the full thing is J_abcd = (pf*Jz)*(Jr)
   subroutine calc_Jr_opt(J,mu,ni,li,nj,lj,nk,lk,nl,ll)

      implicit none
      type(extfield_2bc), intent(inout) :: J
      integer, intent(in) :: mu,ni,li,nj,lj,nk,lk,nl,ll
      real(dp) :: me_p1dp, me_p1d0, me_p1dm
      real(dp) :: me_p2dp, me_p2d0, me_p2dm
      real(dp) :: me_dpp, me_dp0, me_dpm, me_d00, me_d0m, me_dmm, me_del

      if(.not.init) call stop_not_init

      call zero_extfield_2bc(J)

      if (use_p) then
         ! P1 terms (treated same as ddY0)
         me_p1dp = MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu, J%k,+1, 1)
         me_p1d0 = MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu, J%k, 0, 1)
         me_p1dm = MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu, J%k,-1, 1)
         ! P2 terms (treated same as ddY0)
         me_p2dp = MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu, J%k,+1, 2)
         me_p2d0 = MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu, J%k, 0, 2)
         me_p2dm = MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu, J%k,-1, 2)
         ! p1-terms
         J%r0p = me_p1dm
         J%r0z = me_p1d0
         J%r0m = me_p1dp
         ! p2-terms
         J%rp0 = me_p2dm
         J%rz0 = me_p2d0
         J%rm0 = me_p2dp
      end if

      ! Matrix elements depend on the gaussian, so we must loop through terms
      ! and calculate matrix elements for each.
      select case (J%k)
      case default
          call abort(" ERROR, invalid K for extfield 2bc")
      !----------------------------------------------------------------------
      case(+1)

            ! 1st ddY0 terms
            me_dpm  = MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu,  +1,-1, 0)
            me_d00  = MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu,   0, 0, 0)
            me_dp0  = MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu,  +1, 0, 0)
            me_dpp  = MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu,  +1,+1, 0)

            ! 1-spin terms
            J%rkp = me_dpm
            J%rkz = me_dp0
            J%rkm = me_dpp

            ! 2-spin upper triangle
            J%rpz = me_d00
            J%rpm = me_dp0
            J%rzm = me_dpp

            ! 2nd ddY0 terms (only in 2-spin terms)
            J%rpz2 = me_dpm

      !----------------------------------------------------------------------
      case(0)

            ! 1st ddY0 terms
            me_dpm  = MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu,  +1,-1, 0)
            me_d00  = MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu,   0, 0, 0)
            me_d0m  = MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu,   0,-1, 0)
            me_dp0  = MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu,  +1, 0, 0)

            ! 1-spin terms
            J%rkp = me_d0m
            J%rkz = me_d00
            J%rkm = me_dp0

            ! 2-spin upper triangle
            J%rpz = me_d0m
            J%rpm = me_dpm
            J%rzm = me_dp0

            ! no 2nd ddY0 terms for 0th component

      !----------------------------------------------------------------------
      case(-1)

            ! 1st ddY0 terms
            me_dpm  = MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu,  +1,-1, 0)
            me_d00  = MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu,   0, 0, 0)
            me_d0m  = MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu,   0,-1, 0)
            me_dmm  = MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu,  -1,-1, 0)

            ! 1-spin terms
            J%rkp = me_dmm
            J%rkz = me_d0m
            J%rkm = me_dpm

            ! 2-spin upper triangle
            J%rpz = me_dmm
            J%rpm = me_d0m
            J%rzm = me_d00

            ! 2nd ddY0 terms (only in 2-spin terms)
            J%rzm2 = me_dpm

      end select

   end subroutine calc_Jr_opt

   !=========================================================================================
   ! Populate the spatial matrix elements needed for a given spherical component
   ! of the current. Calculates the Z component of the matrix elements.
   subroutine calc_Jz(J,g,ni,nj,nk,nl)

      implicit none
      type(extfield_2bc), intent(inout) :: J
      integer, intent(in) :: g, ni,nj,nk,nl
      real(dp) :: me_p1dp, me_p1d0, me_p1dm
      real(dp) :: me_p2dp, me_p2d0, me_p2dm
      real(dp) :: me_dpp, me_dp0, me_dpm, me_d00, me_d0m, me_dmm, me_del
      integer :: term, mu

      if(.not.init) call stop_not_init

      call zero_extfield_2bc(J)

      term = spatial_tm(g) ! 1st or 2nd term d1d2 Y0
      mu   = spatial_mu(g) ! gaussian index

      ! P1 and P2 terms (treated same as ddY0)
      ! (Only difference is the component of P=K. Do it here to save repeat code
      !  and to easily turn off P terms)
      if (term == 1 .and. use_p) then
         me_p1dp = MatrixElement_ddGz(ni,nj,nk,nl,mu, J%k,+1, 1)
         me_p1d0 = MatrixElement_ddGz(ni,nj,nk,nl,mu, J%k, 0, 1)
         me_p1dm = MatrixElement_ddGz(ni,nj,nk,nl,mu, J%k,-1, 1)
         me_p2dp = MatrixElement_ddGz(ni,nj,nk,nl,mu, J%k,+1, 2)
         me_p2d0 = MatrixElement_ddGz(ni,nj,nk,nl,mu, J%k, 0, 2)
         me_p2dm = MatrixElement_ddGz(ni,nj,nk,nl,mu, J%k,-1, 2)
         ! p1-terms
         J%r0p = +caux_cr2*me_p1dm
         J%r0z = +caux_2  *me_p1d0
         J%r0m = -caux_cr2*me_p1dp
         ! p2-terms
         J%rp0 = +caux_cr2*me_p2dm
         J%rz0 = +caux_2  *me_p2d0
         J%rm0 = -caux_cr2*me_p2dp
      end if

      ! Matrix elements depend on the gaussian, so we must loop through terms
      ! and re-calculate each matrix element for each gaussian.
         select case (J%k)
         case default
             call abort(" ERROR, invalid K for extfield 2bc")
         !----------------------------------------------------------------------
         case(+1)
            select case(term)
            case(0)
               ! Contact term
               me_del = MatrixElement_zZR(ni,nj,nk,nl)

               ! 1-spin and 2-spin terms
               J%rkp = -caux_d1_c4r2*me_del
               J%rpz = +cdel_c2r2*me_del
            case(1)

               ! 1st ddY0 terms
               me_dpm  = MatrixElement_ddGz(ni,nj,nk,nl,mu,  +1,-1, 0)
               me_d00  = MatrixElement_ddGz(ni,nj,nk,nl,mu,   0, 0, 0)
               me_dp0  = MatrixElement_ddGz(ni,nj,nk,nl,mu,  +1, 0, 0)
               me_dpp  = MatrixElement_ddGz(ni,nj,nk,nl,mu,  +1,+1, 0)

               ! 1-spin terms
               J%rkp = -caux_c3_c4r2*me_dpm
               J%rkz = -caux_c3_8   *me_dp0
               J%rkm = +caux_c3_c4r2*me_dpp

               ! 2-spin upper triangle
               J%rpz = +cY0p_c2r2*me_d00
               J%rpm = -cY0p_2   *me_dp0
               J%rzm = -cY0p_c2r2*me_dpp

            case(2)
               ! 2nd ddY0 terms (only in 2-spin terms)
               me_dpm = MatrixElement_ddGz(ni,nj,nk,nl,mu,+1,-1, 0)

               J%rpz = -cY0p_c2r2*me_dpm
            end select
         !----------------------------------------------------------------------
         case(0)
            select case(term)
            case(0)
               ! Contact term
               me_del = MatrixElement_zZR(ni,nj,nk,nl)

               ! 1-spin and 2-spin terms
               J%rkz = +caux_d1_8*me_del
               J%rpm = +cdel_2*me_del
            case(1)

               ! 1st ddY0 terms
               me_dpm  = MatrixElement_ddGz(ni,nj,nk,nl,mu,  +1,-1, 0)
               me_d00  = MatrixElement_ddGz(ni,nj,nk,nl,mu,   0, 0, 0)
               me_d0m  = MatrixElement_ddGz(ni,nj,nk,nl,mu,   0,-1, 0)
               me_dp0  = MatrixElement_ddGz(ni,nj,nk,nl,mu,  +1, 0, 0)
               ! 1-spin terms
               J%rkp = -caux_c3_c4r2*me_d0m
               J%rkz = -caux_c3_8   *me_d00
               J%rkm = +caux_c3_c4r2*me_dp0

               ! 2-spin upper triangle
               J%rpz = +cY0p_c2r2*me_d0m
               J%rpm = -cY0p_4   *me_dpm
               J%rzm = -cY0p_c2r2*me_dp0

            case(2)
                ! no 2nd ddY0 terms for 0th component

            end select
         !----------------------------------------------------------------------
         case(-1)
            select case(term)
            case(0)
               ! Contact term
               me_del = MatrixElement_zZR(ni,nj,nk,nl)

               ! 1-spin and 2-spin terms
               J%rkm = +caux_d1_c4r2*me_del
               J%rzm = -cdel_c2r2*me_del
            case(1)

               ! 1st ddY0 terms
               me_dpm  = MatrixElement_ddGz(ni,nj,nk,nl,mu,  +1,-1, 0)
               me_d00  = MatrixElement_ddGz(ni,nj,nk,nl,mu,   0, 0, 0)
               me_d0m  = MatrixElement_ddGz(ni,nj,nk,nl,mu,   0,-1, 0)
               me_dmm  = MatrixElement_ddGz(ni,nj,nk,nl,mu,  -1,-1, 0)

               ! 1-spin terms
               J%rkp = -caux_c3_c4r2*me_dmm
               J%rkz = -caux_c3_8   *me_d0m
               J%rkm = +caux_c3_c4r2*me_dpm

               ! 2-spin upper triangle
               J%rpz = +cY0p_c2r2*me_dmm
               J%rpm = -cY0p_2   *me_d0m
               J%rzm = -cY0p_c2r2*me_d00

            case(2)
               ! 2nd ddY0 terms (only in 2-spin terms)
               me_dpm = MatrixElement_ddGz(ni,nj,nk,nl,mu,+1,-1, 0)

               J%rzm = +cY0p_c2r2*me_dpm
            end select
         end select

   end subroutine calc_Jz

   ! Populate the spatial matrix elements needed for a given spherical component
   ! of the current. Calculates the radial component of the matrix elements.
   ! NOTE: prefactors are absent here!!!, Jz holds all prefactors, Jr holds just
   ! matrix elements such that the full thing is J_abcd = (pf*Jz)*(Jr)
   subroutine calc_Jr(J,g,ni,li,nj,lj,nk,lk,nl,ll)

      implicit none
      type(extfield_2bc), intent(inout) :: J
      integer, intent(in) :: g,ni,li,nj,lj,nk,lk,nl,ll
      real(dp) :: me_p1dp, me_p1d0, me_p1dm
      real(dp) :: me_p2dp, me_p2d0, me_p2dm
      real(dp) :: me_dpp, me_dp0, me_dpm, me_d00, me_d0m, me_dmm, me_del
      integer :: term, mu

      if(.not.init) call stop_not_init

      call zero_extfield_2bc(J)

      term = spatial_tm(g) ! 1st or 2nd term d1d2 Y0
      mu   = spatial_mu(g) ! gaussian index

      if (term == 1 .and. use_p) then
         ! P1 terms (treated same as ddY0)
         me_p1dp = MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu, J%k,+1, 1)
         me_p1d0 = MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu, J%k, 0, 1)
         me_p1dm = MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu, J%k,-1, 1)
         ! P2 terms (treated same as ddY0)
         me_p2dp = MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu, J%k,+1, 2)
         me_p2d0 = MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu, J%k, 0, 2)
         me_p2dm = MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu, J%k,-1, 2)
         ! p1-terms
         J%r0p = me_p1dm
         J%r0z = me_p1d0
         J%r0m = me_p1dp
         ! p2-terms
         J%rp0 = me_p2dm
         J%rz0 = me_p2d0
         J%rm0 = me_p2dp
      end if

      ! Matrix elements depend on the gaussian, so we must loop through terms
      ! and calculate matrix elements for each.
         select case (J%k)
         case default
             call abort(" ERROR, invalid K for extfield 2bc")
         !----------------------------------------------------------------------
         case(+1)
            select case(term)
            case(0)
               ! Contact term
               me_del = MatrixElement_rZR(ni,li,nj,lj,nk,lk,nl,ll)

               ! 1-spin and 2-spin terms
               J%rkp = me_del
               J%rpz = me_del
            case(1)

               ! 1st ddY0 terms
               me_dpm  = MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu,  +1,-1, 0)
               me_d00  = MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu,   0, 0, 0)
               me_dp0  = MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu,  +1, 0, 0)
               me_dpp  = MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu,  +1,+1, 0)

               ! 1-spin terms
               J%rkp = me_dpm
               J%rkz = me_dp0
               J%rkm = me_dpp

               ! 2-spin upper triangle
               J%rpz = me_d00
               J%rpm = me_dp0
               J%rzm = me_dpp

            case(2)
               ! 2nd ddY0 terms (only in 2-spin terms)
               me_dpm = MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu,+1,-1, 0)

               J%rpz = me_dpm
            end select
         !----------------------------------------------------------------------
         case(0)
            select case(term)
            case(0)
               ! Contact term
               me_del = MatrixElement_rZR(ni,li,nj,lj,nk,lk,nl,ll)

               ! 1-spin and 2-spin terms
               J%rkz = me_del
               J%rpm = me_del
            case(1)

               ! 1st ddY0 terms
               me_dpm  = MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu,  +1,-1, 0)
               me_d00  = MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu,   0, 0, 0)
               me_d0m  = MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu,   0,-1, 0)
               me_dp0  = MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu,  +1, 0, 0)

               ! 1-spin terms
               J%rkp = me_d0m
               J%rkz = me_d00
               J%rkm = me_dp0

               ! 2-spin upper triangle
               J%rpz = me_d0m
               J%rpm = me_dpm
               J%rzm = me_dp0

            case(2)
                ! no 2nd ddY0 terms for 0th component

            end select
         !----------------------------------------------------------------------
         case(-1)
            select case(term)
            case(0)
               ! Contact term
               me_del = MatrixElement_rZR(ni,li,nj,lj,nk,lk,nl,ll)

               ! 1-spin and 2-spin terms
               J%rkm = me_del
               J%rzm = me_del
            case(1)

               ! 1st ddY0 terms
               me_dpm  = MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu,  +1,-1, 0)
               me_d00  = MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu,   0, 0, 0)
               me_d0m  = MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu,   0,-1, 0)
               me_dmm  = MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu,  -1,-1, 0)

               ! 1-spin terms
               J%rkp = me_dmm
               J%rkz = me_d0m
               J%rkm = me_dpm

               ! 2-spin upper triangle
               J%rpz = me_dmm
               J%rpm = me_d0m
               J%rzm = me_d00

            case(2)
               ! 2nd ddY0 terms (only in 2-spin terms)
               me_dpm = MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu,+1,-1, 0)

               J%rzm = me_dpm
            end select
         end select

   end subroutine calc_Jr

   !================================================================================
   function calc_gam(Jd, Je, sa, sc, sd, sb, rhop, rhon, bminus) result(gam)
      implicit none

      type(extfield_2bc), intent(in) :: Jd, Je
      real(dp), intent(in) :: rhop, rhon
      integer, intent(in) :: sa, sc, sd, sb
      logical, intent(in) :: bminus
      real(dp) :: gam

      real(dp) :: Jrsa_d, Jrsx_d, Jrsa_e, Jrsb_e, Jrsx_e

      Jrsa_d = gam_Jrsa_dir(Jd, sa, sc, sd, sb)
      Jrsx_d = gam_Jrsx_dir(Jd, sa, sc, sd, sb) + gam_Jrsxp_dir(Jd, sa, sc, sd, sb)
      Jrsa_e = gam_Jrsa_exc(Je, sa, sc, sd, sb)
      Jrsb_e = gam_Jrsb_exc(Je, sa, sc, sd, sb)
      Jrsx_e = gam_Jrsx_exc(Je, sa, sc, sd, sb) + gam_Jrsxp_exc(Je, sa, sc, sd, sb)

      gam = 0
      if (bminus) then
         gam = ((Jrsa_d + h*Jrsx_d) - (Jrsa_e + h*Jrsx_e))*rhon &
             + ((Jrsa_d - h*Jrsx_d) - (Jrsb_e + h*Jrsx_e))*rhop
      else
         gam = ((Jrsa_d - h*Jrsx_d) - (Jrsb_e + h*Jrsx_e))*rhon &
             + ((Jrsa_d + h*Jrsx_d) - (Jrsa_e + h*Jrsx_e))*rhop
      end if

   end function

   function calc_del(J, sa, sc, sd, sb, kapp, kapn, bminus) result(del)
      implicit none

      type(extfield_2bc), intent(in) :: J
      real(dp), intent(in) :: kapp, kapn
      integer, intent(in) :: sa, sc, sd, sb
      logical, intent(in) :: bminus
      real(dp) :: del

      real(dp) :: Jrsa, Jrsb, Jrsx

      Jrsa = del_Jrsa(J, sa, sc, sd, sb)
      Jrsb = del_Jrsb(J, sa, sc, sd, sb)
      Jrsx = del_Jrsx(J, sa, sc, sd, sb)

      del = 0
      if (bminus) then
         del = (Jrsa + h*Jrsx)*kapn + (Jrsb + h*Jrsx)*kapp
      else
         del = (Jrsa - h*Jrsx)*kapp + (Jrsb - h*Jrsx)*kapn
      end if

   end function

   function calc_gam_sep(Jdn,Jdp,Jen,Jep, sa, sc, sd, sb, bminus) result(gam)
      implicit none

      type(extfield_2bc), intent(in) :: Jdn,Jdp, Jen,Jep
      integer, intent(in) :: sa, sc, sd, sb
      logical, intent(in) :: bminus
      real(dp) :: gam(6), gam1

      type(extfield_2bc) :: Jdn1,Jdp1, Jen1,Jep1
      real(dp) :: Jrsa_dn, Jrsx_dn, Jrsa_en, Jrsb_en, Jrsx_en
      real(dp) :: Jrsa_dp, Jrsx_dp, Jrsa_ep, Jrsb_ep, Jrsx_ep
      real(dp) :: Jrsxp_dn, Jrsxp_dp, Jrsxp_en, Jrsxp_ep

      ! Create a copy
      Jdn1 = Jdn
      Jdp1 = Jdp
      Jen1 = Jen
      Jep1 = Jep

      ! Add in 2nd didjY0 term
      Jdn1%rpz = Jdn1%rpz + Jdn1%rpz2
      Jdp1%rpz = Jdp1%rpz + Jdp1%rpz2
      Jen1%rpz = Jen1%rpz + Jen1%rpz2
      Jep1%rpz = Jep1%rpz + Jep1%rpz2

      Jdn1%rzm = Jdn1%rzm + Jdn1%rzm2
      Jdp1%rzm = Jdp1%rzm + Jdp1%rzm2
      Jen1%rzm = Jen1%rzm + Jen1%rzm2
      Jep1%rzm = Jep1%rzm + Jep1%rzm2

      ! Compute
      Jrsa_dn = gam_Jrsa_dir(Jdn1, sa, sc, sd, sb)
      Jrsx_dn = gam_Jrsx_dir(Jdn1, sa, sc, sd, sb)
      Jrsxp_dn = gam_Jrsxp_dir(Jdn1, sa, sc, sd, sb)

      Jrsa_en = gam_Jrsa_exc(Jen1, sa, sc, sd, sb)
      Jrsb_en = gam_Jrsb_exc(Jen1, sa, sc, sd, sb)
      Jrsx_en = gam_Jrsx_exc(Jen1, sa, sc, sd, sb)
      Jrsxp_en = gam_Jrsxp_exc(Jen1, sa, sc, sd, sb)

      Jrsa_dp = gam_Jrsa_dir(Jdp1, sa, sc, sd, sb)
      Jrsx_dp = gam_Jrsx_dir(Jdp1, sa, sc, sd, sb)
      Jrsxp_dp = gam_Jrsxp_dir(Jdp1, sa, sc, sd, sb)

      Jrsa_ep = gam_Jrsa_exc(Jep1, sa, sc, sd, sb)
      Jrsb_ep = gam_Jrsb_exc(Jep1, sa, sc, sd, sb)
      Jrsx_ep = gam_Jrsx_exc(Jep1, sa, sc, sd, sb)
      Jrsxp_ep = gam_Jrsxp_exc(Jep1, sa, sc, sd, sb)

      gam = 0
      if (bminus) then
         gam(1) = +Jrsa_dn + Jrsa_dp ! c3d
         gam(2) = -Jrsa_en - Jrsb_ep ! c3e

         gam(3) = +h*(Jrsx_dn - Jrsx_dp)   ! c4d
         gam(4) = -h*(Jrsx_en + Jrsx_ep)   ! c4e
         gam(5) = +h*(Jrsxp_dn - Jrsxp_dp) ! cpd
         gam(6) = -h*(Jrsxp_en + Jrsxp_ep) ! cpe

         !gam1 = ((Jrsa_dn + h*(Jrsx_dn+Jrsxp_dn)) - (Jrsa_en + h*(Jrsx_en+Jrsxp_en))) &
         !     + ((Jrsa_dp - h*(Jrsx_dp+Jrsxp_dp)) - (Jrsb_ep + h*(Jrsx_ep+Jrsxp_ep)))
      else
         gam(1) = +Jrsa_dn + Jrsa_dp ! c3d
         gam(2) = -Jrsb_en - Jrsa_ep ! c3e

         gam(3) = -h*(Jrsx_dn - Jrsx_dp)   ! c4d
         gam(4) = -h*(Jrsx_en + Jrsx_ep)   ! c4e
         gam(5) = -h*(Jrsxp_dn - Jrsxp_dp) ! cpd
         gam(6) = -h*(Jrsxp_en + Jrsxp_ep) ! cpe

         !gam1 = ((Jrsa_dn - h*(Jrsx_dn+Jrsxp_dn)) - (Jrsb_en + h*(Jrsx_en+Jrsxp_en))) &
         !     + ((Jrsa_dp + h*(Jrsx_dp+Jrsxp_dp)) - (Jrsa_ep + h*(Jrsx_ep+Jrsxp_ep)))
      end if
      

   end function

   function calc_del_sep(Jn,Jp, sa, sc, sd, sb, bminus) result(del)
      implicit none

      type(extfield_2bc), intent(in) :: Jn,Jp
      integer, intent(in) :: sa, sc, sd, sb
      logical, intent(in) :: bminus
      real(dp) :: del

      type(extfield_2bc) :: Jn1,Jp1
      real(dp) :: Jrsan, Jrsbn, Jrsxn
      real(dp) :: Jrsap, Jrsbp, Jrsxp

      ! Create a copy
      Jn1 = Jn
      Jp1 = Jp

      ! Add in 2nd didjY0 term
      Jn1%rpz = Jn1%rpz + Jn1%rpz2
      Jp1%rpz = Jp1%rpz + Jp1%rpz2

      Jn1%rzm = Jn1%rzm + Jn1%rzm2
      Jp1%rzm = Jp1%rzm + Jp1%rzm2

      ! Compute
      Jrsan = del_Jrsa(Jn1, sa, sc, sd, sb)
      Jrsbn = del_Jrsb(Jn1, sa, sc, sd, sb)
      Jrsxn = del_Jrsx(Jn1, sa, sc, sd, sb)

      Jrsap = del_Jrsa(Jp1, sa, sc, sd, sb)
      Jrsbp = del_Jrsb(Jp1, sa, sc, sd, sb)
      Jrsxp = del_Jrsx(Jp1, sa, sc, sd, sb)

      del = 0
      if (bminus) then
         del = (Jrsan + h*Jrsxn) + (Jrsbp + h*Jrsxp)
      else
         del = (Jrsap - h*Jrsxp) + (Jrsbn - h*Jrsxn)
      end if

   end function

   !================================================================================
   function gam_Jrsa_dir(J, sa, sc, sd, sb) result(Jrs)
      implicit none

      type(extfield_2bc), intent(in) :: J
      integer, intent(in) :: sa, sc, sd, sb

      real(dp) :: Jrs
      real(dp) :: rkp, rkz, rkm
      !real(dp) :: rpz, rpm, rzm
      !real(dp) :: r0p, r0z, r0m
      !real(dp) :: rp0, rz0, rm0

      if(.not.init) call stop_not_init

      ! Shorthand for a given spatial term
      rkp=J%rkp; rkz=J%rkz; rkm=J%rkm
      !rpz=J%rpz; rpm=J%rpm; rzm=J%rzm
      !r0p=J%r0p; r0z=J%r0z; r0m=J%r0m
      !rp0=J%rp0; rz0=J%rz0; rm0=J%rm0

      Jrs = 0

      ! Gam_uu, Gam_dd
      if (sa == sc) then
         ! rho_du = rho_ud
         if (sd+sb == 0) then
            Jrs = 0
         ! rho_dd, rho_uu
         else
            Jrs = (sa*h)*rkz ! sa = +/-1, need spin +/-(1/2)
         end if

      ! Gam_ud
      else if (sa == +1 .and. sc == -1) then
         ! rho_du = rho_ud
         if (sd+sb == 0) then
            Jrs = 0
         ! rho_uu, rho_dd
         else
            Jrs = rkp
         end if

      ! Gam_du
      else if (sa == -1 .and. sc == +1) then
         ! rho_du = rho_ud
         if (sd+sb == 0) then
            Jrs = 0
         ! rho_uu, rho_dd
         else
            Jrs = rkm
         end if

      end if

   end function

   function gam_Jrsb_dir(J, sa, sc, sd, sb) result(Jrs)
      implicit none

      type(extfield_2bc), intent(in) :: J
      integer, intent(in) :: sa, sc, sd, sb

      real(dp) :: Jrs
      real(dp) :: rkp, rkz, rkm
      !real(dp) :: rpz, rpm, rzm
      !real(dp) :: r0p, r0z, r0m
      !real(dp) :: rp0, rz0, rm0


      if(.not.init) call stop_not_init

      ! Shorthand for a given spatial term
      rkp=J%rkp; rkz=J%rkz; rkm=J%rkm
      !rpz=J%rpz; rpm=J%rpm; rzm=J%rzm
      !r0p=J%r0p; r0z=J%r0z; r0m=J%r0m
      !rp0=J%rp0; rz0=J%rz0; rm0=J%rm0

      Jrs = 0

      ! Gam_uu, Gam_dd
      if (sa == sc) then
         ! rho_db = du
         if (sd == -1 .and. sb == +1) then
            Jrs = rkp
         ! rho_db = ud
         else if (sd == +1 .and. sb == -1) then
            Jrs = rkm
         ! rho_dd, rho_uu
         else if (sd+sb > 0) then
            Jrs = +h*rkz
         else
            Jrs = -h*rkz
         end if

      ! Gam_ud, Gam_du
      else
         Jrs = 0
      end if

   end function

   function gam_Jrsx_dir(J, sa, sc, sd, sb) result(Jrs)
      implicit none

      type(extfield_2bc), intent(in) :: J
      integer, intent(in) :: sa, sc, sd, sb

      real(dp) :: Jrs
      !real(dp) :: rkp, rkz, rkm
      real(dp) :: rpz, rpm, rzm
      !real(dp) :: r0p, r0z, r0m
      !real(dp) :: rp0, rz0, rm0


      if(.not.init) call stop_not_init

      ! Shorthand for a given spatial term
      !rkp=J%rkp; rkz=J%rkz; rkm=J%rkm
      rpz=J%rpz; rpm=J%rpm; rzm=J%rzm
      !r0p=J%r0p; r0z=J%r0z; r0m=J%r0m
      !rp0=J%rp0; rz0=J%rz0; rm0=J%rm0

      Jrs = 0

      ! Gam_uu, Gam_dd
      if (sa == sc) then
         ! rho_du = rho_ud
         if (sd == -1 .and. sb == +1) then
            Jrs = (sa*h)*(-rpz)
         else if (sd == +1 .and. sb == -1) then
            Jrs = (sa*h)*(rzm)
         ! rho_uu
         else if (sd+sb > 0) then
            Jrs = 0
         ! rho_dd
         else
            Jrs = 0
         end if

      ! Gam_ud
      else if (sa == +1 .and. sc == -1) then
         ! rho_du = rho_ud
         if (sd == -1 .and. sb == +1) then
            Jrs = 0
         else if (sd == +1 .and. sb == -1) then
            Jrs = rpm
         ! rho_uu
         else if (sd+sb > 0) then
            Jrs = +h*rpz
         ! rho_dd
         else
            Jrs = -h*rpz
         end if

      ! Gam_du
      else if (sa == -1 .and. sc == +1) then
         ! rho_du = rho_ud
         if (sd == -1 .and. sb == +1) then
            Jrs = -rpm
         else if (sd == +1 .and. sb == -1) then
            Jrs = 0
         ! rho_uu
         else if (sd+sb > 0) then
            Jrs = -h*rzm
         ! rho_dd
         else
            Jrs = +h*rzm
         end if

      end if

   end function

   !--------------------------------------------------------------------------------
   function gam_Jrsxp_dir(J, sa, sc, sd, sb) result(Jrs)
      implicit none

      type(extfield_2bc), intent(in) :: J
      integer, intent(in) :: sa, sc, sd, sb

      real(dp) :: Jrs
      !real(dp) :: rkp, rkz, rkm
      !real(dp) :: rpz, rpm, rzm
      real(dp) :: r0p, r0z, r0m
      real(dp) :: rp0, rz0, rm0


      if(.not.init) call stop_not_init

      ! Shorthand for a given spatial term
      !rkp=J%rkp; rkz=J%rkz; rkm=J%rkm
      !rpz=J%rpz; rpm=J%rpm; rzm=J%rzm
      r0p=J%r0p; r0z=J%r0z; r0m=J%r0m
      rp0=J%rp0; rz0=J%rz0; rm0=J%rm0

      Jrs = 0

      ! Gam_uu, Gam_dd
      if (sa == sc) then
         ! rho_du = rho_ud
         if (sd == -1 .and. sb == +1) then
            Jrs = + r0p
         else if (sd == +1 .and. sb == -1) then
            Jrs = + r0m
         ! rho_uu
         else if (sd+sb > 0) then
            Jrs = (sa*h)*rz0 + h*r0z
         ! rho_dd
         else
            Jrs = (sa*h)*rz0 - h*r0z
         end if

      ! Gam_ud
      else if (sa == +1 .and. sc == -1) then
         ! rho_du = rho_ud
         if (sd == -1 .and. sb == +1) then
            Jrs = 0
         else if (sd == +1 .and. sb == -1) then
            Jrs = 0
         ! rho_uu
         else if (sd+sb > 0) then
            Jrs = + rp0
         ! rho_dd
         else
            Jrs = + rp0
         end if

      ! Gam_du
      else if (sa == -1 .and. sc == +1) then
         ! rho_du = rho_ud
         if (sd == -1 .and. sb == +1) then
            Jrs = 0
         else if (sd == +1 .and. sb == -1) then
            Jrs = 0
         ! rho_uu
         else if (sd+sb > 0) then
            Jrs = + rm0
         ! rho_dd
         else
            Jrs = + rm0
         end if

      end if

   end function

   !--------------------------------------------------------------------------------
   function gam_Jrsa_exc(J, sa, sc, sd, sb) result(Jrs)
      implicit none

      type(extfield_2bc), intent(in) :: J
      integer, intent(in) :: sa, sc, sd, sb

      real(dp) :: Jrs
      real(dp) :: rkp, rkz, rkm
      !real(dp) :: rpz, rpm, rzm
      !real(dp) :: r0p, r0z, r0m
      !real(dp) :: rp0, rz0, rm0


      if(.not.init) call stop_not_init

      ! Shorthand for a given spatial term
      rkp=J%rkp; rkz=J%rkz; rkm=J%rkm
      !rpz=J%rpz; rpm=J%rpm; rzm=J%rzm
      !r0p=J%r0p; r0z=J%r0z; r0m=J%r0m
      !rp0=J%rp0; rz0=J%rz0; rm0=J%rm0

      Jrs = 0

      ! Gam_uu
      if (sa == +1 .and. sc == +1) then
         ! rho_du = rho_ud
         if (sd == -1 .and. sb == +1) then
            Jrs = rkp
         else if (sd == +1 .and. sb == -1) then
            Jrs = 0
         ! rho_uu
         else if (sd+sb > 0) then
            Jrs = h*rkz
         ! rho_dd
         else
            Jrs = 0
         end if

      ! Gam_dd
      else if (sa == -1 .and. sc == -1) then
         ! rho_du = rho_ud
         if (sd == -1 .and. sb == +1) then
            Jrs = 0
         else if (sd == +1 .and. sb == -1) then
            Jrs = rkm
         ! rho_uu
         else if (sd+sb > 0) then
            Jrs = 0
         ! rho_dd
         else
            Jrs = -h*rkz
         end if

      ! Gam_ud
      else if (sa == +1 .and. sc == -1) then
         ! rho_du = rho_ud
         if (sd == -1 .and. sb == +1) then
            Jrs = 0
         else if (sd == +1 .and. sb == -1) then
            Jrs = h*rkz
         ! rho_uu
         else if (sd+sb > 0) then
            Jrs = 0
         ! rho_dd
         else
            Jrs = rkp
         end if

      ! Gam_du
      else if (sa == -1 .and. sc == +1) then
         ! rho_du = rho_ud
         if (sd == -1 .and. sb == +1) then
            Jrs = -h*rkz
         else if (sd == +1 .and. sb == -1) then
            Jrs = 0
         ! rho_uu
         else if (sd+sb > 0) then
            Jrs = rkm
         ! rho_dd
         else
            Jrs = 0
         end if

      end if

   end function

   function gam_Jrsb_exc(J, sa, sc, sd, sb) result(Jrs)
      implicit none

      type(extfield_2bc), intent(in) :: J
      integer, intent(in) :: sa, sc, sd, sb

      real(dp) :: Jrs
      real(dp) :: rkp, rkz, rkm
      !real(dp) :: rpz, rpm, rzm
      !real(dp) :: r0p, r0z, r0m
      !real(dp) :: rp0, rz0, rm0


      if(.not.init) call stop_not_init

      ! Shorthand for a given spatial term
      rkp=J%rkp; rkz=J%rkz; rkm=J%rkm
      !rpz=J%rpz; rpm=J%rpm; rzm=J%rzm
      !r0p=J%r0p; r0z=J%r0z; r0m=J%r0m
      !rp0=J%rp0; rz0=J%rz0; rm0=J%rm0

      Jrs = 0

      ! Gam_uu
      if (sa == +1 .and. sc == +1) then
         ! rho_du = rho_ud
         if (sd == -1 .and. sb == +1) then
            Jrs = 0
         else if (sd == +1 .and. sb == -1) then
            Jrs = rkm
         ! rho_uu
         else if (sd+sb > 0) then
            Jrs = h*rkz
         ! rho_dd
         else
            Jrs = 0
         end if

      ! Gam_dd
      else if (sa == -1 .and. sc == -1) then
         ! rho_du = rho_ud
         if (sd == -1 .and. sb == +1) then
            Jrs = rkp
         else if (sd == +1 .and. sb == -1) then
            Jrs = 0
         ! rho_uu
         else if (sd+sb > 0) then
            Jrs = 0
         ! rho_dd
         else
            Jrs = -h*rkz
         end if

      ! Gam_ud
      else if (sa == +1 .and. sc == -1) then
         ! rho_du = rho_ud
         if (sd == -1 .and. sb == +1) then
            Jrs = 0
         else if (sd == +1 .and. sb == -1) then
            Jrs = -h*rkz
         ! rho_uu
         else if (sd+sb > 0) then
            Jrs = rkp
         ! rho_dd
         else
            Jrs = 0
         end if

      ! Gam_du
      else if (sa == -1 .and. sc == +1) then
         ! rho_du = rho_ud
         if (sd == -1 .and. sb == +1) then
            Jrs = +h*rkz
         else if (sd == +1 .and. sb == -1) then
            Jrs = 0
         ! rho_uu
         else if (sd+sb > 0) then
            Jrs = 0
         ! rho_dd
         else
            Jrs = rkm
         end if

      end if

   end function

   function gam_Jrsx_exc(J, sa, sc, sd, sb) result(Jrs)
      implicit none

      type(extfield_2bc), intent(in) :: J
      integer, intent(in) :: sa, sc, sd, sb

      real(dp) :: Jrs
      !real(dp) :: rkp, rkz, rkm
      real(dp) :: rpz, rpm, rzm
      !real(dp) :: r0p, r0z, r0m
      !real(dp) :: rp0, rz0, rm0


      if(.not.init) call stop_not_init

      ! Shorthand for a given spatial term
      !rkp=J%rkp; rkz=J%rkz; rkm=J%rkm
      rpz=J%rpz; rpm=J%rpm; rzm=J%rzm
      !r0p=J%r0p; r0z=J%r0z; r0m=J%r0m
      !rp0=J%rp0; rz0=J%rz0; rm0=J%rm0

      Jrs = 0

      ! Gam_uu
      if (sa == +1 .and. sc == +1) then
         ! rho_du = rho_ud
         if (sd == -1 .and. sb == +1) then
            Jrs = h*rpz
         else if (sd == +1 .and. sb == -1) then
            Jrs = h*rzm
         ! rho_uu
         else if (sd+sb > 0) then
            Jrs = 0
         ! rho_dd
         else
            Jrs = rpm
         end if

      ! Gam_dd
      else if (sa == -1 .and. sc == -1) then
         ! rho_du = rho_ud
         if (sd == -1 .and. sb == +1) then
            Jrs = h*rpz
         else if (sd == +1 .and. sb == -1) then
            Jrs = h*rzm
         ! rho_uu
         else if (sd+sb > 0) then
            Jrs = -rpm
         ! rho_dd
         else
            Jrs = 0
         end if

      ! Gam_ud
      else if (sa == +1 .and. sc == -1) then
         ! rho_du = rho_ud
         if (sd == -1 .and. sb == +1) then
            Jrs = 0
         else if (sd == +1 .and. sb == -1) then
            Jrs = 0
         ! rho_uu
         else if (sd+sb > 0) then
            Jrs = -h*rpz
         ! rho_dd
         else
            Jrs = -h*rpz
         end if

      ! Gam_du
      else if (sa == -1 .and. sc == +1) then
         ! rho_du = rho_ud
         if (sd == -1 .and. sb == +1) then
            Jrs = 0
         else if (sd == +1 .and. sb == -1) then
            Jrs = 0
         ! rho_uu
         else if (sd+sb > 0) then
            Jrs = -h*rzm
         ! rho_dd
         else
            Jrs = -h*rzm
         end if

      end if

   end function

   !--------------------------------------------------------------------------------
   function gam_Jrsxp_exc(J, sa, sc, sd, sb) result(Jrs)
      implicit none

      type(extfield_2bc), intent(in) :: J
      integer, intent(in) :: sa, sc, sd, sb

      real(dp) :: Jrs
      !real(dp) :: rkp, rkz, rkm
      !real(dp) :: rpz, rpm, rzm
      real(dp) :: r0p, r0z, r0m
      real(dp) :: rp0, rz0, rm0


      if(.not.init) call stop_not_init

      ! Shorthand for a given spatial term
      !rkp=J%rkp; rkz=J%rkz; rkm=J%rkm
      !rpz=J%rpz; rpm=J%rpm; rzm=J%rzm
      r0p=J%r0p; r0z=J%r0z; r0m=J%r0m
      rp0=J%rp0; rz0=J%rz0; rm0=J%rm0

      Jrs = 0

      ! Gam_uu
      if (sa == +1 .and. sc == +1) then
         ! rho_du = rho_ud
         if (sd == -1 .and. sb == +1) then
            Jrs = + rp0
         else if (sd == +1 .and. sb == -1) then
            Jrs = + r0m
         ! rho_uu
         else if (sd+sb > 0) then
            Jrs = h*rz0 + h*r0z
         ! rho_dd
         else
            Jrs = 0
         end if

      ! Gam_dd
      else if (sa == -1 .and. sc == -1) then
         ! rho_du = rho_ud
         if (sd == -1 .and. sb == +1) then
            Jrs = + r0p
         else if (sd == +1 .and. sb == -1) then
            Jrs = + rm0
         ! rho_uu
         else if (sd+sb > 0) then
            Jrs = 0
         ! rho_dd
         else
            Jrs = -h*rz0 - h*r0z
         end if

      ! Gam_ud
      else if (sa == +1 .and. sc == -1) then
         ! rho_du = rho_ud
         if (sd == -1 .and. sb == +1) then
            Jrs = 0
         else if (sd == +1 .and. sb == -1) then
            Jrs = h*rz0 - h*r0z
         ! rho_uu
         else if (sd+sb > 0) then
            Jrs = + r0p
         ! rho_dd
         else
            Jrs = + rp0
         end if

      ! Gam_du
      else if (sa == -1 .and. sc == +1) then
         ! rho_du = rho_ud
         if (sd == -1 .and. sb == +1) then
            Jrs = -h*rz0 + h*r0z
         else if (sd == +1 .and. sb == -1) then
            Jrs = 0
         ! rho_uu
         else if (sd+sb > 0) then
            Jrs = + rm0
         ! rho_dd
         else
            Jrs = + r0m
         end if

      end if

   end function

   !--------------------------------------------------------------------------------
   function del_Jrsa(J, sa, sc, sd, sb) result(Jrs)
      implicit none

      type(extfield_2bc), intent(in) :: J
      integer, intent(in) :: sa, sc, sd, sb

      real(dp) :: Jrs
      real(dp) :: rkp, rkz, rkm
      !real(dp) :: rpz, rpm, rzm
      !real(dp) :: r0p, r0z, r0m
      !real(dp) :: rp0, rz0, rm0


      if(.not.init) call stop_not_init

      ! Shorthand for a given spatial term
      rkp=J%rkp; rkz=J%rkz; rkm=J%rkm
      !rpz=J%rpz; rpm=J%rpm; rzm=J%rzm
      !r0p=J%r0p; r0z=J%r0z; r0m=J%r0m
      !rp0=J%rp0; rz0=J%rz0; rm0=J%rm0

      Jrs = 0

      ! Gam_uu
      if (sa == +1 .and. sc == +1) then
         ! rho_du = rho_ud
         if (sd == -1 .and. sb == +1) then
            Jrs = 0
         else if (sd == +1 .and. sb == -1) then
            Jrs = rkp
         ! rho_uu
         else if (sd+sb > 0) then
            Jrs = h*rkz
         ! rho_dd
         else
            Jrs = 0
         end if

      ! Gam_dd
      else if (sa == -1 .and. sc == -1) then
         ! rho_du = rho_ud
         if (sd == -1 .and. sb == +1) then
            Jrs = rkm
         else if (sd == +1 .and. sb == -1) then
            Jrs = 0
         ! rho_uu
         else if (sd+sb > 0) then
            Jrs = 0
         ! rho_dd
         else
            Jrs = -h*rkz
         end if

      ! Gam_ud
      else if (sa == +1 .and. sc == -1) then
         ! rho_du = rho_ud
         if (sd == -1 .and. sb == +1) then
            Jrs = h*rkz
         else if (sd == +1 .and. sb == -1) then
            Jrs = 0
         ! rho_uu
         else if (sd+sb > 0) then
            Jrs = 0
         ! rho_dd
         else
            Jrs = rkp
         end if

      ! Gam_du
      else if (sa == -1 .and. sc == +1) then
         ! rho_du = rho_ud
         if (sd == -1 .and. sb == +1) then
            Jrs = 0
         else if (sd == +1 .and. sb == -1) then
            Jrs = -h*rkz
         ! rho_uu
         else if (sd+sb > 0) then
            Jrs = rkm
         ! rho_dd
         else
            Jrs = 0
         end if

      end if

   end function

   function del_Jrsb(J, sa, sc, sd, sb) result(Jrs)
      implicit none

      type(extfield_2bc), intent(in) :: J
      integer, intent(in) :: sa, sc, sd, sb

      real(dp) :: Jrs
      real(dp) :: rkp, rkz, rkm
      !real(dp) :: rpz, rpm, rzm
      !real(dp) :: r0p, r0z, r0m
      !real(dp) :: rp0, rz0, rm0


      if(.not.init) call stop_not_init

      ! Shorthand for a given spatial term
      rkp=J%rkp; rkz=J%rkz; rkm=J%rkm
      !rpz=J%rpz; rpm=J%rpm; rzm=J%rzm
      !r0p=J%r0p; r0z=J%r0z; r0m=J%r0m
      !rp0=J%rp0; rz0=J%rz0; rm0=J%rm0

      Jrs = 0

      ! Gam_uu
      if (sa == +1 .and. sc == +1) then
         ! rho_du = rho_ud
         if (sd == -1 .and. sb == +1) then
            Jrs = -rkm
         else if (sd == +1 .and. sb == -1) then
            Jrs = 0
         ! rho_uu
         else if (sd+sb > 0) then
            Jrs = -h*rkz
         ! rho_dd
         else
            Jrs = 0
         end if

      ! Gam_dd
      else if (sa == -1 .and. sc == -1) then
         ! rho_du = rho_ud
         if (sd == -1 .and. sb == +1) then
            Jrs = 0
         else if (sd == +1 .and. sb == -1) then
            Jrs = -rkp
         ! rho_uu
         else if (sd+sb > 0) then
            Jrs = 0
         ! rho_dd
         else
            Jrs = +h*rkz
         end if

      ! Gam_ud
      else if (sa == +1 .and. sc == -1) then
         ! rho_du = rho_ud
         if (sd == -1 .and. sb == +1) then
            Jrs = +h*rkz
         else if (sd == +1 .and. sb == -1) then
            Jrs = 0
         ! rho_uu
         else if (sd+sb > 0) then
            Jrs = -rkp
         ! rho_dd
         else
            Jrs = 0
         end if

      ! Gam_du
      else if (sa == -1 .and. sc == +1) then
         ! rho_du = rho_ud
         if (sd == -1 .and. sb == +1) then
            Jrs = 0
         else if (sd == +1 .and. sb == -1) then
            Jrs = -h*rkz
         ! rho_uu
         else if (sd+sb > 0) then
            Jrs = 0
         ! rho_dd
         else
            Jrs = -rkm
         end if

      end if

   end function

   function del_Jrsx(J, sa, sc, sd, sb) result(Jrs)
      implicit none

      type(extfield_2bc), intent(in) :: J
      integer, intent(in) :: sa, sc, sd, sb

      real(dp) :: Jrs
      !real(dp) :: rkp, rkz, rkm
      real(dp) :: rpz, rpm, rzm
      real(dp) :: r0p, r0z, r0m
      real(dp) :: rp0, rz0, rm0


      if(.not.init) call stop_not_init

      ! Shorthand for a given spatial term
      !rkp=J%rkp; rkz=J%rkz; rkm=J%rkm
      rpz=J%rpz; rpm=J%rpm; rzm=J%rzm
      r0p=J%r0p; r0z=J%r0z; r0m=J%r0m
      rp0=J%rp0; rz0=J%rz0; rm0=J%rm0

      Jrs = 0

      ! Gam_uu
      if (sa == +1 .and. sc == +1) then
         ! rho_du = rho_ud
         if (sd == -1 .and. sb == +1) then
            Jrs = -h*rzm - r0m
         else if (sd == +1 .and. sb == -1) then
            Jrs = -h*rpz + rp0
         ! rho_uu
         else if (sd+sb > 0) then
            Jrs = h*rz0 - h*r0z
         ! rho_dd
         else
            Jrs = -rpm
         end if

      ! Gam_dd
      else if (sa == -1 .and. sc == -1) then
         ! rho_du = rho_ud
         if (sd == -1 .and. sb == +1) then
            Jrs = -h*rzm + rm0
         else if (sd == +1 .and. sb == -1) then
            Jrs = -h*rpz - r0p
         ! rho_uu
         else if (sd+sb > 0) then
            Jrs = +rpm
         ! rho_dd
         else
            Jrs = -h*rz0 + h*r0z
         end if

      ! Gam_ud
      else if (sa == +1 .and. sc == -1) then
         ! rho_du = rho_ud
         if (sd == -1 .and. sb == +1) then
            Jrs = h*rz0 + h*r0z
         else if (sd == +1 .and. sb == -1) then
            Jrs = 0
         ! rho_uu
         else if (sd+sb > 0) then
            Jrs = +h*rpz - r0p
         ! rho_dd
         else
            Jrs = +h*rpz + rp0
         end if

      ! Gam_du
      else if (sa == -1 .and. sc == +1) then
         ! rho_du = rho_ud
         if (sd == -1 .and. sb == +1) then
            Jrs = 0
         else if (sd == +1 .and. sb == -1) then
            Jrs = -h*rz0 - h*r0z
         ! rho_uu
         else if (sd+sb > 0) then
            Jrs = +h*rzm + rm0
         ! rho_dd
         else
            Jrs = +h*rzm - r0m
         end if

      end if

   end function

   !--------------------------------------------------------------------------------

   subroutine print_extfield_2bc(J)
      implicit none
      type(extfield_2bc), intent(in) :: J

      if(.not.init) return

      write(*,'(1x,"rkp=",es20.10,", rkz=",es20.10,", rkm=",es20.10)') J%rkp, J%rkz, J%rkm
      write(*,'(1x,"rpz=",es20.10,", rpm=",es20.10,", rzm=",es20.10)') J%rpz, J%rpm, J%rzm
      write(*,'(1x,"r0p=",es20.10,", r0z=",es20.10,", r0m=",es20.10)') J%r0p, J%r0z, J%r0m
      write(*,'(1x,"rp0=",es20.10,", rz0=",es20.10,", rm0=",es20.10)') J%rp0, J%rz0, J%rm0
      write(*,*)

   end subroutine

   function extfield_2bc_is_zero(J) result(zero)
      implicit none
      type(extfield_2bc), intent(in) :: J

      real(dp), parameter :: tol=1e-20
      logical :: zero
      logical :: rkp, rkz, rkm
      logical :: rpz, rpm, rzm
      logical :: r0p, r0z, r0m
      logical :: rp0, rz0, rm0

      rkp = (abs(J%rkp)<tol)
      rpz = (abs(J%rpz)<tol)
      r0p = (abs(J%r0p)<tol)
      rp0 = (abs(J%rp0)<tol)
      rkz = (abs(J%rkz)<tol)
      rpm = (abs(J%rpm)<tol)
      r0z = (abs(J%r0z)<tol)
      rz0 = (abs(J%rz0)<tol)
      rkm = (abs(J%rkm)<tol)
      rzm = (abs(J%rzm)<tol)
      r0m = (abs(J%r0m)<tol)
      rm0 = (abs(J%rm0)<tol)

      zero = (rkp.and.rkz.and.rkm.and.&
              rpz.and.rpm.and.rzm.and.&
              r0p.and.r0z.and.r0m.and.&
              rp0.and.rz0.and.rm0)

   end function

   !!======================================================================
   !!>  Calculates the the necessary two-body potential matrix elements
   !!>  \f$\langle n_{z_i} n_{z_j}|V_z|n_{z_k} n_{z_l} \rangle \f$.
   !!>
   !!>  If the preprocessor variable GOGNY_SYMMETRIES is set to 1 only
   !!>  states with \f$ n_{z_k} \geq n_{z_i} \f$,
   !!>  \f$ n_{z_l} \geq n_{z_j} \f$ are calculated and stored. In any
   !!>  case, states where \f$ n_{z_i}+n_{z_j}+n_{z_k}+n_{z_l} \f$ is not
   !!>  an even number are not calculated since those matrix elements are
   !!>  zero.
   !!======================================================================
   !subroutine CalculateMEJz(op)
   !  use type_extfield
   !  use hfb_solution, only : nzx
   !  implicit none
   !  type(external_field), intent(in) :: op
   !  real(dp) :: Vz1,Vz2,Vz
   !  integer(ipr) :: ii,N,nzi,nzj,nzk,nzl,ig
   !  type(extfield_2bc) :: zero_2bc, Jz


   !  if(.not.init) call stop_not_init

   !  call zero_extfield_2bc(zero_2bc)
   !  call zero_extfield_2bc(Jz)
   !  zero_2bc%k = op%k
   !  Jz%k = op%k

   !  call calculate_Zblock()

   !  n = nzx
   !  ii = Zblock(n,n,n,n); NumVz = ii
   !  if(allocated(MEJz)) deallocate(MEJz)
   !  allocate(MEJz(1:nspatial,0:ii))
   !  MEJz(1:nspatial,0) = zero_2bc

   !  do nzi = 0,n
   !     do nzj = 0,n
   !         do nzk = 0,n
   !           do nzl = 0,n
   !              ! 1) Exchange (a,c) or (b,d) are equal except K=0, the P terms change sign
   !              !    To simplify life, don't use echange symmetry
   !              ! 2) Non-zero MEs for i+j+k+l + (d=0,1,2) even
   !              ii = Zblock(nzi,nzj,nzk,nzl)
   !              do ig = 1,nspatial
   !                 call calc_Jz(Jz,ig, nzi,nzj,nzk,nzl)
   !                 MEJz(ig,ii) = Jz
   !              enddo
   !           enddo
   !        enddo
   !     enddo
   !  enddo

   !end subroutine CalculateMEJz


   !!======================================================================
   !!>  Calculates the necessary two-body potential
   !!>  matrix elements \f$ \langle n_{r_i} \Lambda_i, n_{r_j} \Lambda_j
   !!>  |V_p | n_{r_k} \Lambda_k, n_{r_l} \Lambda_l \rangle\f$.
   !!>
   !!>  If the preprocessor variable GOGNY_SYMMETRIES is set to 1 only
   !!>  the matrix elements that will be used in the calculation of the HFB
   !!>  fields are calculated. In any case, states where \f$ -\Lambda_i
   !!>  - \Lambda_j + \Lambda_k + \Lambda_l \neq 0 \f$ are not calculated
   !!>  since the  matrix element is zero.
   !!======================================================================
   !subroutine CalculateMEJr(op)
   !  use type_extfield
   !  use hfb_solution
   !  implicit none
   !  type(external_field), intent(in) :: op
   !  real(dp) :: Vr
   !  integer(ipr) :: ii,nri,li,nrj,lj,nrk,lk,nrl,ll,ig,im,N
   !  integer(ipr) :: nza,nzb,nzc,nzd,nra,nrb,nrc,nrd,oa,ob,la,lb,lc,ld
   !  integer(ipr) :: nla,nlb,nlc,nld,nsa,nsc,nsac,j1,j2,nsb,nsd,nsdb
   !  integer(ipr) :: ir_abdc,ir_abcd,ir_acbd,jr_abdc,jr_abcd,jr_acbd
   !  integer(ipr) :: ita,itb,itc,itd,iba,ibc,ibd,ibb,ir
   !  integer(ipr) , allocatable :: index_flag(:), nr_flag(:,:), &
   !       nl_flag(:,:),ir_flag(:)
   !  integer(ipr) :: ilauf
   !  type(extfield_2bc) :: zero_2bc, Jr
   !  logical :: gblock
   !  integer(ipr) :: ita_tr, itc_tr

   !   if(.not.init) call stop_not_init

   !   call calculate_Nblock

   !  call zero_extfield_2bc(zero_2bc)
   !  call zero_extfield_2bc(Jr)
   !  zero_2bc%k = op%k
   !  Jr%k = op%k

   !  N  = max(2*nrx,nlx)
   !  im = mod(N,2)

   !  ii = Nblock(n/2,n/2,n/2,n/2,im,-im)+L_block(n/2,n/2,im,-im,-im,n)
   !  NumVr = ii

   !  if(allocated(MEJr)) deallocate(MEJr)
   !  if(.not.allocated(MEJr)) then
   !     allocate(MEJr(1:nspatial,0:ii))
   !     allocate(index_flag(0:ii))
   !     allocate(nr_flag(1:4,0:ii))
   !     allocate(nl_flag(1:4,0:ii))
   !     allocate(ir_flag(0:ii))
   !  endif

   !  call zero_extfield_2bc(zero_2bc)
   !  index_flag = 0
   !  MEJr(1:nspatial,0) = zero_2bc

   !  !We go over the loop that the field calculation will go through in
   !  !order to identify the matrix elements that will be necesary and
   !  !calculate only those on a subsequent loop.
   !  ilauf = 0
   !  ! Loop over just za,zc in [0,1]??? WHY? To get both parities (-1)^{z+l)?
   !  do ita_tr = 1,ntx*2
   !     if (ita_tr > ntx) then
   !        ita = ita_tr - ntx
   !        nra = nr(ita); nza = nz(ita); nla = -nl(ita); nsa = -ns(ita)
   !     else
   !        ita = ita_tr
   !        nra = nr(ita); nza = nz(ita); nla = nl(ita); nsa = ns(ita)
   !     end if
   !     if(nza.gt.1) cycle
   !     iba = ib_itx(ita_tr)
   !     do itc_tr = 1,ntx*2
   !        if (itc_tr > ntx) then
   !           itc = itc_tr - ntx
   !           nrc = nr(itc); nzc = nz(itc); nlc = -nl(itc); nsc = -ns(itc)
   !        else
   !           itc = itc_tr
   !           nrc = nr(itc); nzc = nz(itc); nlc = nl(itc); nsc = ns(itc)
   !        end if
   !        if(nzc.gt.1) cycle
   !        ibc = ib_itx(itc_tr)

   !        ! Are we in a non-zero block for Jeff?
   !        gblock = (op%mat%ir2c(iba) == ibc)
   !        if (.not.gblock) cycle

   !        do itb = 1,nttx
   !           nrb=nrr(itb); nlb=nll(itb); nsb=nss(itb); ibb=noo(itb)
   !           do itd = 1,nttx
   !              nrd=nrr(itd); nld=nll(itd); nsd=nss(itd); ibd=noo(itd)

   !              ! Are we in a non-zero block for densities?
   !              if(ibb.ne.ibd) cycle

   !              if(nla+nlb.ne.nlc+nld+op%k+1&
   !                 .and.nla+nlb.ne.nlc+nld+op%k&
   !                 .and.nla+nlb.ne.nlc+nld+op%k-1) cycle
   !              if(nla-nlb.ne.nlc-nld+op%k+1&
   !                 .and.nla-nlb.ne.nlc-nld+op%k&
   !                 .and.nla-nlb.ne.nlc-nld+op%k-1) cycle

   !              nsdb = nsd + nsb

   !              !do ig = 1,nspatial
   !              !   call calc_Jr(Jr_dir,ig,nra, nla, nrb, nlb, nrc, nlc, nrd, nld)
   !              !   call calc_Jr(Jr_exc,ig,nra, nla, nrb, nlb, nrd, nld, nrc, nlc)
   !              !   call calc_Jr(Jr_dir,ig,nra, nla, nrb,-nlb, nrc, nlc, nrd,-nld)
   !              !   call calc_Jr(Jr_exc,ig,nra, nla, nrb,-nlb, nrd,-nld, nrc, nlc)
   !              !   call calc_Jr(Jr_del,ig,nra, nla, nrb, nlb, nrc, nlc, nrd, nld)
   !              !   call calc_Jr(Jr_del,ig,nra, nla, nrb, nlb, nrc, nlc, nrd, nld)
   !              !   MEJr(ig,ipt) = Jr_dir
   !              !enddo

   !              ! EXCHANGE: ABDC
   !              ir=rindex(nra,nrb,nrd,nrc,nla,nlb,nld,nlc,n)
   !              if(index_flag(ir).eq.0) then
   !                 nr_flag(1,ilauf) = nra
   !                 nr_flag(2,ilauf) = nrb
   !                 nr_flag(3,ilauf) = nrd
   !                 nr_flag(4,ilauf) = nrc
   !                 nl_flag(1,ilauf) = nla
   !                 nl_flag(2,ilauf) = nlb
   !                 nl_flag(3,ilauf) = nld
   !                 nl_flag(4,ilauf) = nlc
   !                 ir_flag(ilauf) = ir
   !                 index_flag(ir) = 1
   !                 ilauf = ilauf + 1
   !              endif

   !              ! DEL: EXCHANGE=ABDC  =  ADBC(B=D)  =  ACBD(exc)  =  A-CB-D(TR)=DEL
   !              !if(nrb.ne.nrd.or.nlb.ne.nld) then
   !                 ir=rindex(nra,nrc,nrb,nrd,nla,-nlc,nlb,-nld,n)
   !                 if(index_flag(ir).eq.0) then
   !                    nr_flag(1,ilauf) = nra
   !                    nr_flag(2,ilauf) = nrc
   !                    nr_flag(3,ilauf) = nrb
   !                    nr_flag(4,ilauf) = nrd
   !                    nl_flag(1,ilauf) = nla
   !                    nl_flag(2,ilauf) =-nlc
   !                    nl_flag(3,ilauf) = nlb
   !                    nl_flag(4,ilauf) =-nld
   !                    ir_flag(ilauf) = ir
   !                    index_flag(ir) = 1
   !                    ilauf = ilauf + 1
   !                 endif
   !              !endif

   !              ! DIRECT: EXCHANGE=ABDC  =  ABCD(C=D)=DIRECT
   !              !     or  DEL=A-CB-D  = ACBD(TR)  =  ABCD(C=B)=DIRECT
   !              !if((nrc.ne.nrd.or.nlc.ne.nld)) then!.and.&
   !              !     !(nrc.ne.nrb.or.nlc.ne.nlb)) then
   !                 ir=rindex(nra,nrb,nrc,nrd,nla,nlb,nlc,nld,n)
   !                 if(index_flag(ir).eq.0) then
   !                    nr_flag(1,ilauf) = nra
   !                    nr_flag(2,ilauf) = nrb
   !                    nr_flag(3,ilauf) = nrc
   !                    nr_flag(4,ilauf) = nrd
   !                    nl_flag(1,ilauf) = nla
   !                    nl_flag(2,ilauf) = nlb
   !                    nl_flag(3,ilauf) = nlc
   !                    nl_flag(4,ilauf) = nld
   !                    ir_flag(ilauf) = ir
   !                    index_flag(ir) = 1
   !                    ilauf = ilauf + 1
   !                 endif
   !              !endif

   !              ! TIME REVERSED: LB,LD --> -LB,-LD, no change unless /= 0
   !              if(nlb.ne.0.or.nld.ne.0) then
   !                 ! EXCHANGE
   !                 ir=rindex(nra,nrb,nrd,nrc,nla,-nlb,-nld,nlc,n)
   !                 if(index_flag(ir).eq.0)then
   !                    nr_flag(1,ilauf) = nra
   !                    nr_flag(2,ilauf) = nrb
   !                    nr_flag(3,ilauf) = nrd
   !                    nr_flag(4,ilauf) = nrc
   !                    nl_flag(1,ilauf) = nla
   !                    nl_flag(2,ilauf) =-nlb
   !                    nl_flag(3,ilauf) =-nld
   !                    nl_flag(4,ilauf) = nlc
   !                    ir_flag(ilauf) = ir
   !                    index_flag(ir) = 1
   !                    ilauf = ilauf+1
   !                 endif
   !                 ! DEL: EXCHANGE=ABDC  =  ADBC(B=D)  =  ACBD(exc)  =  A-CB-D(TR)=DEL
   !                 !if(nrb.ne.nrd .or.nlb.ne.nld) then
   !                    ir=rindex(nra,nrc,nrb,nrd,nla,-nlc,-nlb,nld,n)
   !                    if(index_flag(ir).eq.0) then
   !                       nr_flag(1,ilauf) = nra
   !                       nr_flag(2,ilauf) = nrc
   !                       nr_flag(3,ilauf) = nrb
   !                       nr_flag(4,ilauf) = nrd
   !                       nl_flag(1,ilauf) = nla
   !                       nl_flag(2,ilauf) =-nlc
   !                       nl_flag(3,ilauf) =-nlb
   !                       nl_flag(4,ilauf) = nld
   !                       ir_flag(ilauf) = ir
   !                       index_flag(ir) = 1
   !                       ilauf = ilauf + 1
   !                    endif
   !                 !endif
   !                 ! DIRECT: EXCHANGE=ABDC  =  ABCD(C=D)=DIRECT
   !                 !     or  DEL=A-CB-D  = ACBD(TR)  =  ABCD(C=B)=DIRECT
   !                 !if((nrc.ne.nrb.or.nlb.ne.nlc)) then!.and.&
   !                 !     !(nrc.ne.nrd.or.nlc.ne.nld)) then
   !                    ir=rindex(nra,nrb,nrc,nrd,nla,-nlb,nlc,-nld,n)
   !                    if(index_flag(ir).eq.0) then
   !                       nr_flag(1,ilauf)= nra
   !                       nr_flag(2,ilauf)= nrb
   !                       nr_flag(3,ilauf)= nrc
   !                       nr_flag(4,ilauf)= nrd
   !                       nl_flag(1,ilauf)= nla
   !                       nl_flag(2,ilauf)=-nlb
   !                       nl_flag(3,ilauf)= nlc
   !                       nl_flag(4,ilauf)=-nld
   !                       ir_flag(ilauf)=ir
   !                       index_flag(ir)=1
   !                       ilauf = ilauf + 1
   !                    endif
   !                 !endif
   !              endif ! ld,lb/=0?
   !           enddo !itd
   !        enddo !itb
   !     enddo !itc
   !  enddo !ita

   !  ! The actual calculation of the matrix elements is done below using
   !  ! the flags that were set on the previous loop
   !  do ii = 0,ilauf-1
   !     nra = nr_flag(1,ii)
   !     nrb = nr_flag(2,ii)
   !     nrc = nr_flag(3,ii)
   !     nrd = nr_flag(4,ii)
   !     nla = nl_flag(1,ii)
   !     nlb = nl_flag(2,ii)
   !     nlc = nl_flag(3,ii)
   !     nld = nl_flag(4,ii)
   !     ir_abcd = ir_flag(ii)
   !     do ig = 1,nspatial
   !        call calc_Jr(Jr,ig, nra,nla,nrb,nlb,nrc,nlc,nrd,nld)
   !        MEJr(ig,ir_abcd) = Jr
   !        if (extfield_2bc_is_zero(Jr).and. ig==2) then
   !            print *, "R", ig,nra,nrb,nrc,nrd,nla,nlb,nlc,nld
   !        end if
   !     enddo
   !  enddo
   !end subroutine CalculateMEJr




end module type_extfield_2bc
