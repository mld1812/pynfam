!------------------------------------------------------------------------------
!> This module contains routines for computing matrix elements of
!> charge-changing external fields in single-particle configuration space and
!> storing the properties of the external field in pnfam's type_extfield.
!>
!> @authors M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
!> @authors E.M. Ney, UNC Chapel Hill, 2018-
!------------------------------------------------------------------------------
module pnfam_extfield
    use pnfam_logger
    use pnfam_constants, only : pi
    implicit none
    private
 
    integer,  parameter :: dp = kind(1d0)
 
    character(len=1), dimension(0:1), parameter :: cbeta = (/'+', '-'/)
 
    public :: init_external_field
    public :: number_crossterms
    public :: setup_crossterms
    public :: init_fam_mapping
    public :: deallocate_extfield
    public :: deallocate_crossterms
    public :: set_use_2bc
    public :: tbc_nmlda_da1
    public :: dme_exc
 
 contains
 
    !----------------------------------------------------------------------------
    ! Take in a label and K-projection and create an external_field type.
    ! This is a more consolidated routine than previously: now we assign the
    ! various quantities and ALSO set up the FAM s.p. matrix itself. It is
    ! accessible via op%mat%elem (cf. f%elem previously).
    !----------------------------------------------------------------------------
    subroutine init_external_field(beta_type, label, k, op, use_2bc)
       use type_blockmatrix, only : allocate_blockmatrix, copy_block_structure
       use pnfam_constants, only : translate_uppercase
       use type_extfield
       implicit none
 
       character(len=1),     intent(in)    :: beta_type
       character(len=*),     intent(in)    :: label
       integer,              intent(in)    :: k
       type(external_field), intent(inout) :: op
       integer, optional, intent(in) :: use_2bc
 
       ! Operator label
       op%label = trim(label)
       call translate_uppercase(op%label)
 
       ! Beta minus
       if (trim(beta_type) == '-') then
          op%beta_minus = .true.
       else if (trim(beta_type) == '+') then
          op%beta_minus = .false.
       else
          call error_unknown_operator('beta_type='//trim(beta_type))
       end if
 
       ! Parity
       select case (op%label)
          ! Allowed decays
          case ('F', 'GT')
             op%parity_even = .true.
          ! Forbidden decays
          case ('RS0', 'RS1', 'RS2', 'R', 'P', 'PS0')
             op%parity_even = .false.
          case default
             call error_unknown_operator(trim(op%label))
       end select
 
       ! Angular momentum projection
       op%k = k
 
       ! Rank (needed for rotational energy correction)
       select case (op%label)
          case ('F', 'RS0', 'PS0')
             op%rank = 0
          case ('GT', 'R', 'P', 'RS1')
             op%rank = 1
          case ('RS2')
             op%rank = 2
          case default
             call writelog(" ERROR: Rank not implemented for the operator")
             call writelog(" This should not happen, fix it in extfield.f90")
             call abort
       end select
 
       ! Two body currents
       op%use_2bc = 0
       if (present(use_2bc)) then
          call set_use_2bc(op, use_2bc)
       end if
 
       ! Matrix structure
       !   - 2 body currents: Add additional blockmatrix for pairing field
       !   - 1 body GT is (f11, 0; 0, 0), GT2BC has (f11, f12; -f12*, -f11*)
       call init_fam_mapping(op)
       if (op%use_2bc(3) > 1 .and. op%use_2bc(4) /= 0) then ! Delta needed if 2 or 3
          call allocate_blockmatrix(op%mat12, size(op%mat%elem))
          call copy_block_structure(op%mat, op%mat12)
       end if
 
       call ext_field_operator(op)
 
    end subroutine init_external_field
 
    subroutine ext_field_operator(op)
       use iso_fortran_env, only : error_unit
       use type_extfield
       use type_blockmatrix, only : nb, db, isstart
       use hfb_solution, only : nl, ns, nghl, wf, wfdr, wfdz, y, z
       use hfb_solution, only : hfb_density_coord, nr, nz, npar
       use type_extfield_2bc, only : caux, cd, init_extfield_2bc_type
       use pnfam_constants, only : IT_NEUTRON, IT_PROTON
       ! DME Direct
       !use type_extfield_2bc, only : c3
       !use pnfam_constants, only : Mpi, hbarc
       !use hfb_solution, only : d2rho, drrho, dzrho, get_wf_2nd_derivs
 
       implicit none
       type(external_field), intent(inout) :: op
 
       ! These are assigned from the operator
       integer :: K
       logical :: pty, lpr
       character(len=80) :: label
 
       integer :: ipt, i1, i2, ix1, ix2, ibx1, ibx2, nd1, nd2
       integer :: xl1, xl2, xs1, xs2
       real(dp), dimension(nghl) :: wf_1, wf_2, dr_wf_2, dz_wf_2, dr_wf_1, dz_wf_1, r
       real(dp), dimension(nghl) :: rho_fac, rhon, rhop
       real(dp), dimension(nghl) :: correction_2bc_vector, correction_2bc_axial_charge, correction_2bc_rsL
       ! DME Direct
       !real(dp), dimension(nghl) :: d2z_wf_1, d2r_wf_1, drz_wf_1
       !real(dp), dimension(nghl) :: d2z_wf_2, d2r_wf_2, drz_wf_2
       !real(dp), dimension(nghl) :: r0
       !real(dp) :: dme_c1, me, pre
 
       integer :: debug=0
       ! Define r(:) = 1/(1/rho)
       r(:) = 1.0_dp/y(:)
 
       ! Quantities taken from the operator
       K = op%k
       pty = op%parity_even
       label = op%label
 
       ! Basic tests
       call assert_parity(label=label, pty=pty)
       call assert_K_in_range(label=label, K=K)
 
       ! Compute contact term (sigma tau rho) and/or NM+LDA (sigma tau F(rho))
       rho_fac = 1.0_dp
       if (op%use_2bc(2) /= 0 .and. op%use_2bc(4) /= 0 .and. label == 'GT') then !11/5/23: only work this out if we're looking at GT operator.
          call init_extfield_2bc_type(.false.)
          rho_fac = 0; rhon = 0; rhop = 0
          call hfb_density_coord(IT_NEUTRON,rhon)
          call hfb_density_coord(IT_PROTON, rhop)
 
          ! DME direct
          !dme_c1 = -(caux*4.0_dp*c3)!*((hbarc*hbarc)/(Mpi*Mpi))
          !r0 = 0
          !r0 = rhon + rhop ! N2LO
          !r0 = r0 + (d2rho(:,IT_NEUTRON)+d2rho(:,IT_PROTON))*(hbarc*hbarc)/(Mpi*Mpi) ! N4LO
 
          ! Contact term (same for SNM, ASNM, DME, Full-TBC)
          rho_fac = (caux*2.0_dp*cd)*(rhon + rhop)
          if (op%use_2bc(2) == 2) then
             ! Symmetric nuclear matter calculation at P=0
             rho_fac = rho_fac + tbc_nmlda_da1(rhon, rhop, 0.0_dp, .true.)
          else if (op%use_2bc(2) == 3) then
             ! Asymmetric nuclear matter calculation at P=0
             rho_fac = rho_fac + tbc_nmlda_da1(rhon, rhop, 0.0_dp, .false.)
          else if (op%use_2bc(2) == 4 .or. op%use_2bc(2) == 5) then
             ! DME exchange term
             rho_fac = rho_fac + dme_exc()
          end if
       end if
       ! 3/7/24 Need to call init_extfield_2bc_type to initialize those variables when label = rsL and two body current mode is active
       if (op%use_2bc(4) /= 0 .and. (label == 'RS0' .or. label == 'RS1' .or. label == 'RS2')) then
          call init_extfield_2bc_type(.false.)
       end if
       ! Loop over the FAM block structure to do the calculation
       ipt = 0; op%mat%elem(:) = 0
       do ibx1 = 1, nb
 
          ! Determine the range of particle indices in this block
          ibx2 = op%mat%ir2c(ibx1)
          if (ibx2 == 0) cycle
 
          nd1 = db(ibx1) ; nd2 = db(ibx2) !number of basis states in the blocks ibx1, ibx2.
          if (nd1 == 0 .or. nd2 == 0) cycle
 
          do i2=1, nd2 !loop over basis states in the block ibx2. 
             ! Calculate the true index of the state |2>
             ix2 = i2 + isstart(ibx2) - 1 !correctly index the basis state by adding the index of the first state in block ibx2.
 
             ! Wave functions
             wf_2 = wf(:,ix2) !get the wavefunction data from wf, along with derivatives.
             dr_wf_2 = wfdr(:,ix2)
             dz_wf_2 = wfdz(:,ix2)
 
             ! Quantum numbers
             xl2 = nl(ix2);  xs2 = ns(ix2) !get l and s quantum numbers of the wavefunction.
 
             do i1=1, nd1 !loop over basis states in block ibx1.
                ! Calculate the true index of the state <1|
                ix1 = i1 + isstart(ibx1) - 1 !correctly index as above.
 
                ! Wave functions
                wf_1 = wf(:,ix1)
                dr_wf_1 = wfdr(:,ix1)
                dz_wf_1 = wfdz(:,ix1)
 
                ! Quantum numbers
                xl1 = nl(ix1);  xs1 = ns(ix1)
 
                ! Matrix element identifier
                ipt = ipt + 1
 
                !----------------------------------------------------------------
                ! Definition of the various external field s.p. operators
                !----------------------------------------------------------------
                select case (label)
 
                   ! Unknown operator
                   ! ------------------------------
                   case default
                      call error_unknown_operator(trim(label))
 
 
                   ! Fermi (1)
                   ! ------------------------------
                   case ('F')
                      if ((xl1 == xl2) .and. (xs1 == xs2)) then
                         op%mat%elem(ipt) = dot_product(wf_1(:), wf_2(:)) !just take dot product of the wavefunctions.
                      end if
 
 
                   ! Gamow-Teller (sigma_K)
                   ! ------------------------------
                   case ('GT')
                      lpr = .false.
                      select case (K)
                         case (0)
                            if ((xl1 == xl2) .and. (xs1 == xs2)) then
                               op%mat%elem(ipt) = xs1*dot_product(wf_1(:), rho_fac(:)*wf_2(:))
                               lpr = .true.
                            end if
                         case (1,-1)
                            if ((xl1 == xl2) .and. (xs1 == xs2 + 2*K)) then
                               op%mat%elem(ipt) = -K*sqrt(2.0_dp)*dot_product(wf_1(:), rho_fac(:)*wf_2(:))
                               lpr = .true.
                            end if
                      end select
 
                      !! Gamow-Teller 2BC DME ((sigma.Del)Del(rho + Del^2 rho/m^2)) - IBP twice
                      !! ------------------------------
                      !if (op%use_2bc(2) == 4) then
                      !   me = 0
                      !   call get_wf_2nd_derivs(ix2, d2z_wf_2, d2r_wf_2, drz_wf_2)
                      !   call get_wf_2nd_derivs(ix1, d2z_wf_1, d2r_wf_1, drz_wf_1)
                      !   select case (K)
                      !      case (0)
                      !         pre = dme_c1
                      !         if ((xl1 == xl2) .and. (xs1 == xs2)) then ! Sz
                      !            me = xs1*pre*sum((d2z_wf_1*wf_2 + wf_1*d2z_wf_2 + 2.0_dp*dz_wf_1*dz_wf_2)*r0)
                      !            op%mat%elem(ipt) = op%mat%elem(ipt) + me
                      !         else if (((xl1 == xl2 - 1) .and. (xs1 == xs2 + 2)) .or. &
                      !                  ((xl1 == xl2 + 1) .and. (xs1 == xs2 - 2))) then ! S+ or S-
                      !            me = pre*sum((drz_wf_1*wf_2 + wf_1*drz_wf_2 + dr_wf_1*dz_wf_2 + dz_wf_1*dr_wf_2 &
                      !                          + y*(dz_wf_1*wf_2 + wf_1*dz_wf_2))*r0)
                      !            op%mat%elem(ipt) = op%mat%elem(ipt) + me
                      !         end if
 
                      !      case (1,-1)
                      !         pre = -K*0.5_dp*sqrt(2.0_dp)*dme_c1
                      !         if ((xl1 == xl2 + K) .and. (xs1 == xs2)) then ! Sz
                      !            me = xs1*pre*sum((drz_wf_1*wf_2 + wf_1*drz_wf_2 + dr_wf_1*dz_wf_2 + dz_wf_1*dr_wf_2 &
                      !                              + y*(dz_wf_1*wf_2 + wf_1*dz_wf_2))*r0)
                      !            op%mat%elem(ipt) = op%mat%elem(ipt) + me
                      !         else if ((xl1 == xl2 + K - 1) .and. (xs1 == xs2 + 2)) then ! S+
                      !            me = pre*sum((d2r_wf_1*wf_2 + wf_1*d2r_wf_2 + 2.0_dp*dr_wf_1*dr_wf_2 &
                      !                          + (2 - K)*y*(dr_wf_1*wf_2 + wf_1*dr_wf_2))*r0)
                      !            op%mat%elem(ipt) = op%mat%elem(ipt) + me
                      !         else if ((xl1 == xl2 + K + 1) .and. (xs1 == xs2 - 2)) then ! S-
                      !            me = pre*sum((d2r_wf_1*wf_2 + wf_1*d2r_wf_2 + 2.0_dp*dr_wf_1*dr_wf_2 &
                      !                          + (2 + K)*y*(dr_wf_1*wf_2 + wf_1*dr_wf_2))*r0)
                      !            op%mat%elem(ipt) = op%mat%elem(ipt) + me
                      !         end if
 
                      !   end select
                      !end if
 
                      if (lpr.and.debug<0) then
                      write(*,'(L2,2x,3i4,2x,1e30.20, " |O: ", 2i3, " |L: ", 2i3, " |s: ", 2i3, " |r: ", 2i3, " |z: ", 2i3,&
                          &" |p: ", 2i3, " |diff_O,L,S,P: ",4i3)') (op%mat%ir2c(ibx1) == ibx2), &
                      & ix1, ix2, ipt, op%mat%elem(ipt), 2*xl1+xs1, 2*xl2+xs2, xl1, xl2, xs1, xs2, nr(ix1), nr(ix2),nz(ix1), nz(ix2), &
                      & npar(ix1), npar(ix2),  2*xl1+xs1-(2*xl2+xs2), xl1-xl2, xs1-xs2, npar(ix1)-npar(ix2)
                      end if
 
 
                   ! R (Sqrt(4\pi/3)*rY_1K = r_K)
                   ! ------------------------------
                   case ('R')
                      select case (K)
                         case (0)
                            if ((xl1 == xl2) .and. (xs1 == xs2)) then
                               op%mat%elem(ipt) = dot_product(wf_1(:), z(:)*wf_2(:)) !element multiply wavefunction 2
                            end if
                         case (1,-1)
                            if ((xl1 == xl2 + K) .and. (xs1 == xs2)) then
                               op%mat%elem(ipt) = -K/sqrt(2.0_dp)*dot_product(wf_1(:), r(:)*wf_2(:)) !element-multiply wavefunction 2 by the radius. 
                            end if
                      end select
 
 
                   ! P (p_K/i)
                   ! ------------------------------
                   case ('P')
 
                      select case (K)
                         case (0)
                            if ((xl1 == xl2) .and. (xs1 == xs2)) then
                               if (op%use_2bc(5) == 1) then
                                  correction_2bc_vector = forbidden_2bc_P()
                               else if (op%use_2bc(5) == 2) then !use DME
                                  correction_2bc_vector = dme_vector()
                               end if
                               if (op%use_2bc(1) == 0 .or. op%use_2bc(5) == 0) then !just 1 body.
                                  op%mat%elem(ipt) = -dot_product(wf_1(:), dz_wf_2(:)) !derivative of wavefunction 2 in z direction.
                               else if (op%use_2bc(1) == 1) then !1 body and 2 body.
                                  op%mat%elem(ipt) = -dot_product(wf_1(:) * (1 + correction_2bc_vector(:)), dz_wf_2(:)) 
                               else if (op%use_2bc(1) == 2) then !just 2 body.
                                  op%mat%elem(ipt) = -dot_product(wf_1(:) * correction_2bc_vector(:), dz_wf_2(:)) 
                               end if
                            end if
                         case (1,-1)
                            if ((xl1 == xl2 + K) .and. (xs1 == xs2)) then
                               if (op%use_2bc(5) == 1) then
                                  correction_2bc_vector = forbidden_2bc_P()
                               else if (op%use_2bc(5) == 2) then !use DME
                                  correction_2bc_vector = dme_vector()
                               end if
                               if (op%use_2bc(1) == 0 .or. op%use_2bc(5) == 0) then !just 1 body.
                                  op%mat%elem(ipt) = K/sqrt(2.0_dp) * (dot_product(wf_1(:), dr_wf_2(:))   &
                                  - K*xl2 * dot_product(wf_1(:), y(:)*wf_2(:))) !radial and angular components of gradient. y corresponds to 1/r. 
                               else if (op%use_2bc(1) == 1) then !1 body and 2 body.
                                  op%mat%elem(ipt) = K/sqrt(2.0_dp)* (dot_product(wf_1(:) * (1.0_dp + correction_2bc_vector(:)) , dr_wf_2(:))   &
                                  - K*xl2*dot_product(wf_1(:) * (1.0_dp + correction_2bc_vector(:)), y(:)*wf_2(:)))
                               else if (op%use_2bc(1) == 2) then !just 2 body.
                                  op%mat%elem(ipt) = K/sqrt(2.0_dp)*(dot_product(wf_1(:) * correction_2bc_vector(:), dr_wf_2(:))   &
                                  - K*xl2*dot_product(wf_1(:) * correction_2bc_vector(:), y(:)*wf_2(:))) 
                               end if
                            end if
                         end select
 
 
                   ! RS0 ([rY_1 x sigma]_00/Y_00)
                   ! ------------------------------
                   case ('RS0')
                      if ((xl1 == xl2) .and. (xs1 == xs2)) then
                         if (op%use_2bc(4) == 1) then 
                            correction_2bc_rsL = forbidden_2bc_rsL()
                            op%mat%elem(ipt) = -xs1*dot_product(wf_1(:) * (1.0_dp - correction_2bc_rsL), z(:)*wf_2(:))
                         else
                            op%mat%elem(ipt) = -xs1*dot_product(wf_1(:), z(:)*wf_2(:))
                         end if
                      else if ((xl1 == xl2 + 1) .and. (xs1 == xs2 - 2)) then
                         if (op%use_2bc(4) == 1) then 
                            correction_2bc_rsL = forbidden_2bc_rsL()
                            op%mat%elem(ipt) = -dot_product(wf_1(:) * (1.0_dp - correction_2bc_rsL), r(:)*wf_2(:))
                         else
                            op%mat%elem(ipt) = -dot_product(wf_1(:), r(:)*wf_2(:))
                         end if
                      else if ((xl1 == xl2 - 1) .and. (xs1 == xs2 + 2)) then
                         if (op%use_2bc(4) == 1) then 
                            correction_2bc_rsL = forbidden_2bc_rsL()
                            op%mat%elem(ipt) = -dot_product(wf_1(:) * (1.0_dp - correction_2bc_rsL), r(:)*wf_2(:))
                         else
                            op%mat%elem(ipt) = -dot_product(wf_1(:), r(:)*wf_2(:))
                         end if
                      end if
 
 
                   ! RS1 ([rY_1 x sigma]_1K/Y_00)
                   ! ------------------------------
                   case ('rs1','RS1')
                      select case (K)
                         case (0)
                            if ((xl1 == xl2 - 1) .and. (xs1 == xs2 + 2)) then
                               if (op%use_2bc(4) == 1) then
                                  correction_2bc_rsL = forbidden_2bc_rsL()
                                  op%mat%elem(ipt) = sqrt(3.0_dp/2.0_dp)*dot_product(wf_1(:) * (1.0_dp - correction_2bc_rsL),r(:)*wf_2(:))
                               else
                                  op%mat%elem(ipt) = sqrt(3.0_dp/2.0_dp)*dot_product(wf_1(:),r(:)*wf_2(:))
                               end if
                            else if ((xl1 == xl2 + 1) .and. (xs1 == xs2 - 2)) then
                               if (op%use_2bc(4) == 1) then
                                  correction_2bc_rsL = forbidden_2bc_rsL()
                                  op%mat%elem(ipt) = -sqrt(3.0_dp/2.0_dp)*dot_product(wf_1(:)* (1.0_dp - correction_2bc_rsL),r(:)*wf_2(:))
                               else
                                  op%mat%elem(ipt) = -sqrt(3.0_dp/2.0_dp)*dot_product(wf_1(:),r(:)*wf_2(:))
                               end if 
                            end if
                         case (1,-1)
                            if ((xl1 == xl2) .and. (xs1 == xs2 + 2*K)) then
                               if (op%use_2bc(4) == 1) then
                                  correction_2bc_rsL = forbidden_2bc_rsL()
                                  op%mat%elem(ipt) = sqrt(3.0_dp)*dot_product(wf_1(:)* (1.0_dp - correction_2bc_rsL),z(:)*wf_2(:))
                               else
                                  op%mat%elem(ipt) = sqrt(3.0_dp)*dot_product(wf_1(:), z(:)*wf_2(:))
                               end if
                            else if ((xl1 == xl2 + K) .and. (xs1 == xs2)) then
                               if (op%use_2bc(4) == 1) then
                                  correction_2bc_rsL = forbidden_2bc_rsL()
                                  op%mat%elem(ipt) = -xs1*sqrt(3.0_dp)/2.0_dp*dot_product(wf_1(:)* (1.0_dp - correction_2bc_rsL),r(:)*wf_2(:))
                               else
                                  op%mat%elem(ipt) = -xs1*sqrt(3.0_dp)/2.0_dp*dot_product(wf_1(:), r(:)*wf_2(:))
                               end if
                            end if
                      end select
 
 
                   ! RS2 ([rY_1 x sigma]_2K/Y_00)
                   ! ------------------------------
                   case ('RS2')
                      select case (K)
                         case (0)
                            if ((xl1 == xl2) .and. (xs1 == xs2)) then
                               if (op%use_2bc(4) == 1) then
                                  correction_2bc_rsL = forbidden_2bc_rsL()
                                  op%mat%elem(ipt) = xs1*sqrt(2.0_dp)*dot_product(wf_1(:) * (1.0_dp - correction_2bc_rsL), z(:)*wf_2(:))
                               else
                                  op%mat%elem(ipt) = xs1*sqrt(2.0_dp)*dot_product(wf_1(:), z(:)*wf_2(:))
                               end if 
                            else if ((xl1 == xl2 + 1) .and. (xs1 == xs2 - 2)) then
                               if (op%use_2bc(4) == 1) then
                                  correction_2bc_rsL = forbidden_2bc_rsL()
                                  op%mat%elem(ipt) = -1/sqrt(2.0_dp)*dot_product(wf_1(:) * (1.0_dp - correction_2bc_rsL), r(:)*wf_2(:))
                               else
                                  op%mat%elem(ipt) = -1/sqrt(2.0_dp)*dot_product(wf_1(:), r(:)*wf_2(:))
                               end if 
                            else if ((xl1 == xl2 - 1) .and. (xs1 == xs2 + 2)) then
                               if (op%use_2bc(4) == 1) then
                                  correction_2bc_rsL = forbidden_2bc_rsL()
                                  op%mat%elem(ipt) = -1/sqrt(2.0_dp)*dot_product(wf_1(:) * (1.0_dp - correction_2bc_rsL), r(:)*wf_2(:))
                               else
                                  op%mat%elem(ipt) = -1/sqrt(2.0_dp)*dot_product(wf_1(:), r(:)*wf_2(:))
                               end if
                            end if
                         case (1,-1)
                            if ((xl1 == xl2) .and. (xs1 == xs2 + 2*K)) then
                               if (op%use_2bc(4) == 1) then
                                  correction_2bc_rsL = forbidden_2bc_rsL()
                                  op%mat%elem(ipt) = -K*sqrt(3.0_dp)*dot_product(wf_1(:) * (1.0_dp - correction_2bc_rsL), z(:)*wf_2(:))
                               else
                                  op%mat%elem(ipt) = -K*sqrt(3.0_dp)*dot_product(wf_1(:), z(:)*wf_2(:))
                               end if
                            else if ((xl1 == xl2 + K) .and. (xs1 == xs2)) then
                               if (op%use_2bc(4) == 1) then
                                  correction_2bc_rsL = forbidden_2bc_rsL()
                                  op%mat%elem(ipt) = -K*xs1*sqrt(3.0_dp)/2.0_dp*dot_product(wf_1(:) * (1.0_dp - correction_2bc_rsL), r(:)*wf_2(:))
                               else
                                  op%mat%elem(ipt) = -K*xs1*sqrt(3.0_dp)/2.0_dp*dot_product(wf_1(:), r(:)*wf_2(:))
                               end if
                            end if
                         case (2,-2)
                            if ((xl1 == xl2 + K/2) .and. (xs1 == xs2 + K)) then
                               if (op%use_2bc(4) == 1) then
                                  correction_2bc_rsL = forbidden_2bc_rsL()
                                  op%mat%elem(ipt) = sqrt(3.0_dp)*dot_product(wf_1(:) * (1.0_dp - correction_2bc_rsL), r(:)*wf_2(:))
                               else
                                  op%mat%elem(ipt) = sqrt(3.0_dp)*dot_product(wf_1(:), r(:)*wf_2(:))
                               end if 
                            end if
                      end select
 
 
                      ! PS0 (p.sigma/i)
                      ! ----------------------------------------
                      case ('PS0')
                         if ((xl1 == xl2) .and. (xs1 == xs2)) then
                            if (op%use_2bc(6) == 1) then !assign 2bc correction array
                               correction_2bc_axial_charge = forbidden_2bc_PS0()
                            else if (op%use_2bc(6) == 2) then !use DME
                               correction_2bc_axial_charge = dme_axial()
                            end if
                            if (op%use_2bc(1) == 0 .or. op%use_2bc(6) == 0) then !just one body
                               op%mat%elem(ipt) = -xs1*dot_product(wf_1(:), dz_wf_2(:))
                            else if (op%use_2bc(1) == 1) then !1bc + 2bc
                               op%mat%elem(ipt) = -xs1*dot_product(wf_1(:) * (1 + correction_2bc_axial_charge(:)), dz_wf_2(:))
                            else if (op%use_2bc(1) == 2) then !just 2bc
                               op%mat%elem(ipt) = -xs1*dot_product(wf_1(:) * correction_2bc_axial_charge(:), dz_wf_2(:))
                            end if
                         else if ((xl1 == xl2 - 1) .and. (xs1 == xs2 + 2)) then
                            if (op%use_2bc(6) == 1) then
                               correction_2bc_axial_charge = forbidden_2bc_PS0()
                            else if (op%use_2bc(6) == 2) then !use DME
                               correction_2bc_axial_charge = dme_axial()
                            end if
                            if (op%use_2bc(1) == 0 .or. op%use_2bc(6) == 0) then !just one body
                               op%mat%elem(ipt) = -dot_product(wf_1(:), dr_wf_2(:))      &
                                   - xl2*dot_product(wf_1(:), y(:)*wf_2(:))
                            else if (op%use_2bc(1) == 1) then !1bc + 2bc
                               op%mat%elem(ipt) = -dot_product(wf_1(:) * (1 + correction_2bc_axial_charge(:)), dr_wf_2(:))      &
                                   - xl2*dot_product(wf_1(:) * (1 + correction_2bc_axial_charge(:)), y(:)*wf_2(:))
                            else if (op%use_2bc(1) == 2) then !just 2bc
                               op%mat%elem(ipt) = -dot_product(wf_1(:) * correction_2bc_axial_charge, dr_wf_2(:))      &
                               - xl2*dot_product(wf_1(:) * correction_2bc_axial_charge, y(:)*wf_2(:))
                            end if
                         else if ((xl1 == xl2 + 1) .and. (xs1 == xs2 - 2)) then
                            if (op%use_2bc(6) == 1) then
                               correction_2bc_axial_charge = forbidden_2bc_PS0()
                            else if (op%use_2bc(6) == 2) then !use DME
                               correction_2bc_axial_charge = dme_axial()
                            end if
                            if (op%use_2bc(1) == 0 .or. op%use_2bc(6) == 0) then !just one body
                               op%mat%elem(ipt) = -dot_product(wf_1(:), dr_wf_2(:))      &
                                   + xl2*dot_product(wf_1(:), y(:)*wf_2(:))
                            else if (op%use_2bc(1) == 1) then !1bc + 2bc
                               op%mat%elem(ipt) = -dot_product(wf_1(:) * (1 + correction_2bc_axial_charge(:)), dr_wf_2(:))      &
                                   + xl2*dot_product(wf_1(:) * (1 + correction_2bc_axial_charge(:)), y(:)*wf_2(:))
                            else if (op%use_2bc(1) == 2) then !just 2bc
                               op%mat%elem(ipt) = -dot_product(wf_1(:) * correction_2bc_axial_charge(:), dr_wf_2(:))      &
                                   + xl2*dot_product(wf_1(:) * correction_2bc_axial_charge(:), y(:)*wf_2(:))
                            end if
                         end if
 
                end select ! operator_name
 
             end do ! <1|
          end do ! |2>
       end do ! blocks
 
       if (debug/=0) then
          print *, "END GT 1Body"
       end if
 
    end subroutine ext_field_operator
 
 
    !----------------------------------------------------------------------------
    ! Assertion regarding parity. There has to be a good way to store all these
    ! properties (parity, |K_max|, etc.) of all these operators.
    !----------------------------------------------------------------------------
    subroutine assert_parity(label, pty)
       implicit none
       logical,           intent(in) :: pty
       character(len=80), intent(in) :: label
 
       select case (label)
          ! Even parity - parity_even .eqv. .TRUE.
          case ('F', 'GT')
             if (pty .eqv. .false.) then
                call error_wrong_parity(label)
             end if
          ! Odd parity - parity_even .eqv. .FALSE.
          case  ('RS0', 'RS1', 'RS2', 'R', 'P', 'PS0')
             if (pty .eqv. .true.) then
                call error_wrong_parity(label)
             end if
          case default
             call error_unknown_operator(label)
       end select
    end subroutine assert_parity
 
 
    !----------------------------------------------------------------------------
    ! Assertion regarding |K| ranges.
    !----------------------------------------------------------------------------
    subroutine assert_K_in_range(K, label)
       implicit none
       integer,           intent(in) :: K
       character(len=80), intent(in) :: label
 
       select case (label)
          ! Kmax = 0
          case ('F', 'RS0', 'PS0')
             if (abs(K) > 0) then
                call error_K_out_of_range(label, 0)
             end if
          ! Kmax = 1
          case  ('GT', 'R', 'P', 'RS1')
             if (abs(K) > 1) then
                call error_K_out_of_range(label, 1)
             end if
          ! Kmax = 2
          case ('RS2')
             if (abs(K) > 2) then
                call error_K_out_of_range(label, 2)
             end if
          case default
             call error_unknown_operator(label)
       end select
    end subroutine assert_K_in_range
 
    !----------------------------------------------------------------------------
    ! Assertion regarding use_2bc values
    !----------------------------------------------------------------------------
    subroutine set_use_2bc(op, use_2bc)
       use type_extfield
       use pnfam_constants, only : get_digit
       implicit none
       type(external_field), intent(inout) :: op
       integer, intent(in) :: use_2bc
       integer :: ierr, i1, i2, i3, i4, i5, i6
 
       ierr = 0
       op%use_2bc = 0
       if (use_2bc == 0) return
 
       ! Extract digits
       i1 = get_digit(use_2bc, 5) ! 1=1+2Body, 2=2Body only 
       i2 = get_digit(use_2bc, 4) ! 1=Full Fam, 2=SNM+LDA, 3=ASNM+LDA
       i3 = get_digit(use_2bc, 3) ! 1=Gamma, 2=Delta, 3=Gamma+Delta
       i4 = get_digit(use_2bc, 2) ! if >0, then 2bc GT current is active
       i5 = get_digit(use_2bc, 1) ! if >0, then 2bc P current is active
       i6 = get_digit(use_2bc, 0) ! if >0, then 2bc PS0 current is active
 
       op%use_2bc(1) = i1
       op%use_2bc(2) = i2
       op%use_2bc(3) = i3
       op%use_2bc(4) = i4
       op%use_2bc(5) = i5
       op%use_2bc(6) = i6
       ! Check valid digits
       if ((i1 < 1 .or. i1 > 2) .or. (i2 < 1 .or. i2 > 5) .or. (i3 < 1 .or. i3 > 3)) then
          ierr = 1
       ! Check valid combos: NM+LDA (i1=2,3) has only Gamma (i3=1) ... and DME?
       else if (i1 /= 1 .and. i3 /= 1) then
          ierr = 1
       ! Only Gamma is implemented as of June 2021
       else if (i3 /= 1) then
          ierr = 2
       end if
       if (ierr==1) call abort(" Error: Invalid value supplied for two_body_current_mode.")
       if (ierr==2) call abort(" Error: This two_body_current_mode is not yet operational.")
   
    end subroutine set_use_2bc
 
    !----------------------------------------------------------------------------
    ! Error messages for various cases
    !----------------------------------------------------------------------------
    subroutine error_unknown_operator(name)
       implicit none
       character(len=*), intent(in) :: name
       character(len=200) :: st
 
       write(st,'(3a)') 'Error: Unknown operator "', trim(name), '" in extfield.'
       call abort(st)
    end subroutine error_unknown_operator
 
    subroutine error_wrong_parity(name)
       implicit none
       character(len=*), intent(in) :: name
       character(len=200) :: st
 
       write(st,'(3a)') 'Error: Incorrect parity for operator "', trim(name), '" in extfield.'
       call abort(st)
    end subroutine error_wrong_parity
 
    subroutine error_K_out_of_range(name, K_max)
       implicit none
       integer,          intent(in) :: K_max
       character(len=*), intent(in) :: name
       character(len=200) :: st
 
       write(st,'(3a)') 'Error: K out of range for operator "', trim(name), '" in extfield.'
       call writelog(st)
       write(st,'(a, 1i1)')   'Expected |K| <= ', K_max
       call writelog(st)
       call abort
    end subroutine error_K_out_of_range
    
    
    !-----------------------------------------------------------------------------
    ! This subroutine for initializing the f matrix block structure for general K
    ! and parity;  The HFBTHO basis must be read in before calling this!
    !-----------------------------------------------------------------------------
    subroutine init_fam_mapping(op)
       use type_blockmatrix, only : allocate_blockmatrix, nb, db
       use type_extfield
       use hfb_solution, only : nl, ns, npar
       implicit none
       
       type(external_field), intent(inout) :: op
       
       ! These are assigned from the operator
       integer :: k     ! change in angular momentum projection
       logical :: pty   ! .true. = parity conserved, .false. = parity flipped
       
       integer :: i1, i2, ip1, ip2, i
       integer :: aux1(nb), aux2(nb)
       integer :: nbb                          ! the number of FAM blocks
       integer, allocatable :: ib1(:), ib2(:)  ! the THO basis blocks the FAM block connects <1|FAM|2>
       logical :: pty_blocks
       
       ! K and parity come from the operator itself
       k = op%k
       pty = op%parity_even
       
       ! Since there may be empty (missing) blocks in the THO, we have to count
       ! the number of blocks we'll have the brute way
       pty_blocks = .true.
       ip1 = 1
       do i1=1,nb
          if (any(npar(ip1:ip1+db(i1)-1) /= npar(ip1))) then
             pty_blocks = .false.
             exit
          end if
          ip1 = ip1 + db(i1)
       end do
       nbb = 0
       ip1 = 1
       do i1=1,nb
          ip2 = 1
          do i2=1,nb
             ! Use the first state in the block to determine if Omega and parity match
             if ( (2*nl(ip1)+ns(ip1)-2*nl(ip2)-ns(ip2) == 2*k) .and. &
                  (.not. pty_blocks .or. (pty_blocks .and. (npar(ip1)==npar(ip2) .eqv. pty))) ) then
                ! Found a correspondence: add it to the temporary FAM block list
                ! and exit the inner loop
                nbb = nbb + 1
                aux1(nbb) = i1
                aux2(nbb) = i2
                exit
             end if
             ip2 = ip2 + db(i2)
          end do
          ip1 = ip1 + db(i1)
       end do
       
       ! Now that we know how many FAM blocks exist, we can allocate the actual
       ! FAM block arrays and fill them with the temporary arrays
       if (nbb > 0) then
          allocate(ib1(nbb), ib2(nbb))
          ib1(:) = aux1(1:nbb) ; ib2(:) = aux2(1:nbb)
       else
          call abort(" ERROR: The K and parity are outside the model space. Aborting.")
       end if
       
       ! Initialize the block structure of the matrix f
       call allocate_blockmatrix(op%mat, sum(db(ib1(:))*db(ib2(:))))
       do i = 1, nbb
          op%mat%ir2c(ib1(i)) = ib2(i)
          op%mat%ic2r(ib2(i)) = ib1(i)
       end do
       ip1 = 1
       do i = 1, nb
          if (op%mat%ir2c(i) > 0) then
             op%mat%ir2m(i) = ip1
             op%mat%ic2m(op%mat%ir2c(i)) = ip1
             ip1 = ip1 + db(i)*db(op%mat%ir2c(i))
          end if
       end do
       
    end subroutine
    
    
    subroutine deallocate_extfield(e)
       use type_extfield
       use type_blockmatrix, only : deallocate_blockmatrix
       implicit none
       type(external_field), intent(inout) :: e
       
       call deallocate_blockmatrix(e%mat)
       
    end subroutine deallocate_extfield
    
    
    !----------------------------------------------------------------------------
    ! Return the number of possible cross-terms for a given K and parity.
    !----------------------------------------------------------------------------
    function number_crossterms(op)
       use type_extfield
       implicit none
       type(external_field), intent(in) :: op
       integer :: number_crossterms
       
       number_crossterms = 0
       if (op%parity_even .eqv. .false.) then
          select case (trim(op%label))
             case ('R', 'P', 'RS1')
                number_crossterms = 3
             case ('RS0', 'PS0')
                number_crossterms = 2
             case default
                number_crossterms = 0
          end select
       end if
    end function number_crossterms
    
    
    !----------------------------------------------------------------------------
    ! Define the various cross-terms and labels depending on the intrinsic
    ! J-value of the operator
    ! 11/29/23: Add two body current mode handling for the P operator.
    !----------------------------------------------------------------------------
    subroutine setup_crossterms(op, crossterms, n, use_2bc)
       use type_extfield
       implicit none
       type(external_field), intent(in) :: op
       type(external_field), dimension(:), allocatable, intent(inout) :: crossterms
       integer, intent(out) :: n
       integer, optional, intent(in) :: use_2bc
       
       integer :: nxterms, ixterm, ibeta
       
       ! Count the cross-terms and allocate space
       nxterms = number_crossterms(op=op)
       n = nxterms
       
       ! No cross-terms means no work
       if (nxterms == 0) then
          return
       ! Otherwise, allocate the space and set them up
       else
          call deallocate_crossterms(crossterms)
          allocate(crossterms(1:nxterms))
          
          ! Is the operator beta+ or beta-?
          ! 0 if beta+, 1 if beta- to be the argument for the parameter array cbeta
          ibeta = merge(1, 0, op%beta_minus)
          
          ! Ordering is *VERY* important. If any changes are made here, they
          ! should be compared with the expected cross-term orderings in the
          ! computation of forbidden decay rates!
          ! J = 0
          if (nxterms == 2) then
             ! OP x RS0
             ixterm = 1
             call init_external_field(cbeta(ibeta), label='RS0', k=op%k, op=crossterms(ixterm), use_2bc = use_2bc)
             crossterms(ixterm)%label = trim(op%label)//'x'//trim(crossterms(ixterm)%label)
             ! OP x PS0
             ixterm = 2
             call init_external_field(cbeta(ibeta), label='PS0', k=op%k, op=crossterms(ixterm), use_2bc = use_2bc)
             crossterms(ixterm)%label = trim(op%label)//'x'//trim(crossterms(ixterm)%label)
          ! J = 1
          else if (nxterms == 3) then
             ! OP x R
             ixterm = 1
             call init_external_field(cbeta(ibeta), label='R', k=op%k, op=crossterms(ixterm))
             crossterms(ixterm)%label = trim(op%label)//'x'//trim(crossterms(ixterm)%label)
             ! OP x RS1
             ixterm = 2
             call init_external_field(cbeta(ibeta), label='RS1', k=op%k, op=crossterms(ixterm), use_2bc = use_2bc)
             crossterms(ixterm)%label = trim(op%label)//'x'//trim(crossterms(ixterm)%label)
             ! OP x P
             ixterm = 3
             !if (present(use_2bc)) then
             !   call init_external_field(cbeta(ibeta), label='P', k=op%k, op=crossterms(ixterm), use_2bc = use_2bc)
             !else
             call init_external_field(cbeta(ibeta), label='P', k=op%k, op=crossterms(ixterm), use_2bc = use_2bc)
             !end if 
             crossterms(ixterm)%label = trim(op%label)//'x'//trim(crossterms(ixterm)%label)
          end if
       end if
    end subroutine setup_crossterms
    
    
    subroutine deallocate_crossterms(crossterms)
       use type_extfield
       implicit none
       type(external_field), dimension(:), allocatable, intent(inout) :: crossterms
       
       integer :: i
       
       if (allocated(crossterms)) then
          do i=1,size(crossterms)
             call deallocate_extfield(crossterms(i))
          end do
          deallocate(crossterms)
       end if
       
    end subroutine deallocate_crossterms
 
    !----------------------------------------------------------------------------
    ! TBC NUCLEAR MATTER APPROXIMATION + LDA: PRD 88, 083516 (2013), Eq. 16, A11, A12
    ! Note: Eq. 16 has a typo, there is no 1/4m in the c3 term.
    ! Note: Eq. A11/A12 I cannot reproduce the correct result, there might be a typo
    ! Note: Dimensions are rho = 1/fm^3 , kf = 1/fm
    !----------------------------------------------------------------------------
 
    ! P=p=0 limit, I2=I1==I0 (Klos 2013 Eq.14, Verified by Evan 4/23/21)
    function tbc_nmlda_I0(kf) result(I0)
       use pnfam_constants, only : pi, Mpi, hbarc
       use hfb_solution, only : nghl
       implicit none
       real(dp), intent(in) :: kf
       real(dp) :: kf3, kf2, I0
       real(dp) :: m
 
       m = Mpi/hbarc
       kf2 = kf*kf
       kf3 = kf*kf*kf
 
       I0 = 1.0_dp - 3.0_dp*m*m/(kf2) + 3.0_dp*m*m*m/(kf3)*atan(kf/m)
 
    end function
 
    function tbc_nmlda_I1(kf, Q) result(I1)
       use pnfam_constants, only : pi, Mpi, hbarc
       use hfb_solution, only : nghl
       implicit none
       real(dp), intent(in) :: kf, Q
 
       real(dp) :: kf3, kf2
       real(dp) :: aux1, aux2, aux3, aux4, aux5, I1
       real(dp) :: m
 
       m = Mpi/hbarc
       kf2 = kf*kf
       kf3 = kf*kf*kf
 
       aux1 = 1.0_dp/(512.0_dp*kf3*Q**3)
       aux2 = 8*kf*Q*(48*(kf2 + m*m)**2 + 32*(kf2 - 3*m*m)*Q*Q - 3*Q**4)
       aux3 = 768*m**3*Q**3*(atan((2*kf+Q)/(2*m)) + atan((2*kf-Q)/(2*m))) !atan(2*m*kf/(m*m + 0.25_dp*Q*Q - kf2))
       aux4 = 3*(16*(kf2 + m*m)**2 - 8*(kf2 - 5*m*m)*Q*Q + Q**4)*(4*(kf2 + m*m) - Q*Q)
       aux5 = log((m*m + (kf - Q*0.5_dp)**2)/(m*m + (kf + Q*0.5_dp)**2))
       I1 = aux1*(aux2 + aux3 + aux4*aux5)
 
    end function
 
    function tbc_nmlda_I2(kf, Q) result(I2)
       use pnfam_constants, only : pi, Mpi, hbarc
       use hfb_solution, only : nghl
       implicit none
 
       real(dp), intent(in) :: kf, Q 
 
       real(dp) :: kf3, kf2
       real(dp) :: aux1, aux2, aux3, aux4, aux5, I2
       real(dp) :: m
 
       m = Mpi/hbarc
       kf2 = kf*kf
       kf3 = kf*kf*kf
 
       aux1 = 1.0_dp/(16.0_dp*kf3*Q)
       aux2 = 8*kf*(2*kf2 - 3*m*m)*Q
       aux3 = 24*m**3*Q*(atan((2*kf+Q)/(2*m)) + atan((2*kf-Q)/(2*m))) !atan((2*m*kf)/(m*m + 0.25_dp*Q*Q - kf2))
       aux4 = 3*m*m*(4*kf2 - Q*Q + 4*m*m)
       aux5 = log((m*m + (kf - Q*0.5_dp)**2)/(m*m + (kf + Q*0.5_dp)**2))
       I2 = aux1*(aux2 + aux3 + aux4*aux5)
 
    end function
 
    ! Fermi momentum from density for SNM (rho_snm = rhop+rhon)
    function lda_kf_snm(rho) result(kf)
       use pnfam_constants, only : pi
       use hfb_solution, only : nghl
       implicit none
       real(dp), dimension(nghl), intent(in) :: rho
       real(dp), dimension(nghl) :: kf
 
       kf = 0
       kf = (1.5_dp*pi*pi*rho)**(1.0_dp/3.0_dp)
 
    end function
 
    ! Fermi momentum from density for ASNM (twice that of SNM, rho = rhop or rhon)
    function lda_kf_asnm(rho) result(kf)
       use pnfam_constants, only : pi
       use hfb_solution, only : nghl
       implicit none
       real(dp), dimension(nghl), intent(in) :: rho
       real(dp), dimension(nghl) :: kf
 
       ! Kf is twice Kf of SNM
       kf = 0
       kf = (3.0_dp*pi*pi*rho)**(1.0_dp/3.0_dp)
 
    end function
 
    ! Density dependent coeff of (sigma tau) from effective 2BC
    function tbc_nmlda_da1(rhon, rhop, Q_per_kf, snm) result(da1)
       use pnfam_constants, only : pi, Mpi, hbarc, Fpi, Mn, gA
       use hfb_solution, only : nghl
       use type_extfield_2bc, only : c3, c4, use_p
       implicit none
 
       real(dp), dimension(nghl), intent(in) :: rhon, rhop
       real(dp), intent(in) :: Q_per_kf
       logical, intent(in) :: snm
 
       real(dp) :: caux0, caux3, caux4, ki, ri, Qi, I1, I2
       real(dp), dimension(nghl) :: da1, rho, kf, Q
       integer :: it, i
 
       da1 = 0; rho = 0; kf = 0; Q = 0
 
       ! Factor out gA (note no factor of 1/2 as in TBC module)
       caux0 = (hbarc*hbarc*hbarc)/(Mn*Fpi*Fpi)
       caux3 = -1.0_dp/3.0_dp*c3
       if (use_p) then
          caux3 = caux3 + 1.0_dp/12.0_dp
       end if
       caux4 =  1.0_dp/3.0_dp*(c4 + 0.25_dp)
 
       do it=1,3
          ! Handle loop
          if (snm) then
              if (it/=3) cycle
          else
              if (it==3) cycle
          end if
 
          ! Get isospin specific quantities
          if (it==1) then
             rho = rhon
             kf  = lda_kf_asnm(rho)
          else if (it==2) then
             rho = rhop
             kf  = lda_kf_asnm(rho)
          else
             rho = rhop + rhon
             kf  = lda_kf_snm(rho)
          end if
 
          ! Compute contributions
          do i=1,nghl
             ki = kf(i)
             ri = rho(i)
             Qi = ki*Q_per_kf
             if (abs(Qi) < 1e-3) then
                I1 = tbc_nmlda_I0(ki)
                I2 = I1
             else
                ! I haven't derived this yet.
                !  1. There might be "cross terms", like (I1(kn, Qn) + I1(kn, Qp))*0.5_dp and vice versa
                !  2. I might need the pseudo scalar terms. p=0, but P=Q/=0...
                I1 = tbc_nmlda_I1(ki,Qi) 
                I2 = tbc_nmlda_I2(ki,Qi)
             end if
             da1(i) = da1(i) + caux0*ri*(caux4*(3.0_dp*I2 - I1) + caux3*I1)
          end do
       end do
 
    end function
 
    !tbc_nmlda_I0 function from above, but takes in/returns a full vector of size nghl.
    function tbc_nmlda_I0_vec(kf) result(I0)
       use pnfam_constants, only : pi, Mpi, hbarc
       use hfb_solution, only : nghl
       implicit none
       real(dp), dimension(nghl), intent(in) :: kf
       real(dp), dimension(nghl) :: kf3, kf2, I0
       real(dp) :: m
 
       m = Mpi/hbarc
       kf2 = kf*kf
       kf3 = kf*kf*kf
 
       I0 = 1.0_dp - 3.0_dp*m*m/(kf2) + 3.0_dp*m*m*m/(kf3)*atan(kf/m)
    end function
 
    !3/5/24: 2 body correction term for the rsL operators (expansion of e^-iqr term with two-body Gamow-Teller term).
    !See arXiv:1403.7860v1 "Chiral Two-Body Currents and Neutrinoless Double-Beta Decay" eq 7, "Large-scale nuclear structure calculations for spin-dependent WIMP scattering" eq 14
    function forbidden_2bc_rsL() result (iout)
       use hfb_solution, only : nghl
       use hfb_solution, only : hfb_density_coord
       use type_extfield_2bc, only : caux,c3,c4, cd
       use pnfam_constants, only : Mpi, hbarc, IT_ISOSCALAR
       implicit none
       real(dp), dimension(nghl) :: rho, kf, i0, iout
       real(dp) :: m, ca, ch
       character(200) :: st
       ! Constants
       m = Mpi/hbarc
       ca = 2*caux
       ch = 1._dp/3._dp*(2*c4 - c3 + 0.5_dp)
 
       ! Densities
       call hfb_density_coord(IT_ISOSCALAR, rho)
       kf = lda_kf_snm(rho)
       i0 = tbc_nmlda_I0_vec(kf)
       iout = ca*rho*(cd+ch*i0)
 
       !testing
       !write(st,'(a20,f10.5,a20,f10.5,a20)') 'iout value for kf ', kf(10), 'iout', iout(10), 'in gt2bc.'
       !call writelog(st)
    end function
 
    function dme_u00(kf) result(u00)
       use hfb_solution, only : nghl
       use pnfam_constants, only : hbarc, Mpi
       implicit none
       real(dp), dimension(nghl), intent(in) :: kf
       real(dp), dimension(nghl) :: u00
       real(dp) :: m, limitAtZero
       integer :: loopvar
 
       m = Mpi/hbarc
       u00 = 1.5_dp/(kf*kf)*(1 - m/kf*atan(2*kf/m) + m*m/(4.0_dp*kf*kf)*log(1 + 4*kf*kf/(m*m)))
 
       limitAtZero = 1.0_dp / (m*m) !this function converges to 1/M^2 at zero, but oscillates due to roundoff error around kf = 0.2 MeV.
       do loopvar=1,nghl  
          if (kf(loopvar) < 0.001) then
             u00(loopvar) = limitAtZero
          endif
       end do  
    end function
 
    function dme_u02(kf) result(u02)
       use hfb_solution, only : nghl
       use pnfam_constants, only : hbarc, Mpi
       implicit none
       real(dp), dimension(nghl), intent(in) :: kf
       real(dp), dimension(nghl) :: u02
       real(dp) :: m
       integer :: loopvar
 
       m = Mpi/hbarc
       u02 = 0.75_dp/(kf*kf)*(log(1.0_dp + 4.0_dp*kf*kf/(m*m)) - 4.0_dp*kf*kf/(4.0_dp*kf*kf+m*m))
       do loopvar=1,nghl  
          if (kf(loopvar) < 0.00005) then !0.01 MeV -> 5e-5 fm^-1 is a reasonable cutoff.
             u02(loopvar) = 0
          endif
       end do  
    end function
 
    !need a separate function to handle U_0^2/k_F^2. 
    function dme_u02_kf2(kf) result(u02)
       use hfb_solution, only : nghl
       use pnfam_constants, only : hbarc, Mpi
       implicit none
       real(dp), dimension(nghl), intent(in) :: kf
       real(dp), dimension(nghl) :: u02
       real(dp) :: m, limitAtZero
       integer :: loopvar
 
       m = Mpi/hbarc
       u02 = 0.75_dp/(kf*kf*kf*kf)*(log(1.0_dp + 4.0_dp*kf*kf/(m*m)) - 4.0_dp*kf*kf/(4.0_dp*kf*kf+m*m))
       limitAtZero = 6.0_dp / (m**4)
 
       do loopvar=1,nghl  
          if (kf(loopvar) < 0.001) then !reasonable cutoff, see plot.
             u02(loopvar) = limitAtZero
          endif
       end do  
    end function
 
    !U_2^2 function used in vector current. Since it's divided by kF^2 and multiplied by rho, include all the prefactors.
    function dme_u22(kf) result(u22)
       use hfb_solution, only : nghl
       use pnfam_constants, only : hbarc, Mpi, pi
       implicit none
       real(dp), dimension(nghl), intent(in) :: kf
       real(dp), dimension(nghl) :: u22
       real(dp) :: m
 
       m = Mpi/hbarc
       !fixed definition from Dr Engel's DME currents paper. No oscillation
       u22 = -1.0_dp * kf / (3.0_dp*pi*pi) * (1.0_dp + 8.0_dp*kf*kf/(4.0_dp*kf*kf+m*m) - m*m/(4.0_dp*kf*kf)*log(1.0_dp + 4.0_dp*kf*kf/(m*m)))
       !u22 = 3.0_dp*kf*kf/(8.0_dp*m*m*(4.0_dp*kf*kf+m*m))*(48.0_dp + 20.0_dp*m*m/(kf*kf) &
       !      + m*(4.0_dp*kf*kf + m*m)/(kf*kf*kf)*(m/kf*log(m*m/(4.0_dp*kf*kf + m*m)) - 8.0_dp*atan(2*kf/m)))
 
    end function
 
    function dme_exc() result(rf)
       use hfb_solution, only : nghl
       use hfb_solution, only : hfb_density_coord, d2rho, tau
       use type_extfield_2bc, only : caux, d1,d2,c3,c4
       use pnfam_constants, only : Mpi, hbarc, IT_ISOSCALAR, IT_PROTON, IT_NEUTRON
       implicit none
       real(dp), dimension(nghl) :: rho, kf, u02, u00, u02_kf2, rf1, rf2, rf, d2r0, tau0
       real(dp) :: m, ca, ch, cc
 
       ! Constants
       m = Mpi/hbarc
       ca = 2*caux
       ch = 1._dp/3._dp*(2*c4 - c3 + 0.5_dp)
       cc = -d1 + 2*d2
 
       ! Densities
       call hfb_density_coord(IT_ISOSCALAR, rho)
       d2r0 = d2rho(:,IT_NEUTRON)+d2rho(:,IT_PROTON)
       tau0 = tau(:,IT_NEUTRON)+tau(:,IT_PROTON)
 
       ! DME integrals
       kf = lda_kf_snm(rho) ! 1/fm
       u00 = dme_u00(kf) ! fm^2
       u02 = dme_u02(kf) ! fm^2
       u02_kf2 = dme_u02_kf2(kf)
 
       ! Full expression
       rf1 = ca*(ch*(1 - m*m*u00 - m*m*u02*0.1_dp) + cc)*rho
       rf2 = ca*ch*m*m/(6)*u02_kf2*(0.25_dp*d2r0 - tau0)
       !rf2 = ca*ch*m*m/(6*kf*kf)*u02*(0.25_dp*d2r0 - tau0)
       rf = rf1 - rf2
 
    end function
    
    !no reason to have this in a separate function.
    function forbidden_2bc_P_old(kf, P) result(iout) !compute the integral functions of k_F, P, and Mpi.
       use hfb_solution, only : nghl
       use pnfam_constants, only : hbarc, Mpi
       implicit none
       real(dp), dimension(nghl), intent(in) :: kf
       real(dp), intent(in) :: P
       real(dp), dimension(nghl) :: iout
       real(dp) :: m
 
       m = Mpi/hbarc
       !multiply i0 by 0.5, i2 by 2.
       iout = 0.5_dp / P**3 * (-P * kf * (kf*kf + m*m + P*P) &
       + 0.25_dp * (((m*m + P*P)**2 + kf*kf * (kf*kf - 2 * P*P + 2*m*m)) * log((m*m + (P + kf)**2)/(m*m + (P - kf)**2)))) &
       + 1.0_dp / (8.0_dp * P**3) * ((-kf*kf*(kf*kf + 4.0_dp*m*m - 2.0_dp*P*P) &
       - (3.0_dp*m**4 + 4.0_dp*m*m*P*P + P**4)) * log((m*m + (P + kf)**2)/(m*m + (P - kf)**2)) + 4.0_dp*kf*P*(kf*kf+3.0_dp*m*m+P*P) )
    end function
 
    function forbidden_2bc_P() result (rf)
       use hfb_solution, only : nghl, hfb_density_coord
       use pnfam_constants, only : Mpi, Fpi, Mn, hbarc, gA, pi, IT_ISOSCALAR
       implicit none
       real(dp), dimension(nghl) :: rho, kf, rf, iout
       real(dp) :: P, m
       !character(len=200) :: st
         ! Constants
       !P = 50.0_dp / hbarc !determine an acceptable value for P. 50 MeV? convert to fm^-1. 
       P = sqrt(1.2)*1.361 / 2.0 
       !P = 0.5477_dp * MAXVAL(kf)
       m = Mpi/hbarc
         ! DME integrals
       call hfb_density_coord(IT_ISOSCALAR, rho)
       kf = lda_kf_snm(rho) ! 1/fm
       
       !iout = 0.5_dp / P**3 * (-P * kf * (kf*kf + m*m + P*P) &
       !+ 0.25_dp * (((m*m + P*P)**2 + kf*kf * (kf*kf - 2 * P*P + 2*m*m)) * log((m*m + (P + kf)**2)/(m*m + (P - kf)**2)))) &
       !+ 1.0_dp / (8.0_dp * P**3) * ((-kf*kf*(kf*kf + 4.0_dp*m*m - 2.0_dp*P*P) &
       !- (3.0_dp*m**4 + 4.0_dp*m*m*P*P + P**4)) * log((m*m + (P + kf)**2)/(m*m + (P - kf)**2)) + 4.0_dp*kf*P*(kf*kf+3.0_dp*m*m+P*P) )
 
       !simplified version - see Backup docs for combining 2I_2 + 1/2 I.
       iout = 1.0_dp / (P**3) * (P*M*M*kf - 0.25_dp*M*M*(M*M+P*P+kf*kf)*log((m*m + (P + kf)**2)/(m*m + (P - kf)**2)))
       !iout = forbidden_2bc_P(kf, P) 
       rf = gA*gA*hbarc*Mn/ (4*pi*pi*Fpi*Fpi) * iout !Full correction is gA^2 / (4pi^2 Fpi^2).
       !Multiply out 1/correction factor: -Mn / hbarc.Since we would multiply by two factors of hbarc for the Fpi^2 in denominator, 
       !we end up with Mn * hbarc in numerator. don't need to worry about factoring out a gA since the original P operator has no gA.
       !write(st,'(a13,f16.3,a13, f16.3,a21)') 'max of rf  : ', MAXVAL(rf), 'min of rf  : ', MINVAL(rf), 'in forbidden 2bc exc.'
       !call writelog(st)
    end function
         
    function forbidden_2bc_PS0() result (rf)
       use hfb_solution, only : nghl, hfb_density_coord
       use pnfam_constants, only : Mpi, Fpi, Mn, hbarc, gA, pi, IT_ISOSCALAR
       implicit none
       real(dp), dimension(nghl) :: rho, kf, rf, iout
       real(dp) :: P, m
         ! Constants
       !divide by zero error when using smaller P values..
       P = 275.0 / hbarc
       !P = sqrt(1.2)*1.361 / 2.0 !determine an acceptable value for P. Try using the Fermi gas mean value? (P = sqrt(6/5) * kF^2 / 2)
       !P = 0.5477_dp * MAXVAL(kf)
       m = Mpi / hbarc
         ! DME integrals
       call hfb_density_coord(IT_ISOSCALAR, rho)
       kf = lda_kf_snm(rho) ! 1/fm
 
       iout = 1.0_dp / P**3 * (-P * kf * (kf*kf + m*m + P*P) &
       + 0.25_dp * (((m*m + P*P)**2 + kf*kf * (kf*kf - 2*P*P + 2*m*m)) * log((m*m + (P + kf)**2)/(m*m + (P - kf)**2))))
       rf = -hbarc*Mn/ (8*pi*pi*Fpi*Fpi) * iout
       !factor out the gA since the PS0 ope      !inquire(file="rf_ps0.txt", exist=exists)
       !Reference the dme exc above where gA is factored out of the caux factor.
       !Prefactor here is also -1/Mn, so need to multiply by -Mn. See notes on the signs of the gA, but since this term enhances the one body piece,
       !we can determine the sign accordingly (since iout < 0)
    end function
 
    !U_11 functions in dme current.
    function dme_u11(kf) result(u11)
       use hfb_solution, only : nghl
       use pnfam_constants, only : hbarc, Mpi
       implicit none
       real(dp), dimension(nghl), intent(in) :: kf
       real(dp), dimension(nghl) :: u11
       real(dp) :: m
 
       m = Mpi/hbarc
       u11 = 1.5_dp / (kf * kf) * (atan(2.0_dp*kf/m) - 2.0_dp*kf*m/(4.0_dp*kf*kf + m*m))
    end function
 
    !separate function for U_11/kF^2
    function dme_u11_kf2(kf) result(u11)
       use hfb_solution, only : nghl
       use pnfam_constants, only : hbarc, Mpi
       implicit none
       real(dp), dimension(nghl), intent(in) :: kf
       real(dp), dimension(nghl) :: u11
       real(dp) :: m
       integer :: loopvar
 
       m = Mpi/hbarc
       u11 = 1.5_dp / (kf*kf*kf*kf) * (atan(2.0_dp*kf/m) - 2.0_dp*kf*m/(4.0_dp*kf*kf + m*m))
       do loopvar=1,nghl  
          if (kf(loopvar) < 0.001) then !reasonable cutoff.
             u11(loopvar) = 8.0_dp/(m*m*m*kf(loopvar)) !approximation near zero.
          endif
       end do  
 
    end function
 
    !U_1^0 function. 
    function dme_u10(kf) result(u10)
       use hfb_solution, only : nghl
       use pnfam_constants, only : hbarc, Mpi
       implicit none
       real(dp), dimension(nghl), intent(in) :: kf
       real(dp), dimension(nghl) :: u10
       real(dp) :: m, limitAtZero
       integer :: loopvar
 
       m = Mpi/hbarc
       u10 = 1.5_dp / (kf * kf) * (1.0_dp - m*m/(4.0_dp*kf*kf)*log(1.0_dp+4.0_dp*kf*kf/(m*m)))
       limitAtZero = 3.0_dp / (m*m)
       do loopvar=1,nghl  
          if (kf(loopvar) < 0.001) then !reasonable cutoff.
             u10(loopvar) = limitAtZero !approximation near zero.
          endif
       end do  
    end function
 
    !U_1,-1 function. Since this results in a singularity near 0, it includes an additional factor of kf^3 from multiplying by rho in the axial charge DME>
    function dme_u1_1(kf) result(u1_1)
       use hfb_solution, only : nghl
       use pnfam_constants, only : hbarc, Mpi
       implicit none
       real(dp), dimension(nghl), intent(in) :: kf
       real(dp), dimension(nghl) :: u1_1
       real(dp) :: m
       integer :: loopvar
 
       m = Mpi/hbarc
       u1_1 = m / (2.0_dp) * (1.0_dp + (2.0_dp*kf/m - 1.5_dp*m/kf)*atan(2.0_dp*kf/m) + m*m/(2.0_dp*kf*kf)*log(1.0_dp + 4.0_dp*kf*kf/(m*m)))
       do loopvar=1,nghl  
          if (kf(loopvar) < 0.0001) then !good cutoff to set to 0.
             u1_1(loopvar) = 0
          endif
       end do  
    end function
 
    !U_1^2 function for vector dme.
    function dme_u12(kf) result(u12)
       use hfb_solution, only : nghl
       use pnfam_constants, only : hbarc, Mpi
       implicit none
       real(dp), dimension(nghl), intent(in) :: kf
       real(dp), dimension(nghl) :: u12
       real(dp) :: m
       integer :: loopvar
       m = Mpi/hbarc
       u12 = 3.0_dp / (4.0_dp*kf*kf+m*m) * ((4.0_dp*kf*kf-m*m)/(4.0_dp*kf*kf+m*m) + (1.0_dp + m*m/(4.0_dp*kf*kf))*log(1.0_dp+4.0_dp*kf*kf/(m*m)))
       do loopvar=1,nghl  
          if (kf(loopvar) < 0.0001) then !reasonable cutoff.
             u12(loopvar) = 0 !approximation near zero.
          endif
       end do  
    end function
 
    !need a separate function to handle U_1^2/k_F^2. 
    function dme_u12_kf2(kf) result(u12)
       use hfb_solution, only : nghl
       use pnfam_constants, only : hbarc, Mpi
       implicit none
       real(dp), dimension(nghl), intent(in) :: kf
       real(dp), dimension(nghl) :: u12
       real(dp) :: m, c1, c2
       integer :: loopvar
       m = Mpi/hbarc
       u12 = 3.0_dp / (kf*kf*(4.0_dp*kf*kf+m*m)) * ((4.0_dp*kf*kf-m*m)/(4.0_dp*kf*kf+m*m) + (1.0_dp + m*m/(4.0_dp*kf*kf))*log(1.0_dp+4.0_dp*kf*kf/(m*m)))
       c1 = 30.0_dp / (m**4)
       c2 = -224.0_dp / (m**6)
       do loopvar=1,nghl  
          if (kf(loopvar) < 0.001) then !reasonable cutoff.
             u12(loopvar) = c1 + c2 * kf(loopvar)**2 !use Taylor series approx: u12 = c1 + c2 * kf^2
          endif
       end do  
    end function
 
    !U_2^1 function for vector dme.
    function dme_u21(kf) result(u21)
       use hfb_solution, only : nghl
       use pnfam_constants, only : hbarc, Mpi
       implicit none
       real(dp), dimension(nghl), intent(in) :: kf
       real(dp), dimension(nghl) :: u21
       real(dp) :: m
 
       m = Mpi/hbarc
       u21 = 1.5_dp / (kf*m) * (1.0_dp + (2.0_dp*kf/m - 0.5_dp*m/kf) * atan(2.0_dp*kf/m))
    end function
 
    !U_2^3 function for vector dme.
    function dme_u23(kf) result(u23)
       use hfb_solution, only : nghl
       use pnfam_constants, only : hbarc, Mpi
       implicit none
       real(dp), dimension(nghl), intent(in) :: kf
       real(dp), dimension(nghl) :: u23
       real(dp) :: m
 
       m = Mpi/hbarc
       u23 = 1.5_dp / (m*m) * (3.0_dp*atan(2.0_dp*kf/m) - 2.0_dp * (4.0_dp*m*kf*kf*kf + 3.0_dp*m*m*m*kf)/(4.0_dp*kf*kf + m*m)**2)
    end function
 
    !U_2^4 function for vector dme. Includes additional factors from vector current definition.
    function dme_u24(kf) result(u24)
       use hfb_solution, only : nghl
       use pnfam_constants, only : hbarc, Mpi
       implicit none
       real(dp), dimension(nghl), intent(in) :: kf
       real(dp), dimension(nghl) :: u24
       real(dp) :: m
       m = Mpi / hbarc
       u24 = kf*kf*(80.0_dp*kf**4 + 40.0_dp*kf*kf*m*m - 3.0_dp*m**4)/(4.0_dp*kf*kf+m*m)**3 + 0.75_dp*log(1.0_dp+4.0_dp*kf*kf/(m*m))
 
    end function
 
    !U_2^4 function divided by kf^4. Use taylor series expansion near kf=0. Includes factors from vector DME expression.
    function dme_u24_kf4(kf) result(u24)
       use hfb_solution, only : nghl
       use pnfam_constants, only : hbarc, Mpi
       implicit none
       real(dp), dimension(nghl), intent(in) :: kf
       real(dp), dimension(nghl) :: u24
       real(dp) :: m, c1, c2
       integer :: loopvar
       m = Mpi / hbarc
       u24 = 1.0_dp / (6.0_dp *kf*kf*kf*kf) * (kf*kf*(80.0_dp*kf**4 + 40.0_dp*kf*kf*m*m - 3.0_dp*m**4)/(4.0_dp*kf*kf+m*m)**3 + 0.75_dp*log(1.0_dp+4.0_dp*kf*kf/(m*m)))
       c1 = 35.0_dp / (3.0_dp * m**4)
       c2 = -112.0_dp / (m**6) !Taylor expansion of U_2^4 / kf**4: c1 + c2 kf^2.
       do loopvar=1,nghl  
          if (kf(loopvar) < 0.001) then !reasonable cutoff.
             u24(loopvar) = c1 + c2*kf(loopvar)**2 !use Taylor series approximation
          endif
       end do 
    end function
 
    !DME current for the axial charge term.
    function dme_axial() result(rf)
       use hfb_solution, only : nghl
       use hfb_solution, only : hfb_density_coord, d2rho, tau
       use pnfam_constants, only : Mpi, Mn, Fpi, pi, hbarc, IT_ISOSCALAR, IT_PROTON, IT_NEUTRON
       implicit none
       real(dp), dimension(nghl) :: rho, kf, u11, u11_kf2, u1_1, rf, d2r0, tau0
 
       ! Densities
       call hfb_density_coord(IT_ISOSCALAR, rho)
       d2r0 = d2rho(:,IT_NEUTRON)+d2rho(:,IT_PROTON)
       tau0 = tau(:,IT_NEUTRON)+tau(:,IT_PROTON)
 
       ! DME integrals
       kf = lda_kf_snm(rho) ! 1/fm
       u11 = dme_u11(kf) ! fm^2
       u11_kf2 = dme_u11_kf2(kf)
       u1_1 = dme_u1_1(kf) ! fm^2
 
       !adjust for the additional factor of kf included in U_1^-1 calculation.
       !according to Dr Engel, in constant density nuclear matter, only the first term involving U_1^{-1} is nonzero since the term with tau0 cancels that with 0.6
       rf = hbarc*Mn/(6.0_dp*Fpi*Fpi) * (u1_1 *(2.0_dp/(3.0_dp*pi*pi)) + 0.6_dp*u11*rho + u11_kf2*(0.25_dp*d2r0 - tau0))
    end function
 
    !DME current for the vector current term.
    function dme_vector() result(rf)
       use hfb_solution, only : nghl, y
       use hfb_solution, only : hfb_density_coord, d2rho, tau
       use pnfam_constants, only : Mpi, Mn, Fpi, gA, hbarc, IT_ISOSCALAR, IT_PROTON, IT_NEUTRON, pi
       implicit none
       real(dp), dimension(nghl) :: rho, kf, u10, u12, u12_kf2, u22, u24, u24_kf4, rf, d2r0, tau0, r
       real(dp) :: m
       character(len=200) :: st
       logical :: exists
       integer :: loopvar
 
       m = Mpi/hbarc
       ! Define r(:) = 1/(1/rho)
       r(:) = 1.0_dp/y(:)
       ! Densities
       call hfb_density_coord(IT_ISOSCALAR, rho)
       d2r0 = d2rho(:,IT_NEUTRON)+d2rho(:,IT_PROTON)
       tau0 = tau(:,IT_NEUTRON)+tau(:,IT_PROTON)
 
       rf = 0
 
       ! DME integrals
       kf = lda_kf_snm(rho) ! 1/fm
       u10 = dme_u10(kf)
       u12 = dme_u12(kf)
       u22 = dme_u22(kf)
       u24 = dme_u24(kf)
       u12_kf2 = dme_u12_kf2(kf)
       u24_kf4 = dme_u24_kf4(kf)
 
       !extra factor of kf already included in U_1^-1.
       !rf = gA*gA*Mn*hbarc/(4.0_dp * Fpi*Fpi) * ( (u10 + 0.1_dp*u12 - m*m*u24/(30.0_dp*kf*kf))*rho + u22 &
       
       !rf = gA*gA*Mn*hbarc/(4.0_dp * Fpi*Fpi) * ( (u10 + 0.1_dp*u12)*rho + u22 - u24*kf/(60.0_dp*pi*pi))
       !rf =  (0.25_dp*d2r0 - tau0)
       
       rf = gA*gA*Mn*hbarc/(4.0_dp * Fpi*Fpi) * ( (u10 + 0.1_dp*u12)*rho + u22 - u24*kf/(15.0_dp*pi*pi)&
       + (u12_kf2/6.0_dp - u24_kf4) * (0.25_dp*d2r0 - tau0))
 
       !constant density case: only u10, u22 terms survive, other terms cancel.
       !rf = gA*gA*Mn*hbarc/(4.0_dp * Fpi*Fpi) * (u10*rho + u22)
       !rf = gA*gA*Mn*hbarc/(4.0_dp * Fpi*Fpi) * ((u10 + 0.1_dp*u12)*rho + u22 - u24*kf/(60.0_dp*pi*pi))
       !rf = gA*gA*Mn*hbarc/(4.0_dp * Fpi*Fpi) * ((u10 + 0.1_dp*u12)*rho - kf * u24/(60.0_dp*pi*pi) + u22)
       !rf = gA*gA*Mn*hbarc/(4.0_dp * Fpi*Fpi) * ((u12_kf2/6.0_dp - u24_kf4))
 
       !try writing output to a file. check if exists first - only write once
       !inquire(file="rf_part1_dme_updated.txt", exist=exists)
       !if (exists .eqv. .false.) then
       !   open(12, file = 'rf_part1_dme_updated.txt', status = 'new')
       !   do loopvar=1,nghl  
       !      write(12,*) rf(loopvar)
       !   end do  
       !   close(12)
       !endif
       !write(st,'(a20,f10.5,a20,f10.5,a20)') 'first item in u10', u10(1), 'first item in kf', kf(1), 'in dme.'
       !call writelog(st)
       ! Full expression - multiply by 1/prefactor = -Mn.
    end function
 end module pnfam_extfield