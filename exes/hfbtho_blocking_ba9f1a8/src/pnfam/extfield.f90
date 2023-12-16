!------------------------------------------------------------------------------
! extfield.f90
!
! Charge-changing external fields for the pnFAM code.
!
! M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
!------------------------------------------------------------------------------
module extfield
   use logger
   implicit none
   private
   
   integer,  parameter :: dp = kind(1d0)
   real(dp), parameter :: pi = 3.141592653589793238462643_dp

   character(len=1), dimension(0:1), parameter :: cbeta = (/'+', '-'/)
   
   public :: init_external_field
   public :: number_crossterms
   public :: setup_crossterms
   public :: print_crossterms
   public :: init_fam_mapping
   public :: deallocate_extfield
   public :: deallocate_crossterms
   
contains

   !----------------------------------------------------------------------------
   ! Take in a label and K-projection and create an external_field type.
   ! This is a more consolidated routine than previously: now we assign the
   ! various quantities and ALSO set up the FAM s.p. matrix itself. It is
   ! accessible via op%mat%elem (cf. f%elem previously).
   !----------------------------------------------------------------------------
   subroutine init_external_field(label, k, op)
      use constants, only : translate_uppercase
      use external_field_type
      implicit none

      character(len=*),     intent(in)    :: label
      integer,              intent(in)    :: k
      type(external_field), intent(inout) :: op

      integer :: l
      
      ! Operator label
      l = len_trim(label)
      op%label = label(1:l-1)
      
      ! Normalize to uppercase
      call translate_uppercase(op%label)

      ! Beta minus
      if (label(l:l) == '-') then
         op%beta_minus = .true.
      else
         op%beta_minus = .false.
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
            stop
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
            write(*,'(A)') "ERROR: Rank not implemented for the operator"
            write(*,'(A)') "This should not happen, fix it in extfield.f90"
            stop
      end select
      
      ! Matrix structure
      call init_fam_mapping(op)
      call ext_field_operator(op)
      
   end subroutine init_external_field


   !----------------------------------------------------------------------------
   ! Construct the chosen charge-changing operator in the single-particle basis.
   !----------------------------------------------------------------------------
   subroutine ext_field_operator(op)
      use blockmatrix_type
      use external_field_type
      use hfbtho_basis, only : db, isstart, nl, ns, nghl, wf, wfdr, wfdz, y, z
      use iso_fortran_env, only : error_unit

      implicit none
      type(external_field), intent(inout) :: op
      
      ! These are assigned from the operator
      integer :: K
      logical :: pty
      character(len=80) :: label

      integer :: ipt, i1, i2, ix1, ix2, ibx1, ibx2, nd1, nd2
      integer :: xl1, xl2, xs1, xs2
      real(dp), dimension(nghl) :: wf_1, wf_2, dr_wf_2, dz_wf_2, r

      ! Define r(:) = 1/(1/rho)
      r(:) = 1.0_dp/y(:)
      
      ! Quantities taken from the operator
      K = op%k
      pty = op%parity_even
      label = op%label

      ! Basic tests
      call assert_parity(label=label, pty=pty)
      call assert_K_in_range(label=label, K=K)

      ! Loop over the FAM block structure to do the calculation
      ipt = 0;  op%mat%elem(:) = 0
      do ibx1 = 1, nb

         ! Determine the range of particle indices in this block
         ibx2 = op%mat%ir2c(ibx1)
         if (ibx2 == 0) cycle
         
         nd1 = db(ibx1) ; nd2 = db(ibx2)
         if (nd1 == 0 .or. nd2 == 0) cycle
         
         do i2=1, nd2
            ! Calculate the true index of the state |2>
            ix2 = i2 + isstart(ibx2) - 1

            ! Wave functions
            wf_2 = wf(:,ix2)
            dr_wf_2 = wfdr(:,ix2)
            dz_wf_2 = wfdz(:,ix2)

            ! Quantum numbers
            xl2 = nl(ix2);  xs2 = ns(ix2)

            do i1=1, nd1
               ! Calculate the true index of the state <1|
               ix1 = i1 + isstart(ibx1) - 1

               ! Wave functions
               wf_1 = wf(:,ix1)

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
                     stop


                  ! Fermi (1)
                  ! ------------------------------
                  case ('F')
                     if ((xl1 == xl2) .and. (xs1 == xs2)) then
                        op%mat%elem(ipt) = dot_product(wf_1(:), wf_2(:))
                     end if


                  ! Gamow-Teller (sigma_K)
                  ! ------------------------------
                  case ('GT')
                     select case (K)
                        case (0)
                           if ((xl1 == xl2) .and. (xs1 == xs2)) then
                              op%mat%elem(ipt) = xs1*dot_product(wf_1(:), wf_2(:))
                           end if
                        case (1,-1)
                           if ((xl1 == xl2) .and. (xs1 == xs2 + 2*K)) then
                              op%mat%elem(ipt) = -K*sqrt(2.0_dp)*dot_product(wf_1(:), wf_2(:))
                           end if
                     end select


                  ! R (Sqrt(4\pi/3)*rY_1K = r_K)
                  ! ------------------------------
                  case ('R')
                     select case (K)
                        case (0)
                           if ((xl1 == xl2) .and. (xs1 == xs2)) then
                              op%mat%elem(ipt) = dot_product(wf_1(:), z(:)*wf_2(:))
                           end if
                        case (1,-1)
                           if ((xl1 == xl2 + K) .and. (xs1 == xs2)) then
                              op%mat%elem(ipt) = -K/sqrt(2.0_dp)*dot_product(wf_1(:), r(:)*wf_2(:))
                           end if
                     end select


                  ! P (p_K/i)
                  ! ------------------------------
                  case ('P')
                     select case (K)
                        case (0)
                           if ((xl1 == xl2) .and. (xs1 == xs2)) then
                              op%mat%elem(ipt) = -dot_product(wf_1(:), dz_wf_2(:))
                           end if
                        case (1,-1)
                           if ((xl1 == xl2 + K) .and. (xs1 == xs2)) then
                              op%mat%elem(ipt) = K/sqrt(2.0_dp)*(dot_product(wf_1(:), dr_wf_2(:))   &
                                     - K*xl2*dot_product(wf_1(:), y(:)*wf_2(:)))
                           end if
                     end select


                  ! RS0 ([rY_1 x sigma]_00/Y_00)
                  ! ------------------------------
                  case ('RS0')
                     if ((xl1 == xl2) .and. (xs1 == xs2)) then
                        op%mat%elem(ipt) = -xs1*dot_product(wf_1(:), z(:)*wf_2(:))
                     else if ((xl1 == xl2 + 1) .and. (xs1 == xs2 - 2)) then
                        op%mat%elem(ipt) = -dot_product(wf_1(:), r(:)*wf_2(:))
                     else if ((xl1 == xl2 - 1) .and. (xs1 == xs2 + 2)) then
                        op%mat%elem(ipt) = -dot_product(wf_1(:), r(:)*wf_2(:))
                     end if


                  ! RS1 ([rY_1 x sigma]_1K/Y_00)
                  ! ------------------------------
                  case ('rs1','RS1')
                     select case (K)
                        case (0)
                           if ((xl1 == xl2 - 1) .and. (xs1 == xs2 + 2)) then
                              op%mat%elem(ipt) = sqrt(3.0_dp/2.0_dp)*dot_product(wf_1(:),r(:)*wf_2(:))
                           else if ((xl1 == xl2 + 1) .and. (xs1 == xs2 - 2)) then
                              op%mat%elem(ipt) = -sqrt(3.0_dp/2.0_dp)*dot_product(wf_1(:),r(:)*wf_2(:))
                           end if
                        case (1,-1)
                           if ((xl1 == xl2) .and. (xs1 == xs2 + 2*K)) then
                              op%mat%elem(ipt) = sqrt(3.0_dp)*dot_product(wf_1(:), z(:)*wf_2(:))
                           else if ((xl1 == xl2 + K) .and. (xs1 == xs2)) then
                              op%mat%elem(ipt) = -xs1*sqrt(3.0_dp)/2.0_dp*dot_product(wf_1(:), r(:)*wf_2(:))
                           end if
                     end select


                  ! RS2 ([rY_1 x sigma]_2K/Y_00)
                  ! ------------------------------
                  case ('RS2')
                     select case (K)
                        case (0)
                           if ((xl1 == xl2) .and. (xs1 == xs2)) then
                              op%mat%elem(ipt) = xs1*sqrt(2.0_dp)*dot_product(wf_1(:), z(:)*wf_2(:))
                           else if ((xl1 == xl2 + 1) .and. (xs1 == xs2 - 2)) then
                              op%mat%elem(ipt) = -1/sqrt(2.0_dp)*dot_product(wf_1(:), r(:)*wf_2(:))
                           else if ((xl1 == xl2 - 1) .and. (xs1 == xs2 + 2)) then
                              op%mat%elem(ipt) = -1/sqrt(2.0_dp)*dot_product(wf_1(:), r(:)*wf_2(:))
                           end if
                        case (1,-1)
                           if ((xl1 == xl2) .and. (xs1 == xs2 + 2*K)) then
                              op%mat%elem(ipt) = -K*sqrt(3.0_dp)*dot_product(wf_1(:), z(:)*wf_2(:))
                           else if ((xl1 == xl2 + K) .and. (xs1 == xs2)) then
                              op%mat%elem(ipt) = -K*xs1*sqrt(3.0_dp)/2.0_dp*dot_product(wf_1(:), r(:)*wf_2(:))
                           end if
                        case (2,-2)
                           if ((xl1 == xl2 + K/2) .and. (xs1 == xs2 + K)) then
                              op%mat%elem(ipt) = sqrt(3.0_dp)*dot_product(wf_1(:), r(:)*wf_2(:))
                           end if
                     end select


                     ! PS0 (p.sigma/i)
                     ! ----------------------------------------
                     case ('PS0')
                        if ((xl1 == xl2) .and. (xs1 == xs2)) then
                           op%mat%elem(ipt) = -xs1*dot_product(wf_1(:), dz_wf_2(:))
                        else if ((xl1 == xl2 - 1) .and. (xs1 == xs2 + 2)) then
                           op%mat%elem(ipt) = -dot_product(wf_1(:), dr_wf_2(:))      &
                                  - xl2*dot_product(wf_1(:), y(:)*wf_2(:))
                        else if ((xl1 == xl2 + 1) .and. (xs1 == xs2 - 2)) then
                           op%mat%elem(ipt) = -dot_product(wf_1(:), dr_wf_2(:))      &
                                  + xl2*dot_product(wf_1(:), y(:)*wf_2(:))
                        end if

               end select ! operator_name

            end do ! <1|
         end do ! |2>
      end do ! blocks

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
               stop
            end if
         ! Odd parity - parity_even .eqv. .FALSE.
         case  ('RS0', 'RS1', 'RS2', 'R', 'P', 'PS0')
            if (pty .eqv. .true.) then
               call error_wrong_parity(label)
               stop
            end if
         case default
            call error_unknown_operator(label)
            stop
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
               stop
            end if
         ! Kmax = 1
         case  ('GT', 'R', 'P', 'RS1')
            if (abs(K) > 1) then
               call error_K_out_of_range(label, 1)
               stop
            end if
         ! Kmax = 2
         case ('RS2')
            if (abs(K) > 2) then
               call error_K_out_of_range(label, 2)
               stop
            end if
         case default
            call error_unknown_operator(label)
            stop
      end select
   end subroutine assert_K_in_range


   !----------------------------------------------------------------------------
   ! Error messages for various cases
   !----------------------------------------------------------------------------
   subroutine error_unknown_operator(name)
      implicit none
      character(len=*), intent(in) :: name

      write(*,'(/3a)') 'Error: Unknown operator "', trim(name), '" in extfield.'
   end subroutine error_unknown_operator

   subroutine error_wrong_parity(name)
      implicit none
      character(len=*), intent(in) :: name

      write(*,'(/3a)') 'Error: Incorrect parity for operator "', trim(name), '" in extfield.'
   end subroutine error_wrong_parity

   subroutine error_K_out_of_range(name, K_max)
      implicit none
      integer,          intent(in) :: K_max
      character(len=*), intent(in) :: name

      write(*,'(/3a)') 'Error: K out of range for operator "', trim(name), '" in extfield.'
      write(*,'(a, 1i1)')   'Expected |K| <= ', K_max
   end subroutine error_K_out_of_range
   
   
   !-----------------------------------------------------------------------------
   ! This subroutine for initializing the f matrix block structure for general K
   ! and parity;  The HFBTHO basis must be read in before calling this!
   !-----------------------------------------------------------------------------
   subroutine init_fam_mapping(op)
      use blockmatrix_type
      use external_field_type
      use hfbtho_basis, only : nl, ns, npar
      implicit none
      
      type(external_field), intent(inout) :: op
      
      ! These are assigned from the operator
      integer :: k     ! change in angular momentum projection
      logical :: pty   ! .true. = parity conserved, .false. = parity flipped
      
      integer :: i1, i2, ip1, ip2, i
      integer :: aux1(nb), aux2(nb)
      integer :: nbb                          ! the number of FAM blocks
      integer, allocatable :: ib1(:), ib2(:)  ! the THO basis blocks the FAM block connects <1|FAM|2>
      
      ! K and parity come from the operator itself
      k = op%k
      pty = op%parity_even
      
      ! Since there may be empty (missing) blocks in the THO, we have to count
      ! the number of blocks we'll have the brute way
      nbb = 0
      ip1 = 1
      do i1=1,nb
         ip2 = 1
         do i2=1,nb
            ! Use the first state in the block to determine if Omega and parity match
            if (2*nl(ip1)+ns(ip1)-2*nl(ip2)-ns(ip2) == 2*k .and. &
               ((npar(ip1) == npar(ip2)) .eqv. pty)) then
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
         write(*,*) "The K and parity are outside the model space. Aborting."
         stop
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
      use external_field_type
      use blockmatrix_type
      implicit none
      type(external_field), intent(inout) :: e
      
      call deallocate_blockmatrix(e%mat)
      
   end subroutine deallocate_extfield
   
   
   !----------------------------------------------------------------------------
   ! Return the number of possible cross-terms for a given K and parity.
   !----------------------------------------------------------------------------
   function number_crossterms(op)
      use external_field_type
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
   !----------------------------------------------------------------------------
   subroutine setup_crossterms(op, crossterms, n)
      use external_field_type
      implicit none
      type(external_field), intent(in) :: op
      type(external_field), dimension(:), allocatable, intent(inout) :: crossterms
      integer, intent(out) :: n
      
      integer :: nxterms, ixterm, ibeta
      
      ! Count the cross-terms and allocate space
      nxterms = number_crossterms(op=op)
      n = nxterms
      
      ! No cross-terms means no work
      if (nxterms == 0) then
         return
      ! Otherwise, allocate the space and set them up
      else
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
            call init_external_field(label='RS0'//cbeta(ibeta), k=op%k, op=crossterms(ixterm))
            crossterms(ixterm)%label = trim(op%label)//'x'//trim(crossterms(ixterm)%label)
            ! OP x PS0
            ixterm = 2
            call init_external_field(label='PS0'//cbeta(ibeta), k=op%k, op=crossterms(ixterm))
            crossterms(ixterm)%label = trim(op%label)//'x'//trim(crossterms(ixterm)%label)
         ! J = 1
         else if (nxterms == 3) then
            ! OP x R
            ixterm = 1
            call init_external_field(label='R'//cbeta(ibeta), k=op%k, op=crossterms(ixterm))
            crossterms(ixterm)%label = trim(op%label)//'x'//trim(crossterms(ixterm)%label)
            ! OP x RS1
            ixterm = 2
            call init_external_field(label='RS1'//cbeta(ibeta), k=op%k, op=crossterms(ixterm))
            crossterms(ixterm)%label = trim(op%label)//'x'//trim(crossterms(ixterm)%label)
            ! OP x P
            ixterm = 3
            call init_external_field(label='P'//cbeta(ibeta), k=op%k, op=crossterms(ixterm))
            crossterms(ixterm)%label = trim(op%label)//'x'//trim(crossterms(ixterm)%label)
         end if
      end if
   end subroutine setup_crossterms
   
   
   subroutine deallocate_crossterms(crossterms)
      use external_field_type
      implicit none
      type(external_field), dimension(:), allocatable, intent(inout) :: crossterms
      
      integer :: i
      
      if (allocated(crossterms)) then
         do i=1,3
            call deallocate_extfield(crossterms(i))
         end do
         deallocate(crossterms)
      end if
      
   end subroutine deallocate_crossterms
   
   
   !----------------------------------------------------------------------------
   ! Print the computed cross-terms to STDOUT 
   !----------------------------------------------------------------------------
   subroutine print_crossterms(labels)
      implicit none
      character(len=*), dimension(:), intent(in) :: labels
      character(len=100) :: st
      integer :: i
      
      if (size(labels) > 0) then
         st = 'Cross-terms computed: '//trim(labels(1))
         do i=2, size(labels)
            st = trim(st)//', '//trim(labels(i))
         end do
         call writelog(st)
      end if
   end subroutine print_crossterms

end module extfield
