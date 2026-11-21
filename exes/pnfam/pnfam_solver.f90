!------------------------------------------------------------------------------
!> This module contains the heart of the pnfam solver. It calls the setup routines
!> runs the iterative solver, and writes the output to file.
!>
!> @authors M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
!> @authors E.M. Ney, UNC Chapel Hill, 2018-
!> Added the GMRES solver and Lanczos solver (Dong Min Roh, LBNL)
!------------------------------------------------------------------------------

   !---------------------------------------------------------------------------
   ! The FAM equations
   ! (Ep+En-omega)X + dH^{20} = -F^{20}
   ! (Ep+En+omega)Y + dH^{02} = -F^{02}
   ! are equivalent to the linear response equation
   ! (|A  B | - omega|I  0|)|X| = -|F^{20}|
   ! (|B* A*|        |0 -I|)|Y|    |F^{02}|
   ! and the amplitudes X and Y can be found from an eigenvalue decomposition
   ! |A  B | |X  Y*| = |I  0| |X Y*| |omega  0|
   ! |B* A*| |Y  X*|   |0 -I| |Y X*| |0 -omega|
   !---------------------------------------------------------------------------

module pnfam_solver
   use pnfam_setup
   use type_extfield
   use type_bigblockmatrix

   implicit none
   integer, parameter, private :: dp = kind(1d0)

   logical, private :: use_diagonal_blocks, bminus
   
   !---------------------------------------------------------------------------
   ! 'Green_sign' and 'Greens_power' determine the sign and the power of the
   ! Greens function, respectively.
   ! Greens function is defined as
   !               Greens_sign * (Ep+En-omega)^{Greens_power}
   !
   ! IFAM solver updates X and Y by
   ! X = -(Ep+En-omega)^{-1}(dH^{20}+F^{20})
   ! Y = -(Ep+En+omega)^{-1}(dH^{02}+F^{02})
   ! So, need Greens_sign = -1 and Greens_power = -1
   !
   ! GMRES solver is assigned a preconditioner 
   ! diag((Ep+En-omega)^{-1},(Ep+En+omega)^{-1})
   ! This results in a matvec defined as
   ! X <- X + (Ep+En-omega)^{-1}dH^{20}
   ! Y <- Y + (Ep+En+omega)^{-1}dH^{02}
   ! So, need Greens_sign = 1 and Greens_power = -1
   !
   ! LANCZOS/DOS solver performs matvecs with the matrix of the linear response eq,
   ! which requires access to multiplications with (EP+En-omega) and (Ep+En+omega).
   ! So, need Greens_sign = 1 and Greens_power = 1 
   !---------------------------------------------------------------------------
#ifdef IFAM
   integer, parameter :: Greens_sign = -1
#else
   integer, parameter :: Greens_sign = 1
#endif

#if defined(LANCZOS) || defined(DOS)
   integer, parameter :: Greens_power = 1
#else
   integer, parameter :: Greens_power = -1
#endif

   ! External field
   type(bigblockmatrix), private :: Fsp, Fqp
   type(bigblockmatrix), dimension(:), allocatable, private :: Gsp, Gqp
   ! HFB solution
   type(bigblockmatrix), private :: Wn, Wp
   ! FAM matrices
   type(bigblockmatrix), private :: dHqp_re, dHqp_im, dRqp_re, dRqp_im
   type(bigblockmatrix), private :: dHsp_re, dHsp_im, dRsp_re, dRsp_im
   type(bigblockmatrix), private :: Greens_re, Greens_im, T

   ! DMR
   ! For GMRES
   type(bigblockmatrix), private :: temp_re, temp_im
   type(bigblockmatrix), private :: Init_re, Init_im
   type(bigblockmatrix), private :: Out_re, Out_im

contains
   !---------------------------------------------------------------------------
   ! Subroutine for solving the FAM equations for one value of energy
   !---------------------------------------------------------------------------
   subroutine pnfam_solve
      use pnfam_broyden
      use pnfam_txtoutput
      use pnfam_storage, only : fam_io
      implicit none
      integer :: iexit=1, ierr=0

      ! Initialize FAM inputs (as blockmatrix types, where applicable)
      ! - Logging procedures
      ! - HFB solution (E,U,V)
      ! - External field (F in SP basis)
      ! - Routine for calculating dH
      !---------------------------------------------------------------
      call setup_pnfam
      call openlog(output_txtfile, print_stdout)
      call print_main_header
      call setup_extfield
      !call print_main_header !Testing purposes
      ! Prepare ifam
      !---------------------------------------------------------------
      ! 1. Allocate Block Matices:
      !    - F_QP, G_QP, dH_QP(H20,H02,H11,H11t), dH_SP(h,d),
      !      dR_QP(X,Y,P,Q), dR_SP(rho,kap), Greens(X,Y,P,Q), T(1ff,ff)
      ! 2. Transform F to QP basis
      ! 3. Define block structures (mostly based on F)
      ! 4. Populate statmatrix T and Greens function
      !---------------------------------------------------------------
      call init_pnfam_solver

      ! Initialize fam amplitudes
      !---------------------------------------------------------------
      call fam_io(1,iexit)
      if (iexit > 0) then
         si = 1
         call set_val_bbm(dRqp_re,0d0)
         call set_val_bbm(dRqp_im,0d0)
         call set_broinout_from_dRqp(dRqp_re, dRqp_im,'i')
      else
         call set_dRqp_from_broin(dRqp_re, dRqp_im)
      end if

      ! Run the iterative solver
      !---------------------------------------------------------------
      call log_header
      
      ! Compile directives for calling ifam or gmres or lanczos
      write(st,'(1x,72("-"))'); call writelog(st)
#ifdef IFAM
      write(st,'(30x,a)') 'CALLING IFAM'; call writelog(st)
      call ifam
#elif defined(GMRES)
      write(st,'(30x,a)') 'CALLING GMRES'; call writelog(st)      
      call gmres
#elif defined(DOS)
      write(st,'(30x,a)') 'CALLING DOS'; call writelog(st)
      call dos
#elif defined(LANCZOS)
      write(st,'(30x,a)') 'CALLING LANCZOS'; call writelog(st)
      call lanczos
#endif
      write(st,'(1x,72("-"))'); call writelog(st)

      ! Write the result to footer and store result
      !---------------------------------------------------------------
      call log_result
      call write_amplitudes_table
      call fam_io(2,iexit)
      call closelog

      !call transition_density
      !if (trim(f%label) == 'RS0' .or. trim(f%label) == 'PS0') then
      !   call transition_density_0m
      !end if
      !call blocked_spin_density

   end subroutine pnfam_solve

   !---------------------------------------------------------------------------
   ! Write the converged X (F20) and Y (F02) transition amplitudes to a plain
   ! text file, labeled by the two quasiparticle indices they connect.
   !
   ! Physics reminder (Mustonen et al., PRC 90, 024308 (2014), Eqs. 8-11):
   !   F = sum_{pi,nu} F20_{pi,nu} a^+_pi a^+_nu + F02_{pi,nu} a_nu a_pi
   ! and the FAM amplitudes X_{pi,nu}(w), Y_{pi,nu}(w) solve Eq. (13). For
   ! beta- (bminus=.true.) pi runs over proton quasiparticles and nu over
   ! neutron quasiparticles (and vice-versa for beta+). X lives in the m12
   ! quadrant of dRqp (the F20 block), Y in the m21 quadrant (F02).
   !
   ! Each quasiparticle index is labeled by the block it lives in (an
   ! Omega/parity subspace) and the asymptotic quantum numbers (nr,nz,nl,2*ns)
   ! of the HO basis state at that position in the block, since the
   ! Bogoliubov transformation only mixes states within a single block.
   !---------------------------------------------------------------------------
   subroutine write_amplitudes_table
      use hfb_solution, only : nr, nz, nl, ns, npar
      implicit none

      integer :: lu
      character(len=200) :: ampfile
      character(len=7) :: row_species, col_species

      if (.not. write_amplitudes) return

      ampfile = output_ampfile
      if (trim(ampfile) == "") ampfile = "pnfam_amplitudes.dat"

      if (bminus) then
         row_species = "proton "; col_species = "neutron"
      else
         row_species = "neutron"; col_species = "proton "
      end if

      open(newunit=lu, file=trim(ampfile), status='replace', action='write')

      write(lu,'(a)') '# pnFAM quasiparticle transition amplitudes'
      write(lu,'(a,a,a1,a,i0)') '# Operator: ', trim(operator_name), beta_type, ',  K = ', operator_k
      write(lu,'(a,f12.6,a,f12.6,a)') '# EQRPA = ', real_eqrpa, ' + ', imag_eqrpa, 'i  MeV'
      write(lu,'(5a)') '# row = ', trim(row_species), ' quasiparticle (pi),  col = ', &
         trim(col_species), ' quasiparticle (nu)'
      write(lu,'(a)') '# (nr,nz,nl,2*ns,par) are the quantum numbers of the HO basis state at that&
         & position in the block; the QP mixes states only within its own block.'
      write(lu,'(a)') '# Amplitudes with |X| (or |Y|) below amplitude_cutoff are omitted.'
      write(lu,'(a)') '#'

      write(lu,'(a)') '# ==== X amplitudes (F20) ===='
      call write_one_amplitude_block(lu, dRqp_re%m12, dRqp_im%m12, 'X', row_species, col_species)

      write(lu,'(a)') '#'
      write(lu,'(a)') '# ==== Y amplitudes (F02) ===='
      call write_one_amplitude_block(lu, dRqp_re%m21, dRqp_im%m21, 'Y', row_species, col_species)

      close(lu)

      write(st,'(1x,2a)') " * Wrote transition amplitudes to file: ", trim(ampfile); call writelog(st)
      call writelog("")

   contains

      subroutine write_one_amplitude_block(lu, re_mat, im_mat, name, row_species, col_species)
         implicit none
         integer, intent(in) :: lu
         type(blockmatrix), intent(in) :: re_mat, im_mat
         character(1), intent(in) :: name
         character(7), intent(in) :: row_species, col_species

         integer :: ib, ic, ir, icx, imm, idx, iqp_row, iqp_col
         real(dp) :: re_v, im_v, amp

         write(lu,'(1x,a)') 'row_iqp  blk   st   E_row[MeV]  nr nz nl 2ns par |'// &
                             '  col_iqp  blk   st   E_col[MeV]  nr nz nl 2ns par |'// &
                             '       Re             Im            |amp|'
         write(lu,'(1x,a)') repeat('-',120)

         do ib=1, nb
            ic = re_mat%ir2c(ib)
            if (ic == 0) cycle
            imm = re_mat%ir2m(ib) - 1
            do icx=1, db(ic)
               do ir=1, db(ib)
                  idx = imm + ir + (icx-1)*db(ib)

                  re_v = re_mat%elem(idx)
                  im_v = im_mat%elem(idx)
                  amp  = sqrt(re_v**2 + im_v**2)
                  if (amp < amplitude_cutoff) cycle

                  iqp_row = isstart(ib) - 1 + ir
                  iqp_col = isstart(ic) - 1 + icx

                  write(lu,'(1x,i7,i5,i5,f11.4,4i3,2x,a1,3x,"|",i8,i5,i5,f11.4,4i3,2x,a1,3x,"|",1x,3es15.6)') &
                     iqp_row, ib, ir, qp_energy(iqp_row,.true.), &
                        nr(iqp_row), nz(iqp_row), nl(iqp_row), ns(iqp_row), &
                        merge('+','-', npar(iqp_row)==1), &
                     iqp_col, ic, icx, qp_energy(iqp_col,.false.), &
                        nr(iqp_col), nz(iqp_col), nl(iqp_col), ns(iqp_col), &
                        merge('+','-', npar(iqp_col)==1), &
                     re_v, im_v, amp

               end do
            end do
         end do

      end subroutine write_one_amplitude_block

      !> Energy of the row (is_row=.true.) or column (is_row=.false.) quasiparticle,
      !> which is a proton for beta- and a neutron for beta+ (row), and vice-versa (col).
      real(dp) function qp_energy(iqp, is_row)
         implicit none
         integer, intent(in) :: iqp
         logical, intent(in) :: is_row
         logical :: want_proton

         want_proton = (is_row .eqv. bminus)
         if (want_proton) then
            qp_energy = Ep(iqp)
         else
            qp_energy = En(iqp)
         end if
      end function qp_energy

   end subroutine write_amplitudes_table

   !---------------------------------------------------------------------------
   ! Subroutine for solving the FAM equations for one value of energy. Before
   ! calling this the first time, the module data structures must be initialized
   ! with init_pnfam_solver.
   !---------------------------------------------------------------------------
   subroutine ifam
      use pnfam_constants, only : iu, pi, get_timer
      use pnfam_broyden
      implicit none

      ! Finite temperature
      real(dp)    :: ft_real_factor
      complex(dp) :: ft_cmplx_factor

      integer :: iter, ixterm
      real(dp) :: re_s, im_s
      real(dp) :: time1, time2, timei, timei0
      time1=0; time2=0; timei=0; timei0=0

      ! Iteration header
      write(st,'(1x,72("-"))'); call writelog(st)
      write(st,'(1x,5a)') '    i','               si',&
         '               Re(S)','               Im(S)','      Time'; call writelog(st)
      write(st,'(1x,72("-"))'); call writelog(st)

      ! Main iteration loop
      do iter=0,max_iter

         ! Record the iteration (or show where we are starting from)
         write(st,'(1x, i4,a1, 2x, f15.10, 2x, f18.10, 2x, f18.10, 2x, f8.3)') &
            iter, qrpa_bbroyden, si, re_str(1), im_str(1), timei; call writelog(st)
         if (allocated(g)) then
            do ixterm=1, size(g)
               write(st,'(1x,2x,a10,2x,"Re(X)=",f16.8,2x,"Im(X)=",f16.8)') &
                   trim(g(ixterm)%label), re_str(1+ixterm), im_str(1+ixterm); call writelog(st)
            end do
         end if

         ! Exit condition
         if (si < convergence_epsilon) then
            iter_conv = iter; exit
         end if
         if (iter==max_iter) then
            iter_conv = -1; exit
         end if

         ! Start the iteration
         time1 = get_timer()
         if (iter==0) timei0 = time1

         if (abs(quench_residual_int) < 1e-10) then
            call set_val_bbm(dHqp_im,0d0)
            call set_val_bbm(dHqp_re,0d0)
         else

            ! Compute the perturbed densities in SP space
            if (bminus) then
               call triprod_bbm('n',Wp,'n',dRqp_re,'t',Wn,dRsp_re)
               call triprod_bbm('n',Wp,'n',dRqp_im,'t',Wn,dRsp_im)
            else
               call triprod_bbm('n',Wn,'n',dRqp_re,'t',Wp,dRsp_re)
               call triprod_bbm('n',Wn,'n',dRqp_im,'t',Wp,dRsp_im)
            end if

            call calc_hamiltonian(&
                                  dRsp_re%m11, dRsp_im%m11, dRsp_re%m12, dRsp_im%m12,&
                                  dRsp_re%m22, dRsp_im%m22, dRsp_re%m21, dRsp_im%m21,&
                                  dHsp_re%m11, dHsp_im%m11, dHsp_re%m12, dHsp_im%m12,&
                                  dHsp_re%m22, dHsp_im%m22, dHsp_re%m21, dHsp_im%m21 &
                                  )

            ! Transform the Hamiltonian into quasiparticle basis
            if (bminus) then
               call triprod_bbm('t',Wp,'n',dHsp_re,'n',Wn,dHqp_re)
               call triprod_bbm('t',Wp,'n',dHsp_im,'n',Wn,dHqp_im)
            else
               call triprod_bbm('t',Wn,'n',dHsp_re,'n',Wp,dHqp_re)
               call triprod_bbm('t',Wn,'n',dHsp_im,'n',Wp,dHqp_im)
            end if

            ! Quench the residual interaction
            call scalar_mult_bbm(dHqp_re,quench_residual_int)
            call scalar_mult_bbm(dHqp_im,quench_residual_int)

         end if ! End computation of the residual interaction

         ! Add the external field to perturbed hamiltonian (assuming real Fqp)
         call add_bbm(Fqp, dHqp_re)

         ! Solve the new X/Y, P/Q
         call complex_emult_bbm(Greens_re, Greens_im, dHqp_re, dHqp_im,&
                                         dRqp_re, dRqp_im)
         if (use_diagonal_blocks) then
            call emult_bbm(T, dRqp_re, dRqp_re)
            call emult_bbm(T, dRqp_im, dRqp_im)
         end if

         ! Use Broyden mixing to stabilize convergence
         ! T > 0 / odd-A have twice the number of vectors to handle P and Q
         call qrpa_broyden(iter, dRqp_re, dRqp_im)

         ! Strength function
         call contract_bbm(Fqp,dRqp_re,re_s)
         call contract_bbm(Fqp,dRqp_im,im_s)
         re_str(1) = -re_s/pi
         im_str(1) = -im_s/pi

         ! Cross-terms
         if (allocated(g)) then
            do ixterm=1, size(g)
               call contract_bbm(Gqp(ixterm),dRqp_re,re_s)
               call contract_bbm(Gqp(ixterm),dRqp_im,im_s)
               re_str(1+ixterm) = -re_s/pi
               im_str(1+ixterm) = -im_s/pi
            end do
         end if

         ! Record the iteration time
         time2 = get_timer()
         timei = (time2-time1)

      end do

      write(st,'(1x,72("-"))'); call writelog(st)
      if (iter_conv >= 0 ) then
         write(st,'("  *  FAM iteration converged after   ",i4," steps.  si = ",ES12.5)') iter,si; call writelog(st)
      else
         write(st,'("  *  FAM iteration interrupted after ",i4," steps.  si = ",ES12.5)') iter,si; call writelog(st)
      end if
      write(st,'(1x,72("-"))'); call writelog(st)
      write(st,'(1x,a,es12.3,a)') 'Total CPU time = ',(time2-timei0)/60.0_dp,' minutes' ; call writelog(st)
      call writelog("")

   end subroutine ifam


   !---------------------------------------------------------------------------
   ! GMRES subroutine for solving the FAM equations for one value of energy. 
   ! Before calling this the first time, the module data structures must be 
   ! initialized with init_pnfam_solver.
   !
   ! pnFAM equations for the linear response can be desribed as a linear eq
   !                     HZ = b
   ! where Z contains X and Y, and b is -Fqp.
   ! This GMRES solves a preconditioned version, 
   !                    PHZ = Pb
   ! where the preconditioner P is inv(diag((Ep+En-omega);(Ep+En+omega))),
   ! with an option to restart, i.e., reset the Krylov subspace after few its.
   !---------------------------------------------------------------------------
   subroutine gmres
      use pnfam_constants, only : iu, pi, get_timer
      use pnfam_broyden
      implicit none

      complex(dp), parameter :: cone=(1.0d0, 0.d0), czero=(0.d0, 0.d0) ! declaration of complex one and complex zero
      real(dp), external :: dznrm2 ! LAPACK function for euclidean norm

      complex(dp), allocatable :: basis(:, :), Hmat(:, :) ! Krylov subspace & Hessenberg matrix

      complex(dp), allocatable :: sn(:), cn(:) ! arrays holding Givens rotations
      complex(dp), allocatable :: betavec(:), errorvec(:) ! arrays holding beta & residual norm values
      complex(dp), allocatable :: mres_x(:), Z0(:) ! gmres solution vector & initial vector

      real(dp) :: betanrm0, betanrm, error ! initial residual norm & beta & residual norm
      
      integer :: restart, ncycle ! number of restarts & cycle counter
      integer :: counter ! counter for inner iteration (gmres)

      ! Finite temperature
      real(dp)    :: ft_real_factor
      complex(dp) :: ft_cmplx_factor

      integer :: iter, ixterm, nxy, i, j
      real(dp) :: re_s, im_s
      real(dp) :: time1, time2, timei, timei0
      time1=0; time2=0; timei=0; timei0=0      

      ! Iteration header
      write(st,'(1x,72("-"))'); call writelog(st)
      write(st,'(1x,3a)') '    i', '              Error', '               Time'; call writelog(st)
      write(st,'(1x,72("-"))'); call writelog(st)

      ! Allocate memory
      nxy = size(f%mat%elem)
      allocate(basis(2*nxy, max_iter+1))
      allocate(Hmat(max_iter+1, max_iter))
      allocate(mres_x(2*nxy))
      allocate(Z0(2*nxy))

      allocate(betavec(max_iter+1))
      allocate(cn(max_iter))
      allocate(sn(max_iter))
      allocate(errorvec(max_iter+1))

      !! Set F20=F02
      !do i = 1, nxy
      !   Fqp%m21%elem(i) = Fqp%m12%elem(i)
      !end do

      ! Obtain preconditioned rhs
      ! Set Init_ as Fqp
      Init_re = Fqp
      Init_im = Fqp
      Out_re = Fqp
      Out_im = Fqp
      do i = 1, nxy
         Init_im%m12%elem(i) = 0.0
         Init_im%m21%elem(i) = 0.0
      end do
      ! Apply P * Init_, which is set to Out_
      call complex_emult_bbm(Greens_re, Greens_im, Init_re, Init_im,&
                                         Out_re, Out_im)

      ! Set the number of restarts (NEED THIS ENCODED INTO USER INPUT FILE e.g. GT-k0.in)
      restart = 3
      do ncycle=1,restart

         write(st,'(" ****** At cycle number ", i4, " *****.")') ncycle; call writelog(st)

         ! Start the cycle
         time1 = get_timer()
         if (ncycle==1) timei0 = time1

         basis = czero
         Hmat = czero
         cn = czero
         sn = czero
         betavec = czero
         errorvec = 0.d0

         if (ncycle == 1) then
            ! For fist ever initialization, Z0 is zero and r0 is just b=-Fqp
            Z0 = czero
            do i = 1, nxy
               basis(i, 1) = dcmplx(-Out_re%m12%elem(i), -Out_im%m12%elem(i))
               basis(nxy+i, 1) = dcmplx(-Out_re%m21%elem(i), -Out_im%m21%elem(i))
            end do
            betanrm0 = dznrm2(2*nxy, basis(1, 1), 1) ! use to compute relative residual norm
         else
            ! For restart, need b-H*Z0
            ! Keep Z0 for final solution at the end of cycle
            Z0 = mres_x
            ! Do matvec H*Z0 (dRqp holds mres_x in bigblockmatrix form)
            if (abs(quench_residual_int) < 1e-10) then
               call set_val_bbm(dHqp_im,0d0)
               call set_val_bbm(dHqp_re,0d0)
            else

               ! Compute the perturbed densities in SP space
               if (bminus) then
                  call triprod_bbm('n',Wp,'n',dRqp_re,'t',Wn,dRsp_re)
                  call triprod_bbm('n',Wp,'n',dRqp_im,'t',Wn,dRsp_im)
               else
                  call triprod_bbm('n',Wn,'n',dRqp_re,'t',Wp,dRsp_re)
                  call triprod_bbm('n',Wn,'n',dRqp_im,'t',Wp,dRsp_im)
               end if

               call calc_hamiltonian(& 
                                     dRsp_re%m11, dRsp_im%m11, dRsp_re%m12, dRsp_im%m12,&
                                     dRsp_re%m22, dRsp_im%m22, dRsp_re%m21, dRsp_im%m21,&
                                     dHsp_re%m11, dHsp_im%m11, dHsp_re%m12, dHsp_im%m12,&
                                     dHsp_re%m22, dHsp_im%m22, dHsp_re%m21, dHsp_im%m21 &
                                    )

               ! Transform the Hamiltonian into quasiparticle basis
               if (bminus) then
                  call triprod_bbm('t',Wp,'n',dHsp_re,'n',Wn,dHqp_re)
                  call triprod_bbm('t',Wp,'n',dHsp_im,'n',Wn,dHqp_im)
               else
                  call triprod_bbm('t',Wn,'n',dHsp_re,'n',Wp,dHqp_re)
                  call triprod_bbm('t',Wn,'n',dHsp_im,'n',Wp,dHqp_im)
               end if

               ! Quench the residual interaction
               call scalar_mult_bbm(dHqp_re,quench_residual_int)
               call scalar_mult_bbm(dHqp_im,quench_residual_int)

            end if ! End computation of the residual interaction

            ! Apply preconditioner to dH : P*dH
            temp_re = dRqp_re ! arbitrary
            temp_im = dRqp_im ! arbitrary
            call complex_emult_bbm(Greens_re, Greens_im, dHqp_re, dHqp_im,&
                                            temp_re, temp_im)
            ! Combine Z0 + P*dH
            call add_bbm(temp_re, dRqp_re)
            call add_bbm(temp_im, dRqp_im)
            ! End of matvec

            ! Get b - H*Z0 = -Fqp - H*Z0 = -(Fqp + H*Z0)
            call add_bbm(Out_re, dRqp_re)
            call add_bbm(Out_im, dRqp_im)
            do i = 1, nxy
               basis(i, 1) = dcmplx(-dRqp_re%m12%elem(i), -dRqp_im%m12%elem(i))
               basis(nxy+i, 1) = dcmplx(-dRqp_re%m21%elem(i), -dRqp_im%m21%elem(i))
            end do
         end if

         ! Normalize the first basis vector
         betanrm = dznrm2(2*nxy, basis(1, 1), 1)
         basis(:, 1) = basis(:, 1) / betanrm
         ! Initialize residual
         betavec(1) = betanrm
         errorvec(1) = betanrm

         ! Record the iteration time
         time2 = get_timer()
         timei = (time2-time1)

         ! Record the initial result
         write(st,'(1x, i4, 2x, f20.10, 2x, f15.3)') 0, betanrm, timei; call writelog(st)
         ! End of initialization step

         ! Main gmres loop
         do iter=1,max_iter

            ! Start the iteration
            time1 = get_timer()

            ! Set Z (X and Y)
            do i=1,nxy
               dRqp_re%m12%elem(i) = dble(basis(i, iter))
               dRqp_re%m21%elem(i) = dble(basis(nxy+i, iter))
               dRqp_im%m12%elem(i) = dimag(basis(i, iter))
               dRqp_im%m21%elem(i) = dimag(basis(nxy+i, iter))
            end do 

            ! Do matvec H*Z

            if (abs(quench_residual_int) < 1e-10) then
               call set_val_bbm(dHqp_im,0d0)
               call set_val_bbm(dHqp_re,0d0)
            else

               ! Compute the perturbed densities in SP space
               if (bminus) then
                  call triprod_bbm('n',Wp,'n',dRqp_re,'t',Wn,dRsp_re)
                  call triprod_bbm('n',Wp,'n',dRqp_im,'t',Wn,dRsp_im)
               else
                  call triprod_bbm('n',Wn,'n',dRqp_re,'t',Wp,dRsp_re)
                  call triprod_bbm('n',Wn,'n',dRqp_im,'t',Wp,dRsp_im)
               end if

               call calc_hamiltonian(&
                                     dRsp_re%m11, dRsp_im%m11, dRsp_re%m12, dRsp_im%m12,&
                                     dRsp_re%m22, dRsp_im%m22, dRsp_re%m21, dRsp_im%m21,&
                                     dHsp_re%m11, dHsp_im%m11, dHsp_re%m12, dHsp_im%m12,&
                                     dHsp_re%m22, dHsp_im%m22, dHsp_re%m21, dHsp_im%m21 &
                                     )

               ! Transform the Hamiltonian into quasiparticle basis
               if (bminus) then
                  call triprod_bbm('t',Wp,'n',dHsp_re,'n',Wn,dHqp_re)
                  call triprod_bbm('t',Wp,'n',dHsp_im,'n',Wn,dHqp_im)
               else
                  call triprod_bbm('t',Wn,'n',dHsp_re,'n',Wp,dHqp_re)
                  call triprod_bbm('t',Wn,'n',dHsp_im,'n',Wp,dHqp_im)
               end if

               ! Quench the residual interaction
               call scalar_mult_bbm(dHqp_re,quench_residual_int)
               call scalar_mult_bbm(dHqp_im,quench_residual_int)

            end if ! End computation of the residual interaction
        
            ! Apply preconditioner to dH : P*dH 
            temp_re = dRqp_re ! arbitrary
            temp_im = dRqp_im ! arbitrary
            call complex_emult_bbm(Greens_re, Greens_im, dHqp_re, dHqp_im,&
                                            temp_re, temp_im)
            ! Combine Z + P*dH
            call add_bbm(temp_re, dRqp_re)
            call add_bbm(temp_im, dRqp_im)
            ! End of matvec

            ! Store into Krylov subspace
            do i = 1,nxy
               basis(i, iter+1) = dcmplx(dRqp_re%m12%elem(i), dRqp_im%m12%elem(i))
               basis(nxy+i, iter+1) = dcmplx(dRqp_re%m21%elem(i), dRqp_im%m21%elem(i))
            end do
            ! Orthogonalize w.r.t. previous Krylov basis vectors
            call zgemv('C', 2*nxy, iter, cone, basis, 2*nxy, basis(:, iter+1), 1, czero, Hmat(1:iter, iter), 1)
            call zgemv('N', 2*nxy, iter, -cone, basis, 2*nxy, Hmat(1:iter, iter), 1, cone, basis(:, iter+1), 1)

            ! Normalize current basis vector
            betanrm = dznrm2(2*nxy, basis(1, iter+1), 1)
            basis(:, iter+1) = basis(:, iter+1) / betanrm
            ! Set subdiag of Hessenberg
            Hmat(iter+1, iter) = betanrm

            ! Apply Givens rotations
            do i = 1,iter-1
               call zrot(1, Hmat(i, iter), 1, Hmat(i+1, iter), 1, cn(i), sn(i))   
            end do
            call zrotg(Hmat(iter, iter), Hmat(iter+1, iter), cn(iter), sn(iter))
            Hmat(iter+1, iter) = czero
            ! Residual vector
            call zrot(1, betavec(iter), 1, betavec(iter+1), 1, cn(iter), sn(iter))
            error = abs(betavec(iter+1)) / betanrm0
            errorvec(iter+1) = error

            ! Record the iteration time
            time2 = get_timer()
            timei = (time2-time1)

            ! Record the iteration (or show where we are starting from)
            write(st,'(1x, i4, 2x, f20.10, 2x, f15.3)') iter, error, timei; call writelog(st)

            counter = iter ! Need separate counter for correct number of iterations

            ! Exit condition
            if (error < convergence_epsilon) then
               iter_conv = iter; exit
            end if

         end do ! end inner iteration

         ! Solve triangular problem and get GMRES result mres_x
         call ztrsv('U', 'N', 'N', counter, Hmat, max_iter+1, betavec, 1)
         call zgemv('N', 2*nxy, counter, cone, basis, 2*nxy, betavec, 1, czero, mres_x, 1)
         ! Combine Z0 + mres_x
         mres_x = Z0 + mres_x

         do i=1,nxy
            dRqp_re%m12%elem(i) = dble(mres_x(i))
            dRqp_re%m21%elem(i) = dble(mres_x(nxy+i))
            dRqp_im%m12%elem(i) = imag(mres_x(i))
            dRqp_im%m21%elem(i) = imag(mres_x(nxy+i))
         end do

         ! Strength function
         call contract_bbm(Fqp,dRqp_re,re_s)
         call contract_bbm(Fqp,dRqp_im,im_s)
         re_str(1) = -re_s/pi
         im_str(1) = -im_s/pi

         ! Cross-terms
         if (allocated(g)) then
            do ixterm=1, size(g)
               call contract_bbm(Gqp(ixterm),dRqp_re,re_s)
               call contract_bbm(Gqp(ixterm),dRqp_im,im_s)
               re_str(1+ixterm) = -re_s/pi
               im_str(1+ixterm) = -im_s/pi
            end do
         end if

         time2 = get_timer()

         ! Exit restart if converged
         if (iter_conv > 0) then 
            exit
         end if
         if (ncycle==restart) then
            iter_conv = -1; exit
         end if

      end do !end cycle

      write(st,'(1x,72("-"))'); call writelog(st)
      if (iter_conv >= 0) then
         write(st,'("  *  GMRES iteration converged after   ", i4, " steps.")') (ncycle-1)*max_iter+iter; call writelog(st)
      else
         write(st,'("  *  GMRES iteration did not converge within step   ", i4, " steps.")') (ncycle-1)*max_iter+iter; call writelog(st)
      end if
      write(st,'(1x,72("-"))'); call writelog(st)
      write(st,'(1x,a,es12.3,a)') 'Total CPU time = ',(time2-timei0)/60.0_dp,' minutes' ; call writelog(st)
      call writelog("")

      ! Deallocate
      deallocate(basis)
      deallocate(Hmat)
      deallocate(mres_x)
      deallocate(Z0)

      deallocate(betavec)
      deallocate(cn)
      deallocate(sn)
      deallocate(errorvec)

   end subroutine gmres



   !---------------------------------------------------------------------------
   ! LANCZOS subroutine for computing the response function by solving the
   ! generalized eigenvalue problem.
   ! Before calling this the first time, the module data structures must be
   ! initialized with init_pnfam_solver.
   !
   !                    | A  B| = |X Y||D  ||X Y|^-1
   !                    |-B -A|   |Y X|| -D||Y X|
   !
   ! where X^TX-Y^TY=I and X^Y-Y^X=0, 
   ! and eigenvalues in D are in non-decreasing order.
   !
   ! Solves the reduced eigenvalue problem
   !                         MK = (X-Y)D^2(X-Y)^-1
   ! and
   !                         KM = (X+Y)D^2(X+Y)^-1
   ! where M:=A+B and K:=A-B.
   !
   ! We have access to the matvec with the matrix [A B; B A].
   ! To obtain the matvecs with the matrix M and K, we set the vector as [Z;0]
   ! so that            |A B||Z| = |AZ|
   !                    |B A||0|   |BZ|
   ! then do AZ+BZ for MZ or AZ-BZ for KZ.
   !
   ! Obtains the Ritz values and the first element squared 
   ! of the corresponding eigenvector of MK and KM eigenvalue problems.
   !---------------------------------------------------------------------------
   subroutine lanczos
      use pnfam_constants, only : iu, pi, get_timer
      use pnfam_broyden
      implicit none

      ! Parameters to keep from original code
      integer :: iter, ixterm, nxy, i, j
      real(dp) :: re_s, im_s
      real(dp) :: time1, time2, timei, timei0

      integer :: lu1, lu2, lu3, lu4                                     ! Write integer files for Theta_MK, Theta_KM, Tau2_MK, Tau2_KM
      integer :: ncycle                                                 ! Current number of Lanczos cycle

      complex(dp), parameter :: cone=(1.0d0, 0.d0), czero=(0.d0, 0.d0)  ! declaration of complex one and complex zero
      real(dp), external :: dznrm2                                      ! LAPACK function for euclidean norm

      real(dp), allocatable :: re_rand(:)                               ! Random real parts
      integer :: Nvec                                                   ! Number of Lanczos iterations

      ! Krylov subspace parameters
      complex(dp), allocatable :: basis_MK(:, :), basis_KM(:, :)        ! Krylov subspaces
      complex(dp), allocatable :: Hmat_MK(:, :), Hmat_KM(:, :)          ! Hessenberg matrices
      complex(dp), allocatable :: copyHmat_MK(:, :), copyHmat_KM(:, :)  ! Copy of Hessenberg matrices; for reorthogonalization and zgeev
      real(dp) :: betanrm                                               ! For norm of the basis vectors
      complex(dp), allocatable :: ZZ(:), AZ(:), BZ(:)                   ! For storing intermediate matvecs

      ! zgeev parameters
      complex(dp), allocatable :: leigvecs_MK(:, :), leigvecs_KM(:, :)  ! Left eigenvectors
      complex(dp), allocatable :: reigvecs_MK(:, :), reigvecs_KM(:, :)  ! Right eigenvectors
      complex(dp), allocatable :: workvec_MK(:), workvec_KM(:)          ! Work vectors
      real(dp), allocatable :: rwork_MK(:), rwork_KM(:)                 ! rwork vectors
      integer :: info_MK, info_KM                                       ! Info values (=0 indicates convergence)

      ! Eigen parameters
      complex(dp), allocatable :: Theta_MK(:, :), Theta_KM(:, :)        ! Ritz values
      integer, allocatable :: idx(:)                                    ! For sorting eigenvalues
      integer :: temp_idx, k                                            ! For sorting eigenvalues
      complex(dp) :: temp                                               ! For sorting eigenvalues
      complex(dp), allocatable :: Temp_reigvecs(:, :)                   ! For sorting eigenvalues
      complex(dp), allocatable :: Tau2_MK(:, :), Tau2_KM(:, :)          ! Squared first components of the eigenvectors of Hessenberg matrices
      complex(dp), allocatable :: ritzvecs_MK(:, :), ritzvecs_KM(:, :)  ! Ritz vectors


      ! DMR
      complex(dp), allocatable :: U(:, :), V(:, :)
      complex(dp), allocatable :: F20(:), F02(:)
      complex(dp), allocatable :: VTF20(:, :), UTF02(:, :)



      complex(dp), allocatable :: copyHmat(:, :)
      real(dp), allocatable :: AA(:, :)
      real(dp), allocatable :: BB(:, :)


      complex(dp), allocatable :: Z(:, :)
      complex(dp), allocatable :: KZ(:, :)
      complex(dp), allocatable :: RF0(:, :), RF1(:, :)


      integer :: lu5, lu6, lu7, lu8

      Nvec = 1                                                          ! Shoule be kept as 1
      nxy = size(f%mat%elem)                                            ! Size of A and B

      ! Iteration header
      write(st,'(1x,72("-"))'); call writelog(st)
      write(st,'(1x,3a)') '    i', '              Error', '               Time'; call writelog(st)
      write(st,'(1x,72("-"))'); call writelog(st)

      ! Allocate memory
      allocate(re_rand(nxy))
      allocate(basis_MK(nxy, max_iter+1))
      allocate(basis_KM(nxy, max_iter+1))
      allocate(Hmat_MK(max_iter+1, max_iter))
      allocate(Hmat_KM(max_iter+1, max_iter))
      allocate(copyHmat_MK(max_iter+1, max_iter))
      allocate(copyHmat_KM(max_iter+1, max_iter))
      allocate(ZZ(nxy))
      allocate(AZ(nxy))
      allocate(BZ(nxy))
      allocate(leigvecs_MK(max_iter, max_iter))
      allocate(leigvecs_KM(max_iter, max_iter))
      allocate(reigvecs_MK(max_iter, max_iter))
      allocate(reigvecs_KM(max_iter, max_iter))
      allocate(workvec_MK(2*max_iter+1))
      allocate(workvec_KM(2*max_iter+1))
      allocate(rwork_MK(2*max_iter))
      allocate(rwork_KM(2*max_iter))
      allocate(Theta_MK(max_iter, Nvec))
      allocate(Theta_KM(max_iter, Nvec))
      allocate(idx(max_iter))
      allocate(Temp_reigvecs(max_iter, max_iter))
      allocate(Tau2_MK(max_iter, Nvec))
      allocate(Tau2_KM(max_iter, Nvec))
      allocate(ritzvecs_MK(nxy, max_iter))
      allocate(ritzvecs_KM(nxy, max_iter))


      ! DMR
      allocate(U(nxy, max_iter))
      allocate(V(nxy, max_iter))
      allocate(F20(nxy))
      allocate(F02(nxy))
      allocate(VTF20(max_iter, Nvec))
      allocate(UTF02(max_iter, Nvec))


      allocate(copyHmat(max_iter, max_iter))
      allocate(AA(nxy, nxy))
      allocate(BB(nxy, nxy))




      allocate(Z(nxy, max_iter))
      allocate(KZ(nxy, max_iter))
      allocate(RF0(max_iter, Nvec))
      allocate(RF1(max_iter, Nvec))




      print *, 'nxy is '
      print *, nxy

      time1=0; time2=0; timei=0; timei0=0
      Theta_MK = czero
      Theta_KM = czero
      Tau2_MK = czero
      Tau2_KM = czero

      !! Set F20=F02
      !do i = 1, nxy
      !   Fqp%m21%elem(i) = Fqp%m12%elem(i)
      !end do

      do ncycle=1,Nvec

         write(st,'(" ****** At cycle number ", i4, " *****.")') ncycle; call writelog(st)

         ! Start the cycle
         time1 = get_timer()
         if (ncycle==1) timei0 = time1

         re_rand = 0.d0
         basis_MK = czero
         basis_KM = czero
         Hmat_MK = czero
         Hmat_KM = czero
         copyHmat_MK = czero
         copyHmat_KM = czero
         leigvecs_MK = czero
         leigvecs_KM = czero
         reigvecs_MK = czero
         reigvecs_KM = czero
         workvec_MK = czero
         workvec_KM = czero
         rwork_MK = 0.d0
         rwork_KM = 0.d0
         Temp_reigvecs = czero
         ritzvecs_MK = czero
         ritzvecs_KM = czero

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!! Initialize !!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         do i = 1, nxy
!            basis_MK(i, 1) = Fqp%m12%elem(i)
            basis_MK(i, 1) = Fqp%m12%elem(i) + Fqp%m21%elem(i)
!            basis_KM(i, 1) = Fqp%m12%elem(i) - Fqp%m21%elem(i)
!            
!            ! DMR
!            F20(i) = Fqp%m12%elem(i)
!            F02(i) = Fqp%m21%elem(i)
         end do
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!! End of Initialization !!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!! Normalize the first basis vector !!!!
         !    For MK problem, normalize w.r.t. K  !
         !    For KM problem, normalize w.r.t. M  !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! 
         ! a) MK problem
         !
         ! Set [Z;0]
         do i=1,nxy
            dRqp_re%m12%elem(i) = dble(basis_MK(i, 1))
            dRqp_re%m21%elem(i) = 0.d0
            dRqp_im%m12%elem(i) = dimag(basis_MK(i, 1))
            dRqp_im%m21%elem(i) = 0.d0
         end do
         ! Do matvec and get [AZ;BZ]
         call matvec(dRqp_re, dRqp_im, temp_re, temp_im)
         do i = 1, nxy
            AZ(i) = dcmplx(temp_re%m12%elem(i), temp_im%m12%elem(i))
            BZ(i) = dcmplx(temp_re%m21%elem(i), temp_im%m21%elem(i))
         end do
         ! KZ = AZ - BZ
         ZZ = AZ - BZ
         betanrm = sqrt(real(dot_product(conjg(basis_MK(:, 1)), ZZ)))
         basis_MK(:, 1) = basis_MK(:, 1) / betanrm
!         !
!         ! b) KM problem
!         ! 
!         ! Set [Z;0]
!         do i=1,nxy
!            dRqp_re%m12%elem(i) = dble(basis_KM(i, 1))
!            dRqp_re%m21%elem(i) = 0.d0
!            dRqp_im%m12%elem(i) = dimag(basis_KM(i, 1))
!            dRqp_im%m21%elem(i) = 0.d0
!         end do
!         ! Do matvec and get [AZ; BZ]
!         call matvec(dRqp_re, dRqp_im, temp_re, temp_im)
!         do i = 1, nxy
!            AZ(i) = dcmplx(temp_re%m12%elem(i), temp_im%m12%elem(i))
!            BZ(i) = dcmplx(temp_re%m21%elem(i), temp_im%m21%elem(i))
!         end do
!         ! MZ = AZ + BZ
!         ZZ = AZ + BZ
!         betanrm = sqrt(real(dot_product(conjg(basis_KM(:, 1)), ZZ)))
!         basis_KM(:, 1) = basis_KM(:, 1) / betanrm
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! End of first basis vector normalization!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! Record the iteration time
         time2 = get_timer()
         timei = (time2-time1)

         ! Record the initial result
         write(st,'(1x, i4, 2x, f20.10, 2x, f15.3)') 0, betanrm, timei; call writelog(st)
         ! End of initialization step

         ! Main lanczos loop
         do iter=1,max_iter

            ! Start the iteration
            time1 = get_timer()

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!! Generate Krylov subspace for MK !!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !
            ! a) Matvec with MK
            !
            ! Set [Z;0]
            do i=1,nxy
               dRqp_re%m12%elem(i) = dble(basis_MK(i, iter))
               dRqp_re%m21%elem(i) = 0.d0
               dRqp_im%m12%elem(i) = dimag(basis_MK(i, iter))
               dRqp_im%m21%elem(i) = 0.d0
            end do
            ! Do matvec and get [AZ;BZ]
            call matvec(dRqp_re, dRqp_im, temp_re, temp_im)
            do i = 1, nxy
               AZ(i) = dcmplx(temp_re%m12%elem(i), temp_im%m12%elem(i))
               BZ(i) = dcmplx(temp_re%m21%elem(i), temp_im%m21%elem(i))
            end do
            ! KZ = AZ-BZ
            ZZ = AZ - BZ
            ! Set [KZ;0]
            do i=1,nxy
               dRqp_re%m12%elem(i) = dble(ZZ(i))
               dRqp_re%m21%elem(i) = 0.d0
               dRqp_im%m12%elem(i) = dimag(ZZ(i))
               dRqp_im%m21%elem(i) = 0.d0
            end do
            ! Do matvec and get [AKZ;BKZ]
            call matvec(dRqp_re, dRqp_im, temp_re, temp_im)
            do i = 1, nxy
               AZ(i) = dcmplx(temp_re%m12%elem(i), temp_im%m12%elem(i))
               BZ(i) = dcmplx(temp_re%m21%elem(i), temp_im%m21%elem(i))
            end do
            ! MKZ = AKZ+BKZ
            ZZ = AZ + BZ
            ! Store into basis
            do i = 1,nxy
               basis_MK(i, iter+1) = ZZ(i)
            end do
            !
            ! b) K-inner product of MKZ
            !
            ! Set [MKZ;0]
            do i=1,nxy
               dRqp_re%m12%elem(i) = dble(ZZ(i))
               dRqp_re%m21%elem(i) = 0.d0
               dRqp_im%m12%elem(i) = dimag(ZZ(i))
               dRqp_im%m21%elem(i) = 0.d0
            end do
            ! Do matvec and get [AMKZ;BMKZ]
            call matvec(dRqp_re, dRqp_im, temp_re, temp_im)
            do i = 1, nxy
               AZ(i) = dcmplx(temp_re%m12%elem(i), temp_im%m12%elem(i))
               BZ(i) = dcmplx(temp_re%m21%elem(i), temp_im%m21%elem(i))
            end do
            ! KMKZ = AMKZ-BMKZ
            ZZ = AZ - BZ
            !
            ! c) Orthogonalize w.r.t. previous Krylov basis vectors
            !
            call zgemv('C', nxy, iter, cone, basis_MK, nxy, ZZ, 1, czero, Hmat_MK(1:iter, iter), 1)
            call zgemv('N', nxy, iter, -cone, basis_MK, nxy, Hmat_MK(1:iter, iter), 1, cone, basis_MK(:, iter+1), 1)
            !
            ! d) Reorthogonalize
            !
            !    d)a) K-inner product of f:=basis_MK(:, iter+1)
            !
            ! Set [f;0]
            do i=1,nxy
               dRqp_re%m12%elem(i) = dble(basis_MK(i, iter+1))
               dRqp_re%m21%elem(i) = 0.d0
               dRqp_im%m12%elem(i) = dimag(basis_MK(i, iter+1))
               dRqp_im%m21%elem(i) = 0.d0
            end do
            ! Do matvec and get [Af;Bf]
            call matvec(dRqp_re, dRqp_im, temp_re, temp_im)
            do i = 1, nxy
               AZ(i) = dcmplx(temp_re%m12%elem(i), temp_im%m12%elem(i))
               BZ(i) = dcmplx(temp_re%m21%elem(i), temp_im%m21%elem(i))
            end do
            ! Kf = Af-Bf
            ZZ = AZ - BZ
            !
            !    d)b) Orthogonalize w.r.t. previous Krylov basis vectors
            !
            call zgemv('C', nxy, iter, cone, basis_MK, nxy, ZZ, 1, czero, copyHmat_MK(1:iter, iter), 1)
            call zgemv('N', nxy, iter, -cone, basis_MK, nxy, copyHmat_MK(1:iter, iter), 1, cone, basis_MK(:, iter+1), 1)
            Hmat_MK(1:iter, iter) = Hmat_MK(1:iter, iter) + copyHmat_MK(1:iter, iter)
            
            


            if (iter > 2) then
               do i = 1, iter-2
                  Hmat_MK(i, iter) = czero
               enddo
            end if
            
            
            
            !
            ! e) Normalize current basis vector w.r.t. K
            !
            ! Set [Z;0]
            do i=1,nxy
               dRqp_re%m12%elem(i) = dble(basis_MK(i, iter+1))
               dRqp_re%m21%elem(i) = 0.d0
               dRqp_im%m12%elem(i) = dimag(basis_MK(i, iter+1))
               dRqp_im%m21%elem(i) = 0.d0
            end do
            ! Do matvec and get [AZ;BZ]
            call matvec(dRqp_re, dRqp_im, temp_re, temp_im)
            do i = 1, nxy
               AZ(i) = dcmplx(temp_re%m12%elem(i), temp_im%m12%elem(i))
               BZ(i) = dcmplx(temp_re%m21%elem(i), temp_im%m21%elem(i))
            end do
            ! KZ = AZ-BZ
            ZZ = AZ - BZ
            betanrm = sqrt(real(dot_product(conjg(basis_MK(:, iter+1)), ZZ)))
            basis_MK(:, iter+1) = basis_MK(:, iter+1) / betanrm
            Hmat_MK(iter+1, iter) = betanrm
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! End Krylov subspace generation for MK !!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
            
!            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            !!!!! Generate Krylov subspace for KM !!!!
!            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            !
!            ! a) Matvec with KM
!            !
!            ! Set [Z;0]
!            do i=1,nxy
!               dRqp_re%m12%elem(i) = dble(basis_KM(i, iter))
!               dRqp_re%m21%elem(i) = 0.d0
!               dRqp_im%m12%elem(i) = dimag(basis_KM(i, iter))
!               dRqp_im%m21%elem(i) = 0.d0
!            end do
!            ! Do matvec and get [AZ;BZ]
!            call matvec(dRqp_re, dRqp_im, temp_re, temp_im)
!            do i = 1, nxy
!               AZ(i) = dcmplx(temp_re%m12%elem(i), temp_im%m12%elem(i))
!               BZ(i) = dcmplx(temp_re%m21%elem(i), temp_im%m21%elem(i))
!            end do
!            ! MZ = AZ+BZ
!            ZZ = AZ + BZ
!            ! Set [KZ;0]
!            do i=1,nxy
!               dRqp_re%m12%elem(i) = dble(ZZ(i))
!               dRqp_re%m21%elem(i) = 0.d0
!               dRqp_im%m12%elem(i) = dimag(ZZ(i))
!               dRqp_im%m21%elem(i) = 0.d0
!            end do
!            ! Do matvec and get [AKZ;BKZ]
!            call matvec(dRqp_re, dRqp_im, temp_re, temp_im)
!            do i = 1, nxy
!               AZ(i) = dcmplx(temp_re%m12%elem(i), temp_im%m12%elem(i))
!               BZ(i) = dcmplx(temp_re%m21%elem(i), temp_im%m21%elem(i))
!            end do
!            ! KMZ = AMZ-BMZ
!            ZZ = AZ - BZ
!            ! Store into basis
!            do i = 1,nxy
!               basis_KM(i, iter+1) = ZZ(i)
!            end do
!            !
!            ! b) M-inner product of KMZ
!            !
!            ! Set [KMZ;0]
!            do i=1,nxy
!               dRqp_re%m12%elem(i) = dble(ZZ(i))
!               dRqp_re%m21%elem(i) = 0.d0
!               dRqp_im%m12%elem(i) = dimag(ZZ(i))
!               dRqp_im%m21%elem(i) = 0.d0
!            end do
!            ! Do matvec and get [AKMZ;BKMZ]
!            call matvec(dRqp_re, dRqp_im, temp_re, temp_im)
!            do i = 1, nxy
!               AZ(i) = dcmplx(temp_re%m12%elem(i), temp_im%m12%elem(i))
!               BZ(i) = dcmplx(temp_re%m21%elem(i), temp_im%m21%elem(i))
!            end do
!            ! MKMZ = AKMZ+BKMZ
!            ZZ = AZ + BZ
!            !
!            ! c) Orthogonalize w.r.t. previous Krylov basis vectors
!            !
!            call zgemv('C', nxy, iter, cone, basis_KM, nxy, ZZ, 1, czero, Hmat_KM(1:iter, iter), 1)
!            call zgemv('N', nxy, iter, -cone, basis_KM, nxy, Hmat_KM(1:iter, iter), 1, cone, basis_KM(:, iter+1), 1)
!            !
!            ! d) Reorthogonalize
!            !
!            !    d)a) M-inner product of f:=basis_KM(:, iter+1)
!            !
!            ! Set [f;0]
!            do i=1,nxy
!               dRqp_re%m12%elem(i) = dble(basis_KM(i, iter+1))
!               dRqp_re%m21%elem(i) = 0.d0
!               dRqp_im%m12%elem(i) = dimag(basis_KM(i, iter+1))
!               dRqp_im%m21%elem(i) = 0.d0
!            end do
!            ! Do matvec and get [Af;Bf]
!            call matvec(dRqp_re, dRqp_im, temp_re, temp_im)
!            do i = 1, nxy
!               AZ(i) = dcmplx(temp_re%m12%elem(i), temp_im%m12%elem(i))
!               BZ(i) = dcmplx(temp_re%m21%elem(i), temp_im%m21%elem(i))
!            end do
!            ! Mf = Af+Bf
!            ZZ = AZ + BZ
!            !
!            !    d)b) Orthogonalize w.r.t. previous Krylov basis vectors
!            !
!            call zgemv('C', nxy, iter, cone, basis_KM, nxy, ZZ, 1, czero, copyHmat_KM(1:iter, iter), 1)
!            call zgemv('N', nxy, iter, -cone, basis_KM, nxy, copyHmat_KM(1:iter, iter), 1, cone, basis_KM(:, iter+1), 1)
!            Hmat_KM(1:iter, iter) = Hmat_KM(1:iter, iter) + copyHmat_KM(1:iter, iter)
!            !
!            ! e) Normalize current basis vector w.r.t. M
!            !
!            ! Set [Z;0]
!            do i=1,nxy
!               dRqp_re%m12%elem(i) = dble(basis_KM(i, iter+1))
!               dRqp_re%m21%elem(i) = 0.d0

            time2 = get_timer()
            timei = (time2-time1)
            write(st,'(1x, i4, 2x, f20.10, 2x, f15.3)') iter, betanrm, timei; call writelog(st)

         end do ! end inner iteration

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!! Solve eigenvalue problems !!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         ! a) MK eigenvalue problem
         !
!         call zgeev('N', 'V', max_iter, Hmat_MK, max_iter, Theta_MK(:, ncycle), leigvecs_MK, max_iter, reigvecs_MK, max_iter, &
!                    workvec_MK, 2*max_iter+1, rwork_MK, info_MK)
!         copyHmat_MK (1:max_iter, 1:max_iter)= Hmat_MK (1:max_iter, 1:max_iter)! zgeev subroutine overwrites matrix


         copyHmat = Hmat_MK(1:max_iter, 1:max_iter) ! zgeev subroutine overwrites matrix


         copyHmat = 0.5d0 * (copyHmat + transpose(conjg(copyHmat)))




!         print *, "PRINTING Hmat"
!         do i = 1, max_iter
!            print *, copyHmat(i, :)
!         end do





         call zgeev('N', 'V', max_iter, copyHmat, max_iter, Theta_MK(:, ncycle), leigvecs_MK, max_iter, reigvecs_MK, max_iter, &
                    workvec_MK, 2*max_iter+1, rwork_MK, info_MK)

!         copyHmat_MK = Hmat_MK(1:max_iter, 1:max_iter) ! zgeev subroutine overwrites matrix
!         call zgeev('N', 'V', max_iter, copyHmat_MK, max_iter, Theta_MK(:, ncycle), leigvecs_MK, max_iter, reigvecs_MK, max_iter, &
!                    workvec_MK, 2*max_iter+1, rwork_MK, info_MK)
         ! Initialize indices
         do i = 1, max_iter
            idx(i) = i
         end do
         ! Simple bubble sort
         do i = 1, max_iter-1
            k = i
            do j = i+1, max_iter
               if (real(Theta_MK(j, ncycle)) < real(Theta_MK(k, ncycle))) k = j
            end do
            if (k /= i) then
               ! swap values
               temp = Theta_MK(i, ncycle); Theta_MK(i, ncycle) = Theta_MK(k, ncycle); Theta_MK(k, ncycle) = temp
               ! swap indices
               temp_idx = idx(i); idx(i) = idx(k); idx(k) = temp_idx
            end if
         end do
         Temp_reigvecs = reigvecs_MK
         do i = 1, max_iter
            reigvecs_MK(:, i) = Temp_reigvecs(:, idx(i))
         end do
!         !
!         ! b) KM eigenvalue problem
!         !
!         copyHmat_KM = Hmat_KM(1:max_iter, 1:max_iter) ! zgeev subroutine overwrites matrix
!         call zgeev('N', 'V', max_iter, copyHmat_KM, max_iter, Theta_KM(:, ncycle), leigvecs_KM, max_iter, reigvecs_KM, max_iter, &
!                    workvec_KM, 2*max_iter+1, rwork_KM, info_KM)
!         ! Initialize indices
!         do i = 1, max_iter
!            idx(i) = i
!         end do
!         ! Simple bubble sort
!         do i = 1, max_iter-1
!            k = i
!            do j = i+1, max_iter
!               if (real(Theta_KM(j, ncycle)) < real(Theta_KM(k, ncycle))) k = j
!            end do
!            if (k /= i) then
!               ! swap values
!               temp = Theta_KM(i, ncycle); Theta_KM(i, ncycle) = Theta_KM(k, ncycle); Theta_KM(k, ncycle) = temp
!               ! swap indices
!               temp_idx = idx(i); idx(i) = idx(k); idx(k) = temp_idx
!            end if
!         end do
!         Temp_reigvecs = reigvecs_KM
!         do i = 1, max_iter
!            reigvecs_KM(:, i) = Temp_reigvecs(:, idx(i))
!         end do
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!! End Solve eigenvalue problems !!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!         do i = 1, max_iter
!            Tau2_MK(i, ncycle) = reigvecs_MK(1, i)**2
!         end do
!         do i = 1, max_iter
!            Tau2_KM(i, ncycle) = reigvecs_KM(1, i)**2
!         end do

!         ! Sort
!         print *, 'Printing Ritz values (sorted)'
!         do i = 1, max_iter
!            write(*,'(I5,2X,"(",ES20.12,", ",ES20.12," )")') i, real(Theta(i,ncycle)), aimag(Theta(i,ncycle))
!         end do
!
!         ! Sort
!         print *, 'Printing Ritz values (sorted)'
!         do i = 1, max_iter
!            write(*,'(I5,2X,"(",ES20.12,", ",ES20.12," )")') i, real(Theta2(i,ncycle)), aimag(Theta2(i,ncycle))
!         end do

!         print *, 'Printing Tau2_MK values'
!         do i = 1, max_iter
!            write(*,'(I5,2X,"(",ES20.12,", ",ES20.12," )")') i, real(Tau2_MK(i,ncycle)), aimag(Tau2_MK(i,ncycle))
!         end do
!         print *, 'Printing Tau2_KM values'
!         do i = 1, max_iter
!            write(*,'(I5,2X,"(",ES20.12,", ",ES20.12," )")') i, real(Tau2_KM(i,ncycle)), aimag(Tau2_KM(i,ncycle))
!         end do
         
         ! Compute Ritz vectors
         call zgemm('N', 'N', nxy, max_iter, max_iter, cone, basis_MK, nxy, reigvecs_MK, max_iter, czero, ritzvecs_MK, nxy)
!         call zgemm('N', 'N', nxy, max_iter, max_iter, cone, basis_KM, nxy, reigvecs_KM, max_iter, czero, ritzvecs_KM, nxy)



         ! Want KZ\Omega^{-1/2} and Z\Omega^{1/2}
         ! Z is ritzvecs_MK
         do k = 1, max_iter
            ! Set [Z;0]
            do i=1,nxy
               dRqp_re%m12%elem(i) = dble(ritzvecs_MK(i, k))
               dRqp_re%m21%elem(i) = 0.d0
               dRqp_im%m12%elem(i) = dimag(ritzvecs_MK(i, k))
               dRqp_im%m21%elem(i) = 0.d0
            end do
            ! Do matvec and get [AZ;BZ]
            call matvec(dRqp_re, dRqp_im, temp_re, temp_im)
            do i = 1, nxy
               AZ(i) = dcmplx(temp_re%m12%elem(i), temp_im%m12%elem(i))
               BZ(i) = dcmplx(temp_re%m21%elem(i), temp_im%m21%elem(i))
            end do
            ! KZ = AZ - BZ
            ZZ = AZ - BZ
            KZ(:, k) = ZZ * (real(Theta_MK(k, ncycle)) ** (-0.25))
            Z(:, k) = ritzvecs_MK(:, k) * (real(Theta_MK(k, ncycle)) ** (0.25))
         end do
!         U = (KZ + Z) / 2;
!         V = (kZ - Z) / 2;
         
         ! KZ = U + V
         ! Z = U - V
         do i = 1, nxy
            F20(i) = Fqp%m12%elem(i)
            F02(i) = Fqp%m21%elem(i)
         end do
         call zgemv('C', nxy, max_iter, cone, KZ, nxy, F20+F02, 1, czero, VTF20(:, ncycle), 1)
         call zgemv('C', nxy, max_iter, cone, Z, nxy, F20-F02, 1, czero, UTF02(:, ncycle), 1)
         RF0(:, ncycle) = (VTF20(:, ncycle) + UTF02(:, ncycle)) / 2;
         RF1(:, ncycle) = (VTF20(:, ncycle) - UTF02(:, ncycle)) / 2;






!         ! DMR
!         U = (ritzvecs_KM + ritzvecs_MK) / 2;
!         V = (ritzvecs_KM - ritzvecs_MK) / 2;
!         call zgemv('C', nxy, max_iter, cone, V, nxy, F20, 1, czero, VTF20(:, ncycle), 1)
!         call zgemv('C', nxy, max_iter, cone, U, nxy, F02, 1, czero, UTF02(:, ncycle), 1)

         ! Keeping below because.. original and for writing into file
         ! Strength function
         call contract_bbm(Fqp,dRqp_re,re_s)
         call contract_bbm(Fqp,dRqp_im,im_s)
         re_str(1) = -re_s/pi
         im_str(1) = -im_s/pi

         ! Cross-terms
         if (allocated(g)) then
            do ixterm=1, size(g)
               call contract_bbm(Gqp(ixterm),dRqp_re,re_s)
               call contract_bbm(Gqp(ixterm),dRqp_im,im_s)
               re_str(1+ixterm) = -re_s/pi
               im_str(1+ixterm) = -im_s/pi
            end do
         end if

      end do !end cycle


      time2 = get_timer()




!      do k = 1, nxy
!      do i = 1, nxy
!!         if (i == k) then
!!         dRqp_re%m12%elem(i) = 1.d0
!!         else
!!         dRqp_re%m12%elem(i) = 0.d0
!!         end if
!!         dRqp_re%m21%elem(i) = 0.d0
!         dRqp_re%m12%elem(i) = 0.d0
!         if (i == k) then
!            dRqp_re%m21%elem(i) = 1.d0
!         else
!            dRqp_re%m21%elem(i) = 0.d0
!         endif
!         dRqp_im%m12%elem(i) = 0.d0
!         dRqp_im%m21%elem(i) = 0.d0
!      end do
!      call matvec(dRqp_re, dRqp_im, temp_re, temp_im)
!      do i = 1, nxy
!         BB(i, k) = temp_re%m12%elem(i)
!         AA(i, k) = temp_re%m21%elem(i)
!!         AA(i, k) = temp_re%m12%elem(i)
!!         BB(i, k) = temp_re%m21%elem(i)
!!         AA(i, k) = dcmplx(temp_re%m12%elem(i), temp_im%m12%elem(i))
!!         BB(i, k) = dcmplx(temp_re%m21%elem(i), temp_im%m21%elem(i))
!      end do
!      end do
!      print *, "PRINTING A"
!      do i = 1, nxy
!         print *, AA(i, :)
!      end do
!      print *, "PRINTING B"
!      do i = 1, nxy
!         print *, BB(i, :)
!      end do




!      print *, "PRINTING F20"
!      do i = 1, nxy
!         print *, F20(i)
!      enddo
!      print *, "PRINTING F02"
!      do i = 1, nxy
!         print *, F02(i)
!      enddo


!      print *, "PRINTING THETA"
!      do i = 1, max_iter
!         print *, real(Theta_MK(i,1))
!      enddo
!      print *, "PRINTING RF0"
!      do i = 1, max_iter
!         print *, RF0(i,1)
!      enddo
!      print *, "PRINTING RF1"
!      do i = 1, max_iter
!         print *, RF1(i,1)
!      enddo





      lu1 = 99
      lu2 = 100
      lu3 = 101
      lu4 = 102
      lu5 = 103
      lu6 = 104


      lu7 = 105
      lu8 = 106

      open(unit=lu1, file='theta_mk_values.txt',   status='unknown', action='write')
!      open(unit=lu2, file='theta_km_values.txt',  status='unknown', action='write')
!      open(unit=lu3, file='tau2_mk_values.txt', status='unknown', action='write')
!      open(unit=lu4, file='tau2_km_values.txt', status='unknown', action='write')
!      open(unit=lu5, file='VTF20_values.txt', status='unknown', action='write')
!      open(unit=lu6, file='UTF02_values.txt', status='unknown', action='write')


      open(unit=lu7, file='RF0.txt',   status='unknown', action='write')
      open(unit=lu8, file='RF1.txt',   status='unknown', action='write')

      do i = 1, max_iter
         write(lu1,'(100(1X,ES20.12))') ( real(Theta_MK(i,j)), j=1,ncycle-1 )
!         write(lu2,'(100(1X,ES20.12))') ( real(Theta_KM(i,j)), j=1,ncycle-1 )
!         write(lu3,'(100(1X,ES20.12))') ( real(Tau2_MK(i,j)), j=1,ncycle-1 )
!         write(lu4,'(100(1X,ES20.12))') ( real(Tau2_KM(i,j)), j=1,ncycle-1 )
!         write(lu5,'(100(1X,ES20.12))') ( real(VTF20(i,j)), j=1,ncycle-1 )
!         write(lu6,'(100(1X,ES20.12))') ( real(UTF02(i,j)), j=1,ncycle-1 )
         write(lu7,'(100(1X,ES20.12))') ( RF0(i,j), j=1,ncycle-1 )
         write(lu8,'(100(1X,ES20.12))') ( RF1(i,j), j=1,ncycle-1 )
      end do

      close(lu1)
      close(lu2)
      close(lu3)
      close(lu4)
      close(lu5)
      close(lu6)



      close(lu7)
      close(lu8)



      write(st,'(1x,72("-"))'); call writelog(st)
      if (iter_conv >= 0) then
         write(st,'("  *  LANCZOS iteration converged after   ", i4, " steps.")') (ncycle-2)*max_iter+iter-1; call writelog(st)
      else
         write(st,'("  *  LANCZOS iteration did not converge within step   ", i4, " steps.")') (ncycle-2)*max_iter+iter-1; call writelog(st)
      end if
      write(st,'(1x,72("-"))'); call writelog(st)
      write(st,'(1x,a,es12.3,a)') 'Total CPU time = ',(time2-timei0)/60.0_dp,' minutes' ; call writelog(st)
      call writelog("")

      ! Deallocate memory
      deallocate(re_rand)
      deallocate(basis_MK)
      deallocate(basis_KM)
      deallocate(Hmat_MK)
      deallocate(Hmat_KM)
      deallocate(copyHmat_MK)
      deallocate(copyHmat_KM)
      deallocate(ZZ)
      deallocate(AZ)
      deallocate(BZ)
      deallocate(leigvecs_MK)
      deallocate(leigvecs_KM)
      deallocate(reigvecs_MK)
      deallocate(reigvecs_KM)
      deallocate(workvec_MK)
      deallocate(workvec_KM)
      deallocate(rwork_MK)
      deallocate(rwork_KM)
      deallocate(Theta_MK)
      deallocate(Theta_KM)
      deallocate(Tau2_MK)
      deallocate(Tau2_KM)
      deallocate(idx)
      deallocate(Temp_reigvecs)
      deallocate(ritzvecs_MK)
      deallocate(ritzvecs_KM)

      ! DMR
      deallocate(U)
      deallocate(V)
      deallocate(VTF20)
      deallocate(UTF02)
      deallocate(F20)
      deallocate(F02)



      deallocate(copyHmat)
      deallocate(AA)
      deallocate(BB)



      deallocate(Z)
      deallocate(KZ)
      deallocate(RF0)
      deallocate(RF1)



   end subroutine lanczos



   !---------------------------------------------------------------------------
   ! Like Lanczos subroutine, but solves MK eigenvalue problem only,
   ! and does multiple runs.
   !  
   ! Obtains the Ritz values and the first element squared 
   ! of the corresponding eigenvector of MK eigenvalue problem.
   !---------------------------------------------------------------------------
   subroutine dos
      use pnfam_constants, only : iu, pi, get_timer
      use pnfam_broyden
      implicit none

      ! Parameters to keep from original code
      integer :: iter, ixterm, nxy, i, j
      real(dp) :: re_s, im_s
      real(dp) :: time1, time2, timei, timei0

      integer :: lu1, lu2, lu3, lu4                                     ! Write integer files for Theta_MK, Theta_KM, Tau2_MK, Tau2_KM
      integer :: ncycle                                                 ! Current number of Lanczos cycle

      complex(dp), parameter :: cone=(1.0d0, 0.d0), czero=(0.d0, 0.d0)  ! declaration of complex one and complex zero
      real(dp), external :: dznrm2                                      ! LAPACK function for euclidean norm

      real(dp), allocatable :: re_rand(:)                               ! Random real parts
      integer :: Nvec                                                   ! Number of Lanczos iterations

      ! Krylov subspace parameters
      complex(dp), allocatable :: basis_MK(:, :)                        ! Krylov subspaces
      complex(dp), allocatable :: Hmat_MK(:, :)                         ! Hessenberg matrices
      complex(dp), allocatable :: copyHmat_MK(:, :)                     ! Copy of Hessenberg matrices; for reorthogonalization and zgeev
      real(dp) :: betanrm                                               ! For norm of the basis vectors
      complex(dp), allocatable :: ZZ(:), AZ(:), BZ(:)                   ! For storing intermediate matvecs

      ! zgeev parameters
      complex(dp), allocatable :: leigvecs_MK(:, :)                     ! Left eigenvectors
      complex(dp), allocatable :: reigvecs_MK(:, :)                     ! Right eigenvectors
      complex(dp), allocatable :: workvec_MK(:)                         ! Work vectors
      real(dp), allocatable :: rwork_MK(:)                              ! rwork vectors
      integer :: info_MK                                                ! Info values (=0 indicates convergence)

      ! Eigen parameters
      complex(dp), allocatable :: Theta_MK(:, :)                        ! Ritz values
      integer, allocatable :: idx(:)                                    ! For sorting eigenvalues
      integer :: temp_idx, k                                            ! For sorting eigenvalues
      complex(dp) :: temp                                               ! For sorting eigenvalues
      complex(dp), allocatable :: Temp_reigvecs(:, :)                   ! For sorting eigenvalues
      complex(dp), allocatable :: Tau2_MK(:, :)                         ! Squared first components of the eigenvectors of Hessenberg matrices
      complex(dp), allocatable :: ritzvecs_MK(:, :)                     ! Ritz vectors

      Nvec = 5                                                          ! Number of Lanczos runs
      nxy = size(f%mat%elem)                                            ! Size of A and B

      ! Iteration header
      write(st,'(1x,72("-"))'); call writelog(st)
      write(st,'(1x,3a)') '    i', '              Error', '               Time'; call writelog(st)
      write(st,'(1x,72("-"))'); call writelog(st)

      ! Allocate memory
      allocate(re_rand(nxy))
      allocate(basis_MK(nxy, max_iter+1))
      allocate(Hmat_MK(max_iter+1, max_iter))
      allocate(copyHmat_MK(max_iter+1, max_iter))
      allocate(ZZ(nxy))
      allocate(AZ(nxy))
      allocate(BZ(nxy))
      allocate(leigvecs_MK(max_iter, max_iter))
      allocate(reigvecs_MK(max_iter, max_iter))
      allocate(workvec_MK(2*max_iter+1))
      allocate(rwork_MK(2*max_iter))
      allocate(Theta_MK(max_iter, Nvec))
      allocate(idx(max_iter))
      allocate(Temp_reigvecs(max_iter, max_iter))
      allocate(Tau2_MK(max_iter, Nvec))
      allocate(ritzvecs_MK(nxy, max_iter))

      time1=0; time2=0; timei=0; timei0=0
      Theta_MK = czero
      Tau2_MK = czero

      do ncycle=1,Nvec

         write(st,'(" ****** At cycle number ", i4, " *****.")') ncycle; call writelog(st)

         ! Start the cycle
         time1 = get_timer()
         if (ncycle==1) timei0 = time1

         re_rand = 0.d0
         basis_MK = czero
         Hmat_MK = czero
         copyHmat_MK = czero
         leigvecs_MK = czero
         reigvecs_MK = czero
         workvec_MK = czero
         rwork_MK = 0.d0
         Temp_reigvecs = czero
         ritzvecs_MK = czero

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!! Initialize !!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         call random_seed()
         call random_number(re_rand)
         do i = 1, nxy
            basis_MK(i, 1) = dcmplx(re_rand(i), 0.d0)
         end do
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!! End of Initialization !!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!! Normalize the first basis vector !!!!
         !    For MK problem, normalize w.r.t. K  !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         ! Set [Z;0]
         do i=1,nxy
            dRqp_re%m12%elem(i) = dble(basis_MK(i, 1))
            dRqp_re%m21%elem(i) = 0.d0
            dRqp_im%m12%elem(i) = dimag(basis_MK(i, 1))
            dRqp_im%m21%elem(i) = 0.d0
         end do
         ! Do matvec and get [AZ;BZ]
         call matvec(dRqp_re, dRqp_im, temp_re, temp_im)
         do i = 1, nxy
            AZ(i) = dcmplx(temp_re%m12%elem(i), temp_im%m12%elem(i))
            BZ(i) = dcmplx(temp_re%m21%elem(i), temp_im%m21%elem(i))
         end do
         ! KZ = AZ - BZ
         ZZ = AZ - BZ
         betanrm = sqrt(real(dot_product(conjg(basis_MK(:, 1)), ZZ)))
         basis_MK(:, 1) = basis_MK(:, 1) / betanrm

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! End of first basis vector normalization!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! Record the iteration time
         time2 = get_timer()
         timei = (time2-time1)

         ! Record the initial result
         write(st,'(1x, i4, 2x, f20.10, 2x, f15.3)') 0, betanrm, timei; call writelog(st)
         ! End of initialization step

         ! Main lanczos loop
         do iter=1,max_iter

            ! Start the iteration
            time1 = get_timer()

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!! Generate Krylov subspace for MK !!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !
            ! Set [Z;0]
            do i=1,nxy
               dRqp_re%m12%elem(i) = dble(basis_MK(i, iter))
               dRqp_re%m21%elem(i) = 0.d0
               dRqp_im%m12%elem(i) = dimag(basis_MK(i, iter))
               dRqp_im%m21%elem(i) = 0.d0
            end do
            ! Do matvec and get [AZ;BZ]
            call matvec(dRqp_re, dRqp_im, temp_re, temp_im)
            do i = 1, nxy
               AZ(i) = dcmplx(temp_re%m12%elem(i), temp_im%m12%elem(i))
               BZ(i) = dcmplx(temp_re%m21%elem(i), temp_im%m21%elem(i))
            end do
            ! KZ = AZ-BZ
            ZZ = AZ - BZ
            ! Set [KZ;0]
            do i=1,nxy
               dRqp_re%m12%elem(i) = dble(ZZ(i))
               dRqp_re%m21%elem(i) = 0.d0
               dRqp_im%m12%elem(i) = dimag(ZZ(i))
               dRqp_im%m21%elem(i) = 0.d0
            end do
            ! Do matvec and get [AKZ;BKZ]
            call matvec(dRqp_re, dRqp_im, temp_re, temp_im)
            do i = 1, nxy
               AZ(i) = dcmplx(temp_re%m12%elem(i), temp_im%m12%elem(i))
               BZ(i) = dcmplx(temp_re%m21%elem(i), temp_im%m21%elem(i))
            end do
            ! MKZ = AKZ+BKZ
            ZZ = AZ + BZ
            ! Store into basis
            do i = 1,nxy
               basis_MK(i, iter+1) = ZZ(i)
            end do
            !
            ! b) K-inner product of MKZ
            !
            ! Set [MKZ;0]
            do i=1,nxy
               dRqp_re%m12%elem(i) = dble(ZZ(i))
               dRqp_re%m21%elem(i) = 0.d0
               dRqp_im%m12%elem(i) = dimag(ZZ(i))
               dRqp_im%m21%elem(i) = 0.d0
            end do
            ! Do matvec and get [AMKZ;BMKZ]
            call matvec(dRqp_re, dRqp_im, temp_re, temp_im)
            do i = 1, nxy
               AZ(i) = dcmplx(temp_re%m12%elem(i), temp_im%m12%elem(i))
               BZ(i) = dcmplx(temp_re%m21%elem(i), temp_im%m21%elem(i))
            end do
            ! KMKZ = AMKZ-BMKZ
            ZZ = AZ - BZ
            !
            ! c) Orthogonalize w.r.t. previous Krylov basis vectors
            !
            call zgemv('C', nxy, iter, cone, basis_MK, nxy, ZZ, 1, czero, Hmat_MK(1:iter, iter), 1)
            call zgemv('N', nxy, iter, -cone, basis_MK, nxy, Hmat_MK(1:iter, iter), 1, cone, basis_MK(:, iter+1), 1)
            !
            ! d) Reorthogonalize
            !
            !    d)a) K-inner product of f:=basis_MK(:, iter+1)
            !
            ! Set [f;0]
            do i=1,nxy
               dRqp_re%m12%elem(i) = dble(basis_MK(i, iter+1))
               dRqp_re%m21%elem(i) = 0.d0
               dRqp_im%m12%elem(i) = dimag(basis_MK(i, iter+1))
               dRqp_im%m21%elem(i) = 0.d0
            end do
            ! Do matvec and get [Af;Bf]
            call matvec(dRqp_re, dRqp_im, temp_re, temp_im)
            do i = 1, nxy
               AZ(i) = dcmplx(temp_re%m12%elem(i), temp_im%m12%elem(i))
               BZ(i) = dcmplx(temp_re%m21%elem(i), temp_im%m21%elem(i))
            end do
            ! Kf = Af-Bf
            ZZ = AZ - BZ
            !
            !    d)b) Orthogonalize w.r.t. previous Krylov basis vectors
            !
            call zgemv('C', nxy, iter, cone, basis_MK, nxy, ZZ, 1, czero, copyHmat_MK(1:iter, iter), 1)
            call zgemv('N', nxy, iter, -cone, basis_MK, nxy, copyHmat_MK(1:iter, iter), 1, cone, basis_MK(:, iter+1), 1)
            Hmat_MK(1:iter, iter) = Hmat_MK(1:iter, iter) + copyHmat_MK(1:iter, iter)
            !
            ! e) Normalize current basis vector w.r.t. K
            !
            ! Set [Z;0]
            do i=1,nxy
               dRqp_re%m12%elem(i) = dble(basis_MK(i, iter+1))
               dRqp_re%m21%elem(i) = 0.d0
               dRqp_im%m12%elem(i) = dimag(basis_MK(i, iter+1))
               dRqp_im%m21%elem(i) = 0.d0
            end do
            ! Do matvec and get [AZ;BZ]
            call matvec(dRqp_re, dRqp_im, temp_re, temp_im)
            do i = 1, nxy
               AZ(i) = dcmplx(temp_re%m12%elem(i), temp_im%m12%elem(i))
               BZ(i) = dcmplx(temp_re%m21%elem(i), temp_im%m21%elem(i))
            end do
            ! KZ = AZ-BZ
            ZZ = AZ - BZ
            betanrm = sqrt(real(dot_product(conjg(basis_MK(:, iter+1)), ZZ)))
            basis_MK(:, iter+1) = basis_MK(:, iter+1) / betanrm
            Hmat_MK(iter+1, iter) = betanrm
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! End Krylov subspace generation for MK !!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 

            ! Record the iteration time
            time2 = get_timer()
            timei = (time2-time1)

            ! Record the iteration (or show where we are starting from)
            write(st,'(1x, i4, 2x, f20.10, 2x, f15.3)') iter, betanrm, timei; call writelog(st)

         end do ! end inner iteration

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!! Solve eigenvalue problems !!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         ! a) MK eigenvalue problem
         !
         copyHmat_MK = Hmat_MK(1:max_iter, 1:max_iter) ! zgeev subroutine overwrites matrix
         call zgeev('N', 'V', max_iter, copyHmat_MK, max_iter, Theta_MK(:, ncycle), leigvecs_MK, max_iter, reigvecs_MK, max_iter, &
                    workvec_MK, 2*max_iter+1, rwork_MK, info_MK)
         ! Initialize indices
         do i = 1, max_iter
            idx(i) = i
         end do
         ! Simple bubble sort
         do i = 1, max_iter-1
            k = i
            do j = i+1, max_iter
               if (real(Theta_MK(j, ncycle)) < real(Theta_MK(k, ncycle))) k = j
            end do
            if (k /= i) then
               ! swap values
               temp = Theta_MK(i, ncycle); Theta_MK(i, ncycle) = Theta_MK(k, ncycle); Theta_MK(k, ncycle) = temp
               ! swap indices
               temp_idx = idx(i); idx(i) = idx(k); idx(k) = temp_idx
            end if
         end do
         Temp_reigvecs = reigvecs_MK
         do i = 1, max_iter
            reigvecs_MK(:, i) = Temp_reigvecs(:, idx(i))
         end do
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!! End Solve eigenvalue problems !!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         do i = 1, max_iter
            Tau2_MK(i, ncycle) = reigvecs_MK(1, i)**2
         end do

!         ! Sort
!         print *, 'Printing Ritz values (sorted)'
!         do i = 1, max_iter
!            write(*,'(I5,2X,"(",ES20.12,", ",ES20.12," )")') i, real(Theta(i,ncycle)), aimag(Theta(i,ncycle))
!         end do
!
!         ! Sort
!         print *, 'Printing Ritz values (sorted)'
!         do i = 1, max_iter
!            write(*,'(I5,2X,"(",ES20.12,", ",ES20.12," )")') i, real(Theta2(i,ncycle)), aimag(Theta2(i,ncycle))
!         end do

!         print *, 'Printing Tau2_MK values'
!         do i = 1, max_iter
!            write(*,'(I5,2X,"(",ES20.12,", ",ES20.12," )")') i, real(Tau2_MK(i,ncycle)), aimag(Tau2_MK(i,ncycle))
!         end do
!         print *, 'Printing Tau2_KM values'
!         do i = 1, max_iter
!            write(*,'(I5,2X,"(",ES20.12,", ",ES20.12," )")') i, real(Tau2_KM(i,ncycle)), aimag(Tau2_KM(i,ncycle))
!         end do
         
         ! Compute Ritz vectors
         call zgemm('N', 'N', nxy, max_iter, max_iter, cone, basis_MK, nxy, reigvecs_MK, max_iter, czero, ritzvecs_MK, nxy)


         ! Keeping below because.. original and for writing into file
         ! Strength function
         call contract_bbm(Fqp,dRqp_re,re_s)
         call contract_bbm(Fqp,dRqp_im,im_s)
         re_str(1) = -re_s/pi
         im_str(1) = -im_s/pi

         ! Cross-terms
         if (allocated(g)) then
            do ixterm=1, size(g)
               call contract_bbm(Gqp(ixterm),dRqp_re,re_s)
               call contract_bbm(Gqp(ixterm),dRqp_im,im_s)
               re_str(1+ixterm) = -re_s/pi
               im_str(1+ixterm) = -im_s/pi
            end do
         end if

      end do !end cycle



      lu1 = 99
      lu2 = 100
      lu3 = 101
      lu4 = 102

      open(unit=lu1, file='theta_mk_values.txt',   status='unknown', action='write')
      open(unit=lu3, file='tau2_mk_values.txt', status='unknown', action='write')

      do i = 1, max_iter
         write(lu1,'(100(1X,ES20.12))') ( real(Theta_MK(i,j)), j=1,ncycle-1 )
         write(lu3,'(100(1X,ES20.12))') ( real(Tau2_MK(i,j)), j=1,ncycle-1 )
      end do

      close(lu1)
      close(lu3)



      write(st,'(1x,72("-"))'); call writelog(st)
      if (iter_conv >= 0) then
         write(st,'("  *  DOS iteration converged after   ", i4, " steps.")') (ncycle-2)*max_iter+iter-1; call writelog(st)
      else
         write(st,'("  *  DOS iteration did not converge within step   ", i4, " steps.")') (ncycle-2)*max_iter+iter-1; call writelog(st)
      end if
      write(st,'(1x,72("-"))'); call writelog(st)
      write(st,'(1x,a,es12.3,a)') 'Total CPU time = ',(time2-timei0)/60.0_dp,' minutes' ; call writelog(st)
      call writelog("")

      ! Deallocate memory
      deallocate(re_rand)
      deallocate(basis_MK)
      deallocate(Hmat_MK)
      deallocate(copyHmat_MK)
      deallocate(ZZ)
      deallocate(AZ)
      deallocate(BZ)
      deallocate(leigvecs_MK)
      deallocate(reigvecs_MK)
      deallocate(workvec_MK)
      deallocate(rwork_MK)
      deallocate(Theta_MK)
      deallocate(Tau2_MK)
      deallocate(idx)
      deallocate(Temp_reigvecs)
      deallocate(ritzvecs_MK)


   end subroutine dos





   !---------------------------------------------------------------------------
   ! Allocates and initializes the data structures used by pnfam_solve.
   !---------------------------------------------------------------------------
   subroutine init_pnfam_solver
      use pnfam_constants, only : IT_NEUTRON, IT_PROTON
      use pnfam_broyden, only : qrpa_broin, qrpa_broout, nbroy, broyden_history_size

      implicit none

      integer :: ixterm, i
      logical :: a11,a12,a21,a22

      !---------------------------------------------------------------------------
      ! General initializations
      !---------------------------------------------------------------------------
      ! The complex strength function (array [strength, xterms])
      if (allocated(re_str)) deallocate(re_str, im_str)
      allocate(re_str(1:nxterms+1), im_str(1:nxterms+1))
      re_str(:) = 0; im_str(:) = 0

      ! While the various FAM matrices may have different non-trivial block structures
      ! they still share the same number of elements in the non-trivial blocks
      nxy = size(f%mat%elem)
      bminus = f%beta_minus

      if (ft_active .or. hfb_blo_active) then
         use_diagonal_blocks = .true.
         nbroy = nxy*8
      else
         use_diagonal_blocks = .false.
         nbroy = nxy*4
      end if

      ! Broyden arrays
      if (allocated(qrpa_broin)) deallocate(qrpa_broin,qrpa_broout)
      allocate(qrpa_broin(nbroy),qrpa_broout(nbroy))
      qrpa_broin = 0; qrpa_broout=0
      si = 1
      ! Don't do broyden if no residual interaction, should just need 1 iteration
      if (abs(quench_residual_int) < 1e-10) then
         broyden_history_size = -1
      end if

      !---------------------------------------------------------------------------
      ! Allocate and initialize big-blockmatrix structure
      !---------------------------------------------------------------------------
      ! Which blocks are active
      a11 = .false.; a12 = .true.
      a21 = .true. ; a22 = .false.
      if (use_diagonal_blocks) then
         a11 = .true.; a22 = .true.
      end if

      ! Single particle external field (F_ij*(a^\dagger a) --> m11)
      ! Represents.. : [ f_pn, 0 ; 0, 0 ], for 2bc w/ pairing: [f11/2, f12; -f12*/2, -f11*]
      ! Stores...... : [ f_pn, _ ; _, _ ], for 2bc w/ pairing: [f11/2, f12;  f12/2,   f11 ]
      Fsp%m11 = f%mat
      ! THIS DOES NOT WORK B/C IT SPOILS THE 1 BLOCK PER ROW+COL STRUCTURE IN FQP
      if (allocated(f%mat12%elem)) then
         Fsp%m12 = f%mat12
         Fsp%m21 = f%mat12
         Fsp%m22 = f%mat
         Fsp%m11%elem = f%mat%elem*0.5_dp
         Fsp%m22%elem = f%mat%elem*0.5_dp
         call set_sign_bbm(Fsp, 1d0,1d0,-1d0,-1d0)
      end if
      if (allocated(g)) then
         if (allocated(Gsp)) deallocate(Gsp)
         allocate(Gsp(size(g)), Gqp(size(g)))
         do ixterm=1, size(g)
            Gsp(ixterm)%m11 = g(ixterm)%mat
         end do
      end if

      ! Bogoliubov transform (real matrices)
      ! Represents.. : [ U, V* ; V, U* ]
      ! Stores...... : [ U, V  ; V, U  ]
      Wn%m11 = Un; Wn%m12 = Vn
      Wn%m21 = Vn; Wn%m22 = Un
      Wp%m11 = Up; Wp%m12 = Vp
      Wp%m21 = Vp; Wp%m22 = Up

      ! External fields in QP basis
      ! Represents.. : [ F11, F20 ; -F02, -F11t ]
      ! Stores...... : [ F11, F20 ;  F02,  F11t ]
      call allocate_bbm(Fqp, a11,a12,a21,a22, nxy)
      call set_sign_bbm(Fqp, 1d0,1d0,-1d0,-1d0)
      if (allocated(g)) then
         do ixterm=1, size(g)
            call allocate_bbm(Gqp(ixterm), a11,a12,a21,a22, nxy)
            call set_sign_bbm(Gqp(ixterm), 1d0,1d0,-1d0,-1d0)
         end do
      end if

      ! Perturbed hamiltonian in QP basis
      ! Represents.. : [ H11, H20 ; -H02, -H11t ]
      ! Stores...... : [ H11, H20 ;  H02,  H11t ]
      call allocate_bbm(dHqp_re, a11,a12,a21,a22, nxy)
      call allocate_bbm(dHqp_im, a11,a12,a21,a22, nxy); dHqp_im%imag=.true.
      call set_sign_bbm(dHqp_re, 1d0,1d0,-1d0,-1d0)
      call set_sign_bbm(dHqp_im, 1d0,1d0,-1d0,-1d0)

      ! Perturbed generalized density in QP basis
      ! Represents.. : [ P, X ; -Y, -Q ]
      ! Stores...... : [ P, X ;  Y,  Q ]
      call allocate_bbm(dRqp_re, a11,a12,a21,a22, nxy)
      call allocate_bbm(dRqp_im, a11,a12,a21,a22, nxy); dRqp_im%imag=.true.
      call set_sign_bbm(dRqp_re, 1d0,1d0,-1d0,-1d0)
      call set_sign_bbm(dRqp_im, 1d0,1d0,-1d0,-1d0)

      ! FAM greens function: g = -1/(E1 +/- E2 +/- (re(eqrpa)+im(eqrpa)*i))
      ! Represents.. : diag( gP, gX, gY, gQ )
      ! Stores...... : [ gP, gX ; gY, gQ ]
      call allocate_bbm(Greens_re, a11,a12,a21,a22, nxy)
      call allocate_bbm(Greens_im, a11,a12,a21,a22, nxy); Greens_im%imag=.true.

      ! Statistical occupation matrix: 1ff=(1-f1-f2) and ff=(f2-f1)
      ! Represents.. : diag( T1ff Tff Tff T1ff )
      ! Stores...... : [ T1ff, Tff ; Tff, T1ff ]
      if (use_diagonal_blocks) call allocate_bbm(T, a11,a12,a21,a22, nxy)

      ! Perturbed hamiltonian in SP basis
      ! Represents.. : [ h_pn, Delta(+)_pn ; -Delta(-)*_pn, -(h^T)_pn ]
      ! Stores...... : [ h_pn, Delta(+)_pn ;  Delta(-) _pn,   h_np    ]
      call allocate_bbm(dHsp_re, .true.,.true.,.true.,.true., nxy)
      call allocate_bbm(dHsp_im, .true.,.true.,.true.,.true., nxy); dHsp_im%imag=.true.
      call  set_sign_bbm(dHsp_re,+1d0,+1d0,-1d0,-1d0)
      call set_trans_bbm(dHsp_re,'n', 'n', 'n', 't')
      call  set_sign_bbm(dHsp_im,+1d0,+1d0,+1d0,-1d0)
      call set_trans_bbm(dHsp_im,'n', 'n', 'n', 't')

      ! Perturbed density in SP basis (Note kap is antisymmetric so np<-->pn gives a minus sign)
      ! Represents.. : [ rho_pn,  kap(+)_pn ; -kap(-)*_pn, -(rho^T)_pn ]
      ! Stores...... : [ rho_pn,  kap(+)_np ;  kap(-) _np,   rho_np    ]
      !            ?=? [ rho_pn, rhotilde(+)_pn ; rhotilde(-)_pn, rho_np ]
      call allocate_bbm(dRsp_re, .true.,.true.,.true.,.true., nxy)
      call allocate_bbm(dRsp_im, .true.,.true.,.true.,.true., nxy); dRsp_im%imag=.true.
      call  set_sign_bbm(dRsp_re,+1d0,-1d0,+1d0,-1d0)
      call set_trans_bbm(dRsp_re,'n', 't', 't', 't')
      call  set_sign_bbm(dRsp_im,+1d0,-1d0,-1d0,-1d0)
      call set_trans_bbm(dRsp_im,'n', 't', 't', 't')

      !---------------------------------------------------------------------------
      ! Transform external field to QP basis
      !---------------------------------------------------------------------------
      if (bminus) then
         call triprod_bbm('t',Wp,'n',Fsp,'n',Wn,Fqp)
      else
         call triprod_bbm('t',Wn,'n',Fsp,'n',Wp,Fqp)
      end if
      if (allocated(g)) then
         do ixterm=1, size(g)
            if (bminus) then
               call triprod_bbm('t',Wp,'n',Gsp(ixterm),'n',Wn,Gqp(ixterm))
            else
               call triprod_bbm('t',Wn,'n',Gsp(ixterm),'n',Wp,Gqp(ixterm))
            end if
         end do
      end if

      !---------------------------------------------------------------------------
      ! Initialize individual block structures based on external field
      ! (For arrays whose block structure is not determined by routine triprod)
      !---------------------------------------------------------------------------
      ! dRqp has same structure as Fqp
      call copy_block_structure_bbm(Fqp, dRqp_re)
      call copy_block_structure_bbm(Fqp, dRqp_im)

      ! dHqp has same structure as Fqp
      call copy_block_structure_bbm(Fqp, dHqp_re)
      call copy_block_structure_bbm(Fqp, dHqp_im)

      ! Greens has same structure as dRqp
      call copy_block_structure_bbm(dRqp_re, Greens_re)
      call copy_block_structure_bbm(dRqp_im, Greens_im)

      ! Statmatrix has same structure as Fqp
      if (use_diagonal_blocks) call copy_block_structure_bbm(Fqp, T)

      ! dh_pn has same structure as Fsp
      call copy_block_structure(Fsp%m11, dHsp_re%m11)
      call copy_block_structure(Fsp%m11, dHsp_im%m11)
      ! dDelta(+)_pn has same structure as X
      call copy_block_structure(dRqp_re%m12, dHsp_re%m12)
      call copy_block_structure(dRqp_im%m12, dHsp_im%m12)
      ! dDelta(-)_pn has same structure as Y
      call copy_block_structure(dRqp_re%m21, dHsp_re%m21)
      call copy_block_structure(dRqp_im%m21, dHsp_im%m21)
      ! dh_np has same structure as Fsp
      call copy_block_structure(Fsp%m11, dHsp_re%m22)
      call copy_block_structure(Fsp%m11, dHsp_im%m22)

      !---------------------------------------------------------------------------
      ! Calculate 2QP quantities
      !---------------------------------------------------------------------------
      ! Greens function: -1/(E1 +/- E2 +/- (re(eqrpa)+im(eqrpa)*i))
      ! Calculate the 2QP matrix elements M_{ij} = a*( b*f1_i + c*f2_j + e )**d
      
      ! DMR
      ! Changed the first value to take on Greens_sign value, which is set by preprocessor conditionals
      !
      if (bminus) then
         call matrix_2qp(Greens_sign,+1,+1,Greens_power,-real_eqrpa,-imag_eqrpa,.true., Ep,En,Greens_re%m12)
         call matrix_2qp(Greens_sign,+1,+1,Greens_power,-real_eqrpa,-imag_eqrpa,.false.,Ep,En,Greens_im%m12)
         call matrix_2qp(Greens_sign,+1,+1,Greens_power,+real_eqrpa,+imag_eqrpa,.true., Ep,En,Greens_re%m21)
         call matrix_2qp(Greens_sign,+1,+1,Greens_power,+real_eqrpa,+imag_eqrpa,.false.,Ep,En,Greens_im%m21)
         if (use_diagonal_blocks) then
         call matrix_2qp(Greens_sign,+1,-1,Greens_power,-real_eqrpa,-imag_eqrpa,.true., Ep,En,Greens_re%m11)
         call matrix_2qp(Greens_sign,+1,-1,Greens_power,-real_eqrpa,-imag_eqrpa,.false.,Ep,En,Greens_im%m11)
         call matrix_2qp(Greens_sign,+1,-1,Greens_power,+real_eqrpa,+imag_eqrpa,.true., Ep,En,Greens_re%m22)
         call matrix_2qp(Greens_sign,+1,-1,Greens_power,+real_eqrpa,+imag_eqrpa,.false.,Ep,En,Greens_im%m22)
         end if
      else
         call matrix_2qp(-1,+1,+1,-1,-real_eqrpa,-imag_eqrpa,.true., En,Ep,Greens_re%m12)
         call matrix_2qp(-1,+1,+1,-1,-real_eqrpa,-imag_eqrpa,.false.,En,Ep,Greens_im%m12)
         call matrix_2qp(-1,+1,+1,-1,+real_eqrpa,+imag_eqrpa,.true., En,Ep,Greens_re%m21)
         call matrix_2qp(-1,+1,+1,-1,+real_eqrpa,+imag_eqrpa,.false.,En,Ep,Greens_im%m21)
         if (use_diagonal_blocks) then
         call matrix_2qp(-1,+1,-1,-1,-real_eqrpa,-imag_eqrpa,.true., En,Ep,Greens_re%m11)
         call matrix_2qp(-1,+1,-1,-1,-real_eqrpa,-imag_eqrpa,.false.,En,Ep,Greens_im%m11)
         call matrix_2qp(-1,+1,-1,-1,+real_eqrpa,+imag_eqrpa,.true., En,Ep,Greens_re%m22)
         call matrix_2qp(-1,+1,-1,-1,+real_eqrpa,+imag_eqrpa,.false.,En,Ep,Greens_im%m22)
         end if
      end if

      ! Statistical occupation matrix: (1-f1-f2) and (f2-f1)
      if (use_diagonal_blocks) then
         call set_val_bbm(T,0d0)
         if (bminus) then
            call matrix_2qp(+1,-1,-1,+1, 1.0_dp,0.0_dp, .true., qp_fp,qp_fn, T%m12)
            call matrix_2qp(+1,-1,-1,+1, 1.0_dp,0.0_dp, .true., qp_fp,qp_fn, T%m21)
            call matrix_2qp(+1,-1,+1,+1, 0.0_dp,0.0_dp, .true., qp_fp,qp_fn, T%m11)
            call matrix_2qp(+1,-1,+1,+1, 0.0_dp,0.0_dp, .true., qp_fp,qp_fn, T%m22)
         else
            call matrix_2qp(+1,-1,-1,+1, 1.0_dp,0.0_dp, .true., qp_fn,qp_fp, T%m12)
            call matrix_2qp(+1,-1,-1,+1, 1.0_dp,0.0_dp, .true., qp_fn,qp_fp, T%m21)
            call matrix_2qp(+1,-1,+1,+1, 0.0_dp,0.0_dp, .true., qp_fn,qp_fp, T%m11)
            call matrix_2qp(+1,-1,+1,+1, 0.0_dp,0.0_dp, .true., qp_fn,qp_fp, T%m22)
         end if
      end if

   end subroutine init_pnfam_solver


   subroutine deallocate_pnfam_solver
      implicit none
      integer :: ixterm

      ! External fields in QP basis
      call deallocate_bbm(Fqp)
      if (allocated(g)) then
         do ixterm=1, size(g)
            call deallocate_bbm(Gqp(ixterm))
         end do
      end if

      ! Perturbed hamiltonian in QP basis
      call deallocate_bbm(dHqp_re)
      call deallocate_bbm(dHqp_im)

      ! Perturbed generalized density in QP basis
      call deallocate_bbm(dRqp_re)
      call deallocate_bbm(dRqp_im)

      ! FAM greens function: g = -1/(E1 +/- E2 +/- (re(eqrpa)+im(eqrpa)*i))
      call deallocate_bbm(Greens_re)
      call deallocate_bbm(Greens_im)

      ! Statistical occupation matrix: 1ff=(1-f1-f2) and ff=(f2-f1)
      if (use_diagonal_blocks) call deallocate_bbm(T)

      ! Perturbed hamiltonian in SP basis
      call deallocate_bbm(dHsp_re)
      call deallocate_bbm(dHsp_im)

      ! Perturbed density in SP basis (Note kap is antisymmetric so np<-->pn gives a minus sign)
      call deallocate_bbm(dRsp_re)
      call deallocate_bbm(dRsp_im)

   end subroutine deallocate_pnfam_solver

   !----------------------------------------------------------------------------
   ! Calculate the 2QP matrix elements M_{ij} = a*( b*f1_i + c*f2_j + e )**d
   ! Input arrays are assumed size(f) = dqp (e.g. Ep, En, qp_fp, qp_fn)
   ! - Greens(X): a= -1, b= +1, c= +1, d= -1, e= -(w+hw*i), f1= Ep/En,       f2= En/Ep
   ! - Greens(Y): a= -1, b= +1, c= +1, d= -1, e= +(w+hw*i), f1= Ep/En,       f2= En/Ep
   ! - Greens(P): a= -1, b= +1, c= -1, d= -1, e= -(w+hw*i), f1= Ep/En,       f2= En/Ep
   ! - Greens(Q): a= -1, b= +1, c= -1, d= -1, e= +(w+hw*i), f1= Ep/En,       f2= En/Ep
   ! - T(X/Y)   : a= +1, b= -1, c= -1, d= +1, e= +1,        f1= qp_fp/qp_fn, f2= qp_fn/qp_fp
   ! - T(P/Q)   : a= +1, b= +1, c= -1, d= +1, e=  0,        f1= qp_fp/qp_fn, f2= qp_fn/qp_fp
   !----------------------------------------------------------------------------
   subroutine matrix_2qp(a,b,c,d,ree,ime,re,f1,f2,mat)

      implicit none

      integer, intent(in) :: a,b,c,d
      real(dp), intent(in) :: f1(:), f2(:), ree, ime
      logical, intent(in) :: re
      type(blockmatrix), intent(inout) :: mat

      logical  :: l1, l2
      integer  :: ibr, ibc, ipt, i1, i2, i1a, i1b, i2a, i2b
      complex(dp) :: aux

      mat%elem = 0

      ipt = 1
      do ibr=1, nb
         ibc = mat%ir2c(ibr)
         if (ibc == 0) cycle
            i1a = isstart(ibr) ; i1b = i1a + db(ibr) - 1
            i2a = isstart(ibc) ; i2b = i2a + db(ibc) - 1
            do i2 = i2a,i2b ! true column index
            do i1 = i1a,i1b ! true row index
               aux = a*(b*f1(i1) + c*f2(i2) + cmplx(ree,ime,dp))**d
               if (re) then
                  mat%elem(ipt) = real(aux,dp)
               else
                  mat%elem(ipt) = aimag(aux)
               end if
               ipt = ipt + 1
            end do
            end do
      end do

   end subroutine matrix_2qp


   !----------------------------------------------------------------------------
   ! Does matrix-vector multiplication:
   ! Y_ = H * X_
   !----------------------------------------------------------------------------
   subroutine matvec(X_re, X_im, Y_re, Y_im)
      implicit none

      type(bigblockmatrix), intent(in) :: X_re, X_im
      type(bigblockmatrix), intent(inout) :: Y_re, Y_im
      
      if (abs(quench_residual_int) < 1e-10) then
         call set_val_bbm(dHqp_im,0d0)
         call set_val_bbm(dHqp_re,0d0)
      else
         ! Compute the perturbed densities in SP space
         if (bminus) then
            call triprod_bbm('n',Wp,'n',X_re,'t',Wn,dRsp_re)
            call triprod_bbm('n',Wp,'n',X_im,'t',Wn,dRsp_im)
         else
            call triprod_bbm('n',Wn,'n',X_re,'t',Wp,dRsp_re)
            call triprod_bbm('n',Wn,'n',X_im,'t',Wp,dRsp_im)
         end if
         call calc_hamiltonian(&
                                     dRsp_re%m11, dRsp_im%m11, dRsp_re%m12, dRsp_im%m12,&
                                     dRsp_re%m22, dRsp_im%m22, dRsp_re%m21, dRsp_im%m21,&
                                     dHsp_re%m11, dHsp_im%m11, dHsp_re%m12, dHsp_im%m12,&
                                     dHsp_re%m22, dHsp_im%m22, dHsp_re%m21, dHsp_im%m21 &
                              )
         ! Transform the Hamiltonian into quasiparticle basis
         if (bminus) then
            call triprod_bbm('t',Wp,'n',dHsp_re,'n',Wn,dHqp_re)
            call triprod_bbm('t',Wp,'n',dHsp_im,'n',Wn,dHqp_im)
         else
            call triprod_bbm('t',Wn,'n',dHsp_re,'n',Wp,dHqp_re)
            call triprod_bbm('t',Wn,'n',dHsp_im,'n',Wp,dHqp_im)
         end if
         ! Quench the residual interaction
         call scalar_mult_bbm(dHqp_re,quench_residual_int)
         call scalar_mult_bbm(dHqp_im,quench_residual_int)
      end if ! End computation of the residual interaction
      Y_re = X_re
      Y_im = X_im
      call complex_emult_bbm(Greens_re, Greens_im, X_re, X_im,&
                                            Y_re, Y_im)
      call add_bbm(dHqp_re, Y_re)
      call add_bbm(dHqp_im, Y_im)
   end subroutine matvec



   ! Get the external field objects
   ! (depends on namelist and basis wfs)
   !
   ! init_external_field:
   ! sets: label, beta_minus, parity, K, rank(J)
   ! calls: init_fam_mapping(depends on blockmatrix, K, parity, Kpi-of-blocks)
   ! calls: ext_field_operator (calcs operator, depends on basis, block structure)
   !
   ! setup_cross_terms just calls init_ext_field for the various cross terms
   !----------------------------------------------
   subroutine setup_extfield()
      use pnfam_constants, only : get_timer
      use pnfam_extfield, only : init_external_field, setup_crossterms, set_use_2bc
      use pnfam_storage, only : fam_io
      use pnfam_extfield_2bc, only : effective_2bc_extfield, debug
      use type_extfield_2bc, only : init_extfield_2bc_type
      use type_gamdel_2bc, only : calc_totgamdel
      implicit none
      integer :: iexit, i
      real(dp), allocatable, dimension(:) :: elem_tmp
      real(dp) :: time1, time2

      iexit=1

      ! Calculate 1 body external field
      !11/5/23: since the rest of the code for the two body current mode is designed to handle the GT part, I won't mess with that and 
      !instead just carve out an exception for the P operator with two body current. also only run the two body part if operator = GT.
      !11/29: Add in two body current mode to setup_crossterms. If it's 0 then nothing should change, since set_use_2bc does nothing if it's 0.
      if ((operator_name == 'P' .or. operator_name == 'PS0' .or. operator_name == 'RS0' .or. operator_name == 'RS1' .or. operator_name == 'RS2' .or. operator_name == 'RS0I' .or. operator_name == 'RS1I') .and. two_body_current_mode /= 0) then
         call init_external_field(beta_type=beta_type, label=operator_name, k=operator_k, op=f, use_2bc = two_body_current_mode)
      else
         call init_external_field(beta_type=beta_type, label=operator_name, k=operator_k, op=f)
      end if
      nxterms = 0
      if (compute_crossterms) then
         call setup_crossterms(op=f, crossterms=g, n=nxterms, use_2bc = two_body_current_mode)
      end if

      ! For production runs (no debug) we only calculate/read the 2BC Yukawa part.
      ! The 1BC and 2BC contact part are calculated using the 1BC code and added in.
      ! use_2bc code is:
      ! * 1st digit: 1=1BC+2BC, 2=2BC only
      ! * 2nd digit: 1=Full FAM, 2=SNM+LDA, 3=ASNM+LDA, 4=DME(exc,dir=0), 5=DME(exc)+FAM(dir)
      ! * 3rd digit: 1=Gamma, 2=Delta, 3=Gamma+Delta
      if (two_body_current_mode /= 0 .and. operator_name == 'GT') then
         ! Set the 2BC part of the external field type
         call set_use_2bc(f, two_body_current_mode)
         if (f%use_2bc(2) == 0) return !not calculating 2bc for the GT operator.
         ! Store the 1-body part (if doing 1BC+2BC)
         if (allocated(elem_tmp)) deallocate(elem_tmp)
         allocate(elem_tmp(size(f%mat%elem)))
         elem_tmp = 0
         if (f%use_2bc(1) == 1) then
            elem_tmp = -f%mat%elem ! Park Sign Convention: -sigma tau + 2BC
         end if

         ! Calculate 2-body contact part using (sigma tau rho)  (if not debug).
         ! Also calculates LDA and DME contributions (which have the same contact part as FAM)
         if (abs(debug)<=1) then
            call init_external_field(beta_type=beta_type, label=operator_name, k=operator_k, op=f,&
                    use_2bc=two_body_current_mode)
            do i=1,size(f%mat%elem)
               elem_tmp(i) = elem_tmp(i) + f%mat%elem(i)
            end do
 
            ! If mode needs no FAM calc, we are done: F = -sigma tau + sigma tau f(rho)
            if (f%use_2bc(2) /= 1 .and. f%use_2bc(2) /= 5) then
               f%mat%elem = elem_tmp
               return
            end if
         end if

         ! Finally, read in the Yukawa part or calculate and store it
         f%mat%elem = 0
         call fam_io(-1,iexit)
         if (iexit > 0) then 
            write(st,'(A)') "  * Calculating 2BC matrix elements..."
            call writelog(st)

            time1 = get_timer()
            call effective_2bc_extfield(f)
            time2 = get_timer()

            write(st,'(A, ES12.3)') "  * Two-body current wall time (min): ", (time2-time1)/60.0_dp
            call writelog(st)
            call writelog("")

            ! Use only the direct part of FAM with DME
            if (f%use_2bc(2) == 5) then
                call calc_totgamdel(f%gam, f%mat%elem, 'd')
            end if

            call fam_io(-2,iexit) ! deallocates f%gam/del
         end if

         ! Add 2-body + 1-body and/or contact part
         do i=1,size(f%mat%elem)
            f%mat%elem(i) = f%mat%elem(i) + elem_tmp(i)
         end do

         ! Pairing only - NB I have not checked the contact parts contribution to pairing is correct
         if (f%use_2bc(3) == 2) then
            f%mat12%elem = f%mat%elem
            call deallocate_blockmatrix(f%mat)
         end if

      end if

   end subroutine setup_extfield


   ! NB: Fsp, dRqp, Fqp, imag_eqrpa are all globals in this module
   subroutine transition_density
      use hfb_solution
      use pnfam_constants, only : pi, IT_NEUTRON, IT_PROTON
      use type_blockmatrix
      use type_extfield_2bc, only : init_extfield_2bc_type, caux, cd
      use pnfam_extfield, only : tbc_nmlda_da1, dme_exc
      implicit none
      type(bigblockmatrix) :: dRsp_QRPA, dRqp_QRPA, dR_tmp
      real(dp) :: prefactor, N, Rij, Fij, Sij, amplitude, integral_td1, integral_td2
 
      integer :: ipt, ibx1, ibx2, nd1, nd2, i2, ix2, i1, ix1
      integer :: xl1, xs1, xl2, xs2

      real(dp) :: td(nghl), wf_1(nghl), wf_2(nghl)
      real(dp) :: rhop(nghl), rhon(nghl), rho_fac(nghl)
      character(len=1200) :: denfile, stde

      ! Redundant, but init the couplings
      call init_extfield_2bc_type(.false.)

      ! Transition ampltiude from QP basis solution:
      !   * -1/pi Im(S) = 1/(hw*pi) |<i|F|0>|^2
      !   * |<i|F|0>|^2 = -hw Im(S)
      call contract_bbm(Fqp, dRqp_im, amplitude)
      amplitude = -imag_eqrpa*amplitude ! +/- |<i| F |0>|^2
      amplitude = sqrt(abs(amplitude))  !     |<i| F |0>|

      ! Assume QRPA matrices are real, then:
      !   * Im(X_FAM) = -1/hw (X_QRPA) <i|F|0>
      !   * X_QRPA = (-hw/<i|F|0>) Im(X_FAM)
      call copy_bbm(dRqp_QRPA, dRqp_im) ! dRqp_QRPA = dRqp_im (see copy_bbm for why not "=")
      dRqp_QRPA%imag=.false. ! Declare real
      prefactor = -imag_eqrpa/amplitude
      call scalar_mult_bbm(dRqp_QRPA, prefactor)

      ! Check normalization to see if we are at an eigenvalue:
      !   * N = \sum_{uv} |X_{uv}|^2 - |Y_{uv}|^2
      call copy_bbm(dR_tmp, dRqp_QRPA) ! dR_tmp = dRqp_QRPA
      if (allocated(dR_tmp%m21%elem)) dR_tmp%m21%elem = -dR_tmp%m21%elem
      if (allocated(dR_tmp%m22%elem)) dR_tmp%m22%elem = -dR_tmp%m22%elem
      call contract_bbm(dRqp_QRPA, dR_tmp, N)

      ! Convert X/Y_QRPA to single particle basis
      call copy_bbm(dRsp_QRPA, dRsp_im) ! dRsp_QRPA = dRsp_im
      dRsp_QRPA%imag = .false.! Declare real
      if (bminus) then
         call triprod_bbm('n',Wp,'n',dRqp_QRPA,'t',Wn, dRsp_QRPA)
      else
         call triprod_bbm('n',Wn,'n',dRqp_QRPA,'t',Wp, dRsp_QRPA)
      end if

      ! Calculate density of the transition amplitude
      !   * Since Fsp = (f11, 0; 0, 0), only need to contract with drho_QRPA
      !   * TD = \sum_{ij} rho_{ij} f11_{ij} phi_{i} phi_{j} ??? Also what about spin and lambda?
      ipt = 0; td  = 0
      do ibx1 = 1, nb
         ibx2 = dRsp_QRPA%m11%ir2c(ibx1)
         if (ibx2 == 0) cycle
         nd1 = db(ibx1) ; nd2 = db(ibx2)
         if (nd1 == 0 .or. nd2 == 0) cycle

         do i2=1, nd2
            ! State |2>
            ix2 = i2 + isstart(ibx2) - 1
            wf_2 = wf(:,ix2)
            xl2 = nl(ix2);  xs2 = ns(ix2)

            do i1=1, nd1
               ! State <1|
               ix1 = i1 + isstart(ibx1) - 1
               wf_1 = wf(:,ix1)
               xl1 = nl(ix1);  xs1 = ns(ix1)

               ! Matrix element
               ipt = ipt + 1

               ! Technically non-zero in td, but won't contribute to amplitude
               if (xl1 /= xl2) cycle

               ! Spin matrix element
               select case(f%K)
                  case (0)
                     if (xs1 == xs2) then
                        Sij = xs1
                     end if
                  case (1,-1)
                     if (xs1 == xs2 + 2*f%K) then
                        Sij = -f%K*sqrt(2.0_dp)
                     end if
               end select

               Rij = dRsp_QRPA%m11%elem(ipt)
               td = td + Rij*wf_1*wf_2*Sij

            end do ! i1
         end do ! i2
      end do ! ibx

      ! Get densities and NM density function
      call hfb_density_coord(IT_NEUTRON,rhon)
      call hfb_density_coord(IT_PROTON, rhop)

      rho_fac = (caux*2*cd)*(rhon + rhop)
      if (f%use_2bc(2) == 2) then
         ! Symmetric nuclear matter calculation at P=0
         rho_fac = rho_fac + tbc_nmlda_da1(rhon, rhop, 0.0_dp, .true.)
      else if (f%use_2bc(2) == 3) then
         ! Asymmetric nuclear matter calculation at P=0
         rho_fac = rho_fac + tbc_nmlda_da1(rhon, rhop, 0.0_dp, .false.)
      else if (f%use_2bc(2) == 4) then
         ! DME exchange term
         rho_fac = rho_fac + dme_exc()
      else
         rho_fac = 0
      end if

      integral_td1 = sum(td)
      integral_td2 = sum(td*rho_fac)

      ! Strip integration weights
      td = td*wdcori

      write(denfile,'(a,".den")') trim(fam_output_filename)
      call openlog(denfile, .false.)
      call writelog(" QRPA Transition Amplitude Densities")

      write(stde, '(1x,a20,es20.10e3)') "Re(EQRPA)         = ", real_eqrpa ; call writelog(stde)
      write(stde, '(1x,a20,es20.10e3)') "Im(EQRPA)         = ", imag_eqrpa ; call writelog(stde)
      write(stde, '(1x,a20,es20.10e3)') "QRPA Norm         = ", N ; call writelog(stde)
      write(stde, '(1x,a20,es20.10e3)') "Amplitude         = ", amplitude ; call writelog(stde)
      write(stde, '(1x,a20,es20.10e3)') "Int. trans        = ", integral_td1 ; call writelog(stde)
      write(stde, '(1x,a20,es20.10e3)') "Int. trans*f(rho) = ", integral_td2 ; call writelog(stde)

      write(stde, '(6a20)') "R", "Z", "rhon", "rhop", "trans", "frho_2bc" ; call writelog(stde, .true.)
      do i1 = 1, nghl
         write(stde, '(6ES20.10E3)') 1/y(i1), z(i1), rhon(i1), rhop(i1), td(i1), rho_fac(i1)
         call writelog(stde, .true.)
      end do
      call closelog

   end subroutine transition_density

   !---------------------------------------------------------------------------
   ! Component-wise bigblockmatrix copy, "dest = src".
   !
   ! WORKAROUND: with this project's production flags (-O3 -xCORE-AVX2
   ! -axCORE-AVX512 -fma, ifx 2023.2.0), a plain whole-structure assignment
   ! "dest = src" between two bigblockmatrix variables crashes at runtime with
   !     forrtl: severe (189): LHS and RHS of an assignment statement have
   !     incompatible types
   ! even though the two sides are the exact same type (confirmed with -g
   ! -traceback plus a debugger: the crash is inside the compiler-generated
   ! deep-copy routine, not in anything under our control, and splitting the
   ! assignment component-by-component -- exactly what this routine does --
   ! makes it disappear). This reproduces for both a local-variable
   ! destination and an intent(out) dummy-argument destination, and for both
   ! partially- and fully-allocated source quadrants, so it looks like an ifx
   ! code-generation bug tied to this particular (large, nested-allocatable)
   ! derived type rather than anything about intent/allocation state. Use this
   ! helper instead of "=" for bigblockmatrix values anywhere new code needs
   ! that copy; do not simplify back to "dest = src".
   !---------------------------------------------------------------------------
   subroutine copy_bbm(dest, src)
      implicit none
      type(bigblockmatrix), intent(out) :: dest
      type(bigblockmatrix), intent(in)  :: src

      dest%m11 = src%m11
      dest%m12 = src%m12
      dest%m21 = src%m21
      dest%m22 = src%m22
      dest%t11 = src%t11; dest%t12 = src%t12
      dest%t21 = src%t21; dest%t22 = src%t22
      dest%s11 = src%s11; dest%s12 = src%s12
      dest%s21 = src%s21; dest%s22 = src%s22
      dest%imag = src%imag

   end subroutine copy_bbm


   !---------------------------------------------------------------------------
   ! Shared step 1-3 of the transition-density routines: extract the real QRPA
   ! amplitudes (X^m,Y^m) for the state sitting at the current pole (real_eqrpa
   ! + i*imag_eqrpa) from whichever external field is currently active, and
   ! convert them to the s.p.-basis density-matrix perturbation drho_pn (Eqs.
   ! B1-B4,B8,B9 of Ney, Engel, Schunck, PRC 105, 034349 (2022)). This part is
   ! completely operator-agnostic -- it only depends on the QRPA structure of
   ! the response through Fqp/dRqp_im/dRsp_im/Wp/Wn, not on which field they
   ! were built from. Used by transition_density_0m.
   !
   ! 'amplitude' is |<m|F|0>| for whichever field F is currently active; N is
   ! the QRPA-norm consistency check (should be close to 1 at a good pole).
   !---------------------------------------------------------------------------
   subroutine extract_qrpa_drho(dRsp_QRPA, N, amplitude)
      implicit none
      type(bigblockmatrix), intent(out) :: dRsp_QRPA
      real(dp), intent(out) :: N, amplitude
      type(bigblockmatrix) :: dRqp_QRPA, dR_tmp
      real(dp) :: prefactor

      call contract_bbm(Fqp, dRqp_im, amplitude)
      amplitude = -imag_eqrpa*amplitude ! +/- |<m|F|0>|^2
      amplitude = sqrt(abs(amplitude))  !     |<m|F|0>|

      call copy_bbm(dRqp_QRPA, dRqp_im) ! dRqp_QRPA = dRqp_im (see copy_bbm for why)
      dRqp_QRPA%imag=.false. ! Declare real
      prefactor = -imag_eqrpa/amplitude
      call scalar_mult_bbm(dRqp_QRPA, prefactor)

      call copy_bbm(dR_tmp, dRqp_QRPA) ! dR_tmp = dRqp_QRPA
      if (allocated(dR_tmp%m21%elem)) dR_tmp%m21%elem = -dR_tmp%m21%elem
      if (allocated(dR_tmp%m22%elem)) dR_tmp%m22%elem = -dR_tmp%m22%elem
      call contract_bbm(dRqp_QRPA, dR_tmp, N)

      call copy_bbm(dRsp_QRPA, dRsp_im) ! dRsp_QRPA = dRsp_im (sizes the block structure)
      dRsp_QRPA%imag = .false.! Declare real
      if (bminus) then
         call triprod_bbm('n',Wp,'n',dRqp_QRPA,'t',Wn, dRsp_QRPA)
      else
         call triprod_bbm('n',Wn,'n',dRqp_QRPA,'t',Wp, dRsp_QRPA)
      end if

   end subroutine extract_qrpa_drho


   !---------------------------------------------------------------------------
   ! Project a s.p.-basis density-matrix perturbation drho_pn (from
   ! extract_qrpa_drho) onto the coordinate grid with BOTH the RS0 and the PS0
   ! spatial kernels at once (the two J^pi=0- fields of Mustonen et al., PRC 90,
   ! 024308 (2014), Eqs. A2a,A2d). Kernels copied verbatim from the 'RS0' and
   ! 'PS0' cases of ext_field_operator (pnfam_extfield.f90). RS0's branches are
   !     dLambda=0,  d(2*Sigma)=0   ->  f_pn = -Sigma_1 * <1| z |2>
   !     dLambda=+1, d(2*Sigma)=-2  ->  f_pn = -<1| r |2>
   !     dLambda=-1, d(2*Sigma)=+2  ->  f_pn = -<1| r |2>
   ! PS0's are
   !     dLambda=0,  d(2*Sigma)=0   ->  f_pn = -Sigma_1 <1|d/dz|2>
   !     dLambda=-1, d(2*Sigma)=+2  ->  f_pn = -<1|d/dr|2> - Lambda_2 <1|1/r|2>
   !     dLambda=+1, d(2*Sigma)=-2  ->  f_pn = -<1|d/dr|2> + Lambda_2 <1|1/r|2>
   !
   ! If two-body currents are active for this run (f%use_2bc, set by
   ! set_use_2bc from two_body_current_mode), the same density-dependent
   ! correction factors ext_field_operator multiplies into the RS0/PS0
   ! matrix elements (correction_2bc_rsL, correction_2bc_axial_charge) are
   ! applied here too, using the *same* f%use_2bc dispatch it uses (RS0's
   ! and PS0's crossterm fields are built with the identical raw
   ! two_body_current_mode as the main field -- see setup_crossterms -- so
   ! f%use_2bc is valid for both kernels regardless of which one is the
   ! active operator for this run). Without this, a 2bc-enabled run's
   ! density would silently miss the operator-level 2-body-current piece
   ! (confirmed empirically: bare-kernel integral undershot the true PS0
   ! amplitude by ~25% for a 2bc test case).
   !
   ! Also returns rho_rs0i, the density for RS0I = Ors0 * I(1,1,1,1;r), the
   ! generalized-Fermi-function-weighted variant of RS0 used by the "new"
   ! (UseNewShapeFactors) J^pi=0- shape factor in pynfam/strength/shape_factor.py
   ! (calcSfByJ, else branch). I(1,1,1,1;r) is the closed-form uniform-charge-
   ! sphere approximation (Behrens & Buhring) that ext_field_operator computes
   ! in calculate_I1111_r() -- reproduced verbatim here (it depends only on the
   ! grid and the mass number, not on any FAM solution, so there's nothing to
   ! share/reuse from that module). RS0I picks up the *same* 2bc correction as
   ! RS0 (ext_field_operator gates it on label being RS0/RS1/RS2/RS0I/RS1I
   ! together), so rho_rs0i is just rho_rs0's kernel with wf_1 additionally
   ! weighted by I(r).
   !---------------------------------------------------------------------------
   subroutine project_rs0_ps0_density(dRsp_QRPA, td_rs0, td_ps0, td_rs0i)
      use hfb_solution, only : nl, ns, wf, wfdr, wfdz, y, z, nghl
      use type_blockmatrix
      use type_extfield_2bc, only : init_extfield_2bc_type
      use pnfam_extfield, only : forbidden_2bc_rsL, forbidden_2bc_PS0, dme_axial_2
      implicit none
      type(bigblockmatrix), intent(in) :: dRsp_QRPA
      real(dp), intent(out) :: td_rs0(nghl), td_ps0(nghl), td_rs0i(nghl)

      integer :: ipt, ibx1, ibx2, nd1, nd2, i2, ix2, i1, ix1
      integer :: xl1, xs1, xl2, xs2, Amass
      real(dp) :: wf_1(nghl), wf_2(nghl), dr_wf_2(nghl), dz_wf_2(nghl), r(nghl), radius(nghl)
      real(dp) :: wf_1_rs0(nghl), wf_1_ps0(nghl), wf_1_rs0i(nghl)
      real(dp) :: correction_2bc_rsL(nghl), correction_2bc_axial_charge(nghl)
      real(dp) :: I_func(nghl), R_nuc
      real(dp) :: Rij
      logical :: rsL_snm

      r(:) = 1.0_dp/y(:)
      td_rs0 = 0; td_ps0 = 0; td_rs0i = 0

      ! I(1,1,1,1;r), Behrens-Buhring uniform-sphere approximation -- verbatim
      ! copy of calculate_I1111_r() in pnfam_extfield.f90
      Amass = hfb_npr(3)
      R_nuc = 1.2_dp * real(Amass,dp)**(1.0_dp/3.0_dp)
      radius = sqrt(r**2 + z**2)
      where (radius <= R_nuc)
         I_func = 1.5_dp * (1.0_dp - 0.2_dp*radius**2/R_nuc**2)
      elsewhere
         I_func = 1.5_dp * (R_nuc/radius - 0.2_dp*R_nuc**3/radius**3)
      end where

      ! Precompute the 2bc correction factors (grid functions only, same for
      ! every basis-state pair), exactly as ext_field_operator would for RS0
      ! and PS0 with this run's f%use_2bc.
      correction_2bc_rsL = 0
      if (f%use_2bc(6) > 0) then
         call init_extfield_2bc_type(.false.)
         rsL_snm = (f%use_2bc(6) /= 2)
         correction_2bc_rsL = forbidden_2bc_rsL(rsL_snm, f%use_2bc(2))
      end if

      correction_2bc_axial_charge = 0
      if (f%use_2bc(1) > 0 .and. f%use_2bc(5) > 0) then
         call init_extfield_2bc_type(.false.)
         if (f%use_2bc(5) == 1) then
            correction_2bc_axial_charge = forbidden_2bc_PS0()
         else if (f%use_2bc(5) == 2) then
            correction_2bc_axial_charge = dme_axial_2()
         end if
      end if

      ipt = 0
      do ibx1 = 1, nb
         ibx2 = dRsp_QRPA%m11%ir2c(ibx1)
         if (ibx2 == 0) cycle
         nd1 = db(ibx1) ; nd2 = db(ibx2)
         if (nd1 == 0 .or. nd2 == 0) cycle

         do i2=1, nd2
            ! State |2> = |n>
            ix2 = i2 + isstart(ibx2) - 1
            wf_2 = wf(:,ix2)
            dr_wf_2 = wfdr(:,ix2)
            dz_wf_2 = wfdz(:,ix2)
            xl2 = nl(ix2);  xs2 = ns(ix2)

            do i1=1, nd1
               ! State <1| = <p|
               ix1 = i1 + isstart(ibx1) - 1
               wf_1 = wf(:,ix1)
               xl1 = nl(ix1);  xs1 = ns(ix1)

               ! wf_1, weighted by each operator's 2bc correction (1-body
               ! only when the correction is off, matching ext_field_operator)
               if (f%use_2bc(6) > 0) then
                  wf_1_rs0 = wf_1 * (1.0_dp - correction_2bc_rsL)
               else
                  wf_1_rs0 = wf_1
               end if
               wf_1_rs0i = wf_1_rs0 * I_func
               if (f%use_2bc(1) == 1 .and. f%use_2bc(5) > 0) then
                  wf_1_ps0 = wf_1 * (1.0_dp + correction_2bc_axial_charge)
               else if (f%use_2bc(1) == 2 .and. f%use_2bc(5) > 0) then
                  wf_1_ps0 = wf_1 * correction_2bc_axial_charge
               else
                  wf_1_ps0 = wf_1
               end if

               ipt = ipt + 1
               Rij = dRsp_QRPA%m11%elem(ipt)

               ! RS0 = [r Y_1 x sigma]_00
               if ((xl1 == xl2) .and. (xs1 == xs2)) then
                  td_rs0 = td_rs0 - Rij*xs1*wf_1_rs0*wf_2*z
                  td_rs0i = td_rs0i - Rij*xs1*wf_1_rs0i*wf_2*z
               else if ((xl1 == xl2 + 1) .and. (xs1 == xs2 - 2)) then
                  td_rs0 = td_rs0 - Rij*wf_1_rs0*wf_2*r
                  td_rs0i = td_rs0i - Rij*wf_1_rs0i*wf_2*r
               else if ((xl1 == xl2 - 1) .and. (xs1 == xs2 + 2)) then
                  td_rs0 = td_rs0 - Rij*wf_1_rs0*wf_2*r
                  td_rs0i = td_rs0i - Rij*wf_1_rs0i*wf_2*r
               end if

               ! PS0 = sigma.grad / i
               if ((xl1 == xl2) .and. (xs1 == xs2)) then
                  td_ps0 = td_ps0 - Rij*xs1*wf_1_ps0*dz_wf_2
               else if ((xl1 == xl2 - 1) .and. (xs1 == xs2 + 2)) then
                  td_ps0 = td_ps0 - Rij*(wf_1_ps0*dr_wf_2 + xl2*wf_1_ps0*y*wf_2)
               else if ((xl1 == xl2 + 1) .and. (xs1 == xs2 - 2)) then
                  td_ps0 = td_ps0 - Rij*(wf_1_ps0*dr_wf_2 - xl2*wf_1_ps0*y*wf_2)
               end if

            end do ! i1
         end do ! i2
      end do ! ibx

   end subroutine project_rs0_ps0_density

   !---------------------------------------------------------------------------
   ! Combined RS0/PS0/RS0I transition density for the J^pi=0- nonunique
   ! first-forbidden shape factor. R1(0-) is skipped entirely (it is
   ! generally irrelevant for this analysis); R2(0-) is computed TWICE, once
   ! with each of the two shape-factor conventions in
   ! pynfam/strength/shape_factor.py's calcSfByJ, so they can be compared
   ! side by side:
   !
   !   R2_old -- UseNewShapeFactors=False (Mustonen et al. Eq. A6b): uses
   !             only RS0 and PS0, with the Coulomb/energy dependence
   !             collapsed into the single constant X+ = W0/3 + alpha*Z/(2R)
   !             (Eq. A9), i.e. I(1,1,1,1;r) approximated by its value at the
   !             nuclear surface.
   !   R2_new -- UseNewShapeFactors=True: replaces X+ by the actual
   !             r-dependent generalized Fermi function I(1,1,1,1;r) (the
   !             Behrens-Buhring uniform-sphere form computed in
   !             project_rs0_ps0_density), which multiplies only the
   !             Coulomb/alpha*Z part of X+, not the W0/3 part. This splits
   !             every A^2 or A*(...) term that used to carry a single X+
   !             into an "A" piece (weighted by w0max/3, unchanged) and an
   !             "AI" piece (weighted by alphazed/3, using
   !             AI = <m|Ors0*I(r)|0>, i.e. the RS0I operator, in place of a
   !             second factor of A).
   !
   ! Concretely (alphazed = alpha*Z/R, no factor of 2 -- the 1/2 lives in how
   ! it recombines with w0max/3 below):
   !
   !   R2_old(0-)/gA^2 = (xp^2+1/9)*A^2 + B^2 + 2*xp*A*B,  xp = w0max/3 + alphazed/2
   !   R2_new(0-)/gA^2 =  B^2 + (1/9)*alphazed^2*AI^2 + (w0max^2/9)*A^2
   !                   +  (2/9)*alphazed*w0max*AI*A + (2/3)*alphazed*AI*B
   !                   +  (2/3)*w0max*A*B + (1/9)*A^2
   !
   ! except this is schematic: A, B, AI here are the *physical* operator
   ! matrix elements of Mustonen Eqs. A2a/A2d, which carry dimensionful
   ! prefactors (lambda_e, m_e*c^2/M_n, Theta_K) that pnFAM's bare geometric
   ! kernels (what extract_qrpa_drho/project_rs0_ps0_density actually
   ! return) do not include -- shape_factor.py applies them via its `pre`
   ! dict (pre_RS0=1/lambda_e, pre_PS0=-lambda_e*m_e*c^2/M_n, pre_RS0I=
   ! pre_RS0). Because pre_PS0/pre_RS0 ~ -81, this rescaling is *not* a
   ! minor correction -- it changes which term dominates R2. Every A, B, AI
   ! in the formulas actually implemented below is therefore the bare pnFAM
   ! amplitude multiplied by its own `pre_*` factor at the point of use, to
   ! match shape_factor.py exactly. ec=+1 for both beta-/beta+ (shape_factor.
   ! py's ec=-1 branch is electron capture, not handled here or by bminus).
   !
   ! As before, each bilinear X*Y is spread into a density as the symmetric
   ! combination (X*rho_Y + Y*rho_X)/2 (reducing to X*rho_X for a self term
   ! X=Y), so integrating rho_R2_old/rho_R2_new over the grid reproduces
   ! R2_bare_old/R2_bare_new exactly; rho_rs0i is included as its own column
   ! for the same reason rho_rs0/rho_ps0 were before. See the "R2_bare vs
   ! td_R2" note in pnfam_changes_summary.md for how the two relate.
   !
   ! W0 (Eq. 26/33) and the nuclear radius R follow the same derivation used
   ! throughout this file; gA^2 is still left out (a downstream constant).
   !---------------------------------------------------------------------------
   subroutine transition_density_0m
      use hfb_solution, only : nghl, y, z, wdcori
      use pnfam_constants, only : alpha, mec2, r0, hbar_mec, IT_NEUTRON, IT_PROTON
      use hfb_solution, only : hfb_lambda
      implicit none
      type(bigblockmatrix) :: dRsp_QRPA
      real(dp) :: N, amplitude
      real(dp) :: td_rs0(nghl), td_ps0(nghl), td_rs0i(nghl), td_R2_old(nghl), td_R2_new(nghl)
      real(dp) :: td_rs0xps0(nghl), td_rs0ixps0(nghl)
      real(dp) :: A_rs0, B_ps0, AI_rs0i, R2_bare_old, R2_bare_new
      real(dp) :: omega_max, W0, Rnuc, alphazed, ec, xp
      real(dp) :: pre_rs0, pre_ps0, pre_rs0i
      integer :: Zdaughter, Amass, i1
      real(dp), parameter :: dMnH = 0.78227_dp   ! neutron - H atom mass difference [MeV]
      real(dp), parameter :: Mn_mass = 939.0_dp  ! nucleon mass [MeV], matches pynfam config.py MN
      character(len=1200) :: denfile, stde

      if (trim(f%label) /= 'RS0' .and. trim(f%label) /= 'PS0') then
         call writelog(" transition_density_0m: the external field for this run is neither"// &
            " RS0 nor PS0 -- skipping (rerun pnFAM at a J^pi=0- resonance with either).")
         return
      end if

      call extract_qrpa_drho(dRsp_QRPA, N, amplitude)
      call project_rs0_ps0_density(dRsp_QRPA, td_rs0, td_ps0, td_rs0i)

      A_rs0   = sum(td_rs0)  ! <m|Ors0|0>,        bare (no pre_rs0)
      B_ps0   = sum(td_ps0)  ! <m|Ops0|0>,        bare (no pre_ps0)
      AI_rs0i = sum(td_rs0i) ! <m|Ors0*I(r)|0>,   bare (no pre_rs0i)

      ! Beta-decay kinematics for this state (Eqs. 26,33)
      Amass = hfb_npr(3)
      if (bminus) then
         Zdaughter = hfb_npr(IT_PROTON) + 1
         omega_max = hfb_lambda(IT_NEUTRON) - hfb_lambda(IT_PROTON) + dMnH
      else
         Zdaughter = hfb_npr(IT_PROTON) - 1
         omega_max = hfb_lambda(IT_PROTON) - hfb_lambda(IT_NEUTRON) - dMnH
      end if
      W0       = 1.0_dp + (omega_max - real_eqrpa)/mec2
      Rnuc     = r0 * Amass**(1.0_dp/3.0_dp) / hbar_mec
      alphazed = alpha*real(Zdaughter,dp)/Rnuc
      ec       = 1.0_dp

      ! Operator prefactors (pynfam/strength/shape_factor.py prepStrengths 'pre' dict, K=0)
      pre_rs0  = 1.0_dp/hbar_mec
      pre_ps0  = -hbar_mec*mec2/Mn_mass
      pre_rs0i = pre_rs0

      ! R1(0-) is skipped here -- it is generally irrelevant for this analysis.

      ! R2(0-)/gA^2, OLD shape factor -- shape_factor.py calcSfByJ, J0_R2,
      ! UseNewShapeFactors=False. Only RS0 and PS0 enter (no RS0I); the
      ! Coulomb correction that the new formula carries via RS0I/alphazed is
      ! instead folded into xp here.
      xp = W0/3.0_dp + ec*0.5_dp*alphazed

      R2_bare_old = (xp**2 + 1.0_dp/9.0_dp)*pre_rs0**2*A_rs0**2 &
                  + pre_ps0**2*B_ps0**2 &
                  + 2.0_dp*xp*pre_rs0*pre_ps0*A_rs0*B_ps0

      td_R2_old = (xp**2 + 1.0_dp/9.0_dp)*pre_rs0**2*A_rs0*td_rs0 &
                + pre_ps0**2*B_ps0*td_ps0 &
                + 2.0_dp*xp*pre_rs0*pre_ps0*0.5_dp*(A_rs0*td_ps0 + B_ps0*td_rs0)

      ! R2(0-)/gA^2, NEW shape factor -- shape_factor.py calcSfByJ, J0_R2,
      ! UseNewShapeFactors=True
      R2_bare_new = pre_ps0**2*B_ps0**2 &
              + (1.0_dp/9.0_dp)*alphazed**2*pre_rs0i**2*AI_rs0i**2 &
              + (W0**2/9.0_dp)*pre_rs0**2*A_rs0**2 &
              + (2.0_dp/9.0_dp)*alphazed*W0*pre_rs0i*pre_rs0*AI_rs0i*A_rs0 &
              + (2.0_dp/3.0_dp)*alphazed*pre_rs0i*pre_ps0*AI_rs0i*B_ps0 &
              + (2.0_dp/3.0_dp)*W0*pre_rs0*pre_ps0*A_rs0*B_ps0 &
              + (1.0_dp/9.0_dp)*pre_rs0**2*A_rs0**2

      td_R2_new = pre_ps0**2*B_ps0*td_ps0 &
            + (1.0_dp/9.0_dp)*alphazed**2*pre_rs0i**2*AI_rs0i*td_rs0i &
            + (W0**2/9.0_dp)*pre_rs0**2*A_rs0*td_rs0 &
            + (2.0_dp/9.0_dp)*alphazed*W0*pre_rs0i*pre_rs0*0.5_dp*(AI_rs0i*td_rs0 + A_rs0*td_rs0i) &
            + (2.0_dp/3.0_dp)*alphazed*pre_rs0i*pre_ps0*0.5_dp*(AI_rs0i*td_ps0 + B_ps0*td_rs0i) &
            + (2.0_dp/3.0_dp)*W0*pre_rs0*pre_ps0*0.5_dp*(A_rs0*td_ps0 + B_ps0*td_rs0) &
            + (1.0_dp/9.0_dp)*pre_rs0**2*A_rs0*td_rs0

      ! Bare (no pre_* weighting) crossterm densities, matching the "A*B
      ! (RS0xPS0)"/"AI*B (RS0IxPS0)" scalars above: the usual symmetrized
      ! spread of a bilinear X*Y into a density, (X*rho_Y + Y*rho_X)/2, so
      ! that integrating reproduces A*B / AI*B exactly. There is no unique
      ! pointwise density for a crossterm (any density integrating to the
      ! same number would do) -- this is the same convention used to build
      ! rho_R1/rho_R2 above.
      td_rs0xps0  = 0.5_dp*(A_rs0*td_ps0   + B_ps0*td_rs0)
      td_rs0ixps0 = 0.5_dp*(AI_rs0i*td_ps0 + B_ps0*td_rs0i)

      ! Strip integration weights -> genuine densities on the grid
      td_rs0      = td_rs0     *wdcori
      td_ps0      = td_ps0     *wdcori
      td_rs0i     = td_rs0i    *wdcori
      td_R2_old   = td_R2_old  *wdcori
      td_R2_new   = td_R2_new  *wdcori
      td_rs0xps0  = td_rs0xps0 *wdcori
      td_rs0ixps0 = td_rs0ixps0*wdcori

      write(denfile,'(a,"_0m.den")') trim(fam_output_filename)
      call openlog(denfile, .false.)
      call writelog(" QRPA Transition Amplitude Density -- combined RS0+PS0+RS0I (J^pi = 0-) shape factor")
      call writelog(" R1(0-) is skipped (generally irrelevant); R2(0-) is given both ways:")
      call writelog(" R2_old = shape_factor.py calcSfByJ, J0_R2, UseNewShapeFactors=False (no RS0I)")
      call writelog(" R2_new = shape_factor.py calcSfByJ, J0_R2, UseNewShapeFactors=True  (with RS0I)")

      write(stde, '(1x,a20,es20.10e3)') "Re(EQRPA)      = ", real_eqrpa ; call writelog(stde)
      write(stde, '(1x,a20,es20.10e3)') "Im(EQRPA)      = ", imag_eqrpa ; call writelog(stde)
      write(stde, '(1x,a20,es20.10e3)') "QRPA Norm      = ", N ; call writelog(stde)
      write(stde, '(1x,a20,i0)')        "Z (daughter)   = ", Zdaughter ; call writelog(stde)
      write(stde, '(1x,a20,i0)')        "A (mass)       = ", Amass ; call writelog(stde)
      write(stde, '(1x,a20,es20.10e3)') "W0             = ", W0 ; call writelog(stde)
      write(stde, '(1x,a20,es20.10e3)') "alpha*Z/R      = ", alphazed ; call writelog(stde)
      write(stde, '(1x,a20,es20.10e3)') "A  = <m|Ors0|0>        = ", A_rs0   ; call writelog(stde)
      write(stde, '(1x,a20,es20.10e3)') "B  = <m|Ops0|0>        = ", B_ps0   ; call writelog(stde)
      write(stde, '(1x,a20,es20.10e3)') "AI = <m|Ors0*I(r)|0>   = ", AI_rs0i ; call writelog(stde)
      write(stde, '(1x,a20,es20.10e3)') "A*B (RS0xPS0)          = ", A_rs0*B_ps0     ; call writelog(stde)
      write(stde, '(1x,a20,es20.10e3)') "AI*A (RS0IxRS0)        = ", AI_rs0i*A_rs0   ; call writelog(stde)
      write(stde, '(1x,a20,es20.10e3)') "AI*B (RS0IxPS0)        = ", AI_rs0i*B_ps0   ; call writelog(stde)
      write(stde, '(1x,a20,es20.10e3)') "R2_old(0-)/gA^2        = ", R2_bare_old ; call writelog(stde)
      write(stde, '(1x,a20,es20.10e3)') "R2_new(0-)/gA^2        = ", R2_bare_new ; call writelog(stde)

      ! 'weight' = 1/wdcori is the genuine HFB quadrature weight (incl.
      ! Jacobian) at each grid point, i.e. sum(rho_X*weight) over the grid
      ! reproduces the scalar X_bare exactly -- needed downstream to bin
      ! these point densities into proper volume integrals (e.g. a radial
      ! shell sum ~4*pi*r^2*rho(r)*dr), since the rho_* columns above have
      ! already had this weight divided out (see the "Strip integration
      ! weights" comment above).
      write(stde, '(10a16)') "R", "Z", "rho_rs0", "rho_ps0", "rho_rs0i", "rho_R2_old", "rho_R2_new", &
         "rho_rs0xps0", "rho_rs0ixps0", "weight"
      call writelog(stde, .true.)
      do i1 = 1, nghl
         write(stde, '(10ES16.7E3)') 1/y(i1), z(i1), td_rs0(i1), td_ps0(i1), td_rs0i(i1), td_R2_old(i1), td_R2_new(i1), &
            td_rs0xps0(i1), td_rs0ixps0(i1), 1.0_dp/wdcori(i1)
         call writelog(stde, .true.)
      end do
      call closelog

   end subroutine transition_density_0m


   !---------------------------------------------------------------------------
   ! Spin (z-component) and number density of the odd (blocked) nucleon.
   !
   ! Blocking quasiparticle mu changes the ph density by
   !     drho_ab = U_a,mu U_b,mu - V_a,mu V_b,mu           (real U,V)
   ! whose time-odd part is the odd-nucleon spin density. Within one Omega
   ! block the sigma_z contraction reduces to a per-spin-channel sum, since
   ! each spin channel carries a single Lambda:
   !     s_z(r) = [gU_up^2 - gV_up^2] - [gU_dn^2 - gV_dn^2]
   !     n (r)  = [gU_up^2 - gV_up^2] + [gU_dn^2 - gV_dn^2]
   ! with gX_s(r) = sum_{a in block, spin s} X_a,mu wf_a(r).
   !
   ! Writes one file per blocked isospin (r, z, n_odd, sz_odd). Only the axial
   ! s_z is computed; transverse components couple different Lambda and are not
   ! axially symmetric. Depends only on the HFB solution, so it can be called
   ! any time after setup_pnfam.
   !
   ! The ground-state spin-current |J0|,|J1| (and tau0, d2rho0) are time-even and
   ! thus nonzero for even-even nuclei too; those are always computed. Only the
   ! blocked-nucleon local spin/number density (n_odd, sz_odd, sperp_odd, smag_odd)
   ! requires an odd nucleus -- for even-even it is written as zero (there is no
   ! blocked QP to source it) in a single "*_spin_core.dat" file, rather than
   ! skipping the whole routine.
   !---------------------------------------------------------------------------
   subroutine blocked_spin_density
      use pnfam_constants, only : IT_NEUTRON, IT_PROTON, IT_ISOSCALAR
      use hfb_solution, only : nghl, hfb_density_coord, hfb_spin_current_coord, &
                               hfb_divJ_coord, tau, d2rho
      implicit none

      integer :: it, ib, is
      real(dp) :: intN, intSz, intSm
      real(dp), dimension(nghl) :: n_dens, sz_dens, sp_dens, sm_dens, rho_scalar
      real(dp), dimension(nghl) :: jrn, jzn, jrp, jzp, j0mag, j1mag, tau0, d2rho0
      real(dp), dimension(nghl) :: divjn, divjp, divj0, divj1
      character(len=4)   :: qlab

      ! Isoscalar scalar density (rho_n + rho_p), physical per-volume units,
      ! for the s/rho size comparison against the standard GT term.
      call hfb_density_coord(IT_ISOSCALAR, rho_scalar)

      ! Ground-state vector spin-current (spin-orbit) density: isoscalar |J0| and
      ! isovector |J1|, the first-order (r x J) DME term (time-even -> present in
      ! the even-even core). Also isoscalar kinetic tau0 and Laplacian d2rho0, for
      ! the full GT/J DME coefficient functions (U, I integrals). Always computed,
      ! regardless of whether a QP is blocked.
      call hfb_spin_current_coord(IT_NEUTRON, jrn, jzn)
      call hfb_spin_current_coord(IT_PROTON,  jrp, jzp)
      j0mag  = sqrt((jrn+jrp)**2 + (jzn+jzp)**2)
      j1mag  = sqrt((jrn-jrp)**2 + (jzn-jzp)**2)
      tau0   = tau(:,IT_NEUTRON)   + tau(:,IT_PROTON)
      d2rho0 = d2rho(:,IT_NEUTRON) + d2rho(:,IT_PROTON)

      ! div(J): the crdj (NOT cdrho) coupling partner rho*div(J). Isoscalar
      ! div(J0) and isovector div(J1), same time-even/always-present status as J.
      call hfb_divJ_coord(IT_NEUTRON, divjn)
      call hfb_divJ_coord(IT_PROTON,  divjp)
      divj0 = divjn + divjp
      divj1 = divjn - divjp

      if (.not. hfb_blo_active) then
         call writelog(" blocked_spin_density: no blocked QP (even-even) --"// &
                       " writing spin-current only (core file).")
         n_dens = 0; sz_dens = 0; sp_dens = 0; sm_dens = 0
         intN = 0; intSz = 0; intSm = 0
         call write_spin_density_file("core", 0, 0, n_dens, sz_dens, sp_dens, sm_dens, &
                                      j0mag, j1mag, divj0, divj1, tau0, d2rho0, rho_scalar, &
                                      intN, intSz, intSm)
         return
      end if

      do it = IT_NEUTRON, IT_PROTON
         if (hfb_blo_qp(it) == 0) cycle    ! this isospin is not blocked

         ib = hfb_blo_ib(it)               ! HFBTHO block of the blocked level
         is = hfb_blo_is(it)               ! QP column within that block

         if (it == IT_NEUTRON) then
            call spin_density_of_qp(Un, Vn, ib, is, n_dens, sz_dens, sp_dens, sm_dens, &
                                    intN, intSz, intSm)
            qlab = "neu"
         else
            call spin_density_of_qp(Up, Vp, ib, is, n_dens, sz_dens, sp_dens, sm_dens, &
                                    intN, intSz, intSm)
            qlab = "pro"
         end if

         call write_spin_density_file(qlab, ib, is, n_dens, sz_dens, sp_dens, sm_dens, &
                                      j0mag, j1mag, divj0, divj1, tau0, d2rho0, rho_scalar, &
                                      intN, intSz, intSm)
      end do

   end subroutine blocked_spin_density

   !---------------------------------------------------------------------------
   ! Shared writer for blocked_spin_density: one "<base>_spin_<qlab>.dat" file
   ! with the 13-column (r,z,n_odd,sz_odd,sperp_odd,smag_odd,j0mag,j1mag,
   ! divj0,divj1,tau0,d2rho0,rho_scalar) layout. ib/is are 0 and the odd-nucleon
   ! columns are zero for the even-even "core" case.
   !---------------------------------------------------------------------------
   subroutine write_spin_density_file(qlab, ib, is, n_dens, sz_dens, sp_dens, sm_dens, &
                                      j0mag, j1mag, divj0, divj1, tau0, d2rho0, rho_scalar, &
                                      intN, intSz, intSm)
      use hfb_solution, only : nghl, y, z
      implicit none
      character(len=*), intent(in) :: qlab
      integer, intent(in) :: ib, is
      real(dp), dimension(nghl), intent(in) :: n_dens, sz_dens, sp_dens, sm_dens
      real(dp), dimension(nghl), intent(in) :: j0mag, j1mag, divj0, divj1
      real(dp), dimension(nghl), intent(in) :: tau0, d2rho0, rho_scalar
      real(dp), intent(in) :: intN, intSz, intSm

      integer :: i
      character(len=200) :: fname
      character(len=300) :: stde   ! 13 columns * a20/ES20 = 260 chars; needs > 260

      if (len_trim(fam_output_filename) > 0) then
         write(fname,'(a,"_spin_",a,".dat")') trim(fam_output_filename), trim(qlab)
      else
         write(fname,'("spin_",a,".dat")') trim(qlab)
      end if

      call openlog(fname, .false.)
      write(stde,'(a,a,a,i0,a,i0)') " Blocked-nucleon spin density (",trim(qlab), &
         ")  block=",ib,"  state=",is
      call writelog(stde)
      write(stde,'(1x,a24,es20.10e3)') "Int n   (blocking dN) = ", intN;  call writelog(stde)
      write(stde,'(1x,a24,es20.10e3)') "Int s_z (<sigma_z>)   = ", intSz; call writelog(stde)
      write(stde,'(1x,a24,es20.10e3)') "Int |s| (>= |<sz>|)   = ", intSm; call writelog(stde)
      write(stde,'(13a20)') "r", "z", "n_odd", "sz_odd", "sperp_odd", "smag_odd", &
                           "j0mag", "j1mag", "divj0", "divj1", "tau0", "d2rho0", "rho_scalar"
      call writelog(stde, .true.)
      do i = 1, nghl
         write(stde,'(13ES20.10E3)') 1/y(i), z(i), n_dens(i), sz_dens(i), sp_dens(i), &
                                    sm_dens(i), j0mag(i), j1mag(i), divj0(i), divj1(i), &
                                    tau0(i), d2rho0(i), rho_scalar(i)
         call writelog(stde, .true.)
      end do
      call closelog

   end subroutine write_spin_density_file

   !---------------------------------------------------------------------------
   ! Per-spin-channel contraction giving the number and axial spin density of a
   ! single blocked quasiparticle (block ib, QP column is). U and V share the
   ! same block-diagonal element layout (Vn%ir2m = Un%ir2m), so both columns are
   ! read from the same offset. wf carry sqrt of the integration weights, so the
   ! raw grid sums are the integrals; wdcori strips the weights for the local
   ! (per-volume) density.
   !---------------------------------------------------------------------------
   subroutine spin_density_of_qp(U, V, ib, is, n_dens, sz_dens, sp_dens, sm_dens, &
                                 intN, intSz, intSm)
      use hfb_solution, only : nghl, wf, ns, wdcori
      use type_blockmatrix, only : blockmatrix, db, isstart
      implicit none

      type(blockmatrix), intent(in) :: U, V
      integer, intent(in) :: ib, is
      real(dp), intent(out) :: n_dens(nghl), sz_dens(nghl), sp_dens(nghl), sm_dens(nghl)
      real(dp), intent(out) :: intN, intSz, intSm

      integer :: nd, base, a, ia
      real(dp) :: ua, va
      real(dp), dimension(nghl) :: gU_up, gU_dn, gV_up, gV_dn

      nd   = db(ib)
      base = U%ir2m(ib) + (is-1)*nd     ! column-major within the ib block

      gU_up = 0; gU_dn = 0; gV_up = 0; gV_dn = 0
      do a = 1, nd
         ia = isstart(ib) + a - 1
         ua = U%elem(base + a - 1)
         va = V%elem(base + a - 1)
         if (ns(ia) == 1) then          ! spin up (2*Sigma = +1)
            gU_up = gU_up + ua*wf(:,ia);  gV_up = gV_up + va*wf(:,ia)
         else                           ! spin down (2*Sigma = -1)
            gU_dn = gU_dn + ua*wf(:,ia);  gV_dn = gV_dn + va*wf(:,ia)
         end if
      end do

      ! s_z    : diagonal-in-spin (same Lambda) channel difference
      ! s_perp : |s_x + i s_y| = 2|C|, C the cross-spin-channel (Delta Lambda=1)
      !          product; the transverse vector winds as e^{i phi}, so |s_perp|
      !          is axially symmetric.
      ! s_mag  : sqrt(s_z^2 + s_perp^2), the full local vector magnitude |s(r)|.
      n_dens  = (gU_up**2 - gV_up**2) + (gU_dn**2 - gV_dn**2)
      sz_dens = (gU_up**2 - gV_up**2) - (gU_dn**2 - gV_dn**2)
      sp_dens = 2.0_dp*(gU_up*gU_dn - gV_up*gV_dn)
      sm_dens = sqrt(sz_dens**2 + sp_dens**2)

      ! Integrals (weights already folded into wf^2)
      intN  = sum(n_dens)
      intSz = sum(sz_dens)
      intSm = sum(sm_dens)          ! integral of |s(r)| dV  (>= |<sigma_z>|)

      ! Strip integration weights -> local per-volume densities
      n_dens  = n_dens  * wdcori
      sz_dens = sz_dens * wdcori
      sp_dens = abs(sp_dens) * wdcori    ! report transverse magnitude |s_perp|
      sm_dens = sm_dens * wdcori

   end subroutine spin_density_of_qp

end module pnfam_solver