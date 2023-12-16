!------------------------------------------------------------------------------
!> This modules contains the routines for computing two-body spatial matrix
!> elements of a Gaussian with up to 2 derivatives. The derivatives act either
!> on the Gaussian or on the wavefunctions.
!>
!> The routines are organized as follows, from roughly bottom to top of the module
!> 1. ho_derivative: computes coefficients and quantum numbers for 0,1,2 derivatives
!>    of a 1D HO state.
!> 2. MatrixElement_1d_dho: Wraps (1) to compute 2-body matrix elements of a 1D
!>    gaussian, with one or two derivatives acting on any of the 4 HO states.
!> 3. MatrixElement_1d: Wraps (2) to compute matrix 2-body matrix
!>    elements of a 1D gaussian with derivatives acting on the gaussian, or with
!>    momentum operator acting on wavefunctions to the left and right.
!> 4. MatrixElement_rad*: Wraps (3), using cartesian expansion to compute (r,phi)
!>    component of 2-body matrix elements with CARTESIAN derivatives in x/y.
!> 5. MatrixElement_ddGr*: Wraps (4) to compute (r,phi) component of 2-body matrix
!>    elements involving SPHERICAL components of derivatives.
!> 6. MatrixElement_ddGz: Wraps (3) to compute (z) component same as (5).
!>
!> @authors E.M. Ney, UNC Chapel Hill 2018-
!------------------------------------------------------------------------------
module pnfam_spatial_mtxels
   ! Variables from HFBTHO
   use hfb_solution, only : bp, bz, n00, nlx, nrx, nzx, ntx, nox, hdqp, nttx
   use hfbtho_gogny, only : Cp2c, T_z
   use hfbtho_gogny, only : MatrixElement_z, calculateTz, calculateCpolar2cartesian
   use hfbtho_gogny, only : MatrixElement_zZR, MatrixElement_rZR
   use hfbtho_gogny, only : ib_zrls, i_zrls
   use hfbtho, only : nzz, nzzx, nrr, nll, nss, noo

   implicit none
   private
   integer, parameter, private :: ipr = kind(1), pr  = kind(1.0d0)

   ! Public Functions
   public :: prep_gaussian, &
             MatrixElement_zZR,       MatrixElement_rZR, &
             MatrixElement_z,         MatrixElement_ddGz, &
             MatrixElement_radx,      MatrixElement_ddGr, &
             MatrixElement_rad_cmplx, MatrixElement_ddGr_cmplx, &
             bp, bz, &
             ib_zrls, i_zrls, nzz,nzzx,nrr,nll,nss,noo
             !calculate_Nblock, calculate_Zblock, L_block, zindex, rindex
   public :: n00, MatrixElement_1D ! For debugging

   ! Debug Gaussian
   !integer, parameter,       public :: nY0 = 1 !> n-Gaussians used in Y0 fit
   !real(pr), dimension(nY0), public :: Y0_mu = (/1.2_pr/)
   !real(pr), dimension(nY0), public :: Y0_pf = (/1.0_pr/)

   ! Gaussian fit to Yukawa
   integer, parameter,       public :: nY0 = 6 !> n-Gaussians used in Y0 fit
   real(pr), dimension(nY0), public :: Y0_mu = (/1._pr/sqrt(34.0_pr), 1._pr/sqrt(6.60_pr), 1._pr/sqrt(1.44_pr), &
                                                 1._pr/sqrt(0.38_pr), 1._pr/sqrt(0.15_pr), 1._pr/sqrt(0.13_pr) /)
   real(pr), dimension(nY0), public :: Y0_pf = (/6.79_pr, 2.41_pr, 0.786_pr, 0.241_pr, -0.062_pr, 0.078_pr /)
   integer, parameter,           public :: nspatial = 1 + 2*nY0 !> Number terms to sum over (delta + nG + nG)
   integer, dimension(nspatial), public :: spatial_mu !> Gaussian index for term (delta=0)
   integer, dimension(nspatial), public :: spatial_tm !> type of term (delta=0, ddG1=1, ddG2=2)

   ! Prefactors
   real(pr), parameter :: cr2 = sqrt(2.0_pr)
   real(pr), target :: prefd1z, prefd1r, prefd2z, prefd2r

   ! Arrays for speeding up performance
   real(pr), allocatable, dimension(:,:,:,:,:), private, target :: ME1Dz, ME1Dr
   real(pr), allocatable, dimension(:,:,:,:,:), private, target :: ME1Dr_10, ME1Dr_20, ME1Dr_11, ME1Dr_12, ME1Dr_01, ME1Dr_02

   !! Transforms for storing matrix elements in 1D arrays
   !integer, Public, allocatable :: Zblock(:,:,:,:) !< Indices corresponding to non-zero and different axial matrix elements
   !integer, Public, allocatable :: NBlock(:,:,:,:,:,:) !< Indices corresponding to non-zero and different radial matrix elements


contains

   !------------------------------------------------------------------------------
   !> Initialize arrays and coefficients for computing gaussian matrix elements
   !------------------------------------------------------------------------------
   subroutine prep_gaussian
     use hfbtho, only : finite_range, coulomb_gaussian
     use hfbtho_solver, only : base
     use pnfam_constants, only : Mpi, hbarc, pi

     implicit none

     integer(ipr) :: i, ii, start

     ! Initialize finite_range variables from hfbtho manually here
     if (.not. finite_range .and. .not. coulomb_gaussian) then
        ! This allocates finite range variables that would be allocated in call to thoalloc
        finite_range = .true.
        If(finite_range.or.coulomb_gaussian) Then
           If(Allocated(nrr)) Deallocate(nrr,nll,nss,noo,nzz,nzzx)
           Allocate(nrr(ntx),nll(ntx),nss(ntx),noo(ntx),nzz(ntx,ntx),nzzx(ntx))
        End If
        ! This allocates and populates some quantum number related arrays used by gogny.
        ! NB: It uses some redundant hfbtho_solution arrays, make sure they are still allocated!
        call base(.false.)
     end if

     ! Convert fit f(x)=exp(-x)/x=sum_i[ b_i * exp(-a_i*x^2) ] to Y0(x)=exp(-mpi*x)/(4pi*x)=sum_i[ b_i * exp(-x^2/mu_i^2) ]
     !   - Include the pion mass, Y0(x,mpi) ~ mpi*f(x*mpi)
     !   - Pion mass should be in inverse fermis for hbar=c=1: mpi[1/fm] = mpi[MeV]/hbarc[MeV*fm]
     !   - 1/mu^2 = a*mpi^2, so mu = [1/sqrt(a)]/mpi
     !   - Include extra factor of 1/4pi for definition of Y0
     do i=1,nY0
        Y0_mu(i) = Y0_mu(i)/(Mpi/hbarc)
        Y0_pf(i) = Y0_pf(i)*(Mpi/hbarc)*0.25_pr/pi
     end do

     ! Avoid doing some unnecessary operations under the loops (NB: This must be done before call to calculateME1D)
     prefd1z = 1.0_pr/(bz*cr2)
     prefd1r = 1.0_pr/(bp*cr2)
     prefd2z = 1.0_pr/(bz*bz*2.0_pr)
     prefd2r = 1.0_pr/(bp*bp*2.0_pr)

     ! Calculate the coeffs needed for computing gaussian elements
     ! NB: we need coefficients beyond the basis for derivatives of gaussians (up to d^2|n> ~ |n+/-2>)
     nzx = nzx + 2
     nrx = nrx + 1
     nlx = nlx + 2
     call calculateTz
     call calculateCpolar2cartesian
     nzx = nzx - 2
     nrx = nrx - 1
     nlx = nlx - 2

     ! Calculate 1D gaussian matrix elements (and with derivatives)
     call calculateME1D

     ! NB: I can't sum over Z and R parts of the various spatial terms separately!
     ! I want \sum_i Vi_r*Vi_z,  and (V1_r*V1_z) + (V2_r*V2_z)  /=  (V1_r + V2_r)*(V1_z + V2_z)
     ! So, I need an array with as many terms as summed over in the spatial pieces of the current
     ! I have 12 unique spatial vectors (J_xx) that at most have: delta + d1d2 G_{Y0} + d3d4 G_{Y0}
     ! This means for both V_r and V_z I need an array of size (1 + ny0 + ny0) for each (J_xx)
     start = 1
     do i=0,2
        select case(i)
           case(0) ! Delta
              spatial_tm(start) = 0 ! term type
              spatial_mu(start) = 0 ! gaussian width
              start = start + 1
           case(1,2) ! d1d2 Y0 and d3d4 Y0
              do ii=1, nY0
                 spatial_tm(start-1+ii) = i  ! 1st or 2nd term d1d2 Y0
                 spatial_mu(start-1+ii) = ii ! gaussian index
              end do
              start = start + nY0
        end select
     end do

   end subroutine prep_gaussian

   !------------------------------------------------------------------------------
   !> Populate 1D gaussian matrix elements for every gaussian of width mu, and for
   !> radial and z oscillator lengths, and with all combos of derivatives.
   !------------------------------------------------------------------------------
   subroutine calculateME1D
      implicit none

      integer(ipr) :: ni,nj,nk,nl,n,ig

      ! Bump up the quantum numbers to allow for 2nd derivs
      n = max(nzx+2, max(2*(nrx+1),nlx+2))
      if(allocated(ME1Dz)) deallocate(ME1Dz, ME1Dr)
      allocate(ME1Dz(1:nY0,0:n,0:n,0:n,0:n)); ME1Dz = 0
      allocate(ME1Dr(1:nY0,0:n,0:n,0:n,0:n)); ME1Dr = 0
      do ni = 0,n
         do nj = 0,n
            do nk = 0,n
               do nl = 0,n
                  ! By reflection symmetry these will be zero
                  if(mod(ni+nj+nk+nl,2).ne.0) cycle
                  do ig = 1, nY0
                     ME1Dz(ig,ni,nj,nk,nl) = MatrixElement_z(ni,nj,nk,nl,Y0_mu(ig),bz)*Y0_pf(ig)
                     ME1Dr(ig,ni,nj,nk,nl) = MatrixElement_z(ni,nj,nk,nl,Y0_mu(ig),bp)
                  enddo
               enddo
            enddo
         enddo
      enddo

      ! Use the usual quantum numbers, we're taking derivs here
      n = max(nzx, max(2*nrx,nlx))
      if(allocated(ME1Dr_10)) deallocate(ME1Dr_10, ME1Dr_20, ME1Dr_11, ME1Dr_12, ME1Dr_01, ME1Dr_02)
      allocate(ME1Dr_10(1:nY0,0:n,0:n,0:n,0:n)); ME1Dr_10 = 0
      allocate(ME1Dr_20(1:nY0,0:n,0:n,0:n,0:n)); ME1Dr_20 = 0
      allocate(ME1Dr_11(1:nY0,0:n,0:n,0:n,0:n)); ME1Dr_11 = 0
      allocate(ME1Dr_12(1:nY0,0:n,0:n,0:n,0:n)); ME1Dr_12 = 0
      allocate(ME1Dr_01(1:nY0,0:n,0:n,0:n,0:n)); ME1Dr_01 = 0
      allocate(ME1Dr_02(1:nY0,0:n,0:n,0:n,0:n)); ME1Dr_02 = 0
      do ni = 0,n
         do nj = 0,n
            do nk = 0,n
               do nl = 0,n
                  ! Run over all combos because derivs change n's
                  do ig = 1, nY0
                     ME1Dr_10(ig,ni,nj,nk,nl) = MatrixElement_1d(ni,nj,nk,nl,ig,.false.,1,0)
                     ME1Dr_20(ig,ni,nj,nk,nl) = MatrixElement_1d(ni,nj,nk,nl,ig,.false.,2,0)
                     ME1Dr_11(ig,ni,nj,nk,nl) = MatrixElement_1d(ni,nj,nk,nl,ig,.false.,1,1)
                     ME1Dr_12(ig,ni,nj,nk,nl) = MatrixElement_1d(ni,nj,nk,nl,ig,.false.,1,2)
                     ME1Dr_01(ig,ni,nj,nk,nl) = MatrixElement_1d(ni,nj,nk,nl,ig,.false.,0,1)
                     ME1Dr_02(ig,ni,nj,nk,nl) = MatrixElement_1d(ni,nj,nk,nl,ig,.false.,0,2)
                  enddo
               enddo
            enddo
         enddo
      enddo

   end subroutine

   !------------------------------------------------------------------------------
   !> Calculates z-component of two-body matrix elements of the interaction:
   !> \f[ V = \partial_i \partial_j e^{-(\vec{r}_1 - \vec{r}_2)^2]/\mu^2} \f]
   !> in the axial HO basis. The derivatives d_i are components of the gradient in the
   !> spherical basis, acting on the Gaussian. Optionally, d_i can be switched to 2iP_i,
   !> where \f[ P_i=-i/2(\overrightarrow{\grad}_i - \overleftarrow{\grad}_i) \f] acts on wavefunctions.
   !>
   !> The components of the gradient in the spherical basis are:
   !> \f[ \partial_{+1} = -1/\sqrt{2} (\partial_x + i \partial_y) \f]
   !> \f[ \partial_{0}  = \partial_z \f]
   !> \f[ \partial_{+1} = +1/\sqrt{2} (\partial_x - i \partial_y) \f]
   !------------------------------------------------------------------------------
   function MatrixElement_ddGz(ni,nj,nk,nl,mu,di,dj,dip) result(Vz)
     implicit none
     integer(ipr), intent(in) :: ni !< first  z-component quantum number
     integer(ipr), intent(in) :: nj !< second z-component quantum number
     integer(ipr), intent(in) :: nk !< third  z-component quantum number
     integer(ipr), intent(in) :: nl !< fourth z-component quantum number
     integer(ipr), intent(in) :: mu !< Gaussian width
     integer(ipr), intent(in) :: di,dj !< Component of gradient acting on V
     integer(ipr), intent(in) :: dip !< 0 = ignore, 1,2=change di to 2i pi(x_dip)
     real(pr) :: Vz !< Real matrix element (V_{ijkl})_{z}
     integer(ipr) :: p, check, i

     ! Check inputs
     do i=1,2
        if (i==1) then
           check = di
        else
           check = dj
        end if
        select case(check)
           case default
              write(*,'(A)') 'Invalid components of derivatives supplied.'
              stop
           case(-1,0,1)
              ! valid input
        end select
     end do

     Vz=0
     ! If di --> pi, reduce number d by 1, use dip=(0,1,2) for p part
     if (dip/=0) then
        p = 1
     else
        p = 0
     end if

     ! p0d0 = [1]_r * [pz dz]_z
     if (di==0 .and. dj==0) then
        Vz = MatrixElement_1d(ni,nj,nk,nl,mu,.true.,2-p,dip)
     ! p+d+ = [+1/2(pxdx-pydy+ipxdy+ipydx)]_r * [1]_z
     ! p-d- = [+1/2(pxdx-pydy-ipxdy-ipydx)]_r * [1]_z
     ! p+d- = [-1/2(pxdx+pydy+ipydx-ipxdy)]_r * [1]_z
     ! p-d+ = [-1/2(pxdx+pydy-ipydx+ipxdy)]_r * [1]_z
     else if (abs(di+dj)==2 .or. abs(di+dj)==0) then
        Vz = MatrixElement_1d(ni,nj,nk,nl,mu,.true.,0,0)
     ! p+d0 = [-1/sqrt(2)(px+ipy)]_r * [dz]_z
     ! p-d0 = [+1/sqrt(2)(px-ipy)]_r * [dz]_z
     ! p0d+ = [-1/sqrt(2)(dx+idy)]_r * [pz]_z
     ! p0d- = [+1/sqrt(2)(dx-idy)]_r * [pz]_z
     else if (abs(di+dj)==1) then
        ! ** z is only affected by d-->p if i==0.
        if (di==0) then
           Vz = MatrixElement_1d(ni,nj,nk,nl,mu,.true.,1-p,dip)
        else
           Vz = MatrixElement_1d(ni,nj,nk,nl,mu,.true.,1,0)
        end if
     end if

   end function MatrixElement_ddGz

   !------------------------------------------------------------------------------
   !> This function is similar to MarixElement_ddGrad_cmplx, but uses several properties
   !> to make the routine more efficient:
   !> 1) Matrix elements obey the Wigner-Eckart selection rule  -li-lj+lk+ll = i+j
   !>    where i,j are spherical components of the derivatives (0,+1,-1)
   !> 2) Matrix elements which obey the selection rule are real and can be expressed
   !>    entirely in terms of derivatives in x (vs dx +/- i dy).
   !> This function uses the simplified radial matrix element routine MatrixElement_radx.
   !------------------------------------------------------------------------------
   function MatrixElement_ddGr(ni,li,nj,lj,nk,lk,nl,ll,mu,di,dj,dip) result(Vr)
     implicit none
     integer(ipr),intent(in) :: ni !<first  radial principal quantum number
     integer(ipr),intent(in) :: nj !<second radial principal quantum number
     integer(ipr),intent(in) :: nk !<third  radial principal quantum number
     integer(ipr),intent(in) :: nl !<fourth radial principal quantum number
     integer(ipr),intent(in) :: li !<first  radial orbital quantum number
     integer(ipr),intent(in) :: lj !<second radial orbital quantum number
     integer(ipr),intent(in) :: lk !<third  radial orbital quantum number
     integer(ipr),intent(in) :: ll !<fourth radial orbital quantum number
     integer(ipr),intent(in) :: mu    !< Gaussian width
     integer(ipr),intent(in) :: di,dj !> Component of gradient acting on V
     integer(ipr),intent(in) :: dip   !> 0 = ignore, 1,2=change di to 2i pi(x_dip)
     real(pr) :: Vr !< Real matrix element (V_{ijkl})_{r,phi}

     integer(ipr) :: p, check, i

     ! Check inputs
     do i=1,2
        if (i==1) then
           check = di
        else
           check = dj
        end if
        select case(check)
           case default
              write(*,'(A)') 'Invalid components of derivatives supplied.'
              stop
           case(-1,0,1)
              ! valid input
        end select
     end do

     Vr = 0
     ! Wigner-Eckart selection rule (<i,j| didjG |k,l> ~ <i,j|di,dj,k,l> so i+j = k+l+di+dj)
     if (-li-lj+lk+ll+di+dj /= 0) return

     ! If di --> pi, reduce number d by 1, use dip=(0,1,2) for p part
     if (dip/=0) then
        p = 1
     else
        p = 0
     end if

     ! p0d0 = [1]_r * [pz dz]_z
     if (di==0 .and. dj==0) then
        Vr = MatrixElement_radx(ni,li,nj,lj,nk,lk,nl,ll,mu,0,0)

     ! p+d+ = [1/2(pxdx-pydy+ipxdy+ipydx)]_r * [1]_z
     ! p-d- = [1/2(pxdx-pydy-ipxdy-ipydx)]_r * [1]_z
     ! With selection rule equals 1/2(4*pxdx)
     else if (abs(di+dj)==2) then
        Vr = MatrixElement_radx(ni,li,nj,lj,nk,lk,nl,ll,mu,2-p, dip)*(2.0_pr)

     ! p+d0 = [-1/sqrt(2)(px+ipy)]_r * [dz]_z
     ! d+p0 = [-1/sqrt(2)(dx+idy)]_r * [pz]_z
     ! With selection rule equals -1/sqrt(2)(2*px)
     ! p-d0 = [+1/sqrt(2)(px-ipy)]_r * [dz]_z
     ! d-p0 = [+1/sqrt(2)(dx-idy)]_r * [pz]_z
     ! With selection rule equals +1/sqrt(2)(2*px)
     else if (abs(di+dj)==1) then
        ! ** r is only affected by d-->p if j==0.
        if (dj==0) then
           Vr = MatrixElement_radx(ni,li,nj,lj,nk,lk,nl,ll,mu,1-p, dip)*(-(di+dj)*cr2)
        else
           Vr = MatrixElement_radx(ni,li,nj,lj,nk,lk,nl,ll,mu,1, 0)*(-(di+dj)*cr2)
        end if

     ! p+d- = [-1/2(pxdx+pydy+ipydx-ipxdy)]_r * [1]_z
     ! p-d+ = [-1/2(pxdx+pydy-ipydx+ipxdy)]_r * [1]_z
     ! With selection rule equals -1/2(2*dx^2) or with p, -1/2(2pxdx - 2Im(pydx))
     else if (di+dj==0) then
        Vr = MatrixElement_radx(ni,li,nj,lj,nk,lk,nl,ll,mu,2-p, dip)*(-1.0_pr)
        if (dip/=0) then
           Vr = Vr + real(aimag(MatrixElement_rad_cmplx(ni,li,nj,lj,nk,lk,nl,ll,mu,1,0,0,dip)),pr)
        end if

     end if

   end function MatrixElement_ddGr

   !------------------------------------------------------------------------------
   !> Similar to MatrixElement_rad_cmplx, but restricted to dy,py=0. Results are consequently
   !> all real. This routine is more efficient than its complex counterpart because it skips over
   !> loops than introduce an imaginary component, which we know will be zero.
   !>
   !> ** - This function is an on-the-fly implementation of "HFBTHO_GOGNY.radial_matrix_elements"
   !>      with the additional ability to include derivatives.
   !>    - Here gogny_hyper==1 is required to treat cartesian dervatives.
   !>    - This routine is more general than the one in HFBTHO_GOGNY since it calculates matrix elements
   !>      when -l1-l2+lk+ll /= 0, but less general that MatrixElement_rad_cmplx since it only computes
   !>      with dy,py=0, and therefore matrix elements are all real.
   !>    - Odd numbers of derivatives in y return purely imaginary, other combinations are purely real.
   !> **
   !------------------------------------------------------------------------------
   function MatrixElement_radx(ni,li,nj,lj,nk,lk,nl,ll,mu,dx,px) result(Vr)
     implicit none
     integer(ipr),intent(in) :: ni !<first  radial principal quantum number
     integer(ipr),intent(in) :: nj !<second radial principal quantum number
     integer(ipr),intent(in) :: nk !<third  radial principal quantum number
     integer(ipr),intent(in) :: nl !<fourth radial principal quantum number
     integer(ipr),intent(in) :: li !<first  radial orbital quantum number
     integer(ipr),intent(in) :: lj !<second radial orbital quantum number
     integer(ipr),intent(in) :: lk !<third  radial orbital quantum number
     integer(ipr),intent(in) :: ll !<fourth radial orbital quantum number
     integer(ipr),intent(in) :: dx,px !< power of extra factor x (0,1,2)
     integer(ipr),intent(in) :: mu !< Gaussian width
     real(pr) :: Vr !< Real matrix element (V_{ijkl})_{r,phi}

     real(pr) :: Ci,Cj,Ck,Cl, pf
     integer(ipr) :: nyi,nyj,nyk,nyl,nxi,nxj,nxk,nxl, ip

     real(pr), dimension(:,:,:,:,:), pointer :: ME1Dx

     nullify(ME1Dx)
     if (dx==0.and.px==0) then
       ME1Dx => ME1Dr
     else if (dx==1.and.px==0) then
       ME1Dx => ME1Dr_10
     else if (dx==2.and.px==0) then
       ME1Dx => ME1Dr_20
     else if (dx==1.and.px==1) then
       ME1Dx => ME1Dr_11
     else if (dx==1.and.px==2) then
       ME1Dx => ME1Dr_12
     else if (dx==0.and.px==1) then
       ME1Dx => ME1Dr_01
     else if (dx==0.and.px==2) then
       ME1Dx => ME1Dr_02
     else
        stop "Invalid input"
     end if

     ! Allow computing for -li-lj+lk+ll /=0 because of x-derivatives.
     ! This selection rule is elevated to ddGr.
     Vr=0
     do nyi = 0,2*ni+abs(li)
        nxi = 2*ni+abs(li) - nyi
        Ci = Cp2c(ni,li,nyi)
        do nyj = 0,2*nj+abs(lj)
           nxj = 2*nj+abs(lj)-nyj
           Cj = Cp2c(nj,lj,nyj)
           do nyk = 0,2*nk+abs(lk)
              nxk = 2*nk+abs(lk) - nyk
              Ck = Cp2c(nk,lk,nyk)
              ! Skip values that lead to imaginary prefactor
              do nyl = mod(nyi+nyj+nyk,2),2*nl+abs(ll),2
                 nxl = 2*nl+abs(ll)-nyl
                 Cl = Cp2c(nl,ll,nyl)
                 Vr = Vr + &
                      (-1)**(nyi+nyj+(nyi+nyj+nyk+nyl)/2)*Ci*Cj*Ck*Cl &
                      *ME1Dx(mu,nxi,nxj,nxk,nxl)&
                      *ME1Dr(mu,nyi,nyj,nyk,nyl)
              enddo
           enddo
        enddo
     enddo

   end function MatrixElement_radx

   !------------------------------------------------------------------------------
   !> Calculates radial(r,phi)-component of two-body matrix elements of the interaction:
   !> \f[ V = \partial_i \partial_j e^{-(\vec{r}_1 - \vec{r}_2)^2]/\mu^2} \f]
   !> in the axial HO basis. The derivatives d_i are components of the gradient in the
   !> spherical basis, acting on the Gaussian. Optionally, d_i can be switched to 2iP_i,
   !> where P_i=-i/2(grad_right_i - grad_left_i) acts on wavefunctions.
   !>
   !> The components of the gradient in the spherical basis are:
   !> \f[ \partial_{+1} = -1/\sqrt{2} (\partial_x + i \partial_y) \f]
   !> \f[ \partial_{0}  = \partial_z \f]
   !> \f[ \partial_{+1} = +1/\sqrt{2} (\partial_x - i \partial_y) \f]
   !>
   !> This function uses the full complex radial matrix element routine MatrixElement_rad_cmplx
   !------------------------------------------------------------------------------
   function MatrixElement_ddGr_cmplx(ni,li,nj,lj,nk,lk,nl,ll,mu,di,dj,dip) result(Vrc)
     implicit none
     integer(ipr),intent(in) :: ni !<first  radial principal quantum number
     integer(ipr),intent(in) :: nj !<second radial principal quantum number
     integer(ipr),intent(in) :: nk !<third  radial principal quantum number
     integer(ipr),intent(in) :: nl !<fourth radial principal quantum number
     integer(ipr),intent(in) :: li !<first  radial orbital quantum number
     integer(ipr),intent(in) :: lj !<second radial orbital quantum number
     integer(ipr),intent(in) :: lk !<third  radial orbital quantum number
     integer(ipr),intent(in) :: ll !<fourth radial orbital quantum number
     integer(ipr),intent(in) :: mu    !< Gaussian width
     integer(ipr),intent(in) :: di,dj !> Component of gradient acting on V
     integer(ipr),intent(in) :: dip   !> 0 = ignore, 1,2=change di to 2i pi(x_dip)
     complex(2*pr) :: Vrc, pxdx, pydy, pydx, pxdy, px, py
     integer(ipr) :: p

     Vrc = 0

     ! If di --> pi, reduce number d by 1, use dip=(0,1,2) for p part
     if (dip/=0) then
        p = 1
     else
        p = 0
     end if

     ! p0d0 = [1]_r * [pz dz]_z
     if (di==0 .and. dj==0) then
        Vrc = MatrixElement_rad_cmplx(ni,li,nj,lj,nk,lk,nl,ll,mu,0,0, 0,0)

     ! p+d+ = [1/2(pxdx-pydy+ipxdy+ipydx)]_r * [1]_z
     else if (di==+1 .and. dj==+1) then
        pxdx = MatrixElement_rad_cmplx(ni,li,nj,lj,nk,lk,nl,ll,mu,2-p,0, dip,0)
        pydy = MatrixElement_rad_cmplx(ni,li,nj,lj,nk,lk,nl,ll,mu,0,2-p, 0,dip)
        pxdy = MatrixElement_rad_cmplx(ni,li,nj,lj,nk,lk,nl,ll,mu,1-p,1, dip,0)
        pydx = pxdy
        if (dip/=0) pydx = MatrixElement_rad_cmplx(ni,li,nj,lj,nk,lk,nl,ll,mu,1,1-p, 0,dip)

        Vrc = pxdx - pydy + cmplx(0.,1.,2*pr)*pxdy + cmplx(0.,1.,2*pr)*pydx
        Vrc = 0.5_pr*Vrc

     ! p-d- = [1/2(pxdx-pydy-ipxdy-ipydx)]_r * [1]_z
     else if (di==-1 .and. dj==-1) then
        pxdx = MatrixElement_rad_cmplx(ni,li,nj,lj,nk,lk,nl,ll,mu,2-p,0, dip,0)
        pydy = MatrixElement_rad_cmplx(ni,li,nj,lj,nk,lk,nl,ll,mu,0,2-p, 0,dip)
        pxdy = MatrixElement_rad_cmplx(ni,li,nj,lj,nk,lk,nl,ll,mu,1-p,1, dip,0)
        pydx = pxdy
        if (dip/=0) pydx = MatrixElement_rad_cmplx(ni,li,nj,lj,nk,lk,nl,ll,mu,1,1-p, 0,dip)

        Vrc = pxdx - pydy - cmplx(0.,1.,2*pr)*pxdy - cmplx(0.,1.,2*pr)*pydx
        Vrc = 0.5_pr*Vrc

     ! p+d0 = [-1/sqrt(2)(px+ipy)]_r * [dz]_z
     ! d+p0 = [-1/sqrt(2)(dx+idy)]_r * [pz]_z
     else if ((di==+1 .and. dj==0) .or. (di==0 .and. dj==+1)) then
        if (di==0 .and. dj==+1) then
           px = MatrixElement_rad_cmplx(ni,li,nj,lj,nk,lk,nl,ll,mu,1,0, 0,0)
           py = MatrixElement_rad_cmplx(ni,li,nj,lj,nk,lk,nl,ll,mu,0,1, 0,0)
        else
           px = MatrixElement_rad_cmplx(ni,li,nj,lj,nk,lk,nl,ll,mu,1-p,0, dip,0)
           py = MatrixElement_rad_cmplx(ni,li,nj,lj,nk,lk,nl,ll,mu,0,1-p, 0,dip)
        end if

        Vrc = px + cmplx(0.,1.,2*pr)*py
        Vrc = -Vrc/cr2

     ! p-d0 = [+1/sqrt(2)(px-ipy)]_r * [dz]_z
     ! d-p0 = [+1/sqrt(2)(dx-idy)]_r * [pz]_z
     else if ((di==-1 .and. dj==0) .or. (di==0 .and. dj==-1)) then
        if (di==0 .and. dj==-1) then
           px = MatrixElement_rad_cmplx(ni,li,nj,lj,nk,lk,nl,ll,mu,1,0, 0,0)
           py = MatrixElement_rad_cmplx(ni,li,nj,lj,nk,lk,nl,ll,mu,0,1, 0,0)
        else
           px = MatrixElement_rad_cmplx(ni,li,nj,lj,nk,lk,nl,ll,mu,1-p,0, dip,0)
           py = MatrixElement_rad_cmplx(ni,li,nj,lj,nk,lk,nl,ll,mu,0,1-p, 0,dip)
        end if

        Vrc = px - cmplx(0.,1.,2*pr)*py
        Vrc = +Vrc/cr2

     ! p+d- = [-1/2(pxdx+pydy+ipydx-ipxdy)]_r * [1]_z
     ! d+p- = [-1/2(pxdx+pydy-ipydx+ipxdy)]_r * [1]_z
     else if ((di==+1 .and. dj==-1) .or. (di==-1 .and. dj==+1)) then
        pxdx = MatrixElement_rad_cmplx(ni,li,nj,lj,nk,lk,nl,ll,mu,2-p,0, dip,0)
        pydy = MatrixElement_rad_cmplx(ni,li,nj,lj,nk,lk,nl,ll,mu,0,2-p, 0,dip)
        pydx = 0
        pxdy = 0
        if (dip/=0) then
            pydx = MatrixElement_rad_cmplx(ni,li,nj,lj,nk,lk,nl,ll,mu,1,1-p, 0,dip)
            pxdy = MatrixElement_rad_cmplx(ni,li,nj,lj,nk,lk,nl,ll,mu,1-p,1, dip,0)
            if (di==-1 .and. dj==+1) then
               pydx = - pydx
               pxdy = - pxdy
            end if
        end if

        Vrc = pxdx + pydy + cmplx(0.,1.,2*pr)*pydx - cmplx(0.,1.,2*pr)*pxdy
        Vrc = -0.5_pr*Vrc

     else
        write(*,'(A)') 'Invalid components of derivatives supplied.'
        stop

     end if

   end function MatrixElement_ddGr_cmplx

   !------------------------------------------------------------------------------
   !> Calculates radial(r,phi)-component of two-body matrix elements of Gaussian interaction:
   !> \f[ V = D(dx,dy,px,py) e^{-(\vec{\rho_1} - \vec{\rho_2})^2]/\mu^2} \f]
   !> in the axial HO basis. D stands for an operator that involves 0, 1 or 2
   !> CARTESIAN DERIVATIVES (IN X AND/OR Y) acting either on the Gaussian or on wavefunctions.
   !>
   !> dx,dy = (0,1,2) is the number of cartesian derivatives acting on V (dx=d/dx, dy=d/dy)
   !> px,py = (0,1,2) means 0)no p_xi present, 1)2i p_x1, 2)2i p_x2
   !>         where p_xi=-i/2(Dxi_right - Dxi_left) and Dxi=d/dxi acts on wavefunctions.
   !>
   !> Elements are obtained by transforming the HO wavefunctions
   !> from radial into cartesian coordinates and separating the
   !> integration into a product of two one-dimensional two-body matrix
   !> elements. The transformation is given by
   !> \f[ \langle n_{r_i} \Lambda_i,n_{r_j}\Lambda_j|\hat{V}_p|
   !>     n_{r_k} \Lambda_k,n_{r_l}\Lambda_l \rangle =
   !>     \sum_{n_{y_i}=0}^{2n_{r_i}+|\Lambda_i|}
   !>     \sum_{n_{y_j}=0}^{2n_{r_j}+|\Lambda_j|}
   !>     \sum_{n_{y_k}=0}^{2n_{r_k}+|\Lambda_k|}
   !>     \sum_{n_{y_l}=0}^{2n_{r_l}+|\Lambda_l|}
   !>     C_{n_{x_i} n_{y_i}}^{n_{r_i}\Lambda_i *}
   !>     C_{n_{x_j} n_{y_j}}^{n_{r_j}\Lambda_j *}
   !>     C_{n_{x_k}n_{y_k}}^{n_{r_k}\Lambda_k}
   !>     C_{n_{x_l}n_{y_l}}^{n_{r_l}\Lambda_l}
   !>     \langle n_{x_i}n_{x_j}|\hat{V}_{\rm 1D}|n_{x_k}n_{x_l}\rangle
   !>     \langle n_{y_i}n_{y_j}|\hat{V}_{\rm 1D}|n_{y_k}n_{y_l}\rangle
   !> \f]
   !> But note that we store the real numbers \f[ C_{n_{x_i} n_{y_i}}^{n_{r_i}\Lambda_i } i^{n_{x_i}} \f]
   !>
   !> ** - This function is an on-the-fly implementation of "HFBTHO_GOGNY.radial_matrix_elements"
   !>      with the additional ability to include derivatives in x or y.
   !>    - Here gogny_hyper==1 is required to treat the cartesian derivatives.
   !>    - This routine is more general than the one in HFBTHO_GOGNY since it calculates matrix elements
   !>      when -l1-l2+lk+ll /= 0, and it is not restricted to real results, but it is not as efficient.
   !>    - Odd numbers of derivatives in y return purely imaginary, other combinations are purely real.
   !> **
   !------------------------------------------------------------------------------
   function MatrixElement_rad_cmplx(ni,li,nj,lj,nk,lk,nl,ll,mu,dx,dy,px,py) result(Vr)
     implicit none
     integer(ipr),intent(in) :: ni !<first  radial principal quantum number
     integer(ipr),intent(in) :: nj !<second radial principal quantum number
     integer(ipr),intent(in) :: nk !<third  radial principal quantum number
     integer(ipr),intent(in) :: nl !<fourth radial principal quantum number
     integer(ipr),intent(in) :: li !<first  radial orbital quantum number
     integer(ipr),intent(in) :: lj !<second radial orbital quantum number
     integer(ipr),intent(in) :: lk !<third  radial orbital quantum number
     integer(ipr),intent(in) :: ll !<fourth radial orbital quantum number
     integer(ipr),intent(in) :: dx,px !< power of extra factor x (0,1,2)
     integer(ipr),intent(in) :: dy,py !< power of extra factor y (0,1,2)
     integer(ipr),intent(in) :: mu !< Gaussian width
     complex(2*pr) :: Vr !< Complex matrix element (quad precision) (V_{ijkl})_{r,phi}

     real(pr) :: Ci,Cj,Ck,Cl
     integer(ipr) :: nyi,nyj,nyk,nyl,nxi,nxj,nxk,nxl

     real(pr), dimension(:,:,:,:,:), pointer :: ME1Dx, ME1Dy

     nullify(ME1Dx)
     if (dx==0.and.px==0) then
       ME1Dx => ME1Dr
     else if (dx==1.and.px==0) then
       ME1Dx => ME1Dr_10
     else if (dx==2.and.px==0) then
       ME1Dx => ME1Dr_20
     else if (dx==1.and.px==1) then
       ME1Dx => ME1Dr_11
     else if (dx==1.and.px==2) then
       ME1Dx => ME1Dr_12
     else if (dx==0.and.px==1) then
       ME1Dx => ME1Dr_01
     else if (dx==0.and.px==2) then
       ME1Dx => ME1Dr_02
     else
        stop "Invalid input"
     end if
     nullify(ME1Dy)
     if (dy==0.and.py==0) then
       ME1Dy => ME1Dr
     else if (dy==1.and.py==0) then
       ME1Dy => ME1Dr_10
     else if (dy==2.and.py==0) then
       ME1Dy => ME1Dr_20
     else if (dy==1.and.py==1) then
       ME1Dy => ME1Dr_11
     else if (dy==1.and.py==2) then
       ME1Dy => ME1Dr_12
     else if (dy==0.and.py==1) then
       ME1Dy => ME1Dr_01
     else if (dy==0.and.py==2) then
       ME1Dy => ME1Dr_02
     else
        stop "Invalid input"
     end if

     Vr=0
     do nyi = 0,2*ni+abs(li)
        nxi = 2*ni+abs(li) - nyi
        Ci = Cp2c(ni,li,nyi)
        do nyj = 0,2*nj+abs(lj)
           nxj = 2*nj+abs(lj)-nyj
           Cj = Cp2c(nj,lj,nyj)
           do nyk = 0,2*nk+abs(lk)
              nxk = 2*nk+abs(lk) - nyk
              Ck = Cp2c(nk,lk,nyk)
              do nyl = 0,2*nl+abs(ll)
                 nxl = 2*nl+abs(ll)-nyl
                 Cl = Cp2c(nl,ll,nyl)
                 Vr = Vr + &
                      cmplx(0,1,2*pr)**(nyi+nyj-nyk-nyl)*Ci*Cj*Ck*Cl &
                      *ME1Dx(mu,nxi,nxj,nxk,nxl) &
                      *ME1Dy(mu,nyi,nyj,nyk,nyl)
              enddo
           enddo
        enddo
     enddo

   end function MatrixElement_rad_cmplx

   !------------------------------------------------------------------------------
   !> Calculates 1D-cartesian component of two-body matrix elements of Gaussian interaction:
   !> \f[ V = D(d,p) e^{-(x_1-x_2)^2/\mu^2} \f]
   !> in the axial HO basis. D stands for an operator that involves 0, 1 or 2 derivatives
   !> acting either on the Gaussian or on wavefunctions.
   !>
   !> d = (0,1,2) is the number of derivatives acting on V
   !> p = (0,1,2) means 0)no p_xi present, 1)2i p_x1, 2)2i p_x2
   !>     where p_xi=-i/2(Dxi_right - Dxi_left) and Dxi=d/dxi acts on wavefunctions
   !>
   !> THIS ROUTINE CAN PROBABLY BE MADE MORE ROBUST, WHERE D=IBP, 2iP=R-L.
   !> I JUST WORKED OUT ALL THE COMBINATIONS BY HAND...
   !------------------------------------------------------------------------------
   function MatrixElement_1d(ni,nj,nk,nl,mu,z,d,p) result(V1d)
     implicit none
     integer(ipr), intent(in) :: ni !< first  z-component quantum number
     integer(ipr), intent(in) :: nj !< second z-component quantum number
     integer(ipr), intent(in) :: nk !< third  z-component quantum number
     integer(ipr), intent(in) :: nl !< fourth z-component quantum number
     integer(ipr), intent(in) :: d  !< number of derivatives on V
     integer(ipr), intent(in) :: p  !< momentum operator present (0=no, 1=px1, 2=px2)
     integer(ipr), intent(in) :: mu !< Gaussian width
     logical, intent(in) :: z  !< true = z, false = r
     real(pr) :: V1d !< Real matrix element V_{ijkl}

     real(pr), dimension(:,:,:,:,:), pointer :: ME1D

     nullify(ME1D)
     if (z) then
       ME1D => ME1Dz
     else
       ME1D => ME1Dr
     end if

     V1d = 0

     ! Z-Symmetries (note exchange a,c or b,d is presevered, but changes sign for P terms)
     if (p==0) then
        if(mod(ni+nj+nk+nl+d,2).ne.0) return
     else
        if(mod(ni+nj+nk+nl+d+1,2).ne.0) return
        if(p==1) then
           if (ni==nk) return
        else if(p==2) then
           if (nj==nl) return
        end if
     end if

     ! G=Gaussian
     if (d==0 .and. p==0) then
        V1d = ME1D(mu, ni,nj,nk,nl)

     ! dG = 1 integration by parts on x1
     ! dG = -(i'jkl) - (ijk'l)
     else if (d==1 .and. p==0) then
        V1d =       MatrixElement_1d_dho(ni,nj,nk,nl, 1,0,0,0, mu,z)*(-1.0_pr)
        V1d = V1d + MatrixElement_1d_dho(ni,nj,nk,nl, 0,0,1,0, mu,z)*(-1.0_pr)

     ! dG = 2 integration by parts on x1
     ! d2G = +(i''jkl) + (ijk''l) + 2(i'jk'l)
     else if (d==2 .and. p==0) then
        V1d =       MatrixElement_1d_dho(ni,nj,nk,nl, 2,0,0,0, mu,z)
        V1d = V1d + MatrixElement_1d_dho(ni,nj,nk,nl, 0,0,2,0, mu,z)
        V1d = V1d + MatrixElement_1d_dho(ni,nj,nk,nl, 1,0,1,0, mu,z)*(2.0_pr)

     ! 2ipx1 dG = 1 integration by parts (on x1)  +  (R-L) on x1
     ! 2ipx2 dG = -[ +(i'jkl') + (i'j'kl) - (ijk'l') - (ij'k'l) ]
     !          = -[ +(ijk''l) - (i''jkl) ]
     else if (d==1 .and. p==1) then
        V1d =       MatrixElement_1d_dho(ni,nj,nk,nl, 0,0,2,0, mu,z)*(-1.0_pr)
        V1d = V1d + MatrixElement_1d_dho(ni,nj,nk,nl, 2,0,0,0, mu,z)

     ! 2ipx2 dG = 1 integration by parts (on x1)  +  (R-L) on x2
     ! 2ipx2 dG = -[ +(i'jkl') - (i'j'kl) + (ijk'l') - (ij'k'l) ]
     !          = +[ +(ijkl'') - (ij''kl) ]
     else if (d==1 .and. p==2) then
        V1d =       MatrixElement_1d_dho(ni,nj,nk,nl, 0,0,0,2, mu,z)
        V1d = V1d + MatrixElement_1d_dho(ni,nj,nk,nl, 0,2,0,0, mu,z)*(-1.0_pr)

     ! 2ipx1 G = (R-L) on x1
     ! 2ipx1 G = +(ijk'l) - (i'jkl)
     else if (d==0 .and. p==1) then
        V1d =       MatrixElement_1d_dho(ni,nj,nk,nl, 0,0,1,0, mu,z)
        V1d = V1d + MatrixElement_1d_dho(ni,nj,nk,nl, 1,0,0,0, mu,z)*(-1.0_pr)

     ! 2ipx2 G = (R-L) on x2
     ! 2ipx2 G = +(ijkl') - (ij'kl)
     else if (d==0 .and. p==2) then
        V1d =       MatrixElement_1d_dho(ni,nj,nk,nl, 0,0,0,1, mu,z)
        V1d = V1d + MatrixElement_1d_dho(ni,nj,nk,nl, 0,1,0,0, mu,z)*(-1.0_pr)

     else
        write(*,'(A)') 'Supplied combination of derivatives d,p not implemented.'
        stop

     end if

   end function MatrixElement_1d

   !------------------------------------------------------------------------------
   !> Compute 1D matrix elements with up to 2 derivatives on any HO WF.
   !>
   !> We use negative quantum numbers to indicate combinations to skip. This works
   !> nicely since the new quantum numbers might be n-1 or n-2, which should also be
   !> neglected if that value is negative.
   !>
   !> ** Wraps HFBTHO_GOGNY.MatrixElement_z **
   !------------------------------------------------------------------------------
   function MatrixElement_1d_dho(ni,nj,nk,nl, di,dj,dk,dl, mu,z) result(V1d)
     implicit none
     integer(ipr), intent(in) :: ni,di !< 1st quantum number, n-derivs on WF
     integer(ipr), intent(in) :: nj,dj !< 2nd quantum number, n-derivs on WF
     integer(ipr), intent(in) :: nk,dk !< 3rd quantum number, n-derivs on WF
     integer(ipr), intent(in) :: nl,dl !< 4th quantum number, n-derivs on WF
     integer(ipr), intent(in) :: mu !< Gaussian width
     logical, intent(in) :: z  !< true = z, false = r

     integer(ipr) :: nni
     integer(ipr) :: nnj
     integer(ipr) :: nnk
     integer(ipr) :: nnl
     real(pr) :: V1d, V1d_i !< Real matrix element V_{ijkl}

     integer(ipr), dimension(3) :: nit, njt, nkt, nlt !< temporary array
     real(pr),    dimension(3) :: cit, cjt, ckt, clt !< temporary array
     integer(ipr) :: i, j, k, l, inz, jnz, knz, lnz
     real(pr) :: cci,ccj,cck,ccl

     real(pr), dimension(:,:,:,:,:), pointer :: ME1D


     nullify(ME1D)
     if (z) then
       ME1D => ME1Dz
     else
       ME1D => ME1Dr
     end if

     call ho_derivative(ni,di,z,nit,cit,inz)
     call ho_derivative(nj,dj,z,njt,cjt,jnz)
     call ho_derivative(nk,dk,z,nkt,ckt,knz)
     call ho_derivative(nl,dl,z,nlt,clt,lnz)
     V1d = 0.0_pr
     do l=1,lnz
        nnl = nlt(l)
        ccl = clt(l)
        do k=1,knz
           nnk = nkt(k)
           cck = ckt(k)
           do j=1,jnz
              nnj = njt(j)
              ccj = cjt(j)
              do i=1,inz
                 nni = nit(i)
                 cci = cit(i)
                 ! By reflection symmetry these will be zero
                 if(mod(nni+nnj+nnk+nnl,2).ne.0) cycle
                 V1d_i = ME1D(mu, nni, nnj, nnk, nnl)
                 V1d_i = V1d_i*cci*ccj*cck*ccl
                 V1d = V1d + V1d_i
              end do
           end do
        end do
     end do

   end function MatrixElement_1d_dho

   !------------------------------------------------------------------------------
   !> Derivatives of HO WFs can be expressed as a sum of of HO WFs with different
   !> quantum numbers (and prefactors). This routine computes arrays containing these
   !> quantum numbers (narr) and prefactors (carr) for 0,1,2 derivatives. The results have
   !> (0 derivs = 1 term), (1 derivs = 2 terms), (2 derivs = 3 terms), so all outputs are
   !> size=3. Unused elements of the output arrays are defaulted to carr(i)=0 and
   !> narr(i)=-1. If any new quantum number is narr(i) < 0, carr(i) is forced to 0.
   !------------------------------------------------------------------------------
   subroutine ho_derivative(n,d,z, narr, carr, nnz)

     implicit none
     integer(ipr), intent(in) :: n !< Principle HO quantum number
     integer(ipr), intent(in) :: d !< number of derivatives on HO WF (0,1,2)
     logical, intent(in) :: z !< z or radial part of HO WFs?

     integer(ipr), dimension(3), intent(out) :: narr !< new HO quantum numbers
     real(pr),    dimension(3), intent(out) :: carr !< HO WF coeffs
     integer(ipr), intent(out) :: nnz

     real(pr), pointer :: prefd1, prefd2
     integer(ipr) :: i

     nullify(prefd1,prefd2)
     if (z) then
        prefd1 => prefd1z
        prefd2 => prefd2z
     else
        prefd1 => prefd1r
        prefd2 => prefd2r
     end if

     ! Use negative n to indicate elements to skip (since they're zero)
     carr = 0.0_pr
     narr = -1
     nnz  = 1

     ! Override values based on requested derivative
     select case (d)

        case default
           write(*,'(A)') 'Supplied number of derivatives not implemented.'
           stop

        case (0)
           carr(1) = 1.0_pr
           narr(1) = n

        case (1)
           carr(1) = -sqrt(real(n+1,pr))*prefd1
           narr(1) = n+1
           if (n-1>=0) then
              carr(2) = sqrt(real(n,pr))*prefd1
              narr(2) = n-1
              nnz = 2
           end if

        case (2)
           nnz = 2
           carr(1) = sqrt(real((n+1)*(n+2),pr))*prefd2
           carr(2) = -(2.0_pr*n + 1.0_pr)*prefd2
           narr(1) = n+2
           narr(2) = n
           if (n-2>=0) then
              carr(3) = sqrt(real(n*(n-1),pr))*prefd2
              narr(3) = n-2
              nnz = 3
           end if

     end select

   end subroutine ho_derivative

  !!======================================================================
  !!>  Given four radial and four angular quantum numbers, returns the
  !!>  corresponding index in the one-dimensional array that contains all
  !!>  non-zero and different radial components of two body potential
  !!>  potential matrix elements for a gaussian potential
  !!>
  !!>  If the quantum numbers correspond to a matrix element such that
  !!>  \f$ V^r_{ijkl} = 0\f$, the function returns 0 as the index
  !!=====================================================================
  !function rindex(nri,nrj,nrk,nrl,li,lj,lk,ll,n) result(ir)
  !  implicit none
  !  integer(ipr) :: nri !< First  radial quantum number
  !  integer(ipr) :: nrj !< Second radial quantum number
  !  integer(ipr) :: nrk !< Third  radial quantum number
  !  integer(ipr) :: nrl !< Fourth radial quantum number
  !  integer(ipr) :: li !< First  angular quantum number
  !  integer(ipr) :: lj !< Second angular quantum number
  !  integer(ipr) :: lk !< Third  angular quantum number
  !  integer(ipr) :: ll !< Fourth angular quantum number
  !  integer(ipr) :: n !< Total number of shells
  !  integer(ipr) :: ir
  !  
  !  !if(li+lj.eq.ll+lk) then
  !  !   if(lk.le.-li.and.ll.le.-lj) then
  !  !     ir=nblock(nri,nrj,nrk,nrl,li,lj)+&
  !  !          l_block(nrk,nrl,li,lj,lk,n)
  !  !   elseif(ll.le.-lj) then
  !  !     ir=nblock(nrk,nrj,nri,nrl,-lk,lj)+&
  !  !          l_block(nri,nrl,-lk,lj,-li,n)
  !  !   elseif(lk.le.-li) then
  !  !     ir=nblock(nri,nrl,nrk,nrj,li,-ll)+&
  !  !          l_block(nrk,nrj,li,-ll,lk,n)
  !  !   else
  !  !     ir=nblock(nrk,nrl,nri,nrj,-lk,-ll)+&
  !  !          l_block(nri,nrj,-lk,-ll,-li,n)
  !  !   endif
  !  !else
  !  !   ir = 0
  !  !endif
  !end function rindex

  !!======================================================================
  !!>  Given four axial quantum numbers, returns the corresponding index
  !!>  in the one-dimensional array that contains all non-zero and
  !!>  different axial components.
  !!>
  !!>  If the sum of the four quantum numbers is not an even number, in
  !!>  which case \f$ V^z_{ijkl} = 0\f$, the function returns 0 as the
  !!>  index
  !!======================================================================
  !function zindex(nzi,nzj,nzk,nzl) result(iz)
  !  implicit none
  !  integer(ipr) :: nzi !< First  axial quantum number
  !  integer(ipr) :: nzj !< Second axial quantum number
  !  integer(ipr) :: nzk !< Third  axial quantum number
  !  integer(ipr) :: nzl !< Fourth axial quantum number
  !  integer(ipr) :: iz

  !  ! For K=0, exchange of a,c or b,d changes sign on p terms
  !  ! For K/=0 and non-p-terms are same on exchange of indices
  !  iz = ZBlock(nzi,nzj,nzk,nzl)
  !  ! iz = ZBlock(min(nzi,nzk),min(nzj,nzl),max(nzi,nzk),max(nzj,nzl))
  !end function zindex

  !!======================================================================
  !!>  Calculates and stores in memory what is called the ZBlock. ZBlock
  !!>  is a (N x N x N x N) array that is used to obtain the
  !!>  transformation from the axial component quantum numbers
  !!>  (\f$n_z\f$'s) into the index of an array that contains only the
  !!>  the non-zero axial two body potential matrix elements.
  !!>
  !!>  In particular the ZBlock is
  !!>  \f{eqnarray*}{
  !!>   NB &=&\sum_{i=0}^{n_{z_i}-1} \sum_{j=0}^{N} \sum_{k=i}^{N}
  !!>      \sum_{\substack{l=j\\ i+j+k+l\ {\rm is\ even}}}^{N} 1
  !!>     +\sum_{j=0}^{n_{z_j}-1} \sum_{k=n_{z_i}}^{N}
  !!>      \sum_{\substack{l=j\\ n_{z_i}+j+k+l\ {\rm is\ even}}}^{N} 1 \\
  !!>      & + & \sum_{k=n_{z_i}}^{n_{z_k}-1}
  !!>            \sum_{\substack{l=n_{z_j}\\ n_{z_i}+n_{z_j}+k+
  !!>                  l\ {\rm is\ even}}}^{N} 1
  !!>           +\sum_{\substack{l=n_{z_j}\\ n_{z_i}+n_{z_j}+n_{z_k}
  !!>                  +l\ {\rm is\ even}}}^{n_{z_l}-1} 1  + 1
  !!>  \f}
  !!>  where \f$N\f$ is the number of shells
  !!======================================================================
  !subroutine calculate_Zblock()
  !  implicit none
  !  integer(ipr) :: sni,snj,snk,snl
  !  integer(ipr) :: xni,xnj,xnk,xnl
  !  integer(ipr) :: ii,jj,kk,ll,ix,n
  !  n = nzx
  !  if(allocated(ZBlock)) deallocate(ZBlock)
  !  allocate(ZBlock(0:n,0:n,0:n,0:n))
  !  Zblock = 0
  !  xni = 0
  !  do ii = 0,n
  !     sni = xni
  !     xnj = 0
  !     do jj = 0,n
  !        snj = xnj
  !        xnk = 0
  !        do kk = 0,n!ii,n
  !           snk = xnk
  !           xnl = 0
  !           do ll = 0,n!jj+mod(ii+kk,2),n,2
  !              snl = xnl
  !              ZBlock(ii,jj,kk,ll)=sni+snj+snk+snl+1
  !              xni = xni + 1
  !              xnj = xnj + 1
  !              xnk = xnk + 1
  !              xnl = xnl + 1
  !           enddo
  !        enddo
  !     enddo
  !  enddo
  !end subroutine calculate_Zblock

  !!======================================================================
  !!>  Calculates and stores in memory what is called the NBlock. NBlock
  !!>  is an (N/2 x N/2 x N/2 x N/2 x 2N x 2N) array that is used to
  !!>  obtain the transformation from the radial component quantum numbers
  !!>  (\f$n\f$'s and \f$\Lambda\f$'s) into the index of an array that
  !!>  contains only the non-zero two-body potential matrix elements.
  !!>
  !!>  The Nblock is the part of that index that only depends on the
  !!>  \f$n\f$ quantum numbers and the first two \f$\Lambda\f$ quantum
  !!>  number. This part takes the most time to calculate (and
  !!>  is therefore stored in memory so that it's only calculated once).
  !!>  The other part of that index is calculated with L_block().
  !!>
  !!>  In particular the NBlock is
  !!>  \f{eqnarray*}{
  !!>   NB &=& \sum_{i=0}^{n_{ri} -1} \sum_{j=0}^{N/2} \sum_{k=0}^{N/2}
  !!>          \sum_{l=0}^{N/2} \sum_{p=-N+2i}^{N-2i}\;
  !!>          \sum_{q=-N+2j}^{N-2j}\;
  !!>          \sum_{r=\max(-N+2k,-\min(N-l,-q)+p+q)}
  !!>              ^{\min(\min(N-2k,-p),N-2l+p+q)} 1 \\
  !!>      &+& \sum_{j=0}^{n_{rj}-1} \sum_{k=0}^{N/2}
  !!>          \sum_{l=0}^{N/2} \sum_{p=-N+2n_{ri}}^{N-2n_{ri}}\;
  !!>          \sum_{q=-N+2j}^{N-2j}\;
  !!>          \sum_{r=\max(-N+2k,-\min(N-l,-q)+p+q)}
  !!>              ^{\min(\min(N-2k,-p),N-2l+p+q)} 1 \\
  !!>      &+& \sum_{k=0}^{n_{rk}-1} \sum_{l=0}^{N/2}
  !!>          \sum_{p=-N+2n_{ri}}^{N-2n_{ri}}\;
  !!>          \sum_{q=-N+2n_{rj}}^{N-2n_{rj}}\;
  !!>          \sum_{r=\max(-N+2k,-\min(N-l,-q)+p+q)}
  !!>              ^{\min(\min(N-2k,-p),N-2l+p+q)} 1 \\
  !!>      &+& \sum_{l=0}^{n_{rl}-1} \sum_{p=-N+2n_{ri}}^{N-2n_{ri}}\;
  !!>          \sum_{q=-N+2n_{rj}}^{N-2n_{rj}}\;
  !!>          \sum_{r=\max(-N+2n_{rk},-\min(N-l,-q)+p+q)}
  !!>              ^{\min(\min(N-2k,-p),N-2l+p+q)} 1 \\
  !!>      &+& \sum_{p=-N+2n_{ri}}^{N-2n_{ri}}\;
  !!>          \sum_{q=-N+2n_{rj}}^{N-2n_{rj}}\;
  !!>          \sum_{r=\max(-N+2n_{rk},-\min(N-n_{rl},-q)+p+q)}
  !!>              ^{\min(\min(N-2k,-p),N-2l+p+q)} 1 \\
  !!>      &+& \sum_{q=-N+2n_{rj}}^{N-2n_{rj}}\;
  !!>          \sum_{r=\max(-N+2n_{rk},-\min(N-n_{rl},-q)+\Lambda_i+q)}
  !!>              ^{\min(\min(N-2k,-\Lambda_i),N-2l+\Lambda_i+q)} 1
  !!>  \f}
  !!>  where \f$N\f$ is the number of shells
  !!======================================================================
  !subroutine calculate_Nblock()
  !  implicit none
  !  integer(ipr) :: sni,snj,snk,snl,sli,slj
  !  integer(ipr) :: xni,xnj,xnk,xnl,xli,xlj
  !  integer(ipr) :: ii,jj,kk,ll,i,j,k,ix,n2,n
  !  n = max(2*nrx,nlx)
  !  n2 = n/2
  !  if(allocated(NBlock)) deallocate(NBlock)
  !  allocate(NBlock(0:n2,0:n2,0:n2,0:n2,-n:n,-n:n))
  !  NBlock = 0
  !  xni = 0
  !  do ii = 0,n2
  !     sni = xni
  !     xnj = 0
  !     do jj = 0,n2
  !        snj = xnj
  !        xnk = 0
  !        do kk = 0,n2
  !           snk = xnk
  !           xnl = 0
  !           do ll = 0,n2
  !              snl = xnl
  !              xli = 0
  !              do i = -N+2*ii,N-2*ii
  !                 sli = xli
  !                 xlj = 0
  !                 do j = -N+2*jj,N-2*jj
  !                    slj = xlj
  !                    Nblock(ii,jj,kk,ll,i,j)=sni+snj+snk+snl+sli+slj
  !                    do k = max(-n+2*kk,-min(n-2*ll,-j)+i+j), &
  !                         min(min(n-2*kk,-i),n-2*ll+i+j)
  !                       xni = xni + 1
  !                       xnj = xnj + 1
  !                       xnk = xnk + 1
  !                       xnl = xnl + 1
  !                       xli = xli + 1
  !                       xlj = xlj + 1
  !                    enddo
  !                 enddo
  !              enddo
  !           enddo
  !        enddo
  !     enddo
  !  enddo
  !end subroutine calculate_Nblock

  !!======================================================================
  !!>  Calculates what  is called the LBlock. LBlock
  !!>  is a number that is used to obtain the
  !!>  transformation from the radial component quantum numbers
  !!>  (\f$n\f$'s and \f$\Lambda\f$'s) into the index of an array that
  !!>  contains only the non-zero Two body potential matrix elements.
  !!>
  !!>  The Lblock is the part of that
  !!>  index that depends on five of the radial quantum  numbers
  !!>  (which would require an \f$2 N^5\f$ array to store in memory)
  !!>  but can be calculated on the fly rather quickly. The other part of
  !!>  that index is called the NBlock and is calculated and stored in
  !!>  memory by calculate_nblock().
  !!>
  !!>  @result \f$
  !!>   LB = 1 + \Lambda_k - \max(-N+2n_{rk},-\min(N-2n_{rl},-\Lambda_j)
  !!>                               +\Lambda_i+\Lambda_j) \f$,
  !!>  where \f$N\f$ is the number of shells.
  !!======================================================================
  !function L_block(kn,ln,il,jl,kl,n) result(Lb)
  !  implicit none
  !  integer(ipr), intent(in) :: kn !< third radial quantum number
  !  integer(ipr), intent(in) :: ln !< fourth radial quantum number
  !  integer(ipr), intent(in) :: il !< first orbital quantum number
  !  integer(ipr), intent(in) :: jl !< second orbital quantum number
  !  integer(ipr), intent(in) :: kl !< third orbital quantum number
  !  integer(ipr), intent(in) :: n !< total number of shells
  !  integer(ipr) :: Lb

  !  Lb = 1 + kl - max(-N+2*kn,-min(n-2*ln,-jl)+il+jl)

  !end function L_block


end module pnfam_spatial_mtxels
