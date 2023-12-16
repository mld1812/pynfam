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
!                      FINITE-RANGE GOGNY PACKAGE                      !
!                                                                      !
! ==================================================================== !

!----------------------------------------------------------------------
!>  This module provides all the routines needed to compute the matrix
!>  elements and expectation value of the finite-range Gogny force.
!>  The separability of the matrix elements in axial and radial
!>  components is widely used to signinificantly reduce computation time.
!>  The matrix elements are computed using the Gogny expansions outlined
!>  in \cite younes2009-b
!>
!>  @author
!>  Rodrigo Navarro Perez
!----------------------------------------------------------------------
!  Subroutines: - gogny_matrix_elements
!               - CalculateVzGogny
!               - CalculateVrGogny
!               - calculateTz
!               - calculateTr
!               - calculateCpolar2cartesian
!               - calculateME1D
!               - radial_matrix_elements(ni,li,nj,lj,nk,lk,nl,ll)
!               - calculate_Zblock
!               - calculate_Nblock
!               - trace_product_2(A,B,tr1,tr2)
!               - bdiag_trace(A,tr1,tr2)
!               - bdiag_print(A)
!               - test_HOWF_gauss
!               - LaguerreL(n,alpha,x,Ln,Lnp,Lnm1)
!               - GaussLaguerreWX(alfa,w,x)
!  Functions: - TrCoefficient(n1,k1,n2,k2,n,k)
!             - TrSumTerms(n1,k1,n2,k2,n,k)
!             - MatrixElement_z(ni,nj,nk,nl,mu,b)
!             - Ibarz(m,n,Gz)
!             - MatrixElement_r(ni,li,nj,lj,nk,lk,nl,ll,mu,b)
!             - Ibarr(n1,k1,n2,k2,Gp)
!             - factrl(n)
!             - HyperGeom2F1(a,b,c,x)
!             - binomialco(m,n)
!             - upperfactrl(x,i)
!             - l_block
!             - zindex(nzi,nzj,nzk,nzl)
!             - rindex(nri,nrj,nrk,nrl,li,lj,lk,ll,n)
!             - trace_product(A,B)
!             - MatrixElement_zZR(ni,nj,nk,nl)
!             - MatrixElement_rZR(ni,li,nj,lj,nk,lk,nl,ll)
!             - N_radial(nr,l)
!             - N_axial(nz)
!             - HermiteH(n,x)
!----------------------------------------------------------------------
module hfbtho_gogny
  Use HFBTHO_utilities
  Use HFBTHO
  implicit none
  ! Public variables
  integer(ipr), Public :: NumVz !< Number of matrix elements in the z-direction
  integer(ipr), Public :: NumVr !< Number of matrix elements in the perpendicular direction
  integer(ipr), Public, allocatable :: ib_zrls(:,:,:,:) !<
  integer(ipr), Public, allocatable :: i_zrls(:,:,:,:)
  real(pr), Public, allocatable :: VzGogny(:,:) !< Matrix elements of the Gogny force in the z-direction
                                                !< \f$ \langle n_z m_z | \hat{V} | n'_z m'_z \rangle \f$
  real(pr), Public, allocatable :: VrGogny(:,:) !< Matrix elements of the Gogny force in the perpendicular direction
                                                !< \f$ \langle n_r \Lambda_n m_r \Lambda_m | \hat{V} | n'_r \Lambda'_n m'_r \Lambda'_m \rangle \f$
  real(pr), Public, allocatable :: gamma_g_dir(:,:) !< Matrix of the direct part of the mean field \f$ \Gamma^{(\mathrm{dir.})}_{nm} \f$ for the Gogny force
  real(pr), Public, allocatable :: gamma_g_exc(:,:) !< Matrix of the exchange part of the mean field \f$ \Gamma^{(\mathrm{exc.})}_{nm} \f$ for the Gogny force
  real(pr), Public, allocatable :: delta_g_dir(:,:) !< Matrix of the pairing field \f$ \Delta_{nm} \f$ for the Gogny force
  real(pr), Public, allocatable :: coulf_g_dir(:,:) !< Matrix of the direct part of the mean field \f$ \Gamma^{(\mathrm{dir.})}_{nm} \f$ for the Coulomb force
  real(pr), Public, allocatable :: coulf_g_exc(:,:) !< Matrix of the exchange part of the mean field \f$ \Gamma^{(\mathrm{exc.})}_{nm} \f$ for the Coulomb force
  real(pr), Public, allocatable :: coulf_d_dir(:,:) !< Matrix of the pairing field \f$ \Delta_{nm} \f$ for the Coulomb force
#if(GOGNY_SYMMETRIES==0)
  real(pr), Public, allocatable  :: Vz_Gogny(:,:,:,:,:) !< Matrix elements of the Gogny force in the z-direction
                                                        !< \f$ \langle n_z m_z | \hat{V} | n'_z m'_z \rangle \f$
  real(pr), Public, allocatable  :: Vr_Gogny(:,:,:,:,:,:,:,:,:) !< Matrix elements of the Gogny force in the perpendicular direction
                                                                !< \f$ \langle n_r \Lambda_n m_r \Lambda_m | \hat{V} | n'_r \Lambda'_n m'_r \Lambda'_m \rangle \f$
#endif
  logical, Public :: matrix_elements_calculated
  ! Private variables
  integer(ipr), Private, allocatable :: Zblock(:,:,:,:) !< Indices corresponding to non-zero and different axial matrix elements
  integer(ipr), Private, allocatable :: NBlock(:,:,:,:,:,:) !< Indices corresponding to non-zero and different radial matrix elements
  integer(ipr), Private, allocatable :: gindx(:,:)
  real(pr), public, allocatable, dimension(:,:,:) :: T_z  !< \f$ T_{n_1 n_2}^n \f$= Tz(n1,n2,n) coefficients of the axial Gogny expansion
  real(pr), Private, allocatable, dimension(:,:,:,:,:) :: T_r  !< \f$ T_{n_1,k_1,n_2,k_2}^{n,k1+k2}\f$ = Tr(n1,k1,n2,k2,n) coefficients of the radial Gogny expansion
  real(pr), public, allocatable, dimension(:,:,:) :: Cp2c !< \f$ C^{n_r k}_{n_y}\f$= Cp2c(nr,k,ny) coefficients of the radial to cartesian
                                                           !< coordinates transformation of the Harmonic Oscillator wavefunction
  real(pr), Private, allocatable, dimension(:,:,:,:,:) :: ME1D !< Two body potential matrix elements in one dimension
                                                               !< \f$ \langle n_1 n_2|\hat{V}|n_3 n_4\rangle \f$ = ME1D(n1,n2,n3,n4)
  real(pr), allocatable, dimension(:) :: Vr_ig
!$omp threadprivate(Vr_ig)
contains
  !======================================================================
  !>  Allocates arrays related to finite-range potentials
  !======================================================================
  Subroutine allocate_fr()
    If(Allocated(gamma_g_dir)) Deallocate(gamma_g_dir,gamma_g_exc,delta_g_dir,coulf_g_dir,coulf_g_exc,coulf_d_dir)
    Allocate(gamma_g_dir(ndx**2,2*nbx),gamma_g_exc(ndx**2,2*nbx),delta_g_dir(ndx**2,2*nbx),coulf_g_dir(ndx**2,2*nbx),coulf_g_exc(ndx**2,2*nbx),coulf_d_dir(ndx**2,2*nbx))
  End Subroutine allocate_fr
  !======================================================================
  !>  Calculates the two-body matrix elements \f$ \langle n_1 n_2 |
  !>  V_i(r) | n_3 n_4 \rangle \f$ for a Gaussian potential \f$ V(r) =
  !>  e^{(\mathbf{r}_1-\mathbf{r}_2)^2/\mu_i^2}\f$. The matrix
  !>  elements are separated as the product of axial \f$
  !>  V^z_{n_1n_2n_3n_4}\f$ and radial \f$V^r_{n_1n_2n_3n_4}\f$
  !>  components.
  !>
  !>  If the preprocessor directive GOGNY_HYPER is set to 1,
  !>  \f$V^z_{n_1n_2n_3n_4} \f$ is calculated with an expansion based
  !>  on the hypergeometric function \f$_2F_1\f$. This expansion
  !>  mantains accuracy with an increasing basis size, see
  !>  matrixelement_z() for details. \f$V^r_{n_1n_2n_3n_4} \f$ is
  !>  calculated by expanding each two-dimensional harmonic oscillator
  !>  function as a sum of products of two one-dimensional harmonic
  !>  oscillator functions, see radial_matrix_elements() for more
  !>  details. This expansion allows to calculate \f$V^r\f$ with the
  !>  same subroutine used to calculate \f$V^z\f$ in order to preserve
  !>  accuracy as the basis size increases.
  !>
  !>  For any other value of GOGNY_HYPER a set of more direct
  !>  expansions is used to calculate \f$V^z\f$ and \f$V^r\f$, see
  !>  matrixelement_z() and matrixelement_r() for more details. These
  !>  expansions start to loose accuracy as the basis size increases.
  !>
  !>  The calculation of the radial components is timed and the
  !>  required wall clock time is printed as output
  !======================================================================
  subroutine gogny_matrix_elements
    implicit none
    integer(ipr) ::  t1,t2,countrate,countmax,iw
    if(matrix_elements_calculated) return
    do iw = lout,lfile
       write(iw,*)
       write(iw,*) ' Calculating the finite range matrix elements'
    enddo
    call calculateTz
#if(GOGNY_HYPER==1)
    call calculateCpolar2cartesian
    call calculateME1D
#else
    call calculateTr
#endif
    call calculateVzGogny
    call system_clock(t1,countrate,countmax)
    call calculateVrGogny
    call system_clock(t2,countrate,countmax)
    do iw = lout,lfile
       write(iw,'(a33,f15.4)') '  Matrix elements wall clock time', &
            (t2-t1)/real(countrate,kind=pr)
    enddo
    matrix_elements_calculated = .true.
  end subroutine gogny_matrix_elements

  !======================================================================
  !>  Calculates the the necessary two-body potential matrix elements
  !>  \f$\langle n_{z_i} n_{z_j}|V_z|n_{z_k} n_{z_l} \rangle \f$.
  !>
  !>  If the preprocessor variable GOGNY_SYMMETRIES is set to 1 only
  !>  states with \f$ n_{z_k} \geq n_{z_i} \f$,
  !>  \f$ n_{z_l} \geq n_{z_j} \f$ are calculated and stored. In any
  !>  case, states where \f$ n_{z_i}+n_{z_j}+n_{z_k}+n_{z_l} \f$ is not
  !>  an even number are not calculated since those matrix elements are
  !>  zero.
  !======================================================================
  subroutine CalculateVzGogny()
    implicit none
    real(pr) :: Vz1,Vz2,Vz
    integer(ipr) :: ii,N,nzi,nzj,nzk,nzl,ig
    n = nzx
#if(GOGNY_SYMMETRIES==1)
    call calculate_Zblock
    ii = Zblock(n,n,n,n); NumVz = ii
    if(allocated(VzGogny)) deallocate(VzGogny)
    allocate(VzGogny(1:n_g_all,0:ii))
    VzGogny(1:n_g_all,0) = zero
!$OMP Parallel Default(None) &
!$OMP& SHARED(n,Zblock,VzGogny,mu_g_all,bz,n_g_all) &
!$OMP& PRIVATE(nzi,nzj,nzk,nzl,ii,ig,Vz)
!$OMP DO SCHEDULE(DYNAMIC)
    do nzi = 0,n
       do nzj = 0,n
           do nzk = 0,n
             do nzl = 0,n
                if(nzk.lt.nzi) cycle
                if(nzl.lt.nzj) cycle
                if(mod(nzi+nzj+nzk+nzl,2).ne.0) cycle
                ii = Zblock(nzi,nzj,nzk,nzl)
                do ig = 1,n_g_all
                   Vz = MatrixElement_z(nzi,nzj,nzk,nzl,mu_g_all(ig),bz)
                   VzGogny(ig,ii) = Vz
                enddo
             enddo
          enddo
       enddo
    enddo
!$OMP End Do
!$OMP End Parallel
#else
    if(allocated(Vz_Gogny))  deallocate(Vz_Gogny)
    allocate(Vz_Gogny(1:n_g_all,0:n,0:n,0:n,0:n))
    Vz_gogny = zero
    do nzi = 0,n
       do nzj = 0,n
           do nzk = 0,n
             do nzl = 0,n
                if(mod(nzi+nzj+nzk+nzl,2).ne.0) cycle
                do ig = 1,n_g_all
                   Vz = MatrixElement_z(nzi,nzj,nzk,nzl,mu_g_all(ig),bz)
                   Vz_Gogny(ig,nzi,nzj,nzk,nzl) = Vz
                enddo
             enddo
          enddo
       enddo
    enddo
#endif
  end subroutine CalculateVzGogny


  !======================================================================
  !>  Calculates the necessary two-body potential
  !>  matrix elements \f$ \langle n_{r_i} \Lambda_i, n_{r_j} \Lambda_j
  !>  |V_p | n_{r_k} \Lambda_k, n_{r_l} \Lambda_l \rangle\f$.
  !>
  !>  If the preprocessor variable GOGNY_SYMMETRIES is set to 1 only
  !>  the matrix elements that will be used in the calculation of the HFB
  !>  fields are calculated. In any case, states where \f$ -\Lambda_i
  !>  - \Lambda_j + \Lambda_k + \Lambda_l \neq 0 \f$ are not calculated
  !>  since the  matrix element is zero.
  !======================================================================
  subroutine CalculateVrGogny()
    implicit none
    real(pr) :: Vr
    integer(ipr) :: ii,nri,li,nrj,lj,nrk,lk,nrl,ll,ig,im,N
    integer(ipr) :: nza,nzb,nzc,nzd,nra,nrb,nrc,nrd,oa,ob,la,lb,lc,ld
    integer(ipr) :: nla,nlb,nlc,nld,nsa,nsc,nsac,j1,j2,nsb,nsd,nsdb
    integer(ipr) :: ir_abdc,ir_abcd,ir_acbd,jr_abdc,jr_abcd,jr_acbd
    integer(ipr) :: ita,itb,itc,itd,iba,ibc,ibd,ibb,ir
    integer(ipr) , allocatable :: index_flag(:), nr_flag(:,:), &
         nl_flag(:,:),ir_flag(:)
    integer(ipr) :: ilauf
    !
    N = max(2*nrx,nlx)
    im = mod(N,2)
#if(GOGNY_SYMMETRIES==1)
    call calculate_Nblock
    ii = nblock(n/2,n/2,n/2,n/2,im,-im)+l_block(n/2,n/2,im,-im,-im,n); NumVr = ii
    if(allocated(VrGogny)) deallocate(VrGogny)
    if(.not.allocated(VrGogny)) then
       allocate(VrGogny(1:n_g_all,0:ii))
       allocate(index_flag(0:ii))
       allocate(nr_flag(1:4,0:ii))
       allocate(nl_flag(1:4,0:ii))
       allocate(ir_flag(0:ii))
    endif
    !
    index_flag = 0
    VrGogny(1:n_g_all,0) = zero
    !
    !We go over the loop that the field calculation will go through in
    !order to identify the matrix elements that will be necesary and
    !calculate only those on a subsequent loop.
    ilauf = 0
    if(.not.force_is_dme) then
    do ita = 1,ntx
       nra = nr(ita); nza = nz(ita); nla = nl(ita); nsa = ns(ita)
       if(nza.gt.1) cycle
       iba = ib_zrls(nza,nra,nla,(nsa+1)/2)
       do itc = 1, ita
          nrc = nr(itc); nzc = nz(itc); nlc = nl(itc); nsc = ns(itc)
          if(nzc.gt.1) cycle
          ibc = ib_zrls(nzc,nrc,nlc,(nsc+1)/2)
          nsac = nsa + nsc
          if(ibc.ne.iba) cycle
          do itb = 1,nttx
             nrb=nrr(itb); nlb=nll(itb); nsb=nss(itb); ibb=noo(itb)
             do itd = 1,nttx
                nrd=nrr(itd); nld=nll(itd); nsd=nss(itd); ibd=noo(itd)
                if(ibb.ne.ibd) cycle
                if(nla+nlb.ne.nlc+nld.and.nla-nlb.ne.nlc-nld) cycle
                nsdb = nsd + nsb
                if(nsac.ne.0) then
                   ir=rindex(nra,nrb,nrd,nrc,nla,nlb,nld,nlc,n)
                   if(index_flag(ir).eq.0) then
                      nr_flag(1,ilauf) = nra
                      nr_flag(2,ilauf) = nrb
                      nr_flag(3,ilauf) = nrd
                      nr_flag(4,ilauf) = nrc
                      nl_flag(1,ilauf) = nla
                      nl_flag(2,ilauf) = nlb
                      nl_flag(3,ilauf) = nld
                      nl_flag(4,ilauf) = nlc
                      ir_flag(ilauf) = ir
                      index_flag(ir) = 1
                      ilauf = ilauf + 1
                   endif
                   if(nrb.ne.nrd.or.nlb.ne.nld) then
                      ir=rindex(nra,nrc,nrb,nrd,nla,-nlc,nlb,-nld,n)
                      if(index_flag(ir).eq.0) then
                         nr_flag(1,ilauf) = nra
                         nr_flag(2,ilauf) = nrc
                         nr_flag(3,ilauf) = nrb
                         nr_flag(4,ilauf) = nrd
                         nl_flag(1,ilauf) = nla
                         nl_flag(2,ilauf) =-nlc
                         nl_flag(3,ilauf) = nlb
                         nl_flag(4,ilauf) =-nld
                         ir_flag(ilauf) = ir
                         index_flag(ir) = 1
                         ilauf = ilauf + 1
                      endif
                   endif
                   if((nrc.ne.nrd.or.nlc.ne.nld).and.&
                        (nrc.ne.nrb.or.nlc.ne.nlb)) then
                      ir=rindex(nra,nrb,nrc,nrd,nla,nlb,nlc,nld,n)
                      if(index_flag(ir).eq.0) then
                         nr_flag(1,ilauf) = nra
                         nr_flag(2,ilauf) = nrb
                         nr_flag(3,ilauf) = nrc
                         nr_flag(4,ilauf) = nrd
                         nl_flag(1,ilauf) = nla
                         nl_flag(2,ilauf) = nlb
                         nl_flag(3,ilauf) = nlc
                         nl_flag(4,ilauf) = nld
                         ir_flag(ilauf) = ir
                         index_flag(ir) = 1
                         ilauf = ilauf + 1
                      endif
                   endif
                   if(nlb.ne.0.or.nld.ne.0) then
                      ir=rindex(nra,nrb,nrd,nrc,nla,-nlb,-nld,nlc,n)
                      if(index_flag(ir).eq.0)then
                         nr_flag(1,ilauf) = nra
                         nr_flag(2,ilauf) = nrb
                         nr_flag(3,ilauf) = nrd
                         nr_flag(4,ilauf) = nrc
                         nl_flag(1,ilauf) = nla
                         nl_flag(2,ilauf) =-nlb
                         nl_flag(3,ilauf) =-nld
                         nl_flag(4,ilauf) = nlc
                         ir_flag(ilauf) = ir
                         index_flag(ir) = 1
                         ilauf = ilauf+1
                      endif
                      if(nrb.ne.nrd .or.nlb.ne.nld) then
                         ir=rindex(nra,nrc,nrb,nrd,nla,-nlc,-nlb,nld,n)
                         if(index_flag(ir).eq.0) then
                            nr_flag(1,ilauf) = nra
                            nr_flag(2,ilauf) = nrc
                            nr_flag(3,ilauf) = nrb
                            nr_flag(4,ilauf) = nrd
                            nl_flag(1,ilauf) = nla
                            nl_flag(2,ilauf) =-nlc
                            nl_flag(3,ilauf) =-nlb
                            nl_flag(4,ilauf) = nld
                            ir_flag(ilauf) = ir
                            index_flag(ir) = 1
                            ilauf = ilauf + 1
                         endif
                      endif
                      if((nrc.ne.nrb.or.nlb.ne.nlc).and.&
                           (nrc.ne.nrd.or.nlc.ne.nld)) then
                         ir=rindex(nra,nrb,nrc,nrd,nla,-nlb,nlc,-nld,n)
                         if(index_flag(ir).eq.0) then
                            nr_flag(1,ilauf)= nra
                            nr_flag(2,ilauf)= nrb
                            nr_flag(3,ilauf)= nrc
                            nr_flag(4,ilauf)= nrd
                            nl_flag(1,ilauf)= nla
                            nl_flag(2,ilauf)=-nlb
                            nl_flag(3,ilauf)= nlc
                            nl_flag(4,ilauf)=-nld
                            ir_flag(ilauf)=ir
                            index_flag(ir)=1
                            ilauf = ilauf + 1
                         endif
                      endif
                   endif
                else
                   if(nsa.eq.nsb) then
                      ir=rindex(nra,nrb,nrd,nrc,nla,-nlb,-nld,nlc,n)
                      if(index_flag(ir).eq.0)then
                         nr_flag(1,ilauf) = nra
                         nr_flag(2,ilauf) = nrb
                         nr_flag(3,ilauf) = nrd
                         nr_flag(4,ilauf) = nrc
                         nl_flag(1,ilauf) = nla
                         nl_flag(2,ilauf) =-nlb
                         nl_flag(3,ilauf) =-nld
                         nl_flag(4,ilauf) = nlc
                         ir_flag(ilauf) = ir
                         index_flag(ir) = 1
                         ilauf = ilauf+1
                      endif
                      ir=rindex(nra,nrc,nrb,nrd,nla,-nlc,nlb,-nld,n)
                      if(index_flag(ir).eq.0) then
                         nr_flag(1,ilauf) = nra
                         nr_flag(2,ilauf) = nrc
                         nr_flag(3,ilauf) = nrb
                         nr_flag(4,ilauf) = nrd
                         nl_flag(1,ilauf) = nla
                         nl_flag(2,ilauf) =-nlc
                         nl_flag(3,ilauf) = nlb
                         nl_flag(4,ilauf) =-nld
                         ir_flag(ilauf) = ir
                         index_flag(ir) = 1
                         ilauf = ilauf + 1
                      endif
                   else
                      ir=rindex(nra,nrb,nrd,nrc,nla,nlb,nld,nlc,n)
                      if(index_flag(ir).eq.0) then
                         nr_flag(1,ilauf) = nra
                         nr_flag(2,ilauf) = nrb
                         nr_flag(3,ilauf) = nrd
                         nr_flag(4,ilauf) = nrc
                         nl_flag(1,ilauf) = nla
                         nl_flag(2,ilauf) = nlb
                         nl_flag(3,ilauf) = nld
                         nl_flag(4,ilauf) = nlc
                         ir_flag(ilauf) = ir
                         index_flag(ir) = 1
                         ilauf = ilauf + 1
                      endif
                      ir=rindex(nra,nrc,nrb,nrd,nla,-nlc,-nlb,nld,n)
                      if(index_flag(ir).eq.0) then
                         nr_flag(1,ilauf) = nra
                         nr_flag(2,ilauf) = nrc
                         nr_flag(3,ilauf) = nrb
                         nr_flag(4,ilauf) = nrd
                         nl_flag(1,ilauf) = nla
                         nl_flag(2,ilauf) =-nlc
                         nl_flag(3,ilauf) =-nlb
                         nl_flag(4,ilauf) = nld
                         ir_flag(ilauf) = ir
                         index_flag(ir) = 1
                         ilauf = ilauf + 1
                      endif
                   endif
                endif
             enddo !itd
          enddo !itb
       enddo !itc
    enddo !ita
    else
    do ita = 1,ntx
       nra = nr(ita); nza = nz(ita); nla = nl(ita); nsa = ns(ita)
       if(nza.gt.1) cycle
       iba = ib_zrls(nza,nra,nla,(nsa+1)/2)
       do itc = 1, ita
          nrc = nr(itc); nzc = nz(itc); nlc = nl(itc); nsc = ns(itc)
          if(nzc.gt.1) cycle
          ibc = ib_zrls(nzc,nrc,nlc,(nsc+1)/2)
          nsac = nsa + nsc
          if(ibc.ne.iba) cycle
          do itb = 1,nttx
             nrb=nrr(itb); nlb=nll(itb); nsb=nss(itb); ibb=noo(itb)
             do itd = 1,nttx
                nrd=nrr(itd); nld=nll(itd); nsd=nss(itd); ibd=noo(itd)
                if(ibb.ne.ibd) cycle
                if(nla+nlb.ne.nlc+nld.and.nla-nlb.ne.nlc-nld) cycle
                nsdb = nsd + nsb
                if(nsac.ne.0) then
                   ir=rindex(nra,nrb,nrc,nrd,nla,nlb,nlc,nld,n)
                   if(index_flag(ir).eq.0) then
                      nr_flag(1,ilauf) = nra
                      nr_flag(2,ilauf) = nrb
                      nr_flag(3,ilauf) = nrc
                      nr_flag(4,ilauf) = nrd
                      nl_flag(1,ilauf) = nla
                      nl_flag(2,ilauf) = nlb
                      nl_flag(3,ilauf) = nlc
                      nl_flag(4,ilauf) = nld
                      ir_flag(ilauf) = ir
                      index_flag(ir) = 1
                      ilauf = ilauf + 1
                   endif
                   if(nlb.ne.0.or.nld.ne.0) then
                      ir=rindex(nra,nrb,nrc,nrd,nla,-nlb,nlc,-nld,n)
                      if(index_flag(ir).eq.0) then
                         nr_flag(1,ilauf)= nra
                         nr_flag(2,ilauf)= nrb
                         nr_flag(3,ilauf)= nrc
                         nr_flag(4,ilauf)= nrd
                         nl_flag(1,ilauf)= nla
                         nl_flag(2,ilauf)=-nlb
                         nl_flag(3,ilauf)= nlc
                         nl_flag(4,ilauf)=-nld
                         ir_flag(ilauf)=ir
                         index_flag(ir)=1
                         ilauf = ilauf + 1
                      endif
                   endif
                endif
             enddo !itd
          enddo !itb
       enddo !itc
    enddo !ita
    endif

! The actual calculation of the matrix elements is done below using
! the flags that were set on the previous loop

!$OMP Parallel Default(None) &
!$OMP& SHARED(ilauf,nr_flag,nl_flag,ir_flag,VrGogny,n_g_all) &
!$OMP& PRIVATE(ii,nra,nrb,nrc,nrd,nla,nlb,nlc,nld,ir_abcd,ig)
!$OMP DO SCHEDULE(DYNAMIC)
    do ii = 0,ilauf-1
       nra = nr_flag(1,ii)
       nrb = nr_flag(2,ii)
       nrc = nr_flag(3,ii)
       nrd = nr_flag(4,ii)
       nla = nl_flag(1,ii)
       nlb = nl_flag(2,ii)
       nlc = nl_flag(3,ii)
       nld = nl_flag(4,ii)
       ir_abcd = ir_flag(ii)
       call radial_matrix_elements(nra,nla,nrb,nlb,nrc,nlc,nrd,nld)
       do ig = 1,n_g_all
          VrGogny(ig,ir_abcd) = Vr_ig(ig)
       enddo
    enddo
!$OMP End Do
!$OMP End Parallel
#else
    if(allocated(Vr_Gogny)) deallocate(Vr_Gogny)
    allocate(Vr_Gogny(1:n_g_all,0:n/2,-n:n,0:n/2,-n:n,0:n/2,-n:n,0:n/2,-n:n))
    Vr_gogny = zero
    do nra = 0,n/2
       do nrb = 0,n/2
          do nrc = 0,n/2
             do nrd = 0,n/2
                do nla = -n+2*nra,n-2*nra
                   do nlb = -n+2*nrb,n-2*nrb
                      do nlc = -n+2*nrc,n-2*nrc
                         do nld = -n+2*nrd,n-2*nrd
                            call radial_matrix_elements(nra,nla,nrb,nlb,nrc,nlc,nrd,nld)
                            do ig = 1,n_g_all
                               Vr_gogny(ig,nra,nla,nrb,nlb,nrc,nlc,nrd,nld) = Vr_ig(ig)
                            enddo
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
#endif
  end subroutine CalculateVrGogny


  !======================================================================

  !> Calculates all the non-zero \f$ T_{n_1 n_2}^n \f$
  !> coefficients. These coefficients expand a product of two
  !> one-dimensional harmonic oscillator wave functions (HOWF) into a
  !> linear combination of single one-dimensional HOWF's.

  !> In particular, these coefficients are calculated as
  !> \f[
  !>     T_{n_1 n_2}^n \equiv \frac{\sqrt{n_1!n_2!n!}}
  !>    {\left(\frac{n-n+n}{2}\right)!\left(\frac{n-n+n}{2}\right)!
  !>    \left(\frac{n-n+n}{2}\right)!}.
  !> \f]
  !> The coefficients are stored in the private array Tz(n1,n2,n).
  !>
  !> See \cite younes2009-b for a derivation of these coefficients
  !======================================================================
  subroutine calculateTz()
    implicit none
    integer(ipr) :: n1,n2,n,nz
    real(pr) :: n1f,n2f,nf,m1f,m2f,m3f
#if(GOGNY_HYPER==1)
    nz = max(nzx,max(2*nrx,nlx))
#else
    nz = nzx
#endif
    if(allocated(T_z)) deallocate(T_z)
    allocate(T_z(0:nz,0:nz,0:2*nz))
    T_z = zero
!$OMP Parallel Default(None) &
!$OMP& SHARED(nz,T_z) &
!$OMP& PRIVATE(n1,n2,n,n1f,n2f,nf,m1f,m2f,m3f)
!$OMP DO SCHEDULE(DYNAMIC)
    do n1 = 0, nz
       do n2 = 0, nz
          do n = 0, 2*nz
             if(n2.lt.n1) cycle
             if(n.lt.abs(n2-n1).or.n.gt.n2+n1) cycle
             if(mod(n,2).ne.mod(n2+n1,2)) cycle
             n1f = factrl(n1)
             n2f = factrl(n2)
             nf = factrl(n)
             m1f = factrl((-n1+n2+n)/2)
             m2f = factrl((+n1-n2+n)/2)
             m3f = factrl((+n1+n2-n)/2)
             T_z(n1,n2,n)=(sqrt(n2f)/m1f)*(sqrt(nf)/m2f)*(sqrt(n1f)/m3f)
             T_z(n2,n1,n)=T_z(n1,n2,n)
          enddo
       enddo
    enddo
!$OMP End Do
!$OMP End Parallel
  end subroutine calculateTz

  !======================================================================
  !> Calculates all the nonzero \f$ T_{n_1,k_1,n_2,k_2}^{n,k1+k2}\f$
  !> coefficients to transform a product of two two-dimensional harmonic
  !> oscillator wave functions in radial coordinates into a linear
  !> combination of single two dimensional HO wave functions in radial
  !> coordinates.
  !>
  !> In particular, these coefficients are calculated as
  !> \f[
  !>    T_{n_1,k_1;n_2,k_2}^{n,k_1+k_2} = (-1)^{n_1+n_2-n}
  !>    \sqrt{\frac{n!(n_1+|k_1|)!(n_2+|k_2|)!}{n_1!n_2!(n+|k_1+k_2)!}}
  !>    \sum_{m_1=0}^{n_1} \sum_{m_2=0}^{n_2} (-1)^{m_1+m_2}
  !>    {{n_1}\choose{m_1}} {{n_2}\choose{m_2}}
  !>    {{n_{1,2}-m_1-m_2}\choose{n}}\frac{(n_{1,2}+|k_1+k_2|-m_1-m_2)!}
  !>    {(n_1+|k_1|-m_1)!(n_2+|k_2|-m_2)!},
  !> \f]
  !> where \f$ n_{1,2} = n_1+n_2+\frac{|k_1|+|k_2|-|k_1+k_2|}{2}\f$.
  !> The coefficients are stored on the private array Tr(n1,k1,n2,k2,n)
  !>
  !> See \cite younes2009-b for a derivation of these coefficients
  !======================================================================
  subroutine calculateTr()
    implicit none
    integer(ipr) :: n1,n2,n,k1,k2,n12,nm12,nmax,nm2
    real(pr) :: Tpp,Tmp
    nmax = max(2*nrx,nlx)
    nm2 = nmax/2
    if(allocated(T_r)) deallocate(T_r)
    allocate(T_r(0:nm2,-nmax:nmax,0:nm2,-nmax:nmax,0:nmax))
    do n1 = 0, nm2
       do k1 = 0,nmax-2*n1
          do n2 = n1, nm2
             do k2 = 0,nmax-2*n2
                nm12 = n1 + n2 + (k1+k2-abs(-k1+k2))/2
                do n = 0, nm12
                   Tpp = TrCoefficient(n1, k1,n2,k2,n, k1+k2)
                   Tmp = TrCoefficient(n1,-k1,n2,k2,n,-k1+k2)
                   T_r(n1, k1,n2, k2,n) = Tpp
                   T_r(n2, k2,n1, k1,n) = Tpp
                   T_r(n1,-k1,n2,-k2,n) = Tpp
                   T_r(n2,-k2,n1,-k1,n) = Tpp
                   T_r(n1,-k1,n2, k2,n) = Tmp
                   T_r(n2, k2,n1,-k1,n) = Tmp
                   T_r(n1, k1,n2,-k2,n) = Tmp
                   T_r(n2,-k2,n1, k1,n) = Tmp
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine calculateTr

  !======================================================================
  !> Calculates the \f$ T_{n_1,k_1,n_2,k_2}^{n,k} \f$  coefficient.
  !> Since \f$ T_{n_1,k_1,n_2,k_2}^{n,k}=T_{n,k,n_2,k_2}^{n_1,k_1}=
  !> T_{n_1,k_1,n,k}^{n_2,k_2} \f$, the combination with the least number
  !> of operations is calculated in order to reduce loss of accuracy.
  !> trsumterms() is used to calculate the number of operations for each
  !> combination
  !>
  !> See Eqs. (C.9) and (C.10) of \cite younes2009-b for details
  !>
  !> @result \f$ T_{n_1,k_1,n_2,k_2}^{n,k} = (-1)^{n_1+n_2+n}
  !>         \sqrt{\frac{n!(n_1+|k_1|)!(n_2+|k_2|)!}{n_1!n_2!(n+|k|)!}}
  !>         \sum_{m_1=0}^{n_1} \sum_{m_2=0}^{n_2} (-1)^{m_1+m_2}
  !>         {{n_1}\choose{m_1}} {{n_2}\choose{m_2}}
  !>         {{n_{1,2}-m_1-m_2}\choose{n}}
  !>         \frac{(n_{1,2}+|k|-m_1-m_2)!}{(n_1+|k_1|-m_1)!
  !>         (n_2+|k_2|-m_2)!}\f$, where \f$ n_{1,2} = n_1+n_2+
  !>         \frac{|k_1|+|k_2|-|k|}{2} \f$.
  !======================================================================
  function TrCoefficient(n1,k1,n2,k2,n,k) result(Tr)
    implicit none
    integer(ipr), intent(in) :: n1 !< first  radial principal quantum number
    integer(ipr), intent(in) :: n2 !< second radial principal quantum number
    integer(ipr), intent(in) :: n  !< third  radial principal quantum number
    integer(ipr), intent(in) :: k1 !< first  radial orbital quantum number
    integer(ipr), intent(in) :: k2 !< second radial orbital quantum number
    integer(ipr), intent(in) :: k  !< third  radial orbital quantum number
    real(pr) :: Tr
    integer(ipr), dimension(3) :: Nsums
    integer(ipr) minp,in1,in2,in,ik1,ik2,ik
    integer(ipr):: n12,m1,m2
    real(pr) :: n1f,nk1f,n2f,nk2f,nf,nkf,fac,Ti
    if(mod(k1+k2+k,2).ne.0) then
       stop 'wrong k values in TrCoefficient'
    endif
    Nsums(1) = TrSumTerms(n ,k ,n2,k2,n1,k1)
    Nsums(2) = TrSumTerms(n1,k1,n ,k ,n2,k2)
    NSums(3) = TrSumTerms(n1,k1,n2,k2,n ,k )
    minp = minloc(Nsums,1)
    select case(minp)
    case(1)
       in1 = n
       ik1 = k
       in2 = n2
       ik2 = k2
       in  = n1
       ik  = k1
    case(2)
       in1 = n1
       ik1 = k1
       in2 = n
       ik2 = k
       in  = n2
       ik  = k2
    case(3)
       in1 = n1
       ik1 = k1
       in2 = n2
       ik2 = k2
       in  = n
       ik  = k
    end select
    n12 = in1 + in2 + (abs(ik1)+abs(ik2)-abs(ik))/2
    n1f = factrl(in1)
    nk1f = factrl(in1 + abs(ik1))
    n2f = factrl(in2)
    nk2f = factrl(in2 + abs(ik2))
    nf = factrl(in)
    nkf = factrl(in+abs(ik))
    Tr = 0
    do m1 = 0,min(in1,n12-in)
       do m2 = 0,min(in2,n12-in-m1)
          Ti = (-1)**(m1+m2)*binomialco(in1,m1)&
               *binomialco(in2,m2)*binomialco(n12-m1-m2,in)&
               *(factrl(n12+abs(ik)-m1-m2)&
               /factrl(in1+abs(ik1)-m1))&
               /factrl(in2+abs(ik2)-m2)
          Tr = Tr + Ti
       enddo
    enddo
    Fac = (-1)**(in1+in2-in)*sqrt(nf/nkf)*sqrt(nk1f/n1f)&
         *sqrt(nk2f/n2f)
    Tr = Fac*Tr
  end function TrCoefficient

  !======================================================================
  !> Calculates the number of sumations required to calculate the
  !> \f$ T_{n_1,k_1,n_2,k_2}^{n,k} \f$  coefficient.
  !======================================================================
  function TrSumTerms(n1,k1,n2,k2,n,k) result(NT)
    integer(ipr), intent(in) :: n1 !< first  radial principal quantum number
    integer(ipr), intent(in) :: n2 !< second radial principal quantum number
    integer(ipr), intent(in) :: n  !< third  radial principal quantum number
    integer(ipr), intent(in) :: k1 !< first  radial orbital quantum number
    integer(ipr), intent(in) :: k2 !< second radial orbital quantum number
    integer(ipr), intent(in) :: k  !< third  radial orbital quantum number
    integer(ipr) :: NT
    integer(ipr) :: n12,mu,ab_min,ab_max,a,b
    n12 = n1 + n2 + (abs(k1)+abs(k2)-abs(k))/2
    mu = n12-n
    if(mu.lt.0) then
       NT = 0
       return
    endif
    ab_min = min(n1,n2)
    ab_max = max(n1,n2)
    if(mu.le.ab_min) then
       NT = (mu+1)*(mu+2)/2
    elseif(mu.le.ab_max) then
       NT = (mu+1)*(mu+2)/2 - (mu-ab_min)*(mu-ab_min+1)/2
    else
       if(mu.le.n1+n2) then
          NT = (mu+1)*(mu+2)/2 - (mu-ab_min)*(mu-ab_min+1)/2 &
               - (mu-ab_max)*(mu-ab_max+1)/2
       else
          NT = (n1+1)*(n2+1)
       endif
    endif
  end function TrSumTerms


  !======================================================================
  !>  Calculates all the necessary \f$C_{n_x n_y}^{n k}i^{n_y}\f$
  !>  coefficients with \f$n_x+n_y=2n+|k| \f$ to transform a two-
  !>  dimensional harmonic oscillator wave function (HOWF) in radial
  !>  coordinates into a linear combination of products of
  !>  one-dimensional HOWF in cartesian coordinates.
  !>
  !>  The coefficients are calculated with
  !>  \f[
  !>      C_{n_x n_y}^{n k}i^{n_y} = \frac{2^{-n-|k|/2}(-1)^n
  !>      \sqrt{n!(n+|k|)!}}{\sqrt{(2n+|k|-n_y)!n_y!}}
  !>      \sum_{q=0}^{\min(n_y,n+(|k|-k)/2)}
  !>      {{2n+|k|-n_y}\choose{n-q+\frac{|k|-k}{2}}}
  !>      {{n_y}\choose{q}} (-1)^{n_y-q}
  !>  \f]
  !>  The \f$ i^{n_y}\f$ factor makes the result a real number.
  !>
  !>  These transformation coefficients are only calculated when the
  !>  preprocessor directive GOGNY_HYPER is set to 1. The coefficients
  !>  are stored in the private array Cp2c(nr,k,ny).
  !>
  !>  See \cite younes2009-b for a derivation and details.
  !======================================================================
  subroutine calculateCpolar2cartesian()
    implicit none
    integer(ipr) :: n,k,ny,nx,q,qmax,nshells
    real(pr) :: A,xsum
    nshells = max(2*nrx,nlx)
    if(allocated(Cp2c)) deallocate(Cp2c)
    allocate(cp2c(0:Nshells,-2*Nshells:2*Nshells,0:2*Nshells))
    cp2c = zero
!$OMP Parallel Default(None) &
!$OMP& SHARED(Nshells,Cp2c) &
!$OMP& PRIVATE(n,k,ny,nx,A,xsum,qmax,q)
!$OMP DO SCHEDULE(DYNAMIC)
    do n = 0,Nshells
       do k = -2*nshells,2*nshells
          do ny = 0,2*nshells
             if(abs(k).gt.2*nshells-2*n) cycle
             if(ny.gt.2*n+abs(k)) cycle
             nx = 2*n+abs(k)-ny
             A = (-1)**n*(2**(-n-abs(k)*0.5_pr))*sqrt(factrl(n+abs(k)))&
                  *(sqrt(factrl(n))/(sqrt(factrl(nx))*sqrt(factrl(ny))))
             xsum = zero
             qmax = min(ny,n+(abs(k)-k)/2)
             do q = 0,qmax
                xsum = xsum + BinomialCo(nx,n-q+(abs(k)-k)/2)*BinomialCo(ny,q)*(-1)**(ny-q)
             enddo
             Cp2c(n,k,ny) = A*xsum
          enddo
       enddo
    enddo
!$OMP End Do
!$OMP End Parallel
  end subroutine calculateCpolar2cartesian

  !======================================================================
  !>  Calculates all the non-zero axial two-body matrix elements
  !>  \f$ \langle n_1 n_2|V|n_3 n_4 \rangle \f$ using the radial
  !>  oscillator parameter \f$b_{\perp}\f$.
  !>
  !>  These matrix elements are only calculated when the preprocessor
  !>  directive GOGNY_HYPER is set to 1. The matrix elements are stored
  !>  in the private array ME1D(n1,n2,n3,n4).
  !======================================================================
  subroutine calculateME1D()
    implicit none
    integer(ipr) :: ni,nj,nk,nl,n,ig
    real(pr) :: ME
    n = max(2*nrx,nlx)
    if(allocated(ME1D)) deallocate(ME1D)
    allocate(ME1D(1:n_g_all,0:n,0:n,0:n,0:n))
    ME1D = zero
!$OMP Parallel Default(None) &
!$OMP& SHARED(n,ME1D,mu_g_all,bp,n_g_all) &
!$OMP& PRIVATE(ni,nj,nk,nl,ig,ME)
!$OMP DO SCHEDULE(DYNAMIC)
    do ni = 0,n
       do nj = 0,n
          do nk = 0,n
             do nl = 0,n
                if(nk.lt.ni) cycle
                if(nl.lt.nj) cycle
                if(mod(ni+nj+nk+nl,2).ne.0) cycle
                do ig = 1,n_g_all
                   ME = MatrixElement_z(ni,nj,nk,nl,mu_g_all(ig),bp)
                   ME1D(ig,ni,nj,nk,nl) = ME
                   ME1D(ig,nk,nj,ni,nl) = ME
                   ME1D(ig,ni,nl,nk,nj) = ME
                   ME1D(ig,nk,nl,ni,nj) = ME
                enddo
             enddo
          enddo
       enddo
    enddo
!$OMP End Do
!$OMP End Parallel
  end subroutine calculateME1D

  !======================================================================
  !> Calculates the  axial component of the two body potential matrix
  !> element for a Gaussian potential with axial symmetry.
  !>
  !> If the preproccesor  variable GOGNY_HYPER is set to 1 the matrix
  !> element is calculated using the modified Gogny expansion that uses a
  !> hypergeometric function and preserves accuracy with a large basis
  !> size. See Eq. (D.12) in \cite younes2009-b.
  !> \f[
  !>    \langle n_i n_j| \hat{V}_z|n_k n_l \rangle
  !>    = \frac{\mu}{\sqrt{2\pi^3}b_z}\sum_{n=|n_j-n_l|,2}^{n_j+n_l}
  !>      T_{n_j n_l}^{n}\bar{F}_{n_i n_k}^{n},
  !> \f]
  !> where
  !> \f[
  !>    \bar{F}_{n_i n_k}^{n} = \frac{\Gamma(\xi-n_i)\Gamma(\xi-n_k)
  !>      \Gamma(\xi-n)}{(1+\mu^2/(2b_z^2))^\xi\sqrt{n!n_i!n_k!}}
  !>      \ _2F_1(-n_i,n_k;-\xi+n+1;-\mu^2/(2b_z^2))
  !> \f]
  !> with \f$ \xi = \frac{n_i+n_k+n+1}{2} \f$.
  !>
  !> For any other value of GOGNY_HYPER the matrix element is calculated
  !> using the direct Gogny transformation which looses numerical
  !> accuracy with a large basis size. See Eq. (10) in
  !> \cite younes2009-b.
  !> \f[
  !>    \langle n_i,n_j|V_z|n_k,n_l \rangle = \sqrt{\frac{G_z-1}{G_z+1}}
  !>    \sum_{m=|n_i-n_k|,2}^{n_i+n_k}\sum_{n=|n_j-n_l|,2}^{n_j+n_l}
  !>    T_{n_i,n_k}^{m} T_{n_j,n_l}^{n} \bar{I}(m,n)
  !> \f]
  !======================================================================
  function MatrixElement_z(ni,nj,nk,nl,mu,b) result(Vz)
    implicit none
    integer(ipr), intent(in) :: ni !< first  z-component quantum number
    integer(ipr), intent(in) :: nj !< second z-component quantum number
    integer(ipr), intent(in) :: nk !< third  z-component quantum number
    integer(ipr), intent(in) :: nl !< fourth z-component quantum number
    real(pr), intent(in) :: mu,b
    real(pr) :: Vz
    real(pr) :: Gz
    integer(ipr) :: mz,nz
    integer(ipr), dimension(2) :: Nsums
    integer(ipr) :: minp,ini,inj,ink,inl
    real(pr) Fbar,xi,z
    Vz = 0._pr
#if(GOGNY_HYPER==1)
    Nsums(1) = (nj+nl-abs(nj-nl))/2 + 1
    Nsums(2) = (ni+nk-abs(ni-nk))/2 + 1
    if(nsums(1).eq.nsums(2)) then
       nsums(1) = min(ni,nk)
       nsums(2) = min(nj,nl)
    endif
    minp = minloc(Nsums,1)
    select case(minp)
    case(1)
       ini = ni
       inj = nj
       ink = nk
       inl = nl
    case(2)
       ini = nj
       inj = ni
       ink = nl
       inl = nk
    end select
    z = one+mu**2/(two*b**2)
    Vz = zero
    do nz = abs(inj-inl),inj+inl,2
       if(mod(ini+ink+nz,2).ne.0) exit
       xi = (ini+ink+nz+1)*0.5_pr
       Fbar = gamma(xi-ini)*gamma(xi-ink)*gamma(xi-nz)&
            /(z**xi*sqrt(factrl(ini)*factrl(ink)*factrl(nz)))&
            *HyperGeom2F1(-ini,-ink,nz+1-xi,1-z)
       Vz = Vz + T_z(inj,inl,nz)*Fbar
    enddo
    Vz = mu/(sqrt(two*pi**3)*b)*Vz
#else
    Gz = one + (mu/b)**2
    do mz = abs(ni-nk),ni+nk,2
       do nz = abs(nj-nl),nj+nl,2
          Vz = Vz + T_z(ni,nk,mz)*T_z(nj,nl,nz)*Ibarz(mz,nz,Gz)
       enddo
    enddo
    Vz = sqrt((Gz-1)/(Gz+1))*Vz
#endif
  end function MatrixElement_z

  !======================================================================
  !> Calculates the \f$ \bar{I}(m,n) \f$ coefficient necessary to
  !> calculate the axial component of the two body potential matrix
  !> element for a Gaussian potential with axial symmetry.
  !> See \cite younes2009-b
  !>
  !> @result \f$ \bar{I}(m,n) = \sqrt{\frac{m!n!}{2^{m+n}}} \frac{(-1)
  !> ^{(m-n)/2}}{\left(\frac{m+n}{2}\right)!(1+G_z)^{(m+n)/2}}
  !> {{m+n}\choose{n}}\f$
  !======================================================================
  function Ibarz(m,n,Gz) result(Iz)
    implicit none
    integer(ipr), intent(in) :: m !< An integer
    integer(ipr), intent(in) :: n !< An integer
    real(pr), intent(in) :: Gz
    real(pr) :: Iz
    if(mod(m+n,2).ne.0) then
       Iz =0._pr
       return
    endif
    Iz = (-1)**((m-n)/2)*(sqrt(factrl(m))/factrl((m+n)/2))*&
         sqrt(factrl(n)/(2+2*Gz)**(m+n))*binomialco(m+n,n)
  end function Ibarz



  !======================================================================
  !>  Calculates the radial component of the two-body potential matrix
  !>  element for a given set of radial quantum numbers
  !>
  !>  The matrix elements are calculated for all the Gaussians in the
  !>  finite range functional since each Gaussian has its own range
  !>  parameter in the mu_g_all(:) array.
  !>
  !>  If the preprocessor variable GOGNY_HYPER is set to 1 the matrix
  !>  elements are obtained by transforming the HO wavefunctions
  !>  from radial into cartesian coordinates and separating the
  !>  integration into a product of two one-dimensional two body matrix
  !>  elements. The transformation is given by
  !>  \f[ \langle n_{r_i} \Lambda_i,n_{r_j}\Lambda_j|\hat{V}_p|
  !>      n_{r_k} \Lambda_k,n_{r_l}\Lambda_l \rangle =
  !>      \sum_{n_{y_i}=0}^{2n_{r_i}+|\Lambda_i|}
  !>      \sum_{n_{y_j}=0}^{2n_{r_j}+|\Lambda_j|}
  !>      \sum_{n_{y_k}=0}^{2n_{r_k}+|\Lambda_k|}
  !>      \sum_{n_{y_l}=0}^{2n_{r_l}+|\Lambda_l|}
  !>      C_{n_{x_i} n_{y_i}}^{n_{r_i}\Lambda_i *}
  !>      C_{n_{x_j} n_{y_j}}^{n_{r_j}\Lambda_j *}
  !>      C_{n_{x_k}n_{y_k}}^{n_{r_k}\Lambda_k}
  !>      C_{n_{x_l}n_{y_l}}^{n_{r_l}\Lambda_l}
  !>      \langle n_{x_i}n_{x_j}|\hat{V}_{\rm 1D}|n_{x_k}n_{x_l}\rangle
  !>      \langle n_{y_i}n_{y_j}|\hat{V}_{\rm 1D}|n_{y_k}n_{y_l}\rangle
  !>  \f]
  !>  The one dimentional matrix elements have already been calculated by
  !>  the calculateme1d() subroutine and stored in the ME1D array. Since
  !>  calculateme1d() uses the modified gogny transformation appropriate
  !>  for large basis size, this transformation is also used when
  !>  preserving accuracy with large basis.
  !>
  !>  For any other value of GOGNY_HYPER the function matrixelement_r(),
  !>  which uses the direct Gogny transformation, is used.
  !======================================================================
  subroutine radial_matrix_elements(ni,li,nj,lj,nk,lk,nl,ll)
    implicit none
    integer(ipr),intent(in) :: ni !<first  radial principal quantum number
    integer(ipr),intent(in) :: nj !<second radial principal quantum number
    integer(ipr),intent(in) :: nk !<third  radial principal quantum number
    integer(ipr),intent(in) :: nl !<fourth radial principal quantum number
    integer(ipr),intent(in) :: li !<first  radial orbital quantum number
    integer(ipr),intent(in) :: lj !<second radial orbital quantum number
    integer(ipr),intent(in) :: lk !<third  radial orbital quantum number
    integer(ipr),intent(in) :: ll !<fourth radial orbital quantum number
    real(pr) :: Ci,Cj,Ck,Cl
    integer(ipr) :: nyi,nyj,nyk,nyl,nxi,nxj,nxk,nxl,ig
    if(.not.allocated(Vr_ig)) allocate(Vr_ig(1:n_g_all))
    Vr_ig = zero
    if(-li-lj+lk+ll.ne.0) return
#if(GOGNY_HYPER==1)
    do nyi = 0,2*ni+abs(li)
       nxi = 2*ni+abs(li) - nyi
       Ci = Cp2c(ni,li,nyi)
       do nyj = 0,2*nj+abs(lj)
          nxj = 2*nj+abs(lj)-nyj
          Cj = Cp2c(nj,lj,nyj)
          do nyk = 0,2*nk+abs(lk)
             nxk = 2*nk+abs(lk) - nyk
             Ck = Cp2c(nk,lk,nyk)
             do nyl = mod(nyi+nyj+nyk,2),2*nl+abs(ll),2
                nxl = 2*nl+abs(ll)-nyl
                Cl = Cp2c(nl,ll,nyl)
                do ig = 1,n_g_all
                   Vr_ig(ig) = Vr_ig(ig) + &
                        (-1)**(nyi+nyj+(nyi+nyj+nyk+nyl)/2)*Ci*Cj*Ck*Cl &
                        *ME1D(ig,nxi,nxj,nxk,nxl)*&
                        ME1D(ig,nyi,nyj,nyk,nyl)
                enddo
             enddo
          enddo
       enddo
    enddo
#else
    do ig = 1,n_g_all
       Vr_ig(ig) = MatrixElement_r(ni,li,nj,lj,nk,lk,nl,ll,mu_g_all(ig),bp)
    enddo
#endif
  end subroutine radial_matrix_elements


  !======================================================================
  !> Calculates the radial component of the two body potential matrix
  !> element for a Gaussian potential with axial symmetry using the
  !> direct Gogny transformation.
  !>
  !> This transformation looses accuracy as the basis size increases
  !> and should not be used in production runs.
  !>
  !> See \cite younes2009-b
  !>
  !> @result \f[ <n_{r_i},l_i,n_{r_j},l_j|V_p|n_{r_k},l_k,n_{r_l},l_l> =
  !>         \frac{G_p-1}{G_p+1} \sum_{n_r=0}^{n_{\bar{j},l}}
  !>         \sum_{n=0}^{n_{\bar{i},k}}
  !>         T_{n_{r_i},-l_i,n_{r_k},l_k}^{n,-l_i+l_k}
  !>         T_{n_{r_j},-l_j,n_{r_l},l_l}^{n_r,-l_j+l_l}
  !>         \bar{I}(n_r,-l_j+l_l,n,-l_i+l_k)
  !>         \delta_{-l_j+l_l,l_i+l_k} \f]
  !======================================================================
  function MatrixElement_r(ni,li,nj,lj,nk,lk,nl,ll,mu,b) result(Vr)
    implicit none
    integer(ipr),intent(in) :: ni!< first  radial principal quantum number
    integer(ipr),intent(in) :: nj!< second radial principal quantum number
    integer(ipr),intent(in) :: nk!< third  radial principal quantum number
    integer(ipr),intent(in) :: nl!< fourth radial principal quantum number
    integer(ipr),intent(in) :: li!< first  radial orbital quantum number
    integer(ipr),intent(in) :: lj!< second radial orbital quantum number
    integer(ipr),intent(in) :: lk!< third  radial orbital quantum number
    integer(ipr),intent(in) :: ll!< fourth radial orbital quantum number
    real(pr), intent(in) :: mu!< range of the Gaussian potential
    real(pr), intent(in) :: b!< Harmonic Oscillator length
    real(pr) :: Vr
    real(pr) :: Gp
    integer(ipr) :: mr,nr,nmjl,nmik
    Vr = 0._pr
    if(-li-lj+lk+ll.ne.0) return
    Gp = one + (mu/b)**2
    nmjl = nj+nl+(abs(-lj)+abs(ll)-abs(-lj+ll))/2
    nmik = ni+nk+(abs(-li)+abs(lk)-abs(-li+lk))/2
    do mr = 0,nmik
       do nr = 0,nmjl
          Vr = Vr + T_r(ni,-li,nk,lk,mr)*T_r(nj,-lj,nl,ll,nr)&
               *Ibarr(nr,-lj+ll,mr,-li+lk,Gp)
       enddo
    enddo
    Vr = (Gp-1)/(Gp+1)*Vr
  end function MatrixElement_r

  !======================================================================
  !> Calculates the \f$ \bar{I}(n_1,k_1,n_2,k_2) \f$ coefficient
  !> necessary to calculate the radial component of the two body
  !> potential matrix element component for a Gaussian potential with
  !> axial symmetry
  !>
  !> See \cite younes2009-b
  !>
  !> @result \f$ \bar{I}(n_1,k_1,n_2,k_2) = \sqrt{\frac{(n_1+|k_1|)!n_2!}
  !>         {n_1!(n_2+|k_2|)!}} {{n_1+n_2+|k_1|}\choose{n_2}}
  !>         \frac{\delta_{k_1,-k_2}}{(G_p+1)^{n_1+n_2+|k_1|}}\f$
  !======================================================================
  function Ibarr(n1,k1,n2,k2,Gp) result(Ir)
    implicit none
    integer(ipr),intent(in):: n1 !< first  radial principal quantum number
    integer(ipr),intent(in):: n2 !< second radial principal quantum number
    integer(ipr),intent(in):: k1 !< first  radial orbital quantum number
    integer(ipr),intent(in):: k2 !< second radial orbital quantum number
    real(pr), intent(in) :: Gp
    real(pr) :: Ir
    integer(ipr) :: k
    real(pr) :: n1f,n1kf,n2f,n2kf,n12kf
    if(k1.ne.-k2) then
       Ir =0._pr
       return
    endif
    k = abs(k1)
    n1f = factrl(n1)
    n2f = factrl(n2)
    n1kf = factrl(n1+k)
    n2kf = factrl(n2+k)
    Ir = sqrt(n1kf/n1f)*sqrt(n2f/n2kf)*(binomialco(n1+n2+k,n2)&
         /((Gp+1)**(n1+n2+k)))
  end function Ibarr


  !======================================================================
  !>  Recursively calculates and stores in an array (to avoid
  !>  recalculating in future calls) \f$ 0!, 1!, 2!, \ldots, n! \f$.
  !>
  !>  If \f$n > 170\f$, \f$ n! = \Gamma(n+1) \f$ is used although it will
  !>  overflow on most machines.
  !>
  !>  @result  \f$ n! \f$
  !======================================================================
  function factrl(n) result(fact)
    implicit none
    integer(ipr), intent(in) :: n !< An integer
    real(pr) :: fact
    integer(ipr), save :: ntop = 0
    integer(ipr), parameter :: nmax = 170
    integer(ipr) :: i
    real(pr), dimension(0:nmax), save :: a=0._pr
    if(n.lt.0) then
       write(*,*) 'negative integer in factrl'
       stop
    endif
    if(ntop.eq.0) a(0)=one
    if(n.le.ntop) then
       fact = a(n)
    elseif(n.le.nmax) then
       do i = ntop+1,n
          a(i) = real(i,kind=pr)*a(i-1)
       enddo
       ntop = n
       fact = a(n)
    else
       fact = exp(log_gamma(n+1._pr))
    endif
  end function factrl

  !======================================================================
  !>  Calculates the Hypergeometric function \f$_2F_1(a,b,c;x) \f$ in the
  !>  very particular case where \f$ a \leq 0\f$, \f$ b \leq 0\f$ and
  !>  c is NOT a negative integer such that \f$ |c| \leq \min(|a|,|b|)\f$
  !>
  !>  @result
  !>      \f[
  !>          _2F_1(a,b,c;x) = \sum_{i=0}^d {{d}\choose{i}}{{e}
  !>           \choose{i}} \frac{i!}{(c)_i}x^i
  !>      \f]
  !>      where \f$ d = \min(|a|,|b|)\f$, \f$ e = \max(|a|,|b|)\f$ and
  !>      \f$ (x)_i \f$ is the rising Pochhammer Symbol (also known as
  !>      upper factorial)
  !======================================================================
  function HyperGeom2F1(a,b,c,x) result(HG2F1)
    implicit none
    integer(ipr), intent(in) :: a !<First parameter, has to be negative
    integer(ipr), intent(in) :: b !<Second parameter, has to be negative
    real(pr), intent(in) :: c !< Third parameter, cannot be a negative integer such that \f$ |c| \leq \min(|a|,|b|)\f$
    real(pr), intent(in) :: x !< Value where the function is evaluated
    real(pr) :: HG2F1
    integer(ipr) :: d,e,i
    real(pr) :: HGi
    if(a.gt.0.or.b.gt.0) then
       write(*,*) 'a or b is not a negative integer in HyperGeom2F1'
       stop
    endif
    d = abs(max(a,b))
    e = abs(min(a,b))
    if(c.lt.0._pr.and.abs(c-nint(c)).lt.1.e-14_pr.and.abs(c).le.d) then
       write(*,*) 'c is a negative integer greater or equal than max(a,b)in HyperGeom2F1'
       stop
    endif
    HG2F1 = zero
    do i = 0,d
       HGi = + binomialco(d,i)*binomialco(e,i)*x**i*factrl(i)/&
            upperfactrl(c,i)
       HG2F1 = HG2F1 + HGi
    enddo
  end function HyperGeom2F1

  !======================================================================
  !>  Calculates the binomial coefficient \f$ {{m}\choose{n}} \f$.
  !>
  !>  Overflows are avoided by changing the factorials for logarithms of
  !>  the gamma function and taking the exponential
  !>
  !>  @result  \f$ \displaystyle {{m}\choose{n}}=\frac{m!}{n!(m-n)!)}\f$
  !======================================================================
  function binomialco(m,n) result(bc)
    implicit none
    integer(ipr), intent(in):: m !< A positive integer
    integer(ipr), intent(in):: n !< An integer
    real(pr) :: bc
    if(m.lt.0) then
       write(*,*) 'negative integer in binomialco'
       stop
    endif
    if(n.lt.0.or.n.gt.m) then
       bc = zero
       return
    endif
    if(n.eq.0.or.n.eq.m) then
       bc = one
       return
    endif
    if(n.eq.1.or.n.eq.m-1) then
       bc = real(m,kind=pr)
       return
    endif
    if(m.le.170) then
       if(n.le.m/2) then
          bc = (factrl(m)/factrl(m-n))/factrl(n)
       else
          bc = (factrl(m)/factrl(n))/factrl(m-n)
       endif
    else
       bc = exp(log_gamma(m+1._pr)-log_gamma(n+1._pr)-&
            log_gamma(m-n+1._pr))
    endif
  end function binomialco


  !======================================================================
  !>  Calculates the upper factorial, also knwon as the rising Pochhamer
  !>  symbol \f$ (x)_i \f$
  !>  @result  \f$ (x)_i =  x(x+1)\ldots (x+i-1)\f$
  !======================================================================
  function upperfactrl(x,i) result(upf)
    implicit none
    real(pr), intent(in) :: x !< Value where the upper factorial is evaluated
    integer(ipr), intent(in) :: i !< An integer
    real(pr) :: upf
    integer(ipr) :: j
    if(i.lt.0) then
       write(*,*) 'negative integer in upperfactrl'
       stop
    endif
    upf = one
    if(i.eq.0) return
    do j = 0,i-1
       upf = upf*(x+Real(j,kind=pr))
    enddo
  end function upperfactrl

  !======================================================================
  !>  Calculates and stores in memory what is called the ZBlock. ZBlock
  !>  is a (N x N x N x N) array that is used to obtain the
  !>  transformation from the axial component quantum numbers
  !>  (\f$n_z\f$'s) into the index of an array that contains only the
  !>  the non-zero axial two body potential matrix elements.
  !>
  !>  In particular the ZBlock is
  !>  \f{eqnarray*}{
  !>   NB &=&\sum_{i=0}^{n_{z_i}-1} \sum_{j=0}^{N} \sum_{k=i}^{N}
  !>      \sum_{\substack{l=j\\ i+j+k+l\ {\rm is\ even}}}^{N} 1
  !>     +\sum_{j=0}^{n_{z_j}-1} \sum_{k=n_{z_i}}^{N}
  !>      \sum_{\substack{l=j\\ n_{z_i}+j+k+l\ {\rm is\ even}}}^{N} 1 \\
  !>      & + & \sum_{k=n_{z_i}}^{n_{z_k}-1}
  !>            \sum_{\substack{l=n_{z_j}\\ n_{z_i}+n_{z_j}+k+
  !>                  l\ {\rm is\ even}}}^{N} 1
  !>           +\sum_{\substack{l=n_{z_j}\\ n_{z_i}+n_{z_j}+n_{z_k}
  !>                  +l\ {\rm is\ even}}}^{n_{z_l}-1} 1  + 1
  !>  \f}
  !>  where \f$N\f$ is the number of shells
  !======================================================================
  subroutine calculate_Zblock()
    implicit none
    integer(ipr) :: sni,snj,snk,snl
    integer(ipr) :: xni,xnj,xnk,xnl
    integer(ipr) :: ii,jj,kk,ll,ix,n
    n = nzx
    if(allocated(ZBlock)) deallocate(ZBlock)
    allocate(ZBlock(0:n,0:n,0:n,0:n))
    Zblock = 0
    xni = 0
    do ii = 0,n
       sni = xni
       xnj = 0
       do jj = 0,n
          snj = xnj
          xnk = 0
          do kk = ii,n
             snk = xnk
             xnl = 0
             do ll = jj+mod(ii+kk,2),n,2
                snl = xnl
                ZBlock(ii,jj,kk,ll)=sni+snj+snk+snl+1
                xni = xni + 1
                xnj = xnj + 1
                xnk = xnk + 1
                xnl = xnl + 1
             enddo
          enddo
       enddo
    enddo
  end subroutine calculate_Zblock

  !======================================================================
  !>  Calculates and stores in memory what is called the NBlock. NBlock
  !>  is an (N/2 x N/2 x N/2 x N/2 x 2N x 2N) array that is used to
  !>  obtain the transformation from the radial component quantum numbers
  !>  (\f$n\f$'s and \f$\Lambda\f$'s) into the index of an array that
  !>  contains only the non-zero two-body potential matrix elements.
  !>
  !>  The Nblock is the part of that index that only depends on the
  !>  \f$n\f$ quantum numbers and the first two \f$\Lambda\f$ quantum
  !>  number. This part takes the most time to calculate (and
  !>  is therefore stored in memory so that it's only calculated once).
  !>  The other part of that index is calculated with l_block().
  !>
  !>  In particular the NBlock is
  !>  \f{eqnarray*}{
  !>   NB &=& \sum_{i=0}^{n_{ri} -1} \sum_{j=0}^{N/2} \sum_{k=0}^{N/2}
  !>          \sum_{l=0}^{N/2} \sum_{p=-N+2i}^{N-2i}\;
  !>          \sum_{q=-N+2j}^{N-2j}\;
  !>          \sum_{r=\max(-N+2k,-\min(N-l,-q)+p+q)}
  !>              ^{\min(\min(N-2k,-p),N-2l+p+q)} 1 \\
  !>      &+& \sum_{j=0}^{n_{rj}-1} \sum_{k=0}^{N/2}
  !>          \sum_{l=0}^{N/2} \sum_{p=-N+2n_{ri}}^{N-2n_{ri}}\;
  !>          \sum_{q=-N+2j}^{N-2j}\;
  !>          \sum_{r=\max(-N+2k,-\min(N-l,-q)+p+q)}
  !>              ^{\min(\min(N-2k,-p),N-2l+p+q)} 1 \\
  !>      &+& \sum_{k=0}^{n_{rk}-1} \sum_{l=0}^{N/2}
  !>          \sum_{p=-N+2n_{ri}}^{N-2n_{ri}}\;
  !>          \sum_{q=-N+2n_{rj}}^{N-2n_{rj}}\;
  !>          \sum_{r=\max(-N+2k,-\min(N-l,-q)+p+q)}
  !>              ^{\min(\min(N-2k,-p),N-2l+p+q)} 1 \\
  !>      &+& \sum_{l=0}^{n_{rl}-1} \sum_{p=-N+2n_{ri}}^{N-2n_{ri}}\;
  !>          \sum_{q=-N+2n_{rj}}^{N-2n_{rj}}\;
  !>          \sum_{r=\max(-N+2n_{rk},-\min(N-l,-q)+p+q)}
  !>              ^{\min(\min(N-2k,-p),N-2l+p+q)} 1 \\
  !>      &+& \sum_{p=-N+2n_{ri}}^{N-2n_{ri}}\;
  !>          \sum_{q=-N+2n_{rj}}^{N-2n_{rj}}\;
  !>          \sum_{r=\max(-N+2n_{rk},-\min(N-n_{rl},-q)+p+q)}
  !>              ^{\min(\min(N-2k,-p),N-2l+p+q)} 1 \\
  !>      &+& \sum_{q=-N+2n_{rj}}^{N-2n_{rj}}\;
  !>          \sum_{r=\max(-N+2n_{rk},-\min(N-n_{rl},-q)+\Lambda_i+q)}
  !>              ^{\min(\min(N-2k,-\Lambda_i),N-2l+\Lambda_i+q)} 1
  !>  \f}
  !>  where \f$N\f$ is the number of shells
  !======================================================================
  subroutine calculate_Nblock()
    implicit none
    integer(ipr) :: sni,snj,snk,snl,sli,slj
    integer(ipr) :: xni,xnj,xnk,xnl,xli,xlj
    integer(ipr) :: ii,jj,kk,ll,i,j,k,ix,n2,n
    n = max(2*nrx,nlx)
    n2 = n/2
    if(allocated(NBlock)) deallocate(NBlock)
    allocate(NBlock(0:n2,0:n2,0:n2,0:n2,-n:n,-n:n))
    NBlock = 0
    xni = 0
    do ii = 0,n2
       sni = xni
       xnj = 0
       do jj = 0,n2
          snj = xnj
          xnk = 0
          do kk = 0,n2
             snk = xnk
             xnl = 0
             do ll = 0,n2
                snl = xnl
                xli = 0
                do i = -N+2*ii,N-2*ii
                   sli = xli
                   xlj = 0
                   do j = -N+2*jj,N-2*jj
                      slj = xlj
                      Nblock(ii,jj,kk,ll,i,j)=sni+snj+snk+snl+sli+slj
                      do k = max(-n+2*kk,-min(n-2*ll,-j)+i+j), &
                           min(min(n-2*kk,-i),n-2*ll+i+j)
                         xni = xni + 1
                         xnj = xnj + 1
                         xnk = xnk + 1
                         xnl = xnl + 1
                         xli = xli + 1
                         xlj = xlj + 1
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine calculate_Nblock

  !======================================================================
  !>  Calculates what  is called the LBlock. LBlock
  !>  is a number that is used to obtain the
  !>  transformation from the radial component quantum numbers
  !>  (\f$n\f$'s and \f$\Lambda\f$'s) into the index of an array that
  !>  contains only the non-zero Two body potential matrix elements.
  !>
  !>  The Lblock is the part of that
  !>  index that depends on five of the radial quantum  numbers
  !>  (which would require an \f$2 N^5\f$ array to store in memory)
  !>  but can be calculated on the fly rather quickly. The other part of
  !>  that index is called the NBlock and is calculated and stored in
  !>  memory by calculate_nblock().
  !>
  !>  @result \f$
  !>   LB = 1 + \Lambda_k - \max(-N+2n_{rk},-\min(N-2n_{rl},-\Lambda_j)
  !>                               +\Lambda_i+\Lambda_j) \f$,
  !>  where \f$N\f$ is the number of shells.
  !======================================================================
  function l_block(kn,ln,il,jl,kl,n) result(lb)
    implicit none
    integer(ipr), intent(in) :: kn !< third radial quantum number
    integer(ipr), intent(in) :: ln !< fourth radial quantum number
    integer(ipr), intent(in) :: il !< first orbital quantum number
    integer(ipr), intent(in) :: jl !< second orbital quantum number
    integer(ipr), intent(in) :: kl !< third orbital quantum number
    integer(ipr), intent(in) :: n !< total number of shells
    integer(ipr) :: lb
    ! integer :: i,j,k,mi,mi1,mi2,mi3,mi4,mi5,mi6,mi7,mi8,mi9,ii,mn2
    ! integer :: ma,ma1,ma2,ma3,ma4,ma5,ma6,ma7,ma8,ma9,mi10,m2,mo,sx
    ! integer(ipr) :: i,j,k,dn

    lb = 1 + kl - max(-N+2*kn,-min(n-2*ln,-jl)+il+jl)

    !The code comented bellow is from when the Lblock contained also
    !the sumations depending on in and kn. It is no longer needed
    !but is kept in case it becomes necessary later.

    ! mi  = min(N-2*kn,-il)
    ! mi1 = min(jl-1,-n+2*ln)
    ! ma1 = max(-N+2*jn,mi1+1)
    ! mi2 = min(mi1,2*(kn-ln)-il)
    ! ma2 = max(-N+2*jn,mi2+1)
    ! mi3 = min(jl-1,(2*kn-n-il-mod(abs(2*kn-n-il),2))/2)
    ! ma3 = max(ma1,mi3+1)
    ! mi4 = min(mi2,mi-n+2*ln-il)
    ! ma4 = max(-n+2*jn,mi4+1)
    ! mi5 = min(mi1,mi-n+2*ln-il)
    ! ma5 = max(ma2,mi5+1)
    ! mi6 = min(mi3,mi-n+2*ln-il)
    ! ma6 = max(ma1,mi6+1)
    ! mi7 = min(jl-1,mi-n+2*ln-il)
    ! ma7 = max(ma3,mi7+1)
    ! ma8 = max(-N+2*jn,2*(kn+ln-n)-il)
    ! mi8 = min(mi1,mi+n-2*ln-il)
    ! ma9 = max(ma1,2*(kn+ln-n)-il)
    ! mi9 = min(mi7,n-2*ln)
    ! mi10= min(jl-1,(mi-il-mod(abs(mi-il),2))/2)
    ! if(ma8.le.mi4) then
    !    lb = lb -((-1+ma8-mi4)*(2+2*il-4*kn-4*ln+ma8+mi4+4*n))/2
    ! endif
    ! if(ma4.le.mi2.and.-N+2*kn.le.mi) then
    !    lb = lb + (1 - ma4 + mi2)*(1 - 2*kn + mi + n)
    ! endif
    ! if(ma2.le.mi5.and.-n+2*ln.le.n-2*ln) then
    !    lb = lb + (1 - ma2 + mi5)*(1 - 4*ln + 2*n)
    ! endif
    ! if(ma5.le.mi8) then
    !    lb = lb + ((-1+ma5-mi8)*(-2+2*il+4*ln+ma5-2*mi+mi8-2*n))/2
    ! endif
    ! if(ma9.le.mi6) then
    !    lb = lb -((-1+ma9-mi6)*(2+2*il-4*kn-4*ln+ma9+mi6+4*n))/2
    ! endif
    ! if(ma6.le.mi3.and.-N+2*kn.le.mi) then
    !    lb = lb + (1 - ma6 + mi3)*(1 - 2*kn + mi + n)
    ! endif
    ! if(ma3.le.mi9) then
    !    lb = lb + ((-1+ma3-mi9)*(-2+4*ln+ma3+mi9-2*n))/2
    ! endif
    ! if(ma7.le.mi10) then
    !    lb = lb + (-1 + ma7 - mi10)*(-1 + il + ma7 - mi + mi10)
    ! endif
    ! do i = -n+2*in,il-1
    !    do j = -n+2*jn,n-2*jn
    !       do k = max(-n+2*kn,-min(n-2*ln,-j)+i+j),min(min(n-2*kn,-i),n-2*ln+i+j)
    !          lb = lb + 1
    !       enddo
    !    enddo
    ! enddo
    ! ma = max(-N+2*jn,-N+2*ln+1)
    ! mi1 = min(il-1,-n+2*kn)
    ! ma1 = max(-n+2*in,mi1+1)
    ! mi2 = min(mi1,-n+2*(kn+ln-jn))
    ! ma2 = max(-N+2*in,-n+2*kn)
    ! if(ma2.le.mi2) then
    !    lb = lb -((-1+ma2-mi2)*(6+12*kn**2+4*ma2+ma2**2+5*mi2+&
    !         ma2*mi2+mi2**2+3*(3+ma2+mi2)*n+3*n**2-&
    !         6*kn*(3+ma2+mi2+2*n)))/6 !si1
    ! endif
    ! mi2 = min(N-2*jn,(-n+2*kn-mi1-mod(-n+2*kn-mi1,2))/2)
    ! ma2 = mi2+1
    ! mi3 = min(mi2,2*(ln-kn)-mi1)
    ! ma3 = max(ma,mi3+1)
    ! mi4 = min(mi3,-n+2*(kn+ln-in))
    ! ma4 = max(ma,mi4+1)
    ! ma5 = max(ma,-2*(n-kn-ln)-mi1)
    ! if(ma5.le.mi4) then
    !    lb = lb -((-1+ma5-mi4)*(6+12*kn**2+12*ln**2+4*ma5+ma5**2+&
    !         9*mi1+3*ma5*mi1+3*mi1**2+5*mi4+ma5*mi4+3*mi1*mi4+&
    !         mi4**2+6*kn*(-3+4*ln-ma5-2*mi1-mi4-4*n)+&
    !         6*(3+ma5+2*mi1+mi4)*n+12*n**2-&
    !         6*ln*(3+ma5+2*mi1+mi4+4*n)))/6 !si2
    ! endif
    ! if(ma4.le.mi3.and.-n+2*in.le.mi1) then
    !    lb = lb+((-1+ma4-mi3)*(-1+2*in-mi1-n)*&
    !         (2+2*in-4*kn-4*ln+ma4+mi1+mi3+3*n))/2 !si3
    ! endif
    ! mi4 = min(mi2,N+2*(ln-kn-in))
    ! ma4 = max(ma3,mi4+1)
    ! mi5 = min(mi4,-n+2*(kn+ln-in))
    ! ma5 = max(ma3,mi5+1)
    ! if(ma3.le.mi5) then
    !    lb = lb + (1-ma3+mi5)*(-1+4*kn-2*n)*(-1+2*kn-n) !si4
    ! endif
    ! if(ma5.le.mi4) then
    !    lb = lb +((-1+ma5-mi4)*(-6+12*in**2-36*kn**2-6*ln+12*ln**2+&
    !         ma5-6*ln*ma5+ma5**2+2*mi4-6*ln*mi4+ma5*mi4+mi4**2+&
    !         3*(-5-4*ln+ma5+mi4)*n-9*n**2+&
    !         6*in*(1-4*kn-4*ln+ma5+mi4+2*n)+&
    !         6*kn*(5+4*ln-ma5-mi4+6*n)))/6 !si5
    ! endif
    ! if(ma3.le.mi4) then
    !    lb = lb +((-1+ma3-mi4)*(4*kn-4*ln+ma3+2*mi1+mi4)*&
    !         (-1+4*kn-2*n))/2 !si6
    ! endif
    ! if(ma4.le.mi2.and.-N+2*in.le.mi1) then
    !    lb = lb + (1-ma4+mi2)*(1-2*in+mi1+n)*(1-4*kn+2*n)!si7
    ! endif
    ! mi3 = min(N-2*jn,kn-in)
    ! ma3 = max(ma2,mi3+1)
    ! mi4 = min(mi3,-n+4*kn-2*ln)
    ! ma4 = max(ma2,mi4+1)
    ! mi5 = min(mi4,n+2*(ln-kn-in))
    ! ma5 = max(ma2,mi5+1)
    ! mi6 = min(mi5,-n+2*(kn+ln-in))
    ! ma6 = max(ma2,mi6+1)
    ! if(ma2.le.mi6) then
    !    lb = lb + (1-ma2+mi6)*(-1+4*kn-2*n)*(-1+2*kn-n)!si8
    ! endif
    ! if(ma6.le.mi5) then
    !    lb = lb+((-1+ma6-mi5)*(-6+12*in**2-36*kn**2-6*ln+12*ln**2+&
    !         ma6-6*ln*ma6+ma6**2+2*mi5-6*ln*mi5+ma6*mi5+mi5**2+&
    !         3*(-5-4*ln+ma6+mi5)*n-9*n**2+&
    !         6*in*(1-4*kn-4*ln+ma6+mi5+2*n)+&
    !         6*kn*(5+4*ln-ma6-mi5+6*n)))/6!si9
    ! endif
    ! mi6 = min(mi5,-n+4*kn-2*ln-1)
    ! if(ma2.le.mi6) then
    !    lb = lb -((-1+ma2-mi6)*(-1+4*kn-2*n)*&
    !         (-8*kn+4*ln+ma2+mi6+2*n))/2 !si10
    ! endif
    ! if(ma5.le.mi4) then
    !    lb = lb -((-1+ma5-mi4)*(-1+2*in-2*kn+ma5+mi4)*&
    !         (-1+4*kn-2*n))!si11
    ! endif
    ! mi5 = min(mi3,-n+2*(kn+ln-in))
    ! ma5 = max(ma4,mi5+1)
    ! mi6 = min(mi5,n-2*ln)
    ! if(ma4.le.mi6) then
    !    lb = lb -((-1+ma4-mi6)*(6+12*ln**2-5*ma4+ma4**2-4*mi6+&
    !       ma4*mi6+mi6**2+6*ln*(-3+ma4+mi6-2*n)-&
    !       3*(-3+ma4+mi6)*n+3*n**2))/6 !si12
    ! endif
    ! if(ma5.le.mi3) then
    !    lb = lb +(-1+ma5-mi3)*(-1+2*in-2*kn+ma5+mi3)*&
    !         (1+in-kn-2*ln+n)!si13
    ! endif
    ! mi4 = min(mi3,2*(ln-kn)-mi1)
    ! ma4 = max(ma2,mi4+1)
    ! mi5 = min(mi4,n-2*ln)
    ! if(ma2.le.mi5) then
    !    lb = lb +((-1+ma2-mi5)*(4*ma2**2+3*ma2*mi1+4*ma2*mi5+3*mi1*mi5+&
    !      4*mi5**2-6*kn*(-2+4*ln+ma2+mi5-2*n)-&
    !      3*(ma2+2*mi1+mi5)*n-6*n**2+12*ln*(ma2+mi1+mi5+n)-&
    !      2*(4*ma2+3*mi1+2*mi5+3*n)))/6 !si14
    ! endif
    ! mi5 = min(mi3,-n+4*kn-2*ln)
    ! ma5 = max(ma4,mi5+1)
    ! mi6 = min(mi5,(n-2*kn-mi1-mod(n-2*kn-mi1,2))/2)
    ! ma6 = max(ma4,mi6+1)
    ! if(ma4.le.mi6) then
    !    lb = lb +((-1+ma4-mi6)*(-36*kn**2-5*ma4+4*ma4**2-3*mi1+&
    !         6*ma4*mi1+3*mi1**2-mi6+4*ma4*mi6+6*mi1*mi6+4*mi6**2-&
    !         3*n-6*(ma4+mi1+mi6)*n-9*n**2+&
    !         6*kn*(1+2*ma4+2*mi1+2*mi6+6*n)))/6!si15
    ! endif
    ! if(ma6.le.mi5.and.1.le.2*(n-2*kn)) then
    !    lb = lb +(1-ma6+mi5)*(-1+4*kn-2*n)*(2*kn-n)!si16
    ! endif
    ! mi4 = min(mi3,n-2*ln)
    ! if(ma5.le.mi4) then
    !    lb = lb +((-1+ma5-mi4)*(12*ln**2+(-2+ma5)*ma5-mi4+ma5*mi4+&
    !         mi4**2+6*ln*(-1+ma5+mi4)-&
    !         6*kn*(-2+4*ln+ma5+mi4-2*n)-3*n*(1+n)))/3 !si17
    ! endif
    ! mi6 = min(mi3,(n-2*kn-mi1-mod(n-2*kn-mi1,2))/2)
    ! ma6 = max(ma5,mi6+1)
    ! if(ma5.le.mi6) then
    !    lb = lb +((-1+ma5-mi6)*(4*kn**2-4*ln**2-ma5+ma5**2-mi1+&
    !         2*ma5*mi1+mi1**2+ma5*mi6+2*mi1*mi6+mi6**2+&
    !         kn*(-2+4*ma5+4*mi1+4*mi6-4*n)-&
    !         2*ln*(-1+ma5+mi6-2*n)-(ma5+2*mi1+mi6)*n))/2!si18
    ! endif
    ! mi4 = min(mi3,n-2*ln-1)
    ! if(ma6.le.mi4) then
    !    lb = lb -((-1+ma6-mi4)*(12*ln**2-2*ma6+ma6**2-mi4+ma6*mi4+&
    !         mi4**2+6*ln*(-1+ma6+mi4-2*n)-3*(-1+ma6+mi4)*n+&
    !         3*n**2))/6!si19
    ! endif
    ! mi4 = min(n-2*jn,2*(ln-kn)-mi1)
    ! ma4 = max(ma3,mi4+1)
    ! mi5 = min(mi4,n-2*ln)
    ! if(ma3.le.mi5.and.-N+2*in.le.mi1) then
    !    lb = lb +((-1+ma3-mi5)*(-2+4*ln+ma3+mi5-2*n)*&
    !         (1-2*in+mi1+n))/2!si20
    ! endif
    ! mi5 = min(N-2*jn,N+2*(ln-kn-in))
    ! ma5 = max(ma4,mi5+1)
    ! mi6 = min(mi5,n-2*ln)
    ! if(ma4.le.mi6) then
    !    lb = lb -((-1+ma4-mi6)*(6-24*ln**2-7*ma4+2*ma4**2-5*mi6+&
    !         2*ma4*mi6+2*mi6**2+6*in*(-2+4*ln+ma4+mi6-2*n)+&
    !         6*kn*(-2+4*ln+ma4+mi6-2*n)-6*(-2+ma4+mi6)*n+&
    !         6*n**2))/6 !si21
    ! endif
    ! mi6 = min(mi5,(n-2*kn-mi1-mod(n-2*kn-mi1,2))/2)
    ! ma6 = max(ma4,mi6+1)
    ! if(ma4.le.mi6) then
    !    lb = lb +((-1+ma4-mi6)*(4*kn**2-4*ln**2-ma4+ma4**2-mi1+&
    !         2*ma4*mi1+mi1**2+ma4*mi6+2*mi1*mi6+mi6**2+&
    !         kn*(-2+4*ma4+4*mi1+4*mi6-4*n)-&
    !         2*ln*(-1+ma4+mi6-2*n)-(ma4+2*mi1+mi6)*n))/2!si22
    ! endif
    ! mi6 = min(mi5,n-2*ln-1)
    ! if(ma6.le.mi6) then
    !    lb = lb -((-1+ma6-mi6)*(12*ln**2-2*ma6+ma6**2-mi6+ma6*mi6+&
    !         mi6**2+6*ln*(-1+ma6+mi6-2*n)-3*(-1+ma6+mi6)*n+&
    !         3*n**2))/6!si23
    ! endif
    ! mi6 = min(n-2*jn,(n-2*kn-mi1-mod(n-2*kn-mi1,2))/2)
    ! ma6 = max(ma5,mi6+1)
    ! if(ma5.le.mi6.and.-N+2*in.le.mi1) then
    !    lb = lb -((-1+ma5-mi6)*(-2+2*in+4*kn+2*ma5+mi1+2*mi6-3*n)*&
    !     (-1+2*in-mi1-n))/2!si24
    ! endif
    ! mi6 = N+min(-2*jn,-kn-in)
    ! if(ma6.le.mi6) then
    !    lb = lb -((-1+ma6-mi6)*(6-18*kn-11*ma6-7*mi6+18*n+&
    !       2*(6*in**2+3*in*(-3+4*kn+2*ma6+2*mi6-4*n)+&
    !          2*(3*kn**2+ma6**2+ma6*mi6+mi6**2+&
    !            3*kn*(ma6+mi6-2*n)-3*(ma6+mi6)*n+3*n**2))))/6!si25
    ! endif
    ! mi2 = min(il-1,0)
    ! mi3 = min(mi2,n-4*ln+2*kn)
    ! ma3 = max(ma1,mi3+1)
    ! mi4 = min(mi3,-n+2*(kn+ln-jn))
    ! ma4 = max(ma1,mi4+1)
    ! if(ma1.le.mi4) then
    !    lb = lb -((-1+ma1-mi4)*(6+12*kn**2+4*ma1+ma1**2+5*mi4+&
    !         ma1*mi4+mi4**2+3*(3+ma1+mi4)*n+3*n**2-&
    !         6*kn*(3+ma1+mi4+2*n)))/6!si26
    ! endif
    ! if(ma4.le.mi3.and.jn.le.ln) then
    !    lb = lb +((-1+2*jn-2*ln)*(-1+ma4-mi3)* &
    !         (2+2*jn-4*kn-2*ln+ma4+mi3+2*n))/2!si27
    ! endif
    ! mi4 = min(mi2,n+2*(kn-ln-jn))
    ! ma4 = max(ma3,mi4+1)
    ! mi5 = min(mi4,-n+2*(kn+ln-jn))
    ! ma5 = max(ma3,mi5+1)
    ! if(ma3.le.mi5) then
    !    lb = lb +(1-ma3+mi5)*(-1+4*ln-2*n)*(-1+2*ln-n)!si28
    ! endif
    ! if(ma5.le.mi4) then
    !    lb = lb +((-1+ma5-mi4)*(-6+12*jn**2+12*kn**2+30*ln-36*ln**2+&
    !         ma5-6*ln*ma5+ma5**2+2*mi4-6*ln*mi4+ma5*mi4+mi4**2+&
    !         6*kn*(-1+4*ln-ma5-mi4-2*n)+&
    !         3*(-5+12*ln+ma5+mi4)*n-9*n**2+&
    !         6*jn*(1-4*kn-4*ln+ma5+mi4+2*n)))/6!si29
    ! endif
    ! if(ma3.le.mi4) then
    !    lb = lb +((-1+ma3-mi4)*(-1+4*ln-2*n)*&
    !         (-4*kn+8*ln+ma3+mi4-2*n))/2!si30
    ! endif
    ! if(ma4.le.mi2.and.jn.le.ln) then
    !    lb = lb +(1-2*jn+2*ln)*(1-ma4+mi2)*(1-4*ln+2*n)!si31
    ! endif
    ! mi3 = min(il-1,ln-jn)
    ! ma3 = max(mi2,mi3)+1
    ! mi4 = min(mi3,-n+4*ln-2*kn)
    ! ma4 = max(mi2,mi4)+1
    ! mi5 = min(mi4,n+2*(kn-ln-jn))
    ! ma5 = max(mi2,mi5)+1
    ! mi6 = min(mi5,-n+2*(kn+ln-jn))
    ! ma6 = max(mi2,mi6)+1
    ! if(mi2+1.le.mi6) then
    !    lb = lb +(-mi2+mi6)*(-1+4*ln-2*n)*(-1+2*ln-n)!si32
    ! endif
    ! if(ma6.le.mi5) then
    !    lb = lb +((-1+ma6-mi5)*(-6+12*jn**2+12*kn**2+30*ln-36*ln**2+&
    !         ma6-6*ln*ma6+ma6**2+2*mi5-6*ln*mi5+ma6*mi5+mi5**2+&
    !         6*kn*(-1+4*ln-ma6-mi5-2*n)+&
    !         3*(-5+12*ln+ma6+mi5)*n-9*n**2+&
    !         6*jn*(1-4*kn-4*ln+ma6+mi5+2*n)))/6!si33
    ! endif
    ! mi6 = min(mi5,-n+4*ln-2*kn-1)
    ! if(mi2+1.le.mi6) then
    !    lb = lb -((mi2-mi6)*(-1+4*ln-2*n)*&
    !         (1+4*kn-8*ln+mi2+mi6+2*n))/2!si34
    ! endif
    ! if(ma5.le.mi4) then
    !    lb = lb -((-1+ma5-mi4)*(-1+2*jn-2*ln+ma5+mi4)*&
    !         (-1+4*ln-2*n))!si35
    ! endif
    ! mi5 = min(mi3,-n+2*(kn+ln-jn))
    ! ma5 = max(ma4,mi5+1)
    ! if(ma4.le.mi5) then
    !    lb = lb -((-1+ma4-mi5)*(6+12*kn**2-5*ma4+ma4**2-4*mi5+&
    !         ma4*mi5+mi5**2+6*kn*(-3+ma4+mi5-2*n)-&
    !         3*(-3+ma4+mi5)*n+3*n**2))/6!si36
    ! endif
    ! if(ma5.le.mi3) then
    !    lb = lb +(-1+ma5-mi3)*(-1+2*jn-2*ln+ma5+mi3)*&
    !         (1+jn-2*kn-ln+n)!si37
    ! endif
    ! mi4 = min(mi3,n-4*ln+2*kn)
    ! ma4 = max(mi2,mi4)+1
    ! if(mi2+1.le.mi4) then
    !    lb = lb+((mi2-mi4)*(6*kn*(1+mi2+mi4)+&
    !         2*(-1+mi2**2+mi2*mi4+mi4**2)-3*(1+mi2+mi4)*n))/3!si38
    ! endif
    ! mi5 = min(mi3,-n+4*ln-2*kn)
    ! ma5 = max(ma4,mi5+1)
    ! mi6 = min(mi5,n-2*ln)
    ! ma6 = max(ma4,mi6+1)
    ! if(ma4.le.mi6) then
    !    lb = lb +((-1+ma4-mi6)*(4*ma4**2-mi6+&
    !         ma4*(-5+24*ln+4*mi6-12*n)+4*mi6*(6*ln+mi6-3*n)))/6!si39
    ! endif
    ! if(ma6.le.mi5.and.1.le.2*(n-2*ln)) then
    !    lb = lb + (1 - ma6 + mi5)*(-1 + 4*ln - 2*n)*(2*ln - n)!si40
    ! endif
    ! if(ma5.le.mi3) then
    !    lb = lb +((-1+ma5-mi3)*(12*kn**2+(-2+ma5)*ma5-mi3+ma5*mi3+&
    !      mi3**2+6*kn*(-1-4*ln+ma5+mi3)-&
    !      6*ln*(-2+ma5+mi3-2*n)-3*n*(1+n)))/3!si41
    ! endif
    ! mi6 = min(mi3,n-2*ln)
    ! ma6 = max(ma5,mi6+1)
    ! if(ma5.le.mi6) then
    !    lb = lb +((-1+ma5-mi6)*(-4*kn**2+4*ln*(-1+4*ln)-ma5+ma5**2+&
    !         ma5*mi6+mi6**2-2*kn*(-1+ma5+mi6-2*n)+&
    !         8*ln*(ma5+mi6-2*n)+n-3*(ma5+mi6)*n+3*n**2))/2!si42
    ! endif
    ! if(ma6.le.mi3) then
    !    lb = lb -((-1+ma6-mi3)*(12*kn**2-2*ma6+ma6**2-mi3+ma6*mi3+&
    !         mi3**2+6*kn*(-1+ma6+mi3-2*n)-3*(-1+ma6+mi3)*n+&
    !         3*n**2))/6!si43
    ! endif
    ! mi4 = min(il-1,n-4*ln+2*kn)
    ! ma4 = max(ma3,mi4+1)
    ! if(ma3.le.mi4.and.jn.le.ln) then
    !    lb = lb -((-1+2*jn-2*ln)*(-1+ma3-mi4)*&
    !         (-2+4*kn+ma3+mi4-2*n))/2!si44
    ! endif
    ! mi5 = min(il-1,n+2*(kn-ln-jn))
    ! ma5 = max(ma4,mi5+1)
    ! if(ma4.le.mi5) then
    !    lb = lb -((-1+ma4-mi5)*(6-24*kn**2+24*kn*ln-7*ma4-5*mi5+&
    !         6*jn*(-2+4*kn+ma4+mi5-2*n)+12*n+&
    !         2*(ma4**2+ma4*mi5+mi5**2+3*ln*(-2+ma4+mi5-2*n)-&
    !         3*(ma4+mi5)*n+3*n**2)))/6!si45
    ! endif
    ! mi6 = min(mi5,n-2*ln)
    ! ma6 = max(ma4,mi6+1)
    ! if(ma4.le.mi6) then
    !    lb = lb+((-1+ma4-mi6)*(-4*kn**2+4*ln*(-1+4*ln)-ma4+ma4**2+&
    !         ma4*mi6+mi6**2-2*kn*(-1+ma4+mi6-2*n)+&
    !         8*ln*(ma4+mi6-2*n)+n-3*(ma4+mi6)*n+3*n**2))/2!si46
    ! endif
    ! if(ma6.le.mi5) then
    !    lb = lb -((-1+ma6-mi5)*(12*kn**2-2*ma6+ma6**2-mi5+ma6*mi5+&
    !         mi5**2+6*kn*(-1+ma6+mi5-2*n)-3*(-1+ma6+mi5)*n+&
    !         3*n**2))/6!si47
    ! endif
    ! mi6 = min(il-1,n-2*ln)
    ! ma6 = max(ma5,mi6+1)
    ! if(ma5.le.mi6.and.jn.le.ln) then
    !    lb = lb -((-1+2*jn-2*ln)*(-1+ma5-mi6)*&
    !         (-1+jn+3*ln+ma5+mi6-2*n))!si48
    ! endif
    ! if(ma6.le.il-1) then
    !    lb = lb+((il-ma6)*(17+4*il**2-30*jn-30*ln-15*ma6+&
    !         il*(-15+12*jn+12*ln+4*ma6-12*n)+30*n+&
    !         4*(3*jn**2+6*jn*ln+3*ln**2+3*jn*ma6+3*ln*ma6+ma6**2-&
    !         3*(2*(jn+ln)+ma6)*n+3*n**2)))/6!si49
    ! endif
    ! mi2 = min(il-1,-n+ln+jn)
    ! ma2 = max(ma1,mi2+1)
    ! mi3 = min(N-2*jn,(-n-mi2-mod(abs(-n-mi2),2))/2+kn)
    ! ma3 = max(ma,mi3+1)
    ! if(ma.le.mi3.and.ma1.le.mi2) then
    !    lb = lb +((-1+ma1-mi2)*(-1+ma-mi3)*&
    !         (2-4*kn-4*ln+ma+ma1+mi2+mi3+4*n))/2!si50
    ! endif
    ! mi4 = min(N-2*jn,(-n-ma1-mod(abs(-n-ma1),2))/2+kn)
    ! ma4 = max(ma3,mi4+1)
    ! if(ma3.le.mi4) then
    !    lb = lb -((-1+ma3-mi4)*(-2+2*kn+4*ln-ma1-3*n)*&
    !         (-1-2*kn+ma1+ma3+mi4+n))/2!si51
    ! endif
    ! if(ma3.le.mi4) then
    !    lb = lb +((-1+ma3-mi4)*(4*ma3**2+3*ma3*mi2+4*ma3*mi4+3*mi2*mi4+&
    !         4*mi4**2-6*kn*(-2+4*ln+ma3+mi4-2*n)-&
    !         3*(ma3+2*mi2+mi4)*n-6*n**2+12*ln*(ma3+mi2+mi4+n)-&
    !         2*(4*ma3+3*mi2+2*mi4+3*n)))/6!si52
    ! endif
    ! mi5 = min(N-2*jn,n-2*ln)
    ! if(ma4.le.mi5.and.ma1.le.mi2) then
    !    lb = lb -((-1+ma1-mi2)*(-1+ma4-mi5)*&
    !         (-2+4*ln+ma4+mi5-2*n))/2!si53
    ! endif
    ! mi3 = min(il-1,(-n+2*ln-ma-mod(abs(-n+2*ln-ma),2))/2)
    ! ma3 = max(ma2,mi3+1)
    ! if(ma2.eq.mi3) then
    !    mi4 = min(-n+2*ln-2*ma2,(-n-ma2-mod(abs(-n-ma2),2))/2+kn)
    !    ma4 = max(ma,mi4+1)
    !    if(ma.le.mi4) then
    !       lb = lb -((-1+ma-mi4)*(2-4*kn-4*ln+ma+2*ma2+mi4+4*n))/2!si54
    !    endif
    !    mi5 = min(-n+2*ln-2*ma2,n-2*ln)
    !    if(ma4.le.mi5) then
    !       lb = lb +((-1+ma4-mi5)*(-2+4*ln+ma4+mi5-2*n))/2!si55
    !    endif
    !    mi4 = min(N-2*jn,(-n-ma2-mod(abs(-n-ma2),2))/2+kn)
    !    ma4 = max(-n+2*ln-2*ma2,mi4)+1
    !    if(-n+2*ln-2*ma2+1.le.mi4) then
    !       lb = lb +(1-2*kn-ma2+n)*(-2*ln+2*ma2+mi4+n)!si56
    !    endif
    !    mi5 = min(N-2*jn,-ma2)
    !    if(ma4.le.mi5) then
    !       lb = lb + (-1 + ma4 - mi5)*(-1 + 2*ma2 + ma4 + mi5)!si57
    !    endif
    ! elseif(ma2.lt.mi3) then
    !    if(mod(-mi3,2).eq.1) then
    !       mi4 = min(-n+2*ln-2*mi3,(-n-mi3-mod(abs(-n-mi3),2))/2+kn)
    !       ma4 = max(ma,mi4+1)
    !       if(ma.le.mi4) then
    !          lb = lb -((-1+ma-mi4)*&
    !               (2-4*kn-4*ln+ma+2*mi3+mi4+4*n))/2!si58
    !       endif
    !       mi5 = min(-n+2*ln-2*mi3,n-2*ln)
    !       if(ma4.le.mi5) then
    !          lb = lb +((-1+ma4-mi5)*(-2+4*ln+ma4+mi5-2*n))/2!si59
    !       endif
    !       mi4 = min(n-2*jn,(-n-mi3-mod(abs(-n-mi3),2))/2+kn)
    !       ma4 = max(-n+2*ln-2*mi3+1,mi4+1)
    !       if(-n+2*ln-2*mi3+1.le.mi4) then
    !          lb = lb +(1-2*kn-mi3+n)*(-2*ln+2*mi3+mi4+n)!si60
    !       endif
    !       mi5 = min(N-2*jn,-mi3)
    !       if(ma4.le.mi5) then
    !          lb = lb +(-1 + ma4 - mi5)*(-1 + ma4 + 2*mi3 + mi5)!si61
    !       endif
    !    endif
    !    mi4 = (-ma2-mod(-ma2+1,2))/2
    !    ma4 = (-mi3+mod(-mi3,2))/2
    !    mn2 = (-n-mod(N,2))/2
    !    mi5 = min(mi4,(N-2*ln+kn+mn2-mod(abs(N-2*ln+kn+mn2),3))/3)
    !    ma5 = max(ma4,mi5+1)
    !    if(ma4.le.mi5) then
    !       lb = lb +((-1 + ma4 - mi5)*(-2 + 4*kn + 2*ln - ma - 3*n)*&
    !            (1 + 2*ln - ma + 2*ma4 + 2*mi5 - n))/2 !sx1
    !    endif
    !    mi6 = min(mi4,ma-kn-mn2-1)
    !    ma6 = max(ma5,mi6+1)
    !    mi7 = min(mi6,(n-mod(n,2))/2-ln)
    !    ma7 = max(ma5,mi7+1)
    !    if(ma5.le.mi7) then
    !       lb = lb+((-1+ma5-mi7)*(36*ln**2-3*(-3+ma)*ma+&
    !            2*(-3+ma5*(-7+8*ma5)+mi7+8*ma5*mi7+8*mi7**2)+&
    !            6*ln*(1-2*ma+8*ma5+8*mi7-6*n)-3*n+&
    !            6*(ma - 4*(ma5+mi7))*n+9*n**2))/6 !sx2
    !    endif
    !    if(ma7.le.mi6.and.ma.le.n-2*ln) then
    !       lb = lb+((1-ma7+mi6)*(-2+2*ln+ma-n)*(-1+2*ln+ma-n))/2! sx3
    !    endif
    !    if(ma6.le.mi4) then
    !       lb = lb+((-1+ma6-mi4)*(-2+3*kn**2+ma+ma**2-2*ma*ma6+ma6**2+&
    !            mi4-2*ma*mi4+ma6*mi4+mi4**2-3*mn2+ma6*mn2+mi4*mn2-&
    !            mn2**2+2*ln*(2-2*ma+ma6+mi4+2*mn2)+&
    !            kn*(1+4*ln-4*ma+3*ma6+3*mi4+2*mn2-4*n)+&
    !            2*(-2+2*ma-ma6-mi4-2*mn2)*n))/2!sx4
    !    endif
    !    mi7 = min(mi4,(n-mod(n,2))/2-ln)
    !    ma7 = max(ma6,mi7+1)
    !    if(ma6.le.mi7) then
    !       lb = lb +((-1+ma6-mi7)*(-kn**2+12*ln**2-4*ma6+5*ma6**2+mi7+&
    !            5*ma6*mi7+5*mi7**2+mn2-ma6*mn2-mi7*mn2-mn2**2+&
    !            2*ln*(-1+7*ma6+7*mi7-2*mn2-6*n)-&
    !            kn*(-1+4*ln+ma6+mi7+2*mn2-2*n)+n-7*ma6*n-&
    !            7*mi7*n+2*mn2*n+3*n**2))/2 !sx5
    !    endif
    !    if(ma7.le.mi4) then
    !       lb = lb -((-1+ma7-mi4)*(3*kn**2+12*ln**2-2*ma7+ma7**2-mi4+&
    !            ma7*mi4+mi4**2-3*mn2+3*ma7*mn2+3*mi4*mn2+3*mn2**2+&
    !            6*ln*(-1+ma7+mi4+2*mn2-2*n)+&
    !            3*kn*(-1+4*ln+ma7+mi4+2*mn2-2*n)-&
    !            3*(-1+ma7+mi4+2*mn2)*n+3*n**2))/6!sx6
    !    endif
    !    mi5 = min(mi4,(mn2+kn+n-2*ln-1-mod(abs(mn2+kn+n-2*ln-1),3))/3)
    !    ma5 = max(ma4,mi5+1)
    !    if(ma4.le.mi5) then
    !       lb = lb +((-1+ma4-mi5)*(4*kn**2+4*ln+ma4+5*mi5+&
    !            4*(ma4*(ln+ma4)+(ln+ma4)*mi5+mi5**2)-&
    !            2*kn*(1+4*ln+4*ma4+4*mi5-2*mn2)-2*(ma4+mi5)*mn2+&
    !            (2*kn+4*ln+ma4+mi5-2*mn2)*n-2*n**2-2*(mn2+n)))/2!sx7
    !    endif
    !    mi6 = min(mi5,(n-mod(n,2))/2-jn)
    !    ma6 = max(ma4,mi6+1)
    !    if(ma4.le.mi6) then
    !       lb = lb -((-1+ma4-mi6)*(-ma4+mi6+&
    !            2*(3*kn**2+ma4**2+ma4*mi6+mi6**2-&
    !            3*kn*(ma4+mi6-2*mn2)-3*(ma4+mi6)*mn2+3*mn2**2)))/6!sx8
    !    endif
    !    if(ma6.le.mi5) then
    !       lb = lb +((-1+ma6-mi5)*(8*jn**2-ma6+mi5+8*jn*(ma6+mi5-n)+&
    !            2*(-kn**2+ma6**2+ma6*mi5+mi5**2+kn*(ma6+mi5-2*mn2)+&
    !            ma6*mn2+mi5*mn2-mn2**2-2*(ma6+mi5)*n+n**2)))/2!sx9
    !    endif
    !    mi6 = min(mi4,(n-mod(n,2))/2-jn)
    !    ma6 = max(ma5,mi6+1)
    !    mi7 = min(mi6,(n-1-mod(n-1,2))/2-ln)
    !    if(ma5.le.mi7) then
    !       lb = lb -((-1+ma5-mi7)*(-2*ma5+2*mi7+&
    !            4*(3*ln**2+ma5**2+ma5*mi7+mi7**2+3*ln*(ma5+mi7))-&
    !            6*(2*ln+ma5+mi7)*n+3*n**2))/3!sx10
    !    endif
    !    if(ma6.le.mi4) then
    !       lb = lb +4*(jn-ln)*(-1+ma6-mi4)*(jn+ln+ma6+mi4-n)!sx11
    !    endif
    !    mn2 = (-n+1-mod(N+1,2))/2
    !    mi5 = min(mi4,(N-2*ln+kn+mn2-2-mod(abs(N-2*ln+kn+mn2-2),3))/3)
    !    ma5 = max(ma4,mi5+1)
    !    if(ma4.le.mi5) then
    !       lb = lb +((-1 + ma4 - mi5)*(-2 + 4*kn + 2*ln - ma - 3*n)*&
    !            (3 + 2*ln - ma + 2*ma4 + 2*mi5 - n))/2!sy1
    !    endif
    !    mi6 = min(mi4,ma-kn-mn2-1)
    !    ma6 = max(ma5,mi6+1)
    !    mi7 = min(mi6,(n-1-mod(n-1,2))/2-ln)
    !    ma7 = max(ma5,mi7+1)
    !    if(ma5.le.mi7) then
    !       lb = lb +((-1+ma5-mi7)*(36*ln**2+9*ma-3*ma**2+10*ma5+26*mi7+&
    !            16*(ma5**2+ma5*mi7+mi7**2)-27*n+&
    !            6*(ma-4*(ma5+mi7))*n+9*n**2-&
    !            6*ln*(-9+2*ma-8*ma5-8*mi7+6*n)))/6 !sy2
    !    endif
    !    if(ma7.le.mi6.and.ma.le.n-2*ln) then
    !       lb = lb +((1-ma7+mi6)*(-2+2*ln+ma-n)*(-1+2*ln+ma-n))/2!sy3
    !    endif
    !    if(ma6.le.mi4) then
    !       lb = lb +((-1+ma6-mi4)*(3*kn**2-ma+ma**2+&
    !            ma6-2*ma*ma6+ma6**2+&
    !            2*mi4-2*ma*mi4+ma6*mi4+mi4**2-mn2+ma6*mn2+mi4*mn2-&
    !            mn2**2+2*ln*(2-2*ma+ma6+mi4+2*mn2)+&
    !            kn*(3+4*ln-4*ma+3*ma6+3*mi4+2*mn2-4*n)+&
    !            2*(-2+2*ma-ma6-mi4-2*mn2)*n))/2 !sy4
    !    endif
    !    mi7 = min(mi4,(n-1-mod(n-1,2))/2-ln)
    !    ma7 = max(ma6,mi7+1)
    !    if(ma6.le.mi7) then
    !       lb = lb +((-1+ma6-mi7)*(2-kn**2+12*ln**2+&
    !            4*ma6+5*ma6**2+9*mi7+&
    !            5*ma6*mi7+5*mi7**2+mn2-ma6*mn2-mi7*mn2-mn2**2+&
    !            2*ln*(7+7*ma6+7*mi7-2*mn2-6*n)-&
    !            kn*(-1+4*ln+ma6+mi7+2*mn2-2*n)-&
    !            (7*(1+ma6+mi7)-2*mn2)*n+3*n**2))/2 !sy5
    !    endif
    !    if(ma7.le.mi4) then
    !       lb = lb -((-1+ma7-mi4)*(3*kn**2+12*ln**2-2*ma7+ma7**2-mi4+&
    !            ma7*mi4+mi4**2-3*mn2+3*ma7*mn2+3*mi4*mn2+3*mn2**2+&
    !            6*ln*(-1+ma7+mi4+2*mn2-2*n)+&
    !            3*kn*(-1+4*ln+ma7+mi4+2*mn2-2*n)-&
    !            3*(-1+ma7+mi4+2*mn2)*n+3*n**2))/6!sy6
    !    endif
    !    mi5 = min(mi4,(mn2+kn+n-2*ln-mod(abs(mn2+kn+n-2*ln),3))/3-1)
    !    ma5 = max(ma4,mi5+1)
    !    if(ma4.le.mi5) then
    !       lb = lb +((-1+ma4-mi5)*(4*kn**2+&
    !            4*(2+2*ma4+ma4**2+(3+ma4)*mi5+mi5**2+&
    !            ln*(2+ma4+mi5))-2*(2+ma4+mi5)*mn2+&
    !            (4*ln+ma4+mi5-2*mn2)*n-2*n**2+&
    !            2*kn*(-6-4*ln-4*ma4-4*mi5+2*mn2+n)))/2!sy7
    !    endif
    !    mi6 = min(mi5,(n-1-mod(n-1,2))/2-jn)
    !    ma6 = max(ma4,mi6+1)
    !    if(ma4.le.mi6) then
    !       lb = lb -((-1+ma4-mi6)*(6+6*kn**2+5*ma4+2*ma4**2+7*mi6+&
    !            2*ma4*mi6+2*mi6**2-6*kn*(2+ma4+mi6-2*mn2)-&
    !            6*(2+ma4+mi6)*mn2+6*mn2**2))/6!sy8
    !    endif
    !    if(ma6.le.mi5) then
    !       lb = lb+((-1+ma6-mi5)*(8*jn**2+ma6+3*mi5+4*mn2+&
    !            8*jn*(1+ma6+mi5-n)-4*n+&
    !            2*(-kn**2+ma6**2+ma6*mi5+mi5**2+&
    !            kn*(2+ma6+mi5-2*mn2)+ma6*mn2+mi5*mn2-mn2**2-&
    !            2*(ma6+mi5)*n+n**2)))/2!sy9
    !    endif
    !    mi6 = min(mi4,(n-1-mod(n-1,2))/2-jn)
    !    ma6 = max(ma5,mi6+1)
    !    mi7 = min(mi6,(n-mod(n,2))/2-ln-1)
    !    if(ma5.le.mi7) then
    !       lb = lb -((-1+ma5-mi7)*(3+4*ma5+8*mi7+&
    !            4*(3*ln**2+ma5**2+ma5*mi7+mi7**2+&
    !            3*ln*(1+ma5+mi7))-6*n-6*(2*ln+ma5+mi7)*n+3*n**2))/3!sy10
    !    endif
    !    if(ma6.le.mi4) then
    !       lb = lb +4*(jn-ln)*(-1+ma6-mi4)*(1+jn+ln+ma6+mi4-n)!sy11
    !    endif
    !    if(mod(-ma2,2).eq.0) then
    !       mi4 = min(-n+2*ln-2*ma2,(-n-ma2-mod(abs(-n-ma2),2))/2+kn)
    !       ma4 = max(ma,mi4+1)
    !       if(ma.le.mi4) then
    !          lb = lb -((-1+ma-mi4)*&
    !               (2-4*kn-4*ln+ma+2*ma2+mi4+4*n))/2!si62
    !       endif
    !       mi5 = min(-n+2*ln-2*ma2,n-2*ln)
    !       if(ma4.le.mi5) then
    !          lb = lb +((-1+ma4-mi5)*(-2+4*ln+ma4+mi5-2*n))/2!si63
    !       endif
    !       mi4 = min(n-2*jn,(-n-ma2-mod(abs(-n-ma2),2))/2+kn)
    !       ma4 = max(-n+2*ln-2*ma2,mi4)+1
    !       if(-n+2*ln-2*ma2+1.le.mi4) then
    !          lb = lb +(1-2*kn-ma2+n)*(-2*ln+2*ma2+mi4+n)!si64
    !       endif
    !       mi5 = min(N-2*jn,-ma2)
    !       if(ma4.le.mi5) then
    !          lb = lb +(-1+ma4-mi5)*(-1+2*ma2+ma4+mi5)!si65
    !       endif
    !    endif
    ! endif
    ! mi4 = min(N-2*jn,(-n-il+1-mod(abs(-n-il+1),2))/2+kn)
    ! ma4 = max(ma,mi4+1)
    ! if(ma.le.mi4.and.ma3.le.il-1) then
    !    lb = lb +((il-ma3)*(-1+ma-mi4)*(-3+il+4*kn+ma3-2*n))/2!si66
    ! endif
    ! mi5 = min(N-2*jn,(-n-ma3-mod(abs(-n-ma3),2))/2+kn)
    ! ma5 = max(ma4,mi5+1)
    ! if(ma4.le.mi5) then
    !    lb = lb +((-1+ma4-mi5)*(36*kn**2-3*(-3+ma3)*ma3+ma4+4*ma4**2+&
    !         4*ma4*mi5+mi5*(5+4*mi5)+6*(ma3+2*(ma4+mi5))*n+&
    !         9*n**2-3*(2+n)-6*kn*(-1+2*ma3+4*ma4+4*mi5+6*n)))/6!si67
    ! endif
    ! mi6 = min(mi5,-il+1)
    ! ma6 = max(ma4,mi6+1)
    ! if(ma4.le.mi6) then
    !    lb = lb +(-1+ma4-mi6)*(-1+il+2*kn-n)*(-1+il-2*kn+ma4+mi6+n)!si68
    ! endif
    ! if(ma6.le.mi5) then
    !    lb = lb -((-1+ma6-mi5)*(-ma6+mi5+&
    !         2*(12*kn**2+ma6**2+ma6*mi5+mi5**2+3*(ma6+mi5)*n+&
    !         3*n**2-6*kn*(ma6+mi5+2*n))))/6!si69
    ! endif
    ! mi6 = min(N-2*jn,-il+1)
    ! ma6 = max(ma5,mi6+1)
    ! if(ma5.le.mi6.and.ma3.le.il-1) then
    !    lb = lb +(il-ma3)*(-1+ma5-mi6)*(-2+il+ma3+ma5+mi6)!si70
    ! endif
    ! mi6 = min(N-2*jn,-ma3)
    ! if(ma6.le.mi6) then
    !    lb = lb -((-1+ma6-mi6)*(6-7*ma6-5*mi6+&
    !         2*(3*ma3**2+ma6**2+ma6*mi6+mi6**2+&
    !         3*ma3*(-2+ma6+mi6))))/6!si71
    ! endif
  end function l_block

  !======================================================================
  !>  Given four axial quantum numbers, returns the corresponding index
  !>  in the one-dimensional array that contains all non-zero and
  !>  different axial components.
  !>
  !>  If the sum of the four quantum numbers is not an even number, in
  !>  which case \f$ V^z_{ijkl} = 0\f$, the function returns 0 as the
  !>  index
  !======================================================================
  function zindex(nzi,nzj,nzk,nzl) result(iz)
    implicit none
    integer(ipr) :: nzi !< First  axial quantum number
    integer(ipr) :: nzj !< Second axial quantum number
    integer(ipr) :: nzk !< Third  axial quantum number
    integer(ipr) :: nzl !< Fourth axial quantum number
    integer(ipr) :: iz
    if(mod(nzi+nzj+nzk+nzl,2).eq.0) then
       iz = ZBlock(min(nzi,nzk),min(nzj,nzl),max(nzi,nzk),max(nzj,nzl))
    else
       iz = 0
    endif
  end function zindex

  !======================================================================
  !>  Given four radial and four angular quantum numbers, returns the
  !>  corresponding index in the one-dimensional array that contains all
  !>  non-zero and different radial components of two body potential
  !>  potential matrix elements for a gaussian potential
  !>
  !>  If the quantum numbers correspond to a matrix element such that
  !>  \f$ V^r_{ijkl} = 0\f$, the function returns 0 as the index
  !=====================================================================
  function rindex(nri,nrj,nrk,nrl,li,lj,lk,ll,n) result(ir)
    implicit none
    integer(ipr) :: nri !< First  radial quantum number
    integer(ipr) :: nrj !< Second radial quantum number
    integer(ipr) :: nrk !< Third  radial quantum number
    integer(ipr) :: nrl !< Fourth radial quantum number
    integer(ipr) :: li !< First  angular quantum number
    integer(ipr) :: lj !< Second angular quantum number
    integer(ipr) :: lk !< Third  angular quantum number
    integer(ipr) :: ll !< Fourth angular quantum number
    integer(ipr) :: n !< Total number of shells
    integer(ipr) :: ir
    if(li+lj.eq.ll+lk) then
       if(lk.le.-li.and.ll.le.-lj) then
         ir=nblock(nri,nrj,nrk,nrl,li,lj)+&
              l_block(nrk,nrl,li,lj,lk,n)
       elseif(ll.le.-lj) then
         ir=nblock(nrk,nrj,nri,nrl,-lk,lj)+&
              l_block(nri,nrl,-lk,lj,-li,n)
       elseif(lk.le.-li) then
         ir=nblock(nri,nrl,nrk,nrj,li,-ll)+&
              l_block(nrk,nrj,li,-ll,lk,n)
       else
         ir=nblock(nrk,nrl,nri,nrj,-lk,-ll)+&
              l_block(nri,nrj,-lk,-ll,-li,n)
       endif
    else
       ir = 0
    endif
  end function rindex

  !======================================================================
  !>  Given two equaly sized, square, symmetric and block diagonal
  !>  matrices, returns the trace of the product of such matrices.
  !>
  !>  The matrices are given as two-dimensional arrays. The second index
  !>  indicates the block within the matrix. The first index gives the
  !>  matrix element within each block. In particular the \f$ (i,j) \f$
  !>  pair is transformed into an index given by \f$ i + (j-1)*d \f$,
  !>  where \f$ d \f$ is the size of the block.
  !======================================================================
  function trace_product(A,B) result(tr)
    implicit none
    real(pr),allocatable , dimension(:,:) :: A !< First  Matrix
    real(pr),allocatable , dimension(:,:) :: B !< Second Matrix
    real(pr) :: tr,fac
    integer(ipr) :: ib,ibx,n1,n2,nd,n12,n21
    if(.not.allocated(A).or..not.allocated(B)) then
       tr = 0._pr
       return
    endif
    tr = zero
    do ib = 1,nb
       ibx = ib+nbx
       nd = id(ib)
       do n1 = 1,nd
          do n2 = 1,n1
             if(n1.eq.n2) then
                fac = one
             else
                fac = two
             endif
             n12 = n1+(n2-1)*nd
             tr = tr + fac*A(n12,ib )*B(n12,ib )
             tr = tr + fac*A(n12,ibx)*B(n12,ibx)
          enddo!n2
       enddo!n1
    enddo!ib
  end function trace_product

  !======================================================================
  !>  Given two equaly sized, square, symmetric and block diagonal
  !>  matrices, returns the trace of the product of such matrices.
  !>
  !>  The matrices used with this subroutine are composed of
  !>  two major blocks (one for neutrons and one for protons) which in
  !>  turn are formed by smaller blocks. The subroutine returns a trace
  !>  for each major block and the total trace is the sum of both traces
  !>
  !>  The matrices are given as two dimensional arrays. The second index
  !>  indicates the block within the matrix. The first index gives the
  !>  matrix element within each block. In particular the \f$ (i,j) \f$
  !>  pair is transformed into an index given by \f$ i + (j-1)*d \f$,
  !>  where \f$ d \f$ is the size of the block.
  !======================================================================
  subroutine trace_product_2(A,B,tr1,tr2)
    implicit none
    real(pr),allocatable , dimension(:,:) :: A !< First  Matrix
    real(pr),allocatable , dimension(:,:) :: B !< Second Matrix
    real(pr),intent(out) :: tr1 !< First half of the trace
    real(pr),intent(out) :: tr2 !< Second half of the trace
    real(pr) :: fac
    integer(ipr) :: ib,ibx,n1,n2,nd,n12,n21
    if(.not.allocated(A).or..not.allocated(B)) then
       tr1 = zero
       tr2 = zero
       return
    endif
    tr1 = zero
    tr2 = zero
    do ib = 1,nb
       ibx = ib+nbx
       nd = id(ib)
       do n1 = 1,nd
          do n2 = 1,n1
             if(n1.eq.n2) then
                fac = one
             else
                fac = two
             endif
             n12 = n1+(n2-1)*nd
             tr1 = tr1 + fac*A(n12,ib )*B(n12,ib )
             tr2 = tr2 + fac*A(n12,ibx)*B(n12,ibx)
          enddo!n2
       enddo!n1
    enddo!ib
  end subroutine trace_product_2

  !======================================================================
  !>  Calculates the the trace of a block diagonal matrix.
  !>
  !>  The matrices used with this subroutine are composed of
  !>  two major blocks (one for neutrons and one for protons) which in
  !>  turn are formed by smaller blocks. The subroutine returns a trace
  !>  for each major block and the total trace is the sum of both traces
  !>
  !>  The matrix is given as a two-dimensional array. The second index
  !>  indicates the block within the matrix. The first index gives the
  !>  matrix element within each block. In particular the \f$ (i,j) \f$
  !>  pair is transformed into an index given by
  !>  \f$ i + (j-1)*d \f$, where \f$ d \f$ is the size of the block.
  !======================================================================
  subroutine bdiag_trace(A,tr1,tr2)
    implicit none
    real(pr),allocatable , dimension(:,:), &
         intent(in) :: A !< Block diagonal matrix
    real(pr) , intent(out) :: tr1 !< First  half of the trace
    real(pr) , intent(out) :: tr2 !< Second half of the trace
    integer(ipr) :: ib,ibx,n1,nd,n11
    if(.not.allocated(A)) then
       return
    endif
    tr1 = zero
    tr2 = zero
    do ib = 1,nb
       ibx = ib+nbx
       nd = id(ib)
       do n1 = 1,nd
             n11 = n1+(n1-1)*nd
             tr1 = tr1 + A(n11,ib )
             tr2 = tr2 + A(n11,ibx)
       enddo!n1
    enddo!ib
  end subroutine bdiag_trace

  !======================================================================
  !>  Prints a given block diagonal matrix.
  !>
  !>  The matrix is given as a
  !>  two-dimensional array. The second index indicates the block within
  !>  the matrix. The first index gives the matrix element within each
  !>  block. In particular the \f$ (i,j) \f$ pair is transformed into an
  !>  index given by \f$ i + (j-1)*d \f$, where \f$ d \f$ is the size of
  !>  the block.
  !======================================================================
  subroutine bdiag_print(A)
    implicit none
    real(pr),allocatable , dimension(:,:) :: A !< Block diagonal matrix
    integer(ipr) :: ib,ibx,n1,n2,nd,n12
    if(.not.allocated(A)) then
       return
    endif
    do ib = 1,nb
       ibx = ib+nbx
       nd = id(ib)
       do n1 = 1,nd
          do n2 = 1,nd
             n12 = n2+(n1-1)*nd
             write(*,*) ib ,n1,n2,n12,A(n12,ib )
          enddo
       enddo!n1
    enddo!ib
  end subroutine bdiag_print

  !======================================================================
  !>  Test the orthonormality of the Harmonic oscillator wave functions
  !>  in axial symmetry by printing the integral of all products of two
  !>  wave functions. The integrals are calculated using the appropriate
  !>  gaussian quadrature technique.
  !>
  !>  The two body matrix elements of a contact interaction are also
  !>  printed.
  !>
  !>  This subroutine is used for debugging only and is not used in the
  !>  actual DFT calculations.
  !======================================================================
  subroutine test_HOWF_gauss()
    implicit none
    integer(ipr) :: i,ni,nj,nk,nl,li,lj,lk,ll
    real(pr) :: s,Lni,Lnj,Lnk,Lnl,Nr_i,Nr_j,Nr_k,Nr_l
    real(pr) :: Hni,Hnj,Hnk,Hnl,Nz_i,Nz_j,Nz_k,Nz_l,wher
    real(pr), dimension(1:ngl) :: wlag,xlag
    do ni = 0,n00/2
       do nj = 0,n00/2
          do li = 0,n00-2*ni
             lj = li
             Nr_i = N_radial(ni,li)
             Nr_j = N_radial(nj,lj)
             s = 0._pr
             call GaussLaguerreWX((li+lj)*0.5_pr,wlag,xlag)
             do i = 1,ngl
                call LaguerreL(ni,real(li,kind=pr),xlag(i),Lni)
                call LaguerreL(nj,real(lj,kind=pr),xlag(i),Lnj)
                s = s + wlag(i)*lni*lnj
             enddo
             write(*,'(3i3,f15.8)') ni,nj,li,s*Nr_i*Nr_j*bp**2/2._pr
          enddo
       enddo
    enddo
    write(*,*)
    do ni = 0,n00/2
       do nj = 0,n00/2
          do nk = 0,n00/2
             do nl = 0,n00/2
                do li = -n00+2*ni,n00-2*ni
                   do lj = -n00+2*nj,n00-2*nj
                      do lk = -n00+2*nk,n00-2*nk
                         ll = li+lj-lk
                         if(abs(ll).gt.n00-2*nl) cycle
                         s = MatrixElement_rZR(ni,li,nj,lj,nk,lk,nl,ll)
                         write(*,'(8i3,f15.8)') ni,nj,nk,nl,li,lj,lk&
                              ,ll,s
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    write(*,*)
    do ni = 0,n00
       do nj = 0,n00
          s = 0._pr
          Nz_i = N_axial(ni)
          Nz_j = N_axial(nj)
          do i = 1,ngh
             wher = wh(i)*exp(-xh(i)**2)
             Hni = HermiteH(ni,xh(i))
             Hnj = HermiteH(nj,xh(i))
             s = s + wher*Hni*Hnj*bz
          enddo
          write(*,'(2i3,f15.8)') ni,nj,s*Nz_i*Nz_j
       enddo
    enddo
    write(*,*)
    do ni = 0,n00
       do nj = 0,n00
          do nk = 0,n00
             do nl = 0,n00
                if(mod(ni+nj+nk+nl,2).ne.0) cycle
                write(*,'(4i3,f15.8)') ni,nj,nk,nl,&
                     MatrixElement_zZR(ni,nj,nk,nl)
             enddo
          enddo
       enddo
    enddo
    write(*,*)
  end subroutine test_HOWF_gauss

  !======================================================================
  !>  Given four principal quantum numbers, calculates the two-body
  !>  matrix element of a contact interaction using one-dimensional
  !>  harmonic oscillator wave functions.
  !>
  !>  The matrix element integral is calculated using Gauss-Hermite
  !>  quadrature
  !>
  !>  This function is used for debugging only and is not used in the
  !>  actual DFT calculations.
  !>
  !>  @result \f$ V_{n_i,n_j,n_k,n_l} = \int \phi_{n_i}(z_1;b_z)
  !>          \phi_{n_j}(z_2;b_z) \delta(z_1-z_2) \phi_{n_k}(z_1;b_z)
  !>          \phi_{n_l}(z_2;b_z) dz_1 dz_2 \f$
  !======================================================================
  function MatrixElement_zZR(ni,nj,nk,nl) result(ME)
    implicit none
    integer(ipr), intent(in) :: ni !< first  axial quantum number
    integer(ipr), intent(in) :: nj !< second axial quantum number
    integer(ipr), intent(in) :: nk !< third  axial quantum number
    integer(ipr), intent(in) :: nl !< fourth axial quantum number
    real(pr) :: ME
    integer(ipr) :: i
    real(pr) :: Nz_i,Nz_j,Nz_k,Nz_l,wher,Hni,Hnj,Hnk,Hnl
    ME = 0._pr
    if(mod(ni+nj+nk+nl,2).ne.0) return
    Nz_i = N_axial(ni)
    Nz_j = N_axial(nj)
    Nz_k = N_axial(nk)
    Nz_l = N_axial(nl)
    do i = 1,ngh
       wher = wh(i)*exp(-xh(i)**2)
       Hni = HermiteH(ni,xh(i)/sqrt(2._pr))
       Hnj = HermiteH(nj,xh(i)/sqrt(2._pr))
       Hnk = HermiteH(nk,xh(i)/sqrt(2._pr))
       Hnl = HermiteH(nl,xh(i)/sqrt(2._pr))
       ME = ME + wher*Hni*Hnj*Hnk*Hnl*bz/sqrt(2._pr)
    enddo
    ME = ME*Nz_i*Nz_j*Nz_k*Nz_l
  end function MatrixElement_zZR

  !======================================================================
  !>  Given four radial and four angular quantum numbers, calculates the
  !>  two-body matrix element of a contact interaction using
  !>  two-dimensional harmonic oscillator wave functions in radial
  !>  coordinates
  !>
  !>  The matrix element integral is calculated using Gauss-Laguerre
  !>  quadrature
  !>
  !>  This function is used for debugging only and is not used in the
  !>  actual DFT calculations.
  !>
  !>  @result \f$ V_{n_i,\Lambda_i,n_j,\Lambda_j,n_k,\Lambda_k,n_l,
  !>              \Lambda_l} = \int \phi_{n_i,\Lambda_i}(r_1,\varphi_1;
  !>              b_\perp) \phi_{n_j,\Lambda_j}(r_2,\varphi_2;b_\perp)
  !>              \delta(r_1-r_2) \delta(\varphi_1-\varphi_2)
  !>              \phi_{n_k,\Lambda_k}(r_1,\varphi_1;b_\perp)
  !>              \phi_{n_l,\Lambda_l}(r_2,\varphi_2;b_\perp) dr_1
  !>              d \varphi_2 dr_2 d \varphi_2 \f$
  !======================================================================
  function MatrixElement_rZR(ni,li,nj,lj,nk,lk,nl,ll) result(ME)
    implicit none
    integer(ipr),intent(in) :: ni !<first  radial principal quantum number
    integer(ipr),intent(in) :: nj !<second radial principal quantum number
    integer(ipr),intent(in) :: nk !<third  radial principal quantum number
    integer(ipr),intent(in) :: nl !<fourth radial principal quantum number
    integer(ipr),intent(in) :: li !<first  radial orbital quantum number
    integer(ipr),intent(in) :: lj !<second radial orbital quantum number
    integer(ipr),intent(in) :: lk !<third  radial orbital quantum number
    integer(ipr),intent(in) :: ll !<fourth radial orbital quantum number
    real(pr) :: ME
    integer(ipr) i
    real(pr) :: Nr_i,Nr_j,Nr_k,Nr_l,Lni,Lnj,Lnk,Lnl
    real(pr), dimension(1:ngl) :: wlag,xlag
    ME = 0._pr
    if(li+lj.ne.lk+ll) return
    Nr_i = N_radial(ni,li)
    Nr_j = N_radial(nj,lj)
    Nr_k = N_radial(nk,lk)
    Nr_l = N_radial(nl,ll)
    call GaussLaguerreWX((abs(li)+abs(lj)+abs(lk)+abs(ll))*0.5_pr,&
         wlag,xlag)
    do i = 1,ngl
       call LaguerreL(ni,real(abs(li),kind=pr),xlag(i)/2._pr,Lni)
       call LaguerreL(nj,real(abs(lj),kind=pr),xlag(i)/2._pr,Lnj)
       call LaguerreL(nk,real(abs(lk),kind=pr),xlag(i)/2._pr,Lnk)
       call LaguerreL(nl,real(abs(ll),kind=pr),xlag(i)/2._pr,Lnl)
       ME = ME + wlag(i)*Lni*Lnj*Lnk*Lnl
    enddo
    ME = ME*Nr_i*Nr_j*Nr_k*Nr_l*bp**2*2._pr**(-(abs(li)+abs(lj)+&
         abs(lk)+abs(ll)+4)*0.5_pr)/(2*pi)
  end function MatrixElement_rZR

  !======================================================================
  !>  Given a pair of radial and angular quantum numbers, calculates the
  !>  normalization constant of a  two-dimensional harmonic oscillator
  !>  wave function in radial coordinates
  !>
  !>  This function is used for debugging only and is not used in the
  !>  actual DFT calculations.
  !>
  !>  @result \f$ N_{n_r,\Lambda} = \frac{1}{b_\perp} \left[
  !>              \frac{2n_r!}{(n_r+|\Lambda|)!} \right]^{1/2} \f$
  !======================================================================
  function N_radial(nr,l) result (N_r)
    implicit none
    integer(ipr), intent(in) :: nr !< radial principlar quantum number
    integer(ipr), intent(in) :: l  !< radial orbital quantum number
    real(pr) :: N_r
    real(pr) :: nrf, nrlf
    nrf  = factrl(nr)
    nrlf = factrl(nr+abs(l))
    N_r = sqrt(2*nrf/nrlf)/bp
  end function N_radial

  !======================================================================
  !>  Given a quantum number, calculates the normalization constant of
  !>  a one-dimensional harmonic oscillator wave function
  !>
  !>  This function is used for debugging only and is not used in the
  !>  actual DFT calculations.
  !>
  !>  @result \f$ N_{n_z} = \frac{1}{(b_z\sqrt{\pi}2^{n_z}n_z!)^{1/2}}\f$
  !======================================================================
  function N_axial(nz) result (N_z)
    implicit none
    integer(ipr), intent(in) :: nz !< axial quanutum number
    real(pr) :: N_z
    real(pr) :: nzf
    nzf  = factrl(nz)
    N_z = one/sqrt(bz*sqrt(pi)*(2**nz)*nzf)
  end function N_axial

  !======================================================================
  !>  Calculates the Generalized Laguerre polynomial \f$L_n^\alpha(x)\f$
  !>  along with \f$ \frac{\partial}{\partial x} L_n^\alpha(x) \f$ and
  !>  \f$L_{n-1}^\alpha(x)\f$. The last two are optional
  !======================================================================
  subroutine LaguerreL(n,alpha,x,Ln,Lnp,Lnm1)
    implicit none
    integer(ipr), intent(in) :: n !< Order of the polynomial
    real(pr), intent(in) :: alpha !< &alpha; parameter of the generalized polynomial
    real(pr), intent(in) :: x !< Value where the polynomials are evaluated
    real(pr), intent(out) :: Ln !< Generalized Laguerre Polynomial \f$ L_n^\alpha(x) \f$
    real(pr), optional, intent(out) :: Lnp  !< Derivative of Ln with respect of x
    real(pr), optional, intent(out) :: Lnm1 !< Generalized Laguerre Polynomial \f$ L_{n-1}^\alpha(x) \f$
    real(pr) :: Ljm2, Ljm1, Lj
    integer(ipr) :: j
    Lj = one
    Ljm1 = zero
    do j = 1,n
       Ljm2 = Ljm1
       Ljm1 = Lj
       Lj = ((-x+2*j-1+alpha)*Ljm1-(j-1+alpha)*Ljm2)/real(j,kind=pr)
    end do
    Ln = Lj
    if(present(Lnp)) Lnp = (n*Ln-(n+alpha)*Ljm1)/x
    if(present(Lnm1)) Lnm1 = Ljm1
    return
  end subroutine LaguerreL

  !======================================================================
  !>  Calculates the weights w and nodes x for the Gauss-Laguerre
  !>  quadrature integration with a given parameter \f$ \alpha \f$
  !======================================================================
  subroutine GaussLaguerreWX(alfa,w,x)
    implicit none
    real(pr), intent(in) :: alfa !< \f$\alpha\f$ parameter on the weight function
    real(pr), intent(out), dimension(:) :: x !< Integration nodes
    real(pr), intent(out), dimension(:) :: w !< Integration weights
    integer(ipr), parameter :: maxit = 20
    integer(ipr) :: i,its,n,ai
    real(pr), parameter :: eps = 2.e-14_pr
    real(pr) :: z,z1,Ln,Lnp,Lnm1
    n = size(w)
!    write(*,*) size(w), size(x)
    if(size(w).ne.size(x)) then !Check that both arrays have the same size
       write(*,*) 'Arrays w and x must have the same size in GaussLaguerreWX'
       stop
    end if
    do i = 1,n
       !intial guess for every root of the Laguerre Polynomial
       if(i.eq.1)then
          z=(1+alfa)*(3+0.92_pr*alfa)/(1+2.4_pr*n+1.8_pr*alfa)
       else if(i.eq.2)then
          z=z+(15+6.25_PR*alfa)/(1+0.9_pr*alfa+2.5_pr*n)
       else
          ai=i-2
          z=z+((1+2.55_pr*ai)/(1.9_pr*ai)+1.26_pr*ai*alfa/(1+3.5_pr*ai))&
               *(z-x(i-2))/(1+0.3_pr*alfa)
       endif
       ! Newton method to refine the roots
       do its = 1, maxit
          call LaguerreL(n,alfa,z,Ln,Lnp,Lnm1)
          z1 = z
          z = z1 - Ln/Lnp
!          if(abs(z-z1).le.eps) exit
          if(abs(z-z1).le.eps*z) exit
       enddo
       if(its==maxit+1) write(*,*) 'maxit exceeded in GaussLaguerreWX',its,abs(z-z1)/z,z
       ! Save root and calculate weight
       x(i) = z
       if(alfa.eq.0._pr) then
          w(i) = -1/(n*Lnm1*Lnp)
       else
          w(i) = -gamma(n+alfa)/(factrl(n)*Lnm1*Lnp)
       endif
    end do
  end subroutine GaussLaguerreWX

  !======================================================================
  !> Calculates the Hermite polynomial \f$H_n(x)\f$
  !>
  !> @result \f$ H_n(x) = \left(2x-\frac{d}{dx} \right)^n \cdot 1 \f$
  !======================================================================
  function HermiteH(n,x) result(Hn)
    implicit none
    integer(ipr), intent(in) :: n !< Order of the polynomial
    real(pr), intent(in) :: x !< Value where the polynomial is evaluated
    real(pr) :: Hn  !< Hermite Polynomial \f$ H_n(x) \f$
    real(pr) :: Hjm2, Hjm1, Hj
    integer(ipr) :: j
    Hj = one
    Hjm1 = zero
    do j = 1,n
       Hjm2 = Hjm1
       Hjm1 = Hj
       Hj = two*x*Hjm1-two*Real(j-1,kind=pr)*Hjm2
    end do
    Hn = Hj
    return
  end function HermiteH

end module hfbtho_gogny
