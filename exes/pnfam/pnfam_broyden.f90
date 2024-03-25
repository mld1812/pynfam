!------------------------------------------------------------------------------
!> This modules contains the variables and routines used for applying broyden's
!> method for mixing previous iteration results into the current iteration result.
!> The primary routine "broyden_method" is the same as used in hfbtho v3.00.
!>
!> @authors M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
!> @authors E.M. Ney, UNC Chapel Hill, 2018-
!------------------------------------------------------------------------------
module pnfam_broyden
   use pnfam_logger
   use type_bigblockmatrix, only : bigblockmatrix
   implicit none
   integer, parameter, private :: dp = kind(1d0)
   integer, parameter, private :: pr = kind(1d0)
   integer, parameter, private :: ipr = kind(1)

   character(1) :: qrpa_bbroyden='N'
   integer(ipr) :: broyden_history_size, nbroy
   real(pr)     :: si
   real(pr)     :: qrpa_alphamix = 0.7

   real(dp), allocatable, target :: qrpa_broin(:), qrpa_broout(:)

   integer, allocatable :: sep_points(:) ! Separation points in qrpa_broout for different dRqp
   integer, allocatable :: corr(:) ! From qrpa_broin or qrpa_broout index to dRqp index

contains

   !---------------------------------------------------------------------------
   ! Initialize Broyden's method
   ! dRqp not contributing to dRsp due to zeros in Gre and Gim
   ! resulted from pairing cutoff not be stored in the Broyden vector. 
   !---------------------------------------------------------------------------
   subroutine init_broyden(nxy, Gre, Gim)
      implicit none
      integer, intent(in) :: nxy ! Number of elements in f%mat_n%elem or f%mat_p%elem
      type(bigblockmatrix), intent(in) :: Gre, Gim

      integer :: inc

      if (allocated(qrpa_broin)) deallocate(qrpa_broin,qrpa_broout,sep_points,corr)

      if (allocated(Gre%m11%elem)) then
         allocate(sep_points(8), corr(nxy*8))
      else
         allocate(sep_points(4), corr(nxy*4))
      end if

      ! Initialize corr and sep_points
      inc = 0
      sep_points(:) = 0; corr(:) = 0
      call calc_corr_sep(Gre%m12, Gim%m12, inc, sep_points(1), corr)
      call calc_corr_sep(Gre%m21, Gim%m21, inc, sep_points(2), corr)
      call calc_corr_sep(Gre%m12, Gim%m12, inc, sep_points(3), corr)
      call calc_corr_sep(Gre%m21, Gim%m21, inc, sep_points(4), corr)
      nbroy = sep_points(4)
      if (allocated(Gre%m11%elem)) then
         call calc_corr_sep(Gre%m11, Gim%m11, inc, sep_points(5), corr)
         call calc_corr_sep(Gre%m22, Gim%m22, inc, sep_points(6), corr)
         call calc_corr_sep(Gre%m11, Gim%m11, inc, sep_points(7), corr)
         call calc_corr_sep(Gre%m22, Gim%m22, inc, sep_points(8), corr)
         nbroy = sep_points(8)
      end if

      ! Broyden arrays
      allocate(qrpa_broin(nbroy),qrpa_broout(nbroy))
      qrpa_broin = 0; qrpa_broout=0
      si = 1

      contains

      ! Auxiliary subroutine to establish correspondence
      subroutine calc_corr_sep(Gre, Gim, inc, sep_point, corr)
         use type_bigblockmatrix, only : blockmatrix
         !import, none
         implicit none
         type(blockmatrix), intent(in) :: Gre, Gim
         integer, intent(inout) :: inc, sep_point, corr(:)
         integer :: i
         do i = 1, size(Gre%elem)
            if ((Gre%elem(i) /= 0) .and. (Gim%elem(i) /= 0)) then
               inc = inc + 1
               corr(inc) = i
            end if
         end do
         sep_point = inc
      end subroutine calc_corr_sep

   end subroutine init_broyden

   subroutine set_broinout_from_dRqp(dRqp_re, dRqp_im, io)
      implicit none
      type(bigblockmatrix), intent(in) :: dRqp_re, dRqp_im
      character(1), intent(in) :: io
      ! integer :: nuv, nvec
      real(dp), pointer :: qrpa_bro(:)

      if (.not. allocated(sep_points)) call abort("Error: Broyden not initialized.")

      nullify(qrpa_bro)
      select case(io)
         case('i')
            qrpa_bro => qrpa_broin
         case('o')
            qrpa_bro => qrpa_broout
      end select
      ! nuv   = size(dRqp_re%m12%elem)
      ! nvec  = nbroy/nuv

      ! ! Make sure sizes match up
      ! if (allocated(dRqp_re%m11%elem)) then
      !     if (nvec/=8) call abort
      ! else
      !     if (nvec/=4) call abort
      ! end if

      qrpa_bro(              1:sep_points(1)) = dRqp_re%m12%elem(corr(              1:sep_points(1))) !rex
      qrpa_bro(sep_points(1)+1:sep_points(2)) = dRqp_re%m21%elem(corr(sep_points(1)+1:sep_points(2))) !rey
      qrpa_bro(sep_points(2)+1:sep_points(3)) = dRqp_im%m12%elem(corr(sep_points(2)+1:sep_points(3))) !imx
      qrpa_bro(sep_points(3)+1:sep_points(4)) = dRqp_im%m21%elem(corr(sep_points(3)+1:sep_points(4))) !imy

      if (allocated(dRqp_re%m11%elem)) then
         qrpa_bro(sep_points(4)+1:sep_points(5)) = dRqp_re%m11%elem(corr(sep_points(4)+1:sep_points(5))) !rep
         qrpa_bro(sep_points(5)+1:sep_points(6)) = dRqp_re%m22%elem(corr(sep_points(5)+1:sep_points(6))) !req
         qrpa_bro(sep_points(6)+1:sep_points(7)) = dRqp_im%m11%elem(corr(sep_points(6)+1:sep_points(7))) !imp
         qrpa_bro(sep_points(7)+1:sep_points(8)) = dRqp_im%m22%elem(corr(sep_points(7)+1:sep_points(8))) !imq
      end if

   end subroutine

   subroutine set_dRqp_from_broin(dRqp_re, dRqp_im)
      implicit none
      type(bigblockmatrix), intent(inout) :: dRqp_re, dRqp_im
      ! integer :: nuv, nvec

      if (.not. allocated(sep_points)) call abort("Error: Broyden not initialized.")

      ! nuv   = size(dRqp_re%m12%elem)
      ! nvec  = nbroy/nuv

      ! ! Make sure sizes match up
      ! if (allocated(dRqp_re%m11%elem)) then
      !     if (nvec/=8) call abort
      ! else
      !     if (nvec/=4) call abort
      ! end if

      dRqp_re%m12%elem(corr(              1:sep_points(1))) = qrpa_broin(              1:sep_points(1)) !rex
      dRqp_re%m21%elem(corr(sep_points(1)+1:sep_points(2))) = qrpa_broin(sep_points(1)+1:sep_points(2)) !rey
      dRqp_im%m12%elem(corr(sep_points(2)+1:sep_points(3))) = qrpa_broin(sep_points(2)+1:sep_points(3)) !imx
      dRqp_im%m21%elem(corr(sep_points(3)+1:sep_points(4))) = qrpa_broin(sep_points(3)+1:sep_points(4)) !imy

      if (allocated(dRqp_re%m11%elem)) then
         dRqp_re%m11%elem(corr(sep_points(4)+1:sep_points(5))) = qrpa_broin(sep_points(4)+1:sep_points(5)) !rep
         dRqp_re%m22%elem(corr(sep_points(5)+1:sep_points(6))) = qrpa_broin(sep_points(5)+1:sep_points(6)) !req
         dRqp_im%m11%elem(corr(sep_points(6)+1:sep_points(7))) = qrpa_broin(sep_points(6)+1:sep_points(7)) !imp
         dRqp_im%m22%elem(corr(sep_points(7)+1:sep_points(8))) = qrpa_broin(sep_points(7)+1:sep_points(8)) !imq
      end if

   end subroutine

   !---------------------------------------------------------------------------
   ! A wrapper for the more general modified Broyden mixing routine
   !---------------------------------------------------------------------------
   subroutine qrpa_broyden(niter, dRqp_re, dRqp_im)!nuv, x1, x2, x3, x4, x5, x6, x7, x8)
      implicit none
      integer, intent(in)  :: niter
      type(bigblockmatrix), intent(inout) :: dRqp_re, dRqp_im

      ! The Broyden mixing routine operates on a real array; Therefore we need
      ! to unpack the complex values to a real array, call the Broyden routine,
      ! and then repack the real and imaginary parts to the complex arrays
      call set_broinout_from_dRqp(dRqp_re,dRqp_im,'o')

      call broyden_method(nbroy,qrpa_broout,qrpa_broin,qrpa_alphamix,si,niter,&
          broyden_history_size,qrpa_bbroyden)

      call set_dRqp_from_broin(dRqp_re,dRqp_im)

   end subroutine qrpa_broyden

  !======================================================================================
  ! The following routine has been copy/pasted from HFBTHO v2.00d
  ! Stoitsov et al., Comp. Phys. Comm. 184 (2013) 1592
  !======================================================================================
  Subroutine broyden_method(N,vout,vin,alpha,si,iter,M,bbroyden)
    !---------------------------------------------------------------------
    ! Modified Broyden's method: D.D.Johnson, PRB 38, 12807 (1988)
    ! Adopted from: (C) 2001 PWSCF group
    ! Input :
    !  N      dimension of arrays vin,vout
    !  vin    outpu at previous iteration
    !  vout   output at current iteration
    !  alpha  mixing factor (0 < alpha <= 1)
    !  iter   current iteration number
    !  M      number of iterations in Broyden history
    !         (M=0 Linear mixing)
    ! Output:
    !  si     MaxVal(Abs(vout-vin))
    !  vin    Broyden/Linear mixing result
    !  vout   vout-vin
    !  bbroyden='B' Broyden mixing, curvature>0
    !  bbroyden='L' Linear mixing,  curvature<0
    !---------------------------------------------------------------------
    Implicit None
    Integer(ipr),    Intent(In)     :: N,iter,M
    Real(pr),        Intent(In)     :: alpha
    Real(pr),        Intent(Out)    :: si
    Character(1),    Intent(Out)    :: bbroyden
    Real(pr),        Intent(InOut)  :: vout(N),vin(N)
    Integer(ipr)                    :: iter_used,ipos,inext,info!,i,j
    Integer(ipr), Allocatable,Save  :: iwork(:)
    Real(pr),    Allocatable, Save  :: beta(:,:),beta_store(:,:),work(:)
    Real(pr),    Allocatable, Save  :: df(:,:),dv(:,:),curv(:),gamma(:)
    Real(pr),    Parameter          :: w0=0.010_pr
    Real(pr)                        :: DDOT,DNRM2,normi,sf,curvature
    !
    sf=-1.0_pr; Call DAXPY(N,sf,vin,1,vout,1) ! vout = vout - vin
    si=Maxval(Abs(vout))
    !---------------------------------------------------------------------
    ! No mixing
    !---------------------------------------------------------------------
    If (M < 0) Then
      bbroyden='N'
      sf=+1.0_pr; Call DAXPY(N,sf,vout,1,vin,1) ! vin = vin + vout' = vout
      Return
    End If
    !---------------------------------------------------------------------
    ! Linear mixing
    !---------------------------------------------------------------------
    If(M.Eq.0.Or.iter.Eq.0) Then
       bbroyden='L'; Call DAXPY(N,alpha,vout,1,vin,1) ! vin = vin + alpha*vout
       Return
    End If
    !---------------------------------------------------------------------
    ! Broyden mixing
    !---------------------------------------------------------------------
    bbroyden='B'
    iter_used=Min(iter-1,M)
    ipos=iter-1-((iter-2)/M)*M
    inext=iter-((iter-1)/M)*M
    If (iter.Eq.1) Then
       If(Allocated(curv)) Deallocate(curv,df,dv,beta,work,iwork,gamma)
       Allocate(curv(N),df(N,M),dv(N,M),beta(M,M),beta_store(M,M),work(M),iwork(M),gamma(M))
    Else
       df(:,ipos)=vout(:)-df(:,ipos)
       dv(:,ipos)= vin(:)-dv(:,ipos)
       Normi=1.0_pr/DNRM2(N,df(1,ipos),1)
       Call DSCAL(N,Normi,df(1,ipos),1)
       Call DSCAL(N,Normi,dv(1,ipos),1)
       Call DAXPY(N,alpha,df(1,ipos),1,dv(1,ipos),1) ! Vector u in the paper is stored in dv
       ! construct the upper triangle of beta: only update the row/col related to df(:,ipos)
       Call DGEMV('T',N,iter_used,1.0_pr,df(1,1),N,df(1,ipos),1,0.0_pr,work(1),1)
       beta_store(1:(ipos-1),ipos) = work(1:(ipos-1))
       beta_store(ipos,ipos) = w0*w0  + 1.0_pr
       If (iter-1>M) beta_store(ipos,(ipos+1):iter_used) = work((ipos+1):iter_used)
       beta(:,:) = beta_store(:,:)
    Endif
    Call DGEMV('T',N,iter_used,1.0_pr,df(1,1),N,vout(1),1,0.0_pr,gamma(1),1) ! Vector c in the paper
    Call DSYSV('U',iter_used,1,beta(1,1),M,iwork,gamma(1),M,work,M,info) ! gamma = beta^(-1)*c
    If(info.Ne.0) call abort(' In Broyden: info at DSYSV ')
    curv=alpha*vout
    Call DGEMV('N',N,iter_used,-1.0_pr,dv(1,1),N,gamma(1),1,1.0_pr,curv(1),1) ! curv = curv - gamma*u
    Call DCOPY(N,vout,1,df(1,inext),1)
    Call DCOPY(N,vin ,1,dv(1,inext),1)
   !  curvature=DDOT(N,vout,1,curv,1)
   !  If(.true.) then !curvature.Gt.-1.0_pr) Then
    bbroyden='B'; sf=+1.0_pr; Call DAXPY(N,sf,curv,1,vin,1) ! vin = vin + curv
   !  Else
      !  bbroyden='L'; sf=alpha*0.50_pr; Call DAXPY(N,sf,vout,1,vin,1) ! vin = vin + a/2*vout
   !  End If
  End Subroutine broyden_method

end module pnfam_broyden