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

contains

   subroutine set_broinout_from_dRqp(dRqp_re, dRqp_im, io)
      implicit none
      type(bigblockmatrix), intent(in) :: dRqp_re, dRqp_im
      character(1), intent(in) :: io
      integer :: nuv, nvec
      real(dp), pointer :: qrpa_bro(:)

      nullify(qrpa_bro)
      select case(io)
         case('i')
            qrpa_bro => qrpa_broin
         case('o')
            qrpa_bro => qrpa_broout
      end select
      nuv   = size(dRqp_re%m12%elem)
      nvec  = nbroy/nuv

      ! Make sure sizes match up
      if (allocated(dRqp_re%m11%elem)) then
          if (nvec/=8) call abort
      else
          if (nvec/=4) call abort
      end if

      qrpa_bro(1:nuv)         = dRqp_re%m12%elem !rex
      qrpa_bro(nuv+1:2*nuv)   = dRqp_re%m21%elem !rey
      qrpa_bro(2*nuv+1:3*nuv) = dRqp_im%m12%elem !imx
      qrpa_bro(3*nuv+1:4*nuv) = dRqp_im%m21%elem !imy

      if (allocated(dRqp_re%m11%elem)) then
          qrpa_bro(4*nuv+1:5*nuv) = dRqp_re%m11%elem !rep
          qrpa_bro(5*nuv+1:6*nuv) = dRqp_re%m22%elem !req
          qrpa_bro(6*nuv+1:7*nuv) = dRqp_im%m11%elem !imp
          qrpa_bro(7*nuv+1:8*nuv) = dRqp_im%m22%elem !imq
      end if

   end subroutine

   subroutine set_dRqp_from_broin(dRqp_re, dRqp_im)
      implicit none
      type(bigblockmatrix), intent(inout) :: dRqp_re, dRqp_im
      integer :: nuv, nvec

      nuv   = size(dRqp_re%m12%elem)
      nvec  = nbroy/nuv

      ! Make sure sizes match up
      if (allocated(dRqp_re%m11%elem)) then
          if (nvec/=8) call abort
      else
          if (nvec/=4) call abort
      end if

      dRqp_re%m12%elem = qrpa_broin(1:nuv)         !rex
      dRqp_re%m21%elem = qrpa_broin(nuv+1:2*nuv)   !rey
      dRqp_im%m12%elem = qrpa_broin(2*nuv+1:3*nuv) !imx
      dRqp_im%m21%elem = qrpa_broin(3*nuv+1:4*nuv) !imy

      if (allocated(dRqp_re%m11%elem)) then
          dRqp_re%m11%elem = qrpa_broin(4*nuv+1:5*nuv) !rep
          dRqp_re%m22%elem = qrpa_broin(5*nuv+1:6*nuv) !req
          dRqp_im%m11%elem = qrpa_broin(6*nuv+1:7*nuv) !imp
          dRqp_im%m22%elem = qrpa_broin(7*nuv+1:8*nuv) !imq
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
    Integer(ipr)                    :: i,j,iter_used,ipos,inext,info
    Integer(ipr), Allocatable,Save  :: iwork(:)
    Real(pr),    Allocatable, Save  :: beta(:,:),work(:)
    Real(pr),    Allocatable, Save  :: df(:,:),dv(:,:),curv(:)
    Real(pr),                 Save  :: w0
    Real(pr)                        :: DDOT,DNRM2,normi,gamma,sf,curvature
    !
    sf=-1.0_pr; Call DAXPY(N,sf,vin,1,vout,1) ! vout' = vout - vin
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
       w0=0.010_pr
       If(Allocated(curv)) Deallocate(curv,df,dv,beta,work,iwork)
       Allocate(curv(N),df(N,M),dv(N,M),beta(M,M),work(M),iwork(M))
    Else
       df(:,ipos)=vout(:)-df(:,ipos)
       dv(:,ipos)= vin(:)-dv(:,ipos)
       Normi=1.0_pr/sqrt(DNRM2(N,df(1,ipos),1)**2)
       Call dscal(N,Normi,df(1,ipos),1)
       Call dscal(N,Normi,dv(1,ipos),1)
    Endif
    Do i=1,iter_used
       Do j=i+1,iter_used
          beta(i,j)=DDOT(N,df(1, j),1,df(1,i),1)
       Enddo
       beta(i,i)= w0*w0  + 1.0_pr
    Enddo
    Call DSYTRF('U',iter_used,beta,M,iwork,work,M,info)
    If(info.Ne.0) call abort(' In Broyden: info at DSYTRF ')
    Call DSYTRI('U',iter_used,beta,M,iwork,work,info)
    If(info.Ne.0) call abort(' In Broyden: info at DSYTRI ')
    Do i=1,iter_used
       Do j=i+1,iter_used
          beta(j,i)=beta(i,j)
       Enddo
       work(i)=DDOT(N,df(1,i),1,vout,1)
    Enddo
    curv=alpha*vout
    Do i=1,iter_used
       gamma=0.0_pr
       Do j=1,iter_used
          gamma=gamma+beta(j,i)*work(j)
       Enddo
       curv=curv-gamma*(dv(:,i)+alpha*df(:,i))
    Enddo
    Call DCOPY(N,vout,1,df(1,inext),1)
    Call DCOPY(N,vin ,1,dv(1,inext),1)
    curvature=DDOT(N,vout,1,curv,1)
    If(.true.) then !curvature.Gt.-1.0_pr) Then
       bbroyden='B'; sf=+1.0_pr; Call DAXPY(N,sf,curv,1,vin,1) ! vin = vin + curv
    Else
       bbroyden='L'; sf=alpha*0.50_pr; Call DAXPY(N,sf,vout,1,vin,1) ! vin = vin + a/2*vout
    End If
  End Subroutine broyden_method

end module pnfam_broyden
