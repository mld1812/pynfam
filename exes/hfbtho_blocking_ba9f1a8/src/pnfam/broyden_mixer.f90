module broyden_mixer
   implicit none
   integer, parameter :: pr = kind(1d0)
   integer, parameter :: ipr = kind(1)
   character :: qrpa_bbroyden
   integer(ipr) :: broyden_history_size
   real(pr) :: si
   real(pr) :: qrpa_alphamix = 0.7
contains

  !======================================================================================
  ! The following routine has been copy/pasted from HFBTHO v2.00d
  ! Stoitsov et al., Comp. Phys. Comm. 184 (2013) 1592
  !======================================================================================
  Subroutine qrpa_broyden_method(N,vout,vin,alpha,si,iter,M,bbroyden)
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
    !  bbroyden='B' Broyden mixing, ='L' Linear mixing
    !---------------------------------------------------------------------
    Implicit None
    Integer(ipr),    Intent(In)     :: N,iter,M
    Real(pr),        Intent(In)     :: alpha
    Real(pr),        Intent(Out)    :: si
    Character(1),    Intent(Out)    :: bbroyden
    Real(pr),        Intent(InOut)  :: vout(N),vin(N)
    Integer(ipr)                    :: i,j,iter_used,ipos,inext
    Integer(ipr), Allocatable,Save  :: iwork(:)
    Real(pr),    Allocatable, Save  :: beta(:,:),work(:)
    Real(pr),    Allocatable, Save  :: df(:,:),dv(:,:),curv(:)
    Real(pr),                 Save  :: w0
    Real(pr)                        :: DDOT,DNRM2,normi,gamma,sf
    !
    sf=-1.0_pr; Call DAXPY(N,sf,vin,1,vout,1)         ! vout = vout - vin
    si=Maxval(Abs(vout))
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
    iter_used=Min(iter-1,M); ipos=iter-1-((iter-2)/M)*M; inext=iter-((iter-1)/M)*M
    If (iter.Eq.1) Then
       w0=0.010_pr
       If(Allocated(curv)) Deallocate(curv,df,dv,beta,work,iwork)
       Allocate(curv(N),df(N,M),dv(N,M),beta(M,M),work(M),iwork(M))
    Else
       df(:,ipos)=vout(:)-df(:,ipos)
       dv(:,ipos)= vin(:)-dv(:,ipos)
       Normi=1.0_pr/DNRM2(N,df(1,ipos),1)
       Call dscal(N,Normi,df(1,ipos),1)
       Call dscal(N,Normi,dv(1,ipos),1)
    Endif
    Do i=1,iter_used
       Do j=i+1,iter_used
          beta(i,j)=DDOT(N,df(1, j),1,df(1,i),1)
       Enddo
       beta(i,i)= w0*w0  + 1.0_pr
    Enddo
    Call DSYTRF('U',iter_used,beta,M,iwork,work,M,i)
    If(i.Ne.0) Stop '  In Broyden: info at DSYTRF '
    Call DSYTRI('U',iter_used,beta,M,iwork,work,i)
    If(i.Ne.0) Stop '  In Broyden: info at DSYTRI '
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
    Call dcopy(N,vout,1,df(1,inext),1)
    Call dcopy(N,vin ,1,dv(1,inext),1)
    sf=+1.0_pr; Call DAXPY(N,sf,curv,1,vin,1)
  End Subroutine qrpa_broyden_method

end module
