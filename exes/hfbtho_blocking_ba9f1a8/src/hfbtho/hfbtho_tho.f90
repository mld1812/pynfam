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
!           TRANSFORMED HARMONIC OSCILLATOR BASIS PACKAGE              !
!                                                                      !
! ==================================================================== !

!----------------------------------------------------------------------
!> This module provides the main \theCode DFT solver. It includes
!> routines for the calculation and diagonalization of the HFB matrix;
!> the definition of densities on the quadrature mesh and in configuration
!> space; the self-consistent loop the calculation of expectation values of
!> observables. It also includes a package of a few routines to perform
!> particle number projection.
!>
!>  @author
!>    Mario Stoitsov, Nicolas Schunck
!----------------------------------------------------------------------
!  Subroutines: - f01234(lpr)
!               - thofun(key,r,f,f1,f2,f3,fj,lpr,units)
!               - densitr(it,xr,yr,yrP,msw)
!               - gaupolr(z,x)
!----------------------------------------------------------------------
Module HFBTHO_THO
  Use HFBTHO_utilities
  Use math
  Implicit None
  !
  Character(6), Private :: THO_version='1'
  !
Contains
  !
  !=======================================================================
  !>  Calculates LST-function \f$ f(R) \f$ its 1st through 4th derivatives
  !> 'f1,f2,f3,f4', the Jacobian
  !>    \f[
  !>       fd = \left( f^2(r) \frac{1}{r^2}\frac{df}{dr} \right)^{1/2}
  !>    \f]
  !>  and its first and second derivatives 'fd1' and 'fd2' at the point 'r'.
  !>  The LST function is defined as in routine \ref thofun.
  !>
  !>  All integrations are performed in the dimensionless variables
  !>  'u=xh(ih)' and 'v=xl(il)' being the dimensionless Gauss mesh
  !>  points in cylindrical coordinates. Thus 'r' is found as the inverse
  !>  function of the input variable \f$ qq= \sqrt{xh(ih)^2+xl(il)}', i.e.,
  !>  solving f(r)=qq. The original 'z' and 'rho' variables in 'fm' are
  !>  defined as 'fh(ihli) = zz*bz' and 'fl(ihli) = rr*bp'.
  !>
  !>  All required quantities for evaluating w.f., its first derivative
  !>  and the associated Laplacian are calculated
  !>
  !>  Definitions used:
  !>   - xh(ih),xl(il) : Gauss mesh points
  !>   - bp, bz        : oscillator lengths in fm
  !>   - ilst1         : 0->ho, #0->tho
  !>
  !>  In this way:
  !>   - \delta \rho divided by j^2 is:
  !>         dro = (cz*fs1+cr*fs2+fs3)*qhab*qlab
  !>             +  fs1*qha1b1*qlab
  !>             +  fs2*qhab*qla1b1/(4.*v**2)
  !>             +  fs4*(qha1b+qhab1)*(qla1b+qlab1)/(2.*v)
  !>             +  fs5*(qha1b+qhab1)*qlab
  !>             +  fs6*qhab*(qla1b+qlab1)/(2.*v)
  !>      with
  !>         cr = 1./4. - (nr_a+nr_b+m+1)/(2.*v) + (m/(2.*v))**2
  !>         cz = u*u - (nz_a+nz_b+1)
  !>   - the first (r,z) derivatives are:
  !>         fi'_z = fp1*(qh*ql) + fp2*(qh1*ql) + fp3*(qh*ql1)/(2.*v)
  !>               -> HO -> qh1*ql/bpz
  !>         fi'_r = fp4*(qh*ql) + fp5*(qh1*ql) + fp6*(qh*ql1)/(2.*v)
  !>               -> HO -> qh*ql1/Sqrt(v)/bp
  !=======================================================================
  Subroutine f01234(lpr)
    Use HFBTHO
    Implicit None
    Logical :: lpr
    Integer(ipr), Save :: i,il,ih,ihli,iw,key0=0,key1=1
    Real(pr), Save :: bri,bri2,bzi,bzi2,wv1,u,qq,f,f1,f2,f3,rhoi,fd1,fd2,         &
                      fd12,fdd,r,r2,r3,r4,r5,r6,rr,rr2,rr4,zz,zz2,drr1,drr2,drz1, &
                      drz2,drr12,drz12,g,g1,g2,g3,gg,gg1,gg2,g1g1,uz1,ur1,vz1,vr1,&
                      uz2,ur2,vz2,vr2,ur12,vr12,uz12,vz12
    !
    bri= one/bp; bri2= bri*bri; bzi= one/bz; bzi2= bzi*bzi
    Do il=1,ngl
       wv1   = xl(il)
       Do ih=1,ngh
          ihli = ih+(il-1)*ngh
          u    = xh(ih); qq   = Sqrt(u*u+wv1)
          If(ilst1.Eq.0) Then
             ! ho-case
             r = qq; f = r; f1 = one; f2 = zero; f3 = zero
          Else
             ! tho-case: initial run
             If(ih*il.Eq.1.And.ilst.Lt.0) &
                Call thofun(key0,g,f,f1,f2,f3,g1,.True.,.False.)
             If(iasswrong(3).ne.0) then
                ! reinforce ho results
                r = qq; f = r; f1 = one; f2 = zero; f3 = zero
             else
                ! tho-case: f(r)=qq,f'(r),f''(r),f'''(r), r=Invers_f(r)
                Call thofun(key1,qq,f,f1,f2,f3,r,.False.,.False.)
             End If
          End If
          ! Jacobian calculations
          r2= r*r; r3= r2*r; r4= r2*r2 ;r5= r3*r2 ;r6= r3*r3
          ! fdd=(f(r)^2 f'(r)/r^2)^(1/2),fd1=fdd'/fdd, fd2=fdd''/fdd
          fd1  = f1/f - one/r + half*f2/f1
          fdd  = f*Sqrt(f1)/r
          fd2  = two*(f2/f-fd1/r) - (f2/f1)**2/four + half*f3/f1
          fd12 = fd1**2
          ! g=(f/r)-derivatives
          g    = f/r; g1   =-(f - f1*r)/r2
          g2   = (two*(f - f1*r) + f2*r2)/r3
          g3   = (six*(f1*r - f) - three*f2*r2 + f3*r3)/r4
          ! g4   = (24.0d0*(f-f1*r+half*f2*r2)-four*f3*r3+f4*r4)/r5
          gg   = g*g; gg1= g*g1; gg2= g*g2; g1g1= g1*g1
          ! (rr,zz)-definitions
          rr   = Sqrt(wv1)/g; rhoi = bri/rr; zz   = u/g
          rr2  = rr*rr; rr4= rr2*rr2 ;zz2= zz*zz
          ! (r,z)-derivatives
          drr1  = bri*rr/r;           drz1  = bzi*zz/r
          drr2  = (bri2 - drr1**2)/r; drz2  = (bzi2 - drz1**2)/r
          drr12 = drr1**2;            drz12 = drz1**2
          ! (u,v)-derivatives
          uz1  = bzi*(g+g1*zz2/r);    ur1  = bri*g1*rr*zz/r
          vz1  = bzi*two*gg1*rr2*zz/r
          vr1  = bri*two*rr*(gg+gg1*rr2/r)
          uz2  = bzi2*zz*(three*g1*r2-g1*zz2+g2*r*zz2)/r3
          ur2  = bri2*(g1*r2-g1*rr2+g2*r*rr2)*zz/r3
          vz2  = bzi2*two*rr2*(gg1*r2-gg1*zz2+g1g1*r*zz2+gg2*r*zz2)/r3
          vr2  = gg*r3+five*gg1*r2*rr2-gg1*rr4+g1g1*r*rr4+gg2*r*rr4
          vr2  = vr2*bri2*two/r3; ur12=ur1**2; vr12=vr1**2; uz12=uz1**2
          vz12 = vz1**2
          ! storage
          fh(ihli)= zz*bz; fl(ihli)= rr*bp; fli(ihli)=one/fl(ihli); fd(ihli)= fdd*fdd
          ! for first derivatives
          fp1(ihli)= fd1*drz1; fp2(ihli)= uz1; fp3(ihli)= vz1
          fp4(ihli)= fd1*drr1; fp5(ihli)= ur1; fp6(ihli)= vr1
          ! for the Laplacian
          fs1(ihli) = two*(ur12 + uz12); fs2(ihli) = two*(vr12 + vz12)
          fs3(ihli) = two*(fd1*(drr1*rhoi+drr2+drz2)+ &
                      (fd12+fd2)*(drr12+drz12))
          fs4(ihli) = two*(ur1*vr1 + uz1*vz1)
          fs5(ihli) = four*fd1*(drr1*ur1 + drz1*uz1) + &
                      ur1*rhoi + ur2 + uz2
          fs6(ihli) = four*fd1*(drr1*vr1 + drz1*vz1) + &
                      vr1*rhoi + vr2 + vz2 - (vr12 +vz12)/wv1
       End Do   !ihs
    End Do   !il
    !
    ! Associated (z,r)-weights
    Do il = 1,ngl
       Do ih = 1,ngh
          i = ih + (il-1)*ngh
          wdcor(i) = pi*wh(ih)*wl(il)*bz*bp*bp/fd(i)
          wdcori(i)=one/wdcor(i)
       End Do
    End Do
    !
    If(lpr) Then
       Do iw=lout,lfile
          Write(iw,*)
          If(ilst1.Eq.0) Then
             Write(iw,*) ' ### HO case: wdcor charged'
          Else
             Write(iw,*) ' ### THO case: wdcor charged'
          End If
          Write(iw,*)
       End Do
    End If
    Return
  End Subroutine f01234
  !=======================================================================
  !> Calculates LST-function 'f' and its derivatives 'f1,f2,f3'
  !> at the point 'r' (all dimensionless).
  !=======================================================================
  Subroutine thofun(key,r,f,f1,f2,f3,fj,lpr,units)
    Use HFBTHO
    Implicit None
    Logical :: lpr,units
    Integer(ipr) :: key,msw
    Integer(ipr) :: it,iter,ir,iqq,irmax,irmsit,immho,imm1,&
                    imm2,immm,immmax,imm3
    Real(pr) ::  r,f,f1,f2,f3,fj
    Real(pr), Allocatable :: dsx(:),dsy(:),dsyT(:),dsyi(:),dsyii(:),dsy1(:),         &
                             dsy1i(:),spb0(:),spc0(:),spd0(:),spbi(:),spci(:),spdi(:)
    Real(pr) :: h,hhb,pihhb,c00,snorm,snorm1,assm,asm1,asm2,asm3,               &
                rmsit,rmmho,z1,z10,aaa,bex,rend,fj1,fj2,fj3,                    &
                s,s1,sN,sP,sT,qq,qqup,qqdn,zqq,zqqi,df,zfj1,zfj1i,fjb,aa,bb,yyy,&
                sqsq,rmmm,rmmmb0,z1mmm,rmm1,rmm1b0,z1mm1,rmm2,rmm2b0,z1mm2,     &
                rmmmax,z1mmmax,rmmmaxb0,rmmx,z1mmx,z1mmxx,alaex,aldsy1,decay2
    Real(pr) :: denm1(2),denm2(2),rdenm(2)
    Real(pr) :: epsf=1.0d-14,epsdsy=1.0d-16,epsnorm
    Complex(pr) :: yyy1,bbb1,aac !for the e3rd order equation
    !

    ! Adjustable parameters
    ! Rend=40.0_pr
    epsnorm=0.01_pr

    ! ===========================
    ! KEY=0 INITIAL CALCULATIONS
    ! ===========================
    qq = r
    If(key.Eq.0.And.ilst.Le.0) Then
       write(*,*)
       write(*,*) ' LST transformation...'
       !
       ! steps
       h = 0.01_pr;
       hhb = b0*h; If(units) hhb = h
       pihhb = 4.0_pr*pi*hhb
       !
       ! correct density asymptotic
       !
       !================================
       ! Test neutron/proton asymptotic
       !================================
       Do it=1,2
          ! Neutron/proton density decay constant
          itass=it; decay=ass(itass); decay2=decay**2;
          rmsit=rms(itass)+one; irmsit=Int(rmsit/hhb)
          bb=Real((itass-1)*npr(2),Kind=pr)*Sqrt(1.440d0)/hb0
          ! correct density Rend
          Rend=40.0d0; irmax=Int(Rend/hhb)
          ! Deallocate/Allocate
          if(Allocated(dsx)) Deallocate(dsx,dsy,dsyT,dsyi,dsyii,dsy1,dsy1i,&
                                        spb0,spc0,spd0,spbi,spci,spdi)
          Allocate(dsx(irmax),dsy(irmax),dsyT(irmax),dsyi(irmax),dsyii(irmax),  &
                   dsy1(irmax),dsy1i(irmax),spb0(irmax),spc0(irmax),spd0(irmax),&
                   spbi(irmax),spci(irmax),spdi(irmax))
          ! HFB+HO_{L=0} density 'dsy' and its normalization
          ! integral 'dsyi' at points 'dsx' with step 'hhb=h*b0'
          ! up to the point where 'dsy*dsx*dsx < epsdsy'
          msw=0; snorm=zero
          Do While(Abs(snorm-Real(npr(itass),Kind=pr)).Gt.epsnorm.And.msw.Lt.25)
             msw=msw+6 ! increase for good norm of HFB+HO_{L=0}
             itass = - itass; s1 = zero
             Do ir=1,irmax
                rmmho=hhb*Real(ir,Kind=pr); immho=ir
                ! L=0 component of density for isospin it (s) and isospin 1-it (sT)
                Call densitr(itass,rmmho,sN,sP,msw)
                if(itass.eq.1) then
                   s=sN; sT=sP
                else
                   s=sP; sT=sN
                End If
                s=s*Dnfactor(itass)
                z1=s*rmmho**2; s1=s1+z1
                ! density dsy(ir) at point ir and its integral over r dsyi(ir)
                ! up to that point
                dsyT(ir)=sT; dsy(ir)=s; dsyi(ir)=hhb*s1
                immho=ir !up to the point immho
                If(z1.Lt.epsdsy) Exit
             End Do
             snorm=pihhb*s1  ! HFB+HO_{L=0} norm
          End Do
          ! dsy: density, spb0: first derivative with respect to r
          Call deri(hhb,immho,dsy,spb0)
          ! MIN: Find 'rmm1', the first minimum of Ln(HFB+HO_{L=0})'
          z10 = 1.0d10
          Do ir=irmsit,immho-5
             denm1(it)=dsy(ir); denm2(it)=dsyT(ir)
             z1 = spb0(ir)/dsy(ir)
             If(z1.Le.z10) Then
                imm1 = ir; rmm1 = hhb*Real(ir,Kind=pr);  z1mm1= z1; Else; Exit
             End If
             z10 = z1
          End Do
          rdenm(itass) = rmm1/b0
          If(units) rdenm(itass) = rmm1
          ! no minimum of Ln(HFB+HO_{L=0})'
          If(rmm1.Ge.hhb*Real(immho-5,Kind=pr)) Then
             Write(*,*)
             Write(*,*) '#####################################'
             Write(*,*) 'Please increase Nsh NB!!!(NO THO RUN)'
             Write(*,*) '#####################################'
             Write(*,*)
             iasswrong(itass)=-1
             Stop
             !If(lpr) Then
             !   Open(1110,file='dat0.dat')
             !   Write(1110,*) ' r rhoh Log(rhoh)'
             !   Do ir=5,immho-5
             !      Write(1110,'(14(4x,e13.6))') hhb*Real(ir),dsy(ir),spb0(ir)/dsy(ir)
             !   End Do
             !   Close(1110)
             !End If
             !Return
          End If
       End do
       ! Asymptotics - denm1: density for isospin it, denm2: density for isospin 1-it
       If((denm1(1)-denm2(1))*(denm2(2)-denm1(2)).Le.zero) Then
          itass=1; if(ass(1).gt.ass(2)) itass=2      ! mismatch: use old asymptotic (lower decay)
       Else
          itass=2; If(denm1(1).gt.denm2(1)) itass=1  ! use new asymptotic (higher density)
       End If
       !
       iasswrong(3)=0
       If(iasswrong(itass).ne.0) iasswrong(3)=iasswrong(itass)  ! wrong assymptotic => reinforce HO results
       !
       Write(*,*) '                          min.point         neutron density       proton density'
       Write(*,*) '  1. Neutron min.point ',rdenm(1),denm1(1),denm2(1)
       Write(*,*) '  2. Protons min.point ',rdenm(2),denm2(2),denm1(2)
       Write(*,*) '     Neutron/Proton decay',ass(1),'/',ass(2)
       Write(*,*) '     Chosen Case=',itass
       !
       !===================
       ! Actual asymptotic
       !===================
       ! neutron/proton density decay constant
       decay = ass(itass); decay2=decay**2;
       rmsit= rms(itass)+one;irmsit=Int(rmsit/hhb)
       bb = Real((itass-1)*npr(2),Kind=pr)*Sqrt(1.440d0)/hb0
       ! correct density Rend
       Rend  = 40.0d0; irmax = Int(Rend/hhb)
       ! Deallocate/Allocate
       If(Allocated(dsx)) Deallocate(dsx,dsy,dsyT,dsyi,dsyii,dsy1,dsy1i,&
                                     spb0,spc0,spd0,spbi,spci,spdi)
       Allocate(dsx(irmax),dsy(irmax),dsyi(irmax),dsyii(irmax),   &
                dsy1(irmax),dsy1i(irmax),spb0(irmax),spc0(irmax), &
                spd0(irmax),spbi(irmax),spci(irmax),spdi(irmax))
       ! HFB+HO_{L=0} density 'dsy' and its normalization
       ! integral 'dsyi' at points 'dsx' with step 'hhb=h*b0'
       ! up to the point where 'dsy*dsx*dsx < epsdsy'
       msw = 0; snorm=zero
       Do While(Abs(snorm-Real(npr(itass),Kind=pr)).Gt.0.01.And.msw.Lt.25)
          msw = msw + 6 ! increase for good norm of HFB+HO_{L=0}
          itass = - itass; s1 = zero
          Do ir=1,irmax
             rmmho = hhb*Real(ir,Kind=pr); immho = ir
             Call densitr(itass,rmmho,sN,sP,msw)
             if(itass.eq.1) then
                s=sN; sT=sP
             else
                s=sP; sT=sN
             End If
             s=s*Dnfactor(itass)
             z1 = s*rmmho**2; s1 = s1 + z1
             dsy(ir)= s; dsyi(ir)= hhb*s1 !p-ho density and its integral
             immho = ir                   !up to the point immho
             If(z1.Lt.epsdsy) Exit
          End Do
          snorm = pihhb*s1  ! HFB+HO_{L=0} norm
       End Do
       Call deri(hhb,immho,dsy,spb0)
       ! MIN: Find 'rmm1', the first minimun of Ln(HFB+HO_{L=0})'
       z10 = 1.0d10
       Do ir=irmsit,immho-5
          z1 = spb0(ir)/dsy(ir)
          If(z1.Le.z10) Then
             imm1 = ir; rmm1 = hhb*Real(ir,Kind=pr);  z1mm1= z1; Else; Exit
          End If
          z10 = z1
       End Do
       rmm1b0 = rmm1/b0
       If(units) rmm1b0 = rmm1
       ! MAX: Find 'rmmmax', the first maximum of ln(HFB+HO_{L=0})'
       z10 = z1mm1
       Do ir=imm1,immho-5
          z1 = spb0(ir)/dsy(ir)
          If(z1.Ge.z10) Then
             immmax = ir; rmmmax = hhb*Real(ir,Kind=pr);  z1mmmax= z1; Else; Exit
          End If
          z10 = z1
       End Do
       rmmmaxb0 = rmmmax/b0
       If(units) rmmmaxb0 = rmmmax
       !
       ! END: Find 'rmm2', the last point of ln(HFB+HO_{L=0}) at the level of the first minimum
       z10 = z1mm1
       Do ir=immmax,immho-5
          z1 = spb0(ir)/dsy(ir)
          If(z1.Ge.z10) Then
             imm2 = ir; rmm2 = hhb*Real(ir,Kind=pr);  z1mm2= z1
          End If
       End Do
       rmm2b0 = rmm2/b0
       If(units) rmm2b0 = rmm2
       !
       ! MID: Find 'rmmm', the MinMaX mid point
       z10 = half*(z1mmmax+z1mm1)
       Do ir=imm1,immho
          z1 = spb0(ir)/dsy(ir)
          If(z1.Ge.z10) Exit
          immm = ir; rmmm = hhb*Real(ir,Kind=pr);  z1mmm= z1
       End Do
       rmmmb0 = rmmm/b0
       If(units) rmmmb0 = rmmm
       ! -------------------------------------------------------------
       ! Important points required:
       ! Minimum point         'rmm1'   and its log.density 'z1mm1'
       ! First maximum         'rmmmax' and its log.density 'z1mmmax'
       ! Last acceptable point 'rmm2'   and its log.density 'z1mm2'
       ! Mid point             'rmmm'   and its log.density 'z1mmm'
       ! -------------------------------------------------------------
       ! fit 'aa' from the mid match point 'rmmm','z1mmm'
       If(z1mmm.Ge.-decay) z1mmm= (-decay+z1mm1)/two   !just in case
       ! the 3rd order equation
       sqsq = rmmm*(two*bb+decay2*rmmm)
       bbb1 = one + rmmm*z1mmm
       yyy1 =(2.0d0*bbb1**6 + 18.0d0*bbb1**3*sqsq + 27.0d0*sqsq**2 + &
            3.0d0*Sqrt(3.0d0)*sqsq**1.5d0*Sqrt(4.0d0*bbb1**3 + &
            27.0d0*sqsq))**(1.0d0/3.0d0)
       aa = Real((2.*bbb1**2*yyy1 + 2.**(2.0d0/3.0d0)*yyy1**2 - &
            2.*2.**(1.0d0/3.0d0)*(-bbb1**4 - &
            6.0d0*bbb1*sqsq))/(6.*yyy1),Kind=pr)
       aa = (-4.0d0*bb*rmmm - decay2*rmmm**2 + aa)*0.250d0
       !write(*,*)  ' aa= ',aa !If(aa.Le.zero) aa   = zero  ! in this case take l=0
       !
       ! matching logder at Rmin='rmm1' and 'rmmx'
       ! log density for rmm1<r<rmmx is z1mm1 + aaa*(r-rmm1)**2/r**assm
       ! rmmx  = rmm2     !rmm2 is last acceptable point
       rmmx  = rmmmax      !rmmmax is the firts maximum point
       sqsq  = Sqrt(decay2 + (4.0d0*(aa + bb*rmmx))/rmmx**2)
       z1mmx=-one/rmmx - sqsq - (two*bb + decay2*rmmx)/&
            (4.0d0*aa + rmmx*(4.0d0*bb + decay2*rmmx))
       z1mmxx=(two*(8.0d0*aa**2 + 16.0d0*aa*bb*rmmx + &
            12.0d0*bb**2*rmmx**2 + two*aa*decay2*rmmx**2 + &
            6.0d0*bb*decay2*rmmx**3 + decay2**2*rmmx**4 + &
            rmmx*(two*aa + bb*rmmx)*sqsq*(4.0d0*aa + &
            rmmx*(4.0d0*bb + decay2*rmmx))))/(rmmx**2*(4.0d0*aa + &
            rmmx*(4.0d0*bb + decay2*rmmx))**2)
       assm  = rmmx*(two/(rmmx-rmm1)-z1mmxx/(z1mmx-z1mm1))
       asm1  =  one-assm; asm2= two-assm; asm3= three-assm
       aaa   = -rmmx*rmmx**assm*z1mmxx/ &
            ((rmm1-rmmx)*(assm*(rmm1-rmmx)+two*rmmx))
       alaex = Log(dsy(imm1))-&
            (z1mm1*rmm1 +two*aaa*rmm1**asm3/(asm3*asm2*asm1))
       ! correct density 'dsy1'
       Do ir=1,irmax
          qq = hhb*Real(ir,Kind=pr); dsx(ir) = qq
          If(ir.Le.imm1) Then
             ! inner region (r < rmm1)
             dsy1(ir) = dsy(ir)
          Else
             sqsq=Sqrt(decay2 + 4.0d0*(aa/qq**2 + bb/qq))
             yyy=-qq*sqsq - Log(qq*qq*sqsq) - &
                  (two*bb/decay)*Log((two*bb + decay2*qq+decay*qq*sqsq)/decay)
             ! take complex in the case of negative aa
             aac=aa
             yyy=yyy+two*Sqrt(aac)*&
                  Log((two*aac+bb*qq+Sqrt(aac)*qq*sqsq)/(two*aac**1.50d0*qq))
             If(qq.Le.rmmx) Then
                ! region (rmm1 < r < rmmx)
                aldsy1 = alaex+z1mm1*qq+aaa*qq**asm1*&
                     (asm3*asm2*rmm1**2-two*asm3*asm1*rmm1*qq+asm2*asm1*qq*qq)/ &
                     (asm3*asm2*asm1)
                dsy1(ir)=Exp(aldsy1)
                bex = aldsy1-yyy
             Else
                ! region (r > rmmx)
                dsy1(ir) = Exp(bex + yyy)
             End If
          End If
       End Do
       ! correct density norm
       snorm1 = pihhb*Sum(dsy1*dsx*dsx)
       ! normalized correct density 'dsy1'
       dsy1   = snorm*dsy1/snorm1
       ! zero constant
       c00 = (dsy1(1)/dsy(1))**(1.0d0/3.0d0)
       ! splining correct density and its integral
       s1 = zero
       Do ir=1,irmax
          s1 = s1 +  dsy1(ir)*dsx(ir)**2
          dsy1i(ir) = hhb*s1
       End Do
       !
       ! correct density dsy1 and its integral dsy1i known up to irmax
       Call csplin(irmax,dsx,dsy1 ,spb0,spc0,spd0)
       Call csplin(irmax,dsx,dsy1i,spbi,spci,spdi)
       !
       ! print 'dat1.dat' with HFB+HO and 'correct' densities
       ! and their Log derivatives at lpr=.true.
       If(lpr) Then
          Open(1110,file='density.dat')
          !Write(1110,*) ' r rhoh rhoc Log(rhoh)'' Log(rhoc)'' '
          !Do ir=5,immho-5
          Do ir=5,irmax-5
             ! ho density derivative
             s =(45.0d0*( dsy(ir+1)-dsy(ir-1))-9.0d0*&
                  (dsy(ir+2)-dsy(ir-2))+dsy(ir+3)-dsy(ir-3))/(60.0d0*hhb)/dsy(ir)
             !correct density derivative
             s1=(45.0d0*(dsy1(ir+1)-dsy1(ir-1))-9.0d0*(dsy1(ir+2)-&
                  dsy1(ir-2))+dsy1(ir+3)-dsy1(ir-3))/(60.0d0*hhb)/dsy1(ir)
             Write(1110,'(2(1x,e13.6))') dsx(ir),dsy1(ir)
          End Do
          Close(1110)
          !Open(1111,file='dat1.dat')
          !Write(1111,*) ' Dimensionless_qq  Invers_f Invers_f1 Invers_f2 Invers_f3 '
       End If
       !
       ! =======================================
       ! Calculations at given dimensionless 'qq'
       ! =======================================
       ! f(R->0) = c00*R, therefore Invers_f(R)=R/c00
       fj = h/c00
       bmm3= zero; z1 = zero; ir=0
       Do iqq=1,immho
          qq = Real(iqq,Kind=pr)*h
          ! HFB+HO density and integral at 'b0*qq' or 'qq'
          zqq = dsy(iqq); zqqi = dsyi(iqq)
          ! Iterations to find 'fj = Invers_f(qq)'
          ! NB! f[Invers_f(qq)]=qq
          iter= 0; df= 0.00010d0
          Do While(Abs(df).Ge.epsf.And.iter.Le.500)
             iter = iter + 1
             fjb = fj*b0; If(units) fjb=fj
             Call cseval(irmax,fjb,dsx,dsy1 ,spb0,spc0,spd0,zfj1 )
             Call cseval(irmax,fjb,dsx,dsy1i,spbi,spci,spdi,zfj1i)
             qqup = (Log(zfj1i/zqqi))*zfj1i;
             qqdn = zfj1*b0*fjb**2; If(units) qqdn=zfj1*fjb**2
             ! Secant & Newton
             If(zfj1i.Le.zqqi.And.df.Le.0.0d0) df=-half*df
             If(zfj1i.Gt.zqqi.And.df.Gt.0.0d0) df=-half*df
             If(Abs(qqdn).Gt.Abs(qqup).And.iter.Le.20) df= - qqup/qqdn
             fj = fj + df
          End Do
          fj1 = (zqq*qq*qq)/(zfj1*fj*fj)
          fj2 = (fj1 - z1)/h; z1 =fj1
          If(qq.Gt.rmm2b0) Then
             If(fj1.Ge.bmm3) Then
                bmm3 = fj1; Else; Exit
             End If
          End If
          dsy(iqq)   = fj
          dsyi(iqq)  = fj1
          dsyii(iqq) = fj2
          iqqmax     = iqq
       End Do
       imm3 = iqqmax-50
       rmm3 = Real(imm3,Kind=pr)*h
       amm3 = dsy(imm3)
       bmm3 = dsyi(imm3)
       cmm3 = dsyii(imm3)
       ! second and third derivatives up to 'iqqmax''
       Call deri(h,iqqmax,dsyi,spb0)
       Call deri(h,iqqmax,spb0,spc0)
       !
       If(Allocated(fdsx)) Deallocate(fdsx,fdsy,fdsy1,fdsy2,fdsy3,fspb0,fspc0,  &
                                      fspd0,fspb1,fspc1,fspd1,fspb2,fspc2,fspd2,&
                                      fspb3,fspc3,fspd3)
       Allocate(fdsx(iqqmax),fdsy(iqqmax),fdsy1(iqqmax),fdsy2(iqqmax),  &
                fdsy3(iqqmax),fspb0(iqqmax),fspc0(iqqmax),fspd0(iqqmax),&
                fspb1(iqqmax),fspc1(iqqmax),fspd1(iqqmax),fspb2(iqqmax),&
                fspc2(iqqmax),fspd2(iqqmax),fspb3(iqqmax),fspc3(iqqmax),&
                fspd3(iqqmax))
       Do iqq=1,iqqmax
          fdsx(iqq)  = Real(iqq,Kind=pr)*h
          fdsy(iqq)  = dsy(iqq)
          fdsy1(iqq) = dsyi(iqq)
          fdsy2(iqq) = spb0(iqq)
          fdsy3(iqq) = spc0(iqq)
          !
          ! print 'dat1.dat' with fj=Inverse_f(qq) and its derivatives
          ! fj1..3 at no smoothing when lpr=.true.
          !If(lpr) Then
          !   Write(1111,'(14(4x,e13.6))') fdsx(iqq)*b0,fdsy(iqq),&
          !        fdsy1(iqq),fdsy2(iqq),fdsy3(iqq)
          !End If
       End Do
       !
       If(Allocated(dsx)) Deallocate(dsx,dsy,dsyi,dsyii,dsy1,dsy1i,&
                                     spb0,spc0,spd0,spbi,spci,spdi)
       !
       Call csplin(iqqmax,fdsx,fdsy ,fspb0,fspc0,fspd0)
       Call csplin(iqqmax,fdsx,fdsy1,fspb1,fspc1,fspd1)
       Call csplin(iqqmax,fdsx,fdsy2,fspb2,fspc2,fspd2)
       Call csplin(iqqmax,fdsx,fdsy3,fspb3,fspc3,fspd3)
       !
       Do ir=lout,lfile
          Write(ir,*)
          Write(ir,*) ' Legendre points = ',msw
          Write(ir,*) ' b0, decay=        ',b0,decay
          Write(ir,*) ' h, hhb=           ',h,hhb
          Write(ir,*) ' rms, rmsit=       ',rms(itass),rmsit
          Write(ir,*) ' Rend,   irmax=    ',Rend,irmax
          Write(ir,*) ' HORend, immho=    ',rmmho,immho
          Write(ir,*) ' snorm,snorm1=     ',snorm,snorm1
          Write(ir,*) ' snorm/snorm1=     ',snorm/snorm1
          Write(ir,*) ' min:rmm1,/b0=     ',rmm1,rmm1b0
          Write(ir,*) ' max:rmmmax,/b0=   ',rmmmax,rmmmaxb0
          Write(ir,*) ' last:rmm2,/b0=    ',rmm2,rmm2b0
          Write(ir,*) ' num:rmm3*b0,rmm3= ',rmm3*b0,rmm3
          Write(ir,*) ' rmmho, rmmho/b0=  ',rmmho,rmmho/b0
          Write(ir,*) ' alaex,bex=        ',alaex,bex
          Write(ir,*) ' aa,bb=            ',aa,bb
          Write(ir,*) ' L_eff=            ',(Sqrt(one + 4.0d0*aac)-one)/two
          Write(ir,*) ' amm3, bmm3=       ',amm3,bmm3
          Write(ir,*) ' cmm3, one/c00=    ',cmm3,one/c00
          Write(ir,*) ' z1mm1,z1mmm=      ',z1mm1,z1mmm
          Write(ir,*)
       End Do
       !! print 'dat2..4.dat' after when lpr=.true.
       !! dat2.dat: 'r=b0*qq',correct density
       !! dat3.dat   qq, Invers_f,Invers_f1...3
       !! dat4.dat   qq, f,f1,f2,f3; 'qq' are the Gauss points
       !If(lpr) Then
       !   Close(1111)
       !   Open(1112,file='dat2.dat')
       !   Open(1113,file='dat3.dat')
       !   Open(1114,file='dat4.dat')
       !   Write(1112,*) ' r  den_correct'
       !   Write(1113,*) ' qq Invers_f Invers_f1 Invers_f2 Invers_f3'
       !   Write(1114,*) ' qq f f1 f2 f3'
       !   Do iqq=1,ngh
       !      Do ir=1,ngl
       !         qq   = Sqrt(xh(iqq)**2+xl(ir))
       !         Call densitr(itass,b0*qq,sN,sP,msw)
       !         if(itass.eq.1) then
       !            s=sN; sT=sP
       !         else
       !            s=sP; sT=sN
       !         End If
       !         s=s*Dnfactor(itass)
       !         If(qq.Le.rmm3) Then
       !            Call cseval(iqqmax,qq,fdsx,fdsy ,fspb0,fspc0,fspd0,fj)
       !            Call cseval(iqqmax,qq,fdsx,fdsy1,fspb1,fspc1,fspd1,fj1)
       !            Call cseval(iqqmax,qq,fdsx,fdsy2,fspb2,fspc2,fspd2,fj2)
       !            Call cseval(iqqmax,qq,fdsx,fdsy3,fspb3,fspc3,fspd3,fj3)
       !         Else
       !            fj  = amm3+bmm3*(qq-rmm3)+cmm3*(qq-rmm3)**2/two
       !            fj1 = bmm3+cmm3*(qq-rmm3)
       !            fj2 = cmm3; fj3 = zero
       !         End If
       !         f  = qq; f1 = one/fj1; f2 =-fj2*f1**3
       !         f3 = three*fj2**2*f1**5-fj3*f1**4
       !         s = s*qq*qq/(fj*fj*fj1)
       !         Write(1112,'(14(4x,e13.6))') b0*fj,s
       !         Write(1113,'(14(4x,e13.6))') qq,fj,fj1,fj2,fj3
       !         Write(1114,'(14(4x,e13.6))') qq,f,f1,f2,f3
       !      End Do
       !   End Do
       !   Close(1112); Close(1113);
       !   Close(1114)
       !End If
       ! =========================
    Else  !KEY=1 CALCULATIONS
       ! =========================
       !
       ! Calculations of Invers_f(qq),_f'(qq),_f''(qq),_f'''(qq)
       If(qq.Le.rmm3) Then
          Call cseval(iqqmax,qq,fdsx,fdsy ,fspb0,fspc0,fspd0,fj)
          Call cseval(iqqmax,qq,fdsx,fdsy1,fspb1,fspc1,fspd1,fj1)
          Call cseval(iqqmax,qq,fdsx,fdsy2,fspb2,fspc2,fspd2,fj2)
          Call cseval(iqqmax,qq,fdsx,fdsy3,fspb3,fspc3,fspd3,fj3)
       Else
          fj  = amm3+bmm3*(qq-rmm3)+cmm3*(qq-rmm3)**2/two
          fj1 = bmm3+cmm3*(qq-rmm3); fj2 = cmm3; fj3 = zero
       End If
       ! Calculations of f(fj),f'(fj),f''(fj),f'''(fj)
       f  = qq; f1 = one/fj1; f2 =-fj2*f1**3
       f3 = three*fj2**2*f1**5-fj3*f1**4
    End If
    Return
  End Subroutine thofun
  !=======================================================================
  !> Calculates Legendre decomposition of neutron(proton) 'it=1(2)'
  !> HFB+HO_{L=0}(r) density 'yr' at point 'xr' (in fm)
  !=======================================================================
  Subroutine densitr(it,xr,yr,yrP,msw)
    Use HFBTHO
    Implicit None
    Integer(ipr) :: it,msw
    Integer(ipr), Save :: iw,ik,il,i0,i02,jk,nsa,nsb,nrb,ny,nyy,ib,nd,&
                          n1,n2,n1n2nd,ibit,ibitnb,nzb,mlb,ngh1,ngl1
    ! msw=20 test for protons at the crazy case of U 212 120 92
    ! msw=3  s,a1= 107.665062 80.82396  s/s1= 1.33209
    ! msw=6  s,s1=  88.973061 80.99903  s/s1= 1.09844
    ! msw=12 s,s1=  92.025592 80.99595  s/s1= 1.13617
    ! msw=18 s,s1=  92.000017 80.99595  s/s1= 1.13585
    ! msw=24 s,s1=  92.000010 80.99595  s/s1= 1.13585
    ! Taken up to msw=24
    Real(pr), Allocatable, Save :: xmw(:),yi(:,:)
    Real(pr) :: phy(msw,nzrlx),anl(msw,msw),yl(msw,msw)
    Real(pr), Save :: sl,w,hw,ct2,s,frit,fritP,wdcorin,bzi,bri,ct,st,z,t
    Real(pr) :: xr,yr,yrp
    !
    ngh1=ngh+1; ngl1=ngl+1
    If(it.Lt.0) Then
       it = -it
       If(Allocated(xmw)) Deallocate(xmw,yi)
       Allocate(xmw(msw),yi(msw,msw))
       wdcorin=one/Sqrt(pi*bz*bpp); bzi=one/bz; bri=one/bp
       ! 'msw' mesh-points in 'angle' space
       xmw(1)=zero; hw=half*pi/Real(msw-1,Kind=pr)
       Do il=2,msw
          xmw(il)=hw*Real(il-1,Kind=pr)
       End Do
       ! coefficients for the L-decomposition
       sl=4.0d0
       Do il=1,msw
          sl=0.250d0*sl
          Do ny=1,il
             anl(ny,il) = iv(il-ny)*fak(2*(il+ny-2))*&
                  fi(il-ny)*fi(il+ny-2)*fi(2*ny-2)*sl
          End Do
       End Do
       Do iw=1,msw
          w = xmw(iw); ct2 = Cos(w)**2
          Do il=1,msw
             yi(iw,il) = zero; s = zero
             Do nyy=1,il
                ny = il + 1 - nyy; s = s*ct2 + anl(ny,il)
             End Do
             yl(iw,il) = s*sq(4*il-3)
          End Do
          yi(iw,iw) = one
       End Do
       Call lingd(msw,msw,msw,msw,yl,yi,s,il)
    End If
    ! 'xr/yr' calculations
    Do iw=1,msw
       w = xmw(iw); ct = Cos(w); st = Sin(w)
       z = ct*xr*bzi; t = (st*xr*bri)**2
       Call gaupolr(z,t)
       ik = 0
       Do ib = 1,nb
          nd= id(ib); i0= ia(ib)
          Do n2 = 1,nd
             ik = ik + 1; i02= i0 + n2
             nzb= nz(i02);  nrb= nr(i02); mlb = nl(i02)
             phy(iw,ik) = qh(nzb,ngh1)*ql(nrb,mlb,ngl1)*wdcorin
          End Do
       End Do
    End Do
    ! 'yr' over the blocks
    yr=zero; yrP=zero; ik = 0
    Do ib = 1,nb
       nd= id(ib); i0= ia(ib)
       Do n2 = 1,nd
          jk= ik; ik= ik + 1; i02=i0 + n2; nsb= ns(i02)
          Do n1 = n2,nd
             jk= jk + 1; i02=i0 + n1; nsa= ns(i02)
             If(nsa.Eq.nsb) Then
                ibit=ib; ibitnb= ib+nbx; n1n2nd= n1+(n2-1)*nd
                frit = rk(n1n2nd,ibit); fritP = rk(n1n2nd,ibitnb)
                If(n1.Ne.n2) then
                   frit = two*frit; fritP = two*fritP
                End If
                s = zero
                Do iw=1,msw
                   s = s + yi(1,iw)*phy(iw,ik)*phy(iw,jk)
                End Do
                yr=yr+frit*s; yrP=yrP+fritP*s
             End If
          End Do !n2
       End Do !n1
    End Do !ib
    !
    Return
  End Subroutine densitr
  !=======================================================================
  !
  !=======================================================================
  Subroutine gaupolr(z,x)
    Use HFBTHO
    Implicit None
    Real(pr) :: z,x
    Real(pr) :: w0,w00,w4pii,dsq,d1,d2
    Integer(ipr) :: N,L,NGH1,NGL1
    !
    NGH1=NGH+1; NGL1=NGL+1
    W4PII = PI**(-0.250D0); W0 = W4PII*Exp(-HALF*Z*Z)
    ! W0 = W0*SQRT(Z) NOT MULTIPLIED BY WDCOR
    QH(0,NGH1) = W0; QH(1,NGH1)= SQ(2)*W0*Z
    Do N = 2,NZM
       QH(N,NGH1) = SQI(N)*(SQ(2)*Z*QH(N-1,NGH1)-SQ(N-1)*QH(N-2,NGH1))
    End Do
    W00 = SQ(2)*Exp(-HALF*X)
    Do L = 0,NLM
       If(L.Eq.0) Then
          W0 = W00*Sqrt(HALF)
       Else
          W0 = W00*Sqrt(HALF*X**L)
       End If
       QL(0,L,NGL1) = WFI(L)*W0; QL(1,L,NGL1) = (Real(L+1,Kind=pr)-X)*WFI(L+1)*W0
       Do N = 2,NRM
          DSQ = SQ(N)*SQ(N+L); D1= Real(2*N + L - 1,Kind=pr) - X
          D2  = SQ(N-1)*SQ(N-1+L)
          QL(N,L,NGL1) = (D1*QL(N-1,L,NGL1)-D2*QL(N-2,L,NGL1))/DSQ
       End Do
    End Do
    Return
  End Subroutine gaupolr
  !=======================================================================
  !
  !=======================================================================
End Module HFBTHO_THO
