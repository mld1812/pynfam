!***********************************************************************
!
!    Copyright (c) 2016, Lawrence Livermore National Security, LLC.
!                        Produced at the Lawrence Livermore National
!                        Laboratory.
!                        Written by Nicolas Schunck, schunck1@llnl.gov
!
!    LLNL-CODE-728299 All rights reserved.
!    LLNL-CODE-573953 All rights reserved.
!
!    Copyright 2017, R. Navarro Perez, N. Schunck, R. Lasseri, C. Zhang,
!                    J. Sarich
!    Copyright 2012, M.V. Stoitsov, N. Schunck, M. Kortelainen, H.A. Nam,
!                    N. Michel, J. Sarich, S. Wild
!    Copyright 2005, M.V. Stoitsov, J. Dobaczewski, W. Nazarewicz, P.Ring
!
!    This file is part of HFBTHO.
!
!    HFBTHO is free software: you can redistribute it and/or modify it
!    under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    HFBTHO is distributed in the hope that it will be useful, but
!    WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with HFBTHO. If not, see <http://www.gnu.org/licenses/>.
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
!                    MATHEMATICAL ROUTINES PACKAGE                     !
!                                                                      !
! ==================================================================== !

!-------------------------------------------------------------------
!> This module contains several routines implementing basic linear algebra
!> and mathematical operations
!----------------------------------------------------------------------
!  Subroutines: - lingd(ma,mx,n,m,a,x,d,Ifl)
!               - csplin(n, x, y, b, c, d)
!               - cseval(n,u,x,y,b,c,d,splf0)
!               - deri(h,n,f1,dunl)
!               - gfv()
!               - sdiag(nmax,n,a,d,x,e,is)
!----------------------------------------------------------------------
Module math

  Use HFBTHO_utilities

  Implicit None

Contains
  !=======================================================================
  !> Solves the system of linear equations A*X = B
  !> At the beginning the matrix B is stored in X. During the calculation
  !> it will be overwritten. D is the determinant of A
  !=======================================================================
  Subroutine lingd(ma,mx,n,m,a,x,d,Ifl)
    Integer(ipr) :: ma,mx,n,m,Ifl
    Integer(ipr), Save :: i,j,k,l,k1,n1
    Real(pr) ::  a(ma,m),x(mx,m),d
    Real(pr), Save :: tollim,one,zero,p,q,tol,cp,cq
    Data tollim/1.d-10/,one/1.d0/,zero/0.d0/
    Ifl = 1; p = zero
    Do i=1,n
       q = zero
       Do j=1,n
          q = q + Abs(a(i,j))
       End Do
       If(q.Gt.p) p = q
    End Do
    tol = tollim*p; d   = one
    Do k=1,n
       p = zero
       Do j=k,n
          q = Abs(a(j,k))
          If(q.Lt.p) Cycle
          p = q; i = j
       End Do
       If (p.Lt.tol) Then
          Write (6,200) ('-',j=1,80),tol,i,k,a(i,k),('-',j=1,80)
200     Format (/1x,80a1/' *****  ERROR IN LINGD , TOLERANZ =',e11.4, &
               ' VALUE OF A(',i3,',',i3,') IS ',e11.4/1x,80a1)
          Ifl = -1
          Return
       End If
       cp = one/a(i,k)
       If(i.Ne.k) Then
          d = -d
          Do l=1,m
             cq = x(i,l); x(i,l) = x(k,l); x(k,l) = cq
          End Do
          Do l=k,n
             cq = a(i,l); a(i,l) = a(k,l); a(k,l) = cq
          End Do
       End If
       d = d*a(k,k)
       If(k.Eq.n) Exit
       k1 = k + 1
       Do i=k1,n
          cq=a(i,k)*cp
          Do l=1,m
             x(i,l)=x(i,l)-cq*x(k,l)
          End Do
          Do l=k1,n
             a(i,l)=a(i,l)-cq*a(k,l)
          End Do
       End Do
    End Do
    Do l=1,m
       x(n,l)=x(n,l)*cp
    End Do
    If(n.Eq.1) Return
    n1=n-1
    Do k=1,n1
       cp = one/a(n-k,n-k)
       Do l=1,m
          cq = x(n-k,l)
          Do i=1,k
             cq = cq-a(n-k,n+1-i)*x(n+1-i,l)
          End Do
          x(n-k,l) = cq*cp
       End Do
    End Do
    Return
  End Subroutine lingd
  !=======================================================================
  !> This routine calculates inverse and determinant of a square matrix
  !> 
  !=======================================================================
  Subroutine calculate_inverse(ndim,mat,det,info)
      Implicit None
      Integer(ipr) ndim,info,info1,info2,n
      Real(pr) sgn
      Complex(pr) det
      Integer(ipr), Allocatable :: ipiv(:)
      Complex(pr), Allocatable :: mat(:,:),work(:)
      !
      Allocate(ipiv(ndim),work(ndim))
      ipiv=0; work=Cmplx(0.,0.)
      !
      ! ZGETRF computes an LU factorization of a complex M-by-N matrix A
      ! using partial pivoting with row interchanges.
      Call Zgetrf(ndim,ndim,mat,ndim,ipiv,info1)
      !
      ! Calculating determinant
      !
      det=Cmplx(1.d0,0.d0); sgn=1.d0
      Do n=1,ndim
         det=det*mat(n,n)
         If(ipiv(n).Ne.n) sgn=-sgn
      End Do
      det=sgn*det
      !
      ! DGETRI computes the inverse of a matrix using the LU factorization
      ! computed by DGETRF.
      Call Zgetri(ndim,mat,ndim,ipiv,work,ndim,info2)
      !
      info=info1+info2
      !
      Deallocate(ipiv,work)
      !
  End Subroutine calculate_inverse
  !=======================================================================
  !> The coefficients \f$ b_i, c_i, d_i, i=1,2,\dots,n \f$ are computed
  !> for a cubic interpolating spline
  !>   \f[
  !>      s(x) = y_i + b_i (x-x_i) + c_i (x-x_i)^2 + d_i (x-x_i)^3
  !>   \f]
  !> for \f$ x_i \leq x \leq x_{i+1} \f$
  !> Input
  !>   - n = the number of data points or knots (n.ge.2)
  !>   - x = the abscissas of the knots in strictly increasing order
  !>   - y = the ordinates of the knots
  !> Output
  !>   - b, c, d  = arrays of spline coefficients as defined above.
  !>   \f[
  !>       y_i = s(x_i), \quad\quad
  !>       b_i = \frac{ds}{dx}(x_i),  \quad\quad
  !>       c_i = \frac{1}{2}\frac{d^{2}s}{dx^{2}}(x_i),  \quad\quad
  !>       d_i = \frac{1}{6}\frac{d^{3}s}{dx^{3}}(x_i)
  !>   \f]
  !> The accompanying function subprogram \ref cseval can be used to
  !> evaluate the spline, its derivative or even its 2nd derivative.
  !=======================================================================
  Subroutine csplin(n, x, y, b, c, d)
    Integer(ipr), Save :: nm1,i,ib
    Integer(ipr) :: n
    Real(pr) :: x(n), y(n), b(n), c(n), d(n)
    Real(pr), Save :: t,zero=0.0d0,two=2.0d0,tr=3.0d0
    ! check input for consistency
    If(n.Lt.2) Stop '-n < 2 in csplin call--'
    nm1 = n-1
    Do i = 1, nm1
       If(x(i).Ge.x(i+1)) Stop 'x not strictly ascending in csplin call'
    End Do
    If (n.Ne.2) Then
       ! set up tridiagonal system
       ! b = diagonal, d = offdiagonal, c = right hand side.
       d(1) = x(2) - x(1); c(2) = (y(2) - y(1))/d(1)
       Do i = 2, nm1
          d(i) = x(i+1) - x(i); b(i) = two*(d(i-1) + d(i))
          c(i+1) = (y(i+1) - y(i))/d(i); c(i) = c(i+1) - c(i)
       End Do
       ! end conditions.  third derivatives at  x(1)  and  x(n)
       ! obtained from divided dIfferences
       b(1) = -d(1); b(n) = -d(n-1); c(1) = zero; c(n) = zero
       If (n.Ne.3) Then
          c(1) =  c(3)/(x(4)-x(2))-c(2)/(x(3)-x(1))
          c(n) =  c(n-1)/(x(n)-x(n-2))-c(n-2)/(x(n-1)-x(n-3))
          c(1) =  c(1)*d(1)**2/(x(4)-x(1))
          c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
          ! forward elimination
       Else
          Do i = 2, n
             t = d(i-1)/b(i-1); b(i) = b(i) - t*d(i-1); c(i) = c(i) - t*c(i-1)
          End Do
       End If
       ! back substitution
       c(n) = c(n)/b(n)
       Do ib = 1, nm1
          i = n-ib
          c(i) = (c(i) - d(i)*c(i+1))/b(i)
       End Do
       ! compute polynomial coefficients
       b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + two*c(n))
       Do i = 1, nm1
          b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + two*c(i))
          d(i) = (c(i+1) - c(i))/d(i); c(i) = tr*c(i)
       End Do
       c(n) = tr*c(n); d(n) = d(n-1)
       Return
    Else
       b(1) = (y(2)-y(1))/(x(2)-x(1)); c(1) = zero; d(1) = zero
       Return
    End If
  End Subroutine csplin
  !=======================================================================
  !
  !=======================================================================
  Subroutine cseval(n,u,x,y,b,c,d,splf0)
    Integer(ipr) :: n
    Integer(ipr), Save :: i=1,j,k
    Real(pr) :: x(n),y(n),b(n),c(n),d(n),u,splf0
    Real(pr), Save :: dx
    If(i.Ge.n)      i = 1
    If(u.Lt.x(i))   Go To 10
    If(u.Le.x(i+1)) Go To 30
    ! binary search
10  i = 1
    j = n + 1
20  k = (i+j)/2
    If(u.Lt.x(k)) j = k
    If(u.Ge.x(k)) i = k
    If(j.Gt.i+1) Go To 20
    ! evaluate splf0
30 dx = u - x(i)
    splf0 = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
    Return
  End Subroutine cseval
  !=======================================================================
  !> First derivative of 'f1' if the step is 'h'
  !=======================================================================
  Subroutine deri(h,n,f1,dunl)
    Integer(ipr) :: n
    Integer(ipr), Save :: k
    Real(pr) :: h,f1(n),dunl(n)
    Real(pr), Save :: t60,t12
    Real(pr), Save :: t8=8.0d0,t45=45.0d0,t9=9.0d0
    t60 =1.0d0/(h*60.0d0); t12 =1.0d0/(h*12.0d0)
    !
    dunl(1)  =(t8*f1(2)-f1(3)+f1(1))*t12
    dunl(2)  =(t45*(f1(3)-f1(1))-t9*f1(4)+f1(5)-f1(1))*t60
    dunl(3)  =(t45*(f1(4)-f1(2))-t9*(f1(5)-f1(1))+f1(6))*t60
    dunl(n)  =(-t8*f1(n-1)+f1(n)+f1(n-2))*t12
    dunl(n-1)=(t45*(f1(n)-f1(n-2))+t9*f1(n-3)-f1(n)-f1(n-4))*t60
    dunl(n-2)=(t45*(f1(n-1)-f1(n-3))-t9*(f1(n)-f1(n-4))-f1(n-5))*t60
    Do k=4,n-3
       dunl(k) =(t45*(f1(k+1)-f1(k-1))-t9*(f1(k+2)-f1(k-2))+f1(k+3)-f1(k-3))*t60
    End Do
    Return
  End Subroutine deri
  !=======================================================================
  !> Calculates sign, Sqrt, factorials, etc. of integers and half integers
  !>   - \f$ \mathtt{iv(n) } = (-1)^n    \f$
  !>   - \f$ \mathtt{sq(n) } = \sqrt{n}, \f$
  !>   - \f$ \mathtt{sqi(n)} = 1/\sqrt{n}\f$
  !>   - \f$ \mathtt{fak(n)} = n!        \f$
  !>   - \f$ \mathtt{wf(n) } = \sqrt{n!} \f$
  !>   - \f$ \mathtt{wfi(n)} = 1/\sqrt{n!}\f$
  !=======================================================================
  Subroutine gfv
    Use HFBTHO
    Implicit None
    Integer(ipr) :: i,igfv
    Parameter(igfv=170)               !maximal number for GFV
    If(Allocated(iv)) Deallocate(iv,fak,fi,sq,sqi,wf,wfi)
    Allocate(iv(-igfv:igfv),fak(0:igfv),fi(0:igfv),sq(0:igfv),sqi(0:igfv))
    Allocate(wf(0:igfv),wfi(0:igfv))
    iv(0)=1; sq(0)=zero; sqi(0)=1.0d30
    fak(0)=one; fi(0)=one; wf(0)=one; wfi(0)=one
    Do i=1,igfv
       iv(i)=-iv(i-1)
       iv(-i) = iv(i)
       sq(i)=Sqrt(Real(i,Kind=pr)); sqi(i)=one/sq(i)
       fak(i)=Real(i,Kind=pr)*fak(i-1); fi(i)=one/fak(i)
       wf(i)=sq(i)*wf(i-1); wfi(i)=one/wf(i)
    End Do
  End Subroutine gfv
  !=======================================================================
  !> Diagonalization of a real, symemtric matrix (backup routine if LAPACK fails)
  !=======================================================================
  Subroutine sdiag(nmax,n,a,d,x,e,is)
    !---------------------------------------------------------------------
    ! A   matrix to be diagonalized
    ! D   eigenvalues,  X   eigenvectors, E   auxiliary field
    ! IS=1  eigenvalues are ordered (major component of X is positive)
    ! 0  eigenvalues are not ordered
    !---------------------------------------------------------------------
    Use HFBTHO_utilities, Only: pr,ipr
    Implicit None
    Integer(ipr), Save :: i,j,j1,k,l,im
    Integer(ipr)       :: n,nmax,is
    Real(pr), Save :: f,g,h,hi,s,p,b,r,pra,c
    Real(pr) :: a(nmax,nmax),x(nmax,nmax),e(n),d(n)
    Real(pr), Save :: tol=1.0D-32,eps=9.0D-12,one=1.0_pr,zero=0.0_pr
    !
    If (n.Le.1) Then
       d(1)=a(1,1); x(1,1)=one
       Return
    End If
    Do i=1,n
       Do j=1,i
          x(i,j)=a(i,j)
       End Do
    End Do
    ! householder-reduktion
    i=n
15 Continue
    If (i.Ge.2) Then
       l=i-2
       f=x(i,i-1); g=f; h=zero
       If (l.Gt.0) Then
          Do k=1,l
             h=h+x(i,k)*x(i,k)
          End Do
       End If
       s=h+f*f
       If (s.Lt.tol) Then
          h=zero
          Go To 100
       End If
       If (h.Gt.zero) Then
          l=l+1; g=Sqrt(s)
          If (f.Ge.zero) g=-g
          h=s-f*g; hi=one/h; x(i,i-1)=f-g; f=zero
          If (l.Gt.0) Then
             Do j=1,l
                x(j,i)=x(i,j)*hi
                s=zero
                Do k=1,j
                   s=s+x(j,k)*x(i,k)
                End Do
                j1=j+1
                If (l.Ge.j1) Then
                   Do k=j1,l
                      s=s+x(k,j)*x(i,k)
                   End Do
                End If
                e(j)=s*hi; f=f+s*x(j,i)
             End Do
          End If
          f=f*hi*0.50_pr
          If (l.Gt.0) Then
             Do j=1,l
                s=x(i,j); e(j)=e(j)-f*s; p=e(j)
                Do  k=1,j
                   x(j,k)=x(j,k)-s*e(k)-x(i,k)*p
                End Do
             End Do
          End If
       End If
100  Continue
       d(i)=h; e(i-1)=g; i=i-1
       Go To 15
       ! Bereitstellen der Transformationmatrix
    End If
    d(1)=zero; e(n)=zero; b=zero; f=zero
    Do i=1,n
       l=i-1
       If (d(i).Eq.0.) Go To 221
       If (l.Gt.0) Then
          Do J=1,L
             s=zero
             Do k=1,l
                s=s+x(i,k)*x(k,j)
             End Do
             Do k=1,l
                x(k,j)=x(k,j)-s*x(k,i)
             End Do
          End Do
       End If
221  Continue
       d(i)=x(i,i)
       x(i,i)=one
       If (l.Gt.0) Then
          Do j=1,l
             x(i,j)=zero; x(j,i)=zero
          End Do
       End If
    End Do
    ! Diagonalisieren der Tri-Diagonal-Matrix
    Do l=1,n
       h=eps*(Abs(d(l))+ Abs(e(l)))
       If (h.Gt.b) b=h
       ! Test fuer Splitting
       Do  j=l,n
          If (Abs(e(j)).Le.b) Exit
       End Do
       ! test fuer konvergenz
       If (j.Eq.l) Go To 300
340  p=(d(l+1)-d(l))/(2.0_pr*e(l))
       r=Sqrt(p*p+one); pra=p+r
       If (p.Lt.zero) pra=p-r
       h=d(l)-e(l)/pra
       Do i=l,n
          d(i)=d(i)-h
       End Do
       f=f+h
       ! QR-transformation
       p=d(j); c=one; s=zero; i=j
360  i=i-1
       If (i.Lt.l) Go To 362
       g=c*e(i); h=c*p
       If ( Abs(p).Ge.Abs(e(i))) Then
          c=e(i)/p
          r=Sqrt(c*c+one); e(i+1)=s*p*r; s=c/r; c=one/r
          Go To 365
       End If
       c=p/e(i)
       r=Sqrt(c*c+one); e(i+1)=s*e(i)*r; s=one/r; c=c/r
365  p=c*d(i)-s*g
       d(i+1)=h+s*(c*g+s*d(i))
       Do k=1,n
          h=x(k,i+1); x(k,i+1)=x(k,i)*s+h*c
          x(k,i)=x(k,i)*c-h*s
       End Do
       Go To 360
362  e(l)=s*p
       d(l)=c*p
       If ( Abs(e(l)).Gt.b) Go To 340
       ! konvergenz
300  d(l)=d(l)+f
    End Do
    If (is.Eq.0) Return
    ! ordnen der eigenwerte
    Do i=1,n
       k=i; p=d(i); j1=i+1
       If (j1.Le.n) Then
          Do j=j1,n
             If (d(j).Ge.p) Cycle
             k=j; p=d(j)
          End Do
          If (k.Eq.i) Cycle
          d(k)=d(i); d(i)=p
          Do j=1,n
             p=x(j,i); x(j,i)=x(j,k)
             x(j,k)=p
          End Do
       End If
    End Do
    ! signum
    Do  k=1,n
       s=zero
       Do i=1,n
          h=Abs(x(i,k))
          If (h.Gt.s) Then
             s=h; im=i
          End If
       End Do
       If (x(im,k).Lt.zero) Then
          Do i=1,n
             x(i,k)=-x(i,k)
          End Do
       End If
    End Do
  End Subroutine sdiag
  !=======================================================================
  !
  !=======================================================================
End Module math
