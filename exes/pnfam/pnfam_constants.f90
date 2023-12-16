!------------------------------------------------------------------------------
!> This module contains physics and compuational constants for pnFAM. It also
!> contains some routines for manipulating strings.
!>
!> @authors M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
!> @authors E.M. Ney, UNC Chapel Hill, 2018-
!------------------------------------------------------------------------------
module pnfam_constants
!$ use omp_lib, only : omp_get_wtime
   use pnfam_logger
   implicit none
   public
   
   !----------------------------------------------------------------------------
   ! Precision of the solver
   !----------------------------------------------------------------------------
   integer, parameter, private :: ipr = kind(1)
   integer, parameter, private :: dp  = kind(1.0d0)
   integer, parameter, private :: dc  = kind(cmplx(1.0,1.0,kind=dp))
   
   !----------------------------------------------------------------------------
   ! OpenMP(?)
   !----------------------------------------------------------------------------
   integer :: nthreads=0
   logical :: using_openmp=.false.
   
   !----------------------------------------------------------------------------
   ! PNFAM release version
   !----------------------------------------------------------------------------
   include "pnfam_version.inc"
   
   !----------------------------------------------------------------------------
   ! Numerical constants 
   !----------------------------------------------------------------------------
   complex(dc), parameter :: iu  = cmplx(0.0_dp,1.0_dp,kind=dp)
   real(dp),    parameter :: pi  = 3.141592653589793238462643_dp
   real(dp),    parameter :: ln2 = log(2.0_dp)
   
   !----------------------------------------------------------------------------
   ! Physical constants 
   !----------------------------------------------------------------------------
   ! From NIST CODATA, 2013
   real(dp), parameter :: alpha    = 7.2973525698d-03 ! fine-structure const.
   real(dp), parameter :: hbarc    = 197.3269718_dp   ! hbar c [MeV fm]
   real(dp), parameter :: hbar_mec = 386.15926800_dp  ! (hbar c)/(me c^2) [fm]
   real(dp), parameter :: mec2     = 0.510998928_dp   ! me c^2 [MeV]
   real(dp), parameter :: lambda   = 687_dp           ! breakdown scale of the chiral EFT (MeV)
   
   !! Physics constants used in the Landolt-Bornstein tables
   !! These should be used for any comparisons
   !real(dp), parameter :: alpha    = 1.0_dp/137.0388_dp
   !real(dp), parameter :: hbar_mec = r0/0.42587_dp/alpha
   !real(dp), parameter :: mec2     = 0.511006_dp
   
   ! Beta decay
   real(dp), parameter :: kappa     = 6147.0_dp       ! decay const.
   real(dp), parameter :: ln2_kappa = ln2/kappa       ! multiplying factor for rates   
   real(dp), parameter :: dmassnH   = 0.78227_dp      ! n-H mass difference [MeV]
   real(dp), parameter :: gV        = 1.0_dp          ! default vector coupling constant
   real(dp), parameter :: gA_eff    = 1.0_dp          ! effective abs. value of axial vector coupling
   real(dp), parameter :: gA        = 1.27_dp         ! true abs. value of axial vector coupling
   real(dp), parameter :: Fpi       = 92.4_dp         ! pion decay constant [MeV]
   real(dp), parameter :: Mpi       = 138.04_dp       ! pion mass [MeV]
   
   ! Approximations
   real(dp), parameter :: r0 = 1.2_dp   ! nuclear radius [fm]
   real(dp), parameter :: Mn = 939.0_dp ! nucleon mass [MeV]
   
   !----------------------------------------------------------------------------
   ! Symbols of the chemical elements.  Useful for more human-readable output.
   !----------------------------------------------------------------------------
   character(len=3), parameter :: element_symbols(0:118) = &
      ['?? ',&
       'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ','Na ','Mg ', &
       'Al ','Si ','P  ','S  ','Cl ','Ar ','K  ','Ca ','Sc ','Ti ','V  ','Cr ', &
       'Mn ','Fe ','Co ','Ni ','Cu ','Zn ','Ga ','Ge ','As ','Se ','Br ','Kr ', &
       'Rb ','Sr ','Y  ','Zr ','Nb ','Mo ','Tc ','Ru ','Rh ','Pd ','Ag ','Cd ', &
       'In ','Sn ','Sb ','Te ','I  ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ', &
       'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ','Lu ','Hf ', &
       'Ta ','W  ','Re ','Os ','Ir ','Pt ','Au ','Hg ','Tl ','Pb ','Bi ','Po ', &
       'At ','Rn ','Fr ','Ra ','Ac ','Th ','Pa ','U  ','Np ','Pu ','Am ','Cm ', &
       'Bk ','Cf ','Es ','Fm ','Md ','No ','Lr ','Rf ','Db ','Sg ','Bh ','Hs ', &
       'Mt ','Ds ','Rg ','Cn ','Nh ','Fl ','Mc ','Lv ','Ts ','Og ']
   
   !------------------------------------------------------------------------
   ! Symbols ported over from HFBTHO 
   !------------------------------------------------------------------------
   character(len=*), dimension(2), parameter :: strpar = ['+','-']
   
   !----------------------------------------------------------------------------
   ! Constants to avoid forgetting what different integers mean
   !----------------------------------------------------------------------------
   integer, parameter :: IT_NEUTRON = 1, IT_PROTON = 2, IT_ISOSCALAR = 3
   
   
contains
   
   !----------------------------------------------------------------------------
   ! Converts a string to uppercase (is there really no Fortran standard
   ! procedure for this?)
   !----------------------------------------------------------------------------
   subroutine translate_uppercase(s)
      implicit none
      character(len=*), intent(inout) :: s
      
      integer :: i, cc
      integer, parameter :: ascii_a = iachar('a'), ascii_z = iachar('z'), &
         ascii_diff = iachar('A') - iachar('a')
      
      if (len(s) == 0) return
      do i=1,len(s)
         cc = iachar(s(i:i))
         if (cc >= ascii_a .and. cc <= ascii_z) s(i:i) = achar(cc + ascii_diff)
      end do
      
   end subroutine translate_uppercase

   !----------------------------------------------------------------------------
   ! Auxiliary formatting function: Join a list of strings with comma, trimming each element
   !----------------------------------------------------------------------------
   function txt_clist(labels)
      implicit none
      character(len=*), dimension(:), intent(in) :: labels
      character(len=100) :: txt_clist
      integer :: i
      txt_clist=""
      if (size(labels) > 0) then
         txt_clist = trim(labels(1))
         do i=2, size(labels)
            txt_clist = trim(txt_clist)//', '//trim(labels(i))
         end do
      end if
   end function txt_clist

   !---------------------------------------------------------------------------
   ! Auxiliary formatting function: Pad with spaces on the left
   !---------------------------------------------------------------------------
   function txtr(str, n)
      implicit none
      integer,          intent(in) :: n
      character(len=*), intent(in) :: str
      character(len=n) :: txtr
      txtr = repeat(" ",n-len_trim(str))//trim(str)
   end function txtr

   !---------------------------------------------------------------------------
   ! Auxiliary formatting function: Pad with spaces on the right
   !---------------------------------------------------------------------------
   function txtl(str, n)
      implicit none
      integer,          intent(in) :: n
      character(len=*), intent(in) :: str
      character(len=n) :: txtl
      txtl = trim(str)//repeat(" ",n-len_trim(str))
   end function txtl

   !----------------------------------------------------------------------------
   ! Auxialiary formatting function: Center a string in a field of a fixed width
   !----------------------------------------------------------------------------
   function txtc(text, width)
      implicit none
      character(len=*), intent(in) :: text
      integer,          intent(in) :: width
      character(len=max(width,len_trim(text))) :: txtc
      integer :: lpad

      ! Extra padding on left for odd number
      lpad = (width-len_trim(text))/2 + mod(width-len_trim(text),2)
      if(lpad<0) then
         txtc = text
      else
         txtc = repeat(' ', lpad)//adjustl(text)
      end if
   end function txtc
   
   !----------------------------------------------------------------------------
   ! Convert a "usual" parity integer (+/- 1) to an HFBTHO-style integer
   !----------------------------------------------------------------------------
   function hfbtho_parity(p) result(hfbtho_p)
      implicit none
      integer, intent(in) :: p
      integer :: hfbtho_p
      
      select case (p)
         case ( 1)
            hfbtho_p = 1
         case (-1)
            hfbtho_p = 2
         case default
            call abort(" ERROR in f.hfbtho_parity: p is not +/- 1")
      end select
      
   end function hfbtho_parity

   !----------------------------------------------
   ! Get the wallclock time with or without OpenMP, with high precision
   !----------------------------------------------
   function get_timer() result(t)
      implicit none
      integer(8)  :: it, ir
      real(dp) :: t
      if (using_openmp .eqv. .true.) then
         !$ t = omp_get_wtime()
      else
         call system_clock(count=it, count_rate=ir)
         t = 1.0_dp*it/ir
      end if
   end function get_timer

   !----------------------------------------------
   ! Extract a digit from an integer
   !----------------------------------------------
   function get_digit(n, p) result(d)
      implicit none
      integer, intent(in) :: n, p
      integer :: d

      d = n/10**p - 10*(n/10**(p+1)) ! integer division

   end function get_digit

end module pnfam_constants
