!------------------------------------------------------------------------------
! constants.f90
!
! Physics and compuational constants for pnFAM.
! 
! M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
!------------------------------------------------------------------------------
module constants
   implicit none
   public
   
   !----------------------------------------------------------------------------
   ! Constants in double precision
   !----------------------------------------------------------------------------
   integer, parameter, private :: dp = kind(1.0d0)
   integer, parameter, private :: dc = kind(cmplx(1.0,1.0,kind=dp))
   
   !----------------------------------------------------------------------------
   ! PNFAM release version
   !----------------------------------------------------------------------------
   character(len=*), parameter :: version ='pnFAM version 1.0 [8/13/2015]'
   
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
   real(dp), parameter :: gA        = 1.0_dp          ! default abs. value of axial vector coupling
   
   ! Approximations
   real(dp), parameter :: r0 = 1.2_dp   ! nuclear radius [fm]
   real(dp), parameter :: Mn = 939.0_dp ! nucleon mass [MeV]
   
   !----------------------------------------------------------------------------
   ! Symbols of the chemical elements.  Useful for more human-readable output.
   !----------------------------------------------------------------------------
   character(len=3), parameter :: element_symbols(118) = &
      ['H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ','Na ','Mg ', &
       'Al ','Si ','P  ','S  ','Cl ','Ar ','K  ','Ca ','Sc ','Ti ','V  ','Cr ', &
       'Mn ','Fe ','Co ','Ni ','Cu ','Zn ','Ga ','Ge ','As ','Se ','Br ','Kr ', &
       'Rb ','Sr ','Y  ','Zr ','Nb ','Mo ','Tc ','Ru ','Rh ','Pd ','Ag ','Cd ', &
       'In ','Sn ','Sb ','Te ','I  ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ', &
       'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ','Lu ','Hf ', &
       'Ta ','W  ','Re ','Os ','Ir ','Pt ','Au ','Hg ','Tl ','Pb ','Bi ','Po ', &
       'At ','Rn ','Fr ','Ra ','Ac ','Th ','Pa ','U  ','Np ','Pu ','Am ','Cm ', &
       'Bk ','Cf ','Es ','Fm ','Md ','No ','Lr ','Rf ','Db ','Sg ','Bh ','Hs ', &
       'Mt ','Ds ','Rg ','Cn ','Uut','Fl ','Uup','Lv ','Uus','Uuo']
   
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
            stop "ERROR in f.hfbtho_parity: p is not +/- 1"
      end select
      
   end function hfbtho_parity
   
end module constants
