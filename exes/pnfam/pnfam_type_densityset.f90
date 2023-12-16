!------------------------------------------------------------------------------
!> This module contains the derived type for the set of coordinate-space
!> densities along with its constructor.
!>
!> @authors M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
!> @authors E.M. Ney, UNC Chapel Hill, 2018-
!------------------------------------------------------------------------------
module type_densityset

   implicit none

   integer, parameter, private :: dp = kind(1.0d0)

   type density_set
      ! Time-even densities
      real(dp), dimension(:), allocatable :: rerho, imrho
      real(dp), dimension(:), allocatable :: retau, imtau
      real(dp), dimension(:), allocatable :: retjrr, imtjrr, retjrp, imtjrp, retjrz, imtjrz
      real(dp), dimension(:), allocatable :: retjpr, imtjpr, retjpp, imtjpp, retjpz, imtjpz
      real(dp), dimension(:), allocatable :: retjzr, imtjzr, retjzp, imtjzp, retjzz, imtjzz

      ! Time-odd densities
      real(dp), dimension(:), allocatable :: resr, imsr, resp, imsp, resz, imsz
      real(dp), dimension(:), allocatable :: retr, imtr, retp, imtp, retz, imtz
      real(dp), dimension(:), allocatable :: rejr, imjr, rejp, imjp, rejz, imjz
      real(dp), dimension(:), allocatable :: refr, imfr, refp, imfp, refz, imfz
      real(dp), dimension(:), allocatable :: regs, imgs

      ! Pairing densities
      real(dp), dimension(:), allocatable :: rerb, imrb  ! rho-breve_{10}
      real(dp), dimension(:), allocatable :: resbr, imsbr, resbp, imsbp, resbz, imsbz  ! s-breve

      logical :: alloc=.false. ! every new instance is not allocated

   end type density_set

contains

   subroutine allocate_density(ds, nghl)
      implicit none
      type(density_set), intent(inout) :: ds
      integer, intent(in) :: nghl

      if (allocated(ds%rerho)) then
         deallocate(ds%rerho, ds%imrho)
         deallocate(ds%retau, ds%imtau)
         deallocate(ds%retjrr, ds%imtjrr, ds%retjrp, ds%imtjrp, ds%retjrz, ds%imtjrz)
         deallocate(ds%retjpr, ds%imtjpr, ds%retjpp, ds%imtjpp, ds%retjpz, ds%imtjpz)
         deallocate(ds%retjzr, ds%imtjzr, ds%retjzp, ds%imtjzp, ds%retjzz, ds%imtjzz)
         deallocate(ds%resr, ds%imsr, ds%resp, ds%imsp, ds%resz, ds%imsz)
         deallocate(ds%retr, ds%imtr, ds%retp, ds%imtp, ds%retz, ds%imtz)
         deallocate(ds%rejr, ds%imjr, ds%rejp, ds%imjp, ds%rejz, ds%imjz)
         deallocate(ds%refr, ds%imfr, ds%refp, ds%imfp, ds%refz, ds%imfz)
         deallocate(ds%regs, ds%imgs)
         ! Pairing densities
         deallocate(ds%rerb, ds%imrb)
         deallocate(ds%resbr, ds%imsbr, ds%resbp, ds%imsbp, ds%resbz, ds%imsbz)
      end if

      allocate(ds%rerho(nghl), ds%imrho(nghl))
      allocate(ds%retau(nghl), ds%imtau(nghl))
      allocate(ds%retjrr(nghl), ds%imtjrr(nghl), ds%retjrp(nghl), ds%imtjrp(nghl))
      allocate(ds%retjrz(nghl), ds%imtjrz(nghl), ds%retjpr(nghl), ds%imtjpr(nghl))
      allocate(ds%retjpp(nghl), ds%imtjpp(nghl), ds%retjpz(nghl), ds%imtjpz(nghl))
      allocate(ds%retjzr(nghl), ds%imtjzr(nghl), ds%retjzp(nghl), ds%imtjzp(nghl))
      allocate(ds%retjzz(nghl), ds%imtjzz(nghl))

      allocate(ds%resr(nghl), ds%imsr(nghl), ds%resp(nghl), ds%imsp(nghl))
      allocate(ds%resz(nghl), ds%imsz(nghl))
      allocate(ds%rejr(nghl), ds%imjr(nghl), ds%rejp(nghl), ds%imjp(nghl))
      allocate(ds%rejz(nghl), ds%imjz(nghl))
      allocate(ds%retr(nghl), ds%imtr(nghl), ds%retp(nghl), ds%imtp(nghl))
      allocate(ds%retz(nghl), ds%imtz(nghl))
      allocate(ds%refr(nghl), ds%imfr(nghl), ds%refp(nghl), ds%imfp(nghl))
      allocate(ds%refz(nghl), ds%imfz(nghl))
      allocate(ds%regs(nghl), ds%imgs(nghl))

      ! Pairing densities
      allocate(ds%rerb(nghl), ds%imrb(nghl))
      allocate(ds%resbr(nghl), ds%imsbr(nghl))
      allocate(ds%resbp(nghl), ds%imsbp(nghl))
      allocate(ds%resbz(nghl), ds%imsbz(nghl))

      ds%alloc = .true.

   end subroutine allocate_density
   
   subroutine deallocate_density(ds)
      implicit none
      type(density_set), intent(inout) :: ds

      if (allocated(ds%rerho)) then
         deallocate(ds%rerho, ds%imrho)
         deallocate(ds%retau, ds%imtau)
         deallocate(ds%retjrr, ds%imtjrr, ds%retjrp, ds%imtjrp, ds%retjrz, ds%imtjrz)
         deallocate(ds%retjpr, ds%imtjpr, ds%retjpp, ds%imtjpp, ds%retjpz, ds%imtjpz)
         deallocate(ds%retjzr, ds%imtjzr, ds%retjzp, ds%imtjzp, ds%retjzz, ds%imtjzz)
         deallocate(ds%resr, ds%imsr, ds%resp, ds%imsp, ds%resz, ds%imsz)
         deallocate(ds%retr, ds%imtr, ds%retp, ds%imtp, ds%retz, ds%imtz)
         deallocate(ds%rejr, ds%imjr, ds%rejp, ds%imjp, ds%rejz, ds%imjz)
         deallocate(ds%refr, ds%imfr, ds%refp, ds%imfp, ds%refz, ds%imfz)
         deallocate(ds%regs, ds%imgs)
         deallocate(ds%rerb, ds%imrb, ds%resbz, ds%imsbz, ds%resbr, ds%imsbr, ds%resbp, ds%imsbp)
      end if

      ds%alloc = .false.

   end subroutine deallocate_density

   subroutine zero_density(ds)
      implicit none
      type(density_set), intent(inout) :: ds

      ! Scalar density (rho)
      ds%rerho = 0.0_dp;  ds%imrho = 0.0_dp

      ! Kinetic density (tau)
      ds%retau = 0.0_dp;  ds%imtau = 0.0_dp

      ! Tensor density (J)
      ds%retjrr = 0.0_dp;  ds%imtjrr = 0.0_dp;  ds%retjrp = 0.0_dp;  ds%imtjrp = 0.0_dp
      ds%retjrz = 0.0_dp;  ds%imtjrz = 0.0_dp;  ds%retjpr = 0.0_dp;  ds%imtjpr = 0.0_dp
      ds%retjpp = 0.0_dp;  ds%imtjpp = 0.0_dp;  ds%retjpz = 0.0_dp;  ds%imtjpz = 0.0_dp
      ds%retjzr = 0.0_dp;  ds%imtjzr = 0.0_dp;  ds%retjzp = 0.0_dp;  ds%imtjzp = 0.0_dp
      ds%retjzz = 0.0_dp;  ds%imtjzz = 0.0_dp

      ! Spin density (s)
      ds%resr = 0.0_dp;  ds%imsr = 0.0_dp
      ds%resp = 0.0_dp;  ds%imsp = 0.0_dp
      ds%resz = 0.0_dp;  ds%imsz = 0.0_dp

      ! Current density (j)
      ds%rejr = 0.0_dp;  ds%imjr = 0.0_dp
      ds%rejp = 0.0_dp;  ds%imjp = 0.0_dp
      ds%rejz = 0.0_dp;  ds%imjz = 0.0_dp

      ! Spin-kinetic density (T)
      ds%retr = 0.0_dp;  ds%imtr = 0.0_dp
      ds%retp = 0.0_dp;  ds%imtp = 0.0_dp
      ds%retz = 0.0_dp;  ds%imtz = 0.0_dp

      ! Tensor-kinetic density (F)
      ds%refr = 0.0_dp;  ds%imfr = 0.0_dp
      ds%refp = 0.0_dp;  ds%imfp = 0.0_dp
      ds%refz = 0.0_dp;  ds%imfz = 0.0_dp

      ! Divergence of spin density (Del.s)
      ds%regs = 0.0_dp;  ds%imgs = 0.0_dp

      ! Pairing densities
      ds%rerb  = 0.0_dp;  ds%imrb  = 0.0_dp
      ds%resbr = 0.0_dp;  ds%imsbr = 0.0_dp
      ds%resbp = 0.0_dp;  ds%imsbp = 0.0_dp
      ds%resbz = 0.0_dp;  ds%imsbz = 0.0_dp
      
   end subroutine zero_density

end module type_densityset
