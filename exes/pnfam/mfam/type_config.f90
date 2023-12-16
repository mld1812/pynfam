!------------------------------------------------------------------------------
! hfb_solution_type.f90
!
! Store the Matrix-PNFAM configuration namelist.
!
! M.T. Mustonen and T. Shafer, UNC Chapel Hill 2015
!------------------------------------------------------------------------------
module config_type
   implicit none
   integer, parameter, private :: dp = kind(1.0d0)

   type matrix_config
      ! Namelist FILENAMES and related
      character(len=200) :: file_hfb
      character(len=200) :: file_basename
      character(len=200) :: file_ab
      character(len=200) :: file_xy
      character(len=200) :: file_me_one
      character(len=200) :: file_me_two
      character(len=200) :: file_vv
      ! Namelist BASIS_OPTIONS and related
      integer            :: basis_k
      integer            :: basis_parity
      real(dp)           :: basis_cutoff
      ! Namelist DIAGONALIZATION_OPTIONS and related
      real(dp)           :: diag_ev_limit
      integer            :: diag_blocking
      integer            :: diag_nprow
      ! Namelist MATRIX_ELEMENT_OPTIONS and related
      logical            :: me_check_norm
      logical            :: me_check_complete
      logical            :: me_disable_correction
      ! Namelist ODD_NUCLEUS_OPTIONS and related
      integer            :: odd_block_n(2)
      integer            :: odd_block_p(2)
      character(len=200) :: odd_even_basename
      character(len=200) :: odd_even_file_ab
      character(len=200) :: odd_even_file_xy
   end type matrix_config

contains

   !----------------------------------------------------------------------------
   ! Set sensible null parameters for the config
   !----------------------------------------------------------------------------
   subroutine init_config(config)
      implicit none
      type(matrix_config), intent(inout) :: config

      config%file_hfb              = ""
      config%file_basename         = ""
      config%file_ab               = ""
      config%file_xy               = ""
      config%file_me_one           = ""
      config%file_me_two           = ""
      config%file_vv               = ""
      config%basis_k               = 0
      config%basis_parity          = 1
      config%basis_cutoff          = -1.0_dp
      config%diag_ev_limit         = 120.0_dp
      config%diag_blocking         = 8
      config%diag_nprow            = 1
      config%me_check_norm         = .false.
      config%me_check_complete     = .false.
      config%me_disable_correction = .false.
      config%odd_block_n(:)        = [0, 1]
      config%odd_block_p(:)        = [0, 1]
      config%odd_even_basename     = ""
      config%odd_even_file_ab      = ""
      config%odd_even_file_xy      = ""

   end subroutine init_config

   !----------------------------------------------------------------------------
   ! Read the namelist from file
   !----------------------------------------------------------------------------
   subroutine read_namelist_from_file(fh, config)
      implicit none

      integer,             intent(in)  :: fh
      type(matrix_config), intent(out) :: config

      integer :: ierr

      ! Namelist FILENAMES
      character(len=200) :: hfb_filename = "", storage_basename = ""
      namelist /filenames/ hfb_filename, storage_basename

      ! Namelist BASIS_OPTIONS
      integer  :: k = 0, parity = 1
      real(dp) :: ecutoff = -1.0_dp
      namelist /basis_options/ k, parity, ecutoff, block_n, block_p

      ! Namelist DIAGONALIZATION_OPTIONS
      integer  :: blocking_factor = 8, processes_per_dim = 1
      real(dp) :: eigenvalue_upper_limit = 120.0_dp
      namelist /diagonalization_options/ eigenvalue_upper_limit, blocking_factor, processes_per_dim

      ! Namelist MATRIX_ELEMENT_OPTIONS
      logical :: check_xy_orthonormality = .false., check_xy_completeness = .false., disable_polarization_corr = .false.
      character(len=200) :: even_even_xy_filename = ""
      namelist /matrix_element_options/ check_xy_orthonormality, check_xy_completeness, &
                                        even_even_xy_filename, disable_polarization_corr

      ! Namelist ODD_NUCLEUS_OPTIONS
      integer :: block_n(1:2) = [0, 1], block_p(1:2) = [0, 1]
      character(len=200) :: even_basename = ""
      namelist /odd_nucleus_options/ block_n, block_p, even_basename

      call init_config(config)

      read(fh, nml=filenames, iostat=ierr)
      if (ierr /= 0) stop "ERROR reading namelist in config_type.read_namelist_from_file"

      config%file_hfb          = trim(hfb_filename)
      config%file_basename     = trim(storage_basename)
      config%file_ab           = trim(storage_basename)//".i.ab"
      config%file_xy           = trim(storage_basename)//".xy"
      config%file_me_one       = trim(storage_basename)//".1qp"
      config%file_me_two       = trim(storage_basename)//".2qp"
      config%file_vv           = trim(storage_basename)//".vv"

      read(fh, nml=basis_options)
      if (ierr /= 0) stop "ERROR reading namelist in config_type.read_namelist_from_file"

      config%basis_k           = k
      config%basis_parity      = parity
      config%basis_cutoff      = ecutoff

      read(fh, nml=diagonalization_options)
      if (ierr /= 0) stop "ERROR reading namelist in config_type.read_namelist_from_file"

      config%diag_ev_limit     = eigenvalue_upper_limit
      config%diag_blocking     = blocking_factor
      config%diag_nprow        = processes_per_dim

      read(fh, nml=matrix_element_options)
      if (ierr /= 0) stop "ERROR reading namelist in config_type.read_namelist_from_file"

      config%me_check_norm         = check_xy_orthonormality
      config%me_check_complete     = check_xy_completeness
      config%me_disable_correction = disable_polarization_corr

      read(fh, nml=odd_nucleus_options)
      if (ierr /= 0) stop "ERROR reading namelist in config_type.read_namelist_from_file"

      config%odd_block_n       = block_n(:)
      config%odd_block_p       = block_p(:)
      config%odd_even_basename = trim(even_basename)
      config%odd_even_file_ab  = trim(config%odd_even_basename)//".i.ab"
      config%odd_even_file_xy  = trim(config%odd_even_basename)//".xy"

   end subroutine read_namelist_from_file

end module config_type
