!-------------------------------------------------------------------------------
! pnfam_matrix/aux_odd_information.f90
!
! T. Shafer, UNC Chapel Hill, 2015
!-------------------------------------------------------------------------------
program pnfam_matrix_odd_information
   use logger
   use constants,         only : IT_NEUTRON, IT_PROTON, strpar
   use config_type,       only : matrix_config, read_namelist_from_file
   use hfb_solution_type, only : hfb_solution
   use blockmatrix_type,  only : blockmatrix
   use hfbtho_basis,      only : get_hfbtho_solution, get_hfbtho_occupations, lowest_energy_qp, &
                                 nl, ns, npar
   use pnfam_matrix,      only : input_file_from_cli, write_main_header, qp_reverse, &
                                 find_blocked_level
   implicit none

   integer, parameter :: dp = kind(1.0d0)

   character(len=250)  :: infile
   type(matrix_config) :: config
   type(hfb_solution), target :: hfb
   type(blockmatrix), pointer :: vmat1, vmat2
   integer :: bo_k, bo_parity, bo_isospin, bo_index, bo_index_r, bo_k_sign, bo_final_index
   integer :: twoqp_init_index, twoqp_final_index
   real(dp) :: delta_n, delta_z
   real(dp), pointer :: eqp1(:), eqp2(:)
   real(dp), allocatable :: occ_n(:), occ_p(:)

   !----------------------------------------------------------------------------
   !----------------------------------------------------------------------------

   ! Logging
   openlog  => openlog_stdout
   writelog => writelog_stdout
   closelog => closelog_stdout
   call openlog('')

   ! Read input file and print namelist parameters
   call input_file_from_cli(infile, default='matrix.in')
   open(1, file=trim(infile), status='old')
   call read_namelist_from_file(1, config=config)
   close(1)

   call write_main_header(config, infile, 0, 0, 'Odd-mass Nucleus Information')

   call get_hfbtho_solution(fn=trim(config%file_hfb), ep=hfb%ep, vp=hfb%vp, &
      up=hfb%up, en=hfb%en, vn=hfb%vn, un=hfb%un)

   allocate(occ_n(size(hfb%en)), occ_p(size(hfb%ep)))
   call get_hfbtho_occupations(hfb%vn, occ_n)
   call get_hfbtho_occupations(hfb%vp, occ_p)

   ! If no blocking, exit
   if (config%odd_block_n(1) == 0 .and. config%odd_block_p(1) == 0) then
      stop "No blocking configuration was specified in the input file."
   end if
   if (config%odd_block_n(1) /= 0 .and. config%odd_block_p(1) /= 0) then
      stop "Blocking information specified for both protons and neutrons."
   end if

   write(*,'(a)') "INITIAL-STATE ODD NUCLEON"
   write(*,'(a)') repeat('-', 33)

   ! Find which level to block
   if (config%odd_block_n(1) /= 0) then
      bo_k          = abs(config%odd_block_n(1))
      bo_k_sign     = sign(1, config%odd_block_n(1))
      bo_parity     = config%odd_block_n(2)
      bo_isospin    = IT_NEUTRON
      bo_index      = find_blocked_level(config%odd_block_n(1), bo_parity, bo_isospin, hfb)
      bo_index_r    = qp_reverse(bo_index)
      ! 1 is initial state, 2 is final state
      eqp1  => hfb%en;  eqp2  => hfb%ep
      vmat1 => hfb%vn;  vmat2 => hfb%vp
   else if (config%odd_block_p(1) /= 0) then
      bo_k          = abs(config%odd_block_p(1))
      bo_k_sign     = sign(1, config%odd_block_p(1))
      bo_parity     = config%odd_block_p(2)
      bo_isospin    = IT_PROTON
      bo_index      = find_blocked_level(config%odd_block_p(1), bo_parity, bo_isospin, hfb)
      bo_index_r    = qp_reverse(bo_index)
      eqp1  => hfb%ep;  eqp2  => hfb%en
      vmat1 => hfb%vp;  vmat2 => hfb%vn
   end if

   write(*,'(a,1x,a)')       'Isospin ...... :', trim(merge('Neutron', 'Proton ', bo_isospin == IT_NEUTRON))
   write(*,'(a,1i2,a)')      'K ............ :', 2*nl(bo_index) + ns(bo_index), '/2'
   write(*,'(a,1x,sp,1i0)')  'Parity ....... :', bo_parity
   write(*,'(a,1x,a)')       'PH type ...... :', trim(merge('Particle', 'Hole    ', bo_k_sign == 1))
   write(*,'(a,1f9.6)')      'Lower norm ... :', merge(occ_n(bo_index), occ_p(bo_index), bo_isospin == IT_NEUTRON)
   write(*,'(a,1f9.6,1x,a)') 'QP energy .... :', eqp1(bo_index), 'MeV'

   ! Find the correct final state (p-type vs h-type)
   write(*,*)
   write(*,'(a)') "FINAL-STATE ODD NUCLEON"
   write(*,'(a)') repeat('-', 33)

   bo_final_index = lowest_energy_qp(eqp2, vmat2, 3-bo_isospin, mask_k_positive=.true., &
      mask_occ=(bo_k_sign == -1), mask_unocc=(bo_k_sign == 1))

   write(*,'(a,1x,a)')       'Isospin ...... :', trim(merge('Neutron', 'Proton ', 3-bo_isospin == IT_NEUTRON))
   write(*,'(a,1i2,a)')      'K ............ :', 2*nl(bo_final_index) + ns(bo_final_index), '/2'
   write(*,'(a,1x,a,1i0)')   'Parity ....... :', strpar(npar(bo_final_index)), 1
   write(*,'(a,1x,a)')       'PH type ...... :', trim(merge('Particle', 'Hole    ', bo_k_sign == 1))
   write(*,'(a,1f9.6)')      'Lower norm ... :', merge(occ_n(bo_final_index), occ_p(bo_final_index), 3-bo_isospin == IT_NEUTRON)
   write(*,'(a,1f9.6,1x,a)') 'QP energy .... :', eqp2(bo_final_index), 'MeV'

   write(*,*)
   write(*,'(a)') "ONE-QP CORRECTIONS"
   write(*,'(a)') repeat('-', 33)

   if (bo_isospin == IT_NEUTRON) then
      delta_n = 1-2*occ_n(bo_index)
      delta_z = 0
   else if (bo_isospin == IT_PROTON) then
      delta_n = 0
      delta_z = 1-2*occ_n(bo_index)
   end if

   write(*,'(a,1f11.6,1x,a)') "1-QP energy ...... :", eqp2(bo_final_index)-eqp1(bo_index), 'MeV'
   write(*,'("Neutron number ... :", 1f11.6, " (Initial: ", 1f8.4, ",  Change: ", sp,1f7.4, ")")') &
      sum(occ_n) + delta_n, sum(occ_n), delta_n
   write(*,'("Proton number .... :", 1f11.6, " (Initial: ", 1f8.4, ",  Change: ", sp,1f7.4, ")")') &
      sum(occ_p) + delta_z, sum(occ_p), delta_z
   write(*,'("A final .......... :", 1f11.6, " (Initial: ", 1f8.4, ")")') &
      sum(occ_n) + sum(occ_p) + delta_n + delta_z, sum(occ_n) + sum(occ_p)
   write(*,'("<N-Z> final ...... :", 1f11.6, " (Initial: ", 1f8.4, ")")') &
      sum(occ_n) - sum(occ_p) + delta_n - delta_z, sum(occ_n) - sum(occ_p)

   write(*,*)
   write(*,'(a)') "REF. TWO-QP ENERGY"
   write(*,'(a)') repeat('-', 33)

   ! Calculate lowest 2QP energy in e-e nucleus
   twoqp_init_index  = lowest_energy_qp(hfb%en, hfb%vn, IT_NEUTRON, mask_k_positive=.true., mask_occ=.true.)
   twoqp_final_index = lowest_energy_qp(hfb%ep, hfb%vp, IT_PROTON,  mask_k_positive=.true., mask_unocc=.true.)
   write(*,'(a,1f11.6,1x,a)') "2-QP energy ...... :", hfb%en(twoqp_init_index) + hfb%ep(twoqp_final_index), 'MeV'

end program pnfam_matrix_odd_information
