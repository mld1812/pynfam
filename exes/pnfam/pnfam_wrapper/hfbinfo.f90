!------------------------------------------------------------------------------
! hfbinfo.f90
!
! Print out essential HFB solution information from binary HFBTHO solutions.
! Usage: ./hfbinfo.x [HFBTHO solution filename]
!
! M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2015
!------------------------------------------------------------------------------
program read_HFBTHO
   use constants,        only : strpar, IT_NEUTRON, IT_PROTON
   use blockmatrix_type, only : blockmatrix
   use logger
   use hfb_solution

   implicit none
   
   integer, parameter  :: dp = kind(1.0d0)
   
   character(len=250)    :: hfb_filename
   type(blockmatrix)     :: hfb_vp, hfb_vn, hfb_up, hfb_un
   real(dp), allocatable :: hfb_ep(:), hfb_en(:)
   
   real(dp) :: q, e2qp
   integer  :: k2qp, p2qp, iprot, ineut, iqp, ik, ip, par
   real(dp), allocatable :: occ_n(:), occ_p(:)
   
   openlog  => openlog_stdout
   writelog => writelog_stdout
   closelog => closelog_stdout
   
   !----------------------------------------------------------------------------
   !----------------------------------------------------------------------------
   
   !if (command_argument_count() /= 1) then
   !   write(*,'(a)') "Usage: hfbinfo.x HFBTHO_OUTPUT_FILE"
   !   stop
   !end if
   !call get_command_argument(1, hfb_filename)
   
   call get_hfbtho_solution("", hfb_ep, hfb_vp, hfb_up, hfb_en, hfb_vn, hfb_un)
   
   !! Hack to get at the version number
   !open(11, file=trim(hfb_filename), access='sequential', form='unformatted')
   !read(11) version
   !close(11)
   
   !----------------------------------------------------------------------------
   !----------------------------------------------------------------------------

   ! Results
   write(*,'(1x,a,1i2)') 'HFB STORAGE VERSION:   ', version
   !write(*,'(1x,3a)')    'HFB SOLUTION FILE:    "', trim(hfb_filename), '"'
   
   write(*,*) 
   write(*,'(1x,a,/,1x,a)') 'GLOBAL QUANTITIES', repeat('-',55)
   write(*,'(1x,a,1i12)')         'Particle number (A):   ', hfb_npr(3)
   write(*,'(1x,a, 1f12.6,1x,a)') 'Pairing window (PWI):  ', hfb_pairing_window, 'MeV'
   write(*,'(1x,a, 1f12.6)')      'Total deformation:     ', hfb_beta2
   write(*,'(1x,a, 1f12.6,1x,a)') 'Total energy:          ', hfb_energy,      'MeV'
   
   write(*,*) 
   write(*,*) 
   write(*,'(1x,a,8x,a)') 'NEUTRON/PROTON QUANTITIES', 'Neutrons       Protons'
   write(*,'(1x,a)') repeat('-',55)
   write(*,'(1x,a,2i14)')        'Particle number (N,Z):     ', hfb_npr(1:2)
   write(*,'(1x,a,2f14.6,1x,a)') 'Lambda (ALA):              ', hfb_lambda,  'MeV'
   write(*,'(1x,a,2f14.6,1x,a)') 'Lambda (ALAST):            ', hfb_elast,   'MeV'
   write(*,'(1x,a,2f14.6,1x,a)') 'Avg. pairing gaps (Delta): ', hfb_delta,   'MeV'
   write(*,'(1x,a,2f14.6,1x,a)') 'Pairing strengths:         ', hfb_cpair,   'MeV-fm**(-3)'
   write(*,'(1x,a,2f14.1)')      'Pairing mixture (alpha):   ', hfb_alpha_pair
   
   write(*,*) 
   write(*,*) 
   write(*,'(1x,a,/,1x,a)') 'HFBTHO SOLUTION PROPERTIES', repeat('-',55)
   write(*,'(1x,a,1i10)') 'Num. blocks (nb): ', nb/2
   write(*,'(1x,a,1i10)') 'Num. blocks (nbx):', nbx/2
   write(*,'(1x,a,1i10)') 'Dim. En, Ep:      ', sum(db(1:nb/2))
   write(*,'(1x,a,1i10)',advance='no') 'Active levels (n):', sum(pwi_dim_n(1:nb/2))
   write(*,'(4x,"(", 1i0, " excluded)")') sum(db(1:nb/2))-sum(pwi_dim_n(1:nb/2))
   write(*,'(1x,a,1i10)',advance='no') 'Active levels (p):', sum(pwi_dim_p(1:nb/2))
   write(*,'(4x,"(", 1i0, " excluded)")') sum(db(1:nb/2))-sum(pwi_dim_p(1:nb/2))
   write(*,'(1x,a,1i10)') 'Dim. Vn, Vp:      ', sum(db(1:nb/2)**2)
   write(*,'(1x,a,1i10)') 'NGHL:             ', nghl
   
   ineut = lowest_energy_qp(hfb_en, hfb_vn, IT_NEUTRON, mask_k_positive=.true., mask_occ=.true.,  mask_unocc=.false.)
   iprot = lowest_energy_qp(hfb_ep, hfb_vp, IT_PROTON,  mask_k_positive=.true., mask_occ=.false., mask_unocc=.true.)

   allocate(occ_n(sum(db)), occ_p(sum(db)))
   call get_hfbtho_occupations(hfb_vn, IT_NEUTRON, occ_n)
   call get_hfbtho_occupations(hfb_vp, IT_PROTON, occ_p)

   call get_hfbtho_betadecay_properties(hfb_ep, hfb_en, hfb_vp, hfb_vn, q, e2qp, k2qp, p2qp)
  
   write(*,*)
   write(*,'(1x,a,1f10.6,1x,a)',advance='no') 'Min. 2qp energy:   ', e2qp, 'MeV'
   write(*,'(4x,"(",1i0,a,")")') k2qp, merge('+', '-', p2qp > 0)
   write(*,'(1x,a,1f10.6,1x,a)') 'Computed Q-value:  ', q, 'MeV'
   write(*,'(1x,a,1f10.6,1x,a)') 'Computed EQRPA_max:', q+e2qp, 'MeV'
   
   write(*,*) 
   write(*,*) 
   write(*,'(1x,a,/,1x,a)') 'HFBTHO COUPLING CONSTANTS (S=0, T=1)', repeat('-',55)
   write(*,'(1x,2a)')      'J**2 terms enabled? ', merge('Yes', 'No ', hfb_use_j2terms)
   write(*,'(1x,a,1f24.16)') 'Crho0:', hfb_cr0
   write(*,'(1x,a,1f24.16)') 'Crhor:', hfb_crr
   write(*,'(1x,a,1f24.16)') 'CDrho:', hfb_cdrho
   write(*,'(1x,a,1f24.16)') 'Ctau: ', hfb_ctau
   write(*,'(1x,a,1f24.16)') 'CJ:   ', hfb_ctj
   write(*,'(1x,a,1f24.16)') 'CrDJ: ', hfb_crdj
   
   ! Finite temperature
   write(*,*) 
   write(*,*) 
   write(*,'(1x,a,/,1x,a)') 'HFBTHO FINITE-TEMPERATURE', repeat('-',55)
   write(*,'(1x,2a)')       'Finite temperature active? ', merge('Yes', 'No ', ft_active)
   
   if (ft_active) then
      write(*,'(1x,a,1f10.4,1x,a)') 'Temperature:', ft_temp, 'MeV'
      write(*,'(1x,a,1f14.8)')      'Entropy (n):', ft_entropy(1)
      write(*,'(1x,a,1f14.8)')      'Entropy (p):', ft_entropy(2)
      write(*,'(1x,a,1f14.8)')      'Entropy (t):', ft_entropy(3)
   
      write(*,*)
      write(*,'(1x,a,1f10.6,"  (E_p = ",1f0.4," MeV)")') 'Approx. largest f_p:', &
         maxval(ft_fp), hfb_ep(pwi_qp_p(maxloc(ft_fp(1:sum(pwi_dim_p)))))
      write(*,'(1x,a,1f10.6,"  (E_n = ",1f0.4," MeV)")') 'Approx. largest f_n:', &
         maxval(ft_fn), hfb_en(pwi_qp_n(maxloc(ft_fn(1:sum(pwi_dim_n)))))
   end if
   
   ! Odd nuclei
   write(*,*) 
   write(*,*) 
   write(*,'(1x,a,/,1x,a)') 'HFBTHO ODD-MASS NUCLEI', repeat('-',55)
   write(*,'(1x,2a)') 'Equal-filling active? ', merge('Yes', 'No ', hfb_blo_active)
   
   ! Neutrons
   if (hfb_blo_qp(1) <= 0) then
      write(*,'(1x,a)') 'Neutron odd ...... : No'
   else
      write(*,'(1x,a,"(",1i0,"/2",a,", ",1f8.4," MeV)")') 'Neutron odd ...... : Yes ', &
         2*nl(hfb_blo_qp(1))+ns(hfb_blo_qp(1)), &
         strpar(npar(hfb_blo_qp(1))), &
         hfb_en(hfb_blo_qp(1))
      write(*, '(1x,a,1x,1i0)') 'Blocking Candidate :', keyblo(1)
      write(*,'(1x,a)') 'Ex[B-] (Ep-En) ... :'
      
      do ip=1, 2
         ! Parity is the same (ip=1) or 3-parity (ip=2)
         par = (ip-1)*3 + (-1)**(ip-1)*npar(hfb_blo_qp(1))
         do ik=0, 2
            iqp = lowest_energy_qp(hfb_ep, hfb_vp, IT_PROTON,     &
               mask_k_positive=.true.,                            &
               mask_2k=ik+2*nl(hfb_blo_qp(1))+ns(hfb_blo_qp(1)),  &
               mask_parity=par)
            write(*,'(5x,1i0,a,4x,1f9.6)') ik, strpar(ip), hfb_ep(iqp)-hfb_en(hfb_blo_qp(1))
         end do
      end do
   end if
   
   ! Protons
   if (hfb_blo_qp(2) <= 0) then
      write(*,'(1x,a)') 'Proton odd ....... : No'
   else
      write(*,'(1x,a,"(",1i0,"/2",a,", ",1f8.4," MeV)")') 'Proton odd ....... : Yes ', &
         2*nl(hfb_blo_qp(2))+ns(hfb_blo_qp(2)), &
         strpar(npar(hfb_blo_qp(2))), &
         hfb_ep(hfb_blo_qp(2))
      write(*, '(1x,a,1x,1i0)') 'Blocking Candidate :', keyblo(2)
      write(*,'(1x,a)') 'Ex[B-] (En-Ep) ... :'
      
      do ip=1, 2
         ! Parity is the same (ip=1) or 3-parity (ip=2)
         par = (ip-1)*3 + (-1)**(ip-1)*npar(hfb_blo_qp(2))
         do ik=0, 2
            iqp = lowest_energy_qp(hfb_en, hfb_vn, IT_NEUTRON,    &
               mask_k_positive=.true.,                            &
               mask_2k=ik+2*nl(hfb_blo_qp(2))+ns(hfb_blo_qp(2)),  &
               mask_parity=par)
            write(*,'(5x,1i0,a,4x,1f9.6)') ik, strpar(ip), hfb_en(iqp)-hfb_ep(hfb_blo_qp(2))
         end do
      end do
   end if

end program read_HFBTHO
