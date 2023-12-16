!------------------------------------------------------------------------------
!> This module contains the routines for reading and writing the pnfam solution
!> to and from binary outptu files. It is largely borrowed from HFBTHO but
!> modified for pnfam.
!>
!> @authors E.M. Ney, UNC Chapel Hill, 2018-
!------------------------------------------------------------------------------
module pnfam_storage
  use pnfam_setup
  use pnfam_interaction
  use pnfam_broyden, only : qrpa_broin
  use type_gamdel_2bc
  !use pnfam_hamiltonian, only : rhox
  !use type_densityset

  Implicit None

  integer, parameter, private :: ipr=kind(1), dp=kind(1d0)

  Integer, PUBLIC, SAVE :: VERSION_DATA = 1 !< Version number of the output binary file.
  Integer, PUBLIC, SAVE :: VERSION_READ = 1 !< Version number of the input binary file.
  Character(Len=50), PRIVATE, SAVE :: welfile, tbcfile !< Name of the binary file
  Integer, public, parameter :: lwou=14, lwin=15

Contains
  !=======================================================================
  !> Subroutine that defines a filename
  !=======================================================================
  Subroutine FileLabels(npri,filelabel)
    Implicit None
    Integer(ipr) :: it,nprt(2),npri(3)
    Character(3)  :: snpr(2)
    Character(8)  :: filelabel
    !
    nprt=npri(1:2)
    Do it=1,2
       Write(snpr(it),'(i3.3)') nprt(it)
       If(nprt(it).Lt.10  ) Then
          Write(snpr(it),'(a2,i1)') '00',nprt(it)
       Else
          If(nprt(it).Lt.100 ) Then
             Write(snpr(it),'(a1,i2)') '0',nprt(it)
          Else
             Write(snpr(it),'(i3)') nprt(it)
          End If
       End If
    End Do
    !
    Write(filelabel,'(a3,a1,a3)')  snpr(1),'_',snpr(2)
    Write(welfile,'(a)') trim(output_binfile)
    Write(tbcfile,'(a)') trim(extfield_binfile)
    !
  End Subroutine FileLabels
  !=======================================================================
  !> This subroutine is the central interface to all input/output
  !> operations in HFBTHO. Depending on the value of its input parameter
  !> it will either read a binary file to get input data, or it will write
  !> the current data in a binary file. In the reading phase, the code
  !> tries to read a version number. If this read is successful, it is
  !> assumed the binary file is in version VERSION_DATA
  !=======================================================================
  Subroutine fam_io(is_in,iexit)
    Implicit None
    Integer(ipr), INTENT(IN) :: is_in !> Integer specifying if the file is read (=1) or written (=2)
    Integer(ipr), INTENT(INOUT) :: iexit !> Integer giving the exit status of the routine (0: OK, >0: not OK)
    Integer(ipr) :: checked,checked_again
    Character(50) :: action, filename
    Character(8) :: filelabel
    Integer(ipr) :: is
    is = abs(is_in)
    ! label organization
    Call FileLabels(hfb_npr,filelabel)
    !---------------------------------------------------------------------
    ! Read data to start the calculation
    !---------------------------------------------------------------------
    ! Valid inputs are +/-(1 or 2)
    if (is /= 1 .and. is /= 2) then
       iexit = 1
       return
    end if
    ! Positive values used for fam soln, negatives for tbc
    if (is_in > 0) then
       filename = welfile
    else
       filename = tbcfile
    end if
    If(is.Eq.1) Then
       action = 'Read'
       ! Start from scratch
       If(use_fam_storage == 0 .and. is_in > 0) Then
          iexit=1
       Else
          ! Check the file is valid
          checked = check_file(filename, action, is)
          If(checked==0) Then
             ! Read the version number of the file
             VERSION_READ = version_number(filename)
             ! If version is current, try to read the data
             If(VERSION_READ==VERSION_DATA) Then
                iexit=0
                if (is_in > 0) then
                   Call read_data(iexit)
                else
                   Call read_tbc(iexit)
                end if
             Else
                iexit=1
             End If
          Else
             iexit=1
          End If
       End If
    End If
    !---------------------------------------------------------------------
    ! Write data on disk
    !---------------------------------------------------------------------
    If(is.Eq.2) Then
       action = 'Write'
       If (use_fam_storage == -1 .and. filename /= "" .or. is_in < 0) Then
          ! Check the file is valid
          checked = check_file(filename, action, is)
          If(checked==0) Then
             Call write_version()
             if (is_in > 0) then
                Call write_data()
             else
                Call write_tbc()
             end if
             iexit=0
          Else
             iexit=1
          End If
       Else
          iexit=0
       End If
    End If
  End Subroutine fam_io
  !=======================================================================
  !> Subroutine that just writes the version number on disk
  !=======================================================================
  Subroutine write_version()
    Implicit None
    Write(lwou) VERSION_DATA
  End Subroutine write_version
  !=======================================================================
  !> Function that reads and returns the version number from the input
  !> file. Returns -1 in case the version cannot be read.
  !=======================================================================
  Integer(ipr) Function version_number(filename)
    Implicit None
    Character(50), INTENT(IN) :: filename
    Integer(ipr) :: ierr,version
    version = 1
    Read(lwin,IOSTAT=ierr) version
    If(ierr.NE.0) Then
       version = -1
    End If
    version_number = version
  End Function version_number
  !=======================================================================
  !> Function that checks the status of the file. If the file exists and
  !> can be opened, it is opened and the functions returns an exit status
  !> of 0. If the file should be read but does not exist, the function
  !> returns 1; if the file should be written but does not exist, the
  !> function opens a new file and returns 0.
  !=======================================================================
  Integer(ipr) Function check_file(filename, action, is)
    Implicit None
    Character(50), INTENT(IN) :: filename !> Name of the file to check
    Character(50), INTENT(IN) :: action !> Equal to 'Read' or 'Write', defines what to do if the file does not exist
    Integer(ipr), INTENT(IN) :: is !> 1=read, 2=write pre-hfbdiag quantities, 3=write post-hfbdiag quantities
    !
    Logical :: file_exists,file_opened
    Integer(ipr) :: ierr,iexit
    iexit = 0
    file_exists=.False.; inquire(file=filename, exist=file_exists); ierr=0
    If(Trim(action) == 'Read') Then
       If(file_exists) Then
          file_opened=.False.; inquire(unit=lwin, opened=file_opened)
          If(file_opened) Then
             Close(lwin)
          End If
          Open(lwin,file=filename,status='old',form='unformatted',IOSTAT=ierr)
          If(ierr.NE.0) Then
             iexit = 1
          End If
       Else
          iexit = 1
       End If
    End If
    If(Trim(action) == 'Write') Then
       If(file_exists) Then
          file_opened=.False.; inquire(unit=lwou, opened=file_opened)
          If(file_opened) Then
             Close(lwou)
          End If
          If(is.Eq.3) Then ! Append to existing file
             Open(lwou,file=filename,status='old',form='unformatted',position='append',IOSTAT=ierr)
          Else             ! Overwrite existing file
             Open(lwou,file=filename,status='old',form='unformatted',IOSTAT=ierr)
          End If
          If(ierr.NE.0) Then
             iexit = 1
          End If
       Else
          If(is.Eq.3) Then ! File should exist
              iexit = 1
          Else
              Open(lwou,file=filename,status='new',form='unformatted',IOSTAT=ierr)
              iexit = 0
          End If
       End If
    End If
    check_file = iexit
  End Function check_file
  !=======================================================================
  !> This subroutine reads a binary file containing the results of a
  !> HFBTHO calculation. In the new format (VERSION_DATA=3), the binary
  !> file is structured by keywords and contains the HF and pairing field
  !> on the Gauss quadrature mesh. This new format allows restarts even
  !> when (i) the basis is different (as long as the quadrature mesh is
  !> the same) (ii) the number of constraints is different (iii) metadata
  !> such as number of protons, neutrons, characteristics of the EDF, etc.
  !> are different. The price to pay for this flexibility is an increased
  !> size of the binary file.
  !=======================================================================
  Subroutine read_data(iexit)
    use pnfam_constants, only : iu
    Implicit None
    Integer(ipr), INTENT(INOUT) :: iexit
    Character(Len=8) :: key

    integer :: Z_r,N_r
    integer :: nb_r
    integer, allocatable :: db_r(:)
    integer :: ft_active_r
    integer :: hfb_blo_active_r
    real(dp) :: ft_temp
    integer :: hfb_blo_qpn_r, hfb_blo_qpp_r
    character(len=80) :: interaction_name_r
    integer :: require_gauge_invariance_r
    integer :: require_self_consistency_r
    integer :: force_j2_terms_r
    real(dp) :: vpair_t0_r, vpair_t1_r
    real(dp) :: override_cs0_r, override_csr_r, override_cds_r, override_ct_r
    real(dp) :: override_cf_r, override_cgs_r, override_cj_r, override_csdj_r
    integer :: nghl_r
    real(dp), allocatable :: crho_r(:)
    real(dp) :: cdrho_r, ctau_r, ctj0_r, ctj1_r, ctj2_r, crdj_r
    real(dp), allocatable :: cs_r(:)
    real(dp) :: cds_r, ct_r, cj_r, cgs_r, cf_r, csdj_r
    real(dp), allocatable :: cpair_r(:)
    real(dp), allocatable :: cspair_r(:)
    character(len=80) :: operator_name_r
    integer :: operator_k_r
    integer :: real_eqrpa_r, imag_eqrpa_r
    real(dp) :: quench_residual_int_r
    integer :: nxy_r, ibroy, ibroy_r
    !integer :: ngh_r, ngl_r
    integer :: nxterms, ixterm, nxterms_r
    real(dp), allocatable :: qrpa_broin_r(:)
    ! Populated directly from read: si, dRqp_re%m*%elem, dRqp_im%m*%elem

    !---------------------------------------------------------------------
    ! Read data
    !---------------------------------------------------------------------
    write(st,'(1x,2a)') " * Reading from storage file: ",trim(welfile)
    call writelog(st)
    call writelog("")
    iexit = 0
    ! Loop over all keywords in the file
    Do
       Read(lwin,Err=100,End=99) key
       If(Trim(key)=='hfb_soln') Then
          Read(lwin,Err=100,End=100) Z_r,N_r
          Read(lwin,Err=100,End=100) nb_r
          If(Allocated(db_r)) Deallocate(db_r)
          Allocate(db_r(nb_r))
          Read(lwin,Err=100,End=100) db_r
          ! Enforce same matrix size
          If(nb/=nb_r.or. .not. all(db==db_r)) iexit=1
       End If
       If(Trim(key)=='stat_hfb') Then
          Read(lwin,Err=100,End=100) ft_active_r
          Read(lwin,Err=100,End=100) hfb_blo_active_r
          Read(lwin,Err=100,End=100) ft_temp
          Read(lwin,Err=100,End=100) hfb_blo_qpn_r, hfb_blo_qpp_r
       End If
       If(Trim(key)=='interact') Then
          ! Namelist parameters (overrides used to construct couplings for hamiltonian)
          Read(lwin,Err=100,End=100) interaction_name_r
          Read(lwin,Err=100,End=100) require_gauge_invariance_r
          Read(lwin,Err=100,End=100) require_self_consistency_r
          Read(lwin,Err=100,End=100) force_j2_terms_r
          Read(lwin,Err=100,End=100) vpair_t0_r, vpair_t1_r
          Read(lwin,Err=100,End=100) override_cs0_r, override_csr_r, override_cds_r, override_ct_r
          Read(lwin,Err=100,End=100) override_cf_r, override_cgs_r, override_cj_r, override_csdj_r
          ! Couplings used by hamiltonian (check if nghl consistent?)
          Read(lwin,Err=100,End=100) nghl_r
          If (Allocated(crho_r)) Deallocate(crho_r, cs_r, cpair_r, cspair_r)
          Allocate(crho_r(nghl_r), cs_r(nghl_r), cpair_r(nghl_r), cspair_r(nghl_r))
          ! Time even
          Read(lwin,Err=100,End=100) crho_r
          Read(lwin,Err=100,End=100) cdrho_r, ctau_r, ctj0_r, ctj1_r, ctj2_r, crdj_r
          ! Time odd
          Read(lwin,Err=100,End=100) cs_r
          Read(lwin,Err=100,End=100) cds_r, ct_r, cj_r, cgs_r, cf_r, csdj_r
          ! Pairing
          Read(lwin,Err=100,End=100) cpair_r
          Read(lwin,Err=100,End=100) cspair_r
       End If

       If(Trim(key)=='extfield') Then
          Read(lwin,Err=100,End=100) operator_name_r
          Read(lwin,Err=100,End=100) operator_k_r
          ! Enforce same operator, otherwise block structure and nxy will differ
          call translate_uppercase(operator_name_r)
          call translate_uppercase(operator_name)
          If (trim(operator_name_r)/=trim(operator_name).or.operator_k_r/=operator_k) iexit=1
       End If

       If(Trim(key)=='fam_soln') Then
          Read(lwin,Err=100,End=100) real_eqrpa_r, imag_eqrpa_r
          Read(lwin,Err=100,End=100) quench_residual_int_r
          Read(lwin,Err=100,End=100) si ! True variable

          ! Only read in all the solution if size matches storage
          Read(lwin,Err=100,End=100) nxy_r
          If (nxy_r==nxy) Then
             ibroy=4; ibroy_r=4
             If (ft_active.or.hfb_blo_active) ibroy=8
             If (ft_active_r==1.or.hfb_blo_active_r==1) ibroy_r=8

             If (Allocated(qrpa_broin_r)) Deallocate(qrpa_broin_r)
             Allocate(qrpa_broin_r(ibroy_r*nxy_r))
             Read(lwin,Err=100,End=100) qrpa_broin_r
             ! Allow mismatch in active blocks
             If (ibroy == ibroy_r) Then
                qrpa_broin = qrpa_broin_r
             Else If (ibroy < ibroy_r) Then
                qrpa_broin = qrpa_broin_r(1:ibroy*nxy_r)
             Else
                qrpa_broin(1:ibroy*nxy) = qrpa_broin_r
             End If
             Deallocate(qrpa_broin_r)
          Else
             iexit=1
             Read(lwin,Err=100,End=100)
          End If

          !! THIS IS ONLY PN DENSITIES, ALSO NEED NP
          !Read(lwin,Err=100,End=100) ngh_r, ngl_r
          !call allocate_density(rhox,ngh_r*ngl_r)
          !Read(lwin,Err=100,End=100) rhox%rerho,  rhox%imrho
          !Read(lwin,Err=100,End=100) rhox%retau,  rhox%imtau
          !Read(lwin,Err=100,End=100) rhox%retjrr, rhox%imtjrr, rhox%retjrp, rhox%imtjrp, rhox%retjrz, rhox%imtjrz
          !Read(lwin,Err=100,End=100) rhox%retjpr, rhox%imtjpr, rhox%retjpp, rhox%imtjpp, rhox%retjpz, rhox%imtjpz
          !Read(lwin,Err=100,End=100) rhox%retjzr, rhox%imtjzr, rhox%retjzp, rhox%imtjzp, rhox%retjzz, rhox%imtjzz
          !Read(lwin,Err=100,End=100) rhox%resr,   rhox%imsr,   rhox%resp,   rhox%imsp,   rhox%resz,   rhox%imsz
          !Read(lwin,Err=100,End=100) rhox%retr,   rhox%imtr,   rhox%retp,   rhox%imtp,   rhox%retz,   rhox%imtz
          !Read(lwin,Err=100,End=100) rhox%rejr,   rhox%imjr,   rhox%rejp,   rhox%imjp,   rhox%rejz,   rhox%imjz
          !Read(lwin,Err=100,End=100) rhox%refr,   rhox%imfr,   rhox%refp,   rhox%imfp,   rhox%refz,   rhox%imfz
          !Read(lwin,Err=100,End=100) rhox%regs,   rhox%imgs
          !Read(lwin,Err=100,End=100) rhox%rerb,   rhox%imrb
          !Read(lwin,Err=100,End=100) rhox%resbr,  rhox%imsbr,  rhox%resbp,  rhox%imsbp,  rhox%resbz,  rhox%imsbz
          !If (ngh_r/=ngh .or. ngl_r/=ngl) Then
          !   iexit=1
          !   Write(st,'(1x,a)') " * Inconsistent quadrature mesh."; call writelog(st)
          !   call deallocate_density(rhox)
          !End If

          ! Only read crossterms if nxterms matches storage
          nxterms = size(re_str) - 1
          Read(lwin,Err=100,End=100) re_str(1), im_str(1)
          Read(lwin,Err=100,End=100) nxterms_r
          If (nxterms_r > 0) Then
             Do ixterm=1, nxterms_r
                If (nxterms_r==nxterms) Then
                   Read(lwin,Err=100,End=100) re_str(1+ixterm), im_str(1+ixterm)
                Else
                   Read(lwin,Err=100,End=100)
                End If
             End Do
          End If
       End If
    End Do

 99 Continue
    Close(lwin)
    If (iexit==1) Then
       Write(st,'(1x,a,a,a)')  " * The file '",trim(welfile),"' is incompatible with the current calculation!" ; call writelog(st)
       Write(st,'(1x,a)')      " * Starting calculation from scratch..."; call writelog(st)
       call writelog("")
    End If
    Return
    !
100 Continue
    iexit=1
    Close(lwin)
    Write(st,'(1x,a,a,a)')  " * The file '",trim(welfile),"' is corrupted!" ; call writelog(st)
    Write(st,'(1x,a,a8,a)') " * Problem occurs for key: ",trim(key) ; call writelog(st)
    Write(st,'(1x,a)')      " * Starting calculation from scratch..."; call writelog(st)
    call writelog("")
    !
  End Subroutine read_data
  !=======================================================================
  !> This subroutine writes a binary file containing the results of a
  !> HFBTHO calculation. In the new format (VERSION_DATA=3), the binary
  !> file is structured by keywords and contains the HF and pairing field
  !> on the Gauss quadrature mesh.
  !=======================================================================
  Subroutine write_data()
    Implicit None

    Integer(ipr) :: ixterm, nxterms

    write(lwou) 'hfb_soln'
    write(lwou) hfb_npr(1), hfb_npr(2)
    write(lwou) nb
    write(lwou) db

    write(lwou) 'stat_hfb'
    write(lwou) merge(1,0,ft_active)
    write(lwou) merge(1,0,hfb_blo_active)
    write(lwou) ft_temp
    write(lwou) hfb_blo_qp(1),hfb_blo_qp(2)

    write(lwou) 'interact'
    ! Namelist parameters (overrides used to construct couplings for hamiltonian)
    write(lwou) interaction_name
    write(lwou) merge(1,0,require_gauge_invariance)
    write(lwou) merge(1,0,require_self_consistency)
    write(lwou) merge(1,0,force_j2_terms)
    write(lwou) vpair_t0, vpair_t1
    write(lwou) override_cs0, override_csr, override_cds, override_ct
    write(lwou) override_cf, override_cgs, override_cj, override_csdj
    ! Couplings used by hamiltonian
    write(lwou) size(crho)
    ! Time even
    write(lwou) crho
    write(lwou) cdrho, ctau, ctj0, ctj1, ctj2, crdj
    ! Time odd
    write(lwou) cs
    write(lwou) cds, ct, cj, cgs, cf, csdj
    ! Pairing
    write(lwou) cpair
    write(lwou) cspair

    write(lwou) 'extfield'
    write(lwou) operator_name
    write(lwou) operator_k

    write(lwou) 'fam_soln'
    write(lwou) real_eqrpa, imag_eqrpa
    write(lwou) quench_residual_int
    write(lwou) si
    write(lwou) nxy
    write(lwou) qrpa_broin

    !! THIS IS ONLY PN DENSITIES, ALSO NEED NP
    !write(lwou) ngh, ngl
    !write(lwou) rhox%rerho,  rhox%imrho
    !write(lwou) rhox%retau,  rhox%imtau
    !write(lwou) rhox%retjrr, rhox%imtjrr, rhox%retjrp, rhox%imtjrp, rhox%retjrz, rhox%imtjrz
    !write(lwou) rhox%retjpr, rhox%imtjpr, rhox%retjpp, rhox%imtjpp, rhox%retjpz, rhox%imtjpz
    !write(lwou) rhox%retjzr, rhox%imtjzr, rhox%retjzp, rhox%imtjzp, rhox%retjzz, rhox%imtjzz
    !write(lwou) rhox%resr,   rhox%imsr,   rhox%resp,   rhox%imsp,   rhox%resz,   rhox%imsz
    !write(lwou) rhox%retr,   rhox%imtr,   rhox%retp,   rhox%imtp,   rhox%retz,   rhox%imtz
    !write(lwou) rhox%rejr,   rhox%imjr,   rhox%rejp,   rhox%imjp,   rhox%rejz,   rhox%imjz
    !write(lwou) rhox%refr,   rhox%imfr,   rhox%refp,   rhox%imfp,   rhox%refz,   rhox%imfz
    !write(lwou) rhox%regs,   rhox%imgs
    !write(lwou) rhox%rerb,   rhox%imrb
    !write(lwou) rhox%resbr,  rhox%imsbr,  rhox%resbp,  rhox%imsbp,  rhox%resbz,  rhox%imsbz

    write(lwou) re_str(1), im_str(1)
    nxterms = size(re_str)-1
    write(lwou) nxterms
    If (nxterms > 0) Then
       Do ixterm=1, nxterms
          write(lwou) re_str(1+ixterm), im_str(1+ixterm)
       End Do
    End If

    Close(lwou)

  End Subroutine write_data

  Subroutine write_tbc()
    Implicit None

    Integer(ipr) :: ixterm, fnxy, fnxy12, use_hblas
    real(dp) :: c3used, c4used

    fnxy = size(f%mat%elem)
    fnxy12 = -1
    if (allocated(f%mat12%elem)) fnxy12 = size(f%mat12%elem)
    c3used = two_body_current_lecs(1)
    c4used = two_body_current_lecs(2) + 0.25_dp

    use_hblas = 0
#ifdef USE_HBLAS
#if USE_HBLAS==1
    use_hblas = 1
#endif
#endif

    write(lwou) 'extf_2bc'
    write(lwou) use_hblas
    write(lwou) two_body_current_mode
    write(lwou) merge(1,0,two_body_current_usep)
    write(lwou) two_body_current_lecs
    write(lwou) nb, dqp, fnxy, fnxy12
    write(lwou) f%label
    write(lwou) f%k
    write(lwou) f%rank
    write(lwou) merge(1,0,f%beta_minus)
    write(lwou) merge(1,0,f%parity_even)
    write(lwou) f%use_2bc

    ! Matrix elements (without LECs)
    call strip_gamdel_lecs(f%gam, c3used, c4used)
    write(lwou) f%gam%c3d
    write(lwou) f%gam%c3e
    write(lwou) f%gam%c4d
    write(lwou) f%gam%c4e
    if (two_body_current_usep) then
       write(lwou) f%gam%cpd
       write(lwou) f%gam%cpe
    end if
    if (fnxy12 > 0) then
       call strip_gamdel_lecs(f%del, c3used, c4used)
       write(lwou) f%del%c3d
       write(lwou) f%del%c3e
       write(lwou) f%del%c4d
       write(lwou) f%del%c4e
       if (two_body_current_usep) then
          write(lwou) f%del%cpd
          write(lwou) f%del%cpe
       end if
    end if
    call deallocate_gamdel(f%gam)
    call deallocate_gamdel(f%del)

    write(lwou) nxterms
    !If (nxterms > 0) Then
    !   Do ixterm=1, nxterms
    !      write(lwou) g(ixterm)%label
    !      write(lwou) g(ixterm)%k
    !      write(lwou) g(ixterm)%rank
    !      write(lwou) merge(1,0,g(ixterm)%beta_minus)
    !      write(lwou) merge(1,0,g(ixterm)%parity_even)
    !      write(lwou) g(ixterm)%use_2bc
    !      write(lwou) g(ixterm)%mat%elem
    !      if (fnxy > 0) write(lwou) g(ixterm)%mat12%elem
    !   End Do
    !End If

    Close(lwou)

  End Subroutine write_tbc

  Subroutine read_tbc(iexit)
    use pnfam_constants, only : iu
    use type_gamdel_2bc
    Implicit None
    Integer(ipr), INTENT(INOUT) :: iexit
    Character(Len=8) :: key

    Integer(ipr) :: ixterm, fnxy, fnxy12, use_hblas, i
    Integer(ipr) :: nb_r, dqp_r, nxy_r, nxy12_r
    Character(len=50) :: label_r
    Integer :: k_r, rank_r, beta_minus_r, parity_even_r
    Integer :: use_2bc_r(3)
    Real(dp), allocatable, dimension(:) :: elem_r, elem12_r
    Integer :: tbc_mode_r, tbc_usep_r, use_hblas_r
    Real(dp) :: tbc_lecs_r(3)
    type(gamdel_2bc) :: gam_r, del_r
    Real(dp) :: c3new, c4new

    write(st,'(1x,2a)') " * Reading from extfield file: ",trim(tbcfile)
    call writelog(st)
    call writelog("")
    iexit = 0

    fnxy = size(f%mat%elem)
    fnxy12 = -1
    if (allocated(f%mat12%elem)) fnxy12 = size(f%mat12%elem)
    c3new = two_body_current_lecs(1)
    c4new = two_body_current_lecs(2) + 0.25_dp

    use_hblas = 0
#ifdef USE_HBLAS
#if USE_HBLAS==1
    use_hblas = 1
#endif
#endif

    Do
       Read(lwin,Err=101,End=98) key
       If(Trim(key)=='extf_2bc') Then
          ! 2 body current mode
          !   - tbc_mode, 1st digit doesn't matter b/c 1BC and 2BC are separate
          !   - tbc_mode, 2nd digit doesn't matter b/c 1=file, else=1BC code
          !   - tbc_mode, 3rd digit doesn't matter b/c gam and del are stored separately
          Read(lwin,Err=101,End=101) use_hblas_r ! If mismatch, could just reorder... (not implemented yet)
          Read(lwin,Err=101,End=101) tbc_mode_r  ! If mismatch, doesn't matter, see above
          Read(lwin,Err=101,End=101) tbc_usep_r  ! If mismatch, NEED RERUN only if yes and not in file
          Read(lwin,Err=101,End=101) tbc_lecs_r  ! If mismatch, doesn't matter, uses new ones

          if (use_hblas_r /= use_hblas .or. &
              (tbc_usep_r==0 .and. two_body_current_usep)) then
             iexit = 1
          end if

          ! Matrix structure
          Read(lwin,Err=101,End=101) nb_r, dqp_r, nxy_r, nxy12_r
          if (nb_r /= nb .or. dqp_r /= dqp .or.&
              nxy_r /= fnxy .or. nxy12_r /= fnxy12) then
             iexit = 1
          end if

          ! Operator (skip check on 2bc mode)
          Read(lwin,Err=101,End=101) label_r
          Read(lwin,Err=101,End=101) k_r
          Read(lwin,Err=101,End=101) rank_r
          Read(lwin,Err=101,End=101) beta_minus_r
          Read(lwin,Err=101,End=101) parity_even_r
          Read(lwin,Err=101,End=101) use_2bc_r
          if ((trim(label_r) /= f%label) .or. &
              k_r /= f%k .or. rank_r /= f%rank .or. &
              beta_minus_r /= merge(1,0,f%beta_minus) .or. &
              parity_even_r /= merge(1,0,f%parity_even)) then
             iexit = 1
          end if

          ! Matrix elements
          call allocate_gamdel(gam_r, nxy_r)
          call zero_gamdel(gam_r)
          Read(lwin,Err=101,End=101) gam_r%c3d
          Read(lwin,Err=101,End=101) gam_r%c3e
          Read(lwin,Err=101,End=101) gam_r%c4d ! This is always zero...
          Read(lwin,Err=101,End=101) gam_r%c4e
          if (tbc_usep_r/=0) then
             Read(lwin,Err=101,End=101) gam_r%cpd
             Read(lwin,Err=101,End=101) gam_r%cpe
          end if
          if (nxy_r == fnxy) then
             call mult_gamdel_lecs(gam_r, c3new, c4new)
             if (f%use_2bc(2)==5) then
                call calc_totgamdel(gam_r, f%mat%elem, 'd') ! direct part only
             else
                call calc_totgamdel(gam_r, f%mat%elem)
             end if
          end if

          if (nxy12_r > 0) then
             call allocate_gamdel(del_r, nxy12_r)
             call zero_gamdel(del_r)
             Read(lwin,Err=101,End=101) del_r%c3d
             Read(lwin,Err=101,End=101) del_r%c3e
             Read(lwin,Err=101,End=101) del_r%c4d
             Read(lwin,Err=101,End=101) del_r%c4e
             if (tbc_usep_r/=0) then
                Read(lwin,Err=101,End=101) del_r%cpd
                Read(lwin,Err=101,End=101) del_r%cpe
             end if
             if (fnxy12 > 0 .and. nxy12_r == fnxy12) then
                call mult_gamdel_lecs(del_r, c3new, c4new)
                call calc_totgamdel(del_r, f%mat12%elem)
             end if
          end if
          call deallocate_gamdel(gam_r)
          call deallocate_gamdel(del_r)

          Read(lwin,Err=101,End=101) nxterms
          !If (nxterms > 0) Then
          !   Do ixterm=1, nxterms
          !      ! Operator
          !      Read(lwin,Err=101,End=101) label_r
          !      Read(lwin,Err=101,End=101) k_r
          !      Read(lwin,Err=101,End=101) rank_r
          !      Read(lwin,Err=101,End=101) beta_minus_r
          !      Read(lwin,Err=101,End=101) parity_even_r
          !      Read(lwin,Err=101,End=101) use_2bc_r
          !      if ((trim(label_r) /= g(ixterm)%label) .or. &
          !          k_r /= g(ixterm)%k .or. rank_r /= g(ixterm)%rank .or. &
          !          beta_minus_r /= merge(1,0,g(ixterm)%beta_minus) .or. &
          !          parity_even_r /= merge(1,0,g(ixterm)%parity_even) .or. &
          !          any(use_2bc_r /= g(ixterm)%use_2bc)) then
          !         iexit = 1
          !      end if
          !      ! Matrix elements
          !      If(Allocated(elem_r)) Deallocate(elem_r)
          !      Allocate(elem_r(nxy_r))
          !      Read(lwin,Err=101,End=101) elem_r
          !      if (nxy_r == fnxy) g(ixterm)%mat%elem = elem_r

          !      if (nxy12_r > 0) then
          !         If(Allocated(elem12_r)) Deallocate(elem12_r)
          !         Allocate(elem12_r(nxy12_r))
          !         Read(lwin,Err=101,End=101) elem12_r
          !         if (fnxy12 > 0 .and. nxy12_r == fnxy12) g(ixterm)%mat12%elem = elem12_r
          !      end if
          !   End Do
          !End If

       End If
    End Do

 98 Continue
    Close(lwin)
    If (iexit==1) Then
       Write(st,'(1x,a,a,a)')  " * The file '",trim(tbcfile),"' is incompatible with the current calculation!" ; call writelog(st)
       Write(st,'(1x,a)')      " * Starting calculation from scratch..."; call writelog(st)
       call writelog("")
    End If
    Return
    !
101 Continue
    iexit=1
    Close(lwin)
    Write(st,'(1x,a,a,a)')  " * The file '",trim(tbcfile),"' is corrupted!" ; call writelog(st)
    Write(st,'(1x,a,a8,a)') " * Problem occurs for key: ",trim(key) ; call writelog(st)
    Write(st,'(1x,a)')      " * Starting calculation from scratch..."; call writelog(st)
    call writelog("")
    !
  End Subroutine read_tbc

end module
