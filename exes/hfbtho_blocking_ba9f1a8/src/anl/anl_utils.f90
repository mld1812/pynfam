! anl module for interfacing with HFBTHO using mpi (and still using openmp)

! Expected use:
!    call anl_initialize()
!    call anl_readdata(computation_filename, residual_filename)
!       ... edit hfbtho settings ...
!       ... set up optimization parameters ...
!    call anl_printinfo
!    call anl_computenuclei(x) 
!    call anl_computeresiduals()
!    call anl_printresiduals(file_number, verbose)
!    call anl_finalize()


Module anl
  implicit none 
  integer,parameter :: anl_ENERGY_CLASS  = 1
  integer,parameter :: anl_RADIUS_CLASS  = 18
  integer,parameter :: anl_NOES_CLASS    = 12
  integer,parameter :: anl_POES_CLASS    = 13
  integer,parameter :: anl_MULTI_CLASS   = 101
  integer,parameter :: anl_MAX_RESIDUALS = 300
  integer,parameter :: anl_MAX_NUCLEI = 300
  integer,parameter :: anl_MAX_PARAMETERS = 14
  integer :: anl_rank
  integer :: anl_size
  integer :: anl_comm
  integer :: anl_next  ! kept on proc 0
  integer :: anl_nnuclei ! number of nuclei to run
  integer :: anl_nresiduals
  integer :: anl_nvars ! number of variables 
  integer :: anl_current ! current job running

  double precision :: anl_chi2 ! sum of squares of weighted residuals
  double precision :: anl_lower(anl_MAX_PARAMETERS),anl_upper(anl_MAX_PARAMETERS) ! scaling bounds
  double precision :: anl_xl(anl_MAX_PARAMETERS),anl_xu(anl_MAX_PARAMETERS)       ! unrelaxable bounds
  type anl_global_parameters
     integer     n00      ! number of oscilator shells
     double precision b0  ! basis parameter      
     integer     ilst     ! LST basis control (-1,0,1)
     integer     maxi     ! max iters
     character(len=30) skyrme   ! skyrme force code word ('SLY4', 'SIII', etc.)
     integer     kindhfb  ! 0=VAP,1=PAV w/o LN,2=PAV w/ LN
     integer     ippforce ! 0=no pairing, 1=mixed delta, 2=volume delta
     integer     icstr    ! use deformation constraint  (0=no,1=yes)
     integer     keypj    ! # of gauge angle points Ljp for the VAP and PAV
     integer     iproj    ! projection on different nucleus during VAP (0=no,1=yes)
     integer     npr1pj   ! N+DN=Number of neutrons requested if (iproj==1 && PAV)
     integer     npr2pj   ! Z+DZ=Number of protons requested if (iproj==1 && PAV)
     integer     icminput ! complete cmc calculate but not added
     integer     icrinput ! cranking rc calculated but not added
     logical     usej2
     logical     usecm    ! use center of mass correction
     logical     read_hel
     logical     write_hel
     logical     verbose
     character(Len=90) output_filename
     double precision epsi
  end type anl_global_parameters
  type(anl_global_parameters) anl_globals

  type anl_computation_type
     integer          :: ID
     integer          :: Z
     integer          :: N 
     double precision :: init_deformation
     double precision :: Q2
     integer          :: nk(5)
     integer          :: pk(5)
     character*(8)    :: class_name
     integer          :: shape
     logical          :: necessary ! only true if needed by some residual
     
     ! results
     double precision :: calc_energy
     double precision :: calc_rms
     double precision :: calc_noes
     double precision :: calc_poes
     double precision :: exp_energy
     double precision :: exp_rms
     double precision :: exp_noes
     double precision :: exp_poes
     integer          :: niters
     double precision :: time
     double precision :: si
     integer          :: ierr
     double precision :: deformation

  end type anl_computation_type

  type(anl_computation_type),target,dimension(anl_MAX_NUCLEI) :: anl_computation_list
  type anl_multi_item
     integer :: computation_index
     double precision :: multiplier
  end type anl_multi_item
  type anl_multi
     integer :: nitems
     type(anl_multi_item) :: item(3)
  end type anl_multi
  
  type anl_residual_type
     integer       :: ID
     integer       :: Z
     integer       :: N
     integer       :: res_class ! anl_ENERGY_CLASS, anl_RADIUS_CLASS, etc.
     integer       :: class ! anl_ENERGY_CLASS, anl_RADIUS_CLASS, etc.
     character(len=8) :: class_string
     character(len=50) :: hash
     double precision :: exp_value
     double precision :: calc_value
     double precision :: diff
     double precision :: weight
     double precision :: weighted_diff
     double precision :: exp_energy
     double precision :: exp_rms
     double precision :: exp_noes
     double precision :: exp_poes

     double precision :: calc_energy
     double precision :: calc_rms
     double precision :: calc_noes
     double precision :: calc_poes
     logical          :: is_multi
     integer          :: computation_index ! only valid if not multi
     type(anl_multi)  :: multi
     integer          :: niters
  end type anl_residual_type
  type(anl_residual_type),target,dimension(anl_MAX_RESIDUALS) :: anl_residual_list

  type anl_result_package_type
     integer             :: taskid
     integer             :: flag
     integer             :: niters
     double precision    :: energy
     double precision    :: rms
     double precision    :: noes
     double precision    :: poes
     double precision    :: deformation
     double precision    :: time
     double precision    :: si
     integer             :: ub  ! Used for mpi 
  end type  anl_result_package_type


  
  integer anl_result_package_mpitype

end module anl


subroutine anl_readdata(computation_file,residual_file)
  use anl
  use HFBTHO

  implicit none
  logical keepgoing
  integer i,ioflag
  integer lnk(5),lpk(5)
  integer z,n,nnuclei,nresiduals,nid
  integer pos1,pos2,pos3,pos4,pos5,loc1,loc2
  integer class
  double precision def,expV,weight,nq2
  character*(8)  class_name
  character(len=*) :: computation_file, residual_file
  character(len=800) :: line
  character(len=50) :: hash
  integer :: mf(3),cid(3)
  character :: s(3)
  type(anl_residual_type),pointer :: res

  ! Read in list of nuclei to run simulations on
  nnuclei=0
  keepgoing = .true.
  Open(901,file=computation_file,status='old',IOSTAT=ioflag)
  if (ioflag .ne. 0) then
     print *,'Error reading file ',computation_file
     stop
  end if
  do while (keepgoing) 
     
     read(901,'(A)',IOSTAT=ioflag) line
     if (ioflag .ne. 0 .or. line(1:3) .eq. 'END') then
        keepgoing = .false.
     else if (nnuclei .ge. anl_MAX_RESIDUALS) then
        print *,'ERROR: Number of residuals > MAX_NUCLEI. Adjust and recompile',anl_MAX_NUCLEI
        stop
     else if (line(1:1) .eq. '#') then
        continue
     else
        nnuclei = nnuclei + 1
        pos1 = index(line,'(')+1
        pos2 = index(line,')')-1
        pos3 = index(line,'(',BACK=.TRUE.)+1
        pos4 = index(line,')',BACK=.TRUE.)-1
        pos5 = pos4+3
        backspace(901)
        read(901,*) nid,z,n,nq2,def

        loc1=pos1
        do i=1,5
           if (i .lt. 5) then
              loc2 = loc1+index(line(loc1:),',')-2
           else
              loc2 = pos2
           end if
           if (loc2-loc1 .eq. 0) then
              read(line(loc1:loc2),FMT='(I1)') lpk(i)
           else if (loc2-loc1 .eq.1) then
              read(line(loc1:loc2),FMT='(I2)') lpk(i)
           else if (loc2-loc1 .eq.2) then
              read(line(loc1:loc2),FMT='(I3)') lpk(i)
           else if (loc2-loc1 .eq.3) then
              read(line(loc1:loc2),FMT='(I4)') lpk(i)
           else if (loc2-loc1 .eq.4) then
              read(line(loc1:loc2),FMT='(I5)') lpk(i)
           else if (loc2-loc1 .eq.5) then
              read(line(loc1:loc2),FMT='(I6)') lpk(i)
           else if (loc2-loc1 .ge.6) then
              read(line(loc2-6:loc2),FMT='(I7)') lpk(i)
           endif
           
           loc1 = loc2+2
        enddo

        loc1=pos3
        do i=1,5
           if (i .lt. 5) then
              loc2 = loc1+index(line(loc1:),',')-2
           else
              loc2 = pos4
           end if
           if (loc2-loc1 .eq. 0) then
              read(line(loc1:loc2),FMT='(I1)') lnk(i)
           else if (loc2-loc1 .eq.1) then
              read(line(loc1:loc2),FMT='(I2)') lnk(i)
           else if (loc2-loc1 .eq.2) then
              read(line(loc1:loc2),FMT='(I3)') lnk(i)
           else if (loc2-loc1 .eq.3) then
              read(line(loc1:loc2),FMT='(I4)') lnk(i)
           else if (loc2-loc1 .eq.4) then
              read(line(loc1:loc2),FMT='(I5)') lnk(i)
           else if (loc2-loc1 .eq.5) then
              read(line(loc1:loc2),FMT='(I6)') lnk(i)
           else if (loc2-loc1 .ge.6) then
              read(line(loc2-6:loc2),FMT='(I7)') lnk(i)
           endif
           loc1 = loc2+2
        enddo
        read(line(pos5:), *) class_name
        anl_computation_list(nnuclei)%ID=nid
        anl_computation_list(nnuclei)%Z=z
        anl_computation_list(nnuclei)%N=n
        anl_computation_list(nnuclei)%Q2=nq2
        anl_computation_list(nnuclei)%init_deformation = def
           

        do i=1,5
           anl_computation_list(nnuclei)%pk(i) = lpk(i)
           anl_computation_list(nnuclei)%nk(i) = lnk(i)
        enddo
        anl_computation_list(nnuclei)%class_name = class_name
        if (def .eq. 0.0d0) then
           anl_computation_list(nnuclei)%shape = 1
        else if (def .gt. 0.0d0) then 
           anl_computation_list(nnuclei)%shape = 2
        else 
           anl_computation_list(nnuclei)%shape = 3
        end if
        anl_computation_list%necessary = .False.
     end if
  enddo
  close(901)

  ! read in list of residuals to use in function evaluation
  keepgoing = .true.
  nresiduals = 0
  Open(900 ,file=residual_file,status='old',IOSTAT=ioflag)
  if (ioflag .ne. 0) then
     print *,'Error reading file ',residual_file
     stop
  end if

  do while (keepgoing)
     read(900,'(A)',IOSTAT=ioflag) line
     if (ioflag .ne. 0 .or. line(1:3) .eq. 'END') then
        keepgoing = .false.
     else if (nresiduals .ge. anl_MAX_RESIDUALS) then
        print *,'ERROR: Number of residuals > MAX_RESIDUALS',anl_MAX_RESIDUALS
        stop
     else if (line(1:1) .eq. '#') then
        continue
     else
        !backspace(900)
        read(line,*,IOSTAT=ioflag) nid,z,n,expV,weight,hash,class
        if (ioflag .ne. 0) then
           print *, 'Error reading chi2 file: ioflag=',ioflag,' nid=',nid
           stop
        end if

        nresiduals = nresiduals+1
        res => anl_residual_list(nresiduals)
        res%ID = nid
        res%Z = z
        res%N = n
        res%exp_value = expV
        res%weight = weight
        res%res_class = class
        if (class .eq. anl_ENERGY_CLASS) then
           res%class_string = 'ENERGY'
        else if (class .eq. anl_RADIUS_CLASS) then
           res%class_string = 'RMS'
        else if (class .eq. anl_NOES_CLASS) then
           res%class_string = 'NOES'
        else if (class .eq. anl_POES_CLASS) then
           res%class_string = 'POES'
        end if
        

        ! process hash string
        read(hash,400) mf(1),s(1),cid(1), mf(2),s(2),cid(2), mf(3),s(3),cid(3)
400     format (1X,I1,1X,A1,1X,I3,2X,I1,1X,A1,1X,I3,2X,I1,1X,A1,1X,I3)
        res%hash = hash
        res%multi%nitems = 0
        if (mf(1) .eq. 1 .and. mf(2).eq.0 .and. mf(3).eq.0) then
           res%is_multi = .false.
           res%computation_index = cid(1)
        else
           res%computation_index = -1
           res%is_multi = .true.
           do i=1,3
              if (mf(i) .ne. 0) then
                 res%multi%nitems = res%multi%nitems + 1
                 res%multi%item(i)%computation_index=cid(i)
                 if (s(i) .eq. 'm' .or. s(i) .eq. 'M') then
                    res%multi%item(i)%multiplier = -1.0d0*mf(i)
                 else
                    res%multi%item(i)%multiplier = 1.0d0*mf(i)
                 end if
              end if
           end do
        end if
     end if
  end do
  close(900)
  !print *,'nresiduals = ',nresiduals
  anl_nnuclei = nnuclei
  anl_nresiduals = nresiduals
  call check_residual_dependencies()

end subroutine anl_readdata

! for debugging
subroutine anl_printlists()
  use anl
  implicit none
  integer i
  type(anl_computation_type),pointer :: nuc
  type(anl_residual_type),pointer :: res
  print *,'Nuclei to compute:'

  do i=1,anl_nnuclei
     nuc => anl_computation_list(i)
     if (nuc%necessary) then 
        write (*,715) nuc%id,nuc%z,nuc%n,nuc%q2,nuc%init_deformation,nuc%nk,nuc%pk
     end if
  enddo

  !print *,'Nuclei to ignore:'
  !do i=1,anl_nnuclei
  !   nuc => anl_computation_list(i)
  !   if (.not.nuc%necessary) then 
  !      write (*,715) nuc%id,nuc%z,nuc%n,nuc%q2,nuc%init_deformation,nuc%pk,nuc%nk
  !   end if
  !enddo
715     format(I3,3X,'(',I3,',',I3,') q2=',F16.12,' beta=',F4.1,' (',5I3,')(',5I3,')')


  print *
  print *,'Residuals to compute:',anl_nresiduals
  do i=1,anl_nresiduals
     res => anl_residual_list(i)
     write (*,716) res%id,res%z,res%n,res%exp_value,res%weight,res%res_class
  end do
716 format(I3,3X,'(',I3,',',I3,') exp=',F20.14,' weight=',F5.2,' class=',I3)

end subroutine anl_printlists



! set up mpi structure and default parameters
subroutine anl_initialize(comm)
  use anl
  use HFBTHO_solver
!  use mpi
  implicit none
#if(USE_MPI==1)
  include 'mpif.h'
#endif
  integer :: comm

#if(USE_MPI==1)
  integer :: ierr,i
  integer :: icom_err
  integer :: mpibase, disp(11)
  integer :: blocklen(11),types(11)
  type(anl_result_package_type) :: package
  if (comm .le. 0) then
     call MPI_Init(ierr)
     comm = MPI_COMM_WORLD
     call MPI_Comm_rank(comm,anl_rank,ierr)
     call MPI_Comm_size(comm,anl_size,ierr)
  endif
  anl_comm = comm
  call MPI_GET_ADDRESS(package%taskid,disp(1),icom_err) 
  call MPI_GET_ADDRESS(package%flag,disp(2),icom_err) 
  call MPI_GET_ADDRESS(package%niters,disp(3),icom_err) 
  call MPI_GET_ADDRESS(package%energy,disp(4),icom_err) 
  call MPI_GET_ADDRESS(package%rms,disp(5),icom_err) 
  call MPI_GET_ADDRESS(package%noes,disp(6),icom_err) 
  call MPI_GET_ADDRESS(package%poes,disp(7),icom_err) 
  call MPI_GET_ADDRESS(package%deformation,disp(8),icom_err) 
  call MPI_GET_ADDRESS(package%time,disp(9),icom_err)
  call MPI_GET_ADDRESS(package%si,disp(10),icom_err)
  call MPI_GET_ADDRESS(package%ub,disp(11),icom_err) 
  mpibase = disp(1) 
  do i=1,11
     disp(i) = disp(i) - mpibase
     blocklen(i)=1
  enddo
  types(1) = MPI_INTEGER 
  types(2) = MPI_INTEGER 
  types(3) = MPI_INTEGER 
  types(4) = MPI_DOUBLE_PRECISION 
  types(5) = MPI_DOUBLE_PRECISION 
  types(6) = MPI_DOUBLE_PRECISION 
  types(7) = MPI_DOUBLE_PRECISION 
  types(8) = MPI_DOUBLE_PRECISION 
  types(9) = MPI_DOUBLE_PRECISION
  types(10) = MPI_DOUBLE_PRECISION
  types(11) = MPI_UB
  call MPI_TYPE_STRUCT(11,blocklen,disp,types,anl_result_package_mpitype,icom_err) 
  call MPI_TYPE_COMMIT(anl_result_package_mpitype,icom_err) 
#else
  anl_rank = 0
  anl_size = 1
#endif  
  anl_nvars = 12

  anl_globals%verbose = .false.
  anl_globals%usej2 = .false.
  anl_globals%usecm = .false.
  anl_globals%epsi =  1.0e-5
  anl_globals%n00 = 20
  anl_globals%b0 = -2.319
  anl_globals%ilst = 0
  anl_globals%maxi = 200
  anl_globals%skyrme = 'FITS'
  anl_globals%kindhfb = -1
  anl_globals%ippforce = 1
  anl_globals%icstr = 0
  anl_globals%keypj = 1
  anl_globals%iproj = 0
  anl_globals%npr1pj = 0
  anl_globals%npr2pj = 0
  anl_globals%ICMinput = 1
  anl_globals%ICRinput = 0
  anl_globals%write_hel = .false.
  anl_globals%read_hel = .false.

  anl_lower(1)=     0.14d0; anl_upper(1)=    0.18d0 !rho
  anl_lower(2)=    -17.0d0; anl_upper(2)=   -15.0d0 !e
  anl_lower(3)=    170.0d0; anl_upper(3)=   270.0d0 !k
  anl_lower(4)=     27.0d0; anl_upper(4)=    37.0d0 !ass
  anl_lower(5)=     30.0d0; anl_upper(5)=    70.0d0 !lass
  anl_lower(6)=      0.8d0; anl_upper(6)=     2.0d0 !smass
  anl_lower(7)=   -100.0d0; anl_upper(7)=   -40.0d0 !cdrho0
  anl_lower(8)=   -100.0d0; anl_upper(8)=   100.0d0 !cdrho1
  anl_lower(9)=   -350.0d0; anl_upper(9)=  -150.0d0 !v0(1)
  anl_lower(10)=  -350.0d0; anl_upper(10)= -150.0d0 !v0(2)
  anl_lower(11) = -120.0d0; anl_upper(11)=  -50.0d0 !CdJ0
  anl_lower(12) = -100.0d0; anl_upper(12)=   50.0d0 !CdJ1
  anl_lower(13) = -100.0d0; anl_upper(13)=  100.0d0 !CJ(0)
  anl_lower(14) = -100.0d0; anl_upper(14)=  100.0d0 !CJ(1)

end subroutine anl_initialize

!placeholder
subroutine anl_finalize
#if(USE_MPI==1)  
  !use mpi
  use anl
  implicit none
  include 'mpif.h'
  integer ierr
  call mpi_barrier(anl_comm,ierr)
#endif
end subroutine anl_finalize



subroutine manageworkers
#if(USE_MPI==1)
  use anl
  !use mpi
  implicit none
  include 'mpif.h'

  integer status(100)
  integer izero,minus1
  integer taskid
  integer worker
  
  integer next
  integer tag
  integer ierr
  integer i,finishedtasks,checkedin,necessary_tasks
  integer N,Z
  type(anl_result_package_type) :: incoming
  type(anl_computation_type),pointer :: comp
  logical found_next
  next =0
  izero = 0
  minus1=-1
  tag=1
  finishedtasks=0
  checkedin=0
  necessary_tasks = 0
  do i=1,anl_nnuclei
     if (anl_computation_list(i)%necessary) then
        necessary_tasks = necessary_tasks + 1
     end if
  end do
  do while ((finishedtasks .lt. necessary_tasks) .OR. (checkedin .lt. anl_size-1))
     call MPI_Recv(incoming%taskid,1,anl_result_package_mpitype,mpi_any_source,tag,anl_comm,status,ierr)
     if (ierr .ne. 0) then
        print *,'ierr0=',ierr
     end if
     worker=status(mpi_source)
     taskid = incoming%taskid
     if (taskid .eq. 0) then
        checkedin=checkedin+1
     else
        if (incoming%flag .eq. 0) then
           Z = anl_computation_list(taskid)%Z
           N = anl_computation_list(taskid)%N

           !write (*,301) Z,N,incoming%energy,incoming%rms,incoming%noes,incoming%poes,incoming%niters
301 format('Rcvd: (',I4,',',I4,') E=',f20.13,' r=',f20.13,' no=',f20.13,' po=',f20.13,' it=',i4)
           anl_computation_list(taskid)%ierr = incoming%flag
           anl_computation_list(taskid)%niters = incoming%niters
           anl_computation_list(taskid)%calc_energy = incoming%energy
           anl_computation_list(taskid)%calc_rms = incoming%rms
           anl_computation_list(taskid)%calc_noes = incoming%noes
           anl_computation_list(taskid)%calc_poes = incoming%poes
           anl_computation_list(taskid)%deformation = incoming%deformation
           anl_computation_list(taskid)%time = incoming%time
           anl_computation_list(taskid)%si = incoming%si
        else if (incoming%flag .lt. 0) then
           print *,'Manager: error computing task #',taskid
        else
           print *,'Manager: unnecessary task #',taskid
        end if

        
        finishedtasks = finishedtasks+1
     endif
     ! find next necessary task
     found_next = .False.
     do while (.not. found_next .and. next .lt. anl_nnuclei)
        next = next + 1
        comp => anl_computation_list(next)
        if (comp%necessary) then
           found_next = .True.
        end if
     end do
     if (found_next) then
        call MPI_Send(next,1,MPI_INTEGER,worker,tag,anl_comm,ierr)
        if (ierr .ne. 0) then
           print *,'ierr1=',ierr
        end if
     else
        Call mpi_send(minus1,1,MPI_INTEGER,worker,tag,anl_comm,ierr)
        if (ierr .ne. 0) then
           print *,'ierr2=',ierr
        end if
     end if
  end do
#endif
end subroutine manageworkers



subroutine anl_computenuclei(x)
  use anl
  use HFBTHO
  use HFBTHO_utilities
  use HFBTHO_solver
  implicit none

  double precision x(anl_MAX_PARAMETERS)
  double precision cutoff_tol_init
  double precision starttime,endtime
  integer tag,izero,taskid
  integer feval_flag,info,n,i
  logical idle
  logical keepgoing,singleproc,ismaster
  type(anl_result_package_type) :: outgoing
  type(anl_computation_type),pointer :: computation

  tag=1
  izero=0
  cutoff_tol_init=1.0d-6
  call setparameters(x)
  call printparameters(x)
  
  singleproc=(anl_size .eq. 1)
  ismaster = (anl_rank .eq. 0)
  if (singleproc .or. anl_rank .gt. 0) then
    keepgoing = .true.
    n = 1
    idle = .true.
    feval_flag = 0
    outgoing%taskid=0
    outgoing%flag=0
    outgoing%energy=0.0d0
    outgoing%rms=0.0d0
    outgoing%noes = 0.0d0
    outgoing%poes = 0.0d0
    idle = .false.
    call anl_checkin(outgoing,taskid)
    do while (taskid .gt. 0) 
       computation => anl_computation_list(taskid)
       npr_INI(1)=computation%N
       npr_INI(2)=computation%Z
       do i=1,5
          nkblo_INI(1,i)=computation%nk(i)
          nkblo_INI(2,i)=computation%pk(i)
       enddo

       inin_INI = 1!computation%shape
       if (computation%init_deformation .gt. 0) then
          inin_INI = 2
       end if
       !if (anl_globals%read_hel) then
       !   inin_INI = -inin_INI!computation%shape
       !end if
       q_INI = computation%init_deformation
       b2_0 = computation%init_deformation
       b4_0 = 0.0
       !if (computation%init_deformation .gt. 0.25 .and. computation%init_deformation .lt. 0.35) then
       !   b2_0 = q_INI!computation%init_deformation
       !   b4_0 = q_INI!computation%init_deformation
       !end if
       lambda_active(2) = -1
       expectation_values(2) = computation%Q2
       call CPU_Time(starttime)
       Call HFBTHO_DFT_SOLVER
       call CPU_Time(endtime)

       ierror_flag=0
       info = ierror_flag
       if (ierror_flag .ne. 0) then
          print *,ierror_info(ierror_flag)
       endif
       ! TODO error handling

       outgoing%taskid = taskid
       outgoing%energy = eresu(1)
       outgoing%rms = eresu(18)
       outgoing%noes = eresu(12)+ala2(1)
       outgoing%poes= eresu(13)+ala2(2)
       outgoing%deformation = eresu(3)
       outgoing%niters = iiter
       outgoing%flag = info
       outgoing%time = endtime-starttime
       outgoing%si = si

       
       ! Send result to master and check for next task
       call anl_checkin(outgoing,taskid)
    end do
 else
    call manageworkers()
 endif

#if(USE_MPI==1)
 call mpi_barrier(anl_comm,i)
#endif

end subroutine anl_computenuclei


subroutine check_residual_dependencies()
  use anl
  implicit none   
  type(anl_residual_type),pointer :: res
  type(anl_computation_type),pointer :: comp
  integer i,j
  if (anl_rank .eq. 0) then 
     do i=1,anl_nresiduals
        res => anl_residual_list(i)
        
        if (.not.res%is_multi) then
           comp => anl_computation_list(res%computation_index)
           comp%necessary = .True.
           comp%exp_energy = res%exp_energy
           comp%exp_rms    = res%exp_rms
           comp%exp_noes   = res%exp_noes
           comp%exp_poes   = res%exp_poes
        else
           do j=1,res%multi%nitems
              comp => anl_computation_list(res%multi%item(j)%computation_index)
              comp%necessary = .True.
           end do
        end if
     end do
  end if
end subroutine check_residual_dependencies

subroutine printme(i)
  use anl
  implicit none
  integer i
  type(anl_residual_type),pointer :: res
  res => anl_residual_list(i)
  print *,res%Z,res%N,res%calc_energy,res%exp_energy
end subroutine printme

subroutine anl_computeresiduals()
  use anl
  implicit none
  type(anl_residual_type),pointer :: res
  type(anl_computation_type),pointer :: comp
  integer i,j
  anl_chi2 = 0.0
  if (anl_rank .eq. 0) then
     do i=1,anl_nresiduals
        res => anl_residual_list(i)
        if (res%is_multi) then
           res%calc_value=0.0
           do j=1,res%multi%nitems
              comp => anl_computation_list(res%multi%item(j)%computation_index)
              print *,'spe[',comp%ID,']: ',j,', ',comp%calc_energy ! DEBUG
              res%calc_value = res%calc_value + res%multi%item(j)%multiplier*comp%calc_energy
           end do
        else 
           comp => anl_computation_list(res%computation_index)
           if (res%res_class .eq. anl_ENERGY_CLASS) then
              res%calc_value = comp%calc_energy
           else if (res%res_class .eq. anl_RADIUS_CLASS) then
              res%calc_value = comp%calc_rms
           else if (res%res_class .eq. anl_NOES_CLASS) then
              res%calc_value = comp%calc_noes
           else if (res%res_class .eq. anl_POES_CLASS) then
           ! noes or poes
              res%calc_value = comp%calc_poes
           end if
        end if
        res%diff = res%calc_value - res%exp_value
        res%weighted_diff = res%diff / res%weight
        anl_chi2 = anl_chi2 + res%weighted_diff * res%weighted_diff

     end do
  end if
end subroutine anl_computeresiduals


subroutine anl_printresiduals(fd,verbose)
 use anl
 use,intrinsic :: iso_fortran_env, only : output_unit
 implicit none
 integer fd
 integer i
 logical verbose
 type(anl_residual_type),pointer :: res
 type(anl_computation_type),pointer :: comp
 double precision chi2
 chi2=0
 if (fd .lt. 0) then
    fd = output_unit
 end if
 if (anl_rank .eq. 0) then
    if (verbose) then
       write(fd,*) 'Computed Nuclei:'
       do i=1,anl_nnuclei
          comp => anl_computation_list(i)
          write(fd,100) comp%ID,comp%Z,comp%N,comp%calc_energy,comp%calc_rms,comp%calc_noes,comp%calc_poes,comp%niters,comp%time,comp%si
       end do
    end if
    if (verbose) then
       write(fd,*) 'Weighted difference vector (computed - experimental)/weight'
    end if
    do i=1,anl_nresiduals
       res => anl_residual_list(i)
       if (verbose) then
          write(fd,101) res%Z,res%N,res%class_string,res%calc_value,res%exp_value,res%weight,res%weighted_diff
          chi2 = chi2+res%weighted_diff**2
       else
          write(fd,*) res%weighted_diff
       end if
    end do
    if (verbose) then
       write(fd,102) chi2
    end if
 end if
100 format(I3,I4,I4,F15.8,F11.8,F11.8,F11.8,I4,F12.1,F14.10)
101 format(I3,I4,' ',A8,F24.16,F24.16,F5.2,F24.16)
102 format('||f||^2 = ',f24.16)

end subroutine anl_printresiduals

subroutine setparameters(x)
  use anl
  use HFBTHO
  implicit none
  double precision x(anl_MAX_PARAMETERS)
  character(len=30) :: forcename
  integer i,noforces

  noforces = 0
  forcename = 'FITS'
  call read_unedf_namelist(forcename,noforces)


  if (anl_globals%verbose) then
     n00_INI = anl_globals%n00
  else
     n00_INI = -abs(anl_globals%n00)
  end if
  b0_INI=anl_globals%b0
  ilst_INI=anl_globals%ilst
  max_iter_INI=anl_globals%maxi
  skyrme_INI=anl_globals%skyrme
  kindhfb_INI=anl_globals%kindhfb
  keypj_INI=anl_globals%keypj
  iproj_INI=anl_globals%iproj
  npr1pj_INI=anl_globals%npr1pj
  npr2pj_INI=anl_globals%npr2pj
  epsi_INI=anl_globals%epsi
  write_hel = anl_globals%write_hel
  Add_Pairing_INI=.False.
  set_pairing = .False.
  icou_INI=2
  Parity_INI=.True.
  pwi_INI = 60.0
  cpv1_INI = 0.5
  ngh_INI = 40
  ngl_INI = 40
  nleg_INI = 80
  basis_HFODD_INI = .false.
  nstate_INI = 500
  Print_HFBTHO_Namelist_INI = .False.
  IDEBUG_INI = 0
  switch_on_temperature = .False.
  temper = 0
  !use_Namelist = .False.
  use_INM = .True.
  do i=1,8
     lambda_values(i) = i
     lambda_active(i) = 0
     expectation_values(i) = 0.0
  end do

  !thoout_dat = anl_globals%output_filename
  
  VMASS_NM = 1.249d0  
  VMASS_NM = 1.24983857423227d0 
  CExPar = 1.0 ! fixed
  RHO_NM = x(1)
  E_NM = x(2)
  K_NM = x(3)
  ASS_NM = x(4)
  LASS_NM = x(5)
  SMASS_NM = x(6)
  CRDR(0) =x(7)  !Cdrho0 = x(7)
  CRDR(1) =x(8)  !Cdrho1 = x(8)
  CPV0(0) = x(9)
  CPV0(1) = x(10)
  CRDJ(0) = x(11) !CdJ0 = x(11)
  CRDJ(1) = x(12) !CdJ1 = x(12)

  use_cm_cor = anl_globals%usecm
  if (anl_globals%usej2) then
     CJ(0)=x(13)
     CJ(1)=x(14)
     use_j2terms = .true.
  else
     CJ(0) = 0.0
     CJ(1) = 0.0
     use_j2terms = .false.
  endif

end subroutine setparameters

subroutine anl_printinfo()
  use anl
#if(USE_OPENMP==1)
    Use omp_lib
#endif
  implicit none
  if (anl_rank .eq. 0) then
#if(USE_OPENMP==1)
     print 106,'num_threads',omp_get_max_threads()
#else
     print *, 'openmp disabled, num_threads'
#endif
     print 106,'n00',anl_globals%n00
     print 106,'ilst',anl_globals%ilst
     print 106,'maxi',anl_globals%maxi
     print 106,'kindhfb',anl_globals%kindhfb
     print 106,'ippforce',anl_globals%ippforce
     print 106,'icstr',anl_globals%icstr
     print 106,'keypj',anl_globals%keypj
     print 106,'iproj',anl_globals%iproj
     print 106,'icminput',anl_globals%icminput
     print 106,'icrinput',anl_globals%icrinput
     print 107,'usecm',anl_globals%usecm
     print 107,'usej2',anl_globals%usej2
     print 107,'read_hel',anl_globals%read_hel
     print 107,'write_hel',anl_globals%write_hel
     print 108,'skyrme',anl_globals%skyrme
     print 109,'epsi',anl_globals%epsi
     print 109,'b0',anl_globals%b0
     print *,'Ranges:'
     print 110,anl_lower(1),'rho    ',anl_upper(1)
     print 110,anl_lower(2),'e      ',anl_upper(2)
     print 110,anl_lower(3),'k      ',anl_upper(3)
     print 110,anl_lower(4),'ass    ',anl_upper(4)
     print 110,anl_lower(5),'lass   ',anl_upper(5)
     print 110,anl_lower(6),'smass  ',anl_upper(6)
     print 110,anl_lower(7),'CRDR(0)',anl_upper(7)
     print 110,anl_lower(8),'CRDR(1)',anl_upper(8)
     print 110,anl_lower(9),'CPV0(0)',anl_upper(9)
     print 110,anl_lower(10),'CPV0(1)',anl_upper(10)
     print 110,anl_lower(11),'CRDJ(0)',anl_upper(11)
     print 110,anl_lower(12),'CRDJ(1)',anl_upper(12)
     print 110,anl_lower(13),'CJ(0)  ',anl_upper(13)
     print 110,anl_lower(14),'CJ(1)  ',anl_upper(14)
  end if
106 format(a9,'=',i5)
107 format(a9,'=    ',l1)
108 format(a9,'=    ',a4)
109 format(a9,'=',F20.15)
110 format(f8.2,' < ',a7,' < ', f8.2)
end subroutine anl_printinfo

subroutine printparameters(x)
  use anl
  use HFBTHO
  implicit none
  double precision x(anl_MAX_PARAMETERS)
  integer i
  
  if (anl_rank .eq. 0) then
    print *,'Parameters:'
    print *,'rho   =',RHO_NM
    print *,'e     =',E_NM
    print *,'k     =',K_NM
    print *,'ass   =',ASS_NM
    print *,'lass  =',LASS_NM
    print *,'smass =',SMASS_NM
    print *,'CRDR(0)=',CRDR(0)
    print *,'CRDR(1)=',CRDR(1)
    print *,'CPV0(0)    =',CPV0(0)
    print *,'CPV0(1)    =',CPV0(1)
    print *,'CRDJ(0) = ',CRDJ(0)
    print *,'CRDJ(1) = ',CRDJ(1)
    print *,'CJ(0) =',CJ(0)
    print *,'CJ(1) =',CJ(1)

    print *,'X-transformed:'
    do i=1,anl_nvars
       print 105,i,(x(i)-anl_lower(i))/(anl_upper(i)-anl_lower(i))
    enddo
105 format('x[',i2,']=',f20.13)
    !print 106,'Z',anl_globals%npr
    !print 106,'N',anl_globals%npr2
    !print 109,'q',anl_globals%q
  endif

end subroutine printparameters

subroutine anl_checkin(new_result, next_id)
  use anl
  implicit none
#if(USE_MPI==1)  
  !use mpi
  include 'mpif.h'
  integer mpierr
  integer status(100)
#endif

  type(anl_result_package_type) :: new_result
  integer next_id
  integer next
  type(anl_computation_type),pointer :: comp
  integer taskid,tag
  logical found_next

  tag =1
  if (anl_size .eq. 1) then
     taskid = new_result%taskid
     if (taskid .ne. 0) then
        comp => anl_computation_list(taskid)
        comp%calc_energy = new_result%energy
        comp%calc_rms = new_result%rms
        comp%calc_noes = new_result%noes
        comp%calc_poes = new_result%poes
        comp%niters = new_result%niters
        comp%time = new_result%time
        comp%si = new_result%si
        comp%ierr = new_result%flag
     end if

     found_next = .False.
     next = new_result%taskid + 1
     do while (.not. found_next .and. next .le. anl_nnuclei)
        comp =>anl_computation_list(next)
        if (comp%necessary) then
           found_next = .True.
        else
           next = next + 1
        end if
     end do
     if (.not. found_next) then
        next = -1
     end if
  else
     ! Only used when using mpi
#if(USE_MPI==1)     
     call MPI_Send(new_result%taskid,1,anl_result_package_mpitype,0,tag,anl_comm,mpierr)
     call MPI_Recv(next,1,MPI_INTEGER,mpi_any_source,mpi_any_tag,anl_comm,status,mpierr)
#endif

  end if
  next_id = next
end subroutine anl_checkin



subroutine anl_set_unedf0
  ! Options used when computing unedf0 from sly4
  use anl
  use HFBTHO
  use_cm_cor = .True.
  use_j2terms = .False.
  anl_globals%usej2 = .False.
  anl_nvars = 12
end subroutine anl_set_unedf0

subroutine anl_set_unedf1
  ! Options used when computing unedf1 from unedf0
  use anl
  use HFBTHO
  use_cm_cor = .False.
  use_j2terms = .False.
  anl_globals%usej2 = .False.
  anl_nvars = 12
end subroutine anl_set_unedf1

subroutine anl_set_unedf2
  ! Options used when computing unedf2 from unedf1
  use anl
  use HFBTHO
  use_cm_cor = .False.
  use_j2terms = .True.
  anl_globals%usej2 = .True.
  anl_nvars = 14
end subroutine anl_set_unedf2
  
     
  

  
