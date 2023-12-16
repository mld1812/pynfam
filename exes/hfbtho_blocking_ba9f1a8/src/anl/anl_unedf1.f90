! compute_energy computes the residuals of a set of nuclei at a given parameter set.
! Input files:
!      UNEDF_chi2.dat   --  list of residuals to compute, with experimental data to compare with
!      UNEDF_hfbtho.dat --  list of necessary information to compute each nuclei. 
!                           compute_energy will only compute the nuclei that it deems necessary after
!                           reading the list of residuals to compute
!      x.in             --  list of parameters:
!      

Program compute_energy
  use anl
  use HFBTHO_utilities
  use HFBTHO
  use HFBTHO_solver
  Implicit None
  logical :: singleproc,verbose
  character(len=50) :: comp_filename
  character(len=50) :: chi2_filename
  double precision :: x(15)
  integer :: output,i,ioflag,comm
  
  comm=0
  call anl_initialize(comm)
  anl_globals%write_hel = .False.
  anl_globals%read_hel = .False.
  anl_globals%verbose = .False.
  anl_globals%output_filename ="thoout.dat"
  anl_globals%maxi = 100
  anl_globals%usecm = .False.
  anl_globals%usej2 = .False.
  anl_globals%epsi = 1.0d-5

  Open(801 ,file='x.in',status='old',IOSTAT=ioflag)   ! open the input file  
  if (ioflag .ne. 0) then
     print *,'Error reading file x.in, perhaps it does not exist'
     stop
  end if
  do i=1,anl_nvars
     read(801,*) x(i)
  end do
  close(801)
  comp_filename =  'UNEDF2Argonne_hfbtho.dat'
  chi2_filename =  'UNEDF1_chi2.dat'
  call anl_readdata(comp_filename, chi2_filename)

  !if (anl_rank .eq. 0) then
  !   call anl_printlists()
  !end if

  call anl_printinfo()
  call anl_computenuclei(x)
  call anl_computeresiduals()

  !print to screen
  output = 1
  verbose = .True.
  call anl_printresiduals(output,verbose)

  !write to file
  open(802,file='fval.out',status='unknown')
  output = 802
  verbose=.False.
  call anl_printresiduals(output,verbose)
  close(802)
     
  call anl_finalize()

End Program compute_energy
