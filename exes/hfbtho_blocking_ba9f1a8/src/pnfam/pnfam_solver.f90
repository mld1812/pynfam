!------------------------------------------------------------------------------
! pnfam_solver.f90
!
! M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
!------------------------------------------------------------------------------
module pnfam_solver
   use logger
   use blockmatrix_type
   use external_field_type
   use density_set_type
   implicit none
   integer, parameter, private :: dp = kind(1d0)
   real(dp), parameter :: pi = 3.141592653589793238462643_dp
   
   type(blockmatrix), private :: greenXre, greenXim, greenYre, greenYim, f20, f02,              &
      h20re, h20im, h02re, h02im, rex, imx, rey, imy, rerho_np, imrho_np, rerho_pn, imrho_pn,   &
      rekp, imkp, rekm, imkm, redp, imdp, redm, imdm, reh_np, imh_np, reh_pn, imh_pn
   type(blockmatrix), dimension(:), allocatable, private :: g20, g02
   type(density_set), private :: rhox
   
   ! These blocks are introduced only at finite temperature
   type(blockmatrix), private :: greenPre, greenPim, greenQre, greenQim, f11, f11t, &
         h11re, h11im, h11tre, h11tim, rep, imp, req, imq, &
         ft_fbfa_p, ft_fbfa_q, ft_1fafb_x, ft_1fafb_y
   type(blockmatrix), dimension(:), allocatable, private :: g11, g11t
   
   integer, private :: nxy
   
contains
   
   !---------------------------------------------------------------------------
   ! Subroutine for solving the FAM equations for one value of energy. Before
   ! calling this the first time, the module data structures must be initialized
   ! with init_pnfam_solver.
   !---------------------------------------------------------------------------
   subroutine pnfam_solve(omega, width, Ep, En, vp, up, vn, un, max_iter, &
      convergence_epsilon, skip_residual_interaction, bminus, str, iter_used, cstr)
      use constants, only : iu
      use broyden_mixer
      use hamiltonian
      use hfbtho_basis, only : ft_active, ft_temp
      implicit none
      real(dp), intent(in) :: omega, width, Ep(:), En(:), convergence_epsilon
      integer, intent(in) :: max_iter
      type(blockmatrix), intent(in) :: vp, vn, up, un
      logical, intent(in) :: skip_residual_interaction, bminus
      integer, intent(out) :: iter_used
      real(dp),    dimension(:), intent(out) :: str
      complex(dp), dimension(:), intent(out), optional :: cstr ! complex strength
      
      ! Finite temperature
      real(dp)    :: ft_real_factor
      complex(dp) :: ft_cmplx_factor
      
      integer :: iter, ixterm
      
      ! Initial guess for X and Y
      ! We keep the real and complex part of X and Y separate to allow
      ! optimizing various matrix operations, where often only one of the
      ! matrices is complex and lifting the real ones to be complex is a waste
      ! of computing time
      ReX%elem = 0; ImX%elem = 0 ; ReY%elem = 0 ; ImY%elem = 0
      
      ! Finite temperature
      if (ft_active) then
         ReP%elem = 0; ImP%elem = 0 ; ReQ%elem = 0 ; ImQ%elem = 0
      end if
      
      ! Compute the denominators in the FAM equations with minus sign attached,
      ! e.g. -1/(E_p + E_n - \omega)
      if (bminus) then
         call greenx(Ep, En, -cmplx(omega,width,dp), greenXre, greenXim)
         call greenx(Ep, En,  cmplx(omega,width,dp), greenYre, greenYim)
         
         ! Finite temperature
         if (ft_active) then
            call greenx(Ep, -En, -cmplx(omega,width,dp), greenPre, greenPim)
            call greenx(Ep, -En,  cmplx(omega,width,dp), greenQre, greenQim)            
         end if
      else
         call greenx(En, Ep, -cmplx(omega,width,dp), greenXre, greenXim)
         call greenx(En, Ep,  cmplx(omega,width,dp), greenYre, greenYim)
         
         ! Finite temperature
         if (ft_active) then
            call greenx(En, -Ep, -cmplx(omega,width,dp), greenPre, greenPim)
            call greenx(En, -Ep,  cmplx(omega,width,dp), greenQre, greenQim)            
         end if
      end if
   
      ! The iteration loop
      si = 1
      iter_used = -1  ! in case of no convergence, this value will remain
      
      do iter=0,max_iter
         if (skip_residual_interaction) then
            h20re%elem = f20%elem ; h20im%elem = 0  ! only the external field
            h02re%elem = f02%elem ; h02im%elem = 0
            
            ! Finite temperature
            if (ft_active) then
               h11re%elem  = f11%elem  ; h11im%elem  = 0
               h11tre%elem = f11t%elem ; h11tim%elem = 0
            end if
         else
         ! Compute the perturbed densities
         if (bminus) then
            ! rho_np = -Un X(-)' Vp' + Vn Y(-)' Up'
            call triprod('n', un, 't', rex, 't', vp, -1d0, 0d0, rerho_np)
            call triprod('n', un, 't', imx, 't', vp, -1d0, 0d0, imrho_np)
            call triprod('n', vn, 't', rey, 't', up, 1d0, 1d0, rerho_np)
            call triprod('n', vn, 't', imy, 't', up, 1d0, 1d0, imrho_np)
            
            ! rho_pn = Up X(-) Vn' - Vp Y(-) Un'
            call triprod('n', up, 'n', rex, 't', vn, 1d0, 0d0, rerho_pn)
            call triprod('n', up, 'n', imx, 't', vn, 1d0, 0d0, imrho_pn)
            call triprod('n', vp, 'n', rey, 't', un, -1d0, 1d0, rerho_pn)
            call triprod('n', vp, 'n', imy, 't', un, -1d0, 1d0, imrho_pn)
            
            ! kappa(+)_np = -un X(-)' up' + vn Y(-)' vp'
            call triprod('n', un, 't', rex, 't', up, -1d0, 0d0, rekp)
            call triprod('n', un, 't', imx, 't', up, -1d0, 0d0, imkp)
            call triprod('n', vn, 't', rey, 't', vp, 1d0, 1d0, rekp)
            call triprod('n', vn, 't', imy, 't', vp, 1d0, 1d0, imkp)
            
            ! kappa(-)_np = vn X(-)*' vp' - un Y(-)'* up'
            call triprod('n', vn, 't', rex, 't', vp, 1d0, 0d0, rekm)
            call triprod('n', vn, 't', imx, 't', vp, -1d0, 0d0, imkm)
            call triprod('n', un, 't', rey, 't', up, -1d0, 1d0, rekm)
            call triprod('n', un, 't', imy, 't', up, 1d0, 1d0, imkm)
            
            ! Finite temperature
            if (ft_active) then
               ! rho_pn += up P(-) un' - vp Q(-) vn'
               call triprod('n', up, 'n', rep, 't', un, 1d0, 1d0, rerho_pn)
               call triprod('n', up, 'n', imp, 't', un, 1d0, 1d0, imrho_pn)
               call triprod('n', vp, 'n', req, 't', vn,-1d0, 1d0, rerho_pn)
               call triprod('n', vp, 'n', imq, 't', vn,-1d0, 1d0, imrho_pn)
               
               ! rho_np += Un Q(-)' Up' - Vn P(-)' Vp'
               call triprod('n', un, 't', req, 't', up, 1d0, 1d0, rerho_np)
               call triprod('n', un, 't', imq, 't', up, 1d0, 1d0, imrho_np)
               call triprod('n', vn, 't', rep, 't', vp,-1d0, 1d0, rerho_np)
               call triprod('n', vn, 't', imp, 't', vp,-1d0, 1d0, imrho_np)
               
               ! kappa(+)_np += un Q(-)' vp' - vn P(-)' up'
               call triprod('n', un, 't', req, 't', vp, 1d0, 1d0, rekp)
               call triprod('n', un, 't', imq, 't', vp, 1d0, 1d0, imkp)
               call triprod('n', vn, 't', rep, 't', up,-1d0, 1d0, rekp)
               call triprod('n', vn, 't', imp, 't', up,-1d0, 1d0, imkp)
               
               ! kappa(-)_np += -vn Q(-)'* up' + un P(-)'* vp'
               call triprod('n', vn, 't', req, 't', up,-1d0, 1d0, rekm)
               call triprod('n', vn, 't', imq, 't', up, 1d0, 1d0, imkm)
               call triprod('n', un, 't', rep, 't', vp, 1d0, 1d0, rekm)
               call triprod('n', un, 't', imp, 't', vp,-1d0, 1d0, imkm)
            end if
            
         else
            ! rho_np = Un X(+) Vp' - Vn Y(+) Up'
            call triprod('n', un, 'n', rex, 't', vp, 1d0, 0d0, rerho_np)
            call triprod('n', un, 'n', imx, 't', vp, 1d0, 0d0, imrho_np)
            call triprod('n', vn, 'n', rey, 't', up, -1d0, 1d0, rerho_np)
            call triprod('n', vn, 'n', imy, 't', up, -1d0, 1d0, imrho_np)
            
            ! rho_pn = -Up X(+)' Vn' + Vp Y(+)' Un'
            call triprod('n', up, 't', rex, 't', vn, -1d0, 0d0, rerho_pn)
            call triprod('n', up, 't', imx, 't', vn, -1d0, 0d0, imrho_pn)
            call triprod('n', vp, 't', rey, 't', un, 1d0, 1d0, rerho_pn)
            call triprod('n', vp, 't', imy, 't', un, 1d0, 1d0, imrho_pn)
            
            ! kappa(+)_pn = -up X(+)' un' + vp Y(+)' vn'
            call triprod('n', up, 't', rex, 't', un, -1d0, 0d0, rekp)
            call triprod('n', up, 't', imx, 't', un, -1d0, 0d0, imkp)
            call triprod('n', vp, 't', rey, 't', vn, 1d0, 1d0, rekp)
            call triprod('n', vp, 't', imy, 't', vn, 1d0, 1d0, imkp)
            
            ! kappa(-)_pn = vp X(+)*' vn' - up Y(+)'* un'
            call triprod('n', vp, 't', rex, 't', vn, 1d0, 0d0, rekm)
            call triprod('n', vp, 't', imx, 't', vn, -1d0, 0d0, imkm)
            call triprod('n', up, 't', rey, 't', un, -1d0, 1d0, rekm)
            call triprod('n', up, 't', imy, 't', un, 1d0, 1d0, imkm)
            
            ! Finite temperature
            if (ft_active) then
               ! rho_np += un P(+) up' - vn Q(+) vp'
               call triprod('n', un, 'n', rep, 't', up, 1d0, 1d0, rerho_np)
               call triprod('n', un, 'n', imp, 't', up, 1d0, 1d0, imrho_np)
               call triprod('n', vn, 'n', req, 't', vp,-1d0, 1d0, rerho_np)
               call triprod('n', vn, 'n', imq, 't', vp,-1d0, 1d0, imrho_np)
               
               ! rho_pn += Up Q(+)' Un' - Vp P(+)' Vn'
               call triprod('n', up, 't', req, 't', un, 1d0, 1d0, rerho_pn)
               call triprod('n', up, 't', imq, 't', un, 1d0, 1d0, imrho_pn)
               call triprod('n', vp, 't', rep, 't', vn,-1d0, 1d0, rerho_pn)
               call triprod('n', vp, 't', imp, 't', vn,-1d0, 1d0, imrho_pn)
               
               ! kappa(+)_pn += up Q(+)' vn' - vp P(+)' un'
               call triprod('n', up, 't', req, 't', vn, 1d0, 1d0, rekp)
               call triprod('n', up, 't', imq, 't', vn, 1d0, 1d0, imkp)
               call triprod('n', vp, 't', rep, 't', un,-1d0, 1d0, rekp)
               call triprod('n', vp, 't', imp, 't', un,-1d0, 1d0, imkp)
               
               ! kappa(-)_pn += -vp Q(+)'* un' + up P(+)'* vn'
               call triprod('n', vp, 't', req, 't', un,-1d0, 1d0, rekm)
               call triprod('n', vp, 't', imq, 't', un, 1d0, 1d0, imkm)
               call triprod('n', up, 't', rep, 't', vn, 1d0, 1d0, rekm)
               call triprod('n', up, 't', imp, 't', vn,-1d0, 1d0, imkm)
            end if            
         end if

         ! Compute the perturbed fields (functional derivatives of the EDF)
         ! in the coordinate space and then the mean fields h_pn and h_np
         call density(rerho_pn, imrho_pn, rekp, imkp, rhox)
         call meanfield(rhox, reh_pn, imh_pn)
         call pairingfield(rhox, redp, imdp)    ! Delta+ from rho-tilde+(r)
         
         call density(rerho_np, imrho_np, rekm, imkm, rhox)
         call meanfield(rhox, reh_np, imh_np)
         call pairingfield(rhox, redm, imdm)    ! Delta- from rho-tilde-(r)
         
         ! Transform the Hamiltonian into quasiparticle basis
         if (bminus) then
            ! H20 = up' h(pn) vn - vp' h(np)' un
            !     + up' Delta+(pn) un - vp' Delta-* vn
            call triprod('t', up, 'n', reh_pn, 'n', vn, 1d0, 0d0, h20re)
            call triprod('t', up, 'n', imh_pn, 'n', vn, 1d0, 0d0, h20im)
            call triprod('t', vp, 't', reh_np, 'n', un, -1d0, 1d0, h20re)
            call triprod('t', vp, 't', imh_np, 'n', un, -1d0, 1d0, h20im)
            call triprod('t', up, 'n', redp, 'n', un, 1d0, 1d0, h20re)
            call triprod('t', up, 'n', imdp, 'n', un, 1d0, 1d0, h20im)
            call triprod('t', vp, 'n', redm, 'n', vn, -1d0, 1d0, h20re)
            call triprod('t', vp, 'n', imdm, 'n', vn, 1d0, 1d0, h20im)
            
            ! H02 = up' h(np)' vn - vp' h(pn) un
            !     + up' Delta-* un - vp' Delta+ vn
            call triprod('t', up, 't', reh_np, 'n', vn, 1d0, 0d0, h02re)
            call triprod('t', up, 't', imh_np, 'n', vn, 1d0, 0d0, h02im)
            call triprod('t', vp, 'n', reh_pn, 'n', un, -1d0, 1d0, h02re)
            call triprod('t', vp, 'n', imh_pn, 'n', un, -1d0, 1d0, h02im)
            call triprod('t', up, 'n', redm, 'n', un, 1d0, 1d0, h02re)
            call triprod('t', up, 'n', imdm, 'n', un, -1d0, 1d0, h02im)
            call triprod('t', vp, 'n', redp, 'n', vn, -1d0, 1d0, h02re)
            call triprod('t', vp, 'n', imdp, 'n', vn, -1d0, 1d0, h02im)
            
            ! Finite temperature
            if (ft_active) then
               ! H11  =  up' h(pn) un - vp' h(np)' vn + up' d(+)  vn - vp' d(-)*  un
               call triprod('t', up, 'n', reh_pn, 'n', un, 1d0, 0d0, h11re)
               call triprod('t', up, 'n', imh_pn, 'n', un, 1d0, 0d0, h11im)
               call triprod('t', vp, 't', reh_np, 'n', vn,-1d0, 1d0, h11re)
               call triprod('t', vp, 't', imh_np, 'n', vn,-1d0, 1d0, h11im)
               call triprod('t', up, 'n', redp,   'n', vn, 1d0, 1d0, h11re)
               call triprod('t', up, 'n', imdp,   'n', vn, 1d0, 1d0, h11im)
               call triprod('t', vp, 'n', redm,   'n', un,-1d0, 1d0, h11re)
               call triprod('t', vp, 'n', imdm,   'n', un, 1d0, 1d0, h11im)
               
               ! H11~ = -vp' h(pn) vn + up' h(np)' un - vp' d(+)   un + up' d(-)* vn
               call triprod('t', vp, 'n', reh_pn, 'n', vn,-1d0, 0d0, h11tre)
               call triprod('t', vp, 'n', imh_pn, 'n', vn,-1d0, 0d0, h11tim)
               call triprod('t', up, 't', reh_np, 'n', un, 1d0, 1d0, h11tre)
               call triprod('t', up, 't', imh_np, 'n', un, 1d0, 1d0, h11tim)
               call triprod('t', vp, 'n', redp,   'n', un,-1d0, 1d0, h11tre)
               call triprod('t', vp, 'n', imdp,   'n', un,-1d0, 1d0, h11tim)
               call triprod('t', up, 'n', redm,   'n', vn, 1d0, 1d0, h11tre)
               call triprod('t', up, 'n', imdm,   'n', vn,-1d0, 1d0, h11tim)
            end if
         else
            ! H20 = un' h(np) vn - vn' h(pn)' up
            !     + un' Delta+(np) un - vn' Delta-* vp
            call triprod('t', un, 'n', reh_np, 'n', vp, 1d0, 0d0, h20re)
            call triprod('t', un, 'n', imh_np, 'n', vp, 1d0, 0d0, h20im)
            call triprod('t', vn, 't', reh_pn, 'n', up, -1d0, 1d0, h20re)
            call triprod('t', vn, 't', imh_pn, 'n', up, -1d0, 1d0, h20im)
            call triprod('t', un, 'n', redp, 'n', up, 1d0, 1d0, h20re)
            call triprod('t', un, 'n', imdp, 'n', up, 1d0, 1d0, h20im)
            call triprod('t', vn, 'n', redm, 'n', vp, -1d0, 1d0, h20re)
            call triprod('t', vn, 'n', imdm, 'n', vp, 1d0, 1d0, h20im)
            
            ! H02 = un' h(pn)' vp - vn' h(np) up
            !     + un' Delta-* up - vn' Delta+ vp
            call triprod('t', un, 't', reh_pn, 'n', vp, 1d0, 0d0, h02re)
            call triprod('t', un, 't', imh_pn, 'n', vp, 1d0, 0d0, h02im)
            call triprod('t', vn, 'n', reh_np, 'n', up, -1d0, 1d0, h02re)
            call triprod('t', vn, 'n', imh_np, 'n', up, -1d0, 1d0, h02im)
            call triprod('t', un, 'n', redm, 'n', up, 1d0, 1d0, h02re)
            call triprod('t', un, 'n', imdm, 'n', up, -1d0, 1d0, h02im)
            call triprod('t', vn, 'n', redp, 'n', vp, -1d0, 1d0, h02re)
            call triprod('t', vn, 'n', imdp, 'n', vp, -1d0, 1d0, h02im)
            
            ! Finite temperature
            if (ft_active) then
               ! H11  =  un' h(np) up - vn' h(pn)' vp + un' d(+) vp - vn' d(-)* up
               call triprod('t', un, 'n', reh_np, 'n', up, 1d0, 0d0, h11re)
               call triprod('t', un, 'n', imh_np, 'n', up, 1d0, 0d0, h11im)
               call triprod('t', vn, 't', reh_pn, 'n', vp,-1d0, 1d0, h11re)
               call triprod('t', vn, 't', imh_pn, 'n', vp,-1d0, 1d0, h11im)
               call triprod('t', un, 'n', redp,   'n', vp, 1d0, 1d0, h11re)
               call triprod('t', un, 'n', imdp,   'n', vp, 1d0, 1d0, h11im)
               call triprod('t', vn, 'n', redm,   'n', up,-1d0, 1d0, h11re)
               call triprod('t', vn, 'n', imdm,   'n', up, 1d0, 1d0, h11im)
               
               ! H11~ = -vn' h(np) vp + un' h(pn)' up - vn' d(+) up + un' d(-)* vp
               call triprod('t', vn, 'n', reh_np, 'n', vp,-1d0, 0d0, h11tre)
               call triprod('t', vn, 'n', imh_np, 'n', vp,-1d0, 0d0, h11tim)
               call triprod('t', un, 't', reh_pn, 'n', up, 1d0, 1d0, h11tre)
               call triprod('t', un, 't', imh_pn, 'n', up, 1d0, 1d0, h11tim)
               call triprod('t', vn, 'n', redp,   'n', up,-1d0, 1d0, h11tre)
               call triprod('t', vn, 'n', imdp,   'n', up,-1d0, 1d0, h11tim)
               call triprod('t', un, 'n', redm,   'n', vp, 1d0, 1d0, h11tre)
               call triprod('t', un, 'n', imdm,   'n', vp,-1d0, 1d0, h11tim)
            end if
         end if
         
         
         ! Add the external field (assuming real f20, f02 here)
         h20re%elem = f20%elem + h20re%elem
         h02re%elem = f02%elem + h02re%elem
         
         ! Finite temperature
         if (ft_active) then
            h11re%elem  = f11%elem  + h11re%elem
            h11tre%elem = f11t%elem + h11tre%elem
         end if
         end if  ! End computation of the residual interaction
         
         
         ! Solve the new X and Y
         call complex_multiply(greenXre, greenXim, h20re, h20im, rex, imx)
         call complex_multiply(greenYre, greenYim, h02re, h02im, rey, imy)
         
         ! Finite temperature
         if (ft_active) then
            ! Solve the new P and Q
            call complex_multiply(greenPre, greenPim, h11re,  h11im,  rep, imp)
            call complex_multiply(greenQre, greenQim, h11tre, h11tim, req, imq)
            
            ! There are additional real, T-dependent constants which get multiplied
            rex%elem(:) = ft_1fafb_x%elem(:)*rex%elem(:)
            imx%elem(:) = ft_1fafb_x%elem(:)*imx%elem(:)
            rey%elem(:) = ft_1fafb_y%elem(:)*rey%elem(:)
            imy%elem(:) = ft_1fafb_y%elem(:)*imy%elem(:)
            
            rep%elem(:) = ft_fbfa_p%elem(:)*rep%elem(:)
            imp%elem(:) = ft_fbfa_p%elem(:)*imp%elem(:)
            req%elem(:) = ft_fbfa_q%elem(:)*req%elem(:)
            imq%elem(:) = ft_fbfa_q%elem(:)*imq%elem(:)
         end if
         
         ! Use Broyden mixing to stabilize convergence
         ! Finite temperature has twice the number of vectors to handle P and Q
         if (ft_active) then
            call qrpa_broyden(iter, nxy, rex%elem, rey%elem, imx%elem, imy%elem, &
               rep%elem, imp%elem, req%elem, imq%elem)
         else
            call qrpa_broyden(iter, nxy, rex%elem, rey%elem, imx%elem, imy%elem)            
         end if
         
         
         ! Strength function and cross-terms
         if (.not.ft_active) then
            ! These are the traditional strength functions
            str(1) = -1/pi*(dot_product(f20%elem,imx%elem) + dot_product(f02%elem,imy%elem))
         
            if (present(cstr)) then
               cstr(1) = iu*str(1) - &
                  1/pi*(dot_product(f20%elem,rex%elem) + dot_product(f02%elem,rey%elem))
            end if
            
            ! Cross-terms
            if (allocated(g20)) then
               do ixterm=1, size(g20)
                  str(1+ixterm) = -1/pi * &
                     (dot_product(g20(ixterm)%elem, imx%elem) + &
                      dot_product(g02(ixterm)%elem, imy%elem))

                  if (present(cstr)) then
                     cstr(1+ixterm) = iu*str(1+ixterm) - 1/pi * &
                        (dot_product(g20(ixterm)%elem, rex%elem) + &
                         dot_product(g02(ixterm)%elem, rey%elem))
                  end if
               end do
            end if
         ! Finite temperature
         else if (ft_active) then
            str(1) = -1/pi * &
               (dot_product(f20%elem,imx%elem) + dot_product(f02%elem,imy%elem) + &
                dot_product(f11%elem,imp%elem) + dot_product(f11t%elem,imq%elem))
            
            ! The full complex strength function uses the full complex energy
            ! and the full complex X, Y, P, and Q
            if (present(cstr)) then
               cstr(1) = iu*str(1) - 1/pi * &
                  (dot_product(f20%elem, rex%elem) + &
                   dot_product(f02%elem, rey%elem) + &
                   dot_product(f11%elem, rep%elem) + &
                   dot_product(f11t%elem,req%elem))
            end if
            
            ! Cross-terms
            if (allocated(g20)) then
               do ixterm=1, size(g20)
                  str(1+ixterm) = -1/pi * &
                     (dot_product(g20(ixterm)%elem, imx%elem) + &
                      dot_product(g02(ixterm)%elem, imy%elem) + &
                      dot_product(g11(ixterm)%elem, imp%elem) + &
                      dot_product(g11t(ixterm)%elem,imq%elem))

                  if (present(cstr)) then
                     cstr(1+ixterm) = iu*str(1+ixterm) - 1/pi * &
                        (dot_product(g20(ixterm)%elem, rex%elem) + &
                         dot_product(g02(ixterm)%elem, rey%elem) + &
                         dot_product(g11(ixterm)%elem, rep%elem) + &
                         dot_product(g11t(ixterm)%elem,req%elem))
                  end if
               end do
            end if
            
            ! T-dependent prefactor
            ! Note that this has a pole at E=0, and as a result we need to
            ! abandon the computation of the point E = 0 + i\gamma when the
            ! finite-temperature is active. If we are absolutely at (0,0),
            ! the computation should be dropped due to the pole.
            if (abs(cmplx(omega,width,dp)) < 1.0d-6) then
               write(*,'(a)') 'ERROR: cannot compute the strength function at&
                  & finite temperature for E = 0 + i0 MeV.'
               iter_used = -2
               str(:) = 0
               if (present(cstr)) cstr(:) = 0
               stop
            else if (abs(omega) < 1.0d-6) then
               if (iter == 0) write(*,'(a)') 'Warning: cannot compute the real-valued&
                   & strength function at finite temperature for Re(E) = 0 MeV.'
               str(:) = 0
               if (.not.present(cstr)) then
                  iter_used = -2
                  write(*,'(a)') 'Abandoning this point since a complex strength&
                      & function was not requested.'
                  exit
               end if
            end if
            
            ! Actually compute the prefactors
            ! 1) Compute the prefactor which needs Re(E) only
            ! 2) Compute the full complex prefactor
            ! If Re(E)/T is too large, we have problems. In such a case, fall
            ! back to factor = 1/(1-e(-E/T)) ~ 1/(1-1/e(<huge>)) ~ 1
            if (omega/ft_temp > log(huge(1.0_dp))) then
               ft_real_factor  = 1.0_dp
               ft_cmplx_factor = 1.0_dp
               if (iter == 0) then
                  write(*,'(a)') 'Warning: setting ft_factor = 1 to avoid overflows.'
               end if
            else
               ft_real_factor  = (1-exp(-omega/ft_temp))**(-1)
               ft_cmplx_factor = (1-exp(-cmplx(omega,width,dp)/ft_temp))**(-1)
            end if
            
            ! Multiply in the prefactors
            str(:) = str(:)*ft_real_factor
            if (present(cstr)) cstr(:) = cstr(:)*ft_cmplx_factor
         end if
         
         ! Check for convergence
         if (si < convergence_epsilon) then
            iter_used = iter
            exit
         end if
         
      end do
      
      ! Disabled unless requested --- this routine prints out the important contributors
      ! to each FAM computation (e.g. near an eigenvalue it shows the main contributions
      ! to that transition).
      ! call transition_decomposition(cmplx(omega,width,dp), iter, si, vp, vn)
      
   end subroutine pnfam_solve
   
   
   !---------------------------------------------------------------------------
   ! Allocates and initializes the data structures used by pnfam_solve.
   !---------------------------------------------------------------------------
   subroutine init_pnfam_solver(f, vp, up, vn, un, bminus, g)
      use constants,    only : IT_NEUTRON, IT_PROTON
      use hfbtho_basis, only : ft_active
      
      implicit none
      
      type(blockmatrix), intent(in) :: f, vp, up, vn, un
      type(external_field), dimension(:), allocatable, intent(in) :: g
      logical, intent(in) :: bminus
      
      integer :: ixterm
      
      
      ! While the different FAM matrices may have a different non-trivial block
      ! structure, they still share the same number of elements in the non-trivial
      ! blocks:
      nxy = size(f%elem)
      
      ! Allocate needed arrays
      call allocate_blockmatrix(f20, nxy)
      call allocate_blockmatrix(f02, nxy)
      
      ! Special allocation for the cross-terms
      if (allocated(g)) then
         allocate(g20(1:size(g)), g02(1:size(g)))
         
         do ixterm=1, size(g)       
            call allocate_blockmatrix(g20(ixterm), nxy)
            call allocate_blockmatrix(g02(ixterm), nxy)
         end do
         
         ! Finite temperature
         if (ft_active) then
            allocate(g11(1:size(g)), g11t(1:size(g)))
            
            do ixterm=1, size(g)       
               call allocate_blockmatrix(g11(ixterm),  nxy)
               call allocate_blockmatrix(g11t(ixterm), nxy)
            end do
         end if
         
      end if
   
      call allocate_blockmatrix(h20re, nxy)
      call allocate_blockmatrix(h20im, nxy)
      call allocate_blockmatrix(h02re, nxy)
      call allocate_blockmatrix(h02im, nxy)
   
      call allocate_blockmatrix(rex, nxy)
      call allocate_blockmatrix(imx, nxy)
      call allocate_blockmatrix(rey, nxy)
      call allocate_blockmatrix(imy, nxy)
   
      call allocate_blockmatrix(greenXre, nxy)
      call allocate_blockmatrix(greenXim, nxy)
      call allocate_blockmatrix(greenYre, nxy)
      call allocate_blockmatrix(greenYim, nxy)
   
      call allocate_blockmatrix(reh_np, nxy)
      call allocate_blockmatrix(imh_np, nxy)
      call allocate_blockmatrix(reh_pn, nxy)
      call allocate_blockmatrix(imh_pn, nxy)
      call allocate_blockmatrix(rerho_pn, nxy)
      call allocate_blockmatrix(imrho_pn, nxy)
      call allocate_blockmatrix(rerho_np, nxy)
      call allocate_blockmatrix(imrho_np, nxy)
   
      call allocate_blockmatrix(redp, nxy)
      call allocate_blockmatrix(imdp, nxy)
      call allocate_blockmatrix(redm, nxy)
      call allocate_blockmatrix(imdm, nxy)
      call allocate_blockmatrix(rekp, nxy)
      call allocate_blockmatrix(imkp, nxy)
      call allocate_blockmatrix(rekm, nxy)
      call allocate_blockmatrix(imkm, nxy)
      
      ! Finite temperature
      if (ft_active) then
         call allocate_blockmatrix(f11,        nxy)    ! Diagonal blocks of F
         call allocate_blockmatrix(f11t,       nxy)
         call allocate_blockmatrix(rep,        nxy)    ! Diagonal blocks of \delta R
         call allocate_blockmatrix(imp,        nxy)
         call allocate_blockmatrix(req,        nxy)
         call allocate_blockmatrix(imq,        nxy)
         call allocate_blockmatrix(h11re,      nxy)    ! Diagonal blocks of \delta H
         call allocate_blockmatrix(h11im,      nxy)
         call allocate_blockmatrix(h11tre,     nxy)
         call allocate_blockmatrix(h11tim,     nxy)
         call allocate_blockmatrix(greenPre,   nxy)    ! Diagonal blocks of Green's function
         call allocate_blockmatrix(greenPim,   nxy)
         call allocate_blockmatrix(greenQre,   nxy)
         call allocate_blockmatrix(greenQim,   nxy)
         call allocate_blockmatrix(ft_fbfa_p,  nxy)   ! Statistical factors
         call allocate_blockmatrix(ft_fbfa_q,  nxy)
         call allocate_blockmatrix(ft_1fafb_x, nxy)
         call allocate_blockmatrix(ft_1fafb_y, nxy)
      end if
   
      ! Allocate the coordinate-space densities
      call allocate_density(rhox)
   
      ! Transform the external operator to the quasiparticle basis
      if (bminus) then
         call triprod('t', up, 'n', f, 'n', vn, 1d0, 0d0, f20)       !  up' f- vn
         call triprod('t', vp, 'n', f, 'n', un, -1d0, 0d0, f02)      ! -vp' f- un
         
         ! Finite temperature
         if (ft_active) then
            call triprod('t', up, 'n', f, 'n', un, 1d0, 0d0, f11)    !  up' f- un
            call triprod('t', vp, 'n', f, 'n', vn,-1d0, 0d0, f11t)   ! -vp' f- vn
         end if
      else
         call triprod('t', un, 'n', f, 'n', vp, 1d0, 0d0, f20)       !  un' f+ vp
         call triprod('t', vn, 'n', f, 'n', up, -1d0, 0d0, f02)      ! -vn' f+ un
         
         ! Finite temperature
         if (ft_active) then
            call triprod('t', un, 'n', f, 'n', up, 1d0, 0d0, f11)    !  un' f+ up
            call triprod('t', vn, 'n', f, 'n', vp,-1d0, 0d0, f11t)   ! -vn' f+ vp
         end if
      end if
      
      ! Special transformation for the cross-terms
      if (allocated(g)) then
         do ixterm=1, size(g)
            if (bminus) then
               call triprod('t',up,'n',g(ixterm)%mat,'n',vn, 1d0,0d0,g20(ixterm))
               call triprod('t',vp,'n',g(ixterm)%mat,'n',un,-1d0,0d0,g02(ixterm))
               
               ! Finite temperature
               if (ft_active) then
                  call triprod('t',up,'n',g(ixterm)%mat,'n',un, 1d0,0d0,g11(ixterm))
                  call triprod('t',vp,'n',g(ixterm)%mat,'n',vn,-1d0,0d0,g11t(ixterm))
               end if
            else
               call triprod('t',un,'n',g(ixterm)%mat,'n',vp, 1d0,0d0,g20(ixterm))
               call triprod('t',vn,'n',g(ixterm)%mat,'n',up,-1d0,0d0,g02(ixterm))
               
               ! Finite temperature
               if (ft_active) then
                  call triprod('t',un,'n',g(ixterm)%mat,'n',up, 1d0,0d0,g11(ixterm))
                  call triprod('t',vn,'n',g(ixterm)%mat,'n',vp,-1d0,0d0,g11t(ixterm))
               end if
            end if
         end do
      end if
      
      ! Initialize the block structure of various matrices
      ! This is needed for the block matrices that are not constructed with
      ! the routine triprod
      call copy_block_structure(f20, rex)  ! X has the same structure as F20
      call copy_block_structure(f20, imx)
      call copy_block_structure(f02, rey)  ! Y has the same structure as F02
      call copy_block_structure(f02, imy)
   
      call copy_block_structure(f, reh_pn)  ! h_pn has the same structure as f
      call copy_block_structure(f, imh_pn)
      call copy_block_structure(rex, redp)  ! Delta(+) has the same structure as X
      call copy_block_structure(imx, imdp)
      
      call copy_block_structure(f, reh_np)  ! h_np has the same structure as f
      call copy_block_structure(f, imh_np)
      call copy_block_structure(rey, redm)  ! Delta(-) has the same structure as Y
      call copy_block_structure(imy, imdm)
      
      call copy_block_structure(rex, greenXre)  ! -1/(Ep+En-omega) has the same structure as X
      call copy_block_structure(imx, greenXim)
      call copy_block_structure(rey, greenYre)  ! -1/(Ep+En+omega*) has the same structure as Y
      call copy_block_structure(imy, greenYim)
      
      call copy_block_structure(f20, h20re)  ! delta-H20 has the same structure as F20
      call copy_block_structure(f20, h20im)
      call copy_block_structure(f02, h02re)  ! delta-H02 has the same structure as F02
      call copy_block_structure(f02, h02im)
      
      ! Finite temperature
      if (ft_active) then
         call copy_block_structure(f11,  ReP)
         call copy_block_structure(f11,  ImP)
         call copy_block_structure(f11,  greenPre)
         call copy_block_structure(f11,  greenPim)
         call copy_block_structure(f11,  h11re)
         call copy_block_structure(f11,  h11im)
         call copy_block_structure(f11t, ReQ)
         call copy_block_structure(f11t, ImQ)
         call copy_block_structure(f11t, greenQre)
         call copy_block_structure(f11t, greenQim)
         call copy_block_structure(f11t, h11tre)
         call copy_block_structure(f11t, h11tim)
         
         ! Construct the matrices (f_b-f_a) and (1-f_a-f_b) which (like the
         ! Green's functions) match the P, Q, X, Y matrices in structure
         call copy_block_structure(f11,  ft_fbfa_p)
         call copy_block_structure(f11t, ft_fbfa_q)
         call copy_block_structure(f20,  ft_1fafb_x)
         call copy_block_structure(f02,  ft_1fafb_y)
         
         if (bminus .eqv. .true.) then
            call ft_1fafb_matrix(ta=IT_PROTON, tb=IT_NEUTRON, mat=ft_1fafb_x)
            call ft_1fafb_matrix(ta=IT_PROTON, tb=IT_NEUTRON, mat=ft_1fafb_y)
            call ft_fbfa_matrix (ta=IT_PROTON, tb=IT_NEUTRON, mat=ft_fbfa_p)
            call ft_fbfa_matrix (ta=IT_PROTON, tb=IT_NEUTRON, mat=ft_fbfa_q)
         else
            call ft_1fafb_matrix(ta=IT_NEUTRON, tb=IT_PROTON, mat=ft_1fafb_x)
            call ft_1fafb_matrix(ta=IT_NEUTRON, tb=IT_PROTON, mat=ft_1fafb_y)
            call ft_fbfa_matrix (ta=IT_NEUTRON, tb=IT_PROTON, mat=ft_fbfa_p)
            call ft_fbfa_matrix (ta=IT_NEUTRON, tb=IT_PROTON, mat=ft_fbfa_q)
         end if
      end if
      
   end subroutine init_pnfam_solver
   
   
   subroutine deallocate_pnfam_solver(g)
      use hfbtho_basis, only : ft_active
      implicit none
      type(external_field), dimension(:), allocatable :: g
      
      integer :: ixterm
      
      call deallocate_blockmatrix(f20)
      call deallocate_blockmatrix(f02)
      
      ! Special allocation for the cross-terms
      if (allocated(g)) then
         
         do ixterm=1, size(g)       
            call deallocate_blockmatrix(g20(ixterm))
            call deallocate_blockmatrix(g02(ixterm))
         end do
         
         ! Finite temperature
         if (ft_active) then
            
            do ixterm=1, size(g)       
               call deallocate_blockmatrix(g11(ixterm))
               call deallocate_blockmatrix(g11t(ixterm))
            end do
            deallocate(g11, g11t)
            
         end if
         
         deallocate(g20, g02)
         
      end if
   
      call deallocate_blockmatrix(h20re)
      call deallocate_blockmatrix(h20im)
      call deallocate_blockmatrix(h02re)
      call deallocate_blockmatrix(h02im)
   
      call deallocate_blockmatrix(rex)
      call deallocate_blockmatrix(imx)
      call deallocate_blockmatrix(rey)
      call deallocate_blockmatrix(imy)
   
      call deallocate_blockmatrix(greenXre)
      call deallocate_blockmatrix(greenXim)
      call deallocate_blockmatrix(greenYre)
      call deallocate_blockmatrix(greenYim)
   
      call deallocate_blockmatrix(reh_np)
      call deallocate_blockmatrix(imh_np)
      call deallocate_blockmatrix(reh_pn)
      call deallocate_blockmatrix(imh_pn)
      call deallocate_blockmatrix(rerho_pn)
      call deallocate_blockmatrix(imrho_pn)
      call deallocate_blockmatrix(rerho_np)
      call deallocate_blockmatrix(imrho_np)
   
      call deallocate_blockmatrix(redp)
      call deallocate_blockmatrix(imdp)
      call deallocate_blockmatrix(redm)
      call deallocate_blockmatrix(imdm)
      call deallocate_blockmatrix(rekp)
      call deallocate_blockmatrix(imkp)
      call deallocate_blockmatrix(rekm)
      call deallocate_blockmatrix(imkm)
      
      ! Finite temperature
      if (ft_active) then
         call deallocate_blockmatrix(f11)    ! Diagonal blocks of F
         call deallocate_blockmatrix(f11t)
         call deallocate_blockmatrix(rep)    ! Diagonal blocks of \delta R
         call deallocate_blockmatrix(imp)
         call deallocate_blockmatrix(req)
         call deallocate_blockmatrix(imq)
         call deallocate_blockmatrix(h11re)    ! Diagonal blocks of \delta H
         call deallocate_blockmatrix(h11im)
         call deallocate_blockmatrix(h11tre)
         call deallocate_blockmatrix(h11tim)
         call deallocate_blockmatrix(greenPre)    ! Diagonal blocks of Green's function
         call deallocate_blockmatrix(greenPim)
         call deallocate_blockmatrix(greenQre)
         call deallocate_blockmatrix(greenQim)
         call deallocate_blockmatrix(ft_fbfa_p)   ! Statistical factors
         call deallocate_blockmatrix(ft_fbfa_q)
         call deallocate_blockmatrix(ft_1fafb_x)
         call deallocate_blockmatrix(ft_1fafb_y)
      end if
      
   end subroutine deallocate_pnfam_solver
   
   
   !---------------------------------------------------------------------------
   ! A wrapper for the more general modified Broyden mixing routine
   !---------------------------------------------------------------------------
   subroutine qrpa_broyden(niter, nuv, x1, x2, x3, x4, x5, x6, x7, x8)
      use broyden_mixer, only : qrpa_alphamix, qrpa_broyden_method, si, &
         broyden_history_size, qrpa_bbroyden
      implicit none
      integer, intent(in) :: nuv
      real(dp), intent(inout)           :: x1(nuv), x2(nuv), x3(nuv), x4(nuv)
      real(dp), intent(inout), optional :: x5(nuv), x6(nuv), x7(nuv), x8(nuv)
      integer, intent(in)  :: niter
      integer :: i, nvec
      real(dp), allocatable, save :: qrpa_broin(:), qrpa_broout(:)
      
      ! Calculate the size of the broyden vectors.
      ! 5, 6, 7, 8 are present at finite temperature
      i = 4*nuv
      if (present(x5)) i = i+nuv
      if (present(x6)) i = i+nuv
      if (present(x7)) i = i+nuv
      if (present(x8)) i = i+nuv
      nvec = i/nuv

      ! On the first iteration
      if (niter==0) then
         if (allocated(qrpa_broin)) deallocate(qrpa_broin,qrpa_broout)
         allocate(qrpa_broin(i),qrpa_broout(i))
         qrpa_broin = 0
      end if

      ! The Broyden mixing routine operates on a real array; Therefore we need
      ! to unpack the complex values to a real array, call the Broyden routine,
      ! and then repack the real and imaginary parts to the complex arrays
      qrpa_broout(1:nuv) = x1
      qrpa_broout(nuv+1:2*nuv) = x2
      qrpa_broout(2*nuv+1:3*nuv) = x3
      qrpa_broout(3*nuv+1:4*nuv) = x4
      
      ! nvec should be either 4 (T=0) or 8 (T>0)
      if (nvec /= 4 .and. nvec /= 8) then
         write(*,'(a)') 'Error: nvec /= 4 or 8 in qrpa_broyden'
         stop
      end if
      
      ! Set up the Broyden vectors
      if (nvec == 8) then
         qrpa_broout(4*nuv+1:5*nuv) = x5
         qrpa_broout(5*nuv+1:6*nuv) = x6
         qrpa_broout(6*nuv+1:7*nuv) = x7
         qrpa_broout(7*nuv+1:8*nuv) = x8
      end if
      
      ! Broyden
      call qrpa_broyden_method(i,qrpa_broout,qrpa_broin,qrpa_alphamix,si,niter,broyden_history_size,qrpa_bbroyden)
      if (niter==0) si = 1

      x1 = qrpa_broin(1:nuv)
      x2 = qrpa_broin(nuv+1:2*nuv)
      x3 = qrpa_broin(2*nuv+1:3*nuv)
      x4 = qrpa_broin(3*nuv+1:4*nuv)

      ! Set up the Broyden vectors
      if (nvec == 8) then
         qrpa_broin(4*nuv+1:5*nuv) = x5
         qrpa_broin(5*nuv+1:6*nuv) = x6
         qrpa_broin(6*nuv+1:7*nuv) = x7
         qrpa_broin(7*nuv+1:8*nuv) = x8
      end if

   end subroutine qrpa_broyden


   !---------------------------------------------------------------------------
   ! The "Green functions" in the FAM equations.  The minus sign has been
   ! absorbed here.  This version has been extended for general K and parity.
   !---------------------------------------------------------------------------
   subroutine greenx(E1, E2, omega, greenRe, greenIm)
      use hfbtho_basis
      implicit none
      real(dp), intent(in) :: E1(:), E2(:)
      complex(dp), intent(in) :: omega
      type(blockmatrix), intent(inout) :: greenRe, greenIm

      integer :: ir, ic, ipt, i1, i2, i1a, i1b, i2a, i2b
      
      ipt = 1
      do ir=1,nb   ! block row index
         ic = greenRe%ir2c(ir)   ! block column index
         if (ic == 0) cycle
         i1a = isstart(ir) ; i1b = i1a + db(ir) - 1
         i2a = isstart(ic) ; i2b = i2a + db(ic) - 1
         do i2 = i2a,i2b   ! true column index
         do i1 = i1a,i1b   ! true row index
            greenRe%elem(ipt) = real(-1/(E1(i1) + E2(i2) + omega), dp)
            greenIm%elem(ipt) = aimag(-1/(E1(i1) + E2(i2) + omega))
            ipt = ipt + 1
         end do
         end do
      end do

   end subroutine greenx
   
      
   ! The moment of inertia of the HFB solution in the cranking approximation
   ! Call this separately for protons and neutrons, add the results up
   function moment_of_inertia(E, U, V) result(moi)
      use external_field_type
      use extfield
      use hfbtho_basis
      implicit none
      real(dp), intent(in) :: E(:)
      type(blockmatrix), intent(in) :: U, V  ! qp energies and amplitudes
      real(dp) :: moi
      
      type(blockmatrix) :: M1, M2
      type(external_field) :: Jplus, Jminus
      integer :: ipt, ibx1, ibx2, i1, i2, nd1, nd2, ix1, ix2
      
      ! Finite temperature
      if (ft_active) then
         call writelog('Warning: moment_of_inertia has not been updated for finite-temperature.')
      end if
      
      call writelog('Warning: moment_of_inertia has not been checked for using the PWI.')
      
      ! Because the blockmatrix type requires there to be only one nonzero block
      ! in each row and column, we need to separate the matrix into two:
      ! M1 = U^dagger J_+ V^* - V^dagger J_- U^*
      ! M2 = U^dagger J_- V^* - V^dagger J_+ U^*
      
      ! Construct raising and lowering operators of the total angular momentum
      Jplus%k = 1 ; Jplus%parity_even = .true. ; Jplus%rank = 1
      call init_fam_mapping(Jplus)
      ipt = 1
      do ibx1 = 1,nb
         ibx2 = Jplus%mat%ir2c(ibx1)
         if (ibx2 == 0) cycle
         nd1 = db(ibx1) ; nd2 = db(ibx2)
         if (nd1 == 0 .or. nd2 == 0) cycle
         
         do i2=1, nd2
            ! Calculate the true index of the state |2>
            ix2 = i2 + isstart(ibx2) - 1

            do i1=1, nd1
               ! Calculate the true index of the state <1|
               ix1 = i1 + isstart(ibx1) - 1
               
               !if (nr(ix1)==nr(ix2).and.nz(ix1)==nz(ix2)) then
                  if (nl(ix1)==nl(ix2) .and. ns(ix1)==ns(ix2)+2) then
                     Jplus%mat%elem(ipt) = dot_product(wf(:,ix1), wf(:,ix2))  ! ok
                  end if
                  if (nl(ix1)==nl(ix2)+1 .and. ns(ix1)==ns(ix2)) then
                     Jplus%mat%elem(ipt) = dot_product(-wf(:,ix1)/y(:), wfdz(:,ix2)) &
                        + dot_product(wf(:,ix1)*z(:), wfdr(:,ix2)) &
                        - nl(ix2)*dot_product(wf(:,ix1)*z(:)*y(:), wf(:,ix2))
                  end if
               !end if
               
               ipt = ipt + 1
            end do
         end do
      end do
      
      Jminus%k = -1 ; Jminus%parity_even = .true. ; Jminus%rank = 1
      call init_fam_mapping(Jminus)
      ipt = 1
      do ibx1 = 1,nb
         ibx2 = Jminus%mat%ir2c(ibx1)
         if (ibx2 == 0) cycle
         nd1 = db(ibx1) ; nd2 = db(ibx2)
         if (nd1 == 0 .or. nd2 == 0) cycle
         
         do i2=1, nd2
            ! Calculate the true index of the state |2>
            ix2 = i2 + isstart(ibx2) - 1

            do i1=1, nd1
               ! Calculate the true index of the state <1|
               ix1 = i1 + isstart(ibx1) - 1
               
               !if (nr(ix1)==nr(ix2).and.nz(ix1)==nz(ix2)) then
                  if (nl(ix1)==nl(ix2) .and. ns(ix1)==ns(ix2)-2) then
                     Jminus%mat%elem(ipt) = dot_product(wf(:,ix1)/y(:), wfdz(:,ix2)) &
                        - dot_product(z(:)*wf(:,ix1), wfdr(:,ix2)) &
                        - nl(ix2)*dot_product(wf(:,ix1)*z(:)*y(:), wf(:,ix2))
                  end if
                  if (nl(ix1)==nl(ix2)-1 .and. ns(ix1)==ns(ix2)) then
                     Jminus%mat%elem(ipt) = dot_product(wf(:,ix1), wf(:,ix2))
                  end if
               !end if

               ipt = ipt + 1
            end do
         end do
      end do
      
      ! Prep the matrices M1 and M2
      call allocate_blockmatrix(M1,size(Jplus%mat%elem))
      call allocate_blockmatrix(M2,size(Jplus%mat%elem))
      
      ! Transform them into the qp basis
      call triprod('T', U, 'N', Jplus%mat, 'N', V, 1d0, 0d0, M1)
      call triprod('T', V, 'N', Jminus%mat, 'N', U, -1d0, 1d0, M1)
      call triprod('T', U, 'N', Jminus%mat, 'N', V, 1d0, 0d0, M2)
      call triprod('T', V, 'N', Jplus%mat, 'N', U, -1d0, 1d0, M2)
      
      ! Sum over all qp pairs, block by block
      moi = 0
      ipt = 1
      do ibx1 = 1,nb
         ibx2 = M1%ir2c(ibx1)
         if (ibx2 == 0) cycle
         nd1 = db(ibx1) ; nd2 = db(ibx2)
         if (nd1 == 0 .or. nd2 == 0) cycle
         
         do i2=1, nd2
            ix2 = i2 + isstart(ibx2) - 1
            do i1=1, nd1
               ix1 = i1 + isstart(ibx1) - 1
               moi = moi + M1%elem(ipt)**2/(E(ix1) + E(ix2))
               ipt = ipt + 1
            end do
         end do
      end do
      ipt = 1
      do ibx1 = 1,nb
         ibx2 = M2%ir2c(ibx1)
         if (ibx2 == 0) cycle
         nd1 = db(ibx1) ; nd2 = db(ibx2)
         if (nd1 == 0 .or. nd2 == 0) cycle
         
         do i2=1, nd2
            ix2 = i2 + isstart(ibx2) - 1
            do i1=1, nd1
               ix1 = i1 + isstart(ibx1) - 1
               moi = moi + M2%elem(ipt)**2/(E(ix1) + E(ix2))
               ipt = ipt + 1
            end do
         end do
      end do
      
      moi = moi/4
      
   end function moment_of_inertia
   
   !----------------------------------------------------------------------------
   ! Compute and print out (to '2qp.dat') the most important quasiparticle
   ! levels in any FAM computation. 
   !----------------------------------------------------------------------------
   subroutine transition_decomposition(eqrpa, iters, si, vp, vn)
      use constants,    only : strpar
      use hfbtho_basis, only : ns, nl, npar, isstart, get_hfbtho_occupations
      implicit none
      
      ! Cutoff on when to print this output (B >= strcut) and which levels
      ! to print (relative contribution >= relcut)
      real(dp), parameter :: strcut = -huge(1d0), relcut = 1d-2

      complex(dp), intent(in) :: eqrpa
      real(dp),    intent(in) :: si
      integer,     intent(in) :: iters
      type(blockmatrix), intent(in) :: vp, vn
      
      integer  :: kp, kn, i1, i2, ixy, ifam, ipblock, inblock, ineut, iprot
      real(dp) :: occp(size(vn%elem)), occn(size(vn%elem))
      real(dp) :: total_strength, rel_strength
      
      
      total_strength = (dot_product(imx%elem,f20%elem) + dot_product(imy%elem,f02%elem))/pi
      if (total_strength < strcut) return
      
      
      open(24, file='2qp.dat', status='unknown', position='append')
      write(24,'("EQRPA = (",1f10.6,",",1f10.6,"),   B = ",1f10.6,",   ITER = ",1i0,",   SI = ", 1es12.4)') &
         real(eqrpa,kind=dp), aimag(eqrpa), total_strength, iters, si
      write(24,'(a)') repeat('-',100)
      write(24,'(3a8,2a11,2a14,a)') 'IXY', 'IPROT', 'INEUT', 'KPI PROT', 'KPI NEUT', &
         'OCC PROT', 'OCC NEUT', '   CONTRIB'
      
      
      call get_hfbtho_occupations(vp,occp)
      call get_hfbtho_occupations(vn,occn)
      
      
      ixy = 0
      do ifam=1, nb
         
         ipblock = ifam
         inblock = f20%ir2c(ifam)
         
         if (inblock == 0) cycle
         if (db(ipblock) == 0 .or. db(inblock) == 0) cycle
         
         ! Column-major, <p|O|n>
         do i2=1,db(inblock)
            ineut = isstart(inblock) - 1 + i2
            kn = 2*nl(ineut)+ns(ineut)
            
            do i1=1,db(ipblock)
               iprot = isstart(ipblock) - 1 + i1
               kp = 2*nl(iprot)+ns(iprot)
               
               ixy = ixy + 1
               rel_strength = (imx%elem(ixy)*f20%elem(ixy) + &
                  imy%elem(ixy)*f02%elem(ixy))/pi/total_strength
            
               if (rel_strength >= relcut) then
                  write(24,'(3i8,2(1i8,"/2",1a1),2f14.5,1f10.1,"%")') &
                     ixy, iprot, ineut, abs(kp), strpar(npar(iprot)), abs(kn), &
                     strpar(npar(ineut)), occp(iprot), occn(ineut), 100*rel_strength
               end if
            end do
         end do
      end do
      
      write(24,'(a)') repeat('-',100)
      close(24)
   end subroutine transition_decomposition
   
   
   !----------------------------------------------------------------------------
   ! Compute the finite-temperature statistical factor matrices (1 - f_a - f_b)
   !----------------------------------------------------------------------------
   subroutine ft_1fafb_matrix(ta, tb, mat)
      use constants, only : IT_NEUTRON, IT_PROTON
      use hfbtho_basis, only : pwi_start_n, pwi_start_p, pwi_dim_n, pwi_dim_p, &
         pwi_qp_n, pwi_qp_p, ft_active, ft_fp, ft_fn, isstart
      
      implicit none
      
      integer, intent(in) :: ta, tb
      type(blockmatrix), intent(inout) :: mat
      
      integer :: im, ibr, ibc, ia, ib, ik, ka, kb
      logical :: active_a, active_b
      
      integer,  pointer :: kstarta(:), kstartb(:), kdima(:), kdimb(:), kqpa(:), kqpb(:)
      real(dp), pointer :: fa(:), fb(:)
      
      
      mat%elem(:) = 0
      
      ! Check that finite temperature should be on; if not, exit with the defaults
      if (.not.ft_active) then
         write(*,'(a)') 'WARNING: finite temperature is not active, but ft_1fafb_matrix was called.'
         return
      end if
      
      ! Check isospin integers
      if  ((ta /= IT_NEUTRON .and. ta /= IT_PROTON)   &
      .or. (tb /= IT_NEUTRON .and. tb /= IT_PROTON)) then
         write(*,'(A)') 'ERROR: bad isospin integers in ft_1fafb_matrix.'
         return
      end if
      
      ! Set up pointers
      nullify(kstarta, kdima, kqpa, fa)
      if (ta == IT_NEUTRON) then
         kstarta => pwi_start_n;  kdima => pwi_dim_n;  kqpa => pwi_qp_n;  fa => ft_fn
      else if (ta == IT_PROTON) then
         kstarta => pwi_start_p;  kdima => pwi_dim_p;  kqpa => pwi_qp_p;  fa => ft_fp
      end if
      
      nullify(kstartb, kdimb, kqpb, fb)
      if (tb == IT_NEUTRON) then
         kstartb => pwi_start_n;  kdimb => pwi_dim_n;  kqpb => pwi_qp_n;  fb => ft_fn
      else if (tb == IT_PROTON) then
         kstartb => pwi_start_p;  kdimb => pwi_dim_p;  kqpb => pwi_qp_p;  fb => ft_fp
      end if

      im = 0
      do ibr=1, nb
         
         ibc = mat%ir2c(ibr)
         if (ibc == 0) cycle
         
         do ib=isstart(ibc),isstart(ibc)+db(ibc)-1
            ! Is the q.p. with index ib active?
            active_b = .false.
            do ik=kstartb(ibc),kstartb(ibc)+kdimb(ibc)-1
               if (kqpb(ik) == ib) then
                  active_b = .true.
                  kb = ik
                  exit
               end if
            end do
            
            do ia=isstart(ibr),isstart(ibr)+db(ibr)-1
               ! Is the q.p. with index ia active?
               active_a = .false.
               do ik=kstarta(ibr),kstarta(ibr)+kdima(ibr)-1
                  if (kqpa(ik) == ia) then
                     active_a = .true.
                     ka = ik
                     exit
                  end if
               end do

               im = im + 1
               if (active_a .and. active_b) mat%elem(im) = 1 - fa(ka) - fb(kb)
            end do
         end do
      end do

   end subroutine ft_1fafb_matrix
   
   !----------------------------------------------------------------------------
   ! Compute the finite-temperature statistical factor matrices (f_b - f_a)
   !----------------------------------------------------------------------------
   subroutine ft_fbfa_matrix(ta, tb, mat)
      use constants, only : IT_NEUTRON, IT_PROTON
      use hfbtho_basis, only : pwi_start_n, pwi_start_p, pwi_dim_n, pwi_dim_p, &
         pwi_qp_n, pwi_qp_p, ft_active, ft_fp, ft_fn, isstart
      
      implicit none
      
      integer, intent(in) :: ta, tb
      type(blockmatrix), intent(inout) :: mat
      
      integer :: im, ibr, ibc, ia, ib, ik, ka, kb
      logical :: active_a, active_b
      
      integer,  pointer :: kstarta(:), kstartb(:), kdima(:), kdimb(:), kqpa(:), kqpb(:)
      real(dp), pointer :: fa(:), fb(:)
      
      
      mat%elem(:) = 0
      
      ! Check that finite temperature should be on; if not, exit with the defaults
      if (.not.ft_active) then
         write(*,'(a)') 'WARNING: finite temperature is not active, but ft_1fafb_matrix was called.'
         return
      end if
      
      ! Check isospin integers
      if  ((ta /= IT_NEUTRON .and. ta /= IT_PROTON)   &
      .or. (tb /= IT_NEUTRON .and. tb /= IT_PROTON)) then
         write(*,'(A)') 'ERROR: bad isospin integers in ft_1fafb_matrix.'
         return
      end if
      
      ! Set up pointers
      nullify(kstarta, kdima, kqpa, fa)
      if (ta == IT_NEUTRON) then
         kstarta => pwi_start_n;  kdima => pwi_dim_n;  kqpa => pwi_qp_n;  fa => ft_fn
      else if (ta == IT_PROTON) then
         kstarta => pwi_start_p;  kdima => pwi_dim_p;  kqpa => pwi_qp_p;  fa => ft_fp
      end if
      
      nullify(kstartb, kdimb, kqpb, fb)
      if (tb == IT_NEUTRON) then
         kstartb => pwi_start_n;  kdimb => pwi_dim_n;  kqpb => pwi_qp_n;  fb => ft_fn
      else if (tb == IT_PROTON) then
         kstartb => pwi_start_p;  kdimb => pwi_dim_p;  kqpb => pwi_qp_p;  fb => ft_fp
      end if

      im = 0
      do ibr=1, nb
         
         ibc = mat%ir2c(ibr)
         if (ibc == 0) cycle
         
         do ib=isstart(ibc),isstart(ibc)+db(ibc)-1
            ! Is the q.p. with index ib active?
            active_b = .false.
            do ik=kstartb(ibc),kstartb(ibc)+kdimb(ibc)-1
               if (kqpb(ik) == ib) then
                  active_b = .true.
                  kb = ik
                  exit
               end if
            end do
            
            do ia=isstart(ibr),isstart(ibr)+db(ibr)-1
               ! Is the q.p. with index ia active?
               active_a = .false.
               do ik=kstarta(ibr),kstarta(ibr)+kdima(ibr)-1
                  if (kqpa(ik) == ia) then
                     active_a = .true.
                     ka = ik
                     exit
                  end if
               end do

               im = im + 1
               if (active_a .and. active_b) mat%elem(im) = fb(kb) - fa(ka)
            end do
         end do
      end do
   
   end subroutine ft_fbfa_matrix

end module pnfam_solver
