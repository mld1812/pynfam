!------------------------------------------------------------------------------
!> extfield_2bc.f90
!>
!> Compute the charge-changing external field from two-body weak currents
!> (in single particle configuration space) by contracting the two-body matrix
!> elements of the interaction with the one-body density matrices.
!>
!> E.M. Ney, UNC Chapel Hill 2020
!------------------------------------------------------------------------------
module pnfam_extfield_2bc

   integer,  parameter, private :: dp = kind(1d0), pr = kind(1d0), ipr = kind(1)

    ! 1 = Opt, 2 = Sep, 3 = Notsep, 4 = All
    ! < 0 will print matrix elements, 0 just does opt without print
   integer :: debug = 0

contains

  !> Contract Vr, Vz with the density to calculate the effective one body operator
  !> Fsp_pn = (Gamma_pn, Delta_pn; -Delta_pn*,-Gamma_pn*)
  !
  ! Old Notes: (not sure if they're correct/relevant anymore)
  ! Gamma = J_{abcd} rho_{db} ~ J_{a_c_} --> Like GT+K (a + b = K + c + d)-->(a-c=K)
  ! Delta = J_{acbd} kap_{db} ~ J_{ac__} --> Like ???? (a + c = K + b + d)-->(a+c=K)
  subroutine effective_2bc_extfield(op)
    use type_extfield_2bc
    use type_extfield
    use type_gamdel_2bc
    use hfb_solution
    use type_blockmatrix, only : db, isstart
    use type_blockmatrix, only : allocate_blockmatrix, copy_block_structure, reorder_blockmatrix_basis
    use pnfam_logger

    implicit none

    type(external_field), intent(inout) :: op

    ! Common variables
    integer(ipr) :: ra,rb,rc,rd !> r quantum numbers
    integer(ipr) :: la,lb,lc,ld !> l quantum numbers
    integer(ipr) :: za,zb,zc,zd !> z quantum numbers
    integer(ipr) :: sa,sb,sc,sd !> s quantum numbers
    integer(ipr) :: sac,sdb !> Total spin times 2
    integer :: ita,itb,itc,itd !> State index
    integer :: iba,ibb,ibc,ibd !> block number
    integer :: ida,idb,idc,idd !> Number of states in block
    integer :: isa,isb,isc,isd !> Index of 1st state in block ib, overall
    integer :: ina,inb,inc,ind !> Index of state, within the block
    integer :: nac, ndb !> Column major index within a block
    integer :: img, imd !> Index in vectorized array for 1st element of block
    ! Handle TR states
    integer :: ita_tr,itc_tr
    integer :: a_tr,c_tr, sign_a, sign_c, sign_db
    ! Separable loop variables
    integer :: nzac, npb, npd, nld, nlb, izb, izd, ob, jd, jb, ibbx
    type(extfield_2bc) :: Zrho_n_dir(1:nspatial,0:nrx,0:nrx,0:1,0:1,1:nox+1)
    type(extfield_2bc) :: Zrho_p_dir(1:nspatial,0:nrx,0:nrx,0:1,0:1,1:nox+1)
    type(extfield_2bc) :: Zrho_n_exc(1:nspatial,0:nrx,0:nrx,0:1,0:1,1:nox+1)
    type(extfield_2bc) :: Zrho_p_exc(1:nspatial,0:nrx,0:nrx,0:1,0:1,1:nox+1)
    type(extfield_2bc) :: Zkap_n_dir(1:nspatial,0:nrx,0:nrx,0:1,0:1,1:nox+1)
    type(extfield_2bc) :: Zkap_p_dir(1:nspatial,0:nrx,0:nrx,0:1,0:1,1:nox+1)
    type(extfield_2bc) :: zrho_nd
    type(extfield_2bc) :: zrho_pd
    type(extfield_2bc) :: zrho_ne
    type(extfield_2bc) :: zrho_pe
    type(extfield_2bc) :: zkap_nd
    type(extfield_2bc) :: zkap_pd
    type(extfield_2bc) :: zero_2bc
    ! Special pnfam treatment
    type(extfield_2bc) :: Jz_dir, Jz_exc, Jz_del, Jr_dir, Jr_exc, Jr_del
    logical  :: rblock, kblock, gblock, dblock !> non-zero block indicators
    real(dp) :: rho_p, rho_n, kap_p, kap_n
    real(dp), dimension(6) :: gam, del
    real(dp) :: gam1, del1
    type(gamdel_2bc) :: gam_field, del_field
    integer  :: g, i, nxy
    logical  :: bminus, use_del, use_gam
    ! OpenMP
    integer :: ng, op_k, op_nb
    real(dp), allocatable, dimension(:) :: op_elem, op12_elem
    integer, allocatable, dimension(:) :: op_ir2c, op_ir2m

    ! Desired future behavior for use_2bc:
    ! * I don't want to change operator names anywhere, since 2BC should
    !   just be corrections to 1BC and so should be handled identically
    !   (as far as python is concerned). The treatment of 2BC is entirely
    !   determined by use_2bc value to adjust the extfield upfront.
    ! * 1st digit:
    !      1 = Full FAM, 2 = SNM+LDA, 3 = ASNM+LDA
    ! * 2nd digit:
    !      1 = 1+2Body, 2 = 2Body only
    ! * 3rd digit:
    !      1 = Gamma only, 2 = Delta only, 3 = Gamma+Delta
    ! * 0 = 1Body only
    bminus = op%beta_minus
    use_del = op%use_2bc(3) /= 1
    use_gam = op%use_2bc(3) /= 2

    ! Add additional blockmatrix for pairing field
    ! (Regular GT is (f11, 0; 0, 0), GT2BC has (f11, f12; -f12*, -f11*)
    nxy = size(op%mat%elem)
    if (use_del .and. use_gam) then
       op%mat12%elem = 0
       call allocate_gamdel(op%del, nxy)
       call zero_gamdel(op%del)
       call allocate_gamdel(del_field, nxy)
       call zero_gamdel(del_field)
    end if
    op%mat%elem = 0
    call allocate_gamdel(op%gam, nxy)
    call zero_gamdel(op%gam)
    call allocate_gamdel(gam_field, nxy)
    call zero_gamdel(gam_field)

    ! Populate the arrays needed for calculating Gaussian matrix elements
    call init_extfield_2bc_type(.true.)

    call zero_extfield_2bc(Jz_dir); call zero_extfield_2bc(Jr_dir)
    call zero_extfield_2bc(Jz_exc); call zero_extfield_2bc(Jr_exc)
    call zero_extfield_2bc(Jz_del); call zero_extfield_2bc(Jr_del)
    call zero_extfield_2bc(zero_2bc)
    Jz_dir%k = op%k ; Jr_dir%k = op%k
    Jz_exc%k = op%k ; Jr_exc%k = op%k
    Jz_del%k = op%k ; Jr_del%k = op%k
    zero_2bc%k = op%k

    if (abs(debug)<=1.or.abs(debug)==4) then

       ! All couplings are zero, nothing to do
       if (abs(c4) < 1e-10 .and. abs(c3) < 1e-10 .and. .not. use_p) return

       ! OMP has trouble with allocatable derived types, let's avoid it
       op_nb = size(op%mat%ir2c)

       if (allocated(op_ir2c)) deallocate(op_ir2c, op_ir2m, op_elem)
       allocate(op_ir2c(op_nb), op_ir2m(op_nb), op_elem(nxy))
       op_elem = 0
       op_ir2c = op%mat%ir2c
       op_ir2m = op%mat%ir2m
       op_k = op%k
       ng = nY0
       ! Gam and Del have same block structure (I THINK!) so just need elem
       if (use_del) then
          if (allocated(op12_elem)) deallocate(op12_elem)
          allocate(op12_elem(nxy))
          op12_elem = 0
       end if

!$OMP Parallel Default(NONE)&
!$OMP& SHARED(nzx, zero_2bc, nbx, bminus,&
!$OMP&        nttx, noo, nrr, nll, nss,&
!$OMP&        nzzx, nzz, ib_zrls, i_zrls, db, isstart,&
!$OMP&        rmat, kmat,&
!$OMP&        ng, use_del, use_gam,&
!$OMP&        ntx, hnz, hnr, hnl, hns, ib_itx,&
!$OMP&        op_ir2c, nox, nrx, nrlx,&
!$OMP&        op_k, op_ir2m, op_elem, op12_elem, gam_field, del_field)&
! Privates are declared, uninitialized on each individual thread
!$OMP& PRIVATE(nzac,&
!$OMP&         Zrho_n_dir, Zrho_p_dir, Zrho_n_exc, Zrho_p_exc, Zkap_n_dir, Zkap_p_dir,&
!$OMP&         zc, za,&
!$OMP&         itb, npb, rb, lb, sb, nlb,&
!$OMP&         itd, npd, rd, ld, sd, nld,&
!$OMP&         izb, zb, ibb, jb, idb, isb, ibbx, inb,&
!$OMP&         izd, zd, ibd, jd,      isd,       ind,&
!$OMP&         ndb, rho_n, rho_p, kap_n, kap_p,&
!$OMP&         g, Jz_dir, Jz_exc, Jz_del,&
!$OMP&         ita_tr, ita, a_tr, ra, la, sa, iba, ida, isa, ina, sign_a,&
!$OMP&         itc_tr, itc, c_tr, rc, lc, sc, ibc, idc, isc, inc, sign_c,&
!$OMP&         gblock, gam, del,&
!$OMP&         ob, sdb, sign_db,&
!$OMP&         zrho_nd, zrho_pd, zrho_ne, zrho_pe, zkap_nd, zkap_pd,&
!$OMP&         Jr_dir, Jr_exc, Jr_del,&
!$OMP&         nac, img, imd)

       ! Re-initialized K on each thread
       Jz_dir%k = op_k ; Jr_dir%k = op_k
       Jz_exc%k = op_k ; Jr_exc%k = op_k
       Jz_del%k = op_k ; Jr_del%k = op_k

!$OMP DO SCHEDULE(DYNAMIC)
       ! Loop over za and zc using a single loop, knowing z in [0,nzx]
       ! This ensures the density used in the z contractions matches the one used in the r contraction
       do nzac = 1,(nzx+1)**2
          ! Initialize arrays to have zeroed 2bc types
          Zrho_n_dir = zero_2bc
          Zrho_p_dir = zero_2bc
          Zrho_n_exc = zero_2bc
          Zrho_p_exc = zero_2bc
          Zkap_n_dir = zero_2bc
          Zkap_p_dir = zero_2bc

          ! zc cycles [0,nzx] then za increases by 1, repeat
          zc = mod(nzac-1,nzx+1)   ! Equivalent to: do za=0,nzx
          za = (nzac-zc-1)/(nzx+1) !                   do zc=0,nzx

          ! Loop over (d,b) states with a given (k,r,l,s) < ntx
          do itb = 1,nttx
             npb = noo(itb) ! k = Omega+1/2
             rb = nrr(itb); lb = nll(itb); sb = nss(itb)
             nlb = mod(npb-lb+1,2) ! There are only 2 values of l that give an Omega, l = Omega +/- 1/2
             do itd = 1,nttx
                npd = noo(itd)

                if(npb.ne.npd) cycle ! check kb=kd

                rd = nrr(itd); ld = nll(itd); sd = nss(itd)
                nld = mod(npd-ld+1,2)

                ! Loop over (d,b) states with a given (z,par) for states with a given (k,r,l,s)
                ! nzzx = number of such (z,par) states for i_krls, nzz = nzz(n_krls, n_zpar)
                do izb = 1,nzzx(itb)
                   zb  = nzz(itb,izb)
                   ibb = ib_zrls(zb,rb,lb,(sb+1)/2) ! get standard block from combo of QNs
                   jb  = i_zrls(zb,rb,lb,(sb+1)/2) ! get standard state index from combo of QNs
                   idb = db(ibb)
                   isb = isstart(ibb)-1
                   ibbx= ibb + nbx
                   inb = jb - isb

                   do izd = 1,nzzx(itd)
                      zd  = nzz(itd,izd)
                      ibd = ib_zrls(zd,rd,ld,(sd+1)/2)

                      if(ibb.ne.ibd) cycle ! We checked kb=kd, this checks (kb,parb)=(kd,pard) (for rho AND kap)
                      jd  = i_zrls(zd,rd,ld,(sd+1)/2)
                      isd = isstart(ibd)-1
                      ind = jd - isd

                      ndb = ind + (inb-1)*idb ! Vectorized index within block
                      rho_n = rmat(ndb,ibb)*0.5_dp! rho_db
                      rho_p = rmat(ndb,ibbx)*0.5_dp
                      kap_n = kmat(ndb,ibb) ! -2sb kappa_{d \bar{b}}
                      kap_p = kmat(ndb,ibbx)

                      ! Loop over gaussians and contract over z, populate Z_db(r,l,k)
                      do g = 1, ng
                         if (use_gam) then
                            call calc_Jz_opt(Jz_dir,g,za,zb,zc,zd)
                            if (zc == zd) then
                               Jz_exc = Jz_dir
                            else
                               call calc_Jz_opt(Jz_exc,g,za,zb,zd,zc)
                            end if
                            call add_extfield_2bc(Zrho_n_dir(g,rd,rb,nld,nlb,npb),&
                               scalar_mult_extfield_2bc(Jz_dir, rho_n))
                            call add_extfield_2bc(Zrho_p_dir(g,rd,rb,nld,nlb,npb),&
                               scalar_mult_extfield_2bc(Jz_dir, rho_p))

                            call add_extfield_2bc(Zrho_n_exc(g,rd,rb,nld,nlb,npb),&
                               scalar_mult_extfield_2bc(Jz_exc, rho_n))
                            call add_extfield_2bc(Zrho_p_exc(g,rd,rb,nld,nlb,npb),&
                               scalar_mult_extfield_2bc(Jz_exc, rho_p))
                         end if
                         if (use_del) then
                            call calc_Jz_opt(Jz_del,g,za,zc,zb,zd)
                            call add_extfield_2bc(Zkap_n_dir(g,rd,rb,nld,nlb,npb),&
                               scalar_mult_extfield_2bc(Jz_del, rho_n))
                            call add_extfield_2bc(Zkap_p_dir(g,rd,rb,nld,nlb,npb),&
                               scalar_mult_extfield_2bc(Jz_del, rho_p))
                         end if
                      enddo !ig
                   enddo ! izd
                enddo !izb
             enddo !itd
          enddo !itb

          ! Loop over ALL states (a,c) to populate Gamma_ac (including TR states)
          do ita_tr = 1,ntx*2
             ! Allow us to use the original HFBTHO quantities without explicit TR states
             if (ita_tr > ntx) then
                ita = ita_tr - ntx
                a_tr = -1
             else
                ita = ita_tr
                a_tr = 1
             end if
             if(za.ne.hnz(ita)) cycle
             ! From HFBTHO (no TR states)
             ra = hnr(ita); la = hnl(ita); sa = hns(ita)
             ! Adjust for TR states
             sa = sa*a_tr; la = la*a_tr
             ! From PNFAM (TR states included)
             iba = ib_itx(ita_tr)
             ida = db(iba)
             isa = isstart(iba)-1
             ina = ita_tr-isa
             ! Sign change from spin<0 TR WFs
             sign_a = 1
             if (ita_tr > ntx .and. hns(ita) < 0) sign_a = -1

             do itc_tr = 1,ntx*2
                if (itc_tr > ntx) then
                   itc = itc_tr - ntx
                   c_tr = -1
                else
                   itc = itc_tr
                   c_tr = 1
                end if
                if(zc.ne.hnz(itc)) cycle
                ! From HFBTHO (no TR states)
                rc = hnr(itc); lc = hnl(itc); sc = hns(itc)
                ! Adjust for TR states
                sc = sc*c_tr; lc = lc*c_tr
                ! From PNFAM (TR states included)
                ibc = ib_itx(itc_tr)
                idc = db(ibc)
                isc = isstart(ibc)-1
                inc = itc_tr-isc

                ! Gam and Del have same block structure
                gblock = (op_ir2c(iba) == ibc)
                if (.not.gblock) cycle

                sign_c = 1
                if (itc_tr > ntx .and. hns(itc) < 0) sign_c = -1

                gam = 0.0_dp
                del = 0.0_dp

                ! Loop over (d,b) states (k,l,r,s) for the contraction over r.
                ! There are only 2 ls that give a K, loop over K and 2 ls, which fixes s
                do ob = 0,nox ! nox = max(k-1) = max(Omega - 1/2)
                   ibb = ob + 1 ! This should match npb=k=Omega+1/2, ob=k-1=Omega-1/2
                   do nlb = 0,1
                      lb = ob+nlb
                      sb = -2*nlb+1
                      do nld = 0,1
                         ld = ob+nld

                         ! Currents change nl by K,K+1,K-1 only
                         ! This is already taken care of by block structure??
                         ! This condition is never true in trials...
                         !if(la+lb-lc-ld /= op%k+1 .and.&
                         !   la+lb-lc-ld /= op%k   .and.&
                         !   la+lb-lc-ld /= op%k-1 .and.&
                         !   la-lb-lc+ld /= op%k+1 .and.&
                         !   la-lb-lc+ld /= op%k   .and.&
                         !   la-lb-lc+ld /= op%k-1) then
                         !   cycle
                         !end if

                         sd = -2*nld+1
                         sdb = sd + sb

                         sign_db = 1
                         if (sdb == 0) sign_db = -1

                         do rb = 0,nrx
                            if(2*rb+lb.gt.nrlx) cycle ! Make sure the QNs in the loop match with HO states
                            do rd = 0,nrx
                               if(2*rd+ld.gt.nrlx) cycle

                               ! Loop over gaussians and finish the contraction over (r,l)
                               do g = 1,ng
                                  if (use_gam) then
                                     ! Non TR states
                                     if(la+lb-lc-ld == op_k+1 .or.&
                                        la+lb-lc-ld == op_k   .or.&
                                        la+lb-lc-ld == op_k-1) then
                                        zrho_nd = Zrho_n_dir(g,rd,rb,nld,nlb,ibb)
                                        zrho_pd = Zrho_p_dir(g,rd,rb,nld,nlb,ibb)
                                        zrho_ne = Zrho_n_exc(g,rd,rb,nld,nlb,ibb)
                                        zrho_pe = Zrho_p_exc(g,rd,rb,nld,nlb,ibb)
                                        call calc_Jr_opt(Jr_dir,g,ra,la, rb, lb, rc,lc, rd, ld)
                                        if (rc==rd .and. lc==ld) then
                                           Jr_exc = Jr_dir
                                        else
                                           call calc_Jr_opt(Jr_exc,g,ra,la, rb, lb, rd,ld, rc, lc)
                                        end if
                                        call multiply_extfield_2bc(zrho_nd, Jr_dir)
                                        call multiply_extfield_2bc(zrho_pd, Jr_dir)
                                        call multiply_extfield_2bc(zrho_ne, Jr_exc)
                                        call multiply_extfield_2bc(zrho_pe, Jr_exc)
                                        gam = gam + calc_gam_sep(zrho_nd,zrho_pd,zrho_ne,zrho_pe,sa,sc,sd,sb,bminus)
                                     end if

                                     ! TR states for (d,b) NOTE: (a,c) TR states ARE looped over, (d,b) are NOT
                                     if(la-lb-lc+ld == op_k+1 .or.&
                                        la-lb-lc+ld == op_k   .or.&
                                        la-lb-lc+ld == op_k-1) then
                                        zrho_nd = Zrho_n_dir(g,rd,rb,nld,nlb,ibb)
                                        zrho_pd = Zrho_p_dir(g,rd,rb,nld,nlb,ibb)
                                        zrho_ne = Zrho_n_exc(g,rd,rb,nld,nlb,ibb)
                                        zrho_pe = Zrho_p_exc(g,rd,rb,nld,nlb,ibb)
                                        if ((lb/=0 .or. ld/=0) .or. (la+lb-lc-ld /= op_k+1 .and.&
                                                                     la+lb-lc-ld /= op_k   .and.&
                                                                     la+lb-lc-ld /= op_k+1)) then
                                           call calc_Jr_opt(Jr_dir,g,ra, la, rb,-lb, rc, lc, rd,-ld)
                                        end if
                                        if (rc==rd .and. lc==-ld) then
                                           Jr_exc = Jr_dir
                                        else
                                           call calc_Jr_opt(Jr_exc,g,ra, la, rb,-lb, rd,-ld, rc, lc)
                                        end if
                                        call multiply_extfield_2bc(zrho_nd, Jr_dir)
                                        call multiply_extfield_2bc(zrho_pd, Jr_dir)
                                        call multiply_extfield_2bc(zrho_ne, Jr_exc)
                                        call multiply_extfield_2bc(zrho_pe, Jr_exc)
                                        gam = gam + sign_db*calc_gam_sep(zrho_nd,zrho_pd,zrho_ne,zrho_pe,sa,sc,-sd,-sb,bminus)
                                     end if
                                  end if
                                  if (use_del) then
                                     ! Non TR states
                                     zkap_nd = Zkap_n_dir(g,rd,rb,nld,nlb,ibb)
                                     zkap_pd = Zkap_p_dir(g,rd,rb,nld,nlb,ibb)
                                     call calc_Jr_opt(Jr_del,g,ra,la, rc,-lc, rb,lb, rd,-ld)
                                     call multiply_extfield_2bc(zkap_nd, Jr_del)
                                     call multiply_extfield_2bc(zkap_pd, Jr_del)
                                     del = del + calc_del_sep(zkap_nd,zkap_pd,sa,sc,sd,sb,bminus)
                                     ! TR states (don't think I need to negate kappa cuz its really rho_tilde)
                                     zkap_nd = Zkap_n_dir(g,rd,rb,nld,nlb,ibb)
                                     zkap_pd = Zkap_p_dir(g,rd,rb,nld,nlb,ibb)
                                     call calc_Jr_opt(Jr_del,g,ra,-la, rc, lc, rb,-lb, rd, ld)
                                     call multiply_extfield_2bc(zkap_nd, Jr_del)
                                     call multiply_extfield_2bc(zkap_pd, Jr_del)
                                     del = del + sign_db*calc_del_sep(zkap_nd,zkap_pd,sa,sc,-sd,-sb,bminus)
                                  end if
                               enddo !ig
                            enddo !rd
                         enddo !rb
                      enddo !nld
                   enddo !nlb
                enddo !ob

                nac = ina+(inc-1)*ida
                img = op_ir2m(iba)
                if (use_gam) then
                   gam = gam*sign_a*sign_c
                   op_elem(img-1 + nac) = sum(gam)
                   call fill_gamdel(gam_field, img-1 + nac, gam)
                end if
                if (use_del) then
                   del = del*sign_a*sign_c
                   op12_elem(img-1 + nac) = sum(del)
                   call fill_gamdel(del_field, img-1 + nac, del)
                end if
             enddo !itc
          enddo !ita
       enddo !nzac
!$OMP End Do
!$OMP End Parallel

       if (use_gam .and. use_del) then
          op%mat%elem = op_elem
          op%mat12%elem = op12_elem
          op%gam = gam_field
          op%del = del_field
       else if (use_gam) then
          op%mat%elem = op_elem
          op%gam = gam_field
       else if (use_del) then
          op%mat%elem = op12_elem
          op%del = del_field
       end if
       call deallocate_gamdel(gam_field)
       call deallocate_gamdel(del_field)

       if (debug < 0) then
          print *, "END GT 2BC - Optimized"
          print *, "====================="
       end if

    end if

    if (abs(debug)==2.or.abs(debug)==4) then

       op%mat%elem = 0

       ! Loop over za and zc using a single loop, knowing z in [0,nzx]
       ! This ensures the density used in the z contractions matches the one used in the r contraction
       do nzac = 1,(nzx+1)**2
          ! Initialize arrays to have zeroed 2bc types
          Zrho_n_dir = zero_2bc
          Zrho_p_dir = zero_2bc
          Zrho_n_exc = zero_2bc
          Zrho_p_exc = zero_2bc
          Zkap_n_dir = zero_2bc
          Zkap_p_dir = zero_2bc

          ! zc cycles [0,nzx] then za increases by 1, repeat
          zc = mod(nzac-1,nzx+1)   ! Equivalent to: do za=0,nzx
          za = (nzac-zc-1)/(nzx+1) !                   do zc=0,nzx

          ! Loop over (d,b) states with a given (k,r,l,s) < ntx
          do itb = 1,nttx
             npb = noo(itb) ! k = Omega+1/2
             rb = nrr(itb); lb = nll(itb); sb = nss(itb)
             nlb = mod(npb-lb+1,2) ! There are only 2 values of l that give an Omega, l = Omega +/- 1/2
             do itd = 1,nttx
                npd = noo(itd)

                if(npb.ne.npd) cycle ! check kb=kd

                rd = nrr(itd); ld = nll(itd); sd = nss(itd)
                nld = mod(npd-ld+1,2)

                ! Loop over (d,b) states with a given (z,par) for states with a given (k,r,l,s)
                ! nzzx = number of such (z,par) states for i_krls, nzz = nzz(n_krls, n_zpar)
                do izb = 1,nzzx(itb)
                   zb  = nzz(itb,izb)
                   ibb = ib_zrls(zb,rb,lb,(sb+1)/2) ! get standard block from combo of QNs
                   jb  = i_zrls(zb,rb,lb,(sb+1)/2) ! get standard state index from combo of QNs
                   idb = db(ibb)
                   isb = isstart(ibb)-1
                   ibbx= ibb + nbx
                   inb = jb - isb

                   do izd = 1,nzzx(itd)
                      zd  = nzz(itd,izd)
                      ibd = ib_zrls(zd,rd,ld,(sd+1)/2)

                      if(ibb.ne.ibd) cycle ! We checked kb=kd, this checks (kb,parb)=(kd,pard) (for rho AND kap)
                      jd  = i_zrls(zd,rd,ld,(sd+1)/2)
                      isd = isstart(ibd)-1
                      ind = jd - isd

                      ndb = ind + (inb-1)*idb ! Vectorized index within block
                      rho_n = rmat(ndb,ibb)*0.5_dp! rho_db
                      rho_p = rmat(ndb,ibbx)*0.5_dp
                      kap_n = kmat(ndb,ibb) ! -2sb kappa_{d \bar{b}}
                      kap_p = kmat(ndb,ibbx)

                      ! Loop over gaussians and contract over z, populate Z_db(r,l,k)
                      do g = 1, nspatial
                         call calc_Jz(Jz_dir,g,za,zb,zc,zd)
                         call calc_Jz(Jz_exc,g,za,zb,zd,zc)
                         call calc_Jz(Jz_del,g,za,zc,zb,zd)

                         call add_extfield_2bc(Zrho_n_dir(g,rd,rb,nld,nlb,npb),&
                            scalar_mult_extfield_2bc(Jz_dir, rho_n))
                         call add_extfield_2bc(Zrho_p_dir(g,rd,rb,nld,nlb,npb),&
                            scalar_mult_extfield_2bc(Jz_dir, rho_p))

                         call add_extfield_2bc(Zrho_n_exc(g,rd,rb,nld,nlb,npb),&
                            scalar_mult_extfield_2bc(Jz_exc, rho_n))
                         call add_extfield_2bc(Zrho_p_exc(g,rd,rb,nld,nlb,npb),&
                            scalar_mult_extfield_2bc(Jz_exc, rho_p))

                         call add_extfield_2bc(Zkap_n_dir(g,rd,rb,nld,nlb,npb),&
                            scalar_mult_extfield_2bc(Jz_del, rho_n))
                         call add_extfield_2bc(Zkap_p_dir(g,rd,rb,nld,nlb,npb),&
                            scalar_mult_extfield_2bc(Jz_del, rho_p))
                      enddo !ig
                   enddo ! izd
                enddo !izb
             enddo !itd
          enddo !itb

          ! Loop over ALL states (a,c) to populate Gamma_ac (including TR states)
          do ita_tr = 1,ntx*2
             ! Allow us to use the original HFBTHO quantities without explicit TR states
             if (ita_tr > ntx) then
                ita = ita_tr - ntx
                a_tr = -1
             else
                ita = ita_tr
                a_tr = 1
             end if
             if(za.ne.hnz(ita)) cycle
             ! From HFBTHO (no TR states)
             ra = hnr(ita); la = hnl(ita); sa = hns(ita)
             ! Adjust for TR states
             sa = sa*a_tr; la = la*a_tr
             ! From PNFAM (TR states included)
             iba = ib_itx(ita_tr)
             ida = db(iba)
             isa = isstart(iba)-1
             ina = ita_tr-isa
             ! Sign change from spin<0 TR WFs
             sign_a = 1
             if (ita_tr > ntx .and. hns(ita) < 0) sign_a = -1

             do itc_tr = 1,ntx*2
                if (itc_tr > ntx) then
                   itc = itc_tr - ntx
                   c_tr = -1
                else
                   itc = itc_tr
                   c_tr = 1
                end if
                if(zc.ne.hnz(itc)) cycle
                ! From HFBTHO (no TR states)
                rc = hnr(itc); lc = hnl(itc); sc = hns(itc)
                ! Adjust for TR states
                sc = sc*c_tr; lc = lc*c_tr
                ! From PNFAM (TR states included)
                ibc = ib_itx(itc_tr)
                idc = db(ibc)
                isc = isstart(ibc)-1
                inc = itc_tr-isc

                if (use_del) then
                    dblock = (op%mat12%ir2c(iba) == ibc) ! Are we in a non-zero block?
                else
                    dblock = .false.
                end if
                gblock = (op%mat%ir2c(iba) == ibc) ! Are we in a non-zero block?
                if (.not.gblock .and. .not.dblock) cycle

                sign_c = 1
                if (itc_tr > ntx .and. hns(itc) < 0) sign_c = -1

                gam = 0.0_dp
                del = 0.0_dp

                ! Loop over (d,b) states (k,l,r,s) for the contraction over r.
                ! There are only 2 ls that give a K, loop over K and 2 ls, which fixes s
                do ob = 0,nox ! nox = max(k-1) = max(Omega - 1/2)
                   ibb = ob + 1 ! This should match npb=k=Omega+1/2, ob=k-1=Omega-1/2
                   do nlb = 0,1
                      lb = ob+nlb
                      sb = -2*nlb+1
                      do nld = 0,1
                         ld = ob+nld

                         ! Currents change nl by K,K+1,K-1 only
                         ! This is already taken care of by block structure??
                         ! This condition is never true in trials...
                         !if(la+lb-lc-ld /= op%k+1 .and.&
                         !   la+lb-lc-ld /= op%k   .and.&
                         !   la+lb-lc-ld /= op%k-1 .and.&
                         !   la-lb-lc+ld /= op%k+1 .and.&
                         !   la-lb-lc+ld /= op%k   .and.&
                         !   la-lb-lc+ld /= op%k-1) then
                         !   cycle
                         !end if

                         sd = -2*nld+1
                         sdb = sd + sb

                         sign_db = 1
                         if (sdb == 0) sign_db = -1

                         do rb = 0,nrx
                            if(2*rb+lb.gt.nrlx) cycle ! Make sure the QNs in the loop match with HO states
                            do rd = 0,nrx
                               if(2*rd+ld.gt.nrlx) cycle

                               ! Loop over gaussians and finish the contraction over (r,l)
                               do g = 1,nspatial
                                  if (gblock) then ! always in rblock by construction

                                     ! Non TR states
                                     if(la+lb-lc-ld == op%k+1 .or.&
                                        la+lb-lc-ld == op%k   .or.&
                                        la+lb-lc-ld == op%k-1) then
                                        zrho_nd = Zrho_n_dir(g,rd,rb,nld,nlb,ibb)
                                        zrho_pd = Zrho_p_dir(g,rd,rb,nld,nlb,ibb)
                                        zrho_ne = Zrho_n_exc(g,rd,rb,nld,nlb,ibb)
                                        zrho_pe = Zrho_p_exc(g,rd,rb,nld,nlb,ibb)
                                        call calc_Jr(Jr_dir,g,ra,la, rb, lb, rc,lc, rd, ld)
                                        if (rc==rd .and. lc==ld) then
                                           Jr_exc = Jr_dir
                                        else
                                           call calc_Jr(Jr_exc,g,ra,la, rb, lb, rd,ld, rc, lc)
                                        end if
                                        call multiply_extfield_2bc(zrho_nd, Jr_dir)
                                        call multiply_extfield_2bc(zrho_pd, Jr_dir)
                                        call multiply_extfield_2bc(zrho_ne, Jr_exc)
                                        call multiply_extfield_2bc(zrho_pe, Jr_exc)
                                        gam = gam + calc_gam_sep(zrho_nd,zrho_pd,zrho_ne,zrho_pe,sa,sc,sd,sb,bminus)
                                     end if

                                     ! TR states for (d,b) NOTE: (a,c) TR states ARE looped over, (d,b) are NOT
                                     if(la-lb-lc+ld == op%k+1 .or.&
                                        la-lb-lc+ld == op%k   .or.&
                                        la-lb-lc+ld == op%k-1) then
                                        zrho_nd = Zrho_n_dir(g,rd,rb,nld,nlb,ibb)
                                        zrho_pd = Zrho_p_dir(g,rd,rb,nld,nlb,ibb)
                                        zrho_ne = Zrho_n_exc(g,rd,rb,nld,nlb,ibb)
                                        zrho_pe = Zrho_p_exc(g,rd,rb,nld,nlb,ibb)
                                        if ((lb/=0 .or. ld/=0) .or. (la+lb-lc-ld /= op%k+1 .and.&
                                                                     la+lb-lc-ld /= op%k   .and.&
                                                                     la+lb-lc-ld /= op%k+1)) then
                                           call calc_Jr(Jr_dir,g,ra, la, rb,-lb, rc, lc, rd,-ld)
                                        end if
                                        if (rc==rd .and. lc==-ld) then
                                           Jr_exc = Jr_dir
                                        else
                                           call calc_Jr(Jr_exc,g,ra, la, rb,-lb, rd,-ld, rc, lc)
                                        end if
                                        call multiply_extfield_2bc(zrho_nd, Jr_dir)
                                        call multiply_extfield_2bc(zrho_pd, Jr_dir)
                                        call multiply_extfield_2bc(zrho_ne, Jr_exc)
                                        call multiply_extfield_2bc(zrho_pe, Jr_exc)
                                        gam = gam + sign_db*calc_gam_sep(zrho_nd,zrho_pd,zrho_ne,zrho_pe,sa,sc,-sd,-sb,bminus)
                                     end if

                                  end if
                                  if (dblock .and. use_del) then ! always in kblock by construction
                                     zkap_nd = Zkap_n_dir(g,rd,rb,nld,nlb,ibb)
                                     zkap_pd = Zkap_p_dir(g,rd,rb,nld,nlb,ibb)
                                     call calc_Jr(Jr_del,g,ra,la, rc,-lc, rb,lb, rd,-ld)
                                     call multiply_extfield_2bc(zkap_nd, Jr_del)
                                     call multiply_extfield_2bc(zkap_pd, Jr_del)
                                     del = del + calc_del_sep(zkap_nd,zkap_pd,sa,sc,sd,sb,bminus)
                                     ! TR states
                                     !zkap_nd = scalar_mult_extfield_2bc(Zkap_n_dir(g,rd,rb,nld,nlb,ibb), -1._dp)
                                     !zkap_pd = scalar_mult_extfield_2bc(Zkap_p_dir(g,rd,rb,nld,nlb,ibb), -1._dp)
                                     zkap_nd = Zkap_n_dir(g,rd,rb,nld,nlb,ibb)
                                     zkap_pd = Zkap_p_dir(g,rd,rb,nld,nlb,ibb)
                                     call calc_Jr(Jr_del,g,ra,-la, rc, lc, rb,-lb, rd, ld)
                                     call multiply_extfield_2bc(zkap_nd, Jr_del)
                                     call multiply_extfield_2bc(zkap_pd, Jr_del)
                                     del = del + sign_db*calc_del_sep(zkap_nd,zkap_pd,sa,sc,-sd,-sb,bminus)
                                  end if
                               enddo !ig
                            enddo !rd
                         enddo !rb
                      enddo !nld
                   enddo !nlb
                enddo !ob
                nac = ina+(inc-1)*ida
                if (gblock) then
                   gam = gam*sign_a*sign_c
                   img = op%mat%ir2m(iba)
                   op%mat%elem(img-1 + nac) = sum(gam)
                end if
                if (dblock) then
                   del = del*sign_a*sign_c
                   imd = op%mat12%ir2m(iba)
                   op%mat12%elem(imd-1 + nac) = sum(del)
                end if
                if (debug < 0) then
                   if (abs(sum(gam)) > 1e-16) then
                      nac = ina+(inc-1)*ida
                      img = op%mat%ir2m(iba)
                      write(*,'(L2,2x,3i4,2x,1e30.20, " |O: ", 2i3, " |L: ", 2i3, " |s: ", 2i3, " |r: ", 2i3, " |z: ", 2i3," |p: "&
                      &,2i3," |diff_O,L,S,P: ",4i3)') &
                      &(op%mat%ir2c(iba) == ibc),ita_tr, itc_tr, img-1+nac, sum(gam), 2*la+sa, 2*lc+sc, la, lc, sa, sc, ra,rc, za,zc, &
                      &npar(ita_tr), npar(itc_tr),  2*la+sa - (2*lc+sc), la-lc, sa-sc, npar(ita_tr)-npar(itc_tr)
                   end if
                end if
             enddo !itc
          enddo !ita
       enddo !nzac

       if (debug < 0) then
          print *, "END GT 2BC - Separable"
          print *, "====================="
       end if

    end if
    if (abs(debug)==3.or.abs(debug)==4) then

       op%mat%elem = 0

       ! Run over (a,c) to populate gamma_ac and delta_ac (TR states as well)
       do ita_tr = 1,ntx*2
          if (ita_tr > ntx) then
             ita = ita_tr - ntx
             a_tr = -1
          else
             ita = ita_tr
             a_tr = 1
          end if
          ra = hnr(ita); za = hnz(ita); la = hnl(ita); sa = hns(ita)
          sa = sa*a_tr; la = la*a_tr
          iba = ib_itx(ita_tr)
          ida = db(iba)
          isa = isstart(iba)-1
          ina = ita_tr-isa

          sign_a = 1
          if (ita_tr > ntx .and. hns(ita) < 0) sign_a = -1

          do itc_tr = 1,ntx*2
             if (itc_tr > ntx) then
                itc = itc_tr - ntx
                c_tr = -1
             else
                itc = itc_tr
                c_tr = 1
             end if
             rc = hnr(itc); zc = hnz(itc); lc = hnl(itc); sc = hns(itc)
             sc = sc*c_tr; lc = lc*c_tr
             ibc = ib_itx(itc_tr)
             idc = db(ibc)
             isc = isstart(ibc)-1
             inc = itc_tr-isc

             ! Setting dblock = False should skip everything related to del
             ! except actually calculating the matrix elements.
             if (use_del) then
                 dblock = (op%mat12%ir2c(iba) == ibc) ! Are we in a non-zero block?
             else
                 dblock = .false.
             end if
             gblock = (op%mat%ir2c(iba) == ibc) ! Are we in a non-zero block?
             if (.not.gblock .and. .not.dblock) cycle

             sign_c = 1
             if (itc_tr > ntx .and. hns(itc) < 0) sign_c = -1

             gam1 = 0
             del1 = 0

             ! Run over (b,d) to contract with density_db (Omega > 0)
             do itb = 1,ntx
                rb = hnr(itb); zb = hnz(itb); lb = hnl(itb); sb = hns(itb)
                ibb = ib_itx(itb)
                idb = db(ibb)
                isb = isstart(ibb)-1
                inb = itb-isb
                do itd = 1,ntx
                   rd = hnr(itd); zd = hnz(itd); ld = hnl(itd); sd = hns(itd)
                   ibd = ib_itx(itd)
                   idd = db(ibd)
                   isd = isstart(ibd)-1
                   ind = itd-isd

                   sdb = sb+sd

                   if (ibd > nbx) then
                      kblock = (ibd-nbx == ibb) ! kap is like V
                   else
                      kblock = (ibd+nbx == ibb) ! kap is like V
                   end if
                   rblock = (ibd == ibb) ! rho is block diagonal
                   if (.not.rblock .and. .not.kblock) cycle

                   ! Column major counter WITHIN the block
                   ! Note: idb and idd are equal for kap, since it connects blocks with same |Omega|
                   ! Note: rk = 2*V'V, the 2 due to implicit TR ("m-projection")
                   ! Note: ak = 1/2(-U'V), the 1/2 due to "symmetrization"??, the minus for rho_tilde
                   ndb = ind + (inb-1)*idd

                   rho_n = rmat(ndb,ibd)*0.5_dp! rho_db
                   rho_p = rmat(ndb,ibd+nbx)*0.5_dp
                   kap_n = kmat(ndb,ibd) ! -2sb kappa_{d \bar{b}}
                   kap_p = kmat(ndb,ibd+nbx)

                   sign_db = 1
                   if (sdb == 0) sign_db = -1

                   do g = 1,nspatial
                      if (gblock .and. rblock) then
                         call calc_Jz(Jz_dir,g,za,zb,zc,zd)
                         call calc_Jz(Jz_exc,g,za,zb,zd,zc)
                         call calc_Jr(Jr_dir,g,ra,la, rb, lb, rc,lc, rd, ld)
                         call calc_Jr(Jr_exc,g,ra,la, rb, lb, rd,ld, rc, lc)
                         call multiply_extfield_2bc(Jz_dir, Jr_dir)
                         call multiply_extfield_2bc(Jz_exc, Jr_exc)
                         gam1 = gam1 + calc_gam(Jz_dir,Jz_exc,sa,sc,sd,sb,rho_p,rho_n,bminus)
                         ! TR states for (d,b) NOTE: (a,c) TR states ARE looped over, (d,b) are NOT
                         call calc_Jz(Jz_dir,g,za,zb,zc,zd)
                         call calc_Jz(Jz_exc,g,za,zb,zd,zc)
                         call calc_Jr(Jr_dir,g,ra,la,rb,-lb,rc, lc,rd,-ld)
                         call calc_Jr(Jr_exc,g,ra,la,rb,-lb,rd,-ld,rc, lc)
                         call multiply_extfield_2bc(Jz_dir, Jr_dir)
                         call multiply_extfield_2bc(Jz_exc, Jr_exc)
                         gam1 = gam1 + sign_db*calc_gam(Jz_dir,Jz_exc,sa,sc,-sd,-sb,rho_p,rho_n,bminus)
                      end if
                      if (dblock .and. kblock .and. use_del) then
                         call calc_Jz(Jz_del,g,za,zc,zb,zd)
                         call calc_Jr(Jr_del,g,ra,la, rc,-lc, rb,lb, rd,-ld)
                         call multiply_extfield_2bc(Jz_del, Jr_del)
                         del1 = del1 + calc_del(Jz_del,sa,sc,sd,sb,kap_p,kap_n,bminus)
                         ! TR states for (d,b) NOTE: (a,c) TR states ARE looped over, (d,b) are NOT
                         call calc_Jz(Jz_del,g,za,zc,zb,zd)
                         call calc_Jr(Jr_del,g,ra,la, rc,-lc, rb,-lb, rd,ld)
                         call multiply_extfield_2bc(Jz_del, Jr_del)
                         del1 = del1 + sign_db*calc_del(Jz_del,sa,sc,-sd,-sb,kap_p,kap_n,bminus)
                      end if
                   end do ! spatial terms g
                end do ! d (row)
             end do ! b (col)
             ! inc is block state index, from 1 to ndc
             nac = ina+(inc-1)*ida
             if (gblock) then
                gam1 = gam1*sign_a*sign_c
                img = op%mat%ir2m(iba)
                op%mat%elem(img-1 + nac) = gam1
             end if
             if (dblock) then
                del1 = del1*sign_a*sign_c
                imd = op%mat12%ir2m(iba)
                op%mat12%elem(imd-1 + nac) = del1
             end if

             if (debug < 0) then
                if (abs(gam1) > 1e-16) then
                   nac = ina+(inc-1)*ida
                   img = op%mat%ir2m(iba)
                   write(*,'(L2,2x,3i4,2x,1e30.20, " |O: ", 2i3, " |L: ", 2i3, " |s: ", 2i3, " |r: ", 2i3, " |z: ", 2i3," |p: "&
                   &,2i3," |diff_O,L,S,P: ",4i3)') &
                   &(op%mat%ir2c(iba) == ibc),ita_tr, itc_tr, img-1+nac, gam1, 2*la+sa, 2*lc+sc, la, lc, sa, sc, ra,rc, za,zc, &
                   &npar(ita_tr), npar(itc_tr),  2*la+sa - (2*lc+sc), la-lc, sa-sc, npar(ita_tr)-npar(itc_tr)
                end if
             end if
          end do ! c (col)
       end do ! a (row)

       if (debug < 0) then
          print *, "END GT 2BC - Not Separable"
          print *, "====================="
       end if

    end if ! DEBUG OR NOT

#ifdef USE_HBLAS
#if USE_HBLAS==1
   ! Reorder rows AND cols to match new basis order
   call reorder_blockmatrix_basis(op%mat, new_order, 'a')
   call reorder_gamdel_basis(op%gam, op%mat, new_order, 'a')
   if (use_del) then
      call reorder_blockmatrix_basis(op%mat12, new_order, 'a')
      call reorder_gamdel_basis(op%del, op%mat, new_order, 'a')
   end if
#endif
#endif

  end subroutine effective_2bc_extfield

end module pnfam_extfield_2bc
