module modbulkmicro3_point
  use modprecision, only : field_r
  use modglobal, only :    rdt,rk3step
  use modmicrodata, only : Dv,Kt,delt,nu_a,Sc_num,D_eq, &
                           k_br,k_l,k_r,k_rr,kappa_r, phi, pirhow, &
                           k_1,k_2,k_au,k_c, &
                           x_s,xrmax,xrmin,&
                           l_sb,l_mur_cst,mur_cst,l_rain, &
                           eps0, qcmin
  use modmicrodata3
  implicit none

  private
  public point_processes

  real ::  tmp0     &
          ,ql0      &
          ,qt0      &
          ,esl      &
          ,qvsl     &
          ,qvsi     &
          ,cp_exnf_k &
          ,rhof_k    &
          ,presf_k   &

          ,n_cc, n_ccp, n_ccm  & ! N_{ccn} nr content [ kg^{-1}] of cloud condensation nuclei
          ,n_cl, n_clp, n_clm  & ! N_{c,l} nr content [ kg^{-1}] for liquid cloud droplets,
          ,n_ci, n_cip, n_cim  & ! N_{c,i} nr content [ kg^{-1}] for ice cloud droplets,
          ,n_hr, n_hrp, n_hrm  & ! N_{h,r} nr content [ kg^{-1}] for rain
          ,n_hs, n_hsp, n_hsm  & ! N_{h,s} nr content [ kg^{-1}] for snow
          ,n_hg, n_hgp, n_hgm  & ! N_{h,g} nr content [ kg^{-1}] for graupel
          ,q_cl, q_clp, q_clm  & ! q_{c,l} water content [kg/kg] for liquid cloud droplets,
          ,q_ci, q_cip, q_cim  & ! q_{c,i} water content [kg/kg] for ice cloud droplets,
          ,q_hr, q_hrp, q_hrm  & ! q_{h,r} water content [kg/kg] for rain
          ,q_hs, q_hsp, q_hsm  & ! q_{h,s} water content [kg/kg] for snow
          ,q_hg, q_hgp, q_hgm  & ! q_{h,g} water content [kg/kg] for graupel

          ,x_ci     & ! mean cloud ice size
          ,x_cl     & ! mean cloud water size
          ,x_hs     & ! mean snow size
          ,x_hg     & ! mean graupel size
          ,x_hr     & ! mean raindrop size
          ,D_ci     & ! \bar{D}_ci mean diameter for cloud ice particle
          ,D_cl     & ! \bar{D}_cl mean diameter for cloud water particle
          ,D_hr     & ! \bar{D}_hr mean diameter for raindrops
          ,D_hs     & ! \bar{D}_hs mean diameter for snow particle
          ,D_hg     & ! \bar{D}_hg mean diameter for graupel particle
          ,v_ci     & ! \bar{v}_ci mean velocity for cloud ice particle
          ,v_cl     & ! \bar{v}_cl mean velocity for cloud water droplets
          ,v_hr     & ! \bar{v}_hr mean velocity for raindrops
          ,v_hs     & ! \bar{v}_hs mean velocity for snow particle
          ,v_hg     & ! \bar{v}_hg mean velocity for graupel particle
          ,Dvr,lbdr &

          ,dn_cl_au       &    !< change in number of cloud droplets due to autoconversion
          ,dn_cl_ac       &    !< change in number of cloud droplets due to accretion
          ,dq_ci_dep      &    !< deposition rate for clouds
          ,dq_hs_dep      &    !< deposition rate for snow
          ,dq_hg_dep      &    !< deposition rate for graupel
          ,dq_ci_rime     &    !< riming growth of ice
          ,dn_cl_rime_ci  &    !<  - and impact on n_cl
          ,dq_hs_rime     &    !< riming growth of snow
          ,dn_cl_rime_hs  &    !<  - and impact on n_cl
          ,dq_hg_rime     &    !< riming growth for graupel
          ,dn_cl_rime_hg  &    !<  - and impact on n_cl
          ,dq_hghr_rime   &    !< riming growth for graupel with rain
          ,dq_hr_col_rs   &    !< rain loss from riming of ice+snow->gr
          ,dq_hs_col_rs   &    !< rain number loss from riming of ice+snow->gr
          ,dq_hr_col_ri   &    !< rain loss from riming of ice+rain->gr
          ,dq_ci_col_ri   &    !< ice loss from riming of ice+rain->gr
          ,dn_cl_het      &    !< heterogeneou freezing of cloud water
          ,dn_cl_hom      &    !< homogeneous freezing of cloud water
          ,dn_cl_sc       &    !< cloud self-collection
          ,ret_cc              !< recovery of ccn

  real ::   qtpmcr                    &
           ,thlpmcr

  real, allocatable :: statistics   (:)     &
                      ,tend         (:)

contains
! TODO:  * removed warnings when throwing away too much negative moisture

! NOTES: deposition tendencies are the tendencies after correction from cor_deposit3

  subroutine point_processes(prg,exnf_k_in,rhof_k_in,presf_k_in   &
                            ,sv0,svp,svm,thlpmcr_out,qtpmcr_out   &
                            ,statistics_out,tend_out  )

    use modglobal, only     : cp
    use modmicrodata3, only : in_hr,iq_hr,in_cl,iq_cl,in_cc, &
                              in_ci,iq_ci,in_hs,iq_hs,in_hg,iq_hg
    implicit none
    real(field_r), intent(in)    :: exnf_k_in, rhof_k_in,presf_k_in
    real, intent(in)    :: sv0(ncols),svm(ncols),prg(nprgs)
    real, intent(inout) :: svp(ncols)
    real, intent(out)   :: thlpmcr_out,qtpmcr_out
    real, intent(out)   :: tend_out(ntends),statistics_out(nmphys)

    if (l_statistics) then
      allocate(statistics(nmphys))
    endif
    if (l_tendencies) then
      allocate(tend(ntends))
    endif

    ! base variables
    tmp0 = prg(n_tmp0)
    qt0  = prg(n_qt0)
    ql0  = prg(n_ql0)
    esl  = prg(n_esl)
    qvsl = prg(n_qvsl)
    qvsi = prg(n_qvsi)

    cp_exnf_k = exnf_k_in * cp
    rhof_k = rhof_k_in
    presf_k = presf_k_in

    ! scalars
    n_cc = sv0(in_cc)
    n_cl = sv0(in_cl)
    n_ci = sv0(in_ci)
    n_hr = sv0(in_hr)
    n_hs = sv0(in_hs)
    n_hg = sv0(in_hg)
    q_cl = sv0(iq_cl)
    q_ci = sv0(iq_ci)
    q_hr = sv0(iq_hr)
    q_hs = sv0(iq_hs)
    q_hg = sv0(iq_hg)

    n_ccm = svm(in_cc)
    n_clm = svm(in_cl)
    n_cim = svm(in_ci)
    n_hrm = svm(in_hr)
    n_hsm = svm(in_hs)
    n_hgm = svm(in_hg)
    q_clm = svm(iq_cl)
    q_cim = svm(iq_ci)
    q_hrm = svm(iq_hr)
    q_hsm = svm(iq_hs)
    q_hgm = svm(iq_hg)

    n_ccp  = svp(in_cc)
    n_clp  = svp(in_cl)
    n_cip  = svp(in_ci)
    n_hrp  = svp(in_hr)
    n_hsp  = svp(in_hs)
    n_hgp  = svp(in_hg)
    q_clp  = svp(iq_cl)
    q_cip  = svp(iq_ci)
    q_hrp  = svp(iq_hr)
    q_hsp  = svp(iq_hs)
    q_hgp  = svp(iq_hg)


    ! Reset all values
    ! ----------------

    ! 0 bulk integrals and parameters
    Dvr = 0.
    lbdr = 0.

    x_ci   = 0.0
    x_cl   = 0.0
    x_hr   = 0.0
    x_hs   = 0.0
    x_hg   = 0.0

    D_ci   = 0.0
    D_cl   = 0.0
    D_hr   = 0.0
    D_hs   = 0.0
    D_hg   = 0.0

    v_ci   = 0.0
    v_cl   = 0.0
    v_hr   = 0.0
    v_hs   = 0.0
    v_hg   = 0.0

    ! Initialize the (remaining) global variables
    ! those are calculated in one routines, and used by another.
    ! If we need to run the microphysics with OpenMP, we need to
    ! replace these variables with local variables+arguments.
    ! variable:           used by:
    dn_cl_au       = 0. ! autoconversion3, cloud_self3, recover_cc
    dn_cl_ac       = 0. ! accretion3, recover_cc
    dq_ci_dep      = 0. ! deposit_ice3, cor_deposit3
    dq_hs_dep      = 0. ! deposit_snow3, cor_deposit3
    dq_hg_dep      = 0. ! deposit_graupel3, cor_deposit3
    dq_ci_rime     = 0. ! coll_ici3, ice_multi3, conv_partial3
    dn_cl_rime_ci  = 0. ! coll_ici3, recover_cc
    dq_hs_rime     = 0. ! coll_scs3, ice_multi3, conv_partial3
    dn_cl_rime_hs  = 0. ! coll_scs3, recover_cc
    dq_hg_rime     = 0. ! coll_gcg3, ice_multi3
    dn_cl_rime_hg  = 0. ! coll_gcg3, recover_cc
    dq_hghr_rime   = 0. ! coll_grg3, ice_multi3
    dq_hr_col_rs   = 0. ! coll_rsg3, ice_multi3
    dq_hs_col_rs   = 0. ! coll_rsg3, ice_multi3
    dq_hr_col_ri   = 0. ! coll_rig3, ice_multi3
    dq_ci_col_ri   = 0. ! coll_rig3, ice_multi3
    dn_cl_het      = 0. ! hetfreez3, recover_cc
    dn_cl_hom      = 0. ! homfreez3, recover_cc
    dn_cl_sc       = 0. ! cloud_self3, recover_cc
    ret_cc         = 0. ! evap_rain3, recover_cc

    thlpmcr = 0.
    qtpmcr = 0.


    ! Testing for a noticable amount of rain graupel and snow
    ! -------------------------------------------------------

    ! rain :
    q_hr_mask = (q_hr.gt.q_hr_min).and.(n_hr.gt.0.0)

    ! snow :
    q_hs_mask = (q_hs.gt.q_hs_min).and.(n_hs.gt.0.0)

    ! graupel :
    q_hg_mask = (q_hg.gt.q_hg_min).and.(n_hg.gt.0.0)

    ! liquid clouds
    q_cl_mask = (q_cl.ge.qcliqmin).and.(n_cl.gt.0.0)

    ! ice clouds
    q_ci_mask = (q_ci.ge.qicemin).and.(n_ci.gt.0.0)


    ! calculate Rain DSD integral properties & parameters lbdr
    ! -----------------------------------------------------------------
    call integrals_bulk3

    ! NOTE: not using cloud mask condition, since air can be saturated with respect to ice
    if (tmp0.lt.tmp_inuc) call icenucle3  ! ice nucleation

    ! freezing of water droplets
    ! -----------------------------------------------------------------
    if (tmp0.le.T_3) then
      if (q_cl_mask) then
        call homfreez3   ! homogeneous freezing of cloud droplets
        call hetfreez3   ! heterogeneous freezing
      endif

    ! deposition processes
    ! -----------------------------------------------------------------
      if (q_ci_mask) call deposit_ice3      ! deposition of vapour to cloud ice
      if (q_hs_mask) call deposit_snow3     ! deposition of vapour to snow
      if (q_hg_mask) call deposit_graupel3  ! deposition of vapour to graupel
      call cor_deposit3                     ! correction for deposition
    endif

    ! snow aggregation and self-collection
    ! -----------------------------------------------------------------
    if(q_ci_mask) call ice_aggr3         ! ice selfcollection
    if(q_hs_mask) call snow_self3        ! snow selfcollection - tendency only

    ! collision processes for snow
    ! -----------------------------------------------------------------
    if (q_hs_mask) then
      if (q_ci_mask) call coll_sis3              ! snow selfcollection s+i -> s
      if (q_hg_mask) call coll_gsg3              ! snow selfcollection g+s -> g
    endif
    if (q_cl_mask) then
      if (q_ci_mask) then
        if ((D_cl.gt.D_c_a).and.(D_ci.gt.D_i0)) then
          call coll_ici3                         ! riming i+c -> i
        endif
      endif
      if (q_hs_mask) then
        if ((D_cl.gt.D_c_a).and.(D_hs.gt.D_i0)) then
          call coll_scs3                         ! riming s+c -> s
        endif
      endif
      if (q_hg_mask) then
        if((D_cl.gt.D_c_a).and.(D_hg.gt.D_i0)) then
          call coll_gcg3                         ! riming g+c -> g
        endif
      endif
    endif
    if (q_hg_mask.and.q_hr_mask) call coll_grg3  ! riming g+r -> g

    ! raindrop freezing
    ! -----------------------------------------------------------------
    if (q_hr_mask) then
      if (tmp0.lt.T_3) call rainhetfreez3       ! heterogeneous freezing

    ! collision with conversion
    ! -----------------------------------------------------------------
      if (q_ci_mask) call coll_rig3         ! riming r+i -> g
      if (q_hs_mask) call coll_rsg3         ! riming r+i -> g
    endif

    ! conversions
    ! -----------------------------------------------------------------
    if (tmp0.lt.T_3) then
      call ice_multi3                                      ! ice multiplication of Hallet and Mossop
      if (l_sb_conv_par.and.q_cl_mask) call conv_partial3  ! partial conversion
    else ! tmp0.gt.T_3

    ! melting and evaporation of ice particles
    ! -----------------------------------------------------------------
      call evapmelting3                                    ! melting of ice particles
    endif

    ! basic warm processes
    ! -----------------------------------------------------------------
    if (q_cl_mask) then
      call autoconversion3
      call cloud_self3
    endif
    if (q_hr_mask) then
      call accretion3
      call evap_rain3
    endif

    ! saturation adjustment
    ! -----------------------------------------------------------------
    call satadj3

    ! Recover ccn (except for sedimentation part)
    if(.not.(l_c_ccn)) call recover_cc

    ! Save the results
    ! -----------------------------------------------------------------
    svp(in_cc) = n_ccp
    svp(in_cl) = n_clp
    svp(in_ci) = n_cip
    svp(in_hr) = n_hrp
    svp(in_hs) = n_hsp
    svp(in_hg) = n_hgp
    svp(iq_cl) = q_clp
    svp(iq_ci) = q_cip
    svp(iq_hr) = q_hrp
    svp(iq_hs) = q_hsp
    svp(iq_hg) = q_hgp

    thlpmcr_out = thlpmcr
    qtpmcr_out = qtpmcr

    if (l_statistics) then
      statistics_out = statistics
      deallocate(statistics)
    endif
    if (l_tendencies) then
      tend_out = tend
      deallocate(tend)
    endif
end subroutine point_processes


! calculating rain and other integrals
subroutine integrals_bulk3
  implicit none

  real :: N_r0, mur

  if (q_hr_mask) then
    ! NOTE: code is duplicated in subroutine sedim_rain3
    if (l_rain) then
      ! limiting procedure (as per S&B)
      x_hr    = q_hr/n_hr
      x_hr    = max(xrmin,min(xrmax,x_hr))

      if (l_sb) then
        if(l_sb_classic) then
          D_hr = a_hr *x_hr**b_hr
          v_hr = al_hr*x_hr**be_hr*(rho0s/rhof_k)**ga_hr

          Dvr  = (x_hr/pirhow)**(1./3.)
          N_r0 = rhof_k*n_hr/Dvr
          N_r0 = max(N_0min,min(N_0max,N_r0))

          lbdr = (pirhow*N_r0/(rhof_k*q_hr))**0.25  ! c_lbdr*x_hr**(-mu_hr_cst)
          lbdr = max(lbdr_min, min(lbdr_max,lbdr))
        else ! l_sb_classic
          Dvr = (x_hr/pirhow)**(1./3.)

          D_hr = a_hr *x_hr**b_hr
          v_hr = al_hr*x_hr**be_hr*(rho0s/rhof_k)**ga_hr

          if (l_mur_cst) then
            ! mur = cst
            mur  = mur_cst
          else
            ! mur = f(Dv)
            mur  = min(mur0_G09b,- 1. + c_G09b/ (q_hr*rhof_k)**exp_G09b)  ! G09b
          endif ! l_mur_cst
          lbdr = ((mur+3.)*(mur+2.)*(mur+1.))**(1./3.)/Dvr
        endif !  l_sb_classic
      else !  l_sb
        Dvr = (x_hr/pirhow)**(1./3.)

        D_hr = a_hr *x_hr**b_hr
        v_hr = al_hr*x_hr**be_hr*(rho0s/rhof_k)**ga_hr
      endif ! l_sb
    endif  ! l_rain
  endif ! q_hr_mask

  ! cloud water
  if (q_cl_mask) then
    x_cl = q_cl/n_cl
    x_cl = min(max(x_cl,x_cl_bmin),x_cl_bmax) ! as well - to limit extrme autoconversion

    D_cl = a_cl *x_cl**b_cl
    v_cl = al_cl*x_cl**be_cl*(rho0s/rhof_k)**ga_cl
  endif

  ! cloud ice
  if (q_ci_mask) then
    x_ci = q_ci/n_ci
    x_ci = min(max(x_ci,x_ci_bmin),x_ci_bmax) ! to ensure x is within borders

    D_ci = a_ci *x_ci**b_ci
    v_ci = al_ci*x_ci**be_ci*(rho0s/rhof_k)**ga_ci
  endif

  ! snow
  if (q_hs_mask) then
    x_hs = q_hs/n_hs
    x_hs = min(max(x_hs,x_hs_bmin),x_hs_bmax) ! to ensure x is within borders

    D_hs = a_hs *x_hs**b_hs
    v_hs = al_hs*x_hs**be_hs*(rho0s/rhof_k)**ga_hs
  endif

  ! graupel
  if (q_hg_mask) then
    x_hg = q_hg/n_hg
    x_hg = min(max(x_hg,x_hg_bmin),x_hg_bmax) ! to ensure x is within borders

    D_hg = a_hg *x_hg**b_hg
    v_hg = al_hg*x_hg**be_hg*(rho0s/rhof_k)**ga_hg
  endif
end subroutine integrals_bulk3


! ****************************************************************
! Ice nucleation
!  - based on Seifert & Beheng (2003), p. 53
! ****************************************************************
subroutine icenucle3
  implicit none

  real :: ssice
  real :: n_in,n_tid
  real :: dn_ci_inu    !< ice nucleation rate

  dn_ci_inu = 0.

  ! not always how supersaturated is water vapour, depending on a flag,
  ! it can also include water already in ice particles
  ssice = 0.0

  ! calculate supersaturation with respect to ice
  if (l_sb_inuc_sat) then  ! l_sb_inuc_sat
    ! calculating supersaturation of water vapour only
    ssice = (qt0-q_cl)/qvsi -1.0
    ! ssice = max(0.0, (qt0-q_cl)/qvsi -1.0)
  else  ! l_sb_inuc_sat
    ! ie. cloud ice water is also included supersaturation
    ssice = (qt0-q_cl+q_ci)/qvsi -1.0
    ! ssice = max(0.0, (qt0-q_cl+q_ci)/qvsi -1.0)
  endif ! l_sb_inuc_sat

  if (ssice.gt.ssice_min) then ! condition for nucleation
    if (l_sb_inuc_expl) then ! l_sb_inuc_expl
      ! explicit ice nucleation - not yet included
      dn_ci_inu = 0. ! BUG: should give not-implemented error
    else  ! l_sb_inuc_expl
      ! Meyers et al. (1992)
      n_in = (1.0/rhof_k)*N_inuc*exp( a_M92              &
          +b_M92*min(ssice,ssice_lim))! o: N_inuc*exp( a_M92 + b_M92*ssice)

      ! whether to apply reisnerr correction
      if (l_sb_reisner) then
        ! prepare Reisnerr (1998) correction
        ! w: n_tid = (1.0/rhof_k)*N_inuc_R*exp( - min(tmp0,c_inuc_R)-T_3)
        n_tid = (1.0/rhof_k)*N_inuc_R*exp(b_inuc_R*(T_3- max(tmp0,c_inuc_R)))

        ! performing reisner correction
        n_in = max(a1_inuc_R*n_tid,min(a2_inuc_R*n_tid,n_in))
      endif ! l_sb_reisner

      ! limiting n_in
      n_in = min(n_i_max, n_in)

      ! checking conditions if suitable for nucleation
      if (n_ci.lt.n_in) then ! condition intentionally left this way
        dn_ci_inu = (n_in-n_ci)/delt
      else
        dn_ci_inu = 0.0
        ! note - written this way on purpose
        ! in case of further adjustment of nucleation subroutine
      endif
    endif ! l_sb_inuc_expl

    ! basic correction  - not suitable

    ! update cloud water number
    n_cip = n_cip + dn_ci_inu

    ! update water density [kg kg^{-1}]
    q_cip = q_cip + x_inuc*dn_ci_inu

    ! update liquid water potential temperature
    !  - due to latent heat of melting (freezing in this case)
    qtpmcr = qtpmcr - x_inuc*dn_ci_inu
    thlpmcr = thlpmcr+(rlvi/cp_exnf_k)*x_inuc*dn_ci_inu

    if (l_sb_dbg) then
      if ((qt0-qvsi-q_clm-q_cim)/delt - x_inuc*dn_ci_inu.lt. 0.) then
        write(6,*) 'WARNING: icenucle3 removing too much water' &
                  ,x_inuc*dn_ci_inu, '/', (qt0-qvsi-q_clm-q_cim)/delt
        ! count((qt0-qvsi-q_clm-q_cim)/delt - x_inuc*dn_ci_inu.lt. 0.)
      endif
    endif ! l_sb_dbg
  endif ! ssice.gt.ssice_min

  if (l_tendencies) then
    tend(idn_ci_inu) = dn_ci_inu
  endif
  if (l_statistics) then
    statistics(imphys_dep) = statistics(imphys_dep) + x_inuc*dn_ci_inu
  endif
end subroutine icenucle3


! ****************************************
! Homogeneous freezing
! of cloud droplets
!
! ****************************************
subroutine homfreez3
  use modglobal, only : rlv
  implicit none

  real    :: J_hom   ! freezing rate
  real    :: tmpj    ! adjusted temperature

  ! calculate constants
  real, parameter :: expC_30 = exp(C_30_Cf02)

  real :: dq_cl_hom = 0. !< homogeneous freezing of cloud water

  ! inserting temperature values
  !  -- can be later adjusted to include the effect of chemicals in clouds
  tmpj = tmp0

  if (tmpj>tmp_lim1_Cf02) then
    J_hom = exp(C_Cf02+B_Cf02*(tmpj+CC_Cf02))
  else if(tmpj<tmp_lim2_Cf02) then
    J_hom = expC_30
  else
    J_hom = exp(C_20_Cf02                 &
       + B_21_Cf02*(tmpj-offset_Cf02)     &
       + B_22_Cf02*(tmpj-offset_Cf02)**2  &
       + B_23_Cf02*(tmpj-offset_Cf02)**3  &
       + B_24_Cf02*(tmpj-offset_Cf02)**4  )
  endif

  dn_cl_hom        = -c_mmt_1cl * n_cl * x_cl * J_hom
  dq_cl_hom        = -c_mmt_2cl * q_cl * x_cl * J_hom
  dn_cl_hom        = max(dn_cl_hom,min(0.0,-n_clm/delt-n_clp))
  dq_cl_hom        = max(dq_cl_hom,min(0.0,-q_clm/delt-q_clp))

  ! changes in cloud water
  n_clp = n_clp + dn_cl_hom
  q_clp = q_clp + dq_cl_hom

  ! increase in cloud ice
  n_cip = n_cip - dn_cl_hom
  q_cip = q_cip - dq_cl_hom

  ! changes in the total amount of water
  qtpmcr = qtpmcr + dq_cl_hom

  ! change in th_l due to freezing
  thlpmcr = thlpmcr-((rlv+rlme)/cp_exnf_k)*dq_cl_hom
  if (l_tendencies) then
    tend(idq_cl_hom) = dq_cl_hom
    tend(idn_cl_hom) = dn_cl_hom
  endif
  if (l_statistics) then
    statistics(imphys_freeze) = statistics(imphys_freeze) - dn_cl_hom
  endif
 end subroutine homfreez3


! ****************************************
! Heterogeneous freezing of cloud droplets
! ****************************************
subroutine hetfreez3
  use modglobal, only : rlv
  implicit none

  real :: J_het   ! freezing rate
  real :: dq_cl_het = 0.   !< heterogeneou freezing of cloud water

  J_het = A_het *exp( B_het*(T_3-tmp0) -1)

  dn_cl_het = -c_mmt_1cl * n_cl * x_cl * J_het
  dq_cl_het = -c_mmt_2cl * q_cl * x_cl * J_het

  ! basic correction
  dn_cl_het = max(dn_cl_het,min(0.0,-n_clm/delt-n_clp))
  dq_cl_het = max(dq_cl_het,min(0.0,-q_clm/delt-q_clp))

  ! changes in cloud water
  n_clp = n_clp + dn_cl_het
  q_clp = q_clp + dq_cl_het

  ! increase in cloud ice
  n_cip = n_cip - dn_cl_het
  q_cip = q_cip - dq_cl_het

  ! changes in the total amount of water
  qtpmcr = qtpmcr + dq_cl_het

  ! change in th_l due to freezing
  thlpmcr = thlpmcr-((rlv+rlme)/cp_exnf_k)*(dq_cl_het)
  if (l_tendencies) then
    tend(idq_cl_het) = dq_cl_het
    tend(idn_cl_het) = dn_cl_het
  endif
  if (l_statistics) then
    statistics(imphys_freeze) = statistics(imphys_freeze) - dn_cl_het
  endif
end subroutine hetfreez3


! ****************************************
! Depositional growth of cloud ice particles
!  later tunr into:
!
! Depositional growth of various ice particles
! Inner subroutine - takes input from a wrapper
!  following S&B
! ****************************************
subroutine deposit_ice3
  use modglobal, only : rv,rd,pi
  implicit none

  real :: esi

  real :: F       & ! ventilation factor
         ,q_avail & ! available water for deposition
         ,Si      & ! super or undersaturation with respect to ice
         ,G       & ! G_iv function
         ,Dvic_ci & ! size of ice parti
         ,viic_ci & ! v for ice cloud particles
         ,nrex_ci   ! reynolds number

  q_avail = qt0 - q_cl - qvsi ! NOTE: q_avail < 0 sublimation instead of deposition
  Si = q_avail/qvsi

  ! calculating G_iv
  esi = qvsi*presf_k/(rd/rv+(1.0-rd/rv)*qvsi)
  G = (rv * tmp0) / (Dv*esi) + rlvi/(Kt*tmp0)*(rlvi/(rv*tmp0) -1.)
  G = 1./G

  ! diameter
  Dvic_ci = a_ci*x_ci**b_ci

  ! terminal velocity
  viic_ci = al_ci*((rho0s/rhof_k)**0.5)*x_ci**be_ci

  ! N_re Reynolds number
  nrex_ci = Dvic_ci*viic_ci/nu_a

  ! calculating from prepared ventilation coefficients
  F = aven_1i+bven_1i*Sc_num**(1.0/3.0)*nrex_ci**0.5

  ! depositional growth
  ! k_depos = 4*pi/c_ci  for spherical particles
  ! k_depos = 4*pi/c_ci  for hexagonal plates - included
  dq_ci_dep = (4*pi/c_ci)*n_ci*G*Dvic_ci*F*Si

  if (dq_ci_dep < 0.) then
    ! limiting not to sublimate more than available
    dq_ci_dep = max(dq_ci_dep,-q_cim/delt-q_cip)
  else
    ! limiting not to deposit more than available
    dq_ci_dep = min(dq_ci_dep,q_avail/delt)
  endif

  ! adds to outputs
  q_cip = q_cip + dq_ci_dep

  qtpmcr = qtpmcr - dq_ci_dep
  thlpmcr = thlpmcr + (rlvi/cp_exnf_k)*dq_ci_dep
end subroutine deposit_ice3


! ****************************************
! Depositional growth of snow particles
!  later turn into:
!
! Depositional growth of various ice particles
! Inner subroutine - takes input from a wrapper
!  following S&B
! ****************************************
subroutine deposit_snow3
  use modglobal, only : rv,rd,pi
  implicit none

  real :: esi

  real :: F_hs    & ! ventilation factor
         ,q_avail & ! available water for deposition
         ,Si      & ! super or undersaturation with respect to ice
         ,G       & ! G_iv function
         ,Dvic_hs & ! size of ice parti
         ,viic_hs & ! v for ice cloud particles
         ,nrex_hs   ! reynolds number

  q_avail = qt0 - q_cl - qvsi ! NOTE: q_avail < 0 sublimation instead of deposition
  Si = q_avail/qvsi

  ! calculating G_iv
  esi = qvsi*presf_k/(rd/rv+(1.0-rd/rv)*qvsi)
  G = (rv * tmp0) / (Dv*esi) + rlvi/(Kt*tmp0)*(rlvi/(rv*tmp0) -1.)
  G = 1./G

  ! diameter
  Dvic_hs = a_hs*x_hs**b_hs

  ! terminal velocity
  viic_hs = al_hs*((rho0s/rhof_k)**0.5)*x_hs**be_hs

  ! N_re Reynolds number
  nrex_hs = Dvic_hs*viic_hs/nu_a

  ! calculating from prepared ventilation coefficients
  F_hs = aven_1s+bven_1s*Sc_num**(1.0/3.0)*nrex_hs**0.5

  ! depositional growth
  ! k_depos = 4*pi/c  for spherical particles
  ! k_depos = 4*pi/c  for hexagonal plates - included
  dq_hs_dep = (4*pi/c_hs)*n_hs*G*Dvic_hs*F_hs*Si

  ! limiting not to sublimate more than available
  dq_hs_dep = max(dq_hs_dep,min(0.0,-q_hsm/delt-q_hsp))

  ! limiting not to deposit more than available
  dq_hs_dep = min(dq_hs_dep,max(0.0,q_avail/delt))

  ! adds to outputs
  q_hsp = q_hsp + dq_hs_dep

  qtpmcr = qtpmcr - dq_hs_dep
  thlpmcr = thlpmcr + (rlvi/cp_exnf_k)*dq_hs_dep
end subroutine deposit_snow3


! ****************************************
! Depositional growth of graupel particles
! as well as sublimation
!  later turn into:
!
! Depositional growth of various ice particles
! Inner subroutine - takes input from a wrapper
!  following S&B
! ****************************************
subroutine deposit_graupel3
  use modglobal, only : rv,rd,pi
  implicit none

  real :: esi

  real :: F_hg    & ! ventilation factor
         ,q_avail & ! available water for deposition
         ,Si      & ! super or undersaturation with respect to ice
         ,G       & ! G_iv function
         ,Dvic_hg & ! size of particles
         ,viic_hg & ! v for particles
         ,nrex_hg   ! reynolds number

  q_avail = qt0 - q_cl - qvsi ! NOTE: q_avail < 0 sublimation instead of deposition
  Si = q_avail/qvsi

  ! calculating G_iv
  esi = qvsi*presf_k/(rd/rv+(1.0-rd/rv)*qvsi)
  G = (rv * tmp0) / (Dv*esi) + rlvi/(Kt*tmp0)*(rlvi/(rv*tmp0) -1.)
  G = 1./G

  ! diameter
  Dvic_hg = a_hg*x_hg**b_hg

  ! terminal velocity
  viic_hg = al_hg*((rho0s/rhof_k)**0.5)*x_hg**be_hg

  ! N_re Reynolds number
  nrex_hg = Dvic_hg*viic_hg/nu_a

  ! calculating from prepared ventilation coefficients
  F_hg = aven_1g+bven_1g*Sc_num**(1.0/3.0)*nrex_hg**0.5

  ! depositional growth
  ! k_depos = 4*pi/c  for spherical particles
  ! k_depos = 4*pi/c  for hexagonal plates - included
  dq_hg_dep = (4*pi/c_hg)*n_hg*G*Dvic_hg*F_hg*Si

  ! limiting not to sublimate more than available
  dq_hg_dep = max(dq_hg_dep,min(0.0,-q_hgm/delt-q_hgp))

  ! limiting not to deposit more than available
  dq_hg_dep = min(dq_hg_dep,max(0.0,q_avail/delt))

  ! adds to outputs
  q_hgp = q_hgp + dq_hg_dep

  qtpmcr = qtpmcr - dq_hg_dep
  thlpmcr = thlpmcr + (rlvi/cp_exnf_k)*dq_hg_dep
end subroutine deposit_graupel3


!! ****************************************************************
!!  Limiting condensation and deposition
!!  - to prevent negative values
!!
!!  - call should be located:
!!      - after depositions
!!      - after nucleation correction
!!      - before heterogeneous freezing
!!
!!  ************************************************************
subroutine cor_deposit3
  implicit none

  real    :: tocon,precon,cond_cf
  real    :: cor_dqci_dep,cor_dqhs_dep,cor_dqhg_dep

  ! available water vapour for deposition
  tocon = (qt0-q_clm-qvsi)/delt

  ! consumption of water vapour calculated by nucleation and deposition processes
  precon = dq_ci_dep+dq_hs_dep+dq_hg_dep

  ! run only in oversaturted conditions
  if ((precon.gt.0.0).and.(tocon-precon).lt.0.0) then
    ! preparing additive correctors:
    cond_cf = tocon/precon - 1.0
    cond_cf = max(min(0.0, cond_cf),-1.0)

    ! - corrector for deposition - only if positive deposition
    cor_dqci_dep = cond_cf*max(0.0, dq_ci_dep)
    cor_dqhs_dep = cond_cf*max(0.0, dq_hs_dep)
    cor_dqhg_dep = cond_cf*max(0.0, dq_hg_dep)

    ! and updating values:
    q_cip = q_cip + cor_dqci_dep
    q_hsp = q_hsp + cor_dqhs_dep
    q_hgp = q_hgp + cor_dqhg_dep

    ! - correction for total water
    qtpmcr = qtpmcr -cor_dqhs_dep &
                    -cor_dqhg_dep &
                    -cor_dqci_dep

    ! - and correcting for heat
    thlpmcr = thlpmcr + (rlvi/cp_exnf_k)*cor_dqci_dep  &
                      + (rlvi/cp_exnf_k)*cor_dqhs_dep  &
                      + (rlvi/cp_exnf_k)*cor_dqhg_dep

    ! corrector for process values
    dq_ci_dep = dq_ci_dep + cor_dqci_dep
    dq_hs_dep = dq_hs_dep + cor_dqhs_dep
    dq_hg_dep = dq_hg_dep + cor_dqhg_dep
  endif
  if (l_tendencies) then
    tend(idq_ci_dep) = dq_ci_dep
    tend(idq_hs_dep) = dq_hs_dep
    tend(idq_hg_dep) = dq_hg_dep
  endif
  if (l_statistics) then
    statistics(imphys_dep) = statistics(imphys_dep) + max(0.0,dq_ci_dep)
    statistics(imphys_dep) = statistics(imphys_dep) + max(0.0,dq_hs_dep)
    statistics(imphys_dep) = statistics(imphys_dep) + max(0.0,dq_hg_dep)

    statistics(imphys_sub) = statistics(imphys_sub) + min(0.0,dq_ci_dep)
    statistics(imphys_sub) = statistics(imphys_sub) + min(0.0,dq_hs_dep)
    statistics(imphys_sub) = statistics(imphys_sub) + min(0.0,dq_hg_dep)
  endif
end subroutine cor_deposit3


! ****************************************
! ice selfcollection and aggregation to snow
!  - all collection of ice by ice treated as snow
!
! follows Seifert, 2002
! ****************************************
subroutine ice_aggr3
  use modglobal, only : pi
  implicit none

  real    :: dlt_0aa_i &
            ,dlt_1aa_i &
            ,th_0aa_i  &
            ,th_1aa_i

  real :: dif_D_10, x_minagg_ii, rem_cf_i
  real :: E_ab, E_stick, x_crit_ii

  real :: dq_ci_col_iis = 0. !< self-collection of cloud ice
  real :: dn_ci_col_iis = 0. !< self-collection of cloud ice

  x_crit_ii   = (D_crit_ii/a_ci)**(1.0/b_ci)  ! TODO: this is a constant, precalculate

  if((x_ci.gt.x_crit_ii).and.(q_ci.gt.q_crit_ii)) then
    ! prepare coefficient for remaining water number
    rem_cf_i = (1.0-rem_n_ci_min)/delt

    ! calculate constants
    dlt_0aa_i   = 2*dlt_i0 + dlt_i0i
    dlt_1aa_i   = dlt_i0 + dlt_i1i + dlt_i1
    th_0aa_i    = 2*th_i0 - th_i0i   ! from Seifert, 2002
    th_1aa_i    = th_i0 - th_i1i + th_i1

    ! and the minimal conversion size
    x_minagg_ii = (D_i_b/a_ci)**(1.0/b_ci)
    ! calculating sticking efficiency
    if (l_sb_stickyice) then
      E_stick = c_E_o_s*exp(B_stick *(tmp0+stick_off))
      E_stick = min(c_E_o_s,E_stick)
    else ! l_sb_stickyice
      E_stick = exp(B_stick_ii*(tmp0+stick_off)+C_stick_ii)
      E_stick = min(E_ii_maxst,E_stick)
    endif ! l_sb_stickyice

    ! collision efficiency
    if (l_sb_lim_aggr) then
      ! checking whether sufficient size
      if(D_ci.gt.D_i_a ) then
        if(D_ci.gt.D_i_b ) then
          E_ab = E_ee_m*E_stick
        else
          dif_D_10  = D_i_b-D_i_a  ! difference in ice cloud droplet intervals
          E_ab = (E_ee_m*E_stick /dif_D_10)* (D_ci - D_i_a)
        endif
      else
        E_ab = 0.
      endif
    else ! l_sb_lim_aggr
      E_ab = E_ee_m*E_stick
    endif

    dn_ci_col_iis = -rhof_k*(pi/4)*E_ab                      &
                    *n_ci**2 *dlt_0aa_i*D_ci**2              &
                    *(th_0aa_i *v_ci**2+2*sigma_ci**2)**0.5

    dq_ci_col_iis = -rhof_k*(pi/4)*E_ab                      &
                    *dlt_1aa_i*n_ci*q_ci*D_ci**2             &
                    *(th_1aa_i*v_ci**2+2*sigma_ci**2)**0.5


    ! limiting dq_ci  --   limit x_cv_ii as per ICON 2017
    dq_ci_col_iis = max(dq_ci_col_iis, (-q_cim/delt-q_cip))
    dn_ci_col_iis = max(dn_ci_col_iis, (-rem_cf_i*n_cim-n_cip))
    dn_ci_col_iis = max(dn_ci_col_iis, dq_ci_col_iis/x_minagg_ii)
  endif

  n_cip    = n_cip + dn_ci_col_iis
  q_cip    = q_cip + dq_ci_col_iis
  n_hsp    = n_hsp - dn_ci_col_iis
  q_hsp    = q_hsp - dq_ci_col_iis

  if (l_sb_dbg) then
    if((q_cim+delt*dq_ci_col_iis .lt. 0.0)) then
      write(6,*) 'WARNING: ice_aggr3 too high'
      write(6,*) ' removing more ice than available'
      ! count((q_cim +delt*dq_ci_col_iis).lt. 0.0)
      write(6,*) ' removing too much ice'
      ! count(( q_ci+delt*dq_ci_col_iis).lt. 0.0 )
      write(6,*) ' getting negative q_t'
      ! count(( qt0+delt*q_cip).lt. 0.0 )
    endif

    if((n_cim+delt*dn_ci_col_iis .lt. 0.0)) then
      write(6,*) 'WARNING: ice_aggr3 too high'
      write(6,*) ' removing more ice particles then available'
      ! count((n_cim +delt*dn_ci_col_iis).lt. 0.0)
      write(6,*) ' removing too much ice particles'
      ! count((n_ci+delt*dn_ci_col_iis).lt. 0.0)
    endif
    if (l_tendencies) then
      tend(idq_ci_col_iis) = dq_ci_col_iis
      tend(idn_ci_col_iis) = dn_ci_col_iis
    endif
  endif
end subroutine ice_aggr3


! *************************************************
! snow selfcollection
!
! *************************************************
subroutine snow_self3
  use modglobal, only : pi
  implicit none

  real :: dlt_0aa, th_0aa
  real :: rem_cf_s, E_ab_s
  real :: E_stick
  real :: dn_hs_col_sss = 0.  !< self-collection of snow

  ! adjusting coefficient
  ! prepare coefficient for remaining water number
  rem_cf_s = (1.0-rem_n_hs_min)/delt
  dlt_0aa = 2*dlt_s0 + dlt_s0s
  th_0aa = 2*th_s0 - th_s0s   ! from Seifert, 2002

  ! calculating sticking efficiency
  if (l_sb_stickyice) then
    E_stick = c_E_o_s*exp(B_stick *(tmp0+stick_off))
    E_stick = min(c_E_o_s,E_stick)
  else ! l_sb_stickyice
    E_stick = exp(B_stick_ii*(tmp0+stick_off)+C_stick_ii)
    E_stick = min(E_ss_maxst,E_stick)
  endif ! l_sb_stickyice
  E_ab_s = E_ee_m*E_stick

  dn_hs_col_sss = -rhof_k*(pi/4)*E_ab_s           &
                  *n_hs**2 *dlt_0aa*D_hs**2      &
                  *(th_0aa*v_hs**2+2*sigma_hs**2)**0.5

  dn_hs_col_sss = max(min(0.0,dn_hs_col_sss),-rem_cf_s*n_hsm-n_hsp)

  n_hsp    = n_hsp + dn_hs_col_sss

  if (l_sb_dbg) then
    if(( n_hsm+delt*dn_hs_col_sss .lt. 0.0 )) then
      write(6,*) 'WARNING: snow self-collection too high'
      write(6,*) ' decreasing number of snowflakes below 0'
      ! count((n_hsm+delt*dn_hs_col_sss).lt. 0.0 )
      write(6,*) ' decreasing number of snowflakes too much'
      ! count((n_hs+delt*dn_hs_col_sss).lt. 0.0 )
    endif
  endif
  if (l_tendencies) then
    tend(idn_hs_col_sss) =  dn_hs_col_sss
  endif
end subroutine snow_self3


! ****************************************
! snow collecting cloud ice
!
! s + i -> s
! ****************************************
subroutine coll_sis3
  use modglobal, only : pi
  implicit none

  real :: E_ab, E_stick, rem_cf_i

  real :: dq_hsci_col   = 0.  !< collection s+i - trend in q_hs
  real :: dn_ci_col_hs  = 0.  !< collection s+i - trend in n_ci

  ! adjusting coefficient
  ! prepare coefficient for remaining water number
  rem_cf_i = (1.0-rem_n_ci_min)/delt

  ! calculating sticking efficiency
  E_stick = c_E_o_s*exp(B_stick *(tmp0+stick_off))
  E_stick =min(c_E_o_s,E_stick)
  E_ab = E_ee_m*E_stick

  dq_hsci_col = (rhof_k*pi/4)*E_ab*n_hs                       &
                *q_ci*(dlt_s0*D_hs**2                          &
                +dlt_s1i*D_hs*D_ci+dlt_i1*D_ci**2)             &
                *(th_s0*v_hs**2-th_s1i*v_ci*v_hs               &
                +th_i1*v_ci**2+sigma_hs**2+sigma_ci**2)**0.5

  dn_ci_col_hs = -(rhof_k*pi/4)*E_ab*n_hs                     &
                 *n_ci*(dlt_s0*D_hs**2                         &
                 +dlt_s0i*D_hs*D_ci+dlt_i0*D_ci**2)            &
                 *(th_s0*v_hs**2-th_s0i*v_ci*v_hs              &
                 +th_i0*v_ci**2+sigma_hs**2+sigma_ci**2)**0.5

  dq_hsci_col = min(dq_hsci_col,max(0.0,q_cim/delt+q_cip))   ! following ICON, 2017
  dn_ci_col_hs = max(dn_ci_col_hs,min(0.0,-rem_cf_i*n_cim-n_cip))

  n_cip    = n_cip + dn_ci_col_hs
  q_cip    = q_cip + dq_hsci_col
  q_hsp    = q_hsp - dq_hsci_col

  if (l_sb_dbg) then
    if(q_cim-delt*dq_hsci_col.lt. 0.0) then
      write(6,*) 'WARNING: coll_sis3 too high removing more ice than available'
      write(6,*) ' removing more ice than available'
      ! count((q_cim-delt*dq_hsci_col).lt. 0.0)
      write(6,*) ' removing too much ice'
      ! count((q_ci-delt*dq_hsci_col).lt. 0.0)
      write(6,*) ' getting negative q_t'
      ! count((qt0-delt*q_cip).lt. 0.0 )
    endif

    if(n_cim+delt*dn_ci_col_hs.lt. 0.0) then
      write(6,*) 'WARNING: coll_sis3 too high'
      write(6,*) ' removing more ice particles then available in gridpoints'
      ! count(n_cim+delt*dn_ci_col_hs.lt. 0.0)
      write(6,*) ' removing too many ice particles in gridpoints '
      ! count((n_ci+delt*dn_ci_col_hs).lt. 0.0)
    endif
  endif
  if (l_tendencies) then
    tend(idq_hsci_col ) = dq_hsci_col
    tend(idn_ci_col_hs) = dn_ci_col_hs
  endif
end subroutine coll_sis3


! ****************************************
! graupel collecting snow
!
! g+s -> g
! ****************************************
subroutine coll_gsg3
  use modglobal, only : pi
  implicit none

  real :: E_ab, E_stick, rem_cf_s
  real :: dq_hghs_col  = 0.   !< collection g+s - trend in q_hg
  real :: dn_hs_col_hg = 0.   !< collection g+s - trend in n_hs

  ! adjusting coefficient
  ! prepare coefficient for remaining water number
  rem_cf_s = (1.0-rem_n_hs_min)/delt

  ! calculating sticking efficiency
  E_stick = c_E_o_s*exp(B_stick *(tmp0+stick_off))
  E_stick =min(c_E_o_s,E_stick)
  E_ab = E_ee_m*E_stick

  dq_hghs_col  = (rhof_k*pi/4)*E_ab*n_hg              &
       *q_hs*(dlt_g0*D_hg**2                           &
         +dlt_g1s*D_hg*D_hs+dlt_s1*D_hs**2)            &
       *( th_g0*v_hg**2-th_g1s*v_hs*v_hg               &
         +th_s1*v_hs**2+sigma_hg**2+sigma_hs**2)**0.5

  dn_hs_col_hg = -(rhof_k*pi/4)*E_ab*n_hg             &
       *n_hs*(dlt_g0*D_hg**2                           &
         +dlt_g0s*D_hg*D_hs+dlt_s0*D_hs**2)            &
       *( th_g0*v_hg**2-th_g0s*v_hs*v_hg               &
         +th_s0*v_hs**2+sigma_hg**2+sigma_hs**2)**0.5

  dq_hghs_col = min(dq_hghs_col,&
                    max(0.0,    &
                    q_hsm/delt+q_hsp))   ! following ICON, 2017
  dn_hs_col_hg = max(dn_hs_col_hg,&
                     min(0.0,     &
                     -rem_cf_s*n_hsm-n_hsp))

  n_hsp    = n_hsp + dn_hs_col_hg
  q_hsp    = q_hsp - dq_hghs_col
  q_hgp    = q_hgp + dq_hghs_col

  if (l_sb_dbg) then
    if((q_hsm-delt*dq_hghs_col).lt. 0.0) then
      write(6,*) 'WARNING: coll_gsg3 too high removing more ice than available'
      write(6,*) ' removing more ice than available'
      ! count(q_hsm-delt*dq_hghs_col.lt. 0.0)
      write(6,*) ' removing too much ice'
      ! count(q_hs-delt*dq_hghs_col.lt. 0.0)
      write(6,*) ' getting negative q_t'
      ! count(qt0-delt*q_hsp.lt. 0.0)
    endif

    if(n_hsm+delt*dn_hs_col_hg.lt. 0.0) then
      write(6,*) 'WARNING: coll_gsg3 too high'
      write(6,*) ' removing more ice particles then available'
      ! count(n_hsm+delt*dn_hs_col_hg.lt. 0.0)
      write(6,*) ' removing too many ice particles'
      ! count(n_hs+delt*dn_hs_col_hg.lt. 0.0)
    endif
  endif
  if (l_tendencies) then
    tend(idq_hghs_col ) = dq_hghs_col
    tend(idn_hs_col_hg) = dn_hs_col_hg
  endif
end subroutine coll_gsg3


!****************************************
! Collection of clpoud droplets by ice
!
! - based on Seifert & Beheng (2004)
! - resulting process is:
!    - ice riming by cloud ice
!    - enhanced melting
!
! conditions: (D_cl.gt.D_c_a).and.(D_ci.gt.D_i0)
!    q_cl_mask = .true.
!    q_ci_mask = .true.
!
!****************************************
subroutine coll_ici3
  use modglobal, only : pi
  implicit none

  real :: E_ab
  real :: dif_D_10
  real :: k_enhm
  real :: rem_cf

  real :: dn_col_ici   = 0.
  real :: dq_col_ici   = 0.
  real :: dn_ci_eme_ic = 0.   !< number tendency enhanced melting of cloud ice by cloud water
  real :: dq_ci_eme_ic = 0.   !< mass tendency enhanced melting of cloud ice by cloud water

  if( D_cl.gt.D_c_b ) then
    E_ab = E_i_m
  else
    ! denominator in calculationg collision efficiency
    dif_D_10  = D_c_b-D_c_a

    E_ab = (E_i_m /dif_D_10)* (D_cl - D_c_a)
  endif

  dq_col_ici = (rhof_k*pi/4)*E_ab*n_ci      &
       *q_cl*(dlt_i0*D_ci**2                &
         +dlt_i1c*D_ci*D_cl+dlt_c1*D_cl**2) &
       *( th_i0*v_ci**2-th_i1c*v_cl*v_ci    &
         +th_c1*v_cl**2+sigma_ci**2+sigma_cl**2)**0.5

  dn_col_ici = -(rhof_k*pi/4)*E_ab*n_ci      &
       *n_cl*(dlt_i0*D_ci**2                 &
         +dlt_i0c*D_ci*D_cl+dlt_c0*D_cl**2)  &
       *( th_i0*v_ci**2-th_i0c*v_cl*v_ci     &
         +th_c0*v_cl**2+sigma_ci**2+sigma_cl**2)**0.5

  if(tmp0.lt.T_3) then
    ! riming
    ! adjusting coefficient
    ! prepare coefficient for remaining water number
    rem_cf = (1.0-rem_n_cl_min)/delt

    dq_ci_rime = min(dq_col_ici,max(0.0,q_clm/delt+q_clp))   ! following ICON, 2017
    dn_cl_rime_ci = max(dn_col_ici,min(0.0,-rem_cf*n_clm-n_clp))
  else
    ! scheding and enhanced melting
    k_enhm  = c_water/rlme

    ! calculating the melting
    dq_ci_eme_ic = -k_enhm*(tmp0-T_3)*min(dq_col_ici,q_cl/delt)
    dq_ci_eme_ic = max(dq_ci_eme_ic,min(0.0,-q_cim/delt-q_cip))

    ! calculating number of melted particles
    ! - expected to be proportional to melted mass
    dn_ci_eme_ic = dq_ci_eme_ic*n_ci/q_ci

    ! - but not more than number of interacting particles
    dn_ci_eme_ic = max(dn_ci_eme_ic, max(dn_col_ici,-n_cl/delt))

    ! - and not more than total number of particles
    dn_ci_eme_ic = max(dn_ci_eme_ic,-n_cim/delt-n_cip)
  endif

  n_clp    = n_clp + dn_cl_rime_ci
  q_clp    = q_clp - dq_ci_rime - dq_ci_eme_ic
  n_cip    = n_cip + dn_ci_eme_ic
  q_cip    = q_cip + dq_ci_rime + dq_ci_eme_ic

  qtpmcr = qtpmcr - dq_ci_rime - dq_ci_eme_ic
  thlpmcr = thlpmcr+(rlvi/cp_exnf_k)*(dq_ci_rime + dq_ci_eme_ic)

  if (l_sb_dbg) then
    if(q_clm-delt*dq_col_ici.lt. 0.0) then
      write(6,*) 'WARNING: coll_ici3 too high'
      write(6,*) ' removing more cloud water then available'
      ! count(q_clm-delt*dq_ci_rime.lt. 0.0)
      write(6,*) ' removing too much water'
      ! count(q_cl-delt*dq_ci_rime.lt. 0.0)
      write(6,*) ' getting negative q_t'
      ! count(qt0+delt*q_cip.lt. 0.0)
    endif
    if(n_clm+delt*dn_col_ici.lt. 0.0) then
      write(6,*) 'WARNING: coll_ici3 too high'
      write(6,*) ' removing more droplets then available'
      ! count(n_clm+delt*dn_cl_rime_ci.lt. 0.0)
      write(6,*) ' removing too many droplets'
      ! count(n_cl+delt*dn_cl_rime_ci.lt. 0.0)
    endif
  endif
  if (l_tendencies) then
    tend(idn_ci_eme_ic) = dn_ci_eme_ic
    tend(idq_ci_eme_ic) = dq_ci_eme_ic
    tend(idq_ci_rime)   = dq_ci_rime
    tend(idn_cl_rime_ci)= dn_cl_rime_ci
  endif
  if (l_statistics) then
    statistics(imphys_freeze) = statistics(imphys_freeze) + dq_ci_rime
    statistics(imphys_melt) = statistics(imphys_melt) + dq_ci_eme_ic
  endif
end subroutine coll_ici3


! ****************************************
!  riming of snow
!
! conditions: (D_cl.gt.D_c_a).and.(D_hs.gt.D_i0)
! ****************************************
subroutine coll_scs3
  use modglobal, only : pi
  implicit none

  real :: E_ab
  real :: dif_D_10
  real :: k_enhm
  real :: dn_col_b, dq_col_a
  real :: dn_hs_eme_sc = 0.   !< number tendency enhanced melting of snow by cloud water
  real :: dq_hs_eme_sc = 0.   !< mass tendency enhanced melting of snow by cloud water

  if( D_cl.gt.D_c_b ) then
    E_ab = E_s_m
  else
    ! denominator in calculationg collision efficiency
    dif_D_10  = D_c_b-D_c_a
    E_ab = (E_s_m/dif_D_10)* (D_cl - D_c_a)
  endif

  dq_col_a = (rhof_k*pi/4)*E_ab*n_hs                         &
       *q_cl*(dlt_s0*D_hs**2                                  &
         +dlt_s1c*D_hs*D_cl+dlt_c1*D_cl**2)                   &
       *( th_s0*v_hs**2-th_s1c*v_cl*v_hs                      &
         +th_c1*v_cl**2+sigma_hs**2+sigma_cl**2)**0.5

  dn_col_b = -(rhof_k*pi/4)*E_ab*n_hs                        &
       *n_cl*(dlt_s0*D_hs**2                                  &
         +dlt_s0c*D_hs*D_cl+dlt_c0*D_cl**2)                   &
       *( th_s0*v_hs**2-th_s0c*v_cl*v_hs                      &
         +th_c0*v_cl**2+sigma_hs**2+sigma_cl**2)**0.5

  ! basic limiting
  ! limited by amount of water in droplets
  dq_col_a = min(dq_col_a,q_cl/delt)
  dq_col_a = min(dq_col_a,max(0.0,q_clm/delt+q_clp))
  dn_col_b = max(dn_col_b,min(0.0,-n_clm/delt-n_clp))
  ! limited by number of droplets
  ! then based on temperature

  if(tmp0.lt.T_3) then
    ! riming only
    dq_hs_rime = dq_col_a ! following ICON, 2017

    ! = min(dq_hs_rime,q_clm/delt)
    dn_cl_rime_hs = dn_col_b ! min(dn_cl_rime_hs,-n_clm/delt)

    ! record the change
    ! o dq_hs_rime = dq_col_a
    ! change in the amount of cloud ice
    q_hsp = q_hsp + dq_hs_rime

    ! change in the amount of cloud water
    n_clp = n_clp + dn_cl_rime_hs
    q_clp = q_clp - dq_hs_rime

    ! change in q_t
    qtpmcr = qtpmcr - dq_hs_rime

    ! change in th_l - freezing and removal
    thlpmcr = thlpmcr+ (rlvi/cp_exnf_k)*dq_hs_rime
  else
    ! not riming, but enhanced melting and scheding

    ! enhanced melting
    k_enhm  = c_water/rlme

    ! calculating the melting
    dq_hs_eme_sc = -k_enhm*(tmp0-T_3)*dq_col_a
    dq_hs_eme_sc = max(dq_hs_eme_sc, &
                       min(0.0, &
                           -q_hsm/delt-q_hsp))

    ! calculating number of melted particles
    ! - expected to be proportional to melted mass
    dn_hs_eme_sc = dq_hs_eme_sc*n_hs/q_hs

    ! - not more than number of interacting particles
    ! dn_hs_eme_sc = max(dn_hs_eme_sc, max(dn_col_b,-n_cl/delt))
    ! - and not more than total number of particles
    dn_hs_eme_sc = max(dn_hs_eme_sc, &
                       -min((n_clm/delt+n_clp), &
                            (n_hsm/delt+n_hsp)))

    ! updating tendencies
    ! updating rain
    ! melted snow is turning into rain  RH84
    q_hrp = q_hrp -  dq_hs_eme_sc + dq_col_a ! both mass of snow and droplets
    n_hrp = n_hrp -  dn_hs_eme_sc

    ! snow
    q_hsp = q_hsp + dq_hs_eme_sc
    n_hsp = n_hsp + dn_hs_eme_sc

    ! updating cloud water
    q_clp = q_clp - dq_col_a
    n_clp = n_clp + dn_col_b

    ! updating thermodynamic
    ! qtp : increased by melted water
    qtpmcr = qtpmcr -dq_col_a ! -dq_hs_eme_sc

    ! thlp : melting and adding liquid water
    ! thlpmcr = thlpmcr+(rlvi/cp_exnf_k)*dq_hs_eme_sc
    thlpmcr = thlpmcr+            &
        (rlme/cp_exnf_k)*dq_hs_eme_sc        &
       +(rlvi/cp_exnf_k)*dq_col_a
  endif

  if (l_sb_dbg) then
    if((q_clm-delt*dq_col_a).lt. 0.0) then
      write(6,*) 'WARNING: coll_scs too high'
      write(6,*) ' removing more cloud water than available'
      ! count((q_clm-delt*dq_hs_rime).lt. 0.0)
      write(6,*) ' removing too much cloud water'
      ! count(q_cl-delt*dq_hs_rime.lt. 0.0)
      write(6,*) ' getting negative q_t'
      ! count(qt0+delt*q_clp).lt. 0.0)
    endif

    if(n_clm+delt*dn_col_b.lt. 0.0) then
      write(6,*) 'WARNING: coll_scs too high'
      write(6,*) ' removing more droplets then available'
      ! count(n_clm+delt*dn_cl_rime_hs).lt. 0.0)
      write(6,*) ' removing too many droplets'
      ! count(n_cl+delt*dn_cl_rime_hs.lt. 0.0)
    endif
  endif
  if (l_tendencies) then
    tend(idn_hs_eme_sc ) = dn_hs_eme_sc
    tend(idq_hs_eme_sc ) = dq_hs_eme_sc
    tend(idq_hs_rime   ) = dq_hs_rime
    tend(idn_cl_rime_hs) = dn_cl_rime_hs
  endif
  if (l_statistics) then
    statistics(imphys_freeze) = statistics(imphys_freeze) + dq_hs_rime
    statistics(imphys_melt) = statistics(imphys_melt) + dq_hs_eme_sc
  endif
end subroutine coll_scs3


! ****************************************
!  riming of graupel by cloud droplets
!
!
! conditions: (D_cl.gt.D_c_a).and.(D_hg.gt.D_i0)
!             q_hg_mask = .true.
!             q_cl_mask = .true.
!
! ****************************************
subroutine coll_gcg3
  use modglobal, only : pi
  implicit none

  real :: E_ab
  real :: dif_D_10
  real :: k_enhm
  real :: rem_cf
  real :: dn_col_b, dq_col_a
  real :: dn_hg_eme_gc = 0.  !< number tendency enhanced melting of graupel
  real :: dq_hg_eme_gc = 0.  !< mass tendency enhanced melting of graupel

  ! collision efficiency
  if( D_cl.gt.D_c_b ) then
    E_ab = E_g_m
  else
    dif_D_10 = D_c_b-D_c_a
    E_ab = (E_g_m/dif_D_10)* (D_cl- D_c_a)
  endif

  dq_col_a = (rhof_k*pi/4)*E_ab*n_hg    &
       *q_cl*(dlt_g0*D_hg**2                     &
         +dlt_g1c*D_hg*D_cl+dlt_c1*D_cl**2) &
       *( th_g0*v_hg**2-th_g1c*v_cl*v_hg    &
         +th_c1*v_cl**2+sigma_hg**2+sigma_cl**2)**0.5

  dn_col_b = -(rhof_k*pi/4)*E_ab*n_hg   &
       *n_cl*(dlt_g0*D_hg**2                     &
         +dlt_g0c*D_hg*D_cl+dlt_c0*D_cl**2) &
       *( th_g0*v_hg**2-th_g0c*v_cl*v_hg    &
         +th_c0*v_cl**2+sigma_hg**2+sigma_cl**2)**0.5

  ! initial correction based on amount of cloud water
  rem_cf = (1.0-rem_n_cl_min)/delt
  dq_col_a = min(dq_col_a,max(0.0,q_clm/delt+q_clp)) ! following ICON, 2017
  dn_col_b = max(dn_col_b,min(0.0,-rem_cf*n_clm-n_clp))

  ! then based on temeprature
  if(tmp0.lt.T_3) then
    dq_hg_rime = dq_col_a    ! = min(dq_hg_rime,q_clm/delt)
    dn_cl_rime_hg = dn_col_b ! = min(dn_cl_rime_hg,-n_clm/delt)

    ! change in the amount of graupel
    q_hgp = q_hgp + dq_hg_rime

    ! change in the amount of rain
    n_clp = n_clp + dn_cl_rime_hg
    q_clp = q_clp - dq_hg_rime

    ! no change in q_t
    qtpmcr = qtpmcr - dq_hg_rime

    ! change in th_l - heat release from freezing and removal
    thlpmcr = thlpmcr+ (rlvi/cp_exnf_k)*dq_hg_rime
  else
    ! not riming,but enhanced melting and scheding

    ! enhanced melting
    k_enhm  = c_water/rlme

    ! calculating the melting
    dq_hg_eme_gc = -k_enhm*(tmp0-T_3)*min(dq_col_a,q_cl/delt)
    dq_hg_eme_gc = max(dq_hg_eme_gc,min(0.0,-q_hgm/delt-q_hgp))

    ! calculating number of melted particles
    ! - expected to be proportional to melted mass
    dn_hg_eme_gc = dq_hg_eme_gc*n_hg/q_hg

    ! - but not more than number of interacting particles
    ! dn_hg_eme_gc = max(dn_hg_eme_gc, max(dn_col_b,-n_cl/delt))
    ! - and not more than total number of particles
    dn_hg_eme_gc = max(dn_hg_eme_gc, -min((n_clm/delt+n_clp),(n_hgm/delt+n_hgp)))

    ! updating tendencies
    ! updating rain
    ! based on RH84
    q_hrp = q_hrp - dq_hg_eme_gc + dq_col_a
    n_hrp = n_hrp - dn_col_b

    ! cloud ice
    q_hgp = q_hgp + dq_hg_eme_gc
    n_hgp = n_hgp + dn_hg_eme_gc

    ! updating cloud water
    q_clp = q_clp - dq_col_a !   dq_hg_eme_gc
    n_clp = n_clp + dn_col_b ! + dn_cl_rime_hg

    ! updating thermodynamic
    ! qtp : removal of droplets
    qtpmcr = qtpmcr - dq_col_a ! dq_hg_eme_gc

    ! thlp : melting and adding liquid water
    ! thlpmcr = thlpmcr+(rlvi/cp_exnf_k)*dq_hg_eme_gc
    thlpmcr = thlpmcr + &
      (rlme/cp_exnf_k)*dq_hg_eme_gc+(rlvi/cp_exnf_k)*dq_col_a
  endif

  if (l_sb_dbg) then
    if(q_clm-delt*dq_hg_rime.lt. 0.0) then
      write(6,*) 'WARNING: coll_gcg3 too high'
      write(6,*) ' removing more cloud water then available'
      ! count(q_clm-delt*dq_hg_rime).lt. 0.0)
      write(6,*) ' removing too much water'
      ! count(q_cl-delt*dq_hg_rime.lt. 0.0)
      write(6,*) ' getting negative q_t'
      ! count(qt0+delt*q_clp).lt. 0.0)
    endif

    if(n_clm+delt*dn_cl_rime_hg.lt. 0.0) then
      write(6,*) 'WARNING: coll_gcg3 too high'
      write(6,*) ' removing more droplets then available'
      ! count(n_clm+delt*dn_cl_rime_hg.lt. 0.0)
      write(6,*) ' removing too many droplets'
      ! count(n_cl+delt*dn_cl_rime_hg.lt. 0.0)
    endif
  endif
  if (l_tendencies) then
    tend(idn_hg_eme_gc ) = dn_hg_eme_gc
    tend(idq_hg_eme_gc ) = dq_hg_eme_gc
    tend(idq_hg_rime   ) = dq_hg_rime
    tend(idn_cl_rime_hg) = dn_cl_rime_hg
  endif
  if (l_statistics) then
    statistics(imphys_freeze) = statistics(imphys_freeze) + dq_hg_rime
    statistics(imphys_melt) = statistics(imphys_melt) + dq_hg_eme_gc
  endif
end subroutine coll_gcg3


! ****************************************
!  riming of graupel
!
!
! ****************************************
subroutine coll_grg3
  use modglobal, only : pi
  implicit none

  real :: E_ab
  real :: k_enhm
  real :: rem_cf
  real :: dn_col_b, dq_col_a

  real :: dn_hr_rime_hg = 0.
  real :: dn_hg_eme_gr  = 0.   !< number tendency enhanced melting of graupel by rain
  real :: dq_hg_eme_gr  = 0.   !< mass tendency enhanced melting of graupel by rain

  E_ab =  E_er_m

  dq_col_a = (rhof_k*pi/4)*E_ab*n_hg    &
        *q_hr*(dlt_g0*D_hg**2                     &
          +dlt_g1r*D_hg*D_hr+dlt_r1*D_hr**2) &
        *( th_g0*v_hg**2-th_g1r*v_hr*v_hg    &
          +th_r1*v_hr**2+sigma_hg**2+sigma_hr**2)**0.5

  dn_col_b = -(rhof_k*pi/4)*E_ab*n_hg   &
        *n_hr*(dlt_g0*D_hg**2                     &
          +dlt_g0r*D_hg*D_hr+dlt_r0*D_hr**2) &
        *( th_g0*v_hg**2-th_g0r*v_hr*v_hg    &
          +th_r0*v_hr**2+sigma_hg**2+sigma_hr**2)**0.5

  if (tmp0.lt.T_3) then
    ! riming only
    ! remaining number of particles
    rem_cf = (1.0-rem_n_hr_min)/delt

    dq_hghr_rime = min(dq_col_a,&
                       max(0.0,&
                       q_hrm/delt+q_hrp))   ! following ICON, 2017
    ! = min(dq_hghr_rime,q_hrm/delt)

    dn_hr_rime_hg = max(dn_col_b,&
                        min(0.0,&
                        -rem_cf*n_hrm-n_hrp))
    ! = min(dn_hr_rime_hg,-n_hrm/delt)

    ! change in the amount of graupel
    q_hgp = q_hgp + dq_hghr_rime

    ! change in the amount of rain
    n_hrp = n_hrp + dn_hr_rime_hg
    q_hrp = q_hrp - dq_hghr_rime

    ! no change in q_t
    ! qtpmcr = qtpmcr - dq_col_a

    ! change in th_l - just heat release from freezing
    thlpmcr = thlpmcr+ (rlme/cp_exnf_k)*dq_hghr_rime
  else
    ! not riming,but enhanced melting and scheding

    ! enhanced melting
    k_enhm  = c_water/rlme

    ! calculating the melting
    dq_hg_eme_gr = -k_enhm*(tmp0-T_3)*min(dq_col_a,q_hr/delt)
    dq_hg_eme_gr = max(dq_hg_eme_gr,min(0.0,-q_hgm/delt-q_hgp))

    ! calculating number of melted particles
    ! - expected to be proportional to melted mass
    dn_hg_eme_gr = dq_hg_eme_gr*n_hg/q_hg

    ! - but not more than number of interacting particles
    ! dn_hg_eme_gr = max(dn_hg_eme_gr, max(dn_col_b,-n_hr/delt))
    ! - and not more than total number of particles
    dn_hg_eme_gr = max(dn_hg_eme_gr,&
                       -min((n_hrm/delt+n_hrp),&
                            (n_hgm/delt+n_hgp)))

    ! updating tendencies
    ! updating rain
    q_hrp = q_hrp - dq_hg_eme_gr
    ! n_hrp = n_hrp + dn_hr_rime_hg

    ! graupel
    q_hgp = q_hgp + dq_hg_eme_gr
    n_hgp = n_hgp + dn_hg_eme_gr

    ! updating cloud water
    ! no change - all turned into rain
    ! updating thermodynamic
    ! qt : no change -  increased by melted water goes to rain
    ! qtpmcr = qtpmcr +0.0
    ! thl : melting
    thlpmcr = thlpmcr+(rlme/cp_exnf_k)*dq_hg_eme_gr
  endif

  if (l_sb_dbg) then
    if(q_hrm-delt*dq_col_a.lt. 0.0) then
      write(6,*) 'WARNING: coll_grg3 too high'
      write(6,*) ' removing more rain water than available'
      ! count(q_hrm-delt*dq_hghr_rime.lt. 0.0)
      write(6,*) ' removing too much rain water '
      ! count(q_hr-delt*dq_hghr_rime.lt. 0.0)
      write(6,*) ' getting negative q_t in  '
      ! count(qt0+delt*q_hrp).lt. 0.0)
    endif
    if(n_hrm+delt*dn_col_b.lt. 0.0) then
      write(6,*) 'WARNING: coll_grg too high'
      write(6,*) ' removing more raindrops than available '
      ! count(n_hrm+delt*dn_hr_rime_hg.lt. 0.0)
      write(6,*) ' removing too many raindrops'
      ! count(n_hr+delt*dn_hr_rime_hg.lt. 0.0)
    endif
  endif
  if (l_tendencies) then
    tend(idn_hr_rime_hg) = dn_hr_rime_hg
    tend(idn_hg_eme_gr ) = dn_hg_eme_gr
    tend(idq_hg_eme_gr ) = dq_hg_eme_gr
    tend(idq_hghr_rime ) = dq_hghr_rime
  endif
  if (l_statistics) then
    statistics(imphys_freeze) = statistics(imphys_freeze) + dq_hghr_rime
    statistics(imphys_melt) = statistics(imphys_melt) + dq_hg_eme_gr
  endif
endsubroutine coll_grg3


! **************************************
! Heterogeneous freezing of rain
! **************************************
subroutine rainhetfreez3
  implicit none

  real :: J_het
  real :: dq_hr_het = 0. !< heterogeneou freezing of raindrops
  real :: dn_hr_het = 0. !< heterogeneou freezing of raindrops

  ! maybe only for temperatures below T_3 ?
  J_het = A_het *exp( B_het*(T_3-tmp0) -1)

  dn_hr_het = -c_mmt_1hr * n_hr * x_hr * J_het
  dq_hr_het = -c_mmt_2hr * q_hr * x_hr * J_het

  ! basic correction
  dn_hr_het = max(dn_hr_het,min(0.0,-n_hrm/delt-n_hrp))
  dq_hr_het = max(dq_hr_het,min(0.0,-q_hrm/delt-q_hrp))

  ! decrease in raindrops
  n_hrp = n_hrp + dn_hr_het
  q_hrp = q_hrp + dq_hr_het

  ! increase in graupel
  n_hgp = n_hgp - dn_hr_het
  q_hgp = q_hgp - dq_hr_het

  ! and consumption of aerosols for heterogeneous freezing  ?
  ! n_ccp = n_ccp - dn_hr_het

  ! no changes in the total amount of water
  ! qtpmcr = qtpmcr

  ! change in th_l due to freezing
  thlpmcr = thlpmcr - (rlme/cp_exnf_k)*dq_hr_het

  if (l_tendencies) then
    tend(idq_hr_het) = dq_hr_het
    tend(idn_hr_het) = dn_hr_het
  endif
end subroutine rainhetfreez3


! riming of ice + rain to graupel
! -------------------------------
subroutine coll_rig3
  use modglobal, only : pi
  implicit none

  real :: rem_ci_cf,rem_hr_cf,k_enhm
  real :: E_ab, dn_col_b, dq_col_a, dq_col_b

  real :: dn_ci_col_ri  = 0.  !< ice number loss from riming of ice+rain->gr
  real :: dn_hr_col_ri  = 0.  !< rain number loss from riming of ice+rain->gr
  real :: dn_ci_eme_ri  = 0.  !< number tendency enhanced melting of cloud ice by rain
  real :: dq_ci_eme_ri  = 0.  !< mass tendency enhanced melting of cloud ice  by rain

  ! setting up extra coefficients
  ! remain coefficients
  rem_ci_cf = (1.0-rem_n_ci_min)/delt
  rem_hr_cf = (1.0-rem_n_hr_min)/delt

  ! enhanced melting
  k_enhm  = c_water/rlme

  E_ab =  E_i_m

  dq_col_b = -(rhof_k*pi/4)*E_ab*n_ci           &
       *q_hr*(dlt_i0*D_ci**2                     &
         +dlt_i1r*D_ci*D_hr+dlt_r1*D_hr**2)      &
       *( th_i0*v_ci**2-th_i1r*v_hr*v_ci         &
         +th_r1*v_hr**2+sigma_ci**2+sigma_hr**2)**0.5

  dn_col_b = -(rhof_k*pi/4)*E_ab*n_ci           &
       *n_hr*(dlt_i0*D_ci**2                     &
         +dlt_i0r*D_ci*D_hr+dlt_r0*D_hr**2)      &
       *( th_i0*v_ci**2-th_i0r*v_hr*v_ci         &
         +th_r0*v_hr**2+sigma_ci**2+sigma_hr**2)**0.5

  dq_col_a = (rhof_k*pi/4)*E_ab*q_ci            &
       *n_hr*(dlt_i1*D_ci**2                     &
         +dlt_r1i*D_ci*D_hr+dlt_r0*D_hr**2)      &
       *( th_i1*v_ci**2-th_r1i*v_hr*v_ci         &
         +th_r0*v_hr**2+sigma_ci**2+sigma_hr**2)**0.5

  dq_hr_col_ri =  dq_col_b
  dn_ci_col_ri =  dn_col_b
  dn_hr_col_ri =  dn_col_b
  dq_ci_col_ri = -dq_col_a


  ! first adjustment
  dq_ci_col_ri = max(dq_ci_col_ri,min(0.0,-q_cim/delt-q_cip))   ! following ICON, 2017
  dq_hr_col_ri = max(dq_hr_col_ri,min(0.0,-q_hrm/delt-q_hrp))   ! following ICON, 2017

  ! adjustment of numbers - both ice and water
  dn_ci_col_ri = max(dn_ci_col_ri,min(0.0,-rem_ci_cf*n_cim-n_cip))
  dn_ci_col_ri = max(dn_ci_col_ri,min(0.0,-rem_hr_cf*n_hrm-n_hrp))

  if(tmp0.lt.T_3) then
    ! the collection is just riming
    ! decrease in numeber of raindrops same as decrease in number of ice
    dn_hr_col_ri = dn_ci_col_ri

    ! record the change in cloud ice
    q_cip = q_cip + dq_ci_col_ri
    n_cip = n_cip + dn_ci_col_ri

    ! change in rain
    q_hrp = q_hrp + dq_hr_col_ri
    n_hrp = n_hrp + dn_hr_col_ri

    ! and for graupel
    q_hgp = q_hgp - dq_ci_col_ri - dq_hr_col_ri
    n_hgp = n_hgp - dn_ci_col_ri

    ! change in q_t - decrease in cloud ice
    qtpmcr = qtpmcr + 0.0 ! dq_ci_col_ri

    ! change in th_l - release from freezing and removal of ice
    thlpmcr = thlpmcr - (rlme/cp_exnf_k)*dq_hr_col_ri
  else  ! tmp0.gt.T_3
    ! enhanced melting and graupel formation

    ! calculate the melting
    dq_ci_eme_ri = k_enhm*(tmp0-T_3)*dq_hr_col_ri ! with + due to negative value of dq_hr here

    ! limit melting
    dq_ci_eme_ri = max(dq_ci_eme_ri,min(0.0,-q_cim/delt-q_cip))

    ! calculate how many ice perticles melted
    dn_ci_eme_ri = dq_ci_eme_ri*n_ci/q_ci

    ! limit so it dos not melt more than interacting
    ! also limit so that new graupel not larger that max mean size of source ice ?
    dn_ci_eme_ri =max(dn_ci_eme_ri,dn_ci_col_ri)

    ! update ice
    q_cip = q_cip + dq_ci_eme_ri ! dq_ci_col_ri
    n_cip = n_cip + dn_ci_eme_ri ! dn_ci_col_ri

    ! no collection of raindrops
    dq_hr_col_ri = 0.0
    dq_ci_col_ri = 0.0
    dn_hr_col_ri = 0.0
    dn_ci_col_ri = 0.0

    ! increase in rain mass
    q_hrp = q_hrp - dq_ci_eme_ri

    ! change in thl - heat spent on melting
    thlpmcr = thlpmcr+(rlme/cp_exnf_k)*dq_ci_eme_ri
  endif
  if (l_tendencies) then
    tend(idn_ci_col_ri) = dn_ci_col_ri
    tend(idq_ci_col_ri) = dq_ci_col_ri
    tend(idn_hr_col_ri) = dn_hr_col_ri
    tend(idq_hr_col_ri) = dq_hr_col_ri
    tend(idn_ci_eme_ri) = dn_ci_eme_ri
    tend(idq_ci_eme_ri) = dq_ci_eme_ri
  endif
  if (l_statistics) then
    statistics(imphys_freeze) = statistics(imphys_freeze) - dq_hr_col_ri
    statistics(imphys_melt) = statistics(imphys_melt) + dq_ci_eme_ri
  endif
end subroutine coll_rig3


! riming of rain + snow to graupel
! --------------------------------
subroutine coll_rsg3
  use modglobal, only : pi
  implicit none

  real :: rem_cf, k_enhm
  real :: E_ab, dn_col_b, dq_col_a, dq_col_b

  real :: dn_hr_col_rs  = 0. !< snow loss from riming of ice+snow->gr
  real :: dn_hs_col_rs  = 0. !< snow number loss from riming of ice+snow->gr
  real :: dn_hs_eme_rs  = 0. !< number tendency enhanced melting of snow by rain
  real :: dq_hs_eme_rs  = 0. !< mass tendency enhanced melting of snow by rain

  dn_col_b = 0.0
  dq_col_a = 0.0
  dq_col_b = 0.0

  E_ab =  E_s_m

  dq_col_b = -(rhof_k*pi/4)*E_ab*n_hs           &
       *q_hr*(dlt_s0*D_hs**2                     &
         +dlt_s1r*D_hs*D_hr+dlt_r1*D_hr**2)      &
       *( th_s0*v_hs**2-th_s1r*v_hr*v_hs         &
         +th_r1*v_hr**2+sigma_hs**2+sigma_hr**2)**0.5

  dn_col_b = -(rhof_k*pi/4)*E_ab*n_hs           &
       *n_hr*(dlt_s0*D_hs**2                     &
         +dlt_s0r*D_hs*D_hr+dlt_r0*D_hr**2)      &
       *( th_s0*v_hs**2-th_s0r*v_hr*v_hs         &
         +th_r0*v_hr**2+sigma_hs**2+sigma_hr**2)**0.5

  dq_col_a = (rhof_k*pi/4)*E_ab*q_hs            &
       *n_hr*(dlt_s1*D_hs**2                     &
         +dlt_r1s*D_hs*D_hr+dlt_r0*D_hr**2)      &
       *( th_s1*v_hs**2-th_r1s*v_hr*v_hs         &
         +th_r0*v_hr**2+sigma_hs**2+sigma_hr**2)**0.5

  ! first adjustment following ICON, 2017
  dq_hs_col_rs = max(-dq_col_a,min(0.0,-q_hsm/delt-q_hsp))
  dq_hr_col_rs = max( dq_col_b,min(0.0,-q_hrm/delt-q_hrp))

  ! adjustment of numbers - both ice and snow
  rem_cf = (1.0-rem_n_hs_min)/delt
  dn_hs_col_rs = max(dn_col_b,min(0.0,-rem_cf*n_hsm-n_hsp))

  rem_cf = (1.0-rem_n_hr_min)/delt
  dn_hs_col_rs = max(dn_col_b,min(0.0,-rem_cf*n_hrm-n_hrp))

  if (tmp0.lt.T_3) then
    ! and copying it to the second one
    dn_hr_col_rs = dn_hs_col_rs

    ! record the change in cloud ice
    q_hsp = q_hsp + dq_hs_col_rs
    n_hsp = n_hsp + dn_hs_col_rs

    ! change in rain
    q_hrp = q_hrp + dq_hr_col_rs
    n_hrp = n_hrp + dn_hr_col_rs

    ! and for graupel
    q_hgp = q_hgp - dq_hs_col_rs - dq_hr_col_rs
    n_hgp = n_hgp - dn_hs_col_rs

    ! change in th_l - release from freezing and removal of ice
    thlpmcr = thlpmcr - (rlme/cp_exnf_k)*dq_hr_col_rs
  else  ! tmp0.gt.T_3
    ! enhanced melting and graupel formation
    k_enhm  = c_water/rlme

    ! calculate the melting
    ! with + due to negative value of dq_hr here
    dq_hs_eme_rs = k_enhm*(tmp0-T_3)*dq_hr_col_rs

    ! snow melting
    dq_hs_eme_rs = max(dq_hs_eme_rs,min(0.0,-q_hsm/delt-q_hsp))

    ! calculate how many snow perticles melted
    ! q_hs here is always some small positive number
    dn_hs_eme_rs = dq_hs_eme_rs*n_hs/q_hs

    ! limit so it dos not melt more than interacting
    ! also limit so that new graupel not larger that max mean size of source ice ?
    dn_hs_eme_rs = max(dn_hs_eme_rs,dn_hs_col_rs)

    ! update snow
    q_hsp = q_hsp + dq_hs_eme_rs
    n_hsp = n_hsp + dn_hs_eme_rs

    ! no collection of raindrops
    dq_hr_col_rs = 0.0
    dq_hs_col_rs = 0.0
    dn_hr_col_rs = 0.0
    dn_hs_col_rs = 0.0

    ! increase in rain mass
    q_hrp = q_hrp - dq_hs_eme_rs

    ! change in thl - heat spent on melting
    thlpmcr = thlpmcr+(rlme/cp_exnf_k)*dq_hs_eme_rs
  endif
  if (l_tendencies) then
    tend(idn_hr_col_rs) = dn_hr_col_rs
    tend(idn_hs_col_rs) = dn_hs_col_rs
    tend(idn_hs_eme_rs) = dn_hs_eme_rs
    tend(idq_hs_eme_rs) = dq_hs_eme_rs
    tend(idq_hr_col_rs) = dq_hr_col_rs
    tend(idq_hs_col_rs) = dq_hs_col_rs
  endif
  if (l_statistics) then
    statistics(imphys_freeze) = statistics(imphys_freeze) - dq_hr_col_rs
    statistics(imphys_melt) = statistics(imphys_melt) + dq_hs_eme_rs
  endif
end subroutine coll_rsg3


! Ice multiplication
!   - calls separately H-M process for each of the riming processes
!
! conditions: tmp0.lt.T_3
!
! ------------------
subroutine ice_multi3
  implicit none

  real :: dq_ci_spl, dn_ci_spl, dq_hg_temp
  real :: dn_ci_mul = 0.  !< ice multiplication
  real :: dq_ci_mul = 0.  !< ice multiplication

  dq_hg_temp = 0.0
  dq_ci_spl  = 0.0
  dn_ci_spl  = 0.0

  ! ------------------------------------
  ! calling separate H-M processes
  ! ------------------------------------

  ! i+l -> i
  call hallet_mossop3 (tmp0,dq_ci_rime,q_ci,dq_ci_spl,dn_ci_spl)
  dn_ci_mul = dn_ci_mul + dn_ci_spl

  ! s+l -> s
  call hallet_mossop3 (tmp0,dq_hs_rime,q_hs,dq_ci_spl,dn_ci_spl)
  dn_ci_mul = dn_ci_mul + dn_ci_spl
  dq_ci_mul = dq_ci_mul + dq_ci_spl
  q_hsp     = q_hsp - dq_ci_spl ! effect on snow

  ! g+l -> g
  call hallet_mossop3 (tmp0,dq_hg_rime,q_hg,dq_ci_spl,dn_ci_spl)
  dn_ci_mul = dn_ci_mul + dn_ci_spl
  dq_ci_mul = dq_ci_mul + dq_ci_spl
  q_hgp     = q_hgp - dq_ci_spl ! effect on snow

  ! g+r -> g
  call hallet_mossop3 (tmp0,dq_hghr_rime,q_hg,dq_ci_spl,dn_ci_spl)
  dn_ci_mul = dn_ci_mul + dn_ci_spl
  dq_ci_mul = dq_ci_mul + dq_ci_spl
  q_hgp     = q_hgp - dq_ci_spl  ! effect on snow

  ! s+r -> g  -- does it make sense to include it?
  dq_hg_temp = -dq_hs_col_rs - dq_hr_col_rs
  call hallet_mossop3 (tmp0,dq_hg_temp ,q_hg,dq_ci_spl,dn_ci_spl)
  dn_ci_mul = dn_ci_mul + dn_ci_spl
  dq_ci_mul = dq_ci_mul + dq_ci_spl
  q_hgp     = q_hgp - dq_ci_spl ! effect on snow

  ! i+r-> g   -- does it make sense to include it?
  dq_hg_temp = -dq_ci_col_ri - dq_hr_col_ri
  call hallet_mossop3 (tmp0,dq_hg_temp,q_hg,dq_ci_spl,dn_ci_spl)
  dn_ci_mul = dn_ci_mul + dn_ci_spl
  dq_ci_mul = dq_ci_mul + dq_ci_spl
  q_hgp     = q_hgp - dq_ci_spl  ! effect on snow

  ! ------------------------------------
  !   update
  ! ------------------------------------
  ! add updates to dq_ci, dn_ci
  n_cip = n_cip + dn_ci_mul
  q_cip = q_cip + dq_ci_mul

  if (l_tendencies) then
    tend(idn_ci_mul) = dn_ci_mul
    tend(idq_ci_mul) = dq_ci_mul
  endif
end subroutine ice_multi3


! ****************************************
!  partial conversion to graupel
!
!  - based on Seifgert and Beheng, 2004
!       - the implementation differs from the implementation in ICON
!
!
!  - have to be called after enhanced melting subroutine
!
! ****************************************
subroutine conv_partial3
  use modglobal, only : pi,rhow
  use modmicrodata3, only : rhoeps,al_0ice,al_0snow
  implicit none

  real, parameter  :: pi6rhoe = (pi/6.0)*rhoeps    &
                     ,cc_ci = al_0ice*rhow/rhoeps  &
                     ,cc_hs = al_0snow*rhow/rhoeps

  real :: rem_ci_cf, rem_hs_cf

  real :: dq_ci_cv = 0. !< partial conversion ice -> graupel
  real :: dn_ci_cv = 0.
  real :: dq_hs_cv = 0. !< partial conversion snow-> graupel
  real :: dn_hs_cv = 0.

  ! remain coefficients
  rem_ci_cf = (1.0-rem_n_min_cv)/delt
  rem_hs_cf = (1.0-rem_n_min_cv)/delt

  if (q_ci_mask.and.D_ci.gt.D_mincv_ci) then
    ! remain coefficient
    rem_ci_cf = (1.0-rem_n_min_cv)/delt

    dq_ci_cv = cc_ci*(pi6rhoe*D_ci**3 /x_ci-1.0)
    ! ? not to exceed conversion rate 1 ?
    ! G_ci = min(1.0, G_ci)

    dq_ci_cv = -dq_ci_rime/dq_ci_cv
    dq_ci_cv = max(dq_ci_cv, min(0.0, -q_cim/delt-q_cip))  ! based on ICON, 2017
    ! = max(dq_ci_cv,-q_cim/delt)

    dn_ci_cv = dq_ci_cv/max(x_ci,x_ci_cvmin)
    ! = dq_ci_cv/max(x_ci,x_ci_cvmin)

    dn_ci_cv = max(dn_ci_cv, min(0.0, -rem_ci_cf*n_cim-n_cip))
    ! = max(dn_hs_cv, -n_cim/delt)

    ! change in the amount of graupel
    n_hgp = n_hgp - dn_ci_cv
    q_hgp = q_hgp - dq_ci_cv

    ! change in the amount of cloud ice
    n_cip = n_cip + dn_ci_cv
    q_cip = q_cip + dq_ci_cv
  endif

  ! term for snow conversion
  if (q_hs_mask.and.D_hs.gt.D_mincv_hs) then
    ! remain coefficient
    rem_hs_cf = (1.0-rem_n_min_cv)/delt

    dq_hs_cv = cc_hs*(pi6rhoe*D_hs**3 /x_hs-1)
    ! ? at the sam time, the value should be limited
    ! ? not to exceed conversion rate 1
    ! G_hs = min(1.0, G_hs)
    dq_hs_cv = -dq_hs_rime/dq_hs_cv

    ! correction - not removing more than available
    ! dq_hs_cv = max(dq_hs_cv,-q_hs/delt )

    ! basic correction of the tendency
    dq_hs_cv = max(dq_hs_cv,&
                   min(0.0, &
                   -q_hsm/delt-q_hsp))
    ! dq_hs_cv = max( dq_hs_cv,-q_hsm/delt)

    dn_hs_cv = dq_hs_cv/max(x_hs,x_hs_cvmin)
    ! dn_hs_cv = dq_hs_cv/max(x_hs,x_hs_cvmin)

    dn_hs_cv = max(dn_hs_cv,min(0.0,-rem_hs_cf*n_hsm-n_hsp))
    ! dn_hs_cv = max(dn_hs_cv, -n_hsm/delt)

    ! and the second correction of the q tendency
    ! change in the amount of graupel
    n_hgp = n_hgp - dn_hs_cv
    q_hgp = q_hgp - dq_hs_cv

    ! change in the amount of snow
    n_hsp = n_hsp + dn_hs_cv
    q_hsp = q_hsp + dq_hs_cv
  endif

  ! warnings
  if (l_sb_dbg) then
    if(q_cim+delt*dq_ci_cv.lt. 0) then
      write(6,*) 'WARNING: conv_partial3 removing too much ice', delt*dq_ci_cv, '/', q_cim
    endif
    if(q_hsm+delt*dq_hs_cv.lt. 0) then
      write(6,*) 'WARNING: conv_partial3 removing too much snow', delt*dq_hs_cv, '/', q_hsm
    endif
  endif
  if (l_tendencies) then
    tend(idq_ci_cv) = dq_ci_cv
    tend(idn_ci_cv) = dn_ci_cv
    tend(idq_hs_cv) = dq_hs_cv
    tend(idn_hs_cv) = dn_hs_cv
  endif
end subroutine conv_partial3


! ***************************************************************
! Melting and evaporation of ice particles
!
! - wrapper
!   - calls separately melting processes for each of species
!
! conditions: (tmp0.gt.T_3)
!
! ***************************************************************
subroutine evapmelting3
  use modglobal, only : rlv
  implicit none

  real :: dn_ci_me = 0.   !< number tendency melting of cloud ice
  real :: dq_ci_me = 0.   !< mass tendency melting of cloud ice
  real :: dn_hs_me = 0.   !< number tendency melting of snow
  real :: dq_hs_me = 0.   !< mass tendency melting of snow
  real :: dn_hg_me = 0.   !< number tendency melting of graupel
  real :: dq_hg_me = 0.   !< mass tendency melting of graupel
  real :: dn_ci_ev = 0.   !< number tendency evaporation of cloud ice
  real :: dq_ci_ev = 0.   !< mass tendency evaporation of cloud ice
  real :: dn_hs_ev = 0.   !< number tendency evaporation of snow
  real :: dq_hs_ev = 0.   !< mass tendency evaporation of snow
  real :: dn_hg_ev = 0.   !< number tendency evaporation of graupel
  real :: dq_hg_ev = 0.   !< mass tendency evaporation of graupel

  ! ------------------------------------
  ! calling separate melting processes
  ! ------------------------------------

  ! ice
  if (q_ci.gt.qicemin) then
    call sb_evmelt3(aven_0i,aven_1i,bven_0i,bven_1i,x_ci_bmin        &
                   ,n_ci,n_cip,n_cim,q_ci,q_cip,q_cim                &
                   ,x_ci,D_ci,v_ci,dq_ci_me,dn_ci_me,dq_ci_ev,dn_ci_ev)

    ! loss of ice
    n_cip  = n_cip + dn_ci_me + dn_ci_ev
    q_cip  = q_cip + dq_ci_me + dq_ci_ev

    ! transformed to rain or cloud
    ! for now into rain - improve possibly later
    n_hrp = n_hrp - dn_ci_me
    q_hrp = q_hrp - dq_ci_me

    ! transfomed to water vapour
    qtpmcr = qtpmcr - dq_ci_ev
    ! and heat production : heat spent on melting and evaporation - done lower
  endif

  ! snow
  if (q_hs.gt.qsnowmin) then
    call sb_evmelt3(aven_0s,aven_1s,bven_0s,bven_1s,x_hs_bmin        &
                   ,n_hs,n_hsp,n_hsm,q_hs,q_hsp,q_hsm                &
                   ,x_hs,D_hs,v_hs,dq_hs_me,dn_hs_me,dq_hs_ev,dn_hs_ev)

    ! loss of snow
    n_hsp = n_hsp + dn_hs_me +dn_hs_ev
    q_hsp = q_hsp + dq_hs_me +dq_hs_ev

    ! transformed to rain
    n_hrp = n_hrp - dn_hs_me
    q_hrp = q_hrp - dq_hs_me

    ! transfomed to water vapour
    qtpmcr = qtpmcr-dq_hs_ev
    ! and heat production : heat spent on melting and evaporation - done lower
  endif

  ! graupel
  if (q_hg.gt.qgrmin) then
    call sb_evmelt3(aven_0g,aven_1g,bven_0g,bven_1g,x_hg_bmin         &
                   ,n_hg,n_hgp,n_hgm,q_hg,q_hgp,q_hgm                 &
                   ,x_hg,D_hg,v_hg,dq_hg_me,dn_hg_me,dq_hg_ev,dn_hg_ev)

    ! loss of graupel
    n_hgp = n_hgp + dn_hg_me + dn_hg_ev
    q_hgp = q_hgp + dq_hg_me + dq_hg_ev

    ! transformed to rain
    n_hrp = n_hrp - dn_hg_me
    q_hrp = q_hrp - dq_hg_me

    ! transfomed to water vapour
    qtpmcr = qtpmcr-dq_hg_ev

    ! and heat production : heat spent on melting and evaporation - done lower
    !  - melting goes to rain
    !  - evaporation goes water vapour
    thlpmcr = thlpmcr +                                       &
              (rlme/cp_exnf_k)*(dq_ci_me+dq_hs_me+dq_hg_me) + &
        ((rlv+rlme)/cp_exnf_k)*(dq_ci_ev+dq_hs_ev+dq_hg_ev)
  endif

  if (l_tendencies) then
    tend(idn_ci_me) =  dn_ci_me
    tend(idq_ci_me) =  dq_ci_me
    tend(idn_hs_me) =  dn_hs_me
    tend(idq_hs_me) =  dq_hs_me
    tend(idn_hg_me) =  dn_hg_me
    tend(idq_hg_me) =  dq_hg_me
    tend(idn_ci_ev) =  dn_ci_ev
    tend(idq_ci_ev) =  dq_ci_ev
    tend(idn_hs_ev) =  dn_hs_ev
    tend(idq_hs_ev) =  dq_hs_ev
    tend(idn_hg_ev) =  dn_hg_ev
    tend(idq_hg_ev) =  dq_hg_ev
  endif
  if (l_statistics) then
    statistics(imphys_melt) = statistics(imphys_melt) &
      + dq_ci_me + dq_hs_me + dq_hg_me &
      + dq_ci_ev + dq_hs_ev + dq_hg_ev
    statistics(imphys_ev) = statistics(imphys_ev) + dq_hg_ev
  endif
end subroutine evapmelting3


!> Determine autoconversion rate and adjust qrp and Nrp accordingly
!!
!!   The autoconversion rate is formulated for f(x)=A*x**(nuc)*exp(-Bx),
!!   decaying exponentially for droplet mass x.
!!   It can easily be reformulated for f(x)=A*x**(nuc)*exp(-Bx**(mu)) and
!!   by chosing mu=1/3 one would get a gamma distribution in drop diameter
!!   -> faster rain formation. (Seifert)
subroutine autoconversion3
  use modglobal, only : rlv
  implicit none

  real :: rem_cf
  real :: tau, phi, nuc
  real :: dq_hr_au = 0.   !< change in mass of raindrops due to autoconversion
  real :: dn_hr_au = 0.   !< change in number of raindrops due to autoconversion

  if (l_sb) then
    !
    ! SB autoconversion
    !
    if (l_sb_classic) then  ! l_sb_classic - ie. S&B version
      ! autoconversion coefficient
      k_au = k_cc/(20.0*x_s)

      ! remain coefficient
      rem_cf = (1.0-rem_n_cl_min)/delt

      dq_hr_au = k_au * ( nu_cl_cst+2.) * ( nu_cl_cst +4.) / ( nu_cl_cst +1.)**2.    &
              * (q_cl * x_cl)**2. * rho0s ! *rho**2/rho/rho (= 1)
      tau = 1.0 - q_cl / (q_cl + q_hr) ! NOTE: was qltot

      ! phi_au computation
      phi = k_1 * tau**k_2 * (1.0 -tau**k_2)**3
      dq_hr_au = dq_hr_au * (1.0 + phi/(1.0 -tau)**2)

      ! basic au correction
      dq_hr_au  = min(dq_hr_au,max(0.0,q_clm/delt+q_clp))

      ! and calculate cloud droplet numbers
      dn_cl_au = (-2.0/x_s)*dq_hr_au                                   ! and the droplet number
      dn_cl_au = max(dn_cl_au,min(0.0,-rem_cf*n_clm-n_clp)) ! correct droplet number
      dn_hr_au = -0.5*dn_cl_au                                         ! and the raindrop number

      ! and adding updates
      q_hrp = q_hrp + dq_hr_au
      n_hrp = n_hrp + dn_hr_au
      q_clp = q_clp - dq_hr_au
      n_clp = n_clp + dn_cl_au

      thlpmcr = thlpmcr + (rlv/cp_exnf_k)*dq_hr_au
      qtpmcr  = qtpmcr - dq_hr_au
    else ! l_sb_classic = .false. - original version in DALES
      k_au = k_c/(20.0*x_s)

      nuc = k1nuc*(rhof_k*q_cl*1000.) +k2nuc-1. ! #sb3 G09a
      dq_hr_au = k_au * (nuc+2.) * (nuc+4.) / (nuc+1.)**2.    &
              * (q_cl * x_cl)**2. * rho0s ! *rho**2/rho/rho (= 1)
      tau            = 1.0 - q_cl / (q_cl + q_hr) ! NOTE: was qltot
      phi = k_1 * tau**k_2 * (1.0 -tau**k_2)**3   ! phi_au computation
      dq_hr_au = dq_hr_au * (1.0 + phi/(1.0 -tau)**2)

      ! cloud water numbers
      dn_hr_au = dq_hr_au/x_s
      dn_cl_au = (-2.0/x_s)*dq_hr_au

      ! #sb3 START outputs
      q_hrp = q_hrp + dq_hr_au
      n_hrp = n_hrp + dn_hr_au
      q_clp = q_clp - dq_hr_au
      n_clp = n_clp + dn_cl_au  ! o:  n_clp - (2.0/x_s)*rhof_k*au

      thlpmcr = thlpmcr + (rlv/cp_exnf_k)*dq_hr_au
      qtpmcr  = qtpmcr  - dq_hr_au
    endif ! l_sb_classic
  endif ! l_sb

  if (l_sb_dbg) then
    if (q_cl/delt - dq_hr_au.lt. 0.) then
      write(6,*) 'WARNING: autoconversion too high'
      write(6,*) '  removing more cloud water than available'
      ! count(q_cl/delt - dq_hr_au.lt. 0.)
    end if

    if (n_cl/delt -(2.0/x_s)*dq_hr_au.lt. 0.) then
      write(6,*) 'WARNING: autoconversion too high'
      write(6,*) '  removing more droplets than available'
      ! count(n_cl(2:i1,2:j1,1:kmax)/delt - (2.0/x_s)*dq_hr_au.lt. 0.)
    end if
  endif
  if (l_tendencies) then
    tend(idq_hr_au) = dq_hr_au
    tend(idn_hr_au) = dn_hr_au
    tend(idn_cl_au) = dn_cl_au
  endif
end subroutine autoconversion3


!> Self-collection of cloud droplets
!! written based on S&B to evaluate the self-collection of cloud droplets
subroutine cloud_self3
  use modmicrodata3, only : k_cc,rho0s
  implicit none

  real, parameter :: k_clsc = - k_cc*rho0s

  real    :: rem_cf

  ! calculate constant
  rem_cf = (1.0-rem_n_cl_min)/delt

  dn_cl_sc = k_clsc*(q_cl**2)* &
    (nu_cl_cst+2.0)/(nu_cl_cst+1.0)-dn_cl_au

  ! basic sc collection
  dn_cl_sc = min(0.0,dn_cl_sc)
  dn_cl_sc = max(dn_cl_sc,min(0.0,-rem_cf*n_clm-n_clp))

  ! update
  n_clp = n_clp+dn_cl_sc

  ! no change to q_cl

  if (l_sb_dbg) then
    ! testing for too high values
    if (n_clm/delt - dn_cl_sc.lt. 0.) then
      write(6,*)'WARNING: cloud sc too high'
      write(6,*) '  getting to negative n_cl'
      ! count(n_clm/delt -dn_cl_sc).lt. 0.)
    endif
  endif
  if (l_tendencies) then
    tend(idn_cl_sc) = dn_cl_sc
  endif
end subroutine cloud_self3


!*********************************************************************
! determine accr. + self coll. + br-up rate and adjust qrp and Nrp
! base don S&B
!
! conditions : q_hr_mask
!*********************************************************************
subroutine accretion3
  use modglobal, only : rlv
  implicit none

  real :: phi, Dvrf, tau
  real :: rem_cf
  real :: dq_hr_ac = 0.
  real :: dn_hr_br = 0.
  real :: dn_hr_sc = 0.

  if (l_sb) then

    ! SB accretion
    if (l_sb_classic) then

      if (q_cl_mask) then
        ! since it is forming only where rain
        tau = 1.0 - q_cl/(q_cl + q_hr) ! NOTE: was qltot
        phi = (tau/(tau + k_l))**4.
        dq_hr_ac = k_cr *rhof_k*q_cl*q_hr * phi * &
                         (rho0s/rhof_k)**0.5  ! rho*rho / rho  = rho

        ! basic ac correction
        dq_hr_ac = min(dq_hr_ac,q_cl/delt) ! min(dq_hr_ac,q_clm/delt)

        ! number of cloud droplets
        ! remain coefficient for clouds #sb3
        rem_cf = (1.0-rem_n_cl_min)/delt

        dn_cl_ac = -dq_hr_ac/x_cl
        dn_cl_ac = max(dn_cl_ac,&
                       min(0.0, &
                       -rem_cf*n_clm-n_clp))

        ! update
        q_hrp = q_hrp + dq_hr_ac

        ! no change n_hrp
        ! and changes in water number for clouds
        q_clp = q_clp - dq_hr_ac
        n_clp = n_clp + dn_cl_ac !o: n_clp - rhof_k*ac/x_cl

        qtpmcr = qtpmcr - dq_hr_ac
        thlpmcr = thlpmcr + (rlv/cp_exnf_k)*dq_hr_ac
      endif ! qcl_mask

      ! SB self-collection & Break-up
      dn_hr_sc = -k_rr *rhof_k* q_hr * n_hr  &
        * (1.0 + kappa_r/lbdr)**(-9.)*(rho0s/rhof_k)**0.5

      ! and calculating size of droplets - adjusted
      Dvrf = Dvr ! for now leaving the same
      if (Dvrf.gt.dvrlim) then
        if (Dvrf.gt.dvrbiglim) then
          ! for big drops
          phi = 2.0*exp(kappa_br*(Dvrf-D_eq))-1.0
        else
          ! for smaller drops
          phi = k_br * (Dvrf-D_eq)
        endif
        dn_hr_br = -(phi + 1.) * dn_hr_sc
      endif

    else ! l_sb_classic
      !*********************************************************************
      ! determine accr. rate and adjust qrp and Nrp
      ! accordingly. Break-up : Seifert (2007), a
      !*********************************************************************
      if (q_cl_mask) then
        tau = 1.0 - ql0/(q_cl + q_hr) ! NOTE: was qltot
        phi = (tau/(tau + k_l))**4.

        dq_hr_ac = k_r *rhof_k*ql0 * q_hr * phi * (1.225/rhof_k)**0.5
        q_hrp = q_hrp + dq_hr_ac

        ! no change n_hrp
        q_clp  = q_clp - dq_hr_ac
        qtpmcr  = qtpmcr  - dq_hr_ac
        thlpmcr = thlpmcr + (rlv/cp_exnf_k)*dq_hr_ac

        ! and update on cloud water number
        dn_cl_ac = -dq_hr_ac/x_cl
        n_clp = n_clp + dn_cl_ac
      endif

      dn_hr_sc = -k_rr *rhof_k* q_hr * n_hr  &
        * (1.0 + kappa_r/lbdr*pirhow**(1./3.))**(-9.)*(rho0s/rhof_k)**0.5

      if (Dvr .gt. dvrlim) then
        if (Dvr .gt. dvrbiglim) then
          ! for big drops
          phi = 2.0*exp(kappa_br*(Dvr-D_eq))-1.0
        else
          ! for smaller drops
          phi = k_br* (Dvr-D_eq)
        endif
        dn_hr_br = -(phi + 1.) * dn_hr_sc
      endif
    endif ! l_sb_classic
  endif ! l_sb

  ! correcting for low values
  ! calculating coef for ratio for minimal remaing number of droplets #sb3
  rem_cf = (1.0-rem_n_hr_min)/delt

  n_hrp = n_hrp +max(min(0.0,-rem_cf*n_hrm-n_hrp),dn_hr_sc+dn_hr_br)

  if (l_sb_dbg) then
    if(q_clm/delt-dq_hr_ac.lt. 0.) then
      write(6,*) 'WARNING: accretion removing too much water'
      write(6,*) ' updated below 0'
      ! count(q_clm/delt-dq_hr_ac.lt. 0.0)
    endif

    if(n_hrm/delt+dn_hr_sc+dn_hr_br.lt. 0.) then
      write(6,*) 'WARNING: self-collection of rain too high'
      write(6,*) ' removing more n_hr than available'
      ! count(n_hrm/delt-dn_hr_sc+dn_hr_br.lt. 0.0)
      write(6,*) ' n_hr updated below 0'
      ! count(n_hr/delt-dn_hr_sc+dn_hr_br.lt. 0.0)
    endif

    if (n_clm/delt - dq_hr_ac/x_cl.lt. 0.) then
      write(6,*)'WARNING: ac too large, removing too many droplets'
      ! count(n_clm/delt-dq_hr_ac/x_cl.lt. 0.)
    endif
  endif ! l_sb_dbg
  if (l_tendencies) then
    tend(idn_cl_ac) = dn_cl_ac
    tend(idq_hr_ac) = dq_hr_ac
    tend(idn_hr_br) = dn_hr_br
    tend(idn_hr_sc) = dn_hr_sc
  endif
end subroutine accretion3


!*********************************************************************
! Evaporation of prec. : Seifert & Beheng
! Cond. (S>0.) neglected (all water is condensed on cloud droplets)
!*********************************************************************
subroutine evap_rain3
  use modglobal, only    : rv,rlv,pi
  implicit none

  real :: f0    & ! ventilation factor - moment 0
         ,f1    & ! ventilation factor - moment 1
         ,S     & ! super or undersaturation
         ,G     & ! cond/evap rate of a drop
         ,vihr  & ! mean terminal velocity
         ,nrex  & ! Reynolds number N_re(xr)
         ,x_hrf   ! full x_hr without bounds

  real :: dq_hr_ev = 0.
  real :: dn_hr_ev = 0.

  ! adjusting the calculation for saturation
  S = min(0.,((qt0-q_cl)/qvsl- 1.0))
  G = (rv * tmp0) / (Dv*esl) &
    + rlv/(Kt*tmp0)*(rlv/(rv*tmp0) -1.)
  G = 1./G

  ! terminal velocity  (from mixed scheme)
  vihr = al_hr*((rho0s/rhof_k)**0.5)*x_hr**be_hr

  ! calculating  N_re Reynolds number
  nrex = Dvr*vihr/nu_a

  ! NOTE: eps0 is limiting here.
  ! x_hr it can lead in case of many big raindrops to removal of too many of them.
  x_hrf = q_hr/(n_hr + eps0)

  f0 = aven_0r+bven_0r*Sc_num**(1.0/3.0)*nrex**0.5
  f1 = aven_1r+bven_1r*Sc_num**(1.0/3.0)*nrex**0.5
  f0 = max(0.0,f0)

  dq_hr_ev = 2*pi*n_hr*G*Dvr*f1*S

  dn_hr_ev = 2*pi*n_hr*G*Dvr*f0*S/x_hrf

  ! and limiting it
  dn_hr_ev = min(dn_hr_ev, 0.0)
  dn_hr_ev = max(dn_hr_ev,dq_hr_ev/x_hrf)

  if ((dq_hr_ev+q_hrm/delt.lt.0).or.&
      (dn_hr_ev+n_hrm/delt.lt.0)) then
    dn_hr_ev = - n_hrm/delt
    dq_hr_ev = - q_hrm/delt
  endif

  q_hrp = q_hrp + dq_hr_ev
  n_hrp = n_hrp + dn_hr_ev
  qtpmcr  = qtpmcr - dq_hr_ev
  thlpmcr = thlpmcr + (rlv/cp_exnf_k)*dq_hr_ev

  ! recovery of aerosols ?
  ret_cc = ret_cc - c_ccn_ev_r*min(0.0,dn_hr_ev)

  if (l_tendencies) then
    tend(idq_hr_ev) = dq_hr_ev
    tend(idn_hr_ev) = dn_hr_ev
    tend(iret_cc  ) = ret_cc
  endif
  if (l_statistics) then
    statistics(imphys_ev) = statistics(imphys_ev) + dq_hr_ev
  endif
end subroutine evap_rain3


! =============================
! saturation adjustment
! ============================
! Performs the saturation adjustment to cloud water
!
! Uses split operator:
!   - calculates how much water is nucleated and consumed
!   - calculates remaining amount of water available for condensation
!   - dumps remaining water into cloud
!
!   In addition, check the average size of cloud droplets
!   and evaporates droplets that proportionally
!
! BUG: threshold is probably incorrect, the comments do not match the code
subroutine satadj3
  implicit none

  real :: n_bmax, cogr_max, ql_res
  real :: dq_cl_sa = 0.  !< saturation adjustment
  real :: dn_cl_sa = 0.  !< change in n_cl due to saturation adjustment

  ! NOTE: threshold for adjustment might be lower then threshold for cloud computations, obviously.
  ! The reason is that we want to perform saturation adjustment even for newly nucleated clouds.
  ! However for shorter timesteps (~0.1 s), we should consider a different approach (future plans)
  if (n_clm + n_clp*delt .gt. 0.0) return

  if (l_sb_all_or) then
    !
    ! remaining water =
    !    + condesable water available
    !    - already condensed
    !    - newly condendsed
    !    - removed by mphys processed
    !
    !  calculating amount of available water
    ql_res=(qt0 - qvsl - q_clm)+delt*(qtpmcr - q_clp)
    dq_cl_sa = ql_res/delt
    if ((q_clm+delt*(q_clp+dq_cl_sa)).lt.0.0) then
      dn_cl_sa = -n_clm/delt-n_clp
      dq_cl_sa = -q_clm/delt-q_clp
    endif
  else ! l_sb_all_or
    if (l_sb_dumpall) then
      ! dump all water
      !
      ! remaining water =
      !    + condesable water available
      !    - already condensed
      !    - newly condendsed
      !    - removed by mphys processed
      !

      !  calculating amount of available water
      ql_res=(qt0 - qvsl - q_clm)+delt*(qtpmcr - q_clp)

      ! limiting so it does not remove more water than in clouds
      dq_cl_sa = max((-q_clm/delt)-q_clp,ql_res/delt)

      ! adjusting number of cloud droplets
      ! - calculate min size with this amount of water
      n_bmax = (q_clm+delt*(q_clp+dq_cl_sa))/(0.1*x_cl_bmin)
      n_bmax = max(n_bmax, 0.0)

      ! of course we do not want negative values - but that is alread sorted above
      ! - remove droplets so that mean size in not less than
      dn_cl_sa = min(0.0, (n_bmax-n_clm)/delt-n_clp)

      ! limit change so not in negative numbers
      dn_cl_sa = max((-n_clm/delt)-n_clp,dn_cl_sa)
    else ! l_sb_dumpall
      !
      ! and now if we want to enforce limit on cloud droplet size
      !
      ql_res=(qt0 - qvsl - q_clm)+delt*(qtpmcr - q_clp)

      ! calculate maximal available update
      cogr_max =(1.0/delt)*(x_cogr_max*(n_clm+delt*n_clp)-q_clm)

      ! dump just what is below max size by condensation growth
      dq_cl_sa = min(cogr_max-q_clp,ql_res/delt)

      ! ie. either whole amount, or only that much that droplet size will be: xc_cogr_max
      ! other possibility: require it to be larger than 0
      ! and prevent negative values of svm + delt *svp
      dq_cl_sa = max((-q_clm/delt)-q_clp,dq_cl_sa)

      ! adjusting number of cloud droplets
      ! - calculate min size with this amount of water
      n_bmax = (q_clm+delt*(q_clp+dq_cl_sa))/(0.5*x_cl_bmin)
      n_bmax = max(n_bmax, 0.0)

      ! - remove droplets so that mean size in not by order of magnitude less than x_cl_bmin
      dn_cl_sa = min(0.0, (n_bmax-n_clm)/delt - n_clp)

      ! limit change so not in negative numbers
      dn_cl_sa = max((-n_clm/delt)-n_clp,dn_cl_sa)

      !-------------------------------------
      !! cloud water mixing ratio
      ! + mphys changes in cloud water
      ! + remaining water:
      !    + condesable water available
      !    - already condensed
      !    - newly condendsed
      !      ( {change in cloud water} = {newly condendsed or deposited}
      !                                  {removed by cloud processes}
      !      )
      !      ( - {newly condendsed or deposited} =
      !          - ({change in liquid cloud water}+{change in ice cloud water})
      !          + {removed by cloud processes}
      !      )
      ! -----------------------------
    endif ! l_sb_dumpall
  endif ! l_sb_all_or

  ! and update
  q_clp  = q_clp + dq_cl_sa
  n_clp  = n_clp + dn_cl_sa

  if (l_tendencies) then
    tend(idq_cl_sa) = dq_cl_sa
    tend(idn_cl_sa) = dn_cl_sa
  endif
  if (l_statistics) then
    statistics(imphys_ev) = statistics(imphys_ev) + min(0.0, dq_cl_sa)
    statistics(imphys_cond) = statistics(imphys_cond) + max(0.0, dq_cl_sa)
  endif
end subroutine satadj3


!    recovery of ccn
!
!   - to be later replaced based on advance literature
!   - so far just and easy recovery
!     of ccn based on number of water particles that evaporated, sublimated
!     or got removed with remaining positive n_
! -------------------------------------------------------------------------
subroutine recover_cc
  implicit none

  ! decrease in total amount of potential CCN
  n_ccp = n_ccp         &
        + dn_cl_sc      &
        + dn_cl_au      &
        + dn_cl_ac      &
        + dn_cl_hom     &
        + dn_cl_het     &
        + dn_cl_rime_ci &
        + dn_cl_rime_hs &
        + dn_cl_rime_hg
      ! NOTE: dn_cl_se moved to column processes

  ! recovery of potential CCN
  n_ccp = n_ccp + c_rec_cc*ret_cc
end subroutine recover_cc


! ***************************************************************
! Ice multiplication of Hallet and Mossop (1974)
!
! - written as described in Seifert (2002)
! - implementation similar to the one in ICON model
!
! - returns:
!    dq_i_hm  - change in mass content during the process
!    dn_i_hm  - number of newly produced ice particles
!
! conditions : t0.lt.T_3
!
! ***************************************************************
subroutine hallet_mossop3(t0,dq_rime,q_e,dq_i_hm,dn_i_hm)

  implicit none

  ! inputs
  real, intent(in)  :: t0        ! tmp0 at this gridpoint
  real, intent(in)  :: dq_rime   ! riming rate
  real, intent(in)  :: q_e       ! amount of that ice phase

  ! outputs
  real, intent(out) :: dq_i_hm   ! tendency in q_i by H-M process
  real, intent(out) :: dn_i_hm   ! tendency in n_i by H-M process

  ! constants in calculation
  real, parameter :: c_spl  = c_spl_hm74 &
                    ,c_1_hm = 1.0/(tmp_opt_hm74-tmp_min_hm74) &
                    ,c_2_hm = 1.0/(tmp_opt_hm74-tmp_max_hm74)

  ! local variables
  real :: mult_1,mult_2,mint_1,mint_2    &  ! calculation variables
         ,dn_try,dq_try,rem_cf              ! trial variables

  ! only if riming going on temperature below 0
  if (dq_rime.gt.0) then

     ! setting coefficient for reminder
     rem_cf  = (1.0-rem_q_e_hm)/delt

     ! f_spl calculation following ICON
     mult_1 = c_1_hm * (t0-tmp_min_hm74)   ! positive in the target interval
     mult_2 = c_2_hm * (t0-tmp_max_hm74)   ! positive in the target interval

     ! now for intervals
     mint_1 = max(0.0,min(1.0,mult_1))     !  0 for T<T_min, 1 for T>T_opt
     mint_2 = max(0.0,min(1.0,mult_2))     !  0 for T>T_max, 1 for T<T_opt

     ! calculating prediction for the process
     dn_try = c_spl*mint_1*mint_2*dq_rime
     dq_try = x_ci_spl*dn_try

     ! correcting
     dq_try = min(dq_try,rem_cf*q_e+dq_rime)    !  limit splintering
     dq_try = max(dq_try,0.0)

     ! prepare updates
     dq_i_hm = dq_try
     dn_i_hm = dq_try/x_ci_spl
     !
  else
     dq_i_hm = 0.0
     dn_i_hm = 0.0
  endif
end subroutine hallet_mossop3


! ***************************************************************
! Melting of ice particles
!
! - this is a inner subroutine called from a wrapper
! - written as described in Seifert&Beheng (2004)
!   - based on Pruppacher and Klett (1997)
! - implementation similar to the one in ICON model
!
! NOTE: it is assumed that tmp0.gt.T_3
!                          q_e.gt.qicemin
!
! - returns:
!    dq_me  - mass content melting tendency
!    dn_me  - number content melting tendency
!    dq_ev  - mass content evaporation tendency
!    dn_ev  - number content evaporation tendency
!
! ***************************************************************
subroutine sb_evmelt3(avent0,avent1,bvent0,bvent1,x_bmin,n_e,n_ep,n_em &
                     ,q_e,q_ep,q_em,x_e,D_e,v_e,dq_me,dn_me,dq_ev,dn_ev)

  use modglobal, only : rv,rlv,pi,cp

  implicit none

  ! inputs -------------------------------------------------------------
  real, intent(in)  :: n_e       ! number density of ice particles
  real, intent(in)  :: n_ep      !   and its tendency
  real, intent(in)  :: n_em      !   and the svm value
  real, intent(in)  :: q_e       ! mass density of ice particles
  real, intent(in)  :: q_ep      !   and its tendency
  real, intent(in)  :: q_em      !   and the svm value
  real, intent(in)  :: x_e       ! mean size of the ice particles
  real, intent(in)  :: D_e       ! mean diameter of particles
  real, intent(in)  :: v_e       ! mean terminal velocity of particles
  real, intent(in)  :: avent0    ! ventilation coefficient
  real, intent(in)  :: avent1    ! ventilation coefficient
  real, intent(in)  :: bvent0    ! ventilation coefficient
  real, intent(in)  :: bvent1    ! ventilation coefficient
  real, intent(in)  :: x_bmin    ! minimal size of the hydrometeor
  ! real, intent(in)  :: k_melt  ! depositional growth constant

  ! outputs  -----------------------------------------------------------
  real, intent(out) :: dq_me     ! melting tendency in q_e
  real, intent(out) :: dn_me     ! melting tendency in n_e
  real, intent(out) :: dq_ev     ! evaporation tendency in q_e
  real, intent(out) :: dn_ev     ! evaporation tendency in n_e

  ! local variables ----------------------------------------------------
  real, parameter    :: k_melt = 2*pi/rlme
  real, parameter    :: k_ev   = 2*pi

  real, parameter    :: dvleorv   = Dv*rlvi/rv
  real, parameter    :: ktdtodv   = Kt**2/(cp*rho0s*Dv) ! K_T * D_T / D_v and D_T = K_T/(c_p \rho_s)

  ! allocate fields and fill
  real ::   S      & ! subsaturation
           ,g_me   & ! thermodynamic term for melting
           ,g_ev   & ! thermodynamic term for evaporation
           ,x_er   & ! mean size of particles
           ,f0     & ! ventilation factor
           ,f1     & ! ventilation factor
           ,nrex   & ! reynolds number
           ,me_q   & ! basic melting rate in q
           ,me_n   & ! basic melting rate in n
           ,g_ev_const

  ! - filling
  S           =  0.0
  g_me        =  0.0
  g_ev        =  0.0
  x_er        =  0.0
  f0          =  0.0
  f1          =  0.0
  nrex        =  0.0
  me_q        =  0.0
  me_n        =  0.0

  ! constant evaporation parameter for melting particles
  ! based on: G = (rv * tmp0) / (Dv*esl) + rlv/(Kt*tmp0)*(rlv/(rv*tmp0) -1.)
  g_ev_const = (rv*T_3)/(Dv*eslt3)+rlv/(Kt*T_3)*(rlv/(rv*T_3) -1.)
  g_ev_const = 1.0 / g_ev_const

  ! set tendencies to 0
  dq_me = 0.0
  dn_me = 0.0
  dq_ev = 0.0
  dn_ev = 0.0

  ! preparing calculation
  ! calculating for all cells with the value
  ! calculation of the subsaturation
  S = min(0.0,((qt0-q_cl)/qvsl- 1.0))

  ! calculating the thermodynamic term for evaporation
  g_ev = g_ev_const ! 1.0/g_ev

  ! calculation of the thermodynamic term for melting
  g_me= - k_melt*(ktdtodv*(tmp0-T_3)+dvleorv*(esl/tmp0-eslt3/T_3))

  ! calculating real mean particle mass
  x_er = q_e/n_e

  ! calculating N_re Reynolds number
  nrex= D_e*v_e/nu_a

  ! calculating from prepared ventilation coefficients
  f0  = avent0+bvent0*Sc_num**(1.0/3.0)*nrex**0.5
  f1  = avent1+bvent1*Sc_num**(1.0/3.0)*nrex**0.5

  ! preparing updates for evaporation
  dn_ev = k_ev*g_ev*S*n_e*D_e*f0/max(x_bmin, x_er)
  dq_ev = k_ev*g_ev*S*n_e*D_e*f1

  ! preparing updates for melting
  me_n = g_me*n_e*D_e*f0/max(x_bmin, x_er)
  me_q = g_me*n_e*D_e*f1

  ! and limiting so not removing more than available
  ! basic correction of melting rate
  me_q = min(0.0,max(me_q,-q_em/delt - q_ep))

  ! basic correction for number tendency in melting
  me_n = min(0.0,max(me_n,-n_em/delt - n_ep))

  ! and prevent melting of all particles while leaving mass ?
  dq_ev = min(0.0,max(dq_ev,-q_em/delt - q_ep))

  ! basic correction for number tendency in evaporation
  dn_ev = min(0.0,max(dn_ev,-n_em/delt - n_ep))

  dq_me = me_q
  dn_me = me_n
end subroutine sb_evmelt3

end module modbulkmicro3_point
