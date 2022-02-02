!> \file modmicrodata3.f90
!!  Variables necessary for the microphysics3
!!   based on the microphysics scheme described in:
!!   \see  Seifert (2002)
!!   \see  Seifert and Beheng (Met Atm Phys, 2006)
!>
!!  Variables necessary for the sb3 microphysics
!!  \author Jan Chylik, IGMK
!!  \author Jisk Attema, NLeSC
!>
!  This file is part of DALES.
!
! DALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distrib:uted in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!

  module modmicrodata3

  implicit none
  save

  ! flags used in settings
  logical :: l_sb_classic      = .true.  &  !< .true.  ! whether to use classis S&B setting
            ,l_sb_dumpall      = .true.  &  !< whether to dump all water in cloud droplets
            ,l_sb_all_or       = .true.  &  !< all-or-nothing for number of cloud droplets in saturation adjustment
            ,l_setclouds       = .true.  &  !< whether to start by estimating clouds
            ,l_setccn          = .true.  &  !< whether to set CCN by constant concentration per M
            ,l_corr_neg_qt     = .true.  &  !< whether to adjust qt and thlp in correction for negative hydrometeor updates
            ,l_sb_dbg          = .false. &  !< .true. ! .false.
            ,l_sb_dbg_extra    = .false.     !.true.    ! .true.     !< whether to write debug messages

  ! flags indicating what was already done
  logical :: l_clouds_init     = .false. &  !< whether clouds were already initialised
            ,l_ccn_init        = .false.    !< whether CCN alredy initialised

  ! flags for different parts of microphysical processes
  logical :: l_sb_lim_aggr     = .true.  &  !< whether to start with snow aggregation only after specific size
            ,l_sb_stickyice    = .false. &  !< whether the ice is sticky or the stickyness limited as in Seifert
            ,l_sb_clnuc_first  = .true.  &  !< whether to first correct for nucleated cloud water, then other
            ,l_sb_conv_par     = .true.  &  !< whether to use partial conversion to graupel
            ,l_c_ccn           = .false. &  !< whether to use constant coefficient c_ccn
            ,l_sb_nuc_sat      = .false. &  !< cloud nucleation: whether include q_cl in super-saturation calculatio
            ,l_sb_sat_max      = .true.  &  !< whether to stop nucleation above sat_max
            ,l_sb_nuc_expl     = .true.  &  !< cloud nucleation: whether explicit nucleation
            ,l_sb_nuc_diff     = .true.  &  !< cloud nucleation: whether to use derivation
            ,l_sb_inuc_sat     = .false. &  !< ice nucleation:   whether include q_cl in super-saturation calculation
            ,l_sb_inuc_expl    = .false. &  !< ice nucleation:   whether explicit nucleation
            ,l_sb_reisner      = .true.     !< whether to use Reisner correction in ice nucleation



  logical :: l_hdump          = .false.  &  !< whether to do hydrometeor dumps
            ,l_hbinary        = .false.  &  !<- whether to use binary
            ,l_hdiracc        = .false.  &
            ,l_tendencies     = .true.   &  !<- if to write full tendencies TODO: check if there is an existing flag
            ,l_statistics     = .true.      !<- if to write statistics

  real ::  Nc0             = 70.0e6   &  !<- proposed number of droplet in namelist
             ,xc0_min         = 4.2e-15  &  !<- xcmin  min mean mass of cloud water
             ,Nccn0           = 100.0e6     !<- proposed initial number of cc

  ! addjusting position of species to scalars fields
  ! need to be backward compatible to bulkmicro
  !   (i.e. in_hr = inr, iq_hr = iqr
  !    integer :: inr = 1, iqr=2 )
  integer, parameter :: ncols = 12 ! NOTE: - should be even for the untranspose_svs()
                                   !       - keep the n/q fields next to eachother (performance)
                                   !       - ncols /= nsv, as chemistry could add more scalars
  integer ::  in_hr           =  1       &  ! |<  in_hr : number content [     kg^{-1}] for rain,
             ,iq_hr           =  2       &  ! |<  iq_hr : water content  [ kg  kg^{-1}] for rain,
             ,in_cl           =  3       &  ! |<  in_cl : number content [     kg^{-1}] for cloud droplets,
             ,iq_cl           =  4       &  ! |<  iq_cl : water content  [ kg  kg^{-1}] for cloud droplets,
             ,in_ci           =  5       &  ! |<  in_ci : number content [     kg^{-1}] of ice crystals,
             ,iq_ci           =  6       &  ! |<  iq_ci : mass  content  [ kg  kg^{-1}] of ice crystals,
             ,in_hs           =  7       &  ! |<  in_hs : number content [     kg^{-1}] of snow,
             ,iq_hs           =  8       &  ! |<  iq_hs : water content  [ kg  kg^{-1}] of snow,
             ,in_hg           =  9       &  ! |<  in_hg : number content [     kg^{-1}] of graupels,
             ,iq_hg           = 10       &  ! |<  iq_hg : water content  [ kg  kg^{-1}] of graupels,
             ,in_cc           = 11       &  ! |<  in_cc : number content [     kg^{-1}] of cloud condensation nuclei,
             ,in_tr1          = 12          ! |<  in_tr : number content [ kg^{-1}] kg^{-1}] of a passive tracer,

  ! adding values for setting of microphysiscs
  ! thresholds
  real, parameter ::  q_hr_min     = 1.0e-14    & !<  Rain  specific mixing ratio treshold min
                     ,q_hs_min     = 1.0e-14    & !<  Snow  specific mixing ratio treshold min
                     ,q_hg_min     = 1.0e-14    & !<  Graupel  specific mixing ratio treshold min
                     ,qicemin      = 1.0e-14    & !< ice mixing ratio treshold -have  to be GE q_hs_min, q_hg_min
                     ,qcliqmin     = 1.0e-14    & !< liquid cloud water mixing ratio treshold
                     ,qsnowmin     = 1.0e-14    & !< Snow mixing ratio treshold
                     ,qgrmin       = 1.0e-14    & !< Graupel mixing ratio treshold
                     ,n_c_min      = 1.0e0      & !< cloud droplet number content [ kg^{-1}] treshold
                     ,n_h_min      = 1.0e-4     & !< hydrometeor number content [ kg^{-1}] treshold
                     ,ssice_min    = 1.0e-15    & !< min supersaturation for ice nucleation
                     ,dvrlim       = 0.35e-3    & !< raindrop size threshold for collision breakup
                     ,dvrbiglim    = 0.9e-3     & !< raindrop size threshold for collision breakup
                     ,cc_min_ratio = 0.2        & !< ratio how many ccn has to remain after one step of nucleation
                     ,xc_bmax      = 1.01e-11   & !< max mean size of the cloud droplet
                     ,x_cl_bmax    = xc_bmax    & !< max mean size of cloud droplet
                     ,x_ci_bmax    = 1.0e-7     & ! orig: 1.0e-7 !< max mean size of cloud ice
                     ,x_hr_bmax    = 5.0e-6     & !< max mean size of rain drop
                     ,x_hs_bmax    = 1.0e-6     & ! orig: 1.0e-7 ! !< max mean size of snow
                     ,x_hg_bmax    = 1.0e-4     & !< max mean size of graupel - originally: 1.0e-4
                     ,x_cl_bmin    = 4.20e-15   & !< min mean size of cloud droplet
                     ,x_ci_bmin    = 1.0e-12    & !< min mean size of cloud ice
                     ,x_hr_bmin    = 1.0e-10    & !< min mean size of rain drop
                     ,x_hs_bmin    = 1.73e-9    & !< min mean size of snow
                     ,x_hg_bmin    = 2.6e-10    & !< min mean size of graupel - originally: 1.0e-4
                     ,x_degr_max   = x_ci_bmax  & !< max size of depositional growth
                     ,x_cogr_max   = xc_bmax    & !< max size by condensation growth
                     ,N_0min       = 2.5e5      & !< N [m-3] lower limit in rain integral
                     ,N_0max       = 2.0e7      & !< N [m-3] upper limit in rain integral
                     ,lbdr_min     = 1.0e3      & !< slope lower limit in rain integral
                     ,lbdr_max     = 1.0e4      & !< slope upper limit in rain integral
                     ,ssice_lim    = 0.1        & !< max 10% supersaturation with respect to ice
                     ,split_factor = 1.5        & !< factor in time splitting procedure for precipitation
                     ,rem_n_cl_min = 0.1        & !< ratio of the remaing cloud droplet after sc and br
                     ,rem_n_ci_min = 0.2        & !< ratio of the remaing ice cloud particles after sc
                     ,rem_n_hr_min = 0.1        & !< ratio of the remaing raindrop number after sc and br
                     ,rem_n_hs_min = 0.3        & !< ratio of the minimal snow number after the  self-collection to the original one
                     ,rem_n_min_cv = 0.3       !< ratio of the minimal number remaining after conversion
                     ! ->> continue from here

  ! settings of physics
  real, parameter ::  rlvi       = 2.834e6      & !< latent heat for sublimation ice (value ams glossary)
                     ,rlme       = 3.337e5      & !< latent heat for melting (value ams glossary)
                     ,c_water    = 4.186e3      & !< specific heat of water vapour
                     ,rho0s      = 1.225        & !< reference density
                     ,k1nuc      = 1.58         & !< 1st constant in nuc calculation (G09a)
                     ,k2nuc      = 0.72           !< 2nd constant in nuc calculation (G09a)

  ! settings of physics
  real, parameter ::  k_cr       = 5.25         & !< kernel for accretion
                     ,k_cc       = 4.44e9       & !< kernel cloud-cloud
                     ,kappa_br   = 2.3e3          !< kappa kernel for collision brakeup

  ! parameters in liquid droplet nucleation scheme
  real, parameter ::  x_inuc_SB   = 1.0e-12     & !< mass of nucleated ice particle
                     ,x_cnuc_SB   = 1.0e-12     & !< mass of nucleated liquid droplet
                     ,sat_min_SB  = 1.0e-5      & ! 0.01 ! < min supersaturation for nucleation [%]
                     ,sat_max_SB  = 1.1         & ! 1.1  ! < max supersaturation for nucleation [%]
                     ,n_clmax_marit = 1.5e8     & ! maximum number concentration in case of constant c_ccn
                     ,n_clmax_conti = 1.5e10    & ! maximum number concentration in case of constant c_ccn
                     ,c_ccn_marit = 1.0e8       & !< C_CCN parameter for maritime conditions
                     ,c_ccn_conti = 1.26e9      & !< C_CCN parameter for conitnental conditions
                     ,kappa_marit = 0.462       & !< kappa_ccn parameter for maritime conditions
                     ,kappa_conti = 0.308         !< kappa_ccn parameter for conitnental conditions


  ! additional of standard physical parameters
  real, parameter ::  T_3     = 273.15     & !< T_3 temperature of freezing point for water
                     ,rhoeps  = 900.0        !< turbulent ice density [kg m^{-3}]

  ! setting for statistics
  real, parameter ::  eps_hprec = 1.0e-8     !< threshold for precipitation :,epsqr = 1.0e-8

  ! parameters for cloud particles and hydrometeors
  !   taken directly from S&B, Table 1
  real, parameter ::  a_cl    = 0.124      & !< a param for: cloud
                     ,a_hr    = 0.124      & !<              raindrop
                     ,a_ci    = 0.217      & !<              cloud ice
                     ,a_hs    = 8.156      & !<              snow
                     ,a_hg    = 0.190      & !<              graupel
                     ,b_cl    = 0.33333    & !< b param for: cloud
                     ,b_hr    = 0.33333    & !<              raindrop
                     ,b_ci    = 0.302115   & !<  o: 0.302    cloud ice
                     ,b_hs    = 0.526      & !<              snow
                     ,b_hg    = 0.323      & !<              graupel
                     ,c_cl    = 2.0        & !< capacity param for: cloud
                     ,c_hr    = 2.0        & !<              raindrop
                     ,c_ci    = 3.14159    & !<              cloud ice
                     ,c_hs    = 2.0        & !<              snow
                     ,c_hg    = 2.0        & !<              graupel
                     ,al_cl   = 3.75e5     & !< alpha param for: cloud
                     ,al_hr   = 159.0      & !<              raindrop
                     ,al_ci   = 41.9       & !<  o: 317.0   cloud ice
                     ,al_hs   = 27.7       & !<              snow
                     ,al_hg   = 40.0       & !<              graupel
                     ,be_cl   = 0.66667    & !< beta param for: cloud
                     ,be_hr   = 0.266      & !<              raindrop
                     ,be_ci   = 0.36       & !<  0.363       cloud ice
                     ,be_hs   = 0.216      & !<              snow
                     ,be_hg   = 0.230      & !<              graupel
                     ,ga_cl   = 1.0        & !< gamma param for: cloud
                     ,ga_ci   = 0.5        & !<                  ice
                     ,ga_hr   = 0.5        & !<                  rain
                     ,ga_hs   = 0.5        & !<                  snow
                     ,ga_hg   = 0.5        & !<                  graupel
                     ,nu_cl_cst   = 1.0        & !< nu param for: cloud
                     ,nu_hr_cst   = -0.66667   & !<              raindrop
                     ,nu_ci_cst   = 0.0        & !<    1.0       cloud ice
                     ,nu_hs_cst   = 1.0        & !<              snow
                     ,nu_hg_cst   = 1.0        & !<              graupel
                     ,mu_cl_cst   = 1.0        & !< mu p1.0aram for: cloud
                     ,mu_hr_cst   = 0.33333    & !<              raindrop
                     ,mu_ci_cst   = 0.33333    & !<              cloud ice
                     ,mu_hs_cst   = 0.33333    & !<              snow
                     ,mu_hg_cst   = 0.33333    & !<              graupel
                     ,d_wfallmax_hr = 9.9      & !< default max terminal w for: rain
                     ,d_wfallmax_hg = 11.9     & !< default max terminal w for: graupel
                     ,a_tvsbc     = 9.65       & !<  coeff in terminal velocity param
                     ,b_tvsbc     = 10.3       & !<  coeff in terminal velocity param
                     ,c_tvsbc     = 600.0        !<  coeff in terminal velocity param

  ! rain scheme
  real, parameter ::  mur0_G09b   = 30.0       & !< N in G09b rain scheme
                     ,c_G09b      = 0.008      & !< N in G09b rain scheme
                     ,exp_G09b    = 0.6            !< N in G09b rain scheme

  ! ice scheme
  real, parameter ::  N_M92       = 1.0e3      & !< N in M92 ice nucleation scheme
                     ,a_M92       = -0.639     & !< a in M92 ice nucleation scheme
                     ,b_M92       = 12.96      & !< b in M92 ice nucleation scheme
                     ,tmp_inuc_sb = 268.15     & !< temperature limit below shich is ice nucleation possible
                     ,n_i_max_R019 = 1100.0    & !< basic limiting number of primary ice particles (de Roode 2019
                     ,N_R98       = 0.01       & !< number parameter in Reisner correction R98
                     ,c_R98       = 246.15     & !< temperature parameter in Reisner correction R98
                     ,b_R98       = 0.6        & !< coefficient in Reisner correction
                     ,a1_R98      = 0.1        & !< by order of magnitude below
                     ,a2_R98      = 100.0        !< by order or two orders of magnitude above



  ! parameters for homogeneous freezing - based on Cotton & Field (2002)
  real, parameter :: C_CF02        = -7.36    & !< additive constant in exponential term in CF02
                    ,B_CF02        = -2.996   & !< multiplicative constant in exponential term in CF02
                    ,CC_CF02       = -243.15  & !< additive constant next to T in CF02
                    ,tmp_lim1_CF02 = 243.15   & !< Te boundary for slow freezing
                    ,tmp_lim2_CF02 = 208.15   & !< Te boundary for quick freezing
                    ,offset_CF02   = 273.15   & !< Te offset for calculating in mid freezing
                    ,C_20_CF02     = -243.4   & !< constatnts in mid range
                    ,B_21_CF02     = -14.75   & !< constatnts in mid range
                    ,B_22_CF02     = -0.307   & !< constatnts in mid range
                    ,B_23_CF02     = -0.00287 & !< constatnts in mid range
                    ,B_24_CF02     = -102e-7  & !< constatnts in mid range
                    ,C_30_CF02     = 25.63      !< constant in exponential term in CF02 for very low temp

  ! parameters for heterogeneous freezing and riming
  real, parameter :: A_het      = 0.2         & !< mult const. in heterogeneous freezing
                    ,B_het      = 0.65        & !< mult const. in exp in heterogeneous freezing
                    ,al_0snow   = 0.01        & !< parameter alpha_0,snow in partial conv. s-->g
                    ,al_0ice    = 0.68        & !< parameter alpha_0,snow in partial conv. i-->g
                    ,D_convmin  = 5.0e-4        !< 500 \mi m - min size of ice/snow for partial conv.
                    !-> close it later here

  ! parameters for ventilation coefficients
   real, parameter :: a_v_i      = 0.86       & ! a_v for ice
                     ,a_v_r      = 0.78       & ! a_v for rain
                     ,a_v_s      = 0.78       & ! a_v for snow
                     ,a_v_g      = 0.78       & ! a_v for graupel
                     ,b_v_i      = 0.208      & ! b_v for ice
                     ,b_v_r      = 0.308      & ! b_v for rain
                     ,b_v_s      = 0.308      & ! b_v for snow
                     ,b_v_g      = 0.308        ! b_v for graupel


  ! parameters for ice multiplication
   real, parameter :: c_spl_hm74 = 3.5e8      & !< c_{spl} in Hallet-Mossop process (1974)
                     ,tmp_min_hm74 = 265.0    & !< T_{spl,min} in Hallet-Mossop process (1974)
                     ,tmp_opt_hm74 = 268.0    & !< T_{spl,opt} in Hallet-Mossop process (1974)
                     ,tmp_max_hm74 = 270.0    & !< T_{spl,max} in Hallet-Mossop process (1974)
                     ,x_ci_spl   =x_ci_bmin   & !< size of split ice particles
                     ,rem_q_e_hm   = 0.5        !< how much of the mass of the particle should remain after splintering


  ! parameters for collisions
  real, parameter ::  D_c_a      = 15.0e-6    & !< \bar{D}_c,a : droplet min size for collision effic.
                     ,D_c_b      = 40.0e-6    & !< \bar{D}_c,b : droplet size for high coll. effic.
                     ,D_i_a      = 75.0e-6    & !< \bar{D}_i,a : minimal ice size for coll. effic. - taken from ICON, 2017
                     ,D_i_b      = 398.0e-6   & !< \bar{D}_i,b : ice size for high coll. effic
                     ,D_i0       = 150.0e-6   & !< \bar{D}_i,0 : minimal ice size for coll. effic.
                     ,D_s0       = 150.0e-6   & !< \bar{D}_s,0 : minimal snow size for coll. effic.
                     ,D_g0       = 150.0e-6   & !< \bar{D}_g,0 : minimal graupel size for coll. effic.
                     ,E_i_m      = 0.8        & !< \bar{E}_i(D)  : mean efficiency for ice coll. liq
                     ,E_s_m      = 0.8        & !< \bar{E}_s(D)  : mean efficiency for snow coll. liq
                     ,E_g_m      = 1.0        & !< \bar{E}_g(D)  : mean efficiency for graupel coll. liq
                     ,E_er_m     = 1.0        & !< \bar{E}_er(D) : mean collisions eff. for ice particle - rain
                     ,E_ii_m     = 0.1        & ! 1.0 !< \bar{E}_ee(D) : mean collisions eff. for cloud ice - ice
                     ,E_is_m     = 1.0        & !< \bar{E}_ee(D) : mean collisions eff. for cloud ice - snow
                     ,E_ee_m     = 1.0        & !< \bar{E}_ee(D) : mean collisions eff. for other ice - ice particle
                     ,E_gg_s     = 0.0        & !< \bar{E}_st,gg(D) : mean stick eff. graupel-graupel
                     ,E_gi_s     = 0.0        & !< \bar{E}_st,gi(D) : mean stick eff. graupel-ice
                     ,c_E_o_s    = 1.0        & !< : sticking efficiency for other
                     ,stick_off  = -273.15    & !< temperature offset in calculating sticking eff.
                     ,B_stick    = 0.09       & !< multiplicative constant is sticking efficiency
                     ,B_stick_ii = 0.08059    & !< coefficient (taken from Seifert thesis) 0.035*log(10)
                     ,C_stick_ii = -0.7       & !< additive coefficien in exponent ( -||- )
                     ,E_ii_maxst = 0.2        & !< maximal sticking efficiency for ice-ice
                     ,E_ss_maxst = 0.1        & !< maximal sticking efficiency for snow-snow
                     ,sigma_cl   = 0.0        & !< velocity variance for cloud droplets
                     ,sigma_hr   = 0.0        & !< velocity variance for raindrops
                     ,sigma_ci   = 0.2        & !< velocity variance for cloud ice
                     ,sigma_hs   = 0.2        & !< velocity variance for snow
                     ,sigma_hg   = 0.0        & !< velocity variance for graupels
                     ,D_mincv_ci = 500.0e-6   & !< min size for rime conversion
                     ,D_mincv_hs = 500.0e-6   & !< min size for rime conversion
                     ,x_hs_cvmin = 0.1e-9     & !< sep. size for snow conversion to graupel - from ICON, 2017
                     ,x_ci_cvmin = 0.1e-9     & !< sep. size for ice conversion to graupel
                     ,D_crit_ii  = 100.0e-6   & !< critical size for start of ice aggregation
                     ,q_crit_ii  = 1.0e-6     & !< critical mass for start of ice aggregation
                     ,rime_min   = 1.0e-18    & !< minimal riming rate when considering enhanced melting
                     ,c_ccn_ev_c = 1.0        & !< how often CCN recovered from evaporating cloud droplet
                     ,c_ccn_ev_r = 0.0        & !< how often CCN recovered from evaporating rain drop
                     ,c_rec_cc   = 1.0 ! 1.0    !< ccn recovery - how many evaporated droplets can serve as ccn
                     ! ->later should also include different sticking efficiencies

    ! and now the default value of parameters that can be adjusted  by namelist
    real ::  x_cnuc          = x_cnuc_SB      & !< mass of nucleated liquid droplet
            ,c_ccn           = c_ccn_marit    & !< C_CCN parameter for maritime conditions
            ,n_clmax         = n_clmax_marit  & !< maximum number concentration in case of constant C_CCN
            ,kappa_ccn       = kappa_marit    & !< kappa_ccn parameter for maritime conditions
            ,sat_min         = sat_min_SB     & ! < min supersaturation [%] for nucleation
            ,sat_max         = sat_max_SB     & ! < max supersaturation [%] for nucleation
            ,x_inuc          = x_inuc_SB      & !< mass of nucleated ice particle
            ,N_inuc          = N_M92          & !< N in M92 ice nucleation scheme
            ,a_inuc          = a_M92          & !< a in M92 ice nucleation scheme
            ,b_inuc          = b_M92          & !< b in M92 ice nucleation scheme
            ,tmp_inuc        = tmp_inuc_sb    & !< temperature limit below shich is ice nucleation possible
            ,n_i_max         = n_i_max_R019   & !< basic limiting number of primary ice particles (de Roode 2019
            ,N_inuc_R        = N_R98          & !< number parameter in Reisner correction R98
            ,c_inuc_R        = c_R98          & !< temperature parameter in Reisner correction R98
            ,b_inuc_R        = b_R98          & !< coefficient in Fletcher formula
            ,a1_inuc_R       = a1_R98         & !< by order of magnitude below
            ,a2_inuc_R       = a2_R98           !< by order of magnitude above


  ! parameters for statistics
  real, parameter ::    q_cl_statmin    = 1.0e-5         & !< minimal cloud liqud water content to be counted as cloud
                       ,q_ci_statmin    = 1.0e-5         & !< minimal cloud ice water content to be counted as cloud
                       ,def_hdump_dtav  = 900.             !< default for time outputs

  ! parameters for field outputs
  ! real            ::  cdump_mul          = 1.0e5    !< a factor by which are hydrometeor specie multiplied for recording
  ! integer         ::  cdump_byte         = 4        !< how many byte integer used

  integer, parameter :: nmphys        = 8       ! rounded to power of 2 for alignment (performance)
  integer, parameter :: imphys_freeze = 1     & ! - change in th due to freezing
                       ,imphys_melt   = 2     & ! - change in th due to melting
                       ,imphys_cond   = 3     & ! - change in th due to condensation
                       ,imphys_ev     = 4     & ! - change in th due to evaporation
                       ,imphys_dep    = 5     & ! - change in th due to deposition
                       ,imphys_sub    = 6       ! - change in th due to sublimation

  ! adding other process contributions
  integer, parameter :: ntends = 93 ! The number of tendencies below
  integer, parameter ::  &
     idn_cl_nu       =  1 &      !< droplet nucleation rate
    ,idn_ci_inu      =  2 &      !< ice nucleation rate
    ,idn_cl_au       =  3 &      !< change in number of cloud droplets due to autoconversion
    ,idq_hr_au       =  4 &      !< change in mass of raindrops due to autoconversion
    ,idn_hr_au       =  5 &      !< change in number of raindrops due to autoconversion
    ,idq_hr_ac       =  6 &      !< change in mass of raindrops due to accretion
    ,idn_cl_ac       =  7 &      !< change in number of cloud droplets due to accretion
    ,idn_hr_br       =  8 &      !< change in number of raindrops due to breakup
    ,idn_hr_sc       =  9 &      !< change in number of raindrops due to self-collection
    ,idq_hr_ev       = 10 &      !< change in mass of raindrops due to evaporation
    ,idn_hr_ev       = 11 &      !< change in number of raindrops due to evaporation
    ,idq_ci_dep      = 12 &      !< deposition rate for clouds
    ,idq_hs_dep      = 13 &      !< deposition rate for snow
    ,idq_hg_dep      = 14 &      !< deposition rate for graupel
    ,idq_ci_rime     = 15 &      !< riming growth of ice
    ,idn_cl_rime_ci  = 16 &      !<  - and impact on n_cl
    ,idq_hs_rime     = 17 &      !< riming growth of snow
    ,idn_cl_rime_hs  = 18 &      !<  - and impact on n_cl
    ,idq_hg_rime     = 19 &      !< riming growth for graupel
    ,idn_cl_rime_hg  = 20 &      !<  - and impact on n_cl
    ,idq_hshr_rime   = 21 &      !< riming growth for snow with rain #remove
    ,idn_hr_rime_hs  = 22 &      !<  - and impact on n_hr            # remove
    ,idq_hghr_rime   = 23 &      !< riming growth for graupel with rain
    ,idn_hr_rime_hg  = 24 &      !<  - and impact on n_hr
    ,idq_hr_rime_ri  = 25 &      !< rain loss from riming of ice+rain->gr #remove
    ,idq_ci_rime_ri  = 26 &      !< ice loss from riming of ice+rain->gr #remove
    ,idn_ci_rime_ri  = 27 &      !< ice number loss from riming of ice+rain->gr #remove
    ,idn_hr_rime_ri  = 28 &      !< rain number loss from riming of ice+rain->gr #remove
    ,idq_hr_col_rs   = 29 &      !< rain loss from riming of ice+snow->gr
    ,idq_hs_col_rs   = 30 &      !< rain number loss from riming of ice+snow->gr
    ,idn_hr_col_rs   = 31 &      !< snow loss from riming of ice+snow->gr
    ,idn_hs_col_rs   = 32 &      !< snow number loss from riming of ice+snow->gr
    ,idq_hr_col_ri   = 33 &      !< rain loss from riming of ice+rain->gr
    ,idq_ci_col_ri   = 34 &      !< ice loss from riming of ice+rain->gr
    ,idn_ci_col_ri   = 35 &      !< ice number loss from riming of ice+rain->gr
    ,idn_hr_col_ri   = 36 &      !< rain number loss from riming of ice+rain->gr
    ,idq_cl_het      = 37 &      !< heterogeneou freezing of cloud water
    ,idn_cl_het      = 38 &      !< heterogeneou freezing of cloud water
    ,idq_hr_het      = 39 &      !< heterogeneou freezing of raindrops
    ,idn_hr_het      = 40 &      !< heterogeneou freezing of raindrops
    ,idq_cl_hom      = 41 &      !< homogeneous freezing of cloud water
    ,idn_cl_hom      = 42 &      !< homogeneous freezing of cloud water
    ,idq_ci_col_iis  = 43 &      !< self-collection of cloud ice
    ,idn_ci_col_iis  = 44 &      !< self-collection of cloud ice
    ,idn_hs_col_sss  = 45 &      !< self-collection of snow
    ,idq_hsci_col    = 46 &      !< collection s+i - trend in q_hs
    ,idn_ci_col_hs   = 47 &      !< collection s+i - trend in n_ci
    ,idq_hghs_col    = 48 &      !< collection g+s - trend in q_hg
    ,idn_hs_col_hg   = 49 &      !< collection g+s - trend in n_hs
    ,idq_ci_cv       = 50 &      !< partial conversion ice -> graupel
    ,idn_ci_cv       = 51 &      !< partial conversion ice -> graupel
    ,idq_hs_cv       = 52 &      !< partial conversion snow-> graupel
    ,idn_hs_cv       = 53 &      !< partial conversion snow-> graupel
    ,idn_cl_sc       = 54 &      !< cloud self-collection
    ,idn_ci_mul      = 55 &      !< ice multiplication
    ,idq_ci_mul      = 56 &      !< ice multiplication
    ,idn_ci_me       = 57 &      !< number tendency melting of cloud ice
    ,idq_ci_me       = 58 &      !< mass tendency melting of cloud ice
    ,idn_hs_me       = 59 &      !< number tendency melting of snow
    ,idq_hs_me       = 60 &      !< mass tendency melting of snow
    ,idn_hg_me       = 61 &      !< number tendency melting of graupel
    ,idq_hg_me       = 62 &      !< mass tendency melting of graupel
    ,idn_ci_ev       = 63 &      !< number tendency evaporation of cloud ice
    ,idq_ci_ev       = 64 &      !< mass tendency evaporation of cloud ice
    ,idn_hs_ev       = 65 &      !< number tendency evaporation of snow
    ,idq_hs_ev       = 66 &      !< mass tendency evaporation of snow
    ,idn_hg_ev       = 67 &      !< number tendency evaporation of graupel
    ,idq_hg_ev       = 68 &      !< mass tendency evaporation of graupel
    ,idn_ci_eme_ic   = 69 &      !< number tendency enhanced melting of cloud ice by cloud water
    ,idq_ci_eme_ic   = 70 &      !< mass tendency enhanced melting of cloud ice by cloud water
    ,idn_ci_eme_ri   = 71 &      !< number tendency enhanced melting of cloud ice by rain
    ,idq_ci_eme_ri   = 72 &      !< mass tendency enhanced melting of cloud ice  by rain
    ,idn_hs_eme_sc   = 73 &      !< number tendency enhanced melting of snow by cloud water
    ,idq_hs_eme_sc   = 74 &      !< mass tendency enhanced melting of snow by cloud water
    ,idn_hs_eme_rs   = 75 &      !< number tendency enhanced melting of snow by rain
    ,idq_hs_eme_rs   = 76 &      !< mass tendency enhanced melting of snow by rain
    ,idn_hg_eme_gc   = 77 &      !< number tendency enhanced melting of graupel by liquid clouds
    ,idq_hg_eme_gc   = 78 &      !< mass tendency enhanced melting of graupel by liquid clouds
    ,idn_hg_eme_gr   = 79 &      !< number tendency enhanced melting of graupel by rain
    ,idq_hg_eme_gr   = 80 &      !< mass tendency enhanced melting of graupel by rain
    ,idn_cl_se       = 81 &      !< sedimentation for clouds water - number
    ,idq_cl_se       = 82 &      !<      -||-- mixing ration
    ,idn_ci_se       = 83 &      !< sedimentation for cloud ice - number
    ,idq_ci_se       = 84 &      !<       -||-- mixing ration
    ,idn_hr_se       = 85 &      !< sedimentation for rain - number
    ,idq_hr_se       = 86 &      !<       -||-- mixing ration
    ,idn_hs_se       = 87 &      !< sedimentation for snow - number
    ,idq_hs_se       = 88 &      !<       -||-- mixing ration
    ,idn_hg_se       = 89 &      !< sedimentation for graupel - number
    ,idq_hg_se       = 90 &      !<       -||-- mixing ration
    ,idq_cl_sa       = 91 &      !< saturation adjustment
    ,idn_cl_sa       = 92 &      !< change in n_cl due to saturation adjustment
    ,iret_cc         = 93        !< recovery of ccn

  ! ventilation parameters
  real :: aven_0r,aven_0i,aven_0s,aven_0g                 &
         ,aven_1r,aven_1i,aven_1s,aven_1g                 &
         ,bven_0r,bven_0i,bven_0s,bven_0g                 &
         ,bven_1r,bven_1i,bven_1s,bven_1g

  ! constants in calculation of moments and wall velocities
  real :: c_mmt_1cl, c_mmt_1hr, c_mmt_2cl,c_mmt_2hr       &
         ,c_v_c0, c_v_i0, c_v_s0, c_v_g0                  &
         ,c_v_c1, c_v_i1, c_v_s1, c_v_g1

  ! collision coefficients
  real ::  dlt_c0, dlt_r0, dlt_i0, dlt_s0,dlt_g0          &
          ,dlt_i0c, dlt_i0r, dlt_i0i                      &
          ,dlt_r0i, dlt_r0s, dlt_r0g                      &
          ,dlt_s0c, dlt_s0r, dlt_s0i, dlt_s0s, dlt_s0g    &
          ,dlt_g0c, dlt_g0r, dlt_g0i, dlt_g0s, dlt_g0g    &
          ,th_c0, th_r0, th_i0, th_s0,th_g0               &
          ,th_i0c, th_i0r, th_i0i                         &
          ,th_r0i, th_r0s, th_r0g                         &
          ,th_s0c, th_s0r, th_s0i, th_s0s                 &
          ,th_g0c, th_g0r, th_g0i, th_g0s, th_g0g         &
          ,dlt_c1, dlt_r1, dlt_i1, dlt_s1,dlt_g1          &
          ,dlt_i1c, dlt_i1r, dlt_i1i                      &
          ,dlt_r1i, dlt_r1s, dlt_r1g                      &
          ,dlt_s1c, dlt_s1r, dlt_s1i, dlt_s1s             &
          ,dlt_g1c, dlt_g1r, dlt_g1i, dlt_g1s, dlt_g1g    &
          ,th_c1, th_r1, th_i1, th_s1,th_g1               &
          ,th_i1c, th_i1r, th_i1i                         &
          ,th_r1i, th_r1s, th_r1g                         &
          ,th_s1c, th_s1r, th_s1i, th_s1s                 &
          ,th_g1c, th_g1r, th_g1i, th_g1s, th_g1g

  real :: wfallmax_hr     & !< max sedim w for: rain
         ,wfallmax_ci     & !<              cloud ice
         ,wfallmax_cl     & !<              cloud water
         ,wfallmax_hs     & !<              snow
         ,wfallmax_hg     & !<              graupel
         ,c_lbdr            !< coefficient in lbdr calculation l_sb_classic

  ! adding variables for calculations
  real:: k_mphys   & !< generic microphysics constant
        ,eslt3       !< saturated vapour pressure over water at T_3

  ! pre-processed precipitation fields to output in bulkmicrostat3
  real, allocatable ::  precep_l(:,:) & !< liquid surface precipitation (precep_hr)
                       ,precep_i(:,:)   !< frozen surface precipitation (precep_[ci+hs+hg])

  logical :: q_hr_mask,q_hs_mask,q_hg_mask,q_cl_mask,q_ci_mask

  real, allocatable :: tend_fsum(:,:)           &  ! Sum of individual tendencies
                      ,statistic_mphys(:,:)     &  ! Full sum of selected tendencies
                      ,statistic_sv0_count(:,:) &  ! Count of sv0 > threshold
                      ,statistic_sv0_fsum(:,:)  &  ! Full sum of sv0
                      ,statistic_sv0_csum(:,:)  &  ! Conditional sum: sv0 > threshold
                      ,statistic_svp_fsum(:,:)  &  ! Full sum of svp
                      ,statistic_svp_csum(:,:)     ! Conditional sum: sv0 > threshold

  ! indices in the tranposed array
  integer, parameter :: nprgs     = 8    ! rounded to a power of 2, to help alignment (performance)
  integer, parameter :: n_tmp0    = 1 &
                       ,n_qt0     = 2 &
                       ,n_ql0     = 3 &
                       ,n_esl     = 4 &
                       ,n_qvsl    = 5 &
                       ,n_qvsi    = 6 &
                       ,n_w0      = 7

  end module modmicrodata3
