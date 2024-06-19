!###############################################################################
!
!                       Copyright by
!   National Institute of Public Health and Environment
!                      The Netherlands
!   No part of this software may be used, copied or distributed 
!           without permission of RIVM (2015)
!
!  MODULE             : LE_DryDepos_Gas_DEPAC
!  INTERFACE          : DryDepos_Gas_DEPAC
!  AUTHOR             : Addo van Pul, Jan Willem Erisman, Ferd Sauter, Margreet van Zanten, Roy Wichink Kruit
!  FIRM/INSTITUTE     : RIVM
!  LANGUAGE           : FORTRAN-90
!  DESCRIPTION        : In this subroutine the canopy or surface resistance Rc
!                       is parameterised.
!
! Documentation by Ferd Sauter, Mar 2009.
! Deposition fluxes are computed using one of the following resistance approaches;
! Note that, with the appopriate definitions (see below), B and C are totally equivalent.
!
!           A: classical approach        !     B: with compensation points       ! C: with total            ! D: approximation of C
!                                        !                                       !    compensation point    !    with classical approach
!   zr ---------------------------- Catm !  zr ---------------------------- Catm !  zr --------- Catm       !  zr --------- Catm    
!                 |         |            !                |                      !         |                !         |       
!                 Ra        |            !                Ra                     !         Ra               !         Ra      
!                 |         | F          !                |                      !         |                !         |       
!                 Rb        |            !                Rb                     !         Rb               !         Rb      
!                 |         V            !                |                      !         |                !         |       
!   z0 ---------------------------- Cc   !  z0 ---------------------------- Cc   !  z0 --------- Cc         !  z0 --------- Cc
!       |      |      |     |            !     |        |           |            !         |                !         |              
!       |     Rinc    |     |            !     |       Rinc         |            !        Rc_tot            !        Rc_eff          
!       |      |      |     | F          !     |        |           |            !         |                !         |              
!       Rw   Rsoil   Rstom  |            !    Rw       Rsoil       Rstom         !     --------- Ccomp_tot  !     --------- C=0
!       |      |      |     |            !     |        |           |            !                          !
!       |      |      |     V            !--- Cw ---- Csoil ------ Cstom         !                          !
!   zs ---------------------------- C=0  !
!
!   zr       : reference height (m)
!   z0       : roughness length (m)
!   zs       : surface height (m)
!   Catm     : atmospheric concentration (ug/m3)
!   Cc       : concentration at canopy top (ug/m3)
!   Cw       : external leaf surface water concentration (compensation point) (ug/m3)
!   Csoil    : soil concentration (compensation point) (ug/m3)
!   Cstom    : stomatal concentration (compensation point) (ug/m3)
!   Ccomp_tot: total compensation point (weighed average of separate compensation points) (ug/m3)
!   F        : flux (ug/(m2 s))
!   Ra       : aerodynamic resistance (s/m)
!   Rb       : boundary layer resistance (s/m)
!   Rw       : external leaf resistance (s/m)
!   Rinc     : in canopy resistance (s/m)
!   Rsoil    : soil resistance (s/m)
!   Rsoil_eff: effective soil resistance, Rsoil_eff = Rinc + Rsoil (s/m)
!   Rstom    : stomatal resistance (s/m)
!   Rc_tot   : total canopy or surface resistance (s/m) 1/Rc_tot = 1/Rw + 1/Rsoil_eff + 1/Rstom
!   Rc_eff   : effective total canopy or surface resistance (s/m) rc_eff = ((ra + rb)*ccomp_tot + rc_tot*catm)/(catm-ccomp_tot)
!   vd       : deposition velocity (m/s)
!   ve       : exchange velocity (deposition or emission) (m/s)
!   H        : layer height of numerical model
!
!    A. classical approach             !  B and C. compensation points
!                                      !
!  1     1         1          1        !  1     1         1          1
! -- = ---- + ----------- + -----      ! -- = ---- + ----------- + -----
! Rc     Rw    Rsoil_eff    Rstom      ! Rc    Rw    Rsoil_eff    Rstom
!                                      !
! deposition velocity:                 ! exchange velocity (deposition or emission):
! vd = 1/(Ra + Rb + Rc)                ! ve = 1/(Ra + Rb + Rc)
!                                      !
!                                      ! Separate fluxes over external leaf, soil and stomata (B) equal total flux
!                                      ! between Cc and Ccomp_tot (C):
!                                      !
!                                      !  (Cc - Cw  )    (Cc - Csoil)    (Cc - Cstom)   (Cc - Ccomp_tot)
!                                      !  ------------ + ------------ +  ------------ = ----------------  <=> (Cc-terms cancel)
!                                      !     Rw         Rsoil_eff         Rstom             Rc
!                                      !
!                                      !  1                1              1              1
!                                      ! --- Ccomp_tot = ---- Cw + --------- Csoil + ----- Cstom.
!                                      !  Rc             Rw        Rsoil_eff         Rstom
!                                      !
! F = -vd*Catm                         ! F = -ve*[Catm - Ccomp_tot]
!                                      !
! Mass balance:                        ! Mass balance:
!                                      !
!   dCatm                              !   dCatm
! H ----- = F = -vd*Catm,              ! H ----- = F = -ve*[Catm - Ccomp_tot],
!    dt                                !    dt
!                                      !
! with solution:                       ! with solution (assuming Ccomp_tot constant):
!                                      !
! Catm(t+dt) = Catm(t)*exp(-(vd/H)*dt) ! Catm(t+dt) = Ccomp_tot + [Catm(t) - Ccomp_tot]*exp(-(ve/H)*dt)
!
! Approach D uses the same scheme as C (compensation point), but computes a
! effective total resistance Rc_eff, such that the differential equation of approach A can be used
! (see rc_comp_point_rc_eff).
!
! Note that DEPAC is restricted to the computation of the canopy or surface resistance Rc (A)
! and (optionally) a total compensation point (B/C) and/or an effective Rc (D).
! Ra, Rb, deposition fluxes and new concentrations are computed outside DEPAC.
!
!  Different options available for:
!  -  calculating bulk stomatal resistance
!  -  calculating stomatal compensation point for NH3
!  -  parameterisation of external leaf resistance NH3
! see ipar_ variables below
!
!  -  Wetness indicator
!     nwet=0 dry
!     nwet=1 wet
!     nwet=9 snow
!
!  -  NH3/SO2 ratio
!     ratns=1 low ratio NH3/SO2
!     ratns=2 high ratio NH3/SO2
!     ratns=3 very low ratio NH3/SO2
!    see calling routine:
!     [SO2] < 0.5 then
!        iratns = 1 ! low
!     else
!        [NH3]/[SO2] < 0.02         -> iratns = 3 ! very low
!        [NH3]/[SO2] > 2            -> iratns = 2 ! high
!        0.02 <= [NH3]/[SO2] <= 2   -> iratns = 1 ! low
!     endif
!
!  -  DEPAC is used for the following gaseous components:
!     HNO3, NO, NO2, O3, SO2, NH3.
!     In this version of DEPAC, component numbers are not used anymore, only names
!     (only upper case!).
!
! In order to prevent problems with typing errors in the call to DEPAC or
! when new species are introduced into DEPAC,
! the following CASE construct is used throughout this module; this means that
! if a new species is introduced, all these CASE constructs have to be wended.
! select case(trim(compnam))
! case('HNO3','N2O5','NO3','H2O2')
! case('NO','CO')
! case('NO2')
! case('O3')
! case('SO2')
! case('NH3')
! case default
!    print *, 'error in subroutine ... '
!    print *, 'component ',trim(compnam),' not supported'
!    stop
! end select
!
! For convenience, case entries may be combined, e.g. case('NO','NO2').
! In fact, this type of CASE construct is the only allowed construct where compnam may occur.
!
!****************************************************************************************
!     missing deposition path
!****************************************************************************************
! A missing deposition path (e.g. deposition via stomata for land use "water") is
! represented by parameter values of -999;
! a logical function "missing" is available to check whether a parameter value
! corresponds to a missing deposition path or not.
! The effect of a missing deposition path is that the corresponding conductance is set to 0.
!
! Note that in previous versions of DEPAC, the value -1000 was also (mis)used to
! represent a missing deposition path (probably to avoid rounding error problems
! when checking parameter_value < -999). This use of -1000 has been banned; instead
! the logical function "missing"  uses an EPS (i.e. small value) band to avoid these
! problems with rounding errors.
!
! Missing output data for rc_tot is given a value of -9999; this means that DEPAC was not
! able to calculate a value for the total canopy resistance and that there is no deposition;
! the calling routine should check for this output value of -9999!
!
!****************************************************************************************
!     Land use types
!****************************************************************************************
!     DEPAC is based on RIVM/LBG land use data, 8 classes, plus 1 extra class
!     for special use:
!     1 = grass
!     2 = arable land
!     3 = permanent crops
!     4 = coniferous forest
!     5 = deciduous forest
!     6 = water
!     7 = urban
!     8 = other i.e. short grassy area
!     9 = desert  !used for N-Africa
!
! For the leaf area index LAI and stomatal resistance according to Emberson,
! land uses classes of Emberson are used:
! temp_conif, temp_decid, med_needle, med_broadleaf, temp_crop, med_crop, root_crop, moorland, grass,
! medscrub, wetlands, tundra, desert, water, ice, urban.
!
! Translation tables for leaf area index and for stomatal resistance (Rs) parameters are somewhat different,
! since the Emberson parameterisation of the growing season (needed for the LAI) for temperate crop
! was considered to be too short.
!
!         RIVM/LBG                        LAI Emberson        Rs Emberson
!     1 = grass                           grass               grass
!     2 = arable                          root_crop           temp_crop
!     3 = permanent crops                 root_crop           temp_crop
!     4 = coniferous forest               temp_conif          temp_conif
!     5 = deciduous forest                temp_decid          temp_decid
!     6 = water                           water               water
!     7 = urban                           urban               urban
!     8 = other i.e. short grassy area    grass               grass
!     9 = desert                          desert              desert
!
! Note that in this version of DEPAC all translations have been done already
! and are directly seen in the appropriate values of the parameters used.
!
!  UPDATE HISTORY :
!    1994    , article Erisman & van Pul, Atm. Env.
!    ?       , Franka Loeve (Cap Volmac)
!    Jan 2003, Martien de Haan (ARIS): made single depac module.
!    ?       ,               ? (TNO) : added rc for O3
!    ?       ,               ? (TNO) : separate routines for each species.
!    Nov 2008, Ferd Sauter     (RIVM): v3.0 synthesis of OPS and LOTOS-EUROS versions of DEPAC
!    v3.0      model structure improved; common tasks in separate routines; documentation added
!              names have been changed for readability:
!
!              old name -> new name     : explanation
!              ----------------------------------------
!              rinc     -> rinc         : in canopy resistance (s/m)
!              rso      -> rsoil        : soil resistance, dry soil (s/m)
!              rswet    -> rsoil_wet    : soil resistance, wet soil (s/m)
!              rsfrs    -> rsoil_frozen : soil resistance, frozen soil (s/m)
!              rsoeff   -> rsoil_eff    : effective soil resistance, rsoil_eff = rinc + rsoil (s/m)
!              rext     -> rw           : external leaf resistance (s/m)
!              rs       -> rstom        : stomatal resistance (s/m)
!              rssnow   -> rssnow       : total canopy resistance Rc in case of snow (s/m)
!              rstot    -> rc_tot       : total canopy resistance Rc (s/m)
!              gscu     -> gw           : external leaf (or cuticular) conductance (m/s)
!              gsoeff   -> gsoil_eff    : effective soil conductance (m/s)
!              gs       -> gstom        : stomatal conductance (m/s)
!              gstot    -> gc_tot       : total canopy conductance (m/s)
!
!    10 Dec 2008, Ferd Sauter     (RIVM): v3.1 try-out version with new model structure;
!    v3.1         no calls to separate routines in subroutine depac, but all components
!                 are dealt with in subroutine depac.
!                 This version is NOT developed any further; depacv3.2 is developed from v3.0
!
!    11 Dec 2008, Ferd Sauter     (RIVM): v3.2 new model structure;
!    v3.2         in subroutine depac, calls are made to routines for separate conductances, e.g.
!                 for external, stomatal, soil conductance; the dependence on the components is
!                 placed inside these conductance-routines.
!
!    22 Jan 2009, Ferd Sauter     (RIVM): v3.3 bug fix in season dependency leaf area index;
!    v3.3         see function rc_lai. Older versions of this routine use a wrong numbering of
!                 land use types (no conversion to Olson land use types).
!                 rc_gstom_wes (Wesely) readability improved; routine gives the same results.
!
!    03 Feb 2009, Ferd Sauter     (RIVM): v3.4 Rsoil(NH3,urban) = 100 (was 1000).
!    v3.4
!
!    03 Feb 2009, Ferd Sauter     (RIVM): v3.5 Rinc(grass) = Inf (was 0).
!    v3.5
!
!    03 Feb 2009, Ferd Sauter     (RIVM): v3.6 stomatal compensation point and
!    v3.6         new parameterisation Rw.
!                 New routines:
!                    rc_comp_point (called from depac)
!                    rc_comp_point_rc_eff (called from depac)
!                    rw_nh3_sutton (called from rc_gw)
!                 New option ipar_rw_nh3.    (obsolete in final version MCvZ Nov 2009)
!                 New (optional) arguments of depac: see header of depac.
!
!    02 Mar 2009, Ferd Sauter     (RIVM): v3.7
!    v3.7         - added compensation point for external leaf;                                                  
!                   new parameterisation for Rw (routine rw_nh3_sutton replaces rw_nh3_rwk);                     
!                 - added compensation point for soil; value of compensation point                               
!                   set to zero, due to lack of data.                                                            
!                 - the same parameterisations as in v3.6 for gamma_stom (only new variable name)                
!                 - Wesely for Rstom; note that the parameterisation for Rw is deduced using the Baldocchi       
!                   parameterisation of Rstom; it was decided to change not everything in this version           
!                   but to do it stepwise. See v3.8 for Baldocchi                                                
!
!    10 Mar 2009, Ferd Sauter     (RIVM): v3.8
!    v3.8         - the same as v3.7, but Baldocchi for Rstom
!
!    24 Mar 2009, Ferd Sauter     (RIVM): v3.8.1 LAI in external leaf resistance
!    v3.8.1       gw = (lai/lai(grass)) * gw   (adjusted Oct 2009; lai -> sai and 
!                 sai_grass scaling inside rw_nh3_sutton routine)
!
!    9 Apr 2009, Ferd Sauter     (RIVM): v3.8.2 bug fix in temperature
!    v3.8.2       correction factor Baldocchi BT
!
!    6 July 2009, Ferd Sauter    (RIVM): v3.9 call added to calculate Rstom with Emberson
!    v3.9
!
!    13 Aug 2009, Margreet van Zanten (RIVM): Emberson update, PARshade and PARsun added
!                 (Zhang et al., 2001 Ae 4463-4470 which is update of Norman, 1982
!
!    17 Aug 2009, Margreet van Zanten (RIVM): update gamma parametrization according to
!                 submitted version of Wichink Kruit, 2009  (implemented eq 15 and 14a + eq 12)
!
!    9  Sep 2009, Margreet van Zanten (RIVM): calc of PARdir and PARdiff in Emberson
!                 according to Weiss and Norman 1985
!
!    22 Sep 2009, Ferd Sauter (RIVM): gstom of Emberson scaled with diffc/dO3 (instead of
!    v3.10        erroneously with dwat)
!
!    24 Sep 2009, Margreet van Zanten (RIVM): choices made on lu classes, F_phen set to 1
!                  since described effect is negligible for chosen lu's
!
!    29 Sep 2009, Ferd Sauter (RIVM): Emberson parameterisation of leaf area index (rc_lai);
!                  new subroutine arguments for DEPAC: day_of_year and lat (latitude).
!
!    29 Sep 2009, Arjo Segers (TNO), routine "Rc_GStom_Emb"
!                 Limitted temperature influence to range to avoid 
!                 floating point exceptions.
!
!    2 Oct  2009, Margreet van Zanten (RIVM): v3.10 Merged version of earlier version of 3.10 and 3.9.2
!                 Contains soil compensation point. Currently csoil = 0 for all lu's except lu = ilbg_water
!                 for which a basic parametrization for the calc. of csoil is added.
!                 In this parametrization no disctinction is made between salt and fresh water.
!
!    6 Oct 2009, Margreet van Zanten (RIVM): Distinction made between area index for stom. resistance (LAI)
!                and external resistance Rw (SAI = LAI + area index of stems and branches).
!                rw (external leaf resistance) scaled with fixed sai_grass value (valid for Haarweg data) instead of Gw
!                (external canopy resistance). Only Rw_nh3_sutton is scaled not Rw when it is freezing
!                Old logical vegetetation_present is splitted into LAI_present and SAI_present.
!
!    4 Nov 2009, Margreet van Zanten (RIVM): v3.11, all obsolete code (Wesely, Baldocchi, old LAI parameterisations)
!    v3.11       removed; functionally identical to DEPAC v3.10.
!
!    6 Nov 2009, Margreet van Zanten (RIVM): Differences between LE and OPS version of DEPAC straightened. 
!    v3.11       Choice for either LE or OPS version, guided by description of DEPAC in Erisman et al, 1994.
!                Rw for NO set to -9999. instead of either LE or OPS option
!                Rc calculated for lu = ilbg_other is once more identical to the one for grass (rinc is missing instead of 0)
!
!    6 Nov 2009, version v3.11 renamed in depac.f90
!
!      nov 2009, Arjo Segers (TNO)
!                Added routine "Get_N_S_Regime" to compute N/S ratio.
!
!    25 Nov 2009, Margreet van Zanten (RIVM): bug fix for LE implementation
!    v3.11       ccomp_tot set to zero in rc_special to avoid use of outdated value of ccomp_tot when 
!                calling depac routine for several components in a row (esp. NO after NH3), 
!                ccomp_tot added as optional argument to rc_special routine
!
!   18 Dec 2009, Margreet van Zanten (RIVM): SAI for arable land (lu = ilbg_arable) set to 0.5 outside growing season 
!   v3.13        (instead of 0)
!
!   18 Dec 2009, depac.f90 renamed in depac_2010.f90
!
!   22 Dec 2009, Margreet van Zanten (RIVM): v3.14 -> copy of depac_2010.f90. Write statements for test output added 
!   v3.14        from depac_313. Frozen version
!
!   22 Dec 2009, Margreet van Zanten (RIVM): v3.15 copy of 3.14, two comments related to Rinc added (see 3.13)
!   v3.15         small bug fixed in rw_constant and rc_rctot so that rw can be set to -9999. correctly
!                
!   4 Jan 2010,  Margreet van Zanten (RIVM): depac_f90 and depac_3.11 merged, thus following bug fix dated Nov 25 2009 
!   v3.11        in depac.f90 implemented in 3.11
!                bug fix for LE implementation
!                ccomp_tot set to zero in rc_special to avoid use of outdated value of ccomp_tot when 
!                calling depac routine for several components in a row (esp. NO after NH3), 
!                ccomp_tot added as optional argument to rc_special routine
!
!   04 Jan 2010, Ferd Sauter (RIVM): v3.16 is shell around versions 3.11 ('new' DEPAC for NH3 only) 
!   v3.16        and 3.3 (old DEPAC for other species).       
!                This file is constructed as follows:
!                module m_depac311
!                module m_depac33
!                module m_depac316
!
!   04 Jan 2010, Ferd Sauter (RIVM): iopt_debug -> optional writing of debug output
!   v3.16        added to m_depac311 and m_depac33 in this file 
!
!   04 Jan 2010, Margreet van Zanten(RIVM): frozen version of depac v3.16, renamed in depac_GCN2010
!   depac_GCN2010
!
!   12 Sep 2013, Roy Wichink Kruit(TNO): co-deposition version included and tested in LOTOS-EUROS and extra (optional) output variable added
!             
!   2013-09-17   this version has been derived from the 'hybrid' version 
!   v3.18        depac_GCN2010, which consisted of a shell around version 
!                depac311 (for NH3) and depac33 (for other species).
!                In this version, only depac311 has been retained, with some 
!                bug fixes, FS resulting in version depac318.  
!
!   2015-03-11   synchronized version between LE and OPS DEPAC versions with extra optional in-/output arguments
!   v3.19
!
!   2015-03-11   version with co-deposition of SO2/NH3 for NH3 and extra optional arguments used in LE
!   v3.20
!                
!   2015-03-23   LAI calculation made consistent with EMEP (Simpson et al. 20XX)
!   v3.21
!                
!   2016-08-04   Replaced module variables 'LAI_present', 'SAI_present', and 'icmp'
!                by local variables and arguments to make the routines thread-save;
!                this fix has been introduced before (2015), but was lost somewhere after.
!                Removed suffixes "_3_11", version number is in comment now.
!                (Arjo Segers, TNO)
!
!   2019-03-08   Ruud Janssen (TNO): Added dependency on solubility to calculation of surface resistances
!                to allow for deposition of semi-volatile organic compounds (vbs_cg's in the model).
!                Scaling with Henry's law constants follows Wesely 1989    
!                Also see papers by Hodzic et al., 2014, GRL and Knote et al., 2015, ACP
! 
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################

module LE_DryDepos_Gas_DEPAC

  use GO, only : gol, goPr, goErr
  use modprecision, only : field_r
  
  implicit none
  
  ! --- in/out -----------------------------------

  ! by default all variables and procedures private:
  private
  
  ! main routine:
  public    ::  DryDepos_Gas_DEPAC
  
  ! input function:
  public    ::  Get_N_S_Regime
  public    ::  iratns_default
  
  
  ! --- const ------------------------------------  

  character(len=*), parameter   ::  mname = 'LE_DryDepos_Gas_DEPAC'

! Include constants and landuse dependent parameters from standard input files
! These files are also used by LE_Landuse_data for the rest of the model  
! Parameters include:
!   rsoil   (nlu,ncmp)
!   rsoil_wet   (ncmp)
!   rsoil_frozen(ncmp)
!   diffc       (ncmp)
! Do not pass these as arguments too, this might give strange results ..
#include "depac_lu.inc"

  ! NH3/SO2 ratio regimes:
  integer, parameter  ::  iratns_low      = 1 ! low ratio NH3/SO2
  integer, parameter  ::  iratns_high     = 2 ! high ratio NH3/SO2
  integer, parameter  ::  iratns_very_low = 3 ! very low ratio NH3/SO2
  ! default:
  integer, parameter  ::  iratns_default = iratns_low


contains


!**********************************************************************************
!
! main routine depac
!
!**********************************************************************************

!-------------------------------------------------------------------
! depac: compute total canopy (or surface) resistance Rc for gases
!-------------------------------------------------------------------
subroutine DryDepos_Gas_DEPAC( compnam, day_of_year, lat, t, ust, glrad, sinphi, rh, &
                                 lai, sai, nwet, lu, iratns, &
                                 rc_tot, ccomp_tot, hlaw, react, &
                                 status, &
                                 p, tsea, smi, &
                                 c_ave_prev_nh3, c_ave_prev_so2, catm, gamma_soil_water, &
                                 ra, rb, rc_eff, &
                                 gw_out, gstom_out, gsoil_eff_out, cw_out, cstom_out, csoil_out, &
                                 calc_stom_o3flux, frac_sto_o3_lu, fac_surface_area_2_PLA )

!!DEC$ ATTRIBUTES DLLEXPORT:: DryDepos_Gas_DEPAC

! The last four rows of depac arguments are optional:
!
! A. compute Rc_tot without compensation points (ccomp_tot will be zero):
!     call depac (compnam, day_of_year, lat, t, ust, glrad, sinphi, rh, nwet, lu, iratns, rc_tot, ccomp_tot, [smi])
!
! B. compute Rc_tot with compensation points (used for LOTOS-EUROS):
!     call depac (compnam, day_of_year, lat, t, ust, glrad, sinphi, rh, nwet, lu, iratns, rc_tot, ccomp_tot, [smi], &
!                  c_ave_prev_nh3, c_ave_prev_so2, catm, gamma_soil_water)
!
! C. compute effective Rc based on compensation points (used for OPS):
!     call depac (compnam, day_of_year, lat, t, ust, glrad, sinphi, rh, nwet, lu, iratns, rc_tot, ccomp_tot, [smi], &
!                 c_ave_prev_nh3, c_ave_prev_so2, catm, gamma_soil_water, &
!                 ra, rb, rc_eff)
! X1. Extra (optional) output variables: 
!     call depac (compnam, day_of_year, lat, t, ust, glrad, sinphi, rh, nwet, lu, iratns, rc_tot, ccomp_tot, [smi], &
!                 c_ave_prev_nh3, c_ave_prev_so2, catm, gamma_soil_water, &
!                 ra, rb, rc_eff, &
!                 gw_out, gstom_out, gsoil_eff_out, cw_out, cstom_out, csoil_out, lai_out, sai_out)
! X2. Extra (optional) needed for stomatal ozone flux calculation (only sunlit leaves): 
!     call depac (compnam, day_of_year, lat, t, ust, glrad, sinphi, rh, nwet, lu, iratns, rc_tot, ccomp_tot, [smi], &
!                 c_ave_prev_nh3, c_ave_prev_so2, catm, gamma_soil_water, &
!                 ra, rb, rc_eff, &
!                 gw_out, gstom_out, gsoil_eff_out, cw_out, cstom_out, csoil_out, lai_out, sai_out, &
!                 calc_stom_o3flux, frac_sto_o3_lu, fac_surface_area_2_PLA)

implicit none

character(len=*) , intent(in)  :: compnam          ! component name
                                                   ! 'HNO3','NO','NO2','O3','SO2','NH3'
integer          , intent(in)  :: day_of_year      ! day of year, 1 ... 365 (366)
real             , intent(in)  :: lat              ! latitude Northern hemisphere (degrees) (DEPAC cannot be used for S. hemisphere)
real             , intent(in)  :: t                ! temperature (C)
                                                   ! NB discussion issue is temp T_2m or T_surf or T_leaf?
real             , intent(in)  :: ust              ! friction velocity (m/s)
real(field_r)    , intent(in)  :: glrad            ! global radiation (W/m2)
real             , intent(in)  :: sinphi           ! sin of solar elevation angle
real             , intent(in)  :: rh               ! relative humidity (%)
real             , intent(in)  :: lai              ! one-sidedleaf area index (-)
real             , intent(in)  :: sai              ! surface area index (-) (lai + branches and stems)
integer          , intent(in)  :: nwet             ! wetness indicator; nwet=0 -> dry; nwet=1 -> wet; nwet=9 -> snow
integer          , intent(in)  :: lu               ! land use type, lu = 1,...,nlu
integer          , intent(in)  :: iratns           ! index for NH3/SO2 ratio used for SO2;
                                                   ! iratns = 1: low NH3/SO2
                                                   ! iratns = 2: high NH3/SO2
                                                   ! iratns = 3: very low NH3/SO2
real             , intent(out) :: rc_tot           ! total canopy resistance Rc (s/m)
real             , intent(out) :: ccomp_tot        ! total compensation point (ug/m3) [= 0 for species that don't have a compensation 
                                                   !   point or if arguments (c_ave_prev_nh3,c_ave_prev_so2,catm,gamma_soil_water) are 
                                                   !   missing]
real             , intent(in)  :: hlaw            ! effective Henry's law constant (M atm-1)
real             , intent(in)  :: react           ! reactivity (-)
integer          , intent(out) :: status

real, optional   , intent(in)  :: p                ! pressure (Pa)
real, optional   , intent(in)  :: tsea             ! sea surface temperature (K)
real, optional   , intent(in)  :: smi              ! soil moisture index (there is still a problem with the soil moisture over sea =0 
                                                   ! in the ECMWF data and the interpolation close to the coast)
! optional arguments needed only if compensation points are computed
real, optional   , intent(in)  :: c_ave_prev_nh3   ! air concentration averaged over a previous
                                                   ! period (e.g. previous year or month) (ug/m3)
real, optional   , intent(in)  :: c_ave_prev_so2   ! air concentration averaged over a previous
                                                   ! period (e.g. previous year or month) (ug/m3)
real(field_r), optional, intent(in)  :: catm             ! actual atmospheric concentration (ug/m3)
real, optional   , intent(in)  :: gamma_soil_water ! gammawater from gammawater.nc file
! optional arguments needed only if an effective Rc (based on compensation points) is computed;
! in this case, also the previous three optional arguments are needed
real, optional   , intent(in)  :: ra              ! aerodynamic resistance (s/m)
real, optional   , intent(in)  :: rb              ! boundary layer resistance (s/m)
real, optional   , intent(out) :: rc_eff          ! effective total canopy resistance (s/m)
! optional output
real, optional   , intent(out) :: gw_out          ! external leaf surface conductance (m/s)
real, optional   , intent(out) :: gstom_out       ! stomatal conductance (m/s)
real, optional   , intent(out) :: gsoil_eff_out   ! effective soil conductance (m/s)
real, optional   , intent(out) :: cw_out          ! external leaf surface compensation point (ug/m3)
real, optional   , intent(out) :: cstom_out       ! stomatal compensation point (ug/m3)
real, optional   , intent(out) :: csoil_out       ! soil compensation point (ug/m3)
! optional arguments used for special calculation of stomatal ozone flux
logical, optional, intent(in)  :: calc_stom_o3flux! calculate stomatal ozone flux? (True/False)
real, optional   , intent(out) :: frac_sto_o3_lu  ! contribution of stomatal conductance to total conductance (-) 
real, optional   , intent(out) :: fac_surface_area_2_PLA ! factor to convert from m2 surface area to m2 Projected Leaf Area (PLA) (-)

! constants:
character(len=*), parameter   ::  rname = mname//'/DryDepos_Gas_DEPAC'

! Local variables:
real                           :: laimax          ! maximum leaf area index (-)
logical                        :: ready           ! Rc has been set
                                                  ! = 1 -> constant Rc
                                                  ! = 2 -> temperature dependent Rc
real                           :: gw              ! external leaf conductance (m/s)
real                           :: gstom           ! stomatal conductance (m/s)
real                           :: gmes            ! mesophyll conductance (m/s)
real                           :: gsoil_eff       ! effective soil conductance (m/s)
real                           :: gc_tot          ! total canopy conductance (m/s)
real                           :: cw              ! external leaf surface compensation point (ug/m3)
real                           :: cstom           ! stomatal compensation point (ug/m3)
real                           :: csoil           ! soil compensation point (ug/m3)
 
! Vegetation indicators:
logical            :: LAI_present ! leaves are present for current land use type
logical            :: SAI_present ! vegetation is present for current land use type
 
! Component number taken from component name, paramteres matched with include files
integer           ::  icmp
  
! Define component number
select case ( trim(compnam) )
  
  case ( 'O3', 'o3' )
    icmp = icmp_o3
  
  case ( 'SO2', 'so2' )
    icmp = icmp_so2
  
  case ( 'NO2', 'no2' )
    icmp = icmp_no2
  
  case ( 'NO', 'no' )
    icmp = icmp_no 
  
  case ( 'NH3', 'nh3' )
    icmp = icmp_nh3
  
  case ( 'CO', 'co' )
    icmp = icmp_co
  
  case ( 'NO3', 'no3' )
    icmp = icmp_no3
  
  case ( 'HNO3', 'hno3' )
    icmp = icmp_hno3
  
  case ( 'N2O5', 'n2o5' )
    icmp = icmp_n2o5

  case ( 'H2O2', 'h2o2' )
    icmp = icmp_h2o2

  case ( 'VBS_POG1', 'vbs_pog1' )
    icmp = icmp_vbs_pog1

  case ( 'VBS_POG2', 'vbs_pog2' )
    icmp = icmp_vbs_pog2

  case ( 'VBS_POG3', 'vbs_pog3' )
    icmp = icmp_vbs_pog3

  case ( 'VBS_POG4', 'vbs_pog4' )
    icmp = icmp_vbs_pog4

  case ( 'VBS_POG5', 'vbs_pog5' )
    icmp = icmp_vbs_pog5

  case ( 'VBS_POG6', 'vbs_pog6' )
    icmp = icmp_vbs_pog6

  case ( 'VBS_POG7', 'vbs_pog7' )
    icmp = icmp_vbs_pog7

  case ( 'VBS_POG8', 'vbs_pog8' )
    icmp = icmp_vbs_pog8

  case ( 'VBS_POG9', 'vbs_pog9' )
    icmp = icmp_vbs_pog9

  case ( 'VBS_SISOG1', 'vbs_sisog1' )
    icmp = icmp_vbs_sisog1

  case ( 'VBS_SISOG2', 'vbs_sisog2' )
    icmp = icmp_vbs_sisog2

  case ( 'VBS_SISOG3', 'vbs_sisog3' )
    icmp = icmp_vbs_sisog3

  case ( 'VBS_SISOG4', 'vbs_sisog4' )
    icmp = icmp_vbs_sisog4

  case ( 'VBS_SISOG5', 'vbs_sisog5' )
    icmp = icmp_vbs_sisog5

  case ( 'VBS_SISOG6', 'vbs_sisog6' )
    icmp = icmp_vbs_sisog6

  case ( 'VBS_SISOG7', 'vbs_sisog7' )
    icmp = icmp_vbs_sisog7

  case ( 'VBS_SISOG8', 'vbs_sisog8' )
    icmp = icmp_vbs_sisog8

  case ( 'VBS_ASOG1', 'vbs_asog1' )
    icmp = icmp_vbs_asog1

  case ( 'VBS_ASOG2', 'vbs_asog2' )
    icmp = icmp_vbs_asog2

  case ( 'VBS_ASOG3', 'vbs_asog3' )
    icmp = icmp_vbs_asog3

  case ( 'VBS_ASOG4', 'vbs_asog4' )
    icmp = icmp_vbs_asog4

  case ( 'VBS_ASOG5', 'vbs_asog5' )
    icmp = icmp_vbs_asog5

  case ( 'VBS_ASOG6', 'vbs_asog6' )
    icmp = icmp_vbs_asog6

  case ( 'VBS_BSOG1', 'vbs_bsog1' )
    icmp = icmp_vbs_bsog1

  case ( 'VBS_BSOG2', 'vbs_bsog2' )
    icmp = icmp_vbs_bsog2

  case ( 'VBS_BSOG3', 'vbs_bsog3' )
    icmp = icmp_vbs_bsog3

  case ( 'VBS_BSOG4', 'vbs_bsog4' )
    icmp = icmp_vbs_bsog4

  case ( 'VBS_BSOG5', 'vbs_bsog5' )
    icmp = icmp_vbs_bsog5

  case ( 'VBS_BSOG6', 'vbs_bsog6' )
    icmp = icmp_vbs_bsog6
  case default
    write (gol,'("unsupported component `",a,"`")') trim(compnam); call goErr
    TRACEBACK; status=1; return
  
end select

! inititalize
gw        = 0.
gstom     = 0.
gsoil_eff = 0.
gc_tot    = 0.
cw        = 0.
cstom     = 0.
csoil     = 0.
 
! Compute one-sided leaf area index:
laimax = lai_par(lu)%laimax

! Check whether vegetation is present (in that case the leaf or surface area index > 0):
LAI_present = (lai .gt. 0.0)
SAI_present = (sai .gt. 0.0)

! Set Rc (i.e. rc_tot) in special cases:
!   for NH3, the 'ready' flag is .true. if nwet==9 (snow covered)
call rc_special(icmp,compnam,lu,t,ipar_snow,nwet,rc_tot,ready,ccomp_tot)

! set effective resistance equal to total resistance, this will be changed for compensation points
if ( present(rc_eff) ) then
  rc_eff = rc_tot
end if

! If Rc is not set:
if (.not. ready) then

   ! External conductance:
   call rc_gw(compnam,iratns,t,rh,nwet,SAI_present,sai,gw)         

   ! Stomatal conductance:
   call rc_gstom( icmp, compnam, lu, LAI_present, lai, &
                    glrad, sinphi, t, rh, hlaw, react, gstom, &
                    p=p, smi=smi,calc_stom_o3flux=calc_stom_o3flux )

   ! Mesophyll conductance: 
   call rc_mes( icmp, compnam, hlaw, react, gmes )

   ! Effective soil conductance:
   call rc_gsoil_eff( icmp, compnam, lu, sai, ust, nwet, t, &
                        hlaw, react, icmp_o3, icmp_so2, gsoil_eff, status )
   IF_NOTOK_RETURN(status=1)

   ! Total canopy conductance (gc_tot) and resistance Rc (rc_tot):
   call rc_rctot(gstom,gsoil_eff,gw,gmes,gc_tot,rc_tot)

   ! Compensation points (always returns ccomp_tot; currently ccomp_tot != 0 only for NH3):
   call rc_comp_point( compnam,lu,day_of_year,t,gw,gstom,gsoil_eff,gc_tot,&
             LAI_present, SAI_present, &
             ccomp_tot, status, &
             catm=catm,c_ave_prev_nh3=c_ave_prev_nh3, &
               c_ave_prev_so2=c_ave_prev_so2,gamma_soil_water=gamma_soil_water, &
               tsea=tsea,cw=cw,cstom=cstom,csoil=csoil )
   IF_NOTOK_RETURN(status=1)

   ! Effective Rc based on compensation points:
   if ( present(rc_eff) ) then
      ! check on required arguments:
      if ( (.not. present(catm)) .or. (.not. present(ra)) .or. (.not. present(rb)) ) then
        stop 'output argument rc_eff requires input arguments catm, ra and rb'
      end if
      ! compute rc_eff :
      call rc_comp_point_rc_eff( ccomp_tot, catm, ra, rb, rc_tot, rc_eff, status )
      IF_NOTOK_RETURN(status=1)
   endif
endif

! optional extra output
if ( present( gw_out ) ) then
  gw_out = gw
end if

if ( present( gstom_out ) ) then
  gstom_out = gstom
end if

if ( present( gsoil_eff_out ) ) then
  gsoil_eff_out = gsoil_eff
end if

if ( present( cw_out ) ) then
  cw_out = cw
end if

if ( present( cstom_out ) ) then
  cstom_out = cstom
end if

if ( present( csoil_out ) ) then
  csoil_out = csoil
end if

if ( present ( frac_sto_o3_lu ) ) then
  ! determine contribution of stomatal conductance to total conductance:
  if ( gc_tot .le. 0.0 ) then
    frac_sto_o3_lu = 0.0
  else 
    frac_sto_o3_lu = gstom/gc_tot
  end if
end if  

if ( present( fac_surface_area_2_PLA ) ) then
  ! Stomatal conductance gstom is in m/s per m2 surface area, 
  ! For the stomatal flux calculations (PODy) a flux per Projected Leaf Area (PLA) is required:
  !   1 m2 surface area = LAI m2 total leaf area = 2*LAI m2 PLA
  ! Note that 1 m2 leaf area = 2 * PLA as leaves are considered to be 2-sided
  ! So, 1 m2 PLA = LAI/2
  fac_surface_area_2_PLA = 0.0             ! unit: [-]    
  if ( lai > 0.0 ) then                                                                               
    fac_surface_area_2_PLA = 2.0 / laimax     ! conversion from m2 surface area to m2 PLA.                                     
  else                                                                                                
    fac_surface_area_2_PLA = 0.0                                                                                     
  end if    
end if
!! extra output !!

! ok
status = 0

end subroutine DryDepos_Gas_DEPAC

!************************************************************************************
! help routines
!************************************************************************************


!-------------------------------------------------------------------
! rc_special: compute total canopy resistance in special cases
!-------------------------------------------------------------------
subroutine rc_special(icmp,compnam,lu,t,ipar_snow,nwet,rc_tot,ready,ccomp_tot)

integer         , intent(in)  :: icmp            ! component index
character(len=*), intent(in)  :: compnam         ! component name
integer         , intent(in)  :: lu              ! land use type, lu = 1,...,nlu
real            , intent(in)  :: t               ! temperature (C)
integer         , intent(in)  :: ipar_snow(ncmp) ! parameterisation in case of snow:
                                                 ! = 1 -> constant Rc
                                                 ! = 2 -> temperature dependent Rc
integer         , intent(in)  :: nwet            ! wetness indicator; nwet=0 -> dry; nwet=1 -> wet; nwet=9 -> snow
real            , intent(out) :: rc_tot          ! total canopy resistance Rc (s/m)
logical         , intent(out) :: ready           ! Rc has been set
real            , intent(out) :: ccomp_tot       ! total compensation point (ug/m3)

! rc_tot is not yet set:
ready = .false.

! Default compensation point in special cases = 0:
ccomp_tot = 0.0

select case(trim(compnam))
case('HNO3','N2O5','NO3','H2O2')
   ! No separate resistances for HNO3; just one total canopy resistance:
   if (t .lt. -5.0 .and. nwet .eq. 9) then
     ! T < 5 C and snow:
      rc_tot = 50.
   else
     ! all other circumstances:
      rc_tot = 10.0
   endif
   ready = .true.

case('NO','CO')
      if (lu .eq. ilu_water_sea .or. lu .eq. ilu_water_inland) then       ! water
         rc_tot = 2000.
         ready = .true.
      elseif (nwet .eq. 1) then ! wet
         rc_tot = 2000.
         ready = .true.
      elseif (nwet .eq. 9) then ! snow
         call rc_snow(ipar_snow(icmp),t,rc_tot)
         ready = .true.
      endif

case('NO2','O3','SO2','NH3')
   ! snow surface:
   if (nwet.eq.9) then
      call rc_snow(ipar_snow(icmp),t,rc_tot)
      ready = .true.
   endif

case default
  if ( compnam(1:3) == 'VBS' ) then
    ! snow surface:
    if (nwet.eq.9) then
      call rc_snow(2,t,rc_tot) ! for now, make rc_snow temnperature dependent
      ready = .true.
    end if
  else
   print *, 'error in subroutine rc_special '
   print *, 'component ',trim(compnam),' not supported'
   stop
  end if
end select

end subroutine rc_special

!-------------------------------------------------------------------
! rc_gw: compute external conductance
!-------------------------------------------------------------------
subroutine rc_gw( compnam, iratns,t,rh,nwet, SAI_present, sai, gw )

! Input/output variables:
character(len=*), intent(in)  :: compnam ! component name
integer         , intent(in)  :: iratns  ! index for NH3/SO2 ratio;
                                         ! iratns = 1: low NH3/SO2
                                         ! iratns = 2: high NH3/SO2
                                         ! iratns = 3: very low NH3/SO2
real            , intent(in)  :: t       ! temperature (C)
real            , intent(in)  :: rh      ! relative humidity (%)
integer         , intent(in)  :: nwet    ! wetness indicator; nwet=0 -> dry; nwet=1 -> wet; nwet=9 -> snow
logical         , intent(in)  :: SAI_present
real            , intent(in)  :: sai     ! one-sided leaf area index (-)
real            , intent(out) :: gw      ! external leaf conductance (m/s)

select case(trim(compnam))

!case('HNO3','N2O5','NO3','H2O2') this routine is not called for HNO3

case('NO2')
   call rw_constant(2000.,SAI_present,gw)

case('NO','CO')
   call rw_constant(-9999.,SAI_present,gw)  ! see Erisman et al, 1994 section 3.2.3

case('O3')
   call rw_constant(2500.,SAI_present,gw)

case('SO2')
   call rw_so2( t, nwet, rh, iratns, SAI_present, gw )

case('NH3')
   call rw_nh3_sutton(t,rh,SAI_present,gw)

   ! conversion from leaf resistance to canopy resistance by multiplying with SAI:
   Gw = sai*gw

case default
  if ( compnam(1:3) == 'VBS' ) then
    !rj todo: realistic rw for svoc
    call rw_constant(2000.,SAI_present,gw)
    ! conversion from leaf resistance to canopy resistance by multiplying with SAI:
    Gw = sai*gw
  else
   print *, 'error in subroutine rc_gw '
   print *, 'component ',trim(compnam),' not supported'
   stop
  end if
end select

end subroutine rc_gw

!-------------------------------------------------------------------
! rw_so2: compute external leaf conductance for SO2
!-------------------------------------------------------------------
subroutine rw_so2( t, nwet, rh, iratns, SAI_present, gw )

! Input/output variables:
real   , intent(in)  :: t      ! temperature (C)
integer, intent(in)  :: nwet   ! wetness indicator; nwet=0 -> dry; nwet=1 -> wet; nwet=9 -> snow
real   , intent(in)  :: rh     ! relative humidity (%)
integer, intent(in)  :: iratns ! index for NH3/SO2 ratio;
                               ! iratns = 1: low NH3/SO2
                               ! iratns = 2: high NH3/SO2
                               ! iratns = 3: very low NH3/SO2
logical, intent(in)  :: SAI_present
real   , intent(out) :: gw     ! external leaf conductance (m/s)

! Variables from module:
! SAI_present: vegetation is present

! Local variables:
real                 :: rw     ! external leaf resistance (s/m)

! Check if vegetation present:
if (SAI_present) then

   if (nwet .eq. 0) then
        !--------------------------
        ! dry surface
        !--------------------------
           ! T > -1 C
           if (t .gt. -1.0) then
              if (rh .lt. 81.3) then
                 rw = 25000*exp(-0.0693*rh)
              else
                 rw = 0.58e12*exp(-0.278*rh) + 10.
              endif
           else
              ! -5 C < T <= -1 C
              if (t .gt. -5.0) then
                 rw=200
              else
                 ! T <= -5 C
                 rw=500
              endif
           endif
     else
        !--------------------------
        ! wet surface
        !--------------------------
           rw = 10. !see Table 5, Erisman et al, 1994 Atm. Environment, 0 is impl. as 10
     endif

    ! very low NH3/SO2 ratio:
    if (iratns == iratns_very_low ) rw = rw+50.

    ! Conductance:
    gw = 1./rw
else
   ! no vegetation:
   gw = 0.0
endif

end subroutine rw_so2

!-------------------------------------------------------------------
! rw_nh3_sutton: compute external leaf conductance for NH3,
!                  following Sutton & Fowler, 1993
!-------------------------------------------------------------------
subroutine rw_nh3_sutton(tsurf,rh,SAI_present,gw)

! Input/output variables:
real   , intent(in)  :: tsurf     ! surface temperature (C)
real   , intent(in)  :: rh        ! relative humidity (%)
logical, intent(in)  :: SAI_present
real   , intent(out) :: gw        ! external leaf conductance (m/s)

! Variables from module:
! SAI_present: vegetation is present

! Local variables:
real                 :: rw                ! external leaf resistance (s/m)
real                 :: sai_grass_haarweg ! surface area index at experimental site Haarweg

! Fix SAI_grass at value valid for Haarweg data for which gamma_w parametrization is derived
sai_grass_haarweg = 3.5

!                  100 - rh
!  rw = 2.0 * exp(----------)
!                    12

if (SAI_present) then

   ! External resistance according to Sutton & Fowler, 1993
   rw = 2.0 * exp((100.0 - rh)/12.0)
   rw = sai_grass_haarweg * rw

   ! Frozen soil (from Depac v1):
   if (tsurf .lt. 0) rw = 200

   ! Conductance:
   gw = 1./rw
else
   ! no vegetation:
   gw = 0.0
endif

end subroutine rw_nh3_sutton

!-------------------------------------------------------------------
! rw_constant: compute constant external leaf conductance
!-------------------------------------------------------------------
subroutine rw_constant(rw_val,SAI_present,gw)

! Input/output variables:
real   , intent(in)  :: rw_val ! constant value of Rw
logical, intent(in)  :: SAI_present
real   , intent(out) :: gw     ! wernal leaf conductance (m/s)

! Variables from module:
! SAI_present: vegetation is present

! Compute conductance:
if (SAI_present .and. .not.missing(rw_val)) then
   gw = 1./rw_val
else
   gw = 0.
endif

end subroutine rw_constant

!-------------------------------------------------------------------
! rc_gstom: compute stomatal conductance
!-------------------------------------------------------------------
subroutine rc_gstom( icmp, compnam, lu, LAI_present, lai, &
                       glrad, sinphi, t, rh, henry, react, gstom, &
                        p, smi, calc_stom_o3flux )

! input/output variables:
integer,           intent(in) :: icmp          ! component index
character(len=*),  intent(in) :: compnam       ! component name
integer,           intent(in) :: lu            ! land use type , lu = 1,...,nlu
logical,           intent(in) :: LAI_present
real,              intent(in) :: lai           ! one-sided leaf area index
real(field_r),     intent(in) :: glrad         ! global radiation (W/m2)
real,              intent(in) :: sinphi        ! sin of solar elevation angle
real,              intent(in) :: t             ! temperature (C)
real,              intent(in) :: rh            ! relative humidity (%)
real,              intent(in) :: henry             ! effective Henry's law constant
real,              intent(in) :: react             ! reactivtiy
real,              intent(out):: gstom         ! stomatal conductance (m/s)
real, optional,    intent(in) :: p             ! pressure (Pa)
real, optional,    intent(in) :: smi           ! soil moisture index      
logical, optional, intent(in) :: calc_stom_o3flux  ! calculate stomatal ozone flux??

! variables from module
! LAI_present: vegetation is present
! dwat: diffusion coefficient of water vapour
! dO3 : diffusion coefficient of ozone

! Local variables
real :: vpd  ! vapour pressure deficit (kPa)

select case(trim(compnam))

!case('HNO3','N2O5','NO3','H2O2') this routine is not called for HNO3

case('NO','CO')
   ! for NO stomatal uptake is neglected:
   gstom = 0.0

case('NO2','O3','SO2','NH3')

   ! if vegetation present:
   if (LAI_present) then

      if (glrad .gt. 0.0) then
         call rc_get_vpd(t,rh,vpd)
         call rc_gstom_emb( lu, glrad, t, vpd, LAI_present, lai, sinphi, gstom, p, &
                            smi=smi, calc_stom_o3flux=calc_stom_o3flux )
         gstom = gstom*diffc(icmp)/dO3       ! Gstom of Emberson is derived for ozone 
      else
         gstom = 0.0
      endif
   else
      ! no vegetation; zero conductance (infinite resistance):
      gstom = 0.0
   endif

case default
  if ( compnam(1:3) == 'VBS' ) then
    ! if vegetation present:
    if (LAI_present) then
      if (glrad .gt. 0.0) then
         call rc_get_vpd(t,rh,vpd)
         call rc_gstom_emb( lu, glrad, t, vpd, LAI_present, lai, sinphi, gstom, p, &
                            smi=smi )
         gstom = gstom*diffc(icmp)/dO3        ! Gstom of Emberson is derived for ozone 
         gstom = gstom * (1e-5*henry + react) ! Wesely '89, eq. 7 
      else
         gstom = 0.0
      endif
    else
      ! no vegetation; zero conductance (infinite resistance):
      gstom = 0.0
    endif
  else
    print *, 'error in subroutine rc_gstom '
    print *, 'component ',trim(compnam),' not supported'
    stop
  end if
end select

end subroutine rc_gstom

!-------------------------------------------------------------------
! rc_gstom_emb: stomatal conductance according to Emberson
!-------------------------------------------------------------------
subroutine rc_gstom_emb( lu, glrad, T, vpd, LAI_present, lai, sinp, &
                            Gsto, &
                            p, smi, calc_stom_o3flux )

! History
!   Original code from Lotos-Euros, TNO, M. Schaap (?)
!   2009-08, M.C. van Zanten, Rivm
!     Updated and extended.
!   2009-09, Arjo Segers, TNO
!     Limitted temperature influence to range to avoid
!     floating point exceptions.
!
! Method
!
!   Code based on Emberson et al, 2000, Env. Poll., 403-413
!   Notation conform Unified EMEP Model Description Part 1, ch 8
!
!   In the calculation of F_light the modification of L. Zhang 2001, AE to the PARshade and PARsun
!   parametrizations of Norman 1982 are applied
!
!   F_phen and F_SWP are set to 1
!
!   Land use types DEPAC versus Emberson (Table 5.1, EMEP model description)
!   DEPAC                     Emberson
!     1 = grass                 GR = grassland
!     2 = arable land           TC = temperate crops ( LAI according to RC = rootcrops)
!     3 = permanent crops       TC = temperate crops ( LAI according to RC = rootcrops)
!     4 = coniferous forest     CF = tempareate/boreal coniferous forest
!     5 = deciduous forest      DF = temperate/boreal deciduous forest
!     6 = water                 W  = water
!     7 = urban                 U  = urban
!     8 = other                 GR = grassland
!     9 = desert                DE = desert
!

! Emberson specific declarations

! Input/output variables:
integer,           intent(in) :: lu                ! land use type, lu = 1,...,nlu
real(field_r),     intent(in) :: glrad             ! global radiation (W/m2)
real,              intent(in) :: T                 ! temperature (C)
real,              intent(in) :: vpd               ! vapour pressure deficit (kPa)
logical,           intent(in) :: LAI_present
real,              intent(in) :: lai               ! one-sided leaf area index
real,              intent(in) :: sinp              ! sin of solar elevation angle
real,              intent(out):: Gsto              ! stomatal conductance (m/s)
real, optional,    intent(in) :: p                 ! pressure (Pa)
real, optional,    intent(in) :: smi               ! soil moisture index
logical, optional, intent(in) :: calc_stom_o3flux  ! calculate stomatal ozone flux??

! Local variables:
real                          ::  F_light, F_phen, F_temp, F_vpd, F_swp, bt, F_env
real                          ::  PARdir, PARdiff, PARshade, PARsun, LAIsun, LAIshade, sinphi
real                          ::  pres
real, parameter               ::  p_sealevel = 1.01325e5    ! Pa


! Check whether vegetation is present:
if (LAI_present) then

   ! calculation of correction factors for stomatal conductance
   if (sinp .le. 0.0) then  ! due to the coarse calc. of sinphi, the value is sometimes less than 0.
       sinphi = 0.0001
   else
       sinphi = sinp
   end if
   
   ! ratio between actual and sea-level pressure is used
   ! to correct for height in the computation of PAR ;
   ! should not exceed sea-level pressure therefore ...
   if ( present(p) ) then
       pres = min( p, p_sealevel )
   else 
       pres = p_sealevel
   endif
   
   ! Direct and diffuse PAR, Photoactive (=visible) radiation:
   call par_dir_diff( glrad, sinphi, pres, p_sealevel, PARdir, PARdiff )

   ! PAR for shaded leaves (canopy averaged):
   PARshade = PARdiff*exp(-0.5*LAI**0.7)+0.07*PARdir*(1.1-0.1*LAI)*exp(-sinphi) ! Norman, 1982
   if (glrad .gt. 200 .and. LAI .gt. 2.5) then
      PARshade = PARdiff*exp(-0.5*LAI**0.8)+0.07*PARdir*(1.1-0.1*LAI)*exp(-sinphi) ! Zhang et al., 2001
   end if

   ! PAR for sunlit leaves (canopy averaged):
   ! alpha -> mean angle between leaves and the sun is fixed at 60 deg -> i.e. cos alpha = 0.5
   PARsun = PARdir*0.5/sinphi + PARshade ! Norman, 1982
   if (glrad .gt. 200 .and. LAI .gt. 2.5) then
       PARsun = PARdir**0.8*0.5/sinphi + PARshade ! Zhang et al., 2001
   end if

   ! leaf area index for sunlit and shaded leaves:
   if (sinphi .gt. 0) then
     LAIsun = 2*sinphi*(1-exp(-0.5*LAI/sinphi ))
     LAIshade = LAI -LAIsun
   else
     LAIsun = 0
     LAIshade = LAI
   end if

   ! correction factor for radiation (Emberson):
   if ( present(calc_stom_o3flux) ) then
     ! for calculation of stomatal ozone flux, this is different.
     if ( calc_stom_o3flux ) then
       F_light = 1 - exp(-1.*alpha(lu)*PARsun) 
     else
       F_light = (LAIsun*(1 - exp(-1.*alpha(lu)*PARsun)) + LAIshade*(1 - exp(-1.*alpha(lu)*PARshade)))/LAI 
     end if
   else 
     F_light = (LAIsun*(1 - exp(-1.*alpha(lu)*PARsun)) + LAIshade*(1 - exp(-1.*alpha(lu)*PARshade)))/LAI 
   end if
   
   F_light = max(F_light,F_min(lu))

   !  temperature influence; only non-zero within range [Tmin,Tmax]:
   if ( (Tmin(lu) < T) .and. (T < Tmax(lu)) ) then
     BT = (Tmax(lu)-Topt(lu))/(Topt(lu)-Tmin(lu))
     F_temp = ((T-Tmin(lu))/(Topt(lu)-Tmin(lu))) * ((Tmax(lu)-T)/(Tmax(lu)-Topt(lu)))**BT
   else
     F_temp = 0.0
   end if
   F_temp = max(F_temp,F_min(lu))

   ! vapour pressure deficit influence
   F_vpd = min(1.,((1.-F_min(lu))*(vpd_min(lu)-vpd)/(vpd_min(lu)-vpd_max(lu)) + F_min(lu) ))
   F_vpd = max(F_vpd,F_min(lu))

   ! influence waterpotential PSI > PSI(0).
   if (present(smi)) then
     F_swp = min( 1.0, 2.0*smi)
     !F_swp = 1.0
   else
     F_swp = 1.
   end if
      
   ! influence of phenology on stom. conductance
   ! ignored for now in DEPAC since influence of F_phen on lu classes in use is negligible.
   ! When other EMEP classes (e.g. med. broadleaf) are used f_phen might be too important to ignore
   F_phen = 1.

   ! evaluate total stomatal conductance
   F_env = F_temp*F_vpd*F_swp
   F_env = max(F_env,F_min(lu))
   gsto = G_max(lu) * F_light * F_phen * F_env

   ! gstom expressed per m2 leafarea;
   ! this is converted with LAI to m2 surface.
   Gsto = lai*gsto    ! in m/s

else
  GSto = 0.0
endif

end subroutine rc_gstom_emb

!-------------------------------------------------------------------
! par_dir_diff
!-------------------------------------------------------------------
subroutine par_dir_diff(glrad,sinphi,pres,pres_0,par_dir,par_diff)
!
!     Weiss, A., Norman, J.M. (1985) Partitioning solar radiation into direct and
!     diffuse, visible and near-infrared components. Agric. Forest Meteorol.
!     34, 205-213.
!
!     From a subroutine obtained from Leiming Zhang (via Willem Asman),
!     Meteorological Service of Canada
!     e-mail: leiming.zhang@ec.gc.ca
!
!     Leiming uses solar irradiance. This should be equal to global radiation and
!     Willem Asman set it to global radiation

      real(field_r), intent(in)  :: glrad             ! global radiation (W m-2)
      real,          intent(in)  :: sinphi            ! sine of the solar elevation
      real,          intent(in)  :: pres              ! actual pressure (to correct for height) (Pa)
      real,          intent(in)  :: pres_0            ! pressure at sea level (Pa)
      real,          intent(out) :: par_dir           ! PAR direct : visible (photoactive) direct beam radiation (W m-2)
      real,          intent(out) :: par_diff          ! PAR diffuse: visible (photoactive) diffuse radiation (W m-2)

!     fn              = near-infrared direct beam fraction (dimensionless)
!     fv              = PAR direct beam fraction (dimensionless)
!     ratio           = ratio measured to potential solar radiation (dimensionless)
!     rdm             = potential direct beam near-infrared radiation (W m-2); "potential" means clear-sky
!     rdn             = potential diffuse near-infrared radiation (W m-2)
!     rdu             = visible (PAR) direct beam radiation (W m-2)
!     rdv             = potential visible (PAR) diffuse radiation (W m-2)
!     rn              = near-infrared radiation (W m-2)
!     rv              = visible radiation (W m-2)
!     ww              = water absorption in the near infrared for 10 mm of precipitable water

      real :: rdu,rdv,ww,rdm,rdn,rv,rn,ratio,sv,fv

!     Calculate visible (PAR) direct beam radiation
!     600 W m-2 represents average amount of PAR (400-700 nm wavelength)
!     at the top of the atmosphere; this is roughly 0.45*solar constant (solar constant=1320 Wm-2)
      rdu=600.*exp(-0.185*(pres/pres_0)/sinphi)*sinphi

!     Calculate potential visible diffuse radiation
      rdv=0.4*(600.- rdu)*sinphi

!     Calculate the water absorption in the-near infrared
      ww=1320*10**( -1.195+0.4459*log10(1./sinphi)-0.0345*(log10(1./sinphi))**2 )

!     Calculate potential direct beam near-infrared radiation
      rdm=(720.*exp(-0.06*(pres/pres_0)/sinphi)-ww)*sinphi     !720 = solar constant - 600

!     Calculate potential diffuse near-infrared radiation
      rdn=0.6*(720-rdm-ww)*sinphi

      ! Compute visible and near-infrared radiation
      rv=max(0.1,rdu+rdv)
      rn=max(0.01,rdm+rdn)

      ! Compute ratio between input global radiation and total radiation computed here
      ratio=min(0.9,glrad/(rv+rn))

!     Calculate total visible radiation
      sv=ratio*rv

!     Calculate fraction of PAR in the direct beam
      fv=min(0.99, (0.9-ratio)/0.7)            ! help variable
      fv=max(0.01,rdu/rv*(1.0-fv**0.6667))     ! fraction of PAR in the direct beam

      ! Compute direct and diffuse parts of PAR
      par_dir=fv*sv
      par_diff=sv-par_dir

end subroutine par_dir_diff

!-------------------------------------------------------------------
! rc_get_vpd: get vapour pressure deficit (kPa)
!-------------------------------------------------------------------
subroutine rc_get_vpd(temp, relh,vpd)

implicit none

! Input/output variables:
real    , intent(in)  :: temp   ! temperature (C)
real    , intent(in)  :: relh   ! relative humidity (%)
real    , intent(out) :: vpd    ! vapour pressure deficit (kPa)

! Local variables:
real :: esat

! fit parameters:
real, parameter :: a1 = 6.113718e-1
real, parameter :: a2 = 4.43839e-2
real, parameter :: a3 = 1.39817e-3
real, parameter :: a4 = 2.9295e-5
real, parameter :: a5 = 2.16e-7
real, parameter :: a6 = 3.0e-9

! esat is saturation vapour pressure (kPa) at temp(C) following monteith(1973)
esat = a1 + a2*temp + a3*temp**2 + a4*temp**3 + a5*temp**4 + a6*temp**5
vpd  = esat * (1-relh/100)

end subroutine rc_get_vpd

!-------------------------------------------------------------------
! rc_mes: compute mesophyll conductance
!-------------------------------------------------------------------
subroutine rc_mes( icmp, compnam, henry, react, gmes)

! Input/output variables:
integer, intent(in)           :: icmp    ! component index 
character(len=*), intent(in)  :: compnam ! component name
real, intent(in)              :: henry   ! effective henry's law constant
real, intent(in)              :: react   ! reactivity
real, intent(out)             :: gmes    ! mesophyll conductance (m s-1)

real                          :: rmes    ! mesophyll resistance (s m-1)

if ( compnam(1:3) == 'VBS' ) then
    rmes = 1/(henry/3000.+100.*react)  ! Wesely '89, eq. 6
    gmes = 1/rmes
else
    gmes = 0.0
end if

end subroutine rc_mes

!-------------------------------------------------------------------
! rc_snow: compute total canopy resistance in case of snow
!-------------------------------------------------------------------
subroutine rc_snow(ipar_snow,t,rc_tot)

! Input/output variables:
integer, intent(in)  :: ipar_snow  ! choice of parameterisation in case of snow
                                   ! ipar_snow = 1 : constant parameterisation (=rssnow)
                                   ! ipar_snow = 2 : temperature dependent parameterisation
real   , intent(in)  :: t          ! temperature (C)
real   , intent(out) :: rc_tot     ! total canopy resistance Rc (s/m)

! Local variables:
real, parameter :: rssnow = 2000.  ! constant resistance in case of snow and ipar_snow = 1

! Choose parameterisation with constant or temperature dependent parameterisation:
if (ipar_snow .eq. 1) then
   rc_tot = rssnow
elseif (ipar_snow .eq. 2) then
   if (t .lt. -1.) then
      rc_tot = 500.
   elseif (t .gt.  1.) then
      rc_tot = 70.
   else ! (t .ge. -1. .and. t .le. 1.)
      rc_tot = 70.*(2.-t)
   endif
else
   write(*,*) ' programming error in subroutine rc_snow'
   write(*,*) ' unknown value of ipar_snow: ',ipar_snow
   stop
endif

end subroutine rc_snow

!-------------------------------------------------------------------
! rc_gsoil_eff: compute effective soil conductance
!-------------------------------------------------------------------
subroutine rc_gsoil_eff( icmp, compnam, lu, sai, ust, nwet, t, &
                          henry, react, icmp_o3, icmp_so2, &
                          gsoil_eff, status )

! Input/output variables:
integer, intent(in)  :: icmp           ! component index
character(len=*), intent(in)  :: compnam            ! component name
integer, intent(in)  :: lu             ! land use type, lu = 1,...,nlu
real   , intent(in)  :: sai            ! surface area index
real   , intent(in)  :: ust            ! friction velocity (m/s)
integer, intent(in)  :: nwet           ! index for wetness
                                       ! nwet=0 -> dry; nwet=1 -> wet; nwet=9 -> snow
                                       ! N.B. this routine cannot be called with nwet = 9,
                                       ! nwet = 9 should be handled outside this routine.
real   , intent(in)  :: t              ! temperature (C)
real   , intent(in)  :: henry          ! effective Hnery's law constant (M atm-1)
real   , intent(in)  :: react          ! reactivity (-)
integer, intent(in)  :: icmp_o3        ! index O3
integer, intent(in)  :: icmp_so2       ! index SO2
real   , intent(out) :: gsoil_eff      ! effective soil conductance (m/s)
  integer, intent(out) :: status
  
  character(len=*), parameter   ::  rname = mname//'/rc_gsoil_eff'

! local variables:
real                 :: rinc           ! in canopy resistance  (s/m)
real                 :: rsoil_eff      ! effective soil resistance (s/m)

! Compute in canopy (in crop) resistance:
call rc_rinc(lu,sai,ust,rinc)

! Check for missing deposition path:
if (missing(rinc)) then
  rsoil_eff = -9999.
else
   ! Frozen soil (temperature below 0):
   if (t .lt. 0.0) then
      if (( compnam(1:3) == 'VBS')) then 
         ! Wesely '89, eq. 9
         rsoil_eff = 1/( henry/(1e5*rsoil_frozen(icmp_so2)) + react/rsoil_frozen(icmp_o3) ) + rinc  
      else if (missing(rsoil_frozen(icmp))) then
         rsoil_eff = -9999.
      else
         rsoil_eff = rsoil_frozen(icmp) + rinc
      endif
   else
      ! Non-frozen soil; dry:
      if (nwet .eq. 0) then
         if ( compnam(1:3) == 'VBS' ) then
            ! Wesely '89, eq. 9 
            rsoil_eff = 1/( henry/(1e5*rsoil(lu,icmp_so2)) + react/rsoil(lu,icmp_o3) ) + rinc  
         else if (missing(rsoil(lu,icmp))) then
            rsoil_eff = -9999.
         else
            rsoil_eff = rsoil(lu,icmp) + rinc
         endif

      ! Non-frozen soil; wet:
      elseif (nwet .eq. 1) then
         if ( compnam(1:3) == 'VBS' ) then
            ! Wesely '89, eq. 9
            rsoil_eff = 1/( henry/(1e5*rsoil_wet(icmp_so2)) + react/rsoil_wet(icmp_o3) ) + rinc
         else if (missing(rsoil_wet(icmp))) then
            rsoil_eff = -9999.
         else
            rsoil_eff = rsoil_wet(icmp) + rinc
         endif
      else
           write (gol,'("unsupported nwet value ",i0," should 0 or 1")') nwet; call goErr
           TRACEBACK; status=1; return
      endif
   endif
endif

! Compute conductance:
if (rsoil_eff .gt. 0.0) then
   gsoil_eff = 1./rsoil_eff
else
   gsoil_eff = 0.0
endif

  ! ok
  status = 0

end subroutine rc_gsoil_eff

!-------------------------------------------------------------------
! rc_rinc: compute in canopy (or in crop) resistance
! van Pul and Jacobs, 1993, BLM
!-------------------------------------------------------------------
subroutine rc_rinc(lu,sai,ust,rinc)

!use LE_Landuse_Data, only : b, h

! Input/output variables:
integer, intent(in)  :: lu          ! land use class, lu = 1, ..., nlu
real   , intent(in)  :: sai         ! surface area index
real   , intent(in)  :: ust         ! friction velocity (m/s)
real   , intent(out) :: rinc        ! in canopy resistance (s/m)

!! for in-crop resistance
!! b = empirical constant for computation of rinc (in canopy resistance) (= 14 m-1 or -999 if not applicable)
!! h = vegetation height (m)                  grass arabl  crops conif decid  water  urba  othr  desr  ice   wheat beech spruce semi
!real, dimension(nlu), parameter:: b=       (/ -999,  14,    14,   14,   14,    -999,  -999, -999, -999, -999, 14,   14,   14,    14/)
!real, dimension(nlu), parameter:: h=       (/ -999,  1,     1,    20,   20,    -999,  -999, -999, -999, -999, 1,    20,   20,    1/)  

! Compute Rinc only for arable land, perm. crops, forest; otherwise Rinc = 0:
if (b(lu) .gt. 0.0) then

  ! Check for u* > 0 (otherwise denominator = 0):
  if (ust .gt. 0.0) then
     rinc = b(lu)*h(lu)*sai/ust
  else
     rinc = 1000.0
  endif
else
   if (lu .eq. ilu_grass .or. lu .eq. ilu_other ) then
      rinc = -999. ! no deposition path for grass, other, and semi-natural
   else
      rinc = 0.0   ! no in-canopy resistance
   endif
endif

end subroutine rc_rinc

!-------------------------------------------------------------------
! rc_rctot: compute total canopy (or surface) resistance Rc
!-------------------------------------------------------------------
subroutine rc_rctot(gstom,gsoil_eff,gw,gmes,gc_tot,rc_tot)

! Input/output variables:
real, intent(in)  :: gstom        ! stomatal conductance (s/m)
real, intent(in)  :: gsoil_eff    ! effective soil conductance (s/m)
real, intent(in)  :: gw           ! external leaf conductance (s/m)
real, intent(in)  :: gmes         ! mesophyll conductance (s/m)
real, intent(out) :: gc_tot       ! total canopy conductance (m/s)
real, intent(out) :: rc_tot       ! total canopy resistance Rc (s/m)

! Total conductance:
gc_tot = gstom + gsoil_eff + gw + gmes 

! Total resistance (note: gw can be negative, but no total emission allowed here):
if (gc_tot .le. 0.0 .or. gw .lt. 0.0) then
   rc_tot = -9999.
else
   rc_tot = 1./gc_tot
endif

end subroutine rc_rctot

  !-------------------------------------------------------------------
  ! rc_comp_point: calculate compensation points (stomata, external leaf)
  !-------------------------------------------------------------------
  ! Calculate ccomp, i.e. compensation point for NH3, for different deposition path ways
  ! (external leaf, stomata, soil), according to Roy Wichink Kruit article submitted 2009 Atm. Env.
  !
  !               2.75e15     -1.04e4
  ! ccomp = gamma ------- exp(-------)
  !                 Tk           Tk
  ! with
  ! gamma : [NH4+]/[H+] ratio in apoplast (or leaf)
  ! Tk    : temperature (K)
  !
  ! The [NH4+]/[H+] ratio gamma depends on
  ! 1. for stomata
  !    the average concentration over a previous period:
  !    gamma_stom = gamma_stom_c_fac * c_ave_prev_nh3
  !
  ! 2. for external leaf
  !    actual atmospheric concentration at 4 m. and surface temperature:
  !    gamma_w = -850.+1840.*catm*exp(-0.11*Tsurf) ; Tsurf = surface temperature in C
  !   (multivariate fit on Haarweg data Roy Wichink Kruit 2009)
  !
  ! 3. for soil
  !    unkown yet, except for water surfaces
  !
  ! For components other than NH3, compensation points are 0.

subroutine rc_comp_point( compnam,lu,day_of_year,t,gw,gstom,gsoil_eff,gc_tot,&
                          LAI_present, SAI_present, &
                          ccomp_tot, status, &
                          catm,c_ave_prev_nh3,c_ave_prev_so2,gamma_soil_water,tsea,cw,cstom,csoil )

    ! Input/output variables:
    character(len=*), intent(in)  :: compnam            ! component name
                                                        ! 'HNO3','NO','NO2','O3','SO2','NH3'
    integer, intent(in)           :: lu                 ! land use type, lu = 1,...,9
    integer, intent (in)          :: day_of_year        ! day of year
    real, intent(in)              :: t                  ! temperature (C)
    real, intent(in)              :: gw                 ! external leaf conductance (m/s)
    real, intent(in)              :: gstom              ! stomatal conductance (m/s)
    real, intent(in)              :: gsoil_eff          ! effective soil conductance (m/s)
    real, intent(in)              :: gc_tot             ! total canopy conductance (m/s)
    logical, intent(in)           :: LAI_present
    logical, intent(in)           :: SAI_present
    real, intent(out)             :: ccomp_tot          ! total compensation point (weighted average of
                                                        ! separate compensation points) (ug/m3)
    integer, intent(out)          :: status
    
    real(field_r), optional,  intent(in)   :: catm               ! actual atmospheric concentration (ug/m3)
    real, optional,  intent(in)   :: c_ave_prev_nh3     ! air concentration averaged over a previous
                                                        ! period (e.g. previous year or month) (ug/m3)
    real, optional,  intent(in)   :: c_ave_prev_so2     ! air concentration averaged over a previous
                                                        ! period (e.g. previous year or month) (ug/m3)
    real, optional,  intent(in)   :: gamma_soil_water   ! gammawater from gammawater.nc file
    real, optional,  intent(in)   :: tsea               ! sea surface temperature (K)
    real, optional,  intent(out)  :: cw                 ! external leaf surface compensation point (ug/m3)
    real, optional,  intent(out)  :: cstom              ! stomatal compensation point (ug/m3)
    real, optional,  intent(out)  :: csoil              ! soil compensation point (ug/m3)

    ! constants:
    character(len=*), parameter   ::  rname = mname//'/rc_comp_point'

    !
    ! FOR FUTURE IMPLEMENTATION:
    ! These arrays should be input !
    ! for current parametrization gamma_stom_c_fac is independent of lu
    !
    !                                            grass arable  perm. conif. decid.  water  urban  other  desert     ice    sav   trf    wai    med   semi-   wheat  beech spruce  
    !                                                    land  crops forest forest                                                            scrub  natu
    !                                               1      2      3      4      5      6      7      8       9       10     11    12     13     14     15      16    17     18
    ! 
    ! factor in linear relation between gamma_stom 
    ! and NH3 air concentration; 
    ! gamma_stom = [NH4+]/[H+] ratio in apoplast
    real, parameter :: gamma_stom_c_fac(nlu_ozone_specials) = &
                                                (/  362,   362,   362,   362,   362,  -999,  -999,   362,   -999,  -999,   362,   362,  -999,  -999,   362,   362,   362,   362 /)
    ! factor in linear relation between gamma_soil
    ! and NH3 air concentration; 
    ! gamma_soil = [NH4+]/[H+] ratio in soil
    ! NOTE: for each landuse only one value may be given in either gamma_soil_c_fac(nlu) or gamma_soil(nlu)
    real, parameter :: gamma_soil_c_fac(nlu_ozone_specials)   = &
                                                 (/ -999,  -999,  -999,  -999,  -999,  -999,  -999,  -999,  -999,  -999,  -999,  -999,  -999,  -999,  -999,  -999,  -999,  -999 /) ! for future use
    real, parameter :: gamma_soil_default(nlu_ozone_specials) = &
                                                 (/ -999,  -999,  -999,  -999,  -999,   430,  -999,  -999,  -999,  -999,  -999,  -999,   430,  -999,  -999,  -999,  -999,  -999 /)

    ! Variables from module:
    ! gamma_stom_c_fac: factor in linear relation between gamma_stom and NH3 air concentration.
    ! LAI_present or SAI_present: vegetation is present

    ! Local variables:
    real :: gamma_stom ! [NH4+]/[H+] ratio in apoplast
    real :: gamma_soil ! [NH4+]/[H+] ratio in soil
    real :: gamma_w    ! [NH4+]/[H+] ratio in external leaf surface water
    real :: tk         ! temperature (K)
    real :: tfac       ! temperature factor = (2.75e15/tk)*exp(-1.04e4/tk)
    real :: co_dep_fac ! co-deposition factor 

    select case(trim(compnam))

      !case('HNO3','N2O5','NO3','H2O2')
      !  this routine is not called for HNO3

      case('NO','CO','NO2','O3','SO2')

        ! no compensation points:
        ccomp_tot  = 0.0

      case('NH3')

        ! Temperature factor:
        ! parametrized temperature of surface water including yearly cycle
        ! parametrization based on NL Waterbase data for 2003-2008, for ~25 locations
        ! for which NH4+, PH and temperature measurements are present
        if (lu .eq. ilu_water_sea .or. lu .eq. ilu_water_inland ) then ! water
          !tk = 286.2 + 8.3*sin(day_of_year -113.5)
          tk = tsea ! for LOTOS-EUROS sea surface temperature is available from ECMWF meteo @ Mar11 RWK
        else
           tk = t + 273.15
        endif
        tfac = (2.75e15/tk)*exp(-1.04e4/tk)

        ! Stomatal compensation point:
        if (LAI_present .and. present(c_ave_prev_nh3) ) then
          ! gamma_stom ([NH4+]/[H+] ratio in apoplast) is linearly dependent on an
          ! averaged air concentration in a previous period (stored in soil and leaves):
          gamma_stom = gamma_stom_c_fac(lu)*c_ave_prev_nh3*4.7*exp(-0.071*t)
          ! calculate stomatal compensation point for NH3 in ug/m3:
          cstom = max(0.0,gamma_stom*tfac)
        else
          ! No concentration in previous period or no gamma-c factor:
          cstom = 0.0
        endif

        ! External leaf gamma depends on atmospheric concentration and
        ! surface temperature (assumed to hold for all land use types with vegetation):
        if (SAI_present .and. present(catm) .and. present(c_ave_prev_nh3) .and. present(c_ave_prev_so2) ) then
          gamma_w = -850.+1840.*catm*exp(-0.11*t)
          ! correction gamma_w for co deposition 
          ! xxx documentation to be added
          ! gamma(with SNratio) = [1.12-1.32*SNratio(molair)]  * gamma_original   
          ! where SNratio(molar) = (CSO2longterm/64)/(CNH3longterm/17))
          !   where CSO2longterm and CNH3longterm in ug m-3.
          if (c_ave_prev_nh3 .gt. 0. .and. c_ave_prev_so2 .gt. 0.) then
          co_dep_fac = 1.12 - 1.32 * ((c_ave_prev_so2/64.)/(c_ave_prev_nh3/17.))
          else
            co_dep_fac = 0.0
          endif
          co_dep_fac = max(0.0,co_dep_fac)   
          gamma_w = co_dep_fac * gamma_w 
          cw      = max(0.0,gamma_w*tfac)
        elseif (SAI_present .and. present(catm)) then
          gamma_w = -850.+1840.*catm*exp(-0.11*t)
          cw      = max(0.0,gamma_w*tfac)
        else
          cw = 0.0
        endif

        ! Soil compensation point;
        ! only if soil could take up tracers:
        if ( gamma_soil_c_fac(lu) > 0 ) then
          if (present(c_ave_prev_nh3)) then                                             
            ! gamma_soil ([NH4+]/[H+] ratio in soil) is linearly dependent on an            
            ! averaged air concentration in a previous period (analogue to stomatal gamma;  
            ! can be used to test, but is not effectively used as all other landuses are    
            ! set to -999):                                                                 
            gamma_soil = gamma_soil_c_fac(lu) * c_ave_prev_nh3                              
          else                                                                              
            gamma_soil = 0                                                                  
          end if
        elseif ( gamma_soil_default(lu) > 0 ) then
          ! which landuse ?                                                                 
          if ( lu == ilu_water_sea .or. lu .eq. ilu_water_inland ) then                                                      
            ! gamma_soil for water is determined to be 430 based on Waterbase data          
            ! close to the Dutch coastline, for LOTOS-EUROS gamma values are estimated      
            ! based on annual dry deposition maps with compensation points set to 0.        
            ! See Wichink Kruit et al. (2010) for derivation of gamma_water values.         
            ! These gamma_water values are read from a file (gamma_soil_water).           
            if (present(gamma_soil_water)) then                                           
              gamma_soil = gamma_soil_water                                               
            else                                                                            
              gamma_soil = gamma_soil_default(lu)                                             
            end if                                                                          
          else                                                                              
            gamma_soil = gamma_soil_default(lu)                                             
          end if  ! water or land                                                           
        else
          ! no soil compensation point ...
          gamma_soil = 0.0
        end if  ! soil could take up gas
        ! calculate soil compensation point for NH3 in ug/m3:                             
        csoil = gamma_soil * tfac                                                         

        ! Total compensation point is weighted average of separate compensation points:
        if (gc_tot .gt. 0.0 ) then
          ccomp_tot = (gw/gc_tot)*cw + (gstom/gc_tot)*cstom + (gsoil_eff/gc_tot)*csoil
        else
          ccomp_tot = 0.0
        endif
        
      case default
        if ( compnam(1:3) == 'VBS' ) then
          ! no compensation points:
          ccomp_tot  = 0.0
        else
          write (gol,'("unsupported component `",a,"`")') trim(compnam); call goErr
          TRACEBACK; status=1; return
        end if
    end select
    
    ! ok
    status = 0

end subroutine rc_comp_point

!-------------------------------------------------------------------
! rc_comp_point_rc_eff: calculate the effective resistance Rc
! based on one or more compensation points
!-------------------------------------------------------------------
!
! old name: NH3rc (see depac v3.6 is based on Avero workshop Marc Sutton. p. 173.
! Sutton 1998 AE 473-480)
!
! Documentation by Ferd Sauter, 2008; see also documentation block in header of this module.
! FS 2009-01-29: variable names made consistent with rest of DEPAC module
! FS 2009-03-04: use total compensation point
!
! C: with total compensation point   ! D: approximation of C
!                                    !    with classical approach
!  zr --------- Catm                 !  zr --------- Catm    
!         |                          !         |       
!         Ra                         !         Ra      
!         |                          !         |       
!         Rb                         !         Rb      
!         |                          !         |       
!  z0 --------- Cc                   !  z0 --------- Cc
!         |                          !         |              
!        Rc                          !        Rc_eff          
!         |                          !         |              
!     --------- Ccomp_tot            !     --------- C=0
!
!
! The effective Rc is defined such that instead of using
!
!   F = -vd*[Catm - Ccomp_tot]                                    (1)
!
! we can use the 'normal' flux formula
!
!   F = -vd'*Catm,                                                (2)
!
! with vd' = 1/(Ra + Rb + Rc')                                    (3)
!
! and Rc' the effective Rc (rc_eff).
!                                                (Catm - Ccomp_tot)
! vd'*Catm = vd*(Catm - Ccomp_tot) <=> vd' = vd* ------------------
!                                                      Catm
!
!                                        (Catm - Ccomp_tot)
! 1/(Ra + Rb + Rc') = (1/Ra + Rb + Rc) * ------------------
!                                              Catm
!
!                                          Catm
! (Ra + Rb + Rc') = (Ra + Rb + Rc) * ------------------
!                                     (Catm - Ccomp_tot)
!
!                              Catm
! Rc' = (Ra + Rb + Rc) * ------------------ - Ra - Rb
!                        (Catm - Ccomp_tot)
!
!                        Catm                           Catm
! Rc' = (Ra + Rb) [------------------ - 1 ] + Rc * ------------------
!                  (Catm - Ccomp_tot)              (Catm - Ccomp_tot)
!
! Rc' = [(Ra + Rb)*Ccomp_tot + Rc*Catm ] / (Catm - Ccomp_tot)
!
! This is not the most efficient way to do this;
! in the current LE version, (1) is used directly in the calling routine
! and this routine is not used anymore.
!
! -------------------------------------------------------------------------------------------
! In fact it is the question if this correct; note the difference in differential equations
!
!   dCatm                              !   dCatm
! H ----- = F = -vd'*Catm,             ! H ----- = F = -vd*[Catm - Ccomp_tot],
!    dt                                !    dt
!
! where in the left colum vd' is a function of Catm and not a constant anymore.
!
! Another problem is that if Catm -> 0, we end up with an infinite value of the exchange
! velocity vd.
! -------------------------------------------------------------------------------------------

subroutine rc_comp_point_rc_eff( ccomp_tot, catm, ra, rb, rc_tot, rc_eff, status )

implicit none

! Input/output variables:
real, intent(in)          :: ccomp_tot  ! total compensation point (weighed average of separate compensation points) (ug/m3)
real(field_r), intent(in) :: catm       ! atmospheric concentration (ug/m3)
real, intent(in)          :: ra         ! aerodynamic resistance (s/m)
real, intent(in)          :: rb         ! boundary layer resistance (s/m)
real, intent(in)          :: rc_tot     ! total canopy resistance (s/m)
real, intent(out)         :: rc_eff     ! effective total canopy resistance (s/m)
integer, intent(out)      :: status

! constants:
character(len=*), parameter   ::  rname = mname//'/rc_comp_point_rc_eff'

! Compute effective resistance:
if ( ccomp_tot == 0.0 ) then
  ! trace with no compensiation point ( or compensation point equal to zero)
  rc_eff = rc_tot
    
else if ( ccomp_tot > 0 .and. ( abs( catm - ccomp_tot ) < 1.e-8 ) ) then
  ! surface concentration (almoast) equal to atmospheric concentration
  ! no exchange between surface and atmosphere, infinite RC --> vd=0
  rc_eff = 9999999999.
  
else if ( ccomp_tot > 0 ) then
  ! compensation point available, calculate effective restisance
  rc_eff = ((ra + rb)*ccomp_tot + rc_tot*catm)/(catm-ccomp_tot)
   
else
   rc_eff = -999.
   write (gol,*) 'This should not be possible:'; call goErr
   write (gol,*) 'rc_tot     = ', rc_tot; call goErr
   write (gol,*) 'ccompt_tot = ', ccomp_tot; call goErr
   write (gol,*) 'catm       = ', catm; call goErr
   TRACEBACK; status=1; return
end if

! ok
status = 0

end subroutine rc_comp_point_rc_eff

!-------------------------------------------------------------------
! missing: check for data that correspond with a missing deposition path
!          this data is represented by -999
!-------------------------------------------------------------------

logical function missing(x)

real, intent(in) :: x

! bandwidth for checking (in)equalities of floats
real, parameter :: EPS = 1.0e-5

missing = (abs(x + 999.) .le. EPS)

end function missing


  !**********************************************************************************
  !
  ! NH3/SO2 ratio
  !
  !**********************************************************************************

  !  -  NH3/SO2 ratio
  !     ratns=1 low ratio NH3/SO2
  !     ratns=2 high ratio NH3/SO2
  !     ratns=3 very low ratio NH3/SO2
  !    see calling routine:
  !     [SO2] < 0.5 ppb then
  !        iratns = 1 ! low
  !     else
  !        [NH3]/[SO2] < 0.02         -> iratns = 3 ! very low
  !        [NH3]/[SO2] > 2            -> iratns = 2 ! high
  !        0.02 <= [NH3]/[SO2] <= 2   -> iratns = 1 ! low
  !     endif
  !
  elemental function Get_N_S_Regime( conc_NH3, conc_SO2 ) result ( iratns )

    ! --- in/out ---------------------------------

    real, intent(in)      ::  conc_NH3   ! NH3 concentration (ppb)
    real, intent(in)      ::  conc_SO2   ! SO2 concentration (ppb)
    integer               ::  iratns     ! regime index

    ! --- begin ----------------------------------

    if ( conc_SO2 > 0.5 ) then
      if ( conc_nh3 / conc_so2 < 0.02 ) then
        iratns = iratns_very_low
      else if ( conc_nh3 / conc_so2 > 2.0 ) then
        iratns = iratns_high
      else
        iratns = iratns_low
      end if
    else
      ! default:
      iratns = iratns_low
    end if

  end function Get_N_S_Regime


end module LE_DryDepos_Gas_DEPAC
