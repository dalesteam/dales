!> \file modslurbdata.f90
!  This file is part of DALES.
!
! DALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 2025-2026 Delft University of Technology
! Copyright 2022-2024 University of Helsinki
!
! This file was modified from the original version of the PALM SLUrb model by Sasu Karttunen, which is available at: https://gitlab.palm-model.org/releases/palm_model_system
! authors:
!   Sasu Karttunen <sasu.karttunen@helsinki.fi>
!   André van Ginkel <a.vanginkel@tudelft.nl>

module modslurbdata
    use modprecision, only: field_r
    implicit none
    save
    public

    logical :: lslurb            ! On/off switch LSM

    !-- Default surface description.
    REAL(field_r), DIMENSION(0:45,1:6) ::  building_pars_slurb  !< building default parameters derived from USM
    REAL(field_r), DIMENSION(0:14,1:5) ::  pavement_pars_slurb  !< pavement default parameters derived from LSM
    
    !-- Derived type for the SLUrb model.
    type surf_slurb

        real(field_r), allocatable ::  dt_max(:,:)  !< time step limit for model physical processes (s)

        real(field_r), allocatable ::  dz_road(:,:,:)  !< road layer thickness (m)
        real(field_r), allocatable ::  dz_roof(:,:,:)  !< roof layer thickness (m)
        real(field_r), allocatable ::  dz_wall(:,:,:)  !< wall layer thickness (m)
        real(field_r), allocatable ::  dz_win(:,:,:)   !< window layer thickness (m)
        real(field_r), allocatable ::  zw_win(:,:,:)   !< cumulative window thickness (m)
        !
        !--    Tile-aggregated quantities.
        real(field_r), allocatable ::  albedo_urb(:,:)       !< effective urban albedo (.)
        real(field_r), allocatable ::  emiss_urb(:,:)        !< effective urban emissivity (.)
        real(field_r), allocatable ::  ol_urb(:,:)           !< urban Obukhov length (m)
        real(field_r), allocatable ::  qsws_urb(:,:)         !< total urban latent heat flux (W m^-2)
        real(field_r), allocatable ::  rad_lw_in_urb(:,:)   !< incoming longwave radiation (W m^-2)
        real(field_r), allocatable ::  rad_lw_out_urb(:,:)  !< outgoing longwave radiation (W m^-2)
        real(field_r), allocatable ::  rad_sw_in_urb(:,:)   !< incoming shortwave radiation (W m^-2)
        real(field_r), allocatable ::  rad_sw_out_urb(:,:)  !< outgoing shortwave radiation (W m^-2)
        real(field_r), allocatable ::  ram_urb(:,:)         !< urban aerodynamic resistance for momentum (s m^-1)
        real(field_r), allocatable ::  rib_urb(:,:)         !< urban bulk-Richardson number (.)
        real(field_r), allocatable ::  shf_urb(:,:)         !< total urban sensible heat flux (W m^-2)
        real(field_r), allocatable ::  t_2m_urb(:,:)        !< urban 2-metre temperature (extrapolated) (K)
        real(field_r), allocatable ::  t_c_urb(:,:)         !< complete (area-weighted) urban surface temperature (K)
        real(field_r), allocatable ::  t_h_urb(:,:)         !< effective urban surface temperature (K)
        real(field_r), allocatable ::  thl_rad_urb(:,:)       !< urban radiative surface temperature (K)
        real(field_r), allocatable ::  usws_urb(:,:)        !< urban momentum flux (u-component) (kg m^-1 s^-2)
        real(field_r), allocatable ::  vsws_urb(:,:)        !< urban momentum flux (v-component) (kg m^-1 s^-2)
        real(field_r), allocatable ::  thlskin(:,:)         !< weighted urban roof + canopy liquid water potential temperature (K)
        real(field_r), allocatable ::  qtskin(:,:)          !< weighted urban roof + canopy specific humidity TODOSELF (kg kg^-1)
        !
        !--    Model prognostic variables.
        real(field_r), allocatable ::  m_liq_road_0(:,:)  !< liquid water reservoir on roads (m^3 m^-2)
        real(field_r), allocatable ::  m_liq_road_m(:,:)  !< prev. liquid water reservoir on roads (m^3 m^-2)
        real(field_r), allocatable ::  m_liq_roof_0(:,:)  !< liquid water reservoir on roofs (m^3 m^-2)
        real(field_r), allocatable ::  m_liq_roof_m(:,:)  !< prev. liquid water reservoir on roofs (m^3 m^-2)
        real(field_r), allocatable ::  q_can_0(:,:)       !< canyon mixing ratio (kg kg^-1)
        real(field_r), allocatable ::  q_can_m(:,:)       !< previous canyon mixing ratio (kg kg^-1)
        real(field_r), allocatable ::  t_can_0(:,:)       !< canyon air temperature (K)
        real(field_r), allocatable ::  t_can_m(:,:)       !< prev. canyon temperature (K)

        real(field_r), allocatable ::  t_road_0(:,:,:)      !< road temperature (K)
        real(field_r), allocatable ::  t_road_m(:,:,:)      !< prev. road temperature (K)
        real(field_r), allocatable ::  t_roof_0(:,:,:)      !< roof temperature (K)
        real(field_r), allocatable ::  t_roof_m(:,:,:)      !< prev. roof temperature (K)
        real(field_r), allocatable ::  t_wall_a_0(:,:,:)    !< wall A temperature (K)
        real(field_r), allocatable ::  t_wall_a_m(:,:,:)    !< prev. wall A temperature (K)
        real(field_r), allocatable ::  t_wall_b_0(:,:,:)    !< wall B temperature (K)
        real(field_r), allocatable ::  t_wall_b_m(:,:,:)    !< prev. wall B temperature (K)
        real(field_r), allocatable ::  t_win_a_0(:,:,:)     !< window A temperature (K)
        real(field_r), allocatable ::  t_win_a_m(:,:,:)     !< prev. window A temperature (K)
        real(field_r), allocatable ::  t_win_b_0(:,:,:)     !< window B temperature (K)
        real(field_r), allocatable ::  t_win_b_m(:,:,:)     !< prev. window B temperature (K)
        !
        !--    Tendencies of the prognostic variables.
        real(field_r), allocatable ::  tm_liq_road(:,:)     !< road liquid water reservoir tendency (m^3 m^-2 s^-1)
        real(field_r), allocatable ::  tm_liq_roof(:,:)     !< roof liquid water reservoir tendency (m^3 m^-2 s^-1)
        real(field_r), allocatable ::  tm_roof_runoff(:,:)  !< roof liquid water tendency due to runoff (m^3 m^-2 s^-1)
        real(field_r), allocatable ::  tm_road_runoff(:,:)  !< road liquid water tendency due to runoff (m^3 m^-2 s^-1)
        real(field_r), allocatable ::  tm_roof_precep(:,:)  !< roof liquid water tendency due to precipitation (m^3 m^-2 s^-1)
        real(field_r), allocatable ::  tm_road_precep(:,:)  !< road liquid water tendency due to precipitation (m^3 m^-2 s^-1)
        real(field_r), allocatable ::  tq_can(:,:)          !< canyon mixing ratio tendency (kg kg^-1 s^-1)
        real(field_r), allocatable ::  tt_can(:,:)          !< canyon temperature tendency (K s^-1)

        real(field_r), allocatable ::  tt_road(:,:,:)    !< road temperature tendency (K s^-1)
        real(field_r), allocatable ::  tt_roof(:,:,:)    !< roof temperature tendency (K s^-1)
        real(field_r), allocatable ::  tt_wall_a(:,:,:)  !< wall A temperature tendency (K s^-1)
        real(field_r), allocatable ::  tt_wall_b(:,:,:)  !< wall B temperature tendency (K s^-1)
        real(field_r), allocatable ::  tt_win_a(:,:,:)   !< window A temperature tendency (K s^-1)
        real(field_r), allocatable ::  tt_win_b(:,:,:)   !< window B temperature tendency (K s^-1)
        !
        !--    Diagnostic surface thermodynamic variables.
        real(field_r), allocatable ::  pt_road(:,:)    !< road surface potential temperature (K)
        real(field_r), allocatable ::  pt_roof(:,:)    !< roof surface potential temperature (K)
        real(field_r), allocatable ::  pt_wall_a(:,:)  !< wall A surface potential temperature (K)
        real(field_r), allocatable ::  pt_wall_b(:,:)  !< wall B surface potential temperature (K)
        real(field_r), allocatable ::  pt_win_a(:,:)   !< window A surface potential temperature (K)
        real(field_r), allocatable ::  pt_win_b(:,:)   !< window A surface potential temperature (K)
        real(field_r), allocatable ::  q_road(:,:)     !< road surface mixing ratio (kg kg^-1)
        real(field_r), allocatable ::  q_roof(:,:)     !< roof surface mixing ratio (kg kg^-1)
        real(field_r), allocatable ::  qs_road(:,:)    !< road surface saturation mixing ratio (kg kg^-1)
        real(field_r), allocatable ::  qs_roof(:,:)    !< roof surface saturation mixing ratio (kg kg^-1)
        real(field_r), allocatable ::  vpt_road(:,:)   !< road surface virtual potential temperature (K)
        real(field_r), allocatable ::  vpt_roof(:,:)   !< roof surface virtual potential temperature (K)
        !
        !--    Diagnostic internal sensible heat fluxes.
        real(field_r), allocatable ::  shf_can(:,:)       !< sensible heat flux between the street canyon and the atmosphere (W m^-2)
        real(field_r), allocatable ::  shf_external(:,:)  !< sensible heat flux external to the model (e.g. industry) (W m^-2)
        real(field_r), allocatable ::  shf_road(:,:)      !< road surface sensible heat flux (W m^-2)
        real(field_r), allocatable ::  shf_roof(:,:)      !< roof surface sensible heat flux (W m^-2)
        real(field_r), allocatable ::  shf_traffic(:,:)   !< traffic sensible heat flux (input-only) (W m^-2)
        real(field_r), allocatable ::  shf_wall_a(:,:)    !< wall A sensible heat flux (W m^-2)
        real(field_r), allocatable ::  shf_wall_b(:,:)    !< wall B sensible heat flux (W m^-2)
        real(field_r), allocatable ::  shf_win_a(:,:)     !< window A sensible heat flux (W m^-2)
        real(field_r), allocatable ::  shf_win_b(:,:)     !< window B sensible heat flux (W m^-2)
        !
        !--    Diagnostic internal latent heat fluxes.
        real(field_r), allocatable ::  qsws_can(:,:)       !< latent heat flux between the street canyon and the atmosphere (W m^-2)
        real(field_r), allocatable ::  qsws_external(:,:)  !< latent heat flux external to the model (e.g. industry) (W m^-2)
        real(field_r), allocatable ::  qsws_liq_road(:,:)  !< roof latent heat flux (liquid incl. precipitation) (W m^-2)
        real(field_r), allocatable ::  qsws_liq_roof(:,:)  !< roof latent heat flux (liquid incl. precipitation) (W m^-2)
        real(field_r), allocatable ::  qsws_road(:,:)      !< road latent heat flux (W m^-2)
        real(field_r), allocatable ::  qsws_roof(:,:)      !< roof latent heat flux (W m^-2)
        !
        !--    Liquid water coverages (storages).
        real(field_r), allocatable ::  c_liq_road(:,:)  !< liquid water coverage on road (.)
        real(field_r), allocatable ::  c_liq_roof(:,:)  !< liquid water coverage on roof (.)
        !
        !--    Diagnostic ground heat fluxes.
        real(field_r), allocatable ::  ghf_road(:,:)    !< road ground heat flux (W m^-2)
        real(field_r), allocatable ::  ghf_roof(:,:)    !< roof indoor heat flux (W m^-2)
        real(field_r), allocatable ::  ghf_wall_a(:,:)  !< wall A indoor heat flux (W m^-2)
        real(field_r), allocatable ::  ghf_wall_b(:,:)  !< wall B indoor heat flux (W m^-2)
        real(field_r), allocatable ::  ghf_win_a(:,:)   !< window A indoor heat flux (W m^-2)
        real(field_r), allocatable ::  ghf_win_b(:,:)   !< window B indoor heat flux (W m^-2)
        !
        !--    Model internal radiation fluxes.
        real(field_r), allocatable ::  rad_lw_net_can(:,:)     !< net longwave radiative at canyon top (downwards) (W m^-2)
        real(field_r), allocatable ::  rad_lw_net_road(:,:)    !< net longtwave radiative flux on road (W m^-2)
        real(field_r), allocatable ::  rad_lw_net_roof(:,:)    !< net longwave radiative flux on roof (W m^-2)
        real(field_r), allocatable ::  rad_lw_net_urb(:,:)     !< urban aggegated net longwave radiative flux (W m^-2)
        real(field_r), allocatable ::  rad_lw_net_wall_a(:,:)  !< net longwave radiative flux on wall A (W m^-2)
        real(field_r), allocatable ::  rad_lw_net_wall_b(:,:)  !< net longwave radiative flux wall B (W m^-2)
        real(field_r), allocatable ::  rad_lw_net_win_a(:,:)   !< net longwave radiative flux on wall A (W m^-2)
        real(field_r), allocatable ::  rad_lw_net_win_b(:,:)   !< net longwave radiative flux window B (W m^-2)
        real(field_r), allocatable ::  rad_sw_in_road(:,:)     !< incoming shortwave radiative flux on road (W m^-2)
        real(field_r), allocatable ::  rad_sw_in_win_a(:,:)    !< incoming shortwave radiative flux on window A (W m^-2)
        real(field_r), allocatable ::  rad_sw_in_win_b(:,:)    !< incoming shortwave radiative flux on window B (W m^-2)
        real(field_r), allocatable ::  rad_sw_net_road(:,:)    !< net shortwave radiative flux on road (W m^-2)
        real(field_r), allocatable ::  rad_sw_net_roof(:,:)    !< net shortwave radiative flux on roof (W m^-2)
        real(field_r), allocatable ::  rad_sw_net_urb(:,:)     !< urban aggegated net shortwave radiative flux (W m^-2)
        real(field_r), allocatable ::  rad_sw_net_wall_a(:,:)  !< net shortwave radiative flux on wall A (W m^-2)
        real(field_r), allocatable ::  rad_sw_net_wall_b(:,:)  !< net shortwave radiative flux on wall B (W m^-2)
        real(field_r), allocatable ::  rad_sw_net_win_a(:,:)   !< net shortwave radiative flux on window A (W m^-2)
        real(field_r), allocatable ::  rad_sw_net_win_b(:,:)   !< net shortwave radiative flux on wall B (W m^-2)
        !
        !--    Surface layer model diagnostic variables.
        real(field_r), allocatable ::  ol_can(:,:)      !< canyon top Obukhov length (m)
        real(field_r), allocatable ::  ol_road(:,:)     !< road Obukhov length (m)
        real(field_r), allocatable ::  ol_roof(:,:)     !< roof Obukhov length (m)
        real(field_r), allocatable ::  pt_can(:,:)      !< street canyon virtual potential temperature (K)
        real(field_r), allocatable ::  rib_can(:,:)     !< canyon top bulk Richardson number (.)
        real(field_r), allocatable ::  rib_road(:,:)    !< road bulk Richardson number (.)
        real(field_r), allocatable ::  rib_roof(:,:)    !< roof bulk Richardson number (.)
        real(field_r), allocatable ::  us_can(:,:)      !< friction velocity for canyon resistance calculation (m s^-1)
        real(field_r), allocatable ::  uv_abs_can(:,:)  !< horizontal wind speed in street caynon at half-height (m s^-1)
        real(field_r), allocatable ::  uv_eff_can(:,:)  !< effective horizontal wind speed in street canyon at half-height (m s^-1)
        real(field_r), allocatable ::  vpt_can(:,:)     !< street canyon virtual potential temperature (K)
        !
        !--    Aerodynamic resistances for heat.
        real(field_r), allocatable ::  rah_can(:,:)     !< street canyon air aerodynamic resistance for heat (s m^-1)
        real(field_r), allocatable ::  rah_facade(:,:)  !< wall and window aerodynamic resistance for heat (combined) (s m^-1)
        real(field_r), allocatable ::  rah_road(:,:)    !< road aerodynamic resistance for heat (s m^-1)
        real(field_r), allocatable ::  rah_roof(:,:)    !< roof aerodynamic resistance for heat (s m^-1)
        real(field_r), allocatable ::  rah_wall_a(:,:)  !< wall A aerodynamic resistance for heat (s m^-1)
        real(field_r), allocatable ::  rah_wall_b(:,:)  !< wall B aerodynamic resistance for heat (s m^-1)
        real(field_r), allocatable ::  rah_win_a(:,:)   !< wall A aerodynamic resistance for heat (s m^-1)
        real(field_r), allocatable ::  rah_win_b(:,:)   !< wall B aerodynamic resistance for heat (s m^-1)
        !
        !--    Local friction velocities for roofs and roads.
        real(field_r), allocatable ::  us_road(:,:)  !< friction velocity for roads (m s^-1)
        real(field_r), allocatable ::  us_roof(:,:)  !< friction velocity for roofs (m s^-1)
        !
        !--    Diagnostic variables, defined at the first atmospheric grid level.
        real(field_r), allocatable ::  pt1(:,:)      !< potential temperature (K)
        real(field_r), allocatable ::  q1(:,:)       !< specific humidity (kg kg^-1)
        real(field_r), allocatable ::  us_urb(:,:)   !< friction velocity (m s^-1)
        real(field_r), allocatable ::  uv_abs1(:,:)  !< horizontal wind speed (m s^-1)
        real(field_r), allocatable ::  uv_eff1(:,:)  !< effective horizontal wind speed (m s^-1)
        real(field_r), allocatable ::  vpt1(:,:)     !< virtual potential temperature (K)
        !
        !--    Parameters for the whole urban tile.
        LOGICAL,  allocatable ::  anisotropic_canyon(:,:)  !< boolean flag to mark anisotropic canyon
        real(field_r), allocatable ::  f_bld(:,:)               !< fractional area occupied by buldings (plan area fraction) (.)
        real(field_r), allocatable ::  f_bld_frn(:,:)           !< frontal area fraction of buildings (.)
        real(field_r), allocatable ::  f_win(:,:)               !< window fraction (.)
        real(field_r), allocatable ::  h_bld(:,:)               !< building height (m)
        real(field_r), allocatable ::  hw_can(:,:)              !< canyon aspect ratio (.)
        real(field_r), allocatable ::  svf_road(:,:)            !< sky-view factor for road (.)
        real(field_r), allocatable ::  svf_wall(:,:)            !< sky-view-factor for walls (.)
        real(field_r), allocatable ::  theta_can(:,:)           !< canyon orientation / road direction in radians (.)
        real(field_r), allocatable ::  z0_urb(:,:)              !< aerodynamic roughness length of the urban surface (m)
        !
        !--    Material properties.
        real(field_r), allocatable ::  albedo_road(:,:)         !< albedo of the road (.)
        real(field_r), allocatable ::  albedo_roof(:,:)         !< albedo of the roof (.)
        real(field_r), allocatable ::  albedo_wall(:,:)         !< albedo of the wall (.)
        real(field_r), allocatable ::  albedo_wall_win(:,:)     !< weighted average of wall and window albedos for reflections (.)
        real(field_r), allocatable ::  albedo_win(:,:)          !< albedo of the window (.)
        real(field_r), allocatable ::  emiss_road(:,:)          !< emissivity of the road (.)
        real(field_r), allocatable ::  emiss_roof(:,:)          !< emissivity of the roof (.)
        real(field_r), allocatable ::  emiss_wall(:,:)          !< emissivity of the wall (.)
        real(field_r), allocatable ::  emiss_win(:,:)           !< emissivity of the window (.)
        real(field_r), allocatable ::  transmissivity_win(:,:)  !< transmissivity of the window layers (.)
        real(field_r), allocatable ::  z0_road(:,:)             !< aerodynamic roughness length for momentum for roads (m)
        real(field_r), allocatable ::  z0_roof(:,:)             !< aerodynamic roughness length for momentum of roofs (m)
        real(field_r), allocatable ::  z0_wall(:,:)             !< aerodynamic roughness length for walls and windows (m)
        real(field_r), allocatable ::  z0h_road(:,:)            !< aerodynamic roughness length for heat for roads (m)
        real(field_r), allocatable ::  z0h_roof(:,:)            !< aerodynamic roughness length for heat for roofs (m)

        real(field_r), allocatable ::  absorption_win(:,:,:)  !< fraction of absorbed shortwave radiation over glass sheet (.)
        real(field_r), allocatable ::  c_road(:,:,:)          !< total (specific c * layer depth) heat capacity of the road (J m^-2 K^-1)
        real(field_r), allocatable ::  c_roof(:,:,:)          !< total (specific c * layer depth) heat capacity of the roof (J m^-2 K^-1)
        real(field_r), allocatable ::  c_wall(:,:,:)          !< total (specific c * layer depth) heat capacity of the wall (J m^-2 K^-1)
        real(field_r), allocatable ::  c_win(:,:,:)           !< total (specific c * layer depth) heat heat capacity of the window (J m^-2 K^-1)
        real(field_r), allocatable ::  lambda_road(:,:,:)     !< thermal conductivity of the road (W m^-1 K^-1)
        real(field_r), allocatable ::  lambda_roof(:,:,:)     !< thermal conductivity of the roof (W m^-1 K^-1)
        real(field_r), allocatable ::  lambda_wall(:,:,:)     !< thermal conductivity of the wall (W m^-1 K^-1)
        real(field_r), allocatable ::  lambda_win(:,:,:)      !< effective thermal conductivity of the window (W m^-1 K^-1)
        !
        !--    Building indoor parameters.
        real(field_r), allocatable ::  t_indoor(:,:)  !< building indoor temperature (K)
        !
        !--    Soil parameters.
        real(field_r), allocatable ::  t_soil(:,:)  !< fixed soil top temperature (K)
        !
        !--    Pre-computed total layer conductivities.
        real(field_r), allocatable ::  conductivity_road(:,:,:)  !< total conductivity bewtween road layers (lambda_h / dz) (W m^-2 K^-1)
        real(field_r), allocatable ::  conductivity_roof(:,:,:)  !< total conductivity between roof layers (lambda_h / dz) (W m^-2 K^-1)
        real(field_r), allocatable ::  conductivity_wall(:,:,:)  !< total conductivity between wall layers (lambda_h / dz) (W m^-2 K^-1)
        real(field_r), allocatable ::  conductivity_win(:,:,:)   !< total conductivity between window layers (lambda_h / dz) (W m^-2 K^-1)
        !
        !--    Pre-computed variables and coefficients.
        real(field_r), allocatable ::  sw_ref_denom(:,:)      !< SW radiation reflection denominator (.)
        real(field_r), allocatable ::  uv_abs_can_coef(:,:)   !< coefficient for the canyon wind speed (.)
        real(field_r), allocatable ::  wall_hor_a_ratio(:,:)  !< wall-to-horizontal area ratio (unused) (.)
        real(field_r), allocatable ::  z_mo(:,:)              !< reference height for MOST for the atmosphere (m)
        real(field_r), allocatable ::  z_mo_can(:,:)          !< canyon reference height for MOST (canyon half-height) (m)

        real(field_r), allocatable ::  lw_road_coef(:,:,:)  !< LW radiation coefficients for roads (W m^-2 K^-4 & . & W m^-2 K^-4 ...)
        real(field_r), allocatable ::  lw_roof_coef(:,:,:)  !< LW radiation coefficients for roofs (W m^-2 K^-4 & . & W m^-2 K^-4 ...)
        real(field_r), allocatable ::  lw_wall_coef(:,:,:)  !< LW radiation coefficients for walls (W m^-2 K^-4 & . & W m^-2 K^-4 ...)
        real(field_r), allocatable ::  lw_win_coef(:,:,:)   !< LW radiation coefficients for walls (W m^-2 K^-4 & . & W m^-2 K^-4 ...)

    end type surf_slurb

    type(surf_slurb) :: slurb_tile

    real(field_r), allocatable ::  ln_z_z0_roof(:,:)   !< temporary array to store logarithm ZELFTODO (.)
    real(field_r), allocatable ::  ln_z_z0h_roof(:,:)  !< temporary array to store logarithm (.)
    real(field_r), allocatable ::  ln_z_z0_urb(:,:)    !< temporary array to store logarithm (.)
    real(field_r), allocatable ::  pt_surface(:,:)     !< temporary array to store weighted temperature (K)
    
    real(field_r), allocatable ::  ln_z_z0_road(:,:)   !< temporary array to store logarithm ZELFTODO (.)
    real(field_r), allocatable ::  ln_z_z0h_road(:,:)  !< temporary array to store logarithm (.)

    REAL(field_r) ::  dt_slurb = HUGE( 1.0_field_r )  !< maximum allowed timestep of SLUrb

    !
    !-- Model constants.
    REAL(field_r) ::  drho_l_lv  !< 1/(rho_l * l_v) (J^-1 m^3)
    REAL(field_r) ::  rho_lv     !< rho_surface * l_v (J m^-3)

    !
    !-- Parameter defaults.
    REAL(field_r), PARAMETER ::  m_liq_max_road = 1.0E-3_field_r  !< maximum capacity of the liquid water reservoir on roads (i,j) (m^3 m^-2)
    REAL(field_r), PARAMETER ::  m_liq_max_roof = 1.0E-3_field_r  !< maximum capacity of the liquid water reservoir on roofs (i,j) (m^3 m^-2)
    ! REAL(field_r), PARAMETER ::  m_liq_max_road = 100_field_r  !< maximum capacity of the liquid water reservoir on roads (i,j) (m^3 m^-2)
    ! REAL(field_r), PARAMETER ::  m_liq_max_roof = 100_field_r  !< maximum capacity of the liquid water reservoir on roofs (i,j) (m^3 m^-2)
    REAL(field_r), PARAMETER ::  rah_max   = 1.0E6_field_r        !< maximum aerodynamic resistance for scalars (s m^-1)
    REAL(field_r), PARAMETER ::  rah_min   = 1.0_field_r          !< minimum aerodynamic resistance for scalars (s m^-1)
    REAL(field_r), PARAMETER ::  ram_min   = 1.0_field_r          !< minimum aerodynamic resistance for momentum (s m^-1) (TODOSELF)
    REAL(field_r), PARAMETER ::  urb_thres = 1.0E-2_field_r       !< minimum urban fraction to consider (1%) (.)
    REAL(field_r), PARAMETER ::  us_min    = 1.0E-8_field_r       !< minimum friction velocity (m s^-1)
    REAL(field_r), PARAMETER ::  zeta_min  = 1.0E-3_field_r       !< minimum stability parameter absolute value (neutral limit) (.)
    !
    !-- slurb_parameters namelist defaults.
    CHARACTER(LEN=20) ::  aero_roughness_heat = 'kanda'                !< SLURrb namelist parameter
    CHARACTER(LEN=20) ::  facade_resistance_parametrization = 'doe-2'  !< SLURrb namelist parameter
    CHARACTER(LEN=20) ::  street_canyon_wspeed_factor = 'surfex'       !< SLURrb namelist parameter

    integer ::  building_type = 2     !< SLURrb namelist parameter
    integer ::  n_layers_roads = 4    !< SLURrb namelist parameter
    integer ::  n_layers_roofs = 4    !< SLURrb namelist parameter
    integer ::  n_layers_walls = 4    !< SLURrb namelist parameter
    integer ::  n_layers_windows = 4  !< SLURrb namelist parameter
    integer ::  pavement_type = 2     !< SLURrb namelist parameter

    LOGICAL ::  anisotropic_street_canyons = .FALSE.  !< SLURrb namelist parameter
    LOGICAL ::  moist_physics = .true.                !< SLURrb namelist parameter
    logical ::  lread_from_netcdf = .true.             !< SLURrb namelist parameter

    REAL(field_r) ::  building_frontal_area_fraction = -9999.0_field_r  !< SLURrb namelist parameter (.)
    REAL(field_r) ::  building_height = -9999.0_field_r                 !< SLURrb namelist parameter (m)
    REAL(field_r) ::  building_indoor_temperature =  -9999.0_field_r    !< SLURrb namelist parameter (K)
    REAL(field_r) ::  building_plan_area_fraction = -9999.0_field_r     !< SLURrb namelist parameter (.)
    REAl(field_r) ::  deep_soil_temperature = -9999.0_field_r           !< SLURrb namelist parameter (K)
    REAL(field_r) ::  qsws_external = 0.0_field_r                       !< SLURrb namelist parameter (W m^-2 s^-1)
    REAL(field_r) ::  shf_external = 0.0_field_r                        !< SLURrb namelist parameter (W m^-2 s^-1)
    REAL(field_r) ::  shf_traffic = 0.0_field_r                         !< SLURrb namelist parameter (W m^-2 s^-1)
    REAL(field_r) ::  street_canyon_aspect_ratio = -9999.0_field_r      !< SLURrb namelist parameter (.)
    REAL(field_r) ::  street_canyon_orientation = -9999.0_field_r       !< SLURrb namelist parameter (.)
    REAL(field_r) ::  urban_fraction = -9999.0_field_r                  !< SLURrb namelist parameter (.)
    REAL(field_r) ::  urban_roughness_length = -9999.0_field_r          !< SLURrb namelist parameter (m)
    REAL(field_r) ::  window_fraction = -9999.0_field_r                 !< SLURrb namelist parameter (.)


    REAL(field_r), PARAMETER ::  ol_max   = 1.0E6_field_r   !< allowed absolute maximum value Obukhov length (m)
    REAL(field_r), PARAMETER ::  ol_min   = 1.0E-6_field_r  !< allowed absolute minimum value Obukhov length (m)
    REAL(field_r), PARAMETER ::  ol_tol   = 1.0E-4_field_r  !< convergence limit for Obukhov length, relative tolerance (m)
    REAL(field_r), PARAMETER ::  rib_max  = 1.0E1_field_r   !< maximum bulk Richardson number (absolute value) (.)

    !-- Internal logical switches for character-based namelist settings.
    !TODO ADD CHECKS
    LOGICAL ::  facade_rah_doe       = .TRUE.  !< facade resistance parameterization using DOE-2
    LOGICAL ::  facade_rah_kray      = .FALSE.  !< facade resistance parameterization using Krayenhoff&Voogt (2007)
    LOGICAL ::  facade_rah_rowley    = .FALSE.  !< facade resistance parameterization using Rowley (1932)
    LOGICAL ::  roughness_kanda      = .FALSE.  !< roughness parameterization of horizontal surfaces using Kanda et al. (2007)
    LOGICAL ::  uv_can_factor_kray   = .FALSE.  !< street canyon wind speed factor following Krayenhoff&Voogt (2007)
    LOGICAL ::  uv_can_factor_masson = .FALSE.  !< street canyon wind speed factor following Masson (2000)
    LOGICAL ::  uv_can_factor_surfex = .TRUE.  !< street canyon wind speed factor following the SURFEX model

    !-- Default subsurface layer configuration.
    INTEGER ::  nzt_wall  !< top of the wall model (outer surface)
    INTEGER ::  nzb_wall  !< bottom of the wall model (inside surface)
    INTEGER ::  nzt_win   !< top of the window model (outer surface)
    INTEGER ::  nzb_win   !< bottom of the window model (inside surface)
    INTEGER ::  nzt_roof  !< top of the roof model (outer surface)
    INTEGER ::  nzb_roof  !< bottom the roof model (inside surface)
    INTEGER ::  nzt_road  !< top of the road model
    INTEGER ::  nzb_road  !< bottom of the road model

    real(field_r), allocatable :: fraction_slurb(:,:) !< (.)


    real(field_r) :: output_fill_value = -99999.0_field_r
    logical :: data_output_raw = .false.
    logical :: spinup = .false.
    logical :: calc_t_2m = .true.
    logical :: calc_t_c = .true.
    logical :: calc_t_h = .true.

    logical :: enable_slurb = .false. !< switch to enable slurb model. Is set to true if any slurb tile found.

contains
!--------------------------------------------------------------------------------------------------!
!   MODLULE PREDEFINED PARAMETERS
!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Default parameters for the building types. These are based on the PALM urban surface mod. (urban_surface_mod.f90)
!> These values can only be considered valid for german buildings, and are not necessarily representative for other regions.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE slurb_default_pars

!
!-- Residential, < 1950.
    building_pars_slurb(:,1) = (/                                                                  &
       0.18_field_r,        &   !< parameter 0   - [-] window fraction
       0.02_field_r,        &   !< parameter 1   - [m] 1st roof layer thickness (outside)
       0.04_field_r,        &   !< parameter 2   - [m] 2nd roof layer thickness
       0.02_field_r,        &   !< parameter 3   - [m] 3rd roof layer thickness
       0.02_field_r,        &   !< parameter 4   - [m] 4th roof layer thickness (inside)
       1.51200E6_field_r,   &   !< parameter 5   - [J/(m3*K)] specific heat capacity 1st roof layer (outside)
       0.70965E6_field_r,   &   !< parameter 6   - [J/(m3*K)] specific heat capacity 2nd roof layer
       0.70965E6_field_r,   &   !< parameter 7   - [J/(m3*K)] specific heat capacity 3rd roof layer
       1.52600E6_field_r,   &   !< parameter 8   - [J/(m3*K)] specific heat capacity 4th roof layer (inside)
       0.520_field_r,       &   !< parameter 9   - [W/(m*K)] thermal conductivity 1st roof layer (outside)
       0.120_field_r,       &   !< parameter 10  - [W/(m*K)] thermal conductivity 2nd roof layer
       0.120_field_r,       &   !< parameter 11  - [W/(m*K)] thermal conductivity 3rd roof layer
       0.700_field_r,       &   !< parameter 12  - [W/(m*K)] thermal conductivity 4th roof layer (inside)
       0.15_field_r,        &   !< parameter 13  - [m] z0 roughness length for momentum
       0.17_field_r,        &   !< parameter 14  - [-] albedo
       0.90_field_r,        &   !< parameter 15  - [-] emissivity
       0.02_field_r,        &   !< parameter 16  - [m] 1st wall layer thickness (outside)
       0.18_field_r,        &   !< parameter 17  - [m] 2nd wall layer thickness
       0.18_field_r,        &   !< parameter 18  - [m] 3rd wall layer thickness
       0.02_field_r,        &   !< parameter 19  - [m] 4th wall layer thickness
       1.5200E6_field_r,    &   !< parameter 20  - [J/(m3*K)] specific heat capacity 1st wall layer (outside)
       1.5120E6_field_r,    &   !< parameter 21  - [J/(m3*K)] specific heat capacity 2nd wall layer
       1.5120E6_field_r,    &   !< parameter 22  - [J/(m3*K)] specific heat capacity 3rd wall layer
       1.5260E6_field_r,    &   !< parameter 23  - [J/(m3*K)] specific heat capacity 4th wall layer (inside)
       0.930_field_r,       &   !< parameter 24  - [W/(m*K)] thermal conductivity 1st wall layer (outside)
       0.810_field_r,       &   !< parameter 25  - [W/(m*K)] thermal conductivity 2nd wall layer
       0.810_field_r,       &   !< parameter 26  - [W/(m*K)] thermal conductivity 3rd wall layer
       0.700_field_r,       &   !< parameter 27  - [W/(m*K)] thermal conductivity 4th wall layer (inside)
       0.001_field_r,       &   !< parameter 28  - [m] z0 roughness length for momentum
       0.30_field_r,        &   !< parameter 29  - [-] albedo
       0.93_field_r,        &   !< parameter 30  - [-] emissivity
       0.02_field_r,        &   !< parameter 31  - [m] 1st window layer thickness (glass sheet + air total) (outside)
       0.02_field_r,        &   !< parameter 32  - [m] 2rd window layer thickness
       0.02_field_r,        &   !< parameter 33  - [m] 3rd window layer thickness
       0.02_field_r,        &   !< parameter 34  - [m] 4th window layer thickness (inside)
       1.736E6_field_r,     &   !< parameter 35  - [J/(m3*K)] specific heat capacity 1st window layer (outside)
       1.736E6_field_r,     &   !< parameter 36  - [J/(m3*K)] specific heat capacity 2nd window layer
       1.736E6_field_r,     &   !< parameter 37  - [J/(m3*K)] specific heat capacity 3rd window layer
       1.736E6_field_r,     &   !< parameter 38  - [J/(m3*K)] specific heat capacity 4th window layer (inside)
       0.45_field_r,        &   !< parameter 39  - [W/(m*K)] thermal conductivity 1st window layer (outside)
       0.45_field_r,        &   !< parameter 40  - [W/(m*K)] thermal conductivity 2nd window layer
       0.45_field_r,        &   !< parameter 41  - [W/(m*K)] thermal conductivity 3rd window layer
       0.45_field_r,        &   !< parameter 42  - [W/(m*K)] thermal conductivity 4th window layer (inside)
       0.70_field_r,        &   !< parameter 43  - [-] transmissivity
       0.12_field_r,        &   !< parameter 44  - [-] albedo
       0.91_field_r         &   !< parameter 45  - [-] emissivity
    /)

!
!-- Residential, 1950 - 2000.
    building_pars_slurb(:,2) = (/                                                                  &
       0.25_field_r,        &   !< parameter 0   - [-] window fraction
       0.02_field_r,        &   !< parameter 1   - [m] 1st roof layer thickness (outside)
       0.15_field_r,        &   !< parameter 2   - [m] 2nd roof layer thickness
       0.20_field_r,        &   !< parameter 3   - [m] 3rd roof layer thickness
       0.02_field_r,        &   !< parameter 4   - [m] 4th roof layer thickness (inside)
       1.70000E6_field_r,   &   !< parameter 5   - [J/(m3*K)] specific heat capacity 1st roof layer (outside)
       0.07920E6_field_r,   &   !< parameter 6   - [J/(m3*K)] specific heat capacity 2nd roof layer
       2.11200E6_field_r,   &   !< parameter 7   - [J/(m3*K)] specific heat capacity 3rd roof layer
       1.52600E6_field_r,   &   !< parameter 8   - [J/(m3*K)] specific heat capacity 4th roof layer (inside)
       0.160_field_r,       &   !< parameter 9   - [W/(m*K)] thermal conductivity 1st roof layer (outside)
       0.046_field_r,       &   !< parameter 10  - [W/(m*K)] thermal conductivity 2nd roof layer
       2.100_field_r,       &   !< parameter 11  - [W/(m*K)] thermal conductivity 3rd roof layer
       0.700_field_r,       &   !< parameter 12  - [W/(m*K)] thermal conductivity 4th roof layer (inside)
       0.15_field_r,        &   !< parameter 13  - [m] z0 roughness length for momentum
       0.10_field_r,        &   !< parameter 14  - [-] albedo
       0.95_field_r,        &   !< parameter 15  - [-] emissivity
       0.02_field_r,        &   !< parameter 16  - [m] 1st wall layer thickness (outside)
       0.06_field_r,        &   !< parameter 17  - [m] 2nd wall layer thickness
       0.24_field_r,        &   !< parameter 18  - [m] 3rd wall layer thickness
       0.02_field_r,        &   !< parameter 19  - [m] 4th wall layer thickness
       1.5200E6_field_r,    &   !< parameter 20  - [J/(m3*K)] specific heat capacity 1st wall layer (outside)
       0.0792E6_field_r,    &   !< parameter 21  - [J/(m3*K)] specific heat capacity 2nd wall layer
       2.1120E6_field_r,    &   !< parameter 22  - [J/(m3*K)] specific heat capacity 3rd wall layer
       1.5260E6_field_r,    &   !< parameter 23  - [J/(m3*K)] specific heat capacity 4th wall layer (inside)
       0.930_field_r,       &   !< parameter 24  - [W/(m*K)] thermal conductivity 1st wall layer (outside)
       0.046_field_r,       &   !< parameter 25  - [W/(m*K)] thermal conductivity 2nd wall layer
       2.100_field_r,       &   !< parameter 26  - [W/(m*K)] thermal conductivity 3rd wall layer
       0.700_field_r,       &   !< parameter 27  - [W/(m*K)] thermal conductivity 4th wall layer (inside)
       0.001_field_r,       &   !< parameter 28  - [m] z0 roughness length for momentum
       0.30_field_r,        &   !< parameter 29  - [-] albedo
       0.93_field_r,        &   !< parameter 30  - [-] emissivity
       0.02_field_r,        &   !< parameter 31  - [m] 1st window layer thickness (glass sheet + air total) (outside)
       0.02_field_r,        &   !< parameter 32  - [m] 2rd window layer thickness
       0.02_field_r,        &   !< parameter 33  - [m] 3rd window layer thickness
       0.02_field_r,        &   !< parameter 34  - [m] 4th window layer thickness (inside)
       1.736E6_field_r,     &   !< parameter 35  - [J/(m3*K)] specific heat capacity 1st window layer (outside)
       1.736E6_field_r,     &   !< parameter 36  - [J/(m3*K)] specific heat capacity 2nd window layer
       1.736E6_field_r,     &   !< parameter 37  - [J/(m3*K)] specific heat capacity 3rd window layer
       1.736E6_field_r,     &   !< parameter 38  - [J/(m3*K)] specific heat capacity 4th window layer (inside)
       0.18_field_r,        &   !< parameter 39  - [W/(m*K)] thermal conductivity 1st window layer (outside)
       0.18_field_r,        &   !< parameter 40  - [W/(m*K)] thermal conductivity 2nd window layer
       0.18_field_r,        &   !< parameter 41  - [W/(m*K)] thermal conductivity 3rd window layer
       0.18_field_r,        &   !< parameter 42  - [W/(m*K)] thermal conductivity 4th window layer (inside)
       0.65_field_r,        &   !< parameter 43  - [-] transmissivity
       0.15_field_r,        &   !< parameter 44  - [-] albedo
       0.87_field_r         &   !< parameter 45  - [-] emissivity
    /)

!
!-- Residential, > 2000.
    building_pars_slurb(:,3) = (/                                                                  &
       0.29_field_r,        &   !< parameter 0   - [-] window fraction
       0.02_field_r,        &   !< parameter 1   - [m] 1st roof layer thickness (outside)
       0.04_field_r,        &   !< parameter 2   - [m] 2nd roof layer thickness
       0.30_field_r,        &   !< parameter 3   - [m] 3rd roof layer thickness
       0.02_field_r,        &   !< parameter 4   - [m] 4th roof layer thickness (inside)
       3.75360E6_field_r,   &   !< parameter 5   - [J/(m3*K)] specific heat capacity 1st roof layer (outside)
       0.70965E6_field_r,   &   !< parameter 6   - [J/(m3*K)] specific heat capacity 2nd roof layer
       0.07920E6_field_r,   &   !< parameter 7   - [J/(m3*K)] specific heat capacity 3rd roof layer
       1.52600E6_field_r,   &   !< parameter 8   - [J/(m3*K)] specific heat capacity 4th roof layer (inside)
       0.520_field_r,       &   !< parameter 9   - [W/(m*K)] thermal conductivity 1st roof layer (outside)
       0.120_field_r,       &   !< parameter 10  - [W/(m*K)] thermal conductivity 2nd roof layer
       0.035_field_r,       &   !< parameter 11  - [W/(m*K)] thermal conductivity 3rd roof layer
       0.700_field_r,       &   !< parameter 12  - [W/(m*K)] thermal conductivity 4th roof layer (inside)
       0.15_field_r,        &   !< parameter 13  - [m] z0 roughness length for momentum
       0.17_field_r,        &   !< parameter 14  - [-] albedo
       0.92_field_r,        &   !< parameter 15  - [-] emissivity
       0.02_field_r,        &   !< parameter 16  - [m] 1st wall layer thickness (outside)
       0.20_field_r,        &   !< parameter 17  - [m] 2nd wall layer thickness
       0.36_field_r,        &   !< parameter 18  - [m] 3rd wall layer thickness
       0.02_field_r,        &   !< parameter 19  - [m] 4th wall layer thickness
       1.5200E6_field_r,    &   !< parameter 20  - [J/(m3*K)] specific heat capacity 1st wall layer (outside)
       0.0792E6_field_r,    &   !< parameter 21  - [J/(m3*K)] specific heat capacity 2nd wall layer
       1.3400E6_field_r,    &   !< parameter 22  - [J/(m3*K)] specific heat capacity 3rd wall layer
       1.5260E6_field_r,    &   !< parameter 23  - [J/(m3*K)] specific heat capacity 4th wall layer (inside)
       0.930_field_r,       &   !< parameter 24  - [W/(m*K)] thermal conductivity 1st wall layer (outside)
       0.035_field_r,       &   !< parameter 25  - [W/(m*K)] thermal conductivity 2nd wall layer
       0.680_field_r,       &   !< parameter 26  - [W/(m*K)] thermal conductivity 3rd wall layer
       0.700_field_r,       &   !< parameter 27  - [W/(m*K)] thermal conductivity 4th wall layer (inside)
       0.001_field_r,       &   !< parameter 28  - [m] z0 roughness length for momentum
       0.37_field_r,        &   !< parameter 29  - [-] albedo
       0.93_field_r,        &   !< parameter 30  - [-] emissivity
       0.02_field_r,        &   !< parameter 31  - [m] 1st window layer thickness (glass sheet + air total) (outside)
       0.02_field_r,        &   !< parameter 32  - [m] 2rd window layer thickness
       0.02_field_r,        &   !< parameter 33  - [m] 3rd window layer thickness
       0.02_field_r,        &   !< parameter 34  - [m] 4th window layer thickness (inside)
       1.736E6_field_r,     &   !< parameter 35  - [J/(m3*K)] specific heat capacity 1st window layer (outside)
       1.736E6_field_r,     &   !< parameter 36  - [J/(m3*K)] specific heat capacity 2nd window layer
       1.736E6_field_r,     &   !< parameter 37  - [J/(m3*K)] specific heat capacity 3rd window layer
       1.736E6_field_r,     &   !< parameter 38  - [J/(m3*K)] specific heat capacity 4th window layer (inside)
       0.11_field_r,        &   !< parameter 39  - [W/(m*K)] thermal conductivity 1st window layer (outside)
       0.11_field_r,        &   !< parameter 40  - [W/(m*K)] thermal conductivity 2nd window layer
       0.11_field_r,        &   !< parameter 41  - [W/(m*K)] thermal conductivity 3rd window layer
       0.11_field_r,        &   !< parameter 42  - [W/(m*K)] thermal conductivity 4th window layer (inside)
       0.57_field_r,        &   !< parameter 43  - [-] transmissivity
       0.18_field_r,        &   !< parameter 44  - [-] albedo
       0.80_field_r         &   !< parameter 45  - [-] emissivity
    /)

!
!-- Office, < 1950.
    building_pars_slurb(:,4) = (/                                                                  &
       0.18_field_r,        &   !< parameter 0   - [-] window fraction
       0.02_field_r,        &   !< parameter 1   - [m] 1st roof layer thickness (outside)
       0.04_field_r,        &   !< parameter 2   - [m] 2nd roof layer thickness
       0.02_field_r,        &   !< parameter 3   - [m] 3rd roof layer thickness
       0.02_field_r,        &   !< parameter 4   - [m] 4th roof layer thickness (inside)
       1.51200E6_field_r,   &   !< parameter 5   - [J/(m3*K)] specific heat capacity 1st roof layer (outside)
       0.70965E6_field_r,   &   !< parameter 6   - [J/(m3*K)] specific heat capacity 2nd roof layer
       0.70965E6_field_r,   &   !< parameter 7   - [J/(m3*K)] specific heat capacity 3rd roof layer
       1.52600E6_field_r,   &   !< parameter 8   - [J/(m3*K)] specific heat capacity 4th roof layer (inside)
       0.520_field_r,       &   !< parameter 9   - [W/(m*K)] thermal conductivity 1st roof layer (outside)
       0.120_field_r,       &   !< parameter 10  - [W/(m*K)] thermal conductivity 2nd roof layer
       0.120_field_r,       &   !< parameter 11  - [W/(m*K)] thermal conductivity 3rd roof layer
       0.700_field_r,       &   !< parameter 12  - [W/(m*K)] thermal conductivity 4th roof layer (inside)
       0.15_field_r,        &   !< parameter 13  - [m] z0 roughness length for momentum
       0.17_field_r,        &   !< parameter 14  - [-] albedo
       0.90_field_r,        &   !< parameter 15  - [-] emissivity
       0.02_field_r,        &   !< parameter 16  - [m] 1st wall layer thickness (outside)
       0.18_field_r,        &   !< parameter 17  - [m] 2nd wall layer thickness
       0.18_field_r,        &   !< parameter 18  - [m] 3rd wall layer thickness
       0.02_field_r,        &   !< parameter 19  - [m] 4th wall layer thickness
       1.5200E6_field_r,    &   !< parameter 20  - [J/(m3*K)] specific heat capacity 1st wall layer (outside)
       1.5120E6_field_r,    &   !< parameter 21  - [J/(m3*K)] specific heat capacity 2nd wall layer
       1.5120E6_field_r,    &   !< parameter 22  - [J/(m3*K)] specific heat capacity 3rd wall layer
       1.5260E6_field_r,    &   !< parameter 23  - [J/(m3*K)] specific heat capacity 4th wall layer (inside)
       0.930_field_r,       &   !< parameter 24  - [W/(m*K)] thermal conductivity 1st wall layer (outside)
       0.810_field_r,       &   !< parameter 25  - [W/(m*K)] thermal conductivity 2nd wall layer
       0.810_field_r,       &   !< parameter 26  - [W/(m*K)] thermal conductivity 3rd wall layer
       0.700_field_r,       &   !< parameter 27  - [W/(m*K)] thermal conductivity 4th wall layer (inside)
       0.001_field_r,       &   !< parameter 28  - [m] z0 roughness length for momentum
       0.30_field_r,        &   !< parameter 29  - [-] albedo
       0.93_field_r,        &   !< parameter 30  - [-] emissivity
       0.02_field_r,        &   !< parameter 31  - [m] 1st window layer thickness (glass sheet + air total) (outside)
       0.02_field_r,        &   !< parameter 32  - [m] 2rd window layer thickness
       0.02_field_r,        &   !< parameter 33  - [m] 3rd window layer thickness
       0.02_field_r,        &   !< parameter 34  - [m] 4th window layer thickness (inside)
       1.736E6_field_r,     &   !< parameter 35  - [J/(m3*K)] specific heat capacity 1st window layer (outside)
       1.736E6_field_r,     &   !< parameter 36  - [J/(m3*K)] specific heat capacity 2nd window layer
       1.736E6_field_r,     &   !< parameter 37  - [J/(m3*K)] specific heat capacity 3rd window layer
       1.736E6_field_r,     &   !< parameter 38  - [J/(m3*K)] specific heat capacity 4th window layer (inside)
       0.45_field_r,        &   !< parameter 39  - [W/(m*K)] thermal conductivity 1st window layer (outside)
       0.45_field_r,        &   !< parameter 40  - [W/(m*K)] thermal conductivity 2nd window layer
       0.45_field_r,        &   !< parameter 41  - [W/(m*K)] thermal conductivity 3rd window layer
       0.45_field_r,        &   !< parameter 42  - [W/(m*K)] thermal conductivity 4th window layer (inside)
       0.70_field_r,        &   !< parameter 43  - [-] transmissivity
       0.12_field_r,        &   !< parameter 44  - [-] albedo
       0.91_field_r         &   !< parameter 45  - [-] emissivity
    /)

!
!-- Office, 1950 - 2000.
    building_pars_slurb(:,5) = (/                                                                  &
       0.25_field_r,        &   !< parameter 0   - [-] window fraction
       0.02_field_r,        &   !< parameter 1   - [m] 1st roof layer thickness (outside)
       0.15_field_r,        &   !< parameter 2   - [m] 2nd roof layer thickness
       0.20_field_r,        &   !< parameter 3   - [m] 3rd roof layer thickness
       0.02_field_r,        &   !< parameter 4   - [m] 4th roof layer thickness (inside)
       1.70000E6_field_r,   &   !< parameter 5   - [J/(m3*K)] specific heat capacity 1st roof layer (outside)
       0.07920E6_field_r,   &   !< parameter 6   - [J/(m3*K)] specific heat capacity 2nd roof layer
       2.11200E6_field_r,   &   !< parameter 7   - [J/(m3*K)] specific heat capacity 3rd roof layer
       1.52600E6_field_r,   &   !< parameter 8   - [J/(m3*K)] specific heat capacity 4th roof layer (inside)
       0.160_field_r,       &   !< parameter 9   - [W/(m*K)] thermal conductivity 1st roof layer (outside)
       0.046_field_r,       &   !< parameter 10  - [W/(m*K)] thermal conductivity 2nd roof layer
       2.100_field_r,       &   !< parameter 11  - [W/(m*K)] thermal conductivity 3rd roof layer
       0.700_field_r,       &   !< parameter 12  - [W/(m*K)] thermal conductivity 4th roof layer (inside)
       0.15_field_r,        &   !< parameter 13  - [m] z0 roughness length for momentum
       0.10_field_r,        &   !< parameter 14  - [-] albedo
       0.95_field_r,        &   !< parameter 15  - [-] emissivity
       0.02_field_r,        &   !< parameter 16  - [m] 1st wall layer thickness (outside)
       0.06_field_r,        &   !< parameter 17  - [m] 2nd wall layer thickness
       0.24_field_r,        &   !< parameter 18  - [m] 3rd wall layer thickness
       0.02_field_r,        &   !< parameter 19  - [m] 4th wall layer thickness
       1.5200E6_field_r,    &   !< parameter 20  - [J/(m3*K)] specific heat capacity 1st wall layer (outside)
       0.0792E6_field_r,    &   !< parameter 21  - [J/(m3*K)] specific heat capacity 2nd wall layer
       2.1120E6_field_r,    &   !< parameter 22  - [J/(m3*K)] specific heat capacity 3rd wall layer
       1.5260E6_field_r,    &   !< parameter 23  - [J/(m3*K)] specific heat capacity 4th wall layer (inside)
       0.930_field_r,       &   !< parameter 24  - [W/(m*K)] thermal conductivity 1st wall layer (outside)
       0.046_field_r,       &   !< parameter 25  - [W/(m*K)] thermal conductivity 2nd wall layer
       2.100_field_r,       &   !< parameter 26  - [W/(m*K)] thermal conductivity 3rd wall layer
       0.700_field_r,       &   !< parameter 27  - [W/(m*K)] thermal conductivity 4th wall layer (inside)
       0.001_field_r,       &   !< parameter 28  - [m] z0 roughness length for momentum
       0.30_field_r,        &   !< parameter 29  - [-] albedo
       0.93_field_r,        &   !< parameter 30  - [-] emissivity
       0.02_field_r,        &   !< parameter 31  - [m] 1st window layer thickness (glass sheet + air total) (outside)
       0.02_field_r,        &   !< parameter 32  - [m] 2rd window layer thickness
       0.02_field_r,        &   !< parameter 33  - [m] 3rd window layer thickness
       0.02_field_r,        &   !< parameter 34  - [m] 4th window layer thickness (inside)
       1.736E6_field_r,     &   !< parameter 35  - [J/(m3*K)] specific heat capacity 1st window layer (outside)
       1.736E6_field_r,     &   !< parameter 36  - [J/(m3*K)] specific heat capacity 2nd window layer
       1.736E6_field_r,     &   !< parameter 37  - [J/(m3*K)] specific heat capacity 3rd window layer
       1.736E6_field_r,     &   !< parameter 38  - [J/(m3*K)] specific heat capacity 4th window layer (inside)
       0.18_field_r,        &   !< parameter 39  - [W/(m*K)] thermal conductivity 1st window layer (outside)
       0.18_field_r,        &   !< parameter 40  - [W/(m*K)] thermal conductivity 2nd window layer
       0.18_field_r,        &   !< parameter 41  - [W/(m*K)] thermal conductivity 3rd window layer
       0.18_field_r,        &   !< parameter 42  - [W/(m*K)] thermal conductivity 4th window layer (inside)
       0.65_field_r,        &   !< parameter 43  - [-] transmissivity
       0.15_field_r,        &   !< parameter 44  - [-] albedo
       0.87_field_r         &   !< parameter 45  - [-] emissivity
    /)

!
!-- Office, > 2000.
    building_pars_slurb(:,6) = (/                                                                  &
       0.29_field_r,        &   !< parameter 0   - [-] window fraction
       0.02_field_r,        &   !< parameter 1   - [m] 1st roof layer thickness (outside)
       0.04_field_r,        &   !< parameter 2   - [m] 2nd roof layer thickness
       0.30_field_r,        &   !< parameter 3   - [m] 3rd roof layer thickness
       0.02_field_r,        &   !< parameter 4   - [m] 4th roof layer thickness (inside)
       3.75360E6_field_r,   &   !< parameter 5   - [J/(m3*K)] specific heat capacity 1st roof layer (outside)
       0.70965E6_field_r,   &   !< parameter 6   - [J/(m3*K)] specific heat capacity 2nd roof layer
       0.07920E6_field_r,   &   !< parameter 7   - [J/(m3*K)] specific heat capacity 3rd roof layer
       1.52600E6_field_r,   &   !< parameter 8   - [J/(m3*K)] specific heat capacity 4th roof layer (inside)
       0.520_field_r,       &   !< parameter 9   - [W/(m*K)] thermal conductivity 1st roof layer (outside)
       0.120_field_r,       &   !< parameter 10  - [W/(m*K)] thermal conductivity 2nd roof layer
       0.035_field_r,       &   !< parameter 11  - [W/(m*K)] thermal conductivity 3rd roof layer
       0.700_field_r,       &   !< parameter 12  - [W/(m*K)] thermal conductivity 4th roof layer (inside)
       0.15_field_r,        &   !< parameter 13  - [m] z0 roughness length for momentum
       0.17_field_r,        &   !< parameter 14  - [-] albedo
       0.92_field_r,        &   !< parameter 15  - [-] emissivity
       0.02_field_r,        &   !< parameter 16  - [m] 1st wall layer thickness (outside)
       0.20_field_r,        &   !< parameter 17  - [m] 2nd wall layer thickness
       0.36_field_r,        &   !< parameter 18  - [m] 3rd wall layer thickness
       0.02_field_r,        &   !< parameter 19  - [m] 4th wall layer thickness
       1.5200E6_field_r,    &   !< parameter 20  - [J/(m3*K)] specific heat capacity 1st wall layer (outside)
       0.0792E6_field_r,    &   !< parameter 21  - [J/(m3*K)] specific heat capacity 2nd wall layer
       1.3400E6_field_r,    &   !< parameter 22  - [J/(m3*K)] specific heat capacity 3rd wall layer
       1.5260E6_field_r,    &   !< parameter 23  - [J/(m3*K)] specific heat capacity 4th wall layer (inside)
       0.930_field_r,       &   !< parameter 24  - [W/(m*K)] thermal conductivity 1st wall layer (outside)
       0.035_field_r,       &   !< parameter 25  - [W/(m*K)] thermal conductivity 2nd wall layer
       0.680_field_r,       &   !< parameter 26  - [W/(m*K)] thermal conductivity 3rd wall layer
       0.700_field_r,       &   !< parameter 27  - [W/(m*K)] thermal conductivity 4th wall layer (inside)
       0.001_field_r,       &   !< parameter 28  - [m] z0 roughness length for momentum
       0.37_field_r,        &   !< parameter 29  - [-] albedo
       0.93_field_r,        &   !< parameter 30  - [-] emissivity
       0.02_field_r,        &   !< parameter 31  - [m] 1st window layer thickness (glass sheet + air total) (outside)
       0.02_field_r,        &   !< parameter 32  - [m] 2rd window layer thickness
       0.02_field_r,        &   !< parameter 33  - [m] 3rd window layer thickness
       0.02_field_r,        &   !< parameter 34  - [m] 4th window layer thickness (inside)
       1.736E6_field_r,     &   !< parameter 35  - [J/(m3*K)] specific heat capacity 1st window layer (outside)
       1.736E6_field_r,     &   !< parameter 36  - [J/(m3*K)] specific heat capacity 2nd window layer
       1.736E6_field_r,     &   !< parameter 37  - [J/(m3*K)] specific heat capacity 3rd window layer
       1.736E6_field_r,     &   !< parameter 38  - [J/(m3*K)] specific heat capacity 4th window layer (inside)
       0.11_field_r,        &   !< parameter 39  - [W/(m*K)] thermal conductivity 1st window layer (outside)
       0.11_field_r,        &   !< parameter 40  - [W/(m*K)] thermal conductivity 2nd window layer
       0.11_field_r,        &   !< parameter 41  - [W/(m*K)] thermal conductivity 3rd window layer
       0.11_field_r,        &   !< parameter 42  - [W/(m*K)] thermal conductivity 4th window layer (inside)
       0.57_field_r,        &   !< parameter 43  - [-] transmissivity
       0.18_field_r,        &   !< parameter 44  - [-] albedo
       0.80_field_r         &   !< parameter 45  - [-] emissivity
    /)

!
!-- Asphalt concrete mix (I-II), stone aggregate(III), gravel and soil(IV), PALM-LSM default.
    pavement_pars_slurb(:,1) = (/                                                                  &
       0.01_field_r,      &   !< parameter 0   - [m] 1st pavement layer thickness (top)
       0.04_field_r,      &   !< parameter 1   - [m] 2nd pavement layer thickness
       0.20_field_r,      &   !< parameter 2   - [m] 3rd pavement layer thickness
       1.00_field_r,      &   !< parameter 3   - [m] 4th pavement layer thickness (bottom)
       2.00E6_field_r,    &   !< parameter 4   - [J/(m3*K)] heat capacity 1st pavement layer (top)
       2.00E6_field_r,    &   !< parameter 5   - [J/(m3*K)] heat capacity 2nd pavement layer
       2.00E6_field_r,    &   !< parameter 6   - [J/(m3*K)] heat capacity 3rd pavement layer
       1.40E6_field_r,    &   !< parameter 7   - [J/(m3*K)] heat capacity 4th pavement layer (bottom)
       1.00_field_r,      &   !< parameter 8   - [W/(m*K)] thermal conductivity 1st pavement layer (top)
       1.00_field_r,      &   !< parameter 9   - [W/(m*K)] thermal conductivity 2nd pavement layer
       2.10_field_r,      &   !< parameter 10  - [W/(m*K)] thermal conductivity 3rd pavement layer
       0.40_field_r,      &   !< parameter 11  - [W/(m*K)] thermal conductivity 4th pavement layer (bottom)
       5.0E-2_field_r,    &   !< parameter 12  - [m] z0 roughness length for momentum
       0.17_field_r,      &   !< parameter 13  - [-] albedo
       0.93_field_r       &   !< parameter 14  - [-] emissivity
    /)

!
!-- Asphalt concrete (I-II), stone aggregate (III), gravel and soil (IV), Masson et al. (2002).
    pavement_pars_slurb(:,2) = (/                                                                  &
       0.01_field_r,      &   !< parameter 0   - [m] 1st pavement layer thickness (top)
       0.04_field_r,      &   !< parameter 1   - [m] 2nd pavement layer thickness
       0.20_field_r,      &   !< parameter 2   - [m] 3rd pavement layer thickness
       1.00_field_r,      &   !< parameter 3   - [m] 4th pavement layer thickness (bottom)
       1.74E6_field_r,    &   !< parameter 4   - [J/(m3*K)] heat capacity 1st pavement layer (top)
       1.74E6_field_r,    &   !< parameter 5   - [J/(m3*K)] heat capacity 2nd pavement layer
       2.00E6_field_r,    &   !< parameter 6   - [J/(m3*K)] heat capacity 3rd pavement layer
       1.40E6_field_r,    &   !< parameter 7   - [J/(m3*K)] heat capacity 4th pavement layer (bottom)
       0.82_field_r,      &   !< parameter 8   - [W/(m*K)] thermal conductivity 1st pavement layer (top)
       0.82_field_r,      &   !< parameter 9   - [W/(m*K)] thermal conductivity 2nd pavement layer
       2.10_field_r,      &   !< parameter 10  - [W/(m*K)] thermal conductivity 3rd pavement layer
       0.40_field_r,      &   !< parameter 11  - [W/(m*K)] thermal conductivity 4th pavement layer (bottom)
       5.0E-2_field_r,    &   !< parameter 12  - [m] z0 roughness length for momentum
       0.10_field_r,      &   !< parameter 13  - [-] albedo
       0.95_field_r       &   !< parameter 14  - [-] emissivity
    /)

!
!-- Concrete (Portland concrete, I-II), stone aggregate (III), gravel and soil (IV),
!-- Masson et al. (2002) and Yaghoobian et al. (2009).
    pavement_pars_slurb(:,3) = (/                                                                  &
       0.01_field_r,      &   !< parameter 0   - [m] 1st pavement layer thickness (top)
       0.04_field_r,      &   !< parameter 1   - [m] 2nd pavement layer thickness
       0.20_field_r,      &   !< parameter 2   - [m] 3rd pavement layer thickness
       1.00_field_r,      &   !< parameter 3   - [m] 4th pavement layer thickness (bottom)
       2.11E6_field_r,    &   !< parameter 4   - [J/(m3*K)] heat capacity 1st pavement layer (top)
       2.11E6_field_r,    &   !< parameter 5   - [J/(m3*K)] heat capacity 2nd pavement layer
       2.00E6_field_r,    &   !< parameter 6   - [J/(m3*K)] heat capacity 3rd pavement layer
       1.40E6_field_r,    &   !< parameter 7   - [J/(m3*K)] heat capacity 4th pavement layer (bottom)
       1.51_field_r,      &   !< parameter 8   - [W/(m*K)] thermal conductivity 1st pavement layer (top)
       1.51_field_r,      &   !< parameter 9   - [W/(m*K)] thermal conductivity 2nd pavement layer
       2.10_field_r,      &   !< parameter 10  - [W/(m*K)] thermal conductivity 3rd pavement layer
       0.40_field_r,      &   !< parameter 11  - [W/(m*K)] thermal conductivity 4th pavement layer (bottom)
       5.0E-2_field_r,    &   !< parameter 12  - [m] z0 roughness length for momentum
       0.30_field_r,      &   !< parameter 13  - [-] albedo
       0.90_field_r       &   !< parameter 14  - [-] emissivity
    /)

!
!-- Sett (I-II), stone aggregate (III), gravel and soil (IV),Masson et al. (2002), Oke (1987)
!-- and Mandanici et al. (2016).
    pavement_pars_slurb(:,4) = (/                                                                  &
       0.01_field_r,      &   !< parameter 0   - [m] 1st pavement layer thickness (top)
       0.04_field_r,      &   !< parameter 1   - [m] 2nd pavement layer thickness
       0.20_field_r,      &   !< parameter 2   - [m] 3rd pavement layer thickness
       1.00_field_r,      &   !< parameter 3   - [m] 4th pavement layer thickness (bottom)
       2.25E6_field_r,    &   !< parameter 4   - [J/(m3*K)] heat capacity 1st pavement layer (top)
       2.25E6_field_r,    &   !< parameter 5   - [J/(m3*K)] heat capacity 2nd pavement layer
       2.00E6_field_r,    &   !< parameter 6   - [J/(m3*K)] heat capacity 3rd pavement layer
       1.40E6_field_r,    &   !< parameter 7   - [J/(m3*K)] heat capacity 4th pavement layer (bottom)
       2.19_field_r,      &   !< parameter 8   - [W/(m*K)] thermal conductivity 1st pavement layer (top)
       2.19_field_r,      &   !< parameter 9   - [W/(m*K)] thermal conductivity 2nd pavement layer
       2.10_field_r,      &   !< parameter 10  - [W/(m*K)] thermal conductivity 3rd pavement layer
       0.40_field_r,      &   !< parameter 11  - [W/(m*K)] thermal conductivity 4th pavement layer (bottom)
       5.0E-2_field_r,    &   !< parameter 12  - [m] z0 roughness length for momentum
       0.17_field_r,      &   !< parameter 13  - [-] albedo
       0.95_field_r       &   !< parameter 14  - [-] emissivity
    /)

!
!-- Pavement stones (I-II), stone aggregate (III), gravel and soil (IV),
!-- Masson et al. (2002), Oke (1987) and Göttsche & Hulley (2012).
    pavement_pars_slurb(:,5) = (/                                                                  &
       0.01_field_r,      &   !< parameter 0   - [m] 1st pavement layer thickness (top)
       0.04_field_r,      &   !< parameter 1   - [m] 2nd pavement layer thickness
       0.20_field_r,      &   !< parameter 2   - [m] 3rd pavement layer thickness
       1.00_field_r,      &   !< parameter 3   - [m] 4th pavement layer thickness (bottom)
       2.25E6_field_r,    &   !< parameter 4   - [J/(m3*K)] heat capacity 1st pavement layer (top)
       2.25E6_field_r,    &   !< parameter 5   - [J/(m3*K)] heat capacity 2nd pavement layer
       2.00E6_field_r,    &   !< parameter 6   - [J/(m3*K)] heat capacity 3rd pavement layer
       1.40E6_field_r,    &   !< parameter 7   - [J/(m3*K)] heat capacity 4th pavement layer (bottom)
       2.19_field_r,      &   !< parameter 8   - [W/(m*K)] thermal conductivity 1st pavement layer (top)
       2.19_field_r,      &   !< parameter 9   - [W/(m*K)] thermal conductivity 2nd pavement layer
       2.10_field_r,      &   !< parameter 10  - [W/(m*K)] thermal conductivity 3rd pavement layer
       0.40_field_r,      &   !< parameter 11  - [W/(m*K)] thermal conductivity 4th pavement layer (bottom)
       5.0E-2_field_r,    &   !< parameter 12  - [m] z0 roughness length for momentum
       0.17_field_r,      &   !< parameter 13  - [-] albedo
       0.93_field_r       &   !< parameter 14  - [-] emissivity
    /)

 END SUBROUTINE slurb_default_pars
end module modslurbdata