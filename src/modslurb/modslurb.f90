!> \file modslurb.f90
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
module modslurb
    use netcdf
    use modprecision, only : field_r
    use ieee_arithmetic, only: ieee_is_nan
    use modslurbdata
    use modslurb_energybalance, only: slurb_energy_balance_model
    use modslurb_radiationmodel, only: slurb_radiation_model
    use modslurb_resistance_stability, only: calc_canyon_resistances, calc_urban_resistances
    use modslurbhelpers, only: slurb_set_previous_timestep, slurb_read_namelist, magnus, calc_1d_heat_equation
    use modglobal, only: g => grav, cp => cp, kappa => fkar, pi
    use modfields, only: rho_air_zw => rhobf, exnf
    use modlogging, only: finish, warning

    implicit none

    character(len=*), parameter :: modname = 'modslurb'
    private
    
    public :: slurb_radiation_model
    public :: slurb_canyon_model
    public :: slurb_energy_balance_model
    public :: calc_canyon_resistances
    public :: calc_urban_resistances
    public :: slurb_urban_aggregation_model
    public :: slurb_update_external_vars
    public :: slurb_read_namelist
    public :: slurb_set_previous_timestep
    public :: fraction_slurb
    public :: slurb_tile
    public :: initslurb
    public :: exitslurb
    public :: enable_slurb
    !
    !-- Arrays for output temporal averaging.
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  albedo_urb_av         !< road liquid water coverage 
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  c_liq_road_av         !< road liquid water coverage
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  c_liq_roof_av         !< roof liquid water coverage
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  emiss_urb_av          !< road liquid water coverage
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  ghf_road_av           !< heat flux between the road bottom layer and soil
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  ghf_roof_av           !< heat flux between the roof bottom layer and indoor air
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  ghf_wall_a_av         !< heat flux between the wall a bottom layer and indoor air
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  ghf_wall_b_av         !< heat flux between the wall b bottom layer and indoor air
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  ghf_win_a_av          !< heat flux between the window a bottom layer and indoor air
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  ghf_win_b_av          !< heat flux between the window b bottom layer and indoor air
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  m_liq_road_av         !< liquid water reservoir on roads
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  m_liq_roof_av         !< liquid water reservoir on roofs
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  ol_can_av             !< street canyon top obukhov length
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  ol_road_av            !< road obukhov length
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  ol_roof_av            !< roof obukhov length
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  ol_urb_av             !< urban obukhov length for momentum flux
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  pt_can_av             !< street canyon air potential temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  pt_road_av            !< road surface potential temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  pt_roof_av            !< roof surface potential temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  pt_wall_a_av          !< wall a surface potential temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  pt_wall_b_av          !< wall b surface potential temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  pt_win_a_av           !< window a surface potential temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  pt_win_b_av           !< window b surface potential temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  q_can_av              !< street canyon water vapour mixing ratio
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  q_road_av             !< road surface mixing ratio
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  q_roof_av             !< roof surface mixing ratio
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  qs_road_av            !< road surface saturation mixing ratio
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  qs_roof_av            !< roof surface saturation mixing ratio
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  qsws_can_av           !< latent heat flux between the street canyon and the atmosphere
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  qsws_external_av      !< latent heat flux external to the model (e.g. from industry)
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  qsws_road_av          !< latent heat flux between the road and the street canyon air
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  qsws_roof_av          !< latent heat flux between the roof and the atmosphere
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  qsws_urb_av           !< urban aggregated latent heat flux
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rad_lw_net_road_av    !< road surface net longwave radiative flux
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rad_lw_net_roof_av    !< roof surface net longwave radiative flux
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rad_lw_net_urb_av     !< urban aggergated net longwave radiative flux
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rad_lw_net_wall_a_av  !< wall a surface net longwave radiative flux
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rad_lw_net_wall_b_av  !< wall b surface net longwave radiative flux
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rad_lw_net_win_a_av   !< window a surface net longwave radiative flux
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rad_lw_net_win_b_av   !< window b surface net longwave radiative flux
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rad_sw_net_road_av    !< road surface net shortwave radiative flux
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rad_sw_net_roof_av    !< roof surface net shortwave radiative flux
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rad_sw_net_urb_av     !< aggegated urban surface net shortwave radiative flux
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rad_sw_net_wall_a_av  !< wall a surface net shortwave radiative flux
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rad_sw_net_wall_b_av  !< wall b surface net shortwave radiative flux
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rad_sw_net_win_a_av   !< window a surface net shortwave radiative flux
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rad_sw_net_win_b_av   !< window b surface net shortwave radiative flux
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rad_sw_tr_win_a_av    !< window a surface transmitted shortwave radiative flux
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rad_sw_tr_win_b_av    !< window b surface transmitted shortwave radiative flux
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rah_can_av            !< street canyon aerodynamic resistance for heat
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rah_road_av           !< road aerodynamic resistance for heat
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rah_roof_av           !< roof aerodynamic resistance for heat
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rah_wall_a_av         !< wall A aerodynamic resistance for heat (DOE-2)
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rah_wall_b_av         !< wall B aerodynamic resistance for heat (DOE-2)
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rah_win_a_av          !< window A aerodynamic resistance for heat (DOE-2)
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rah_win_b_av          !< window B aerodynamic resistance for heat (DOE-2)
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rah_facade_av         !< wall and window aerodynamic resistance for heat (combined)
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  ram_urb_av            !< urban aerodynamic resistance for momentum
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rib_can_av            !< street canyon top bulk richardson number
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rib_road_av           !< road bulk richardson number
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rib_roof_av           !< roof bulk richardson number
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  shf_can_av            !< sensible heat flux between the street canyon and the atmosphere
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  shf_external_av       !< sensible heat flux external to the model (e.g. from industry)
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  shf_road_av           !< sensible heat flux between the road and the street canyon air
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  shf_roof_av           !< sensible heat flux between the roof and the atmosphere
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  shf_traffic_av        !< sensible heat flux from traffic to the canyon air
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  shf_urb_av            !< urban aggregated sensible heat flux
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  shf_wall_a_av         !< sensible heat flux between the wall a and the canyon air
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  shf_wall_b_av         !< sensible heat flux between the wall b and the canyon air
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  shf_win_a_av          !< sensible heat flux between the window a and the canyon air
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  shf_win_b_av          !< sensible heat flux between the window b and the canyon air
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  t_2m_urb_av           !< extrapolated 2-metre urban surface temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  t_c_urb_av            !< complete urban surface temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  t_can_av              !< street canyon air temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  t_h_urb_av            !< effective urban surface temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  t_rad_urb_av          !< effective urban surface radiative temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  t_surf_road_av        !< road surface temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  t_surf_roof_av        !< roof surface temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  t_surf_wall_a_av      !< wall a surface temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  t_surf_wall_b_av      !< wall b surface temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  t_surf_win_a_av       !< window a surface temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  t_surf_win_b_av       !< window b surface temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  us_can_av             !< friction velocity for street canyons
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  us_road_av            !< friction velocity for roads
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  us_roof_av            !< friction velocity for roofs
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  us_urb_av             !< urban friction velocity
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  usws_urb_av           !< urban momentum flux (u-component)
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  uv_abs_can_av         !< street canyon wind speed
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  uv_eff_can_av         !< street canyon effective wind speed
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  vpt_can_av            !< street canyon air virtual potential temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  vpt_road_av           !< road surface virtual potential temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  vpt_roof_av           !< roof surface virtual potential temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  vsws_urb_av           !< urban momentum flux (v-component
    !
    !-- Arrays for output of unmodified LSM fluxes (2D).
    ! REAL(field_r), DIMENSION(:,:,:), ALLOCATABLE ::  shf_lsm_av   !< sensible heat flux from lsm surfaces
    ! REAL(field_r), DIMENSION(:,:,:), ALLOCATABLE ::  qsws_lsm_av  !< latent heat flux from lsm surfaces

    !
    !-- Arrays for output temporal averaging (2D).
    ! REAL(field_r), DIMENSION(:,:,:), ALLOCATABLE ::  t_road_av    !< road temperature (all layers)
    ! REAL(field_r), DIMENSION(:,:,:), ALLOCATABLE ::  t_roof_av    !< roof temperature (all layers)
    ! REAL(field_r), DIMENSION(:,:,:), ALLOCATABLE ::  t_wall_a_av  !< wall a temperature (all layers)
    ! REAL(field_r), DIMENSION(:,:,:), ALLOCATABLE ::  t_wall_b_av  !< wall b temperature (all layers)
    ! REAL(field_r), DIMENSION(:,:,:), ALLOCATABLE ::  t_win_a_av   !< window a temperature (all layers)
    ! REAL(field_r), DIMENSION(:,:,:), ALLOCATABLE ::  t_win_b_av   !< window b temperature (all layers)



    

contains


subroutine initslurb
    use modglobal,   only : i1, j1
    use modlsmdata, only : ilu, nlu, tile
    use modslurbhelpers, only: slurb_bulk_allocations
    use modchecksim, only: check_array
    implicit none
    integer, parameter :: rkind = kind( 1.0_field_r )

    integer :: i, j, slurb_ilu
    character(len=*), parameter :: routine = modname//'/initslurb'
    !-- Initialize bounds for subsurface layers.
    nzt_road = 1
    nzb_road = n_layers_roads
    nzt_roof = 1
    nzb_roof = n_layers_roofs
    nzt_wall = 1
    nzb_wall = n_layers_walls
    nzt_win  = 1
    nzb_win  = n_layers_windows

    

    do ilu=1,nlu
        if (tile(ilu)%lushort == "slb") then
            ! make sure only one tile is slb, enable_slurb is false by default, set to true if any tile found.
            if (enable_slurb) then
                call finish(routine, "Only one SLUrb tile at a time supported. Make sure there's only one tile with lushort=slb.")
            endif
            enable_slurb = .true.
            slurb_ilu = ilu
        endif
    end do
    if (.not. enable_slurb) then
        return
    end if
    call warning(routine, "SLUrb module enabled. Keep in mind that calculation of effective albedo is not implemented yet!")
    call warning(routine, "SLUrb module enabled. Keep in mind that different building drag parametrizations have not been tested yet!")
    call warning(routine, "SLUrb module enabled. Keep in mind that moist_physics=false has not been tested yet!")
    call warning(routine, "SLUrb module enabled. Only rrtmgp radiation has been tested with SLUrb!")

    call check_array([deep_soil_temperature],"deep_soil_temperature", routine, &
    threshold=[real(100.0_field_r,rkind),real( 400.0_field_r,rkind)], &
    stop_if_invalid=.true.)
    call check_array([building_indoor_temperature],"building_indoor_temperature", &
    routine, threshold=[real(100.0_field_r,rkind),real( 400.0_field_r,rkind)], &
    stop_if_invalid=.true.)

    call slurb_bulk_allocations
    do j=2,j1
        do i=2,i1
            fraction_slurb(i,j) = tile(slurb_ilu)%frac(i,j) ! (.)
        enddo
    enddo


    slurb_tile%dt_max(:,:) = HUGE( 1.0_field_r ) ! (s)

    

    call process_surface_parameters

    call precompute_latent_variables

    call init_slurb_variables

    ! call radiation
    
    



end subroutine initslurb

subroutine exitslurb
    use modslurbhelpers, only: slurb_bulk_deallocations
    implicit none
    if (.not. enable_slurb) then
        return
    end if
    call slurb_bulk_deallocations

end subroutine exitslurb

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Updates model external variables, e.g. the variables defined at the first atmospheric level based
    !> on atmospheric simulation state as well as the temporally dynamic SLUrb input variables.
    !--------------------------------------------------------------------------------------------------!
 SUBROUTINE slurb_update_external_vars

    use modglobal, only : cp, rlv, cu, cv, i1, j1
    use modfields, only : ql0, u0, v0, qt0, exnf, thl0
    implicit none
    INTEGER ::  i      !< loop index
    INTEGER::  j      !< loop index
    INTEGER ::  k_atm  !< k index of the first atmospheric level

    ! REAL(field_r) ::  fac_dt  !< factor for linear interpolation between timesteps
    REAL(field_r) ::  vtws    !< buoyancy flux (m K s^-1)
    REAL(field_r) ::  ws      !< free convection velocity scale

    real(field_r) :: du, dv

    real :: rho_cp !< cp * rho (J m^-3 K^-1)
    rho_cp = cp * rho_air_zw(1) !TODOSELF


    do j=2,j1
      do i=2,i1
        k_atm = 1

    !        k_atm = topo_top_ind(j,i,0) + 1


            ! K = K + J/kg /(J/kg K^-1) * (kg/kg)
            slurb_tile%pt1(i,j)  = thl0(i, j, k_atm) + (rlv/(cp * exnf(k_atm)))  * ql0(i,j,k_atm)
            slurb_tile%q1(i,j)   = qt0(i, j, k_atm) - ql0(i, j, k_atm) !TODOSELF BUG
            slurb_tile%vpt1(i,j) = slurb_tile%pt1(i,j) * ( 1.0_field_r + 0.61_field_r * slurb_tile%q1(i,j) )

            du = 0.5*(u0(i,j,1) + u0(i+1,j,1)) + cu
            dv = 0.5*(v0(i,j,1) + v0(i,j+1,1)) + cv
            slurb_tile%uv_abs1(i,j) = sqrt(du**2 + dv**2)
            ! slurb_tile%uv_abs1(i,j) = max(0.1, sqrt(du**2 + dv**2)) DALES VERSION




            !--    Calculate surface-parallel absolute velocity uv_eff1 at cell center using
            !--    free convection scale (w_star, for unstable cases).
            ! m K s^-1 = (J^-1 kg K)(kg^-1 m^3) J s^-1 m^-2 + (J/kg)/(J/kg/K)*W/m^2
            ! m K s^-1 = (J^-1 kg K)(kg^-1 m^3) J s^-1 m^-2 + (J kg^-1) ( J K^-1 kg^-1)^-1 W/m^2
            ! m K s^-1 = (J^-1 kg K)(kg^-1 m^3) J s^-1 m^-2 + (J kg^-1) J^-1 K kg W/m^2
            ! m K s^-1 = (J^-1 kg K)(kg^-1 m^3) J s^-1 m^-2 +  K W m^-2

            ! m K s^-1 = (J^-1 kg K)(kg^-1 m^3) J s^-1 m^-2 + (J/kg)/(J kg^-1 K^-1 kg m^-3)*W/m^2
            ! m K s^-1 = (J^-1 kg K)(kg^-1 m^3) J s^-1 m^-2 + (J K^-1 m^-3)^-1 W/m^2
            ! m K s^-1 = (J^-1 kg K)(kg^-1 m^3) J s^-1 m^-2 + J^-1 K W m^-2 m^3
            ! m K s^-1 = (J^-1 kg K)(kg^-1 m^3) J s^-1 m^-2 + K s^-1 m^-2 m^3
            ! m K s^-1 = (J^-1 kg K)(kg^-1 m^3) J s^-1 m^-2 + K s^-1 m
            vtws = (1/(rho_cp)) * slurb_tile%shf_urb(i,j) + (1 / rho_cp) * slurb_tile%qsws_urb(i,j) ! (m K s^-1)
            !
            !--    No scaling for stable cases:
            vtws = MERGE( vtws, 0.0_field_r, vtws > 0.0_field_r )
            ! m/s = (m s^-2 K^-1 m ?)^(1/3)
            ! m/s = (m^2 s^-2 K^-1 (m K s^-1))^(1/3)
            ! m/s = (m^3 s^-3)^(1/3)
            ws = ( g / slurb_tile%pt1(i,j) * slurb_tile%z_mo(i,j) * vtws )**( 1.0_field_r / 3.0_field_r )  ! (m s^-1)

            slurb_tile%uv_eff1(i, j) = sqrt(du**2 + dv**2 + ws**2)
        enddo
    enddo

    ! IF ( slurb_dynamic%ntime > 0 )  CALL update_dynamic_inputs

end subroutine slurb_update_external_vars




 !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Process parameters dependent on the building/pavement type and properties.
    !--------------------------------------------------------------------------------------------------!
 SUBROUTINE process_surface_parameters
    use modglobal, only: i1, j1, i2, j2, cexpnr, zf, imax, jmax, itot, jtot
    use modmpi,      only : myid, myidx, myidy
    use modslurbdata, only : slurb_default_pars, building_pars_slurb, pavement_pars_slurb
    use modstat_nc, only : read_nc_field, nchandle_error
    use modchecksim, only: check_array
    implicit none
    character(len=*), parameter :: routine = modname//'/process_surface_parameters'

    real(field_r), allocatable :: canyon_orientation_tmp(:,:)
    integer, parameter :: rkind = kind( 1.0_field_r )
    INTEGER, DIMENSION(:,:), ALLOCATABLE ::  type_tmp  !< array to contain building type temporarily
    integer i,j,k, ncid
    
    if (lread_from_netcdf) then
        call nchandle_error(nf90_open('inslurb.'//cexpnr//'.nc', NF90_NOWRITE, ncid))
    endif



    !
    !-- Process variables related to urban form.
    !
    !-- Internally, f_bld refers to building plan area fraction of the urban surface. However, it is
    !-- more common to report the building plan area fraction as a fraction of the total surface, e.g.
    !-- in the case of LCZs. We want the user input correspond to the latter, thus the scaling.
    slurb_tile%f_bld(:,:) = -9999.0_field_r
    call read_nc_field(ncid, 'f_bld', slurb_tile%f_bld(2:i1,2:j1), fillvalue=0.5_field_r, requirefill=.true., &
                                start = (/1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/imax, jmax/) )
    call check_array(slurb_tile%f_bld(2:i1,2:j1), "f_bld", routine, threshold=[real( TINY( 1.0_field_r ),rkind),real(  1.0_field_r ,rkind)], stop_if_invalid=.true.)
    do j=2,j1
      do i=2,i1
        if ( fraction_slurb(i,j) /= 0) then
            slurb_tile%f_bld(i, j) = slurb_tile%f_bld(i, j) / fraction_slurb(i, j)
        endif
        enddo
    enddo

    slurb_tile%f_bld_frn(:,:) = -9999.0_field_r

    call read_nc_field(ncid, 'f_bld_frn', slurb_tile%f_bld_frn(2:i1,2:j1), fillvalue=0.2_field_r, requirefill=.true.,&
                                start = (/1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/imax, jmax/) )


    slurb_tile%h_bld(:,:) = -9999.0_field_r
    call read_nc_field(ncid, 'h_bld', slurb_tile%h_bld(2:i1,2:j1), fillvalue=building_height, requirefill=.true.,&
                                start = (/1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/imax, jmax/) )
    call check_array(slurb_tile%h_bld(2:i1,2:j1), "h_bld", routine, threshold=[real( 0.1_field_r,rkind),real(  1000.0_field_r ,rkind)], stop_if_invalid=.true.)

    !
    !-- Urban surface and street canyon MOST heights.
    do j=2,j1
      do i=2,i1
       slurb_tile%z_mo(i,j) = 0.5_field_r * (zf(2) - zf(1)) ! (m)
    !    slurb_tile%z_mo(i,j) = 0.5_field_r *  dzw(topo_top_ind(j,i,0)+1)
       slurb_tile%z_mo_can(i,j) = 0.5_field_r * slurb_tile%h_bld(i,j) ! (m)
      enddo
    enddo
    !
    !-- Process canyon direction information if anisotropic street canyons are enabled.
    IF ( anisotropic_street_canyons )  THEN
       ALLOCATE( canyon_orientation_tmp(i2,j2) )
       canyon_orientation_tmp(:,:) = -9999.0_field_r
        IF ( street_canyon_orientation /= -9999.0_field_r  )  THEN
            canyon_orientation_tmp(:,:) = street_canyon_orientation
        ENDIF

        call read_nc_field(ncid, 'street_canyon_orientation', canyon_orientation_tmp(2:i1,2:j1), requirefill=.false.,&
                                start = (/1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/imax, jmax/) )
    !
    !--    In order to make it possible to have a mix of isotropic and anisotropic tiles, use
    !--    anisotropic canyons only in the case a canyon orientation has been given either for all tiles
    !--    in the namelist or per-patch basis in input file.
        do j=2,j1
            do i=2,i1
                IF ( canyon_orientation_tmp(i,j) /= -9999.0_field_r )  THEN
                    slurb_tile%anisotropic_canyon(i,j) = .TRUE.
                    slurb_tile%theta_can(i,j) = canyon_orientation_tmp(i,j) * ( pi / 180.0_field_r )
                ELSE
                    slurb_tile%anisotropic_canyon(i,j) = .FALSE.
                    slurb_tile%theta_can(i,j) = -9999.0_field_r
                ENDIF
            enddo
        enddo
       DEALLOCATE( canyon_orientation_tmp )
    ELSE
        do j=2,j1
            do i=2,i1
                slurb_tile%anisotropic_canyon(i,j) = .FALSE.
                slurb_tile%theta_can(i,j) = -9999.0_field_r
            enddo
        enddo
    ENDIF

    call read_nc_field(ncid, 'hw_can', slurb_tile%hw_can(2:i1,2:j1), fillvalue=street_canyon_aspect_ratio, requirefill=.true.,&
                                start = (/1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/imax, jmax/) )
    call check_array(slurb_tile%hw_can(2:i1,2:j1), "hw_can", routine, threshold=[real( TINY(0.0_field_r),rkind),real(  1000.0_field_r ,rkind)], stop_if_invalid=.true.) !TODOSELF REALISTIC FALUE

    call read_nc_field(ncid, 'z0_urb', slurb_tile%z0_urb(2:i1,2:j1), fillvalue=urban_roughness_length, requirefill=.true.,&
                                start = (/1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/imax, jmax/) )
    call check_array(slurb_tile%z0_urb(2:i1,2:j1), "z0_urb", routine, threshold=[real( TINY(0.0_field_r),rkind),real(  1000.0_field_r ,rkind)], stop_if_invalid=.true.) !TODOSELF REALISTIC FALUE


    ALLOCATE( type_tmp(i2, j2) )
    type_tmp(:,:) = 2 ! init building type at 2

    CALL slurb_default_pars

    call read_nc_field(ncid, "building_type", type_tmp, fillvalue=2, requirefill=.true.,&
                                start = (/1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/imax, jmax/) )
    call check_array(type_tmp, "building_type", routine, threshold=[1,6], stop_if_invalid=.true.)

    do j=2,j1
      do i=2,i1
        slurb_tile%f_win(i,j) = building_pars_slurb(0,type_tmp(i,j))

        IF ( n_layers_roofs == 4 )  THEN
            slurb_tile%dz_roof(1,i,j) = building_pars_slurb(1,type_tmp(i,j))
            slurb_tile%dz_roof(2,i,j) = building_pars_slurb(2,type_tmp(i,j))
            slurb_tile%dz_roof(3,i,j) = building_pars_slurb(3,type_tmp(i,j))
            slurb_tile%dz_roof(4,i,j) = building_pars_slurb(4,type_tmp(i,j))

            slurb_tile%c_roof(1,i,j) = building_pars_slurb(5,type_tmp(i,j))
            slurb_tile%c_roof(2,i,j) = building_pars_slurb(6,type_tmp(i,j))
            slurb_tile%c_roof(3,i,j) = building_pars_slurb(7,type_tmp(i,j))
            slurb_tile%c_roof(4,i,j) = building_pars_slurb(8,type_tmp(i,j))

            slurb_tile%lambda_roof(1,i,j) = building_pars_slurb(9,type_tmp(i,j))
            slurb_tile%lambda_roof(2,i,j) = building_pars_slurb(10,type_tmp(i,j))
            slurb_tile%lambda_roof(3,i,j) = building_pars_slurb(11,type_tmp(i,j))
            slurb_tile%lambda_roof(4,i,j) = building_pars_slurb(12,type_tmp(i,j))
        ENDIF

        slurb_tile%z0_roof(i,j)     = building_pars_slurb(13,type_tmp(i,j))
        slurb_tile%z0h_roof(i,j)    = building_pars_slurb(13,type_tmp(i,j)) * 1.0E-2
        slurb_tile%albedo_roof(i,j) = building_pars_slurb(14,type_tmp(i,j))
        slurb_tile%emiss_roof(i,j)  = building_pars_slurb(15,type_tmp(i,j))

        IF ( n_layers_walls == 4 )  THEN
            slurb_tile%dz_wall(1,i,j) = building_pars_slurb(16,type_tmp(i,j))
            slurb_tile%dz_wall(2,i,j) = building_pars_slurb(17,type_tmp(i,j))
            slurb_tile%dz_wall(3,i,j) = building_pars_slurb(18,type_tmp(i,j))
            slurb_tile%dz_wall(4,i,j) = building_pars_slurb(19,type_tmp(i,j))

            slurb_tile%c_wall(1,i,j) = building_pars_slurb(20,type_tmp(i,j))
            slurb_tile%c_wall(2,i,j) = building_pars_slurb(21,type_tmp(i,j))
            slurb_tile%c_wall(3,i,j) = building_pars_slurb(22,type_tmp(i,j))
            slurb_tile%c_wall(4,i,j) = building_pars_slurb(23,type_tmp(i,j))

            slurb_tile%lambda_wall(1,i,j) = building_pars_slurb(24,type_tmp(i,j))
            slurb_tile%lambda_wall(2,i,j) = building_pars_slurb(25,type_tmp(i,j))
            slurb_tile%lambda_wall(3,i,j) = building_pars_slurb(26,type_tmp(i,j))
            slurb_tile%lambda_wall(4,i,j) = building_pars_slurb(27,type_tmp(i,j))
        ENDIF

        slurb_tile%z0_wall(i,j)     = building_pars_slurb(28,type_tmp(i,j))
        slurb_tile%albedo_wall(i,j) = building_pars_slurb(29,type_tmp(i,j))
        slurb_tile%emiss_wall(i,j)  = building_pars_slurb(30,type_tmp(i,j))

        IF ( n_layers_windows == 4 )  THEN
            slurb_tile%dz_win(1,i,j) = building_pars_slurb(31,type_tmp(i,j))
            slurb_tile%dz_win(2,i,j) = building_pars_slurb(32,type_tmp(i,j))
            slurb_tile%dz_win(3,i,j) = building_pars_slurb(33,type_tmp(i,j))
            slurb_tile%dz_win(4,i,j) = building_pars_slurb(34,type_tmp(i,j))

            slurb_tile%c_win(1,i,j) = building_pars_slurb(35,type_tmp(i,j))
            slurb_tile%c_win(2,i,j) = building_pars_slurb(36,type_tmp(i,j))
            slurb_tile%c_win(3,i,j) = building_pars_slurb(37,type_tmp(i,j))
            slurb_tile%c_win(4,i,j) = building_pars_slurb(38,type_tmp(i,j))

            slurb_tile%lambda_win(1,i,j) = building_pars_slurb(39,type_tmp(i,j))
            slurb_tile%lambda_win(2,i,j) = building_pars_slurb(40,type_tmp(i,j))
            slurb_tile%lambda_win(3,i,j) = building_pars_slurb(41,type_tmp(i,j))
            slurb_tile%lambda_win(4,i,j) = building_pars_slurb(42,type_tmp(i,j))
        ENDIF

        slurb_tile%transmissivity_win(i,j) = building_pars_slurb(43,type_tmp(i,j))

        slurb_tile%albedo_win(i,j) = building_pars_slurb(44,type_tmp(i,j))
        slurb_tile%emiss_win(i,j)  = building_pars_slurb(45,type_tmp(i,j))

        enddo
    enddo
    if (lread_from_netcdf) then
        call read_nc_field(ncid, "pavement_type", type_tmp, fillvalue=2, requirefill=.true.,&
                                start = (/1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/imax, jmax/) )
        call check_array(type_tmp, "pavement_type", routine, threshold=[1,5], stop_if_invalid=.true.)
    endif
    !
    !-- Process pavement type.

    do j=2,j1
      do i=2,i1
        IF ( n_layers_roads == 4 )  THEN
            slurb_tile%dz_road(1,i,j) = pavement_pars_slurb(0,type_tmp(i,j))
            slurb_tile%dz_road(2,i,j) = pavement_pars_slurb(1,type_tmp(i,j))
            slurb_tile%dz_road(3,i,j) = pavement_pars_slurb(2,type_tmp(i,j))
            slurb_tile%dz_road(4,i,j) = pavement_pars_slurb(3,type_tmp(i,j))

            slurb_tile%c_road(1,i,j) = pavement_pars_slurb(4,type_tmp(i,j))
            slurb_tile%c_road(2,i,j) = pavement_pars_slurb(5,type_tmp(i,j))
            slurb_tile%c_road(3,i,j) = pavement_pars_slurb(6,type_tmp(i,j))
            slurb_tile%c_road(4,i,j) = pavement_pars_slurb(7,type_tmp(i,j))

            slurb_tile%lambda_road(1,i,j) = pavement_pars_slurb(8,type_tmp(i,j))
            slurb_tile%lambda_road(2,i,j) = pavement_pars_slurb(9,type_tmp(i,j))
            slurb_tile%lambda_road(3,i,j) = pavement_pars_slurb(10,type_tmp(i,j))
            slurb_tile%lambda_road(4,i,j) = pavement_pars_slurb(11,type_tmp(i,j))

            slurb_tile%z0_road(i,j)     = pavement_pars_slurb(12,type_tmp(i,j))
            slurb_tile%z0h_road(i,j)    = pavement_pars_slurb(12,type_tmp(i,j)) * 1.0E-2
            slurb_tile%albedo_road(i,j) = pavement_pars_slurb(13,type_tmp(i,j))
            slurb_tile%emiss_road(i,j)  = pavement_pars_slurb(14,type_tmp(i,j))
        ENDIF
        enddo
    enddo
    DEALLOCATE( type_tmp )




    !
    !-- Process material layer information such as thickness, heat capacities, if given.
    !-- By default, use information provided on building type.
    if (lread_from_netcdf) then
        call read_nc_field(ncid, 'albedo_roof', slurb_tile%albedo_roof(2:i1,2:j1), requirefill=.false.,&
                                start = (/1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/imax, jmax/) )
        call check_array(slurb_tile%albedo_roof(2:i1,2:j1), "albedo_roof", routine, threshold=[real( 0.0_field_r,rkind),real(  1.0_field_r ,rkind)], stop_if_invalid=.true.)

        call read_nc_field(ncid, 'dz_roof', slurb_tile%dz_roof(:,2:i1,2:j1), requirefill=.false.,&
                                start = (/1, 1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/n_layers_roofs, imax, jmax/) )
        call check_array(slurb_tile%dz_roof(:,2:i1,2:j1), "dz_roof", routine, threshold=[real( TINY( 1.0_field_r ),rkind),real(  HUGE( 1.0_field_r ) ,rkind)], stop_if_invalid=.true.)

        call read_nc_field(ncid, 'emiss_roof', slurb_tile%emiss_roof(2:i1,2:j1), requirefill=.false.,&
                                start = (/1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/imax, jmax/) )
        call check_array(slurb_tile%emiss_roof(2:i1,2:j1), "emiss_roof", routine, threshold=[real( 0.0_field_r,rkind),real(  1.0_field_r ,rkind)], stop_if_invalid=.true.)

        call read_nc_field(ncid, 'c_roof', slurb_tile%c_roof(:,2:i1,2:j1), requirefill=.false.,&
                                start = (/1, 1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/n_layers_roofs, imax, jmax/) )
        call check_array(slurb_tile%c_roof(:,2:i1,2:j1), "c_roof", routine, threshold=[real( TINY( 1.0_field_r ),rkind),real(  HUGE( 1.0_field_r ) ,rkind)], stop_if_invalid=.true.)

        call read_nc_field(ncid, 'z0_roof', slurb_tile%z0_roof(2:i1,2:j1), requirefill=.false.,&
                                start = (/1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/imax, jmax/) )
        call check_array(slurb_tile%z0_roof(2:i1,2:j1), "z0_roof", routine, threshold=[real( TINY( 1.0_field_r ),rkind),real(0.5_field_r * MINVAL( slurb_tile%z_mo(2:i1,2:j1)),rkind) ], stop_if_invalid=.true.)

        call read_nc_field(ncid, 'z0h_roof', slurb_tile%z0h_roof(2:i1,2:j1), requirefill=.false.,&
                                start = (/1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/imax, jmax/) )
        call check_array(slurb_tile%z0h_roof(2:i1,2:j1), "z0h_roof", routine, threshold=[real( TINY( 1.0_field_r ),rkind),real(0.5_field_r * MINVAL( slurb_tile%z_mo(2:i1,2:j1)),rkind) ], stop_if_invalid=.true.)

        call read_nc_field(ncid, 'lambda_roof', slurb_tile%lambda_roof(:,2:i1,2:j1), requirefill=.false.,&
                                start = (/1, 1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/n_layers_roofs, imax, jmax/) )
        call check_array(slurb_tile%lambda_roof(:,2:i1,2:j1), "lambda_roof", routine, threshold=[real( TINY( 1.0_field_r ),rkind),real( HUGE( 1.0_field_r ) ,rkind)], stop_if_invalid=.true.)

        call read_nc_field(ncid, 'albedo_wall', slurb_tile%albedo_wall(2:i1,2:j1), requirefill=.false.,&
                                start = (/1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/imax, jmax/) )
        call check_array(slurb_tile%albedo_wall(2:i1,2:j1), "albedo_wall", routine, threshold=[real( 0.0_field_r,rkind),real(  1.0_field_r ,rkind)], stop_if_invalid=.true.)

        call read_nc_field(ncid, 'dz_wall', slurb_tile%dz_wall(:,2:i1,2:j1), requirefill=.false.,&
                                start = (/1, 1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/n_layers_walls, imax, jmax/) )
        call check_array(slurb_tile%dz_wall(:,2:i1,2:j1), "dz_wall", routine, threshold=[real( TINY( 1.0_field_r ),rkind),real(  HUGE( 1.0_field_r ) ,rkind)], stop_if_invalid=.true.)

        call read_nc_field(ncid, 'emiss_wall', slurb_tile%emiss_wall(2:i1,2:j1), requirefill=.false.,&
                                start = (/1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/imax, jmax/) )
        call check_array(slurb_tile%emiss_wall(2:i1,2:j1), "emiss_wall", routine, threshold=[real( 0.0_field_r,rkind),real(  1.0_field_r ,rkind)], stop_if_invalid=.true.)

        call read_nc_field(ncid, 'c_wall', slurb_tile%c_wall(:,2:i1,2:j1), requirefill=.false.,&
                                start = (/1, 1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/n_layers_walls, imax, jmax/) )
        call check_array(slurb_tile%c_wall(:,2:i1,2:j1), "c_wall", routine, threshold=[real( TINY( 1.0_field_r ),rkind),real(  HUGE( 1.0_field_r ) ,rkind)], stop_if_invalid=.true.)

        call read_nc_field(ncid, 'z0_wall', slurb_tile%z0_wall(2:i1,2:j1), requirefill=.false.,&
                                start = (/1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/imax, jmax/) )
        call check_array(slurb_tile%z0_wall(2:i1,2:j1), "z0_wall", routine, threshold=[real( TINY( 1.0_field_r ),rkind),real(  1.0_field_r ,rkind)], stop_if_invalid=.true.)

        call read_nc_field(ncid, 'lambda_wall', slurb_tile%lambda_wall(:,2:i1,2:j1), requirefill=.false.,&
                                start = (/1, 1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/n_layers_walls, imax, jmax/) )
        call check_array(slurb_tile%lambda_wall(:,2:i1,2:j1), "lambda_wall", routine, threshold=[real( TINY( 1.0_field_r ),rkind),real( HUGE( 1.0_field_r ) ,rkind)], stop_if_invalid=.true.)

        call read_nc_field(ncid, 'albedo_win', slurb_tile%albedo_win(2:i1,2:j1), requirefill=.false.,&
                                start = (/1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/imax, jmax/) )
        call check_array(slurb_tile%albedo_win(2:i1,2:j1), "albedo_win", routine, threshold=[real( 0.0_field_r,rkind),real(  1.0_field_r ,rkind)], stop_if_invalid=.true.)
        
        call read_nc_field(ncid, 'dz_win', slurb_tile%dz_win(:,2:i1,2:j1), requirefill=.false.,&
                                start = (/1, 1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/n_layers_windows, imax, jmax/) )
        call check_array(slurb_tile%dz_win(:,2:i1,2:j1), "dz_win", routine, threshold=[real( TINY( 1.0_field_r ),rkind),real(  HUGE( 1.0_field_r ) ,rkind)], stop_if_invalid=.true.)

        call read_nc_field(ncid, 'emiss_win', slurb_tile%emiss_win(2:i1,2:j1), requirefill=.false.,&
                                start = (/1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/imax, jmax/) )
        call check_array(slurb_tile%emiss_win(2:i1,2:j1), "emiss_win", routine, threshold=[real( 0.0_field_r,rkind),real(  1.0_field_r ,rkind)], stop_if_invalid=.true.)

        call read_nc_field(ncid, 'c_win', slurb_tile%c_win(:,2:i1,2:j1), requirefill=.false.,&
                                start = (/1, 1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/n_layers_windows, imax, jmax/) )
        call check_array(slurb_tile%c_win(:,2:i1,2:j1), "c_win", routine, threshold=[real( TINY( 1.0_field_r ),rkind),real(  HUGE( 1.0_field_r ) ,rkind)], stop_if_invalid=.true.)

        ! slurb_tile%c_win = slurb_tile%c_win * slurb_tile%dz_win
        call read_nc_field(ncid, 'lambda_win', slurb_tile%lambda_win(:,2:i1,2:j1), requirefill=.false.,&
                                start = (/1, 1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/n_layers_windows, imax, jmax/) )
        call check_array(slurb_tile%lambda_win(:,2:i1,2:j1), "lambda_win", routine, threshold=[real( TINY( 1.0_field_r ),rkind),real( HUGE( 1.0_field_r ) ,rkind)], stop_if_invalid=.true.)

        call read_nc_field(ncid, 'transmissivity_win', slurb_tile%transmissivity_win(2:i1,2:j1), requirefill=.false.,&
                                start = (/1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/imax, jmax/) )
        call check_array(slurb_tile%transmissivity_win(2:i1,2:j1), "transmissivity_win", routine, threshold=[real( 0.0_field_r,rkind),real( 1.0_field_r ,rkind)], stop_if_invalid=.true.)

        call read_nc_field(ncid, 'f_win', slurb_tile%f_win(2:i1,2:j1), requirefill=.false.,&
                                start = (/1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/imax, jmax/) )
        call check_array(slurb_tile%f_win(2:i1,2:j1), "f_win", routine, threshold=[real( 0.0_field_r,rkind),real(  1.0_field_r ,rkind)], stop_if_invalid=.true.)


        call read_nc_field(ncid, 'albedo_road', slurb_tile%albedo_road(2:i1,2:j1), requirefill=.false.,&
                                start = (/1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/imax, jmax/) )
        call check_array(slurb_tile%albedo_road(2:i1,2:j1), "albedo_road", routine, threshold=[real( 0.0_field_r,rkind),real(  1.0_field_r ,rkind)], stop_if_invalid=.true.)

        call read_nc_field(ncid, 'dz_road', slurb_tile%dz_road(:,2:i1,2:j1), requirefill=.false.,&
                                start = (/1, 1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/n_layers_roads, imax, jmax/) )
        call check_array(slurb_tile%dz_road(:,2:i1,2:j1), "dz_road", routine, threshold=[real( TINY( 1.0_field_r ),rkind),real(  HUGE( 1.0_field_r ) ,rkind)], stop_if_invalid=.true.)

        call read_nc_field(ncid, 'emiss_road', slurb_tile%emiss_road(2:i1,2:j1), requirefill=.false.,&
                                start = (/1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/imax, jmax/) )
        call check_array(slurb_tile%emiss_road(2:i1,2:j1), "emiss_road", routine, threshold=[real( 0.0_field_r,rkind),real(  1.0_field_r ,rkind)], stop_if_invalid=.true.)

        call read_nc_field(ncid, 'c_road', slurb_tile%c_road(:,2:i1,2:j1), requirefill=.false.,&
                                start = (/1, 1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/n_layers_roads, imax, jmax/) )
        call check_array(slurb_tile%c_road(:,2:i1,2:j1), "c_road", routine, threshold=[real( TINY( 1.0_field_r ),rkind),real(  HUGE( 1.0_field_r ) ,rkind)], stop_if_invalid=.true.)

        call read_nc_field(ncid, 'z0_road', slurb_tile%z0_road(2:i1,2:j1), requirefill=.false.,&
                                start = (/1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/imax, jmax/) )
        call check_array(slurb_tile%z0_road(2:i1,2:j1), "z0_road", routine, threshold=[real( TINY( 1.0_field_r ),rkind),real(  1.0_field_r ,rkind)], stop_if_invalid=.true.)

        call read_nc_field(ncid, 'z0h_road', slurb_tile%z0h_road(2:i1,2:j1), requirefill=.false.,&
                                start = (/1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/imax, jmax/) )
        call check_array(slurb_tile%z0h_road(2:i1,2:j1), "z0h_road", routine, threshold=[real( TINY( 1.0_field_r ),rkind),real(  1.0_field_r ,rkind)], stop_if_invalid=.true.)

        call read_nc_field(ncid, 'lambda_road', slurb_tile%lambda_road(:,2:i1,2:j1), requirefill=.false.,&
                                start = (/1, 1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/n_layers_roads, imax, jmax/) )
        call check_array(slurb_tile%lambda_road(:,2:i1,2:j1), "lambda_road", routine, threshold=[real( TINY( 1.0_field_r ),rkind),real( HUGE( 1.0_field_r ) ,rkind)], stop_if_invalid=.true.)

    endif

        !-- check if z_mo and z_mo_can are realistic compared to z0 and z0h
    do j=2,j1
      do i=2,i1
        if (slurb_tile%z0_road(i,j) >= slurb_tile%z_mo_can(i,j)) then
            call finish(routine, 'i,j=',i,',',j,' z0_road=', slurb_tile%z0_road(i,j), ' z_mo_can=', slurb_tile%z_mo_can(i,j))
        endif
      enddo
    enddo

    !
    !-- SLUrb uses the total layer heat capacity instead of specific heat capacity,
    !-- so multiply c_roof by dz_roof.

    do j=2,j1
      do i=2,i1
        do k=1,4
            slurb_tile%c_roof(k,i,j) = slurb_tile%c_roof(k,i,j) * slurb_tile%dz_roof(k,i,j)
        enddo
      enddo
    enddo
    !
    !-- SLUrb uses the total layer heat capacity instead of specific heat capacity,
    !-- so multiply c_wall by dz_wall.
    do j=2,j1
      do i=2,i1
        do k=1,4
            slurb_tile%c_wall(k,i,j) = slurb_tile%c_wall(k,i,j) * slurb_tile%dz_wall(k,i,j)
        enddo
      enddo
    enddo
    !
    !-- SLUrb uses the total layer heat capacity instead of specific heat capacity, so multiply c_win
    !-- by dz_win.
    do j=2,j1
      do i=2,i1
        do k=1,4
            slurb_tile%c_win(k,i,j) = slurb_tile%c_win(k,i,j) * slurb_tile%dz_win(k,i,j)
        enddo
      enddo
    enddo
    !
    !-- SLUrb uses the total layer heat capacity instead of specific heat capacity,
    !-- so multiply c_road by dz_road.

    do j=2,j1
      do i=2,i1
        do k=1,4
            slurb_tile%c_road(k,i,j) = slurb_tile%c_road(k,i,j) * slurb_tile%dz_road(k,i,j)
        enddo
      enddo
    enddo

    !
    !-- Compute weighted wall-window albedo.
    do j=2,j1
      do i=2,i1
       slurb_tile%albedo_wall_win(i,j) = ( 1.0_field_r - slurb_tile%f_win(i,j) ) * slurb_tile%albedo_wall(i,j) +                &
                                 slurb_tile%f_win(i,j) * slurb_tile%albedo_win(i,j)
      enddo
    enddo
    !
    !-- Compute the cumulative layer thickness zw for windows.
    do j=2,j1
        do i=2,i1
            slurb_tile%zw_win(nzt_win,i,j) = slurb_tile%dz_win(nzt_win,i,j)
            DO  k = nzt_win+1, nzb_win
                slurb_tile%zw_win(k,i,j) = slurb_tile%zw_win(k-1,i,j) + slurb_tile%dz_win(k,i,j)
            enddo
        enddo
    enddo

 END SUBROUTINE process_surface_parameters



    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Initializes SLUrb model variables.
    !--------------------------------------------------------------------------------------------------!
 SUBROUTINE init_slurb_variables
    use modfields, only : thl0, ql0, qt0, u0, v0, exnf
    use modglobal, only : cp, rlv, cu, cv, i1, j1, ep
    use modsurface, only : ps
    REAL(field_r) ::  bc_atm  !< initial atmospheric boundary condition for temperature
    REAL(field_r) ::  e_s     !< initial water vapor saturation pressure
    real du,dv
    integer i,j,k_atm,k_topo

    do j=2,j1
        do i=2,i1

        k_atm = 1 !TODOSELF check vertical levels for density calculations
        k_topo = 1

        ! TODOSELF is ql0 the correct liquid water?
        ! in PALM pt=liquid water potential temperature
        ! this implies slurb_tile%pt1 = pt+L/cpexn ql0, slurb_tile%q1 = q - ql, vpt1 = pt1 * (1+0.61q1)
        slurb_tile%pt1(i,j)  = thl0(i, j, k_atm) + (rlv/(cp * exnf(k_atm)))  * ql0(i,j,k_atm)
        slurb_tile%q1(i,j)   = qt0(i, j, k_atm) - ql0(i, j, k_atm)
        slurb_tile%vpt1(i,j) = slurb_tile%pt1(i,j) * ( 1.0_field_r + 0.61_field_r * slurb_tile%q1(i,j) )

        

        du = 0.5*(u0(i,j,1) + u0(i+1,j,1)) + cu
        dv = 0.5*(v0(i,j,1) + v0(i,j+1,1)) + cv
        slurb_tile%uv_abs1(i,j) = sqrt(du**2 + dv**2)
        ! slurb_tile%uv_abs1(i,j) = max(0.1, sqrt(du**2 + dv**2)) DALES VERSION
        slurb_tile%uv_eff1(i,j) = slurb_tile%uv_abs1(i,j)


        !--       If spinup is enabled for current run, use diurnal mean spinup pt as the initial
        !--       atmospheric boundary condition. Otherwise, use the first atmospheric grid level.
        !--       Does DALES have a spinup option?
        ! IF ( spinup )  THEN
        !     bc_atm = spinup_pt_mean * exnf(k_topo)
        ! ELSE
        bc_atm = slurb_tile%pt1(i,j) * exnf(k_topo)
        ! ENDIF

        slurb_tile%us_urb(i,j)  = 1.0_field_r
        slurb_tile%us_can(i,j)  = 1.0_field_r
        slurb_tile%us_roof(i,j) = 1.0_field_r
        slurb_tile%us_road(i,j) = 1.0_field_r

        slurb_tile%uv_abs_can(i,j) = slurb_tile%uv_abs_can_coef(i,j) * slurb_tile%uv_abs1(i,j)
        slurb_tile%uv_eff_can(i,j) = slurb_tile%uv_abs_can(i,j)

        slurb_tile%shf_urb(i,j)  = 0.0_field_r
        slurb_tile%qsws_urb(i,j) = 0.0_field_r

        slurb_tile%t_can_0(i,j)   = bc_atm
        slurb_tile%t_can_m(i,j) = slurb_tile%t_can_0(i,j)

        slurb_tile%t_indoor(i,j) = building_indoor_temperature
        slurb_tile%t_soil(i,j) = deep_soil_temperature
        slurb_tile%shf_external(i,j) = shf_external
        slurb_tile%qsws_external(i,j) = qsws_external

    !
    !--       For subsurface temps, a steady-state 1D heat equation solution will be used as the
    !--       initial temperature profile. This might or might not speed up the spinup process.
    !--       In case of windowless facade, set window temps to fill value to prevent meaningless
    !--       output values. Vice versa for the opposite case.
        IF ( slurb_tile%f_win(i,j) < 1.0_field_r )  THEN
            slurb_tile%t_wall_a_0(:,i,j) = calc_1d_heat_equation( SIZE( slurb_tile%t_wall_a_0, 1 ), bc_atm,         &
                                                        slurb_tile%t_indoor(i,j),                         &
                                                        slurb_tile%conductivity_wall(:,i,j) )
            slurb_tile%t_wall_a_m(:,i,j) = slurb_tile%t_wall_a_0(:,i,j)
            slurb_tile%t_wall_b_0(:,i,j)   = slurb_tile%t_wall_a_0(:,i,j)
            slurb_tile%t_wall_b_m(:,i,j) = slurb_tile%t_wall_a_0(:,i,j)
        ELSE
            ! PALM extra debug option. If data_output_raw=false, set wall temperatures to fill value to prevent meaningless output values.
            IF ( .NOT. data_output_raw )  THEN
                slurb_tile%t_wall_a_0(:,i,j)   = output_fill_value
                slurb_tile%t_wall_a_m(:,i,j) = output_fill_value
            ENDIF
        ENDIF

        IF ( slurb_tile%f_win(i,j) > 0.0_field_r )  THEN
            slurb_tile%t_win_a_0(:,i,j) = calc_1d_heat_equation( SIZE( slurb_tile%t_win_a_0, 1 ), bc_atm,           &
                                                        slurb_tile%t_indoor(i,j),                          &
                                                        slurb_tile%conductivity_win(:,i,j) )
            slurb_tile%t_win_a_m(:,i,j) = slurb_tile%t_win_a_0(:,i,j)
            slurb_tile%t_win_b_0(:,i,j)   = slurb_tile%t_win_a_0(:,i,j)
            slurb_tile%t_win_b_m(:,i,j) = slurb_tile%t_win_a_0(:,i,j)
        ELSE
            ! PALM extra debug option. If data_output_raw=false, set wall temperatures to fill value to prevent meaningless output values.
            IF ( .NOT. data_output_raw )  THEN
                slurb_tile%t_win_a_0(:,i,j)   = output_fill_value
                slurb_tile%t_win_a_m(:,i,j) = output_fill_value
            ENDIF
        ENDIF

        slurb_tile%t_roof_0(:,i,j) = calc_1d_heat_equation( SIZE( slurb_tile%t_roof_0, 1 ), bc_atm,                &
                                                    slurb_tile%t_indoor(i,j), slurb_tile%conductivity_roof(:,i,j) )
        slurb_tile%t_roof_m(:,i,j) = slurb_tile%t_roof_0(:,i,j)

        slurb_tile%t_road_0(:,i,j) = calc_1d_heat_equation( SIZE( slurb_tile%t_road_0, 1 ), bc_atm,                &
                                                    slurb_tile%t_soil(i,j), slurb_tile%conductivity_road(:,i,j) )
        slurb_tile%t_road_m(:,i,j) = slurb_tile%t_road_0(:,i,j)

        IF ( moist_physics )  THEN
            slurb_tile%vpt_can(i,j) = 0.0_field_r

            slurb_tile%q_can_0(i,j)      = slurb_tile%q1(i,j)
            slurb_tile%q_can_m(i,j)      = slurb_tile%q_can_0(i,j)
            slurb_tile%m_liq_roof_m(i,j)   = 0.0_field_r
            slurb_tile%m_liq_roof_0(i,j) = slurb_tile%m_liq_roof_m(i,j)
            slurb_tile%m_liq_road_m(i,j)   = 0.0_field_r
            slurb_tile%m_liq_road_0(i,j) = slurb_tile%m_liq_road_m(i,j)

            slurb_tile%q_roof(i,j) = slurb_tile%q1(i,j)
            slurb_tile%q_road(i,j) = slurb_tile%q_can_0(i,j)
        ENDIF

        slurb_tile%ol_roof(i,j) = slurb_tile%z_mo(i,j)     / zeta_min
        slurb_tile%ol_road(i,j) = slurb_tile%z_mo_can(i,j) / zeta_min
        slurb_tile%ol_can(i,j)  = slurb_tile%z_mo(i,j)     / zeta_min
        slurb_tile%ol_urb(i,j)  = slurb_tile%z_mo(i,j)     / zeta_min

        !
        !--    Init potential temperatures and virtual potential temperatures. These need to be computed
        !--    also for the restart case, as d_exner is not yet available when rrd routines are called.
        slurb_tile%pt_can(i,j) = slurb_tile%t_can_0(i,j) / exnf(k_topo)

        IF ( slurb_tile%f_win(i,j) < 1.0_field_r )  THEN
            slurb_tile%pt_wall_a(i,j) = slurb_tile%t_wall_a_0(nzt_wall,i,j) / exnf(k_topo)
            slurb_tile%pt_wall_b(i,j) = slurb_tile%t_wall_b_0(nzt_wall,i,j) / exnf(k_topo)
        ELSE
            ! PALM extra debug option. If data_output_raw=false, set wall temperatures to fill value to prevent meaningless output values.
            IF ( .NOT. data_output_raw )  THEN
                slurb_tile%pt_wall_a(i,j) = output_fill_value
                slurb_tile%pt_wall_b(i,j) = output_fill_value
            ENDIF
        ENDIF

        IF ( slurb_tile%f_win(i,j) > 0.0_field_r )  THEN
            slurb_tile%pt_win_a(i,j) = slurb_tile%t_win_a_0(nzt_win,i,j) / exnf(k_topo)
            slurb_tile%pt_win_b(i,j) = slurb_tile%t_win_b_0(nzt_win,i,j) / exnf(k_topo)
        ELSE
            ! PALM extra debug option. If data_output_raw=false, set wall temperatures to fill value to prevent meaningless output values.
            IF ( .NOT. data_output_raw )  THEN
                slurb_tile%pt_win_a(i,j) = output_fill_value
                slurb_tile%pt_win_b(i,j) = output_fill_value
            ENDIF
        ENDIF

        slurb_tile%pt_roof(i,j) = slurb_tile%t_roof_0(nzt_roof,i,j) / exnf(k_topo)
        slurb_tile%pt_road(i,j) = slurb_tile%t_road_0(nzt_road,i,j) / exnf(k_topo)

        IF ( moist_physics )  THEN
            slurb_tile%vpt_can(i,j)  = slurb_tile%pt_can(i,j)  * ( 1.0_field_r + 0.61_field_r * slurb_tile%q_can_0(i,j)  )
            slurb_tile%vpt_roof(i,j) = slurb_tile%pt_roof(i,j) * ( 1.0_field_r + 0.61_field_r * slurb_tile%q_roof(i,j) )
            slurb_tile%vpt_road(i,j) = slurb_tile%pt_road(i,j) * ( 1.0_field_r + 0.61_field_r * slurb_tile%q_road(i,j) )
        ENDIF

    !
    !--    Initialize tendencies to zero.
        slurb_tile%tt_can(i,j)      = 0.0_field_r
        slurb_tile%tt_wall_a(:,i,j) = 0.0_field_r
        slurb_tile%tt_wall_b(:,i,j) = 0.0_field_r
        slurb_tile%tt_win_a(:,i,j)  = 0.0_field_r
        slurb_tile%tt_win_b(:,i,j)  = 0.0_field_r
        slurb_tile%tt_roof(:,i,j)   = 0.0_field_r
        slurb_tile%tt_road(:,i,j)   = 0.0_field_r
        IF ( moist_physics )  THEN
            slurb_tile%tq_can(i,j)      = 0.0_field_r
            slurb_tile%tm_liq_roof(i,j) = 0.0_field_r
            slurb_tile%tm_liq_road(i,j) = 0.0_field_r
            slurb_tile%tm_roof_runoff(i,j) = 0.0_field_r
            slurb_tile%tm_road_runoff(i,j) = 0.0_field_r
            slurb_tile%tm_roof_precep(i,j) = 0.0_field_r
            slurb_tile%tm_road_precep(i,j) = 0.0_field_r
        ENDIF

    !
    !--    Initialize model variables which are not used prior to an assignment in the model itself.
    !--    Thus these initializations should not end up being used in code, but as this is not
    !--    guaranteed with e.g. future changes, initialize them nevertheless. For the same reason,
    !--    these are not included in the restart data. But if in the future there is an usage prior to
    !--    proper assignment by the model, the respective variable should be added to restart routines,
    !--    and given a proper intialization.
        slurb_tile%albedo_urb(i,j)     = 0.0_field_r
        slurb_tile%emiss_urb(i,j)      = 1.0_field_r
        slurb_tile%rad_lw_in_urb(i,j)  = 0.0_field_r
        slurb_tile%rad_lw_out_urb(i,j) = 0.0_field_r
        slurb_tile%rad_sw_in_urb(i,j)  = 0.0_field_r
        slurb_tile%rad_sw_out_urb(i,j) = 0.0_field_r
        slurb_tile%ram_urb(i,j)        = 1E3_field_r
        slurb_tile%rib_urb(i,j)        = 0.0_field_r
        slurb_tile%t_2m_urb(i,j)       = 0.0_field_r
        slurb_tile%t_c_urb(i,j)        = 0.0_field_r
        slurb_tile%t_h_urb(i,j)        = 0.0_field_r
        slurb_tile%t_rad_urb(i,j)      = 0.0_field_r
        slurb_tile%usws_urb(i,j)       = 0.0_field_r
        slurb_tile%vsws_urb(i,j)       = 0.0_field_r

        slurb_tile%shf_can(i,j)    = 0.0_field_r
        slurb_tile%shf_road(i,j)   = 0.0_field_r
        slurb_tile%shf_roof(i,j)   = 0.0_field_r
        slurb_tile%shf_wall_a(i,j) = 0.0_field_r
        slurb_tile%shf_wall_b(i,j) = 0.0_field_r
        slurb_tile%shf_win_a(i,j)  = 0.0_field_r
        slurb_tile%shf_win_b(i,j)  = 0.0_field_r

        slurb_tile%ghf_road(i,j)   = 0.0_field_r
        slurb_tile%ghf_roof(i,j)   = 0.0_field_r
        slurb_tile%ghf_wall_a(i,j) = 0.0_field_r
        slurb_tile%ghf_wall_b(i,j) = 0.0_field_r
        slurb_tile%ghf_win_a(i,j)  = 0.0_field_r
        slurb_tile%ghf_win_b(i,j)  = 0.0_field_r

        slurb_tile%rad_lw_net_can(i,j)    = 0.0_field_r
        slurb_tile%rad_lw_net_road(i,j)   = 0.0_field_r
        slurb_tile%rad_lw_net_roof(i,j)   = 0.0_field_r
        slurb_tile%rad_lw_net_urb(i,j)    = 0.0_field_r
        slurb_tile%rad_lw_net_wall_a(i,j) = 0.0_field_r
        slurb_tile%rad_lw_net_wall_b(i,j) = 0.0_field_r
        slurb_tile%rad_lw_net_win_a(i,j)  = 0.0_field_r
        slurb_tile%rad_lw_net_win_b(i,j)  = 0.0_field_r
        slurb_tile%rad_sw_in_road(i,j)    = 0.0_field_r
        slurb_tile%rad_sw_in_win_a(i,j)   = 0.0_field_r
        slurb_tile%rad_sw_in_win_b(i,j)   = 0.0_field_r
        slurb_tile%rad_sw_net_road(i,j)   = 0.0_field_r
        slurb_tile%rad_sw_net_roof(i,j)   = 0.0_field_r
        slurb_tile%rad_sw_net_urb(i,j)    = 0.0_field_r
        slurb_tile%rad_sw_net_wall_a(i,j) = 0.0_field_r
        slurb_tile%rad_sw_net_wall_b(i,j) = 0.0_field_r
        slurb_tile%rad_sw_net_win_a (i,j) = 0.0_field_r
        slurb_tile%rad_sw_net_win_b (i,j) = 0.0_field_r

        slurb_tile%rib_can(i,j) = 0.0_field_r
        slurb_tile%rib_road(i,j) = 0.0_field_r
        slurb_tile%rib_roof(i,j) = 0.0_field_r

        slurb_tile%rah_can(i,j)  = 1E3_field_r
        slurb_tile%rah_road(i,j) = 1E3_field_r
        slurb_tile%rah_roof    = 1E3_field_r

        IF ( facade_rah_doe )  THEN
            slurb_tile%rah_wall_a(i,j) = 1E3_field_r
            slurb_tile%rah_wall_b(i,j) = 1E3_field_r
            slurb_tile%rah_win_a(i,j)  = 1E3_field_r
            slurb_tile%rah_win_b(i,j)  = 1E3_field_r
        ELSE
            slurb_tile%rah_facade(i,j) = 1E3_field_r
        ENDIF

        IF ( moist_physics )  THEN
            slurb_tile%qsws_can(i,j)      = 0.0_field_r
            slurb_tile%qsws_liq_road(i,j) = 0.0_field_r
            slurb_tile%qsws_liq_roof(i,j) = 0.0_field_r
            slurb_tile%qsws_road(i,j)     = 0.0_field_r
            slurb_tile%qsws_roof(i,j)     = 0.0_field_r

            e_s = magnus( MIN( slurb_tile%t_road_0(nzt_road,i,j), 333.15_field_r ) )
            slurb_tile%qs_road(i,j) = ep * e_s / ( ps - e_s )
            e_s = magnus( MIN( slurb_tile%t_roof_0(nzt_roof,i,j), 333.15_field_r ) )
            slurb_tile%qs_roof(i,j) = ep * e_s / ( ps - e_s )

            slurb_tile%c_liq_road(i,j)    = MIN( 1.0_field_r, ( slurb_tile%m_liq_road_m(i,j) / m_liq_max_road )**0.67 )
            slurb_tile%c_liq_roof(i,j)    = MIN( 1.0_field_r, ( slurb_tile%m_liq_roof_m(i,j) / m_liq_max_roof )**0.67 )
        ENDIF


        enddo
    enddo

    !
    !-- Calculate logarithms of ratio z/z0.
    !>  TODO: Since the ratios do not change during the simulation, they can be stored once at the
    !>        and stored in surf_slurb, like it is done for the other surface types, too.
            

    do j=2,j1
      do i=2,i1
       ln_z_z0_roof(i,j)  = LOG( slurb_tile%z_mo(i,j) / slurb_tile%z0_roof(i,j)  )
       ln_z_z0h_roof(i,j) = LOG( slurb_tile%z_mo(i,j) / slurb_tile%z0h_roof(i,j) )
       ln_z_z0_urb(i,j)   = LOG( slurb_tile%z_mo(i,j) / slurb_tile%z0_urb(i,j)   )
      enddo
    enddo

END SUBROUTINE init_slurb_variables

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Computes latent variables which can be inferred from the inputs.
    !--------------------------------------------------------------------------------------------------!
 SUBROUTINE precompute_latent_variables

    use modfields, only: rhof
    use modglobal, only: rhow, rlv, i1, j1, boltz
    use modmpi, only: comm3d

    use modmpi, only: comm3d, mpierr,mpi_min, D_MPI_ALLREDUCE
    
    implicit none

    character(len=*), parameter :: routine = modname//'/precompute_latent_variables'

    REAL(field_r) ::  emiss_facade       !< aggregated facade emissivity
    REAL(field_r) ::  f_wall             !< wall fraction
    REAL(field_r) ::  f_win              !< window fraction
    REAL(field_r) ::  wake               !< wake parameter for U_can parametrization in SURFEX
    REAL(field_r) ::  win_nonrefl_1side  !< 1-side nonreflected radiation (for windows)
    REAL(field_r) ::  win_absorp         !< window absorption coefficient
    real(field_r) :: dt_slurb_individual

    integer i,j,k
    !
    !-- Precompute model constants.
    rho_lv = rlv * rhof(1)
    drho_l_lv = 1.0_field_r / (rhow * rlv)

    !
    !-- Precompute layer total conductivities from layer thicknesses and thermal conductivities.
    do j=2,j1
        do i=2,i1
            DO  k = nzt_roof, nzb_roof-1
                slurb_tile%conductivity_roof(k,i,j) = 2.0_field_r / ( slurb_tile%dz_roof(k,i,j)   / slurb_tile%lambda_roof(k,i,j) +   &
                                                        slurb_tile%dz_roof(k+1,i,j) / slurb_tile%lambda_roof(k+1,i,j) )
            ENDDO
            slurb_tile%conductivity_roof(nzb_roof,i,j) = 2.0_field_r * slurb_tile%lambda_roof(nzb_roof,i,j) /                &
                                                    slurb_tile%dz_roof(nzb_roof,i,j)

            DO  k = nzt_wall, nzb_wall-1
                slurb_tile%conductivity_wall(k,i,j) = 2.0_field_r / ( slurb_tile%dz_wall(k,i,j)   / slurb_tile%lambda_wall(k,i,j) +   &
                                                        slurb_tile%dz_wall(k+1,i,j) / slurb_tile%lambda_wall(k+1,i,j) )
            ENDDO
            slurb_tile%conductivity_wall(nzb_wall,i,j) = 2.0_field_r * slurb_tile%lambda_wall(nzb_wall,i,j) /                &
                                                    slurb_tile%dz_wall(nzb_wall,i,j)

            DO  k = nzt_win, nzb_win-1
                slurb_tile%conductivity_win(k,i,j) = 2.0_field_r / ( slurb_tile%dz_win(k,i,j)   / slurb_tile%lambda_win(k,i,j) +      &
                                                        slurb_tile%dz_win(k+1,i,j) / slurb_tile%lambda_win(k+1,i,j) )
            ENDDO
            slurb_tile%conductivity_win(nzb_wall,i,j) = 2.0_field_r * slurb_tile%lambda_wall(nzb_win,i,j) /                  &
                                                slurb_tile%dz_wall(nzb_win,i,j)

        !
        !--    For the road, the last conductance depends on the soil conductance, so we need to
        !--    compute it during time-stepping (as it depends on soil moisture).
            DO  k = nzt_road, nzb_road-1
                slurb_tile%conductivity_road(k,i,j) = 2.0_field_r / ( slurb_tile%dz_road(k,i,j)   / slurb_tile%lambda_road(k,i,j) +   &
                                                        slurb_tile%dz_road(k+1,i,j) / slurb_tile%lambda_road(k+1,i,j) )
            ENDDO
            slurb_tile%conductivity_road(nzb_road,i,j) = 2.0_field_r * slurb_tile%lambda_road(nzb_road,i,j) /                &
                                                    slurb_tile%dz_road(nzb_road,i,j)
        enddo
    ENDDO

    !
    !-- Precompute sky-view factors.
    slurb_tile%svf_road(:,:) = 0.0_field_r
    slurb_tile%svf_wall(:,:) = 0.0_field_r

    do j=2,j1
        do i=2,i1
            slurb_tile%svf_road(i,j) = SQRT( (slurb_tile%hw_can(i,j))**2 + 1.0_field_r ) - slurb_tile%hw_can(i,j)

            slurb_tile%svf_wall(i,j) = ( slurb_tile%hw_can(i,j) + 1.0_field_r - SQRT( slurb_tile%hw_can(i,j)**2 + 1.0_field_r ) ) /       &
                                ( 2.0_field_r * slurb_tile%hw_can(i,j) )
        enddo
    ENDDO

    !
    !-- Precompute urban emissivity based on SVFs.
    do j=2,j1
        do i=2,i1
            slurb_tile%emiss_urb(i,j) = slurb_tile%f_bld(i,j) * slurb_tile%emiss_roof(i,j)                                      &
                                + ( 1.0_field_r - slurb_tile%f_bld(i,j) )                                            &
                                    * ( slurb_tile%svf_road(i,j) * slurb_tile%emiss_road(i,j)                            &
                                        + slurb_tile%svf_wall(i,j) * slurb_tile%hw_can(i,j)                                &
                                        * ( ( 1.0_field_r - slurb_tile%f_win(i,j) ) * 2.0_field_r * slurb_tile%emiss_wall(i,j)     &
                                            + slurb_tile%f_win(i,j) * 2.0_field_r * slurb_tile%emiss_win(i,j) ) )
        enddo
    ENDDO

    !
    !-- Preompute the longwave interaction coefficients for surface elements as these are
    !-- static in time. Based on Johnson et al. (1991) general formula. Absorption from reflected
    !-- radiation is taken into account only after first reflection. The first reflections contribute
    !-- around 5% of the total LW budget, while higher order reflections would contribute only <0.5%.
    !-- Coefficients are grouped per variable, so they can be effectively used in time-stepping
    !-- without wasting too much computational time or memory.
    do j=2,j1
        do i=2,i1
        !
        !--    Compute aggregated facade emissivity for simplification of reflections.
            f_win  = slurb_tile%f_win(i,j)
            f_wall = ( 1.0_field_r - f_win )
            emiss_facade = f_wall * slurb_tile%emiss_wall(i,j) + f_win * slurb_tile%emiss_win(i,j)
        !
        !--    Roof.
        !--    To be multiplied by t_roof**4 in the LW budget:
            slurb_tile%lw_roof_coef(1,i,j) = -slurb_tile%emiss_roof(i,j) * boltz
        !
        !--    To be multiplied by lw_rad_in_urb in the LW budget:
            slurb_tile%lw_roof_coef(2,i,j) = slurb_tile%emiss_roof(i,j)
        !
        !--    Roads.
        !--    To be multiplied by t_road**4 in the LW budget:
            slurb_tile%lw_road_coef(1,i,j) = ( - slurb_tile%emiss_road(i,j)                                             &
                                        + slurb_tile%emiss_road(i,j)**2 * ( 1 - emiss_facade )                   &
                                            * ( 1.0_field_r - slurb_tile%svf_road(i,j) ) * slurb_tile%svf_wall(i,j)             &
                                        ) * boltz
        !
        !--    To be multiplied by lw_rad_in in the LW budget:
            slurb_tile%lw_road_coef(2,i,j) = slurb_tile%emiss_road(i,j) * slurb_tile%svf_road(i,j)                              &
                                    - slurb_tile%emiss_road(i,j) * ( 1.0_field_r - emiss_facade ) * slurb_tile%svf_wall(i,j) &
                                        * ( 1.0_field_r - slurb_tile%svf_road(i,j) )
        !
        !--    To be multiplied by (t_wall_a**4 + t_wall_b**4) in the LW budget:
            slurb_tile%lw_road_coef(3,i,j) = ( slurb_tile%emiss_road(i,j) * slurb_tile%emiss_wall(i,j)                          &
                                            * ( 1.0_field_r - slurb_tile%svf_road(i,j) )                                &
                                        + slurb_tile%emiss_road(i,j) * slurb_tile%emiss_wall(i,j)                        &
                                            * ( 1.0_field_r - emiss_facade )                                    &
                                            * ( 1.0_field_r - slurb_tile%svf_road(i,j) )                                &
                                            * ( 1.0_field_r - 2.0_field_r * slurb_tile%svf_wall(i,j) )                       &
                                        ) * 0.5_field_r * f_wall * boltz
            slurb_tile%lw_road_coef(4,i,j) = ( slurb_tile%emiss_road(i,j) * slurb_tile%emiss_win(i,j)                           &
                                            * ( 1.0_field_r - slurb_tile%svf_road(i,j) )                                &
                                        + slurb_tile%emiss_road(i,j) * slurb_tile%emiss_win(i,j)                         &
                                            * ( 1.0_field_r - emiss_facade )                                    &
                                            * ( 1.0_field_r - slurb_tile%svf_road(i,j) )                                &
                                            * ( 1.0_field_r - 2.0_field_r * slurb_tile%svf_wall(i,j) )                       &
                                        ) * 0.5_field_r * f_win * boltz
        !
        !--    Walls.
        !--    To be multiplied by t_wall_a**4 in the LW budget:
            slurb_tile%lw_wall_coef(1,i,j) = ( - slurb_tile%emiss_wall(i,j)                                             &
                                        + 0.5_field_r * f_wall * slurb_tile%emiss_wall(i,j)**2                         &
                                        * ( 1.0_field_r - slurb_tile%emiss_road(i,j) ) * slurb_tile%svf_wall(i,j)            &
                                        * ( 1.0_field_r - slurb_tile%svf_road(i,j) )                                 &
                                        + f_wall * slurb_tile%emiss_wall(i,j)**2                                  &
                                        * ( 1.0_field_r - emiss_facade )                                     &
                                        * ( 1.0_field_r - 2.0_field_r * slurb_tile%svf_wall(i,j) )**2                     &
                                        ) * boltz
        !
        !--    To be multiplied by lw_rad_in in the LW budget:
            slurb_tile%lw_wall_coef(2,i,j) = slurb_tile%emiss_wall(i,j) * slurb_tile%svf_wall(i,j)                              &
                                        + slurb_tile%emiss_wall(i,j) * ( 1.0_field_r - slurb_tile%emiss_road(i,j) )             &
                                        * slurb_tile%svf_wall(i,j) * slurb_tile%svf_road(i,j)                            &
                                        + slurb_tile%emiss_wall(i,j)                                               &
                                        * ( 1.0_field_r - emiss_facade )                                      &
                                        * slurb_tile%svf_wall(i,j) * slurb_tile%svf_road(i,j)                            &
                                        + slurb_tile%emiss_wall(i,j)                                               &
                                        * ( 1.0_field_r - emiss_facade )                                      &
                                        * slurb_tile%svf_wall(i,j)                                               &
                                        * ( 1.0_field_r - 2.0_field_r * slurb_tile%svf_wall(i,j) )
        !
        !--    To be multiplied by t_wall_b**4 in the LW budget:
            slurb_tile%lw_wall_coef(3,i,j) = ( 0.5_field_r * f_wall * slurb_tile%emiss_wall(i,j)**2                          &
                                        * ( 1.0_field_r - slurb_tile%emiss_road(i,j) )                                &
                                        * slurb_tile%svf_wall(i,j)                                               &
                                        * ( 1.0_field_r - slurb_tile%svf_road(i,j) )                                  &
                                            + f_wall * slurb_tile%emiss_wall(i,j)**2                               &
                                            * ( 1.0_field_r - 2.0_field_r * slurb_tile%svf_wall(i,j) )                     &
                                        ) * boltz
        !
        !--    To be multiplied by t_win_a**4 in the LW budget:
            slurb_tile%lw_wall_coef(4,i,j) = ( 0.5_field_r * f_win * slurb_tile%emiss_wall(i,j)                              &
                                            * slurb_tile%emiss_win(i,j)                                            &
                                            * ( 1.0_field_r - slurb_tile%emiss_road(i,j) )                              &
                                            * slurb_tile%svf_wall(i,j)                                             &
                                            * ( 1.0_field_r - slurb_tile%svf_road(i,j) )                                &
                                        + f_win * slurb_tile%emiss_wall(i,j)                                     &
                                            * slurb_tile%emiss_win(i,j)                                           &
                                            * ( 1.0_field_r - emiss_facade )                                   &
                                            * ( 1.0_field_r - 2.0_field_r * slurb_tile%svf_wall(i,j) )**2                   &
                                        ) * boltz
        !
        !--    To be multiplied by t_win_b**4 in the LW budget:
            slurb_tile%lw_wall_coef(5,i,j) = ( 0.5_field_r * f_win * slurb_tile%emiss_wall(i,j)                              &
                                            * slurb_tile%emiss_win(i,j)                                           &
                                            * ( 1.0_field_r - slurb_tile%emiss_road(i,j) )                             &
                                            * slurb_tile%svf_wall(i,j)                                            &
                                            * ( 1.0_field_r - slurb_tile%svf_road(i,j) )                               &
                                        + f_win * slurb_tile%emiss_wall(i,j)                                     &
                                            * slurb_tile%emiss_win(i,j)                                           &
                                            * ( 1.0_field_r - 2.0_field_r * slurb_tile%svf_wall(i,j) )                      &
                                        ) * boltz
        !
        !--    To be multiplied by t_road**4 in the LW budget:
            slurb_tile%lw_wall_coef(6,i,j) = ( slurb_tile%emiss_wall(i,j) * slurb_tile%emiss_road(i,j)                          &
                                            * slurb_tile%svf_wall(i,j)                                            &
                                        + slurb_tile%emiss_wall(i,j) * slurb_tile%emiss_road(i,j)                        &
                                            * ( 1.0_field_r - emiss_facade )                                   &
                                            * slurb_tile%svf_wall(i,j)                                            &
                                            * ( 1.0_field_r - 2.0_field_r * slurb_tile%svf_wall(i,j) )                      &
                                        ) * boltz
        !
        !--    Windows.
        !--    To be multiplied by t_wall_a**4 in the LW budget:
            slurb_tile%lw_win_coef(1,i,j) = ( - slurb_tile%emiss_win(i,j)                                               &
                                        + 0.5_field_r * f_win * slurb_tile%emiss_win(i,j)**2                           &
                                        * ( 1.0_field_r - slurb_tile%emiss_road(i,j) )                               &
                                        * slurb_tile%svf_wall(i,j)                                              &
                                        * ( 1.0_field_r - slurb_tile%svf_road(i,j) )                                 &
                                        + f_win * slurb_tile%emiss_win(i,j)**2                                    &
                                        * ( 1.0_field_r - emiss_facade )                                     &
                                        * ( 1.0_field_r - 2.0_field_r * slurb_tile%svf_wall(i,j) )**2                     &
                                    ) * boltz
        !
        !--    To be multiplied by lw_rad_in in the LW budget:
            slurb_tile%lw_win_coef(2,i,j) = slurb_tile%emiss_win(i,j) * slurb_tile%svf_wall(i,j)                                &
                                    + slurb_tile%emiss_win(i,j)                                                 &
                                        * ( 1.0_field_r - slurb_tile%emiss_road(i,j) )                                 &
                                        * slurb_tile%svf_wall(i,j) * slurb_tile%svf_road(i,j)                             &
                                    + slurb_tile%emiss_win(i,j)                                                 &
                                        * ( 1.0_field_r - emiss_facade )                                       &
                                        * slurb_tile%svf_wall(i,j) * slurb_tile%svf_road(i,j)                             &
                                    + slurb_tile%emiss_win(i,j)                                                 &
                                        * ( 1.0_field_r - emiss_facade )                                       &
                                        * slurb_tile%svf_wall(i,j)                                                &
                                        * ( 1.0_field_r - 2.0_field_r * slurb_tile%svf_wall(i,j) )
        !
        !--    To be multiplied by t_win_b**4 in the LW budget:
            slurb_tile%lw_win_coef(3,i,j) = ( 0.5_field_r * f_win * slurb_tile%emiss_win(i,j)**2                             &
                                        * ( 1.0_field_r - slurb_tile%emiss_road(i,j) )                               &
                                        * slurb_tile%svf_wall(i,j)                                              &
                                        * ( 1.0_field_r - slurb_tile%svf_road(i,j) )                                 &
                                        + f_win * slurb_tile%emiss_win(i,j)**2                                    &
                                        * ( 1.0_field_r - 2.0_field_r * slurb_tile%svf_wall(i,j) )                        &
                                    ) * boltz
        !
        !--    To be multiplied by t_wall_a**4 in the LW budget:
            slurb_tile%lw_win_coef(4,i,j) = ( 0.5_field_r * f_wall * slurb_tile%emiss_win(i,j)                               &
                                        * slurb_tile%emiss_wall(i,j)                                            &
                                        * ( 1.0_field_r - slurb_tile%emiss_road(i,j) )                               &
                                        * slurb_tile%svf_wall(i,j)                                              &
                                        * ( 1.0_field_r - slurb_tile%svf_road(i,j) )                                 &
                                        + f_wall * slurb_tile%emiss_win(i,j)                                      &
                                        * slurb_tile%emiss_wall(i,j)                                            &
                                        * ( 1.0_field_r - emiss_facade )                                     &
                                        * ( 1.0_field_r - 2.0_field_r * slurb_tile%svf_wall(i,j) )**2                     &
                                    ) * boltz
        !
        !--    To be multiplied by t_wall_b**4 in the LW budget:
            slurb_tile%lw_win_coef(5,i,j) = ( 0.5_field_r * f_wall * slurb_tile%emiss_win(i,j)                               &
                                        * slurb_tile%emiss_wall(i,j)                                            &
                                        * ( 1.0_field_r - slurb_tile%emiss_road(i,j) )                               &
                                        * slurb_tile%svf_wall(i,j)                                              &
                                        * ( 1.0_field_r - slurb_tile%svf_road(i,j) )                                 &
                                        + f_wall * slurb_tile%emiss_win(i,j)                                      &
                                        * slurb_tile%emiss_wall(i,j)                                            &
                                        * ( 1.0_field_r - 2.0_field_r * slurb_tile%svf_wall(i,j) )                        &
                                    ) * boltz
        !
        !--    To be multiplied by t_road**4 in the LW budget:
            slurb_tile%lw_win_coef(6,i,j) = ( slurb_tile%emiss_win(i,j) * slurb_tile%emiss_road(i,j)                            &
                                        * slurb_tile%svf_wall(i,j)                                              &
                                        + slurb_tile%emiss_win(i,j) * slurb_tile%emiss_road(i,j)                          &
                                        * ( 1.0_field_r - emiss_facade )                                     &
                                        * slurb_tile%svf_wall(i,j)                                              &
                                        * ( 1.0_field_r - 2.0_field_r * slurb_tile%svf_wall(i,j) )                        &
                                    ) * boltz
        enddo
    ENDDO

    !
    !-- Precompute shortwave radiation reflection denominator.
    do j=2,j1
        do i=2,i1
            slurb_tile%sw_ref_denom(i,j) = 1.0_field_r - slurb_tile%albedo_road(i,j)                                         &
                                    * slurb_tile%albedo_wall_win(i,j)                                             &
                                    * slurb_tile%svf_wall(i,j) * ( 1.0_field_r - slurb_tile%svf_road(i,j) )                    &
                                    - slurb_tile%albedo_wall_win(i,j) * ( 1.0_field_r - 2.0_field_r * slurb_tile%svf_wall(i,j) )
        enddo
    ENDDO

    !
    !-- Compute window layer shortwave absorption based on USM documentation.
    !-- @todo This computation needs checking. Now sw_transmitted is not simply equal to
    !-- sw_net_win*transmissivity, as 1.0_field_r - SUM(absorption(:,i,j)) != transmissivity(i,j). This is
    !-- mitigated for now at the output side. No side effects for the model, as the transmitted
    !-- radiation is purely an output.
    do j=2,j1
        do i=2,i1
            win_nonrefl_1side = 1.0 - (slurb_tile%albedo_win(i,j) + slurb_tile%transmissivity_win(i,j)                  &
                                + 1.0_field_r  - SQRT( ( slurb_tile%albedo_win(i,j)                                  &
                                + slurb_tile%transmissivity_win(i,j) + 1.0_field_r )**2                              &
                                - 4.0_field_r * slurb_tile%albedo_win(i,j) ) ) / 2.0_field_r

            win_absorp = -LOG( ( slurb_tile%transmissivity_win(i,j) + slurb_tile%albedo_win(i,j)                        &
                                    - 1.0_field_r + win_nonrefl_1side ) / win_nonrefl_1side                     &
                                ) / slurb_tile%zw_win(nzb_win,i,j)

            DO  k = nzt_win, nzb_win
                IF ( k /= nzt_win)  THEN
        !
        !--          The absorbed fraction is difference between cumulative absorption over the layer.
                    slurb_tile%absorption_win(k,i,j) = win_nonrefl_1side                                          &
                                                * ( EXP( -win_absorp * slurb_tile%zw_win(k-1,i,j) ) -              &
                                                    EXP( -win_absorp * slurb_tile%zw_win(k,i,j)   ) )
                ELSE
        !
        !--          For the first layer, it is the cumulative absorption so far.
                    slurb_tile%absorption_win(k,i,j) = win_nonrefl_1side *                                        &
                                                ( 1.0_field_r - EXP( -win_absorp * slurb_tile%zw_win(k,i,j) ) )
                ENDIF
            ENDDO
        enddo
    ENDDO

    !
    !-- Coefficient for the canyon wind speed Krayenhoff & Voogt (2007) Eq. (9).
    IF ( uv_can_factor_kray )  THEN
        do j=2,j1
            do i=2,i1
                slurb_tile%uv_abs_can_coef(i,j) = LOG( slurb_tile%h_bld(i,j)  / ( 3.0_field_r * slurb_tile%z0_urb(i,j) ) ) /          &
                                    LOG( ( slurb_tile%z_mo(i,j) + slurb_tile%h_bld(i,j) / 3.0_field_r ) / slurb_tile%z0_urb(i,j) ) * &
                                    EXP( -slurb_tile%f_bld_frn(i,j) / ( 2.0_field_r * ( 1.0_field_r - slurb_tile%f_bld(i,j) ) ) )
            enddo
       ENDDO
    !
    !-- Coefficient for the canyon windspeed as derived in Masson (2000) (original TEB)
    ELSEIF ( uv_can_factor_masson )  THEN
        do j=2,j1
            do i=2,i1
                slurb_tile%uv_abs_can_coef(i,j) = LOG( slurb_tile%h_bld(i,j)  / ( 3.0_field_r * slurb_tile%z0_urb(i,j) ) ) /          &
                                    LOG( ( slurb_tile%z_mo(i,j) + slurb_tile%h_bld(i,j) / 3.0_field_r ) / slurb_tile%z0_urb(i,j) ) * &
                                    EXP( -slurb_tile%hw_can(i,j) / 4.0_field_r )
            enddo
       ENDDO
    !
    !-- Coefficient for the canyon windspeed as implemented in SURFEX v8.1
    ELSEIF ( uv_can_factor_surfex )  THEN
        do j=2,j1
            do i=2,i1
                wake = 1.0_field_r + ( 2.0_field_r / pi - 1.0_field_r ) * 2.0 * ( slurb_tile%hw_can(i,j) - 0.5_field_r )
                wake = MAX( MIN( wake, 1.0_field_r ), 2.0_field_r / pi )
                slurb_tile%uv_abs_can_coef(i,j) = wake * EXP( - slurb_tile%hw_can(i,j) / 4.0_field_r ) *                      &
                                            LOG( 2.0_field_r * slurb_tile%h_bld(i,j) / ( 3.0_field_r * slurb_tile%z0_urb(i,j) ) ) /  &
                                            LOG( ( slurb_tile%z_mo(i,j) + 2.0_field_r * slurb_tile%h_bld(i,j) ) /               &
                                                ( 3.0_field_r * slurb_tile%z0_urb(i,j) ) )
        enddo
       ENDDO
    ELSE
        call finish(routine, 'no canyon coefficient set for calculating absolute canyon velocity, this might have unintended consequences.')
    endif
    !
    !-- Compute minimum timestep based on SLUrb internal diffusivities.
    do j=2,j1
        do i=2,i1

        !
        !--    Criterion based on subsurface heat diffusion. Heat capacities are already multiplied with
        !--    layer thickness and conductivity is already divided with it. Thus, no need to multiply with
        !--    dz**2 like in usm and lsm. Dimension analysis:
        !--    c [J m-2 K-1] and layer conductivity [W m-2 K-1] -> [J/W] -> [s]
            slurb_tile%dt_max(i,j) = MIN( slurb_tile%dt_max(i,j),                                                       &
                                    MINVAL( slurb_tile%c_roof(:,i,j) / slurb_tile%conductivity_roof(:,i,j) ) )

            slurb_tile%dt_max(i,j) = MIN( slurb_tile%dt_max(i,j),                                                       &
                                    MINVAL( slurb_tile%c_road(:,i,j) / slurb_tile%conductivity_road(:,i,j) ) )

            slurb_tile%dt_max(i,j) = MIN( slurb_tile%dt_max(i,j),                                                       &
                                    MINVAL( slurb_tile%c_wall(:,i,j) / slurb_tile%conductivity_wall(:,i,j) ) )

            IF ( slurb_tile%f_win(i,j) > 0.0_field_r )  THEN
                slurb_tile%dt_max(i,j) = MIN( slurb_tile%dt_max(i,j),                                                    &
                                        MINVAL( slurb_tile%c_win(:,i,j) / slurb_tile%conductivity_win(:,i,j) ) )
            ENDIF
        enddo

    ENDDO
    !
    !-- Consider a pre-factor (1/8) for the diffusion criterion.
    dt_slurb_individual = MINVAL( slurb_tile%dt_max(2:i1,2:j1) ) * 0.125_field_r
    ! IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL D_MPI_ALLREDUCE( dt_slurb, dt_slurb_individual, 1, mpi_min, comm3d, mpierr )

    call warning(routine, 'Maximum timestep for SLUrb estimated to be ', dt_slurb, ' seconds. This is probably calculated incorrectly, and this is NOT applied to tstep yet!!')
    ! #endif

 END SUBROUTINE precompute_latent_variables




 !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Compute the dynamical conditions (wind speed, pt, q, vpt) in the street canyon.
    !--------------------------------------------------------------------------------------------------!
SUBROUTINE slurb_canyon_model
    use modglobal, only : i1, j1, cp, rlv, rk3step, rdt
    INTEGER ::  i       !< loop index (x-direction)
    INTEGER ::  j       !< loop index (y-direction)
    INTEGER ::  k_topo  !< k index of topography
    INTEGER ::  k_atm   !< k index of the first atmospheric level

    LOGICAL  ::  runge_l  !< timestep scheme switch for vectorization

    REAL(field_r) ::  c           !< total heat capacity of canyon air column per square metre (J K^-1 m^-2)
    REAL(field_r) ::  coef_1      !< coefficient A for the prognostic equation (W m^-2)
    REAL(field_r) ::  coef_2      !< coefficient B for the prognostic equation (W m^-2 K^-1)
    REAL(field_r) ::  f_shf       !< factor for the sensible heat flux  (W m^-2 K^-1)
    REAL(field_r) ::  f_qsws      !< factor for the latent heat flux (W m^-2)
    REAL(field_r) ::  qsws_surf   !< aggregated latent heat flux from canyon surfaces per unit area (W m^-2)
    REAL(field_r) ::  shf_surf    !< aggregated sensible heat flux from canyon surfaces per unit area (W m^-2)
    REAL(field_r) ::  tq_new      !< mixing ratio tendency for the new RK3 time step
    REAL(field_r) ::  tt_new      !< temperature tendency for the new RK3 time step
    REAL(field_r) ::  vtws        !< buoyancy flux (m K s^-1)
    REAL(field_r) ::  ws          !< free-convection scale (m/s)
    REAL(field_r) ::  q_can_p_imp !< mixing ratio for the new RK3 time step
    REAL(field_r) ::  t_can_p_imp !< temperature for the new RK3 time step

    real :: rk3coef !< (s)
    real :: rho_cp !< cp * rho (J m^-3 K^-1)

    rk3coef = rdt / (4. - dble(rk3step))

    k_topo = 1
    k_atm = 1
    rho_cp = cp * rho_air_zw(k_topo)

    ! runge_l = ( timestep_scheme(1:5) == 'runge' )
    runge_l = .true.

        do j=2,j1
      do i=2,i1

   !
   !--    Index offset of surface element point with respect to adjoining atmospheric grid point.
         !  k_topo = topo_top_ind(j,i,0)
         !  k_atm  = topo_top_ind(j,i,0) + 1


         f_shf  = rho_cp / slurb_tile%rah_can(i,j)
   !
   !--    Consider total air mass column within the street canyon.
         c = rho_cp * slurb_tile%h_bld(i,j)
   !
   !--    In canyon temperature prognostic equation, we use already computed fluxes from surfaces
   !--    in order to ensure consistency and conservation of energy. Thus, only the fluxes between
   !--    canyon air and the atmosphere are linearized.
   !
   !--    Aggregated sensible heat flux from canyon surfaces (per unit area).
         shf_surf = slurb_tile%hw_can(i,j) * ( ( 1.0_field_r - slurb_tile%f_win(i,j) ) *                                  &
                                       ( slurb_tile%shf_wall_a(i,j) + slurb_tile%shf_wall_b(i,j) ) +                 &
                                       slurb_tile%f_win(i,j) * ( slurb_tile%shf_win_a(i,j) + slurb_tile%shf_win_b(i,j) )     &
                                    ) + slurb_tile%shf_road(i,j)
   !
   !--    Aggregated flux doesn't contain c_p yet.
         shf_surf = shf_surf

   !
   !--    Coefficients for the prognostic equation of street canyon temperature.
         coef_1 = f_shf * slurb_tile%pt1(i,j) + (shf_surf)
         coef_2 = f_shf * (1 / exnf(k_topo))

        t_can_p_imp = (coef_1 * rk3coef + c * slurb_tile%t_can_0(i,j)) / (c + coef_2 * rk3coef)
        slurb_tile%tt_can(i,j) = (t_can_p_imp - slurb_tile%t_can_0(i,j)) / rk3coef
        slurb_tile%t_can_0(i,j) = t_can_p_imp

         !
         !--    Calculate new pt and shf from canyon to atmosphere.
         slurb_tile%pt_can(i,j) = slurb_tile%t_can_0(i,j) * (1 / exnf(k_topo))
         slurb_tile%shf_can(i,j) = -f_shf * ( slurb_tile%pt1(i,j) - slurb_tile%pt_can(i,j) )

   !
   !--    Compute prognostic street canyon mixing ratio.
         IF ( moist_physics )  THEN

            f_qsws = rho_lv / slurb_tile%rah_can(i,j)
   !
   !--       Same for the latent heat flux. Currently only the roads, walls are always dry.
   !--       This is a placeholder aggregation for street canyon vegetation,
   !--       e.g. green walls, low vegetation etc.
   !
   !--       Compute new prognostic canyon mixing ratio.
            coef_1 = f_qsws * slurb_tile%q1(i,j) + slurb_tile%qsws_road(i,j)
            coef_2 = f_qsws


            !--       Prevent negative mixing ratios due to temporal discretization. This is done before
            !--       the computation of tq_new in order to conserve energy.
            !--       Here our "latent heat capacity" is the canyon air column total mass.
            c = rho_lv * slurb_tile%h_bld(i,j)

            q_can_p_imp = (( coef_1 * rk3coef + c * slurb_tile%q_can_0(i,j) ) / ( c + coef_2 * rk3coef ))
            IF ( q_can_p_imp < 0.0_field_r )  q_can_p_imp = 0.0_field_r
            slurb_tile%tq_can(i,j) = (q_can_p_imp - slurb_tile%q_can_0(i,j)) / rk3coef
            slurb_tile%q_can_0(i,j) = q_can_p_imp

            slurb_tile%vpt_can(i,j) = slurb_tile%pt_can(i,j) * ( 1.0_field_r + 0.61_field_r * slurb_tile%q_can_0(i,j) )
            slurb_tile%qsws_can(i,j) = - f_qsws * ( slurb_tile%q1(i,j) - slurb_tile%q_can_0(i,j) )

         ENDIF
   !
   !--    Compute the canyon horizontal wind speed, Eq. (9), Krayenhoff & Voogt (2007).
         slurb_tile%uv_abs_can(i,j) = slurb_tile%uv_abs_can_coef(i,j) * slurb_tile%uv_abs1(i,j)
   !
   !--    Calculate the canyon effective wind speed taking into account turbulent processes
   !--    Lemonsu et al. (2004) Eqs. (2-3).
   !
   !--    Free convection scale (wstar) at building roof height (for unstable cases)
   !--    In case of moist physics, use virtual temperature (buoyancy) flux in free-convection scale.
        !  write(*,*) "shf_can"
        !  write (*,*) i,j,slurb_tile%shf_can(i,j)
        !  write(*,*) "qsws_can"
        !  write (*,*) i,j,slurb_tile%qsws_can(i,j)
         IF ( moist_physics )  THEN
            vtws =  (1/(rho_cp)) * slurb_tile%shf_can(i,j) + (1/(rho_cp)) * slurb_tile%qsws_can(i,j)
         ELSE
            vtws =  (1/(rho_cp)) * slurb_tile%shf_can(i,j)
         ENDIF
   !
   !--    No scaling for stable cases:
        !  write(*,*) "VTWS"
        !  write (*,*) i,j,vtws
         vtws = MERGE( vtws, 0.0_field_r, vtws > 0.0_field_r )
         ws = ( g / slurb_tile%pt_can(i,j) * slurb_tile%z_mo_can(i,j) * vtws )**( 1.0_field_r / 3.0_field_r )
   !
   !--    Canyon effective wind speed taking account both mean and turbulent wind.
         slurb_tile%uv_eff_can(i,j) = SQRT( slurb_tile%uv_abs_can(i,j)**2 + ( slurb_tile%us_can(i,j) + ws )**2 )

      enddo
    enddo

    ! IF ( debug_output_timestep )  THEN
    !    WRITE( debug_string, * ) 'slurb_canyon_model'
    !    CALL debug_message( debug_string, 'end' )
    ! ENDIF

 END SUBROUTINE slurb_canyon_model


    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> SLUrb's internal model to model urban surface - atmosphere coupling.
    !--------------------------------------------------------------------------------------------------!
 SUBROUTINE slurb_urban_aggregation_model
    use modglobal, only : i1, j1, rlv, cp
    use modfields, only : u0, v0, rhof, qt0, thl0
    INTEGER ::  i       !< running index
    INTEGER ::  j       !< running index
    INTEGER ::  k_topo  !< k index of topography
    INTEGER ::  k_atm   !< k index of the first atmospheric level

    LOGICAL ::  runge_l  !< flag for timestep scheme to allow vectorization
    real :: rhocp_i, rholv_i


    do j=2,j1
      do i=2,i1
      !  k_topo = topo_top_ind(j,i,0)
      !  k_atm = topo_top_ind(j,i,0) + 1
       k_topo = 1
       k_atm = 1

        !
        !--    For shf and qsws, use direct aggregation.
       ! must be W m^-2
       slurb_tile%shf_urb(i,j) = slurb_tile%f_bld(i,j) * slurb_tile%shf_roof(i,j) +                                        &
                         ( 1.0_field_r - slurb_tile%f_bld(i,j) ) * slurb_tile%shf_can(i,j) + slurb_tile%shf_external(i,j)

       IF ( moist_physics )  THEN
          slurb_tile%qsws_urb(i,j) = slurb_tile%f_bld(i,j) * slurb_tile%qsws_roof(i,j) +                                   &
                             ( 1.0_field_r - slurb_tile%f_bld(i,j) ) * slurb_tile%qsws_can(i,j) + slurb_tile%qsws_external(i,j)
       ENDIF

    !
    !--    Calculate momentum flux for horizontal wind components.
       slurb_tile%usws_urb(i,j) = -u0(i,j,k_atm) / slurb_tile%ram_urb(i,j) * rho_air_zw(k_topo)
       slurb_tile%vsws_urb(i,j) = -v0(i,j,k_atm) / slurb_tile%ram_urb(i,j) * rho_air_zw(k_topo)

    !
    !--    Aggregate radiative fluxes. Note that this aggregation is done here rather than in the
    !--    slurb_radiation_model on purpose to include the longwave term dependent on
    !--    the surface temperature of given surface.
    !
    !--    Compute the net LW radiation flux at canyon top and for urban surface.
       slurb_tile%rad_lw_net_can(i,j) = slurb_tile%rad_lw_net_road(i,j) + slurb_tile%hw_can(i,j) *                         &
                                (   ( 1.0_field_r - slurb_tile%f_win(i,j) ) *                                   &
                                    ( slurb_tile%rad_lw_net_wall_a(i,j) + slurb_tile%rad_lw_net_wall_b(i,j) )      &
                                  + slurb_tile%f_win(i,j) *                                                &
                                    ( slurb_tile%rad_lw_net_win_a(i,j)  + slurb_tile%rad_lw_net_win_b(i,j) )       &
                                )

       slurb_tile%rad_lw_net_urb(i,j) = slurb_tile%f_bld(i,j) * slurb_tile%rad_lw_net_roof(i,j) +                          &
                                ( 1.0_field_r - slurb_tile%f_bld(i,j) ) * slurb_tile%rad_lw_net_can(i,j)
    !
    !--    Outgoing LW flux.
       slurb_tile%rad_lw_out_urb(i,j) = slurb_tile%rad_lw_in_urb(i,j) - slurb_tile%rad_lw_net_urb(i,j)
    !
    !--    Calculate urban aggregated surface temperatures.
       CALL calc_urban_aggregated_temperatures

        rhocp_i = 1. / (rhof(1) * cp)
        rholv_i = 1. / (rhof(1) * rlv)
        ! Calculate surface values
        slurb_tile%thlskin(i,j) = thl0(i,j,1) + (slurb_tile%shf_urb(i,j) * rhocp_i) * slurb_tile%ram_urb(i,j)
        slurb_tile%qtskin (i,j) = qt0(i,j,1) + (slurb_tile%qsws_urb(i,j) * rholv_i) * slurb_tile%ram_urb(i,j)

      enddo
    enddo

 CONTAINS


    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Update the aggregated urban surface temperatures. These are diagnostic outputs and are not
    !> prognostic model variables. Four aggregated urban surface temperatures are computed:
    !> 1) Effective surface temperature T_H derived from conservation of heat flux contributions
    !> 2) Radiative surface temperature T_rad derived from the outgoing LW radiation
    !> 3) Complete surface temperature T_C which is an area-weighted temperature of all facets
    !> 4) Theoretical temperature at 2 m height extrapolated using stability-corrected log profile
    !> For 1-3 formulations of Kanda et al. 2005, adapted for SLUrb configuration, are used.
    !< Note that prognostic temperatures (suffix _p) are used. These are the ones that are output
    !> for current time step, as the timelevel is swapped right after the prognostic equation calls.
    !--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_urban_aggregated_temperatures
    use modglobal, only : boltz, rlv, cp
    use modslurb_resistance_stability, only: psi_h, psi_m
    REAL(field_r) ::  c_h_roof    !< bulk heat transfer coefficient for roof (J kg^-1 K^-1)
    REAL(field_r) ::  c_h_wall_a  !< bulk heat transfer coefficient for wall a (J kg^-1 K^-1)
    REAL(field_r) ::  c_h_wall_b  !< bulk heat transfer coefficient for wall b (J kg^-1 K^-1)
    REAL(field_r) ::  c_h_win_a   !< bulk heat transfer coefficient for window a (J kg^-1 K^-1)
    REAL(field_r) ::  c_h_win_b   !< bulk heat transfer coefficient for window b (J kg^-1 K^-1)
    REAL(field_r) ::  c_h_road    !< bulk heat transfer coefficient for road (J kg^-1 K^-1)
    REAL(field_r) ::  ts          !< scaling temperature (K)
    REAL(field_r) ::  vtws        !< virtual potential temperature flux (buoyancy flux) (m K s^-1)
    real :: rho_cp !< (J m^-3 K^-1)
    real :: rho !< (kg m^-3)
    rho_cp = cp * rho_air_zw(k_topo)
    rho = rho_air_zw(k_topo)
    !
    !-- 1) Effective surface temperature T_H.
    !-- First, compute the bulk heat transfer coefficients.
    IF ( calc_t_h )  THEN
       c_h_roof = ABS( slurb_tile%shf_roof(i,j) / ( rho * slurb_tile%uv_eff1(i,j) *                             &
                       ( slurb_tile%t_roof_0(nzt_roof,i,j) - slurb_tile%pt1(i,j) * exnf(k_atm) ) ) )

       c_h_wall_a = ABS( slurb_tile%shf_wall_a(i,j) / ( rho * slurb_tile%uv_eff1(i,j) *                         &
                         ( slurb_tile%t_wall_a_0(nzt_wall,i,j) - slurb_tile%pt1(i,j) * exnf(k_atm) ) ) )

       c_h_wall_b = ABS( slurb_tile%shf_wall_b(i,j) / ( rho * slurb_tile%uv_eff1(i,j) *                         &
                         ( slurb_tile%t_wall_b_0(nzt_wall,i,j) - slurb_tile%pt1(i,j) * exnf(k_atm) ) ) )

       c_h_win_a = ABS( slurb_tile%shf_win_a(i,j) / ( rho * slurb_tile%uv_eff1(i,j) *                           &
                        ( slurb_tile%t_win_a_0(nzt_win,i,j) - slurb_tile%pt1(i,j) * exnf(k_atm) ) ) )

       c_h_win_b = ABS( slurb_tile%shf_win_b(i,j) / ( rho * slurb_tile%uv_eff1(i,j) *                           &
                        ( slurb_tile%t_win_b_0(nzt_win,i,j) - slurb_tile%pt1(i,j) * exnf(k_atm) ) ) )

       c_h_road = ABS( slurb_tile%shf_road(i,j) / ( rho * slurb_tile%uv_eff1(i,j) *                             &
                       ( slurb_tile%t_road_0(nzt_road,i,j) - slurb_tile%pt1(i,j) * exnf(k_atm) ) ) )

       slurb_tile%t_h_urb(i,j) = ( ( 1.0_field_r - slurb_tile%f_bld(i,j) ) *                                            &
                           ( slurb_tile%hw_can(i,j) * (                                                    &
                                                ( 1.0_field_r - slurb_tile%f_win(i,j) ) *                       &
                                                ( c_h_wall_a * slurb_tile%t_wall_a_0(nzt_wall,i,j)         &
                                                + c_h_wall_b * slurb_tile%t_wall_b_0(nzt_wall,i,j) )       &
                                              + slurb_tile%f_win(i,j) *                                    &
                                                ( c_h_win_a * slurb_tile%t_win_a_0(nzt_win,i,j)            &
                                                + c_h_win_b * slurb_tile%t_win_b_0(nzt_win,i,j) )          &
                                              )                                                    &
                           + c_h_road * slurb_tile%t_road_0(nzt_road,i,j)                                  &
                           )                                                                       &
                         + slurb_tile%f_bld(i,j) * c_h_roof * slurb_tile%t_roof_0(nzt_roof,i,j)                    &
                         ) /                                                                       &
                         ( ( 1.0_field_r - slurb_tile%f_bld(i,j) ) *                                            &
                           ( slurb_tile%hw_can(i,j) * (                                                    &
                                                ( 1.0_field_r - slurb_tile%f_win(i,j) ) *                       &
                                                ( c_h_wall_a + c_h_wall_b )                        &
                                              + slurb_tile%f_win(i,j) *                                    &
                                                ( c_h_win_a + c_h_win_b )                          &
                                              )                                                    &
                           + c_h_road                                                              &
                           )                                                                       &
                         + slurb_tile%f_bld(i,j) * c_h_roof + 1E-10_field_r                                     &
                         )
    ENDIF

    !
    !-- 2) Radiative surface temperature T_rad.
    slurb_tile%t_rad_urb(i,j) = SQRT( SQRT( slurb_tile%rad_lw_out_urb(i,j) / ( slurb_tile%emiss_urb(i,j) * boltz ) ) )

    !
    !-- 3) Complete surface temperature T_C, similarly to T_H but without the C_h weighting.
    IF ( calc_t_c )  THEN
       slurb_tile%t_c_urb(i,j) = ( ( 1.0_field_r - slurb_tile%f_bld(i,j) ) *                                            &
                           ( slurb_tile%hw_can(i,j) * ( ( 1.0_field_r - slurb_tile%f_win(i,j) ) *                       &
                                     ( slurb_tile%t_wall_a_0(nzt_wall,i,j) + slurb_tile%t_wall_b_0(nzt_wall,i,j) ) &
                                     + slurb_tile%f_win(i,j) *                                             &
                                     ( slurb_tile%t_win_a_0(nzt_win,i,j)   + slurb_tile%t_win_b_0(nzt_win,i,j)   ) &
                                              )                                                    &
                           + slurb_tile%t_road_0(nzt_road,i,j)                                             &
                           )                                                                       &
                           + slurb_tile%f_bld(i,j) * slurb_tile%t_roof_0(nzt_roof,i,j)                             &
                         ) /                                                                       &
                         ( ( 1.0_field_r - slurb_tile%f_bld(i,j) ) * ( 2.0_field_r * slurb_tile%hw_can(i,j) + 1.0_field_r )       &
                           + slurb_tile%f_bld(i,j)                                                         &
                         )
    ENDIF

    !
    !-- 4) Theoretical 2 m temperature extrapolated using MOST.
    IF ( calc_t_2m )  THEN
       IF ( moist_physics )  THEN
          vtws =  (1/(rho_cp)) * slurb_tile%shf_can(i,j) + (1/(rho_cp)) * slurb_tile%qsws_can(i,j)
       ELSE
          vtws =  (1/(rho_cp)) * slurb_tile%shf_can(i,j)
          ! m K s^-1 = (kg m^-3 J kg^-1 K^-1)^-1 W m^-2
          ! m K s^-1 = m^3 J^-1 K W m^-2
          ! m K s^-1 = m K s^-1
          ! for the moist case:
          ! m K s^-1 = J kg^-1 J^-1 kg K W m^-2
          ! m K s^-1 = K J s^-1 m^-2
          ! need to multiply by J^-1 m^3
       ENDIF
       ts = -vtws / slurb_tile%us_urb(i,j)

       slurb_tile%t_2m_urb(i,j) = ts / kappa *                                                             &
                          ( LOG( 2.0_field_r / ( slurb_tile%z_mo(i,j) + slurb_tile%h_bld(i,j) ) ) -                     &
                            psi_h( 2.0_field_r / slurb_tile%ol_urb(i,j) ) +                                     &
                            psi_h( ( slurb_tile%z_mo(i,j) + slurb_tile%h_bld(i,j) ) / slurb_tile%ol_urb(i,j) )             &
                          ) + slurb_tile%pt1(i,j) * exnf(k_atm)
    ENDIF


 END SUBROUTINE calc_urban_aggregated_temperatures

 END SUBROUTINE slurb_urban_aggregation_model

end module modslurb