module modslurbcrosssection

  use modlogging,         only : warning
  use modnetcdf_file_t,   only : cross_section_file_t, slurb_3d_file_t
  use modprecision,       only : field_r
  use modstat_nc_files,   only : add_output_file, is_sampling_timestep

  implicit none

  private

  public :: initslurbcrosssection, slurbcrosssection

  type(cross_section_file_t) :: slurb_urb_file
  type(slurb_3d_file_t)      :: slurb_roof_file
  type(slurb_3d_file_t)      :: slurb_road_file
  type(slurb_3d_file_t)      :: slurb_wall_file
  type(slurb_3d_file_t)      :: slurb_win_file

  integer :: slurb_urb_file_id = 0
  integer :: slurb_roof_file_id = 0
  integer :: slurb_road_file_id = 0
  integer :: slurb_wall_file_id = 0
  integer :: slurb_win_file_id = 0

  logical :: slurb_urb_enabled = .false.
  logical :: slurb_roof_enabled = .false.
  logical :: slurb_road_enabled = .false.
  logical :: slurb_wall_enabled = .false.
  logical :: slurb_win_enabled = .false.

contains

  subroutine write_slurb_constants()
    use modglobal,    only : i1, j1, itot, jtot
    use modslurbdata, only : slurb_tile, &
                             nzt_roof, nzb_roof, nzt_road, nzb_road, nzt_wall, nzb_wall, nzt_win, nzb_win

    implicit none

    type(cross_section_file_t) :: slurb_constants_2d_file
    type(slurb_3d_file_t)      :: slurb_constants_roof_file
    type(slurb_3d_file_t)      :: slurb_constants_road_file
    type(slurb_3d_file_t)      :: slurb_constants_wall_file
    type(slurb_3d_file_t)      :: slurb_constants_win_file

    real(field_r), pointer :: emiss_urb(:,:), z0_urb(:,:), theta_can(:,:), svf_wall(:,:), hw_can(:,:), f_bld(:,:), &
                  f_bld_frn(:,:), f_win(:,:), h_bld(:,:), &
                              t_indoor(:,:), t_soil(:,:), z_mo(:,:), z_mo_can(:,:), wall_hor_a_ratio(:,:), uv_abs_can_coef(:,:)
    real(field_r), pointer :: dz_roof(:,:,:), dz_road(:,:,:), dz_wall(:,:,:), dz_win(:,:,:)
    real(field_r), pointer :: c_roof(:,:,:), c_road(:,:,:), c_wall(:,:,:), c_win(:,:,:), absorption_win(:,:,:)

    slurb_constants_roof_file = slurb_3d_file_t('slurbcross_constants_roof', nzs=nzb_roof-nzt_roof+1, lgpu=.false.)
    call slurb_constants_roof_file%add_var('dz_roof', 'roof layer thickness', 'm', 'tttts_slurb')
    call slurb_constants_roof_file%add_var('c_roof', 'roof layer heat capacity', 'J/m^2/K', 'tttts_slurb')

    slurb_constants_road_file = slurb_3d_file_t('slurbcross_constants_road', nzs=nzb_road-nzt_road+1, lgpu=.false.)
    call slurb_constants_road_file%add_var('dz_road', 'road layer thickness', 'm', 'tttts_slurb')
    call slurb_constants_road_file%add_var('c_road', 'road layer heat capacity', 'J/m^2/K', 'tttts_slurb')

    slurb_constants_wall_file = slurb_3d_file_t('slurbcross_constants_wall', nzs=nzb_wall-nzt_wall+1, lgpu=.false.)
    call slurb_constants_wall_file%add_var('dz_wall', 'wall layer thickness', 'm', 'tttts_slurb')
    call slurb_constants_wall_file%add_var('c_wall', 'wall layer heat capacity', 'J/m^2/K', 'tttts_slurb')

    slurb_constants_win_file = slurb_3d_file_t('slurbcross_constants_win', nzs=nzb_win-nzt_win+1, lgpu=.false.)
    call slurb_constants_win_file%add_var('dz_win', 'window layer thickness', 'm', 'tttts_slurb')
    call slurb_constants_win_file%add_var('c_win', 'window layer heat capacity', 'J/m^2/K', 'tttts_slurb')
    call slurb_constants_win_file%add_var('absorption_win', 'absorption_win', '', 'tttts_slurb')

    slurb_constants_2d_file = cross_section_file_t('slurbcross_constants_2d', nx=itot, ny=jtot, lgpu=.false.)
    call slurb_constants_2d_file%add_var('z_mo', 'reference height for MOST for the atmosphere', 'm', 'tt0t')
    call slurb_constants_2d_file%add_var('z_mo_can', 'canyon reference height for MOST', 'm', 'tt0t')
    call slurb_constants_2d_file%add_var('wall_hor_a_ratio', 'wall-to-horizontal area ratio', '-', 'tt0t')
    call slurb_constants_2d_file%add_var('uv_abs_can_coef', 'coefficient for the canyon wind speed', '', 'tt0t')
    call slurb_constants_2d_file%add_var('z0_urb', 'aerodynamic roughness length of the urban surface', 'm', 'tt0t')
    call slurb_constants_2d_file%add_var('theta_can', 'canyon orientation / road direction in radians', 'rad', 'tt0t')
    call slurb_constants_2d_file%add_var('svf_wall', 'sky-view-factor for walls', '-', 'tt0t')
    call slurb_constants_2d_file%add_var('hw_can', 'canyon aspect ratio', '-', 'tt0t')
    call slurb_constants_2d_file%add_var('f_bld', 'fractional area occupied by buildings', '-', 'tt0t')
    call slurb_constants_2d_file%add_var('f_bld_frn', 'frontal area fraction of buildings', '-', 'tt0t')
    call slurb_constants_2d_file%add_var('f_win', 'window fraction', '-', 'tt0t')
    call slurb_constants_2d_file%add_var('h_bld', 'building height', 'm', 'tt0t')
    call slurb_constants_2d_file%add_var('emiss_urb', 'effective urban emissivity', '-', 'tt0t')
    call slurb_constants_2d_file%add_var('t_indoor', 'building indoor temperature', 'K', 'tt0t')
    call slurb_constants_2d_file%add_var('t_soil', 'fixed soil top temperature', 'K', 'tt0t')

    call slurb_constants_roof_file%open()
    call slurb_constants_road_file%open()
    call slurb_constants_wall_file%open()
    call slurb_constants_win_file%open()
    call slurb_constants_2d_file%open()

    call slurb_constants_roof_file%get_pointer('dz_roof', dz_roof)
    call slurb_constants_roof_file%get_pointer('c_roof', c_roof)
    call slurb_constants_road_file%get_pointer('dz_road', dz_road)
    call slurb_constants_road_file%get_pointer('c_road', c_road)
    call slurb_constants_wall_file%get_pointer('dz_wall', dz_wall)
    call slurb_constants_wall_file%get_pointer('c_wall', c_wall)
    call slurb_constants_win_file%get_pointer('dz_win', dz_win)
    call slurb_constants_win_file%get_pointer('c_win', c_win)
    call slurb_constants_win_file%get_pointer('absorption_win', absorption_win)

    call slurb_constants_2d_file%get_pointer('z_mo', z_mo)
    call slurb_constants_2d_file%get_pointer('z_mo_can', z_mo_can)
    call slurb_constants_2d_file%get_pointer('wall_hor_a_ratio', wall_hor_a_ratio)
    call slurb_constants_2d_file%get_pointer('uv_abs_can_coef', uv_abs_can_coef)
    call slurb_constants_2d_file%get_pointer('z0_urb', z0_urb)
    call slurb_constants_2d_file%get_pointer('theta_can', theta_can)
    call slurb_constants_2d_file%get_pointer('svf_wall', svf_wall)
    call slurb_constants_2d_file%get_pointer('hw_can', hw_can)
    call slurb_constants_2d_file%get_pointer('f_bld', f_bld)
    call slurb_constants_2d_file%get_pointer('f_bld_frn', f_bld_frn)
    call slurb_constants_2d_file%get_pointer('f_win', f_win)
    call slurb_constants_2d_file%get_pointer('h_bld', h_bld)
    call slurb_constants_2d_file%get_pointer('emiss_urb', emiss_urb)
    call slurb_constants_2d_file%get_pointer('t_indoor', t_indoor)
    call slurb_constants_2d_file%get_pointer('t_soil', t_soil)

    dz_roof(:,:,:) = slurb_tile%dz_roof(:,2:i1,2:j1)
    dz_road(:,:,:) = slurb_tile%dz_road(:,2:i1,2:j1)
    dz_wall(:,:,:) = slurb_tile%dz_wall(:,2:i1,2:j1)
    dz_win(:,:,:) = slurb_tile%dz_win(:,2:i1,2:j1)
    c_roof(:,:,:) = slurb_tile%c_roof(:,2:i1,2:j1)
    c_road(:,:,:) = slurb_tile%c_road(:,2:i1,2:j1)
    c_wall(:,:,:) = slurb_tile%c_wall(:,2:i1,2:j1)
    c_win(:,:,:) = slurb_tile%c_win(:,2:i1,2:j1)
    absorption_win(:,:,:) = slurb_tile%absorption_win(:,2:i1,2:j1)
    z_mo(:,:) = slurb_tile%z_mo(2:i1,2:j1)
    z_mo_can(:,:) = slurb_tile%z_mo_can(2:i1,2:j1)
    wall_hor_a_ratio(:,:) = 0.0_field_r
    uv_abs_can_coef(:,:) = slurb_tile%uv_abs_can_coef(2:i1,2:j1)
    z0_urb(:,:) = slurb_tile%z0_urb(2:i1,2:j1)
    theta_can(:,:) = slurb_tile%theta_can(2:i1,2:j1)
    svf_wall(:,:) = slurb_tile%svf_wall(2:i1,2:j1)
    hw_can(:,:) = slurb_tile%hw_can(2:i1,2:j1)
    f_bld(:,:) = slurb_tile%f_bld(2:i1,2:j1)
    f_bld_frn(:,:) = slurb_tile%f_bld_frn(2:i1,2:j1)
    f_win(:,:) = slurb_tile%f_win(2:i1,2:j1)
    h_bld(:,:) = slurb_tile%h_bld(2:i1,2:j1)
    emiss_urb(:,:) = slurb_tile%emiss_urb(2:i1,2:j1)
    t_indoor(:,:) = slurb_tile%t_indoor(2:i1,2:j1)
    t_soil(:,:) = slurb_tile%t_soil(2:i1,2:j1)

    call slurb_constants_roof_file%write()
    call slurb_constants_road_file%write()
    call slurb_constants_wall_file%write()
    call slurb_constants_win_file%write()
    call slurb_constants_2d_file%write()

    call slurb_constants_roof_file%close()
    call slurb_constants_road_file%close()
    call slurb_constants_wall_file%close()
    call slurb_constants_win_file%close()
    call slurb_constants_2d_file%close()
  end subroutine write_slurb_constants

  subroutine initslurbcrosssection
    use modglobal,    only : itot, jtot
    use modslurbdata, only : enable_slurb, dtav_slurb, output_slurb_bc, output_slurb_constants, slurb_cross_output, &
                             slurb_cross_output_roof, slurb_cross_output_road, slurb_cross_output_wall_win, &
                             slurb_cross_output_tendencies, slurb_cross_output_radiation, &
                             nzt_roof, nzb_roof, nzt_road, nzb_road, nzt_wall, nzb_wall, nzt_win, nzb_win

    implicit none

    if (.not. enable_slurb) return
    if (.not. slurb_cross_output) return


    if (output_slurb_constants) then
      call write_slurb_constants()
    end if

    slurb_urb_file = cross_section_file_t('slurbcross_urb', nx=itot, ny=jtot, lgpu=.false.)
    call add_output_file(slurb_urb_file, dtav_slurb, slurb_urb_file_id)
    slurb_urb_enabled = .true.

    call slurb_urb_file%add_var('albedo_urb', 'effective urban albedo', '-', 'tt0t')
    call slurb_urb_file%add_var('ol_urb', 'urban Obukhov length', 'L', 'tt0t')
    call slurb_urb_file%add_var('qsws_urb', 'total urban latent heat flux', 'W/m^2', 'tt0t')
    call slurb_urb_file%add_var('ram_urb', 'urban aerodynamic resistance for momentum', 's/m', 'tt0t')
    call slurb_urb_file%add_var('rib_urb', 'urban bulk-Richardson number', '-', 'tt0t')
    call slurb_urb_file%add_var('shf_urb', 'total urban sensible heat flux', 'W/m^2', 'tt0t')
    call slurb_urb_file%add_var('t_2m_urb', 'urban 2-metre temperature', 'K', 'tt0t')
    call slurb_urb_file%add_var('t_c_urb', 'complete urban surface temperature', 'K', 'tt0t')
    call slurb_urb_file%add_var('t_h_urb', 'effective urban surface temperature', 'K', 'tt0t')
    call slurb_urb_file%add_var('thl_rad_urb', 'urban radiative surface liquid water potential temperature', 'K', 'tt0t')
    call slurb_urb_file%add_var('usws_urb', 'urban momentum flux u-component', 'm^2/s^2', 'tt0t')
    call slurb_urb_file%add_var('vsws_urb', 'urban momentum flux v-component', 'm^2/s^2', 'tt0t')
    call slurb_urb_file%add_var('thlskin', 'urban skin liquid water potential temperature', 'K', 'tt0t')
    call slurb_urb_file%add_var('qtskin', 'urban skin specific humidity', 'kg/kg', 'tt0t')
    call slurb_urb_file%add_var('q_can_0', 'canyon mixing ratio', 'kg/kg', 'tt0t')
    call slurb_urb_file%add_var('t_can_0', 'canyon air temperature', 'K', 'tt0t')
    call slurb_urb_file%add_var('shf_can', 'sensible heat flux between the street canyon and the atmosphere', 'W/m^2', 'tt0t')
    call slurb_urb_file%add_var('shf_external', 'sensible heat flux external to the model', 'W/m^2', 'tt0t')
    call slurb_urb_file%add_var('shf_traffic', 'traffic sensible heat flux', 'W/m^2', 'tt0t')
    call slurb_urb_file%add_var('qsws_can', 'latent heat flux between the street canyon and the atmosphere', 'W/m^2', 'tt0t')
    call slurb_urb_file%add_var('qsws_external', 'latent heat flux external to the model', 'W/m^2', 'tt0t')
    call slurb_urb_file%add_var('ol_can', 'canyon top Obukhov length', 'm', 'tt0t')
    call slurb_urb_file%add_var('pt_can', 'street canyon virtual potential temperature', 'K', 'tt0t')
    call slurb_urb_file%add_var('rib_can', 'canyon top bulk Richardson number', '-', 'tt0t')
    call slurb_urb_file%add_var('us_can', 'friction velocity for canyon resistance calculation', 'm/s', 'tt0t')
    call slurb_urb_file%add_var('uv_abs_can', 'horizontal wind speed in street canyon at half-height', 'm/s', 'tt0t')
    call slurb_urb_file%add_var('uv_eff_can', 'effective horizontal wind speed in street canyon at half-height', 'm/s', 'tt0t')
    call slurb_urb_file%add_var('vpt_can', 'street canyon virtual potential temperature', 'K', 'tt0t')
    call slurb_urb_file%add_var('rah_can', 'street canyon air aerodynamic resistance for heat', 's/m', 'tt0t')
    call slurb_urb_file%add_var('rah_facade', 'wall and window aerodynamic resistance for heat combined', 's/m', 'tt0t')
    call slurb_urb_file%add_var('us_urb', 'friction velocity', 'm/s', 'tt0t')
    call slurb_urb_file%add_var('wall_hor_a_ratio', 'wall-to-horizontal area ratio', '-', 'tt0t')

    if (output_slurb_bc) then
      call slurb_urb_file%add_var('pt1', 'potential temperature', 'K', 'tt0t')
      call slurb_urb_file%add_var('q1', 'specific humidity', 'kg/kg', 'tt0t')
      call slurb_urb_file%add_var('uv_abs1', 'horizontal wind speed', 'm/s', 'tt0t')
      call slurb_urb_file%add_var('uv_eff1', 'effective horizontal wind speed', 'm/s', 'tt0t')
      call slurb_urb_file%add_var('vpt1', 'virtual potential temperature', 'K', 'tt0t')
    end if

    if (slurb_cross_output_tendencies) then
      call slurb_urb_file%add_var('tq_can', 'canyon mixing ratio tendency', 'kg/kg/s', 'tt0t')
      call slurb_urb_file%add_var('tt_can', 'canyon temperature tendency', 'W/m^3', 'tt0t')
    end if

    if (slurb_cross_output_radiation) then
      call slurb_urb_file%add_var('rad_lw_in_urb', 'incoming longwave radiation', 'W/m^2', 'tt0t')
      call slurb_urb_file%add_var('rad_lw_out_urb', 'outgoing longwave radiation', 'W/m^2', 'tt0t')
      call slurb_urb_file%add_var('rad_sw_in_urb', 'incoming shortwave radiation', 'W/m^2', 'tt0t')
      call slurb_urb_file%add_var('rad_sw_out_urb', 'outgoing shortwave radiation', 'W/m^2', 'tt0t')
      call slurb_urb_file%add_var('rad_lw_net_can', 'net longwave radiative at canyon top', 'W/m^2', 'tt0t')
      call slurb_urb_file%add_var('rad_lw_net_urb', 'urban aggegated net longwave radiative flux', 'W/m^2', 'tt0t')
      call slurb_urb_file%add_var('rad_sw_net_urb', 'urban aggegated net shortwave radiative flux', 'W/m^2', 'tt0t')
      call slurb_urb_file%add_var('sw_ref_denom', 'SW radiation reflection denominator', '', 'tt0t')
    end if

    if (slurb_cross_output_roof) then
      slurb_roof_file = slurb_3d_file_t('slurbcross_roof', nzs=nzb_roof-nzt_roof+1, lgpu=.false.)
      call add_output_file(slurb_roof_file, dtav_slurb, slurb_roof_file_id)
      slurb_roof_enabled = .true.
      call slurb_roof_file%add_var('t_roof', 'temperature roof', 'K', 'tttts_slurb')
      if (slurb_cross_output_tendencies) call slurb_roof_file%add_var('tt_roof', 'tendency roof', 'K/s', 'tttts_slurb')
    end if

    if (slurb_cross_output_road) then
      slurb_road_file = slurb_3d_file_t('slurbcross_road', nzs=nzb_road-nzt_road+1, lgpu=.false.)
      call add_output_file(slurb_road_file, dtav_slurb, slurb_road_file_id)
      slurb_road_enabled = .true.
      call slurb_road_file%add_var('t_road', 'temperature road', 'K', 'tttts_slurb')
      if (slurb_cross_output_tendencies) call slurb_road_file%add_var('tt_road', 'tendency road', 'K/s', 'tttts_slurb')
    end if

    if (slurb_cross_output_wall_win) then
      slurb_wall_file = slurb_3d_file_t('slurbcross_wall', nzs=nzb_wall-nzt_wall+1, lgpu=.false.)
      call add_output_file(slurb_wall_file, dtav_slurb, slurb_wall_file_id)
      slurb_wall_enabled = .true.
      call slurb_wall_file%add_var('t_wall_a', 'temperature wall a', 'K', 'tttts_slurb')
      call slurb_wall_file%add_var('t_wall_b', 'temperature wall b', 'K', 'tttts_slurb')
      if (slurb_cross_output_tendencies) then
        call slurb_wall_file%add_var('tt_wall_a', 'tendency wall a', 'K/s', 'tttts_slurb')
        call slurb_wall_file%add_var('tt_wall_b', 'tendency wall b', 'K/s', 'tttts_slurb')
      end if

      slurb_win_file = slurb_3d_file_t('slurbcross_win', nzs=nzb_win-nzt_win+1, lgpu=.false.)
      call add_output_file(slurb_win_file, dtav_slurb, slurb_win_file_id)
      slurb_win_enabled = .true.
      call slurb_win_file%add_var('t_win_a', 'temperature win a', 'K', 'tttts_slurb')
      call slurb_win_file%add_var('t_win_b', 'temperature win b', 'K', 'tttts_slurb')
      if (slurb_cross_output_tendencies) then
        call slurb_win_file%add_var('tt_win_a', 'tendency win a', 'K/s', 'tttts_slurb')
        call slurb_win_file%add_var('tt_win_b', 'tendency win b', 'K/s', 'tttts_slurb')
      end if
    end if
  end subroutine initslurbcrosssection

  subroutine slurbcrosssection
    use modglobal,    only : i1, j1, cp
    use modfields,    only : rhof
    use modslurbdata, only : slurb_tile, facade_rah_doe, enable_slurb, output_slurb_bc, slurb_cross_output, &
                             slurb_cross_output_tendencies, slurb_cross_output_radiation
    use modstat_nc,   only : lnetcdf

    implicit none

    real(field_r), pointer :: albedo_urb(:,:), ol_urb(:,:), qsws_urb(:,:), rad_lw_in_urb(:,:), &
                              rad_lw_out_urb(:,:), rad_sw_in_urb(:,:), rad_sw_out_urb(:,:), ram_urb(:,:), rib_urb(:,:), &
                              shf_urb(:,:), t_2m_urb(:,:), t_c_urb(:,:), t_h_urb(:,:), thl_rad_urb(:,:), usws_urb(:,:), &
                              vsws_urb(:,:), thlskin(:,:), qtskin(:,:), q_can_0(:,:), t_can_0(:,:), &
                              tq_can(:,:), tt_can(:,:), shf_can(:,:), shf_external(:,:), shf_traffic(:,:), &
                              qsws_can(:,:), qsws_external(:,:), rad_lw_net_can(:,:), rad_lw_net_urb(:,:), rad_sw_net_urb(:,:), &
                              ol_can(:,:), pt_can(:,:), rib_can(:,:), us_can(:,:), uv_abs_can(:,:), uv_eff_can(:,:), &
                              vpt_can(:,:), rah_can(:,:), rah_facade(:,:), pt1(:,:), q1(:,:), us_urb(:,:), uv_abs1(:,:), &
                              uv_eff1(:,:), vpt1(:,:), sw_ref_denom(:,:), &
                              wall_hor_a_ratio(:,:), z_mo(:,:)
    real(field_r), pointer :: tt_roof(:,:,:), t_roof(:,:,:)
    real(field_r), pointer :: tt_road(:,:,:), t_road(:,:,:)
    real(field_r), pointer :: tt_wall_a(:,:,:), tt_wall_b(:,:,:), t_wall_a(:,:,:), t_wall_b(:,:,:)
    real(field_r), pointer :: tt_win_a(:,:,:), tt_win_b(:,:,:), &
                              t_win_a(:,:,:), t_win_b(:,:,:)

    if (.not. (lnetcdf .and. enable_slurb .and. slurb_cross_output)) return

    if (slurb_urb_enabled .and. is_sampling_timestep(slurb_urb_file_id)) then
      call slurb_urb_file%get_pointer('albedo_urb', albedo_urb)
      call slurb_urb_file%get_pointer('ol_urb', ol_urb)
      call slurb_urb_file%get_pointer('qsws_urb', qsws_urb)
      call slurb_urb_file%get_pointer('ram_urb', ram_urb)
      call slurb_urb_file%get_pointer('rib_urb', rib_urb)
      call slurb_urb_file%get_pointer('shf_urb', shf_urb)
      call slurb_urb_file%get_pointer('t_2m_urb', t_2m_urb)
      call slurb_urb_file%get_pointer('t_c_urb', t_c_urb)
      call slurb_urb_file%get_pointer('t_h_urb', t_h_urb)
      call slurb_urb_file%get_pointer('thl_rad_urb', thl_rad_urb)
      call slurb_urb_file%get_pointer('usws_urb', usws_urb)
      call slurb_urb_file%get_pointer('vsws_urb', vsws_urb)
      call slurb_urb_file%get_pointer('thlskin', thlskin)
      call slurb_urb_file%get_pointer('qtskin', qtskin)
      call slurb_urb_file%get_pointer('q_can_0', q_can_0)
      call slurb_urb_file%get_pointer('t_can_0', t_can_0)
      call slurb_urb_file%get_pointer('shf_can', shf_can)
      call slurb_urb_file%get_pointer('shf_external', shf_external)
      call slurb_urb_file%get_pointer('shf_traffic', shf_traffic)
      call slurb_urb_file%get_pointer('qsws_can', qsws_can)
      call slurb_urb_file%get_pointer('qsws_external', qsws_external)
      call slurb_urb_file%get_pointer('ol_can', ol_can)
      call slurb_urb_file%get_pointer('pt_can', pt_can)
      call slurb_urb_file%get_pointer('rib_can', rib_can)
      call slurb_urb_file%get_pointer('us_can', us_can)
      call slurb_urb_file%get_pointer('uv_abs_can', uv_abs_can)
      call slurb_urb_file%get_pointer('uv_eff_can', uv_eff_can)
      call slurb_urb_file%get_pointer('vpt_can', vpt_can)
      call slurb_urb_file%get_pointer('rah_can', rah_can)
      call slurb_urb_file%get_pointer('rah_facade', rah_facade)
      call slurb_urb_file%get_pointer('us_urb', us_urb)
            if (output_slurb_bc) then
              call slurb_urb_file%get_pointer('pt1', pt1)
              call slurb_urb_file%get_pointer('q1', q1)
              call slurb_urb_file%get_pointer('uv_abs1', uv_abs1)
              call slurb_urb_file%get_pointer('uv_eff1', uv_eff1)
              call slurb_urb_file%get_pointer('vpt1', vpt1)
            end if

      if (slurb_cross_output_tendencies) then
        call slurb_urb_file%get_pointer('tq_can', tq_can)
        call slurb_urb_file%get_pointer('tt_can', tt_can)
      end if

      if (slurb_cross_output_radiation) then
        call slurb_urb_file%get_pointer('rad_lw_in_urb', rad_lw_in_urb)
        call slurb_urb_file%get_pointer('rad_lw_out_urb', rad_lw_out_urb)
        call slurb_urb_file%get_pointer('rad_sw_in_urb', rad_sw_in_urb)
        call slurb_urb_file%get_pointer('rad_sw_out_urb', rad_sw_out_urb)
        call slurb_urb_file%get_pointer('rad_lw_net_can', rad_lw_net_can)
        call slurb_urb_file%get_pointer('rad_lw_net_urb', rad_lw_net_urb)
        call slurb_urb_file%get_pointer('rad_sw_net_urb', rad_sw_net_urb)
        call slurb_urb_file%get_pointer('sw_ref_denom', sw_ref_denom)
      end if

      albedo_urb(:,:) = slurb_tile%albedo_urb(2:i1,2:j1)
      ol_urb(:,:) = slurb_tile%ol_urb(2:i1,2:j1)
      qsws_urb(:,:) = slurb_tile%qsws_urb(2:i1,2:j1)
      ram_urb(:,:) = slurb_tile%ram_urb(2:i1,2:j1)
      rib_urb(:,:) = slurb_tile%rib_urb(2:i1,2:j1)
      shf_urb(:,:) = slurb_tile%shf_urb(2:i1,2:j1)
      t_2m_urb(:,:) = slurb_tile%t_2m_urb(2:i1,2:j1)
      t_c_urb(:,:) = slurb_tile%t_c_urb(2:i1,2:j1)
      t_h_urb(:,:) = slurb_tile%t_h_urb(2:i1,2:j1)
      thl_rad_urb(:,:) = slurb_tile%thl_rad_urb(2:i1,2:j1)
      usws_urb(:,:) = slurb_tile%usws_urb(2:i1,2:j1)
      vsws_urb(:,:) = slurb_tile%vsws_urb(2:i1,2:j1)
      thlskin(:,:) = slurb_tile%thlskin(2:i1,2:j1)
      qtskin(:,:) = slurb_tile%qtskin(2:i1,2:j1)
      q_can_0(:,:) = slurb_tile%q_can_0(2:i1,2:j1)
      t_can_0(:,:) = slurb_tile%t_can_0(2:i1,2:j1)
      shf_can(:,:) = slurb_tile%shf_can(2:i1,2:j1)
      shf_external(:,:) = slurb_tile%shf_external(2:i1,2:j1)
      shf_traffic(:,:) = 0.0_field_r
      qsws_can(:,:) = slurb_tile%qsws_can(2:i1,2:j1)
      qsws_external(:,:) = slurb_tile%qsws_external(2:i1,2:j1)
      ol_can(:,:) = slurb_tile%ol_can(2:i1,2:j1)
      pt_can(:,:) = slurb_tile%pt_can(2:i1,2:j1)
      rib_can(:,:) = slurb_tile%rib_can(2:i1,2:j1)
      us_can(:,:) = slurb_tile%us_can(2:i1,2:j1)
      uv_abs_can(:,:) = slurb_tile%uv_abs_can(2:i1,2:j1)
      uv_eff_can(:,:) = slurb_tile%uv_eff_can(2:i1,2:j1)
      vpt_can(:,:) = slurb_tile%vpt_can(2:i1,2:j1)
      rah_can(:,:) = slurb_tile%rah_can(2:i1,2:j1)
      if (facade_rah_doe) then
        rah_facade(:,:) = 0.0_field_r
      else
        rah_facade(:,:) = slurb_tile%rah_facade(2:i1,2:j1)
      end if
      us_urb(:,:) = slurb_tile%us_urb(2:i1,2:j1)

      if (output_slurb_bc) then
        pt1(:,:) = slurb_tile%pt1(2:i1,2:j1)
        q1(:,:) = slurb_tile%q1(2:i1,2:j1)
        uv_abs1(:,:) = slurb_tile%uv_abs1(2:i1,2:j1)
        uv_eff1(:,:) = slurb_tile%uv_eff1(2:i1,2:j1)
        vpt1(:,:) = slurb_tile%vpt1(2:i1,2:j1)
      end if

      if (slurb_cross_output_tendencies) then
        tq_can(:,:) = slurb_tile%tq_can(2:i1,2:j1) * rhof(1)
        tt_can(:,:) = slurb_tile%tt_can(2:i1,2:j1) * cp * rhof(1)
      end if

      if (slurb_cross_output_radiation) then
        rad_lw_in_urb(:,:) = slurb_tile%rad_lw_in_urb(2:i1,2:j1)
        rad_lw_out_urb(:,:) = slurb_tile%rad_lw_out_urb(2:i1,2:j1)
        rad_sw_in_urb(:,:) = slurb_tile%rad_sw_in_urb(2:i1,2:j1)
        rad_sw_out_urb(:,:) = slurb_tile%rad_sw_out_urb(2:i1,2:j1)
        rad_lw_net_can(:,:) = slurb_tile%rad_lw_net_can(2:i1,2:j1)
        rad_lw_net_urb(:,:) = slurb_tile%rad_lw_net_urb(2:i1,2:j1)
        rad_sw_net_urb(:,:) = slurb_tile%rad_sw_net_urb(2:i1,2:j1)
        sw_ref_denom(:,:) = slurb_tile%sw_ref_denom(2:i1,2:j1)
      end if
    end if

    if (slurb_roof_enabled .and. is_sampling_timestep(slurb_roof_file_id)) then
      call slurb_roof_file%get_pointer('t_roof', t_roof)

      t_roof(:,:,:) = slurb_tile%t_roof_0(:,2:i1,2:j1)

      if (slurb_cross_output_tendencies) then
        call slurb_roof_file%get_pointer('tt_roof', tt_roof)
        tt_roof(:,:,:) = slurb_tile%tt_roof(:,2:i1,2:j1)
      end if
    end if

    if (slurb_road_enabled .and. is_sampling_timestep(slurb_road_file_id)) then
      call slurb_road_file%get_pointer('t_road', t_road)

      t_road(:,:,:) = slurb_tile%t_road_0(:,2:i1,2:j1)

      if (slurb_cross_output_tendencies) then
        call slurb_road_file%get_pointer('tt_road', tt_road)
        tt_road(:,:,:) = slurb_tile%tt_road(:,2:i1,2:j1)
      end if
    end if

    if (slurb_wall_enabled .and. is_sampling_timestep(slurb_wall_file_id)) then
      call slurb_wall_file%get_pointer('t_wall_a', t_wall_a)
      call slurb_wall_file%get_pointer('t_wall_b', t_wall_b)

      t_wall_a(:,:,:) = slurb_tile%t_wall_a_0(:,2:i1,2:j1)
      t_wall_b(:,:,:) = slurb_tile%t_wall_b_0(:,2:i1,2:j1)

      if (slurb_cross_output_tendencies) then
        call slurb_wall_file%get_pointer('tt_wall_a', tt_wall_a)
        call slurb_wall_file%get_pointer('tt_wall_b', tt_wall_b)
        tt_wall_a(:,:,:) = slurb_tile%tt_wall_a(:,2:i1,2:j1)
        tt_wall_b(:,:,:) = slurb_tile%tt_wall_b(:,2:i1,2:j1)
      end if
    end if

    if (slurb_win_enabled .and. is_sampling_timestep(slurb_win_file_id)) then
      call slurb_win_file%get_pointer('t_win_a', t_win_a)
      call slurb_win_file%get_pointer('t_win_b', t_win_b)

      t_win_a(:,:,:) = slurb_tile%t_win_a_0(:,2:i1,2:j1)
      t_win_b(:,:,:) = slurb_tile%t_win_b_0(:,2:i1,2:j1)

      if (slurb_cross_output_tendencies) then
        call slurb_win_file%get_pointer('tt_win_a', tt_win_a)
        call slurb_win_file%get_pointer('tt_win_b', tt_win_b)
        tt_win_a(:,:,:) = slurb_tile%tt_win_a(:,2:i1,2:j1)
        tt_win_b(:,:,:) = slurb_tile%tt_win_b(:,2:i1,2:j1)
      end if

    end if
  end subroutine slurbcrosssection

end module modslurbcrosssection
