!> \file modslurbhelpers.f90
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
module modslurbhelpers
    use modslurbdata

    contains

subroutine slurb_read_namelist(nml_filename)
    use modglobal,   only : ifnamopt, checknamelisterror
    use modmpi,      only : myid, comm3d, mpierr, D_MPI_BCAST
    use fortran_support,       only: nnml_output
    implicit none

    character(len=*), intent(in) :: nml_filename


    integer :: ierr

    ! Namelist definition
    namelist /NAMSLURB/ &
        urban_fraction, urban_roughness_length, building_plan_area_fraction, building_frontal_area_fraction, building_height, window_fraction,&
        street_canyon_aspect_ratio, building_type, pavement_type, anisotropic_street_canyons, street_canyon_orientation, deep_soil_temperature,building_indoor_temperature,shf_external,qsws_external

    ! Read namelist
    if (myid == 0) then
        open(ifnamopt, file=nml_filename, status='old', iostat=ierr)
        read(ifnamopt, NAMSLURB, iostat=ierr)
        call checknamelisterror(ierr, ifnamopt, 'NAMSLURB')
        write(nnml_output, NAMSLURB)
        close(ifnamopt)
    end if

    ! Broadcast namelist values to all MPI tasks
    call D_MPI_BCAST(urban_fraction,  1, 0, comm3d, mpierr)
    call D_MPI_BCAST(urban_roughness_length,   1, 0, comm3d, mpierr)
    call D_MPI_BCAST(building_plan_area_fraction,            1, 0, comm3d, mpierr)
    call D_MPI_BCAST(building_frontal_area_fraction,             1, 0, comm3d, mpierr)
    call D_MPI_BCAST(building_height,       1, 0, comm3d, mpierr)
    call D_MPI_BCAST(window_fraction,   1, 0, comm3d, mpierr)
    call D_MPI_BCAST(street_canyon_aspect_ratio,       1, 0, comm3d, mpierr)
    call D_MPI_BCAST(building_type, 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(pavement_type, 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(anisotropic_street_canyons, 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(street_canyon_orientation, 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(deep_soil_temperature, 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(building_indoor_temperature, 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(shf_external, 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(qsws_external, 1, 0, comm3d, mpierr)
end subroutine slurb_read_namelist

subroutine slurb_bulk_allocations
    use modglobal, only: i2, j2
    use, intrinsic :: IEEE_ARITHMETIC

    allocate(fraction_slurb(i2,j2))

    allocate(ln_z_z0_roof(i2,j2))
    allocate(ln_z_z0h_roof(i2,j2))
    allocate(ln_z_z0_urb(i2,j2))
    allocate(pt_surface(i2,j2))
    allocate(ln_z_z0_road(i2,j2))
    allocate(ln_z_z0h_road(i2,j2))
    !-- Bulk allocation
    ALLOCATE( slurb_tile%dz_roof(nzt_roof:nzb_roof,i2,j2) )
    ALLOCATE( slurb_tile%dz_wall(nzt_wall:nzb_wall,i2,j2) )
    ALLOCATE( slurb_tile%dz_road(nzt_road:nzb_road,i2,j2) )
    ALLOCATE( slurb_tile%dz_win(nzt_win:nzb_win,i2,j2) )
    ALLOCATE( slurb_tile%zw_win(nzt_win:nzb_win,i2,j2) )

    ALLOCATE( slurb_tile%t_c_urb(i2,j2) )
    ALLOCATE( slurb_tile%thl_rad_urb(i2,j2) )
    ALLOCATE( slurb_tile%t_h_urb(i2,j2) )
    ALLOCATE( slurb_tile%t_2m_urb(i2,j2) )
    ALLOCATE( slurb_tile%shf_urb(i2,j2) )
    ALLOCATE( slurb_tile%qsws_urb(i2,j2) )
    ALLOCATE( slurb_tile%ol_urb(i2,j2) )
    ALLOCATE( slurb_tile%rib_urb(i2,j2) )
    ALLOCATE( slurb_tile%ram_urb(i2,j2) )
    ALLOCATE( slurb_tile%usws_urb(i2,j2) )
    ALLOCATE( slurb_tile%vsws_urb(i2,j2) )
    ALLOCATE( slurb_tile%thlskin(i2,j2) )
    ALLOCATE( slurb_tile%qtskin(i2,j2) )

    ALLOCATE( slurb_tile%albedo_urb(i2,j2) )
    ALLOCATE( slurb_tile%emiss_urb(i2,j2) )

    ALLOCATE( slurb_tile%t_indoor(i2,j2) )
    ALLOCATE( slurb_tile%t_soil(i2,j2) )

    ALLOCATE( slurb_tile%tt_can(i2,j2) )
    ALLOCATE( slurb_tile%tt_wall_a(nzt_wall:nzb_wall,i2,j2) )
    ALLOCATE( slurb_tile%tt_wall_b(nzt_wall:nzb_wall,i2,j2) )
    ALLOCATE( slurb_tile%tt_win_a(nzt_win:nzb_win,i2,j2) )
    ALLOCATE( slurb_tile%tt_win_b(nzt_win:nzb_win,i2,j2) )
    ALLOCATE( slurb_tile%tt_roof(nzt_roof:nzb_roof,i2,j2) )
    ALLOCATE( slurb_tile%tt_road(nzt_road:nzb_road,i2,j2) )

    ALLOCATE( slurb_tile%pt_wall_a(i2,j2) )
    ALLOCATE( slurb_tile%pt_wall_b(i2,j2) )
    ALLOCATE( slurb_tile%pt_win_a(i2,j2) )
    ALLOCATE( slurb_tile%pt_win_b(i2,j2) )
    ALLOCATE( slurb_tile%pt_roof(i2,j2) )
    ALLOCATE( slurb_tile%pt_road(i2,j2) )

    ALLOCATE( slurb_tile%shf_can(i2,j2) )
    ALLOCATE( slurb_tile%shf_roof(i2,j2) )
    ALLOCATE( slurb_tile%shf_road(i2,j2) )
    ALLOCATE( slurb_tile%shf_wall_a(i2,j2) )
    ALLOCATE( slurb_tile%shf_wall_b(i2,j2) )
    ALLOCATE( slurb_tile%shf_win_a(i2,j2) )
    ALLOCATE( slurb_tile%shf_win_b(i2,j2) )

    ALLOCATE( slurb_tile%shf_external(i2,j2) )
    ALLOCATE( slurb_tile%shf_traffic(i2,j2) )

    ALLOCATE( slurb_tile%ghf_road(i2,j2) )
    ALLOCATE( slurb_tile%ghf_roof(i2,j2) )
    ALLOCATE( slurb_tile%ghf_wall_a(i2,j2) )
    ALLOCATE( slurb_tile%ghf_wall_b(i2,j2) )
    ALLOCATE( slurb_tile%ghf_win_a(i2,j2) )
    ALLOCATE( slurb_tile%ghf_win_b(i2,j2) )

    ALLOCATE( slurb_tile%rad_lw_in_urb(i2,j2) )
    ALLOCATE( slurb_tile%rad_sw_in_urb(i2,j2) )
    ALLOCATE( slurb_tile%rad_lw_out_urb(i2,j2) )
    ALLOCATE( slurb_tile%rad_sw_out_urb(i2,j2) )

    ALLOCATE( slurb_tile%rad_lw_net_urb(i2,j2) )
    ALLOCATE( slurb_tile%rad_sw_net_urb(i2,j2) )

    ALLOCATE( slurb_tile%rad_lw_net_can(i2,j2) )

    ALLOCATE( slurb_tile%rad_lw_net_roof(i2,j2) )
    ALLOCATE( slurb_tile%rad_sw_net_roof(i2,j2) )
    ALLOCATE( slurb_tile%rad_lw_net_road(i2,j2) )
    ALLOCATE( slurb_tile%rad_sw_net_road(i2,j2) )
    ALLOCATE( slurb_tile%rad_sw_in_road(i2,j2) )
    ALLOCATE( slurb_tile%rad_lw_net_wall_a(i2,j2) )
    ALLOCATE( slurb_tile%rad_sw_net_wall_a(i2,j2) )
    ALLOCATE( slurb_tile%rad_lw_net_wall_b(i2,j2) )
    ALLOCATE( slurb_tile%rad_sw_net_wall_b(i2,j2) )
    ALLOCATE( slurb_tile%rad_lw_net_win_a(i2,j2) )
    ALLOCATE( slurb_tile%rad_sw_net_win_a(i2,j2) )
    ALLOCATE( slurb_tile%rad_sw_in_win_a(i2,j2) )
    ALLOCATE( slurb_tile%rad_lw_net_win_b(i2,j2) )
    ALLOCATE( slurb_tile%rad_sw_net_win_b(i2,j2) )
    ALLOCATE( slurb_tile%rad_sw_in_win_b(i2,j2) )

    ALLOCATE( slurb_tile%pt_can(i2,j2) )
    ALLOCATE( slurb_tile%uv_abs_can(i2,j2) )
    ALLOCATE( slurb_tile%uv_eff_can(i2,j2) )
    ALLOCATE( slurb_tile%us_can(i2,j2) )
    ALLOCATE( slurb_tile%rib_can(i2,j2) )
    ALLOCATE( slurb_tile%ol_can(i2,j2) )

    ALLOCATE( slurb_tile%rib_roof(i2,j2) )
    ALLOCATE( slurb_tile%ol_roof(i2,j2) )
    ALLOCATE( slurb_tile%rib_road(i2,j2) )
    ALLOCATE( slurb_tile%ol_road(i2,j2) )

    ALLOCATE( slurb_tile%us_roof(i2,j2) )
    ALLOCATE( slurb_tile%us_road(i2,j2) )

    ALLOCATE( slurb_tile%hw_can(i2,j2) )
    ALLOCATE( slurb_tile%anisotropic_canyon(i2,j2) )
    ALLOCATE( slurb_tile%theta_can(i2,j2) )
    ALLOCATE( slurb_tile%h_bld(i2,j2) )
    ALLOCATE( slurb_tile%f_bld(i2,j2) )
    ALLOCATE( slurb_tile%f_bld_frn(i2,j2) )
    ALLOCATE( slurb_tile%f_win(i2,j2) )
    ALLOCATE( slurb_tile%svf_road(i2,j2) )
    ALLOCATE( slurb_tile%svf_wall(i2,j2) )
    ALLOCATE( slurb_tile%z0_urb(i2,j2) )

    ALLOCATE( slurb_tile%rah_roof(i2,j2) )
    ALLOCATE( slurb_tile%rah_road(i2,j2) )
    ALLOCATE( slurb_tile%rah_can(i2,j2) )

    IF ( facade_rah_doe )  THEN
       ALLOCATE( slurb_tile%rah_wall_a(i2,j2) )
       ALLOCATE( slurb_tile%rah_wall_b(i2,j2) )
       ALLOCATE( slurb_tile%rah_win_a(i2,j2) )
       ALLOCATE( slurb_tile%rah_win_b(i2,j2) )
    ELSE
       ALLOCATE( slurb_tile%rah_facade(i2,j2) )
    ENDIF

    ALLOCATE( slurb_tile%lambda_roof(nzt_roof:nzb_roof,i2,j2) )
    ALLOCATE( slurb_tile%c_roof(nzt_roof:nzb_roof,i2,j2) )
    ALLOCATE( slurb_tile%albedo_roof(i2,j2) )
    ALLOCATE( slurb_tile%emiss_roof(i2,j2) )
    ALLOCATE( slurb_tile%z0_roof(i2,j2) )
    ALLOCATE( slurb_tile%z0h_roof(i2,j2) )
    ALLOCATE( slurb_tile%lambda_wall(nzt_wall:nzb_wall,i2,j2) )
    ALLOCATE( slurb_tile%c_wall(nzt_wall:nzb_wall,i2,j2) )
    ALLOCATE( slurb_tile%albedo_wall(i2,j2) )
    ALLOCATE( slurb_tile%emiss_wall(i2,j2) )
    ALLOCATE( slurb_tile%z0_wall(i2,j2) )
    ALLOCATE( slurb_tile%lambda_win(nzt_win:nzb_win,i2,j2) )
    ALLOCATE( slurb_tile%c_win(nzt_win:nzb_win,i2,j2) )
    ALLOCATE( slurb_tile%albedo_wall_win(i2,j2) )
    ALLOCATE( slurb_tile%albedo_win(i2,j2) )
    ALLOCATE( slurb_tile%emiss_win(i2,j2) )
    ALLOCATE( slurb_tile%transmissivity_win(i2,j2) )
    ALLOCATE( slurb_tile%absorption_win(nzt_win:nzb_win,i2,j2) )
    ALLOCATE( slurb_tile%lambda_road(nzt_road:nzb_road,i2,j2) )
    ALLOCATE( slurb_tile%c_road(nzt_road:nzb_road,i2,j2) )
    ALLOCATE( slurb_tile%albedo_road(i2,j2) )
    ALLOCATE( slurb_tile%emiss_road(i2,j2) )
    ALLOCATE( slurb_tile%z0_road(i2,j2) )
    ALLOCATE( slurb_tile%z0h_road(i2,j2) )

    ALLOCATE( slurb_tile%conductivity_roof(nzt_roof:nzb_roof,i2,j2) )
    ALLOCATE( slurb_tile%conductivity_wall(nzt_wall:nzb_wall,i2,j2) )
    ALLOCATE( slurb_tile%conductivity_win(nzt_win:nzb_win,i2,j2) )
    ALLOCATE( slurb_tile%conductivity_road(nzt_road:nzb_road,i2,j2) )

    ALLOCATE( slurb_tile%z_mo(i2,j2) )
    ALLOCATE( slurb_tile%z_mo_can(i2,j2) )
    ALLOCATE( slurb_tile%uv_abs_can_coef(i2,j2) )
    ALLOCATE( slurb_tile%wall_hor_a_ratio(i2,j2) )


    ALLOCATE( slurb_tile%lw_roof_coef(1:2,i2,j2) )
    ALLOCATE( slurb_tile%lw_road_coef(1:4,i2,j2) )
    ALLOCATE( slurb_tile%lw_wall_coef(1:6,i2,j2) )
    ALLOCATE( slurb_tile%lw_win_coef(1:6,i2,j2) )
    ALLOCATE( slurb_tile%sw_ref_denom(i2,j2) )

    ALLOCATE( slurb_tile%us_urb(i2,j2) )
    ALLOCATE( slurb_tile%uv_eff1(i2,j2) )
    ALLOCATE( slurb_tile%uv_abs1(i2,j2) )
    ALLOCATE( slurb_tile%pt1(i2,j2) )

    ALLOCATE( slurb_tile%t_can_0(i2,j2) )
    ALLOCATE( slurb_tile%t_can_m(i2,j2) )
    ALLOCATE( slurb_tile%t_wall_a_0(nzt_wall:nzb_wall,i2,j2) )
    ALLOCATE( slurb_tile%t_wall_a_m(nzt_wall:nzb_wall,i2,j2) )
    ALLOCATE( slurb_tile%t_wall_b_0(nzt_wall:nzb_wall,i2,j2) )
    ALLOCATE( slurb_tile%t_wall_b_m(nzt_wall:nzb_wall,i2,j2) )
    ALLOCATE( slurb_tile%t_win_a_0(nzt_win:nzb_win,i2,j2) )
    ALLOCATE( slurb_tile%t_win_a_m(nzt_win:nzb_win,i2,j2) )
    ALLOCATE( slurb_tile%t_win_b_0(nzt_win:nzb_win,i2,j2) )
    ALLOCATE( slurb_tile%t_win_b_m(nzt_win:nzb_win,i2,j2) )
    ALLOCATE( slurb_tile%t_roof_0(nzt_roof:nzb_roof,i2,j2) )
    ALLOCATE( slurb_tile%t_roof_m(nzt_roof:nzb_roof,i2,j2) )
    ALLOCATE( slurb_tile%t_road_0(nzt_road:nzb_road,i2,j2) )
    ALLOCATE( slurb_tile%t_road_m(nzt_road:nzb_road,i2,j2) )

    IF ( moist_physics )  THEN
       ALLOCATE( slurb_tile%tq_can(i2,j2))
       ALLOCATE( slurb_tile%tm_liq_roof(i2,j2) )
       ALLOCATE( slurb_tile%tm_liq_road(i2,j2) )
       ALLOCATE( slurb_tile%tm_roof_runoff(i2,j2) )
       ALLOCATE( slurb_tile%tm_road_runoff(i2,j2) )
       ALLOCATE( slurb_tile%tm_roof_precep(i2,j2) )
       ALLOCATE( slurb_tile%tm_road_precep(i2,j2) )

       ALLOCATE( slurb_tile%vpt_roof(i2,j2) )
       ALLOCATE( slurb_tile%vpt_road(i2,j2) )

       ALLOCATE( slurb_tile%q_roof(i2,j2) )
       ALLOCATE( slurb_tile%q_road(i2,j2) )
       ALLOCATE( slurb_tile%qs_roof(i2,j2) )
       ALLOCATE( slurb_tile%qs_road(i2,j2) )

       ALLOCATE( slurb_tile%qsws_can(i2,j2) )
       ALLOCATE( slurb_tile%qsws_roof(i2,j2) )
       ALLOCATE( slurb_tile%qsws_road(i2,j2) )
       ALLOCATE( slurb_tile%qsws_liq_roof(i2,j2) )
       ALLOCATE( slurb_tile%qsws_liq_road(i2,j2) )

       ALLOCATE( slurb_tile%c_liq_roof(i2,j2) )
       ALLOCATE( slurb_tile%c_liq_road(i2,j2) )

       ALLOCATE( slurb_tile%vpt_can(i2,j2) )

       ALLOCATE( slurb_tile%q1(i2,j2) )
       ALLOCATE( slurb_tile%vpt1(i2,j2) )

       ALLOCATE( slurb_tile%qsws_external(i2,j2) )

       ALLOCATE( slurb_tile%q_can_0(i2,j2) )
       ALLOCATE( slurb_tile%q_can_m(i2,j2) )
       ALLOCATE( slurb_tile%m_liq_roof_0(i2,j2) )
       ALLOCATE( slurb_tile%m_liq_roof_m(i2,j2) )
       ALLOCATE( slurb_tile%m_liq_road_0(i2,j2) )
       ALLOCATE( slurb_tile%m_liq_road_m(i2,j2) )

    ENDIF

    ALLOCATE( slurb_tile%dt_max(i2,j2) )



#ifndef __FUJITSU ! Fujitsu compiler doesn't like these initializations (March 2026)
    fraction_slurb(:,:) = ieee_value(fraction_slurb,ieee_signaling_nan)

    ln_z_z0_roof(:,:) = ieee_value(ln_z_z0_roof,ieee_signaling_nan)
    ln_z_z0h_roof(:,:) = ieee_value(ln_z_z0h_roof,ieee_signaling_nan)
    ln_z_z0_urb(:,:) = ieee_value(ln_z_z0_urb,ieee_signaling_nan)
    pt_surface(:,:) = ieee_value(pt_surface,ieee_signaling_nan)
    ln_z_z0_road(:,:) = ieee_value(ln_z_z0_road,ieee_signaling_nan)
    ln_z_z0h_road(:,:) = ieee_value(ln_z_z0h_road,ieee_signaling_nan)
    !-- Bulk allocation
    slurb_tile%dz_roof(:,:,:) = ieee_value(slurb_tile%dz_roof,ieee_signaling_nan)
    slurb_tile%dz_wall(:,:,:) = ieee_value(slurb_tile%dz_wall,ieee_signaling_nan)
    slurb_tile%dz_road(:,:,:) = ieee_value(slurb_tile%dz_road,ieee_signaling_nan)
    slurb_tile%dz_win(:,:,:) = ieee_value(slurb_tile%dz_win,ieee_signaling_nan)
    slurb_tile%zw_win(:,:,:) = ieee_value(slurb_tile%zw_win,ieee_signaling_nan)

    slurb_tile%t_c_urb(:,:) = ieee_value(slurb_tile%t_c_urb,ieee_signaling_nan)
    slurb_tile%thl_rad_urb(:,:) = ieee_value(slurb_tile%thl_rad_urb,ieee_signaling_nan)
    slurb_tile%t_h_urb(:,:) = ieee_value(slurb_tile%t_h_urb,ieee_signaling_nan)
    slurb_tile%t_2m_urb(:,:) = ieee_value(slurb_tile%t_2m_urb,ieee_signaling_nan)
    slurb_tile%shf_urb(:,:) = ieee_value(slurb_tile%shf_urb,ieee_signaling_nan)
    slurb_tile%qsws_urb(:,:) = ieee_value(slurb_tile%qsws_urb,ieee_signaling_nan)
    slurb_tile%ol_urb(:,:) = ieee_value(slurb_tile%ol_urb,ieee_signaling_nan)
    slurb_tile%rib_urb(:,:) = ieee_value(slurb_tile%rib_urb,ieee_signaling_nan)
    slurb_tile%ram_urb(:,:) = ieee_value(slurb_tile%ram_urb,ieee_signaling_nan)
    slurb_tile%usws_urb(:,:) = ieee_value(slurb_tile%usws_urb,ieee_signaling_nan)
    slurb_tile%vsws_urb(:,:) = ieee_value(slurb_tile%vsws_urb,ieee_signaling_nan)
    slurb_tile%thlskin(:,:) = ieee_value(slurb_tile%thlskin,ieee_signaling_nan)
    slurb_tile%qtskin(:,:) = ieee_value(slurb_tile%qtskin,ieee_signaling_nan)

    slurb_tile%albedo_urb(:,:) = ieee_value(slurb_tile%albedo_urb,ieee_signaling_nan)
    slurb_tile%emiss_urb(:,:) = ieee_value(slurb_tile%emiss_urb,ieee_signaling_nan)

    slurb_tile%t_indoor(:,:) = ieee_value(slurb_tile%t_indoor,ieee_signaling_nan)
    slurb_tile%t_soil(:,:) = ieee_value(slurb_tile%t_soil,ieee_signaling_nan)

    slurb_tile%tt_can(:,:) = ieee_value(slurb_tile%tt_can,ieee_signaling_nan)
    slurb_tile%tt_wall_a(:,:,:) = ieee_value(slurb_tile%tt_wall_a,ieee_signaling_nan)
    slurb_tile%tt_wall_b(:,:,:) = ieee_value(slurb_tile%tt_wall_b,ieee_signaling_nan)
    slurb_tile%tt_win_a(:,:,:) = ieee_value(slurb_tile%tt_win_a,ieee_signaling_nan)
    slurb_tile%tt_win_b(:,:,:) = ieee_value(slurb_tile%tt_win_b,ieee_signaling_nan)
    slurb_tile%tt_roof(:,:,:) = ieee_value(slurb_tile%tt_roof,ieee_signaling_nan)
    slurb_tile%tt_road(:,:,:) = ieee_value(slurb_tile%tt_road,ieee_signaling_nan)

    slurb_tile%pt_wall_a(:,:) = ieee_value(slurb_tile%pt_wall_a,ieee_signaling_nan)
    slurb_tile%pt_wall_b(:,:) = ieee_value(slurb_tile%pt_wall_b,ieee_signaling_nan)
    slurb_tile%pt_win_a(:,:) = ieee_value(slurb_tile%pt_win_a,ieee_signaling_nan)
    slurb_tile%pt_win_b(:,:) = ieee_value(slurb_tile%pt_win_b,ieee_signaling_nan)
    slurb_tile%pt_roof(:,:) = ieee_value(slurb_tile%pt_roof,ieee_signaling_nan)
    slurb_tile%pt_road(:,:) = ieee_value(slurb_tile%pt_road,ieee_signaling_nan)

    slurb_tile%shf_can(:,:) = ieee_value(slurb_tile%shf_can,ieee_signaling_nan)
    slurb_tile%shf_roof(:,:) = ieee_value(slurb_tile%shf_roof,ieee_signaling_nan)
    slurb_tile%shf_road(:,:) = ieee_value(slurb_tile%shf_road,ieee_signaling_nan)
    slurb_tile%shf_wall_a(:,:) = ieee_value(slurb_tile%shf_wall_a,ieee_signaling_nan)
    slurb_tile%shf_wall_b(:,:) = ieee_value(slurb_tile%shf_wall_b,ieee_signaling_nan)
    slurb_tile%shf_win_a(:,:) = ieee_value(slurb_tile%shf_win_a,ieee_signaling_nan)
    slurb_tile%shf_win_b(:,:) = ieee_value(slurb_tile%shf_win_b,ieee_signaling_nan)

    slurb_tile%shf_external(:,:) = ieee_value(slurb_tile%shf_external,ieee_signaling_nan)
    slurb_tile%shf_traffic(:,:) = ieee_value(slurb_tile%shf_traffic,ieee_signaling_nan)

    slurb_tile%ghf_road(:,:) = ieee_value(slurb_tile%ghf_road,ieee_signaling_nan)
    slurb_tile%ghf_roof(:,:) = ieee_value(slurb_tile%ghf_roof,ieee_signaling_nan)
    slurb_tile%ghf_wall_a(:,:) = ieee_value(slurb_tile%ghf_wall_a,ieee_signaling_nan)
    slurb_tile%ghf_wall_b(:,:) = ieee_value(slurb_tile%ghf_wall_b,ieee_signaling_nan)
    slurb_tile%ghf_win_a(:,:) = ieee_value(slurb_tile%ghf_win_a,ieee_signaling_nan)
    slurb_tile%ghf_win_b(:,:) = ieee_value(slurb_tile%ghf_win_b,ieee_signaling_nan)

    slurb_tile%rad_lw_in_urb(:,:) = ieee_value(slurb_tile%rad_lw_in_urb,ieee_signaling_nan)
    slurb_tile%rad_sw_in_urb(:,:) = ieee_value(slurb_tile%rad_sw_in_urb,ieee_signaling_nan)
    slurb_tile%rad_lw_out_urb(:,:) = ieee_value(slurb_tile%rad_lw_out_urb,ieee_signaling_nan)
    slurb_tile%rad_sw_out_urb(:,:) = ieee_value(slurb_tile%rad_sw_out_urb,ieee_signaling_nan)

    slurb_tile%rad_lw_net_urb(:,:) = ieee_value(slurb_tile%rad_lw_net_urb,ieee_signaling_nan)
    slurb_tile%rad_sw_net_urb(:,:) = ieee_value(slurb_tile%rad_sw_net_urb,ieee_signaling_nan)

    slurb_tile%rad_lw_net_can(:,:) = ieee_value(slurb_tile%rad_lw_net_can,ieee_signaling_nan)

    slurb_tile%rad_lw_net_roof(:,:) = ieee_value(slurb_tile%rad_lw_net_roof,ieee_signaling_nan)
    slurb_tile%rad_sw_net_roof(:,:) = ieee_value(slurb_tile%rad_sw_net_roof,ieee_signaling_nan)
    slurb_tile%rad_lw_net_road(:,:) = ieee_value(slurb_tile%rad_lw_net_road,ieee_signaling_nan)
    slurb_tile%rad_sw_net_road(:,:) = ieee_value(slurb_tile%rad_sw_net_road,ieee_signaling_nan)
    slurb_tile%rad_sw_in_road(:,:) = ieee_value(slurb_tile%rad_sw_in_road,ieee_signaling_nan)
    slurb_tile%rad_lw_net_wall_a(:,:) = ieee_value(slurb_tile%rad_lw_net_wall_a,ieee_signaling_nan)
    slurb_tile%rad_sw_net_wall_a(:,:) = ieee_value(slurb_tile%rad_sw_net_wall_a,ieee_signaling_nan)
    slurb_tile%rad_lw_net_wall_b(:,:) = ieee_value(slurb_tile%rad_lw_net_wall_b,ieee_signaling_nan)
    slurb_tile%rad_sw_net_wall_b(:,:) = ieee_value(slurb_tile%rad_sw_net_wall_b,ieee_signaling_nan)
    slurb_tile%rad_lw_net_win_a(:,:) = ieee_value(slurb_tile%rad_lw_net_win_a,ieee_signaling_nan)
    slurb_tile%rad_sw_net_win_a(:,:) = ieee_value(slurb_tile%rad_sw_net_win_a,ieee_signaling_nan)
    slurb_tile%rad_sw_in_win_a(:,:) = ieee_value(slurb_tile%rad_sw_in_win_a,ieee_signaling_nan)
    slurb_tile%rad_lw_net_win_b(:,:) = ieee_value(slurb_tile%rad_lw_net_win_b,ieee_signaling_nan)
    slurb_tile%rad_sw_net_win_b(:,:) = ieee_value(slurb_tile%rad_sw_net_win_b,ieee_signaling_nan)
    slurb_tile%rad_sw_in_win_b(:,:) = ieee_value(slurb_tile%rad_sw_in_win_b,ieee_signaling_nan)

    slurb_tile%pt_can(:,:) = ieee_value(slurb_tile%pt_can,ieee_signaling_nan)
    slurb_tile%uv_abs_can(:,:) = ieee_value(slurb_tile%uv_abs_can,ieee_signaling_nan)
    slurb_tile%uv_eff_can(:,:) = ieee_value(slurb_tile%uv_eff_can,ieee_signaling_nan)
    slurb_tile%us_can(:,:) = ieee_value(slurb_tile%us_can,ieee_signaling_nan)
    slurb_tile%rib_can(:,:) = ieee_value(slurb_tile%rib_can,ieee_signaling_nan)
    slurb_tile%ol_can(:,:) = ieee_value(slurb_tile%ol_can,ieee_signaling_nan)

    slurb_tile%rib_roof(:,:) = ieee_value(slurb_tile%rib_roof,ieee_signaling_nan)
    slurb_tile%ol_roof(:,:) = ieee_value(slurb_tile%ol_roof,ieee_signaling_nan)
    slurb_tile%rib_road(:,:) = ieee_value(slurb_tile%rib_road,ieee_signaling_nan)
    slurb_tile%ol_road(:,:) = ieee_value(slurb_tile%ol_road,ieee_signaling_nan)

    slurb_tile%us_roof(:,:) = ieee_value(slurb_tile%us_roof,ieee_signaling_nan)
    slurb_tile%us_road(:,:) = ieee_value(slurb_tile%us_road,ieee_signaling_nan)

    slurb_tile%hw_can(:,:) = ieee_value(slurb_tile%hw_can,ieee_signaling_nan)
    ! slurb_tile%anisotropic_canyon(:,:) = ieee_value(slurb_tile%anisotropic_canyon,ieee_signaling_nan)
    slurb_tile%theta_can(:,:) = ieee_value(slurb_tile%theta_can,ieee_signaling_nan)
    slurb_tile%h_bld(:,:) = ieee_value(slurb_tile%h_bld,ieee_signaling_nan)
    slurb_tile%f_bld(:,:) = ieee_value(slurb_tile%f_bld,ieee_signaling_nan)
    slurb_tile%f_bld_frn(:,:) = ieee_value(slurb_tile%f_bld_frn,ieee_signaling_nan)
    slurb_tile%f_win(:,:) = ieee_value(slurb_tile%f_win,ieee_signaling_nan)
    slurb_tile%svf_road(:,:) = ieee_value(slurb_tile%svf_road,ieee_signaling_nan)
    slurb_tile%svf_wall(:,:) = ieee_value(slurb_tile%svf_wall,ieee_signaling_nan)
    slurb_tile%z0_urb(:,:) = ieee_value(slurb_tile%z0_urb,ieee_signaling_nan)

    slurb_tile%rah_roof(:,:) = ieee_value(slurb_tile%rah_roof,ieee_signaling_nan)
    slurb_tile%rah_road(:,:) = ieee_value(slurb_tile%rah_road,ieee_signaling_nan)
    slurb_tile%rah_can(:,:) = ieee_value(slurb_tile%rah_can,ieee_signaling_nan)

    IF ( facade_rah_doe )  THEN
       slurb_tile%rah_wall_a(:,:) = ieee_value(slurb_tile%rah_wall_a,ieee_signaling_nan)
       slurb_tile%rah_wall_b(:,:) = ieee_value(slurb_tile%rah_wall_b,ieee_signaling_nan)
       slurb_tile%rah_win_a(:,:) = ieee_value(slurb_tile%rah_win_a,ieee_signaling_nan)
       slurb_tile%rah_win_b(:,:) = ieee_value(slurb_tile%rah_win_b,ieee_signaling_nan)
    ELSE
       slurb_tile%rah_facade(:,:) = ieee_value(slurb_tile%rah_facade,ieee_signaling_nan)
    ENDIF

    slurb_tile%lambda_roof(:,:,:) = ieee_value(slurb_tile%lambda_roof,ieee_signaling_nan)
    slurb_tile%c_roof(:,:,:) = ieee_value(slurb_tile%c_roof,ieee_signaling_nan)
    slurb_tile%albedo_roof(:,:) = ieee_value(slurb_tile%albedo_roof,ieee_signaling_nan)
    slurb_tile%emiss_roof(:,:) = ieee_value(slurb_tile%emiss_roof,ieee_signaling_nan)
    slurb_tile%z0_roof(:,:) = ieee_value(slurb_tile%z0_roof,ieee_signaling_nan)
    slurb_tile%z0h_roof(:,:) = ieee_value(slurb_tile%z0h_roof,ieee_signaling_nan)
    slurb_tile%lambda_wall(:,:,:) = ieee_value(slurb_tile%lambda_wall,ieee_signaling_nan)
    slurb_tile%c_wall(:,:,:) = ieee_value(slurb_tile%c_wall,ieee_signaling_nan)
    slurb_tile%albedo_wall(:,:) = ieee_value(slurb_tile%albedo_wall,ieee_signaling_nan)
    slurb_tile%emiss_wall(:,:) = ieee_value(slurb_tile%emiss_wall,ieee_signaling_nan)
    slurb_tile%z0_wall(:,:) = ieee_value(slurb_tile%z0_wall,ieee_signaling_nan)
    slurb_tile%lambda_win(:,:,:) = ieee_value(slurb_tile%lambda_win,ieee_signaling_nan)
    slurb_tile%c_win(:,:,:) = ieee_value(slurb_tile%c_win,ieee_signaling_nan)
    slurb_tile%albedo_wall_win(:,:) = ieee_value(slurb_tile%albedo_wall_win,ieee_signaling_nan)
    slurb_tile%albedo_win(:,:) = ieee_value(slurb_tile%albedo_win,ieee_signaling_nan)
    slurb_tile%emiss_win(:,:) = ieee_value(slurb_tile%emiss_win,ieee_signaling_nan)
    slurb_tile%transmissivity_win(:,:) = ieee_value(slurb_tile%transmissivity_win,ieee_signaling_nan)
    slurb_tile%absorption_win(:,:,:) = ieee_value(slurb_tile%absorption_win,ieee_signaling_nan)
    slurb_tile%lambda_road(:,:,:) = ieee_value(slurb_tile%lambda_road,ieee_signaling_nan)
    slurb_tile%c_road(:,:,:) = ieee_value(slurb_tile%c_road,ieee_signaling_nan)
    slurb_tile%albedo_road(:,:) = ieee_value(slurb_tile%albedo_road,ieee_signaling_nan)
    slurb_tile%emiss_road(:,:) = ieee_value(slurb_tile%emiss_road,ieee_signaling_nan)
    slurb_tile%z0_road(:,:) = ieee_value(slurb_tile%z0_road,ieee_signaling_nan)
    slurb_tile%z0h_road(:,:) = ieee_value(slurb_tile%z0h_road,ieee_signaling_nan)

    slurb_tile%conductivity_roof(:,:,:) = ieee_value(slurb_tile%conductivity_roof,ieee_signaling_nan)
    slurb_tile%conductivity_wall(:,:,:) = ieee_value(slurb_tile%conductivity_wall,ieee_signaling_nan)
    slurb_tile%conductivity_win(:,:,:) = ieee_value(slurb_tile%conductivity_win,ieee_signaling_nan)
    slurb_tile%conductivity_road(:,:,:) = ieee_value(slurb_tile%conductivity_road,ieee_signaling_nan)

    slurb_tile%z_mo(:,:) = ieee_value(slurb_tile%z_mo,ieee_signaling_nan)
    slurb_tile%z_mo_can(:,:) = ieee_value(slurb_tile%z_mo_can,ieee_signaling_nan)
    slurb_tile%uv_abs_can_coef(:,:) = ieee_value(slurb_tile%uv_abs_can_coef,ieee_signaling_nan)
    slurb_tile%wall_hor_a_ratio(:,:) = ieee_value(slurb_tile%wall_hor_a_ratio,ieee_signaling_nan)


    slurb_tile%lw_roof_coef(:,:,:) = ieee_value(slurb_tile%lw_roof_coef,ieee_signaling_nan)
    slurb_tile%lw_road_coef(:,:,:) = ieee_value(slurb_tile%lw_road_coef,ieee_signaling_nan)
    slurb_tile%lw_wall_coef(:,:,:) = ieee_value(slurb_tile%lw_wall_coef,ieee_signaling_nan)
    slurb_tile%lw_win_coef(:,:,:) = ieee_value(slurb_tile%lw_win_coef,ieee_signaling_nan)
    slurb_tile%sw_ref_denom(:,:) = ieee_value(slurb_tile%sw_ref_denom,ieee_signaling_nan)

    slurb_tile%us_urb(:,:) = ieee_value(slurb_tile%us_urb,ieee_signaling_nan)
    slurb_tile%uv_eff1(:,:) = ieee_value(slurb_tile%uv_eff1,ieee_signaling_nan)
    slurb_tile%uv_abs1(:,:) = ieee_value(slurb_tile%uv_abs1,ieee_signaling_nan)
    slurb_tile%pt1(:,:) = ieee_value(slurb_tile%pt1,ieee_signaling_nan)

    slurb_tile%t_can_0(:,:) = ieee_value(slurb_tile%t_can_0,ieee_signaling_nan)
    slurb_tile%t_can_m(:,:) = ieee_value(slurb_tile%t_can_m,ieee_signaling_nan)
    slurb_tile%t_wall_a_0(:,:,:) = ieee_value(slurb_tile%t_wall_a_0,ieee_signaling_nan)
    slurb_tile%t_wall_a_m(:,:,:) = ieee_value(slurb_tile%t_wall_a_m,ieee_signaling_nan)
    slurb_tile%t_wall_b_0(:,:,:) = ieee_value(slurb_tile%t_wall_b_0,ieee_signaling_nan)
    slurb_tile%t_wall_b_m(:,:,:) = ieee_value(slurb_tile%t_wall_b_m,ieee_signaling_nan)
    slurb_tile%t_win_a_0(:,:,:) = ieee_value(slurb_tile%t_win_a_0,ieee_signaling_nan)
    slurb_tile%t_win_a_m(:,:,:) = ieee_value(slurb_tile%t_win_a_m,ieee_signaling_nan)
    slurb_tile%t_win_b_0(:,:,:) = ieee_value(slurb_tile%t_win_b_0,ieee_signaling_nan)
    slurb_tile%t_win_b_m(:,:,:) = ieee_value(slurb_tile%t_win_b_m,ieee_signaling_nan)
    slurb_tile%t_roof_0(:,:,:) = ieee_value(slurb_tile%t_roof_0,ieee_signaling_nan)
    slurb_tile%t_roof_m(:,:,:) = ieee_value(slurb_tile%t_roof_m,ieee_signaling_nan)
    slurb_tile%t_road_0(:,:,:) = ieee_value(slurb_tile%t_road_0,ieee_signaling_nan)
    slurb_tile%t_road_m(:,:,:) = ieee_value(slurb_tile%t_road_m,ieee_signaling_nan)

    ! ! initialise pointers
    ! slurb_tile%q_can_0 => slurb_tile%q_can_m; slurb_tile%q_can_m => slurb_tile%q_can_0
    ! slurb_tile%t_can_0 => slurb_tile%t_can_m; slurb_tile%t_can_m => slurb_tile%t_can_0
    ! slurb_tile%m_liq_road_m => slurb_tile%m_liq_road_m; slurb_tile%m_liq_road_0 => slurb_tile%m_liq_road_0
    ! slurb_tile%m_liq_roof_m => slurb_tile%m_liq_roof_m; slurb_tile%m_liq_roof_0 => slurb_tile%m_liq_roof_0
    ! slurb_tile%t_wall_a_0 => slurb_tile%t_wall_a_m; slurb_tile%t_wall_a_m => slurb_tile%t_wall_a_0
    ! slurb_tile%t_wall_b_0 => slurb_tile%t_wall_b_m; slurb_tile%t_wall_b_m => slurb_tile%t_wall_b_0
    ! slurb_tile%t_win_a_0 => slurb_tile%t_win_a_m; slurb_tile%t_win_a_m => slurb_tile%t_win_a_0
    ! slurb_tile%t_win_b_0 => slurb_tile%t_win_b_m; slurb_tile%t_win_b_m => slurb_tile%t_win_b_0
    ! slurb_tile%t_roof_0 => slurb_tile%t_roof_m; slurb_tile%t_roof_m => slurb_tile%t_roof_0
    ! slurb_tile%t_road_0 => slurb_tile%t_road_m; slurb_tile%t_road_m => slurb_tile%t_road_0

    IF ( moist_physics )  THEN
       slurb_tile%tq_can(:,:) = ieee_value(slurb_tile%tq_can,ieee_signaling_nan)
       slurb_tile%tm_liq_roof(:,:) = ieee_value(slurb_tile%tm_liq_roof,ieee_signaling_nan)
       slurb_tile%tm_liq_road(:,:) = ieee_value(slurb_tile%tm_liq_road,ieee_signaling_nan)
       slurb_tile%tm_roof_runoff(:,:) = ieee_value(slurb_tile%tm_roof_runoff,ieee_signaling_nan)
       slurb_tile%tm_road_runoff(:,:) = ieee_value(slurb_tile%tm_road_runoff,ieee_signaling_nan)
       slurb_tile%tm_roof_precep(:,:) = ieee_value(slurb_tile%tm_roof_precep,ieee_signaling_nan)
       slurb_tile%tm_road_precep(:,:) = ieee_value(slurb_tile%tm_road_precep,ieee_signaling_nan)

       slurb_tile%vpt_roof(:,:) = ieee_value(slurb_tile%vpt_roof,ieee_signaling_nan)
       slurb_tile%vpt_road(:,:) = ieee_value(slurb_tile%vpt_road,ieee_signaling_nan)

       slurb_tile%q_roof(:,:) = ieee_value(slurb_tile%q_roof,ieee_signaling_nan)
       slurb_tile%q_road(:,:) = ieee_value(slurb_tile%q_road,ieee_signaling_nan)
       slurb_tile%qs_roof(:,:) = ieee_value(slurb_tile%qs_roof,ieee_signaling_nan)
       slurb_tile%qs_road(:,:) = ieee_value(slurb_tile%qs_road,ieee_signaling_nan)

       slurb_tile%qsws_can(:,:) = ieee_value(slurb_tile%qsws_can,ieee_signaling_nan)
       slurb_tile%qsws_roof(:,:) = ieee_value(slurb_tile%qsws_roof,ieee_signaling_nan)
       slurb_tile%qsws_road(:,:) = ieee_value(slurb_tile%qsws_road,ieee_signaling_nan)
       slurb_tile%qsws_liq_roof(:,:) = ieee_value(slurb_tile%qsws_liq_roof,ieee_signaling_nan)
       slurb_tile%qsws_liq_road(:,:) = ieee_value(slurb_tile%qsws_liq_road,ieee_signaling_nan)

       slurb_tile%c_liq_roof(:,:) = ieee_value(slurb_tile%c_liq_roof,ieee_signaling_nan)
       slurb_tile%c_liq_road(:,:) = ieee_value(slurb_tile%c_liq_road,ieee_signaling_nan)

       slurb_tile%vpt_can(:,:) = ieee_value(slurb_tile%vpt_can,ieee_signaling_nan)

       slurb_tile%q1(:,:) = ieee_value(slurb_tile%q1,ieee_signaling_nan)
       slurb_tile%vpt1(:,:) = ieee_value(slurb_tile%vpt1,ieee_signaling_nan)

       slurb_tile%qsws_external(:,:) = ieee_value(slurb_tile%qsws_external,ieee_signaling_nan)

       slurb_tile%q_can_0(:,:) = ieee_value(slurb_tile%q_can_0,ieee_signaling_nan)
       slurb_tile%q_can_m(:,:) = ieee_value(slurb_tile%q_can_m,ieee_signaling_nan)
       slurb_tile%m_liq_roof_0(:,:) = ieee_value(slurb_tile%m_liq_roof_0,ieee_signaling_nan)
       slurb_tile%m_liq_roof_m(:,:) = ieee_value(slurb_tile%m_liq_roof_m,ieee_signaling_nan)
       slurb_tile%m_liq_road_0(:,:) = ieee_value(slurb_tile%m_liq_road_0,ieee_signaling_nan)
       slurb_tile%m_liq_road_m(:,:) = ieee_value(slurb_tile%m_liq_road_m,ieee_signaling_nan)

    ENDIF

    slurb_tile%dt_max(:,:) = ieee_value(slurb_tile%dt_max,ieee_signaling_nan)
#endif
end subroutine slurb_bulk_allocations

subroutine slurb_bulk_deallocations
    DEALLOCATE(fraction_slurb)

    DEALLOCATE(ln_z_z0_roof)
    DEALLOCATE(ln_z_z0h_roof)
    DEALLOCATE(ln_z_z0_urb)
    DEALLOCATE(pt_surface)
    DEALLOCATE(ln_z_z0_road)
    DEALLOCATE(ln_z_z0h_road)
    !-- Bulk allocation
    DEALLOCATE( slurb_tile%dz_roof)
    DEALLOCATE( slurb_tile%dz_wall)
    DEALLOCATE( slurb_tile%dz_road)
    DEALLOCATE( slurb_tile%dz_win)
    DEALLOCATE( slurb_tile%zw_win)

    DEALLOCATE( slurb_tile%t_c_urb)
    DEALLOCATE( slurb_tile%thl_rad_urb)
    DEALLOCATE( slurb_tile%t_h_urb)
    DEALLOCATE( slurb_tile%t_2m_urb)
    DEALLOCATE( slurb_tile%shf_urb)
    DEALLOCATE( slurb_tile%qsws_urb)
    DEALLOCATE( slurb_tile%ol_urb)
    DEALLOCATE( slurb_tile%rib_urb)
    DEALLOCATE( slurb_tile%ram_urb)
    DEALLOCATE( slurb_tile%usws_urb)
    DEALLOCATE( slurb_tile%vsws_urb)
    DEALLOCATE( slurb_tile%thlskin)
    DEALLOCATE( slurb_tile%qtskin)

    DEALLOCATE( slurb_tile%albedo_urb)
    DEALLOCATE( slurb_tile%emiss_urb)

    DEALLOCATE( slurb_tile%t_indoor)
    DEALLOCATE( slurb_tile%t_soil)

    DEALLOCATE( slurb_tile%tt_can)
    DEALLOCATE( slurb_tile%tt_wall_a)
    DEALLOCATE( slurb_tile%tt_wall_b)
    DEALLOCATE( slurb_tile%tt_win_a)
    DEALLOCATE( slurb_tile%tt_win_b)
    DEALLOCATE( slurb_tile%tt_roof)
    DEALLOCATE( slurb_tile%tt_road)

    DEALLOCATE( slurb_tile%pt_wall_a)
    DEALLOCATE( slurb_tile%pt_wall_b)
    DEALLOCATE( slurb_tile%pt_win_a)
    DEALLOCATE( slurb_tile%pt_win_b)
    DEALLOCATE( slurb_tile%pt_roof)
    DEALLOCATE( slurb_tile%pt_road)

    DEALLOCATE( slurb_tile%shf_can)
    DEALLOCATE( slurb_tile%shf_roof)
    DEALLOCATE( slurb_tile%shf_road)
    DEALLOCATE( slurb_tile%shf_wall_a)
    DEALLOCATE( slurb_tile%shf_wall_b)
    DEALLOCATE( slurb_tile%shf_win_a)
    DEALLOCATE( slurb_tile%shf_win_b)

    DEALLOCATE( slurb_tile%shf_external)
    DEALLOCATE( slurb_tile%shf_traffic)

    DEALLOCATE( slurb_tile%ghf_road)
    DEALLOCATE( slurb_tile%ghf_roof)
    DEALLOCATE( slurb_tile%ghf_wall_a)
    DEALLOCATE( slurb_tile%ghf_wall_b)
    DEALLOCATE( slurb_tile%ghf_win_a)
    DEALLOCATE( slurb_tile%ghf_win_b)

    DEALLOCATE( slurb_tile%rad_lw_in_urb)
    DEALLOCATE( slurb_tile%rad_sw_in_urb)
    DEALLOCATE( slurb_tile%rad_lw_out_urb)
    DEALLOCATE( slurb_tile%rad_sw_out_urb)

    DEALLOCATE( slurb_tile%rad_lw_net_urb)
    DEALLOCATE( slurb_tile%rad_sw_net_urb)

    DEALLOCATE( slurb_tile%rad_lw_net_can)

    DEALLOCATE( slurb_tile%rad_lw_net_roof)
    DEALLOCATE( slurb_tile%rad_sw_net_roof)
    DEALLOCATE( slurb_tile%rad_lw_net_road)
    DEALLOCATE( slurb_tile%rad_sw_net_road)
    DEALLOCATE( slurb_tile%rad_sw_in_road)
    DEALLOCATE( slurb_tile%rad_lw_net_wall_a)
    DEALLOCATE( slurb_tile%rad_sw_net_wall_a)
    DEALLOCATE( slurb_tile%rad_lw_net_wall_b)
    DEALLOCATE( slurb_tile%rad_sw_net_wall_b)
    DEALLOCATE( slurb_tile%rad_lw_net_win_a)
    DEALLOCATE( slurb_tile%rad_sw_net_win_a)
    DEALLOCATE( slurb_tile%rad_sw_in_win_a)
    DEALLOCATE( slurb_tile%rad_lw_net_win_b)
    DEALLOCATE( slurb_tile%rad_sw_net_win_b)
    DEALLOCATE( slurb_tile%rad_sw_in_win_b)

    DEALLOCATE( slurb_tile%pt_can)
    DEALLOCATE( slurb_tile%uv_abs_can)
    DEALLOCATE( slurb_tile%uv_eff_can)
    DEALLOCATE( slurb_tile%us_can)
    DEALLOCATE( slurb_tile%rib_can)
    DEALLOCATE( slurb_tile%ol_can)

    DEALLOCATE( slurb_tile%rib_roof)
    DEALLOCATE( slurb_tile%ol_roof)
    DEALLOCATE( slurb_tile%rib_road)
    DEALLOCATE( slurb_tile%ol_road)

    DEALLOCATE( slurb_tile%us_roof)
    DEALLOCATE( slurb_tile%us_road)

    DEALLOCATE( slurb_tile%hw_can)
    DEALLOCATE( slurb_tile%anisotropic_canyon)
    DEALLOCATE( slurb_tile%theta_can)
    DEALLOCATE( slurb_tile%h_bld)
    DEALLOCATE( slurb_tile%f_bld)
    DEALLOCATE( slurb_tile%f_bld_frn)
    DEALLOCATE( slurb_tile%f_win)
    DEALLOCATE( slurb_tile%svf_road)
    DEALLOCATE( slurb_tile%svf_wall)
    DEALLOCATE( slurb_tile%z0_urb)

    DEALLOCATE( slurb_tile%rah_roof)
    DEALLOCATE( slurb_tile%rah_road)
    DEALLOCATE( slurb_tile%rah_can)

    IF ( facade_rah_doe )  THEN
       DEALLOCATE( slurb_tile%rah_wall_a)
       DEALLOCATE( slurb_tile%rah_wall_b)
       DEALLOCATE( slurb_tile%rah_win_a)
       DEALLOCATE( slurb_tile%rah_win_b)
    ELSE
       DEALLOCATE( slurb_tile%rah_facade)
    ENDIF

    DEALLOCATE( slurb_tile%lambda_roof)
    DEALLOCATE( slurb_tile%c_roof)
    DEALLOCATE( slurb_tile%albedo_roof)
    DEALLOCATE( slurb_tile%emiss_roof)
    DEALLOCATE( slurb_tile%z0_roof)
    DEALLOCATE( slurb_tile%z0h_roof)
    DEALLOCATE( slurb_tile%lambda_wall)
    DEALLOCATE( slurb_tile%c_wall)
    DEALLOCATE( slurb_tile%albedo_wall)
    DEALLOCATE( slurb_tile%emiss_wall)
    DEALLOCATE( slurb_tile%z0_wall)
    DEALLOCATE( slurb_tile%lambda_win)
    DEALLOCATE( slurb_tile%c_win)
    DEALLOCATE( slurb_tile%albedo_wall_win)
    DEALLOCATE( slurb_tile%albedo_win)
    DEALLOCATE( slurb_tile%emiss_win)
    DEALLOCATE( slurb_tile%transmissivity_win)
    DEALLOCATE( slurb_tile%absorption_win)
    DEALLOCATE( slurb_tile%lambda_road)
    DEALLOCATE( slurb_tile%c_road)
    DEALLOCATE( slurb_tile%albedo_road)
    DEALLOCATE( slurb_tile%emiss_road)
    DEALLOCATE( slurb_tile%z0_road)
    DEALLOCATE( slurb_tile%z0h_road)

    DEALLOCATE( slurb_tile%conductivity_roof)
    DEALLOCATE( slurb_tile%conductivity_wall)
    DEALLOCATE( slurb_tile%conductivity_win)
    DEALLOCATE( slurb_tile%conductivity_road)

    DEALLOCATE( slurb_tile%z_mo)
    DEALLOCATE( slurb_tile%z_mo_can)
    DEALLOCATE( slurb_tile%uv_abs_can_coef)
    DEALLOCATE( slurb_tile%wall_hor_a_ratio)


    DEALLOCATE( slurb_tile%lw_roof_coef)
    DEALLOCATE( slurb_tile%lw_road_coef)
    DEALLOCATE( slurb_tile%lw_wall_coef)
    DEALLOCATE( slurb_tile%lw_win_coef)
    DEALLOCATE( slurb_tile%sw_ref_denom)

    DEALLOCATE( slurb_tile%us_urb)
    DEALLOCATE( slurb_tile%uv_eff1)
    DEALLOCATE( slurb_tile%uv_abs1)
    DEALLOCATE( slurb_tile%pt1)

    DEALLOCATE( slurb_tile%t_can_0)
    DEALLOCATE( slurb_tile%t_can_m)
    DEALLOCATE( slurb_tile%t_wall_a_0)
    DEALLOCATE( slurb_tile%t_wall_a_m)
    DEALLOCATE( slurb_tile%t_wall_b_0)
    DEALLOCATE( slurb_tile%t_wall_b_m)
    DEALLOCATE( slurb_tile%t_win_a_0)
    DEALLOCATE( slurb_tile%t_win_a_m)
    DEALLOCATE( slurb_tile%t_win_b_0)
    DEALLOCATE( slurb_tile%t_win_b_m)
    DEALLOCATE( slurb_tile%t_roof_0)
    DEALLOCATE( slurb_tile%t_roof_m)
    DEALLOCATE( slurb_tile%t_road_0)
    DEALLOCATE( slurb_tile%t_road_m)

    IF ( moist_physics )  THEN
       DEALLOCATE( slurb_tile%tq_can)
       DEALLOCATE( slurb_tile%tm_liq_roof)
       DEALLOCATE( slurb_tile%tm_liq_road)
       DEALLOCATE( slurb_tile%tm_roof_runoff)
       DEALLOCATE( slurb_tile%tm_road_runoff)
       DEALLOCATE( slurb_tile%tm_roof_precep)
       DEALLOCATE( slurb_tile%tm_road_precep)
       DEALLOCATE( slurb_tile%vpt_roof)
       DEALLOCATE( slurb_tile%vpt_road)

       DEALLOCATE( slurb_tile%q_roof)
       DEALLOCATE( slurb_tile%q_road)
       DEALLOCATE( slurb_tile%qs_roof)
       DEALLOCATE( slurb_tile%qs_road)

       DEALLOCATE( slurb_tile%qsws_can)
       DEALLOCATE( slurb_tile%qsws_roof)
       DEALLOCATE( slurb_tile%qsws_road)
       DEALLOCATE( slurb_tile%qsws_liq_roof)
       DEALLOCATE( slurb_tile%qsws_liq_road)

       DEALLOCATE( slurb_tile%c_liq_roof)
       DEALLOCATE( slurb_tile%c_liq_road)

       DEALLOCATE( slurb_tile%vpt_can)

       DEALLOCATE( slurb_tile%q1)
       DEALLOCATE( slurb_tile%vpt1)

       DEALLOCATE( slurb_tile%qsws_external)

       DEALLOCATE( slurb_tile%q_can_0)
       DEALLOCATE( slurb_tile%q_can_m)
       DEALLOCATE( slurb_tile%m_liq_roof_0)
       DEALLOCATE( slurb_tile%m_liq_roof_m)
       DEALLOCATE( slurb_tile%m_liq_road_0)
       DEALLOCATE( slurb_tile%m_liq_road_m)
    ENDIF
end subroutine slurb_bulk_deallocations


    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    ! Swap timelevel of the SLUrb model.
    !--------------------------------------------------------------------------------------------------!
 SUBROUTINE slurb_set_previous_timestep()
    use modglobal, only: rk3step

    if (rk3step == 1) then
        slurb_tile%t_wall_a_m(:,:,:) = slurb_tile%t_wall_a_0(:,:,:)
        slurb_tile%t_wall_b_m(:,:,:) = slurb_tile%t_wall_b_0(:,:,:)
        slurb_tile%t_win_a_m(:,:,:) = slurb_tile%t_win_a_0(:,:,:)
        slurb_tile%t_win_b_m(:,:,:) = slurb_tile%t_win_b_0(:,:,:)
        slurb_tile%t_roof_m(:,:,:) = slurb_tile%t_roof_0(:,:,:)
        slurb_tile%t_road_m(:,:,:) = slurb_tile%t_road_0(:,:,:)
        slurb_tile%t_can_m(:,:) = slurb_tile%t_can_0(:,:)
        slurb_tile%q_can_m(:,:) = slurb_tile%q_can_0(:,:)
        slurb_tile%m_liq_roof_m(:,:) = slurb_tile%m_liq_roof_0(:,:)
        slurb_tile%m_liq_road_m(:,:) = slurb_tile%m_liq_road_0(:,:)
    endif

 END SUBROUTINE slurb_set_previous_timestep

  !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> This function computes the magnus formula (Press et al., 1992).
    !> The magnus formula is needed to calculate the saturation vapor pressure.
    !--------------------------------------------------------------------------------------------------!
 FUNCTION magnus( t )
    !$ACC ROUTINE SEQ

    IMPLICIT NONE

    REAL(field_r), INTENT(IN) ::  t  !< temperature (K)

    REAL(field_r) ::  magnus

    !
    !-- Saturation vapor pressure for a specific temperature:
    magnus =  611.2_field_r * EXP( 17.62_field_r * ( t - 273.15_field_r ) / ( t - 29.65_field_r  ) )

 END FUNCTION magnus

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Computes a steady-state solution for 1D heat equation using Gauss-Seidel iteration.
    !> This is used to initialize the material temperatures for roofs, walls, windows and roads,
    !> shortening the time required for the spinup. For windows, SW absorption is not considered.
    !--------------------------------------------------------------------------------------------------!
 !TODOSELF PURE
 FUNCTION calc_1d_heat_equation( result_size , t_bc_1, t_bc_2, lambda ) RESULT( t_result )

    INTEGER, INTENT(IN) ::  result_size  !< output target size

    REAL(field_r), INTENT(IN) ::  t_bc_1  !< outer t boundary condition
    REAL(field_r), INTENT(IN) ::  t_bc_2  !< inner t boundary condition

    REAL(field_r), DIMENSION(:), INTENT(IN) ::  lambda  !< total layer heat conductivity

    INTEGER ::  ix  !< iteration counter
    INTEGER ::  kx  !< layer running index

    REAL(field_r), PARAMETER ::  omega = 1.0_field_r   !< relaxation to control convergence
    REAL(field_r), PARAMETER ::  tol = 1.0E-6_field_r  !< maximum residual for convergence

    REAL(field_r) ::  res    !< iteration residual for convergence check
    REAL(field_r) ::  t_old  !< previous t of layer for convergence check

    REAL(field_r), DIMENSION(result_size) ::  t_result  !< result t profile
    REAL(field_r), DIMENSION(1:result_size+1) ::  t     !< intermediate t array containing also the BCs


    !
    !-- Set boundary conditions for the iteration array. The layer against the atmosphere will have a
    !-- constant boundary condition and the other boundary is treated similarly to the inner layer
    !-- boundary condition in prognostic equations. Thus, an extra layer is neeeded for the inner
    !-- temperature array for iteration.
    t(LBOUND( t, 1 )) = t_bc_1
    t(UBOUND( t, 1 )) = t_bc_2

    !
    !-- Set initial guess for temperature for subsurface layers.
    t(LBOUND( t, 1 )+1:UBOUND( t, 1 )-1) = (t_bc_1 + t_bc_2) / 2.0_field_r
 
    !
    !-- Gauss-Seidel iteration.
    DO  ix = 1, 1000
       DO  kx = LBOUND( t, 1 )+1, UBOUND( t, 1 )-1
          t_old = t(kx)
          res = 0.0_field_r
          t(kx) = ( lambda(kx) * ( t(kx+1) - t(kx) ) + lambda(kx-1) * ( t(kx-1) - t(kx) ) ) /      &
                  ( lambda(kx) + lambda(kx-1) ) * omega + t(kx)
          res = MAX( res, ABS( t(kx) - t_old ) )
       ENDDO
    !
    !--    Check for convergence using the stored maximum residual.
       IF ( res < tol )  EXIT
    ENDDO
    !
    !-- Return the solution.
    t_result(:) = t(LBOUND( t, 1 ):UBOUND( t, 1 )-1)

 END FUNCTION calc_1d_heat_equation

end module modslurbhelpers
