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
    allocate( slurb_tile%dz_roof(nzt_roof:nzb_roof,i2,j2) )
    allocate( slurb_tile%dz_wall(nzt_wall:nzb_wall,i2,j2) )
    allocate( slurb_tile%dz_road(nzt_road:nzb_road,i2,j2) )
    allocate( slurb_tile%dz_win(nzt_win:nzb_win,i2,j2) )
    allocate( slurb_tile%zw_win(nzt_win:nzb_win,i2,j2) )

    allocate( slurb_tile%t_c_urb(i2,j2) )
    allocate( slurb_tile%thl_rad_urb(i2,j2) )
    allocate( slurb_tile%t_h_urb(i2,j2) )
    allocate( slurb_tile%t_2m_urb(i2,j2) )
    allocate( slurb_tile%shf_urb(i2,j2) )
    allocate( slurb_tile%qsws_urb(i2,j2) )
    allocate( slurb_tile%ol_urb(i2,j2) )
    allocate( slurb_tile%rib_urb(i2,j2) )
    allocate( slurb_tile%ram_urb(i2,j2) )
    allocate( slurb_tile%usws_urb(i2,j2) )
    allocate( slurb_tile%vsws_urb(i2,j2) )
    allocate( slurb_tile%thlskin(i2,j2) )
    allocate( slurb_tile%qtskin(i2,j2) )

    allocate( slurb_tile%albedo_urb(i2,j2) )
    allocate( slurb_tile%emiss_urb(i2,j2) )

    allocate( slurb_tile%t_indoor(i2,j2) )
    allocate( slurb_tile%t_soil(i2,j2) )

    allocate( slurb_tile%tt_can(i2,j2) )
    allocate( slurb_tile%tt_wall_a(nzt_wall:nzb_wall,i2,j2) )
    allocate( slurb_tile%tt_wall_b(nzt_wall:nzb_wall,i2,j2) )
    allocate( slurb_tile%tt_win_a(nzt_win:nzb_win,i2,j2) )
    allocate( slurb_tile%tt_win_b(nzt_win:nzb_win,i2,j2) )
    allocate( slurb_tile%tt_roof(nzt_roof:nzb_roof,i2,j2) )
    allocate( slurb_tile%tt_road(nzt_road:nzb_road,i2,j2) )

    allocate( slurb_tile%pt_wall_a(i2,j2) )
    allocate( slurb_tile%pt_wall_b(i2,j2) )
    allocate( slurb_tile%pt_win_a(i2,j2) )
    allocate( slurb_tile%pt_win_b(i2,j2) )
    allocate( slurb_tile%pt_roof(i2,j2) )
    allocate( slurb_tile%pt_road(i2,j2) )

    allocate( slurb_tile%shf_can(i2,j2) )
    allocate( slurb_tile%shf_roof(i2,j2) )
    allocate( slurb_tile%shf_road(i2,j2) )
    allocate( slurb_tile%shf_wall_a(i2,j2) )
    allocate( slurb_tile%shf_wall_b(i2,j2) )
    allocate( slurb_tile%shf_win_a(i2,j2) )
    allocate( slurb_tile%shf_win_b(i2,j2) )

    allocate( slurb_tile%shf_external(i2,j2) )
    allocate( slurb_tile%shf_traffic(i2,j2) )

    allocate( slurb_tile%ghf_road(i2,j2) )
    allocate( slurb_tile%ghf_roof(i2,j2) )
    allocate( slurb_tile%ghf_wall_a(i2,j2) )
    allocate( slurb_tile%ghf_wall_b(i2,j2) )
    allocate( slurb_tile%ghf_win_a(i2,j2) )
    allocate( slurb_tile%ghf_win_b(i2,j2) )

    allocate( slurb_tile%rad_lw_in_urb(i2,j2) )
    allocate( slurb_tile%rad_sw_in_urb(i2,j2) )
    allocate( slurb_tile%rad_lw_out_urb(i2,j2) )
    allocate( slurb_tile%rad_sw_out_urb(i2,j2) )

    allocate( slurb_tile%rad_lw_net_urb(i2,j2) )
    allocate( slurb_tile%rad_sw_net_urb(i2,j2) )

    allocate( slurb_tile%rad_lw_net_can(i2,j2) )

    allocate( slurb_tile%rad_lw_net_roof(i2,j2) )
    allocate( slurb_tile%rad_sw_net_roof(i2,j2) )
    allocate( slurb_tile%rad_lw_net_road(i2,j2) )
    allocate( slurb_tile%rad_sw_net_road(i2,j2) )
    allocate( slurb_tile%rad_sw_in_road(i2,j2) )
    allocate( slurb_tile%rad_lw_net_wall_a(i2,j2) )
    allocate( slurb_tile%rad_sw_net_wall_a(i2,j2) )
    allocate( slurb_tile%rad_lw_net_wall_b(i2,j2) )
    allocate( slurb_tile%rad_sw_net_wall_b(i2,j2) )
    allocate( slurb_tile%rad_lw_net_win_a(i2,j2) )
    allocate( slurb_tile%rad_sw_net_win_a(i2,j2) )
    allocate( slurb_tile%rad_sw_in_win_a(i2,j2) )
    allocate( slurb_tile%rad_lw_net_win_b(i2,j2) )
    allocate( slurb_tile%rad_sw_net_win_b(i2,j2) )
    allocate( slurb_tile%rad_sw_in_win_b(i2,j2) )

    allocate( slurb_tile%pt_can(i2,j2) )
    allocate( slurb_tile%uv_abs_can(i2,j2) )
    allocate( slurb_tile%uv_eff_can(i2,j2) )
    allocate( slurb_tile%us_can(i2,j2) )
    allocate( slurb_tile%rib_can(i2,j2) )
    allocate( slurb_tile%ol_can(i2,j2) )

    allocate( slurb_tile%rib_roof(i2,j2) )
    allocate( slurb_tile%ol_roof(i2,j2) )
    allocate( slurb_tile%rib_road(i2,j2) )
    allocate( slurb_tile%ol_road(i2,j2) )

    allocate( slurb_tile%us_roof(i2,j2) )
    allocate( slurb_tile%us_road(i2,j2) )

    allocate( slurb_tile%hw_can(i2,j2) )
    allocate( slurb_tile%anisotropic_canyon(i2,j2) )
    allocate( slurb_tile%theta_can(i2,j2) )
    allocate( slurb_tile%h_bld(i2,j2) )
    allocate( slurb_tile%f_bld(i2,j2) )
    allocate( slurb_tile%f_bld_frn(i2,j2) )
    allocate( slurb_tile%f_win(i2,j2) )
    allocate( slurb_tile%svf_road(i2,j2) )
    allocate( slurb_tile%svf_wall(i2,j2) )
    allocate( slurb_tile%z0_urb(i2,j2) )

    allocate( slurb_tile%rah_roof(i2,j2) )
    allocate( slurb_tile%rah_road(i2,j2) )
    allocate( slurb_tile%rah_can(i2,j2) )

    if ( facade_rah_doe )  then
       allocate( slurb_tile%rah_wall_a(i2,j2) )
       allocate( slurb_tile%rah_wall_b(i2,j2) )
       allocate( slurb_tile%rah_win_a(i2,j2) )
       allocate( slurb_tile%rah_win_b(i2,j2) )
    else
       allocate( slurb_tile%rah_facade(i2,j2) )
    endif

    allocate( slurb_tile%lambda_roof(nzt_roof:nzb_roof,i2,j2) )
    allocate( slurb_tile%c_roof(nzt_roof:nzb_roof,i2,j2) )
    allocate( slurb_tile%albedo_roof(i2,j2) )
    allocate( slurb_tile%emiss_roof(i2,j2) )
    allocate( slurb_tile%z0_roof(i2,j2) )
    allocate( slurb_tile%z0h_roof(i2,j2) )
    allocate( slurb_tile%lambda_wall(nzt_wall:nzb_wall,i2,j2) )
    allocate( slurb_tile%c_wall(nzt_wall:nzb_wall,i2,j2) )
    allocate( slurb_tile%albedo_wall(i2,j2) )
    allocate( slurb_tile%emiss_wall(i2,j2) )
    allocate( slurb_tile%z0_wall(i2,j2) )
    allocate( slurb_tile%lambda_win(nzt_win:nzb_win,i2,j2) )
    allocate( slurb_tile%c_win(nzt_win:nzb_win,i2,j2) )
    allocate( slurb_tile%albedo_wall_win(i2,j2) )
    allocate( slurb_tile%albedo_win(i2,j2) )
    allocate( slurb_tile%emiss_win(i2,j2) )
    allocate( slurb_tile%transmissivity_win(i2,j2) )
    allocate( slurb_tile%absorption_win(nzt_win:nzb_win,i2,j2) )
    allocate( slurb_tile%lambda_road(nzt_road:nzb_road,i2,j2) )
    allocate( slurb_tile%c_road(nzt_road:nzb_road,i2,j2) )
    allocate( slurb_tile%albedo_road(i2,j2) )
    allocate( slurb_tile%emiss_road(i2,j2) )
    allocate( slurb_tile%z0_road(i2,j2) )
    allocate( slurb_tile%z0h_road(i2,j2) )

    allocate( slurb_tile%conductivity_roof(nzt_roof:nzb_roof,i2,j2) )
    allocate( slurb_tile%conductivity_wall(nzt_wall:nzb_wall,i2,j2) )
    allocate( slurb_tile%conductivity_win(nzt_win:nzb_win,i2,j2) )
    allocate( slurb_tile%conductivity_road(nzt_road:nzb_road,i2,j2) )

    allocate( slurb_tile%z_mo(i2,j2) )
    allocate( slurb_tile%z_mo_can(i2,j2) )
    allocate( slurb_tile%uv_abs_can_coef(i2,j2) )
    allocate( slurb_tile%wall_hor_a_ratio(i2,j2) )


    allocate( slurb_tile%lw_roof_coef(1:2,i2,j2) )
    allocate( slurb_tile%lw_road_coef(1:4,i2,j2) )
    allocate( slurb_tile%lw_wall_coef(1:6,i2,j2) )
    allocate( slurb_tile%lw_win_coef(1:6,i2,j2) )
    allocate( slurb_tile%sw_ref_denom(i2,j2) )

    allocate( slurb_tile%us_urb(i2,j2) )
    allocate( slurb_tile%uv_eff1(i2,j2) )
    allocate( slurb_tile%uv_abs1(i2,j2) )
    allocate( slurb_tile%pt1(i2,j2) )

    allocate( slurb_tile%t_can_0(i2,j2) )
    allocate( slurb_tile%t_can_m(i2,j2) )
    allocate( slurb_tile%t_wall_a_0(nzt_wall:nzb_wall,i2,j2) )
    allocate( slurb_tile%t_wall_a_m(nzt_wall:nzb_wall,i2,j2) )
    allocate( slurb_tile%t_wall_b_0(nzt_wall:nzb_wall,i2,j2) )
    allocate( slurb_tile%t_wall_b_m(nzt_wall:nzb_wall,i2,j2) )
    allocate( slurb_tile%t_win_a_0(nzt_win:nzb_win,i2,j2) )
    allocate( slurb_tile%t_win_a_m(nzt_win:nzb_win,i2,j2) )
    allocate( slurb_tile%t_win_b_0(nzt_win:nzb_win,i2,j2) )
    allocate( slurb_tile%t_win_b_m(nzt_win:nzb_win,i2,j2) )
    allocate( slurb_tile%t_roof_0(nzt_roof:nzb_roof,i2,j2) )
    allocate( slurb_tile%t_roof_m(nzt_roof:nzb_roof,i2,j2) )
    allocate( slurb_tile%t_road_0(nzt_road:nzb_road,i2,j2) )
    allocate( slurb_tile%t_road_m(nzt_road:nzb_road,i2,j2) )

    if ( moist_physics )  then
       allocate( slurb_tile%tq_can(i2,j2))
       allocate( slurb_tile%tm_liq_roof(i2,j2) )
       allocate( slurb_tile%tm_liq_road(i2,j2) )
       allocate( slurb_tile%tm_roof_runoff(i2,j2) )
       allocate( slurb_tile%tm_road_runoff(i2,j2) )
       allocate( slurb_tile%tm_roof_precep(i2,j2) )
       allocate( slurb_tile%tm_road_precep(i2,j2) )

       allocate( slurb_tile%vpt_roof(i2,j2) )
       allocate( slurb_tile%vpt_road(i2,j2) )

       allocate( slurb_tile%q_roof(i2,j2) )
       allocate( slurb_tile%q_road(i2,j2) )
       allocate( slurb_tile%qs_roof(i2,j2) )
       allocate( slurb_tile%qs_road(i2,j2) )

       allocate( slurb_tile%qsws_can(i2,j2) )
       allocate( slurb_tile%qsws_roof(i2,j2) )
       allocate( slurb_tile%qsws_road(i2,j2) )
       allocate( slurb_tile%qsws_liq_roof(i2,j2) )
       allocate( slurb_tile%qsws_liq_road(i2,j2) )

       allocate( slurb_tile%c_liq_roof(i2,j2) )
       allocate( slurb_tile%c_liq_road(i2,j2) )

       allocate( slurb_tile%vpt_can(i2,j2) )

       allocate( slurb_tile%q1(i2,j2) )
       allocate( slurb_tile%vpt1(i2,j2) )

       allocate( slurb_tile%qsws_external(i2,j2) )

       allocate( slurb_tile%q_can_0(i2,j2) )
       allocate( slurb_tile%q_can_m(i2,j2) )
       allocate( slurb_tile%m_liq_roof_0(i2,j2) )
       allocate( slurb_tile%m_liq_roof_m(i2,j2) )
       allocate( slurb_tile%m_liq_road_0(i2,j2) )
       allocate( slurb_tile%m_liq_road_m(i2,j2) )

    endif

    allocate( slurb_tile%dt_max(i2,j2) )



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

    if ( facade_rah_doe )  then
       slurb_tile%rah_wall_a(:,:) = ieee_value(slurb_tile%rah_wall_a,ieee_signaling_nan)
       slurb_tile%rah_wall_b(:,:) = ieee_value(slurb_tile%rah_wall_b,ieee_signaling_nan)
       slurb_tile%rah_win_a(:,:) = ieee_value(slurb_tile%rah_win_a,ieee_signaling_nan)
       slurb_tile%rah_win_b(:,:) = ieee_value(slurb_tile%rah_win_b,ieee_signaling_nan)
    else
       slurb_tile%rah_facade(:,:) = ieee_value(slurb_tile%rah_facade,ieee_signaling_nan)
    endif

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

    if ( moist_physics )  then
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

    endif

    slurb_tile%dt_max(:,:) = ieee_value(slurb_tile%dt_max,ieee_signaling_nan)
#endif
end subroutine slurb_bulk_allocations

subroutine slurb_bulk_deallocations
    deallocate(fraction_slurb)

    deallocate(ln_z_z0_roof)
    deallocate(ln_z_z0h_roof)
    deallocate(ln_z_z0_urb)
    deallocate(pt_surface)
    deallocate(ln_z_z0_road)
    deallocate(ln_z_z0h_road)
    !-- Bulk allocation
    deallocate( slurb_tile%dz_roof)
    deallocate( slurb_tile%dz_wall)
    deallocate( slurb_tile%dz_road)
    deallocate( slurb_tile%dz_win)
    deallocate( slurb_tile%zw_win)

    deallocate( slurb_tile%t_c_urb)
    deallocate( slurb_tile%thl_rad_urb)
    deallocate( slurb_tile%t_h_urb)
    deallocate( slurb_tile%t_2m_urb)
    deallocate( slurb_tile%shf_urb)
    deallocate( slurb_tile%qsws_urb)
    deallocate( slurb_tile%ol_urb)
    deallocate( slurb_tile%rib_urb)
    deallocate( slurb_tile%ram_urb)
    deallocate( slurb_tile%usws_urb)
    deallocate( slurb_tile%vsws_urb)
    deallocate( slurb_tile%thlskin)
    deallocate( slurb_tile%qtskin)

    deallocate( slurb_tile%albedo_urb)
    deallocate( slurb_tile%emiss_urb)

    deallocate( slurb_tile%t_indoor)
    deallocate( slurb_tile%t_soil)

    deallocate( slurb_tile%tt_can)
    deallocate( slurb_tile%tt_wall_a)
    deallocate( slurb_tile%tt_wall_b)
    deallocate( slurb_tile%tt_win_a)
    deallocate( slurb_tile%tt_win_b)
    deallocate( slurb_tile%tt_roof)
    deallocate( slurb_tile%tt_road)

    deallocate( slurb_tile%pt_wall_a)
    deallocate( slurb_tile%pt_wall_b)
    deallocate( slurb_tile%pt_win_a)
    deallocate( slurb_tile%pt_win_b)
    deallocate( slurb_tile%pt_roof)
    deallocate( slurb_tile%pt_road)

    deallocate( slurb_tile%shf_can)
    deallocate( slurb_tile%shf_roof)
    deallocate( slurb_tile%shf_road)
    deallocate( slurb_tile%shf_wall_a)
    deallocate( slurb_tile%shf_wall_b)
    deallocate( slurb_tile%shf_win_a)
    deallocate( slurb_tile%shf_win_b)

    deallocate( slurb_tile%shf_external)
    deallocate( slurb_tile%shf_traffic)

    deallocate( slurb_tile%ghf_road)
    deallocate( slurb_tile%ghf_roof)
    deallocate( slurb_tile%ghf_wall_a)
    deallocate( slurb_tile%ghf_wall_b)
    deallocate( slurb_tile%ghf_win_a)
    deallocate( slurb_tile%ghf_win_b)

    deallocate( slurb_tile%rad_lw_in_urb)
    deallocate( slurb_tile%rad_sw_in_urb)
    deallocate( slurb_tile%rad_lw_out_urb)
    deallocate( slurb_tile%rad_sw_out_urb)

    deallocate( slurb_tile%rad_lw_net_urb)
    deallocate( slurb_tile%rad_sw_net_urb)

    deallocate( slurb_tile%rad_lw_net_can)

    deallocate( slurb_tile%rad_lw_net_roof)
    deallocate( slurb_tile%rad_sw_net_roof)
    deallocate( slurb_tile%rad_lw_net_road)
    deallocate( slurb_tile%rad_sw_net_road)
    deallocate( slurb_tile%rad_sw_in_road)
    deallocate( slurb_tile%rad_lw_net_wall_a)
    deallocate( slurb_tile%rad_sw_net_wall_a)
    deallocate( slurb_tile%rad_lw_net_wall_b)
    deallocate( slurb_tile%rad_sw_net_wall_b)
    deallocate( slurb_tile%rad_lw_net_win_a)
    deallocate( slurb_tile%rad_sw_net_win_a)
    deallocate( slurb_tile%rad_sw_in_win_a)
    deallocate( slurb_tile%rad_lw_net_win_b)
    deallocate( slurb_tile%rad_sw_net_win_b)
    deallocate( slurb_tile%rad_sw_in_win_b)

    deallocate( slurb_tile%pt_can)
    deallocate( slurb_tile%uv_abs_can)
    deallocate( slurb_tile%uv_eff_can)
    deallocate( slurb_tile%us_can)
    deallocate( slurb_tile%rib_can)
    deallocate( slurb_tile%ol_can)

    deallocate( slurb_tile%rib_roof)
    deallocate( slurb_tile%ol_roof)
    deallocate( slurb_tile%rib_road)
    deallocate( slurb_tile%ol_road)

    deallocate( slurb_tile%us_roof)
    deallocate( slurb_tile%us_road)

    deallocate( slurb_tile%hw_can)
    deallocate( slurb_tile%anisotropic_canyon)
    deallocate( slurb_tile%theta_can)
    deallocate( slurb_tile%h_bld)
    deallocate( slurb_tile%f_bld)
    deallocate( slurb_tile%f_bld_frn)
    deallocate( slurb_tile%f_win)
    deallocate( slurb_tile%svf_road)
    deallocate( slurb_tile%svf_wall)
    deallocate( slurb_tile%z0_urb)

    deallocate( slurb_tile%rah_roof)
    deallocate( slurb_tile%rah_road)
    deallocate( slurb_tile%rah_can)

    if ( facade_rah_doe )  then
       deallocate( slurb_tile%rah_wall_a)
       deallocate( slurb_tile%rah_wall_b)
       deallocate( slurb_tile%rah_win_a)
       deallocate( slurb_tile%rah_win_b)
    else
       deallocate( slurb_tile%rah_facade)
    endif

    deallocate( slurb_tile%lambda_roof)
    deallocate( slurb_tile%c_roof)
    deallocate( slurb_tile%albedo_roof)
    deallocate( slurb_tile%emiss_roof)
    deallocate( slurb_tile%z0_roof)
    deallocate( slurb_tile%z0h_roof)
    deallocate( slurb_tile%lambda_wall)
    deallocate( slurb_tile%c_wall)
    deallocate( slurb_tile%albedo_wall)
    deallocate( slurb_tile%emiss_wall)
    deallocate( slurb_tile%z0_wall)
    deallocate( slurb_tile%lambda_win)
    deallocate( slurb_tile%c_win)
    deallocate( slurb_tile%albedo_wall_win)
    deallocate( slurb_tile%albedo_win)
    deallocate( slurb_tile%emiss_win)
    deallocate( slurb_tile%transmissivity_win)
    deallocate( slurb_tile%absorption_win)
    deallocate( slurb_tile%lambda_road)
    deallocate( slurb_tile%c_road)
    deallocate( slurb_tile%albedo_road)
    deallocate( slurb_tile%emiss_road)
    deallocate( slurb_tile%z0_road)
    deallocate( slurb_tile%z0h_road)

    deallocate( slurb_tile%conductivity_roof)
    deallocate( slurb_tile%conductivity_wall)
    deallocate( slurb_tile%conductivity_win)
    deallocate( slurb_tile%conductivity_road)

    deallocate( slurb_tile%z_mo)
    deallocate( slurb_tile%z_mo_can)
    deallocate( slurb_tile%uv_abs_can_coef)
    deallocate( slurb_tile%wall_hor_a_ratio)


    deallocate( slurb_tile%lw_roof_coef)
    deallocate( slurb_tile%lw_road_coef)
    deallocate( slurb_tile%lw_wall_coef)
    deallocate( slurb_tile%lw_win_coef)
    deallocate( slurb_tile%sw_ref_denom)

    deallocate( slurb_tile%us_urb)
    deallocate( slurb_tile%uv_eff1)
    deallocate( slurb_tile%uv_abs1)
    deallocate( slurb_tile%pt1)

    deallocate( slurb_tile%t_can_0)
    deallocate( slurb_tile%t_can_m)
    deallocate( slurb_tile%t_wall_a_0)
    deallocate( slurb_tile%t_wall_a_m)
    deallocate( slurb_tile%t_wall_b_0)
    deallocate( slurb_tile%t_wall_b_m)
    deallocate( slurb_tile%t_win_a_0)
    deallocate( slurb_tile%t_win_a_m)
    deallocate( slurb_tile%t_win_b_0)
    deallocate( slurb_tile%t_win_b_m)
    deallocate( slurb_tile%t_roof_0)
    deallocate( slurb_tile%t_roof_m)
    deallocate( slurb_tile%t_road_0)
    deallocate( slurb_tile%t_road_m)

    if ( moist_physics )  then
       deallocate( slurb_tile%tq_can)
       deallocate( slurb_tile%tm_liq_roof)
       deallocate( slurb_tile%tm_liq_road)
       deallocate( slurb_tile%tm_roof_runoff)
       deallocate( slurb_tile%tm_road_runoff)
       deallocate( slurb_tile%tm_roof_precep)
       deallocate( slurb_tile%tm_road_precep)
       deallocate( slurb_tile%vpt_roof)
       deallocate( slurb_tile%vpt_road)

       deallocate( slurb_tile%q_roof)
       deallocate( slurb_tile%q_road)
       deallocate( slurb_tile%qs_roof)
       deallocate( slurb_tile%qs_road)

       deallocate( slurb_tile%qsws_can)
       deallocate( slurb_tile%qsws_roof)
       deallocate( slurb_tile%qsws_road)
       deallocate( slurb_tile%qsws_liq_roof)
       deallocate( slurb_tile%qsws_liq_road)

       deallocate( slurb_tile%c_liq_roof)
       deallocate( slurb_tile%c_liq_road)

       deallocate( slurb_tile%vpt_can)

       deallocate( slurb_tile%q1)
       deallocate( slurb_tile%vpt1)

       deallocate( slurb_tile%qsws_external)

       deallocate( slurb_tile%q_can_0)
       deallocate( slurb_tile%q_can_m)
       deallocate( slurb_tile%m_liq_roof_0)
       deallocate( slurb_tile%m_liq_roof_m)
       deallocate( slurb_tile%m_liq_road_0)
       deallocate( slurb_tile%m_liq_road_m)
    endif
end subroutine slurb_bulk_deallocations


    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    ! Swap timelevel of the SLUrb model.
    !--------------------------------------------------------------------------------------------------!
 subroutine slurb_set_previous_timestep()
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

 end subroutine slurb_set_previous_timestep

  !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> This function computes the magnus formula (Press et al., 1992).
    !> The magnus formula is needed to calculate the saturation vapor pressure.
    !--------------------------------------------------------------------------------------------------!
 function magnus( t )
    !$ACC ROUTINE SEQ

    IMPLICIT NONE
!$omp declare target

    real(field_r), intent(in) ::  t  !< temperature (K)

    real(field_r) ::  magnus

    !
    !-- Saturation vapor pressure for a specific temperature:
    magnus =  611.2_field_r * EXP( 17.62_field_r * ( t - 273.15_field_r ) / ( t - 29.65_field_r  ) )

 end function magnus

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Computes a steady-state solution for 1D heat equation using Gauss-Seidel iteration.
    !> This is used to initialize the material temperatures for roofs, walls, windows and roads,
    !> shortening the time required for the spinup. For windows, SW absorption is not considered.
    !--------------------------------------------------------------------------------------------------!
 !TODOSELF pure
 function calc_1d_heat_equation( result_size , t_bc_1, t_bc_2, lambda ) RESULT( t_result )

    integer, intent(in) ::  result_size  !< output target size

    real(field_r), intent(in) ::  t_bc_1  !< outer t boundary condition
    real(field_r), intent(in) ::  t_bc_2  !< inner t boundary condition

    real(field_r), dimension(:), intent(in) ::  lambda  !< total layer heat conductivity

    integer ::  ix  !< iteration counter
    integer ::  kx  !< layer running index

    real(field_r), PARAMETER ::  omega = 1.0_field_r   !< relaxation to control convergence
    real(field_r), PARAMETER ::  tol = 1.0E-6_field_r  !< maximum residual for convergence

    real(field_r) ::  res    !< iteration residual for convergence check
    real(field_r) ::  t_old  !< previous t of layer for convergence check

    real(field_r), dimension(result_size) ::  t_result  !< result t profile
    real(field_r), dimension(1:result_size+1) ::  t     !< intermediate t array containing also the BCs


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
    do  ix = 1, 1000
       do  kx = LBOUND( t, 1 )+1, UBOUND( t, 1 )-1
          t_old = t(kx)
          res = 0.0_field_r
          t(kx) = ( lambda(kx) * ( t(kx+1) - t(kx) ) + lambda(kx-1) * ( t(kx-1) - t(kx) ) ) /      &
                  ( lambda(kx) + lambda(kx-1) ) * omega + t(kx)
          res = MAX( res, ABS( t(kx) - t_old ) )
       enddo
    !
    !--    Check for convergence using the stored maximum residual.
       if ( res < tol )  EXIT
    enddo
    !
    !-- Return the solution.
    t_result(:) = t(LBOUND( t, 1 ):UBOUND( t, 1 )-1)

 end function calc_1d_heat_equation

end module modslurbhelpers
