!> \file modslurb_energybalance.f90
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
module modslurb_energybalance
    use modslurbdata
    use modfields, only : exnf, rhobf, ql0, rhof, rho_air_zw => rhobf
    use modmicrodata, only : precep, imicro !TODOSELF TEST MICRO
    use modsurface,  only : ps
    use modslurbhelpers, only: magnus
    real:: rho_cp  !< cp * rho (J m^-3 K^-1)

    contains


    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Surface and subsurface energy balance computations of roofs, walls, windows and roads.
    !--------------------------------------------------------------------------------------------------!
 subroutine slurb_energy_balance_model
   use modglobal, only : i1, j1, cp, rlv, rhow, rk3step, rdt, ep

   implicit none
    integer ::  i       !< loop index (x-direction)
    integer ::  j       !< loop index (y-direction)
    integer ::  k_topo  !< k index of topography
    integer ::  k_atm   !< k index of the first atmospheric level




    k_topo = 1
    k_atm = 1
    rho_cp = cp * rho_air_zw(k_topo) ! TODO, check if this is the right density to use for roof/wall/window/road energy balance calculations, or if we should use a different height level..
   do j=2,j1
      do i=2,i1
      !  k_topo = topo_top_ind(j,i,0)
      !  k_atm = topo_top_ind(j,i,0) + 1


    !
    !--    Call specific models for all the facets.
       call roof_model
       call wall_model
       if ( slurb_tile%f_win(i,j) /= 0.0_field_r )  call window_model
       call road_model
      enddo
   enddo

   contains


    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Computes the new surface prognostic temperature for current time step using RK3.
    !--------------------------------------------------------------------------------------------------!
 subroutine calc_surf_t_p ( t_m, t_0, tt_current, coef_1, coef_2, c )
    use modglobal, only : rk3step, rdt
    implicit none

    real(field_r), intent(in) ::  c       !< total layer heat capacity (J m^-2 K^-1)
    real(field_r), intent(in) ::  coef_1  !< coefficient A in the prognostic equation (W m^-2)
    real(field_r), intent(in) ::  coef_2  !< coefficient B in the prognostic equation (W m^-2 K^-1)
    real(field_r), intent(in) ::  t_m       !< current layer temperature (K)

    ! real(field_r) :: tend
    real(field_r), intent(inout) ::  t_0  !< new layer temperature (K)

    real(field_r), intent(inout) ::  tt_current  !< current temperature tendency (K s^-1)

    real(field_r) ::  tt_new  !< new temperature tendency (K s^-1)
    real(field_r) :: t_new_implicit

    real :: rk3coef
    real :: rdt3 

    rdt3 = rdt / 3

    rk3coef = rdt / (4. - dble(rk3step))
    !-- Compute the RK3 tendency for next time step.
    if ( c /= 0.0_field_r )  then
        t_new_implicit = ( ( coef_1 * (rk3coef) + c * t_0 )  / ( c + coef_2 * (rk3coef)  ))
        tt_new = (t_new_implicit - t_0) / rk3coef
        t_0 = t_new_implicit

        tt_current = tt_new


    endif

 end subroutine calc_surf_t_p
    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Computes the new layer prognostic temperature by solving the Fourier diffusion equation.
    !--------------------------------------------------------------------------------------------------!
 subroutine calc_heat_diffusion ( t_m, t_0, tt_current, c, lambda, t_bc, sw_in, phi )
    use modglobal, only : rk3step, rdt
    real(field_r), intent(in) ::  t_bc  !< temperature boundary condition (K)

    real(field_r), intent(in), OPTIONAL ::  sw_in  !< incoming shortwave radiation for windows

    real(field_r), dimension(:), intent(in) ::  c       !< total heat capacity of the layer (J m^-2 K^-1)
    real(field_r), dimension(:), intent(in) ::  lambda  !< total heat conductivity between layers (W m^-2 K^-1)
    real(field_r), dimension(:), intent(in) ::  t_m       !< current time level temperature (K)

    real(field_r), dimension(:), intent(in), OPTIONAL ::  phi  !< fraction of incoming shortwave radiation absorbed at window layer

    real(field_r), dimension(:), intent(inout) ::  t_0  !< new layer temperature (K)

    real(field_r), dimension(:), intent(inout) ::  tt_current  !< current temperature tendency (K s^-1)

    integer ::  k  !< material layer loop index

    real(field_r) ::  tt_new  !<  new temperature tendency (K s^-1)

    real, allocatable :: temp1(:)

    real :: rk3coef

    rk3coef = rdt / (4. - dble(rk3step))

    !
    !-- Loop through non-boundary layers of the material.
    !-- @todo Split loop into three to move IFs out for better vecotrization.
    do  k = LBOUND( t_0, 1 ) + 1, UBOUND( t_0, 1 )
    !
    !--    New prognostic layer temperature.
    !--    Compute the t between neighbouring layers.
      if ( k /= UBOUND( t_0 , 1 ) )  then
         tt_new = ( 1.0_field_r / c(k) ) * ( lambda(k) * ( t_0(k+1) - t_0(k) ) +                           &
                  lambda(k-1) * ( t_0(k-1) - t_0(k) ) )
      else
   !
   !--    Use a constant value boundary condition (skin temperature) for the innermost layer.
         tt_new = ( 1.0_field_r / c(k) ) * ( lambda(k) * ( t_bc - t_0(k) ) +                             &
                  lambda(k-1) * ( t_0(k-1) - t_0(k) ) )
      endif
      !
      !--    Add tendency from absorbed shortwave radiation.
      if ( PRESENT( sw_in ) )  then
          tt_new = tt_new + ( 1.0_field_r / c(k) ) * sw_in * phi(k)
      endif

      t_0(k) = t_m(k) + (rk3coef) * ( tt_new )

      tt_current(k) = tt_new

    enddo

 end subroutine calc_heat_diffusion



    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Models the surface energy balance and subsurface heat diffusion for roofs.
    !--------------------------------------------------------------------------------------------------!
 subroutine roof_model
   use modglobal, only : rk3step, rdt
    real(field_r) ::  coef_1              !< coefficient A of the prognostic equation
    real(field_r) ::  coef_2              !< coefficient B of the prognostic equation
    real(field_r) ::  dq_s_dt             !< water vapour mixing ratio tendency
    real(field_r) ::  e_s                 !< saturation water vapour pressure
    real(field_r) ::  e_s_dt              !< saturation water vapour pressure tendency
    real(field_r) ::  f_shf               !< factor for the roof sensible heat flux (W m^-2 K^-1)
    real(field_r) ::  f_qsws_liq          !< factor for the latent heat flux from/to liquid water reservoir (W m^-2)
    real(field_r) ::  tm_new              !< new liquid water reservoir tendency (m s^-1)
    real(field_r) ::  tm_new_limited      !< new liquid water reservoir tendency limited by the max reservoir (m s^-1)
    real :: rk3coef

    rk3coef = rdt / (4. - dble(rk3step))


    !
    !-- Surface sensible heat flux factor.
    f_shf = rho_cp / slurb_tile%rah_roof(i,j)

    !
    !-- Compute the nominator and denominator coefficients in
    !-- the prognostic equation for the moist case.
    if ( moist_physics )  then
    !
    !--    Computation of factor for the latent heat flux due to
    !--    liquid water reservoir evaporation/condensation.
       e_s = ps * slurb_tile%qs_roof(i,j) / ( slurb_tile%qs_roof(i,j) + ep )

    !
    !--    In case of evaporation, evaporate only for the liquid water coverage area,
    !--    in case of condensation, use the total surface.
       if ( slurb_tile%qs_roof(i,j) > slurb_tile%q1(i,j) )  then
          f_qsws_liq = rho_lv * slurb_tile%c_liq_roof(i,j) / slurb_tile%rah_roof(i,j)
       else
          f_qsws_liq = rho_lv / slurb_tile%rah_roof(i,j)
       endif

       e_s_dt = e_s * ( 17.62_field_r / ( slurb_tile%t_roof_0(nzt_roof,i,j) -  29.65_field_r ) -                       &
                        17.62_field_r * ( slurb_tile%t_roof_0(nzt_roof,i,j) - 273.15_field_r ) /                       &
                        ( slurb_tile%t_roof_0(nzt_roof,i,j) - 29.65_field_r )**2                                  &
                      )

       dq_s_dt = ep * e_s_dt / ( ps - e_s_dt )

    !
    !--    The coefficients for the moist prognostic equation for temperature.
       coef_1 = slurb_tile%rad_sw_net_roof(i,j) + slurb_tile%rad_lw_net_roof(i,j)                                  &
                - 3.0_field_r * slurb_tile%lw_roof_coef(1,i,j) * slurb_tile%t_roof_0(nzt_roof,i,j)**4                     &
                + f_shf * slurb_tile%pt1(i,j)                                                              &
                + f_qsws_liq * ( slurb_tile%q1(i,j) - slurb_tile%qs_roof(i,j)                                      &
                                 + dq_s_dt * slurb_tile%t_roof_0(nzt_roof,i,j) )                             &
                + slurb_tile%conductivity_roof(nzt_roof,i,j) * slurb_tile%t_roof_0(nzt_roof+1,i,j)

       coef_2 = -4.0_field_r * slurb_tile%lw_roof_coef(1,i,j) * slurb_tile%t_roof_0(nzt_roof,i,j)**3                      &
                + f_shf * (1 / exnf(k_topo))                                                          &
                + f_qsws_liq * dq_s_dt                                                             &
                + slurb_tile%conductivity_roof(nzt_roof,i,j)

    else
    !
    !-- The coefficients for the dry prognostic equation for temperature.
       coef_1 = slurb_tile%rad_sw_net_roof(i,j) + slurb_tile%rad_lw_net_roof(i,j)                                  &
                -3.0_field_r * slurb_tile%lw_roof_coef(1,i,j) * slurb_tile%t_roof_0(nzt_roof,i,j)**4                      &
                + f_shf * slurb_tile%pt1(i,j)                                                              &
                + slurb_tile%conductivity_roof(nzt_roof,i,j) * slurb_tile%t_roof_0(nzt_roof+1,i,j)

       coef_2 = -4.0_field_r * slurb_tile%lw_roof_coef(1,i,j) * slurb_tile%t_roof_0(nzt_roof,i,j)**3                      &
                + f_shf * (1 / exnf(k_topo))                                                          &
                + slurb_tile%conductivity_roof(nzt_roof,i,j)
    endif

    call calc_surf_t_p( slurb_tile%t_roof_m(nzt_roof,i,j), slurb_tile%t_roof_0(nzt_roof,i,j),                        &
                        slurb_tile%tt_roof(nzt_roof,i,j), coef_1, coef_2, slurb_tile%c_roof(nzt_roof,i,j) )

    !
    !-- Explicit solution of the Fourier heat equation for the subsurface layers.
    call calc_heat_diffusion( slurb_tile%t_roof_m(:,i,j), slurb_tile%t_roof_0(:,i,j), slurb_tile%tt_roof(:,i,j),             &
                              slurb_tile%c_roof(:,i,j), slurb_tile%conductivity_roof(:,i,j), slurb_tile%t_indoor(i,j) )

    !
    !-- Compute the diagnostic fluxes for the roof surface.
    slurb_tile%ghf_roof(i,j) = slurb_tile%conductivity_roof(nzb_roof,i,j) *                                        &
                       ( slurb_tile%t_roof_0(nzb_roof,i,j) - slurb_tile%t_indoor(i,j) )

    slurb_tile%pt_roof(i,j) = slurb_tile%t_roof_0(nzt_roof,i,j) * (1 / exnf(k_topo))

    slurb_tile%shf_roof(i,j) = -f_shf * ( slurb_tile%pt1(i,j) - slurb_tile%pt_roof(i,j) )

    !
    !-- Update longwave radiative flux following linearization.
    slurb_tile%rad_lw_net_roof(i,j) = slurb_tile%rad_lw_net_roof(i,j)                                              &
                              + slurb_tile%lw_roof_coef(1,i,j) * slurb_tile%t_roof_m(nzt_roof,i,j)**4                &
                              - 4.0_field_r * slurb_tile%lw_roof_coef(1,i,j) * slurb_tile%t_roof_0(nzt_roof,i,j)**3       &
                              * ( slurb_tile%t_roof_m(nzt_roof,i,j) - slurb_tile%t_roof_0(nzt_roof,i,j) )

    !
    !-- Compute the water vapor flux from/to liquid water reservoir and the prognostic reservoir level.
    if ( moist_physics )  then
       slurb_tile%qsws_liq_roof(i,j) = -f_qsws_liq * ( slurb_tile%q1(i,j) - slurb_tile%qs_roof(i,j) +                      &
                                               dq_s_dt * slurb_tile%t_roof_m(nzt_roof,i,j) -                 &
                                               dq_s_dt * slurb_tile%t_roof_0(nzt_roof,i,j)                 &
                                             )

       slurb_tile%qsws_roof(i,j) = slurb_tile%qsws_liq_roof(i,j)
    !
    !
    !--    Modification due to precipitiation. If the liquid reservoir is full, the liquid water
    !--    is assumed to be drained into the drainage system (liquid water is not conserved).
    !--    The precipitation flux is not included in the surface-atmosphere latent heat flux (qsws).
       if (imicro == 0 .or. imicro == 1) then
            slurb_tile%qsws_liq_roof(i,j) = slurb_tile%qsws_roof(i,j)
        else
            !   if ( slurb_tile%m_liq_roof_0(i,j) < m_liq_max_roof )  then
            slurb_tile%tm_roof_precep(i,j) = precep(i,j,k_atm)
            slurb_tile%qsws_liq_roof(i,j) = (slurb_tile%qsws_roof(i,j) - slurb_tile%tm_roof_precep(i,j) * rhof(k_atm) * rlv)
  

          !todoself even morme assume precipitation
       endif
          !todoself assume precipitation
    !
    !--    Compute the total latent heat flux.
       slurb_tile%qsws_roof(i,j) = slurb_tile%qsws_roof(i,j)
    !
    !--    Compute the prognostic liquid water reservoir.
       tm_new = - slurb_tile%qsws_liq_roof(i,j) * drho_l_lv
       slurb_tile%m_liq_roof_0(i,j)  = slurb_tile%m_liq_roof_m(i,j) + rk3coef * tm_new
    !
    !--    Check if the liquid water reservoir is overfull. If so, drain excess to the
    !--    assumed drainage system (water is not conserved here).



    if ((slurb_tile%m_liq_roof_0(i,j) > m_liq_max_roof)) then
        ! tm_new = (m_liq_max_roof - slurb_tile%m_liq_roof_m(i,j)) / rk3coef
        slurb_tile%tm_roof_runoff(i,j) = (slurb_tile%m_liq_roof_0(i,j) - m_liq_max_roof) / rk3coef
        slurb_tile%m_liq_roof_0(i,j) = m_liq_max_roof
    else
        slurb_tile%tm_roof_runoff(i,j) = 0
    endif

    !
    !--    Check for negative water reservoir. @todo store the removed water as runoff for output.
       slurb_tile%m_liq_roof_0(i,j) = MAX( slurb_tile%m_liq_roof_0(i,j), 0.0_field_r )
    !
       
       slurb_tile%tm_liq_roof(i,j) = tm_new
    !
    !--    Compute the new liquid water coverage.
       slurb_tile%c_liq_roof(i,j) = MIN( 1.0_field_r, ( slurb_tile%m_liq_roof_0(i,j) / m_liq_max_roof )**0.67 )
    !
    !--    Compute new saturation mixing ratio.
       e_s = magnus( MIN( slurb_tile%t_roof_0(nzt_roof,i,j), 333.15_field_r ) )
       slurb_tile%qs_roof(i,j) = ep * e_s / ( ps - e_s )
    !
    !--    Calculate new mixing ratio and vpt at roof surface.
       slurb_tile%q_roof(i,j) = q_surf( slurb_tile%qs_roof(i,j), slurb_tile%rah_roof(i,j), slurb_tile%q1(i,j), f_qsws_liq )
       slurb_tile%vpt_roof(i,j) = slurb_tile%pt_roof(i,j) * ( 1.0_field_r + 0.61_field_r * slurb_tile%q_roof(i,j) )

    endif

 end subroutine roof_model


    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Models the surface energy balance and subsurface heat diffusion for roads.
    !--------------------------------------------------------------------------------------------------!
 subroutine road_model
   use modglobal, only : rk3step, rdt
    real(field_r) ::  coef_1              !< coefficient A of the prognostic equation
    real(field_r) ::  coef_2              !< coefficient B of the prognostic equation
    real(field_r) ::  dq_s_dt             !< water vapour mixing ratio tendency
    real(field_r) ::  e_s                 !< saturation water vapour pressure
    real(field_r) ::  e_s_dt              !< saturation water vapour pressure tendency
    real(field_r) ::  f_shf               !< factor for the road sensible heat flux (W m^-2 K^-1)
    real(field_r) ::  f_qsws_liq          !< factor for the latent heat flux from/to liquid water reservoir (W m^-2)
    real(field_r) ::  tm_new              !< new liquid water reservoir tendency (m s^-1)
    real(field_r) ::  tm_new_limited      !< new liquid water reservoir tendency limited by the max reservoir (m s^-1)
    real :: rk3coef

    rk3coef = rdt / (4. - dble(rk3step))

    !
    !-- Surface sensible heat flux factor.
    f_shf = rho_cp / slurb_tile%rah_road(i,j)
    !
    !-- Compute the nominator and denominator coefficients in
    !-- the prognostic equation for the moist case.
    if ( moist_physics )  then
    !
    !--    Computation of factor for the latent heat flux due to
    !--    liquid water reservoir evaporation/condensation.
       e_s = ps * slurb_tile%qs_road(i,j) / ( slurb_tile%qs_road(i,j) + ep )

    !
    !--    In case of evaporation, evaporate only for the liquid water coverage area,
    !--    in case of condensation, use the total surface.
       if ( slurb_tile%qs_road(i,j) > slurb_tile%q_can_0(i,j) )  then
          f_qsws_liq = rho_lv * slurb_tile%c_liq_road(i,j) / slurb_tile%rah_road(i,j)
       else
          f_qsws_liq = rho_lv / slurb_tile%rah_road(i,j)
       endif

       e_s_dt = e_s * ( 17.62_field_r / ( slurb_tile%t_road_0(nzt_road,i,j) - 29.65_field_r ) -                        &
                        17.62_field_r * ( slurb_tile%t_road_0(nzt_road,i,j) - 273.15_field_r ) /                       &
                        ( slurb_tile%t_road_0(nzt_road,i,j) - 29.65_field_r)**2                                   &
                      )

       dq_s_dt = ep * e_s_dt / ( ps - e_s_dt )

    !
    !--    The coefficients for the moist prognostic equation for temperature. For the longwave balance,
    !--    both direct emission and the effect of backreflection are linearized.
       coef_1 = slurb_tile%rad_sw_net_road(i,j) + slurb_tile%rad_lw_net_road(i,j)                                  &
                -3.0_field_r * slurb_tile%lw_road_coef(1,i,j) * slurb_tile%t_road_0(nzt_road,i,j)**4                      &
                + f_shf * slurb_tile%t_can_0(i,j)                                                            &
                + f_qsws_liq * ( slurb_tile%q_can_0(i,j) - slurb_tile%qs_road(i,j)                                   &
                                 + dq_s_dt * slurb_tile%t_road_0(nzt_road,i,j) )                             &
                + slurb_tile%conductivity_road(nzt_road,i,j) * slurb_tile%t_road_0(nzt_road+1,i,j)

       coef_2 = -4.0_field_r * slurb_tile%lw_road_coef(1,i,j) * slurb_tile%t_road_0(nzt_road,i,j)**3                      &
                + f_shf                                                                            &
                + f_qsws_liq * dq_s_dt                                                             &
                + slurb_tile%conductivity_road(nzt_road,i,j)

    else
    !
    !--    The coefficients for the dry prognostic equation for temperature.
       coef_1 = slurb_tile%rad_sw_net_road(i,j) + slurb_tile%rad_lw_net_road(i,j)                                  &
                -3.0_field_r * slurb_tile%lw_road_coef(1,i,j) * slurb_tile%t_road_0(nzt_road,i,j)**4                      &
                + f_shf * slurb_tile%t_can_0(i,j)                                                            &
                + slurb_tile%conductivity_road(nzt_road,i,j) * slurb_tile%t_road_0(nzt_road+1,i,j)

       coef_2 = -4.0_field_r * slurb_tile%lw_road_coef(1,i,j) * slurb_tile%t_road_0(nzt_road,i,j)**3                      &
                + f_shf                                                                            &
                + slurb_tile%conductivity_road(nzt_road,i,j)
    endif

    call calc_surf_t_p( slurb_tile%t_road_m(nzt_road,i,j), slurb_tile%t_road_0(nzt_road,i,j),                        &
                        slurb_tile%tt_road(nzt_road,i,j), coef_1, coef_2, slurb_tile%c_road(nzt_road,i,j) )

    !
    !-- Heat diffusion through subsurface layers.
    call calc_heat_diffusion( slurb_tile%t_road_m(:,i,j), slurb_tile%t_road_0(:,i,j), slurb_tile%tt_road(:,i,j),             &
                              slurb_tile%c_road(:,i,j), slurb_tile%conductivity_road(:,i,j), slurb_tile%t_soil(i,j) )


    slurb_tile%shf_road(i,j) = -f_shf * ( slurb_tile%t_can_m(i,j) - slurb_tile%t_road_0(nzt_road,i,j) )

    slurb_tile%pt_road(i,j)  = slurb_tile%t_road_0(nzt_road,i,j) * (1 / exnf(k_topo))

    slurb_tile%ghf_road(i,j) = slurb_tile%conductivity_road(nzb_road,i,j) *                                        &
                       ( slurb_tile%t_road_0(nzb_road,i,j) - slurb_tile%t_soil(i,j) )

    !
    !-- Update longwave radiative flux following linearization.
    slurb_tile%rad_lw_net_road(i,j) = slurb_tile%rad_lw_net_road(i,j)                                              &
                              + slurb_tile%lw_road_coef(1,i,j) * slurb_tile%t_road_m(nzt_road,i,j)**4                &
                              - 4.0_field_r * slurb_tile%lw_road_coef(1,i,j) * slurb_tile%t_road_m(nzt_road,i,j)**3       &
                              * ( slurb_tile%t_road_m(nzt_road,i,j) - slurb_tile%t_road_0(nzt_road,i,j) )

    !
    !-- Compute the water vapor flux from/to liquid water reservoir and the prognostic reservoir level.
    if ( moist_physics )  then
       slurb_tile%qsws_liq_road(i,j) = -f_qsws_liq * ( slurb_tile%q_can_0(i,j) - slurb_tile%qs_road(i,j) +                   &
                                               dq_s_dt * slurb_tile%t_road_m(nzt_road,i,j) -                 &
                                               dq_s_dt * slurb_tile%t_road_0(nzt_road,i,j)                 &
                                             )

       slurb_tile%qsws_road(i,j) = slurb_tile%qsws_liq_road(i,j)
    !
    !--    Modification due to precipitiation. If the liquid reservoir is full, the liquid water
    !--    is assumed to be drained into the drainage system (liquid water is not conserved).
    !--    The precipitation flux is not included in the surface-atmosphere latent heat flux (qsws).
    if (imicro == 0 .or. imicro == 1) then ! this should be the same as if PRECIPITATION
        slurb_tile%qsws_liq_road(i,j) = slurb_tile%qsws_road(i,j)
    else
        ! if ( slurb_tile%m_liq_road_0(i,j) < m_liq_max_road )  then
            slurb_tile%tm_road_precep(i,j) = precep(i,j,k_atm)
            slurb_tile%qsws_liq_road(i,j) = (slurb_tile%qsws_road(i,j) - slurb_tile%tm_road_precep(i,j) * rhof(k_atm) * rlv)
        ! endif

          !todoself even morme assume precipitation
    endif
       ! liquid water reservoir is in m^3/m^2, rain rate in m/s (m^3/m^2 /s)
    !
    !--    Compute the total latent heat flux.
       slurb_tile%qsws_road(i,j) = slurb_tile%qsws_road(i,j)
    !
    !--    Compute the prognostic liquid water reservoir.
       tm_new = - slurb_tile%qsws_liq_road(i,j) * drho_l_lv
       slurb_tile%m_liq_road_0(i,j) = slurb_tile%m_liq_road_m(i,j) + rk3coef * tm_new
    !
    !--    Check if the liquid water reservoir is overfull. If so, drain excess to the
    !--    assumed drainage system (water is not conserved here).
    if ((slurb_tile%m_liq_road_0(i,j) > m_liq_max_road)) then
        slurb_tile%tm_road_runoff(i,j) = (slurb_tile%m_liq_road_0(i,j) - m_liq_max_road) / rk3coef
        slurb_tile%m_liq_road_0(i,j) = m_liq_max_road
    else
        slurb_tile%tm_road_runoff(i,j) = 0
    endif
    !
    !--    Check for negative water reservoir. Should we adjust qsws_road accordingly?
       slurb_tile%m_liq_road_0(i,j) = MAX( slurb_tile%m_liq_road_0(i,j), 0.0_field_r )
    !
    !--    Compute RK3 tendency
       slurb_tile%tm_liq_road(i,j) = tm_new
    !
    !--    Compute the new liquid water coverage.
       slurb_tile%c_liq_road(i,j) = MIN( 1.0_field_r, ( slurb_tile%m_liq_road_0(i,j) / m_liq_max_road )**0.67 )
    !
    !--    Compute new saturation mixing ratio.
       e_s = magnus( MIN( slurb_tile%t_road_0(nzt_road,i,j), 333.15_field_r ) )
       slurb_tile%qs_road(i,j) = ep * e_s / ( ps - e_s )
    !
    !--    Calculate new mixing ratio and vpt at road surface.
       slurb_tile%q_road(i,j) = q_surf( slurb_tile%qs_road(i,j), slurb_tile%rah_road(i,j), slurb_tile%q1(i,j), f_qsws_liq )
       slurb_tile%vpt_road(i,j) = slurb_tile%pt_road(i,j) * ( 1.0_field_r + 0.61_field_r * slurb_tile%q_road(i,j) )

    endif

 end subroutine road_model


    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Models the surface energy balance and subsurface heat diffusion for both walls.
    !--------------------------------------------------------------------------------------------------!
 subroutine wall_model
   

   real(field_r) ::  coef_1   !< coefficient A of the prognostic equation
   real(field_r) ::  coef_2   !< coefficient B of the prognostic equation
   real(field_r) ::  f_shf_a  !< factor for the wall surface heat flux (W m^-2 K^-1)
   real(field_r) ::  f_shf_b  !< factor for the wall surface heat flux (W m^-2 K^-1)


    if ( facade_rah_doe )  then
       f_shf_a = rho_cp / slurb_tile%rah_wall_a(i,j)
       if ( slurb_tile%anisotropic_canyon(i,j) )  f_shf_b = rho_cp / slurb_tile%rah_wall_b(i,j)
    else
       f_shf_a = rho_cp / slurb_tile%rah_facade(i,j)
       if ( slurb_tile%anisotropic_canyon(i,j) )  f_shf_b = f_shf_a
    endif

    !
    !-- The coefficients for the moist prognostic equation for temperature. For the longwave balance,
    !-- both direct emission and the effect of backreflection are linearized. The linearization depends
    !-- if the canyon is isotropic or not as an average backreflection is used for isotropic canyons.
    !-- We consider the walls are dry in all cases, so moist physical processes are not considered.
    if ( slurb_tile%anisotropic_canyon(i,j) )  then
       coef_1 = slurb_tile%rad_sw_net_wall_a(i,j) + slurb_tile%rad_lw_net_wall_a(i,j)                              &
                - 3.0_field_r * slurb_tile%lw_wall_coef(1,i,j) * slurb_tile%t_wall_a_0(nzt_wall,i,j)**4                   &
                + f_shf_a * slurb_tile%t_can_0(i,j)                                                          &
                + slurb_tile%conductivity_wall(nzt_wall,i,j) * slurb_tile%t_wall_a_0(nzt_wall+1,i,j)

       coef_2 = -4.0_field_r * slurb_tile%lw_wall_coef(1,i,j) * slurb_tile%t_wall_a_0(nzt_wall,i,j)**3                    &
                + f_shf_a                                                                          &
                + slurb_tile%conductivity_wall(nzt_wall,i,j)

       call calc_surf_t_p(slurb_tile%t_wall_a_m(nzt_wall,i,j), slurb_tile%t_wall_a_0(nzt_wall,i,j),                  &
                          slurb_tile%tt_wall_a(nzt_wall,i,j), coef_1, coef_2, slurb_tile%c_wall(nzt_wall,i,j) )

       coef_1 = slurb_tile%rad_sw_net_wall_b(i,j) + slurb_tile%rad_lw_net_wall_b(i,j)                              &
                - 3.0_field_r * slurb_tile%lw_wall_coef(1,i,j) * slurb_tile%t_wall_b_0(nzt_wall,i,j)**4                   &
                + f_shf_b * slurb_tile%t_can_0(i,j)                                                          &
                + slurb_tile%conductivity_wall(nzt_wall,i,j) * slurb_tile%t_wall_b_0(nzt_wall+1,i,j)

       coef_2 = -4.0_field_r * slurb_tile%lw_wall_coef(1,i,j) * slurb_tile%t_wall_b_0(nzt_wall,i,j)**3                    &
                + f_shf_b                                                                          &
                + slurb_tile%conductivity_wall(nzt_wall,i,j)

       call calc_surf_t_p( slurb_tile%t_wall_b_m(nzt_wall,i,j), slurb_tile%t_wall_b_0(nzt_wall,i,j),                 &
                           slurb_tile%tt_wall_b(nzt_wall,i,j), coef_1, coef_2, slurb_tile%c_wall(nzt_wall,i,j) )
    else
    !
    !--    In case of isotropic canyon, wall A and B temperatures are averaged, and thus the prognostic
    !--    equation for t_wall_a is representative of both of the walls. Thus both the terms for
    !--    t_wall_a as well as for t_wall_b in the longwave radiation balance has a dependency on
    !--    the surface temperature.
       coef_1 = slurb_tile%rad_sw_net_wall_a(i,j) + slurb_tile%rad_lw_net_wall_a(i,j)                              &
                - 3.0_field_r * ( slurb_tile%lw_wall_coef(1,i,j) + slurb_tile%lw_wall_coef(3,i,j) )                     &
                   * slurb_tile%t_wall_a_0(nzt_wall,i,j)**4                                                  &
                + f_shf_a * slurb_tile%t_can_0(i,j)                                                          &
                + slurb_tile%conductivity_wall(nzt_wall,i,j) * slurb_tile%t_wall_a_0(nzt_wall+1,i,j)

       coef_2 = -4.0_field_r * ( slurb_tile%lw_wall_coef(1,i,j) + slurb_tile%lw_wall_coef(3,i,j) )                      &
                   * slurb_tile%t_wall_a_0(nzt_wall,i,j)**3                                                  &
                + f_shf_a                                                                          &
                + slurb_tile%conductivity_wall(nzt_wall,i,j)

       call calc_surf_t_p( slurb_tile%t_wall_a_m(nzt_wall,i,j), slurb_tile%t_wall_a_0(nzt_wall,i,j),                 &
                           slurb_tile%tt_wall_a(nzt_wall,i,j), coef_1, coef_2, slurb_tile%c_wall(nzt_wall,i,j) )
    endif
    slurb_tile%pt_wall_a(i,j)  = slurb_tile%t_wall_a_0(nzt_wall,i,j) * (1 / exnf(k_topo))
    slurb_tile%shf_wall_a(i,j) = -f_shf_a * ( slurb_tile%t_can_m(i,j) - slurb_tile%t_wall_a_0(nzt_wall,i,j) )

    !
    !-- Heat diffusion through subsurface layers.
    call calc_heat_diffusion( slurb_tile%t_wall_a_m(:,i,j), slurb_tile%t_wall_a_0(:,i,j), slurb_tile%tt_wall_a(:,i,j),       &
                              slurb_tile%c_wall(:,i,j), slurb_tile%conductivity_wall(:,i,j), slurb_tile%t_indoor(i,j) )

    slurb_tile%ghf_wall_a(i,j) = slurb_tile%conductivity_wall(nzb_wall,i,j) *                                      &
                         ( slurb_tile%t_wall_a_0(nzb_wall,i,j) - slurb_tile%t_indoor(i,j) )

    !
    !-- Same treatment for wall B if this is an anisotropic canyon, otherwise copy.
    if ( slurb_tile%anisotropic_canyon(i,j) )  then
       slurb_tile%pt_wall_b(i,j)  = slurb_tile%t_wall_b_0(nzt_wall,i,j) * (1 / exnf(k_topo))
       slurb_tile%shf_wall_b(i,j) = -f_shf_b * ( slurb_tile%t_can_m(i,j) - slurb_tile%t_wall_b_0(nzt_wall,i,j) )


       slurb_tile%ghf_wall_b(i,j) = slurb_tile%conductivity_wall(nzb_wall,i,j) *                                   &
                            ( slurb_tile%t_wall_b_0(nzb_wall,i,j) - slurb_tile%t_indoor(i,j) )

    call calc_heat_diffusion( slurb_tile%t_wall_b_m(:,i,j), slurb_tile%t_wall_b_0(:,i,j), slurb_tile%tt_wall_b(:,i,j),    &
        slurb_tile%c_wall(:,i,j), slurb_tile%conductivity_wall(:,i,j), slurb_tile%t_indoor(i,j) )
    !
    !--    Update longwave radiative fluxes following linearization.
       slurb_tile%rad_lw_net_wall_a(i,j) = slurb_tile%rad_lw_net_wall_a(i,j)                                       &
                                   + slurb_tile%lw_wall_coef(1,i,j) * slurb_tile%t_wall_a_m(nzt_wall,i,j)**4         &
                                   - 4.0_field_r * slurb_tile%lw_wall_coef(1,i,j)                               &
                                      * slurb_tile%t_wall_a_m(nzt_wall,i,j)**3                               &
                                   * ( slurb_tile%t_wall_a_m(nzt_wall,i,j) - slurb_tile%t_wall_a_0(nzt_wall,i,j) )

       slurb_tile%rad_lw_net_wall_b(i,j) = slurb_tile%rad_lw_net_wall_b(i,j)                                       &
                                   + slurb_tile%lw_wall_coef(1,i,j) * slurb_tile%t_wall_b_m(nzt_wall,i,j)**4         &
                                   - 4.0_field_r * slurb_tile%lw_wall_coef(1,i,j)                               &
                                      * slurb_tile%t_wall_b_m(nzt_wall,i,j)**3                               &
                                   * ( slurb_tile%t_wall_b_m(nzt_wall,i,j) - slurb_tile%t_wall_b_0(nzt_wall,i,j) )
    else
    !
    !--    Copy all layers including the surface for wall B.
       slurb_tile%t_wall_b_0(:,i,j) = slurb_tile%t_wall_a_0(:,i,j)
       slurb_tile%tt_wall_b(:,i,j)  = slurb_tile%tt_wall_a(:,i,j)
       slurb_tile%pt_wall_b(i,j)    = slurb_tile%pt_wall_a(i,j)
       slurb_tile%shf_wall_b(i,j)   = slurb_tile%shf_wall_a(i,j)
       slurb_tile%ghf_wall_b(i,j)   = slurb_tile%ghf_wall_a(i,j)

    !
    !--    For longwave radiative flux, we need to add terms for both the t_wall_a and the t_wall_b
    !--    in the longwave balance, thus coefficients 1 and 3 are summed here.
       slurb_tile%rad_lw_net_wall_a(i,j) = slurb_tile%rad_lw_net_wall_a(i,j)                                       &
                                   + ( slurb_tile%lw_wall_coef(1,i,j) + slurb_tile%lw_wall_coef(3,i,j) )           &
                                      * slurb_tile%t_wall_a_m(nzt_wall,i,j)**4                               &
                                   - 4.0_field_r * ( slurb_tile%lw_wall_coef(1,i,j)                             &
                                                + slurb_tile%lw_wall_coef(3,i,j) )                         &
                                   * slurb_tile%t_wall_a_m(nzt_wall,i,j)**3                                  &
                                   * ( slurb_tile%t_wall_a_m(nzt_wall,i,j) - slurb_tile%t_wall_a_0(nzt_wall,i,j) )
       slurb_tile%rad_lw_net_wall_b(i,j) = slurb_tile%rad_lw_net_wall_a(i,j)
    endif

 end subroutine wall_model


    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Models the surface energy balance, SW transmission and subsurface heat diffusion for windows.
    !--------------------------------------------------------------------------------------------------!
 subroutine window_model

    real(field_r) ::  coef_1   !< coefficient A of the prognostic equation
    real(field_r) ::  coef_2   !< coefficient B of the prognostic equation
    real(field_r) ::  f_shf_a  !< factor for the window surface heat flux
    real(field_r) ::  f_shf_b  !< factor for the window surface heat flux


    if ( facade_rah_doe )  then
       f_shf_a = rho_cp / slurb_tile%rah_win_a(i,j)
       if ( slurb_tile%anisotropic_canyon(i,j) )  f_shf_b = rho_cp / slurb_tile%rah_win_b(i,j)
    else
       f_shf_a = rho_cp / slurb_tile%rah_facade(i,j)
       if ( slurb_tile%anisotropic_canyon(i,j) )  f_shf_b = f_shf_a
    endif

    !
    !-- Computation of the prognostic equation similarly to the walls, with exception of added
    !-- shortwave transmission component for surface and subsurface layers. Explanatory comments
    !-- are not repeated from the wall model, comments reflect differences specific to windows.
    if ( slurb_tile%anisotropic_canyon(i,j) )  then
    !
    !--    For windows, some of the incoming shortwave radiation is transmitted through the material.
       coef_1 = slurb_tile%rad_sw_net_win_a(i,j) * slurb_tile%absorption_win(nzt_win,i,j)                          &
                + slurb_tile%rad_lw_net_win_a(i,j)                                                         &
                - 3.0_field_r * slurb_tile%lw_win_coef(1,i,j) * slurb_tile%t_win_a_0(nzt_win,i,j)**4                      &
                + f_shf_a * slurb_tile%t_can_0(i,j)                                                          &
                + slurb_tile%conductivity_win(nzt_win,i,j) * slurb_tile%t_win_a_0(nzt_win+1,i,j)

       coef_2 = -4.0_field_r * slurb_tile%lw_win_coef(1,i,j) * slurb_tile%t_win_a_0(nzt_win,i,j)**3                       &
                + f_shf_a                                                                          &
                + slurb_tile%conductivity_win(nzt_win,i,j)

       call calc_surf_t_p( slurb_tile%t_win_a_m(nzt_win,i,j), slurb_tile%t_win_a_0(nzt_win,i,j),                     &
                           slurb_tile%tt_win_a(nzt_win,i,j), coef_1, coef_2, slurb_tile%c_win(nzt_win,i,j) )

       coef_1 = slurb_tile%rad_sw_net_win_b(i,j) * slurb_tile%absorption_win(nzt_win,i,j)                          &
                + slurb_tile%rad_lw_net_win_b(i,j)                                                         &
                - 3.0_field_r * slurb_tile%lw_win_coef(1,i,j) * slurb_tile%t_win_b_0(nzt_win,i,j)**4                      &
                + f_shf_b * slurb_tile%t_can_0(i,j)                                                          &
                + slurb_tile%conductivity_win(nzt_win,i,j) * slurb_tile%t_win_b_0(nzt_win+1,i,j)

       coef_2 = -4.0_field_r * slurb_tile%lw_win_coef(1,i,j) * slurb_tile%t_win_b_0(nzt_win,i,j)**3                       &
                + f_shf_b                                                                          &
                + slurb_tile%conductivity_win(nzt_win,i,j)

       call calc_surf_t_p( slurb_tile%t_win_b_m(nzt_win,i,j), slurb_tile%t_win_b_0(nzt_win,i,j),                     &
                           slurb_tile%tt_win_b(nzt_win,i,j), coef_1, coef_2, slurb_tile%c_win(nzt_win,i,j) )

    else
       coef_1 = slurb_tile%rad_sw_net_win_a(i,j) * slurb_tile%absorption_win(nzt_win,i,j)                          &
                + slurb_tile%rad_lw_net_win_a(i,j)                                                         &
                - 3.0_field_r * ( slurb_tile%lw_win_coef(1,i,j) + slurb_tile%lw_win_coef(3,i,j) )                       &
                   * slurb_tile%t_win_a_0(nzt_win,i,j)**4                                                    &
                + f_shf_a * slurb_tile%t_can_0(i,j)                                                          &
                + slurb_tile%conductivity_win(nzt_win,i,j) * slurb_tile%t_win_a_0(nzt_win+1,i,j)

       coef_2 = -4.0_field_r * ( slurb_tile%lw_win_coef(1,i,j) + slurb_tile%lw_win_coef(3,i,j) )                        &
                   * slurb_tile%t_win_a_0(nzt_win,i,j)**3                                                    &
                + f_shf_a                                                                          &
                + slurb_tile%conductivity_win(nzt_win,i,j)

       call calc_surf_t_p( slurb_tile%t_win_a_m(nzt_win,i,j), slurb_tile%t_win_a_0(nzt_win,i,j),                     &
                           slurb_tile%tt_win_a(nzt_win,i,j), coef_1, coef_2, slurb_tile%c_win(nzt_win,i,j) )
    endif
    !
    !-- The transmitted shortwave radiation is included also in the prognostic equations for material
    !-- subsurface temperatures.
    slurb_tile%pt_win_a(i,j)  = slurb_tile%t_win_a_0(nzt_win,i,j) * (1 / exnf(k_topo))
    slurb_tile%shf_win_a(i,j) = -f_shf_a * ( slurb_tile%t_can_m(i,j) - slurb_tile%t_win_a_0(nzt_win,i,j) )



    call calc_heat_diffusion( slurb_tile%t_win_a_m(:,i,j), slurb_tile%t_win_a_0(:,i,j),                              &
                              slurb_tile%tt_win_a(:,i,j), slurb_tile%c_win(:,i,j),                                 &
                              slurb_tile%conductivity_win(:,i,j), slurb_tile%t_indoor(i,j),                        &
                              slurb_tile%rad_sw_net_win_a(i,j), slurb_tile%absorption_win(:,i,j) )


    slurb_tile%ghf_win_a(i,j) = slurb_tile%conductivity_win(nzb_win,i,j) *                                         &
                        ( slurb_tile%t_win_a_0(nzb_win,i,j) - slurb_tile%t_indoor(i,j) )

    if ( slurb_tile%anisotropic_canyon(i,j) )  then
       slurb_tile%pt_win_b(i,j)  = slurb_tile%t_win_b_0(nzt_win,i,j) * (1 / exnf(k_topo))
       slurb_tile%shf_win_b(i,j) = -f_shf_b * ( slurb_tile%t_can_m(i,j) - slurb_tile%t_win_b_0(nzt_win,i,j) )


       call calc_heat_diffusion( slurb_tile%t_win_b_m(:,i,j), slurb_tile%t_win_b_0(:,i,j),                           &
                                 slurb_tile%tt_win_b(:,i,j), slurb_tile%c_win(:,i,j),                              &
                                 slurb_tile%conductivity_win(:,i,j), slurb_tile%t_indoor(i,j),                     &
                                 slurb_tile%rad_sw_net_win_b(i,j), slurb_tile%absorption_win(:,i,j) )

       slurb_tile%ghf_win_b(i,j) = slurb_tile%conductivity_win(nzb_win,i,j) *                                      &
                           ( slurb_tile%t_win_b_0(nzb_win,i,j) - slurb_tile%t_indoor(i,j) )

       slurb_tile%rad_lw_net_win_a(i,j) = slurb_tile%rad_lw_net_win_a(i,j)                                         &
                                  + slurb_tile%lw_win_coef(1,i,j) * slurb_tile%t_win_a_m(nzt_win,i,j)**4             &
                                  - 4.0_field_r * slurb_tile%lw_win_coef(1,i,j) * slurb_tile%t_win_a_m(nzt_win,i,j)**3    &
                                     * ( slurb_tile%t_win_a_m(nzt_win,i,j) - slurb_tile%t_win_a_0(nzt_win,i,j) )

       slurb_tile%rad_lw_net_win_b(i,j) = slurb_tile%rad_lw_net_win_b(i,j)                                         &
                                  + slurb_tile%lw_win_coef(1,i,j) * slurb_tile%t_win_b_m(nzt_win,i,j)**4             &
                                  - 4.0_field_r * slurb_tile%lw_win_coef(1,i,j) * slurb_tile%t_win_b_m(nzt_win,i,j)**3    &
                                     * ( slurb_tile%t_win_b_m(nzt_win,i,j) - slurb_tile%t_win_b_0(nzt_win,i,j) )
    else
       slurb_tile%t_win_b_0(:,i,j) = slurb_tile%t_win_a_0(:,i,j)
       slurb_tile%tt_win_b(:,i,j)  = slurb_tile%tt_win_a(:,i,j)
       slurb_tile%pt_win_b(i,j)    = slurb_tile%pt_win_a(i,j)
       slurb_tile%shf_win_b(i,j)   = slurb_tile%shf_win_a(i,j)
       slurb_tile%ghf_win_b(i,j)   = slurb_tile%ghf_win_a(i,j)

       slurb_tile%rad_lw_net_win_a(i,j) = slurb_tile%rad_lw_net_win_a(i,j)                                         &
                                  + ( slurb_tile%lw_win_coef(1,i,j) + slurb_tile%lw_win_coef(3,i,j) )              &
                                     * slurb_tile%t_win_a_m(nzt_win,i,j)**4                                  &
                                  - 4.0_field_r * ( slurb_tile%lw_win_coef(1,i,j) + slurb_tile%lw_win_coef(3,i,j) )     &
                                  * slurb_tile%t_win_a_m(nzt_win,i,j)**3                                     &
                                  * ( slurb_tile%t_win_a_m(nzt_win,i,j) - slurb_tile%t_win_a_0(nzt_win,i,j) )
       slurb_tile%rad_lw_net_win_b(i,j) = slurb_tile%rad_lw_net_win_a(i,j)
    endif

 end subroutine window_model


    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Calculate surface mixing ratio using resistance weighting.
    !--------------------------------------------------------------------------------------------------!
 pure function q_surf( q_s, rah, q_a, f_qsws )

    real(field_r), intent(in) ::  f_qsws  !< factor for the latent heat flux (W m^-2)
    real(field_r), intent(in) ::  q_a     !< mixing ratio of adjacent air
    real(field_r), intent(in) ::  q_s     !< saturation mixing ratio at the surface
    real(field_r), intent(in) ::  rah     !< aerodynamic resistance for heat (and for water vapor)

    real(field_r) ::  q_surf  !< mixing ratio for the surface
    real(field_r) ::  res     !< total surface resistance


    !
    !-- Total surface resistance.
    res = rah / ( rah + ABS( rho_lv / ( f_qsws + 1.0E-20_field_r ) - rah ) )

    !--    Assume equal liquid water content in canyon as in air above.
    q_surf = res * q_s + ( 1.0_field_r - res ) * ( q_a - ql0(i,j,k_atm) )


 end function q_surf
end subroutine slurb_energy_balance_model
end module modslurb_energybalance