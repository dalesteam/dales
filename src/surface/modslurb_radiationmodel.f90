!> \file modslurb_radiationmodel.f90
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
module modslurb_radiationmodel
    use modprecision, only: field_r
    use modglobal, only: pi
    use modslurbdata
    real(field_r) ::  azimuth        !< solar azimuth angle
    real(field_r) ::  tan_zenith     !< tangent of the solar zenith angle
    real ::  cos_zenith              !< cosine of solar zenith angle
    real ::  zenith                  !< solar zenith angle
    real :: sun_dir_lon, sun_dir_lat
    contains

 !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    ! Shortwave and longwave radiation parametrisations of the model.
    !--------------------------------------------------------------------------------------------------!
 subroutine slurb_radiation_model
    use modglobal, only : i1,j1,xtime,rtimee,xday,xlat,xlon
    use modraddata, only : zenith_lon_lat
    integer ::  i            !< loop index
    integer ::  j            !< loop index
    integer ::  k_topo       !< k index of topography top
    integer ::  k_atm        !< k index of the first atmospheric level





      !
      !-- Compute the solar zenith and azimuth angles for the current time and location. These are needed
      !-- for the shortwave radiation calculations
    call zenith_lon_lat(xtime*3600_field_r + rtimee, xday, xlat, xlon, cos_zenith, sun_dir_lon, sun_dir_lat)
    azimuth = ATAN2( sun_dir_lon, sun_dir_lat )
    zenith = ACOS( cos_zenith )
    !
    !-- Split the incoming SW radiation into direct and diffuse parts.
    !-- Direct-diffuse SW split is quite weirdly done in the radiation mod if radiation
    !-- interactions are enabled. However, we do need it here even without interactions.
   do j=2,j1
      do i=2,i1

       k_topo = 1
       k_atm = 1

    !
    !--    Update SLUrb internal radiative fluxes based on the new surface temperatures
    !--    Compute the internal longwave radiation interactions at every timestep.
       call calc_rad_lw

    !
    !--    Compute the SW radiation fluxesd.
    !--    Do this only if the radiation model has updated SW fluxes at previous timestep,
    !--    as otherwise the computation would just yield the same fluxes.
    !    if ( radiation_called  .OR.  first_call )  call calc_rad_sw TODOSELF
       call calc_rad_sw
    enddo
   enddo

    !
    !-- Private functions and subroutines of slurb_radiation_model.
    contains


    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Computes the LW radiative fluxes and their differentials for the time step.
    !--------------------------------------------------------------------------------------------------!
 subroutine calc_rad_lw
   use modglobal, only : boltz
   use modraddata, only : lwd
    real(field_r) ::  t_rad_sky  !< Radiative temperature of the sky


    !
    !-- Compute the effective radiative temperature of the incoming LW radiation.
    slurb_tile%rad_lw_in_urb(i,j) = abs(lwd(i,j,1)) !TODO check which level to use, and check which radiation model is allowed.
    ! lwd is positive in DALES
    t_rad_sky = SQRT( SQRT( slurb_tile%rad_lw_in_urb(i,j) / boltz ) )

    !
    !-- Computation of net LW fluxes based on Lemonsu et al. 2012 Eqs. (1-3) (+ windows).
    !-- Note that these are NOT YET the final net longwave fluxes for the surfaces, as the term
    !-- dependent on the surface's own surface temperature (coef=1) is omitted at this stage.
    !-- This term is added after computing the prognostic equation for the surface temperature,
    !-- as it is included in the prognostic equations in an linearized form.
    slurb_tile%rad_lw_net_roof(i,j) = slurb_tile%lw_roof_coef(2,i,j) * slurb_tile%rad_lw_in_urb(i,j)


    slurb_tile%rad_lw_net_road(i,j) = slurb_tile%lw_road_coef(2,i,j) * slurb_tile%rad_lw_in_urb(i,j) +                     &
                              slurb_tile%lw_road_coef(3,i,j) * slurb_tile%t_wall_a_0(nzt_wall,i,j)**4 +              &
                              slurb_tile%lw_road_coef(3,i,j) * slurb_tile%t_wall_b_0(nzt_wall,i,j)**4 +              &
                              slurb_tile%lw_road_coef(4,i,j) * slurb_tile%t_win_a_0(nzt_win,i,j)**4 +                &
                              slurb_tile%lw_road_coef(4,i,j) * slurb_tile%t_win_b_0(nzt_win,i,j)**4

    !
    !-- The term dependent on t_wall_b is omitted at this stage, as for isotropic canyons the mean wall
    !-- temperature is used, including both wall A and B interactions. Thus, the terms for both
    !-- t_wall_a and t_wall_b have to be included in linearization. For anisotropic canyons there is
    !-- no direct dependence, so it can be directly added (see below).
    slurb_tile%rad_lw_net_wall_a(i,j) = slurb_tile%lw_wall_coef(2,i,j) * slurb_tile%rad_lw_in_urb(i,j) +                   &
                                slurb_tile%lw_wall_coef(4,i,j) * slurb_tile%t_win_a_0(nzt_win,i,j)**4 +              &
                                slurb_tile%lw_wall_coef(5,i,j) * slurb_tile%t_win_b_0(nzt_win,i,j)**4 +              &
                                slurb_tile%lw_wall_coef(6,i,j) * slurb_tile%t_road_0(nzt_road,i,j)**4

    if ( slurb_tile%f_win(i,j) > 0.0_field_r )  then
       slurb_tile%rad_lw_net_win_a(i,j) = slurb_tile%lw_win_coef(2,i,j) * slurb_tile%rad_lw_in_urb(i,j) +                  &
                                  slurb_tile%lw_win_coef(4,i,j) * slurb_tile%t_wall_a_0(nzt_wall,i,j)**4 +           &
                                  slurb_tile%lw_win_coef(5,i,j) * slurb_tile%t_wall_b_0(nzt_wall,i,j)**4 +           &
                                  slurb_tile%lw_win_coef(6,i,j) * slurb_tile%t_road_0(nzt_road,i,j)**4
    endif

    !
    !-- Inverse for facade B, if anisotropic canyons are used. If not, copy.
    if ( slurb_tile%anisotropic_canyon(i,j) )  then
    !
    !--    In case of anisotropic canyons, t_wall_b doesn't have dependency on t_wall_a in the
    !--    prognostic equation, and thus it's contribution to longwave balance can be directly added
    !--    to the net longwave radiation before prognostic equations. Vice versa for t_wall_b.
       slurb_tile%rad_lw_net_wall_a(i,j) = slurb_tile%rad_lw_net_wall_a(i,j) +                                     &
                                   slurb_tile%lw_wall_coef(3,i,j) * slurb_tile%t_wall_b_0(nzt_wall,i,j)**4

    !
    !--    Note that for wall (and window) B the coefficients 4 and 5 are also swapped.
       slurb_tile%rad_lw_net_wall_b(i,j) = slurb_tile%lw_wall_coef(2,i,j) * slurb_tile%rad_lw_in_urb(i,j) +                &
                                   slurb_tile%lw_wall_coef(3,i,j) * slurb_tile%t_wall_a_0(nzt_wall,i,j)**4 +         &
                                   slurb_tile%lw_wall_coef(4,i,j) * slurb_tile%t_win_b_0(nzt_win,i,j)**4 +           &
                                   slurb_tile%lw_wall_coef(5,i,j) * slurb_tile%t_win_a_0(nzt_win,i,j)**4 +           &
                                   slurb_tile%lw_wall_coef(6,i,j) * slurb_tile%t_road_0(nzt_road,i,j)**4

       if ( slurb_tile%f_win(i,j) > 0.0_field_r )  then
          slurb_tile%rad_lw_net_win_a(i,j) = slurb_tile%rad_lw_net_win_a(i,j) +                                    &
                                     slurb_tile%lw_win_coef(3,i,j) * slurb_tile%t_win_b_0(nzt_win,i,j)**4

          slurb_tile%rad_lw_net_win_b(i,j) = slurb_tile%lw_win_coef(2,i,j) * slurb_tile%rad_lw_in_urb(i,j) +               &
                                     slurb_tile%lw_win_coef(3,i,j) * slurb_tile%t_win_a_0(nzt_win,i,j)**4 +          &
                                     slurb_tile%lw_win_coef(4,i,j) * slurb_tile%t_wall_b_0(nzt_wall,i,j)**4 +        &
                                     slurb_tile%lw_win_coef(5,i,j) * slurb_tile%t_wall_a_0(nzt_wall,i,j)**4 +        &
                                     slurb_tile%lw_win_coef(6,i,j) * slurb_tile%t_road_0(nzt_road,i,j)**4
       endif
    else
       slurb_tile%rad_lw_net_wall_b(i,j) = slurb_tile%rad_lw_net_wall_a(i,j)
       slurb_tile%rad_lw_net_win_b(i,j)  = slurb_tile%rad_lw_net_win_a(i,j)
    endif

 end subroutine calc_rad_lw


    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Computes the SW radiative fluxes for the time step.
    !--------------------------------------------------------------------------------------------------!
 subroutine calc_rad_sw
   use modraddata, only : swdir, swdif

    real(field_r) ::  rad_sw_diff_road      !< incoming diffuse shortwave radiation on road
    real(field_r) ::  rad_sw_diff_wall_a    !< incoming diffuse shortwave radiation on wall A
    real(field_r) ::  rad_sw_diff_wall_b    !< incoming diffuse shortwave radiation on wall B
    real(field_r) ::  rad_sw_dir_road       !< incoming direct shortwave radiation on road
    real(field_r) ::  rad_sw_dir_wall_a     !< incoming direct shortwave radiation on wall A
    real(field_r) ::  rad_sw_dir_wall_b     !< incoming direct shortwave radiation on wall B
    real(field_r) ::  rad_sw_ref_nomin      !< nominator of the sum of reflections at infinity.
    real(field_r) ::  rad_sw_wall_modifier  !< modifier term for anisotropic walls
    real(field_r) ::  theta0                !< critical canyon orientation for road illumination
    real(field_r) ::  w_inf                 !< mean wall reflection at infinity



    !
    !-- Check if there is any shortwave radiation to take care of in the first place.
    if ( .NOT. ( cos_zenith > tiny(cos_zenith) ) )  then
       slurb_tile%rad_sw_in_urb(i,j)     = 0.0_field_r
       slurb_tile%rad_sw_net_urb(i,j)    = 0.0_field_r
       slurb_tile%rad_sw_net_roof(i,j)   = 0.0_field_r
       slurb_tile%rad_sw_net_road(i,j)   = 0.0_field_r
       slurb_tile%rad_sw_net_wall_a(i,j) = 0.0_field_r
       slurb_tile%rad_sw_net_wall_b(i,j) = 0.0_field_r
       slurb_tile%albedo_urb(i,j)        = 0.1_field_r
       RETURN
    endif

    ! whatever radiation model we use, shortwave DOWN will always be positive, so we ensure that by taking absolute value.
    slurb_tile%rad_sw_in_urb(i,j) = abs(swdir(i,j,1)) + abs(swdif(i,j,1))

    !
    !-- Compute the net shortwave radiation for roofs, which is the simplest case.
    slurb_tile%rad_sw_net_roof(i,j) = ( 1.0_field_r - slurb_tile%albedo_roof(i,j) ) * slurb_tile%rad_sw_in_urb(i,j)

    !
    !-- Next, compute then et shortwave radiation within the street canyon. This is quite complex,
    !-- including the effect of shading and within-canyon reflections. See Lemonsu et al. (2012)
    !-- for reference.

    !
    !-- Calculate tangent of the zenith angle, with limiters and safety margins applied to prevent
    !-- floating point overflows and division by zero. Shouldn't affect the physics too much.
    if ( ABS( 0.5_field_r * pi - zenith ) < 1.0E-6_field_r )  then
       if ( 0.5_field_r * pi - zenith >  0.0_field_r )  tan_zenith = TAN( 0.5_field_r * pi - 1.0E-6_field_r )
       if ( 0.5_field_r * pi - zenith <= 0.0_field_r )  tan_zenith = TAN( 0.5_field_r * pi + 1.0E-6_field_r )
    ELSEIF ( ABS( zenith ) < 1.0E-6_field_r )  then
       tan_zenith = SIGN(1.0, zenith) * TAN( 1.0E-6_field_r )
    else
       tan_zenith = TAN( zenith )
    endif

    !
    !-- Direct SW radiation received by the walls (and windows), the road and vegetation.
    if ( slurb_tile%anisotropic_canyon(i,j) )  then
    !
    !--    Lemonsu et al. (2012) Eq. (A1)
    !--    @note There is an error in this equation in the article. It should be that
    !--    the direct radiation on road should decrease when difference between the sun azimuth
    !--    angles increase, not vice versa.
       rad_sw_dir_road = abs(swdir(i,j,1)) * MAX( 0.0_field_r, 1.0_field_r - slurb_tile%hw_can(i,j) *               &
                         tan_zenith *  SIN( ABS( azimuth - slurb_tile%theta_can(i,j) ) ) )

    !
    !--    Lemonsu et al. (2012) Eqs. (A2-A4)
       rad_sw_dir_wall_a = ( abs(swdir(i,j,1)) - rad_sw_dir_road ) * 0.5_field_r / slurb_tile%hw_can(i,j)

       if ( SIN( azimuth - slurb_tile%theta_can(i,j) ) > 0.0_field_r )  then
          rad_sw_dir_wall_a = 2.0_field_r * rad_sw_dir_wall_a
          rad_sw_dir_wall_b = 0.0_field_r
       else
          rad_sw_dir_wall_b = 2.0_field_r * rad_sw_dir_wall_a
          rad_sw_dir_wall_a = 0.0_field_r
       endif

    else
    !
    !--    Revert to the anisotropic integrated solution by Masson (2000).
    !
    !--    Calculate the critical canyon orientation theta0 for anisotropic street canyons.
       theta0 = ASIN( MIN( 1.0_field_r / ( tan_zenith * slurb_tile%hw_can(i,j) ), 1.0_field_r ) )

    !
    !--    Masson (2000) Eqs. (13-15)
       rad_sw_dir_road = abs(swdir(i,j,1)) * ( 2.0_field_r * theta0 / pi -                             &
                         2.0_field_r * tan_zenith / pi * slurb_tile%hw_can(i,j) * ( 1.0_field_r - COS( theta0 ) ) )

       rad_sw_dir_wall_a = ( abs(swdir(i,j,1)) - rad_sw_dir_road ) * 0.5_field_r / slurb_tile%hw_can(i,j)

       rad_sw_dir_wall_b = rad_sw_dir_wall_a

   endif

    !
    !-- Diffuse (from sky) solar radiation received by the surfaces.
    rad_sw_diff_road   = abs(swdif(i,j,1)) * slurb_tile%svf_road(i,j)
    rad_sw_diff_wall_a = abs(swdif(i,j,1)) * slurb_tile%svf_wall(i,j)
    rad_sw_diff_wall_b = rad_sw_diff_wall_a

    !
    !-- Canyon internal scattering based on both Masson (2000) Eqs. (16-20) and
    !-- Lemonsu et al. (2012) Appendix A2. This has been modified to include windows: the weighted
    !-- average reflection from walls and windows is taken into account by using weighted average
    !-- albedo. The wall and window surfaces are assumed to be uniformly distributed.

    !
    !-- Nominator of the sum of reflections at infinity.
    rad_sw_ref_nomin = slurb_tile%albedo_wall_win(i,j) * ( rad_sw_dir_wall_a + rad_sw_diff_wall_a +        &
                       rad_sw_dir_wall_b + rad_sw_diff_wall_b ) / 2.0_field_r +                         &
                       slurb_tile%albedo_wall_win(i,j) * slurb_tile%svf_wall(i,j) * slurb_tile%albedo_road(i,j) *          &
                       rad_sw_dir_road

    !
    !-- Sum of refelctions at infinity.
    w_inf = rad_sw_ref_nomin / slurb_tile%sw_ref_denom(i,j)

    !
    !-- Total solar radiation absorbed after infinite reflections.
    slurb_tile%rad_sw_in_road(i,j) = rad_sw_dir_road + rad_sw_diff_road +                                  &
                             ( 1.0_field_r - slurb_tile%svf_road(i,j) ) * w_inf
    slurb_tile%rad_sw_net_road(i,j) = ( 1.0_field_r - slurb_tile%albedo_road(i,j) ) * slurb_tile%rad_sw_in_road(i,j)

    slurb_tile%rad_sw_net_wall_a(i,j) = ( 1.0_field_r - slurb_tile%albedo_wall(i,j) ) *                                 &
                                ( 0.5_field_r * ( rad_sw_dir_wall_a + rad_sw_diff_wall_a +              &
                                             rad_sw_dir_wall_b + rad_sw_diff_wall_b )              &
                                + slurb_tile%albedo_road(i,j) * slurb_tile%svf_wall(i,j) *                         &
                                  ( rad_sw_dir_road + rad_sw_diff_road )                           &
                                + slurb_tile%albedo_road(i,j) * slurb_tile%svf_wall(i,j) *                         &
                                  ( 1.0_field_r - slurb_tile%svf_road(i,j) ) * w_inf                            &
                                + ( 1.0_field_r - 2.0_field_r * slurb_tile%svf_wall(i,j) ) * w_inf                   &
                                )

    slurb_tile%rad_sw_net_wall_b(i,j) = slurb_tile%rad_sw_net_wall_a(i,j)

    if ( slurb_tile%f_win(i,j) /= 0.0_field_r  )  then
       slurb_tile%rad_sw_in_win_a(i,j) =   0.5_field_r * ( rad_sw_dir_wall_a + rad_sw_diff_wall_a               &
                                            + rad_sw_dir_wall_b + rad_sw_diff_wall_b )             &
                                 + slurb_tile%albedo_road(i,j) * slurb_tile%svf_wall(i,j) *                        &
                                   ( rad_sw_dir_road + rad_sw_diff_road )                          &
                                 + slurb_tile%albedo_road(i,j) * slurb_tile%svf_wall(i,j) *                        &
                                   ( 1.0_field_r - slurb_tile%svf_road(i,j) ) * w_inf                           &
                                 + ( 1.0_field_r - 2.0_field_r * slurb_tile%svf_wall(i,j) ) * w_inf

       slurb_tile%rad_sw_net_win_a(i,j) = ( 1.0_field_r - slurb_tile%albedo_win(i,j) ) * slurb_tile%rad_sw_in_win_a(i,j)

       slurb_tile%rad_sw_in_win_b(i,j)  = slurb_tile%rad_sw_in_win_a(i,j)
       slurb_tile%rad_sw_net_win_b(i,j) = slurb_tile%rad_sw_net_win_a(i,j)
    endif

    !
    !-- Modification of reflected solar radiation for anisotropic street canyons.
    if ( slurb_tile%anisotropic_canyon(i,j) )  then
       rad_sw_wall_modifier = ( 1.0_field_r + slurb_tile%albedo_wall_win(i,j) *                                 &
                                ( 1.0_field_r - 2.0_field_r * slurb_tile%svf_wall(i,j) ) /                           &
                                ( 1.0_field_r + slurb_tile%albedo_wall_win(i,j) *                               &
                                  ( 1.0_field_r - 2.0_field_r * slurb_tile%svf_wall(i,j) ) )                         &
                              ) *                                                                  &
                              0.5_field_r * ( ( rad_sw_dir_wall_a + rad_sw_diff_wall_a )                &
                                       - ( rad_sw_dir_wall_b + rad_sw_diff_wall_b ) )

       slurb_tile%rad_sw_net_wall_a(i,j) = slurb_tile%rad_sw_net_wall_a(i,j) +                                     &
                                   ( 1.0_field_r - slurb_tile%albedo_wall(i,j) ) * rad_sw_wall_modifier

       slurb_tile%rad_sw_net_wall_b(i,j) = slurb_tile%rad_sw_net_wall_b(i,j) -                                     &
                                   ( 1.0_field_r - slurb_tile%albedo_wall(i,j) ) * rad_sw_wall_modifier

       if ( slurb_tile%f_win(i,j) /= 0.0_field_r )  then
          slurb_tile%rad_sw_in_win_a(i,j)  = slurb_tile%rad_sw_in_win_a(i,j) + rad_sw_wall_modifier
          slurb_tile%rad_sw_net_win_a(i,j) = slurb_tile%rad_sw_in_win_a(i,j) * ( 1.0_field_r - slurb_tile%albedo_win(i,j) )
          slurb_tile%rad_sw_in_win_b(i,j)  = slurb_tile%rad_sw_in_win_b(i,j) - rad_sw_wall_modifier
          slurb_tile%rad_sw_net_win_b(i,j) = slurb_tile%rad_sw_in_win_b(i,j) * ( 1.0_field_r - slurb_tile%albedo_win(i,j) )
       endif
    endif

    !
    !-- The upward shortwave radiation is computed as residual of absorbed radiation per uniturban
    !-- area. Aggregated effective albedo of urban surface is computed so that the raditaiton models end
    !-- up with the same figure for outgoing shortwave radiation.
    slurb_tile%rad_sw_out_urb(i,j) = slurb_tile%rad_sw_in_urb(i,j) -                                               &
                             ( ( 1.0_field_r - slurb_tile%f_bld(i,j) ) *                                        &
                               ( slurb_tile%hw_can(i,j) * ( ( 1.0_field_r - slurb_tile%f_win(i,j) ) *                   &
                                       ( slurb_tile%rad_sw_net_wall_a(i,j) + slurb_tile%rad_sw_net_wall_b(i,j) )   &
                                       + slurb_tile%f_win(i,j) *                                           &
                                       ( slurb_tile%rad_sw_net_win_a(i,j)  + slurb_tile%rad_sw_net_win_b(i,j)  )   &
                                                  )                                                &
                               + slurb_tile%rad_sw_net_road(i,j)                                           &
                               )                                                                   &
                             + slurb_tile%f_bld(i,j) * slurb_tile%rad_sw_net_roof(i,j)                             &
                             )

    !
    !-- Compute the net SW flux for diagnostics and output.
    slurb_tile%rad_sw_net_urb(i,j) = slurb_tile%rad_sw_in_urb(i,j) - slurb_tile%rad_sw_out_urb(i,j)

    !
    !-- Save effective albedo for the radiation model.
    if (slurb_tile%rad_sw_in_urb(i,j) /= 0.0) then
     slurb_tile%albedo_urb(i,j) = slurb_tile%rad_sw_out_urb(i,j) / slurb_tile%rad_sw_in_urb(i,j)
    endif

 end subroutine calc_rad_sw

 end subroutine slurb_radiation_model
end module modslurb_radiationmodel