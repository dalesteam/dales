!> \file modslurbresistance_stability.f90
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
module modslurb_resistance_stability
    use modprecision, only: field_r
    use modslurbdata
    use modglobal, only: g => grav, cp, kappa => fkar, pi
    use modfields, only: rho_air_zw => rhobf

    contains
!--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Computes the heat and momentum fluxes between the atmosphere and the urban surface.
    !--------------------------------------------------------------------------------------------------!
 subroutine calc_urban_resistances
    use modglobal, only : i1, j1
    implicit none
    integer i, j


    !
    !-- Compute friction velocity and aerodynamic resistance for momentum for the whole urban surface.
    !-- As SLUrb doesn't explicitly compute the momentum flux for each individual surface, and as
    !-- pressure drag needs to be included in the total urban drag, us_urb and rah_urb are computed
    !-- using given roughness length for whole urban fabric (z0_urb, user input). To compute the MOST
    !-- stability corrections, we follow the SURFEX implementation where weighted pt/vpt from canyons
    !-- and roofs is used to represent the pt/vpt at roof level. For urban heat fluxes, aggregated
    !-- values from roofs and canyons are directly used, so rah_urb is not needed.
        if ( moist_physics )  then

        do j=2,j1
            do i=2,i1
                pt_surface(i,j) = slurb_tile%f_bld(i,j)              * slurb_tile%vpt_roof(i,j) +                          &
                                ( 1.0_field_r - slurb_tile%f_bld(i,j) ) * slurb_tile%vpt_can(i,j)
                call calc_rib( slurb_tile%vpt1(i,j), pt_surface(i,j), slurb_tile%rib_urb(i,j), slurb_tile%uv_eff1(i,j), slurb_tile%z_mo(i,j) )
            enddo
        enddo
       

    else

        do j=2,j1
            do i=2,i1
                pt_surface(i,j) = slurb_tile%f_bld(i,j)              * slurb_tile%pt_roof(i,j) +                           &
                                ( 1.0_field_r - slurb_tile%f_bld(i,j) ) * slurb_tile%pt_can(i,j)
                call calc_rib( slurb_tile%pt1(i,j), pt_surface(i,j), slurb_tile%rib_urb(i,j), slurb_tile%uv_eff1(i,j), slurb_tile%z_mo(i,j) )
            enddo
        enddo

    endif

    do j=2,j1
        do i=2,i1
            call calc_ol( ln_z_z0_urb(i,j), ln_z_z0_urb(i,j), slurb_tile%ol_urb(i,j), slurb_tile%rib_urb(i,j), slurb_tile%z0_urb(i,j),       &
                  slurb_tile%z0_urb(i,j), slurb_tile%z_mo(i,j) )
        enddo
    enddo



    do j=2,j1
        do i=2,i1
       slurb_tile%us_urb(i,j) = kappa * slurb_tile%uv_eff1(i,j) /                                                  &
                        ( LOG( slurb_tile%z_mo(i,j) / slurb_tile%z0_urb(i,j) ) -                                   &
                          psi_m( slurb_tile%z_mo(i,j) / slurb_tile%ol_urb(i,j) ) +                                 &
                          psi_m( slurb_tile%z0_urb(i,j) / slurb_tile%ol_urb(i,j) ) )
        enddo
    enddo

    !
    !-- Ensure physical friction velocity (might be needed due to instabilities in e.g. initialization)
    do j=2,j1
        do i=2,i1
       if ( slurb_tile%us_urb(i,j) <= us_min ) slurb_tile%us_urb(i,j) = us_min

       slurb_tile%ram_urb(i,j) = 1.0_field_r / ( kappa * slurb_tile%us_urb(i,j) ) *                                     &
                         ( LOG( slurb_tile%z_mo(i,j) / slurb_tile%z0_urb(i,j) ) -                                  &
                           psi_m( slurb_tile%z_mo(i,j) / slurb_tile%ol_urb(i,j) ) +                                &
                           psi_m( slurb_tile%z0_urb(i,j) / slurb_tile%ol_urb(i,j) ) )

       if ( slurb_tile%ram_urb(i,j) < ram_min )  slurb_tile%ram_urb(i,j) = ram_min
        enddo
    enddo
    !
    !-- For street canyons, effective mixing between canyon half-height and roof height is assumed, thus
    !-- z_mo is used as as reference height when considering atmosphere-street canyon air mixing.
    !-- This is equivalent to mixing of canyon air between the roof top level and the first atm grid
    !-- level. For canyons, roughness length for the whole urban fabric (z0_urb) is used instead of
    !-- local z0m/z0h. This is based on an assumption that turbulence can effectively mix the two air
    !-- masses (canyon air and the atmospheric air). The same assumption is used in TEB/SURFEX.
    !-- Using e.g. the Kanda et al. (2007) parametrization or any other surface parametrization for
    !-- canyon z0h would yield unrealistically low mixing.
    !
    !-- Update z0h for roofs following Kanda et al. (2007) parametrization if enabled.
    if ( roughness_kanda )  then
        do j=2,j1
            do i=2,i1
                slurb_tile%z0h_roof(i,j) = slurb_tile%z0_roof(i,j) * 7.4_field_r *                                            &
                                    EXP( -1.29_field_r * SQRT( SQRT( slurb_tile%z0_roof(i,j) * slurb_tile%us_roof(i,j) /       &
                                                                1.461E-5_field_r) ) )
                ln_z_z0h_roof(i,j) = LOG( slurb_tile%z_mo(i,j) / slurb_tile%z0h_roof(i,j) )
            enddo
        enddo
    endif

    if ( moist_physics )  then
        do j=2,j1
            do i=2, i1
                call calc_rib( slurb_tile%vpt1(i,j), slurb_tile%vpt_roof(i,j), slurb_tile%rib_roof(i,j), slurb_tile%uv_eff1(i,j), slurb_tile%z_mo(i,j) )
                call calc_rib( slurb_tile%vpt1(i,j), slurb_tile%vpt_can(i,j),  slurb_tile%rib_can(i,j),  slurb_tile%uv_eff1(i,j), slurb_tile%z_mo(i,j) )
        enddo
    enddo
    else
        do j=2,j1
            do i=2, i1
                call calc_rib( slurb_tile%pt1(i,j), slurb_tile%pt_roof(i,j), slurb_tile%rib_roof(i,j), slurb_tile%uv_eff1(i,j), slurb_tile%z_mo(i,j) )
                call calc_rib( slurb_tile%pt1(i,j), slurb_tile%pt_can(i,j),  slurb_tile%rib_can(i,j),  slurb_tile%uv_eff1(i,j), slurb_tile%z_mo(i,j) )
            enddo
        enddo
    endif
    do j=2,j1
        do i=2, i1
            ! write(*,*), "now_ol"
            call calc_ol( ln_z_z0_roof(i,j), ln_z_z0h_roof(i,j), slurb_tile%ol_roof(i,j), slurb_tile%rib_roof(i,j), slurb_tile%z0_roof(i,j), &
                        slurb_tile%z0h_roof(i,j), slurb_tile%z_mo(i,j) )
            ! write(*,*), i,j
            call calc_ol( ln_z_z0_urb(i,j), ln_z_z0_urb(i,j), slurb_tile%ol_can(i,j), slurb_tile%rib_can(i,j), slurb_tile%z0_urb(i,j),       &
                        slurb_tile%z0_urb(i,j), slurb_tile%z_mo(i,j) )
            ! write(*,*), "alldone"
        enddo
    enddo

    !
    !-- Compute the local friction velocity for roof and canyon.
    do j=2,j1
      do i=2,i1
       slurb_tile%us_roof(i,j) = kappa * slurb_tile%uv_eff1(i,j) /                                                 &
                         ( LOG( slurb_tile%z_mo(i,j) / slurb_tile%z0_roof(i,j) ) -                                 &
                           psi_m( slurb_tile%z_mo(i,j) / slurb_tile%ol_roof(i,j) ) +                               &
                           psi_m( slurb_tile%z0_roof(i,j) / slurb_tile%ol_roof(i,j) ) )

    !
    !--    For canyons, use urban roughness length (assume the air mixes efficiently
    !--    between the canyon air and atmosphere).
       slurb_tile%us_can(i,j) = kappa * slurb_tile%uv_eff1(i,j) /                                                  &
                        ( LOG( slurb_tile%z_mo(i,j) / slurb_tile%z0_urb(i,j) ) -                                   &
                          psi_m( slurb_tile%z_mo(i,j) / slurb_tile%ol_can(i,j) ) +                                 &
                          psi_m( slurb_tile%z0_urb(i,j) / slurb_tile%ol_can(i,j) ) )

    !
    !--    Ensure physical friction velocity.
       if ( slurb_tile%us_roof(i,j) <= us_min )  slurb_tile%us_roof(i,j) = us_min
      enddo
    enddo

    !
    !-- Compute the aerodynamic resistances for heat.
    do j=2,j1
      do i=2,i1
       slurb_tile%rah_roof(i,j) = 1.0_field_r / ( kappa * slurb_tile%us_roof(i,j) ) *                                   &
                          ( LOG( slurb_tile%z_mo(i,j) / slurb_tile%z0h_roof(i,j) ) -                               &
                            psi_h( slurb_tile%z_mo(i,j) / slurb_tile%ol_roof(i,j) ) +                              &
                            psi_h( slurb_tile%z0h_roof(i,j) / slurb_tile%ol_roof(i,j) ) )

       slurb_tile%rah_can(i,j) = 1.0_field_r / ( kappa * slurb_tile%us_can(i,j) ) *                                     &
                         ( LOG( slurb_tile%z_mo(i,j) / slurb_tile%z0_urb(i,j) ) -                                  &
                           psi_h( slurb_tile%z_mo(i,j) / slurb_tile%ol_can(i,j) ) +                                &
                           psi_h( slurb_tile%z0_urb(i,j) / slurb_tile%ol_can(i,j) ) )

       if ( slurb_tile%rah_roof(i,j) < rah_min )  slurb_tile%rah_roof(i,j) = rah_min
       if ( slurb_tile%rah_roof(i,j) > rah_max )  slurb_tile%rah_roof(i,j) = rah_max
    !
    !--    Use ram_min for canyon air as turbulence is able to mix the air.
       if ( slurb_tile%rah_can(i,j) < ram_min )  slurb_tile%rah_can(i,j) = ram_min
       if ( slurb_tile%rah_can(i,j) > rah_max )  slurb_tile%rah_can(i,j) = rah_max
      enddo
    enddo

 end subroutine calc_urban_resistances



    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Model for the surface resistances within the street canyon.
    !--------------------------------------------------------------------------------------------------!
 subroutine calc_canyon_resistances
    use modglobal, only : i1, j1
    
    implicit none

    integer ::  i       !< loop index x-direction
    integer ::  j       !< loop index y-direction
    integer ::  k_topo  !< k-index of topography
    ! integer ::  m       !< running index of surface tiles


    !
    !-- Calculate logarithms of ratio z/z0.
    !>  TODO: Since the ratios do not change during the simulation, they can be stored once at the
    !>        and stored in surf_slurb, like it is done for the other surface types, too.
    do j=2,j1
      do i=2,i1
       ln_z_z0_road(i,j)  = LOG( slurb_tile%z_mo_can(i,j) / slurb_tile%z0_road(i,j)  )
       ln_z_z0h_road(i,j) = LOG( slurb_tile%z_mo_can(i,j) / slurb_tile%z0h_road(i,j) )
       ln_z_z0_roof(i,j)  = LOG( slurb_tile%z_mo(i,j)     / slurb_tile%z0_roof(i,j)  )
       ln_z_z0h_roof(i,j) = LOG( slurb_tile%z_mo(i,j)     / slurb_tile%z0h_roof(i,j) )
      enddo
    enddo

    !
    !-- Update z0h for roads following Kanda et al. (2007) parametrization if necessary.
    if ( roughness_kanda )  then
    do j=2,j1
      do i=2,i1
            slurb_tile%z0h_road(i,j) = slurb_tile%z0_road(i,j) * 7.4_field_r *                                            &
                                EXP( -1.29_field_r * SQRT( SQRT( slurb_tile%z0_road(i,j) * slurb_tile%us_road(i,j) /       &
                                                            1.461E-5_field_r) ) )
            ln_z_z0h_road(i,j) = LOG( slurb_tile%z_mo_can(i,j) / slurb_tile%z0h_road(i,j) )
             enddo
        enddo
    endif

    !
    !-- Compute the new Obukhov length for road.
    if ( moist_physics )  then
        do j=2,j1
            do i=2,i1
                call calc_rib( slurb_tile%vpt_can(i,j), slurb_tile%vpt_road(i,j), slurb_tile%rib_road(i,j), slurb_tile%uv_eff_can(i,j),        &
                             slurb_tile%z_mo_can(i,j) )
            enddo
        enddo
    else
        do j=2,j1
            do i=2,i1
                call calc_rib( slurb_tile%pt_can(i,j), slurb_tile%pt_road(i,j), slurb_tile%rib_road(i,j), slurb_tile%uv_eff_can(i,j),          &
                               slurb_tile%z_mo_can(i,j) )
            enddo
        enddo
    endif

    do j=2,j1
      do i=2,i1
        call calc_ol( ln_z_z0_road(i,j), ln_z_z0h_road(i,j), slurb_tile%ol_road(i,j), slurb_tile%rib_road(i,j), slurb_tile%z0_road(i,j), &
                    slurb_tile%z0h_road(i,j), slurb_tile%z_mo_can(i,j) )
      enddo
    enddo

    !
    !-- Compute the local friction velocity for roads.
    do j=2,j1
      do i=2,i1
        slurb_tile%us_road(i,j) = kappa * slurb_tile%uv_eff_can(i,j) /                                              &
                            ( LOG( slurb_tile%z_mo_can(i,j) / slurb_tile%z0_road(i,j) ) -                             &
                            psi_m( slurb_tile%z_mo_can(i,j) / slurb_tile%ol_road(i,j) ) +                           &
                            psi_m( slurb_tile%z0_road(i,j) / slurb_tile%ol_road(i,j) ) )
        enddo
    enddo

    !
    !-- The resistance between the street canyon air and facades (walls and windows).
    if ( facade_rah_doe )  then

        do j=2,j1
            do i=2,i1
                ! k_topo = topo_top_ind(j,i,0)  ! ZELFTODO
                k_topo = 1
                slurb_tile%rah_wall_a(i,j) = rah_doe2( k_topo, slurb_tile%t_can_0(i,j), slurb_tile%t_wall_a_0(nzt_wall,i,j),         &
                                                slurb_tile%uv_eff_can(i,j), .TRUE. )
                if ( slurb_tile%rah_wall_a(i,j) < rah_min )  slurb_tile%rah_wall_a(i,j) = rah_min
                if ( slurb_tile%rah_wall_a(i,j) > rah_max )  slurb_tile%rah_wall_a(i,j) = rah_max
                if ( slurb_tile%f_win(i,j) /= 0.0_field_r )  then
                    slurb_tile%rah_win_a(i,j) = rah_doe2( k_topo, slurb_tile%t_can_0(i,j), slurb_tile%t_win_a_0(nzt_win,i,j),         &
                                                slurb_tile%uv_eff_can(i,j), .FALSE. )
                    if ( slurb_tile%rah_win_a(i,j) < rah_min )  slurb_tile%rah_win_a(i,j) = rah_min
                    if ( slurb_tile%rah_win_a(i,j) > rah_max )  slurb_tile%rah_win_a(i,j) = rah_max
                endif

                if ( slurb_tile%anisotropic_canyon(i,j) )  then
                    slurb_tile%rah_wall_b(i,j) = rah_doe2( k_topo, slurb_tile%t_can_0(i,j), slurb_tile%t_wall_b_0(nzt_wall,i,j),      &
                                                    slurb_tile%uv_eff_can(i,j), .TRUE. )
                    if ( slurb_tile%rah_wall_b(i,j) < rah_min )  slurb_tile%rah_wall_b(i,j) = rah_min
                    if ( slurb_tile%rah_wall_b(i,j) > rah_max )  slurb_tile%rah_wall_b(i,j) = rah_max
                    if ( slurb_tile%f_win(i,j) /= 0.0_field_r )  then
                        slurb_tile%rah_win_b(i,j) = rah_doe2( k_topo, slurb_tile%t_can_0(i,j), slurb_tile%t_win_b_0(nzt_win,i,j),      &
                                                    slurb_tile%uv_eff_can(i,j), .FALSE. )
                        if ( slurb_tile%rah_win_b(i,j) < rah_min )  slurb_tile%rah_win_b(i,j) = rah_min
                        if ( slurb_tile%rah_win_b(i,j) > rah_max )  slurb_tile%rah_win_b(i,j) = rah_max
                    endif
                endif
             enddo
        enddo

    ELSEIF ( facade_rah_kray )  then

        do j=2,j1
            do i=2,i1
                ! k_topo = topo_top_ind(j,i,0) ! ZELFTODO
                k_topo = 1
                slurb_tile%rah_facade(i,j) = rah_kray( k_topo, slurb_tile%z0_wall(i,j), slurb_tile%uv_eff_can(i,j) )
                if ( slurb_tile%rah_facade(i,j) < rah_min )  slurb_tile%rah_facade(i,j) = rah_min
                if ( slurb_tile%rah_facade(i,j) > rah_max )  slurb_tile%rah_facade(i,j) = rah_max
            enddo
        enddo

    ELSEIF ( facade_rah_rowley )  then
    !
    !--    Rowley et al. (1930) , Cole and Sturrock (1977)  Mills (1993).
        do j=2,j1
            do i=2,i1
                ! k_topo = topo_top_ind(j,i,0) ! ZELFTODO
                k_topo = 1
                slurb_tile%rah_facade(i,j) = cp * rho_air_zw(k_topo) / ( 11.8_field_r + 4.2_field_r * slurb_tile%uv_eff_can(i,j) )
                if ( slurb_tile%rah_facade(i,j) < rah_min )  slurb_tile%rah_facade(i,j) = rah_min
                if ( slurb_tile%rah_facade(i,j) > rah_max )  slurb_tile%rah_facade(i,j) = rah_max
            enddo
        enddo
    endif

    do j=2,j1
        do i=2,i1
            slurb_tile%rah_road(i,j) = 1.0_field_r / ( kappa * slurb_tile%us_road(i,j) ) *                                    &
                                ( ln_z_z0h_road(i,j) -                                                     &
                                    psi_h( slurb_tile%z_mo_can(i,j) / slurb_tile%ol_road(i,j) ) +                          &
                                    psi_h( slurb_tile%z0h_road(i,j) / slurb_tile%ol_road(i,j) ) )

            if ( slurb_tile%rah_road(i,j) < rah_min )  slurb_tile%rah_road(i,j) = rah_min
            if ( slurb_tile%rah_road(i,j) > rah_max )  slurb_tile%rah_road(i,j) = rah_max
        enddo
    enddo

 end subroutine calc_canyon_resistances



    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Calculate the Obukhov length (L).
    !--------------------------------------------------------------------------------------------------!
 subroutine calc_ol(ln_z_z0, ln_z_z0h, ol, rib, z0, z0h, z_mo )

    IMPLICIT NONE

    real(field_r), intent(in)    ::  ln_z_z0   !< logarithm (z/z0)
    real(field_r), intent(in)    ::  ln_z_z0h  !< logarithm (z/z0h)
    real(field_r), intent(inout) ::  ol        !< Obukhov length
    real(field_r), intent(in)    ::  rib       !< Richardson flux number
    real(field_r), intent(in)    ::  z0        !< rougness length for momentum
    real(field_r), intent(in)    ::  z0h       !< rougness length for scalar quantities
    real(field_r), intent(in)    ::  z_mo      !< constant flux layer height

    integer ::  iter  !< Newton iteration step

    ! LOGICAL ::  convergence_reached  !< convergence switch for vectorization

    real(field_r) ::  f        !< function for Newton iteration: f = Ri - [...]/[...]^2 = 0
    real(field_r) ::  f_d_ol   !< derivative of f
    real(field_r) ::  ol_l     !< lower bound of L for Newton iteration
    real(field_r) ::  ol_m     !< previous value of L for Newton iteration
    real(field_r) ::  ol_prev  !< previous time step value of L
    real(field_r) ::  ol_u     !< upper bound of L for Newton iteration

    ! real ::  ol_prev_vec  !< temporary array required for vectorization

    !
    !-- Calculate the Obukhov length using Newton iteration.

    !
    !--       Store current value in case the Newton iteration fails.
            ol_prev = ol
    !
    !--       Flip the sign of the initial Obukhov length if the stability has changed from stable to
    !--       unstable or vice versa and set it to a moderate value. A moderate value is also chosen,
    !--       if the Obukhov length from the last time step reached the maximum threshold value.
            if ( rib * ol < 0.0_field_r  .OR.  ABS( ol ) == ol_max )  then
                if ( rib > 0.0_field_r )  ol =  100.0_field_r
                if ( rib < 0.0_field_r )  ol = -100.0_field_r
            endif
    !
    !--       Iteration to find Obukhov length.
            iter = 0
            do
                iter = iter + 1
    !
    !--          In case of divergence, use the value of the previous time step.
                if ( iter > 1000 )  then
                ol = ol_prev
                EXIT
                endif

    !
    !--          Calculate step size for central difference.
                ol_m = ol
                ol_l = ol_m - 0.001_field_r * ol_m
                ol_u = ol_m + 0.001_field_r * ol_m

    !
    !--             Calculate f = Ri - [...]/[...]^2 = 0.
                f = rib - ( z_mo / ol_m ) * ( ln_z_z0h - psi_h( z_mo / ol_m )          &
                                                                + psi_h( z0h  / ol_m ) )        &
                                                / ( ln_z_z0  - psi_m( z_mo / ol_m )          &
                                                                + psi_m( z0   / ol_m ) )**2
    !
    !--             Calculate df/dL.
                f_d_ol = ( - ( z_mo / ol_u ) * ( ln_z_z0h - psi_h( z_mo / ol_u )          &
                                                                + psi_h( z0h  / ol_u ) )        &
                                                / ( ln_z_z0  - psi_m( z_mo / ol_u )          &
                                                                + psi_m( z0   / ol_u ) )**2     &
                            + ( z_mo / ol_l ) * ( ln_z_z0h - psi_h( z_mo / ol_l )          &
                                                                + psi_h( z0h  / ol_l ) )        &
                                                / ( ln_z_z0  - psi_m( z_mo / ol_l )          &
                                                                + psi_m( z0   / ol_l ) )**2     &
                            ) / ( ol_u - ol_l )
    !
    !--          Calculate new L.
                ol = ol_m - f / f_d_ol
    !
    !--          Ensure that the bulk Richardson number and the Obukhov length have the same sign and
    !--          ensure convergence. If the sign is not the same, the above calculated Obukhov length
    !--          obviously overshooted to the opposite side, so the next iteration should start with
    !--          a smaller value.
                if ( ol * ol_m < 0.0_field_r )  ol = ol_m * 0.5_field_r
    !
    !--          In the deep neutral zone, set L to the maximum allowed value.
                if ( ABS( ol ) > ol_max )  then
                ol = SIGN( ol_max, ol )
                EXIT
                endif
    !
    !--          Assure that Obukhov length does not become zero.
                if ( ABS( ol ) < ol_min )  then
                ol = SIGN( ol_min, ol )
                EXIT
                endif
    !
    !--          Check for convergence.
                if ( ABS( ( ol - ol_m ) /  ol ) < ol_tol )  EXIT

            enddo

 end subroutine calc_ol



    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Calculate the bulk Richardson number for given surface (z0) temperature.
    !--------------------------------------------------------------------------------------------------!
 subroutine calc_rib( pt1, pt_surface, rib, uvw_abs, z_mo )
    implicit none

    real(field_r), intent(in)  ::  pt1          !< potential temperature at first grid level
    real(field_r), intent(in)  ::  pt_surface   !< skin-surface potential temperature
    real(field_r), intent(out) ::  rib          !< Richardson flux number
    real(field_r), intent(in)  ::  uvw_abs      !< absolute surface-parallel velocity on grid center
    real(field_r), intent(in)  ::  z_mo         !< constant flux layer height

    !-- Evaluate bulk Richardson number.
    rib = g * z_mo * ( pt1 - pt_surface ) / ( uvw_abs**2 * pt1 + 1.0E-20_field_r )

    !
    !-- For the SLUrb model, limit to |rib| < |rib_max| to dampen possible instabilities during
    !-- initialization.
    if ( ABS( rib ) > rib_max )  rib = SIGN( rib_max, rib )

 end subroutine calc_rib



    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Integrated stability function for momentum.
    !--------------------------------------------------------------------------------------------------!
 pure function psi_m( zeta )

    IMPLICIT NONE

    real(field_r), intent(in) ::  zeta   !< Stability parameter z/L

    real(field_r) ::  psi_m  !< Integrated similarity function result
    real(field_r) ::  x      !< dummy variable

    real(field_r), PARAMETER ::  a = 1.0_field_r            !< constant
    real(field_r), PARAMETER ::  b = 0.66666666666_field_r  !< constant
    real(field_r), PARAMETER ::  c = 5.0_field_r            !< constant
    real(field_r), PARAMETER ::  d = 0.35_field_r           !< constant
    real(field_r), PARAMETER ::  c_d_d = c / d         !< constant
    real(field_r), PARAMETER ::  bc_d_d = b * c / d    !< constant


    if ( zeta < 0.0_field_r )  then
       x = SQRT( SQRT( 1.0_field_r  - 16.0_field_r * zeta ) )
       psi_m = pi * 0.5_field_r - 2.0_field_r * ATAN( x ) + LOG( ( 1.0_field_r + x )**2                           &
               * ( 1.0_field_r + x**2 ) * 0.125_field_r )
    else

       psi_m = - b * ( zeta - c_d_d ) * EXP( -d * zeta ) - a * zeta - bc_d_d
    !
    !--    Old version for stable conditions (only valid for z/L < 0.5) psi_m = - 5.0_field_r * zeta

    endif

 end function psi_m


    !--------------------------------------------------------------------------------------------------!
    ! Description:
    !------------
    !> Integrated stability function for heat and moisture.
    !--------------------------------------------------------------------------------------------------!
 pure function psi_h( zeta )

    IMPLICIT NONE

    real(field_r), intent(in) ::  zeta   !< stability parameter z/L

    real(field_r) ::  psi_h  !< integrated similarity function result
    real(field_r) ::  x      !< dummy variable

    real(field_r), PARAMETER ::  a = 1.0_field_r            !< constant
    real(field_r), PARAMETER ::  b = 0.66666666666_field_r  !< constant
    real(field_r), PARAMETER ::  c = 5.0_field_r            !< constant
    real(field_r), PARAMETER ::  d = 0.35_field_r           !< constant
    real(field_r), PARAMETER ::  c_d_d = c / d         !< constant
    real(field_r), PARAMETER ::  bc_d_d = b * c / d    !< constant


    if ( zeta < 0.0_field_r )  then
       x = SQRT( 1.0_field_r  - 16.0_field_r * zeta )
       psi_h = 2.0_field_r * LOG( (1.0_field_r + x ) / 2.0_field_r )
    else
       psi_h = - b * ( zeta - c_d_d ) * EXP( -d * zeta ) - (1.0_field_r                                 &
               + 0.66666666666_field_r * a * zeta )**1.5_field_r - bc_d_d + 1.0_field_r
    !
    !--    Old version for stable conditions (only valid for z/L < 0.5)
    !--    psi_h = - 5.0_field_r * zeta
    endif

 end function psi_h


    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Calculates stability function for momentum
    !>
    !> @author Hauke Wurps
    !--------------------------------------------------------------------------------------------------!
 pure function phi_m( zeta )

    IMPLICIT NONE

    real(field_r), intent(in) ::  zeta   !< stability parameter z/L

    real(field_r) ::  phi_m  !< value of the function

    real(field_r), PARAMETER ::  a = 16.0_field_r  !< constant
    real(field_r), PARAMETER ::  c = 5.0_field_r   !< constant

    if ( zeta < 0.0_field_r )  then
       phi_m = 1.0_field_r / SQRT( SQRT( 1.0_field_r - a * zeta ) )
    else
       phi_m = 1.0_field_r + c * zeta
    endif

 end function phi_m



    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Compute aerodynamic resistance for heat for vertical surfaces following DOE-2 parametrization,
    !> which takes natural convection into account. Average of leeward and windward sides.
    !> Source: EnegyPlus 23.2.0 Engineering Reference p.68.
    !--------------------------------------------------------------------------------------------------!
 pure function rah_doe2( k_topo, t_air, t_surf, u_eff, rough )

    LOGICAL, intent(in) ::  rough  !< flag for rough surface, true for walls, false for windows

    integer, intent(in) ::  k_topo  !< k-index of topography

    real(field_r), intent(in) ::  t_air   !< temperature of adjacent air
    real(field_r), intent(in) ::  t_surf  !< surface temperature
    real(field_r), intent(in) ::  u_eff   !< effective wind speed

    real(field_r), PARAMETER ::  r_f = 1.52_field_r  !< surface roughness multiplier

    real(field_r) ::  chtcn       !< convective heat transfer coefficient for natural convection
    real(field_r) ::  chtcs       !< convective heat transfer coefficient for smooth surface
    real(field_r) ::  chtcs_lee   !< convective heat transfer coefficient for smooth surface (leeward)
    real(field_r) ::  chtcs_wind  !< convective heat transfer coefficient for smooth surface (windward)
    real(field_r) ::  rah_doe2    !< resulting resistance


    chtcn = 1.31_field_r * ABS( t_air - t_surf )**0.33333_field_r

    chtcs_lee  = SQRT( chtcn**2 + ( 2.86_field_r * u_eff**0.617_field_r )**2 )
    chtcs_wind = SQRT( chtcn**2 + ( 2.38_field_r * u_eff**0.89_field_r  )**2 )

    chtcs = 0.5 * ( chtcs_lee + chtcs_wind )

    if ( rough )  then
       rah_doe2 = cp * rho_air_zw(k_topo) / ( chtcn + r_f * ( chtcs - chtcn ) )
    else
       rah_doe2 = cp * rho_air_zw(k_topo) / chtcs
    endif

 end function rah_doe2


    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Compute aerodynamic resistance for heat for vertical surfaces following
    !> Krayenhoff & Voogt (2007).
    !--------------------------------------------------------------------------------------------------!
 pure function rah_kray( k_topo, z0, u_eff )

    integer, intent(in) ::  k_topo  !< k-index of topography

    real(field_r), intent(in) ::  u_eff  !< effective wind speed
    real(field_r), intent(in) ::  z0     !< roughness length for momentum

    real(field_r) ::  kray_coeff  !< denominator for the parametrization
    real(field_r) ::  rah_kray    !< resulting resistance


    !
    !-- Compute denominator first, ensuring it is a positive number.
    kray_coeff = MAX( z0 * 1000.0_field_r * ( 11.8_field_r + 4.2_field_r * u_eff ) - 4.0_field_r, 1.0E-3_field_r )

    rah_kray = cp * rho_air_zw(k_topo) / kray_coeff

 end function
end module modslurb_resistance_stability