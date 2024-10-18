!
! Copyright (c) 2020-2020 Wageningen University and Research (WUR)
!
! This file is part of DALES
!
! DALES is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with DALES.  If not, see <http://www.gnu.org/licenses/>.
!

module modlsm
    use netcdf
    use modprecision, only : field_r
    implicit none

    public :: initlsm, lsm, exitlsm, init_lsm_tiles, lags, llsm

    logical :: llsm            ! On/off switch LSM
    logical :: lfreedrainage   ! Free drainage bottom BC for soil moisture
    logical :: lags            ! Switch for A-Gs scheme

    ! Interpolation types soil from full to half level
    integer :: iinterp_t, iinterp_theta
    integer, parameter :: iinterp_amean = 1  ! val = 0.5*(v1+v2)
    integer, parameter :: iinterp_gmean = 2  ! val = sqrt(v1*v2)
    integer, parameter :: iinterp_hmean = 3  ! val = ((dz1+dz2)*v1*v2)/(dz1*v2+dz2*v1)
    integer, parameter :: iinterp_max   = 4  ! val = max(a,b)

    ! Soil grid
    integer :: kmax_soil = -1
    real :: z_size_soil = -1
    real, allocatable :: z_soil(:), zh_soil(:)
    real, allocatable :: dz_soil(:), dzh_soil(:)
    real, allocatable :: dzi_soil(:), dzhi_soil(:)

    ! Soil index in `van_genuchten_parameters.nc` lookup table:
    integer, allocatable :: soil_index(:,:,:)

    ! Source term in soil water budget
    real, allocatable :: phiw_source(:,:,:)

    ! Precipitation, interception et al.
    real, allocatable :: throughfall(:,:)
    real, allocatable :: interception(:,:)
    real, allocatable :: wl_max(:,:)

    ! Dependency factor canopy resistance on VPD (high veg only)
    real, allocatable :: gD(:,:)

    ! Reduction functions canopy resistance
    real, allocatable :: f1(:,:), f2b(:,:)

    ! Random
    real, allocatable :: du_tot(:,:), thv_1(:,:), land_frac(:,:), cveg(:,:)

    ! A-Gs
    real, allocatable :: an_co2(:,:), resp_co2(:,:)
    integer :: co2_index = -1

    ! Data structure for sub-grid tiles
    type T_lsm_tile
    ! Fixed LU properties
        ! Land use name
        character(len=64) :: luname
        character(len=3)  :: lushort
        ! Check if LU type is vegetation
        logical           :: lveg
        ! Check if LU type is water
        logical           :: laqu
        ! Static properties:
        real, allocatable :: z0m(:,:), z0h(:,:)
        ! Base tile fraction (i.e. without liquid water)
        real, allocatable :: base_frac(:,:)
        ! Conductivity skin layer:
        real, allocatable :: lambda_stable(:,:), lambda_unstable(:,:)
        ! Vegetation properties:
        real, allocatable :: lai(:,:)
        ! Surface properties:
        real, allocatable :: rs_min(:,:)
        ! Root fraction parameters
        real, allocatable :: a_r(:,:), b_r(:,:)

    ! Dynamic land surface properties
        ! Dynamic tile fraction:
        real, allocatable :: frac(:,:)
        ! Monin-obukhov / surface layer:
        real, allocatable :: obuk(:,:), ustar(:,:), ra(:,:)
        ! Surface fluxes:
        real, allocatable :: H(:,:), LE(:,:), G(:,:)
        real, allocatable :: wthl(:,:), wqt(:,:)
        ! Surface (potential) temperature and humidity:
        real, allocatable :: tskin(:,:), thlskin(:,:), qtskin(:,:)
        ! Buoyancy difference surface - atmosphere
        real, allocatable :: db(:,:)
        ! Vegetation properties:
        real, allocatable :: rs(:,:)
        real, allocatable :: f2(:,:), f3(:,:), gD(:,:)
        ! Root fraction in soil
        real, allocatable :: root_frac(:,:,:)
        ! Root fraction weighted mean soil water content
        real, allocatable :: phiw_mean(:,:)

    ! LU dependent deposition parameters
        ! In-canopy resistance parameters
!        real, allocatable :: R_inc_b(:,:), R_inc_h(:,:)
        ! SAI = SAI_a * LAI + SAI_b
        real, allocatable :: SAI_a(:,:), SAI_b(:,:)
        ! Minimum correction factor for stomatal resistance
!        real, allocatable :: fmin(:,:)
        ! Alpha value for light correction of stomatal resistance
!        real, allocatable :: alpha(:,:)
        ! Min, optimum and max emperatures for temperature correction of stomatal resistance
!        real, allocatable :: Tmin(:,:), Topt(:,:), Tmax(:,:)
        ! Maximum leaf stomatal conductance for ozone
!        real, allocatable :: gs_max(:,:)
        ! Minimum and maximum vapour pressure deficit parameters
!        real, allocatable :: vpd_min(:,:), vpd_max(:,:)
        ! Gamma parameter for calculating stomata compensation point
!        real, allocatable :: gamma_stom(:,:)
        ! Gamma correction factor for calculating soil compensation point
!        real, allocatable :: gamma_soil_c_fac(:,:)
        ! Gamma parameter for calculating soil compensation point
!        real, allocatable :: gamma_soil_default(:,:)

    end type T_lsm_tile

    !Tiles for all LU types
    integer :: ilu, nlu, ilu_ws
    type(T_lsm_tile), allocatable :: tile(:)

    ! Land-surface / van Genuchten parameters from NetCDF input table.
    real, allocatable :: &
        theta_res(:), theta_wp(:), theta_fc(:), theta_sat(:), &
        gamma_theta_sat(:), vg_a(:), vg_l(:), vg_n(:)
    ! Derived soil parameters
    real, allocatable :: &
        vg_m(:), &
        lambda_theta_min(:), lambda_theta_max(:), &
        gamma_theta_min(:), gamma_theta_max(:), &
        gamma_t_dry(:), rho_C(:)

contains

subroutine lsm
    implicit none

    if (.not. llsm) return

    ! Calculate dynamic tile fractions,
    ! based on the amount of liquid water on vegetation.
    call calc_tile_fractions

    ! Calculate root fraction weighted mean soil water content.
    do ilu=1,nlu
      if (tile(ilu)%lveg) then
        call calc_theta_mean(tile(ilu))
      end if
    end do

    ! Calculate canopy/soil resistances.
    if (lags) then
        call calc_canopy_resistance_ags
    else
        call calc_canopy_resistance_js
    endif

    ! Calculate aerodynamic resistance (and u*, obuk).
    call calc_stability

    ! Set grid point averaged boundary conditions (thls, qts, gradients, ..)
    call calc_bulk_bcs

    ! Calculate soil tendencies
    ! Calc diffusivity heat:
    call calc_thermal_properties
    ! Solve diffusion equation:
    call integrate_t_soil

    ! Calc diffusivity and conductivity soil moisture:
    call calc_hydraulic_properties
    ! Calculate tendency due to root water extraction
    call calc_root_water_extraction
    ! Solve diffusion equation:
    call integrate_theta_soil

    ! Update liquid water reservoir
    call calc_liquid_reservoir

end subroutine lsm

!
! Calculate dynamic tile fractions, based on liquid water on vegetation
!
subroutine calc_tile_fractions
    use modglobal, only : i1, j1, eps1
    use modsurfdata, only : wl
    implicit none

    integer :: i, j
    real :: c_liq
    real :: base_frac_sum
    real :: sum_frac, sum_basefrac

    ! Water tile fraction is not variable
    do ilu=1,nlu
      if (tile(ilu)%laqu) then
        tile(ilu)%frac(:,:) = tile(ilu)%base_frac(:,:)
      end if
    end do

    do j=2, j1
      do i=2, i1
        c_liq = min(1., wl(i,j)/wl_max(i,j))
        base_frac_sum = 0
        do ilu=1,nlu
          if (tile(ilu)%laqu) then
            cycle
          else if (ilu == ilu_ws) then
            tile(ilu)%frac(i,j) = land_frac(i,j) * c_liq
          else
            tile(ilu)%frac(i,j) = (1.-c_liq) * tile(ilu)%base_frac(i,j)
          end if
          base_frac_sum = base_frac_sum + tile(ilu)%base_frac(i,j)
        end do
      end do
    end do

    do j=2, j1
      do i=2, i1
        sum_basefrac=0
        sum_frac=0
        do ilu=1,nlu
          sum_basefrac = sum_basefrac + tile(ilu)%base_frac(i,j)
          sum_frac = sum_frac + tile(ilu)%frac(i,j)
        end do

        if (sum_frac < 1.-1e-5) then
          print*,'ERROR: sum of LU fractions below 1:',sum_frac
          stop
        endif
        if (sum_basefrac < 1.-1e-5) then
          print*,'ERROR: sum of LU base fractions below 1:',sum_basefrac
          stop
        endif
      end do
    end do
end subroutine calc_tile_fractions

!
! Calculate changes in the liquid water reservoir
!
subroutine calc_liquid_reservoir
    use modglobal,    only : rk3step, rdt, i1, j1, rhow, rlv
    use modsurfdata,  only : wl, wlm
    use modmicrodata, only : imicro, precep
    use modfields,    only : rhof
    implicit none

    integer :: i, j
    real :: tend, rk3coef, rainrate, wl_tend_max, wl_tend_min
    real :: wl_tend_liq, wl_tend_dew, wl_tend_precip, wl_tend_sum, wl_tend_lim

    real, parameter :: intercept_eff = 0.5
    real, parameter :: to_ms  = 1./(rhow*rlv)

    rk3coef = rdt / (4. - dble(rk3step))
    if(rk3step == 1) wlm(:,:) = wl(:,:)

    do j=2, j1
        do i=2, i1
            wl_tend_dew = 0
            wl_tend_liq = 0

            ! Max and min possible tendencies
            wl_tend_min = -wlm(i,j) / rk3coef
            wl_tend_max = (wl_max(i,j) - wlm(i,j)) / rk3coef

            do ilu=1,nlu
              if (tile(ilu)%laqu) then
                cycle
              end if

              ! Tendency due to evaporation from liquid water reservoir/tile.
              !if (trim(tile(ilu)%lushort) == 'ws') then
              if (ilu == ilu_ws) then
                wl_tend_liq = wl_tend_liq -max(0., tile(ilu)%frac(i,j) * tile(ilu)%LE(i,j) * to_ms)
              end if

              ! Tendency due to dewfall into vegetation/soil/liquid water tiles
              wl_tend_dew = wl_tend_dew &
                -( min(0., tile(ilu)%frac(i,j) * tile(ilu)%LE(i,j) * to_ms) )
            end do

            ! Tendency due to interception of precipitation by vegetation
            if (imicro == 0) then
                rainrate = 0.
            else
               !rainrate = -sed_qr(i,j,1)/rhow
                rainrate = precep(i,j,1)*rhof(1)/rhow ! positive for a downward rain flux
            end if

            wl_tend_precip = intercept_eff * cveg(i,j) * rainrate

            ! Total and limited tendencies
            wl_tend_sum = wl_tend_liq + wl_tend_dew + wl_tend_precip
            wl_tend_lim = min(wl_tend_max, max(wl_tend_min,  wl_tend_sum))

            ! Diagnose interception and throughfall. Upward flux is positive.
            throughfall(i,j) = &
                -(1.-cveg(i,j)) * rainrate &
                -(1.-intercept_eff) * cveg(i,j) * rainrate &
                + min(0., wl_tend_lim - wl_tend_sum)
            interception(i,j) = max(0., wl_tend_lim)

            ! Integrate
            wl(i,j) = wlm(i,j) + rk3coef * wl_tend_lim

        end do
    end do

end subroutine calc_liquid_reservoir


!
! Calculate root fraction weighted mean soil water content
!
subroutine calc_theta_mean(tile)
    use modglobal,   only : i1, j1
    use modsurfdata, only : phiw
    implicit none

    type(T_lsm_tile), intent(inout) :: tile
    integer :: i, j, k, si
    real :: theta_lim

    tile%phiw_mean(:,:) = 0.

    do k=1, kmax_soil
        do j=2,j1
            do i=2,i1
                si = soil_index(i,j,k)
                theta_lim = max(phiw(i,j,k), theta_wp(si))
                tile%phiw_mean(i,j) = tile%phiw_mean(i,j) + tile%root_frac(i,j,k) * &
                    (theta_lim - theta_wp(si)) / (theta_fc(si) - theta_wp(si))
            end do
        end do
    end do

end subroutine calc_theta_mean

!
! Calculate canopy and soil resistances using Jarvis-Stewart method.
!
subroutine calc_canopy_resistance_js
    use modglobal,   only : i1, j1
    use modfields,   only : thl0, qt0, exnf, presf
    use modsurface,  only : ps
    use modraddata,  only : swd
    use modsurfdata, only : phiw
    implicit none

    integer :: i, j, k, si
    real :: swd_pos, T, esat, e, theta_min, theta_rel

    ! Constants f1 calculation:
    real, parameter :: a_f1 = 0.81
    real, parameter :: b_f1 = 0.004
    real, parameter :: c_f1 = 0.05

    k = kmax_soil
    do j=2,j1
        do i=2,i1
            si = soil_index(i,j,k)

            ! f1: reduction vegetation resistance as f(sw_in):
            swd_pos = max(0._field_r, -swd(i,j,1))
            f1(i,j) = 1./min(1., (b_f1*swd_pos + c_f1) / (a_f1 * (b_f1*swd_pos + 1.)))

            ! f2: reduction vegetation resistance as f(theta):
            do ilu=1,nlu
              if (tile(ilu)%lveg) then
                tile(ilu)%f2(i,j) = 1./min(1., max(1.e-9, tile(ilu)%phiw_mean(i,j)))
              endif
            enddo

            ! f3: reduction vegetation resistance as f(VPD) (high veg only):
            T    = thl0(i,j,1) * exnf(1)
            esat = 0.611e3 * exp(17.2694 * (T - 273.16) / (T - 35.86))
            e    = qt0(i,j,1) * presf(1) / 0.622

            do ilu=1,nlu
              if (tile(ilu)%lveg) then
                tile(ilu)%f3(i,j) = 1./exp(-tile(ilu)%gD(i,j) * (esat-e))
              endif
            enddo

            ! f2b: reduction soil resistance as f(theta)
            theta_min = cveg(i,j) * theta_wp(si) + (1.-cveg(i,j)) * theta_res(si);
            theta_rel = (phiw(i,j,k) - theta_min) / (theta_fc(si) - theta_min);
            f2b(i,j)  = 1./min(1., max(1.e-9, theta_rel))

            ! Calculate canopy and soil resistance
            do ilu=1,nlu
              if (tile(ilu)%lveg) then
                tile(ilu)%rs(i,j) = tile(ilu)%rs_min(i,j) / tile(ilu)%lai(i,j) * tile(ilu)%f2(i,j) * tile(ilu)%f3(i,j) * f1(i,j)
              else if (trim(tile(ilu)%lushort) == 'bs' .or. trim(tile(ilu)%lushort) == 'brn') then !TODO; special function for bare soil
                tile(ilu)%rs(i,j) = tile(ilu)%rs_min(i,j) / f2b(i,j)
              else if (ilu == ilu_ws) then
                tile(ilu)%rs(i,j) = 0
              else
                tile(ilu)%rs(i,j) = tile(ilu)%rs_min(i,j)
              endif
            enddo

        end do
    end do

end subroutine calc_canopy_resistance_js

!
! Calculate canopy and soil resistances using A-Gs (plant physiology).
! In addition, this calculates/sets the surface CO2 fluxes...
!
subroutine calc_canopy_resistance_ags

end subroutine calc_canopy_resistance_ags

!
! Calculate Obukhov length, ustar, and aerodynamic resistance, for all tiles
!
subroutine calc_stability
  use modglobal, only : i1, j1, cu, cv, rv, rd
  use modfields, only : u0, v0, thl0, qt0
  implicit none

  real, parameter :: du_min = 0.1
  real :: du, dv
  integer :: i, j

  ! Calculate properties shared by all tiles:
  ! Absolute wind speed difference, and virtual potential temperature atmosphere
  do j=2,j1
      do i=2,i1
          du = 0.5*(u0(i,j,1) + u0(i+1,j,1)) + cu
          dv = 0.5*(v0(i,j,1) + v0(i,j+1,1)) + cv
          du_tot(i,j) = max(0.1, sqrt(du**2 + dv**2))

          thv_1(i,j) = thl0(i,j,1)  * (1.+(rv/rd-1.)*qt0(i,j,1))
      end do
  end do

  do ilu=1, nlu
    call calc_obuk_ustar_ra(tile(ilu))
  end do
end subroutine calc_stability

!
! Calculate Obukhov length and ustar, for single tile
!
subroutine calc_obuk_ustar_ra(tile)
    use modglobal, only : i1, j1, rd, rv, grav, zf
    use modfields, only : u0, v0
    implicit none

    type(T_lsm_tile), intent(inout) :: tile
    integer :: i, j
    real :: thvs

    do j=2,j1
        do i=2,i1
            if (tile%frac(i,j) > 0) then
                ! Buoyancy difference surface - atmosphere
                thvs = tile%thlskin(i,j) * (1.+(rv/rd-1.)*tile%qtskin(i,j))
                tile%db(i,j) = grav/thvs * (thv_1(i,j) - thvs)

                ! Iteratively find Obukhov length
                tile%obuk(i,j) = calc_obuk_dirichlet( &
                    tile%obuk(i,j), du_tot(i,j), tile%db(i,j), real(zf(1), 8), tile%z0m(i,j), tile%z0h(i,j))
            end if
        end do
    end do

    do j=2,j1
        do i=2,i1
            if (tile%frac(i,j) > 0) then
                ! Calculate friction velocity and aerodynamic resistance
                tile%ustar(i,j) = du_tot(i,j) * fm(real(zf(1), 8), tile%z0m(i,j), tile%obuk(i,j))
                tile%ra(i,j)    = 1./(tile%ustar(i,j) * fh(real(zf(1), 8), tile%z0h(i,j), tile%obuk(i,j)))
            end if
        end do
    end do

end subroutine calc_obuk_ustar_ra

!
! Calculate surface temperature for each tile, and calculate
! surface fluxes (H, LE, G0, wthl, wqt) and values (thlskin, qtskin)
!
subroutine calc_tile_bcs(tile)
    use modglobal,   only : i1, j1, cp, rlv, boltz
    use modfields,   only : exnh, exnf, presh, thl0, qt0, rhof
    use modraddata,  only : swd, swu, lwd, lwu
    use modsurfdata, only : ps, tsoil
    implicit none

    type(T_lsm_tile), intent(inout) :: tile
    integer :: i, j
    real :: Ts, thvs, esats, qsats, desatdTs, dqsatdTs, &
        rs_lim, fH, fLE, fG, num, denom, Ta, qsat_new, &
        rhocp_i, rholv_i, Qnet

    rhocp_i = 1. / (rhof(1) * cp)
    rholv_i = 1. / (rhof(1) * rlv)

    do j=2, j1
        do i=2, i1
            if (tile%frac(i,j) > 0) then

                ! Disable canopy resistance in case of dew fall
                Ts    = tile%thlskin(i,j) * exnh(1)
                esats = 0.611e3 * exp(17.2694 * (Ts - 273.16) / (Ts - 35.86))
                qsats = 0.622 * esats / ps

                if (qsats < qt0(i,j,1)) then
                    rs_lim = 0.
                else
                    rs_lim = tile%rs(i,j)
                end if

                ! Calculate recuring terms
                ! NOTE: this should use the surface density, not rhof(1).
                !       Not sure if rhoh is available somewhere...
                fH  = rhof(1) * cp  / tile%ra(i,j)
                fLE = rhof(1) * rlv / (tile%ra(i,j) + rs_lim)

                if (tile%db(i,j) > 0) then
                    fG = tile%lambda_stable(i,j)
                else
                    fG = tile%lambda_unstable(i,j)
                end if

                ! Net radiation; negative sign = net input of radiation at surface
                Qnet = swd(i,j,1) + swu(i,j,1) + lwd(i,j,1) + lwu(i,j,1)

                ! Solve new skin temperature from SEB
                desatdTs = esats * (17.2694 / (Ts - 35.86) - 17.2694 * (Ts - 273.16) / (Ts - 35.86)**2.)
                dqsatdTs = 0.622 * desatdTs / ps
                Ta = thl0(i,j,1) * exnf(1)

                num = -(Qnet - lwu(i,j,1) &
                      - fH * Ta + (qsats - dqsatdTs * Ts - qt0(i,j,1)) * fLE &
                      - fG * tsoil(i,j,kmax_soil) - 3.*boltz * Ts**4)
                denom = (fH + fLE * dqsatdTs + fG + 4.*boltz * Ts**3)
                tile%tskin(i,j) = num / denom

                ! Update qsat with linearised relation, to make sure that the SEB closes
                qsat_new = qsats + dqsatdTs * (tile%tskin(i,j) - Ts)

                ! Calculate energetic surface fluxes
                tile%H (i,j) = fH  * (tile%tskin(i,j) - Ta)
                tile%LE(i,j) = fLE * (qsat_new - qt0(i,j,1))
                tile%G (i,j) = fG  * (tsoil(i,j,kmax_soil) - tile%tskin(i,j))

                ! Calculate kinematic surface fluxes
                tile%wthl(i,j) = tile%H (i,j) * rhocp_i
                tile%wqt (i,j) = tile%LE(i,j) * rholv_i

                ! Calculate surface values
                tile%thlskin(i,j) = thl0(i,j,1) + tile%wthl(i,j) * tile%ra(i,j)
                tile%qtskin (i,j) = qt0 (i,j,1) + tile%wqt (i,j) * tile%ra(i,j)
            end if
        end do
    end do

end subroutine calc_tile_bcs

!
! Calculate BCs and fluxes for the water tile
!
!subroutine calc_water_bcs
subroutine calc_water_bcs(tile)
    use modglobal,   only : i1, j1, cp, rlv
    use modfields,   only : exnh, thl0, qt0, rhof
    use modsurfdata, only : ps
    implicit none

    type(T_lsm_tile), intent(inout) :: tile
    integer :: i, j
    real :: esats

    do j=2, j1
      do i=2, i1
        if (tile%frac(i,j) > 0) then
            !! Calculate BCs. `tile_aq%tskin` is fixed in time (for now...)
            !tile_aq%thlskin(i,j) = tile_aq%tskin(i,j) / exnh(1)
            !esats = 0.611e3 * exp(17.2694 * (tile_aq%tskin(i,j) - 273.16) / (tile_aq%tskin(i,j) - 35.86))
            !tile_aq%qtskin(i,j) = 0.622 * esats / ps
            !
            !! Calculate kinematic fluxes
            !tile_aq%wthl(i,j) = 1./tile_aq%ra(i,j) * (tile_aq%thlskin(i,j) - thl0(i,j,1))
            !tile_aq%wqt (i,j) = 1./tile_aq%ra(i,j) * (tile_aq%qtskin (i,j) - qt0 (i,j,1))
            !
            !! Calculate energetic fluxes
            !tile_aq%H (i,j) = tile_aq%wthl(i,j) * rhof(1) * cp
            !tile_aq%LE(i,j) = tile_aq%wqt (i,j) * rhof(1) * rlv
            !tile_aq%G (i,j) = 0.
            ! Calculate BCs. `tile_aq%tskin` is fixed in time (for now...)

            tile%thlskin(i,j) = tile%tskin(i,j) / exnh(1)
            esats = 0.611e3 * exp(17.2694 * (tile%tskin(i,j) - 273.16) / (tile%tskin(i,j) - 35.86))
            tile%qtskin(i,j) = 0.622 * esats / ps

            ! Calculate kinematic fluxes
            tile%wthl(i,j) = 1./tile%ra(i,j) * (tile%thlskin(i,j) - thl0(i,j,1))
            tile%wqt (i,j) = 1./tile%ra(i,j) * (tile%qtskin (i,j) - qt0 (i,j,1))

            ! Calculate energetic fluxes
            tile%H (i,j) = tile%wthl(i,j) * rhof(1) * cp
            tile%LE(i,j) = tile%wqt (i,j) * rhof(1) * rlv
            tile%G (i,j) = 0.
        end if
      end do
    end do

end subroutine calc_water_bcs

!
! Set grid point mean boundary conditions, as used by
! the diffusion scheme, thermodynamics, ...
!
subroutine calc_bulk_bcs
    use modglobal,   only : i1, j1, i2, j2, cp, rlv, fkar, zf, cu, cv, grav, rv, rd, lopenbc,lboundary,lperiodic, eps1
    use modfields,   only : rhof, thl0, qt0, u0, v0, thvh
    use modsurface,  only : phim, phih
    use modmpi,      only : excjs
    use modopenboundary, only : openboundary_excjs
    use modsurfdata, only : &
        H, LE, G0, tskin, qskin, thlflux, qtflux, dthldz, dqtdz, &
        dudz, dvdz, ustar, obl, cliq, ra, rsveg, rssoil
    implicit none

    integer :: i, j
    real :: rhocp_i, rholv_i, ucu, vcv, bflux
    real, pointer :: ustar_3D(:,:,:)

    rhocp_i = 1. / (rhof(1) * cp)
    rholv_i = 1. / (rhof(1) * rlv)

    ! Calculate surface temperature for each tile, and calculate
    ! surface fluxes (H, LE, G0, wthl, wqt) and values (thlskin, qtskin)
    do ilu=1,nlu
      if (tile(ilu)%laqu) then
        call calc_water_bcs(tile(ilu))
      else
        call calc_tile_bcs(tile(ilu))
      endif
    enddo

    do j=2,j1
        do i=2,i1
            H(i,j) = 0
            LE(i,j) = 0
            G0(i,j) = 0
            ustar(i,j) = 0
            tskin(i,j) = 0
            qskin(i,j) = 0
            rsveg(i,j) = 0
            rssoil(i,j) = 0
            do ilu=1,nlu
              H(i,j)      = H(i,j)     + tile(ilu)%frac(i,j) * tile(ilu)%H(i,j)
              LE(i,j)     = LE(i,j)    + tile(ilu)%frac(i,j) * tile(ilu)%LE(i,j)
              G0(i,j)     = G0(i,j)    + tile(ilu)%frac(i,j) * tile(ilu)%G(i,j)
              ustar(i,j)  = ustar(i,j) + tile(ilu)%frac(i,j) * tile(ilu)%ustar(i,j)
              tskin(i,j)  = tskin(i,j) + tile(ilu)%frac(i,j) * tile(ilu)%thlskin(i,j)
              qskin(i,j)  = qskin(i,j) + tile(ilu)%frac(i,j) * tile(ilu)%qtskin(i,j)
            enddo

            ! Kinematic surface fluxes
            thlflux(i,j) =  H(i,j)  * rhocp_i
            qtflux (i,j) =  LE(i,j) * rholv_i

            ! Calculate mean Obukhov length from mean fluxes
            bflux = grav/thvh(1) * (thlflux(i,j) * (1.-(1.-rv/rd)*qskin(i,j)) - &
                (1.-rv/rd)*tskin(i,j)*qtflux(i,j))
            obl(i,j) = -ustar(i,j)**3 / (fkar * bflux)

            ! MO gradients
            dthldz(i,j) = -thlflux(i,j) / (fkar * zf(1) * ustar(i,j)) * phih(zf(1)/obl(i,j))
            dqtdz (i,j) = -qtflux (i,j) / (fkar * zf(1) * ustar(i,j)) * phih(zf(1)/obl(i,j))

            ! NOTE: dudz, dvdz are at the grid center (full level), not the velocity locations.
            ucu = 0.5*(u0(i,j,1) + u0(i+1,j,1))+cu
            vcv = 0.5*(v0(i,j,1) + v0(i,j+1,1))+cv

            dudz(i,j) = ustar(i,j) / (fkar * zf(1)) * phim(zf(1)/obl(i,j)) * (ucu/du_tot(i,j))
            dvdz(i,j) = ustar(i,j) / (fkar * zf(1)) * phim(zf(1)/obl(i,j)) * (vcv/du_tot(i,j))

            ! Cyclic BCs where needed.
            ustar_3D(1:i2,1:j2,1:1) => ustar

            if(lopenbc) then ! Only use periodicity for non-domain boundaries when openboundaries are used
              call openboundary_excjs(ustar_3D, 2,i1,2,j1,1,1,1,1, &
                   & (.not.lboundary(1:4)).or.lperiodic(1:4))
           else
              !call excjs(ustar,2,i1,2,j1,1,1,1,1)
              call excjs(ustar_3D,2,i1,2,j1,1,1,1,1)
           endif

            ! Just for diagnostics (modlsmcrosssection)
            do ilu=1,nlu
              !if (trim(tile(ilu)%lushort) == 'ws') then
              if (ilu == ilu_ws) then
                cliq(i,j) = tile(ilu)%frac(i,j) / land_frac(i,j)
              endif
            enddo

            ! Calculate ra consistent with mean flux and temperature difference:
            ra(i,j) = (tskin(i,j) - thl0(i,j,1)) / thlflux(i,j)

            if (cveg(i,j) > 0) then
                do ilu=1,nlu
                  if (tile(ilu)%lveg) then
                    rsveg(i,j) = rsveg(i,j) + tile(ilu)%frac(i,j) * tile(ilu)%rs(i,j)
                    ! if ( tile(ilu)%frac(i,j) <= 0.0 ) then
                    !   cycle
                    ! endif
                    ! rsveg(i,j) = 1/(rsveg(i,j)**(-1) + (tile(ilu)%rs(i,j) * tile(ilu)%frac(i,j))**(-1))
                  endif
                enddo
                rsveg(i,j) = rsveg(i,j) / cveg(i,j)
            else
                rsveg(i,j) = 0.
            end if
            do ilu=1,nlu
            ! TODO: flexible solution for LU types with soil resistance
              if ( .not. (tile(ilu)%lveg .or. tile(ilu)%laqu .or. (ilu == ilu_ws) ) ) then
                  rssoil(i,j) = rssoil(i,j) + tile(ilu)%rs(i,j) * tile(ilu)%frac(i,j)
                  ! if ( tile(ilu)%frac(i,j) <= 0.0 ) then
                  !   cycle
                  ! endif
                  ! rssoil(i,j) = 1/(rssoil(i,j)**(-1) + (tile(ilu)%rs(i,j) * tile(ilu)%frac(i,j))**(-1))
              end if
            end do
        end do
    end do

end subroutine calc_bulk_bcs

!
! Interpolation soil from full to half levels,
! using various interpolation methods
!
subroutine interpolate_soil(fieldh, field, iinterp)
    use modglobal, only : i1, j1
    implicit none
    real, intent(inout) :: fieldh(:,:,:)
    real, intent(in)    :: field(:,:,:)
    integer, intent(in) :: iinterp
    integer i, j, k

    if (iinterp == iinterp_amean) then
        do k=2,kmax_soil
            do j=2,j1
                do i=2,i1
                    fieldh(i,j,k) = 0.5*(field(i,j,k-1) + field(i,j,k))
                end do
            end do
        end do
    else if (iinterp == iinterp_gmean) then
        do k=2,kmax_soil
            do j=2,j1
                do i=2,i1
                    fieldh(i,j,k) = sqrt(field(i,j,k-1) * field(i,j,k))
                end do
            end do
        end do
    else if (iinterp == iinterp_hmean) then
        do k=2,kmax_soil
            do j=2,j1
                do i=2,i1
                    fieldh(i,j,k) = ((dz_soil(k-1)+dz_soil(k))*field(i,j,k-1)*field(i,j,k)) / &
                            (dz_soil(k-1)*field(i,j,k) + dz_soil(k)*field(i,j,k-1))
                end do
            end do
        end do
    else if (iinterp == iinterp_max) then
        do k=2,kmax_soil
            do j=2,j1
                do i=2,i1
                    fieldh(i,j,k) = max(field(i,j,k-1), field(i,j,k))
                end do
            end do
        end do
    else
        print*,'ERROR: invalid soil interpolation type iinterp=',iinterp
        stop
    end if

end subroutine interpolate_soil

!
! Calculate temperature diffusivity soil at full and half levels
!
subroutine calc_thermal_properties
    use modglobal, only : i1, j1, gamma_t_matrix, gamma_t_water
    use modsurfdata, only : lambda, lambdah, phiw
    implicit none

    integer :: i, j, k, si
    real :: lambda_t_sat, kersten, gamma_t

    ! Calculate diffusivity heat
    do k=1,kmax_soil
        do j=2,j1
            do i=2,i1
                si = soil_index(i,j,k)

                ! Heat conductivity at saturation (from IFS code..)
                lambda_T_sat =   gamma_T_matrix**(1.-theta_sat(si)) &
                               * gamma_T_water**phiw(i,j,k) &
                               * 2.2**(theta_sat(si) - phiw(i,j,k))

                ! Kersten number for fine soils [IFS eq 8.64] (-)
                kersten = log10(max(0.1, phiw(i,j,k) / theta_sat(si))) + 1.

                ! Heat conductivity soil [IFS eq 8.62] (W m-1 K-1)
                gamma_t = kersten * (lambda_T_sat - gamma_t_dry(si)) + gamma_t_dry(si)

                ! Heat diffusivity (m2 s-1)
                lambda(i,j,k) = gamma_t / rho_C(si)
            end do
        end do
    end do

    ! Interpolate to half levels
    call interpolate_soil(lambdah, lambda, iinterp_t)

end subroutine calc_thermal_properties

!
! Calculate soil moisture diffusivity and conductivity,
! using van Genuchten parameterisation.
!
subroutine calc_hydraulic_properties
    use modglobal, only : i1, j1
    use modsurfdata, only : phiw, lambdas, lambdash, gammas, gammash
    implicit none

    integer :: i, j, k, si
    real :: theta_lim, theta_norm

    ! Calculate diffusivity and conductivity soil moisture
    do k=1,kmax_soil
        do j=2,j1
            do i=2,i1
                si = soil_index(i,j,k)

                ! Limit soil moisture just above the residual soil moisture content
                ! and just below saturation
                theta_lim = min(max(phiw(i,j,k), 1.001*theta_res(si)), theta_sat(si))

                ! Dimensionless soil water content
                theta_norm = (theta_lim - theta_res(si)) / (theta_sat(si) - theta_res(si))

                ! Calculate & limit the diffusivity
                lambdas(i,j,k) = calc_diffusivity_vg( &
                        theta_norm, vg_a(si), vg_l(si), vg_m(si), &
                        gamma_theta_sat(si), theta_sat(si), theta_res(si))
                lambdas(i,j,k) = max(min(lambdas(i,j,k), lambda_theta_max(si)), lambda_theta_min(si))

                ! Calculate & limit the conductivity
                gammas(i,j,k) = calc_conductivity_vg( &
                        theta_norm, vg_l(si), vg_m(si), gamma_theta_sat(si))
                gammas(i,j,k) = max(min(gammas(i,j,k), gamma_theta_max(si)), gamma_theta_min(si))
            end do
        end do
    end do

    ! Interpolate to half levels
    call interpolate_soil(lambdash, lambdas, iinterp_theta)
    call interpolate_soil(gammash,  gammas,  iinterp_theta)

    ! Optionally, set free drainage bottom BC
    if (lfreedrainage) then
        gammash(:,:,1) = gammash(:,:,2)
    else
        gammash(:,:,1) = 0.
    end if

end subroutine calc_hydraulic_properties

!
! Calculate soil moisture tendency from root water extraction
!
subroutine calc_root_water_extraction
    use modglobal, only : i1, j1, rhow, rlv
    use modsurfdata, only : phiw
    implicit none

    integer :: i, j, k
    real :: phiw_rf, phi_frac, LE
    real, parameter :: fac = 1./(rhow * rlv)

    phiw_source = 0
    do j=2, j1
      do i=2, i1
        do ilu=1,nlu
          if (.not. tile(ilu)%lveg) then
            cycle
          else
            LE = tile(ilu)%frac(i,j) * tile(ilu)%LE(i,j)
            phiw_rf = 0.

            do k=1, kmax_soil
                phiw_rf = phiw_rf + tile(ilu)%root_frac(i,j,k) * phiw(i,j,k)
            end do

            do k=1, kmax_soil
                phi_frac = tile(ilu)%root_frac(i,j,k) * phiw(i,j,k) / phiw_rf

                phiw_source(i,j,k) = phiw_source(i,j,k) &
                                     - max(0., LE) * fac * dzi_soil(k) * phi_frac

            end do

          end if
        end do
      end do
    end do

end subroutine calc_root_water_extraction

!
! Solve diffusion equation (explicit) for soil temperature.
! Top flux = G0/rho_C, bottom flux = zero.
!
subroutine integrate_t_soil
    use modglobal, only : rk3step, rdt, i1, j1
    use modsurfdata, only : tsoil, tsoilm, lambdah, G0

    implicit none
    integer :: i, j, k, si
    real :: tend, rk3coef, flux_top

    rk3coef = rdt / (4. - dble(rk3step))
    if(rk3step == 1) tsoilm(:,:,:) = tsoil(:,:,:)

    ! Top soil layer
    k = kmax_soil
    do j=2,j1
        do i=2,i1
            si = soil_index(i,j,k)
            flux_top = G0(i,j) / rho_C(si)
            tend = (-flux_top - (lambdah(i,j,k) * (tsoil(i,j,k) - tsoil(i,j,k-1)) * dzhi_soil(k)))*dzi_soil(k)

            tsoil(i,j,k) = tsoilm(i,j,k) + rk3coef * tend
        end do
    end do

    ! Bottom soil layer
    k = 1
    do j=2,j1
        do i=2,i1
            tend = ((lambdah(i,j,k+1) * (tsoil(i,j,k+1) - tsoil(i,j,k)) * dzhi_soil(k+1)))*dzi_soil(k)

            tsoil(i,j,k) = tsoilm(i,j,k) + rk3coef * tend
        end do
    end do

    ! Interior
    do k=2,kmax_soil-1
        do j=2,j1
            do i=2,i1
                tend = ((lambdah(i,j,k+1) * (tsoil(i,j,k+1) - tsoil(i,j,k  )) * dzhi_soil(k+1)) &
                      - (lambdah(i,j,k  ) * (tsoil(i,j,k  ) - tsoil(i,j,k-1)) * dzhi_soil(k  ))) * dzi_soil(k)

                tsoil(i,j,k) = tsoilm(i,j,k) + rk3coef * tend
            end do
        end do
    end do

end subroutine integrate_t_soil

!
! Solve diffusion equation (explicit) for soil moisture,
! including a source term for root water extraction.
!
subroutine integrate_theta_soil
    use modglobal, only : rk3step, rdt, i1, j1, rhow, rlv
    use modsurfdata, only : phiw, phiwm, lambdash, gammash
    use modmpi, only : myidx, myidy
    implicit none
    integer :: i, j, k, si
    real :: tend, rk3coef, flux_top, fac

    rk3coef = rdt / (4. - dble(rk3step))
    if(rk3step == 1) phiwm(:,:,:) = phiw(:,:,:)

    fac = 1./(rhow * rlv)

     ! Top soil layer
    k = kmax_soil
    do j=2,j1
        do i=2,i1
          do ilu=1,nlu
            !if (trim(tile(ilu)%lushort) == 'bs') then
            if (trim(tile(ilu)%lushort) == 'bs' .or. trim(tile(ilu)%lushort) == 'brn') then !TODO; special function for bare soil
              flux_top = tile(ilu)%frac(i,j) * tile(ilu)%LE(i,j) * fac + throughfall(i,j)
              tend = (-flux_top - (lambdash(i,j,k) * (phiw(i,j,k) - phiw(i,j,k-1)) * dzhi_soil(k)))*dzi_soil(k) &
                    - gammash(i,j,k) * dzi_soil(k) + phiw_source(i,j,k)
              phiw(i,j,k) = phiwm(i,j,k) + rk3coef * tend
            end if
          end do
        end do
    end do

    ! Bottom soil layer
    k = 1
    do j=2,j1
        do i=2,i1
            tend = ((lambdash(i,j,k+1) * (phiw(i,j,k+1) - phiw(i,j,k)) * dzhi_soil(k+1)))*dzi_soil(k) &
                   + (gammash(i,j,k+1) - gammash(i,j,k)) * dzi_soil(k) + phiw_source(i,j,k)

            phiw(i,j,k) = phiwm(i,j,k) + rk3coef * tend
        end do
    end do

    ! Interior
    do k=2,kmax_soil-1
        do j=2,j1
            do i=2,i1
                tend = ((lambdash(i,j,k+1) * (phiw(i,j,k+1) - phiw(i,j,k  )) * dzhi_soil(k+1)) &
                      - (lambdash(i,j,k  ) * (phiw(i,j,k  ) - phiw(i,j,k-1)) * dzhi_soil(k  ))) * dzi_soil(k) &
                      + (gammash(i,j,k+1) - gammash(i,j,k)) * dzi_soil(k) + phiw_source(i,j,k)

                phiw(i,j,k) = phiwm(i,j,k) + rk3coef * tend
            end do
        end do
    end do

    ! Range check of phiw
    if (maxval(phiw) > 1 .or. minval(phiw) < 0) then
       write(*,*) 'phiw out or range 0...1'
       write(*,*) 'myid{x,y} =', myidx, myidy
       write(*,*) 'max', maxval(phiw), 'at', maxloc(phiw)
       write(*,*) 'min', minval(phiw), 'at', minloc(phiw)
       !stop
    end if

end subroutine integrate_theta_soil

!
! Initialise the land-surface model
!
subroutine initlsm
    use modglobal,   only : ifnamopt, fname_options, checknamelisterror, lwarmstart
    use modmpi,      only : myid, comm3d, mpierr, D_MPI_BCAST
    use modsurfdata, only : isurf
    use modemisdata, only : l_emission

    implicit none

    integer :: ierr
    logical :: lheterogeneous

    ! Namelist definition
    namelist /NAMLSM/ &
        lheterogeneous, lfreedrainage, lags, dz_soil, iinterp_t, iinterp_theta, co2_index, nlu
    llsm = (isurf == 11)

    if (llsm) then

        allocate(dz_soil(kmax_soil))
        ! Read namelist
        if (myid == 0) then
            open(ifnamopt, file=fname_options, status='old', iostat=ierr)
            read(ifnamopt, NAMLSM, iostat=ierr)
            call checknamelisterror(ierr, ifnamopt, 'NAMLSM')
            write(6, NAMLSM)
            close(ifnamopt)
        end if

        ! Broadcast namelist values to all MPI tasks
        call D_MPI_BCAST(lheterogeneous,  1, 0, comm3d, mpierr)
        call D_MPI_BCAST(lfreedrainage,   1, 0, comm3d, mpierr)
        call D_MPI_BCAST(lags,            1, 0, comm3d, mpierr)
        call D_MPI_BCAST(nlu,             1, 0, comm3d, mpierr)
        call D_MPI_BCAST(iinterp_t,       1, 0, comm3d, mpierr)
        call D_MPI_BCAST(iinterp_theta,   1, 0, comm3d, mpierr)
        call D_MPI_BCAST(co2_index,       1, 0, comm3d, mpierr)
        call D_MPI_BCAST(dz_soil, kmax_soil, 0, comm3d, mpierr)

        ! check if nlu_file==nlu-1 ('wet skin' is not in file)
        if (.not. lheterogeneous) then
          write(6,"(A100, i3)") "Homogeneous land use; Note that 1 additional LU type (ws) is added on runtime. &
                                    Include it in nlu  ", nlu
        end if

        allocate(tile(nlu), stat=ierr)
        if (ierr/=0) stop

        ! Checks on input
        if (lags .and. .not. l_emission .and. (co2_index == -1)) then
            print*,'With A-Gs enabled and without the emission module, the `co2_index`'
            print*,'in the scalar array should be specified in the `NAMLSM` group.'
            stop
        endif

        ! Create/calculate soil grid properties
        call create_soil_grid

        ! Allocate required fields in modsurfacedata, and arrays / tiles from this module:
        call allocate_fields

        if (lheterogeneous) then
            ! Initialise heterogeneous LSM from external input
            !call init_heterogeneous
            call init_heterogeneous_nc
        else
            ! Initialise homogeneous LSM from namelist input
            call init_homogeneous
        end if

        ! Read the soil parameter table
        call read_soil_table

        ! Calculate derived soil properties
        call calc_soil_properties

        ! Calculate root fractions soil (low and high veg)
        call calc_root_fractions

    end if

end subroutine initlsm

!
! Cleanup (deallocate) the land-surface model
!
subroutine exitlsm
    use modsurfdata, only : &
        tsoil, tsoilm, phiw, phiwm, &
        lambda, lambdah, lambdas, lambdash, gammas, gammash, &
        H, LE, G0, Qnet, wl, wlm, &
        tendskin, cliq, rsveg, rssoil
    implicit none

    if (.not. llsm) return

    ! Allocated from `read_soil_table`:
    deallocate( theta_res, theta_wp, theta_fc, theta_sat, gamma_theta_sat, vg_a, vg_l, vg_n )

    ! Allocated from `allocate_fields`:
    deallocate( soil_index, tsoil, tsoilm, phiw, phiwm, phiw_source )
    deallocate( wl, wlm, wl_max )
    deallocate( throughfall, interception )
    deallocate( f1, f2b )
    deallocate( du_tot, thv_1, land_frac, cveg )
    deallocate( lambda, lambdah, lambdas, lambdash, gammas, gammash )
    deallocate( Qnet, H, LE, G0 )
    deallocate( tendskin, cliq, rsveg, rssoil )

    ! Allocated from `create_soil_grid`:
    deallocate( z_soil, dz_soil, dzi_soil, zh_soil, dzh_soil, dzhi_soil )

    if (lags) deallocate(an_co2, resp_co2)

    ! Tiles, allocated from `allocate_tile`:
    do ilu=1,nlu
      call deallocate_tile(tile(ilu))
    enddo

end subroutine exitlsm

!
! Calculate soil grid properties
!
subroutine create_soil_grid
    implicit none
    integer :: k

    allocate(z_soil   (kmax_soil  ))
    allocate(dzi_soil (kmax_soil  ))
    allocate(zh_soil  (kmax_soil+1))
    allocate(dzh_soil (kmax_soil+1))
    allocate(dzhi_soil(kmax_soil+1))

    !
    ! Full level is in the middle of two half levels
    !
    z_size_soil = -sum(dz_soil)

    ! Half level height
    zh_soil(1) = z_size_soil
    zh_soil(kmax_soil+1) = 0.
    do k=2, kmax_soil
        zh_soil(k) = zh_soil(k-1) + dz_soil(k-1)
    end do

    ! Full level height
    do k=1, kmax_soil
        z_soil(k) = 0.5*(zh_soil(k) + zh_soil(k+1))
    end do

    do k=2, kmax_soil
        dzh_soil(k) = z_soil(k) - z_soil(k-1)
    end do

    dzh_soil(1) = 2*(z_soil(1) - zh_soil(1))
    dzh_soil(kmax_soil+1) = 2*(-z_soil(kmax_soil))

    ! Inverse grid spacings
    dzi_soil(:) = 1./dz_soil(:)
    dzhi_soil(:) = 1./dzh_soil(:)

    !
    ! Half level is in the middle of two full levels
    !
    !! Half level heights
    !zh_soil(1) = z_size_soil
    !zh_soil(kmax_soil+1) = 0.

    !do k=2, kmax_soil
    !    zh_soil(k) = 0.5*(z_soil(k) + z_soil(k-1))
    !end do

    !! Grid spacing full and half levels
    !do k=1, kmax_soil
    !    dz_soil(k) = zh_soil(k+1) - zh_soil(k)
    !end do

    !do k=2, kmax_soil
    !    dzh_soil(k) = z_soil(k) - z_soil(k-1)
    !end do

    !dzh_soil(1) = 2*(z_soil(1) - zh_soil(1))
    !dzh_soil(kmax_soil+1) = 2*(-z_soil(kmax_soil))

    !! Inverse grid spacings
    !dzi_soil(:) = 1./dz_soil(:)
    !dzhi_soil(:) = 1./dzh_soil(:)

end subroutine create_soil_grid

!
! Allocate all LSM fields
!
subroutine allocate_fields
    use modglobal, only : i2, j2
    use modsurfdata, only : &
        tsoil, tsoilm, phiw, phiwm, &
        lambda, lambdah, lambdas, lambdash, gammas, gammash, &
        H, LE, G0, Qnet, wl, wlm, &
        tendskin, cliq, rsveg, rssoil
    implicit none

    ! Allocate soil variables
    allocate(soil_index(i2, j2, kmax_soil))

    allocate(tsoil (i2, j2, kmax_soil))
    allocate(tsoilm(i2, j2, kmax_soil))
    allocate(phiw  (i2, j2, kmax_soil))
    allocate(phiwm (i2, j2, kmax_soil))

    allocate(wl (i2, j2))
    allocate(wlm(i2, j2))

    allocate(wl_max(i2, j2))

    allocate(phiw_source(i2, j2, kmax_soil))

    allocate(throughfall(i2, j2))
    allocate(interception(i2, j2))

    allocate(f1(i2, j2))
    allocate(f2b(i2, j2))

    allocate(du_tot(i2, j2))
    allocate(thv_1(i2, j2))
    allocate(land_frac(i2, j2))
    allocate(cveg(i2, j2))

    ! NOTE: names differ from what is described in modsurfdata!
    ! Diffusivity temperature:
    allocate(lambda (i2, j2, kmax_soil  ))
    allocate(lambdah(i2, j2, kmax_soil+1))

    ! Diffusivity theta:
    allocate(lambdas (i2, j2, kmax_soil))
    allocate(lambdash(i2, j2, kmax_soil+1))

    ! Conductivity theta:
    allocate(gammas (i2, j2, kmax_soil))
    allocate(gammash(i2, j2, kmax_soil+1))

    ! Tile averaged surface fluxes
    allocate(Qnet (i2, j2))
    allocate(H    (i2, j2))
    allocate(LE   (i2, j2))
    allocate(G0   (i2, j2))

    ! TMP, to prevent segfaults in modlsmcrosssection.f90
    allocate(tendskin(i2, j2))
    allocate(cliq(i2, j2))
    allocate(rsveg(i2, j2))
    allocate(rssoil(i2, j2))

    ! A-Gs
    if (lags) then
        allocate(an_co2(i2, j2))
        allocate(resp_co2(i2, j2))
    endif

    ! Allocate the tiled variables
    do ilu=1,nlu
      call allocate_tile(tile(ilu))
    enddo

    ilu_ws = nlu

end subroutine allocate_fields

!
! Allocate all fields of a LSM tile
!
subroutine allocate_tile(tile)
    use modglobal, only : i2, j2
    implicit none
    type(T_lsm_tile), intent(inout) :: tile

    ! Static properties:
    allocate(tile % z0m(i2, j2))
    allocate(tile % z0h(i2, j2))

    ! Base and dynamic tile fraction:
    allocate(tile % base_frac(i2, j2))
    allocate(tile % frac(i2, j2))

    ! Monin-obukhov / surface layer:
    allocate(tile % obuk(i2, j2))
    allocate(tile % ustar(i2, j2))
    allocate(tile % ra(i2, j2))

    ! Conductivity skin layer:
    allocate(tile % lambda_stable(i2, j2))
    allocate(tile % lambda_unstable(i2, j2))

    ! Surface fluxes:
    allocate(tile % H(i2, j2))
    allocate(tile % LE(i2, j2))
    allocate(tile % G(i2, j2))
    allocate(tile % wthl(i2, j2))
    allocate(tile % wqt(i2, j2))

    ! Surface temperature and humidity:
    allocate(tile % tskin(i2, j2))
    allocate(tile % thlskin(i2, j2))
    allocate(tile % qtskin(i2, j2))

    ! Buoyancy difference surface-atmosphere
    allocate(tile % db(i2, j2))

    ! Vegetation properties:
    allocate(tile % lai(i2, j2))
    allocate(tile % rs_min(i2, j2))
    allocate(tile % rs(i2, j2))

    ! Root fraction in soil:
    allocate(tile % a_r(i2, j2))
    allocate(tile % b_r(i2, j2))
    allocate(tile % root_frac(i2, j2, kmax_soil))

    ! Root fraction weigthed mean soil water content
    allocate(tile % phiw_mean(i2, j2))

    ! Reduction functions canopy resistance
    allocate(tile % f2(i2, j2))
    allocate(tile % f3(i2, j2))
    allocate(tile % gD(i2, j2))

    ! In-canopy resistance parameters
!    allocate(tile % R_inc_b(i2, j2))
!    allocate(tile % R_inc_h(i2, j2))
    ! SAI = SAI_a * LAI + SAI_b
    allocate(tile % SAI_a(i2, j2))
    allocate(tile % SAI_b(i2, j2))
    ! Minimum correction factor for stomatal resistance
!    allocate(tile % fmin(i2, j2))
    ! Alpha value for light correction of stomatal resistance
!    allocate(tile % alpha(i2, j2))
    ! Min, optimum and max emperatures for temperature correction of stomatal resistance
!    allocate(tile % Tmin(i2, j2))
!    allocate(tile % Topt(i2, j2))
!    allocate(tile % Tmax(i2, j2))
    ! Maximum leaf stomatal conductance for ozone
!    allocate(tile % gs_max(i2, j2))
    ! Minimum and maximum vapour pressure deficit parameters
!    allocate(tile % vpd_min(i2, j2))
!    allocate(tile % vpd_max(i2, j2))
    ! Gamma parameter for calculating stomata compensation point
!    allocate(tile % gamma_stom(i2, j2))
    ! Gamma correction factor for calculating soil compensation point
!    allocate(tile % gamma_soil_c_fac(i2, j2))
    ! Gamma parameter for calculating soil compensation point
!    allocate(tile % gamma_soil_default(i2, j2))

end subroutine allocate_tile

!
! Deallocate all fields of a LSM tile
!
subroutine deallocate_tile(tile)
    implicit none
    type(T_lsm_tile), intent(inout) :: tile
    deallocate( tile%z0m, tile%z0h, tile%base_frac, tile%frac )
    deallocate( tile%obuk, tile%ustar, tile%ra )
    deallocate( tile%lambda_stable, tile%lambda_unstable )
    deallocate( tile%H, tile%LE, tile%G, tile%wthl, tile%wqt)
    deallocate( tile%tskin, tile%thlskin, tile%qtskin)
    deallocate( tile%db, tile%lai, tile%rs_min, tile%rs )
    deallocate( tile%a_r, tile%b_r, tile%root_frac, tile%phiw_mean )
    deallocate( tile%f2, tile%f3, tile%gD )
!   deallocate( tile%R_inc_b, tile%R_inc_h)
    deallocate( tile%SAI_a, tile%SAI_b)
!   deallocate( tile%fmin, tile%alpha)
!   deallocate( tile%Tmin, tile%Topt, tile%Tmax)
!   deallocate( tile%gs_max, tile%vpd_min, tile%vpd_max)
!   deallocate( tile%gamma_stom, tile%gamma_soil_c_fac, tile%gamma_soil_default)

end subroutine deallocate_tile

!
! Init some of the tiled variables, in case of cold start.
! Called from modstartup -> readinitfiles()
!
subroutine init_lsm_tiles
    use modfields, only : thlprof, qtprof
    implicit none

    do ilu=1,nlu
      tile(ilu) % thlskin(:,:) = thlprof(1)
      tile(ilu) % qtskin (:,:) = qtprof(1)
      tile(ilu) % obuk   (:,:) = -0.1
    end do

end subroutine init_lsm_tiles

!
! Initialise the LSM homogeneous from namelist input
!
subroutine init_homogeneous
    use modglobal,   only : ifnamopt, fname_options, checknamelisterror, lwarmstart, eps1
    use modmpi,      only : myid, comm3d, mpierr, D_MPI_BCAST
    use modsurfdata, only : tsoil, tsoilm, phiw, phiwm, wl, wlm, wmax
    implicit none

    integer :: ierr, k
    integer :: ilu_lv, ilu_hv, ilu_bs, ilu_aq, ilu_ap, ilu_ws
    real :: c_low, c_high, c_bare, c_water, c_asph
    real :: z0m_low, z0m_high, z0m_bare, z0m_water, z0m_asph
    real :: z0h_low, z0h_high, z0h_bare, z0h_water, z0h_asph
    real :: lambda_s_low, lambda_s_high, lambda_s_bare, lambda_s_asph
    real :: lambda_us_low, lambda_us_high, lambda_us_bare, lambda_us_asph
    real :: lai_low, lai_high
    real :: rs_min_low, rs_min_high, rs_min_bare, rs_min_asph
    real :: ar_low, br_low, ar_high, br_high
    real :: gD_high
    real :: tskin_water

    real, allocatable :: t_soil_p(:), theta_soil_p(:)
    integer, allocatable :: soil_index_p(:)

    ! Read namelist
    namelist /NAMLSM_HOMOGENEOUS/ &
        c_low, c_high, c_bare, c_water, c_asph, &
        z0m_low, z0m_high, z0m_bare, z0m_water, z0m_asph, &
        z0h_low, z0h_high, z0h_bare, z0h_water, z0h_asph, &
        lambda_s_low, lambda_s_high, lambda_s_bare, lambda_s_asph, &
        lambda_us_low, lambda_us_high, lambda_us_bare, lambda_us_asph, &
        lai_low, lai_high, &
        rs_min_low, rs_min_high, rs_min_bare, rs_min_asph, &
        t_soil_p, theta_soil_p, soil_index_p, &
        ar_low, br_low, ar_high, br_high, &
        gD_high, tskin_water

    allocate(t_soil_p(kmax_soil), theta_soil_p(kmax_soil))
    allocate(soil_index_p(kmax_soil))
  !todo : check if len(soil_index_p)==len(t_soil_p)==len(theta_soil_p)==len(kmax_soil)
  !if not, throw error
    if (myid == 0) then
        open(ifnamopt, file=fname_options, status='old', iostat=ierr)
        read(ifnamopt, NAMLSM_HOMOGENEOUS, iostat=ierr)
        call checknamelisterror(ierr, ifnamopt, 'NAMLSM_HOMOGENEOUS')
        write(6, NAMLSM_HOMOGENEOUS)
        close(ifnamopt)
    end if

    ! Broadcast to all MPI tasks
    call D_MPI_BCAST(c_low,   1, 0, comm3d, mpierr)
    call D_MPI_BCAST(c_high,  1, 0, comm3d, mpierr)
    call D_MPI_BCAST(c_bare,  1, 0, comm3d, mpierr)
    call D_MPI_BCAST(c_water, 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(c_asph,  1, 0, comm3d, mpierr)

    call D_MPI_BCAST(z0m_low,   1, 0, comm3d, mpierr)
    call D_MPI_BCAST(z0m_high,  1, 0, comm3d, mpierr)
    call D_MPI_BCAST(z0m_bare,  1, 0, comm3d, mpierr)
    call D_MPI_BCAST(z0m_water, 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(z0m_asph,  1, 0, comm3d, mpierr)

    call D_MPI_BCAST(z0h_low,   1, 0, comm3d, mpierr)
    call D_MPI_BCAST(z0h_high,  1, 0, comm3d, mpierr)
    call D_MPI_BCAST(z0h_bare,  1, 0, comm3d, mpierr)
    call D_MPI_BCAST(z0h_water, 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(z0h_asph,  1, 0, comm3d, mpierr)

    call D_MPI_BCAST(lambda_s_low,  1, 0, comm3d, mpierr)
    call D_MPI_BCAST(lambda_s_high, 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(lambda_s_bare, 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(lambda_s_asph, 1, 0, comm3d, mpierr)

    call D_MPI_BCAST(lambda_us_low,  1, 0, comm3d, mpierr)
    call D_MPI_BCAST(lambda_us_high, 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(lambda_us_bare, 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(lambda_us_asph, 1, 0, comm3d, mpierr)

    call D_MPI_BCAST(lai_low,  1,  0, comm3d, mpierr)
    call D_MPI_BCAST(lai_high, 1,  0, comm3d, mpierr)

    call D_MPI_BCAST(rs_min_low,  1, 0, comm3d, mpierr)
    call D_MPI_BCAST(rs_min_high, 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(rs_min_bare, 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(rs_min_asph, 1, 0, comm3d, mpierr)

    call D_MPI_BCAST(t_soil_p,     kmax_soil, 0, comm3d, mpierr)
    call D_MPI_BCAST(theta_soil_p, kmax_soil, 0, comm3d, mpierr)
    call D_MPI_BCAST(soil_index_p, kmax_soil, 0, comm3d, mpierr)

    call D_MPI_BCAST(ar_low,  1, 0, comm3d, mpierr)
    call D_MPI_BCAST(br_low,  1, 0, comm3d, mpierr)
    call D_MPI_BCAST(ar_high, 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(br_high, 1, 0, comm3d, mpierr)

    call D_MPI_BCAST(gD_high, 1, 0, comm3d, mpierr)

    call D_MPI_BCAST(tskin_water, 1, 0, comm3d, mpierr)

    ! hack to set LU indices
    ilu_lv = 1
    ilu_hv = 2
    ilu_bs = 3
    ilu_aq = 4
    ilu_ap = 5
    ilu_ws = 6

    ! hack to set LU names & type
    ! todo: read from namoptions
    tile(ilu_lv) % lushort = 'lv'
    tile(ilu_hv) % lushort = 'hv'
    tile(ilu_bs) % lushort = 'bs'
    tile(ilu_aq) % lushort = 'aq'
    tile(ilu_ap) % lushort = 'ap'
    tile(ilu_ws) % lushort = 'ws'

    tile(ilu_lv)%luname = 'low vegetation'
    tile(ilu_hv)%luname = 'high vegetation'
    tile(ilu_bs)%luname = 'bare soil'
    tile(ilu_aq)%luname = 'water'
    tile(ilu_ap)%luname = 'asphalt'
    tile(ilu_ws)%luname = 'water on leaf'

    tile(ilu_lv)%lveg = .True.
    tile(ilu_hv)%lveg = .True.
    tile(ilu_bs)%lveg = .False.
    tile(ilu_aq)%lveg = .False.
    tile(ilu_ap)%lveg = .False.
    tile(ilu_ws)%lveg = .False.

    tile(ilu_lv)%laqu = .False.
    tile(ilu_hv)%laqu = .False.
    tile(ilu_bs)%laqu = .False.
    tile(ilu_aq)%laqu = .True.
    tile(ilu_ap)%laqu = .False.
    tile(ilu_ws)%laqu = .False.

    ! Set values
    tile(ilu_lv) % base_frac(:,:) = c_low
    tile(ilu_hv) % base_frac(:,:) = c_high
    tile(ilu_bs) % base_frac(:,:) = c_bare
    tile(ilu_aq) % base_frac(:,:) = c_water
    tile(ilu_ap) % base_frac(:,:) = c_asph
    tile(ilu_ws) % base_frac(:,:) = 0.

    tile(ilu_lv) % z0m(:,:) = z0m_low
    tile(ilu_hv) % z0m(:,:) = z0m_high
    tile(ilu_bs) % z0m(:,:) = z0m_bare
    tile(ilu_aq) % z0m(:,:) = z0m_water
    tile(ilu_ap) % z0m(:,:) = z0m_asph

    tile(ilu_lv) % z0h(:,:) = z0h_low
    tile(ilu_hv) % z0h(:,:) = z0h_high
    tile(ilu_bs) % z0h(:,:) = z0h_bare
    tile(ilu_aq) % z0h(:,:) = z0h_water
    tile(ilu_ap) % z0h(:,:) = z0h_asph

    tile(ilu_lv) % lambda_stable(:,:) = lambda_s_low
    tile(ilu_hv) % lambda_stable(:,:) = lambda_s_high
    tile(ilu_bs) % lambda_stable(:,:) = lambda_s_bare
    tile(ilu_ap) % lambda_stable(:,:) = lambda_s_asph

    tile(ilu_lv) % lambda_unstable(:,:) = lambda_us_low
    tile(ilu_hv) % lambda_unstable(:,:) = lambda_us_high
    tile(ilu_bs) % lambda_unstable(:,:) = lambda_us_bare
    tile(ilu_ap) % lambda_unstable(:,:) = lambda_us_asph

    tile(ilu_lv) % lai(:,:) = lai_low
    tile(ilu_hv) % lai(:,:) = lai_high

    tile(ilu_lv) % rs_min(:,:) = rs_min_low
    tile(ilu_hv) % rs_min(:,:) = rs_min_high
    tile(ilu_bs) % rs_min(:,:) = rs_min_bare
    tile(ilu_ap) % rs_min(:,:) = rs_min_asph

    tile(ilu_lv) % a_r(:,:) = ar_low
    tile(ilu_lv) % b_r(:,:) = br_low
    tile(ilu_hv) % a_r(:,:) = ar_high
    tile(ilu_hv) % b_r(:,:) = br_high

    tile(ilu_hv) % gD(:,:) = gD_high

    tile(ilu_aq) % tskin(:,:) = tskin_water

    if (.not. lwarmstart) then
        ! Init prognostic variables in case of cold start.
        ! For a warm start, these are read from restartfiles
        ! in modstartup.
        do k=1, kmax_soil
            tsoil(:,:,k) = t_soil_p(k)
            phiw (:,:,k) = theta_soil_p(k)
        end do

        tsoilm(:,:,:) = tsoil(:,:,:)
        phiwm (:,:,:) = phiw (:,:,:)

        wl(:,:)  = 0.
        wlm(:,:) = 0.
    end if

    do k=1, kmax_soil
        soil_index(:,:,k) = soil_index_p(k)
    end do

    ! Set properties wet skin tile
    tile(ilu_ws) % z0m(:,:) = c_low*z0m_low + c_high*z0m_high + c_bare*z0m_bare
    tile(ilu_ws) % z0h(:,:) = c_low*z0h_low + c_high*z0h_high + c_bare*z0h_bare

    tile(ilu_ws) % lambda_stable(:,:)   = c_low*lambda_s_low  + c_high*lambda_s_high  + c_bare*lambda_s_bare + c_asph*lambda_s_asph
    tile(ilu_ws) % lambda_unstable(:,:) = c_low*lambda_us_low + c_high*lambda_us_high + c_bare*lambda_us_bare + c_asph*lambda_us_asph

    ! Calculate land fraction, and limit to prevent div/0's
    land_frac(:,:) = 1.-c_water
    where (land_frac == 0) land_frac = eps1

    ! Calculate vegetation fraction, and limit to prevent div/0's
    cveg(:,:) = tile(ilu_lv)%base_frac(:,:) + tile(ilu_hv)%base_frac(:,:)
    where (cveg == 0) cveg = eps1


    ! Max liquid water per grid point, accounting for LAI
    wl_max(:,:) = wmax * ( &
            tile(ilu_lv)%base_frac(:,:) * tile(ilu_lv)%lai(:,:) + &
            tile(ilu_hv)%base_frac(:,:) * tile(ilu_hv)%lai(:,:) + &
            tile(ilu_ap)%base_frac(:,:) + &
            tile(ilu_bs)%base_frac(:,:)) / land_frac(:,:)
    where (wl_max == 0) wl_max = eps1

    ! Cleanup!
    deallocate(t_soil_p, theta_soil_p, soil_index_p)

end subroutine init_homogeneous

!
! Init LSM properties from external input: read from NetCDF
!
subroutine init_heterogeneous_nc

    use netcdf
    use modglobal,   only : i1, j1, lwarmstart, iexpnr, eps1
    use modglobal,   only : i2, j2
    use modmpi,      only : myid, myidx, myidy
    use modglobal,   only : imax, jmax, itot, jtot, ldrydep

    use modsurfdata, only : tsoil, phiw, wl, wlm, wmax
    implicit none

    !integer       :: ilu_lv, ilu_hv, ilu_aq, ilu_ap, ilu_ws, ilu_bs
    integer       :: ilu_ws
    real          :: lai_tmp(i2,j2)
    character(32) :: input_file = 'lsm.inp_xxx.nc'
    integer       :: ncid , varid
    character     :: lvegs(nlu-1)
    character     :: laqus(nlu-1)
    integer       :: nlu_file, len_x, len_y

    write(input_file(9:11), '(i3.3)') iexpnr

    write(6,"(A18, A32)") "Reading LSM input: ", input_file
    call check( nf90_open(input_file, nf90_nowrite, ncid) )

    if (myid==0) then
      call check( nf90_inq_dimid(ncid, 'nlu', varid) )
      call check( nf90_inquire_dimension(ncid, varid, len=nlu_file) )

      ! check if nlu_file==nlu-1 ('wet skin' is not in file)
      if (nlu-1 /= nlu_file) then
        write(6,"(A58, i3, A3, i3)") "STOPPED. Number of LU types in file differs from nlu-1:  ", nlu-1, " /=", nlu_file
        stop
      end if

      ! check if dimensions of lsm.inp_xxx.nc agree with the DALES domain
      call check( nf90_inq_dimid(ncid, 'x', varid) )
      call check( nf90_inquire_dimension(ncid, varid, len=len_x) )
      if (len_x /= itot) then
        write(6,"(A62, i3, A3, i3)") "STOPPED. x-dimension of lsm.inp differs from DALES domain: ", len_x, " /=", itot
        stop
      end if
      call check( nf90_inq_dimid(ncid, 'y', varid) )
      call check( nf90_inquire_dimension(ncid, varid, len=len_y) )
      if (len_y /= jtot) then
        write(6,"(A62, i3, A3, i3)") "STOPPED. y-dimension of lsm.inp differs from DALES domain: ", len_y, " /=", jtot
        stop
      end if
    end if

    ! set land use names and indices, and whether it's vegetation or not
    do ilu=1,nlu-1
      call check( nf90_inq_varid( ncid, 'luname', varid) )
      call check( nf90_get_var(ncid, varid, tile(ilu)%luname , &
                               start = (/ 1 , ilu /) , &
                               count = (/ 32, 1   /)) )

      call check( nf90_inq_varid( ncid, 'lushort', varid) )
      call check( nf90_get_var(ncid, varid, tile(ilu)%lushort , &
                              start = (/ 1, ilu /) , &
                              count = (/ 3, 1   /)) )

      call check( nf90_inq_varid( ncid, 'lveg', varid) )
      call check( nf90_get_var(ncid, varid, lvegs(ilu), &
                                start = (/ 1, ilu /), &
                                count = (/ 1, 1   /) ) )
      if (trim(lvegs(ilu)) .eq. "T") then
        tile(ilu)%lveg = .true.
      else if (trim(lvegs(ilu)) .eq. "F") then
        tile(ilu)%lveg = .false.
      else
        print *,'logical lveg not defined'
        stop
      end if

      call check( nf90_inq_varid( ncid, 'laqu', varid) )
      call check( nf90_get_var(ncid, varid, laqus(ilu), &
                                start = (/ 1, ilu /), &
                                count = (/ 1, 1   /) ) )
      if (trim(laqus(ilu)) .eq. "T") then
        tile(ilu)%laqu = .true.
      else if (trim(laqus(ilu)) .eq. "F") then
        tile(ilu)%laqu = .false.
      else
        print *,'logical laqu not defined'
        stop
      end if

    end do

    ! wet skin is always the last LU type; set explicitly
    ilu_ws            = nlu
    tile(nlu)%luname  = 'wet skin'
    tile(nlu)%lushort = 'ws '
    tile(nlu)%lveg    = .false.
    tile(nlu)%laqu    = .false.

    ! 2D surface fields
    do ilu=1,nlu-1
      write(*,*) 'reading variables for LU type: ', trim(tile(ilu)%lushort)
      ! LU cover
      call check( nf90_inq_varid( ncid, 'cover_'//trim(tile(ilu)%lushort), varid) )
      call check( nf90_get_var(ncid, varid, tile(ilu)%base_frac(2:i1, 2:j1) , &
                                start = (/1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/imax, jmax/) ) )
      ! z0m
      call check( nf90_inq_varid( ncid, 'z0m_'//trim(tile(ilu)%lushort), varid) )
      call check( nf90_get_var(ncid, varid, tile(ilu)%z0m(2:i1, 2:j1) , &
                                start = (/1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/imax, jmax/) ) )
      ! z0h
      call check( nf90_inq_varid( ncid, 'z0h_'//trim(tile(ilu)%lushort), varid) )
      call check( nf90_get_var(ncid, varid, tile(ilu)%z0h(2:i1, 2:j1) , &
                                start = (/1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/imax, jmax/) ) )
      ! lambda stable
      call check( nf90_inq_varid( ncid, 'lambda_s_'//trim(tile(ilu)%lushort), varid) )
      call check( nf90_get_var(ncid, varid, tile(ilu)%lambda_stable(2:i1, 2:j1) , &
                                start = (/1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/imax, jmax/) ) )
      ! lambda unstable
      call check( nf90_inq_varid( ncid, 'lambda_us_'//trim(tile(ilu)%lushort), varid) )
      call check( nf90_get_var(ncid, varid, tile(ilu)%lambda_unstable(2:i1, 2:j1) , &
                                start = (/1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/imax, jmax/) ) )
      ! rs_min
      call check( nf90_inq_varid( ncid, 'rs_min_'//trim(tile(ilu)%lushort), varid) )
      call check( nf90_get_var(ncid, varid, tile(ilu)%rs_min(2:i1, 2:j1) , &
                                start = (/1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/imax, jmax/) ) )
      ! lai
      call check( nf90_inq_varid( ncid, 'lai_'//trim(tile(ilu)%lushort), varid) )
      call check( nf90_get_var(ncid, varid, tile(ilu)%lai(2:i1, 2:j1) , &
                                start = (/1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/imax, jmax/) ) )
      ! a_r
      call check( nf90_inq_varid( ncid, 'ar_'//trim(tile(ilu)%lushort), varid) )
      call check( nf90_get_var(ncid, varid, tile(ilu)%a_r(2:i1, 2:j1) , &
                                start = (/1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/imax, jmax/) ) )
      ! b_r
      call check( nf90_inq_varid( ncid, 'br_'//trim(tile(ilu)%lushort), varid) )
      call check( nf90_get_var(ncid, varid, tile(ilu)%b_r(2:i1, 2:j1) , &
                                start = (/1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/imax, jmax/) ) )
      ! gD
      call check( nf90_inq_varid( ncid, 'gD_'//trim(tile(ilu)%lushort), varid) )
      call check( nf90_get_var(ncid, varid, tile(ilu)%gD(2:i1, 2:j1) , &
                                start = (/1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/imax, jmax/) ) )
    ! tskin
      call check( nf90_inq_varid( ncid, 'tskin_'//trim(tile(ilu)%lushort), varid) )
      call check( nf90_get_var(ncid, varid, tile(ilu)%tskin(2:i1, 2:j1) , &
                              start = (/1 + myidx * imax, 1 + myidy * jmax/), &
                              count = (/imax, jmax/) ) )

    !!! deposition parameters
    if (ldrydep) then
      ! ! R_inc_b
      ! call check( nf90_inq_varid( ncid, 'R_inc_b_'//trim(tile(ilu)%lushort), varid) )
      ! call check( nf90_get_var(ncid, varid, tile(ilu)%R_inc_b(2:i1, 2:j1) , &
      !                           start = (/1 + myidx * imax, 1 + myidy * jmax/), &
      !                           count = (/imax, jmax/) ) )
      ! ! R_inc_h
      ! call check( nf90_inq_varid( ncid, 'R_inc_h_'//trim(tile(ilu)%lushort), varid) )
      ! call check( nf90_get_var(ncid, varid, tile(ilu)%R_inc_h(2:i1, 2:j1) , &
      !                           start = (/1 + myidx * imax, 1 + myidy * jmax/), &
      !                           count = (/imax, jmax/) ) )
      ! SAI_a
      call check( nf90_inq_varid( ncid, 'SAI_a_'//trim(tile(ilu)%lushort), varid) )
      call check( nf90_get_var(ncid, varid, tile(ilu)%SAI_a(2:i1, 2:j1) , &
                                start = (/1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/imax, jmax/) ) )
      ! SAI_b
      call check( nf90_inq_varid( ncid, 'SAI_b_'//trim(tile(ilu)%lushort), varid) )
      call check( nf90_get_var(ncid, varid, tile(ilu)%SAI_b(2:i1, 2:j1) , &
                                start = (/1 + myidx * imax, 1 + myidy * jmax/), &
                                count = (/imax, jmax/) ) )
     !  ! fmin
     !  call check( nf90_inq_varid( ncid, 'fmin_'//trim(tile(ilu)%lushort), varid) )
     !  call check( nf90_get_var(ncid, varid, tile(ilu)%fmin(2:i1, 2:j1) , &
     !                            start = (/1 + myidx * imax, 1 + myidy * jmax/), &
     !                            count = (/imax, jmax/) ) )
     !  ! alpha
     !  call check( nf90_inq_varid( ncid, 'alpha_'//trim(tile(ilu)%lushort), varid) )
     !  call check( nf90_get_var(ncid, varid, tile(ilu)%alpha(2:i1, 2:j1) , &
     !                            start = (/1 + myidx * imax, 1 + myidy * jmax/), &
     !                            count = (/imax, jmax/) ) )
     !  ! Tmin
     !  call check( nf90_inq_varid( ncid, 'Tmin_'//trim(tile(ilu)%lushort), varid) )
     !  call check( nf90_get_var(ncid, varid, tile(ilu)%Tmin(2:i1, 2:j1) , &
     !                            start = (/1 + myidx * imax, 1 + myidy * jmax/), &
     !                            count = (/imax, jmax/) ) )
     !  ! Topt
     !  call check( nf90_inq_varid( ncid, 'Topt_'//trim(tile(ilu)%lushort), varid) )
     !  call check( nf90_get_var(ncid, varid, tile(ilu)%Topt(2:i1, 2:j1) , &
     !                            start = (/1 + myidx * imax, 1 + myidy * jmax/), &
     !                            count = (/imax, jmax/) ) )
     !  ! Tmax
     !  call check( nf90_inq_varid( ncid, 'Tmax_'//trim(tile(ilu)%lushort), varid) )
     !  call check( nf90_get_var(ncid, varid, tile(ilu)%Tmax(2:i1, 2:j1) , &
     !                            start = (/1 + myidx * imax, 1 + myidy * jmax/), &
     !                            count = (/imax, jmax/) ) )
     ! ! gs_max
     !  call check( nf90_inq_varid( ncid, 'gs_max_'//trim(tile(ilu)%lushort), varid) )
     !  call check( nf90_get_var(ncid, varid, tile(ilu)%gs_max(2:i1, 2:j1) , &
     !                            start = (/1 + myidx * imax, 1 + myidy * jmax/), &
     !                            count = (/imax, jmax/) ) )
     ! ! vpd_min
     !  call check( nf90_inq_varid( ncid, 'vpd_min_'//trim(tile(ilu)%lushort), varid) )
     !  call check( nf90_get_var(ncid, varid, tile(ilu)%vpd_min(2:i1, 2:j1) , &
     !                            start = (/1 + myidx * imax, 1 + myidy * jmax/), &
     !                            count = (/imax, jmax/) ) )
     ! ! vpd_max
     !  call check( nf90_inq_varid( ncid, 'vpd_max_'//trim(tile(ilu)%lushort), varid) )
     !  call check( nf90_get_var(ncid, varid, tile(ilu)%vpd_max(2:i1, 2:j1) , &
     !                            start = (/1 + myidx * imax, 1 + myidy * jmax/), &
     !                            count = (/imax, jmax/) ) )
     ! ! gamma_stom
     !  call check( nf90_inq_varid( ncid, 'gamma_stom_'//trim(tile(ilu)%lushort), varid) )
     !  call check( nf90_get_var(ncid, varid, tile(ilu)%gamma_stom(2:i1, 2:j1) , &
     !                            start = (/1 + myidx * imax, 1 + myidy * jmax/), &
     !                            count = (/imax, jmax/) ) )
     ! ! gamma_soil_c_fac
     !  call check( nf90_inq_varid( ncid, 'gamma_soil_c_fac_'//trim(tile(ilu)%lushort), varid) )
     !  call check( nf90_get_var(ncid, varid, tile(ilu)%gamma_soil_c_fac(2:i1, 2:j1) , &
     !                            start = (/1 + myidx * imax, 1 + myidy * jmax/), &
     !                            count = (/imax, jmax/) ) )
     ! ! gamma_soil_default
     !  call check( nf90_inq_varid( ncid, 'gamma_soil_default_'//trim(tile(ilu)%lushort), varid) )
     !  call check( nf90_get_var(ncid, varid, tile(ilu)%gamma_soil_default(2:i1, 2:j1) , &
     !                            start = (/1 + myidx * imax, 1 + myidy * jmax/), &
     !                            count = (/imax, jmax/) ))
      end if
    end do

    ! 3D soil fields
    ! soil index
    call check( nf90_inq_varid( ncid, 'index_soil', varid) )
    call check( nf90_get_var(ncid, varid, soil_index(2:i1, 2:j1, 1:kmax_soil) , &
                              start = (/1 + myidx * imax, 1 + myidy * jmax, 1/), &
                              count = (/imax, jmax, kmax_soil/) ) )

    if (.not. lwarmstart) then
      ! soil temperature
      call check( nf90_inq_varid( ncid, 't_soil', varid) )
      call check( nf90_get_var(ncid, varid, tsoil(2:i1, 2:j1, 1:kmax_soil) , &
                                start = (/1 + myidx * imax, 1 + myidy * jmax, 1/), &
                                count = (/imax, jmax, kmax_soil/) ) )

      ! soil moisture
      call check( nf90_inq_varid( ncid, 'theta_soil', varid) )
      call check( nf90_get_var(ncid, varid, phiw(2:i1, 2:j1, 1:kmax_soil) , &
                                start = (/1 + myidx * imax, 1 + myidy * jmax, 1/), &
                                count = (/imax, jmax, kmax_soil/) ) )

      ! liquid water
      wl(:,:)  = 0.
      wlm(:,:) = 0.
    end if

    call check( nf90_close(ncid) )

    ! Derived quantities
    tile(ilu_ws)%base_frac(:,:) = 0.

    ! Set properties wet skin tile
    tile(ilu_ws)%z0m(:,:) = 0
    tile(ilu_ws)%z0h(:,:) = 0
    tile(ilu_ws)%lambda_stable(:,:) = 0
    tile(ilu_ws)%lambda_unstable(:,:) = 0

    do ilu=1,nlu
      if (tile(ilu)%laqu) then
        cycle
      else
        tile(ilu_ws)%z0m(:,:) = tile(ilu_ws)%z0m(:,:) + tile(ilu)%base_frac(:,:)*tile(ilu)%z0m(:,:)
        tile(ilu_ws)%z0h(:,:) = tile(ilu_ws)%z0h(:,:) + tile(ilu)%base_frac(:,:)*tile(ilu)%z0h(:,:)
        tile(ilu_ws)%lambda_stable(:,:) = tile(ilu_ws)%lambda_stable(:,:) + tile(ilu)%base_frac(:,:)*tile(ilu)%lambda_stable(:,:)
        tile(ilu_ws)%lambda_unstable(:,:) = tile(ilu_ws)%lambda_unstable(:,:) + tile(ilu)%base_frac(:,:)*tile(ilu)%lambda_unstable(:,:)
      end if
    end do

    ! Calculate vegetation fraction, and limit to prevent div/0's
    do ilu=1,nlu
      if (tile(ilu)%lveg) then
        cveg(:,:) = cveg(:,:) + tile(ilu)%base_frac(:,:)
      end if
    end do
    where (cveg == 0) cveg = eps1

    ! Calculate land fraction, and limit to prevent div/0's
    do ilu=1,nlu
      if (tile(ilu)%laqu) then
        land_frac(:,:) = 1.-tile(ilu)%base_frac(:,:)
      end if
    end do
    where (land_frac == 0) land_frac = eps1

    ! Max liquid water per grid point, accounting for LAI
    wl_max = 0
    do ilu=1,nlu
      if (tile(ilu)%laqu) then
        cycle
      endif
      if (tile(ilu)%lveg) then
        lai_tmp(:,:) = tile(ilu)%lai(:,:)
        wl_max(:,:)  = wl_max(:,:) + tile(ilu)%base_frac(:,:) * lai_tmp(:,:)
      else
        wl_max(:,:)  = wl_max(:,:) + tile(ilu)%base_frac(:,:)
      endif
    enddo
    wl_max(:,:) = wl_max(:,:) * wmax/land_frac(:,:)
    where (wl_max == 0) wl_max = eps1

    ! !! debugging: some checks
    ! write(*,*) '...checking LU inputs: mean, min and max, accounting for ghost cells'
    ! do ilu=1,nlu
    !   write(*,*) 'lushort       ', tile(ilu)%lushort
    !   write(*,*) 'base_frac     ', sum(tile(ilu)%base_frac(2:i1,2:j1))/size(tile(ilu)%base_frac(2:i1,2:j1)), minval(tile(ilu)%base_frac(2:i1,2:j1)), maxval(tile(ilu)%base_frac(2:i1,2:j1))
    !   write(*,*) 'z0h           ', sum(tile(ilu)%z0h(2:i1,2:j1))/size(tile(ilu)%z0h(2:i1,2:j1)), minval(tile(ilu)%z0h(2:i1,2:j1)), maxval(tile(ilu)%z0h(2:i1,2:j1))
    !   write(*,*) 'z0m           ', sum(tile(ilu)%z0m(2:i1,2:j1))/size(tile(ilu)%z0m(2:i1,2:j1)), minval(tile(ilu)%z0m(2:i1,2:j1)), maxval(tile(ilu)%z0m(2:i1,2:j1))
    !   write(*,*) 'rs_min        ', sum(tile(ilu)%rs_min(2:i1,2:j1))/size(tile(ilu)%rs_min(2:i1,2:j1)), minval(tile(ilu)%rs_min(2:i1,2:j1)), maxval(tile(ilu)%rs_min(2:i1,2:j1))
    !   write(*,*) 'gD            ', sum(tile(ilu)%gD(2:i1,2:j1))/size(tile(ilu)%gD(2:i1,2:j1)), minval(tile(ilu)%gD(2:i1,2:j1)), maxval(tile(ilu)%gD(2:i1,2:j1))

    !   write(*,*) 'lambda_stable ', sum(tile(ilu)%lambda_stable(2:i1,2:j1))/size(tile(ilu)%lambda_stable(2:i1,2:j1)), minval(tile(ilu)%lambda_stable(2:i1,2:j1)), maxval(tile(ilu)%lambda_stable(2:i1,2:j1))
    !   write(*,*) 'lambda_unstab ', sum(tile(ilu)%lambda_unstable(2:i1,2:j1))/size(tile(ilu)%lambda_unstable(2:i1,2:j1)), minval(tile(ilu)%lambda_unstable(2:i1,2:j1)), maxval(tile(ilu)%lambda_unstable(2:i1,2:j1))

    !   write(*,*) 'lai           ', sum(tile(ilu)%lai(2:i1,2:j1))/size(tile(ilu)%lai(2:i1,2:j1)), minval(tile(ilu)%lai(2:i1,2:j1)), maxval(tile(ilu)%lai(2:i1,2:j1))
    !   write(*,*) 'a_r           ', sum(tile(ilu)%a_r(2:i1,2:j1))/size(tile(ilu)%a_r(2:i1,2:j1)), minval(tile(ilu)%a_r(2:i1,2:j1)), maxval(tile(ilu)%a_r(2:i1,2:j1))
    !   write(*,*) 'b_r           ', sum(tile(ilu)%b_r(2:i1,2:j1))/size(tile(ilu)%b_r(2:i1,2:j1)), minval(tile(ilu)%b_r(2:i1,2:j1)), maxval(tile(ilu)%b_r(2:i1,2:j1))
    !   write(*,*) 'tskin         ', sum(tile(ilu)%tskin(2:i1,2:j1))/size(tile(ilu)%tskin(2:i1,2:j1)), minval(tile(ilu)%tskin(2:i1,2:j1)), maxval(tile(ilu)%tskin(2:i1,2:j1))
    !   write(*,*) 'R_inc_b       ', sum(tile(ilu)%R_inc_b(2:i1,2:j1))/size(tile(ilu)%R_inc_b(2:i1,2:j1)), minval(tile(ilu)%R_inc_b(2:i1,2:j1)), maxval(tile(ilu)%R_inc_b(2:i1,2:j1))
    !   write(*,*) ''
    ! end do

    ! write(*,*) 'soil_index top', sum(soil_index(2:i1,2:j1,4))/size(soil_index(2:i1,2:j1,4)), minval(soil_index(2:i1,2:j1,4)), maxval(soil_index(2:i1,2:j1,4))

    ! write(*,*) 'tsoil top     ', sum(tsoil(2:i1,2:j1,4))/size(tsoil(2:i1,2:j1,4)), minval(tsoil(2:i1,2:j1,4)), maxval(tsoil(2:i1,2:j1,4))
    ! write(*,*) 'phiw top      ', sum(phiw(2:i1,2:j1,4))/size(phiw(2:i1,2:j1,4)), minval(phiw(2:i1,2:j1,4)), maxval(phiw(2:i1,2:j1,4))
    ! write(*,*) 'cveg     ', sum(cveg)/size(cveg), minval(cveg), maxval(cveg)
    ! write(*,*) 'land_frac', sum(land_frac)/size(land_frac), minval(land_frac), maxval(land_frac)
    ! write(*,*) 'wl_max   ', sum(wl_max)/size(wl_max), minval(wl_max), maxval(wl_max)
    ! write(*,*) 'wmax     ', wmax
    ! !call flush()

end subroutine init_heterogeneous_nc

!
! Read the input table with the (van Genuchten) soil parameters
!
subroutine read_soil_table
    use modmpi, only : myid, comm3d, mpierr, D_MPI_BCAST
    implicit none
    integer :: table_size, ncid, dimid, varid

    if (myid == 0) then
        ! Open the NetCDF file and read the table size
        print*,'Reading "van_genuchten_parameters.nc"'
        call check( nf90_open('van_genuchten_parameters.nc', nf90_nowrite, ncid) )
        call check( nf90_inq_dimid(ncid, 'index', dimid) )
        call check( nf90_inquire_dimension(ncid, dimid, len=table_size) )
    end if

    call D_MPI_BCAST(table_size, 1, 0, comm3d, mpierr)

    ! Allocate variables
    allocate( theta_res(table_size), theta_wp(table_size), theta_fc(table_size) )
    allocate( theta_sat(table_size), gamma_theta_sat(table_size) )
    allocate( vg_a(table_size), vg_l(table_size), vg_n(table_size) )

    if (myid == 0) then
        ! Read variables
        call check( nf90_inq_varid(ncid, 'theta_res', varid) )
        call check( nf90_get_var(ncid, varid, theta_res) )

        call check( nf90_inq_varid(ncid, 'theta_wp', varid) )
        call check( nf90_get_var(ncid, varid, theta_wp) )

        call check( nf90_inq_varid(ncid, 'theta_fc', varid) )
        call check( nf90_get_var(ncid, varid, theta_fc) )

        call check( nf90_inq_varid(ncid, 'theta_sat', varid) )
        call check( nf90_get_var(ncid, varid, theta_sat) )

        call check( nf90_inq_varid(ncid, 'gamma_sat', varid) )
        call check( nf90_get_var(ncid, varid, gamma_theta_sat) )

        call check( nf90_inq_varid(ncid, 'alpha', varid) )
        call check( nf90_get_var(ncid, varid, vg_a) )

        call check( nf90_inq_varid(ncid, 'l', varid) )
        call check( nf90_get_var(ncid, varid, vg_l) )

        call check( nf90_inq_varid(ncid, 'n', varid) )
        call check( nf90_get_var(ncid, varid, vg_n) )

        call check( nf90_close(ncid) )
    end if

    ! Broadcast to other MPI tasks
    call D_MPI_BCAST(theta_res,       table_size, 0, comm3d, mpierr)
    call D_MPI_BCAST(theta_wp,        table_size, 0, comm3d, mpierr)
    call D_MPI_BCAST(theta_fc,        table_size, 0, comm3d, mpierr)
    call D_MPI_BCAST(theta_sat,       table_size, 0, comm3d, mpierr)
    call D_MPI_BCAST(gamma_theta_sat, table_size, 0, comm3d, mpierr)
    call D_MPI_BCAST(vg_a,            table_size, 0, comm3d, mpierr)
    call D_MPI_BCAST(vg_l,            table_size, 0, comm3d, mpierr)
    call D_MPI_BCAST(vg_n,            table_size, 0, comm3d, mpierr)

end subroutine read_soil_table

!
! Calculate derived (tabulated) soil properties
!
subroutine calc_soil_properties
    use modglobal, only : rho_solid_soil, rho_c_matrix, rho_c_water, eps1
    implicit none

    integer :: i, table_size
    real :: theta_norm_min, theta_norm_max, rho_dry

    table_size = size(vg_n)
    allocate(vg_m(table_size))
    allocate(lambda_theta_min(table_size))
    allocate(lambda_theta_max(table_size))
    allocate(gamma_theta_min(table_size))
    allocate(gamma_theta_max(table_size))
    allocate(gamma_t_dry(table_size))
    allocate(rho_C(table_size))

    do i=1, table_size
        ! van Genuchten parameter `m`
        vg_m(i) = (1. - (1. / vg_n(i)))

        ! Min/max values diffusivity soil moisture
        theta_norm_min = (1.001 * theta_res(i) - theta_res(i)) / (theta_sat(i) - theta_res(i)) + eps1
        theta_norm_max = (0.999 * theta_sat(i) - theta_res(i)) / (theta_sat(i) - theta_res(i))

        lambda_theta_min(i) = calc_diffusivity_vg( &
                theta_norm_min, vg_a(i), vg_l(i), vg_m(i), gamma_theta_sat(i), &
                theta_sat(i), theta_res(i))
        lambda_theta_max(i) = calc_diffusivity_vg( &
                theta_norm_max, vg_a(i), vg_l(i), vg_m(i), gamma_theta_sat(i), &
                theta_sat(i), theta_res(i))

        ! Min/max values conductivity soil moisture
        gamma_theta_min(i) = 0.
        gamma_theta_max(i) = gamma_theta_sat(i)

        ! Conductivity temperature
        rho_dry = (1. - theta_sat(i)) * rho_solid_soil  ! Density of soil (kg m-3)
        gamma_t_dry(i) = (0.135 * rho_dry + 64.7) / (rho_solid_soil - 0.947 * rho_dry)
        rho_C(i) = (1. - theta_sat(i)) * rho_C_matrix + theta_fc(i) * rho_C_water
    end do

end subroutine calc_soil_properties

!
! Calculate root fraction for the low and high vegetation tiles
!
subroutine calc_root_fractions
    use modglobal, only : i1, j1
    implicit none
    real :: root_sum(nlu)
    integer i, j, k

    do j=2,j1
      do i=2,i1
        root_sum = 0
        do k=2, kmax_soil
          do ilu=1,nlu
            if (tile(ilu)%lveg) then
              tile(ilu)%root_frac(i,j,k) = 0.5 * (&
                  exp(tile(ilu)%a_r(i,j) * zh_soil(k+1)) + &
                  exp(tile(ilu)%b_r(i,j) * zh_soil(k+1)) - &
                  exp(tile(ilu)%a_r(i,j) * zh_soil(k  )) - &
                  exp(tile(ilu)%b_r(i,j) * zh_soil(k  )))
              root_sum(ilu) = root_sum(ilu) + tile(ilu)%root_frac(i,j,k)
            end if
          end do
        end do

        ! Make sure that the fractions sum to one.
        do ilu=1,nlu
          if (tile(ilu)%lveg) then
            tile(ilu)%root_frac(i,j,1) = 1. - root_sum(ilu)
          end if
        end do
      end do
    end do

end subroutine calc_root_fractions

!
! Iterative Rib -> Obukhov length solver
!
function calc_obuk_dirichlet(L_in, du, db_in, zsl, z0m, z0h) result(res)
    use modglobal, only : fkar
    implicit none
    real, intent(in) :: L_in, du, db_in, zsl, z0m, z0h

    integer :: m, n, nlim
    real :: res, L, db, Lmax, Ri, L0, Lstart, Lend, fx0, fxdif

    m = 0
    nlim = 10
    Lmax = 1e10
    L = L_in
    db = db_in

    ! The solver does not have a solution for large Ri numbers,
    ! i.e. the `fx` equation below has no zero crossing.
    ! The limit of 0.13 typically results in a minimum (positive) L of ~1.
    ! IFS caps z/L at 5, so for a typical zsl at ~L=2.
    !Ri = fkar * db * zsl / du**2
    !if (Ri > 0.13) then
    !    print*,'WARNING: Ri out of range, returning L=1'
    !    res = 1
    !    return
    !end if

    ! Avoid buoyancy difference of zero:
    if (db >= 0) then
        db = max(db, 1e-9)
    else
        db = min(db, -1e-9)
    end if

    ! Allow for one restart of iterative procedure:
    do while (m <= 1)
        ! if L and db are of different sign, or the last calculation did not converge,
        ! the stability has changed and the procedure needs to be reset
        if (L*db <= 0) then
            nlim = 200
            if (db >= 0) then
                L = 1e-9
            else
                L = -1e-9
            end if
        end if

        ! Make sure the iteration starts
        if (db >= 0) then
            L0 = 1e30
        else
            L0 = -1e30
        end if

        ! Exit on convergence or on iteration count
        n = 0
        fxdif = 1
        do while (abs((L - L0) / L0) > 0.001 .and. n < nlim .and. abs(L) < Lmax)
            L0     = L
            Lstart = L - 0.001*L
            Lend   = L + 0.001*L

            fx0    = fx(zsl, L0, du, db, z0h, z0m)
            fxdif  = (fx(zsl, Lend, du, db, z0h, z0m) - fx(zsl, Lstart, du, db, z0h, z0m)) / (Lend - Lstart)
            L      = L - fx0/fxdif
            n      = n+1
        end do

        if (n < nlim .and. abs(L) < Lmax) then
            ! Convergence has been reached
            res = L
            return
        else
            ! Convergence has not been reached, procedure restarted once
            L = 1e-9
            m = m+1
            nlim = 200
        end if
    end do

    if (m > 1) then
        print*,'WARNING: convergence has not been reached in Obukhov length iteration'
        print*,'Input: ', L_in, du, db_in, zsl, z0m, z0h
        !stop
        res = 1e-9
        return
    end if

end function calc_obuk_dirichlet

pure function fx(zsl, L, du, db, z0h, z0m) result(res)
    implicit none
    real, intent(in) :: zsl, L, du, db, z0h, z0m
    real :: res, fkar
    fkar = 0.4

    res = zsl/L - fkar * zsl * db * fh(zsl, z0h, L) / (du * fm(zsl, z0m, L))**2
end function fx

pure function fm(zsl, z0m, L) result(res)
    use modglobal, only : fkar
    use modsurface, only : psim
    implicit none
    real, intent(in) :: zsl, z0m, L
    real :: res

    res = fkar / (log(zsl/z0m) - psim(zsl/L) + psim(z0m/L))
end function fm

pure function fh(zsl, z0h, L) result(res)
    use modglobal, only : fkar
    use modsurface, only : psih
    implicit none
    real, intent(in) :: zsl, z0h, L
    real :: res

    res = fkar / (log(zsl/z0h) - psih(zsl/L) + psih(z0h/L))
end function fh

!
! Convert soil hydraulic head to soil water content, using van Genuchten parameterisation.
!
pure function psi_to_theta(theta_res, theta_sat, vg_a, vg_n, vg_m, psi) result(res)
    implicit none
    real, intent(in) :: theta_res, theta_sat, vg_a, vg_n, vg_m, psi
    real :: res

    res = theta_res + (theta_sat - theta_res) * (1. / (1.+ abs(vg_a * psi)**vg_n))**vg_m
end function psi_to_theta

!
! Convert soil water content to hydraulic head, using van Genuchten parameterisation.
!
pure function theta_to_psi(theta_res, theta_sat, vg_a, vg_n, vg_m, theta) result(res)
    implicit none
    real, intent(in) :: theta_res, theta_sat, vg_a, vg_n, vg_m, theta
    real :: res

    res = -(vg_a**(-vg_n) * (-1. + ((theta_res - theta_sat)/(theta_res - theta))**(1./vg_m)))**(1./vg_n)
end function theta_to_psi

!
! Calculate hydraulic diffusivity using van Genuchten parameterisation.
!
pure function calc_diffusivity_vg( &
        theta_norm, vg_a, vg_l, vg_m, lambda_sat, theta_sat, theta_res) result(res)
    implicit none
    real, intent(in) :: theta_norm, vg_a, vg_l, vg_m, lambda_sat, theta_sat, theta_res
    real :: res

    res = (1.-vg_m)*lambda_sat / (vg_a * vg_m * (theta_sat-theta_res)) * theta_norm**(vg_l-(1./vg_m)) * &
             (  (1.-theta_norm**(1./vg_m))**(-vg_m) + (1.-theta_norm**(1./vg_m))**vg_m - 2. )
end function calc_diffusivity_vg

!
! Calculate hydraulic conductivity using van Genuchten parameterisation.
!
pure function calc_conductivity_vg(theta_norm, vg_l, vg_m, gamma_sat) result(res)
    implicit none
    real, intent(in) :: theta_norm, vg_l, vg_m, gamma_sat
    real :: res

    res = gamma_sat * theta_norm**vg_l * ( 1.- (1.-theta_norm**(1./vg_m))**vg_m )**2.

end function calc_conductivity_vg

!
! Check NetCDF calls
!
subroutine check(status)
    integer, intent (in) :: status
    if(status /= nf90_noerr) then
        print *,'NetCDF error: ', trim(nf90_strerror(status))
        call abort
    end if
end subroutine check

end module modlsm
