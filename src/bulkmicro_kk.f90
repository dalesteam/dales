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
!  Copyright 1993-2026 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!
!> Kernels for Khairoutdinov-Kogan microphysics.
module bulkmicro_kk
  use modglobal,         only: i1, ih, j1, jh, k1, rlv, cp, pi, rv, pirhow
  use modmicrodata,      only: Nc_0, delt
  use modbulkmicro_data, only: qrmin, qcmin, l_mur_cst, mur_cst
  use modprecision,      only: field_r
  use modtimer,          only: timer_tic, timer_toc

  implicit none

  private

  character(len=*), parameter :: modname = 'bulkmicro_kk'

  public :: autoconversion_kk
  public :: accretion_kk
  public :: evaporation_kk
  public :: sedimentation_rain_kk
  public :: calc_sed_qr_kk
  public :: calc_sed_nr_kk
  public :: xrmin, xrmax

  real(field_r), parameter :: &
    c_evap = 0.87,  & !< Coefficient for evaporation.
    D0 = 50e-6,     & !< Diameter separating cloud and precipitation parts of the DSD.
    Dv = 2.4e-5,    & !< Diffusivity of water vapor [m2/s].
    Kt = 2.5e-2,    & !< Conductivity of heat [J/(sKm)].
    wfallmax = 9.9, & !< Terminal fall velocity.
    xrmin = 0.0,    &
    xrmax = 5.2e-7, & !< Max mean mass of pw.
    eps = 1e-18

contains

  include 'microphysics.inc'

  !> Calculate the autoconversion term.
  !!
  !! \param ql0 Liquid water mixing ratio.
  !! \param rhof Density at full levels.
  !! \param exnf Exner function at full levels.
  !! \param qcbase Lowest level with cloud.
  !! \param qcroof Highest level with cloud.
  !! \param thlpmcr Tendency of $\theta_l$.
  !! \param qtpmcr Tendency of $\q_t$.
  !! \param qrp Tendency of rain water mixing ratio.
  !! \param Nrp Tendency of rain drop number concentration.
  subroutine autoconversion_kk(ql0, Nc, rhof, exnf, qcbase, qcroof, thlpmcr, &
                            qtpmcr, qrp, Nrp, Ncp)
    real(field_r), intent(in)    :: ql0(2-ih:i1+ih,2-jh:j1+jh,1:k1)
    real(field_r), intent(in)    :: nc(2:,2:,:)
    real(field_r), intent(in)    :: rhof(1:k1)
    real(field_r), intent(in)    :: exnf(1:k1)

    integer,       intent(in)    :: qcbase, qcroof

    real(field_r), intent(inout) :: thlpmcr(2:i1,2:j1,1:k1)
    real(field_r), intent(inout) :: qtpmcr(2-ih:i1+ih,2-jh:j1+jh,1:k1)
    real(field_r), intent(inout) :: qrp(2:i1,2:j1,1:k1)
    real(field_r), intent(inout) :: Nrp(2:i1,2:j1,1:k1)

    real(field_r), intent(inout), optional :: Ncp(2:i1,2:j1,1:k1)

    character(len=*), parameter :: routine = modname//'/autoconversion_kk'

    integer       :: i, j, k
    real(field_r) :: &
      au, &
      xc

    call timer_tic(routine, 1)

    if (qcbase > qcroof) then
      call timer_toc(routine)
      return
    end if

    !$acc parallel loop collapse(3) default(present) private(au)
!!$omp target teams loop private(au) collapse(3)&
!!$omp defaultmap(present:aggregate) defaultmap(present:allocatable)
    do k = qcbase, qcroof
      do j = 2, j1
        do i = 2, i1
           if (ql0(i,j,k) > qcmin) then
              au = 1350 * ql0(i,j,k)**(2.47_field_r) &
                   * (nc(i,j,k) / 1E6)**(-1.79_field_r)
              au = min(ql0(i,j,k) / delt, au)
              qrp(i,j,k) = qrp(i,j,k) + au
              qtpmcr(i,j,k) = qtpmcr(i,j,k) - au
              thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlv / (cp * exnf(k))) * au
              Nrp(i,j,k) = Nrp(i,j,k) + au * rhof(k) / (pirhow * D0**3)
              if (present(Ncp)) then
                xc = rhof(k) * ql0(i,j,k) / (Nc(i,j,k) + eps)
                Ncp(i,j,k) = Ncp(i,j,k) - au / xc * rhof(k)
              end if
           endif
        enddo
      enddo
    enddo

    call timer_toc(routine)

  end subroutine autoconversion_kk

  !> Calculate the accretion term.
  !!
  !! \param ql0 Liquid water mixing ratio.
  !! \param qr Rain water mixing ratio.
  !! \param exnf Exner function at full levels.
  !! \param qcbase Lowest level with cloud.
  !! \param qcroof Highest level with cloud.
  !! \param qrbase Lowest level with rain.
  !! \param qrroof Highest level with rain.
  !! \param thlpmcr Tendency of $\theta_l$.
  !! \param qtpmcr Tendency of total water mixing ratio.
  !! \param qrp Tendency of rain water mixing ratio.
  subroutine accretion_kk(ql0, Nc, qr, rhof, exnf, qcbase, qcroof, qrbase, qrroof, &
                          thlpmcr, qtpmcr, qrp, Ncp)
    real(field_r), intent(in)    :: ql0(2-ih:i1+ih,2-jh:j1+jh,1:k1)
    real(field_r), intent(in)    :: Nc(2:i1,2:j1,1:k1)
    real(field_r), intent(in)    :: qr(2:i1,2:j1,1:k1)
    real(field_r), intent(in)    :: rhof(1:k1)
    real(field_r), intent(in)    :: exnf(1:k1)

    integer,       intent(in)    :: qcbase, qcroof
    integer,       intent(in)    :: qrbase, qrroof

    real(field_r), intent(inout) :: thlpmcr(2:i1,2:j1,1:k1)
    real(field_r), intent(inout) :: qtpmcr(2-ih:i1+ih,2-jh:j1+jh,1:k1)
    real(field_r), intent(inout) :: qrp(2:i1,2:j1,1:k1)

    real(field_r), intent(inout), optional :: Ncp(2:i1,2:j1,1:k1)

    character(len=*), parameter :: routine = modname//'/accretion_kk'

    integer       :: i, j, k
    real(field_r) :: &
      ac, &
      xc

    call timer_tic(routine, 1)

    if (max(qrbase, qcbase) > min(qcroof, qcroof)) then
      call timer_toc(routine)
      return
    end if

    !$acc parallel loop collapse(3) default(present) private(ac)
!!$omp target teams loop private(ac) collapse(3)&
!!$omp defaultmap(present:aggregate) defaultmap(present:allocatable)
    do k = max(qrbase, qcbase), min(qcroof, qrroof)
      do j = 2, j1
        do i = 2, i1
          if (ql0(i,j,k) > qcmin .and. qr(i,j,k) > qrmin) then
            ac = 67 * (ql0(i,j,k) * qr(i,j,k))**1.15_field_r
            qrp(i,j,k) = qrp(i,j,k) + ac
            qtpmcr(i,j,k) = qtpmcr(i,j,k) - ac
            thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlv / (cp * exnf(k))) * ac
            if (present(Ncp)) then
              xc = rhof(k) * ql0(i,j,k) / (Nc(i,j,k) + eps)
              Ncp(i,j,k) = Ncp(i,j,k) - ac / xc * rhof(k)
            end if
          endif
        enddo
      enddo
    enddo

    call timer_toc(routine)

  end subroutine accretion_kk

  !> Calculate the evaporation term.
  !!
  !! \param ql0 Liquid water mixing ratio.
  !! \param qt0 Total water mixing ratio.
  !! \param qvsl
  !! \param esl
  !! \param tmp0 Temperature.
  !! \param qrm Rain water mixing ratio at previous time step.
  !! \param Nrm Rain drop number concentration at previous time step.
  !! \param Nr Rain drop number concentration.
  !! \param rhof Density at full levels.
  !! \param exnf Exner function at full levels.
  !! \param qrbase Lowest level with rain.
  !! \param qrroof Highest level with rain.
  !! \param Dvr Rain water mean diameter.
  !! \param xr Mean mass of rain drops.
  !! \param delt Time step size.
  !! \param thlpmcr Tendency of $\theta_l$.
  !! \param qtpmcr Tendency of total water mixing ratio.
  !! \param qrp Tendency of rain water mixing ratio.
  !! \param Nrp Tendency of rain drop number concentration.
  subroutine evaporation_kk(ql0, qt0, qvsl, esl, tmp0, qrm, Nrm, Nr, qr, rhof, exnf, &
                            qrbase, qrroof, delt, thlpmcr, qtpmcr, qrp, Nrp)
    real(field_r), intent(in)    :: ql0(2-ih:i1+ih,2-jh:j1+jh,1:k1)
    real(field_r), intent(in)    :: qt0(2-ih:i1+ih,2-jh:j1+jh,1:k1)
    real(field_r), intent(in)    :: qvsl(2-ih:i1+ih,2-jh:j1+jh,1:k1)
    real(field_r), intent(in)    :: esl(2-ih:i1+ih,2-jh:j1+jh,1:k1)
    real(field_r), intent(in)    :: tmp0(2-ih:i1+ih,2-jh:j1+jh,1:k1)
    real(field_r), intent(in)    :: qrm(2-ih:i1+ih,2-jh:j1+jh,1:k1)
    real(field_r), intent(in)    :: Nrm(2-ih:i1+ih,2-jh:j1+jh,1:k1)
    real(field_r), intent(in)    :: Nr(2:i1,2:j1,1:k1)
    real(field_r), intent(in)    :: qr(2:i1,2:j1,1:k1)

    real(field_r), intent(in)    :: rhof(1:k1)
    real(field_r), intent(in)    :: exnf(1:k1)

    integer,       intent(in)    :: qrbase, qrroof

    real(field_r), intent(in)    :: delt

    real(field_r), intent(inout) :: thlpmcr(2:i1,2:j1,1:k1)
    real(field_r), intent(inout) :: qtpmcr(2-ih:i1+ih,2-jh:j1+jh,1:k1)
    real(field_r), intent(inout) :: qrp(2:i1,2:j1,1:k1)
    real(field_r), intent(inout) :: Nrp(2:i1,2:j1,1:k1)

    character(len=*), parameter :: routine = modname//'/evaporation_kk'

    integer       :: i, j, k
    real(field_r) :: S, G
    real(field_r) :: evap, Nevap
    real(field_r) :: xr, dvr

    call timer_tic(routine, 1)

    if (qrbase > qrroof) then
      call timer_toc(routine)
      return
    end if

    !$acc parallel loop collapse(3) default(present) private(S, G, evap, Nevap)
!!$omp target teams loop private(s,g,evap,nevap) collapse(3)&
!!$omp defaultmap(present:aggregate) defaultmap(present:allocatable)
    do k = qrbase, qrroof
      do j = 2, j1
        do i = 2, i1
          if (qr(i,j,k) > qrmin) then
            xr = calc_xr(rhof(k), qr(i,j,k), nr(i,j,k), xrmin, xrmax)
            dvr = calc_dvr(xr)

            S = min(0.0_field_r, (qt0(i,j,k) - ql0(i,j,k)) / qvsl(i,j,k) - 1)
            G = (Rv * tmp0(i,j,k)) / (Dv * esl(i,j,k)) + rlv / &
                (Kt * tmp0(i,j,k)) * (rlv / (Rv * tmp0(i,j,k)) - 1)
            G = 1 / G

            evap = c_evap * 2 * pi * Dvr * G * S * Nr(i,j,k) / rhof(k)
            Nevap = evap * rhof(k) / xr

            if (evap < - qrm(i,j,k) / delt) then
              Nevap = - Nrm(i,j,k) / delt
              evap  = - qrm(i,j,k) / delt
            endif

            qrp(i,j,k) = qrp(i,j,k) + evap
            Nrp(i,j,k) = Nrp(i,j,k) + Nevap

            qtpmcr(i,j,k) = qtpmcr(i,j,k) - evap
            thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlv / (cp * exnf(k))) * evap
          endif
        enddo
      enddo
    enddo

    call timer_toc(routine)

  end subroutine evaporation_kk

  !> Calculate the sedimentation rate of the rain water content.
  !!
  !! @param[in] qr Rain water content.
  !! @param[in] nr Rain droplet number concentration.
  !! @param[in] rho Air density.
  !!
  !! @returns sedimentation rate of nr.
  elemental function calc_sed_qr_kk(qr, nr, rho) result(sed_qr)
!!$omp declare target

    real(field_r), intent(in) :: qr, nr, rho
    
    real(field_r) :: xr, dvr, wfall_qr, sed_qr

    !$acc routine seq

    xr = calc_xr(rho, qr, nr, xrmin, xrmax)
    dvr = calc_dvr(xr)

    wfall_qr = max(0._field_r, 0.006_field_r * 1E6 * dvr - 0.2_field_r)

    sed_qr  = wfall_qr * qr * rho ! m/s * kg/m3

  end function calc_sed_qr_kk

  !> Calculate the sedimentation rate of the rain water content.
  !!
  !! @param[in] qr Rain water content.
  !! @param[in] nr Rain droplet number concentration.
  !! @param[in] rho Air density.
  !!
  !! @returns sedimentation rate of nr.
  elemental function calc_sed_nr_kk(qr, nr, rho) result(sed_nr)
!!$omp declare target

    real(field_r), intent(in) :: qr, nr, rho
    
    real(field_r) :: xr, dvr, wfall_nr, sed_nr

    !$acc routine seq

    xr = calc_xr(rho, qr, nr, xrmin, xrmax)
    dvr = calc_dvr(xr)

    wfall_nr = max(0._field_r, 0.0035_field_r * 1E6 * dvr - 0.1_field_r)

    sed_nr  = wfall_nr * nr

  end function calc_sed_nr_kk

  !> Calculate the sedimentation term.
  !!
  !! \param qr Rain water mixing ratio.
  !! \param Nr Rain drop number concentration.
  !! \param rhof Density at full levels.
  !! \param dzf Thickness of vertical levels.
  !! \param qrbase Lowest level with rain.
  !! \param qrroof Highest level with rain.
  !! \param delt Time step size.
  !! \param Dvr Rain water mean diameter.
  !! \param xr Mean mass of rain drops.
  !! \param qrp Tendency of rain water mixing ratio.
  !! \param Nrp Tendency of rain drop number concentration.
  !! \param precep Precipitation.
#ifndef DALES_GPU
  subroutine sedimentation_rain_kk(qr, Nr, rhof, dzf, qrbase, qrroof, &
                                delt, qrp, Nrp, precep)
    real(field_r), intent(in)    :: qr(2:i1,2:j1,1:k1)
    real(field_r), intent(in)    :: Nr(2:i1,2:j1,1:k1)
    real(field_r), intent(in)    :: rhof(1:k1)
    real(field_r), intent(in)    :: dzf(1:k1)

    integer,       intent(inout) :: qrbase
    integer,       intent(in)    :: qrroof

    real(field_r), intent(in)    :: delt

    real(field_r), intent(inout) :: qrp(2:i1,2:j1,1:k1)
    real(field_r), intent(inout) :: Nrp(2:i1,2:j1,1:k1)
    real(field_r), intent(out)   :: precep(2:i1,2:j1,1:k1)

    character(len=*), parameter :: routine = modname//'/sedimentation_rain_kk'

    integer       :: i, j, k, jn
    integer       :: n_spl      !<  sedimentation time splitting loop
    real(field_r) :: sed_qr
    real(field_r) :: sed_Nr
    real(field_r) :: xr, dvr

    real(field_r), allocatable :: qr_spl(:,:,:), Nr_spl(:,:,:)

    real(field_r) :: dt_spl

    precep(:,:,:) = 0 ! zero the precipitation flux field
                      ! the update below is not always performed

    call timer_tic(routine, 1)

    if (qrbase > qrroof) then
      call timer_toc(routine)
      return
    end if

    allocate(qr_spl(2:i1,2:j1,1:k1))
    allocate(Nr_spl(2:i1,2:j1,1:k1))

    n_spl = ceiling(wfallmax * delt / minval(dzf))
    dt_spl = delt / real(n_spl, kind=field_r)

    do jn = 1, n_spl ! time splitting loop
      if (jn == 1) then
        qr_spl(:,:,:) = qr(:,:,:)
        Nr_spl(:,:,:) = Nr(:,:,:)
      else
        ! lower the rain base by one level to include the rain fall
        ! from the previous step
        qrbase = max(1, qrbase - 1)
      end if

      do k = qrbase, qrroof
        do j = 2, j1
          do i = 2, i1
            if (qr_spl(i,j,k) > qrmin) then
              xr = calc_xr(rhof(k), qr_spl(i,j,k), nr(i,j,k), xrmin, xrmax)
              dvr = calc_dvr(xr)

              sed_qr = max(0.0_field_r, 0.006_field_r * 1E6_field_r * Dvr - 0.2_field_r) * qr_spl(i,j,k) * rhof(k)
              sed_Nr = max(0.0_field_r, 0.0035_field_r * 1E6_field_r * Dvr - 0.1_field_r) * Nr_spl(i,j,k)

              qr_spl(i,j,k) = qr_spl(i,j,k) - sed_qr * dt_spl / (dzf(k) * rhof(k))
              Nr_spl(i,j,k) = Nr_spl(i,j,k) - sed_Nr * dt_spl / dzf(k)

              if (k > 1) then
                qr_spl(i,j,k-1) = qr_spl(i,j,k-1) + sed_qr * dt_spl / (dzf(k-1) * rhof(k-1))
                Nr_spl(i,j,k-1) = Nr_spl(i,j,k-1) + sed_Nr * dt_spl / dzf(k-1)
              endif
              if (jn==1) then
                precep(i,j,k) = sed_qr / rhof(k)   ! kg kg-1 m s-1
              endif
            endif
          enddo
        enddo
      enddo
    end do ! time splitting loop

    ! the last time splitting step lowered the base level
    ! and we still need to adjust for it
    qrbase = max(1, qrbase - 1)

    Nrp(:,:,qrbase:qrroof) = Nrp(:,:,qrbase:qrroof) + &
      (Nr_spl(:,:,qrbase:qrroof) - Nr(:,:,qrbase:qrroof))/delt

    qrp(:,:,qrbase:qrroof) = qrp(:,:,qrbase:qrroof) + &
      (qr_spl(:,:,qrbase:qrroof) - qr(:,:,qrbase:qrroof))/delt

    deallocate(qr_spl, Nr_spl)

    call timer_toc(routine)

  end subroutine sedimentation_rain_kk
#else
  subroutine sedimentation_rain_kk(qr, Nr, rhof, dzf, qrbase, qrroof, &
                                    delt, qrp, Nrp, precep)
    real(field_r), intent(in)    :: qr(2:i1,2:j1,1:k1)
    real(field_r), intent(in)    :: Nr(2:i1,2:j1,1:k1)
    real(field_r), intent(in)    :: rhof(1:k1)
    real(field_r), intent(in)    :: dzf(1:k1)

    integer,       intent(inout) :: qrbase
    integer,       intent(in)    :: qrroof

    real(field_r), intent(in)    :: delt

    real(field_r), intent(inout) :: qrp(2:i1,2:j1,1:k1)
    real(field_r), intent(inout) :: Nrp(2:i1,2:j1,1:k1)
    real(field_r), intent(out)   :: precep(2:i1,2:j1,1:k1)

    character(len=*), parameter :: routine = modname//'/sedimentation_rain_kk'

    integer       :: i, j, k, jn, sedimbase
    integer       :: n_spl      !<  sedimentation time splitting loop
    real(field_r) :: sed_qr
    real(field_r) :: sed_Nr
    real(field_r) :: dt_spl
    real(field_r) :: delt_inv
    real(field_r) :: xr, dvr

    real(field_r), allocatable :: qr_spl(:,:,:), Nr_spl(:,:,:)
    real(field_r), allocatable :: qr_tmp(:,:,:), Nr_tmp(:,:,:)

    !$acc parallel loop collapse(3) default(present)
!!$omp target teams loop collapse(3) defaultmap(present:aggregate)&
!!$omp defaultmap(present:allocatable)
    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          precep(i,j,k) = 0.0
        end do
      end do
    end do

    call timer_tic(routine, 1)

    if (qrbase > qrroof) then
      call timer_toc(routine)
      return
    end if

    allocate(qr_spl(2:i1,2:j1,1:k1))
    allocate(Nr_spl(2:i1,2:j1,1:k1))
    allocate(qr_tmp(2:i1,2:j1,1:k1))
    allocate(Nr_tmp(2:i1,2:j1,1:k1))

    !$acc enter data create(qr_spl, Nr_spl, qr_tmp, Nr_tmp)
!$omp target enter data map(alloc:qr_spl,nr_spl,qr_tmp,nr_tmp)

    n_spl = ceiling(wfallmax * delt / minval(dzf))
    dt_spl = delt / real(n_spl, kind=field_r)

    do jn = 1, n_spl ! time splitting loop
      if (jn == 1) then
        !$acc parallel loop collapse(3) default(present)
!!$omp target teams loop collapse(3) defaultmap(present:aggregate)&
!!$omp defaultmap(present:allocatable)
        do k = 1, k1
          do j = 2, j1
            do i = 2, i1
              qr_spl(i,j,k) = qr(i,j,k)
              Nr_spl(i,j,k) = Nr(i,j,k)
              qr_tmp(i,j,k) = qr(i,j,k)
              Nr_tmp(i,j,k) = Nr(i,j,k)
            end do
          end do
        end do
      else
        !Copy from tmp into spl
        !$acc parallel loop collapse(3) default(present)
!!$omp target teams loop collapse(3) defaultmap(present:aggregate)&
!!$omp defaultmap(present:allocatable)
        do k = 1, k1
          do j = 2, j1
            do i = 2, i1
              qr_spl(i,j,k) = qr_tmp(i,j,k)
              Nr_spl(i,j,k) = Nr_tmp(i,j,k)
            end do
          end do
        end do

        ! lower the rain base by one level to include the rain fall
        ! from the previous step
        qrbase = max(1, qrbase - 1)
      end if

      ! Compute precep
      if (jn == 1) then
        !$acc parallel loop collapse(3) default(present)
!!$omp target teams loop collapse(3) defaultmap(present:aggregate)&
!!$omp defaultmap(present:allocatable)
        do k = qrbase, qrroof
          do j = 2, j1
            do i = 2, i1
              if (qr_spl(i,j,k) > qrmin) then
                xr = calc_xr(rhof(k), qr_spl(i,j,k), nr(i,j,k), xrmin, xrmax)
                dvr = calc_dvr(xr)
                precep(i,j,k) = max(0.0_field_r, 0.006_field_r * 1E6_field_r * dvr - 0.2_field_r) * qr_spl(i,j,k)
              endif
            enddo
          enddo
        enddo
      end if ! jn == 1

      sedimbase = qrbase

      ! k qrbase if == 1
      if (qrbase == 1) then
        sedimbase = sedimbase + 1
        k = 1

        !$acc parallel loop collapse(2) default(present) private(sed_qr, sed_Nr)
!!$omp target teams loop private(sed_qr,sed_nr) collapse(2)&
!!$omp defaultmap(present:aggregate) defaultmap(present:allocatable)
        do j = 2, j1
          do i = 2, i1
            if (qr_spl(i,j,k) > qrmin) then
              xr = calc_xr(rhof(k), qr_spl(i,j,k), nr(i,j,k), xrmin, xrmax)
              dvr = calc_dvr(xr)

              sed_qr = max(0.0_field_r, 0.006_field_r *1E6_field_r * dvr - 0.2_field_r) * qr_spl(i,j,k) * rhof(k)
              sed_Nr = max(0.0_field_r, 0.0035_field_r *1E6_field_r * dvr - 0.1_field_r) * Nr_spl(i,j,k)

              qr_tmp(i,j,k) = qr_tmp(i,j,k) - sed_qr * dt_spl / (dzf(k) * rhof(k))
              Nr_tmp(i,j,k) = Nr_tmp(i,j,k) - sed_Nr * dt_spl / dzf(k)
            endif
          enddo
        enddo
      end if ! qrbase == 1

      !$acc parallel loop collapse(3) default(present) private(sed_qr, sed_Nr)
!!$omp target teams loop private(sed_qr,sed_nr) collapse(3)&
!!$omp defaultmap(present:aggregate) defaultmap(present:allocatable)
      do k = sedimbase, qrroof
        do j = 2, j1
          do i = 2, i1
            if (qr_spl(i,j,k) > qrmin) then
              xr = calc_xr(rhof(k), qr_spl(i,j,k), nr(i,j,k), xrmin, xrmax)
              dvr = calc_dvr(xr)

              sed_qr = max(0.0_field_r, 0.006_field_r *1E6_field_r * dvr - 0.2_field_r) * qr_spl(i,j,k) * rhof(k)
              sed_Nr = max(0.0_field_r, 0.0035_field_r *1E6_field_r * dvr - 0.1_field_r) * Nr_spl(i,j,k)

              !$acc atomic update
!!$omp atomic update
              qr_tmp(i,j,k) = qr_tmp(i,j,k) - sed_qr*dt_spl/(dzf(k)*rhof(k))
              !$acc atomic update
!!$omp atomic update
              Nr_tmp(i,j,k) = Nr_tmp(i,j,k) - sed_Nr*dt_spl/dzf(k)

              !$acc atomic update
!!$omp atomic update
              qr_tmp(i,j,k-1) = qr_tmp(i,j,k-1) + sed_qr*dt_spl/(dzf(k-1)*rhof(k-1))
              !$acc atomic update
!!$omp atomic update
              Nr_tmp(i,j,k-1) = Nr_tmp(i,j,k-1) + sed_Nr*dt_spl/dzf(k-1)
            endif
          enddo
        enddo
      enddo
    end do ! time splitting loop

    ! the last time splitting step lowered the base level
    ! and we still need to adjust for it
    qrbase = max(1, qrbase - 1)

    delt_inv = 1 / delt

    !$acc parallel loop collapse(3) default(present)
!!$omp target teams loop collapse(3) defaultmap(present:aggregate)&
!!$omp defaultmap(present:allocatable)
    do k = qrbase, qrroof
      do j = 2, j1
        do i = 2, i1
          Nrp(i,j,k) = Nrp(i,j,k) + (Nr_tmp(i,j,k) - Nr(i,j,k)) * delt_inv
          qrp(i,j,k) = qrp(i,j,k) + (qr_tmp(i,j,k) - qr(i,j,k)) * delt_inv
        end do
      end do
    end do

    !$acc exit data delete(qr_spl, Nr_spl, qr_tmp, Nr_tmp)
!$omp target exit data map(delete:qr_spl,nr_spl,qr_tmp,nr_tmp)

    deallocate(qr_spl, Nr_spl, qr_tmp, Nr_tmp)

    call timer_toc(routine)

  end subroutine sedimentation_rain_kk
#endif

end module bulkmicro_kk
