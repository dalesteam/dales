module modaerosol

  use, intrinsic :: iso_fortran_env

  use modaerosol_mode_t, only: aerosol_mode_t, hydrometeor_mode_t, mode_t, &
                               mode_container_t, connect_modes
  use modaerosol_common, only: maxspecies, maxmodes, iNUS, iAIS, iACS, iCOS, &
                               iAII, iACI, iCOI, iINC, iINR, iSO4, iSS, iPOM, &
                               iBC, iDU
  use modglobal,         only: ifnamopt, fname_options, checknamelisterror, &
                               cexpnr, i1, j1, k1, ih, jh, pi, nsv, rhow, kmax, &
                               rk3step, rd, pirhow
  use modfields,         only: sv0, svp
  use modmicrodata,      only: qcmin, delt
  use modmpi,            only: myid, D_MPI_BCAST, commwrld, mpierr
  use modprecision,      only: field_r
  use modtimer,          only: timer_tic, timer_toc
  use modbulkmicro_data, only: l_sb, qrmin, l_mur_cst, mur_cst
  use bulkmicro_sb,      only: calc_sed_qr_sb, calc_sed_nr_sb
  use bulkmicro_kk,      only: calc_sed_nr_kk, calc_sed_qr_kk
  use modstat_nc

  implicit none

  private

  character(len=*), parameter :: modname = 'modaerosol'

  public :: laerosol
  public :: aerosol_read_namelist
  public :: init_aerosol
  public :: aerosol_prepare
  public :: aerosol_finish
  public :: aerosol_activation
  public :: aerosol_cloud_to_rain
  public :: aerosol_resuspend_rain
  public :: aerosol_resuspend_cloud
  public :: aerosol_sedimentation_rain
  public :: aerosol_scavenging_rain

  interface erfcinv
    module procedure :: erfcinv_real32
    module procedure :: erfcinv_real64
  end interface

  ! Public variables
  logical :: &
    laerosol = .false., & !< Switch for enabling aerosol scheme
    lso4 = .false.,     & !< Switch for enabling sulphur
    lss = .false.,      & !< Switch for enabling sea salt
    lpom = .false.,     & !< Switch for enabling particulate organic matter
    lbc = .false.,      & !< Switch for enabling black carbon
    ldu = .false.         !< Switch for enabling dust

  real(field_r), allocatable :: &
    sed_qr(:,:,:), & ! Rain sedimentation rate, needed for scavenging.
    qlm(:,:,:)       ! Cloud water at previous time step, needed for resuspension.

  type(mode_container_t) :: &
    modes(9)           ! List of all modes.

  type(aerosol_mode_t), target :: &
    modes_f(7)         ! List of free aerosol modes.

  type(hydrometeor_mode_t), target :: &
    modes_h(iINC:iINR) ! List of in-hydrometeor modes. (currently only cloud and rain)

contains

  include 'erfcinv.inc'

  subroutine aerosol_read_namelist(nml_filename)

    character(len=*), intent(in) :: &
      nml_filename

    integer :: &
      ierr

    namelist /NAMAEROSOL/ laerosol, lso4, lss, lpom, lbc, ldu

    if (myid == 0) then
      ! Namelist
      open(ifnamopt, file=nml_filename, status='old', iostat=ierr)
      read(ifnamopt, NAMAEROSOL, iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMAEROSOL')
      close(ifnamopt)
    end if

    call d_mpi_bcast(laerosol, 1, 0, commwrld, mpierr)
    call d_mpi_bcast(lso4, 1, 0, commwrld, mpierr)
    call d_mpi_bcast(lss, 1, 0, commwrld, mpierr)
    call d_mpi_bcast(lpom, 1, 0, commwrld, mpierr)
    call d_mpi_bcast(lbc, 1, 0, commwrld, mpierr)
    call d_mpi_bcast(ldu, 1, 0, commwrld, mpierr)

  end subroutine aerosol_read_namelist

  !> Initialize aerosol module.
  subroutine init_aerosol()

    character(len=*), parameter :: &
      routine = modname//"/initaerosol"

    integer :: &
      imod

    logical :: &
      mode_config(maxspecies,maxmodes) = .false.

    if (.not. laerosol) return

    do imod = 1, 7
      modes(imod)%p => modes_f(imod)
    end do

    do imod = iINC, iINR
      modes(imod)%p => modes_h(imod)
    end do

    ! Configure the modes.

    !                | NUS | AIS | ACS | COS | AII | ACI | COI | INC | INR |
    ! SO4            |  x  |  x  |  x  |  x  |     |     |     |  x  |  x  |
    ! Sea Salt       |     |     |  x  |  x  |     |     |     |  x  |  x  |
    ! Organic Matter |     |  x  |  x  |  x  |  x  |     |     |  x  |  x  |
    ! Black Carbon   |     |  x  |  x  |  x  |  x  |     |     |  x  |  x  |
    ! Dust           |     |     |  x  |  x  |     |  x  |  x  |  x  |  x  |

    if (lso4) mode_config(iSO4,:) = [.true., .true., .true., .true., .false., &
                                     .false., .false., .true., .true.]
    if (lss)  mode_config(iSS,:)  = [.false., .false., .true., .true., &
                                     .false., .false., .false., .true., .true.]
    if (lpom) mode_config(iPOM,:) = [.false., .true., .true., .true., .true., &
                                     .false., .false., .true., .true.]
    if (lbc)  mode_config(iBC,:)  = [.false., .true., .true., .true., .true., &
                                     .false., .false., .true., .true.]
    if (ldu)  mode_config(iDU,:)  = [.false., .false., .true., .true., &
                                     .false., .true., .true., .true., .true.]

    do imod = 1, maxmodes
      call modes(imod)%p%init(imod, mode_config(:,imod))
    end do

    ! Connect free aerosol modes to in-hydrometeor modes
    do imod = 1, maxmodes
      select type(mode => modes(imod)%p)
        class is (aerosol_mode_t)
          call connect_modes(mode, modes(iINC)%p, mode%to_hydro)
          print *, mode%name, " ", mode%to_hydro%cnct
      end select
    end do

    allocate(sed_qr(2:i1,2:j1,k1), qlm(2:i1,2:j1,k1))

  end subroutine init_aerosol

  !> Prepares aerosol fields for microphysics calculations.
  subroutine aerosol_prepare()

    character(len=*), parameter :: &
      routine = modname//'/aerosol_prepare'

    integer :: &
      imod

    call timer_tic(routine, 1)

    do imod = 1, maxmodes
      call modes(imod)%p%prepare(sv0)
    end do

    call timer_toc(routine)

  end subroutine aerosol_prepare

  !> Copy out tendencies.
  subroutine aerosol_finish()

    character(len=*), parameter :: &
      routine = modname//'/aerosol_finish'

    integer :: &
      imod

    call timer_tic(routine, 1)

    do imod = 1, maxmodes
      call modes(imod)%p%finish(svp)
    end do

    call timer_toc(routine)

  end subroutine aerosol_finish

  pure function calc_mean_rho(q, rho) result(rho_m)
    
    real(field_r), intent(in) :: q(:), rho(:)

    real(field_r) :: &
      m,             &
      rho_m

    integer :: &
      s

    m = 0
    rho_m = 0

    do s = 1, size(q)
      m = m + q(s)
      rho_m = rho_m + q(s) / rho(s)
    end do

    m = max(0.0_field_r, m)
    rho_m = max(0.0_field_r, m / (rho_m + 1E-16))

  end function calc_mean_rho

  !> Compute the median diameter of a log-normal distribution.
  !!
  !! @param[in] n Number concentration.
  !! @param[in] q Mass concentration (dim=nspecies).
  !! @param[in] rho Aerosol densities (dim=nspecies).
  !! @param[in] sig_g Geometric standard deviation of the distribution.
  pure function calc_median_diameter(n, q, rho, sig_g) result(dm)

    real(field_r), intent(in) :: n, q(:), rho(:), sig_g

    !$acc routine seq

    real(field_r) :: &
      m,     & ! Total aerosol mass.
      rho_m, & ! Mean density.
      dm       ! Median diameter.

    integer :: &
      s ! Loop index

    m = 0
    rho_m = 0

    do s = 1, size(q)
      m = m + q(s)
      rho_m = rho_m + q(s) / rho(s)
    end do

    m = max(0.0_field_r, m)
    rho_m = max(0.0_field_r, m / (rho_m + 1E-16))

    dm = ((6 * m) / (pi * n * rho_m + 1E-16))**(1.0_field_r / 3) &
         * exp(- 0.5_field_r * 3 * log(sig_g) * log(sig_g))

    dm = max(0.0_field_r, dm)

  end function calc_median_diameter

  !> Aerosol activation based on updraft velocity.
  !!
  !! @see https://doi.org/10.5194/acp-15-9217-2015
  !!
  !! @param[in] ql Cloud water content.
  !! @param[in] w Vertical velocity.
  !! @param[in] nc Cloud droplet number concentration.
  !! @param[in] delt Time step size.
  !! @param[inout] ncp Tendency of cloud droplet number concentration.
  subroutine aerosol_activation(ql, w, nc, delt, ncp)

    real(field_r), intent(in) :: &
      ql(2-ih:,2-jh:,:),         &
      w(2-ih:,2-jh:,:),          &
      nc(2:,2:,:),               &
      delt

    real(field_r), intent(inout) :: &
      ncp(2:,2:,:)

    character(*),  parameter :: &
      routine = modname//"/aerosol_activation"

    real(field_r), parameter :: &
      r_crit = 35E-9 ! Critical radius for activation.

    class(aerosol_mode_t), pointer :: &
      m_ais, & ! Soluble Aitken mode.
      m_acs, & ! Soluble accumulation mode.
      m_cos    ! Soluble coarse mode.

    class(hydrometeor_mode_t), pointer :: &
      m_inc ! In-cloud mode.

    integer :: &
      i, j, k, s, m, & ! Loop indices.
      st               ! Target index in in-cloud mode.

    real(field_r) :: &
      dm,            & ! Median diameter.
      n_act,         & ! Number of activated particles.
      dncdt,         & ! Potential activation tendency.
      fn,            & ! Activated fraction of number concentration.
      fm,            & ! Activated fraciton of mass concentrattion.
      tend_n,        & ! Real tendency of number concentration.
      tend_m,        & ! Real tendency of mass concentration.
      w0               ! Updraft velocity.

    call timer_tic(routine, 1)

    m_ais => modes_f(iAIS)
    m_acs => modes_f(iACS)
    m_cos => modes_f(iCOS)
    m_inc => modes_h(iINC)

    do k = 1, kmax
      do j = 2, j1
        do i = 2, i1
          if (ql(i,j,k) > qcmin) then
            if (m_ais%nspecies > 0) then
              dm = calc_median_diameter(m_ais%n(i,j,k), m_ais%q(i,j,k,:), &
                                        m_ais%rho, m_ais%sig_g)
              if (dm > 0) then
              fn = 1 - 0.5_field_r * erfc(-log(2 * r_crit / &
                   dm + 1E-30) / (sqrt(2.0_field_r) * log(m_ais%sig_g)))
              end if
            end if

            n_act = 1E-6 * (m_acs%n(i,j,k) + m_cos%n(i,j,k) + &
                            fn * m_ais%n(i,j,k))

            w0 = max(0.0_field_r, w(i,j,k))
            dncdt = 1E6 / delt * &
                    (0.1 * (w0 * 100 * n_act / &
                     (w0 * 100 + 0.023_field_r * n_act + 1E-20)))**1.27_field_r &
                    - 1E-6 * Nc(i,j,k)
            dncdt = max(dncdt, 0.0_field_r)

            fn = dncdt * delt / (m_cos%n(i,j,k) + 1E-20)
            fn = max(min(fn, 1.0_field_r), 0.0_field_r)
            fm = 1 - 0.5_field_r * erfc(erfcinv(2 * fn) &
                 - 3 * log(m_cos%sig_g) / sqrt(2.0_field_r))

            tend_n = fn * m_cos%n(i,j,k) / delt
            tend_n = max(0.0_field_r, tend_n)
            m_cos%np(i,j,k) = m_cos%np(i,j,k) - tend_n
            Ncp(i,j,k) = Ncp(i,j,k) + tend_n

            do s = 1, m_cos%nspecies
              tend_m = max(0.0_field_r, fm * m_cos%q(i,j,k,s) / delt)
              st = m_cos%to_hydro%cnct(2,s)
              m_cos%qp(i,j,k,s) = m_cos%qp(i,j,k,s) - tend_m
              m_inc%qp(i,j,k,st) = m_inc%qp(i,j,k,st) + tend_m
            end do

            dncdt = dncdt - tend_n
            dncdt = max(0.0_field_r, dncdt)

            if (dncdt > 0) then
              fn = dncdt * delt / (m_acs%n(i,j,k) + 1E-20)
              fn = max(min(fn, 1.0_field_r), 0.0_field_r)
              fm = 1 - 0.5_field_r * erfc(erfcinv(2 * fn) &
                - 3 * log(m_acs%sig_g) / sqrt(2.0_field_r))

              tend_n = fn * m_acs%n(i,j,k) / delt
              tend_n = max(0.0_field_r, tend_n)

              m_acs%np(i,j,k) = m_acs%np(i,j,k) - tend_n
              Ncp(i,j,k) = Ncp(i,j,k) + tend_n

              do s = 1, m_acs%nspecies
                tend_m = max(0.0_field_r, fm * m_acs%q(i,j,k,s) / delt)
                st = m_acs%to_hydro%cnct(2,s)
                m_acs%qp(i,j,k,s) = m_acs%qp(i,j,k,s) - tend_m
                m_inc%qp(i,j,k,st) = m_inc%qp(i,j,k,st) + tend_m
              end do

              dncdt = dncdt - tend_n
              dncdt = max(0.0_field_r, dncdt)
            end if

            if (m_ais%nspecies > 0 .and. dncdt > 0) then
              fn = dncdt * delt / (m_ais%n(i,j,k) + 1E-20)
              fn = max(min(fn, 1.0_field_r), 0.0_field_r)
              fm = 1 - 0.5_field_r * &
                erfc(erfcinv(2 * fn) - 3 * log(m_ais%sig_g) / sqrt(2.0_field_r))
              fm = merge(1.0_field_r, fm, fn > 1.0_field_r)

              tend_n = fn * m_ais%n(i,j,k) / delt
              tend_n = max(0.0_field_r, tend_n)

              m_ais%np(i,j,k) = m_ais%np(i,j,k) - tend_n
              Ncp(i,j,k) = Ncp(i,j,k) + tend_n

              do s = 1, m_ais%nspecies
                tend_m = max(0.0_field_r, fm * m_ais%q(i,j,k,s) / delt)
                st = m_ais%to_hydro%cnct(2,s)
                m_ais%qp(i,j,k,s) = m_ais%qp(i,j,k,s) - tend_m
                m_inc%qp(i,j,k,st) = m_inc%qp(i,j,k,st) + tend_m
              end do
            end if
          end if
        end do
      end do
    end do
 
    call timer_toc(routine)

  end subroutine aerosol_activation

  !> Move aerosol from in-cloud to in-rain mode based on the tendency of some
  !! rain generation process.
  !!
  !! @param[in] qc Cloud water content.
  !! @param[in] qrp Tendency of rain water content.
  subroutine aerosol_cloud_to_rain(qc, qrp)

    real(field_r), intent(in) :: &
      qc(2-ih:,2-jh:,:),         &
      qrp(2:,2:,:)

    character(len=*), parameter :: &
      routine = modname//'/aerosol_cloud_to_rain'

    class(hydrometeor_mode_t), pointer :: &
      m_inc, & ! In-cloud mode.
      m_inr    ! In-rain mode.

    integer :: &
      i, j, k, s ! Loop indices

    real(field_r) :: &
      dqadt ! Tendency of in-rain aerosol

    call timer_tic(routine, 1)

    m_inc => modes_h(iINC)
    m_inr => modes_h(iINR)

    !$acc parallel loop collapse(4) default(present) private(dqadt)
    do s = 1, m_inc%nspecies
      do k = 1, kmax
        do j = 2, j1
          do i = 2, i1
            if (qrp(i,j,k) > 0) then
              dqadt = qrp(i,j,k) / qc(i,j,k) * m_inc%q(i,j,k,s)
              m_inc%qp(i,j,k,s) = m_inc%qp(i,j,k,s) - dqadt
              m_inr%qp(i,j,k,s) = m_inr%qp(i,j,k,s) + dqadt
            end if
          end do
        end do
      end do
    end do

    call timer_toc(routine)

  end subroutine aerosol_cloud_to_rain

  !> Computes the resuspension of aerosol particles by evaporating rain droplets.
  !!
  !! Aerosols are resuspended over the ACS and COS modes. Only rain drops that 
  !! fully evaporate should resuspend an aerosol particle, so the resuspended
  !! aerosol mass is corrected using a correction factor.
  !!
  !! @see https://doi.org/10.1016/j.atmosres.2005.10.012
  !!
  !! @param[in] qr Rain water content.
  !! @param[in] qrp Tendency of rain water content from evaporation.
  !! @param[in] nrp Tendency of rain number concentration from evaporation.
  !! @param[in] delt Time step size.
  subroutine aerosol_resuspend_rain(qr, qrp, nrp, delt)

    real(field_r), intent(in) :: &
      qr(2:,2:,:),               &
      qrp(2:,2:,:),              &
      nrp(2:,2:,:),              &
      delt

    character(len=*), parameter :: routine = &
      modname//'/aerosol_resuspend_rain'

    real(field_r), parameter :: &
      Dc = 1E-9 ! Diameter separating the accumulation and coarse modes.

    class(aerosol_mode_t), pointer :: &
      m_acs, & ! Soluble accumulation mode.
      m_cos    ! Soluble coarse mode.

    class(hydrometeor_mode_t), pointer :: &
      m_inr ! In-rain mode

    integer :: &
      i, j, k, s, & ! Loop indices.
      itype,      & ! Aerosol type.
      target_idx    ! Index in target mode.

    real(field_r) :: &
      f_evp,             & ! Fraction of evaporated rain water.
      eps,               & ! Correction factor.
      evapm(maxspecies), & ! Mass of evaporated aerosol.
      evapn,             & ! Number of evaporated aerosol.
      dn,                & ! Number median diameter of evaporated aerosol.
      dm,                & ! Mass median diameter of evaporated aerosol.
      fn,                & ! Number fraction of aerosol resuspended in ACS mode.
      fm                   ! Mass fraction of aerosol resuspended in COS mode.

    call timer_tic(routine, 1)

    m_acs => modes_f(iACS)
    m_cos => modes_f(iCOS)
    m_inr => modes_h(iINR)

    do k = 1, kmax
      do j = 2, j1
        do i = 2, i1
          if (qr(i,j,k) > 0) then
            f_evp = (-qrp(i,j,k) * delt) / (qr(i,j,k) + 1E-40)
            f_evp = max(min(f_evp, 1.0_field_r), 0.0_field_r)

            ! Correction factor from Gong et al. (2006).
            ! Evaluates to 1 for f_evp = 1
            eps = (1 - exp(-2 * sqrt(f_evp)) * (1 + 2 * sqrt(f_evp) &
                  + 2 * f_evp + (4.0_field_r/3) * f_evp**(3.0_field_r/2))) &
                  * (1 - f_evp) + f_evp * f_evp
          
            evapm(:) = eps * f_evp * m_inr%q(i,j,k,:) / delt
            evapn = max(0.0_field_r, -1 * nrp(i,j,k))

            ! Compute the median diameter of the resuspended aerosol.
            dn = calc_median_diameter(evapn, evapm(:), m_inr%rho, 1.5_field_r)
            dm = dn * exp(3 * log(1.5_field_r)**2)

            fn = 0.5_field_r * erfc(-log(dc/(dn + 1E-40)) &
                                    / (log(1.5_field_r) * sqrt(2.0_field_r)))
            fm = 0.5_field_r * erfc(-log(dc/(dm + 1E-40)) &
                                    / (log(1.5_field_r) * sqrt(2.0_field_r)))

            m_acs%np(i,j,k) = m_acs%np(i,j,k) + fn * evapn
            m_cos%np(i,j,k) = m_cos%np(i,j,k) + (1 - fn) * evapn

            do s = 1, m_inr%nspecies
              m_inr%qp(i,j,k,s) = m_inr%qp(i,j,k,s) - evapm(s)
              m_acs%qp(i,j,k,s) = m_acs%qp(i,j,k,s) + fm * evapm(s)
              m_cos%qp(i,j,k,s) = m_cos%qp(i,j,k,s) + (1 - fm) * evapm(s)
            end do
          end if
        end do
      end do
    end do

    call timer_toc(routine)

  end subroutine aerosol_resuspend_rain

  !> Computes the resuspension of aerosol particles by evaporating cloud droplets.
  !!
  !! Uses a similar approach as aerosol_resuspend_rain(). Additionally computes
  !! tendency of cloud droplet number concentration due to evaporation.
  !!
  !! @param[in] ql Cloud water at current time step.
  !! @param[in] nc Cloud droplet number concentration.
  !! @param[in] delt Time step size.
  !! @param[inout] ncp Tendency of cloud droplet number concentration.
  subroutine aerosol_resuspend_cloud(ql, nc, delt, ncp)

    real(field_r), intent(in) :: &
      ql(2-ih:,2-jh:,:),         &
      nc(2:,2:,:),               &
      delt

    real(field_r), intent(inout) :: &
      ncp(2:,2:,:)

    character(len=*), parameter :: &
      routine = modname//'/aero_resuspend_cloud'

    real(field_r), parameter :: &
      Dc = 1E-9 ! Diameter separating the accumulation and coarse modes.

    class(aerosol_mode_t), pointer :: &
      m_acs, & ! Soluble accumulation mode.
      m_cos    ! Soluble coarse mode.

    class(hydrometeor_mode_t), pointer :: &
      m_inc ! In-cloud mode

    integer :: &
      i, j, k, s, & ! Loop indices.
      itype,      & ! Aerosol type.
      target_idx    ! Index in target mode.

    real(field_r) :: &
      f_evp,             & ! Fraction of evaporated rain water.
      eps,               & ! Correction factor.
      evapm(maxspecies), & ! Mass of evaporated aerosol.
      evapn,             & ! Number of evaporated aerosol.
      dn,                & ! Number median diameter of evaporated aerosol.
      dm,                & ! Mass median diameter of evaporated aerosol.
      fn,                & ! Number fraction of aerosol resuspended in ACS mode.
      fm                   ! Mass fraction of aerosol resuspended in COS mode.

    call timer_tic(routine, 1)

    m_acs => modes_f(iACS)
    m_cos => modes_f(iCOS)
    m_inc => modes_h(iINC)

    do k = 1, kmax
      do j = 2, j1
        do i = 2, i1
          f_evp = (qlm(i,j,k) - ql(i,j,k)) / (qlm(i,j,k) + 1E-40)
          f_evp = max(min(f_evp, 1.0_field_r), 0.0_field_r)

          ! If no more cloud is in this grid cell, evaporate all aerosol.
          f_evp = merge(f_evp, 1.0_field_r, ql(i,j,k) > 0.0_field_r)

          ! Correction factor from Gong et al. (2006).
          ! Evaluates to 1 for f_evp = 1
          eps = (1 - exp(-2 * sqrt(f_evp)) * (1 + 2 * sqrt(f_evp) &
                + 2 * f_evp + (4.0_field_r/3) * f_evp**(3.0_field_r/2))) &
                * (1 - f_evp) + f_evp * f_evp
          
          evapm(:) = eps * f_evp * m_inc%q(i,j,k,:) / delt
          evapn = f_evp * nc(i,j,k) / delt

          ! Compute the median diameter of the resuspended aerosol.
          dn = calc_median_diameter(evapn, evapm(:), m_inc%rho, 1.5_field_r)
          dm = dn * exp(3 * log(1.5_field_r)**2)

          fn = 0.5_field_r * erfc(-log(dc/(dn + 1E-40)) &
                                  / (log(1.5_field_r) * sqrt(2.0_field_r)))
          fm = 0.5_field_r * erfc(-log(dc/(dm + 1E-40)) &
                                  / (log(1.5_field_r) * sqrt(2.0_field_r)))

          ncp(i,j,k) = ncp(i,j,k) - evapn

          m_acs%np(i,j,k) = m_acs%np(i,j,k) + fn * evapn
          m_cos%np(i,j,k) = m_cos%np(i,j,k) + (1 - fn) * evapn

          do s = 1, m_inc%nspecies
            m_inc%qp(i,j,k,s) = m_inc%qp(i,j,k,s) - evapm(s)
            m_acs%qp(i,j,k,s) = m_acs%qp(i,j,k,s) + fm * evapm(s)
            m_cos%qp(i,j,k,s) = m_cos%qp(i,j,k,s) + (1 - fm) * evapm(s)
          end do
        end do
      end do
    end do

    if (rk3step == 3) then
      do k = 1, k1
        do j = 2, j1
          do i = 2, i1
            qlm(i,j,k) = ql(i,j,k)
          end do
        end do
      end do
    end if

    call timer_toc(routine)

  end subroutine aerosol_resuspend_cloud

  !> Compute flux of in-rain aerosols due to sedimentation of rain drops.
  !!
  !! @param[in] qr Rain water content.
  !! @param[in] nr Rain number concentration.
  !! @param[in] rho Air density.
  !! @param[in] dzf Thickness of vertical levels.
  !! @param[in] qrbase Lowest level with rain.
  !! @param[in] qrroof Highest level with rain.
  !! @param[in] delt Time step size.
  subroutine aerosol_sedimentation_rain(qr, nr, rho, dzf, qrbase, qrroof, delt)

    real(field_r), intent(in)  :: &
      qr(2:,2:,:),                &
      nr(2:,2:,:),                &
      rho(:),                     &
      dzf(:),                     &
      delt

    integer, intent(in)  :: &
      qrbase,               &
      qrroof

    character(len=*), parameter :: &
      routine = modname//'/aerosol_sedimentation_rain'

    class(hydrometeor_mode_t), pointer :: &
      m_inr

    integer ::    &
      i, j, k, s, & ! Loop indices.
      ts,         & ! Time index.
      n_spl         ! Number of sub-timesteps.

    real(field_r) :: &
      dt_spl,        & ! Sub-timestep size.
      sed_nr           ! Sedimentation rate of number concentration.

    real(field_r), pointer :: &
      qr_spl(:,:,:),          & ! Rain water content at sub-timesteps.
      nr_spl(:,:,:),          & ! Rain number concentration at sub-timesteps.
      qa_spl(:,:,:,:)           ! Aerosol mass at sub-timesteps.

    call timer_tic(routine, 1)

    m_inr => modes_h(iINR)

    allocate(qr_spl(2:i1,2:j1,1:k1), nr_spl(2:i1,2:j1,1:k1), &
             qa_spl(1:m_inr%nspecies,2:i1,2:j1,1:k1))

    n_spl = ceiling(9.9 * delt / minval(dzf))
    dt_spl = delt / real(n_spl, kind=field_r)

    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          qr_spl(i,j,k) = qr(i,j,k)
          nr_spl(i,j,k) = nr(i,j,k)
        end do
      end do
    end do

    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          do s = 1, m_inr%nspecies
            qa_spl(s,i,j,k) = m_inr%q(i,j,k,s)
          end do
        end do
      end do
    end do

    do ts = 1, n_spl
      do k = qrbase, qrroof
        do j = 2, j1 
          do i = 2, i1
            if (qr_spl(i,j,k) > qrmin .and. nr_spl(i,j,k) > 0) then
              if (l_sb) then
                sed_qr(i,j,k) = calc_sed_qr_sb(qr_spl(i,j,k), nr_spl(i,j,k), &
                                               rho(k))
                sed_nr = calc_sed_nr_sb(qr_spl(i,j,k), nr_spl(i,j,k), rho(k))
              else
                sed_qr(i,j,k) = calc_sed_qr_kk(qr_spl(i,j,k), nr_spl(i,j,k), &
                                               rho(k))
                sed_nr = calc_sed_nr_kk(qr_spl(i,j,k), nr_spl(i,j,k), rho(k))
              end if
              qr_spl(i,j,k) = qr_spl(i,j,k) - sed_qr(i,j,k) * dt_spl &
                              / (dzf(k) * rho(k))
              nr_spl(i,j,k) = nr_spl(i,j,k) - sed_nr * dt_spl / dzf(k)
              if (k > 1) then
                qr_spl(i,j,k-1) = qr_spl(i,j,k-1) + sed_qr(i,j,k) * dt_spl &
                                  / (dzf(k-1) * rho(k-1))
                nr_spl(i,j,k-1) = nr_spl(i,j,k-1) + sed_nr * dt_spl / dzf(k-1)
              end if
              do s = 1, m_inr%nspecies
                qa_spl(s,i,j,k) = qa_spl(s,i,j,k) - sed_qr(i,j,k) / qr_spl(i,j,k) &
                                  * qa_spl(s,i,j,k) * dt_spl / (dzf(k) * rho(k))
                if (k > 1) then
                  qa_spl(s,i,j,k-1) = qa_spl(s,i,j,k-1) + sed_qr(i,j,k) &
                                      / qr_spl(i,j,k) * qa_spl(s,i,j,k) &
                                      * dt_spl / (dzf(k-1) * rho(k-1))
                end if
              end do
            end if
          end do
        end do
      end do
    end do

    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          do s = 1, m_inr%nspecies
            m_inr%qp(i,j,k,s) = m_inr%qp(i,j,k,s) + &
                               (qa_spl(s,i,j,k) - m_inr%q(i,j,k,s)) / delt
          end do
        end do
      end do
    end do

    deallocate(qr_spl, nr_spl, qa_spl)
    
    call timer_toc(routine)

  end subroutine aerosol_sedimentation_rain

  include 'microphysics.inc'

  subroutine aerosol_scavenging_rain(nr, qr, thl, rho, exn, p, delt)

    real(field_r), intent(in) :: &
      nr(2:,2:,:),               &
      qr(2:,2:,:),               &
      thl(2-ih:,2-jh:,:),        &
      rho(:),                    &
      exn(:),                    &
      p(:),                      &
      delt

    real(field_r), parameter :: &
      mu_a = 1.409E-5,          & ! Viscosity of air.
      mu_w = 1.8E-3,            & ! Viscosity of water.
      kb = 1.38065E-23            ! Boltzmann constant.

    integer :: &
      i, j, k, s, imod, ts

    class(aerosol_mode_t), pointer :: &
      mode

    class(hydrometeor_mode_t), pointer :: &
      m_inr

    real(field_r) :: &
      re,            & ! Reynolds number.
      xr,            & ! Mean mass of rain drops.
      dvr,           & ! Mean rain drop diameter.
      dm,            & ! Median diameter of aerosol particle.
      dn,            & ! Decay rate.
      cnacc,         &
      v,             & ! Terminal velocity of rain drops.
      T,             & ! Temperature.
      E,             & ! Collision efficiency.
      lbd,           & ! Mean free path of air.
      gamma,         & ! Scavenging coefficient.
      csc,           & ! Cunningham slip correction factor.
      sc,            & ! Aerosol Schmidt number.
      db,            & ! Brownian diffusion coefficient.
      tau,           & ! Particle relaxation factor.
      st,            & ! Stokes number.
      phi,           & ! Particle diameter ratio.
      rho_p,         & ! Particle density.
      mav,           &
      mtot,          &
      s_st             ! Critical Stokes number.

    m_inr => modes_h(iINR)

    do imod = 1, size(modes_f)
      mode => modes_f(imod) 
      if (mode%nspecies > 0) then
        do k = 1, kmax
          do j = 2, j1
            do i = 2, i1
              if (sed_qr(i,j,k) > 0.0_field_r) then
                xr = calc_xr(rho(k), qr(i,j,k), nr(i,j,k), 2.6E-10, 5.0E-6) 
                dvr = calc_dvr(xr)
                T = thl(i,j,k) * exn(k)
                dm = calc_median_diameter(mode%n(i,j,k), mode%q(i,j,k,:), &
                                          mode%rho, mode%sig_g)
                v = 9.65 - 9.8 * exp(-600 * dvr)
                rho_p = calc_mean_rho(mode%q(i,j,k,:), mode%rho)

                lbd = 2 * mu_a / (p(k) * sqrt(8 / (pi * rd * T)))
                csc = 1 + 2 * lbd / dvr * (1.257 + 0.4 * exp(-0.55 * dvr / lbd))

                re = 0.5 * rho(k) * dvr * v / mu_a
                db = kb * T * csc / (2 * pi * mu_a * dm)
                sc = mu_a / (rho(k) * db)
                tau = (rho_p - rho(k)) * dm**2 * csc / (18 * mu_a)
                st = 2 * tau * v / dvr
                phi = dm / dvr

                s_st = (1.2 + (log(1 + re) / 12)) / (1 + log(1 + re))

                E = 4 / (re * sc) * (1 + 0.4 * sqrt(re) * sc**(1.0_field_r / 3) + &
                                     0.16 * sqrt(re) * sqrt(sc)) &
                    + 4 * phi * (mu_a / mu_w + (1 + 2 * sqrt(re)) * phi) &
                    + (max(st - s_st, 0.0_field_r) &
                       / (st - s_st + 2.0_field_r / 3))**(1.5) * sqrt(rho_p / rhow)
              
                ! Tost et al. 2016
                E = max(min(E, 1.0_field_r), 0.0_field_r)
                gamma = 1.5 * E / (0.5 * dvr * 1E3) * (sed_qr(i,j,k) * rho(k))
                cnacc = 1 - exp(-delt*gamma)

                do s = 1, mode%nspecies
                  ts = mode%to_hydro%cnct(s,2)
                  mode%qp(i,j,k,s) = mode%qp(i,j,k,s) &
                                     - mode%q(i,j,k,s) * cnacc / delt
                  m_inr%qp(i,j,k,ts) = m_inr%qp(i,j,k,ts) &
                                     + mode%q(i,j,k,s) * cnacc / delt
                end do

                ! This is not right:
                mode%np(i,j,k) = mode%np(i,j,k) - cnacc * mode%n(i,j,k) / delt
              end if
            end do
          end do
        end do
      end if
    end do

  end subroutine aerosol_scavenging_rain

end module modaerosol