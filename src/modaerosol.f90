module modaerosol

  use, intrinsic :: iso_fortran_env

  use modaerosol_mode_t,     only: aerosol_mode_t, hydrometeor_mode_t, mode_t, &
                                   mode_container_t, connect_modes
  use modaerosol_common,     only: maxspecies, maxmodes, iNUS, iAIS, iACS, &
                                   iCOS, iAII, iACI, iCOI, iINC, iINR, iSO4, &
                                   iSS, iPOM, iBC, iDU, calc_median_diameter
  use modaerosol_scavenging, only: init_scavenging, &
                                   aerosol_scavenging_rain_lut, &
                                   aerosol_scavenging_cloud_lut
  use modglobal,             only: ifnamopt, fname_options, &
                                   checknamelisterror, cexpnr, i1, j1, k1, ih, &
                                   jh, pi, nsv, rhow, kmax, rk3step, rd, pirhow, &
                                   timee, rk3step, cp, rlv
  use modfields,             only: sv0, svp, svm, ql0
  use modmicrodata,          only: qcmin, delt
  use modmpi,                only: myid, D_MPI_BCAST, commwrld, mpierr
  use modprecision,          only: field_r
  use modtimer,              only: timer_tic, timer_toc
  use modbulkmicro_data,     only: l_sb, qrmin, l_mur_cst, mur_cst
  use bulkmicro_sb,          only: calc_sed_qr_sb, calc_sed_nr_sb
  use bulkmicro_kk,          only: calc_sed_nr_kk, calc_sed_qr_kk
  use fortran_support,       only: nnml_output
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
  public :: aerosol_scavenging_cloud

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

  real(field_r), pointer :: &
    qr_spl(:,:,:),          & ! Rain water content at sub-timesteps.
    nr_spl(:,:,:),          & ! Rain number concentration at sub-timesteps.
    qa_spl(:,:,:,:)           ! Aerosol mass at sub-timesteps.

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
      write(nnml_output, NAMAEROSOL)
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
!      call modes(imod)%p%init(imod, mode_config(:,imod))
    end do

    ! Connect free aerosol modes to in-hydrometeor modes
    do imod = 1, maxmodes
      select type(mode => modes(imod)%p)
        class is (aerosol_mode_t)
          call connect_modes(mode, modes(iINC)%p, mode%to_hydro)
          !$acc enter data copyin(mode%to_hydro%cnct)
      end select
    end do

    allocate(sed_qr(2:i1,2:j1,k1), qlm(2:i1,2:j1,k1))

    allocate(qr_spl(2:i1,2:j1,1:k1), nr_spl(2:i1,2:j1,1:k1), &
             qa_spl(1:modes_h(iINR)%nspecies,2:i1,2:j1,1:k1))

    !$acc enter data create(qr_spl(2:i1,2:j1,1:k1), nr_spl(2:i1,2:j1,1:k1), &
    !$acc                   qa_spl(1:modes_h(iINR)%nspecies,2:i1,2:j1,1:k1))

    sed_qr(:,:,:) = 0
    qlm(:,:,:) = 0

    !$acc enter data copyin(sed_qr(2:i1,2:j1,1:k1), qlm(2:i1,2:j1,1:k1))

    call init_scavenging()

  end subroutine init_aerosol

  !> Prepares aerosol fields for microphysics calculations.
  subroutine aerosol_prepare()

    character(len=*), parameter :: &
      routine = modname//'/aerosol_prepare'

    integer :: i, j, k, imod

    call timer_tic(routine, 2)

    !$acc wait

    do imod = 1, maxmodes
!      call modes(imod)%p%prepare(sv0)
    end do

    if (rk3step == 3 .or. timee < 0.01) then
      !$acc parallel loop collapse(3) default(present)
      do k = 1, k1
        do j = 2, j1
          do i = 2, i1
            qlm(i,j,k) = ql0(i,j,k)
          end do
        end do
      end do
    end if

    !$acc wait

    call timer_toc(routine)

  end subroutine aerosol_prepare

  !> Copy out tendencies.
  subroutine aerosol_finish()

    character(len=*), parameter :: &
      routine = modname//'/aerosol_finish'

    integer :: &
      imod

    call timer_tic(routine, 2)

    !$acc wait

    do imod = 1, maxmodes

    !call modes(imod)%p%finish(svp, svm, delt)

    end do

    !$acc wait

    call timer_toc(routine)

  end subroutine aerosol_finish

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

    call timer_tic(routine, 2)

    m_ais => modes_f(iAIS)
    m_acs => modes_f(iACS)
    m_cos => modes_f(iCOS)
    m_inc => modes_h(iINC)
    
    !$acc parallel loop collapse(3) default(present) &
    !$acc private(dm, fn, n_act, w0, dncdt, fm, tend_n, tend_m, st)
    do k = 1, kmax
      do j = 2, j1
        do i = 2, i1
          if (ql(i,j,k) > qcmin) then
            if (m_ais%nspecies > 0) then
              dm = calc_median_diameter(m_ais%n(i,j,k), m_ais%q(:,i,j,k), &
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
              tend_m = max(0.0_field_r, fm * m_cos%q(s,i,j,k) / delt)
              st = m_cos%to_hydro%cnct(2,s)
              m_cos%qp(s,i,j,k) = m_cos%qp(s,i,j,k) - tend_m
              m_inc%qp(st,i,j,k) = m_inc%qp(st,i,j,k) + tend_m
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
                tend_m = max(0.0_field_r, fm * m_acs%q(s,i,j,k) / delt)
                st = m_acs%to_hydro%cnct(2,s)
                m_acs%qp(s,i,j,k) = m_acs%qp(s,i,j,k) - tend_m
                m_inc%qp(st,i,j,k) = m_inc%qp(st,i,j,k) + tend_m
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
                tend_m = max(0.0_field_r, fm * m_ais%q(s,i,j,k) / delt)
                st = m_ais%to_hydro%cnct(2,s)
                m_ais%qp(s,i,j,k) = m_ais%qp(s,i,j,k) - tend_m
                m_inc%qp(st,i,j,k) = m_inc%qp(st,i,j,k) + tend_m
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

    real(field_r) :: frac

    call timer_tic(routine, 2)

    m_inc => modes_h(iINC)
    m_inr => modes_h(iINR)

    !$acc parallel loop collapse(4) default(present) private(dqadt)
    do k = 1, kmax
      do j = 2, j1
        do i = 2, i1
          do s = 1, m_inc%nspecies
            if (qrp(i,j,k) > 0) then
              frac = min(max(qrp(i,j,k) / qc(i,j,k), 0.0_field_r), 1.0_field_r)
              dqadt = qrp(i,j,k) / qc(i,j,k) * m_inc%q(s,i,j,k)
              m_inc%qp(s,i,j,k) = m_inc%qp(s,i,j,k) - dqadt
              m_inr%qp(s,i,j,k) = m_inr%qp(s,i,j,k) + dqadt
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
      Dc = 1E-6 ! Diameter separating the accumulation and coarse modes.

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

    call timer_tic(routine, 2)

    m_acs => modes_f(iACS)
    m_cos => modes_f(iCOS)
    m_inr => modes_h(iINR)

    !$acc parallel loop collapse(3) default(present) &
    !$acc private(f_evp, eps, evapm, evapn, dn, dm, fn, fm)
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
          
            evapm(:) = eps * f_evp * m_inr%q(:,i,j,k) / delt
            evapn = max(0.0_field_r, -1 * nrp(i,j,k))

            ! Compute the median diameter of the resuspended aerosol.
            dn = calc_median_diameter(evapn, evapm(:), m_inr%rho, 1.5_field_r)
            dm = dn * exp(3 * log(1.5_field_r)**2)

            fn = 0.5_field_r * erfc(-log(dc/(dn + 1E-40)) &
                                    / (log(1.5_field_r) * sqrt(2.0_field_r)))
            fm = 0.5_field_r * erfc(-log(dc/(dm + 1E-40)) &
                                    / (log(1.5_field_r) * sqrt(2.0_field_r)))


            fn = min(max(fn, 0.0_field_r), 1.0_field_r)
            fm = min(max(fm, 0.0_field_r), 1.0_field_r)

            m_acs%np(i,j,k) = m_acs%np(i,j,k) + fn * evapn
            m_cos%np(i,j,k) = m_cos%np(i,j,k) + (1 - fn) * evapn

            do s = 1, m_inr%nspecies
              m_inr%qp(s,i,j,k) = m_inr%qp(s,i,j,k) - evapm(s)
              m_acs%qp(s,i,j,k) = m_acs%qp(s,i,j,k) + fm * evapm(s)
              m_cos%qp(s,i,j,k) = m_cos%qp(s,i,j,k) + (1 - fm) * evapm(s)
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
      Dc = 1E-6 ! Diameter separating the accumulation and coarse modes.

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

    call timer_tic(routine, 2)

    m_acs => modes_f(iACS)
    m_cos => modes_f(iCOS)
    m_inc => modes_h(iINC)

    !$acc parallel loop collapse(3) default(present) &
    !$acc private(f_evp, eps, evapm, evapn, dn, dm, fn, fm)
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
          
          evapm(:) = eps * f_evp * m_inc%q(:,i,j,k) / delt
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
            m_inc%qp(s,i,j,k) = m_inc%qp(s,i,j,k) - evapm(s)
            m_acs%qp(s,i,j,k) = m_acs%qp(s,i,j,k) + fm * evapm(s)
            m_cos%qp(s,i,j,k) = m_cos%qp(s,i,j,k) + (1 - fm) * evapm(s)
          end do
        end do
      end do
    end do

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

    call timer_tic(routine, 2)

    m_inr => modes_h(iINR)


    n_spl = ceiling(9.9 * delt / minval(dzf))
    dt_spl = delt / real(n_spl, kind=field_r)

    !$acc parallel loop collapse(3) default(present)
    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          qr_spl(i,j,k) = qr(i,j,k)
          nr_spl(i,j,k) = nr(i,j,k)
        end do
      end do
    end do

    !$acc parallel loop collapse(4) default(present)
    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          do s = 1, m_inr%nspecies
            qa_spl(s,i,j,k) = m_inr%q(s,i,j,k)
          end do
        end do
      end do
    end do

    do ts = 1, n_spl
      !$acc parallel loop collapse(3) default(present) private(sed_nr)
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
              !$acc atomic update
              qr_spl(i,j,k) = qr_spl(i,j,k) - sed_qr(i,j,k) * dt_spl &
                              / (dzf(k) * rho(k))
              !$acc atomic update
              nr_spl(i,j,k) = nr_spl(i,j,k) - sed_nr * dt_spl / dzf(k)
              if (k > 1) then
                !$acc atomic update
                qr_spl(i,j,k-1) = qr_spl(i,j,k-1) + sed_qr(i,j,k) * dt_spl &
                                  / (dzf(k-1) * rho(k-1))
                !$acc atomic update
                nr_spl(i,j,k-1) = nr_spl(i,j,k-1) + sed_nr * dt_spl / dzf(k-1)
              end if
              do s = 1, m_inr%nspecies
                !$acc atomic update
                qa_spl(s,i,j,k) = qa_spl(s,i,j,k) - sed_qr(i,j,k) / qr_spl(i,j,k) &
                                  * qa_spl(s,i,j,k) * dt_spl / (dzf(k) * rho(k))
                if (k > 1) then
                  !$acc atomic update
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

    !$acc parallel loop collapse(4) default(present)
    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          do s = 1, m_inr%nspecies
            m_inr%qp(s,i,j,k) = m_inr%qp(s,i,j,k) + &
                               (qa_spl(s,i,j,k) - m_inr%q(s,i,j,k)) / delt
          end do
        end do
      end do
    end do
    
    call timer_toc(routine)

  end subroutine aerosol_sedimentation_rain
  
  !> Compute scavenging of aerosols by rain drops.
  subroutine aerosol_scavenging_rain(qr, nr, rho, delt)

    real(field_r), intent(in) :: qr(2:,2:,:) !< Rain water content [kg kg-1].
    real(field_r), intent(in) :: nr(2:,2:,:) !< Rain number concentration [m-3].
    real(field_r), intent(in) :: rho(:)      !< Air density [kg m-3].
    real(field_r), intent(in) :: delt        !< Time step size [s].

    integer :: imod

    do imod = 1, size(modes_f)
      call aerosol_scavenging_rain_lut(qr, nr, rho, delt, modes_f(imod), &
                                       modes_h(iINR))
    end do

  end subroutine aerosol_scavenging_rain

  !> Compute scavenging of aerosols by cloud droplets.
  subroutine aerosol_scavenging_cloud(ql, nc, rho, delt)

    real(field_r), intent(in) :: ql(2-ih:,2-jh:,:) !< Cloud water content [kg kg-1].
    real(field_r), intent(in) :: nc(2:,2:,:)       !< Cloud droplet number concentration [m-3].
    real(field_r), intent(in) :: rho(:)            !< Air density [kg m-3].
    real(field_r), intent(in) :: delt              !< Time step size [s].

    integer :: imod

    do imod = 1, size(modes_f)
      call aerosol_scavenging_cloud_lut(ql, nc, rho, delt, modes_f(imod), modes_h(iINC))
    end do

  end subroutine aerosol_scavenging_cloud

end module modaerosol
