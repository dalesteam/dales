!> Adapter routines for connecting DALES host data to LCM-owned interfaces.
#ifdef USE_LCM
#if FIELD_PRECISION != 64
#error "LCM Option A pointer attachment requires DALES FIELD_PRECISION=64"
#endif
module modlcm_adapter
  use iso_fortran_env, only : real64
  use lcm_host_interface, only : lcm_config_t, lcm_grid_t, lcm_parallel_t,   &
                                 lcm_init_config, lcm_init_grid,              &
                                 lcm_init_parallel, lcm_set_config,           &
                                 lcm_set_grid, lcm_set_parallel,              &
                                 lcm_configure_runtime, lcm_attach_fields,    &
                                 lcm_attach_thermodynamic_history, lcm_initialize, lcm_advance, &
                                 lcm_get_liquid_water_mass, lcm_get_condensation_mass_change, &
                                 lcm_pass_timer_routines
  use modboundary, only : lboundopen, ksp, tsc
  use modfields, only : u0, v0, w0, tmp0, dse0, qv0, qt0, qtm, thlm, ql0, &
                         presf, rhof, rhobf, sv0, svm, svp, sv0av, exnf, thlp, &
                         u0av, v0av, whls, dthldxls, dthldyls, dthldtls, &
                         dqtdxls, dqtdyls, dqtdtls
  use modglobal, only : imax, jmax, kmax, itot, jtot, i1, j1, k1, &
                        ih, jh, kh, dx, dy, dzf, zf, zh, rdt, cp, lopenbc, lwarmstart, &
                        iadv_sv, iadv_null, nsv, rlv, grav, ladaptive, dtmax, igrw_damp
  use modlcm_namelist, only : lcm_apply_namelist_config, lcm_piggyback
  use modsurfdata, only : thlflux, qtflux, svflux
  use modnudge, only : lnudge
  use ieee_arithmetic, only : ieee_is_finite
  use modlogging, only : finish
  use modmpi, only : comm3d, myidx, myidy, excjs, D_MPI_ALLREDUCE, MPI_SUM
  use modtracers, only : add_tracer
  use modtimer, only : timer_tic, timer_toc

  implicit none

  private

  public :: configure_lcm_runtime
  public :: init_lcm
  public :: advance_lcm_step
  public :: prepare_lcm_piggyback, force_lcm_piggyback

  type(lcm_grid_t), save :: lcm_grid
  type(lcm_parallel_t), save :: lcm_parallel
  type(lcm_config_t), save :: lcm_config
  real(real64), allocatable, target, save :: lcm_static_energy(:,:,:)
  real(real64), allocatable, target, save :: lcm_w(:,:,:)
  real(real64), allocatable, target, save :: lcm_vapor_mixing_ratio(:,:,:)
  real(real64), allocatable, target, save :: lcm_dry_air_density(:,:,:)
  integer, save :: i_t_prev, i_p_prev, i_qv_prev
  integer, save :: i_theta, i_vapor
  real(real64), allocatable, target, save :: lcm_temperature(:,:,:)
  real(real64), allocatable, save :: liquid_mass(:,:,:), phase_change(:,:,:)
  logical, save :: piggyback_initialized = .false.

contains

  subroutine configure_lcm_runtime()
    integer :: nsv_before

    if (lopenbc .or. lwarmstart) then
      call finish('modlcm_adapter', 'LCM history tracers require a periodic cold start.')
    end if
    if (iadv_sv == iadv_null) then
      call finish('modlcm_adapter', 'LCM history tracers require scalar advection.')
    end if
    if (lcm_piggyback) then
#ifdef DALES_GPU
      call finish('modlcm_adapter', 'LCM passenger field updates currently require a CPU build.')
#endif
      if (ladaptive) call finish('modlcm_adapter', 'LCM piggybacking requires a constant timestep.')
      nsv_before = nsv
      call add_tracer('lcm_theta', 'LCM passenger potential temperature', 'K', isv=i_theta)
      call add_tracer('lcm_qv', 'LCM passenger vapor per reference-air mass', 'kg kg-1', isv=i_vapor)
      if (nsv /= nsv_before + 2) &
        call finish('modlcm_adapter', 'The lcm_theta and lcm_qv tracer names are reserved.')
    end if
    nsv_before = nsv
    call add_tracer('lcm_t_prev', 'LCM transported previous temperature', 'K', isv=i_t_prev)
    call add_tracer('lcm_p_prev', 'LCM transported previous pressure', 'Pa', isv=i_p_prev)
    call add_tracer('lcm_qv_prev', 'LCM transported previous vapor mixing ratio', 'kg kg-1', isv=i_qv_prev)
    if (nsv /= nsv_before + 3) then
      call finish('modlcm_adapter', 'The lcm_t_prev, lcm_p_prev, and lcm_qv_prev tracer names are reserved.')
    end if

    call lcm_init_config(lcm_config)
    call lcm_apply_namelist_config(lcm_config)
    if (lcm_piggyback) then
      if (lcm_config%switch_lem_micro .or. lcm_config%switch_ice .or. &
          lcm_config%switch_collection .or. lcm_config%switch_sedimentation .or. &
          lcm_config%switch_particle_source .or. lcm_config%switch_sgs_velocities) &
        call finish('modlcm_adapter', &
          'LCM piggybacking requires cell-mean liquid growth without sources, collection, sedimentation or SGS motion.')
    end if
    call lcm_set_config(lcm_config)

    ! LCM derives the complete Cartesian topology from DALES's existing
    ! communicator; the adapter does not duplicate neighbor metadata.
    call lcm_init_parallel(lcm_parallel, comm3d)
    call lcm_set_parallel(lcm_parallel)

    call lcm_init_grid(                                                   &
      grid=lcm_grid,                                                      &
      nx_local=imax, ny_local=jmax, nz=kmax,                              &
      nx_global=itot, ny_global=jtot,                                     &
      i_start=2, i_end=i1, j_start=2, j_end=j1,                           &
      dx=real(dx, kind=real64), dy=real(dy, kind=real64),                 &
      z_center=real(zf(1:kmax), kind=real64),                             &
      z_face=real(zh(1:kmax+1), kind=real64),                             &
      dz=real(dzf(1:kmax), kind=real64),                                  &
      halo_x=ih, halo_y=jh, halo_z=kh,                                    &
      global_i_start=myidx * imax + 1, global_j_start=myidy * jmax + 1,   &
      y_dimension=1)

    call lcm_set_grid(lcm_grid)
    call lcm_configure_runtime()
  end subroutine configure_lcm_runtime

  subroutine init_lcm()
    call lcm_pass_timer_routines(timer_tic, timer_toc)

    if (lboundopen) call finish('modlcm_adapter', 'LCM history tracers require periodic scalar boundaries.')
    call update_lcm_atmosphere()

    call lcm_attach_fields(                                             &
      u=u0(2:i1+1, 2:j1, 1:kmax),                                       &
      v=v0(2:i1, 2:j1+1, 1:kmax),                                       &
      w=lcm_w,                                                         &
      static_energy=lcm_static_energy,                                   &
      temperature=lcm_temperature(:,:,1:kmax),                           &
      qv=lcm_vapor_mixing_ratio(:,:,1:kmax),                            &
      pressure=presf(1:kmax),                                           &
      density=rhof(1:kmax),                                            &
      dry_air_density=lcm_dry_air_density)
    call refresh_lcm_history()
    call lcm_attach_thermodynamic_history( &
      temperature=sv0(2:i1, 2:j1, 1:kmax, i_t_prev), &
      pressure=sv0(2:i1, 2:j1, 1:kmax, i_p_prev), &
      qv=sv0(2:i1, 2:j1, 1:kmax, i_qv_prev))
    call lcm_initialize()
    if (lcm_piggyback) then
      allocate(liquid_mass(imax,jmax,kmax),phase_change(imax,jmax,kmax))
      call initialize_piggyback_state()
    end if
  end subroutine init_lcm

  subroutine advance_lcm_step()
    ! Called once after the final RK integration, boundaries and thermodynamics
    ! by lcm_microphysics; rdt is the full atmospheric timestep.
    call update_lcm_atmosphere()
    if (lcm_piggyback) then
      if (abs(rdt-dtmax)>1.e-12_real64) call finish('modlcm_adapter', 'LCM did not receive the fixed timestep.')
      ! Particle migration changes local loading, so refresh humidity units and
      ! dry-air mass after migration, immediately before condensation.
      call lcm_advance(real(rdt,real64), refresh_atmosphere=update_piggyback_atmosphere)
      call apply_piggyback_condensation()
    else
      call lcm_advance(real(rdt,real64))
    end if
    call refresh_lcm_history()
  end subroutine advance_lcm_step

  subroutine refresh_lcm_history()
    integer :: k, n
    integer :: history_tracers(3)

    ! The saved atmosphere is the starting state for the next RK timestep.
    sv0(2:i1, 2:j1, 1:k1, i_t_prev) = lcm_temperature
    sv0(2:i1, 2:j1, 1:k1, i_qv_prev) = lcm_vapor_mixing_ratio
    do k = 1, k1
      sv0(2:i1, 2:j1, k, i_p_prev) = presf(k)
    end do
    if (lcm_piggyback .and. piggyback_initialized) then
      call sync_tracer(i_theta)
      call sync_tracer(i_vapor)
    end if
    history_tracers = [i_t_prev, i_p_prev, i_qv_prev]
    do n = 1, size(history_tracers)
      call sync_tracer(history_tracers(n))
    end do
    if (lcm_piggyback) call update_tracer_means([i_theta, i_vapor, i_t_prev, i_p_prev, i_qv_prev])
  end subroutine refresh_lcm_history

  subroutine update_lcm_atmosphere()
    integer :: k

    if (.not. allocated(lcm_static_energy)) then
      allocate(lcm_static_energy(imax, jmax, kmax))
      allocate(lcm_temperature(imax,jmax,k1),lcm_w(imax,jmax,k1))
      allocate(lcm_vapor_mixing_ratio(imax, jmax, k1))
      allocate(lcm_dry_air_density(imax, jmax, kmax))
    end if

    ! Prescribed large-scale motion is applied separately from w0 to DALES
    ! scalars. Give particles the same motion at the shared vertical faces.
    do k = 1, kmax+1
      lcm_w(:,:,k) = w0(2:i1, 2:j1, k) + whls(k)
    end do
    if (lcm_piggyback .and. piggyback_initialized) then
      call update_piggyback_atmosphere()
      return
    end if
    lcm_temperature = tmp0(2:i1,2:j1,1:k1)

    ! DALES water fields are per unit moist-air mass; LCM vapor is per unit
    ! dry-air mass. Use total water, including liquid, for the dry-air fraction.
    lcm_vapor_mixing_ratio = qv0(2:i1, 2:j1, 1:k1) / &
                            (1.0_real64 - qt0(2:i1, 2:j1, 1:k1))
    do k = 1, kmax
      lcm_dry_air_density(:,:,k) = rhof(k) * &
                                 (1.0_real64 - qt0(2:i1, 2:j1, k))
    end do

    ! LCM follows SAM-LCM and stores static energy in temperature units.
    lcm_static_energy = real(dse0(2:i1, 2:j1, 1:kmax), real64) / &
                        real(cp, real64)
  end subroutine update_lcm_atmosphere


  !> Same imposed surface/radiative forcing; no bulk microphysics tendencies.
  subroutine prepare_lcm_piggyback()
    if (.not. lcm_piggyback) return
    ! initnudge reads NAMNUDGE after init_lcm, so check the loaded setting here.
    if (lnudge) call finish('modlcm_adapter', 'Passenger nudging is not implemented; disable lnudge.')
    ! Called after radiation and surface, before advection/diffusion. At this
    ! point thlp contains the imposed radiative temperature tendency only.
    svflux(:,:,i_theta)=thlflux
    svflux(:,:,i_vapor)=qtflux
    svp(2:i1,2:j1,1:kmax,i_theta)=svp(2:i1,2:j1,1:kmax,i_theta)+thlp(2:i1,2:j1,1:kmax)
  end subroutine prepare_lcm_piggyback

  subroutine force_lcm_piggyback()
    integer :: k
    if (.not. lcm_piggyback) return
    ! DALES lstend already applies subsidence to each tracer's own profile.
    ! Add only the common prescribed horizontal and time-dependent sources.
    do k=1,kmax
      svp(2:i1,2:j1,k,i_theta)=svp(2:i1,2:j1,k,i_theta) &
        -u0av(k)*dthldxls(k)-v0av(k)*dthldyls(k)+dthldtls(k)
      svp(2:i1,2:j1,k,i_vapor)=svp(2:i1,2:j1,k,i_vapor) &
        -u0av(k)*dqtdxls(k)-v0av(k)*dqtdyls(k)+dqtdtls(k)
    end do
    ! The driver damps thermal perturbations in its upper sponge layer;
    ! ordinary passive tracers do not receive this operator automatically.
    if (igrw_damp>=1 .and. igrw_damp<=3) then
      do k=ksp,kmax
        svp(2:i1,2:j1,k,i_theta)=svp(2:i1,2:j1,k,i_theta) &
          -(sv0(2:i1,2:j1,k,i_theta)-sv0av(k,i_theta))*tsc(k)
        svp(2:i1,2:j1,k,i_vapor)=svp(2:i1,2:j1,k,i_vapor) &
          -(sv0(2:i1,2:j1,k,i_vapor)-sv0av(k,i_vapor))*tsc(k)
      end do
    end if
  end subroutine force_lcm_piggyback

  ! Passenger scalar budgets use the same anelastic reference-air mass
  ! (rhobf*cell_volume) as DALES advection and subgrid diffusion.
  subroutine initialize_piggyback_state()
    integer :: k
    real(real64) :: air_mass
    call lcm_get_liquid_water_mass(liquid_mass)
    do k=1,k1
      sv0(2:i1,2:j1,k,i_theta)=tmp0(2:i1,2:j1,k)/exnf(k)
      sv0(2:i1,2:j1,k,i_vapor)=qv0(2:i1,2:j1,k)
    end do
    do k=1,kmax
      air_mass=rhobf(k)*dx*dy*dzf(k)
      ! Match the bulk driver's water and moist enthalpy, including the
      ! initialized particle water. Condensation changes theta, not theta_l.
      sv0(2:i1,2:j1,k,i_vapor)=qt0(2:i1,2:j1,k)-liquid_mass(:,:,k)/air_mass
      sv0(2:i1,2:j1,k,i_theta)=(tmp0(2:i1,2:j1,k) &
        +rlv/cp*(liquid_mass(:,:,k)/air_mass-ql0(2:i1,2:j1,k)))/exnf(k)
    end do
    piggyback_initialized=.true.
    call update_piggyback_atmosphere()
    call refresh_lcm_history()
    ! A cold start can perturb DALES's two RK levels independently. Match
    ! each level rather than replacing the passenger base state by level 0.
    do k=1,kmax
      air_mass=rhobf(k)*dx*dy*dzf(k)
      svm(2:i1,2:j1,k,i_vapor)=qtm(2:i1,2:j1,k)-liquid_mass(:,:,k)/air_mass
      svm(2:i1,2:j1,k,i_theta)=thlm(2:i1,2:j1,k)+rlv/cp/exnf(k)*liquid_mass(:,:,k)/air_mass
      svm(2:i1,2:j1,k,i_t_prev)=svm(2:i1,2:j1,k,i_theta)*exnf(k)
      svm(2:i1,2:j1,k,i_qv_prev)=svm(2:i1,2:j1,k,i_vapor)/(1.0_real64-qtm(2:i1,2:j1,k))
    end do
    call excjs(svm(:,:,:,i_theta),2,i1,2,j1,1,k1,ih,jh)
    call excjs(svm(:,:,:,i_vapor),2,i1,2,j1,1,k1,ih,jh)
    call excjs(svm(:,:,:,i_t_prev),2,i1,2,j1,1,k1,ih,jh)
    call excjs(svm(:,:,:,i_qv_prev),2,i1,2,j1,1,k1,ih,jh)
  end subroutine initialize_piggyback_state

  subroutine update_piggyback_atmosphere()
    integer :: k
    real(real64) :: air_mass
    real(real64) :: dry_fraction(imax,jmax)
    call lcm_get_liquid_water_mass(liquid_mass)
    do k=1,k1
      lcm_temperature(:,:,k)=sv0(2:i1,2:j1,k,i_theta)*exnf(k)
      if(k<=kmax)then
        air_mass=rhobf(k)*dx*dy*dzf(k)
        dry_fraction=1.0_real64-sv0(2:i1,2:j1,k,i_vapor)-liquid_mass(:,:,k)/air_mass
        lcm_dry_air_density(:,:,k)=rhobf(k)*dry_fraction
      else
        dry_fraction=1.0_real64-sv0(2:i1,2:j1,k,i_vapor)
      endif
      if (.not. all(ieee_is_finite(dry_fraction)) .or. &
          any(dry_fraction <= 0.0_real64) .or. any(sv0(2:i1,2:j1,k,i_vapor) < 0.0_real64)) &
        call finish('modlcm_adapter','Invalid passenger water state; no clipping applied.')
      lcm_vapor_mixing_ratio(:,:,k)=sv0(2:i1,2:j1,k,i_vapor)/dry_fraction
      if(k<=kmax)lcm_static_energy(:,:,k)=lcm_temperature(:,:,k)+grav/cp*zf(k)
    end do
    if (.not. all(ieee_is_finite(lcm_temperature(:,:,1:kmax))) .or. &
        any(lcm_temperature(:,:,1:kmax) <= 0.0_real64)) &
      call finish('modlcm_adapter', 'Invalid passenger absolute temperature.')
  end subroutine update_piggyback_atmosphere

  subroutine apply_piggyback_condensation()
    integer :: k
    real(real64) :: air_mass
    call lcm_get_condensation_mass_change(phase_change)
    do k=1,kmax
      air_mass=rhobf(k)*dx*dy*dzf(k)
      sv0(2:i1,2:j1,k,i_vapor)=sv0(2:i1,2:j1,k,i_vapor)-phase_change(:,:,k)/air_mass
      sv0(2:i1,2:j1,k,i_theta)=sv0(2:i1,2:j1,k,i_theta) &
        +rlv/cp/exnf(k)*phase_change(:,:,k)/air_mass
    end do
    call update_piggyback_atmosphere()
  end subroutine apply_piggyback_condensation

  subroutine sync_tracer(tracer)
    integer, intent(in) :: tracer
    call excjs(sv0(:,:,:,tracer), 2, i1, 2, j1, 1, k1, ih, jh)
    svm(:,:,:,tracer) = sv0(:,:,:,tracer)
  end subroutine sync_tracer

  subroutine update_tracer_means(tracers)
    integer, intent(in) :: tracers(:)
    integer :: k, n, ierr
    real(real64) :: local_mean(k1,size(tracers)), global_mean(k1,size(tracers))

    ! Feedback occurs after DALES computed its slab averages. Refresh all five
    ! means for the next top-boundary/sponge call, using one MPI reduction.
    do n = 1, size(tracers)
      do k = 1, k1
        local_mean(k,n) = sum(sv0(2:i1,2:j1,k,tracers(n))) / real(itot*jtot,real64)
      end do
    end do
    call D_MPI_ALLREDUCE(local_mean, global_mean, size(local_mean), MPI_SUM, comm3d, ierr)
    do n = 1, size(tracers)
      sv0av(:,tracers(n)) = global_mean(:,n)
    end do
  end subroutine update_tracer_means

end module modlcm_adapter
#endif
