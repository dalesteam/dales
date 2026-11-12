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
                                 lcm_attach_thermodynamic_history, lcm_initialize, lcm_advance
  use modboundary, only : lboundopen
  use modfields, only : u0, v0, w0, tmp0, dse0, qv0, qt0, presf, rhof, sv0, svm
  use modglobal, only : imax, jmax, kmax, itot, jtot, i1, j1, k1, &
                        ih, jh, kh, dx, dy, dzf, zf, zh, rdt, cp, lopenbc, lwarmstart, iadv_sv, iadv_null, nsv
  use modlcm_namelist, only : lcm_apply_namelist_config
  use modlogging, only : finish
  use modmpi, only : comm3d, myidx, myidy, excjs
  use modtracers, only : add_tracer
#ifdef LCM_VALIDATION_CHECKS
  use modlcm_validation_debug, only : print_lcm_validation_checks
#endif

  implicit none

  private

  public :: configure_lcm_runtime
  public :: init_lcm
  public :: advance_lcm_step

  type(lcm_grid_t), save :: lcm_grid
  type(lcm_parallel_t), save :: lcm_parallel
  type(lcm_config_t), save :: lcm_config
  real(real64), allocatable, target, save :: lcm_static_energy(:,:,:)
  real(real64), allocatable, target, save :: lcm_vapor_mixing_ratio(:,:,:)
  real(real64), allocatable, target, save :: lcm_dry_air_density(:,:,:)
  integer, save :: i_t_prev, i_p_prev, i_qv_prev

contains

  subroutine configure_lcm_runtime()
    integer :: nsv_before

    if (lopenbc .or. lwarmstart) then
      call finish('modlcm_adapter', 'LCM history tracers require a periodic cold start.')
    end if
    if (iadv_sv == iadv_null) then
      call finish('modlcm_adapter', 'LCM history tracers require scalar advection.')
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
    if (lboundopen) call finish('modlcm_adapter', 'LCM history tracers require periodic scalar boundaries.')
    call update_lcm_atmosphere()

    call lcm_attach_fields(                                             &
      u=u0(2:i1+1, 2:j1, 1:kmax),                                       &
      v=v0(2:i1, 2:j1+1, 1:kmax),                                       &
      w=w0(2:i1, 2:j1, 1:kmax+1),                                       &
      static_energy=lcm_static_energy,                                   &
      temperature=tmp0(2:i1, 2:j1, 1:kmax),                             &
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
#ifdef LCM_VALIDATION_CHECKS
    call print_lcm_validation_checks(lcm_config)
#endif
  end subroutine init_lcm

  subroutine advance_lcm_step()
    ! Called once after the final RK integration, boundaries and thermodynamics
    ! by lcm_microphysics; rdt is the full atmospheric timestep.
    call update_lcm_atmosphere()
    call lcm_advance(real(rdt, kind=real64))
    call refresh_lcm_history()
  end subroutine advance_lcm_step

  subroutine refresh_lcm_history()
    integer :: k, n, tracer
    integer :: history_tracers(3)

    ! The saved atmosphere is the starting state for the next RK timestep.
    sv0(2:i1, 2:j1, 1:k1, i_t_prev) = tmp0(2:i1, 2:j1, 1:k1)
    sv0(2:i1, 2:j1, 1:k1, i_qv_prev) = lcm_vapor_mixing_ratio
    do k = 1, k1
      sv0(2:i1, 2:j1, k, i_p_prev) = presf(k)
    end do
    history_tracers = [i_t_prev, i_p_prev, i_qv_prev]
    do n = 1, size(history_tracers)
      tracer = history_tracers(n)
      call excjs(sv0(:,:,:,tracer), 2, i1, 2, j1, 1, k1, ih, jh)
      svm(:,:,:,tracer) = sv0(:,:,:,tracer)
    end do
  end subroutine refresh_lcm_history

  subroutine update_lcm_atmosphere()
    integer :: k

    if (.not. allocated(lcm_static_energy)) then
      allocate(lcm_static_energy(imax, jmax, kmax))
      allocate(lcm_vapor_mixing_ratio(imax, jmax, k1))
      allocate(lcm_dry_air_density(imax, jmax, kmax))
    end if

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

end module modlcm_adapter
#endif
