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
                                 lcm_initialize, lcm_advance
  use modfields, only : u0, v0, w0, tmp0, dse0, qv0, presf, rhof
  use modglobal, only : imax, jmax, kmax, itot, jtot, i1, j1, &
                        ih, jh, kh, dx, dy, dzf, zf, zh, rdt, cp
  use modlcm_namelist, only : lcm_apply_namelist_config
  use modmpi, only : comm3d, myidx, myidy
#ifdef LCM_VALIDATION_CHECKS
  use modlcm_validation_debug, only : print_lcm_validation_checks
#endif

  implicit none

  private

  public :: configure_lcm_runtime
  public :: init_lcm
  public :: lcm_microphysics

  type(lcm_grid_t), save :: lcm_grid
  type(lcm_parallel_t), save :: lcm_parallel
  type(lcm_config_t), save :: lcm_config
  real(real64), allocatable, target, save :: lcm_static_energy(:,:,:)

contains

  subroutine configure_lcm_runtime()
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
    call update_lcm_static_energy_temperature_units()

    call lcm_attach_fields(                                             &
      u=u0(2:i1+1, 2:j1, 1:kmax),                                       &
      v=v0(2:i1, 2:j1+1, 1:kmax),                                       &
      w=w0(2:i1, 2:j1, 1:kmax+1),                                       &
      static_energy=lcm_static_energy,                                   &
      temperature=tmp0(2:i1, 2:j1, 1:kmax),                             &
      qv=qv0(2:i1, 2:j1, 1:kmax),                                       &
      pressure=presf(1:kmax),                                           &
      density=rhof(1:kmax))
    call lcm_initialize()
#ifdef LCM_VALIDATION_CHECKS
    call print_lcm_validation_checks(lcm_config)
#endif
  end subroutine init_lcm

  subroutine lcm_microphysics()
    ! Called once after the final RK integration, boundaries and thermodynamics
    ! by lcm_after_dynamics; rdt is the full atmospheric timestep.
    call update_lcm_static_energy_temperature_units()
    call lcm_advance(real(rdt, kind=real64))
  end subroutine lcm_microphysics

  subroutine update_lcm_static_energy_temperature_units()
    if (.not. allocated(lcm_static_energy)) then
      allocate(lcm_static_energy(imax, jmax, kmax))
    end if

    ! LCM follows SAM-LCM and stores static energy in temperature units.
    lcm_static_energy = real(dse0(2:i1, 2:j1, 1:kmax), real64) / &
                        real(cp, real64)
  end subroutine update_lcm_static_energy_temperature_units

end module modlcm_adapter
#endif
