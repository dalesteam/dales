!> Adapter routines for connecting DALES host data to LCM-owned interfaces.
#ifdef USE_LCM
#if FIELD_PRECISION != 64
#error "LCM Option A pointer attachment requires DALES FIELD_PRECISION=64"
#endif
module modlcm_adapter
  use iso_fortran_env, only : real64
  use lcm_host_interface, only : lcm_config_t, lcm_grid_t, lcm_init_config, &
                                 lcm_init_grid, lcm_set_config, lcm_set_grid, &
                                 lcm_attach_fields, lcm_initialize, lcm_advance
  use modfields, only : u0, v0, w0, tmp0, dse0, qv0, presf, rhof
  use modglobal, only : imax, jmax, kmax, itot, jtot, i1, j1, &
                        ih, jh, kh, dx, dy, dzf, zf, zh, rdt
  use modlcm_namelist, only : lcm_apply_namelist_config
  use modmpi, only : myidx, myidy, nprocx, nprocy, nbrwest,  &
                     nbreast, nbrsouth, nbrnorth, periods

  implicit none

  private

  public :: init_lcm
  public :: lcm_microphysics

  type(lcm_grid_t), save :: lcm_grid
  type(lcm_config_t), save :: lcm_config

contains

  subroutine init_lcm()
    call lcm_init_config(lcm_config)
    call lcm_apply_namelist_config(lcm_config)
    call lcm_set_config(lcm_config)

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
      subdomain_x_index=myidx, subdomain_y_index=myidy,                   &
      n_subdomains_x=nprocx, n_subdomains_y=nprocy,                       &
      global_i_start=myidx * imax + 1, global_j_start=myidy * jmax + 1,   &
      west_rank=nbrwest, east_rank=nbreast,                               &
      south_rank=nbrsouth, north_rank=nbrnorth,                           &
      periodic_x=periods(1), periodic_y=periods(2),                       &
      y_dimension=1)

    call lcm_set_grid(lcm_grid)
    call lcm_attach_fields(                                             &
      u=u0(2:i1, 2:j1, 1:kmax),                                         &
      v=v0(2:i1, 2:j1, 1:kmax),                                         &
      w=w0(2:i1, 2:j1, 1:kmax),                                         &
      static_energy=dse0(2:i1, 2:j1, 1:kmax),                           &
      temperature=tmp0(2:i1, 2:j1, 1:kmax),                             &
      qv=qv0(2:i1, 2:j1, 1:kmax),                                       &
      pressure=presf(1:kmax),                                           &
      density=rhof(1:kmax))
    call lcm_initialize()
  end subroutine init_lcm

  subroutine lcm_microphysics()
    call lcm_advance(real(rdt, kind=real64))
  end subroutine lcm_microphysics

end module modlcm_adapter
#endif
