!> Adapter routines for connecting DALES host data to LCM-owned interfaces.
module modlcm_adapter
  use iso_fortran_env, only : real64
  use lcm_types, only : lcm_grid_t, lcm_init_grid

  implicit none

  private

  public :: build_lcm_grid_from_dales

contains

  subroutine build_lcm_grid_from_dales(grid)
    use modglobal, only : imax, jmax, kmax, itot, jtot, i1, j1, &
                          ih, jh, kh, dx, dy, dzf, zf, zh
    use modmpi, only : myidx, myidy, nprocx, nprocy, nbrwest,  &
                       nbreast, nbrsouth, nbrnorth, periods

    type(lcm_grid_t), intent(out) :: grid

    call lcm_init_grid(                                                   &
      grid=grid,                                                          &
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
      periodic_x=periods(1), periodic_y=periods(2))
  end subroutine build_lcm_grid_from_dales

end module modlcm_adapter
