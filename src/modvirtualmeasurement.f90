!> Virtual measurement output at user-selected (x,y) points in one NetCDF file.
module modvirtualmeasurement

  use fortran_support,   only: nnml_output
  use modglobal,         only: ifnamopt, fname_options, checknamelisterror, &
                               itot, jtot, imax, jmax, &
                               x0, y0, dx, dy, cu, cv, cp, rlv, fkar, dtav_glob
  use modmpi,            only: myid, myidx, myidy, comm3d, mpierr, d_mpi_bcast, &
                               nprocs, D_MPI_ALLREDUCE, mpi_sum
  use modlogging,        only: finish
  use modlsm,            only: f1, f2b, tile, nlu
  use modnetcdf_file_t,  only: multi_timeseries_file_t
  use modprecision,      only: field_r
  use modsurfdata,       only: isurf
  use modstat_nc_files,  only: add_output_file, is_sampling_timestep

  implicit none

  private

  character(len=*), parameter :: modname = 'modvirtualmeasurement'
  integer, parameter :: max_points = 999 !< Max number of virtual measurement points.

  public :: initvirtualmeasurement
  public :: virtualmeasurement
  public :: exitvirtualmeasurement

  logical :: lvirtualmeasurement = .false.
  real    :: dtav !< Sampling interval [s].
  integer :: npoints = 0
  integer :: x_idx(max_points) = 0
  integer :: y_idx(max_points) = 0

  type(multi_timeseries_file_t) :: ofile
  integer :: ofile_id
  integer, allocatable :: i_local(:), j_local(:), global_index(:)

contains

  subroutine initvirtualmeasurement

    implicit none

    integer :: ierr, ip, il, nlocal, rank_offset
    integer :: x_idx_global, y_idx_global
    real(field_r) :: locx_point, locy_point
    real(field_r), allocatable :: locx(:), locy(:), locx_all(:), locy_all(:)
    integer, allocatable :: nlocal_send(:), nlocal_all_arr(:)

    character(len=*), parameter :: routine = modname//'/initvirtualmeasurement'

    namelist /NAMVIRTUALMEASUREMENT/ lvirtualmeasurement, dtav, npoints, x_idx, y_idx

    npoints = 0
    x_idx = 0
    y_idx = 0
    dtav = dtav_glob

    if (myid == 0) then
      open(ifnamopt, file=fname_options, status='old', iostat=ierr)
      read(ifnamopt, NAMVIRTUALMEASUREMENT, iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMVIRTUALMEASUREMENT')
      write(nnml_output, NAMVIRTUALMEASUREMENT)
      close(ifnamopt)
    end if

    call d_mpi_bcast(lvirtualmeasurement, 1, 0, comm3d, mpierr)
    call d_mpi_bcast(dtav, 1, 0, comm3d, mpierr)
    call d_mpi_bcast(npoints, 1, 0, comm3d, mpierr)
    call d_mpi_bcast(x_idx, max_points, 0, comm3d, mpierr)
    call d_mpi_bcast(y_idx, max_points, 0, comm3d, mpierr)

    if (.not. lvirtualmeasurement) return

    if (npoints < 1 .or. npoints > max_points) then
      call finish(routine, 'NAMVIRTUALMEASUREMENT: npoints must be between 1 and ', max_points)
    end if

    if (dtav <= 0.0_field_r) then
      call finish(routine, 'NAMVIRTUALMEASUREMENT: dtav must be > 0')
    end if

    do ip = 1, npoints
      if (x_idx(ip) < 1 .or. x_idx(ip) > itot) then
        call finish(routine, 'NAMVIRTUALMEASUREMENT: x_idx out of range at entry ', ip)
      end if
      if (y_idx(ip) < 1 .or. y_idx(ip) > jtot) then
        call finish(routine, 'NAMVIRTUALMEASUREMENT: y_idx out of range at entry ', ip)
      end if
    end do

    nlocal = 0
    do ip = 1, npoints
      x_idx_global = x_idx(ip)
      y_idx_global = y_idx(ip)
      if (x_idx_global >= myidx*imax + 1 .and. x_idx_global <= myidx*imax + imax .and. &
          y_idx_global >= myidy*jmax + 1 .and. y_idx_global <= myidy*jmax + jmax) then
        nlocal = nlocal + 1
      end if
    end do

    allocate(i_local(nlocal), j_local(nlocal), global_index(nlocal))
    allocate(locx(nlocal), locy(nlocal))
    allocate(locx_all(npoints), locy_all(npoints))

    ! here we determine the global index of each local point, and the local x,y coordinates of each point
    nlocal = 0
    do ip = 1, npoints
      x_idx_global = x_idx(ip)
      y_idx_global = y_idx(ip)
      if (x_idx_global >= myidx*imax + 1 .and. x_idx_global <= myidx*imax + imax .and. &
          y_idx_global >= myidy*jmax + 1 .and. y_idx_global <= myidy*jmax + jmax) then
        nlocal = nlocal + 1
        i_local(nlocal) = x_idx_global - myidx*imax + 1
        j_local(nlocal) = y_idx_global - myidy*jmax + 1
        global_index(nlocal) = ip
        locx_point = real(x0 + dx * (real(x_idx_global, kind=field_r) - 0.5_field_r), kind=field_r)
        locy_point = real(y0 + dy * (real(y_idx_global, kind=field_r) - 0.5_field_r), kind=field_r)
        locx(nlocal) = locx_point
        locy(nlocal) = locy_point
      end if
    end do

    ! we want to write each measurement into a contiguous slice in a collective array in the output file.
    ! to do this, we need to know the offset of each rank's local points in the global array.    
    ! we create an array of length nprocs. then in index myid we store the number of local points.
    ! myid starts at 0, so we add 1 for proper indicing
    allocate(nlocal_send(nprocs), nlocal_all_arr(nprocs))
    nlocal_send = 0
    nlocal_send(myid + 1) = nlocal

    ! then we do an allreduce sum to get the total number of points on each rank.
    call D_MPI_ALLREDUCE(nlocal_send, nlocal_all_arr, nprocs, MPI_SUM, comm3d, mpierr)
    ! we sum the array up to myid-1 to get the offset for this rank UP UNTIL THIS ID.
    ! that way, we have counted all measurement points before our own in the global array
    rank_offset = sum(nlocal_all_arr(1:myid))
    deallocate(nlocal_send, nlocal_all_arr)

    ! we construct the contiguous slice of points in the global array
    global_index = [(rank_offset + il, il = 1, nlocal)]

    ! set all global x,y coordinates to 0
    locx_all = 0.0_field_r
    locy_all = 0.0_field_r

    ! hopefully, at this point all ranks are properly separated
    do il = 1, nlocal
      locx_all(rank_offset + il) = locx(il)
      locy_all(rank_offset + il) = locy(il)
    end do
    ! when we sum, and we're properly separated, each PE will be left with the correct coordinate array
    call D_MPI_ALLREDUCE(locx_all, npoints, MPI_SUM, comm3d, mpierr)
    call D_MPI_ALLREDUCE(locy_all, npoints, MPI_SUM, comm3d, mpierr)

    ofile = multi_timeseries_file_t('virtualmeasurement', npoints, lgpu=.false., &
                                      locx=locx_all, locy=locy_all, point_ids=global_index)
    call add_output_file(ofile, dtav, ofile_id)

    call ofile%add_var('tskin', 'skin temperature', 'K', 'it','tt0t')
    call ofile%add_var('thlskin', 'skin liquid water potential temperature', 'K', 'it','tt0t')
    call ofile%add_var('qskin', 'skin specific humidity', 'kg/kg', 'it','tt0t')
    call ofile%add_var('swd', 'surface shortwave downward radiation', 'W/m^2', 'it','tt0t')
    call ofile%add_var('swu', 'surface shortwave upward radiation', 'W/m^2', 'it','tt0t')
    call ofile%add_var('lwd', 'surface longwave downward radiation', 'W/m^2', 'it','tt0t')
    call ofile%add_var('lwu', 'surface longwave upward radiation', 'W/m^2', 'it','tt0t')
    call ofile%add_var('swnet', 'surface net shortwave radiation', 'W/m^2', 'it','tt0t')
    call ofile%add_var('lwnet', 'surface net longwave radiation', 'W/m^2', 'it','tt0t')
    call ofile%add_var('t1', 'air temperature at first atmospheric level', 'K', 'it','tt0t')
    call ofile%add_var('q1', 'water vapor specific humidity at first atmospheric level', 'kg/kg', 'it','tt0t')
    call ofile%add_var('ql1', 'liquid water specific humidity at first atmospheric level', 'kg/kg', 'it','tt0t')
    call ofile%add_var('Qnet', 'surface net radiation', 'W/m^2', 'it','tt0t')
    call ofile%add_var('H', 'sensible heat flux', 'W/m^2', 'it','tt0t')
    call ofile%add_var('LE', 'latent heat flux', 'W/m^2', 'it','tt0t')
    call ofile%add_var('G0', 'ground heat flux', 'W/m^2', 'it','tt0t')
    call ofile%add_var('ustar', 'friction velocity', 'm/s', 'it','tt0t')
    call ofile%add_var('obuk', 'Obukhov length', 'm', 'it','tt0t')
    call ofile%add_var('ra', 'aerodynamic resistance', 's/m', 'it','tt0t')
    call ofile%add_var('wind1', 'wind speed at first atmospheric level', 'm/s', 'it','tt0t')

    select case (isurf)
    case (11)
      call ofile%add_var('wind10m', 'tile-averaged MOST-derived 10m wind speed', 'm/s', 'it','tt0t')
      call ofile%add_var('t2m', 'tile-averaged MOST-derived 2m temperature', 'K', 'it','tt0t')
      call ofile%add_var('q2m', 'tile-averaged MOST-derived 2m water vapor specific humidity', 'kg/kg', 'it','tt0t')
      call ofile%add_var('ql2m', 'tile-averaged MOST-derived 2m liquid water specific humidity', 'kg/kg', 'it','tt0t')
      call ofile%add_var('rsveg', 'vegetation resistance', 's/m', 'it','tt0t')
      call ofile%add_var('rssoil', 'soil evaporation resistance', 's/m', 'it','tt0t')
      call ofile%add_var('cliq', 'fraction of vegetated surface covered with liquid water', '-', 'it','tt0t')
      call ofile%add_var('Wl', 'liquid water reservoir', 'm', 'it','tt0t')
      call ofile%add_var('f1', 'f1(SWD) function vegetation resistance', '-', 'it','tt0t')
      call ofile%add_var('f2_b', 'f2(theta) function soil resistance', '-', 'it','tt0t')
    case default
      call ofile%add_var('wind10m', 'MOST-derived 10m wind speed', 'm/s', 'it','tt0t')
      call ofile%add_var('t2m', 'MOST-derived 2m temperature', 'K', 'it','tt0t')
      call ofile%add_var('q2m', 'MOST-derived 2m water vapor specific humidity', 'kg/kg', 'it','tt0t')
      call ofile%add_var('ql2m', 'MOST-derived 2m liquid water specific humidity', 'kg/kg', 'it','tt0t')
      call ofile%add_var('rs', 'surface resistance', 's/m', 'it','tt0t')
    end select

    deallocate(locx, locy, locx_all, locy_all)

  end subroutine initvirtualmeasurement

  pure function most_wind_speed(ustar_in, obuk_in, z0m_in, z_out, fallback_value) result(wind_out)
    use modsurface, only: psim

    implicit none

    real(field_r), intent(in) :: ustar_in, obuk_in, z0m_in, z_out, fallback_value

    real(field_r) :: log_term, wind_out

    if (ustar_in > 1.0e-8_field_r .and. abs(obuk_in) > 1.0e-8_field_r .and. z0m_in > 0.0_field_r .and. z_out > 1.01 * z0m_in) then
      log_term = log(z_out / z0m_in) - psim(z_out / obuk_in) + psim(z0m_in / obuk_in)
      wind_out = ustar_in / fkar * log_term
    else
      wind_out = fallback_value
    end if
  end function most_wind_speed

  pure function most_scalar_value(surface_value, scalar_flux, ustar_in, obuk_in, z0h_in, z_out, fallback_value) result(scalar_out)
    use modsurface, only: psih

    implicit none

    real(field_r), intent(in) :: surface_value, scalar_flux, ustar_in, obuk_in, z0h_in, z_out, fallback_value

    real(field_r) :: log_term, scalar_out

    if (ustar_in > 1.0e-8_field_r .and. abs(obuk_in) > 1.0e-8_field_r .and. z0h_in > 0.0_field_r .and. z_out > 1.01 * z0h_in) then
      log_term = log(z_out / z0h_in) - psih(z_out / obuk_in) + psih(z0h_in / obuk_in)
      scalar_out = surface_value - scalar_flux * log_term / (fkar * ustar_in)
    else
      scalar_out = fallback_value
    end if
  end function most_scalar_value

  subroutine virtualmeasurement
    use modfields,         only: u0, v0, thl0, qt0, ql0, exnf, exnh, presf, rhof
    use modslurb,          only: slurb_tile
    use modsurfdata,       only: isurf, tskin, qskin, Qnet, G0, H, LE, ustar, obl, thlflux, qtflux, z0h, z0m, ra, rs, &
                                 rsveg, rssoil, cliq, Wl
    use modraddata,        only: swd, swu, lwd, lwu
    use modthermodynamics, only: calc_qsat

    implicit none

    integer :: i, j, point_index, ilu
    real(field_r) :: current_obuk, current_qflux, current_qskin, current_tflux, current_tskin
    real(field_r) :: current_ustar, current_z0h, current_z0m
    real(field_r) :: qsat_2m, qt1, qt2m_total, upcu, urban_wqt, urban_wthl, vpcv, tile_wind10m, tile_weight_sum

    real(field_r), pointer :: p_tskin, p_thlskin, p_qskin, p_qnet
    real(field_r), pointer :: p_swd, p_swu, p_lwd, p_lwu, p_swnet, p_lwnet
    real(field_r), pointer :: p_wind1, p_wind10m, p_t1, p_q1, p_ql1
    real(field_r), pointer :: p_t2m, p_q2m, p_ql2m
    real(field_r), pointer :: p_ustar, p_obuk, p_ra, p_rs, p_h, p_le, p_g0
    real(field_r), pointer :: p_rsveg, p_rssoil, p_cliq, p_wl, p_f1, p_f2b

    

    if (.not. lvirtualmeasurement) return
    if (.not. is_sampling_timestep(ofile_id)) return

    do point_index = 1, size(i_local)

      i = i_local(point_index)
      j = j_local(point_index)

      call ofile%get_pointer('tskin', global_index(point_index), p_tskin)
      call ofile%get_pointer('thlskin', global_index(point_index), p_thlskin)
      call ofile%get_pointer('qskin', global_index(point_index), p_qskin)
      call ofile%get_pointer('swd', global_index(point_index), p_swd)
      call ofile%get_pointer('swu', global_index(point_index), p_swu)
      call ofile%get_pointer('lwd', global_index(point_index), p_lwd)
      call ofile%get_pointer('lwu', global_index(point_index), p_lwu)
      call ofile%get_pointer('swnet', global_index(point_index), p_swnet)
      call ofile%get_pointer('lwnet', global_index(point_index), p_lwnet)
      call ofile%get_pointer('t1', global_index(point_index), p_t1)
      call ofile%get_pointer('q1', global_index(point_index), p_q1)
      call ofile%get_pointer('ql1', global_index(point_index), p_ql1)
      call ofile%get_pointer('Qnet', global_index(point_index), p_qnet)
      call ofile%get_pointer('ustar', global_index(point_index), p_ustar)
      call ofile%get_pointer('obuk', global_index(point_index), p_obuk)
      call ofile%get_pointer('ra', global_index(point_index), p_ra)
      call ofile%get_pointer('H', global_index(point_index), p_h)
      call ofile%get_pointer('LE', global_index(point_index), p_le)
      call ofile%get_pointer('G0', global_index(point_index), p_g0)
      call ofile%get_pointer('wind1', global_index(point_index), p_wind1)
      call ofile%get_pointer('wind10m', global_index(point_index), p_wind10m)
      call ofile%get_pointer('t2m', global_index(point_index), p_t2m)
      call ofile%get_pointer('q2m', global_index(point_index), p_q2m)
      call ofile%get_pointer('ql2m', global_index(point_index), p_ql2m)

      p_tskin = tskin(i,j) * exnh(1)
      p_thlskin = tskin(i,j)
      p_qskin = qskin(i,j)
      p_swd = swd(i,j,1)
      p_swu = swu(i,j,1)
      p_lwd = lwd(i,j,1)
      p_lwu = lwu(i,j,1)
      p_swnet = swd(i,j,1) + swu(i,j,1)
      p_lwnet = lwd(i,j,1) + lwu(i,j,1)
      p_qnet = Qnet(i,j)
      p_t1 = (thl0(i,j,1) + rlv*ql0(i,j,1)/(cp*exnf(1))) * exnf(1)
      qt1 = qt0(i,j,1)
      p_q1 = qt0(i,j,1) - ql0(i,j,1)
      p_ql1 = ql0(i,j,1)
      p_ustar = ustar(i,j)
      p_obuk = obl(i,j)
      p_ra = ra(i,j)
      p_h = H(i,j)
      p_le = LE(i,j)
      p_g0 = G0(i,j)

      upcu = 0.5_field_r * (u0(i,j,1) + u0(i+1,j,1)) + cu
      vpcv = 0.5_field_r * (v0(i,j,1) + v0(i,j+1,1)) + cv
      p_wind1 = max(sqrt(upcu**2 + vpcv**2), 0.1_field_r)

      select case (isurf)
      case (11)
        call ofile%get_pointer('rsveg', global_index(point_index), p_rsveg)
        call ofile%get_pointer('rssoil', global_index(point_index), p_rssoil)
        call ofile%get_pointer('cliq', global_index(point_index), p_cliq)
        call ofile%get_pointer('Wl', global_index(point_index), p_wl)
        call ofile%get_pointer('f1', global_index(point_index), p_f1)
        call ofile%get_pointer('f2_b', global_index(point_index), p_f2b)

        tile_wind10m = 0.0_field_r
        p_t2m = 0.0_field_r
        qt2m_total = 0.0_field_r
        tile_weight_sum = 0.0_field_r
        ! we loop over LSM tiles, as each tile has a different z0(m/h), so it's hard to define which z0 to use. could also possible average z0 in future...
        do ilu = 1, nlu
          if (tile(ilu)%frac(i,j) <= 0.0_field_r) cycle

          if (tile(ilu)%lushort == 'slb') then
            urban_wthl = slurb_tile%shf_urb(i,j) / (rhof(1) * cp)
            urban_wqt = slurb_tile%qsws_urb(i,j) / (rhof(1) * rlv)

            current_ustar = slurb_tile%us_urb(i,j)
            current_obuk = slurb_tile%ol_urb(i,j)
            current_z0m = slurb_tile%z0_urb(i,j)
            current_z0h = slurb_tile%z0_urb(i,j)
            current_tskin = slurb_tile%thlskin(i,j) * exnh(1)
            current_qskin = slurb_tile%qtskin(i,j)
            current_tflux = urban_wthl * exnh(1)
            current_qflux = urban_wqt
          else
            current_ustar = tile(ilu)%ustar(i,j)
            current_obuk = tile(ilu)%obuk(i,j)
            current_z0m = tile(ilu)%z0m(i,j)
            current_z0h = tile(ilu)%z0h(i,j)
            current_tskin = tile(ilu)%thlskin(i,j) * exnh(1)
            current_qskin = tile(ilu)%qtskin(i,j)
            current_tflux = tile(ilu)%wthl(i,j) * exnh(1)
            current_qflux = tile(ilu)%wqt(i,j)
          end if

          tile_wind10m = tile_wind10m + tile(ilu)%frac(i,j) * most_wind_speed(current_ustar, current_obuk, current_z0m, 10.0_field_r, -999.0_field_r)
          ! we do MOST on the real temperature, so we need to convert back to potential temperature for the output as we don't have exn(2m)
          p_t2m = p_t2m + tile(ilu)%frac(i,j) * most_scalar_value(current_tskin, current_tflux, current_ustar, current_obuk, current_z0h, 2.0_field_r, -999.0_field_r)
          qt2m_total = qt2m_total + tile(ilu)%frac(i,j) * most_scalar_value(current_qskin, current_qflux, current_ustar, current_obuk, current_z0h, 2.0_field_r, -999.0_field_r)

          tile_weight_sum = tile_weight_sum + tile(ilu)%frac(i,j)
        end do

        if (tile_weight_sum > 0.0_field_r) then
          p_wind10m = tile_wind10m / tile_weight_sum
          p_t2m = p_t2m / tile_weight_sum
          qt2m_total = qt2m_total / tile_weight_sum
        else
          p_wind10m = -999.0_field_r
          p_t2m = -999.0_field_r
          qt2m_total = -999.0_field_r
        end if

        p_rsveg = rsveg(i,j)
        p_rssoil = rssoil(i,j)
        p_cliq = cliq(i,j)
        p_wl = Wl(i,j)
        p_f1 = f1(i,j)
        p_f2b = f2b(i,j)
      case default
        call ofile%get_pointer('rs', global_index(point_index), p_rs)

        p_rs = rs(i,j)

        p_wind10m = most_wind_speed(ustar(i,j), obl(i,j), z0m(i,j), 10.0_field_r, -999.0_field_r)
        ! we do MOST on the real temperature, so we need to convert back to potential temperature for the output as we don't have exn(2m)
        p_t2m = most_scalar_value(tskin(i,j) * exnh(1), thlflux(i,j) * exnh(1), ustar(i,j), obl(i,j), z0h(i,j), 2.0_field_r, -999.0_field_r)
        qt2m_total = most_scalar_value(qskin(i,j), qtflux(i,j), ustar(i,j), obl(i,j), z0h(i,j), 2.0_field_r, -999.0_field_r)
      end select

      qsat_2m = calc_qsat(p_t2m, presf(1))
      p_ql2m = max(qt2m_total - qsat_2m, 0.0_field_r)
      p_q2m = max(qt2m_total - p_ql2m, 0.0_field_r)
    end do

  end subroutine virtualmeasurement

  subroutine exitvirtualmeasurement
    if (allocated(i_local)) deallocate(i_local)
    if (allocated(j_local)) deallocate(j_local)
    if (allocated(global_index)) deallocate(global_index)
  end subroutine exitvirtualmeasurement

end module modvirtualmeasurement
