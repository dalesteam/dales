!> Subroutines for aerosol scavenging.
module modaerosol_scavenging

  use bulkmicro_sb,      only: calc_sed_qr_sb
  use modaerosol_common, only: calc_median_diameter
  use modaerosol_mode_t, only: aerosol_mode_t, hydrometeor_mode_t
  use modbulkmicro_data, only: qrmin, qcmin
  use modglobal,         only: i1, j1, kmax, ih, jh, rhow, pi
  use modmpi,            only: myid, d_mpi_bcast
  use modprecision,      only: field_r
  use modstat_nc,        only: nchandle_error

  use netcdf

  implicit none

  private

  public :: init_scavenging
  public :: aerosol_scavenging_rain_lut
  public :: aerosol_scavenging_cloud_lut

  real(field_r), allocatable :: log_rp_inc(:) !< Log of aerosol particle radii, for in-cloud LUT.
  real(field_r), allocatable :: log_rp_blc(:) !< Log of aerosol particle radii, for below-cloud LUT.

  real(field_r), allocatable :: log_rr(:) !< Log of the rain rates [log(mm hr-1)].
  real(field_r), allocatable :: log_rc(:) !< Log of the cloud droplet radii [log(10-6 m)].

  ! TODO: it makes no sense to store these luts in FP64; the original tables
  ! only provide 5 digits of accuracy! Better to default to FP32, or even FP16
  ! on GPU.
  real(field_r), allocatable :: gamma_inc_m(:,:) !< Scavenging rates for mass, in-cloud [s-1].
  real(field_r), allocatable :: gamma_inc_n(:,:) !< Scavenging rates for number, in-cloud [s-1].
  real(field_r), allocatable :: gamma_blc_m(:,:) !< Scavenging rates for mass, below-cloud [s-1].
  real(field_r), allocatable :: gamma_blc_n(:,:) !< Scavenging rates for number, below-cloud [s-1].

contains

  !> Read the scavenging lookup tables.
  subroutine init_scavenging()

    integer :: istat, ncid, varid, dimid
    integer :: dims_inc(2) !< Dimensions of the in-cloud LUT.
    integer :: dims_blc(2) !< DImensions of the below-cloud LUT.

    call nchandle_error(nf90_open("scavenging_lut.nc", NF90_NOWRITE, ncid))

    call nchandle_error(nf90_inq_dimid(ncid, "rp_blc", dimid))
    call nchandle_error(nf90_inquire_dimension(ncid, dimid, len=dims_blc(2)))

    call nchandle_error(nf90_inq_dimid(ncid, "log_rr", dimid))
    call nchandle_error(nf90_inquire_dimension(ncid, dimid, len=dims_blc(1)))

    call nchandle_error(nf90_inq_dimid(ncid, "rp_inc", dimid))
    call nchandle_error(nf90_inquire_dimension(ncid, dimid, len=dims_inc(2)))

    call nchandle_error(nf90_inq_dimid(ncid, "log_rc", dimid))
    call nchandle_error(nf90_inquire_dimension(ncid, dimid, len=dims_inc(1)))

    allocate(log_rp_inc(dims_inc(2)), &
             log_rp_blc(dims_blc(2)), &
             log_rr(dims_blc(1)), &
             log_rc(dims_inc(1)), &
             gamma_inc_m(dims_inc(1),dims_inc(2)), &
             gamma_inc_n(dims_inc(1),dims_inc(2)), &
             gamma_blc_m(dims_blc(1),dims_blc(2)), &
             gamma_blc_n(dims_blc(1),dims_blc(2)))

    call nchandle_error(nf90_inq_varid(ncid, "rp_blc", varid))
    call nchandle_error(nf90_get_var(ncid, varid, log_rp_blc))

    call nchandle_error(nf90_inq_varid(ncid, "log_rr", varid))
    call nchandle_error(nf90_get_var(ncid, varid, log_rr))

    call nchandle_error(nf90_inq_varid(ncid, "rp_inc", varid))
    call nchandle_error(nf90_get_var(ncid, varid, log_rp_inc))

    call nchandle_error(nf90_inq_varid(ncid, "log_rc", varid))
    call nchandle_error(nf90_get_var(ncid, varid, log_rc))

    call nchandle_error(nf90_inq_varid(ncid, "gamma_inc_m", varid))
    call nchandle_error(nf90_get_var(ncid, varid, gamma_inc_m))

    call nchandle_error(nf90_inq_varid(ncid, "gamma_inc_n", varid))
    call nchandle_error(nf90_get_var(ncid, varid, gamma_inc_n))

    call nchandle_error(nf90_inq_varid(ncid, "gamma_blc_m", varid))
    call nchandle_error(nf90_get_var(ncid, varid, gamma_blc_m))

    call nchandle_error(nf90_inq_varid(ncid, "gamma_blc_n", varid))
    call nchandle_error(nf90_get_var(ncid, varid, gamma_blc_n))

    call nchandle_error(nf90_close(ncid))

  end subroutine init_scavenging

  !> Compute washout of aerosols by precipitation (below-cloud scavenging).
  !!
  !! Method taken from Croft et al. (2009). Eq (1):
  !!
  !! \f[
  !!   dC/dt = Cf(RF)
  !! \f]
  !!
  !! Where \f$C\f$ is the ambient mixing ratio of tracer, \f$f\f$ the cloud
  !! fraction, \f$R\f$ the normalized scavenging coeff. and \f$F\f$ the
  !! precipitation flux.
  !!
  !! However:
  !! 1) The values in the lookuptable (L) are applied as R*F as stated on page
  !!    4656 of Croft et al.
  !! 2) In DALES a gridbox is either cover by cloud or not. So cloud fraction
  !!    is not used, i.e. equals 1 if  method is applied.
  subroutine aerosol_scavenging_rain_lut(qr, nr, rho, delt, f_mode, r_mode)

    real(field_r), intent(in) :: qr(2:,2:,:) !< Rain water content [kg kg-1].
    real(field_r), intent(in) :: nr(2:,2:,:) !< Rain number concentration [m-3].
    real(field_r), intent(in) :: rho(:)      !< Air density [kg m-3].
    real(field_r), intent(in) :: delt        !< Time step size [s].

    class(aerosol_mode_t),     intent(inout) :: f_mode !< Free mode.
    class(hydrometeor_mode_t), intent(inout) :: r_mode !< In-rain mode.

    real(field_r) :: sed_qr  !< Sedimentation rate (= rain rate?) [mm hr-1].
    real(field_r) :: rm      !< Median aerosol radius [microns].
    real(field_r) :: gamma_n !< Number scavenging rate [s-1].
    real(field_r) :: gamma_m !< Mass scavenging rate [s-1].

    integer :: i, j, k, s, st

    do k = 1, kmax
      do j = 2, j1
        do i = 2, i1
          sed_qr = calc_sed_qr_sb(qr(i,j,k), nr(i,j,k), rho(k)) * 3600
          if (qr(i,j,k) > qrmin .and. sed_qr > 0.01) then
            sed_qr = min(max(sed_qr, 0.01001_field_r), 99.999_field_r)

            rm = calc_median_diameter(f_mode%n(i,j,k), f_mode%q(i,j,k,:), &
                                      f_mode%rho, f_mode%sig_g) * 0.5 * 1E6
            rm = min(max(1.5E-3_field_r, rm), 0.9999E3_field_r)

            gamma_n = interpolate_lut(gamma_blc_n, log_rr, log_rp_blc, &
                                      [log(sed_qr), log(rm)])

            gamma_m = interpolate_lut(gamma_blc_m, log_rr, log_rp_blc, &
                                      [log(sed_qr), log(rm)])

            f_mode%np(i,j,k) = f_mode%np(i,j,k) - gamma_n * f_mode%n(i,j,k)

            do s = 1, f_mode%nspecies
              st = f_mode%to_hydro%cnct(2,s)
              f_mode%qp(i,j,k,s) = f_mode%qp(i,j,k,s) &
                                   - gamma_m * f_mode%q(i,j,k,s)
              r_mode%qp(i,j,k,s) = r_mode%qp(i,j,k,s) &
                                   + gamma_m * f_mode%q(i,j,k,s)
            end do
          end if
        end do
      end do
    end do

  end subroutine aerosol_scavenging_rain_lut

  !> Compute washout of aerosols by cloud droplets (in-cloud scavenging).
  !!
  !! Apply the interpolation to find efficiency coefficient
  !! Calculate tendency and transfer to corresponding arrays
  !!
  !! Note that:
  !! Stated in pers. comm. from Betty Croft:
  !! "Data files containing in-cloud scavenging coefficients for a cloud
  !! droplet number concentration (CDNC) of 1 cm^-3."
  !! Also:
  !! "Coefficients are given as a function of mode radius and can be used
  !! directly by multiplying by the CDNC"
  !! So we multiply by the CDNC in units of cm^-1, i.e. Nc(i,j,k)*1e-6
  subroutine aerosol_scavenging_cloud_lut(qc, nc, rho, delt, f_mode, c_mode)

    real(field_r), intent(in) :: qc(2-ih:,2-jh:,:) !< Cloud water content [kg kg-1].
    real(field_r), intent(in) :: nc(2:,2:,:)       !< Cloud number concentration [m-3].
    real(field_r), intent(in) :: rho(:)            !< Air density [kg m-3].
    real(field_r), intent(in) :: delt              !< Time step size [s].

    class(aerosol_mode_t),     intent(inout) :: f_mode !< Free aerosol mode.
    class(hydrometeor_mode_t), intent(inout) :: c_mode !< In-cloud mode.

    logical       :: limit   !< Limit the scavenging rate.
    real(field_r) :: rc      !< Mean cloud drop radius [microns].
    real(field_r) :: rm      !< Mean aerosol radius [cm].
    real(field_r) :: gamma_n !< Number scavenging rate [s-1].
    real(field_r) :: gamma_m !< Mass scavenging rate [s-1].

    integer :: i, j, k, s, st

    do k = 1, kmax
      do j = 2, j1
        do i = 2, i1
          if (qc(i,j,k) > qcmin) then
            rc = 1E6 * (3 * qc(i,j,k) * rho(k) &
                  / (4 * pi * nc(i,j,k) * rhow + 1E-16))**(1.0_field_r / 3)
            rc = min(max(rc, 5.001), 49.999)

            rm = calc_median_diameter(f_mode%n(i,j,k), f_mode%q(i,j,k,:), &
                                      f_mode%rho, f_mode%sig_g) * 0.5 * 100
            rm = min(max(rm, 1.0E-8), 8.0E-3)

            gamma_n = interpolate_lut(gamma_inc_n, log_rc, log_rp_inc, &
                                      [log(rc), log(rm)]) * nc(i,j,k) * 1E-6

            gamma_m = interpolate_lut(gamma_inc_m, log_rc, log_rp_inc, &
                                      [log(rc), log(rm)]) * nc(i,j,k) * 1E-6

            limit = gamma_n * delt > 1.0_field_r &
                    .or. gamma_m * delt > 1.0_field_r
            gamma_n = merge(1 / delt, gamma_n, limit)
            gamma_m = merge(1 / delt, gamma_m, limit)

            f_mode%np(i,j,k) = f_mode%np(i,j,k) - gamma_n * f_mode%n(i,j,k)

            do s = 1, f_mode%nspecies
              st = f_mode%to_hydro%cnct(2,s)
              f_mode%qp(i,j,k,s) = f_mode%qp(i,j,k,s) &
                                   - gamma_m * f_mode%q(i,j,k,s)
              c_mode%qp(i,j,k,s) = c_mode%qp(i,j,k,s) &
                                   + gamma_m * f_mode%q(i,j,k,s)
            end do
          end if
        end do
      end do
    end do

  end subroutine aerosol_scavenging_cloud_lut

  !> Find the index of a value in a given array, where arr(i) < val < array(i+1)
  pure function binary_search(array, value) result(idx)

    real(field_r), intent(in) :: array(:) !< Array to search.
    real(field_r), intent(in) :: value    !< Value to match.

    integer :: idx
    integer :: left, right

    left = 1
    right = size(array)
    idx = right / 2

    do while (left <= right)
      if (value >= array(idx) .and. value < array(idx + 1)) then
        exit
      else if (value < array(idx)) then
        right = idx
        idx = (right + left) / 2
      else if (value > array(idx + 1)) then
        left = idx
        idx = (right + left) / 2
      end if
    end do

  end function binary_search

  !> Interpolates a lookup table at a given location.
  pure function interpolate_lut(lut, xc, yc, loc) result(f)

    real(field_r), intent(in) :: lut(:,:) !< Lookup table.
    real(field_r), intent(in) :: xc(:)    !< Coordinates of the first dimension.
    real(field_r), intent(in) :: yc(:)    !< Coordinates of the second dimension.
    real(field_r), intent(in) :: loc(2)   !< Location [x,y] to interpolate at.

    integer       :: i, j
    real(field_r) :: x, y
    real(field_r) :: x1, x2
    real(field_r) :: y1, y2
    real(field_r) :: f

    x = loc(1)
    y = loc(2)

    i = binary_search(xc, x)
    j = binary_search(yc, y)

    x1 = xc(i)
    x2 = xc(i+1)

    y1 = yc(j)
    y2 = yc(j+1)

    f = lut(i,j) * (x2 - x) * (y2 - y) + lut(i+1,j) * (x - x1) * (y2 - y) &
        + lut(i,j+1) * (x2 - x) * (y - y1) + lut(i+1,j+1) * (x - x1) * (y - y1)

    f = f / ((x2 - x1) * (y2 - y1))

  end function interpolate_lut

end module modaerosol_scavenging