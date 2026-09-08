!> Subroutines for aerosol scavenging.
module modaerosol_scavenging

  use bulkmicro_sb,      only: calc_sed_qr_sb
  use modaerosol_common, only: calc_median_diameter, maxspecies
  use modaerosol_mode_t, only: aerosol_mode_t, hydrometeor_mode_t
  use modbulkmicro_data, only: qrmin, qcmin
  use modglobal,         only: i1, j1, kmax, ih, jh, rhow, pi
  use modmpi,            only: myid, d_mpi_bcast
  use modprecision,      only: field_r
  use modstat_nc,        only: nchandle_error
  use modtimer,          only: timer_tic, timer_toc
  use fortran_support,   only: finish

  use netcdf

  implicit none

  private

  character(len=*), parameter :: modname = 'modaerosol_scavenging'

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

    character(len=*), parameter :: routine =  modname//'/init_scavenging'

    integer :: istat, ncid, varid, dimid
    integer :: dims_inc(2) !< Dimensions of the in-cloud LUT.
    integer :: dims_blc(2) !< DImensions of the below-cloud LUT.
    logical :: file_exists

    inquire(file="scavenging_lut.nc", exist=file_exists)

    if (.not. file_exists) then
      call finish(routine, "scavenging lookup table file 'scavenging_lut.nc' &
        &not found.")
    end if

    call nchandle_error(nf90_open("scavenging_lut.nc", NF90_NOWRITE, ncid))

    call nchandle_error(nf90_inq_dimid(ncid, "log_rp_blc", dimid))
    call nchandle_error(nf90_inquire_dimension(ncid, dimid, len=dims_blc(2)))

    call nchandle_error(nf90_inq_dimid(ncid, "log_rr", dimid))
    call nchandle_error(nf90_inquire_dimension(ncid, dimid, len=dims_blc(1)))

    call nchandle_error(nf90_inq_dimid(ncid, "log_rp_inc", dimid))
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

    call nchandle_error(nf90_inq_varid(ncid, "log_rp_blc", varid))
    call nchandle_error(nf90_get_var(ncid, varid, log_rp_blc))

    call nchandle_error(nf90_inq_varid(ncid, "log_rr", varid))
    call nchandle_error(nf90_get_var(ncid, varid, log_rr))

    call nchandle_error(nf90_inq_varid(ncid, "log_rp_inc", varid))
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

    !$acc enter data copyin(log_rp_inc(1:dims_inc(2)), &
    !$acc                   log_rp_blc(1:dims_blc(2)), &
    !$acc                   log_rr(1:dims_blc(1)), &
    !$acc                   log_rc(1:dims_inc(1)), &
    !$acc                   gamma_inc_m(1:dims_inc(1),1:dims_inc(2)), &
    !$acc                   gamma_inc_n(1:dims_inc(1),1:dims_inc(2)), &
    !$acc                   gamma_blc_m(1:dims_blc(1),1:dims_blc(2)), &
    !$acc                   gamma_blc_n(1:dims_blc(1),1:dims_blc(2)))

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

    character(len=*), parameter :: routine = &
      modname//'/aerosol_scavenging_rain_lut'

    real(field_r) :: sed_qr  !< Sedimentation rate (= rain rate?) [mm hr-1].
    real(field_r) :: rm      !< Median aerosol radius [microns].
    real(field_r) :: gamma_n !< Number scavenging rate [s-1].
    real(field_r) :: gamma_m !< Mass scavenging rate [s-1].

    integer :: i, j, k, s, st

    call timer_tic(routine, 2)

    !$acc parallel loop collapse(3) default(present) &
    !$acc private(sed_qr, rm, gamma_n, gamma_m, st)
    do k = 1, kmax
      do j = 2, j1
        do i = 2, i1
          sed_qr = calc_sed_qr_sb(qr(i,j,k), nr(i,j,k), rho(k)) * 3600
          if (qr(i,j,k) > qrmin .and. sed_qr > 0.01 .and. f_mode%nspecies > 0 .and. nr(i,j,k) > 1E3) then
            sed_qr = log(sed_qr)
            sed_qr = min(max(sed_qr, -4.60517_field_r), 4.60517_field_r)

            rm = calc_median_diameter(f_mode%n(i,j,k), f_mode%q(:,i,j,k), &
                                      f_mode%rho, f_mode%sig_g) * 0.5 * 1E6
            rm = log(rm + 1E-16)
            rm = min(max(rm, -6.907755_field_r), 6.907755_field_r)

            gamma_n = interpolate_lut(gamma_blc_n, log_rr, log_rp_blc, &
                                      [sed_qr, rm])

            gamma_m = interpolate_lut(gamma_blc_m, log_rr, log_rp_blc, &
                                      [sed_qr, rm])

            f_mode%np(i,j,k) = f_mode%np(i,j,k) - gamma_n * f_mode%n(i,j,k)

            do s = 1, f_mode%nspecies
              st = f_mode%to_hydro%cnct(2,s)
              f_mode%qp(s,i,j,k) = f_mode%qp(s,i,j,k) &
                                   - gamma_m * f_mode%q(s,i,j,k)
              r_mode%qp(st,i,j,k) = r_mode%qp(st,i,j,k) &
                                   + gamma_m * f_mode%q(s,i,j,k)
            end do
          end if
        end do
      end do
    end do

    call timer_toc(routine)

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

    character(len=*), parameter :: routine = &
      modname//'/aerosol_scavenging_cloud_lut'

    logical       :: limit   !< Limit the scavenging rate.
    real(field_r) :: rc      !< Mean cloud drop radius [microns].
    real(field_r) :: rm      !< Mean aerosol radius [cm].
    real(field_r) :: gamma_n !< Number scavenging rate [s-1].
    real(field_r) :: gamma_m !< Mass scavenging rate [s-1].

    integer :: i, j, k, s, st

    call timer_tic(routine, 2)

    !$acc parallel loop collapse(3) default(present) &
    !$acc private(rc, rm, gamma_n, gamma_m, limit, st)
    do k = 1, kmax
      do j = 2, j1
        do i = 2, i1
          if (qc(i,j,k) > qcmin .and. f_mode%nspecies > 0 .and. nc(i,j,k) > 1E3) then
            rc = 1E6 * (3 * qc(i,j,k) * rho(k) &
                  / (4 * pi * nc(i,j,k) * rhow + 1E-16))**(1.0_field_r / 3)
            rc = log(rc)
            rc = min(max(rc, 1.609438_field_r), 3.912023_field_r)

            rm = calc_median_diameter(f_mode%n(i,j,k), f_mode%q(:,i,j,k), &
                                      f_mode%rho, f_mode%sig_g) * 0.5 * 100
            rm = log(rm + 1E-16)
            rm = min(max(rm, -18.42068_field_r), -4.788786_field_r)

            gamma_n = interpolate_lut(gamma_inc_n, log_rc, log_rp_inc, &
                                      [rc, rm]) * nc(i,j,k) * 1E-6

            gamma_m = interpolate_lut(gamma_inc_m, log_rc, log_rp_inc, &
                                      [rc, rm]) * nc(i,j,k) * 1E-6

            limit = gamma_n * delt > 1.0_field_r &
                    .or. gamma_m * delt > 1.0_field_r
            gamma_n = merge(1 / delt, gamma_n, limit)
            gamma_m = merge(1 / delt, gamma_m, limit)

            f_mode%np(i,j,k) = f_mode%np(i,j,k) - gamma_n * f_mode%n(i,j,k)

            do s = 1, f_mode%nspecies
              st = f_mode%to_hydro%cnct(2,s)
              f_mode%qp(s,i,j,k) = f_mode%qp(s,i,j,k) &
                                   - gamma_m * f_mode%q(s,i,j,k)
              c_mode%qp(st,i,j,k) = c_mode%qp(st,i,j,k) &
                                   + gamma_m * f_mode%q(s,i,j,k)
            end do
          end if
        end do
      end do
    end do

    call timer_toc(routine)

  end subroutine aerosol_scavenging_cloud_lut

  !> Find the index of a value in a given array, where arr(i) < val < array(i+1)
  pure function binary_search(array, value) result(idx)

    real(field_r), intent(in) :: array(:) !< Array to search.
    real(field_r), intent(in) :: value    !< Value to match.

    integer :: idx
    integer :: left, right

    left = 1
    right = size(array)
    idx = (right + left) / 2

    if (value <= array(left)) then
      idx = left
    else if (value >= array(right)) then
      idx = right
    else
      do while (left <= right)
        if (value >= array(idx) .and. value <= array(idx + 1)) then
          exit
        else if (value < array(idx)) then
          right = idx
          idx = (right + left) / 2
        else if (value > array(idx + 1)) then
          left = idx
          idx = (right + left) / 2
        end if
      end do
    end if

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

    if (i == size(xc)) then
      i = i - 1
    end if

    if (j == size(yc)) then
      j = j - 1
    end if

    x1 = xc(i)
    x2 = xc(i+1)
    y1 = yc(j)
    y2 = yc(j+1)

    f = lut(i,j) * (x2 - x) * (y2 - y) + lut(i+1,j) * (x - x1) * (y2 - y) &
        + lut(i,j+1) * (x2 - x) * (y - y1) + lut(i+1,j+1) * (x - x1) * (y - y1)

    f = f / ((x2 - x1) * (y2 - y1))

  end function interpolate_lut

end module modaerosol_scavenging
