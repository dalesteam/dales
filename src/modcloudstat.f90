!> Computes various cloud statistics.
module modcloudstat

  use fortran_support,   only: nnml_output
  use modglobal,         only: i1, j1, kmax, dzf, zf, imax, jmax, itot, jtot, &
                               ifnamopt, checknamelisterror, output_prefix
  use modfields,         only: ql0, qt0, sv0, rhobf, ql0av, exnf, thvf, w0, thl0
  use modprecision,      only: field_r
  use modthermodynamics, only: calc_virt_pot_temp
  use modtimer,          only: timer_tic, timer_toc
  use modtracers,        only: get_tracer_index
  use modnetcdf_file_t,  only: cross_section_file_t
  use modstat_nc_files,  only: add_output_file, is_sampling_timestep
  use modmpi,            only: cmyidy, cmyidx, myidx, myidy, nprocx, nprocy, &
                               comm3d, d_mpi_bcast, myid

  implicit none

  private

  character(len=*), parameter :: modname = 'modcloudstat'

  public :: cloudstat_read_namelist
  public :: init_cloudstat
  public :: do_cloudstat

  type(cross_section_file_t) :: ofile
  integer                    :: ofile_id

  logical :: lcloudstat
  real    :: dtav = 60

  integer :: iqr
  integer :: fileid

contains

  !> Read cloudstat namelist.
  subroutine cloudstat_read_namelist(nml_filename)

    character(len=*), intent(in) :: nml_filename

    integer :: ierr

    namelist /cloudstat/ lcloudstat, dtav

    if (myid == 0) then
      open(ifnamopt, file=nml_filename, status='old', iostat=ierr)
      read(ifnamopt, cloudstat, iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'cloudstat')
      write(nnml_output, cloudstat)
      close(ifnamopt)
    end if

    call d_mpi_bcast(lcloudstat, 1, 0, comm3d, ierr)
    call d_mpi_bcast(dtav, 1, 0, comm3d, ierr)

  end subroutine cloudstat_read_namelist

  subroutine init_cloudstat

    if (lcloudstat) then
      ! Make a new NetCDF file
      ofile = cross_section_file_t(trim(output_prefix)//'cloudstat.nc', nx=itot, ny=jtot)

      ! Add the file to the list of output files
      call add_output_file(ofile, dtav, ofile_id)

      call ofile%add_var('lwp', 'liquid water path', 'kg/m2', 'tt0t')
      call ofile%add_var('twp', 'total water path', 'kg/m2', 'tt0t')

      ! Check if we have rain
      iqr = get_tracer_index('qr')
      if(iqr > 0) call ofile%add_var('rwp', 'rain water path', 'kg/m2', 'tt0t')

      call ofile%add_var('thlcb', 'thl at cloudbase', 'K', 'tt0t')
      call ofile%add_var('buoycb', 'buoyancy at cloudbase', 'K', 'tt0t')
      call ofile%add_var('qtcb', 'qt at cloudbase', 'kg/kg', 'tt0t')
      call ofile%add_var('qlcb', 'ql at cloudbase', 'kg/kg', 'tt0t')
      call ofile%add_var('wcb', 'w at cloudbase', 'm/s', 'tt0t')
      call ofile%add_var('hw2cb', '1/2 W^2 at the top of the subcloud layer', &
                         'm^2/s^2', 'tt0t')
      call ofile%add_var('cldtop', 'cloud top height', 'm', 'tt0t')
      call ofile%add_var('buoymax', 'maximum buoyancy', 'K', 'tt0t')
      call ofile%add_var('hw2max', 'highest 1/2 W^2 at the top of the subcloud layer', &
                         'm^2/s^2', 'tt0t')
    end if

  end subroutine init_cloudstat

  subroutine do_cloudstat

    character(len=*), parameter :: routine = modname//'/do_cloudstat'

    integer :: i, j, k

    integer :: kcb !< Index of cloud base height.

    real(field_r) :: ql_at_cb
    real(field_r) :: qt_at_cb
    real(field_r) :: thl_at_cb
    real(field_r) :: thv_at_cb
    real(field_r) :: w_at_cb
    real(field_r) :: thv

    real(field_r), pointer :: lwp(:,:)    !< Liquid water path [kg/m2]
    real(field_r), pointer :: twp(:,:)    !< Total water path [kg/m2]
    real(field_r), pointer :: rwp(:,:)    !< Rain water path [kg/m2]
    real(field_r), pointer :: thlcb(:,:)  !< Liquid potential temperature at cloud base [K]
    real(field_r), pointer :: buoycb(:,:) !< Buoyancy at cloud base [K]
    real(field_r), pointer :: qtcb(:,:)   !< Total water specific humidity at cloud base [kg/kg]
    real(field_r), pointer :: qlcb(:,:)   !< Liquid water specific humditiy at cloud base [kg/kg]
    real(field_r), pointer :: wcb(:,:)    !< Vertical velocity at cloud base [m/s]
    real(field_r), pointer :: hw2cb(:,:)  !< Highest 1/2 w^2 at cloud base [m2/s2]
    real(field_r), pointer :: cldtop(:,:)
    real(field_r), pointer :: buoymax(:,:)
    real(field_r), pointer :: hw2max(:,:)

    if (lcloudstat) then
      if (is_sampling_timestep(ofile_id)) then

        call timer_tic(routine)

        ! ------------------------------------------------------------------------
        ! Liquid, total and rain water paths
        ! ------------------------------------------------------------------------

        call ofile%get_pointer('lwp', lwp)
        call ofile%get_pointer('twp', twp)

        !$acc parallel loop collapse(3) default(present)
        do k = 1, kmax
          do j = 2, j1
            do i = 2, i1
              !$acc atomic update
              lwp(i,j) = lwp(i,j) + rhobf(k) * ql0(i,j,k) * dzf(k)
              !$acc atomic update
              twp(i,j) = twp(i,j) + rhobf(k) * qt0(i,j,k) * dzf(k)
            end do
          end do
        end do

        if (iqr > 0) then
          call ofile%get_pointer('rwp', rwp)

          !$acc parallel loop collapse(3) default(present)
          do k = 1, kmax
            do j = 2, j1
              do i = 2, i1
                !$acc atomic update
                rwp(i,j) = rwp(i,j) + rhobf(k) * sv0(i,j,k,iqr) * dzf(k)
              end do
            end do
          end do

        end if

        ! ------------------------------------------------------------------------
        ! Compute cloud base height (highest level below which it is non-cloudy)
        ! ------------------------------------------------------------------------

        !$acc serial default(present) copyout(kcb)
        do k = 2, kmax
          if (ql0av(k-1) < 0.001) then
            kcb = k
            exit
          end if
        end do
        !$acc end serial

        ! ------------------------------------------------------------------------
        ! Variables at cloud-base
        ! ------------------------------------------------------------------------

        call ofile%get_pointer('thlcb', thlcb)
        call ofile%get_pointer('buoycb', buoycb)
        call ofile%get_pointer('qtcb', qtcb)
        call ofile%get_pointer('qlcb', qlcb)

        !$acc parallel loop collapse(2) default(present) &
        !$acc private(thl_at_cb, ql_at_cb, thv_at_cb)
        do j = 2, j1
          do i = 2, i1
            thl_at_cb = thl0(i,j,kcb)
            qt_at_cb = qt0(i,j,kcb)
            ql_at_cb = ql0(i,j,kcb)
            thv_at_cb = calc_virt_pot_temp(thl_at_cb, qt_at_cb, ql_at_cb, &
                                           exnf(kcb))
            thlcb(i,j) = thl_at_cb
            buoycb(i,j) = thv_at_cb - thvf(kcb)
            qtcb(i,j) = qt_at_cb
            qlcb(i,j) = ql_at_cb
          end do
        end do

        call ofile%get_pointer('wcb', wcb)
        call ofile%get_pointer('hw2cb', hw2cb)

        !$acc parallel loop collapse(2) default(present) &
        !$acc private(w_at_cb)
        do j = 2, j1
          do i = 2, i1
            w_at_cb = (w0(i,j,kcb) + w0(i,j,kcb+1)) / 2
            wcb(i,j) = w_at_cb
            hw2cb(i,j) = 0.5_field_r * w_at_cb * abs(w_at_cb)
          end do
        end do

        ! ------------------------------------------------------------------------
        ! Max values
        ! ------------------------------------------------------------------------

        call ofile%get_pointer('cldtop', cldtop)
        call ofile%get_pointer('buoymax', buoymax)

        !$acc parallel loop seq default(present)
        do k = 1, kmax
          !$acc loop collapse(2)
          do j = 2, j1
            do i = 2, i1
              thv = calc_virt_pot_temp(thl0(i,j,k), qt0(i,j,k), ql0(i,j,k), &
                                       exnf(k))
              if (ql0(i,j,k) > 1.0E-10_field_r) cldtop(i,j) = zf(k)
              if (thv - thvf(k) > buoymax(i,j)) buoymax(i,j) = thv - thvf(k)
            end do
          end do
        end do

        call ofile%get_pointer('hw2max', hw2max)

        !$acc parallel loop seq default(present)
        do k = 1, kmax
          !$acc loop collapse(2)
          do j = 2, j1
            do i = 2, i1
              if (w0(i,j,k)**2 > hw2max(i,j)) then
                hw2max(i,j) = 0.5 * w0(i,j,k) * abs(w0(i,j,k))
              end if
            end do
          end do
        end do

        !$acc wait

        call timer_toc(routine)

      end if
    end if

  end subroutine do_cloudstat

end module modcloudstat
