!> Computes various cloud statistics.
module modcloudstat

  use modglobal,         only: i1, j1, kmax, dzf, zf
  use modfields,         only: ql0, qt0, sv0, rhobf, ql0av, exnf, thvf, w0, thl0
  use modstat_2d,        only: add_slice, get_slice, do_stats
  use modprecision,      only: field_r
  use modthermodynamics, only: calc_virt_pot_temp
  use modtimer,          only: timer_tic, timer_toc
  use modtracers,        only: get_tracer_index

  implicit none

  private

  character(len=*), parameter :: modname = 'modcloudstat'

  public :: init_cloudstat
  public :: do_cloudstat

  integer :: iqr

contains

  subroutine init_cloudstat

    call add_slice('lwp', 'liquid water path', 'kg/m2', 'tt0t')
    call add_slice('twp', 'total water path', 'kg/m2', 'tt0t')

    ! Check if we have rain
    iqr = get_tracer_index('qr')
    if(iqr > 0) call add_slice('rwp', 'rain water path', 'kg/m2', 'tt0t')

    call add_slice('thlcb', 'thl at cloudbase', 'K', 'tt0t')
    call add_slice('buoycb', 'buoyancy at cloudbase', 'K', 'tt0t')
    call add_slice('qtcb', 'qt at cloudbase', 'kg/kg', 'tt0t')
    call add_slice('qlcb', 'ql at cloudbase', 'kg/kg', 'tt0t')
    call add_slice('wcb', 'w at cloudbase', 'm/s', 'tt0t')
    call add_slice('hw2cb', '1/2 W^2 at the top of the subcloud layer', &
                   'm^2/s^2', 'tt0t')
    call add_slice('cldtop', 'cloud top height', 'm', 'tt0t')
    call add_slice('buoymax', 'maximum buoyancy', 'K', 'tt0t')
    call add_slice('hw2max', 'highest 1/2 W^2 at the top of the subcloud layer', &
                   'm^2/s^2', 'tt0t')

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

    if (do_stats) then

      call timer_tic(routine)

      ! ------------------------------------------------------------------------
      ! Liquid, total and rain water paths
      ! ------------------------------------------------------------------------

      lwp => get_slice('lwp')
      twp => get_slice('twp')

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
        rwp => get_slice('rwp')

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

      ! ------------------------------------------------------------------------
      ! Variables at cloud-base
      ! ------------------------------------------------------------------------

      thlcb => get_slice('thlcb')
      buoycb => get_slice('buoycb')
      qtcb => get_slice('qtcb')
      qlcb => get_slice('qlcb')

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

      wcb => get_slice('wcb')
      hw2cb => get_slice('hw2cb')

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

      cldtop => get_slice('cldtop')
      buoymax => get_slice('buoymax')

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

      hw2max => get_slice('hw2max')

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

  end subroutine do_cloudstat

end module modcloudstat