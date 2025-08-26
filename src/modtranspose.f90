!> Transposition routines for Poisson solver.
module modtranspose

  use modglobal,    only: itot, jtot, imax, jmax, kmax, i1, j1
  use modmpi,       only: D_MPI_ALLTOALL, commrow, commcol, nprocs, nprocx, &
                          nprocy
  use modprecision, only: pois_r, longint

  implicit none

  private

  public :: transpose_get_buffer_size
  public :: init_transpose
  public :: transpose_z_to_x
  public :: transpose_x_to_z
  public :: transpose_x_to_y
  public :: transpose_y_to_x
  public :: transpose_y_to_z
  public :: transpose_z_to_y

  public :: iony
  public :: jonx
  public :: konx

  integer, protected :: &
    iony,    & !< Number of cells in the x-direction on y-contiguous pencils.
    jonx,    & !< Number of cells in the y-direction on x-contiguous pencils.
    konx       !< Number of cells in the z-direction on x/y-contiguous pencils.

contains

  !> Compute the minimum size of the workspace for transposing.
  !!
  !! @return minimum buffer size for transposing.
  function transpose_get_buffer_size() result(size)
    integer(longint) :: size_x, size_y, size_z, size

    ! x-contiguous pencils
    size_x = itot * jmax * konx

    ! y-contiguous pencils
    size_y = iony * jtot * konx

    ! z-contiguous pencils
    size_z = imax * jmax * kmax

    size = max(size_x, size_y, size_z)

  end function transpose_get_buffer_size

  !> Initialize transposes
  !!
  !! @param[out] size Required size of the work space for transposes.
  subroutine init_transpose(size)

    integer(longint), intent(out) :: size

    konx = kmax / nprocx
    if (mod(kmax, nprocx) > 0) konx = konx + 1

    iony = itot / nprocy
    if (mod(itot, nprocy) > 0) iony = iony + 1

    jonx = jtot / nprocx
    if (mod(jtot, nprocx) > 0) jonx = jonx + 1

    size = transpose_get_buffer_size()

  end subroutine init_transpose

  ! CJ: In cuDecomp, the "vertical" direction is labeled as Y.
  ! Hence, z_to_x is actually y_to_x in cuDecomp.

  !> Tranpose z-contiguous pencils to x-contiguous pencils.
  !!
  !! @param[in] pz Input data.
  !! @param[out] px Output data.
  !! @param[out] buffer Buffer for transposing.
  subroutine transpose_z_to_x(pz, px, buffer)

    real(pois_r), intent(in)  :: pz(:,:,:)
    real(pois_r), intent(out) :: px(:,:,:)
    real(pois_r), intent(out) :: buffer(:)

    integer :: i, j, k, n, ii
    integer :: mpierr

    if (nprocs == 1) then
      !$acc parallel loop collapse(3) default(present)
      do k=1,kmax
        do j=1,jtot
          do i=1,itot
            px(i,j,k) = pz(i+1,j+1,k)
          end do
        end do
      end do
    else
      !$acc parallel loop collapse(4) default(present) private(ii)
      do n = 0, nprocx-1
        do k = 1, konx
          do j = 2, j1
            do i = 2, i1
              ii = (i-1) + (j-2)*imax + (k-1)*imax*jmax + n*imax*jmax*konx
              if (k+n*konx <= kmax) buffer(ii) = pz(i,j,k+n*konx) 
            end do
          end do
        end do
      end do

      call D_MPI_ALLTOALL(buffer, imax*jmax*konx, commrow, mpierr, lacc=.true.)

      !$acc parallel loop collapse(4) default(present) private(ii)
      do n = 0, nprocx-1
        do k = 1, konx
          do j = 1, jmax
            do i = 1, imax
              ii = i + (j-1)*imax + (k-1)*imax*jmax + n*imax*jmax*konx
              px(i+n*imax,j,k) = buffer(ii)
            end do
          end do
        end do
      end do
    end if

  end subroutine transpose_z_to_x

  !> Tranpose x-contiguous pencils to z-contiguous pencils.
  !!
  !! @param[in] px Input data.
  !! @param[out] pz Output data.
  !! @param[out] buffer Buffer for transposing.
  subroutine transpose_x_to_z(px, pz, buffer)

    real(pois_r), intent(in)  :: px(:,:,:)
    real(pois_r), intent(out) :: pz(:,:,:)
    real(pois_r), intent(out) :: buffer(:)

    integer :: i, j, k, n, ii
    integer :: mpierr

    if (nprocs == 1) then
      !$acc parallel loop collapse(3) default(present)
      do k = 1, kmax
        do j = 1, jtot
          do i = 1, itot
            pz(i+1,j+1,k) = px(i,j,k)
          end do
        end do
      end do
    else
      !$acc parallel loop collapse(4) default(present) private(ii)
      do n = 0, nprocx-1
        do k = 1, konx
          do j = 1, jmax
            do i = 1, imax
              ii = i + (j-1)*imax + (k-1)*imax*jmax + n*imax*jmax*konx
              buffer(ii) = px(i+n*imax,j,k)
            end do
          end do
        end do
      end do

      call D_MPI_ALLTOALL(buffer, imax*jmax*konx, commrow, mpierr, lacc=.true.)

      !$acc parallel loop collapse(4) default(present) private(ii)
      do n = 0, nprocx-1
        do k = 1, konx
          do j = 2, j1
            do i = 2, i1
              ii = (i-1) + (j-2)*imax + (k-1)*imax*jmax + n*imax*jmax*konx
              if (k+n*konx <= kmax) pz(i,j,k+n*konx) = buffer(ii)
            end do
          end do
        end do
      end do
    end if

  end subroutine transpose_x_to_z

  !> Tranpose x-contiguous pencils to y-contiguous pencils.
  !!
  !! @param[in] px Input data.
  !! @param[out] py Output data.
  !! @param[out] buffer Buffer for transposing.
  subroutine transpose_x_to_y(px, py, buffer)

    real(pois_r), intent(in)  :: px(:,:,:)
    real(pois_r), intent(out) :: py(:,:,:)
    real(pois_r), intent(out) :: buffer(:)

    integer :: i, j, k, n, ii
    integer :: mpierr

    if (nprocs == 1) then
      !$acc parallel loop collapse(3) default(present) private(ii)
      do k = 1, kmax
        do j = 1, jtot
          do i = 1, itot
            ii = i + (j-1)*itot + (k-1)*itot*jtot
            buffer(ii) = px(i,j,k)
          end do
        end do
      end do

      !$acc parallel loop collapse(3) default(present) private(ii)
      do k = 1, kmax
        do j = 1, jtot
         do i = 1, itot
            ii = i + (j-1)*itot + (k-1)*itot*jtot
            py(j,k,i) = buffer(ii)
          end do
        end do
      end do
    else
      !$acc parallel loop collapse(4) default(present) private(ii)
      do n = 0, nprocy-1
        do k = 1, konx
          do j = 1, jmax
            do i = 1, iony
              ii = i + (j-1)*iony + (k-1)*iony*jmax + n*iony*jmax*konx
              if (i <= itot) buffer(ii) = px(i+n*iony,j,k)
            end do
          end do
        end do
      end do

      call D_MPI_ALLTOALL(buffer, iony*jmax*konx, commcol, mpierr, lacc=.true.)


      !$acc parallel loop collapse(4) default(present) private(ii)
      do n = 0, nprocy-1
        do k = 1, konx
          do i = 1, iony
            do j = 1, jmax
              ii = i + (j-1)*iony + (k-1)*iony*jmax + n*iony*jmax*konx
              py(j+n*jmax,k,i) = buffer(ii)
            end do
          end do
        end do
      end do

    end if

  end subroutine transpose_x_to_y

  !> Tranpose y-contiguous pencils to x-contiguous pencils.
  !!
  !! @param[in] py Input data.
  !! @param[out] px Output data.
  !! @param[out] buffer Buffer for transposing.
  subroutine transpose_y_to_x(py, px, buffer)

    real(pois_r), intent(in)  :: py(:,:,:)
    real(pois_r), intent(out) :: px(:,:,:)
    real(pois_r), intent(out) :: buffer(:)

    integer :: i, j, k, n, ii
    integer :: mpierr

    if (nprocs == 1) then
      !$acc parallel loop collapse(3) default(present) private(ii)
      do k = 1, kmax
        do j = 1, jtot
          do i = 1, itot
            ii = j + (i-1)*jtot + (k-1)*itot*jtot
            buffer(ii) = py(j,k,i)
          end do
        end do
      end do

      !$acc parallel loop collapse(3) default(present) private(ii)
      do k = 1, kmax
        do j = 1, jtot
          do i = 1, itot
            ii = j + (i-1)*jtot + (k-1)*itot*jtot
            px(i,j,k) = buffer(ii)
          end do
        end do
      end do
    else
      !$acc parallel loop collapse(4) default(present) private(ii)
      do n = 0, nprocy-1
        do k = 1, konx
          do i = 1, iony
            do j = 1, jmax
              ii = i + (j-1)*iony + (k-1)*iony*jmax + n*iony*jmax*konx
              buffer(ii) = py(j+n*jmax,k,i)
            end do
          end do
        end do
      end do

      call D_MPI_ALLTOALL(buffer, iony*jmax*konx, commcol, mpierr, lacc=.true.)

      !$acc parallel loop collapse(4) default(present) private(ii)
      do n = 0, nprocy-1
        do k = 1, konx
          do j = 1, jmax
            do i = 1, iony
              ii = i + (j-1)*iony + (k-1)*iony*jmax + n*iony*jmax*konx
              if (i+n*iony <= itot) px(i+n*iony,j,k) = buffer(ii)
            end do
          end do
        end do
      end do
    end if

  end subroutine transpose_y_to_x

  !> Tranpose y-contiguous pencils to z-contiguous pencils.
  !!
  !! @param[in] py Input data.
  !! @param[out] pz Output data.
  !! @param[out] buffer Buffer for transposing.
  subroutine transpose_y_to_z(py, pz, buffer)

    real(pois_r), intent(in)  :: py(:,:,:)
    real(pois_r), intent(out) :: pz(:,:,:)
    real(pois_r), intent(out) :: buffer(:)

    integer :: i, j, k, n, ii
    integer :: mpierr

    if (nprocs == 1) then
      !$acc parallel loop collapse(3) default(present)
      do k = 1, kmax
        do j = 1, jtot
          do i = 1, itot
            pz(i+1,j+1,k) = py(j,k,i)
          end do
        end do
      end do
    else
      !$acc parallel loop collapse(4) default(present) private(ii)
      do n = 0, nprocx-1
        do k = 1, konx
          do i = 1, iony
            do j = 1, jonx
              ii = j + (i-1)*jonx + (k-1)*iony*jonx + n*iony*jonx*konx
              if (j+n*jonx <= jtot) buffer(ii) = py(j+n*jonx,k,i)
            end do
          end do
        end do
      end do

      call D_MPI_ALLTOALL(buffer, iony*jonx*konx, commrow, mpierr, lacc=.true.)

      !$acc parallel loop collapse(4) default(present) private(ii)
      do n = 0, nprocx-1
        do k = 1, konx
          do j = 1, jonx
            do i = 1, iony
              ii = j + (i-1)*jonx + (k-1)*iony*jonx + n*iony*jonx*konx
              if (k+n*konx <= kmax) pz(i,j,k+n*konx) = buffer(ii)
            end do
          end do
        end do
      end do
    end if

  end subroutine transpose_y_to_z

  !> Tranpose y-contiguous pencils to x-contiguous pencils.
  !!
  !! @param[in] py Input data.
  !! @param[out] px Output data.
  !! @param[out] buffer Buffer for transposing.
  subroutine transpose_z_to_y(pz, py, buffer)

    real(pois_r), intent(in)  :: pz(:,:,:)
    real(pois_r), intent(out) :: py(:,:,:)
    real(pois_r), intent(out) :: buffer(:)

    integer :: i, j, k, n, ii
    integer :: mpierr

    if (nprocs == 1) then
      !$acc parallel loop collapse(3) default(present)
      do k=1,kmax
        do j=1,jtot
          do i=1,itot
            py(j,k,i) = pz(i+1,j+1,k)
          end do
        end do
      end do
    else
      !$acc parallel loop collapse(4) default(present) private(ii)
      do n = 0, nprocx-1
        do k = 1, konx
          do j = 1, jonx
            do i = 1, iony
              ii = j + (i-1)*jonx + (k-1)*iony*jonx + n*iony*jonx*konx
              if (k+n*konx <= kmax) buffer(ii) = pz(i,j,k+n*konx)
            end do
          end do
        end do
      end do

      call D_MPI_ALLTOALL(buffer, iony*jonx*konx, commrow, mpierr, lacc=.true.)

      !$acc parallel loop collapse(4) default(present) private(ii)
      do n = 0, nprocx-1
        do k = 1, konx
          do i = 1, iony
            do j = 1, jonx
              ii = j + (i-1)*jonx + (k-1)*iony*jonx + n*iony*jonx*konx
              if (j+n*jonx <= jtot) py(j+n*jonx,k,i) = buffer(ii)
            end do
          end do
        end do
      end do
    end if

  end subroutine transpose_z_to_y

end module modtranspose