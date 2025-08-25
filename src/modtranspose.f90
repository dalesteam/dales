!> Transposition routines for Poisson solver.
module modtranspose

  use modglobal,    only: itot, jtot, imax, jmax, kmax, i1, j1
  use modmpi,       only: D_MPI_ALLTOALL, commrow, commcol, nprocs, nprocx, nprocy
  use modprecision, only: pois_r

  implicit none

  private

  public :: init_transpose
  public :: transpose_z_to_x
  public :: transpose_x_to_z
  public :: transpose_x_to_y
  public :: transpose_y_to_x
  public :: transpose_y_to_z
  public :: transpose_z_to_y

  integer :: iony, jonx, konx
  integer :: mpierr

contains

  subroutine init_transpose(iony_, jonx_, konx_)
    integer, intent(in) :: iony_, jonx_, konx_

    iony = iony_
    jonx = jonx_
    konx = konx_

  end subroutine init_transpose

  subroutine transpose_z_to_x(p, px, buffer)

    real(pois_r), intent(in)  :: p(:,:,:)
    real(pois_r), intent(out) :: px(:,:,:)
    real(pois_r), intent(out) :: buffer(:)

    integer :: i, j, k, n, ii

    if (nprocs == 1) then
      !$acc parallel loop collapse(3) default(present)
      do k=1,kmax
        do j=1,jtot
          do i=1,itot
            px(i,j,k) = p(i+1,j+1,k)
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
              if (k+n*konx <= kmax) buffer(ii) = p(i,j,k+n*konx) 
            end do
          end do
        end do
      end do

      call D_MPI_ALLTOALL(buffer, imax*jmax*konx, &
                          commrow, mpierr, lacc=.true.)

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

  subroutine transpose_x_to_z(p, px, buffer)

    real(pois_r), intent(in)  :: px(:,:,:)
    real(pois_r), intent(out) :: p(:,:,:)
    real(pois_r), intent(out) :: buffer(:)

    integer :: i, j, k, n, ii

    if (nprocs == 1) then
      !$acc parallel loop collapse(3) default(present)
      do k = 1, kmax
        do j = 1, jtot
          do i = 1, itot
            p(i+1,j+1,k) = px(i,j,k)
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

      call D_MPI_ALLTOALL(buffer, imax*jmax*konx, &
                          commrow, mpierr, lacc=.true.)

      !$acc parallel loop collapse(4) default(present) private(ii)
      do n = 0, nprocx-1
        do k = 1, konx
          do j = 2, j1
            do i = 2, i1
              ii = (i-1) + (j-2)*imax + (k-1)*imax*jmax + n*imax*jmax*konx
              if (k+n*konx <= kmax) p(i,j,k+n*konx) = buffer(ii)
            end do
          end do
        end do
      end do
    end if

  end subroutine transpose_x_to_z

  subroutine transpose_x_to_y(px, py, buffer)

    real(pois_r), intent(in)  :: px(:,:,:)
    real(pois_r), intent(out) :: py(:,:,:)
    real(pois_r), intent(out) :: buffer(:)

    integer :: i, j, k, n, ii

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

      call D_MPI_ALLTOALL(buffer, iony*jmax*konx, &
                          commcol, mpierr, lacc=.true.)


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

  subroutine transpose_y_to_x(px, py, buffer)

    real(pois_r), intent(in)  :: py(:,:,:)
    real(pois_r), intent(out) :: px(:,:,:)
    real(pois_r), intent(out) :: buffer(:)

    integer :: i, j, k, n, ii

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

      call D_MPI_ALLTOALL(buffer, iony*jmax*konx, &
                          commcol, mpierr, lacc=.true.)

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

  subroutine transpose_y_to_z(py, Fp, buffer)

    real(pois_r), intent(in)  :: py(:,:,:)
    real(pois_r), intent(out) :: Fp(:,:,:)
    real(pois_r), intent(out) :: buffer(:)

    integer :: i, j, k, n, ii

    if (nprocs == 1) then
      !$acc parallel loop collapse(3) default(present)
      do k = 1, kmax
        do j = 1, jtot
          do i = 1, itot
            Fp(i+1,j+1,k) = py(j,k,i)
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

      call D_MPI_ALLTOALL(buffer, iony*jonx*konx, &
                          commrow, mpierr, lacc=.true.)

      !$acc parallel loop collapse(4) default(present) private(ii)
      do n = 0, nprocx-1
        do k = 1, konx
          do j = 1, jonx
            do i = 1, iony
              ii = j + (i-1)*jonx + (k-1)*iony*jonx + n*iony*jonx*konx
              if (k+n*konx <= kmax) Fp(i,j,k+n*konx) = buffer(ii)
            end do
          end do
        end do
      end do

    end if

  end subroutine transpose_y_to_z

  subroutine transpose_z_to_y(py, Fp, buffer)

    real(pois_r), intent(in)  :: Fp(:,:,:)
    real(pois_r), intent(out) :: py(:,:,:)
    real(pois_r), intent(out) :: buffer(:)

    integer :: i, j, k, n, ii

    if (nprocs == 1) then
      !$acc parallel loop collapse(3) default(present)
      do k=1,kmax
        do j=1,jtot
          do i=1,itot
            py(j,k,i) = Fp(i+1,j+1,k)
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
              if (k+n*konx <= kmax) buffer(ii) = Fp(i,j,k+n*konx)
            end do
          end do
        end do
      end do

      call D_MPI_ALLTOALL(buffer, iony*jonx*konx, &
                          commrow, mpierr, lacc=.true.)

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