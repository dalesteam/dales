!> Transposition routines for Poisson solver.
module modtranspose

  use modglobal,    only: itot, jtot, imax, jmax, kmax, i1, j1
  use modmpi,       only: MPI_ALLTOALL, commrow, commcol, nprocs, nprocx, &
                          nprocy, MPI_REAL4, MPI_REAL8, MPI_DATATYPE, &
                          MPI_IN_PLACE
  use modprecision, only: pois_r, longint
  use modtimer,     only: ltimer, timer_tic, timer_toc

  implicit none

  private

  character(len=*), parameter :: modname = 'modtranspose'

  public :: t_transposer

  interface t_transposer
    procedure :: build_transposer
  end interface t_transposer

  ! Derived type that wraps the transpose routines.
  type :: t_transposer
    !> Pencil dimensions.
    integer :: iony, jonx, konx
  contains
    procedure :: get_buffer_size => transpose_get_buffer_size
    procedure :: z_to_x => transpose_z_to_x
    procedure :: x_to_z => transpose_x_to_z
    procedure :: x_to_y => transpose_x_to_y
    procedure :: y_to_x => transpose_y_to_x
    procedure :: y_to_z => transpose_y_to_z
    procedure :: z_to_y => transpose_z_to_y
  end type t_transposer

#if POIS_PRECISION == 64
  type(MPI_DATATYPE), parameter :: MPI_DTYPE = MPI_REAL8
#else
  type(MPI_DATATYPE), parameter :: MPI_DTYPE = MPI_REAL4
#endif

contains

  !> Initialize transpose object.
  pure function build_transposer() result(self)

    type(t_transposer) :: self

    self%iony = itot / nprocy
    if (mod(itot, nprocy) > 0) self%iony = self%iony + 1

    self%jonx = jtot / nprocx
    if (mod(jtot, nprocx) > 0) self%jonx = self%jonx + 1

    self%konx = kmax / nprocx
    if (mod(kmax, nprocx) > 0) self%konx = self%konx + 1

  end function build_transposer

  !> Compute the minimum size of the workspace for transposing.
  !!
  !! @return minimum buffer size for transposing.
  pure function transpose_get_buffer_size(self) result(size)
    
    class(t_transposer), intent(in) :: self

    integer(longint) :: size_x, size_y, size_z, size

    ! x-contiguous pencils
    size_x = itot * jmax * self%konx

    ! y-contiguous pencils
    size_y = self%iony * jtot * self%konx

    ! z-contiguous pencils
    size_z = imax * jmax * kmax

    size = max(size_x, size_y, size_z)

  end function transpose_get_buffer_size

  !> Tranpose z-contiguous pencils to x-contiguous pencils.
  !!
  !! @param[in] pz Input data.
  !! @param[out] px Output data.
  !! @param[out] buffer Buffer for transposing.
  subroutine transpose_z_to_x(self, pz, px, buffer)

    class(t_transposer), intent(in)  :: self
    real(pois_r),        intent(in)  :: pz(:,:,:)
    real(pois_r),        intent(out) :: px(:,:,:)
    real(pois_r),        intent(out) :: buffer(:)

    character(len=*), parameter :: routine = modname//'/transpose_z_to_x'

    integer :: i, j, k, n, ii
    integer :: mpierr

    if (ltimer) call timer_tic(routine, 2)

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
        do k = 1, self%konx
          do j = 2, j1
            do i = 2, i1
              ii = (i-1) + (j-2)*imax + (k-1)*imax*jmax + n*imax*jmax*self%konx
              if (k+n*self%konx <= kmax) buffer(ii) = pz(i,j,k+n*self%konx) 
            end do
          end do
        end do
      end do

      !$acc host_data use_device(buffer)
      call MPI_ALLTOALL(MPI_IN_PLACE, 0, MPI_DTYPE, &
                        buffer, imax*jmax*self%konx, MPI_DTYPE, commrow, mpierr)
      !$acc end host_data

      !$acc parallel loop collapse(4) default(present) private(ii)
      do n = 0, nprocx-1
        do k = 1, self%konx
          do j = 1, jmax
            do i = 1, imax
              ii = i + (j-1)*imax + (k-1)*imax*jmax + n*imax*jmax*self%konx
              px(i+n*imax,j,k) = buffer(ii)
            end do
          end do
        end do
      end do
    end if

    if(ltimer) call timer_toc(routine)

  end subroutine transpose_z_to_x

  !> Tranpose x-contiguous pencils to z-contiguous pencils.
  !!
  !! @param[in] px Input data.
  !! @param[out] pz Output data.
  !! @param[out] buffer Buffer for transposing.
  subroutine transpose_x_to_z(self, px, pz, buffer)

    class(t_transposer), intent(in)  :: self
    real(pois_r),        intent(in)  :: px(:,:,:)
    real(pois_r),        intent(out) :: pz(:,:,:)
    real(pois_r),        intent(out) :: buffer(:)

    character(len=*), parameter :: routine = modname//'/transpose_x_to_z'

    integer :: i, j, k, n, ii
    integer :: mpierr

    if (ltimer) call timer_tic(routine, 2)

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
        do k = 1, self%konx
          do j = 1, jmax
            do i = 1, imax
              ii = i + (j-1)*imax + (k-1)*imax*jmax + n*imax*jmax*self%konx
              buffer(ii) = px(i+n*imax,j,k)
            end do
          end do
        end do
      end do

      !$acc host_data use_device(buffer)
      call MPI_ALLTOALL(MPI_IN_PLACE, 0, MPI_DTYPE, &
                        buffer, imax*jmax*self%konx, MPI_DTYPE, commrow, mpierr)
      !$acc end host_data

      !$acc parallel loop collapse(4) default(present) private(ii)
      do n = 0, nprocx-1
        do k = 1, self%konx
          do j = 2, j1
            do i = 2, i1
              ii = (i-1) + (j-2)*imax + (k-1)*imax*jmax + n*imax*jmax*self%konx
              if (k+n*self%konx <= kmax) pz(i,j,k+n*self%konx) = buffer(ii)
            end do
          end do
        end do
      end do
    end if

    if (ltimer) call timer_toc(routine)

  end subroutine transpose_x_to_z

  !> Tranpose x-contiguous pencils to y-contiguous pencils.
  !!
  !! @param[in] px Input data.
  !! @param[out] py Output data.
  !! @param[out] buffer Buffer for transposing.
  subroutine transpose_x_to_y(self, px, py, buffer)

    class(t_transposer), intent(in)  :: self
    real(pois_r),        intent(in)  :: px(:,:,:)
    real(pois_r),        intent(out) :: py(:,:,:)
    real(pois_r),        intent(out) :: buffer(:)

    character(len=*), parameter :: routine = modname//'/transpose_x_to_y'

    integer :: i, j, k, n, ii
    integer :: mpierr

    if (ltimer) call timer_tic(routine, 2)

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
        do k = 1, self%konx
          do j = 1, jmax
            do i = 1, self%iony
              ii = i + (j-1)*self%iony + (k-1)*self%iony*jmax + n*self%iony*jmax*self%konx
              if (i <= itot) buffer(ii) = px(i+n*self%iony,j,k)
            end do
          end do
        end do
      end do

      !$acc host_data use_device(buffer)
      call MPI_ALLTOALL(MPI_IN_PLACE, 0, MPI_DTYPE, &
                        buffer, self%iony*jmax*self%konx, MPI_DTYPE, &
                        commcol, mpierr)
      !$acc end host_data

      !$acc parallel loop collapse(4) default(present) private(ii)
      do n = 0, nprocy-1
        do k = 1, self%konx
          do i = 1, self%iony
            do j = 1, jmax
              ii = i + (j-1)*self%iony + (k-1)*self%iony*jmax + n*self%iony*jmax*self%konx
              py(j+n*jmax,k,i) = buffer(ii)
            end do
          end do
        end do
      end do
    end if

    if (ltimer) call timer_toc(routine)

  end subroutine transpose_x_to_y

  !> Tranpose y-contiguous pencils to x-contiguous pencils.
  !!
  !! @param[in] py Input data.
  !! @param[out] px Output data.
  !! @param[out] buffer Buffer for transposing.
  subroutine transpose_y_to_x(self, py, px, buffer)

    class(t_transposer), intent(in)  :: self
    real(pois_r),        intent(in)  :: py(:,:,:)
    real(pois_r),        intent(out) :: px(:,:,:)
    real(pois_r),        intent(out) :: buffer(:)

    character(len=*), parameter :: routine = modname//'/transpose_y_to_x'

    integer :: i, j, k, n, ii
    integer :: mpierr

    if (ltimer) call timer_tic(routine, 2)

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
        do k = 1, self%konx
          do i = 1, self%iony
            do j = 1, jmax
              ii = i + (j-1)*self%iony + (k-1)*self%iony*jmax + n*self%iony*jmax*self%konx
              buffer(ii) = py(j+n*jmax,k,i)
            end do
          end do
        end do
      end do

      !$acc host_data use_device(buffer)
      call MPI_ALLTOALL(MPI_IN_PLACE, 0, MPI_DTYPE, &
                        buffer, self%iony*jmax*self%konx, MPI_DTYPE, &
                        commcol, mpierr)
      !$acc end host_data

      !$acc parallel loop collapse(4) default(present) private(ii)
      do n = 0, nprocy-1
        do k = 1, self%konx
          do j = 1, jmax
            do i = 1, self%iony
              ii = i + (j-1)*self%iony + (k-1)*self%iony*jmax + n*self%iony*jmax*self%konx
              if (i+n*self%iony <= itot) px(i+n*self%iony,j,k) = buffer(ii)
            end do
          end do
        end do
      end do
    end if

    if (ltimer) call timer_toc(routine)

  end subroutine transpose_y_to_x

  !> Tranpose y-contiguous pencils to z-contiguous pencils.
  !!
  !! @param[in] py Input data.
  !! @param[out] pz Output data.
  !! @param[out] buffer Buffer for transposing.
  subroutine transpose_y_to_z(self, py, pz, buffer)

    class(t_transposer), intent(in)  :: self
    real(pois_r),        intent(in)  :: py(:,:,:)
    real(pois_r),        intent(out) :: pz(:,:,:)
    real(pois_r),        intent(out) :: buffer(:)

    character(len=*), parameter :: routine = modname//'/transpose_y_to_z'

    integer :: i, j, k, n, ii
    integer :: mpierr

    if (ltimer) call timer_tic(routine, 2)

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
        do k = 1, self%konx
          do i = 1, self%iony
            do j = 1, self%jonx
              ii = j + (i-1)*self%jonx + (k-1)*self%iony*self%jonx + n*self%iony*self%jonx*self%konx
              if (j+n*self%jonx <= jtot) buffer(ii) = py(j+n*self%jonx,k,i)
            end do
          end do
        end do
      end do

      !$acc host_data use_device(buffer)
      call MPI_ALLTOALL(MPI_IN_PLACE, 0, MPI_DTYPE, &
                        buffer, self%iony*self%jonx*self%konx, MPI_DTYPE, &
                        commrow, mpierr)
      !$acc end host_data

      !$acc parallel loop collapse(4) default(present) private(ii)
      do n = 0, nprocx-1
        do k = 1, self%konx
          do j = 1, self%jonx
            do i = 1, self%iony
              ii = j + (i-1)*self%jonx + (k-1)*self%iony*self%jonx + n*self%iony*self%jonx*self%konx
              if (k+n*self%konx <= kmax) pz(i,j,k+n*self%konx) = buffer(ii)
            end do
          end do
        end do
      end do
    end if

    if (ltimer) call timer_toc(routine)

  end subroutine transpose_y_to_z

  !> Tranpose y-contiguous pencils to x-contiguous pencils.
  !!
  !! @param[in] py Input data.
  !! @param[out] px Output data.
  !! @param[out] buffer Buffer for transposing.
  subroutine transpose_z_to_y(self, pz, py, buffer)

    class(t_transposer), intent(in)  :: self
    real(pois_r),        intent(in)  :: pz(:,:,:)
    real(pois_r),        intent(out) :: py(:,:,:)
    real(pois_r),        intent(out) :: buffer(:)

    character(len=*), parameter :: routine = modname//'/transpose_z_to_y'

    integer :: i, j, k, n, ii
    integer :: mpierr

    if (ltimer) call timer_tic(routine, 2)

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
        do k = 1, self%konx
          do j = 1, self%jonx
            do i = 1, self%iony
              ii = j + (i-1)*self%jonx + (k-1)*self%iony*self%jonx + n*self%iony*self%jonx*self%konx
              if (k+n*self%konx <= kmax) buffer(ii) = pz(i,j,k+n*self%konx)
            end do
          end do
        end do
      end do

      !$acc host_data use_device(buffer)
      call MPI_ALLTOALL(MPI_IN_PLACE, 0, MPI_DTYPE, &
                        buffer, self%iony*self%jonx*self%konx, MPI_DTYPE, &
                        commrow, mpierr)
      !$acc end host_data

      !$acc parallel loop collapse(4) default(present) private(ii)
      do n = 0, nprocx-1
        do k = 1, self%konx
          do i = 1, self%iony
            do j = 1, self%jonx
              ii = j + (i-1)*self%jonx + (k-1)*self%iony*self%jonx + n*self%iony*self%jonx*self%konx
              if (j+n*self%jonx <= jtot) py(j+n*self%jonx,k,i) = buffer(ii)
            end do
          end do
        end do
      end do
    end if

    if (ltimer) call timer_toc(routine)

  end subroutine transpose_z_to_y

end module modtranspose