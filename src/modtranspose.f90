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
!!$omp declare mapper (t_transposer::x) map ( &
!!$omp  x%iony &
!!$omp , x%jonx &
!!$omp , x%konx &
!!$omp )

  type(MPI_DATATYPE) :: MPI_DTYPE

contains

  !> Initialize transpose object.
  pure function build_transposer() result(this)

    type(t_transposer) :: this

    this%iony = itot / nprocy
    if (mod(itot, nprocy) > 0) this%iony = this%iony + 1

    this%jonx = jtot / nprocx
    if (mod(jtot, nprocx) > 0) this%jonx = this%jonx + 1

    this%konx = kmax / nprocx
    if (mod(kmax, nprocx) > 0) this%konx = this%konx + 1

  end function build_transposer

  !> Compute the minimum size of the workspace for transposing.
  !!
  !! @return minimum buffer size for transposing.
  pure function transpose_get_buffer_size(this) result(size)
    
    class(t_transposer), intent(in) :: this

    integer(longint) :: size_x, size_y, size_z, size_z2, size

    ! x-contiguous pencils
    size_x = itot * jmax * this%konx

    ! y-contiguous pencils
    size_y = this%iony * jtot * this%konx

    ! z-contiguous pencils
    size_z = imax * jmax * kmax

    ! z-contiguious pencils, other flavour
    size_z2 = this%iony * this%jonx * this%konx * nprocx

    size = max(size_x, size_y, size_z, size_z2)

  end function transpose_get_buffer_size

  !> Tranpose z-contiguous pencils to x-contiguous pencils.
  !!
  !! @param[in] pz Input data.
  !! @param[out] px Output data.
  !! @param[out] buffer Buffer for transposing.
  subroutine transpose_z_to_x(this, pz, px, buffer)

    class(t_transposer),          intent(in)  :: this
    real(pois_r),        pointer, intent(in)  :: pz(:,:,:)
    real(pois_r),        pointer, intent(in)  :: px(:,:,:)
    real(pois_r),                 intent(out) :: buffer(:)

    character(len=*), parameter :: routine = modname//'/transpose_z_to_x'

    integer :: i, j, k, n, ii
    integer :: n1, n2, n3
    integer :: mpierr

#if POIS_PRECISION == 64
     MPI_DTYPE = MPI_REAL8
#else
     MPI_DTYPE = MPI_REAL4
#endif

    if (ltimer) call timer_tic(routine, 2)

    if (nprocs == 1) then
      !$acc parallel loop collapse(3) default(present)
!!$omp target teams loop collapse(3) defaultmap(present:aggregate)&
!!$omp defaultmap(present:allocatable)
      do k=1,kmax
        do j=1,jtot
          do i=1,itot
            px(i,j,k) = pz(i+1,j+1,k)
          end do
        end do
      end do
    else

      n1 = imax
      n2 = jmax
      n3 = this%konx
      
      !$acc parallel loop collapse(4) default(present) private(ii)
!!$omp target teams loop private(ii) collapse(4)&
!!$omp defaultmap(present:aggregate) defaultmap(present:allocatable)
      do n = 0, nprocx-1
        do k = 1, n3
          do j = 1, n2
            do i = 1, n1
              ii = i + (j-1)*n1 + (k-1)*n1*n2 + n*n1*n2*n3
              if (k+n*n3 <= kmax) buffer(ii) = pz(i+1,j+1,k+n*n3) 
            end do
          end do
        end do
      end do

      !$acc host_data use_device(buffer)
!!$omp target update from(buffer)
      call MPI_ALLTOALL(MPI_IN_PLACE, 0, MPI_DTYPE, &
                        buffer, n1*n2*n3, MPI_DTYPE, commrow, mpierr)
      !$acc end host_data
!$omp target update to(buffer)

      !$acc parallel loop collapse(4) default(present) private(ii)
!!$omp target teams loop private(ii) collapse(4)&
!!$omp defaultmap(present:aggregate) defaultmap(present:allocatable)
      do n = 0, nprocx-1
        do k = 1, n3
          do j = 1, n2
            do i = 1, n1
              ii = i + (j-1)*n1 + (k-1)*n1*n2 + n*n1*n2*n3
              px(i+n*n1,j,k) = buffer(ii)
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
  subroutine transpose_x_to_z(this, px, pz, buffer)

    class(t_transposer),          intent(in)  :: this
    real(pois_r),        pointer, intent(in)  :: px(:,:,:)
    real(pois_r),        pointer, intent(in)  :: pz(:,:,:)
    real(pois_r),                 intent(out) :: buffer(:)

    character(len=*), parameter :: routine = modname//'/transpose_x_to_z'

    integer :: i, j, k, n, ii
    integer :: n1, n2, n3
    integer :: mpierr

    if (ltimer) call timer_tic(routine, 2)

    if (nprocs == 1) then
      !$acc parallel loop collapse(3) default(present)
!!$omp target teams loop collapse(3) defaultmap(present:aggregate)&
!!$omp defaultmap(present:allocatable)
      do k = 1, kmax
        do j = 1, jtot
          do i = 1, itot
            pz(i+1,j+1,k) = px(i,j,k)
          end do
        end do
      end do
    else
      
      n1 = imax
      n2 = jmax
      n3 = this%konx

      !$acc parallel loop collapse(4) default(present) private(ii)
!!$omp target teams loop private(ii) collapse(4)&
!!$omp defaultmap(present:aggregate) defaultmap(present:allocatable)
      do n = 0, nprocx-1
        do k = 1, n3
          do j = 1, n2
            do i = 1, n1
              ii = i + (j-1)*n1 + (k-1)*n1*n2 + n*n1*n2*n3
              buffer(ii) = px(i+n*n1,j,k)
            end do
          end do
        end do
      end do

      !$acc host_data use_device(buffer)
!!$omp target update from(buffer)
      call MPI_ALLTOALL(MPI_IN_PLACE, 0, MPI_DTYPE, &
                        buffer, n1*n2*n3, MPI_DTYPE, commrow, mpierr)
      !$acc end host_data
!$omp target update to(buffer)

      !$acc parallel loop collapse(4) default(present) private(ii)
!!$omp target teams loop private(ii) collapse(4)&
!!$omp defaultmap(present:aggregate) defaultmap(present:allocatable)
      do n = 0, nprocx-1
        do k = 1, n3
          do j = 1, n2
            do i = 1, n1
              ii = i + (j-1)*n1 + (k-1)*n1*n2 + n*n1*n2*n3
              if (k+n*n3 <= kmax) pz(i+1,j+1,k+n*n3) = buffer(ii)
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
  subroutine transpose_x_to_y(this, px, py, buffer)

    class(t_transposer),          intent(in)  :: this
    real(pois_r),        pointer, intent(in)  :: px(:,:,:)
    real(pois_r),        pointer, intent(in)  :: py(:,:,:)
    real(pois_r),                 intent(out) :: buffer(:)

    character(len=*), parameter :: routine = modname//'/transpose_x_to_y'

    integer :: i, j, k, n, ii
    integer :: n1, n2, n3
    integer :: mpierr

    if (ltimer) call timer_tic(routine, 2)

    if (nprocs == 1) then
      !$acc parallel loop collapse(3) default(present) private(ii)
!!$omp target teams loop private(ii) collapse(3)&
!!$omp defaultmap(present:aggregate) defaultmap(present:allocatable)
      do k = 1, kmax
        do j = 1, jtot
          do i = 1, itot
            ii = i + (j-1)*itot + (k-1)*itot*jtot
            buffer(ii) = px(i,j,k)
          end do
        end do
      end do

      !$acc parallel loop collapse(3) default(present) private(ii)
!!$omp target teams loop private(ii) collapse(3)&
!!$omp defaultmap(present:aggregate) defaultmap(present:allocatable)
      do k = 1, kmax
        do j = 1, jtot
         do i = 1, itot
            ii = i + (j-1)*itot + (k-1)*itot*jtot
            py(j,k,i) = buffer(ii)
          end do
        end do
      end do
    else

      n1 = this%iony
      n2 = jmax
      n3 = this%konx

      !$acc parallel loop collapse(4) default(present) private(ii)
!!$omp target teams loop private(ii) collapse(4)&
!!$omp defaultmap(present:aggregate) defaultmap(present:allocatable)
      do n = 0, nprocy-1
        do k = 1, n3
          do j = 1, n2
            do i = 1, n1
              ii = i + (j-1)*n1 + (k-1)*n1*n2 + n*n1*n2*n3
              if (i+n*n1 <= itot) buffer(ii) = px(i+n*n1,j,k)
            end do
          end do
        end do
      end do

      !$acc host_data use_device(buffer)
!!$omp target update from(buffer)
      call MPI_ALLTOALL(MPI_IN_PLACE, 0, MPI_DTYPE, &
                        buffer, n1*n2*n3, MPI_DTYPE, &
                        commcol, mpierr)
      !$acc end host_data
!$omp target update to(buffer)

      !$acc parallel loop collapse(4) default(present) private(ii)
!!$omp target teams loop private(ii) collapse(4)&
!!$omp defaultmap(present:aggregate) defaultmap(present:allocatable)
      do n = 0, nprocy-1
        do k = 1, n3
          do i = 1, n1
            do j = 1, n2
              ii = i + (j-1)*n1 + (k-1)*n1*n2 + n*n1*n2*n3
              py(j+n*n2,k,i) = buffer(ii)
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
  subroutine transpose_y_to_x(this, py, px, buffer)

    class(t_transposer),          intent(in)  :: this
    real(pois_r),        pointer, intent(in)  :: py(:,:,:)
    real(pois_r),        pointer, intent(in)  :: px(:,:,:)
    real(pois_r),                 intent(out) :: buffer(:)

    character(len=*), parameter :: routine = modname//'/transpose_y_to_x'

    integer :: i, j, k, n, ii
    integer :: n1, n2, n3
    integer :: mpierr

    if (ltimer) call timer_tic(routine, 2)

    if (nprocs == 1) then
      !$acc parallel loop collapse(3) default(present) private(ii)
!!$omp target teams loop private(ii) collapse(3)&
!!$omp defaultmap(present:aggregate) defaultmap(present:allocatable)
      do k = 1, kmax
        do j = 1, jtot
          do i = 1, itot
            ii = j + (i-1)*jtot + (k-1)*itot*jtot
            buffer(ii) = py(j,k,i)
          end do
        end do
      end do

      !$acc parallel loop collapse(3) default(present) private(ii)
!!$omp target teams loop private(ii) collapse(3)&
!!$omp defaultmap(present:aggregate) defaultmap(present:allocatable)
      do k = 1, kmax
        do j = 1, jtot
          do i = 1, itot
            ii = j + (i-1)*jtot + (k-1)*itot*jtot
            px(i,j,k) = buffer(ii)
          end do
        end do
      end do
    else

      n1 = this%iony
      n2 = jmax
      n3 = this%konx

      !$acc parallel loop collapse(4) default(present) private(ii)
!!$omp target teams loop private(ii) collapse(4)&
!!$omp defaultmap(present:aggregate) defaultmap(present:allocatable)
      do n = 0, nprocy-1
        do k = 1, n3
          do i = 1, n1
            do j = 1, n2
              ii = i + (j-1)*n1 + (k-1)*n1*n2 + n*n1*n2*n3
              buffer(ii) = py(j+n*n2,k,i)
            end do
          end do
        end do
      end do

      !$acc host_data use_device(buffer)
!!$omp target update from(buffer)
      call MPI_ALLTOALL(MPI_IN_PLACE, 0, MPI_DTYPE, &
                        buffer, n1*n2*n3, MPI_DTYPE, &
                        commcol, mpierr)
      !$acc end host_data
!$omp target update to(buffer)

      !$acc parallel loop collapse(4) default(present) private(ii)
!!$omp target teams loop private(ii) collapse(4)&
!!$omp defaultmap(present:aggregate) defaultmap(present:allocatable)
      do n = 0, nprocy-1
        do k = 1, n3
          do j = 1, n2
            do i = 1, n1
              ii = i + (j-1)*n1 + (k-1)*n1*n2 + n*n1*n2*n3
              if (i+n*n1 <= itot) px(i+n*n1,j,k) = buffer(ii)
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
  subroutine transpose_y_to_z(this, py, pz, buffer)

    class(t_transposer),          intent(in)  :: this
    real(pois_r),        pointer, intent(in)  :: py(:,:,:)
    real(pois_r),        pointer, intent(in)  :: pz(:,:,:)
    real(pois_r),                 intent(out) :: buffer(:)

    character(len=*), parameter :: routine = modname//'/transpose_y_to_z'

    integer :: i, j, k, n, ii
    integer :: n1, n2, n3
    integer :: mpierr

    if (ltimer) call timer_tic(routine, 2)

    if (nprocs == 1) then
      !$acc parallel loop collapse(3) default(present)
!!$omp target teams loop collapse(3) defaultmap(present:aggregate)&
!!$omp defaultmap(present:allocatable)
      do k = 1, kmax
        do j = 1, jtot
          do i = 1, itot
            pz(i+1,j+1,k) = py(j,k,i)
          end do
        end do
      end do
    else

      n1 = this%jonx
      n2 = this%iony
      n3 = this%konx

      !$acc parallel loop collapse(4) default(present) private(ii)
!!$omp target teams loop private(ii) collapse(4)&
!!$omp defaultmap(present:aggregate) defaultmap(present:allocatable)
      do n = 0, nprocx-1
        do k = 1, n3
          do i = 1, n2
            do j = 1, n1
              ii = j + (i-1)*n1 + (k-1)*n1*n2 + n*n1*n2*n3
              if (j+n*n1 <= jtot) buffer(ii) = py(j+n*n1,k,i)
            end do
          end do
        end do
      end do

      !$acc host_data use_device(buffer)
!!$omp target update from(buffer)
      call MPI_ALLTOALL(MPI_IN_PLACE, 0, MPI_DTYPE, &
                        buffer, n1*n2*n3, MPI_DTYPE, &
                        commrow, mpierr)
      !$acc end host_data
!$omp target update to(buffer)

      !$acc parallel loop collapse(4) default(present) private(ii)
!!$omp target teams loop private(ii) collapse(4)&
!!$omp defaultmap(present:aggregate) defaultmap(present:allocatable)
      do n = 0, nprocx-1
        do k = 1, n3
          do j = 1, n1
            do i = 1, n2
              ii = j + (i-1)*n1 + (k-1)*n1*n2 + n*n1*n2*n3
              if (k+n*n3 <= kmax) pz(i,j,k+n*n3) = buffer(ii)
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
  subroutine transpose_z_to_y(this, pz, py, buffer)

    class(t_transposer),          intent(in)  :: this
    real(pois_r),        pointer, intent(in)  :: pz(:,:,:)
    real(pois_r),        pointer, intent(in)  :: py(:,:,:)
    real(pois_r),                 intent(out) :: buffer(:)

    character(len=*), parameter :: routine = modname//'/transpose_z_to_y'

    integer :: i, j, k, n, ii
    integer :: n1, n2, n3
    integer :: mpierr

    if (ltimer) call timer_tic(routine, 2)

    if (nprocs == 1) then
      !$acc parallel loop collapse(3) default(present)
!!$omp target teams loop collapse(3) defaultmap(present:aggregate)&
!!$omp defaultmap(present:allocatable)
      do k=1,kmax
        do j=1,jtot
          do i=1,itot
            py(j,k,i) = pz(i+1,j+1,k)
          end do
        end do
      end do
    else

      n1 = this%jonx
      n2 = this%iony
      n3 = this%konx

      !$acc parallel loop collapse(4) default(present) private(ii)
!!$omp target teams loop private(ii) collapse(4)&
!!$omp defaultmap(present:aggregate) defaultmap(present:allocatable)
      do n = 0, nprocx-1
        do k = 1, n3
          do j = 1, n1
            do i = 1, n2
              ii = j + (i-1)*n1 + (k-1)*n1*n2 + n*n1*n2*n3
              if (k+n*n3 <= kmax) buffer(ii) = pz(i,j,k+n*n3)
            end do
          end do
        end do
      end do

      !$acc host_data use_device(buffer)
!!$omp target update from(buffer)
      call MPI_ALLTOALL(MPI_IN_PLACE, 0, MPI_DTYPE, &
                        buffer, n1*n2*n3, MPI_DTYPE, &
                        commrow, mpierr)
      !$acc end host_data
!$omp target update to(buffer)

      !$acc parallel loop collapse(4) default(present) private(ii)
!!$omp target teams loop private(ii) collapse(4)&
!!$omp defaultmap(present:aggregate) defaultmap(present:allocatable)
      do n = 0, nprocx-1
        do k = 1, n3
          do i = 1, n2
            do j = 1, n1
              ii = j + (i-1)*n1 + (k-1)*n1*n2 + n*n1*n2*n3
              if (j+n*n1 <= jtot) py(j+n*n1,k,i) = buffer(ii)
            end do
          end do
        end do
      end do
    end if

    if (ltimer) call timer_toc(routine)

  end subroutine transpose_z_to_y

end module modtranspose
