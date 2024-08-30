module modtranspose

  use modprecision, only: pois_r
  use modglobal,    only: itot, jtot, kmax, i1, j1, k1, ih, jh, imax, jmax
  use modmpi,       only: D_MPI_ALLTOALL, commrow, commcol, nprocs, &
    &                     nprocx, nprocy, mpierr

  implicit none

  private

  public :: transpose_a1, transpose_a1inv
  public :: transpose_a2, transpose_a2inv
  public :: transpose_a3, transpose_a3inv

contains

! Transpose functions:
!   Data is stored such that a whole column, 1:kmax, is on one node (pencil layout).
!   These funtions redistribute the data over the MPI nodes so that a whole dimension
!   is on one processor.
!   To facilitate processing after the MPI transpose, data is further transposed locally
!   such that the complete dimension, one of itot, jtot, kmax, is consecutive (ie. fastest).
!
!   The figure below shows the transpose between dimension 0 and 1, or 'k' and 'i':
!
!
!         /-------------/|                                /-------------/|
!        /../          / |                               /             / |
!  kmax  |------------|  |                               |------------|  |
!        |..|         |  |            ==>                |            |  |
!        |..|         |  |                         konx  |____________|. |
!        |..|         | /    jmax                        |............|./    jmax
!   1    |--|--|--|---|/   1                        1    |------------|/   1
!
!        1  imax                                         1            itot
!
!  transpose_a1: dimensions 0 and 2:
!            p012(2-ih:i1+ih,2-jh:j1+jh,kmax) <=> p210(itot,jmax,konx)
!
!  transpose_a2: dimensions 1 and 2:
!            p210(      itot,      jmax,konx) <=> p201(jtot,konx,iony)
!
!  transpose_a3: dimensions 0 and 2:
!            p201(      jtot,      konx,iony) <=> p102(iony,jonx,kmax)
!

  subroutine transpose_a1(p, px, konx, bufin, bufout)
  
    real(pois_r), pointer, intent(in)  :: p(:,:,:)
    real(pois_r), pointer, intent(out) :: px(:,:,:)
    integer,               intent(in)  :: konx
    real(pois_r),          intent(inout)  :: bufin(:)
    real(pois_r),          intent(inout)  :: bufout(:)
  
    integer :: i, j, k, n, ii
  
    if (nprocs == 1) then
      !$acc parallel loop collapse(3) default(present)
      do k = 1, kmax
        do j = 1, jtot
          do i = 1, itot
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
              if (k+n*konx <= kmax) bufin(ii) = p(i,j,k+n*konx) 
            end do
          end do
        end do
      end do
  
      !$acc host_data use_device(workspace_0, workspace_1)
      call D_MPI_ALLTOALL(bufin, imax*jmax*konx, &
                          bufout, imax*jmax*konx, &
                          commrow, mpierr)
      !$acc end host_data
  
      !$acc parallel loop collapse(4) default(present) private(ii)
      do n = 0, nprocx-1
        do k = 1, konx
          do j = 1, jmax
            do i = 1, imax
              ii = i + (j-1)*imax + (k-1)*imax*jmax + n*imax*jmax*konx
              px(i+n*imax,j,k) = bufout(ii)
            end do
          end do
        end do
      end do
    end if
  
  end subroutine transpose_a1
  
  subroutine transpose_a1inv(p, px, konx, bufin, bufout)
  
    real(pois_r), pointer, intent(in)  :: px(:,:,:)
    real(pois_r), pointer, intent(out) :: p(:,:,:)
    integer,               intent(in)  :: konx
    real(pois_r),          intent(inout)  :: bufin(:)
    real(pois_r),          intent(inout)  :: bufout(:)
  
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
              bufin(ii) = px(i+n*imax,j,k)
            end do
          end do
        end do
      end do
  
      !$acc host_data use_device(workspace_0, workspace_1)
      call D_MPI_ALLTOALL(bufin, imax*jmax*konx, &
                          bufout, imax*jmax*konx, &
                          commrow, mpierr)
      !$acc end host_data
  
      !$acc parallel loop collapse(4) default(present) private(ii)
      do n = 0, nprocx-1
        do k = 1, konx
          do j = 2, j1
            do i = 2, i1
              ii = (i-1) + (j-2)*imax + (k-1)*imax*jmax + n*imax*jmax*konx
              if (k+n*konx <= kmax) p(i,j,k+n*konx) = bufout(ii)
            end do
          end do
        end do
      end do
    end if
  
  end subroutine transpose_a1inv
  
  subroutine transpose_a2(px, py, iony, konx, bufin, bufout)
  
    real(pois_r), pointer, intent(in)  :: px(:,:,:)
    real(pois_r), pointer, intent(out) :: py(:,:,:)
    integer,               intent(in)  :: iony
    integer,               intent(in)  :: konx
    real(pois_r),          intent(inout)  :: bufin(:)
    real(pois_r),          intent(inout)  :: bufout(:)
  
    integer :: i, j, k, n, ii
  
    if (nprocs == 1) then
      !$acc parallel loop collapse(3) default(present) private(ii)
      do k = 1, kmax
        do j = 1, jtot
          do i = 1, itot
            ii = i + (j-1)*itot + (k-1)*itot*jtot
            bufin(ii) = px(i,j,k)
          end do
        end do
      end do
  
      !$acc parallel loop collapse(3) default(present) private(ii)
      do k = 1, kmax
        do j = 1, jtot
         do i = 1, itot
            ii = i + (j-1)*itot + (k-1)*itot*jtot
            py(j,k,i) = bufout(ii)
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
              if (i <= itot) bufin(ii) = px(i+n*iony,j,k)
            end do
          end do
        end do
      end do
  
      !$acc host_data use_device(workspace_0, workspace_1)
      call D_MPI_ALLTOALL(bufin, iony*jmax*konx, &
                          bufout, iony*jmax*konx, &
                          commcol, mpierr)
      !$acc end host_data
  
      !$acc parallel loop collapse(4) default(present) private(ii)
      do n = 0, nprocy-1
        do k = 1, konx
          do i = 1, iony
            do j = 1, jmax
              ii = i + (j-1)*iony + (k-1)*iony*jmax + n*iony*jmax*konx
              py(j+n*jmax,k,i) = bufout(ii)
            end do
          end do
        end do
      end do
              
    end if
  
  end subroutine transpose_a2
  
  subroutine transpose_a2inv(px, py, iony, konx, bufin, bufout)

    real(pois_r), pointer, intent(in)  :: px(:,:,:)
    real(pois_r), pointer, intent(out) :: py(:,:,:)
    integer,               intent(in)  :: iony
    integer,               intent(in)  :: konx
    real(pois_r),          intent(inout)  :: bufin(:)
    real(pois_r),          intent(inout)  :: bufout(:)
  
    integer :: i, j, k, n, ii
  
    if (nprocs == 1) then
      !$acc parallel loop collapse(3) default(present) private(ii)
      do k = 1, kmax
        do j = 1, jtot
          do i = 1, itot
            ii = j + (i-1)*jtot + (k-1)*itot*jtot
            bufin(ii) = py(j,k,i)
          end do
        end do
      end do
  
      !$acc parallel loop collapse(3) default(present) private(ii)
      do k = 1, kmax
        do j = 1, jtot
          do i = 1, itot
            ii = j + (i-1)*jtot + (k-1)*itot*jtot
            px(i,j,k) = bufin(ii)
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
              bufin(ii) = py(j+n*jmax,k,i)
            end do
          end do
        end do
      end do
  
      !$acc host_data use_device(workspace_0, workspace_1)
      call D_MPI_ALLTOALL(bufin, iony*jmax*konx, &
                          bufout, iony*jmax*konx, &
                          commcol, mpierr)
      !$acc end host_data
  
      !$acc parallel loop collapse(4) default(present) private(ii)
      do n = 0, nprocy-1
        do k = 1, konx
          do j = 1, jmax
            do i = 1, iony
              ii = i + (j-1)*iony + (k-1)*iony*jmax + n*iony*jmax*konx
              if (i+n*iony <= itot) px(i+n*iony,j,k) = bufout(ii)
            end do
          end do
        end do
      end do
    end if
  
  end subroutine transpose_a2inv
  
  subroutine transpose_a3(py, Fp, iony, jonx, konx, bufin, bufout)

    real(pois_r), pointer, intent(in)  :: py(:,:,:)
    real(pois_r), pointer, intent(out) :: Fp(:,:,:)
    integer,               intent(in)  :: iony
    integer,               intent(in)  :: jonx
    integer,               intent(in)  :: konx
    real(pois_r),          intent(inout)  :: bufin(:)
    real(pois_r),          intent(inout)  :: bufout(:)
  
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
              if (j+n*jonx <= jtot) bufin(ii) = py(j+n*jonx,k,i)
            end do
          end do
        end do
      end do
  
      !$acc host_data use_device(workspace_0, workspace_1)
      call D_MPI_ALLTOALL(bufin, iony*jonx*konx, &
                          bufout, iony*jonx*konx, &
                          commrow, mpierr)
      !$acc end host_data
  
      !$acc parallel loop collapse(4) default(present) private(ii)
      do n = 0, nprocx-1
        do k = 1, konx
          do j = 1, jonx
            do i = 1, iony
              ii = j + (i-1)*jonx + (k-1)*iony*jonx + n*iony*jonx*konx
              if (k+n*konx <= kmax) Fp(i,j,k+n*konx) = bufout(ii)
            end do
          end do
        end do
      end do
  
    end if
  
  end subroutine transpose_a3
  
  subroutine transpose_a3inv(py, Fp, iony, jonx, konx, bufin, bufout)
  
    real(pois_r), pointer, intent(in)  :: Fp(:,:,:)
    real(pois_r), pointer, intent(out) :: py(:,:,:)
    integer,               intent(in)  :: iony
    integer,               intent(in)  :: jonx
    integer,               intent(in)  :: konx
    real(pois_r),          intent(inout)  :: bufin(:)
    real(pois_r),          intent(inout)  :: bufout(:)
  
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
              if (k+n*konx <= kmax) bufin(ii) = Fp(i,j,k+n*konx)
            end do
          end do
        end do
      end do
  
      !$acc host_data use_device(workspace_0, workspace_1)
      call D_MPI_ALLTOALL(bufin, iony*jonx*konx, &
                          bufout, iony*jonx*konx, &
                          commrow, mpierr)
      !$acc end host_data
  
      !$acc parallel loop collapse(4) default(present) private(ii)
      do n = 0, nprocx-1
        do k = 1, konx
          do i = 1, iony
            do j = 1, jonx
              ii = j + (i-1)*jonx + (k-1)*iony*jonx + n*iony*jonx*konx
              if (j+n*jonx <= jtot) py(j+n*jonx,k,i) = bufout(ii)
            end do
          end do
        end do
      end do
    end if
  
  end subroutine transpose_a3inv

end module modtranspose