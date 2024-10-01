#define NX 4
#define NY 4
module modtranspose

  use modprecision, only: pois_r
  use modglobal,    only: itot, jtot, kmax, i1, j1, k1, ih, jh, imax, jmax
  use modmpi,       only: D_MPI_ALLTOALL, commrow, commcol, nprocs, &
    &                     nprocx, nprocy, mpierr, myidx, myidy
  use mpi_f08

  implicit none

  private

  public :: inittranspose
  public :: transpose_a1, transpose_a1inv
  public :: transpose_a2, transpose_a2inv
  public :: transpose_a3, transpose_a3inv

  real(pois_r), target, allocatable :: buffer(:) !< Buffer for transposes
  integer                           :: iony, jonx, konx
  integer                           :: iony_me, jonx_me, konx_me

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
  subroutine inittranspose

    integer :: sz

    ! setup the matrix transposes.
    ! konx is the number of vertical (k) points per processor in the x direction.
    ! iony and jonx are the number of i (j) points per processor in the y (x) direction.
    ! it is of course best if the points can be distributed equally, but if not
    ! just let one row or column of processors do less points

    konx = kmax / nprocx
    if ( mod(kmax, nprocx) > 0 ) then
      konx = konx + 1
    endif

    iony = itot / nprocy
    if ( mod(itot, nprocy) > 0 ) then
      iony = iony + 1
    endif

    jonx = jtot / nprocx
    if ( mod(jtot, nprocx) > 0 ) then
      jonx = jonx + 1
    endif

    ! how many data elements are in use in the current process?
    ! if the number of elements is not divisible by nprocx or nprocy,
    ! one process may have fewer elements, and some processes may have *no* elements
    konx_me = max(min(konx, kmax - konx*myidx), 0)
    iony_me = max(min(iony, itot - iony*myidy), 0)
    jonx_me = max(min(jonx, jtot - jonx*myidx), 0)    

    sz = max(imax * jmax * konx * nprocx, &
             iony * jmax * konx * nprocy, &
             iony * jonx * konx * nprocx)

    allocate(buffer(sz))

  end subroutine inittranspose

  subroutine transpose_a1(p, px)
  
    real(pois_r), pointer, intent(in)  :: p(:,:,:)
    real(pois_r), pointer, intent(out) :: px(:,:,:)
  
    integer :: i, j, k, n, ii
  
    if (nprocs == 1) then
      ! Already coalesced, no need to optimize further
      !$acc parallel loop collapse(3) default(present) private(ii)
      do k = 1, kmax
        do j = 1, jtot
          do i = 1, itot
            ii = i + (j-1)*itot + (k-1)*itot*jtot
            buffer(ii) = p(i+1,j+1,k)
          end do
        end do
      end do

      !$acc parallel loop collapse(3) default(present) private(ii)
      do k = 1, kmax
        do j = 1, jtot
          do i = 1, itot
            ii = i + (j-1)*itot + (k-1)*itot*jtot
            px(i,j,k) = buffer(ii)
          end do
        end do
      end do
    else
      !$acc parallel loop collapse(4) default(present) private(ii)
      do n = 0, nprocx - 1
        do k = 1, konx_me
          do j = 1, jmax
            do i = 1, imax
              if (k + n * konx_me <= kmax) then
                ii = i + (j-1)*imax + (k-1)*imax*jmax + n*imax*jmax*konx
                buffer(ii) = p(i+1,j+1,k+n*konx_me)
              end if
            end do
          end do
        end do
      end do
  
      ! CJ: updated to use in-place all-to-all so we only need to allocate one buffer.
      ! TODO: make generic interface for this one
#if POIS_PRECISION==32
      call MPI_ALLTOALL(MPI_IN_PLACE, 0, MPI_REAL4, buffer, imax*jmax*konx, MPI_REAL4, commrow, mpierr)
#else
      call MPI_ALLTOALL(MPI_IN_PLACE, 0, MPI_REAL8, buffer, imax*jmax*konx, MPI_REAL8, commrow, mpierr)
#endif
  
      !$acc parallel loop collapse(3) default(present) private(ii)
      do k = 1, konx
        do j = 1, jmax
          do i = 1, itot
            ii = i + (j-1)*itot + (k-1)*itot*jmax
            px(i,j,k) = buffer(ii)
          end do
        end do
      end do
    end if
  
  end subroutine transpose_a1
  
  subroutine transpose_a1inv(p, px)
  
    real(pois_r), pointer, intent(in)  :: px(:,:,:)
    real(pois_r), pointer, intent(out) :: p(:,:,:)
  
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
            p(i+1,j+1,k) = buffer(ii) 
          end do
        end do
      end do
    else
      !$acc parallel loop collapse(3) default(present) private(ii)
      do k = 1, konx
        do j = 1, jmax
          do i = 1, itot
            ii = i + (j-1)*itot + (k-1)*itot*jmax
            buffer(ii) = px(i,j,k)
          end do
        end do
      end do
  
      ! CJ: updated to use in-place all-to-all so we only need to allocate one buffer.
      ! TODO: make generic interface for this one
#if POIS_PRECISION==32
      call MPI_ALLTOALL(MPI_IN_PLACE, 0, MPI_REAL4, buffer, imax*jmax*konx, MPI_REAL4, commrow, mpierr)
#else
      call MPI_ALLTOALL(MPI_IN_PLACE, 0, MPI_REAL8, buffer, imax*jmax*konx, MPI_REAL8, commrow, mpierr)
#endif
  
      !$acc parallel loop collapse(4) default(present) private(ii)
      do n = 0, nprocx - 1
        do k = 1, konx_me
          do j = 1, jmax
            do i = 1, imax
              if (k + n * konx_me <= kmax) then
                ii = i + (j-1)*imax + (k-1)*imax*jmax + n*imax*jmax*konx
                p(i+1,j+1,k+n*konx_me) = buffer(ii)
              end if
            end do
          end do
        end do
      end do
    end if
  
  end subroutine transpose_a1inv
  
  subroutine transpose_a2(px, py)
  
    real(pois_r), pointer, intent(in)   :: px(:,:,:)
    real(pois_r), pointer, intent(out)  :: py(:,:,:)
 
    real(pois_r) :: shmem(NX, NY)
    integer      :: i, j, k, n, ii, iB, jB

    if (nprocs == 1) then
      !$acc parallel loop gang collapse(3) default(present) &
      !$acc& private(shmem) vector_length(NX*NY)
      do k = 1, kmax
        do j = 0, jtot - 1, NY
          do i = 0, itot, NX
            !$acc cache(shmem(1:NX,1:NY))
            !$acc loop vector collapse(2)
            do jB = 1, NY
              do iB = 1, NX
                if (i + iB <= itot .and. j + jB <= jtot) then
                  shmem(iB,jB) = px(i,j,k)
                end if
              end do
            end do
            !$acc loop vector collapse(2) private(ii)
            do iB = 1, NX
              do jB = 1, NY
                if (i + iB <= itot .and. j + jB <= jtot) then
                  ii = (j+jB) + ((i+iB)-1)*jtot + (k-1)*itot*jtot
                  buffer(ii) = shmem(jB,iB)
                end if
              end do
            end do
          end do
        end do
      end do

      !$acc parallel loop collapse(3) private(ii)
      do k = 1, kmax
        do i = 1, itot
          do j = 1, jtot
            ii = j + (i-1)*jtot + (k-1)*itot*jtot
            py(j,k,i) = buffer(ii)
          end do
        end do
      end do
    else
      !$acc parallel loop collapse(4) default(present) private(ii)
      do n = 0, nprocy - 1
        do k = 1, konx_me
          do j = 1, jmax
            do i = 1, iony_me
              if (i + n * iony_me <= itot) then
                ii = i + (j-1)*iony + (k-1)*iony*jmax + n*iony*jmax*konx
                buffer(ii) = px(i+n*iony_me,j,k)
              end if 
            end do
          end do
        end do
      end do
  
#if POIS_PRECISION==32
      call MPI_ALLTOALL(MPI_IN_PLACE, 0, MPI_REAL4, buffer, iony*jmax*konx, MPI_REAL4, commcol, mpierr)
#else
      call MPI_ALLTOALL(MPI_IN_PLACE, 0, MPI_REAL8, buffer, iony*jmax*konx, MPI_REAL8, commcol, mpierr)
#endif

      !$acc parallel loop gang collapse(4) default(present) &
      !$acc& private(shmem) vector_length(NX*NY)
      do n = 0, nprocy - 1
        do k = 1, konx_me
          do i = 0, iony_me - 1, NX
            do j = 0, jmax - 1, NY
              !$acc cache(shmem(1:NX,1:NY))
              !$acc loop vector collapse(2) private(ii)
              do jB = 1, NY
                do iB = 1, NX
                  if (i + iB <= iony_me .and. j + jB <= jmax) then
                    ii = (i+iB) + ((j+jB)-1)*iony + (k-1)*iony*jmax + n*iony*jmax*konx
                    shmem(iB,jB) = buffer(ii)
                  end if
                end do
              end do
              !$acc loop vector collapse(2) private(ii)
              do iB = 1, NX
                do jB = 1, NY
                  if (i + iB <= iony_me .and. j + jB <= jmax) then
                    py(j+jB+n*jmax,k,i+iB) = shmem(iB,jB)
                  end if
                end do
              end do
            end do
          end do
        end do
      end do
    end if
  
  end subroutine transpose_a2
  
  subroutine transpose_a2inv(px, py)

    real(pois_r), pointer, intent(in)  :: px(:,:,:)
    real(pois_r), pointer, intent(out) :: py(:,:,:)
  
    real(pois_r) :: shmem(NY,NX)
    integer      :: i, j, k, n, ii, iB, jB
  
    if (nprocs == 1) then
      !$acc parallel loop gang collapse(3) default(present) &
      !$acc& private(shmem(1:NY,1:NX)) vector_length(NX*NY)
      do k = 1, kmax
        do j = 0, jtot - 1, NY
          do i = 0, itot, NX
            !$acc cache(shmem(1:NY,1:NX))
            !$acc loop vector collapse(2)
            do iB = 1, NX
              do jB = 1, NY
                if (i + iB <= itot .and. j + jB <= jtot) then
                  shmem(jB,iB) = py(j,k,i)
                end if
              end do
            end do
            !$acc loop vector collapse(2) private(ii)
            do jB = 1, NY
              do iB = 1, NX
                if (i + iB <= itot .and. j + jB <= jtot) then
                  ii = j + (i-1)*jtot + (k-1)*itot*jtot
                  buffer(ii) = shmem(jB,iB)
                end if
              end do
            end do
          end do
        end do
      end do

      !$acc parallel loop collapse(3) private(ii)
      do k = 1, kmax
        do j = 1, jtot
          do i = 1, itot
            ii = j + (i-1)*jtot + (k-1)*itot*jtot
            px(i,j,k) = buffer(ii)
          end do
        end do
      end do
    else
      !$acc parallel loop gang collapse(4) default(present)
      !$acc& private(shmem(1:NY,1:NX)) vector_length(NX*NY)
      do n = 0, nprocy - 1
        do k = 1, konx_me
          do i = 0, iony_me - 1, NX
            do j = 0, jmax - 1, NY
              !$acc cache(shmem(1:NY,1:NX))
              !$acc loop vector collapse(2)
              do iB = 1, NX
                do jB = 1, NY
                  if (j + jB <= jmax .and. i + iB <= iony_me) then
                    shmem(jB,iB) = py(j+jB+n*jmax,k,i+iB)
                  end if
                end do
              end do
              !$acc loop vector collapse(2)
              do jB = 1, NY
                do iB = 1, NX
                  if (j + jB <= jmax .and. i + iB <= iony_me) then
                    ii = (i+iB) + ((j+jB)-1)*iony_me + (k-1)*iony_me*jmax + n*iony_me*jmax*konx_me
                    buffer(ii) = shmem(jB,iB)
                  end if
                end do
              end do
            end do
          end do
        end do
      end do

#if POIS_PRECISION==32
      call MPI_ALLTOALL(MPI_IN_PLACE, 0, MPI_REAL4, buffer, iony*jmax*konx, MPI_REAL4, commcol, mpierr)
#else
      call MPI_ALLTOALL(MPI_IN_PLACE, 0, MPI_REAL8, buffer, iony*jmax*konx, MPI_REAL8, commcol, mpierr)
#endif
  
      !$acc parallel loop collapse(4) default(present) private(ii)
      do n = 0, nprocy-1
        do k = 1, konx
          do j = 1, jmax
            do i = 1, iony
              if (i + n * iony_me <= itot) then
                ii = i + (j-1)*iony + (k-1)*iony*jmax + n*iony*jmax*konx
                px(i+n*iony,j,k) = buffer(ii)
              end if
            end do
          end do
        end do
      end do
    end if
  
  end subroutine transpose_a2inv
  
  subroutine transpose_a3(py, Fp)

    real(pois_r), pointer, intent(in)  :: py(:,:,:)
    real(pois_r), pointer, intent(out) :: Fp(:,:,:)
  
    real(pois_r) :: shmem(NY,NX)
    integer      :: i, j, k, n, ii, iB, jB
  
    if (nprocs == 1) then
      !$acc parallel loop gang collapse(3) default(present)
      !$acc& private(shmem(1:NY,1:NX)) vector_length(NX*NY)
      do k = 1, kmax
        do j = 1, jtot
          do i = 1, itot
            !$acc cache(shmem(1:NY,1:NX))
            !$acc loop vector collapse(2)
            do iB = 1, NX
              do jB = 1, NY
                if (i + iB <= itot .and. j + jB <= jtot) then
                  shmem(jB,iB) = py(j+jB,k,iB)
                end if
              end do
            end do
            !$acc loop vector collapse(2) private(ii)
            do jB = 1, NY
              do iB = 1, NX
                if (i + iB <= itot .and. j + jB <= jtot) then
                  ii = (i+iB) + ((j+jB)-1)*itot + (k-1)*itot*jtot
                  buffer(ii) = shmem(j+jB,i+iB)
                end if
              end do
            end do
          end do
        end do
      end do

      !$acc parallel loop collapse(3) default(present)
      do k = 1, kmax
        do j = 1, jtot
          do i = 1, itot
            ii = i + (j-1)*itot + (k-1)*itot*jtot
            Fp(i+1,j+1,k) = buffer(ii)
          end do
        end do
      end do
    else
      !$acc parallel loop collapse(4) default(present) private(ii)
      do n = 0, nprocx-1
        do k = 1, konx_me
          do i = 1, iony_me
            do j = 1, jonx_me
              if (j+n*jonx_me <= jtot) then
                ii = j + (i-1)*jonx_me + (k-1)*iony_me*jonx_me + n*iony_me*jonx_me*konx_me
                buffer(ii) = py(j+n*jonx_me,k,i)
              end if
            end do
          end do
        end do
      end do
  
#if POIS_PRECISION==32
      call MPI_ALLTOALL(MPI_IN_PLACE, 0, MPI_REAL4, buffer, iony_me*jonx_me*konx_me, MPI_REAL4, commrow, mpierr)
#else
      call MPI_ALLTOALL(MPI_IN_PLACE, 0, MPI_REAL8, buffer, iony_me*jonx_me*konx_me, MPI_REAL8, commrow, mpierr)
#endif
  
      !$acc parallel loop gang collapse(4) default(present) private(ii)
      !$acc& private(shmem(1:NY,1:NX)) vector_length(NX*NY)
      do n = 0, nprocx-1
        do k = 1, konx_me
          do j = 0, jonx_me - 1, NY
            do i = 0, iony_me - 1, NX
              !$acc cache(shmem(1:NY,1:NX))
              !$acc loop collapse(2) private(ii)
              do iB = 1, NX
                do jB = 1, NY
                  if (i + iB <= iony_me .and. j + jB <= jonx_me) then
                    ii = (j+jb) + ((i+iB)-1)*jonx_me + (k-1)*iony_me*jonx_me + n*iony_me*jonx_me*konx_me
                    shmem(jB,iB) = buffer(ii)
                  end if
                end do
              end do
              !$acc loop collapse(2)
              do jB = 1, NY
                do iB = 1, NX
                  if (i + iB <= iony_me .and. j + jB <= jonx_me .and. k + n * konx_me <= kmax) then
                    Fp(i+iB,j+jB,k+n*konx_me) = shmem(jB,iB)
                  end if
                end do
              end do
            end do
          end do
        end do
      end do
  
    end if
  
  end subroutine transpose_a3
  
  subroutine transpose_a3inv(py, Fp)
  
    real(pois_r), pointer, intent(in)  :: Fp(:,:,:)
    real(pois_r), pointer, intent(out) :: py(:,:,:)
  
    real(pois_r) :: shmem(NX,NY)
    integer      :: i, j, k, n, ii, iB, jB
  
    if (nprocs == 1) then
      !$acc parallel loop gang collapse(3) default(present)
      !$acc& private(shmem(1:NX,1:NY)) vector_length(NX*NY)
      do k = 1, kmax
        do j = 0, jtot - 1, NY
          do i = 0, itot - 1, NX
            !$acc cache(shmem(1:NX,1:NY))
            !$acc loop vector collapse(2)
            do jB = 1, NY
              do iB = 1, NX
                if (i + iB <= itot .and. j + jB <= jtot) then
                  shmem(iB,jB) = Fp(i+iB+1,j+jB+1,k)
                end if
              end do
            end do
            !$acc loop vector collapse(2) private(ii)
            do iB = 1, NX
              do jB = 1, NY
                if (i + iB <= itot .and. j + jB <= jtot) then
                  ii = (j+jB) + ((i+iB)-1)*jtot + (k-1)*itot*jtot
                  buffer(ii) = shmem(iB,jB)
                end if
              end do
            end do
          end do
        end do
      end do

      !$acc parallel loop collapse(3) default(present) private(ii)
      do k = 1, kmax
        do i = 1, itot
          do j = 1, jtot
            ii = j + (i-1)*jtot + (k-1)*itot*jtot
            py(j,k,i) = buffer(ii)
          end do
        end do
      end do
    else
      !$acc parallel loop gang collapse(4) default(present) private(ii)
      !$acc& private(shmem(1:NX,1:NY)) vector_length(NX*NY)
      do n = 0, nprocx-1
        do k = 1, konx_me
          do j = 0, jonx_me - 1, NY
            do i = 0, iony_me - 1, NX
              !$acc cache(shmem(1:NX,1:NY))
              !$acc loop vector collapse(2)
              do jB = 1, NY
                do iB = 1, NX
                  if (i + iB <= iony_me .and. j + jB <= jonx_me .and. k + n * konx_me <= kmax) then
                    shmem(iB,jB) = Fp(i+iB,j+jB,k+n*konx_me)
                  end if
                end do
              end do
              !$acc loop vector collapse(2) private(ii)
              do iB = 1, NX
                do jB = 1, NY
                  if (i + iB <= iony_me .and. j + jB <= jonx_me) then
                    ii = (j+jb) + ((i+iB)-1)*jonx_me + (k-1)*iony_me*jonx_me + n*iony_me*jonx_me*konx_me
                    buffer(ii) = shmem(iB,jB)
                  end if
                end do
              end do
            end do
          end do
        end do
      end do
  
#if POIS_PRECISION==32
      call MPI_ALLTOALL(MPI_IN_PLACE, 0, MPI_REAL4, buffer, iony_me*jonx_me*konx_me, MPI_REAL4, commrow, mpierr)
#else
      call MPI_ALLTOALL(MPI_IN_PLACE, 0, MPI_REAL8, buffer, iony_me*jonx_me*konx_me, MPI_REAL8, commrow, mpierr)
#endif
  
      !$acc parallel loop collapse(4) default(present) private(ii)
      do n = 0, nprocx - 1
        do k = 1, konx_me
          do i = 1, iony_me
            do j = 1, jonx_me
              if (j + n * jonx_me <= jtot) then
                ii = j + (i-1)*jonx_me + (k-1)*iony_me*jonx_me + n*iony_me*jonx_me*konx_me
                py(j+n*jonx_me,k,i) = buffer(ii)
              end if
            end do
          end do
        end do
      end do
    end if
  
  end subroutine transpose_a3inv

end module modtranspose
