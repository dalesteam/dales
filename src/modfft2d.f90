!> \file modfft2d.f90
!!  Perfoms a 2d Fourier transform on a NPROCX x NPROCY mesh of MPI nodes.

!>
!!  Perfoms the Fourier transform over x and y needed for solving the Poisson equation
!>
!!  \author Jisk Attema, Netherlands eScience Center
!  This file is part of DALES.
!
! DALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 2014 Netherlands eScience Center
!
module modfft2d

implicit none

save
  integer                             :: nkonx, nkony
  real, dimension(:),     allocatable :: winew, wjnew
  real, dimension(:,:,:), allocatable :: worka, workb
  real, dimension(:),     allocatable :: bufin, bufout

contains

  subroutine fft2dinit()

    use modmpi, only   : nprocx, nprocy
    use modglobal, only: itot, jtot, imax, jmax, kmax

    integer :: sz

! setup the matrix rotation.
! nkonx and nkony are the number of vertical (k) points per processor in the x and y direction.
! it is of course best if the kmax points can be distributed equally, but if not
! just let one row or column of processors do less points (controlled by 'ke' in the transpose functions)

    nkonx = kmax / nprocx
    if ( mod(kmax, nprocx) > 0 ) then
      nkonx = nkonx + 1
    endif

    nkony = kmax / nprocy
    if ( mod(kmax, nprocy) > 0 ) then
      nkony = nkony + 1
    endif

! Allocate communication buffers for the transpose functions
    sz = max( imax * jmax * nkonx * nprocx, &
              imax * jmax * nkony * nprocy  )

    allocate(bufin (sz))
    allocate(bufout(sz))

! Allocate temporary arrays to hold the rotated matrix
    allocate(worka(itot,jmax,nkonx))
    allocate(workb(jtot,imax,nkony))

! Prepare 1d FFT transforms
    allocate(winew(2*itot+15),wjnew(2*jtot+15))

    CALL rffti(itot,winew)
    CALL rffti(jtot,wjnew)

  end subroutine 


  subroutine fft2dexit()

    deallocate(bufin,bufout)
    deallocate(worka,workb)
    deallocate(winew,wjnew)

  end subroutine 


!
! Transpose functions:
!   Data is stored such that a whole column, 1:kmax, is on one node.
!   These funtions redistribute the matrix over the MPI nodes so that
!   a whole dimension (x or y) is on one processor. This enables us to do
!   a Fourier transform over that dimension.
!   
!
!
! Rotation a:  p(imax,jmax,kmax)                      ptrans(itot,jmax,nkonx)
!                  x     y  full                              full  y     x
!
!         /-------------/|                                /-------------/|
!        /../          / |                               /             / |
!  kmax  |------------|  |                               |------------|  |
!        |..|         |  |            ==>                |            |  |         
!        |..|         |  |                        nkonx  |____________|. |
!        |..|         | /    jmax                        |............|./    jmax
!   1    |--|--|--|---|/   1                        1    |------------|/   1
!                                                
!        1  imax                                         1            itot
!
!
!
! Rotation b:  p(imax,jmax,kmax)                      ptrans(itot,jmax,nkonx)
!                  x     y  full                              full  y     x
!
!         /-------------/|                                /-------------/|
!        /../          / |                               /             / |
!  kmax  |------------|  |                               |------------|  |
!        |..|         |  |            ==>                |..../       |  |         
!        |..|         |  |                        nkony  |___/________|. |
!        |..|         | /    jmax                        |...|        |./    jtot
!   1    |--|--|--|---|/   1                        1    |---|--------|/   1
!                                                
!        1  imax                                         1  imax          
!
!


  subroutine transpose_a(p,ptrans,ih,jh)
! data are on a single processor in the k-direction for p
! data are on a single processor in the i-direction for ptrans

    use mpi
    use modmpi, only : commrow, mpierr, my_real, nprocx
    use modglobal, only : i1,j1, itot, imax,jmax, kmax
    implicit none

    integer, intent(in)  :: ih,jh
    real, intent(in)  ::   p(2-ih:i1+ih,2-jh:j1+jh,kmax)
    real, intent(out) ::   ptrans(itot,jmax,nkonx)

    integer :: n, i,j,k, ii

    ii = 0
    do n=0,nprocx-1
    do k=n*nkonx + 1, (n+1)*nkonx
    do j=2,j1
    do i=2,i1
      ii = ii + 1
      if( k <= kmax ) then
        bufin(ii) = p(i,j,k)
      endif
    enddo
    enddo
    enddo
    enddo

    call MPI_ALLTOALL(bufin,   (imax*jmax*nkonx),my_real, &
                      bufout,  (imax*jmax*nkonx),my_real, &
                      commrow,mpierr)

    ii = 0
    do n=0,nprocx-1
    do k=1,nkonx
    do j=1,jmax
    do i=n*imax + 1, (n+1)*imax
        ii = ii + 1
        ptrans(i,j,k) = bufout(ii)
    enddo
    enddo
    enddo
    enddo

  end subroutine

  subroutine transpose_ainv(p,ptrans,ih,jh)
! data are on a single processor in the k-direction for p
! data are on a single processor in the i-direction for ptrans

    use mpi
    use modmpi, only : commrow, mpierr, my_real, nprocx
    use modglobal, only : i1,j1, itot, imax,jmax, kmax
    implicit none

    integer, intent(in)  :: ih,jh
    real, intent(inout)  :: p(2-ih:i1+ih,2-jh:j1+jh,kmax)
    real, intent(in)     :: ptrans(itot,jmax,nkonx)

    integer :: n, i,j,k, ii

    ii = 0
    do n=0,nprocx-1
    do k=1,nkonx
    do j=1,jmax
    do i=n*imax + 1, (n+1)*imax
      ii = ii + 1
      bufin(ii) = ptrans(i,j,k)
    enddo
    enddo
    enddo
    enddo

    call MPI_ALLTOALL(bufin,   (imax*jmax*nkonx),my_real, &
                      bufout,  (imax*jmax*nkonx),my_real, &
                      commrow,mpierr)

    ii = 0
    do n=0,nprocx-1
    do k=n*nkonx + 1, (n+1)*nkonx
    do j=2,j1
    do i=2,i1
      ii = ii + 1
      if( k <= kmax ) then
        p(i,j,k) = bufout(ii)
      endif
    enddo
    enddo
    enddo
    enddo

  end subroutine

  subroutine transpose_b(p,ptrans,ih,jh)
! data are on a single processor in the k-direction for p
! data are on a single processor in the i-direction for ptrans

    use mpi
    use modmpi, only : commcol, mpierr, nprocy, my_real
    use modglobal, only : i1,j1, jtot, imax,jmax, kmax
    implicit none

    integer, intent(in)  :: ih,jh
    real, intent(in)  ::   p(2-ih:i1+ih,2-jh:j1+jh,kmax)
    real, intent(out) ::   ptrans(jtot,imax,nkony)

    integer :: n, i,j,k, ii

    ii = 0
    do n=0,nprocy-1
    do k=n*nkony + 1, (n+1)*nkony
    do i=2,i1
    do j=2,j1
      ii = ii + 1
      if( k <= kmax ) then
         bufin(ii) = p(i,j,k)
      endif
    enddo
    enddo
    enddo
    enddo

    call MPI_ALLTOALL(bufin,   (imax*jmax*nkony),my_real, &
                      bufout,  (imax*jmax*nkony),my_real, &
                      commcol,mpierr)


    ii = 0
    do n=0,nprocy-1
    do k=1,nkony
    do i=1,imax
    do j=n*jmax + 1, (n+1)*jmax
      ii = ii + 1
      ptrans(j,i,k) = bufout(ii)
    enddo
    enddo
    enddo
    enddo

  end subroutine

  subroutine transpose_binv(p,ptrans,ih,jh)
! data are on a single processor in the k-direction for p
! data are on a single processor in the i-direction for ptrans

    use mpi
    use modmpi, only : commcol, mpierr, nprocy, my_real
    use modglobal, only : i1,j1, jtot, imax,jmax, kmax
    implicit none

    integer, intent(in)  :: ih,jh
    real, intent(inout)  ::   p(2-ih:i1+ih,2-jh:j1+jh,kmax)
    real, intent(in)     ::   ptrans(jtot,imax,nkony)

    integer :: n, i,j,k, ii

    ii = 0
    do n=0,nprocy-1
    do k=1,nkony
    do i=1,imax
    do j=n*jmax + 1, (n+1)*jmax
      ii = ii + 1
      bufin(ii) = ptrans(j,i,k) 
    enddo
    enddo
    enddo
    enddo

    call MPI_ALLTOALL(bufin,   (imax*jmax*nkony),my_real, &
                      bufout,  (imax*jmax*nkony),my_real, &
                      commcol,mpierr)

    ii = 0
    do n=0,nprocy-1
    do k=n*nkony + 1, (n+1)*nkony
    do i=2,i1
    do j=2,j1
      ii = ii + 1
      if( k <= kmax ) then
         p(i,j,k) = bufout(ii)
      endif
    enddo
    enddo
    enddo
    enddo

  end subroutine


  subroutine fft2df(p,ih,jh)

    use modglobal, only : imax, jmax, kmax, itot, jtot, ijtot, i1, j1
    use modmpi, only    : myidx,myidy

    integer, intent(in) :: ih,jh
    real, intent(inout) :: p(2-ih:i1+ih,2-jh:j1+jh,kmax)

    integer :: i,j,k, ke

! fft over i

    ke = min(kmax - myidx * nkonx, nkonx)

    call transpose_a(p, worka, ih, jh)
    do k=1,ke
    do j=1,jmax
      call rfftf(itot,worka(1,j,k),winew)
    enddo
    enddo
    call transpose_ainv(p, worka, ih, jh)

! fft over j

    ke = min(kmax - myidy * nkony, nkony)

    call transpose_b(p, workb, ih, jh)
    do k=1,ke
    do i=1,imax
      call rfftf(jtot,workb(1,i,k),wjnew)
    enddo
    enddo
    call transpose_binv(p, workb, ih, jh)            

    p(:,:,:) = p(:,:,:) / sqrt(ijtot)
  end subroutine


  subroutine fft2db(p,ih,jh)

    use modglobal, only : imax, jmax, kmax, itot, jtot, ijtot, i1, j1
    use modmpi, only    : myidx,myidy

    integer, intent(in) :: ih,jh
    real, intent(inout) :: p(2-ih:i1+ih,2-jh:j1+jh,kmax)

    integer :: i,j,k, ke

! inverse fft over j

    ke = min(kmax - myidy * nkony, nkony)

    call transpose_b(p, workb, ih, jh)
    do k=1,ke
    do i=1,imax
      call rfftb(jtot,workb(1,i,k),wjnew)
    enddo
    enddo
    call transpose_binv(p, workb, ih, jh)

! inverse fft over i

    ke = min(kmax - myidx * nkonx, nkonx)

    call transpose_a(p, worka, ih, jh)
    do k=1,ke
    do j=1,jmax
      call rfftb(itot,worka(1,j,k),winew)
    enddo
    enddo
    call transpose_ainv(p, worka, ih, jh)

    p(:,:,:) = p(:,:,:) / sqrt(ijtot)
  end subroutine

end module modfft2d
