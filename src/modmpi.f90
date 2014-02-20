!> \file modmpi.f90
!!  Layer to deal with the parallelization.

!>
!!  Layer to deal with the parallelization.
!>
!!  \author Matthieu Pourquie, TU Delft
!!  \par Revision list
!!  \todo Documentation
!!  \todo 2D/3D parallelization
!!  \todo Include interfaces for MPI_ALLREDUCE, MPI_ALLTOALL, MPI_BCAST,
!! MPI_SENDRECV to get rid of pure mpi calls in the code
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
!  Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!



module modmpi
use mpi
implicit none
save
  integer comm3d, commrow, commcol
  integer nbrnorth
  integer nbrsouth
  integer nbreast
  integer nbrwest
  integer myid
  integer myidx, myidy
  integer nprocs
  integer nprocx
  integer nprocy
  integer mpierr
  integer my_real
  real    CPU_program    !end time
  real    CPU_program0   !start time

   character(3) :: cmyid

contains
  subroutine initmpi
    implicit none
    integer dims(2)
    integer coords(2)
    logical periods(2)

    call MPI_INIT(mpierr)
    MY_REAL = MPI_DOUBLE_PRECISION
    call MPI_COMM_RANK( MPI_COMM_WORLD, myid, mpierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, mpierr )

! Specify the # procs in each direction.
! specifying a 0 means that MPI will try to find a useful # procs in
! the corresponding  direction,

    dims(1) = 0
    dims(2) = 0

! directions 1 and 2 are chosen periodic

    periods(1) = .true.
    periods(2) = .true.

! find suitable # procs in each direction

    call MPI_DIMS_CREATE( nprocs, 2, dims, mpierr )

    nprocx = dims(0)
    nprocy = dims(1)

! create the Cartesian communicator, denoted by the integer comm3d

    call MPI_CART_CREATE(MPI_COMM_WORLD, 2, dims, periods, .true., &
                        comm3d, mpierr )

! Get my processor number in this communicator

    call MPI_COMM_RANK( comm3d, myid, mpierr )

! when applying boundary conditions, we need to know which processors
! are neighbours in all 3 directions
! these are determined with the aid of the MPI routine MPI_CART_SHIFT,

    call MPI_CART_SHIFT( comm3d, 0,  1, nbrwest,  nbreast ,   mpierr )
    call MPI_CART_SHIFT( comm3d, 1,  1, nbrsouth, nbrnorth,   mpierr )

! Setup the row- and column- communicators
    call MPI_Cart_sub( comm3d, (/.TRUE.,.FALSE./), commrow, mpierr )
    call MPI_Cart_sub( comm3d, (/.FALSE.,.TRUE./), commcol, mpierr )

! Get the processors number in these communicators
    call MPI_COMM_RANK( commrow, myidx, mpierr )
    call MPI_COMM_RANK( commcol, myidy, mpierr )

! determine some useful MPI datatypes for sending/receiving data

     write(cmyid,'(i3.3)') myid

    if(myid==0)then
      CPU_program0 = MPI_Wtime()
    end if

    write(*,*)'nprocs = ', nprocs
    write(*,*)'myid = ', myid, coords(1), coords(2)

  end subroutine initmpi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine exitmpi
    implicit none


    if(myid==0)then
      CPU_program = MPI_Wtime() - CPU_program0
      write(6,*)'TOTAL CPU time = ', CPU_program
    end if

    call MPI_Comm_free( comm3d, mpierr )
    call MPI_FINALIZE(mpierr)
  end subroutine exitmpi

  subroutine barrou()
    implicit none
    call MPI_BARRIER(comm3d,mpierr)

  return
  end subroutine barrou



  subroutine excj( a, sx, ex, sy, ey, sz,ez)
    implicit none

  integer sx, ex, sy, ey, sz, ez
  real a(sx:ex, sy:ey, sz:ez)
  integer iiget, status(MPI_STATUS_SIZE)
  integer ii, i, j, k
  integer nssize, ewsize
  real,allocatable, dimension(:) :: nssend, nsrecv
  real,allocatable, dimension(:) :: ewsend, ewrecv
  nssize = (ex - sx + 1)*(ez - sz + 1)
  ewsize = (ey - sy + 1)*(ez - sz + 1)

  allocate( nssend(nssize),&
            nsrecv(nssize),&
            ewsend(ewsize),&
            ewrecv(ewsize))

!   communicate north/south
  
  if(nbrnorth/=MPI_PROC_NULL .AND. nbrsouth/=MPI_PROC_NULL)then
    ii = 0
    do k=sz,ez
    do i=sx,ex
      ii = ii + 1
      nssend(ii) = a(i,ey-1,k)
    enddo
    enddo

    call MPI_SENDRECV(  nssend,  nssize, MY_REAL, nbrnorth, 4, &
                        nsrecv,  iiget , MY_REAL, nbrsouth, 4, &
                        comm3d,  status, mpierr )

    ii = 0
    do k=sz,ez
    do i=sx,ex
      ii = ii + 1
      a(i,sy,k) = nsrecv(ii)

      nssend(ii) = a(i,sy+1,k)
    enddo
    enddo

    call MPI_SENDRECV(  nssend,  nssize, MY_REAL, nbrsouth, 5, &
                        nsrecv,  iiget , MY_REAL, nbrnorth, 5, &
                        comm3d,  status, mpierr )

    ii = 0
    do k=sz,ez
    do i=sx,ex
      ii = ii + 1
      a(i,ey,k) = nsrecv(ii)
    enddo
    enddo
  endif

!   communicate east/west

  if(nbreast/=MPI_PROC_NULL .AND. nbrwest/=MPI_PROC_NULL)then
    ii = 0
    do k=sz,ez
    do i=sy,ey
      ii = ii + 1
      ewsend(ii) = a(ex-1,i,k)
    enddo
    enddo

    call MPI_SENDRECV(  ewsend,  ewsize, MY_REAL, nbreast, 6, &
                        ewrecv,  iiget , MY_REAL, nbrwest, 6, &
                        comm3d,  status, mpierr )
    ii = 0
    do k=sz,ez
    do i=sy,ey
      ii = ii + 1
      a(sx,i,k) = ewrecv(ii)

      ewsend(ii) = a(sx+1,i,k)
    enddo
    enddo

    call MPI_SENDRECV(  ewsend,  ewsize, MY_REAL, nbrwest, 7, &
                        ewrecv,  iiget , MY_REAL, nbreast, 7, &
                        comm3d,  status, mpierr )
    ii = 0
    do k=sz,ez
    do i=sy,ey
      ii = ii + 1
      a(ex,i,k) = ewrecv(ii)
    enddo
    enddo
  endif


  deallocate (nssend,nsrecv,ewsend,ewrecv)

  return
  end subroutine excj

  subroutine excjs( a, sx, ex, sy, ey, sz,ez,ih,jh)
    implicit none
  integer sx, ex, sy, ey, sz, ez, ih, jh
  real a(sx-ih:ex+ih, sy-jh:ey+jh, sz:ez)
  integer status(MPI_STATUS_SIZE), iiget
  integer ii, i, j, k
  integer nssize, ewsize
  real,allocatable, dimension(:) :: nssend,nsrecv
  real,allocatable, dimension(:) :: ewsend,ewrecv
  nssize = jh*(ex - sx + 1 + 2*ih)*(ez - sz + 1)
  ewsize = ih*(ey - sy + 1 + 2*jh)*(ez - sz + 1)

  allocate( nssend(nssize),&
            nsrecv(nssize),&
            ewsend(ewsize),&
            ewrecv(ewsize))

!   Communicate north/south
  if(nbrnorth/=MPI_PROC_NULL .AND. nbrsouth/=MPI_PROC_NULL)then
    ii = 0
    do j=1,jh
    do k=sz,ez
    do i=sx-ih,ex+ih
      ii = ii + 1
      nssend(ii) = a(i,ey-j+1,k)
    enddo
    enddo
    enddo

    call MPI_SENDRECV(  nssend,  nssize, MY_REAL, nbrnorth, 4, &
                        nsrecv,  iiget , MY_REAL, nbrsouth, 4, &
                        comm3d,  status, mpierr )
    ii = 0
    do j=1,jh
    do k=sz,ez
    do i=sx-ih,ex+ih
      ii = ii + 1
      a(i,sy-j,k) = nsrecv(ii)

      nssend(ii) = a(i,sy+j-1,k)
    enddo
    enddo
    enddo

    call MPI_SENDRECV(  nssend,  nssize, MY_REAL, nbrsouth, 5, &
                        nsrecv,  iiget , MY_REAL, nbrnorth, 5, &
                        comm3d,  status, mpierr )

    ii = 0
    do j=1,jh
    do k=sz,ez
    do i=sx-ih,ex+ih
      ii = ii + 1
      a(i,ey+j,k) = nsrecv(ii)
    enddo
    enddo
    enddo
  endif

!   Communicate east/west
  if(nbreast/=MPI_PROC_NULL .AND. nbrwest/=MPI_PROC_NULL)then
    ii = 0
    do i=1,ih
    do k=sz,ez
    do j=sy-jh,ey+jh
      ii = ii + 1
      nssend(ii) = a(ex-i+1,j,k)
    enddo
    enddo
    enddo

    call MPI_SENDRECV(  ewsend,  ewsize, MY_REAL, nbreast, 6, &
                        ewrecv,  iiget , MY_REAL, nbrwest, 6, &
                        comm3d,  status, mpierr )
    ii = 0
    do i=1,ih
    do k=sz,ez
    do j=sy-jh,ey+jh
      ii = ii + 1
      a(sx-i,j,k) = nsrecv(ii)

      nssend(ii) = a(sx+i-1,j,k)
    enddo
    enddo
    enddo

    call MPI_SENDRECV(  ewsend,  ewsize, MY_REAL, nbrwest, 7, &
                        ewrecv,  iiget , MY_REAL, nbreast, 7, &
                        comm3d,  status, mpierr )

    ii = 0
    do i=1,jh
    do k=sz,ez
    do j=sy-jh,ey+jh
      ii = ii + 1
      a(ey+i,j,k) = nsrecv(ii)
    enddo
    enddo
    enddo
  endif

  deallocate (nssend,nsrecv,ewsend,ewrecv)

  return
  end subroutine excjs

  subroutine slabsum(aver,ks,kf,var,ib,ie,jb,je,kb,ke,ibs,ies,jbs,jes,kbs,kes)
    implicit none

    integer :: ks,kf
    integer :: ib,ie,jb,je,kb,ke,ibs,ies,jbs,jes,kbs,kes
    real    :: aver(ks:kf)
    real    :: var (ib:ie,jb:je,kb:ke)
    real    :: averl(ks:kf)
    real    :: avers(ks:kf)
    integer :: k

    averl       = 0.
    avers       = 0.

    do k=kbs,kes
      averl(k) = sum(var(ibs:ies,jbs:jes,k))
    enddo

    call MPI_ALLREDUCE(averl, avers, kf-ks+1,  MY_REAL, &
                          MPI_SUM, comm3d,mpierr)

    aver = aver + avers

    return
  end subroutine slabsum


end module
