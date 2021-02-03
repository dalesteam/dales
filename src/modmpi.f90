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
  integer  :: commwrld, comm3d, commrow, commcol
  logical  :: libmode !Library mode: skip finalize, assumed to be called externally
  integer  :: nbrnorth
  integer  :: nbrsouth
  integer  :: nbreast
  integer  :: nbrwest
  integer  :: myid
  integer  :: myidx, myidy
  integer  :: nprocs
  integer  :: nprocx = 1
  integer  :: nprocy = 0
  integer  :: mpierr
  integer  :: my_real = MPI_DOUBLE_PRECISION

  real     :: CPU_program    !end time
  real     :: CPU_program0   !start time

  character(8) :: cmyid
  character(3) :: cmyidx, cmyidy

contains

  ! Initializes the world communicator within dales. Optionally this communicator is passed from an external caller.
  ! TODO: Handle errors correctly.
  subroutine initmpicomm(comm)
    implicit none
    integer, intent(in),optional  :: comm
    logical                       :: init

    call MPI_INITIALIZED(init,mpierr)

    if(.not.init) then
        call MPI_INIT(mpierr)
    endif

    MY_REAL = MPI_DOUBLE_PRECISION

    if(present(comm)) then
        libmode=.true.
        if(comm==MPI_COMM_WORLD) then
            commwrld=comm
        else
            call MPI_COMM_DUP(comm,commwrld,mpierr)
        endif
    else
        libmode=.false.
        commwrld=MPI_COMM_WORLD
    endif

    call MPI_COMM_RANK( commwrld, myid, mpierr )
    call MPI_COMM_SIZE( commwrld, nprocs, mpierr )
  end subroutine initmpicomm

  
! This routine does the setup of the MPI mesh
! NPROCS
!        is the number of processors, this is set at run time, ie. mpirun -np 10
! NPROCX, NPROCY
!        are the number of processors in the x and y-direction. This set in the MPIOPT namelist.
!        A value of 0 lets MPI try to determine a suitable value
!        The old 'slab' parallelisation is equal to nprocx=1 and nprocy=0
!
! When using a large number of processors it is recommended to set NPROCX=NPROCY=0
! Otherwise, set NPROCX=1 and NPROCY=0 is probably faster (default)
!
! NOTE: the code is not symmetrical in NPROCX and NPROCY and NPROCX=0 NPROCY=1 will be
!       slower than the default.
!

  subroutine initmpi
    implicit none
    integer dims(2)
    logical periods(2)

! Specify the # procs in each direction.
! specifying a 0 means that MPI will try to find a useful # procs in
! the corresponding direction

    dims(1) = nprocx
    dims(2) = nprocy

! directions 1 and 2 are chosen periodic

    periods(1) = .true.
    periods(2) = .true.

! find suitable # procs in each direction

    call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, mpierr)
    call MPI_DIMS_CREATE( nprocs, 2, dims, mpierr )

    nprocx = dims(1)
    nprocy = dims(2)

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

! Get the processors ranks in these communicators
    call MPI_COMM_RANK( commrow, myidx, mpierr )
    call MPI_COMM_RANK( commcol, myidy, mpierr )

    if(myid==0)then
      CPU_program0 = MPI_Wtime()
      write(*,*) 'MPI mesh nprocx, nprocy: ', nprocx, nprocy
    end if

    !write(*,*)'myid, myidx, myidy, n, e, s, w = ', myid, myidx, myidy, nbrnorth, nbreast, nbrsouth, nbrwest
    write(cmyid,'(a,i3.3,a,i3.3)') 'x', myidx, 'y', myidy
    write(cmyidx,'(i3.3)') myidx
    write(cmyidy,'(i3.3)') myidy

  end subroutine initmpi


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine exitmpi
   implicit none
    logical :: mpifin

    if(myid==0)then
      CPU_program = MPI_Wtime() - CPU_program0
      write(6,*)'TOTAL wall time = ', CPU_program
    end if

    call MPI_Comm_free( comm3d, mpierr )

    if(commwrld/=MPI_COMM_WORLD .and. myid==0) then
        call MPI_COMM_FREE(commwrld,mpierr)
    endif

    call MPI_FINALIZED(mpifin,mpierr)
    if(.not.mpifin .and. .not.libmode) then
        call MPI_FINALIZE(mpierr)
    endif
  end subroutine exitmpi

  subroutine excjs(a,sx,ex,sy,ey,sz,ez,ih,jh)
  implicit none
  integer sx, ex, sy, ey, sz, ez, ih, jh
  real a(sx-ih:ex+ih, sy-jh:ey+jh, sz:ez)
  integer status(MPI_STATUS_SIZE)
  integer ii, i, j, k
  integer reqn, reqs, reqe, reqw
  integer nssize, ewsize
  real,allocatable, dimension(:) :: sendn,recvn
  real,allocatable, dimension(:) :: sends,recvs
  real,allocatable, dimension(:) :: sende,recve
  real,allocatable, dimension(:) :: sendw,recvw

!   Calculate buffer size
  nssize = jh*(ex - sx + 1 + 2*ih)*(ez - sz + 1)
  ewsize = ih*(ey - sy + 1 + 2*jh)*(ez - sz + 1)

!   Allocate send / receive buffers
  allocate(sendn(nssize),sends(nssize))
  allocate(sende(ewsize),sendw(ewsize))

  allocate(recvn(nssize),recvs(nssize))
  allocate(recve(ewsize),recvw(ewsize))

  if(nprocy .gt. 1)then
    !   Send north/south
    ii = 0
    do j=1,jh
    do k=sz,ez
    do i=sx-ih,ex+ih
      ii = ii + 1
      sendn(ii) = a(i,ey-j+1,k)
      sends(ii) = a(i,sy+j-1,k)
    enddo
    enddo
    enddo

    call MPI_ISEND(sendn, nssize, MY_REAL, nbrnorth, 4, comm3d, reqn, mpierr)
    call MPI_ISEND(sends, nssize, MY_REAL, nbrsouth, 5, comm3d, reqs, mpierr)

    !   Receive south/north
    call MPI_RECV(recvs, nssize, MY_REAL, nbrsouth, 4, comm3d, status, mpierr)
    call MPI_RECV(recvn, nssize, MY_REAL, nbrnorth, 5, comm3d, status, mpierr)

    ii = 0
    do j=1,jh
    do k=sz,ez
    do i=sx-ih,ex+ih
      ii = ii + 1
      a(i,sy-j,k) = recvs(ii)
      a(i,ey+j,k) = recvn(ii)
    enddo
    enddo
    enddo
  else
    ! Single processor, make sure the field is periodic
    do j=1,jh
    do k=sz,ez
    do i=sx-ih,ex+ih
      a(i,sy-j,k) = a(i,ey-j+1,k)
      a(i,ey+j,k) = a(i,sy+j-1,k)
    enddo
    enddo
    enddo
  endif

  if(nprocx .gt. 1)then
    !   Send east/west
    ii = 0
    do i=1,ih
    do k=sz,ez
    do j=sy-jh,ey+jh
      ii = ii + 1
      sende(ii) = a(ex-i+1,j,k)
      sendw(ii) = a(sx+i-1,j,k)
    enddo
    enddo
    enddo

    call MPI_ISEND(sende, ewsize, MY_REAL, nbreast, 6, comm3d, reqe, mpierr)
    call MPI_ISEND(sendw, ewsize, MY_REAL, nbrwest, 7, comm3d, reqw, mpierr)

    !   Receive west/east
    call MPI_RECV(recvw, ewsize, MY_REAL, nbrwest, 6, comm3d, status, mpierr)
    call MPI_RECV(recve, ewsize, MY_REAL, nbreast, 7, comm3d, status, mpierr)

    ii = 0
    do i=1,ih
    do k=sz,ez
    do j=sy-jh,ey+jh
      ii = ii + 1
      a(sx-i,j,k) = recvw(ii)
      a(ex+i,j,k) = recve(ii)
    enddo
    enddo
    enddo
  else
    ! Single processor, make sure the field is periodic
    do i=1,ih
    do k=sz,ez
    do j=sy-jh,ey+jh
      a(sx-i,j,k) = a(ex-i+1,j,k)
      a(ex+i,j,k) = a(sx+i-1,j,k)
    enddo
    enddo
    enddo
  endif

  if(nprocx.gt.1)then
    call MPI_WAIT(reqe, status, mpierr)
    call MPI_WAIT(reqw, status, mpierr)
  endif

  if(nprocy.gt.1)then
    call MPI_WAIT(reqn, status, mpierr)
    call MPI_WAIT(reqs, status, mpierr)
  endif

  deallocate (sendn, sends, sende, sendw)
  deallocate (recvn, recvs, recve, recvw)

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
  
  subroutine mpi_get_time(val)
   real, intent(out) :: val
 
   val = MPI_Wtime()
   call MPI_BCAST(val,1,MY_REAL   ,0,comm3d,mpierr)

  end subroutine mpi_get_time

  ! Gather a variable l(imax,jmax) along a row (ie. constant myidy)
  ! into              g(itot,jmax) at the processor with myix=0

  subroutine gatherrow(l,g,imax,jmax,itot)
    implicit none

    integer, intent(in) :: itot,imax,jmax
    real, intent(in)    :: l(imax,jmax)
    real, intent(out)   :: g(itot,jmax)

    integer  :: n,i,j, ii
    real     :: sbuffer(imax * jmax)
    real     :: rbuffer(itot * jmax)

    ii = 0
    do j=1,jmax
    do i=1,imax
       ii = ii + 1
       sbuffer(ii) = l(i,j)
    enddo
    enddo

    call MPI_GATHER(sbuffer,jmax*imax,MY_REAL, &
                    rbuffer,jmax*imax,MY_REAL, &
                    0, commrow,mpierr)

    if(myidx == 0) THEN
      ii = 0
      do n=0,nprocx-1
      do j=1,jmax
      do i=1 + n*imax,(n+1)*imax
          ii = ii + 1
          g(i,j) = rbuffer(ii)
      enddo
      enddo
      enddo
    endif

  end subroutine gatherrow

end module
