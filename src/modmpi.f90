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
use iso_fortran_env, only : real32,real64,int32
use iso_c_binding  , only : c_loc, c_f_pointer
use mpi_f08
implicit none
save
  type(MPI_COMM) :: commwrld, comm3d, commrow, commcol
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

  real     :: CPU_program    !end time
  real     :: CPU_program0   !start time

  character(8) :: cmyid
  character(3) :: cmyidx, cmyidy

  !------------------------------------------------------------------------------
  ! D_MPI_INTERFACE
  !
  ! Victor Azizi Escience Center 2021
  !
  ! purpose
  ! -------
  ! These interfaces determine the correct KIND for the mpi procedures from the
  ! KIND of the input arguments, this means that code that calls this functions 
  ! does not have to worry about changing the MPI_TYPE with a changing KIND
  ! 
  ! The argument list is the same as the corresponding MPI_ functions, but with
  ! the MPI_TYPE omitted
  !------------------------------------------------------------------------------
  interface D_MPI_ISEND
    procedure :: D_MPI_ISEND_REAL32
    procedure :: D_MPI_ISEND_REAL64
  end interface
  interface D_MPI_IRECV
    procedure :: D_MPI_IRECV_REAL32
    procedure :: D_MPI_IRECV_REAL64
  end interface
  interface D_MPI_RECV
    procedure :: D_MPI_RECV_REAL32
    procedure :: D_MPI_RECV_REAL64
  end interface
  interface D_MPI_BCAST
    procedure :: D_MPI_BCAST_REAL64_S
    procedure :: D_MPI_BCAST_INT32_S
    procedure :: D_MPI_BCAST_LOGICAL_S
    procedure :: D_MPI_BCAST_REAL32
    procedure :: D_MPI_BCAST_REAL64_R1
    procedure :: D_MPI_BCAST_REAL64_R2
    procedure :: D_MPI_BCAST_INT32_R1
    procedure :: D_MPI_BCAST_INT32_R2
    procedure :: D_MPI_BCAST_LOGICAL_R1
    procedure :: D_MPI_BCAST_STRING
  end interface
  interface D_MPI_ALLREDUCE
    procedure :: D_MPI_ALLREDUCE_REAL32_S
    procedure :: D_MPI_ALLREDUCE_REAL64_S
    procedure :: D_MPI_ALLREDUCE_INT32_S
    procedure :: D_MPI_ALLREDUCE_REAL32_R1
    procedure :: D_MPI_ALLREDUCE_REAL64_R1
    procedure :: D_MPI_ALLREDUCE_REAL64_R2
    procedure :: D_MPI_ALLREDUCE_REAL64_R3
    procedure :: D_MPI_ALLREDUCE_INT32_R1
    procedure :: D_MPI_ALLREDUCE_INT32_R2
    procedure :: D_MPI_ALLREDUCE_REAL32_IP
    procedure :: D_MPI_ALLREDUCE_REAL64_IP
  end interface
  interface D_MPI_ALLTOALL
    procedure :: D_MPI_ALLTOALL_REAL32
    procedure :: D_MPI_ALLTOALL_REAL64
  end interface
  interface D_MPI_REDUCE
    procedure :: D_MPI_REDUCE_REAL32
    procedure :: D_MPI_REDUCE_REAL64_R1
    procedure :: D_MPI_REDUCE_REAL64_R2
    procedure :: D_MPI_REDUCE_REAL64_R3
    procedure :: D_MPI_REDUCE_REAL32_IP
    procedure :: D_MPI_REDUCE_REAL64_IP_R1
    procedure :: D_MPI_REDUCE_REAL64_IP_R2
  end interface
  interface D_MPI_GATHER
    procedure :: D_MPI_GATHER_REAL32
    procedure :: D_MPI_GATHER_REAL64
  end interface

  interface excjs
    procedure :: excjs_real32
    procedure :: excjs_real64
  end interface
  interface slabsum
    procedure :: slabsum_real32
    procedure :: slabsum_real64
  end interface

contains
  ! Initializes the world communicator within dales. Optionally this communicator is passed from an external caller.
  ! TODO: Handle errors correctly.
  subroutine initmpicomm(comm)
    implicit none
    type(MPI_COMM), intent(in),optional  :: comm
    logical                              :: init
    integer                              :: ierr

    call MPI_INITIALIZED(init,mpierr)

    if(.not.init) then
        call MPI_INIT(mpierr)
    endif

    call MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN, ierr)
    if (ierr /= MPI_SUCCESS) then
      call abort
    end if

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

  subroutine excjs_real32(a,sx,ex,sy,ey,sz,ez,ih,jh)
  implicit none
  integer sx, ex, sy, ey, sz, ez, ih, jh
  real(real32) a(sx-ih:ex+ih, sy-jh:ey+jh, sz:ez)
  type(MPI_STATUS)  :: status
  integer :: xl, yl, zl
  type(MPI_REQUEST) :: reqn, reqs, reqe, reqw
  type(MPI_REQUEST) :: reqrn, reqrs, reqre, reqrw
  integer nssize, ewsize
  real(real32),allocatable, dimension(:) :: sendn,recvn &
                                          , sends,recvs &
                                          , sende,recve &
                                          , sendw,recvw

! Calulate buffer lengths
  xl = size(a,1)
  yl = size(a,2)
  zl = size(a,3)

!   Calculate buffer size
  nssize = xl*jh*zl
  ewsize = ih*yl*zl


  if(nprocy .gt. 1)then

    !   Allocate send / receive buffers
    allocate(sendn(nssize),sends(nssize),recvn(nssize),recvs(nssize))

    sendn = reshape(a(:,ey-jh+1:ey,:),(/nssize/))
    sends = reshape(a(:,sy:sy+jh-1,:),(/nssize/))

    !   Send north/south
    call D_MPI_ISEND(sendn, nssize, nbrnorth, 4, comm3d, reqn, mpierr)
    call D_MPI_ISEND(sends, nssize, nbrsouth, 5, comm3d, reqs, mpierr)

    !   Receive south/north
    call D_MPI_IRECV(recvs, nssize, nbrsouth, 4, comm3d, reqrs, mpierr)
    call D_MPI_IRECV(recvn, nssize, nbrnorth, 5, comm3d, reqrn, mpierr)

    ! Wait until data is received
    call MPI_WAIT(reqrs, status, mpierr)
    call MPI_WAIT(reqrn, status, mpierr)

    ! Write back buffers
    a(:,sy-jh:sy-1,:) = reshape(recvs,(/xl,jh,zl/))
    a(:,ey+1:ey+jh,:) = reshape(recvn,(/xl,jh,zl/))

  else

    ! Single processor, make sure the field is periodic
    a(:,sy-jh:sy-1,:) = a(:,ey-jh+1:ey,:)
    a(:,ey+1:ey+jh,:) = a(:,sy:sy+jh-1,:)
    
  endif

  if(nprocx .gt. 1)then

    !   Allocate send / receive buffers
    allocate(sende(ewsize),sendw(ewsize),recve(ewsize),recvw(ewsize))

    sende = reshape(a(ex-ih+1:ex,:,:),(/ewsize/))
    sendw = reshape(a(sx:sx+ih-1,:,:),(/ewsize/))

    !   Send east/west
    call D_MPI_ISEND(sende, ewsize, nbreast, 6, comm3d, reqe, mpierr)
    call D_MPI_ISEND(sendw, ewsize, nbrwest, 7, comm3d, reqw, mpierr)

    !   Receive west/east
    call D_MPI_IRECV(recvw, ewsize, nbrwest, 6, comm3d, reqrw, mpierr)
    call D_MPI_IRECV(recve, ewsize, nbreast, 7, comm3d, reqre, mpierr)

    ! Wait until data is received
    call MPI_WAIT(reqrw, status, mpierr)
    call MPI_WAIT(reqre, status, mpierr)

    ! Write back buffers
    a(sx-ih:sx-1,:,:) = reshape(recvw,(/ih,yl,zl/))
    a(ex+1:ex+ih,:,:) = reshape(recve,(/ih,yl,zl/))

  else

    ! Single processor, make sure the field is periodic
    a(sx-ih:sx-1,:,:) = a(ex-ih+1:ex,:,:)
    a(ex+1:ex+ih,:,:) = a(sx:sx+ih-1,:,:)

  endif

  if(nprocy.gt.1)then

    ! Make sure data is sent
    call MPI_WAIT(reqn, status, mpierr)
    call MPI_WAIT(reqs, status, mpierr)
  
    deallocate (sendn, sends)
    deallocate (recvn, recvs)

  endif

  if(nprocx.gt.1)then

    ! Make sure data is sent
    call MPI_WAIT(reqe, status, mpierr)
    call MPI_WAIT(reqw, status, mpierr)

    ! Deallocate buffers
    deallocate (sende, sendw)
    deallocate (recve, recvw)

  endif

  return
  end subroutine excjs_real32

  subroutine excjs_real64(a,sx,ex,sy,ey,sz,ez,ih,jh)
  implicit none
  integer sx, ex, sy, ey, sz, ez, ih, jh
  real(real64) a(sx-ih:ex+ih, sy-jh:ey+jh, sz:ez)
  type(MPI_STATUS)  :: status
  integer :: xl, yl, zl
  type(MPI_REQUEST) :: reqn, reqs, reqe, reqw
  type(MPI_REQUEST) :: reqrn, reqrs, reqre, reqrw
  integer nssize, ewsize
  real(real64),allocatable, dimension(:) :: sendn,recvn &
                                          , sends,recvs &
                                          , sende,recve &
                                          , sendw,recvw

! Calulate buffer lengths
  xl = size(a,1)
  yl = size(a,2)
  zl = size(a,3)

!   Calculate buffer size
  nssize = xl*jh*zl
  ewsize = ih*yl*zl


  if(nprocy .gt. 1)then

    !   Allocate send / receive buffers
    allocate(sendn(nssize),sends(nssize),recvn(nssize),recvs(nssize))

    sendn = reshape(a(:,ey-jh+1:ey,:),(/nssize/))
    sends = reshape(a(:,sy:sy+jh-1,:),(/nssize/))

    !   Send north/south
    call D_MPI_ISEND(sendn, nssize, nbrnorth, 4, comm3d, reqn, mpierr)
    call D_MPI_ISEND(sends, nssize, nbrsouth, 5, comm3d, reqs, mpierr)

    !   Receive south/north
    call D_MPI_IRECV(recvs, nssize, nbrsouth, 4, comm3d, reqrs, mpierr)
    call D_MPI_IRECV(recvn, nssize, nbrnorth, 5, comm3d, reqrn, mpierr)

    ! Wait until data is received
    call MPI_WAIT(reqrs, status, mpierr)
    call MPI_WAIT(reqrn, status, mpierr)


    ! Write back buffers
    a(:,sy-jh:sy-1,:) = reshape(recvs,(/xl,jh,zl/))
    a(:,ey+1:ey+jh,:) = reshape(recvn,(/xl,jh,zl/))

  else

    ! Single processor, make sure the field is periodic
    a(:,sy-jh:sy-1,:) = a(:,ey-jh+1:ey,:)
    a(:,ey+1:ey+jh,:) = a(:,sy:sy+jh-1,:)
    
  endif

  if(nprocx .gt. 1)then

    !   Allocate send / receive buffers
    allocate(sende(ewsize),sendw(ewsize),recve(ewsize),recvw(ewsize))

    sende = reshape(a(ex-ih+1:ex,:,:),(/ewsize/))
    sendw = reshape(a(sx:sx+ih-1,:,:),(/ewsize/))

    !   Send east/west
    call D_MPI_ISEND(sende, ewsize, nbreast, 6, comm3d, reqe, mpierr)
    call D_MPI_ISEND(sendw, ewsize, nbrwest, 7, comm3d, reqw, mpierr)

    !   Receive west/east
    call D_MPI_IRECV(recvw, ewsize, nbrwest, 6, comm3d, reqrw, mpierr)
    call D_MPI_IRECV(recve, ewsize, nbreast, 7, comm3d, reqre, mpierr)

    ! Wait until data is received
    call MPI_WAIT(reqrw, status, mpierr)
    call MPI_WAIT(reqre, status, mpierr)

    ! Write back buffers
    a(sx-ih:sx-1,:,:) = reshape(recvw,(/ih,yl,zl/))
    a(ex+1:ex+ih,:,:) = reshape(recve,(/ih,yl,zl/))

  else

    ! Single processor, make sure the field is periodic
    a(sx-ih:sx-1,:,:) = a(ex-ih+1:ex,:,:)
    a(ex+1:ex+ih,:,:) = a(sx:sx+ih-1,:,:)

  endif

  if(nprocy.gt.1)then

    ! Make sure data is sent
    call MPI_WAIT(reqn, status, mpierr)
    if (mpierr /= MPI_SUCCESS) call abort
    call MPI_WAIT(reqs, status, mpierr)
    if (mpierr /= MPI_SUCCESS) call abort
  
    deallocate (sendn, sends)
    deallocate (recvn, recvs)

  endif

  if(nprocx.gt.1)then

    ! Make sure data is sent
    call MPI_WAIT(reqe, status, mpierr)
    if (mpierr /= MPI_SUCCESS) call abort
    call MPI_WAIT(reqw, status, mpierr)
    if (mpierr /= MPI_SUCCESS) call abort

    ! Deallocate buffers
    deallocate (sende, sendw)
    deallocate (recve, recvw)

  endif

  return
  end subroutine excjs_real64

  subroutine slabsum_real32(aver,ks,kf,var,ib,ie,jb,je,kb,ke,ibs,ies,jbs,jes,kbs,kes)
    implicit none

    integer      :: ks,kf
    integer      :: ib,ie,jb,je,kb,ke,ibs,ies,jbs,jes,kbs,kes
    real(real32) :: aver(ks:kf)
    real(real32) :: var (ib:ie,jb:je,kb:ke)
    real(real32) :: averl(ks:kf)
    real(real32) :: avers(ks:kf)
    integer      :: k

    averl       = 0.
    avers       = 0.

    do k=kbs,kes
      averl(k) = sum(var(ibs:ies,jbs:jes,k))
    enddo

    call MPI_ALLREDUCE(averl, avers, kf-ks+1,  MPI_REAL4, &
                       MPI_SUM, comm3d,mpierr)

    aver = aver + avers

    return
  end subroutine slabsum_real32

  subroutine slabsum_real64(aver,ks,kf,var,ib,ie,jb,je,kb,ke,ibs,ies,jbs,jes,kbs,kes)
    implicit none

    integer      :: ks,kf
    integer      :: ib,ie,jb,je,kb,ke,ibs,ies,jbs,jes,kbs,kes
    real(real64) :: aver(ks:kf)
    real(real64) :: var (ib:ie,jb:je,kb:ke)
    real(real64) :: averl(ks:kf)
    real(real64) :: avers(ks:kf)
    integer      :: k

    averl       = 0.
    avers       = 0.

    do k=kbs,kes
      averl(k) = sum(var(ibs:ies,jbs:jes,k))
    enddo

    call MPI_ALLREDUCE(averl, avers, kf-ks+1,  MPI_REAL8, &
                       MPI_SUM, comm3d,mpierr)

    aver = aver + avers

    return
  end subroutine slabsum_real64

  subroutine mpi_get_time(val)
   real(real32), intent(out) :: val
 
   val = MPI_Wtime()
   call MPI_BCAST(val,1,MPI_REAL4,0,comm3d,mpierr)

  end subroutine mpi_get_time

  ! Gather a variable l(imax,jmax) along a row (ie. constant myidy)
  ! into              g(itot,jmax) at the processor with myix=0

  subroutine gatherrow(l,g,imax,jmax,itot)
    implicit none

    integer, intent(in) :: itot,imax,jmax
    real, intent(in)    :: l(imax,jmax)
    real, intent(out)   :: g(itot,jmax)

    integer      :: n,i,j, ii
    real :: sbuffer(imax * jmax)
    real :: rbuffer(itot * jmax)

    ii = 0
    do j=1,jmax
    do i=1,imax
       ii = ii + 1
       sbuffer(ii) = l(i,j)
    enddo
    enddo

    call D_MPI_GATHER(sbuffer,jmax*imax, &
                      rbuffer,jmax*imax, &
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

  !MPI interfaces instantations for the various types
  subroutine D_MPI_ISEND_REAL32(buf, count, dest, tag, comm, request, ierror)
    implicit none
    real(real32), target, contiguous, asynchronous, intent(inout) ::   buf(:)
    integer       ::   count, dest, tag, ierror
    type(MPI_COMM):: comm
    type(MPI_REQUEST) :: request
    ! This cast is necessary for the fujitsu fortran compiler, hopefully it
    ! will not be neccessary anymore sometime in the future (then please remove the
    ! cast and pass the buf directly!
    ! Furthermore, the target attribute of the buf argument should then be
    ! changed to asynchronous, when the subroutine is non-blocking (MPI_I prefix)
    ! To match the actual MPI_XXX interface of mpi_f08 better
    BYTE, pointer :: buf_b
    call c_f_pointer(c_loc(buf),buf_b)
    call MPI_ISEND(buf_b,count,MPI_REAL4,dest,tag,comm,request,ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ISEND_REAL32
  subroutine D_MPI_ISEND_REAL64(buf, count, dest, tag, comm, request, ierror)
    implicit none
    real(real64), contiguous, asynchronous, intent(inout) ::   buf(:)
    integer       ::   count, dest, tag, ierror
    type(MPI_COMM):: comm
    type(MPI_REQUEST) :: request
    call MPI_ISEND(buf,count,MPI_REAL8,dest,tag,comm,request,ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ISEND_REAL64

  subroutine D_MPI_IRECV_REAL32(buf, count, source, tag, comm, request, ierror)
    implicit none
    real(real32), target, asynchronous,contiguous, intent(inout)  ::   buf(:)
    integer        :: count, source, tag, ierror
    type(MPI_COMM) :: comm
    type(MPI_REQUEST) :: request
    BYTE, pointer :: buf_b
    call c_f_pointer(c_loc(buf),buf_b)
    call MPI_IRECV(buf_b,count,MPI_REAL4,source,tag,comm,request,ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_IRECV_REAL32
  subroutine D_MPI_IRECV_REAL64(buf, count, source, tag, comm, request, ierror)
    implicit none
    real(real64), contiguous, asynchronous, intent(inout)  ::   buf(:)
    integer        :: count, source, tag, ierror
    type(MPI_COMM) :: comm
    type(MPI_REQUEST) :: request
    call MPI_IRECV(buf,count,MPI_REAL8,source,tag,comm,request,ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_IRECV_REAL64
  
  subroutine D_MPI_RECV_REAL32(buf, count, source, tag, comm, status, ierror)
    implicit none
    real(real32), target, contiguous, intent(inout)  ::   buf(:)
    integer        :: count, source, tag, ierror
    type(MPI_COMM) :: comm
    type(MPI_STATUS) :: status
    BYTE, pointer :: buf_b
    call c_f_pointer(c_loc(buf),buf_b)
    call MPI_RECV(buf_b,count,MPI_REAL4,source,tag,comm,status,ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_RECV_REAL32
  subroutine D_MPI_RECV_REAL64(buf, count, source, tag, comm, status, ierror)
    implicit none
    real(real64), contiguous, intent(inout)  ::   buf(:)
    integer        :: count, source, tag, ierror
    type(MPI_COMM) :: comm
    type(MPI_STATUS) :: status
    call MPI_RECV(buf,count,MPI_REAL8,source,tag,comm,status,ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_RECV_REAL64

  subroutine D_MPI_BCAST_REAL64_S(buffer, count, root, comm, ierror)
    implicit none
    real(real64), intent(inout)  ::  buffer
    integer        :: count, root, ierror
    type(MPI_COMM) :: comm
    call MPI_BCAST(buffer, count, MPI_REAL8, root, comm, ierror)
  end subroutine D_MPI_BCAST_REAL64_S
  subroutine D_MPI_BCAST_INT32_S(buffer, count, root, comm, ierror)
    implicit none
    integer(int32), intent(inout) ::  buffer
    integer        :: count, root, ierror
    type(MPI_COMM) :: comm
    call MPI_BCAST(buffer, count, MPI_INTEGER4, root, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_BCAST_INT32_S
  subroutine D_MPI_BCAST_LOGICAL_S(buffer, count, root, comm, ierror)
    implicit none
    logical, intent(inout)        :: buffer
    integer        :: count, root, ierror
    type(MPI_COMM) :: comm
    call MPI_BCAST(buffer, count, MPI_LOGICAL, root, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_BCAST_LOGICAL_S

  subroutine D_MPI_BCAST_REAL32(buffer, count, root, comm, ierror)
    implicit none
    real(real32), target, contiguous, intent(inout)   ::  buffer(:)
    integer        :: count, root, ierror
    type(MPI_COMM) :: comm
    BYTE, pointer :: buf_b
    call c_f_pointer(c_loc(buffer),buf_b)
    call MPI_BCAST(buf_b, count, MPI_REAL4, root, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_BCAST_REAL32
  subroutine D_MPI_BCAST_REAL64_R1(buffer, count, root, comm, ierror)
    implicit none
    real(real64), contiguous, intent(inout)  ::  buffer(:)
    integer        :: count, root, ierror
    type(MPI_COMM) :: comm
    call MPI_BCAST(buffer, count, MPI_REAL8, root, comm, ierror)
  end subroutine D_MPI_BCAST_REAL64_R1
  subroutine D_MPI_BCAST_REAL64_R2(buffer, count, root, comm, ierror)
    implicit none
    real(real64), contiguous, intent(inout)  ::  buffer(:,:)
    integer        :: count, root, ierror
    type(MPI_COMM) :: comm
    call MPI_BCAST(buffer, count, MPI_REAL8, root, comm, ierror)
  end subroutine D_MPI_BCAST_REAL64_R2
  subroutine D_MPI_BCAST_INT32_R2(buffer, count, root, comm, ierror)
    implicit none
    integer(int32), contiguous, intent(inout) ::  buffer(:,:)
    integer        :: count, root, ierror
    type(MPI_COMM) :: comm
    call MPI_BCAST(buffer, count, MPI_INTEGER4, root, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_BCAST_INT32_R2
  subroutine D_MPI_BCAST_INT32_R1(buffer, count, root, comm, ierror)
    implicit none
    integer(int32), contiguous, intent(inout) ::  buffer(:)
    integer        :: count, root, ierror
    type(MPI_COMM) :: comm
    call MPI_BCAST(buffer, count, MPI_INTEGER4, root, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_BCAST_INT32_R1
  subroutine D_MPI_BCAST_LOGICAL_R1(buffer, count, root, comm, ierror)
    implicit none
    logical, contiguous, intent(inout)        :: buffer(:)
    integer        :: count, root, ierror
    type(MPI_COMM) :: comm
    call MPI_BCAST(buffer, count, MPI_LOGICAL, root, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_BCAST_LOGICAL_R1
  subroutine D_MPI_BCAST_STRING(buffer, count, root, comm, ierror)
    implicit none
    character(len = *) :: buffer
    integer        :: count, root, ierror
    type(MPI_COMM) :: comm
    call MPI_BCAST(buffer, count, MPI_CHARACTER, root, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_BCAST_STRING

  subroutine D_MPI_ALLREDUCE_REAL32_S(sendbuf, recvbuf, count, op, comm, ierror)
    implicit none
    real(real32), target, intent(inout)   :: sendbuf, recvbuf
    integer        :: count, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    BYTE, pointer :: bufs_b, bufr_b
    call c_f_pointer(c_loc(sendbuf),bufs_b)
    call c_f_pointer(c_loc(recvbuf),bufr_b)
    call MPI_ALLREDUCE(bufs_b, bufr_b, count, MPI_REAL4, op, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ALLREDUCE_REAL32_S
  subroutine D_MPI_ALLREDUCE_REAL64_S(sendbuf, recvbuf, count, op, comm, ierror)
    implicit none
    real(real64), intent(inout)   :: sendbuf, recvbuf
    integer        :: count, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    call MPI_ALLREDUCE(sendbuf, recvbuf, count, MPI_REAL8, op, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ALLREDUCE_REAL64_S
  subroutine D_MPI_ALLREDUCE_INT32_S(sendbuf, recvbuf, count, op, comm, ierror)
    implicit none
    integer(int32), intent(inout) :: sendbuf, recvbuf
    integer        :: count, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    call MPI_ALLREDUCE(sendbuf, recvbuf, count, MPI_INTEGER4, op, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ALLREDUCE_INT32_S
  subroutine D_MPI_ALLREDUCE_REAL32_R1(sendbuf, recvbuf, count, op, comm, ierror)
    implicit none
    real(real32), target, contiguous, intent(inout)   :: sendbuf(:), recvbuf(:)
    integer        :: count, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    BYTE, pointer :: bufs_b, bufr_b
    call c_f_pointer(c_loc(sendbuf),bufs_b)
    call c_f_pointer(c_loc(recvbuf),bufr_b)
    call MPI_ALLREDUCE(bufs_b, bufr_b, count, MPI_REAL4, op, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ALLREDUCE_REAL32_R1
  subroutine D_MPI_ALLREDUCE_REAL64_R1(sendbuf, recvbuf, count, op, comm, ierror)
    implicit none
    real(real64), contiguous, intent(inout)   :: sendbuf(:), recvbuf(:)
    integer        :: count, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    call MPI_ALLREDUCE(sendbuf, recvbuf, count, MPI_REAL8, op, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ALLREDUCE_REAL64_R1
  subroutine D_MPI_ALLREDUCE_REAL64_R2(sendbuf, recvbuf, count, op, comm, ierror)
    implicit none
    real(real64), contiguous, intent(inout)   :: sendbuf(:,:), recvbuf(:,:)
    integer        :: count, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    call MPI_ALLREDUCE(sendbuf, recvbuf, count, MPI_REAL8, op, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ALLREDUCE_REAL64_R2
  subroutine D_MPI_ALLREDUCE_REAL64_R3(sendbuf, recvbuf, count, op, comm, ierror)
    implicit none
    real(real64), contiguous, intent(inout)   :: sendbuf(:,:,:), recvbuf(:,:,:)
    integer        :: count, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    call MPI_ALLREDUCE(sendbuf, recvbuf, count, MPI_REAL8, op, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ALLREDUCE_REAL64_R3
  subroutine D_MPI_ALLREDUCE_INT32_R2(sendbuf, recvbuf, count, op, comm, ierror)
    implicit none
    integer(int32), contiguous, intent(inout) :: sendbuf(:,:), recvbuf(:,:)
    integer        :: count, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    call MPI_ALLREDUCE(sendbuf, recvbuf, count, MPI_INTEGER4, op, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ALLREDUCE_INT32_R2
  subroutine D_MPI_ALLREDUCE_INT32_R1(sendbuf, recvbuf, count, op, comm, ierror)
    implicit none
    integer(int32), contiguous, intent(inout) :: sendbuf(:), recvbuf(:)
    integer        :: count, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    call MPI_ALLREDUCE(sendbuf, recvbuf, count, MPI_INTEGER4, op, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ALLREDUCE_INT32_R1

  subroutine D_MPI_ALLREDUCE_REAL32_IP(recvbuf, count, op, comm, ierror)
    implicit none
    real(real32), target, contiguous, intent(inout)   :: recvbuf(:)
    integer        :: count, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    BYTE, pointer  :: bufr_b
    call c_f_pointer(c_loc(recvbuf),bufr_b)
    call MPI_ALLREDUCE(MPI_IN_PLACE, bufr_b, count, MPI_REAL4, op, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ALLREDUCE_REAL32_IP
  subroutine D_MPI_ALLREDUCE_REAL64_IP(recvbuf, count, op, comm, ierror)
    implicit none
    real(real64), contiguous, intent(inout)   :: recvbuf(:)
    integer        :: count, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    call MPI_ALLREDUCE(MPI_IN_PLACE, recvbuf, count, MPI_REAL8, op, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ALLREDUCE_REAL64_IP

  subroutine D_MPI_ALLTOALL_REAL32(sendbuf, sendcount, recvbuf, recvcount, comm, ierror)
    implicit none
    real(real32), target, contiguous, intent(inout)   :: sendbuf(:), recvbuf(:)
    integer        :: sendcount, recvcount, ierror
    type(MPI_COMM) :: comm
    BYTE, pointer :: bufs_b, bufr_b
    call c_f_pointer(c_loc(sendbuf),bufs_b)
    call c_f_pointer(c_loc(recvbuf),bufr_b)
    call MPI_ALLTOALL(bufs_b, sendcount, MPI_REAL4, bufr_b, recvcount, MPI_REAL4, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ALLTOALL_REAL32
  subroutine D_MPI_ALLTOALL_REAL64(sendbuf, sendcount, recvbuf, recvcount, comm, ierror)
    implicit none
    real(real64), contiguous, intent(inout)   :: sendbuf(:), recvbuf(:)
    integer        :: sendcount, recvcount, ierror
    type(MPI_COMM) :: comm
    call MPI_ALLTOALL(sendbuf, sendcount, MPI_REAL8, recvbuf, recvcount, MPI_REAL8, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ALLTOALL_REAL64

  subroutine D_MPI_REDUCE_REAL32(sendbuf, recvbuf, count, op, root, comm, ierror)
    implicit none
    real(real32), target, contiguous, intent(inout)   :: sendbuf(:), recvbuf(:)
    integer        :: count, root, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    BYTE, pointer :: bufs_b, bufr_b
    call c_f_pointer(c_loc(sendbuf),bufs_b)
    call c_f_pointer(c_loc(recvbuf),bufr_b)
    call MPI_REDUCE(bufs_b, bufr_b, count, MPI_REAL4, op, root, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_REDUCE_REAL32
  subroutine D_MPI_REDUCE_REAL64_R1(sendbuf, recvbuf, count, op, root, comm, ierror)
    implicit none
    real(real64), contiguous, intent(inout)   :: sendbuf(:), recvbuf(:)
    integer        :: count, root, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    call MPI_REDUCE(sendbuf, recvbuf, count, MPI_REAL8, op, root, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_REDUCE_REAL64_R1
  subroutine D_MPI_REDUCE_REAL64_R2(sendbuf, recvbuf, count, op, root, comm, ierror)
    implicit none
    real(real64), contiguous, intent(inout)   :: sendbuf(:,:), recvbuf(:,:)
    integer        :: count, root, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    call MPI_REDUCE(sendbuf, recvbuf, count, MPI_REAL8, op, root, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_REDUCE_REAL64_R2
  subroutine D_MPI_REDUCE_REAL64_R3(sendbuf, recvbuf, count, op, root, comm, ierror)
    implicit none
    real(real64), contiguous, intent(inout)   :: sendbuf(:,:,:), recvbuf(:,:,:)
    integer        :: count, root, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    call MPI_REDUCE(sendbuf, recvbuf, count, MPI_REAL8, op, root, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_REDUCE_REAL64_R3
  subroutine D_MPI_REDUCE_REAL32_IP(recvbuf, count, op, root, comm, ierror)
    implicit none
    real(real32), target, contiguous, intent(inout)   :: recvbuf(:)
    integer        :: count, root, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    BYTE, pointer :: bufr_b
    call c_f_pointer(c_loc(recvbuf),bufr_b)
    call MPI_REDUCE(MPI_IN_PLACE, bufr_b, count, MPI_REAL4, op, root, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_REDUCE_REAL32_IP
  subroutine D_MPI_REDUCE_REAL64_IP_R1(recvbuf, count, op, root, comm, ierror)
    implicit none
    real(real64), contiguous, intent(inout)  :: recvbuf(:)
    integer        :: count, root, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    call MPI_REDUCE(MPI_IN_PLACE, recvbuf, count, MPI_REAL8, op, root, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_REDUCE_REAL64_IP_R1
  subroutine D_MPI_REDUCE_REAL64_IP_R2(recvbuf, count, op, root, comm, ierror)
    implicit none
    real(real64), contiguous, intent(inout)  :: recvbuf(:,:)
    integer        :: count, root, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    call MPI_REDUCE(MPI_IN_PLACE, recvbuf, count, MPI_REAL8, op, root, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_REDUCE_REAL64_IP_R2

  subroutine D_MPI_GATHER_REAL32(sendbuf, sendcount, recvbuf, recvcount, root, comm, ierror)
    implicit none
    real(real32), target, contiguous, intent(inout)   :: sendbuf(:), recvbuf(:)
    integer        :: sendcount, recvcount, root, ierror
    type(MPI_COMM) :: comm
    BYTE, pointer :: bufs_b, bufr_b
    call c_f_pointer(c_loc(sendbuf),bufs_b)
    call c_f_pointer(c_loc(recvbuf),bufr_b)
    call MPI_GATHER( bufs_b, sendcount, MPI_REAL4 &
                   , bufr_b, recvcount, MPI_REAL4 &
                   , root, comm, ierror )
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_GATHER_REAL32
  subroutine D_MPI_GATHER_REAL64(sendbuf, sendcount, recvbuf, recvcount, root, comm, ierror)
    implicit none
    real(real64), contiguous, intent(inout)   :: sendbuf(:), recvbuf(:)
    integer        :: sendcount, recvcount, root, ierror
    type(MPI_COMM) :: comm
    call MPI_GATHER( sendbuf, sendcount, MPI_REAL8 &
                   , recvbuf, recvcount, MPI_REAL8 &
                   , root, comm, ierror )
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_GATHER_REAL64

end module
