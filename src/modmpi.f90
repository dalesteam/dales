!> \file modmpi.f90
!!  Layer to deal with the parallelization.

!>
!!  Layer to deal with the parallelization.
!>
!!  \author Matthieu Pourquie, TU Delft
!!  \author Jisk Attema
!!  \author Victor Azizi
!!  \author Fredrik Jansson
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
use modmpiinterface
#if defined(_OPENACC)
use openacc
#endif
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
  logical  :: periods(2) = .true.

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
  ! The implementations can be found in modmpiinterface
  !
  ! The argument list is the same as the corresponding MPI_ functions, but with
  ! the MPI_TYPE omitted
  !------------------------------------------------------------------------------
interface D_MPI_ISEND
    procedure :: D_MPI_ISEND_REAL32_R1
    procedure :: D_MPI_ISEND_REAL64_R1
    procedure :: D_MPI_ISEND_LOGICAL_R1
  end interface
  interface D_MPI_IRECV
    procedure :: D_MPI_IRECV_REAL32_R1
    procedure :: D_MPI_IRECV_REAL64_R1
    procedure :: D_MPI_IRECV_LOGICAL_R1
  end interface
  interface D_MPI_RECV
    procedure :: D_MPI_RECV_REAL32_R1
    procedure :: D_MPI_RECV_REAL64_R1
  end interface
  interface D_MPI_BCAST
    procedure :: D_MPI_BCAST_LOGICAL_S
    procedure :: D_MPI_BCAST_REAL64_S
    procedure :: D_MPI_BCAST_REAL32_S
    procedure :: D_MPI_BCAST_INT32_S
    procedure :: D_MPI_BCAST_REAL32_R1
    procedure :: D_MPI_BCAST_REAL32_R2
    procedure :: D_MPI_BCAST_REAL32_R3
    procedure :: D_MPI_BCAST_REAL64_R1
    procedure :: D_MPI_BCAST_REAL64_R2
    procedure :: D_MPI_BCAST_REAL64_R3
    procedure :: D_MPI_BCAST_INT32_R1
    procedure :: D_MPI_BCAST_INT32_R2
    procedure :: D_MPI_BCAST_INT32_R3
    procedure :: D_MPI_BCAST_LOGICAL_R1
    procedure :: D_MPI_BCAST_STRING
    procedure :: D_MPI_BCAST_STRING_R1
  end interface
  interface D_MPI_ALLREDUCE
    procedure :: D_MPI_ALLREDUCE_REAL32_S
    procedure :: D_MPI_ALLREDUCE_REAL64_S
    procedure :: D_MPI_ALLREDUCE_INT32_S
    procedure :: D_MPI_ALLREDUCE_REAL32_R1
    procedure :: D_MPI_ALLREDUCE_REAL32_R2
    procedure :: D_MPI_ALLREDUCE_REAL32_R3
    procedure :: D_MPI_ALLREDUCE_REAL64_R1
    procedure :: D_MPI_ALLREDUCE_REAL64_R2
    procedure :: D_MPI_ALLREDUCE_REAL64_R3
    procedure :: D_MPI_ALLREDUCE_INT32_R1
    procedure :: D_MPI_ALLREDUCE_INT32_R2
    procedure :: D_MPI_ALLREDUCE_REAL32_IP
    procedure :: D_MPI_ALLREDUCE_REAL64_IP
  end interface
  interface D_MPI_ALLTOALL
    procedure :: D_MPI_ALLTOALL_REAL32_R1
    procedure :: D_MPI_ALLTOALL_REAL64_R1
  end interface
  interface D_MPI_REDUCE
    procedure :: D_MPI_REDUCE_REAL32_R1
    procedure :: D_MPI_REDUCE_REAL32_R2
    procedure :: D_MPI_REDUCE_REAL32_R3
    procedure :: D_MPI_REDUCE_REAL64_R1
    procedure :: D_MPI_REDUCE_REAL64_R2
    procedure :: D_MPI_REDUCE_REAL64_R3
    procedure :: D_MPI_REDUCE_REAL32_IP_R1
    procedure :: D_MPI_REDUCE_REAL32_IP_R2
    procedure :: D_MPI_REDUCE_REAL64_IP_R1
    procedure :: D_MPI_REDUCE_REAL64_IP_R2
  end interface
  interface D_MPI_GATHER
    procedure :: D_MPI_GATHER_REAL32_R1
    procedure :: D_MPI_GATHER_REAL64_R1
  end interface
  interface excjs
    procedure :: excjs_real32
    procedure :: excjs_real64
    procedure :: excjs_logical
  end interface excjs

  interface openboundary_excjs
    procedure :: openboundary_excjs_real32
    procedure :: openboundary_excjs_real64
  end interface openboundary_excjs

  interface slabsum
    procedure :: slabsum_real32
    procedure :: slabsum_real64
  end interface
  interface slabsum_multi
    procedure :: slabsum_real32_5fields
    procedure :: slabsum_real32_4fields
    procedure :: slabsum_real32_3fields
    procedure :: slabsum_real32_2fields
    procedure :: slabsum_real64_5fields
    procedure :: slabsum_real64_4fields
    procedure :: slabsum_real64_3fields
    procedure :: slabsum_real64_2fields
  end interface


contains

! Subroutine for detecting and reporting namelist errors.
! Prints the last line read before failiure, as debugging help.
  subroutine checkmpierror (mpierr, location)
    implicit none
    integer, intent(in) :: mpierr
    character(*), intent(in) :: location
    integer len, err
    character(len = MPI_MAX_ERROR_STRING) :: str

    if (mpierr /= MPI_SUCCESS) then
       print *, 'MPI error', mpierr, 'in ', location
       ! look up the meaning of the error code
       call mpi_error_string(mpierr, str, len, err)
       if (err /= MPI_SUCCESS) then
          print *, 'Another error occurred when looking up the error code', err
          STOP
       endif
       print *, trim(str)
       STOP
    endif
  end subroutine checkmpierror

  ! Initializes the world communicator within dales. Optionally this communicator is passed from an external caller.
  subroutine initmpicomm(comm)
    implicit none
    type(MPI_COMM), intent(in),optional  :: comm
    logical                              :: init

    call MPI_INITIALIZED(init,mpierr)
    call checkmpierror(mpierr, 'MPI_INITIALIZED')

    if(.not.init) then
        call MPI_INIT(mpierr)
        call checkmpierror(mpierr, 'MPI_INIT')
     endif

    call MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN, mpierr)
    call checkmpierror(mpierr, 'MPI_Comm_set_errhandler')

    if(present(comm)) then
        libmode=.true.
        if(comm==MPI_COMM_WORLD) then
            commwrld=comm
        else
            call MPI_COMM_DUP(comm,commwrld,mpierr)
            call checkmpierror(mpierr, 'MPI_COMM_DUP')
        endif
    else
        libmode=.false.
        commwrld=MPI_COMM_WORLD
    endif

    call MPI_COMM_RANK( commwrld, myid, mpierr )
    call checkmpierror(mpierr, 'MPI_COMM_RANK')
    call MPI_COMM_SIZE( commwrld, nprocs, mpierr )
    call checkmpierror(mpierr, 'MPI_COMM_SIZE')
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
 ! if either nprocx = 0 or nprocy = 0 a value is computed automatically
 ! considering the total number of processors but not the itot,jtot grid size
    call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, mpierr)

    call checkmpierror(mpierr, 'MPI_COMM_SIZE')

    call MPI_DIMS_CREATE( nprocs, 2, dims, mpierr )
    if (mpierr /= MPI_SUCCESS) then
       if (myid == 0) then
          print *, 'MPI grid setup failed. '
          print *, '  nprocx', nprocx
          print *, '  nprocy', nprocy
          print *, '  nprocs', nprocs
          print *, 'nprocx * nprocy = nprocs is required but could not be achieved.'
       endif
    endif
    call checkmpierror(mpierr, 'MPI_DIMS_CREATE')

    nprocx = dims(1)
    nprocy = dims(2)

! create the Cartesian communicator, denoted by the integer comm3d

    call MPI_CART_CREATE(MPI_COMM_WORLD, 2, dims, periods, .true., &
                         comm3d, mpierr )
    call checkmpierror(mpierr, 'MPI_CART_CREATE')

! Get my processor number in this communicator

    call MPI_COMM_RANK( comm3d, myid, mpierr )
    call checkmpierror(mpierr, 'MPI_COMM_RANK')

! when applying boundary conditions, we need to know which processors
! are neighbours in all 3 directions
! these are determined with the aid of the MPI routine MPI_CART_SHIFT,

    call MPI_CART_SHIFT( comm3d, 0,  1, nbrwest,  nbreast ,   mpierr )
    call checkmpierror(mpierr, 'MPI_CART_SHIFT')
    call MPI_CART_SHIFT( comm3d, 1,  1, nbrsouth, nbrnorth,   mpierr )
    call checkmpierror(mpierr, 'MPI_CART_SHIFT')

! Setup the row- and column- communicators
    call MPI_Cart_sub( comm3d, (/.TRUE.,.FALSE./), commrow, mpierr )
    call checkmpierror(mpierr, 'MPI_Cart_sub')
    call MPI_Cart_sub( comm3d, (/.FALSE.,.TRUE./), commcol, mpierr )
    call checkmpierror(mpierr, 'MPI_Cart_sub')

! Get the processors ranks in these communicators
    call MPI_COMM_RANK( commrow, myidx, mpierr )
    call checkmpierror(mpierr, 'MPI_COMM_RANK')
    call MPI_COMM_RANK( commcol, myidy, mpierr )
    call checkmpierror(mpierr, 'MPI_COMM_RANK')

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
    call checkmpierror(mpierr, 'MPI_Comm_free')

    if(commwrld/=MPI_COMM_WORLD .and. myid==0) then
        call MPI_COMM_FREE(commwrld,mpierr)
        call checkmpierror(mpierr, 'MPI_COMM_FREE')
    endif

    call MPI_FINALIZED(mpifin,mpierr)
    call checkmpierror(mpierr, 'MPI_FINALIZED')
    if(.not.mpifin .and. .not.libmode) then
        call MPI_FINALIZE(mpierr)
        call checkmpierror(mpierr, 'MPI_FINALIZE')
    endif
  end subroutine exitmpi

  subroutine excjs_real32(a,sx,ex,sy,ey,sz,ez,ih,jh)
  implicit none
  integer :: sx, ex, sy, ey, sz, ez, ih, jh
  real(real32) :: a(sx-ih:ex+ih, sy-jh:ey+jh, sz:ez)
  type(MPI_STATUS) :: status
  integer :: xl, yl, zl
  type(MPI_REQUEST) :: reqn, reqs, reqe, reqw
  type(MPI_REQUEST) :: reqrn, reqrs, reqre, reqrw
  integer :: nssize, ewsize
  integer :: i, j, k, ii
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
    !$acc enter data copyin(sendn, sends, recvn, recvs)

    !$acc parallel loop collapse(3) default(present) private(ii)
    do k = 1, zl
      do j = 1, jh
        do i = 1, xl
          ii = i + (j-1)*xl + (k-1)*xl*jh
          sendn(ii) = a(sx-ih+i-1,ey-jh+j,k)
          sends(ii) = a(sx-ih+i-1,sy+j-1,k)
        end do
      end do
    end do

    !   Send north/south
    call D_MPI_ISEND(sendn, nssize, nbrnorth, 4, comm3d, reqn, mpierr, ondevice=.true.)
    call D_MPI_ISEND(sends, nssize, nbrsouth, 5, comm3d, reqs, mpierr, ondevice=.true.)

    !   Receive south/north
    call D_MPI_IRECV(recvs, nssize, nbrsouth, 4, comm3d, reqrs, mpierr, ondevice=.true.)
    call D_MPI_IRECV(recvn, nssize, nbrnorth, 5, comm3d, reqrn, mpierr, ondevice=.true.)

    ! Wait until data is received
    call MPI_WAIT(reqrs, status, mpierr)
    call MPI_WAIT(reqrn, status, mpierr)

    ! Write back buffers
    !$acc parallel loop collapse(3) default(present) private(ii)
    do k = 1, zl
      do j = 1, jh
        do i = 1, xl
          ii = i + (j-1)*xl + (k-1)*xl*jh
          a(sx-ih+i-1,ey+j,k) = recvn(ii)
          a(sx-ih+i-1,sy-jh+(j-1),k) = recvs(ii)
        end do
      end do
    end do

  else

    ! Single processor, make sure the field is periodic
    !$acc kernels default(present) async
    a(:,sy-jh:sy-1,:) = a(:,ey-jh+1:ey,:)
    a(:,ey+1:ey+jh,:) = a(:,sy:sy+jh-1,:)
    !$acc end kernels

  endif

  if(nprocx .gt. 1)then

    !   Allocate send / receive buffers
    allocate(sende(ewsize),sendw(ewsize),recve(ewsize),recvw(ewsize))
    !$acc enter data copyin(sende, sendw, recve, recvw)

    !$acc parallel loop collapse(3) default(present) private(ii)
    do k = 1, zl
      do j = 1, yl
        do i = 1, ih
          ii = i + (j-1)*ih + (k-1)*ih*yl
          sende(ii) = a(ex-ih+i,sy-jh+j-1,k)
          sendw(ii) = a(sx+i-1,sy-jh+j-1,k)
        end do
      end do
    end do

    !   Send east/west
    call D_MPI_ISEND(sende, ewsize, nbreast, 6, comm3d, reqe, mpierr, ondevice=.true.)
    call D_MPI_ISEND(sendw, ewsize, nbrwest, 7, comm3d, reqw, mpierr, ondevice=.true.)

    !   Receive west/east
    call D_MPI_IRECV(recvw, ewsize, nbrwest, 6, comm3d, reqrw, mpierr, ondevice=.true.)
    call D_MPI_IRECV(recve, ewsize, nbreast, 7, comm3d, reqre, mpierr, ondevice=.true.)

    ! Wait until data is received
    call MPI_WAIT(reqrw, status, mpierr)
    call MPI_WAIT(reqre, status, mpierr)

    ! Write back buffers
    !$acc parallel loop collapse(3) default(present) private(ii)
    do k = 1, zl
      do j = 1, yl
        do i = 1, ih
          ii = i + (j-1)*ih + (k-1)*ih*yl
          a(sx-jh+(i-1),sy-jh+j-1,k) = recvw(ii)
          a(ex+i,sy-jh+j-1,k) = recve(ii)
        end do
      end do
    end do

  else

    ! Single processor, make sure the field is periodic
    !$acc kernels default(present) async
    a(sx-ih:sx-1,:,:) = a(ex-ih+1:ex,:,:)
    a(ex+1:ex+ih,:,:) = a(sx:sx+ih-1,:,:)
    !$acc end kernels

  endif

  if(nprocy.gt.1)then

    ! Make sure data is sent
    call MPI_WAIT(reqn, status, mpierr)
    call MPI_WAIT(reqs, status, mpierr)

    !$acc exit data delete(sendn, sends, recvn, recvs)
    deallocate (sendn, sends)
    deallocate (recvn, recvs)

  endif

  if(nprocx.gt.1)then

    ! Make sure data is sent
    call MPI_WAIT(reqe, status, mpierr)
    call MPI_WAIT(reqw, status, mpierr)

    ! Deallocate buffers
    !$acc exit data delete(sende, sendw, recve, recvw)
    deallocate (sende, sendw)
    deallocate (recve, recvw)

  endif

  !$acc wait

  end subroutine excjs_real32

  subroutine openboundary_excjs_real32(a,sx,ex,sy,ey,sz,ez,ih,jh,switch)
  implicit none
  integer sx, ex, sy, ey, sz, ez, ih, jh
  real(real32) a(sx-ih:ex+ih, sy-jh:ey+jh, sz:ez)
  logical, dimension(4), intent(in) :: switch ! flags for which boundaries need exchange
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

    if(switch(4)) sendn = reshape(a(:,ey-jh+1:ey,:),(/nssize/))
    if(switch(3)) sends = reshape(a(:,sy:sy+jh-1,:),(/nssize/))

    !   Send north/south
    if(switch(4)) call D_MPI_ISEND(sendn, nssize, nbrnorth, 4, comm3d, reqn, mpierr)
    if(switch(3)) call D_MPI_ISEND(sends, nssize, nbrsouth, 5, comm3d, reqs, mpierr)

    !   Receive south/north
    if(switch(3)) call D_MPI_IRECV(recvs, nssize, nbrsouth, 4, comm3d, reqrs, mpierr)
    if(switch(4)) call D_MPI_IRECV(recvn, nssize, nbrnorth, 5, comm3d, reqrn, mpierr)

    ! Wait until data is received
    if(switch(3)) call MPI_WAIT(reqrs, status, mpierr)
    if(switch(4)) call MPI_WAIT(reqrn, status, mpierr)

    ! Write back buffers
    if(switch(3)) a(:,sy-jh:sy-1,:) = reshape(recvs,(/xl,jh,zl/))
    if(switch(4)) a(:,ey+1:ey+jh,:) = reshape(recvn,(/xl,jh,zl/))

  else

    ! Single processor, make sure the field is periodic
    if(switch(3)) a(:,sy-jh:sy-1,:) = a(:,ey-jh+1:ey,:)
    if(switch(3)) a(:,ey+1:ey+jh,:) = a(:,sy:sy+jh-1,:)

  endif

  if(nprocx .gt. 1)then

    !   Allocate send / receive buffers
    allocate(sende(ewsize),sendw(ewsize),recve(ewsize),recvw(ewsize))

    if(switch(2)) sende = reshape(a(ex-ih+1:ex,:,:),(/ewsize/))
    if(switch(1)) sendw = reshape(a(sx:sx+ih-1,:,:),(/ewsize/))

    !   Send east/west
    if(switch(2)) call D_MPI_ISEND(sende, ewsize, nbreast, 6, comm3d, reqe, mpierr)
    if(switch(1)) call D_MPI_ISEND(sendw, ewsize, nbrwest, 7, comm3d, reqw, mpierr)

    !   Receive west/east
    if(switch(1)) call D_MPI_IRECV(recvw, ewsize, nbrwest, 6, comm3d, reqrw, mpierr)
    if(switch(2)) call D_MPI_IRECV(recve, ewsize, nbreast, 7, comm3d, reqre, mpierr)

    ! Wait until data is received
    if(switch(1)) call MPI_WAIT(reqrw, status, mpierr)
    if(switch(2)) call MPI_WAIT(reqre, status, mpierr)

    ! Write back buffers
    if(switch(1)) a(sx-ih:sx-1,:,:) = reshape(recvw,(/ih,yl,zl/))
    if(switch(2)) a(ex+1:ex+ih,:,:) = reshape(recve,(/ih,yl,zl/))

  else

    ! Single processor, make sure the field is periodic
    if(switch(1)) a(sx-ih:sx-1,:,:) = a(ex-ih+1:ex,:,:)
    if(switch(1)) a(ex+1:ex+ih,:,:) = a(sx:sx+ih-1,:,:)

  endif

  if(nprocy.gt.1)then

    ! Make sure data is sent
    if(switch(4)) call MPI_WAIT(reqn, status, mpierr)
    if(switch(3)) call MPI_WAIT(reqs, status, mpierr)

    deallocate (sendn, sends)
    deallocate (recvn, recvs)

  endif

  if(nprocx.gt.1)then

    ! Make sure data is sent
    if(switch(2)) call MPI_WAIT(reqe, status, mpierr)
    if(switch(1)) call MPI_WAIT(reqw, status, mpierr)

    ! Deallocate buffers
    deallocate (sende, sendw)
    deallocate (recve, recvw)

  endif
  end subroutine openboundary_excjs_real32


  subroutine openboundary_excjs_real64(a,sx,ex,sy,ey,sz,ez,ih,jh,switch)
  implicit none
  integer sx, ex, sy, ey, sz, ez, ih, jh
  real(real64) a(sx-ih:ex+ih, sy-jh:ey+jh, sz:ez)
  logical, dimension(4), intent(in) :: switch ! flags for which boundaries need exchange
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

    if(switch(4)) sendn = reshape(a(:,ey-jh+1:ey,:),(/nssize/))
    if(switch(3)) sends = reshape(a(:,sy:sy+jh-1,:),(/nssize/))

    !   Send north/south
    if(switch(4)) call D_MPI_ISEND(sendn, nssize, nbrnorth, 4, comm3d, reqn, mpierr)
    if(switch(3)) call D_MPI_ISEND(sends, nssize, nbrsouth, 5, comm3d, reqs, mpierr)

    !   Receive south/north
    if(switch(3)) call D_MPI_IRECV(recvs, nssize, nbrsouth, 4, comm3d, reqrs, mpierr)
    if(switch(4)) call D_MPI_IRECV(recvn, nssize, nbrnorth, 5, comm3d, reqrn, mpierr)

    ! Wait until data is received
    if(switch(3)) call MPI_WAIT(reqrs, status, mpierr)
    if(switch(4)) call MPI_WAIT(reqrn, status, mpierr)

    ! Write back buffers
    if(switch(3)) a(:,sy-jh:sy-1,:) = reshape(recvs,(/xl,jh,zl/))
    if(switch(4)) a(:,ey+1:ey+jh,:) = reshape(recvn,(/xl,jh,zl/))

  else

    ! Single processor, make sure the field is periodic
    if(switch(3)) a(:,sy-jh:sy-1,:) = a(:,ey-jh+1:ey,:)
    if(switch(3)) a(:,ey+1:ey+jh,:) = a(:,sy:sy+jh-1,:)

  endif

  if(nprocx .gt. 1)then

    !   Allocate send / receive buffers
    allocate(sende(ewsize),sendw(ewsize),recve(ewsize),recvw(ewsize))

    if(switch(2)) sende = reshape(a(ex-ih+1:ex,:,:),(/ewsize/))
    if(switch(1)) sendw = reshape(a(sx:sx+ih-1,:,:),(/ewsize/))

    !   Send east/west
    if(switch(2)) call D_MPI_ISEND(sende, ewsize, nbreast, 6, comm3d, reqe, mpierr)
    if(switch(1)) call D_MPI_ISEND(sendw, ewsize, nbrwest, 7, comm3d, reqw, mpierr)

    !   Receive west/east
    if(switch(1)) call D_MPI_IRECV(recvw, ewsize, nbrwest, 6, comm3d, reqrw, mpierr)
    if(switch(2)) call D_MPI_IRECV(recve, ewsize, nbreast, 7, comm3d, reqre, mpierr)

    ! Wait until data is received
    if(switch(1)) call MPI_WAIT(reqrw, status, mpierr)
    if(switch(2)) call MPI_WAIT(reqre, status, mpierr)

    ! Write back buffers
    if(switch(1)) a(sx-ih:sx-1,:,:) = reshape(recvw,(/ih,yl,zl/))
    if(switch(2)) a(ex+1:ex+ih,:,:) = reshape(recve,(/ih,yl,zl/))

  else

    ! Single processor, make sure the field is periodic
    if(switch(1)) a(sx-ih:sx-1,:,:) = a(ex-ih+1:ex,:,:)
    if(switch(1)) a(ex+1:ex+ih,:,:) = a(sx:sx+ih-1,:,:)

  endif

  if(nprocy.gt.1)then

    ! Make sure data is sent
    if(switch(4)) call MPI_WAIT(reqn, status, mpierr)
    if(switch(3)) call MPI_WAIT(reqs, status, mpierr)

    deallocate (sendn, sends)
    deallocate (recvn, recvs)

  endif

  if(nprocx.gt.1)then

    ! Make sure data is sent
    if(switch(2)) call MPI_WAIT(reqe, status, mpierr)
    if(switch(1)) call MPI_WAIT(reqw, status, mpierr)

    ! Deallocate buffers
    deallocate (sende, sendw)
    deallocate (recve, recvw)

  endif
  end subroutine openboundary_excjs_real64


  subroutine excjs_real64(a,sx,ex,sy,ey,sz,ez,ih,jh)
  implicit none
  integer :: sx, ex, sy, ey, sz, ez, ih, jh
  real(real64) :: a(sx-ih:ex+ih, sy-jh:ey+jh, sz:ez)
  type(MPI_STATUS) :: status
  integer :: xl, yl, zl
  type(MPI_REQUEST) :: reqn, reqs, reqe, reqw
  type(MPI_REQUEST) :: reqrn, reqrs, reqre, reqrw
  integer nssize, ewsize
  integer :: i, j, k, ii
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
    ! TODO: allocate these once
    allocate(sendn(nssize),sends(nssize),recvn(nssize),recvs(nssize))
    !$acc enter data copyin(sendn, sends, recvn, recvs)

    !$acc parallel loop collapse(3) default(present) private(ii)
    do k = 1, zl
      do j = 1, jh
        do i = 1, xl
          ii = i + (j-1)*xl + (k-1)*xl*jh
          sendn(ii) = a(sx-ih+i-1,ey-jh+j,k)
          sends(ii) = a(sx-ih+i-1,sy+j-1,k)
        end do
      end do
    end do

    !   Send north/south
    call D_MPI_ISEND(sendn, nssize, nbrnorth, 4, comm3d, reqn, mpierr, ondevice=.true.)
    call D_MPI_ISEND(sends, nssize, nbrsouth, 5, comm3d, reqs, mpierr, ondevice=.true.)

    !   Receive south/north
    call D_MPI_IRECV(recvs, nssize, nbrsouth, 4, comm3d, reqrs, mpierr, ondevice=.true.)
    call D_MPI_IRECV(recvn, nssize, nbrnorth, 5, comm3d, reqrn, mpierr, ondevice=.true.)

    ! Wait until data is received
    call MPI_WAIT(reqrs, status, mpierr)
    call MPI_WAIT(reqrn, status, mpierr)


    ! Write back buffers
    !$acc parallel loop collapse(3) default(present) private(ii)
    do k = 1, zl
      do j = 1, jh
        do i = 1, xl
          ii = i + (j-1)*xl + (k-1)*xl*jh
          a(sx-ih+i-1,ey+j,k) = recvn(ii)
          a(sx-ih+i-1,sy-jh+(j-1),k) = recvs(ii)
        end do
      end do
    end do
    
  else

    ! Single processor, make sure the field is periodic
    !$acc kernels default(present) async
    a(:,sy-jh:sy-1,:) = a(:,ey-jh+1:ey,:)
    a(:,ey+1:ey+jh,:) = a(:,sy:sy+jh-1,:)
    !$acc end kernels

  endif

  if(nprocx .gt. 1)then

    !   Allocate send / receive buffers
    allocate(sende(ewsize),sendw(ewsize),recve(ewsize),recvw(ewsize))
    !$acc enter data copyin(sende, sendw, recve, recvw)

    !$acc parallel loop collapse(3) default(present) private(ii)
    do k = 1, zl
      do j = 1, yl
        do i = 1, ih
          ii = i + (j-1)*ih + (k-1)*ih*yl
          sende(ii) = a(ex-ih+i,sy-jh+j-1,k)
          sendw(ii) = a(sx+i-1,sy-jh+j-1,k)
        end do
      end do
    end do

    !   Send east/west
    call D_MPI_ISEND(sende, ewsize, nbreast, 6, comm3d, reqe, mpierr, ondevice=.true.)
    call D_MPI_ISEND(sendw, ewsize, nbrwest, 7, comm3d, reqw, mpierr, ondevice=.true.)

    !   Receive west/east
    call D_MPI_IRECV(recvw, ewsize, nbrwest, 6, comm3d, reqrw, mpierr, ondevice=.true.)
    call D_MPI_IRECV(recve, ewsize, nbreast, 7, comm3d, reqre, mpierr, ondevice=.true.)

    ! Wait until data is received
    call MPI_WAIT(reqrw, status, mpierr)
    call MPI_WAIT(reqre, status, mpierr)

    ! Write back buffers
    !$acc parallel loop collapse(3) default(present) private(ii)
    do k = 1, zl
      do j = 1, yl
        do i = 1, ih
          ii = i + (j-1)*ih + (k-1)*ih*yl
          a(sx-ih+(i-1),sy-jh+j-1,k) = recvw(ii)
          a(ex+i,sy-jh+j-1,k) = recve(ii)
        end do
      end do
    end do

  else

    ! Single processor, make sure the field is periodic
    !$acc kernels default(present) async
    a(sx-ih:sx-1,:,:) = a(ex-ih+1:ex,:,:)
    a(ex+1:ex+ih,:,:) = a(sx:sx+ih-1,:,:)
    !$acc end kernels
  endif

  if(nprocy.gt.1)then

    ! Make sure data is sent
    call MPI_WAIT(reqn, status, mpierr)
    if (mpierr /= MPI_SUCCESS) call abort
    call MPI_WAIT(reqs, status, mpierr)
    if (mpierr /= MPI_SUCCESS) call abort

    !$acc exit data delete(sendn, sends, recvn, recvs)
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
    !$acc exit data delete(sende, sendw, recve, recvw)
    deallocate (sende, sendw)
    deallocate (recve, recvw)

  endif

  !$acc wait

  end subroutine excjs_real64

  subroutine excjs_logical(a,sx,ex,sy,ey,sz,ez,ih,jh)
  implicit none
  integer sx, ex, sy, ey, sz, ez, ih, jh
  logical a(sx-ih:ex+ih, sy-jh:ey+jh, sz:ez)
  type(MPI_STATUS)  :: status
  integer :: xl, yl, zl
  type(MPI_REQUEST) :: reqn, reqs, reqe, reqw
  type(MPI_REQUEST) :: reqrn, reqrs, reqre, reqrw
  integer nssize, ewsize
  logical,allocatable, dimension(:) :: sendn,recvn &
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
  end subroutine excjs_logical

  subroutine slabsum_real32(aver,ks,kf,var,ib,ie,jb,je,kb,ke,ibs,ies,jbs,jes,kbs,kes,on_gpu)
    implicit none

    integer           :: ks,kf
    integer           :: ib,ie,jb,je,kb,ke,ibs,ies,jbs,jes,kbs,kes
    real(real32)      :: aver(ks:kf)
    real(real32)      :: var (ib:ie,jb:je,kb:ke)
    logical, optional :: on_gpu
    real(real32)      :: averl(ks:kf)
    real(real32)      :: avers(ks:kf)
    integer           :: k


    if (present(on_gpu)) then
      !$acc kernels default(present)
      do k = kbs, kes
        aver(k) = aver(k) + sum(var(ibs:ies, jbs:jes, k))
      end do
      !$acc end kernels
      call MPI_ALLREDUCE(MPI_IN_PLACE, aver, kf-ks+1, MPI_REAL4, MPI_SUM, comm3d, mpierr)
    else
      averl       = 0.
      avers       = 0.
      do k=kbs,kes
        averl(k) = sum(var(ibs:ies,jbs:jes,k))
      enddo
      call MPI_ALLREDUCE(averl, avers, kf-ks+1,  MPI_REAL4, &
                         MPI_SUM, comm3d,mpierr)
      
      aver = aver + avers
    endif

    return
  end subroutine slabsum_real32

  subroutine slabsum_real64(aver,ks,kf,var,ib,ie,jb,je,kb,ke,ibs,ies,jbs,jes,kbs,kes,on_gpu)
    implicit none

    integer           :: ks,kf
    integer           :: ib,ie,jb,je,kb,ke,ibs,ies,jbs,jes,kbs,kes
    real(real64)      :: aver(ks:kf)
    real(real64)      :: var (ib:ie,jb:je,kb:ke)
    logical, optional :: on_gpu
    real(real64)      :: averl(ks:kf)
    real(real64)      :: avers(ks:kf)
    integer           :: k

    
    if (present(on_gpu)) then
      !$acc kernels default(present)
      do k = kbs, kes
        aver(k) = aver(k) + sum(var(ibs:ies, jbs:jes, k))
      end do
      !$acc end kernels
      call MPI_ALLREDUCE(MPI_IN_PLACE, aver, kf-ks+1, MPI_REAL8, MPI_SUM, comm3d, mpierr)
    else
      averl       = 0.
      avers       = 0.
      do k=kbs,kes
        averl(k) = sum(var(ibs:ies,jbs:jes,k))
      enddo
      call MPI_ALLREDUCE(averl, avers, kf-ks+1,  MPI_REAL8, &
                         MPI_SUM, comm3d,mpierr)

      aver = aver + avers
    endif

    return
  end subroutine slabsum_real64

  subroutine slabsum_real32_5fields(aver1,ks,kf,var1,ib,ie,jb,je,kb,ke,ibs,ies,jbs,jes,kbs,kes, &
                                    aver2,      var2, &
                                    aver3,      var3, &
                                    aver4,      var4, &
                                    aver5,      var5, &
                                    on_gpu)
    implicit none

    integer :: ks, kf
    integer :: ib, ie, jb, je, kb, ke, ibs, ies, jbs, jes, kbs, kes
    real(real32) :: aver1(ks:kf), aver2(ks:kf), aver3(ks:kf), &
                    aver4(ks:kf), aver5(ks:kf)
    real(real32) :: var1(ib:ie, jb:je, kb:ke)
    real(real32) :: var2(ib:ie, jb:je, kb:ke)
    real(real32) :: var3(ib:ie, jb:je, kb:ke)
    real(real32) :: var4(ib:ie, jb:je, kb:ke)
    real(real32) :: var5(ib:ie, jb:je, kb:ke)
    logical, optional :: on_gpu
    real(real32), dimension(:,:), allocatable :: sum2d
    real(real32) :: sum_lcl1, sum_lcl2, sum_lcl3, sum_lcl4, sum_lcl5
    integer :: i,j,k

    allocate(sum2d(kf-ks+1,5))
    if (present(on_gpu)) then
      !$acc enter data create(sum2d)
      !$acc parallel loop gang default(present) private(sum_lcl1, sum_lcl2, sum_lcl3, sum_lcl4, sum_lcl5)
      do k = kbs, kes
        sum_lcl1 = 0.0
        sum_lcl2 = 0.0
        sum_lcl3 = 0.0
        sum_lcl4 = 0.0
        sum_lcl5 = 0.0
        !$acc loop collapse(2) reduction(+:sum_lcl1, sum_lcl2, sum_lcl3, sum_lcl4, sum_lcl5)
        do j = jbs, jes
          do i = ibs, ies
            sum_lcl1 = sum_lcl1 + var1(i, j, k)
            sum_lcl2 = sum_lcl2 + var2(i, j, k)
            sum_lcl3 = sum_lcl3 + var3(i, j, k)
            sum_lcl4 = sum_lcl4 + var4(i, j, k)
            sum_lcl5 = sum_lcl5 + var5(i, j, k)
          end do
        end do
        sum2d(k,1) = sum_lcl1
        sum2d(k,2) = sum_lcl2
        sum2d(k,3) = sum_lcl3
        sum2d(k,4) = sum_lcl4
        sum2d(k,5) = sum_lcl5
      end do
      !$acc exit data copyout(sum2d)
    else
      sum2d = 0
      do k = kbs, kes
        do j = jbs, jes
          do i = ibs, ies
            sum2d(k,1) = sum2d(k,1) + var1(i, j, k)
            sum2d(k,2) = sum2d(k,2) + var2(i, j, k)
            sum2d(k,3) = sum2d(k,3) + var3(i, j, k)
            sum2d(k,4) = sum2d(k,4) + var4(i, j, k)
            sum2d(k,5) = sum2d(k,5) + var5(i, j, k)
          end do
        end do
      end do
    endif

    call MPI_ALLREDUCE(MPI_IN_PLACE, sum2d, (kf-ks+1)*5, MPI_REAL4, MPI_SUM, comm3d, mpierr)
    aver1(:) = sum2d(:,1)
    aver2(:) = sum2d(:,2)
    aver3(:) = sum2d(:,3)
    aver4(:) = sum2d(:,4)
    aver5(:) = sum2d(:,5)

    if (present(on_gpu)) then
      !$acc exit data delete(sum2d)
      deallocate(sum2d)
    else
      deallocate(sum2d)
    endif

  end subroutine slabsum_real32_5fields

  subroutine slabsum_real32_4fields(aver1,ks,kf,var1,ib,ie,jb,je,kb,ke,ibs,ies,jbs,jes,kbs,kes, &
                                    aver2,      var2, &
                                    aver3,      var3, &
                                    aver4,      var4, &
                                    on_gpu)
    implicit none

    integer :: ks, kf
    integer :: ib, ie, jb, je, kb, ke, ibs, ies, jbs, jes, kbs, kes
    real(real32) :: aver1(ks:kf), aver2(ks:kf), aver3(ks:kf), &
                    aver4(ks:kf)
    real(real32) :: var1(ib:ie, jb:je, kb:ke)
    real(real32) :: var2(ib:ie, jb:je, kb:ke)
    real(real32) :: var3(ib:ie, jb:je, kb:ke)
    real(real32) :: var4(ib:ie, jb:je, kb:ke)
    logical, optional :: on_gpu
    real(real32), dimension(:,:), allocatable :: sum2d
    real(real32) :: sum_lcl1, sum_lcl2, sum_lcl3, sum_lcl4
    integer :: i,j,k

    allocate(sum2d(kf-ks+1,4))
    if (present(on_gpu)) then
      !$acc enter data create(sum2d)
      !$acc parallel loop gang default(present) private(sum_lcl1, sum_lcl2, sum_lcl3, sum_lcl4)
      do k = kbs, kes
        sum_lcl1 = 0.0
        sum_lcl2 = 0.0
        sum_lcl3 = 0.0
        sum_lcl4 = 0.0
        !$acc loop collapse(2) reduction(+:sum_lcl1, sum_lcl2, sum_lcl3, sum_lcl4)
        do j = jbs, jes
          do i = ibs, ies
            sum_lcl1 = sum_lcl1 + var1(i, j, k)
            sum_lcl2 = sum_lcl2 + var2(i, j, k)
            sum_lcl3 = sum_lcl3 + var3(i, j, k)
            sum_lcl4 = sum_lcl4 + var4(i, j, k)
          end do
        end do
        sum2d(k,1) = sum_lcl1
        sum2d(k,2) = sum_lcl2
        sum2d(k,3) = sum_lcl3
        sum2d(k,4) = sum_lcl4
      end do
      !$acc exit data copyout(sum2d)
    else
      sum2d = 0
      do k = kbs, kes
        do j = jbs, jes
          do i = ibs, ies
            sum2d(k,1) = sum2d(k,1) + var1(i, j, k)
            sum2d(k,2) = sum2d(k,2) + var2(i, j, k)
            sum2d(k,3) = sum2d(k,3) + var3(i, j, k)
            sum2d(k,4) = sum2d(k,4) + var4(i, j, k)
          end do
        end do
      end do
    endif

    call MPI_ALLREDUCE(MPI_IN_PLACE, sum2d, (kf-ks+1)*4, MPI_REAL4, MPI_SUM, comm3d, mpierr)
    aver1(:) = sum2d(:,1)
    aver2(:) = sum2d(:,2)
    aver3(:) = sum2d(:,3)
    aver4(:) = sum2d(:,4)

    if (present(on_gpu)) then
      !$acc exit data delete(sum2d)
      deallocate(sum2d)
    else
      deallocate(sum2d)
    endif

  end subroutine slabsum_real32_4fields

  subroutine slabsum_real32_3fields(aver1,ks,kf,var1,ib,ie,jb,je,kb,ke,ibs,ies,jbs,jes,kbs,kes, &
                                    aver2,      var2, &
                                    aver3,      var3, &
                                    on_gpu)
    implicit none

    integer :: ks, kf
    integer :: ib, ie, jb, je, kb, ke, ibs, ies, jbs, jes, kbs, kes
    real(real32) :: aver1(ks:kf), aver2(ks:kf), aver3(ks:kf)
    real(real32) :: var1(ib:ie, jb:je, kb:ke)
    real(real32) :: var2(ib:ie, jb:je, kb:ke)
    real(real32) :: var3(ib:ie, jb:je, kb:ke)
    logical, optional :: on_gpu
    real(real32), dimension(:,:), allocatable :: sum2d
    real(real32) :: sum_lcl1, sum_lcl2, sum_lcl3
    integer :: i,j,k

    allocate(sum2d(kf-ks+1,3))

    if (present(on_gpu)) then
      !$acc enter data create(sum2d)
      !$acc parallel loop gang default(present) private(sum_lcl1, sum_lcl2, sum_lcl3)
      do k = kbs, kes
        sum_lcl1 = 0.0
        sum_lcl2 = 0.0
        sum_lcl3 = 0.0
        !$acc loop collapse(2) reduction(+:sum_lcl1, sum_lcl2, sum_lcl3)
        do j = jbs, jes
          do i = ibs, ies
            sum_lcl1 = sum_lcl1 + var1(i, j, k)
            sum_lcl2 = sum_lcl2 + var2(i, j, k)
            sum_lcl3 = sum_lcl3 + var3(i, j, k)
          end do
        end do
        sum2d(k,1) = sum_lcl1
        sum2d(k,2) = sum_lcl2
        sum2d(k,3) = sum_lcl3
      end do
      !$acc exit data copyout(sum2d)
    else
      sum2d = 0
      do k = kbs, kes
        do j = jbs, jes
          do i = ibs, ies
            sum2d(k,1) = sum2d(k,1) + var1(i, j, k)
            sum2d(k,2) = sum2d(k,2) + var2(i, j, k)
            sum2d(k,3) = sum2d(k,3) + var3(i, j, k)
          end do
        end do
      end do
    endif

    call MPI_ALLREDUCE(MPI_IN_PLACE, sum2d, (kf-ks+1)*3, MPI_REAL4, MPI_SUM, comm3d, mpierr)
    aver1(:) = sum2d(:,1)
    aver2(:) = sum2d(:,2)
    aver3(:) = sum2d(:,3)

    if (present(on_gpu)) then
      !$acc exit data delete(sum2d)
      deallocate(sum2d)
    else
      deallocate(sum2d)
    endif

  end subroutine slabsum_real32_3fields

  subroutine slabsum_real32_2fields(aver1,ks,kf,var1,ib,ie,jb,je,kb,ke,ibs,ies,jbs,jes,kbs,kes, &
                                    aver2,      var2, &
                                    on_gpu)
    implicit none

    integer :: ks, kf
    integer :: ib, ie, jb, je, kb, ke, ibs, ies, jbs, jes, kbs, kes
    real(real32) :: aver1(ks:kf), aver2(ks:kf)
    real(real32) :: var1(ib:ie, jb:je, kb:ke)
    real(real32) :: var2(ib:ie, jb:je, kb:ke)
    logical, optional :: on_gpu
    real(real32), dimension(:,:), allocatable :: sum2d
    real(real32) :: sum_lcl1, sum_lcl2
    integer :: i,j,k

    allocate(sum2d(kf-ks+1,2))

    if (present(on_gpu)) then
      !$acc enter data create(sum2d)
      !$acc parallel loop gang default(present) private(sum_lcl1, sum_lcl2)
      do k = kbs, kes
        sum_lcl1 = 0.0
        sum_lcl2 = 0.0
        !$acc loop collapse(2) reduction(+:sum_lcl1, sum_lcl2)
        do j = jbs, jes
          do i = ibs, ies
            sum_lcl1 = sum_lcl1 + var1(i, j, k)
            sum_lcl2 = sum_lcl2 + var2(i, j, k)
          end do
        end do
        sum2d(k,1) = sum_lcl1
        sum2d(k,2) = sum_lcl2
      end do
      !$acc exit data copyout(sum2d)
    else
      sum2d = 0
      do k = kbs, kes
        do j = jbs, jes
          do i = ibs, ies
            sum2d(k,1) = sum2d(k,1) + var1(i, j, k)
            sum2d(k,2) = sum2d(k,2) + var2(i, j, k)
          end do
        end do
      end do
    endif

    call MPI_ALLREDUCE(MPI_IN_PLACE, sum2d, (kf-ks+1)*2, MPI_REAL4, MPI_SUM, comm3d, mpierr)
    aver1(:) = sum2d(:,1)
    aver2(:) = sum2d(:,2)

    if (present(on_gpu)) then
      !$acc exit data delete(sum2d)
      deallocate(sum2d)
    else
      deallocate(sum2d)
    endif

  end subroutine slabsum_real32_2fields

  subroutine slabsum_real64_5fields(aver1,ks,kf,var1,ib,ie,jb,je,kb,ke,ibs,ies,jbs,jes,kbs,kes, &
                                    aver2,      var2, &
                                    aver3,      var3, &
                                    aver4,      var4, &
                                    aver5,      var5, &
                                    on_gpu)
    implicit none

    integer :: ks, kf
    integer :: ib, ie, jb, je, kb, ke, ibs, ies, jbs, jes, kbs, kes
    real(real64) :: aver1(ks:kf), aver2(ks:kf), aver3(ks:kf), &
                    aver4(ks:kf), aver5(ks:kf)
    real(real64) :: var1(ib:ie, jb:je, kb:ke)
    real(real64) :: var2(ib:ie, jb:je, kb:ke)
    real(real64) :: var3(ib:ie, jb:je, kb:ke)
    real(real64) :: var4(ib:ie, jb:je, kb:ke)
    real(real64) :: var5(ib:ie, jb:je, kb:ke)
    logical, optional :: on_gpu
    real(real64), dimension(:,:), allocatable :: sum2d
    real(real64) :: sum_lcl1, sum_lcl2, sum_lcl3, sum_lcl4, sum_lcl5
    integer :: i,j,k

    allocate(sum2d(kf-ks+1,5))

    if (present(on_gpu)) then
      !$acc enter data create(sum2d)
      !$acc parallel loop gang default(present) private(sum_lcl1, sum_lcl2, sum_lcl3, sum_lcl4, sum_lcl5)
      do k = kbs, kes
        sum_lcl1 = 0.0
        sum_lcl2 = 0.0
        sum_lcl3 = 0.0
        sum_lcl4 = 0.0
        sum_lcl5 = 0.0
        !$acc loop collapse(2) reduction(+:sum_lcl1, sum_lcl2, sum_lcl3, sum_lcl4, sum_lcl5)
        do j = jbs, jes
          do i = ibs, ies
            sum_lcl1 = sum_lcl1 + var1(i, j, k)
            sum_lcl2 = sum_lcl2 + var2(i, j, k)
            sum_lcl3 = sum_lcl3 + var3(i, j, k)
            sum_lcl4 = sum_lcl4 + var4(i, j, k)
            sum_lcl5 = sum_lcl5 + var5(i, j, k)
          end do
        end do
        sum2d(k,1) = sum_lcl1
        sum2d(k,2) = sum_lcl2
        sum2d(k,3) = sum_lcl3
        sum2d(k,4) = sum_lcl4
        sum2d(k,5) = sum_lcl5
      end do
      !$acc exit data copyout(sum2d)
    else
      sum2d = 0
      do k = kbs, kes
        do j = jbs, jes
          do i = ibs, ies
            sum2d(k,1) = sum2d(k,1) + var1(i, j, k)
            sum2d(k,2) = sum2d(k,2) + var2(i, j, k)
            sum2d(k,3) = sum2d(k,3) + var3(i, j, k)
            sum2d(k,4) = sum2d(k,4) + var4(i, j, k)
            sum2d(k,5) = sum2d(k,5) + var5(i, j, k)
          end do
        end do
      end do
    endif

    call MPI_ALLREDUCE(MPI_IN_PLACE, sum2d, (kf-ks+1)*5, MPI_REAL8, MPI_SUM, comm3d, mpierr)
    aver1(:) = sum2d(:,1)
    aver2(:) = sum2d(:,2)
    aver3(:) = sum2d(:,3)
    aver4(:) = sum2d(:,4)
    aver5(:) = sum2d(:,5)

    if (present(on_gpu)) then
      !$acc exit data delete(sum2d)
      deallocate(sum2d)
    else
      deallocate(sum2d)
    endif

  end subroutine slabsum_real64_5fields

  subroutine slabsum_real64_4fields(aver1,ks,kf,var1,ib,ie,jb,je,kb,ke,ibs,ies,jbs,jes,kbs,kes, &
                                    aver2,      var2, &
                                    aver3,      var3, &
                                    aver4,      var4, &
                                    on_gpu)
    implicit none

    integer :: ks, kf
    integer :: ib, ie, jb, je, kb, ke, ibs, ies, jbs, jes, kbs, kes
    real(real64) :: aver1(ks:kf), aver2(ks:kf), aver3(ks:kf), &
                    aver4(ks:kf)
    real(real64) :: var1(ib:ie, jb:je, kb:ke)
    real(real64) :: var2(ib:ie, jb:je, kb:ke)
    real(real64) :: var3(ib:ie, jb:je, kb:ke)
    real(real64) :: var4(ib:ie, jb:je, kb:ke)
    logical, optional :: on_gpu
    real(real64), dimension(:,:), allocatable :: sum2d
    real(real64) :: sum_lcl1, sum_lcl2, sum_lcl3, sum_lcl4
    integer :: i,j,k

    allocate(sum2d(kf-ks+1,4))

    if (present(on_gpu)) then
      !$acc enter data create(sum2d)
      !$acc parallel loop gang default(present) private(sum_lcl1, sum_lcl2, sum_lcl3, sum_lcl4)
      do k = kbs, kes
        sum_lcl1 = 0.0
        sum_lcl2 = 0.0
        sum_lcl3 = 0.0
        sum_lcl4 = 0.0
        !$acc loop collapse(2) reduction(+:sum_lcl1, sum_lcl2, sum_lcl3, sum_lcl4)
        do j = jbs, jes
          do i = ibs, ies
            sum_lcl1 = sum_lcl1 + var1(i, j, k)
            sum_lcl2 = sum_lcl2 + var2(i, j, k)
            sum_lcl3 = sum_lcl3 + var3(i, j, k)
            sum_lcl4 = sum_lcl4 + var4(i, j, k)
          end do
        end do
        sum2d(k,1) = sum_lcl1
        sum2d(k,2) = sum_lcl2
        sum2d(k,3) = sum_lcl3
        sum2d(k,4) = sum_lcl4
      end do
      !$acc exit data copyout(sum2d)
    else
      sum2d = 0
      do k = kbs, kes
        do j = jbs, jes
          do i = ibs, ies
            sum2d(k,1) = sum2d(k,1) + var1(i, j, k)
            sum2d(k,2) = sum2d(k,2) + var2(i, j, k)
            sum2d(k,3) = sum2d(k,3) + var3(i, j, k)
            sum2d(k,4) = sum2d(k,4) + var4(i, j, k)
          end do
        end do
      end do
    endif

    call MPI_ALLREDUCE(MPI_IN_PLACE, sum2d, (kf-ks+1)*4, MPI_REAL8, MPI_SUM, comm3d, mpierr)
    aver1(:) = sum2d(:,1)
    aver2(:) = sum2d(:,2)
    aver3(:) = sum2d(:,3)
    aver4(:) = sum2d(:,4)

    if (present(on_gpu)) then
      !$acc exit data delete(sum2d)
      deallocate(sum2d)
    else
      deallocate(sum2d)
    endif

  end subroutine slabsum_real64_4fields

  subroutine slabsum_real64_3fields(aver1,ks,kf,var1,ib,ie,jb,je,kb,ke,ibs,ies,jbs,jes,kbs,kes, &
                                    aver2,      var2, &
                                    aver3,      var3, &
                                    on_gpu)
    implicit none

    integer :: ks, kf
    integer :: ib, ie, jb, je, kb, ke, ibs, ies, jbs, jes, kbs, kes
    real(real64) :: aver1(ks:kf), aver2(ks:kf), aver3(ks:kf)
    real(real64) :: var1(ib:ie, jb:je, kb:ke)
    real(real64) :: var2(ib:ie, jb:je, kb:ke)
    real(real64) :: var3(ib:ie, jb:je, kb:ke)
    logical, optional :: on_gpu
    real(real64), dimension(:,:), allocatable :: sum2d
    real(real64) :: sum_lcl1, sum_lcl2, sum_lcl3
    integer :: i,j,k

    allocate(sum2d(kf-ks+1,3))

    if (present(on_gpu)) then
      !$acc enter data create(sum2d)
      !$acc parallel loop gang default(present) private(sum_lcl1, sum_lcl2, sum_lcl3)
      do k = kbs, kes
        sum_lcl1 = 0.0
        sum_lcl2 = 0.0
        sum_lcl3 = 0.0
        !$acc loop collapse(2) reduction(+:sum_lcl1, sum_lcl2, sum_lcl3)
        do j = jbs, jes
          do i = ibs, ies
            sum_lcl1 = sum_lcl1 + var1(i, j, k)
            sum_lcl2 = sum_lcl2 + var2(i, j, k)
            sum_lcl3 = sum_lcl3 + var3(i, j, k)
          end do
        end do
        sum2d(k,1) = sum_lcl1
        sum2d(k,2) = sum_lcl2
        sum2d(k,3) = sum_lcl3
      end do
      !$acc exit data copyout(sum2d)
    else
      sum2d = 0
      do k = kbs, kes
        do j = jbs, jes
          do i = ibs, ies
            sum2d(k,1) = sum2d(k,1) + var1(i, j, k)
            sum2d(k,2) = sum2d(k,2) + var2(i, j, k)
            sum2d(k,3) = sum2d(k,3) + var3(i, j, k)
          end do
        end do
      end do
    endif

    call MPI_ALLREDUCE(MPI_IN_PLACE, sum2d, (kf-ks+1)*3, MPI_REAL8, MPI_SUM, comm3d, mpierr)
    aver1(:) = sum2d(:,1)
    aver2(:) = sum2d(:,2)
    aver3(:) = sum2d(:,3)

    if (present(on_gpu)) then
      !$acc exit data delete(sum2d)
      deallocate(sum2d)
    else
      deallocate(sum2d)
    endif

  end subroutine slabsum_real64_3fields

  subroutine slabsum_real64_2fields(aver1,ks,kf,var1,ib,ie,jb,je,kb,ke,ibs,ies,jbs,jes,kbs,kes, &
                                    aver2,      var2, &
                                    on_gpu)
    implicit none

    integer :: ks, kf
    integer :: ib, ie, jb, je, kb, ke, ibs, ies, jbs, jes, kbs, kes
    real(real64) :: aver1(ks:kf), aver2(ks:kf)
    real(real64) :: var1(ib:ie, jb:je, kb:ke)
    real(real64) :: var2(ib:ie, jb:je, kb:ke)
    logical, optional :: on_gpu
    real(real64), dimension(:,:), allocatable :: sum2d
    real(real64) :: sum_lcl1, sum_lcl2
    integer :: i,j,k

    allocate(sum2d(kf-ks+1,2))

    if (present(on_gpu)) then
      !$acc enter data create(sum2d)
      !$acc parallel loop gang default(present) private(sum_lcl1, sum_lcl2)
      do k = kbs, kes
        sum_lcl1 = 0.0
        sum_lcl2 = 0.0
        !$acc loop collapse(2) reduction(+:sum_lcl1, sum_lcl2)
        do j = jbs, jes
          do i = ibs, ies
            sum_lcl1 = sum_lcl1 + var1(i, j, k)
            sum_lcl2 = sum_lcl2 + var2(i, j, k)
          end do
        end do
        sum2d(k,1) = sum_lcl1
        sum2d(k,2) = sum_lcl2
      end do
      !$acc exit data copyout(sum2d)
    else
      sum2d = 0
      do k = kbs, kes
        do j = jbs, jes
          do i = ibs, ies
            sum2d(k,1) = sum2d(k,1) + var1(i, j, k)
            sum2d(k,2) = sum2d(k,2) + var2(i, j, k)
          end do
        end do
      end do
    endif

    call MPI_ALLREDUCE(MPI_IN_PLACE, sum2d, (kf-ks+1)*2, MPI_REAL8, MPI_SUM, comm3d, mpierr)
    aver1(:) = sum2d(:,1)
    aver2(:) = sum2d(:,2)

    if (present(on_gpu)) then
      !$acc exit data delete(sum2d)
      deallocate(sum2d)
    else
      deallocate(sum2d)
    endif

  end subroutine slabsum_real64_2fields

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

  subroutine closeboundaries(field, sx,ex,hx, sy,ey,hy, sz,ez, fillvalue_in)

    ! Clears halo cells of domain edges
    ! Author: Marco de Bruine (VU), 2020
    use modprecision, only: field_r
    implicit none
    integer sx,ex,hx, sy,ey,hy, sz,ez
    real(field_r) :: field(sx-hx:ex+hx,sy-hy:ey+hy,sz:ez)
    real :: fillvalue = 0
    real(field_r),optional :: fillvalue_in
    if(present(fillvalue_in)) fillvalue=fillvalue_in

    if      (myidx == 0) then        ! WEST boundary
      field(sx-hx:sx-1,  sy-hy:ey+hy, sz:ez) = fillvalue
    else if (myidx == nprocx-1) then ! EAST boundary
      field(ex+1 :ex+hx, sy-hy:ey+hy, sz:ez) = fillvalue
    end if

    if      (myidy == 0) then        ! SOUTH boundary
      field(sx-hx:ex+hx, sy-hy:sy-1,  sz:ez) = fillvalue
    else if (myidy == nprocy-1) then ! NORTH boundary
      field(sx-hx:ex+hx, ey+1:ey+hy, sz:ez) = fillvalue
    end if

  end subroutine closeboundaries
end module
