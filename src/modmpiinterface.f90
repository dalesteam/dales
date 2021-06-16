!> \file modmpiinterface.f90
!
! \author: Victor Azizi
!
! This file contains the separate interfaces for communicating with MPI
! The seperate interfaces are necessary to send the correct number of bytes
!
! Not all compilers accept (..) as an argument, therefore we write all possibly
! used combinations explicitly (_S _R1 _R2 _R3) where needed

! Some compilers (fujitsu) do not accept REAL32 arrays, so we cast them to bytes

! _IP stands for in  place
! _S stands for scalar
! _R1 stands for an array with rank 1
! _R2 stands for an array with rank 2
! _R3 stands for an array with rank 3
module modmpiinterface
use mpi_f08
use iso_c_binding  , only : c_loc, c_f_pointer, c_ptr
use iso_fortran_env, only : real32,real64,int32
implicit none

contains

!>D_MPI_ISEND
  !MPI interfaces instantations for the various types
  subroutine D_MPI_ISEND_REAL32_R1(buf, count, dest, tag, comm, request, ierror)
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
    type(c_ptr) :: c_p
    c_p = c_loc(buf)
    call c_f_pointer(c_p,buf_b)
    call MPI_ISEND(buf_b,count,MPI_REAL4,dest,tag,comm,request,ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ISEND_REAL32_R1
  subroutine D_MPI_ISEND_REAL64_R1(buf, count, dest, tag, comm, request, ierror)
    implicit none
    real(real64), contiguous, asynchronous, intent(inout) ::   buf(:)
    integer       ::   count, dest, tag, ierror
    type(MPI_COMM):: comm
    type(MPI_REQUEST) :: request
    call MPI_ISEND(buf,count,MPI_REAL8,dest,tag,comm,request,ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ISEND_REAL64_R1

!>D_MPI_IRECV
  subroutine D_MPI_IRECV_REAL32_R1(buf, count, source, tag, comm, request, ierror)
    implicit none
    real(real32), target, asynchronous,contiguous, intent(inout)  ::   buf(:)
    integer        :: count, source, tag, ierror
    type(MPI_COMM) :: comm
    type(MPI_REQUEST) :: request
    BYTE, pointer :: buf_b
    type(c_ptr) :: c_p
    c_p = c_loc(buf)
    call c_f_pointer(c_p,buf_b)
    call MPI_IRECV(buf_b,count,MPI_REAL4,source,tag,comm,request,ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_IRECV_REAL32_R1
  subroutine D_MPI_IRECV_REAL64_R1(buf, count, source, tag, comm, request, ierror)
    implicit none
    real(real64), contiguous, asynchronous, intent(inout)  ::   buf(:)
    integer        :: count, source, tag, ierror
    type(MPI_COMM) :: comm
    type(MPI_REQUEST) :: request
    call MPI_IRECV(buf,count,MPI_REAL8,source,tag,comm,request,ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_IRECV_REAL64_R1
  
!>D_MPI_RECV
  subroutine D_MPI_RECV_REAL32_R1(buf, count, source, tag, comm, status, ierror)
    implicit none
    real(real32), target, contiguous, intent(inout)  ::   buf(:)
    integer        :: count, source, tag, ierror
    type(MPI_COMM) :: comm
    type(MPI_STATUS) :: status
    BYTE, pointer :: buf_b
    type(c_ptr) :: c_p
    c_p = c_loc(buf)
    call c_f_pointer(c_p,buf_b)
    call MPI_RECV(buf_b,count,MPI_REAL4,source,tag,comm,status,ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_RECV_REAL32_R1
  subroutine D_MPI_RECV_REAL64_R1(buf, count, source, tag, comm, status, ierror)
    implicit none
    real(real64), contiguous, intent(inout)  ::   buf(:)
    integer        :: count, source, tag, ierror
    type(MPI_COMM) :: comm
    type(MPI_STATUS) :: status
    call MPI_RECV(buf,count,MPI_REAL8,source,tag,comm,status,ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_RECV_REAL64_R1

!>D_MPI_BCAST
  subroutine D_MPI_BCAST_REAL32_S(buffer, count, root, comm, ierror)
    implicit none
    real(real32), intent(inout)  ::  buffer
    integer        :: count, root, ierror
    type(MPI_COMM) :: comm
    call MPI_BCAST(buffer, count, MPI_REAL4, root, comm, ierror)
  end subroutine D_MPI_BCAST_REAL32_S
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
  subroutine D_MPI_BCAST_REAL32_R1(buffer, count, root, comm, ierror)
    implicit none
    real(real32), target, contiguous, intent(inout)   ::  buffer(:)
    integer        :: count, root, ierror
    type(MPI_COMM) :: comm
    BYTE, pointer :: buf_b
    type(c_ptr) :: c_p
    c_p = c_loc(buffer)
    call c_f_pointer(c_p,buf_b)
    call MPI_BCAST(buf_b, count, MPI_REAL4, root, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_BCAST_REAL32_R1
  subroutine D_MPI_BCAST_REAL32_R2(buffer, count, root, comm, ierror)
    implicit none
    real(real32), target, contiguous, intent(inout)   ::  buffer(:,:)
    integer        :: count, root, ierror
    type(MPI_COMM) :: comm
    BYTE, pointer :: buf_b
    type(c_ptr) :: c_p
    c_p = c_loc(buffer)
    call c_f_pointer(c_p,buf_b)
    call MPI_BCAST(buf_b, count, MPI_REAL4, root, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_BCAST_REAL32_R2
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

!>D_MPI_ALLREDUCE
  subroutine D_MPI_ALLREDUCE_REAL32_S(sendbuf, recvbuf, count, op, comm, ierror)
    implicit none
    real(real32), target, intent(inout)   :: sendbuf, recvbuf
    integer        :: count, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    BYTE, pointer :: bufs_b, bufr_b
    type(c_ptr) :: c_p
    c_p = c_loc(sendbuf)
    call c_f_pointer(c_p,bufs_b)
    c_p = c_loc(recvbuf)
    call c_f_pointer(c_p,bufr_b)
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
    type(c_ptr) :: c_p
    c_p = c_loc(sendbuf)
    call c_f_pointer(c_p,bufs_b)
    c_p = c_loc(recvbuf)
    call c_f_pointer(c_p,bufr_b)
    call MPI_ALLREDUCE(bufs_b, bufr_b, count, MPI_REAL4, op, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ALLREDUCE_REAL32_R1
  subroutine D_MPI_ALLREDUCE_REAL32_R2(sendbuf, recvbuf, count, op, comm, ierror)
    implicit none
    real(real32), target, contiguous, intent(inout)   :: sendbuf(:,:), recvbuf(:,:)
    integer        :: count, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    BYTE, pointer :: bufs_b, bufr_b
    type(c_ptr) :: c_p
    c_p = c_loc(sendbuf)
    call c_f_pointer(c_p,bufs_b)
    c_p = c_loc(recvbuf)
    call c_f_pointer(c_p,bufr_b)
    call MPI_ALLREDUCE(bufs_b, bufr_b, count, MPI_REAL4, op, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ALLREDUCE_REAL32_R2
  subroutine D_MPI_ALLREDUCE_REAL32_R3(sendbuf, recvbuf, count, op, comm, ierror)
    implicit none
    real(real32), target, contiguous, intent(inout)   :: sendbuf(:,:,:), recvbuf(:,:,:)
    integer        :: count, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    BYTE, pointer :: bufs_b, bufr_b
    type(c_ptr) :: c_p
    c_p = c_loc(sendbuf)
    call c_f_pointer(c_p,bufs_b)
    c_p = c_loc(recvbuf)
    call c_f_pointer(c_p,bufr_b)
    call MPI_ALLREDUCE(bufs_b, bufr_b, count, MPI_REAL4, op, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ALLREDUCE_REAL32_R3
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
    type(c_ptr) :: c_p
    c_p = c_loc(recvbuf)
    call c_f_pointer(c_p,bufr_b)
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

!>D_MPI_ALLTOALL
  subroutine D_MPI_ALLTOALL_REAL32_R1(sendbuf, sendcount, recvbuf, recvcount, comm, ierror)
    implicit none
    real(real32), target, contiguous, intent(inout)   :: sendbuf(:), recvbuf(:)
    integer        :: sendcount, recvcount, ierror
    type(MPI_COMM) :: comm
    BYTE, pointer :: bufs_b, bufr_b
    type(c_ptr) :: c_p
    c_p = c_loc(sendbuf)
    call c_f_pointer(c_p,bufs_b)
    c_p = c_loc(recvbuf)
    call c_f_pointer(c_p,bufr_b)
    call MPI_ALLTOALL(bufs_b, sendcount, MPI_REAL4, bufr_b, recvcount, MPI_REAL4, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ALLTOALL_REAL32_R1
  subroutine D_MPI_ALLTOALL_REAL64_R1(sendbuf, sendcount, recvbuf, recvcount, comm, ierror)
    implicit none
    real(real64), contiguous, intent(inout)   :: sendbuf(:), recvbuf(:)
    integer        :: sendcount, recvcount, ierror
    type(MPI_COMM) :: comm
    call MPI_ALLTOALL(sendbuf, sendcount, MPI_REAL8, recvbuf, recvcount, MPI_REAL8, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ALLTOALL_REAL64_R1

!>D_MPI_REDUCE
  subroutine D_MPI_REDUCE_REAL32_R1(sendbuf, recvbuf, count, op, root, comm, ierror)
    implicit none
    real(real32), target, contiguous, intent(inout)   :: sendbuf(:), recvbuf(:)
    integer        :: count, root, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    BYTE, pointer :: bufs_b, bufr_b
    type(c_ptr) :: c_p
    c_p = c_loc(sendbuf)
    call c_f_pointer(c_p,bufs_b)
    c_p = c_loc(recvbuf)
    call c_f_pointer(c_p,bufr_b)
    call MPI_REDUCE(bufs_b, bufr_b, count, MPI_REAL4, op, root, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_REDUCE_REAL32_R1
  subroutine D_MPI_REDUCE_REAL32_R2(sendbuf, recvbuf, count, op, root, comm, ierror)
    implicit none
    real(real32), target, contiguous, intent(inout)   :: sendbuf(:,:), recvbuf(:,:)
    integer        :: count, root, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    BYTE, pointer :: bufs_b, bufr_b
    type(c_ptr) :: c_p
    c_p = c_loc(sendbuf)
    call c_f_pointer(c_p,bufs_b)
    c_p = c_loc(recvbuf)
    call c_f_pointer(c_p,bufr_b)
    call MPI_REDUCE(bufs_b, bufr_b, count, MPI_REAL4, op, root, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_REDUCE_REAL32_R2
  subroutine D_MPI_REDUCE_REAL32_R3(sendbuf, recvbuf, count, op, root, comm, ierror)
    implicit none
    real(real32), target, contiguous, intent(inout)   :: sendbuf(:,:,:), recvbuf(:,:,:)
    integer        :: count, root, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    BYTE, pointer :: bufs_b, bufr_b
    type(c_ptr) :: c_p
    c_p = c_loc(sendbuf)
    call c_f_pointer(c_p,bufs_b)
    c_p = c_loc(recvbuf)
    call c_f_pointer(c_p,bufr_b)
    call MPI_REDUCE(bufs_b, bufr_b, count, MPI_REAL4, op, root, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_REDUCE_REAL32_R3
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
  subroutine D_MPI_REDUCE_REAL32_IP_R1(recvbuf, count, op, root, comm, ierror)
    implicit none
    real(real32), target, contiguous, intent(inout)   :: recvbuf(:)
    integer        :: count, root, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    BYTE, pointer :: bufr_b
    type(c_ptr) :: c_p
    c_p = c_loc(recvbuf)
    call c_f_pointer(c_p,bufr_b)
    call MPI_REDUCE(MPI_IN_PLACE, bufr_b, count, MPI_REAL4, op, root, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_REDUCE_REAL32_IP_R1
  subroutine D_MPI_REDUCE_REAL32_IP_R2(recvbuf, count, op, root, comm, ierror)
    implicit none
    real(real32), target, contiguous, intent(inout)   :: recvbuf(:,:)
    integer        :: count, root, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    BYTE, pointer :: bufr_b
    type(c_ptr) :: c_p
    c_p = c_loc(recvbuf)
    call c_f_pointer(c_p,bufr_b)
    call MPI_REDUCE(MPI_IN_PLACE, bufr_b, count, MPI_REAL4, op, root, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_REDUCE_REAL32_IP_R2
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

!>D_MPI_GATHER
  subroutine D_MPI_GATHER_REAL32_R1(sendbuf, sendcount, recvbuf, recvcount, root, comm, ierror)
    implicit none
    real(real32), target, contiguous, intent(inout)   :: sendbuf(:), recvbuf(:)
    integer        :: sendcount, recvcount, root, ierror
    type(MPI_COMM) :: comm
    BYTE, pointer :: bufs_b, bufr_b
    type(c_ptr) :: c_p
    c_p = c_loc(sendbuf)
    call c_f_pointer(c_p,bufs_b)
    c_p = c_loc(recvbuf)
    call c_f_pointer(c_p,bufr_b)
    call MPI_GATHER( bufs_b, sendcount, MPI_REAL4 &
                   , bufr_b, recvcount, MPI_REAL4 &
                   , root, comm, ierror )
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_GATHER_REAL32_R1
  subroutine D_MPI_GATHER_REAL64_R1(sendbuf, sendcount, recvbuf, recvcount, root, comm, ierror)
    implicit none
    real(real64), contiguous, intent(inout)   :: sendbuf(:), recvbuf(:)
    integer        :: sendcount, recvcount, root, ierror
    type(MPI_COMM) :: comm
    call MPI_GATHER( sendbuf, sendcount, MPI_REAL8 &
                   , recvbuf, recvcount, MPI_REAL8 &
                   , root, comm, ierror )
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_GATHER_REAL64_R1

end module
