!> \file modmpiinterface.f90
!
! \author: Victor Azizi
! \author: Caspar Jungbacker, TU Delft
!
! This file contains the separate interfaces for communicating with GPU-aware MPI
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
module modgpumpiinterface
use mpi_f08
use iso_c_binding  , only : c_loc, c_f_pointer, c_ptr
use iso_fortran_env, only : real32,real64,int32
implicit none

#if defined(_OPENACC)

contains

!>D_MPI_ISEND
  !MPI interfaces instantations for the various types
  subroutine D_MPI_ISEND_REAL32_R1_GPU(buf, count, dest, tag, comm, request, ierror)
    implicit none
    real(real32), device, contiguous, intent(inout) ::   buf(:)
    integer       ::   count, dest, tag, ierror
    type(MPI_COMM):: comm
    type(MPI_REQUEST) :: request
    call MPI_ISEND(buf,count,MPI_REAL4,dest,tag,comm,request,ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ISEND_REAL32_R1_GPU
  subroutine D_MPI_ISEND_REAL64_R1_GPU(buf, count, dest, tag, comm, request, ierror)
    implicit none
    real(real64), device, contiguous, intent(inout) ::   buf(:)
    integer       ::   count, dest, tag, ierror
    type(MPI_COMM):: comm
    type(MPI_REQUEST) :: request
    call MPI_ISEND(buf,count,MPI_REAL8,dest,tag,comm,request,ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ISEND_REAL64_R1_GPU
  subroutine D_MPI_ISEND_LOGICAL_R1_GPU(buf, count, dest, tag, comm, request, ierror)
    implicit none
    logical, device, contiguous, intent(inout) ::   buf(:)
    integer       ::   count, dest, tag, ierror
    type(MPI_COMM):: comm
    type(MPI_REQUEST) :: request
    call MPI_ISEND(buf,count,MPI_LOGICAL,dest,tag,comm,request,ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ISEND_LOGICAL_R1_GPU

!>D_MPI_IRECV
  subroutine D_MPI_IRECV_REAL32_R1_GPU(buf, count, source, tag, comm, request, ierror)
    implicit none
    real(real32), device, contiguous, intent(inout)  ::   buf(:)
    integer        :: count, source, tag, ierror
    type(MPI_COMM) :: comm
    type(MPI_REQUEST) :: request
    call MPI_IRECV(buf,count,MPI_REAL4,source,tag,comm,request,ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_IRECV_REAL32_R1_GPU
  subroutine D_MPI_IRECV_REAL64_R1_GPU(buf, count, source, tag, comm, request, ierror)
    implicit none
    real(real64), device, contiguous, intent(inout)  ::   buf(:)
    integer        :: count, source, tag, ierror
    type(MPI_COMM) :: comm
    type(MPI_REQUEST) :: request
    call MPI_IRECV(buf,count,MPI_REAL8,source,tag,comm,request,ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_IRECV_REAL64_R1_GPU
  subroutine D_MPI_IRECV_LOGICAL_R1_GPU(buf, count, source, tag, comm, request, ierror)
    implicit none
    logical, device, contiguous, intent(inout)  ::   buf(:)
    integer        :: count, source, tag, ierror
    type(MPI_COMM) :: comm
    type(MPI_REQUEST) :: request
    call MPI_IRECV(buf,count,MPI_LOGICAL,source,tag,comm,request,ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_IRECV_LOGICAL_R1_GPU
  
!>D_MPI_RECV
  subroutine D_MPI_RECV_REAL32_R1_GPU(buf, count, source, tag, comm, status, ierror)
    implicit none
    real(real32), device, contiguous, intent(inout)  ::   buf(:)
    integer        :: count, source, tag, ierror
    type(MPI_COMM) :: comm
    type(MPI_STATUS) :: status
    call MPI_RECV(buf,count,MPI_REAL4,source,tag,comm,status,ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_RECV_REAL32_R1_GPU
  subroutine D_MPI_RECV_REAL64_R1_GPU(buf, count, source, tag, comm, status, ierror)
    implicit none
    real(real64), device, contiguous, intent(inout)  ::   buf(:)
    integer        :: count, source, tag, ierror
    type(MPI_COMM) :: comm
    type(MPI_STATUS) :: status
    call MPI_RECV(buf,count,MPI_REAL8,source,tag,comm,status,ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_RECV_REAL64_R1_GPU

!>D_MPI_BCAST
  subroutine D_MPI_BCAST_REAL32_S_GPU(buffer, count, root, comm, ierror)
    implicit none
    real(real32), device, intent(inout)  ::  buffer
    integer        :: count, root, ierror
    type(MPI_COMM) :: comm
    call MPI_BCAST(buffer, count, MPI_REAL4, root, comm, ierror)
  end subroutine D_MPI_BCAST_REAL32_S_GPU
  subroutine D_MPI_BCAST_REAL64_S_GPU(buffer, count, root, comm, ierror)
    implicit none
    real(real64), device, intent(inout)  ::  buffer
    integer        :: count, root, ierror
    type(MPI_COMM) :: comm
    call MPI_BCAST(buffer, count, MPI_REAL8, root, comm, ierror)
  end subroutine D_MPI_BCAST_REAL64_S_GPU
  subroutine D_MPI_BCAST_INT32_S_GPU(buffer, count, root, comm, ierror)
    implicit none
    integer(int32), device, intent(inout) ::  buffer
    integer        :: count, root, ierror
    type(MPI_COMM) :: comm
    call MPI_BCAST(buffer, count, MPI_INTEGER4, root, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_BCAST_INT32_S_GPU
  subroutine D_MPI_BCAST_LOGICAL_S_GPU(buffer, count, root, comm, ierror)
    implicit none
    logical, device, intent(inout)        :: buffer
    integer        :: count, root, ierror
    type(MPI_COMM) :: comm
    call MPI_BCAST(buffer, count, MPI_LOGICAL, root, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_BCAST_LOGICAL_S_GPU
  subroutine D_MPI_BCAST_REAL32_R1_GPU(buffer, count, root, comm, ierror)
    implicit none
    real(real32), device, contiguous, intent(inout)   ::  buffer(:)
    integer        :: count, root, ierror
    type(MPI_COMM) :: comm
    call MPI_BCAST(buffer, count, MPI_REAL4, root, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_BCAST_REAL32_R1_GPU
  subroutine D_MPI_BCAST_REAL32_R2_GPU(buffer, count, root, comm, ierror)
    implicit none
    real(real32), device, contiguous, intent(inout)   ::  buffer(:,:)
    integer        :: count, root, ierror
    type(MPI_COMM) :: comm
    call MPI_BCAST(buffer, count, MPI_REAL4, root, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_BCAST_REAL32_R2_GPU
  subroutine D_MPI_BCAST_REAL64_R1_GPU(buffer, count, root, comm, ierror)
    implicit none
    real(real64), device, contiguous, intent(inout)  ::  buffer(:)
    integer        :: count, root, ierror
    type(MPI_COMM) :: comm
    call MPI_BCAST(buffer, count, MPI_REAL8, root, comm, ierror)
  end subroutine D_MPI_BCAST_REAL64_R1_GPU
  subroutine D_MPI_BCAST_REAL64_R2_GPU(buffer, count, root, comm, ierror)
    implicit none
    real(real64), device, contiguous, intent(inout)  ::  buffer(:,:)
    integer        :: count, root, ierror
    type(MPI_COMM) :: comm
    call MPI_BCAST(buffer, count, MPI_REAL8, root, comm, ierror)
  end subroutine D_MPI_BCAST_REAL64_R2_GPU
  subroutine D_MPI_BCAST_INT32_R2_GPU(buffer, count, root, comm, ierror)
    implicit none
    integer(int32), device, contiguous, intent(inout) ::  buffer(:,:)
    integer        :: count, root, ierror
    type(MPI_COMM) :: comm
    call MPI_BCAST(buffer, count, MPI_INTEGER4, root, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_BCAST_INT32_R2_GPU
  subroutine D_MPI_BCAST_INT32_R1_GPU(buffer, count, root, comm, ierror)
    implicit none
    integer(int32), device, contiguous, intent(inout) ::  buffer(:)
    integer        :: count, root, ierror
    type(MPI_COMM) :: comm
    call MPI_BCAST(buffer, count, MPI_INTEGER4, root, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_BCAST_INT32_R1_GPU
  subroutine D_MPI_BCAST_LOGICAL_R1_GPU(buffer, count, root, comm, ierror)
    implicit none
    logical, device, contiguous, intent(inout)        :: buffer(:)
    integer        :: count, root, ierror
    type(MPI_COMM) :: comm
    call MPI_BCAST(buffer, count, MPI_LOGICAL, root, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_BCAST_LOGICAL_R1_GPU

!>D_MPI_ALLREDUCE
  subroutine D_MPI_ALLREDUCE_REAL32_S_GPU(sendbuf, recvbuf, count, op, comm, ierror)
    implicit none
    real(real32), device, intent(inout)   :: sendbuf, recvbuf
    integer        :: count, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    call MPI_ALLREDUCE(sendbuf, recvbuf, count, MPI_REAL4, op, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ALLREDUCE_REAL32_S_GPU
  subroutine D_MPI_ALLREDUCE_REAL64_S_GPU(sendbuf, recvbuf, count, op, comm, ierror)
    implicit none
    real(real64), device, intent(inout)   :: sendbuf, recvbuf
    integer        :: count, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    call MPI_ALLREDUCE(sendbuf, recvbuf, count, MPI_REAL8, op, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ALLREDUCE_REAL64_S_GPU
  subroutine D_MPI_ALLREDUCE_INT32_S_GPU(sendbuf, recvbuf, count, op, comm, ierror)
    implicit none
    integer(int32), device, intent(inout) :: sendbuf, recvbuf
    integer        :: count, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    call MPI_ALLREDUCE(sendbuf, recvbuf, count, MPI_INTEGER4, op, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ALLREDUCE_INT32_S_GPU
  subroutine D_MPI_ALLREDUCE_REAL32_R1_GPU(sendbuf, recvbuf, count, op, comm, ierror)
    implicit none
    real(real32), device, contiguous, intent(inout)   :: sendbuf(:), recvbuf(:)
    integer        :: count, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    call MPI_ALLREDUCE(sendbuf, recvbuf, count, MPI_REAL4, op, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ALLREDUCE_REAL32_R1_GPU
  subroutine D_MPI_ALLREDUCE_REAL32_R2_GPU(sendbuf, recvbuf, count, op, comm, ierror)
    implicit none
    real(real32), device, contiguous, intent(inout)   :: sendbuf(:,:), recvbuf(:,:)
    integer        :: count, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    call MPI_ALLREDUCE(sendbuf, recvbuf, count, MPI_REAL4, op, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ALLREDUCE_REAL32_R2_GPU
  subroutine D_MPI_ALLREDUCE_REAL32_R3_GPU(sendbuf, recvbuf, count, op, comm, ierror)
    implicit none
    real(real32), device, contiguous, intent(inout)   :: sendbuf(:,:,:), recvbuf(:,:,:)
    integer        :: count, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    call MPI_ALLREDUCE(sendbuf,recvbuf, count, MPI_REAL4, op, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ALLREDUCE_REAL32_R3_GPU
  subroutine D_MPI_ALLREDUCE_REAL64_R1_GPU(sendbuf, recvbuf, count, op, comm, ierror)
    implicit none
    real(real64), device, contiguous, intent(inout)   :: sendbuf(:), recvbuf(:)
    integer        :: count, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    call MPI_ALLREDUCE(sendbuf, recvbuf, count, MPI_REAL8, op, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ALLREDUCE_REAL64_R1_GPU
  subroutine D_MPI_ALLREDUCE_REAL64_R2_GPU(sendbuf, recvbuf, count, op, comm, ierror)
    implicit none
    real(real64), device, contiguous, intent(inout)   :: sendbuf(:,:), recvbuf(:,:)
    integer        :: count, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    call MPI_ALLREDUCE(sendbuf, recvbuf, count, MPI_REAL8, op, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ALLREDUCE_REAL64_R2_GPU
  subroutine D_MPI_ALLREDUCE_REAL64_R3_GPU(sendbuf, recvbuf, count, op, comm, ierror)
    implicit none
    real(real64), device, contiguous, intent(inout)   :: sendbuf(:,:,:), recvbuf(:,:,:)
    integer        :: count, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    call MPI_ALLREDUCE(sendbuf, recvbuf, count, MPI_REAL8, op, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ALLREDUCE_REAL64_R3_GPU
  subroutine D_MPI_ALLREDUCE_INT32_R2_GPU(sendbuf, recvbuf, count, op, comm, ierror)
    implicit none
    integer(int32), device, contiguous, intent(inout) :: sendbuf(:,:), recvbuf(:,:)
    integer        :: count, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    call MPI_ALLREDUCE(sendbuf, recvbuf, count, MPI_INTEGER4, op, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ALLREDUCE_INT32_R2_GPU
  subroutine D_MPI_ALLREDUCE_INT32_R1_GPU(sendbuf, recvbuf, count, op, comm, ierror)
    implicit none
    integer(int32), device, contiguous, intent(inout) :: sendbuf(:), recvbuf(:)
    integer        :: count, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    call MPI_ALLREDUCE(sendbuf, recvbuf, count, MPI_INTEGER4, op, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ALLREDUCE_INT32_R1_GPU
  subroutine D_MPI_ALLREDUCE_REAL32_IP_GPU(recvbuf, count, op, comm, ierror)
    implicit none
    real(real32), device, contiguous, intent(inout)   :: recvbuf(:)
    integer        :: count, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    call MPI_ALLREDUCE(MPI_IN_PLACE, recvbuf, count, MPI_REAL4, op, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ALLREDUCE_REAL32_IP_GPU
  subroutine D_MPI_ALLREDUCE_REAL64_IP_GPU(recvbuf, count, op, comm, ierror)
    implicit none
    real(real64), device, contiguous, intent(inout)   :: recvbuf(:)
    integer        :: count, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    call MPI_ALLREDUCE(MPI_IN_PLACE, recvbuf, count, MPI_REAL8, op, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ALLREDUCE_REAL64_IP_GPU

!>D_MPI_ALLTOALL
  subroutine D_MPI_ALLTOALL_REAL32_R1_GPU(sendbuf, sendcount, recvbuf, recvcount, comm, ierror)
    implicit none
    real(real32), device, contiguous, intent(inout)   :: sendbuf(:), recvbuf(:)
    integer        :: sendcount, recvcount, ierror
    type(MPI_COMM) :: comm
    call MPI_ALLTOALL(sendbuf, sendcount, MPI_REAL4, recvbuf, recvcount, MPI_REAL4, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ALLTOALL_REAL32_R1_GPU
  subroutine D_MPI_ALLTOALL_REAL64_R1_GPU(sendbuf, sendcount, recvbuf, recvcount, comm, ierror)
    implicit none
    real(real64), device, contiguous, intent(inout)   :: sendbuf(:), recvbuf(:)
    integer        :: sendcount, recvcount, ierror
    type(MPI_COMM) :: comm
    call MPI_ALLTOALL(sendbuf, sendcount, MPI_REAL8, recvbuf, recvcount, MPI_REAL8, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_ALLTOALL_REAL64_R1_GPU

!>D_MPI_REDUCE
  subroutine D_MPI_REDUCE_REAL32_R1_GPU(sendbuf, recvbuf, count, op, root, comm, ierror)
    implicit none
    real(real32), device, contiguous, intent(inout)   :: sendbuf(:), recvbuf(:)
    integer        :: count, root, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    call MPI_REDUCE(sendbuf, recvbuf, count, MPI_REAL4, op, root, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_REDUCE_REAL32_R1_GPU
  subroutine D_MPI_REDUCE_REAL32_R2_GPU(sendbuf, recvbuf, count, op, root, comm, ierror)
    implicit none
    real(real32), device, contiguous, intent(inout)   :: sendbuf(:,:), recvbuf(:,:)
    integer        :: count, root, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    call MPI_REDUCE(sendbuf, recvbuf, count, MPI_REAL4, op, root, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_REDUCE_REAL32_R2_GPU
  subroutine D_MPI_REDUCE_REAL32_R3_GPU(sendbuf, recvbuf, count, op, root, comm, ierror)
    implicit none
    real(real32), device, contiguous, intent(inout)   :: sendbuf(:,:,:), recvbuf(:,:,:)
    integer        :: count, root, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    call MPI_REDUCE(sendbuf, recvbuf, count, MPI_REAL4, op, root, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_REDUCE_REAL32_R3_GPU
  subroutine D_MPI_REDUCE_REAL64_R1_GPU(sendbuf, recvbuf, count, op, root, comm, ierror)
    implicit none
    real(real64), device, contiguous, intent(inout)   :: sendbuf(:), recvbuf(:)
    integer        :: count, root, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    call MPI_REDUCE(sendbuf, recvbuf, count, MPI_REAL8, op, root, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_REDUCE_REAL64_R1_GPU
  subroutine D_MPI_REDUCE_REAL64_R2_GPU(sendbuf, recvbuf, count, op, root, comm, ierror)
    implicit none
    real(real64), device, contiguous, intent(inout)   :: sendbuf(:,:), recvbuf(:,:)
    integer        :: count, root, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    call MPI_REDUCE(sendbuf, recvbuf, count, MPI_REAL8, op, root, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_REDUCE_REAL64_R2_GPU
  subroutine D_MPI_REDUCE_REAL64_R3_GPU(sendbuf, recvbuf, count, op, root, comm, ierror)
    implicit none
    real(real64), device, contiguous, intent(inout)   :: sendbuf(:,:,:), recvbuf(:,:,:)
    integer        :: count, root, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    call MPI_REDUCE(sendbuf, recvbuf, count, MPI_REAL8, op, root, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_REDUCE_REAL64_R3_GPU
  subroutine D_MPI_REDUCE_REAL32_IP_R1_GPU(recvbuf, count, op, root, comm, ierror)
    implicit none
    real(real32), device, contiguous, intent(inout)   :: recvbuf(:)
    integer        :: count, root, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    call MPI_REDUCE(MPI_IN_PLACE, recvbuf, count, MPI_REAL4, op, root, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_REDUCE_REAL32_IP_R1_GPU
  subroutine D_MPI_REDUCE_REAL32_IP_R2_GPU(recvbuf, count, op, root, comm, ierror)
    implicit none
    real(real32), device, contiguous, intent(inout)   :: recvbuf(:,:)
    integer        :: count, root, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    call MPI_REDUCE(MPI_IN_PLACE, recvbuf, count, MPI_REAL4, op, root, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_REDUCE_REAL32_IP_R2_GPU
  subroutine D_MPI_REDUCE_REAL64_IP_R1_GPU(recvbuf, count, op, root, comm, ierror)
    implicit none
    real(real64), device, contiguous, intent(inout)  :: recvbuf(:)
    integer        :: count, root, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    call MPI_REDUCE(MPI_IN_PLACE, recvbuf, count, MPI_REAL8, op, root, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_REDUCE_REAL64_IP_R1_GPU
  subroutine D_MPI_REDUCE_REAL64_IP_R2_GPU(recvbuf, count, op, root, comm, ierror)
    implicit none
    real(real64), device, contiguous, intent(inout)  :: recvbuf(:,:)
    integer        :: count, root, ierror
    type(MPI_OP)   :: op
    type(MPI_COMM) :: comm
    call MPI_REDUCE(MPI_IN_PLACE, recvbuf, count, MPI_REAL8, op, root, comm, ierror)
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_REDUCE_REAL64_IP_R2_GPU

!>D_MPI_GATHER
  subroutine D_MPI_GATHER_REAL32_R1_GPU(sendbuf, sendcount, recvbuf, recvcount, root, comm, ierror)
    implicit none
    real(real32), device, contiguous, intent(inout)   :: sendbuf(:), recvbuf(:)
    integer        :: sendcount, recvcount, root, ierror
    type(MPI_COMM) :: comm
    call MPI_GATHER( sendbuf, sendcount, MPI_REAL4 &
                   , recvbuf, recvcount, MPI_REAL4 &
                   , root, comm, ierror )
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_GATHER_REAL32_R1_GPU
  subroutine D_MPI_GATHER_REAL64_R1_GPU(sendbuf, sendcount, recvbuf, recvcount, root, comm, ierror)
    implicit none
    real(real64), device, contiguous, intent(inout)   :: sendbuf(:), recvbuf(:)
    integer        :: sendcount, recvcount, root, ierror
    type(MPI_COMM) :: comm
    call MPI_GATHER( sendbuf, sendcount, MPI_REAL8 &
                   , recvbuf, recvcount, MPI_REAL8 &
                   , root, comm, ierror )
    if (ierror /= MPI_SUCCESS) call abort
  end subroutine D_MPI_GATHER_REAL64_R1_GPU

#endif

end module
