! -
!
! SPDX-FileCopyrightText: Copyright (c) 2022 Pedro Costa. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
!
! a simple timer, see https://github.com/p-costa/first-timer
!
module modtimer

  use, intrinsic :: iso_fortran_env, only: dp => real64

  use modmpi,          only: myid, d_mpi_bcast, comm3d, nprocs, &
                             d_mpi_allreduce, mpi_min, mpi_max, mpi_sum, &
                             mpi_wtime
  use modglobal,       only: checknamelisterror, ifnamopt, fname_options
  use fortran_support, only: nnml_output, nout, t_table, add_table_column, &
                             set_table_entry, print_table, initialize_table, &
                             find_next_free_unit, real2string, int2string, &
                             finish
#if defined(USE_NVTX)
  use modnvtx
#endif

  implicit none

  private

  public :: timer_read_namelist
  public :: timer_init
  public :: timer_tic
  public :: timer_toc
  public :: timer_cleanup
  public :: ltimer
  public :: output_timings

  character(len=*), parameter :: modname = "modtimer"

  logical, parameter :: GPU_DEFAULT_SYNC = .true.
  integer, parameter :: max_name_len = 50

  character(max_name_len), allocatable :: timer_names(:)
  integer,                 allocatable :: timer_counts(:)
  integer,                 allocatable :: timer_counter(:)
  real(dp),                allocatable :: timer_tictoc(:)
  real(dp),                allocatable :: timer_elapsed_acc(:)
  real(dp),                allocatable :: timer_elapsed_min(:)
  real(dp),                allocatable :: timer_elapsed_max(:)
  logical,                 allocatable :: timer_is_nvtx(:)

  integer :: ntimers = 0            !< Number of timers.

  logical :: ltimer = .false.       !< Switch for enabling/disabling timings.
  logical :: ltimer_print = .true.  !< Switch for printing timing results to std out.
  logical :: ltimer_write = .false. !< Switch for writing timing results to a csv file.
  logical :: lverbose = .false.     !< Switch for printing per-rank statistics.

contains

  subroutine timer_read_namelist(nml_filename)

    character(len=*), intent(in) :: nml_filename

    integer :: ierr

    namelist /timer/ ltimer, ltimer_print, ltimer_write

    if (myid == 0) then
      open(ifnamopt, file=nml_filename, status="old", iostat=ierr)
      read(ifnamopt, timer, iostat=ierr)
      call checknamelisterror(ierr, ifnamopt , "timer")
      write(nnml_output, timer)
      close(ifnamopt)
    end if

    call D_MPI_BCAST(ltimer, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(ltimer_print, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(ltimer_write, 1, 0, comm3d, ierr)

  end subroutine timer_read_namelist

  subroutine timer_init()

    allocate(timer_names(0), &
             timer_counts(0), &
             timer_counter(0), &
             timer_tictoc(0), &
             timer_elapsed_acc(0), &
             timer_elapsed_min(0), &
             timer_elapsed_max(0), &
             timer_is_nvtx(0))

  end subroutine timer_init

  subroutine output_timings()

    character(len=*), parameter :: routine = modname//"/output_timings"

    character(len=64), parameter :: column_names(12) = [character(len=64) :: &
      "label", &
      "calls", &
      "total elapsed time [s]", &
      "avg elapsed time per call [s]", &
      "min time per call [s]", &
      "max time per call [s]", &
      "avg elapsed time per call (minimum across tasks) [s]", &
      "avg elapsed time per call (maximum across tasks) [s]", &
      "min time per call (minimum across tasks) [s]", &
      "min time per call (maximum across tasks) [s]", &
      "max time per call (minimum across tasks) [s]", &
      "max time per call (maximum across tasks) [s]" &
    ]

    integer :: itable_out
    integer :: ierr
    integer :: i

    real(dp), allocatable :: timing_results_acc(:,:)
    real(dp), allocatable :: timing_results_min(:,:)
    real(dp), allocatable :: timing_results_max(:,:)

    type(t_table) :: table

    if (.not. ltimer) return

    ! Accumulate timing results across all MPI tasks
    allocate(timing_results_acc(ntimers,3), &
             timing_results_min(ntimers,3), &
             timing_results_max(ntimers,3))

    call D_MPI_ALLREDUCE(timer_elapsed_acc(:), timing_results_acc(:,1), ntimers, MPI_MIN, comm3d, &
                         ierr)
    call D_MPI_ALLREDUCE(timer_elapsed_acc(:), timing_results_acc(:,2), ntimers, MPI_MAX, comm3d, &
                         ierr)
    call D_MPI_ALLREDUCE(timer_elapsed_acc(:), timing_results_acc(:,3), ntimers, MPI_SUM, comm3d, &
                         ierr)
    timing_results_acc(:,3) = timing_results_acc(:,3) / nprocs

    call D_MPI_ALLREDUCE(timer_elapsed_min(:), timing_results_min(:,1), ntimers, MPI_MIN, comm3d, &
                         ierr)
    call D_MPI_ALLREDUCE(timer_elapsed_min(:), timing_results_min(:,2), ntimers, MPI_MAX, comm3d, &
                         ierr)
    call D_MPI_ALLREDUCE(timer_elapsed_min(:), timing_results_min(:,3), ntimers, MPI_SUM, comm3d, &
                         ierr)
    timing_results_min(:,3) = timing_results_min(:,3) / nprocs

    call D_MPI_ALLREDUCE(timer_elapsed_max(:), timing_results_max(:,1), ntimers, MPI_MIN, comm3d, &
                         ierr)
    call D_MPI_ALLREDUCE(timer_elapsed_max(:), timing_results_max(:,2), ntimers, MPI_MAX, comm3d, &
                         ierr)
    call D_MPI_ALLREDUCE(timer_elapsed_max(:), timing_results_max(:,3), ntimers, MPI_SUM, comm3d, &
                         ierr)
    timing_results_max(:,3) = timing_results_max(:,3) / nprocs

    ! Format the timing results as a nice table
    call initialize_table(table)

    do i = 1, ntimers
      call set_table_entry(table, i, column_names(1), timer_names(i))
      call set_table_entry(table, i, column_names(2), int2string(timer_counts(i)))
      call set_table_entry(table, i, column_names(3), real2string(timing_results_acc(i,3)))
      call set_table_entry(table, i, column_names(4), &
                           real2string(timing_results_acc(i,3) / timer_counts(i)))
      call set_table_entry(table, i, column_names(5), real2string(timing_results_min(i,3)))
      call set_table_entry(table, i, column_names(6), real2string(timing_results_max(i,3)))
      if (lverbose) then
        call set_table_entry(table, i, column_names(7), real2string(timing_results_acc(i,1)))
        call set_table_entry(table, i, column_names(8), real2string(timing_results_acc(i,2)))
        call set_table_entry(table, i, column_names(9), real2string(timing_results_min(i,1)))
        call set_table_entry(table, i, column_names(10), real2string(timing_results_min(i,2)))
        call set_table_entry(table, i, column_names(11), real2string(timing_results_max(i,1)))
        call set_table_entry(table, i, column_names(12), real2string(timing_results_max(i,2)))
      end if
      end do

    itable_out = find_next_free_unit(10, 20)

    ! Write to file
    open(unit=itable_out, file="timings.txt", status="replace", action="write", iostat=ierr)

    if (ierr /= 0) then
      call finish(routine, "could not open timings file for writing")
    end if

    call print_table(table, opt_dstfile=itable_out)

  end subroutine output_timings

  subroutine timer_tic(timer_name,nvtx_id_fix,nvtx_color,nvtx_id_inc,nvtx_gpu_stream)

    character(len=*), intent(in) :: timer_name
    integer         , intent(in   ), optional :: nvtx_id_fix     ! if <= 0, only label and no color
    character(len=1), intent(in   ), optional :: nvtx_color      ! g/b/y/m/c/r/w following matplotlib's convention
    integer         , intent(inout), optional :: nvtx_id_inc     ! to increment the id, e.g.: call timer_tic(name,nvtx_id_inc=i_nvtx)
    integer         , intent(in   ), optional :: nvtx_gpu_stream ! to optionally sync host/device over a stream/queue (asynchronous if < 0)
    integer :: idx,nvtx_id
    logical :: is_nvtx,is_gpu_sync

    if (.not. ltimer) return

    !
    idx = findloc(timer_names, timer_name, dim=1)
    if (idx <= 0) then
      ntimers = ntimers + 1
      call concatenate_c(timer_names,timer_name)
      timer_counts      = [timer_counts     ,0          ]
      timer_counter     = [timer_counter    ,0          ]
      timer_tictoc      = [timer_tictoc     ,0._dp      ]
      timer_elapsed_acc = [timer_elapsed_acc,0._dp      ]
      timer_elapsed_min = [timer_elapsed_min,huge(0._dp)]
      timer_elapsed_max = [timer_elapsed_max,tiny(0._dp)]
      timer_is_nvtx     = [timer_is_nvtx    ,.false.    ]
      idx = ntimers
    end if
    timer_counter(idx)     = timer_counter(idx) + 1
    timer_tictoc(idx) = MPI_WTIME()
#if defined(USE_NVTX)
    is_nvtx = .false.
    if(     present(nvtx_id_inc)) then
      nvtx_id = nvtx_id_inc
      if(nvtx_id == huge(1)) nvtx_id_inc = 0 ! avoid overflow
      nvtx_id_inc = nvtx_id_inc + 1
      is_nvtx = .true.
    else if(present(nvtx_id_fix)) then
      nvtx_id = nvtx_id_fix
      is_nvtx = .true.
    else if(present(nvtx_color )) then
      is_nvtx = .true.
    end if
    if(is_nvtx) then
      is_gpu_sync = GPU_DEFAULT_SYNC
      if(present(nvtx_gpu_stream)) then
        if(nvtx_gpu_stream < 0) then
          is_gpu_sync = .false.
        end if
      end if
      if(is_gpu_sync) then
        if(.not.present(nvtx_gpu_stream)) then
          !$acc wait
        else
          !$acc wait(nvtx_gpu_stream)
        end if
      end if
      if(     present(nvtx_color)) then
        call nvtxStartRange(trim(timer_name),color=nvtx_color)
      else if(nvtx_id > 0        ) then
          call nvtxStartRange(trim(timer_name),id=nvtx_id)
      else
        call nvtxStartRange(trim(timer_name))
      end if
      timer_is_nvtx(idx) = .true.
    end if
#endif
  end subroutine timer_tic
  subroutine timer_toc(timer_name,nvtx_gpu_stream,ierror)
    character(*), intent(in) :: timer_name
    integer, intent(in), optional :: nvtx_gpu_stream
    integer, intent(out), optional :: ierror
    integer :: idx
    logical :: is_gpu_sync

    if (.not. ltimer) return
    
    if(present(ierror)) ierror = 0
    idx = findloc(timer_names, timer_name, dim=1)
    if (idx > 0) then
      timer_tictoc(idx)      = MPI_WTIME() - timer_tictoc(idx)
      timer_elapsed_acc(idx) =    (timer_elapsed_acc(idx)+timer_tictoc(idx))
      timer_elapsed_min(idx) = min(timer_elapsed_min(idx),timer_tictoc(idx))
      timer_elapsed_max(idx) = max(timer_elapsed_max(idx),timer_tictoc(idx))
      timer_counts(idx)      = timer_counts(idx) + 1
      timer_counter(idx)     = timer_counter(idx) - 1
      if(timer_is_nvtx(idx)) then
        is_gpu_sync = GPU_DEFAULT_SYNC
        if(present(nvtx_gpu_stream)) then
          if(nvtx_gpu_stream < 0) then
            is_gpu_sync = .false.
          end if
        end if
        if(is_gpu_sync) then
          if(.not.present(nvtx_gpu_stream)) then
            !$acc wait
          else
            !$acc wait(nvtx_gpu_stream)
          end if
        end if
#if defined(USE_NVTX)
        call nvtxEndRange
#endif
      end if
    else
      if(present(ierror)) ierror = 1
    end if

  end subroutine timer_toc

  subroutine timer_cleanup

    integer :: i
    do i = 1,ntimers
      if (timer_counter(i) .ne. 0) then
        print*,'WARNING: malformed timer: ', timer_names(i)
      end if
    end do
    if (allocated(timer_names)) then
      deallocate(timer_names,timer_counts,timer_counter,timer_elapsed_acc,timer_elapsed_min,timer_elapsed_max)
    end if

  end subroutine timer_cleanup

  subroutine concatenate_c(arr,val)

    character(*), intent(inout), allocatable, dimension(:) :: arr
    character(*), intent(in   ) :: val
    character(:), allocatable, dimension(:) :: arr_tmp
    integer :: n
    n = size(arr)
    allocate(arr_tmp,source=arr)
    deallocate(arr); allocate(arr(n+1))
    arr(1:n) = arr_tmp(:); arr(n+1) = val

  end subroutine concatenate_c

end module modtimer
