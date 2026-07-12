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
                             finish, warning
#if defined(USE_CUDA)
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

  integer, parameter :: max_name_len = 50

  character(max_name_len), allocatable :: timer_names(:)       !< Names of the timers.
  integer,                 allocatable :: timer_counts(:)      !< Number of calls to each timer.
  integer,                 allocatable :: timer_counter(:)     !< Nesting level of timer.
  real(dp),                allocatable :: timer_tictoc(:)      !< Elapsed time for the current tic/toc call [s].
  real(dp),                allocatable :: timer_elapsed_acc(:) !< Total elapsed time accumulated for each timer [s].
  real(dp),                allocatable :: timer_elapsed_min(:) !< Minimum elapsed time for each timer [s].
  real(dp),                allocatable :: timer_elapsed_max(:) !< Maximum elapsed time for each timer [s].

  integer :: ntimers = 0 !< Number of timers.
  integer :: level = 0   !< Current level of nested timers.

  logical, protected :: ltimer = .false.   !< Switch for enabling/disabling timings.
  logical            :: lverbose = .false. !< Switch for printing per-rank statistics.
  logical            :: lnvtx = .true.     !< Switch for enabling NVTX regions.
  integer            :: max_level = 0      !< Maximum level of nested timers (if < 1, then no maximum).

contains

  !> Read the timer namelist.
  subroutine timer_read_namelist(nml_filename)

    character(len=*), intent(in) :: nml_filename !< Name of the namelist file.

    integer :: ierr

    namelist /timer/ ltimer, lnvtx, max_level

    if (myid == 0) then
      open(ifnamopt, file=nml_filename, status="old", iostat=ierr)
      read(ifnamopt, timer, iostat=ierr)
      call checknamelisterror(ierr, ifnamopt , "timer")
      write(nnml_output, timer)
      close(ifnamopt)
    end if

    call D_MPI_BCAST(ltimer, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(lnvtx, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(max_level, 1, 0, comm3d, ierr)

  end subroutine timer_read_namelist

  !> Initialize the timer module and allocate arrays.
  subroutine timer_init()

    allocate(timer_names(0), &
             timer_counts(0), &
             timer_counter(0), &
             timer_tictoc(0), &
             timer_elapsed_acc(0), &
             timer_elapsed_min(0), &
             timer_elapsed_max(0))

  end subroutine timer_init

  !> Accumulate timing results across all MPI tasks and print to file.
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

  !> Start the timer and record the start time.
  subroutine timer_tic(timer_name, opt_nvtx_id)

    character(len=*), intent(in) :: timer_name !< Name of the timer to start.

    integer, intent(in), optional :: opt_nvtx_id !< Optional NVTX range ID.

    integer :: idx
    integer :: nvtx_id
    real    :: wtime

    if (.not. ltimer) return
    if (myid == 0) then
       wtime = MPI_Wtime()
       write(*, *) wtime, '*', timer_name
    end if
    level = level + 1

    if (level > max_level .and. max_level > 0) return

    nvtx_id = -1
    if (present(opt_nvtx_id)) nvtx_id = opt_nvtx_id

    idx = findloc(timer_names, timer_name, dim=1)

    if (idx <= 0) then
      ntimers = ntimers + 1
      call concatenate_c(timer_names,timer_name)
      timer_counts      = [timer_counts,      0          ]
      timer_counter     = [timer_counter,     0          ]
      timer_tictoc      = [timer_tictoc,      0._dp      ]
      timer_elapsed_acc = [timer_elapsed_acc, 0._dp      ]
      timer_elapsed_min = [timer_elapsed_min, huge(0._dp)]
      timer_elapsed_max = [timer_elapsed_max, tiny(0._dp)]
      idx = ntimers
    end if
    timer_counter(idx) = timer_counter(idx) + 1
    timer_tictoc(idx) = MPI_WTIME()

    !$acc wait

#if defined(USE_CUDA)
    if (lnvtx) then
      if (nvtx_id > 0) then
        call nvtxStartRange(trim(timer_name), id=nvtx_id)
      else
        call nvtxStartRange(trim(timer_name))
      end if
    end if
#endif

  end subroutine timer_tic

  !> Stop the timer and accumulate the elapsed time.
  subroutine timer_toc(timer_name)

    character(len=*), intent(in) :: timer_name !< Name of the timer to stop.

    character(len=*), parameter :: routine = modname//"/timer_toc"

    integer :: idx

    if (.not. ltimer) return

    if (level > max_level .and. max_level > 0) then
      level = level - 1
      return
    else
      level = level - 1

      idx = findloc(timer_names, timer_name, dim=1)

      if (idx > 0) then
        timer_tictoc(idx)      = MPI_WTIME() - timer_tictoc(idx)
        timer_elapsed_acc(idx) =    (timer_elapsed_acc(idx) + timer_tictoc(idx))
        timer_elapsed_min(idx) = min(timer_elapsed_min(idx), timer_tictoc(idx))
        timer_elapsed_max(idx) = max(timer_elapsed_max(idx), timer_tictoc(idx))
        timer_counts(idx)      = timer_counts(idx) + 1
        timer_counter(idx)     = timer_counter(idx) - 1

        !$acc wait
#if defined(USE_CUDA)
        if (lnvtx) call nvtxEndRange
#endif
      else
        call finish(routine, "timer " // trim(timer_name) // " not found")
      end if
    end if

  end subroutine timer_toc

  !> Check for unbalanced tic/toc calls and deallocate arrays.
  subroutine timer_cleanup()

    character(len=*), parameter :: routine = modname//"/timer_cleanup"

    integer :: i

    do i = 1,ntimers
      if (timer_counter(i) /= 0) then
        call warning(routine, "timer " // trim(timer_names(i)) // " has unbalanced tic/toc calls")
      end if
    end do
    if (allocated(timer_names)) then
      deallocate(timer_names, timer_counts, timer_counter, timer_elapsed_acc, timer_elapsed_min, &
                 timer_elapsed_max)
    end if

  end subroutine timer_cleanup

  !> Concatenate a character array with a new value.
  subroutine concatenate_c(arr, val)

    character(len=max_name_len), intent(inout), allocatable :: arr(:) !< Character array to concatenate to.

    character(len=*),            intent(in)    :: val !< New value to concatenate.

    character(len=max_name_len), allocatable :: arr_tmp(:)

    integer :: n

    n = size(arr)
    allocate(arr_tmp, source=arr)
    deallocate(arr)
    allocate(arr(n + 1))
    arr(1:n) = arr_tmp(:)
    arr(n+1) = val

  end subroutine concatenate_c

end module modtimer
