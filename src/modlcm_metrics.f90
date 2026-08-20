!> End-of-experiment reporting for optional standalone LCM metrics.
#ifdef USE_LCM
module modlcm_metrics
  use iso_fortran_env, only : int64, real64
  use lcm_host_interface, only : lcm_get_particle_exchange_metrics,       &
                                 lcm_get_particle_exchange_timing_summary,&
                                 lcm_get_particle_storage_metrics,        &
                                 lcm_particle_exchange_metrics_t,         &
                                 lcm_particle_exchange_timing_summary_t,  &
                                 lcm_particle_storage_metrics_t,          &
                                 lcm_rank_timing_summary_t,               &
                                 lcm_reset_particle_exchange_metrics,     &
                                 lcm_set_particle_exchange_metrics_enabled
  use modglobal, only : cexpnr
  use modmpi, only : checkmpierror, comm3d, MPI_Datatype, MPI_MAX,        &
                     MPI_MIN, MPI_REDUCE, MPI_SUM, MPI_TYPECLASS_INTEGER,&
                     MPI_TYPE_MATCH_SIZE, myid, nprocs

  implicit none

  private

  integer, parameter :: exchange_metric_count = 23
  integer, parameter :: storage_metric_count = 9
  integer, parameter :: integer_metric_count = exchange_metric_count +   &
                                                     storage_metric_count

  logical, save :: metrics_are_active = .false.

  public :: finish_lcm_metrics
  public :: start_lcm_metrics

contains

  !> Reset and enable instrumentation only for explicitly requested runs.
  subroutine start_lcm_metrics(enabled)
    logical, intent(in) :: enabled

    metrics_are_active = enabled
    call lcm_set_particle_exchange_metrics_enabled(enabled)
    if (enabled) call lcm_reset_particle_exchange_metrics()
  end subroutine start_lcm_metrics

  !> Collect rank-local metrics and let rank zero write one experiment CSV.
  subroutine finish_lcm_metrics()
    type(lcm_particle_exchange_metrics_t) :: exchange_metrics
    type(lcm_particle_exchange_timing_summary_t) :: timing_summary
    type(lcm_particle_storage_metrics_t) :: storage_metrics
    type(MPI_Datatype) :: int64_datatype
    integer(int64) :: local_values(integer_metric_count)
    integer(int64) :: maximum_values(integer_metric_count)
    integer(int64) :: minimum_values(integer_metric_count)
    integer(int64) :: summed_values(integer_metric_count)
    integer :: ierr

    if (.not. metrics_are_active) return

    ! This summary is collective. Every DALES rank reaches this routine from
    ! exitmicrophysics before comm3d is released by exitmpi.
    call lcm_get_particle_exchange_timing_summary(timing_summary)
    call lcm_get_particle_exchange_metrics(exchange_metrics)
    call lcm_get_particle_storage_metrics(storage_metrics)
    call pack_integer_metrics(exchange_metrics, storage_metrics, local_values)

    call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_INTEGER,                     &
                             storage_size(local_values(1)) / 8,          &
                             int64_datatype, ierr)
    call checkmpierror(ierr, 'LCM metrics integer datatype')
    call MPI_REDUCE(local_values, summed_values, integer_metric_count,   &
                    int64_datatype, MPI_SUM, 0, comm3d, ierr)
    call checkmpierror(ierr, 'LCM metrics sums')
    call MPI_REDUCE(local_values, minimum_values, integer_metric_count,  &
                    int64_datatype, MPI_MIN, 0, comm3d, ierr)
    call checkmpierror(ierr, 'LCM metrics minima')
    call MPI_REDUCE(local_values, maximum_values, integer_metric_count,  &
                    int64_datatype, MPI_MAX, 0, comm3d, ierr)
    call checkmpierror(ierr, 'LCM metrics maxima')

    if (myid == 0) then
      call write_metrics_csv(timing_summary, summed_values, minimum_values, &
                             maximum_values)
    end if

    call lcm_set_particle_exchange_metrics_enabled(.false.)
    metrics_are_active = .false.
  end subroutine finish_lcm_metrics

  subroutine pack_integer_metrics(exchange, storage, values)
    type(lcm_particle_exchange_metrics_t), intent(in) :: exchange
    type(lcm_particle_storage_metrics_t), intent(in) :: storage
    integer(int64), intent(out) :: values(integer_metric_count)

    values(1) = exchange%x_particles_examined
    values(2) = exchange%y_particles_examined
    values(3) = exchange%x_particles_outgoing
    values(4) = exchange%y_particles_outgoing
    values(5) = exchange%x_count_exchange_calls
    values(6) = exchange%y_count_exchange_calls
    values(7) = exchange%x_payload_exchange_calls
    values(8) = exchange%y_payload_exchange_calls
    values(9) = exchange%east_particles_sent
    values(10) = exchange%west_particles_sent
    values(11) = exchange%south_particles_sent
    values(12) = exchange%north_particles_sent
    values(13) = exchange%east_payload_bytes_sent
    values(14) = exchange%west_payload_bytes_sent
    values(15) = exchange%south_payload_bytes_sent
    values(16) = exchange%north_payload_bytes_sent
    values(17) = exchange%x_received_insertion_calls
    values(18) = exchange%y_received_insertion_calls
    values(19) = exchange%x_particles_inserted
    values(20) = exchange%y_particles_inserted
    values(21) = exchange%ownership_update_calls
    values(22) = exchange%ownership_particle_examinations
    values(23) = exchange%ownership_particles_moved

    values(24) = storage%cell_count
    values(25) = storage%allocated_cell_count
    values(26) = storage%empty_allocated_cell_count
    values(27) = storage%active_particle_count
    values(28) = storage%reserved_particle_capacity
    values(29) = storage%unused_reserved_particle_capacity
    values(30) = storage%particle_payload_bytes_per_slot
    values(31) = storage%active_particle_payload_bytes
    values(32) = storage%reserved_particle_payload_bytes
  end subroutine pack_integer_metrics

  subroutine write_metrics_csv(timing, sums, minima, maxima)
    type(lcm_particle_exchange_timing_summary_t), intent(in) :: timing
    integer(int64), intent(in) :: sums(integer_metric_count)
    integer(int64), intent(in) :: minima(integer_metric_count)
    integer(int64), intent(in) :: maxima(integer_metric_count)
    character(len=64), parameter :: exchange_names(exchange_metric_count) = [ &
      character(len=64) ::                                                   &
      'x_particles_examined', 'y_particles_examined',                       &
      'x_particles_outgoing', 'y_particles_outgoing',                       &
      'x_count_exchange_calls', 'y_count_exchange_calls',                   &
      'x_payload_exchange_calls', 'y_payload_exchange_calls',               &
      'east_particles_sent', 'west_particles_sent',                         &
      'south_particles_sent', 'north_particles_sent',                       &
      'east_payload_bytes_sent', 'west_payload_bytes_sent',                 &
      'south_payload_bytes_sent', 'north_payload_bytes_sent',               &
      'x_received_insertion_calls', 'y_received_insertion_calls',           &
      'x_particles_inserted', 'y_particles_inserted',                       &
      'ownership_update_calls', 'ownership_particle_examinations',          &
      'ownership_particles_moved' ]
    character(len=16), parameter :: exchange_units(exchange_metric_count) = [ &
      character(len=16) ::                                                   &
      'particles', 'particles', 'particles', 'particles',                   &
      'calls', 'calls', 'calls', 'calls',                                   &
      'particles', 'particles', 'particles', 'particles',                   &
      'bytes', 'bytes', 'bytes', 'bytes',                                   &
      'calls', 'calls', 'particles', 'particles',                           &
      'calls', 'particles', 'particles' ]
    character(len=64), parameter :: storage_names(storage_metric_count) = [ &
      character(len=64) ::                                                  &
      'cell_count', 'allocated_cell_count', 'empty_allocated_cell_count',   &
      'active_particle_count', 'reserved_particle_capacity',                &
      'unused_reserved_particle_capacity', 'particle_payload_bytes_per_slot',&
      'active_particle_payload_bytes', 'reserved_particle_payload_bytes' ]
    character(len=16), parameter :: storage_units(storage_metric_count) = [ &
      character(len=16) ::                                                  &
      'cells', 'cells', 'cells', 'particles', 'slots', 'slots',             &
      'bytes', 'bytes', 'bytes' ]
    character(len=128) :: filename
    integer :: i
    integer :: iostat
    integer :: unit

    filename = 'lcm_particle_metrics.' // trim(cexpnr) // '.csv'
    open(newunit=unit, file=trim(filename), status='replace', action='write', &
         iostat=iostat)
    if (iostat /= 0) error stop 'DALES failed to open the LCM metrics CSV'

    write(unit, '(a)') 'experiment,category,metric,unit,aggregation,value,rank_count'
    call write_timing_rows(unit, timing)

    do i = 1, exchange_metric_count
      call write_integer_rows(unit, 'exchange_counter', exchange_names(i), &
                              exchange_units(i), sums(i), minima(i),        &
                              maxima(i))
    end do
    do i = 1, storage_metric_count
      call write_integer_rows(unit, 'storage', storage_names(i),           &
                              storage_units(i), sums(exchange_metric_count+i),&
                              minima(exchange_metric_count+i),             &
                              maxima(exchange_metric_count+i))
    end do

    close(unit)
    write(6, '(a)') 'LCM particle metrics written to ' // trim(filename)
  end subroutine write_metrics_csv

  subroutine write_timing_rows(unit, timing)
    integer, intent(in) :: unit
    type(lcm_particle_exchange_timing_summary_t), intent(in) :: timing

    call write_one_timing(unit, 'x_outgoing_preparation',                 &
                          timing%x_outgoing_preparation)
    call write_one_timing(unit, 'y_outgoing_preparation',                 &
                          timing%y_outgoing_preparation)
    call write_one_timing(unit, 'x_count_exchange', timing%x_count_exchange)
    call write_one_timing(unit, 'y_count_exchange', timing%y_count_exchange)
    call write_one_timing(unit, 'x_payload_exchange',                     &
                          timing%x_payload_exchange)
    call write_one_timing(unit, 'y_payload_exchange',                     &
                          timing%y_payload_exchange)
    call write_one_timing(unit, 'x_received_insertion',                   &
                          timing%x_received_insertion)
    call write_one_timing(unit, 'y_received_insertion',                   &
                          timing%y_received_insertion)
    call write_one_timing(unit, 'ownership_update', timing%ownership_update)
    call write_one_timing(unit, 'advance_total', timing%advance_total)
    call write_one_timing(unit, 'prepare_particle_velocities',             &
                          timing%prepare_particle_velocities)
    call write_one_timing(unit, 'refine_particle_velocities',              &
                          timing%refine_particle_velocities)
    call write_one_timing(unit, 'update_particle_positions',               &
                          timing%update_particle_positions)
    call write_one_timing(unit, 'particle_boundaries',                      &
                          timing%particle_boundaries)
    call write_one_timing(unit, 'migration_total', timing%migration_total)
    call write_one_timing(unit, 'migration_buffer_clear',                  &
                          timing%migration_buffer_clear)
    call write_one_timing(unit, 'migration_plan_allocate',                 &
                          timing%migration_plan_allocate)
    call write_one_timing(unit, 'migration_plan_classification',           &
                          timing%migration_plan_classification)
    call write_one_timing(unit, 'migration_plan_apply',                    &
                          timing%migration_plan_apply)
    call write_one_timing(unit, 'migration_plan_deallocate',               &
                          timing%migration_plan_deallocate)
    call write_one_timing(unit, 'particle_exchange_total',                 &
                          timing%particle_exchange_total)
    call write_one_timing(unit, 'particle_sources', timing%particle_sources)
  end subroutine write_timing_rows

  subroutine write_one_timing(unit, name, timing)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: name
    type(lcm_rank_timing_summary_t), intent(in) :: timing

    call write_real_row(unit, 'timing', name, 'seconds', 'rank_minimum',  &
                        timing%minimum_seconds)
    call write_real_row(unit, 'timing', name, 'seconds', 'rank_maximum',  &
                        timing%maximum_seconds)
    call write_real_row(unit, 'timing', name, 'seconds', 'rank_mean',     &
                        timing%mean_seconds)
  end subroutine write_one_timing

  subroutine write_integer_rows(unit, category, name, metric_unit,        &
                                total, minimum, maximum)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: category
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: metric_unit
    integer(int64), intent(in) :: total
    integer(int64), intent(in) :: minimum
    integer(int64), intent(in) :: maximum

    call write_integer_row(unit, category, name, metric_unit, 'global_sum', &
                           total)
    call write_integer_row(unit, category, name, metric_unit, 'rank_minimum',&
                           minimum)
    call write_integer_row(unit, category, name, metric_unit, 'rank_maximum',&
                           maximum)
    call write_real_row(unit, category, name, metric_unit, 'rank_mean',    &
                        real(total, real64) / real(nprocs, real64))
  end subroutine write_integer_rows

  subroutine write_integer_row(unit, category, name, metric_unit,         &
                               aggregation, value)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: category
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: metric_unit
    character(len=*), intent(in) :: aggregation
    integer(int64), intent(in) :: value

    write(unit, '(a,4(",",a),",",i0,",",i0)') trim(cexpnr),           &
      trim(category), trim(name), trim(metric_unit), trim(aggregation),   &
      value, nprocs
  end subroutine write_integer_row

  subroutine write_real_row(unit, category, name, metric_unit,            &
                            aggregation, value)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: category
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: metric_unit
    character(len=*), intent(in) :: aggregation
    real(real64), intent(in) :: value

    write(unit, '(a,4(",",a),",",es24.16e3,",",i0)') trim(cexpnr),    &
      trim(category), trim(name), trim(metric_unit), trim(aggregation),   &
      value, nprocs
  end subroutine write_real_row

end module modlcm_metrics
#endif
