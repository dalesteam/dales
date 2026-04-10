!> Output file manager. Does the timekeeping for all output files.
module modstat_nc_files

  use modlogging,       only: finish
  use modglobal,        only: timee, rtimee, rk3step, dt_lim, tres, ladaptive, &
                              dtmax
  use modnetcdf_file_t, only: netcdf_file_t
  use modprecision,     only: field_r, longint

  implicit none

  private

  character(len=*), parameter :: modname = 'modstat_nc_files'

  public :: is_sampling_timestep
  public :: add_output_file
  public :: init_output_files
  public :: write_output_files
  public :: close_output_files
  public :: stats_limit_timestep

  integer, parameter :: MAX_FILES = 10 !< Max number of NetCDF files.

  type netcdf_file_list_entry_t
    class(netcdf_file_t), pointer :: file => null()
    integer(longint) :: dt_sample
    integer(longint) :: dt_write
  end type netcdf_file_list_entry_t

  integer :: nfiles = 0 !< Number of active output files

  type(netcdf_file_list_entry_t) :: file_list(MAX_FILES) !< List of output files

contains

  !> Check if sampling is due for a given file.
  function is_sampling_timestep(id) result(do_sample)

    integer, intent(in) :: id !< File identifier.

    logical :: do_sample !< True if sampling should be done for this file.

    if (rk3step == 3 .and. timee > 0.001) then
      do_sample = mod(timee, file_list(id)%dt_sample) == 0.0_field_r
    else
      do_sample = .false.
    end if

  end function is_sampling_timestep

  !> Check if writing is due for a given file.
  function is_writing_timestep(id) result(do_write)

    integer, intent(in) :: id !< File identifier.

    logical :: do_write !< True if writing should be done for this file.

    if (rk3step == 3 .and. timee > 0.001) then
      do_write = mod(timee, file_list(id)%dt_write) == 0.0_field_r
    else
      do_write = .false.
    end if

  end function is_writing_timestep

  !> Add a new file to the list of files.
  subroutine add_output_file(file, dt_sample, id, dt_write)

    class(netcdf_file_t), target, intent(in) :: file
    real,                         intent(in) :: dt_sample

    integer, intent(out) :: id

    real, intent(in), optional :: dt_write

    character(len=*), parameter :: routine = modname//'/add_output_file'

    if (dt_sample < 0.0_field_r) then
      call finish(routine, 'dt_sample cannot be negative (file: '&
                  //trim(file%filename)//')')
    end if

    if (.not. ladaptive .and. mod(dt_sample, dtmax) > 1.0E-4_field_r) then
      call finish(routine, 'adaptive time stepping is disabled, so dt_sample&
        & should be an integer multiple of dtmax (file: '&
        //trim(file%filename)//')')
    end if

    nfiles = nfiles + 1
    id = nfiles

    file_list(id)%file => file 
    file_list(id)%dt_sample = int(dt_sample / tres, kind=longint)

    if (present(dt_write)) then
      if (mod(dt_write, dt_sample) > 1.0E-4_field_r) then
        call finish(routine, 'dt_write should be an integer multiple of&
          & dt_sample (file: '//trim(file%filename)//')')
      end if
      file_list(id)%dt_write = int(dt_write / tres, kind=longint)
    else
      file_list(id)%dt_write = file_list(id)%dt_sample
    end if

  end subroutine add_output_file

  !> Calls the initialization subroutine for all files in the file list
  subroutine init_output_files()

    integer :: ifile

    do ifile = 1, nfiles
      call file_list(ifile)%file%open
    end do

  end subroutine init_output_files

  !> Loop over the files and write those that are due.
  subroutine write_output_files()

    integer :: ifile

    do ifile = 1, nfiles
      if (is_writing_timestep(ifile)) then
        call file_list(ifile)%file%write
      end if
    end do

  end subroutine write_output_files

  !> Limit the time step if needed for sampling or writing
  subroutine stats_limit_timestep()

    integer          :: ifile

    integer(longint) :: dts       !< Delta t for sampling
    integer(longint) :: dtw       !< Delta t for writing
    integer(longint) :: time_left !< Time left before sampling or writing needs to be done

    if (rk3step == 3) then
      do ifile = 1, nfiles
        dts = file_list(ifile)%dt_sample
        dtw = file_list(ifile)%dt_write
        time_left = min(dts - mod(timee, dts), &
                        dtw - mod(timee, dtw))
        dt_lim = min(dt_lim, time_left)
      end do
    end if

  end subroutine stats_limit_timestep

  !> Closes all files.
  subroutine close_output_files()

    integer :: ifile

    do ifile = 1, nfiles
      call file_list(ifile)%file%close()
    end do

  end subroutine close_output_files

end module modstat_nc_files