!> Output file manager.
module modstat_nc_files

  use modglobal,        only: timee, tres, dt_lim, rk3step
  use modmpi,           only: mpi_comm, cmyid, myid
  use modnetcdf_file_t, only: netcdf_file_t
  use modprecision,     only: longint, field_r
  use modstat_nc,       only: NC_HAVE_PARALLEL

  implicit none

  private

  public :: make_netcdf_file
  public :: is_sampling_timestep
  public :: add_variable
  public :: init_nc_files
  public :: write_nc_files
  public :: stats_limit_timestep
  public :: get_pointer

  integer, parameter :: MAX_FILES = 10 !< Max number of NetCDF files.

  integer :: nfiles = 0

  type(netcdf_file_t), target :: file_list(MAX_FILES)

  ! Time keeping
  integer(longint) :: sampling_dts(MAX_FILES)
  integer(longint) :: writing_dts(MAX_FILES)

  interface get_pointer
    procedure :: get_pointer_0d
    procedure :: get_pointer_1d
    procedure :: get_pointer_2d
    procedure :: get_pointer_3d
  end interface

contains

  !> Define a new NetCDF file.
  function make_netcdf_file(filename, dimension_lengths, dt_sample, dt_write, &
                            nprocs, ranks, comm) result(id)

    character(len=*),  intent(in) :: filename             !< Base file name
    integer,           intent(in) :: dimension_lengths(5) !< Length of each dimension
    real,              intent(in) :: dt_sample            !< Interval of sampling

    integer,        intent(in), optional :: dt_write  !< Interval of writing
    integer,        intent(in), optional :: nprocs(5) !< Number of processes in each dimension
    integer,        intent(in), optional :: ranks(5)  !< Rank of the calling process in each dimension
    type(mpi_comm), intent(in), optional :: comm      !< MPI communicator containing all ranks that participate in this file

    integer           :: idim

    integer           :: id                !< File identifier (NOT the ncid!)
    integer           :: offsets(5)        !< Offset in each dimension
    integer           :: nvals(5)          !< Number of values per dim that each rank will write
    integer           :: my_dim_lengths(5) !< Actual dimension lengths
    character(len=80) :: suffix = ''       !< File suffix: <filename>.<suffix>.nc

    nfiles = nfiles + 1
    id = nfiles

    if (present(nprocs)) then
      ! Multiple processes write to this file
      if (any(nprocs > 1)) then
        if (present(comm) .and. NC_HAVE_PARALLEL) then
          ! All ranks write to the same file
          do idim = 1, size(dimension_lengths)
            if (dimension_lengths(idim) > 0) then
              nvals(idim) = dimension_lengths(idim) / nprocs(idim)
              offsets(idim) = nvals(idim) * ranks(idim) + 1
            else
              nvals(idim) = 0
              offsets(idim) = 0
            end if
            my_dim_lengths(idim) = dimension_lengths(idim)
          end do
        else
          ! One file per rank
          do idim = 1, size(dimension_lengths)
            if (dimension_lengths(idim) > 0) then
              nvals(idim) = dimension_lengths(idim) / nprocs(idim)
              offsets(idim) = 1
            else
              nvals(idim) = 0
              offsets(idim) = 0
            end if
            my_dim_lengths(:) = nvals(:)
            ! Set suffix of each file to MPI id
            suffix = cmyid
          end do
        end if
      end if
    end if

    file_list(id) = netcdf_file_t(filename//trim(suffix), my_dim_lengths, &
                                  nvals, offsets, comm)

    sampling_dts(id) = int(dt_sample / tres, kind=longint)

    if (present(dt_write)) then
      writing_dts(id) = int(dt_write / tres, kind=longint)
    else
      writing_dts(id) = sampling_dts(id)
    end if

  end function make_netcdf_file

  !> Check if sampling is due for a given file.
  function is_sampling_timestep(id) result(do_sample)

    integer, intent(in) :: id !< File identifier

    logical :: do_sample !< True if sampling should be done for this file

    if (rk3step == 3) then
      do_sample = mod(timee, sampling_dts(id)) == 0
    else
      do_sample = .false.
    end if

  end function is_sampling_timestep

  function is_writing_timestep(id) result(do_write)

    integer, intent(in) :: id !< File identifier

    logical :: do_write !< True if writing should be done for this file

    if (rk3step == 3) then
      do_write = mod(timee, writing_dts(id)) == 0
    else
      do_write = .false.
    end if

  end function is_writing_timestep

  !> Add a variable to a file.
  subroutine add_variable(id, name, long_name, unit, dim)

    integer,          intent(in) :: id        !< File identifier
    character(len=*), intent(in) :: name      !< Name of the variable
    character(len=*), intent(in) :: long_name !< Long name of the variable
    character(len=*), intent(in) :: unit      !< Unit of the variable
    character(len=*), intent(in) :: dim       !< Dimension of the variable

    call file_list(id)%add_var(name, long_name, unit, dim)

  end subroutine add_variable

  !> Calls the initialization subroutine for all files in the file list
  subroutine init_nc_files()

    integer :: ifile

    do ifile = 1, nfiles
      call file_list(ifile)%init
    end do

  end subroutine init_nc_files

  !> Loop over the files and write those that are due.
  subroutine write_nc_files()

    integer :: ifile

    do ifile = 1, nfiles
      if (is_writing_timestep(ifile)) then
        call file_list(ifile)%write
      end if
    end do

  end subroutine write_nc_files

  !> Limit the time step if needed for sampling or writing
  subroutine stats_limit_timestep()

    integer          :: ifile

    integer(longint) :: dts       !< Delta t for sampling
    integer(longint) :: dtw       !< Delta t for writing
    integer(longint) :: time_left !< Time left before sampling or writing needs to be done

    if (rk3step == 3) then
      do ifile = 1, nfiles
        dts = sampling_dts(ifile)
        dtw = writing_dts(ifile)
        time_left = min(dts - mod(timee, dts), &
                        dtw - mod(timee, dtw))
        dt_lim = min(dt_lim, time_left)
      end do
    end if

  end subroutine stats_limit_timestep

  !> Setup pointer to the buffer of a scalar variable in a given file
  subroutine get_pointer_0d(id, name, ptr)

    integer,          intent(in) :: id   !< File identifier
    character(len=*), intent(in) :: name !< Name of the variable

    real(field_r), pointer, intent(out) :: ptr !< Pointer to the variable buffer

    integer :: var_id

    var_id = file_list(id)%get_var_id(name)

    if (var_id /= -1) then
      ptr => file_list(id)%data_0d(var_id)
    else
      ptr => null()
    end if

  end subroutine get_pointer_0d

  !> Setup pointer to the buffer of a 1D variable in a given file
  subroutine get_pointer_1d(id, name, ptr, lbound)

    integer,          intent(in) :: id   !< File identifier
    character(len=*), intent(in) :: name !< Name of the variable

    real(field_r), pointer, intent(out) :: ptr(:) !< Pointer to the variable buffer

    integer, intent(in), optional :: lbound !< Lower bound of the pointer

    integer :: var_id
    integer :: start = 1

    if (present(lbound)) start = lbound

    var_id = file_list(id)%get_var_id(name)

    if (var_id /= -1) then
      ptr(start:) => file_list(id)%data_1d(:,var_id)
    else
      ptr => null()
    end if

  end subroutine get_pointer_1d

  !> Setup pointer to the buffer of a 2D variable in a given file
  subroutine get_pointer_2d(id, name, ptr, lbound)

    integer,          intent(in) :: id   !< File identifier
    character(len=*), intent(in) :: name !< Name of the variable

    real(field_r), pointer, intent(out) :: ptr(:,:) !< Pointer to the variable buffer

    integer, intent(in), optional :: lbound(2) !< Lower bounds of the pointer

    integer :: var_id
    integer :: start(2)

    if (present(lbound)) then
      start(:) = lbound(:)
    else
      start(:) = 1
    end if

    var_id = file_list(id)%get_var_id(name)

    if (var_id /= -1) then
      ptr(start(1):,start(2):) => file_list(id)%data_2d(:,:,var_id)
    else
      ptr => null()
    end if

  end subroutine get_pointer_2d

  !> Setup pointer to the buffer of a 3D variable in a given file
  subroutine get_pointer_3d(id, name, ptr, lbound)

    integer,          intent(in) :: id   !< File identifier
    character(len=*), intent(in) :: name !< Name of the variable

    real(field_r), pointer, intent(out) :: ptr(:,:,:) !< Pointer to the variable buffer

    integer, intent(in), optional :: lbound(3) !< Lower bounds of the pointer

    integer :: var_id
    integer :: start(3)

    if (present(lbound)) then
      start(:) = lbound(:)
    else
      start(:) = 1
    end if

    var_id = file_list(id)%get_var_id(name)

    if (var_id /= -1) then
      ptr(start(1):,start(2):,start(3):) => file_list(id)%data_3d(:,:,:,var_id)
    else
      ptr => null()
    end if

  end subroutine get_pointer_3d

end module modstat_nc_files
