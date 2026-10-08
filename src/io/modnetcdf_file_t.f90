module modnetcdf_file_t

  use fortran_support, only: finish
  use modglobal,       only: cexpnr
  use modstat_nc,      only: exitstat_nc

  implicit none

  private

  character(len=*), parameter :: modname = 'modnetcdf_file_t'

  type, abstract :: netcdf_file_t
    character(len=64)  :: filename
    integer            :: ncid = 0
    integer            :: nvar = 0
    integer            :: nrec = 0
    character(len=128) :: timeinfo(1,4)
    character(len=128), allocatable :: names(:,:)
  contains
    procedure :: add_var => netcdf_file_add_var
    procedure :: get_var_id => netcdf_file_get_var_id
    procedure :: set_filename => netcdf_file_set_filename
    procedure :: close => netcdf_file_close
    procedure(netcdf_file_open),  deferred :: open
    procedure(netcdf_file_write), deferred :: write
  end type netcdf_file_t

  abstract interface
    subroutine netcdf_file_open(this)
      import :: netcdf_file_t
      class(netcdf_file_t), intent(inout) :: this
    end subroutine netcdf_file_open
  end interface

  abstract interface
    subroutine netcdf_file_write(this)
      import :: netcdf_file_t
      class(netcdf_file_t), intent(inout) :: this
    end subroutine netcdf_file_write
  end interface

contains

  !> Add a variable to a NetCDF file.
  subroutine netcdf_file_add_var(this, name, long_name, unit, dim)

    class(netcdf_file_t), intent(inout) :: this

    character(len=*), intent(in) :: name       !< Name of variable.
    character(len=*), intent(in) :: long_name  !< Long name of variable.
    character(len=*), intent(in) :: unit       !< Unit of variable.
    character(len=*), intent(in) :: dim        !< Dimensions of variable.

    character(len=*), parameter :: routine = modname//'/netcdf_file_add_var'

    integer :: id

    character(len=128), allocatable :: tmp(:,:)

    ! Check if this variable already exists
    id = this%get_var_id(name)

    if (id > 0) then
      call finish(routine, 'variable '//trim(name)//' already exists in file ' &
                            //trim(this%filename))
    else
      ! Grow the list of names
      this%nvar = this%nvar + 1
      allocate(tmp(this%nvar,4)) 
      if (allocated(this%names)) tmp(1:this%nvar-1,:) = this%names(:,:)
      call move_alloc(tmp, this%names)

      ! Add info of the new variable
      this%names(this%nvar,1) = trim(name)
      this%names(this%nvar,2) = trim(long_name)
      this%names(this%nvar,3) = trim(unit)
      this%names(this%nvar,4) = trim(dim)
    end if

  end subroutine netcdf_file_add_var

  !> Get the internal index of variable by name.
  function netcdf_file_get_var_id(this, name) result(id)

    class(netcdf_file_t), intent(in) :: this
    character(len=*),     intent(in) :: name !< Name of the variable.

    integer :: id !< Index of variable in variable list of this file

    !we explicitly don't set. found=.false. at the declaration as this will persist the result of a search across calls...
    logical :: found
    found = .false.

    do id = 1, this%nvar
      if (trim(name) == trim(this%names(id,1))) then
        found = .true.
        exit
      end if
    end do

    if (.not. found) id = -1

  end function netcdf_file_get_var_id

  !> Sets the filename, takes care of the .nc suffix and any other suffix.
  subroutine netcdf_file_set_filename(this, filename, suffix)

    class(netcdf_file_t), intent(inout) :: this

    character(len=*), intent(in) :: filename

    character(len=*), intent(in), optional :: suffix

    integer           :: strlen
    character(len=64) :: the_filename

    the_filename = adjustr(filename)
    strlen = len(trim(the_filename))

    ! Find out if the given filename already ends in '.nc'
    ! If so, strip it
    if (the_filename(strlen-3+1:) == '.nc') then
      the_filename(strlen-3+1:) = ' '
      the_filename = trim(the_filename)
    end if

    if (present(suffix)) the_filename = trim(the_filename)//'.'//trim(suffix)

    ! Add the experiment ID and the .nc suffix
    the_filename = trim(the_filename)//'.'//cexpnr//'.nc'

    this%filename = the_filename

  end subroutine netcdf_file_set_filename

  !> Close the NetCDF file.
  subroutine netcdf_file_close(this)

    class(netcdf_file_t), intent(inout) :: this

    if (this%ncid /= 0) call exitstat_nc(this%ncid)

  end subroutine netcdf_file_close

end module modnetcdf_file_t