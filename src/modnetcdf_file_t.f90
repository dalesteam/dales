!> Type definitions for various NetCDF file types.
module modnetcdf_file_t

  use fortran_support, only: finish
  use modglobal,       only: imax, jmax, kmax, itot, jtot, rtimee, cexpnr
  use modmpi,          only: comm3d, myidx, myidy, cmyid, nprocx, nprocy, &
                             mpi_comm, commrow, commcol
  use modprecision,    only: field_r
  use modstat_nc

  implicit none

  private

  character(len=*), parameter :: modname = 'modnetcdf_file_t'

  public :: netcdf_file_t
  public :: time_series_file_t
  public :: profiles_file_t
  public :: cross_section_file_t
  public :: field_dump_file_t

  interface time_series_file_t
    procedure :: time_series_file_open
  end interface time_series_file_t

  interface profiles_file_t
    procedure :: profiles_file_open
  end interface profiles_file_t

  interface cross_section_file_t
    procedure :: cross_section_file_open
  end interface cross_section_file_t

  interface field_dump_file_t
    procedure :: field_dump_file_open
  end interface field_dump_file_t

  !> Base NetCDF file type.
  type, abstract :: netcdf_file_t
    character(len=80) :: filename      !< Name of the file.
    integer           :: ncid = 0      !< NetCDF file ID.
    integer           :: nvar = 0      !< Number of variables.
    integer           :: nrec = 0      !< Number of records.
    character(len=80) :: timeinfo(1,4) !< Metadata of time dimension.
    character(len=80), allocatable :: names(:,:) !< Variable metadata.
  contains
    procedure :: add_var => netcdf_file_add_var
    procedure :: get_var_id => netcdf_file_get_var_id
    procedure :: set_filename => netcdf_file_set_filename
    procedure(netcdf_file_init),  deferred :: init
    procedure(netcdf_file_write), deferred :: write
  end type netcdf_file_t

  abstract interface
    subroutine netcdf_file_init(this)
      import :: netcdf_file_t
      class(netcdf_file_t), intent(inout) :: this
    end subroutine netcdf_file_init
  end interface

  abstract interface
    subroutine netcdf_file_write(this)
      import :: netcdf_file_t
      class(netcdf_file_t), intent(inout) :: this
    end subroutine netcdf_file_write
  end interface

  !> File containing time series data.
  type, extends(netcdf_file_t) :: time_series_file_t
    private
    real(field_r), allocatable :: buffer(:) !< Memory for variable data.
  contains
    procedure :: init => time_series_file_init
    procedure :: write => time_series_file_write
    procedure :: get_pointer => time_series_file_get_pointer
  end type time_series_file_t

  !> File containing vertical profiles.
  type, extends(netcdf_file_t) :: profiles_file_t
    private
    integer :: nz = 0  !< Number of vertical levels.
    integer :: nzs = 0 !< Number of vertical levels in soil grid.
    real(field_r), allocatable :: buffer(:,:) !< Memory for variable data.
  contains
    procedure :: init => profiles_file_init
    procedure :: write => profiles_file_write
    procedure :: get_pointer => profiles_file_get_pointer
  end type profiles_file_t

  !> File containing cross sections.
  type, extends(netcdf_file_t) :: cross_section_file_t
    private
    integer :: nx = 0      !< Number of cells in the x-direction.
    integer :: ny = 0      !< Number of cells in the y-direction.
    integer :: nz = 0      !< Number of vertical layers.
    integer :: nzs = 0     !< Number of vertical layers in the soil grid.
    integer :: nslices = 0 !< Number of cross sections.
    integer :: x_start = 0 !< Writing offset in the x-direction.
    integer :: y_start = 0 !< Writing offset in the y-direction.
    integer :: nvals_x = 0 !< Number of values that will be written by calling process in the x-direction.
    integer :: nvals_y = 0 !< Number of values that will be written by calling process in the y-direction.
    real(field_r), allocatable :: buffer(:,:,:) !< Memory for variable data.
  contains
    procedure :: init => cross_section_file_init
    procedure :: write => cross_section_file_write
    procedure :: get_pointer => cross_section_file_get_pointer
  end type cross_section_file_t

  !> File containing 3D field dumps.
  type, extends(netcdf_file_t) :: field_dump_file_t
    private
    integer :: nx = 0      !< Number of cells in the x-direction.
    integer :: ny = 0      !< Number of cells in the y-direction.
    integer :: nz = 0      !< Number of vertical levels.
    integer :: nzs = 0     !< Number of vertical levels in the soil grid.
    integer :: ncoarse = 1 !< Coarse graining factor.
    integer :: klo = 1     !< Vertical lower bound.
    integer :: khi = 1     !< Vertical upper bound.
    integer :: x_start = 0 !< Writing offset in the x-direction.
    integer :: y_start = 0 !< Writing offset in the y-direction.
    integer :: nvals_x = 0 !< Number of values that will be written by calling process in the x-direction.
    integer :: nvals_y = 0 !< Number of values that will be written by calling process in the y-direction.
    real(field_r), allocatable :: buffer(:,:,:,:) !< Memory for variable data.
  contains
    procedure :: init => field_dump_file_init
    procedure :: write => field_dump_file_write
    procedure :: get_pointer => field_dump_file_get_pointer
  end type field_dump_file_t

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

    character(len=80), allocatable :: tmp(:,:)

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

    logical :: found = .false.

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
    character(len=80) :: the_filename

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

  !> Open a NetCDF file containing time series data.
  function time_series_file_open(filename) result(this)

    character(len=*), intent(in)  :: filename !< Name of the file.

    type(time_series_file_t) :: this !< New time series file object.

    call this%set_filename(filename)
    call open_nc(this%filename, this%ncid, this%nrec)

  end function time_series_file_open

  !> Define dimensions and allocate memory.
  subroutine time_series_file_init(this)

    class(time_series_file_t), intent(inout) :: this

    integer :: ivar

    call nctiminfo(this%timeinfo)

    if (this%nrec == 0) then
      call define_nc(this%ncid, 1, this%timeinfo)
      call writestat_dims_nc(this%ncid)
    end if

    call define_nc(this%ncid, this%nvar, this%names)

    allocate(this%buffer(this%nvar))

    do ivar = 1, this%nvar
      this%buffer(ivar) = 0.0_field_r
    end do

  end subroutine time_series_file_init

  !> Write data to disk.
  subroutine time_series_file_write(this)

    class(time_series_file_t), intent(inout) :: this

    integer :: ivar

    call writestat_nc(this%ncid, this%nvar, this%names, this%buffer, &
                      this%nrec, .true.)

    do ivar = 1, this%nvar
      this%buffer(ivar) = 0.0_field_r
    end do

  end subroutine time_series_file_write

  !> Setup a pointer to the buffer of a variable.
  subroutine time_series_file_get_pointer(this, name, ptr)

    class(time_series_file_t), target, intent(in) :: this
    character(len=*),                  intent(in) :: name !< Name of the variable.

    real(field_r), pointer, intent(out) :: ptr !< Pointer to the variable's buffer.

    character(len=*), parameter :: routine = &
      modname//'time_series_file_get_pointer'

    integer :: id

    id = this%get_var_id(name)

    if (id < 0) then
      call finish(routine, 'variable '//trim(name)//' not found')
    end if

    ptr => this%buffer(id)

  end subroutine time_series_file_get_pointer

  !> Open a NetCDF file containing vertical profiles.
  function profiles_file_open(filename, nz, nzs) result(this)

    character(len=*), intent(in) :: filename !< Name of the file.

    integer, intent(in), optional :: nz  !< Number of vertical levels.
    integer, intent(in), optional :: nzs !< Number of vertical levels in the soil grid.

    type(profiles_file_t) :: this !< New profiles file object.

    character(len=*), parameter :: routine = modname//'/profiles_file_open'

    if (.not. any([present(nz), present(nzs)], dim=1)) then
      call finish(routine, 'no vertical dimension length given')
    end if

    call this%set_filename(filename)
    
    if (present(nz)) this%nz = nz
    if (present(nzs)) this%nzs = nzs

    call open_nc(this%filename, this%ncid, this%nrec, n3=nz, ns=nzs)

  end function profiles_file_open

  !> Define dimensions and allocate memory.
  subroutine profiles_file_init(this)

    class(profiles_file_t), intent(inout) :: this

    integer :: k, ivar

    call nctiminfo(this%timeinfo)

    if (this%nrec == 0) then
      call define_nc(this%ncid, 1, this%timeinfo)
      call writestat_dims_nc(this%ncid)
    end if

    call define_nc(this%ncid, this%nvar, this%names)

    if (this%nz > 0) then
      allocate(this%buffer(this%nz,this%nvar))
    else if (this%nzs > 0) then
      allocate(this%buffer(this%nzs,this%nvar))
    end if

    do k = 1, size(this%buffer, dim=1)
      do ivar = 1, this%nvar
        this%buffer(k,ivar) = 0.0_field_r
      end do
    end do

  end subroutine profiles_file_init

  !> Write data to disk.
  subroutine profiles_file_write(this)

    class(profiles_file_t), intent(inout) :: this

    integer :: k, ivar

    call writestat_nc(this%ncid, 1, this%timeinfo, [rtimee], this%nrec, &
                      lraise=.true.)
    call writestat_nc(this%ncid, this%nvar, this%names, this%buffer, &
                      this%nrec, dim1=size(this%buffer, dim=1))

    do k = 1, size(this%buffer, dim=1) 
      do ivar = 1, this%nvar
        this%buffer(k,ivar) = 0.0_field_r
      end do
    end do

  end subroutine profiles_file_write

  !> Setup a pointer to the buffer of a variable.
  subroutine profiles_file_get_pointer(this, name, ptr)

    class(profiles_file_t), target, intent(in) :: this
    character(len=*),               intent(in) :: name !< Name of the variable.

    real(field_r), pointer, intent(out) :: ptr(:) !< Pointer to the variable's buffer.

    character(len=*), parameter :: routine = &
      modname//'/profiles_file_get_pointer'

    integer :: id

    id = this%get_var_id(name)

    if (id < 0) then
      call finish(routine, 'variable '//trim(name)//' not found')
    end if

    ptr => this%buffer(:,id)

  end subroutine profiles_file_get_pointer

  !> Open a NetCDF file containing cross sections.
  function cross_section_file_open(filename, nx, ny, nz, nzs) result(this)

    character(len=*), intent(in) :: filename !< Name of the file.

    integer, intent(in), optional :: nx  !< Number of cells in the x-direction.
    integer, intent(in), optional :: ny  !< Number of cells in the y-direction.
    integer, intent(in), optional :: nz  !< Number of vertical levels.
    integer, intent(in), optional :: nzs !< Number of vertical levels in the soil grid.

    type(cross_section_file_t) :: this !< New cross section file object.

    character(len=*), parameter :: routine = modname//'/cross_section_file_open'

    type(mpi_comm), pointer :: comm

    if (count([present(nx), present(ny), present(nz), present(nzs)], dim=1) > 2) then
      call finish(routine, 'cross section file can contain only two dimensions')
    end if
    
    if (present(nz)) this%nz = nz
    if (present(nzs)) this%nzs = nzs

    if (NC_HAVE_PARALLEL) then
      if (present(nx)) then
        this%nx = itot
        this%x_start = myidx * imax + 1
        this%nvals_x = imax
      end if
      if (present(ny)) then
        this%ny = jtot
        this%y_start = myidy * jmax + 1
        this%nvals_y = jmax
      end if

      call this%set_filename(filename)

      ! Determine which MPI communicator to use
      if (present(nx) .and. present(ny)) then
        comm => comm3d
      else if (present(nx)) then 
        comm => commrow
      else
        comm => commcol
      end if

      call open_nc(this%filename, this%ncid, this%nrec, n1=this%nx, &
                   n2=this%ny, n3=nz, ns=nzs, comm=comm)
    else
      if (present(nx)) then
        this%nx = imax
        this%x_start = 1
        this%nvals_x = imax
      end if
      if (present(ny)) then
        this%ny = imax
        this%y_start = 1
        this%nvals_y = jmax
      end if

      ! Every rank writes to its own file, distinguish with cmyid
      call this%set_filename(filename, suffix=cmyid)

      call open_nc(this%filename, this%ncid, this%nrec, n1=this%nvals_x, &
                   n2=this%nvals_y, n3=nz, ns=nzs)
    end if

  end function cross_section_file_open

  !> Define dimensions and allocate memory.
  subroutine cross_section_file_init(this)

    class(cross_section_file_t), intent(inout) :: this

    integer :: n1, n2, ivar

    call nctiminfo(this%timeinfo(1,:))

    if (this%nrec == 0) then
      call define_nc(this%ncid, 1, this%timeinfo, lcollective=.true.)
      call writestat_dims_nc(this%ncid, offset_x=this%x_start, &
                             offset_y=this%y_start)
    end if

    call define_nc(this%ncid, this%nvar, this%names, lcollective=.true.)

    if (this%nvals_x > 0 .and. this%nvals_y > 0) then
      n1 = this%nvals_x
      n2 = this%nvals_y
    else
      if (this%nvals_x > 0) then
        n1 = this%nvals_x
      else
        n1 = this%nvals_y
      end if
      if (this%nz > 0) then
        n2 = this%nz
      else
        n2 = this%nzs
      end if
    end if

    allocate(this%buffer(n1,n2,this%nvar))

    do ivar = 1, this%nvar
      do n2 = 1, size(this%buffer, dim=2)
        do n1 = 1, size(this%buffer, dim=1)
          this%buffer(n1,n2,ivar) = 0.0_field_r
        end do
      end do
    end do

  end subroutine cross_section_file_init

  !> Write data to disk.
  subroutine cross_section_file_write(this)

    class(cross_section_file_t), intent(inout) :: this

    integer :: n1, n2, ivar
    integer :: offsets(2)

    if (this%x_start > 0 .and. this%y_start > 0) then
      offsets = [this%x_start, this%y_start]
    else if (this%x_start > 0) then
      offsets = [this%x_start, 1]
    else if (this%y_start > 0) then
      offsets = [this%y_start, 1]
    else
      offsets = [1, 1]
    end if

    call writestat_nc(this%ncid, 1, this%timeinfo, [rtimee], this%nrec, &
                      lraise=.true.)
    call writestat_nc(this%ncid, this%nvar, this%names, this%buffer, &
                      this%nrec, dim1=size(this%buffer, dim=1), &
                      dim2=size(this%buffer, dim=2), offsets=offsets)

    do ivar = 1, this%nvar
      do n2 = 1, size(this%buffer, dim=2)
        do n1 = 1, size(this%buffer, dim=1)
          this%buffer(n1,n2,ivar) = 0.0_field_r
        end do
      end do
    end do

  end subroutine cross_section_file_write

  !> Setup a pointer to the buffer of a variable.
  !!
  !! @note Lower bounds of the pointer will be set to 2 for horizontal directions.
  subroutine cross_section_file_get_pointer(this, name, ptr)

    class(cross_section_file_t), target, intent(in) :: this
    character(len=*),                    intent(in) :: name

    real(field_r), pointer, intent(out) :: ptr(:,:)

    character(len=*), parameter :: routine = &
      modname//'/cross_section_file_get_pointer'

    integer :: id
    integer :: lbound_1 !< Lower bound of first pointer dimension
    integer :: lbound_2 !< Lower bound of second pointer dimension

    id = this%get_var_id(name)

    ! Set lower bounds to 2 so that we can easily access pointers and
    ! fields in conjunction
    lbound_1 = 2
    lbound_2 = 2

    ! For vertical slices, the second dimension (= vertical) starts from 1
    if (.not. (this%nx > 0 .and. this%ny > 0)) lbound_2 = 1

    ptr(lbound_1:,lbound_2:) => this%buffer(:,:,id)

  end subroutine cross_section_file_get_pointer

  !> Open a NetCDF file containing 3D field dumps.
  function field_dump_file_open(filename, nz, nzs, ncoarse, klo, khi) &
    result(this)

    character(len=*), intent(in) :: filename !< Name of the file.

    integer, intent(in), optional :: nz      !< Number of vertical levels.
    integer, intent(in), optional :: nzs     !< Number of vertical levels in the soil grid.
    integer, intent(in), optional :: ncoarse !< Coarse graining factor (e.g.: 2 means saving only half of the data).
    integer, intent(in), optional :: klo     !< Vertical lower bound.
    integer, intent(in), optional :: khi     !< Vertical upper bound.

    type(field_dump_file_t) :: this !< New field dump file object.

    character(len=*), parameter :: routine = modname//'/field_dump_file_open'

    if (.not. present(nz) .and. .not. present(nzs)) then
      call finish(routine, 'field dump file has to have at least one&
        & vertical dimension')
    end if

    if (present(ncoarse)) this%ncoarse = ncoarse

    if (present(klo)) then
      this%klo = klo
    else
      this%klo = 1
    end if
    
    if (present(khi)) then
      this%khi = khi
    else
      this%khi = kmax
    end if

    if (NC_HAVE_PARALLEL) then
      call this%set_filename(filename)
      this%nx = itot
      this%ny = jtot
      this%x_start = myidx * imax + 1
      this%y_start = myidy * jmax + 1
    else
      call this%set_filename(filename, suffix=cmyid)
      this%nx = imax
      this%ny = jmax
      this%x_start = 1
      this%y_start = 1
    end if

    this%nx = this%nx / ncoarse
    this%ny = this%ny / ncoarse
    this%nvals_x = imax / ncoarse
    this%nvals_y = jmax / ncoarse

    if (present(nz)) this%nz = khi - klo + 1 ! Passed value of nz not actually used...
    
    if (NC_HAVE_PARALLEL) then
      call open_nc(this%filename, this%ncid, this%nrec, n1=this%nvals_x, &
                   n2=this%nvals_y, n3=nz, ns=nzs, comm=comm3d)
    else
      call open_nc(this%filename, this%ncid, this%nrec, n1=this%nvals_x, &
                   n2=this%nvals_y, n3=nz, ns=nzs)
    end if

  end function field_dump_file_open

  !> Define dimensions and allocate memory.
  subroutine field_dump_file_init(this)

    class(field_dump_file_t), intent(inout) :: this

    integer :: n1, n2, n3, ivar

    call nctiminfo(this%timeinfo)

    if (this%nrec == 0) then
      call define_nc(this%ncid, 1, this%timeinfo, lcollective=.true.)
      call writestat_dims_nc(this%ncid, ncoarse=this%ncoarse, klow=this%klo, &
                             offset_x=this%x_start, offset_y=this%y_start)
    end if

    call define_nc(this%ncid, this%nvar, this%names, lcollective=.true.)

    if (this%nz > 0) then
      n3 = this%nz
    else
      n3 = this%nzs
    end if

    allocate(this%buffer(this%nvals_x,this%nvals_y,n3,this%nvar))

    do ivar = 1, this%nvar
      do n3 = 1, size(this%buffer, dim=3)
        do n2 = 1, size(this%buffer, dim=2)
          do n1 = 1, size(this%buffer, dim=1)
            this%buffer(n1,n2,n3,ivar) = 0.0_field_r
          end do
        end do
      end do
    end do

  end subroutine field_dump_file_init

  !> Write data to disk.
  subroutine field_dump_file_write(this)

    class(field_dump_file_t), intent(inout) :: this

    integer :: n1, n2, n3, ivar

    call writestat_nc(this%ncid, 1, this%timeinfo, [rtimee], this%nrec, &
                      lraise=.true.)
    call writestat_nc(this%ncid, this%nvar, this%names, this%buffer, &
                      this%nrec, dim1=this%nvals_x, dim2=this%nvals_y, &
                      dim3=size(this%buffer, dim=3), &
                      offsets=[this%x_start, this%y_start, 1])

    do ivar = 1, this%nvar
      do n3 = 1, size(this%buffer, dim=3)
        do n2 = 1, size(this%buffer, dim=2)
          do n1 = 1, size(this%buffer, dim=1)
            this%buffer(n1,n2,n3,ivar) = 0.0_field_r
          end do
        end do
      end do
    end do

  end subroutine field_dump_file_write

  !> Setup a pointer to the buffer of a variable.
  !!
  !! @note Lower bounds of the pointer will be set to 2 for horizontal directions.
  subroutine field_dump_file_get_pointer(this, name, ptr)

    class(field_dump_file_t), target, intent(in) :: this
    character(len=*),                 intent(in) :: name

    real(field_r), pointer, intent(out) :: ptr(:,:,:)

    character(len=*), parameter :: routine = &
      modname//'/field_dump_file_get_pointer'

    integer :: id

    id = this%get_var_id(name)

    if (id < 0) then
      call finish(routine, 'variable '//trim(name)//' not found')
    end if

    ptr(2:,2:,1:) => this%buffer(:,:,:,id)

  end subroutine field_dump_file_get_pointer

end module modnetcdf_file_t
