!> NetCDF file type.
module modnetcdf_file_t

  use fortran_support, only: finish, message
  use modglobal,       only: rtimee
  use modmpi,          only: mpi_comm
  use modprecision,    only: field_r
  use modstat_nc

  implicit none

  private

  character(len=*), parameter :: modname = 'modnetcdf_file_t'

  public :: netcdf_file_t

  interface netcdf_file_t
    procedure :: netcdf_file_open
  end interface netcdf_file_t

  type :: netcdf_file_t
    ! Metadata
    character(len=80) :: filename = '' !< File name
    character(len=80) :: suffix = ''   !< File suffix: <filename>.<suffix>.nc
    integer           :: ncid = -1     !< NCID
    integer           :: nrec = -1     !< Number of records (= time steps)
    integer           :: nvar = 0      !< Number of variables
    character(len=80) :: timeinfo(1,4) !< Metadata for time
    integer           :: ndims         !< Number of used dimensions
    ! Information about the dimensions (local to each MPI rank)
    integer, allocatable :: offsets(:) !< Starting index of each dimension
    integer, allocatable :: nvals(:)   !< Number of values in each dimension
    ! Variable data
    character(len=80), allocatable :: names(:,:) !< dim 1 = name, dim 2 = long name, dim 3 = unit, dim 4 = dimension
    real(field_r),     allocatable :: data_0d(:)
    real(field_r),     allocatable :: data_1d(:,:)
    real(field_r),     allocatable :: data_2d(:,:,:)
    real(field_r),     allocatable :: data_3d(:,:,:,:)
  contains
    procedure :: add_var => netcdf_file_add_var
    procedure :: close => netcdf_file_close
    procedure :: write => netcdf_file_write
    procedure :: init => netcdf_file_init
    procedure :: get_var_id => netcdf_file_get_var_id
  end type netcdf_file_t

contains

  !> Open a new NetCDF file.
  function netcdf_file_open(filename, dimension_lengths, nvals, offsets, comm, &
                            suffix) result(this)

    character(len=*), intent(in) :: filename             !< File name
    integer,          intent(in) :: dimension_lengths(5) !< Lengths of the dimensions [n1, n2, n3, ns, nq].

    integer,          intent(in), optional :: nvals(5)   !< How much values will be written per dimension. 
    integer,          intent(in), optional :: offsets(5) !< Starting indices in each dimension.
    type(mpi_comm),   intent(in), optional :: comm       !< MPI communicator containing all ranks that write to this file.
    character(len=*), intent(in), optional :: suffix     !< File suffix

    type(netcdf_file_t) :: this !< New NetCDF file object

    integer        :: idim

    integer        :: ndims !< Number of dimensions with length greater than 0

    this%filename = filename

    ! Check how many dimensions we have to define
    ndims = count(dimension_lengths > 0, dim=1)
    this%ndims = ndims
    allocate(this%offsets(ndims), this%nvals(ndims))

    if (present(suffix)) then
      this%suffix = '.'//trim(suffix)
    end if

    ! Set offsets and number of values for this process
    do idim = 1, 5
      if (dimension_lengths(idim) > 0) then 
        if (present(offsets)) then
          this%offsets(idim) = offsets(idim)
        else
          this%offsets(idim) = 1
        end if

        if(present(nvals)) then
          this%nvals(idim) = nvals(idim)
        else
          this%nvals(idim) = 1
        end if
      end if
    end do

    ! Finally, open the file
    call open_nc(trim(this%filename)//trim(this%suffix)//'.nc', this%ncid, &
                 this%nrec, dimension_lengths(1), &
                 dimension_lengths(2), dimension_lengths(3), &
                 dimension_lengths(4), dimension_lengths(5), comm)

  end function netcdf_file_open

  !> Add a variable to the file.
  subroutine netcdf_file_add_var(this, name, long_name, unit, dim)

    class(netcdf_file_t), intent(inout) :: this

    character(len=*), intent(in) :: name      !< Name of variable
    character(len=*), intent(in) :: long_name !< Long name of variable
    character(len=*), intent(in) :: unit      !< Unit of variable
    character(len=*), intent(in) :: dim       !< Dimensions of variable

    character(len=*), parameter :: routine = modname//'/netcdf_file_add_var'

    integer :: i

    logical :: found = .false.

    character(len=80), allocatable :: tmp(:,:)

    ! Check if the new variable already exists in this file
    do i = 1, this%nvar
      if (trim(name) == trim(this%names(i,1))) then
        found = .true.
      end if
    end do

    if (found) then
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

      call message(routine, 'variable '//trim(name)//' added to file ' &
                            //trim(this%filename))
    end if

  end subroutine netcdf_file_add_var

  !> Ends the definition stage and allocates memory for the variables.
  subroutine netcdf_file_init(this)

    class(netcdf_file_t), intent(inout) :: this

    integer :: buffer_len

    select case (this%ndims)
      case (0)
        allocate(this%data_0d(this%nvar))
        this%data_0d(:) = 0.0_field_r
      case (1)
        allocate(this%data_1d(this%nvals(1), this%nvar))
        this%data_1d(:,:) = 0.0_field_r
      case (2)
        allocate(this%data_2d(this%nvals(1), this%nvals(2), this%nvar))
        this%data_2d(:,:,:) = 0.0_field_r
      case (3)
        allocate(this%data_3d(this%nvals(1), this%nvals(2), this%nvals(3), &
                              this%nvar))
        this%data_3d(:,:,:,:) = 0.0_field_r
      case default
    end select

    call nctiminfo(this%timeinfo)

    if (this%nrec == 0) then
      call define_nc(this%ncid, 1, this%timeinfo)
      call writestat_dims_nc(this%ncid)
    end if

    call define_nc(this%ncid, this%nvar, this%names)

  end subroutine netcdf_file_init

  !> Write variable data to disk and reset the buffer.
  subroutine netcdf_file_write(this)

    class(netcdf_file_t), intent(inout) :: this
    
    integer :: i, j, k, ivar

    call writestat_nc(this%ncid, 1, this%timeinfo, [rtimee], this%nrec, .true.)
    
    select case (this%ndims)
      case (0)
        call writestat_nc(this%ncid, this%nvar, this%names, this%data_0d, &
                          this%nrec, .false.)
        
        do ivar = 1, this%nvar
          this%data_0d(ivar) = 0
        end do

      case (1)
        call writestat_nc(this%ncid, this%nvar, this%names(:,:), this%data_1d, &
                          this%nrec, this%nvals(1))

        do ivar = 1, this%nvar
          do k = 1, size(this%data_1d, dim=1)
            this%data_1d(k,ivar) = 0
          end do
        end do

      case (2)
        call writestat_nc(this%ncid, this%nvar, this%names, this%data_2d, &
                          this%nrec, this%nvals(1), this%nvals(2), &
                          offsets=this%offsets)

        do ivar = 1, this%nvar
          do j = 1, size(this%data_2d, dim=2)
            do i = 1, size(this%data_2d, dim=1)
              this%data_2d(i,j,ivar) = 0
            end do
          end do
        end do

      case (3)
        call writestat_nc(this%ncid, this%nvar, this%names, this%data_3d, &
                          this%nrec, this%nvals(1), this%nvals(2), &
                          this%nvals(3), offsets=this%offsets) 

        do ivar = 1, this%nvar
          do k = 1, size(this%data_3d, dim=3)
            do j = 1, size(this%data_3d, dim=2)
              do i = 1, size(this%data_3d, dim=1)
                this%data_3d(i,j,k,ivar) = 0
              end do
            end do
          end do
        end do

      case default
    end select

  end subroutine netcdf_file_write

  !> Close NetCDF file.
  subroutine netcdf_file_close(this)

    class(netcdf_file_t), intent(inout) :: this

    call exitstat_nc(this%ncid)
  
  end subroutine netcdf_file_close

  !> Get the internal index of variable by name.
  function netcdf_file_get_var_id(this, name) result(id)

    class(netcdf_file_t), intent(in) :: this
    character(len=*),     intent(in) :: name !< Name of the variable

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

end module modnetcdf_file_t