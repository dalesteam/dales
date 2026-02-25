module modstat_2d

  use fortran_support, only: finish, message, nnml_output
  use modglobal,       only: ifnamopt, checknamelisterror, imax, jmax, i1, j1, &
                             tres, rtimee, dtav_glob, timeav_glob, btime, &
                             timee, dt_lim, rk3step, fname_options
  use modmpi,          only: cmyidx, cmyidy, myid, d_mpi_bcast, comm3d
  use modprecision,    only: field_r, longint
  use modstat_nc,      only: nctiminfo, open_nc, define_nc, writestat_dims_nc, &
                             exitstat_nc, writestat_nc, ncinfo

  implicit none

  private

  character(len=*), parameter :: modname = 'modstat_2d'

  public :: stat_2d_read_namelist
  public :: add_slice
  public :: init_stat_2d
  public :: write_2d
  public :: get_slice
  public :: sample_2d
  public :: do_stats

  integer          :: nvar = 0 !< Number of statistical variables.
  real             :: dtav     !< Sampling interval [s].
  integer(longint) :: idtav    !< Integer sampling interval.
  integer(longint) :: tnext    !< Next sampling time step.

  character(len=80) :: fname !< Filename.
  integer           :: ncid  !< NetCDF file ID.
  integer           :: nrec  !< Number of records in file (= n time steps).

  character(len=80), allocatable :: ncname(:,:)
  character(len=80)              :: tncname(1,4)

  logical :: do_stats
  logical :: write_stats
  integer :: nsamples

  real(field_r), allocatable, target :: output_data(:,:,:) !< Array containing all 2D slices.

contains

  subroutine stat_2d_read_namelist(nml_filename)

    character(len=*), intent(in) :: nml_filename

    integer :: ierr

    namelist /output_2d/ dtav

    dtav = dtav_glob

    if (myid == 0) then
      open(ifnamopt, file=fname_options, status='old', iostat=ierr)
      read(ifnamopt, output_2d, iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'output_2d')
      write(nnml_output, output_2d)
      close(ifnamopt)
    end if

    call d_mpi_bcast(dtav, 1, 0, comm3d, ierr)

  end subroutine stat_2d_read_namelist

  subroutine init_stat_2d()

    idtav = int(dtav / tres, kind=longint)

    tnext = idtav + btime

    do_stats = .false.

    fname = 'dales-output-2d.nc'

    if (nvar > 0) then

      allocate(output_data(2:i1,2:j1,1:nvar))

      output_data(:,:,:) = 0
      
      !$acc enter data copyin(output_data)

      call nctiminfo(tncname(1,:))
      call open_nc(fname, ncid, nrec, n1=imax, n2=jmax, lparallel=.true.)

      if (nrec == 0) then
        call define_nc(ncid, 1, tncname)
        call writestat_dims_nc(ncid, lparallel=.true.)
      end if

      call define_nc(ncid, nvar, ncname)

    end if

  end subroutine init_stat_2d

  function find_index(name) result(index)

    character(len=*), intent(in) :: name

    integer :: index

    do index = 1, nvar
      if (trim(name) == trim(ncname(index,1))) return
    end do

    index = 0

  end function find_index

  !> Add a 2D slice to the output and setup a pointer to the buffer.
  subroutine add_slice(name, long_name, unit, dim)

    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: long_name
    character(len=*), intent(in) :: unit
    character(len=*), intent(in) :: dim

    character(len=*), parameter :: routine = modname//'/add_slice'

    character(len=80), allocatable :: tmp_ncname(:,:)

    ! Allocate array for metadata
    if (.not. allocated(ncname)) then
      allocate(ncname(1,4))
    else
      ! Check if given name already exists. For the long name, we don't care.
      if (find_index(name) /= 0) then
        call finish(routine, 'variable '//trim(name)//' already exists')
      else
        ! If already allocated, grow in size by 1
        allocate(tmp_ncname(size(ncname, dim=1) + 1, 4))
        tmp_ncname(1:nvar,:) = ncname(1:nvar,:)
        call move_alloc(tmp_ncname, ncname)
      end if
    end if

    nvar = nvar + 1

    call ncinfo(ncname(nvar,:), name, long_name, unit, dim)

  end subroutine add_slice

  !> Get a pointer to the output buffer by variable name.
  function get_slice(name) result(ptr)

    character(len=*), intent(in) :: name

    character(len=*), parameter :: routine = modname//'/get_slice'

    real(field_r), pointer :: ptr(:,:)

    integer :: i
    logical :: found = .false.

    do i = 1, nvar
      if (trim(name) == trim(ncname(i,1))) then
        found = .true.
        exit
      end if
    end do

    if (.not. found) then
      call finish(routine, 'variable '//trim(name)//' not found')
    end if

    ptr(2:,2:) => output_data(2:,2:,i)

  end function get_slice

  subroutine sample_2d()

    ! Reset switch
    do_stats = .false.

    if (nvar < 0) return
    if (rk3step /= 3) return

    if (timee < tnext) then
      dt_lim = minval([dt_lim, tnext - timee])
    else
      do_stats = .true.
      tnext = tnext + idtav
    end if

  end subroutine sample_2d

  subroutine write_2d()

    integer :: i, j, n

    if (do_stats) then

      !$acc update host(output_data)

      call writestat_nc(ncid, 1, tncname, [rtimee], nrec, .true.)
      call writestat_nc(ncid, nvar, ncname, output_data(:,:,:), nrec, imax, jmax, lparallel=.true.)

      !$acc parallel loop collapse(3) default(present)
      do n = 1, nvar
        do j = 2, j1
          do i = 2, i1
            output_data(i,j,n) = 0
          end do
        end do
      end do

    end if

  end subroutine write_2d


end module modstat_2d