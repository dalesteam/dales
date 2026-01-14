module modstat_profiles

  use modglobal,      only: kmax, fname_options, ifnamopt, checknamelisterror, &
                            i1, j1, k1, kmax, ijtot, btime, ih, dt_lim, timee, &
                            tres, rtimee, imax, cexpnr, dtav_glob, timeav_glob, &
                            rk3step
  use modmpi,         only: d_mpi_bcast, comm3d, mpierr, print_info_stderr, &
                            cmyidx, cmyidy
  use modstat_nc
  use modprecision,   only: field_r, longint
  use modslabaverage, only: slabavg
  use modlogging, only: finish
  use fortran_support, only: nnml_output

  implicit none

  private

  character(len=*), parameter :: modname = 'modstat_profiles'

  public :: is_sampling_timestep
  public :: is_writing_timestep
  public :: do_procblock

  public :: add_profile
  public :: init_profiles
  public :: sample_profiles
  public :: sample_field
  public :: write_profiles
  public :: exit_profiles

  interface sample_field
    module procedure sample_field
    module procedure sample_field_masked
  end interface sample_field

  ! Namelist options
  logical :: lstat
  logical :: lprocblock = .false.
  real    :: dtav, timeav

  ! NetCDF file variables
  character(len=80) :: fname
  integer           :: ncid, nrec

  ! NetCDF data
  character(len=80), allocatable :: ncname(:,:)
  character(len=80)              :: tncname(1,4)
  real,              allocatable :: profiles(:,:)

  integer          :: nvar = 0
  integer(longint) :: idtav, itimeav, tnext, tnextwrite
  integer          :: nsamples
  logical          :: do_stats
  logical          :: write_stats
  logical          :: my_task_writes

  real(field_r), allocatable :: slab_average(:)

contains

  logical function is_sampling_timestep()

    is_sampling_timestep = do_stats

  end function is_sampling_timestep

  logical function is_writing_timestep()

    is_writing_timestep = write_stats

  end function is_writing_timestep

  logical function do_procblock()

    do_procblock = lprocblock

  end function do_procblock

  function find_index(name) result(index)

    character(len=*), intent(in) :: name

    integer :: index

    do index = 1, nvar
      if (trim(name) == trim(ncname(index,1))) return
    end do

    index = 0

  end function find_index

  subroutine add_profile(name, long_name, unit, dim)

    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: long_name
    character(len=*), intent(in) :: unit
    character(len=*), intent(in) :: dim

    character(len=*), parameter :: routine = modname//'/add_profile'

    character(len=80), allocatable :: tmp_ncname(:,:)
    integer :: idx

    ! Allocate array for metadata
    if (.not. allocated(ncname)) then
      allocate(ncname(1,4))
    else
      ! Check if given name already exists. For the long name, we don't care.
      if (find_index(name) /= 0) then
        call finish(routine, 'profile ', trim(name), ' already exists')
      else
        ! If already allocated, grow in size by 1
        allocate(tmp_ncname(size(ncname, dim=1) + 1, 4))
        tmp_ncname(1:nvar,:) = ncname(1:nvar,:)
        call move_alloc(tmp_ncname, ncname)
      end if
    end if

    nvar = nvar + 1

    call ncinfo(ncname(nvar,:), name, long_name, unit, dim)

  end subroutine add_profile

  ! Read namelist and allocate memory
  subroutine init_profiles

    integer :: ierr

    namelist /NAMOUT1D/ lstat, lprocblock, dtav, timeav

    dtav = dtav_glob
    timeav = timeav_glob

    ! Read namelist
    if (myid == 0) then
      open(ifnamopt, file=fname_options, status='old', iostat=ierr)
      read(ifnamopt, NAMOUT1D, iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMOUT1D')
      write(nnml_output, NAMOUT1D)
      close(ifnamopt)
    end if

    call d_mpi_bcast(lstat, 1, 0, comm3d, mpierr)
    call d_mpi_bcast(lprocblock, 1, 0, comm3d, mpierr)
    call d_mpi_bcast(dtav, 1, 0, comm3d, mpierr)
    call d_mpi_bcast(timeav, 1, 0, comm3d, mpierr)
    call d_mpi_bcast(nvar, 1, 0, comm3d, mpierr)

    idtav = int(dtav / tres, kind=longint)
    itimeav = int(timeav / tres, kind=longint)

    tnext = idtav + btime
    tnextwrite = itimeav + btime
    nsamples = int(itimeav / idtav)

    do_stats = .false.
    write_stats = .false.

    if (.not. lstat) return

    allocate(profiles(kmax,nvar))
    allocate(slab_average(1:k1))

    profiles = 0

    !$acc enter data copyin(profiles) create(slab_average) async

    if (lprocblock) then
      my_task_writes = .true. ! All MPI ranks write to a file
      fname = 'new-profiles.x'//cmyidx//'.y'//cmyidy//'.'//cexpnr//'.nc'
    else
      fname = 'new-profiles.'//cexpnr//'.nc'
      if (myid == 0) then ! Only the root MPI rank writes to a file
        my_task_writes = .true.
      else
        my_task_writes = .false.
      end if
    end if

    ! Initialize the NetCDF file
    if (my_task_writes) then
      call nctiminfo(tncname(1,:))
      call open_nc(fname, ncid, nrec, n3=kmax)

      if (nrec == 0) then
        call define_nc(ncid, 1, tncname)
        call writestat_dims_nc(ncid)
      end if

      call define_nc(ncid, nvar, ncname)
    end if

  end subroutine init_profiles


  subroutine exit_profiles

    if (my_task_writes .and. lstat) call exitstat_nc(ncid)

  end subroutine exit_profiles


  !> Bookkeeping routine
  !!
  !! Determines if profiles need to be sampled/written this timestep, and
  !! limits the time step.
  subroutine sample_profiles

    ! Reset switch
    do_stats = .false.
    write_stats = .false.

    if (.not. lstat) return
    if (rk3step /= 3) return

    if (timee < tnext .and. timee < tnextwrite) then
      dt_lim = minval([dt_lim, tnext - timee, tnextwrite - timee])
      return
    end if

    if (timee >= tnext) then
      do_stats = .true.
      tnext = tnext + idtav
    end if

    if (timee >= tnextwrite) then
      write_stats = .true.
      tnextwrite = tnextwrite + itimeav
    end if

  end subroutine sample_profiles

  subroutine sample_field(name, field)

    character(len=*), intent(in) :: name
    real(field_r),    intent(in) :: field(:,:,:)

    character(len=*), parameter :: routine = modname//'/sample_profile'

    integer :: k, idx
    integer :: nh

    if (do_stats) then

      ! TODO: a hash is probably more efficient here
      ! Find location in the list of profiles
      idx = find_index(name)

      if (idx == 0) then
        call finish(routine, 'profile ', trim(name), ' not found')
      end if

      ! A bit hacky maybe: figure out if the given field has ghost cells
      nh = (size(field, dim=1) - imax) / 2

      call slabavg(field, nh, slab_average, local=lprocblock)

      !$acc parallel loop vector default(present) async
      do k = 1, kmax
        profiles(k,idx) = profiles(k,idx) + slab_average(k)
      end do

    end if

  end subroutine sample_field

  subroutine sample_field_masked(name, field, mask)

    character(len=*), intent(in) :: name
    real(field_r),    intent(in) :: field(:,:,:)
    logical,          intent(in) :: mask(:,:,:)

    character(len=*), parameter :: routine = modname//'/sample_profile_masked'

    integer :: k, idx
    integer :: nh

    if (do_stats) then

      ! Find location in the list of profiles
      idx = find_index(name)

      if (idx == 0) then
        call finish(routine, 'profile ', trim(name), ' not found')
      end if

      ! A bit hacky maybe: figure out if the given field has ghost cells
      nh = (size(field, dim=1) - imax) / 2

      call slabavg(field, mask, nh, slab_average, local=lprocblock)

      !$acc parallel loop vector default(present) async
      do k = 1, kmax
        profiles(k,idx) = profiles(k,idx) + slab_average(k)
      end do

    end if

  end subroutine sample_field_masked

  subroutine write_profiles

    integer :: k, n

    if (write_stats) then

      !$acc parallel loop collapse(2) default(present)
      do n = 1, nvar
        do k = 1, kmax
          profiles(k,n) = profiles(k,n) / nsamples
        end do
      end do

      !$acc update host(profiles)

      if (my_task_writes) then
        call writestat_nc(ncid, 1, tncname, [rtimee], nrec, .true.)
        call writestat_nc(ncid, nvar, ncname, profiles(1:kmax,:), nrec, kmax)
      end if

      ! Reset averages
      !$acc parallel loop collapse(2) default(present) async
      do n = 1, nvar
        do k = 1, kmax
          profiles(k,n) = 0
        end do
      end do

    end if

  end subroutine write_profiles

end module modstat_profiles
