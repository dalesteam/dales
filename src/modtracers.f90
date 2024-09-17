#define DEBUG
!> \file modtracers.f90
!! Definitions and functions for passive and reactive tracers

!>
!!  \author Ruud Janssen, TNO
!  This file is part of DALES.
!
! DALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 1993-2023 Delft University of Technology, Wageningen
!  University, Utrecht University, KNMI, TNO
!

module modtracers

  use modglobal, only : nsv
  use modprecision, only: field_r
  use modtracer_type, only: T_tracer

  implicit none
  save
  
  public :: inittracers
  public :: read_tracer_props
  public :: assign_tracer_props
  public :: add_tracer

  public :: defined_tracers
  public :: tracer_prop
  public :: findval

  integer, parameter:: max_tracs  =  31 !<  Max. number of tracers that can be defined
  character(len=10),dimension(max_tracs) :: tracname_short = "NA" ! Short tracer name
  character(len=32),dimension(max_tracs) :: tracname_long  = "NA" ! Long tracer name
  character(len=10),dimension(max_tracs) :: tracer_unit = "NA" ! Unit of tracer
  real(field_r)     :: molar_mass(max_tracs)       = -1.0 ! Molar mass of tracer (g mol-1)
  logical  :: tracer_is_emitted(max_tracs) = .false. ! Tracer is emitted (T/F)
  logical  :: tracer_is_reactive(max_tracs) = .false. ! Tracer is reactive (T/F)
  logical  :: tracer_is_deposited(max_tracs) = .false. ! Tracer is deposited (T/F)
  logical  :: tracer_is_photosynth(max_tracs) = .false. ! Tracer is photosynthesized (T/F)
  logical  :: tracer_is_microphys(max_tracs) = .false. ! Tracer is involved in cloud microphysics (T/F)
  logical  :: tracer_is_aerosol(max_tracs) = .false.
  character(3) :: tracer_mode(max_tracs) = "NA"
  real(field_r) :: tracer_kappa(max_tracs) = -1.0
  real(field_r) :: tracer_rho(max_tracs) = -1.0

  integer :: defined_tracers

  integer  :: iname
  logical  :: ltracers = .false. ! tracer switch
  character(len = 6), dimension(200) :: &
              tracernames = (/ ('      ', iname=1, 200) /) ! list with scalar names,
  logical  :: lnc = .true. ! Read properties from NetCDF file

  ! Trace type for each tracer
  integer :: isv
  type(T_tracer), allocatable :: tracer_prop(:)

  interface findval
    module procedure findval_real
    module procedure findval_integer
    module procedure findval_logical
    module procedure findval_character
  end interface findval

contains


  !> Initialize tracer definition.
  !!
  !! Read the namelist NAMTRACERS from the namoptions file, distribute
  !! the parameters to all processes and allocate the tracer (SV) arrays.
  subroutine inittracers
    ! read namelist    
    ! init tracer type

    use modglobal,        only : ifnamopt, fname_options, checknamelisterror
    use modmpi,           only : myid, comm3d, mpierr, d_mpi_bcast

    implicit none

    ! Auxiliary variables
    integer  :: ierr

    ! Namelist definition
    namelist /NAMTRACERS/ &
        ltracers, tracernames, lnc

    ! Read namelist
    if (myid == 0) then
        open(ifnamopt, file=fname_options, status='old', iostat=ierr)
        read(ifnamopt, NAMTRACERS, iostat=ierr)
        call checknamelisterror(ierr, ifnamopt, 'NAMTRACERS')
        write(6, NAMTRACERS)
        close(ifnamopt)
    end if

    ! Broadcast namelist values to all MPI tasks
    call d_mpi_bcast(ltracers,             1, 0, comm3d, ierr)
    call d_mpi_bcast(tracernames(1:200), 200, 0, comm3d, ierr)
    call d_mpi_bcast(lnc, 1, 0, comm3d, ierr)

    allocate(tracer_prop(nsv), stat=ierr)
    if (ierr/=0) stop

    if (lnc) then
      call read_tracer_props_from_netcdf
    else
      call read_tracer_props
      call assign_tracer_props
    end if
    
  end subroutine inittracers

  !
  ! Cleanup (deallocate) the tracers
  ! 
  subroutine exittracers
    ! use ..., only :
    implicit none

    if (.not. ltracers) return

    ! Tracer properties:
    deallocate(tracer_prop)

  end subroutine exittracers

  subroutine add_tracer(name, long_name, unit, molar_mass, lemis, lreact, &
                        ldep, lags, lmicro, laero, rho, kappa, mode, iINC, &
                        iINR, isv)

    use modmpi, only : myid

    character(len=*), intent(in)           :: name
    character(len=*), intent(in), optional :: long_name
    character(len=*), intent(in), optional :: unit
    character(len=3), intent(in), optional :: mode
    real(field_r),    intent(in), optional :: molar_mass
    real(field_r),    intent(in), optional :: rho
    real(field_r),    intent(in), optional :: kappa
    logical,          intent(in), optional :: lemis
    logical,          intent(in), optional :: lreact
    logical,          intent(in), optional :: ldep
    logical,          intent(in), optional :: lags
    logical,          intent(in), optional :: lmicro
    logical,          intent(in), optional :: laero
    integer,          intent(in), optional :: iINC
    integer,          intent(in), optional :: iINR
    integer,          intent(out)          :: isv

    defined_tracers = defined_tracers + 1

    if (defined_tracers > max_tracs) then
      if (myid == 0) write(*,*) "Maximum number of tracers reached!"
      stop
    end if
     
    ! Check if tracer already exists

    isv = defined_tracers
    
    tracer_prop(isv)%tracname = name
    tracer_prop(isv)%trac_idx = isv
    if (present(long_name)) tracer_prop(isv)%traclong = trim(long_name)
    if (present(unit)) tracer_prop(isv)%unit = unit
    if (present(mode)) tracer_prop(isv)%mode = mode
    if (present(molar_mass)) tracer_prop(isv)%molar_mass = molar_mass
    if (present(rho)) tracer_prop(isv)%kappa = kappa
    if (present(kappa)) tracer_prop(isv)%rho = rho
    if (present(lemis)) tracer_prop(isv)%lemis = lemis
    if (present(lreact)) tracer_prop(isv)%lreact = lreact
    if (present(ldep)) tracer_prop(isv)%ldep = ldep
    if (present(lags)) tracer_prop(isv)%lags = lags
    if (present(lmicro)) tracer_prop(isv)%lmicro = lmicro
    if (present(laero)) tracer_prop(isv)%laero = laero
    if (present(iINC)) tracer_prop(isv)%iINC = iINC
    if (present(iINR)) tracer_prop(isv)%iINR = iINR
        
  end subroutine add_tracer

  subroutine read_tracer_props_from_netcdf

    use go,                   only : goSplitString_s
    use modglobal,            only : cexpnr
    use modmicrodata_aerosol, only : modenames, modes
    use modmpi,               only : myid
    use modnetcdf

    integer              :: ncid, ierr
    integer, allocatable :: varid(:)
    integer              :: itrac, ntrac, imod, nmodes, tdx, defined_tracers, idx
    logical              :: laero
    character(27)        :: modes_s
    character(3)         :: my_modes(9)
    integer              :: defined_aerosols = 0

    character(len=16) :: tracname
    character(len=32) :: traclong
    character(len=16) :: unit     
    real              :: the_molar_mass
    integer           :: trac_idx
    logical           :: lemis
    logical           :: lreact
    logical           :: ldep
    logical           :: lags   
    logical           :: lmicro     
    real(field_r)     :: rho, kappa
    
    if (myid == 0) then
      call check(nf90_open("tracers."//cexpnr//".nc", nf90_nowrite, ncid), __FILE__, __LINE__)

      call check(nf90_inquire(ncid, nVariables=ntrac), __FILE__, __LINE__)
      
      allocate(varid(ntrac))

      ! Get the ID's of all tracers
      call check(nf90_inq_varids(ncid, ntrac, varid), __FILE__, __LINE__)
      
      tdx = 0
      
      do itrac = 1, ntrac
        call check(nf90_inquire_variable(ncid, varid(itrac), name=tracname), __FILE__, __LINE__)

        ! Try to find other properties, using defaults if they're not found
        call getattr(ncid, varid(itrac), "long_name", traclong, "dummy")
        call getattr(ncid, varid(itrac), "unit", unit, "-")
        call getattr(ncid, varid(itrac), "molar_mass", the_molar_mass, -999.)
        call getattr(ncid, varid(itrac), "is_emitted", lemis, .false.)
        call getattr(ncid, varid(itrac), "is_reactive", lreact, .false.)
        call getattr(ncid, varid(itrac), "is_deposited", ldep, .false.)
        call getattr(ncid, varid(itrac), "is_photosynth", lags, .false.)
        call getattr(ncid, varid(itrac), "is_microphys", lmicro, .false.)
        call getattr(ncid, varid(itrac), "kappa", kappa, -1.0_field_r)
        call getattr(ncid, varid(itrac), "rho", rho, -1.0_field_r)
        call getattr(ncid, varid(itrac), "is_aerosol", laero, .false.)
          
        if (laero) then
          defined_aerosols = defined_aerosols + 1

          ! Figure out the modes
          call check(nf90_get_att(ncid, varid(itrac), "modes", modes_s), __FILE__, __LINE__)
          call goSplitString_s(trim(modes_s), nmodes, my_modes, ierr, sep=',')

          ! Check if the provided modes are valid
          do imod = 1, nmodes
            idx = findloc(modenames, my_modes(imod), dim=1)
            if (idx == 0) then
              write(*,*) "Mode "//trim(my_modes(imod))//" is not valid!" 
              stop
            end if
          end do
          
          ! Now add a tracer for each mass concentration, we will read this later from init.XXX.nc
          do imod = 1, nmodes
            call add_tracer(trim(tracname)//"_"//my_modes(imod), trim(traclong)//", "//trim(my_modes(imod))//" mode", &
              &             unit, the_molar_mass, lemis, lreact, ldep, lags, lmicro, laero, rho, &
              &             kappa, my_modes(imod), defined_aerosols, defined_aerosols, tdx)
          end do
          
          ! Manually add in-cloud and in-rain modes
          call add_tracer(trim(tracname)//"_"//"inc", trim(traclong)//", in-cloud mode", &
            &             unit, the_molar_mass, lemis, lreact, ldep, lags, lmicro, laero, rho, &
            &             kappa, "inc", isv=tdx)
          call add_tracer(trim(tracname)//"_"//"inr", trim(traclong)//", in-rain mode", &
            &             unit, the_molar_mass, lemis, lreact, ldep, lags, lmicro, laero, rho, &
            &             kappa, "inr", isv=tdx)
            
        else ! laero
          call add_tracer(tracname, trim(traclong), unit, the_molar_mass, lemis, lreact, ldep, lags, &
            &             lmicro, laero, rho, kappa, "N/A", isv=tdx)          
        end if ! laero
        
      end do ! ntrac

      deallocate(varid)
      defined_tracers = tdx

      ! Print properties
      do itrac = 1, defined_tracers
        call tracer_prop(itrac) % print_properties
      end do

    end if ! myid == 0

  end subroutine read_tracer_props_from_netcdf

  !! Read the list of available tracers from tracerdata.inp
  !! and their properties
  subroutine read_tracer_props

    use modglobal,  only : ifinput, cexpnr
    use modmpi,     only : myid
    use modnetcdf

    implicit none

    integer :: tdx, ierr, defined_tracers
    character(len=1500) :: readbuffer

    defined_tracers = 0

    ! CJ: Kept for compatibility, for TNO people :)
    open (ifinput,file='tracerdata.inp')
    ierr = 0
    do while (ierr == 0)
      read(ifinput, '(A)', iostat=ierr) readbuffer

      if (ierr == 0) then   !So no end of file is encountered
        if (readbuffer(1:1)=='#') then
          if (myid == 0)   print *, trim(readbuffer)
        else
          if (myid == 0)   print *, trim(readbuffer)
          defined_tracers = defined_tracers + 1
          tdx = defined_tracers
          read(readbuffer, *, iostat=ierr) tracname_short(tdx), tracname_long(tdx), tracer_unit(tdx), molar_mass(tdx), &
                                           tracer_is_emitted(tdx), tracer_is_reactive(tdx), tracer_is_deposited(tdx), &
                                           tracer_is_photosynth(tdx), tracer_is_microphys(tdx)

          if ( (molar_mass(tdx) .lt. 0) .and. (molar_mass(tdx) .gt. -1.1) ) then
            if (myid == 0) stop "MODTRACERS: a molar mass value is not set in the tracer input file"
          end if
        endif
      endif
    enddo

    close(ifinput)

  end subroutine read_tracer_props


  !! For each tracer in 'tracernames', get the properties as defined in 'tracerdata.inp' and
  !! read by 'read_tracer_props'
  !! Fill the 'tracer_prop' data structure with these properties
  subroutine assign_tracer_props
    !!! ! use modtracdata ! arrays with tracer properties

    implicit none

    ! CJ: Why not do this in read_tracer_props?

    do isv=1,nsv
      ! First assign tracer index values. They are equal to the sv index by default
      tracer_prop(isv) % trac_idx = isv

      ! match species by short name and 
      ! look up species props in modtracdata arrays
      tracer_prop(isv) % tracname = trim(tracernames(isv))
      tracer_prop(isv) % traclong = trim(findval_character(tracernames(isv), tracname_short, &
                                      tracname_long, defltvalue='dummy longname'))  ! Default is 'dummy '
      tracer_prop(isv) % unit     = trim(findval_character(tracernames(isv), tracname_short, &
                                      tracer_unit, defltvalue='dummy unit'))  ! Default is 'dummy unit'
      tracer_prop(isv) % molar_mass = findval_real(tracernames(isv), tracname_short, &
                                      molar_mass, defltvalue=-999.)  ! Default is -999.
      tracer_prop(isv) % lemis    = findval_logical(tracernames(isv), tracname_short, &
                                      tracer_is_emitted, defltvalue=.false.)  ! Default is False
      tracer_prop(isv) % lreact   = findval_logical(tracernames(isv), tracname_short, &
                                      tracer_is_reactive, defltvalue=.false.)  ! Default is False
      tracer_prop(isv) % ldep     = findval_logical(tracernames(isv), tracname_short, &
                                      tracer_is_deposited, defltvalue=.false.)  ! Default is False
      tracer_prop(isv) % lags     = findval_logical(tracernames(isv), tracname_short, &
                                      tracer_is_photosynth, defltvalue=.false.)  ! Default is False
      tracer_prop(isv) % lmicro   = findval_logical(tracernames(isv), tracname_short, &
                                      tracer_is_microphys, defltvalue=.false.)  ! Default is False

      ! write(*,*) 'tracer props ', tracer_prop(isv) % trac_idx, tracer_prop(isv) % tracname, tracer_prop(isv) % traclong, tracer_prop(isv) % unit, tracer_prop(isv) % lemis, tracer_prop(isv) % ldep, tracer_prop(isv) % lreact, tracer_prop(isv) % lags, tracer_prop(isv) % lmicro
    end do

  end subroutine assign_tracer_props

  ! > Find a value in an array of values based on a key in an array of keys
  ! !
  ! ! Lookup a value in one array (`values`) at the index position of `key` in the array `keys`.
  ! ! Provide the optional default value (`defltvalue`), that is returned when `key` isn't found.
  ! ! When no default is given, zero is returned.
  ! !
  ! ! @param[in] key The key to be looked up
  ! ! @param[in] keys The array of keys
  ! ! @param[in] values The array of values
  ! ! @param[in] defltvalue The default value (optional, default 0.0)
  ! ! @return The value corresponding to the key
  pure real function findval_real(key, keys, values, defltvalue)
    implicit none
    character(*), intent(in) :: key 
    character(*), intent(in) :: keys(:)
    real, intent(in) :: values(:)
    real, intent(in), optional :: defltvalue
    character(6) :: fkey
    integer :: idx(1)

    fkey = trim(key)

    idx = findloc(keys, fkey)
    if (idx(1) /= 0 .and. idx(1) <= size(values)) then
      findval_real = values(idx(1))
    else
      if (present(defltvalue)) then
        findval_real = defltvalue 
      else
        findval_real = 0.0 
      end if
    end if
  end function findval_real

  ! > Find a value in an array of values based on a key in an array of keys
  ! !
  ! ! Lookup a value in one array (`values`) at the index position of `key` in the array `keys`.
  ! ! Provide the optional default value (`defltvalue`), that is returned when `key` isn't found.
  ! ! When no default is given, zero is returned.
  ! !
  ! ! @param[in] key The key to be looked up
  ! ! @param[in] keys The array of keys
  ! ! @param[in] values The array of values
  ! ! @param[in] defltvalue The default value (optional, default -1)
  ! ! @return The value corresponding to the key
  pure integer function findval_integer(key, keys, values, defltvalue)
    implicit none
    character(*), intent(in) :: key 
    character(*), intent(in) :: keys(:)
    integer, intent(in) :: values(:)
    integer, intent(in), optional :: defltvalue
    character(6) :: fkey
    integer :: idx(1)

    fkey = trim(key)

    idx = findloc(keys, fkey)
    if (idx(1) /= 0 .and. idx(1) <= size(values)) then
      findval_integer = values(idx(1))
    else
      if (present(defltvalue)) then
        findval_integer = defltvalue 
      else
        findval_integer = -1 
      end if
    end if
  end function findval_integer

  ! > Find a value in an array of values based on a key in an array of keys
  ! !
  ! ! Lookup a value in one array (`values`) at the index position of `key` in the array `keys`.
  ! ! Provide the optional default value (`defltvalue`), that is returned when `key` isn't found.
  ! ! When no default is given, zero is returned.
  ! !
  ! ! @param[in] key The key to be looked up
  ! ! @param[in] keys The array of keys
  ! ! @param[in] values The array of values
  ! ! @param[in] defltvalue The default value (optional, default .False.)
  ! ! @return The value corresponding to the key
  pure logical function findval_logical(key, keys, values, defltvalue)
    implicit none
    character(*), intent(in) :: key 
    character(*), intent(in) :: keys(:)
    logical, intent(in) :: values(:)
    logical, intent(in), optional :: defltvalue
    character(6) :: fkey
    integer :: idx(1)

    fkey = trim(key)

    idx = findloc(keys, fkey)
    if (idx(1) /= 0 .and. idx(1) <= size(values)) then
      findval_logical = values(idx(1))
    else
      if (present(defltvalue)) then
        findval_logical = defltvalue 
      else
        findval_logical = .False. 
      end if
    end if
  end function findval_logical

  ! > Find a value in an array of values based on a key in an array of keys
  ! !
  ! ! Lookup a value in one array (`values`) at the index position of `key` in the array `keys`.
  ! ! Provide the optional default value (`defltvalue`), that is returned when `key` isn't found.
  ! ! When no default is given, zero is returned.
  ! !
  ! ! @param[in] key The key to be looked up
  ! ! @param[in] keys The array of keys
  ! ! @param[in] values The array of values
  ! ! @param[in] defltvalue The default value (optional, default .False.)
  ! ! @return The value corresponding to the key
  pure character(32) function findval_character(key, keys, values, defltvalue)
    implicit none
    character(*), intent(in) :: key 
    character(*), intent(in) :: keys(:)
    character(*), intent(in) :: values(:)
    character(6), intent(in), optional :: defltvalue
    character(6) :: fkey
    integer :: idx(1)

    fkey = trim(key)
    idx = findloc(keys, fkey)
    if (idx(1) /= 0 .and. idx(1) <= size(values)) then
      findval_character = trim(values(idx(1)))
    else
      if (present(defltvalue)) then
        findval_character = defltvalue 
      else
        findval_character = 'dummy'
      end if
    end if
  end function findval_character

end module modtracers
