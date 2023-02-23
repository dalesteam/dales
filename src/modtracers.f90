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

  use modfields, only : svm, sv0, svp  
  use modglobal, only : nsv

  implicit none

  save
  public  :: inittracers, assign_tracer_props

  integer  :: iname
  logical  :: ltracers = .false. ! tracer switch
  character(len = 6), dimension(100) :: & 
              tracernames = (/ ('      ', iname=1, 100) /) ! list with scalar names,

  ! Data structure for tracers
  type T_tracer
  ! Fixed tracer properties
      ! Tracer name
      character(len=16) :: tracname
      ! Tracer long name
      character(len=64) :: traclong 
      ! Tracer unit
      character(len=16) :: unit     
      ! Tracer index in sv0, svm, svp
      integer           :: trac_idx 
      ! Boolean if tracer is emitted 
      logical           :: lemis
      ! Boolean if tracer is reactive
      logical           :: lreact
      ! Boolean if tracer is deposited
      logical           :: ldep
      ! Boolean if in A-gs
      logical           :: lags   
      ! Boolean if in cloud microphysics
      logical           :: lmicro     
      ! ! Static tracer properties:
      ! real :: diffusivity
  end type T_tracer

  ! Trace type for each tracer
  integer :: isv
  type(T_tracer), allocatable :: tracer_prop(:)

contains

  !> Initialize tracer definition.
  !!
  !! Read the namelist NAMTRACERS from the namoptions file, distribute
  !! the parameters to all processes and allocate the tracer (SV) arrays.
  subroutine inittracers
    ! read namelist    
    ! init tracer type

    use modglobal,        only : ifnamopt, fname_options, checknamelisterror
    use modmpi,           only : myid, comm3d, mpierr, mpi_logical, mpi_integer, my_real, mpi_character

    implicit none

    ! Auxiliary variables
    integer  :: ierr

    ! Namelist definition
    namelist /NAMTRACERS/ &
        ltracers, tracernames

    ! Read namelist
    if (myid == 0) then
        open(ifnamopt, file=fname_options, status='old', iostat=ierr)
        read(ifnamopt, NAMTRACERS, iostat=ierr)
        call checknamelisterror(ierr, ifnamopt, 'NAMTRACERS')
        write(6, NAMTRACERS)
        close(ifnamopt)
    end if

    ! Broadcast namelist values to all MPI tasks
    call MPI_BCAST(ltracers,             1, mpi_logical, 0, comm3d, mpierr)
    call mpi_bcast(tracernames(1:100), 100, mpi_character, 0, comm3d, ierr)

    allocate(tracer_prop(nsv), stat=ierr)
    if (ierr/=0) stop

    call assign_tracer_props

    ! do isv = 1, nsv
    !   write(6,"(A17, A6)") "species: ", trim(tracer_prop(isv) % tracname)
    ! enddo

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

  !! For each tracer in 'tracernames', get the properties as defined in 'modtracdata.f90'
  !! Fill the 'tracer_prop' data structure with these properties
  subroutine assign_tracer_props
    use modtracdata ! arrays with tracer properties

    implicit none

    ! First assign tracer index values
    ! They are equal to the sv index
    do isv=1,nsv
      ! tracer_prop(isv) % tracname = 'Dummy name'
      ! tracer_prop(isv) % traclong = 'Dummy long name'
      ! tracer_prop(isv) % unit     = 'Dummy unit'
      tracer_prop(isv) % trac_idx = isv
    end do

    do isv=1,nsv
      ! match species by short name and 
      ! look up species props in modtracdata arrays
      tracer_prop(isv) % tracname = trim(tracernames(isv))
      tracer_prop(isv) % traclong = trim(findval_character(tracernames(isv), tracname_short, &
                                      tracname_long, defltvalue='dummy longname'))  ! Default is 'dummy '
      ! TODO: write general find_character function to find string of of arbitrary lenght
      tracer_prop(isv) % unit     = trim(findval_character(tracernames(isv), tracname_short, &
                                      tracer_unit, defltvalue='dummy unit'))  ! Default is 'dummy unit'
      tracer_prop(isv) % lemis    = findval_logical(tracernames(isv), tracname_short, &
                                      tracer_is_emitted, defltvalue=.false.)  ! Default is False
      tracer_prop(isv) % lreact   = findval_logical(tracernames(isv), tracname_short, &
                                      tracer_is_reactive, defltvalue=.false.)  ! Default is False
      tracer_prop(isv) % ldep     = findval_logical(tracernames(isv), tracname_short, &
                                      tracer_is_deposited, defltvalue=.false.)  ! Default is False
      tracer_prop(isv) % lags     = findval_logical(tracernames(isv), tracname_short, &
                                      tracer_is_photosynthesized, defltvalue=.false.)  ! Default is False
      tracer_prop(isv) % lmicro   = findval_logical(tracernames(isv), tracname_short, &
                                      tracer_is_microphys, defltvalue=.false.)  ! Default is False

      write(*,*) 'tracer props ', tracer_prop(isv) % trac_idx, tracer_prop(isv) % tracname, tracer_prop(isv) % traclong, tracer_prop(isv) % unit, tracer_prop(isv) % lemis, tracer_prop(isv) % ldep, tracer_prop(isv) % lreact, tracer_prop(isv) % lags, tracer_prop(isv) % lmicro
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