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
  use modtracer_type, only: T_tracer

  implicit none

  save
  public  :: inittracers, read_tracer_props, assign_tracer_props

  integer, parameter:: max_tracs  =  31 !<  Max. number of tracers that can be defined
  character(len=10),dimension(max_tracs) :: tracname_short = "NA" ! Short tracer name
  character(len=32),dimension(max_tracs) :: tracname_long  = "NA" ! Long tracer name
  character(len=10),dimension(max_tracs) :: tracer_unit = "NA" ! Unit of tracer
  real     :: molar_mass(max_tracs)       = -1.0 ! Molar mass of tracer (g mol-1)
  logical  :: tracer_is_emitted(max_tracs) = .false. ! Tracer is emitted (T/F)
  logical  :: tracer_is_reactive(max_tracs) = .false. ! Tracer is reactive (T/F)
  logical  :: tracer_is_deposited(max_tracs) = .false. ! Tracer is deposited (T/F)
  logical  :: tracer_is_photosynth(max_tracs) = .false. ! Tracer is photosynthesized (T/F)
  logical  :: tracer_is_microphys(max_tracs) = .false. ! Tracer is involved in cloud microphysics (T/F)

  integer  :: iname
  logical  :: ltracers = .false. ! tracer switch
  character(len = 6), dimension(200) :: &
              tracernames = (/ ('      ', iname=1, 200) /) ! list with scalar names,


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
    use modmpi,           only : myid, comm3d, mpierr, d_mpi_bcast

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
    call d_mpi_bcast(ltracers,             1, 0, comm3d, ierr)
    call d_mpi_bcast(tracernames(1:200), 200, 0, comm3d, ierr)

    allocate(tracer_prop(nsv), stat=ierr)
    if (ierr/=0) stop

    call read_tracer_props
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

  !! Read the list of available tracers from tracerdata.inp
  !! and their properties
  subroutine read_tracer_props

    use modglobal,  only : ifinput
    use modmpi,     only : myid

    implicit none

    integer   :: tdx, ierr, defined_tracers
    character(len=1500) :: readbuffer

    defined_tracers = 0
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

end module modtracers
