!> \file modtracers.f90
!! Definitions and functions for passive and reactive tracers

!>
!!  \author Ruud Janssen, TNO
!!  \author Caspar Jungbacker, TU Delft
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

  use modglobal,      only: nsv, i1, ih, j1, jh, k1, kmax, cexpnr, iinput, &
                            input_ascii
  use modprecision,   only: field_r
  use modfields,      only: svm, sv0, svp, sv0av, svprof
  use modmpi,         only: myid, comm3d, d_mpi_bcast, print_info_stderr
  use go,             only: goSplitString_s
  use modstat_nc
  use utils

  implicit none

  private

  save

  character(len=*), parameter :: modname = 'modtracers'

  public :: inittracers
  public :: add_tracer
  public :: allocate_tracers
  public :: exittracers
  public :: tracer_profs_from_netcdf

  public :: tracer_prop
  public :: nsv_user

  type T_tracer
    character(len=16) :: tracname           !< Tracer name
    character(len=64) :: traclong           !< Tracer long name
    character(len=16) :: unit = '-'         !< Tracer unit
    real(field_r)     :: molar_mass = -999. !< Molecular mass of tracer (g mol-1)
    integer           :: trac_idx = 1       !< Tracer index in sv0, svm, svp
    logical           :: lemis = .false.    !< Boolean if tracer is emitted 
    logical           :: lreact = .false.   !< Boolean if tracer is reactive
    logical           :: ldep = .false.     !< Boolean if tracer is deposited
    logical           :: lags = .false.     !< Boolean if in A-gs
    logical           :: lmicro = .false.   !< Boolean if in cloud microphysics
    real(field_r)     :: wsvsurf = 0        !< Kinematic surface flux (- m/s)
  end type T_tracer

  integer, protected :: nsv_user !< Number of user-provided tracers

  type(T_tracer), allocatable, protected :: tracer_prop(:) !< List of tracers

contains

  !> Initialize tracer definition.
  subroutine inittracers
    character(len=*), parameter :: routine = modname//'::inittracers'

    character(len=128) :: file_profs
    logical            :: file_exists

    ! First, make sure that the tracer input file exists before calling
    ! tracer_props_from_xxx
    if (iinput == input_ascii) then
      file_profs = 'scalar.inp.'//cexpnr
    else
      file_profs = 'tracers.'//cexpnr//'.nc'
    end if

    inquire(file=file_profs, exist=file_exists)

    if (file_exists) then
      if (iinput == input_ascii) then
        call tracer_props_from_ascii(file_profs, 'tracerdata.inp')
      else
        call tracer_props_from_netcdf(file_profs)
      end if
    else
      call print_info_stderr(routine, trim(file_profs)//' not found')
      nsv_user = 0
    end if

  end subroutine inittracers

  !> \brief Define a new tracer
  !!
  !! \param name Short name of the tracer.
  !! \param long_name Full name of the tracer.
  !! \param unit Unit.
  !! \param molar_mass Molar mass (g mol^-1)
  !! \param lemis Tracer is emitted.
  !! \param lreact Tracer is reactive.
  !! \param ldep Tracer is deposited.
  !! \param lags Tracer is photosynthesized.
  !! \param lmicro Tracer is involved in cloud microphysics.
  !! \param wsvsurf Kinematic surface flux (- m/s).
  !! \note All tracers should be added before readinitfiles is called!
  subroutine add_tracer(name, long_name, unit, molar_mass, lemis, lreact, &
                        ldep, lags, lmicro, wsvsurf, isv)
    character(len=*), intent(in)            :: name
    character(len=*), intent(in),  optional :: long_name
    character(len=*), intent(in),  optional :: unit
    real(field_r),    intent(in),  optional :: molar_mass
    logical,          intent(in),  optional :: lemis
    logical,          intent(in),  optional :: lreact
    logical,          intent(in),  optional :: ldep
    logical,          intent(in),  optional :: lags
    logical,          intent(in),  optional :: lmicro
    real(field_r),    intent(in),  optional :: wsvsurf
    integer,          intent(out), optional :: isv

    character(len=*), parameter :: routine = modname//'::add_tracer'

    integer                     :: s
    character(len=1024)         :: message = ''
    type(T_tracer), allocatable :: tmp(:)

    ! Check if we have already allocated memory
    if (allocated(sv0)) then
      call print_info_stderr(routine, 'adding new tracers after memory is &
        &allocated is not allowed (tracer: '//trim(name)//')')
      error stop
    end if

    ! Check if the tracer already exists. If so, don't add a new one.
    if (nsv > 0) then
      do s = 1, nsv
        if (trim(to_lower(name)) == &
            trim(to_lower(tracer_prop(s) % tracname))) then
          write(message, '(a,a,a)') 'tracer ', trim(name), ' already defined'
          call print_info_stderr(routine, message)
          if(present(isv)) isv = s
          return
        end if
      end do
    end if

    nsv = nsv + 1

    ! Expand the list
    if (.not. allocated(tracer_prop)) then
      allocate(tracer_prop(nsv))
    else
      allocate(tmp(nsv))
      tmp(1:nsv-1) = tracer_prop(1:nsv-1)
      call move_alloc(tmp, tracer_prop)
    end if

    tracer_prop(nsv) % tracname = name
    tracer_prop(nsv) % trac_idx = nsv

    ! Long name. Use short name if not provided
    if (present(long_name)) then
      tracer_prop(nsv) % traclong = trim(long_name)
    else
      tracer_prop(nsv) % traclong = name
    end if

    if (present(unit)) tracer_prop(nsv) % unit = unit
    if (present(molar_mass)) tracer_prop(nsv) % molar_mass = molar_mass
    if (present(lemis)) tracer_prop(nsv) % lemis = lemis
    if (present(lreact)) tracer_prop(nsv) % lreact = lreact
    if (present(ldep)) tracer_prop(nsv) % ldep = ldep
    if (present(lags)) tracer_prop(nsv) % lags = lags
    if (present(lmicro)) tracer_prop(nsv) % lmicro = lmicro
    if (present(wsvsurf)) tracer_prop(nsv) % wsvsurf = wsvsurf

    if (present(isv)) isv = nsv

  end subroutine add_tracer

  !> Allocates all tracer fields
  subroutine allocate_tracers
    integer        :: isv
    type(T_tracer) :: tracer

    ! At this point, all tracers should be defined.

    ! Print tracer properties
    if (myid == 0) then
      write(6, '(a17,a17,a7,a9,a10,a11,a11)') &
        'Tracer           ', &
        'Unit             ', &
        'Index  ', &
        'Emitted  ', &
        'Reactive  ', &
        'Deposited  ', &
        'Surf. Flux '
      write(6, '(a)') repeat('-', 81)
      do isv = 1, nsv
        tracer = tracer_prop(isv)
        write(6, '(a,x,a,x,i3,4x,l3,6x,l3,7x,l3,8x,e10.4,x)') & ! Ugh
          tracer%tracname, &
          tracer%unit, &
          tracer%trac_idx, &
          tracer%lemis, &
          tracer%lreact, &
          tracer%ldep, &
          tracer%wsvsurf
      end do
    end if

    allocate(svm(2-ih:i1+ih,2-jh:j1+jh,k1,nsv), &
             sv0(2-ih:i1+ih,2-jh:j1+jh,k1,nsv), &
             svp(2-ih:i1+ih,2-jh:j1+jh,k1,nsv), &
             sv0av(k1,nsv), svprof(k1,nsv))

    svm(:,:,:,:) = 0
    sv0(:,:,:,:) = 0
    svp(:,:,:,:) = 0
    sv0av(:,:) = 0
    svprof(:,:) = 0

    !$acc enter data copyin(svm(2-ih:i1+ih,2-jh:j1+jh,1:k1,1:nsv), &
    !$acc&                  sv0(2-ih:i1+ih,2-jh:j1+jh,1:k1,1:nsv), &
    !$acc&                  svp(2-ih:i1+ih,2-jh:j1+jh,1:k1,1:nsv), &
    !$acc&                  sv0av(1:k1,1:nsv), svprof(1:k1,1:nsv))

  end subroutine allocate_tracers

  !> Deallocates all tracers fields
  subroutine exittracers

    !$acc exit data delete(svm(2-ih:i1+ih,2-jh:j1+jh,1:k1,1:nsv), &
    !$acc&                 sv0(2-ih:i1+ih,2-jh:j1+jh,1:k1,1:nsv), &
    !$acc&                 svp(2-ih:i1+ih,2-jh:j1+jh,1:k1,1:nsv), &
    !$acc&                 sv0av(1:k1,1:nsv), svprof(1:k1,1:nsv))

    if (nsv > 0) then
      deallocate(tracer_prop)
      deallocate(svm, sv0, svp, sv0av, svprof)
    end if

  end subroutine exittracers

  !> \brief Setup tracers from ASCII input files
  !!
  !! \param file_profiles Name of file containing profiles of tracers. From 
  !! this file, the number of tracers is determined.
  !! \param file_properties Name of file containing tracer attributes.
  subroutine tracer_props_from_ascii(file_profiles, file_properties)
    character(len=*), intent(in) :: file_profiles
    character(len=*), intent(in) :: file_properties

    character(len=*), parameter :: routine = &
      modname//'::tracer_props_from_ascii'
    integer,          parameter :: max_tracs = 100 !< Max. number of tracers that can be defined

    character(len=512) :: line
    character(len=7)   :: headers(max_tracs)
    integer            :: ierr
    integer            :: nheader
    integer            :: isv, n

    ! Buffers for tracer properties
    character(len=10) :: tracname_short(max_tracs) = 'NA'
    character(len=32) :: tracname_long(max_tracs)  = 'NA'
    character(len=10) :: tracer_unit(max_tracs) = 'NA'
    real(field_r)     :: molar_mass(max_tracs) = -1.0
    logical           :: tracer_is_emitted(max_tracs) = .false.
    logical           :: tracer_is_reactive(max_tracs) = .false.
    logical           :: tracer_is_deposited(max_tracs) = .false.
    logical           :: tracer_is_photosynth(max_tracs) = .false.
    logical           :: tracer_is_microphys(max_tracs) = .false.
    real(field_r)     :: wsvsurf(max_tracs) = 0.0

    open(1, file=file_profiles, status='old', iostat=ierr)

    if (ierr /= 0) then
      call print_info_stderr(routine, 'Error opening '//trim(file_profiles))
      error stop
    end if

    read(1, '(a512)') line
    read(1, '(a512)') line

    ! Determine the number of tracers from the header
    call goSplitString_s(line, nheader, headers, ierr, sep=' ')

    nsv_user = nheader - 1

    close(1)

    ! Read table of tracer properties, and put into the respective buffers
    open(1, file=file_properties, status='old', iostat=ierr)

    if (ierr /= 0) then
      call print_info_stderr(routine, 'Error opening '//trim(file_properties))
    else
      ierr = 0
      isv = 0
      do while (ierr == 0)
        read(1, '(A)', iostat=ierr) line

        if (ierr == 0) then ! So no end of file is encountered
          if (line(1:1)=='#') then
             cycle
          else
             isv = isv + 1
             read(line, *, iostat=ierr) tracname_short(isv), &
                                        tracname_long(isv), &
                                        tracer_unit(isv), &
                                        molar_mass(isv), &
                                        tracer_is_emitted(isv), &
                                        tracer_is_reactive(isv), &
                                        tracer_is_deposited(isv), &
                                        tracer_is_photosynth(isv), &
                                        tracer_is_microphys(isv), &
                                        wsvsurf(isv)
          end if
        end if
      end do
      close(1)
    end if


    ! For every tracer, find the properties
    do n = 2, nheader ! Skip the first column, containing the heights
      call add_tracer( &
        name=trim(headers(n)), &
        long_name=trim(findval(headers(n), tracname_short, tracname_long, &
                               defltvalue='dummy longname')), & ! Default is 'dummy '
        unit=trim(findval(headers(n), tracname_short, &
                    tracer_unit, defltvalue='dummy unit')), & ! Default is 'dummy unit'
        molar_mass=findval(headers(n), tracname_short, &
                     molar_mass, defltvalue=-999._field_r), & ! Default is -999.
        lemis=findval(headers(n), tracname_short, &
                tracer_is_emitted, defltvalue=.false.), & ! Default is False
        lreact=findval(headers(n), tracname_short, &
                 tracer_is_reactive, defltvalue=.false.), & ! Default is False
        ldep=findval(headers(n), tracname_short, &
               tracer_is_deposited, defltvalue=.false.), & ! Default is False
        lags=findval(headers(n), tracname_short, &
               tracer_is_photosynth, defltvalue=.false.), & ! Default is False
        lmicro=findval(headers(n), tracname_short, &
                 tracer_is_microphys, defltvalue=.false.), & ! Default is False
        wsvsurf=findval(headers(n), tracname_short, &
                  wsvsurf, defltvalue=0.0_field_r) &
      )
    end do
    
  end subroutine tracer_props_from_ascii

  !> \brief Read tracer properties from tracers.XXX.nc
  !!
  !! \param filename Name of the input file.
  subroutine tracer_props_from_netcdf(filename)
    character(len=*), intent(in)  :: filename

    character(len=*), parameter :: routine = &
      modname//'::tracer_props_from_netcdf'

    integer              :: ncid, nvars, ivar
    logical              :: file_exists
    integer, allocatable :: varids(:)

    ! Tracer attributes
    character(len=NF90_MAX_NAME) :: name
    character(len=32)            :: long_name
    character(len=16)            :: unit
    real(field_r)                :: molar_mass
    logical                      :: lemis, lreact, ldep, lags

    call nchandle_error(nf90_open(filename, NF90_NOWRITE, ncid))
    call nchandle_error(nf90_inquire(ncid, nVariables=nvars))

    allocate(varids(nvars))

    call nchandle_error(nf90_inq_varids(ncid, nvars, varids))

    nsv_user = nvars

    do ivar = 1, nvars
      call nchandle_error(nf90_inquire_variable(ncid, varids(ivar), name=name))

      ! Read attributes
      call read_nc_attribute(ncid, varids(ivar), 'long_name', long_name, &
                             default='dummy')
      call read_nc_attribute(ncid, varids(ivar), 'unit', unit, default='---')
      call read_nc_attribute(ncid, varids(ivar), 'molar_mass', molar_mass, &
                             default=-999._field_r)
      call read_nc_attribute(ncid, varids(ivar), 'lemis', lemis, &
                             default=.false.)
      call read_nc_attribute(ncid, varids(ivar), 'lreact', lreact, &
                             default=.false.)
      call read_nc_attribute(ncid, varids(ivar), 'ldep', ldep, default=.false.)
      call read_nc_attribute(ncid, varids(ivar), 'lags', lags, default=.false.)

      ! Setup tracer
      call add_tracer(trim(name), long_name=trim(long_name), unit=unit, &
                      molar_mass=molar_mass, lemis=lemis, lreact=lreact, &
                      ldep=ldep, lags=lags, lmicro=.false.)
    end do

    deallocate(varids)

  end subroutine tracer_props_from_netcdf

  !> \brief Read initial profiles for tracers
  !!
  !! \param filename Name of the input file.
  !! \param tracers List of tracers to read input for.
  !! \param nsv Number of tracers.
  !! \param svprof 2D array (z,s) to place initial profiles in.
  subroutine tracer_profs_from_netcdf(filename, tracers, svprof)
    character(len=*),   intent(in)  :: filename
    type(T_tracer),     intent(in)  :: tracers(:)
    real(field_r),      intent(out) :: svprof(:,:)

    integer :: ncid, ivar

    call nchandle_error(nf90_open(filename, NF90_NOWRITE, ncid))

    do ivar = 1, nsv
      call read_nc_field(ncid, trim(tracers(ivar) % tracname), svprof(:,ivar), &
                         start=1, count=kmax, fillvalue=0._field_r)
    end do

  end subroutine tracer_profs_from_netcdf

end module modtracers
