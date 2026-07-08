!> \file modemission.f90
!!  (Anthropogenic) emissions

!>
!!  \author Marco de Bruine, VU
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
!  Copyright 1993-2026 Delft University of Technology, Wageningen
!  University, Utrecht University, KNMI
!

module modemission
use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
use modprecision, only: field_r
use ieee_arithmetic, only: ieee_is_nan
use modemisdata
use modtracers,       only : tracer_prop
use modlogging, only: finish, warning, message

implicit none

character(len=*), parameter :: modname = 'modemission'

contains

  subroutine initemission

    use modglobal,    only : i2, j2,kmax, nsv, ifnamopt, fname_options, checknamelisterror
    use modmpi,       only : myid, comm3d, d_mpi_bcast
    use moddatetime,  only : datex, prevday, nextday
    use fortran_support, only : nnml_output
    use modlogging, only: profile_output

    implicit none

    character(len=*), parameter :: routine = modname//"/initemission"

    ! Auxiliary variables
    integer :: ierr, l

    ! --- Read & broadcast namelist EMISSION -----------------------------------
    ! namelist/NAMEMISSION/ l_emission, kemis, svskip, emisnames, svco2sum
    namelist/NAMEMISSION/ l_emission, l_points, explicit_plume_rise, kemis, nemis, emisnames, l_scale, scalefactor

    if (myid == 0) then

      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMEMISSION,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMEMISSION')
      write(nnml_output, NAMEMISSION)
      close(ifnamopt)

    endif

    call d_mpi_bcast(l_emission,           1,  0, comm3d, ierr)
    call d_mpi_bcast(l_points,             1,  0, comm3d, ierr)
    call d_mpi_bcast(explicit_plume_rise,  1,  0, comm3d, ierr)
    call d_mpi_bcast(kemis,                1,  0, comm3d, ierr)
    call d_mpi_bcast(emisnames(1:100),   100,  0, comm3d, ierr)
    call d_mpi_bcast(nemis,                1,  0, comm3d, ierr)
    call d_mpi_bcast(l_scale,              1,  0, comm3d, ierr)
    call d_mpi_bcast(scalefactor(1:100), 100,  0, comm3d, ierr)

    ! -- Interaction with AGs   ----------------------------------------------------
    if (.not. (l_emission)) return

#if defined(DALES_GPU)
    call finish(routine, "emissions are not supported on GPU")
#endif
    !allocate(co2fields(nsv))

    !co2fields = 0
    ! co2fields(svskip+1:nsv) = index(emisnames(1:nsv-svskip), "co2")
    !co2fields = index(emisnames, "co2")
    
    !svco2sum = findloc(emisnames, value = "co2sum", dim = 1)

    ! svco2ags = findloc(emisnames(1:nsv-svskip), value = "co2ags", dim = 1)
    ! svco2ags = svco2ags + svskip
    !svco2ags = findloc(emisnames, value = "co2ags", dim = 1)

    ! svco2veg = findloc(emisnames(1:nsv-svskip), value = "co2veg", dim = 1)
    ! svco2veg = svco2veg + svskip
    !svco2veg = findloc(emisnames, value = "co2veg", dim = 1)
    
    
    ! Find CO2 index and set conversions
    do l = 1, nsv
        if (trim(tracer_prop(l)%tracname) == 'co2ags') then                       
            svco2ags= tracer_prop(l)%trac_idx
        else if (trim(tracer_prop(l)%tracname) == 'co2veg') then 
            svco2veg= tracer_prop(l)%trac_idx
        else if (trim(tracer_prop(l)%tracname) == 'co2sum') then
            svco2sum = tracer_prop(l)%trac_idx
        end if 
    enddo


    if (myid == 0) then
      write(profile_output,*) 'modemission: co2fields (scalar fields with CO2 0=no, 1=yes)'
      write(profile_output,*) co2fields
      write(profile_output,*) 'modemission: svco2ags (scalar field number for AGS emissions)'
      write(profile_output,*) svco2ags
      write(profile_output,*) 'modemission: svco2veg (scalar field number for AGS emissions)'
      write(profile_output,*) svco2veg
      write(profile_output,*) 'modemission: svco2sum (scalar field number for AGS emissions)'
      write(profile_output,*) svco2sum
      write(profile_output,*) 'number of emitted species'
      write(profile_output,*) nemis
    endif

    ! --- Local pre-calculations and settings
    if (kemis == -1) kemis = kmax

    ! --- Read emission files for first time step ----------------------------------

    ! Two hourly emission fields are loaded at all times:
    ! (1) before model time,   t_field < t_model, "in-the-past"
    ! (2) ahead of model time, t_field > t_model, "in-the-future"
    ! allocate(emisfield(i2, j2, kemis, svskip+1:nsv, 2))
    allocate(emisfield(i2, j2, kemis, nemis, 2))


    if (datex(5) >= 30) then
      call reademission(    datex(1),   datex(2),   datex(3),   datex(4), emisfield(:,:,:,:,1))

      if (datex(4) == 23) then
        call reademission(nextday(1), nextday(2), nextday(3),          0, emisfield(:,:,:,:,2))
      else
        call reademission(  datex(1),   datex(2),   datex(3), datex(4)+1, emisfield(:,:,:,:,2))
      endif

    else
      call reademission(    datex(1),   datex(2),   datex(3),   datex(4), emisfield(:,:,:,:,2))

      if (datex(4) == 0) then
        call reademission(prevday(1), prevday(2), prevday(3),         23, emisfield(:,:,:,:,1))
      else
        call reademission(  datex(1),   datex(2),   datex(3), datex(4)-1, emisfield(:,:,:,:,1))
      endif

    endif

    if (l_points) then
    
        ! Ensure point_sources is allocated on all ranks
        if (.not. allocated(point_sources)) then
            allocate(point_sources(nsv))
        endif
    
        if (myid == 0) then
            ! First, inquire how many points exist for each time/tracer
            call inquirepoints(datex(1), datex(2), datex(3), datex(4))
            if (datex(5) >= 30) then
                if (datex(4) == 23) then
                    call inquirepoints(nextday(1), nextday(2), nextday(3), 0)
                else
                    call inquirepoints(datex(1), datex(2), datex(3), datex(4)+1)
                endif
            else
                if (datex(4) == 0) then
                    call inquirepoints(prevday(1), prevday(2), prevday(3), 23)
                else
                    call inquirepoints(datex(1), datex(2), datex(3), datex(4)-1)
                endif
            endif
        
            ! Root reads point source data
            if (datex(5) >= 30) then
                call readpoints(datex(1), datex(2), datex(3), datex(4), 1)  ! t1 = current
                if (datex(4) == 23) then
                    call readpoints(nextday(1), nextday(2), nextday(3), 0, 2)
                else
                    call readpoints(datex(1), datex(2), datex(3), datex(4)+1, 2)
                endif
            else
                call readpoints(datex(1), datex(2), datex(3), datex(4), 2)
                if (datex(4) == 0) then
                    call readpoints(prevday(1), prevday(2), prevday(3), 23, 1)
                else
                    call readpoints(datex(1), datex(2), datex(3), datex(4)-1, 1)
                endif
            endif
        endif

        ! Distribute point source data to all ranks
        do l = 1, nsv
            if (tracer_prop(l)%lemis) then
                call distributepoints(l)
            endif
        end do
    endif

  end subroutine initemission

  subroutine reademission(iyear, imonth, iday, ihour, emisfield)

  ! ----------------------------------------------------------------------
  ! Reading of emission files
  ! Multiple/all tracers
  ! ----------------------------------------------------------------------

    use netcdf
    use modmpi,      only : myidx, myidy
    use modglobal,   only : i1, j1, i2, j2, imax, jmax, nsv

    implicit none

    character(len=*), parameter :: routine = modname//'/reademission'

    integer, intent(in)  :: iyear, imonth, iday, ihour
    ! real, intent(out)    :: emisfield(i2, j2, kemis, 1+svskip:nsv)
    real, intent(out)    :: emisfield(i2, j2, kemis, nemis)

    integer              :: ncid, varid
    integer              :: isv, iem
    integer              :: unitlength = 64
    character(len=12)    :: sdatetime
    character (len = 80) :: unit

    ! Create string from given date
    write(sdatetime, "(I0.4,2I0.2,2I0.2)") iyear, imonth, iday, ihour, 0

    call message(routine, "Reading emission: ", sdatetime)

    iem = 1
    do isv = 1, nsv
      if (tracer_prop(isv)%lemis) then
        ! check tracer unit
        ! give warning when emission file is not available for a species which is emitted
        if (iem > nemis) then
          call warning(routine, "More emitted species than declared in NAMEMISSION: ", iem, " > ", nemis)
        endif
        call message(routine, "Reading tracer: ", tracer_prop(isv)%trac_idx, " ", trim(tracer_prop(isv)%tracname))
        call check( nf90_open( 'emissions/'//trim(tracer_prop(isv)%tracname)//'_emis_'//sdatetime//'_3d.nc', NF90_NOWRITE, ncid))
        call check( nf90_inq_varid( ncid, tracer_prop(isv)%tracname, varid) )
        call check( nf90_get_var  ( ncid, varid, emisfield(2:i1,2:j1,1:kemis,iem), &
                                    start = (/1 + myidx * imax, 1 + myidy * jmax, 1, 1/), &
                                    count = (/imax, jmax, kemis, 1/) ) )
        call check( nf90_inquire_attribute(ncid, varid, 'units',  len = unitlength) )
        call check( nf90_get_att( ncid, varid, 'units', unit) )
        call check( nf90_close( ncid ) )
        ! write(6,"(A22, A22)") "Reading tracer unit: ", trim(unit)
        if ( trim(unit) /= 'kg hour-1' ) then
          !!! TODO: make this an ERROR after updating the emission pre-processor
          call warning(routine, "WARNING: emission units do not match: " , trim(unit), " /= kg hour-1")
        endif
        iem = iem + 1
      else
        call warning(routine,"Tracer not emitted: ", tracer_prop(isv)%trac_idx, trim(tracer_prop(isv)%tracname))
      endif
    end do

  contains

  subroutine check(status)
    integer, intent(in) :: status

    if(status /= nf90_noerr) then
      call finish(routine, 'NetCDF error in modemission. See outputfile for more information. Error ', trim(nf90_strerror(status)))
    end if
  end subroutine check

  end subroutine reademission

  subroutine emission
  ! ----------------------------------------------------------------------
  ! Read appropriate emission fields, interpolate and transfer to svp
  !
  ! NOTES
  ! 1. Emission files (currently) in kg per gridbox per hour!
  !    What results from this routine now is ug/g, i.e. we scale for time,
  !    gridbox size and air density AND apply a factor of 1e6.
  !
  !    Note that svp is tracer tendency in ug g-1 s-1
  !
  ! 2. R. Janssen 2023/06/29
  !    Convention applied to read emissions in kg hour-1 (per grid cell).
  !    In this routine, we convert to mixing ratios (i.e. ppm or ppb),
  !    because this is the unit that the chemistry scheme needs.
  !
  !    Emitted tracers now align properly with "non-emitted" tracers, e.g.
  !    cloud scalars and secondary chemical components
  
  ! 3. A. Doyennel 2025/05
  !    Separate model input for point sources with/without explicit simulation of plume rise
  !    (Now: supports multiple chemical tracers, having different point sources number )
  !
  ! ----------------------------------------------------------------------

    use modfields,   only : svm, svp
    use modglobal,   only : i1, j1, nsv, &
                            rdt, rtimee, rk3step, &
                            dzf, dx, dy
    use modfields,   only : rhof
    use moddatetime, only : datex, nextday
    use modlsm,      only : lags
    use modmpi,      only : myid

    implicit none

    character(len=*), parameter :: routine = modname//'/emission'

    integer         :: i, j, k, l, iem

    real            :: emistime_s, emistime_e ! Emission timers
    real, parameter :: div3600 = 1./3600.     ! Quick division
    real            :: tend
    real            :: conv_factor, factor, sf
    real, parameter :: MW_air = 28.97

    if (.not. (l_emission)) return

    ! --------------------------------------------------------------------------
    ! Interpolate and apply emission
    ! --------------------------------------------------------------------------
    emistime_s = mod(rtimee +       1800., 3600.)*div3600

    ! MdB NOTE : Better way to do this? Problem is the broadcasting of 1D arrays
    ! rhof and dzf to emisfield. For now, loop over k.
    ! BvS NOTE: I wrote out the loop, to prevent needing a temporary 2D field to store `tend`.

    do k = 1, kemis
      do i = 2, i1
        do j = 2, j1
          ! do l = svskip+1, nsv
          iem = 1
          do l = 1, nsv
            if (.not. tracer_prop(l)%lemis)  cycle
            tend = ((1. - emistime_s)*emisfield(i,j,k,iem,1) + &
                          emistime_s *emisfield(i,j,k,iem,2))

            ! old unit conversion: from kg/hour to ug/g
            ! conv_factor = 1/(3600.*rhof(k)*dzf(k)*dx*dy*1e-6)

            ! new unit conversion: from kg hour-1 to ppb or ppm
            if ( trim(tracer_prop(l)%unit) == 'ppb' ) then
              factor = 1.e9
            elseif ( trim(tracer_prop(l)%unit) == 'ppm' ) then
              factor = 1.e6
            elseif ( trim(tracer_prop(l)%unit) == 'kg m-3') then
              factor = 1.0 ! direct mass concentration, no scale factor here 
            else
              call finish(routine, 'factor not defined for this unit: ', trim(tracer_prop(l)%unit))
            endif

            if (tracer_prop(l)%molar_mass < 0.0 .and. .not. ( trim(tracer_prop(l)%unit) == 'kg m-3')) then
              call finish(routine, 'molar mass not defined for this tracer: ', trim(tracer_prop(l)%tracname))
            endif
            
            if  ( trim(tracer_prop(l)%unit) == 'kg m-3') then
              conv_factor = 1.0 / (dzf(k)*dx*dy) * div3600  ! [kg m-3 s-1]
            else
              conv_factor = 1/(rhof(k)*dzf(k)*dx*dy) * div3600 * MW_air/tracer_prop(l)%molar_mass * factor
            endif
            
            if (l_scale) then
              sf = scalefactor(iem)
            else
              sf = 1.
            endif
            
            if (lags .and. tracer_prop(l)%tracname == 'co2') then
              ! Add tendency to CO2 sum field
              if (trim(tracer_prop(l)%tracname) == 'co2sum') then
                  svp(i,j,k,svco2sum) = svp(i,j,k,svco2sum) + tend * conv_factor * sf
              end if 
            endif

            ! Add tendency to tracer field
            svp(i,j,k,tracer_prop(l)%trac_idx) = svp(i,j,k,tracer_prop(l)%trac_idx) + tend * conv_factor * sf
            !if (i==10 .and. j==10 .and. k==1) then
            ! write(6,"(A18, I2, A7)") "applying species: ", tracer_prop(l)%trac_idx, trim(tracer_prop(l)%tracname)
            ! write(*,*) 'indices   ', i,j,k,tracer_prop(l)%trac_idx
            ! write(*,*) 'emisfield ', emisfield(i,j,k,iem,1)
            ! write(*,*) 'tend      ', tend
            ! write(*,*) 'svp       ', svp(i,j,k,tracer_prop(l)%trac_idx)
            !endif
            
            iem = iem + 1
          end do
        end do
      end do
    end do

    ! -----
    ! Point sources
    ! Intra-hour interpolation is applied
    iem = 1
    do l = 1, nsv
        if (.not. tracer_prop(l)%lemis) cycle
    
        ! Check if current tracer has point sources
        if (l_points .and. (point_sources(l)%npoints > 0)) then
            call applypoints(l, iem)  ! Pass the tracer index to applypoints
        end if
        iem=iem+1
    end do

    ! --------------------------------------------------------------------------
    ! Read emission files when neccesary, i.e. simulation reaches half hour mark
    ! after current timestep
    ! --------------------------------------------------------------------------

    if ( rk3step == 3 ) then
        emistime_e = mod(rtimee + rdt + 1800., 3600.)*div3600

        if ( emistime_e < emistime_s ) then
            ! Transfer data from 'ahead-of-modeltime' field to 'past-modeltime' field
            emisfield(:,:,:,:,1) = emisfield(:,:,:,:,2)

            ! Read new 'ahead-of-modeltime' emission field
            if ( datex(4) == 23 ) then
                call reademission(nextday(1), nextday(2), nextday(3),          0, emisfield(:,:,:,:,2))
            else
                call reademission(  datex(1),   datex(2),   datex(3), datex(4)+1, emisfield(:,:,:,:,2))
            endif
        endif


        ! --------------------------------------------------------------------------
        ! Point sources, with inter-hour interpolation
        ! --------------------------------------------------------------------------
        ! Check if you need to handle point sources

        if (l_points .and. emistime_e < emistime_s) then

            ! Ensure point_sources is allocated on all ranks
          if (.not. allocated(point_sources)) then
              allocate(point_sources(nsv))
          endif
          
          do l = 1, nsv
            if (.not. tracer_prop(l)%lemis)  cycle
            
            if (myid == 0) then
                if (datex(4) == 23) then
                    if (point_sources(l)%npoints == 0) then
                        call inquirepoints(nextday(1), nextday(2), nextday(3), 0)
                    endif
                
                    call readpoints(nextday(1), nextday(2), nextday(3), 0, 1)
                else
                    if (point_sources(l)%npoints == 0) then
                        call inquirepoints(datex(1), datex(2), datex(3), datex(4)+1)
                    else
                        point_sources(l)%data(:,:,1) = point_sources(l)%data(:,:,2)
                    endif
                
                    call readpoints(datex(1), datex(2), datex(3), datex(4)+1, 2)
                endif
            endif
            
            ! Distribute point source data to all ranks
            call distributepoints(l)
          
          end do
      endif

    endif

  end subroutine emission

  ! --------------------------------------------------------------------------
  ! Cleanup after run.
  ! --------------------------------------------------------------------------
  subroutine exitemission
    implicit none
    integer         :: l
    
    if (.not. (l_emission)) return
    deallocate(emisfield)
    !deallocate(co2fields)

    if (l_points) then
        do l = 1, size(point_sources)
            if (allocated(point_sources(l)%data)) deallocate(point_sources(l)%data)
        end do
        if (allocated(point_sources)) deallocate(point_sources)
    endif

  end subroutine exitemission

  subroutine inquirepoints(iyear, imonth, iday, ihour)
    !A. Doyennel 2025/05
    !Check how many point sources per tracer in the simulation domain 
    
    use netcdf
    use modmpi,    only: myidx, myidy
    use modglobal, only: nsv

    implicit none

    character(len=*), parameter :: routine = modname//'/inquirepoints'

    integer, intent(in) :: iyear, imonth, iday, ihour
    integer :: ncid, ndimid, np, l
    logical :: points_exist
    character(512) :: fullpath
    character(16)  :: tracname
    character(256) :: filename
    
    do l = 1, nsv
        if (.not. tracer_prop(l)%lemis)  cycle
        tracname = trim(tracer_prop(l)%tracname)

        filename = 'pointsources.____________.' // trim(tracname) // '.nc'
        write(filename(14:17), '(i4.4)') iyear
        write(filename(18:19), '(i2.2)') imonth
        write(filename(20:21), '(i2.2)') iday
        write(filename(22:23), '(i2.2)') ihour
        write(filename(24:25), '(i2.2)') 0  ! minutes

        fullpath = 'emissions/' // trim(filename)

        inquire(file=trim(fullpath), exist=points_exist)

        if (points_exist) then
            call check(nf90_open(trim(fullpath), NF90_NOWRITE, ncid))
            call check(nf90_inq_dimid(ncid, "n", ndimid))
            call check(nf90_inquire_dimension(ncid, ndimid, len=np))
            point_sources(l)%npoints = np
            call check(nf90_close(ncid))
            call message(routine, 'Filename: ', trim(filename), ' Tracer:', trim(tracname), ' npoints=', np)
        else
            point_sources(l)%npoints = 0
            call message(routine, 'Filename: ', trim(filename), ' Tracer: ', trim(tracname), ' has no point sources.')
        end if
    end do

  contains

    subroutine check(status)
        integer, intent(in) :: status
        if (status /= nf90_noerr) then
            call finish(routine, 'NetCDF error in inquirepoints: ', trim(nf90_strerror(status)))
        end if
    end subroutine check

  end subroutine inquirepoints

  subroutine readpoints(iyear, imonth, iday, ihour, itime)
    
    !A. Doyennel 2025/05
    !Read point sources per tracer in the simulation domain 
    
    use netcdf
    use modmpi,    only: myidx, myidy
    use modglobal, only: nsv
    implicit none

    character(len=*), parameter :: routine = modname//'/readpoints'

    integer, intent(in) :: iyear, imonth, iday, ihour, itime

    integer :: ncid, varid, l, np
    character(256) :: filename, fullpath
    character(16)  :: tracname

    do l = 1, nsv
    
        if (.not. tracer_prop(l)%lemis)  cycle
        np = point_sources(l)%npoints
        if (np == 0) cycle

        tracname = trim(tracer_prop(l)%tracname)

        ! Construct filename
        filename = 'pointsources.____________.' // trim(tracname) // '.nc'
        write(filename(14:17), '(i4.4)') iyear
        write(filename(18:19), '(i2.2)') imonth
        write(filename(20:21), '(i2.2)') iday
        write(filename(22:23), '(i2.2)') ihour
        write(filename(24:25), '(i2.2)') 0  ! minutes
        fullpath = 'emissions/' // trim(filename)

        ! Allocate storage
        if (.not. allocated(point_sources(l)%data)) then
            allocate(point_sources(l)%data(np,7,2))
        endif

        ! Read NetCDF
        call check(nf90_open(trim(fullpath), NF90_NOWRITE, ncid))

        call check(nf90_inq_varid(ncid, "x_idx", varid)) !Global domain indexes!
        call check(nf90_get_var(ncid, varid, point_sources(l)%data(:,1,itime)))
        call check(nf90_inq_varid(ncid, "y_idx", varid)) !Global domain indexes!
        call check(nf90_get_var(ncid, varid, point_sources(l)%data(:,2,itime)))
        call check(nf90_inq_varid(ncid, "height", varid))
        call check(nf90_get_var(ncid, varid, point_sources(l)%data(:,3,itime)))
        call check(nf90_inq_varid(ncid, "temperature", varid))
        call check(nf90_get_var(ncid, varid, point_sources(l)%data(:,4,itime)))
        call check(nf90_inq_varid(ncid, "volume", varid))
        call check(nf90_get_var(ncid, varid, point_sources(l)%data(:,5,itime)))
        call check(nf90_inq_varid(ncid, "stack_exit_area", varid))
        call check(nf90_get_var(ncid, varid, point_sources(l)%data(:,6,itime)))
        call check(nf90_inq_varid(ncid, "emission", varid))
        call check(nf90_get_var(ncid, varid, point_sources(l)%data(:,7,itime)))
        
        call check(nf90_close(ncid))

        call message(routine, 'Read ', np, ' point sources', ' from ', trim(filename), ' for tracer: ', trim(tracname))
    end do

  contains

    subroutine check(status)
        integer, intent(in) :: status
        if (status /= nf90_noerr) then
            call finish(routine, 'NetCDF error in readpoints.', trim(nf90_strerror(status)))
        end if
    end subroutine check

  end subroutine readpoints
  
  subroutine distributepoints(l)
  
    !A. Doyennel 2025/05
    !Distribute point sources for each rank
    
    use modmpi, only: myid, comm3d, D_MPI_BCAST_INT32_R1, D_MPI_BCAST_REAL64_R3
    implicit none

    integer, intent(in) :: l
    integer :: npoints, ierr
    integer :: npoints_array(1)  ! Temporary array to hold npoints

    !------------------------------------------------------------
    ! Broadcast number of point sources for tracer `l` from rank 0
    !------------------------------------------------------------
    if (myid == 0) then
        npoints = point_sources(l)%npoints
    endif

    ! Broadcast integer npoints using the specific subroutine
    npoints_array(1) = npoints
    call D_MPI_BCAST_INT32_R1(npoints_array, 1, 0, comm3d, ierr)  ! Broadcast npoints in an array
    npoints = npoints_array(1)  ! Retrieve npoints from the array

    point_sources(l)%npoints = npoints

    !------------------------------------------------------------
    ! Allocate point_sources(l)%data if not already done (on non-root)
    !------------------------------------------------------------
    if (myid /= 0 .and. npoints > 0) then
        if (.not. allocated(point_sources(l)%data)) then
            allocate(point_sources(l)%data(npoints, 7, 2))
        endif
    endif

    !------------------------------------------------------------
    ! Broadcast the actual point source data (7 vars, 2 times)
    !------------------------------------------------------------
    if (npoints > 0) then
        ! Broadcast the actual point source data using the specific subroutine for 3D REAL64
        call D_MPI_BCAST_REAL64_R3(point_sources(l)%data, npoints*7*2, 0, comm3d, ierr)
    endif
    
  end subroutine distributepoints

  subroutine applypoints(l, iem)
  
    ! A. Doyennel 2025/05
    ! ----------------------------------------------------------------------  !
    ! Purpose:
    !   Applies point source emission tendencies to svp.
    !
    !   Handles plume injection based on two approaches:
    !
    !   a) Explicit plume rise simulation (flag explicit_plume_rise in namoptions):
    !      - Alters potential temperature and momentum tendencies directly (by heat and vertical velocity from the plume)
    !        emission applied around stack heights index
    !      - Suitable for high-resolution LES, where plume dynamics are resolved (<=50m horizontal resolutions recommended)
    !
    !   b) Parameterized modeling of plume rise:
    !      - Uses Briggs’ empirical model formulation
    !      - Useful for coarse LES grid setups where plume rise cannot be properly resolved explicitly 
    !
    !  !   Notes:
    !     - Approach (b) (Briggs) does not modify any model atmospheric thermodynamic fields.
    !       It only estimates plume rise height and distributes emissions accordingly in svp.
    !     - In contrast, approach (a) directly modifies potential temperature and vertical momentum mechanically
    !       to simulate buoyant plume dynamics.
    !   
    !     - Supports application of point sources to multiple scalar tracers
    !       (accessed by outer loop through l)
    ! ----------------------------------------------------------------------
    
    use modfields, only: svp, u0, v0, tmp0, rhof
    use modglobal, only : kmax, dx, dy, dzf,rdt, rtimee, nsv, zh, zf, imax, jmax
    use modmpi, only : myidx, myidy

    implicit none

    character(len=*), parameter :: routine = modname//'/applypoints'

    integer, intent(in) :: l, iem
    integer :: ipoint, ix, iy, iz, isv, izt, izb, iheight, i, j, k
    real    :: emis_b,emis_a, emis_top, emis_bot, emis_in_between
    real    :: plume_top_fraction, plume_bottom_fraction, plumefactor
    real    ::  hmax, ztop, zbottom

    real            :: emistime_s, emistime_e ! Emission timers
    real, parameter :: div3600 = 1./3600.     ! Quick division
    real            :: tend
    real            :: factor, sf
    real, parameter :: MW_air = 28.97
    
    integer :: ixg, iyg         ! Global indices
    integer :: istart, jstart
    logical :: point_is_local
    
    ! For Gaussian vertical distribution (implicit)
    logical :: use_gaussian
    real    :: dz_plume, mean_dz, plume_sigma
    real    :: z_layer_center, z_plume_center
    real    :: weight, sum_weight, emis_k, denom
    real, dimension(kmax) :: plume_shape
    integer :: nlevels
    
    ! Local variables for Gaussian emission injection (explicit)
    real :: sigma_h, sigma_z       ! horizontal and vertical Gaussian spread [m]
    real :: z_center, z_layer      ! plume center height and current layer height [m]
    integer :: k_center, klow, khigh
    real :: weight_xy, weight_z    ! Gaussian weights in x-y and z
    real :: total_weight           ! sum of all weights for normalization
    real :: r2                     ! squared horizontal distance [m^2]
    real :: dv                     ! fraction of emission assigned to a cell
    real, parameter :: MIN_SIGMA_H = 5.0
    real, parameter :: MIN_SIGMA_Z = 2.0
    integer :: i_min,i_max,j_min,j_max
    
    ! Calculate start of this rank’s local domain (DALES global indices start at 2!)
    istart = myidx * imax + 2
    jstart = myidy * jmax + 2                
     
    ! Loop over each point source for the current tracer
    do ipoint = 1, point_sources(l)%npoints
                
            ! Read global grid indices (+2) from data array
            ixg = int(point_sources(l)%data(ipoint, 1, 1)+0.1)  !  full grid (+2) x-index 
            iyg = int(point_sources(l)%data(ipoint, 2, 1)+0.1)  !  full grid (+2) y-index
            
            ! Check if this point falls within local subdomain
            point_is_local = ixg >= istart .and. ixg < istart + imax .and. &
                 iyg >= jstart .and. iyg < jstart + jmax

            if (point_is_local) then
                ! Convert to local indices
                ix = ixg - istart + 2  
                iy = iyg - jstart + 2
                 
                ! Emission values for past and ahead model time
                emis_b = point_sources(l)%data(ipoint, 7, 1)  ! Emission for 'past-modeltime'
                emis_a = point_sources(l)%data(ipoint, 7, 2)  ! Emission for 'ahead-of-modeltime'
                                
                ! --------------------------------------------------------------------------
                ! Interpolate emission (now, the temporal interpolation is in the same way as for area emissions)
                ! --------------------------------------------------------------------------
                emistime_s = mod(rtimee +       1800., 3600.)*div3600
                tend=((1. - emistime_s)*emis_b + emistime_s *emis_a)

                ! old unit conversion: from kg/hour to ug/g
                ! conv_factor = 1/(3600.*rhof(k)*dzf(k)*dx*dy*1e-6)

                ! new unit conversion: from kg hour-1 to ppb or ppm

                if ( trim(tracer_prop(l)%unit) == 'ppb' ) then
                    factor = 1.e9
                elseif ( trim(tracer_prop(l)%unit) == 'ppm' ) then
                    factor = 1.e6
                elseif ( trim(tracer_prop(l)%unit) == 'kg m-3') then
                    factor = 1.0 ! direct mass concentration, no scale factor here 
                else
                    call finish(routine, 'factor not defined for this unit', trim(tracer_prop(l)%unit))
                endif
                
                if (tracer_prop(l)%molar_mass < 0.0 .and. .not. ( trim(tracer_prop(l)%unit) == 'kg m-3')) then
                    call finish(routine, 'molar mass not defined for this tracer', trim(tracer_prop(l)%tracname))
                endif

                sf = merge(scalefactor(iem), 1.0, l_scale) !Analog of if-else statement
                
                if (explicit_plume_rise) then !flag for explicit LES simulation of emission plume rise:
                
                    ! ----------------------------------------------------------------------------
                    ! Inject heat into the potential temperature tendency field (thlp)
                    ! ----------------------------------------------------------------------------
                    call inject_heat_source(ix, iy, &
                            point_sources(l)%data(ipoint, 3, 1), &  ! Stack height [m]
                            point_sources(l)%data(ipoint, 4, 1), &  ! Exhaust temp Ts [K]
                            point_sources(l)%data(ipoint, 5, 1), &  ! Volumetric flow rate Vs [m³/s]
                            use_gaussian = .true.)

                    ! ----------------------------------------------------------------------------
                    ! Inject momentum into the vertical velocity tendency field (wp)
                    ! ----------------------------------------------------------------------------                    
                    call inject_momentum_source(ix, iy, &
                            point_sources(l)%data(ipoint, 3, 1), &  ! Stack height [m]
                            point_sources(l)%data(ipoint, 4, 1), &  ! Exhaust temp Ts [K]
                            point_sources(l)%data(ipoint, 5, 1), &  ! Volumetric flow rate Vs [m³/s]
                            point_sources(l)%data(ipoint, 6, 1), &  ! Stack exit area [m²]
                            use_gaussian = .true.)  

                    ! ----------------------------------------------------------------------------
                    ! Inject emission tendencies into tracer source (svp) (Gaussian spread)
                    ! ----------------------------------------------------------------------------

                    ! === Plume center and grid indices ===
                    k_center = minloc(abs(zf - point_sources(l)%data(ipoint,3,1)), dim=1)
                    z_center = zf(k_center)

                    ! Horizontal and vertical spread (PALM-like)
                    sigma_h = max(1.5*dx, MIN_SIGMA_H)
                    sigma_z = max(1.5*dzf(k_center), MIN_SIGMA_Z)

                    ! Plume grid (~3x3x3 cells) (tune me, now the 1x1x1 is used!)
                    !i_min = max(1, ix-1); i_max = min(size(svp,1), ix+1)
                    !j_min = max(1, iy-1); j_max = min(size(svp,2), iy+1)
                    !izb  = max(1, k_center-1); izt = min(size(zf), k_center+1)
                    
                    ! Plume grid indices (~1x1x1)
                    i_min = ix; i_max = ix
                    j_min = iy; j_max = iy
                    izb  = k_center; izt = k_center

                    ! === Compute total Gaussian weight over plume grid ===
                    total_weight = 0.0
                    do k = izb, izt
                        z_layer = zf(k)
                        weight_z = exp(-((z_layer - z_center)**2)/(2.0*sigma_z**2))
                        do i = i_min,i_max
                            do j = j_min,j_max
                                r2 = ((real(i-ix)*dx)**2 + (real(j-iy)*dy)**2)
                                weight_xy = exp(-r2/(2.0*sigma_h**2))
                                total_weight = total_weight + weight_z * weight_xy
                            end do
                        end do
                    end do

                    ! === Apply emission tendency (distributed over plume grid) ===
                    do k = izb, izt
                        z_layer = zf(k)
                        weight_z = exp(-((z_layer - z_center)**2)/(2.0*sigma_z**2))
                        do i = i_min,i_max
                            do j = j_min,j_max
                                r2 = ((real(i-ix)*dx)**2 + (real(j-iy)*dy)**2)
                                weight_xy = exp(-r2/(2.0*sigma_h**2))
                                ! Normalized 3D Gaussian
                                dv = (weight_z * weight_xy) / total_weight
                                svp(i,j,k,tracer_prop(l)%trac_idx) = svp(i,j,k,tracer_prop(l)%trac_idx) &
                                                 + compute_tendency_func(tend * dv, k, factor, l, sf)
                            end do
                        end do
                    end do
              
                else
                
                    ! ===Briggs empirical model (for coarse resolution simulations): Compute vertical range and inject into izb to izt
                    ! Call briggs subroutine to calculate plume parameters
                    call briggs(tmp0(ix, iy, 1:kmax), &                             ! Temperature profile
                        sqrt(v0(ix, iy, 1:kmax)**2 + u0(ix, iy, 1:kmax)**2), & ! Total horizontal windspeed profile

                        point_sources(l)%data(ipoint, 4, 1), &    ! Source temperature
                        point_sources(l)%data(ipoint, 5, 1), &    ! Source volumetric flow rate
                        point_sources(l)%data(ipoint, 3, 1), &    ! Source stack height

                        izt, plume_top_fraction, &               ! Full level index for plume top
                        izb, plume_bottom_fraction, hmax, ztop, zbottom)  

                    !-------------------------------------------------------------------------------------------------

                    ! Emissions are per source, so refactor to emission per gridbox
                    ! ALSO: Emissions are per source, per hour so refactor to account for pressure, gridboxsize and seconds below:

                    if (izt - izb > 1) then

                        denom = (izt - izb - 1) + plume_bottom_fraction + plume_top_fraction

                        if (denom > 0.0) then
                            emis_bot = tend * plume_bottom_fraction / denom
                            emis_top = tend * plume_top_fraction    / denom
                            emis_in_between = tend / denom
                        else
                            ! Fallback: should not happen, but keep safe
                            emis_bot = 0.0
                            emis_top = 0.0
                            emis_in_between = 0.0
                        end if

                        if ((plume_top_fraction>1) .OR. (plume_bottom_fraction>1) .OR. (plume_top_fraction<0) .OR. (plume_bottom_fraction<0)) then
                            print*,'plume_top_fraction, plume_bottom_fraction', plume_top_fraction, plume_bottom_fraction
                        endif

                    end if

                    ! Apply interpolated-in-time point source emissions:
                    !-------------------------------------------------------------------------------------------------

                    if (izb == izt) then

                        svp(ix, iy, izb, tracer_prop(l)%trac_idx) = svp(ix, iy, izb, tracer_prop(l)%trac_idx) + compute_tendency_func(tend, izb,factor,l, sf)
                    
                    else if (izt - izb == 1) then

                        svp(ix, iy, izb, tracer_prop(l)%trac_idx) = svp(ix, iy, izb, tracer_prop(l)%trac_idx) + compute_tendency_func(tend/2, izb,factor,l, sf)
                        svp(ix, iy, izt, tracer_prop(l)%trac_idx) = svp(ix, iy, izt, tracer_prop(l)%trac_idx) + compute_tendency_func(tend/2, izt,factor,l, sf)

                    else if (izt - izb > 1) then
                        ! Compute mean dz in affected layers
                        dz_plume = 0.0
                        nlevels = izt - izb + 1

                        do k = izb, izt
                            dz_plume = dz_plume + dzf(k)
                        end do

                        mean_dz = dz_plume / real(nlevels)

                        ! Determine whether to use Gaussian
                        use_gaussian = (hmax > 20.0) .and. ((hmax / mean_dz) >= 3.0)

                        if (use_gaussian) then
                            ! Gaussian vertical distribution
                            do k = izb, izt
                                ! Compute layer center height relative to plume center
                                z_layer_center = zf(k)
                                z_plume_center = 0.5 * (ztop + zbottom)
                                plume_sigma = max((ztop - zbottom) / 4.0, 0.5 * mean_dz)  ! Stddev = 1/4 of plume height range

                                ! Gaussian weight (unnormalized)
                                weight = exp(- ((z_layer_center - z_plume_center)**2) / (2.0 * plume_sigma**2))
                                plume_shape(k) = weight
                            end do

                            ! Normalize weights so total = 1
                            sum_weight = 0.0
                            do k = izb, izt
                                sum_weight = sum_weight + plume_shape(k)
                            end do

                            do k = izb, izt
                                emis_k = tend * (plume_shape(k) / sum_weight)
                                svp(ix, iy, k, tracer_prop(l)%trac_idx) = svp(ix, iy, k, tracer_prop(l)%trac_idx) + &
                                compute_tendency_func(emis_k, k, factor, l, sf)
                            end do

                        else
                            ! Uniform linear interpolation
                            svp(ix, iy, izb, tracer_prop(l)%trac_idx) = svp(ix, iy, izb, tracer_prop(l)%trac_idx) + &
                                compute_tendency_func(emis_bot, izb, factor, l, sf)
                            svp(ix, iy, izt, tracer_prop(l)%trac_idx) = svp(ix, iy, izt, tracer_prop(l)%trac_idx) + &
                                compute_tendency_func(emis_top, izt, factor, l, sf)

                            do k = izb+1, izt-1
                                svp(ix, iy, k, tracer_prop(l)%trac_idx) = svp(ix, iy, k, tracer_prop(l)%trac_idx) + &
                                compute_tendency_func(emis_in_between, k, factor, l, sf)
                            end do
                        end if
                    end if
            end if
        end if
    enddo

    contains

    function compute_tendency_func(etend, k, factor, l, sf) result(tendency)
        integer, intent(in) :: k, l
        real, intent(in) :: etend, factor, sf
        real :: tendency
        real :: tmp
        logical :: is_valid
        real, parameter :: div3600 = 1.0 / 3600.0
        real, parameter :: MW_air = 28.97
        
        if ( trim(tracer_prop(l)%unit) == 'kg m-3') then
           tmp = etend * ((1.0 / (dzf(k)*dx*dy)) * div3600) * sf
        else 
           tmp = etend * ((1/(rhof(k)*dzf(k)*dx*dy)) * div3600 * MW_air / tracer_prop(l)%molar_mass * factor) * sf
        endif
        
        is_valid = (tmp == tmp .and. abs(tmp) <= 1.0e6)

        if (.not. is_valid) then
            print*,'Warning: point source emission tendency value is invalid (set to 0.0). Raw value = ', tmp, k
        end if

        tendency = merge(tmp, 0.0, is_valid)

    end function compute_tendency_func

  end subroutine applypoints

  subroutine briggs(Ta, U, Ts, Vs, hs, iztop, ztop_frac, izbottom, zbottom_frac, hmax, ztop, zbottom)

    !Briggs algorithm to calculate the vertical plume rise above the stack height
    !The detail description can be found in Gordon et al., (2017) and Akingunola et al., (2018)

    use modglobal,   only : zh, zf, dzf, kmax, pi, cp, grav
    !----------
    ! Ta Atmospheric temperature, K
    ! U  Total horizontal wind speed sqrt(v0ˆ2 + u0ˆ2), m/s
    ! Ts Emission temperature K
    ! Vs Emission volumetric flow rate m³/s
    ! hs Emission stack height, m
    ! tzh Atmospheric temperature at half-level grid
    ! uzh Wind speed at half-level grid
    ! ths Atmospheric temperature at stack height
    ! uhs Wind speed at stack height
    ! hmax plume rise height (calculated relative to hs)
    ! S the stability parameter
    ! Fb buoyancy flux at stack height
    ! F1 Residual buoyancy flux ( iz +1)
    ! F0 Residual buoyancy flux ( iz)
    ! F0_old Residual buoyancy flux ( iz-1)
    ! iz half-grid index starts from the top of "stack" layer
    ! ---------

    implicit none

    real, intent(in)  :: Ts, hs, vs
    real(field_r), dimension(kmax), intent(in) :: Ta, U
    integer, intent(out) :: iztop, izbottom
    real,    intent(out) :: ztop_frac, zbottom_frac, hmax, ztop, zbottom

    integer :: iz, i, kbelow, kabove, ieq, iz0
    real    :: F0, F0_old, F1, Fb, dT, dU, ths, uhs
    real             :: gradT, S, S_eff, upperh, lowerh, lowerw
    real, dimension(kmax+1) :: tzh, uzh
    real, parameter :: min_plume_thickness = 10.0 !tune it if needed
    
    !============================================================
    ! 1. Compute half-level (interface) values for T and U
    !============================================================

    do i = 1, kmax-1
        tzh(i+1) = 0.5 * (Ta(i) + Ta(i+1))
        uzh(i+1) = 0.5 * (U(i) + U(i+1))
    end do
    
    tzh(1) = Ta(1)
    uzh(1) = U(1)
    
    !============================================================
    ! 2. Interpolate T and U at stack height (between interfaces)
    !============================================================
    kbelow = maxloc(zh, dim=1, mask=zh <= hs)
    kabove = minloc(zh, dim=1, mask=zh >  hs)

    if (kbelow < 1) kbelow = 1
    if (kabove < kbelow) kabove = kbelow + 1

    ths = tzh(kbelow) + (tzh(kabove) - tzh(kbelow)) * &
          (hs - zh(kbelow)) / max(zh(kabove) - zh(kbelow), 1.0e-6)

    uhs = uzh(kbelow) + (uzh(kabove) - uzh(kbelow)) * &
          (hs - zh(kbelow)) / max(zh(kabove) - zh(kbelow), 1.0e-6)

    !============================================================
    ! 3. Compute initial buoyancy flux at stack height
    !============================================================
    Fb = 0.0
    F1 = 0.0
    
    if (Ts > ths) Fb = (grav / pi) * Vs * (Ts - ths) / Ts
    
    F1     = Fb
    F0     = Fb
    F0_old = Fb
    hmax   = 0.0

    ! Layer integration
    iz = kabove
    iz0 = iz
    S = 0.0

    do while ((F1 > 0.0) .and. (iz <= kmax))

        if (iz == iz0) then
            gradT = (tzh(iz) - ths) / (zh(iz) - hs)
            lowerh = 0.0
            upperh = zh(iz) - hs
            lowerw = 0.5 * (uhs + uzh(iz))
        else
            gradT = (tzh(iz) - tzh(iz - 1)) / (zh(iz) - zh(iz - 1))
            lowerh = zh(iz - 1) - hs
            upperh = zh(iz) - hs
            lowerw = 0.5 * (uzh(iz) + uzh(iz - 1))
        end if

        S = grav / tzh(iz) * (gradT + grav / cp)
        S_eff = max(S, 0.0)

        F0_old = F0
        F0 = F1

        if (S_eff > 0.0) then
            F1 = min( &
             F0 - 0.015 * S_eff * max(F0_old,1.0e-6)**(1.0/3.0) * &
                   (upperh**(8.0/3.0) - lowerh**(8.0/3.0)), &
             F0 - 0.053 * S_eff * lowerw * (upperh**3 - lowerh**3) )
        else
             F1 = F0
        end if

        if (F1 <= 0.0) then
            if (abs(S_eff) > 1.0e-6 .and. F0_old > 0.0 .and. lowerw > 0.0) then
                hmax = min( &
                    ((F0 / (0.015 * S_eff * max(F0_old,1e-6)**(1.0/3.0)))**(3.0/8.0) + lowerh), &
                    ((F0 / (0.053 * S_eff * lowerw))**(1.0/3.0) + lowerh) )
                        
            else
                hmax = max(upperh, 0.0)
            end if

            exit
        end if

        iz = iz + 1
    end do

    ! Final plume geometry
    zbottom = hs
    ztop    = hs + max(hmax, min_plume_thickness)

    if (ztop < zbottom + min_plume_thickness) ztop = zbottom + min_plume_thickness

    ! Map plume bounds to DALES grid
    izbottom = minloc(zh, dim=1, mask=zh >= zbottom) - 1
    iztop    = minloc(zh, dim=1, mask=zh >= ztop)    - 1
    izbottom = max(izbottom, 1)
    iztop    = max(iztop, izbottom + 1)

    zbottom_frac = (zh(izbottom+1) - zbottom) / dzf(izbottom)
    ztop_frac    = (ztop - zh(iztop)) / dzf(iztop)
    
    zbottom_frac = max(0.0, min(1.0, zbottom_frac))
    ztop_frac    = max(0.0, min(1.0, ztop_frac))

  end subroutine
  
  subroutine inject_heat_source(ix, iy, hs, Ts, Vs, use_gaussian)
    use modglobal, only: rdt, dzf, zf, dx, dy, cp
    use modfields, only: tmp0, rhof, thlp, exnf
    implicit none

    integer, intent(in) :: ix, iy
    real,    intent(in) :: hs       ! Stack height [m]
    real,    intent(in) :: Ts       ! Stack exit temperature [K]
    real,    intent(in) :: Vs       ! Volumetric flow rate [m³/s]
    logical, intent(in), optional :: use_gaussian

    ! === Physical constants and limits ===
    real, parameter :: MAX_DELTA_THETA = 1.5 ! Max allowed heating [K per timestep] [tune me!]
    !TODO: cup buoyancy instead of theta_tend
    real, parameter :: MIN_SIGMA_H = 5.0     ! minimum horizontal spread [m]
    real, parameter :: MIN_SIGMA_Z = 2.0     ! minimum vertical spread [m]

    ! === Locals ===
    integer :: i,j,k, k_center, kmin,kmax, i_min,i_max,j_min,j_max
    real :: Ta, rho_air, emission_power, heat_tend, rho_stack
    real :: sigma_h, sigma_z, z_center, z_layer, weight_xy, weight_z, total_weight
    real :: volume, max_heat_tend, r2

    ! Find vertical index nearest to stack height
    k_center = minloc(abs(zf - hs), dim=1)
    Ta = tmp0(ix, iy, k_center)
    rho_air = max(1.0e-6, rhof(k_center))
    rho_stack = rho_air * Ta / Ts !ideal gas stack density
    
    ! === Input safety checks ===
    if (Ts <= 0.0 .or. Ta <= 0.0 .or. Vs <= 0.0 .or. rho_air <= 0.0) then
        print *, 'WARNING: Bad inputs in inject_heat_source @', ix, iy, &
                 ' Ts:', Ts, ' Ta:', Ta, ' Vs:', Vs, ' rho:', rho_air, ' hs:', hs, ' k_center:', k_center
        return
    end if

    ! === Total heat emission power ===
    ! Units: [kg/m³] * [m³/s] * [J/kg·K] * [K] = [W] = [J/s]
    emission_power = rho_stack * Vs * cp * (Ts - Ta)

    ! Adaptive Gaussian spread based on grid resolution (PALM-like)
    sigma_h = max(1.5*dx, MIN_SIGMA_H)
    sigma_z = max(1.5*dzf(k_center), MIN_SIGMA_Z)
    z_center = zf(k_center)

    ! Determine plume indices (~3x3x3 cells) (tune me, now the 1x1x1 is used!)
    !i_min = max(1, ix - 1); i_max = min(size(thlp,1), ix + 1)
    !j_min = max(1, iy - 1); j_max = min(size(thlp,2), iy + 1)
    !kmin  = max(1, k_center - 1); kmax = min(size(zf), k_center + 1)
    
    ! Plume grid indices (~1x1x1)
    i_min = ix; i_max = ix
    j_min = iy; j_max = iy
    kmin  = k_center; kmax = k_center

    ! Compute total Gaussian weight
    total_weight = 0.0
    do k = kmin,kmax
        z_layer = zf(k)
        weight_z = exp(-((z_layer - z_center)**2)/(2.0*sigma_z**2))
        do i=i_min,i_max
            do j=j_min,j_max
                r2 = ((real(i-ix)*dx)**2 + (real(j-iy)*dy)**2)
                weight_xy = exp(-r2/(2.0*sigma_h**2))
                total_weight = total_weight + weight_z*weight_xy
            end do
        end do
    end do

    ! === Apply distributed heating tendency [K/s] ===
    do k = kmin,kmax
        z_layer = zf(k)
        weight_z = exp(-((z_layer - z_center)**2)/(2.0*sigma_z**2))
        ! === Max allowed temperature tendency [K/s] 
        max_heat_tend = (MAX_DELTA_THETA * exnf(k)) / rdt
        do i=i_min,i_max
            do j=j_min,j_max
                r2 = ((real(i-ix)*dx)**2 + (real(j-iy)*dy)**2)
                weight_xy = exp(-r2/(2.0*sigma_h**2))
                rho_air = max(1.0e-6, rhof(k))   ! Update for each level
                volume = dx*dy*dzf(k)
                ! [K/s] = [W] * [unitless] / ([kg/m³] * [m³] * [J/kg·K])
                heat_tend = emission_power * weight_z * weight_xy / (total_weight*volume*rho_air*cp)
                if (heat_tend > max_heat_tend) then
                  print *, 'Capping plume theta tendency at', ix, iy, k, ':', heat_tend, '→', max_heat_tend
                  heat_tend = max_heat_tend
                end if
                !print *, ' theta_tend value at', ix, iy, k, ':', heat_tend / exnf(k), '[K/s]'
                thlp(i,j,k) = thlp(i,j,k) + (heat_tend / exnf(k)) ! [K/s] exnf Convert temperature tendency to potential temperature tendency
            end do
        end do
    end do

  end subroutine inject_heat_source

  !=====================================================

  subroutine inject_momentum_source(ix, iy, hs, Ts, Vs, As, use_gaussian)
    !Mechanical injection of the vertical wind velosity
    use modglobal, only: rdt, dzf, zf, dx, dy, pi
    use modfields, only: tmp0, rhof, wp ! wp is vertical wind tendency [m/s²]
    implicit none

    integer, intent(in) :: ix, iy
    real,    intent(in) :: hs       ! Stack height [m]
    real,    intent(in) :: Ts       ! Stack exit temperature [K]
    real,    intent(in) :: Vs       ! Volumetric flow rate [m³/s]
    real,    intent(in) :: As       ! Stack exit area [m²]
    logical, intent(in), optional :: use_gaussian

    real, parameter :: MAX_DVZ_PER_STEP = 0.7 ! Maximum change in vertical velocity [m/s] allowed per timestep [tune me!]
    real, parameter :: MIN_SIGMA_H = 5.0
    real, parameter :: MIN_SIGMA_Z = 2.0

    integer :: i,j,k, k_center, kmin,kmax, i_min,i_max,j_min,j_max
    real :: w_exit, rho_air, dvz, sigma_h, sigma_z, Ta, rho_stack !,D,r,A
    real :: z_center, z_layer, weight_xy, weight_z, total_weight, volume, r2
    real :: MAX_W_TEND ! Max dvz injection rate [m/s²]

    !  ===  Max dvz injection rate assumption [m/s²]  ===
    MAX_W_TEND = MAX_DVZ_PER_STEP / rdt

    ! === Geometry assumptions (rough estimate) ===
    !D = hs / 10.0            ! Effective stack diameter [m] [tune me!]
    !r = D / 2.0              ! Stack radius [m]
    !A = pi * r**2            ! Stack exit area [m²]
    ! Note: Stack exit area [m²] now is used from the input
    
    ! === Safety check ===
    if (As <= 0.0 .or. Vs <= 0.0) then
        print *, 'WARNING: Invalid stack geometry or zero flow at ix=', ix, 'iy=', iy
        return
    end if

    ! === Compute exit velocity ===
    w_exit = Vs / As          ! Stack exit velocity [m/s]

    ! === Find vertical index closest to stack height ===
    k_center = minloc(abs(zf - hs), dim=1)
    
    Ta = tmp0(ix, iy, k_center)
    rho_air = max(1.0e-6, rhof(k_center))
    rho_stack = rho_air * Ta / Ts !ideal gas stack density

    ! Adaptive Gaussian spread
    sigma_h = max(2.0*dx, MIN_SIGMA_H)
    sigma_z = max(1.5*dzf(k_center), MIN_SIGMA_Z)
    z_center = zf(k_center)

    ! Plume grid indices (~3x3x3) (tune me, now the 1x1x1 is used!)
    !i_min = max(1, ix-1); i_max = min(size(wp,1), ix+1)
    !j_min = max(1, iy-1); j_max = min(size(wp,2), iy+1)
    !kmin  = max(1, k_center-1); kmax = min(size(zf), k_center+1)
    
    ! Plume grid indices (~1x1x1)
    i_min = ix; i_max = ix
    j_min = iy; j_max = iy
    kmin  = k_center; kmax = k_center
    
    ! First pass: compute weighted sum including volume and density
    total_weight = 0.0
    do k=kmin,kmax
      rho_air = max(1.0e-6, rhof(k))
      do i=i_min,i_max
        do j=j_min,j_max
            volume = dx*dy*dzf(k)
            weight_z = exp(-((zf(k) - z_center)**2)/(2.0*sigma_z**2))
            r2 = ((real(i-ix)*dx)**2 + (real(j-iy)*dy)**2)
            weight_xy = exp(-r2/(2.0*sigma_h**2))
            total_weight = total_weight + weight_z * weight_xy
        end do
      end do
    end do

    ! Apply momentum tendency (mass flux based) purely mechanical jet!
    do k=kmin,kmax
      rho_air = max(1.0e-6, rhof(k)) ! Update for each layer
      volume = dx*dy*dzf(k) ! Grid cell volume [m³]
      do i=i_min,i_max
        do j=j_min,j_max
            weight_z = exp(-((zf(k) - z_center)**2)/(2.0*sigma_z**2))
            r2 = ((real(i-ix)*dx)**2 + (real(j-iy)*dy)**2)
            weight_xy = exp(-r2/(2.0*sigma_h**2)) 
            dvz = (weight_z * weight_xy / total_weight) * (rho_stack * Vs * w_exit) / (rho_air * volume)  ! Vertical velocity tendency [m/s²]
            if (dvz > MAX_W_TEND) then
                  print *, 'Capping plume wind tendency at', ix, iy, k, ':', dvz, '→', MAX_W_TEND
                  dvz = MAX_W_TEND
            end if
            !print *, ' dvz value at', ix, iy, k, ':', dvz, '[m/s²]'
            wp(i,j,k) = wp(i,j,k) + dvz
        end do
      end do
    end do

  end subroutine inject_momentum_source

end module modemission
