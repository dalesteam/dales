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
!  Copyright 1993-2020 Delft University of Technology, Wageningen
!  University, Utrecht University, KNMI
!

module modemission

use modemisdata
use modtracers,       only : tracer_prop

implicit none

! integer :: nemis

contains

  subroutine initemission

    use modglobal,    only : i2, j2,kmax, nsv, ifnamopt, fname_options, checknamelisterror
    use modmpi,       only : myid, comm3d, d_mpi_bcast
    use moddatetime,  only : datex, prevday, nextday

    implicit none

    ! Auxiliary variables
    integer :: ierr

    ! --- Read & broadcast namelist EMISSION -----------------------------------
    ! namelist/NAMEMISSION/ l_emission, kemis, svskip, emisnames, svco2sum
    namelist/NAMEMISSION/ l_emission, kemis, nemis, emisnames, l_scale, scalefactor

    if (myid == 0) then

      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMEMISSION,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMEMISSION')
      write(6, NAMEMISSION)
      close(ifnamopt)

    endif

    call d_mpi_bcast(l_emission,           1,  0, comm3d, ierr)
    call d_mpi_bcast(kemis,                1,  0, comm3d, ierr)
    call d_mpi_bcast(emisnames(1:100),   100,  0, comm3d, ierr)
    call d_mpi_bcast(nemis,                1,  0, comm3d, ierr)
    call d_mpi_bcast(l_scale,              1,  0, comm3d, ierr)
    call d_mpi_bcast(scalefactor(1:100), 100,  0, comm3d, ierr)

    ! -- Interaction with AGs   ----------------------------------------------------
    if (.not. (l_emission)) return
    allocate(co2fields(nsv))

    co2fields = 0
    ! co2fields(svskip+1:nsv) = index(emisnames(1:nsv-svskip), "co2")
    co2fields = index(emisnames, "co2")


    ! svco2ags = findloc(emisnames(1:nsv-svskip), value = "co2ags", dim = 1)
    ! svco2ags = svco2ags + svskip
    svco2ags = findloc(emisnames, value = "co2ags", dim = 1)

    ! svco2veg = findloc(emisnames(1:nsv-svskip), value = "co2veg", dim = 1)
    ! svco2veg = svco2veg + svskip
    svco2veg = findloc(emisnames, value = "co2veg", dim = 1)

    if (myid == 0) then
      write(6,*) 'modemission: co2fields (scalar fields with CO2 0=no, 1=yes)'
      write(6,*) co2fields
      write(6,*) 'modemission: svco2ags (scalar field number for AGS emissions)'
      write(6,*) svco2ags
      write(6,*) 'modemission: svco2veg (scalar field number for AGS emissions)'
      write(6,*) svco2veg
      write(6,*) 'number of emitted species'
      write(6,*) nemis
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

    write(6,"(A18, A12)") "Reading emission: ", sdatetime

    iem = 1
    do isv = 1, nsv
      if (tracer_prop(isv)%lemis) then
        ! check tracer unit
        ! give warning when emission file is not available for a species which is emitted  
        if (iem > nemis) then
          write(6,"(A52, I2, A3, I2)") "More emitted species than declared in NAMEMISSION: ", iem, " > ", nemis
        endif
        write(6,"(A17, I2, A7)") "Reading tracer: ", tracer_prop(isv)%trac_idx, trim(tracer_prop(isv)%tracname)
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
          write(6,"(A38, A36, A14)") "WARNING: emission units do not match: " , trim(unit), " /= kg hour-1"
        endif
        iem = iem + 1
      else 
        write(6,"(A20, I2, A7)") "Tracer not emitted: ", tracer_prop(isv)%trac_idx, trim(tracer_prop(isv)%tracname)  
      endif
    end do

  contains

  subroutine check(status)
    integer, intent(in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop 'NetCDF error in modemission. See outputfile for more information.'
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
  ! ----------------------------------------------------------------------

    use modfields,   only : svp
    use modglobal,   only : i1, j1, nsv, &
                            rdt, rtimee, rk3step, &
                            dzf, dx, dy
    use modfields,   only : rhof
    use moddatetime, only : datex, nextday
    use modlsm,      only : lags

    implicit none

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
            else
              print *, trim(tracer_prop(l)%unit)
              STOP 'factor not defined for this unit'
            endif

            if (tracer_prop(l)%molar_mass < 0. ) then
              print *, trim(tracer_prop(l)%tracname)
              STOP 'molar mass not defined for this tracer'
            endif
            conv_factor = 1/(rhof(k)*dzf(k)*dx*dy) * div3600 * MW_air/tracer_prop(l)%molar_mass * factor

            if (l_scale) then
              sf = scalefactor(iem)
            else
              sf = 1.
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
            if (lags) then
              ! Add tendency to CO2 sum field
              svp(i,j,k,svco2sum) = svp(i,j,k,svco2sum) + tend * conv_factor * sf
            endif
            iem = iem + 1
          end do
        end do
      end do
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

    endif

  end subroutine emission

  ! --------------------------------------------------------------------------
  ! Cleanup after run.
  ! --------------------------------------------------------------------------
  subroutine exitemission

    implicit none

    if (.not. (l_emission)) return

    deallocate(emisfield)
    deallocate(co2fields)

  end subroutine exitemission

end module modemission
