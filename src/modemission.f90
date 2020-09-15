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

implicit none

contains

  subroutine initemission

    use modglobal,    only : i2, j2,kmax, nsv, ifnamopt, fname_options
    use modmpi,       only : myid, comm3d, mpi_logical, mpi_integer, mpi_character
    use moddatetime,  only : datex, prevday, nextday

    implicit none
  
    ! Auxiliary variables
    integer :: ierr

    ! --- Read & broadcast namelist EMISSION -----------------------------------
    namelist/NAMEMISSION/ l_emission, kemis, sv_skip, svlist 

    if (myid == 0) then

      open(ifnamopt, file=fname_options, status='old', iostat=ierr)
      read(ifnamopt, NAMEMISSION, iostat=ierr)

      if (ierr > 0) then
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMEMISSION'
      endif

      write(6, NAMEMISSION)
      close(ifnamopt)

    endif

    call mpi_bcast(l_emission,    1,   mpi_logical,   0, comm3d, ierr)
    call mpi_bcast(kemis,         1,   mpi_integer,   0, comm3d, ierr)
    call mpi_bcast(sv_skip,       1,   mpi_integer,   0, comm3d, ierr)
    call mpi_bcast(svlist(1:100), 100, mpi_character, 0, comm3d, ierr)

    ! --- Local pre-calculations and settings
    if (.not. (l_emission)) return

    if (kemis == -1) kemis = kmax

    ! --- Read emission files for first time step ----------------------------------

    ! Two hourly emission fields are loaded at all times: 
    ! (1) before model time,   t_field < t_model, "in-the-past"
    ! (2) ahead of model time, t_field > t_model, "in-the-future"
    allocate(svemis(i2, j2, kemis, sv_skip+1:nsv, 2))

    if (datex(5) >= 30) then
      call reademission(    datex(1),   datex(2),   datex(3),   datex(4), svemis(:,:,:,:,1))
      
      if (datex(4) == 23) then  
        call reademission(nextday(1), nextday(2), nextday(3),          0, svemis(:,:,:,:,2))
      else
        call reademission(  datex(1),   datex(2),   datex(3), datex(4)+1, svemis(:,:,:,:,2))
      endif

    else
      call reademission(    datex(1),   datex(2),   datex(3),   datex(4), svemis(:,:,:,:,2)) 
     
      if (datex(4) == 0) then
        call reademission(prevday(1), prevday(2), prevday(3),         23, svemis(:,:,:,:,1))
      else
        call reademission(  datex(1),   datex(2),   datex(3), datex(4)-1, svemis(:,:,:,:,1))
      endif

    endif 
  end subroutine initemission

  subroutine reademission(iyear, imonth, iday, ihour, emisfield)

  ! ----------------------------------------------------------------------
  ! Reading of emission files
  ! Multiple/all tracers
  ! ----------------------------------------------------------------------

    use mpi,         only : MPI_INFO_NULL
    use netcdf
    use modmpi,      only : myid, myidx, myidy, comm3d
    use modglobal,   only : i1, j1, i2, j2, imax, jmax, nsv
    use moddatetime, only : datex

    implicit none

    integer, intent(in)  :: iyear, imonth, iday, ihour     
    real, intent(out)    :: emisfield(i2, j2, kemis, 1+sv_skip:nsv)

    integer, parameter   :: ndim = 3 
    integer              :: start(ndim), count(ndim)
    integer              :: ncid, varid
    integer              :: isv
 
    character(len=12)    :: sdatetime

    ! Create string from given date
    write(sdatetime, "(I0.4,2I0.2,2I0.2)") iyear, imonth, iday, ihour, 0 

    write(6,"(A18, A12)") "Reading emission: ", sdatetime

    do isv = 1+sv_skip, nsv

      call check( nf90_open( trim(svlist(isv-sv_skip))//'_emis_'//sdatetime//'_3d.nc', IOR(NF90_NOWRITE, NF90_MPIIO), &
                              ncid, comm = comm3d, info = MPI_INFO_NULL) )
      call check( nf90_inq_varid( ncid, svlist(isv-sv_skip), varid) )
      call check( nf90_get_var  ( ncid, varid, emisfield(2:i1,2:j1,1:kemis,isv), & 
                                  start = (/1 + myidx * imax, 1 + myidy * jmax, 1, 1/), &
                                  count = (/imax, jmax, kemis, 1/) ) )
      call check( nf90_close( ncid ) )

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
  !    gridbx size and air density AND apply a factor of 1e6.
  ! 
  ! TODO
  ! 1. MDB Align properly with "non-emitted" tracers, i.e. cloud scalars from e.g.
  ! microphysics/chemistry
  ! ----------------------------------------------------------------------

    use modfields,   only : svp
    use modglobal,   only : i1, j1, ih, jh, nsv, &
                            rdt, rtimee, rk3step, &
                            dzf, dx, dy    
    use modfields,   only : rhof
    use moddatetime, only : datex, nextday
    
    implicit none

    integer         :: k
    
    real            :: emistime_s, emistime_e ! Emission timers
    real, parameter :: div3600 = 1./3600.     ! Quick division
 
    if (.not. (l_emission)) return
 
    ! --------------------------------------------------------------------------
    ! Interpolate and apply emission
    ! --------------------------------------------------------------------------
    emistime_s = mod(rtimee +       1800., 3600.)*div3600

    ! MdB NOTE : Better way to do this? Problem is the broadcasting of 1D arrays
    ! rhof and dzf to svemis. For now, loop over k.

    do k = 1,kemis
      svp(2:i1, 2:j1, k, sv_skip+1:nsv) = svp(2:i1, 2:j1, k, sv_skip+1:nsv)    + &
            ((1. - emistime_s)*svemis(2:i1, 2:j1, k, sv_skip+1:nsv, 1) + &
                   emistime_s *svemis(2:i1, 2:j1, k, sv_skip+1:nsv, 2))/(3600.*rhof(k)*dzf(k)*dx*dy*1e-6) 
    end do

    ! --------------------------------------------------------------------------
    ! Read emission files when neccesary, i.e. simulation reaches half hour mark
    ! after current timestep
    ! --------------------------------------------------------------------------

    if ( rk3step == 3 ) then
      emistime_e = mod(rtimee + rdt + 1800., 3600.)*div3600

      if ( emistime_e < emistime_s ) then
        ! Transfer data from 'ahead-of-modeltime' field to 'past-modeltime' field
        svemis(:,:,:,:,1) = svemis(:,:,:,:,2)
        
        ! Read new 'ahead-of-modeltime' emission field
        if ( datex(4) == 23 ) then
          call reademission(nextday(1), nextday(2), nextday(3),          0, svemis(:,:,:,:,2))
        else
          call reademission(  datex(1),   datex(2),   datex(3), datex(4)+1, svemis(:,:,:,:,2))
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

    deallocate(svemis)

  end subroutine exitemission

end module modemission
