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

    use modglobal,    only : i2, j2, nsv
    use moddatetime,  only : datex, prevday, nextday

    implicit none

    allocate(svemis_a(i2,j2,nsv))
    allocate(svemis_b(i2,j2,nsv))

    if (datex(5) >= 30) then
      call reademission(    datex(1),   datex(2),   datex(3),   datex(4), svemis_a)
      
      if (datex(4) == 23) then  
        call reademission(nextday(1), nextday(2), nextday(3),          0, svemis_b)
      else
        call reademission(  datex(1),   datex(2),   datex(3), datex(4)+1, svemis_b)
      endif

    else
      call reademission(    datex(1),   datex(2),   datex(3),   datex(4), svemis_b) 
     
      if (datex(4) == 0) then
        call reademission(prevday(1), prevday(2), prevday(3),         23, svemis_a)
      else
        call reademission(  datex(1),   datex(2),   datex(3), datex(4)-1, svemis_a)
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
    real, intent(out)    :: emisfield(i2,j2,nsv)

    integer, parameter   :: ndim = 3 
    integer              :: start(ndim), count(ndim)
    integer              :: ncid, varid
    integer              :: isv
 
    character(len=12)    :: sdatetime

    ! Create string from given date
    write(sdatetime, "(I0.4,2I0.2,2I0.2)") iyear, imonth, iday, ihour, 0 

    write(6,"(A18, A12)") "Reading emission: ", sdatetime

    do isv = 1, nsv

      if (emislist(isv)) then
        call check( nf90_open( svlist(isv)//'_emis_'//sdatetime//'.nc', IOR(NF90_NOWRITE, NF90_MPIIO), &
                                ncid, comm = comm3d, info = MPI_INFO_NULL) )
        call check( nf90_inq_varid( ncid, svlist(isv), varid) )
        call check( nf90_get_var  ( ncid, varid, emisfield(2:i1,2:j1,isv), & 
                                    start = (/1 + myidx * imax, 1 + myidy * jmax, 1/), &
                                    count = (/imax, jmax, 1/) ) )
        call check( nf90_close( ncid ) )
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
  !    gridbx size and air density AND apply a factor of 1e6.
  ! 
  ! TODO
  ! 1. MDB Align properly with non-chem tracers, i.e. cloud scalars from e.g.
  ! microphysics.
  ! ----------------------------------------------------------------------

    use modfields,   only : svp
    use modglobal,   only : i1, j1, ih, jh, nsv, &
                            rdt, rtimee, rk3step, &
                            dzf, dx, dy    
    use modfields,   only : rhof
    use moddatetime, only : datex, nextday
    
    implicit none

    real            :: emistime_s, emistime_e ! Emission timers
    real, parameter :: div3600 = 1./3600.     ! Quick division

    ! --------------------------------------------------------------------------
    ! Interpolate and apply emission
    ! --------------------------------------------------------------------------
    emistime_s = mod(rtimee +       1800., 3600.)*div3600

    svp(2:i1,2:j1,1,1:nsv) = svp(2:i1,2:j1,1,1:nsv) + &
                             ((1. - emistime_s)*svemis_a(2:i1, 2:j1, 1:nsv) + &
                                    emistime_s *svemis_b(2:i1, 2:j1, 1:nsv))/(3600.*rhof(1)*dzf(1)*dx*dy*1e-6) 

    ! --------------------------------------------------------------------------
    ! Read emission files when neccesary, i.e. simulation reaches half hour mark
    ! after current timestep
    ! --------------------------------------------------------------------------

    if ( rk3step == 3 ) then
      emistime_e = mod(rtimee + rdt + 1800., 3600.)*div3600

      if ( emistime_e < emistime_s ) then
        ! Transfer data from 'ahead-of-modeltime' field to 'past-modeltime' field
        svemis_a = svemis_b
        
        ! Read new 'ahead-of-modeltime' emission field
        if ( datex(4) == 23 ) then
          call reademission(nextday(1), nextday(2), nextday(3),          0, svemis_b)
        else
          call reademission(  datex(1),   datex(2),   datex(3), datex(4)+1, svemis_b)
        endif
      endif

    endif

  end subroutine emission

  ! --------------------------------------------------------------------------
  ! Cleanup after run.
  ! --------------------------------------------------------------------------
  subroutine exitemission
  
    implicit none

    deallocate(svemis_a, svemis_b)

  end subroutine exitemission

end module modemission
