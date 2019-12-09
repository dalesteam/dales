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

  subroutine initemis

  ! ----------------------------------------------------------------------
  ! Reading of emission files
  ! EXPAND
  ! ----------------------------------------------------------------------

    use mpi
    use netcdf
    use modmpi,      only : myidx, myidy
    use modglobal,   only : i1, j1, i2, j2, imax, jmax, nsv

    implicit none

    character (len = *), parameter :: fname = "emission_input.nc"

    integer, parameter :: ndim = 3 
    integer            :: start(ndim), count(ndim)
    integer            :: ncid, varid
    integer            :: n

    allocate(svemis(i2,j2,nsv))

    call check( nf90_open     ( fname, IOR(NF90_NOWRITE, NF90_MPIIO), ncid, &
                                comm = MPI_COMM_WORLD, info = MPI_INFO_NULL) )
    call check( nf90_inq_varid( ncid, "emission", varid) )
    call check( nf90_get_var  ( ncid, varid, svemis(2:i1,2:i1,1), & 
                                start = (/1 + myidx * imax, 1 + myidy * jmax, 1/), &
                                count = (/imax, jmax, 1/) ) )
    call check( nf90_close    ( ncid ) )

  contains

    subroutine check(status)
      integer, intent(in) :: status
    
      if(status /= nf90_noerr) then 
        print *, trim(nf90_strerror(status))
        stop 2
      end if
    end subroutine check

  end subroutine initemis

  ! ----
  ! Do calculations on emission data and transfer to svp array
  ! ----
  subroutine emission

  use modfields, only : svp
  use modglobal, only : i1, j1, ih, jh, nsv
  use modfields, only : rhof

    implicit none

    svp(2:i1,2:j1,1,1:nsv) = svp(2:i1,2:j1,1,1:nsv) + svemis(2:i1,2:j1,1:nsv)/rhof(1) 

  end subroutine emission

  ! ----
  ! Cleanup after run.
  ! ----
  subroutine exitemission
  
    implicit none

    deallocate(svemis)

  end subroutine exitemission

end module modemission
