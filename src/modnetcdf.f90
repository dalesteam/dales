!> \file modnetcdf.f90
!!  Convenience functions for working with NetCDF.
!>
!!  \author Caspar Jungbacker, Delft University of Technology
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
! Copyright 2024 Delft University of Technology
!
module modnetcdf

  use netcdf

  implicit none

  public
  
  public :: check

contains
 
  !> Checks return code of netcdf calls
  !!
  !! If an error occurs, stops the program and print information about
  !! the location of the error
  subroutine check(status, file, line)

    use modmpi, only: myid

    integer,      intent(in) :: status, line
    character(*), intent(in) :: file

    if (status /= nf90_noerr) then
      if (myid == 0) then
        write(*,*) "NetCDF error in: ", file, " on line: ", line
        write(*,*) trim(nf90_strerror(status))
      end if
      stop
    end if

  end subroutine check

end module modnetcdf
