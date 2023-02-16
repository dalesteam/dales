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

  implicit none

  save
  public  :: inittracers

  integer  :: iname
  logical  :: ltracers = .false. ! tracer switch
  character(len = 6), dimension(100) :: & 
              tracernames = (/ ('      ', iname=1, 100) /) ! list with scalar names,

contains

  !> Initialize tracer definition.
  !!
  !! Read the namelist NAMTRACERS from the namoptions file, distribute
  !! the parameters to all processes and allocate the tracer (SV) arrays.
  subroutine inittracers
    ! read namelist    
    ! init tracer array

    use modglobal, only : ifnamopt, fname_options, checknamelisterror
    use modmpi,    only : myid, comm3d, mpierr, mpi_logical, mpi_integer, my_real, mpi_character

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
    ! call MPI_BCAST(nsv,                  1, mpi_integer, 0, comm3d, mpierr)
    call mpi_bcast(tracernames(1:100), 100, mpi_character, 0, comm3d, ierr)

  end subroutine inittracers


  !! this function is inserted here as an example of how to use "findval"
  !> Retrieve the three letter LU class acronym from the tile and return the corresponding index in DEPAC arrays
  !!
  !! This function is an interface between the current LSM implementation in DALES and the one used in DEPAC.
  !! Once both DALES and DEPAC incorporate an LSM model that is not hard linked to data, this function is no longer needed.
  !!
  !! @param[in] luclass The LU class acronym from the tile (`lushort`)
  !! @returns The corresponding index in DEPAC arrays defined in `depac_lu.inc`  
  ! pure integer function get_depac_luindex(luclass)
  !   implicit none
  !   character(len=3), intent(in) :: luclass
  !   integer :: idx(1)
  !   character(len=3) :: depac_lu

  !   select case (luclass)
  !     case ('fcd', 'fce')
  !       depac_lu = 'cnf'
  !     case ('fbd', 'fbe')
  !       depac_lu = 'dec'
  !     case ('aqu')
  !       depac_lu = 'wai'
  !     case ('sem')
  !       depac_lu = 'grs'
  !     case ('brn')
  !       depac_lu = 'dsr'
  !     case default
  !       depac_lu = luclass
  !   end select

  !   idx = findloc(lu_name_abbr, depac_lu)
  !   get_depac_luindex = idx(1)
  ! end function get_depac_luindex

  !> Find a value in an array of values based on a key in an array of keys
  !!
  !! Lookup a value in one array (`values`) at the index position of `key` in the array `keys`.
  !! Provide the optional default value (`defltvalue`), that is returned when `key` isn't found.
  !! When no default is given, zero is returned.
  !!
  !! @param[in] key The key to be looked up
  !! @param[in] keys The array of keys
  !! @param[in] values The array of values
  !! @param[in] defltvalue The default value (optional, default 0.0)
  !! @return The value corresponding to the key
  ! pure real function findval(key, keys, values, defltvalue)
  !   implicit none
  !   character(*), intent(in) :: key 
  !   character(*), intent(in) :: keys(:)
  !   real, intent(in) :: values(:)
  !   real, intent(in), optional :: defltvalue
  !   character(6) :: fkey
  !   integer :: idx(1)

  !   fkey = trim(key)

  !   idx = findloc(keys, fkey)
  !   if (idx(1) /= 0 .and. idx(1) <= size(values)) then
  !     findval = values(idx(1))
  !   else
  !     if (present(defltvalue)) then
  !       findval = defltvalue 
  !     else
  !       findval = 0.0 
  !     end if
  !   end if
  ! end function findval

end module modtracers