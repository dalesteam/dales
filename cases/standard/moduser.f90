!> \file moduser.f90
!! A dummy file for cases where one wants additional forcings
!----------------------------------------------------------------------------
! This file is part of DALES.
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
! Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!----------------------------------------------------------------------------
!
!
module moduser
use modlogging, only: finish

character(len=*), parameter :: modname = 'moduser'

contains
subroutine initsurf_user
  implicit none
  character(len=*), parameter :: routine = modname//'/initsurf_user'
  call finish(routine, "surf_user invoked (isurf = 10) but a case-specific moduser.f90 was not compiled in.")
end subroutine initsurf_user

subroutine force_user
  implicit none
  character(len=*), parameter :: routine = modname//'/force_user'
  call finish(routine, "force_user invoked but a case-specific moduser.f90 was not compiled in.")
end subroutine force_user

subroutine rad_user
  implicit none
  character(len=*), parameter :: routine = modname//'/rad_user'
  call finish(routine, "rad_user invoked but a case-specific moduser.f90 was not compiled in.")
end subroutine rad_user

subroutine micro_user
  implicit none
  character(len=*), parameter :: routine = modname//'/micro_user'
    call finish(routine, "micro_user invoked but a case-specific moduser.f90 was not compiled in.")
end subroutine micro_user

subroutine surf_user
  implicit none
  character(len=*), parameter :: routine = modname//'/surf_user'
  call finish(routine, "surf_user invoked (isurf = 10) but a case-specific moduser.f90 was not compiled in.")
end subroutine surf_user

end module moduser
