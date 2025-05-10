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
! Copyright 1993-2025 The DALES team.
!

!> Microphysics abstraction layer.
!!  \author Hans Cuijpers, IMAU
!!  \author Thijs Heus,MPI-M
!!  \author Steef B\"oing, TU Delft
module modmicrophysics

  use modglobal,     only: ifnamopt, checknamelisterror
  use moddrizzle,    only: drizzle
  use modbulkmicro,  only: initbulkmicro, bulkmicro_read_namelist, &
                           exitbulkmicro, bulkmicro
  use modbulkmicro3, only: initbulkmicro3, bulkmicro3_read_namelist, &
                           exitbulkmicro3, bulkmicro3
  use modmicrodata,  only: imicro, lstat
  use modsimpleice,  only: initsimpleice, simpleice_read_namelist, &
                           exitsimpleice, simpleice
  use modsimpleice2, only: initsimpleice2, exitsimpleice2, simpleice2
  use modmpi,        only: myid, D_MPI_BCAST, comm3d, print_info_stderr
  use modtimer,      only: timer_tic, timer_toc
  use moduser,       only: micro_user

  implicit none

  private

  public :: microphysics_read_namelist
  public :: initmicrophysics
  public :: microphysics
  public :: exitmicrophysics

  character(len=*), parameter :: modname = 'modmicrophysics'

  integer, parameter :: &
    imicro_none = 0,    & !< No microphysics.
    imicro_drizzle = 1, & !< Drizzle microphyics.
    imicro_bulk = 2,    & !< Double-moment warm microphysics.
    imicro_sice = 5,    & !< Single-moment mixed-phase microphysics.
    imicro_sice2 = 6,   & !< Single-moment mixed-phase microphysics (alternative implementation).
    imicro_user = 10,   & !< User-provided microphysics.
    imicro_bulk3 = 11     !< Double-moment mixed-phase microphysics.

  namelist /nammicrophysics/ imicro, lstat

contains

  !> Read microphysics namelist entry and broadcast settings.
  subroutine microphysics_read_namelist(nml_filename)

    character(len=*), intent(in) :: nml_filename

    integer :: ierr

    if (myid == 0) then
      open(ifnamopt, file=nml_filename, status='old', iostat=ierr)
      read(ifnamopt, nammicrophysics, iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'nammicrophysics')
      close(ifnamopt)
    end if

    call D_MPI_BCAST(imicro, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(lstat, 1, 0, comm3d, ierr)

    ! Read namelist of selected microphysics scheme
    select case(imicro)
      case(imicro_none, imicro_drizzle, imicro_user)
        ! Do nothing
      case(imicro_bulk)
        call bulkmicro_read_namelist(nml_filename)
      case(imicro_sice, imicro_sice2) ! simpleice2 uses the same namelist
        call simpleice_read_namelist(nml_filename)
      case(imicro_bulk3)
        call bulkmicro3_read_namelist(nml_filename) 
      case default
        call print_info_stderr(modname, 'invalid option selected for imicro')
        error stop
    end select

  end subroutine microphysics_read_namelist

  !> Call the initialization routine of the selected microphysical scheme.
  subroutine initmicrophysics()

    select case(imicro)
      case(imicro_bulk)
        call initbulkmicro
      case(imicro_sice)
        call initsimpleice
      case(imicro_sice2)
        call initsimpleice2
      case(imicro_bulk3)
        call initbulkmicro3
    end select

  end subroutine initmicrophysics

  !> Do the microphysics.
  subroutine microphysics

    character(len=*), parameter :: routine = modname//'/microphysics'

    call timer_tic(routine, 0)

    select case (imicro)
      case(imicro_drizzle)
        call drizzle
      case(imicro_bulk)
        call bulkmicro
      case(imicro_sice)
         call simpleice
      case(imicro_sice2)
        call simpleice2
      case(imicro_bulk3)
        call bulkmicro3
      case(imicro_user)
        call micro_user
    end select

    call timer_toc(routine)

  end subroutine microphysics

  !> Calls the clean-up routine for the selected microphysical scheme.
  subroutine exitmicrophysics

    select case (imicro)
      case(imicro_bulk)
        call exitbulkmicro
      case(imicro_sice)
        call exitsimpleice
      case(imicro_sice2)
        call exitsimpleice2
      case(imicro_bulk3)
        call exitbulkmicro3
    end select

  end subroutine exitmicrophysics

end module modmicrophysics