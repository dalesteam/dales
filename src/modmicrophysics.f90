!> \file modmicrophysics.f90
!!  Microphysics abstraction layer.

!>
!!  Microphysics abstraction layer.
!>
!!  Also provides the drizzle routine
!!  \author Hans Cuijpers, IMAU
!!  \author Thijs Heus,MPI-M
!!  \author Steef B\"oing, TU Delft
!!  \todo Documentation
!!  \par Revision list
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
!  Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!



module modmicrophysics
  use modglobal,     only: ifnamopt, checknamelisterror
  use modbulkmicro,  only: initbulkmicro, bulkmicro_read_namelist
  use modbulkmicro3, only: initbulkmicro3, bulkmicro3_read_namelist
  use modmicrodata,  only: imicro, lstat
  use modsimpleice,  only: initsimpleice, simpleice_read_namelist
  use modsimpleice2, only: initsimpleice2
  use modmpi,        only: myid, D_MPI_BCAST, comm3d, print_info_stderr

implicit none

  private

  public :: microphysics_read_namelist
  public :: initmicrophysics
  public :: microsources
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
      case(imicro_bulk)
        call bulkmicro_read_namelist(nml_filename)
      case(imicro_sice, imicro_sice2)
        call simpleice_read_namelist(nml_filename)
      case(imicro_bulk3)
        call bulkmicro3_read_namelist(nml_filename) 
    end select

  end subroutine microphysics_read_namelist

  subroutine initmicrophysics()

    character(len=*), parameter :: routine = modname//'/init_microphysics'

    select case(imicro)
    case(imicro_bulk)
       call initbulkmicro
    case(imicro_sice)
       call initsimpleice
    case(imicro_sice2)
       call initsimpleice2
      case(imicro_bulk3)
        call initbulkmicro3
      case default
        call print_info_stderr(modname, 'invalid option selected for imicro')
        error stop
    end select

  end subroutine initmicrophysics


  subroutine microphysics
    use modmicrodata, only : imicro, imicro_none, imicro_drizzle, imicro_bulk, imicro_bin, imicro_user
    implicit none
    select case (imicro)
    case(imicro_none)
    case(imicro_drizzle)
    case(imicro_bulk)
!       call bulkmicro
    case(imicro_bin)
!       call binmicro
    case(imicro_user)
    end select
  end subroutine microphysics

  subroutine microsources
   use moduser,      only : micro_user
   use modbulkmicro, only : bulkmicro
   use modsimpleice, only : simpleice
   use modsimpleice2, only : simpleice2
   use modmicrodata, only : imicro, imicro_drizzle, imicro_bulk, imicro_bin, &
                            imicro_sice, imicro_sice2, imicro_user, imicro_none, &
                            imicro_bulk3
   use modbulkmicro3, only : bulkmicro3 !#sb3
   use modtimer
!     use modbinmicro,  only : binmicrosources
    implicit none

    call timer_tic('modmicrophysics/microsources', 0)

    select case (imicro)
    case(imicro_none)
    case(imicro_drizzle)
      call drizzle
    case(imicro_bulk)
      call bulkmicro
    case(imicro_bin)
!       call binmicrosources
    case(imicro_sice)
       call simpleice
    case(imicro_sice2)
      call simpleice2
    case(imicro_bulk3)  !#sb3
      call bulkmicro3   !#sb3
    case(imicro_user)
      call micro_user
    end select

    call timer_toc('modmicrophysics/microsources')

  end subroutine microsources

  subroutine exitmicrophysics
    use modbulkmicro, only : exitbulkmicro
    use modsimpleice, only : exitsimpleice
    use modsimpleice2, only : exitsimpleice2
    use modmicrodata, only : imicro, imicro_none, imicro_drizzle, imicro_bin, &
                             imicro_user, imicro_bulk, imicro_sice, imicro_sice2, &
                             imicro_bulk3
    use modbulkmicro3, only : exitbulkmicro3 !#sb3
 !     use modbinmicro,  only : exitbinmicro
    implicit none

     select case (imicro)
     case(imicro_none)
     case(imicro_drizzle)
     case(imicro_bulk)
!       call exitbulkmicro
     case(imicro_bin)
!       call exitbinmicro
     case(imicro_user)
     case(imicro_sice)
        call exitsimpleice
     case(imicro_sice2)
      call exitsimpleice2
     case(imicro_bulk3)    !#sb3
      call exitbulkmicro3  !#sb3
  end select
  end subroutine exitmicrophysics

 subroutine drizzle

!-----------------------------------------------------------------|
!                                                                 |
!      Hans Cuijpers   I.M.A.U.  23 May 1995                      |
!                                                                 |
!     purpose.                                                    |
!     --------                                                    |
!                                                                 |
!      Calculates gravitational settling (or rainfall rate)       |
!                                                                 |
!**   interface.                                                  |
!     ----------                                                  |
!                                                                 |
!     *drizzle* is called from *program*.                         |
!                                                                 |
!-----------------------------------------------------------------|

  use modglobal, only : i1,j1,kmax,rlv,cp,dzf,pi
  use modfields, only : qtp,ql0,thlp,rhof,exnf
  use modmicrodata, only : c_st, Nc_0, rhow, sig_g
  implicit none
  real :: sedc,csed
  integer :: i, j, k
    csed = c_St*(3./(4.*pi*rhow))**(2./3.)*exp(5.*log(sig_g)**2.)*Nc_0**(-2./3.)
  sedc = 0.

  do k=1,kmax
  do j=2,j1
  do i=2,i1
  if (ql0(i,j,k)>0.0) then
    sedc= csed*((ql0(i,j,k+1)*rhof(k+1))**(5./3.)-(ql0(i,j,k)*rhof(k))**(5./3.))/(dzf(k)*rhof(k))
    qtp(i,j,k) = qtp(i,j,k) + sedc
    thlp(i,j,k) = thlp(i,j,k) - (rlv/(cp*exnf(k)))*sedc
  endif
  enddo
  enddo
  enddo
  return
  end subroutine drizzle

end module modmicrophysics