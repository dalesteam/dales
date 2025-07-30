!> \file modmicrodata.f90
!!  Variables necessary for the microphysics

!>
!!  Variables necessary for the microphysics
!>
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

module modmicrodata

  use modprecision, only: field_r

  implicit none

  ! Microphysics scheme options
  integer, parameter :: &
    imicro_none = 0,    & !< No microphysics.
    imicro_drizzle = 1, & !< Drizzle microphyics.
    imicro_bulk = 2,    & !< Double-moment warm microphysics.
    imicro_bin = 3,     & !< Bin microphysics. (why is this still here?)
    imicro_sice = 5,    & !< Single-moment mixed-phase microphysics.
    imicro_sice2 = 6,   & !< Single-moment mixed-phase microphysics (alternative implementation).
    imicro_user = 10,   & !< User-provided microphysics.
    imicro_bulk3 = 11     !< Double-moment mixed-phase microphysics.

  ! Thresholds for statistics
  real(field_r), parameter :: &
    epscloud = 0.01e-3,       &
    epsprec = 3.65e-5,        &
    epsqr = 1.0e-8

  ! User settings
  integer :: imicro = 0     !< Selected scheme.
  logical ::        &
    lstat = .true., & !< Compute intermediate statistics.
    l_rain = .true.   !< Switch for rain calculations.

  ! Cloud settings
  ! Might be used without microphysics (e.g.: radiation), so kept here.
  real(field_r) :: &
    Nc_0 = 70e6,   & !< Cloud droplet number concentration [1/m^3].
    sig_g = 1.34     !< Std. dev. of cloud droplet size distribution.

  ! Indices of rain-related tracers in tracer array. Kept here, because some
  ! statistics use them. Better to switch to using get_tracer_index in the future.
  integer ::  &
    inr = -1, & !< Rain droplet number concentration.
    iqr = -1    !< Rain water mixing ratio.

  real(field_r) :: delt !< Time step size. Better to put this in tstep.f90?

  ! Temp arrays for qt and thl tendencies, and precip
  real(field_r), allocatable :: &
    qtpmcr(:,:,:),              & !< Qt tendency.
    thlpmcr(:,:,:),             & !< Thl tendency.
    precep(:,:,:)                 !< Precipitation [m/s].

end module modmicrodata