!> \file modraddata.f90
!! Variable definitions and auxilary routines for radiation

!>
!! Variable definitions and auxilary routines for radiation.
!>
!! This routine should have no dependency on any other routine, save perhaps modglobal or modfields.
!!  \author Thijs Heus, MPI-M
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



module modraddata

! implicit none
  use modglobal, only : longint

SAVE

  integer, parameter :: irad_none  = 0   !< 0=no radiation
  integer, parameter :: irad_full  = 1   !< 1=full radiation
  integer, parameter :: irad_par   = 2   !< 2=parameterized radiation
  integer, parameter :: irad_lsm   = 3   !< 3=simple surface radiation for land surface model
  integer, parameter :: irad_user  = 10  !< 10=user specified radiation

  logical :: rad_ls      = .true.   !< prescribed radiative forcing
  logical :: rad_longw   = .true.   !< parameterized longwave radiative forcing
  logical :: rad_shortw  = .true.   !< parameterized shortwave radiative forcing
  logical :: rad_smoke   = .false.  !< longwave divergence for smoke cloud
  logical :: useMcICA    = .true.   !< Use the Monte Carlo Independent Column Approach

  real              :: timerad = 0 !<  timescale of the radiation scheme
  integer(kind=longint)           :: itimerad = 0 !<  timescale of the radiation scheme
  integer (kind=longint)          :: tnext   = 0 !<  time of the first upcoming call of the radiation scheme
  real :: rka        = 130.   !< extinction coefficient in radpar scheme
  real :: dlwtop     = 74.    !< longwave radiative flux divergence at top of domain
  real :: dlwbot     = 0.     !< longwave radiative flux divergence near the surface
  real :: sw0        = 1100.0 !< direct component at top of the cloud (W/m^2), diffuse not possible
  real :: gc         = 0.85   !< asymmetry factor of droplet scattering angle distribution
  real :: reff       = 1.e-5  !< cloud droplet effective radius (m)
  integer :: isvsmoke = 1     !< number of passive scalar to be used for optical depth calculation
  integer :: iradiation = irad_none !< Selection parameter for type of radiation scheme
  integer :: irad    = -1  !< Deprecated selection parameter for the type of radiation scheme


  real mu                    !< cosine of the solar zenith angle

  real, allocatable :: thlprad(:,:,:)!<   the radiative tendencies
  real, allocatable :: swd(:,:,:)    !<   shortwave downward radiative flux
  real, allocatable :: swu(:,:,:)    !<   shortwave upward radiative flux
  real, allocatable :: lwd(:,:,:)    !<   longwave downward radiative flux
  real, allocatable :: lwu(:,:,:)    !<   longwave upward radiative flux
  real, allocatable :: swuToA(:,:),swdToA(:,:),lwuToA(:,:),lwdToA(:,:) !< Top of the atmosphere radiative fluxes

contains
!< Calculation of the cosine of the zenith angle
!< \param time UTC Time of the simulation
!< \param xday Day at the start of the simulation
!< \param xlat Latitude of the domain
!< \param xlon Longitude of the domain
  real function zenith(time, xday, xlat,xlon)
    use modglobal, only : pi
!     implicit none
    real, intent(in) :: time, xday, xlat, xlon
    real :: phi,el,obliq,xlam,declin,hora
    real :: day,daytime
    day    = xday + floor(time/86400.)
    daytime= mod(time,86400.)

    phi    = xlat * pi/180.
    el     = xlon * pi/180.
    obliq  = 23.45 * pi/180.
    xlam   = 4.88 + 0.0172 * day
    declin = asin(sin(obliq)*sin(xlam))
    hora   = el-pi + 2.*pi*(daytime/86400.)
    zenith = max(0.,sin(declin)*sin(phi)+cos(declin)*cos(phi)* &
                                                         cos(hora))
  end function zenith

!< Calculation of the albedo at sea
!< From the RRTMG scheme interface by P. Blossey.
  subroutine par_albedo(coszrs,albdir)!,albdif)
  !-----------------------------------------------------------------------
  ! Computes surface albedos over ocean 
  ! and the surface (added by Marat Khairoutdinov, remove by Johan)
  !
  !  Uses solar zenith angle to compute albedo for direct
  !  radiation; diffuse radiation values constant; albedo
  !  independent of spectral interval and other physical
  !  factors such as ocean surface wind speed.
  !
  ! For more details , see Briegleb, Bruce P., 1992: Delta-Eddington
  ! Approximation for Solar Radiation in the NCAR Community Climate Model,
  ! Journal of Geophysical Research, Vol 97, D7, pp7603-7612).
  !-----------------------------------------------------------------------
   implicit none

   real, intent(in)            :: coszrs     ! Cosine of the solar zenith angle
   real, intent(out)           :: albdir     ! Srf alb for direct rad 
!   real, intent(out), optional :: albdif     ! Srf alb for diffuse rad

   real :: eps = 1e-5

   if (coszrs <= eps) then
     albdir = 0.
   else
     albdir = ( .026 / (coszrs**1.7 + .065)) + &
               (.15*(coszrs - 0.10) * (coszrs - 0.50) * (coszrs - 1.00) )
   endif
!   albdif = 0.06

  end subroutine par_albedo

end module modraddata
