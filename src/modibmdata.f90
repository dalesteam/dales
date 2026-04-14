!> \file modibmdata.f90
!! Provides variable and paramater values for the grid-conforming Immersed Boundary Method (IBM) 

!>
!!  \author Michael Koene, Delft University of Technology, 2018-2019
!!  \author Stephan de Roode, Delft University of Technology, 2018-2024
!!  \author Steven van der Linden, Delft University of Technology, 2025-
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
! Copyright 2025 Delft University of Technology
!

module modibmdata
  
  use modprecision, only: field_r
  implicit none
  save

  logical :: lapply_ibm     = .false.        !< Switch to enable immersed boundary method 
  logical :: lwallheat      = .false.        !< Switch to apply lateral heat flux from buildings
  logical :: lpoislast      = .true.         !< Switch to use the Poisson solver after the Immersed boundary method
                                             !  .false. will set the order to: ZeroVelocity -> PoissonSolver -> IBM

  real(field_r)    :: thlwall        = 293.           !< Wall temperature at the sides of the buildings [K]
  real(field_r)    :: qtwall         = 0.             !< Wall specific humidity [kg/kg]
  real(field_r)    :: thlroof        = 293.           !< Obstacle roof (top) temperature [K]
  real(field_r)    :: qtroof         = 0.             !< Obstacle roof specific humidity [kg/kg]
  real(field_r)    :: thlibm         = 293            !< Interior potential temperature of obstacle [K]
  real(field_r)    :: qtibm          = 0.             !< Interior specific humidity of obstacle [kg/kg]

  real(field_r)    :: z0m_wall       = 0.03           !< Roughness length for momentum at walls [m]
  real(field_r)    :: z0h_wall       = 0.03           !< Roughness length for heat/scalars at walls [m]
                                                      ! compare with 0.03 m for open flat terrain, grass, few isolated obstacles
end module modibmdata
