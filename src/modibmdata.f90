!> \file modibmdata.f90
!!  Provides variable and paramater values for the grid-conforming Immersed Boundary Method (IBM) 

!> References:
!! (1) Pourquie, M., W.-P. Breugem, and B. J. Boersma, 2009: Some issues related to the use of immersed boundary methods
!! to represent square obstacles. International Journal for Multiscale Computational Engineering, 7 (6), 509–522.
!! (2) Tomas, J., 2016: Obstacle-resolving large-eddy simulation of dispersion in urban environments: Effects of stability and roughness geometry, 
!! Delft University of Technology, Delft, The Netherlands. https://doi.org/10.4233/uuid:5d93a697-be49-4f63-b871-b763bc327139
!>
!!  \author Michael Koene, Delft University of Technology, 2018-2019
!!  \author Stephan de Roode, Delft University of Technology, 2018-2024
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

!! By Michael Koene (MK), TU Delft, section Atmospheric Physics, 28 January 2019
!! cstep: subroutine airslabsum  moved from modmpi to here to avoid mutual dependencies

module modibmdata
  use modprecision, only: field_r
  implicit none
  save

  logical :: lapply_ibm     = .false.        !< Switch to enable immersed boundary method 
  logical :: lwallheat      = .false.        !< Switch to apply lateral heat flux from buildings
  logical :: lpoislast      = .true.         !< Switch to use the Poisson solver after the Immersed boundary method
                                             !  .false. will set the order to: ZeroVelocity -> PoissonSolver -> IBM

  real(field_r)    :: thlwall        = 293.           !< Wall temperature for temperature flux at the sides of the buildings, needed for lateral flux
  real(field_r)    :: thlroof        = 293.           !< Obstacle roof (top) temperature
  real(field_r)    :: qtroof         = 0.             !< Obstacle roof specific humidity
  real(field_r)    :: thlibm         = 293            !< Interior potential temperature of obstacle
  real(field_r)    :: qtibm          = 0.             !< Wall specific humidity for the latent heat flux at the sides and top of the buildings
                                             !< In modsurface it will be set to the saturation value (but this needs to be adapted)
  real(field_r)    :: z0m_wall       = 0.03           !< compare with 0.03 m for open flat terrain, grass, few isolated obstacles
  real(field_r)    :: z0h_wall       = 0.03           !< compare with 0.03 m for open flat terrain, grass, few isolated obstacles

end module modibmdata
