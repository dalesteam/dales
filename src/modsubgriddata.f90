!!> \file modsubdata.f90
!!!  Provides variable definitions for Calculates and applies the Sub Filter Scale diffusion
!
!>
!!  Calculates and applies the Sub Filter Scale diffusion
!>
!!  \author Pier Siebesma, K.N.M.I.
!!  \author Stephan de Roode,TU Delft
!!  \author Chiel van Heerwaarden, Wageningen U.R.
!!  \author Thijs Heus,MPI-M
!!  \par Revision list
!!  \todo Documentation
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

module modsubgriddata
implicit none
save
! private
  logical :: ldelta       = .false. !<  switch for subgrid length formulation (on/off)
  logical :: lmason       = .false. !<  switch for decreased length scale near the surface
  logical :: lD80R        = .false. !<  switch for D80-R scheme
  logical :: lsmagorinsky = .false. !<  switch for smagorinsky subgrid scheme
  real :: cf      = 2.5  !< filter constant
  real :: Rigc    = 0.25 !< critical Richardson number
  real :: Prandtl = (1.0/3.0)
  real :: cm      = 0.12
  real :: cn      = 0.76
  real :: ch1     = 1.
  real :: ch2     = 2.
  real :: ce1     = 0.19
  real :: ce2     = 0.51
  real :: cs      = -1
  real :: nmason  = 2.   !< exponent in Mason correction function
  real :: alpha_kolm  = 1.5     !< factor in Kolmogorov expression for spectral energy
  real :: beta_kolm   = 1.      !< factor in Kolmogorov relation for temperature spectrum
  logical :: sgs_surface_fix = .false.  !< which fix to apply to coupling of SGSTKE to surface


  real, allocatable :: ekm(:,:,:)   !< k-coefficient for momentum
  real, allocatable :: ekh(:,:,:)   !< k-coefficient for heat and q_tot
  real, allocatable :: sbdiss(:,:,:)!< dissiation
  real, allocatable :: sbshr(:,:,:) !< shear production
  real, allocatable :: sbbuo(:,:,:) !< buoyancy production / destruction
  real, allocatable :: zlt(:,:,:)   !< filter width

  real, allocatable :: csz(:)       !< Smagorinsky constant

end module modsubgriddata

