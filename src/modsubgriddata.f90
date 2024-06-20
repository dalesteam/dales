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
use modprecision, only: field_r
implicit none
save
! private
  logical :: ldelta       = .false. !<  switch for subgrid length formulation (on/off)
  logical :: lmason       = .false. !<  switch for decreased length scale near the surface
  logical :: lsmagorinsky = .false. !<  switch for smagorinsky subgrid scheme
  logical :: lanisotrop   = .false. !<  switch for anisotropic diffusion

  real(field_r) :: cf      = 2.5  !< filter constant
  real(field_r) :: Rigc    = 0.25 !< critical Richardson number
  real(field_r) :: Prandtl = (1.0/3.0)
  real(field_r) :: cm      = 0.12
  real(field_r) :: cn      = 0.76
  real(field_r) :: ch1     = 1.
  real(field_r) :: ch2     = 2.
  real(field_r) :: ch                   ! ch = ch1 + ch2
  real(field_r) :: ce1     = 0.19
  real(field_r) :: ce2     = 0.51
  real(field_r) :: cs      = -1
  real(field_r) :: nmason  = 2.   !< exponent in Mason correction function
  real(field_r) :: alpha_kolm  = 1.5     !< factor in Kolmogorov expression for spectral energy
  real(field_r) :: beta_kolm   = 1.      !< factor in Kolmogorov relation for temperature spectrum
  logical :: sgs_surface_fix = .false.  !< which fix to apply to coupling of SGSTKE to surface


  real(field_r), allocatable :: ekm(:,:,:)   !< k-coefficient for momentum
  real(field_r), allocatable :: ekh(:,:,:)   !< k-coefficient for heat and q_tot
  real(field_r), allocatable :: sbdiss(:,:,:)!< dissiation
  real(field_r), allocatable :: sbshr(:,:,:) !< shear production
  real(field_r), allocatable :: sbbuo(:,:,:) !< buoyancy production / destruction
  real(field_r), allocatable :: zlt(:,:,:)   !< filter width

  real(field_r), allocatable :: csz(:)       !< Smagorinsky constant

  real(field_r), allocatable :: anis_fac(:)  !< grid anisotropy factor 
end module modsubgriddata

