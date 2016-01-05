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
! Arnold Moene added 14-1-2015
  logical :: sgs_surface_fix = .true.  !< which fix to apply to coupling of SGSTKE to surface


  real, allocatable :: ekm(:,:,:)   !< k-coefficient for momentum
  real, allocatable :: ekh(:,:,:)   !< k-coefficient for heat and q_tot
  real, allocatable :: sbdiss(:,:,:)!< dissiation
  real, allocatable :: sbshr(:,:,:) !< shear production
  real, allocatable :: sbbuo(:,:,:) !< buoyancy production / destruction
  real, allocatable :: zlt(:,:,:)   !< filter width
  
  real, allocatable :: csz(:)       !< Smagorinsky constant

  !victor: onderstaande variabelen zijn allemaal nodig voor het subgrid model van Sullivan
  
  real, allocatable :: suu(:,:,:) !< element 11 of Strain tensor = 1/2(du_i/dx_j + du_j/dx_i) 
  real, allocatable :: svv(:,:,:) !< element 22 of Strain tensor = 1/2(du_i/dx_j + du_j/dx_i)
  real, allocatable :: sww(:,:,:) !< element 33 of Strain tensor = 1/2(du_i/dx_j + du_j/dx_i)  
  real, allocatable :: suv(:,:,:) !< element 12 of Strain tensor = 1/2(du_i/dx_j + du_j/dx_i)
  real, allocatable :: suw(:,:,:) !< element 13 of Strain tensor = 1/2(du_i/dx_j + du_j/dx_i)
  real, allocatable :: svw(:,:,:) !< element 23 of Strain tensor = 1/2(du_i/dx_j + du_j/dx_i)
  
  real, allocatable :: suum(:) !< element 11 of slab averaged strain
  real, allocatable :: svvm(:) !< element 22 of slab averaged strain
  real, allocatable :: swwm(:) !< element 33 of slab averaged strain  
  real, allocatable :: suvm(:) !< element 12 of slab averaged strain
  real, allocatable :: suwm(:) !< element 13 of slab averaged strain
  real, allocatable :: svwm(:) !< element 23 of slab averaged strain
   
  real, allocatable :: svaruu(:) !< variance of 11 component of slabstrain
  real, allocatable :: svarvv(:)
  real, allocatable :: svarww(:)
  real, allocatable :: svaruv(:)
  real, allocatable :: svaruw(:)
  real, allocatable :: svarvw(:)
  
  real, allocatable :: svaruul(:,:,:) !< local variance of 11 component of slabstrain
  real, allocatable :: svarvvl(:,:,:)
  real, allocatable :: svarwwl(:,:,:)
  real, allocatable :: svaruvl(:,:,:)
  real, allocatable :: svaruwl(:,:,:)
  real, allocatable :: svarvwl(:,:,:)
  
  
  
  real, allocatable :: sacc(:) !< S' , defined as sqrt( 2 <Sij-<Sij>><Sij-<Sij>>)
  real, allocatable :: sslab(:) !< <S>, defined by sqrt(2<Sij><Sij>)
  real, allocatable :: gamma(:) ! < isotropy factor, defined by S'/(S'+<S>)
  real, allocatable :: sdemean(:,:,:) !<  (Sij-<Sij>)^2, as used in the Sullivan subgrid shear production term
  
  real, allocatable :: ekm2(:)
  real, allocatable :: uwr(:,:)
  real, allocatable :: vwr(:,:)
  
  real:: ustar2
  real:: phimzf2
  
  real:: zisullivan

end module modsubgriddata
