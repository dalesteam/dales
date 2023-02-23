!> \file modemisdata.f90
!! Variable definitions and auxilary routines for emissions

!>
!! Variable definitions and auxilary routines for emission
!>
!! This routine should have no dependency on any other routine, save
!perhaps modglobal or modfields.
!!  \author Marco de Bruine, VU
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
!  Copyright 1993-2020 Delft University of Technology, Wageningen
!  University, Utrecht University, KNMI, Vrije Universiteit Amsterdam
!
! TODO INTEGRATION WITH MODCHEM, possibly switch?
! For example, nchem, firstchem etc, but also structure for species
! 'location', i.e. switch scalar field represents which species?

module modemisdata
  
  implicit none
  save

  integer  :: iname
  ! ---------------------------------------------------------------------!
  ! Namelist variables                                                   !
  ! ---------------------------------------------------------------------!
   
  logical  :: l_emission = .false. ! scalar emission switch
  integer  :: kemis    = -1, &     ! no. of layers to include for emission
              svskip   =  0, &     ! no. scalars to exclude for emission
              nemis    = 0         ! no. of emitted scalars  

  character(len = 6), dimension(100) :: & 
              emisnames = (/ ('      ', iname=1, 100) /) ! list with scalar names,
                          ! each name must(!) be 6 characters long for now  

  ! Interaction with AGs ------------------------------------------------
  integer :: svco2ags = -1       ! Scalar field number for AGs soil respiration
  integer :: svco2veg = -1       ! Scalar field number for AGs net CO2 flux
  integer :: svco2sum = -1       ! Scalar field which holds the sum of CO2
  integer, allocatable :: co2fields(:) ! Array defnining co2 fields for AGs

  ! ---------------------------------------------------------------------!
  ! Main variables                                                       !
  ! ---------------------------------------------------------------------!

  real, allocatable :: emisfield(:,:,:,:,:) !array for emission fields

end module modemisdata
