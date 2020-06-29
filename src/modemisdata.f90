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
!  University, Utrecht University, KNMI
!


module modemisdata
  
  implicit none
  save

  real, allocatable :: &
    svemis_a  (:,:,:), & 
    svemis_b  (:,:,:) 

  ! TODO INTEGRATION WITH MODCHEM, possibly switch?
  ! For example, nchem, firstchem etc, but also structure for species
  ! 'location', i.e. switch scalar field represents which species?
  
  character (len = 3), dimension(2) :: svlist   = (/'nr',    'qr',    'co2',  'ch4'  /)
  logical,             dimension(2) :: emislist = (/.false., .false., .true., .false./)

end module modemisdata
