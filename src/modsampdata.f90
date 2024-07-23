!> \file modsampdata.f90
!! Variable definitions for the modsampling and modsamptend routines

!>
!! This routine should have no dependency on any other routine, save perhaps modglobal or modfields.
!!  \author Huug Ouwersloot
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



module modsampdata

! implicit none

SAVE
  real    :: dtav
  real    :: timeav
  logical :: lsampall   = .true.  !< switch for sampling (on/off)
  logical :: lsampcl    = .false. !< switch for conditional sampling cloud (on/off)
  logical :: lsampco    = .false. !< switch for conditional sampling core (on/off)
  logical :: lsampup    = .false. !< switch for conditional sampling updraft (on/off)
  logical :: lsampbuup  = .false. !< switch for conditional sampling buoyant updraft (on/off)
  logical :: lsampcldup = .false. !< switch for condtional sampling cloudy updraft (on/off)
  logical :: lsamptend  = .false. !< switch to also sample tendencies
  logical :: lsamptendu = .false. !< switch to sample u tendencies
  logical :: lsamptendv = .false. !< switch to sample v tendencies
  logical :: lsamptendw = .false. !< switch to sample w tendencies
  logical :: lsamptendthl = .false. !< switch to sample thl tendencies
  logical :: lsamptendqt = .false. !< switch to sample qt tendencies
  logical :: lsamptendqr = .false. !< switch to sample qr tendencies
  logical :: lsamptendnr = .false. !< switch to sample nr tendencies
  logical :: lprocblock = .false. !< switch to write (so far only tendencies) per processor block instead of domain-averaged
  logical :: ltenddec = .false. !< switch to get variables needed to scale-decompose (processor-averaged) advective tendencies

end module modsampdata
