!> \file modsampsurfdata.f90
!! Variable definitions for the modsamplingsurf routine

!>
!! This routine should have no dependency on any other routine, save perhaps modglobal or modfields.
!!  \author Xabier Pedruzo
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



module modsampsurfdata

! implicit none

SAVE
  real    :: dtav
  real    :: timeav
  logical :: lsampallsurf   = .true.  !< switch for sampling (on/off)
  logical :: lsampupsurf    = .false. !< switch for conditional sampling updraft (on/off)
  logical :: lsampbuupsurf  = .false. !< switch for conditional sampling buoyant updraft (on/off)
  logical :: lsampclsurf    = .false. !< switch for conditional sampling cloud (on/off)
  logical :: lsampclearsurf    = .false. !< switch for conditional sampling cloud (on/off)
  logical :: lsampcl0o0_0o5surf    = .false. !< switch for conditional sampling shallow (0<cumulative tau<0.5) cloud (on/off)
  logical :: lsampcl0o5_1o0surf    = .false. !< switch for conditional sampling shallow (0<cumulative tau<0.5) cloud (on/off)
  logical :: lsampcl1o0_1o5surf    = .false. !< switch for conditional sampling shallow (0<cumulative tau<0.5) cloud (on/off)
  logical :: lsampcl1o5_2o0surf    = .false. !< switch for conditional sampling shallow (0<cumulative tau<0.5) cloud (on/off)
  logical :: lsampcl2o0_2o5surf    = .false. !< switch for conditional sampling shallow (0<cumulative tau<0.5) cloud (on/off)
  logical :: lsampcl2o5_3o0surf    = .false. !< switch for conditional sampling shallow (0<cumulative tau<0.5) cloud (on/off)
  logical :: lsampcl3o0_3o5surf    = .false. !< switch for conditional sampling shallow (0<cumulative tau<0.5) cloud (on/off)
  logical :: lsampcl3o5_4o0surf    = .false. !< switch for conditional sampling shallow (0<cumulative tau<0.5) cloud (on/off)
  logical :: lsampcl4o0_4o5surf    = .false. !< switch for conditional sampling shallow (0<cumulative tau<0.5) cloud (on/off)
  logical :: lsampcl4o5_5o0surf    = .false. !< switch for conditional sampling shallow (0<cumulative tau<0.5) cloud (on/off)
  logical :: lsampcl5_6surf    = .false. !< switch for conditional sampling shallow (0<cumulative tau<0.5) cloud (on/off)
  logical :: lsampcl6_7surf    = .false. !< switch for conditional sampling shallow (0<cumulative tau<0.5) cloud (on/off)
  logical :: lsampcl7_8surf    = .false. !< switch for conditional sampling shallow (0<cumulative tau<0.5) cloud (on/off)
  logical :: lsampcl8_9surf    = .false. !< switch for conditional sampling shallow (0<cumulative tau<0.5) cloud (on/off)
  logical :: lsampcl9_10surf    = .false. !< switch for conditional sampling shallow (0<cumulative tau<0.5) cloud (on/off)
  logical :: lsampcl10_15surf    = .false. !< switch for conditional sampling shallow (2<cumulative tau<5) cloud (on/off)
  logical :: lsampcl15_20surf    = .false. !< switch for conditional sampling shallow (15<cumulative tau<20) cloud (on/off)
  logical :: lsampclO20surf    = .false. !< switch for conditional sampling deep (cumulative tau>20) cloud (on/off)

end module modsampsurfdata
