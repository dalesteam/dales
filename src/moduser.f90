!> \file moduser.f90
!! A dummy file for cases where one wants additional forcings
!----------------------------------------------------------------------------
! This file is part of DALES.
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
! Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!----------------------------------------------------------------------------
!
!
module moduser

contains
subroutine force_user
  implicit none
  end subroutine force_user

subroutine rad_user
  implicit none
end subroutine rad_user

subroutine micro_user
  implicit none
end subroutine micro_user

subroutine surf_user
  implicit none
end subroutine surf_user

subroutine idealisedfluxes
!    thls is not nudged
!    ps is not given
use modsurfdata, only: wqsurf,wtsurf
use modglobal, only : timee,xtime, pi
use modfields, only : rhobh ! use base state since transports are considered

real :: xfact

xfact=cos(pi/2.*(12.75-(xtime+timee/3600.))/5.25)
if(xfact.le.0.) xfact=0.
wqsurf =(554.*xfact**1.3)/(rlv*rhobh(1))
wtsurf =(270.*xfact**1.5)/(cpd*rhobh(1))

end subroutine idealisedfluxes

end module moduser
