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
subroutine initsurf_user
  implicit none
end subroutine initsurf_user

subroutine force_user
  implicit none
end subroutine force_user

subroutine force_testbed
  use modglobal, only : i1,j1,k1,kmax
  use modfields, only : qt0,thlp,qtp,up,vp, &
                        dqtdtls,dthldtls,dudtls,dvdtls   
  implicit none

  integer i,j,k
  !real dqt    (i1,j1,k1)
  !real idqtmax(i1,j1)


  !dqt    (:,:,:) = 0.
  !idqtmax(:,:)   = 0


  !--- search for level with greatest qt jump ---
  !do k=1,k1-1
  !  dqt(:,:,k) = qt0 (:,:,k) - qt0 (:,:,k+1)
  !enddo
  !idqtmax = maxloc( dqt, 3)   !find the index with the maximum jump

  !do k=1,kmax
  !  print *, "   qt-profile (2,2): ", k, qt0(2,2,k), dqt(2,2,k)
  !end do
  !print *, "         maxloc (2,2):", idqtmax(2,2)


  !--- add LS tendencies to budget ---
  do k=1,kmax
    do j=2,j1
      do i=2,i1

        thlp(i,j,k) = thlp(i,j,k) + dthldtls(k)
        qtp (i,j,k) = qtp (i,j,k) + dqtdtls(k)
        up  (i,j,k) = up  (i,j,k) + dudtls(k)
        vp  (i,j,k) = vp  (i,j,k) + dvdtls(k)

      enddo
    enddo
  enddo

  return
end subroutine force_testbed

subroutine rad_user
  implicit none
end subroutine rad_user

subroutine micro_user
  implicit none
end subroutine micro_user

subroutine surf_user
  implicit none
end subroutine surf_user

end module moduser
