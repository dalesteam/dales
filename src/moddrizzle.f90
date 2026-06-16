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
! Copyright 1993-2026 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!

! CJ: To me, this looks like it's doing exactly the same as sedimentation_cloud
! in the bulk microphysics. If that is the case, we should get rid of this one.

!> "Drizzle" microphysics.
module moddrizzle

  use modglobal,         only: i1, j1, kmax, rlv, cp, dzf, pi, rhow
  use modfields,         only: qtp, ql0, thlp, rhof, exnf
  use modmicrodata,      only: Nc_0, sig_g
  use modbulkmicro_data, only: c_st
  use modprecision,      only: field_r

contains
  !-----------------------------------------------------------------|
  !                                                                 |
  !      Hans Cuijpers   I.M.A.U.  23 May 1995                      |
  !                                                                 |
  !     purpose.                                                    |
  !     --------                                                    |
  !                                                                 |
  !      Calculates gravitational settling (or rainfall rate)       |
  !                                                                 |
  !**   interface.                                                  |
  !     ----------                                                  |
  !                                                                 |
  !     *drizzle* is called from *program*.                         |
  !                                                                 |
  !-----------------------------------------------------------------|
  subroutine drizzle

    integer :: &
      i, j, k
    real(field_r) :: &
      sedc, &
      csed

    csed = c_St*(3./(4.*pi*rhow))**(2./3.)*exp(5.*log(sig_g)**2.)
    sedc = 0.

    do k=1,kmax
      do j=2,j1
        do i=2,i1
          if (ql0(i,j,k)>0.0) then
            sedc= csed*((ql0(i,j,k+1)*rhof(k+1))**(5./3.)-(ql0(i,j,k)*rhof(k))**(5./3.))/(dzf(k)*rhof(k))
            qtp(i,j,k) = qtp(i,j,k) + sedc
            thlp(i,j,k) = thlp(i,j,k) - (rlv/(cp*exnf(k)))*sedc
          end if
        end do
      end do
    end do

  end subroutine drizzle

end module moddrizzle
