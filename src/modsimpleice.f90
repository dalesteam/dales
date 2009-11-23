!> \file modsimpleice.f90

!>
!!  Ice microphysics.
!>
!! Calculates ice microphysics in a cheap scheme without prognostic CDC
!  simpleice is called from *modmicrophysics*
!! \see  Grabowski, 1998, JAS
!!  \author Steef B\"oing, TU Delft
!!  \par Revision list
!! \todo test
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


module modsimpleice
  use modmicrodata

  implicit none
  save
  contains

!> Initializes and allocates the arrays
  subroutine initsimpleice
    use modglobal, only : ih,i1,jh,j1,k1
    use modmpi,    only : myid
    implicit none

    allocate (qr(2-ih:i1+ih,2-jh:j1+jh,k1)        & ! qr is converted from a scalar variable
             ,qrp(2-ih:i1+ih,2-jh:j1+jh,k1)       & ! qr tendency due to microphysics only, for statistics
             ,thlpmcr(2-ih:i1+ih,2-jh:j1+jh,k1)   & ! thl tendency due to microphysics only, for statistics
             ,qtpmcr(2-ih:i1+ih,2-jh:j1+jh,k1)    & ! qt tendency due to microphysics only, for statistics
             ,sed_qr(2-ih:i1+ih,2-jh:j1+jh,k1)    & ! sedimentation rain droplets mixing ratio
             ,qr_spl(2-ih:i1+ih,2-jh:j1+jh,k1)    & ! time-splitting substep qr
             ,ilratio(2-ih:i1+ih,2-jh:j1+jh,k1)     & ! partition ratio liquid/solid
             ,lambdal(2-ih:i1+ih,2-jh:j1+jh,k1)   & ! slope parameter for liquid phase
             ,lambdai(2-ih:i1+ih,2-jh:j1+jh,k1))    ! slope parameter for ice phase

    allocate (Nr(2-ih:i1+ih,2-jh:j1+jh,k1)        & ! CDC, only for statistics
             ,qc(2-ih:i1+ih,2-jh:j1+jh,k1))         ! = ql, also for statistics

    allocate (qrmask(2-ih:i1+ih,2-jh:j1+jh,k1)    & ! mask for rain water
             ,qcmask(2-ih:i1+ih,2-jh:j1+jh,k1))     ! mask for cloud water

  end subroutine initsimpleice

!> Cleaning up after the run
  subroutine exitsimpleice
    implicit none
  end subroutine exitsimpleice

!> Calculates the microphysical source term.
  subroutine simpleice
    implicit none

    if (l_rain) then
!     call bulkmicrotend
      call autoconvert
!     call bulkmicrotend
      call accrete
!     call bulkmicrotend
      call evaposite
!     call bulkmicrotend
      call precipitate
!     call bulkmicrotend
    endif

  end subroutine simpleice

  subroutine autoconvert
    implicit none
  end subroutine autoconvert

  subroutine accrete
    implicit none
  end subroutine accrete

  subroutine precipitate
    implicit none
  end subroutine precipitate

  subroutine evaposite
    implicit none
  end subroutine evaposite

end module modsimpleice


