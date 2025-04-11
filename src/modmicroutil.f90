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
! Copyright 1993-2025, The DALES Team.
!

!> Common utility routines for microphysics.
module modmicroutil

  use modmicrodata, only: l_mur_cst, mur_cst, pirhow
  use modprecision, only: field_r

  implicit none

  public

contains
  !> Compute Mu parameter of the Gamma distribution.
  !! \param qr Rain water mixing ratio.
  !! \param rho Air density.
  !! \return Mu.
  function calc_mur(qr, rho) result(mu)

    real(field_r), intent(in) :: qr, rho

    real(field_r) :: mu
    !$acc routine seq

    if (.not. l_mur_cst) then
      mu = min(30.0_field_r, -1 + 0.008_field_r / (qr * rho**0.6_field_r))
    else
      mu = mur_cst
    end if

  end function calc_mur

  !> Compute mean droplet mass
  !! \param rho Air density.
  !! \param qr Rain water mixing ratio.
  !! \param nr Droplet number concentration.
  !! \param xrmin Lower bound of mass.
  !! \param xrmax Upper bound of mass.
  !! \return Droplet mass.
  function calc_xr(rho, qr, nr, xrmin, xrmax) result(xr)

    real(field_r), intent(in) :: rho, qr, nr, xrmin, xrmax

    real(field_r) :: xr
    !$acc routine seq

    xr = rho * qr / nr
    xr = min(max(xr, xrmin), xrmax)

  end function calc_xr

  !> Compute mean droplet diameter.
  !! \param xr Mass of droplet.
  !! \return Droplet diameter.
  function calc_dvr(xr) result(d)

    real(field_r), intent(in) :: xr

    real(field_r) :: d
    !$acc routine seq

    d = (xr / pirhow)**(1._field_r/3)

  end function calc_dvr

  !> Compute lambda parameter of the Gamma distribution.
  !! \param mur Mu parameter of the Gamma distribution.
  !! \param dr Mean droplet diameter.
  !! \return Lambda.
  function calc_lbdr(mur, dr) result(lambda)

    real(field_r), intent(in) :: mur, dr

    real(field_r) :: lambda
    !$acc routine seq

    lambda = ((mur + 3) * (mur + 2) * (mur + 1))**(1.0_field_r/3) / dr

  end function calc_lbdr

end module modmicroutil