!
! Copyright (c) 2020-2020 Wageningen University and Research (WUR)
!
! This file is part of DALES
!
! DALES is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with DALES.  If not, see <http://www.gnu.org/licenses/>.
!

module modlsm
    use netcdf
    implicit none
    public :: initlsm, lsm, exitlsm

    ! Land-surface / van Genuchten parameters from NetCDF input table.
    real, allocatable :: &
        theta_res(:), theta_wp(:), theta_fc(:), theta_sat(:), gamma_sat(:), vg_a(:), vg_l(:), vg_n(:)
contains

!
! Initialise the land-surface model
!
subroutine initlsm
    implicit none

    ! Read the soil parameter table
    call read_soil_table
end subroutine initlsm

subroutine lsm
    implicit none
end subroutine lsm

!
! Cleanup (deallocate) the land-surface model
!
subroutine exitlsm
    implicit none

    deallocate( theta_res, theta_wp, theta_fc, theta_sat, gamma_sat, vg_a, vg_l, vg_n )
end subroutine exitlsm

!
! Read the input table with the (van Genuchten) soil parameters
!
subroutine read_soil_table
    implicit none
    integer :: n, ncid, dimid, varid

    ! Open the NetCDF file and read the table size
    print*,'Reading "van_genuchten_parameters.nc"'
    call check( nf90_open('van_genuchten_parameters.nc', nf90_nowrite, ncid) )
    call check( nf90_inq_dimid(ncid, 'index', dimid) )
    call check( nf90_inquire_dimension(ncid, dimid, len=n) )

    ! Allocate variables
    allocate( &
        theta_res(n), theta_wp(n), theta_fc(n), theta_sat(n), &
        gamma_sat(n), vg_a(n), vg_l(n), vg_n(n) )

    ! Read variables
    call check( nf90_inq_varid(ncid, 'theta_res', varid) )
    call check( nf90_get_var(ncid, varid, theta_res) )

    call check( nf90_inq_varid(ncid, 'theta_wp', varid) )
    call check( nf90_get_var(ncid, varid, theta_wp) )

    call check( nf90_inq_varid(ncid, 'theta_fc', varid) )
    call check( nf90_get_var(ncid, varid, theta_fc) )

    call check( nf90_inq_varid(ncid, 'theta_sat', varid) )
    call check( nf90_get_var(ncid, varid, theta_sat) )

    call check( nf90_inq_varid(ncid, 'gamma_sat', varid) )
    call check( nf90_get_var(ncid, varid, gamma_sat) )

    call check( nf90_inq_varid(ncid, 'alpha', varid) )
    call check( nf90_get_var(ncid, varid, vg_a) )

    call check( nf90_inq_varid(ncid, 'l', varid) )
    call check( nf90_get_var(ncid, varid, vg_l) )

    call check( nf90_inq_varid(ncid, 'n', varid) )
    call check( nf90_get_var(ncid, varid, vg_n) )

    call check( nf90_close(ncid) )

end subroutine read_soil_table

!
! Convert soil hydraulic head to soil water content, using van Genuchten parameterisation.
!
pure function psi_to_theta(theta_res, theta_sat, vg_a, vg_n, vg_m, psi) result(res)
    implicit none
    real, intent(in) :: theta_res, theta_sat, vg_a, vg_n, vg_m, psi
    real :: res

    res = theta_res + (theta_sat - theta_res) * (1. / (1.+ abs(vg_a * psi)**vg_n))**vg_m
end function psi_to_theta

!
! Convert soil water content to hydraulic head, using van Genuchten parameterisation.
!
pure function theta_to_psi(theta_res, theta_sat, vg_a, vg_n, vg_m, theta) result(res)
    implicit none
    real, intent(in) :: theta_res, theta_sat, vg_a, vg_n, vg_m, theta
    real :: res

    res = -(vg_a**(-vg_n) * (-1. + ((theta_res - theta_sat)/(theta_res - theta))**(1./vg_m)))**(1./vg_n)
end function theta_to_psi

!
! Calculate hydraulic diffusivity using van Genuchten parameterisation.
!
pure function calc_diffusivity_vg( &
        theta_norm, vg_a, vg_l, vg_m, lambda_sat, theta_sat, theta_res) result(res)
    implicit none
    real, intent(in) :: theta_norm, vg_a, vg_l, vg_m, lambda_sat, theta_sat, theta_res
    real :: res

    res = (1.-vg_m)*lambda_sat / (vg_a * vg_m * (theta_sat-theta_res)) * theta_norm**(vg_l-(1./vg_m)) * &
             (  (1.-theta_norm**(1./vg_m))**(-vg_m) + (1.-theta_norm**(1/vg_m))**vg_m - 2. )
end function calc_diffusivity_vg

!
! Calculate hydraulic conductivity using van Genuchten parameterisation.
!
pure function calc_conductivity_vg(theta_norm, vg_l, vg_m, lambda_sat) result(res)
    implicit none
    real, intent(in) :: theta_norm, vg_l, vg_m, lambda_sat
    real :: res

    res = lambda_sat * theta_norm**vg_l * ( 1.- (1.-theta_norm**(1./vg_m))**vg_m )**2.
end function calc_conductivity_vg

!
! Check NetCDF calls
!
subroutine check(status)
    integer, intent (in) :: status
    if(status /= nf90_noerr) then
        print *,'NetCDF error: ', trim(nf90_strerror(status))
        stop
    end if
end subroutine check

end module modlsm
