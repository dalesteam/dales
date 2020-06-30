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

    logical :: sw_lsm   ! On/off switch LSM

    ! Soil grid
    integer :: kmax_soil
    real :: z_size_soil
    real, allocatable :: z_soil(:), zh_soil(:)
    real, allocatable :: dz_soil(:), dzh_soil(:)
    real, allocatable :: dzi_soil(:), dzhi_soil(:)

    ! Soil properties
    integer, allocatable :: soil_index(:,:,:)

    ! Data structure for sub-grid tiles
    type lsm_tile
        ! Static properties:
        real, allocatable :: z0m(:,:), z0h(:,:)
        ! Dynamic tile fraction:
        real, allocatable :: frac(:,:)
        ! Monin-obukhov / surface layer:
        real, allocatable :: obuk(:,:), ustar(:,:), ra(:,:)
        ! Conductivity skin layer:
        real, allocatable :: lambda_stable(:,:), lambda_unstable(:,:)
        ! Surface fluxes:
        real, allocatable :: H(:,:), LE(:,:), G(:,:), wthl(:,:), wqt(:,:)
        ! Surface temperature and humidity:
        real, allocatable :: tskin(:,:), qskin(:,:)
        ! Vegetation properties:
        real, allocatable :: lai(:,:), rs_min(:,:)
    end type lsm_tile

    ! Tiles for low veg (lv), high veg (hv), bare soil (bs), wet skin (ws):
    type(lsm_tile) tile_lv, tile_hv, tile_bs, tile_ws

    ! Land-surface / van Genuchten parameters from NetCDF input table.
    real, allocatable :: &
        theta_res(:), theta_wp(:), theta_fc(:), theta_sat(:), &
        gamma_sat(:), vg_a(:), vg_l(:), vg_n(:)

contains

subroutine lsm
    use modsurfdata, only : thlflux, qtflux
    implicit none

    thlflux(:,:) = 0.
    qtflux(:,:) = 0.

end subroutine lsm

!
! Initialise the land-surface model
!
subroutine initlsm
    use modglobal,   only : ifnamopt, fname_options, checknamelisterror
    use modmpi,      only : myid, comm3d, mpierr, mpi_logical, my_real
    use modsurfdata, only : isurf, zsoil
    implicit none

    integer :: ierr
    logical :: sw_homogeneous

    ! Read namelist
    namelist /NAMLSM/ &
        sw_lsm, sw_homogeneous, zsoil, z_size_soil

    allocate(zsoil(kmax_soil))

    if (myid == 0) then
        open(ifnamopt, file=fname_options, status='old', iostat=ierr)
        read(ifnamopt, NAMLSM, iostat=ierr)
        call checknamelisterror(ierr, ifnamopt, 'NAMLSM')
        write(6, NAMLSM)
        close(ifnamopt)
    end if

    ! Broadcast namelist values to all MPI tasks
    call MPI_BCAST(sw_lsm, 1, mpi_logical, 0, comm3d, mpierr)
    call MPI_BCAST(sw_homogeneous, 1, mpi_logical, 0, comm3d, mpierr)

    call MPI_BCAST(z_soil, kmax_soil, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(z_size_soil, 1, my_real, 0, comm3d, mpierr)

    if (sw_lsm) then
        ! Create/calculate soil grid properties
        call create_soil_grid

        ! Allocate required fields in modsurfacedata, and arrays / tiles from this module:
        call allocate_fields

        if (sw_homogeneous) then
            ! Initialise homogeneous LSM from namelist input
            call init_homogeneous
        else
            ! Initialise heterogeneous LSM from external input
            stop 'Heterogeneous LSM not (yet) implemented!'
        end if

        ! Read the soil parameter table
        call read_soil_table
    end if

end subroutine initlsm

!
! Cleanup (deallocate) the land-surface model
!
subroutine exitlsm
    implicit none

    deallocate( theta_res, theta_wp, theta_fc, theta_sat, gamma_sat, vg_a, vg_l, vg_n )
end subroutine exitlsm

!
! Calculate soil grid properties
!
subroutine create_soil_grid
    implicit none
    integer :: k

    allocate(dz_soil  (kmax_soil  ))
    allocate(dzi_soil (kmax_soil  ))
    allocate(zh_soil  (kmax_soil+1))
    allocate(dzh_soil (kmax_soil+1))
    allocate(dzhi_soil(kmax_soil+1))

    ! Half level heights
    zh_soil(1) = z_size_soil
    zh_soil(kmax_soil+1) = 0.

    do k=2, kmax_soil
        zh_soil(k) = 0.5*(z_soil(k) + z_soil(k-1))
    end do

    ! Grid spacing full and half levels
    do k=1, kmax_soil
        dz_soil(k) = zh_soil(k+1) - zh_soil(k)
    end do

    do k=2, kmax_soil
        dzh_soil(k) = z_soil(k) - z_soil(k-1)
    end do

    dzh_soil(1) = 2*(z_soil(1) - zh_soil(1))
    dzh_soil(kmax_soil+1) = 2*(-z_soil(kmax_soil))

    ! Inverse grid spacings
    dzi_soil(:) = 1./dz_soil(:)
    dzhi_soil(:) = 1./dzh_soil(:)

end subroutine create_soil_grid

!
! Allocate all LSM fields
!
subroutine allocate_fields
    use modglobal, only : i2, j2
    use modsurfdata, only : &
        tsoil, phiw
    implicit none

    ! Allocate soil variables
    allocate(soil_index(i2, j2, kmax_soil))
    allocate(tsoil     (i2, j2, kmax_soil))
    allocate(phiw      (i2, j2, kmax_soil))

    ! Allocate the tiled variables
    call allocate_tile(tile_lv)
    call allocate_tile(tile_hv)
    call allocate_tile(tile_bs)
    call allocate_tile(tile_ws)

end subroutine allocate_fields

!
! Allocate all fields of a LSM tile
!
subroutine allocate_tile(tile)
    use modglobal, only : i2, j2
    implicit none
    type(lsm_tile), intent(inout) :: tile

    ! Static properties:
    allocate(tile % z0m(i2, j2))
    allocate(tile % z0h(i2, j2))

    ! Dynamic tile fraction:
    allocate(tile % frac(i2, j2))

    ! Monin-obukhov / surface layer:
    allocate(tile % obuk(i2, j2))
    allocate(tile % ustar(i2, j2))
    allocate(tile % ra(i2, j2))

    ! Conductivity skin layer:
    allocate(tile % lambda_stable(i2, j2))
    allocate(tile % lambda_unstable(i2, j2))

    ! Surface fluxes:
    allocate(tile % H(i2, j2))
    allocate(tile % LE(i2, j2))
    allocate(tile % G(i2, j2))
    allocate(tile % wthl(i2, j2))
    allocate(tile % wqt(i2, j2))

    ! Surface temperature and humidity:
    allocate(tile % tskin(i2, j2))
    allocate(tile % qskin(i2, j2))

    ! Vegetation properties:
    allocate(tile % lai(i2, j2))
    allocate(tile % rs_min(i2, j2))

end subroutine allocate_tile

!
! Initialise the LSM homogeneous from namelist input
!
subroutine init_homogeneous
    use modglobal,   only : ifnamopt, fname_options, checknamelisterror
    use modmpi,      only : myid, comm3d, mpierr, mpi_logical, my_real
    use modsurfdata, only : tsoil, phiw
    implicit none

    integer :: ierr, k
    real :: z0m_low, z0m_high, z0m_bare
    real :: z0h_low, z0h_high, z0h_bare
    real :: lambda_s_low, lambda_s_high, lambda_s_bare
    real :: lambda_us_low, lambda_us_high, lambda_us_bare
    real :: lai_low, lai_high
    real :: rs_min_low, rs_min_high

    real, allocatable :: t_soil_p(:), theta_soil_p(:)

    ! Read namelist
    namelist /NAMLSM_HOMOGENEOUS/ &
        z0m_low, z0m_high, z0m_bare, &
        z0h_low, z0h_high, z0h_bare, &
        lambda_s_low, lambda_s_high, lambda_s_bare, &
        lambda_us_low, lambda_us_high, lambda_us_bare, &
        lai_low, lai_high, &
        rs_min_low, rs_min_high, &
        t_soil_p, theta_soil_p

    allocate(t_soil_p(kmax_soil), theta_soil_p(kmax_soil))

    if (myid == 0) then
        open(ifnamopt, file=fname_options, status='old', iostat=ierr)
        read(ifnamopt, NAMLSM_HOMOGENEOUS, iostat=ierr)
        call checknamelisterror(ierr, ifnamopt, 'NAMLSM_HOMOGENEOUS')
        write(6, NAMLSM_HOMOGENEOUS)
        close(ifnamopt)
    end if

    ! Broadcast to all MPI tasks
    call MPI_BCAST(z0m_low,  1, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(z0m_high, 1, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(z0m_bare, 1, my_real, 0, comm3d, mpierr)

    call MPI_BCAST(z0h_low,  1, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(z0h_high, 1, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(z0h_bare, 1, my_real, 0, comm3d, mpierr)

    call MPI_BCAST(lambda_s_low,  1, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(lambda_s_high, 1, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(lambda_s_bare, 1, my_real, 0, comm3d, mpierr)

    call MPI_BCAST(lambda_us_low,  1, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(lambda_us_high, 1, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(lambda_us_bare, 1, my_real, 0, comm3d, mpierr)

    call MPI_BCAST(lai_low,  1, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(lai_high, 1, my_real, 0, comm3d, mpierr)

    call MPI_BCAST(rs_min_low,  1, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(rs_min_high, 1, my_real, 0, comm3d, mpierr)

    call MPI_BCAST(t_soil_p,     kmax_soil, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(theta_soil_p, kmax_soil, my_real, 0, comm3d, mpierr)

    ! Set values
    tile_lv % z0m(:,:) = z0m_low
    tile_hv % z0m(:,:) = z0m_high
    tile_bs % z0m(:,:) = z0m_bare

    tile_lv % z0h(:,:) = z0h_low
    tile_hv % z0h(:,:) = z0h_high
    tile_bs % z0h(:,:) = z0h_bare

    tile_lv % lambda_stable(:,:) = lambda_s_low
    tile_hv % lambda_stable(:,:) = lambda_s_high
    tile_bs % lambda_stable(:,:) = lambda_s_bare

    tile_lv % lambda_unstable(:,:) = lambda_us_low
    tile_hv % lambda_unstable(:,:) = lambda_us_high
    tile_bs % lambda_unstable(:,:) = lambda_us_bare

    tile_lv % lai(:,:) = lai_low
    tile_hv % lai(:,:) = lai_high

    tile_lv % rs_min(:,:) = rs_min_low
    tile_hv % rs_min(:,:) = rs_min_high

    do k=1, kmax_soil
        tsoil(:,:,k) = t_soil_p(k)
        phiw (:,:,k) = theta_soil_p(k)
    end do

    ! Cleanup!
    deallocate(t_soil_p, theta_soil_p)

end subroutine init_homogeneous

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
