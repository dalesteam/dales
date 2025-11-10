
module modbulkmicro_stat

  use bulkmicro_sb,      only: xrmin_sb => xrmin, xrmax_sb => xrmax
  use bulkmicro_kk,      only: xrmin_kk => xrmin, xrmax_kk => xrmax
  use modfields,         only: ql0, rhof
  use modglobal,         only: i1, j1, k1, pirhow
  use modmicrodata,      only: precep, l_rain
  use modbulkmicro_data, only: qr, nr, epscloud, epsqr, epsprec, l_sb, &
                               l_sedc, l_mur_cst, mur_cst
  use modstat_profiles,  only: add_profile, sample_field, is_sampling_timestep
  use modprecision,      only: field_r

  implicit none

  private

  public :: init_bulkmicro_stat
  public :: bulkmicro_stat

contains

  include 'microphysics.inc'

  subroutine init_bulkmicro_stat

    call add_profile('cfrac', 'Cloud fraction', '-', 'tt')

    if (l_rain) then
      call add_profile('rainrate', 'Echo rain rate', 'W/m^2', 'tt')
      call add_profile('preccount', 'Precipitation flux area fraction', '-', 'tt')
      call add_profile('nrrain', 'Rain droplet number concentration', '#/m3', 'tt')
      call add_profile('raincount', 'Rain water content area fraction', '-', 'tt')
      call add_profile('precmn', 'Rain rate', 'W/m^2', 'tt')
      call add_profile('dvrmn', 'Precipitation mean diameter', 'm', 'tt')
      call add_profile('qrmn', 'Precipitation specific humidity', 'kg/kg', 'tt')
      call add_profile('npauto', 'Autoconversion rain drop tendency', '#/m3/s', 'tt')
      call add_profile('npaccr', 'Accretion rain drop tendency', '#/m3/s', 'tt')
      call add_profile('npsed', 'Sedimentation rain drop tendency', '#/m3/s', 'tt')
      call add_profile('npevap', 'Evaporation rain drop tendency', '#/m3/s', 'tt')
      call add_profile('npclip', 'Rain drop tendency due to clipping', '#/m3/s', 'tt')
      call add_profile('nptot', 'Total rain drop tendency', '#/m3/s',  'tt')
      call add_profile('qrpauto', 'Autoconversion rain water content tendency', 'kg/kg/s', 'tt')
      call add_profile('qrpaccr', 'Accretion rain water content tendency', 'kg/kg/s', 'tt')
      call add_profile('qrpsed', 'Sedimentation rain water content tendency', 'kg/kg/s', 'tt')
      call add_profile('qrpevap', 'Evaporation rain water content tendency', 'kg/kg/s', 'tt')
      call add_profile('qrpclip', 'Rain water content tendency due to clipping', 'kg/kg/s', 'tt')
      call add_profile('qrptot', 'Total rain water content tendency', 'kg/kg/s', 'tt')
    end if

    if (l_sedc) then
      call add_profile('qtpsedc', 'Sedimentation total water content tendency', 'kg/kg/s', 'tt')
    end if

  end subroutine init_bulkmicro_stat

  subroutine bulkmicro_stat

    integer       :: i, j, k
    real(field_r) :: xr, xrmin, xrmax

    real(field_r), allocatable :: is_cloud(:,:,:)
    real(field_r), allocatable :: is_rain(:,:,:)
    real(field_r), allocatable :: is_precip(:,:,:)
    real(field_r), allocatable :: dvr(:,:,:)

    if (is_sampling_timestep()) then

      ! TODO: implement some temp field allocator
      allocate(is_cloud(2:i1,2:j1,1:k1))
      allocate(is_rain(2:i1,2:j1,1:k1))
      allocate(is_precip(2:i1,2:j1,1:k1))
      allocate(dvr(2:i1,2:j1,1:k1))

      !$acc data create(is_cloud, is_rain, is_precip, dvr)

      if (l_sb) then
        xrmin = xrmin_sb
        xrmax = xrmax_sb
      else
        xrmin = xrmin_kk
        xrmax = xrmax_kk
      end if

      !$acc parallel loop collapse(3) default(present) private(xr)
      do k = 1, k1
        do j = 2, j1
          do i = 2, i1
            is_cloud(i,j,k) = merge(1, 0, ql0(i,j,k) > epscloud)
            is_rain(i,j,k) = merge(1, 0, qr(i,j,k) > epsqr)
            is_precip(i,j,k) = merge(1, 0, precep(i,j,k) > epsprec)
            xr = calc_xr(rhof(k), qr(i,j,k), nr(i,j,k), xrmin, xrmax)
            dvr(i,j,k) = merge(calc_dvr(xr), 0.0_field_r, qr(i,j,k) > epsqr)
          end do
        end do
      end do

      call sample_field('cfrac', is_cloud)
      call sample_field('raincount', is_rain)
      call sample_field('preccount', is_precip)
      call sample_field('dvrmn', dvr)
      call sample_field('nrrain', nr)
      call sample_field('qrmn', qr)
      call sample_field('precmn', precep)

      !$acc end data

      deallocate(is_cloud, is_rain, is_precip, dvr)

    end if

  end subroutine bulkmicro_stat

end module modbulkmicro_stat