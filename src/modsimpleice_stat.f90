!> Statistics routines for single-moment ice microphysics scheme.
module modsimpleice_stat

  use modglobal,         only: i1, j1, k1, kmax, ih, jh
  use modmicrodata,      only: l_rain
  use modmicrodata,      only: epscloud, epsqr, epsprec
  use modstat_profiles,  only: add_profile, sample_field, do_procblock, &
                               is_sampling_timestep
  use modprecision,      only: field_r

  implicit none

  private

  public :: init_simpleice_stat
  public :: simpleice_stat

contains

  subroutine init_simpleice_stat()

    character(len=*), parameter :: routine = 'init_simpleice_stat'

    character(len=4) :: dim

    if (do_procblock()) then
      dim = 'tttt'
    else
      dim = 'tt'
    end if

    call add_profile('cfrac', 'Cloud fraction', '-', 'tt')
    call add_profile('rainrate', 'Echo rain rate', 'W/m2', dim)
    call add_profile('preccount', 'Precipitation flux area fraction', '-', dim)
    call add_profile('raincount', 'Rain water content area fraction', '-', dim)
    call add_profile('precmn', 'Rain rate', 'W/m2', dim)
    call add_profile('qrmn', 'Precipitation specific humidity', 'kg/kg', dim)
    call add_profile('qrpaccr', 'Accretion rain water content tendency', &
                     'kg/kg/s', dim)
    call add_profile('qrpauto', 'Autoconversion rain water content tendency', &
                     'kg/kg/s', dim)
    call add_profile('qrpsed', 'Sedimentation rain water content tendency', &
                     'kg/kg/s', dim)
    call add_profile('qrpevap', 'Evaporation rain water content tendency', &
                     'kg/kg/s', dim)
    call add_profile('qrpclip', 'Rain water content tendency due to clipping', &
                     'kg/kg/s', dim)
    call add_profile('qrptot', 'Total rain water content tendency', &
                     'kg/kg/s', dim)

  end subroutine init_simpleice_stat

  subroutine simpleice_stat(ql, qr, precip)

    real(field_r), intent(in) :: ql(2-ih:,2-jh:,:) !< Cloud liquid water content [-]
    real(field_r), intent(in) :: qr(2:,2:,:)       !< Rain water content [-]
    real(field_r), intent(in) :: precip(2:,2:,:)   !< Precipitation flux [kg/m2/s]

    real(field_r), allocatable :: is_cloud(:,:,:)
    real(field_r), allocatable :: is_rain(:,:,:)
    real(field_r), allocatable :: is_precip(:,:,:)

    integer :: i, j, k

    if (is_sampling_timestep()) then

      allocate(is_cloud(2:i1,2:j1,1:k1))
      allocate(is_rain(2:i1,2:j1,1:k1))
      allocate(is_precip(2:i1,2:j1,1:k1))

      !$acc data create(is_cloud, is_rain, is_precip)
!!$omp target data map(alloc:is_cloud,is_rain,is_precip)

      !$acc parallel loop collapse(3) default(present)
!!$omp target teams loop collapse(3) defaultmap(present:aggregate)&
!!$omp defaultmap(present:allocatable)
      do k = 1, kmax
        do j = 2, j1
          do i = 2, i1
            is_cloud(i,j,k) = merge(1, 0, ql(i,j,k) > epscloud)
            is_rain(i,j,k) = merge(1, 0, qr(i,j,k) > epsqr)
            is_precip(i,j,k) = merge(1, 0, precip(i,j,k) > epsprec)
          end do
        end do
      end do

      call sample_field('cfrac', is_cloud)
      call sample_field('raincount', is_rain)
      call sample_field('preccount', is_precip)
      call sample_field('qrmn', qr)
      call sample_field('precmn', precip)

      !$acc end data
!!$omp end target data

      deallocate(is_cloud, is_rain, is_precip)

    end if

  end subroutine simpleice_stat

end module modsimpleice_stat
