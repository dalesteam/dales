!> Common data and routines for aerosol microphysics.
module modaerosol_common

  use modprecision, only: field_r
  use modglobal,    only: pi

  implicit none

  public

  integer, parameter :: &
    maxspecies = 5,     &
    maxmodes = 9

  character(len=*), parameter :: &
    aerosol_names(maxspecies) = [character(len=3) :: &
      'so4', 'ss', 'pom', 'bc', 'du'],               &
    aerosol_stdnames(maxspecies) = [character(len=26) :: &
      'sulfate',                                         &
      'sea_salt',                                        &
      'particulate_organic_matter',                      &
      'black_carbon',                                    &
      'dust'],                                           &
    mode_names(maxmodes) = &
      ['nus', 'ais', 'acs', 'cos', 'aii', 'aci', 'coi', 'inc', 'inr'], &
    mode_longnames(maxmodes) = [character(len=27) :: &
      'soluble_nucleation_mode',     &
      'soluble_Aitken_mode',         &
      'soluble_accumulation_mode',   &
      'soluble_coarse_mode',         &
      'insoluble_Aitken_mode',       &
      'insoluble_accumulation_mode', &
      'insoluble_coarse_mode',       &
      'in_cloud_mode',               &
      'in_rain_mode']

  real(field_r), parameter :: &
    aerosol_densities(maxspecies) = [1841, 2165, 1800, 1300, 2560], &
    sigma_g_modes(maxmodes) = [1.59, 1.59, 1.59, 2.0, 1.59, 1.59, 2.0, 1.5, 1.5]
  
  ! Indices of modes in mode list.
  integer, parameter :: &
    iNUS = 1, & ! Nucleation mode.
    iAIS = 2, & ! Aitken soluble mode.
    iACS = 3, & ! Accumulation soluble mode.
    iCOS = 4, & ! Coarse soluble mode.
    iAII = 5, & ! Aitken insoluble mode.
    iACI = 6, & ! Accumulation insoluble mode.
    iCOI = 7, & ! Coarse insoluble mode.
    iINC = 8, & ! In-cloud mode.
    iINR = 9    ! In-rain mode.

  ! Aerosol identifiers.
  integer, parameter :: &
    iSO4 = 1, & ! Sulfate.
    iSS = 2,  & ! Sea salt.
    iPOM = 3, & ! Primary organic matter.
    iBC = 4,  & ! Black carbon.
    iDU = 5     ! Dust.

contains

  !> Compute the median diameter of a log-normally distributed mode.
!NVF$ INLINE
  pure function calc_median_diameter(n, q, rho, sig_g) result(dm)

    real(field_r), intent(in) :: n      !< Aerosol number concentration [m-3].
    real(field_r), intent(in) :: q(:)   !< Aerosol mass concentrations [kg kg-1].
    real(field_r), intent(in) :: rho(:) !< Aerosol densities [kg m-3].
    real(field_r), intent(in) :: sig_g  !< Geometric standard deviation.

    real(field_r) :: m     !< Total aerosol mass concentration [kg kg-1].
    real(field_r) :: rho_m !< Mean aerosol density [kg m-3].
    real(field_r) :: dm    !< Median diameter [m].

    integer :: s !< Loop index

    m = 0
    rho_m = 0

    do s = 1, size(q)
      if (rho(s) > 0.0_field_r) then
        m = m + q(s)
        rho_m = rho_m + q(s) / rho(s)
      end if
    end do

    m = max(0.0_field_r, m)
    rho_m = max(0.0_field_r, m / (rho_m + 1E-16))

    dm = ((6 * m) / (pi * n * rho_m + 1E-16))**(1.0_field_r / 3) &
         * exp(- 0.5_field_r * 3 * log(sig_g) * log(sig_g))

    dm = max(0.0_field_r, dm)

  end function calc_median_diameter

end module modaerosol_common