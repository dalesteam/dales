!> Common data and routines for aerosol microphysics.
module modaerosol_common

  use modprecision, only: field_r

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

end module modaerosol_common