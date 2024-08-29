module modmicrodata_aerosol

  use modprecision, only: field_r
  use modmode_type, only: mode_t
  use modglobal,    only: i1, j1, k1, ih, jh

  implicit none

  integer,      parameter :: maxmodes = 9 !< M7 + in-rain + in-cloud
  character(3), parameter :: modenames(9) = (/ 'nus', 'ais', 'acs', 'cos', 'aii', 'aci', 'coi', 'inr', 'inc'/)
  integer,      parameter :: iNUS = 1, iAIS = 2, iACS = 3, iCOS = 4, iAII = 5, iACI = 6, &
    &                        iCOI = 7, iINR = 8, iINC = 9

  type(mode_t) :: modes(maxmodes)

  integer :: nmodes                    !< Number of modes in use
  logical :: lkohler = .true.
  real(field_r) :: ssat

  integer, parameter :: ncld = 10, &
    &                   nrai = 5, &
    &                   naer_blc = 100, &
    &                   naer_inc = 60

  real(field_r), parameter :: ncmin = 1.0e3_field_r

  real(field_r), allocatable :: Nc(:,:,:)
  real(field_r), allocatable :: sed_qr_dup(:,:,:)

  real(field_r), dimension(naer_inc) :: aerrad,     logaerrad
  real(field_r), dimension(naer_blc) :: aerrad_blc, logaerrad_blc
  real(field_r), dimension(ncld)     :: cldrad,     logcldrad
  real(field_r), dimension(nrai)     :: rainrate,   lograinrate

  real(field_r), dimension(ncld,naer_inc) :: incmass, incnumb
  real(field_r), dimension(nrai,naer_blc) :: blcmass, blcnumb

end module modmicrodata_aerosol
