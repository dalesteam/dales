module modmicrodata_aerosol

  use modprecision, only: field_r
  use modmode_type, only: mode_t

  implicit none

  integer, parameter :: maxmodes = 9 !< M7 + in-rain + in-cloud

  type(mode_t) :: modes(maxmodes)

  integer :: nmodes                    !< Number of modes in use
  logical :: lkohler = .true.

end module modmicrodata_aerosol
