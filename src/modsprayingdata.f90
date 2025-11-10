!> \file modspraying.f90
!! Stephan de Roode and Annelot Broerze

module modsprayingdata

  implicit none
  save

  logical :: lwater_spraying   = .false.     !< Switch to enable water and sea salt spraying
  logical :: lsalt_spraying    = .false.     !< Switch to enable sea salt spraying
  logical :: lsalt_sponge      = .false.     !< Switch to enable nudging of salt to 0 at the boundary
  !< Set default value for water spraying
  !< Set default location

  integer :: i_glob_spray = 2  !these are global grid points (numbering in the whole domain)
  integer :: j_glob_spray = 2
  integer :: k_glob_spray = 2

  integer :: i_loc_spray = -999
  integer :: j_loc_spray = -999
  integer :: k_loc_spray = -999

  real :: water_spray_rate = 1.    ! kg/sec water spraying excluding salt
  real :: salt_spray_rate  = 0.030 ! kg/sec salt spraying

  real :: dqldt_spraying = 0.   ! convert water_spray_rate to local value in LES grid
  real :: dsvdt_spraying = 0.   ! convert salt spray rate to local value in LES grid

  character(20) :: tracer = "salt" ! name of the sprayed scalar
  integer :: isv_salt = -1

  real :: salinity = 0.03 ! this definition assumes 1 kg of sea water contains a mass of salt equal to salinity kg

end module modsprayingdata
