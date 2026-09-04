!> Common data for the spraying module.
!!
!! @author Stephan de Roode
!! @author Annelot Broerze
module modsprayingdata

  use modprecision, only: field_r

  implicit none

  public

  logical :: lwater_spraying = .false. !< Switch to enable water and sea salt spraying
  logical :: lsalt_spraying  = .false. !< Switch to enable sea salt spraying
  logical :: lsalt_sponge    = .false. !< Switch to enable nudging of salt to 0 at the boundary
  logical :: lcoupled        = .false. !< Enable coupling to aerosol module

  integer :: i_glob_spray = 2 !< Global i index of spraying point.
  integer :: j_glob_spray = 2 !< Global j index of spraying point.
  integer :: k_glob_spray = 2 !< Global k index of spraying point.

  integer :: i_spray = -999 !< Local i index of spraying point
  integer :: j_spray = -999 !< Local j index of spraying point
  integer :: k_spray = -999 !< Local k index of spraying point

  real(field_r) :: water_spray_rate = 1.    !< Water spray rate [kg/s]
  real(field_r) :: salt_spray_rate  = 0.030 !< Salt spray rate [kg/s]
  real(field_r) :: salinity = 0.03          !< Salinity of sprayed water [kg of salt per kg of seawater]

  character(len=20) :: tracer = "salt"     !< Name of the sprayed scalar (only used if lcoupled is false)
  integer           :: isv_salt = -1       !< Tracer index for salt mass concentration
  integer           :: isv_salt_n = -1     !< Tracer index for salt number concentration (only used if lcoupled is true)
  character(len=3)  :: target_mode = 'acs' !< Aerosol mode to spray in (acs or cos), only used if lcoupled is true

  logical :: my_process_sprays = .false. !< Whether this process should apply spraying (i.e. whether the spraying point is located on this process)

end module modsprayingdata
