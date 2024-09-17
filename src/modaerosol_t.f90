module modaerosol_t
  use modtracer_type, only: T_tracer

  implicit none

  private

  type, public :: aerosol_t
    character(32) :: name
    integer       :: mode_idx(9)
  end type aerosol_t

contains

end module modaerosol_t