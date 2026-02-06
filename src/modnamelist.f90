!> Read namelists and check settings.
module modnamelist

  use modaerosol,        only: aerosol_read_namelist
  use modchecksim,       only: checksim_read_namelist
  use modlatsponge,      only: lateral_sponge_read_namelist
  use modmicrophysics,   only: microphysics_read_namelist
  use modpois,           only: poisson_solver_read_namelist
  use modsurface,        only: surface_read_namelist
  use modthermodynamics, only: thermodynamics_read_namelist

  implicit none

  private

  public :: read_namelists

contains

  !> Read all namelists.
  subroutine read_namelists(nml_filename)

    character(len=*), intent(in) :: nml_filename

    ! Core modules
    call thermodynamics_read_namelist(nml_filename)
    call surface_read_namelist(nml_filename)
    call poisson_solver_read_namelist(nml_filename)
    call checksim_read_namelist(nml_filename)

    ! Add-on modules
    call aerosol_read_namelist(nml_filename)
    call lateral_sponge_read_namelist(nml_filename)
    call microphysics_read_namelist(nml_filename)
  
  end subroutine read_namelists

end module modnamelist