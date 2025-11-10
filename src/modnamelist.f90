!> Read namelists and check settings.
module modnamelist

  use modchecksim,     only: checksim_read_namelist
  use modmicrophysics, only: microphysics_read_namelist

  implicit none

  private

  public :: read_namelists

contains

  !> Read all namelists.
  subroutine read_namelists(nml_filename)

    character(len=*), intent(in) :: nml_filename

    ! Core modules
    call checksim_read_namelist(nml_filename)

    ! Add-on modules
    call microphysics_read_namelist(nml_filename)
  
  end subroutine read_namelists

end module modnamelist