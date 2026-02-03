module modaerosol_sources

  use modaerosol_mode_t, only: aerosol_mode_t
  use modaerosol_common, only: calc_median_diameter, aerosol_densities
  use modglobal,         only: pi
  use modprecision,      only: field_r
  use fortran_support,   only: finish

  implicit none

  private

  public :: aerosol_point_source
  public :: NUMBER_SOURCE, MASS_SOURCE

  character(len=*), parameter :: modname = 'modaerosol_sources'

  integer, parameter :: NUMBER_SOURCE = 0
  integer, parameter :: MASS_SOURCE = 1

contains

  !> Apply a point source of aerosols to a given mode.
  !!
  !! Note: since we don't have any mechanism to transfer aerosols between modes due to
  !! growth of particles, we assume here that the source aerosol follows the 
  !! same distribution as the target mode. I.e., the source aerosol follows a 
  !! lognormal distribution with the same median diameter and geometric
  !! standard deviation as the target mode.
  subroutine aerosol_point_source(mode, ispecies, source_strength, location, &
                                  source_type)

    type(aerosol_mode_t), intent(inout) :: mode            !< Aerosol mode to which the source is applied
    integer,              intent(in)    :: ispecies        !< Aerosol species index
    real(field_r),        intent(in)    :: source_strength !< Strength of the source [kg/kg/s or #/kg/s]
    integer,              intent(in)    :: location(3)     !< Grid location (i,j,k) where the source is applied
    integer,              intent(in)    :: source_type     !< Type of source: NUMBER_SOURCE or MASS_SOURCE

    character(len=*), parameter :: routine = modname//'/aerosol_point_source'

    integer :: i, j, k

    integer       :: idx(1) !< Index of species in mode
    real(field_r) :: rho    !< Aerosol density of source [kg/m3]
    real(field_r) :: dm     !< Median diameter of mode [m]
    real(field_r) :: fac    !< Conversion factor between number and mass source

    i = location(1)
    j = location(2)
    k = location(3)

    rho = aerosol_densities(ispecies)
    dm = calc_median_diameter(mode%n(i,j,k), mode%q(:,i,j,k), mode%rho, &
                              mode%sig_g)
    fac = pi / 6 * rho * dm**3 * exp(4.5 * log(mode%sig_g)**2)

    idx = findloc(mode%itype, ispecies)

    select case (source_type)
      case (NUMBER_SOURCE)
        mode%qp(idx(1),i,j,k) = mode%qp(idx(1),i,j,k) + source_strength * fac
      case (MASS_SOURCE)
        mode%np(i,j,k) = mode%np(i,j,k) + source_strength / fac
      case default
        call finish(routine, 'Invalid source_type')
    end select

  end subroutine aerosol_point_source

end module modaerosol_sources