module modtracer_type

  use modprecision, only: field_r
 
  implicit none
  
  type T_tracer
  ! Fixed tracer properties
      ! Tracer name
      character(len=16) :: tracname
      ! Tracer long name
      character(len=32) :: traclong
      ! Tracer unit
      character(len=16) :: unit     
      ! Moleculare mass of tracer (g mol-1)
      real              :: molar_mass
      ! Tracer index in sv0, svm, svp
      integer           :: trac_idx
      ! Boolean if tracer is emitted 
      logical           :: lemis
      ! Boolean if tracer is reactive
      logical           :: lreact
      ! Boolean if tracer is deposited
      logical           :: ldep
      ! Boolean if in A-gs
      logical           :: lags   
      ! Boolean if in cloud microphysics
      logical           :: lmicro     
      ! Boolean if tracer is an aerosol
      logical           :: laero
      ! ! Static tracer properties:
      ! real :: diffusivity
      real(field_r) :: rho
      real(field_r) :: kappa
      character(3) :: mode !< Comma-separated list of the modes

  contains
    procedure :: print_properties => tracer_print_properties
  end type T_tracer

contains

  subroutine tracer_print_properties(self)

    class(T_tracer), intent(in) :: self

    write(*,*) "Tracer: ", self%tracname
    write(*,*) "  long name: ", trim(self%traclong)
    write(*,*) "  is aerosol: ", self%laero

    if (self%laero) then
      write(*,*) "  mode: ", trim(self%mode)
      write(*,'(a12,f6.1)') "  density: ", self%rho
      write(*,'(a19,f4.2)') "  hygroscopicity: ", self%kappa
    end if

  end subroutine tracer_print_properties

end module modtracer_type
