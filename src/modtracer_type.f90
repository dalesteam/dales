module modtracer_type

  use modprecision, only: field_r
 
  implicit none
  
  type T_tracer
  ! Fixed tracer properties
      ! Tracer name
      character(len=16) :: tracname
      ! Tracer long name
      character(len=64) :: traclong="dummy long name"
      ! Tracer unit
      character(len=16) :: unit="dummy unit"
      ! Moleculare mass of tracer (g mol-1)
      real(field_r)     :: molar_mass=-999.
      ! Tracer index in sv0, svm, svp
      integer           :: trac_idx=-1
      ! Boolean if tracer is emitted 
      logical           :: lemis=.false.
      ! Boolean if tracer is reactive
      logical           :: lreact=.false.
      ! Boolean if tracer is deposited
      logical           :: ldep=.false.
      ! Boolean if in A-gs
      logical           :: lags=.false.
      ! Boolean if in cloud microphysics
      logical           :: lmicro=.false.
      ! ! Static tracer properties:
      ! real :: diffusivity

  contains
    procedure :: print_properties => tracer_print_properties
  end type T_tracer

contains

  subroutine tracer_print_properties(self)

    class(T_tracer), intent(in) :: self

    write(*,*) "Tracer: ", self%tracname
    write(*,*) "  long name  : ", trim(self%traclong)
    write(*,*) "  unit       : ", trim(self%unit)
    write(*,*) "  molar mass : ", self%molar_mass
    write(*,*) "  index      : ", self%trac_idx
    write(*,*) "  lemis      : ", self%lemis
    write(*,*) "  lreact     : ", self%lreact
    write(*,*) "  ldep       : ", self%ldep
    write(*,*) "  lags       : ", self%lags
    write(*,*) "  lmicro     : ", self%lmicro

  end subroutine tracer_print_properties

end module modtracer_type
