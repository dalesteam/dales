module modtracer_type

  use modprecision, only: field_r

  implicit none

  type T_tracer
    character(len=16) :: tracname           !< Tracer name
    character(len=64) :: traclong = "dummy long name" !< Tracer long name
    character(len=16) :: unit = "dummy unit" !< Tracer unit
    real(field_r)     :: molar_mass = -999. !< Molecular mass of tracer (g mol-1)
    integer           :: trac_idx = -1      !< Tracer index in sv0, svm, svp
    logical           :: lemis = .false.    !< Boolean if tracer is emitted
    logical           :: lreact = .false.   !< Boolean if tracer is reactive
    logical           :: ldep = .false.     !< Boolean if tracer is deposited
    logical           :: lags = .false.     !< Boolean if in A-gs
    logical           :: lmicro = .false.   !< Boolean if in cloud microphysics
    real(field_r)     :: wsvsurf = 0        !< Kinematic surface flux (- m/s)
    logical           :: lnudge = .false.   !< Boolean if tracer is nudged
    logical           :: lnudge_init = .true. !< Nudge to initial conditions (true) or to a profile (false).
    integer           :: inudge_loc = 0     !< Location of nudging (see modnudge for options)
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
    write(*,*) "  wsvsurf    : ", self%wsvsurf
    write(*,*) "  lnudge     : ", self%lnudge
    write(*,*) "  lnudge_init: ", self%lnudge_init
    write(*,*) "  inudge_loc : ", self%inudge_loc
  end subroutine tracer_print_properties

end module modtracer_type