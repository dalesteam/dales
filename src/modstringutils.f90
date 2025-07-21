module modstringutils

  use, intrinsic :: iso_fortran_env, only: real32, real64

  implicit none

  private

  interface number2string
    module procedure :: int2string
    module procedure :: real2string
    module procedure :: double2string
  end interface number2string

  character(len=*), parameter :: modname = 'modstringutils'

  public :: number2string

contains

  function int2string(n, opt_fmt) result(int_string)
    integer,          intent(in) :: n
    character(len=*), intent(in), optional :: opt_fmt

    character(len=11) :: int_string
    character(len=11) :: fmt

    if (present(opt_fmt)) then
      fmt = opt_fmt
    else
      fmt = '(i11)'
    end if
    write(int_string, fmt) n
    int_string = adjustl(int_string)

  end function int2string

  function real2string(n, opt_fmt) result(real_string)
    real(real32),     intent(in) :: n
    character(len=*), intent(in), optional :: opt_fmt

    character(len=32) :: real_string
    character(len=11) :: fmt

    if (present(opt_fmt)) then
      fmt = opt_fmt
    else
      fmt = '(g32.5)'
    end if
    write(real_string, fmt) n
    real_string = adjustl(real_string)

  end function real2string

  function double2string(n, opt_fmt) result(double_string)
    real(real64),     intent(in) :: n
    character(len=*), intent(in), optional :: opt_fmt

    character(len=32) :: double_string
    character(len=11) :: fmt

    if (present(opt_fmt)) then
      fmt = opt_fmt
    else
      fmt = '(g32.8)'
    end if
    write(double_string, fmt) n
    double_string = adjustl(double_string)

  end function double2string

end module modstringutils