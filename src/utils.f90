!> Module with various handy functions/subroutines
module utils

  implicit none

  private

  interface findval
    module procedure :: findval_character
    module procedure :: findval_integer
    module procedure :: findval_logical
    module procedure :: findval_real
  end interface findval
  
  public :: findval

contains

  !> Find a value in an array of values based on a key in an array of keys
  !!
  !! Lookup a value in one array (`values`) at the index position of `key` in the array `keys`.
  !! Provide the optional default value (`defltvalue`), that is returned when `key` isn't found.
  !! When no default is given, zero is returned.
  !!
  !! @param[in] key The key to be looked up
  !! @param[in] keys The array of keys
  !! @param[in] values The array of values
  !! @param[in] defltvalue The default value (optional, default 0.0)
  !! @return The value corresponding to the key
  pure real function findval_real(key, keys, values, defltvalue)
    character(*), intent(in)           :: key 
    character(*), intent(in)           :: keys(:)
    real,         intent(in)           :: values(:)
    real,         intent(in), optional :: defltvalue

    character(6) :: fkey
    integer      :: idx(1)

    fkey = trim(key)

    idx = findloc(keys, fkey)
    if (idx(1) /= 0 .and. idx(1) <= size(values)) then
      findval_real = values(idx(1))
    else
      if (present(defltvalue)) then
        findval_real = defltvalue 
      else
        findval_real = 0.0 
      end if
    end if
  end function findval_real

  !> Find a value in an array of values based on a key in an array of keys
  !!
  !! Lookup a value in one array (`values`) at the index position of `key` in the array `keys`.
  !! Provide the optional default value (`defltvalue`), that is returned when `key` isn't found.
  !! When no default is given, zero is returned.
  !!
  !! @param[in] key The key to be looked up
  !! @param[in] keys The array of keys
  !! @param[in] values The array of values
  !! @param[in] defltvalue The default value (optional, default -1)
  !! @return The value corresponding to the key
  pure integer function findval_integer(key, keys, values, defltvalue)
    character(*), intent(in)           :: key 
    character(*), intent(in)           :: keys(:)
    integer,      intent(in)           :: values(:)
    integer,      intent(in), optional :: defltvalue

    character(6) :: fkey
    integer      :: idx(1)

    fkey = trim(key)

    idx = findloc(keys, fkey)
    if (idx(1) /= 0 .and. idx(1) <= size(values)) then
      findval_integer = values(idx(1))
    else
      if (present(defltvalue)) then
        findval_integer = defltvalue 
      else
        findval_integer = -1 
      end if
    end if
  end function findval_integer

  !> Find a value in an array of values based on a key in an array of keys
  !!
  !! Lookup a value in one array (`values`) at the index position of `key` in the array `keys`.
  !! Provide the optional default value (`defltvalue`), that is returned when `key` isn't found.
  !! When no default is given, zero is returned.
  !!
  !! @param[in] key The key to be looked up
  !! @param[in] keys The array of keys
  !! @param[in] values The array of values
  !! @param[in] defltvalue The default value (optional, default .False.)
  !! @return The value corresponding to the key
  pure logical function findval_logical(key, keys, values, defltvalue)
    character(*), intent(in)           :: key 
    character(*), intent(in)           :: keys(:)
    logical,      intent(in)           :: values(:)
    logical,      intent(in), optional :: defltvalue

    character(6) :: fkey
    integer      :: idx(1)

    fkey = trim(key)

    idx = findloc(keys, fkey)
    if (idx(1) /= 0 .and. idx(1) <= size(values)) then
      findval_logical = values(idx(1))
    else
      if (present(defltvalue)) then
        findval_logical = defltvalue 
      else
        findval_logical = .False. 
      end if
    end if
  end function findval_logical

  !> Find a value in an array of values based on a key in an array of keys
  !!
  !! Lookup a value in one array (`values`) at the index position of `key` in the array `keys`.
  !! Provide the optional default value (`defltvalue`), that is returned when `key` isn't found.
  !! When no default is given, zero is returned.
  !!
  !! @param[in] key The key to be looked up
  !! @param[in] keys The array of keys
  !! @param[in] values The array of values
  !! @param[in] defltvalue The default value (optional, default .False.)
  !! @return The value corresponding to the key
  pure character(32) function findval_character(key, keys, values, defltvalue)
    character(*), intent(in)           :: key 
    character(*), intent(in)           :: keys(:)
    character(*), intent(in)           :: values(:)
    character(6), intent(in), optional :: defltvalue

    character(6) :: fkey
    integer      :: idx(1)

    fkey = trim(key)
    idx = findloc(keys, fkey)
    if (idx(1) /= 0 .and. idx(1) <= size(values)) then
      findval_character = trim(values(idx(1)))
    else
      if (present(defltvalue)) then
        findval_character = defltvalue 
      else
        findval_character = 'dummy'
      end if
    end if
  end function findval_character

end module utils