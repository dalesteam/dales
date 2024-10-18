!> Module with various handy functions/subroutines
module utils

  implicit none

  private

  interface findval
    module procedure :: findval_character
    module procedure :: findval_integer
    module procedure :: findval_logical
    module procedure :: findval_real4
    module procedure :: findval_real8
  end interface findval
  
  public :: findval
  public :: to_lower
  public :: to_upper

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
  pure real(4) function findval_real4(key, keys, values, defltvalue)
    character(*), intent(in)           :: key 
    character(*), intent(in)           :: keys(:)
    real(4),      intent(in)           :: values(:)
    real(4),      intent(in), optional :: defltvalue

    character(6) :: fkey
    integer      :: idx(1)

    fkey = trim(key)

    idx = findloc(keys, fkey)
    if (idx(1) /= 0 .and. idx(1) <= size(values)) then
      findval_real4 = values(idx(1))
    else
      if (present(defltvalue)) then
        findval_real4 = defltvalue 
      else
        findval_real4 = 0.0 
      end if
    end if
  end function findval_real4

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
  pure real(8) function findval_real8(key, keys, values, defltvalue)
    character(*), intent(in)           :: key 
    character(*), intent(in)           :: keys(:)
    real(8),      intent(in)           :: values(:)
    real(8),      intent(in), optional :: defltvalue

    character(6) :: fkey
    integer      :: idx(1)

    fkey = trim(key)

    idx = findloc(keys, fkey)
    if (idx(1) /= 0 .and. idx(1) <= size(values)) then
      findval_real8 = values(idx(1))
    else
      if (present(defltvalue)) then
        findval_real8 = defltvalue 
      else
        findval_real8 = 0.0 
      end if
    end if
  end function findval_real8

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

      !> Transform lowercase letters to uppercase ones in a string.
    !!
    !! Leave capitals, numbers and punctuation untouched
    !!
    !! @param[in] strIn String to convert to uppercase
    !! @returns Uppercase string
    !!
    !! @see https://stackoverflow.com/questions/10759375/how-can-i-write-a-to-upper-or-to-lower-function-in-f90
  function to_upper(strIn) result(strOut)
    ! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
    ! Original author: Clive Page

    implicit none

    character(len=*), intent(in) :: strIn
    character(len=len(strIn)) :: strOut
    integer :: i,j

    do i = 1, len(strIn)
        j = iachar(strIn(i:i))
        if (j>= iachar("a") .and. j<=iachar("z") ) then
            strOut(i:i) = achar(iachar(strIn(i:i))-32)
        else
            strOut(i:i) = strIn(i:i)
        end if
    end do

  end function to_upper

  !> Transform uppercase letters to lowercase ones in a string.
  !!
  !! Leave lowercase, numbers and punctuation untouched
  !!
  !! @param[in] strIn String to convert to lowercase
  !! @returns lowercase string
  !!
  !! @see https://stackoverflow.com/questions/10759375/how-can-i-write-a-to-upper-or-to-lower-function-in-f90
  function to_lower(strIn) result(strOut)
      ! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
      ! Original author: Clive Page

      implicit none

      character(len=*), intent(in) :: strIn
      character(len=len(strIn)) :: strOut
      integer :: i,j

      do i = 1, len(strIn)
          j = iachar(strIn(i:i))
          if (j>= iachar("A") .and. j<=iachar("Z") ) then
              strOut(i:i) = achar(iachar(strIn(i:i))+32)
          else
              strOut(i:i) = strIn(i:i)
          end if
      end do

  end function to_lower

end module utils