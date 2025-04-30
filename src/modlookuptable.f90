!> A Fortran 90 module for creating lookup tables. These tables can be used to
!> efficiently interpolate one or more values.
!>
!> Author: Jannis Teunissen
module modlookuptable
  use modprecision, only: field_r
  implicit none
  private

  ! The precision of the real numbers used in the tables
  integer, parameter :: dp = field_r

  ! Table spacing
  integer, parameter, public :: LT_xspacing_linear = 1
  integer, parameter, public :: LT_xspacing_quadratic = 2
  integer, parameter, public :: LT_xspacing_cubic = 3

  ! ** Routines for finding indices in sorted lists **
  public :: find_index_linear
  public :: find_index_bsearch
  public :: find_index_adaptive

  !> The lookup table type. There can be one or more columns, for which values
  !> can be looked up for a given 'x-coordinate'.
  type LT_t
     integer  :: n_points !< The number of points
     integer  :: n_cols   !< The number of columns
     integer  :: xspacing  !< Type of table spacing
     real(dp) :: x_min    !< The minimum lookup coordinate
     real(dp) :: inv_fac  !< The inverse x-spacing
     logical  :: extrapolate_above !< Linearly extrapolate above x_max

     ! The table is stored in two ways, to speed up different types of lookups.
     real(dp), allocatable :: x(:) !< The x values in the table
     real(dp), allocatable :: cols_rows(:, :)
     real(dp), allocatable :: rows_cols(:, :)
  end type LT_t

  !> The 2D lookup table type
  type LT2_t
     integer               :: n_points(2) !< The size of the table
     integer               :: n_cols      !< The number of columns/variables
     integer               :: xspacing(2)  !< Type of table spacing
     real(dp)              :: x_min(2)    !< The minimum lookup coordinate
     real(dp)              :: inv_fac(2)  !< The inverse x-spacing
     real(dp), allocatable :: x1(:)       !< The x values in the table
     real(dp), allocatable :: x2(:)       !< The x values in the table
     real(dp), allocatable :: rows_cols(:, :, :)
  end type LT2_t

  !> Type to indicate a location in the lookup table, which can be used to speed
  !> up multiple lookups of different columns.
  type LT_loc_t
     integer  :: low_ix   !< The x-value lies between low_ix and low_ix+1
     real(dp) :: low_frac !< The distance from low_ix (up to low_ix+1), given
                          !< as a real number between 0 and 1.
  end type LT_loc_t

  !> Type to indicate a location in a 2D lookup table
  type LT2_loc_t
     !> The x-value lies between low_ix and low_ix+1
     integer  :: low_ix(2)
     !> The distance from low_ix (up to low_ix+1), given as a real number
     !> between 0 and 1.
     real(dp) :: low_frac(2)
  end type LT2_loc_t

  ! Public types
  public :: LT_t
  public :: LT_loc_t
  public :: LT2_t
  public :: LT2_loc_t

  ! Generic methods
  public :: LT_get_spaced_data  ! Convert values to regularly spaced
  public :: LT_lin_interp_list  ! Linearly interpolate a list

  ! Public methods
  public :: LT_create          ! Create a new lookup table
  public :: LT_set_col         ! Set one table column
  public :: LT_set_col_data    ! Set one table column
  public :: LT_add_col         ! Add a column
  public :: LT_get_loc         ! Get the index (row) of a value
  public :: LT_get_col         ! Interpolate one column
  public :: LT_get_mcol        ! Interpolate multiple columns
  public :: LT_get_col_at_loc  ! Get one column at location
  public :: LT_get_mcol_at_loc ! Get multiple columns at location
  public :: LT_to_file         ! Store lookup table in file
  public :: LT_from_file       ! Restore lookup table from file

  ! Public methods for 2D tables
  public :: LT2_create         ! Create a new lookup table
  public :: LT2_set_col        ! Set one table column
  public :: LT2_set_col_data   ! Set one table column
  public :: LT2_get_loc        ! Get the index (row) of a value
  public :: LT2_get_col        ! Interpolate one column
  public :: LT2_get_col_at_loc ! Get one column at location
  public :: LT2_get_col_inline ! Interpolate one column, inlined version

contains

  ! ** Routines for finding indices **

  !> Linear search of sorted list for the smallest ix such that list(ix) >= val.
  !> On failure, returns size(list)+1
  pure function find_index_linear(list, val) result(ix)
    real(dp), intent(in) :: list(:) !< Sorted list
    real(dp), intent(in) :: val     !< Value to search for
    integer              :: ix

    do ix = 1, size(list)
       if (list(ix) >= val) exit
    end do
  end function find_index_linear

  !> Binary search of sorted list for the smallest ix such that list(ix) >= val.
  !> On failure, returns size(list)+1
  pure function find_index_bsearch(list, val) result(ix)
    real(dp), intent(in) :: list(:) !< Sorted list
    real(dp), intent(in) :: val     !< Value to search for
    integer              :: ix, i_min, i_max, i_middle

    i_min = 1
    i_max = size(list)

    do while (i_min < i_max)
       ! This safely performs: i_middle = (i_max + i_min) / 2
       i_middle = i_min + ishft(i_max - i_min, -1)

       if (list(i_middle) >= val) then
          i_max = i_middle
       else
          i_min = i_middle + 1
       end if
    end do

    ix = i_min
    if (val > list(ix)) ix = size(list) + 1
  end function find_index_bsearch

  !> Adaptive search (combination of linear and binary search) of sorted list
  !> for the smallest ix such that list(ix) >= val. On failure, returns
  !> size(list)+1
  pure function find_index_adaptive(list, val) result(ix)
    real(dp), intent(in) :: list(:) !< Sorted list
    real(dp), intent(in) :: val     !< Value to search for
    integer              :: ix
    integer, parameter   :: binary_search_limit = 40

    if (size(list) < binary_search_limit) then
       ix = find_index_linear(list, val)
    else
       ix = find_index_bsearch(list, val)
    end if
  end function find_index_adaptive

  !> Compute by use of linear interpolation the value in the middle
  ! of a domain D = [x_list(1) , x_list(size(x_list))].
  ! If x_value is left of domain  D,
  ! then the value becomes the value at the left side of D,
  ! if x_value is right of domain D,
  ! then the value becomes the value at the rigth side of D
  pure subroutine LT_lin_interp_list(x_list, y_list, x_value, y_value)
    real(dp), intent(in)  :: x_list(:), y_list(:)
    real(dp), intent(in)  :: x_value
    real(dp), intent(out) :: y_value

    integer               :: ix, iMin, iMax
    real(dp)              :: temp

    iMin = 1
    iMax = size(x_list)

    if (x_value <= x_list(iMin)) then
       y_value = y_list(iMin)
    else if (x_value >= x_list(iMax)) then
       y_value = y_list(iMax)
    else
       ix = find_index_adaptive(x_list, x_value)
       temp = (x_value - x_list(ix-1)) / (x_list(ix) - x_list(ix-1))
       y_value = (1 - temp) * y_list(ix-1) + temp * y_list(ix)
    end if
  end subroutine LT_lin_interp_list

  ! ** 1D lookup table routines **

  !> This function returns a new lookup table
  function LT_create(x_min, x_max, n_points, n_cols, xspacing, &
       extrapolate_above) result(my_lt)
    real(dp), intent(in) :: x_min  !< Minimum x-coordinate
    real(dp), intent(in) :: x_max  !< Maximum x-coordinate
    integer, intent(in)  :: n_points !< How many x-values to store
    integer, intent(in)  :: n_cols !< Number of variables that will be looked up
    integer, intent(in), optional :: xspacing !< Spacing of data
    !> Linearly extrapolate above x_max
    logical, intent(in), optional :: extrapolate_above
    type(LT_t)           :: my_lt

    if (x_max <= x_min) error stop "x_max should be > x_min"
    if (n_points <= 1) error stop "n_points should be bigger than 1"

    my_lt%n_points = n_points
    my_lt%n_cols   = n_cols
    my_lt%x_min    = x_min
    my_lt%xspacing  = LT_xspacing_linear
    if (present(xspacing)) my_lt%xspacing = xspacing
    my_lt%extrapolate_above = .false.
    if (present(extrapolate_above)) my_lt%extrapolate_above = extrapolate_above

    allocate(my_lt%x(n_points))
    call table_set_x(n_points, my_lt%xspacing, x_min, x_max, &
         my_lt%x, my_lt%inv_fac)

    allocate(my_lt%cols_rows(n_cols, n_points))
    allocate(my_lt%rows_cols(n_points, n_cols))
    my_lt%cols_rows = 0
    my_lt%rows_cols = 0
  end function LT_create

  subroutine table_set_x(n_points, xspacing, x_min, x_max, x, inv_fac)
    integer, intent(in)   :: n_points
    integer, intent(in)   :: xspacing
    real(dp), intent(in)  :: x_min, x_max
    real(dp), intent(out) :: x(n_points)
    real(dp), intent(out) :: inv_fac

    select case (xspacing)
    case (LT_xspacing_linear)
       inv_fac = (n_points - 1)/(x_max - x_min)
    case (LT_xspacing_quadratic)
       inv_fac = (n_points - 1.0_dp)**2/(x_max - x_min)
    case (LT_xspacing_cubic)
       inv_fac = (n_points - 1.0_dp)**3/(x_max - x_min)
    case default
       error stop "Unknown spacing"
    end select

    x = get_x(x_min, x_max, n_points, xspacing)
  end subroutine table_set_x

  !> Linearly interpolate the (x, y) input data to the new_x coordinates
  function LT_get_spaced_data(in_x, in_y, new_x) result(out_yy)
    real(dp), intent(in) :: in_x(:), in_y(:), new_x(:)
    real(dp)             :: out_yy(size(new_x))
    integer              :: ix, n

    n = size(in_x)
    if (n < 2) error stop "size(in_x) < 2"
    if (size(in_x) /= size(in_y)) error stop "in_x and in_y not of same size"
    if (minval(in_x(2:) - in_x(1:n-1)) <= 0) &
         error stop "in_x should strictly increase"

    do ix = 1, size(new_x)
       call LT_lin_interp_list(in_x, in_y, new_x(ix), out_yy(ix))
    end do
  end function LT_get_spaced_data

  !> Fill the column with index col_ix after linearly interpolating
  subroutine LT_set_col(my_lt, col_ix, x, y)
    type(LT_t), intent(inout) :: my_lt
    integer, intent(in)       :: col_ix
    real(dp), intent(in)      :: x(:), y(:)

    if (col_ix < 0 .or. col_ix > my_lt%n_cols) &
         error stop "should have 1 <= col_ix <= n_cols"

    my_lt%cols_rows(col_ix, :) = LT_get_spaced_data(x, y, my_lt%x)
    my_lt%rows_cols(:, col_ix) = my_lt%cols_rows(col_ix, :)
  end subroutine LT_set_col

  !> Fill the column with index col_ix with y data
  subroutine LT_set_col_data(my_lt, col_ix, y)
    type(LT_t), intent(inout) :: my_lt
    integer, intent(in)       :: col_ix
    real(dp), intent(in)      :: y(:)

    if (col_ix < 0 .or. col_ix > my_lt%n_cols) &
         error stop "should have 1 <= col_ix <= n_cols"
    if (size(y) /= my_lt%n_points) error stop "size(y) /= number of rows"

    my_lt%cols_rows(col_ix, :) = y
    my_lt%rows_cols(:, col_ix) = y
  end subroutine LT_set_col_data

  !> Add a new column by linearly interpolating the (x, y) data
  subroutine LT_add_col(my_lt, x, y)
    type(LT_t), intent(inout) :: my_lt
    real(dp), intent(in)      :: x(:), y(:)
    type(LT_t)                :: temp_lt

    temp_lt = my_lt
    deallocate(my_lt%cols_rows)
    deallocate(my_lt%rows_cols)
    allocate(my_lt%rows_cols(my_lt%n_points, my_lt%n_cols+1))
    allocate(my_lt%cols_rows(my_lt%n_cols+1, my_lt%n_points))

    my_lt%cols_rows(1:my_lt%n_cols, :) = temp_lt%cols_rows
    my_lt%rows_cols(:, 1:my_lt%n_cols) = temp_lt%rows_cols
    my_lt%n_cols                       = my_lt%n_cols + 1
    my_lt%cols_rows(my_lt%n_cols, :)   = LT_get_spaced_data(x, y, my_lt%x)
    my_lt%rows_cols(:, my_lt%n_cols)   = my_lt%cols_rows(my_lt%n_cols, :)
  end subroutine LT_add_col

  !> Returns the x-coordinates of the lookup table
  pure function get_x(x_min, x_max, n_points, xspacing) result(x)
    real(dp), intent(in) :: x_min, x_max
    integer, intent(in)  :: n_points, xspacing
    real(dp)             :: x(n_points), tmp
    integer              :: ix

    tmp = 1.0_dp / (n_points-1)

    select case (xspacing)
    case (LT_xspacing_linear)
       do ix = 1, n_points
          x(ix) = (ix-1) * tmp
       end do
    case (LT_xspacing_quadratic)
       do ix = 1, n_points
          x(ix) = ((ix-1) * tmp)**2
       end do
    case (LT_xspacing_cubic)
       do ix = 1, n_points
          x(ix) = ((ix-1) * tmp)**3
       end do
    end select

    x = x_min + x * (x_max - x_min)
  end function get_x

  !> Get a location in the lookup table
  elemental function LT_get_loc(my_lt, x) result(my_loc)
    type(LT_t), intent(in) :: my_lt
    real(dp), intent(in)   :: x
    type(LT_loc_t)         :: my_loc
    real(dp)               :: frac
    real(dp), parameter    :: one_third = 1/3.0_dp

    frac = (x - my_lt%x_min) * my_lt%inv_fac

    select case (my_lt%xspacing)
    case (LT_xspacing_quadratic)
       if (frac > 0) frac = sqrt(frac)
    case (LT_xspacing_cubic)
       if (frac > 0) frac = frac**one_third
    end select

    ! Check bounds
    if (frac <= 0) then
       my_loc%low_ix   = 1
       my_loc%low_frac = 1
    else if (frac >= my_lt%n_points - 1) then
       my_loc%low_ix   = my_lt%n_points - 1
       if (my_lt%extrapolate_above) then
          my_loc%low_frac = my_loc%low_ix - frac
       else
          my_loc%low_frac = 0
       end if
    else
       my_loc%low_ix   = ceiling(frac)
       my_loc%low_frac = my_loc%low_ix - frac
    end if

  end function LT_get_loc

  !> Get the values of all columns at x
  pure function LT_get_mcol(my_lt, x) result(col_values)
    type(LT_t), intent(in) :: my_lt
    real(dp), intent(in)   :: x
    real(dp)               :: col_values(my_lt%n_cols)
    type(LT_loc_t)         :: loc

    loc        = LT_get_loc(my_lt, x)
    col_values = LT_get_mcol_at_loc(my_lt, loc)
  end function LT_get_mcol

  !> Get the value of a single column at x
  elemental function LT_get_col(my_lt, col_ix, x) result(col_value)
    type(LT_t), intent(in) :: my_lt
    integer, intent(in)    :: col_ix
    real(dp), intent(in)   :: x
    real(dp)               :: col_value
    type(LT_loc_t)         :: loc

    loc       = LT_get_loc(my_lt, x)
    col_value = LT_get_col_at_loc(my_lt, col_ix, loc)
  end function LT_get_col

  !> Get the values of all columns at a location
  pure function LT_get_mcol_at_loc(my_lt, loc) result(col_values)
    type(LT_t), intent(in)     :: my_lt
    type(LT_loc_t), intent(in) :: loc
    real(dp)                   :: col_values(my_lt%n_cols)

    col_values = loc%low_frac * my_lt%cols_rows(:, loc%low_ix) + &
         (1-loc%low_frac) * my_lt%cols_rows(:, loc%low_ix+1)
  end function LT_get_mcol_at_loc

  !> Get the value of a single column at a location
  elemental function LT_get_col_at_loc(my_lt, col_ix, loc) result(col_value)
    type(LT_t), intent(in)     :: my_lt
    integer, intent(in)        :: col_ix
    type(LT_loc_t), intent(in) :: loc
    real(dp)                   :: col_value

    col_value = loc%low_frac * my_lt%rows_cols(loc%low_ix, col_ix) + &
         (1-loc%low_frac) * my_lt%rows_cols(loc%low_ix+1, col_ix)
  end function LT_get_col_at_loc

  !> Write the lookup table to file (in binary, potentially unportable)
  subroutine LT_to_file(my_lt, filename)
    type(LT_t), intent(in)       :: my_lt
    character(len=*), intent(in) :: filename
    integer                      :: my_unit

    open(newunit=my_unit, file=trim(filename), form='UNFORMATTED', &
         access='STREAM', status='REPLACE')
    write(my_unit) my_lt%n_points, my_lt%n_cols
    write(my_unit) my_lt%x_min, my_lt%inv_fac, my_lt%xspacing
    write(my_unit) my_lt%x, my_lt%cols_rows
    close(my_unit)
  end subroutine LT_to_file

  !> Read the lookup table from file (in binary, potentially unportable)
  subroutine LT_from_file(my_lt, filename)
    type(LT_t), intent(inout)    :: my_lt
    character(len=*), intent(in) :: filename
    integer                      :: my_unit

    open(newunit=my_unit, file=trim(filename), form='UNFORMATTED', &
         access='STREAM', status='OLD')
    read(my_unit) my_lt%n_points, my_lt%n_cols
    read(my_unit) my_lt%x_min, my_lt%inv_fac, my_lt%xspacing

    allocate(my_lt%x(my_lt%n_points))
    allocate(my_lt%cols_rows(my_lt%n_cols, my_lt%n_points))
    allocate(my_lt%rows_cols(my_lt%n_points, my_lt%n_cols))

    read(my_unit) my_lt%x, my_lt%cols_rows
    my_lt%rows_cols = transpose(my_lt%cols_rows)

    close(my_unit)
  end subroutine LT_from_file

  ! ** 2D lookup table routines **

  !> This function returns a new lookup table
  function LT2_create(x_min, x_max, n_points, n_cols, xspacing) result(my_lt)
    real(dp), intent(in)          :: x_min(2)    !< Minimum coordinate
    real(dp), intent(in)          :: x_max(2)    !< Maximum coordinate
    integer, intent(in)           :: n_points(2) !< How many values to store
    integer, intent(in)           :: n_cols      !< Number of variables that will be looked up
    integer, intent(in), optional :: xspacing(2)  !< Spacing of data
    type(LT2_t)                   :: my_lt

    if (any(x_max <= x_min)) stop "LT2_create error: x_max <= x_min"
    if (any(n_points <= 1)) stop "LT2_create error: n_points <= 1"

    !$acc enter data create(my_lt) async

    my_lt%n_points = n_points
    my_lt%n_cols   = n_cols
    my_lt%x_min    = x_min
    my_lt%xspacing  = LT_xspacing_linear
    if (present(xspacing)) my_lt%xspacing = xspacing

    allocate(my_lt%x1(n_points(1)))
    allocate(my_lt%x2(n_points(2)))

    call table_set_x(n_points(1), my_lt%xspacing(1), x_min(1), x_max(1), &
         my_lt%x1, my_lt%inv_fac(1))
    call table_set_x(n_points(2), my_lt%xspacing(2), x_min(2), x_max(2), &
         my_lt%x2, my_lt%inv_fac(2))

    allocate(my_lt%rows_cols(n_points(1), n_points(2), n_cols))
    my_lt%rows_cols = 0

    !$acc update device(my_lt%x_min, my_lt%xmax, my_lt%n_points, my_lt%n_cols) &
    !$acc update device(xspacing) if(present(xspacing))
    !$acc enter data copyin(my_lt%rows_cols)

  end function LT2_create

  !> Fill the column with index col_ix using linearly interpolated data
  subroutine LT2_set_col(my_lt, col_ix, x1, x2, y)
    type(LT2_t), intent(inout) :: my_lt
    integer, intent(in)        :: col_ix
    real(dp), intent(in)       :: x1(:), x2(:), y(:, :)
    real(dp), allocatable      :: tmp(:, :)
    integer                    :: ix

    allocate(tmp(my_lt%n_points(1), size(x2)))

    ! Interpolate along first coordinate
    do ix = 1, size(x2)
       tmp(:, ix) = LT_get_spaced_data(x1, y(:, ix), my_lt%x1)
    end do

    ! Interpolate along second coordinate
    do ix = 1, my_lt%n_points(1)
       my_lt%rows_cols(ix, :, col_ix) = &
            LT_get_spaced_data(x2, tmp(ix, :), my_lt%x2)
    end do

    !$acc update device(my_lt%rows_cols)

  end subroutine LT2_set_col

  !> Fill the column with index col_ix with y data
  subroutine LT2_set_col_data(my_lt, col_ix, y)
    type(LT2_t), intent(inout) :: my_lt
    integer, intent(in)        :: col_ix
    real(dp), intent(in)       :: y(:, :)

    if (col_ix < 0 .or. col_ix > my_lt%n_cols) &
         error stop "should have 1 <= col_ix <= n_cols"
    if (any(shape(y) /= my_lt%n_points)) error stop "shape(y) /= n_points"

    my_lt%rows_cols(:, :, col_ix) = y
  end subroutine LT2_set_col_data

  !> Get a location in the lookup table
  elemental function LT2_get_loc(my_lt, x1, x2) result(my_loc)
    !$acc routine seq
    type(LT2_t), intent(in) :: my_lt
    real(dp), intent(in)    :: x1, x2
    type(LT2_loc_t)         :: my_loc
    real(dp)                :: frac(2)

    frac            = ([x1, x2] - my_lt%x_min) * my_lt%inv_fac
    my_loc%low_ix   = ceiling(frac)
    my_loc%low_frac = my_loc%low_ix - frac

    ! Check bounds
    where (my_loc%low_ix < 1)
       my_loc%low_ix   = 1
       my_loc%low_frac = 1
    end where

    where (my_loc%low_ix >= my_lt%n_points - 1)
       my_loc%low_ix   = my_lt%n_points - 1
       my_loc%low_frac = 0
    end where
  end function LT2_get_loc

  !> Get the value of a single column at x
  pure function LT2_get_col(my_lt, col_ix, x1, x2) result(col_value)
    !$acc routine seq
    type(LT2_t), intent(in) :: my_lt
    integer, intent(in)     :: col_ix
    real(dp), intent(in)    :: x1, x2
    real(dp)                :: col_value
    type(LT2_loc_t)         :: loc

    loc       = LT2_get_loc(my_lt, x1, x2)
    col_value = LT2_get_col_at_loc(my_lt, col_ix, loc)
  end function LT2_get_col

  !> Get the value of a single column at a location
  pure function LT2_get_col_at_loc(my_lt, col_ix, loc) result(col_value)
    !$acc routine seq
    type(LT2_t), intent(in)     :: my_lt
    integer, intent(in)         :: col_ix
    type(LT2_loc_t), intent(in) :: loc
    integer                     :: ix(2)
    real(dp)                    :: w(2, 2)
    real(dp)                    :: col_value

    ! Bilinear interpolation
    w(1, 1) = loc%low_frac(1) * loc%low_frac(2)
    w(2, 1) = (1 - loc%low_frac(1)) * loc%low_frac(2)
    w(1, 2) = loc%low_frac(1) * (1 - loc%low_frac(2))
    w(2, 2) = (1 - loc%low_frac(1)) * (1 - loc%low_frac(2))
    ix = loc%low_ix

    col_value = sum(w * my_lt%rows_cols(ix(1):ix(1)+1, &
         ix(2):ix(2)+1, col_ix))
  end function LT2_get_col_at_loc

  !> Manually inlined version, a bit easier on the routine seq's
  elemental function LT2_get_col_inline(my_lt, col_ix, x1, x2) result(col_value)
    !$acc routine seq
    type(LT2_t), intent(in) :: my_lt
    integer,     intent(in) :: col_ix
    real(dp),    intent(in) :: x1, x2

    real(dp)        :: col_value
    real(dp)        :: w(2,2)
    real(dp)        :: frac(2)
    integer         :: ix(2)
    type(LT2_loc_t) :: loc

    frac = ([x1, x2] - my_lt%x_min) * my_lt%inv_fac
    loc%low_ix = ceiling(frac)
    loc%low_frac = loc%low_ix - frac

    ! Check bounds
    where (loc%low_ix < 1)
       loc%low_ix   = 1
       loc%low_frac = 1
    end where

    where (loc%low_ix >= my_lt%n_points - 1)
       loc%low_ix   = my_lt%n_points - 1
       loc%low_frac = 0
    end where

    ! Bilinear interpolation
    w(1, 1) = loc%low_frac(1) * loc%low_frac(2)
    w(2, 1) = (1 - loc%low_frac(1)) * loc%low_frac(2)
    w(1, 2) = loc%low_frac(1) * (1 - loc%low_frac(2))
    w(2, 2) = (1 - loc%low_frac(1)) * (1 - loc%low_frac(2))
    ix = loc%low_ix

    col_value = sum(w * my_lt%rows_cols(ix(1):ix(1)+1, &
         ix(2):ix(2)+1, col_ix))
  end function LT2_get_col_inline

end module modlookuptable