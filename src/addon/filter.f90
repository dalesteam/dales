subroutine filter(v2f, ndx, ghostx, ghosty)
  use modglobal, only : i1,ih,i2,j1,jh,j2
  ! top hat filter
  ! for now, only filter in horizontal
  implicit none

  integer, intent(in)    :: ndx, ghostx, ghosty
  real,    intent(inout) :: v2f(2-ghostx:i1+ghostx,2-ghosty:j1+ghosty)
  real                   :: v2fin(2-ghostx:i1+ghostx,2-ghosty:j1+ghosty)
  integer                :: i,j,m,n
  real                   :: weight
  integer                :: starti, startj, endi, endj

  v2fin(:,:) = v2f(:,:)
  v2f(:,:)   = 0.
  
  weight   = ndx**2.

  if(mod(ndx,2) == 0) then
    ! Filter width is even
    starti = 2  - ghostx + ndx / 2
    startj = 2  - ghosty + ndx / 2
    endi   = i1 + ghostx - ndx / 2
    endj   = j1 + ghosty - ndx / 2

    if((starti > 2 ) .or. (startj > 2)) stop 'ERROR: Filter larger than array'

    do j = startj, endj
      do i = starti, endi
        do m = - ndx / 2, ndx / 2
          do n = - ndx / 2, ndx / 2
            ! Check if we are on the left edge or right edge
            if(m == - ndx / 2 .or. m == ndx / 2) then
              ! Check if we are on the top edge or bottom edge
              if(n == - ndx / 2 .or. n == ndx / 2) then
                v2f(i,j) = v2f(i,j) + 0.25 * v2fin(i+m,j+n)
              else
                v2f(i,j) = v2f(i,j) + 0.5 * v2fin(i+m,j+n)
              end if
            ! We are not on left or right edge
            else
              ! Check if we are on top edge or bottom edge
              if(n == - ndx / 2 .or. n == ndx / 2) then
                v2f(i,j) = v2f(i,j) + 0.5 * v2fin(i+m,j+n)
              ! We are not on any edge
              else
                v2f(i,j) = v2f(i,j) + 1.0 * v2fin(i+m,j+n)
              end if
            end if
          end do
        end do
        v2f(i,j) = v2f(i,j) / weight
      end do
    end do

  else
    ! Filter width is odd
    starti = 2  - ghostx + (ndx - 1) / 2
    startj = 2  - ghosty + (ndx - 1) / 2
    endi   = i1 + ghostx - (ndx - 1) / 2
    endj   = j1 + ghosty - (ndx - 1) / 2
    
    if((starti > 2 ) .or. (startj > 2)) stop 'ERROR: Filter larger than array'

    do j = startj, endj
      do i = starti, endi
        do m = - (ndx - 1) / 2, (ndx - 1) / 2
          do n = - (ndx - 1) / 2, (ndx - 1) / 2
             v2f(i,j) = v2f(i,j) + v2fin(i+m,j+n)
          end do
        end do
        v2f(i,j) = v2f(i,j) / weight
      end do
    end do
  end if

  return
end subroutine filter
