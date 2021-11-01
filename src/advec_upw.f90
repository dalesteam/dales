!> \file advec_upw.f90
!!  Does advection with a 1st order upwind scheme.
!! \par Revision list
!! \par Authors
!  This file is part of DALES.
!
! DALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!

module advec_upw
use modprecision, only : field_r
contains
!> Advection at cell center
subroutine advecc_upw(putin,putout)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi,dyi,dzi
  use modfields, only : u0, v0, w0, rhobf
  implicit none

  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: putin !< Input: the cell centered field
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: putout !< Output: the tendency

  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1) :: put
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1) :: rhoputin
  integer :: i,j,k


  do k=1,k1
    do j=2-jh,j1+jh
      do i=2-ih,i1+ih
      rhoputin(i,j,k)=rhobf(k)*putin(i,j,k)
      end do
    end do
  end do

  do k=1,k1
    do j=2,j1
      do i=2,i1+1
       if( u0(i,j,k) > 0 ) then
           put(i,j,k) = rhoputin(i-1,j,k)
       else
           put(i,j,k) = rhoputin(i,j,k)
       endif
      enddo
    enddo
  enddo

  do k=1,k1
    do j=2,j1
      do i=2,i1
        putout(i,j,k) = putout(i,j,k) - &
            (1./rhobf(k))*(u0(i+1,j,k)*put(i+1,j,k)-u0(i,j,k)*put(i,j,k))*dxi
      enddo
    enddo
  enddo

  do k=1,k1
    do j=2,j1+1
      do i=2,i1
       if( v0(i,j,k) > 0 ) then
           put(i,j,k) = rhoputin(i,j-1,k)
       else
           put(i,j,k) = rhoputin(i,j,k)
       endif
      enddo
    enddo
  enddo
  do k=1,k1
    do j=2,j1
      do i=2,i1
        putout(i,j,k) = putout(i,j,k) - &
            (1./rhobf(k))*(v0(i,j+1,k)*put(i,j+1,k)-v0(i,j,k)*put(i,j,k))*dyi
      enddo
    enddo
  enddo

  put(2:i1,2:j1, 1) = 0
  put(2:i1,2:j1,k1) = 0
  do k=2,kmax
    do j=2,j1
      do i=2,i1
       if( w0(i,j,k) > 0 ) then
           put(i,j,k) = rhoputin(i,j,k-1)
       else
           put(i,j,k) = rhoputin(i,j,k)
       endif
      enddo
    enddo
  enddo
  do k=1,kmax
    do j=2,j1
      do i=2,i1
        putout(i,j,k) = putout(i,j,k) - &
            (1./rhobf(k))*(w0(i,j,k+1)*put(i,j,k+1)-w0(i,j,k)*put(i,j,k))*dzi
      enddo
    enddo
  enddo

end subroutine advecc_upw
end module advec_upw
