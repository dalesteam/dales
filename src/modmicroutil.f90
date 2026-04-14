! This file is part of DALES.
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
! Copyright 1993-2026, The DALES Team.
!

!> Common utility routines for microphysics.
module modmicroutil

  use modprecision, only: field_r

  implicit none

  public

contains

  subroutine zero_field(field)

    real(field_r), intent(out) :: field(:,:,:)

    integer :: i, j, k

    !$acc parallel loop collapse(3) default(present)
!!$omp target teams loop collapse(3) defaultmap(present:aggregate)&
!!$omp defaultmap(present:allocatable)
    do k = 1, size(field, dim=3)
      do j = 1, size(field, dim=2)
        do i = 1, size(field, dim=1)
          field(i,j,k) = 0
        end do
      end do
    end do

  end subroutine zero_field

  subroutine sum_fields(field, other)

    real(field_r), intent(in)  :: field(:,:,:)
    real(field_r), intent(out) :: other(:,:,:)

    integer :: i, j, k

    !$acc parallel loop collapse(3) default(present)
!!$omp target teams loop collapse(3) defaultmap(present:aggregate)&
!!$omp defaultmap(present:allocatable)
    do k = 1, size(field, dim=3)
      do j = 1, size(field, dim=2)
        do i = 1, size(field, dim=1)
          other(i,j,k) = other(i,j,k) + field(i,j,k)
        end do
      end do
    end do

  end subroutine sum_fields

end module modmicroutil
