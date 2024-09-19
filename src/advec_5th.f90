!> \file advec_5th.f90
!!  Does advection with a 5th order upwind scheme.
!! \par Revision list
!! \par Authors
!! \see Wicker and Scamarock 2002
!!
!! By adding a small dissipative term to the sixth order flux, a fifth order
!! scheme is created that is nearly monotone:
!! \latexonly
!! \begin{eqnarray}
!!  F_{i-\frac{1}{2}}^{5th} &=& F_{i-\frac{1}{2}}^{6th} -
!! \left|\frac{\fav{u}_{i-\frac{1}{2}}}{60}\right|\left[10(\phi_i-\phi_{i-1})\right
!! . \nonumber\\\\
!! &&-\left.5(\phi_{i+1}-\phi_{i-2})+(\phi_{i+2}-\phi_{i-3})\right].
!! \end{eqnarray}
!! \endlatexonly
!!
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

module advec_5th
use modprecision, only : field_r
contains
!> Horizontal advection at cell center
subroutine hadvecc_5th(a_in, a_out)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi,dyi
  use modfields, only : u0, v0

  implicit none

  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: a_in !< Input: the cell centered field
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: a_out !< Output: the tendency

  integer :: i,j,k

  !$acc parallel loop collapse(3) default(present)
  do k = 1, kmax
    do j = 2, j1
      do i = 2, i1
        a_out(i,j,k)  = a_out(i,j,k) - (  &
              ( &
                  u0(i+1,j,k)/60&
                  *(37*(a_in(i+1,j,k)+a_in(i,j,k))-8*(a_in(i+2,j,k)+a_in(i-1,j,k))+(a_in(i+3,j,k)+a_in(i-2,j,k)))&
                  -sign(1._field_r,u0(i+1,j,k))*u0(i+1,j,k)/60&
                  *(10*(a_in(i+1,j,k)-a_in(i,j,k))-5*(a_in(i+2,j,k)-a_in(i-1,j,k))+(a_in(i+3,j,k)-a_in(i-2,j,k)))&
                  -u0(i,j,k)/60&
                  *(37*(a_in(i,j,k)+a_in(i-1,j,k))-8*(a_in(i+1,j,k)+a_in(i-2,j,k))+(a_in(i+2,j,k)+a_in(i-3,j,k)))&
                  +sign(1._field_r,u0(i,j,k))*u0(i,j,k)/60&
                  *(10*(a_in(i,j,k)-a_in(i-1,j,k))-5*(a_in(i+1,j,k)-a_in(i-2,j,k))+(a_in(i+2,j,k)-a_in(i-3,j,k)))&
              )* dxi &
            +( &
                  v0(i,j+1,k)/60&
                  *(37*(a_in(i,j+1,k)+a_in(i,j,k))-8*(a_in(i,j+2,k)+a_in(i,j-1,k))+(a_in(i,j+3,k)+a_in(i,j-2,k)))&
                  -sign(1._field_r,v0(i,j+1,k))*v0(i,j+1,k)/60&
                  *(10*(a_in(i,j+1,k)-a_in(i,j,k))-5*(a_in(i,j+2,k)-a_in(i,j-1,k))+(a_in(i,j+3,k)-a_in(i,j-2,k)))&
                  -v0(i,j,k)/60&
                  *(37*(a_in(i,j,k)+a_in(i,j-1,k))-8*(a_in(i,j+1,k)+a_in(i,j-2,k))+(a_in(i,j+2,k)+a_in(i,j-3,k)))&
                  +sign(1._field_r,v0(i,j,k))*v0(i,j,k)/60&
                  *(10*(a_in(i,j,k)-a_in(i,j-1,k))-5*(a_in(i,j+1,k)-a_in(i,j-2,k))+(a_in(i,j+2,k)-a_in(i,j-3,k)))&
              )* dyi &
              )
      end do
    end do
  end do
end subroutine hadvecc_5th

!> Vertical advection at cell center
subroutine vadvecc_5th(a_in, a_out)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dzfi
  use modfields, only : w0, rhobf

  implicit none

  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: a_in !< Input: the cell centered field
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: a_out !< Output: the tendency

  integer :: i,j,k

  k = 1
  !$acc parallel loop collapse(2) default(present) async(1)
  do j = 2, j1
    do i = 2, i1
      a_out(i,j,k) = a_out(i,j,k) - ( (1/rhobf(1)) * ( &
                    w0(i,j,k+1) * (rhobf(k+1)*a_in(i,j,k+1) + rhobf(k)*a_in(i,j,k)) &
              ) * ( 0.5_field_r * dzfi(k) ) &
            )
    end do
  end do

  ! CvH do 2nd order for bottom and top
  k = 2
  !$acc parallel loop collapse(2) default(present) async(1)
  do j = 2, j1
    do i = 2, i1
      a_out(i,j,k) = a_out(i,j,k) - ( (1/rhobf(k)) * ( &
                    w0(i,j,k+1) * (rhobf(k+1)*a_in(i,j,k+1)+rhobf(k)*a_in(i,j,k)) &
                  - w0(i,j,k  ) * (rhobf(k-1)*a_in(i,j,k-1)+rhobf(k)*a_in(i,j,k)) &
              ) * ( 0.5_field_r * dzfi(k) ) &
            )
    end do
  end do

  k = kmax-1
  !$acc parallel loop collapse(2) default(present) async(1)
  do j = 2, j1
    do i = 2, i1
      a_out(i,j,k) = a_out(i,j,k) - ( (1/rhobf(k)) * ( &
                    w0(i,j,k+1) * (rhobf(k+1)*a_in(i,j,k+1)+rhobf(k)*a_in(i,j,k)) &
                  - w0(i,j,k  ) * (rhobf(k-1)*a_in(i,j,k-1)+rhobf(k)*a_in(i,j,k)) &
              ) * ( 0.5_field_r * dzfi(k) ) &
            )
    end do
  end do

  k = kmax
  !$acc parallel loop collapse(2) default(present) async(1)
  do j = 2, j1
    do i = 2, i1
      a_out(i,j,k) = a_out(i,j,k) - ( (1/rhobf(k)) * ( &
                    w0(i,j,k+1) * (rhobf(k+1)*a_in(i,j,k+1)+rhobf(k)*a_in(i,j,k)) &
                  - w0(i,j,k  ) * (rhobf(k-1)*a_in(i,j,k-1)+rhobf(k)*a_in(i,j,k)) &
              ) * ( 0.5_field_r * dzfi(k) ) &
            )
    end do
  end do

  k = 3
  !$acc parallel loop collapse(2) default(present) async(1)
  do j = 2, j1
    do i = 2, i1
      a_out(i,j,k) = a_out(i,j,k) - ( (1/rhobf(k)) * ( &
                    w0(i,j,k+1) / 60 &
                    * ( 37*(rhobf(k+1)*a_in(i,j,k+1)+rhobf(k  )*a_in(i,j,k  )) &
                       - 8*(rhobf(k+2)*a_in(i,j,k+2)+rhobf(k-1)*a_in(i,j,k-1)) &
                       +   (rhobf(k+3)*a_in(i,j,k+3)+rhobf(k-2)*a_in(i,j,k-2)) ) &
                  - sign(1._field_r,w0(i,j,k+1))*w0(i,j,k+1) / 60 &
                    * ( 10*(rhobf(k+1)*a_in(i,j,k+1)-rhobf(k  )*a_in(i,j,k  )) &
                       - 5*(rhobf(k+2)*a_in(i,j,k+2)-rhobf(k-1)*a_in(i,j,k-1)) &
                       +   (rhobf(k+3)*a_in(i,j,k+3)-rhobf(k-2)*a_in(i,j,k-2)) ) &
                  - w0(i,j,k) / 2 &
                    *      (rhobf(k-1)*a_in(i,j,k-1)+rhobf(k  )*a_in(i,j,k  )) &
              ) * dzfi(k) &
            )
    end do
  end do

  k = kmax-2
  !$acc parallel loop collapse(2) default(present) async(1)
  do j = 2, j1
    do i = 2, i1
      a_out(i,j,k) = a_out(i,j,k) - ( (1/rhobf(k)) * ( &
                    w0(i,j,k+1) / 2 &
                    *      (rhobf(k+1)*a_in(i,j,k+1)+rhobf(k  )*a_in(i,j,k  )) &
                  - w0(i,j,k  ) / 60 &
                    * ( 37*(rhobf(k  )*a_in(i,j,k  )+rhobf(k-1)*a_in(i,j,k-1)) &
                       - 8*(rhobf(k+1)*a_in(i,j,k+1)+rhobf(k-2)*a_in(i,j,k-2)) &
                       +   (rhobf(k+2)*a_in(i,j,k+2)+rhobf(k-3)*a_in(i,j,k-3)) ) &
                  + sign(1._field_r,w0(i,j,k))*w0(i,j,k) / 60 &
                    * ( 10*(rhobf(k  )*a_in(i,j,k  )-rhobf(k-1)*a_in(i,j,k-1)) &
                       - 5*(rhobf(k+1)*a_in(i,j,k+1)-rhobf(k-2)*a_in(i,j,k-2)) &
                       +   (rhobf(k+2)*a_in(i,j,k+2)-rhobf(k-3)*a_in(i,j,k-3)) ) &
              ) * dzfi(k) &
            )
    end do
  end do


  !$acc parallel loop collapse(3) default(present) async(2)
  do k = 4, kmax-3
    do j = 2, j1
      do i = 2, i1
        a_out(i,j,k) = a_out(i,j,k) - ( (1/rhobf(k)) * ( &
                      w0(i,j,k+1) / 60 &
                      * ( 37*(rhobf(k+1)*a_in(i,j,k+1)+rhobf(k  )*a_in(i,j,k  )) &
                         - 8*(rhobf(k+2)*a_in(i,j,k+2)+rhobf(k-1)*a_in(i,j,k-1)) &
                         +   (rhobf(k+3)*a_in(i,j,k+3)+rhobf(k-2)*a_in(i,j,k-2)) ) &
                    - sign(1._field_r,w0(i,j,k+1))*w0(i,j,k+1) / 60 &
                      * ( 10*(rhobf(k+1)*a_in(i,j,k+1)-rhobf(k  )*a_in(i,j,k  )) &
                         - 5*(rhobf(k+2)*a_in(i,j,k+2)-rhobf(k-1)*a_in(i,j,k-1)) &
                         +   (rhobf(k+3)*a_in(i,j,k+3)-rhobf(k-2)*a_in(i,j,k-2)) ) &
                    - w0(i,j,k) / 60 &
                      * ( 37*(rhobf(k  )*a_in(i,j,k  )+rhobf(k-1)*a_in(i,j,k-1)) &
                         - 8*(rhobf(k+1)*a_in(i,j,k+1)+rhobf(k-2)*a_in(i,j,k-2)) &
                         +   (rhobf(k+2)*a_in(i,j,k+2)+rhobf(k-3)*a_in(i,j,k-3)) ) &
                    + sign(1._field_r,w0(i,j,k))*w0(i,j,k) / 60 &
                      * ( 10*(rhobf(k  )*a_in(i,j,k  )-rhobf(k-1)*a_in(i,j,k-1)) &
                         - 5*(rhobf(k+1)*a_in(i,j,k+1)-rhobf(k-2)*a_in(i,j,k-2)) &
                         +   (rhobf(k+2)*a_in(i,j,k+2)-rhobf(k-3)*a_in(i,j,k-3)) ) &
                ) * dzfi(k) &
              )
      end do
    end do
  end do
  !$acc wait
end subroutine vadvecc_5th

!> Horizontal advection at the u point.
subroutine hadvecu_5th(a_in,a_out)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi5,dyi5
  use modfields, only : u0, v0
  implicit none

  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: a_in !< Input: the u field
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: a_out !< Output: the tendency

  integer :: i,j,k

  !$acc parallel loop collapse(3) default(present) async(1)
  do k = 1, kmax
    do j = 2, j1
      do i = 2, i1
        a_out(i,j,k)  = a_out(i,j,k)- ( &
              ( &
                  (u0(i+1,j,k)+u0(i,j,k))/60&
                  *(37*(a_in(i+1,j,k)+a_in(i,j,k))-8*(a_in(i+2,j,k)+a_in(i-1,j,k))+(a_in(i+3,j,k)+a_in(i-2,j,k)))&
                  -sign(1._field_r,(u0(i+1,j,k)+u0(i,j,k)))*(u0(i+1,j,k)+u0(i,j,k))/60&
                  *(10*(a_in(i+1,j,k)-a_in(i,j,k))-5*(a_in(i+2,j,k)-a_in(i-1,j,k))+(a_in(i+3,j,k)-a_in(i-2,j,k)))&
                  -(u0(i,j,k)+u0(i-1,j,k))/60&
                  *(37*(a_in(i,j,k)+a_in(i-1,j,k))-8*(a_in(i+1,j,k)+a_in(i-2,j,k))+(a_in(i+2,j,k)+a_in(i-3,j,k)))&
                  +sign(1._field_r,(u0(i,j,k)+u0(i-1,j,k)))*(u0(i,j,k)+u0(i-1,j,k))/60&
                  *(10*(a_in(i,j,k)-a_in(i-1,j,k))-5*(a_in(i+1,j,k)-a_in(i-2,j,k))+(a_in(i+2,j,k)-a_in(i-3,j,k)))&
              )*dxi5&
            +(&
                  (v0(i,j+1,k)+v0(i-1,j+1,k))/60&
                  *(37*(a_in(i,j+1,k)+a_in(i,j,k))-8*(a_in(i,j+2,k)+a_in(i,j-1,k))+(a_in(i,j+3,k)+a_in(i,j-2,k)))&
                  -sign(1._field_r,(v0(i,j+1,k)+v0(i-1,j+1,k)))*(v0(i,j+1,k)+v0(i-1,j+1,k))/60&
                  *(10*(a_in(i,j+1,k)-a_in(i,j,k))-5*(a_in(i,j+2,k)-a_in(i,j-1,k))+(a_in(i,j+3,k)-a_in(i,j-2,k)))&
                  -(v0(i,j,k)+v0(i-1,j,k))/60&
                  *(37*(a_in(i,j,k)+a_in(i,j-1,k))-8*(a_in(i,j+1,k)+a_in(i,j-2,k))+(a_in(i,j+2,k)+a_in(i,j-3,k)))&
                  +sign(1._field_r,(v0(i,j,k)+v0(i-1,j,k)))*(v0(i,j,k)+v0(i-1,j,k))/60&
                  *(10*(a_in(i,j,k)-a_in(i,j-1,k))-5*(a_in(i,j+1,k)-a_in(i,j-2,k))+(a_in(i,j+2,k)-a_in(i,j-3,k)))&
              )* dyi5 &
              )
      end do
    end do
  end do
end subroutine hadvecu_5th

!> Vertical advection at the u point.
subroutine vadvecu_5th(a_in,a_out)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dzfi
  use modfields, only : w0, rhobf
  implicit none

  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: a_in !< Input: the u field
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: a_out !< Output: the tendency

  integer :: i,j,k

  k = 1
  !$acc parallel loop collapse(2) default(present) async(1)
  do j = 2, j1
    do i = 2, i1
      a_out(i,j,k) = a_out(i,j,k) - ( &
          (1/rhobf(1))*( &
            ( rhobf(k+1)*a_in(i,j,k+1) + rhobf(k)*a_in(i,j,k)) *(w0(i,j,k+1)+ w0(i-1,j,k+1)) &
            ) * (0.25_field_r * dzfi(k)) &
            )
    end do
  end do

  k = 2
  !$acc parallel loop collapse(2) default(present) async(1)
  do j = 2, j1
    do i = 2, i1
      a_out(i,j,k) = a_out(i,j,k) - ( &
          (1/rhobf(k))*( &
            (rhobf(k)*a_in(i,j,k)+rhobf(k+1)*a_in(i,j,k+1) )*(w0(i,j,k+1)+w0(i-1,j,k+1)) &
          - (rhobf(k)*a_in(i,j,k)+rhobf(k-1)*a_in(i,j,k-1) )*(w0(i,j,k  )+w0(i-1,j,k  )) &
            ) * (0.25_field_r * dzfi(k)) &
            )
    end do
  end do

  k = kmax-1
  !$acc parallel loop collapse(2) default(present) async(1)
  do j = 2, j1
    do i = 2, i1
      a_out(i,j,k) = a_out(i,j,k) - ( &
          (1/rhobf(k))*( &
            (rhobf(k)*a_in(i,j,k)+rhobf(k+1)*a_in(i,j,k+1) )*(w0(i,j,k+1)+w0(i-1,j,k+1)) &
          - (rhobf(k)*a_in(i,j,k)+rhobf(k-1)*a_in(i,j,k-1) )*(w0(i,j,k  )+w0(i-1,j,k  )) &
            ) * (0.25_field_r * dzfi(k)) &
            )
    end do
  end do

  k = kmax
  !$acc parallel loop collapse(2) default(present) async(1)
  do j = 2, j1
    do i = 2, i1
      a_out(i,j,k) = a_out(i,j,k) - ( &
          (1/rhobf(k))*( &
            (rhobf(k)*a_in(i,j,k)+rhobf(k+1)*a_in(i,j,k+1) )*(w0(i,j,k+1)+w0(i-1,j,k+1)) &
          - (rhobf(k)*a_in(i,j,k)+rhobf(k-1)*a_in(i,j,k-1) )*(w0(i,j,k  )+w0(i-1,j,k  )) &
            ) * (0.25_field_r * dzfi(k)) &
            )
    end do
  end do

  k = 3
  !$acc parallel loop collapse(2) default(present) async(1)
  do j = 2, j1
    do i = 2, i1
      a_out(i,j,k) = a_out(i,j,k) - ( (1/rhobf(k)) * ( &
                    ( w0(i,j,k+1) + w0(i-1,j,k+1) ) / 60 &
                    * ( 37*(rhobf(k+1)*a_in(i,j,k+1)+rhobf(k  )*a_in(i,j,k  )) &
                       - 8*(rhobf(k+2)*a_in(i,j,k+2)+rhobf(k-1)*a_in(i,j,k-1)) &
                       +   (rhobf(k+3)*a_in(i,j,k+3)+rhobf(k-2)*a_in(i,j,k-2)) )&
                  - sign(1._field_r,(w0(i,j,k+1)+w0(i-1,j,k+1)))*(w0(i,j,k+1)+w0(i-1,j,k+1)) / 60 &
                    * ( 10*(rhobf(k+1)*a_in(i,j,k+1)-rhobf(k  )*a_in(i,j,k  )) &
                       - 5*(rhobf(k+2)*a_in(i,j,k+2)-rhobf(k-1)*a_in(i,j,k-1)) &
                       +   (rhobf(k+3)*a_in(i,j,k+3)-rhobf(k-2)*a_in(i,j,k-2)) )&
                  - ( w0(i,j,k  ) + w0(i-1,j,k  ) ) / 2 &
                    *      (rhobf(k  )*a_in(i,j,k  )+rhobf(k-1)*a_in(i,j,k-1) ) &
            ) * (0.5_field_r * dzfi(k)) &
            )
    end do
  end do

  k = kmax-2
  !$acc parallel loop collapse(2) default(present) async(1)
  do j = 2, j1
    do i = 2, i1
      a_out(i,j,k) = a_out(i,j,k) - ( (1/rhobf(k)) * ( &
                    ( w0(i,j,k+1) + w0(i-1,j,k+1)) / 2 &
                    * ( rhobf(k)*a_in(i,j,k)+rhobf(k+1)*a_in(i,j,k+1) ) &
                  - ( w0(i,j,k  ) + w0(i-1,j,k  )) / 60 &
                    * ( 37*(rhobf(k  )*a_in(i,j,k  )+rhobf(k-1)*a_in(i,j,k-1)) &
                       - 8*(rhobf(k+1)*a_in(i,j,k+1)+rhobf(k-2)*a_in(i,j,k-2)) &
                       +   (rhobf(k+2)*a_in(i,j,k+2)+rhobf(k-3)*a_in(i,j,k-3)) ) &
                  + sign(1._field_r,(w0(i,j,k)+w0(i-1,j,k)))*(w0(i,j,k)+w0(i-1,j,k)) / 60 &
                    * ( 10*(rhobf(k  )*a_in(i,j,k  )-rhobf(k-1)*a_in(i,j,k-1)) &
                       - 5*(rhobf(k+1)*a_in(i,j,k+1)-rhobf(k-2)*a_in(i,j,k-2)) &
                       +   (rhobf(k+2)*a_in(i,j,k+2)-rhobf(k-3)*a_in(i,j,k-3)) ) &
            ) * (0.5_field_r * dzfi(k)) &
            )
    end do
  end do

  !$acc parallel loop collapse(3) default(present) async(2)
  do k = 4, kmax-3
    do j = 2, j1
      do i = 2, i1
        a_out(i,j,k) = a_out(i,j,k) - ( (1/rhobf(k)) * ( &
                      ( w0(i,j,k+1) + w0(i-1,j,k+1) ) / 60 &
                      * ( 37*(rhobf(k+1)*a_in(i,j,k+1)+rhobf(k  )*a_in(i,j,k  )) &
                         - 8*(rhobf(k+2)*a_in(i,j,k+2)+rhobf(k-1)*a_in(i,j,k-1)) &
                         +   (rhobf(k+3)*a_in(i,j,k+3)+rhobf(k-2)*a_in(i,j,k-2)) ) &
                    - sign(1._field_r,(w0(i,j,k+1)+w0(i-1,j,k+1))) * (w0(i,j,k+1)+w0(i-1,j,k+1)) / 60 &
                      * ( 10*(rhobf(k+1)*a_in(i,j,k+1)-rhobf(k  )*a_in(i,j,k  )) &
                         - 5*(rhobf(k+2)*a_in(i,j,k+2)-rhobf(k-1)*a_in(i,j,k-1)) &
                         +   (rhobf(k+3)*a_in(i,j,k+3)-rhobf(k-2)*a_in(i,j,k-2)) ) &
                    - ( w0(i,j,k  ) + w0(i-1,j,k) ) / 60 &
                      * ( 37*(rhobf(k  )*a_in(i,j,k  )+rhobf(k-1)*a_in(i,j,k-1)) &
                         - 8*(rhobf(k+1)*a_in(i,j,k+1)+rhobf(k-2)*a_in(i,j,k-2)) &
                         +   (rhobf(k+2)*a_in(i,j,k+2)+rhobf(k-3)*a_in(i,j,k-3)) ) &
                    + sign(1._field_r,(w0(i,j,k  )+w0(i-1,j,k  ))) * (w0(i,j,k  )+w0(i-1,j,k  )) / 60 &
                      * ( 10*(rhobf(k  )*a_in(i,j,k  )-rhobf(k-1)*a_in(i,j,k-1)) &
                         - 5*(rhobf(k+1)*a_in(i,j,k+1)-rhobf(k-2)*a_in(i,j,k-2)) &
                         +   (rhobf(k+2)*a_in(i,j,k+2)-rhobf(k-3)*a_in(i,j,k-3)) ) &
                ) * (0.5_field_r * dzfi(k)) &
                )
      end do
    end do
  end do
end subroutine vadvecu_5th

!> Horizontal advection at the v point.
subroutine hadvecv_5th(a_in, a_out)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi5,dyi5
  use modfields, only : u0, v0
  implicit none

  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: a_in !< Input: the v field
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: a_out !< Output: the tendency
  integer :: i,j,k

  !$acc parallel loop collapse(3) default(present) async(2)
  do k = 1, kmax
    do j = 2, j1
      do i = 2, i1
        a_out(i,j,k)  = a_out(i,j,k) - ( &
              ( &
                  (u0(i+1,j,k)+u0(i+1,j-1,k))/60&
                  *(37*(a_in(i+1,j,k)+a_in(i,j,k))-8*(a_in(i+2,j,k)+a_in(i-1,j,k))+(a_in(i+3,j,k)+a_in(i-2,j,k)))&
                  -sign(1._field_r,(u0(i+1,j,k)+u0(i+1,j-1,k)))*(u0(i+1,j,k)+u0(i+1,j-1,k))/60&
                  *(10*(a_in(i+1,j,k)-a_in(i,j,k))-5*(a_in(i+2,j,k)-a_in(i-1,j,k))+(a_in(i+3,j,k)-a_in(i-2,j,k)))&
                  -(u0(i,j,k)+u0(i,j-1,k))/60&
                  *(37*(a_in(i,j,k)+a_in(i-1,j,k))-8*(a_in(i+1,j,k)+a_in(i-2,j,k))+(a_in(i+2,j,k)+a_in(i-3,j,k)))&
                  +sign(1._field_r,(u0(i,j,k)+u0(i,j-1,k)))*(u0(i,j,k)+u0(i,j-1,k))/60&
                  *(10*(a_in(i,j,k)-a_in(i-1,j,k))-5*(a_in(i+1,j,k)-a_in(i-2,j,k))+(a_in(i+2,j,k)-a_in(i-3,j,k)))&
                )*dxi5&
              +(&
                  (v0(i,j+1,k)+v0(i,j,k))/60&
                  *(37*(a_in(i,j+1,k)+a_in(i,j,k))-8*(a_in(i,j+2,k)+a_in(i,j-1,k))+(a_in(i,j+3,k)+a_in(i,j-2,k)))&
                  -sign(1._field_r,(v0(i,j+1,k)+v0(i,j,k)))*(v0(i,j+1,k)+v0(i,j,k))/60&
                  *(10*(a_in(i,j+1,k)-a_in(i,j,k))-5*(a_in(i,j+2,k)-a_in(i,j-1,k))+(a_in(i,j+3,k)-a_in(i,j-2,k)))&
                  -(v0(i,j,k)+v0(i,j-1,k))/60&
                  *(37*(a_in(i,j,k)+a_in(i,j-1,k))-8*(a_in(i,j+1,k)+a_in(i,j-2,k))+(a_in(i,j+2,k)+a_in(i,j-3,k)))&
                  +sign(1._field_r,(v0(i,j,k)+v0(i,j-1,k)))*(v0(i,j,k)+v0(i,j-1,k))/60&
                  *(10*(a_in(i,j,k)-a_in(i,j-1,k))-5*(a_in(i,j+1,k)-a_in(i,j-2,k))+(a_in(i,j+2,k)-a_in(i,j-3,k)))&
                )* dyi5 &
                )
      end do
    end do
  end do

end subroutine hadvecv_5th

!> Vertical advection at the v point.
subroutine vadvecv_5th(a_in, a_out)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dzfi
  use modfields, only : w0, rhobf
  implicit none

  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: a_in !< Input: the v field
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: a_out !< Output: the tendency

  integer :: i,j,k

  k = 1
  !$acc parallel loop collapse(2) default(present) async(3)
  do j = 2, j1
    do i = 2, i1
      a_out(i,j,k) = a_out(i,j,k) - ( &
            (1/rhobf(1))*( &
              (w0(i,j,k+1)+w0(i,j-1,k+1)) *(rhobf(k+1)*a_in(i,j,k+1)+rhobf(k)*a_in(i,j,k)) &
              ) * (0.25_field_r * dzfi(k)) &
            )
    end do
  end do

  k = 2
  !$acc parallel loop collapse(2) default(present) async(3)
  do j = 2, j1
    do i = 2, i1
      a_out(i,j,k) = a_out(i,j,k) - ( &
            (1/rhobf(k))*( &
              (w0(i,j,k+1)+w0(i,j-1,k+1))*(rhobf(k+1)*a_in(i,j,k+1)+rhobf(k)*a_in(i,j,k)) &
             -(w0(i,j,k  )+w0(i,j-1,k  ))*(rhobf(k-1)*a_in(i,j,k-1)+rhobf(k)*a_in(i,j,k)) &
              ) * (0.25_field_r * dzfi(k)) &
            )
    end do
  end do

  k = kmax-1
  !$acc parallel loop collapse(2) default(present) async(3)
  do j = 2, j1
    do i = 2, i1
      a_out(i,j,k) = a_out(i,j,k) - ( &
            (1/rhobf(k))*( &
              (w0(i,j,k+1)+w0(i,j-1,k+1))*(rhobf(k+1)*a_in(i,j,k+1)+rhobf(k)*a_in(i,j,k)) &
             -(w0(i,j,k  )+w0(i,j-1,k  ))*(rhobf(k-1)*a_in(i,j,k-1)+rhobf(k)*a_in(i,j,k)) &
              ) * (0.25_field_r * dzfi(k)) &
            )
    end do
  end do

  k = kmax
  !$acc parallel loop collapse(2) default(present) async(3)
  do j = 2, j1
    do i = 2, i1
      a_out(i,j,k) = a_out(i,j,k) - ( &
            (1/rhobf(k))*( &
              (w0(i,j,k+1)+w0(i,j-1,k+1))*(rhobf(k+1)*a_in(i,j,k+1)+rhobf(k)*a_in(i,j,k)) &
             -(w0(i,j,k  )+w0(i,j-1,k  ))*(rhobf(k-1)*a_in(i,j,k-1)+rhobf(k)*a_in(i,j,k)) &
              ) * (0.25_field_r * dzfi(k)) &
            )
    end do
  end do

  k = 3
  !$acc parallel loop collapse(2) default(present) async(3)
  do j = 2, j1
    do i = 2, i1
      a_out(i,j,k) = a_out(i,j,k)- ( (1/rhobf(k)) * ( &
                    ( w0(i,j,k+1) + w0(i,j-1,k+1) ) / 60 &
                    * ( 37*(rhobf(k+1)*a_in(i,j,k+1)+rhobf(k  )*a_in(i,j,k  )) &
                       - 8*(rhobf(k+2)*a_in(i,j,k+2)+rhobf(k-1)*a_in(i,j,k-1)) &
                       +   (rhobf(k+3)*a_in(i,j,k+3)+rhobf(k-2)*a_in(i,j,k-2)) ) &
                  - sign(1._field_r,(w0(i,j,k+1)+w0(i,j-1,k+1)))*(w0(i,j,k+1)+w0(i,j-1,k+1)) / 60 &
                    * ( 10*(rhobf(k+1)*a_in(i,j,k+1)-rhobf(k  )*a_in(i,j,k  )) &
                       - 5*(rhobf(k+2)*a_in(i,j,k+2)-rhobf(k-1)*a_in(i,j,k-1)) &
                       +   (rhobf(k+3)*a_in(i,j,k+3)-rhobf(k-2)*a_in(i,j,k-2)) ) &
                  - ( w0(i,j,k  ) + w0(i,j-1,k  ) ) / 2 &
                    * (rhobf(k-1)*a_in(i,j,k-1)+rhobf(k)*a_in(i,j,k)) &
              ) * (0.5_field_r * dzfi(k)) &
            )
    end do
  end do

  k = kmax-2
  !$acc parallel loop collapse(2) default(present) async(3)
  do j = 2, j1
    do i = 2, i1
      a_out(i,j,k) = a_out(i,j,k) - ( (1/rhobf(k)) * ( &
                    ( w0(i,j,k+1) + w0(i,j-1,k+1) ) / 2 &
                     * (rhobf(k+1)*a_in(i,j,k+1)+rhobf(k)*a_in(i,j,k)) &
                  - ( w0(i,j,k  ) + w0(i,j-1,k  ) ) / 60 &
                     * ( 37*(rhobf(k  )*a_in(i,j,k  )+rhobf(k-1)*a_in(i,j,k-1)) &
                        - 8*(rhobf(k+1)*a_in(i,j,k+1)+rhobf(k-2)*a_in(i,j,k-2)) &
                        +   (rhobf(k+2)*a_in(i,j,k+2)+rhobf(k-3)*a_in(i,j,k-3)) ) &
                  + sign(1._field_r,(w0(i,j,k)+w0(i,j-1,k)))*(w0(i,j,k)+w0(i,j-1,k)) / 60 &
                     * ( 10*(rhobf(k  )*a_in(i,j,k  )-rhobf(k-1)*a_in(i,j,k-1)) &
                        - 5*(rhobf(k+1)*a_in(i,j,k+1)-rhobf(k-2)*a_in(i,j,k-2)) &
                        +   (rhobf(k+2)*a_in(i,j,k+2)-rhobf(k-3)*a_in(i,j,k-3)) ) &
              ) * (0.5_field_r * dzfi(k)) &
            )
    end do
  end do

  !$acc parallel loop collapse(3) default(present) async(4)
  do k = 4, kmax-3
    do j = 2, j1
      do i = 2, i1
        a_out(i,j,k) = a_out(i,j,k) - ( (1/rhobf(k)) * ( &
                      ( w0(i,j,k+1) + w0(i,j-1,k+1) ) / 60 &
                      * ( 37*(rhobf(k+1)*a_in(i,j,k+1)+rhobf(k  )*a_in(i,j,k  )) &
                         - 8*(rhobf(k+2)*a_in(i,j,k+2)+rhobf(k-1)*a_in(i,j,k-1)) &
                         +   (rhobf(k+3)*a_in(i,j,k+3)+rhobf(k-2)*a_in(i,j,k-2)) ) &
                    - sign(1._field_r,(w0(i,j,k+1)+w0(i,j-1,k+1)))*(w0(i,j,k+1)+w0(i,j-1,k+1)) / 60 &
                      * ( 10*(rhobf(k+1)*a_in(i,j,k+1)-rhobf(k  )*a_in(i,j,k  )) &
                         - 5*(rhobf(k+2)*a_in(i,j,k+2)-rhobf(k-1)*a_in(i,j,k-1)) &
                         +   (rhobf(k+3)*a_in(i,j,k+3)-rhobf(k-2)*a_in(i,j,k-2)) ) &
                    - ( w0(i,j,k) + w0(i,j-1,k) ) / 60&
                      * ( 37*(rhobf(k  )*a_in(i,j,k  )+rhobf(k-1)*a_in(i,j,k-1)) &
                         - 8*(rhobf(k+1)*a_in(i,j,k+1)+rhobf(k-2)*a_in(i,j,k-2)) &
                         +   (rhobf(k+2)*a_in(i,j,k+2)+rhobf(k-3)*a_in(i,j,k-3)) ) &
                    + sign(1._field_r,(w0(i,j,k)+w0(i,j-1,k)))*(w0(i,j,k)+w0(i,j-1,k)) / 60 &
                      * ( 10*(rhobf(k  )*a_in(i,j,k  )-rhobf(k-1)*a_in(i,j,k-1)) &
                         - 5*(rhobf(k+1)*a_in(i,j,k+1)-rhobf(k-2)*a_in(i,j,k-2)) &
                         +   (rhobf(k+2)*a_in(i,j,k+2)-rhobf(k-3)*a_in(i,j,k-3)) ) &
                ) * (0.5_field_r * dzfi(k)) &
              )
      end do
    end do
  end do
end subroutine vadvecv_5th

subroutine hadvecw_5th(a_in, a_out)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi5,dyi5
  use modfields, only : u0, v0
  implicit none

  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: a_in !< Input: the w field
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: a_out !< Output: the tendency

  integer :: i,j,k

  !$acc parallel loop collapse(3) default(present) async(3)
  do k = 2, kmax
    do j = 2, j1
      do i = 2, i1
        a_out(i,j,k) = a_out(i,j,k) - ( &
              (&
                  (u0(i+1,j,k)+u0(i+1,j,k-1))/60&
                  *(37*(a_in(i+1,j,k)+a_in(i,j,k))-8*(a_in(i+2,j,k)+a_in(i-1,j,k))+(a_in(i+3,j,k)+a_in(i-2,j,k)))&
                  -sign(1._field_r,(u0(i+1,j,k)+u0(i+1,j,k-1)))*(u0(i+1,j,k)+u0(i+1,j,k-1))/60&
                  *(10*(a_in(i+1,j,k)-a_in(i,j,k))-5*(a_in(i+2,j,k)-a_in(i-1,j,k))+(a_in(i+3,j,k)-a_in(i-2,j,k)))&
                  -(u0(i,j,k)+u0(i,j,k-1))/60&
                  *(37*(a_in(i,j,k)+a_in(i-1,j,k))-8*(a_in(i+1,j,k)+a_in(i-2,j,k))+(a_in(i+2,j,k)+a_in(i-3,j,k)))&
                  +sign(1._field_r,(u0(i,j,k)+u0(i,j,k-1)))*(u0(i,j,k)+u0(i,j,k-1))/60&
                  *(10*(a_in(i,j,k)-a_in(i-1,j,k))-5*(a_in(i+1,j,k)-a_in(i-2,j,k))+(a_in(i+2,j,k)-a_in(i-3,j,k)))&
              )*dxi5&
            + (&
                  (v0(i,j+1,k)+v0(i,j+1,k-1))/60&
                  *(37*(a_in(i,j+1,k)+a_in(i,j,k))-8*(a_in(i,j+2,k)+a_in(i,j-1,k))+(a_in(i,j+3,k)+a_in(i,j-2,k)))&
                  -sign(1._field_r,(v0(i,j+1,k)+v0(i,j+1,k-1)))*(v0(i,j+1,k)+v0(i,j+1,k-1))/60&
                  *(10*(a_in(i,j+1,k)-a_in(i,j,k))-5*(a_in(i,j+2,k)-a_in(i,j-1,k))+(a_in(i,j+3,k)-a_in(i,j-2,k)))&
                  -(v0(i,j,k)+v0(i,j,k-1))/60&
                  *(37*(a_in(i,j,k)+a_in(i,j-1,k))-8*(a_in(i,j+1,k)+a_in(i,j-2,k))+(a_in(i,j+2,k)+a_in(i,j-3,k)))&
                  +sign(1._field_r,(v0(i,j,k)+v0(i,j,k-1)))*(v0(i,j,k)+v0(i,j,k-1))/60&
                  *(10*(a_in(i,j,k)-a_in(i,j-1,k))-5*(a_in(i,j+1,k)-a_in(i,j-2,k))+(a_in(i,j+2,k)-a_in(i,j-3,k)))&
              )* dyi5 &
              )
      end do
    end do
  end do
end subroutine hadvecw_5th

!> Vertical advection at the w point.
subroutine vadvecw_5th(a_in, a_out)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dzhi
  use modfields, only : w0, rhobh
  implicit none

  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: a_in !< Input: the w field
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: a_out !< Output: the tendency

  integer :: i,j,k

  k = 2
  !$acc parallel loop collapse(2) default(present) async(5)
  do j = 2, j1
    do i = 2, i1
      a_out(i,j,k) = a_out(i,j,k) - (  (1/rhobh(k)) * ( &
                     ( rhobh(k)*a_in(i,j,k)+rhobh(k+1)*a_in(i,j,k+1) )*(w0(i,j,k) + w0(i,j,k+1)) &
                   - ( rhobh(k)*a_in(i,j,k)+rhobh(k-1)*a_in(i,j,k-1) )*(w0(i,j,k) + w0(i,j,k-1)) &
                  ) * (0.25_field_r * dzhi(k)) &
            )
    end do
  end do

  k = kmax-1
  !$acc parallel loop collapse(2) default(present) async(5)
  do j = 2, j1
    do i = 2, i1
      a_out(i,j,k) = a_out(i,j,k) - (  (1/rhobh(k)) * ( &
                     ( rhobh(k)*a_in(i,j,k)+rhobh(k+1)*a_in(i,j,k+1) )*(w0(i,j,k) + w0(i,j,k+1)) &
                   - ( rhobh(k)*a_in(i,j,k)+rhobh(k-1)*a_in(i,j,k-1) )*(w0(i,j,k) + w0(i,j,k-1)) &
                  ) * (0.25_field_r * dzhi(k)) &
            )
    end do
  end do

  k = kmax
  !$acc parallel loop collapse(2) default(present) async(5)
  do j = 2, j1
    do i = 2, i1
      a_out(i,j,k) = a_out(i,j,k) - (  (1/rhobh(k)) * ( &
                     ( rhobh(k)*a_in(i,j,k)+rhobh(k+1)*a_in(i,j,k+1) )*(w0(i,j,k) + w0(i,j,k+1)) &
                   - ( rhobh(k)*a_in(i,j,k)+rhobh(k-1)*a_in(i,j,k-1) )*(w0(i,j,k) + w0(i,j,k-1)) &
                  ) * (0.25_field_r * dzhi(k)) &
            )
    end do
  end do

  k = 3
  !$acc parallel loop collapse(2) default(present) async(5)
  do j = 2, j1
    do i = 2, i1
      a_out(i,j,k) = a_out(i,j,k) - ( (1/rhobh(k)) * ( &
                    ( w0(i,j,k) + w0(i,j,k+1) ) / 60 &
                    * ( 37*(rhobh(k+1)*a_in(i,j,k+1)+rhobh(k  )*a_in(i,j,k  )) &
                       - 8*(rhobh(k+2)*a_in(i,j,k+2)+rhobh(k-1)*a_in(i,j,k-1)) &
                       +   (rhobh(k+3)*a_in(i,j,k+3)+rhobh(k-2)*a_in(i,j,k-2)) ) &
                  - sign(1._field_r,(w0(i,j,k)+w0(i,j,k+1)))*(w0(i,j,k)+w0(i,j,k+1)) / 60 &
                    * ( 10*(rhobh(k+1)*a_in(i,j,k+1)-rhobh(k  )*a_in(i,j,k  )) &
                       - 5*(rhobh(k+2)*a_in(i,j,k+2)-rhobh(k-1)*a_in(i,j,k-1)) &
                       +   (rhobh(k+3)*a_in(i,j,k+3)-rhobh(k-2)*a_in(i,j,k-2)) ) &
                  - (w0(i,j,k) + w0(i,j,k-1)) / 2 &
                    *      (rhobh(k  )*a_in(i,j,k  )+rhobh(k-1)*a_in(i,j,k-1)) &
              ) * (0.5_field_r * dzhi(k)) &
            )
    end do
  end do

  k = kmax-2
  !$acc parallel loop collapse(2) default(present) async(5)
  do j = 2, j1
    do i = 2, i1
      a_out(i,j,k) = a_out(i,j,k) - ( (1/rhobh(k)) * ( &
                    ( w0(i,j,k) + w0(i,j,k+1) ) / 2 &
                     *      (rhobh(k  )*a_in(i,j,k  )+rhobh(k+1)*a_in(i,j,k+1)) &
                  - ( w0(i,j,k) + w0(i,j,k-1) ) / 60 &
                     * ( 37*(rhobh(k  )*a_in(i,j,k  )+rhobh(k-1)*a_in(i,j,k-1)) &
                        - 8*(rhobh(k+1)*a_in(i,j,k+1)+rhobh(k-2)*a_in(i,j,k-2)) &
                        +   (rhobh(k+2)*a_in(i,j,k+2)+rhobh(k-3)*a_in(i,j,k-3)) ) &
                  + sign(1._field_r,(w0(i,j,k)+w0(i,j,k-1)))*(w0(i,j,k)+w0(i,j,k-1)) / 60 &
                     * ( 10*(rhobh(k  )*a_in(i,j,k  )-rhobh(k-1)*a_in(i,j,k-1)) &
                        - 5*(rhobh(k+1)*a_in(i,j,k+1)-rhobh(k-2)*a_in(i,j,k-2)) &
                        +   (rhobh(k+2)*a_in(i,j,k+2)-rhobh(k-3)*a_in(i,j,k-3)) ) &
              ) * (0.5_field_r * dzhi(k)) &
            )
    end do
  end do

  !$acc parallel loop collapse(3) default(present) async(6)
  do k = 4, kmax-3
    do j = 2, j1
      do i = 2, i1
         a_out(i,j,k) = a_out(i,j,k) - ( (1/rhobh(k)) * ( &
                       ( w0(i,j,k) + w0(i,j,k+1) ) / 60 &
                        * ( 37*(rhobh(k+1)*a_in(i,j,k+1)+rhobh(k  )*a_in(i,j,k  )) &
                           - 8*(rhobh(k+2)*a_in(i,j,k+2)+rhobh(k-1)*a_in(i,j,k-1)) &
                           +   (rhobh(k+3)*a_in(i,j,k+3)+rhobh(k-2)*a_in(i,j,k-2)) ) &
                     - sign(1._field_r,(w0(i,j,k)+w0(i,j,k+1)))*(w0(i,j,k)+w0(i,j,k+1)) / 60 &
                        * ( 10*(rhobh(k+1)*a_in(i,j,k+1)-rhobh(k  )*a_in(i,j,k  )) &
                           - 5*(rhobh(k+2)*a_in(i,j,k+2)-rhobh(k-1)*a_in(i,j,k-1)) &
                           +   (rhobh(k+3)*a_in(i,j,k+3)-rhobh(k-2)*a_in(i,j,k-2)) ) &
                     - ( w0(i,j,k) + w0(i,j,k-1) ) / 60 &
                        * ( 37*(rhobh(k  )*a_in(i,j,k  )+rhobh(k-1)*a_in(i,j,k-1)) &
                           - 8*(rhobh(k+1)*a_in(i,j,k+1)+rhobh(k-2)*a_in(i,j,k-2)) &
                           +   (rhobh(k+2)*a_in(i,j,k+2)+rhobh(k-3)*a_in(i,j,k-3)) ) &
                     + sign(1._field_r,(w0(i,j,k)+w0(i,j,k-1)))*(w0(i,j,k)+w0(i,j,k-1)) / 60 &
                        * ( 10*(rhobh(k  )*a_in(i,j,k  )-rhobh(k-1)*a_in(i,j,k-1)) &
                           - 5*(rhobh(k+1)*a_in(i,j,k+1)-rhobh(k-2)*a_in(i,j,k-2)) &
                           +   (rhobh(k+2)*a_in(i,j,k+2)-rhobh(k-3)*a_in(i,j,k-3)) ) &
                 ) * (0.5_field_r * dzhi(k)) &
               )
      end do
    end do
  end do
end subroutine vadvecw_5th
end module advec_5th
