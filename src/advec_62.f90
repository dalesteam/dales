!> \file advec_62.f90
!!  Does advection with a 6th order horizontal central differencing scheme.
!!  Vertical advection is done with 2nd order scheme.
!! \par Revision list
!! \par Authors
!! \see Wicker and Skamarock 2002
!!
!! A higher-order accuracy in the calculation of the advection is reached with a
!! sixth order central differencing scheme.
!! \latexonly
!!!! \endlatexonly
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

module advec_62
use modprecision, only : field_r
contains
!> Advection at cell center
subroutine advecc_62(a_in, a_out)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi,dyi,dzf,dzh,leq
  use modfields, only : u0, v0, w0,rhobf

  implicit none

  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: a_in !< Input: the cell centered field
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: a_out !< Output: the tendency
  real(field_r)                                       :: inv2dzfk, rhobf_p, rhobf_m

  integer :: i,j,k

  if (leq) then

    k = 1
    inv2dzfk = 1./(2. * dzf(k))
    rhobf_p = rhobf(k+1)/rhobf(k)
    !$acc parallel loop collapse(2) default(present) async(1)
    do j = 2, j1
      do i = 2, i1
        a_out(i,j,k)  = a_out(i,j,k)- ( &
             ( &
             u0(i+1,j,k)/60 &
             *(37*(a_in(i+1,j,k)+a_in(i,j,k))-8*(a_in(i+2,j,k)+a_in(i-1,j,k))+(a_in(i+3,j,k)+a_in(i-2,j,k)))&
             -u0(i,j,k)/60 &
             *(37*(a_in(i,j,k)+a_in(i-1,j,k))-8*(a_in(i+1,j,k)+a_in(i-2,j,k))+(a_in(i+2,j,k)+a_in(i-3,j,k)))&
             )*dxi&
             +(&
             v0(i,j+1,k)/60 &
             *(37*(a_in(i,j+1,k)+a_in(i,j,k))-8*(a_in(i,j+2,k)+a_in(i,j-1,k))+(a_in(i,j+3,k)+a_in(i,j-2,k)))&
             -v0(i,j,k)/60 &
             *(37*(a_in(i,j,k)+a_in(i,j-1,k))-8*(a_in(i,j+1,k)+a_in(i,j-2,k))+(a_in(i,j+2,k)+a_in(i,j-3,k)))&
             )* dyi &
             + ( &
             w0(i,j,k+1) * (rhobf_p * a_in(i,j,k+1) + a_in(i,j,k)) &
             ) * inv2dzfk &
             )
      end do
    end do

    !$acc parallel loop collapse(3) default(present) async(2)
    do k = 2, kmax
      do j = 2, j1
        do i = 2, i1
          inv2dzfk = 1./(2. * dzf(k))
          rhobf_p = rhobf(k+1)/rhobf(k)
          rhobf_m = rhobf(k-1)/rhobf(k)
          a_out(i,j,k)  = a_out(i,j,k)- (  &
              ( &
                  u0(i+1,j,k)/60. &
                  *(37.*(a_in(i+1,j,k)+a_in(i,j,k))-8.*(a_in(i+2,j,k)+a_in(i-1,j,k))+(a_in(i+3,j,k)+a_in(i-2,j,k)))&
                  -u0(i,j,k)/60. &
                  *(37.*(a_in(i,j,k)+a_in(i-1,j,k))-8.*(a_in(i+1,j,k)+a_in(i-2,j,k))+(a_in(i+2,j,k)+a_in(i-3,j,k)))&
              )*dxi&
            +(&
                  v0(i,j+1,k)/60. &
                  *(37.*(a_in(i,j+1,k)+a_in(i,j,k))-8.*(a_in(i,j+2,k)+a_in(i,j-1,k))+(a_in(i,j+3,k)+a_in(i,j-2,k)))&
                  -v0(i,j,k)/60. &
                  *(37.*(a_in(i,j,k)+a_in(i,j-1,k))-8.*(a_in(i,j+1,k)+a_in(i,j-2,k))+(a_in(i,j+2,k)+a_in(i,j-3,k)))&
              )* dyi &
            + ( &
              w0(i,j,k+1) * (rhobf_p * a_in(i,j,k+1) + a_in(i,j,k)) &
              -w0(i,j,k)  * (rhobf_m * a_in(i,j,k-1) + a_in(i,j,k)) &
              ) * inv2dzfk &
              )
        end do
      end do
    end do

  else ! non-equidistant grid
    k = 1
    inv2dzfk = 1./(2. * dzf(k))
    rhobf_p = rhobf(k+1)/rhobf(k)
    !$acc parallel loop collapse(2) default(present) async(1)
    do j = 2, j1
      do i = 2, i1
        a_out(i,j,k)  = a_out(i,j,k) - ( &
             ( &
             u0(i+1,j,k)/60 &
             *(37*(a_in(i+1,j,k)+a_in(i,j,k))-8*(a_in(i+2,j,k)+a_in(i-1,j,k))+(a_in(i+3,j,k)+a_in(i-2,j,k)))&
             -u0(i,j,k)/60 &
             *(37*(a_in(i,j,k)+a_in(i-1,j,k))-8*(a_in(i+1,j,k)+a_in(i-2,j,k))+(a_in(i+2,j,k)+a_in(i-3,j,k)))&
             )*dxi&
             +(&
             v0(i,j+1,k)/60 &
             *(37*(a_in(i,j+1,k)+a_in(i,j,k))-8*(a_in(i,j+2,k)+a_in(i,j-1,k))+(a_in(i,j+3,k)+a_in(i,j-2,k)))&
             -v0(i,j,k)/60 &
             *(37*(a_in(i,j,k)+a_in(i,j-1,k))-8*(a_in(i,j+1,k)+a_in(i,j-2,k))+(a_in(i,j+2,k)+a_in(i,j-3,k)))&
             )* dyi &
             + ( &
             w0(i,j,k+1) * (rhobf_p * a_in(i,j,k+1) * dzf(k) +  a_in(i,j,k) * dzf(k+1) ) / dzh(k+1) &
             ) * inv2dzfk &
             )
      end do
    end do

    !$acc parallel loop collapse(3) default(present) async(2)
    do k = 2, kmax
      do j = 2, j1
        do i = 2, i1
          inv2dzfk = 1./(2. * dzf(k))
          rhobf_p = rhobf(k+1)/rhobf(k)
          rhobf_m = rhobf(k-1)/rhobf(k)
          a_out(i,j,k)  = a_out(i,j,k)- (  &
              ( &
                  u0(i+1,j,k)/60. &
                  *(37.*(a_in(i+1,j,k)+a_in(i,j,k))-8.*(a_in(i+2,j,k)+a_in(i-1,j,k))+(a_in(i+3,j,k)+a_in(i-2,j,k)))&
                  -u0(i,j,k)/60. &
                  *(37.*(a_in(i,j,k)+a_in(i-1,j,k))-8.*(a_in(i+1,j,k)+a_in(i-2,j,k))+(a_in(i+2,j,k)+a_in(i-3,j,k)))&
              )*dxi&
            +(&
                  v0(i,j+1,k)/60. &
                  *(37.*(a_in(i,j+1,k)+a_in(i,j,k))-8.*(a_in(i,j+2,k)+a_in(i,j-1,k))+(a_in(i,j+3,k)+a_in(i,j-2,k)))&
                  -v0(i,j,k)/60. &
                  *(37.*(a_in(i,j,k)+a_in(i,j-1,k))-8.*(a_in(i,j+1,k)+a_in(i,j-2,k))+(a_in(i,j+2,k)+a_in(i,j-3,k)))&
              )* dyi &
              + ( &
               w0(i,j,k+1) * (rhobf_p * a_in(i,j,k+1) * dzf(k) +  a_in(i,j,k) * dzf(k+1) ) / dzh(k+1) &
              -w0(i,j,k )  * (rhobf_m * a_in(i,j,k-1) * dzf(k) +  a_in(i,j,k) * dzf(k-1) ) / dzh(k) &
              ) * inv2dzfk &
              )
        end do
      end do
    end do
  end if
  !$acc wait(1,2)
end subroutine advecc_62

!> Advection at the u point.
subroutine advecu_62(a_in,a_out)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi5,dyi5,dzf,dzh,leq
  use modfields, only : u0, v0, w0,rhobf

  implicit none

  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: a_in !< Input: the u field
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: a_out !< Output: the tendency
  real(field_r)                                       :: inv4dzfk, rhobf_p, rhobf_m

  integer :: i,j,k

  if (leq) then
    k = 1
    inv4dzfk = 1./(4. * dzf(k))
    rhobf_p = rhobf(k+1)/rhobf(k)
    !$acc parallel loop collapse(2) default(present) async(1)
    do j = 2, j1
      do i = 2, i1
        a_out(i,j,k) = a_out(i,j,k)- ( &
             ( &
             (u0(i+1,j,k)+u0(i,j,k))/60 &
             *(37*(a_in(i+1,j,k)+a_in(i,j,k))-8*(a_in(i+2,j,k)+a_in(i-1,j,k))+(a_in(i+3,j,k)+a_in(i-2,j,k)))&
             -(u0(i,j,k)+u0(i-1,j,k))/60 &
             *(37*(a_in(i,j,k)+a_in(i-1,j,k))-8*(a_in(i+1,j,k)+a_in(i-2,j,k))+(a_in(i+2,j,k)+a_in(i-3,j,k)))&
             )*dxi5&
             +(&
             (v0(i,j+1,k)+v0(i-1,j+1,k))/60 &
             *(37*(a_in(i,j+1,k)+a_in(i,j,k))-8*(a_in(i,j+2,k)+a_in(i,j-1,k))+(a_in(i,j+3,k)+a_in(i,j-2,k)))&
             -(v0(i,j,k)+v0(i-1,j,k))/60 &
             *(37*(a_in(i,j,k)+a_in(i,j-1,k))-8*(a_in(i,j+1,k)+a_in(i,j-2,k))+(a_in(i,j+2,k)+a_in(i,j-3,k)))&
             )* dyi5 &
             + ( &
             (rhobf_p * a_in(i,j,k+1) + a_in(i,j,k)) *(w0(i,j,k+1)+ w0(i-1,j,k+1)) &
             ) * inv4dzfk &
             )
      end do
    end do

    !$acc parallel loop collapse(3) default(present) async(1)
    do k = 2, kmax
      do j = 2, j1
        do i = 2, i1
          inv4dzfk = 1./(4. * dzf(k))
          rhobf_p = rhobf(k+1)/rhobf(k)
          rhobf_m = rhobf(k-1)/rhobf(k)
          a_out(i,j,k)  = a_out(i,j,k)- ( &
                (&
                    (u0(i+1,j,k)+u0(i,j,k))/60 &
                    *(37*(a_in(i+1,j,k)+a_in(i,j,k))-8*(a_in(i+2,j,k)+a_in(i-1,j,k))+(a_in(i+3,j,k)+a_in(i-2,j,k)))&
                    -(u0(i,j,k)+u0(i-1,j,k))/60 &
                    *(37*(a_in(i,j,k)+a_in(i-1,j,k))-8*(a_in(i+1,j,k)+a_in(i-2,j,k))+(a_in(i+2,j,k)+a_in(i-3,j,k)))&
                )*dxi5&
              +(&
                    (v0(i,j+1,k)+v0(i-1,j+1,k))/60 &
                    *(37*(a_in(i,j+1,k)+a_in(i,j,k))-8*(a_in(i,j+2,k)+a_in(i,j-1,k))+(a_in(i,j+3,k)+a_in(i,j-2,k)))&
                    -(v0(i,j,k)+v0(i-1,j,k))/60 &
                    *(37*(a_in(i,j,k)+a_in(i,j-1,k))-8*(a_in(i,j+1,k)+a_in(i,j-2,k))+(a_in(i,j+2,k)+a_in(i,j-3,k)))&
                )* dyi5 &
              + ( &
                     (a_in(i,j,k) + rhobf_p * a_in(i,j,k+1) )*(w0(i,j,k+1)+w0(i-1,j,k+1)) &
                    -(a_in(i,j,k) + rhobf_m * a_in(i,j,k-1) )*(w0(i,j,k  )+w0(i-1,j,k  )) &
                ) * inv4dzfk &
                )
        end do
      end do
    end do

  else ! non-equidistant grid
    k = 1
    inv4dzfk = 1./(4. * dzf(k))
    rhobf_p = rhobf(k+1)/rhobf(k)
    !$acc parallel loop collapse(2) default(present) async(1)
    do j = 2, j1
      do i = 2, i1
        a_out(i,j,k)  = a_out(i,j,k)- ( &
             ( &
             (u0(i+1,j,k)+u0(i,j,k))/60 &
             *(37*(a_in(i+1,j,k)+a_in(i,j,k))-8*(a_in(i+2,j,k)+a_in(i-1,j,k))+(a_in(i+3,j,k)+a_in(i-2,j,k)))&
             -(u0(i,j,k)+u0(i-1,j,k))/60 &
             *(37*(a_in(i,j,k)+a_in(i-1,j,k))-8*(a_in(i+1,j,k)+a_in(i-2,j,k))+(a_in(i+2,j,k)+a_in(i-3,j,k)))&
             )*dxi5&
             +(&
             (v0(i,j+1,k)+v0(i-1,j+1,k))/60 &
             *(37*(a_in(i,j+1,k)+a_in(i,j,k))-8*(a_in(i,j+2,k)+a_in(i,j-1,k))+(a_in(i,j+3,k)+a_in(i,j-2,k)))&
             -(v0(i,j,k)+v0(i-1,j,k))/60 &
             *(37*(a_in(i,j,k)+a_in(i,j-1,k))-8*(a_in(i,j+1,k)+a_in(i,j-2,k))+(a_in(i,j+2,k)+a_in(i,j-3,k)))&
             )* dyi5 &
             + ( &
             (a_in(i,j,k) * dzf(k+1) + rhobf_p * a_in(i,j,k+1)*dzf(k)) * (w0(i,j,k+1)+w0(i-1,j,k+1)) / dzh(k+1) &
             ) * inv4dzfk &
             )
      end do
    end do

    !$acc parallel loop collapse(3) default(present) async(1)
    do k = 2, kmax
      do j = 2, j1
        do i = 2, i1
          inv4dzfk = 1./(4. * dzf(k))
          rhobf_p = rhobf(k+1)/rhobf(k)
          rhobf_m = rhobf(k-1)/rhobf(k)
          a_out(i,j,k)  = a_out(i,j,k)- ( &
                (&
                    (u0(i+1,j,k)+u0(i,j,k))/60 &
                    *(37*(a_in(i+1,j,k)+a_in(i,j,k))-8*(a_in(i+2,j,k)+a_in(i-1,j,k))+(a_in(i+3,j,k)+a_in(i-2,j,k)))&
                    -(u0(i,j,k)+u0(i-1,j,k))/60 &
                    *(37*(a_in(i,j,k)+a_in(i-1,j,k))-8*(a_in(i+1,j,k)+a_in(i-2,j,k))+(a_in(i+2,j,k)+a_in(i-3,j,k)))&
                )*dxi5&
              +(&
                    (v0(i,j+1,k)+v0(i-1,j+1,k))/60 &
                    *(37*(a_in(i,j+1,k)+a_in(i,j,k))-8*(a_in(i,j+2,k)+a_in(i,j-1,k))+(a_in(i,j+3,k)+a_in(i,j-2,k)))&
                    -(v0(i,j,k)+v0(i-1,j,k))/60 &
                    *(37*(a_in(i,j,k)+a_in(i,j-1,k))-8*(a_in(i,j+1,k)+a_in(i,j-2,k))+(a_in(i,j+2,k)+a_in(i,j-3,k)))&
                )* dyi5 &
              + ( &
                     (a_in(i,j,k) * dzf(k+1) + rhobf_p * a_in(i,j,k+1) *dzf(k)) * (w0(i,j,k+1)+w0(i-1,j,k+1)) / dzh(k+1) &
                    -(a_in(i,j,k) * dzf(k-1) + rhobf_m * a_in(i,j,k-1) *dzf(k)) * (w0(i,j,k  )+w0(i-1,j,k  )) / dzh(k)   &
                ) * inv4dzfk &
                )
        end do
      end do
    end do
  end if
end subroutine advecu_62

!> Advection at the v point.
subroutine advecv_62(a_in, a_out)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi5,dyi5,dzf,dzh,leq
  use modfields, only : u0, v0, w0,rhobf
  implicit none

  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: a_in !< Input: the v field
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: a_out !< Output: the tendency
  real(field_r)                                       :: inv4dzfk, rhobf_p, rhobf_m

  integer :: i,j,k

  if (leq) then
    k = 1
    inv4dzfk = 1./(4. * dzf(k))
    rhobf_p = rhobf(k+1)/rhobf(k)

    !$acc parallel loop collapse(2) default(present) async(2)
    do j = 2, j1
      do i = 2, i1
        a_out(i,j,k)  = a_out(i,j,k)- ( &
             ( &
             (u0(i+1,j,k)+u0(i+1,j-1,k))/60 &
             *(37*(a_in(i+1,j,k)+a_in(i,j,k))-8*(a_in(i+2,j,k)+a_in(i-1,j,k))+(a_in(i+3,j,k)+a_in(i-2,j,k)))&
             -(u0(i,j,k)+u0(i,j-1,k))/60 &
             *(37*(a_in(i,j,k)+a_in(i-1,j,k))-8*(a_in(i+1,j,k)+a_in(i-2,j,k))+(a_in(i+2,j,k)+a_in(i-3,j,k)))&
             )*dxi5&
             +(&
             (v0(i,j+1,k)+v0(i,j,k))/60 &
             *(37*(a_in(i,j+1,k)+a_in(i,j,k))-8*(a_in(i,j+2,k)+a_in(i,j-1,k))+(a_in(i,j+3,k)+a_in(i,j-2,k)))&
             -(v0(i,j,k)+v0(i,j-1,k))/60 &
             *(37*(a_in(i,j,k)+a_in(i,j-1,k))-8*(a_in(i,j+1,k)+a_in(i,j-2,k))+(a_in(i,j+2,k)+a_in(i,j-3,k)))&
             )* dyi5 &
             +( &
             (w0(i,j,k+1)+w0(i,j-1,k+1)) *(rhobf_p * a_in(i,j,k+1)+a_in(i,j,k)) &
             ) * inv4dzfk  &
             )
      end do
    end do

    !$acc parallel loop collapse(3) default(present) async(2)
    do k = 2, kmax
      do j = 2, j1
        do i = 2, i1
          inv4dzfk = 1./(4. * dzf(k))
          rhobf_p = rhobf(k+1)/rhobf(k)
          rhobf_m = rhobf(k-1)/rhobf(k)
          a_out(i,j,k)  = a_out(i,j,k)- ( &
                ( &
                    (u0(i+1,j,k)+u0(i+1,j-1,k))/60 &
                    *(37*(a_in(i+1,j,k)+a_in(i,j,k))-8*(a_in(i+2,j,k)+a_in(i-1,j,k))+(a_in(i+3,j,k)+a_in(i-2,j,k)))&
                    -(u0(i,j,k)+u0(i,j-1,k))/60 &
                    *(37*(a_in(i,j,k)+a_in(i-1,j,k))-8*(a_in(i+1,j,k)+a_in(i-2,j,k))+(a_in(i+2,j,k)+a_in(i-3,j,k)))&
                 )*dxi5&
                +(&
                    (v0(i,j+1,k)+v0(i,j,k))/60 &
                    *(37*(a_in(i,j+1,k)+a_in(i,j,k))-8*(a_in(i,j+2,k)+a_in(i,j-1,k))+(a_in(i,j+3,k)+a_in(i,j-2,k)))&
                    -(v0(i,j,k)+v0(i,j-1,k))/60 &
                    *(37*(a_in(i,j,k)+a_in(i,j-1,k))-8*(a_in(i,j+1,k)+a_in(i,j-2,k))+(a_in(i,j+2,k)+a_in(i,j-3,k)))&
                  )* dyi5 &
                + ( &
                (w0(i,j,k+1)+w0(i,j-1,k+1))*(rhobf_p * a_in(i,j,k+1) + a_in(i,j,k))&
                -(w0(i,j,k) +w0(i,j-1,k))  *(rhobf_m * a_in(i,j,k-1) + a_in(i,j,k))&
                  ) * inv4dzfk  &
                  )
        end do
      end do
    end do

  else ! non-equidistant grid
    k = 1
    inv4dzfk = 1./(4. * dzf(k))
    rhobf_p = rhobf(k+1)/rhobf(k)
    !$acc parallel loop collapse(2) default(present) async(2)
    do j = 2, j1
      do i = 2, i1
         a_out(i,j,k)  = a_out(i,j,k)- ( &
              ( &
              (u0(i+1,j,k)+u0(i+1,j-1,k))/60. &
              *(37.*(a_in(i+1,j,k)+a_in(i,j,k))-8.*(a_in(i+2,j,k)+a_in(i-1,j,k))+(a_in(i+3,j,k)+a_in(i-2,j,k)))&
              -(u0(i,j,k)+u0(i,j-1,k))/60. &
              *(37.*(a_in(i,j,k)+a_in(i-1,j,k))-8.*(a_in(i+1,j,k)+a_in(i-2,j,k))+(a_in(i+2,j,k)+a_in(i-3,j,k)))&
              )*dxi5&
              +(&
              (v0(i,j+1,k)+v0(i,j,k))/60. &
              *(37.*(a_in(i,j+1,k)+a_in(i,j,k))-8.*(a_in(i,j+2,k)+a_in(i,j-1,k))+(a_in(i,j+3,k)+a_in(i,j-2,k)))&
              -(v0(i,j,k)+v0(i,j-1,k))/60. &
              *(37.*(a_in(i,j,k)+a_in(i,j-1,k))-8.*(a_in(i,j+1,k)+a_in(i,j-2,k))+(a_in(i,j+2,k)+a_in(i,j-3,k)))&
              )* dyi5 &
              +( &
              (w0(i,j,k+1)+w0(i,j-1,k+1)) * (rhobf_p * a_in(i,j,k+1)*dzf(k) + a_in(i,j,k)*dzf(k+1)) / dzh(k+1)&
              ) * inv4dzfk  &
              )
      end do
    end do

    !$acc parallel loop collapse(3) default(present) async(2)
    do k = 2, kmax
      do j = 2, j1
        do i = 2, i1
          inv4dzfk = 1./(4. * dzf(k))
          rhobf_p = rhobf(k+1)/rhobf(k)
          rhobf_m = rhobf(k-1)/rhobf(k)
          a_out(i,j,k)  = a_out(i,j,k)- ( &
                ( &
                    (u0(i+1,j,k)+u0(i+1,j-1,k))/60 &
                    *(37*(a_in(i+1,j,k)+a_in(i,j,k))-8*(a_in(i+2,j,k)+a_in(i-1,j,k))+(a_in(i+3,j,k)+a_in(i-2,j,k)))&
                    -(u0(i,j,k)+u0(i,j-1,k))/60 &
                    *(37*(a_in(i,j,k)+a_in(i-1,j,k))-8*(a_in(i+1,j,k)+a_in(i-2,j,k))+(a_in(i+2,j,k)+a_in(i-3,j,k)))&
                 )*dxi5&
                +(&
                    (v0(i,j+1,k)+v0(i,j,k))/60 &
                    *(37*(a_in(i,j+1,k)+a_in(i,j,k))-8*(a_in(i,j+2,k)+a_in(i,j-1,k))+(a_in(i,j+3,k)+a_in(i,j-2,k)))&
                    -(v0(i,j,k)+v0(i,j-1,k))/60 &
                    *(37*(a_in(i,j,k)+a_in(i,j-1,k))-8*(a_in(i,j+1,k)+a_in(i,j-2,k))+(a_in(i,j+2,k)+a_in(i,j-3,k)))&
                  )* dyi5 &
                  + ( &
                  (w0(i,j,k+1)+w0(i,j-1,k+1)) * (rhobf_p * a_in(i,j,k+1) * dzf(k) + a_in(i,j,k) * dzf(k+1)) / dzh(k+1)&
                 -(w0(i,j,k)  +w0(i,j-1,k))   * (rhobf_m * a_in(i,j,k-1) * dzf(k) + a_in(i,j,k) * dzf(k-1)) / dzh(k)  &
                  ) * inv4dzfk  &
                  )
        end do
      end do
    end do
  end if
end subroutine advecv_62

!> Advection at the w point.
subroutine advecw_62(a_in, a_out)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi5,dyi5,dzh
  use modfields, only : u0, v0, w0,rhobh
  implicit none

  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: a_in !< Input: the w field
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: a_out !< Output: the tendency
  real(field_r)                                       :: inv4dzhk, rhobh_p, rhobh_m

  integer :: i,j,k

  !$acc parallel loop collapse(3) default(present) async(3)
  do k = 2, kmax
    do j = 2, j1
      do i = 2, i1
        inv4dzhk = 1./(4. * dzh(k))
        rhobh_p = rhobh(k+1)/rhobh(k)
        rhobh_m = rhobh(k-1)/rhobh(k)
        a_out(i,j,k)  = a_out(i,j,k)- ( &
              (&
                  (u0(i+1,j,k)+u0(i+1,j,k-1))/60. &
                  *(37.*(a_in(i+1,j,k)+a_in(i,j,k))-8.*(a_in(i+2,j,k)+a_in(i-1,j,k))+(a_in(i+3,j,k)+a_in(i-2,j,k)))&
                  -(u0(i,j,k)+u0(i,j,k-1))/60. &
                  *(37.*(a_in(i,j,k)+a_in(i-1,j,k))-8.*(a_in(i+1,j,k)+a_in(i-2,j,k))+(a_in(i+2,j,k)+a_in(i-3,j,k)))&
              )*dxi5&
             +(&
                  (v0(i,j+1,k)+v0(i,j+1,k-1))/60. &
                  *(37.*(a_in(i,j+1,k)+a_in(i,j,k))-8.*(a_in(i,j+2,k)+a_in(i,j-1,k))+(a_in(i,j+3,k)+a_in(i,j-2,k)))&
                  -(v0(i,j,k)+v0(i,j,k-1))/60. &
                  *(37.*(a_in(i,j,k)+a_in(i,j-1,k))-8.*(a_in(i,j+1,k)+a_in(i,j-2,k))+(a_in(i,j+2,k)+a_in(i,j-3,k)))&
               )* dyi5 &
             + ( &
                   (a_in(i,j,k)+rhobh_p * a_in(i,j,k+1) )*(w0(i,j,k) + w0(i,j,k+1)) &
                  -(a_in(i,j,k)+rhobh_m * a_in(i,j,k-1) )*(w0(i,j,k) + w0(i,j,k-1)) &
               )*inv4dzhk &
               )
      end do
    end do
  end do
  !Advection for u, v and w called sequentially in modavection. Only sync here.
  !$acc wait(1,2,3)
end subroutine advecw_62

end module advec_62
