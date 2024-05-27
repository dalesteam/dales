!> \file advec_52.f90
!!  Does advection with a 5th order upwind scheme and 2nd order in the vertical
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

module advec_52
use modprecision, only : field_r
contains
!> Advection at cell center
subroutine advecc_52(a_in, a_out)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi,dyi,dzf,dzh,leq
  use modfields, only : u0, v0, w0,rhobf
  implicit none

  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: a_in !< Input: the cell centered field
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: a_out !< Output: the tendency
  real                                       :: inv2dzfk, rhobf_p, rhobf_m

  integer :: i,j,k

  if (leq) then

    k = 1
    inv2dzfk = 1./(2. * dzf(k))
    rhobf_p = rhobf(k+1)/rhobf(k)

    !$acc parallel loop collapse(2) default(present) async(1)
    do j = 2, j1
      do i = 2, i1
        a_out(i,j,k) = a_out(i,j,k) - ( &
             ( &
             u0(i+1,j,k)/60&
             *(37*(a_in(i+1,j,k)+a_in(i,j,k))-8*(a_in(i+2,j,k)+a_in(i-1,j,k))+(a_in(i+3,j,k)+a_in(i-2,j,k)))&
             -abs(u0(i+1,j,k))/60&
             *(10*(a_in(i+1,j,k)-a_in(i,j,k))-5*(a_in(i+2,j,k)-a_in(i-1,j,k))+(a_in(i+3,j,k)-a_in(i-2,j,k)))&
             -u0(i,j,k)/60&
             *(37*(a_in(i,j,k)+a_in(i-1,j,k))-8*(a_in(i+1,j,k)+a_in(i-2,j,k))+(a_in(i+2,j,k)+a_in(i-3,j,k)))&
             +abs(u0(i,j,k))/60&
             *(10*(a_in(i,j,k)-a_in(i-1,j,k))-5*(a_in(i+1,j,k)-a_in(i-2,j,k))+(a_in(i+2,j,k)-a_in(i-3,j,k)))&
             )*dxi&
             +(&
             v0(i,j+1,k)/60&
             *(37*(a_in(i,j+1,k)+a_in(i,j,k))-8*(a_in(i,j+2,k)+a_in(i,j-1,k))+(a_in(i,j+3,k)+a_in(i,j-2,k)))&
             -abs(v0(i,j+1,k))/60&
             *(10*(a_in(i,j+1,k)-a_in(i,j,k))-5*(a_in(i,j+2,k)-a_in(i,j-1,k))+(a_in(i,j+3,k)-a_in(i,j-2,k)))&
             -v0(i,j,k)/60&
             *(37*(a_in(i,j,k)+a_in(i,j-1,k))-8*(a_in(i,j+1,k)+a_in(i,j-2,k))+(a_in(i,j+2,k)+a_in(i,j-3,k)))&
             +abs(v0(i,j,k))/60&
             *(10*(a_in(i,j,k)-a_in(i,j-1,k))-5*(a_in(i,j+1,k)-a_in(i,j-2,k))+(a_in(i,j+2,k)-a_in(i,j-3,k))) &
             )* dyi &
             + ( &
             w0(i,j,k+1) * (rhobf_p * a_in(i,j,k+1) + a_in(i,j,k)) &
             ) * inv2dzfk  &
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
                    u0(i+1,j,k)/60&
                    *(37*(a_in(i+1,j,k)+a_in(i,j,k))-8*(a_in(i+2,j,k)+a_in(i-1,j,k))+(a_in(i+3,j,k)+a_in(i-2,j,k)))&
                    -abs(u0(i+1,j,k))/60&
                    *(10*(a_in(i+1,j,k)-a_in(i,j,k))-5*(a_in(i+2,j,k)-a_in(i-1,j,k))+(a_in(i+3,j,k)-a_in(i-2,j,k)))&
                    -u0(i,j,k)/60&
                    *(37*(a_in(i,j,k)+a_in(i-1,j,k))-8*(a_in(i+1,j,k)+a_in(i-2,j,k))+(a_in(i+2,j,k)+a_in(i-3,j,k)))&
                    +abs(u0(i,j,k))/60&
                    *(10*(a_in(i,j,k)-a_in(i-1,j,k))-5*(a_in(i+1,j,k)-a_in(i-2,j,k))+(a_in(i+2,j,k)-a_in(i-3,j,k)))&
                )*dxi&
              +(&
                    v0(i,j+1,k)/60&
                    *(37*(a_in(i,j+1,k)+a_in(i,j,k))-8*(a_in(i,j+2,k)+a_in(i,j-1,k))+(a_in(i,j+3,k)+a_in(i,j-2,k)))&
                    -abs(v0(i,j+1,k))/60&
                    *(10*(a_in(i,j+1,k)-a_in(i,j,k))-5*(a_in(i,j+2,k)-a_in(i,j-1,k))+(a_in(i,j+3,k)-a_in(i,j-2,k)))&
                    -v0(i,j,k)/60&
                    *(37*(a_in(i,j,k)+a_in(i,j-1,k))-8*(a_in(i,j+1,k)+a_in(i,j-2,k))+(a_in(i,j+2,k)+a_in(i,j-3,k)))&
                    +abs(v0(i,j,k))/60&
                    *(10*(a_in(i,j,k)-a_in(i,j-1,k))-5*(a_in(i,j+1,k)-a_in(i,j-2,k))+(a_in(i,j+2,k)-a_in(i,j-3,k)))&
                )* dyi &
              + ( &
                w0(i,j,k+1) * (rhobf_p * a_in(i,j,k+1) +  a_in(i,j,k)) &
                -w0(i,j,k)  * (rhobf_m * a_in(i,j,k-1) +  a_in(i,j,k)) &
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
        a_out(i,j,k) = a_out(i,j,k) - ( &
             ( &
             u0(i+1,j,k)/60&
             *(37*(a_in(i+1,j,k)+a_in(i,j,k))-8*(a_in(i+2,j,k)+a_in(i-1,j,k))+(a_in(i+3,j,k)+a_in(i-2,j,k)))&
             -abs(u0(i+1,j,k))/60&
             *(10*(a_in(i+1,j,k)-a_in(i,j,k))-5*(a_in(i+2,j,k)-a_in(i-1,j,k))+(a_in(i+3,j,k)-a_in(i-2,j,k)))&
             -u0(i,j,k)/60&
             *(37*(a_in(i,j,k)+a_in(i-1,j,k))-8*(a_in(i+1,j,k)+a_in(i-2,j,k))+(a_in(i+2,j,k)+a_in(i-3,j,k)))&
             +abs(u0(i,j,k))/60&
             *(10*(a_in(i,j,k)-a_in(i-1,j,k))-5*(a_in(i+1,j,k)-a_in(i-2,j,k))+(a_in(i+2,j,k)-a_in(i-3,j,k)))&
             )*dxi&
             +(&
             v0(i,j+1,k)/60&
             *(37*(a_in(i,j+1,k)+a_in(i,j,k))-8*(a_in(i,j+2,k)+a_in(i,j-1,k))+(a_in(i,j+3,k)+a_in(i,j-2,k)))&
             -abs(v0(i,j+1,k))/60&
             *(10*(a_in(i,j+1,k)-a_in(i,j,k))-5*(a_in(i,j+2,k)-a_in(i,j-1,k))+(a_in(i,j+3,k)-a_in(i,j-2,k)))&
             -v0(i,j,k)/60&
             *(37*(a_in(i,j,k)+a_in(i,j-1,k))-8*(a_in(i,j+1,k)+a_in(i,j-2,k))+(a_in(i,j+2,k)+a_in(i,j-3,k)))&
             +abs(v0(i,j,k))/60&
             *(10*(a_in(i,j,k)-a_in(i,j-1,k))-5*(a_in(i,j+1,k)-a_in(i,j-2,k))+(a_in(i,j+2,k)-a_in(i,j-3,k))) &
             )* dyi &
             + ( &
             w0(i,j,k+1) * (rhobf_p * a_in(i,j,k+1) * dzf(k) +  a_in(i,j,k) * dzf(k+1) ) / dzh(k+1) &
             ) * inv2dzfk  &
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
          a_out(i,j,k) = a_out(i,j,k) - (  &
                ( &
                    u0(i+1,j,k)/60&
                    *(37*(a_in(i+1,j,k)+a_in(i,j,k))-8*(a_in(i+2,j,k)+a_in(i-1,j,k))+(a_in(i+3,j,k)+a_in(i-2,j,k)))&
                    -abs(u0(i+1,j,k))/60&
                    *(10*(a_in(i+1,j,k)-a_in(i,j,k))-5*(a_in(i+2,j,k)-a_in(i-1,j,k))+(a_in(i+3,j,k)-a_in(i-2,j,k)))&
                    -u0(i,j,k)/60&
                    *(37*(a_in(i,j,k)+a_in(i-1,j,k))-8*(a_in(i+1,j,k)+a_in(i-2,j,k))+(a_in(i+2,j,k)+a_in(i-3,j,k)))&
                    +abs(u0(i,j,k))/60&
                    *(10*(a_in(i,j,k)-a_in(i-1,j,k))-5*(a_in(i+1,j,k)-a_in(i-2,j,k))+(a_in(i+2,j,k)-a_in(i-3,j,k)))&
                )*dxi&
              +(&
                    v0(i,j+1,k)/60&
                    *(37*(a_in(i,j+1,k)+a_in(i,j,k))-8*(a_in(i,j+2,k)+a_in(i,j-1,k))+(a_in(i,j+3,k)+a_in(i,j-2,k)))&
                    -abs(v0(i,j+1,k))/60&
                    *(10*(a_in(i,j+1,k)-a_in(i,j,k))-5*(a_in(i,j+2,k)-a_in(i,j-1,k))+(a_in(i,j+3,k)-a_in(i,j-2,k)))&
                    -v0(i,j,k)/60&
                    *(37*(a_in(i,j,k)+a_in(i,j-1,k))-8*(a_in(i,j+1,k)+a_in(i,j-2,k))+(a_in(i,j+2,k)+a_in(i,j-3,k)))&
                    +abs(v0(i,j,k))/60&
                    *(10*(a_in(i,j,k)-a_in(i,j-1,k))-5*(a_in(i,j+1,k)-a_in(i,j-2,k))+(a_in(i,j+2,k)-a_in(i,j-3,k)))&
                )* dyi &
              + ( &
              w0(i,j,k+1) * (rhobf_p * a_in(i,j,k+1) * dzf(k) +  a_in(i,j,k) * dzf(k+1) ) / dzh(k+1) &
              -w0(i,j,k ) * (rhobf_m * a_in(i,j,k-1) * dzf(k) +  a_in(i,j,k) * dzf(k-1) ) / dzh(k) &
                ) * inv2dzfk &
                )
        end do
      end do
    end do
  end if
  !$acc wait(1,2)

end subroutine advecc_52


!> Advection at the u point.
subroutine advecu_52(a_in,a_out)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi5,dyi5,dzf,dzh,leq
  use modfields, only : u0, v0, w0,rhobf
  implicit none

  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: a_in !< Input: the u field
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: a_out !< Output: the tendency

  integer :: i,j,k

  if (leq) then
    k = 1
    !$acc parallel loop collapse(2) default(present) async(1)
    do j = 2, j1
      do i = 2, i1
        a_out(i,j,k)  = a_out(i,j,k)- ( &
             (&
             (u0(i+1,j,k)+u0(i,j,k))/60&
             *(37*(a_in(i+1,j,k)+a_in(i,j,k))-8*(a_in(i+2,j,k)+a_in(i-1,j,k))+(a_in(i+3,j,k)+a_in(i-2,j,k)))&
             -sign(1._field_r,(u0(i+1,j,k)+u0(i,j,k)))*(u0(i+1,j,k)+u0(i,j,k))/60&
             *(10*(a_in(i+1,j,k)-a_in(i,j,k))-5*(a_in(i+2,j,k)-a_in(i-1,j,k))+(a_in(i+3,j,k)-a_in(i-2,j,k)))&
             -(u0(i,j,k)+u0(i-1,j,k))/60&
             *(37*(a_in(i,j,k)+a_in(i-1,j,k))-8*(a_in(i+1,j,k)+a_in(i-2,j,k))+(a_in(i+2,j,k)+a_in(i-3,j,k)))&
             +sign(1._field_r,(u0(i,j,k)+u0(i-1,j,k)))*(u0(i,j,k)+u0(i-1,j,k))/60&
             *(10*(a_in(i,j,k)-a_in(i-1,j,k))-5*(a_in(i+1,j,k)-a_in(i-2,j,k))+(a_in(i+2,j,k)-a_in(i-3,j,k)))&
             )*dxi5 &
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
             +(1./rhobf(k))*( &
             ( rhobf(k+1)*a_in(i,j,k+1) + rhobf(k) * a_in(i,j,k)) *(w0(i,j,k+1)+ w0(i-1,j,k+1)) &
             ) / (4.*dzf(k)) &
             )
      end do
    end do

    !$acc parallel loop collapse(3) default(present) async(1)
    do k = 2, kmax
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
              +(1./rhobf(k))*( &
                (rhobf(k) * a_in(i,j,k) + rhobf(k+1) * a_in(i,j,k+1) )*(w0(i,j,k+1)+w0(i-1,j,k+1)) &
               -(rhobf(k) * a_in(i,j,k) + rhobf(k-1) * a_in(i,j,k-1) )*(w0(i,j,k  )+w0(i-1,j,k  )) &
                ) / (4. * dzf(k)) &
                )
        end do
      end do
    end do

  else ! non-equidistant grid

    k = 1
    !$acc parallel loop collapse(2) default(present) async(1)
    do j = 2, j1
      do i = 2, i1
        a_out(i,j,k) = a_out(i,j,k) - ( &
               (&
               (u0(i+1,j,k)+u0(i,j,k))/60&
               *(37*(a_in(i+1,j,k)+a_in(i,j,k))-8*(a_in(i+2,j,k)+a_in(i-1,j,k))+(a_in(i+3,j,k)+a_in(i-2,j,k)))&
               -sign(1._field_r,(u0(i+1,j,k)+u0(i,j,k)))*(u0(i+1,j,k)+u0(i,j,k))/60&
               *(10*(a_in(i+1,j,k)-a_in(i,j,k))-5*(a_in(i+2,j,k)-a_in(i-1,j,k))+(a_in(i+3,j,k)-a_in(i-2,j,k)))&
               -(u0(i,j,k)+u0(i-1,j,k))/60&
               *(37*(a_in(i,j,k)+a_in(i-1,j,k))-8*(a_in(i+1,j,k)+a_in(i-2,j,k))+(a_in(i+2,j,k)+a_in(i-3,j,k)))&
               +sign(1._field_r,(u0(i,j,k)+u0(i-1,j,k)))*(u0(i,j,k)+u0(i-1,j,k))/60&
               *(10*(a_in(i,j,k)-a_in(i-1,j,k))-5*(a_in(i+1,j,k)-a_in(i-2,j,k))+(a_in(i+2,j,k)-a_in(i-3,j,k)))&
               )*dxi5 &
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
               +(1./rhobf(k))*( &
               ( rhobf(k+1) * a_in(i,j,k+1)*dzf(k)   + rhobf(k)   * a_in(i,j,k)  *dzf(k+1) ) / dzh(k+1)  *( w0(i,j,k+1)+ w0(i-1,j,k+1) ) &
               ) / (4.*dzf(k)) &
               )
       end do
    end do

    !$acc parallel loop collapse(3) default(present) async(1)
    do k = 2, kmax
      do j = 2, j1
        do i = 2, i1
          a_out(i,j,k) = a_out(i,j,k) - ( &
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
                +(1./rhobf(k))*( &
                ( rhobf(k+1) * a_in(i,j,k+1)*dzf(k)   + rhobf(k)   * a_in(i,j,k)  *dzf(k+1) ) / dzh(k+1)  *( w0(i,j,k+1)+ w0(i-1,j,k+1) ) &
                -( rhobf(k)  * a_in(i,j,k)  *dzf(k-1) + rhobf(k-1) * a_in(i,j,k-1)*dzf(k)   ) / dzh(k)    *( w0(i,j,k)  + w0(i-1,j,k)   ) &
                ) / (4. * dzf(k)) &
                )
        end do
      end do
    end do
  end if
end subroutine advecu_52

!> Advection at the v point.
subroutine advecv_52(a_in, a_out)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi5,dyi5,dzf,dzh,leq
  use modfields, only : u0, v0, w0,rhobf
  implicit none

  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: a_in !< Input: the v field
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: a_out !< Output: the tendency
  integer :: i,j,k

  if (leq) then

    k = 1
    !$acc parallel loop collapse(2) default(present) async(2)
    do j = 2, j1
      do i = 2, i1
        a_out(i,j,k)  = a_out(i,j,k)- ( &
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
             +(1/rhobf(k))*( &
             (w0(i,j,k+1)+w0(i,j-1,k+1)) *(rhobf(k+1) * a_in(i,j,k+1) + rhobf(k) * a_in(i,j,k)) &
             ) / (4 * dzf(k)) &
             )
      end do
    end do

    !$acc parallel loop collapse(3) default(present) async(2)
    do k = 2, kmax
      do j = 2, j1
        do i = 2, i1
          a_out(i,j,k)  = a_out(i,j,k)- ( &
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
                +(1./rhobf(k))*( &
                  (w0(i,j,k+1)+w0(i,j-1,k+1))*(rhobf(k+1) * a_in(i,j,k+1) + rhobf(k) * a_in(i,j,k)) &
                  -(w0(i,j,k) +w0(i,j-1,k))  *(rhobf(k-1) * a_in(i,j,k-1) + rhobf(k) * a_in(i,j,k)) &
                  ) / (4. * dzf(k)) &
                  )

        end do
      end do
    end do

  else ! non-equidistant grid
    k = 1
    !$acc parallel loop collapse(2) default(present) async(2)
    do j = 2, j1
      do i = 2, i1
        a_out(i,j,k) = a_out(i,j,k) - ( &
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
               +(1./rhobf(k))*( &
               (w0(i,j,k+1)+w0(i,j-1,k+1)) * (rhobf(k+1) * a_in(i,j,k+1)*dzf(k) + rhobf(k) * a_in(i,j,k)*dzf(k+1)) / dzh(k+1) &
               ) / (4. * dzf(k)) &
               )
       end do
    end do

    !$acc parallel loop collapse(3) default(present) async(2)
    do k = 2, kmax
      do j = 2, j1
        do i = 2, i1
          a_out(i,j,k)  = a_out(i,j,k)- ( &
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
                    -sign(1._field_r,(v0(i,j+1,k)+v0(i,j,k)))*(v0(i,j+1,k)+v0(i,j+1,k))/60&
                    *(10*(a_in(i,j+1,k)-a_in(i,j,k))-5*(a_in(i,j+2,k)-a_in(i,j-1,k))+(a_in(i,j+3,k)-a_in(i,j-2,k)))&
                    -(v0(i,j,k)+v0(i,j-1,k))/60&
                    *(37*(a_in(i,j,k)+a_in(i,j-1,k))-8*(a_in(i,j+1,k)+a_in(i,j-2,k))+(a_in(i,j+2,k)+a_in(i,j-3,k)))&
                    +sign(1._field_r,(v0(i,j,k)+v0(i,j-1,k)))*(v0(i,j,k)+v0(i,j-1,k))/60&
                    *(10*(a_in(i,j,k)-a_in(i,j-1,k))-5*(a_in(i,j+1,k)-a_in(i,j-2,k))+(a_in(i,j+2,k)-a_in(i,j-3,k)))&
                  )* dyi5 &
                +(1./rhobf(k))*( &
                (w0(i,j,k+1)+w0(i,j-1,k+1)) * (rhobf(k+1) * a_in(i,j,k+1)*dzf(k) + rhobf(k) * a_in(i,j,k)*dzf(k+1)) / dzh(k+1) &
               -(w0(i,j,k)  +  w0(i,j-1,k)) * (rhobf(k-1) * a_in(i,j,k-1)*dzf(k) + rhobf(k) * a_in(i,j,k)*dzf(k-1)) / dzh(k) &
                ) / (4. * dzf(k)) &
                )
          ! note advec_2nd had rhobf(k) instead of rhobf(k+1) on top row, which seems wrong
          ! fixed here
        end do
      end do
    end do
  end if
end subroutine advecv_52

!> Advection at the w point.
subroutine advecw_52(a_in, a_out)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi5,dyi5,dzh,dzh,leq
  use modfields, only : u0, v0, w0,rhobh
  implicit none

  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: a_in !< Input: the w field
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: a_out !< Output: the tendency

  integer :: i,j,k
  ! if (leq) then
  ! FJ: judging from advec_2nd, equidistant and non-equidistant cases are similar

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
            + (1./rhobh(k))*( &
            (rhobh(k) * a_in(i,j,k) + rhobh(k+1) * a_in(i,j,k+1) )*(w0(i,j,k) + w0(i,j,k+1)) &
              -(rhobh(k) * a_in(i,j,k) + rhobh(k-1) * a_in(i,j,k-1) )*(w0(i,j,k) + w0(i,j,k-1)) &
              )/ (4. * dzh(k)) &
              )
      end do
    end do
  end do
  !Advection for u, v and w called sequentially in modavection. Only sync here.
  !$acc wait(1,2,3)
end subroutine advecw_52

end module advec_52
