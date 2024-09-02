!> \file advec_2nd.f90
!!  Does advection with a 2nd order central differencing scheme.
!! \par Revision list
!! \par Authors
!! Second order central differencing can be used for variables where neither very
!! high accuracy nor strict monotonicity is necessary.
!! \latexonly
!!\begin{eqnarray}
!! F_{i-\frac{1}{2}}^{2nd} &=&
!!\fav{u}_{i-\frac{1}{2}}\frac{\phi_{i}+\phi_{i-1}}{2},
!!\end{eqnarray}
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

module advec_2nd
use modprecision, only : field_r
contains
!> Horizontal advection at cell center
subroutine hadvecc_2nd(a_in,a_out,istart,iend,jstart,jend)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi5,dyi5
  use modfields, only : u0, v0
  implicit none

  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: a_in !< Input: the cell centered field
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: a_out !< Output: the tendency
  integer, intent(in) :: istart,iend,jstart,jend !< Input: start and end indices for advection routine

  integer :: i,j,k

  !$acc parallel loop collapse(3) default(present) async(1)
  do k = 1, kmax
    do j = jstart, jend
      do i = istart, iend
        a_out(i,j,k)  = a_out(i,j,k)- (  &
              ( &
              u0(i+1,j,k) * ( a_in(i+1,j,k) + a_in(i,j,k) ) &
             -u0(i ,j,k) * ( a_in(i-1,j,k) + a_in(i,j,k) ) &
              )* dxi5 &
            +( &
              v0(i,j+1,k) * ( a_in(i,j+1,k) + a_in(i,j,k) ) &
             -v0(i,j ,k) * ( a_in(i,j-1,k) + a_in(i,j,k) ) &
              )* dyi5 )
      end do
    end do
  end do
  !$acc wait

end subroutine hadvecc_2nd

!> Vertical advection at cell center
subroutine vadvecc_2nd(a_in,a_out,istart,iend,jstart,jend)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dzi5,dzf,dzfi,dzhi,leq
  use modfields, only : w0, rhobf
  implicit none

  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: a_in !< Input: the cell centered field
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: a_out !< Output: the tendency
  integer, intent(in) :: istart,iend,jstart,jend !< Input: start and end indices for advection routine

  integer :: i,j,k

  if (leq) then ! equidistant grid
    !$acc parallel loop collapse(2) default(present) wait(1) async(2)
    do j = jstart, jend
      do i = istart, iend
        a_out(i,j,1) = a_out(i,j,1)- (1/rhobf(1))*( &
                       w0(i,j,2) * (rhobf(2) * a_in(i,j,2) + rhobf(1) * a_in(i,j,1) ) &
                       ) * dzi5
      end do
    end do

    !$acc parallel loop collapse(3) default(present) wait(1) async(3)
    do k = 2, kmax
      do j = jstart, jend
        do i = istart, iend
          a_out(i,j,k) = a_out(i,j,k)- (1/rhobf(k))*( &
                         w0(i,j,k+1) * (rhobf(k+1) * a_in(i,j,k+1) + rhobf(k) * a_in(i,j,k)) &
                         -w0(i,j,k)   * (rhobf(k-1) * a_in(i,j,k-1)+ rhobf(k) * a_in(i,j,k)) &
                         )*dzi5
        end do
      end do
    end do

  else   ! non-equidistant grid
    !$acc parallel loop collapse(2) default(present) wait(1) async(2)
    do j = jstart, jend
      do i = istart, iend
        a_out(i,j,1) = a_out(i,j,1)- (1/rhobf(1))*( &
                       w0(i,j,2) * (rhobf(2) * a_in(i,j,2) * dzf(1) + rhobf(1) * a_in(i,j,1) * dzf(2) ) * (0.5_field_r * dzhi(2)) &
                       ) * dzfi(1)
      end do
    end do

    !$acc parallel loop collapse(3) default(present) wait(1) async(3)
    do k = 2, kmax
      do j = jstart, jend
        do i = istart, iend
          a_out(i,j,k) = a_out(i,j,k)- (1/rhobf(k))*( &
                         w0(i,j,k+1) * (rhobf(k+1) * a_in(i,j,k+1) * dzf(k) + rhobf(k) * a_in(i,j,k) * dzf(k+1) ) * dzhi(k+1) &
                        -w0(i,j,k ) * (rhobf(k-1) * a_in(i,j,k-1) * dzf(k) + rhobf(k) * a_in(i,j,k) * dzf(k-1) ) * dzhi(k) &
                         ) * (0.5_field_r * dzfi(k))
        end do
      end do
    end do
  end if
  !$acc wait

end subroutine vadvecc_2nd


!> Horizontal advection at the u point.
subroutine hadvecu_2nd(a_in, a_out,istart,iend,jstart,jend)
  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxiq,dyiq
  use modfields, only : u0, v0, up, vp

  implicit none

  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in) :: a_in
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: a_out
  integer, intent(in) :: istart,iend,jstart,jend !< Input: start and end indices for advection routine

  integer :: i,j,k,ip,im,jp,jm,kp,km

  !$acc parallel loop collapse(3) default(present) async(1)
  do k = 1, kmax
    do j = jstart, jend
      do i = istart, iend
        a_out(i,j,k)  = a_out(i,j,k)- ( &
                ( &
                (a_in(i,j,k)+a_in(i+1,j,k))*(u0(i,j,k)+u0(i+1,j,k)) &
                -(a_in(i,j,k)+a_in(i-1,j,k))*(u0(i,j,k)+u0(i-1,j,k)) &
                )*dxiq &
                +(  &
                (a_in(i,j,k)+a_in(i,j+1,k))*(v0(i,j+1,k)+v0(i-1,j+1 ,k)) &
                -(a_in(i,j,k)+a_in(i,j-1,k))*(v0(i,j  ,k)+v0(i-1,j  ,k)) &
                )*dyiq )
      end do
    end do
  end do

end subroutine hadvecu_2nd

!> Vertical advection at the u point.
subroutine vadvecu_2nd(a_in, a_out,istart,iend,jstart,jend)
  use modglobal, only : i1,ih,j1,jh,k1,kmax,dziq,dzf,dzfi,dzhi,leq
  use modfields, only : w0, rhobf

  implicit none

  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in) :: a_in
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: a_out
  integer, intent(in) :: istart,iend,jstart,jend !< Input: start and end indices for advection routine

  integer :: i,j,k,ip,im,jp,jm,kp,km

  if (leq) then

    !$acc parallel loop collapse(2) default(present) async(1)
    do j = jstart, jend
      do i = istart, iend
        a_out(i,j,1) = a_out(i,j,1)-(1/rhobf(1))*( &
                       ( rhobf(2) * a_in(i,j,2) + rhobf(1) * a_in(i,j,1))*( w0(i,j,2)+ w0(i-1,j,2) ) &
                       ) *dziq
      end do
    end do

    !$acc parallel loop collapse(3) default(present) async(1)
    do k = 2, kmax
      do j = jstart, jend
        do i = istart, iend
          a_out(i,j,k) = a_out(i,j,k)- (1/rhobf(k))*( &
                         (rhobf(k) * a_in(i,j,k) + rhobf(k+1) * a_in(i,j,k+1) )*(w0(i,j,k+1)+w0(i-1,j,k+1)) &
                         -(rhobf(k) * a_in(i,j,k) + rhobf(k-1) * a_in(i,j,k-1) )*(w0(i,j,k )+w0(i-1,j,k )) &
                             )*dziq
        end do
      end do
    end do

  else

    !$acc parallel loop collapse(2) default(present) async(1)
    do j = jstart, jend
      do i = istart, iend
        a_out(i,j,1) = a_out(i,j,1)- (1/rhobf(1))*( &
                       ( rhobf(2) * a_in(i,j,2)*dzf(1) + rhobf(1) * a_in(i,j,1)*dzf(2) ) * dzhi(2) &
                         *( w0(i,j,2)+ w0(i-1,j,2) )) * (0.25_field_r * dzfi(1))
      end do
    end do

    !$acc parallel loop collapse(3) default(present) async(1)
    do k = 2, kmax
      do j = jstart, jend
        do i = istart, iend
          a_out(i,j,k)  = a_out(i,j,k)- (1/rhobf(k))*( &
                ( rhobf(k+1) * a_in(i,j,k+1)*dzf(k) + rhobf(k) * a_in(i,j,k)*dzf(k+1) ) * dzhi(k+1) &
                  *( w0(i,j,k+1)+ w0(i-1,j,k+1) ) &
               -( rhobf(k) * a_in(i,j,k)*dzf(k-1) + rhobf(k-1) * a_in(i,j,k-1)*dzf(k) ) * dzhi(k) &
                  *( w0(i,j,k)  + w0(i-1,j,k)   ) &
                ) * (0.25_field_r*dzfi(k))
        end do
      end do
    end do
  end if
end subroutine vadvecu_2nd

!> Horizontal advection at the v point.
subroutine hadvecv_2nd(a_in, a_out,istart,iend,jstart,jend)
  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxiq,dyiq
  use modfields, only : u0, v0
  implicit none


  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: a_in !< Input: the v-field
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: a_out !< Output: the tendency
  integer, intent(in) :: istart,iend,jstart,jend !< Input: start and end indices for advection routine

  integer :: i,j,k,ip,im,jp,jm,kp,km

  !$acc parallel loop collapse(3) default(present) async(2)
  do k = 1, kmax
    do j = jstart, jend
      do i = istart, iend
        a_out(i,j,k)  = a_out(i,j,k)- ( &
              ( &
              ( u0(i+1,j,k)+u0(i+1,j-1,k))*(a_in(i,j,k)+a_in(i+1,j,k)) &
              -(u0(i ,j,k)+u0(i ,j-1,k))*(a_in(i,j,k)+a_in(i-1,j,k)) &
              )*dxiq &
              +( &
              ( v0(i,j+1,k)+v0(i,j,k))*(a_in(i,j,k)+a_in(i,j+1,k)) &
              -(v0(i,j-1,k)+v0(i,j,k))*(a_in(i,j,k)+a_in(i,j-1,k)) &
              )*dyiq )
      end do
    end do
  end do

end subroutine hadvecv_2nd

!> Vertical advection at the v point.
subroutine vadvecv_2nd(a_in, a_out,istart,iend,jstart,jend)
  use modglobal, only : i1,ih,j1,jh,k1,kmax,dziq,dzf,dzfi,dzhi,leq
  use modfields, only : w0, rhobf
  implicit none


  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: a_in !< Input: the v-field
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: a_out !< Output: the tendency
  integer, intent(in) :: istart,iend,jstart,jend !< Input: start and end indices for advection routine

  integer :: i,j,k,ip,im,jp,jm,kp,km

  if (leq) then
    !$acc parallel loop collapse(2) default(present) async(2)
    do j = jstart, jend
      do i = istart, iend
        a_out(i,j,1)  = a_out(i,j,1)- (1/rhobf(1))*( &
           (w0(i,j,2)+w0(i,j-1,2))*(rhobf(2) * a_in(i,j,2)+rhobf(1) * a_in(i,j,1)) &
            )*dziq
      end do
    end do

    !$acc parallel loop collapse(3) default(present) async(2)
    do k = 2, kmax
      do j = jstart, jend
        do i = istart, iend
          a_out(i,j,k)  = a_out(i,j,k)- (1/rhobf(k))*( &
                ( w0(i,j,k+1)+w0(i,j-1,k+1))*(rhobf(k+1) * a_in(i,j,k+1) + rhobf(k) * a_in(i,j,k)) &
                -(w0(i,j,k) +w0(i,j-1,k)) *(rhobf(k-1) * a_in(i,j,k-1) + rhobf(k) * a_in(i,j,k)) &
                )*dziq
        end do
      end do
    end do

  else
    !$acc parallel loop collapse(2) default(present) async(2)
    do j = jstart, jend
      do i = istart, iend
        a_out(i,j,1)  = a_out(i,j,1)- (1/rhobf(1))*( &
          (w0(i,j,2)+w0(i,j-1,2)) &
          *(rhobf(2) * a_in(i,j,2)*dzf(1) + rhobf(1) * a_in(i,j,1)*dzf(2) ) * dzhi(2) &
          ) * (0.25_field_r * dzfi(1))
      end do
    end do

    !$acc parallel loop collapse(3) default(present) async(2)
    do k = 2, kmax
      do j = jstart, jend
        do i = istart, iend
          a_out(i,j,k)  = a_out(i,j,k)- (1/rhobf(k))*( &
            (w0(i,j,k+1)+w0(i,j-1,k+1)) &
            *(rhobf(k+1) * a_in(i,j,k+1)*dzf(k) + rhobf(k) * a_in(i,j,k)*dzf(k+1) ) * dzhi(k+1) &
            -(w0(i,j,k)+w0(i,j-1,k)) &
            *(rhobf(k-1) * a_in(i,j,k-1)*dzf(k) + rhobf(k) * a_in(i,j,k)*dzf(k-1)) * dzhi(k) &
            ) * (0.25_field_r * dzfi(k))
        end do
      end do
    end do
  end if

end subroutine vadvecv_2nd

!> Horiozntal advection at the w point.
subroutine hadvecw_2nd(a_in,a_out,istart,iend,jstart,jend)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxiq,dyiq,dziq,dzf,dzhi,leq
  use modfields, only : u0, v0, rhobh
  implicit none

  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: a_in !< Input: the w-field
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: a_out !< Output: the tendency
  integer, intent(in) :: istart,iend,jstart,jend !< Input: start and end indices for advection routine

  integer :: i,j,k,ip,im,jp,jm,kp,km

  if (leq) then

    !$acc parallel loop collapse(3) default(present) async(3)
    do k = 2, kmax
      do j = jstart, jend
        do i = istart, iend
          a_out(i,j,k)  = a_out(i,j,k)- ( &
                ( &
                (a_in(i+1,j,k)+a_in(i,j,k))*(u0(i+1,j,k)+u0(i+1,j,k-1)) &
               -(a_in(i-1,j,k)+a_in(i,j,k))*(u0(i  ,j,k)+u0(i  ,j,k-1)) &
                )*dxiq &
              + &
                ( &
                (a_in(i,j+1,k)+a_in(i,j,k))*(v0(i,j+1,k)+v0(i,j+1,k-1)) &
               -(a_in(i,j-1,k)+a_in(i,j,k))*(v0(i,j  ,k)+v0(i,j  ,k-1)) &
                )*dyiq &
                )
        end do
      end do
    end do
  else
    !$acc parallel loop collapse(3) default(present) async(3)
    do k = 2, kmax
      do j = jstart, jend
        do i = istart, iend
          a_out(i,j,k)  = a_out(i,j,k) - (1/rhobh(k))*( &
                ( &
                ( rhobh(k) * a_in(i+1,j,k) + rhobh(k) * a_in(i,j,k) ) &
              *( dzf(k-1)*u0(i+1,j,k) + dzf(k)*u0(i+1,j,k-1) ) &
              -( rhobh(k) * a_in(i,j,k) + rhobh(k) * a_in(i-1,j,k) ) &
              *( dzf(k-1)*u0(i,j,k)+dzf(k)*u0(i ,j,k-1) ) &
                )*dxiq * dzhi(k) &
              + &
                ( &
                ( rhobh(k) * a_in(i,j+1,k) + rhobh(k) * a_in(i,j,k) ) &
              *( dzf(k-1)*v0(i,j+1,k) + dzf(k)*v0(i,j+1,k-1) ) &
              -( rhobh(k) * a_in(i,j,k) + rhobh(k) * a_in(i,j-1,k) ) &
              *( dzf(k-1)*v0(i,j,k) + dzf(k)*v0(i,j,k-1) ) &
                ) *dyiq * dzhi(k) &
                )
        end do
      end do
    end do
  end if
end subroutine hadvecw_2nd

!> Vertical advection at the w point.
subroutine vadvecw_2nd(a_in,a_out,istart,iend,jstart,jend)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dziq,dzf,dzhi,leq
  use modfields, only : w0, rhobh
  implicit none

  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: a_in !< Input: the w-field
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: a_out !< Output: the tendency
  integer, intent(in) :: istart,iend,jstart,jend !< Input: start and end indices for advection routine

  integer :: i,j,k,ip,im,jp,jm,kp,km

  if (leq) then

    !$acc parallel loop collapse(3) default(present) async(3)
    do k = 2, kmax
      do j = jstart, jend
        do i = istart, iend
          a_out(i,j,k)  = a_out(i,j,k)- ( &
                (1/rhobh(k))*( &
                (rhobh(k) * a_in(i,j,k) + rhobh(k+1) * a_in(i,j,k+1) )*(w0(i,j,k) + w0(i,j,k+1)) &
               -(rhobh(k) * a_in(i,j,k) + rhobh(k-1) * a_in(i,j,k-1) )*(w0(i,j,k) + w0(i,j,k-1)) &
                )*dziq &
                )
        end do
      end do
    end do
  else
    !$acc parallel loop collapse(3) default(present) async(3)
    do k = 2, kmax
      do j = jstart, jend
        do i = istart, iend
          a_out(i,j,k)  = a_out(i,j,k) - (1/rhobh(k))*( &
                ( &
                ( rhobh(k) * a_in(i,j,k) + rhobh(k+1) * a_in(i,j,k+1) ) * (w0(i,j,k) + w0(i,j,k+1) ) &
               -( rhobh(k) * a_in(i,j,k) + rhobh(k-1) * a_in(i,j,k-1) ) * (w0(i,j,k) + w0(i,j,k-1) ) &
                ) * (0.25_field_r * dzhi(k) ) &
                )
        end do
      end do
    end do
  end if
  !$acc wait(1,2,3)
end subroutine vadvecw_2nd

end module advec_2nd
