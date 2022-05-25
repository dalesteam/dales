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
!> Advection at cell center
subroutine advecc_2nd(a_in,a_out)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi5,dyi5,dzi5,dzf,dzh,leq
  use modfields, only : u0, v0, w0, rhobf
  implicit none

  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: a_in !< Input: the cell centered field
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: a_out !< Output: the tendency
!  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1) :: rho_a_in

  integer :: i,j,k

!  do k=1,k1
!    do j=2-jh,j1+jh
!      do i=2-ih,i1+ih
!      rho_a_in(i,j,k)=rhobf(k)*a_in(i,j,k)
!      end do
!    end do
!  end do

  do k=1,kmax
    do j=2,j1
      do i=2,i1
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

  if (leq) then ! equidistant grid

    do j=2,j1
      do i=2,i1
        a_out(i,j,1)  = a_out(i,j,1)- (1./rhobf(1))*( &
                w0(i,j,2) * (rhobf(2) * a_in(i,j,2) + rhobf(1) * a_in(i,j,1) ) &
                ) * dzi5
      end do
    end do

    do j=2,j1
    do k=2,kmax         
       do i=2,i1
          a_out(i,j,k)  = a_out(i,j,k)- (1./rhobf(k))*( &
                w0(i,j,k+1) * (rhobf(k+1) * a_in(i,j,k+1) + rhobf(k) * a_in(i,j,k)) &
                -w0(i,j,k)   * (rhobf(k-1) * a_in(i,j,k-1)+ rhobf(k) * a_in(i,j,k)) &
                )*dzi5
        end do
      end do
    end do

  else   ! non-equidistant grid

    do j=2,j1
      do i=2,i1
        a_out(i,j,1)  = a_out(i,j,1)- (1./rhobf(1))*( &
                w0(i,j,2) * (rhobf(2) * a_in(i,j,2) * dzf(1) + rhobf(1) * a_in(i,j,1) * dzf(2) ) / (2.*dzh(2)) &
                ) / dzf(1)
      end do
    end do

    do j=2,j1
    do k=2,kmax
       do i=2,i1
          a_out(i,j,k)  = a_out(i,j,k)- (1./rhobf(k))*( &
                w0(i,j,k+1) * (rhobf(k+1) * a_in(i,j,k+1) * dzf(k) + rhobf(k) * a_in(i,j,k) * dzf(k+1) ) / dzh(k+1) &
               -w0(i,j,k ) * (rhobf(k-1) * a_in(i,j,k-1) * dzf(k) + rhobf(k) * a_in(i,j,k) * dzf(k-1) ) / dzh(k) &
                )/ (2. * dzf(k))
        end do
      end do
    end do

  end if

end subroutine advecc_2nd


!> Advection at the u point.
subroutine advecu_2nd(a_in, a_out)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxiq,dyiq,dziq,dzf,dzh,leq
  use modfields, only : u0, v0, w0, rhobf
  implicit none

  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in) :: a_in
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: a_out
!  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1) :: rho_a_in

  integer :: i,j,k,ip,im,jp,jm,kp,km

!  do k=1,k1
!    do j=2-jh,j1+jh
!      do i=2-ih,i1+ih
!      rho_a_in(i,j,k)=rhobf(k)*a_in(i,j,k)
!      end do
!    end do
!  end do

  do k=1,kmax
  km=k-1
  kp=k+1
    do j=2,j1
    jm=j-1
    jp=j+1
      do i=2,i1
      im=i-1
      ip=i+1
        a_out(i,j,k)  = a_out(i,j,k)- ( &
                ( &
                (a_in(i,j,k)+a_in(ip,j,k))*(u0(i,j,k)+u0(ip,j,k)) &
                -(a_in(i,j,k)+a_in(im,j,k))*(u0(i,j,k)+u0(im,j,k)) &
                )*dxiq &
                +(  &
                (a_in(i,j,k)+a_in(i,jp,k))*(v0(i,jp,k)+v0(im,jp ,k)) &
                -(a_in(i,j,k)+a_in(i,jm,k))*(v0(i,j  ,k)+v0(im,j  ,k)) &
                )*dyiq )

      end do
    end do
  end do

  if (leq) then

    do j=2,j1
    jm=j-1
    jp=j+1
      do i=2,i1
      im=i-1
      ip=i+1
        a_out(i,j,1)  = a_out(i,j,1)-(1./rhobf(1))*( &
            ( rhobf(2) * a_in(i,j,2) + rhobf(1) * a_in(i,j,1))*( w0(i,j,2)+ w0(im,j,2) ) &
            ) *dziq
      end do
    end do

    do j=2,j1
    jm=j-1
    jp=j+1
    do k=2,kmax
       km=k-1
       kp=k+1
       do i=2,i1
          im=i-1
          ip=i+1
          a_out(i,j,k)  = a_out(i,j,k)- (1./rhobf(k))*( &
              (rhobf(k) * a_in(i,j,k) + rhobf(kp) * a_in(i,j,kp) )*(w0(i,j,kp)+w0(im,j,kp)) &
              -(rhobf(k) * a_in(i,j,k) + rhobf(km) * a_in(i,j,km) )*(w0(i,j,k )+w0(im,j,k )) &
                  )*dziq
        end do
      end do
    end do

  else

    do j=2,j1
    jm=j-1
    jp=j+1
      do i=2,i1
      im=i-1
      ip=i+1
        a_out(i,j,1)  = a_out(i,j,1)- (1./rhobf(1))*( &
              ( rhobf(2) * a_in(i,j,2)*dzf(1) + rhobf(1) * a_in(i,j,1)*dzf(2) ) / dzh(2) &
                *( w0(i,j,2)+ w0(im,j,2) ))/ (4.*dzf(1))
      end do
    end do

    do j=2,j1
    jm=j-1
    jp=j+1
    do k=2,kmax
       km=k-1
       kp=k+1
       do i=2,i1
          im=i-1
          ip=i+1
          a_out(i,j,k)  = a_out(i,j,k)- (1./rhobf(k))*( &
                ( rhobf(kp) * a_in(i,j,kp)*dzf(k) + rhobf(k) * a_in(i,j,k)*dzf(kp) ) / dzh(kp) &
                  *( w0(i,j,kp)+ w0(im,j,kp) ) &
               -( rhobf(k) * a_in(i,j,k)*dzf(km) + rhobf(km) * a_in(i,j,km)*dzf(k) ) / dzh(k) &
                  *( w0(i,j,k)  + w0(im,j,k)   ) &
                )/ (4.*dzf(k))
        end do
      end do
    end do




  end if

end subroutine advecu_2nd


!> Advection at the v point.
subroutine advecv_2nd(a_in, a_out)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxiq,dyiq,dziq,dzf,dzh,leq
  use modfields, only : u0, v0, w0, rhobf
  implicit none

  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: a_in !< Input: the v-field
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: a_out !< Output: the tendency
!  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1) :: rho_a_in

  integer :: i,j,k,ip,im,jp,jm,kp,km

!  do k=1,k1
!    do j=2-jh,j1+jh
!      do i=2-ih,i1+ih
!      rho_a_in(i,j,k)=rhobf(k)*a_in(i,j,k)
!      end do
!    end do
!  end do

  do k=1,kmax
  km=k-1
  kp=k+1
    do j=2,j1
    jm=j-1
    jp=j+1
      do i=2,i1
      im=i-1
      ip=i+1

        a_out(i,j,k)  = a_out(i,j,k)- ( &
              ( &
              ( u0(ip,j,k)+u0(ip,jm,k))*(a_in(i,j,k)+a_in(ip,j,k)) &
              -(u0(i ,j,k)+u0(i ,jm,k))*(a_in(i,j,k)+a_in(im,j,k)) &
              )*dxiq &
              +( &
              ( v0(i,jp,k)+v0(i,j,k))*(a_in(i,j,k)+a_in(i,jp,k)) &
              -(v0(i,jm,k)+v0(i,j,k))*(a_in(i,j,k)+a_in(i,jm,k)) &
              )*dyiq )
      end do
    end do
  end do

  if (leq) then

    do j=2,j1
    jm=j-1
    jp=j+1
      do i=2,i1
      im=i-1
      ip=i+1
        a_out(i,j,1)  = a_out(i,j,1)- (1./rhobf(1))*( &
           (w0(i,j,2)+w0(i,jm,2))*(rhobf(2) * a_in(i,j,2)+rhobf(1) * a_in(i,j,1)) &
            )*dziq
      end do
    end do

    do j=2,j1
    jm=j-1
    jp=j+1
    do k=2,kmax
       km=k-1
       kp=k+1
       do i=2,i1
          im=i-1
          ip=i+1
          a_out(i,j,k)  = a_out(i,j,k)- (1./rhobf(k))*( &
                ( w0(i,j,kp)+w0(i,jm,kp))*(rhobf(kp) * a_in(i,j,kp) + rhobf(k) * a_in(i,j,k)) &
                -(w0(i,j,k) +w0(i,jm,k)) *(rhobf(km) * a_in(i,j,km) + rhobf(k) * a_in(i,j,k)) &
                )*dziq
        end do
      end do
    end do

  else
    do j=2,j1
    jm=j-1
    jp=j+1
      do i=2,i1
        im=i-1
        ip=i+1
        a_out(i,j,1)  = a_out(i,j,1)- (1./rhobf(1))*( &
          (w0(i,j,2)+w0(i,jm,2)) &
          *(rhobf(2) * a_in(i,j,2)*dzf(1) + rhobf(1) * a_in(i,j,1)*dzf(2) )/ dzh(2) &
          ) / (4. * dzf(1))
      end do
    end do

    do j=2,j1
    jm=j-1
    jp=j+1
    do k=2,kmax
       km=k-1
       kp=k+1
       do i=2,i1
          im=i-1
          ip=i+1
          a_out(i,j,k)  = a_out(i,j,k)- (1./rhobf(k))*( &
            (w0(i,j,kp)+w0(i,jm,kp)) &
            *(rhobf(kp) * a_in(i,j,kp)*dzf(k) + rhobf(k) * a_in(i,j,k)*dzf(kp) )/ dzh(kp) &
            -(w0(i,j,k)+w0(i,jm,k)) &
            *(rhobf(km) * a_in(i,j,km)*dzf(k) + rhobf(k) * a_in(i,j,k)*dzf(km)) / dzh(k) &
            ) / (4. * dzf(k))
        end do
      end do
    end do

  end if

end subroutine advecv_2nd



!> Advection at the w point.
subroutine advecw_2nd(a_in,a_out)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxiq,dyiq,dziq,dzf,dzh,leq
  use modfields, only : u0, v0, w0, rhobh
  implicit none

  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: a_in !< Input: the w-field
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: a_out !< Output: the tendency
!  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1) :: rho_a_in

  integer :: i,j,k,ip,im,jp,jm,kp,km

!  do k=1,k1
!    do j=2-jh,j1+jh
!      do i=2-ih,i1+ih
!      rho_a_in(i,j,k)=rhobh(k)*a_in(i,j,k)
!      end do
!    end do
!  end do

  if (leq) then


    do k=2,kmax
    km=k-1
    kp=k+1
      do j=2,j1
      jm=j-1
      jp=j+1
        do i=2,i1
        im=i-1
        ip=i+1

          a_out(i,j,k)  = a_out(i,j,k)- ( &
                ( &
                (a_in(ip,j,k)+a_in(i,j,k))*(u0(ip,j,k)+u0(ip,j,km)) &
              -(a_in(im,j,k)+a_in(i,j,k))*(u0(i  ,j,k)+u0(i  ,j,km)) &
                )*dxiq &
              + &
                ( &
                (a_in(i,jp,k)+a_in(i,j,k))*(v0(i,jp,k)+v0(i,jp,km)) &
              -(a_in(i,jm,k)+a_in(i,j,k))*(v0(i,j  ,k)+v0(i,j  ,km)) &
                )*dyiq &
              + &
                (1./rhobh(k))*( &
                (rhobh(k) * a_in(i,j,k) + rhobh(kp) * a_in(i,j,kp) )*(w0(i,j,k) + w0(i,j,kp)) &
               -(rhobh(k) * a_in(i,j,k) + rhobh(km) * a_in(i,j,km) )*(w0(i,j,k) + w0(i,j,km)) &
                )*dziq &
                )

        end do
      end do
    end do
  else
    do k=2,kmax
    km=k-1
    kp=k+1
      do j=2,j1
        jm=j-1
        jp=j+1
        do i=2,i1
        im=i-1
        ip=i+1

          a_out(i,j,k)  = a_out(i,j,k) - (1./rhobh(k))*( &
                ( &
                ( rhobh(k) * a_in(ip,j,k) + rhobh(k) * a_in(i,j,k) ) &
              *( dzf(km)*u0(ip,j,k) + dzf(k)*u0(ip,j,km) ) &
              -( rhobh(k) * a_in(i,j,k) + rhobh(k) * a_in(im,j,k) ) &
              *( dzf(km)*u0(i,j,k)+dzf(k)*u0(i ,j,km) ) &
                )*dxiq / dzh(k) &
              + &
                ( &
                ( rhobh(k) * a_in(i,jp,k) + rhobh(k) * a_in(i,j,k) ) &
              *( dzf(km)*v0(i,jp,k) + dzf(k)*v0(i,jp,km) ) &
              -( rhobh(k) * a_in(i,j,k) + rhobh(k) * a_in(i,j-1,k) ) &
              *( dzf(km)*v0(i,j,k) + dzf(k)*v0(i,j,km) ) &
                ) *dyiq / dzh(k) &
              + &
                ( &
                ( rhobh(k) * a_in(i,j,k) + rhobh(kp) * a_in(i,j,kp) ) * (w0(i,j,k) + w0(i,j,kp) ) &
               -( rhobh(k) * a_in(i,j,k) + rhobh(km) * a_in(i,j,km) ) * (w0(i,j,k) + w0(i,j,km) ) &
                ) / (4. *dzh(k) ) &
                )

        end do
      end do
    end do

  end if


end subroutine advecw_2nd
end module advec_2nd
