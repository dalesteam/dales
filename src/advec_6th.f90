!> \file advec_6th.f90
!!  Does advection with a 6th order central differencing scheme.
!! \par Revision list
!! \par Authors
!! \see Wicker and Scamarock 2002
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

module advec_6th
use modprecision, only : field_r
contains
!> Advection at cell center
subroutine advecc_6th(a_in, a_out, istart, iend, jstart, jend)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi,dyi,dzf
  use modfields, only : u0, v0, w0, rhobf

  implicit none

  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in) :: a_in
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: a_out
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1) :: rho_a_in
  integer, intent(in) :: istart, iend, jstart ,jend

  integer :: i,j,k

!   if (leq) then

  do k=1,k1
    do j=jstart-jh,jend+jh
      do i=istart-ih,iend+ih
      rho_a_in(i,j,k)=rhobf(k)*a_in(i,j,k)
      end do
    end do
  end do

  do k=1,kmax
    do j=jstart,jend
      do i=istart,iend

        if(k==1) then

          a_out(i,j,k)  = a_out(i,j,k)- ( &
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
                +(1./rhobf(1))*( &
                w0(i,j,k+1) * (rho_a_in(i,j,k+1) + rho_a_in(i,j,k)) / 2. &
                ) / dzf(k) &
                )
       
        elseif(k==2) then
        !CvH do 2nd order for influx and 4th order for outflux

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
                +(1./rhobf(k))*( &
                      w0(i,j,k+1) / 12. &
                      *(7.*(rho_a_in(i,j,k+1)+rho_a_in(i,j,k))-(rho_a_in(i,j,k+2)+rho_a_in(i,j,k-1)))&
                      - w0(i,j,k)  * (rho_a_in(i,j,k-1)+rho_a_in(i,j,k)) / 2. &
                  ) / dzf(k) &
                  )

        elseif(k == 3) then
        !CvH do 6th order for outflux and 4th for influx

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
                +(1./rhobf(k))*( &
                      w0(i,j,k+1)/60. &
                      *(37.*(rho_a_in(i,j,k+1)+rho_a_in(i,j,k))-8.*(rho_a_in(i,j,k+2)+rho_a_in(i,j,k-1))&
                           +(rho_a_in(i,j,k+3)+rho_a_in(i,j,k-2)))&
                      -w0(i,j,k)/12. &
                      *(7.*(rho_a_in(i,j,k)+rho_a_in(i,j,k-1))-(rho_a_in(i,j,k+1)+rho_a_in(i,j,k-2)))&
                  )/dzf(k) &
                  )

        elseif(k==kmax-1) then
        !CvH do 6th order for influx and 4th order for outflux

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
                +(1./rhobf(k))*( &
                      w0(i,j,k+1)/12. &
                      *(7.*(rho_a_in(i,j,k+1)+rho_a_in(i,j,k))-(rho_a_in(i,j,k+2)+rho_a_in(i,j,k-1)))&
                      -w0(i,j,k)/60. &
                      *(37.*(rho_a_in(i,j,k)+rho_a_in(i,j,k-1))-8.*(rho_a_in(i,j,k+1)+rho_a_in(i,j,k-2))&
                          +(rho_a_in(i,j,k+2)+rho_a_in(i,j,k-3)))&
                  )/dzf(k) &
                  )
 
        elseif(k==kmax) then
        !CvH do 4th order for influx and 2nd order for outflux
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
                +(1./rhobf(k))*( &
                     w0(i,j,k+1) * (rho_a_in(i,j,k+1)+rho_a_in(i,j,k)) / 2. &
                     -w0(i,j,k)/12. &
                     *(7.*(rho_a_in(i,j,k)+rho_a_in(i,j,k-1))-(rho_a_in(i,j,k+1)+rho_a_in(i,j,k-2)))&
                  )/dzf(k) &
                  )

        else

            a_out(i,j,k)  = a_out(i,j,k)- (  &
                  ( &
                      u0(i+1,j,k)/60. &
                      *(37.*(a_in(i+1,j,k)+a_in(i,j,k))-8.*(a_in(i+2,j,k)+a_in(i-1,j,k))+(a_in(i+3,j,k)+a_in(i-2,j,k)))&
                      -u0(i,j,k)/60. &
                      *(37.*(a_in(i,j,k)+a_in(i-1,j,k))-8.*(a_in(i+1,j,k)+a_in(i-2,j,k))+(a_in(i+2,j,k)+a_in(i-3,j,k)))&
                  )* dxi &
                +( &
                      v0(i,j+1,k)/60. &
                      *(37.*(a_in(i,j+1,k)+a_in(i,j,k))-8.*(a_in(i,j+2,k)+a_in(i,j-1,k))+(a_in(i,j+3,k)+a_in(i,j-2,k)))&
                      -v0(i,j,k)/60. &
                      *(37.*(a_in(i,j,k)+a_in(i,j-1,k))-8.*(a_in(i,j+1,k)+a_in(i,j-2,k))+(a_in(i,j+2,k)+a_in(i,j-3,k)))&
                  )*dyi&
                +(1./rhobf(k))*(&
                      w0(i,j,k+1)/60. &
                      *(37.*(rho_a_in(i,j,k+1)+rho_a_in(i,j,k))-8.*(rho_a_in(i,j,k+2)+rho_a_in(i,j,k-1))&
                          +(rho_a_in(i,j,k+3)+rho_a_in(i,j,k-2)))&
                      -w0(i,j,k)/60. &
                      *(37.*(rho_a_in(i,j,k)+rho_a_in(i,j,k-1))-8.*(rho_a_in(i,j,k+1)+rho_a_in(i,j,k-2))&
                          +(rho_a_in(i,j,k+2)+rho_a_in(i,j,k-3)))&
                  )/dzf(k) &
                  )
        
        end if

      end do
    end do
  end do


!   end if

end subroutine advecc_6th


!> Advection at the u point.
subroutine advecu_6th(a_in,a_out,istart,iend,jstart,jend)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi5,dyi5,dzf
  use modfields, only : u0, v0, w0, rhobf

  implicit none
  integer, intent(in) :: istart, iend, jstart, jend
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: a_in !< Input: the u field
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: a_out !< Output: the tendency
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1) :: rho_a_in

  integer :: i,j,k


  do k=1,k1
    do j=jstart-jh,jend+jh
      do i=istart-ih,iend+ih
      rho_a_in(i,j,k)=rhobf(k)*a_in(i,j,k)
      end do
    end do
  end do

!   if (leq) then

    do k=1,kmax
      do j=jstart,jend
        do i=istart,iend

        if(k==1) then

          a_out(i,j,k)  = a_out(i,j,k)- ( &
                ( &
                    (u0(i+1,j,k)+u0(i,j,k))/60. &
                    *(37.*(a_in(i+1,j,k)+a_in(i,j,k))-8.*(a_in(i+2,j,k)+a_in(i-1,j,k))+(a_in(i+3,j,k)+a_in(i-2,j,k)))&
                    -(u0(i,j,k)+u0(i-1,j,k))/60. &
                    *(37.*(a_in(i,j,k)+a_in(i-1,j,k))-8.*(a_in(i+1,j,k)+a_in(i-2,j,k))+(a_in(i+2,j,k)+a_in(i-3,j,k)))&
                )*dxi5&
              +(&
                    (v0(i,j+1,k)+v0(i-1,j+1,k))/60. &
                    *(37.*(a_in(i,j+1,k)+a_in(i,j,k))-8.*(a_in(i,j+2,k)+a_in(i,j-1,k))+(a_in(i,j+3,k)+a_in(i,j-2,k)))&
                    -(v0(i,j,k)+v0(i-1,j,k))/60. &
                    *(37.*(a_in(i,j,k)+a_in(i,j-1,k))-8.*(a_in(i,j+1,k)+a_in(i,j-2,k))+(a_in(i,j+2,k)+a_in(i,j-3,k)))&
                )* dyi5 &
              +(1./rhobf(1))*( &
                    (rho_a_in(i,j,k+1) + rho_a_in(i,j,k)) *(w0(i,j,k+1)+w0(i-1,j,k+1)) / 2. &
                )/(2.*dzf(k)) &
                )

        elseif(k==2) then

          a_out(i,j,k)  = a_out(i,j,k)- ( &
                (&
                    (u0(i+1,j,k)+u0(i,j,k))/60. &
                    *(37.*(a_in(i+1,j,k)+a_in(i,j,k))-8.*(a_in(i+2,j,k)+a_in(i-1,j,k))+(a_in(i+3,j,k)+a_in(i-2,j,k)))&
                    -(u0(i,j,k)+u0(i-1,j,k))/60. &
                    *(37.*(a_in(i,j,k)+a_in(i-1,j,k))-8.*(a_in(i+1,j,k)+a_in(i-2,j,k))+(a_in(i+2,j,k)+a_in(i-3,j,k)))&
                )*dxi5&
              +(&
                    (v0(i,j+1,k)+v0(i-1,j+1,k))/60. &
                    *(37.*(a_in(i,j+1,k)+a_in(i,j,k))-8.*(a_in(i,j+2,k)+a_in(i,j-1,k))+(a_in(i,j+3,k)+a_in(i,j-2,k)))&
                    -(v0(i,j,k)+v0(i-1,j,k))/60. &
                    *(37.*(a_in(i,j,k)+a_in(i,j-1,k))-8.*(a_in(i,j+1,k)+a_in(i,j-2,k))+(a_in(i,j+2,k)+a_in(i,j-3,k)))&
                )* dyi5 &
              +(1./rhobf(k))*( &
                    (w0(i,j,k+1)+w0(i-1,j,k+1))/12. &
                    *(7.*(rho_a_in(i,j,k+1)+rho_a_in(i,j,k))-(rho_a_in(i,j,k+2)+rho_a_in(i,j,k-1)))&
                    -(rho_a_in(i,j,k)+rho_a_in(i,j,k-1))*(w0(i,j,k)+w0(i-1,j,k)) / 2. &
                ) / (2.*dzf(k)) &
                )

        elseif(k==3) then

          a_out(i,j,k)  = a_out(i,j,k)- ( &
                (&
                    (u0(i+1,j,k)+u0(i,j,k))/60. &
                    *(37.*(a_in(i+1,j,k)+a_in(i,j,k))-8.*(a_in(i+2,j,k)+a_in(i-1,j,k))+(a_in(i+3,j,k)+a_in(i-2,j,k)))&
                    -(u0(i,j,k)+u0(i-1,j,k))/60. &
                    *(37.*(a_in(i,j,k)+a_in(i-1,j,k))-8.*(a_in(i+1,j,k)+a_in(i-2,j,k))+(a_in(i+2,j,k)+a_in(i-3,j,k)))&
                )*dxi5&
              +(&
                    (v0(i,j+1,k)+v0(i-1,j+1,k))/60. &
                    *(37.*(a_in(i,j+1,k)+a_in(i,j,k))-8.*(a_in(i,j+2,k)+a_in(i,j-1,k))+(a_in(i,j+3,k)+a_in(i,j-2,k)))&
                    -(v0(i,j,k)+v0(i-1,j,k))/60. &
                    *(37.*(a_in(i,j,k)+a_in(i,j-1,k))-8.*(a_in(i,j+1,k)+a_in(i,j-2,k))+(a_in(i,j+2,k)+a_in(i,j-3,k)))&
                )* dyi5 &
              +(1./rhobf(k))*( &
                    (w0(i,j,k+1)+w0(i-1,j,k+1))/60. &
                    *(37.*(rho_a_in(i,j,k+1)+rho_a_in(i,j,k))-8.*(rho_a_in(i,j,k+2)+rho_a_in(i,j,k-1))&
                        +(rho_a_in(i,j,k+3)+rho_a_in(i,j,k-2)))&
                    -(w0(i,j,k)+w0(i-1,j,k))/12. &
                    *(7.*(rho_a_in(i,j,k)+rho_a_in(i,j,k-1))-(rho_a_in(i,j,k+1)+rho_a_in(i,j,k-2)))&
                ) / (2.*dzf(k)) &
                )

        elseif(k==kmax-1) then

          a_out(i,j,k)  = a_out(i,j,k)- ( &
                (&
                    (u0(i+1,j,k)+u0(i,j,k))/60. &
                    *(37.*(a_in(i+1,j,k)+a_in(i,j,k))-8.*(a_in(i+2,j,k)+a_in(i-1,j,k))+(a_in(i+3,j,k)+a_in(i-2,j,k)))&
                    -(u0(i,j,k)+u0(i-1,j,k))/60. &
                    *(37.*(a_in(i,j,k)+a_in(i-1,j,k))-8.*(a_in(i+1,j,k)+a_in(i-2,j,k))+(a_in(i+2,j,k)+a_in(i-3,j,k)))&
                )*dxi5&
              +(&
                    (v0(i,j+1,k)+v0(i-1,j+1,k))/60. &
                    *(37.*(a_in(i,j+1,k)+a_in(i,j,k))-8.*(a_in(i,j+2,k)+a_in(i,j-1,k))+(a_in(i,j+3,k)+a_in(i,j-2,k)))&
                    -(v0(i,j,k)+v0(i-1,j,k))/60. &
                    *(37.*(a_in(i,j,k)+a_in(i,j-1,k))-8.*(a_in(i,j+1,k)+a_in(i,j-2,k))+(a_in(i,j+2,k)+a_in(i,j-3,k)))&
                )* dyi5 &
              +(1./rhobf(k))*( &
                    (w0(i,j,k+1)+w0(i-1,j,k+1))/12. &
                    *(7.*(rho_a_in(i,j,k+1)+rho_a_in(i,j,k))-(rho_a_in(i,j,k+2)+rho_a_in(i,j,k-1)))&
                    -(w0(i,j,k)+w0(i-1,j,k))/60. &
                    *(37.*(rho_a_in(i,j,k)+rho_a_in(i,j,k-1))-8.*(rho_a_in(i,j,k+1)+rho_a_in(i,j,k-2))&
                        +(rho_a_in(i,j,k+2)+rho_a_in(i,j,k-3)))&
                ) / (2.*dzf(k)) &
                )
 
        elseif(k==kmax) then

          a_out(i,j,k)  = a_out(i,j,k)- ( &
                (&
                    (u0(i+1,j,k)+u0(i,j,k))/60. &
                    *(37.*(a_in(i+1,j,k)+a_in(i,j,k))-8.*(a_in(i+2,j,k)+a_in(i-1,j,k))+(a_in(i+3,j,k)+a_in(i-2,j,k)))&
                    -(u0(i,j,k)+u0(i-1,j,k))/60. &
                    *(37.*(a_in(i,j,k)+a_in(i-1,j,k))-8.*(a_in(i+1,j,k)+a_in(i-2,j,k))+(a_in(i+2,j,k)+a_in(i-3,j,k)))&
                )*dxi5&
              +(&
                    (v0(i,j+1,k)+v0(i-1,j+1,k))/60. &
                    *(37.*(a_in(i,j+1,k)+a_in(i,j,k))-8.*(a_in(i,j+2,k)+a_in(i,j-1,k))+(a_in(i,j+3,k)+a_in(i,j-2,k)))&
                    -(v0(i,j,k)+v0(i-1,j,k))/60. &
                    *(37.*(a_in(i,j,k)+a_in(i,j-1,k))-8.*(a_in(i,j+1,k)+a_in(i,j-2,k))+(a_in(i,j+2,k)+a_in(i,j-3,k)))&
                )* dyi5 &
              +(1./rhobf(k))*( &
                    (rho_a_in(i,j,k)+rho_a_in(i,j,k+1))*(w0(i,j,k+1)+w0(i-1,j,k+1))/2. &
                    -(w0(i,j,k)+w0(i-1,j,k))/12. &
                    *(7.*(rho_a_in(i,j,k)+rho_a_in(i,j,k-1))-(rho_a_in(i,j,k+1)+rho_a_in(i,j,k-2)))&
                ) / (2.*dzf(k)) &
                )
       
        else

          a_out(i,j,k)  = a_out(i,j,k)- ( &
                  ( &
                      (u0(i+1,j,k)+u0(i,j,k))/60. &
                      *(37.*(a_in(i+1,j,k)+a_in(i,j,k))-8.*(a_in(i+2,j,k)+a_in(i-1,j,k))+(a_in(i+3,j,k)+a_in(i-2,j,k)))&
                      -(u0(i,j,k)+u0(i-1,j,k))/60. &
                      *(37.*(a_in(i,j,k)+a_in(i-1,j,k))-8.*(a_in(i+1,j,k)+a_in(i-2,j,k))+(a_in(i+2,j,k)+a_in(i-3,j,k)))&
                  )*dxi5&
                +(&
                      (v0(i,j+1,k)+v0(i-1,j+1,k))/60. &
                      *(37.*(a_in(i,j+1,k)+a_in(i,j,k))-8.*(a_in(i,j+2,k)+a_in(i,j-1,k))+(a_in(i,j+3,k)+a_in(i,j-2,k)))&
                      -(v0(i,j,k)+v0(i-1,j,k))/60. &
                      *(37.*(a_in(i,j,k)+a_in(i,j-1,k))-8.*(a_in(i,j+1,k)+a_in(i,j-2,k))+(a_in(i,j+2,k)+a_in(i,j-3,k)))&
                  )*dyi5&
                + (1./rhobf(k))*(&
                      (w0(i,j,k+1)+w0(i-1,j,k+1))/60. &
                      *(37.*(rho_a_in(i,j,k+1)+rho_a_in(i,j,k))-8.*(rho_a_in(i,j,k+2)+rho_a_in(i,j,k-1))&
                          +(rho_a_in(i,j,k+3)+rho_a_in(i,j,k-2)))&
                      -(w0(i,j,k)+w0(i-1,j,k))/60. &
                      *(37.*(rho_a_in(i,j,k)+rho_a_in(i,j,k-1))-8.*(rho_a_in(i,j,k+1)+rho_a_in(i,j,k-2))&
                          +(rho_a_in(i,j,k+2)+rho_a_in(i,j,k-3)))&
                  ) / (2.*dzf(k)) &
                  )
        end if

        end do
      end do
    end do

!   end if

end subroutine advecu_6th

!> Advection at the v point.
subroutine advecv_6th(a_in, a_out, istart, iend, jstart, jend)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi5,dyi5,dzf
  use modfields, only : u0, v0, w0, rhobf
  implicit none
  integer, intent(in) :: istart, iend, jstart, jend
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: a_in !< Input: the v field
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: a_out !< Output: the tendency
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1) :: rho_a_in

  integer :: i,j,k

  do k=1,k1
    do j=jstart-jh,jend+jh
      do i=istart-ih,iend+ih
      rho_a_in(i,j,k)=rhobf(k)*a_in(i,j,k)
      end do
    end do
  end do

!   if (leq) then

    do k=1,kmax
      do j=jstart,jend
        do i=istart,iend

        if(k==1) then

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
                +(1./rhobf(1))*( &
                    (w0(i,j,k+1)+w0(i,j-1,k+1)) *(rho_a_in(i,j,k+1)+rho_a_in(i,j,k)) / 2.&
                  ) / (2.*dzf(k)) &
                  )

        elseif(k==2) then

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
                +(1./rhobf(k))*( &
                    (w0(i,j,k+1)+w0(i,j-1,k+1))/12. &
                    *(7.*(rho_a_in(i,j,k+1)+rho_a_in(i,j,k))-(rho_a_in(i,j,k+2)+rho_a_in(i,j,k-1)))&
                    -(w0(i,j,k)  +w0(i,j-1,k))  *(rho_a_in(i,j,k-1)+rho_a_in(i,j,k)) / 2.&
                  )/(2.*dzf(k)) &
                  )

        elseif(k==3) then

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
                +(1./rhobf(k))*( &
                    (w0(i,j,k+1)+w0(i,j-1,k+1))/60. &
                    *(37.*(rho_a_in(i,j,k+1)+rho_a_in(i,j,k))-8.*(rho_a_in(i,j,k+2)+rho_a_in(i,j,k-1))&
                        +(rho_a_in(i,j,k+3)+rho_a_in(i,j,k-2)))&
                    -(w0(i,j,k)+w0(i,j-1,k))/12. &
                    *(7.*(rho_a_in(i,j,k)+rho_a_in(i,j,k-1))-(rho_a_in(i,j,k+1)+rho_a_in(i,j,k-2)))&
                  )/(2.*dzf(k)) &
                  )

        elseif(k==kmax-1) then

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
                +(1./rhobf(k))*( &
                    (w0(i,j,k+1)+w0(i,j-1,k+1))/12. &
                    *(7.*(rho_a_in(i,j,k+1)+rho_a_in(i,j,k))-(rho_a_in(i,j,k+2)+rho_a_in(i,j,k-1)))&
                    -(w0(i,j,k)+w0(i,j-1,k))/60. &
                    *(37.*(rho_a_in(i,j,k)+rho_a_in(i,j,k-1))-8.*(rho_a_in(i,j,k+1)+rho_a_in(i,j,k-2))&
                        +(rho_a_in(i,j,k+2)+rho_a_in(i,j,k-3)))&
                  ) / (2.*dzf(k)) &
                  )

        elseif(k==kmax) then

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
                +(1./rhobf(k))*( &
                    (w0(i,j,k+1)+w0(i,j-1,k+1))*(rho_a_in(i,j,k+1)+rho_a_in(i,j,k)) / 2.&
                    -(w0(i,j,k)+w0(i,j-1,k))/12. &
                    *(7.*(rho_a_in(i,j,k)+rho_a_in(i,j,k-1))-(rho_a_in(i,j,k+1)+rho_a_in(i,j,k-2)))&
                  ) / (2.*dzf(k)) &
                  )

        else

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
                  )*dyi5&
                +(1./rhobf(k))*(&
                      (w0(i,j,k+1)+w0(i,j-1,k+1))/60. &
                      *(37.*(rho_a_in(i,j,k+1)+rho_a_in(i,j,k))-8.*(rho_a_in(i,j,k+2)+rho_a_in(i,j,k-1))&
                          +(rho_a_in(i,j,k+3)+rho_a_in(i,j,k-2)))&
                      -(w0(i,j,k)+w0(i,j-1,k))/60. &
                      *(37.*(rho_a_in(i,j,k)+rho_a_in(i,j,k-1))-8.*(rho_a_in(i,j,k+1)+rho_a_in(i,j,k-2))&
                          +(rho_a_in(i,j,k+2)+rho_a_in(i,j,k-3)))&
                  ) / (2.*dzf(k)) &
                  )

        end if

        end do
      end do
    end do


!   end if

end subroutine advecv_6th


!> Advection at the w point.
subroutine advecw_6th(a_in, a_out, istart, iend, jstart, jend)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi5,dyi5,dzh
  use modfields, only : u0, v0, w0, rhobh
  implicit none

  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: a_in !< Input: the w field
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: a_out !< Output: the tendency
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1) :: rho_a_in
  integer, intent(in) :: istart, iend, jstart, jend

  integer :: i,j,k


  do k=1,k1
    do j=jstart-jh,jend+jh
      do i=istart-ih,iend+ih
      rho_a_in(i,j,k)=rhobh(k)*a_in(i,j,k)
      end do
    end do
  end do

!   if (leq) then

    do k=2,kmax
      do j=jstart,jend
        do i=istart,iend

          if(k==2) then
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
                + (1./rhobh(k))*( &
                     (w0(i,j,k)+w0(i,j,k+1))/12. &
                     *(7.*(rho_a_in(i,j,k+1)+rho_a_in(i,j,k))-(rho_a_in(i,j,k+2)+rho_a_in(i,j,k-1)))&
                     -(rho_a_in(i,j,k)+rho_a_in(i,j,k-1))*(w0(i,j,k) + w0(i,j,k-1)) / 2. &
                  )/(2.*dzh(k)) &
                  )
 
          elseif(k==3) then
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
                + (1./rhobh(k))*( &
                     (w0(i,j,k)+w0(i,j,k+1))/60. &
                     *(37.*(rho_a_in(i,j,k+1)+rho_a_in(i,j,k))-8.*(rho_a_in(i,j,k+2)+rho_a_in(i,j,k-1))&
                         +(rho_a_in(i,j,k+3)+rho_a_in(i,j,k-2)))&
                     -(w0(i,j,k)+w0(i,j,k-1))/12. &
                     *(7.*(rho_a_in(i,j,k)+rho_a_in(i,j,k-1))-(rho_a_in(i,j,k+1)+rho_a_in(i,j,k-2))) &
                  )/(2.*dzh(k)) &
                  )

          elseif(k==kmax-1) then
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
                + (1./rhobh(k))*( &
                     (w0(i,j,k)+w0(i,j,k+1))/12. &
                     *(7.*(rho_a_in(i,j,k+1)+rho_a_in(i,j,k))-(rho_a_in(i,j,k+2)+rho_a_in(i,j,k-1))) &
                     -(w0(i,j,k)+w0(i,j,k-1))/60. &
                     *(37.*(rho_a_in(i,j,k)+rho_a_in(i,j,k-1))-8.*(rho_a_in(i,j,k+1)+rho_a_in(i,j,k-2))&
                         +(rho_a_in(i,j,k+2)+rho_a_in(i,j,k-3)))&
                  )/(2.*dzh(k)) &
                  )

          elseif(k==kmax) then
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
                + (1./rhobh(k))*( &
                     (rho_a_in(i,j,k)+rho_a_in(i,j,k+1) )*(w0(i,j,k) + w0(i,j,k+1)) / 2. &
                     -(w0(i,j,k)+w0(i,j,k-1))/12. &
                     *(7.*(rho_a_in(i,j,k)+rho_a_in(i,j,k-1))-(rho_a_in(i,j,k+1)+rho_a_in(i,j,k-2))) &
                  ) / (2.*dzh(k)) &
                  )
 
          else

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
                  )*dyi5&
                + (1./rhobh(k))*(&
                      (w0(i,j,k)+w0(i,j,k+1))/60. &
                      *(37.*(rho_a_in(i,j,k+1)+rho_a_in(i,j,k))-8.*(rho_a_in(i,j,k+2)+rho_a_in(i,j,k-1))&
                          +(rho_a_in(i,j,k+3)+rho_a_in(i,j,k-2)))&
                      -(w0(i,j,k)+w0(i,j,k-1))/60. &
                      *(37.*(rho_a_in(i,j,k)+rho_a_in(i,j,k-1))-8.*(rho_a_in(i,j,k+1)+rho_a_in(i,j,k-2))&
                          +(rho_a_in(i,j,k+2)+rho_a_in(i,j,k-3)))&
                  ) / (2.*dzh(k)) &
                  )
          end if
        end do
      end do
    end do

end subroutine advecw_6th
end module advec_6th
