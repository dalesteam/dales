!> \file advec_weno.f90
!  Does advection with the 5th order advection scheme that was already present in DALES, except
!  around locations were discontinuities arise. There, the 5th order WENO scheme is used.
!  Discontinuities should be preserved, while the damping of high wavenumber components of the flow,
!  common in pure WENO solutions.
!  The scheme is more or less equal to the scheme proposed by Hill and Pullin (2004) but using 
!  the fifth order advection scheme instead of their tuned one.
!  JvdD (2011)
!
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
module advec_hybrid
implicit none
contains

subroutine advecc_hybrid(pin,pout)
  use modglobal, only : ih,i1,jh,j1,kmax,k1 &
                       ,dxi,dyi,dzf,lambda_max
  use modfields, only : u0,v0,w0
  use advec_weno,only : ip_weno
  implicit none

  real,dimension(2-ih:i1+ih,2-jh:j1+jh,k1),intent(in) :: pin  !< Input: the cell centered field (qt,thetal,sv etc)
  real,dimension(2-ih:i1+ih,2-jh:j1+jh,k1),intent(out):: pout !< Output: the tendency for the input field (qtp,thetalp,svp etc)

  ! local variables
  real,dimension(2-ih:i1+ih,2-jh:j1+jh,k1) :: x_smooth,y_smooth,z_smooth
  integer :: i,j,k
  integer :: kp3,kp2,kp1,km1,km2,km3
  integer :: jp3,jp2,jp1,jm1,jm2,jm3
  integer :: ip3,ip2,ip1,im1,im2,im3
  
  ! calculate the smoothness indicator in all 3 directions.
  x_smooth = smoothness(pin,1)
  y_smooth = smoothness(pin,2)
  z_smooth = smoothness(pin,3)

  do j=2,j1
    jp3=j+3;jp2=j+2;jp1=j+1;jm1=j-1;jm2=j-2;jm3=j-3
    do i=2,i1
      ip3=i+3;ip2=i+2;ip1=i+1;im1=i-1;im2=i-2;im3=i-3
      do k=1,kmax
        ! Use 5th order WENO in the horizontal directions
!        pout(i,j,k) = pout(i,j,k) - ( &
!                    (u0(ip1,j,k)*ip_weno(pin(ip3,j,k),pin(ip2,j,k),pin(ip1,j,k),pin(i,j,k),pin(im1,j,k),pin(im2,j,k),u0(ip1,j,k)>=0.) -      &
!                     u0(i,j,k)  *ip_weno(pin(ip2,j,k),pin(ip1,j,k),pin(i,j,k),pin(im1,j,k),pin(im2,j,k),pin(im3,j,k),u0(i,j,k)  >=0.))*dxi   &
!                   +(v0(i,jp1,k)*ip_weno(pin(i,jp3,k),pin(i,jp2,k),pin(i,jp1,k),pin(i,j,k),pin(i,jm1,k),pin(i,jm2,k),v0(i,jp1,k)>=0.) -      &
!                     v0(i,j,k)  *ip_weno(pin(i,jp2,k),pin(i,jp1,k),pin(i,j,k),pin(i,jm1,k),pin(i,jm2,k),pin(i,jm3,k),v0(i,j,k)  >=0.))*dyi   &
!                                    )

        ! Use original 5th order scheme in the horizontal directions, WENO is probably only really useful in the vertical direction at the inversion,
        ! although shallow cumulus convection could also profit.
!        pout(i,j,k) = pout(i,j,k) - ( &
!                    (u0(ip1,j,k)*ip_5th(pin(ip3,j,k),pin(ip2,j,k),pin(ip1,j,k),pin(i,j,k),pin(im1,j,k),pin(im2,j,k),sign(1.,u0(ip1,j,k))) -   &
!                     u0(i,j,k)  *ip_5th(pin(ip2,j,k),pin(ip1,j,k),pin(i,j,k),pin(im1,j,k),pin(im2,j,k),pin(im3,j,k),sign(1.,u0(i,j,k)) ))*dxi &
!                   +(v0(i,jp1,k)*ip_5th(pin(i,jp3,k),pin(i,jp2,k),pin(i,jp1,k),pin(i,j,k),pin(i,jm1,k),pin(i,jm2,k),sign(1.,v0(i,jp1,k))) -   &
!                     v0(i,j,k)  *ip_5th(pin(i,jp2,k),pin(i,jp1,k),pin(i,j,k),pin(i,jm1,k),pin(i,jm2,k),pin(i,jm3,k),sign(1.,v0(i,j,k)) ))*dyi &
!                                    )
        ! Do advection in the x-direction
        if (x_smooth(i,j,k) >= lambda_max) then ! not smooth, use weno
          pout(i,j,k) = pout(i,j,k) - ( &
                      (u0(ip1,j,k)*ip_weno(pin(ip3,j,k),pin(ip2,j,k),pin(ip1,j,k),pin(i,j,k),pin(im1,j,k),pin(im2,j,k),u0(ip1,j,k)>=0.) -      &
                       u0(i,j,k)  *ip_weno(pin(ip2,j,k),pin(ip1,j,k),pin(i,j,k),pin(im1,j,k),pin(im2,j,k),pin(im3,j,k),u0(i,j,k)  >=0.))*dxi   &
                                      )
        else ! smooth, use 5th order scheme
          pout(i,j,k) = pout(i,j,k) - ( &
                      (u0(ip1,j,k)*ip_5th(pin(ip3,j,k),pin(ip2,j,k),pin(ip1,j,k),pin(i,j,k),pin(im1,j,k),pin(im2,j,k),sign(1.,u0(ip1,j,k))) -   &
                       u0(i,j,k)  *ip_5th(pin(ip2,j,k),pin(ip1,j,k),pin(i,j,k),pin(im1,j,k),pin(im2,j,k),pin(im3,j,k),sign(1.,u0(i,j,k)) ))*dxi &
                                      )
        end if
        ! Do advection in the y-direction 
        if (y_smooth(i,j,k) >= lambda_max) then ! not smooth, use weno
          pout(i,j,k) = pout(i,j,k) - ( &
                      (v0(i,jp1,k)*ip_weno(pin(i,jp3,k),pin(i,jp2,k),pin(i,jp1,k),pin(i,j,k),pin(i,jm1,k),pin(i,jm2,k),v0(i,jp1,k)>=0.) -      &
                       v0(i,j,k)  *ip_weno(pin(i,jp2,k),pin(i,jp1,k),pin(i,j,k),pin(i,jm1,k),pin(i,jm2,k),pin(i,jm3,k),v0(i,j,k)  >=0.))*dyi   &
                                      )
        else ! smooth, use 5th order
          pout(i,j,k) = pout(i,j,k) - ( &
                      (v0(i,jp1,k)*ip_5th(pin(i,jp3,k),pin(i,jp2,k),pin(i,jp1,k),pin(i,j,k),pin(i,jm1,k),pin(i,jm2,k),sign(1.,v0(i,jp1,k))) -   &
                       v0(i,j,k)  *ip_5th(pin(i,jp2,k),pin(i,jp1,k),pin(i,j,k),pin(i,jm1,k),pin(i,jm2,k),pin(i,jm3,k),sign(1.,v0(i,j,k)) ))*dyi &
                                      )
        end if

        ! Do advection in the z-direction
        if (k==1) then
          ! Forward difference at k==1
          pout(i,j,k) = pout(i,j,k) - ( &
                       w0(i,j,k+1) * (pin(i,j,k+1)+pin(i,j,k)) &
                                      )/(2.*dzf(k))
        elseif (k==2 .or. k==kmax-1 .or. k==kmax) then
          ! 2nd order difference at k==2 and k==3
          pout(i,j,k) = pout(i,j,k) - ( &
                       w0(i,j,k+1) * (pin(i,j,k+1)+pin(i,j,k)) &
                      -w0(i,j,k)   * (pin(i,j,k-1)+pin(i,j,k)) &
                                      )/(2.*dzf(k))
        elseif (k==3) then
          ! 5th order at top, 2nd order below
          pout(i,j,k) = pout(i,j,k) - ( &
                       w0(i,j,k+1)*ip_5th(pin(i,j,k+3),pin(i,j,k+2),pin(i,j,k+1),pin(i,j,k),pin(i,j,k-1),pin(i,j,k-2),sign(1.,w0(i,j,k+1))) &
                      -w0(i,j,k)      * (pin(i,j,k-1)+pin(i,j,k))/2. &
                                      )/dzf(k)
        else ! No domain top or bottom level
          ! 5th order WENO for all layers between k=4 and k=kmax-2
          kp3=k+3;kp2=k+2;kp1=k+1;km1=k-1;km2=k-2;km3=k-3
          if (z_smooth(i,j,k) >= lambda_max) then ! field is locally not very smooth
            pout(i,j,k) = pout(i,j,k) - ( &
                        (w0(i,j,kp1)*ip_weno(pin(i,j,kp3),pin(i,j,kp2),pin(i,j,kp1),pin(i,j,k),pin(i,j,km1),pin(i,j,km2),w0(i,j,kp1)>=0.) - &
                         w0(i,j,k)  *ip_weno(pin(i,j,kp2),pin(i,j,kp1),pin(i,j,k),pin(i,j,km1),pin(i,j,km2),pin(i,j,km3),w0(i,j,k)  >=0.) ) &
                                        )/dzf(k)
          else ! field is smooth and no weno method is necessary
            pout(i,j,k) = pout(i,j,k) - ( &
                        (w0(i,j,kp1)*ip_5th(pin(i,j,kp3),pin(i,j,kp2),pin(i,j,kp1),pin(i,j,k),pin(i,j,km1),pin(i,j,km2),sign(1.,w0(i,j,kp1))) - &
                         w0(i,j,k)  *ip_5th(pin(i,j,kp2),pin(i,j,kp1),pin(i,j,k),pin(i,j,km1),pin(i,j,km2),pin(i,j,km3),sign(1.,w0(i,j,k  ))) ) &
                                        )/dzf(k)
          end if !Checking for non-smooth locations
        end if !Selection of bottom/top levels
      end do !Loop over k
    end do !Loop over i
  end do !Loop over j


end subroutine advecc_hybrid

!=======================================================================================
! Function that interpolates cell centered values of a variable to the appropriate cell face
! using a six point stencil. Fifth order accurate, because of the diffusive term included.
!=======================================================================================
elemental function ip_5th(vp2,vp1,v,vm1,vm2,vm3,sgn)
  implicit none
  real            :: ip_5th
  real,intent(in) :: vp2,vp1,v,vm1,vm2,vm3
  real,intent(in) :: sgn

  ip_5th = (37.*(v+vm1)-8.*(vp1+vm2)+(vp2+vm3)-sgn*(10.*(v-vm1)-5.*(vp1-vm2)+(vp2-vm3)))/60.

end function ip_5th

!=======================================================================================
! Subroutine that calculates lambda, which is a so-called smoothness parameter, defined by
! Blossey and Durran (2008). Could also be calculated using method of Hill and Pullin, although
! that is more computationally demanding.
!=======================================================================================

function smoothness(pin,dir)
  use modglobal, only : kmax,ih,i1,jh,j1,k1
  use modfields, only : u0,v0,w0
  implicit none
  real,intent(in),dimension(2-ih:i1+ih,2-jh:j1+jh,k1) :: pin
  integer,intent(in) :: dir
  real,dimension(2-ih:i1+ih,2-jh:j1+jh,k1) :: smoothness
  real,dimension(3) :: gam=0.
  real    :: eps_hybrid,phi_tilde
  integer :: i,j,k,n
  integer :: ip3,ip2,ip1,im1,im2,im3
  integer :: jp3,jp2,jp1,jm1,jm2,jm3
  integer :: kp3,kp2,kp1,km1,km2,km3

  smoothness(:,:,:)=0.
  !eps_hybrid indicates the minimum in jumps in the scalar that are worthy of attention
  !depends on specific scalar phi_tilde={1e-3 for qt, 1 for theta, 1e3 for N}
  !(from email correspondence with Peter Blossey)
  if (any(pin>=1.e5)) then ! probably number density
    phi_tilde = 1.e3
  elseif (any(pin>=1.)) then ! probably (potential) temperature
    phi_tilde = 1.
  else ! probably qt
    phi_tilde = 1.e-3
  end if
  
  ! Calculate 'small' value to keep from dividing by zero. Value is relative to a
  ! representative value of the variable being advected (=phi_tilde where phi \in {qt,thl,sv,u,v,w,tke})
  eps_hybrid = 1.e-8*phi_tilde**2

  select case (dir)
  case (1) ! x-direction
    do k=1,k1
      do j=2,j1
        do i=2,i1
          ip3=i+3;ip2=i+2;ip1=i+1;im1=i-1;im2=i-2;im3=i-3
          ! Calculate the smoothness for each stencil in an upwind configuration
          if (u0(i,j,k).ge.0.) then
            gam(:) = (pin(im1:ip1,j,k)-pin(im2:i,j,k))**2 + &
                     (pin(im2:i,j,k)-pin(im3:im1,j,k))**2
          else
            gam(:) = (pin(i:ip2,j,k)-pin(im1:ip1,j,k))**2 + &
                     (pin(im1:ip1,j,k)-pin(im2:i,j,k))**2
          end if
          smoothness(i,j,k) = maxval(gam)/(minval(gam)+eps_hybrid)
        end do
      end do
    end do
  case (2) ! y-direction
    do k=1,k1
      do j=2,j1
        jp3=j+3;jp2=j+2;jp1=j+1;jm1=j-1;jm2=j-2;jm3=j-3
        do i=2,i1
          ! Calculate the smoothness for each stencil in an upwind configuration
          if (u0(i,j,k).ge.0.) then
            gam(:) = (pin(i,jm1:jp1,k)-pin(i,jm2:j,k))**2 + &
                     (pin(i,jm2:j,k)-pin(i,jm3:jm1,k))**2
          else
            gam(:) = (pin(i,j:jp2,k)-pin(i,jm1:jp1,k))**2 + &
                     (pin(i,jm1:jp1,k)-pin(i,jm2:j,k))**2
          end if
          smoothness(i,j,k) = maxval(gam)/(minval(gam)+eps_hybrid)
        end do
      end do
    end do
  case (3) ! z-direction
    do k=4,kmax-2 ! Do not analyse bottom and top levels, because WENO cannot be used there anyway
      kp3=k+3;kp2=k+2;kp1=k+1;km1=k-1;km2=k-2;km3=k-3
      do j=2,j1
        do i=2,i1
          ! Calculate the smoothness for each stencil in an upwind configuration
          if (w0(i,j,k).ge.0.) then
            gam(:) = (pin(i,j,km1:kp1)-pin(i,j,km2:k))**2 + &
                     (pin(i,j,km2:k)-pin(i,j,km3:km1))**2
          else
            gam(:) = (pin(i,j,k:kp2)-pin(i,j,km1:kp1))**2 + &
                     (pin(i,j,km1:kp1)-pin(i,j,km2:k))**2
          end if
          smoothness(i,j,k) = maxval(gam)/(minval(gam)+eps_hybrid)
        end do
      end do
    end do
  case default
    stop 'ERROR: incorrect direction selected'
  end select
end function smoothness

end module advec_hybrid
