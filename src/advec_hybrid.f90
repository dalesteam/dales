!> \file advec_hybrid.f90
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
use modprecision, only: field_r
implicit none
contains

subroutine advecc_hybrid(pin,pout)
  use modglobal, only : ih,i1,jh,j1,kmax,k1 &
                       ,dxi,dyi,dzf,lambda_crit
  use modfields, only : u0,v0,w0,rhobf
  implicit none

  ! input and output variables
  real(field_r),dimension(2-ih:i1+ih,2-jh:j1+jh,k1),intent(in)   :: pin  !< Input: the cell centered field (qt,thetal,sv etc)
  real(field_r),dimension(2-ih:i1+ih,2-jh:j1+jh,k1),intent(inout):: pout !< Output: the tendency for the input field (qtp,thetalp,svp etc)

  real(field_r),dimension(2-ih:i1+ih,2-jh:j1+jh,k1) :: rhopin  !< 3D density profile * input

  ! local variables
  logical,dimension(2:i1+1,2:j1+1,k1) :: lsmx,lsmy,lsmz
  real,dimension(2:i1+1,2:j1+1,k1) :: pfacex,pfacey,pfacez ! face values, defined at the same interfaces as u0,v0 and w0 respectively
  integer :: i,j,k
  integer :: kp2,km3
  integer :: jp2,jm3
  integer :: ip2,im3

  ! Anelastic approx.
  do k=1,k1
    do j=2-jh,j1+jh
      do i=2-ih,i1+ih
        rhopin(i,j,k)=rhobf(k)*pin(i,j,k)
      end do
    end do
  end do

  ! Initialize face values
  pfacex=0.;pfacey=0.;pfacez=0.

  ! calculate the smoothness indicator in all 3 directions and immediately check if it exceeds the critical value
  ! note that the smoothness indicator is defined at cell [i]faces[/i].
  lsmx = smoothness(pin,1)<lambda_crit
  lsmy = smoothness(pin,2)<lambda_crit
  lsmz = smoothness(pin,3)<lambda_crit

  ! Do interpolation and save cell face values in a matrix. Avoid calculating cell face values twice.
  ! In the old configuration, it was possible that a value at a face was calculated twice with different methods!
  ! This was because the smoothness metric was supposed to be defined at cell centers, while it is actually at a face.

  ! Levels outside loop over k, to avoid if then else inside do-loop
  ! Lowest level (k=1) is at the surface, z=0. No flux is calculated with it
  ! Simple difference to determine second and third half levels and highest two half levels
  ! The smoothness is not checked at these levels, because 5th order (WENO) advection cannot be used anyway.
  pfacez(2:i1,2:j1,2:3)     = ( rhopin(2:i1,2:j1,2:3)+rhopin(2:i1,2:j1,1:2) )/2    !note that pin and pfacez do not have the same dimensions!
  pfacez(2:i1,2:j1,kmax:k1) = ( rhopin(2:i1,2:j1,kmax:k1)+rhopin(2:i1,2:j1,kmax-1:kmax) )/2

  ! Start looping over the remaining points
  do j=2,j1+1
    jp2=j+2;jm3=j-3
    do i=2,i1+1
      ip2=i+2;im3=i-3
      ! Loop over first two height levels to do horizontal interpolation
      do k=1,3
        pfacex(i,j,k) = ip_hybrid(pin(im3:ip2,j,k),u0(i,j,k)>=0.,lsmx(i,j,k))
        pfacey(i,j,k) = ip_hybrid(pin(i,jm3:jp2,k),v0(i,j,k)>=0.,lsmy(i,j,k))
      end do
      ! Loop over last two height levels to do horizontal interpolation
      do k=kmax,k1
        pfacex(i,j,k) = ip_hybrid(pin(im3:ip2,j,k),u0(i,j,k)>=0.,lsmx(i,j,k))
        pfacey(i,j,k) = ip_hybrid(pin(i,jm3:jp2,k),v0(i,j,k)>=0.,lsmy(i,j,k))
      end do
      ! Loop over rest of levels, for horizontal and vertical faces
      do k=4,kmax-1
        kp2=k+2;km3=k-3
        pfacex(i,j,k) = ip_hybrid(pin(im3:ip2,j,k),u0(i,j,k)>=0.,lsmx(i,j,k))
        pfacey(i,j,k) = ip_hybrid(pin(i,jm3:jp2,k),v0(i,j,k)>=0.,lsmy(i,j,k))
        pfacez(i,j,k) = ip_hybrid(rhopin(i,j,km3:kp2),w0(i,j,k)>=0.,lsmz(i,j,k))
      end do !Loop over k
    end do !Loop over i
  end do !Loop over j

  ! Calculate actual tendencies by multiplying matrices, accept in the vertical, since dzf(k)
  ! does not have the appropriate dimensions.
  do k=1,kmax
    pout(2:i1,2:j1,k) = pout(2:i1,2:j1,k) - ( &
                              (u0(3:i1+1,2:j1,k)*pfacex(3:i1+1,2:j1,k) -    &
                               u0(2:i1,2:j1,k)*pfacex(2:i1,2:j1,k) )*dxi    &
                             +(v0(2:i1,3:j1+1,k)*pfacey(2:i1,3:j1+1,k) -    &
                               v0(2:i1,2:j1,k)*pfacey(2:i1,2:j1,k) )*dyi    &
                             +(1./rhobf(k))*(w0(2:i1,2:j1,k+1)*pfacez(2:i1,2:j1,k+1) -    &
                               w0(2:i1,2:j1,k)*pfacez(2:i1,2:j1,k) )/dzf(k) &
                                            )
  end do

end subroutine advecc_hybrid

!=======================================================================================
! Function that interpolates cell centered values of a variable to the appropriate cell face
! using a six point stencil. Fifth order accurate, because of the diffusive term included.
!=======================================================================================
!elemental function ip_5th(vp2,vp1,v,vm1,vm2,vm3,sgn)
!  implicit none
!  real            :: ip_5th
!  real,intent(in) :: vp2,vp1,v,vm1,vm2,vm3
!  real,intent(in) :: sgn
!
!  ip_5th = (37.*(v+vm1)-8.*(vp1+vm2)+(vp2+vm3)-sgn*(10.*(v-vm1)-5.*(vp1-vm2)+(vp2-vm3)))/60.
!
!end function ip_5th

!=======================================================================================
! This function checks whether fifth order advection should be used, or WENO method, based
! on the smoothness of different stencils
!=======================================================================================
function ip_hybrid(vin,lpos,lsmooth)
  implicit none
  real :: ip_hybrid
  real(field_r),intent(in),dimension(-3:2) :: vin   ! Subset of a variable
  logical,intent(in) :: lpos               ! Positive of negative velocity
  logical,intent(in) :: lsmooth            ! Locally smooth or non-smooth

  !local variables
  real,parameter    :: c1=13./12.,c2=1./4. ! Multiplication constants
  real,dimension(3) :: wgtOpt=(/.1,.6,.3/) ! Optimal weights (see eg Hill and Pullin)
  real,dimension(3) :: beta,             & ! Smoothness measure
                       wgt,              & ! Weighting factor
                       varFace             ! Interpolated value at cell face
  real              :: wgtfac              ! Normalization factor for the weights
  integer,parameter :: pweno=1             ! Exponent used in WENO scheme
                                           ! Values 1,2 or 3 should work
  real,parameter    :: epsWeno=1e-12       ! Small value set to keep from dividing by zero
  integer           :: sgn                 ! 1 if velocity is positive, -1 if velocity is negative

  sgn=-1

  if (lsmooth) then ! field around this location is smooth -> use regular 5th order (upwind)
    if (lpos) sgn=1 ! set sgn, to account for different wind directions
    ip_hybrid = (37.*(vin(0)+vin(-1))-8.*(vin(1)+vin(-2))+(vin(2)+vin(-3)) &
                  -sgn*(10.*(vin(0)-vin(-1))-5.*(vin(1)-vin(-2))+(vin(2)-vin(-3))))/60.
  else ! field around this location is non-smooth -> use weno
    if (lpos) then !Positive velocity at cell face
      !compute smoothness indicators for each of the stencils
      beta(1) = c1*(vin(-3)-2*vin(-2)+vin(-1))**2 + c2*(vin(-3)-4*vin(-2)+3*vin(-1))**2
      beta(2) = c1*(vin(-2)-2*vin(-1)+vin(0) )**2 + c2*(vin(-2)-vin(0))**2
      beta(3) = c1*(vin(-1)-2*vin(0) +vin(1) )**2 + c2*(3*vin(-1)-4*vin(0)+vin(1))**2

      !interpolated values of the variable at the cell faces using each of the stencils
      varFace(1) = (2*vin(-3)- 7*vin(-2)+ 11*vin(-1))/6
      varFace(2) = (- vin(-2)+ 5*vin(-1)+  2*vin(0) )/6
      varFace(3) = (2*vin(-1)+ 5*vin(0) -    vin(1) )/6
    else !Negative velocity at cell face
      !compute smoothness indicators for each of the stencils
      !the following is found by mirroring the equations for positive velocity
      beta(1) = c1*(vin(0) -2*vin(1) +vin(2))**2 + c2*(3*vin(0)-4*vin(1)+vin(2))**2
      beta(2) = c1*(vin(-1)-2*vin(0) +vin(1))**2 + c2*(vin(-1)-vin(1))**2
      beta(3) = c1*(vin(-2)-2*vin(-1)+vin(0))**2 + c2*(vin(-2)-4*vin(-1)+3*vin(0))**2

      !interpolated values of the variable at the cell faces using each of the stencils
      varFace(1) = (11*vin(0) -7*vin(1) + 2*vin(2))/6
      varFace(2) = ( 2*vin(-1)+5*vin(0) -   vin(1))/6
      varFace(3) = ( - vin(-2)+5*vin(-1)+ 2*vin(0))/6
    end if

    !compute weights
    wgt = wgtOpt*(epsWeno+beta)**(-pweno)
    wgtfac = sum(wgt)**(-1)

    ! compute interpolated value
    ip_hybrid = sum(wgt(:)*varFace(:))*wgtfac
  end if

end function ip_hybrid

!=======================================================================================
! Subroutine that calculates lambda, which is a so-called smoothness parameter, defined by
! Blossey and Durran (2008). Could also be calculated using method of Hill and Pullin, although
! that is more computationally demanding.
!=======================================================================================

function smoothness(pin,dir)
  use modglobal, only : kmax,ih,i1,jh,j1,k1
  use modfields, only : u0,v0,w0
  implicit none
  real(field_r),intent(in),dimension(2-ih:i1+ih,2-jh:j1+jh,k1) :: pin
  integer,intent(in) :: dir
  real,dimension(2:i1+1,2:j1+1,k1) :: smoothness
  real,dimension(3) :: gam=0.
  real    :: eps_hybrid,phi_tilde
  integer :: i,j,k
  integer :: ip2,ip1,im1,im2,im3
  integer :: jp2,jp1,jm1,jm2,jm3
  integer :: kp2,kp1,km1,km2,km3

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
      do j=2,j1+1
        do i=2,i1+1
          ip2=i+2;ip1=i+1;im1=i-1;im2=i-2;im3=i-3
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
      do j=2,j1+1
        jp2=j+2;jp1=j+1;jm1=j-1;jm2=j-2;jm3=j-3
        do i=2,i1+1
          ! Calculate the smoothness for each stencil in an upwind configuration
          if (v0(i,j,k).ge.0.) then
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
      kp2=k+2;kp1=k+1;km1=k-1;km2=k-2;km3=k-3
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
