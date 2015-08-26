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
implicit none

contains

subroutine advecc_hybrid(pin,pout,phi_tilde)
  use modglobal, only : ih,i1,i2,jh,j1,j2,kmax,k1,dxi,dyi,rdt
  use modfields, only : u0,v0,w0,rhobf,rhodzi
  implicit none

  ! input and output variables
  real,dimension(2-ih:,2-jh:,:),intent(in)   :: pin  !< Input: the cell centered field (qt,thetal,sv etc)
  real,dimension(2-ih:,2-jh:,:),intent(inout):: pout !< Output: the tendency for the input field (qtp,thetalp,svp etc)
  real,intent(in) :: phi_tilde

  ! local variables
  real,dimension(2-ih:i1+ih,2-jh:j1+jh,k1) :: rhopin  !< 3D density profile * input
  real,dimension(2:i2,2:j2,k1) :: pfacex,pfacey,pfacez
  real    :: eps_hybrid
  integer :: i,j,k
  integer :: ip,jp,kp

  ! Anelastic approx.
  do k=1,k1
    rhopin(:,:,k)=rhobf(k)*pin(:,:,k)
  end do

  ! Calculate 'small' value to keep from dividing by zero. Value is relative to a
  ! representative value of the variable being advected (=phi_tilde where phi \in {qt,thl,sv,tke})
  eps_hybrid = 1.e-8*phi_tilde**2

  ! calculate the smoothness indicator in all 3 directions and immediately check if it exceeds the critical value
  ! note that the smoothness indicator is defined at cell [i]faces[/i].
  call interp_x(   pin,pfacex,eps_hybrid)
  call interp_y(   pin,pfacey,eps_hybrid)
  call interp_z(rhopin,pfacez,eps_hybrid)

  ! Levels outside loop over k, to avoid if then else inside do-loop
  ! Lowest level (k=1) is at the surface, z=0. No flux is calculated with it
  ! Simple difference to determine second and third half levels and highest two half levels
  ! The smoothness is not checked at these levels, because 5th order (WENO) advection cannot be used anyway.
  pfacez(2:i1,2:j1,2:3)     = .5*( rhopin(2:i1,2:j1,2:3)    +rhopin(2:i1,2:j1,1:2)         )    !note that pin and pfacez do not have the same dimensions!
  pfacez(2:i1,2:j1,kmax:k1) = .5*( rhopin(2:i1,2:j1,kmax:k1)+rhopin(2:i1,2:j1,kmax-1:kmax) )
  
  ! Calculate actual tendencies by multiplying matrices, accept in the vertical, since rhodzi(k)
  ! does not have the appropriate dimensions.
  do k=1,kmax
    kp=k+1
    do j=2,j1
      jp=j+1
      do i=2,i1
        ip=i+1
        pout(i,j,k) = pout(i,j,k) - ( &
                            (u0(ip,j,k)*pfacex(ip,j,k) -        &
                             u0(i,j,k)*pfacex(i,j,k))*dxi +     &
                            (v0(i,jp,k)*pfacey(i,jp,k) -        &
                             v0(i,j,k)*pfacey(i,j,k))*dyi +     &
                            (w0(i,j,kp)*pfacez(i,j,kp) -        &
                             w0(i,j,k)*pfacez(i,j,k))*rhodzi(k) &
                                    )
      end do
    end do
  end do

end subroutine advecc_hybrid

!=======================================================================================
! Subroutine that determines whether a stencil is smooth or not. If it is, then
! 5th order advection is applied. If not, then WENO interpolation is used. The
! smoothness function of Blossey and Durran (2008) is used. Could also be calculated using 
! method of Hill and Pullin, although that is more computationally demanding.
!=======================================================================================
subroutine interp_x(pin,pout,eps_hybrid)
  use modglobal, only : kmax,ih,i1,i2,jh,j1,j2,k1,lambda_crit,rdt
  use modfields, only : u0,smoothx
  implicit none
  real,intent(in),dimension(2-ih:,2-jh:,:) :: pin
  real,intent(out),dimension(2:,2:,:) :: pout
  real,intent(in) :: eps_hybrid
  real,dimension(3) :: gam
  real    :: lambda
  integer :: i,j,k
  integer :: ip2,ip1,im1,im2,im3

  do k=1,k1
    do j=2,j2
      do i=2,i2
      ip2=i+2;ip1=i+1;im1=i-1;im2=i-2;im3=i-3
        ! Calculate the smoothness for each stencil in an upwind configuration
        if (u0(i,j,k)>=0.) then
          gam(:) = (pin(im1:ip1,j,k)-pin(im2:i,j,k))**2 + &
                   (pin(im2:i,j,k)-pin(im3:im1,j,k))**2
        else
          gam(:) = (pin(i:ip2,j,k)-pin(im1:ip1,j,k))**2 + &
                   (pin(im1:ip1,j,k)-pin(im2:i,j,k))**2
        end if
        lambda = maxval(gam)/(minval(gam)+eps_hybrid)
        ! Select which interpolation schemes should be used based on smoothness metric lambda
        if (lambda > lambda_crit) then
!          if (abs(eps_hybrid-1.e-14)<1e-18) smoothx(k) = smoothx(k) + rdt
          pout(i,j,k) = ip_weno(pin(im3:ip2,j,k),u0(i,j,k)>=0.)
        else
          pout(i,j,k) = ip_5th (pin(im3:ip2,j,k),sign(1.,u0(i,j,k)))
        end if
      end do
    end do
  end do
end subroutine interp_x

subroutine interp_y(pin,pout,eps_hybrid)
  use modglobal, only : kmax,ih,i1,i2,jh,j1,j2,k1,lambda_crit,rdt
  use modfields, only : v0,smoothy
  implicit none
  real,intent(in),dimension(2-ih:,2-jh:,:) :: pin
  real,intent(out),dimension(2:,2:,:) :: pout
  real,intent(in) :: eps_hybrid
  real,dimension(3) :: gam
  real    :: lambda
  integer :: i,j,k
  integer :: jp2,jp1,jm1,jm2,jm3
 
  do k=1,k1
    do j=2,j2
      jp2=j+2;jp1=j+1;jm1=j-1;jm2=j-2;jm3=j-3
      do i=2,i2
        ! Calculate the smoothness for each stencil in an upwind configuration
        if (v0(i,j,k)>=0.) then
          gam(:) = (pin(i,jm1:jp1,k)-pin(i,jm2:j,k))**2 + &
                   (pin(i,jm2:j,k)-pin(i,jm3:jm1,k))**2
        else
          gam(:) = (pin(i,j:jp2,k)-pin(i,jm1:jp1,k))**2 + &
                   (pin(i,jm1:jp1,k)-pin(i,jm2:j,k))**2
        end if
        lambda = maxval(gam)/(minval(gam)+eps_hybrid)
        ! Select which interpolation schemes should be used based on smoothness metric lambda
        if (lambda > lambda_crit) then
!          if (abs(eps_hybrid-1.e-14)<1e-18) smoothy(k) = smoothy(k) + rdt
          pout(i,j,k) = ip_weno(pin(i,jm3:jp2,k),v0(i,j,k)>=0.)
        else
          pout(i,j,k) = ip_5th (pin(i,jm3:jp2,k),sign(1.,v0(i,j,k)))
        end if
      end do
    end do
  end do
end subroutine interp_y

subroutine interp_z(pin,pout,eps_hybrid)
  use modglobal, only : kmax,ih,i1,i2,jh,j1,j2,k1,lambda_crit,rdt
  use modfields, only : w0,smoothz
  implicit none
  real,intent(in),dimension(2-ih:,2-jh:,:) :: pin
  real,intent(out),dimension(2:,2:,:) :: pout
  real,intent(in) :: eps_hybrid
  real,dimension(3) :: gam
  real    :: lambda
  integer :: i,j,k
  integer :: kp2,kp1,km1,km2,km3

  do k=4,kmax-1 ! Do not analyse bottom and top levels, because WENO cannot be used there anyway
    kp2=k+2;kp1=k+1;km1=k-1;km2=k-2;km3=k-3
    do j=2,j1
      do i=2,i1
        ! Calculate the smoothness for each stencil in an upwind configuration
        if (w0(i,j,k)>=0.) then
          gam(:) = (pin(i,j,km1:kp1)-pin(i,j,km2:k))**2 + &
                   (pin(i,j,km2:k)-pin(i,j,km3:km1))**2
        else
          gam(:) = (pin(i,j,k:kp2)-pin(i,j,km1:kp1))**2 + &
                   (pin(i,j,km1:kp1)-pin(i,j,km2:k))**2
        end if
        lambda = maxval(gam)/(minval(gam)+eps_hybrid)
        ! Select which interpolation schemes should be used based on smoothness metric lambda
        if (lambda > lambda_crit) then
!          if (abs(eps_hybrid-1.e-14)<1e-18) smoothz(k) = smoothz(k) + rdt
          pout(i,j,k) = ip_weno(pin(i,j,km3:kp2),w0(i,j,k)>=0.)
        else
          pout(i,j,k) = ip_5th (pin(i,j,km3:kp2),sign(1.,w0(i,j,k)))
        end if
      end do
    end do
  end do
end subroutine interp_z

!=======================================================================================
! Function that interpolates cell centered values of a variable to the appropriate cell face
! using a six point stencil. Fifth order accurate, because of the diffusive term included.
!=======================================================================================
real function ip_5th(vp,sgn) result(face)
  implicit none
  real,dimension(:),intent(in) :: vp
  real,intent(in) :: sgn

  face = (37*(vp(4)+vp(3))-8*(vp(5)+vp(2))+(vp(6)+vp(1))-sgn*(10*(vp(4)-vp(3))-5*(vp(5)-vp(2))+(vp(6)-vp(1))))/60

end function ip_5th

!=======================================================================================
! The value at the cell face is determined by doing WENO interpolation (Hill and
! Pullin 2004)
! 
!=======================================================================================
real function ip_weno(vp,lpos) result(face)
  implicit none
  logical,intent(in) :: lpos               ! Positive or negative velocity
  real,dimension(:),intent(in) :: vp
  !local variables
  real,parameter    :: c1=13./12.,c2=1./4. ! Multiplication constants
  real,parameter,dimension(3) :: wgtOpt=(/.1,.6,.3/) ! Optimal weights (see eg Hill and Pullin)
  real,dimension(3) :: beta,             & ! Smoothness measure
                       wgt,              & ! Weighting factor
                       varFace             ! Interpolated value at cell face
  real              :: wgtfac              ! Normalization factor for the weights
  integer,parameter :: pweno=1             ! Exponent used in WENO scheme
                                           ! Values 1,2 or 3 should work
  real,parameter    :: epsWeno=1.e-12      ! Small value set to keep from dividing by zero
 
  ! field around this location is non-smooth -> use weno
  if (lpos) then !Positive velocity at cell face
    !compute smoothness indicators for each of the stencils
    beta(1) = c1*(vp(1)-2*vp(2)+vp(3))**2 + c2*(  vp(1)-4*vp(2)+3*vp(3))**2
    beta(2) = c1*(vp(2)-2*vp(3)+vp(4))**2 + c2*(  vp(2)        -  vp(4))**2
    beta(3) = c1*(vp(3)-2*vp(4)+vp(5))**2 + c2*(3*vp(3)-4*vp(4)+  vp(5))**2

    !interpolated values of the variable at the cell faces using each of the stencils
    varFace(1) = (2*vp(1)-7*vp(2)+11*vp(3))/6
    varFace(2) = (- vp(2)+5*vp(3)+ 2*vp(4))/6
    varFace(3) = (2*vp(3)+5*vp(4)-   vp(5))/6
  else !Negative velocity at cell face
    !compute smoothness indicators for each of the stencils
    !the following is found by mirroring the equations for positive velocity
    beta(1) = c1*(vp(4)-2*vp(5)+vp(6))**2 + c2*(3*vp(4)-4*vp(5)+  vp(6))**2
    beta(2) = c1*(vp(3)-2*vp(4)+vp(5))**2 + c2*(  vp(3)        -  vp(5))**2
    beta(3) = c1*(vp(2)-2*vp(3)+vp(4))**2 + c2*(  vp(2)-4*vp(3)+3*vp(4))**2

    !interpolated values of the variable at the cell faces using each of the stencils
    varFace(1) = (11*vp(4)-7*vp(5)+2*vp(6))/6
    varFace(2) = ( 2*vp(3)+5*vp(4)-  vp(5))/6
    varFace(3) = ( - vp(2)+5*vp(3)+2*vp(4))/6
  end if

  !compute weights
  wgt = wgtOpt*(epsWeno+beta)**(-pweno)
  wgtfac = sum(wgt)**(-1)

  ! compute interpolated value 
  face = sum(wgt(:)*varFace(:))*wgtfac

end function ip_weno

end module advec_hybrid
