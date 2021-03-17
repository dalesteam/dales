!> \file advec_hybrid.f90
!  Does advection with the 5th order advection scheme that was already present in DALES, except
!  around locations were discontinuities arise. There, the 5th order WENO scheme is used.
!  Discontinuities should be preserved, while the damping of high wavenumber components of the flow,
!  common in pure WENO solutions.
!  The scheme is more or less equal to the scheme proposed by Hill and Pullin (2004)
!  [https://doi.org/10.1016/j.jcp.2003.07.032] but using
!  the fifth order advection scheme instead of their tuned one.
! 
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
module advec_hybrid_f
use modprecision, only : field_r
implicit none
contains

subroutine advecc_hybrid_f(pin, pout, phi_tilde_in)
  use modglobal, only : ih,i1,jh,j1,kmax,k1,dxi,dyi,dzf,lambda_crit
  use modfields, only : u0,v0,w0,rhobf

  implicit none
  
  real(field_r),dimension(2-ih:i1+ih,2-jh:j1+jh,k1),intent(in)   :: pin  !< Input: the cell centered field (qt,thetal,sv etc)
  real(field_r),dimension(2-ih:i1+ih,2-jh:j1+jh,k1),intent(inout):: pout !< Output: the tendency for the input field (qtp,thetalp,svp etc)
  real,optional,intent(in) :: phi_tilde_in   !< Order of magnitude of the field, used in the smoothness criterion. Optional. 

  real :: phi_tilde
  logical :: lsmx,lsmy,lsmz                   ! smoothness flags
  real,dimension(2:i1+1,2:j1+1,k1) :: pfacex,pfacey,pfacez  ! face values, defined at the same interfaces as u0,v0 and w0 respectively
  real,dimension(3)    :: gam                 ! used for smoothness test
  integer              :: i,j,k
  real                 :: eps_hybrid
  real,dimension(-3:2) :: vin                 ! Subset of the field to be advected
  real                 :: sgn                 ! 1 if velocity is positive, -1 if velocity is negative
  real,parameter       :: c1=13./12.,c2=1./4. ! Multiplication constants
  integer,parameter    :: pweno=1             ! Exponent used in WENO scheme
                                              ! Values 1,2 or 3 should work
  real,parameter       :: epsWeno=1e-12       ! Small value set to keep from dividing by zero
  real,dimension(3)    :: wgtOpt=(/.1,.6,.3/) ! Optimal weights (see eg Hill and Pullin)
  real,dimension(3)    :: beta,          &    ! Smoothness measure
                          wgt,         &      ! Weighting factor
                          varFace             ! Interpolated value at cell face
  real                 :: wgtfac              ! Normalization factor for the weights


  ! phi_tilde is some kind of order-of-magnitude, used to calculate eps_hybrid
  ! it's unclear to me if it's necessary for this scheme or just used to avoid dividing by 0
  ! but for now it's here, to give the same results as the original routine
  ! If phi_tilde is passed as a function argument, use it, otherwise determine heuristically as in the original
  ! e12 may be on either side of 1,
  ! other fields like T, qt should always end up with the same phi_tilde.
  if(.not. present(phi_tilde_in)) then
     if (any(pin>=1.e5)) then   ! probably number density
        phi_tilde = 1.e3
     elseif (any(pin>=1.)) then ! probably (potential) temperature
        phi_tilde = 1.
     else                         ! probably qt
        phi_tilde = 1.e-3
     end if
  else
     phi_tilde = phi_tilde_in
  end if
  eps_hybrid = 1.e-8*phi_tilde**2
   
  do k=1,k1
     do j=2,j1+1
        do i=2,i1+1
           ! determine smoothness lsmx
           if (u0(i,j,k).ge.0.) then
              gam(:) = (pin(i-1:i+1,j,k)-pin(i-2:i,j,k))**2 + &
                   (pin(i-2:i,j,k)-pin(i-3:i-1,j,k))**2
           else
              gam(:) = (pin(i:i+2,j,k)-pin(i-1:i+1,j,k))**2 + &
                   (pin(i-1:i+1,j,k)-pin(i-2:i,j,k))**2
           end if
           !lsmx = maxval(gam)/(minval(gam)+eps_hybrid) < lambda_crit
           lsmx = maxval(gam) < lambda_crit * (minval(gam)+eps_hybrid)
           
           ! determine smoothness lsmy
           if (v0(i,j,k).ge.0.) then
              gam(:) = (pin(i,j-1:j+1,k)-pin(i,j-2:j,k))**2 + &
                   (pin(i,j-2:j,k)-pin(i,j-3:j-1,k))**2
           else
              gam(:) = (pin(i,j:j+2,k)-pin(i,j-1:j+1,k))**2 + &
                   (pin(i,j-1:j+1,k)-pin(i,j-2:j,k))**2
           end if
           !lsmy = maxval(gam)/(minval(gam)+eps_hybrid) < lambda_crit
           lsmy = maxval(gam) < lambda_crit * (minval(gam)+eps_hybrid)


           ! advection in x
           vin(-3:2) = pin(i-3:i+2,j,k)
           if (lsmx) then ! field around this location is smooth -> use regular 5th order (upwind)
              sgn = sign(1.0_field_r,u0(i,j,k))  ! set sgn, to account for different wind directions
              pfacex(i,j,k) = (37.*(vin(0)+vin(-1))-8.*(vin(1)+vin(-2))+(vin(2)+vin(-3)) &
                   -sgn*(10.*(vin(0)-vin(-1))-5.*(vin(1)-vin(-2))+(vin(2)-vin(-3))))/60.
           else ! field around this location is non-smooth -> use weno
              if (u0(i,j,k) >= 0) then !Positive velocity at cell face
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
              pfacex(i,j,k) = sum(wgt(:)*varFace(:))*wgtfac
              
           end if

            ! advection in y
           vin(-3:2) = pin(i,j-3:j+2,k)
           if (lsmy) then ! field around this location is smooth -> use regular 5th order (upwind)
              sgn = sign(1.0_field_r,v0(i,j,k))  ! set sgn, to account for different wind directions
              pfacey(i,j,k) = (37.*(vin(0)+vin(-1))-8.*(vin(1)+vin(-2))+(vin(2)+vin(-3)) &
                   -sgn*(10.*(vin(0)-vin(-1))-5.*(vin(1)-vin(-2))+(vin(2)-vin(-3))))/60.
           else ! field around this location is non-smooth -> use weno
              if (v0(i,j,k) >= 0) then !Positive velocity at cell face
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
              pfacey(i,j,k) = sum(wgt(:)*varFace(:))*wgtfac
           end if
           
           ! advection in z
           if (k < 4 .or. k >= kmax) then
              ! special treatment of top and bottom layers
              if (k == 1) then
                 pfacez(i,j,k) = 0
              else
                 pfacez(i,j,k) =  ( rhobf(k)*pin(i,j,k) + rhobf(k-1)*pin(i,j,k-1) ) * .5
              end if
           else
              
              lsmz = .true.
              if (k < kmax-1) then   ! the original scheme considers k=kmax-1 fully smooth
                 ! determine smoothness lsmz
                 if (w0(i,j,k).ge.0.) then
                    gam(:) = (pin(i,j,k-1:k+1)-pin(i,j,k-2:k))**2 + &
                         (pin(i,j,k-2:k)-pin(i,j,k-3:k-1))**2
                 else
                    gam(:) = (pin(i,j,k:k+2)-pin(i,j,k-1:k+1))**2 + &
                         (pin(i,j,k-1:k+1)-pin(i,j,k-2:k))**2
                 end if
                 ! lsmz = maxval(gam)/(minval(gam)+eps_hybrid) < lambda_crit
                 lsmz = maxval(gam) < lambda_crit * (minval(gam)+eps_hybrid)
              end if
              
              vin(-3:2) = pin(i,j,k-3:k+2) * rhobf(k-3:k+2)
              if (lsmz) then ! field around this location is smooth -> use regular 5th order (upwind)
                 sgn = sign(1.0_field_r,w0(i,j,k))  ! set sgn, to account for different wind directions
                 pfacez(i,j,k) = (37.*(vin(0)+vin(-1))-8.*(vin(1)+vin(-2))+(vin(2)+vin(-3)) &
                      -sgn*(10.*(vin(0)-vin(-1))-5.*(vin(1)-vin(-2))+(vin(2)-vin(-3))))/60.
              else ! field around this location is non-smooth -> use weno
                 if (w0(i,j,k) >= 0) then !Positive velocity at cell face
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
                 pfacez(i,j,k) = sum(wgt(:)*varFace(:))*wgtfac
              end if
           end if
           
           !kp2=k+2;km3=k-3
           !pfacex(i,j,k) = ip_hybrid(pin(im3:ip2,j,k),u0(i,j,k)>=0.,lsmx(i,j,k))
           !pfacey(i,j,k) = ip_hybrid(pin(i,jm3:jp2,k),v0(i,j,k)>=0.,lsmy(i,j,k))
           !pfacez(i,j,k) = ip_hybrid(rhopin(i,j,km3:kp2),w0(i,j,k)>=0.,lsmz(i,j,k))
        end do
     end do
  end do !Loop over k

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
end subroutine advecc_hybrid_f


end module advec_hybrid_f
