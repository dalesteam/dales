!> \file advec_weno.f90
!  Does advection with a 5th order weighted essentially non-oscillatory advection scheme
!  Should be able to remove spurious overshoots (mainly at the inversion at this point) present
!  in other schemes (also in 5th order scheme, not in Kappa, but that scheme is very diffusive).
!  Scheme is basically the same as in Jiang and Shu (1996)
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
module advec_weno
implicit none
contains

subroutine advecc_weno(pin,pout)

  use modglobal, only : ih,i1,jh,j1,kmax,k1,leq &
                        ,dx,dxi,dy,dyi,dzi,dzf
  use modfields, only : u0,v0,w0

  implicit none

  real,dimension(2-ih:i1+ih,2-jh:j1+jh,k1),intent(in) :: pin  !< Input: the cell centered field (qt,thetal,sv etc)
  real,dimension(2-ih:i1+ih,2-jh:j1+jh,k1),intent(out):: pout !< Output: the tendency for the input field (qtp,thetalp,svp etc)

  ! local variables
  integer :: i,j,k
  integer :: kp3,kp2,kp1,km1,km2,km3
  integer :: jp3,jp2,jp1,jm1,jm2,jm3
  integer :: ip3,ip2,ip1,im1,im2,im3

  do j=2,j1
    jp3=j+3;jp2=j+2;jp1=j+1;jm1=j-1;jm2=j-2;jm3=j-3
    do i=2,i1
      ip3=i+3;ip2=i+2;ip1=i+1;im1=i-1;im2=i-2;im3=i-3
      do k=1,kmax
        ! Use 5th order WENO in the horizontal directions
        pout(i,j,k) = pout(i,j,k) - ( &
                    (u0(ip1,j,k)*ip_weno(pin(ip3,j,k),pin(ip2,j,k),pin(ip1,j,k),pin(i,j,k),pin(im1,j,k),pin(im2,j,k),u0(ip1,j,k)>=0.) -       &
                     u0(i,j,k)  *ip_weno(pin(ip2,j,k),pin(ip1,j,k),pin(i,j,k),pin(im1,j,k),pin(im2,j,k),pin(im3,j,k),u0(i,j,k)  >=0.))*dxi    &
                   +(v0(i,jp1,k)*ip_weno(pin(i,jp3,k),pin(i,jp2,k),pin(i,jp1,k),pin(i,j,k),pin(i,jm1,k),pin(i,jm2,k),v0(i,jp1,k)>=0.) -       &
                     v0(i,j,k)  *ip_weno(pin(i,jp2,k),pin(i,jp1,k),pin(i,j,k),pin(i,jm1,k),pin(i,jm2,k),pin(i,jm3,k),v0(i,j,k)  >=0.))*dyi    &
                                    )
        if (k==1) then
          ! Forward difference at k==1
          pout(i,j,k) = pout(i,j,k) - ( &
                          w0(i,j,k+1) * (pin(i,j,k+1) + pin(i,j,k)) &
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
                      w0(i,j,k+1)/60.&
                      *(37.*(pin(i,j,k+1)+pin(i,j,k))-8.*(pin(i,j,k+2)+pin(i,j,k-1))+(pin(i,j,k+3)+pin(i,j,k-2)))&
                      -sign(1.,w0(i,j,k+1))*w0(i,j,k+1)/60.&
                      *(10.*(pin(i,j,k+1)-pin(i,j,k))-5.*(pin(i,j,k+2)-pin(i,j,k-1))+(pin(i,j,k+3)-pin(i,j,k-2)))&
                      -w0(i,j,k)      * (pin(i,j,k-1)+pin(i,j,k))/2. &
                                      )/dzf(k)
        else
          ! 5th order WENO for all layers between k=4 and k=kmax-2
          kp3=k+3;kp2=k+2;kp1=k+1;km1=k-1;km2=k-2;km3=k-3

          pout(i,j,k) = pout(i,j,k) - ( &
                      (w0(i,j,kp1)*ip_weno(pin(i,j,kp3),pin(i,j,kp2),pin(i,j,kp1),pin(i,j,k),pin(i,j,km1),pin(i,j,km2),w0(i,j,kp1)>=0.) -       &
                       w0(i,j,k)  *ip_weno(pin(i,j,kp2),pin(i,j,kp1),pin(i,j,k),pin(i,j,km1),pin(i,j,km2),pin(i,j,km3),w0(i,j,k)  >=0.))/dzf(k) &
                                      )
        end if
      end do !Loop over k
    end do !Loop over i
  end do !Loop over j

end subroutine

! Function that does the interpolation (ip) of the cell centered field values to the appropriate cell edges.
function ip_weno(vp2,vp1,v,vm1,vm2,vm3,lpos)
  implicit none
  real              :: ip_weno
  logical,intent(in):: lpos
  real,intent(in)   :: vp2,vp1,v,vm1,vm2,vm3
  ! NOTE: vp2 is not used if velocity is in the forward direction (i.e. wh >0)

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

  if (lpos) then !Positive velocity at cell face
    !compute smoothness indicators for each of the stencils
    beta(1) = c1*(vm3-2*vm2+vm1)**2 + c2*(vm3-4*vm2+3*vm1)**2
    beta(2) = c1*(vm2-2*vm1+v  )**2 + c2*(vm2-v)**2
    beta(3) = c1*(vm1-2*v+vp1  )**2 + c2*(3*vm1-4*v+vp1)**2

    !interpolated values of the variable at the cell faces using each of the stencils
    varFace(1) = (2*vm3- 7*vm2+ 11*vm1)/6
    varFace(2) = (- vm2+ 5*vm1+  2*v  )/6
    varFace(3) = (2*vm1+ 5*v  -    vp1)/6
  else !Negative velocity at cell face
    !compute smoothness indicators for each of the stencils
    !the following is found by mirroring the equations for positive velocity
    beta(1) = c1*(v  -2*vp1+vp2)**2 + c2*(3*v-4*vp1+vp2)**2
    beta(2) = c1*(vm1-2*v  +vp1)**2 + c2*(vm1-vp1)**2
    beta(3) = c1*(vm2-2*vm1+v  )**2 + c2*(vm2-4*vm1+3*v)**2

    !interpolated values of the variable at the cell faces using each of the stencils
    varFace(1) = (11*v   -7*vp1+ 2*vp2)/6
    varFace(2) = ( 2*vm1 +5*v  -   vp1)/6
    varFace(3) = ( - vm2 +5*vm1+ 2*v  )/6
  end if

  !compute weights
  wgt = wgtOpt*(epsWeno+beta)**(-pweno)
  wgtfac = sum(wgt)**(-1)

  ! compute interpolated value 
  ip_weno = sum(wgt(:)*varFace(:))*wgtfac

end function ip_weno

end module
