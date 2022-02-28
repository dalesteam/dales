!> \file advec_kappa.f90
!!  Does advection with a kappa limiter scheme.
!! \par Revision list
!! \par Authors
!! \see Hundsdorfer et al 1995
!!
!! For advection of scalars that need to be strictly monotone (for example chemically reacting species)
!! the kappa scheme has been implemented:
!! \latexonly
!! \begin{eqnarray}
!!  F_{i-\frac{1}{2}}^{\kappa} &=& \fav{u}_{i-\frac{1}{2}}
!!  \left[\phi_{i-1}+\frac{1}{2}\kappa_{i-\frac{1}{2}}\left(\phi_{i-1}-\phi_{i-2}\right)\right],
!! \end{eqnarray}
!! in case $\fav{u}>0$. $\kappa_{i-\smfrac{1}{2}}$ serves as a switch between higher order advection and
!! first order upwind in case of strong upwind gradients of $\phi$.
!! \endlatexonly
!! This makes the scheme monotone, but also rather dissipative.
!!
!! limiter phi(r) = max(0, min(2*r, 2, K(r)))   (20) in Hundsdorfer 1995
!! K(r) = 1./3.+2./3.*r here -> kappa = 1/3 -> third-order upwind-biased scheme
!!
!! Changes 2020 by Jisk Attema and Fredrik Jansson:
!! - support for non-uniform vertical grid by replacing dzi by 1/dzf(k).
!!   neither the gradient nor the limiter has been modified.
!! - vectorization
!!   - rlim function inlined and rewritten without division and eps1
!!   - both branches of if uvw0 > 0 calculated, then selected
!!     in order to enable vectorization
!!   - merge k-loops of the x,y,z advection steps for better cache efficiency
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
module advec_kappa
use modprecision, only : field_r
contains

subroutine advecc_kappa(a,p)
  use modglobal, only : i1,i2,ih,j1,j2,jh,k1,kmax,dxi,dyi,dzf
  use modfields, only : u0, v0, w0, rhobf
  implicit none
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in) :: a
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: p
  
  real      d1,d2,cf
  real :: d1m, d2m, d1p, cfm, cfp, work
  integer   i,j,k,k_low,k_high

  ! find the lowest and highest k level with a non-zero value in a 
  k_low = -1
  do k=1,k1
     if (any(a(:,:,k).ne.0.)) then
        k_low = k
        exit
     endif
  enddo
  if (k_low == -1) then
     ! a == zero
     return
  endif

  do k=k1,1,-1
     if (any(a(:,:,k).ne.0.)) then
        k_high = k
        exit
     endif
  enddo

  if (k_low <= 3) then
     ! vertical advection from layer 1 to 2, special case. k=2
     do j=2,j1
        do i=2,i1 ! YES
           d1m = 0
           d2m = rhobf(2) * a(i,j,2) - rhobf(1) * a(i,j,1)
           cfm = rhobf(1) * a(i,j,1)
           d1p = rhobf(2) * a(i,j,2) - rhobf(3) * a(i,j,3)
           cfp = rhobf(2) * a(i,j,2)

           if (w0(i,j,2) > 0) then
              d1 = d1m
              d2 = d2m
              cf = cfm
           else
              d1 = d1p
              d2 = -d2m
              cf = cfp
           end if

           work = cf + &
                min(abs(d1), abs(d2), abs((d1/6.0) + (d2/3.0))) * &
                (sign(0.5, d1) + sign(0.5, d2))

           work = work * w0(i,j,2)
           p(i,j,1) = p(i,j,1) - (1./(rhobf(1)*dzf(1)))*work
           p(i,j,2) = p(i,j,2) + (1./(rhobf(2)*dzf(2)))*work
        end do
     end do
  end if

  !do k=1,kmax
  do k= max(k_low-1,1), min(k_high+2, kmax) ! loop accesses k-2, k-1, k, k+1
     do j=2,j1
        do i=2,i2 ! YES
           d2m = a(i  ,j,k)-a(i-1,j,k)
           d1m = a(i-1,j,k)-a(i-2,j,k)
           d1p = a(i  ,j,k)-a(i+1,j,k)
           cfm = a(i-1,j,k)
           cfp = a(i  ,j,k)

           if (u0(i,j,k) > 0) then
              d1 = d1m
              d2 = d2m
              cf = cfm
           else
              d1 = d1p
              d2 = -d2m
              cf = cfp
           end if

           work = cf + &
                min(abs(d1), abs(d2), abs((d1/6.0) + (d2/3.0))) * &
                (sign(0.5, d1) + sign(0.5, d2))

           work = work * u0(i,j,k) * dxi
           p(i-1,j,k) = p(i-1,j,k) - work
           p(i,j,k)   = p(i,j,k)   + work
        end do
     end do
     !  end do

     !  do k=1,kmax
     do j=2,j2
        do i=2,i1 ! YES
           d1m = a(i,j-1,k)-a(i,j-2,k)
           d1p = a(i,j  ,k)-a(i,j+1,k)
           d2m = a(i,j  ,k)-a(i,j-1,k)
           cfm = a(i,j-1,k)
           cfp = a(i,j  ,k)

           if (v0(i,j,k) > 0) then
              d1 = d1m
              d2 = d2m
              cf = cfm
           else
              d1 = d1p
              d2 = -d2m
              cf = cfp
           end if

           work = cf + &
                min(abs(d1), abs(d2), abs((d1/6.0) + (d2/3.0))) * &
                (sign(0.5, d1) + sign(0.5, d2))

           work = work * v0(i,j,k) * dyi
           p(i,j-1,k) = p(i,j-1,k) - work
           p(i,j,k)   = p(i,j,k)   + work
        end do
     end do
     !  end do

     !  do k=3,kmax
     if (k >= 3) then
        do j=2,j1
           do i=2,i1 ! YES
              d1m = rhobf(k-1) * a(i,j,k-1) - rhobf(k-2) * a(i,j,k-2)
              d2m = rhobf(k)   * a(i,j,k  ) - rhobf(k-1) * a(i,j,k-1)
              d1p = rhobf(k)   * a(i,j,k  ) - rhobf(k+1) * a(i,j,k+1)
              cfm = rhobf(k-1) * a(i,j,k-1)
              cfp = rhobf(k)   * a(i,j,k  )

              if (w0(i,j,k) > 0) then
                 d1 = d1m
                 d2 = d2m
                 cf = cfm
              else
                 d1 = d1p
                 d2 = -d2m
                 cf = cfp
              end if

              work = cf + &
                   min(abs(d1), abs(d2), abs((d1/6.0) + (d2/3.0))) * &
                   (sign(0.5, d1) + sign(0.5, d2))

              work = work * w0(i,j,k)
              p(i,j,k-1) = p(i,j,k-1) - (1./(rhobf(k-1)*dzf(k-1)))*work
              p(i,j,k)   = p(i,j,k)   + (1./(rhobf(k)  *dzf(k)  ))*work
           end do
        end do
     end if
  end do
end subroutine advecc_kappa

subroutine  halflev_kappa(a,p)

  use modglobal, only : i1,ih,j1,jh,k1
    use modfields, only : w0, rhobf, rhobh
    implicit none
    real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in) :: a
    real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: p
    real(field_r)      d1,d2,cf
    integer   i,j,k


    do k=3,k1
      do j=2,j1
        do i=2,i1
          if (w0(i,j,k)>=0) then
            d1 = rhobf(k-1) * a(i,j,k-1) - rhobf(k-2) * a(i,j,k-2)
            d2 = rhobf(k) * a(i,j,k  )   - rhobf(k-1) * a(i,j,k-1)
            cf = rhobf(k-1) * a(i,j,k-1)
          else
            d1 = rhobf(k)   * a(i,j,k  ) - rhobf(k+1) * a(i,j,k+1)
            d2 = rhobf(k-1) * a(i,j,k-1) - rhobf(k)   * a(i,j,k  )
            cf = rhobf(k)   * a(i,j,k  )
          end if
          p(i,j,k) = (1./rhobh(k))*(cf + rlim(d1,d2))
        end do
      end do
    end do

    do j=2,j1
      do i=2,i1
        if (w0(i,j,2)>=0) then
          d1 = 0
          d2 = rhobf(2) * a(i,j,2) - rhobf(1) * a(i,j,1)
          cf = rhobf(1) * a(i,j,1)
        else
          d1 = rhobf(2) * a(i,j,2) - rhobf(3) * a(i,j,3)
          d2 = rhobf(1) * a(i,j,1) - rhobf(2) * a(i,j,2)
          cf = rhobf(2) * a(i,j,2)
        end if
        p(i,j,2) = (1./rhobh(2))*(cf + rlim(d1,d2))
      end do
    end do

  end subroutine halflev_kappa

!> Determination of the limiter function
  real function rlim(d1,d2)
    use modglobal, only : eps1
    implicit none
    real(field_r), intent(in) :: d1 !< Scalar flux at 1.5 cells upwind
    real(field_r), intent(in) :: d2 !< Scalar flux at 0.5 cells upwind

    real(field_r) ri,phir

    ri    = (d2+eps1)/(d1+eps1)
    phir  = max(0.,min(2.*ri,min(1./3.+2./3.*ri,2.)))
    rlim  = 0.5*phir*d1
    end function rlim



end module advec_kappa
