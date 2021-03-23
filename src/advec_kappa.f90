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
contains

subroutine advecc_kappa(putin,putout)
  use modglobal, only : i1,i2,ih,j1,j2,jh,k1,kmax,dxi,dyi,dzf
  use modfields, only : u0, v0, w0, rhobf
  implicit none
  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in) :: putin
  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: putout
  
  real      d1,d2,cf
  real :: d1m, d2m, d1p, cfm, cfp, work
  integer   i,j,k,k_low,k_high

  ! find the lowest and highest k level with a non-zero value in putin 
  k_low = -1
  do k=1,k1
     if (any(putin(:,:,k).ne.0.)) then
        k_low = k
        exit
     endif
  enddo
  if (k_low == -1) then
     ! putin == zero
     return
  endif

  k_high = kmax
  do k=k1,1,-1
     if (any(putin(:,:,k).ne.0.)) then
        k_high = k
        exit
     endif
  enddo

  if (k_low <= 3) then
     ! from layer 1 to 2, special case. k=2
     do j=2,j1
        do i=2,i1 ! YES
           d1m = 0
           d2m = rhobf(2) * putin(i,j,2) - rhobf(1) * putin(i,j,1)
           cfm = rhobf(1) * putin(i,j,1)
           d1p = rhobf(2) * putin(i,j,2) - rhobf(3) * putin(i,j,3)
           !d2p = rhobf(1) * putin(i,j,1) - rhobf(2) * putin(i,j,2  ) ! d2p = -d2m
           cfp = rhobf(2) * putin(i,j,2)

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
           putout(i,j,1) = putout(i,j,1) - (1./(rhobf(1)*dzf(1)))*work
           putout(i,j,2) = putout(i,j,2) + (1./(rhobf(2)*dzf(2)))*work
        end do
     end do
  end if

  !do k=1,kmax
  do k= max(k_low-1,1), min(k_high+2, kmax) ! loop accesses k-2, k-1, k, k+1
     do j=2,j1
        do i=2,i2 ! YES
           d2m =  putin(i  ,j,k) -putin(i-1,j,k)
           ! d2p = -putin(i  ,j,k) +putin(i-1,j,k) ! d2p = -d2m
           d1m = putin(i-1,j,k)-putin(i-2,j,k)
           d1p = putin(i  ,j,k)-putin(i+1,j,k)
           cfm = putin(i-1,j,k)
           cfp = putin(i  ,j,k)

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
           putout(i-1,j,k) = putout(i-1,j,k) - work
           putout(i,j,k)   = putout(i,j,k)   + work
        end do
     end do
     !  end do

     !  do k=1,kmax
     do j=2,j2
        do i=2,i1 ! YES
           d1m = putin(i,j-1,k)-putin(i,j-2,k)
           d1p = putin(i,j  ,k)-putin(i,j+1,k)
           d2m = putin(i,j  ,k)-putin(i,j-1,k)
           !d2p = putin(i,j-1,k)-putin(i,j  ,k) ! d2p = -d2m
           cfm = putin(i,j-1,k)
           cfp = putin(i,j  ,k)

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
           putout(i,j-1,k) = putout(i,j-1,k) - work
           putout(i,j,k)   = putout(i,j,k)   + work
        end do
     end do
     !  end do

     !  do k=3,kmax
     if (k >= 3) then
        do j=2,j1
           do i=2,i1 ! YES
              d1m = rhobf(k-1) * putin(i,j,k-1) - rhobf(k-2) * putin(i,j,k-2)
              d2m = rhobf(k)   * putin(i,j,k  ) - rhobf(k-1) * putin(i,j,k-1)
              d1p = rhobf(k)   * putin(i,j,k  ) - rhobf(k+1) * putin(i,j,k+1)
              ! d2p = rhobf(k-1) * putin(i,j,k-1) - rhobf(k)   * putin(i,j,k  ) ! d2p = -d2m
              cfm = rhobf(k-1) * putin(i,j,k-1)
              cfp = rhobf(k)   * putin(i,j,k  )

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
              putout(i,j,k-1) = putout(i,j,k-1) - (1./(rhobf(k-1)*dzf(k-1)))*work
              putout(i,j,k)   = putout(i,j,k)   + (1./(rhobf(k)  *dzf(k)  ))*work
           end do
        end do
     end if
  end do
end subroutine advecc_kappa


subroutine advecc_kappa_old(putin,putout)

  use modglobal, only : i1,i2,ih,j1,j2,jh,k1,kmax,dxi,dyi,dzf
  use modfields, only : u0, v0, w0, rhobf
  implicit none
  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in) :: putin
  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: putout
  real      d1,d2,cf

  integer   i,j,k

  do k=1,kmax
    do j=2,j1
      do i=2,i2
        if (u0(i,j,k)>0) then
          d1 = putin(i-1,j,k)-putin(i-2,j,k)
          d2 = putin(i  ,j,k)-putin(i-1,j,k)
          cf = putin(i-1,j,k)
        else
          d1 = putin(i  ,j,k)-putin(i+1,j,k)
          d2 = putin(i-1,j,k)-putin(i  ,j,k)
          cf = putin(i  ,j,k)
        end if
        cf = cf + rlim(d1,d2)
        putout(i-1,j,k) = putout(i-1,j,k) - cf * u0(i,j,k) * dxi
        putout(i,j,k)   = putout(i,j,k)   + cf * u0(i,j,k) * dxi
      end do
    end do
  end do

  do k=1,kmax
    do j=2,j2
      do i=2,i1
        if (v0(i,j,k)>0) then
          d1 = putin(i,j-1,k)-putin(i,j-2,k)
          d2 = putin(i,j  ,k)-putin(i,j-1,k)
          cf = putin(i,j-1,k)
        else
          d1 = putin(i,j  ,k)-putin(i,j+1,k)
          d2 = putin(i,j-1,k)-putin(i,j  ,k)
          cf = putin(i,j  ,k)
        end if
        cf = cf + rlim(d1,d2)
        putout(i,j-1,k) = putout(i,j-1,k) - cf * v0(i,j,k) * dyi
        putout(i,j,k)   = putout(i,j,k)   + cf * v0(i,j,k) * dyi
      end do
    end do
  end do

  do k=3,kmax
    do j=2,j1
      do i=2,i1
        if (w0(i,j,k)>0) then
          d1 = rhobf(k-1) * putin(i,j,k-1) - rhobf(k-2) * putin(i,j,k-2)
          d2 = rhobf(k)   * putin(i,j,k  ) - rhobf(k-1) * putin(i,j,k-1)
          cf = rhobf(k-1) * putin(i,j,k-1)
        else
          d1 = rhobf(k)   * putin(i,j,k  ) - rhobf(k+1) * putin(i,j,k+1)
          d2 = rhobf(k-1) * putin(i,j,k-1) - rhobf(k)   * putin(i,j,k  )
          cf = rhobf(k)   * putin(i,j,k  )
        end if
        cf = cf + rlim(d1,d2)
        putout(i,j,k-1) = putout(i,j,k-1) - (1./rhobf(k-1))*cf * w0(i,j,k) * (1.0 / dzf(k-1))
        putout(i,j,k)   = putout(i,j,k)   + (1./rhobf(k))*cf   * w0(i,j,k) * (1.0 / dzf(k))
      end do
    end do
  end do

  !special case for vertical advection at k=2
  do j=2,j1
    do i=2,i1
      if (w0(i,j,2)>0) then
        d1 = 0
        d2 = rhobf(2) * putin(i,j,2)   - rhobf(1) * putin(i,j,1)
        cf = rhobf(1) * putin(i,j,1)
      else
        d1 = rhobf(2) * putin(i,j,2)   - rhobf(3) * putin(i,j,3)
        d2 = rhobf(1) * putin(i,j,1)   - rhobf(2) * putin(i,j,2)
        cf = rhobf(2) * putin(i,j,2)
      end if
      cf = cf + rlim(d1,d2)
      putout(i,j,1) = putout(i,j,1) - (1./rhobf(1))*cf * w0(i,j,2) * (1.0 / dzf(1))
      putout(i,j,2) = putout(i,j,2) + (1./rhobf(2))*cf * w0(i,j,2) * (1.0 / dzf(2))
    end do
  end do

  end subroutine advecc_kappa_old

subroutine  halflev_kappa(putin,putout)

  use modglobal, only : i1,ih,j1,jh,k1
    use modfields, only : w0, rhobf, rhobh
    implicit none
    real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in) :: putin
    real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: putout
    real      d1,d2,cf
    integer   i,j,k


    do k=3,k1
      do j=2,j1
        do i=2,i1
          if (w0(i,j,k)>=0) then
            d1 = rhobf(k-1) * putin(i,j,k-1) - rhobf(k-2) * putin(i,j,k-2)
            d2 = rhobf(k) * putin(i,j,k  )   - rhobf(k-1) * putin(i,j,k-1)
            cf = rhobf(k-1) * putin(i,j,k-1)
          else
            d1 = rhobf(k)   * putin(i,j,k  ) - rhobf(k+1) * putin(i,j,k+1)
            d2 = rhobf(k-1) * putin(i,j,k-1) - rhobf(k)   * putin(i,j,k  )
            cf = rhobf(k)   * putin(i,j,k  )
          end if
          putout(i,j,k) = (1./rhobh(k))*(cf + rlim(d1,d2))
        end do
      end do
    end do

    do j=2,j1
      do i=2,i1
        if (w0(i,j,2)>=0) then
          d1 = 0
          d2 = rhobf(2) * putin(i,j,2) - rhobf(1) * putin(i,j,1)
          cf = rhobf(1) * putin(i,j,1)
        else
          d1 = rhobf(2) * putin(i,j,2) - rhobf(3) * putin(i,j,3)
          d2 = rhobf(1) * putin(i,j,1) - rhobf(2) * putin(i,j,2)
          cf = rhobf(2) * putin(i,j,2)
        end if
        putout(i,j,2) = (1./rhobh(2))*(cf + rlim(d1,d2))
      end do
    end do

  end subroutine halflev_kappa

!> Determination of the limiter function
  real function rlim(d1,d2)
    use modglobal, only : eps1
    implicit none
    real, intent(in) :: d1 !< Scalar flux at 1.5 cells upwind
    real, intent(in) :: d2 !< Scalar flux at 0.5 cells upwind

    real ri,phir

    ri    = (d2+eps1)/(d1+eps1)
    phir  = max(0.,min(2.*ri,min(1./3.+2./3.*ri,2.)))
    rlim  = 0.5*phir*d1
    end function rlim



end module advec_kappa
