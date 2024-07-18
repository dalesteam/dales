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

!< Horizontal advection with the kappa scheme

subroutine hadvecc_kappa(a_in,a_out)
  use modglobal, only : i1,i2,ih,j1,j2,jh,k1,kmax,dxi,dyi
  use modfields, only : u0, v0
  implicit none
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in) :: a_in
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: a_out

  real(field_r) :: d1, d2, cf
  real(field_r) :: d1m, d2m, d1p, cfm, cfp, work
  integer :: i, j, k, k_low, k_high, kbeg, kend

  ! find the lowest and highest k level with a non-zero value in a_in

#if defined(DALES_GPU)
  ! This version is faster with OpenACC acceleration as it enables collapse(3)
  k_low = k1 + 1

  !$acc parallel loop collapse(3) default(present) reduction(min:k_low)
  do k = 1, k1
    do j = 2, j1
      do i = 2, i1
        if (a_in(i,j,k).ne.0.) then
          k_low = min(k,k_low)
        endif
      enddo
    enddo
  enddo

  ! a_in == zero
  if (k_low == (k1 + 1)) then
    return
  endif

  k_high = 0
  !$acc parallel loop collapse(3) default(present) reduction(max:k_high)
  do k = k_low, k1
    do j = 2, j1
      do i = 2, i1
        if ((a_in(i,j,k).ne.0.)) then
          k_high = max(k,k_high)
        endif
      enddo
    enddo
  enddo

#else
  k_low = -1
  do k=1,k1
     if (any(a_in(:,:,k).ne.0.)) then
        k_low = k
        exit
     endif
  enddo
  if (k_low == -1) then
     ! a_in == zero
     return
  endif

  do k=k1,1,-1
     if (any(a_in(:,:,k).ne.0.)) then
        k_high = k
        exit
     endif
  enddo
#endif

  ! loop accesses k-2, k-1, k, k+1
  kbeg = max(k_low-1,1)
  kend = min(k_high+2, kmax)

  !$acc parallel loop collapse(3) default(present) async(2)
  do k = kbeg, kend
    do j = 2, j1
      do i = 2, i2
        d2m = a_in(i  ,j,k)-a_in(i-1,j,k)
        d1m = a_in(i-1,j,k)-a_in(i-2,j,k)
        d1p = a_in(i  ,j,k)-a_in(i+1,j,k)
        cfm = a_in(i-1,j,k)
        cfp = a_in(i  ,j,k)

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
             min(abs(d1), abs(d2), abs((d1/6.0_field_r) + (d2/3.0_field_r))) * &
             (sign(0.5_field_r, d1) + sign(0.5_field_r, d2))

        work = work * u0(i,j,k) * dxi
        !$acc atomic update
        a_out(i-1,j,k) = a_out(i-1,j,k) - work
        !$acc atomic update
        a_out(i,j,k)   = a_out(i,j,k)   + work
      end do
    end do
  end do

  !$acc parallel loop collapse(3) default(present) async(3)
  do k = kbeg, kend
    do j = 2, j2
      do i = 2, i1
        d1m = a_in(i,j-1,k)-a_in(i,j-2,k)
        d1p = a_in(i,j  ,k)-a_in(i,j+1,k)
        d2m = a_in(i,j  ,k)-a_in(i,j-1,k)
        cfm = a_in(i,j-1,k)
        cfp = a_in(i,j  ,k)

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
             min(abs(d1), abs(d2), abs((d1/6.0_field_r) + (d2/3.0_field_r))) * &
             (sign(0.5_field_r, d1) + sign(0.5_field_r, d2))

        work = work * v0(i,j,k) * dyi
        !$acc atomic update
        a_out(i,j-1,k) = a_out(i,j-1,k) - work
        !$acc atomic update
        a_out(i,j,k)   = a_out(i,j,k)   + work
      end do
    end do
  end do

end subroutine hadvecc_kappa

!< Vertical advection with the kappa scheme

subroutine vadvecc_kappa(a_in,a_out)
  use modglobal, only : i1,i2,ih,j1,j2,jh,k1,kmax,dzf
  use modfields, only : w0, rhobf
  implicit none
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in) :: a_in
  real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: a_out

  real(field_r) :: d1, d2, cf
  real(field_r) :: d1m, d2m, d1p, cfm, cfp, work
  integer :: i, j, k, k_low, k_high, kend

  ! find the lowest and highest k level with a non-zero value in a_in
  ! TODO: Unnecessary to do this twice every time step, should do it when hadvecc is called and then reuse

#if defined(DALES_GPU)
  ! This version is faster with OpenACC acceleration as it enables collapse(3)
  k_low = k1 + 1

  !$acc parallel loop collapse(3) default(present) reduction(min:k_low)
  do k = 1, k1
    do j = 2, j1
      do i = 2, i1
        if (a_in(i,j,k).ne.0.) then
          k_low = min(k,k_low)
        endif
      enddo
    enddo
  enddo

  ! a_in == zero
  if (k_low == (k1 + 1)) then
    return
  endif

  k_high = 0
  !$acc parallel loop collapse(3) default(present) reduction(max:k_high)
  do k = k_low, k1
    do j = 2, j1
      do i = 2, i1
        if ((a_in(i,j,k).ne.0.)) then
          k_high = max(k,k_high)
        endif
      enddo
    enddo
  enddo

#else
  k_low = -1
  do k=1,k1
     if (any(a_in(:,:,k).ne.0.)) then
        k_low = k
        exit
     endif
  enddo
  if (k_low == -1) then
     ! a_in == zero
     return
  endif

  do k=k1,1,-1
     if (any(a_in(:,:,k).ne.0.)) then
        k_high = k
        exit
     endif
  enddo
#endif

  ! vertical advection from layer 1 to 2, special case. k=2
  if (k_low <= 3) then
    !$acc parallel loop collapse(2) default(present) async(1)
    do j = 2, j1
      do i = 2, i1
        d1m = 0
        d2m = rhobf(2) * a_in(i,j,2) - rhobf(1) * a_in(i,j,1)
        cfm = rhobf(1) * a_in(i,j,1)
        d1p = rhobf(2) * a_in(i,j,2) - rhobf(3) * a_in(i,j,3)
        cfp = rhobf(2) * a_in(i,j,2)

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
             min(abs(d1), abs(d2), abs((d1/6.0_field_r) + (d2/3.0_field_r))) * &
             (sign(0.5_field_r, d1) + sign(0.5_field_r, d2))

        work = work * w0(i,j,2)
        !$acc atomic update
        a_out(i,j,1) = a_out(i,j,1) - (1/(rhobf(1)*dzf(1)))*work
        !$acc atomic update
        a_out(i,j,2) = a_out(i,j,2) + (1/(rhobf(2)*dzf(2)))*work
      end do
    end do
  end if

  kend = min(k_high+2, kmax)

  !$acc parallel loop collapse(3) default(present) async(4)
  do k = 3, kend
    do j = 2, j1
      do i = 2, i1
        d1m = rhobf(k-1) * a_in(i,j,k-1) - rhobf(k-2) * a_in(i,j,k-2)
        d2m = rhobf(k)   * a_in(i,j,k  ) - rhobf(k-1) * a_in(i,j,k-1)
        d1p = rhobf(k)   * a_in(i,j,k  ) - rhobf(k+1) * a_in(i,j,k+1)
        cfm = rhobf(k-1) * a_in(i,j,k-1)
        cfp = rhobf(k)   * a_in(i,j,k  )

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
             min(abs(d1), abs(d2), abs((d1/6.0_field_r) + (d2/3.0_field_r))) * &
             (sign(0.5_field_r, d1) + sign(0.5_field_r, d2))

        work = work * w0(i,j,k)
        !$acc atomic update
        a_out(i,j,k-1) = a_out(i,j,k-1) - (1/(rhobf(k-1)*dzf(k-1)))*work
        !$acc atomic update
        a_out(i,j,k)   = a_out(i,j,k)   + (1/(rhobf(k)  *dzf(k)  ))*work
      end do
    end do
  end do
  !$acc wait

end subroutine vadvecc_kappa

subroutine  halflev_kappa(a_in,a_out)

    use modglobal, only : i1,ih,j1,jh,k1
    use modfields, only : w0, rhobf, rhobh
    implicit none
    real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in) :: a_in
    real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: a_out
    real(field_r)      d1,d2,cf
    integer   i,j,k

    !$acc parallel loop collapse(3) private(d1,d2,cf) default(present) async(1)
    do k = 3, k1
      do j = 2, j1
        do i = 2, i1
          if (w0(i,j,k)>=0) then
            d1 = rhobf(k-1) * a_in(i,j,k-1) - rhobf(k-2) * a_in(i,j,k-2)
            d2 = rhobf(k)   * a_in(i,j,k  ) - rhobf(k-1) * a_in(i,j,k-1)
            cf = rhobf(k-1) * a_in(i,j,k-1)
          else
            d1 = rhobf(k)   * a_in(i,j,k  ) - rhobf(k+1) * a_in(i,j,k+1)
            d2 = rhobf(k-1) * a_in(i,j,k-1) - rhobf(k)   * a_in(i,j,k  )
            cf = rhobf(k)   * a_in(i,j,k  )
          end if
          a_out(i,j,k) = (1/rhobh(k))*(cf + rlim(d1,d2))
        end do
      end do
    end do

    !$acc parallel loop collapse(2) private(d1,d2,cf) default(present) async(2)
    do j = 2, j1
      do i = 2, i1
        if (w0(i,j,2)>=0) then
          d1 = 0
          d2 = rhobf(2) * a_in(i,j,2) - rhobf(1) * a_in(i,j,1)
          cf = rhobf(1) * a_in(i,j,1)
        else
          d1 = rhobf(2) * a_in(i,j,2) - rhobf(3) * a_in(i,j,3)
          d2 = rhobf(1) * a_in(i,j,1) - rhobf(2) * a_in(i,j,2)
          cf = rhobf(2) * a_in(i,j,2)
        end if
        a_out(i,j,2) = (1/rhobh(2))*(cf + rlim(d1,d2))
      end do
    end do
    !$acc wait(1,2)

  end subroutine halflev_kappa

!> Determination of the limiter function
  real function rlim(d1,d2)
    use modglobal, only : eps1
    implicit none
    real(field_r), intent(in) :: d1 !< Scalar flux at 1.5 cells upwind
    real(field_r), intent(in) :: d2 !< Scalar flux at 0.5 cells upwind

    real(field_r) ri,phir

    ri    = (d2+eps1)/(d1+eps1)
    phir  = max(0._field_r,min(2._field_r*ri,min(1._field_r/3._field_r + &
                                         2._field_r/3._field_r*ri , 2._field_r)))
    rlim  = 0.5_field_r*phir*d1
    end function rlim

end module advec_kappa
