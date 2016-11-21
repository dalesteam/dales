!> \file modpois.f90
!!  Solves the Poisson equation for the pressure fluctuations

!>
!!  Solves the Poisson equation for the pressure fluctuations
!>
!!  \author Harm Jonker, TU Delft
!!  \author Hans Cuijpers, IMAU
!!  \todo documentation
!!  \par Revision list
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
!            2014      Netherlands eScience Center
!

module modpois

implicit none
private
public :: initpois,poisson,exitpois,p
save

  real,allocatable :: p(:,:,:)                            ! pressure fluctuations
  real,allocatable :: xyrt(:,:)                           ! constant factors in poisson equation

contains
  subroutine initpois
    use modglobal, only : i1,j1,kmax,imax,jmax,itot,jtot,dxi,dyi,pi,ih,jh
    use modfft2d,  only : fft2dinit
    use modmpi, only    : myidx, myidy

    implicit none

    integer :: i,j,iv,jv
    real    :: fac
    real    :: xrt(itot), yrt(jtot)

    allocate(   p(2-ih:i1+ih,2-jh:j1+jh,kmax))
    allocate(xyrt(2-ih:i1+ih,2-jh:j1+jh)     )

  ! Generate Eigenvalues xrt and yrt

  !  I --> direction

    fac = 1./(2.*itot)
    do i=3,itot,2
      xrt(i-1)=-4.*dxi*dxi*(sin(float((i-1))*pi*fac))**2
      xrt(i)  = xrt(i-1)
    end do
    xrt(1    ) = 0.
    xrt(itot ) = -4.*dxi*dxi

  !  J --> direction

    fac = 1./(2.*jtot)
    do j=3,jtot,2
      yrt(j-1)=-4.*dyi*dyi*(sin(float((j-1))*pi*fac))**2
      yrt(j  )= yrt(j-1)
    end do
    yrt(1    ) = 0.
    yrt(jtot ) = -4.*dyi*dyi

  ! Combine I and J directions
  ! Note that:
  ! 1. MPI starts counting at 0 so it should be myidy * jmax
  ! 2. real data, ie. no halo, starts at index 2 in the array xyrt(2,2) <-> xrt(1), yrt(1)

    do j=2,j1
      jv = j + myidy*jmax - 1
      do i=2,i1
        iv = i + myidx*imax - 1
        xyrt(i,j)=(xrt(iv)+yrt(jv)) !!! LH
      end do
    end do

    call fft2dinit()

  end subroutine initpois

  subroutine poisson
    implicit none

    call fillps
    call solmpj
    call tderive

  end subroutine poisson

  subroutine exitpois
    use modfft2d,  only : fft2dexit
    implicit none

    deallocate(p,xyrt)

    call fft2dexit()

  end subroutine exitpois
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fillps

  ! Chiel van Heerwaarden,  19 June 2007
  ! Adapted fillps for RK3 time loop


    use modfields, only : up, vp, wp, um, vm, wm, rhobf,rhobh
    use modglobal, only : rk3step, i1,j1,kmax,k1, dx,dy,dzf,rdt
    use modmpi,    only : excj
    implicit none
    real,allocatable :: pup(:,:,:), pvp(:,:,:), pwp(:,:,:)
    integer i,j,k
    real rk3coef

    allocate(pup(2-1:i1+1,2-1:j1+1,k1))
    allocate(pvp(2-1:i1+1,2-1:j1+1,k1))
    allocate(pwp(2-1:i1+1,2-1:j1+1,k1))

    rk3coef = rdt / (4. - dble(rk3step))

    do k=1,kmax
      do j=2,j1
        do i=2,i1
          pup(i,j,k) = up(i,j,k) + um(i,j,k) / rk3coef
          pvp(i,j,k) = vp(i,j,k) + vm(i,j,k) / rk3coef
          pwp(i,j,k) = wp(i,j,k) + wm(i,j,k) / rk3coef
        end do
      end do
    end do


  !****************************************************************

  !     Fill the right hand for the poisson solver.
  !     Call excj to set the values in the halo zone.
  !     Also we take wp(i,j,1) and wp(i,j,k1) equal to zero.

  !**************************************************************

    do j=2,j1
      do i=2,i1
        pwp(i,j,1)  = 0.
        pwp(i,j,k1) = 0.
      end do
    end do

    call excj( pup, 2-1, i1+1, 2-1, j1+1, 1, k1 )
    call excj( pvp, 2-1, i1+1, 2-1, j1+1, 1, k1 )
    call excj( pwp, 2-1, i1+1, 2-1, j1+1, 1, k1 )

    do k=1,kmax
      do j=2,j1
        do i=2,i1
          p(i,j,k)  =  rhobf(k)*(( pup(i+1,j,k)-pup(i,j,k) ) / dx &
                                +( pvp(i,j+1,k)-pvp(i,j,k) ) / dy ) &
                      +( pwp(i,j,k+1)*rhobh(k+1)-pwp(i,j,k)*rhobh(k)) / dzf(k)
        end do
      end do
    end do

    deallocate( pup,pvp,pwp )

  end subroutine fillps

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine solmpj
! version: working version, barrou's removed,
!          correct timing fft's
!          AAPC with MPI-provided routines
!          uses only 2 AAPC, using MPI-rovided routines,
!          to pre-distribute the arrays s.t. complete
!          2-D planes are present on each processor
!          uses ALLTOALL instead of ALLTOALLV
!          ONLY distribution in j-direction allowed

! NOTE: input array p is supposed to have the ip1ray distribution,
!       i.e. the entire range of the first index must be present on
!       each processor

!******************************************************************
!********************  FAST POISSON SOLVER ************************
!*****                                                        *****
!***               P_xx + P_yy + P_zz  =f(x,y,z)                ***
!****                                                         *****
!******************************************************************
!   FOURIER TRANSFORMS IN X AND Y DIRECTION   GIVE:
!   a^2 P + b^2 P + P_zz =  F(x,y,z) = FFT_i [ FTT_j (f(x,y,z))]

!   where a and b are the KNOWN eigenvalues, and P_zz is

!   P_zz =[ P_{i,j,k+1} - 2 P_{i,j,k} +P_{i,j,k-1} ] / (dz * dz)

!   a^2 P + b^2 +P_zz =
!   [P_{i,j,k+1}-(2+a^2+ b^2) P_{i,j,k}+P_{i,j,k-1}]/(dz*dz)=F( x,y,z)

!   The equation above results in a tridiagonal system in k which
!   can be solved with Gaussian elemination --> P
!   The P we have found with the Gaussian elemination is still in
!   the Fourier Space and 2 backward FFTS are necessary to compute
!   the physical P
!******************************************************************
!******************************************************************
!******************************************************************
!****   Programmer: Bendiks Jan Boersma                      ******
!****               email : b.j.boersma@wbmt.tudelft.nl      ******
!****                                                        ******
!****   USES      :  VFFTPACK   (netlib)                     ******
!****             :  FFTPACK    (netlib)                     ******
!****                (B.J. Boersma & L.J.P. Timmermans)      ******
!****                                                        ******
!******************************************************************
!******************************************************************

! mpi-version, no master region for timing
!              copy times all included

    use modfft2d,  only : fft2df, fft2db
    use modglobal, only : kmax,dzf,dzh,i1,j1,ih,jh
    use modfields, only : rhobf, rhobh
    implicit none

    real :: a(kmax),b(kmax),c(kmax)

  ! allocate d in the same shape as p and xyrt
    real, allocatable :: d(:,:,:)

    real    z,ak,bk,bbk
    integer i, j, k

    allocate(d(2-ih:i1+ih,2-jh:j1+jh,kmax))
    
  ! Forward FFT
    call fft2df(p,ih,jh)

  ! Generate tridiagonal matrix

    do k=1,kmax
      ! SB fixed the coefficients
      a(k)=rhobh(k)  /(dzf(k)*dzh(k  ))
      c(k)=rhobh(k+1)/(dzf(k)*dzh(k+1))
      b(k)=-(a(k)+c(k))
    end do

    b(1   )=b(1)+a(1)
    a(1   )=0.
    b(kmax)=b(kmax)+c(kmax)
    c(kmax)=0.

    ! SOLVE TRIDIAGONAL SYSTEMS WITH GAUSSIAN ELEMINATION
    do j=2,j1
      do i=2,i1
        z        = 1./(b(1)+rhobf(1)*xyrt(i,j))
        d(i,j,1) = c(1)*z
        p(i,j,1) = p(i,j,1)*z
      end do
    end do

    do k=2,kmax-1
      do  j=2,j1
        do  i=2,i1
          bbk      = b(k)+rhobf(k)*xyrt(i,j)
          z        = 1./(bbk-a(k)*d(i,j,k-1))
          d(i,j,k) = c(k)*z
          p(i,j,k) = (p(i,j,k)-a(k)*p(i,j,k-1))*z
        end do
      end do
    end do

    ak =a(kmax)
    bk =b(kmax)
    do j=2,j1
      do i=2,i1
        bbk = bk + rhobf(kmax)*xyrt(i,j)
        z        = bbk-ak*d(i,j,kmax-1)
        if(z/=0.) then
          p(i,j,kmax) = (p(i,j,kmax)-ak*p(i,j,kmax-1))/z
        else
          p(i,j,kmax) =0.
        end if
      end do
    end do

    do k=kmax-1,1,-1
      do j=2,j1
        do i=2,i1
          p(i,j,k) = p(i,j,k)-d(i,j,k)*p(i,j,k+1)
        end do
      end do
    end do

    ! Backward FFT
    call fft2db(p,ih,jh)

    deallocate(d)
    return    
  end subroutine solmpj
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine tderive

!-----------------------------------------------------------------|
!                                                                 |
!*** *tderive*  read input fields for initialisation              |
!                                                                 |
!      Hans Cuijpers   I.M.A.U.     06/01/1995                    |
!                                                                 |
!     purpose.                                                    |
!     --------                                                    |
!                                                                 |
!     Refill array p with pressure values.                        |
!                                                                 |
!     Further we set cyclic boundary conditions for the pressure- |
!     fluctuations in the x-y plane.                              |
!                                                                 |
!**   interface.                                                  |
!     ----------                                                  |
!                                                                 |
!             *tderive* is called from *program*.                 |
!                                                                 |
!-----------------------------------------------------------------|

    use modfields, only : up, vp, wp
    use modglobal, only : i1,j1,kmax,dx,dy,dzh,ih,jh
    use modmpi,    only : excjs
    implicit none
    integer i,j,k

  ! **  Cyclic boundary conditions **************
  ! **  set by the commcart communication in excj
  call excjs( p, 2,i1,2,j1,1,kmax,ih,jh)

  !*****************************************************************
  ! **  Calculate time-derivative for the velocities with known ****
  ! **  pressure gradients.  ***************************************
  !*****************************************************************

    do k=1,kmax
    do j=2,j1
    do i=2,i1
      up(i,j,k) = up(i,j,k)-(p(i,j,k)-p(i-1,j,k))/dx
      vp(i,j,k) = vp(i,j,k)-(p(i,j,k)-p(i,j-1,k))/dy
    end do
    end do
    end do

    do k=2,kmax
    do j=2,j1
    do i=2,i1
      wp(i,j,k) = wp(i,j,k)-(p(i,j,k)-p(i,j,k-1))/dzh(k)
    end do
    end do
    end do

    return
  end subroutine tderive

end module modpois
