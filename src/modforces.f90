!> \file modforces.f90
!!  Calculates the other forces and sources in the equations.

!>
!!  Calculates the other forces and sources in the equations.
!>
!!  This includes the large scale forcings, the coriolis and the subsidence
!!  \author Pier Siebesma, K.N.M.I.
!!  \author Stephan de Roode,TU Delft
!!  \author Chiel van Heerwaarden, Wageningen U.R.
!!  \author Thijs Heus,MPI-M
!!  \author Steef Böing, TU Delft
!!  \par Revision list
!!  \todo Documentation
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



module modforces
!Calculates additional forces and large scale tendencies
implicit none
save
private
public :: forces, coriolis, lstend,lforce_user
logical :: lforce_user = .false.
contains
  subroutine forces

!-----------------------------------------------------------------|
!                                                                 |
!      Hans Cuijpers   I.M.A.U.                                   |
!      Pier Siebesma   K.N.M.I.     06/01/1995                    |
!                                                                 |
!     purpose.                                                    |
!     --------                                                    |
!                                                                 |
!      Calculates all other terms in the N-S equation,            |
!      except for the diffusion and advection terms.              |
!                                                                 |
!**   interface.                                                  |
!     ----------                                                  |
!                                                                 |
!     *forces* is called from *program*.                          |
!                                                                 |
!-----------------------------------------------------------------|

  use modglobal, only : i1,j1,kmax,dzh,dzf,grav, lpressgrad
  use modfields, only : sv0,up,vp,wp,thv0h,dpdxl,dpdyl,thvh
  use moduser,   only : force_user
  use modmicrodata, only : imicro, imicro_bulk, imicro_bin, imicro_sice,iqr
  implicit none

  integer i, j, k, jm, jp, km, kp

  if (lforce_user) call force_user

  if((imicro==imicro_sice).or.(imicro==imicro_bulk).or.(imicro==imicro_bin)) then
    do k=2,kmax
      kp=k+1
      km=k-1
    do j=2,j1
      jp=j+1
      jm=j-1
    do i=2,i1
    
    if (lpressgrad) then
      up(i,j,k) = up(i,j,k) - dpdxl(k)      !RN LS pressure gradient force in x,y directions;
      vp(i,j,k) = vp(i,j,k) - dpdyl(k)
    end if
      wp(i,j,k) = wp(i,j,k) + grav*(thv0h(i,j,k)-thvh(k))/thvh(k) - &
                  grav*(sv0(i,j,k,iqr)*dzf(k-1)+sv0(i,j,k-1,iqr)*dzf(k))/(2.0*dzh(k))
    end do
    end do
    end do
  else
    do k=2,kmax
      kp=k+1
      km=k-1
    do j=2,j1
      jp=j+1
      jm=j-1
    do i=2,i1
      up(i,j,k) = up(i,j,k) - dpdxl(k)
      vp(i,j,k) = vp(i,j,k) - dpdyl(k)
      wp(i,j,k) = wp(i,j,k) + grav*(thv0h(i,j,k)-thvh(k))/thvh(k)
    end do
    end do
    end do
  end if

!     --------------------------------------------
!     special treatment for lowest full level: k=1
!     --------------------------------------------

  do j=2,j1
    jp = j+1
    jm = j-1
  do i=2,i1

    if (lpressgrad) then
      up(i,j,1) = up(i,j,1) - dpdxl(1)
      vp(i,j,1) = vp(i,j,1) - dpdyl(1)
    end if

    wp(i,j,1) = 0.0

  end do
  end do
!     ----------------------------------------------end i,j-loop


  return
  end subroutine forces
  subroutine coriolis

!-----------------------------------------------------------------|
!                                                                 |
!      Thijs Heus TU Delft                                        |
!                                                                 |
!     purpose.                                                    |
!     --------                                                    |
!                                                                 |
!      Calculates the Coriolis force.                             |
!                                                                 |
!**   interface.                                                  |
!     ----------                                                  |
!                                                                 |
!     *coriolis* is called from *program*.                          |
!                                                                 |
!-----------------------------------------------------------------|

  use modglobal, only : i1,j1,kmax,dzh,dzf,cu,cv,om22,om23,lcoriol
  use modfields, only : u0,v0,w0,up,vp,wp
  implicit none

  integer i, j, k, jm, jp, km, kp

  if (lcoriol .eqv. .false.) return
  
  do k=2,kmax
    kp=k+1
    km=k-1
  do j=2,j1
    jp=j+1
    jm=j-1
  do i=2,i1

    up(i,j,k) = up(i,j,k)+ cv*om23 &
          +(v0(i,j,k)+v0(i,jp,k)+v0(i-1,j,k)+v0(i-1,jp,k))*om23*0.25 &
          -(w0(i,j,k)+w0(i,j,kp)+w0(i-1,j,kp)+w0(i-1,j,k))*om22*0.25

    vp(i,j,k) = vp(i,j,k)  - cu*om23 &
          -(u0(i,j,k)+u0(i,jm,k)+u0(i+1,jm,k)+u0(i+1,j,k))*om23*0.25


    wp(i,j,k) = wp(i,j,k) + cu*om22 +( (dzf(km) * (u0(i,j,k)  + u0(i+1,j,k) )    &
                +    dzf(k)  * (u0(i,j,km) + u0(i+1,j,km))  ) / dzh(k) ) &
                * om22*0.25
  end do
  end do
!     -------------------------------------------end i&j-loop
  end do
!     -------------------------------------------end k-loop

!     --------------------------------------------
!     special treatment for lowest full level: k=1
!     --------------------------------------------

  do j=2,j1
    jp = j+1
    jm = j-1
  do i=2,i1

    up(i,j,1) = up(i,j,1)  + cv*om23 &
          +(v0(i,j,1)+v0(i,jp,1)+v0(i-1,j,1)+v0(i-1,jp,1))*om23*0.25 &
          -(w0(i,j,1)+w0(i,j ,2)+w0(i-1,j,2)+w0(i-1,j ,1))*om22*0.25

    vp(i,j,1) = vp(i,j,1) - cu*om23 &
          -(u0(i,j,1)+u0(i,jm,1)+u0(i+1,jm,1)+u0(i+1,j,1))*om23*0.25

    wp(i,j,1) = 0.0

  end do
  end do
!     ----------------------------------------------end i,j-loop


  return
  end subroutine coriolis

  subroutine lstend

!-----------------------------------------------------------------|
!                                                                 |
!*** *lstend*  calculates large-scale tendencies                  |
!                                                                 |
!      Pier Siebesma   K.N.M.I.     06/01/1995                    |
!                                                                 |
!     purpose.                                                    |
!     --------                                                    |
!                                                                 |
!     calculates and adds large-scale tendencies due to           |
!     large scale advection and subsidence.                       |
!                                                                 |
!**   interface.                                                  |
!     ----------                                                  |
!                                                                 |
!             *lstend* is called from *program*.                  |
!                                                                 |
!-----------------------------------------------------------------|

  use modglobal, only : i1,j1,kmax,dzh,nsv,lmomsubs
  use modfields, only : up,vp,thlp,qtp,svp,&
                        whls, u0av,v0av,thl0,qt0,sv0,u0,v0,&
                        dudxls,dudyls,dvdxls,dvdyls,dthldxls,dthldyls,dqtdxls,dqtdyls, &
                        dqtdtls, dthldtls, dudtls, dvdtls
  implicit none

  integer k,kp,km

!     1. DETERMINE LARGE SCALE TENDENCIES
!        --------------------------------

!     1.1 lowest model level above surface : only downward component
!     1.2 other model levels twostream

  do k=1,kmax
    kp=k+1
    km=k-1

    if (whls(kp).lt.0) then   !downwind scheme for subsidence
       thlp(2:i1,2:j1,k) = thlp(2:i1,2:j1,k) - whls(kp) * (thl0(2:i1,2:j1,kp) - thl0(2:i1,2:j1,k))/dzh(kp)
       qtp (2:i1,2:j1,k) = qtp (2:i1,2:j1,k) - whls(kp) * (qt0 (2:i1,2:j1,kp) - qt0 (2:i1,2:j1,k))/dzh(kp)
       if (lmomsubs) then
          up(2:i1,2:j1,k) = up(2:i1,2:j1,k) - whls(kp) * (u0(2:i1,2:j1,kp) - u0(2:i1,2:j1,k))/dzh(kp)
          vp(2:i1,2:j1,k) = vp(2:i1,2:j1,k) - whls(kp) * (v0(2:i1,2:j1,kp) - v0(2:i1,2:j1,k))/dzh(kp)
       endif

       svp(2:i1,2:j1,k,:) = svp(2:i1,2:j1,k,:) - whls(kp) * (sv0(2:i1,2:j1,kp,:) - sv0(2:i1,2:j1,k,:))/dzh(kp)

    else !downwind scheme for mean upward motions
       if (k > 1) then !neglect effect of mean ascending on tendencies at the lowest full level
          thlp(2:i1,2:j1,k) = thlp(2:i1,2:j1,k) - whls(k) * (thl0(2:i1,2:j1,k) - thl0(2:i1,2:j1,km))/dzh(k)
          qtp (2:i1,2:j1,k) = qtp (2:i1,2:j1,k) - whls(k) * (qt0 (2:i1,2:j1,k) - qt0 (2:i1,2:j1,km))/dzh(k)
          if (lmomsubs) then
             up(2:i1,2:j1,k) = up(2:i1,2:j1,k) - whls(k) * (u0(2:i1,2:j1,k) - u0(2:i1,2:j1,km))/dzh(k)
             vp(2:i1,2:j1,k) = vp(2:i1,2:j1,k) - whls(k) * (v0(2:i1,2:j1,k) - v0(2:i1,2:j1,km))/dzh(k)
          endif
          svp(2:i1,2:j1,k,:) = svp(2:i1,2:j1,k,:)-whls(k) * (sv0(2:i1,2:j1,k,:) - sv0(2:i1,2:j1,km,:))/dzh(k)
       endif
    endif

    thlp(2:i1,2:j1,k) = thlp(2:i1,2:j1,k)-u0av(k)*dthldxls(k)-v0av(k)*dthldyls(k) + dthldtls(k)
    qtp (2:i1,2:j1,k) = qtp (2:i1,2:j1,k)-u0av(k)*dqtdxls (k)-v0av(k)*dqtdyls (k) + dqtdtls(k)
    up  (2:i1,2:j1,k) = up  (2:i1,2:j1,k)-u0av(k)*dudxls  (k)-v0av(k)*dudyls  (k) + dudtls(k)
    vp  (2:i1,2:j1,k) = vp  (2:i1,2:j1,k)-u0av(k)*dvdxls  (k)-v0av(k)*dvdyls  (k) + dvdtls(k)
    
  enddo

  return
  end subroutine lstend

end module modforces
