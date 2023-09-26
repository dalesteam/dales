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
  use modmicrodata, only : imicro, imicro_bulk, imicro_bin, imicro_sice, imicro_sice2, iqr
  implicit none

  integer i, j, k

  if (lforce_user) call force_user

  if (lpressgrad) then
     !$acc kernels default(present) async(1)
     do k=1,kmax
        up(:,:,k) = up(:,:,k) - dpdxl(k)      !RN LS pressure gradient force in x,y directions;
        vp(:,:,k) = vp(:,:,k) - dpdyl(k)
     end do
     !$acc end kernels
  end if

  if((imicro==imicro_sice).or.(imicro==imicro_sice2).or.(imicro==imicro_bulk).or.(imicro==imicro_bin)) then
    !$acc kernels default(present) async(2)
    do k=2,kmax
       wp(:,:,k) = wp(:,:,k) + grav*(thv0h(:,:,k)-thvh(k))/thvh(k) - &
                  grav*(sv0(:,:,k,iqr)*dzf(k-1)+sv0(:,:,k-1,iqr)*dzf(k))/(2.0*dzh(k))
    end do
    !$acc end kernels
  else
    !$acc kernels default(present) async(2)
    do k=2,kmax
      wp(:,:,k) = wp(:,:,k) + grav*(thv0h(:,:,k)-thvh(k))/thvh(k)
    end do
    !$acc end kernels
  end if

!     --------------------------------------------
!     special treatment for lowest full level: k=1
!     --------------------------------------------
  !$acc kernels default(present) async(3)
  wp(:,:,1) = 0.0
  !$acc end kernels

  !$acc wait

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

  integer i, j, k

  if (lcoriol .eqv. .false.) return

  !$acc parallel loop collapse(3) default(present) async
  do k=2,kmax
    do j=2,j1
      do i=2,i1

        up(i,j,k) = up(i,j,k)+ cv*om23 &
              +(v0(i,j,k)+v0(i,j+1,k)+v0(i-1,j,k)+v0(i-1,j+1,k))*om23*0.25 &
              -(w0(i,j,k)+w0(i,j,k+1)+w0(i-1,j,k+1)+w0(i-1,j,k))*om22*0.25

        vp(i,j,k) = vp(i,j,k)  - cu*om23 &
              -(u0(i,j,k)+u0(i,j-1,k)+u0(i+1,j-1,k)+u0(i+1,j,k))*om23*0.25


        wp(i,j,k) = wp(i,j,k) + cu*om22 +( (dzf(k-1) * (u0(i,j,k)  + u0(i+1,j,k) )    &
                    +    dzf(k)  * (u0(i,j,k-1) + u0(i+1,j,k-1))  ) / dzh(k) ) &
                    * om22*0.25
      end do
    end do
!     -------------------------------------------end i&j-loop
  end do
!     -------------------------------------------end k-loop

!     --------------------------------------------
!     special treatment for lowest full level: k=1
!     --------------------------------------------
  !$acc parallel loop collapse(2) default(present) async
  do j=2,j1
    do i=2,i1

      up(i,j,1) = up(i,j,1)  + cv*om23 &
            +(v0(i,j,1)+v0(i,j+1,1)+v0(i-1,j,1)+v0(i-1,j+1,1))*om23*0.25 &
            -(w0(i,j,1)+w0(i,j ,2)+w0(i-1,j,2)+w0(i-1,j ,1))*om22*0.25

      vp(i,j,1) = vp(i,j,1) - cu*om23 &
            -(u0(i,j,1)+u0(i,j-1,1)+u0(i+1,j-1,1)+u0(i+1,j,1))*om23*0.25

      wp(i,j,1) = 0.0

    end do
  end do
!     ----------------------------------------------end i,j-loop
  !$acc wait
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

  integer i, j, k, n

!     1. DETERMINE LARGE SCALE TENDENCIES
!        --------------------------------

!     1.1 lowest model level above surface : only downward component
!     1.2 other model levels twostream

  !$acc parallel loop gang default(present) async(1)
  do k=1,kmax
    if (whls(k+1) < 0) then   !downwind scheme for subsidence
      !$acc loop collapse(2)
      do j=2,j1
        do i=2,i1
          thlp(i,j,k) = thlp(i,j,k) - whls(k+1) * (thl0(i,j,k+1) - thl0(i,j,k))/dzh(k+1)
          qtp (i,j,k) = qtp (i,j,k) - whls(k+1) * (qt0 (i,j,k+1) - qt0 (i,j,k))/dzh(k+1)
        end do
      end do
    else !downwind scheme for mean upward motions
      if (k > 1) then !neglect effect of mean ascending on tendencies at the lowest full level
        !$acc loop collapse(2)
        do j=2,j1
          do i=2,i1
            thlp(i,j,k) = thlp(i,j,k) - whls(k) * (thl0(i,j,k) - thl0(i,j,k-1))/dzh(k)
            qtp (i,j,k) = qtp (i,j,k) - whls(k) * (qt0 (i,j,k) - qt0 (i,j,k-1))/dzh(k)
          end do
        end do
      endif
    endif
  end do

  if (lmomsubs) then
    !$acc parallel loop gang default(present) async(2)
    do k=1,kmax
      if (whls(k+1) < 0) then
        !$acc loop collapse(2)
        do j=1,j1
          do i=1,i1
            up(i,j,k) = up(i,j,k) - whls(k+1) * (u0(i,j,k+1) - u0(i,j,k))/dzh(k+1)
            vp(i,j,k) = vp(i,j,k) - whls(k+1) * (v0(i,j,k+1) - v0(i,j,k))/dzh(k+1)
          end do
        end do
      else
        if (k > 1) then
          !$acc loop collapse(2)
          do j=1,j1
            do i=1,i1 
              up(i,j,k) = up(i,j,k) - whls(k) * (u0(i,j,k) - u0(i,j,k-1))/dzh(k)
              vp(i,j,k) = vp(i,j,k) - whls(k) * (v0(i,j,k) - v0(i,j,k-1))/dzh(k)
            end do
          end do
        end if
      end if
    end do
  end if 
  
  !$acc parallel loop collapse(3) async wait(1,2)
  do k=1,kmax
    do j=2,j1
      do i=2,i1
        thlp(i,j,k) = thlp(i,j,k)-u0av(k)*dthldxls(k)-v0av(k)*dthldyls(k) + dthldtls(k)
        qtp (i,j,k) = qtp (i,j,k)-u0av(k)*dqtdxls (k)-v0av(k)*dqtdyls (k) + dqtdtls(k)
        up  (i,j,k) = up  (i,j,k)-u0av(k)*dudxls  (k)-v0av(k)*dudyls  (k) + dudtls(k)
        vp  (i,j,k) = vp  (i,j,k)-u0av(k)*dvdxls  (k)-v0av(k)*dvdyls  (k) + dvdtls(k)
      end do
    end do
  end do
  
  ! Only do above for scalars if there are any scalar fields
  if (nsv > 0) then
    do n=1,nsv
      !$acc parallel loop gang default(present) async
      do k=1,kmax
        if (whls(k+1).lt.0) then
          !$acc loop collapse(2)
          do j=1,j1
            do i=1,i1
              svp(i,j,k,n) = svp(i,j,k,n) - whls(k+1) * (sv0(i,j,k+1,n) - sv0(i,j,k,n))/dzh(k+1)
            end do
          end do
        else
          if (k > 1) then
            !$acc loop collapse(2)
            do j=1,j1
              do i=1,i1
                svp(i,j,k,n) = svp(i,j,k,n)-whls(k) * (sv0(i,j,k,n) - sv0(i,j,k-1,n))/dzh(k)
              end do
            end do
          end if
        end if
      end do
    end do
  end if
  
  !$acc wait

  return
  end subroutine lstend

end module modforces
