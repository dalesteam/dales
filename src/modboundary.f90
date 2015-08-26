!> \file modboundary.f90
!!  Takes care of all the boundaries, except for the surface
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

!>
!! Modboundary takes care of all the boundaries, except for the surface
!>
!! This module takes care of the periodic boundary conditions in x- and y, of
!! the top inlet conditions, of the gravity wave damping.
!! \par Revision list
!! \par Authors
!!
module modboundary


implicit none
save
private
public :: initboundary, boundary, exitboundary,grwdamp, ksp,cyclich
  integer :: ksp = -1                 !<    lowest level of sponge layer
  real,allocatable :: tsc(:)          !<   damping coefficients to be used in grwdamp.
  real,allocatable :: thl_nudge (:), &  !<   nudging profiles for the different variables
                      qt_nudge  (:), &  !    used for igrw_damp=(4)
                      u_nudge   (:), &
                      v_nudge   (:)
  real :: rnu0 = 1./3600.
  logical :: isInitSponge = .false.   !<   Switch for initialization of nudge profiles in
                                      !    the sponge layer
  logical :: isInitExtra  = .false.

contains
!>
!! Initializing Boundary; specifically the sponge layer
!>
  subroutine initboundary
    use modglobal, only : k1,kmax,pi,zf
    implicit none

    real    :: zspb, zspt
    integer :: k
    allocate(tsc(k1))
    allocate(thl_nudge(k1),qt_nudge(k1),u_nudge(k1),v_nudge(k1))
! Sponge layer
    if (ksp==-1) then
      ksp  = min(3*kmax/4,kmax - 15)
    end if

    zspb    = zf(ksp)
    zspt    = zf(kmax)

    tsc(1:ksp-1) = 0.0
    do k=ksp,kmax
      tsc(k) = rnu0*sin(0.5*pi*(zf(k)-zspb)/(zspt-zspb))**2
!      tsc(k)  = 0.5*rnu0*(1-cos(pi*(zf(k)-zspb)/(zspt-zspb))) ! CGILS nudging toward initial conditions
    end do
    tsc(k1)=tsc(kmax)
  end subroutine initboundary

!>
!! Execute boundary conditions
!>
!! The boundary conditions at the top of the domain and in the horizontal
!! directions are relatively straightforward. In the horizontal directions,
!! periodic boundary conditions are applied.
!! \latexonly
!! At the top of the domain, we take:
!! \begin{equation}
!!  \derr{\fav{u}}{z} = \derr{\fav{v}}{z} = 0;\\ \fav{w} = 0;\\
!! \derr{\fav{\varphi}}{z} =  \mr{cst}.
!! \end{equation}
!! \endlatexonly
  subroutine boundary
    use modthermodynamics, only : avgProfs
  implicit none

    call avgProfs ! Calculates domain averaged vertical profiles of prognostic variables
    call cyclicm
    call cyclich
    call topm
    call toph
    call avgProfs ! Calculates domain averaged vertical profiles of prognostic variables
  end subroutine boundary
!> Cleans up after the run
  subroutine exitboundary
  implicit none
    deallocate(tsc)
  end subroutine exitboundary

!> Sets lateral periodic boundary conditions for the scalars
 subroutine cyclich

  use modglobal, only : i1,i2,ih,j1,jh,k1,nsv
  use modfields, only : thl0,qt0,sv0
  use modmpi,    only : excjs
  integer n,m

  do m = 1,ih
    thl0(2-m,:,:)   = thl0(i2-m,:,:)
    thl0(i1+m,:,:)  = thl0(1+m,:,:)
    qt0(2-m,:,:)    = qt0(i2-m,:,:)
    qt0(i1+m,:,:)   = qt0(1+m,:,:)
    sv0(2-m,:,:,:)  = sv0(i2-m,:,:,:)
    sv0(i1+m,:,:,:) = sv0(1+m,:,:,:)
  end do

  call excjs( thl0           , 2,i1,2,j1,1,k1,ih,jh)
  call excjs( qt0            , 2,i1,2,j1,1,k1,ih,jh)

  do n=1,nsv
    call excjs( sv0(:,:,:,n)   , 2,i1,2,j1,1,k1,ih,jh)
  enddo

  return
  end subroutine cyclich

!>set lateral periodic boundary conditions for momentum
 subroutine cyclicm

  use modglobal, only : i1,i2,ih,j1,jh,k1
  use modfields, only : u0,v0,w0,e120
  use modmpi,    only : excjs

  integer m

  do m = 1,ih

    u0(2-m,:,:)    = u0(i2-m,:,:)
    u0(i1+m,:,:)   = u0(1+m,:,:)
    v0(2-m,:,:)    = v0(i2-m,:,:)
    v0(i1+m,:,:)   = v0(1+m,:,:)
    w0(2-m,:,:)    = w0(i2-m,:,:)
    w0(i1+m,:,:)   = w0(1+m,:,:)

    e120(2-m,:,:)  = e120(i2-m,:,:)
    e120(i1+m,:,:) = e120(1+m,:,:)

  end do

  call excjs( u0  , 2,i1,2,j1,1,k1,ih,jh)
  call excjs( v0  , 2,i1,2,j1,1,k1,ih,jh)
  call excjs( w0  , 2,i1,2,j1,1,k1,ih,jh)
  call excjs( e120, 2,i1,2,j1,1,k1,ih,jh)

  return
  end subroutine cyclicm

!>
!! grwdamp damps gravity waves in the upper part of the domain.
!>
!! The lower limit of the damping region is set by ksp
!! Horizontal fluctuations at the top of the domain (for instance gravity waves)
!! are damped out by a sponge layer through an additional forcing/source term.
!! \latexonly
!! \begin{eqnarray}
!! \force{i}{sp}(z) &=& -\frac{1}{t^{\mr{sp}}}\left(\xav{\fav{u_i}}-\fav{u_i}\right), \\\\
!!  \source{\varphi}{sp}(z) &=& -\frac{1}{t^{\mr{sp}}}\left(\xav{\varphi}-\varphi\right),
!! \end{eqnarray}
!! with $t^{\mr{sp}}$ a relaxation time scale that goes from
!! $t^{\mr{sp}}_0=1/(2.75\times10^{-3})\mr{s}\approx 6$min at the top of the domain
!! to infinity at the bottom of the sponge layer.
!! \endlatexonly
 subroutine grwdamp
  use modglobal, only : i1,j1,kmax,cu,cv,lcoriol,igrw_damp,geodamptime,nsv
  use modfields, only : up,vp,wp,thlp,qtp,u0,v0,w0,thl0,qt0,sv0,ug,vg & 
                        ,thl0av,qt0av,sv0av,u0av,v0av
  implicit none

  integer k,n

  select case(igrw_damp)
  case(0) !do nothing
  case(1)
    do k=ksp,kmax
      up(:,:,k)  = up(:,:,k)-(u0(:,:,k)-(u0av(k)-cu))*tsc(k)
      vp(:,:,k)  = vp(:,:,k)-(v0(:,:,k)-(v0av(k)-cv))*tsc(k)
      wp(:,:,k)  = wp(:,:,k)-w0(:,:,k)*tsc(k)
      thlp(:,:,k)= thlp(:,:,k)-(thl0(:,:,k)-thl0av(k))*tsc(k)
      qtp(:,:,k) = qtp(:,:,k)-(qt0(:,:,k)-qt0av(k))*tsc(k)
    end do
    if(lcoriol) then
    do k=ksp,kmax
      up(:,:,k)  = up(:,:,k)-(u0(:,:,k)-(ug(k)-cu))*((1./(geodamptime*rnu0))*tsc(k))
      vp(:,:,k)  = vp(:,:,k)-(v0(:,:,k)-(vg(k)-cv))*((1./(geodamptime*rnu0))*tsc(k))
    end do
    end if
  case(2)
    do k=ksp,kmax
      up(:,:,k)  = up(:,:,k)-(u0(:,:,k)-(ug(k)-cu))*tsc(k)
      vp(:,:,k)  = vp(:,:,k)-(v0(:,:,k)-(vg(k)-cv))*tsc(k)
      wp(:,:,k)  = wp(:,:,k)-w0(:,:,k)*tsc(k)
      thlp(:,:,k)= thlp(:,:,k)-(thl0(:,:,k)-thl0av(k))*tsc(k)
      qtp(:,:,k) = qtp(:,:,k)-(qt0(:,:,k)-qt0av(k))*tsc(k)
    end do
  case(3)
    do k=ksp,kmax
      up(:,:,k)  = up(:,:,k)-(u0(:,:,k)-(u0av(k)-cu))*tsc(k)
      vp(:,:,k)  = vp(:,:,k)-(v0(:,:,k)-(v0av(k)-cv))*tsc(k)
      wp(:,:,k)  = wp(:,:,k)-w0(:,:,k)*tsc(k)
      thlp(:,:,k)= thlp(:,:,k)-(thl0(:,:,k)-thl0av(k))*tsc(k)
      qtp(:,:,k) = qtp(:,:,k)-(qt0(:,:,k)-qt0av(k))*tsc(k)
    end do
  case(4)
    ! Option for nudging thl, qt, u and v to their initial profiles
    if (.not.isInitSponge) then
      thl_nudge(ksp:kmax) = thl0av(ksp:kmax)
      qt_nudge (ksp:kmax) = qt0av (ksp:kmax)
      u_nudge  (ksp:kmax) = u0av  (ksp:kmax)
      v_nudge  (ksp:kmax) = v0av  (ksp:kmax)
      isInitSponge = .true.
    end if

    do k=ksp,kmax
      up(:,:,k)  = up(:,:,k)  - (u0(:,:,k)-(u_nudge(k)-cu))*tsc(k)
      vp(:,:,k)  = vp(:,:,k)  - (v0(:,:,k)-(v_nudge(k)-cv))*tsc(k)
      wp(:,:,k)  = wp(:,:,k)  -  w0(:,:,k)*tsc(k)
      thlp(:,:,k)= thlp(:,:,k)- (thl0(:,:,k)-thl_nudge(k))*tsc(k)
      qtp(:,:,k) = qtp(:,:,k) - (qt0(:,:,k) -qt_nudge (k))*tsc(k)
    end do
  case default
    stop "no gravity wave damping option selected"
  end select

  
  ! Additional to gravity wave damping, set qt, thl and sv0(:) equal to slabaverage
  ! at level kmax.
  ! Originally done in subroutine tqaver, now using averages from modthermodynamics
  ! These lines make the new time integration scheme crash. Furthermore, it is not
  ! good practice to manually set boundary values inside the domain.
  ! Note that k1 values are set through extrapolation in subroutine toph

  !thl0(:,:,kmax) = thl0av(kmax)
  !qt0 (:,:,kmax) = qt0av(kmax)
  !do n=1,nsv
  !  sv0(:,:,kmax,n) = sv0av(kmax,n)
  !end do

  return
  end subroutine grwdamp

!> Sets top boundary conditions for scalars
  subroutine toph

    use modglobal, only : i1,j1,kmax,k1,nsv,dtheta,dqt,dsv,dzh
    use modfields, only : thl0,qt0,sv0,thl0av,qt0av,sv0av
    implicit none
    integer :: n,i,j
    integer,parameter :: kav=5

    if (.not. isInitExtra) then
      ! **  Top conditions :
      ! Calculate new gradient over several of the top levels, to be used
      ! to extrapolate thl and qt to level k1 !JvdD
      !dtheta = sum((thl0av(kmax-kav+1:kmax)-thl0av(kmax-kav:kmax-1))/ &
      !           dzh(kmax-kav+1:kmax))/kav
      !dqt    = sum((qt0av (kmax-kav+1:kmax)-qt0av (kmax-kav:kmax-1))/ &
      !           dzh(kmax-kav+1:kmax))/kav
      dtheta = (thl0av(kmax)-thl0av(kmax-1))/dzh(kmax)
      dqt    = (qt0av (kmax)-qt0av (kmax-1))/dzh(kmax)
      do n=1,nsv
        !dsv(n) = sum((sv0av(kmax-kav+1:kmax,n)-sv0av(kmax-kav:kmax-1,n))/ &
        !           dzh(kmax-kav:kmax-1))/kav
        dsv(n) = (sv0av(kmax,n)-sv0av(kmax-1,n))/ &
                   dzh(kmax)
      enddo
      isInitExtra = .true.
    end if

!    thl0(2:i1,2:j1,k1) = thl0(2:i1,2:j1,kmax) + dtheta*dzh(k1)
!    qt0(2:i1,2:j1,k1)  = qt0 (2:i1,2:j1,kmax) + dqt*dzh(k1)
!    do n=1,nsv
!      sv0(2:i1,2:j1,k1,n) = sv0(2:i1,2:j1,kmax,n) + dsv(n)*dzh(k1)
!    enddo
!    do i=2,i1
!      do j=2,j1
        thl0(:,:,k1) = thl0(:,:,kmax) + dtheta*dzh(k1)
        qt0(:,:,k1)  = qt0 (:,:,kmax) + dqt*dzh(k1)
        do n=1,nsv
          sv0(:,:,k1,n) = sv0(:,:,kmax,n) + dsv(n)*dzh(k1)
        enddo
!      end do
!    end do

    return
  end subroutine toph
!> Sets top boundary conditions for momentum
  subroutine topm

    use modglobal, only : i1,j1,kmax,k1,e12min
    use modfields, only : u0,v0,w0,e120
    implicit none
    u0(2:i1,2:j1,k1)   = u0(2:i1,2:j1,kmax)
    v0(2:i1,2:j1,k1)   = v0(2:i1,2:j1,kmax)
    w0(2:i1,2:j1,k1)   = 0.
    e120(2:i1,2:j1,k1) = e12min

  return
  end subroutine topm

end module
