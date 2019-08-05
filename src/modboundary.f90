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
  real :: rnu0 = 2.75e-3
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
! Sponge layer
    if (ksp==-1) then
      ksp  = min(3*kmax/4,kmax - 15)
    end if

    zspb    = zf(ksp)
    zspt    = zf(kmax)

    tsc(1:ksp-1) = 0.0
    do k=ksp,kmax
      tsc(k) = rnu0*sin(0.5*pi*(zf(k)-zspb)/(zspt-zspb))**2
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
  implicit none

    call cyclicm
    call cyclich
    call topm
    call toph
  end subroutine boundary
!> Cleans up after the run
  subroutine exitboundary
  implicit none
    deallocate(tsc)
  end subroutine exitboundary

!> Sets lateral periodic boundary conditions for the scalars
 subroutine cyclich

  use modglobal, only : i1,ih,j1,jh,k1,nsv
  use modfields, only : thl0,thlm,qt0,qtm,sv0,svm
  use modmpi,    only : excjs

  integer n

  call excjs( thl0           , 2,i1,2,j1,1,k1,ih,jh,.true.)
  call excjs( qt0            , 2,i1,2,j1,1,k1,ih,jh,.true.)
  call excjs( thlm           , 2,i1,2,j1,1,k1,ih,jh,.true.)
  call excjs( qtm            , 2,i1,2,j1,1,k1,ih,jh,.true.)

  do n=1,nsv
    call excjs( sv0(:,:,:,n)   , 2,i1,2,j1,1,k1,ih,jh,.true.)
    call excjs( svm(:,:,:,n)   , 2,i1,2,j1,1,k1,ih,jh,.true.)
  enddo

  return
  end subroutine cyclich

!>set lateral periodic boundary conditions for momentum
 subroutine cyclicm

  use modglobal, only : i1,ih,j1,jh,k1
  use modfields, only : u0,um,v0,vm,w0,wm,e120,e12m
  use modmpi,    only : excjs

  call excjs( u0  , 2,i1,2,j1,1,k1,ih,jh,.true.)
  call excjs( v0  , 2,i1,2,j1,1,k1,ih,jh,.true.)
  call excjs( w0  , 2,i1,2,j1,1,k1,ih,jh,.true.)
  call excjs( um  , 2,i1,2,j1,1,k1,ih,jh,.true.)
  call excjs( vm  , 2,i1,2,j1,1,k1,ih,jh,.true.)
  call excjs( wm  , 2,i1,2,j1,1,k1,ih,jh,.true.)
  call excjs( e120  , 2,i1,2,j1,1,k1,ih,jh,.true.)
  call excjs( e12m  , 2,i1,2,j1,1,k1,ih,jh,.true.)

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
  use modglobal, only : i1,j1,kmax,cu,cv,lcoriol,igrw_damp,geodamptime,nsv,rdt,unudge,dzf
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
  case(-1)
    up(:,:,:) = up(:,:,:) - unudge * ( sum((u0av(1:kmax) - ug(1:kmax)) * dzf(1:kmax)) / sum(dzf(1:kmax)) ) / rdt
    vp(:,:,:) = vp(:,:,:) - unudge * ( sum((v0av(1:kmax) - vg(1:kmax)) * dzf(1:kmax)) / sum(dzf(1:kmax)) ) / rdt
  case default
    stop "no gravity wave damping option selected"
  end select

  ! Additional to gravity wave damping, set qt, thl and sv0(:) equal to slabaverage
  ! at level kmax.
  ! Originally done in subroutine tqaver, now using averages from modthermodynamics

  thl0(2:i1,2:j1,kmax) = thl0av(kmax)
  qt0 (2:i1,2:j1,kmax) = qt0av(kmax)
  do n=1,nsv
    sv0(2:i1,2:j1,kmax,n) = sv0av(kmax,n)
  end do

  return
  end subroutine grwdamp

!> Sets top boundary conditions for scalars
  subroutine toph

  use modglobal, only : kmax,k1,nsv,dtheta,dqt,dsv,dzh
  use modfields, only : thl0,thlm,qt0,qtm,sv0,svm &
                       ,thl0av,qt0av,sv0av
  implicit none
  integer :: n
  integer,parameter :: kav=5

! **  Top conditions :
  ! Calculate new gradient over several of the top levels, to be used
  ! to extrapolate thl and qt to level k1 !JvdD
  dtheta = sum((thl0av(kmax-kav+1:kmax)-thl0av(kmax-kav:kmax-1))/ &
             dzh(kmax-kav+1:kmax))/kav
  dqt    = sum((qt0av (kmax-kav+1:kmax)-qt0av (kmax-kav:kmax-1))/ &
             dzh(kmax-kav+1:kmax))/kav
  do n=1,nsv
    dsv(n) = sum((sv0av(kmax-kav+1:kmax,n)-sv0av(kmax-kav:kmax-1,n))/ &
               dzh(kmax-kav:kmax-1))/kav
  enddo

  thl0(:,:,k1) = thl0(:,:,kmax) + dtheta*dzh(k1)
  qt0(:,:,k1)  = qt0 (:,:,kmax) + dqt*dzh(k1)

  thlm(:,:,k1) = thlm(:,:,kmax) + dtheta*dzh(k1)
  qtm(:,:,k1)  = qtm (:,:,kmax) + dqt*dzh(k1)
  do n=1,nsv
    sv0(:,:,k1,n) = sv0(:,:,kmax,n) + dsv(n)*dzh(k1)
    svm(:,:,k1,n) = svm(:,:,kmax,n) + dsv(n)*dzh(k1)
  enddo

  return
  end subroutine toph
!> Sets top boundary conditions for momentum
  subroutine topm

    use modglobal, only : kmax,k1,e12min,lrigidlid
    use modfields, only : u0,v0,w0,e120,um,vm,wm,e12m
    implicit none
    u0(:,:,k1)   = u0(:,:,kmax)
    v0(:,:,k1)   = v0(:,:,kmax)
    w0(:,:,k1)   = 0.0
    e120(:,:,k1) = e12min
    if (lrigidlid) e120(:,:,k1) = e120(:,:,kmax)

    um(:,:,k1)   = um(:,:,kmax)
    vm(:,:,k1)   = vm(:,:,kmax)
    wm(:,:,k1)   = 0.0
    e12m(:,:,k1) = e12min
    if (lrigidlid) e12m(:,:,k1) = e12m(:,:,kmax)

  return
  end subroutine topm

!!>Set thl, qt and sv(n) equal to slab average at level kmax
! Functionality added to subroutine 'grwdamp' !JvdD
!
!  subroutine tqaver
!
!  use modmpi,    only : comm3d,mpierr,my_real, mpi_sum
!  use modglobal, only : i1,j1,kmax,nsv,ijtot
!  use modfields, only : thl0,qt0,sv0
!  implicit none
!
!  real thl0a, qt0a
!  real thl0al, qt0al
!  integer n
!  real,allocatable, dimension(:) :: sv0al, sv0a
!  allocate (sv0al(nsv),sv0a(nsv))
!
!  thl0al=sum(thl0(2:i1,2:j1,kmax))
!  qt0al =sum(qt0(2:i1,2:j1,kmax))
!
!  do n=1,nsv
!    sv0al(n) = sum(sv0(2:i1,2:j1,kmax,n))
!  enddo
!
!  call MPI_ALLREDUCE(thl0al, thl0a, 1,    MY_REAL, &
!                         MPI_SUM, comm3d,mpierr)
!  call MPI_ALLREDUCE(qt0al, qt0a , 1,     MY_REAL, &
!                         MPI_SUM, comm3d,mpierr)
!  if(nsv > 0) then
!    call MPI_ALLREDUCE(sv0al, sv0a , nsv,   MY_REAL, &
!                           MPI_SUM, comm3d,mpierr)
!  end if
!
!
!  thl0a=thl0a/ijtot
!  qt0a =qt0a/ijtot
!  sv0a = sv0a/ijtot
!
!  thl0(2:i1,2:j1,kmax)=thl0a
!  qt0(2:i1,2:j1,kmax) =qt0a
!  do n=1,nsv
!    sv0(2:i1,2:j1,kmax,n) = sv0a(n)
!  enddo
!  deallocate (sv0al,sv0a)
!
!  return
!  end subroutine tqaver



end module
