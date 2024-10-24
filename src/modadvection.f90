!> \file advection.f90
!!  Advection management

!>
!!  Advection management
!! \par Revision list
!! Thijs Heus, Chiel van Heerwaarden, 15 June 2007
!! \par Authors
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

module modadvection
use modtimer
contains
!> Advection redirection function
subroutine advection

  use modglobal,      only : lmoist, nsv, iadv_mom,iadv_tke,iadv_thl,iadv_qt,iadv_sv, &
                             iadv_cd2,iadv_5th,iadv_52,iadv_cd6,iadv_62,iadv_kappa,&
                             iadv_upw,iadv_hybrid,iadv_hybrid_f,iadv_null,leq,&
                             lopenbc,lboundary,lperiodic,i1,j1
  use modfields,      only : u0,up,v0,vp,w0,wp,e120,e12p,thl0,thlp,qt0,qtp,sv0,svp
  use modsubgrid,     only : lsmagorinsky
  use modsamptend,    only : samptend, tend_hadv,tend_vadv
  use advec_2nd,      only : hadvecu_2nd, vadvecu_2nd, hadvecv_2nd, vadvecv_2nd, &
                             hadvecw_2nd, vadvecw_2nd, hadvecc_2nd, vadvecc_2nd
  use advec_5th,      only : hadvecu_5th, vadvecu_5th, hadvecv_5th, vadvecv_5th, &
                             hadvecw_5th, vadvecw_5th, hadvecc_5th, vadvecc_5th
  use advec_6th,      only : hadvecu_6th, vadvecu_6th, hadvecv_6th, vadvecv_6th, &
                             hadvecw_6th, vadvecw_6th, hadvecc_6th, vadvecc_6th
  use advec_hybrid,   only : hadvecc_hybrid,vadvecc_hybrid
  use advec_hybrid_f, only : hadvecc_hybrid_f,vadvecc_hybrid_f
  use advec_kappa,    only : hadvecc_kappa, vadvecc_kappa
  use advec_upw,      only : hadvecc_upw, vadvecc_upw
  use modopenboundary,only : advecc_2nd_boundary_buffer,advecu_2nd_boundary_buffer,advecv_2nd_boundary_buffer,advecw_2nd_boundary_buffer
  implicit none
  integer :: n,istart,iend,jstart,jend,ibuffer,jbuffer

  call timer_tic('modadvection/advection', 0)

  ! leq = .false. ! for testing that the non-uniform advection routines agree with the uniform ones
                  ! when the grid is uniform

  ! Horizontal advection
  istart = 2; iend = i1; jstart = 2; jend = j1
  select case(iadv_mom)
    case(iadv_cd2)
      call hadvecu_2nd(u0,up,merge(istart+1, istart, (lboundary(1) .and. .not. lperiodic(1))),iend,jstart,jend)
      call hadvecv_2nd(v0,vp,istart,iend,merge(jstart+1, jstart, (lboundary(3) .and. .not. lperiodic(3))),jend)
      call hadvecw_2nd(w0,wp,istart,iend,jstart,jend)
      !$acc wait(1,2,3)
    case(iadv_5th)
      !if (.not. leq) stop "advec_5th does not support a non-uniform vertical grid."
      ibuffer = 2; jbuffer = 2
      istart = 2; iend = i1; jstart = 2; jend = j1
      if(lopenbc) call advecu_2nd_boundary_buffer(u0,up,istart,iend,jstart,jend,ibuffer,jbuffer)
      call hadvecu_5th(u0,up,istart,iend,jstart,jend)
      istart = 2; iend = i1; jstart = 2; jend = j1
      if(lopenbc) call advecv_2nd_boundary_buffer(v0,vp,istart,iend,jstart,jend,ibuffer,jbuffer)
      call hadvecv_5th(v0,vp,istart,iend,jstart,jend)
      istart = 2; iend = i1; jstart = 2; jend = j1
      if(lopenbc) call advecw_2nd_boundary_buffer(w0,wp,istart,iend,jstart,jend,ibuffer,jbuffer)
      call hadvecw_5th(w0,wp,istart,iend,jstart,jend)
      !$acc wait(1,2,3)
    case(iadv_52)
      ibuffer = 2; jbuffer = 2
      istart = 2; iend = i1; jstart = 2; jend = j1
      if(lopenbc) call advecu_2nd_boundary_buffer(u0,up,istart,iend,jstart,jend,ibuffer,jbuffer)
      call hadvecu_5th(u0,up,istart,iend,jstart,jend)
      istart = 2; iend = i1; jstart = 2; jend = j1
      if(lopenbc) call advecv_2nd_boundary_buffer(v0,vp,istart,iend,jstart,jend,ibuffer,jbuffer)
      call hadvecv_5th(v0,vp,istart,iend,jstart,jend)
      istart = 2; iend = i1; jstart = 2; jend = j1
      if(lopenbc) call advecw_2nd_boundary_buffer(w0,wp,istart,iend,jstart,jend,ibuffer,jbuffer)
      call hadvecw_5th(w0,wp,istart,iend,jstart,jend)
      !$acc wait(1,2,3)
    case(iadv_cd6)
      !if (.not. leq) stop "advec_6th does not support a non-uniform vertical grid."
      ibuffer = 2; jbuffer = 2
      istart = 2; iend = i1; jstart = 2; jend = j1
      if(lopenbc) call advecu_2nd_boundary_buffer(u0,up,istart,iend,jstart,jend,ibuffer,jbuffer)
      call hadvecu_6th(u0,up,istart,iend,jstart,jend)
      istart = 2; iend = i1; jstart = 2; jend = j1
      if(lopenbc) call advecv_2nd_boundary_buffer(v0,vp,istart,iend,jstart,jend,ibuffer,jbuffer)
      call hadvecv_6th(v0,vp,istart,iend,jstart,jend)
      istart = 2; iend = i1; jstart = 2; jend = j1
      if(lopenbc) call advecw_2nd_boundary_buffer(w0,wp,istart,iend,jstart,jend,ibuffer,jbuffer)
      call hadvecw_6th(w0,wp,istart,iend,jstart,jend)
      !$acc wait(1,2,3)
    case(iadv_62)
      ibuffer = 2; jbuffer = 2
      istart = 2; iend = i1; jstart = 2; jend = j1
      if(lopenbc) call advecu_2nd_boundary_buffer(u0,up,istart,iend,jstart,jend,ibuffer,jbuffer)
      call hadvecu_6th(u0,up,istart,iend,jstart,jend)
      istart = 2; iend = i1; jstart = 2; jend = j1
      if(lopenbc) call advecv_2nd_boundary_buffer(v0,vp,istart,iend,jstart,jend,ibuffer,jbuffer)
      call hadvecv_6th(v0,vp,istart,iend,jstart,jend)
      istart = 2; iend = i1; jstart = 2; jend = j1
      if(lopenbc) call advecw_2nd_boundary_buffer(w0,wp,istart,iend,jstart,jend,ibuffer,jbuffer)
      call hadvecw_6th(w0,wp,istart,iend,jstart,jend)
      !$acc wait(1,2,3)
    case(iadv_hybrid)
      !if (.not. leq) stop "advec_5th does not support a non-uniform vertical grid."
      call hadvecu_5th(u0,up,istart,iend,jstart,jend)
      call hadvecv_5th(v0,vp,istart,iend,jstart,jend)
      call hadvecw_5th(w0,wp,istart,iend,jstart,jend)
      !$acc wait(1,2,3)
    case(iadv_hybrid_f)
      !if (.not. leq) stop "advec_5th does not support a non-uniform vertical grid."
      call hadvecu_5th(u0,up,istart,iend,jstart,jend)
      call hadvecv_5th(v0,vp,istart,iend,jstart,jend)
      call hadvecw_5th(w0,wp,istart,iend,jstart,jend)
      !$acc wait(1,2,3)
    case(iadv_null)
      ! null advection scheme
      stop "Null advection scheme selected for iadv_mom - probably a bad idea."
    case default
      stop "Unknown advection scheme "
  end select

istart = 2; iend = i1; jstart = 2; jend = j1
  select case(iadv_thl)
    case(iadv_cd2)
      call hadvecc_2nd(thl0,thlp,istart,iend,jstart,jend)
    case(iadv_5th)
      ibuffer = 2; jbuffer = 2
      if(lopenbc) call advecc_2nd_boundary_buffer(thl0,thlp,istart,iend,jstart,jend,ibuffer,jbuffer)
      call hadvecc_5th(thl0,thlp,istart,iend,jstart,jend)
    case(iadv_52)
      ibuffer = 2; jbuffer = 2
      if(lopenbc) call advecc_2nd_boundary_buffer(thl0,thlp,istart,iend,jstart,jend,ibuffer,jbuffer)
      call hadvecc_5th(thl0,thlp,istart,iend,jstart,jend)
    case(iadv_cd6)
      !if (.not. leq) stop "advec_6th does not support a non-uniform vertical grid."
      ibuffer = 2; jbuffer = 2
      if(lopenbc) call advecc_2nd_boundary_buffer(thl0,thlp,istart,iend,jstart,jend,ibuffer,jbuffer)
      call hadvecc_6th(thl0,thlp,istart,iend,jstart,jend)
    case(iadv_62)
      ibuffer = 2; jbuffer = 2
      if(lopenbc) call advecc_2nd_boundary_buffer(thl0,thlp,istart,iend,jstart,jend,ibuffer,jbuffer)
      call hadvecc_6th(thl0,thlp,istart,iend,jstart,jend)
    case(iadv_kappa)
      call hadvecc_kappa(thl0,thlp)
    case(iadv_upw)
      if (.not. leq) stop "advec_upw does not support a non-uniform vertical grid."
      call hadvecc_upw(thl0,thlp)
    case(iadv_hybrid)
       !if (.not. leq) stop "advec_hybrid does not support a non-uniform vertical grid."
      call hadvecc_hybrid(thl0,thlp)
    case(iadv_hybrid_f)
      !if (.not. leq) stop "advec_hybrid_f does not support a non-uniform vertical grid."
      call hadvecc_hybrid_f(thl0,thlp,1.0)
    case(iadv_null)
      ! null advection scheme
      stop "Null advection scheme selected for iadv_thl - probably a bad idea."
    case default
      stop "Unknown advection scheme "
  end select

  if (lmoist) then
    istart = 2; iend = i1; jstart = 2; jend = j1
    select case(iadv_qt)
      case(iadv_cd2)
        call hadvecc_2nd(qt0,qtp,istart,iend,jstart,jend)
      case(iadv_5th)
        ibuffer = 2; jbuffer = 2
        if(lopenbc) call advecc_2nd_boundary_buffer(qt0,qtp,istart,iend,jstart,jend,ibuffer,jbuffer)
        call hadvecc_5th(qt0,qtp,istart,iend,jstart,jend)
      case(iadv_52)
        ibuffer = 2; jbuffer = 2
        if(lopenbc) call advecc_2nd_boundary_buffer(qt0,qtp,istart,iend,jstart,jend,ibuffer,jbuffer)
        call hadvecc_5th(qt0,qtp,istart,iend,jstart,jend)
      case(iadv_cd6)
        !if (.not. leq) stop "advec_6th does not support a non-uniform vertical grid."
        ibuffer = 2; jbuffer = 2
        if(lopenbc) call advecc_2nd_boundary_buffer(qt0,qtp,istart,iend,jstart,jend,ibuffer,jbuffer)
        call hadvecc_6th(qt0,qtp,istart,iend,jstart,jend)
      case(iadv_62)
        ibuffer = 2; jbuffer = 2
        if(lopenbc) call advecc_2nd_boundary_buffer(qt0,qtp,istart,iend,jstart,jend,ibuffer,jbuffer)
        call hadvecc_6th(qt0,qtp,istart,iend,jstart,jend)
      case(iadv_kappa)
        call hadvecc_kappa(qt0,qtp)
      case(iadv_upw)
        if (.not. leq) stop "advec_upw does not support a non-uniform vertical grid."
        call hadvecc_upw(qt0,qtp)
      case(iadv_hybrid)
        !if (.not. leq) stop "advec_hybrid does not support a non-uniform vertical grid."
        call hadvecc_hybrid(qt0,qtp)
      case(iadv_hybrid_f)
        !if (.not. leq) stop "advec_hybrid_f does not support a non-uniform vertical grid."
        call hadvecc_hybrid_f(qt0,qtp,1e-3)
      case(iadv_null)
        ! null advection scheme
        stop "Null advection scheme selected for iadv_qt - probably a bad idea."
      case default
        stop "Unknown advection scheme "
    end select
  end if

  do n=1,nsv
    istart = 2; iend = i1; jstart = 2; jend = j1
    select case(iadv_sv)
    case(iadv_cd2)
      call hadvecc_2nd(sv0(:,:,:,n),svp(:,:,:,n),istart,iend,jstart,jend)
    case(iadv_5th)
      ibuffer = 2; jbuffer = 2
      if(lopenbc) call advecc_2nd_boundary_buffer(sv0(:,:,:,n),svp(:,:,:,n),istart,iend,jstart,jend,ibuffer,jbuffer)
      call hadvecc_5th(sv0(:,:,:,n),svp(:,:,:,n),istart,iend,jstart,jend)
    case(iadv_52)
      ibuffer = 2; jbuffer = 2
      if(lopenbc) call advecc_2nd_boundary_buffer(sv0(:,:,:,n),svp(:,:,:,n),istart,iend,jstart,jend,ibuffer,jbuffer)
      call hadvecc_5th(sv0(:,:,:,n),svp(:,:,:,n),istart,iend,jstart,jend)
    case(iadv_cd6)
      !if (.not. leq) stop "advec_6th does not support a non-uniform vertical grid."
      ibuffer = 2; jbuffer = 2
      if(lopenbc) call advecc_2nd_boundary_buffer(sv0(:,:,:,n),svp(:,:,:,n),istart,iend,jstart,jend,ibuffer,jbuffer)
      call hadvecc_6th(sv0(:,:,:,n),svp(:,:,:,n),istart,iend,jstart,jend)
    case(iadv_62)
      ibuffer = 2; jbuffer = 2
      if(lopenbc) call advecc_2nd_boundary_buffer(sv0(:,:,:,n),svp(:,:,:,n),istart,iend,jstart,jend,ibuffer,jbuffer)
      call hadvecc_6th(sv0(:,:,:,n),svp(:,:,:,n),istart,iend,jstart,jend)
    case(iadv_kappa)
      call hadvecc_kappa(sv0(:,:,:,n),svp(:,:,:,n))
    case(iadv_upw)
      if (.not. leq) stop "advec_upw does not support a non-uniform vertical grid."
      call hadvecc_upw(sv0(:,:,:,n),svp(:,:,:,n))
    case(iadv_hybrid)
      !if (.not. leq) stop "advec_hybrid does not support a non-uniform vertical grid."
      call hadvecc_hybrid(sv0(:,:,:,n),svp(:,:,:,n))
    case(iadv_hybrid_f)
      !if (.not. leq) stop "advec_hybrid_f does not support a non-uniform vertical grid."
      call hadvecc_hybrid_f(sv0(:,:,:,n),svp(:,:,:,n))
    case(iadv_null)
       ! null advection scheme - do nothing
    case default
      stop "Unknown advection scheme "
    end select
  end do
  !$acc wait
  call samptend(tend_hadv)

! Vertical advection
  istart = 2; iend = i1; jstart = 2; jend = j1
  select case(iadv_mom)
    case(iadv_cd2)
      ! Horizontal advection
      call vadvecu_2nd(u0,up,merge(istart+1, istart, (lboundary(1) .and. .not. lperiodic(1))),iend,jstart,jend)
      call vadvecv_2nd(v0,vp,istart,iend,merge(jstart+1, jstart, (lboundary(3) .and. .not. lperiodic(3))),jend)
      call vadvecw_2nd(w0,wp,istart,iend,jstart,jend)
    case(iadv_5th)
!       !if (.not. leq) stop "advec_5th does not support a non-uniform vertical grid."
      call vadvecu_5th(u0,up,merge(istart+1, istart, (lboundary(1) .and. .not. lperiodic(1))),iend,jstart,jend)
      call vadvecv_5th(v0,vp,istart,iend,merge(jstart+1, jstart, (lboundary(3) .and. .not. lperiodic(3))),jend)
      call vadvecw_5th(w0,wp,istart,iend,jstart,jend)
    case(iadv_52)
      call vadvecu_2nd(u0,up,merge(istart+1, istart, (lboundary(1) .and. .not. lperiodic(1))),iend,jstart,jend)
      call vadvecv_2nd(v0,vp,istart,iend,merge(jstart+1, jstart, (lboundary(3) .and. .not. lperiodic(3))),jend)
      call vadvecw_2nd(w0,wp,istart,iend,jstart,jend)
    case(iadv_cd6)
!       !if (.not. leq) stop "advec_6th does not support a non-uniform vertical grid."
      call vadvecu_6th(u0,up,merge(istart+1, istart, (lboundary(1) .and. .not. lperiodic(1))),iend,jstart,jend)
      call vadvecv_6th(v0,vp,istart,iend,merge(jstart+1, jstart, (lboundary(3) .and. .not. lperiodic(3))),jend)
      call vadvecw_6th(w0,wp,istart,iend,jstart,jend)
    case(iadv_62)
      call vadvecu_2nd(u0,up,merge(istart+1, istart, (lboundary(1) .and. .not. lperiodic(1))),iend,jstart,jend)
      call vadvecv_2nd(v0,vp,istart,iend,merge(jstart+1, jstart, (lboundary(3) .and. .not. lperiodic(3))),jend)
      call vadvecw_2nd(w0,wp,istart,iend,jstart,jend)
    case(iadv_hybrid)
      !if (.not. leq) stop "advec_5th does not support a non-uniform vertical grid."
      call vadvecu_5th(u0,up,istart,iend,jstart,jend)
      call vadvecv_5th(v0,vp,istart,iend,jstart,jend)
      call vadvecw_5th(w0,wp,istart,iend,jstart,jend)
    case(iadv_hybrid_f)
      !if (.not. leq) stop "advec_5th does not support a non-uniform vertical grid."
      call vadvecu_5th(u0,up,istart,iend,jstart,jend)
      call vadvecv_5th(v0,vp,istart,iend,jstart,jend)
      call vadvecw_5th(w0,wp,istart,iend,jstart,jend)
  end select
  !$acc wait

  istart = 2; iend = i1; jstart = 2; jend = j1
  select case(iadv_thl)
    case(iadv_cd2)
      call vadvecc_2nd(thl0,thlp,istart,iend,jstart,jend)
    case(iadv_5th)
!       !if (.not. leq) stop "advec_5th does not support a non-uniform vertical grid."
      call vadvecc_5th(thl0,thlp,istart,iend,jstart,jend)
    case(iadv_52)
      call vadvecc_2nd(thl0,thlp,istart,iend,jstart,jend)
    case(iadv_cd6)
!       !if (.not. leq) stop "advec_6th does not support a non-uniform vertical grid."
      call vadvecc_6th(thl0,thlp,istart,iend,jstart,jend)
    case(iadv_62)
      call vadvecc_2nd(thl0,thlp,istart,iend,jstart,jend)
    case(iadv_kappa)
      call vadvecc_kappa(thl0,thlp)
    case(iadv_upw)
      call vadvecc_upw(thl0,thlp)
    case(iadv_hybrid)
       !if (.not. leq) stop "advec_hybrid does not support a non-uniform vertical grid."
      call vadvecc_hybrid(thl0,thlp)
    case(iadv_hybrid_f)
      !if (.not. leq) stop "advec_hybrid_f does not support a non-uniform vertical grid."
      call vadvecc_hybrid_f(thl0,thlp,1.0)
  end select
  if (lmoist) then
    istart = 2; iend = i1; jstart = 2; jend = j1
    select case(iadv_qt)
      case(iadv_cd2)
        call vadvecc_2nd(qt0,qtp,istart,iend,jstart,jend)
      case(iadv_5th)
!         !if (.not. leq) stop "advec_5th does not support a non-uniform vertical grid."
        call vadvecc_5th(qt0,qtp,istart,iend,jstart,jend)
      case(iadv_52)
        call vadvecc_2nd(qt0,qtp,istart,iend,jstart,jend)
      case(iadv_cd6)
!         !if (.not. leq) stop "advec_6th does not support a non-uniform vertical grid."
        call vadvecc_6th(qt0,qtp,istart,iend,jstart,jend)
      case(iadv_62)
        call vadvecc_2nd(qt0,qtp,istart,iend,jstart,jend)
      case(iadv_kappa)
        call vadvecc_kappa(qt0,qtp)
      case(iadv_upw)
        call vadvecc_upw(qt0,qtp)
      case(iadv_hybrid)
        !if (.not. leq) stop "advec_hybrid does not support a non-uniform vertical grid."
        call vadvecc_hybrid(qt0,qtp)
      case(iadv_hybrid_f)
        !if (.not. leq) stop "advec_hybrid_f does not support a non-uniform vertical grid."
        call vadvecc_hybrid_f(qt0,qtp,1e-3)
    end select
  end if

  do n=1,nsv
    istart = 2; iend = i1; jstart = 2; jend = j1
    select case(iadv_sv)
    case(iadv_cd2)
      call vadvecc_2nd(sv0(:,:,:,n),svp(:,:,:,n),istart,iend,jstart,jend)
    case(iadv_5th)
!       !if (.not. leq) stop "advec_5th does not support a non-uniform vertical grid."
      call vadvecc_5th(sv0(:,:,:,n),svp(:,:,:,n),istart,iend,jstart,jend)
    case(iadv_52)
      call vadvecc_2nd(sv0(:,:,:,n),svp(:,:,:,n),istart,iend,jstart,jend)
    case(iadv_cd6)
!       !if (.not. leq) stop "advec_6th does not support a non-uniform vertical grid."
      call vadvecc_6th(sv0(:,:,:,n),svp(:,:,:,n),istart,iend,jstart,jend)
    case(iadv_62)
      call vadvecc_2nd(sv0(:,:,:,n),svp(:,:,:,n),istart,iend,jstart,jend)
    case(iadv_kappa)
      call vadvecc_kappa(sv0(:,:,:,n),svp(:,:,:,n))
    case(iadv_upw)
      call vadvecc_upw(sv0(:,:,:,n),svp(:,:,:,n))
    case(iadv_hybrid)
      !if (.not. leq) stop "advec_hybrid does not support a non-uniform vertical grid."
      call vadvecc_hybrid(sv0(:,:,:,n),svp(:,:,:,n))
    case(iadv_hybrid_f)
      !if (.not. leq) stop "advec_hybrid_f does not support a non-uniform vertical grid."
      call vadvecc_hybrid_f(sv0(:,:,:,n),svp(:,:,:,n))
    end select
  end do

  ! Horizontal and vertical tke advection together
  if (.not. lsmagorinsky) then
    istart = 2; iend = i1; jstart = 2; jend = j1
    select case(iadv_tke)
      case(iadv_cd2)
        call hadvecc_2nd(e120,e12p,istart,iend,jstart,jend)
        call vadvecc_2nd(e120,e12p,istart,iend,jstart,jend)
      case(iadv_5th)
        !if (.not. leq) stop "advec_5th does not support a non-uniform vertical grid."
        ibuffer = 2; jbuffer = 2
        if(lopenbc) call advecc_2nd_boundary_buffer(e120,e12p,istart,iend,jstart,jend,ibuffer,jbuffer)
        call hadvecc_5th(e120,e12p,istart,iend,jstart,jend)
        istart = 2; iend = i1; jstart = 2; jend = j1
        call vadvecc_5th(e120,e12p,istart,iend,jstart,jend)
      case(iadv_52)
        ibuffer = 2; jbuffer = 2
        if(lopenbc) call advecc_2nd_boundary_buffer(e120,e12p,istart,iend,jstart,jend,ibuffer,jbuffer)
        call hadvecc_5th(e120,e12p,istart,iend,jstart,jend)
        istart = 2; iend = i1; jstart = 2; jend = j1
        call vadvecc_2nd(e120,e12p,istart,iend,jstart,jend)
      case(iadv_cd6)
        !if (.not. leq) stop "advec_6th does not support a non-uniform vertical grid."
        ibuffer = 2; jbuffer = 2
        if(lopenbc) call advecc_2nd_boundary_buffer(e120,e12p,istart,iend,jstart,jend,ibuffer,jbuffer)
        call hadvecc_6th(e120,e12p,istart,iend,jstart,jend)
        istart = 2; iend = i1; jstart = 2; jend = j1
        call vadvecc_6th(e120,e12p,istart,iend,jstart,jend)
      case(iadv_62)
        ibuffer = 2; jbuffer = 2
        if(lopenbc) call advecc_2nd_boundary_buffer(e120,e12p,istart,iend,jstart,jend,ibuffer,jbuffer)
        call hadvecc_6th(e120,e12p,istart,iend,jstart,jend)
        istart = 2; iend = i1; jstart = 2; jend = j1
        call vadvecc_2nd(e120,e12p,istart,iend,jstart,jend)
      case(iadv_kappa)
        call hadvecc_kappa(e120,e12p)
        call vadvecc_kappa(e120,e12p)
      case(iadv_hybrid)
        !if (.not. leq) stop "advec_hybrid does not support a non-uniform vertical grid."
        call hadvecc_hybrid(e120,e12p)
        call vadvecc_hybrid(e120,e12p)
      case(iadv_hybrid_f)
        !if (.not. leq) stop "advec_hybrid_f does not support a non-uniform vertical grid."
        call hadvecc_hybrid_f(e120,e12p)
        call vadvecc_hybrid_f(e120,e12p)
      case(iadv_null)
        ! null advection scheme
        stop "Null advection scheme selected for iadv_tke - probably a bad idea."
      case default
        stop "Unknown advection scheme "
    end select
  end if
  call samptend(tend_vadv)

  !$acc wait
  call timer_toc('modadvection/advection')
end subroutine advection
end module modadvection
