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
contains
!> Advection redirection function
subroutine advection

  use modglobal,      only : lmoist, nsv, iadv_mom,iadv_tke,iadv_thl,iadv_qt,iadv_sv, &
                             iadv_cd2,iadv_5th,iadv_52,iadv_cd6,iadv_62,iadv_kappa,iadv_upw,iadv_hybrid,iadv_hybrid_f,iadv_null,leq
  use modfields,      only : u0,up,v0,vp,w0,wp,e120,e12p,thl0,thlp,qt0,qtp,sv0,svp
  use modsubgrid,     only : lsmagorinsky
  use advec_2nd,      only : advecu_2nd, advecv_2nd, advecw_2nd, advecc_2nd
  use advec_52,       only : advecu_52, advecv_52, advecw_52, advecc_52
  use advec_5th,      only : advecu_5th, advecv_5th, advecw_5th, advecc_5th
  use advec_62,       only : advecu_62, advecv_62, advecw_62, advecc_62
  use advec_6th,      only : advecu_6th, advecv_6th, advecw_6th, advecc_6th
  use advec_hybrid,   only : advecc_hybrid
  use advec_hybrid_f, only : advecc_hybrid_f
  use advec_kappa,    only : advecc_kappa
  use advec_upw,      only : advecc_upw
  implicit none
  integer :: n

  ! leq = .false. ! for testing that the non-uniform advection routines agree with the uniform ones
                  ! when the grid is uniform
  
  select case(iadv_mom)
    case(iadv_cd2)
      call advecu_2nd(u0,up)
      call advecv_2nd(v0,vp)
      call advecw_2nd(w0,wp)
    case(iadv_5th)
      !if (.not. leq) stop "advec_5th does not support a non-uniform vertical grid."
      call advecu_5th(u0,up)
      call advecv_5th(v0,vp)
      call advecw_5th(w0,wp)
    case(iadv_52)
      call advecu_52(u0,up)
      call advecv_52(v0,vp)
      call advecw_52(w0,wp)
    case(iadv_cd6)
      !if (.not. leq) stop "advec_6th does not support a non-uniform vertical grid."
      call advecu_6th(u0,up)
      call advecv_6th(v0,vp)
      call advecw_6th(w0,wp)
    case(iadv_62)
      call advecu_62(u0,up)
      call advecv_62(v0,vp)
      call advecw_62(w0,wp)
    case(iadv_hybrid)
      !if (.not. leq) stop "advec_5th does not support a non-uniform vertical grid."
      call advecu_5th(u0,up)
      call advecv_5th(v0,vp)
      call advecw_5th(w0,wp)
    case(iadv_hybrid_f)
      !if (.not. leq) stop "advec_5th does not support a non-uniform vertical grid."
      call advecu_5th(u0,up)
      call advecv_5th(v0,vp)
      call advecw_5th(w0,wp)
    case(iadv_null)
      ! null advection scheme 
      stop "Null advection scheme selected for iadv_mom - probably a bad idea."
    case default
      stop "Unknown advection scheme "
  end select

  if (.not. lsmagorinsky) then
    select case(iadv_tke)
      case(iadv_cd2)
        call advecc_2nd(e120,e12p)
      case(iadv_5th)
        !if (.not. leq) stop "advec_5th does not support a non-uniform vertical grid."
        call advecc_5th(e120,e12p)
      case(iadv_52)
        call advecc_52(e120,e12p)
      case(iadv_cd6)
        !if (.not. leq) stop "advec_6th does not support a non-uniform vertical grid."
        call advecc_6th(e120,e12p)
      case(iadv_62)
        call advecc_62(e120,e12p)
      case(iadv_kappa)
        call advecc_kappa(e120,e12p)
      case(iadv_hybrid)
        !if (.not. leq) stop "advec_hybrid does not support a non-uniform vertical grid."
        call advecc_hybrid(e120,e12p)
      case(iadv_hybrid_f)
        !if (.not. leq) stop "advec_hybrid_f does not support a non-uniform vertical grid."
        call advecc_hybrid_f(e120,e12p)         
      case(iadv_null)
        ! null advection scheme 
        stop "Null advection scheme selected for iadv_tke - probably a bad idea."
      case default
        stop "Unknown advection scheme "
    end select
  end if

  select case(iadv_thl)
    case(iadv_cd2)
      call advecc_2nd(thl0,thlp)
    case(iadv_5th)
      !if (.not. leq) stop "advec_5th does not support a non-uniform vertical grid."
      call advecc_5th(thl0,thlp)
    case(iadv_52)
      call advecc_52(thl0,thlp)
    case(iadv_cd6)
      !if (.not. leq) stop "advec_6th does not support a non-uniform vertical grid."
      call advecc_6th(thl0,thlp)
    case(iadv_62)
      call advecc_62(thl0,thlp)
    case(iadv_kappa)
      call advecc_kappa(thl0,thlp)
    case(iadv_upw)
      if (.not. leq) stop "advec_upw does not support a non-uniform vertical grid."
      call advecc_upw(thl0,thlp)
    case(iadv_hybrid)
       !if (.not. leq) stop "advec_hybrid does not support a non-uniform vertical grid."
      call advecc_hybrid(thl0,thlp)
    case(iadv_hybrid_f)
      !if (.not. leq) stop "advec_hybrid_f does not support a non-uniform vertical grid."
      call advecc_hybrid_f(thl0,thlp,1.0)
    case(iadv_null)
      ! null advection scheme 
      stop "Null advection scheme selected for iadv_thl - probably a bad idea." 
    case default
      stop "Unknown advection scheme "
  end select
  if (lmoist) then
    select case(iadv_qt)
      case(iadv_cd2)
        call advecc_2nd(qt0,qtp)
      case(iadv_5th)
        !if (.not. leq) stop "advec_5th does not support a non-uniform vertical grid."
        call advecc_5th(qt0,qtp)
      case(iadv_52)
        call advecc_52(qt0,qtp)
      case(iadv_cd6)
        !if (.not. leq) stop "advec_6th does not support a non-uniform vertical grid."
        call advecc_6th(qt0,qtp)
      case(iadv_62)
        call advecc_62(qt0,qtp)
      case(iadv_kappa)
        call advecc_kappa(qt0,qtp)
      case(iadv_upw)
        if (.not. leq) stop "advec_upw does not support a non-uniform vertical grid."
        call advecc_upw(qt0,qtp)
      case(iadv_hybrid)
        !if (.not. leq) stop "advec_hybrid does not support a non-uniform vertical grid."
        call advecc_hybrid(qt0,qtp)
      case(iadv_hybrid_f)
        !if (.not. leq) stop "advec_hybrid_f does not support a non-uniform vertical grid."
        call advecc_hybrid_f(qt0,qtp,1e-3)
      case(iadv_null)
        ! null advection scheme 
        stop "Null advection scheme selected for iadv_qt - probably a bad idea."
      case default
        stop "Unknown advection scheme "
    end select
  end if
  do n=1,nsv
    select case(iadv_sv(n))
    case(iadv_cd2)
      call advecc_2nd(sv0(:,:,:,n),svp(:,:,:,n))
    case(iadv_5th)
      !if (.not. leq) stop "advec_5th does not support a non-uniform vertical grid."
      call advecc_5th(sv0(:,:,:,n),svp(:,:,:,n))
    case(iadv_52)
      call advecc_52(sv0(:,:,:,n),svp(:,:,:,n))
    case(iadv_cd6)
      !if (.not. leq) stop "advec_6th does not support a non-uniform vertical grid."
      call advecc_6th(sv0(:,:,:,n),svp(:,:,:,n))
    case(iadv_62)
      call advecc_62(sv0(:,:,:,n),svp(:,:,:,n))
    case(iadv_kappa)
      call advecc_kappa(sv0(:,:,:,n),svp(:,:,:,n))
    case(iadv_upw)
      if (.not. leq) stop "advec_upw does not support a non-uniform vertical grid."
      call advecc_upw(sv0(:,:,:,n),svp(:,:,:,n))
    case(iadv_hybrid)
      !if (.not. leq) stop "advec_hybrid does not support a non-uniform vertical grid."
      call advecc_hybrid(sv0(:,:,:,n),svp(:,:,:,n))
    case(iadv_hybrid_f)
      !if (.not. leq) stop "advec_hybrid_f does not support a non-uniform vertical grid."
      call advecc_hybrid_f(sv0(:,:,:,n),svp(:,:,:,n))      
    case(iadv_null)
       ! null advection scheme - do nothing
    case default
      stop "Unknown advection scheme "
    end select
  end do

end subroutine advection
end module modadvection
