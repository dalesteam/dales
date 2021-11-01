!> \file modfields.f90
!!  Declares, allocates and initializes the 3D fields

!>
!!  Declares, allocates and initializes the 3D fields
!>

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

module modfields
use modprecision

implicit none
save

  ! Prognostic variables

  real(field_r), allocatable :: um(:,:,:)        !<   x-component of velocity at time step t-1
  real(field_r), allocatable :: vm(:,:,:)        !<   y-component of velocity at time step t-1
  real(field_r), allocatable :: wm(:,:,:)        !<   z-component of velocity at time step t-1
  real(field_r), allocatable :: thlm(:,:,:)      !<   liq. water pot. temperature at time step t-1
  real(field_r), allocatable :: e12m(:,:,:)      !<   square root of turb. kin. energy at time step t-1
  real(field_r), allocatable :: qtm(:,:,:)       !<   total specific humidity at time step t
  real(field_r), allocatable :: u0(:,:,:)        !<   x-component of velocity at time step t
  real(field_r), allocatable :: v0(:,:,:)        !<   y-component of velocity at time step t
  real(field_r), allocatable :: w0(:,:,:)        !<   z-component of velocity at time step t
  real(field_r), allocatable :: thl0(:,:,:)      !<   liq. water pot. temperature at time step t
  real(field_r), allocatable :: thl0h(:,:,:)     !<  3d-field of theta_l at half levels for kappa scheme
  real(field_r), allocatable :: qt0h(:,:,:)      !<  3d-field of q_tot   at half levels for kappa scheme
  real(field_r), allocatable :: e120(:,:,:)      !<   square root of turb. kin. energy at time step t
  real(field_r), allocatable :: qt0(:,:,:)       !<   total specific humidity at time step t

  real(field_r), allocatable :: up(:,:,:)        !<   tendency of um
  real(field_r), allocatable :: vp(:,:,:)        !<   tendency of vm
  real(field_r), allocatable :: wp(:,:,:)        !<   tendency of wm
  real(field_r), allocatable :: thlp(:,:,:)      !<   tendency of thlm
  real(field_r), allocatable :: e12p(:,:,:)      !<   tendency of e12m
  real(field_r), allocatable :: qtp(:,:,:)       !<   tendency of qtm

  real(field_r), allocatable :: svm(:,:,:,:)   !<  scalar sv(n) at time step t-1
  real(field_r), allocatable :: sv0(:,:,:,:)   !<  scalar sv(n) at time step t
  real(field_r), allocatable :: svp(:,:,:,:)   !<  tendency of sv(n)

  ! Base state variables
  real(field_r), allocatable :: rhobf(:)       !<   Base state density, full level
  real(field_r), allocatable :: rhobh(:)       !<   Base state density, half level

  real(field_r), allocatable :: drhobdzf(:)       !<   Base state density, derivative at full level
  real(field_r), allocatable :: drhobdzh(:)       !<   Base state density, derivative at half level

  ! Cloud edge variables
  real(field_r), allocatable :: ql0(:,:,:)  !<   liquid water content
  real(field_r), allocatable :: tmp0(:,:,:) !<   temperature at full level
  real(field_r), allocatable :: thv0h(:,:,:)!<   theta_v at half level

  real(field_r), allocatable :: whls(:)                       !<   large scale vert velocity at half levels

  real(field_r), allocatable :: presf(:)                      !<   hydrostatic pressure at full level
  real(field_r), allocatable :: presh(:)                      !<   hydrostatic pressure at half level
  real(field_r), allocatable :: initial_presf(:)              !<   initial hydrostatic pressure at full level
  real(field_r), allocatable :: initial_presh(:)              !<   initial hydrostatic pressure at half level
  real(field_r), allocatable :: exnf(:)                       !<   hydrostatic exner function at full level
  real(field_r), allocatable :: exnh(:)                       !<   hydrostatic exner function at half level
  real(field_r), allocatable :: thvf(:)                       !<   hydrostatic thetav at full level
  real(field_r), allocatable :: thvh(:)                       !<   hydrostatic thetav at half level
  real(field_r), allocatable :: rhof(:)                       !<   slab averaged density at full level
  real(field_r), allocatable :: qt0av(:)                      !<   slab averaged q_tot
  real(field_r), allocatable :: ql0av(:)                      !<   slab averaged q_liq

  real(field_r), allocatable :: thl0av(:)                     !<   slab averaged th_liq
  real(field_r), allocatable :: u0av(:)                       !<   slab averaged u
  real(field_r), allocatable :: v0av(:)                       !<   slab averaged v
  real(field_r), allocatable :: ug(:)                       !<   geostrophic u-wind
  real(field_r), allocatable :: vg(:)                       !<   geostrophic v-wind

  real(field_r), allocatable :: dpdxl(:)                      !<   large scale pressure x-gradient
  real(field_r), allocatable :: dpdyl(:)                      !<   large scale pressure y-gradient

  real(field_r), allocatable :: dthldxls(:)                   !<   large scale x-gradient of th_liq
  real(field_r), allocatable :: dthldyls(:)                   !<   large scale y-gradient of th_liq
  real(field_r), allocatable :: dthldtls(:)                   !<   large scale tendency of thl

  real(field_r), allocatable :: dqtdxls(:)                    !<   large scale x-gradient of q_tot
  real(field_r), allocatable :: dqtdyls(:)                    !<   large scale y-gradient of q_tot
  real(field_r), allocatable :: dqtdtls(:)                    !<   large scale tendency of q_tot

  real(field_r), allocatable :: dudxls(:)                     !<   large scale x-gradient of u
  real(field_r), allocatable :: dudyls(:)                     !<   large scale y-gradient of u
  real(field_r), allocatable :: dudtls(:)                     !<   large scale tendency of u

  real(field_r), allocatable :: dvdxls(:)                     !<   large scale x-gradient of v
  real(field_r), allocatable :: dvdyls(:)                     !<   large scale y-gradient of v
  real(field_r), allocatable :: dvdtls(:)                     !<   large scale tendency of v

  real(field_r), allocatable :: wfls  (:)                     !<   large scale vertical velocity
  real(field_r), allocatable :: ql0h(:,:,:)
  real(field_r), allocatable :: dthvdz(:,:,:)!<   theta_v at half level

  real(field_r), allocatable :: thlprof(:)                    !<   initial thl-profile
  real(field_r), allocatable :: qtprof(:)                     !<   initial qt-profile
  real(field_r), allocatable :: uprof(:)                      !<   initial u-profile
  real(field_r), allocatable :: vprof(:)                      !<   initial v-profile
  real(field_r), allocatable :: e12prof(:)                    !<   initial subgrid sqrt(TKE) profile
  real(field_r), allocatable :: sv0av(:,:)                  !<   slab average of sv(n)
  real(field_r), allocatable :: svprof(:,:)                 !<   initial sv(n)-profile

  real(field_r), allocatable :: thlpcar(:)                    !< prescribed radiatively forced thl tendency
  real(field_r), allocatable :: SW_up_TOA(:,:), SW_dn_TOA(:,:), LW_up_TOA(:,:), LW_dn_TOA(:,:)
  real(field_r), allocatable :: qvsl(:,:,:)
  real(field_r), allocatable :: qvsi(:,:,:)
  real(field_r), allocatable :: esl(:,:,:)

  real(field_r), allocatable :: qsat(:,:,:)
  real(field_r), allocatable :: surf_rain(:,:)               !< integrated surface rain 

contains
!> Allocate and initialize the prognostic variables
subroutine initfields

    use modglobal, only : i1,ih,j1,jh,k1,nsv
    ! Allocation of prognostic variables
    implicit none

    allocate(um   (2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(vm   (2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(wm   (2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(thlm (2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(e12m (2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(qtm  (2-ih:i1+ih,2-jh:j1+jh,k1))

    allocate(u0   (2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(v0   (2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(w0   (2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(thl0 (2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(e120 (2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(qt0  (2-ih:i1+ih,2-jh:j1+jh,k1))

    allocate(up   (2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(vp   (2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(wp   (2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(thlp (2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(e12p (2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(qtp  (2-ih:i1+ih,2-jh:j1+jh,k1))

    allocate(svm  (2-ih:i1+ih,2-jh:j1+jh,k1,nsv))
    allocate(sv0  (2-ih:i1+ih,2-jh:j1+jh,k1,nsv))
    allocate(svp  (2-ih:i1+ih,2-jh:j1+jh,k1,nsv))

    allocate(thl0h(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(qt0h (2-ih:i1+ih,2-jh:j1+jh,k1))

    ! Allocation of base state variables
    allocate(rhobf   (k1))
    allocate(rhobh   (k1))
    allocate(drhobdzf(k1))
    allocate(drhobdzh(k1))

    ! Allocation of diagnostic variables
    allocate(ql0   (2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(ql0h  (2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(tmp0  (2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(thv0h (2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(dthvdz(2-ih:i1+ih,2-jh:j1+jh,k1))

    allocate(whls         (k1))
    allocate(presf        (k1))
    allocate(presh        (k1))
    allocate(initial_presf(k1))
    allocate(initial_presh(k1))
    allocate(exnf         (k1))
    allocate(exnh         (k1))
    allocate(thvf         (k1))
    allocate(thvh         (k1))
    allocate(rhof         (k1))
    allocate(qt0av        (k1))
    allocate(ql0av        (k1))
    allocate(thl0av       (k1))
    allocate(u0av         (k1))
    allocate(v0av         (k1))
    allocate(ug           (k1))
    allocate(vg           (k1))
    allocate(dpdxl        (k1))
    allocate(dpdyl        (k1))

    allocate(dthldxls(k1))
    allocate(dthldyls(k1))
    allocate(dthldtls(k1))

    allocate(dqtdxls(k1))
    allocate(dqtdyls(k1))
    allocate(dqtdtls(k1))

    allocate(dudxls(k1))
    allocate(dudyls(k1))
    allocate(dudtls(k1))

    allocate(dvdxls(k1))
    allocate(dvdyls(k1))
    allocate(dvdtls(k1))

    allocate(wfls   (k1))
    allocate(thlprof(k1))
    allocate(qtprof (k1))
    allocate(uprof  (k1))
    allocate(vprof  (k1))
    allocate(e12prof(k1))
    allocate(sv0av  (k1,nsv))
    allocate(svprof (k1,nsv))
    allocate(thlpcar(k1))

    allocate(SW_up_TOA(2-ih:i1+ih,2-jh:j1+jh))
    allocate(SW_dn_TOA(2-ih:i1+ih,2-jh:j1+jh))
    allocate(LW_up_TOA(2-ih:i1+ih,2-jh:j1+jh))
    allocate(LW_dn_TOA(2-ih:i1+ih,2-jh:j1+jh))

    allocate (qvsl(2-ih:i1+ih,2-jh:j1+jh,k1)    & ! qv-liquid
             ,qvsi(2-ih:i1+ih,2-jh:j1+jh,k1)    & ! qv-ice
             ,esl (2-ih:i1+ih,2-jh:j1+jh,k1)    & ! es-liquid
             ,qsat(2-ih:i1+ih,2-jh:j1+jh,k1))

    allocate(surf_rain(2-ih:i1+ih,2-jh:j1+jh))

    um=0.;u0=0.;up=0.
    vm=0.;v0=0.;vp=0.
    wm=0.;w0=0.;wp=0.
    thlm=0.;thl0=0.;thlp=0.
    qtm=0.;qt0=0.;qtp=0.
    e12m=0.;e120=0.;e12p=0.
    svm=0.;sv0=0.;svp=0.

    rhobf=0.;rhobh=0.;drhobdzf=0.;drhobdzh=0.
    ql0=0.;tmp0=0.;ql0h=0.;thv0h=0.;thl0h=0.;qt0h=0.
    presf=0.;presh=0.;exnf=0.;exnh=0.;thvh=0.;thvf=0.;rhof=0.    ! OG
    qt0av=0.;ql0av=0.;thl0av=0.;u0av=0.;v0av=0.;sv0av=0.
    thlprof=0.;qtprof=0.;uprof=0.;vprof=0.;e12prof=0.;svprof=0.
    ug=0.;vg=0.;dpdxl=0.;dpdyl=0.;wfls=0.;whls=0.
    thlpcar = 0.
    dthldxls=0.;dthldyls=0.;dthldtls=0.
    dqtdxls=0.;dqtdyls=0.;dqtdtls=0.
    dudxls=0.;dudyls=0.;dudtls=0.
    dvdxls=0.;dvdyls=0.;dvdtls=0.
    dthvdz=0.
    SW_up_TOA=0.;SW_dn_TOA=0.;LW_up_TOA=0.;LW_dn_TOA=0.
    qvsl=0.;qvsi=0.;esl=0.
    qsat=0.

    surf_rain = 0
  end subroutine initfields

!> Deallocate the fields
  subroutine exitfields
  implicit none
    deallocate(um,vm,wm,thlm,e12m,qtm,u0,v0,w0,thl0,thl0h,qt0h,e120,qt0)
    deallocate(up,vp,wp,thlp,e12p,qtp)
    deallocate(svm,sv0,svp)
    deallocate(rhobf,rhobh)
    deallocate(drhobdzf,drhobdzh)
    deallocate(ql0,tmp0,ql0h,thv0h,dthvdz,whls,presf,presh,initial_presf,initial_presh,exnf,exnh,thvh,thvf,rhof,qt0av,ql0av,thl0av,u0av,v0av)
    deallocate(ug,vg,dpdxl,dpdyl,wfls)
    deallocate(dthldxls,dthldyls,dthldtls,dqtdxls,dqtdyls,dqtdtls)
    deallocate(dudxls,dudyls,dudtls,dvdxls,dvdyls,dvdtls)
    deallocate(thlprof,qtprof,uprof,vprof,e12prof,sv0av,svprof)
    deallocate(thlpcar)
    deallocate(SW_up_TOA,SW_dn_TOA,LW_up_TOA,LW_dn_TOA)
    deallocate(qvsl,qvsi,esl)
    deallocate(qsat)
    deallocate(surf_rain)
    end subroutine exitfields

end module modfields
