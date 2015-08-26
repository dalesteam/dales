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

implicit none
save

  ! Prognostic variables

  real, allocatable :: u0(:,:,:)        !<   x-component of velocity at time step t
  real, allocatable :: v0(:,:,:)        !<   y-component of velocity at time step t
  real, allocatable :: w0(:,:,:)        !<   z-component of velocity at time step t
  real, allocatable :: thl0(:,:,:)      !<   liq. water pot. temperature at time step t
  real, allocatable :: thl0h(:,:,:)     !<  3d-field of theta_l at half levels for kappa scheme
  real, allocatable :: qt0h(:,:,:)      !<  3d-field of q_tot   at half levels for kappa scheme
  real, allocatable :: e120(:,:,:)      !<   turb. kin. energy at time step t
  real, allocatable :: qt0(:,:,:)       !<   total specific humidity at time step t

  ! Subgrid variables, moved from subgriddata !JvdD
  real, allocatable :: ekm(:,:,:)       !< k-coefficient for momentum
  real, allocatable :: ekh(:,:,:)       !< k-coefficient for heat and q_tot

  ! Prognostic fields at previous timestep
  !  Only used for the Wicker and Skamarock RK3 scheme
  !  Unnecessary when using the Williamson low storage RK3 scheme
  real, allocatable :: um(:,:,:)        !<   x-component of velocity at time step t-1
  real, allocatable :: vm(:,:,:)        !<   y-component of velocity at time step t-1
  real, allocatable :: wm(:,:,:)        !<   z-component of velocity at time step t-1
  real, allocatable :: thlm(:,:,:)      !<   liq. water pot. temperature at time step t-1
  real, allocatable :: qtm(:,:,:)       !<   total specific humidity at time step t
  real, allocatable :: e12m(:,:,:)      !<   turb. kin. energy at time step t-1
  real, allocatable :: svm(:,:,:,:)     !<  scalar sv(n) at time step t-1

  ! Tendencies
  real, allocatable :: up(:,:,:)        !<   tendency of um
  real, allocatable :: vp(:,:,:)        !<   tendency of vm
  real, allocatable :: wp(:,:,:)        !<   tendency of wm
  real, allocatable :: wp_store(:,:,:)  !<   tendency of wm, dummy variable for w-budget sampling
  real, allocatable :: thlp(:,:,:)      !<   tendency of thlm
  real, allocatable :: e12p(:,:,:)      !<   tendency of e12m
  real, allocatable :: qtp(:,:,:)       !<   tendency of qtm

  real, allocatable :: sv0(:,:,:,:)   !<  scalar sv(n) at time step t
  real, allocatable :: svp(:,:,:,:)   !<  tendency of sv(n)

  ! Base state variables
  real, allocatable :: rhobf(:)       !<   Base state density, full level
  real, allocatable :: rhobh(:)       !<   Base state density, half level

  real, allocatable :: drhobdzf(:)       !<   Base state density, derivative at full level
  real, allocatable :: drhobdzh(:)       !<   Base state density, derivative at half level

  real, allocatable :: rhodzi(:)      !<   Inverse of rhobf*dzf

  ! Surface gradient of thl and qt
  ! Moved these from surfdata because they are also used in thermodynamics and are therefore not optional
  real, allocatable :: dthldz(:,:)
  real, allocatable :: dqtdz(:,:)

  ! Following fields are included in restarts. It is more clear to add them to
  ! modfields to avoid unallocated arrays when doing a restart.
  real, allocatable, dimension(:,:) :: tskin,   &
                                       qskin,   &
                                       obl,     &
                                       ustar,   &
                                       thlflux, &
                                       qtflux
  real, allocatable                 :: svflux(:,:,:)
  ! Patch variables that are also in the restart files
  real, allocatable, dimension(:,:) :: ps_patch,   &
                                       thls_patch, &
                                       qts_patch,  &
                                       thvs_patch, &
                                       oblpatch

  ! Cloud edge variables
  real, allocatable :: cloudarea(:,:,:)
  real, allocatable :: cloudnr(:,:,:)
  real, allocatable :: cloudnrold(:,:,:)
  real, allocatable :: distcld(:,:,:)
  real, allocatable :: distcr(:,:,:)
  real, allocatable :: distqr(:,:)
  real, allocatable :: distdiv(:,:)
  real, allocatable :: distcon(:,:)
  real, allocatable :: distbuoy(:,:)
  real, allocatable :: distw(:,:)

  real, allocatable :: ql0(:,:,:)  !<   liquid water content
  real, allocatable :: tmp0(:,:,:) !<   temperature at full level
  real, allocatable :: thv0h(:,:,:)!<   theta_v at half level

  real, allocatable :: whls(:)                       !<   large scale vert velocity at half levels

  real, allocatable :: presf(:)                      !<   hydrostatic pressure at full level
  real, allocatable :: presh(:)                      !<   hydrostatic pressure at half level
  real, allocatable :: exnf(:)                       !<   hydrostatic exner function at full level
  real, allocatable :: exnh(:)                       !<   hydrostatic exner function at half level
  real, allocatable :: thvf(:)                       !<   hydrostatic thetav at full level
  real, allocatable :: thvh(:)                       !<   hydrostatic thetav at half level
  real, allocatable :: rhof(:)                       !<   slab averaged density at full level
  real, allocatable :: qt0av(:)                      !<   slab averaged q_tot
  real, allocatable :: ql0av(:)                      !<   slab averaged q_liq

  real, allocatable :: thl0av(:)                     !<   slab averaged th_liq
  real, allocatable :: u0av(:)                       !<   slab averaged u
  real, allocatable :: v0av(:)                       !<   slab averaged v
  real, allocatable :: ug(:)                       !<   geostrophic u-wind
  real, allocatable :: vg(:)                       !<   geostrophic v-wind

  real, allocatable :: dpdxl(:)                      !<   large scale pressure x-gradient
  real, allocatable :: dpdyl(:)                      !<   large scale pressure y-gradient

  real, allocatable :: dthldxls(:)                   !<   large scale x-gradient of th_liq
  real, allocatable :: dthldyls(:)                   !<   large scale y-gradient of th_liq
  real, allocatable :: dqtdxls(:)                    !<   large scale x-gradient of q_tot
  real, allocatable :: dqtdyls(:)                    !<   large scale y-gradient of q_tot
  real, allocatable :: dqtdtls(:)                    !<   large scale y-gradient of q_tot
  real, allocatable :: dudxls(:)                     !<   large scale x-gradient of u

  real, allocatable :: dudyls(:)                     !<   large scale y-gradient of u
  real, allocatable :: dvdxls(:)                     !<   large scale x-gradient of v
  real, allocatable :: dvdyls(:)                     !<   large scale y-gradient of v
  real, allocatable :: wfls  (:)                     !<   large scale y-gradient of v
  real, allocatable :: ql0h(:,:,:)
  real, allocatable :: dthvdz(:,:,:)!<   theta_v at half level

  real, allocatable :: thlprof(:)                    !<   initial thl-profile
  real, allocatable :: qtprof(:)                     !<   initial qt-profile
  real, allocatable :: uprof(:)                      !<   initial u-profile
  real, allocatable :: vprof(:)                      !<   initial v-profile
  real, allocatable :: e12prof(:)                    !<   initial subgrid TKE profile
  real, allocatable :: sv0av(:,:)                  !<   slab average of sv(n)
  real, allocatable :: svprof(:,:)                 !<   initial sv(n)-profile

  real, allocatable :: thlpcar(:)                    !< prescribed radiatively forced thl tendency
  real, allocatable :: qvsl(:,:,:)
  real, allocatable :: qvsi(:,:,:)
  real, allocatable :: esl(:,:,:)

  ! Store (smooth or not) in a 1d array
  real, allocatable :: smoothx(:),smoothy(:),smoothz(:) !< Smoothnesses are all based on qt!!

contains
!> Allocate and initialize the prognostic variables
subroutine initfields

    use modglobal, only : i1,i2,ih,j1,j2,jh,kmax,k1,nsv,iTimeInt,iTimeWicker,iTimeTVD
    use modsurfdata, only : lhetero,xpatches,ypatches
    ! Allocation of prognostic variable
    implicit none

    allocate(u0(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(v0(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(w0(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(thl0(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(thl0h(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(qt0h(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(e120(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(qt0(2-ih:i1+ih,2-jh:j1+jh,k1))

    allocate(ekm(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(ekh(2-ih:i1+ih,2-jh:j1+jh,k1))

    ! Only allocate these fields when required JvdD
    if (iTimeInt==iTimeWicker .or. iTimeInt==iTimeTVD) then
      allocate(um(2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(vm(2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(wm(2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(thlm(2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(e12m(2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(qtm(2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(svm(2-ih:i1+ih,2-jh:j1+jh,k1,nsv))
      um=0;vm=0;wm=0;thlm=0;e12m=0;qtm=0;svm=0
    end if

    allocate(up(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(vp(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(wp(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(wp_store(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(thlp(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(e12p(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(qtp(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(sv0(2-ih:i1+ih,2-jh:j1+jh,k1,nsv))
    allocate(svp(2-ih:i1+ih,2-jh:j1+jh,k1,nsv))

    allocate(dqtdz (i2,j2))
    allocate(dthldz(i2,j2))
    allocate(tskin(i2,j2),   &
             qskin(i2,j2),   &
             obl(i2,j2),     &
             ustar(i2,j2),   &
             thlflux(i2,j2), &
             qtflux(i2,j2)   )
    allocate(svflux(i2,j2,nsv))
    if (lhetero) then
      allocate(ps_patch  (xpatches,ypatches), &
               thls_patch(xpatches,ypatches), &
               qts_patch (xpatches,ypatches), &
               thvs_patch(xpatches,ypatches), &
               oblpatch  (xpatches,ypatches)  )
    end if

    allocate(cloudarea(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(cloudnr(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(cloudnrold(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(distcld(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(distcr(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(distqr(2-ih:i1+ih,2-jh:j1+jh))
    allocate(distdiv(2-ih:i1+ih,2-jh:j1+jh))
    allocate(distcon(2-ih:i1+ih,2-jh:j1+jh))
    allocate(distbuoy(2-ih:i1+ih,2-jh:j1+jh))
    allocate(distw(2-ih:i1+ih,2-jh:j1+jh))

    ! Allocation of base state variables
    allocate(rhobf(k1))
    allocate(rhobh(k1))
    allocate(drhobdzf(k1))
    allocate(drhobdzh(k1))
    allocate(rhodzi(k1))

    ! Allocation of diagnostic variables
    allocate(ql0(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(tmp0(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(thv0h(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(whls(k1))
    allocate(presf(k1))
    allocate(presh(k1))
    allocate(exnf(k1))
    allocate(exnh(k1))
    allocate(thvf(k1))
    allocate(thvh(k1))
    allocate(rhof(k1))
    allocate(qt0av(k1))
    allocate(ql0av(k1))
    allocate(thl0av(k1))
    allocate(u0av(k1))
    allocate(v0av(k1))
    allocate(ug(k1))
    allocate(vg(k1))
    allocate(dpdxl(k1))
    allocate(dpdyl(k1))
    allocate(dthldxls(k1))
    allocate(dthldyls(k1))
    allocate(dqtdxls(k1))
    allocate(dqtdyls(k1))
    allocate(dqtdtls(k1))
    allocate(dudxls(k1))
    allocate(dudyls(k1))
    allocate(dvdxls(k1))
    allocate(dvdyls(k1))
    allocate(wfls  (k1))
    allocate(ql0h(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(dthvdz(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(thlprof(k1))
    allocate(qtprof(k1))
    allocate(uprof(k1))
    allocate(vprof(k1))
    allocate(e12prof(k1))
    allocate(sv0av(k1,nsv))
    allocate(svprof(k1,nsv))
    allocate(thlpcar(k1))

    allocate(smoothx(k1),smoothy(k1),smoothz(k1))

    allocate (qvsl(2-ih:i1+ih,2-jh:j1+jh,k1)    & ! qv-liquid
             ,qvsi(2-ih:i1+ih,2-jh:j1+jh,k1)    & ! qv ice
             ,esl(2-ih:i1+ih,2-jh:j1+jh,k1))     ! es-liquid

    u0=0.;up=0.
    v0=0.;vp=0.
    w0=0.;wp=0.;wp_store=0.
    thl0=0.;thlp=0.
    qt0=0.;qtp=0.
    e120=0.;e12p=0.
    sv0=0.;svp=0.

    ekm=0.;ekh=0.

    rhobf=0.;rhobh=0.;drhobdzf=0.;drhobdzh=0.; rhodzi=0.
    ql0=0.;tmp0=0.;ql0h=0.;thv0h=0.;thl0h=0.;qt0h=0.
    presf=0.;presh=0.;exnf=0.;exnh=0.;thvh=0.;thvf=0.;rhof=0.    ! OG
    qt0av=0.;ql0av=0.;thl0av=0.;u0av=0.;v0av=0.;sv0av=0.
    thlprof=0.;qtprof=0.;uprof=0.;vprof=0.;e12prof=0.;svprof=0.
    ug=0.;vg=0.;dpdxl=0.;dpdyl=0.;wfls=0.;whls=0.;thlpcar = 0.
    dthldxls=0.;dthldyls=0.;dqtdxls=0.;dqtdyls=0.;dudxls=0.;dudyls=0.;dvdxls=0.;dvdyls=0.
    dthvdz=0.;dthldz=0.;dqtdz=0.
    tskin=0.;qskin=0.;obl=0.;ustar=0.;thlflux=0.;qtflux=0.;svflux=0.
    if (lhetero) then
      ps_patch=-1.;thls_patch=-1.;qts_patch=-1.;thvs_patch=-1.;oblpatch=-.1
    end if
    qvsl=0.;qvsi=0.;esl=0.

    smoothx=0.; smoothy=0.; smoothz=0.

    cloudarea=0.;cloudnr=0.;cloudnrold=0.;distcld=0.;distcr=0.;distqr=0.;distdiv=0.;distcon=0.;distbuoy=0.;distw=0.

  end subroutine initfields

!> Deallocate the fields
  subroutine exitfields
  use modglobal, only : iTimeInt,iTimeWicker,iTimeTVD
  use modsurfdata, only : lhetero,xpatches,ypatches
  implicit none
    if (iTimeInt==iTimeWicker .or. iTimeInt==iTimeTVD) deallocate(um,vm,wm,thlm,e12m,qtm,svm)
    if (lhetero) deallocate(ps_patch,thls_patch,qts_patch,thvs_patch,oblpatch)
    deallocate(u0,v0,w0,thl0,thl0h,qt0h,e120,qt0)
    deallocate(ekm,ekh)
    deallocate(up,vp,wp,wp_store,thlp,e12p,qtp)
    deallocate(sv0,svp)
    deallocate(dthldz,dqtdz)
    deallocate(tskin,qskin,obl,ustar,thlflux,qtflux)
    deallocate(svflux)
    deallocate(rhobf,rhobh)
    deallocate(drhobdzf,drhobdzh)
    deallocate(rhodzi)
    deallocate(ql0,tmp0,ql0h,thv0h,dthvdz,whls,presf,presh,exnf,exnh,thvh,thvf,rhof,qt0av,ql0av,thl0av,u0av,v0av)
    deallocate(ug,vg,dpdxl,dpdyl,dthldxls,dthldyls,dqtdxls,dqtdyls,dqtdtls,dudxls,dudyls,dvdxls,dvdyls,wfls)
    deallocate(thlprof,qtprof,uprof,vprof,e12prof,sv0av,svprof)
    deallocate(thlpcar)
    deallocate(cloudarea,cloudnr,cloudnrold,distcld,distcr,distqr,distdiv,distcon,distbuoy,distw)
    deallocate(qvsl,qvsi,esl)
    deallocate(smoothx,smoothy,smoothz)
    end subroutine exitfields

end module modfields
