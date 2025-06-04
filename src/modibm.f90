!> \file modibm.f90
!!  Grid-conforming Immersed Boundary Method (IBM)

!> References:
!! (1) Pourquie, M., W.-P. Breugem, and B. J. Boersma, 2009: Some issues related to the use of immersed boundary methods
!! to represent square obstacles. International Journal for Multiscale Computational Engineering, 7 (6), 509–522.
!! (2) Tomas, J., 2016: Obstacle-resolving large-eddy simulation of dispersion in urban environments: Effects of stability and roughness geometry,
!! Delft University of Technology, Delft, The Netherlands. https://doi.org/10.4233/uuid:5d93a697-be49-4f63-b871-b763bc327139
!>
!!  \author Michael Koene, Delft University of Technology, 2018-2019
!!  \author Stephan de Roode, Delft University of Technology, 2018-2024
!!  \author Steven van der Linden, Delft University of Technology, 2025-
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
! Copyright 2025 Delft University of Technology
!

module modibm
  use modglobal,    only : rd, rv, grav, ijtot
  use modprecision, only : field_r
  use modsurface,   only : psim, psih
  use modsurfdata,  only : thvs
  use modibmdata,   only : lapply_ibm,lpoislast, lwallheat, &
                            thlwall, thlroof, qtroof, thlibm, qtibm, &
                            z0m_wall, z0h_wall
  use modtimer
  implicit none
  save
  private

  public :: ixw_p, ixw_m, iyw_p, iyw_m, izw_p, iobst

  !< Later extendable to put different wall types in fortran types
  integer :: Nxwalls_plus, Nywalls_plus, Nzwalls_plus, Nxwalls_min, Nywalls_min, Nzwalls_min, Nobst
  integer, allocatable :: ixw_p(:,:), ixw_m(:,:), iyw_p(:,:), iyw_m(:,:), izw_p(:,:)!, izw_m(:,:)
  integer, allocatable :: iobst(:,:)
  logical, allocatable :: fluid_mask(:,:,:) !< Logical which is .false. for internal building points

  !< Additional parameters
  real(field_r) :: dx_half, dy_half, Cm_xwall, Cm_ywall, Cd_xwall, Cd_ywall, Cm_zwall, Cd_zwall, z_MO

  public :: initibm, exitibm, applyibm, zerowallvelocity, fluid_mask

contains

  subroutine initibm
    use modglobal,      only : zh, zf, itot, jtot, ih, i1, jh, j1, k1, imax, jmax, kmax, cexpnr, ifnamopt, ifinput, &
                               fname_options, nsv, cu, cv, ijtot, &
                               iadv_mom,iadv_tke,iadv_thl,iadv_qt,iadv_sv,iadv_cd2, &
                               ibas_prf, &
                               dx,dy,fkar,&
                               checknamelisterror
    use modmpi,         only : myid, comm3d, mpierr, myidx, myidy, d_mpi_bcast, excjs, D_MPI_ALLREDUCE, &
                               mpi_max, mpi_sum
    use modsurface,     only : lmostlocal
    use modsubgriddata, only : lanisotrop, lsmagorinsky

    implicit none

    integer         :: i, j, k, ierr, no = 0
    integer         :: advarr(5)
    character(100)  :: readstring

    ! Temporary fields and profiles for the processing of the IBM input
    real(field_r),  allocatable :: bc_height(:,:) !< Height of immersed boundary at cell center i,j
    integer,        allocatable :: tixw_p(:,:), tixw_m(:,:), tiyw_p(:,:), tiyw_m(:,:), tizw_p(:,:)!, izw_m(:,:)

    ! Read in NAMOPTIONS parameters related to IBM
    namelist/NAMIBM/ &
      lapply_ibm, lwallheat, thlwall, thlibm, thlroof, qtibm, lpoislast, z0m_wall, z0h_wall

    call timer_tic('modibm/initibm',0)

    if( myid==0 ) then

      open (ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMIBM,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMIBM')
      write (6,NAMIBM)
      close (ifnamopt)

      ! Do some checks for conflicting settings and warn/stop further execution
      if ( lapply_ibm ) then
        if (abs(cu) > 0) stop 'Domain translation not allowed with IBM, set cu to zero'
        if (abs(cv) > 0) stop 'Domain translation not allowed with IBM, set cv to zero'

        if (ibas_prf .ne. 2) then
          ibas_prf = 2
          write (6,*) 'ibas_pr is overwritten to 2 (Boussinesq approximation with constant density)'
          write (6,*) 'height dependent density gives probles with correction of vertical advective'
          write (6,*) 'tendencies at the top of obstacles'
        end if

        !< check for use of 2nd order advection. Ideally, IBM should work with kappa advection for tracers, IMPLEMENT later.
        advarr = (/iadv_mom,iadv_tke,iadv_thl,iadv_qt,iadv_sv/)
        if (any(advarr/=iadv_cd2)) then
          write (6,*) 'Current IBM implementation only works with 2nd order advection'
          write (6,*) 'Proper check for kappa advection of scalars to be implemented'
          write (6,*) 'iostat error: ', ierr
          stop 'ERROR: Problem in namoptions NAMIBM'
        end if

        if ( lanisotrop ) then
          write (6,*) 'WARNING: you are using IBM with anisotropic grids in x,y-direction'
          write (6,*) 'while possible, this may cause unexpected results (blending effects)'
        end if

        if ( .not.(lsmagorinsky) ) then
          write (6,*) 'WARNING: subgrid TKE (e120) is not (yet) explicitly corrected for walls'
          write (6,*) 'This includes production, destruction and dissipation terms in sgs-tke budget'
        end if

        !! Check doesn't work currently, related to order in startup routine
        ! if ( .not. lmostlocal ) then
        !   stop 'ERROR: you must use local monin-obukhov with IBM'
        ! end if

      end if

    end if

    !> Broadcast IBM settings to all ranks
    call D_MPI_BCAST(lapply_ibm         ,    1, 0, comm3d, mpierr)
    call D_MPI_BCAST(lwallheat          ,    1, 0, comm3d, mpierr)  !zero (Cd->0) or nonzero heat flux from wall
    call D_MPI_BCAST(thlwall            ,    1, 0, comm3d, mpierr)
    call D_MPI_BCAST(thlibm             ,    1, 0, comm3d, mpierr)
    call D_MPI_BCAST(qtibm              ,    1, 0, comm3d, mpierr)
    call D_MPI_BCAST(thlroof            ,    1, 0, comm3d, mpierr)
    call D_MPI_BCAST(lpoislast          ,    1, 0, comm3d, mpierr)
    call D_MPI_BCAST(z0m_wall           ,    1, 0, comm3d, mpierr)
    call D_MPI_BCAST(z0h_wall           ,    1, 0, comm3d, mpierr)

    !< Step out of further subroutine when IBM is switched off
    if (.not. (lapply_ibm)) return

    ! TODO change to allow for wall dependent roughness
    ! Calculate law-of-wall coefficients for vertical walls (constant; no stability correction on vertical walls)
    dx_half = 0.5 * dx
    dy_half = 0.5 * dy
    Cm_xwall = (fkar/(log(dx_half/z0m_wall)))**2
    Cm_ywall = (fkar/(log(dy_half/z0m_wall)))**2

    ! TODO extend later as well with more options: constant flux, reactive wall temperature
    ! Set thermal boundary conditions for vertical walls to either Dirichlet (lwallheat=.true.) or zero flux
    if (lwallheat) then
      Cd_xwall = fkar**2 / log(dx_half/z0m_wall) / log(dx_half/z0h_wall)
      Cd_ywall = fkar**2 / log(dy_half/z0m_wall) / log(dy_half/z0h_wall)
    else
      Cd_xwall = 0.
      Cd_ywall = 0.
    end if

    ! Allocate fields and profiles to be used by IBM
    allocate(bc_height(itot+1,jtot+1))             ! use itot+1, jtot+1 and start writing at index = 2 to conform to field indices
    allocate(fluid_mask(2-ih:i1+ih,2-jh:j1+jh,k1))

    bc_height(:,:)  = 0.
    fluid_mask (:,:,:)    = .true.

    ! Definition of obstacles
    if (myid==0) then

      write (6,*) 'Reading inputfile ibm.inp.',cexpnr

      open (ifinput,file='ibm.inp.'//cexpnr)
        do k=1,7
          read (ifinput,'(a100)') readstring
          write (6,*) readstring
        end do

        do j=jtot+1,2,-1
          do i=2,itot+1
            read(ifinput,'(F6.1)') bc_height(i,j)
          end do
        end do

      close(ifinput)

      write(6,*) 'Succesfully read inputfile in modibm'

    end if

    !> Broadcast building heights to all ranks
    call D_MPI_BCAST(bc_height, (itot+1)*(jtot+1), 0, comm3d, mpierr)

    !> Determine obstacle cells. Use obstacle height is above midpoint of vertical cell (= full levels). Corresponds to >50% of cell being filled.
    Nobst = 0
    do i=2,i1
      do j=2,j1
         do k=1,kmax
            if (zf(k) <= bc_height(i+myidx*imax,j+myidy*jmax)) then
              fluid_mask(i,j,k) = .false.     ! Set grid cell to obstacle
              Nobst             = Nobst + 1   ! Increase counter of obstacles by one
              ! ksfc  (i,j)   = k + 1         ! Store index of top roof (at half level; w-flux)
           end if
        end do
      end do
    end do

    !> Set ghost cells for fluid mask
    call excjs(fluid_mask  , 2,i1,2,j1,1,k1,ih,jh)

    if( myid == 0 ) then
      write(6,*) 'Start determination of wall positions'
    end if

    !> Identify sidewalls based on fluid_mask-field
    !!  u-positions (with index i) are to the left of grid center (with index i)
    !!  walls with normal in minus x-direction will therefore be at same index, walls with normal in plus direction at i+1
    !!  see visual example below
    !!
    !! building positions   :   0     X     X     X     0
    !! grid centered index  :  i-2   i-1    i    i+1   i+2
    !! fluid_mask           :   F     T     T     T     F
    !! xwall min yes/no     :   F     T     F     F     F
    !! xwall plus yes/no    :   F     F     F     F     T

    !> Preset (temporary) arrays for wall indices (Nobst [+i1 (or j1)*k1] provides upper bound)
    allocate(iobst (Nobst,3))                                                                                             !< can be immediately set to correct size
    allocate(tixw_p(Nobst+i1*k1,3), tixw_m(Nobst+i1*k1,3), tiyw_p(Nobst+j1*k1,3), tiyw_m(Nobst+j1*k1,3), tizw_p(Nobst,3)) !< for x- and y-walls in positive and negative directions

    Nxwalls_plus = 0; Nxwalls_min = 0; Nywalls_plus = 0; Nywalls_min = 0; Nzwalls_plus = 0
    !! Potential rework to get .not. statements out -> fill plus walls first
    no = 0;
    do k=1,kmax
      do j=2,j1
        do i=2,i1

          if ( .not.(fluid_mask(i,j,k)) ) then
            no          = no + 1                  ! local counter for obstacle points
            iobst(no,1) = i
            iobst(no,2) = j
            iobst(no,3) = k
          end if

          !> Because we always compare 1 grid index lower, no double counts of walls occur
          if (fluid_mask(i,j,k) .neqv. fluid_mask(i-1,j,k)) then  ! signals wall in x-direction
            if ( .not.(fluid_mask(i,j,k)) ) then                  ! wall with normal in minus x-direction
              Nxwalls_min            = Nxwalls_min  + 1
              tixw_m(Nxwalls_min, 1) = i
              tixw_m(Nxwalls_min, 2) = j
              tixw_m(Nxwalls_min, 3) = k
            else                                                  ! normal in plus x-direction
              Nxwalls_plus           = Nxwalls_plus + 1
              tixw_p(Nxwalls_plus,1) = i
              tixw_p(Nxwalls_plus,2) = j
              tixw_p(Nxwalls_plus,3) = k
            end if
          end if

          if (fluid_mask(i,j,k) .neqv. fluid_mask(i,j-1,k)) then  ! signals wall in y-direction
            if ( .not.(fluid_mask(i,j,k)) ) then                  ! wall with normal in minus y-direction
              Nywalls_min            = Nywalls_min  + 1
              tiyw_m(Nywalls_min, 1) = i
              tiyw_m(Nywalls_min, 2) = j
              tiyw_m(Nywalls_min, 3) = k
            else                                                  ! normal in plus y-direction
              Nywalls_plus           = Nywalls_plus + 1
              tiyw_p(Nywalls_plus,1) = i
              tiyw_p(Nywalls_plus,2) = j
              tiyw_p(Nywalls_plus,3) = k
            end if
          end if

          if (k == 1)  cycle                                      ! if near surface, go to next loop iteration (as fluid_mask(i,j,0) doesn't exist)

          if (fluid_mask(i,j,k) .neqv. fluid_mask(i,j,k-1)) then  ! signals wall in z-direction
            ! if ( fluid_mask(i,j,k) ) then                       ! wall with normal in plus y-direction
            Nzwalls_plus            = Nzwalls_plus  + 1
            tizw_p(Nzwalls_plus, 1) = i
            tizw_p(Nzwalls_plus, 2) = j
            tizw_p(Nzwalls_plus, 3) = k
            !! Current impementation does not allow overhanging obstacles (i.e. walls in negative z)
            ! else                                                ! normal in minus z-direction
            !   Nzwalls_min           = Nzwalls_min + 1
            !   tiz_m(Nzwalls_min,1) = i
            !   tizw_m(Nzwalls_min,2) = j
            !   tizw_m(Nzwalls_min,3) = k
            ! end if
          end if

        end do
      end do
    end do

    !> Should be same.
    if (no /= Nobst) then
      write (6,*) 'Number of identified obstacle points during'
      write (6,*) 'wall determination and prior stage do not match!'
      stop 'ERROR: Problem in subroutine initibm'
    end if

    !> Copy temporary arrays with indices into final ones
    allocate(ixw_p(Nxwalls_plus,3))
    allocate(ixw_m(Nxwalls_min, 3))
    allocate(iyw_p(Nywalls_plus,3))
    allocate(iyw_m(Nywalls_min ,3))
    allocate(izw_p(Nzwalls_plus,3))

    ixw_p(1:Nxwalls_plus,:) = tixw_p(1:Nxwalls_plus,:)
    ixw_m(1:Nxwalls_min ,:) = tixw_m(1:Nxwalls_min ,:)
    iyw_p(1:Nywalls_plus,:) = tiyw_p(1:Nywalls_plus,:)
    iyw_m(1:Nywalls_min ,:) = tiyw_m(1:Nywalls_min ,:)
    izw_p(1:Nzwalls_plus,:) = tizw_p(1:Nzwalls_plus,:)

    !> Free memory of temperorary fields and arrays
    deallocate(bc_height)
    deallocate(tixw_p, tixw_m, tiyw_p, tiyw_m, tizw_p)

    !> Finally, copy data to GPU
    !$acc enter data copyin(fluid_mask, iobst, ixw_p, ixw_m, iyw_p, iyw_m, izw_p)

    call timer_toc('modibm/initibm')

    return
  end subroutine initibm

  subroutine exitibm
    !< Step out of further subroutine when IBM is switched off
    if (.not. (lapply_ibm)) return

    !$acc exit data delete(fluid_mask, iobst, ixw_p, ixw_m, iyw_p, iyw_m, izw_p)

    deallocate(iobst)
    deallocate(ixw_p)
    deallocate(ixw_m)
    deallocate(iyw_p)
    deallocate(iyw_m)
    deallocate(izw_p)

    deallocate(fluid_mask)

    return
  end subroutine exitibm

  subroutine applyibm
    use modfields,      only : um, vm, wm, thlm, qtm, e12m, svm, &
                               u0, v0, w0, thl0, qt0, e120, sv0, &
                               up, vp, wp, thlp, qtp, e12p, svp, &
                               thl0av, qt0av, rhobf, rhobh
    use modglobal,      only : rk3step, kmax, i1, j1, k1, ih, jh, rdt, timee, dx, dy, dx2i, dy2i, dzh, dzhi, dzf, dzfi, zf, zh, nsv, e12min, fkar
    use modsurface,     only : lneutral, lmostlocal
    use modsubgriddata, only : ekm, ekh
    use modmpi,         only : excjs

    implicit none

    integer           :: i, j, k, nn, nc
    real(field_r)     :: rk3coef, rk3coefi
    real(field_r)     :: emmo, empo, emom, emop, eomm
    real(field_r)     :: u_at_v_min, u_at_v_plus, v_at_u_min, v_at_u_plus
    real(field_r)     :: w_at_v_min, w_at_v_plus, w_at_u_min, w_at_u_plus
    real(field_r)     :: u_at_w_min, u_at_w_plus, v_at_w_min, v_at_w_plus
    real(field_r)     :: uspeed, z_MO
    real(field_r)     :: tau_vu_plus, tau_vu_min, tau_wu_min, tau_wu_plus, tau_uv_min, tau_uv_plus, tau_wv_min, tau_wv_plus
    real              :: Lob
    
    if (.not. lapply_ibm) return

    call timer_tic('modibm/applyibm',0)
    ! Commented lines below taken from Stephan's version (needed for reproduction of his results)
    ! thlibm = thl0av(1) !assumes inside air has the same temperature as the air outside, at near-surface level
    ! qtibm  = qt0av (1)

    rk3coef = rdt / (4. - dble(rk3step))
    rk3coefi = 1. / rk3coef

    !> Start by setting enforcing cyclic boundary conditions for everything, as correction may go over MPI bound
    ! call excjs( u0  , 2,i1,2,j1,1,k1,ih,jh)
    ! call excjs( v0  , 2,i1,2,j1,1,k1,ih,jh)
    ! call excjs( w0  , 2,i1,2,j1,1,k1,ih,jh)
    call excjs( e120  , 2,i1,2,j1,1,k1,ih,jh)
    call excjs( up    , 2,i1,2,j1,1,k1,ih,jh)
    call excjs( vp    , 2,i1,2,j1,1,k1,ih,jh)
    call excjs( wp    , 2,i1,2,j1,1,k1,ih,jh)
    call excjs( e12p  , 2,i1,2,j1,1,k1,ih,jh)
    call excjs( thlp  , 2,i1,2,j1,1,k1,ih,jh)
    call excjs( qtp   , 2,i1,2,j1,1,k1,ih,jh)

    !> Set tendencies inside obstacles (i.e., correct for any drift from previous integration step)
    !$acc parallel loop gang vector default(present)
    do nn = 1,Nobst
      i = iobst(nn,1)
      j = iobst(nn,2)
      k = iobst(nn,3)

      ! Correction of velocities also corrects walls w. normals in negative x,y-directions (due to staggered grid arrangement)
      ! still do these again later for consistency with respect to positive x,y-directions
      up(i,j,k)   = -u0(i,j,k)*rk3coefi
      vp(i,j,k)   = -v0(i,j,k)*rk3coefi
      wp(i,j,k)   = -w0(i,j,k)*rk3coefi
      thlp(i,j,k) = (thlibm - thl0(i,j,k))*rk3coefi
      qtp (i,j,k) = (qtibm  - qt0(i,j,k) )*rk3coefi
      e12p(i,j,k) = (e12min - e120(i,j,k))*rk3coefi   !< Maintain e12min to prevent ekm/ekh going to NaN
                                                      !< Maybe better to explicitly set eddy diffusivities to tiny value
      do nc=1,nsv
          svp(i,j,k,nc) = - sv0(i,j,k,nc)*rk3coefi
      end do

    end do

    !> All corrections consist of 2 parts: cancel the "wrong" diffusion imposed by modsubgrid at fluid points at the wall, and add wall friction/flux.
    !> In this framework, we assume no correction for advection is needed, as velocities at walls should already be close to zero (may better enforce later)

    !> Correct tendencies for walls in positive z-direction (only works when k>1, which should be the case for vertical walls [see initibm])
    !$acc parallel loop gang vector default(present)
    do nn = 1,Nzwalls_plus
      i = izw_p(nn,1)
      j = izw_p(nn,2)
      k = izw_p(nn,3)

      !> First set the normal velocity at the wall (correct for any drift)
      wp(i,j,k)   = -w0(i,j,k)*rk3coefi

      !> Remove "wrong" diffusive tendencies
      ! for u(i,j,k):
      emom = ( dzf(k-1) * ( ekm(i,j,k  ) + ekm(i+1,j,k  ) )  + &
            dzf(k  ) * ( ekm(i,j,k-1) + ekm(i+1,j,k-1) ) ) * &
            ( 0.25_field_r * dzhi(k) )
      up(i,j,k) = up(i,j,k) + 0.5_field_r * rhobh(k)/rhobf(k) * emom * ( (u0(i,j,k) - u0(i,j,k-1)) *dzhi(k)) *dzfi(k)

      ! for u(i+1,j,k):
      emom = ( dzf(k-1) * ( ekm(i,j,k  ) + ekm(i-1,j,k  ) )  + &
            dzf(k  ) * ( ekm(i,j,k-1) + ekm(i-1,j,k-1) ) ) * &
            ( 0.25_field_r * dzhi(k) )
      up(i+1,j,k) = up(i+1,j,k) + 0.5_field_r * rhobh(k)/rhobf(k) * emom * ( (u0(i+1,j,k) - u0(i+1,j,k-1)) *dzhi(k)) *dzfi(k)

      ! for v(i,j,k):
      eomm = ( dzf(k-1) * ( ekm(i,j,k  ) + ekm(i,j-1,k ) ) + &
            dzf(k  ) * ( ekm(i,j,k-1) + ekm(i,j-1,k-1) ) ) * &
            ( .25_field_r * dzhi(k) )
      vp(i,j,k) = vp(i,j,k) + 0.5_field_r * rhobh(k)/rhobf(k) * eomm * ( (v0(i,j,k) - v0(i,j,k-1)) * dzhi(k)) * dzfi(k)

      ! for v(i,j+1,k):
      eomm = ( dzf(k-1) * ( ekm(i,j,k  ) + ekm(i,j+1,k ) ) + &
            dzf(k  ) * ( ekm(i,j,k-1) + ekm(i,j+1,k-1) ) ) * &
            ( .25_field_r * dzhi(k) )
      vp(i,j+1,k) = vp(i,j+1,k) + 0.5_field_r * rhobh(k)/rhobf(k) * eomm * ( (v0(i,j+1,k)-v0(i,j+1,k-1)) * dzhi(k)) * dzfi(k)

      !> Include correct Monin-Obukhov wall drag (with stability correction)
      uspeed  = 0.5_field_r * sqrt( ( u0(i,j,k) + u0(i+1,j,k) )**2 + ( v0(i,j,k) + v0(i,j+1,k) )**2 )
      z_MO    = zf(k) - zh(k)

      !> If locally no wind, automatically no flux, so default drag coefficient to zero REWORK THIS!!!
      if (uspeed < 0.01) then
        Cm_zwall = 0
        Cd_zwall = 0
      else if (lneutral) then
        Cm_zwall = fkar**2 / (log(z_MO / z0m_wall))** 2
        Cd_zwall = fkar**2 / (log(z_MO / z0m_wall)) / (log(z_MO / z0h_wall))
      else
        Lob      = getobl_local(uspeed, thl0(i,j,k), qt0(i,j,k), thlroof, qtroof, z_MO, z0m_wall, z0h_wall)
        Cm_zwall = fkar**2 / (log(z_MO / z0m_wall) - psim(z_MO / Lob) + psim(z0m_wall / Lob))** 2
        Cd_zwall = fkar**2 / (log(z_MO / z0m_wall) - psim(z_MO / Lob) + psim(z0m_wall / Lob)) / (log(z_MO / z0h_wall) - psih(z_MO / Lob) + psih(z0h_wall / Lob))
      end if

      up(i  ,j,k) = up(i  ,j,k) - 0.25_field_r * rhobh(k)/rhobf(k) * Cm_zwall * ( u0(i,j,k) + u0(i+1,j,k) ) * uspeed * dzfi(k)
      up(i+1,j,k) = up(i+1,j,k) - 0.25_field_r * rhobh(k)/rhobf(k) * Cm_zwall * ( u0(i,j,k) + v0(i+1,j,k) ) * uspeed * dzfi(k)
      vp(i,j  ,k) = vp(i,j  ,k) - 0.25_field_r * rhobh(k)/rhobf(k) * Cm_zwall * ( v0(i,j,k) + v0(i,j+1,k) ) * uspeed * dzfi(k)
      vp(i,j+1,k) = vp(i,j+1,k) - 0.25_field_r * rhobh(k)/rhobf(k) * Cm_zwall * ( v0(i,j,k) + v0(i,j+1,k) ) * uspeed * dzfi(k)

      ! SvdL: tentative fix for vertical diffusion of temperature over z-walls. Do theck later. Also here, 0.5 comes from interpolation of ekm.
      thlp(i,j,k) = thlp(i,j,k) + 0.5_field_r * rhobh(k)/rhobf(k) * ( ( ( dzf(k-1) * ekm(i,j,k) ) + ( dzf(k ) * ekm(i,j,k-1)) )* dzhi(k)  ) * ( thl0(i,j,k) - thl0(i-1,j,k) ) * dzhi(k) * dzfi(k)
      
      !> not sure yet how to properly include zero flux conditions (... copy relevant isurf parts from modsurface.f90 here)
      ! if flux = constant, it is just the flux to feed into here. Yet, roof temperatures should still be diagnosed somehow (requiring energy balance)
      ! doubting wheter rhobh(k)/rhobf(k) are truly correct here, i.e., at right positions..
      thlp(i,j,k) = thlp(i,j,k) - rhobh(k)/rhobf(k) * Cd_zwall * ( thl0(i,j,k) - thlroof ) * uspeed * dzfi(k)
      qtp(i,j,k)  =  qtp(i,j,k) - rhobh(k)/rhobf(k) * Cd_zwall * ( qt0(i,j,k) - qtroof )   * uspeed * dzfi(k)
      ! no flux of e12 from wall into flow

      !> not sure how to include constant, i.e., fixed emission fluxes for tracers yet (... confer with Caspar)
      do nc=1,nsv
        svp(i,j,k,nc) = svp(i,j,k,nc) - 0 !rhobh(k)/rhobf(k) * Cd_zwall * ( sv0(i,j,k,nc) -svroof(nc) ) * uspeed * dzfi(k) !< currently set to zero
      end do

      ! slight differences in implementation w.r.t. regular diffusion (as the interpolations there are different)
      ! old option below was interpolating always to momentum location.. yet somehow seems less consistent with general framework
          ! u_at_v_min  = 0.25_field_r * ( u0(i,j,k) + u0(i+1,j,k) + u0(i,j-1,k) + u0(i+1,j-1,k) )
          ! u_at_v_plus = 0.25_field_r * ( u0(i,j,k) + u0(i+1,j,k) + u0(i,j+1,k) + u0(i+1,j+1,k) )
          ! v_at_u_min  = 0.25_field_r * ( v0(i,j,k) + v0(i-1,j,k) + v0(i,j+1,k) + v0(i-1,j+1,k) )
          ! v_at_u_plus = 0.25_field_r * ( v0(i,j,k) + v0(i+1,j,k) + v0(i,j+1,k) + v0(i+1,j+1,k) )
          ! uspeed      = 0.5_field_r * sqrt( ( u0(i,j,k) + u0(i+1,j,k) )**2 + ( v0(i,j,k) + v0(i,j+1,k) )**2 )

          ! up(i  ,j,k) = up(i  ,j,k) - 0.5_field_r * rhobh(k)/rhobf(k) * Cm_zwall * u0(i  ,j,k) * sqrt( u0(i  ,j,k)**2 +  v_at_u_min**2 ) * dzfi(k)
          ! up(i+1,j,k) = up(i+1,j,k) - 0.5_field_r * rhobh(k)/rhobf(k) * Cm_zwall * u0(i+1,j,k) * sqrt( u0(i+1,j,k)**2 + v_at_u_plus**2 ) * dzfi(k)
          ! vp(i,j  ,k) = vp(i,j  ,k) - 0.5_field_r * rhobh(k)/rhobf(k) * Cm_zwall * v0(i,j  ,k) * sqrt( v0(i,j  ,k)**2 +  u_at_v_min**2 ) * dzfi(k)
          ! vp(i,j+1,k) = vp(i,j+1,k) - 0.5_field_r * rhobh(k)/rhobf(k) * Cm_zwall * v0(i,j+1,k) * sqrt( v0(i,j+1,k)**2 + u_at_v_plus**2 ) * dzfi(k)

    end do

    !> Correct tendencies for walls in positive x-direction
    !$acc parallel loop gang vector default(present)
    do nn = 1,Nxwalls_plus
      i = ixw_p(nn,1)
      j = ixw_p(nn,2)
      k = ixw_p(nn,3)

      !> First set the normal velocity at the wall (correct for any drift)
      up(i,j,k)   = -u0(i,j,k)*rk3coefi

      !> Remove "wrong" diffusive tendencies and replace with wall drag
      ! for v(i,j,k):
      emmo = 0.25_field_r * ( ekm(i,j,k) + ekm(i,j-1,k) + ekm(i-1,j,k) + ekm(i-1,j-1,k) )
      w_at_v_plus = 0.25_field_r * ( w0(i,j,k) + w0(i,j,k+1) + w0(i,j-1,k) + w0(i,j-1,k+1) )  !at v(i,j,k)
      tau_vu_plus = log_wallaw(v0(i,j,k), w_at_v_plus, Cm_xwall)

      vp(i,j  ,k) = vp(i,j  ,k) + 0.5_field_r * emmo * ( (v0(i,j,k) - v0(i-1,j,k) ) / dx) / dx - 0.5_field_r * tau_vu_plus / dx ! factor 0.5 originates to avoid double correction
      ! NOTE: if x-plane next to current one is also a wall, part will be corrected via v(i,j+1,k) in that plane. if plane is not a wall, we do want to allow for partial contribution to diffusive flux

      ! for v(i,j+1,k):
      empo = 0.25_field_r * ( ekm(i,j,k) + ekm(i,j+1,k) + ekm(i-1,j,k) + ekm(i-1,j+1,k) )
      w_at_v_plus = 0.25_field_r * (w0(i,j+1,k) + w0(i,j+1,k+1) + w0(i,j,k) + w0(i,j,k+1) )
      tau_vu_plus = log_wallaw(v0(i,j+1,k), w_at_v_plus, Cm_xwall)

      vp(i,j+1,k) = vp(i,j+1,k) + 0.5_field_r * empo * ( (v0(i,j+1,k) - v0(i-1,j+1,k) ) / dx) / dx - 0.5_field_r * tau_vu_plus / dx

      ! for w(i,j,k):
      if (k /= 1) then ! not correctable when at surface (k = 1)
        emom = ( dzf(k-1) * ( ekm(i,j,k)  + ekm(i-1,j,k)  )  + &
                dzf(k)  * ( ekm(i,j,k-1) + ekm(i-1,j,k-1) ) ) / &
                ( 4.0_field_r * dzh(k) )
        v_at_w_plus = 0.25_field_r * ( v0(i,j,k-1) + v0(i,j,k) + v0(i,j+1,k-1) + v0(i,j+1,k) )
        tau_wu_plus = log_wallaw(w0(i  ,j,k),v_at_w_plus,Cm_xwall)

        wp(i,j,k  ) = wp(i,j,k  ) + 0.5_field_r * emom * ( (w0(i,j,k  ) - w0(i-1,j,  k)) / dx ) / dx - 0.5_field_r * tau_wu_plus / dx

      end if

      ! for w(i,j,k+1):
      emop = ( dzf(k) * ( ekm(i,j,k+1)  + ekm(i-1,j,k+1)  )  + &
              dzf(k+1)  * ( ekm(i,j,k) + ekm(i-1,j,k) ) ) / &
              ( 4.0_field_r * dzh(k+1) )
      v_at_w_plus = 0.25_field_r * ( v0(i,j,k) + v0(i,j,k+1) + v0(i,j+1,k) + v0(i,j+1,k+1) )
      tau_wu_plus = log_wallaw(w0(i  ,j,k+1), v_at_w_plus, Cm_xwall)

      wp(i,j,k+1) = wp(i,j,k+1) + 0.5_field_r * emop * ( (w0(i,j,k+1) - w0(i-1,j,k+1) ) / dx ) / dx - 0.5_field_r * tau_wu_plus / dx

      !> Enforcing zero flux by correction of tendencies of temperature, moisture and other scalars (by negating diffusion term)

      ! scalars at (i,j,k) are "to the right" of the x-walls, so correct on s(i,j,k)
      ! "+" because we reflect it back, and 0.5_field_r comes from interpolation of ekh
      thlp(i,j,k) = thlp(i,j,k) + 0.5_field_r * ( ekh(i,j,k) + ekh(i-1,j,k) ) * ( thl0(i,j,k) - thl0(i-1,j,k) ) * dx2i
      qtp(i,j,k)  =  qtp(i,j,k) + 0.5_field_r * ( ekh(i,j,k) + ekh(i-1,j,k) ) * (  qt0(i,j,k) -  qt0(i-1,j,k) ) * dx2i

      do nc=1,nsv
        svp(i,j,k,nc) = svp(i,j,k,nc) + 0.5_field_r * ( ekh(i,j,k) + ekh(i-1,j,k) ) * ( sv0(i,j,k,nc) - sv0(i-1,j,k,nc) ) * dx2i
      end do
      !call xwalle12(i,j,k) ! correction is ignored assuming u,v,w,subgrid TKE are near zero inside buildings

      !> Finally, set correct heat flux from wall to fluid
      uspeed = 0.5_field_r * ( ( v0(i,j,k) + v0(i,j+1,k) )**2 + ( w0(i,j,k) + w0(i,j,k+1) )**2 )**0.5_field_r
      thlp(i,j,k) = thlp(i,j,k) + Cd_xwall * uspeed * (thlwall - thl0(i,j,k)) / dx

    end do

    !> Correct tendencies for walls in negative x-direction
    !$acc parallel loop gang vector default(present)
    do nn = 1,Nxwalls_min
      i = ixw_m(nn,1)
      j = ixw_m(nn,2)
      k = ixw_m(nn,3)

      ! First set the normal velocity at the wall (correct for any drift)
      ! note: these should actually already be corrected in loop over iobst (due to staggered grid arrangement)
      up(i,j,k)   = -u0(i,j,k)*rk3coefi

      !> Remove "wrong" diffusive tendencies and replace with wall drag
      ! for v(i-1,j,k):
      emmo = 0.25_field_r * ( ekm(i,j,k) + ekm(i,j-1,k) + ekm(i-1,j,k) + ekm(i-1,j-1,k) )
      w_at_v_min  = 0.25_field_r * ( w0(i-1,j,k) + w0(i-1,j,k+1) + w0(i-1,j-1,k) + w0(i-1,j-1,k+1) )  !at v(i-1,j,k)
      tau_vu_min = log_wallaw(v0(i-1,j,k), w_at_v_min , Cm_xwall) !if v0 > 0, tau > 0, minus sign in tendency enforces opposing friction

      vp(i-1,j  ,k) = vp(i-1,j  ,k) - 0.5_field_r * emmo * ( (v0(i,j,k) - v0(i-1,j,k) ) / dx) / dx  - 0.5_field_r * tau_vu_min / dx ! factor 0.5 originates to avoid double correction

      ! for v(i-1,j+1,k):
      empo = 0.25_field_r * ( ekm(i,j,k) + ekm(i,j+1,k) + ekm(i-1,j,k) + ekm(i-1,j+1,k) )
      w_at_v_min  = 0.25_field_r * ( w0(i-1,j+1,k) + w0(i-1,j+1,k+1) + w0(i-1,j,k) + w0(i-1,j,k+1) )
      tau_vu_min = log_wallaw(v0(i-1,j+1,k), w_at_v_min , Cm_xwall)

      vp(i-1,j+1,k) = vp(i-1,j+1,k) - 0.5_field_r * empo * ( (v0(i,j+1,k) - v0(i-1,j+1,k) ) / dx) / dx - 0.5_field_r * tau_vu_min / dx

      ! for w(i-1,j,k):
      if (k /= 1) then ! not correctable when at surface (k = 1)
        emom = ( dzf(k-1) * ( ekm(i,j,k)  + ekm(i-1,j,k)  )  + &
                dzf(k)  * ( ekm(i,j,k-1) + ekm(i-1,j,k-1) ) ) / &
                ( 4.0_field_r * dzh(k) )
        v_at_w_min = 0.25_field_r * (v0(i-1,j,k-1)+v0(i-1,j,k)+v0(i-1,j+1,k-1)+v0(i-1,j+1,k) )
        tau_wu_min = log_wallaw(w0(i-1,j,k), v_at_w_min , Cm_xwall)

        wp(i-1,j,k  ) = wp(i-1,j,k  ) - 0.5_field_r * emom * ( (w0(i,j,k  ) - w0(i-1,j,  k)) / dx ) / dx - 0.5_field_r * tau_wu_min / dx
      end if

      ! for w(i-1,j,k+1):
      emop = ( dzf(k) * ( ekm(i,j,k+1)  + ekm(i-1,j,k+1)  )  + &
              dzf(k+1)  * ( ekm(i,j,k) + ekm(i-1,j,k) ) ) / &
              ( 4.0_field_r * dzh(k+1) )
      v_at_w_min  = 0.25_field_r * (v0(i-1,j,k)+v0(i-1,j,k+1)+v0(i-1,j+1,k)+v0(i-1,j+1,k+1) )
      tau_wu_min = log_wallaw(w0(i-1,j,k+1), v_at_w_min , Cm_xwall)

      wp(i-1,j,k+1) = wp(i-1,j,k+1) - 0.5_field_r * emop * ( (w0(i,j,k+1) - w0(i-1,j,k+1) ) / dx ) / dx - 0.5_field_r * tau_wu_min / dx

      !> Enforcing zero flux by correction of tendencies of temperature, moisture and other scalars (by negating diffusion term)

      ! scalars at (i,j,k) are "to the right" of the x-walls, so correct on s(i-1,j,k) which is on the left for minus x-walls
      ! "-" because we reflect it back, and 0.5_field_r comes from interpolation of ekh
      thlp(i-1,j,k) = thlp(i-1,j,k) - 0.5_field_r * ( ekh(i,j,k) + ekh(i-1,j,k) ) * ( thl0(i,j,k) - thl0(i-1,j,k) ) * dx2i
      qtp(i-1,j,k)  =  qtp(i-1,j,k) - 0.5_field_r * ( ekh(i,j,k) + ekh(i-1,j,k) ) * (  qt0(i,j,k) -  qt0(i-1,j,k) ) * dx2i

      do nc=1,nsv
        svp(i,j,k,nc) = svp(i,j,k,nc) - 0.5_field_r * ( ekh(i,j,k) + ekh(i-1,j,k) ) * ( sv0(i,j,k,nc) - sv0(i-1,j,k,nc) ) * dx2i
      end do

      !> Finally, set correct heat flux from wall to fluid
      uspeed = 0.5_field_r * ( ( v0(i-1,j,k) + v0(i-1,j+1,k) )**2 + ( w0(i-1,j,k) + w0(i-1,j,k+1) )**2 )**0.5_field_r
      thlp(i-1,j,k) = thlp(i-1,j,k) + Cd_xwall * uspeed * (thlwall - thl0(i-1,j,k)) / dx

    end do

    !> Correct tendencies for walls in positive y-direction
    !$acc parallel loop gang vector default(present)
    do nn = 1,Nywalls_plus
      i = iyw_p(nn,1)
      j = iyw_p(nn,2)
      k = iyw_p(nn,3)

      ! First set the normal velocity at the wall (correct for any drift)
      vp(i,j,k)   = -v0(i,j,k)*rk3coefi

      !> Remove "wrong" diffusive tendencies and replace with wall drag
      ! for u(i,j,k):
      emmo = 0.25_field_r * ( ekm(i,j,k) + ekm(i,j-1,k) + ekm(i-1,j,k) + ekm(i-1,j-1,k) )
      w_at_u_plus = 0.25_field_r * ( w0(i,j,k) + w0(i,j,k+1) + w0(i-1,j,k) + w0(i-1,j,k+1) )  !at u(i,j,k)
      tau_uv_plus = log_wallaw(u0(i,j,k), w_at_u_plus, Cm_ywall)

      up(i,j  ,k) = up(i,j  ,k) + 0.5_field_r * emmo * ( (u0(i,j,k) - u0(i,j-1,k) ) / dy) / dy - 0.5_field_r * tau_uv_plus / dy

      ! for u(i+1,j,k):
      empo = 0.25_field_r * ( ekm(i,j,k) + ekm(i+1,j+1,k) + ekm(i,j,k) + ekm(i,j+1,k) )
      w_at_u_plus = 0.25_field_r * (w0(i+1,j,k) + w0(i+1,j,k+1) + w0(i,j,k) + w0(i,j,k+1) )
      tau_uv_plus = log_wallaw(u0(i+1,j,k), w_at_u_plus, Cm_ywall)

      up(i+1,j,k) = up(i+1,j,k) + 0.5_field_r * empo * ( (u0(i+1,j,k) - u0(i,j,k) ) / dy) / dy - 0.5_field_r * tau_uv_plus / dy

      ! for w(i,j,k):
      if (k /= 1) then ! not correctable when at surface (k = 1)
        emom = ( dzf(k-1) * ( ekm(i,j,k)  + ekm(i,j-1,k)  )  + &
                dzf(k)  * ( ekm(i,j,k-1) + ekm(i,j-1,k-1) ) ) / &
                ( 4.0_field_r * dzh(k) )
        u_at_w_plus = 0.25_field_r * ( u0(i,j,k-1) + u0(i,j,k) + u0(i+1,j,k-1) + u0(i+1,j,k) )
        tau_wv_plus = log_wallaw(w0(i  ,j,k), u_at_w_plus, Cm_ywall)

        wp(i,j,k  ) = wp(i,j,k  ) + 0.5_field_r * emom * ( (w0(i,j,k  ) - w0(i,j-1,  k)) / dy ) / dy - 0.5_field_r * tau_wv_plus / dy

      end if

      ! for w(i,j,k+1):
      emop = ( dzf(k) * ( ekm(i,j,k+1)  + ekm(i,j-1,k+1)  )  + &
              dzf(k+1)  * ( ekm(i,j,k) + ekm(i,j-1,k) ) ) / &
              ( 4.0_field_r * dzh(k+1) )
      u_at_w_plus = 0.25_field_r * ( u0(i,j,k) + u0(i,j,k+1) + u0(i+1,j,k) + u0(i+1,j,k+1) )
      tau_wv_plus = log_wallaw(w0(i  ,j,k+1), u_at_w_plus, Cm_ywall)

      wp(i,j,k+1) = wp(i,j,k+1) + 0.5_field_r * emop * ( (w0(i,j,k+1) - w0(i,j-1,k+1) ) / dy ) / dy - 0.5_field_r * tau_wv_plus / dy

      !> Enforcing zero flux by correction of tendencies of temperature, moisture and other scalars (by negating diffusion term)

      ! scalars at (i,j,k) are "to the right" of the y-walls, so correct on s(i,j,k)
      ! "+" because we reflect it back, and 0.5_field_r comes from interpolation of ekh
      thlp(i,j,k) = thlp(i,j,k) + 0.5_field_r * ( ekh(i,j,k) + ekh(i,j-1,k) ) * ( thl0(i,j,k) - thl0(i,j-1,k) ) * dy2i
      qtp(i,j,k)  =  qtp(i,j,k) + 0.5_field_r * ( ekh(i,j,k) + ekh(i,j-1,k) ) * (  qt0(i,j,k) -  qt0(i,j-1,k) ) * dy2i

      do nc=1,nsv
        svp(i,j,k,nc) = svp(i,j,k,nc) + 0.5_field_r * ( ekh(i,j,k) + ekh(i,j-1,k) ) * ( sv0(i,j,k,nc) - sv0(i,j-1,k,nc) ) * dy2i
      end do
      !call xwalle12(i,j,k) ! correction is ignored assuming u,v,w,subgrid TKE are near zero inside buildings

      !> Finally, set correct heat flux from wall to fluid
      uspeed = 0.5_field_r * ( ( u0(i,j,k) + u0(i+1,j,k) )**2 + ( w0(i,j,k) + w0(i,j,k+1) )**2 )**0.5_field_r
      thlp(i,j,k) = thlp(i,j,k) + Cd_ywall * uspeed * (thlwall - thl0(i,j,k)) / dy

    end do

    !> Correct tendencies for walls in negative y-direction
    !$acc parallel loop gang vector default(present)
    do nn = 1,Nywalls_min
      i = iyw_m(nn,1)
      j = iyw_m(nn,2)
      k = iyw_m(nn,3)

      ! First set the normal velocity at the wall (correct for any drift)
      ! note: these should actually already be corrected in loop over iobst (due to staggered grid arrangement)
      vp(i,j,k)   = -v0(i,j,k)*rk3coefi

      !> Remove "wrong" diffusive tendencies and replace with wall drag
      ! for u(i,j-1,k):
      emmo = 0.25_field_r * ( ekm(i,j,k) + ekm(i,j-1,k) + ekm(i-1,j,k) + ekm(i-1,j-1,k) )
      w_at_u_min  = 0.25_field_r * ( w0(i,j-1,k) + w0(i,j-1,k+1) + w0(i-1,j-1,k) + w0(i-1,j-1,k+1) )  !at u(i,j-1,k)
      tau_uv_min = log_wallaw(u0(i,j-1,k), w_at_u_min , Cm_ywall) !if v0 > 0, tau > 0, minus sign in tendency enforces opposing friction

      up(i,j-1,k) = up(i,j-1  ,k) - 0.5_field_r * emmo * ( (u0(i,j,k) - u0(i,j-1,k) ) / dy) / dy  - 0.5_field_r * tau_uv_min / dy ! factor 0.5 originates to avoid double correction

      ! for u(i+1,j-1,k):
      empo = 0.25_field_r * ( ekm(i,j,k) + ekm(i,j-1,k) + ekm(i+1,j,k) + ekm(i+1,j-1,k) )
      w_at_u_min = 0.25_field_r * ( w0(i+1,j-1,k) + w0(i+1,j-1,k+1) + w0(i,j-1,k) + w0(i,j-1,k+1) )
      tau_uv_min = log_wallaw(u0(i+1,j-1,k), w_at_u_min , Cm_ywall)

      up(i+1,j-1,k) = up(i+1,j-1,k) - 0.5_field_r * empo * ( (u0(i+1,j,k) - u0(i+1,j-1,k) ) / dy) / dy - 0.5_field_r * tau_uv_min / dy

      ! for w(i,j-1,k):
      if (k /= 1) then ! not correctable when at surface (k = 1)
        emom = ( dzf(k-1) * ( ekm(i,j,k)  + ekm(i,j-1,k)  )  + &
                dzf(k)  * ( ekm(i,j,k-1) + ekm(i,j-1,k-1) ) ) / &
                ( 4.0_field_r * dzh(k) )
        u_at_w_min = 0.25_field_r * ( u0(i,j-1,k-1) + u0(i,j-1,k) + u0(i+1,j-1,k-1) + u0(i+1,j-1,k) )
        tau_wv_min = log_wallaw(w0(i,j-1,k), v_at_w_min , Cm_ywall)

        wp(i,j-1,k  ) = wp(i,j-1,k  ) - 0.5_field_r * emom * ( (w0(i,j,k  ) - w0(i,j-1,  k)) / dy ) / dy - 0.5_field_r * tau_wv_min / dy
      end if

      ! for w(i,j-1,k+1):
      emop = ( dzf(k) * ( ekm(i,j,k+1)  + ekm(i,j-1,k+1)  )  + &
              dzf(k+1)  * ( ekm(i,j,k) + ekm(i,j-1,k) ) ) / &
              ( 4.0_field_r * dzh(k+1) )
      u_at_w_min = 0.25_field_r * ( u0(i,j-1,k) + u0(i,j-1,k+1) + u0(i+1,j-1,k) + u0(i+1,j-1,k+1) )
      tau_wv_min = log_wallaw(w0(i,j-1,k+1), u_at_w_min , Cm_ywall)

      wp(i,j-1,k+1) = wp(i,j-1,k+1) - 0.5_field_r * emop * ( (w0(i,j,k+1) - w0(i,j-1,k+1) ) / dy ) / dy - 0.5_field_r * tau_wv_min / dy

      !> Enforcing zero flux by correction of tendencies of temperature, moisture and other scalars (by negating diffusion term)

      ! scalars at (i,j,k) are "to the right" of the y-walls, so correct on s(i,j-1,k) which is on the left for minus y-walls
      ! "-" because we reflect it back, and 0.5_field_r comes from interpolation of ekh
      thlp(i,j-1,k) = thlp(i,j-1,k) - 0.5_field_r * ( ekh(i,j,k) + ekh(i,j-1,k) ) * ( thl0(i,j,k) - thl0(i,j-1,k) ) * dy2i
      qtp(i,j-1,k)  =  qtp(i,j-1,k) - 0.5_field_r * ( ekh(i,j,k) + ekh(i,j-1,k) ) * (  qt0(i,j,k) -  qt0(i,j-1,k) ) * dy2i

      do nc=1,nsv
        svp(i,j-1,k,nc) = svp(i,j-1,k,nc) - 0.5_field_r * ( ekh(i,j,k) + ekh(i,j-1,k) ) * ( sv0(i,j,k,nc) - sv0(i,j-1,k,nc) ) * dy2i
      end do

      !> Finally, set correct heat flux from wall to fluid
      uspeed = 0.5_field_r * ( ( u0(i,j-1,k) + u0(i+1,j-1,k) )**2 + ( w0(i,j-1,k) + w0(i,j-1,k+1) )**2 )**0.5_field_r
      thlp(i,j-1,k) = thlp(i,j-1,k) + Cd_ywall * uspeed * (thlwall - thl0(i,j-1,k)) / dy

    end do

    call timer_toc('modibm/applyibm')
    return
  end subroutine applyibm

  subroutine zerowallvelocity !<- MK: Force velocity at the immersed boundaries to 0 for a better interaction with the poissonsolver

    use modfields,      only : um, vm, wm, up, vp, wp, u0, v0, w0
    use modglobal,      only : rk3step, kmax, i1, j1, k1, ih, jh, rdt
    use modmpi,         only : excjs

    implicit none
    integer  :: i, j, k, nn
    real     :: rk3coef,rk3coefi

    rk3coef = rdt / (4. - dble(rk3step))
    rk3coefi = 1. / rk3coef

    !> Set tendencies inside obstacled (i.e., correct for any drift from previous integration step)
    do nn = 1,Nobst
      i = iobst(nn,1)
      j = iobst(nn,2)
      k = iobst(nn,3)

      ! Correction of velocities also corrects walls w. normals in negative x,y-directions (due to staggered grid arrangement)
      u0(i,j,k)   = 0!-um(i,j,k) * rk3coefi
      v0(i,j,k)   = 0!-vm(i,j,k) * rk3coefi
      ! w0(i,j,k)   = 0!-wm(i,j,k) * rk3coefi

      ! Do trick: moving one index up by 1 in each direction also corrects walls with normals pointing in positive direction
      u0(i+1,j,k) = 0!-um(i+1,j,k) * rk3coefi
      v0(i,j+1,k) = 0!-vm(i,j+1,k) * rk3coefi
      w0(i,j,k+1) = 0!-wm(i,j,k+1) * rk3coefi

    end do

    ! call excjs( up  , 2,i1,2,j1,1,k1,ih,jh)
    ! call excjs( vp  , 2,i1,2,j1,1,k1,ih,jh)
    ! call excjs( wp  , 2,i1,2,j1,1,k1,ih,jh)
    call excjs( u0  , 2,i1,2,j1,1,k1,ih,jh)
    call excjs( v0  , 2,i1,2,j1,1,k1,ih,jh)
    call excjs( w0  , 2,i1,2,j1,1,k1,ih,jh)
    return
  end subroutine zerowallvelocity

  !> Calculates the Obukhov length iteratively (modified from modsurface.f90 implementation)
  function getobl_local(uspeed,thl,qt,thlroof,qtroof,z_MO,z0m_wall,z0h_wall) result (Lob)
    !$acc routine seq
    real(field_r), intent(in) :: uspeed, thl, qt, thlroof, qtroof, z_MO, z0m_wall, z0h_wall

    real                :: Lob

    ! local internal variables
    integer             :: i, j, iter
    real                :: Rib, fx, fxdif, Lold, Lstart, Lend, thv, thvsl, thvroof, horv2

    thv   = thl * (1. + (rv/rd - 1.) * qt)
    thvsl = thlroof * (1. + (rv/rd - 1.) * qtroof)
    horv2 = max(uspeed**2, 0.01)

    Rib = grav / thvs * z_MO * (thv - thvsl) / horv2 !! WAAR KOMT THVS vandaan!!!!!!!!

    if (Rib == 0) then
        ! Rib can be 0 if there is no surface flux
        ! L is capped at 1e6 below, so use the same cap here
        Lob = 1e6
        write(*,*) 'Obukhov length: Rib = 0 -> setting Lob=1e6'
    else
        iter = 0
        ! L = obl(i,j) ! previous value is best guess for new value, yet we don't have that saved currently.. consider saving later..

        if(Rib > 0) Lob = 0.01
        if(Rib < 0) Lob = -0.01

        do while (.true.)
          iter    = iter + 1
          Lold    = Lob

          fx     = Rib - z_MO / Lob * (log(z_MO / z0m_wall) - psih(z_MO / Lob) + psih(z0m_wall / Lob)) / &
                    (log(z_MO / z0m_wall) - psim(z_MO / Lob) + psim(z0m_wall / Lob))**2.
          Lstart = Lob - 0.001*Lob
          Lend   = Lob + 0.001*Lob

          fxdif  = ( (- z_MO / Lstart * (log(z_MO / z0m_wall) - psih(z_MO / Lstart) + psih(z0m_wall / Lstart)) /&
                  (log(z_MO / z0m_wall) - psim(z_MO / Lstart) + psim(z0m_wall / Lstart)) ** 2.) - (-z_MO / Lend * &
                  (log(z_MO / z0m_wall) - psih(z_MO / Lend) + psih(z0m_wall / Lend)) / (log(z_MO / z0m_wall) - psim(z_MO / Lend)&
                  + psim(z0m_wall / Lend)) ** 2.) ) / (Lstart - Lend)

          Lob = Lob - fx / fxdif
          if(Rib * Lob < 0. .or. abs(Lob) == 1e5) then
              if(Rib > 0) Lob = 0.01
              if(Rib < 0) Lob = -0.01
          end if
          if(abs((Lob - Lold)/Lob) < 1e-4) exit
          if(iter > 1000) then
            stop 'Obukhov length calculation does not converge in IBM!'
          end if
        end do

        if (abs(Lob)>1e6) Lob = sign(1.0e6,Lob)
    end if

    return

  end function getobl_local

  function log_wallaw(u1,u2,Cm_hor_wall) result(tau)
    !$acc routine seq
    real(field_r), intent(in) :: u1,u2,Cm_hor_wall

    real(field_r)             :: tau

    tau  = Cm_hor_wall * sqrt(u1**2 + u2**2) * u1   !not a minus sign here but in the subroutine above, where it ensures force the direction of the wind
                                                    !similar as michael who states "give tau the same sign as utan"
  end function log_wallaw

  ! SvdL: wall corection for e12 currenly neglected. To be tested in detail later.
  ! So, we now allow for diffusive flux of e12 over the walls, while e12 and ekm/ekh should be insignificant inside buildings. Just outside ekm is probably too large to accurately determine such flux at all -> so we should consider proper cancellation.

  ! subroutine xwalle12(i,j,k)

  !   use modglobal,      only : dx2i, dx, dy, dzh
  !   use modsubgriddata, only : ekm
  !   use modfields,      only : e12p, e120, u0, v0, w0

  !   implicit none

  !   integer, intent(in)    :: i,j,k

  !   if(.not.(k==1)) then
  !     if(.not. (e12p(i,j,k)==0)) then
  !       e12p(i,j,k)   = e12p(i,j,k)   - (ekm(i,j,k)+ekm(i-1,j,k))*(e120(i,j,k)-e120(i-1,j,k))*dx2i &
  !                                 + ekm(i,j,k)/(2*e120(i,j,k))* (&  !source terms
  !                                    -((w0(i,j,k+1)-w0(i-1,j,k+1))  / dx             + &
  !                                      (u0(i,j,k+1)-u0(i,j,k))      / dzh(k+1) )**2  + &
  !                                    +(2.*(w0(i,j,k+1))             / dx             + &
  !                                      (u0(i,j,k+1)-u0(i,j,k))      / dzh(k+1) )**2  + &

  !                                    -((w0(i,j,k)-w0(i-1,j,k))      / dx             + &
  !                                      (u0(i,j,k)-u0(i,j,k-1))      / dzh(k)   )**2  + &
  !                                    +(2.*(w0(i,j,k))               / dx             + &
  !                                      (u0(i,j,k)-u0(i,j,k-1))      / dzh(k)   )**2  + &

  !                                    -((u0(i,j+1,k)-u0(i,j,k))      / dy             + &
  !                                      (v0(i,j+1,k)-v0(i-1,j+1,k))  / dx       )**2  + &
  !                                    +((u0(i,j+1,k)-u0(i,j,k))      / dy             + &
  !                                      (2.*v0(i,j+1,k))             / dx       )**2  + &

  !                                    -((u0(i,j,k)-u0(i,j-1,k))      / dy             + &
  !                                      (v0(i,j,k)-v0(i-1,j,k))      / dx       )**2  + &
  !                                    +((u0(i,j,k)-u0(i,j-1,k))      / dy             + &
  !                                      (2.*v0(i,j,k))               / dx       )**2    &
  !                                   )
  !     elseif(.not. (e12p(i-1,j,k)==0)) then
  !       e12p(i-1,j,k) = e12p(i-1,j,k) + (ekm(i,j,k)+ekm(i-1,j,k))*(e120(i,j,k)-e120(i-1,j,k))*dx2i &
  !                                 + ekm(i-1,j,k)/(2*e120(i-1,j,k))* (&  !source terms
  !                                      -((w0(i,j,k)-w0(i-1,j,k))    / dx             + &
  !                                        (u0(i,j,k)-u0(i,j,k-1))    / dzh(k)   )**2  + &
  !                                      +((-2.*w0(i-1,j,k))          / dx             + &
  !                                        (u0(i,j,k)-u0(i,j,k-1))    / dzh(k)   )**2  + &

  !                                      -((w0(i,j,k+1)-w0(i-1,j,k+1))/ dx             + &
  !                                        (u0(i,j,k+1)-u0(i,j,k))    / dzh(k+1) )**2  + &
  !                                      +((-2.*w0(i-1,j,k+1))        / dx             + &
  !                                        (u0(i,j,k+1)-u0(i,j,k))    / dzh(k+1) )**2  + &

  !                                      -((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
  !                                        (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &
  !                                      +((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
  !                                        (-2.*v0(i-1,j,k))          / dx       )**2  + &

  !                                      -((u0(i,j+1,k)-u0(i,j,k))    / dy             + &
  !                                        (v0(i,j+1,k)-v0(i-1,j+1,k))/ dx       )**2  + &
  !                                      +((u0(i,j+1,k)-u0(i,j,k))    / dy             + &
  !                                        (-2.*v0(i-1,j+1,k))        / dx       )**2    &
  !                                   )
  !   end if
  !   else !Special treatment for the lowest full level: k=1
  !     if(.not. (e12p(i,j,k)==0)) then
  !       e12p(i,j,k)   = e12p(i,j,k)   - (ekm(i,j,k)+ekm(i-1,j,k))*(e120(i,j,k)-e120(i-1,j,k))*dx2i &
  !                                 + ekm(i,j,k)/(2*e120(i,j,k))* (&  !source terms
  !                                    -((u0(i,j+1,k)-u0(i,j,k))      / dy             + &
  !                                      (v0(i,j+1,k)-v0(i-1,j+1,k))  / dx       )**2  + &
  !                                    +((u0(i,j+1,k)-u0(i,j,k))      / dy             + &
  !                                      (2.*v0(i,j+1,k))             / dx       )**2  + &

  !                                    -((u0(i,j,k)-u0(i,j-1,k))      / dy             + &
  !                                      (v0(i,j,k)-v0(i-1,j,k))      / dx       )**2  + &
  !                                    +((u0(i,j,k)-u0(i,j-1,k))      / dy             + &
  !                                      (2.*v0(i,j,k))               / dx       )**2    &
  !                                   )
  !     elseif(.not. (e12p(i-1,j,k)==0)) then
  !       e12p(i-1,j,k) = e12p(i-1,j,k) + (ekm(i,j,k)+ekm(i-1,j,k))*(e120(i,j,k)-e120(i-1,j,k))*dx2i &
  !                                 + ekm(i-1,j,k)/(2*e120(i-1,j,k))* (&  !source terms

  !                                      -((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
  !                                        (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &
  !                                      +((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
  !                                        (-2.*v0(i-1,j,k))          / dx       )**2  + &

  !                                      -((u0(i,j+1,k)-u0(i,j,k))    / dy             + &
  !                                        (v0(i,j+1,k)-v0(i-1,j+1,k))/ dx       )**2  + &
  !                                      +((u0(i,j+1,k)-u0(i,j,k))    / dy             + &
  !                                        (-2.*v0(i-1,j+1,k))        / dx       )**2    &
  !                                   )
  !     end if
  !   end if
  ! end subroutine xwalle12

  ! subroutine ywalle12(i,j,k)

  !   use modglobal,      only : dy2i, dx, dy, dzh
  !   use modsubgriddata, only : ekm
  !   use modfields,      only : e12p, e120, u0, v0, w0

  !   implicit none

  !   integer, intent(in)    :: i,j,k

  !   if(.not.(k==1)) then
  !     if(.not. (e12p(i,j,k)==0)) then
  !       e12p(i,j,k)   = e12p(i,j,k)   - (ekm(i,j,k)+ekm(i,j-1,k))*(e120(i,j,k)-e120(i,j-1,k))*dy2i &
  !                                 + ekm(i,j,k)/(2.*e120(i,j,k))* (&  !source terms
  !                                      -((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
  !                                        (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &
  !                                      +((2.*u0(i,j,k))             / dy             + &
  !                                        (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &

  !                                      -((u0(i+1,j,k)-u0(i+1,j-1,k))/ dy             + &
  !                                        (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2  + &
  !                                      +((2.*u0(i+1,j,k))           / dy             + &
  !                                        (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2  + &

  !                                      -((v0(i,j,k+1)-v0(i,j,k))    / dzh(k+1)       + &
  !                                        (w0(i,j,k+1)-w0(i,j-1,k+1))/ dy       )**2  + &
  !                                      +((v0(i,j,k+1)-v0(i,j,k))    / dzh(k+1)       + &
  !                                        (2.*w0(i,j,k+1))           / dy       )**2  + &

  !                                      -((v0(i,j,k)-v0(i,j,k-1))    / dzh(k)         + &
  !                                        (w0(i,j,k)-w0(i,j-1,k))    / dy       )**2  + &
  !                                      +((v0(i,j,k)-v0(i,j,k-1))    / dzh(k)         + &
  !                                        (2.*w0(i,j,k))             / dy       )**2    &
  !                                   )
  !     elseif(.not. (e12p(i,j-1,k)==0)) then
  !       e12p(i,j-1,k) = e12p(i,j-1,k) + (ekm(i,j,k)+ekm(i,j-1,k))*(e120(i,j,k)-e120(i,j-1,k))*dy2i &
  !                                 + ekm(i,j-1,k)/(2.*e120(i,j-1,k))* (&  !source terms
  !                                      -((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
  !                                        (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &
  !                                      +((-2.*u0(i,j-1,k))          / dy             + &
  !                                        (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &

  !                                      -((u0(i+1,j,k)-u0(i+1,j-1,k))/ dy             + &
  !                                        (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2  + &
  !                                      +((-2.*u0(i+1,j-1,k))        / dy             + &
  !                                        (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2  + &

  !                                      -((v0(i,j,k)-v0(i,j,k-1))    / dzh(k)         + &
  !                                        (w0(i,j,k)-w0(i,j-1,k))    / dy       )**2  + &
  !                                      +((v0(i,j,k)-v0(i,j,k-1))    / dzh(k)         + &
  !                                        (-2.*w0(i,j-1,k))          / dy       )**2  + &

  !                                      -((v0(i,j,k+1)-v0(i,j,k))    / dzh(k+1)       + &
  !                                        (w0(i,j,k+1)-w0(i,j-1,k+1))/ dy       )**2  + &
  !                                      +((v0(i,j,k+1)-v0(i,j,k))    / dzh(k+1)       + &
  !                                        (-2.*w0(i,j-1,k+1))        / dy       )**2    &
  !                                   )
  !     end if
  !   else !Special treatment for the lowest full level: k=1
  !     if(.not. (e12p(i,j,k)==0)) then
  !       e12p(i,j,k)   = e12p(i,j,k)   - (ekm(i,j,k)+ekm(i,j-1,k))*(e120(i,j,k)-e120(i,j-1,k))*dy2i &
  !                                 + ekm(i,j,k)/(2.*e120(i,j,k))* (&  !source terms
  !                                      -((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
  !                                        (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &
  !                                      +((2.*u0(i,j,k))             / dy             + &
  !                                        (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &

  !                                      -((u0(i+1,j,k)-u0(i+1,j-1,k))/ dy             + &
  !                                        (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2  + &
  !                                      +((2.*u0(i+1,j,k))           / dy             + &
  !                                        (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2    &
  !                                   )

  !     elseif(.not. (e12p(i,j-1,k)==0)) then
  !       e12p(i,j-1,k) = e12p(i,j-1,k) + (ekm(i,j,k)+ekm(i,j-1,k))*(e120(i,j,k)-e120(i,j-1,k))*dy2i &
  !                                 + ekm(i,j-1,k)/(2.*e120(i,j-1,k))* (&  !source terms
  !                                      -((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
  !                                        (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &
  !                                      +((-2.*u0(i,j-1,k))          / dy             + &
  !                                        (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &

  !                                      -((u0(i+1,j,k)-u0(i+1,j-1,k))/ dy             + &
  !                                        (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2  + &
  !                                      +((-2.*u0(i+1,j-1,k))        / dy             + &
  !                                        (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2    &
  !                                   )
  !     end if
  !   end if
  ! end subroutine ywalle12

  ! subroutine bulk_wall_temp(uspeed,thl,Cd,dx,thlp)
  !   implicit none

  !   real(field_r),intent(in) :: uspeed
  !   real,intent(in) :: Cd,dx
  !   real(field_r),intent(in) :: thl
  !   real(field_r),intent(out) :: thlp

  !   thlp = Cd * uspeed * (thlwall - thl) / dx

  !   return
  ! end subroutine bulk_wall_temp

end module modibm
