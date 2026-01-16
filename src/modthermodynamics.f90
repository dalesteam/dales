!>\file modthermodynamics.f90
!! Do the thermodynamics

!>
!! Do the thermodynamics
!>
!! Timeseries of the most relevant parameters. Written to tmser1.expnr and tmsurf.expnr
!! If netcdf is true, this module leads the tmser.expnr.nc output
!!  \author Pier Siebesma, K.N.M.I.
!!  \author Stephan de Roode,TU Delft
!!  \author Thijs Heus,MPI-M
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
!

module modthermodynamics
  use modglobal,    only: checknamelisterror, ifnamopt, i1, j1, k1, ih, jh, &
                          rv, rlv, cp, rd
  use modmpi,       only: myid, d_mpi_bcast, commwrld
  use modprecision, only : field_r
  use modtimer
  use modlogging, only: finish
  use fortran_support, only: nnml_output
  implicit none
  character(len=*), parameter :: modname = 'modthermodynamics'
!   private
  public :: thermodynamics,calc_halflev
  public :: ttab
  public :: esatltab
  public :: esatitab
  public :: esatmtab
  public :: calc_qsat
  public :: thermodynamics_read_namelist

  logical :: lmoist = .true.       !< Switch to calculate moisture fields.
  logical :: lnoclouds = .false.   !< Switch to enable/disable thl calculations.
  logical :: lconstexner = .false. !< Switch to use the initial pressure profile in the exner function.
  logical :: lbaseexner = .false.  !< Switch to use the base pressure profile in the exner function.

  real, allocatable :: th0av(:)
  real(field_r), allocatable :: thv0(:,:,:)
  real :: chi_half=0.5  !< set wet, dry or intermediate (default) mixing over the cloud edge
  real, allocatable :: thetah(:), qth(:), qlh(:)

  real(field_r), protected :: ttab(1:2000)
  real(field_r), protected :: esatltab(1:2000)
  real(field_r), protected :: esatitab(1:2000)
  real(field_r), protected :: esatmtab(1:2000)

  !$acc declare create(ttab, esatltab, esatitab, esatmtab)

contains

  !> Read thermodynamics namelist.
  subroutine thermodynamics_read_namelist(nml_filename)

    character(len=*), intent(in) :: nml_filename

    integer :: ierr

    namelist /thermodynamics/ lmoist, chi_half, lconstexner, lbaseexner, lnoclouds

    if (myid == 0) then
      open(ifnamopt, file=nml_filename, status='old', action='read', &
           iostat=ierr)
      read(ifnamopt, nml=thermodynamics, iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'thermodynamics')
      write(nnml_output, thermodynamics)
      close(ifnamopt)
    end if

    call d_mpi_bcast(lmoist, 1, 0, commwrld, ierr)
    call d_mpi_bcast(chi_half, 1, 0, commwrld, ierr)
    call d_mpi_bcast(lconstexner, 1, 0, commwrld, ierr)
    call d_mpi_bcast(lbaseexner, 1, 0, commwrld, ierr)

  end subroutine thermodynamics_read_namelist

!> Allocate and initialize arrays
  subroutine initthermodynamics
    use modglobal, only : ih,i1,jh,j1,k1,tdn,tup
    use modmicrodata, only: imicro,imicro_bulk3
    implicit none
    real :: ilratio
    integer :: m

    allocate(th0av(k1))
    allocate(thv0(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(thetah(k1), qth(k1), qlh(k1))

    th0av(:) = 0.

    !$acc enter data copyin(th0av, thv0, thetah, qth, qlh)

    ! esatltab(m) gives the saturation vapor pressure over water at T corresponding to m
    ! esatitab(m) is the same over ice
    ! esatmtab(m) is interpolated between the ice and liquid values with ilratio
    ! http://www.radiativetransfer.org/misc/atmlabdoc/atmlab/h2o/thermodynamics/e_eq_water_mk.html
    ! Murphy and Koop 2005 parameterization formula.
    do m=1,2000
       ttab(m)=150.+0.2*m
       esatltab(m)=exp(54.842763-6763.22/ttab(m)-4.21*log(ttab(m))+0.000367*ttab(m)+&
            tanh(0.0415*(ttab(m)-218.8))*(53.878-1331.22/ttab(m)-9.44523*log(ttab(m))+ 0.014025*ttab(m)))

       esatitab(m)=exp(9.550426-5723.265/ttab(m)+3.53068*log(ttab(m))-0.00728332*ttab(m))
       ilratio = max(0.,min(1.,(ttab(m)-tdn)/(tup-tdn)))
       if(imicro.eq.imicro_bulk3) then
          ! bulkmicro3 thermodynamics is for liquid only, ice is explicitely accounted for separately.
          esatmtab(m) = esatltab(m)
       else
          ! for all other microphysics, saturation is w.r.t. liquid and ice
          esatmtab(m) = ilratio*esatltab(m) + (1-ilratio)*esatitab(m)
       end if
    end do

    !$acc update device(ttab, esatltab, esatitab, esatmtab)

  end subroutine initthermodynamics

!> Do moist thermodynamics.
!! Calculate the liquid water content, do the microphysics, calculate the mean hydrostatic pressure,
!! calculate the fields at the half levels, and finally calculate the virtual potential temperature.
  subroutine thermodynamics
    use modglobal,  only : timee,k1,i1,j1,ih,jh,rd,rv,ijtot,cp,rlv
    use modfields,  only : thl0, qt0, ql0, presf, exnf, thvh, thv0h, qt0av, ql0av, thvf, rhof, ql0h, thl0h, qt0h, presh, exnh
    use modmpi,     only : slabsum
    use modibm,     only : fluid_mask
    use modibmdata, only : lapply_ibm
    use modslabaverage, only : slabavg
    implicit none
    integer:: i, j, k

    real(field_r) :: T
    logical :: too_hot, too_cold

    call timer_tic('modthermodynamics/thermodynamics', 0)

    if (timee < 0.01) then
      call diagfld
    end if
    if (lmoist .and. (.not. lnoclouds)) then
      
      ! Before we do the saturation adjustment, check if 150 K < T < 550 K
      ! If this is not the case, we will read outside of the bounds of
      ! esatmtab

      !$acc parallel loop collapse(3) default(present) async(1) private(T)
      do k = 1, k1
        do j = 2 , j1
          do i = 2, i1
            T = thl0(i,j,k) * exnf(k)
            if (T < 150) then
              !$acc atomic write
              too_cold = .true.
            else if (T > 550) then
              !$acc atomic write
              too_hot = .true.
            end if
          end do
        end do
      end do

      call saturation_adjustment(qt0, thl0, presf, exnf, ql0)
      call calc_dry_tmp ! tmp0 is used in statistics
                         ! can consider calculating it only when needed
    end if
    call diagfld

    call calc_halflev !calculate halflevel values of qt0 and thl0

    if (lmoist .and. (.not. lnoclouds)) then
      call saturation_adjustment(qt0h, thl0h, presh, exnh, ql0h)
    end if

    ! recalculate thv and rho on the basis of results
    call calthv

    !$acc parallel loop collapse(3) default(present) async(2)
    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          thv0(i,j,k) = (thl0(i,j,k)+rlv*ql0(i,j,k)/(cp*exnf(k))) &
                      * (1+(rv/rd-1)*qt0(i,j,k)-rv/rd*ql0(i,j,k))
        end do
      end do
    end do

    !$acc kernels default(present)
    thvh(:) = 0.0
    thvf(:) = 0.0
    !$acc end kernels

    if (.not. lapply_ibm) then
      call slabsum(thvh,1,k1,thv0h,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1, on_gpu=.true.) ! redefine halflevel thv using calculated thv
      call slabsum(thvf,1,k1,thv0,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1, on_gpu=.true.)
    else
      call slabavg(thv0h,fluid_mask,ih,thvh)
      call slabavg(thv0,fluid_mask,ih,thvf)
    end if

    !$acc kernels default(present) async(1)
    if (.not. lapply_ibm) then
      thvh(:) = thvh(:)/ijtot
      thvf(:) = thvf(:)/ijtot
    end if
    thvh(1) = th0av(1)*(1+(rv/rd-1)*qt0av(1)-rv/rd*ql0av(1)) ! override first level
    !$acc end kernels

    !$acc parallel loop default(present) async(1)
    do k = 1, k1
      rhof(k) = presf(k)/(rd*thvf(k)*exnf(k))
    end do

    !$acc wait
    call timer_toc('modthermodynamics/thermodynamics')
  end subroutine thermodynamics

!> Cleans up after the run
  subroutine exitthermodynamics
    implicit none
    !$acc exit data delete(th0av, thv0, thetah, qth, qlh)
    deallocate(th0av, thv0, thetah, qth, qlh)
  end subroutine exitthermodynamics

  !> Calculate real temperature tmp0 from thl0, for the dry case i.e. ql=0
  subroutine calc_dry_tmp
    use modglobal, only : i1,j1,k1
    use modfields, only : thl0,exnf
    use modfields, only : tmp0

    implicit none
    integer :: i, j, k

    !$acc parallel loop collapse(3) default(present) async
    do k = 1,k1
       do j = 2,j1
          do i = 2,i1
             tmp0(i,j,k) = exnf(k)*thl0(i,j,k)
          end do
       end do
    end do

  end subroutine calc_dry_tmp

!> Calculate thetav and dthvdz
  subroutine calthv
    use modglobal, only : i1,j1,k1,kmax,zf,dzh,rlv,rd,rv,cp,eps1
    use modfields, only : thl0,thl0h,ql0,ql0h,qt0,qt0h,exnf,exnh,thv0h,dthvdz
    use modsurfdata,only : dthldz,dqtdz
    implicit none

    integer i, j, k
    real(field_r)    qs
    real(field_r)    a_surf,b_surf,dq,dth,dthv,temp
    real(field_r)    a_dry, b_dry, a_moist, b_moist, c_liquid, epsilon, eps_I, chi_sat, chi
    real(field_r)    del_thv_sat, del_thv_dry

    call timer_tic('modthermodynamics/calthv', 1)

    dthvdz = 0

    if (lmoist) then
      !$acc parallel loop collapse(3) default(present) async(1)
      do k = 2, k1
        do j = 2, j1
          do i = 2, i1
            thv0h(i,j,k) = (thl0h(i,j,k)+rlv*ql0h(i,j,k)/(cp*exnh(k))) &
                          *(1+(rv/rd-1)*qt0h(i,j,k)-rv/rd*ql0h(i,j,k))
          end do
        end do
      end do

      !TODO: fix the branching in this loop
      !$acc parallel loop collapse(3) default(present) &
      !$acc& private(a_dry, b_dry, a_moist, b_moist, c_liquid, epsilon, eps_I, chi_sat, chi, dthv, del_thv_dry, del_thv_sat, temp, qs, dq, dth) &
      !$acc& async(2)
      do k = 2, kmax
        do j = 2 , j1
          do i = 2, i1
!
!         default thv jump computed unsaturated
!
            epsilon = rd/rv
            eps_I = 1/epsilon - 1  !cstep approx 0.608

            a_dry = 1 + eps_I * qt0(i,j,k)
            b_dry = eps_I * thl0(i,j,k)

            dth = thl0(i,j,k+1)-thl0(i,j,k-1)
            dq  = qt0(i,j,k+1)-qt0(i,j,k-1)

            del_thv_dry = a_dry   * dth + b_dry * dq

            dthv = del_thv_dry

            if  (ql0(i,j,k)> 0) then  !include moist thermodynamics

               temp = thl0(i,j,k)*exnf(k)+(rlv/cp)*ql0(i,j,k)
               qs   = qt0(i,j,k) - ql0(i,j,k)

               a_moist = (1-qt0(i,j,k)+qs/epsilon*(1+rlv/(rv*temp))) &
                        /(1+rlv**2*qs/(cp*rv*temp**2))
               b_moist = a_moist*rlv/cp-temp
               c_liquid = a_dry * rlv / cp - thl0(i,j,k) / epsilon

               del_thv_sat = a_moist * dth + b_moist * dq

               chi     = 2*chi_half*(zf(k) - zf(k-1))/(dzh(k)+dzh(k+1))
               chi_sat = c_liquid * ql0(i,j,k) / (del_thv_dry - del_thv_sat)

               if (chi < chi_sat) then  !mixed parcel is saturated
                 dthv = del_thv_sat
              end if
            end if

            dthvdz(i,j,k) = dthv/(dzh(k+1)+dzh(k))
          end do
        end do
      end do

      !$acc parallel loop collapse(2) default(present) private(temp, qs, a_surf, b_surf) async(3)
      do j=2,j1
        do i=2,i1
          if(ql0(i,j,1)>0) then
            temp = thl0(i,j,1)*exnf(1)+(rlv/cp)*ql0(i,j,1)
            qs   = qt0(i,j,1) - ql0(i,j,1)
            a_surf   = (1-qt0(i,j,1)+rv/rd*qs*(1+rlv/(rv*temp))) &
                      /(1+rlv**2*qs/(cp*rv*temp**2))
            b_surf   = a_surf*rlv/(temp*cp)-1

          else
            a_surf = 1+(rv/rd-1)*qt0(i,j,1)
            b_surf = rv/rd-1

          end if
          dthvdz(i,j,1) = a_surf*dthldz(i,j) + b_surf*thl0(i,j,1)*dqtdz(i,j)
        end do
      end do

    else
      !$acc parallel loop collapse(3) default(present)
      do k = 2, k1
        do j = 2, j1
          do i = 2, i1
            thv0h(i,j,k)  = thl0h(i,j,k)
          end do
        end do
      end do

      !$acc parallel loop collapse(3) default(present)
      do k = 2, kmax
        do j = 2, j1
          do i = 2, i1
            dthvdz(i,j,k) = (thl0(i,j,k+1)-thl0(i,j,k-1))/(dzh(k+1)+dzh(k))
          end do
        end do
      end do

      !$acc parallel loop collapse(2) default(present)
      do j = 2, j1
        do i = 2, i1
          dthvdz(i,j,1) = dthldz(i,j)
        end do
      end do
    end if

    !$acc parallel loop collapse(3) default(present) async wait(2, 3)
    do k = 1, kmax
      do j = 2, j1
        do i = 2, i1
          if(abs(dthvdz(i,j,k)) < eps1) then
            dthvdz(i,j,k) = sign(eps1, dthvdz(i,j,k))
          end if
        end do
      end do
    end do

    !$acc wait

    call timer_toc('modthermodynamics/calthv')

  end subroutine calthv
!> Calculate diagnostic slab averaged fields.
!!     Calculates slab averaged fields assuming
!!     hydrostatic equilibrium for: u,v,theta_l,theta_v,
!!     qt,ql,exner,pressure and the density
!! \author      Pier Siebesma   K.N.M.I.     06/01/1995
  subroutine diagfld
  use modglobal,  only : i1,ih,j1,jh,k1,nsv,zh,zf,cu,cv,ijtot,grav,rlv,cp,rd,rv,pref0,timee
  use modfields,  only : u0,v0,thl0,qt0,ql0,sv0,u0av,v0av,thl0av,qt0av,ql0av,sv0av, &
                        presf,presh,exnf,exnh,rhof,thvf
  use modsurfdata,only : thls,ps
  use modmpi,     only : slabsum
  use modibm,     only : fluid_mask
  use modibmdata, only : lapply_ibm
  use modslabaverage, only : slabavg
  implicit none

  integer :: k,n

  call timer_tic('modthermodynamics/diagfld', 1)


!*********************************************************
!  1.0   calculate average profiles of u,v,thl,qt and ql *
!        assuming hydrostatic equilibrium                *
!*********************************************************

! initialise local MPI arrays

  !$acc kernels default(present)
  u0av = 0.0
  v0av = 0.0
  thl0av = 0.0
  th0av  = 0.0
  qt0av  = 0.0
  ql0av  = 0.0
  sv0av = 0.
  !$acc end kernels

  !CvH changed momentum array dimensions to same value as scalars!
  if (.not. lapply_ibm) then
    call slabsum(u0av  ,1,k1,u0  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1, on_gpu=.true.)
    call slabsum(v0av  ,1,k1,v0  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1, on_gpu=.true.)
    call slabsum(thl0av,1,k1,thl0,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1, on_gpu=.true.)
    call slabsum(qt0av ,1,k1,qt0 ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1, on_gpu=.true.)
    call slabsum(ql0av ,1,k1,ql0 ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1, on_gpu=.true.)
    do n=1,nsv
      call slabsum(sv0av(1:1,n),1,k1,sv0(:,:,:,n),2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1, on_gpu=.true.)
    end do
  else
    call slabavg(u0,fluid_mask,ih,u0av)
    call slabavg(v0,fluid_mask,ih,v0av)
    call slabavg(thl0,fluid_mask,ih,thl0av)
    call slabavg(qt0,fluid_mask,ih,qt0av)
    call slabavg(ql0,fluid_mask,ih,ql0av)
    do n=1,nsv
      call slabavg(sv0(:,:,:,n),fluid_mask,ih,sv0av(:,n))
    end do
  end if

  if (.not. lapply_ibm) then
    !$acc kernels default(present)
    u0av   = u0av  /ijtot + cu
    v0av   = v0av  /ijtot + cv
    thl0av = thl0av/ijtot
    qt0av  = qt0av /ijtot
    ql0av  = ql0av /ijtot
    sv0av  = sv0av /ijtot
    !$acc end kernels
  end if
  if ((timee < 0.01 .or. .not. lconstexner) .and. .not. lbaseexner) then
    !$acc kernels default(present)
    exnf   = 1-grav*zf/(cp*thls)
    exnh   = 1-grav*zh/(cp*thls)
    !$acc end kernels
  endif
  !$acc kernels default(present)
  th0av  = thl0av+ (rlv/cp)*ql0av/exnf
  !$acc end kernels

!***********************************************************
!  2.0   calculate average profile of pressure at full and *
!        half levels, assuming hydrostatic equilibrium.    *
!***********************************************************

!    2.1 Use first guess of theta, then recalculate theta

   call fromztop

   !$acc kernels default(present)
   th0av = thl0av + (rlv/cp)*ql0av/exnf
   if ((timee < 0.01 .or. .not. lconstexner) .and. .not. lbaseexner) then
      exnf = (presf/pref0)**(rd/cp)
   endif
   !$acc end kernels

!    2.2 Use new updated value of theta for determination of pressure

   call fromztop

!***********************************************************
!  3.0   Construct density profiles and exner function     *
!       for further use in the program                     *
!***********************************************************

!  3.1 determine exner
   if ((timee < 0.01 .or. .not. lconstexner) .and. .not. lbaseexner) then
     !$acc serial default(present) async(1)
     exnh(1) = (ps/pref0)**(rd/cp)
     exnf(1) = (presf(1)/pref0)**(rd/cp)
     !$acc end serial

     !$acc parallel loop default(present) async(2)
     do k=2,k1
       exnf(k) = (presf(k)/pref0)**(rd/cp)
       exnh(k) = (presh(k)/pref0)**(rd/cp)
     end do
   endif

!  3.2 determine rho
   !$acc parallel loop default(present) async wait(1, 2)
   do k=1,k1
     thvf(k) = th0av(k)*exnf(k)*(1+(rv/rd-1)*qt0av(k)-rv/rd*ql0av(k))
     rhof(k) = presf(k)/(rd*thvf(k))
   end do
   !$acc wait

   call timer_toc('modthermodynamics/diagfld')

   return
  end subroutine diagfld

!> Calculates slab averaged pressure
!!      Input :  zf,zh,theta and qt profile
!!      Output:  pressure profile at full and
!!               half levels
!!
!!      Method: Using hydrostatic equilibrium
!!
!!                              -g*pref0**(rd/cp)
!! =====>       dp**(rd/cp)/dz = --------------
!!                                 cp*thetav
!! \author Pier Siebesma   K.N.M.I.     06/01/1995
  subroutine fromztop

  use modglobal, only : k1,dzf,dzh,rv,rd,cp,zf,grav,pref0
  use modfields, only : qt0av,ql0av,presf,presh,thvh,thvf
  use modsurfdata,only : ps
  implicit none

  integer   k
  real(field_r)  rdocp

  call timer_tic('modthermodynamics/fromztop', 1)

  rdocp = rd/cp

!**************************************************
!    1.0 Determine theta and qt at half levels    *
!**************************************************

  !$acc parallel loop default(present)
  do k=2,k1
    thetah(k) = (th0av(k)*dzf(k-1) + th0av(k-1)*dzf(k))/(2*dzh(k))
    qth   (k) = (qt0av(k)*dzf(k-1) + qt0av(k-1)*dzf(k))/(2*dzh(k))
    qlh   (k) = (ql0av(k)*dzf(k-1) + ql0av(k-1)*dzf(k))/(2*dzh(k))
  end do

!**************************************************
!     2.1  calculate pressures at full levels     *
!          assuming hydrostatic equilibrium       *
!**************************************************

!     1: lowest level: use first level value for safety!

  !$acc update self(thetah, qth, qlh, th0av, qt0av, ql0av)

  thvh(1) = th0av(1)*(1+(rv/rd-1)*qt0av(1)-rv/rd*ql0av(1))
  presf(1) = ps**rdocp - grav*(pref0**rdocp)*zf(1) /(cp*thvh(1))
  presf(1) = presf(1)**(1/rdocp)

!     2: higher levels

  do k=2,k1
    thvh(k)  = thetah(k)*(1+(rv/rd-1)*qth(k)-rv/rd*qlh(k))
    presf(k) = presf(k-1)**rdocp - &
                   grav*(pref0**rdocp)*dzh(k) /(cp*thvh(k))
    presf(k) = presf(k)**(1/rdocp)
  end do

!**************************************************
!     2.2   calculate pressures at half levels    *
!           assuming hydrostatic equilibrium      *
!**************************************************

  presh(1) = ps
  thvf(1) = th0av(1)*(1+(rv/rd-1)*qt0av(1)-rv/rd*ql0av(1))

  do k=2,k1
    thvf(k)  = th0av(k)*(1+(rv/rd-1)*qt0av(k)-rv/rd*ql0av(k))
    presh(k) = presh(k-1)**rdocp - &
                   grav*(pref0**rdocp)*dzf(k-1) / (cp*thvf(k-1))
    presh(k) = presh(k)**(1/rdocp)
  end do

  !$acc update device(thvh, presf, thvf, presh)
  call timer_toc('modthermodynamics/fromztop')

  return
  end subroutine fromztop

!> Magnus formulas for q_sat over liquid and ice
!> from Huang 2018 https://doi.org/10.1175/JAMC-D-17-0334.
!> Warning: for performance, check that rd/rv etc are pre-computed
  pure function qsat_magnus(T, p) result(qsat)
    use modglobal, only : rd,rv,tup,tdn
    implicit none
    real(field_r), intent(in) :: T, p
    real :: qsat
    real ilratio, TC, esl, esi, es
    ilratio = max(0._field_r,min(1._field_r,(T-tdn)/(tup-tdn)))

    TC = T - 273.15 ! in Celcius
    esl = 610.94_field_r * exp( (17.625_field_r*TC) / (TC+243.04_field_r) ) ! Magnus
    esi = 611.21_field_r * exp( (22.587_field_r*TC) / (TC+273.86_field_r) ) ! Magnus

    ! interpolated saturation vapor pressure
    es = ilratio*esl + (1-ilratio)*esi

    ! convert saturation vapor pressure to saturation humidity
    qsat = (rd/rv) * es / (p - (1-rd/rv)*es)
  end function qsat_magnus

!> Huang's formulas for q_sat over liquid and ice
!> from Huang 2018 https://doi.org/10.1175/JAMC-D-17-0334.
!> should be more accurate than Magnus, at the cost of more divisions
!> Warning: for performance, check that rd/rv etc are pre-computed
  pure function qsat_huang(T, p) result(qsat)
    use modglobal, only : rd,rv,tup,tdn
    implicit none
    real(field_r), intent(in) :: T, p
    real :: qsat
    real ilratio, TC, esl, esi, es
    ilratio = max(0._field_r,min(1._field_r,(T-tdn)/(tup-tdn)))

    TC = T - 273.15_field_r ! in Celcius
    esl = exp(34.494_field_r - 4924.99_field_r / (TC  + 237.1_field_r)) /  (TC+105)**1.57_field_r  ! Huang
    esi = exp(43.494_field_r - 6545.8_field_r/(TC+278)) / (TC+868)**2              ! Huang

    ! interpolated saturation vapor pressure
    es = ilratio*esl + (1-ilratio)*esi

    ! convert saturation vapor pressure to saturation humidity
    qsat = (rd/rv) * es / (p - (1-rd/rv)*es)
  end function qsat_huang

  ! return esat for ice-liquid mix using table
  pure function esat_tab(T) result(es)

    implicit none
    !$acc routine seq
    real(field_r), intent(in) :: T
    integer :: tlonr
    real(field_r) :: tlo, thi, es

    ! interpolated ice-liquid saturation vapor pressure from table
    ! note if imicto==imicro_bulk3, the table is for liquid only
    tlonr=int((T-150)*5)
    tlo = 150 + 0.2_field_r*tlonr
    thi = tlo + 0.2_field_r
    es = (thi-T)*5*esatmtab(tlonr)+(T-tlo)*5*esatmtab(tlonr+1)
  end function esat_tab

!> q_sat over liquid and ice, using interpolation in a table created in modglobal.
!> seems to be faster than the Magnus formula (on CPU)
  pure function qsat_tab(T, p) result(qsat)
    use modglobal, only : rd,rv

    implicit none
    !$acc routine seq
    real(field_r), intent(in) :: T, p
    real(field_r) :: qsat
    integer :: tlonr
    real(field_r) :: tlo, thi, es

    ! interpolated ice-liquid saturation vapor pressure from table
    tlonr=int((T-150)*5)
    tlo = 150 + 0.2_field_r*tlonr
    thi = tlo + 0.2_field_r
    es = (thi-T)*5*esatmtab(tlonr)+(T-tlo)*5*esatmtab(tlonr+1)

    ! convert saturation vapor pressure to saturation humidity
    qsat = (rd/rv) * es / (p - (1-rd/rv)*es)
  end function qsat_tab

  !> Compute the saturation specific humidity
  !!
  !! This is just a wrapper around qsat_tab, but that can be changed to any of
  !! the other qsat functions if needed.
  pure function calc_qsat(T, p) result(qsat)

    real(field_r), intent(in) :: T !< Temperature [K]
    real(field_r), intent(in) :: p !< Pressure [Pa]

    real(field_r) :: qsat !< Saturation specific humidity [kg/kg]

    qsat = qsat_tab(T, p)

  end function calc_qsat

  !> Compute the cloud water content via the saturation adjustment method.
  subroutine saturation_adjustment(qt, thl, pres, exn, ql, opt_stream)

    real(field_r), intent(in) :: qt(2-ih:,2-jh:,:)  !< Total water specific humidity [kg/kg]
    real(field_r), intent(in) :: thl(2-ih:,2-jh:,:) !< Liquid water potential temperature [K]
    real(field_r), intent(in) :: pres(:)            !< Pressure [Pa]
    real(field_r), intent(in) :: exn(:)             !< Exner function [-]

    real(field_r), intent(inout) :: ql(2-ih:,2-jh:,:) !< Liquid water specific humidity [kg/kg].

    integer, optional, intent(in) :: opt_stream !< (Optional) OpenACC stream ID.

    character(len=*), parameter :: routine = modname//'saturation_adjustment'

    integer :: i, j, k
    integer :: stream = 1 !< OpenACC stream ID.

    real(field_r) :: b      !< Factor in the equation for saturation specific humidity [-]
    real(field_r) :: qli    !< Intermediate value of liquid water specific humidity [kg/kg]
    real(field_r) :: qsat   !< Saturation specific humidity [kg/kg]
    real(field_r) :: qti    !< Intermediate value of total water specific humditiy [kg/kg]
    real(field_r) :: Tl     !< Liquid water temperature [K]
    real(field_r) :: Tl_min !< Minimum value of the liquid water temperature [K]
    real(field_r) :: qt_max !< Maximum value of the total water specific humidity [kg/kg]

    if (present(opt_stream)) stream = opt_stream

    call timer_tic(routine, 1)

    !$acc parallel loop gang vector collapse(3) default(present) async(stream) &
    !$acc private(b, qli, qsat, qti, Tl)
    do k = 1, k1
      ! Find lowest thl and highest qt in the slab.
      ! If they in combination are not saturated, the whole slab is below saturation.
      !
      ! TODO: on GPU, test if it's cheaper to just do the computation instead.
      TL_min = minval(thl(2:i1,2:j1,k)) * exn(k)
      qt_max = maxval(qt(2:i1,2:j1,k))
      qsat = qsat_tab(TL_min, pres(k))
      if (qt_max > qsat) then
        do j = 2, j1
          do i = 2, i1
            qti = qt(i,j,k)

            ! First step
            Tl = exn(k) * thl(i,j,k)
            qsat = qsat_tab(Tl, pres(k))
            b = rlv**2 / (rv * cp * Tl**2)
            qsat = qsat * (1 + b * qti) / (1 + b * qsat)

            ! Update the starting point
            qli = max(qti - qsat, 0.0_field_r)
            Tl = Tl + (rlv / cp) * qli
            qti = qti - qli

            ! Second step
            qsat = qsat_tab(Tl, pres(k))
            b = rlv**2 / (rv * cp * Tl**2)
            qsat = qsat * (1 + b * qti) / (1 + b * qsat)

            ql(i,j,k) = max(qt(i,j,k) - qsat, 0.0_field_r)
          end do
        end do
      end if
    end do
    
    call timer_toc(routine)

  end subroutine saturation_adjustment

  !> Diagnose saturation specific humidities over liquid and ice.
  subroutine calc_saturation_humidities(qt, ql, thl, pres, exn, esl, qvsl, qvsi)

    real(field_r), intent(in) :: qt(2-ih:,2-jh:,:)  !< Total water specific humidity [kg/kg]
    real(field_r), intent(in) :: ql(2-ih:,2-jh:,:)  !< Liquid water specific humidity [kg/kg]
    real(field_r), intent(in) :: thl(2-ih:,2-jh:,:) !< Liquid water potential temperature [K]
    real(field_r), intent(in) :: pres(:)            !< Pressure [Pa]
    real(field_r), intent(in) :: exn(:)             !< Exner function [-]

    real(field_r), intent(out) :: esl(2-ih:,2-jh:,:)  !< Liquid water saturation pressure [Pa]
    real(field_r), intent(out) :: qvsl(2-ih:,2-jh:,:) !< Liquid water saturation humidity [kg/kg]
    real(field_r), intent(out) :: qvsi(2-ih:,2-jh:,:) !< Ice saturation humidity [kg/kg]

    character(len=*), parameter :: routine = &
      modname//'calc_saturation_humidities'

    integer :: i, j, k
    
    real(field_r) :: esi   !< Saturation vapor pressure for ice (not stored) [Pa]
    real(field_r) :: qsat  !< Saturation specific humidity [kg/kg]
    real(field_r) :: T     !< Temperature [K]
    real(field_r) :: thi   !< Upper bound temperature for interpolation [K]
    real(field_r) :: tlo   !< Lower bound temperature for interpolation [K]
    real(field_r) :: tlonr !< Index of temperature in esat lookuptable

    !$acc parallel loop collapse(3) default(present) async(1) &
    !$acc private(T, tlonr, tlo, thi, esi)
    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          qsat = max(qt(i,j,k) - ql(i,j,k), 1.0_field_r)
          T = exn(k) * thl(i,j,k) + (rlv / cp) * ql(i,j,k)
          tlonr = int((T - 150) * 5)
          tlo = 150 + 0.2_field_r * tlonr
          thi = tlo + 0.2_field_r

          ! Liquid
          esl(i,j,k) = (thi - T) * 5 * esatltab(tlonr) &
                       + (T - tlo) * 5 * esatltab(tlonr + 1)
          qvsl(i,j,k) = rd / rv * esl(i,j,k) &
                        / (pres(k) - (1 - rd / rv) * esl(i,j,k))

          ! Ice
          esi = (thi - T) * 5 * esatitab(tlonr) &
                + (T - tlo) * 5 * esatitab(tlonr + 1)
          qvsi(i,j,k) = rd / rv * esi / (pres(k) - (1 - rd / rv) * esi)
        end do
      end do
    end do

  end subroutine calc_saturation_humidities

!> Calculates the scalars at half levels.
!! If the kappa advection scheme is active, interpolation needs to be done consistently.
  subroutine calc_halflev
    use modglobal, only : i1, j1, k1, dzf, dzhi, iadv_thl, iadv_qt, iadv_kappa
    use modfields, only : thl0, thl0h, qt0, qt0h
    use modsurfdata,only: qts, thls
    use advec_kappa,only: halflev_kappa
    implicit none

    integer :: i, j, k

    call timer_tic('modthermodynamics/calc_halflev', 1)

    if (iadv_thl==iadv_kappa) then
      call halflev_kappa(thl0,thl0h)
    else
      !$acc parallel loop collapse(3) default(present) async(1)
      do k = 2, k1
        do j = 2 ,j1
          do i = 2 ,i1
            thl0h(i,j,k) = (thl0(i,j,k)*dzf(k-1)+thl0(i,j,k-1)*dzf(k)) * (0.5_field_r * dzhi(k))
          end do
        end do
      end do
    end if

    !$acc parallel loop collapse(2) default(present) async(1)
    do j = 2, j1
      do i = 2, i1
        thl0h(i,j,1) = thls
      end do
    end do

    if (iadv_qt==iadv_kappa) then
      call halflev_kappa(qt0,qt0h)
    else
      !$acc parallel loop collapse(3) default(present) async(1)
      do k = 2, k1
        do j = 2, j1
          do i = 2, i1
            qt0h(i,j,k)  = (qt0(i,j,k)*dzf(k-1)+qt0(i,j,k-1)*dzf(k)) * (0.5_field_r * dzhi(k))
          end do
        end do
      end do

      !$acc parallel loop collapse(2) default(present) async(1)
      do j = 2, j1
        do i = 2, i1
          qt0h(i,j,1) = qts
        end do
      end do
    end if

    !$acc wait(1)
    call timer_toc('modthermodynamics/calc_halflev')
  end subroutine calc_halflev

end module modthermodynamics
