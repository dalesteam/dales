!>\file modthermodynamics.f90
!! Do the thermodynamics
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
!  Copyright 1993-2026 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!

!> Module for thermodynamics calculations.
module modthermodynamics
  use modglobal,       only: checknamelisterror, ifnamopt, i1, j1, k1, ih, jh, &
                             rv, rlv, cp, rd, dzf, dzhi, iadv_kappa, iadv_qt, &
                             iadv_thl, tdn, tup, timee, ijtot, kmax, zf, dzh, &
                             eps1, cu, cv, grav, pref0, nsv, zh
  use modfields,       only: qt0, thl0, qt0h, thl0h, ql0, presf, exnf, thvh, &
                             thv0h, qt0av, ql0av, thvf, rhof, ql0h, presh, exnh, &
                             u0, v0, sv0, u0av, v0av, thl0av, ql0av, sv0av, &
                             tmp0, dthvdz, thl0h, qt0h, esl, qvsl, qvsi
  use modsurfdata,     only: qts, thls, ps, dthldz, dqtdz
  use modmpi,          only: myid, d_mpi_bcast, commwrld, slabsum
  use modmicrodata,    only: imicro, imicro_bulk3, imicro_none
  use modibmdata,      only: lapply_ibm, fluid_mask
  use modslabaverage,  only: slabavg
  use modslabaverage,  only: slabavg_gpu
  use advec_kappa,     only: halflev_kappa
  use modprecision,    only: field_r
  use modtimer,        only: timer_tic, timer_toc
  use fortran_support, only: nnml_output, finish
  use modgpu
  implicit none
  character(len=*), parameter :: modname = 'modthermodynamics'
!   private
  public :: thermodynamics
  public :: ttab
  public :: esatltab
  public :: esatitab
  public :: esatmtab
  public :: calc_qsat
  public :: thermodynamics_read_namelist
  public :: calc_halflev

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
  !$omp declare target (ttab,esatltab,esatitab,esatmtab)

contains

  !> Calculate the virtual potential temperature.
  elemental function calc_virt_pot_temp(thl, qt, ql, exn) result(thv)
    !$omp declare target

    real(field_r), intent(in) :: thl  !< Liquid water potential temperature [K]
    real(field_r), intent(in) :: qt   !< Total water specific humidity [kg/kg]
    real(field_r), intent(in) :: ql   !< Liquid water specific humidity [kg/kg]
    real(field_r), intent(in) :: exn  !< Exner function [-]

    real(field_r) :: thv !< Virtual potential temperature [K]

    thv = (thl + rlv * ql / (cp * exn)) * (1 + (rv / rd - 1) * qt - rv / rd * ql)

  end function calc_virt_pot_temp

  !> Read thermodynamics namelist.
  subroutine thermodynamics_read_namelist(nml_filename)

    character(len=*), intent(in) :: nml_filename

    integer :: ierr
    logical :: lqlnr = .true. !< deprecated and ignored, kept for compatibility

    namelist /thermodynamics/ lmoist, chi_half, lconstexner, lbaseexner, &
                              lnoclouds, lqlnr

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
    real :: ilratio
    integer :: m

    allocate(th0av(k1))
    allocate(thv0(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(thetah(k1), qth(k1), qlh(k1))

    th0av(:) = 0.

    !$acc enter data copyin(th0av, thv0, thetah, qth, qlh)
    !$omp target enter data map(to:th0av,thv0,thetah,qth,qlh)

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
    !$omp target update to(ttab,esatltab,esatitab,esatmtab)

  end subroutine initthermodynamics

  !> Main thermodynamics subroutine.
  !!
  !! Calculates liquid water content, mean hydrostatic pressure and virtual 
  !! potential temperature.
  subroutine thermodynamics
    character(len=*), parameter :: routine = modname//'/thermodynamics'

    integer:: i, j, k

    real(field_r) :: T
    logical :: too_hot, too_cold

    call timer_tic(routine, 0)

    !$omp target update to(u0, v0, thl0, qt0, ql0)
    !$omp target update to(presf, presh, exnf, exnh, ql0h)
    !$omp target update to(zf, zh, dzf, dzhi, dthldz)

    if (timee < 0.01) then
      call diagfld
    end if

    if (lmoist .and. (.not. lnoclouds)) then

      ! Before we do the saturation adjustment, check if 150 K < T < 550 K
      ! If this is not the case, we will read outside of the bounds of
      ! esatmtab

      too_cold = .false.
      too_hot = .false.

      !$acc parallel loop collapse(3) default(present) async(1) private(T) &
      !$acc firstprivate(too_cold, too_hot)
      !$omp target teams loop private(T) collapse(3) reduction(.or.: too_cold,&
      !$omp too_hot) defaultmap(present:aggregate)&
      !$omp defaultmap(present:allocatable)
      do k = 1, k1
        do j = 2, j1
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

      !$acc wait
      if (too_cold) then
        call finish(routine, 'temperature below 150 K encountered!')
      else if (too_hot) then
        call finish(routine, 'temperature above 550 K encountered!')
      end if

      ! Do the saturation adjustment on the full levels
#if defined(DALES_GPU)
      call saturation_adjustment_gpu(qt0, thl0, presf, exnf, ql0, opt_stream=1)
#else
      call saturation_adjustment(qt0, thl0, presf, exnf, ql0, opt_stream=1)
#endif

      call diagfld

      ! Interpolate thl and qt to the half levels
      call calc_halflev(thl0, dzf, dzhi, thls, iadv_thl == iadv_kappa, thl0h)
      call calc_halflev(qt0, dzf, dzhi, qts, iadv_qt == iadv_kappa, qt0h)

      ! Do saturation adjustment again on the half levels
#if defined(DALES_GPU)
      call saturation_adjustment_gpu(qt0h, thl0h, presh, exnh, ql0h, opt_stream=1)
#else
      call saturation_adjustment(qt0h, thl0h, presh, exnh, ql0h, opt_stream=1)
#endif

      if (imicro /= imicro_none) then
        call calc_saturation_humidities(qt0, ql0, thl0, presf, exnf, esl, &
                                        qvsl, qvsi)
      end if
    else
      call calc_dry_tmp ! tmp0 is used in statistics
                         ! can consider calculating it only when needed
      call diagfld

      ! Interpolate thl and qt to the half levels
      call calc_halflev(thl0, dzf, dzhi, thls, iadv_thl == iadv_kappa, thl0h)
      call calc_halflev(qt0, dzf, dzhi, qts, iadv_qt == iadv_kappa, qt0h)

    end if

    ! recalculate thv and rho on the basis of results
    call calthv

    !$acc parallel loop collapse(3) default(present) async(1)
    !$omp target teams loop collapse(3) defaultmap(present:aggregate)&
    !$omp defaultmap(present:allocatable)
    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          thv0(i,j,k) = calc_virt_pot_temp(thl0(i,j,k), qt0(i,j,k), &
                                           ql0(i,j,k), exnf(k))
        end do
      end do
    end do

    !$acc parallel loop gang(static:1) default(present) async wait(1)
    !$omp target teams loop defaultmap(present:aggregate)&
    !$omp defaultmap(present:allocatable)
    do k = 1, k1
      thvh(k) = 0.0_field_r
    end do

    !$acc parallel loop gang(static:1) default(present) async wait(1)
    !$omp target teams loop defaultmap(present:aggregate)&
    !$omp defaultmap(present:allocatable)
    do k = 1, k1
      thvf(k) = 0.0_field_r
    end do

    !$acc wait

    if (.not. lapply_ibm) then
      call slabavg_gpu(thv0h, ih, thvh)
      call slabavg_gpu(thv0, ih, thvf)
    else
      call slabavg(thv0h,fluid_mask,ih,thvh)
      call slabavg(thv0,fluid_mask,ih,thvf)
    end if

    !$acc serial default(present) async(1)
    !$omp target defaultmap(present:allocatable)
    thvh(1) = th0av(1)*(1+(rv/rd-1)*qt0av(1)-rv/rd*ql0av(1))
    !$acc end serial
    !$omp end target

    !$acc parallel loop default(present) async(1)
    !$omp target teams loop defaultmap(present:aggregate)&
    !$omp defaultmap(present:allocatable)
    do k = 1, k1
      rhof(k) = presf(k)/(rd*thvf(k)*exnf(k))
    end do

    !$acc wait
    !$omp target update from(thv0h,tmp0,dthvdz)
    !$omp target update from(esl, qvsl, qvsi)
    !$omp target update from(thl0h, qt0h, ql0h, ql0)
    !$omp target update from(th0av, exnf, exnh, thvf, rhof, dzh)
    !$omp target update from(thv0, thvh, thvf, rhof)

    call timer_toc(routine)

  end subroutine thermodynamics

  !> Cleans up after the run
  subroutine exitthermodynamics
    !$acc exit data delete(th0av, thv0, thetah, qth, qlh)
    !$omp target exit data map(delete:th0av,thv0,thetah,qth,qlh)
    deallocate(th0av, thv0, thetah, qth, qlh)
  end subroutine exitthermodynamics

  !> Calculate real temperature tmp0 from thl0, for the dry case i.e. ql=0
  subroutine calc_dry_tmp
    integer :: i, j, k

    !$acc parallel loop collapse(3) default(present) async(1)
    !$omp target teams loop collapse(3) defaultmap(present:aggregate)&
    !$omp defaultmap(present:allocatable)
    do k = 1,k1
       do j = 2,j1
          do i = 2,i1
             tmp0(i,j,k) = exnf(k)*thl0(i,j,k)
          end do
       end do
    end do

  end subroutine calc_dry_tmp

  !> Calculate thetav and dthvdz
  !> also calculates and stores tmp0 for statistics
  subroutine calthv

    character(len=*), parameter :: routine = modname//'/calthv'

    integer i, j, k
    real(field_r)    qs
    real(field_r)    a_surf,b_surf,dq,dth,dthv,temp
    real(field_r)    a_dry, b_dry, a_moist, b_moist, c_liquid, epsilon, eps_I, chi_sat, chi
    real(field_r)    del_thv_sat, del_thv_dry

    call timer_tic(routine, 1)

    if (lmoist) then
      !$acc parallel loop collapse(3) default(present) async(1)
      !$omp target teams loop collapse(3) defaultmap(present:aggregate)&
      !$omp defaultmap(present:allocatable)
      do k = 2, k1
        do j = 2, j1
          do i = 2, i1
            thv0h(i,j,k) = calc_virt_pot_temp(thl0h(i,j,k), qt0h(i,j,k), &
                                              ql0h(i,j,k), exnh(k))
          end do
        end do
      end do

      !TODO: fix the branching in this loop
      !$acc parallel loop collapse(3) default(present) &
      !$acc private(a_dry, b_dry, a_moist, b_moist, c_liquid, epsilon, eps_I, &
      !$acc         chi_sat, chi, dthv, del_thv_dry, del_thv_sat, temp, qs, dq, dth) &
      !$acc async(1)
      !$omp target teams loop private(a_dry,b_dry,a_moist,b_moist,c_liquid,&
      !$omp epsilon,eps_i,chi_sat,chi,dthv,del_thv_dry,del_thv_sat,temp,qs,&
      !$omp dq,dth) collapse(3) defaultmap(present:aggregate)&
      !$omp defaultmap(present:allocatable)
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
               tmp0(i,j,k) = temp !stored for statistics
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
            else
                tmp0(i,j,k) = thl0(i,j,k)*exnf(k) !stored for statistics
            end if

            dthvdz(i,j,k) = dthv/(dzh(k+1)+dzh(k))
          end do
        end do
      end do

      !$acc parallel loop collapse(2) default(present) private(temp, qs, a_surf, b_surf) async(1)
      !$omp target teams loop private(temp,qs,a_surf,b_surf) collapse(2)&
      !$omp defaultmap(present:aggregate) defaultmap(present:allocatable)
      do j=2,j1
        do i=2,i1
          if(ql0(i,j,1)>0) then
            temp = thl0(i,j,1)*exnf(1)+(rlv/cp)*ql0(i,j,1)
            tmp0(i,j,1) = temp !stored for statistics
            qs   = qt0(i,j,1) - ql0(i,j,1)
            a_surf   = (1-qt0(i,j,1)+rv/rd*qs*(1+rlv/(rv*temp))) &
                      /(1+rlv**2*qs/(cp*rv*temp**2))
            b_surf   = a_surf*rlv/(temp*cp)-1

          else
            tmp0(i,j,1) = thl0(i,j,1)*exnf(1) !stored for statistics
            a_surf = 1+(rv/rd-1)*qt0(i,j,1)
            b_surf = rv/rd-1

          end if
          dthvdz(i,j,1) = a_surf*dthldz(i,j) + b_surf*thl0(i,j,1)*dqtdz(i,j)
        end do
      end do

    else
      !$acc parallel loop collapse(3) default(present) async(1)
       !$omp target teams loop collapse(3) defaultmap(present:aggregate)&
       !$omp defaultmap(present:allocatable)
      do k = 2, k1
        do j = 2, j1
          do i = 2, i1
            thv0h(i,j,k)  = thl0h(i,j,k)
          end do
        end do
      end do

      !$acc parallel loop collapse(3) default(present) async(1)
      !$omp target teams loop collapse(3) defaultmap(present:aggregate)&
      !$omp defaultmap(present:allocatable)
      do k = 2, kmax
        do j = 2, j1
          do i = 2, i1
            dthvdz(i,j,k) = (thl0(i,j,k+1)-thl0(i,j,k-1))/(dzh(k+1)+dzh(k))
          end do
        end do
      end do

      !$acc parallel loop collapse(2) default(present) async(1)
      !$omp target teams loop collapse(2) defaultmap(present:aggregate)&
      !$omp defaultmap(present:allocatable)
      do j = 2, j1
        do i = 2, i1
          dthvdz(i,j,1) = dthldz(i,j)
        end do
      end do
    end if

    !$acc parallel loop collapse(3) default(present) async(1)
    !$omp target teams loop collapse(3) defaultmap(present:aggregate)&
    !$omp defaultmap(present:allocatable)
    do k = 1, kmax
      do j = 2, j1
        do i = 2, i1
          if(abs(dthvdz(i,j,k)) < eps1) then
            dthvdz(i,j,k) = sign(eps1, dthvdz(i,j,k))
          end if
        end do
      end do
    end do

    call timer_toc(routine)

  end subroutine calthv

  !> Diagnones slab averaged fields assuming hydrostatic equilibrium.
  subroutine diagfld

    character(len=*), parameter :: routine = modname//'/diagfld'

    integer :: k,n

    call timer_tic(routine, 1)

    ! 1. Compute slab averaged fields

    !$acc parallel loop gang(static:1) default(present) async wait(1)
    !$omp target teams loop defaultmap(present:aggregate)&
    !$omp defaultmap(present:allocatable)
    do k = 1, k1
      u0av(k) = 0.0_field_r
    end do

    !$acc parallel loop gang(static:1) default(present) async wait(1)
    !$omp target teams loop defaultmap(present:aggregate)&
    !$omp defaultmap(present:allocatable)
    do k = 1, k1
      v0av(k) = 0.0_field_r
    end do

    !$acc parallel loop gang(static:1) default(present) async wait(1)
    !$omp target teams loop defaultmap(present:aggregate)&
    !$omp defaultmap(present:allocatable)
    do k = 1, k1
      thl0av(k) = 0.0_field_r
    end do

    !$acc parallel loop gang(static:1) default(present) async wait(1)
    !$omp target teams loop defaultmap(present:aggregate)&
    !$omp defaultmap(present:allocatable)
    do k = 1, k1
      th0av(k) = 0.0_field_r
    end do

    !$acc parallel loop gang(static:1) default(present) async wait(1)
    !$omp target teams loop defaultmap(present:aggregate)&
    !$omp defaultmap(present:allocatable)
    do k = 1, k1
      qt0av(k) = 0.0_field_r
    end do

    !$acc parallel loop gang(static:1) default(present) async wait(1)
    !$omp target teams loop defaultmap(present:aggregate)&
    !$omp defaultmap(present:allocatable)
    do k = 1, k1
      ql0av(k) = 0.0_field_r
    end do

    !$acc parallel loop gang vector collapse(2) default(present) async wait(1)
    !$omp target teams loop collapse(2)
    do k = 1, k1
      do n = 1, nsv
        sv0av(k,n) = 0.0_field_r
      end do
    end do

    !$acc wait

    ! If the IBM is enabled, exclude the building cells from the averages
    if (.not. lapply_ibm) then
      call slabavg_gpu(u0,ih,u0av)
      call slabavg_gpu(v0,ih,v0av)
      call slabavg_gpu(thl0,ih,thl0av)
      call slabavg_gpu(qt0,ih,qt0av)
      call slabavg_gpu(ql0,ih,ql0av)
      do n=1,nsv
        call slabavg_gpu(sv0(:,:,:,n),ih,sv0av(:,n))
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

    if ((timee < 0.01 .or. .not. lconstexner) .and. .not. lbaseexner) then
      !$acc parallel loop gang(static:1) default(present)
      !$omp target teams loop defaultmap(present:aggregate)&
      !$omp defaultmap(present:allocatable)
      do k = 1, k1
        exnf(k) = 1 - grav * zf(k) / (cp * thls)
        exnh(k) = 1 - grav * zh(k) / (cp * thls)
      end do
    endif

    !$acc parallel loop gang(static:1) default(present) async(1)
    !$omp target teams loop defaultmap(present:aggregate)&
    !$omp defaultmap(present:allocatable)
    do k = 1, k1
      th0av(k) = thl0av(k) + (rlv / cp) * ql0av(k) / exnf(k)
    end do

    ! 2. Calculate the pressure profiles assuming hydrostatic equilibrium.

    ! 2.1 Use first guess of theta, then recalculate theta

    !$omp target update to(dzf,dzh)
    call fromztop

    !$acc parallel loop gang(static:1) default(present) async(1)
    !$omp target teams loop defaultmap(present:allocatable)
    do k = 1, k1
      th0av(k) = thl0av(k) + (rlv / cp) * ql0av(k) / exnf(k)
    end do

    if ((timee < 0.01 .or. .not. lconstexner) .and. .not. lbaseexner) then
      !$acc parallel loop gang(static:1) default(present) async(1)
      !$omp target teams loop defaultmap(present:allocatable)
      do k = 1, k1
        exnf(k) = (presf(k) / pref0)**(rd / cp)
      end do
    end if

    ! 2.2 Use new updated value of theta for determination of pressure
    call fromztop

    ! 3. Construct density profiles and exner function

    if ((timee < 0.01 .or. .not. lconstexner) .and. .not. lbaseexner) then
      !$acc serial default(present) async(1)
      !$omp target defaultmap(present:allocatable)
      exnh(1) = (ps/pref0)**(rd/cp)
      exnf(1) = (presf(1)/pref0)**(rd/cp)
      !$acc end serial
      !$omp end target

      !$acc parallel loop default(present) async(1)
      !$omp target teams loop defaultmap(present:aggregate)&
      !$omp defaultmap(present:allocatable)
      do k=2,k1
        exnf(k) = (presf(k)/pref0)**(rd/cp)
        exnh(k) = (presh(k)/pref0)**(rd/cp)
      end do
    endif

    !$acc parallel loop default(present) async(1)
    !$omp target teams loop defaultmap(present:aggregate)&
    !$omp defaultmap(present:allocatable)
    do k=1,k1
      thvf(k) = th0av(k)*exnf(k)*(1+(rv/rd-1)*qt0av(k)-rv/rd*ql0av(k))
      rhof(k) = presf(k)/(rd*thvf(k))
    end do

    call timer_toc(routine)

  end subroutine diagfld

  !> Calculates slab averaged pressure.
  subroutine fromztop

    character(len=*), parameter :: routine = modname//'/fromztop'

    integer   k
    real(field_r)  rdocp

    call timer_tic(routine, 2)

    rdocp = rd/cp

    ! Interpolate theta and qt to half levels

    !$acc parallel loop default(present) async(1)
    !$omp target teams loop defaultmap(present:allocatable)
    do k=2,k1
      thetah(k) = (th0av(k)*dzf(k-1) + th0av(k-1)*dzf(k))/(2*dzh(k))
      qth   (k) = (qt0av(k)*dzf(k-1) + qt0av(k-1)*dzf(k))/(2*dzh(k))
      qlh   (k) = (ql0av(k)*dzf(k-1) + ql0av(k-1)*dzf(k))/(2*dzh(k))
    end do

    ! Calculate pressures at full levels
    ! Do this on the CPU for now; these loops are serial so GPU is very slow!

    !$acc update self(thetah, qth, qlh, th0av, qt0av, ql0av) async(1)
    !$omp target update from(thetah, qth, qlh, th0av, qt0av, ql0av)
    !$acc wait

    thvh(1) = th0av(1)*(1+(rv/rd-1)*qt0av(1)-rv/rd*ql0av(1))
    presf(1) = ps**rdocp - grav*(pref0**rdocp)*zf(1) /(cp*thvh(1))
    presf(1) = presf(1)**(1/rdocp)

    do k=2,k1
      thvh(k)  = thetah(k)*(1+(rv/rd-1)*qth(k)-rv/rd*qlh(k))
      presf(k) = presf(k-1)**rdocp - &
                     grav*(pref0**rdocp)*dzh(k) /(cp*thvh(k))
      presf(k) = presf(k)**(1/rdocp)
    end do

    ! Calculate pressures at half levels

    presh(1) = ps
    thvf(1) = th0av(1)*(1+(rv/rd-1)*qt0av(1)-rv/rd*ql0av(1))

    do k=2,k1
      thvf(k)  = th0av(k)*(1+(rv/rd-1)*qt0av(k)-rv/rd*ql0av(k))
      presh(k) = presh(k-1)**rdocp - &
                     grav*(pref0**rdocp)*dzf(k-1) / (cp*thvf(k-1))
      presh(k) = presh(k)**(1/rdocp)
    end do

    !$acc update device(thvh, presf, thvf, presh) async(1)
    !$omp target update to(thvh, presf, thvf, presh)

    call timer_toc(routine)

  end subroutine fromztop

  !> Magnus formulas for q_sat over liquid and ice.
  !!
  !! from Huang 2018 https://doi.org/10.1175/JAMC-D-17-0334.
  !! Warning: for performance, check that rd/rv etc are pre-computed
  pure function qsat_magnus(T, p) result(qsat)
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

  !> Huang's formulas for q_sat over liquid and ice.
  !!
  !! from Huang 2018 https://doi.org/10.1175/JAMC-D-17-0334.
  !! should be more accurate than Magnus, at the cost of more divisions
  !! Warning: for performance, check that rd/rv etc are pre-computed
  pure function qsat_huang(T, p) result(qsat)
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

  !> Compute the saturation vapor pressure via table lookup.
  pure function esat_tab(T) result(es)
    !$omp declare target
    !$acc routine seq
    real(field_r), intent(in) :: T
    integer :: tlo
    real(field_r) :: es, interp_w

    ! interpolated ice-liquid saturation vapor pressure from table
    ! note if imicto==imicro_bulk3, the table is for liquid only
    interp_w = (T - 150) * 5
    tlo = int(interp_w)
    interp_w = interp_w - tlo
    es = (1 - interp_w) * esatmtab(tlo) + interp_w * esatmtab(tlo+1)
  end function esat_tab

  !> Computes the saturation specific humidity via table lookup.
  pure function qsat_tab(T, p) result(qsat)
    !$omp declare target
    !$acc routine seq
    real(field_r), intent(in) :: T, p
    real(field_r) :: qsat
    integer :: tlo
    real(field_r) :: es, interp_w

    interp_w = (T - 150) * 5
    tlo = int(interp_w)
    interp_w = interp_w - tlo
    es = (1 - interp_w) * esatmtab(tlo) + interp_w * esatmtab(tlo+1)

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

    !$acc parallel loop gang default(present) async(stream) &
    !$acc private(b, qli, qsat, qti, Tl)
    !$omp target teams loop private(b,qli,qsat,qti,tl)&
    !$omp defaultmap(present:aggregate) defaultmap(present:allocatable)
    do k = 1, k1
      ! Find lowest thl and highest qt in the slab.
      ! If they in combination are not saturated, the whole slab is below saturation.
      !
      ! TODO: on GPU, test if it's cheaper to just do the computation instead.
      TL_min = minval(thl(2:i1,2:j1,k)) * exn(k)
      qt_max = maxval(qt(2:i1,2:j1,k))
      qsat = qsat_tab(TL_min, pres(k))
      if (qt_max > qsat) then
        !$acc loop vector collapse(2)
        !$omp loop collapse(2)
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
      else
        ql(:,:,k) = 0
      end if
    end do

    call timer_toc(routine)

  end subroutine saturation_adjustment

#if defined(DALES_GPU)
  !> Compute the cloud water content via the saturation adjustment method.
  subroutine saturation_adjustment_gpu(qt, thl, pres, exn, ql, opt_stream)

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
    ! FIXME: GPU divergence
    !!$omp target teams loop private(b,qli,qsat,qti,tl) collapse(3)&
    !!$omp defaultmap(present:aggregate) defaultmap(present:allocatable)
    do k = 1, k1
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
    end do
    !$omp target update to(ql)

    call timer_toc(routine)

  end subroutine saturation_adjustment_gpu
#endif

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

    real(field_r) :: esi      !< Saturation vapor pressure for ice (not stored) [Pa]
    real(field_r) :: qsat     !< Saturation specific humidity [kg/kg]
    real(field_r) :: T        !< Temperature [K]
    real(field_r) :: interp_w !< Interpolation temperature [K]
    integer       :: tlo      !< Index of temperature in esat lookuptable

    !$acc parallel loop collapse(3) default(present) async(1) &
    !$acc private(qsat, T, interp_w, tlo, esi)
    !$omp target teams loop private(qsat,t,interp_w,tlo,esi) collapse(3)&
    !$omp defaultmap(present:aggregate) defaultmap(present:allocatable)
    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          qsat = max(qt(i,j,k) - ql(i,j,k), 1.0_field_r)
          T = exn(k) * thl(i,j,k) + (rlv / cp) * ql(i,j,k)

          interp_w = (T - 150) * 5
          tlo = int(interp_w)
          interp_w = interp_w - tlo

          ! Liquid
          esl(i,j,k) = (1 - interp_w) * esatltab(tlo) &
                       + interp_w * esatltab(tlo + 1)
          qvsl(i,j,k) = rd / rv * esl(i,j,k) &
                        / (pres(k) - (1 - rd / rv) * esl(i,j,k))

          ! Ice
          esi = (1 - interp_w) * esatitab(tlo) &
                + interp_w * esatitab(tlo + 1)
          qvsi(i,j,k) = rd / rv * esi / (pres(k) - (1 - rd / rv) * esi)
        end do
      end do
    end do

  end subroutine calc_saturation_humidities

  !> Interpolate a scalar field to the half levels.
  !!
  !! CJ: this is a very general subroutine which could be in its own module.
  !!     It's used only here and in modstartup, so I'll leave it here for now.
  subroutine calc_halflev(phi, dzf, dzhi, phi_surf, lkappa, phi_half, &
                          opt_stream)

    real(field_r), intent(in) :: phi(2-ih:,2-jh:,:) !< Scalar field to interpolate [*]
    real(field_r), intent(in) :: dzf(:)             !< Thickness of full model levels [m]
    real(field_r), intent(in) :: dzhi(:)            !< Inverse of thickness of half model levels [m]
    real(field_r), intent(in) :: phi_surf           !< Value to use at the surface level [*]
    logical,       intent(in) :: lkappa             !< Interpolate using Kappa scheme

    real(field_r), intent(out) :: phi_half(2-ih:,2-jh:,:) !< Scalar field on half levels [*]

    integer, intent(in), optional :: opt_stream !< (Optional) OpenACC stream ID.

    integer :: i, j, k
    integer :: stream = 1

    if (present(opt_stream)) stream = opt_stream

    if (lkappa) then
      call halflev_kappa(phi, phi_half)
    else
      !$acc parallel loop collapse(3) default(present) async(stream)
      !$omp target teams loop collapse(3) defaultmap(present:aggregate)&
      !$omp defaultmap(present:allocatable)
      do k = 2, k1
        do j = 2, j1
          do i = 2, i1
            phi_half(i,j,k) = (phi(i,j,k) * dzf(k-1) + phi(i,j,k-1) * dzf(k)) &
                              * (0.5_field_r * dzhi(k))
          end do
        end do
      end do
    end if

    !$acc parallel loop collapse(2) default(present) async(stream)
    !$omp target teams loop collapse(2) defaultmap(present:aggregate)&
    !$omp defaultmap(present:allocatable)
    do j = 2, j1
      do i = 2, i1
        phi_half(i,j,1) = phi_surf
      end do
    end do

  end subroutine calc_halflev

end module modthermodynamics
