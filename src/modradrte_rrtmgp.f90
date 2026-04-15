!> \file modradrte_rrtmgp.f90
!!  Interfaces with the radiation library RTE-RRTMGP from Earth System Radiation group

!>
!!  Interfaces with the radiation library RTE-RRTMGP from Earth System Radiation group
!>
!!  \author Laurent Soucasse, Netherlands eScience Center
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
!  Copyright 2023 Netherlands eScience Center
!
module modradrte_rrtmgp
  use modraddata
  use modprecision, only : field_r
  use modtimer
  use modlogging, only: finish
  use modrrtmgp_utils,        only: load_gas_optics, load_cloud_optics, &
                                    stop_on_err
  ! RTE-RRTMGP modules
  use mo_optical_props,       only: ty_optical_props, &
                                    ty_optical_props_arry, &
                                    ty_optical_props_1scl, ty_optical_props_2str
  use mo_gas_optics_rrtmgp,   only: ty_gas_optics_rrtmgp
  use mo_cloud_optics_rrtmgp, only: ty_cloud_optics_rrtmgp
  use mo_source_functions,    only: ty_source_func_lw
  use mo_fluxes,              only: ty_fluxes_broadband
  use mo_gas_concentrations,  only: ty_gas_concs
  use mo_rte_kind,            only: wl

  implicit none

  character(len=*), parameter :: modname = 'modradrte_rrtmgp'

  private
  ! RRTMGP variables
  type(ty_gas_concs)                        :: gas_concs
  type(ty_source_func_lw), save             :: sources_lw
  type(ty_gas_optics_rrtmgp)                :: k_dist_lw, k_dist_sw
  type(ty_cloud_optics_rrtmgp)              :: cloud_optics_lw, cloud_optics_sw
  class(ty_optical_props_arry), allocatable :: atmos_lw, atmos_sw, clouds_lw, clouds_sw
  type(ty_fluxes_broadband)                 :: fluxes_lw, fluxes_sw, fluxes_cs_lw, fluxes_cs_sw
  real(kind=kind_rb), dimension(:,:), allocatable :: inc_sw_flux, sfc_alb_dir, sfc_alb_dif
  !Specify gas names, the first five (h2o, o3, co2, ch4 and n2o) are mandatory as they are major absorbers
  integer, parameter                        :: ngas = 10
  character(len=5), dimension(ngas)         :: gas_names = ['h2o  ', 'o3   ', 'co2  ', 'ch4  ', 'n2o  ', 'o2   ', 'cfc11', 'cfc12', 'cfc22', 'ccl4 ']
  integer                                   :: nlay, nlev, ncol, nbndlw, nbndsw, ngptsw
  logical                                   :: initialized = .false.

  public :: radrte_rrtmgp, exit_radrte_rrtmgp

contains

  subroutine init_radrte_rrtmgp
    use mo_rte_config,         only: rte_config_checks

    ! DALES modules
    use modradrrtmg,           only: readSounding, readTraceProfs
    use modmpi,                only: myid
    use modfields,             only: initial_presh, initial_presf
    use modglobal,             only: imax, jmax, kmax, k1
    implicit none

    character(len=*), parameter :: routine = modname//'/init_radrte_rrtmgp'

    integer                 :: k, npatch, ierr(3)=0
    character(len=256)      :: k_dist_file_lw = "rrtmgp-gas-lw-g128.nc"
    character(len=256)      :: k_dist_file_sw = "rrtmgp-gas-sw-g112.nc"
    character(len=256)      :: cloud_optics_file_lw = "rrtmgp-clouds-lw-bnd.nc"
    character(len=256)      :: cloud_optics_file_sw = "rrtmgp-clouds-sw-bnd.nc"

    ! Reading sounding (patch above Dales domain), only once
    call readSounding(initial_presh(k1)/100.,npatch_start,npatch_end)

    if(npatch_end.ne.npatch_start) then
      npatch = npatch_end - npatch_start + 1
    else
      call finish(routine,  'ERROR: No valid radiation sounding found above the LES domain, check sounding input file')
    end if

    !old notations nlay-1=kradmax=nzrad, nlay=krad1, nlev=krad2
    nlay = kmax + npatch + 1
    nlev = nlay + 1 ! necessary?
    !the indices below are necessary for the readTraceProfs routine
    kradmax=nlay-1
    krad1=nlay
    krad2=nlev

    !Set the default value of nbatch if not provided in nameoptions
    if(nbatch==0) nbatch = jmax
    !Check if jmax is a mutliple of nbatch, if user provided
    if(mod(jmax,nbatch)/=0) call finish(routine, 'ERROR: Wrong batch number specified')
    ncol = imax*jmax/nbatch

    ! Allocating working variables
    allocate(layerP(ncol,nlay), &
             layerT(ncol,nlay), &
             h2ovmr(ncol,nlay), &
             tracevmr(ncol,nlay), &
             liquidRe(ncol,nlay), &
             iceRe(ncol,nlay), &
             LWP_slice(ncol,nlay), &
             IWP_slice(ncol,nlay), &
             tg_slice(ncol), &
             presf_input(nlay-1), &
             solarZenithAngleCos(ncol), &
             STAT=ierr(1))
    allocate(interfaceP(ncol,nlay+1), &
             interfaceT(ncol,nlay+1), &
             lwUp_slice(ncol,nlay+1), &
             lwDown_slice(ncol,nlay+1), &
             swUp_slice(ncol,nlay+1), &
             swDown_slice(ncol,nlay+1), &
             swDownDir_slice(ncol,nlay+1), &
             presh_input(nlay), &
             STAT=ierr(2))
    if(doclearsky) then
      allocate(lwUpCS_slice(ncol,nlay+1), &
               lwDownCS_slice(ncol,nlay+1), &
               swUpCS_slice(ncol,nlay+1), &
               swDownCS_slice(ncol,nlay+1), &
               STAT=ierr(3))
    endif
    if(any(ierr(:)/=0)) then
      call finish(routine, 'ERROR: Could not allocate input/output arrays for radiation variables')
    end if

    ! Pressure, trace gases and sounding (patch above the DALES domain) initialization
    ! Patch sounding profile pressures above domain pressures (convert to hPa!)
    presf_input(1:kmax)   = initial_presf(1:kmax)  /100.
    presh_input(1:k1)     = initial_presh(1:k1)/100.
    if(npatch>0) then
      presf_input(k1  :kradmax) = psnd(npatch_start:npatch_end)
      presh_input(k1+1:kradmax) = 0.5*(psnd(npatch_start:npatch_end-1) &
                                     + psnd(npatch_start+1:npatch_end))
      presh_input(krad1) = max(0.5*psnd(npatch_end),            &
                               1.5*psnd(npatch_end)-0.5*psnd(npatch_end-1))
    end if

    ! Set up pressure layer and interface values (pressures in SI units, i.e. Pa)
    do k=1,nlay-1
      layerP(:,k) = presf_input(k)*100.00
    enddo
    do k=1,nlay
      interfaceP(:,k) = presh_input(k)*100.00
    enddo
    layerP(:,nlay) = 0.5*presh_input(nlay)*100.00
    interfaceP(:, nlay+1) = min(1.e-4_kind_rb , 0.25*layerP(1,nlay))

    ! Set up temperature and h2o layer values above the DALES domain
    do k=1,nlay-kmax-1
      layerT(:,kmax+k) = tsnd(npatch_start+k-1)
      h2ovmr(:,kmax+k) = mwdry/mwh2o * qsnd(npatch_start+k-1)
    enddo
    h2ovmr(:,nlay) = h2ovmr(:,nlay-1)
    layerT(:,nlay) = 2.*layerT(:,nlay-1)-layerT(:, nlay-2)

    ! Reading Trace Profiles
    call readTraceProfs
    if(myid==0) write(*,*) 'Trace gas profile have been read'

    ! Specific RRTMGP initialization
    call stop_on_err(gas_concs%init(gas_names))
    !setup trace gases concentration once for all
    !it seems the array used by the set_vmr function has to be on the GPU... to be tested
    !$acc data create(tracevmr)
!!$omp target data map(alloc:tracevmr)
    do k=1,nlay; tracevmr(:,k) = o3(k); enddo
    !$acc update device(tracevmr)
!!$omp target update to(tracevmr)
    call stop_on_err(gas_concs%set_vmr(trim(gas_names(2)), tracevmr))
    do k=1,nlay; tracevmr(:,k) = co2(k); enddo
    !$acc update device(tracevmr)
!!$omp target update to(tracevmr)
    call stop_on_err(gas_concs%set_vmr(trim(gas_names(3)), tracevmr))
    do k=1,nlay; tracevmr(:,k) = ch4(k); enddo
    !$acc update device(tracevmr)
!!$omp target update to(tracevmr)
    call stop_on_err(gas_concs%set_vmr(trim(gas_names(4)), tracevmr))
    do k=1,nlay; tracevmr(:,k) = n2o(k); enddo
    !$acc update device(tracevmr)
!!$omp target update to(tracevmr)
    call stop_on_err(gas_concs%set_vmr(trim(gas_names(5)), tracevmr))
    do k=1,nlay; tracevmr(:,k) = o2(k); enddo
    !$acc update device(tracevmr)
!!$omp target update to(tracevmr)
    call stop_on_err(gas_concs%set_vmr(trim(gas_names(6)), tracevmr))
    do k=1,nlay; tracevmr(:,k) = cfc11(k); enddo
    !$acc update device(tracevmr)
!!$omp target update to(tracevmr)
    call stop_on_err(gas_concs%set_vmr(trim(gas_names(7)), tracevmr))
    do k=1,nlay; tracevmr(:,k) = cfc12(k); enddo
    !$acc update device(tracevmr)
!!$omp target update to(tracevmr)
    call stop_on_err(gas_concs%set_vmr(trim(gas_names(8)), tracevmr))
    do k=1,nlay; tracevmr(:,k) = cfc22(k); enddo
    !$acc update device(tracevmr)
!!$omp target update to(tracevmr)
    call stop_on_err(gas_concs%set_vmr(trim(gas_names(9)), tracevmr))
    do k=1,nlay; tracevmr(:,k) = ccl4(k); enddo
    !$acc update device(tracevmr)
!!$omp target update to(tracevmr)
    call stop_on_err(gas_concs%set_vmr(trim(gas_names(10)), tracevmr))
    !$acc end data
!!$omp end target data
    deallocate(tracevmr)

    !$acc enter data copyin(layerP, layerT, interfaceP, h2ovmr)
!$omp target enter data map(to:layerp,layert,interfacep,h2ovmr)
    !$acc enter data create(interfaceT, tg_slice)
!$omp target enter data map(alloc:interfacet,tg_slice)
    !$acc enter data create(lwUp_slice, lwDown_slice)
!$omp target enter data map(alloc:lwup_slice,lwdown_slice)
    !$acc enter data create(swUp_slice, swDown_slice, swDownDir_slice)
!$omp target enter data map(alloc:swup_slice,swdown_slice,&
!$omp swdowndir_slice)
    !$acc enter data create(solarZenithAngleCos)
!$omp target enter data map(alloc:solarzenithanglecos)
    !$acc enter data create(liquidRe, iceRe, LWP_slice, IWP_slice)
!$omp target enter data map(alloc:liquidre,icere,lwp_slice,iwp_slice)
    if(doclearsky) then
      !$acc enter data create(lwUpCS_slice, lwDownCS_slice)
!$omp target enter data map(alloc:lwupcs_slice,lwdowncs_slice)
      !$acc enter data create(swUpCS_slice, swDownCS_slice)
!$omp target enter data map(alloc:swupcs_slice,swdowncs_slice)
    endif

    ! Longwave init
    if(rad_longw) then

      ! Load k distributions
      call load_gas_optics(k_dist_lw, k_dist_file_lw, gas_concs)
      if(.not. k_dist_lw%source_is_internal()) &
        call finish(routine, "k-distribution file isn't LW")
      nbndlw = k_dist_lw%get_nband()

      ! Initialize gas optical properties
      allocate(ty_optical_props_1scl::atmos_lw)
      select type(atmos_lw)
        class is (ty_optical_props_1scl)
          call stop_on_err(atmos_lw%alloc_1scl(ncol, nlay, k_dist_lw))
          !$acc enter data copyin(atmos_lw) create(atmos_lw%tau)
!$omp target enter data map(to:atmos_lw) map(alloc:atmos_lw%tau)
      end select

      ! Load cloud property data
      call load_cloud_optics(cloud_optics_lw, cloud_optics_file_lw)
      call stop_on_err(cloud_optics_lw%set_ice_roughness(2))

      ! Initialize cloud optical properties
      allocate(ty_optical_props_1scl::clouds_lw)
      call stop_on_err(clouds_lw%init(k_dist_lw%get_band_lims_wavenumber()))
      select type(clouds_lw)
        class is (ty_optical_props_1scl)
          call stop_on_err(clouds_lw%alloc_1scl(ncol, nlay))
          !$acc enter data copyin(clouds_lw) create(clouds_lw%tau)
!$omp target enter data map(to:clouds_lw) map(alloc:clouds_lw%tau)
      end select

      ! Allocate source term and define emissivity
      call stop_on_err(sources_lw%alloc(ncol, nlay, k_dist_lw))
      allocate(emis(nbndlw,ncol))
      ! we set emis from the surface data now, so emis is initialized later on
      !$acc enter data copyin(emis, sources_lw)
!$omp target enter data map(to:emis,sources_lw)
      !$acc enter data create(sources_lw%lay_source, sources_lw%lev_source, &
      !$acc&                  sources_lw%sfc_source, sources_lw%sfc_source_Jac)
!$omp target enter data map(alloc:sources_lw%lay_source,&
!$omp sources_lw%lev_source,sources_lw%sfc_source,&
!$omp sources_lw%sfc_source_jac)

      ! Define lw fluxes pointers
      fluxes_lw%flux_up => lwUp_slice(:,:)
      fluxes_lw%flux_dn => lwDown_slice(:,:)
      if(doclearsky) then
        fluxes_cs_lw%flux_up => lwUpCS_slice(:,:)
        fluxes_cs_lw%flux_dn => lwDownCS_slice(:,:)
      endif
    endif

    ! Shortwave init
    if(rad_shortw) then

      ! Load k distributions
      call load_gas_optics(k_dist_sw, k_dist_file_sw, gas_concs)
      if(k_dist_sw%source_is_internal()) &
        call finish(routine, "k-distribution file isn't SW")
      nbndsw = k_dist_sw%get_nband()
      ngptsw = k_dist_sw%get_ngpt()

      ! Initialize gas optical properties
      allocate(ty_optical_props_2str::atmos_sw)
      select type(atmos_sw)
        class is (ty_optical_props_2str)
          call stop_on_err(atmos_sw%alloc_2str(ncol, nlay, k_dist_sw))
          !$acc enter data copyin(atmos_sw) create(atmos_sw%tau, atmos_sw%ssa, atmos_sw%g)
!$omp target enter data map(to:atmos_sw) map(alloc:atmos_sw%tau,&
!$omp atmos_sw%ssa,atmos_sw%g)
      end select

      ! Load cloud property data
      call load_cloud_optics(cloud_optics_sw, cloud_optics_file_sw)
      call stop_on_err(cloud_optics_sw%set_ice_roughness(2))

      ! Initialize cloud optical properties
      allocate(ty_optical_props_2str::clouds_sw)
      call stop_on_err(clouds_sw%init(k_dist_sw%get_band_lims_wavenumber()))
      select type(clouds_sw)
        class is (ty_optical_props_2str)
          call stop_on_err(clouds_sw%alloc_2str(ncol, nlay))
          !$acc enter data copyin(clouds_sw) create(clouds_sw%tau, clouds_sw%ssa, clouds_sw%g)
!$omp target enter data map(to:clouds_sw) map(alloc:clouds_sw%tau,&
!$omp clouds_sw%ssa,clouds_sw%g)
      end select

      ! Define boundary conditions
      allocate(inc_sw_flux(ncol,ngptsw))
      allocate(sfc_alb_dir(nbndsw,ncol), sfc_alb_dif(nbndsw,ncol))
      !$acc enter data create(inc_sw_flux, sfc_alb_dir, sfc_alb_dif)
!$omp target enter data map(alloc:inc_sw_flux,sfc_alb_dir,sfc_alb_dif)

      fluxes_sw%flux_up => swUp_slice(:,:)
      fluxes_sw%flux_dn => swDown_slice(:,:)
      fluxes_sw%flux_dn_dir => swDownDir_slice(:,:)
      if(doclearsky) then
        fluxes_cs_sw%flux_up => swUpCS_slice(:,:)
        fluxes_cs_sw%flux_dn => swDownCS_slice(:,:)
      endif
    endif

    call rte_config_checks(.false._wl)

    initialized = .true.

  end subroutine init_radrte_rrtmgp

  subroutine radrte_rrtmgp
    use mo_rte_lw,             only: rte_lw
    use mo_rte_sw,             only: rte_sw
    implicit none

    logical                 :: sunUp = .false.
    integer                 :: ibatch

    if(.not.initialized) call init_radrte_rrtmgp

    do ibatch = 1, nbatch

      call setupColumnProfiles(ibatch)

      if(rad_longw) then
        call setupEmis(ibatch)
        ! Compute optical properties and source
        call timer_tic('modradrte_rrtmgp/lwgasoptics', 0)
        call stop_on_err(k_dist_lw%gas_optics(layerP, interfaceP, & ! p_lay, p_lev (in, Pa)
                                              layerT, tg_slice, & ! t_lay, t_sfc (in, K)
                                              gas_concs, & ! gas volume mixing ratios (in)
                                              atmos_lw, & ! Optical properties (inout)
                                              sources_lw, & ! Planck source (inout)
                                              tlev = interfaceT)) ! t_lev (optional input, K)
        call timer_toc('modradrte_rrtmgp/lwgasoptics')

        ! Solve clear sky radiation transport if required
        if(doclearsky) then
          call stop_on_err(rte_lw(atmos_lw, & ! optical properties (in)
                                  sources_lw, & ! source function (in)
                                  emis, & ! emissivity at surface (in)
                                  fluxes_cs_lw)) ! fluxes (W/m2, inout)
        endif

        ! Compute and add cloud properties
        call timer_tic('modradrte_rrtmgp/lwcloudsoptics', 0)
        call stop_on_err(cloud_optics_lw%cloud_optics(LWP_slice, & ! cloud liquid water path (in, g/m2)
                                                      IWP_slice, & ! cloud ice water path (in, g/m2)
                                                      liquidRe, & ! cloud liquid particle effective size (in, microns)
                                                      iceRe, & ! cloud ice particle effective radius (in, microns)
                                                      clouds_lw)) ! cloud optical properties lw (inout)
        call stop_on_err(clouds_lw%increment(atmos_lw))
        call timer_toc('modradrte_rrtmgp/lwcloudsoptics')

        ! Solve radiation transport
        call timer_tic('modradrte_rrtmgp/lwrtesolve', 0)
        call stop_on_err(rte_lw(atmos_lw, & ! optical properties (in)
                                sources_lw, & ! source function (in)
                                emis, & ! emissivity at surface (in)
                                fluxes_lw)) ! fluxes (W/m2, inout)
        call timer_toc('modradrte_rrtmgp/lwrtesolve')
      endif

      if(rad_shortw) then

        ! setup incoming flux and albedo as a function of the zenith angle
        call setupSW(sunUp, ibatch)

        if(sunUp) then
          ! Compute optical properties and incoming shortwave flux
          call timer_tic('modradrte_rrtmgp/swgasoptics', 0)
          call stop_on_err(k_dist_sw%gas_optics(layerP, interfaceP, & ! p_lay, p_lev (in, Pa)
                                                layerT, & ! t_lay (in, K)
                                                gas_concs, & ! gas volume mixing ratios (in)
                                                atmos_sw, & ! Optical properties (inout)
                                                inc_sw_flux)) ! Incoming shortwave flux (inout)
          call timer_toc('modradrte_rrtmgp/swgasoptics')

          ! Solve clear sky radiation transport if required
          if(doclearsky) then
            call stop_on_err(rte_sw(atmos_sw, & ! optical properties (in)
                                    solarZenithAngleCos, & ! cosine of the solar zenith angle (in)
                                    inc_sw_flux, & ! solar incoming flux (in)
                                    sfc_alb_dir, sfc_alb_dif, & ! surface albedos, direct and diffuse (in)
                                    fluxes_cs_sw)) ! fluxes (inout, W/m2)
          endif

          ! Compute and add cloud properties
          call timer_tic('modradrte_rrtmgp/swcloudsoptics', 0)
          call stop_on_err(cloud_optics_sw%cloud_optics(LWP_slice, & ! cloud liquid water path (in, g/m2)
                                                        IWP_slice, & ! cloud ice water path (in, g/m2)
                                                        liquidRe, & ! cloud liquid particle effective size (in, microns)
                                                        iceRe, & ! cloud ice particle effective radius (in, microns)
                                                        clouds_sw)) ! cloud optical properties sw (inout)
          call stop_on_err(clouds_sw%delta_scale())
          call stop_on_err(clouds_sw%increment(atmos_sw))
          call timer_toc('modradrte_rrtmgp/swcloudsoptics')

          ! Solve radiation transport
          call timer_tic('modradrte_rrtmgp/swrtesolve', 0)
          call stop_on_err(rte_sw(atmos_sw, & ! optical properties (in)
                                  solarZenithAngleCos, & ! cosine of the solar zenith angle (in)
                                  inc_sw_flux, & ! solar incoming flux (in)
                                  sfc_alb_dir, sfc_alb_dif, & ! surface albedos, direct and diffuse (in)
                                  fluxes_sw)) ! fluxes (inout, W/m2)
          call timer_toc('modradrte_rrtmgp/swrtesolve')

        else
          call zero_SW_for_sunDown()
        end if

      endif

      call getFluxProfiles(ibatch)

    enddo

  end subroutine radrte_rrtmgp

  subroutine exit_radrte_rrtmgp
    implicit none

    if (rad_longw) then
      !$acc exit data delete(sources_lw)
!$omp target exit data map(delete:sources_lw)
      !$acc exit data delete(emis)
!$omp target exit data map(delete:emis)
      deallocate(emis)

      !$acc exit data delete(atmos_lw%tau) delete(atmos_lw)
!$omp target exit data map(delete:atmos_lw%tau,atmos_lw)
      !$acc exit data delete(clouds_lw%tau) delete(clouds_lw)
!$omp target exit data map(delete:clouds_lw%tau,clouds_lw)

    endif

    if (rad_shortw) then
      !$acc exit data delete(inc_sw_flux, sfc_alb_dir, sfc_alb_dif)
!$omp target exit data map(delete:inc_sw_flux,sfc_alb_dir,sfc_alb_dif)
      deallocate(inc_sw_flux, sfc_alb_dir, sfc_alb_dif)

      !$acc exit data delete(atmos_sw%tau) delete(atmos_sw)
!$omp target exit data map(delete:atmos_sw%tau,atmos_sw)
      !$acc exit data delete(clouds_sw%tau) delete(clouds_sw)
!$omp target exit data map(delete:clouds_sw%tau,clouds_sw)
    endif

    !$acc exit data delete(liquidRe, iceRe, LWP_slice, IWP_slice)
!$omp target exit data map(delete:liquidre,icere,lwp_slice,iwp_slice)
    !$acc exit data delete(solarZenithAngleCos)
!$omp target exit data map(delete:solarzenithanglecos)
    !$acc exit data delete(lwUp_slice, lwDown_slice)
!$omp target exit data map(delete:lwup_slice,lwdown_slice)
    !$acc exit data delete(swUp_slice, swDown_slice, swDownDir_slice)
!$omp target exit data map(delete:swup_slice,swdown_slice,&
!$omp swdowndir_slice)
    !$acc exit data delete(layerP, layerT, interfaceP, interfaceT, tg_slice, h2ovmr)
!$omp target exit data map(delete:layerp,layert,interfacep,interfacet,&
!$omp tg_slice,h2ovmr)
    if(doclearsky) then
      !$acc exit data delete(lwUpCS_slice, lwDownCS_slice, &
      !$acc&                 swUpCS_slice, swDownCS_slice)
!$omp target exit data map(delete:lwupcs_slice,lwdowncs_slice,&
!$omp swupcs_slice,swdowncs_slice)
    endif

    if(isAllocated_RadInputsOutputs) then
      deallocate(layerP, &
                 layerT, &
                 h2ovmr, &
                 tracevmr, &
                 liquidRe, &
                 iceRe, &
                 LWP_slice, &
                 IWP_slice, &
                 tg_slice, &
                 presf_input, &
                 solarZenithAngleCos)
      deallocate(interfaceP, &
                 interfaceT, &
                 lwUp_slice, &
                 lwDown_slice, &
                 swUp_slice, &
                 swDown_slice, &
                 swDownDir_slice, &
                 presh_input)
      if(doclearsky) then
        deallocate(lwUpCS_slice, &
                   lwDownCS_slice, &
                   swUpCS_slice, &
                   swDownCS_slice)
      endif
    end if

  end subroutine exit_radrte_rrtmgp
  subroutine zero_SW_for_sunDown
    implicit none
    integer :: i, k

    ! Make sure the SW output is 0 if sun is down.
    ! If the sun is down, rrtmgp does not initialize the sw fluxes arrays, so we need to set them to zero here to avoid uninitialized values in the output.
    !$acc parallel loop collapse(2) default(present)
!!$omp target teams loop collapse(2) defaultmap(present:aggregate)&
!!$omp defaultmap(present:allocatable)
    do i=1,ncol
      do k=1,nlay+1
        swUp_slice(i,k) = 0.
        swDown_slice(i,k) = 0.
        swDownDir_slice(i,k) = 0.
      enddo
    enddo

    if(doclearsky) then
    !$acc parallel loop collapse(2) default(present)
!!$omp target teams loop collapse(2) defaultmap(present:aggregate)&
!!$omp defaultmap(present:allocatable)
     do i=1,ncol
      do k=1,nlay+1
          swUpCS_slice(i,k) = 0.
          swDownCS_slice(i,k) = 0.
      enddo
     enddo
    end if

  end subroutine zero_SW_for_sunDown
  subroutine setupColumnProfiles(ibatch)

    use modglobal,   only: imax, jmax, kmax, i1, grav, kind_rb, rlv, cp, rd, pref0, tup, tdn
    use modfields,   only: thl0, qt0, ql0, exnf, rhof, sv0
    use modsurfdata, only: tskin, ps
    use modmicrodata, only : Nc_0,sig_g
    use modtracers, only: get_tracer_index

    implicit none

    integer, intent(in) :: ibatch
    integer :: jstart, jend
    integer :: i, j, k, icol
    integer :: inc
    real, parameter :: pi = 3.14159265358979
    real, parameter :: rho_liq = 1000., IWC0=50e-3 ! both in kg/m3

    real(kind=kind_rb) :: exners, reff_factor, ilratio, layerMass, qci, qcl, B_function
    real(kind_rb), allocatable :: nc_slice(:,:)

    allocate(nc_slice(ncol,nlay+1))
    !$acc enter data create(nc_slice)
!$omp target enter data map(alloc:nc_slice)

    exners = (ps/pref0)**(rd/cp)
    !reff_factor = 1e6*(3. /(4.*pi*Nc_0*rho_liq) )**(1./3.) * exp(log(sig_g)**2 )

    ! Set up j indices to be treated
    jstart = (ibatch-1) * jmax/nbatch + 2
    jend   =  ibatch    * jmax/nbatch + 1

    ! Set up layer values within the DALES domain
    !$acc parallel loop collapse(3) default(present) private(icol)
!!$omp target teams loop private(icol) collapse(3)&
!!$omp defaultmap(present:aggregate) defaultmap(present:allocatable)
    do k=1,kmax
      do j=jstart, jend
        do i=2,i1 !i1=imax+1
          icol=i-1+(j-jstart)*imax
          layerT(icol,k) = thl0(i,j,k) * exnf(k) + (rlv / cp) * ql0(i,j,k)
          h2ovmr(icol,k) = mwdry/mwh2o * max(qt0(i,j,k)-ql0(i,j,k),1e-10_field_r) !avoid negative values
        enddo
      enddo
    enddo

    call stop_on_err(gas_concs%set_vmr(trim(gas_names(1)), h2ovmr))

    ! Set up temperature interface values
    !$acc parallel loop collapse(3) default(present) private(icol)
!!$omp target teams loop private(icol) collapse(3)&
!!$omp defaultmap(present:aggregate) defaultmap(present:allocatable)
    do k=2,nlay
      do j=jstart, jend
        do i=2,i1 !i1=imax+1
          icol=i-1+(j-jstart)*imax
          interfaceT(icol,k) = (layerT(icol,k-1) + layerT(icol, k))/2.
        enddo
      enddo
    enddo
    !$acc parallel loop collapse(2) default(present) private(icol)
!!$omp target teams loop private(icol) collapse(2)&
!!$omp defaultmap(present:aggregate) defaultmap(present:allocatable)
    do j=jstart, jend
      do i=2,i1 !i1=imax+1
        icol=i-1+(j-jstart)*imax
        tg_slice(icol) = tskin(i,j)*exners
        interfaceT(icol, 1) = tg_slice(icol) !enforce ground temperature
        interfaceT(icol, nlay+1) = 2.*layerT(icol, nlay) - interfaceT(icol, nlay)
      enddo
    enddo

    ! Setup cloud properties (above the DALES domain everyhting is set to zero)
    !$acc kernels default(present)
!!$omp target defaultmap(present:aggregate)&
!!$omp defaultmap(present:allocatable) defaultmap(tofrom:scalar)
    LWP_slice = 0.0
    IWP_slice = 0.0
    liquidRe = 0.
    iceRe = 0.
    nc_slice = 0.0
    !$acc end kernels
!!$omp end target

    inc = get_tracer_index('Nc')

    if (inc > 0) then
       !$acc parallel loop collapse(3) default(present) private(icol)
!!$omp target teams loop private(icol) collapse(3)&
!!$omp defaultmap(present:aggregate) defaultmap(present:allocatable)
       do k=1,kmax
          do j=jstart, jend
             do i=2,i1
                icol=i-1+(j-jstart)*imax
                nc_slice(icol,k) = sv0(i,j,k,iNc)
             end do
          end do
       end do
    else
       !$acc parallel loop collapse(3) default(present) private(icol)
!!$omp target teams loop private(icol) collapse(3)&
!!$omp defaultmap(present:aggregate) defaultmap(present:allocatable)
       do k=1,kmax
          do j=jstart, jend
             do i=2,i1
                icol=i-1+(j-jstart)*imax
                nc_slice(icol,k) = Nc_0
             end do
          end do
       end do
    endif

    !$acc parallel loop collapse(3) default(present) private(icol,ilratio,layerMass,qcl,qci,B_function)
!!$omp target teams loop private(icol,ilratio,layermass,qcl,qci,&
!!$omp b_function) collapse(3) defaultmap(present:aggregate)&
!!$omp defaultmap(present:allocatable)
    do k=1,kmax
      do j=jstart, jend
        do i=2,i1 !i1=imax+1
          icol=i-1+(j-jstart)*imax
          ! set up working variables
          ilratio  = max(0.,min(1.,(layerT(icol,k)-tdn)/(tup-tdn)))! cloud water vs cloud ice partitioning
          layerMass = (interfaceP(icol,k)-interfaceP(icol,k+1))/grav !kg/m2
          qcl = ql0(i,j,k) * ilratio
          qci = ql0(i,j,k) * (1.0-ilratio)

          LWP_slice(icol,k) = qcl * layerMass*1e3 !g/m2
          IWP_slice(icol,k) = qci * layerMass*1e3 !g/m2

          if (LWP_slice(icol,k).gt.0.) then
            liquidRe(icol, k) = 1.e6*( 3.*( 1.e-3*LWP_slice(icol,k)/layerMass ) &
                              /(4.*pi*nc_slice(icol,k)*rho_liq) )**(1./3.) * exp(log(sig_g)**2 )
            !cstep: equation above contains function of many constants, are now absorbed in reff_factor
            !liquidRe(icol, k) = reff_factor  * qcl**(1./3.)

            if(liquidRe(icol,k).lt.2.5) liquidRe(icol,k) = 2.5
            if(liquidRe(icol,k).gt.20.) liquidRe(icol,k) = 20.
          endif

          if (IWP_slice(icol,k).gt.0) then
             !cstep Ou Liou: tempC = layerT(icol,k)--tmelt
             !cstep Ou Liou  iceRe(icol,k) = 326.3 + 12.42 * tempC + 0.197 * tempC**2 + 0.0012 * tempC**3  !cstep : Ou Liou 1995
             B_function =  -2 + 0.001 *(273.-layerT(icol,k))**1.5 * log10(qci*rhof(k)/IWC0) !Eq. 14 Wyser 1998
             iceRe(icol,k) = 377.4 + 203.3 * B_function + 37.91 * B_function**2 + 2.3696 * B_function**3 !micrometer, Wyser 1998, Eq. 35

            ! Ice optical properties in RRTMGP LUTs are given as a function of effective diameter
            !
            iceRe(icol,k) = iceRe(icol,k) * 2

             if(iceRe(icol,k).lt.10.) iceRe(icol,k) = 10.
             if(iceRe(icol,k).gt.180.) iceRe(icol,k) = 180.

          endif

        enddo
      enddo
    enddo

    !$acc exit data delete(nc_slice)
!$omp target exit data map(delete:nc_slice)
    deallocate(nc_slice)

  end subroutine setupColumnProfiles

  subroutine getFluxProfiles(ibatch)

    use modglobal,   only: i1, k1, imax, jmax, kmax, cp, dzf
    use modfields,   only: exnf, rhof

    implicit none

    integer, intent(in) :: ibatch
    integer :: jstart, jend
    integer :: i,j,k,icol

    ! Set up j indices to be treated
    jstart = (ibatch-1) * jmax/nbatch + 2
    jend   =  ibatch    * jmax/nbatch + 1

    !$acc parallel loop collapse(3) default(present) private(icol)
!!$omp target teams loop private(icol) collapse(3)&
!!$omp defaultmap(present:aggregate) defaultmap(present:allocatable)
    do k=1,k1
      do j=jstart, jend
        do i=2,i1 !i1=imax+1
          icol=i-1+(j-jstart)*imax
          lwu(i,j,k) = lwUp_slice(icol,k)
          lwd(i,j,k) =-lwDown_slice(icol,k)
          swu(i,j,k) = swUp_slice(icol,k)
          swd(i,j,k) =-swDown_slice(icol,k)
          swdir(i,j,k) = -swDownDir_slice(icol,k)
          swdif(i,j,k) = -(swDown_slice(icol,k) - swDownDir_slice(icol,k))
        enddo
      enddo
    enddo
    !$acc parallel loop collapse(2) default(present) private(icol)
!!$omp target teams loop private(icol) collapse(2)&
!!$omp defaultmap(present:aggregate) defaultmap(present:allocatable)
    do j=jstart, jend
      do i=2,i1 !i1=imax+1
        icol=i-1+(j-jstart)*imax
        LW_up_TOA(i,j) = lwUp_slice(icol,nlay+1)
        LW_dn_TOA(i,j) =-lwDown_slice(icol,nlay+1)
        SW_up_TOA(i,j) = swUp_slice(icol,nlay+1)
        SW_dn_TOA(i,j) =-swDown_slice(icol,nlay+1)
      enddo
    enddo

    if(doclearsky) then
      !$acc parallel loop collapse(3) default(present) private(icol)
!!$omp target teams loop private(icol) collapse(3)&
!!$omp defaultmap(present:aggregate) defaultmap(present:allocatable)
      do k=1,k1
        do j=jstart, jend
          do i=2,i1 !i1=imax+1
            icol=i-1+(j-jstart)*imax
            lwuca(i,j,k) = lwUpCS_slice(icol,k)
            lwdca(i,j,k) = -lwDownCS_slice(icol,k)
            swuca(i,j,k) =  swUpCS_slice(icol,k)
            swdca(i,j,k) = -swDownCS_slice(icol,k)
          enddo
        enddo
      enddo
      !$acc parallel loop collapse(2) default(present) private(icol)
!!$omp target teams loop private(icol) collapse(2)&
!!$omp defaultmap(present:aggregate) defaultmap(present:allocatable)
      do j=jstart, jend
        do i=2,i1 !i1=imax+1
          icol=i-1+(j-jstart)*imax
          SW_up_ca_TOA(i,j) = swUpCS_slice(icol,k)
          SW_dn_ca_TOA(i,j) =-swDownCS_slice(icol,k)
          LW_up_ca_TOA(i,j) = lwUpCS_slice(icol,k)
          LW_dn_ca_TOA(i,j) =-lwDownCS_slice(icol,k)
        enddo
      enddo
    endif

    !$acc parallel loop collapse(3) default(present)
!!$omp target teams loop collapse(3) defaultmap(present:aggregate)&
!!$omp defaultmap(present:allocatable)
    do k=1,kmax
      do j=jstart, jend
        do i=2,i1
          thlprad(i,j,k) = thlprad(i,j,k)-(lwd(i,j,k+1)-lwd(i,j,k)+lwu(i,j,k+1)-lwu(i,j,k)&
                                         +swd(i,j,k+1)-swd(i,j,k)+swu(i,j,k+1)-swu(i,j,k)) &
                                          /(rhof(k)*cp*exnf(k)*dzf(k))
        end do
      end do
    end do


  end subroutine

  subroutine setupSW(sunUp, ibatch)

    use modglobal,   only : xday,xlat,xlon,xtime,rtimee,imax,jmax,i1
    use shr_orb_mod, only : shr_orb_decl
    use modsurfdata, only : albedo

    implicit none

    integer, intent(in) :: ibatch
    integer :: jstart, jend
    integer :: i, j, icol
    logical,intent(out) :: sunUp
    real                :: dayForSW

    ! Set up j indices to be treated. We need to set the albedo correctly.
    jstart = (ibatch-1) * jmax/nbatch + 2
    jend   =  ibatch    * jmax/nbatch + 1

    if(doseasons) then
      ! The diurnal cycle of insolation will vary
      ! according to time of year of the current day.
      dayForSW = xday + (xtime + rtimee/3600) / 24
    end if

    call shr_orb_decl(dayForSW) ! Saves some orbital values to modraddata
    solarZenithAngleCos(:) =  &
         zenith(xtime*3600 + rtimee, xday, xlat, xlon) ! Used function in modraddata
    !$acc update device(solarZenithAngleCos)
!!$omp target update to(solarzenithanglecos)

    sunUp = .false.
    ! if all values in solarZenithAngleCos are >= its smallest positive, non-zero element
    if (all(solarZenithAngleCos(:) >= tiny(solarZenithAngleCos))) then
      sunUp = .true.

      ! Constant albedo for now
      ! Albedos can be computed as a function of solarZenithAngleCos,
      ! so it makes sense to keep the init here
      !$acc parallel loop collapse(2) default(present) private(icol)
!!$omp target teams loop private(icol) collapse(2)&
!!$omp defaultmap(present:aggregate) defaultmap(present:allocatable)
      do j=jstart, jend
        do i=2,i1 !i1=imax+1
          icol=i-1+(j-jstart)*imax
          sfc_alb_dif(:,icol) = albedo(i,j)
          sfc_alb_dir(:,icol) = albedo(i,j)
        enddo
      enddo
      !

    end if

  end subroutine setupSW

  subroutine setupEmis(ibatch)
    use modglobal,   only : imax,jmax,i1
    use modsurfdata, only : emissivity

    implicit none

    integer, intent(in) :: ibatch
    integer :: jstart, jend
    integer :: i, j, icol
    ! Set up j indices to be treated. We need to set the emissivity correctly.
    jstart = (ibatch-1) * jmax/nbatch + 2
    jend   =  ibatch    * jmax/nbatch + 1
      !$acc parallel loop collapse(2) default(present) private(icol)
!!$omp target teams loop private(icol) collapse(2)&
!!$omp defaultmap(present:aggregate) defaultmap(present:allocatable)
      do j=jstart, jend
        do i=2,i1 !i1=imax+1
          icol=i-1+(j-jstart)*imax
          emis(:,icol) = emissivity(i,j)
        enddo
      enddo
  end subroutine setupEmis

end module modradrte_rrtmgp
