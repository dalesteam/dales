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
  ! RTE-RRTMGP modules
  use mo_optical_props,      only: ty_optical_props, &
                                   ty_optical_props_arry, ty_optical_props_1scl, ty_optical_props_2str
  use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
  use mo_cloud_optics,       only: ty_cloud_optics
  use mo_source_functions,   only: ty_source_func_lw
  use mo_fluxes,             only: ty_fluxes_broadband
  use mo_gas_concentrations, only: ty_gas_concs

  implicit none

  private
  ! RRTMGP variables
  type(ty_gas_concs)                        :: gas_concs
  type(ty_source_func_lw), save             :: sources_lw
  type(ty_gas_optics_rrtmgp)                :: k_dist_lw, k_dist_sw
  type(ty_cloud_optics)                     :: cloud_optics_lw, cloud_optics_sw
  class(ty_optical_props_arry), allocatable :: atmos_lw, atmos_sw, clouds_lw, clouds_sw
  type(ty_fluxes_broadband)                 :: fluxes_lw, fluxes_sw, fluxes_cs_lw, fluxes_cs_sw
  real(kind=kind_rb), dimension(:,:), allocatable :: inc_sw_flux, sfc_alb_dir, sfc_alb_dif
  !Specify gas names, the first five (h2o, o3, co2, ch4 and n2o) are mandatory as they are major absorbers
  integer, parameter                        :: ngas = 10
  character(len=5), dimension(ngas)         :: gas_names = ['h2o  ', 'o3   ', 'co2  ', 'ch4  ', 'n2o  ', 'o2   ', 'cfc11', 'cfc12', 'cfc22', 'ccl4 ']
  integer                                   :: nlay, nlev, ncol, nbatch, nbndlw, nbndsw, ngptsw
  logical                                   :: doclearsky = .false.

  public :: radrte_rrtmgp

contains

  subroutine stop_on_err(error_msg)
    use iso_fortran_env, only : error_unit
    implicit none

    character(len=*), intent(in) :: error_msg

    if(error_msg /= "") then
      write (error_unit,*) trim(error_msg)
      write (error_unit,*) "modradrte_rrtmgp stopped"
      error stop 1
    end if
  end subroutine stop_on_err

  subroutine radrte_rrtmgp
    use mo_rte_lw,             only: rte_lw
    use mo_rte_sw,             only: rte_sw
    use mo_load_coefficients,  only: load_and_init
    use mo_load_cloud_coefficients, &
                               only: load_cld_lutcoeff, load_cld_padecoeff
    ! DALES modules
    use modradrrtmg,           only: readSounding, readTraceProfs
    use modmpi,                only: myid
    use modglobal,             only: imax, jmax, kmax, k1 
    use modfields,             only: initial_presh, initial_presf
    implicit none

    logical                 :: top_at_1 = .false., sunUp = .false.
    integer                 :: npatch, ibatch, ierr(3)=0
    character(len=256)      :: k_dist_file_lw = "rrtmgp-data-lw-g128-210809.nc"
    character(len=256)      :: k_dist_file_sw = "rrtmgp-data-sw-g112-210809.nc"
    character(len=256)      :: cloud_optics_file_lw = "rrtmgp-cloud-optics-coeffs-lw.nc"
    character(len=256)      :: cloud_optics_file_sw = "rrtmgp-cloud-optics-coeffs-reordered-sw.nc"
    integer                 :: ilay, icol, k

    ! Reading sounding (patch above Dales domain), only once
    if(.not.isReadSounding) then
      call readSounding(initial_presh(k1)/100.,npatch_start,npatch_end)

      if(npatch_end.ne.npatch_start) then
        npatch = npatch_end - npatch_start + 1
      else
        if(myid==0) write(*,*) 'No sounding levels above the LES domain, check sounding input file'
        stop 'ERROR: No valid radiation sounding found (modradrte_rrtmgp.f90)'
      end if

      !old notations nlay-1=kradmax=nzrad, nlay=krad1, nlev=krad2
      nlay = kmax + npatch + 1
      nlev = nlay + 1 ! necessary?
      !the two indices below are necessary for the readTraceProfs routine
      krad1=nlay
      krad2=nlev

      !Set the number of batch, between 1 and jmax, depending on memory available
      nbatch = jmax
      !Check if jmax is a mutliple of nbatch
      if(mod(jmax,nbatch)/=0) stop 'ERROR: Wrong batch number specified in modradrte_rrtmgp.f90'
      ncol = imax*jmax/nbatch

      isReadSounding = .true.
    end if

    ! Allocating working variables, only once
    if(.not.isAllocated_RadInputsOutputs) then
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
        if(myid==0) write(*,*) 'Could not allocate input/output arrays in modradrte_rrtmgp'
        stop 'ERROR: Radiation variables could not be allocated in modradrte_rrtmgp.f90'
      else
        isAllocated_RadInputsOutputs = .true.
      end if
    end if       

    ! Reading Trace Profiles
    if(.not.isReadTraceProfiles) then
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
      call readTraceProfs

      if(myid==0) write(*,*) 'Trace gas profile have been read'
      isReadTraceProfiles = .true.
    end if

    ! Specific RRTMGP initialization
    if(.not.isInitializedRrtmg) then

      call stop_on_err(gas_concs%init(gas_names))
      !setup trace gases concentration once for all
      do k=1,nlay; tracevmr(:,k) = o3(k); enddo
      call stop_on_err(gas_concs%set_vmr(trim(gas_names(2)), tracevmr))
      do k=1,nlay; tracevmr(:,k) = co2(k); enddo
      call stop_on_err(gas_concs%set_vmr(trim(gas_names(3)), tracevmr))
      do k=1,nlay; tracevmr(:,k) = ch4(k); enddo
      call stop_on_err(gas_concs%set_vmr(trim(gas_names(4)), tracevmr))
      do k=1,nlay; tracevmr(:,k) = n2o(k); enddo
      call stop_on_err(gas_concs%set_vmr(trim(gas_names(5)), tracevmr))
      do k=1,nlay; tracevmr(:,k) = o2(k); enddo
      call stop_on_err(gas_concs%set_vmr(trim(gas_names(6)), tracevmr))
      do k=1,nlay; tracevmr(:,k) = cfc11(k); enddo
      call stop_on_err(gas_concs%set_vmr(trim(gas_names(7)), tracevmr))
      do k=1,nlay; tracevmr(:,k) = cfc12(k); enddo
      call stop_on_err(gas_concs%set_vmr(trim(gas_names(8)), tracevmr))
      do k=1,nlay; tracevmr(:,k) = cfc22(k); enddo
      call stop_on_err(gas_concs%set_vmr(trim(gas_names(9)), tracevmr))
      do k=1,nlay; tracevmr(:,k) = ccl4(k); enddo
      call stop_on_err(gas_concs%set_vmr(trim(gas_names(10)), tracevmr))
      deallocate(tracevmr)

      ! Longwave init
      if(rad_longw) then

        ! Load k distributions
        call load_and_init(k_dist_lw, k_dist_file_lw, gas_concs)
        if(.not. k_dist_lw%source_is_internal()) &
          stop "modradrte_rrtmgp: k-distribution file isn't LW"
        nbndlw = k_dist_lw%get_nband()

        ! Initialize gas optical properties
        allocate(ty_optical_props_1scl::atmos_lw)
        select type(atmos_lw)
          class is (ty_optical_props_1scl)
            call stop_on_err(atmos_lw%alloc_1scl(ncol, nlay, k_dist_lw))
        end select

        ! Load cloud property data
        call load_cld_lutcoeff (cloud_optics_lw, cloud_optics_file_lw)
        !call load_cld_padecoeff(cloud_optics_lw, cloud_optics_file_lw)
        call stop_on_err(cloud_optics_lw%set_ice_roughness(2))

        ! Initialize cloud optical properties
        allocate(ty_optical_props_1scl::clouds_lw)
        call stop_on_err(clouds_lw%init(k_dist_lw%get_band_lims_wavenumber()))
        select type(clouds_lw)
          class is (ty_optical_props_1scl)
            call stop_on_err(clouds_lw%alloc_1scl(ncol, nlay))
        end select

        ! Allocate source term and define emissivity
        call stop_on_err(sources_lw%alloc(ncol, nlay, k_dist_lw))
        allocate(emis(nbndlw,ncol))
        emis=0.95

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
        call load_and_init(k_dist_sw, k_dist_file_sw, gas_concs)
        if(k_dist_sw%source_is_internal()) &
          stop "modradrte_rrtmgp: k-distribution file isn't SW"
        nbndsw = k_dist_sw%get_nband()
        ngptsw = k_dist_sw%get_ngpt()

        ! Initialize gas optical properties
        allocate(ty_optical_props_2str::atmos_sw)
        select type(atmos_sw)
          class is (ty_optical_props_2str)
            call stop_on_err(atmos_sw%alloc_2str(ncol, nlay, k_dist_sw))
        end select

        ! Load cloud property data
        call load_cld_lutcoeff (cloud_optics_sw, cloud_optics_file_sw)
        !call load_cld_padecoeff(cloud_optics_sw, cloud_optics_file_sw)
        call stop_on_err(cloud_optics_sw%set_ice_roughness(2))

        ! Initialize cloud optical properties
        allocate(ty_optical_props_2str::clouds_sw)
        call stop_on_err(clouds_sw%init(k_dist_sw%get_band_lims_wavenumber()))
        select type(clouds_sw)
          class is (ty_optical_props_2str)
            call stop_on_err(clouds_sw%alloc_2str(ncol, nlay))
        end select

        ! Define boundary conditions
        allocate(inc_sw_flux(ncol,ngptsw))
        allocate(sfc_alb_dir(nbndsw,ncol), sfc_alb_dif(nbndsw,ncol))

        fluxes_sw%flux_up => swUp_slice(:,:)
        fluxes_sw%flux_dn => swDown_slice(:,:)
        fluxes_sw%flux_dn_dir => swDownDir_slice(:,:)
        if(doclearsky) then
          fluxes_cs_sw%flux_up => swUpCS_slice(:,:)
          fluxes_cs_sw%flux_dn => swDownCS_slice(:,:)
        endif
      endif

      isInitializedRrtmg = .true.

    end if

    do ibatch = 1, nbatch

      call setupColumnProfiles(ibatch)

      if(rad_longw) then
        ! Compute optical properties and source
        call stop_on_err(k_dist_lw%gas_optics(layerP, interfaceP, & ! p_lay, p_lev (in, Pa)
                                              layerT, tg_slice, & ! t_lay, t_sfc (in, K)
                                              gas_concs, & ! gas volume mixing ratios (in)
                                              atmos_lw, & ! Optical properties (inout)
                                              sources_lw, & ! Planck source (inout)
                                              tlev = interfaceT)) ! t_lev (optional input, K)

        ! Solve clear sky radiation transport if required
        if(doclearsky) then
          call stop_on_err(rte_lw(atmos_lw, & ! optical properties (in)
                                  top_at_1, & ! Is the top of the domain at index 1? (in)
                                  sources_lw, & ! source function (in)
                                  emis, & ! emissivity at surface (in)
                                  fluxes_cs_lw)) ! fluxes (W/m2, inout)
        endif

        ! Compute and add cloud properties
        call stop_on_err(cloud_optics_lw%cloud_optics(LWP_slice, & ! cloud liquid water path (in, g/m2)
                                                      IWP_slice, & ! cloud ice water path (in, g/m2)
                                                      liquidRe, & ! cloud liquid particle effective size (in, microns)
                                                      iceRe, & ! cloud ice particle effective radius (in, microns)
                                                      clouds_lw)) ! cloud optical properties lw (inout)
        call stop_on_err(clouds_lw%increment(atmos_lw))

        ! Solve radiation transport
        call stop_on_err(rte_lw(atmos_lw, & ! optical properties (in)
                                top_at_1, & ! Is the top of the domain at index 1? (in)
                                sources_lw, & ! source function (in)
                                emis, & ! emissivity at surface (in)
                                fluxes_lw)) ! fluxes (W/m2, inout)
      endif

      if(rad_shortw) then

        ! setup incoming flux and albedo as a function of the zenith angle
        call setupSW(sunUp)

        if(sunUp) then
          ! Compute optical properties and incoming shortwave flux
          call stop_on_err(k_dist_sw%gas_optics(layerP, interfaceP, & ! p_lay, p_lev (in, Pa)
                                                layerT, & ! t_lay (in, K)
                                                gas_concs, & ! gas volume mixing ratios (in)
                                                atmos_sw, & ! Optical properties (inout)
                                                inc_sw_flux)) ! Incoming shortwave flux (inout)

          ! Solve clear sky radiation transport if required
          if(doclearsky) then
            call stop_on_err(rte_sw(atmos_sw, & ! optical properties (in)
                                    top_at_1, & ! Is the top of the domain at index 1? (in)
                                    solarZenithAngleCos, & ! cosine of the solar zenith angle (in)
                                    inc_sw_flux, & ! solar incoming flux (in)
                                    sfc_alb_dir, sfc_alb_dif, & ! surface albedos, direct and diffuse (in)
                                    fluxes_cs_sw)) ! fluxes (inout, W/m2)
          endif

          ! Compute and add cloud properties
          call stop_on_err(cloud_optics_sw%cloud_optics(LWP_slice, & ! cloud liquid water path (in, g/m2)
                                                        IWP_slice, & ! cloud ice water path (in, g/m2)
                                                        liquidRe, & ! cloud liquid particle effective size (in, microns)
                                                        iceRe, & ! cloud ice particle effective radius (in, microns)
                                                        clouds_sw)) ! cloud optical properties sw (inout)
          call stop_on_err(clouds_sw%delta_scale())
          call stop_on_err(clouds_sw%increment(atmos_sw))

          ! Solve radiation transport
          call stop_on_err(rte_sw(atmos_sw, & ! optical properties (in)
                                  top_at_1, & ! Is the top of the domain at index 1? (in)
                                  solarZenithAngleCos, & ! cosine of the solar zenith angle (in)
                                  inc_sw_flux, & ! solar incoming flux (in)
                                  sfc_alb_dir, sfc_alb_dif, & ! surface albedos, direct and diffuse (in)
                                  fluxes_sw)) ! fluxes (inout, W/m2)

        endif

      endif

      call getFluxProfiles(ibatch)

    enddo

  end subroutine radrte_rrtmgp

  subroutine setupColumnProfiles(ibatch)

    use modglobal,   only: imax, jmax, kmax, i1, j1, grav, kind_rb, rlv, cp, rd, pref0, tup, tdn
    use modfields,   only: thl0, qt0, ql0, exnf, rhof
    use modsurfdata, only: tskin, ps
    use modmicrodata, only : Nc_0,sig_g

    implicit none

    integer, intent(in) :: ibatch
    integer :: jstart, jend
    integer :: i, j, k, icol
    real(SHR_KIND_R4), parameter :: pi = 3.14159265358979
    real, parameter :: rho_liq = 1000., IWC0=50e-3 ! both in kg/m3

    real(kind=kind_rb) :: exners, reff_factor, ilratio, layerMass, qci, qcl, B_function

    exners = (ps/pref0)**(rd/cp)
    reff_factor = 1e6*(3. /(4.*pi*Nc_0*rho_liq) )**(1./3.) * exp(log(sig_g)**2 )

    ! Set up j indices to be treated
    jstart = (ibatch-1) * jmax/nbatch + 2
    jend   =  ibatch    * jmax/nbatch + 1

    ! Set up layer values within the DALES domain
    do j=jstart, jend
      do i=2,i1 !i1=imax+1
        icol=i-1+(j-jstart)*imax
        tg_slice(icol) = tskin(i,j)*exners
        do k=1,kmax
          layerP(icol,k) = presf_input(k)
          layerT(icol,k) = thl0(i,j,k) * exnf(k) + (rlv / cp) * ql0(i,j,k)
          h2ovmr(icol,k) = mwdry/mwh2o * max(qt0(i,j,k)-ql0(i,j,k),1e-10) !avoid negative values
        enddo
      enddo
    enddo

    ! Set up layer values above the DALES domain
    do j=jstart, jend
      do i=2,i1 !i1=imax+1
        icol=i-1+(j-jstart)*imax
        do k=1,nlay-kmax-1
          layerP(icol,kmax+k) = presf_input(kmax+k)
          layerT(icol,kmax+k) = tsnd(npatch_start+k-1)
          h2ovmr(icol,kmax+k) = mwdry/mwh2o * qsnd(npatch_start+k-1)
        enddo
        h2ovmr(icol,nlay) = h2ovmr(icol,nlay-1)
        layerP(icol,nlay) = 0.5*presh_input(nlay)
        layerT(icol,nlay) = 2.*layerT(icol,nlay-1)-layerT(icol, nlay-2)
      enddo
    enddo

    call stop_on_err(gas_concs%set_vmr(trim(gas_names(1)), h2ovmr))

    ! Set up interface values // use table assignment?
    do j=jstart, jend
      do i=2,i1 !i1=imax+1
        icol=i-1+(j-jstart)*imax
        interfaceT(icol, 1) = tg_slice(icol) !enforce ground temperature
        interfaceP(icol, 1) = presh_input(1)
        do k=2,nlay
          interfaceT(icol,k) = (layerT(icol,k-1) + layerT(icol, k))/2.
          interfaceP(icol,k) = presh_input(k)
        enddo
        interfaceT(icol, nlay+1) = 2.*layerT(icol, nlay) - interfaceT(icol, nlay)
        interfaceP(icol, nlay+1) = min(1.e-4_kind_rb , 0.25*layerP(1,nlay))
      enddo
    enddo

    ! pressures are in SI unit, i.e. Pa
    layerP(:,:)=layerP(:,:)*100.00
    interfaceP(:,:)=interfaceP(:,:)*100.00

    ! Setup cloud properties (above the DALES domain everyhting is set to zero)
    LWP_slice = 0.0
    IWP_slice = 0.0
    liquidRe = 0.
    iceRe = 0.
    do j=jstart, jend
      do i=2,i1 !i1=imax+1
        icol=i-1+(j-jstart)*imax
        do k=1,kmax
          ! set up working variables
          ilratio  = max(0.,min(1.,(layerT(icol,k)-tdn)/(tup-tdn)))! cloud water vs cloud ice partitioning
          layerMass = (interfaceP(icol,k)-interfaceP(icol,k+1))/grav !kg/m2
          qcl = ql0(i,j,k) * ilratio
          qci = ql0(i,j,k) * (1.0-ilratio)

          LWP_slice(icol,k) = qcl * layerMass*1e3 !g/m2
          IWP_slice(icol,k) = qci * layerMass*1e3 !g/m2


          if (LWP_slice(icol,k).gt.0.) then
            !cstep liquidRe(icol, k) = 1.e6*( 3.*( 1.e-3*LWP_slice(icol,k)/layerMass ) &
            !cstep                  /(4.*pi*Nc_0*rho_liq) )**(1./3.) * exp(log(sig_g)**2 )
            !cstep: equation above contains function of many constants, are now absorbed in reff_factor
            liquidRe(icol, k) = reff_factor  * qcl**(1./3.)

            if(liquidRe(icol,k).lt.2.5) liquidRe(icol,k) = 2.5
            if(liquidRe(icol,k).gt.20.) liquidRe(icol,k) = 20.
          endif

          if (IWP_slice(icol,k).gt.0) then
             !cstep Ou Liou: tempC = layerT(icol,k)--tmelt
             !cstep Ou Liou  iceRe(icol,k) = 326.3 + 12.42 * tempC + 0.197 * tempC**2 + 0.0012 * tempC**3  !cstep : Ou Liou 1995
             B_function =  -2 + 0.001 *(273.-layerT(icol,k))**1.5 * log10(qci*rhof(k)/IWC0) !Eq. 14 Wyser 1998
             iceRe(icol,k) = 377.4 + 203.3 * B_function + 37.91 * B_function**2 + 2.3696 * B_function**3 !micrometer, Wyser 1998, Eq. 35
             if (isnan(iceRe(icol,k))) then
                write (*,*) "B", B_function, "iceRe", iceRe(icol,k), "qci", qci, "layerT", layerT(icol,k)
                stop "modradrte_rrtmgp: iceRe is nan."
             endif

             if(iceRe(icol,k).lt.10.) iceRe(icol,k) = 10.
             if(iceRe(icol,k).gt.180.) iceRe(icol,k) = 180.
          endif

        enddo
      enddo
    enddo

  end subroutine setupColumnProfiles

  subroutine getFluxProfiles(ibatch)

    use modglobal,   only: i1, j1, k1, imax, jmax, kmax, cp, dzf
    use modfields,   only: exnf, rhof

    implicit none

    integer, intent(in) :: ibatch
    integer :: jstart, jend
    integer :: i,j,k,icol

    ! Set up j indices to be treated
    jstart = (ibatch-1) * jmax/nbatch + 2
    jend   =  ibatch    * jmax/nbatch + 1

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

  subroutine setupSW(sunUp)

    use modglobal,   only : xday,xlat,xlon,imax,xtime,rtimee
    use shr_orb_mod, only : shr_orb_decl
    use modsurfdata, only : albedoav

    implicit none

    logical,intent(out) :: sunUp
    real                :: dayForSW

    if(doseasons) then
      ! The diurnal cycle of insolation will vary
      ! according to time of year of the current day.
      dayForSW = xday + (xtime + rtimee/3600) / 24
    end if

    call shr_orb_decl(dayForSW) ! Saves some orbital values to modraddata
    solarZenithAngleCos(:) =  &
         zenith(xtime*3600 + rtimee, xday, xlat, xlon) ! Used function in modraddata

    sunUp = .false.
    ! if all values in solarZenithAngleCos are >= its smallest positive, non-zero element
    if (all(solarZenithAngleCos(:) >= tiny(solarZenithAngleCos))) then
      sunUp = .true.

      ! Constant albedo for now
      ! Albedos can be computed as a function of solarZenithAngleCos,
      ! so it makes sense to keep the init here
      sfc_alb_dir=albedoav
      sfc_alb_dif=albedoav

    end if

  end subroutine setupSW

end module modradrte_rrtmgp
