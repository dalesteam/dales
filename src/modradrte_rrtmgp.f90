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
  implicit none

  private
  integer :: nlay, nlev, ncol, nbnblw
  public :: radrte_rrtmgp

contains

  subroutine stop_on_err(error_msg)
    use iso_fortran_env, only : error_unit
    character(len=*), intent(in) :: error_msg

    if(error_msg /= "") then
      write (error_unit,*) trim(error_msg)
      write (error_unit,*) "modradrte_rrtmgp stopped"
      error stop 1
    end if
  end subroutine stop_on_err

  subroutine radrte_rrtmgp
    ! RTE-RRTMGP modules
    use mo_optical_props,      only: ty_optical_props, &
                                     ty_optical_props_arry, ty_optical_props_1scl, ty_optical_props_2str
    use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
    !use mo_cloud_optics,       only: ty_cloud_optics
    use mo_gas_concentrations, only: ty_gas_concs
    use mo_source_functions,   only: ty_source_func_lw
    use mo_fluxes,             only: ty_fluxes_broadband
    use mo_rte_lw,             only: rte_lw
    use mo_rte_sw,             only: rte_sw
    use mo_load_coefficients,  only: load_and_init  
    ! DALES modules
    use modradrrtmg,           only: readSounding
    use modmpi,                only: myid
    use modglobal,             only: imax, jmax, kmax, k1 
    use modfields,             only: initial_presh, initial_presf

    type(ty_source_func_lw), save             :: sources_lw
    type(ty_gas_optics_rrtmgp)                :: k_dist_lw!, k_dist_sw
    !type(ty_cloud_optics)                     :: cloud_optics
    type(ty_gas_concs)                        :: gas_concs
    class(ty_optical_props_arry), allocatable :: atmos_lw!, clouds_lw
    type(ty_fluxes_broadband)                 :: fluxes_lw

    logical                                   :: top_at_1 = .false.
    integer                                   :: nbndlw, npatch, ierr(2)=0
    integer, parameter                        :: ngas = 1
    character(len=256)                        :: k_dist_file_lw, k_dist_file_sw
    character(len=3), dimension(ngas)         :: gas_names = ['h2o']!, 'co2', 'o3 ', 'n2o', 'co ', 'ch4', 'o2 ', 'n2 ']

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
      ncol = imax*jmax

      isReadSounding = .true.
    end if

    ! Allocating working variables, only once
    if(.not.isAllocated_RadInputsOutputs) then
      allocate(layerP(ncol,nlay), &
               layerT(ncol,nlay), &
               h2ovmr(ncol,nlay), &
               tg_slice(ncol), &
               presf_input(nlay-1), &
               STAT=ierr(1))
      allocate(interfaceP(ncol,nlay+1), &
               interfaceT(ncol,nlay+1), &
               lwUp_slice(ncol,nlay+1), &
               lwDown_slice(ncol,nlay+1), &
               presh_input(nlay), &
               STAT=ierr(2))
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
      !call readTraceProfs ! Disable it for now, only H2O!!!

      if(myid==0) write(*,*) 'Trace gas profile have been read'
      isReadTraceProfiles = .true.
    end if

    ! Specific RRTMGP initialization
    if(.not.isInitializedRrtmg) then
      call stop_on_err(gas_concs%init(gas_names))

      if (rad_longw) then

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
          class default
            call stop_on_err("modradrte_rrtmgp.f90: Don't recognize the kind of optical properties ")
        end select
 
        ! Initialize cloud optical properties
        !allocate(ty_optical_props_1scl::clouds)

        ! Allocate source term and define BC
        call stop_on_err(sources_lw%alloc(ncol, nlay, k_dist_lw))
        allocate(emis(nbnblw,ncol))
        emis=0.95

      endif
      if (rad_shortw) then
        ! Some sw init 
      endif
      isInitializedRrtmg = .true.
    end if

    call setupColumnProfiles()
    call stop_on_err(gas_concs%set_vmr(trim(gas_names(1)), h2ovmr))

    if(rad_longw) then
 
      fluxes_lw%flux_up => lwUp_slice(:,:)
      fluxes_lw%flux_dn => lwDown_slice(:,:)

      ! Compute optical properties and source
      call stop_on_err(k_dist_lw%gas_optics(layerP, interfaceP, & ! p_lay, p_lev (in)
                                            layerT, tg_slice, & ! t_lay, t_sfc (in)
                                            gas_concs, & ! gas volume mixing ratios
                                            atmos_lw, & ! Optical properties (inout)
                                            sources_lw, & ! Planck source (inout)
                                            tlev = interfaceT)) !t_lev (optional input)

      ! Add cloud properties
      !call stop_on_err(clouds_lw%increment(atmos_lw))

      ! Solve radiation transport
      call stop_on_err(rte_lw(atmos_lw, & ! optical properties
                              top_at_1, & ! Is the top of the domain at index 1?
                              sources_lw, & ! source function
                              emis, & ! emissivity at surface
                              fluxes_lw)) ! fluxes
    endif

    call getFluxProfiles()

  end subroutine radrte_rrtmgp

  subroutine setupColumnProfiles

    use modglobal,   only: imax, jmax, kmax, i1, j1, kind_rb, rlv, cp, rd, pref0
    use modfields,   only: thl0, qt0, ql0, exnf
    use modsurfdata, only: tskin, ps

    implicit none

    integer :: i, j, k, icol
    real(kind=kind_rb) :: exners

    exners = (ps/pref0)**(rd/cp)

    ! Set up layer values within the DALES domain
    do i=2,i1 !i1=imax+1
      do j=2,j1 !j1=jmax+1
        icol=j-1+(i-2)*jmax
        tg_slice(icol) = tskin(i,j)*exners
        do k=1,kmax
          layerP(icol,k) = presf_input(k)
          layerT(icol,k) = thl0(i,j,k) * exnf(k) + (rlv / cp) * ql0(i,j,k)
          h2ovmr(icol,k) = mwdry/mwh2o * max(qt0(i,j,k)-ql0(i,j,k),1e-10) !avoid negative values
        enddo
      enddo
    enddo

    ! Set up layer values above the DALES domain
    do i=2,i1 !i1=imax+1
      do j=2,j1 !j1=jmax+1
        icol=j-1+(i-2)*jmax
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

    ! Set up interface values // use table assignment?
    do i=2,i1 !i1=imax+1
      do j=2,j1 !j1=jmax+1
        icol=j-1+(i-2)*jmax
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

  end subroutine setupColumnProfiles

  subroutine getFluxProfiles

    use modglobal,   only: i1, j1, k1, jmax, kmax, cp, dzf
    use modfields,   only: exnf, rhof

    implicit none

    integer :: i,j,k,icol

    do i=2,i1
      do j=2,j1
        icol=j-1+(i-2)*jmax
        do k=1,k1
          lwu(i,j,k) = lwUp_slice(icol,k)
          lwd(i,j,k) =-lwDown_slice(icol,k)
        enddo
        LW_up_TOA(i,j) =  lwUp_slice(icol,nlay+1)
        LW_dn_TOA(i,j) = -lwDown_slice(icol,nlay+1)
      enddo
    enddo

    do k=1,kmax
      do j=2,j1
        do i=2,i1
          thlprad(i,j,k) = thlprad(i,j,k)-(lwd(i,j,k+1)-lwd(i,j,k)+lwu(i,j,k+1)-lwu(i,j,k))&
                                         !+swd(i,j,k+1)-swd(i,j,k)+swu(i,j,k+1)-swu(i,j,k)) &
                                          /(rhof(k)*cp*exnf(k)*dzf(k))
        end do
      end do
    end do


  end subroutine

end module modradrte_rrtmgp
