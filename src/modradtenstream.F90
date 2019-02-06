module modradtenstream


#ifndef HAVE_TENSTREAM
  contains
    subroutine dales_tenstream
      stop 'Your build does not support the TenStream solver for 3D Radiative Transfer... please reconfigure with -DWITH_TENSTREAM=ON'
    end subroutine

#else
  use m_data_parameters, only: ireals, iintegers, mpiint, one, zero,  init_mpi_data_parameters, pi, default_str_len
  use m_dyn_atm_to_rrtmg, only: t_tenstr_atm, setup_tenstr_atm, destroy_tenstr_atm
  use m_pprts_rrtmg, only : pprts_rrtmg, destroy_pprts_rrtmg
  use m_pprts_base, only : t_solver_3_10,t_solver_8_16
  use m_helper_functions, only : imp_allgather_int_inplace,rad2deg,deg2rad,CHKERR,reorder_mpi_comm
  use modfields, only : &
    qt0, &    ! (:,:,:) total specific humidity at time step t
    ql0, &    ! (:,:,:) liquid water content
    tmp0, &   ! (:,:,:) temperature at full level
    thl0, &   ! (:,:,:) liq. water pot. temperature at time step t
    thv0h, &  ! (:,:,:) theta_v at half level
    rhof, &   ! (:)     slab averaged density at full level
    exnf, &   ! (:)     hydrostatic exner function at full level
    presf, &  ! (:)     hydrostatic pressure at full level
    presh     ! (:)     hydrostatic pressure at half level

  use modglobal, only : imax,jmax,kmax,i1,j1,k1, &
      rlv, cp, dzf

  use modraddata, only : &
      mwdry, mwh2o, & !< molecular weights
      lCnstAlbedo,  & !< use of a constant albedo
      zenith , & !<   computes the cosine of the solar zenith angle
      thlprad,thlprSW,thlprLW, & !<   the radiative tendencies
      swd    , & !<   shortwave downward radiative flux
      swdir  , & !<   Direct shortwave downward radiative flux
      swdif  , & !<   Difuse shortwave downward radiative flux
      swu    , & !<   shortwave upward radiative flux
      lwd    , & !<   longwave downward radiative flux
      lwu    , & !<   longwave upward radiative flux
      SW_up_TOA, SW_dn_TOA, LW_up_TOA, LW_dn_TOA !< Top of the atmosphere radiative fluxes

  implicit none
  save

  private
  public :: dales_tenstream

  logical,parameter :: ldebug=.False.
!  logical,parameter :: ldebug=.True.
  type(t_solver_3_10) :: pprts_solver
 

contains

  subroutine dales_tenstream
    use modmpi, only : comm3d, myid, nprocx, nprocy
    use modglobal, only : dx, dy, xday, xlat, xlon,cp,Rd, xtime, rtimee,pref0
    use modmicrodata, only : Nc_0,sig_g
    use mpi, only : mpi_barrier
    use modsurfdata, only : albedoav,tskin,ps

    character(len=default_str_len),parameter :: atm_filename='afglus_100m.dat'
    real(ireals),allocatable, dimension(:,:,:) :: edir,edn,eup,abso ! [nlev_merged(-1), nxp, nyp]

    real(ireals) :: phi0=180, theta0
    real(ireals) :: albedo_thermal, albedo_solar

    real(ireals),dimension(1:k1,  2:i1, 2:j1),target :: d_plev, d_tlev
    real(ireals),dimension(1:kmax,2:i1, 2:j1),target :: d_tlay, d_h2ovmr
    real(ireals),dimension(1:kmax,2:i1, 2:j1),target :: d_lwc, d_reliq ! [g/kg], [micron]
    real(ireals), pointer, dimension(:,:) :: pplev, ptlev, ptlay,plwc,ph2ovmr, preliq

    real(ireals) :: mu, modeltime
    integer(mpiint) :: inp_comm,mpierr
    integer(iintegers) :: i, j, k, kk
    integer(iintegers), allocatable :: nxproc(:), nyproc(:)

    real(ireals), parameter :: solar_min_sza=85 ! minimum solar zenith angle -- below, dont compute solar rad
    real(ireals), parameter :: rho_liq = 1000
  
    type(t_tenstr_atm) :: atm
  
    mu = real(zenith(xtime*3600 + rtimee, xday, xlat, xlon), ireals)
    theta0 = rad2deg(acos(mu))
  
    albedo_thermal = 0
    if (lCnstAlbedo) then
      albedo_solar = real(albedoav, ireals)
    else
      call CHKERR(1_mpiint, 'Tenstream currently only supports lCnstAlbedo=.True.')
    end if

    if(ldebug .and. myid.eq.0) then
      print *,'Calling DALES TenStream Wrapper with:'
      print *,'dx/dy', dx,dy
      print *,'Domain x1', 2,k1,':', 2,i1,':', 2,j1
      print *,'Domain xmax', 2,kmax,':', 2,imax,':', 2,jmax
      print *,'presh', presh,' : (', shape(presh),')'
      print *,'presf', presf,' : (', shape(presf),')'
      print *,'mu',xtime*3600, rtimee, xday, xlat, xlon, '::', mu, '::', theta0, 'deg'
      print *,'allocated sw flxs',allocated(swdir), allocated(swdif), allocated(swd), allocated(swu)
      print *,'shape sw flxs',shape(swdir), shape(swdif), shape(swd), shape(swu)
      print *,'allocated lw flxs',allocated(lwd), allocated(lwu)
      print *,'shape lw flxs', shape(lwd), shape(lwu)
      print *,'allocated thlprad', allocated(thlprad), shape(thlprad)
    endif

    !inp_comm = comm3d
    call reorder_mpi_comm(comm3d, nprocx,nprocy, inp_comm)
    allocate(nxproc(nprocx), source=imax)
    allocate(nyproc(nprocy), source=jmax)
    
    modeltime = real(rtimee, ireals)

    do j=2,j1
      do i=2,i1
        d_plev(:, i, j) = real(presf(:), ireals)/100
        do k=1,kmax
          d_tlay(k,i,j) = real(thl0(i,j,k) * exnf(k) + (rlv / cp) * ql0(i,j,k), ireals)
          d_h2ovmr(k,i,j) = real(mwdry/mwh2o * max((qt0(i,j,k) - ql0(i,j,k)),1e-18), ireals)

          d_lwc(k,i,j) = real(ql0(i,j,k) * 1e3, ireals)

          d_reliq(k,i,j) = real(1.e6*( 3.*( 1e-3 * d_lwc(k,i,j) ) &
                            /(4.*pi*Nc_0*rho_liq) )**(1./3.) * exp(log(sig_g)**2 ), ireals)
          d_reliq(k,i,j) = min(max(d_reliq(k,i,j), 2.5_ireals), 60._ireals)

        enddo
      enddo
    enddo
    do k=2,kmax
      d_tlev(k,:,:) = real(d_tlay(k-1,:,:) + d_tlay(k,:,:), ireals)*.5_ireals
    enddo
    d_tlev(k1,:,:) = 2*d_tlay(kmax,:,:) - d_tlev(kmax,:,:) ! linear extrapolation towards top
    !d_tlev(1,:,:) = 2*d_tlay(1,:,:) - d_tlev(2,:,:) ! same for bottom layer TODO: this needs something better, at least when we have an interactive surface
    d_tlev(1,:,:) = tskin(2:i1,2:j1) * (ps/pref0) ** (rd/cp) !use skin temperature from surface scheme
    if(ldebug .and. myid.eq.0) then
      print *,'shape for d_plev', shape(d_plev)
      print *,'shape for d_tlay', shape(d_tlay)
      print *,'d_tlay(:,2,2)', d_tlay(:,2,2)
      print *,'d_tlev(:,2,2)', d_tlev(:,2,2)
      print *,'d_lwc(:,2,2)', d_lwc(:,2,2)
      print *,'d_reliq(:,2,2)', d_reliq(:,2,2)
    endif

    call init_mpi_data_parameters(inp_comm)


    pplev(1:size(d_plev,1),1:size(d_plev,2)*size(d_plev,3)) => d_plev
    ptlev(1:size(d_tlev,1),1:size(d_tlev,2)*size(d_tlev,3)) => d_tlev
    ptlay(1:size(d_tlay,1),1:size(d_tlay,2)*size(d_tlay,3)) => d_tlay
    ph2ovmr(1:size(d_h2ovmr,1),1:size(d_h2ovmr,2)*size(d_h2ovmr,3)) => d_h2ovmr   
    plwc (1:size(d_lwc ,1),1:size(d_lwc ,2)*size(d_lwc ,3)) => d_lwc
    preliq(1:size(d_reliq,1),1:size(d_reliq,2)*size(d_reliq,3)) => d_reliq

    call setup_tenstr_atm(inp_comm,.False.,atm_filename, &
      pplev,ptlev,atm,ptlay,d_h2ovmr=ph2ovmr,d_lwc=plwc,d_reliq=preliq)





    ! Thermal RT
    call pprts_rrtmg(inp_comm, pprts_solver,atm,imax,jmax, &
            real(dx, ireals), real(dy, ireals),           &
            phi0, theta0,                                 &
            albedo_thermal, albedo_solar,   &
            .true., .false.,                              &
            edir,edn,eup,abso,                            &
            !d_plev=d_plev, d_tlev=d_tlev,                 &
           ! d_h2ovmr=d_h2ovmr,                            &
            !d_lwc=d_lwc, d_reliq=d_reliq,                 &
            nxproc=nxproc, nyproc=nyproc, opt_time=modeltime)

    if(ldebug .and. myid.eq.0) then
      print *,'shape abso', shape(abso)
      do k=lbound(abso,1), ubound(abso,1)
        print *,k,'ediff',edn(k,1,1), eup(k,1,1),'::',abso(k,1,1)
      enddo
      k = ubound(edn,1)
      print *,k,'ediff',edn(k,1,1), eup(k,1,1)
    endif

    do k=1,k1
      kk = ubound(edn,1) - (k-1) ! TenStream providing fluxes from TOA to Srfc and on the bigger grid...
      lwd(2:i1, 2:j1, k) = -edn(kk, :, :)
      lwu(2:i1, 2:j1, k) = eup(kk, :, :)
    enddo
    LW_dn_TOA(2:i1, 2:j1) = -edn(1, :, :)
    LW_up_TOA(2:i1, 2:j1) = eup(1, :, :)

    do k=kmax,1,-1
      kk = ubound(abso,1) - (k-1) ! TenStream providing fluxes from TOA to Srfc and on the bigger grid...
      thlprad(2:i1, 2:j1, k) = thlprad(2:i1, 2:j1, k) + abso(kk, :, :) / (rhof(k)*cp*exnf(k))
      thlprLW(2:i1, 2:j1, k) = thlprLW(2:i1, 2:j1, k) + abso(kk, :, :) / (rhof(k)*cp*exnf(k))
      if(myid.eq.0 .and. ldebug) then
        print *,'thermal thlprad before',k,'::',kk,'::',thlprad(2, 2, k),'::', abso(kk,2,2) / (rhof(k)*cp*exnf(k)) *3600*24, 'K/day'
      endif
    enddo


    if(mu.gt.cos(deg2rad(solar_min_sza))) then
      ! Solar RT
      call pprts_rrtmg(inp_comm, pprts_solver,atm,imax,jmax,&
              real(dx,ireals), real(dy,ireals),             &
              phi0, theta0,                                 &
              albedo_thermal, albedo_solar,   &
              .false., .true.,                              &
              edir,edn,eup,abso,                            &
          !    d_plev=d_plev, d_tlev=d_tlev,                 &
          !    d_h2ovmr=d_h2ovmr,                            &
          !    d_lwc=d_lwc, d_reliq=d_reliq,                 &
              nxproc=nxproc, nyproc=nyproc, opt_time=modeltime)

      if(ldebug .and. myid.eq.0) then
        print *,'shape abso', shape(abso)
        do k=lbound(abso,1), ubound(abso,1)
          print *,k,'edir',edir(k,1,1),'ediff',edn(k,1,1), eup(k,1,1),'::',abso(k,1,1)
        enddo
        k = ubound(edir,1)
        print *,k,'edir',edir(k,1,1),'ediff',edn(k,1,1), eup(k,1,1)
      endif
      
      do k=1,k1
        kk = ubound(edn,1) - (k-1) ! TenStream providing fluxes from TOA to Srfc and on the bigger grid...
        swdir(2:i1, 2:j1, k) = edir(kk, :, :)
        swdif(2:i1, 2:j1, k) = -edn(kk, :, :)
        swd(2:i1, 2:j1, k) = -(edir(kk, :, :) + edn(kk, :, :))
        swu(2:i1, 2:j1, k) = eup(kk, :, :)
      enddo
      SW_dn_TOA(2:i1, 2:j1) = -(edir(1, :, :) + edn(1, :, :))
      SW_up_TOA(2:i1, 2:j1) = eup(1, :, :)

      do k=kmax,1,-1
        kk = ubound(abso,1) - (k-1) ! TenStream providing fluxes from TOA to Srfc and on the bigger grid...
        thlprad(2:i1, 2:j1, k) = thlprad(2:i1, 2:j1, k) + abso(kk, :, :) / (rhof(k)*cp*exnf(k))
        thlprSW(2:i1, 2:j1, k) = thlprSW(2:i1, 2:j1, k) + abso(kk, :, :) / (rhof(k)*cp*exnf(k))
       if(myid.eq.0 .and. ldebug) then
          print *,'solar thlprad before',k,'::',kk,'::',thlprad(2, 2, k)*3600*24,'::', &
            abso(kk,2,2), &
            (rhof(k)*cp*exnf(k)), &
            abso(kk,2,2) / (rhof(k)*cp*exnf(k)) *3600*24, 'K/day'
        endif
      enddo

    else ! sun is not up
      swdir(2:i1, 2:j1, :) = 0
      swdif(2:i1, 2:j1, :) = 0
      swd(2:i1, 2:j1, :) = 0
      swu(2:i1, 2:j1, :) = 0
      SW_dn_TOA(2:i1, 2:j1) = 0
      SW_up_TOA(2:i1, 2:j1) = 0
    endif

  end subroutine
#endif
end module modradtenstream
