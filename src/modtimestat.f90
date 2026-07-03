!> \file modtimestat.f90
!!  Timestat calculates timeseries of several variables

!>
!! Timestat calculates timeseries of several variables
!>
!! Timeseries of the most relevant parameters. Written to tmser1.expnr and tmsurf.expnr
!! If netcdf is true, this module leads the tmser.expnr.nc output
!!  \author Pier Siebesma, K.N.M.I.
!!  \author Stephan de Roode, TU Delft
!!  \author Chiel van Heerwaarden, Wageningen U.R.
!!  \author Thijs Heus, MPI-M
!!  \author Fredrik Jansson, TU Delft
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


module modtimestat
  use modtimer
  use modprecision, only : longint, field_r
  use modlogging, only: finish

implicit none
character(len=*), parameter :: modname = 'modtimestat'
! private
! PUBLIC :: inittimestat, timestat
save
!NetCDF variables

  integer :: nvar
  integer :: ncid,nrec = 0
  integer :: ivar_rad = 0  ! starting index for radiation quantities
  character(80) :: fname = 'tmser.xxx.nc'
  !character(80),dimension(nvar,4) :: ncname
  character(80), allocatable, dimension(:,:)    :: ncname
  character(40) :: name

  real    :: dtav
  integer(kind=longint) :: idtav,tnext
  logical :: ltimestat= .false. !<switch for timestatistics (on/off)
  real    :: zi,ziold=-1, we
  integer, parameter :: iblh_flux = 1, iblh_grad = 2, iblh_thres = 3
  integer, parameter :: iblh_thv = -1, iblh_thl = -2, iblh_qt = -3
  integer :: iblh_meth = iblh_grad, iblh_var = iblh_thv
  integer :: blh_nsamp = 4
  real    :: blh_thres=-1 ,blh_sign=1.0
  real(field_r) :: zbaseav, ztopav, ztopmax,zbasemin
  real(field_r) :: qlintav, qlintmax, tke_tot
  real(field_r) :: prav, pravl
  real(field_r) :: qtintav, qrintav
  real(field_r) :: cc, wmax, qlmax
  real(field_r) :: qlint, qtint, qrint
  logical:: store_zi = .false.
  real(field_r), allocatable, dimension(:) :: profile, gradient, dgrad
  real(field_r), allocatable, dimension(:,:,:) :: blh_fld
  real(field_r), allocatable,dimension(:,:,:) :: sv0h

contains
!> Initializing Timestat. Read out the namelist, initializing the variables
  subroutine inittimestat
    use modmpi,    only : myid,comm3d,mpierr,D_MPI_BCAST
    use modglobal, only : ifnamopt, fname_options,cexpnr,dtmax,ifoutput,dtav_glob,tres,&
                          ladaptive,k1,kmax,rd,rv,dt_lim,btime,i1,j1,lwarmstart,checknamelisterror, &
                          ih ,jh
    use modfields, only : thlprof,qtprof,svprof
    use modsurfdata, only : isurf
    use modstat_nc, only : lnetcdf, open_nc, define_nc, ncinfo, nctiminfo
    use modraddata, only : iradiation
    use modlsm, only : lags
    use fortran_support, only: nnml_output
    implicit none

    character(len=*), parameter :: routine = modname//'/inittimestat'

    integer :: ierr,k,location = 1
    integer :: i,j

    namelist/NAMTIMESTAT/ & !< namelist
    dtav,ltimestat,blh_thres,iblh_meth,iblh_var,blh_nsamp !! namelist contents


    dtav=dtav_glob
    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMTIMESTAT,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMTIMESTAT')
      write(nnml_output ,NAMTIMESTAT)
      close(ifnamopt)
    end if

    call D_MPI_BCAST(dtav     ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(ltimestat  ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(blh_thres,1,0,comm3d,mpierr)
    call D_MPI_BCAST(iblh_meth  ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(iblh_var   ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(blh_nsamp  ,1,0,comm3d,mpierr)
    idtav = int(dtav / tres, kind=kind(idtav))

    tnext = idtav+btime

    nvar = 24
    ivar_rad = 25
    if(isurf == 1) then
       nvar = nvar + 11
       ivar_rad = ivar_rad + 11
    else if (isurf  == 11) then
      nvar = nvar + 7
      ivar_rad = ivar_rad + 7
      if (lags) then
        nvar = nvar + 2
        ivar_rad = ivar_rad + 2
      end if
    end if
    if (iradiation /= 0) then
       nvar = nvar + 19
    end if

    if(.not.(ltimestat)) return

    call timer_tic('modtimestat/inittimestat', 0)

    dt_lim = min(dt_lim,tnext)

    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      call finish(routine, 'TIMESTAT: dtav should be a integer multiple of dtmax')
    end if

    allocate(blh_fld(2-ih:i1+ih,2-jh:j1+jh,k1),sv0h(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(profile(k1),gradient(k1),dgrad(k1))

    gradient = 0.
    profile = 0.

    select case (iblh_var)
    case(iblh_qt)
      profile = qtprof
    case(iblh_thl)
      profile = thlprof
    case(iblh_thv)
      do k=1,k1
        profile(k) = thlprof(k)*(1+(rv/rd-1)*qtprof(k))
      end do
    case(1:)
      profile = svprof(:,iblh_var)
    end select
    blh_sign = sign(1.0_field_r,profile(kmax)-profile(1))

    select case(iblh_meth)
    case (iblh_flux)
    case (iblh_grad)
    case (iblh_thres)
      if (blh_thres<0) then
        do k=kmax,2,-1
          if (blh_sign*(profile(k+1) - profile(k-1)) > gradient(k)) then
            location = k
            gradient = blh_sign*(profile(k+1) - profile(k-1))
          endif
        enddo
        blh_thres=profile(location)
        if (myid==0) write (*,*) 'TIMESTAT: blh_tres =',blh_thres
      end if
    case default
      call finish(routine, 'TIMESTAT: Incorrect iblh_meth')
    end select

    if(myid==0) then
       if (.not. lwarmstart) then
          !tmser1
          open (ifoutput,file='tmser1.'//cexpnr,status='replace',position='append')
          write(ifoutput,'(2a)') &
               '#  time      cc     z_cbase    z_ctop_avg  z_ctop_max      zi         we', &
               '   <<ql>>  <<ql>>_max   w_max   tke     ql_max'
          close(ifoutput)
          !tmsurf
          open (ifoutput,file='tmsurf.'//cexpnr,status='replace',position='append')
          write(ifoutput,'(2a)') &
               '#  time        ust        tst        qst         obukh', &
               '      thls        z0        wthls      wthvs      wqls '
          close(ifoutput)
          if(isurf == 1) then
             open (ifoutput,file='tmlsm.'//cexpnr,status='replace',position='append')
             write(ifoutput,'(4a)') &
                  '#     time      Qnet        H          LE         G0  ', &
                  '   tendskin     rs         ra        tskin        cliq  ', &
                  '    Wl          rssoil     rsveg       Resp       wco2         An', &
                  '    gcco2'
             write(ifoutput,'(4a)') &
                  '#      [s]     [W/m2]     [W/m2]     [W/m2]     [W/m2]', &
                  '   [W/m2]      [s/m]       [s/m]     [K]          [-]   ', &
                  '   [m]          [s/m]      [s/m]   [mgCm2/s]               [mgCm2/s]',&
                  '   [m/s]  '
             close(ifoutput)
          end if

       endif

       if (lnetcdf) then
        allocate(ncname(nvar,4))

        fname(7:9) = cexpnr
        call nctiminfo(ncname(1,:))
        call ncinfo(ncname( 2,:),'cfrac','Cloud fraction','-','time')
        call ncinfo(ncname( 3,:),'zb','Cloud-base height','m','time')
        call ncinfo(ncname( 4,:),'zc_av','Average Cloud-top height','m','time')
        call ncinfo(ncname( 5,:),'zc_max','Maximum Cloud-top height','m','time')
        call ncinfo(ncname( 6,:),'zi','Boundary layer height','m','time')
        call ncinfo(ncname( 7,:),'we','Entrainment velocity','m/s','time')
        call ncinfo(ncname( 8,:),'lwp_bar','Liquid-water path','kg/m^2','time')
        call ncinfo(ncname( 9,:),'lwp_max','Maximum Liquid-water path','kg/m^2','time')
        call ncinfo(ncname(10,:),'wmax','Maximum vertical velocity','m/s','time')
        call ncinfo(ncname(11,:),'vtke','Vertical integral of total TKE','kg/s^2','time')
        call ncinfo(ncname(12,:),'lmax','Maximum liquid water specific humidity','kg/kg','time')
        call ncinfo(ncname(13,:),'ustar','Surface friction velocity','m/s','time')
        call ncinfo(ncname(14,:),'tstr','Turbulent temperature scale','K','time')
        call ncinfo(ncname(15,:),'qtstr','Turbulent humidity scale','K','time')
        call ncinfo(ncname(16,:),'obukh','Obukhov Length','m','time')
        call ncinfo(ncname(17,:),'thlskin','Surface liquid water potential temperature','K','time')
        call ncinfo(ncname(18,:),'z0','Roughness height','m','time')
        call ncinfo(ncname(19,:),'wtheta','Surface kinematic temperature flux','K m/s','time')
        call ncinfo(ncname(20,:),'wthetav','Surface kinematic virtual temperature flux','K m/s','time')
        call ncinfo(ncname(21,:),'wq','Surface kinematic moisture flux','kg/kg m/s','time')
        call ncinfo(ncname(22,:),'twp_bar','Total water path','kg/m^2','time')
        call ncinfo(ncname(23,:),'rwp_bar','Rain water path','kg/m^2','time')
        call ncinfo(ncname(24,:),'pr','surface precipitation rate','kg/m^2/s','time')

        if(isurf==1) then
          call ncinfo(ncname(25,:),'Qnet','Net radiation','W/m^2','time')
          call ncinfo(ncname(26,:),'H','Sensible heat flux','W/m^2','time')
          call ncinfo(ncname(27,:),'LE','Latent heat flux','W/m^2','time')
          call ncinfo(ncname(28,:),'G0','Ground heat flux','W/m^2','time')
          call ncinfo(ncname(29,:),'tendskin','Skin tendency','W/m^2','time')
          call ncinfo(ncname(30,:),'rs','Surface resistance','s/m','time')
          call ncinfo(ncname(31,:),'ra','Aerodynamic resistance','s/m','time')
          call ncinfo(ncname(32,:),'cliq','Fraction of vegetated surface covered with liquid water','-','time')
          call ncinfo(ncname(33,:),'Wl','Liquid water reservoir','m','time')
          call ncinfo(ncname(34,:),'rssoil','Soil evaporation resistance','s/m','time')
          call ncinfo(ncname(35,:),'rsveg','Vegitation resistance','s/m','time')
        else if (isurf == 11) then
          call ncinfo(ncname(25,:),'Qnet','Net radiation','W/m^2','time')
          call ncinfo(ncname(26,:),'H','Sensible heat flux','W/m^2','time')
          call ncinfo(ncname(27,:),'LE','Latent heat flux','W/m^2','time')
          call ncinfo(ncname(28,:),'G','Ground heat flux','W/m^2','time')
          call ncinfo(ncname(29,:),'f1','Reduction canopy resistance f(swd)','-','time')
          call ncinfo(ncname(30,:),'f2b','Reduction soil resistance f(theta)','-','time')
          call ncinfo(ncname(31,:),'wl','Liquid water reservoir','m','time')

          if (lags) then
            call ncinfo(ncname(32,:),'an_co2','Net CO2 assimilation','ppb m s-1','time')
            call ncinfo(ncname(33,:),'resp_co2','CO2 respiration soil','ppb m s-1','time')
          end if
        end if

        if (iradiation /= 0) then
          call ncinfo(ncname(ivar_rad+ 0,:),'rlds',   'surface downwelling longwave flux','W/m^2','time')
          call ncinfo(ncname(ivar_rad+ 1,:),'rlus',   'surface upwelling longwave flux','W/m^2','time')
          call ncinfo(ncname(ivar_rad+ 2,:),'rsds',   'surface downwelling shortwave flux','W/m^2','time')
          call ncinfo(ncname(ivar_rad+ 3,:),'rsus',   'surface upwelling shortwave flux','W/m^2','time')
          call ncinfo(ncname(ivar_rad+ 4,:),'rsdscs', 'surface downwelling shortwave flux - clear sky','W/m^2','time')
          call ncinfo(ncname(ivar_rad+ 5,:),'rsuscs', 'surface upwelling shortwave flux - clear sky','W/m^2','time')
          call ncinfo(ncname(ivar_rad+ 6,:),'rldscs', 'surface downwelling longwave flux - clear sky','W/m^2','time')
          call ncinfo(ncname(ivar_rad+ 7,:),'rluscs', 'surface upwelling longwave flux - clear sky','W/m^2','time')

          call ncinfo(ncname(ivar_rad+ 8,:),'rsdt',   'TOA incoming shortwave flux','W/m^2','time')
          call ncinfo(ncname(ivar_rad+ 9,:),'rsut',   'TOA outgoing shortwave flux','W/m^2','time')
          call ncinfo(ncname(ivar_rad+10,:),'rlut',   'TOA outgoing longwave flux','W/m^2','time')
          call ncinfo(ncname(ivar_rad+11,:),'rsutcs', 'TOA outgoing shortwave flux -clear sky','W/m^2','time')
          call ncinfo(ncname(ivar_rad+12,:),'rlutcs', 'TOA outgoing longwave flux -clear sky','W/m^2','time')

          call ncinfo(ncname(ivar_rad+13,:),'rsdtm',  'TOM incoming shortwave flux','W/m^2','time')
          call ncinfo(ncname(ivar_rad+14,:),'rldtm',  'TOM incoming longwave flux','W/m^2','time')
          call ncinfo(ncname(ivar_rad+15,:),'rsutm',  'TOM outgoing shortwave flux','W/m^2','time')
          call ncinfo(ncname(ivar_rad+16,:),'rlutm',  'TOM outgoing longwave flux','W/m^2','time')
          call ncinfo(ncname(ivar_rad+17,:),'rsutmcs','TOM outgoing shortwave flux -clear sky','W/m^2','time')
          call ncinfo(ncname(ivar_rad+18,:),'rlutmcs','TOM outgoing longwave flux -clear sky','W/m^2','time')
        end if

        call open_nc(fname,  ncid,nrec)
        if(nrec==0) call define_nc( ncid, NVar, ncname)
      end if
    end if

    !$acc enter data copyin(blh_fld, sv0h, profile, gradient, dgrad)

    call timer_toc('modtimestat/inittimestat')

  end subroutine inittimestat

!>Run timestat. Calculate and write the statistics
  subroutine timestat

    use modglobal,  only : i1,j1, k1,kmax,zf,dzf,cu,cv,rv,rd,eps1, &
                          ijtot,timee,rtimee,dt_lim,rk3step,cexpnr,ifoutput
    use modmicrodata, only : imicro, imicro_sice, imicro_sice2, imicro_bulk, imicro_bulk3, precep
    use modfields,  only : e120,qt0,ql0,u0av,v0av,rhobf,rhof,u0,v0,w0,sv0
    use modsurfdata,only : wtsurf, wqsurf, isurf,ustar,thlflux,qtflux,z0,oblav,qts,thls,&
                           Qnet, H, LE, G0, rs, ra, tskin, tendskin, &
                           cliq,rsveg,rssoil,Wl, &
                           obl, wco2av, Anav, Respav,gcco2av
    use modmpi,     only : mpi_sum,mpi_max,mpi_min,comm3d,mpierr,myid, D_MPI_ALLREDUCE
    use modstat_nc,  only : lnetcdf, writestat_nc,nc_fillvalue
    use modlsm,     only : tile, f1, f2b, nlu, lags, an_co2, resp_co2
#if defined(_OPENACC)
    use modgpu, only: update_host
#endif
    use modraddata, only :  lwd,lwu,swd,swu,lwdca,lwuca,swdca,swuca, &
                            iradiation, doclearsky
    use modtracers, only : get_tracer_index
    implicit none

    real(field_r)   :: zbaseavl, ztopavl, ztopmaxl, ztop, zbaseminl
    real(field_r)   :: qlintavl, qlintmaxl, tke_totl
    real(field_r)   :: qrintavl, qtintavl
    real(field_r)   :: ccl, wmaxl, qlmaxl
    real(field_r)   :: ust,tst,qst,ustl,tstl,qstl,thlfluxl,qtfluxl
    real(field_r)   :: usttst, ustqst
    real(field_r)   :: wts, wqls,wthvs
    real(field_r)   :: c1,c2 !Used to calculate wthvs
    real,dimension(nvar) :: vars

    ! lsm variables
    real   :: Qnetavl, Havl, LEavl, G0avl, tendskinavl, rsavl, raavl, tskinavl,Wlavl,cliqavl,rsvegavl,rssoilavl
    real   :: Qnetav, Hav, LEav, G0av, tendskinav, rsav, raav, tskinav,Wlav,cliqav,rsvegav,rssoilav

    ! LSM tiled variables
    real   :: obuk_av(nlu)
    real   :: ustar_av(nlu)
    real   :: ra_av(nlu)
    real   :: f1_av, f2_av(nlu), f3_av(nlu), f2b_av
    real   :: rs_av(nlu)
    real   :: c_av(nlu)
    real   :: H_av(nlu)
    real   :: LE_av(nlu)
    real   :: G_av(nlu)
    real   :: thlskin_av(nlu)
    real   :: qtskin_av(nlu)
    real   :: an_co2_av, resp_co2_av

    ! Radiation variables for reductions
    real(field_r) :: &
      s_lwd_surf,    & !< Surface downwelling longwave flux
      s_lwu_surf,    & !< Surface upwelling longwave flux
      s_swd_surf,    & !< Surface downwelling shortwave flux
      s_swu_surf,    & !< Surface upwelling shortwave flux
      s_swd_toa,     & !< TOA downwelling shortwave flux
      s_swu_toa,     & !< TOA upwelling shortwave flux
      s_lwu_toa,     & !< TOA upwelling longwave flux
      s_swd_tom,     & !< TOM downwelling shortwave flux
      s_lwd_tom,     & !< TOM downwelling longwave flux
      s_swu_tom,     & !< TOM upwelling shortwave flux
      s_lwu_tom,     & !< TOM upwelling longwave flux
      s_swd_surf_ca, & !< Surface downwelling shorwave flux, clear sky
      s_swu_surf_ca, & !< Surface upwelling shortwave flux, clear sky
      s_lwd_surf_ca, & !< Surface downwelling longwave flux, clear sky
      s_lwu_surf_ca, & !< Surface upwelling longwave flux, clear sky
      s_swu_toa_ca,  & !< TOA upwelling shortwave flux, clear sky
      s_lwu_toa_ca,  & !< TOA upwelling longwave flux, clear sky
      s_swu_tom_ca,  & !< TOM upwelling shortwave flux, clear sky
      s_lwu_tom_ca     !< TOM upwelling longwave flux, clear sky

    integer:: i, j, k, ilu, iqr

    if (.not.(ltimestat)) return
    if (rk3step/=3) return
    if(timee<tnext) then
      dt_lim = min(dt_lim,tnext-timee)
      return
    end if

    call timer_tic('modtimestat/timestat', 0)

    tnext = tnext+idtav
    dt_lim = minval((/dt_lim,tnext-timee/))

    !      -----------------------------------------------------------
  !     1     EVALUATION OF CLOUD COVER, CLOUD BASE, ETC.
  !    -----------------------------------------------------------

  !     -----------------------------------------------------
  !     1.   Set A:  entrainment and time evolution
  !     -----------------------------------------------------

    zbaseavl = 0.0
    ztopavl = 0.0
    zbaseminl = zf(kmax)
    store_zi = .true.

    call calcblheight

    store_zi = .false.

  !     --------------------------------------------------------------
  !     9.2  liq. waterpath, cloudcover, cloudbase and cloudtop
  !     --------------------------------------------------------------

    ccl      = 0.0
    qlintavl = 0.0
    qtintavl = 0.0
    qrintavl = 0.0
    qlintmaxl= 0.0
    tke_totl = 0.0

    ! Make qlint 2D array to get rid of the if statement?
    !$acc parallel loop collapse(2) default(present) reduction(+:ccl, qlintavl, qtintavl) &
    !$acc& reduction(max: qlintmaxl) private(qlint, qtint) async
    do j = 2, j1
      do i = 2, i1
        qlint = 0.
        qtint = 0.
        !$acc loop reduction(+: qlint, qtint)
        do k = 1, kmax
          qlint = qlint + ql0(i,j,k)*rhof(k)*dzf(k)
          qtint = qtint + qt0(i,j,k)*rhof(k)*dzf(k)
        end do
        if (qlint > 0.) then
          ccl = ccl + 1
          qlintavl = qlintavl + qlint
          qlintmaxl = max(qlint, qlintmaxl)
        end if
        qtintavl = qtintavl + qtint
      end do
    end do

    if (imicro == imicro_sice .or. imicro == imicro_sice2 .or. imicro == imicro_bulk .or. imicro == imicro_bulk3) then
       iqr = get_tracer_index("qr")
       if (iqr == 0) then
          iqr = get_tracer_index("qhr")
       endif
       !$acc parallel loop collapse(2) default(present) reduction(+:qrintavl) &
      !$acc& private(qrint) async
      do j = 2, j1
        do i = 2, i1
          qrint = 0.0
          !$acc loop reduction(+: qrint)
          do k = 1, kmax
            qrint = qrint + sv0(i, j, k, iqr) * rhof(k) * dzf(k)
          end do
          qrintavl = qrintavl + qrint
        end do
      end do
    end if

    !$acc parallel loop collapse(2) default(present) reduction(+: zbaseavl) reduction(min: zbaseminl) async
    do j = 2, j1
      do i = 2, i1
        !$acc loop seq
        do k = 1, kmax
          if (ql0(i,j,k) > 0.) then
            zbaseavl = zbaseavl + zf(k)
            zbaseminl = min(zf(k),zbaseminl)
            exit
          end if
        end do
      end do
    end do

  !     ---------------------------------------
  !     9.3  determine maximum ql_max and w_max
  !     ---------------------------------------

    wmaxl  = 0.0
    qlmaxl = 0.0
    ztopavl = 0.0
    ztopmaxl = 0.0

    !$acc parallel loop collapse(2) default(present) reduction(+:ztopavl) private(ztop)
    do j = 2, j1
      do i = 2, i1
        ztop = 0.0
        !$acc loop seq
        do k = 1, kmax
          if (ql0(i,j,k) > 0.0) then
            ztop = zf(k)
          endif
          wmaxl = max(w0(i,j,k), wmaxl)
          qlmaxl = max(ql0(i,j,k), qlmaxl)
        end do
        ztopavl = ztopavl + ztop
        if (ztop > ztopmaxl) ztopmaxl = ztop
      end do
    end do

  !     -------------------------
  !     9.5  Domain Averaged TKE
  !     -------------------------

    !$acc parallel loop collapse(3) default(present) reduction(+: tke_totl) async
    do k = 1, kmax
      do j = 2, j1
        do i = 2, i1
          tke_totl = tke_totl +(0.5*( &
                               (0.5*(u0(i,j,k)+u0(i+1,j,k))+cu-u0av(k))**2 &
                              +(0.5*(v0(i,j,k)+v0(i,j+1,k))+cv-v0av(k))**2 &
                              +(0.5*(w0(i,j,k)+w0(i,j,k+1))           )**2 &
                                    ) + e120(i,j,k)**2 ) * dzf(k) * rhof(k)
        end do
      end do
    end do



!     -------------------------
!     9.6  Horizontally  Averaged ustar, tstar and obl
!     -------------------------

    ustl = 0
    tstl = 0
    qstl = 0
    !$acc parallel loop collapse(2) default(present) reduction(+:ustl,tstl,qstl) async
    do j = 2, j1
       do i = 2, i1
          ustl = ustl + ustar(i,j)
          tstl = tstl - thlflux(i,j) / ustar(i,j)
          qstl = qstl - qtflux (i,j) / ustar(i,j)
       end do
    end do

    if(isurf < 3) then
       thlfluxl = 0
       qtfluxl  = 0
       !$acc parallel loop collapse(2) default(present) reduction(+:thlfluxl,qtfluxl) async
       do j = 2, j1
          do i = 2, i1
             thlfluxl = thlfluxl + thlflux(i, j)
             qtfluxl  = qtfluxl  + qtflux (i, j)
          end do
       end do
    end if
    ! note ! ACC wait is far below

  ! -----------------------------------
  ! 9.7 Communication and normalisation
  ! -----------------------------------

    call D_MPI_ALLREDUCE(ccl   , cc   , 1,       &
                          MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(qlintavl, qlintav, 1  , &
                          MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(qtintavl, qtintav, 1  , &
                          MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(qrintavl, qrintav, 1  , &
                          MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(qlintmaxl, qlintmax, 1, &
                          MPI_MAX, comm3d,mpierr)
    call D_MPI_ALLREDUCE(zbaseavl, zbaseav, 1,   &
                          MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(zbaseminl, zbasemin, 1, &
                          MPI_MIN, comm3d,mpierr)
    prav = 0
    if (imicro == imicro_sice .or. imicro == imicro_sice2 .or. imicro == imicro_bulk) then
       pravl = 0
       !$acc parallel loop collapse(2) default(present) reduction(+:pravl)
       do j = 2, j1
          do i = 2, i1
             pravl = pravl + precep(i,j,1)
          end do
       end do

       call D_MPI_ALLREDUCE(pravl, prav, 1, MPI_SUM, comm3d,mpierr)
    end if

    call D_MPI_ALLREDUCE(wmaxl   , wmax   , 1,   &
                          MPI_MAX, comm3d,mpierr)
    call D_MPI_ALLREDUCE(qlmaxl, qlmax, 1,       &
                          MPI_MAX, comm3d,mpierr)
    call D_MPI_ALLREDUCE(ztopavl, ztopav, 1,     &
                          MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(ztopmaxl, ztopmax, 1,   &
                          MPI_MAX, comm3d,mpierr)

    if (cc > 0.0) then
      zbaseav = zbaseav / cc
      ztopav  = ztopav / cc
    else
      zbaseav = 0.0
      ztopav = 0.0
    end if

    cc      = cc / ijtot
    qlintav = qlintav / ijtot !domain averaged liquid water path
    qtintav = qtintav / ijtot !domain averaged total water path
    qrintav = qrintav / ijtot !domain averaged rain water path
    prav    = prav*rhobf(1) / ijtot !domain averaged precipitation rate

    call D_MPI_ALLREDUCE(tke_totl, tke_tot, 1,   &
                          MPI_SUM, comm3d,mpierr)

    tke_tot = tke_tot / ijtot

    !$acc wait ! wait for sum of ustl etc and thlfluxl,qtfluxl
    call D_MPI_ALLREDUCE(ustl, ust, 1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(tstl, tst, 1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(qstl, qst, 1, MPI_SUM, comm3d,mpierr)

    if (isurf < 3) then
      call D_MPI_ALLREDUCE(thlfluxl, usttst, 1, MPI_SUM, comm3d,mpierr)
      call D_MPI_ALLREDUCE(qtfluxl,  ustqst, 1, MPI_SUM, comm3d,mpierr)
      usttst = -usttst / ijtot
      ustqst = -ustqst / ijtot
    end if

    ust = ust / ijtot
    tst = tst / ijtot
    qst = qst / ijtot

    !Constants c1 and c2
    c1   = 1.+(rv/rd-1)*qts
    c2   = (rv/rd-1)

    if(isurf >= 3) then
      wts  = wtsurf
      wqls = wqsurf
      wthvs = c1*wts + c2*thls*wqls
    else
      wts  = -usttst
      wqls = -ustqst
      wthvs = c1*wts + c2*thls*wqls
    end if

  !  9.8  Create statistics for the land surface scheme
    if(isurf == 1) then
      Qnetavl      = sum(Qnet(2:i1,2:j1))
      Havl         = sum(H(2:i1,2:j1))
      LEavl        = sum(LE(2:i1,2:j1))
      G0avl        = sum(G0(2:i1,2:j1))
      tendskinavl  = sum(tendskin(2:i1,2:j1))
      rsavl        = sum(rs(2:i1,2:j1))
      raavl        = sum(ra(2:i1,2:j1))
      cliqavl      = sum(cliq(2:i1,2:j1))
      Wlavl        = sum(wl(2:i1,2:j1))
      rsvegavl     = sum(rsveg(2:i1,2:j1))
      rssoilavl    = sum(rssoil(2:i1,2:j1))
      tskinavl     = sum(tskin(2:i1,2:j1))

      call D_MPI_ALLREDUCE(Qnetavl,     Qnetav,     1,  MPI_SUM, comm3d,mpierr)
      call D_MPI_ALLREDUCE(Havl,        Hav,        1,  MPI_SUM, comm3d,mpierr)
      call D_MPI_ALLREDUCE(LEavl,       LEav,       1,  MPI_SUM, comm3d,mpierr)
      call D_MPI_ALLREDUCE(G0avl,       G0av,       1,  MPI_SUM, comm3d,mpierr)
      call D_MPI_ALLREDUCE(tendskinavl, tendskinav, 1,  MPI_SUM, comm3d,mpierr)
      call D_MPI_ALLREDUCE(rsavl,       rsav,       1,  MPI_SUM, comm3d,mpierr)
      call D_MPI_ALLREDUCE(raavl,       raav,       1,  MPI_SUM, comm3d,mpierr)
      call D_MPI_ALLREDUCE(cliqavl,     cliqav,     1,  MPI_SUM, comm3d,mpierr)
      call D_MPI_ALLREDUCE(wlavl,       wlav,       1,  MPI_SUM, comm3d,mpierr)
      call D_MPI_ALLREDUCE(rsvegavl,    rsvegav,    1,  MPI_SUM, comm3d,mpierr)
      call D_MPI_ALLREDUCE(rssoilavl,   rssoilav,   1,  MPI_SUM, comm3d,mpierr)
      call D_MPI_ALLREDUCE(tskinavl,    tskinav,    1,  MPI_SUM, comm3d,mpierr)

      Qnetav        = Qnetav      / ijtot
      Hav           = Hav         / ijtot
      LEav          = LEav        / ijtot
      G0av          = G0av        / ijtot
      tendskinav    = tendskinav  / ijtot
      rsav          = rsav        / ijtot
      raav          = raav        / ijtot
      cliqav        = cliqav      / ijtot
      wlav          = wlav        / ijtot
      rsvegav       = rsvegav     / ijtot
      rssoilav      = rssoilav    / ijtot
      tskinav       = tskinav     / ijtot

    else if (isurf == 11) then
      Qnet(2:i1,2:j1) = swd(2:i1,2:j1,1) + swu(2:i1,2:j1,1) + lwd(2:i1,2:j1,1) + lwu(2:i1,2:j1,1)

      ! TODO: replace mean_2d with slabsum?
      Qnetav = mean_2d(Qnet)
      Hav    = mean_2d(H)
      LEav   = mean_2d(LE)
      G0av   = mean_2d(G0)
      oblav  = mean_2d(obl)

      ! Tiled variables
      obuk_av = 0
      ustar_av = 0
      ra_av = 0
      f2_av = 0
      f3_av = 0
      do ilu=1,nlu
        !skip for ws and slb
        if (trim(tile(ilu)%lushort) == 'ws'.or. &
            trim(tile(ilu)%lushort) == 'slb') cycle
        obuk_av(ilu)  = mean_2d(tile(ilu)%obuk)
        ustar_av(ilu) = mean_2d(tile(ilu)%ustar)
        ra_av(ilu)    = mean_2d(tile(ilu)%ra)
        f2_av(ilu)    = mean_2d(tile(ilu)%f2)
        f3_av(ilu)    = mean_2d(tile(ilu)%f3)
      end do

      do ilu=1,nlu
        !skip for ws, aq and slb
        if (trim(tile(ilu)%lushort) == 'ws' .or. &
            trim(tile(ilu)%lushort) == 'aq' .or. &
            trim(tile(ilu)%lushort) == 'slb') cycle 
        rs_av(ilu)    = mean_2d(tile(ilu)%rs)
      end do

      f1_av    = mean_2d(f1)
      f2b_av   = mean_2d(f2b)

      do ilu=1,nlu
        ! skip for slb
        if (trim(tile(ilu)%lushort) == 'slb') cycle
        c_av(ilu)       = mean_2d(tile(ilu)%frac)
        H_av(ilu)       = mean_2d(tile(ilu)%H)
        LE_av(ilu)      = mean_2d(tile(ilu)%LE)
        thlskin_av(ilu) = mean_2d(tile(ilu)%thlskin) 
        qtskin_av(ilu)  = mean_2d(tile(ilu)%qtskin)  !TODO urb skin roof/can
      end do

      wlav = mean_2d(wl)

      do ilu=1,nlu
        !skip for aq and slb
        if (trim(tile(ilu)%lushort) == 'aq' .or. &
            trim(tile(ilu)%lushort) == 'slb') cycle
        G_av(ilu)    = mean_2d(tile(ilu)%G)
      end do

      if (lags) then
        an_co2_av   = mean_2d(an_co2)
        resp_co2_av = mean_2d(resp_co2)
      endif
     end if

    ! calculate radiation fluxes at surface, TOM, TOA
    if (iradiation /= 0) then

      s_lwd_surf = 0
      s_lwu_surf = 0
      s_swd_surf = 0
      s_swu_surf = 0
      s_swd_toa = 0
      s_swu_toa = 0
      s_lwu_toa = 0
      s_swd_tom = 0
      s_lwd_tom = 0
      s_swu_tom = 0
      s_lwu_tom = 0
      s_swd_surf_ca = 0
      s_swu_surf_ca = 0
      s_lwd_surf_ca = 0
      s_lwu_surf_ca = 0
      s_swu_toa_ca = 0
      s_lwu_toa_ca = 0
      s_swu_tom_ca = 0
      s_lwu_tom_ca = 0

      ! Surface fluxes

      !$acc parallel loop gang vector collapse(2) default(present) &
      !$acc reduction(+: s_swd_surf, s_swu_surf, s_lwd_surf, s_lwu_surf) async
      do j = 2, j1
        do i = 2, i1
          s_swd_surf = s_swd_surf + swd(i,j,1)
          s_swu_surf = s_swu_surf + swu(i,j,1)
          s_lwd_surf = s_lwd_surf + lwd(i,j,1)
          s_lwu_surf = s_lwu_surf + lwu(i,j,1)
        end do
      end do

      ! Top of atmosphere fluxes

      !$acc parallel loop gang vector collapse(2) default(present) &
      !$acc reduction(+: s_swd_toa, s_swu_toa, s_lwu_toa) async
      do j = 2, j1
        do i = 2, i1
          s_swd_toa = s_swd_toa + swd(i,j,k1)
          s_swu_toa = s_swu_toa + swu(i,j,k1)
          s_lwu_toa = s_lwu_toa + lwu(i,j,k1)
        end do
      end do

      ! Top of model fluxes

      !$acc parallel loop gang vector collapse(2) default(present) &
      !$acc reduction(+: s_swd_tom, s_swu_tom, s_lwd_tom, s_lwu_tom) async
      do j = 2, j1
        do i = 2, i1
          s_swd_tom = s_swd_tom + swd(i,j,kmax)
          s_swu_tom = s_swu_tom + swu(i,j,kmax)
          s_lwd_tom = s_lwd_tom + lwd(i,j,kmax)
          s_lwu_tom = s_lwu_tom + lwu(i,j,kmax)
        end do
      end do

      if (doclearsky) then

        ! Surface fluxes

        !$acc parallel loop gang vector collapse(2) default(present) &
        !$acc reduction(+: s_swd_surf_ca, s_swu_surf_ca, &
        !$acc              s_lwd_surf_ca, s_lwu_surf_ca) async
        do j = 2, j1
          do i = 2, i1
            s_swd_surf_ca = s_swd_surf_ca + swdca(i,j,1)
            s_swu_surf_ca = s_swu_surf_ca + swuca(i,j,1)
            s_lwd_surf_ca = s_lwd_surf_ca + lwdca(i,j,1)
            s_lwu_surf_ca = s_lwu_surf_ca + lwuca(i,j,1)
          end do
        end do

        ! Top of atmosphere fluxes

        !$acc parallel loop gang vector collapse(2) default(present) &
        !$acc reduction(+: s_swu_toa_ca, s_lwu_toa_ca) async
        do j = 2, j1
          do i = 2, i1
            s_swu_toa_ca = s_swu_toa_ca + swuca(i,j,k1)
            s_lwu_toa_ca = s_lwu_toa_ca + lwuca(i,j,k1)
          end do
        end do

        ! Top of model fluxes

        !$acc parallel loop gang vector collapse(2) default(present) &
        !$acc reduction(+: s_swu_tom_ca, s_lwu_tom_ca) async
        do j = 2, j1
          do i = 2, i1
            s_swu_tom_ca = s_swu_tom_ca + swuca(i,j,kmax)
            s_lwu_tom_ca = s_lwu_tom_ca + lwuca(i,j,kmax)
          end do
        end do
      end if

      !$acc wait

      vars(ivar_rad+ 0) = abs(s_lwd_surf)  !'rlds',   'surface downwelling longwave flux'
      vars(ivar_rad+ 1) = abs(s_lwu_surf)  !'rlus',   'surface upwelling longwave flux'
      vars(ivar_rad+ 2) = abs(s_swd_surf)  !'rsds',   'surface downwelling shortwave flux'
      vars(ivar_rad+ 3) = abs(s_swu_surf)  !'rsus',   'surface upwelling shortwave flux'

      vars(ivar_rad+ 4) = abs(s_swd_surf_ca)  !'rsdscs', 'surface downwelling shortwave flux - clear sky'
      vars(ivar_rad+ 5) = abs(s_swu_surf_ca)  !'rsuscs', 'surface upwelling shortwave flux - clear sky'
      vars(ivar_rad+ 6) = abs(s_lwd_surf_ca)  !'rldscs', 'surface downwelling longwave flux - clear sky'
      vars(ivar_rad+ 7) = abs(s_lwu_surf_ca)  !'rluscs', 'surface upwelling longwave flux - clear sky'

      vars(ivar_rad+ 8) = abs(s_swd_toa)    !'rsdt',   'TOA incoming shortwave flux','W/m^2'
      vars(ivar_rad+ 9) = abs(s_swu_toa)    !'rsut',   'TOA outgoing shortwave flux','W/m^2'
      vars(ivar_rad+10) = abs(s_lwu_toa)    !'rlut',   'TOA outgoing longwave flux','W/m^2'
      vars(ivar_rad+11) = abs(s_swu_toa_ca) !'rsutcs', 'TOA outgoing shortwave flux -clear sky'
      vars(ivar_rad+12) = abs(s_lwu_toa_ca) !'rlutcs', 'TOA outgoing longwave flux -clear sky'

      vars(ivar_rad+13) = abs(s_swd_tom)   !'rsdtm',  'TOM incoming shortwave flux'
      vars(ivar_rad+14) = abs(s_lwd_tom)   !'rsdtm',  'TOM incoming longwave flux'
      vars(ivar_rad+15) = abs(s_swu_tom)   !'rsutm',  'TOM outgoing shortwave flux'
      vars(ivar_rad+16) = abs(s_lwu_tom)   !'rlutm',  'TOM outgoing longwave flux'
      vars(ivar_rad+17) = abs(s_swu_tom_ca) !'rsutmcs','TOM outgoing shortwave flux -clear sky'
      vars(ivar_rad+18) = abs(s_lwu_tom_ca) !'rlutmcs','TOM outgoing longwave flux -clear sky'

      call D_MPI_ALLREDUCE(vars(ivar_rad:ivar_rad+18), 19,  MPI_SUM, comm3d,mpierr)

      vars(ivar_rad:ivar_rad+18) = vars(ivar_rad:ivar_rad+18) / ijtot
    end if


  !  9.8  write the results to output file
  !     ---------------------------------------

    if(myid==0)then
       !tmser1
      open (ifoutput,file='tmser1.'//cexpnr,position='append')
      write( ifoutput,'(f10.2,f6.3,4f12.3,f10.4,5f9.3)') &
          rtimee, &
          cc, &
          zbaseav, &
          ztopav, &
          ztopmax, &
          zi, &
          we, &
          qlintav*1000., &
          qlintmax*1000., &
          wmax, &
          tke_tot, &
          qlmax*1000.
      close(ifoutput)

      !tmsurf
      open (ifoutput,file='tmsurf.'//cexpnr,position='append')
      write( ifoutput,'(f10.2,4e11.3,f11.3,4e11.3)') &
          rtimee   ,&
          ust     ,&
          tst     ,&
          qst     ,&
          oblav   ,&
          thls    ,&
          z0      ,&
          wts     ,&
          wthvs    ,&
          wqls
      close(ifoutput)

      if (isurf == 1) then
        !tmlsm
        open (ifoutput,file='tmlsm.'//cexpnr,position='append')
        write(ifoutput,'(f10.2,9f11.3,e13.3, 5f11.3,e13.3)') &
            rtimee       ,&
            Qnetav      ,&
            Hav         ,&
            LEav        ,&
            G0av        ,&
            tendskinav  ,&
            rsav        ,&
            raav        ,&
            tskinav     ,&
            cliqav      ,&
            wlav        ,&
            rssoilav    ,&
            rsvegav     ,&
            Respav      ,&
            wco2av      ,&
            Anav        ,&
            gcco2av
        close(ifoutput)
      end if
      if (lnetcdf) then
        vars( 1) = rtimee
        vars( 2) = cc
        vars( 3) = zbaseav
        if (vars(3)<eps1) vars(3) = nc_fillvalue
        vars( 4) = ztopav
        if (vars(4)<eps1) vars(4) = nc_fillvalue
        vars( 5) = ztopmax
        if (vars(5)<eps1) vars(5) = nc_fillvalue
        vars( 6) = zi
        vars( 7) = we
        vars( 8) = qlintav
        vars( 9) = qlintmax
        vars(10) = wmax
        vars(11) = tke_tot
        vars(12) = qlmax
        vars(13) = ust
        vars(14) = tst
        vars(15) = qst
        vars(16) = oblav
        vars(17) = thls
        vars(18) = z0
        vars(19) = wts
        vars(20) = wthvs
        vars(21) = wqls
        vars(22) = qtintav
        vars(23) = qrintav
        vars(24) = prav

        if (isurf == 1) then
          vars(25) = Qnetav
          vars(26) = Hav
          vars(27) = LEav
          vars(28) = G0av
          vars(29) = tendskinav
          vars(30) = rsav
          vars(31) = raav
          vars(32) = cliqav
          vars(33) = wlav
          vars(34) = rssoilav
          vars(35) = rsvegav
        else if (isurf == 11) then
          vars(25) = Qnetav
          vars(26) = Hav
          vars(27) = LEav
          vars(28) = G0av
          vars(29) = f1_av
          vars(30) = f2b_av
          vars(31) = wlav

          if (lags) then
            vars(32) = an_co2_av
            vars(33) = resp_co2_av
          end if
        end if

        call writestat_nc(ncid,nvar,ncname,vars,nrec,.true.)
      end if
    end if

    call timer_toc('modtimestat/timestat')

  end subroutine timestat

  function mean_2d(var_2d) result(res)
    use modglobal, only : i1, j1, ijtot
    use modmpi, only : mpi_sum, comm3d, mpierr, d_mpi_allreduce
    implicit none

    real, intent(in) :: var_2d(:,:)
    real :: res, var_sum_l, var_sum

    var_sum_l = sum(var_2d(2:i1, 2:j1))
    call d_mpi_allreduce(var_sum_l, var_sum, 1, mpi_sum, comm3d, mpierr)
    res = var_sum / ijtot
  end function mean_2d
    

!>Calculate the boundary layer height
!!
!! There are 3 available ways to calculate the boundary layer height:
!! - By determining the minimum flux in some scalar, e.g. buoyancy
!! - By determining the minimum local gradient of some scalar, averaged over a definable number of columns
!! - By monitoring a threshold value of some scalar, averaged over a definable number of columns
  subroutine calcblheight

    use modglobal,  only : i1,j1,kmax,k1,cp,rlv,imax,rd,zh,dzh,zf,dzf,rv,ijtot,iadv_sv,iadv_kappa
    use modfields,  only : w0,qt0,qt0h,ql0,thl0,thl0h,thv0h,sv0,exnf,whls
    use modsurfdata,only : svs
    use modmpi,     only : mpierr, comm3d,mpi_sum, D_MPI_ALLREDUCE
    use advec_kappa,only : halflev_kappa

    implicit none

    real    :: zil, dhdt, locval, oldlocval
    integer :: location, i, j, k, nsamp, stride
    real, allocatable,dimension(:,:,:) :: blh_fld2



    zil = 0.0
    !$acc kernels default(present)
    gradient(:) = 0.0
    dgrad(:) = 0.0
    !$acc end kernels

    select case (iblh_meth)
      case (iblh_flux)
        select case (iblh_var)
          case(iblh_qt)
            !$acc kernels default(present)
            blh_fld(:,:,:) = w0(:,:,:)*qt0h(:,:,:)
            !$acc end kernels
          case(iblh_thl)
            !$acc kernels default(present)
            blh_fld(:,:,:) = w0(:,:,:)*thl0h(:,:,:)
            !$acc end kernels
          case(iblh_thv)
            !$acc kernels default(present)
            blh_fld(:,:,:) = w0(:,:,:)*thv0h(:,:,:)
            !$acc end kernels
          case(1:)
            if (iadv_sv == iadv_kappa) then
              call halflev_kappa(sv0(:,:,:,iblh_var),sv0h)
              !$acc kernels default(present)
              sv0h(2:i1,2:j1,1) = svs(iblh_var)
              blh_fld(:,:,:) = w0(:,:,:)*sv0h(:,:,:)
              !$acc end kernels
            else
              !$acc kernels default(present)
              do k = 2, k1
                do j = 2, j1
                  do i = 2, i1
                    sv0h(i,j,k) = (sv0(i,j,k,iblh_var)*dzf(k-1)+sv0(i,j,k-1,iblh_var)*dzf(k))/(2*dzh(k))
                  enddo
                enddo
              enddo
              sv0h(2:i1,2:j1,1) = svs(iblh_var)
              blh_fld(:,:,:) = w0(:,:,:)*sv0h(:,:,:)
              !$acc end kernels
            end if
        end select

      case (iblh_grad,iblh_thres)
        select case (iblh_var)
          case(iblh_qt)
            !$acc kernels default(present)
            blh_fld(:,:,:) = qt0(:,:,:)
            !$acc end kernels
          case(iblh_thl)
            !$acc kernels default(present)
            blh_fld(:,:,:) = thl0(:,:,:)
            !$acc end kernels
          case(iblh_thv)
            !$acc parallel loop collapse(3) default(present) async
            do k = 1, k1
              do j = 2, j1
                do i = 2, i1
                  blh_fld(i,j,k) = (thl0(i,j,k)+rlv*ql0(i,j,k)/(cp*exnf(k))) &
                                   *(1+(rv/rd-1)*qt0(i,j,k)-rv/rd*ql0(i,j,k))
                end do
              end do
            end do
          case(1:)
            !$acc kernels default(present)
            blh_fld(:,:,:) = sv0(2:i1,2:j1,1:k1,iblh_var)
            !$acc end kernels
        end select

    end select

    select case (iblh_meth)
      case (iblh_flux)
        stride = ceiling(real(imax)/real(blh_nsamp))
        !$acc parallel loop collapse(2) default(present) reduction(+:zil) async
        do i = 2, stride+1
          do j = 2, j1
            nsamp =  ceiling(real(i1-i+1)/real(stride))
            zil = zil + nsamp*zh(minloc(sum(blh_fld(i:i1:stride,j,:),1),1))
          end do
        end do

      case (iblh_grad)
        stride = ceiling(real(imax)/real(blh_nsamp))
        !$acc parallel loop collapse(2) default(present) reduction(+:zil)
        do i = 2, stride+1
          do j = 2, j1
            nsamp =  ceiling(real(i1-i+1)/real(stride))
            profile(:) = sum(blh_fld(i:i1:stride,j,:),1)
            select case (iblh_var)
              case(iblh_qt) !Water vapour gradients near the inversion layer can be either positive or negative
                gradient(2:k1) = abs(profile(2:k1) - profile(1:kmax))/dzh(2:k1)
              case(iblh_thl,iblh_thv) !temperature jumps near the inversion layer are always positive
                gradient(2:k1) = (profile(2:k1) - profile(1:kmax))/dzh(2:k1)
              case default
                gradient(2:k1) = (profile(2:k1) - profile(1:kmax))/dzh(2:k1)
            end select
            dgrad(2:kmax)    = (gradient(3:k1) - gradient(2:kmax))/dzf(2:kmax)
            location = maxloc(gradient,1)
            zil  = zil + nsamp*(zh(location-1) - dzh(location)*dgrad(location-1)/(dgrad(location)-dgrad(location-1) + 1.e-8))
          enddo
        enddo

      case (iblh_thres)
        stride = ceiling(real(imax)/real(blh_nsamp))
        do i=2,stride+1
          nsamp =  ceiling(real(i1-i+1)/real(stride))
          do j=2,j1
            locval = 0.0
            do k=kmax,1,-1
              oldlocval = locval
              locval = blh_sign*sum(blh_fld(i:i1:stride,j,k))/nsamp
              if (locval < blh_sign*blh_thres) then
                zil = zil + nsamp *(zf(k) +  (blh_sign*blh_thres-locval) &
                          *dzh(k+1)/(oldlocval-locval))
                exit
              endif
            enddo
          enddo
        enddo

    end select

    call D_MPI_ALLREDUCE(zil, zi, 1, MPI_SUM, comm3d,mpierr)
    zi = zi / ijtot

    if (ziold< 0) ziold = zi
    dhdt = (zi-ziold)/dtav
    if(store_zi) ziold = zi

    k=2
    do while (zh(k)<zi .and. k < kmax)
      k=k+1
    end do
    we = dhdt - whls (k)   !include for large-scale vertical velocity

    !$acc wait

  end subroutine calcblheight

!> Clean up when leaving the run
  subroutine exittimestat
    use modmpi, only : myid
    use modstat_nc, only : exitstat_nc,lnetcdf
    implicit none

    if(ltimestat .and. lnetcdf .and. myid==0) call exitstat_nc(ncid)
    if(.not.ltimestat) return

    !$acc exit data delete(blh_fld, sv0h, profile, gradient, dgrad)

    deallocate(blh_fld,sv0h)
    deallocate(profile,gradient,dgrad)

  end subroutine exittimestat

end module modtimestat
