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
!  Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!


module modtimestat
  use modtimer
  use modprecision, only : longint, field_r

implicit none
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

  !Variables for heterogeneity
  real(field_r), allocatable :: u0av_patch (:,:)     ! patch averaged um    at full level
  real(field_r), allocatable :: v0av_patch (:,:)     ! patch averaged vm    at full level
  real(field_r), allocatable :: w0av_patch (:,:)     ! patch averaged wm    at full level
  real(field_r),allocatable, dimension(:,:) :: zbase_field, ztop_field, cc_field, qlint_field, tke_tot_field
  real(field_r),allocatable, dimension(:,:) :: zbase_patch, ztop_patch, zbasemin_patch, zbasemin_patchl
  real(field_r),allocatable, dimension(:,:) :: cc_patch, qlint_patch, qlintmax_patch, qlintmax_patchl, tke_tot_patch
  real(field_r),allocatable, dimension(:,:) :: wmax_patch, wmax_patchl, qlmax_patch, qlmax_patchl, ztopmax_patch, ztopmax_patchl
  real(field_r),allocatable, dimension(:,:) :: ust_patch, qst_patch, tst_patch, wthls_patch, wqls_patch, wthvs_patch
  !In combination with isurf = 1
  real(field_r),allocatable, dimension(:,:) :: Qnet_patch, H_patch, LE_patch, G0_patch, tendskin_patch,rs_patch,ra_patch
  real(field_r),allocatable, dimension(:,:) :: cliq_patch, wl_patch, rsveg_patch, rssoil_patch, tskin_patch, obl_patch
  real(field_r),allocatable, dimension(:,:) :: zi_patch,ziold_patch,we_patch, zi_field

contains
!> Initializing Timestat. Read out the namelist, initializing the variables
  subroutine inittimestat
    use modmpi,    only : myid,comm3d,mpierr,D_MPI_BCAST
    use modglobal, only : ifnamopt, fname_options,cexpnr,dtmax,ifoutput,dtav_glob,tres,&
                          ladaptive,k1,kmax,rd,rv,dt_lim,btime,i1,j1,lwarmstart,checknamelisterror, &
                          ih ,jh
    use modfields, only : thlprof,qtprof,svprof
    use modsurfdata, only : isurf, lhetero, xpatches, ypatches
    use modstat_nc, only : lnetcdf, open_nc, define_nc, ncinfo, nctiminfo
    use modraddata, only : iradiation
    use modlsm, only : lags
    implicit none
    integer :: ierr,k,location = 1
    integer :: i,j

    namelist/NAMTIMESTAT/ & !< namelist
    dtav,ltimestat,blh_thres,iblh_meth,iblh_var,blh_nsamp !! namelist contents

    call timer_tic('modtimestat/inittimestat', 0)

    dtav=dtav_glob
    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMTIMESTAT,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMTIMESTAT')
      write(6 ,NAMTIMESTAT)
      close(ifnamopt)
    end if

    call D_MPI_BCAST(dtav     ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(ltimestat  ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(blh_thres,1,0,comm3d,mpierr)
    call D_MPI_BCAST(iblh_meth  ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(iblh_var   ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(blh_nsamp  ,1,0,comm3d,mpierr)
    idtav = dtav/tres

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
    dt_lim = min(dt_lim,tnext)

    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'TIMESTAT: dtav should be a integer multiple of dtmax'
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
      stop 'TIMESTAT: Incorrect iblh_meth'
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

          if(lhetero) then
             do i=1,xpatches
                do j=1,ypatches
                   name = 'tmser1patchiiixjjj.'//cexpnr
                   write (name(12:14),'(i3.3)') i
                   write (name(16:18),'(i3.3)') j
                   open (ifoutput,file=name,status='replace',position='append')
                   write(ifoutput,'(2a)') &
                        '#  time      cc     z_cbase    z_ctop_avg  z_ctop_max      zi         we', &
                        '   <<ql>>  <<ql>>_max   w_max   tke     ql_max'
                   close(ifoutput)

                   name = 'tmsurfpatchiiixjjj.'//cexpnr
                   write (name(12:14),'(i3.3)') i
                   write (name(16:18),'(i3.3)') j
                   open (ifoutput,file=name,status='replace',position='append')
                   write(ifoutput,'(2a)') &
                        '#  time        ust        tst        qst         obukh', &
                        '      thls        z0        wthls      wthvs      wqls '
                   close(ifoutput)

                   if(isurf == 1) then
                      name = 'tmlsmpatchiiixjjj.'//cexpnr
                      write (name(11:13),'(i3.3)') i
                      write (name(15:17),'(i3.3)') j
                      open (ifoutput,file=name,status='replace',position='append')
                      write(ifoutput,'(3a)') &
                           '#     time      Qnet        H          LE         G0  ', &
                           '   tendskin     rs         ra        tskin        cliq  ', &
                           '    Wl          rssoil     rsveg'
                      write(ifoutput,'(3a)') &
                           '#      [s]     [W/m2]     [W/m2]     [W/m2]     [W/m2]', &
                           '   [W/m2]      [s/m]       [s/m]     [K]          [-]   ', &
                           '   [m]          [s/m]      [s/m]'
                      close(ifoutput)
                   endif
                enddo
             enddo
          endif
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

    if (lhetero) then
      allocate(zbase_field  (2:i1,2:j1))
      allocate(ztop_field   (2:i1,2:j1))
      allocate(cc_field     (2:i1,2:j1))
      allocate(qlint_field  (2:i1,2:j1))
      allocate(tke_tot_field(2:i1,2:j1))

      allocate(zbase_patch    (xpatches,ypatches))
      allocate(ztop_patch     (xpatches,ypatches))
      allocate(zbasemin_patch (xpatches,ypatches))
      allocate(zbasemin_patchl(xpatches,ypatches))
      allocate(cc_patch       (xpatches,ypatches))
      allocate(qlint_patch    (xpatches,ypatches))
      allocate(qlintmax_patch (xpatches,ypatches))
      allocate(qlintmax_patchl(xpatches,ypatches))
      allocate(tke_tot_patch  (xpatches,ypatches))
      allocate(wmax_patch     (xpatches,ypatches))
      allocate(wmax_patchl    (xpatches,ypatches))
      allocate(qlmax_patch    (xpatches,ypatches))
      allocate(qlmax_patchl   (xpatches,ypatches))
      allocate(ztopmax_patch  (xpatches,ypatches))
      allocate(ztopmax_patchl (xpatches,ypatches))

      allocate(u0av_patch(xpatches,ypatches))
      allocate(v0av_patch(xpatches,ypatches))
      allocate(w0av_patch(xpatches,ypatches))

      allocate(ust_patch (xpatches,ypatches))
      allocate(qst_patch (xpatches,ypatches))
      allocate(tst_patch (xpatches,ypatches))
      allocate(wthls_patch (xpatches,ypatches))
      allocate(wqls_patch(xpatches,ypatches))
      allocate(wthvs_patch(xpatches,ypatches))

      allocate(Qnet_patch    (xpatches,ypatches))
      allocate(H_patch       (xpatches,ypatches))
      allocate(LE_patch      (xpatches,ypatches))
      allocate(G0_patch      (xpatches,ypatches))
      allocate(tendskin_patch(xpatches,ypatches))
      allocate(rs_patch      (xpatches,ypatches))
      allocate(ra_patch      (xpatches,ypatches))
      allocate(cliq_patch    (xpatches,ypatches))
      allocate(wl_patch      (xpatches,ypatches))
      allocate(rsveg_patch   (xpatches,ypatches))
      allocate(rssoil_patch  (xpatches,ypatches))
      allocate(tskin_patch   (xpatches,ypatches))
      allocate(obl_patch     (xpatches,ypatches))

      allocate(zi_patch      (xpatches,ypatches))
      allocate(zi_field      (2:i1,2:j1))
      allocate(ziold_patch   (xpatches,ypatches))
      ziold_patch = -1
      allocate(we_patch      (xpatches,ypatches))
    endif

    !$acc enter data copyin(blh_fld, sv0h, profile, gradient, dgrad)

    call timer_toc('modtimestat/inittimestat')

  end subroutine inittimestat

!>Run timestat. Calculate and write the statistics
  subroutine timestat

    use modglobal,  only : i1,j1,kmax,zf,dzf,cu,cv,rv,rd,eps1, &
                          ijtot,timee,rtimee,dt_lim,rk3step,cexpnr,ifoutput
    use modmicrodata, only : imicro, iqr, imicro_sice, imicro_sice2, imicro_bulk, precep
    use modfields,  only : e120,qt0,ql0,u0av,v0av,rhobf,rhof,u0,v0,w0,sv0
    use modsurfdata,only : wtsurf, wqsurf, isurf,ustar,thlflux,qtflux,z0,oblav,qts,thls,&
                           Qnet, H, LE, G0, rs, ra, tskin, tendskin, &
                           cliq,rsveg,rssoil,Wl, &
                           lhetero, xpatches, ypatches, qts_patch, wt_patch, wq_patch, thls_patch,obl,z0mav_patch, wco2av, Anav, Respav,gcco2av
    use modsurface, only : patchxnr,patchynr
    use modmpi,     only : mpi_sum,mpi_max,mpi_min,comm3d,mpierr,myid, D_MPI_ALLREDUCE
    use modstat_nc,  only : lnetcdf, writestat_nc,nc_fillvalue
    use modlsm,     only : tile, f1, f2b, nlu, lags, an_co2, resp_co2
#if defined(_OPENACC)
    use modgpu, only: update_host
#endif
    use modraddata, only :  lwd,lwu,swd,swu,lwdca,lwuca,swdca,swuca,SW_up_TOA,SW_dn_TOA,LW_up_TOA,&
                            SW_up_ca_TOA,LW_up_ca_TOA,iradiation
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

    ! heterogeneity variables
    integer:: patchx, patchy

    integer:: i, j, k, ilu

    if (.not.(ltimestat)) return
    if (rk3step/=3) return
    if(timee<tnext) then
      dt_lim = min(dt_lim,tnext-timee)
      return
    end if

    call timer_tic('modtimestat/timestat', 0)

    tnext = tnext+idtav
    dt_lim = minval((/dt_lim,tnext-timee/))

    if (lhetero) then
      zbase_field    = 0
      ztop_field     = 0
      cc_field       = 0
      qlint_field    = 0
      tke_tot_field  = 0

      zbase_patch    = 0
      ztop_patch     = 0
      zbasemin_patch = zf(kmax)
      zbasemin_patchl= zf(kmax)
      cc_patch       = 0
      qlint_patch    = 0
      qlintmax_patch = 0
      qlintmax_patchl= 0
      tke_tot_patch  = 0
      wmax_patch     = 0
      wmax_patchl    = 0
      qlmax_patch    = 0
      qlmax_patchl   = 0
      ztopmax_patch  = 0
      ztopmax_patchl = 0

      ust_patch      = 0
      qst_patch      = 0
      tst_patch      = 0
      wthls_patch      = 0
      wqls_patch     = 0
      wthvs_patch     = 0

      Qnet_patch     = 0
      H_patch        = 0
      LE_patch       = 0
      G0_patch       = 0
      tendskin_patch = 0
      rs_patch       = 0
      ra_patch       = 0
      cliq_patch     = 0
      wl_patch       = 0
      rsveg_patch    = 0
      rssoil_patch   = 0
      tskin_patch    = 0
      obl_patch      = 0

      zi_patch       = 0
      zi_field       = 0
      we_patch       = 0
    endif

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
    if (lhetero) then
      do j = 2, j1
        do i = 2, i1
          patchy = patchynr(j)
          patchx = patchxnr(i)
          qlint = 0.
          qtint = 0.
          do k = 1, kmax
            qlint = qlint + ql0(i,j,k)*rhof(k)*dzf(k)
            qtint = qtint + qt0(i,j,k)*rhof(k)*dzf(k)
          end do
          if (qlint > 0.) then
            ccl = ccl + 1
            qlintavl = qlintavl + qlint
            qlintmaxl = max(qlint, qlintmaxl)
            cc_field(i,j)                  = 1.0
            qlint_field(i,j)               = qlint
            qlintmax_patchl(patchx,patchy) = max(qlintmax_patchl(patchx,patchy),qlint)
          end if
          qtintavl = qtintavl + qtint
          qrintavl = qrintavl + qrint
        end do
      end do
    else
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
    end if

    if (imicro > 0) then
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

    if (lhetero) then
      do j = 2, j1
        do i = 2, i1
          patchy = patchynr(j)
          patchx = patchxnr(i)
          do k = 1, kmax
            if (ql0(i,j,k) > 0.) then
              zbaseavl = zbaseavl + zf(k)
              zbaseminl = min(zf(k),zbaseminl)
              zbase_field(i,j) = zf(k)
              zbasemin_patchl(patchx,patchy) = min(zbasemin_patchl(patchx,patchy),zf(k))
              exit
            end if
          end do
        end do
      end do
    else
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
    end if

  !     ---------------------------------------
  !     9.3  determine maximum ql_max and w_max
  !     ---------------------------------------

    wmaxl  = 0.0
    qlmaxl = 0.0
    ztopavl = 0.0
    ztopmaxl = 0.0

    if (lhetero) then
      do j = 2, j1
        do i = 1, i1
          patchy = patchynr(j)
          patchx = patchxnr(i)
          ztop = 0.0
          do k = 1, kmax
            if (ql0(i,j,k) > 0.0) then
              ztop = zf(k)
              ztop_field(i,j) = zf(k)
            end if
            wmaxl = max(w0(i,j,k),wmaxl)
            qlmaxl = max(ql0(i,j,k),qlmaxl)
            wmax_patchl(patchx,patchy)  = max(wmax_patchl (patchx,patchy),w0 (i,j,k))
            qlmax_patchl(patchx,patchy) = max(qlmax_patchl(patchx,patchy),ql0(i,j,k))
          end do
          ztopavl = ztopavl + ztop
          if (ztop> ztopmaxl) ztopmaxl = ztop
          ztop_field = ztop
          if (ztop > ztopmax_patchl(patchx,patchy)) ztopmax_patchl(patchx,patchy) = ztop
        end do
      end do
    else
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
    end if

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

    if (lhetero) then
      do k = 1, kmax
        u0av_patch = patchsum_1level_field_r(u0(2:i1,2:j1,k)) * (xpatches*ypatches/ijtot)
        v0av_patch = patchsum_1level_field_r(v0(2:i1,2:j1,k)) * (xpatches*ypatches/ijtot)
        w0av_patch = patchsum_1level_field_r(w0(2:i1,2:j1,k)) * (xpatches*ypatches/ijtot)
        do j = 2, j1
          do i = 2, i1
            patchy = patchynr(j)
            patchx = patchxnr(i)

            tke_tot_field(i,j) = tke_tot_field(i,j) + (0.5*( &
                                     (0.5*(u0(i,j,k)+u0(i+1,j,k))+cu-u0av_patch(patchx,patchy))**2 + &
                                     (0.5*(v0(i,j,k)+v0(i,j+1,k))+cv-v0av_patch(patchx,patchy))**2 + &
                                     (0.5*(w0(i,j,k)+w0(i,j,k+1))   -w0av_patch(patchx,patchy))**2 &
                                         ) + e120(i,j,k)**2 ) * dzf(k) * rhof(k)
          end do
        end do
      end do
      tke_tot_patch = patchsum_1level_field_r(tke_tot_field) * (xpatches*ypatches/ijtot)
    end if



!     -------------------------
!     9.6  Horizontally  Averaged ustar, tstar and obl
!     -------------------------
    !$acc kernels default(present) async
    ustl=sum(ustar(2:i1,2:j1))
    tstl=sum(- thlflux(2:i1,2:j1) / ustar(2:i1,2:j1))
    qstl=sum(- qtflux (2:i1,2:j1) / ustar(2:i1,2:j1))
    !$acc end kernels

    if(isurf < 3) then
      !$acc kernels default(present) async
      thlfluxl = sum(thlflux(2:i1, 2:j1))
      qtfluxl  = sum(qtflux (2:i1, 2:j1))
      !$acc end kernels
    end if

  ! -----------------------------------
  ! 9.7 Communication and normalisation
  ! -----------------------------------
    !$acc wait
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
    if (imicro == imicro_sice .or. imicro == imicro_sice2 .or. imicro == imicro_bulk) then
       pravl = sum(precep(2:i1,2:j1,1))
       call D_MPI_ALLREDUCE(pravl, prav, 1, MPI_SUM, comm3d,mpierr)
    end if
    if (lhetero) then
      cc_patch    = patchsum_1level_field_r(cc_field   )
      qlint_patch = patchsum_1level_field_r(qlint_field)
      zbase_patch = patchsum_1level_field_r(zbase_field)
      call D_MPI_ALLREDUCE(qlintmax_patchl, qlintmax_patch, xpatches*ypatches, MPI_MAX, comm3d,mpierr)
      call D_MPI_ALLREDUCE(zbasemin_patchl, zbasemin_patch, xpatches*ypatches, MPI_MIN, comm3d,mpierr)
    endif

    call D_MPI_ALLREDUCE(wmaxl   , wmax   , 1,   &
                          MPI_MAX, comm3d,mpierr)
    call D_MPI_ALLREDUCE(qlmaxl, qlmax, 1,       &
                          MPI_MAX, comm3d,mpierr)
    call D_MPI_ALLREDUCE(ztopavl, ztopav, 1,     &
                          MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(ztopmaxl, ztopmax, 1,   &
                          MPI_MAX, comm3d,mpierr)

    if (lhetero) then
      call D_MPI_ALLREDUCE(wmax_patchl ,      wmax_patch, xpatches*ypatches, MPI_MAX, comm3d,mpierr)
      call D_MPI_ALLREDUCE(qlmax_patchl,     qlmax_patch, xpatches*ypatches, MPI_MAX, comm3d,mpierr)
      call D_MPI_ALLREDUCE(ztopmax_patchl, ztopmax_patch, xpatches*ypatches, MPI_MAX, comm3d,mpierr)
      ztop_patch = patchsum_1level_field_r(ztop_field)
    endif

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

    if (lhetero) then
      do j=1,ypatches
         do i=1,xpatches
           if (cc_patch(i,j) > 0.0) then
             zbase_patch(i,j) = zbase_patch(i,j)/cc_patch(i,j)
             ztop_patch (i,j) = ztop_patch(i,j) /cc_patch(i,j)
           else
             zbase_patch(i,j) = 0.0
             ztop_patch (i,j) = 0.0
           endif
           cc_patch    = cc_patch    * (xpatches*ypatches/ijtot)
           qlint_patch = qlint_patch * (xpatches*ypatches/ijtot)
        enddo
      enddo
    endif

    call D_MPI_ALLREDUCE(tke_totl, tke_tot, 1,   &
                          MPI_SUM, comm3d,mpierr)

    tke_tot = tke_tot / ijtot

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

    if (lhetero) then
      ust_patch = patchsum_1level(ustar(2:i1,2:j1)) * (xpatches*ypatches/ijtot)
      tst_patch = patchsum_1level(- thlflux(2:i1,2:j1) / ustar(2:i1,2:j1)) * (xpatches*ypatches/ijtot)
      qst_patch = patchsum_1level(-  qtflux(2:i1,2:j1) / ustar(2:i1,2:j1)) * (xpatches*ypatches/ijtot)
    endif

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

    if (lhetero) then
      if(isurf < 3) then
        wthls_patch  = patchsum_1level(thlflux(2:i1, 2:j1)) * (xpatches*ypatches/ijtot)
        wqls_patch = patchsum_1level( qtflux(2:i1, 2:j1)) * (xpatches*ypatches/ijtot)
      else
        wthls_patch  = wt_patch
        wqls_patch = wq_patch
      endif
      wthvs_patch = (1.+(rv/rd-1)*qts_patch) * wthls_patch + c2 * (thls_patch) * wq_patch
      obl_patch  = patchsum_1level(obl(2:i1, 2:j1)) * (xpatches*ypatches/ijtot)
    endif

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

      if (lhetero) then
        Qnet_patch     = patchsum_1level(Qnet    (2:i1, 2:j1)) * (xpatches*ypatches/ijtot)
        H_patch        = patchsum_1level(H       (2:i1, 2:j1)) * (xpatches*ypatches/ijtot)
        LE_patch       = patchsum_1level(LE      (2:i1, 2:j1)) * (xpatches*ypatches/ijtot)
        G0_patch       = patchsum_1level(G0      (2:i1, 2:j1)) * (xpatches*ypatches/ijtot)
        tendskin_patch = patchsum_1level(tendskin(2:i1, 2:j1)) * (xpatches*ypatches/ijtot)
        rs_patch       = patchsum_1level(rs      (2:i1, 2:j1)) * (xpatches*ypatches/ijtot)
        ra_patch       = patchsum_1level(ra      (2:i1, 2:j1)) * (xpatches*ypatches/ijtot)
        cliq_patch     = patchsum_1level(cliq    (2:i1, 2:j1)) * (xpatches*ypatches/ijtot)
        wl_patch       = patchsum_1level(wl      (2:i1, 2:j1)) * (xpatches*ypatches/ijtot)
        rsveg_patch    = patchsum_1level(rsveg   (2:i1, 2:j1)) * (xpatches*ypatches/ijtot)
        rssoil_patch   = patchsum_1level(rssoil  (2:i1, 2:j1)) * (xpatches*ypatches/ijtot)
        tskin_patch    = patchsum_1level(tskin   (2:i1, 2:j1)) * (xpatches*ypatches/ijtot)
      endif
    else if (isurf == 11) then
      Qnet(:,:) = swd(i,j,1) + swu(i,j,1) + lwd(i,j,1) + lwu(i,j,1)

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
        !skip for ws
        if (trim(tile(ilu)%lushort) == 'ws') cycle 
        obuk_av(ilu)  = mean_2d(tile(ilu)%obuk)
        ustar_av(ilu) = mean_2d(tile(ilu)%ustar)
        ra_av(ilu)    = mean_2d(tile(ilu)%ra)
        f2_av(ilu)    = mean_2d(tile(ilu)%f2)
        f3_av(ilu)    = mean_2d(tile(ilu)%f3)
      end do

      do ilu=1,nlu
        !skip for ws and aq
        if (trim(tile(ilu)%lushort) == 'ws' .or. &
            trim(tile(ilu)%lushort) == 'aq') cycle 
        rs_av(ilu)    = mean_2d(tile(ilu)%rs)
      end do

      f1_av    = mean_2d(f1)
      f2b_av   = mean_2d(f2b)

      do ilu=1,nlu
        c_av(ilu)       = mean_2d(tile(ilu)%frac)
        H_av(ilu)       = mean_2d(tile(ilu)%H)
        LE_av(ilu)      = mean_2d(tile(ilu)%LE)
        thlskin_av(ilu) = mean_2d(tile(ilu)%thlskin) 
        qtskin_av(ilu)  = mean_2d(tile(ilu)%qtskin) 
      end do

      wlav = mean_2d(wl)

      do ilu=1,nlu
        !skip for aq
        if (trim(tile(ilu)%lushort) == 'aq') cycle
        G_av(ilu)    = mean_2d(tile(ilu)%G)
      end do

      if (lags) then
        an_co2_av   = mean_2d(an_co2)
        resp_co2_av = mean_2d(resp_co2)
      endif
     end if

    ! calculate radiation fluxes at surface, TOM, TOA
    if (iradiation /= 0) then
       vars(ivar_rad+ 0) = abs(sum(lwd(2:i1,2:j1,1)))  !'rlds',   'surface downwelling longwave flux'
       vars(ivar_rad+ 1) = abs(sum(lwu(2:i1,2:j1,1)))  !'rlus',   'surface upwelling longwave flux'
       vars(ivar_rad+ 2) = abs(sum(swd(2:i1,2:j1,1)))  !'rsds',   'surface downwelling shortwave flux'
       vars(ivar_rad+ 3) = abs(sum(swu(2:i1,2:j1,1)))  !'rsus',   'surface upwelling shortwave flux'

       vars(ivar_rad+ 4) = abs(sum(swdca(2:i1,2:j1,1)))  !'rsdscs', 'surface downwelling shortwave flux - clear sky'
       vars(ivar_rad+ 5) = abs(sum(swuca(2:i1,2:j1,1)))  !'rsuscs', 'surface upwelling shortwave flux - clear sky'
       vars(ivar_rad+ 6) = abs(sum(lwdca(2:i1,2:j1,1)))  !'rldscs', 'surface downwelling longwave flux - clear sky'
       vars(ivar_rad+ 7) = abs(sum(lwuca(2:i1,2:j1,1)))  !'rluscs', 'surface upwelling longwave flux - clear sky'

       vars(ivar_rad+ 8) = abs(sum(SW_dn_TOA(2:i1,2:j1)))    !'rsdt',   'TOA incoming shortwave flux','W/m^2'
       vars(ivar_rad+ 9) = abs(sum(SW_up_TOA(2:i1,2:j1)))    !'rsut',   'TOA outgoing shortwave flux','W/m^2'
       vars(ivar_rad+10) = abs(sum(LW_up_TOA(2:i1,2:j1)))    !'rlut',   'TOA outgoing longwave flux','W/m^2'
       vars(ivar_rad+11) = abs(sum(SW_up_ca_TOA(2:i1,2:j1))) !'rsutcs', 'TOA outgoing shortwave flux -clear sky'
       vars(ivar_rad+12) = abs(sum(LW_up_ca_TOA(2:i1,2:j1))) !'rlutcs', 'TOA outgoing longwave flux -clear sky'

       vars(ivar_rad+13) = abs(sum(swd(2:i1,2:j1,kmax)))   !'rsdtm',  'TOM incoming shortwave flux'
       vars(ivar_rad+14) = abs(sum(lwd(2:i1,2:j1,kmax)))   !'rsdtm',  'TOM incoming longwave flux'
       vars(ivar_rad+15) = abs(sum(swu(2:i1,2:j1,kmax)))   !'rsutm',  'TOM outgoing shortwave flux'
       vars(ivar_rad+16) = abs(sum(lwu(2:i1,2:j1,kmax)))   !'rlutm',  'TOM outgoing longwave flux'
       vars(ivar_rad+17) = abs(sum(swuca(2:i1,2:j1,kmax))) !'rsutmcs','TOM outgoing shortwave flux -clear sky'
       vars(ivar_rad+18) = abs(sum(lwuca(2:i1,2:j1,kmax))) !'rlutmcs','TOM outgoing longwave flux -clear sky'

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

      if(lhetero) then
        do i=1,xpatches
          do j=1,ypatches
            name = 'tmser1patchiiixjjj.'//cexpnr
            write (name(12:14),'(i3.3)') i
            write (name(16:18),'(i3.3)') j
            open (ifoutput,file=name,position='append')
            write( ifoutput,'(f10.2,f6.3,4f12.3,f10.4,5f9.3)') &
              rtimee, &
              cc_patch(i,j), &
              zbase_patch(i,j), &
              ztop_patch(i,j), &
              ztopmax_patch(i,j), &
              zi_patch(i,j), &
              we_patch(i,j), &
              qlint_patch(i,j)*1000., &
              qlintmax_patch(i,j)*1000., &
              wmax_patch(i,j), &
              tke_tot_patch(i,j), &
              qlmax_patch(i,j)*1000.
            close(ifoutput)

            name = 'tmsurfpatchiiixjjj.'//cexpnr
            write (name(12:14),'(i3.3)') i
            write (name(16:18),'(i3.3)') j
            open (ifoutput,file=name,position='append')
            write( ifoutput,'(f10.2,4e11.3,f11.3,4e11.3)') &
              rtimee      ,&
              ust_patch(i,j)   ,&
              tst_patch(i,j)   ,&
              qst_patch(i,j)   ,&
              obl_patch(i,j)   ,&
              thls_patch(i,j)  ,&
              z0mav_patch(i,j) ,&
              wthls_patch(i,j)   ,&
              wthvs_patch(i,j)  ,&
              wqls_patch(i,j)
            close(ifoutput)

            if(isurf == 1) then
              name = 'tmlsmpatchiiixjjj.'//cexpnr
              write (name(11:13),'(i3.3)') i
              write (name(15:17),'(i3.3)') j
              open (ifoutput,file=name,position='append')
              write(ifoutput,'(f10.2,9f11.3,e13.3, 2f11.3)') &
                rtimee           ,&
                Qnet_patch(i,j)       ,&
                H_patch(i,j)          ,&
                LE_patch(i,j)         ,&
                G0_patch(i,j)         ,&
                tendskin_patch(i,j)   ,&
                rs_patch(i,j)         ,&
                ra_patch(i,j)         ,&
                tskin_patch(i,j)      ,&
                cliq_patch(i,j)       ,&
                wl_patch(i,j)         ,&
                rssoil_patch(i,j)     ,&
                rsveg_patch(i,j)
              close(ifoutput)
            endif
          enddo
        enddo
      endif

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

    use modglobal,  only : ih,i1,jh,j1,kmax,k1,cp,rlv,imax,rd,zh,dzh,zf,dzf,rv,ijtot,iadv_sv,iadv_kappa
    use modfields,  only : w0,qt0,qt0h,ql0,thl0,thl0h,thv0h,sv0,exnf,whls
    use modsurfdata,only : svs, lhetero, xpatches, ypatches
    use modsurface, only : patchxnr,patchynr
    use modmpi,     only : mpierr, comm3d,mpi_sum, D_MPI_ALLREDUCE
    use advec_kappa,only : halflev_kappa

    implicit none

    real    :: zil, dhdt, locval, oldlocval
    integer :: location, i, j, k, nsamp, stride
    real, allocatable,dimension(:,:,:) :: blh_fld2


    if (lhetero) then
      allocate(blh_fld2(2:i1,2:j1,k1))
    endif

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

    if (lhetero) then !Not using stride, instead use adjacent grid points in x-direction (to prevent processor communication)

      select case (iblh_meth)
        case (iblh_flux)
          blh_fld2(2,2:j1,:)        = blh_fld(i1,2:j1,:)       + blh_fld(2,2:j1,:)        + blh_fld(3,2:j1,:)
          blh_fld2(3:(i1-1),2:j1,:) = blh_fld(2:(i1-2),2:j1,:) + blh_fld(3:(i1-1),2:j1,:) + blh_fld(4:i1,2:j1,:)
          blh_fld2(i1,2:j1,:)       = blh_fld(i1-1,2:j1,:)     + blh_fld(i1,2:j1,:)       + blh_fld(2,2:j1,:)

          do i=2,i1
            do j=2,j1
              zi_field(i,j) = zh(minloc(blh_fld2(i,j,:),1))
            end do
          end do

        case (iblh_grad)
          blh_fld2(2,2:j1,:)        = blh_fld(i1,2:j1,:)       + blh_fld(2,2:j1,:)        + blh_fld(3,2:j1,:)
          blh_fld2(3:(i1-1),2:j1,:) = blh_fld(2:(i1-2),2:j1,:) + blh_fld(3:(i1-1),2:j1,:) + blh_fld(4:i1,2:j1,:)
          blh_fld2(i1,2:j1,:)       = blh_fld(i1-1,2:j1,:)     + blh_fld(i1,2:j1,:)       + blh_fld(2,2:j1,:)

          do i=2,i1
            do j=2,j1
              profile  = blh_fld2(i,j,:)
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
              zi_field(i,j) = (zh(location-1) - dzh(location)*dgrad(location-1)/(dgrad(location)-dgrad(location-1) + 1.e-8))
            enddo
          enddo

        case (iblh_thres)
          blh_fld2(2,2:j1,:)        = blh_fld(i1,2:j1,:)       + blh_fld(2,2:j1,:)        + blh_fld(3,2:j1,:)
          blh_fld2(3:(i1-1),2:j1,:) = blh_fld(2:(i1-2),2:j1,:) + blh_fld(3:(i1-1),2:j1,:) + blh_fld(4:i1,2:j1,:)
          blh_fld2(i1,2:j1,:)       = blh_fld(i1-1,2:j1,:)     + blh_fld(i1,2:j1,:)       + blh_fld(2,2:j1,:)

          do i=2,i1
            do j=2,j1
              locval = 0.0
              do k=kmax,1,-1
                oldlocval = locval
                locval = blh_sign*blh_fld2(i,j,k)/3
                if (locval < blh_sign*blh_thres) then
                  zi_field(i,j) = (zf(k) +  (blh_sign*blh_thres-locval) * dzh(k+1)/(oldlocval-locval))
                  exit
                endif
              enddo
            enddo
          enddo

      end select

    endif

    call D_MPI_ALLREDUCE(zil, zi, 1, MPI_SUM, comm3d,mpierr)
    zi = zi / ijtot

    if (lhetero) then
      zi_patch = patchsum_1level_field_r(zi_field) * (xpatches*ypatches/ijtot)
    endif

    if (ziold< 0) ziold = zi
    dhdt = (zi-ziold)/dtav
    if(store_zi) ziold = zi

    k=2
    do while (zh(k)<zi .and. k < kmax)
      k=k+1
    end do
    we = dhdt - whls (k)   !include for large-scale vertical velocity

    if (lhetero) then
      do j=1,ypatches
        do i=1,xpatches
          if (ziold_patch(i,j)<0) ziold_patch(i,j) = zi_patch(i,j)
          k=2
          do while (zh(k)<zi_patch(i,j) .and. k < kmax)
            k=k+1
          end do
          we_patch(i,j) = ((zi_patch(i,j)-ziold_patch(i,j))/dtav) - whls(k)
        enddo
      enddo
      if(store_zi) ziold_patch = zi_patch
    endif

    if (lhetero) then
      deallocate(blh_fld2)
    endif

    !$acc wait

  end subroutine calcblheight

!> Clean up when leaving the run
  subroutine exittimestat
    use modmpi, only : myid
    use modstat_nc, only : exitstat_nc,lnetcdf
    use modsurfdata,only :lhetero
    implicit none

    if(ltimestat .and. lnetcdf .and. myid==0) call exitstat_nc(ncid)
    if(.not.ltimestat) return

    !$acc exit data delete(blh_fld, sv0h, profile, gradient, dgrad)

    deallocate(blh_fld,sv0h)
    deallocate(profile,gradient,dgrad)

    if (lhetero) then
      deallocate(zbase_field  )
      deallocate(ztop_field   )
      deallocate(cc_field     )
      deallocate(qlint_field  )
      deallocate(tke_tot_field)

      deallocate(zbase_patch    )
      deallocate(ztop_patch     )
      deallocate(zbasemin_patch )
      deallocate(zbasemin_patchl)
      deallocate(cc_patch       )
      deallocate(qlint_patch    )
      deallocate(qlintmax_patch )
      deallocate(qlintmax_patchl)
      deallocate(tke_tot_patch  )
      deallocate(wmax_patch     )
      deallocate(wmax_patchl    )
      deallocate(qlmax_patch    )
      deallocate(qlmax_patchl   )
      deallocate(ztopmax_patch  )
      deallocate(ztopmax_patchl )

      deallocate(u0av_patch)
      deallocate(v0av_patch)
      deallocate(w0av_patch)

      deallocate(ust_patch )
      deallocate(qst_patch )
      deallocate(tst_patch )
      deallocate(wthls_patch )
      deallocate(wqls_patch)
      deallocate(wthvs_patch)

      deallocate(Qnet_patch    )
      deallocate(H_patch       )
      deallocate(LE_patch      )
      deallocate(G0_patch      )
      deallocate(tendskin_patch)
      deallocate(rs_patch      )
      deallocate(ra_patch      )
      deallocate(cliq_patch    )
      deallocate(wl_patch      )
      deallocate(rsveg_patch   )
      deallocate(rssoil_patch  )
      deallocate(tskin_patch   )
      deallocate(obl_patch     )

      deallocate(zi_patch      )
      deallocate(zi_field      )
      deallocate(ziold_patch   )
      deallocate(we_patch      )
    endif
  end subroutine exittimestat

 function patchsum_1level(x)
   use modglobal,  only : imax,jmax
   use modsurface, only : patchxnr, patchynr
   use modsurfdata,only : xpatches,ypatches
   use modmpi,     only : mpierr,comm3d,mpi_sum, D_MPI_ALLREDUCE
   implicit none
   real                :: patchsum_1level(xpatches,ypatches),xl(xpatches,ypatches)
   real, intent(in)    :: x(imax,jmax)
   integer             :: i,j,iind,jind

   patchsum_1level = 0
   xl              = 0

   do j=1,jmax
     jind  = patchynr(j)

     do i=1,imax
       iind  = patchxnr(i)

       xl(iind,jind) = xl(iind,jind) + x(i,j)
     enddo
   enddo

  call D_MPI_ALLREDUCE(xl,patchsum_1level, xpatches*ypatches,MPI_SUM, comm3d,mpierr)

  end function

!> It might be the same as the normal one, it might have a different kind ...
 function patchsum_1level_field_r(x)
   use modglobal,  only : imax,jmax
   use modsurface, only : patchxnr, patchynr
   use modsurfdata,only : xpatches,ypatches
   use modmpi,     only : mpierr,comm3d,mpi_sum, D_MPI_ALLREDUCE
   implicit none
   real                :: patchsum_1level_field_r(xpatches,ypatches),xl(xpatches,ypatches)
   real(field_r), intent(in) :: x(imax,jmax)
   integer             :: i,j,iind,jind

   patchsum_1level_field_r = 0
   xl              = 0

   do j=1,jmax
     jind  = patchynr(j)

     do i=1,imax
       iind  = patchxnr(i)

       xl(iind,jind) = xl(iind,jind) + x(i,j)
     enddo
   enddo

  call D_MPI_ALLREDUCE(xl,patchsum_1level_field_r, xpatches*ypatches,MPI_SUM, comm3d,mpierr)

  end function

end module modtimestat
