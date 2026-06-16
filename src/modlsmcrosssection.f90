!> \file modlsmcrosssection.f90
!!   Dumps an instantenous lsmcrosssection of the field
!! lsmcrosssections in the yz-plane and in the xy-plane            |
!        of u,v,w,thl,thv,qt,ql. Written to movv_*.expnr and movh_*.expnr
!! If netcdf is true, this module leads the cross.myid.expnr.nc output

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
module modlsmcrosssection


  use modglobal, only : longint
  use modsurfdata,only : ksoilmax
  use modlogging, only: finish

implicit none
character(len=*), parameter :: modname = 'modlsmcrosssection'
private
PUBLIC :: initlsmcrosssection, lsmcrosssection,exitlsmcrosssection
save
!NetCDF variables
  integer,parameter :: nvar=2
  integer :: nvar3
  integer :: nvar4
  integer :: nvars3d4
  integer :: ncid1 = 0
  integer :: ncid2 = 0
  integer :: ncid3 = 0
  integer :: ncid4 = 0
  integer :: nrec1 = 0
  integer :: nrec2 = 0
  integer :: nrec3 = 0
  integer :: nrec4 = 0
  integer :: crossheight
  character(4) :: cheight
  character(80) :: fname1 = 'lsmcrossxz.xxxxyxxx.xxx.nc'
  character(80) :: fname2 = 'lsmcrossxy.xxxx.xxxxyxxx.xxx.nc'
  character(80) :: fname3 = 'surfcross.xxxxyxxx.xxx.nc'
  character(80) :: fname4 = 'slrbcross.xxxxyxxx.xxx.nc'
  character(80), dimension(nvar,4) :: ncname1
  character(80), dimension(1,4) :: tncname1
  character(80), dimension(nvar,4) :: ncname2
  character(80), dimension(1,4) :: tncname2
  character(80), allocatable, dimension(:,:) :: ncname3
  character(80), allocatable, dimension(:,:) :: ncname4
  character(80), dimension(1,4) :: tncname3
  character(80), dimension(1,4) :: tncname4

  real    :: dtav
  integer(kind=longint) :: idtav, tnext
  logical :: lcross = .false.     !< switch for doing the lsmcrosssection (on/off)
  logical :: lcrosssoil = .false. !< switch for doing vertical soil crosssection (on/off)
  integer :: crossplane = 2       !< Location of the xz lsmcrosssection

contains
!> Initializing lsmcrosssection. Read out the namelist, initializing the variables
  subroutine initlsmcrosssection
    use modmpi,      only : myid,myidy,mpierr,comm3d,cmyid,D_MPI_BCAST
    use modglobal,   only : imax,jmax,ifnamopt,fname_options,dtmax,dtav_glob,ladaptive,j1,dt_lim,cexpnr,tres,btime,checknamelisterror,&
                            output_prefix
    use modstat_nc,  only : lnetcdf,open_nc, define_nc,ncinfo,nctiminfo,writestat_dims_nc
    use modsurfdata, only : isurf
    use modlsm,      only : lags
    use modslurb,    only : enable_slurb
    use fortran_support, only: nnml_output
    implicit none

    character(len=*), parameter :: routine = modname//'/initlsmcrossection'

    integer :: ierr

    namelist/NAMLSMCROSSSECTION/ &
    lcross, lcrosssoil, dtav, crossheight, crossplane

    crossheight=2
    ncid2=2
    nrec2=0

    dtav = dtav_glob
    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMLSMCROSSSECTION,iostat=ierr)
       call checknamelisterror(ierr, ifnamopt, 'NAMLSMCROSSSECTION')
      write(nnml_output ,NAMLSMCROSSSECTION)
      close(ifnamopt)
    end if

    if (lcross .and. .not. (isurf == 1 .or. isurf == 2 .or. isurf == 11)) then
       lcross = .FALSE.
       write (6,*) "Ignoring lcross, lsmcrossection currently implemented only for isurf==1, 2, or 11."
    endif

    if (lcrosssoil .and. .not. (isurf == 1 .or. isurf == 11)) then
       lcrosssoil = .FALSE.
       write (6,*) "Ignoring lcrosssoil, lsm soil crossection currently implemented only for isurf==1 or 11."
    endif

    call D_MPI_BCAST(dtav       ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(lcross     ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(lcrosssoil ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(crossheight,1,0,comm3d,mpierr)
    call D_MPI_BCAST(crossplane ,1,0,comm3d,mpierr)

    idtav = int(dtav / tres, kind=kind(idtav))
    tnext   = idtav+btime
    if(.not.(lcross .or. lcrosssoil)) return
    dt_lim = min(dt_lim,tnext)

    if (lcrosssoil) then
    if((crossheight>ksoilmax) .or. crossplane>j1) then
      call finish(routine, 'lsmcrosssection: lsmcrosssection out of range')
    end if
    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      call finish(routine, 'lsmcrosssection: dtav should be a integer multiple of dtmax')
    end if
    if (lnetcdf) then
      if (myidy==0) then
        fname1(12:19) = cmyid
        fname1(21:23) = cexpnr
        call nctiminfo(tncname1(1,:))
        call ncinfo(ncname1( 1,:),'tsoil', 'xz crosssection of the Soil temperature','K','t0tts')
        call ncinfo(ncname1( 2,:),'phiw', 'xz crosssection of the Soil moisture','m3/m3','t0tts')
        call open_nc(trim(output_prefix)//fname1,  ncid1,nrec1,n1=imax,ns=ksoilmax)
        if (nrec1==0) then
          call define_nc( ncid1, 1, tncname1)
          call writestat_dims_nc(ncid1)
          call define_nc( ncid1, NVar, ncname1)
        end if
      end if
        write(cheight,'(i4.4)') crossheight
        fname2(12:15) = cheight
        fname2(17:24) = cmyid
        fname2(26:28) = cexpnr
        call nctiminfo(tncname2(1,:))
        call ncinfo(ncname2( 1,:),'tsoil', 'xy crosssection of the Soil temperature','K','tt0t')
        call ncinfo(ncname2( 2,:),'phiw', 'xy crosssection of the Soil moisture','m3/m3','tt0t')
        call open_nc(trim(output_prefix)//fname2,  ncid2,nrec2,n1=imax,n2=jmax)
        if (nrec2==0) then
          call define_nc( ncid2, 1, tncname2)
          call writestat_dims_nc(ncid2)
          call define_nc( ncid2, NVar, ncname2)
        end if
     end if
  end if


  if (lnetcdf .and. lcross) then
!       ! Surface values
        fname3(11:18) = cmyid
        fname3(20:22) = cexpnr
        if (isurf == 1) then
            nvar3 = 12
         else if (isurf == 2) then
            nvar3 = 8
        else if (isurf == 11) then
            nvar3 = 14
            if (lags) nvar3 = nvar3 + 2
        end if

        allocate(ncname3(nvar3,4))

        if (isurf == 1) then
            call nctiminfo(tncname3(1,:))
            call ncinfo(ncname3( 1,:),'Qnet','Net radiation','W/m^2','tt0t')
            call ncinfo(ncname3( 2,:),'H','Sensible heat flux','W/m^2','tt0t')
            call ncinfo(ncname3( 3,:),'LE','Latent heat flux','W/m^2','tt0t')
            call ncinfo(ncname3( 4,:),'G0','Ground heat flux','W/m^2','tt0t')
            call ncinfo(ncname3( 5,:),'tskin','Skin temperature','K','tt0t')
            call ncinfo(ncname3( 6,:),'tendskin','Skin tendency','W/m^2','tt0t')
            call ncinfo(ncname3( 7,:),'rs','Surface resistance','s/m','tt0t')
            call ncinfo(ncname3( 8,:),'ra','Aerodynamic resistance','s/m','tt0t')
            call ncinfo(ncname3( 9,:),'cliq','Fraction of vegetated surface covered with liquid water','-','tt0t')
            call ncinfo(ncname3(10,:),'Wl','Liquid water reservoir','m','tt0t')
            call ncinfo(ncname3(11,:),'rssoil','Soil evaporation resistance','s/m','tt0t')
            call ncinfo(ncname3(12,:),'rsveg','Vegetation resistance','s/m','tt0t')
            call open_nc(trim(output_prefix)//fname3,  ncid3,nrec3,n1=imax,n2=jmax)
            if (nrec3==0) then
              call define_nc(ncid3, 1, tncname3)
              call writestat_dims_nc(ncid3)
            end if
            call define_nc(ncid3, nvar3, ncname3)
        else if (isurf == 2) then
            call nctiminfo(tncname3(1,:))
            call ncinfo(ncname3( 1,:),'hfss','Surface upward sensible heat flux','W/m^2','tt0t')
            call ncinfo(ncname3( 2,:),'hfls','Surface upward latent heat flux','W/m^2','tt0t')
            call ncinfo(ncname3( 3,:),'obuk', 'Obukhov length', 'm', 'tt0t')
            call ncinfo(ncname3( 4,:),'ustar', 'Friction velocity', 'm/s^-1', 'tt0t')
            call ncinfo(ncname3( 5,:),'Cs', 'Drag coefficient for scalars', '-', 'tt0t')
            call ncinfo(ncname3( 6,:),'Cm', 'Drag coefficient for momentum','-', 'tt0t')
            call ncinfo(ncname3( 7,:),'z0h','Surface roughness length for heat','m', 'tt0t')
            call ncinfo(ncname3( 8,:),'z0m','Surface roughness length for momentum','m', 'tt0t')
            call open_nc(trim(output_prefix)//fname3,  ncid3,nrec3,n1=imax,n2=jmax)
            if (nrec3==0) then
               call define_nc(ncid3, 1, tncname3)
               call writestat_dims_nc(ncid3)
            end if
            call define_nc(ncid3, nvar3, ncname3)
        else if (isurf == 11) then
            call nctiminfo(tncname3(1,:))
            call ncinfo(ncname3( 1,:),'H', 'Sensible heat flux', 'W/m^2', 'tt0t')
            call ncinfo(ncname3( 2,:),'LE', 'Latent heat flux', 'W/m^2', 'tt0t')
            call ncinfo(ncname3( 3,:),'G0', 'Ground heat flux', 'W/m^2', 'tt0t')
            call ncinfo(ncname3( 4,:),'tskin', 'Skin temperature', 'K', 'tt0t')
            call ncinfo(ncname3( 5,:),'obuk', 'Obukhov length', 'm', 'tt0t')
            call ncinfo(ncname3( 6,:),'ustar', 'Friction velocity', 'm/s^-1', 'tt0t')
            call ncinfo(ncname3( 7,:),'cliq', 'Fraction of vegetated surface covered with liquid water', '-', 'tt0t')
            call ncinfo(ncname3( 8,:),'wl', 'Liquid water reservoir', 'm', 'tt0t')
            call ncinfo(ncname3( 9,:),'ra', 'Aerodynamic resistance', 's/m', 'tt0t')
            call ncinfo(ncname3(10,:),'rssoil', 'Soil evaporation resistance', 's/m', 'tt0t')
            call ncinfo(ncname3(11,:),'rsveg', 'Vegetation resistance', 's/m', 'tt0t')
            call ncinfo(ncname3(12,:),'f1', 'f1(SWD) function vegetation resistance', 's/m', 'tt0t')
            call ncinfo(ncname3(13,:),'f2_b', 'f2(theta) function soil resistance', 's/m', 'tt0t')
            call ncinfo(ncname3(14,:),'Qnet', 'Net radiation', 'W/m^2', 'tt0t')
            if (lags) then
              call ncinfo(ncname3(15,:),'an_co2', 'Net CO2 assimilation', 'ppm m s-1', 'tt0t')
              call ncinfo(ncname3(16,:),'resp_co2', 'CO2 respiration soil + plant', 'ppm m s-1', 'tt0t')
            end if

            call open_nc(trim(output_prefix)//fname3,  ncid3,nrec3,n1=imax,n2=jmax)
            if (nrec3==0) then
              call define_nc(ncid3, 1, tncname3)
              call writestat_dims_nc(ncid3)
            end if
            call define_nc( ncid3, nvar3, ncname3)
        end if

        if (enable_slurb) then
        !
        ! Surface values
        fname4(11:18) = cmyid
        fname4(20:22) = cexpnr
        
        nvar4 = 145
        nvars3d4 = 21
        nvar4 = nvar4 + nvars3d4
        nrec4=0

        allocate(ncname4(nvar4,4))

        call nctiminfo(tncname4(1,:))
        call ncinfo(ncname4( 1,:),'albedo_urb', 'effective urban albedo', '', 'tt0t')
        call ncinfo(ncname4( 2,:),'emiss_urb', 'effective urban emissivity', '', 'tt0t')
        call ncinfo(ncname4( 3,:),'ol_urb', 'urban Obukhov length', '', 'tt0t')
        call ncinfo(ncname4( 4,:),'qsws_urb', 'total urban latent heat flux', '', 'tt0t')
        call ncinfo(ncname4( 5,:),'rad_lw_in_urb', 'incoming longwave radiation', '', 'tt0t')
        call ncinfo(ncname4( 6,:),'rad_lw_out_urb', 'outgoing longwave radiation', '', 'tt0t')
        call ncinfo(ncname4( 7,:),'rad_sw_in_urb', 'incoming shortwave radiation', '', 'tt0t')
        call ncinfo(ncname4( 8,:),'rad_sw_out_urb', 'outgoing shortwave radiation', '', 'tt0t')
        call ncinfo(ncname4( 9,:),'ram_urb', 'urban aerodynamic resistance for momentum', '', 'tt0t')
        call ncinfo(ncname4( 10,:),'rib_urb', 'urban bulk-Richardson number', '', 'tt0t')
        call ncinfo(ncname4( 11,:),'shf_urb', 'total urban sensible heat flux', '', 'tt0t')
        call ncinfo(ncname4( 12,:),'t_2m_urb', 'urban 2-metre temperature ', 'extrapolated) (K', 'tt0t')
        call ncinfo(ncname4( 13,:),'t_c_urb', 'complete ', 'area-weighted) urban surface temperature (K', 'tt0t')
        call ncinfo(ncname4( 14,:),'t_h_urb', 'effective urban surface temperature ', 'K', 'tt0t')
        call ncinfo(ncname4( 15,:),'thl_rad_urb', 'urban radiative surface liquid water potential temperature ', 'K', 'tt0t')
        call ncinfo(ncname4( 16,:),'usws_urb', 'urban momentum flux ', 'u-component', 'tt0t')
        call ncinfo(ncname4( 17,:),'vsws_urb', 'urban momentum flux ', 'v-component', 'tt0t')
        call ncinfo(ncname4( 18,:),'thlskin', 'urban skin liquid water potential temperature', '', 'tt0t')
        call ncinfo(ncname4( 19,:),'qtskin', 'urban skin specific humidity TODOSELF', '', 'tt0t')
        call ncinfo(ncname4( 20,:),'m_liq_road_0', 'liquid water reservoir on roads', '', 'tt0t')
        call ncinfo(ncname4( 21,:),'m_liq_road_m', 'prev. liquid water reservoir on roads', '', 'tt0t')
        call ncinfo(ncname4( 22,:),'m_liq_roof_0', 'liquid water reservoir on roofs', '', 'tt0t')
        call ncinfo(ncname4( 23,:),'m_liq_roof_m', 'prev. liquid water reservoir on roofs', '', 'tt0t')
        call ncinfo(ncname4( 24,:),'q_can_0', 'canyon mixing ratio ', 'kg/kg', 'tt0t')
        call ncinfo(ncname4( 25,:),'q_can_m', 'previous canyon mixing ratio ', 'kg/kg', 'tt0t')
        call ncinfo(ncname4( 26,:),'t_can_0', 'canyon air temperature ', 'K', 'tt0t')
        call ncinfo(ncname4( 27,:),'t_can_m', 'prev. canyon temperature ', 'K', 'tt0t')
        call ncinfo(ncname4( 28,:),'tm_liq_road', 'road liquid water reservoir tendency', '', 'tt0t')
        call ncinfo(ncname4( 29,:),'tm_liq_roof', 'roof liquid water reservoir tendency', '', 'tt0t')
        call ncinfo(ncname4( 30,:),'tq_can', 'canyon mixing ratio tendency ', 'kg/kg/s', 'tt0t')
        call ncinfo(ncname4( 31,:),'tt_can', 'canyon temperature tendency ', 'K/s', 'tt0t')
        call ncinfo(ncname4( 32,:),'pt_road', 'road surface potential temperature', '', 'tt0t')
        call ncinfo(ncname4( 33,:),'pt_roof', 'roof surface potential temperature', '', 'tt0t')
        call ncinfo(ncname4( 34,:),'pt_wall_a', 'wall A surface potential temperature', '', 'tt0t')
        call ncinfo(ncname4( 35,:),'pt_wall_b', 'wall B surface potential temperature', '', 'tt0t')
        call ncinfo(ncname4( 36,:),'pt_win_a', 'window A surface potential temperature', '', 'tt0t')
        call ncinfo(ncname4( 37,:),'pt_win_b', 'window A surface potential temperature', '', 'tt0t')
        call ncinfo(ncname4( 38,:),'q_road', 'road surface mixing ratio', '', 'tt0t')
        call ncinfo(ncname4( 39,:),'q_roof', 'roof surface mixing ratio', '', 'tt0t')
        call ncinfo(ncname4( 40,:),'qs_road', 'road surface saturation mixing ratio', '', 'tt0t')
        call ncinfo(ncname4( 41,:),'qs_roof', 'roof surface saturation mixing ratio', '', 'tt0t')
        call ncinfo(ncname4( 42,:),'vpt_road', 'road surface virtual potential temperature', '', 'tt0t')
        call ncinfo(ncname4( 43,:),'vpt_roof', 'roof surface virtual potential temperature', '', 'tt0t')
        call ncinfo(ncname4( 43,:),'shf_can', 'sensible heat flux between the street canyon and the atmosphere', '', 'tt0t')
        call ncinfo(ncname4( 44,:),'shf_external', 'sensible heat flux external to the model ', 'e.g. industry', 'tt0t')
        call ncinfo(ncname4( 45,:),'shf_road', 'road surface sensible heat flux', '', 'tt0t')
        call ncinfo(ncname4( 46,:),'shf_roof', 'roof surface sensible heat flux', '', 'tt0t')
        call ncinfo(ncname4( 47,:),'shf_traffic', 'traffic sensible heat flux ', 'input-only', 'tt0t')
        call ncinfo(ncname4( 48,:),'shf_wall_a', 'wall A sensible heat flux', '', 'tt0t')
        call ncinfo(ncname4( 49,:),'shf_wall_b', 'wall B sensible heat flux', '', 'tt0t')
        call ncinfo(ncname4( 50,:),'shf_win_a', 'window A sensible heat flux', '', 'tt0t')
        call ncinfo(ncname4( 51,:),'shf_win_b', 'window B sensible heat flux', '', 'tt0t')
        call ncinfo(ncname4( 52,:),'qsws_can', 'latent heat flux between the street canyon and the atmosphere', '', 'tt0t')
        call ncinfo(ncname4( 53,:),'qsws_external', 'latent heat flux external to the model ', 'e.g. industry', 'tt0t')
        call ncinfo(ncname4( 54,:),'qsws_liq_road', 'roof latent heat flux ', 'liquid incl. precipitation', 'tt0t')
        call ncinfo(ncname4( 55,:),'qsws_liq_roof', 'roof latent heat flux ', 'liquid incl. precipitation', 'tt0t')
        call ncinfo(ncname4( 56,:),'qsws_road', 'road latent heat flux', '', 'tt0t')
        call ncinfo(ncname4( 57,:),'qsws_roof', 'roof latent heat flux', '', 'tt0t')
        call ncinfo(ncname4( 58,:),'c_liq_road', 'liquid water coverage on road', '', 'tt0t')
        call ncinfo(ncname4( 59,:),'c_liq_roof', 'liquid water coverage on roof', '', 'tt0t')
        call ncinfo(ncname4( 60,:),'ghf_road', 'road ground heat flux', '', 'tt0t')
        call ncinfo(ncname4( 61,:),'ghf_roof', 'roof indoor heat flux', '', 'tt0t')
        call ncinfo(ncname4( 62,:),'ghf_wall_a', 'wall A indoor heat flux', '', 'tt0t')
        call ncinfo(ncname4( 63,:),'ghf_wall_b', 'wall B indoor heat flux', '', 'tt0t')
        call ncinfo(ncname4( 64,:),'ghf_win_a', 'window A indoor heat flux', '', 'tt0t')
        call ncinfo(ncname4( 65,:),'ghf_win_b', 'window B indoor heat flux', '', 'tt0t')
        call ncinfo(ncname4( 66,:),'rad_lw_net_can', 'net longwave radiative at canyon top ', 'downwards', 'tt0t')
        call ncinfo(ncname4( 67,:),'rad_lw_net_road', 'net longtwave radiative flux on road', '', 'tt0t')
        call ncinfo(ncname4( 68,:),'rad_lw_net_roof', 'net longwave radiative flux on roof', '', 'tt0t')
        call ncinfo(ncname4( 69,:),'rad_lw_net_urb', 'urban aggegated net longwave radiative flux', '', 'tt0t')
        call ncinfo(ncname4( 70,:),'rad_lw_net_wall_a', 'net longwave radiative flux on wall A', '', 'tt0t')
        call ncinfo(ncname4( 71,:),'rad_lw_net_wall_b', 'net longwave radiative flux wall B', '', 'tt0t')
        call ncinfo(ncname4( 72,:),'rad_lw_net_win_a', 'net longwave radiative flux on wall A', '', 'tt0t')
        call ncinfo(ncname4( 73,:),'rad_lw_net_win_b', 'net longwave radiative flux window B', '', 'tt0t')
        call ncinfo(ncname4( 74,:),'rad_sw_in_road', 'incoming shortwave radiative flux on road', '', 'tt0t')
        call ncinfo(ncname4( 75,:),'rad_sw_in_win_a', 'incoming shortwave radiative flux on window A', '', 'tt0t')
        call ncinfo(ncname4( 76,:),'rad_sw_in_win_b', 'incoming shortwave radiative flux on window B', '', 'tt0t')
        call ncinfo(ncname4( 77,:),'rad_sw_net_road', 'net shortwave radiative flux on road', '', 'tt0t')
        call ncinfo(ncname4( 78,:),'rad_sw_net_roof', 'net shortwave radiative flux on roof', '', 'tt0t')
        call ncinfo(ncname4( 79,:),'rad_sw_net_urb', 'urban aggegated net shortwave radiative flux', '', 'tt0t')
        call ncinfo(ncname4( 80,:),'rad_sw_net_wall_a', 'net shortwave radiative flux on wall A', '', 'tt0t')
        call ncinfo(ncname4( 81,:),'rad_sw_net_wall_b', 'net shortwave radiative flux on wall B', '', 'tt0t')
        call ncinfo(ncname4( 82,:),'rad_sw_net_win_a', 'net shortwave radiative flux on window A', '', 'tt0t')
        call ncinfo(ncname4( 83,:),'rad_sw_net_win_b', 'net shortwave radiative flux on wall B', '', 'tt0t')
        call ncinfo(ncname4( 84,:),'ol_can', 'canyon top Obukhov length', '', 'tt0t')
        call ncinfo(ncname4( 85,:),'ol_road', 'road Obukhov length', '', 'tt0t')
        call ncinfo(ncname4( 86,:),'ol_roof', 'rroof Obukhov length', '', 'tt0t')
        call ncinfo(ncname4( 87,:),'pt_can', 'street canyon virtual potential temperature ', 'K', 'tt0t')
        call ncinfo(ncname4( 88,:),'rib_can', 'canyon top bulk Richardson number', '', 'tt0t')
        call ncinfo(ncname4( 89,:),'rib_road', 'road bulk Richardson number', '', 'tt0t')
        call ncinfo(ncname4( 90,:),'rib_roof', 'roof bulk Richardson number', '', 'tt0t')
        call ncinfo(ncname4( 91,:),'us_can', 'friction velocity for canyon resistance calculation', '', 'tt0t')
        call ncinfo(ncname4( 92,:),'uv_abs_can', 'horizontal wind speed in street caynon at half-height', '', 'tt0t')
        call ncinfo(ncname4( 93,:),'uv_eff_can', 'effective horizontal wind speed in street canyon at half-height', '', 'tt0t')
        call ncinfo(ncname4( 94,:),'vpt_can', 'street canyon virtual potential temperature ', 'K', 'tt0t')
        call ncinfo(ncname4( 95,:),'rah_can', 'street canyon air aerodynamic resistance for heat', '', 'tt0t')
        call ncinfo(ncname4( 96,:),'rah_facade', 'wall and window aerodynamic resistance for heat ', 'combined', 'tt0t')
        call ncinfo(ncname4( 97,:),'rah_road', 'road aerodynamic resistance for heat', '', 'tt0t')
        call ncinfo(ncname4( 98,:),'rah_roof', 'roof aerodynamic resistance for heat', '', 'tt0t')
        call ncinfo(ncname4( 99,:),'rah_wall_a', 'wall A aerodynamic resistance for heat', '', 'tt0t')
        call ncinfo(ncname4( 100,:),'rah_wall_b', 'wall B aerodynamic resistance for heat', '', 'tt0t')
        call ncinfo(ncname4( 101,:),'rah_win_a', 'wall A aerodynamic resistance for heat', '', 'tt0t')
        call ncinfo(ncname4( 102,:),'rah_win_b', 'wall B aerodynamic resistance for heat', '', 'tt0t')
        call ncinfo(ncname4( 103,:),'us_road', 'friction velocity for roads', '', 'tt0t')
        call ncinfo(ncname4( 104,:),'us_roof', 'friction velocity for roofs', '', 'tt0t')
        call ncinfo(ncname4( 105,:),'pt1', 'potential temperature', '', 'tt0t')
        call ncinfo(ncname4( 106,:),'q1', 'specific humidity', '', 'tt0t')
        call ncinfo(ncname4( 107,:),'us_urb', 'friction velocity', '', 'tt0t')
        call ncinfo(ncname4( 108,:),'uv_abs1', 'horizontal wind speed', '', 'tt0t')
        call ncinfo(ncname4( 109,:),'uv_eff1', 'effective horizontal wind speed', '', 'tt0t')
        call ncinfo(ncname4( 110,:),'vpt1', 'virtual potential temperature', '', 'tt0t')
        call ncinfo(ncname4( 111,:),'f_bld', 'fractional area occupied by buldings ', 'plan area fraction', 'tt0t')
        call ncinfo(ncname4( 112,:),'f_bld_frn', 'frontal area fraction of buildings', '', 'tt0t')
        call ncinfo(ncname4( 113,:),'f_win', 'window fraction', '', 'tt0t')
        call ncinfo(ncname4( 114,:),'h_bld', 'building height', '', 'tt0t')
        call ncinfo(ncname4( 115,:),'hw_can', 'canyon aspect ratio', '', 'tt0t')
        call ncinfo(ncname4( 116,:),'svf_road', 'sky-view factor for road', '', 'tt0t')
        call ncinfo(ncname4( 117,:),'svf_wall', 'sky-view-factor for walls', '', 'tt0t')
        call ncinfo(ncname4( 118,:),'theta_can', 'canyon orientation / road direction in radians', '', 'tt0t')
        call ncinfo(ncname4( 119,:),'z0_urb', 'aerodynamic roughness length of the urban surface', '', 'tt0t')
        call ncinfo(ncname4( 120,:),'albedo_road', 'albedo of the road', '', 'tt0t')
        call ncinfo(ncname4( 121,:),'albedo_roof', 'albedo of the roof', '', 'tt0t')
        call ncinfo(ncname4( 122,:),'albedo_wall', 'albedo of the wall', '', 'tt0t')
        call ncinfo(ncname4( 123,:),'albedo_wall_win', 'weighted average of wall and window albedos for reflections', '', 'tt0t')
        call ncinfo(ncname4( 124,:),'albedo_win', 'albedo of the window', '', 'tt0t')
        call ncinfo(ncname4( 125,:),'emiss_road', 'emissivity of the road', '', 'tt0t')
        call ncinfo(ncname4( 126,:),'emiss_roof', 'emissivity of the roof', '', 'tt0t')
        call ncinfo(ncname4( 127,:),'emiss_wall', 'emissivity of the wall', '', 'tt0t')
        call ncinfo(ncname4( 128,:),'emiss_win', 'emissivity of the window', '', 'tt0t')
        call ncinfo(ncname4( 129,:),'transmissivity_win', 'transmissivity of the window layers', '', 'tt0t')
        call ncinfo(ncname4( 130,:),'z0_road', 'aerodynamic roughness length for momentum for roads', '', 'tt0t')
        call ncinfo(ncname4( 131,:),'z0_roof', 'aerodynamic roughness length for momentum of roofs', '', 'tt0t')
        call ncinfo(ncname4( 132,:),'z0_wall', 'aerodynamic roughness length for walls and windows', '', 'tt0t')
        call ncinfo(ncname4( 133,:),'z0h_road', 'aerodynamic roughness length for heat for roads', '', 'tt0t')
        call ncinfo(ncname4( 134,:),'z0h_roof', 'aerodynamic roughness length for heat for roofs', '', 'tt0t')
        call ncinfo(ncname4( 135,:),'t_indoor', 'building indoor temperature ', 'K', 'tt0t')
        call ncinfo(ncname4( 136,:),'t_soil', 'fixed soil top temperature ', 'K', 'tt0t')
        call ncinfo(ncname4( 137,:),'sw_ref_denom', 'SW radiation reflection denominator', '', 'tt0t')
        call ncinfo(ncname4( 138,:),'uv_abs_can_coef', 'coefficient for the canyon wind speed', '', 'tt0t')
        call ncinfo(ncname4( 139,:),'wall_hor_a_ratio', 'wall-to-horizontal area ratio', '', 'tt0t')
        call ncinfo(ncname4( 140,:),'z_mo', 'reference height for MOST for the atmosphere', '', 'tt0t')
        call ncinfo(ncname4( 141,:),'z_mo_can', 'canyon reference height for MOST ', 'canyon half-height', 'tt0t')
        call ncinfo(ncname4( 142,:),'tm_roof_runoff', 'tm_roof_runoff', 'm s^-1', 'tt0t')
        call ncinfo(ncname4( 143,:),'tm_road_runoff', 'tm_road_runoff', 'm s^-1', 'tt0t')
        call ncinfo(ncname4( 144,:),'tm_roof_precep', 'tm_roof_precep', 'm s^-1', 'tt0t')
        call ncinfo(ncname4( 145,:),'tm_road_precep', 'tm_road_precep', 'm s^-1', 'tt0t')

        call ncinfo(ncname4( nvar4-nvars3d4+1,:),'tt_wall_a', 'tendency wall a ', 'tt_wall_a', 'tttts_slurb')
        call ncinfo(ncname4( nvar4-nvars3d4+2,:),'tt_wall_b', 'tendency wall b ', 'tt_wall_b', 'tttts_slurb')
        call ncinfo(ncname4( nvar4-nvars3d4+3,:),'tt_roof', 'tendency roof a ', 'tt_roof', 'tttts_slurb')
        ! call ncinfo(ncname4( 145,:),'tt_roof_b', 'tendency roof b ', 'tt_roof_b', 'tt0t')
        call ncinfo(ncname4( nvar4-nvars3d4+4,:),'tt_win_a', 'tendency win a ', 'tt_win_a', 'tttts_slurb')
        call ncinfo(ncname4( nvar4-nvars3d4+5,:),'tt_win_b', 'tendency win b ', 'tt_win_b', 'tttts_slurb')
        call ncinfo(ncname4( nvar4-nvars3d4+6,:),'tt_road', 'tendency road a ', 'tt_road', 'tttts_slurb')
        ! call ncinfo(ncname4( 149,:),'tt_road_b', 'tendency road b ', 'tt_road_b', 'tt0t')
        call ncinfo(ncname4( nvar4-nvars3d4+7,:),'c_wall', 'c_wall ', 'c_wall', 'tttts_slurb')
        call ncinfo(ncname4( nvar4-nvars3d4+8,:),'c_roof', 'c_roof ', 'c_roof', 'tttts_slurb')
        call ncinfo(ncname4( nvar4-nvars3d4+9,:),'c_win', 'c_win ', 'c_win', 'tttts_slurb')
        call ncinfo(ncname4( nvar4-nvars3d4+10,:),'c_road', 'c_road ', 'c_road', 'tttts_slurb')
        call ncinfo(ncname4( nvar4-nvars3d4+11,:),'absorption_win', 'absorption_win ', 'absorption_win', 'tttts_slurb')
        call ncinfo(ncname4( nvar4-nvars3d4+12,:),'dz_wall', 'dz_wall ', 'dz_wall', 'tttts_slurb')
        call ncinfo(ncname4( nvar4-nvars3d4+13,:),'dz_roof', 'dz_roof ', 'dz_roof', 'tttts_slurb')
        call ncinfo(ncname4( nvar4-nvars3d4+14,:),'dz_win', 'dz_win ', 'dz_win', 'tttts_slurb')
        call ncinfo(ncname4( nvar4-nvars3d4+15,:),'dz_road', 'dz_road ', 'dz_road', 'tttts_slurb')
        call ncinfo(ncname4( nvar4-nvars3d4+16,:),'t_wall_a', 'temperature wall a ', 't_wall_a', 'tttts_slurb')
        call ncinfo(ncname4( nvar4-nvars3d4+17,:),'t_wall_b', 'temperature wall b ', 't_wall_b', 'tttts_slurb')
        call ncinfo(ncname4( nvar4-nvars3d4+18,:),'t_roof', 'temperature roof a ', 't_roof', 'tttts_slurb')
        ! call ncinfo(ncname4( 145,:),'tt_roof_b', 'tendency roof b ', 'tt_roof_b', 'tt0t')
        call ncinfo(ncname4( nvar4-nvars3d4+19,:),'t_win_a', 'temperature win a ', 't_win_a', 'tttts_slurb')
        call ncinfo(ncname4( nvar4-nvars3d4+20,:),'t_win_b', 'temperature win b ', 't_win_b', 'tttts_slurb')
        call ncinfo(ncname4( nvar4-nvars3d4+21,:),'t_road', 'temperature road a ', 't_road', 'tttts_slurb')

        


        call open_nc(trim(output_prefix)//fname4,  ncid4,nrec4,n1=imax,n2=jmax,ns=4)
        if (nrec4==0) then
          call define_nc(ncid4, 1, tncname4)
          call writestat_dims_nc(ncid4)
        end if
        call define_nc( ncid4, nvar4, ncname4)
    end if
  end if
  end subroutine initlsmcrosssection
!>Run lsmcrosssection. Mainly timekeeping
  subroutine lsmcrosssection
    use modglobal, only : rk3step,timee,dt_lim
    use modstat_nc, only : writestat_nc
    implicit none


    if (.not. (lcross .or. lcrosssoil)) return
    if (rk3step/=3) return
    if(timee<tnext) then
      dt_lim = min(dt_lim,tnext-timee)
      return
    end if
    tnext = tnext+idtav
    dt_lim = minval((/dt_lim,tnext-timee/))

    if (lcrosssoil) then
       call wrtvert
       call wrthorz
    end if
    if (lcross) then
       call wrtsurf
       call wrtslurb
    end if
  end subroutine lsmcrosssection


!> Do the xz lsmcrosssections and dump them to file
  subroutine wrtvert
  use modglobal, only : imax,i1,cexpnr,ifoutput,rtimee
  use modsurfdata, only : tsoil, phiw
  use modmpi,    only : myidy
  use modstat_nc, only : lnetcdf, writestat_nc
  implicit none

  integer i,k

  real, allocatable :: vars(:,:,:)

  if( myidy /= 0 ) return

    open(ifoutput,file='movv_tsoil.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((tsoil(i,crossplane,k),i=2,i1),k=1,ksoilmax)
    close(ifoutput)

    open(ifoutput,file='movv_phiw.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((phiw(i,crossplane,k),i=2,i1),k=1,ksoilmax)
    close(ifoutput)
    if (lnetcdf) then
      allocate(vars(1:imax,1:ksoilmax,2))
      vars(:,:,1) = tsoil(2:i1,crossplane,1:ksoilmax)
      vars(:,:,2) = phiw(2:i1,crossplane,1:ksoilmax)
      call writestat_nc(ncid1,1,tncname1,(/rtimee/),nrec1,.true.)
      call writestat_nc(ncid1,2,ncname1(1:2,:),vars,nrec1,imax,ksoilmax)
      deallocate(vars)
    end if

  end subroutine wrtvert

!> Do the xy lsmcrosssections and dump them to file
  subroutine wrthorz
    use modglobal, only : imax,jmax,i1,j1,cexpnr,ifoutput,rtimee
    use modsurfdata, only : tsoil,phiw
    use modstat_nc, only : lnetcdf, writestat_nc
    implicit none

    ! LOCAL
    integer i,j
    real, allocatable :: vars(:,:,:)

    write(cheight,'(i4.4)') crossheight
    open(ifoutput,file='movh_tsoil.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((tsoil(i,j,crossheight),i=2,i1),j=2,j1)
    close(ifoutput)

    open(ifoutput,file='movh_phiw.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((phiw(i,j,crossheight),i=2,i1),j=2,j1)
    close(ifoutput)

    if (lnetcdf) then
        allocate(vars(1:imax,1:jmax,2))
        vars(:,:,1) = tsoil(2:i1,2:j1,crossheight)
        vars(:,:,2) = phiw(2:i1,2:j1,crossheight)
        call writestat_nc(ncid2,1,tncname2,(/rtimee/),nrec2,.true.)
        call writestat_nc(ncid2,2,ncname2(1:2,:),vars,nrec2,imax,jmax)
        deallocate(vars)
    end if

  end subroutine wrthorz

  !> Do the xy lsmcrosssections and dump them to file
  subroutine wrtsurf
    use modglobal, only : imax,jmax,i1,j1,rtimee,cp,rlv
    use modfields, only : rhof
    use modsurfdata, only : Qnet, H, LE, G0, rs, ra, tskin, tendskin, &
                            cliq, rsveg, rssoil, Wl, isurf, obl, ustar, &
                            Cs, Cm, z0h, z0m, qtflux, thlflux
    use modlsm, only : f1, f2b, lags, an_co2, resp_co2
    use modstat_nc, only : lnetcdf, writestat_nc
    implicit none

    ! LOCAL
    real, allocatable :: vars(:,:,:)
    if (lnetcdf) then
        allocate(vars(1:imax,1:jmax,nvar3))

        if (isurf == 1) then
            vars(:,:,1) = qnet(2:i1,2:j1)
            vars(:,:,2) = h(2:i1,2:j1)
            vars(:,:,3) = le(2:i1,2:j1)
            vars(:,:,4) = g0(2:i1,2:j1)
            vars(:,:,5) = tskin(2:i1,2:j1)
            vars(:,:,6) = tendskin(2:i1,2:j1)
            vars(:,:,7) = rs(2:i1,2:j1)
            vars(:,:,8) = ra(2:i1,2:j1)
            vars(:,:,9) = cliq(2:i1,2:j1)
            vars(:,:,10) = Wl(2:i1,2:j1)
            vars(:,:,11) = rssoil(2:i1,2:j1)
            vars(:,:,12) = rsveg(2:i1,2:j1)
        else if (isurf == 2) then
            vars(:,:, 1) = rhof(1) * cp * thlflux(2:i1,2:j1)  ! H not allocated for isurf=2
            vars(:,:, 2) = rhof(1) * rlv * qtflux (2:i1,2:j1) ! LE not allocated for isurf=2
            vars(:,:, 3) = obl(2:i1,2:j1)
            vars(:,:, 4) = ustar(2:i1,2:j1)
            vars(:,:, 5) = Cs(2:i1,2:j1)
            vars(:,:, 6) = Cm(2:i1,2:j1)
            vars(:,:, 7) = z0h(2:i1,2:j1)
            vars(:,:, 8) = z0m(2:i1,2:j1)
        else if (isurf == 11) then
            !$acc update host(H, LE, G0, tskin, obl, ustar, cliq, Wl, ra, &
            !$acc& rssoil, rsveg, f1, f2b)
            vars(:,:, 1) = H(2:i1,2:j1)
            vars(:,:, 2) = LE(2:i1,2:j1)
            vars(:,:, 3) = G0(2:i1,2:j1)
            vars(:,:, 4) = tskin(2:i1,2:j1)
            vars(:,:, 5) = obl(2:i1,2:j1)
            vars(:,:, 6) = ustar(2:i1,2:j1)
            vars(:,:, 7) = cliq(2:i1,2:j1)
            vars(:,:, 8) = Wl(2:i1,2:j1)
            vars(:,:, 9) = ra(2:i1,2:j1)
            vars(:,:,10) = rssoil(2:i1,2:j1)
            vars(:,:,11) = rsveg(2:i1,2:j1)
            vars(:,:,12) = f1(2:i1,2:j1)
            vars(:,:,13) = f2b(2:i1,2:j1)
            vars(:,:,14) = Qnet(2:i1,2:j1)
            if (lags) then
              vars(:,:,15) = an_co2(2:i1,2:j1)
              vars(:,:,16) = resp_co2(2:i1,2:j1)
            endif
        end if

        call writestat_nc(ncid3, 1, tncname3, (/rtimee/), nrec3, .true.)
        call writestat_nc(ncid3, nvar3, ncname3(1:nvar3,:), vars, nrec3, imax, jmax)

        deallocate(vars)
    end if


  end subroutine wrtsurf
  
  subroutine wrtslurb
    use modglobal, only : imax,jmax,i1,j1,rtimee,rlv,cp
    use modfields, only : rhof
    use modslurbdata, only : slurb_tile, facade_rah_doe, enable_slurb
    use modstat_nc, only : lnetcdf, writestat_nc
    implicit none

    ! LOCAL
    real, allocatable :: vars(:,:,:)
    real, allocatable :: vars3d(:,:,:,:)

    if (lnetcdf.and.enable_slurb) then
        allocate(vars(1:imax,1:jmax,nvar4-nvars3d4))
        allocate(vars3d(4,1:imax,1:jmax,nvars3d4))

        vars(:,:,1) = slurb_tile%albedo_urb(2:i1,2:j1)
        vars(:,:,2) = slurb_tile%emiss_urb(2:i1,2:j1)
        vars(:,:,3) = slurb_tile%ol_urb(2:i1,2:j1)
        vars(:,:,4) = slurb_tile%qsws_urb(2:i1,2:j1)
        vars(:,:,5) = slurb_tile%rad_lw_in_urb(2:i1,2:j1)
        vars(:,:,6) = slurb_tile%rad_lw_out_urb(2:i1,2:j1)
        vars(:,:,7) = slurb_tile%rad_sw_in_urb(2:i1,2:j1)
        vars(:,:,8) = slurb_tile%rad_sw_out_urb(2:i1,2:j1)
        vars(:,:,9) = slurb_tile%ram_urb(2:i1,2:j1)
        vars(:,:,10) = slurb_tile%rib_urb(2:i1,2:j1)
        vars(:,:,11) = slurb_tile%shf_urb(2:i1,2:j1)
        vars(:,:,12) = slurb_tile%t_2m_urb(2:i1,2:j1)
        vars(:,:,13) = slurb_tile%t_c_urb(2:i1,2:j1)
        vars(:,:,14) = slurb_tile%t_h_urb(2:i1,2:j1)
        vars(:,:,15) = slurb_tile%thl_rad_urb(2:i1,2:j1)
        vars(:,:,16) = slurb_tile%usws_urb(2:i1,2:j1)
        vars(:,:,17) = slurb_tile%vsws_urb(2:i1,2:j1)
        vars(:,:,18) = slurb_tile%thlskin(2:i1,2:j1)
        vars(:,:,19) = slurb_tile%qtskin(2:i1,2:j1)
        vars(:,:,20) = slurb_tile%m_liq_road_0(2:i1,2:j1)
        vars(:,:,21) = slurb_tile%m_liq_road_m(2:i1,2:j1)
        vars(:,:,22) = slurb_tile%m_liq_roof_0(2:i1,2:j1)
        vars(:,:,23) = slurb_tile%m_liq_roof_m(2:i1,2:j1)
        vars(:,:,24) = slurb_tile%q_can_0(2:i1,2:j1)
        vars(:,:,25) = slurb_tile%q_can_m(2:i1,2:j1)
        vars(:,:,26) = slurb_tile%t_can_0(2:i1,2:j1)
        vars(:,:,27) = slurb_tile%t_can_m(2:i1,2:j1)
        vars(:,:,28) = slurb_tile%tm_liq_road(2:i1,2:j1)
        vars(:,:,29) = slurb_tile%tm_liq_roof(2:i1,2:j1)
        vars(:,:,30) = slurb_tile%tq_can(2:i1,2:j1) * rhof(1)
        ! K s^-1 * J K^-1 kg^-1 * kg m^-3
        !  W m^-3
        vars(:,:,31) = slurb_tile%tt_can(2:i1,2:j1) * cp * rhof(1)
        vars(:,:,32) = slurb_tile%pt_road(2:i1,2:j1)
        vars(:,:,33) = slurb_tile%pt_roof(2:i1,2:j1)
        vars(:,:,34) = slurb_tile%pt_wall_a(2:i1,2:j1)
        vars(:,:,35) = slurb_tile%pt_wall_b(2:i1,2:j1)
        vars(:,:,36) = slurb_tile%pt_win_a(2:i1,2:j1)
        vars(:,:,37) = slurb_tile%pt_win_b(2:i1,2:j1)
        vars(:,:,38) = slurb_tile%q_road(2:i1,2:j1)
        vars(:,:,39) = slurb_tile%q_roof(2:i1,2:j1)
        vars(:,:,40) = slurb_tile%qs_road(2:i1,2:j1)
        vars(:,:,41) = slurb_tile%qs_roof(2:i1,2:j1)
        vars(:,:,42) = slurb_tile%vpt_road(2:i1,2:j1)
        vars(:,:,43) = slurb_tile%vpt_roof(2:i1,2:j1)
        vars(:,:,43) = slurb_tile%shf_can(2:i1,2:j1)
        vars(:,:,44) = slurb_tile%shf_external(2:i1,2:j1)
        vars(:,:,45) = slurb_tile%shf_road(2:i1,2:j1)
        vars(:,:,46) = slurb_tile%shf_roof(2:i1,2:j1)
        vars(:,:,47) = 0!slurb_tile%shf_traffic(2:i1,2:j1)
        vars(:,:,48) = slurb_tile%shf_wall_a(2:i1,2:j1)
        vars(:,:,49) = slurb_tile%shf_wall_b(2:i1,2:j1)
        vars(:,:,50) = slurb_tile%shf_win_a(2:i1,2:j1)
        vars(:,:,51) = slurb_tile%shf_win_b(2:i1,2:j1)
        vars(:,:,52) = slurb_tile%qsws_can(2:i1,2:j1)
        vars(:,:,53) = slurb_tile%qsws_external(2:i1,2:j1)
        vars(:,:,54) = slurb_tile%qsws_liq_road(2:i1,2:j1)
        vars(:,:,55) = slurb_tile%qsws_liq_roof(2:i1,2:j1)
        vars(:,:,56) = slurb_tile%qsws_road(2:i1,2:j1)
        vars(:,:,57) = slurb_tile%qsws_roof(2:i1,2:j1)
        vars(:,:,58) = slurb_tile%c_liq_road(2:i1,2:j1)
        vars(:,:,59) = slurb_tile%c_liq_roof(2:i1,2:j1)
        vars(:,:,60) = slurb_tile%ghf_road(2:i1,2:j1)
        vars(:,:,61) = slurb_tile%ghf_roof(2:i1,2:j1)
        vars(:,:,62) = slurb_tile%ghf_wall_a(2:i1,2:j1)
        vars(:,:,63) = slurb_tile%ghf_wall_b(2:i1,2:j1)
        vars(:,:,64) = slurb_tile%ghf_win_a(2:i1,2:j1)
        vars(:,:,65) = slurb_tile%ghf_win_b(2:i1,2:j1)
        vars(:,:,66) = slurb_tile%rad_lw_net_can(2:i1,2:j1)
        vars(:,:,67) = slurb_tile%rad_lw_net_road(2:i1,2:j1)
        vars(:,:,68) = slurb_tile%rad_lw_net_roof(2:i1,2:j1)
        vars(:,:,69) = slurb_tile%rad_lw_net_urb(2:i1,2:j1)
        vars(:,:,70) = slurb_tile%rad_lw_net_wall_a(2:i1,2:j1)
        vars(:,:,71) = slurb_tile%rad_lw_net_wall_b(2:i1,2:j1)
        vars(:,:,72) = slurb_tile%rad_lw_net_win_a(2:i1,2:j1)
        vars(:,:,73) = slurb_tile%rad_lw_net_win_b(2:i1,2:j1)
        vars(:,:,74) = slurb_tile%rad_sw_in_road(2:i1,2:j1)
        vars(:,:,75) = slurb_tile%rad_sw_in_win_a(2:i1,2:j1)
        vars(:,:,76) = slurb_tile%rad_sw_in_win_b(2:i1,2:j1)
        vars(:,:,77) = slurb_tile%rad_sw_net_road(2:i1,2:j1)
        vars(:,:,78) = slurb_tile%rad_sw_net_roof(2:i1,2:j1)
        vars(:,:,79) = slurb_tile%rad_sw_net_urb(2:i1,2:j1)
        vars(:,:,80) = slurb_tile%rad_sw_net_wall_a(2:i1,2:j1)
        vars(:,:,81) = slurb_tile%rad_sw_net_wall_b(2:i1,2:j1)
        vars(:,:,82) = slurb_tile%rad_sw_net_win_a(2:i1,2:j1)
        vars(:,:,83) = slurb_tile%rad_sw_net_win_b(2:i1,2:j1)
        vars(:,:,84) = slurb_tile%ol_can(2:i1,2:j1)
        vars(:,:,85) = slurb_tile%ol_road(2:i1,2:j1)
        vars(:,:,86) = slurb_tile%ol_roof(2:i1,2:j1)
        vars(:,:,87) = slurb_tile%pt_can(2:i1,2:j1)
        vars(:,:,88) = slurb_tile%rib_can(2:i1,2:j1)
        vars(:,:,89) = slurb_tile%rib_road(2:i1,2:j1)
        vars(:,:,90) = slurb_tile%rib_roof(2:i1,2:j1)
        vars(:,:,91) = slurb_tile%us_can(2:i1,2:j1)
        vars(:,:,92) = slurb_tile%uv_abs_can(2:i1,2:j1)
        vars(:,:,93) = slurb_tile%uv_eff_can(2:i1,2:j1)
        vars(:,:,94) = slurb_tile%vpt_can(2:i1,2:j1)
        vars(:,:,95) = slurb_tile%rah_can(2:i1,2:j1)
        vars(:,:,97) = slurb_tile%rah_road(2:i1,2:j1)
        vars(:,:,98) = slurb_tile%rah_roof(2:i1,2:j1)
        if (facade_rah_doe) then
          vars(:,:,96) = 0
          vars(:,:,99) = slurb_tile%rah_wall_a(2:i1,2:j1)
          vars(:,:,100) = slurb_tile%rah_wall_b(2:i1,2:j1)
          vars(:,:,101) = slurb_tile%rah_win_a(2:i1,2:j1)
          vars(:,:,102) = slurb_tile%rah_win_b(2:i1,2:j1)
        else
          vars(:,:,96) = slurb_tile%rah_facade(2:i1,2:j1)
          vars(:,:,99) = 0
          vars(:,:,100) = 0
          vars(:,:,101) =0
          vars(:,:,102) =0
        endif
        vars(:,:,103) = slurb_tile%us_road(2:i1,2:j1)
        vars(:,:,104) = slurb_tile%us_roof(2:i1,2:j1)
        vars(:,:,105) = slurb_tile%pt1(2:i1,2:j1)
        vars(:,:,106) = slurb_tile%q1(2:i1,2:j1)
        vars(:,:,107) = slurb_tile%us_urb(2:i1,2:j1)
        vars(:,:,108) = slurb_tile%uv_abs1(2:i1,2:j1)
        vars(:,:,109) = slurb_tile%uv_eff1(2:i1,2:j1)
        vars(:,:,110) = slurb_tile%vpt1(2:i1,2:j1)
        vars(:,:,111) = slurb_tile%f_bld(2:i1,2:j1)
        vars(:,:,112) = slurb_tile%f_bld_frn(2:i1,2:j1)
        vars(:,:,113) = slurb_tile%f_win(2:i1,2:j1)
        vars(:,:,114) = slurb_tile%h_bld(2:i1,2:j1)
        vars(:,:,115) = slurb_tile%hw_can(2:i1,2:j1)
        vars(:,:,116) = slurb_tile%svf_road(2:i1,2:j1)
        vars(:,:,117) = slurb_tile%svf_wall(2:i1,2:j1)
        vars(:,:,118) = slurb_tile%theta_can(2:i1,2:j1)
        vars(:,:,119) = slurb_tile%z0_urb(2:i1,2:j1)
        vars(:,:,120) = slurb_tile%albedo_road(2:i1,2:j1)
        vars(:,:,121) = slurb_tile%albedo_roof(2:i1,2:j1)
        vars(:,:,122) = slurb_tile%albedo_wall(2:i1,2:j1)
        vars(:,:,123) = slurb_tile%albedo_wall_win(2:i1,2:j1)
        vars(:,:,124) = slurb_tile%albedo_win(2:i1,2:j1)
        vars(:,:,125) = slurb_tile%emiss_road(2:i1,2:j1)
        vars(:,:,126) = slurb_tile%emiss_roof(2:i1,2:j1)
        vars(:,:,127) = slurb_tile%emiss_wall(2:i1,2:j1)
        vars(:,:,128) = slurb_tile%emiss_win(2:i1,2:j1)
        vars(:,:,129) = slurb_tile%transmissivity_win(2:i1,2:j1)
        vars(:,:,130) = slurb_tile%z0_road(2:i1,2:j1)
        vars(:,:,131) = slurb_tile%z0_roof(2:i1,2:j1)
        vars(:,:,132) = slurb_tile%z0_wall(2:i1,2:j1)
        vars(:,:,133) = slurb_tile%z0h_road(2:i1,2:j1)
        vars(:,:,134) = slurb_tile%z0h_roof(2:i1,2:j1)
        vars(:,:,135) = slurb_tile%t_indoor(2:i1,2:j1)
        vars(:,:,136) = slurb_tile%t_soil(2:i1,2:j1)
        vars(:,:,137) = slurb_tile%sw_ref_denom(2:i1,2:j1)
        vars(:,:,138) = slurb_tile%uv_abs_can_coef(2:i1,2:j1)
        vars(:,:,139) = 0!slurb_tile%wall_hor_a_ratio(2:i1,2:j1)
        vars(:,:,140) = slurb_tile%z_mo(2:i1,2:j1)
        vars(:,:,141) = slurb_tile%z_mo_can(2:i1,2:j1)
        vars(:,:,142) = slurb_tile%tm_roof_runoff(2:i1,2:j1)
        vars(:,:,143) = slurb_tile%tm_road_runoff(2:i1,2:j1)
        vars(:,:,144) = slurb_tile%tm_roof_precep(2:i1,2:j1) * rhof(1)
        vars(:,:,145) = slurb_tile%tm_road_precep(2:i1,2:j1) * rhof(1)

        vars3d(:,:,:,1) = slurb_tile%tt_wall_a(:,2:i1,2:j1)
        vars3d(:,:,:,2) = slurb_tile%tt_wall_b(:,2:i1,2:j1)
        vars3d(:,:,:,3) = slurb_tile%tt_roof(:,2:i1,2:j1)
        ! vars(:,:,145) = slurb_tile%tt_roof_b(1,2:i1,2:j1)
        vars3d(:,:,:,4) = slurb_tile%tt_win_a(:,2:i1,2:j1)
        vars3d(:,:,:,5) = slurb_tile%tt_win_b(:,2:i1,2:j1)
        vars3d(:,:,:,6) = slurb_tile%tt_road(:,2:i1,2:j1)
        ! vars(:,:,149) = slurb_tile%tt_road_b(1,2:i1,2:j1)
        vars3d(:,:,:,7) = slurb_tile%c_wall(:,2:i1,2:j1)
        vars3d(:,:,:,8) = slurb_tile%c_roof(:,2:i1,2:j1)
        vars3d(:,:,:,9) = slurb_tile%c_win(:,2:i1,2:j1)
        vars3d(:,:,:,10) = slurb_tile%c_road(:,2:i1,2:j1)
        vars3d(:,:,:,11) = slurb_tile%absorption_win(:,2:i1,2:j1)
        vars3d(:,:,:,12) = slurb_tile%dz_wall(:,2:i1,2:j1)
        vars3d(:,:,:,13) = slurb_tile%dz_roof(:,2:i1,2:j1)
        vars3d(:,:,:,14) = slurb_tile%dz_win(:,2:i1,2:j1)
        vars3d(:,:,:,15) = slurb_tile%dz_road(:,2:i1,2:j1)
        vars3d(:,:,:,16) = slurb_tile%t_wall_a_0(:,2:i1,2:j1)
        vars3d(:,:,:,17) = slurb_tile%t_wall_b_0(:,2:i1,2:j1)
        vars3d(:,:,:,18) = slurb_tile%t_roof_0(:,2:i1,2:j1)
        ! vars(:,:,145) = slurb_tile%tt_roof_b(1,2:i1,2:j1)
        vars3d(:,:,:,19) = slurb_tile%t_win_a_0(:,2:i1,2:j1)
        vars3d(:,:,:,20) = slurb_tile%t_win_b_0(:,2:i1,2:j1)
        vars3d(:,:,:,21) = slurb_tile%t_road_0(:,2:i1,2:j1)



        call writestat_nc(ncid4, 1, tncname4, (/rtimee/), nrec4, .true.)
        call writestat_nc(ncid4, nvar4-nvars3d4, ncname4(1:nvar4-nvars3d4,:), vars, nrec4, imax, jmax)
        call writestat_nc(ncid4, nvars3d4, ncname4(nvar4-nvars3d4+1:nvar4,:), vars3d, nrec4, 4, imax, jmax)

        deallocate(vars)
        deallocate(vars3d)
    end if
  end subroutine wrtslurb
!> Clean up when leaving the run
  subroutine exitlsmcrosssection
    use modstat_nc, only : exitstat_nc,lnetcdf
    use modmpi, only : myidy
    use modslurb, only : enable_slurb
    implicit none

    if(lcrosssoil .and. lnetcdf) then
      if (myidy==0) then
        call exitstat_nc(ncid1) ! xz soil
      end if
      call exitstat_nc(ncid2)   ! xy soil
    end if
    if(lcross .and. lnetcdf) then
       call exitstat_nc(ncid3) ! surface
       deallocate(ncname3)
       if (enable_slurb) then
         deallocate(ncname4)
         call exitstat_nc(ncid4)
       end if
    end if
  end subroutine exitlsmcrosssection

end module modlsmcrosssection
