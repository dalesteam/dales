!> \file modtimestat.f90
!!  Timestat calculates timeseries of several variables
!>
!! Timestat calculates timeseries of several variables
!>
!! Timeseries of the most relevant parameters. Written to tmser1.expnr and tmsurf.expnr
!! If netcdf is true, this module leads the tmser.expnr.nc output
!!  \author Pier Siebesma, K.N.M.I.
!!  \author Stephan de Roode,TU Delft
!!  \author Chiel van Heerwaarden, Wageningen U.R.
!!  \author Thijs Heus,MPI-M
!!  \rewritten by Johan van der Dussen TU Delft
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module modtimestat

use modglobal, only   : longint

implicit none
save

! Initialize module-wide variables
logical              :: ltimestat                             ! On/off switch time statistics
real                 :: dtav                                  ! Averaging timestep (in s)
integer(kind=longint):: idtav                , &              ! Averaging timestep (in ms)
                        tnext                                 ! Time for next entry in timeseries (in s)
integer,parameter    :: iblh_flux = 1        , &              ! Flux method of BL height calculation
                        iblh_grad = 2        , &              ! Gradient method
                        iblh_thres= 3        , &              ! Threshold method
                        iblh_thv  =-1        , &              ! Virtual potential temperature is used
                        iblh_thl  =-2        , &              ! Liquid water potential temperature
                        iblh_qt   =-3                         ! Total water content
integer              :: iblh_meth = iblh_grad, &              ! Selection of BL height calculation method
                        iblh_var  = iblh_thl                  ! Selection of variable to base BL height on 
real                 :: blh_thres = -1                        ! Threshold used if iblh_meth=iblh_thres=3
integer              :: blh_nsamp = 4        , &              ! Number of samples for fancy BL height calculation
                        ncid                 , &              ! ID of NetCDF file
                        nRec                                  ! Record number
integer,parameter    :: nVarMax=50                            ! Maximum number of variables (arbitrary, increase if necessary)
real                 :: vars(nVarMax)                         ! Holds the variables to be written
real,parameter       :: missVal=-999.                         ! Missing value (useful in netcdf files)
! Initialize output variables
real                 :: cc                   , &              ! Cloud cover [-]
                        zbaseav              , &              ! Average cloud base height [m]
                        ztopav               , &              ! Average cloud top height [m]
                        ztopmax              , &              ! Maximum cloud top height [m]
                        zi                   , &              ! Inversion height [m]
                        we                   , &              ! Entrainment rate [m/s]
                        qlintav              , &              ! Average liquid water path [kg/m^2]
                        qlintmax             , &              ! Maximum liquid water path [kg/m^2]
                        wmax                 , &              ! Maximum updraft velocity [m/s]
                        tkeint               , &              ! Average vertically integrated TKE [m^3/s^2]
                        qlmax                , &              ! Maximum liquid water content [kg/kg]
                        ust                  , &              ! Friction velocity [m/s]
                        tst                  , &              ! Turbulent temperature scale [K]
                        qst                  , &              ! Turbulent humidity scale [kg/kg]
                        wts                  , &              ! Surface kinematic temperature flux [K m/s]
                        wtvs                 , &              ! Surface kinematic virtual temp. flux [K m/s]
                        wqts                                  ! Surface humidity flux [kg m/kg s]
! Surface scheme variables
real                 :: Qnetav,Hav,LEav,G0av,tendskinav, &
                        rsav,raav,tskinav,cliqav,wlav,   &
                        rssoilav,rsvegav
! User specified variables
real                 :: zbasemin             , &              ! Minimum cloud base height [m]
                        qlintvar             , &              ! Liquid water path variance [kg^2/m^4]
                        alb_sfc              , &              ! Surface albedo [-]
                        alb_tod              , &              ! Albedo just above cloud top [-]
                        alb_toa              , &              ! Albedo at the top of the atmosphere [-]
                        qrintav              , &              ! Rain water path [kg/m^2]
                        Nc_av                , &              ! Average cloud droplet number dens. [/cm^3]
                        Nr_av                , &              ! Average rain droplet number dens. [/cm^3]
                        prec_srf             , &              ! Mean precipitation flux at the surface [kg m/kg s]
                        prec_fracsrf                          ! Fraction of boxes with precip/total at surface [-]

!====== Start of subroutines ========================================================================;
contains
  subroutine inittimestat
    use modglobal,   only : ifnamopt,fname_options,dtav_glob,tres,btime,dt_lim,dtmax,ladaptive
    use modmpi,      only : myid,MY_REAL,MPI_LOGICAL,MPI_INTEGER,comm3d,mpierr
    use modstat_nc,  only : lnetcdf
    implicit none
    integer          :: ierr                                  ! Status specifier

    namelist/NAMTIMESTAT/ &                                   ! Read specific timestat namelist
    ltimestat,dtav,iblh_meth,iblh_var,blh_thres,blh_nsamp     ! Overwrites standard values

    dtav = dtav_glob                                          ! Set equal to global averaging timestep 
    if (myid==0) then
      open(ifnamopt, file=fname_options, status='old', iostat=ierr)
      read(ifnamopt, NAMTIMESTAT                     , iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions NAMTIMESTAT'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMTIMESTAT'
      endif
      write(6 ,NAMTIMESTAT)
      close(ifnamopt)
    end if

    ! Broadcast the variables read in from the namelist to all cores,
    ! since only the first one (with myid=0) read them.
    call MPI_BCAST(dtav     , 1, MY_REAL     ,0,comm3d,mpierr)
    call MPI_BCAST(ltimestat, 1, MPI_LOGICAL ,0,comm3d,mpierr)
    call MPI_BCAST(blh_thres, 1, MY_REAL     ,0,comm3d,mpierr)
    call MPI_BCAST(iblh_meth, 1, MPI_INTEGER ,0,comm3d,mpierr)
    call MPI_BCAST(iblh_var , 1, MPI_INTEGER ,0,comm3d,mpierr)
    call MPI_BCAST(blh_nsamp, 1, MPI_INTEGER ,0,comm3d,mpierr)

    if (.not.(ltimestat)) return                              ! Leave subroutine if ltimestat=.false.
    
    idtav  = dtav  / tres                                     ! Change time from s to ms
    tnext  = idtav + btime                                    ! Absolute time for next entry in timeseries
    dt_lim = min(dt_lim,tnext)                                ! Limits the models timestep if necessary

    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'TIMESTAT: dtav should be a integer multiple of dtmax'
    end if

    call calcBLHeight                                         !
    call timestatAscii                                        ! First time is only for initialization
    if (lnetcdf) then                                         !
      call timestatNC                                         ! 
    end if

  end subroutine inittimestat

!====================================================================================================;

  subroutine timestat
    use modglobal,   only : rk3step,timee,dt_lim
    use modstat_nc,  only : lnetcdf
    implicit none

    integer              :: n

    if (.not.(ltimestat)) return
    if (rk3step/=3) return
    if(timee<tnext) then
      dt_lim = min(dt_lim,tnext-timee)
      return
    end if
    tnext = tnext+idtav
    dt_lim = minval((/dt_lim,tnext-timee/))
    
    ! Start actual calculations !
    call calcBLHeight
    call calcTimestat(n)
    call calcTimestatUser(n)

    ! Write variables to files !
    call timestatAscii
    if (lnetcdf) then
      call timestatNC
    end if

  end subroutine timestat

!====================================================================================================;

  subroutine calcTimestat(n)
    use modglobal,    only : i1,j1,kmax,dzf,zf,rslabs,cu,cv,ep2,rtimee
    use modfields,    only : rhof,ql0,um,vm,wm,u0av,v0av,e12m
    use modsurfdata,  only : ustar,thlflux,qtflux,isurf,wtsurf,wqsurf,qts,thls,z0, &
                             Qnet,H,LE,G0,tendskin,rs,ra,tskin,cliq,wl,rssoil,rsveg,oblav
    use modmicrodata, only : imicro,imicro_bulk,qr
    use modmpi,       only : MY_REAL,MPI_SUM,MPI_MAX,MPI_MIN,comm3d,mpierr
    implicit none

    integer, intent(out) :: n
    ! Every variable that ends with an 'l' is local to the CPU, MPI communication is needed to
    ! find the domain wide value.
    real                 :: qlint            , &              ! Temporary value LWP
                            qrint            , &              ! Temporary value RWP
                            qlint2           , &              ! Temporary value LWP^2
                            zbase            , &              ! Temporary value cloud base
                            ztop             , &              ! Temporary value cloud top
                            ccl              , &              ! Cloud cover
                            qlintavl         , &              ! Sum of LWPs
                            qlintmaxl        , &              ! Maximum LWP
                            qlint2l          , &              ! Sum of LWP variance
                            qrintavl         , &              ! Sum of RWPs
                            zbaseavl         , &              ! Sum of cloud bases
                            zbaseminl        , &              ! Sum of minimum cloud bases
                            ztopavl          , &              ! Sum of cloud tops 
                            ztopmaxl         , &              ! Sum of maximum cloud tops
                            wmaxl            , &              ! Maximum vertical velocity w
                            qlmaxl           , &              ! Maximum liquid water content ql
                            tkeintl                           ! Sum of vertically integrated tke
    real                 :: c1,c2                             ! Constants used for calculation wtvs
    integer              :: i,j,k                             ! Indices

    ! Set all CPU local values to zero or another appropriate value
    ccl = 0; qlintavl = 0; qlintmaxl = 0.; qlint2l = 0.; qrintavl = 0.; zbaseavl = 0.
    zbaseminl = zf(kmax); ztopavl = 0.; ztopmaxl = 0.; wmaxl = 0.; qlmaxl = 0.; tkeintl = 0.

    ! Start do loop over the horizontal domain i,j:
    do j=2,j1
      do i=2,i1
        qlint     = 0.
        qrint     = 0.
        ztop      = 0.

        ! Calculate LWP and RWP
        do k=1,kmax
          if (imicro == imicro_bulk) then
            qrint = qrint + qr (i,j,k)*rhof(k)*dzf(k)
          end if
          qlint   = qlint + ql0(i,j,k)*rhof(k)*dzf(k)
        end do

        ! Calculate LWP (variance) and cloud cover
        qlint2l = qlint2l + qlint**2
        if (qlint>0.) then
          ccl      = ccl      + 1.0
          qlintavl = qlintavl + qlint
          qlintmaxl = max(qlint,qlintmaxl)
        end if

        ! Calculate average RWP
        if (qrint>0.) then
          qrintavl = qrintavl + qrint
        end if

        ! Find cloud base variables
        do k=1,kmax
          if (ql0(i,j,k) > 0.) then
            zbase    = zf(k)
            zbaseavl = zbaseavl + zbase
            if (zbase < zbaseminl) zbaseminl = zbase
            exit
          end if
        end do

        ! Find cloud top variables
        do  k=1,kmax
          if (ql0(i,j,k) > 0) ztop = zf(k)
        end do

        if (ztop.ne.0) then  ! Non-cloudy columns not taken into account
           ztopavl = ztopavl + ztop
           if (ztop > ztopmaxl) ztopmaxl = ztop
        endif

        ! Find maximum vertical velocity and maximum ql in the domain
        wmaxl  = maxval( wm(i,j,1:kmax))
        qlmaxl = maxval(ql0(i,j,1:kmax))

        ! Calculate the vertical integral of the turbulent kinetic energy
        do  k=1,kmax
          tkeintl   = tkeintl + ( 0.5*( &
                           (0.5*(um(i,j,k)+um(i+1,j,k))+cu-u0av(k))**2 &
                          +(0.5*(vm(i,j,k)+vm(i,j+1,k))+cv-v0av(k))**2 &
                          +(0.5*(wm(i,j,k)+wm(i,j,k+1))           )**2 &
                                      ) + e12m(i,j,k)**2 ) * dzf(k)
        end do

      end do
    end do
    ! End do loop over the horizontal domain i,j
    
    ! Perform sums that don't require a do-loop
    call horAverage(ustar       ,ust)
    call horAverage(thlflux/ust ,tst)
    call horAverage(qtflux /ust ,qst)

    !+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
    ! Communication between CPU to find domain wide value for each variable
    call MPI_ALLREDUCE(ccl      , cc      , 1, MY_REAL, MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(qlintavl , qlintav , 1, MY_REAL, MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(qlintmaxl, qlintmax, 1, MY_REAL, MPI_MAX, comm3d,mpierr)
    call MPI_ALLREDUCE(qlint2l  , qlint2  , 1, MY_REAL, MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(qrintavl , qrintav , 1, MY_REAL, MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(zbaseavl , zbaseav , 1, MY_REAL, MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(zbaseminl, zbasemin, 1, MY_REAL, MPI_MIN, comm3d,mpierr)
    call MPI_ALLREDUCE(ztopavl  , ztopav  , 1, MY_REAL, MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(ztopmaxl , ztopmax , 1, MY_REAL, MPI_MAX, comm3d,mpierr)
    call MPI_ALLREDUCE(wmaxl    , wmax    , 1, MY_REAL, MPI_MAX, comm3d,mpierr)
    call MPI_ALLREDUCE(qlmaxl   , qlmax   , 1, MY_REAL, MPI_MAX, comm3d,mpierr)
    call MPI_ALLREDUCE(tkeintl  , tkeint  , 1, MY_REAL, MPI_MAX, comm3d,mpierr)
    !+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    ! Normalize fields if necessary
    if (cc > 0.0) then
      zbaseav = zbaseav/cc
      ztopav  = ztopav /cc
    else
      zbaseav   = missVal
      ztopav    = missVal
      zbasemin  = missVal
      ztopmax   = missVal
    end if
    
    cc      = cc      / rslabs
    qlintav = qlintav / rslabs
    qlint2  = qlint2  / rslabs
    qrintav = qrintav / rslabs
    tkeint  = tkeint  / rslabs
    
    ! Calculation of some extra variables, more convenient when at the end
    !+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
    qlintvar = (qlint2-(qlintav)**2)
    ! Surface fluxes of theta,thetav,qt; depends on the surface scheme used.
    if (isurf < 3) then
      call horAverage(thlflux, wts )
      call horAverage(qtflux , wqts)
    else
      wts  = wtsurf
      wqts = wqsurf
    end if
    ! Calculation of wtvs. NOTE: not valid in case of ql at surface
    c1   = 1.+ep2*qts
    c2   =    ep2*thls
    wtvs = c1*wts + c2*wqts

    ! Calculation of some variables specific for land-surface scheme
    if (isurf == 1) then
      call horAverage(Qnet     ,Qnetav    )
      call horAverage(H        ,Hav       )
      call horAverage(LE       ,LEav      )
      call horAverage(G0       ,G0av      )
      call horAverage(tendskin ,tendskinav)
      call horAverage(rs       ,rsav      )
      call horAverage(ra       ,raav      )
      call horAverage(tskin    ,tskinav   )
      call horAverage(cliq     ,cliqav    )
      call horAverage(wl       ,wlav      )
      call horAverage(rssoil   ,rssoilav  )
      call horAverage(rsveg    ,rsvegav   )
    end if
   
    ! Make sure the variables are in the same order as in ncName!
    n=1
    vars(n) = rtimee;       n=n+1
    vars(n) = cc;           n=n+1
    vars(n) = zbaseav;      n=n+1
    vars(n) = ztopav;       n=n+1
    vars(n) = ztopmax;      n=n+1
    vars(n) = zi;           n=n+1
    vars(n) = we;           n=n+1
    vars(n) = qlintav;      n=n+1
    vars(n) = qlintmax;     n=n+1
    vars(n) = wmax;         n=n+1
    vars(n) = tkeint;       n=n+1
    vars(n) = qlmax;        n=n+1
    vars(n) = ust;          n=n+1
    vars(n) = tst;          n=n+1
    vars(n) = qst;          n=n+1
    vars(n) = oblav;        n=n+1
    vars(n) = thls;         n=n+1
    vars(n) = z0;           n=n+1
    vars(n) = wts;          n=n+1
    vars(n) = wtvs;         n=n+1
    vars(n) = wqts;         n=n+1

    if (isurf == 1) then
      vars(n) = Qnetav;     n=n+1
      vars(n) = Hav;        n=n+1
      vars(n) = LEav;       n=n+1
      vars(n) = G0av;       n=n+1
      vars(n) = tendskinav; n=n+1
      vars(n) = rsav;       n=n+1
      vars(n) = raav;       n=n+1
      vars(n) = tskinav;    n=n+1
      vars(n) = cliqav;     n=n+1
      vars(n) = wlav;       n=n+1
      vars(n) = rssoilav;   n=n+1
      vars(n) = rsvegav;    n=n+1
    end if

  end subroutine calcTimestat

!===========================================================================================================!
! In this subroutine, user specific variables can be added to the timestat output. Three steps are          !
! necessary, namely: 1) adding the information to ncName by using ncinfo                                    !
!                    2) calculation of the variable                                                         !
!                    3) adding it to vars                                                                   !
! All these steps can be done locally, although it is also possible to do calculations in the calcTimeStat  !
! routine if that's more convenient. This increases the chance of merging warnings and errors however.      !
!===========================================================================================================!

  subroutine calcTimestatUser(n,ncName)
    use modstat_nc,  only           : ncinfo,lnetcdf
    use modglobal,   only           : kmax,i1,j1,rslabs
    use modsurfdata, only           : albedo
    use modraddata,  only           : swu,swd,swuToA,swdToA
    use modmicrodata,only           : qc,qr,precep,Nc,Nr,epscloud,epsqr,epsprec, &
                                      imicro, imicro_bulk
    use modmpi,      only           : MPI_INTEGER,MY_REAL,MPI_SUM,comm3d,mpierr,myid

    implicit none
    integer, intent(inout)         :: n
    character(80), dimension(:,:), &
      intent(inout), optional      :: ncName
    logical, save                  :: isInitTimestatUser = .false. ! Check for initialization of user timestat
    real                           :: prec_cntl          , &
                                      Nc_avl             , &
                                      Nr_avl                 
    integer                        :: cld_boxesl         , & 
                                      cld_boxes          , &
                                      rain_boxesl        , &
                                      rain_boxes
    real                           :: swuav              , &       ! Average shortwave up at domain top (kmax)
                                      swdav                        ! Average shortwave down at domain top (kmax)
    integer                        :: i,j,k
    
    != Add the names of the variables you want to add here
    if (lnetcdf .and. .not. isInitTimestatUser) then
      call ncinfo(ncName(n,:),'zb_min','Minimum cloud base height','m','time');                            n=n+1
      call ncinfo(ncName(n,:),'lwpvar','Variance of liquid water path','kg^2/m^4','time');                 n=n+1
      call ncinfo(ncName(n,:),'alb_sfc','Surface albedo','-','time');                                      n=n+1
      call ncinfo(ncName(n,:),'alb_tod','Top of LES domain shortwave albedo','-','time');                  n=n+1
      call ncinfo(ncName(n,:),'alb_toa','Top of atmosphere shortwave albedo','-','time');                  n=n+1
      call ncinfo(ncName(n,:),'rwp','Average rain water path','kg/m^2','time');                            n=n+1
      call ncinfo(ncName(n,:),'Nc','Mean cloud drop conc. (cloudy cells only)','m^-3','time');             n=n+1
      call ncinfo(ncName(n,:),'Nr','Mean rain drop conc. (rainy cells only)','m^-3','time');               n=n+1
      call ncinfo(ncName(n,:),'prec_srf','Mean precipitation flux at the surface','kg/kg m/s','time');     n=n+1
      call ncinfo(ncName(n,:),'prec_fracsrf','Fraction of boxes with precip/total at surface','-','time'); n=n+1
      isInitTimestatUser = .true.
      return
    end if

    != Calculate the variables that you want to add here, or somewhere in 'calcTimestat' if that is more convenient
    call horAverage( albedo, alb_sfc )
    
    ! Albedo at the top of the domain
    call horAverage(  swu(1:i1,1:j1,kmax), swuav  )
    call horAverage( -swd(1:i1,1:j1,kmax), swdav  )
    if (swdav > 0.) then
      alb_tod = swuav / swdav
    else
      alb_tod = missVal
    end if
    
    ! Albedo at the top of the atmosphere
    call horAverage(  swuToA(1:i1,1:j1), swuav  )
    call horAverage( -swdToA(1:i1,1:j1), swdav  )
    if (swdav > 0.) then
      alb_toa = swuav / swdav
    else
      alb_toa = missVal
    end if

    ! Some microphysics variables
    if (imicro == imicro_bulk) then
      do k=1,kmax
      do j=2,j1
      do i=2,i1

        if ( qc(i,j,k) > epscloud ) then
          Nc_avl     = Nc_avl + Nc(i,j,k)
          cld_boxesl = cld_boxesl + 1
        end if

        if ( qr(i,j,k) > epsqr ) then
          Nr_avl      = Nr_avl + Nr(i,j,k)
          rain_boxesl = rain_boxesl + 1
        end if

      end do
      end do
      end do

      call MPI_ALLREDUCE(cld_boxesl , cld_boxes , 1, MPI_INTEGER, MPI_SUM, comm3d, mpierr)
      call MPI_ALLREDUCE(rain_boxesl, rain_boxes, 1, MPI_INTEGER, MPI_SUM, comm3d, mpierr)
      call MPI_ALLREDUCE(Nc_avl     , Nc_av     , 1, MY_REAL    , MPI_SUM, comm3d, mpierr)
      call MPI_ALLREDUCE(Nr_avl     , Nr_av     , 1, MY_REAL    , MPI_SUM, comm3d, mpierr)

      ! Failsafe divide by zero Nc; -1. is converted to missing values
      if (cld_boxes == 0) then
        Nc_av  = missVal
      else
        Nc_av  = Nc_av / cld_boxes
      end if

      ! Failsafe divide by zero Nr; -1. is converted to missing values
      if (rain_boxes == 0) then
        Nr_av  = missVal
      else
        Nr_av  = Nr_av / rain_boxes
      end if

      ! Calculate surface precipitation and fraction of gridboxes with precipitation
      call horAverage( precep(:,:,1) , prec_srf      )
      
      prec_cntl       = count ( precep  (2:i1,2:j1,1) > epsprec )
      call MPI_ALLREDUCE(prec_cntl, prec_fracsrf, 1, MY_REAL , MPI_SUM, comm3d, mpierr)
      prec_fracsrf   = prec_fracsrf / rslabs

    end if ! imicro == imicrobulk
    
    != Put the variables in 'vars', which is written to the file later
    vars(n) = zbasemin;     n=n+1
    vars(n) = qlintvar;     n=n+1
    vars(n) = alb_sfc;      n=n+1
    vars(n) = alb_tod;      n=n+1
    vars(n) = alb_toa;      n=n+1
    vars(n) = qrintav;      n=n+1
    vars(n) = Nc_av;        n=n+1
    vars(n) = Nr_av;        n=n+1
    vars(n) = prec_srf;     n=n+1
    vars(n) = prec_fracsrf; n=n+1

  end subroutine calcTimestatUser

!====================================================================================================;

  subroutine timestatNC
    use modsurfdata, only : isurf,thls,z0
    use modglobal,   only : cexpnr,rtimee
    use modstat_nc,  only : ncinfo,open_nc,define_nc,writestat_nc
    use modmpi,      only : myid
    implicit none
    
    character(80)                  :: fname = 'tmser.xxx.nc'   ! Name of the NetCDF output file
    character(80),allocatable,  &                              ! Variable that holds: name, description,
                   dimension(:,:)  :: ncName                   ! unit and dimension of variables
    logical, save                  :: isInitNCOut = .false.    ! Check for initialization of NC output
    integer                        :: n=1                  , & ! Variable number
                                      nVar                     ! Total number of variables

    if (.not.isInitNCOut) then
      allocate(ncName(nVarMax,4))
      call ncinfo(ncName(n,:),'time','Time','s','time');                                                   n=n+1
      call ncinfo(ncName(n,:),'cfrac','Cloud fraction','-','time');                                        n=n+1
      call ncinfo(ncName(n,:),'zb','Cloud-base height','m','time');                                        n=n+1
      call ncinfo(ncName(n,:),'zc_av','Average Cloud-top height','m','time');                              n=n+1
      call ncinfo(ncName(n,:),'zc_max','Maximum Cloud-top height','m','time');                             n=n+1
      call ncinfo(ncName(n,:),'zi','Boundary layer height','m','time');                                    n=n+1
      call ncinfo(ncName(n,:),'we','Entrainment velocity','m/s','time');                                   n=n+1
      call ncinfo(ncName(n,:),'lwp_bar','Average liquid-water path','kg/m^2','time');                      n=n+1
      call ncinfo(ncName(n,:),'lwp_max','Maximum liquid-water path','kg/m^2','time');                      n=n+1
      call ncinfo(ncName(n,:),'wmax','Maximum vertical velocity','m/s','time');                            n=n+1
      call ncinfo(ncName(n,:),'vtke','Vertical integral of total TKE','m^3/s^2','time');                   n=n+1
      call ncinfo(ncName(n,:),'lmax','Maximum liquid water mixing ratio','kg/kg','time');                  n=n+1
      call ncinfo(ncName(n,:),'ustar','Surface friction velocity','m/s','time');                           n=n+1
      call ncinfo(ncName(n,:),'tstr','Turbulent temperature scale','K','time');                            n=n+1
      call ncinfo(ncName(n,:),'qtstr','Turbulent humidity scale','K','time');                              n=n+1
      call ncinfo(ncName(n,:),'obukh','Obukhov Length','m','time');                                        n=n+1
      call ncinfo(ncName(n,:),'thlskin','Surface liquid water potential temperature','K','time');          n=n+1
      call ncinfo(ncName(n,:),'z0','Roughness height','m','time');                                         n=n+1
      call ncinfo(ncName(n,:),'wtheta','Surface kinematic temperature flux','K m/s','time');               n=n+1
      call ncinfo(ncName(n,:),'wthetav','Surface kinematic virtual temperature flux','K m/s','time');      n=n+1
      call ncinfo(ncName(n,:),'wq','Surface kinematic moisture flux','kg/kg m/s','time');                  n=n+1

      if (isurf==1) then
        call ncinfo(ncName(n,:),'Qnet','Net radiation','W/m^2','time');                                    n=n+1
        call ncinfo(ncName(n,:),'H','Sensible heat flux','W/m^2','time');                                  n=n+1
        call ncinfo(ncName(n,:),'LE','Latent heat flux','W/m^2','time');                                   n=n+1
        call ncinfo(ncName(n,:),'G0','Ground heat flux','W/m^2','time');                                   n=n+1
        call ncinfo(ncName(n,:),'tendskin','Skin tendency','W/m^2','time');                                n=n+1
        call ncinfo(ncName(n,:),'rs','Surface resistance','s/m','time');                                   n=n+1
        call ncinfo(ncName(n,:),'ra','Aerodynamic resistance','s/m','time');                               n=n+1
        call ncinfo(ncName(n,:),'tskin','Skin temperature','W/m^2','time');                                n=n+1
        call ncinfo(ncName(n,:),'cliq','Fraction of vegetated sfc covered with liquid water','-','time');  n=n+1
        call ncinfo(ncName(n,:),'Wl','Liquid water reservoir','m','time');                                 n=n+1
        call ncinfo(ncName(n,:),'rssoil','Soil evaporation resistance','s/m','time');                      n=n+1
        call ncinfo(ncName(n,:),'rsveg','Vegitation resistance','s/m','time');                             n=n+1
      end if

      call calcTimestatUser(n,ncName=ncName)

      nVar = n - 1

      fname(7:9) = cexpnr
      if (myid==0) then
        call open_nc(fname, ncid, nRec)
        call define_nc(ncid, nVar, ncName(1:nVar,:))
      end if

      isInitNCOut = .true.
      return
    end if
    
    if (myid == 0) then
      call writestat_nc(ncid,nVar,ncName,vars(1:nVar), nRec,.true.)
    end if

  end subroutine timestatNC

!====================================================================================================;

  subroutine timestatAscii
    use modglobal,   only : ifoutput,cexpnr,rtimee
    use modsurfdata, only : isurf
    use modmpi,      only : myid
    implicit none

    logical, save         :: isInitAsciiOut = .false. ! Check for initialization of ASCII output

    if (.not.isInitAsciiOut) then

      ! Write headers for the ASCII output files. Use only 1 CPU: myid=0
      if(myid==0) then

        ! Header tmser1.[expnr]
        open (ifoutput,file='tmser1.'//cexpnr,status='replace',position='append')
        write(ifoutput,'(2a)') &
               '#  time      cc     z_cbase_avg  z_ctop_avg  z_ctop_max    ', &
               ' zi       we    <<ql>>  <<ql>>_max   w_max   tke     ql_max'
        close(ifoutput)

        ! Header tmsurf.[expnr]
        open (ifoutput,file='tmsurf.'//cexpnr,status='replace',position='append')
        write(ifoutput,'(2a)') &
               '#  time        ust        tst        qst         obukh     ', &
               ' thls        z0        wthls      wthvs      wqts '
        close(ifoutput)

        if(isurf == 1) then
          ! Header tmlsm.[expnr]
          open (ifoutput,file='tmlsm.'//cexpnr,status='replace',position='append')
          write(ifoutput,'(3a)') &
                 '#     time      Qnet        H          LE         G0     ', &
                 'tendskin     rs         ra        tskin        cliq      ', &
                 'Wl          rssoil     rsveg'
          write(ifoutput,'(3a)') &
                 '#      [s]     [W/m2]     [W/m2]     [W/m2]     [W/m2]   ', &
                 '[W/m2]      [s/m]       [s/m]     [K]          [-]      ' , &
                 '[m]          [s/m]      [s/m]'
          close(ifoutput)
        end if

      end if ! myid == 0

      isInitAsciiOut = .true.
      return
    end if

    if (myid == 0) then
      open (ifoutput,file='tmser1.'//cexpnr,position='append')
      write( ifoutput,'(f10.2,f5.2,4f12.3,f10.4,5f9.3)') &
          vars(1:7)      , &
          vars(8:9)*1000., &
          vars(10:11)    , &
          vars(12)*1000.
      close(ifoutput)

      open (ifoutput,file='tmsurf.'//cexpnr,position='append')
      write( ifoutput,'(f10.2,3e11.3,2f11.3,4e11.3)   ') &
          vars(1)        , &
          vars(13:21)
      close(ifoutput)

      if (isurf == 1) then
        open (ifoutput,file='tmlsm.'//cexpnr,position='append')
        write(ifoutput,'(f10.2,9f11.3,e13.3, 2f11.3)') &
            vars(1)      , &
            vars(22:33)
        close(ifoutput)
      end if

    end if

  end subroutine timestatAscii

!=================================================================================================================;
! Subroutine that calculates the boundary layer height and the entrainment                                        ;
!=================================================================================================================;

  subroutine calcBLHeight
    use modglobal,  only            : ih,i1,jh,j1,kmax,k1,cp,rlv,imax,rd,zh,dzh,zf,dzf,rv,rslabs,iadv_sv,iadv_kappa
    use modfields,  only            : thlprof,qtprof,svprof,w0,qt0,qt0h,ql0,thl0,thl0h,thv0h,sv0,exnf,whls
    use modsurfdata,only            : svs
    use modmpi,     only            : myid,MPI_SUM,MY_REAL,comm3d,mpierr
    implicit none

    logical, save                  :: isInitBLHeight = .false. ! Check for initialization of BL height calculations
    real,allocatable,dimension(:,:,:) :: sv0h,  &              ! Half-level concentration of passive scalar
                                      blh_fld                  ! Field that is used to calculate Bl height
    real,allocatable,dimension(:)  :: profile,  &              ! Temporary variable to store qt/thl/thv profiles
                                      dprof,    &              ! First order derivative (wrt height) of profile
                                      ddprof                   ! Second order derivative (wrt height) of profile                    
    real                           :: zil,      &              ! Inversion height (CPU local value)
                                      ziold=-1, &              ! Temporary value of the average inversion height
                                      gradient, &              ! Used to find the maximum gradient (initialization)
                                      dhdt,     &              ! Change in inversion height (with time)
                                      locval,   &              ! -
                                      oldlocval,&              ! -
                                      blh_sign                 ! -1/+1, to use with gradient
    integer                        :: location, &              ! Location of the strongest gradient in prof
                                      nsamp,    &              ! Number of samples ??
                                      stride,   &              ! ??
                                      i,j,k

    if (.not. isInitBLHeight) then
      allocate(profile(1:k1))

      select case (iblh_var)
      case(iblh_qt)
        profile = qtprof
      case(iblh_thl)
        profile = thlprof
      case(iblh_thv)
        do k=1,kmax
          profile(k) = thlprof(k)*(1+(rv/rd-1)*qtprof(k))
        end do
      case(1:)
        profile = svprof(:,iblh_var)
      end select

      blh_sign = sign(1.0,profile(kmax)-profile(1))

      select case(iblh_meth)
      case (iblh_flux)
      case (iblh_grad)
      case (iblh_thres)
        if (blh_thres<0) then
          do k=kmax,2,-1
            if (blh_sign*(profile(k+1) - profile(k-1)) > gradient) then
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

      deallocate(profile)
      isInitBLHeight = .true.

      return
    end if

    != End of initialization, start actual calculations =!

    allocate(blh_fld(2-ih:i1+ih,2-jh:j1+jh,k1), &
             sv0h   (2-ih:i1+ih,2-jh:j1+jh,k1)  )
    allocate(profile (k1),                      &
             dprof(k1)   ,                      &
             ddprof(k1)                         )

    zil      = 0.
    dprof    = 0.
    ddprof   = 0.
    select case (iblh_meth)
    case (iblh_flux)
      select case (iblh_var)
      case(iblh_qt)
        blh_fld = w0*qt0h
      case(iblh_thl)
        blh_fld = w0*thl0h
      case(iblh_thv)
        blh_fld = w0*thv0h
      case(1:)
        if (iadv_sv(iblh_var)==iadv_kappa) then
          call halflev_kappa(sv0(:,:,:,iblh_var),sv0h)
        else
          do  j=2,j1
          do  i=2,i1
          do  k=2,kmax
            sv0h(i,j,k) = (sv0(i,j,k,iblh_var)*dzf(k-1)+sv0(i,j,k-1,iblh_var)*dzf(k))/(2*dzh(k))
          enddo
          enddo
          enddo
          sv0h(2:i1,2:j1,1) = svs(iblh_var)
          blh_fld = w0*sv0h
        end if
      end select
    case (iblh_grad,iblh_thres)
      select case (iblh_var)
      case(iblh_qt)
        blh_fld = qt0
      case(iblh_thl)
        blh_fld = thl0
      case(iblh_thv)
        do k=1,kmax
          blh_fld(2:i1,2:j1,k) = (thl0(2:i1,2:j1,k)+rlv*ql0(2:i1,2:j1,k)/(cp*exnf(k))) &
                      *(1+(rv/rd-1)*qt0(2:i1,2:j1,k)-rv/rd*ql0(2:i1,2:j1,k))
        end do
      case(1:)
        blh_fld = sv0(2:i1,2:j1,1:k1,iblh_var)
      end select
    end select

    select case (iblh_meth)
    case (iblh_flux)
      stride = ceiling(real(imax)/real(blh_nsamp))
      do i=2,stride+1
        nsamp =  ceiling(real(i1-i+1)/real(stride))
        do j=2,j1
          zil = zil + nsamp*zh(minloc(sum(blh_fld(i:i1:stride,j,:),1),1))
        end do
      end do
    case (iblh_grad)
      stride = ceiling(real(imax)/real(blh_nsamp))
      do i=2,stride+1
        nsamp =  ceiling(real(i1-i+1)/real(stride))
        do j=2,j1
          profile  = sum(blh_fld(i:i1:stride,j,:),1)
          select case (iblh_var)
          case(iblh_qt) !Water vapour gradients near the inversion layer can be either positive or negative
            dprof(2:kmax) = abs(profile(2:kmax) - profile(1:kmax-1))/dzh(2:kmax)
          case(iblh_thl,iblh_thv) !temperature jumps near the inversion layer are always positive
            dprof(2:kmax) = (profile(2:kmax) - profile(1:kmax-1))/dzh(2:kmax)
          case default
            dprof(2:kmax) = (profile(2:kmax) - profile(1:kmax-1))/dzh(2:kmax)
          end select
          ddprof(2:kmax-1)    = (dprof(3:kmax) - dprof(2:kmax-1))/dzf(2:kmax-1)
          location = maxloc(dprof,1)
          zil  = zil + nsamp*(zh(location-1) - dzh(location)*ddprof(location-1)/(ddprof(location)-ddprof(location-1) + 1.e-8))
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
    call MPI_ALLREDUCE(zil, zi, 1, MY_REAL, MPI_SUM, comm3d,mpierr)
    zi = zi / rslabs

    if (ziold < 0) ziold = zi                     ! Set initial value
    dhdt = (zi-ziold)/dtav
    ziold = zi

    k=2
    do while (zh(k)<zi .and. k < kmax)
      k=k+1
    end do
    we = dhdt - whls (k)   !include for large-scale vertical velocity
    deallocate(blh_fld,sv0h        )
    deallocate(profile,dprof,ddprof)

  end subroutine calcBLHeight

!====================================================================================================;
! This subroutine is used to calculate a domain average for a variable

  subroutine horAverage(varl,varav)
    use modglobal, only : rslabs,i1,j1
    use modmpi,    only : MY_REAL,MPI_SUM,comm3d,mpierr
    implicit none

    real, dimension(i1,j1), intent(in) :: varl
    real, intent(out)                  :: varav
    real                               :: varavl
    varav = 0.; varavl = 0.

    varavl = sum(varl(2:i1,2:j1))
    call MPI_ALLREDUCE(varavl, varav, 1, MY_REAL,MPI_SUM,comm3d,mpierr)
    varav  = varav/rslabs

  end subroutine horAverage

!====================================================================================================;

  subroutine exittimestat
    use modmpi,     only : myid
    use modstat_nc, only : exitstat_nc,lnetcdf
    implicit none

    if(ltimestat .and. lnetcdf .and. myid==0) call exitstat_nc(ncid)

  end subroutine exittimestat

end module modtimestat
