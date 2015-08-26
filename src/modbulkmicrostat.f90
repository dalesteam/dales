!> \file modbulkmicrostat.f90
!!  Calculates profiles coming from the bulkmicrophysics


!>
!!  Calculates profiles coming from the bulkmicrophysics
!>
!! Profiles coming from the bulkmicrophysics. Written to precep.expnr for the
!! rain rates etc., and to qlptend.expnr, nptend.expnr and qtptend.expnr for the
!! tendencies is rain water content, droplet number, and total water content,
!! respectively.
!! If netcdf is true, this module also writes in the profiles.expnr.nc output
!!  \author Olivier Geoffroy, KNMI
!!  \author Johan van de Dussen, TU Delft
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
module modbulkmicrostat
  use modglobal, only : longint

implicit none
private
PUBLIC  :: initbulkmicrostat, bulkmicrostat, exitbulkmicrostat,bulkmicrotend,lmicrostat
save
!NetCDF variables
  integer,parameter :: nvar = 26
  integer :: ncid,nrec = 0
  character(80) :: fname = 'bulkmicrostat.xxx.nc'
  character(80),dimension(nvar,4) :: ncname
  character(80),dimension(1,4) :: tncname

  real          :: dtav, timeav
  integer(kind=longint):: idtav, itimeav, tnext, tnextwrite
  integer          :: nsamples
  logical          :: lmicrostat = .false.
  integer, parameter      :: nrfields = 5  , &
                             iauto    = 2  , &
                             iaccr    = 3  , &
                             ievap    = 4  , &
                             ised     = 5
  real, allocatable, dimension(:,:)  :: Npav    , &
                                        Npmn    , &
                                        qlpav  , &
                                        qlpmn  , &
                                        qtpav  , &
                                        qtpmn
  real, allocatable, dimension(:)    :: precavl  , &
               precav  , &
               precmn  , &
               sedavl , &
               sedav  ,&
               sedmn  , &
               preccountavl  , &
               preccountav  , &
               preccountmn  , &
               prec_prcavl  , &
               prec_prcav  , &
               prec_prcmn  , &
               prec_fracavl  , &
               prec_fracav  , &
               prec_fracmn  , &
               cloudcountavl, &
               cloudcountav  , &
               cloudcountmn  , &
               raincountavl  , &
               raincountav  , &
               raincountmn  , &
               Nrrainavl  , &
               Nrrainav  , &
               Nrrainmn  , &
               Nrcloudavl  , &
               Nrcloudav  , &
               Nrcloudmn  , &
               qravl  , &
               qrav    , &
               qrmn    , &
               Dvravl  , &
               Dvrav  , &
               Dvrmn

contains
!> Initialization routine, reads namelists and inits variables
subroutine initbulkmicrostat
    use modmpi,    only  : myid, mpi_logical, my_real, comm3d, mpierr
    use modglobal, only  : ifnamopt, fname_options, cexpnr, ifoutput, &
              dtav_glob, timeav_glob, ladaptive, k1,kmax, dtmax,btime,tres,&
                           i1,j1,k1,ih,jh
    use modstat_nc, only : lnetcdf, open_nc,define_nc,redefine_nc,ncinfo,writestat_dims_nc
    use modgenstat, only : idtav_prof=>idtav, itimeav_prof=>itimeav,ncid_prof=>ncid
    use modmicrodata,only: imicro, imicro_bulk,Dvr,precep,qrevap,qrsed,qrsrc,l_rain
    implicit none
    integer      :: ierr

    namelist/NAMBULKMICROSTAT/ &
    lmicrostat, dtav, timeav

    !if (imicro /=imicro_bulk) return

    dtav  = dtav_glob
    timeav  = timeav_glob
    if(myid==0)then
      open (ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMBULKMICROSTAT,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions NAMBULKMICROSTAT'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMBULKMICROSTAT'
      endif
      write(6,NAMBULKMICROSTAT)
      close(ifnamopt)
    end if

    call MPI_BCAST(lmicrostat  ,1,MPI_LOGICAL  ,0,comm3d,mpierr)
    call MPI_BCAST(dtav    ,1,MY_REAL  ,0,comm3d,mpierr)
    call MPI_BCAST(timeav    ,1,MY_REAL  ,0,comm3d,mpierr)
    idtav = dtav/tres
    itimeav = timeav/tres

    tnext      = idtav   +btime
    tnextwrite = itimeav +btime
    nsamples = itimeav/idtav

    if (l_rain .or. lmicrostat) then
      allocate( Dvr   (2-ih:i1+ih,2-jh:j1+jh,k1), &
                precep(2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate( qrevap(k1), &
                qrsed(k1) , &
                qrsrc(k1) )
      Dvr=0.; precep=0.
      qrevap=0.; qrsed=0.; qrsrc=0.
    end if
    
    if (.not. lmicrostat) return
    if (abs(timeav/dtav - nsamples) > 1e-4) then
      stop 'timeav must be an integer multiple of dtav (NAMBULKMICROSTAT)'
    end if
    if (.not. ladaptive .and. abs(dtav/dtmax - nint(dtav/dtmax)) > 1e-4) then
      stop 'dtav must be an integer multiple of dtmax (NAMBULKMICROSTAT)'
    end if

    allocate(Npav    (k1, nrfields)  , &
       Npmn    (k1, nrfields)  , &
       qlpav    (k1, nrfields)  , &
       qlpmn    (k1, nrfields)  , &
       qtpav    (k1, nrfields)  , &
       qtpmn    (k1, nrfields)  )
    allocate(precavl  (k1)    , &
       precav    (k1)    , &
       precmn    (k1)    , &
       sedavl (k1)      , &
       sedav (k1)       ,&
       sedmn    (k1)    , &
       preccountavl  (k1)    , &
       preccountav  (k1)    , &
       preccountmn  (k1)    , &
       prec_prcavl  (k1)    , &
       prec_prcav  (k1)    , &
       prec_prcmn  (k1)    , &
       prec_fracavl  (k1)    , &
       prec_fracav  (k1)    , &
       prec_fracmn  (k1)    , &
       cloudcountavl  (k1)    , &
       cloudcountav  (k1)    , &
       cloudcountmn  (k1)    , &
       raincountavl  (k1)    , &
       raincountav  (k1)    , &
       raincountmn  (k1)    , &
       Nrrainavl  (k1)    , &
       Nrrainav  (k1)    , &
       Nrrainmn  (k1)    , &
       Nrcloudavl  (k1)    , &
       Nrcloudav  (k1)    , &
       Nrcloudmn  (k1)    , &
       qravl    (k1)    , &
       qrav    (k1)    , &
       qrmn    (k1)    , &
       Dvravl    (k1)    , &
       Dvrav    (k1)    , &
       Dvrmn    (k1)    )
    Npmn    = 0.0
    qlpmn    = 0.0
    qtpmn    = 0.0
    precmn    = 0.0
    sedmn    = 0.0
    preccountmn  = 0.0
    prec_prcmn  = 0.0
    prec_fracmn  = 0.0
    cloudcountmn  = 0.0
    raincountmn  = 0.0
    Nrrainmn  = 0.0
    Nrcloudmn  = 0.0
    qrmn    = 0.0
    Dvrmn    = 0.0

    if (myid == 0) then
      open (ifoutput,file = 'precep.'//cexpnr ,status = 'replace')
      close(ifoutput)
      open (ifoutput,file = 'nptend.'//cexpnr ,status = 'replace')
      close(ifoutput)
      open (ifoutput,file = 'qlptend.'//cexpnr,status = 'replace')
      close(ifoutput)
      open (ifoutput,file = 'qtptend.'//cexpnr,status = 'replace')
      close(ifoutput)
    end if
    if (lnetcdf) then
      idtav = idtav_prof
      itimeav = itimeav_prof
      tnext      = idtav+btime
      tnextwrite = itimeav+btime
      nsamples = itimeav/idtav
      if (myid==0) then
        fname(15:17) = cexpnr
        call ncinfo(tncname(1,:),'time','Time','s','time')
        call ncinfo(ncname( 1,:),'qr','Average rain water content','kg/kg','tt')
        call ncinfo(ncname( 2,:),'cfrac','Cloud fraction','-','tt')
        call ncinfo(ncname( 3,:),'raincount','Rainy gridcell fraction','-','tt')
        call ncinfo(ncname( 4,:),'preccount','Fraction of gridcells with precip>eps_prec','-','mt')
        call ncinfo(ncname( 5,:),'prec_frac','Fraction of precip. cells with precip>epsprec','-','tt')
        call ncinfo(ncname( 6,:),'prec','Average precipitation flux','kg/kg m/s','mt')
        call ncinfo(ncname( 7,:),'sed','Average sedimentation flux','kg/kg m/s','mt')
        call ncinfo(ncname( 8,:),'prec_prc','Average precipitation flux (prec. gridboxes only)','kg/kg m/s','tt')
        call ncinfo(ncname( 9,:),'dvrmn','Average precip. drop diameter (rainy cells only)','m','tt')
        call ncinfo(ncname(10,:),'nrcloud','Cloud droplet number density (cloudy cells only)','m^-3','mt')
        call ncinfo(ncname(11,:),'nrrain','Rain droplet number density (rainy cells only)','m^-3','mt')

        call ncinfo(ncname(12,:),'npauto','Autoconversion rain drop tendency','#/m3/s','tt')
        call ncinfo(ncname(13,:),'npaccr','Accretion rain drop tendency','#/m3/s','tt')
        call ncinfo(ncname(14,:),'npsed','Sedimentation rain drop tendency','#/m3/s','tt')
        call ncinfo(ncname(15,:),'npevap','Evaporation rain drop tendency','#/m3/s','tt')
        call ncinfo(ncname(16,:),'nptot','Total rain drop tendency','#/m3/s','tt')
        call ncinfo(ncname(17,:),'qrpauto','Autoconversion rain water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(18,:),'qrpaccr','Accretion rain water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(19,:),'qrpsed','Sedimentation rain water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(20,:),'qrpevap','Evaporation rain water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(21,:),'qrptot','Total rain water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(22,:),'qtpauto','Autoconversion total water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(23,:),'qtpaccr','Accretion total water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(24,:),'qtpsed','Sedimentation total water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(25,:),'qtpevap','Evaporation total water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(26,:),'qtptot','Total total water content tendency','kg/kg/s','tt')
!        call open_nc(fname,ncid,n3=kmax)
!        call define_nc( ncid, 1, tncname)
!        call writestat_dims_nc(ncid)
!        call redefine_nc(ncid_prof)
        call define_nc( ncid_prof, NVar, ncname)
      end if

   end if

  end subroutine initbulkmicrostat

!------------------------------------------------------------------------------!
!> General routine, does the timekeeping
  subroutine bulkmicrostat
    use modmpi,    only  : myid
    use modglobal, only  : rkStep,rkMaxStep, timee, dt_lim
    implicit none
    if (.not. lmicrostat)  return
    if (rkStep /= rkMaxStep)  return
    if (timee == 0)    return
    if (timee < tnext .and. timee < tnextwrite) then
      dt_lim  = minval((/dt_lim, tnext - timee, tnextwrite - timee/))
      return
    end if
    if (timee >= tnext) then
      tnext = tnext + idtav
      call dobulkmicrostat
    end if
    if (timee >= tnextwrite) then
      tnextwrite = tnextwrite + itimeav
      call writebulkmicrostat
    end if

  end subroutine bulkmicrostat

!------------------------------------------------------------------------------!
!> Performs the calculations for rainrate etc.
  subroutine dobulkmicrostat
    use modmpi,        only  : my_real, mpi_sum, comm3d, mpierr
    use modglobal,     only  : i1, j1, k1, rslabs
    use modmicrodata,  only  : precep,sedc,l_sedc,Dvr,epscloud,epsqr,epsprec,iqr,iNr,Nc_0
    use modfields,     only  : ql0,sv0
    implicit none

    integer      :: k

    cloudcountav = 0.0
    raincountav  = 0.0
    preccountav  = 0.0
    prec_fracav  = 0.0
    prec_prcav   = 0.0
    Dvrav        = 0.0
    Nrcloudav    = 0.0
    Nrrainav     = 0.0
    precav       = 0.0
    sedav        = 0.0
    qrav         = 0.0

    do k = 1,k1
      cloudcountavl(k)  = count(ql0     (2:i1,2:j1,k) > epscloud)    ! Number of cloudy gridboxes
      raincountavl (k)  = count(sv0 (2:i1,2:j1,k,iqr) > epsqr   )    ! Number of rainy gridboxes
      preccountavl (k)  = count(precep  (2:i1,2:j1,k) > epsprec )    ! Number of precipitating gridboxes > epsprec
      prec_fracavl (k)  = count(precep  (2:i1,2:j1,k) > 1e-8    )    ! Number of gridboxes with precip > 1e-8
      prec_prcavl  (k)  = sum  (precep  (2:i1,2:j1,k) , precep(2:i1,2:j1,k) > epsprec ) ! Sum of precip fluxes for boxes with precip > eps 
      Dvravl       (k)  = sum  (Dvr     (2:i1,2:j1,k) , sv0(2:i1,2:j1,k,iqr) > epsqr   ) ! Sum of prec water mean diam for boxes with qr > eps
      Nrcloudavl   (k)  = cloudcountavl(k)*Nc_0 ! Sum of cloud drop number dens for boxes with ql0 > eps
      Nrrainavl    (k)  = sum  (sv0 (2:i1,2:j1,k,iNr) , sv0(2:i1,2:j1,k,iqr) > epsqr) ! Sum of rain drop number dens for boxes with qr > eps
      precavl      (k)  = sum  (precep  (2:i1,2:j1,k))               ! Total sum of precipitation fluxes
      if (l_sedc) then
        sedavl       (k)  = sum  (sedc    (2:i1,2:j1,k))               ! Total sum of cloud water sedimentation fluxes
      end if
      qravl        (k)  = sum  (sv0 (2:i1,2:j1,k,iqr))               ! Total sum of rain water content
    end do

    call MPI_ALLREDUCE(cloudcountavl,cloudcountav ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(raincountavl ,raincountav  ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(preccountavl ,preccountav  ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(prec_fracavl ,prec_fracav  ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(prec_prcavl  ,prec_prcav   ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(Dvravl       ,Dvrav        ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(Nrcloudavl   ,Nrcloudav    ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(Nrrainavl    ,Nrrainav     ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(precavl      ,precav       ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(sedavl       ,sedav        ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(qravl        ,qrav         ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)

    cloudcountmn = cloudcountmn +  cloudcountav /rslabs   ! fraction of cloudy gridcells
    raincountmn  = raincountmn  +  raincountav  /rslabs   ! fraction of rainy gridcells 
    preccountmn  = preccountmn  +  preccountav  /rslabs   ! fraction of precipitating gridcells
    ! Average prec_prc, Nrcloud and Nrrain over precip/cloudy/rainy gridcells only !
    do k=1,k1
      if (cloudcountav(k) > 0) then
        Nrcloudmn(k)  = Nrcloudmn(k) + Nrcloudav(k) / cloudcountav(k)
      end if
      if (raincountav(k) > 0) then
        Nrrainmn(k)   = Nrrainmn(k) + Nrrainav(k) / raincountav(k)
        Dvrmn(k)      = Dvrmn(k)    + Dvrav(k)    / raincountav(k)
      end if
      if (preccountav(k) > 0) then
        prec_prcmn(k) = prec_prcmn(k) +  prec_prcav(k) / preccountav(k)
      end if
      if (prec_fracav(k) > 0) then
        ! fraction of precipitating cells with prec >= epsprec
        prec_fracmn(k)  = prec_fracmn(k)  +  preccountav(k) / prec_fracav(k)
      end if
    end do
    precmn       = precmn       +  precav   /rslabs
    sedmn        = sedmn        +  sedav    /rslabs
    qrmn         = qrmn         +  qrav     /rslabs

  end subroutine dobulkmicrostat

!------------------------------------------------------------------------------!
!> Performs the calculations for the tendencies etc.
  subroutine bulkmicrotend
    use modmpi,    only  : slabsum
    use modglobal,    only  : rkStep,rkMaxStep, timee, dt_lim, k1, ih, i1, jh, j1, rslabs
    use modfields,    only : qtp
    use modmicrodata,  only  : qrp, Nrp,l_sedc
    implicit none

    real, dimension(:), allocatable  :: avfield
    integer        :: ifield = 0

    if (.not. lmicrostat)  return
    if (rkStep /= rkMaxStep)  return
    if (timee == 0)    return
    if (timee < tnext .and. timee < tnextwrite) then
      dt_lim  = minval((/dt_lim, tnext - timee, tnextwrite - timee/))
      return
    end if
!    tnext = tnext+dtav

    allocate(avfield(k1))

    ifield    = mod(ifield, nrfields) + 1

    avfield    = 0.0
    call slabsum(avfield  ,1,k1,Nrp  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    Npav(:,ifield)  = avfield - sum(Npav  (:,1:ifield-1),2)

    avfield    = 0.0
    call slabsum(avfield  ,1,k1,qrp  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    qlpav(:,ifield) = avfield - sum(qlpav  (:,1:ifield-1),2)

    avfield    = 0.0
    call slabsum(avfield  ,1,k1,-qrp ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    qtpav(:,ifield) = avfield - sum(qtpav  (:,1:ifield-1),2)

    if (ifield == nrfields) then
      Npmn    = Npmn  + Npav  /nsamples/rslabs
      qlpmn    = qlpmn  + qlpav /nsamples/rslabs
      qtpmn    = qtpmn  + qtpav /nsamples/rslabs
      Npav    = 0.0
      qlpav    = 0.0
      qtpav    = 0.0
    end if

    deallocate(avfield)
  end subroutine bulkmicrotend

!------------------------------------------------------------------------------!
!> Write the stats to file
  subroutine writebulkmicrostat
    use modmpi,    only  : myid
    use modglobal,    only  : rtimee, ifoutput, cexpnr, k1,kmax, &
              rlv, zf,rslabs
    use modfields,    only  : presf,rhof
    use modstat_nc, only: lnetcdf, writestat_nc
    use modgenstat, only: ncid_prof=>ncid,nrec_prof=>nrec
    use modmicrodata, only : qrsrc,qrsed,qrevap,tint

    implicit none
    real,dimension(k1,nvar) :: vars

    integer    :: nsecs, nhrs, nminut
    integer    :: k

    nsecs    = nint(rtimee)
    nhrs    = int (nsecs/3600)
    nminut    = int (nsecs/60)-nhrs*60
    nsecs    = mod (nsecs,60)

    cloudcountmn    = cloudcountmn  /nsamples
    raincountmn     = raincountmn   /nsamples
    preccountmn     = preccountmn   /nsamples
    prec_fracmn     = prec_fracmn   /nsamples
    prec_prcmn      = prec_prcmn    /nsamples
    Dvrmn           = Dvrmn         /nsamples
    Nrrainmn        = Nrrainmn      /nsamples
    Nrcloudmn       = Nrcloudmn     /nsamples
    precmn          = precmn        /nsamples
    sedmn           = sedmn         /nsamples
    qrmn            = qrmn          /nsamples

!!! Already done !!!
!    where (raincountmn > 0.)
!      Dvrmn        = Dvrmn / raincountmn
!    elsewhere
!      Dvrmn = 0.0
!    end where
!    where (preccountmn > 0.)
!      prec_prcmn = prec_prcmn/preccountmn
!    elsewhere
!      prec_prcmn = 0.0
!    end where

    if (myid == 0) then
    open (ifoutput,file='precep.'//cexpnr,position='append')
    write(ifoutput,'(//2A,/A,F5.0,A,I4,A,I2,A,I2,A)')         &
      '#-------------------------------------------------------------'   &
      ,'---------------------------------)'           &
      ,'#',(timeav),'--- AVERAGING TIMESTEP --- '         &
      ,nhrs,':',nminut,':',nsecs             &
      ,'   HRS:MIN:SEC AFTER INITIALIZATION '
    write (ifoutput,'(2A/A/A/2A/2A/2A)')             &
      '#------------------------------------------------------------'     &
      ,'------------'                 &
      ,'#               --------   PRECIPITATION ------    '       &
      ,'#                                                           '     &
      ,'# LEV HEIGHT   RHO(k)  PRES  |CLOUDCOVER  ECHORAINRATE  PRECCOUNT '   &
      ,'    NRRAIN      RAINCOUNT     PREC(k)     <Dvr(k)>     <qr(k)>'   &
      ,'#      (M)             (MB)  |----------  ---W/M2----   --------- '   &
      ,'    ------      ---------     -------     --------    ---------'   &
      ,'#-----------------------------------------------------------------'   &
      ,'---------------------------------------------------------------'
    write(ifoutput,'(I4,F8.2,F8.3,F7.1,8E13.5)') &
      (k          , &
      zf    (k)      , &
      rhof    (k)      , &
      presf    (k)/100.    , &
      cloudcountmn  (k)      , &
      prec_prcmn  (k)*rhof(k)*rlv  , &
      preccountmn  (k)      , &
      Nrrainmn  (k)      , &
      raincountmn  (k)      , &
      precmn    (k)*rhof(k)*rlv  , &
      Dvrmn    (k)      , &
      qrmn    (k)      , &
      k=1,kmax)
    close(ifoutput)

    open (ifoutput,file='nptend.'//cexpnr,position='append')
    write(ifoutput,'(//2A,/A,F5.0,A,I4,A,I2,A,I2,A)')         &
      '#-------------------------------------------------------------'   &
      ,'---------------------------------)'           &
      ,'#',(timeav),'--- AVERAGING TIMESTEP --- '         &
      ,nhrs,':',nminut,':',nsecs             &
      ,'   HRS:MIN:SEC AFTER INITIALIZATION '
    write (ifoutput,'(2A/A/A/2A/A/A)')             &
      '#------------------------------------------------------------'     &
      , '------------'               &
      ,'#               --------   T E N D E N C I E S NRAIN ------    '     &
      ,'#                                                           '     &
      ,'# LEV HEIGHT   PRES  |  AUTO         ACCR          SEDIM    '     &
      ,'     EVAP         TOT '             &
      ,'#      (M)   (MB)  |  ---------   (#/M3/S)      ----------'     &
      ,'#-----------------------------------------------------------'
    write(ifoutput,'(I4,F8.2,F7.1,5E13.5)') &
      (k          , &
      zf    (k)      , &
      presf    (k)/100.    , &
      Npmn    (k,iauto)    , &
      Npmn    (k,iaccr)    , &
      Npmn    (k,ised)    , &
      Npmn    (k,ievap)    , &
      sum(Npmn  (k,2:nrfields))    , &
      k=1,kmax)
    close(ifoutput)

    open (ifoutput,file='qlptend.'//cexpnr,position='append')
    write(ifoutput,'(//2A,/A,F5.0,A,I4,A,I2,A,I2,A)')         &
      '#-------------------------------------------------------------'   &
      ,'---------------------------------)'           &
      ,'#',(timeav),'--- AVERAGING TIMESTEP --- '         &
      ,nhrs,':',nminut,':',nsecs             &
      ,'   HRS:MIN:SEC AFTER INITIALIZATION '
    write (ifoutput,'(2A/A/A/2A/A/A)')             &
      '#------------------------------------------------------------'     &
      , '------------'               &
      ,'#               --------   T E N D E N C I E S QRAIN ------    '   &
      ,'#                                                           '     &
      ,'# LEV HEIGHT   PRES  |  AUTO         ACCR          SEDIM    '     &
      ,'     EVAP         TOT '             &
      ,'#      (M)   (MB)  |  ---------   (KG/KG/S)      ----------'     &
      ,'#-----------------------------------------------------------'
    write(ifoutput,'(I4,F8.2,F7.1,5E13.5)') &
      (k          , &
      zf    (k)      , &
      presf    (k)/100.    , &
      qlpmn    (k,iauto)    , &
      qlpmn    (k,iaccr)    , &
      qlpmn    (k,ised)    , &
      qlpmn    (k,ievap)    , &
      sum(qlpmn  (k,2:nrfields))    , &
                        k=1,kmax)
    close(ifoutput)

    open (ifoutput,file='qtptend.'//cexpnr,position='append')
    write(ifoutput,'(//2A,/A,F5.0,A,I4,A,I2,A,I2,A)')         &
      '#-------------------------------------------------------------'   &
      ,'---------------------------------)'           &
      ,'#',(timeav),'--- AVERAGING TIMESTEP --- '         &
      ,nhrs,':',nminut,':',nsecs             &
      ,'   HRS:MIN:SEC AFTER INITIALIZATION '
    write (ifoutput,'(2A/A/A/2A/A/A)')             &
      '#------------------------------------------------------------'     &
      , '------------'               &
      ,'#               --------   T E N D E N C I E S QTP ------    '   &
      ,'#                                                           '     &
      ,'# LEV HEIGHT   PRES  |  AUTO         ACCR          SEDIM    '     &
      ,'     EVAP         TOT '             &
      ,'#      (M)   (MB)  |  ---------   (KG/KG/S)      ----------'     &
      ,'#-----------------------------------------------------------'
    write(ifoutput,'(I4,F8.2,F7.1,5E13.5)') &
      (k          , &
      zf    (k)      , &
      presf    (k)/100.    , &
      qtpmn    (k,iauto)    , &
      qtpmn    (k,iaccr)    , &
      qtpmn    (k,ised)    , &
      qtpmn    (k,ievap)    , &
      sum    (qtpmn(k,2:nrfields))  , &
      k=1,kmax)
      close(ifoutput)
      if (lnetcdf) then
        vars(:, 1) = qrmn
        vars(:, 2) = cloudcountmn    ! = cloudfraction
        vars(:, 3) = raincountmn
        vars(:, 4) = preccountmn
        vars(:, 5) = prec_fracmn
!        vars(:, 6) = precmn    (:)*rhoz(2,2,:)*rlv
        vars(:, 6) = precmn      ! now in kg/kg m/s
!        vars(:, 7) = prec_prcmn  (:)*rhoz(2,2,:)*rlv
        vars(:, 7) = sedmn
        vars(:, 8) = prec_prcmn  ! now in kg/kg m/s
        vars(:, 9) = Dvrmn
        vars(:,10) = Nrcloudmn
        vars(:,11) = Nrrainmn

        vars(:,12) =Npmn     (:,iauto)
        vars(:,13) =Npmn     (:,iaccr)
        vars(:,14) =Npmn     (:,ised)
        vars(:,15) =Npmn     (:,ievap)
        vars(:,16) =sum(Npmn (:,2:nrfields),2)
        vars(:,17) =qrsrc    (:)/(rslabs*tint)
        vars(:,18) =0.
        vars(:,19) =qrsed    (:)/(rslabs*tint)
        vars(:,20) =qrevap   (:)/(rslabs*tint)
        vars(:,21) =(qrsrc+qrsed+qrevap)/(rslabs*tint)
        vars(:,22) =-qrsrc    (:)/(rslabs*tint)
        vars(:,23) =0.
        vars(:,24) =-qrsed    (:)/(rslabs*tint)
        vars(:,25) =-qrevap   (:)/(rslabs*tint)
        vars(:,26) =-(qrsrc+qrsed+qrevap)/(rslabs*tint)

        call writestat_nc(ncid_prof,nvar,ncname,vars(1:kmax,:),nrec_prof,kmax)
      end if

    end if

    qrsrc = 0.
    qrsed = 0.
    qrevap = 0.
    tint = 0.

    cloudcountmn = 0.0
    raincountmn  = 0.0
    preccountmn  = 0.0
    prec_prcmn   = 0.0
    prec_fracmn  = 0.0
    Dvrmn        = 0.0
    Nrrainmn     = 0.0
    Nrcloudmn    = 0.0
    precmn       = 0.0
    sedmn        = 0.0
    qrmn         = 0.0
    Npmn         = 0.0
    qlpmn        = 0.0
    qtpmn        = 0.0

  end subroutine writebulkmicrostat

!------------------------------------------------------------------------------!

  subroutine exitbulkmicrostat
    use modmicrodata, only : Dvr,precep,qrevap,qrsed,qrsrc
    implicit none
    if (.not. lmicrostat)  return

    deallocate( Dvr,precep)
    deallocate( qrevap, qrsed, qrsrc )
    deallocate(Npav      , &
         Npmn      , &
         qlpav    , &
         qlpmn    , &
         qtpav    , &
         qtpmn    )
    deallocate(precavl    , &
         precav    , &
         precmn    , &
         sedavl   , &
         sedav    , &
         sedmn    , &
         preccountavl    , &
         preccountav    , &
         preccountmn    , &
         prec_prcavl    , &
         prec_prcav    , &
         prec_prcmn    , &
         prec_fracavl    , &
         prec_fracav    , &
         prec_fracmn    , &
         cloudcountavl  , &
         cloudcountav    , &
         cloudcountmn    , &
         raincountavl    , &
         raincountav    , &
         raincountmn    , &
         Nrrainavl    , &
         Nrrainav    , &
         Nrrainmn    , &
         Nrcloudavl    , &
         Nrcloudav    , &
         Nrcloudmn    , &
         qravl    , &
         qrav      , &
         qrmn      , &
         Dvravl    , &
         Dvrav    , &
         Dvrmn    )
  end subroutine exitbulkmicrostat

!------------------------------------------------------------------------------!

end module
