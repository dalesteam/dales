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
PUBLIC  :: initbulkmicrostat, bulkmicrostat, exitbulkmicrostat, bulkmicrotend, bulkaertend
save
!NetCDF variables
  integer,parameter :: nvar = 74
  character(80),dimension(nvar,4) :: ncname
  character(80),dimension(1,4) :: tncname
  real          :: dtav, timeav
  integer(kind=longint):: idtav, itimeav, tnext, tnextwrite
  integer          :: nsamples
  logical          :: lmicrostat = .false.
  integer, parameter      :: &
               nrfields = 5, &
               iauto = 2, &
               iaccr = 3, &
               ievpr = 4, &
               isedr = 5, &

               aerfields = 10, & 
               jacti = 1, &
               jscvc = 2, &
               jevpc = 3, &
               jscvr = 4, &
               jevpr = 5, &
               jauto = 6, &
               jaccr = 7, &
               jslfc = 8, &
               jslfr = 9, &
               jsedr = 10

  real, allocatable, dimension(:,:)  :: &
                 Npav,   Npmn, &
                qlpav,  qlpmn, &
                qtpav,  qtpmn, &
                 ncav,   ncmn, &
               so4cav, so4cmn, &
                bccav,  bccmn, &
               pomcav, pomcmn, & 
                sscav,  sscmn, &
                ducav,  ducmn, &
                 nrav,   nrmn, &
               so4rav, so4rmn, &
                bcrav,  bcrmn, &
               pomrav, pomrmn, & 
                ssrav,  ssrmn, &
                durav,  durmn
 
  real, allocatable, dimension(:)    :: precavl  , &
               precav  , &
               precmn  , &
               preccountavl  , &
               preccountav  , &
               preccountmn  , &
               prec_prcavl  , &
               prec_prcav  , &
               prec_prcmn  , &
               cloudcountavl, &
               cloudcountav  , &
               cloudcountmn  , &
               raincountavl  , &
               raincountav  , &
               raincountmn  , &
               Nrrainavl  , &
               Nrrainav  , &
               Nrrainmn  , &
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
              dtav_glob, timeav_glob, ladaptive, k1, dtmax,btime,tres
    use modstat_nc, only : lnetcdf,define_nc,ncinfo,writestat_dims_nc
    use modgenstat, only : idtav_prof=>idtav, itimeav_prof=>itimeav,ncid_prof=>ncid
    use modmicrodata,only: imicro, imicro_bulk, imicro_sice
    implicit none
    integer      :: ierr

    namelist/NAMBULKMICROSTAT/ &
    lmicrostat, dtav, timeav

    if ((imicro /=imicro_bulk) .and. (imicro /= imicro_sice)) return

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
    nsamples   = itimeav/idtav

    if (.not. lmicrostat) return
    if (abs(timeav/dtav - nsamples) > 1e-4) then
      stop 'timeav must be an integer multiple of dtav (NAMBULKMICROSTAT)'
    end if
    if (.not. ladaptive .and. abs(dtav/dtmax - nint(dtav/dtmax)) > 1e-4) then
      stop 'dtav must be an integer multiple of dtmax (NAMBULKMICROSTAT)'
    end if

    allocate( &
         Npav (k1, nrfields),   Npmn (k1, nrfields), &
        qlpav (k1, nrfields),  qlpmn (k1, nrfields), &
        qtpav (k1, nrfields),  qtpmn (k1, nrfields), &  

         ncav (k1, aerfields),   ncmn (k1, aerfields), &
       so4cav (k1, aerfields), so4cmn (k1, aerfields), &
        bccav (k1, aerfields),  bccmn (k1, aerfields), &
       pomcav (k1, aerfields), pomcmn (k1, aerfields), &
        sscav (k1, aerfields),  sscmn (k1, aerfields), &
        ducav (k1, aerfields),  ducmn (k1, aerfields), &
         nrav (k1, aerfields),   nrmn (k1, aerfields), &
       so4rav (k1, aerfields), so4rmn (k1, aerfields), &
        bcrav (k1, aerfields),  bcrmn (k1, aerfields), &
       pomrav (k1, aerfields), pomrmn (k1, aerfields), &
        ssrav (k1, aerfields),  ssrmn (k1, aerfields), &
        durav (k1, aerfields),  durmn (k1, aerfields)  )

    allocate(precavl (k1), &
       precav        (k1), &
       precmn        (k1), &
       preccountavl  (k1), &
       preccountav   (k1), &
       preccountmn   (k1), &
       prec_prcavl   (k1), &
       prec_prcav    (k1), &
       prec_prcmn    (k1), &
       cloudcountavl (k1), &
       cloudcountav  (k1), &
       cloudcountmn  (k1), &
       raincountavl  (k1), &
       raincountav   (k1), &
       raincountmn   (k1), &
       Nrrainavl     (k1), &
       Nrrainav      (k1), &
       Nrrainmn      (k1), &
       qravl         (k1), &
       qrav          (k1), &
       qrmn          (k1), &
       Dvravl        (k1), &
       Dvrav         (k1), &
       Dvrmn         (k1))

    Npmn   = 0.0
    qlpmn  = 0.0
    qtpmn  = 0.0

    ncmn   = 0.0 
    so4cmn = 0.0
    bccmn  = 0.0
    pomcmn = 0.0
    sscmn  = 0.0
    ducmn  = 0.0
    nrmn   = 0.0 
    so4rmn = 0.0
    bcrmn  = 0.0
    pomrmn = 0.0
    ssrmn  = 0.0
    durmn  = 0.0

    precmn       = 0.0
    preccountmn  = 0.0
    prec_prcmn   = 0.0
    cloudcountmn = 0.0
    raincountmn  = 0.0
    Nrrainmn     = 0.0
    qrmn         = 0.0
    Dvrmn        = 0.0

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
      idtav      = idtav_prof
      itimeav    = itimeav_prof
      tnext      = idtav+btime
      tnextwrite = itimeav+btime
      nsamples   = itimeav/idtav

      if (myid==0) then
        call ncinfo(tncname(1,:),'time','Time',              's',  'time')
        call ncinfo(ncname( 1,:),'cfrac','Cloud fraction',   '-'    ,'tt')
        call ncinfo(ncname( 2,:),'rainrate','Echo Rain Rate','W/m^2','tt')
        call ncinfo(ncname( 3,:),'preccount','Preccount',    'W/m^2','tt')
        call ncinfo(ncname( 4,:),'nrrain','nrrain',          'W/m^2','tt')
        call ncinfo(ncname( 5,:),'raincount','raincount',    'W/m^2','tt')
        call ncinfo(ncname( 6,:),'precmn','precmn',          'W/m^2','tt')
        call ncinfo(ncname( 7,:),'dvrmn','dvrmn',            'W/m^2','tt')
        call ncinfo(ncname( 8,:),'qrmn','qrmn',              'W/m^2','tt')
        call ncinfo(ncname( 9,:),'npauto','Autoconversion rain drop tendency',           '#/m3/s', 'tt')
        call ncinfo(ncname(10,:),'npaccr','Accretion rain drop tendency',                '#/m3/s', 'tt')
        call ncinfo(ncname(11,:),'npsed', 'Sedimentation rain drop tendency',            '#/m3/s', 'tt')
        call ncinfo(ncname(12,:),'npevap','Evaporation rain drop tendency',              '#/m3/s', 'tt')
        call ncinfo(ncname(13,:),'nptot', 'Total rain drop tendency',                    '#/m3/s', 'tt')
        call ncinfo(ncname(14,:),'qrpauto','Autoconversion rain water content tendency', 'kg/kg/s','tt')
        call ncinfo(ncname(15,:),'qrpaccr','Accretion rain water content tendency',      'kg/kg/s','tt')
        call ncinfo(ncname(16,:),'qrpsed', 'Sedimentation rain water content tendency',  'kg/kg/s','tt')
        call ncinfo(ncname(17,:),'qrpevap','Evaporation rain water content tendency',    'kg/kg/s','tt')
        call ncinfo(ncname(18,:),'qrptot', 'Total rain water content tendency',          'kg/kg/s','tt')
        call ncinfo(ncname(19,:),'qtpauto','Autoconversion total water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(20,:),'qtpaccr','Accretion total water content tendency',     'kg/kg/s','tt')
        call ncinfo(ncname(21,:),'qtpsed', 'Sedimentation total water content tendency', 'kg/kg/s','tt')
        call ncinfo(ncname(22,:),'qtpevap','Evaporation total water content tendency',   'kg/kg/s','tt')
        call ncinfo(ncname(23,:),'qtptot', 'Total total water content tendency',         'kg/kg/s','tt')

        ! New aerosol statistics
        ! TODO: Missing, decrease in free aerosol number resulting from scavenging
        !       Reason is that we now use cloud/rain variables, but number scavenging
        !       only applies to free aerosol tracers
        call ncinfo(ncname(24,:),  'ncacti','Activation           cloud drop tendency','?/?/s','tt')        
        call ncinfo(ncname(25,:),  'ncslfc','Cloud selfcollection cloud drop tendency','?/?/s','tt')        
        call ncinfo(ncname(26,:),  'ncauto','Autoconversion       cloud drop tendency','?/?/s','tt')        
        call ncinfo(ncname(27,:),  'ncaccr','Accretion            cloud drop tendency','?/?/s','tt')        
        call ncinfo(ncname(28,:),  'ncevpc','Cloud evaporation    cloud drop tendency','?/?/s','tt')        

        call ncinfo(ncname(29,:),'so4cacti','Activation           in-cloud sulphate tendency','?/?/s','tt')
        call ncinfo(ncname(30,:),'so4cauto','Autoconversion       in-cloud sulphate tendency','?/?/s','tt')
        call ncinfo(ncname(31,:),'so4caccr','Accretion            in-cloud sulphate tendency','?/?/s','tt')
        call ncinfo(ncname(32,:),'so4cscvc','In-cloud scavenging  in-cloud sulphate tendency','?/?/s','tt')
        call ncinfo(ncname(33,:),'so4cevpc','Cloud evaporation    in-cloud sulphate tendency','?/?/s','tt')

        call ncinfo(ncname(34,:), 'bccacti','Activation           in-cloud black carbon tendency','?/?/s','tt')
        call ncinfo(ncname(35,:), 'bccauto','Autoconversion       in-cloud black carbon tendency','?/?/s','tt')
        call ncinfo(ncname(36,:), 'bccaccr','Accretion            in-cloud black carbon tendency','?/?/s','tt')
        call ncinfo(ncname(37,:), 'bccscvc','In-cloud scavenging  in-cloud black carbon tendency','?/?/s','tt')
        call ncinfo(ncname(38,:), 'bccevpc','Cloud evaporation    in-cloud black carbon tendency','?/?/s','tt')

        call ncinfo(ncname(39,:),'pomcacti','Activation           in-cloud organic matter tendency','?/?/s','tt')
        call ncinfo(ncname(40,:),'pomcauto','Autoconversion       in-cloud organic matter tendency','?/?/s','tt')
        call ncinfo(ncname(41,:),'pomcaccr','Accretion            in-cloud organic matter tendency','?/?/s','tt')
        call ncinfo(ncname(42,:),'pomcscvc','In-cloud scavenging  in-cloud organic matter tendency','?/?/s','tt')
        call ncinfo(ncname(43,:),'pomcevpc','Cloud evaporation    in-cloud organic matter tendency','?/?/s','tt')

        call ncinfo(ncname(44,:), 'sscacti','Activation           in-cloud sea salt tendency','?/?/s','tt')
        call ncinfo(ncname(45,:), 'sscauto','Autoconversion       in-cloud sea salt tendency','?/?/s','tt')
        call ncinfo(ncname(46,:), 'sscaccr','Accretion            in-cloud sea salt tendency','?/?/s','tt')
        call ncinfo(ncname(47,:), 'sscscvc','In-cloud scavenging  in-cloud sea salt tendency','?/?/s','tt')
        call ncinfo(ncname(48,:), 'sscevpc','Cloud evaporation    in-cloud sea salt tendency','?/?/s','tt')

        call ncinfo(ncname(49,:), 'ducacti','Activation           in-cloud dust tendency','?/?/s','tt')
        call ncinfo(ncname(50,:), 'ducauto','Autoconversion       in-cloud dust tendency','?/?/s','tt')
        call ncinfo(ncname(51,:), 'ducaccr','Accretion            in-cloud dust tendency','?/?/s','tt')
        call ncinfo(ncname(52,:), 'ducscvc','In-cloud scavenging  in-cloud dust tendency','?/?/s','tt')
        call ncinfo(ncname(53,:), 'ducevpc','Cloud evaporation    in-cloud dust tendency','?/?/s','tt')

        call ncinfo(ncname(54,:), 'nrauto','Autoconversion          rain drop tendency','?/?/s','tt')
        call ncinfo(ncname(55,:), 'nrevpr','Rain evaporation        rain drop tendency','?/?/s','tt')
        call ncinfo(ncname(56,:), 'nrsedr','Rain sedimentation      rain drop tendency','?/?/s','tt')
        call ncinfo(ncname(57,:), 'nrslfr','Rain selfcollection     rain drop tendency','?/?/s','tt')

        call ncinfo(ncname(58,:),'so4rscvr','Below-cloud scavenging in-rain sulphate tendency','?/?/s','tt')
        call ncinfo(ncname(59,:),'so4revpr','Rain evaporation       in-rain sulphate tendency','?/?/s','tt')
        call ncinfo(ncname(60,:),'so4rsedr','Rain sedimentation     in-rain sulphate tendency','?/?/s','tt')

        call ncinfo(ncname(61,:), 'bcrscvr','Below-cloud scavenging in-rain black carbon tendency','?/?/s','tt')
        call ncinfo(ncname(62,:), 'bcrevpr','Rain evaporation       in-rain black carbon tendency','?/?/s','tt')
        call ncinfo(ncname(63,:), 'bcrsedr','Rain sedimentation     in-rain black carbon tendency','?/?/s','tt')

        call ncinfo(ncname(64,:),'pomrscvr','Below-cloud scavenging in-rain organic matter tendency','?/?/s','tt')
        call ncinfo(ncname(65,:),'pomrevpr','Rain evaporation       in-rain organic matter tendency','?/?/s','tt')
        call ncinfo(ncname(66,:),'pomrsedr','Rain sedimentation     in-rain organic matter tendency','?/?/s','tt')

        call ncinfo(ncname(67,:), 'ssrscvr','Below-cloud scavenging in-rain sea salt tendency','?/?/s','tt')
        call ncinfo(ncname(68,:), 'ssrevpr','Rain evaporation       in-rain sea salt tendency','?/?/s','tt')
        call ncinfo(ncname(69,:), 'ssrsedr','Rain sedimentation     in-rain sea salt tendency','?/?/s','tt')

        call ncinfo(ncname(70,:), 'durscvr','Below-cloud scavenging in-rain dust tendency','?/?/s','tt')
        call ncinfo(ncname(71,:), 'durevpr','Rain evaporation       in-rain dust tendency','?/?/s','tt')
        call ncinfo(ncname(72,:), 'dursedr','Rain sedimentation     in-rain dust tendency','?/?/s','tt')
        
        call ncinfo(ncname(73,:),  'ncscvc','In-cloud number scavenging tendency','#/kg/s','tt')        
        call ncinfo(ncname(74,:),  'nrscvr','In-rain  number scavenging tendency','#/kg/s','tt')        

        call define_nc( ncid_prof, nvar, ncname)
      end if
   end if

  end subroutine initbulkmicrostat

!------------------------------------------------------------------------------!
!> General routine, does the timekeeping
  subroutine bulkmicrostat
    use modglobal,    only  : rk3step, timee, dt_lim
    implicit none
    if (.not. lmicrostat)  return
    if (rk3step /= 3)  return
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
    use modmpi,       only : my_real, mpi_sum, comm3d, mpierr
    use modglobal,    only : i1, j1, k1, ijtot
    use modmicrodata, only : qr, precep, Dvr, Nr, epscloud, epsqr, epsprec,imicro, imicro_bulk
    use modfields,    only : ql0

    implicit none

    integer :: k

    precav       = 0.0
    preccountav  = 0.0
    prec_prcav   = 0.0
    cloudcountav = 0.0
    raincountav  = 0.0
    Nrrainav     = 0.0
    qrav         = 0.0
    Dvrav        = 0.0

    do k = 1,k1
      cloudcountavl(k)  = count(ql0   (2:i1,2:j1,k) > epscloud)
      raincountavl (k)  = count(qr    (2:i1,2:j1,k) > epsqr)
      preccountavl (k)  = count(precep(2:i1,2:j1,k) > epsprec)
      prec_prcavl  (k)  = sum  (precep(2:i1,2:j1,k), precep(2:i1,2:j1,k) > epsprec)
      Nrrainavl    (k)  = sum  (Nr    (2:i1,2:j1,k))
      precavl      (k)  = sum  (precep(2:i1,2:j1,k))
      qravl        (k)  = sum  (qr    (2:i1,2:j1,k))
      if (imicro==imicro_bulk) then
        Dvravl     (k)  = sum  (Dvr   (2:i1,2:j1,k), qr    (2:i1,2:j1,k) > epsqr)
      end if
    end do

    call MPI_ALLREDUCE(cloudcountavl, cloudcountav, k1, MY_REAL, MPI_SUM, comm3d, mpierr)
    call MPI_ALLREDUCE(raincountavl , raincountav,  k1, MY_REAL, MPI_SUM, comm3d, mpierr)
    call MPI_ALLREDUCE(preccountavl , preccountav,  k1, MY_REAL, MPI_SUM, comm3d, mpierr)
    call MPI_ALLREDUCE(prec_prcavl  , prec_prcav,   k1, MY_REAL, MPI_SUM, comm3d, mpierr)
    call MPI_ALLREDUCE(Dvravl       , Dvrav,        k1, MY_REAL, MPI_SUM, comm3d, mpierr)
    call MPI_ALLREDUCE(Nrrainavl    , Nrrainav,     k1, MY_REAL, MPI_SUM, comm3d, mpierr)
    call MPI_ALLREDUCE(precavl      , precav,       k1, MY_REAL, MPI_SUM, comm3d, mpierr)
    call MPI_ALLREDUCE(qravl        , qrav,         k1, MY_REAL, MPI_SUM, comm3d, mpierr)

    cloudcountmn = cloudcountmn + cloudcountav/ijtot
    raincountmn  = raincountmn  + raincountav /ijtot
    preccountmn  = preccountmn  + preccountav /ijtot
    prec_prcmn   = prec_prcmn   + prec_prcav  /ijtot
    Dvrmn        = Dvrmn        + Dvrav       /ijtot
    Nrrainmn     = Nrrainmn     + Nrrainav    /ijtot
    precmn       = precmn       + precav      /ijtot
    qrmn         = qrmn         + qrav        /ijtot

  end subroutine dobulkmicrostat

!------------------------------------------------------------------------------!
!> Performs the calculations for the tendencies etc.
  subroutine bulkmicrotend
    use modmpi,       only : slabsum
    use modglobal,    only : rk3step, timee, dt_lim, k1, ih, i1, jh, j1, ijtot
    use modfields,    only : qtp
    use modmicrodata, only : qrp, aer_tend, inr
   
    implicit none

    real, dimension(:), allocatable  :: avfield
    integer                          :: ifield = 0

    if (.not. lmicrostat)  return
    if (rk3step /= 3)  return
    if (timee == 0)    return
    if (timee < tnext .and. timee < tnextwrite) then
      dt_lim  = minval((/dt_lim, tnext - timee, tnextwrite - timee/))
      return
    end if

    allocate(avfield(k1))

    ifield    = mod(ifield, nrfields) + 1

    avfield    = 0.0
    call slabsum(avfield,1,k1,aer_tend(:,:,:,inr),2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    Npav(:,ifield)  = avfield - sum(Npav  (:,1:ifield-1),2)

    avfield    = 0.0
    call slabsum(avfield  ,1,k1,qrp  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    qlpav(:,ifield) = avfield - sum(qlpav  (:,1:ifield-1),2)

    avfield    = 0.0
    call slabsum(avfield  ,1,k1,qtp  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    qtpav(:,ifield) = avfield - sum(qtpav  (:,1:ifield-1),2)

    if (ifield == nrfields) then
        Npmn =   Npmn +   Npav/nsamples/ijtot;   Npav = 0.0
       qlpmn =  qlpmn +  qlpav/nsamples/ijtot;  qlpav = 0.0
       qtpmn =  qtpmn +  qtpav/nsamples/ijtot;  qtpav = 0.0
    end if

    deallocate(avfield)
  end subroutine bulkmicrotend

  ! --------------------------------------------------------------------------
  ! Performs the calculations for the aerosol tendencies statistics
  subroutine bulkaertend
     use modmpi,       only : slabsum
     use modglobal,    only : rk3step, timee, dt_lim, k1, ih, i1, jh, j1, ijtot
     use modmicrodata, only : naer, &
                              inus_n, iais_n, iacs_n, icos_n, iaii_n, iaci_n, icoi_n, &  
                              inc, iso4cld, ibccld, ipomcld, isscld, iducld, &
                              inr, iso4rai, ibcrai, ipomrai, issrai, idurai, & 
                              aer_acti, aer_scvc, aer_evpc, &
                              aer_scvr, aer_evpr, & 
                              aer_auto, aer_accr, &
                              aer_slfc, aer_slfr, aer_sedr 
     implicit none

     real, dimension(:),       allocatable  :: avfield
     real, dimension(:,:,:,:), allocatable  :: tend

     integer :: ifield

     if (.not. lmicrostat)  return
     if (rk3step /= 3)  return
     if (timee == 0)    return
     if (timee < tnext .and. timee < tnextwrite) then
       dt_lim  = minval((/dt_lim, tnext - timee, tnextwrite - timee/))
       return
     end if

     allocate(avfield(k1))
     allocate(tend(2-ih:i1+ih,2-jh:j1+jh,k1,naer))

     do ifield = 1,aerfields
        select case (ifield)
           case(jacti); tend = aer_acti(:,:,:,:)
           case(jscvc); tend = aer_scvc(:,:,:,:)
           case(jevpc); tend = aer_evpc(:,:,:,:)
           case(jscvr); tend = aer_scvr(:,:,:,:)
           case(jevpr); tend = aer_evpr(:,:,:,:)
           case(jauto); tend = aer_auto(:,:,:,:)
           case(jaccr); tend = aer_accr(:,:,:,:)
           case(jslfc); tend = aer_slfc(:,:,:,:)
           case(jslfr); tend = aer_slfr(:,:,:,:)
           case(jsedr); tend = aer_sedr(:,:,:,:)
        end select
                
        avfield = 0.0
        ! In-cloud --------------------------------- ------------------------------------------!

        if (ifield == jscvc) then
           call slabsum(avfield,1,k1,tend(:,:,:,inus_n)+tend(:,:,:,iais_n)+tend(:,:,:,iacs_n)+tend(:,:,:,icos_n)+ &
                                     tend(:,:,:,iaii_n)+tend(:,:,:,iaci_n)+tend(:,:,:,icoi_n),2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
           ncmn(:,ifield)   = ncmn(:,ifield)   + avfield/nsamples/ijtot
        else
           call slabsum(avfield,1,k1,tend(:,:,:,inc    ),2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
           ncmn(:,ifield)   = ncmn(:,ifield)   + avfield/nsamples/ijtot
        endif

        avfield = 0.0
        call slabsum(avfield,1,k1,tend(:,:,:,iso4cld),2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
        so4cmn(:,ifield) = so4cmn(:,ifield) + avfield/nsamples/ijtot
    
        avfield = 0.0
        call slabsum(avfield,1,k1,tend(:,:,:,ibccld ),2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
        bccmn(:,ifield)  = bccmn(:,ifield)  + avfield/nsamples/ijtot
    
        avfield = 0.0
        call slabsum(avfield,1,k1,tend(:,:,:,ipomcld),2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
        pomcmn(:,ifield) = pomcmn(:,ifield) + avfield/nsamples/ijtot  
    
        avfield = 0.0
        call slabsum(avfield,1,k1,tend(:,:,:,isscld ),2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
        sscmn(:,ifield)  = sscmn(:,ifield)  + avfield/nsamples/ijtot        

        avfield = 0.0
        call slabsum(avfield,1,k1,tend(:,:,:,iducld ),2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
        ducmn(:,ifield)  = ducmn(:,ifield)  + avfield/nsamples/ijtot
    
        ! In-rain ----------------------------------------------------------------------------!
        avfield = 0.0
        
        if (ifield == jscvr) then
           call slabsum(avfield,1,k1,tend(:,:,:,inus_n)+tend(:,:,:,iais_n)+tend(:,:,:,iacs_n)+tend(:,:,:,icos_n)+ & 
                                     tend(:,:,:,iaii_n)+tend(:,:,:,iaci_n)+tend(:,:,:,icoi_n),2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
           nrmn(:,ifield)   = nrmn(:,ifield)   + avfield/nsamples/ijtot
        else
           call slabsum(avfield,1,k1,tend(:,:,:,inr   ),2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
           nrmn(:,ifield)   = nrmn(:,ifield)   + avfield/nsamples/ijtot
        end if

        avfield = 0.0
        call slabsum(avfield,1,k1,tend(:,:,:,iso4rai),2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
        so4rmn(:,ifield) = so4rmn(:,ifield) + avfield/nsamples/ijtot
    
        avfield = 0.0
        call slabsum(avfield,1,k1,tend(:,:,:,ibcrai ),2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
        bcrmn(:,ifield)  = bcrmn(:,ifield)  + avfield/nsamples/ijtot
    
        avfield = 0.0
        call slabsum(avfield,1,k1,tend(:,:,:,ipomrai),2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
        pomrmn(:,ifield) = pomrmn(:,ifield) + avfield/nsamples/ijtot  
        
        avfield = 0.0
        call slabsum(avfield,1,k1,tend(:,:,:,issrai ),2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
        ssrmn(:,ifield)  = ssrmn(:,ifield)  + avfield/nsamples/ijtot        

        avfield = 0.0
        call slabsum(avfield,1,k1,tend(:,:,:,idurai ),2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
        durmn(:,ifield)  = durmn(:,ifield)  + avfield/nsamples/ijtot
    
     enddo
       
     deallocate(avfield) 
     deallocate(tend)
   end subroutine bulkaertend

!------------------------------------------------------------------------------!
!> Write the stats to file
  subroutine writebulkmicrostat
    use modmpi,    only  : myid
    use modglobal,    only  : rtimee, ifoutput, cexpnr, k1,kmax, &
              rlv, zf
    use modfields,    only  : presf,rhof
    use modstat_nc, only: lnetcdf, writestat_nc
    use modgenstat, only: ncid_prof=>ncid,nrec_prof=>nrec

    implicit none
    real,dimension(k1,nvar) :: vars

    integer    :: nsecs, nhrs, nminut
    integer    :: k

    nsecs  = nint(rtimee)
    nhrs   = int (nsecs/3600)
    nminut = int (nsecs/60)-nhrs*60
    nsecs  = mod (nsecs,60)

    cloudcountmn  = cloudcountmn/nsamples
    raincountmn   = raincountmn /nsamples
    preccountmn   = preccountmn /nsamples
    prec_prcmn    = prec_prcmn  /nsamples
    Dvrmn         = Dvrmn       /nsamples
    Nrrainmn      = Nrrainmn    /nsamples
    precmn        = precmn      /nsamples
    qrmn          = qrmn        /nsamples

    where (raincountmn > 0.)
      Dvrmn        = Dvrmn / raincountmn
    elsewhere
      Dvrmn = 0.0
    end where
    where (preccountmn > 0.)
      prec_prcmn = prec_prcmn/preccountmn
    elsewhere
      prec_prcmn = 0.0
    end where

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
    write(ifoutput,'(I4,F10.2,F8.3,F7.1,8E13.5)') &
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
    write(ifoutput,'(I4,F10.2,F7.1,5E13.5)') &
      (k          , &
      zf    (k)      , &
      presf    (k)/100.    , &
      Npmn    (k,iauto)    , &
      Npmn    (k,iaccr)    , &
      Npmn    (k,isedr)    , &
      Npmn    (k,ievpr)    , &
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
    write(ifoutput,'(I4,F10.2,F7.1,5E13.5)') &
      (k          , &
      zf    (k)      , &
      presf    (k)/100.    , &
      qlpmn    (k,iauto)    , &
      qlpmn    (k,iaccr)    , &
      qlpmn    (k,isedr)    , &
      qlpmn    (k,ievpr)    , &
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
    write(ifoutput,'(I4,F10.2,F7.1,5E13.5)') &
      (k          , &
      zf    (k)      , &
      presf    (k)/100.    , &
      qtpmn    (k,iauto)    , &
      qtpmn    (k,iaccr)    , &
      qtpmn    (k,isedr)    , &
      qtpmn    (k,ievpr)    , &
      sum    (qtpmn(k,2:nrfields))  , &
      k=1,kmax)
      close(ifoutput)

      if (lnetcdf) then
        vars(:, 1) = cloudcountmn
        vars(:, 2) = prec_prcmn  (:)*rhof(:)*rlv
        vars(:, 3) = preccountmn  (:)
        vars(:, 4) = Nrrainmn  (:)
        vars(:, 5) = raincountmn  (:)
        vars(:, 6) = precmn    (:)*rhof(:)*rlv
        vars(:, 7) = Dvrmn    (:)
        vars(:, 8) = qrmn    (:)
        vars(:, 9) = Npmn (:,iauto)
        vars(:,10) = Npmn (:,iaccr)
        vars(:,11) = Npmn (:,isedr)
        vars(:,12) = Npmn (:,ievpr)
        do k=1,k1
        vars(k,13) =sum(Npmn  (k,2:nrfields))
        enddo
        vars(:,14) = qlpmn (:,iauto)
        vars(:,15) = qlpmn (:,iaccr)
        vars(:,16) = qlpmn (:,isedr)
        vars(:,17) = qlpmn (:,ievpr)
        do k=1,k1
        vars(k,18) =sum(qlpmn  (k,2:nrfields))
        enddo
        vars(:,19) = qtpmn (:,iauto)
        vars(:,20) = qtpmn (:,iaccr)
        vars(:,21) = qtpmn (:,isedr)
        vars(:,22) = qtpmn (:,ievpr)
        do k=1,k1
        vars(k,23) =sum(qtpmn  (k,2:nrfields))
        enddo
        !New aerosol statistics
        vars(:,24) =   ncmn(:,jacti)
        vars(:,25) =   ncmn(:,jslfc)
        vars(:,26) =   ncmn(:,jauto)
        vars(:,27) =   ncmn(:,jaccr)
        vars(:,28) =   ncmn(:,jevpc)

        vars(:,29) = so4cmn(:,jacti)
        vars(:,30) = so4cmn(:,jauto)
        vars(:,31) = so4cmn(:,jaccr)
        vars(:,32) = so4cmn(:,jscvc)
        vars(:,33) = so4cmn(:,jevpc)

        vars(:,34) =  bccmn(:,jacti)
        vars(:,35) =  bccmn(:,jauto)
        vars(:,36) =  bccmn(:,jaccr)
        vars(:,37) =  bccmn(:,jscvc)
        vars(:,38) =  bccmn(:,jevpc)

        vars(:,39) = pomcmn(:,jacti)
        vars(:,40) = pomcmn(:,jauto)
        vars(:,41) = pomcmn(:,jaccr)
        vars(:,42) = pomcmn(:,jscvc)
        vars(:,43) = pomcmn(:,jevpc)

        vars(:,44) =  sscmn(:,jacti)
        vars(:,45) =  sscmn(:,jauto)
        vars(:,46) =  sscmn(:,jaccr)
        vars(:,47) =  sscmn(:,jscvc)
        vars(:,48) =  sscmn(:,jevpc)

        vars(:,49) =  ducmn(:,jacti)
        vars(:,50) =  ducmn(:,jauto)
        vars(:,51) =  ducmn(:,jaccr)
        vars(:,52) =  ducmn(:,jscvc)
        vars(:,53) =  ducmn(:,jevpc)

        vars(:,54) =   nrmn(:,jauto)      
        vars(:,55) =   nrmn(:,jevpr)
        vars(:,56) =   nrmn(:,jsedr)
        vars(:,57) =   nrmn(:,jslfr)

        vars(:,58) = so4rmn(:,jscvr)
        vars(:,59) = so4rmn(:,jevpr)
        vars(:,60) = so4rmn(:,jsedr)
                
        vars(:,61) =  bcrmn(:,jscvr)
        vars(:,62) =  bcrmn(:,jevpr)
        vars(:,63) =  bcrmn(:,jsedr)
        
        vars(:,64) = pomrmn(:,jscvr)
        vars(:,65) = pomrmn(:,jevpr)
        vars(:,66) = pomrmn(:,jsedr)

        vars(:,67) =  ssrmn(:,jscvr)
        vars(:,68) =  ssrmn(:,jevpr)
        vars(:,69) =  ssrmn(:,jsedr)
        
        vars(:,70) =  durmn(:,jscvr)
        vars(:,71) =  durmn(:,jevpr)
        vars(:,72) =  durmn(:,jsedr)

        vars(:,73) =   ncmn(:,jscvc)
        vars(:,74) =   nrmn(:,jscvr)

        call writestat_nc(ncid_prof,nvar,ncname,vars(1:kmax,:),nrec_prof,kmax)
      end if

    end if

    cloudcountmn = 0.0
    raincountmn  = 0.0
    preccountmn  = 0.0
    prec_prcmn   = 0.0
    Dvrmn        = 0.0
    Nrrainmn     = 0.0
    precmn       = 0.0
    qrmn         = 0.0

    Npmn   = 0.0
    qlpmn  = 0.0
    qtpmn  = 0.0

    ncmn   = 0.0
    so4cmn = 0.0
    bccmn  = 0.0
    pomcmn = 0.0
    sscmn  = 0.0
    ducmn  = 0.0
    nrmn   = 0.0
    so4rmn = 0.0
    bcrmn  = 0.0
    pomrmn = 0.0
    ssrmn  = 0.0
    durmn  = 0.0

  end subroutine writebulkmicrostat

!------------------------------------------------------------------------------!

  subroutine exitbulkmicrostat
    implicit none
    if (.not. lmicrostat)  return

    deallocate(  Npav,   Npmn, &
                qlpav,  qlpmn, &
                qtpav,  qtpmn, &
                 ncav,   ncmn, &
               so4cav, so4cmn, &
                bccav,  bccmn, &
               pomcav, pomcmn, &
                sscav,  sscmn, &
                ducav,  ducmn, &
                 nrav,   nrmn, &
               so4rav, so4rmn, &
                bcrav,  bcrmn, &
               pomrav, pomrmn, &
                ssrav,  ssrmn, &
                durav,  durmn  )

    deallocate(precavl, &
         precav       , &
         precmn       , &
         preccountavl , &
         preccountav  , &
         preccountmn  , &
         prec_prcavl  , &
         prec_prcav   , &
         prec_prcmn   , &
         cloudcountavl, &
         cloudcountav , &
         cloudcountmn , &
         raincountavl , &
         raincountav  , &
         raincountmn  , &
         Nrrainavl    , &
         Nrrainav     , &
         Nrrainmn     , &
         qravl        , &
         qrav         , &
         qrmn         , &
         Dvravl       , &
         Dvrav        , &
         Dvrmn)

  end subroutine exitbulkmicrostat

!------------------------------------------------------------------------------!

end module
