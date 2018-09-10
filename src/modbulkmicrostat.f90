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
PUBLIC  :: initbulkmicrostat, bulkmicrostat, exitbulkmicrostat, bulkmicrotend
save
!NetCDF variables
  integer,parameter :: nvar = 74
  character(80),dimension(nvar,4) :: ncname
  character(80),dimension(1,4) :: tncname
  real          :: dtav, timeav
  integer(kind=longint):: idtav, itimeav, tnext, tnextwrite
  integer          :: nsamples
  logical          :: lmicrostat = .false.
  integer, parameter      :: nrfields = 10 , &
               iacti = 2 , &  
               islfc = 3 , & 
               isedc = 4 , &
               iauto = 5 , &
               iaccr = 6 , &
               ievpr = 7 , &
               isedr = 8 , &
               iscav = 9 , &
               ievpc = 10      
!               iauto  = 2 , &
!               iaccr  = 3 , &
!               ievap  = 4 , &
!               ised   = 5
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
    nsamples = itimeav/idtav

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
         ncav (k1, nrfields),   ncmn (k1, nrfields), &  
       so4cav (k1, nrfields), so4cmn (k1, nrfields), &  
        bccav (k1, nrfields),  bccmn (k1, nrfields), &  
       pomcav (k1, nrfields), pomcmn (k1, nrfields), &  
        sscav (k1, nrfields),  sscmn (k1, nrfields), &  
        ducav (k1, nrfields),  ducmn (k1, nrfields), &  
       so4rav (k1, nrfields), so4rmn (k1, nrfields), &  
        bcrav (k1, nrfields),  bcrmn (k1, nrfields), &  
       pomrav (k1, nrfields), pomrmn (k1, nrfields), &  
        ssrav (k1, nrfields),  ssrmn (k1, nrfields), &  
        durav (k1, nrfields),  durmn (k1, nrfields)  )  

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
        call ncinfo(ncname(26,:),  'ncsedc','Cloud sedimentation  cloud drop tendency','?/?/s','tt')        
        call ncinfo(ncname(27,:),  'ncauto','Autoconversion       cloud drop tendency','?/?/s','tt')        
        call ncinfo(ncname(28,:),  'ncaccr','Accretion            cloud drop tendency','?/?/s','tt')        
        call ncinfo(ncname(29,:),  'ncevpc','Cloud evaporation    cloud drop tendency','?/?/s','tt')        

        call ncinfo(ncname(30,:),'so4cacti','Activation           in-cloud sulphate tendency','?/?/s','tt')
        call ncinfo(ncname(31,:),'so4csedc','Cloud sedimentation  in-cloud sulphate tendency','?/?/s','tt')
        call ncinfo(ncname(32,:),'so4cauto','Autoconversion       in-cloud sulphate tendency','?/?/s','tt')
        call ncinfo(ncname(33,:),'so4caccr','Accretion            in-cloud sulphate tendency','?/?/s','tt')
        call ncinfo(ncname(34,:),'so4cscvc','In-cloud scavenging  in-cloud sulphate tendency','?/?/s','tt')
        call ncinfo(ncname(35,:),'so4cevpc','Cloud evaporation    in-cloud sulphate tendency','?/?/s','tt')

        call ncinfo(ncname(36,:), 'bccacti','Activation           in-cloud black carbon tendency','?/?/s','tt')
        call ncinfo(ncname(37,:), 'bccsedc','Cloud sedimentation  in-cloud black carbon tendency','?/?/s','tt')
        call ncinfo(ncname(38,:), 'bccauto','Autoconversion       in-cloud black carbon tendency','?/?/s','tt')
        call ncinfo(ncname(39,:), 'bccaccr','Accretion            in-cloud black carbon tendency','?/?/s','tt')
        call ncinfo(ncname(40,:), 'bccscvc','In-cloud scavenging  in-cloud black carbon tendency','?/?/s','tt')
        call ncinfo(ncname(41,:), 'bccevpc','Cloud evaporation    in-cloud black carbon tendency','?/?/s','tt')

        call ncinfo(ncname(42,:),'pomcacti','Activation           in-cloud organic matter tendency','?/?/s','tt')
        call ncinfo(ncname(43,:),'pomcsedc','Cloud sedimentation  in-cloud organic matter tendency','?/?/s','tt')
        call ncinfo(ncname(44,:),'pomcauto','Autoconversion       in-cloud organic matter tendency','?/?/s','tt')
        call ncinfo(ncname(45,:),'pomcaccr','Accretion            in-cloud organic matter tendency','?/?/s','tt')
        call ncinfo(ncname(46,:),'pomcscvc','In-cloud scavenging  in-cloud organic matter tendency','?/?/s','tt')
        call ncinfo(ncname(47,:),'pomcevpc','Cloud evaporation    in-cloud organic matter tendency','?/?/s','tt')

        call ncinfo(ncname(48,:), 'sscacti','Activation           in-cloud sea salt tendency','?/?/s','tt')
        call ncinfo(ncname(49,:), 'sscsedc','Cloud sedimentation  in-cloud sea salt tendency','?/?/s','tt')
        call ncinfo(ncname(50,:), 'sscauto','Autoconversion       in-cloud sea salt tendency','?/?/s','tt')
        call ncinfo(ncname(51,:), 'sscaccr','Accretion            in-cloud sea salt tendency','?/?/s','tt')
        call ncinfo(ncname(52,:), 'sscscvc','In-cloud scavenging  in-cloud sea salt tendency','?/?/s','tt')
        call ncinfo(ncname(53,:), 'sscevpc','Cloud evaporation    in-cloud sea salt tendency','?/?/s','tt')

        call ncinfo(ncname(54,:), 'ducacti','Activation           in-cloud dust tendency','?/?/s','tt')
        call ncinfo(ncname(55,:), 'ducsedc','Cloud sedimentation  in-cloud dust tendency','?/?/s','tt')
        call ncinfo(ncname(56,:), 'ducauto','Autoconversion       in-cloud dust tendency','?/?/s','tt')
        call ncinfo(ncname(57,:), 'ducaccr','Accretion            in-cloud dust tendency','?/?/s','tt')
        call ncinfo(ncname(58,:), 'ducscvc','In-cloud scavenging  in-cloud dust tendency','?/?/s','tt')
        call ncinfo(ncname(59,:), 'ducevpc','Cloud evaporation    in-cloud dust tendency','?/?/s','tt')

        call ncinfo(ncname(60,:),'so4revpr','Rain evaporation       in-rain sulphate tendency','?/?/s','tt')
        call ncinfo(ncname(61,:),'so4rsedr','Rain sedimentation     in-rain sulphate tendency','?/?/s','tt')
        call ncinfo(ncname(62,:),'so4rscvr','Below-cloud scavenging in-rain sulphate tendency','?/?/s','tt')

        call ncinfo(ncname(63,:), 'bcrevpr','Rain evaporation       in-rain black carbon tendency','?/?/s','tt')
        call ncinfo(ncname(64,:), 'bcrsedr','Rain sedimentation     in-rain black carbon tendency','?/?/s','tt')
        call ncinfo(ncname(65,:), 'bcrscvr','Below-cloud scavenging in-rain black carbon tendency','?/?/s','tt')

        call ncinfo(ncname(66,:),'pomrevpr','Rain evaporation       in-rain organic matter tendency','?/?/s','tt')
        call ncinfo(ncname(67,:),'pomrsedr','Rain sedimentation     in-rain organic matter tendency','?/?/s','tt')
        call ncinfo(ncname(68,:),'pomrscvr','Below-cloud scavenging in-rain organic matter tendency','?/?/s','tt')

        call ncinfo(ncname(69,:), 'ssrevpr','Rain evaporation       in-rain sea salt tendency','?/?/s','tt')
        call ncinfo(ncname(70,:), 'ssrsedr','Rain sedimentation     in-rain sea salt tendency','?/?/s','tt')
        call ncinfo(ncname(71,:), 'ssrscvr','Below-cloud scavenging in-rain sea salt tendency','?/?/s','tt')

        call ncinfo(ncname(72,:), 'durevpr','Rain evaporation       in-rain dust tendency','?/?/s','tt')
        call ncinfo(ncname(73,:), 'dursedr','Rain sedimentation     in-rain dust tendency','?/?/s','tt')
        call ncinfo(ncname(74,:), 'durscvr','Below-cloud scavenging in-rain dust tendency','?/?/s','tt')
        
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
    use modmpi,    only  : my_real, mpi_sum, comm3d, mpierr
    use modglobal,    only  : i1, j1, k1, ijtot
    use modmicrodata,  only  : qr, precep, Dvr, Nr, epscloud, epsqr, epsprec,imicro, imicro_bulk
    use modfields,  only  : ql0
    implicit none

    integer      :: k

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
    use modmicrodata, only : qrp, aer_tend, &
                             inc, iso4cld, ibccld, ipomcld, isscld, iducld, &
                             inr, iso4rai, ibcrai, ipomrai, issrai, idurai
   
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

!    write(6,"(I1, A28)") ifield, 'Start gathering tendencies:' 
    avfield    = 0.0
!   call slabsum(avfield,1,k1,Nrp                                   ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    call slabsum(avfield,1,k1,aer_tend(:,:,:,inr),2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    Npav(:,ifield)  = avfield - sum(Npav  (:,1:ifield-1),2)
!    write(6,"(A3)") "Nr"

    avfield    = 0.0
    call slabsum(avfield  ,1,k1,qrp  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    qlpav(:,ifield) = avfield - sum(qlpav  (:,1:ifield-1),2)
!    write(6,"(A3)") "qr"

    avfield    = 0.0
    call slabsum(avfield  ,1,k1,qtp  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    qtpav(:,ifield) = avfield - sum(qtpav  (:,1:ifield-1),2)
!    write(6,"(A3)") "qt"

    ! -------------------------------------------------------------------------------------------------------------
    ! Cloud droplets (inc)
    avfield    = 0.0
    call slabsum(avfield,1,k1,aer_tend(:,:,:,inc    ),2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    ncav(:,ifield)  = avfield - sum(ncav  (:,1:ifield-1),2)
!    write(6,"(A3)") "nc"

    ! Rain droplets (inr)

    ! In-cloud sulphate (iso4cld)
    avfield    = 0.0
    call slabsum(avfield,1,k1,aer_tend(:,:,:,iso4cld),2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    so4cav(:,ifield)  = avfield - sum(so4cav  (:,1:ifield-1),2)
!    write(6,"(A5)") "so4c"

    ! In-cloud black carbon (ibccld)
    avfield    = 0.0
    call slabsum(avfield,1,k1,aer_tend(:,:,:,ibccld ),2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    bccav (:,ifield)  = avfield - sum(bccav  (:,1:ifield-1),2)
!    write(6,"(A5)") "so4c"

    ! In-cloud organic matter (ipomcld)
    avfield    = 0.0
    call slabsum(avfield,1,k1,aer_tend(:,:,:,ipomcld),2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    pomcav(:,ifield)  = avfield - sum(pomcav  (:,1:ifield-1),2)
!    write(6,"(A5)") "pomc"

    ! In-cloud sea salt (isscld)
    avfield    = 0.0
    call slabsum(avfield,1,k1,aer_tend(:,:,:,isscld),2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    sscav (:,ifield)  = avfield - sum(sscav  (:,1:ifield-1),2)
!    write(6,"(A5)") "ssc"

    ! In-cloud mineral dust (iducld)
    avfield    = 0.0
    call slabsum(avfield,1,k1,aer_tend(:,:,:,iducld),2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    ducav (:,ifield)  = avfield - sum(ducav  (:,1:ifield-1),2)
!    write(6,"(A5)") "duc"

    ! In-rain sulphate (iso4rai) 
    avfield    = 0.0
    call slabsum(avfield,1,k1,aer_tend(:,:,:,iso4rai),2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    so4rav(:,ifield)  = avfield - sum(so4rav  (:,1:ifield-1),2)
!    write(6,"(A5)") "so4r"

    ! In-rain black carbon (ibcrai)
    avfield    = 0.0
    call slabsum(avfield,1,k1,aer_tend(:,:,:,ibcrai ),2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    bcrav (:,ifield)  = avfield - sum(bcrav  (:,1:ifield-1),2)
 !   write(6,"(A5)") "bcr"

    ! In-rain organic matter (ipomrai)
    avfield    = 0.0
    call slabsum(avfield,1,k1,aer_tend(:,:,:,ipomrai),2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    pomrav(:,ifield)  = avfield - sum(pomrav  (:,1:ifield-1),2)
!    write(6,"(A5)") "pomr"

    ! In-rain sea salt (issrai)
    avfield    = 0.0
!    write(6,"(A15,I3)") "In-rain salt", issrai    
    call slabsum(avfield,1,k1,aer_tend(:,:,:,issrai ),2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
!    write(6,"(A5,E11.4)") 'ssr1', sum(avfield)    
    ssrav(:,ifield)  = avfield - sum(ssrav  (:,1:ifield-1),2)
!    write(6,"(A5, 2E11.4)") "ssr2", sum(ssrav(:,ifield)), sum(sum(ssrav(:,1:ifield-1),2))

    ! In-rain mineral dust (idurai)
    avfield    = 0.0
!    write(6,"(A15,I3,E11.4)") "In-rain dust", idurai, SUM(aer_tend(:,:,:,idurai))
    call slabsum(avfield,1,k1,aer_tend(:,:,:,idurai ),2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
!    write(6,"(A5,E11.4)") 'dur1', sum(avfield)    
    durav(:,ifield)  = avfield - sum(durav  (:,1:ifield-1),2)
!    write(6,"(A5, 2E11.4)") "dur2", sum(ssrav(:,ifield)), sum(sum(ssrav(:,1:ifield-1),2))

    if (ifield == nrfields) then
        Npmn =   Npmn +   Npav/nsamples/ijtot;   Npav = 0.0
       qlpmn =  qlpmn +  qlpav/nsamples/ijtot;  qlpav = 0.0
       qtpmn =  qtpmn +  qtpav/nsamples/ijtot;  qtpav = 0.0

        ncmn =   ncmn +   ncav/nsamples/ijtot;   ncav = 0.0
      so4cmn = so4cmn + so4cav/nsamples/ijtot; so4cav = 0.0
       bccmn =  bccmn +  bccav/nsamples/ijtot;  bccav = 0.0
      pomcmn = pomcmn + pomcav/nsamples/ijtot; pomcav = 0.0
       sscmn =  sscmn +  sscav/nsamples/ijtot;  sscav = 0.0
       ducmn =  ducmn +  ducav/nsamples/ijtot;  ducav = 0.0
      so4rmn = so4rmn + so4rav/nsamples/ijtot; so4rav = 0.0
       bcrmn =  bcrmn +  bcrav/nsamples/ijtot;  bcrav = 0.0
      pomrmn = pomrmn + pomrav/nsamples/ijtot; pomrav = 0.0
       ssrmn =  ssrmn +  ssrav/nsamples/ijtot;  ssrav = 0.0
       durmn =  durmn +  durav/nsamples/ijtot;  durav = 0.0
    end if

    deallocate(avfield)
  end subroutine bulkmicrotend

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
        vars(:,24) =   ncmn(:,iacti)
        vars(:,25) =   ncmn(:,islfc)
        vars(:,26) =   ncmn(:,isedc)
        vars(:,27) =   ncmn(:,iauto)       
        vars(:,28) =   ncmn(:,iaccr)
        vars(:,29) =   ncmn(:,ievpc)

        vars(:,30) = so4cmn(:,iacti)
        vars(:,31) = so4cmn(:,isedc)
        vars(:,32) = so4cmn(:,iauto)
        vars(:,33) = so4cmn(:,iaccr)
        vars(:,34) = so4cmn(:,iscav)
        vars(:,35) = so4cmn(:,ievpc)

        vars(:,36) =  bccmn(:,iacti)
        vars(:,37) =  bccmn(:,isedc)
        vars(:,38) =  bccmn(:,iauto)
        vars(:,39) =  bccmn(:,iaccr)
        vars(:,40) =  bccmn(:,iscav)
        vars(:,41) =  bccmn(:,ievpc)

        vars(:,42) = pomcmn(:,iacti)
        vars(:,43) = pomcmn(:,isedc)
        vars(:,44) = pomcmn(:,iauto)
        vars(:,45) = pomcmn(:,iaccr)
        vars(:,46) = pomcmn(:,iscav)
        vars(:,47) = pomcmn(:,ievpc)

        vars(:,48) =  sscmn(:,iacti)
        vars(:,49) =  sscmn(:,isedc)
        vars(:,50) =  sscmn(:,iauto)
        vars(:,51) =  sscmn(:,iaccr)
        vars(:,52) =  sscmn(:,iscav)
        vars(:,53) =  sscmn(:,ievpc)

        vars(:,54) =  ducmn(:,iacti)
        vars(:,55) =  ducmn(:,isedc)
        vars(:,56) =  ducmn(:,iauto)
        vars(:,57) =  ducmn(:,iaccr)
        vars(:,58) =  ducmn(:,iscav)
        vars(:,59) =  ducmn(:,ievpc)

        vars(:,60) = so4rmn(:,ievpr)
        vars(:,61) = so4rmn(:,isedr)
        vars(:,62) = so4rmn(:,iscav)

        vars(:,63) =  bcrmn(:,ievpr)
        vars(:,64) =  bcrmn(:,isedr)
        vars(:,65) =  bcrmn(:,iscav)

        vars(:,66) = pomrmn(:,ievpr)
        vars(:,67) = pomrmn(:,isedr)
        vars(:,68) = pomrmn(:,iscav)

        vars(:,69) =  ssrmn(:,ievpr)
        vars(:,70) =  ssrmn(:,isedr)
        vars(:,71) =  ssrmn(:,iscav)

        vars(:,72) =  durmn(:,ievpr)
        vars(:,73) =  durmn(:,isedr)
        vars(:,74) =  durmn(:,iscav)

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
               so4rav, so4rmn, &
                bcrav,  bcrmn, &
               pomrav, pomrmn, &
                ssrav,  ssrmn, &
                durav,  durmn  )

    deallocate(precavl    , &
         precav    , &
         precmn    , &
         preccountavl    , &
         preccountav    , &
         preccountmn    , &
         prec_prcavl    , &
         prec_prcav    , &
         prec_prcmn    , &
         cloudcountavl  , &
         cloudcountav    , &
         cloudcountmn    , &
         raincountavl    , &
         raincountav    , &
         raincountmn    , &
         Nrrainavl    , &
         Nrrainav    , &
         Nrrainmn    , &
         qravl    , &
         qrav      , &
         qrmn      , &
         Dvravl    , &
         Dvrav    , &
         Dvrmn)

  end subroutine exitbulkmicrostat

!------------------------------------------------------------------------------!

end module
