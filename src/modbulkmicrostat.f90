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
  use modprecision, only : longint, field_r

implicit none
private
PUBLIC  :: initbulkmicrostat, bulkmicrostat, exitbulkmicrostat, bulkmicrotend
save
!NetCDF variables
  integer,parameter :: nvar = 23
  character(80),dimension(nvar,4) :: ncname
  character(80),dimension(1,4) :: tncname
  real          :: dtav, timeav
  integer(kind=longint):: idtav, itimeav, tnext, tnextwrite
  integer          :: nsamples
  logical          :: lmicrostat = .false.
  integer, parameter      :: nrfields = 5  , &
                 iauto    = 2 , &
                  iaccr    = 3 , &
               ievap    = 4 , &
               ised      = 5
  real, allocatable, dimension(:,:)  :: Npav    , &
               Npmn    , &
               qlpav  , &
               qlpmn  , &
               qtpav  , &
               qtpmn
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
    use modmpi,    only  : myid, mpi_logical, comm3d, mpierr, D_MPI_BCAST
    use modglobal, only  : ifnamopt, fname_options, cexpnr, ifoutput, &
         dtav_glob, timeav_glob, ladaptive, k1, dtmax,btime,tres,lwarmstart,checknamelisterror
    use modstat_nc, only : lnetcdf,define_nc,ncinfo,nctiminfo,writestat_dims_nc
    use modgenstat, only : idtav_prof=>idtav, itimeav_prof=>itimeav,ncid_prof=>ncid
    use modmicrodata,only: imicro, imicro_bulk, imicro_sice, imicro_bulk3 !#sb3
    use modbulkmicrostat3,only:initbulkmicrostat3  ! #sb3
    implicit none
    integer      :: ierr

    namelist/NAMBULKMICROSTAT/ &
    lmicrostat, dtav, timeav
    
    ! #sb3 START
    if (imicro.eq.imicro_bulk3) then
       call initbulkmicrostat3
       return
    endif
    ! #sb3 END

    if ((imicro /=imicro_bulk) .and. (imicro /= imicro_sice)) return

    
    
    dtav  = dtav_glob
    timeav  = timeav_glob
    if(myid==0)then
      open (ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMBULKMICROSTAT,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMBULKMICROSTAT')
      write(6,NAMBULKMICROSTAT)
      close(ifnamopt)
    end if

    call D_MPI_BCAST(lmicrostat,1,0,comm3d,mpierr)
    call D_MPI_BCAST(dtav      ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(timeav    ,1,0,comm3d,mpierr)
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

    allocate(Npav    (k1, nrfields)  , &
       Npmn    (k1, nrfields)  , &
       qlpav    (k1, nrfields)  , &
       qlpmn    (k1, nrfields)  , &
       qtpav    (k1, nrfields)  , &
       qtpmn    (k1, nrfields)  )
    allocate(precavl  (k1)    , &
       precav    (k1)    , &
       precmn    (k1)    , &
       preccountavl  (k1)    , &
       preccountav  (k1)    , &
       preccountmn  (k1)    , &
       prec_prcavl  (k1)    , &
       prec_prcav  (k1)    , &
       prec_prcmn  (k1)    , &
       cloudcountavl  (k1)    , &
       cloudcountav  (k1)    , &
       cloudcountmn  (k1)    , &
       raincountavl  (k1)    , &
       raincountav  (k1)    , &
       raincountmn  (k1)    , &
       Nrrainavl  (k1)    , &
       Nrrainav  (k1)    , &
       Nrrainmn  (k1)    , &
       qravl    (k1)    , &
       qrav    (k1)    , &
       qrmn    (k1)    , &
       Dvravl    (k1)    , &
       Dvrav    (k1)    , &
       Dvrmn    (k1))
    Npmn    = 0.0
    qlpmn    = 0.0
    qtpmn    = 0.0
    precmn    = 0.0
    preccountmn  = 0.0
    prec_prcmn  = 0.0
    cloudcountmn  = 0.0
    raincountmn  = 0.0
    Nrrainmn  = 0.0
    qrmn    = 0.0
    Dvrmn    = 0.0


    if (myid == 0 .and. .not. lwarmstart) then
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
        call nctiminfo(tncname(1,:))
        call ncinfo(ncname( 1,:),'cfrac','Cloud fraction','-','tt')
        call ncinfo(ncname( 2,:),'rainrate','Echo Rain Rate','W/m^2','tt')
        call ncinfo(ncname( 3,:),'preccount','Preccount','W/m^2','tt')
        call ncinfo(ncname( 4,:),'nrrain','nrrain','W/m^2','tt')
        call ncinfo(ncname( 5,:),'raincount','raincount','W/m^2','tt')
        call ncinfo(ncname( 6,:),'precmn','precmn','W/m^2','tt')
        call ncinfo(ncname( 7,:),'dvrmn','dvrmn','W/m^2','tt')
        call ncinfo(ncname( 8,:),'qrmn','qrmn','W/m^2','tt')
        call ncinfo(ncname( 9,:),'npauto','Autoconversion rain drop tendency','#/m3/s','tt')
        call ncinfo(ncname(10,:),'npaccr','Accretion rain drop tendency','#/m3/s','tt')
        call ncinfo(ncname(11,:),'npsed','Sedimentation rain drop tendency','#/m3/s','tt')
        call ncinfo(ncname(12,:),'npevap','Evaporation rain drop tendency','#/m3/s','tt')
        call ncinfo(ncname(13,:),'nptot','Total rain drop tendency','#/m3/s','tt')
        call ncinfo(ncname(14,:),'qrpauto','Autoconversion rain water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(15,:),'qrpaccr','Accretion rain water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(16,:),'qrpsed','Sedimentation rain water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(17,:),'qrpevap','Evaporation rain water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(18,:),'qrptot','Total rain water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(19,:),'qtpauto','Autoconversion total water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(20,:),'qtpaccr','Accretion total water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(21,:),'qtpsed','Sedimentation total water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(22,:),'qtpevap','Evaporation total water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(23,:),'qtptot','Total total water content tendency','kg/kg/s','tt')
        call define_nc( ncid_prof, NVar, ncname)
      end if

   end if

  end subroutine initbulkmicrostat

!------------------------------------------------------------------------------!
!> General routine, does the timekeeping
  subroutine bulkmicrostat
    use modglobal,    only  : rk3step, timee, dt_lim
    use modmicrodata,only: imicro, imicro_bulk3             ! #sb3
    implicit none
    
    if (imicro.eq.imicro_bulk3) return  ! #sb3 treated separately
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
    use modmpi,    only  :  mpi_sum, comm3d, mpierr, D_MPI_ALLREDUCE
    use modglobal,    only  : i1, j1, k1, ijtot
    use modmicrodata,  only  : qr,precep,Dvr,Nr,epscloud,epsqr,epsprec,imicro,imicro_bulk
    use modfields,  only  : ql0
    implicit none

    integer :: k

    precav = 0.0
    preccountav = 0.0
    prec_prcav = 0.0
    cloudcountav = 0.0
    raincountav = 0.0
    Nrrainav = 0.0
    qrav = 0.0
    Dvrav = 0.0

    do k = 1,k1
      cloudcountavl(k)  = count(ql0     (2:i1,2:j1,k) > epscloud)
      raincountavl (k)  = count(qr      (2:i1,2:j1,k) > epsqr)
      preccountavl (k)  = count(precep  (2:i1,2:j1,k) > epsprec)
      prec_prcavl  (k)  = sum  (precep  (2:i1,2:j1,k) , precep(2:i1,2:j1,k) > epsprec)
      Nrrainavl    (k)  = sum  (Nr      (2:i1,2:j1,k))
      precavl      (k)  = sum  (precep  (2:i1,2:j1,k))
      qravl        (k)  = sum  (qr      (2:i1,2:j1,k))
      if (imicro==imicro_bulk) then
        Dvravl     (k)  = sum  (Dvr     (2:i1,2:j1,k) , qr(2:i1,2:j1,k) > epsqr)
      end if
    end do

    call D_MPI_ALLREDUCE(cloudcountavl,cloudcountav,k1,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(raincountavl ,raincountav ,k1,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(preccountavl ,preccountav ,k1,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(prec_prcavl  ,prec_prcav  ,k1,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(Dvravl       ,Dvrav       ,k1,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(Nrrainavl    ,Nrrainav    ,k1,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(precavl      ,precav      ,k1,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(qravl        ,qrav        ,k1,MPI_SUM,comm3d,mpierr)

    cloudcountmn  = cloudcountmn  +  cloudcountav  /ijtot
    raincountmn  = raincountmn  +  raincountav  /ijtot
    preccountmn  = preccountmn  +  preccountav  /ijtot
    prec_prcmn  = prec_prcmn  +  prec_prcav  /ijtot
    Dvrmn    = Dvrmn    +  Dvrav    /ijtot
    Nrrainmn  = Nrrainmn  +  Nrrainav  /ijtot
    precmn    = precmn  +  precav    /ijtot
    qrmn    = qrmn    +  qrav    /ijtot

  end subroutine dobulkmicrostat

!------------------------------------------------------------------------------!
!> Performs the calculations for the tendencies etc.
  subroutine bulkmicrotend
    use modmpi,    only  : slabsum
    use modglobal,    only  : rk3step, timee, dt_lim, k1, ih, i1, jh, j1, ijtot
    use modfields,    only  : qtp
    use modmicrodata,  only  : qrp, Nrp
    implicit none

    real(field_r), dimension(:), allocatable  :: avfield
    integer        :: ifield = 0

    if (.not. lmicrostat)  return
    if (rk3step /= 3)  return
    if (timee == 0)    return
    if (timee < tnext .and. timee < tnextwrite) then
      dt_lim  = minval((/dt_lim, tnext - timee, tnextwrite - timee/))
      return
    end if
!    tnext = tnext+dtav

    allocate(avfield(k1))

    ifield    = mod(ifield, nrfields) + 1

    avfield    = 0.0
    call slabsum(avfield  ,1,k1,Nrp  ,2,i1,2,j1,1,k1,2,i1,2,j1,1,k1)
    Npav(:,ifield)  = avfield - sum(Npav  (:,1:ifield-1),2)

    avfield    = 0.0
    call slabsum(avfield  ,1,k1,qrp  ,2,i1,2,j1,1,k1,2,i1,2,j1,1,k1)
    qlpav(:,ifield) = avfield - sum(qlpav  (:,1:ifield-1),2)

    avfield    = 0.0
    call slabsum(avfield  ,1,k1,qtp  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    qtpav(:,ifield) = avfield - sum(qtpav  (:,1:ifield-1),2)

    if (ifield == nrfields) then
      Npmn    = Npmn  + Npav  /nsamples/ijtot
      qlpmn    = qlpmn  + qlpav /nsamples/ijtot
      qtpmn    = qtpmn  + qtpav /nsamples/ijtot
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
              rlv, zf
    use modfields,    only  : presf,rhof
    use modstat_nc, only: lnetcdf, writestat_nc
    use modgenstat, only: ncid_prof=>ncid,nrec_prof=>nrec

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
                prec_prcmn      = prec_prcmn    /nsamples
                Dvrmn           = Dvrmn         /nsamples
                Nrrainmn        = Nrrainmn      /nsamples
                precmn          = precmn        /nsamples
                qrmn            = qrmn          /nsamples

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
    write(ifoutput,'(I4,F10.2,F7.1,5E13.5)') &
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
    write(ifoutput,'(I4,F10.2,F7.1,5E13.5)') &
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
        vars(:, 1) = cloudcountmn
        vars(:, 2) = prec_prcmn  (:)*rhof(:)*rlv
        vars(:, 3) = preccountmn  (:)
        vars(:, 4) = Nrrainmn  (:)
        vars(:, 5) = raincountmn  (:)
        vars(:, 6) = precmn    (:)*rhof(:)*rlv
        vars(:, 7) = Dvrmn    (:)
        vars(:, 8) = qrmn    (:)
        vars(:, 9) =Npmn    (:,iauto)
        vars(:,10) =Npmn    (:,iaccr)
        vars(:,11) =Npmn    (:,ised)
        vars(:,12) =Npmn    (:,ievap)
        do k=1,k1
        vars(k,13) =sum(Npmn  (k,2:nrfields))
        enddo
        vars(:,14) =qlpmn    (:,iauto)
        vars(:,15) =qlpmn    (:,iaccr)
        vars(:,16) =qlpmn    (:,ised)
        vars(:,17) =qlpmn    (:,ievap)
        do k=1,k1
        vars(k,18) =sum(qlpmn  (k,2:nrfields))
        enddo
        vars(:,19) =qtpmn    (:,iauto)
        vars(:,20) =qtpmn    (:,iaccr)
        vars(:,21) =qtpmn    (:,ised)
        vars(:,22) =qtpmn    (:,ievap)
        do k=1,k1
        vars(k,23) =sum(qtpmn  (k,2:nrfields))
        enddo
        call writestat_nc(ncid_prof,nvar,ncname,vars(1:kmax,:),nrec_prof,kmax)
      end if

    end if

    cloudcountmn    = 0.0
    raincountmn    = 0.0
    preccountmn    = 0.0
    prec_prcmn    = 0.0
    Dvrmn      = 0.0
    Nrrainmn    = 0.0
    precmn      = 0.0
    qrmn      = 0.0
    Npmn      = 0.0
    qlpmn      = 0.0
    qtpmn      = 0.0

  end subroutine writebulkmicrostat

!------------------------------------------------------------------------------!

  subroutine exitbulkmicrostat
    use modmicrodata,only: imicro, imicro_bulk3     ! #sb3
    use modbulkmicrostat3, only:exitbulkmicrostat3  ! #sb3
    implicit none
    ! #sb3 START
    if (imicro.eq.imicro_bulk3) then
       call exitbulkmicrostat3
       return
    endif
    ! #sb3 END
    
    if (.not. lmicrostat)  return

    deallocate(Npav      , &
         Npmn      , &
         qlpav    , &
         qlpmn    , &
         qtpav    , &
         qtpmn    )
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
