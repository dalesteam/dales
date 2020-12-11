!> \file modcanstat.f90
!!  Calculates the canopy statistics XPB


!>
!!  Calculates the canopy statistics
!>
!! Profiles of the radiative statistics. Written to radstat.expnr
!! If netcdf is true, this module also writes in the profiles.expnr.nc output
!!  \author Stephan de Roode, TU Delft
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
module modcanstat

  use modglobal, only : longint

implicit none
!private
PUBLIC :: initcanstat, canstat, exitcanstat
save
!NetCDF variables
  integer,parameter :: nvar = 17
  character(80),dimension(nvar,4) :: ncname

  real    :: dtav, timeav
  integer(kind=longint) :: idtav,itimeav,tnext,tnextwrite
  integer :: nsamples
  logical :: lstat= .false. !< switch to enable the canopy statistics (on/off)

!     ------

!   --------------
  real, allocatable :: shcanav  (:)   
  real, allocatable :: lecanav  (:)  
  real, allocatable :: fco2canav(:)       
  real, allocatable :: sthetaav (:)  
  real, allocatable :: sqtav    (:) 
  real, allocatable :: sco2av   (:)  
  real, allocatable :: t_leafshadav   (:)  
  real, allocatable :: t_leafsunav    (:)  
  real, allocatable :: gcc_leafshadav (:)  
  real, allocatable :: gcc_leafsunav  (:)  
  real, allocatable :: absSWleaf_shadav(:)
  real, allocatable :: absSWleaf_allsunav(:)
  real, allocatable :: absSWlayerav(:)
  real, allocatable :: cfSLav   (:)  
  
  real, allocatable :: shcanmn  (:)   
  real, allocatable :: lecanmn  (:)  
  real, allocatable :: fco2canmn(:)       
  real, allocatable :: sthetamn (:)  
  real, allocatable :: sqtmn    (:) 
  real, allocatable :: sco2mn   (:)  
  real, allocatable :: t_leafshadmn   (:)  
  real, allocatable :: t_leafsunmn    (:)  
  real, allocatable :: gcc_leafshadmn (:)  
  real, allocatable :: gcc_leafsunmn  (:)  
  real, allocatable :: absSWleaf_shadmn(:)  
  real, allocatable :: absSWleaf_allsunmn(:)  
  real, allocatable :: absSWlayermn(:)  
  real, allocatable :: cfSLmn   (:)  

contains
!> Initialization routine, reads namelists and inits variables
  subroutine initcanstat
    use modmpi,    only : myid,mpierr, comm3d,my_real, mpi_logical
    use modglobal, only : dtmax, ifnamopt,fname_options, ifoutput,&
                          cexpnr,dtav_glob,timeav_glob,ladaptive,dt_lim,btime,tres
    use modstat_nc, only : lnetcdf,define_nc,ncinfo
    use modgenstat, only : idtav_prof=>idtav, itimeav_prof=>itimeav,ncid_prof=>ncid
    use modcanopy, only : ncanopy

    implicit none

    integer ierr
    namelist/NAMCANSTAT/ &
    dtav,timeav,lstat

    dtav=dtav_glob;timeav=timeav_glob
    lstat = .false.

    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMCANSTAT,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions NAMCANSTAT'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMCANSTAT'
      endif
      write(6 ,NAMCANSTAT)
      close(ifnamopt)
    end if

    call MPI_BCAST(timeav     ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(dtav       ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(lstat   ,1,MPI_LOGICAL,0,comm3d,mpierr)
    idtav = dtav/tres
    itimeav = timeav/tres

    tnext      = idtav   +btime
    tnextwrite = itimeav +btime
    nsamples = itimeav/idtav


    if(.not.(lstat)) return
    dt_lim = min(dt_lim,tnext)

    if (abs(timeav/dtav-nsamples)>1e-4) then
      stop 'timeav must be a integer multiple of dtav'
    end if
    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'dtav should be a integer multiple of dtmax'
    end if
    
    allocate(shcanav  (ncanopy))
    allocate(lecanav  (ncanopy))
    allocate(fco2canav(ncanopy))
    allocate(sthetaav (ncanopy))
    allocate(sqtav    (ncanopy))
    allocate(sco2av   (ncanopy))
    allocate(t_leafshadav   (ncanopy))
    allocate(t_leafsunav    (ncanopy))
    allocate(gcc_leafshadav (ncanopy))
    allocate(gcc_leafsunav  (ncanopy))
    allocate(absSWleaf_shadav    (ncanopy))
    allocate(absSWleaf_allsunav  (ncanopy))
    allocate(absSWlayerav  (ncanopy))
    allocate(cfSLav   (ncanopy))

    allocate(shcanmn  (ncanopy))
    allocate(lecanmn  (ncanopy))
    allocate(fco2canmn(ncanopy))
    allocate(sthetamn (ncanopy))
    allocate(sqtmn    (ncanopy))
    allocate(sco2mn   (ncanopy))
    allocate(t_leafshadmn   (ncanopy))
    allocate(t_leafsunmn    (ncanopy))
    allocate(gcc_leafshadmn (ncanopy))
    allocate(gcc_leafsunmn  (ncanopy))
    allocate(absSWleaf_shadmn    (ncanopy))
    allocate(absSWleaf_allsunmn  (ncanopy))
    allocate(absSWlayermn  (ncanopy))
    allocate(cfSLmn   (ncanopy))

    shcanmn   = 0.0
    lecanmn   = 0.0
    fco2canmn = 0.0
    sthetamn  = 0.0
    sqtmn     = 0.0
    sco2mn    = 0.0
    t_leafshadmn   = 0.0
    t_leafsunmn    = 0.0
    gcc_leafshadmn = 0.0
    gcc_leafsunmn  = 0.0
    absSWleaf_shadmn  = 0.0
    absSWleaf_allsunmn  = 0.0
    absSWlayermn  = 0.0
    cfSLmn    = 0.0

    if(myid==0)then
      open (ifoutput,file='canstat.'//cexpnr,status='replace')
      close (ifoutput)
    end if
    if (lnetcdf) then
      idtav = idtav_prof
      itimeav = itimeav_prof
      tnext      = idtav+btime
      tnextwrite = itimeav+btime
      nsamples = itimeav/idtav

      if (myid==0) then
        call ncinfo(ncname( 1,:),'padf','Plant area density at full level','m^2/m^3','tt')
        call ncinfo(ncname( 2,:),'pai','Plant area index at full level','m^2/m^2','tt')
        call ncinfo(ncname( 3,:),'paih','Plant area index at half level','m^2/m^2','tt')
        call ncinfo(ncname( 4,:),'cfSL','Fraction of sunlit leaves in canopy','-','tt')
        call ncinfo(ncname( 5,:),'shcan','Canopy sensible heat source','W/m^3','tt')
        call ncinfo(ncname( 6,:),'lecan','Canopy latent heat source','W/m^3','tt')
        call ncinfo(ncname( 7,:),'fco2can','Canopy CO2 source','mg C/(s m^3)','tt')
        call ncinfo(ncname( 8,:),'sthetamn','Canopy temperature tendency','K/s','tt')
        call ncinfo(ncname( 9,:),'sqtmn','Canopy specific humidity tendency','Kg_w/Kg_a/s','tt')
        call ncinfo(ncname( 10,:),'sco2mn','Canopy CO2 tendency','ppb/s','tt')
        call ncinfo(ncname( 11,:),'t_leafshadmn','Shaded leave temperature','K','tt')
        call ncinfo(ncname( 12,:),'t_leafsunmn','Sunlit leaf temperature','K','tt')
        call ncinfo(ncname( 13,:),'gcc_leafshadmn','Shaded leave stomatal CO2 conductance','','tt')
        call ncinfo(ncname( 14,:),'gcc_leafsunmn','Sunlit leave stomatal CO2 conductance','','tt')
        call ncinfo(ncname( 15,:),'absSWleaf_shadmn','Absorbed SW by sahded leaves','W/m2','tt')
        call ncinfo(ncname( 16,:),'absSWleaf_allsunmn','Absorbed SW by sunlit leaves','W/m2','tt')
        call ncinfo(ncname( 17,:),'absSWlayermn','Absorbed SW in the layer','W/m2','tt') ! weighing sunlit and shaded leave fraction

        call define_nc( ncid_prof, NVar, ncname)
      end if

   end if

  end subroutine initcanstat
!> General routine, does the timekeeping
  subroutine canstat
    use modglobal, only : rk3step,timee,dt_lim
    implicit none
    if (.not. lstat) return
    if (rk3step/=3) return
    if(timee<tnext .and. timee<tnextwrite) then
      dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))
      return
    end if
    if (timee>=tnext) then
      tnext = tnext+idtav
      call do_canstat
    end if
    if (timee>=tnextwrite) then
      tnextwrite = tnextwrite+itimeav
      call writecanstat
    end if
    dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))

  end subroutine canstat

!> Calculates the statistics
  subroutine do_canstat

    use modmpi,    only : slabsum
    use modglobal, only : ijtot,i1,j1,ih,jh
    use modcanopy, only : S_theta,S_qt,S_co2,sh_can,le_can,Fco2_can,ncanopy,cfSL,t_leafshad,&
                          t_leafsun,gcc_leafshad,gcc_leafsun,absSWleaf_shad,absSWleaf_allsun,absSWlayer
    implicit none

    shcanav   = 0.
    lecanav   = 0.
    fco2canav = 0.
    sthetaav  = 0.
    sqtav     = 0.
    sco2av    = 0.
    t_leafshadav    = 0.
    t_leafsunav     = 0.
    gcc_leafshadav  = 0.
    gcc_leafsunav   = 0.
    absSWleaf_shadav   = 0.
    absSWleaf_allsunav   = 0.
    absSWlayerav   = 0.
    cfSLav    = 0.

    call slabsum(shcanav  ,1,ncanopy, sh_can  ,2-ih,i1+ih,2-jh,j1+jh,1,ncanopy,2,i1,2,j1,1,ncanopy)
    call slabsum(lecanav  ,1,ncanopy, le_can  ,2-ih,i1+ih,2-jh,j1+jh,1,ncanopy,2,i1,2,j1,1,ncanopy)
    call slabsum(fco2canav,1,ncanopy, Fco2_can,2-ih,i1+ih,2-jh,j1+jh,1,ncanopy,2,i1,2,j1,1,ncanopy)
    call slabsum(sthetaav ,1,ncanopy, S_theta ,2-ih,i1+ih,2-jh,j1+jh,1,ncanopy,2,i1,2,j1,1,ncanopy)
    call slabsum(sqtav    ,1,ncanopy, S_qt    ,2-ih,i1+ih,2-jh,j1+jh,1,ncanopy,2,i1,2,j1,1,ncanopy)
    call slabsum(sco2av   ,1,ncanopy, S_co2   ,2-ih,i1+ih,2-jh,j1+jh,1,ncanopy,2,i1,2,j1,1,ncanopy)
    call slabsum(t_leafshadav    ,1,ncanopy, t_leafshad   ,2-ih,i1+ih,2-jh,j1+jh,1,ncanopy,2,i1,2,j1,1,ncanopy)
    call slabsum(t_leafsunav     ,1,ncanopy, t_leafsun    ,2-ih,i1+ih,2-jh,j1+jh,1,ncanopy,2,i1,2,j1,1,ncanopy)
    call slabsum(gcc_leafshadav  ,1,ncanopy, gcc_leafshad ,2-ih,i1+ih,2-jh,j1+jh,1,ncanopy,2,i1,2,j1,1,ncanopy)
    call slabsum(gcc_leafsunav   ,1,ncanopy, gcc_leafsun  ,2-ih,i1+ih,2-jh,j1+jh,1,ncanopy,2,i1,2,j1,1,ncanopy)
    call slabsum(absSWleaf_shadav     ,1,ncanopy, absSWleaf_shad  (:,:,1:ncanopy)  ,2-ih,i1+ih,2-jh,j1+jh,1,ncanopy,2,i1,2,j1,1,ncanopy)
    call slabsum(absSWleaf_allsunav   ,1,ncanopy, absSWleaf_allsun(:,:,1:ncanopy)  ,2-ih,i1+ih,2-jh,j1+jh,1,ncanopy,2,i1,2,j1,1,ncanopy)
    call slabsum(absSWlayerav         ,1,ncanopy, absSWlayer      (:,:,1:ncanopy)  ,2-ih,i1+ih,2-jh,j1+jh,1,ncanopy,2,i1,2,j1,1,ncanopy)

    cfSLav = cfSL  ! no slab average needed as it only depends on time, not on space
 !    ADD SLAB AVERAGES TO TIME MEAN

    shcanmn      = shcanmn     + shcanav     / ijtot
    lecanmn      = lecanmn     + lecanav     / ijtot
    fco2canmn    = fco2canmn   + fco2canav   / ijtot
    sthetamn     = sthetamn    + sthetaav    / ijtot
    sqtmn        = sqtmn       + sqtav       / ijtot
    sco2mn       = sco2mn      + sco2av      / ijtot
    t_leafshadmn        = t_leafshadmn         + t_leafshadav        / ijtot
    t_leafsunmn         = t_leafsunmn          + t_leafsunav         / ijtot
    gcc_leafshadmn      = gcc_leafshadmn       + gcc_leafshadav      / ijtot
    gcc_leafsunmn       = gcc_leafsunmn        + gcc_leafsunav       / ijtot
    absSWleaf_shadmn   = absSWleaf_shadmn    + absSWleaf_shadav   / ijtot
    absSWleaf_allsunmn = absSWleaf_allsunmn  + absSWleaf_allsunav / ijtot
    absSWlayermn = absSWlayermn  + absSWlayerav / ijtot
    cfSLmn       = cfSLmn      + cfSLav

  end subroutine do_canstat

!> Write the statistics to file
  subroutine writecanstat
      use modmpi,    only : myid
      use modglobal, only : cexpnr,ifoutput,zf,rtimee,zh
      use modstat_nc, only: lnetcdf, writestat_nc
      use modgenstat, only: ncid_prof=>ncid,nrec_prof=>nrec
      use modcanopy, only : ncanopy,pai,padf,paih
      implicit none
      real,dimension(ncanopy,nvar) :: vars
      integer nsecs, nhrs, nminut,k


      nsecs   = nint(rtimee)
      nhrs    = int(nsecs/3600)
      nminut  = int(nsecs/60)-nhrs*60
      nsecs   = mod(nsecs,60)

      sthetamn   = sthetamn   /nsamples
      sqtmn      = sqtmn      /nsamples
      sco2mn     = sco2mn     /nsamples
      t_leafshadmn   = t_leafshadmn    /nsamples
      t_leafsunmn    = t_leafsunmn     /nsamples
      gcc_leafshadmn = gcc_leafshadmn  /nsamples
      gcc_leafsunmn  = gcc_leafsunmn   /nsamples
      absSWleaf_shadmn   = absSWleaf_shadmn     /nsamples
      absSWleaf_allsunmn = absSWleaf_allsunmn   /nsamples
      absSWlayermn       = absSWlayermn         /nsamples
      shcanmn    = shcanmn    /nsamples
      lecanmn    = lecanmn    /nsamples
      fco2canmn  = fco2canmn  /nsamples
      cfSLmn     = cfSLmn     /nsamples
  !     ----------------------
  !     2.0  write the fields
  !           ----------------

    if(myid==0)then
      open (ifoutput,file='canstat.'//cexpnr,position='append')
      write(ifoutput,'(//A,/A,F5.0,A,I4,A,I2,A,I2,A)') &
      '#--------------------------------------------------------'      &
      ,'#',(timeav),'--- AVERAGING TIMESTEP --- '      &
      ,nhrs,':',nminut,':',nsecs      &
      ,'   HRS:MIN:SEC AFTER INITIALIZATION '
      write (ifoutput,'(A/3A/3A)') &
          '#--------------------------------------------------------------------------' &
          ,'#LEV     HEIGHT    HEIGHT_H    PAD         PAI          PAI_H       cfSL_H'&
          ,'      SH_CAN       LE_CAN       FCO2_CAN      S_THETA       S_QT       S_CO2     ' &
          ,'T_LSHAD   T_LSUN   GCC_LSHAD    GCC_LSUN  absSW_LSHAD absSW_LSUN absSW_LAYER' &
          ,'#          (M)       (M)    (M^2/M^3)    (M^2/M^2)    (M^2/M^2)      (-)'&
          ,'        (W/M^3)      (W/M^3)  (mg C / S / M^3)  (K/S)      (KG/(KG S)   (PPB/S)  '&
          ,'   (K)     (K)       ()            ()      (W/M^2)   (W/M^2)   (W/M^2)    '
      do k=1,ncanopy
        write(ifoutput,'(I4,2F10.2,10E13.4,2F9.3,2E13.4,3F9.3)') &
            k,zf(k),zh(k),padf(k),pai(k),paih(k),&
            cfSLmn(k),            & 
            shcanmn(k),           &
            lecanmn(k),           &
            fco2canmn(k),         &
            sthetamn(k),          & 
            sqtmn(k),             &  
            sco2mn(k),            &  
            t_leafshadmn(k),      &
            t_leafsunmn(k),       &
            gcc_leafshadmn(k),    &
            gcc_leafsunmn(k),     &
            absSWleaf_shadmn(k),  &
            absSWleaf_allsunmn(k),&
            absSWlayermn(k)      
      end do
      close (ifoutput)

      if (lnetcdf) then
        vars(:, 1) = padf
        vars(:, 2) = pai
        vars(:, 3) = paih
        vars(:, 4) = cfSLmn
        vars(:, 5) = shcanmn
        vars(:, 6) = lecanmn
        vars(:, 7) = fco2canmn
        vars(:, 8) = sthetamn
        vars(:, 9) = sqtmn
        vars(:, 10) = sco2mn
        vars(:, 11) = t_leafshadmn
        vars(:, 12) = t_leafsunmn
        vars(:, 13) = gcc_leafshadmn
        vars(:, 14) = gcc_leafsunmn
        vars(:, 15) = absSWleaf_shadmn
        vars(:, 16) = absSWleaf_allsunmn
        vars(:, 17) = absSWlayermn
       call writestat_nc(ncid_prof,nvar,ncname,vars(1:ncanopy,:),nrec_prof,ncanopy)
      end if
    end if !

    shcanmn   = 0.0
    cfSLmn    = 0.0
    lecanmn   = 0.0
    fco2canmn = 0.0
    sthetamn  = 0.0
    sqtmn    = 0.0
    sco2mn    = 0.0
    t_leafshadmn    = 0.0
    t_leafsunmn     = 0.0
    gcc_leafshadmn  = 0.0
    gcc_leafsunmn   = 0.0
    absSWleaf_shadmn   = 0.0
    absSWleaf_allsunmn  = 0.0
    absSWlayermn  = 0.0

  end subroutine writecanstat

!> Cleans up after the run
  subroutine exitcanstat
    implicit none

    !deallocate variables that are needed in modcaniation

    if(.not.(lstat)) return

    deallocate(shcanav,lecanav,fco2canav,sthetaav,sqtav,sco2av)
    deallocate(t_leafshadav,t_leafsunav,gcc_leafshadav,gcc_leafsunav,absSWleaf_shadav,absSWleaf_allsunav,absSWlayerav)
    deallocate(shcanmn,lecanmn,fco2canmn,sthetamn,sqtmn,sco2mn,cfSLmn)
    deallocate(t_leafshadmn,t_leafsunmn,gcc_leafshadmn,gcc_leafsunmn,absSWleaf_shadmn,absSWleaf_allsunmn,absSWlayermn)

  end subroutine exitcanstat


end module modcanstat
