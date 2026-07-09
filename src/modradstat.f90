!> \file modradstat.f90
!!  Calculates the radiative statistics


!>
!!  Calculates the radiative statistics
!>
!! This module is activated by setting lstat=.true. in the NAMRADSTAT namelist.
!! Profiles of the radiative statistics are written to radstat.expnr
!! If lnetcdf is true, this module also writes in the profiles.expnr.nc output.
!!  \author Stephan de Roode, TU Delft
!
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
module modradstat

  use modglobal, only : longint
  use modprecision, only: field_r
  use modlogging, only: finish

implicit none
character(len=*), parameter :: modname = 'modradstat'
!private
PUBLIC :: initradstat, radstat, exitradstat
save
!NetCDF variables
  integer,parameter :: nvar = 14
  character(80),dimension(nvar,4) :: ncname

  real    :: dtav, timeav
  integer(kind=longint) :: idtav,itimeav,tnext,tnextwrite
  integer :: nsamples
  logical :: lstat= .false. !< switch to enable the radiative statistics (on/off)
  logical :: lradclearair= .false. !< switch to enable the radiative statistics (on/off)

!     ------

!   --------------
  real(field_r), allocatable :: thltendav(:)
  real(field_r), allocatable :: thllwtendav(:)
  real(field_r), allocatable :: thlswtendav(:)
  real(field_r), allocatable :: thllwtendcaav(:)
  real(field_r), allocatable :: thlswtendcaav(:)
  real(field_r), allocatable :: lwuav(:)
  real(field_r), allocatable :: lwdav(:)
  real(field_r), allocatable :: swdav(:)
  real(field_r), allocatable :: swdirav(:)
  real(field_r), allocatable :: swdifav(:)
  real(field_r), allocatable :: swuav(:)
  real(field_r), allocatable :: lwucaav(:)
  real(field_r), allocatable :: lwdcaav(:)
  real(field_r), allocatable :: swdcaav(:)
  real(field_r), allocatable :: swucaav(:)

!
  real, allocatable :: thltendmn(:)
  real, allocatable :: thllwtendmn(:)
  real, allocatable :: thlswtendmn(:)
  real, allocatable :: thllwtendcamn(:)
  real, allocatable :: thlswtendcamn(:)
  real, allocatable :: lwumn(:)
  real, allocatable :: lwdmn(:)
  real, allocatable :: swdmn(:)
  real, allocatable :: swdirmn(:)
  real, allocatable :: swdifmn(:)
  real, allocatable :: swumn(:)
  real, allocatable :: lwucamn(:)
  real, allocatable :: lwdcamn(:)
  real, allocatable :: swdcamn(:)
  real, allocatable :: swucamn(:)
  real, allocatable :: thlradlsmn(:)

contains
!> Initialization routine, reads namelists and inits variables
  subroutine initradstat
    use modmpi,    only : myid,mpierr, comm3d, D_MPI_BCAST
    use modglobal, only : dtmax, k1, ifnamopt,fname_options, ifoutput,&
                          cexpnr,dtav_glob,timeav_glob,ladaptive,dt_lim,btime,tres,lwarmstart,checknamelisterror
    use modstat_nc, only : lnetcdf,define_nc,ncinfo
    use fortran_support, only: nnml_output
    use modgenstat, only : idtav_prof=>idtav, itimeav_prof=>itimeav,ncid_prof=>ncid

    implicit none

    character(len=*), parameter :: routine = modname//'/initradstat'

    integer ierr
    namelist/NAMRADSTAT/ &
    dtav,timeav,lstat,lradclearair

    dtav=dtav_glob;timeav=timeav_glob

    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMRADSTAT,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMRADSTAT')
      write(nnml_output ,NAMRADSTAT)
      close(ifnamopt)
    end if

    call D_MPI_BCAST(timeav      ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(dtav        ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(lstat       ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(lradclearair,1,0,comm3d,mpierr)
    idtav = int(dtav / tres, kind=kind(idtav))
    itimeav = int(timeav / tres, kind=kind(itimeav))

    tnext      = idtav   +btime
    tnextwrite = itimeav +btime
    nsamples = int(itimeav / idtav)


    if(.not.(lstat)) return
    dt_lim = min(dt_lim,tnext)

    if (abs(timeav/dtav-nsamples)>1e-4) then
      call finish(routine, 'timeav must be a integer multiple of dtav')
    end if
    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      call finish(routine, 'dtav should be a integer multiple of dtmax')
    end if

    allocate(lwuav(k1))
    allocate(lwdav(k1))
    allocate(swdav(k1))
    allocate(swdirav(k1))
    allocate(swdifav(k1))
    allocate(swuav(k1))
    allocate(lwucaav(k1))
    allocate(lwdcaav(k1))
    allocate(swdcaav(k1))
    allocate(swucaav(k1))
    allocate(thllwtendav(k1))
    allocate(thltendav(k1))
    allocate(thlswtendav(k1))
    allocate(thllwtendcaav(k1))
    allocate(thlswtendcaav(k1))

    allocate(lwumn(k1))
    allocate(lwdmn(k1))
    allocate(swdmn(k1))
    allocate(swdirmn(k1))
    allocate(swdifmn(k1))
    allocate(swumn(k1))
    allocate(lwucamn(k1))
    allocate(lwdcamn(k1))
    allocate(swdcamn(k1))
    allocate(swucamn(k1))
    allocate(thllwtendmn(k1))
    allocate(thltendmn(k1))
    allocate(thlswtendmn(k1))
    allocate(thlradlsmn(k1))
    allocate(thllwtendcamn(k1))
    allocate(thlswtendcamn(k1))

    lwumn = 0.0
    lwdmn = 0.0
    swdmn = 0.0
    swdirmn = 0.0
    swdifmn = 0.0
    swumn = 0.0
    lwucamn = 0.0
    lwdcamn = 0.0
    swdcamn = 0.0
    swucamn = 0.0
    thltendmn = 0.0
    thllwtendmn = 0.0
    thlswtendmn = 0.0
    thlradlsmn  = 0.0
    thllwtendcamn = 0.0
    thlswtendcamn = 0.0

    if(myid==0 .and. .not. lwarmstart)then
      open (ifoutput,file='radstat.'//cexpnr,status='replace')
      close (ifoutput)
    end if
    if (lnetcdf) then
      idtav = idtav_prof
      itimeav = itimeav_prof
      tnext      = idtav+btime
      tnextwrite = itimeav+btime
      nsamples = int(itimeav / idtav)

      if (myid==0) then
        call ncinfo(ncname( 1,:),'thltend','Total radiative tendency','K/s','tt')
        call ncinfo(ncname( 2,:),'thllwtend','Long wave radiative tendency','K/s','tt')
        call ncinfo(ncname( 3,:),'thlswtend','Short wave radiative tendency','K/s','tt')
        call ncinfo(ncname( 4,:),'thlradls','Prescribed large scale radiative tendency','K/s','tt')
        call ncinfo(ncname( 5,:),'lwu','Long wave upward radiative flux','W/m^2','mt')
        call ncinfo(ncname( 6,:),'lwd','Long wave downward radiative flux','W/m^2','mt')
        call ncinfo(ncname( 7,:),'swu','Short wave upward radiative flux','W/m^2','mt')
        call ncinfo(ncname( 8,:),'swd','Short wave downward radiative flux','W/m^2','mt')
        call ncinfo(ncname( 9,:),'lwuca','Long wave clear air upward radiative flux','W/m^2','mt')
        call ncinfo(ncname(10,:),'lwdca','Long wave clear air downward radiative flux','W/m^2','mt')
        call ncinfo(ncname(11,:),'swuca','Short wave clear air upward radiative flux','W/m^2','mt')
        call ncinfo(ncname(12,:),'swdca','Short wave clear air downward radiative flux','W/m^2','mt')
        call ncinfo(ncname(13,:),'thllwtendca','Long wave clear air radiative tendency','K/s','tt')
        call ncinfo(ncname(14,:),'thlswtendca','Short wave clear air radiative tendency','K/s','tt')

        call define_nc( ncid_prof, NVar, ncname)
      end if

   end if

   !$acc enter data copyin(lwuav, lwdav, swdav, swdirav, swdifav, swuav, lwucaav, lwdcaav, swdcaav, swucaav, &
   !$acc&                  thllwtendav, thltendav,thlswtendav, thllwtendcaav, thlswtendcaav, &
   !$acc&                  lwumn, lwdmn, swdmn, swdirmn, swdifmn, swumn, lwucamn, lwdcamn, swdcamn, swucamn, &
   !$acc&                  thllwtendmn, thltendmn, thlswtendmn, thlradlsmn, thllwtendcamn, thlswtendcamn)
   !$omp target enter data map(to:lwuav,lwdav,swdav,swdirav,swdifav,swuav,&
   !$omp lwucaav,lwdcaav,swdcaav,swucaav,thllwtendav,thltendav,&
   !$omp thlswtendav,thllwtendcaav,thlswtendcaav,lwumn,lwdmn,swdmn,&
   !$omp swdirmn,swdifmn,swumn,lwucamn,lwdcamn,swdcamn,swucamn,&
   !$omp thllwtendmn,thltendmn,thlswtendmn,thlradlsmn,thllwtendcamn,&
   !$omp thlswtendcamn)

  end subroutine initradstat
!> General routine, does the timekeeping
  subroutine radstat
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
      call do_radstat
    end if
    if (timee>=tnextwrite) then
      tnextwrite = tnextwrite+itimeav
      call writeradstat
    end if
    dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))

  end subroutine radstat

!> Calculates the statistics
  subroutine do_radstat

    use modmpi,    only :  slabsum
    use modglobal, only : kmax,ijtot,cp,dzf,i1,j1,k1,ih,jh
    use modfields, only : thlpcar,rhof,exnf
    use modraddata, only : lwd,lwu,swd,swdir,swdif,swu,thlprad

    implicit none
    integer :: k

    !$acc parallel loop default(present)
    !$omp target teams loop defaultmap(present:aggregate)&
    !$omp defaultmap(present:allocatable)
    do k=1,k1
       lwdav(k)  = 0.
       lwuav(k)  = 0.
       swdav(k)  = 0.
       swdirav(k) = 0.
       swdifav(k) = 0.
       swuav(k)  = 0.
       thltendav(k) = 0.
       thllwtendav(k) = 0.
       thlswtendav(k) = 0.
       thltendav(k) = 0.
       thllwtendcaav(k) = 0.
       thlswtendcaav(k) = 0.
    end do

    call slabsum(lwdav ,1,k1,lwd ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    call slabsum(lwuav ,1,k1,lwu ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    call slabsum(swdav ,1,k1,swd ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    call slabsum(swdirav ,1,k1,swdir ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    call slabsum(swdifav ,1,k1,swdif ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    call slabsum(swuav ,1,k1,swu ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    call slabsum(thltendav ,1,k1,thlprad ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)

    !$acc parallel loop default(present)
    !$omp target teams loop defaultmap(present:aggregate)&
    !$omp defaultmap(present:allocatable)
    do k=1,kmax
       thllwtendav(k) = (abs(lwdav(k+1)) - abs(lwuav(k+1)) - abs(lwdav(k)) + abs(lwuav(k)) )/(rhof(k)*exnf(k)*cp*dzf(k))
       thlswtendav(k) = (abs(swdav(k+1)) - abs(swuav(k+1)) - abs(swdav(k)) + abs(swuav(k)) )/(rhof(k)*exnf(k)*cp*dzf(k))
       ! absolute values here to handle the different sign conventions in different radiation schemes.
    end do

 !    ADD SLAB AVERAGES TO TIME MEAN
    !$acc parallel loop default(present)
    !$omp target teams loop defaultmap(present:aggregate)&
    !$omp defaultmap(present:allocatable)
    do k=1,k1
       lwumn(k)       = lwumn(k)       + lwuav(k)       / ijtot
       lwdmn(k)       = lwdmn(k)       + lwdav(k)       / ijtot
       swdmn(k)       = swdmn(k)       + swdav(k)       / ijtot
       swdirmn(k)     = swdirmn(k)     + swdirav(k)     / ijtot
       swdifmn(k)     = swdifmn(k)     + swdifav(k)     / ijtot
       swumn(k)       = swumn(k)       + swuav(k)       / ijtot
       thltendmn(k)   = thltendmn(k)   + thltendav(k)   / ijtot
       thllwtendmn(k) = thllwtendmn(k) + thllwtendav(k) / ijtot
       thlswtendmn(k) = thlswtendmn(k) + thlswtendav(k) / ijtot
       thlradlsmn(k)  = thlradlsmn(k)  + thlpcar(k)
    end do


    if (lradclearair) then
        call radclearair
        !$acc parallel loop default(present)
        !$omp target teams loop defaultmap(present:aggregate)&
        !$omp defaultmap(present:allocatable)
        do k=1,kmax
          thllwtendcaav(k) = (-lwdcaav(k+1) - lwucaav(k+1) + lwdcaav(k) + lwucaav(k))/(rhof(k)*exnf(k)*cp*dzf(k))
          thlswtendcaav(k) = (-swdcaav(k+1) - swucaav(k+1) + swdcaav(k) + swucaav(k))/(rhof(k)*exnf(k)*cp*dzf(k))
        enddo

        !$acc parallel loop default(present)
        !$omp target teams loop defaultmap(present:aggregate)&
        !$omp defaultmap(present:allocatable)
        do k=1,k1
           thllwtendcamn(k) = thllwtendcamn(k) + thllwtendcaav(k) / ijtot
           thlswtendcamn(k) = thlswtendcamn(k) + thlswtendcaav(k) / ijtot
        end do
    endif


  end subroutine do_radstat

  subroutine radclearair
    use modradfull,    only : d4stream
    use modglobal,    only : i1,ih,j1,jh,kmax,k1,cp,rlv,rd,pref0,ijtot
    use modfields,    only : rhof, exnf, thl0,qt0,ql0
    use modsurfdata,  only : albedo, tskin, qskin, thvs, ps
    use modmicrodata, only : Nc_0
    use modmpi,    only :  slabsum
    use modraddata, only: irad_full, iradiation, swdca,swuca,lwdca,lwuca
    implicit none
    real, dimension(k1)  :: rhof_b, exnf_b
    real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1) :: temp_b, qv_b, ql_b
    integer :: i,j,k

    real :: exnersurf
    !$acc parallel loop default(present)
    !$omp target teams loop defaultmap(present:aggregate)&
    !$omp defaultmap(present:allocatable)
    do k=1,k1
       lwdcaav(k)  = 0.
       lwucaav(k)  = 0.
       swdcaav(k)  = 0.
       swucaav(k)  = 0.
    end do

    !irad_full case is not ported to GPU, see modradiation.f90
    if (iradiation.eq.irad_full) then   !rrtmg has calculated lwdca already
       !take care of UCLALES z-shift for thermo variables.
      do k=1,kmax
        rhof_b(k+1)     = rhof(k)
        exnf_b(k+1)     = exnf(k)
        do j=2,j1
          do i=2,i1
            qv_b(i,j,k+1)   = qt0(i,j,k) - ql0(i,j,k)
            ql_b(i,j,k+1)   = 0.
            temp_b(i,j,k+1) = thl0(i,j,k)*exnf(k)+(rlv/cp)*ql0(i,j,k)
          end do
        end do
      end do

      !take care of the surface boundary conditions
      !CvH edit, extrapolation creates instability in surface scheme
      exnersurf = (ps/pref0) ** (rd/cp)
      rhof_b(1) = ps / (rd * thvs * exnersurf)
      exnf_b(1) = exnersurf

      !rhof_b(1) = rhof(1) + 2*zf(1)/dzf(1)*(rhof(1)-rhof(2))
      !exnf_b(1) = exnh(1) + 0.5*dzf(1)*(exnh(1)-exnf(1))

      do j=2,j1
        do i=2,i1
          ql_b(i,j,1)   = 0.! CvH, no ql at surface
          qv_b(i,j,1)   = qskin(i,j) !CvH, no ql at surface thus qv = qt
          temp_b(i,j,1) = tskin(i,j)*exnersurf
        end do
      end do

      call d4stream(i1,ih,j1,jh,k1,tskin,albedo,real(Nc_0),rhof_b,exnf_b*cp,temp_b,qv_b,ql_b,swdca,swuca,lwdca,lwuca)
    end if

    call slabsum(lwdcaav ,1,k1,lwdca ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    call slabsum(lwucaav ,1,k1,lwuca ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    call slabsum(swdcaav ,1,k1,swdca ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    call slabsum(swucaav ,1,k1,swuca ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)

 !    ADD SLAB AVERAGES TO TIME MEAN

    !$acc parallel loop default(present)
    !$omp target teams loop defaultmap(present:aggregate)&
    !$omp defaultmap(present:allocatable)
    do k=1,k1
       lwucamn(k) = lwucamn(k) + lwucaav(k)/ijtot
       lwdcamn(k) = lwdcamn(k) + lwdcaav(k)/ijtot
       swdcamn(k) = swdcamn(k) + swdcaav(k)/ijtot
       swucamn(k) = swucamn(k) + swucaav(k)/ijtot
    end do
  end subroutine radclearair

!> Write the statistics to file
  subroutine writeradstat
      use modmpi,    only : myid
      use modglobal, only : cexpnr,ifoutput,kmax,k1,zf,zh,rtimee
      use modstat_nc, only: lnetcdf, writestat_nc
      use modgenstat, only: ncid_prof=>ncid,nrec_prof=>nrec
      use modraddata, only : iradiation,irad_par,irad_rrtmg,irad_rte_rrtmgp
      implicit none
      real,dimension(k1,nvar) :: vars
      integer nsecs, nhrs, nminut,k


      nsecs   = nint(rtimee)
      nhrs    = int(nsecs/3600)
      nminut  = int(nsecs/60)-nhrs*60
      nsecs   = mod(nsecs,60)

      !$acc parallel loop default(present)
      !$omp target teams loop defaultmap(present:aggregate)&
      !$omp defaultmap(present:allocatable)
      do k=1,k1
         lwumn(k)   = lwumn(k)    /nsamples
         lwdmn(k)   = lwdmn(k)    /nsamples
         swdmn(k)   = swdmn(k)    /nsamples
         swdirmn(k)   = swdirmn(k)    /nsamples
         swdifmn(k)   = swdifmn(k)    /nsamples
         swumn(k)   = swumn(k)    /nsamples
         lwucamn(k)   = lwucamn(k)    /nsamples
         lwdcamn(k)   = lwdcamn(k)    /nsamples
         swdcamn(k)   = swdcamn(k)    /nsamples
         swucamn(k)   = swucamn(k)    /nsamples
         thllwtendmn(k) = thllwtendmn(k) /nsamples
         thlswtendmn(k) = thlswtendmn(k) /nsamples
         thllwtendcamn(k) = thllwtendcamn(k) /nsamples
         thlswtendcamn(k) = thlswtendcamn(k) /nsamples
         thlradlsmn(k)  = thlradlsmn(k)  /nsamples
         thltendmn(k)   = thltendmn(k)   /nsamples
      end do

      !$acc update self(lwumn, lwdmn, swdmn, swumn, thllwtendmn, thlswtendmn, &
      !$acc&            lwucamn, lwdcamn, swucamn, swdcamn, swdirmn, swdifmn, &
      !$acc&            thltendmn, thlradlsmn, thllwtendcamn, thlswtendcamn)
      !$omp target update from(lwumn,lwdmn,swdmn,swumn,thllwtendmn,&
      !$omp thlswtendmn,lwucamn,lwdcamn,swucamn,swdcamn,swdirmn,swdifmn,&
      !$omp thltendmn,thlradlsmn,thllwtendcamn,thlswtendcamn)

  !     ----------------------
  !     2.0  write the fields
  !           ----------------

    if(myid==0)then
      open (ifoutput,file='radstat.'//cexpnr,position='append')
      write(ifoutput,'(//A,/A,F5.0,A,I4,A,I2,A,I2,A)') &
      '#--------------------------------------------------------'      &
      ,'#',(timeav),'--- AVERAGING TIMESTEP --- '      &
      ,nhrs,':',nminut,':',nsecs      &
      ,'   HRS:MIN:SEC AFTER INITIALIZATION '
      write (ifoutput,'(A/2A/2A)') &
          '#--------------------------------------------------------------------------' &
          ,'#LEV RAD_FLX_HGHT  THL_HGHT  LW_UP        LW_DN        SW_UP       SW_DN       ' &
          ,'TL_LW_TEND   TL_SW_TEND   TL_LS_TEND   TL_TEND' &
          ,'#    (M)    (M)      (W/M^2)      (W/M^2)      (W/M^2)      (W/M^2)      ' &
          ,'(K/H)         (K/H)        (K/H)        (K/H)'
      do k=1,kmax
        write(ifoutput,'(I4,2F10.2,12E13.4)') &
            k,zh(k), zf(k),&
            lwumn(k),&
            lwdmn(k),&
            swumn(k),&
            swdmn(k),&
            thllwtendmn(k)*3600,&
            thlswtendmn(k)*3600,&
            thlradlsmn(k) *3600,&
            thltendmn(k)  *3600,&
            lwucamn(k),&
            lwdcamn(k),&
            swucamn(k),&
            swdcamn(k)
      end do
      close (ifoutput)

     if(iradiation == irad_par .or. iradiation ==irad_rrtmg .or. iradiation == irad_rte_rrtmgp) then ! delta eddington or RRTMG)
      open (ifoutput,file='radsplitstat.'//cexpnr,position='append')
      write(ifoutput,'(//A,/A,F5.0,A,I4,A,I2,A,I2,A)') &
      '#--------------------------------------------------------'      &
      ,'#',(timeav),'--- AVERAGING TIMESTEP --- '      &
      ,nhrs,':',nminut,':',nsecs      &
      ,'   HRS:MIN:SEC AFTER INITIALIZATION '
      write (ifoutput,'(A/2A/2A)') &
          '#--------------------------------------------------------------------------' &
          ,'#LEV  HGHT        LW_UP       LW_DN         SW_UP      SW_DIR_DN    ' &
          ,'SW_DIF_DN     SW_DN      TL_SW_TEND' &
          ,'#  (M)           (W/M^2)     (W/M^2)       (W/M^2)      (W/M^2)     ' &
          ,'(W/M^2)      (W/M^2)       (K/DAY)'
      do k=1,kmax
        write(ifoutput,'(I4,F10.2,7E13.4)') &
            k,zh(k),&
            lwumn(k),&
            lwdmn(k),&
            swumn(k),&
            swdirmn(k),&
            swdifmn(k),&
            swdmn(k),&
            thlswtendmn(k)*3600*24
      end do
      close (ifoutput)
      endif

      if (lnetcdf) then
        vars(:, 1) = thltendmn
        vars(:, 2) = thllwtendmn
        vars(:, 3) = thlswtendmn
        vars(:, 4) = thlradlsmn
        vars(:, 5) = abs(lwumn)
        vars(:, 6) = abs(lwdmn)
        vars(:, 7) = abs(swumn)
        vars(:, 8) = abs(swdmn)
        vars(:, 9) = abs(lwucamn)
        vars(:,10) = abs(lwdcamn)
        vars(:,11) = abs(swucamn)
        vars(:,12) = abs(swdcamn)
        vars(:,13) = thllwtendcamn
        vars(:,14) = thlswtendcamn

       call writestat_nc(ncid_prof,nvar,ncname,vars(1:kmax,:),nrec_prof,kmax)
      end if
    end if ! end if(myid==0)

    !$acc parallel loop default(present)
    !$omp target teams loop defaultmap(present:aggregate)&
    !$omp defaultmap(present:allocatable)
    do k=1,k1
       lwumn(k) = 0.0
       lwdmn(k) = 0.0
       swdmn(k) = 0.0
       swdirmn(k) = 0.0 ! not in netCDF yet
       swdifmn(k) = 0.0 ! not in netCDF yet
       swumn(k) = 0.0
       lwucamn(k) = 0.0
       lwdcamn(k) = 0.0
       swdcamn(k) = 0.0
       swucamn(k) = 0.0
       thllwtendmn(k) = 0.0
       thlswtendmn(k) = 0.0
       thlradlsmn(k)  = 0.0
       thltendmn(k)  = 0.0
       thllwtendcamn(k) = 0.0
       thlswtendcamn(k) = 0.0
    end do

  end subroutine writeradstat

!> Cleans up after the run
  subroutine exitradstat
    implicit none

    !deallocate variables that are needed in modradiation

    if(.not.(lstat)) return

    !$acc exit data delete(lwuav, lwdav, swdav, swdirav, swdifav, swuav, lwucaav, lwdcaav, swdcaav, swucaav, &
    !$acc&                 thllwtendav, thltendav,thlswtendav, thllwtendcaav, thlswtendcaav, &
    !$acc&                 lwumn, lwdmn, swdmn, swdirmn, swdifmn, swumn, lwucamn, lwdcamn, swdcamn, swucamn, &
    !$acc&                 thllwtendmn, thltendmn, thlswtendmn, thlradlsmn, thllwtendcamn, thlswtendcamn)
    !$omp target exit data map(delete:lwuav,lwdav,swdav,swdirav,swdifav,&
    !$omp swuav,lwucaav,lwdcaav,swdcaav,swucaav,thllwtendav,thltendav,&
    !$omp thlswtendav,thllwtendcaav,thlswtendcaav,lwumn,lwdmn,swdmn,&
    !$omp swdirmn,swdifmn,swumn,lwucamn,lwdcamn,swdcamn,swucamn,&
    !$omp thllwtendmn,thltendmn,thlswtendmn,thlradlsmn,thllwtendcamn,&
    !$omp thlswtendcamn)

    deallocate(lwuav,lwdav,swdav,swdirav,swdifav,swuav)
    deallocate(lwucaav, lwdcaav, swucaav, swdcaav)
    deallocate(thllwtendav,thlswtendav, thltendav)
    deallocate(lwumn,lwdmn,swdmn,swdirmn,swdifmn,swumn)
    deallocate(lwucamn, lwdcamn, swucamn, swdcamn)
    deallocate(thllwtendmn,thlswtendmn,thlradlsmn, thltendmn)
    deallocate(thllwtendcaav,thlswtendcaav)
    deallocate(thllwtendcamn,thlswtendcamn)

  end subroutine exitradstat


end module modradstat
