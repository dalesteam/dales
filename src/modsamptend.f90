!> \file modsamptend.f90
!!  Calculates the tendencies of the main fields


!>
!!  Calculates the tendencies of the main fields
!>
!! Profiles of the individual terms of the prognostic equations.  Written to *tend.expnr
!! If netcdf is true, this module also writes in the profiles.expnr.nc output
!!  \author Thijs Heus, MPI
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
module modsamptend
  use modglobal, only : longint
  implicit none
  private
  public :: initsamptend, samptend, exitsamptend, leibniztend
  save
!NetCDF variables
  integer, parameter :: nvar = 66
  integer, dimension(10) :: nrec5
  integer, dimension(10) :: ncid5
  character(80), dimension(10) :: fname5
  character(80),dimension(nvar,4) :: ncname5
  character(80),dimension(1,4) :: tncname5
  real    :: dtav, timeav
  integer(kind=longint) :: idtav,itimeav,tnext,tnextwrite
  integer,public,parameter :: tend_tot=1,tend_start=1,tend_adv=2,tend_subg=3,tend_force=4,&
                       tend_rad=5,tend_ls=6,tend_micro=7, tend_topbound=8,tend_pois=9,tend_addon=10, tend_coriolis=11, tend_totlb=12
  integer,parameter :: nrfields = 12
  character(20),dimension(10) :: samplname,longsamplname
  integer :: nsamples,isamp,isamptot
  logical :: lsampcl  = .false. !< switch for conditional sampling cloud (on/off)
  logical :: lsampco  = .false. !< switch for conditional sampling core (on/off)
  logical :: lsampup  = .false. !< switch for conditional sampling updraft (on/off)
  logical :: lsampbuup  = .false. !< switch for conditional sampling buoyant updraft (on/off)
  logical :: lsampall = .false. !< switch for sampling all data (on/off)
  logical :: ldosamptendwrite = .false. !< write tendencies after leibniz terms have been detemined
  logical :: ldosamptendleib = .false. !< determine leibniz terms
  real :: lastrk3coef

  real, allocatable :: uptm(:,:,:),vptm(:,:,:),wptm(:,:,:),thlptm(:,:,:),qtptm(:,:,:),qrptm(:,:,:),nrptm(:,:,:)
  real, allocatable :: upav(:,:,:),vpav(:,:,:),wpav(:,:,:),thlpav(:,:,:),qtpav(:,:,:),qrpav(:,:,:),nrpav(:,:,:)
  real, allocatable :: upmn(:,:,:),vpmn(:,:,:),wpmn(:,:,:),thlpmn(:,:,:),qtpmn(:,:,:),qrpmn(:,:,:),nrpmn(:,:,:)
  real, allocatable :: ust(:,:),vst(:,:),wst(:,:),thlst(:,:),qtst(:,:),qrst(:,:),nrst(:,:)
  logical, allocatable :: tendmask(:,:,:,:)
  integer, allocatable :: nrsamptot(:,:),nrsamp(:,:),nrsamplast(:,:),nrsampnew(:,:)

contains
!> Initialization routine, reads namelists and inits variables
subroutine initsamptend
    use modmpi,   only : mpierr,my_real,mpi_logical,comm3d,myid,cmyid
    use modglobal,only : cexpnr,dtmax,imax,jmax,kmax,ifnamopt,fname_options,k1,dtav_glob,timeav_glob,ladaptive, dt_lim,btime,kmax,tres,ifoutput,cexpnr,j1,jh,i1,ih
    use modstat_nc, only : open_nc,define_nc,redefine_nc,ncinfo,writestat_dims_nc,lnetcdf

    implicit none
    integer :: ierr,i

    namelist/NAMSAMPLING/ &
    dtav,timeav,lsampcl,lsampco,lsampup,lsampbuup,lsampall

    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMSAMPLING,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions NAMSAMPLING'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMSAMPLING'
      endif
      write(6 ,NAMSAMPLING)
      close(ifnamopt)
    end if

    call MPI_BCAST(dtav,1,MY_REAL,0,comm3d,mpierr)
    call MPI_BCAST(timeav,1,MY_REAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampcl,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampco,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampup,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampall,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampbuup,1,MPI_LOGICAL,0,comm3d,mpierr)

    isamptot = 0
    if (lsampall) then
      isamptot = isamptot + 1
      samplname (isamptot) = 'all'
      longsamplname(isamptot) = 'All '
    endif
    if (lsampup) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'upd'
      longsamplname(isamptot) = 'Updraft '
    end if
    if (lsampbuup) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'buupd'
      longsamplname(isamptot) = 'Buoyant Updraft '
    end if
    if (lsampcl) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'cld'
      longsamplname(isamptot) = 'Cloud '
    end if
    if (lsampco) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'cldcr'
      longsamplname(isamptot) = 'Cloud Core '
    end if

    if(isamptot == 0) return
    if(.not.(lnetcdf)) return !only in netcdf at the moment

    idtav = dtav/tres
    itimeav = timeav/tres
    tnext      = idtav   +btime
    tnextwrite = itimeav +btime
    nsamples = itimeav/idtav
    dt_lim = min(dt_lim,tnext)

    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'dtav should be a integer multiple of dtmax'
    end if
    if (abs(timeav/dtav-nsamples)>1e-4) then
      stop 'timeav should be a integer multiple of dtav'
    end if

    allocate (uptm(k1,nrfields,isamptot),vptm(k1,nrfields,isamptot),wptm(k1,nrfields,isamptot),thlptm(k1,nrfields,isamptot),qtptm(k1,nrfields,isamptot),qrptm(k1,nrfields,isamptot),nrptm(k1,nrfields,isamptot))
    allocate (upmn(k1,nrfields,isamptot),vpmn(k1,nrfields,isamptot),wpmn(k1,nrfields,isamptot),thlpmn(k1,nrfields,isamptot),qtpmn(k1,nrfields,isamptot),qrpmn(k1,nrfields,isamptot),nrpmn(k1,nrfields,isamptot))
    allocate (upav(k1,nrfields,isamptot),vpav(k1,nrfields,isamptot),wpav(k1,nrfields,isamptot),thlpav(k1,nrfields,isamptot),qtpav(k1,nrfields,isamptot),qrpav(k1,nrfields,isamptot),nrpav(k1,nrfields,isamptot))
    allocate (tendmask(2-ih:i1+ih,2-jh:j1+jh,k1,isamptot))
    allocate (nrsamptot(k1,isamptot),nrsamp(k1,isamptot),nrsamplast(k1,isamptot),nrsampnew(k1,isamptot))
    allocate (ust(k1,isamptot),vst(k1,isamptot),wst(k1,isamptot),thlst(k1,isamptot),qtst(k1,isamptot),qrst(k1,isamptot),nrst(k1,isamptot))
    uptm = 0.
    vptm = 0.
    wptm = 0.
    thlptm = 0.
    qtptm = 0.
    qrptm = 0.
    nrptm = 0.
    upmn = 0.
    vpmn = 0.
    wpmn = 0.
    thlpmn = 0.
    qtpmn = 0.
    qrpmn = 0.
    nrpmn = 0.
    upav = 0.
    vpav = 0.
    wpav = 0.
    thlpav = 0.
    qtpav = 0.
    qrpav = 0.
    nrpav = 0.
    ust = 0.
    vst = 0.
    wst = 0.
    thlst = 0.
    qtst = 0.
    qrst = 0.
    nrst = 0.

    tendmask=.false.
    nrsamp=0
    nrsamptot=0
    nrsamplast=0
    nrsampnew=0

    nrec5=0
    ncid5=0
    do i=1,10
    fname5(i)=''
    enddo
    do isamp=1,isamptot
      ncid5(isamp)=isamp
    enddo

    idtav = dtav/tres
    itimeav = timeav/tres
    tnext      = idtav   +btime
    tnextwrite = itimeav +btime
    nsamples = itimeav/idtav

    if (myid==0) then
      do isamp=1,isamptot
        fname5(isamp)='samptend.'//trim(samplname(isamp))//'.'//cexpnr//'.nc'
        call ncinfo(tncname5(1,:),'time','Time','s','time')
        call ncinfo(ncname5( 1,:),'utendadv','U advective tendency','m/s^2','tt')
        call ncinfo(ncname5( 2,:),'utenddif','U diffusive tendency','m/s^2','tt')
        call ncinfo(ncname5( 3,:),'utendfor','U tendency due to other forces','m/s^2','tt')
        call ncinfo(ncname5( 4,:),'utendcor','U coriolis tendency','m/s^2','tt')
        call ncinfo(ncname5( 5,:),'utendls','U large scale tendency','m/s^2','tt')
        call ncinfo(ncname5( 6,:),'utendtop','U top boundary tendency','m/s^2','tt')
        call ncinfo(ncname5( 7,:),'utendpois','U pressure gradient tendency','m/s^2','tt')
        call ncinfo(ncname5( 8,:),'utendaddon','U in addons tendency','m/s^2','tt')
        call ncinfo(ncname5( 9,:),'utendtot','U total tendency','m/s^2','tt')
        call ncinfo(ncname5(10,:),'utendleib','U total tendency with leibniz terms','m/s^2','tt')
        call ncinfo(ncname5(11,:),'vtendadv','V advective tendency','m/s^2','tt')
        call ncinfo(ncname5(12,:),'vtenddif','V diffusive tendency','m/s^2','tt')
        call ncinfo(ncname5(13,:),'vtendfor','V tendency due to other forces','m/s^2','tt')
        call ncinfo(ncname5(14,:),'vtendcor','V coriolis tendency','m/s^2','tt')
        call ncinfo(ncname5(15,:),'vtendls','V large scale tendency','m/s^2','tt')
        call ncinfo(ncname5(16,:),'vtendtop','V top boundary tendency','m/s^2','tt')
        call ncinfo(ncname5(17,:),'vtendpois','V pressure gradient tendency','m/s^2','tt')
        call ncinfo(ncname5(18,:),'vtendaddon','V in addons tendency','m/s^2','tt')
        call ncinfo(ncname5(19,:),'vtendtot','V total tendency','m/s^2','tt')
        call ncinfo(ncname5(20,:),'vtendleib','V total tendency with leibniz terms','m/s^2','tt')
        call ncinfo(ncname5(21,:),'wtendadv','W advective tendency','m/s^2','mt')
        call ncinfo(ncname5(22,:),'wtenddif','W diffusive tendency','m/s^2','mt')
        call ncinfo(ncname5(23,:),'wtendfor','W tendency due to other forces','m/s^2','mt')
        call ncinfo(ncname5(24,:),'wtendcor','W coriolis tendency','m/s^2','mt')
        call ncinfo(ncname5(25,:),'wtendls','W large scale tendency','m/s^2','mt')
        call ncinfo(ncname5(26,:),'wtendtop','W top boundary tendency','m/s^2','mt')
        call ncinfo(ncname5(27,:),'wtendpois','W pressure gradient tendency','m/s^2','mt')
        call ncinfo(ncname5(28,:),'wtendaddon','W in addons tendency','m/s^2','mt')
        call ncinfo(ncname5(29,:),'wtendtot','W total tendency','m/s^2','mt')
        call ncinfo(ncname5(30,:),'wtendleib','W total tendency with leibniz terms','m/s^2','tt')
        call ncinfo(ncname5(31,:),'tltendadv','theta_l advective tendency','K/s','tt')
        call ncinfo(ncname5(32,:),'tltenddif','theta_l diffusive tendency','K/s','tt')
        call ncinfo(ncname5(33,:),'tltendrad','theta_l radiative tendency','K/s','tt')
        call ncinfo(ncname5(34,:),'tltendmicro','theta_l microphysical tendency','K/s','tt')
        call ncinfo(ncname5(35,:),'tltendls','theta_l large scale tendency','K/s','tt')
        call ncinfo(ncname5(36,:),'tltendtop','theta_l  top boundary tendency','K/s','tt')
        call ncinfo(ncname5(37,:),'tltendaddon','theta_l in addons tendency','K/s','tt')
        call ncinfo(ncname5(38,:),'tltendtot','theta_l total tendency','K/s','tt')
        call ncinfo(ncname5(39,:),'tltendleib','theta_l total tendency with leibniz terms','K/s','tt')
        call ncinfo(ncname5(40,:),'qttendadv','total water content advective tendency','kg/kg/s','tt')
        call ncinfo(ncname5(41,:),'qttenddif','total water content diffusive tendency','kg/kg/s','tt')
        call ncinfo(ncname5(42,:),'qttendrad','total water content radiative tendency','kg/kg/s','tt')
        call ncinfo(ncname5(43,:),'qttendmicro','total water content microphysical tendency','kg/kg/s','tt')
        call ncinfo(ncname5(44,:),'qttendls','total water content large scale tendency','kg/kg/s','tt')
        call ncinfo(ncname5(45,:),'qttendtop','total water content  top boundary tendency','kg/kg/s','tt')
        call ncinfo(ncname5(46,:),'qttendaddon','total water content in addons tendency','kg/kg/s','tt')
        call ncinfo(ncname5(47,:),'qttendtot','total water content total tendency','kg/kg/s','tt')
        call ncinfo(ncname5(48,:),'qttendleib','total water content total tendency with leibniz terms','kg/kg/s','tt')
        call ncinfo(ncname5(49,:),'qrtendadv','total water content advective tendency','kg/kg/s','tt')
        call ncinfo(ncname5(50,:),'qrtenddif','total water content diffusive tendency','kg/kg/s','tt')
        call ncinfo(ncname5(51,:),'qrtendrad','total water content radiative tendency','kg/kg/s','tt')
        call ncinfo(ncname5(52,:),'qrtendmicro','total water content microphysical tendency','kg/kg/s','tt')
        call ncinfo(ncname5(53,:),'qrtendls','total water content large scale tendency','kg/kg/s','tt')
        call ncinfo(ncname5(54,:),'qrtendtop','total water content  top boundary tendency','kg/kg/s','tt')
        call ncinfo(ncname5(55,:),'qrtendaddon','total water content in addons tendency','kg/kg/s','tt')
        call ncinfo(ncname5(56,:),'qrtendtot','total water content total tendency','kg/kg/s','tt')
        call ncinfo(ncname5(57,:),'qrtendleib','total water content total tendency with leibniz terms','kg/kg/s','tt')
        call ncinfo(ncname5(58,:),'nrtendadv','RDNC advective tendency','/kg/s','tt')
        call ncinfo(ncname5(59,:),'nrtenddif','RDNC diffusive tendency','/kg/s','tt')
        call ncinfo(ncname5(60,:),'nrtendrad','RDNC radiative tendency','/kg/s','tt')
        call ncinfo(ncname5(61,:),'nrtendmicro','RDNC microphysical tendency','/kg/s','tt')
        call ncinfo(ncname5(62,:),'nrtendls','RDNC large scale tendency','/kg/s','tt')
        call ncinfo(ncname5(63,:),'nrtendtop','RDNC top boundary tendency','/kg/s','tt')
        call ncinfo(ncname5(64,:),'nrtendaddon','RDNC addons tendency','/kg/s','tt')
        call ncinfo(ncname5(65,:),'nrtendtot','RDNC total tendency','/kg/s','tt')
        call ncinfo(ncname5(66,:),'nrtendleib','RDNC total tendency with leibniz terms','/kg/s','tt')
        call open_nc(fname5(isamp),ncid5(isamp),n3=kmax)
        call define_nc(ncid5(isamp),1,tncname5)
        call writestat_dims_nc(ncid5(isamp))
        call redefine_nc(ncid5(isamp))
        call define_nc(ncid5(isamp),nvar,ncname5)
      enddo
    end if

  end subroutine initsamptend

!> Performs the statistics, keeps track of what the tendencies were last time, and what they are this time.
  subroutine samptend(tendterm,firstterm,lastterm)
    use modmpi,    only : myid,slabsum
    use modglobal, only : i1,i2,j1,j2,kmax,k1,ih,jh,&
                          cp,rv,rlv,rd,rslabs,&
                          grav,om22,cu,timee,rk3step,dt_lim,rslabs,btime,nsv,rdt
    use modfields, only : up,vp,wp,thlp,qtp,svp,w0,thl0,ql0,exnf,qt0,u0,v0,sv0
    use modmicrodata, only : iqr,inr
    use modstat_nc, only : lnetcdf
    implicit none
    integer, intent(in)           :: tendterm !< name of the term to write down
    logical, intent(in), optional :: lastterm !< true if this is the last term of the equations; the write routine is entered.
    logical, intent(in), optional :: firstterm !< true if this is the first term of the equations
    real, allocatable, dimension(:,:,:) :: w0f,wpf
    real, allocatable, dimension(:,:,:) :: thv0
    real, allocatable, dimension(:) :: thvav
    integer :: i,j,k

    if(isamptot == 0) return
    if(.not.(lnetcdf)) return !only in netcdf at the moment
    if (rk3step/=3) return
    if(timee<tnext) then
      dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))
      return
    end if

    IF (present(firstterm)) THEN
    IF (firstterm) THEN
      nrsamplast=0
      tendmask=.false.
      uptm = 0.
      vptm = 0.
      wptm = 0.
      thlptm = 0.
      qtptm = 0.
      qrptm = 0.
      nrptm = 0.
      ust = 0.
      vst = 0.
      wst = 0.
      thlst = 0.
      qtst = 0.
      qrst = 0.
      nrst = 0.

      allocate(thv0(2-ih:i1+ih,2-jh:j1+jh,k1),&
                w0f(2-ih:i1+ih,2-jh:j1+jh,k1))   
      allocate(thvav(k1))
  
      do k=1,k1
        thv0(2:i1,2:j1,k) = (thl0(2:i1,2:j1,k)+rlv*ql0(2:i1,2:j1,k)/(cp*exnf(k))) &
                    *(1+(rv/rd-1)*qt0(2:i1,2:j1,k)-rv/rd*ql0(2:i1,2:j1,k))
      enddo
      do k=1,kmax
        w0f (2:i1,2:j1,k) = 0.5*(w0 (2:i1,2:j1,k) + w0  (2:i1,2:j1,k+1))
      end do
  
      thvav = 0.0
      call slabsum(thvav,1,k1,thv0,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      thvav = thvav/rslabs

      do isamp=1,isamptot
        select case (samplname(isamp))
        case ('upd')
          do i=2,i1
          do j=2,j1
          do k=1,kmax
            if (w0f(i,j,k)>0.) then
                tendmask(i,j,k,isamp) = .true.
            endif
          enddo
          enddo
          enddo
        case ('buup')
          do i=2,i1
          do j=2,j1
          do k=1,kmax
            if ((w0f(i,j,k)>0.0).and.(thv0(i,j,k) > thvav(k))) then
                tendmask(i,j,k,isamp) = .true.
            endif
          enddo
          enddo
          enddo
        case ('cld')
          do i=2,i1
          do j=2,j1
          do k=1,kmax
            if (ql0(i,j,k)>epsilon(1.0)) then
                tendmask(i,j,k,isamp) = .true.
            endif
          enddo
          enddo
          enddo
        case ('cldcr')
          do i=2,i1
          do j=2,j1
          do k=1,kmax
            if (ql0(i,j,k)>epsilon(1.0).and.thv0(i,j,k) > thvav(k)) then
                tendmask(i,j,k,isamp) = .true.
            endif
          enddo
          enddo
          enddo
        case ('all')
            tendmask(:,:,:,isamp)  = .true.
        end select
        do k=1,kmax
          nrsamp(k,isamp)= nrsamp(k,isamp)+count(tendmask(2:i1,2:j1,k,isamp))
          nrsamplast(k,isamp)= count(tendmask(2:i1,2:j1,k,isamp))
        end do
      enddo

      deallocate(thv0,w0f)
      deallocate(thvav)

      do isamp=1,isamptot
      do k=1,kmax
        ust(k,isamp) = sum(u0(2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))
        vst(k,isamp) = sum(v0(2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))
        wst(k,isamp) = sum(w0(2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))
        thlst(k,isamp) = sum(thl0(2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))
        qtst(k,isamp) = sum(qt0(2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))
        if(nsv>1) then
        qrst(k,isamp) = sum(sv0(2:i1,2:j1,k,iqr),tendmask(2:i1,2:j1,k,isamp))
        nrst(k,isamp) = sum(sv0(2:i1,2:j1,k,inr),tendmask(2:i1,2:j1,k,isamp))
        endif
      end do
      end do

      ldosamptendleib=.true.
      lastrk3coef = rdt / (4. - dble(rk3step))

    ENDIF
    ENDIF

    allocate(wpf(2-ih:i1+ih,2-jh:j1+jh,k1))  
    
    do k=1,kmax
      wpf (2:i1,2:j1,k) = 0.5*(wp (2:i1,2:j1,k) + wp  (2:i1,2:j1,k+1))
    end do

    do isamp=1,isamptot
    do k=1,kmax
      uptm(k,tendterm,isamp) = sum(up (2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))-uptm (k,tend_tot,isamp)
      vptm(k,tendterm,isamp) = sum(vp (2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))-vptm (k,tend_tot,isamp)
      wptm(k,tendterm,isamp) = sum(wpf (2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))-wptm (k,tend_tot,isamp)
      thlptm(k,tendterm,isamp) = sum(thlp (2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))-thlptm (k,tend_tot,isamp)
      qtptm(k,tendterm,isamp) = sum(qtp (2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))-qtptm (k,tend_tot,isamp)
      if(nsv>1) then
      qrptm(k,tendterm,isamp) = sum(svp (2:i1,2:j1,k,iqr),tendmask(2:i1,2:j1,k,isamp))-qrptm (k,tend_tot,isamp)
      nrptm(k,tendterm,isamp) = sum(svp (2:i1,2:j1,k,inr),tendmask(2:i1,2:j1,k,isamp))-nrptm (k,tend_tot,isamp)
      endif
      uptm(k,tend_tot,isamp) = sum(up (2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))
      vptm(k,tend_tot,isamp) = sum(vp (2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))
      wptm(k,tend_tot,isamp) = sum(wpf (2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))
      thlptm(k,tend_tot,isamp) = sum(thlp (2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))
      qtptm(k,tend_tot,isamp) = sum(qtp (2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))
      qrptm(k,tend_tot,isamp) = sum(svp (2:i1,2:j1,k,iqr),tendmask(2:i1,2:j1,k,isamp))
      nrptm(k,tend_tot,isamp) = sum(svp (2:i1,2:j1,k,inr),tendmask(2:i1,2:j1,k,isamp))
      upav(k,tendterm,isamp) = upav(k,tendterm,isamp)+uptm(k,tendterm,isamp)
      vpav(k,tendterm,isamp) = vpav(k,tendterm,isamp)+vptm(k,tendterm,isamp)
      wpav(k,tendterm,isamp) = wpav(k,tendterm,isamp)+wptm(k,tendterm,isamp)
      thlpav(k,tendterm,isamp) = thlpav(k,tendterm,isamp)+thlptm(k,tendterm,isamp)
      qtpav(k,tendterm,isamp) = qtpav(k,tendterm,isamp)+qtptm(k,tendterm,isamp)
      qrpav(k,tendterm,isamp) = qrpav(k,tendterm,isamp)+qrptm(k,tendterm,isamp)
      nrpav(k,tendterm,isamp) = nrpav(k,tendterm,isamp)+nrptm(k,tendterm,isamp)
    end do
    end do

    deallocate(wpf)

    IF (present(lastterm)) THEN
    IF (lastterm) THEN
      do isamp=1,isamptot
      do k=1,kmax
        upav(k,tend_tot,isamp) = upav(k,tend_tot,isamp)+uptm(k,tend_tot,isamp)
        vpav(k,tend_tot,isamp) = vpav(k,tend_tot,isamp)+vptm(k,tend_tot,isamp)
        wpav(k,tend_tot,isamp) = wpav(k,tend_tot,isamp)+wptm(k,tend_tot,isamp)
        thlpav(k,tend_tot,isamp) = thlpav(k,tend_tot,isamp)+thlptm(k,tend_tot,isamp)
        qtpav(k,tend_tot,isamp) = qtpav(k,tend_tot,isamp)+qtptm(k,tend_tot,isamp)
        qrpav(k,tend_tot,isamp) = qrpav(k,tend_tot,isamp)+qrptm(k,tend_tot,isamp)
        nrpav(k,tend_tot,isamp) = nrpav(k,tend_tot,isamp)+nrptm(k,tend_tot,isamp)
      enddo
      enddo
      tnext = tnext+idtav

      if (timee>=tnextwrite) then
        tnextwrite = tnextwrite+itimeav
        ldosamptendwrite = .true.
      end if
      dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))
    END IF
    END IF
    
  end subroutine samptend

  subroutine leibniztend
    use modmpi,    only : myid,slabsum
    use modglobal, only : i1,i2,j1,j2,kmax,k1,ih,jh,&
                          cp,rv,rlv,rd,rslabs,&
                          grav,om22,cu,timee,rk3step,dt_lim,rslabs,btime,nsv,rdt
    use modfields, only : up,vp,wp,thlp,qtp,svp,w0,thl0,ql0,exnf,qt0,u0,v0,sv0
    use modmicrodata, only : iqr,inr
    implicit none
    real, allocatable, dimension(:,:,:) :: w0f,wpf
    real, allocatable, dimension(:,:,:) :: thv0
    real, allocatable, dimension(:) :: thvav
    integer :: i,j,k

      if(.not.(ldosamptendleib)) return
      ldosamptendleib=.false.
      tendmask=.false.
      nrsampnew=0

      allocate(thv0(2-ih:i1+ih,2-jh:j1+jh,k1),&
                w0f(2-ih:i1+ih,2-jh:j1+jh,k1))   
      allocate(thvav(k1))
  
      do k=1,k1
        thv0(2:i1,2:j1,k) = (thl0(2:i1,2:j1,k)+rlv*ql0(2:i1,2:j1,k)/(cp*exnf(k))) &
                    *(1+(rv/rd-1)*qt0(2:i1,2:j1,k)-rv/rd*ql0(2:i1,2:j1,k))
      enddo
      do k=1,kmax
        w0f (2:i1,2:j1,k) = 0.5*(w0 (2:i1,2:j1,k) + w0  (2:i1,2:j1,k+1))
      end do
  
      thvav = 0.0
      call slabsum(thvav,1,k1,thv0,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      thvav = thvav/rslabs

      do isamp=1,isamptot
        select case (samplname(isamp))
        case ('upd')
          do i=2,i1
          do j=2,j1
          do k=1,kmax
            if (w0f(i,j,k)>0.) then
                tendmask(i,j,k,isamp) = .true.
            endif
          enddo
          enddo
          enddo
        case ('buup')
          do i=2,i1
          do j=2,j1
          do k=1,kmax
            if ((w0f(i,j,k)>0.0).and.(thv0(i,j,k) > thvav(k))) then
                tendmask(i,j,k,isamp) = .true.
            endif
          enddo
          enddo
          enddo
        case ('cld')
          do i=2,i1
          do j=2,j1
          do k=1,kmax
            if (ql0(i,j,k)>epsilon(1.0)) then
                tendmask(i,j,k,isamp) = .true.
            endif
          enddo
          enddo
          enddo
        case ('cldcr')
          do i=2,i1
          do j=2,j1
          do k=1,kmax
            if (ql0(i,j,k)>epsilon(1.0).and.thv0(i,j,k) > thvav(k)) then
                tendmask(i,j,k,isamp) = .true.
            endif
          enddo
          enddo
          enddo
        case ('all')
            tendmask(:,:,:,isamp)  = .true.
        end select
        do k=1,kmax
          nrsampnew(k,isamp)= count(tendmask(2:i1,2:j1,k,isamp))
        end do
      enddo

      deallocate(thv0,w0f)
      deallocate(thvav)

      do isamp=1,isamptot
      do k=1,kmax
        if((nrsampnew(k,isamp)>0).and.(nrsamplast(k,isamp)>0)) then! only do if sampling can be performed at both points in time
        upav(k,tend_totlb,isamp) = upav(k,tend_totlb,isamp)+(ust(k,isamp)-sum(u0(2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))*nrsamplast(k,isamp)/nrsampnew(k,isamp))/lastrk3coef
        vpav(k,tend_totlb,isamp) = vpav(k,tend_totlb,isamp)+(vst(k,isamp)-sum(v0(2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))*nrsamplast(k,isamp)/nrsampnew(k,isamp))/lastrk3coef
        wpav(k,tend_totlb,isamp) = wpav(k,tend_totlb,isamp)+(wst(k,isamp)-sum(w0(2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))*nrsamplast(k,isamp)/nrsampnew(k,isamp))/lastrk3coef
        thlpav(k,tend_totlb,isamp) = thlpav(k,tend_totlb,isamp)+(thlst(k,isamp)-sum(thl0(2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))*nrsamplast(k,isamp)/nrsampnew(k,isamp))/lastrk3coef
        qtpav(k,tend_totlb,isamp) = qtpav(k,tend_totlb,isamp)+(qtst(k,isamp)-sum(qt0(2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))*nrsamplast(k,isamp)/nrsampnew(k,isamp))/lastrk3coef
        if(nsv>1) then
        qrpav(k,tend_totlb,isamp) = qrpav(k,tend_totlb,isamp)+(qrst(k,isamp)-sum(sv0(2:i1,2:j1,k,iqr),tendmask(2:i1,2:j1,k,isamp))*nrsamplast(k,isamp)/nrsampnew(k,isamp))/lastrk3coef
        nrpav(k,tend_totlb,isamp) = nrpav(k,tend_totlb,isamp)+(nrst(k,isamp)-sum(sv0(2:i1,2:j1,k,inr),tendmask(2:i1,2:j1,k,isamp))*nrsamplast(k,isamp)/nrsampnew(k,isamp))/lastrk3coef
        endif
        endif
      enddo
      enddo

      if(ldosamptendwrite) then
        ldosamptendwrite=.false.
        call writesamptend
        upav = 0.
        vpav = 0.
        wpav = 0.
        thlpav = 0.
        qtpav = 0.
        qrpav = 0.
        nrpav = 0.
        nrsamp = 0
      endif

  end subroutine leibniztend

!> Write the statistics to file
  subroutine writesamptend
    use modglobal, only : kmax,k1,rtimee
    use modmpi,    only : mpi_allreduce,mpi_integer,myid,comm3d,mpierr,my_real,mpi_sum
    use modstat_nc, only: lnetcdf,writestat_nc
    implicit none
    integer :: field,k
    real, allocatable :: vars(:,:)
    allocate(vars(1:k1,nvar))
    upmn = 0.
    vpmn = 0.
    wpmn = 0.
    thlpmn = 0.
    qtpmn = 0.
    qrpmn = 0.
    nrpmn = 0.
    nrsamptot=0

    call MPI_ALLREDUCE(nrsamp   ,nrsamptot ,k1*isamptot,MPI_INTEGER,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(upav     ,upmn      ,k1*nrfields*isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(vpav     ,vpmn      ,k1*nrfields*isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(wpav     ,wpmn      ,k1*nrfields*isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(thlpav   ,thlpmn    ,k1*nrfields*isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(qtpav    ,qtpmn     ,k1*nrfields*isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(qrpav    ,qrpmn     ,k1*nrfields*isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(nrpav    ,nrpmn     ,k1*nrfields*isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)

    do field=1,nrfields
    do isamp=1,isamptot
    do k=1,k1
      if (nrsamptot(k,isamp)>0) then
        upmn  (k,field,isamp) = upmn (k,field,isamp)/nrsamptot(k,isamp)
        vpmn  (k,field,isamp) = vpmn (k,field,isamp)/nrsamptot(k,isamp)
        wpmn  (k,field,isamp) = wpmn (k,field,isamp)/nrsamptot(k,isamp)
        thlpmn(k,field,isamp) = thlpmn (k,field,isamp)/nrsamptot(k,isamp)
        qtpmn (k,field,isamp) = qtpmn (k,field,isamp)/nrsamptot(k,isamp)
        qrpmn (k,field,isamp) = qrpmn (k,field,isamp)/nrsamptot(k,isamp)
        nrpmn (k,field,isamp) = nrpmn (k,field,isamp)/nrsamptot(k,isamp)
      endif
    enddo
    enddo
    enddo

    if(myid == 0) then 
      if (lnetcdf) then
        do isamp=1,isamptot
          vars=0.
          vars(:, 1) = upmn(:,tend_adv,isamp)
          vars(:, 2) = upmn(:,tend_subg,isamp)
          vars(:, 3) = upmn(:,tend_force,isamp)
          vars(:, 4) = upmn(:,tend_coriolis,isamp)
          vars(:, 5) = upmn(:,tend_ls,isamp)
          vars(:, 6) = upmn(:,tend_topbound,isamp)
          vars(:, 7) = upmn(:,tend_pois,isamp)
          vars(:, 8) = upmn(:,tend_addon,isamp)
          vars(:, 9) = upmn(:,tend_tot,isamp)
          vars(:,10) = upmn(:,tend_totlb,isamp)
          vars(:,11) = vpmn(:,tend_adv,isamp)
          vars(:,12) = vpmn(:,tend_subg,isamp)
          vars(:,13) = vpmn(:,tend_force,isamp)
          vars(:,14) = vpmn(:,tend_coriolis,isamp)
          vars(:,15) = vpmn(:,tend_ls,isamp)
          vars(:,16) = vpmn(:,tend_topbound,isamp)
          vars(:,17) = vpmn(:,tend_pois,isamp)
          vars(:,18) = vpmn(:,tend_addon,isamp)
          vars(:,19) = vpmn(:,tend_tot,isamp)
          vars(:,20) = vpmn(:,tend_totlb,isamp)
          vars(:,21) = wpmn(:,tend_adv,isamp)
          vars(:,22) = wpmn(:,tend_subg,isamp)
          vars(:,23) = wpmn(:,tend_force,isamp)
          vars(:,24) = wpmn(:,tend_coriolis,isamp)
          vars(:,25) = wpmn(:,tend_ls,isamp)
          vars(:,26) = wpmn(:,tend_topbound,isamp)
          vars(:,27) = wpmn(:,tend_pois,isamp)
          vars(:,28) = wpmn(:,tend_addon,isamp)
          vars(:,29) = wpmn(:,tend_tot,isamp)
          vars(:,30) = wpmn(:,tend_totlb,isamp)
          vars(:,31) = thlpmn(:,tend_adv,isamp)
          vars(:,32) = thlpmn(:,tend_subg,isamp)
          vars(:,33) = thlpmn(:,tend_rad,isamp)
          vars(:,34) = thlpmn(:,tend_micro,isamp)
          vars(:,35) = thlpmn(:,tend_ls,isamp)
          vars(:,36) = thlpmn(:,tend_topbound,isamp)
          vars(:,37) = thlpmn(:,tend_addon,isamp)
          vars(:,38) = thlpmn(:,tend_tot,isamp)
          vars(:,39) = thlpmn(:,tend_totlb,isamp)
          vars(:,40) = qtpmn(:,tend_adv,isamp)
          vars(:,41) = qtpmn(:,tend_subg,isamp)
          vars(:,42) = qtpmn(:,tend_rad,isamp)
          vars(:,43) = qtpmn(:,tend_micro,isamp)
          vars(:,44) = qtpmn(:,tend_ls,isamp)
          vars(:,45) = qtpmn(:,tend_topbound,isamp)
          vars(:,46) = qtpmn(:,tend_addon,isamp)
          vars(:,47) = qrpmn(:,tend_tot,isamp)
          vars(:,48) = qrpmn(:,tend_totlb,isamp)
          vars(:,49) = qrpmn(:,tend_adv,isamp)
          vars(:,50) = qrpmn(:,tend_subg,isamp)
          vars(:,51) = qrpmn(:,tend_rad,isamp)
          vars(:,52) = qrpmn(:,tend_micro,isamp)
          vars(:,53) = qrpmn(:,tend_ls,isamp)
          vars(:,54) = qrpmn(:,tend_topbound,isamp)
          vars(:,55) = qrpmn(:,tend_addon,isamp)
          vars(:,56) = qrpmn(:,tend_tot,isamp)
          vars(:,57) = qrpmn(:,tend_totlb,isamp)
          vars(:,58) = nrpmn(:,tend_adv,isamp)
          vars(:,59) = nrpmn(:,tend_subg,isamp)
          vars(:,60) = nrpmn(:,tend_rad,isamp)
          vars(:,61) = nrpmn(:,tend_micro,isamp)
          vars(:,62) = nrpmn(:,tend_ls,isamp)
          vars(:,63) = nrpmn(:,tend_topbound,isamp)
          vars(:,64) = nrpmn(:,tend_addon,isamp)
          vars(:,65) = nrpmn(:,tend_tot,isamp)
          vars(:,66) = nrpmn(:,tend_totlb,isamp)
          call writestat_nc(ncid5(isamp),1,tncname5,(/rtimee/),nrec5(isamp),.true.)
          call writestat_nc(ncid5(isamp),nvar,ncname5,vars(1:kmax,:),nrec5(isamp),kmax)
        enddo
      end if
    end if
    deallocate(vars)

  end subroutine writesamptend

!> Cleans up after the run
  subroutine exitsamptend
    use modstat_nc, only: lnetcdf
  implicit none

    if(isamptot == 0) return
    if(.not.(lnetcdf)) return
    deallocate (uptm,vptm,wptm,thlptm,qtptm,qrptm,nrptm)
    deallocate (upmn,vpmn,wpmn,thlpmn,qtpmn,qrpmn,nrpmn)
    deallocate (upav,vpav,wpav,thlpav,qtpav,qrpav,nrpav)
    deallocate (ust,vst,wst,thlst,qtst,qrst,nrst)
    deallocate (tendmask)
    deallocate (nrsamptot,nrsamp,nrsamplast,nrsampnew)

  end subroutine exitsamptend

end module
