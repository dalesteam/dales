!> \file modgenstat.f90
!!  Genstat calculates slab averages of several variables

!>
!!  Genstat calculates slab averages of several variables
!>
!!  Written to fields.expnr, moments.expnr and flux1.expnr and flux2.expnr
!! If netcdf is true, this module leads the profiles.expnr.nc output
!!  \author Hans Cuijpers, K.N.M.I.
!!  \author Pier Siebesma, K.N.M.I.
!!  \author Stephan de Roode,TU Delft
!!  \author Chiel van Heerwaarden, Wageningen U.R.
!!  \author Thijs Heus,MPI-M
!!  \par Revision list
!!  \todo Documentation
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



module modgenstat

    !-----------------------------------------------------------------|
    !                                                                 |
    !*** *stattend*  calculates generic slabaveraged statistics       |
    !                                                                 |
    !      Pier Siebesma   K.N.M.I.     12/05/1995                    |
    !      Hans Cuijpers   I.M.A.U.     21/06/1995                    |
    !                                                                 |
    !     purpose.                                                    |
    !     --------                                                    |
    !                                                                 |
    !     stattend.f calculates:                                      |
    !                                                                 |
    !     * time averaged fieldsof theta_l, theta, theta_v, qt, qv   |
    !        ql, u and v                                              |
    !     * time averaged tendencies  of theta_l, theta, qt, qv, ql,  |
    !        u and v                                                  |
    !     * time averaged turbulent fluxes of theta_l, theta_v,       |
    !       theta, qt, qv, ql, u and v                                |
    !     * variances of qt, w, u, theta_l and theta_v                |
    !     * skewness of qt and w                                      |
    !*** *genstat*  calculates timeseries of several variables       |
    !                                                                 |
    !____________________SETTINGS_AND_SWITCHES________________________|
    !                     IN &NAMTIMESTAT                             |
    !                                                                 |
    !    dtav           SAMPLING INTERVAL                             |
    !                                                                 |
    !    timeav         INTERVAL OF WRITING                           |
    !                                                                 |
    !    lstat      SWITCH TO ENABLE TIMESERIES                       |
    !-----------------------------------------------------------------|
  use modprecision

  implicit none
  ! private
  PUBLIC :: initgenstat, genstat, exitgenstat
  save

!NetCDF variables
  integer :: nvar = 40
  integer :: ncid,nrec = 0
  character(80) :: fname = 'profiles.xxx.nc'
  character(80),allocatable, dimension(:,:) :: ncname
  character(80),dimension(1,4) :: tncname

  real    :: dtav, timeav
  integer(kind=longint) :: idtav,itimeav,tnext,tnextwrite
  logical :: lstat= .false. ! switch for conditional sampling cloud (on/off)
  integer :: nsamples
!     ----  total fields  ---

  real, allocatable  :: umn   (:)       ,vmn   (:)
  real, allocatable  :: thlmn (:)       ,thvmn (:)
  real, allocatable  :: qtmn  (:)       ,qlmn  (:),  qlhmn(:),cfracmn(:)

! real, allocatable  ::     --- fluxes (resolved, subgrid and total) ---
  real, allocatable  :: wthlsmn (:),wthlrmn (:),wthltmn(:)
  real, allocatable  :: wthvsmn (:),wthvrmn (:),wthvtmn(:)
  real, allocatable  :: wqlsmn (:),wqlrmn (:),wqltmn(:)
  real, allocatable  :: wqtsmn (:),wqtrmn (:),wqttmn(:)
  real, allocatable  :: wsvsmn (:,:),wsvrmn(:,:),wsvtmn(:,:)
  real, allocatable  :: uwtmn  (:),vwtmn  (:) !total    uw, vw
  real, allocatable  :: uwsmn  (:),vwsmn  (:) !resolved uw, vw
  real, allocatable  :: uwrmn  (:),vwrmn  (:) !subgrid  uw, vw
! real, allocatable  ::     --- various moments ---

!   real, allocatable  :: rmn        (:), r2mn   (:), r3mn (:), rhmn (:)
  real, allocatable  :: w2mn       (:), skewmn (:)
  real, allocatable  :: w2submn    (:)
  real, allocatable  :: u2mn       (:), v2mn  (:),     qt2mn(:)
  real, allocatable  :: thl2mn     (:), thv2mn(:),     th2mn(:),     ql2mn(:)
!   real, allocatable  :: qs2mn      (:), qsmn  (:)
  real, allocatable  :: svmmn(:,:),svptmn(:,:),svplsmn(:,:),svpmn(:,:)
  real, allocatable  :: sv2mn(:,:)

  real(field_r), allocatable :: umav (:)     ! slab averaged ql_0    at full level
  real(field_r), allocatable :: vmav (:)     ! slab averaged ql_0    at full level
  real(field_r), allocatable :: thlmav (:)     ! slab averaged ql_0    at full level
  real(field_r), allocatable :: thmav (:)     ! slab averaged ql_0    at full level
  real(field_r), allocatable :: qtmav (:)     ! slab averaged ql_0    at full level
  real(field_r), allocatable :: qlmav (:)     ! slab averaged ql_0    at full level
  real, allocatable :: cfracav (:)     ! slab averaged ql_0    at full level
  real(field_r), allocatable :: svmav (:,:)     ! slab averaged ql_0    at full level
  real, allocatable :: svpav(:,:)                  !  slab average total tendency of sv(n)
  real, allocatable :: svptav(:,:)                 !  slab average tendency of sv(n) due to turb.

  real, allocatable :: uptav(:)                      !  slab averaged tendency for u
  real, allocatable :: vptav(:)                      !  slab averaged tendency for v

  real, allocatable :: thptav(:)     ! slab averaged turbulence tendency of theta
  real, allocatable :: qlptav(:)     ! slab averaged turbulence tendency of q_liq
  real, allocatable :: uwtot (:)     ! slab averaged tot w-u flux at half levels
  real, allocatable :: vwtot (:)     ! slab averaged tot w-v flux at half levels
  real, allocatable :: uwsub (:)     ! slab averaged sgs w-u flux at half levels
  real, allocatable :: vwsub (:)     ! slab averaged sgs w-v flux at half levels
  real, allocatable :: uwres (:)     ! slab averaged res w-u flux at half levels
  real, allocatable :: vwres (:)     ! slab averaged res w-v flux at half levels


  real, allocatable :: thvhav(:)
  real, allocatable :: th0av(:)
  real, allocatable :: wthlsub(:)     ! slab averaged sub w-theta_l flux at half levels
  real, allocatable :: wthlres(:)     ! slab averaged res w-theta_l flux at half levels
  real, allocatable :: wthltot(:)     ! slab averaged tot w-theta_l flux at half levels

  real, allocatable :: wqtsub(:)    ! slab averaged sub w-qtot    flux at half levels
  real, allocatable :: wqtres(:)    ! slab averaged res w-qtot    flux at half levels
  real, allocatable :: wqttot(:)    ! slab averaged tot w-qtot    flux at half levels

  real, allocatable :: wqlsub(:)    ! slab averaged sub w-ql      flux at half levels
  real, allocatable :: wqlres(:)    ! slab averaged tot w-ql      flux at half levels
  real, allocatable :: wqltot(:)    ! slab averaged tot w-ql      flux at half levels

  real, allocatable :: wthvtot(:)    !  slab averaged total wthv-flux
  real, allocatable :: wthvres(:)    !  slab averaged res   wthv-flux
  real, allocatable :: wthvsub(:)    !  slab averaged sub   wthv-flux

  real, allocatable :: wsvsub(:,:)! slab averaged sub w-sv(n)  flux
  real, allocatable :: wsvres(:,:)! slab averaged res w-sv(n)  flux
  real, allocatable :: wsvtot(:,:)! slab averaged tot w-sv(n)  flux

  real, allocatable :: cszmn(:)    ! Smagorinsky constant
  real, allocatable :: cszav(:)    ! Smagorinsky constant

  real, allocatable :: qlmnlast(:)
  real, allocatable :: wthvtmnlast(:)

    ! Local fields

  real,allocatable, dimension(:):: qlptavl   ! slab averaged turbulence tendency of q_liq
  real,allocatable, dimension(:):: uwsubl
  real,allocatable, dimension(:):: vwsubl
  real,allocatable, dimension(:):: uwresl
  real,allocatable, dimension(:):: vwresl
  real,allocatable, dimension(:):: qlhav
  real,allocatable, dimension(:):: u2av    , &
            v2av    , &
            w2av    , &
            w3av    , &
            w2subav , &
            qt2av   , &
            thl2av  , &
            thv2av  , &
            th2av   , &
            ql2av
  real, allocatable, dimension(:,:) :: sv2av
  real(field_r),allocatable, dimension(:,:,:)::  thv0
  real(field_r),allocatable, dimension(:)::   thvmav
  real(field_r),allocatable, dimension(:,:,:):: sv0h

contains

  subroutine initgenstat
    use modmpi, only : myid,mpierr, comm3d, mpi_logical, D_MPI_BCAST
    use modglobal, only : i1, ih, j1, jh, kmax, k1, nsv, ifnamopt, fname_options, ifoutput, &
                          cexpnr, dtav_glob, timeav_glob, dt_lim, btime, tres, &
                          lwarmstart, checknamelisterror
    use modstat_nc, only : lnetcdf, open_nc, ncinfo, define_nc, nctiminfo, writestat_dims_nc
    use modsurfdata, only : isurf, ksoilmax

    implicit none

    integer n, ierr
    character(40) :: name
    character(3) :: csvname
    namelist/NAMGENSTAT/ &
    dtav,timeav,lstat

    dtav=dtav_glob;timeav=timeav_glob

    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMGENSTAT,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMGENSTAT')
      write(6 ,NAMGENSTAT)
      close(ifnamopt)
    end if

    call D_MPI_BCAST(timeav     ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(dtav       ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(lstat      ,1,0,comm3d,mpierr)
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

    allocate(umn(k1)       ,vmn   (k1))
    allocate(thlmn (k1)       ,thvmn (k1))
    allocate(qtmn  (k1)       ,qlmn  (k1),  qlhmn(k1),cfracmn(k1))
    allocate(wthlsmn (k1),wthlrmn (k1),wthltmn(k1))
    allocate(wthvsmn (k1),wthvrmn (k1),wthvtmn(k1))
    allocate(wqlsmn (k1),wqlrmn (k1),wqltmn(k1))
    allocate(wqtsmn (k1),wqtrmn (k1),wqttmn(k1))
    allocate(wsvsmn (k1,nsv),wsvrmn(k1,nsv),wsvtmn(k1,nsv))
    allocate(uwtmn  (k1),vwtmn  (k1),uwrmn  (k1),  vwrmn(k1),uwsmn(k1),vwsmn(k1))
!     allocate(rmn(k1), r2mn   (k1), r3mn (k1), rhmn (k1))
    allocate(w2mn       (k1), skewmn (k1))
    allocate(w2submn    (k1))
    allocate(u2mn       (k1), v2mn  (k1),     qt2mn(k1))
    allocate(thl2mn     (k1), thv2mn(k1),     th2mn(k1),     ql2mn(k1))
!     allocate(qs2mn      (k1), qsmn  (k1))
    allocate(svmmn(k1,nsv),svptmn(k1,nsv),svplsmn(k1,nsv),svpmn(k1,nsv))
    allocate(sv2mn(k1,nsv))

    allocate(umav (k1))
    allocate(vmav (k1))
    allocate(thlmav (k1))
    allocate(thmav (k1))
    allocate(qtmav (k1))
    allocate(qlmav (k1))
    allocate(cfracav(k1))
    allocate(svmav (k1,nsv))
    allocate(uptav(k1))
    allocate(vptav(k1))

    allocate(thptav(k1))
    allocate(qlptav(k1))
    allocate(uwtot (k1))
    allocate(vwtot (k1))
    allocate(uwsub (k1))
    allocate(vwsub (k1))
    allocate(uwres (k1))
    allocate(vwres (k1))
    allocate(wthlsub(k1))
    allocate(wthlres(k1))
    allocate(wthltot(k1))
    allocate(wqtsub(k1))
    allocate(wqtres(k1))
    allocate(wqttot(k1))
    allocate(wqlsub(k1))
    allocate(wqlres(k1))
    allocate(wqltot(k1))
    allocate(wthvtot(k1))
    allocate(wthvres(k1))
    allocate(wthvsub(k1))
    allocate(wsvsub(k1,nsv))
    allocate(wsvres(k1,nsv))
    allocate(wsvtot(k1,nsv))
    allocate(thvhav(k1))
    allocate(th0av(k1))
    allocate(svptav(k1,nsv))
    allocate(svpav(k1,nsv))

    allocate(cszmn(k1), cszav(k1))

    allocate(qlmnlast(k1))
    allocate(wthvtmnlast(k1))

    allocate( qlhav(k1))
    allocate( u2av    (k1), &
              v2av    (k1), &
              w2av    (k1), &
              w3av    (k1), &
              w2subav (k1), &
              qt2av   (k1), &
              thl2av  (k1), &
              thv2av  (k1), &
              th2av   (k1), &
              ql2av   (k1), &
              sv2av   (k1,nsv))
    allocate(thv0(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(thvmav(k1))
    allocate(sv0h(2-ih:i1+ih,2-jh:j1+jh,k1))

    umn      = 0.
    vmn      = 0.
    thlmn    = 0.
    thvmn    = 0.
    qtmn     = 0.
    qlmn     = 0.
    qlhmn    = 0.
    cfracmn  = 0.

    wthlsmn =  0.
    wthlrmn =  0.
    wthltmn =  0.

    wthvsmn =  0.
    wthvrmn =  0.
    wthvtmn =  0.

    wqtsmn =  0.
    wqtrmn =  0.
    wqttmn =  0.

    wqlsmn =  0.
    wqlrmn =  0.
    wqltmn =  0.

    uwtmn  = 0.
    vwtmn  = 0.
    uwrmn  = 0.
    vwrmn  = 0.
    uwsmn  = 0.
    vwsmn  = 0.


    u2mn     = 0.
    v2mn     = 0.
    w2mn     = 0.
    w2submn  = 0.
    skewmn   = 0.
    qt2mn    = 0.
    thl2mn   = 0.
    thv2mn   = 0.
    th2mn    = 0.
    ql2mn    = 0.
!     qs2mn    = 0.
!     qsmn     = 0.
!     rhmn     = 0.
!     rmn      = 0.
!     r2mn     = 0.
!     r3mn     = 0.

    svmmn   = 0.
    svpmn   = 0.
    svpav   = 0.
!     svplsmn = 0.
    svptmn  = 0.
    svptav  = 0.

    sv2mn = 0.

    wsvsmn = 0.
    wsvrmn = 0.
    wsvtmn = 0.

    cszav = 0.
    cszmn = 0.

    qlmnlast = 0.
    wthvtmnlast = 0.

    if(myid==0)then
       if (.not. lwarmstart) then
          open (ifoutput,file='field.'//cexpnr,status='replace')
          close (ifoutput)
          open (ifoutput,file='flux1.'//cexpnr,status='replace')
          close (ifoutput)
          open (ifoutput,file='flux2.'//cexpnr,status='replace')
          close (ifoutput)
          open (ifoutput,file='moments.'//cexpnr,status='replace')
          close (ifoutput)
          do n=1,nsv
             name = 'svnnnfld.'//cexpnr
             write (name(3:5),'(i3.3)') n
             open (ifoutput,file=name,status='replace')
             close (ifoutput)
          end do
          do n=1,nsv
             name = 'svnnnflx.'//cexpnr
             write (name(3:5),'(i3.3)') n
             open (ifoutput,file=name,status='replace')
             close (ifoutput)
          end do
       end if
       
    if (lnetcdf) then
      fname(10:12) = cexpnr
      nvar = nvar + 7*nsv
      allocate(ncname(nvar,4))
      call nctiminfo(tncname(1,:))
      call ncinfo(ncname( 1,:),'rhof','Full level slab averaged density','kg/m^3','tt')
      call ncinfo(ncname( 2,:),'rhobf','Full level base-state density','kg/m^3','tt')
      call ncinfo(ncname( 3,:),'rhobh','Half level base-state density','kg/m^3','mt')
      call ncinfo(ncname( 4,:),'presh','Pressure at cell center','Pa','tt')
      call ncinfo(ncname( 5,:),'u','West-East velocity','m/s','tt')
      call ncinfo(ncname( 6,:),'v','South-North velocity','m/s','tt')
      call ncinfo(ncname( 7,:),'thl','Liquid water potential temperature','K','tt')
      call ncinfo(ncname( 8,:),'thv','Virtual potential temperature','K','tt')
      call ncinfo(ncname( 9,:),'qt','Total water specific humidity','kg/kg','tt')
      call ncinfo(ncname(10,:),'ql','Liquid water specific humidity','kg/kg','tt')
      call ncinfo(ncname(11,:),'wthls','SFS-Theta_l flux','Km/s','mt')
      call ncinfo(ncname(12,:),'wthlr','Resolved Theta_l flux','Km/s','mt')
      call ncinfo(ncname(13,:),'wthlt','Total Theta_l flux','Km/s','mt')
      call ncinfo(ncname(14,:),'wthvs','SFS-buoyancy flux','Km/s','mt')
      call ncinfo(ncname(15,:),'wthvr','Resolved buoyancy flux','Km/s','mt')
      call ncinfo(ncname(16,:),'wthvt','Total buoyancy flux','Km/s','mt')
      call ncinfo(ncname(17,:),'wqts','SFS-moisture flux','kg/kg m/s','mt')
      call ncinfo(ncname(18,:),'wqtr','Resolved moisture flux','kg/kg m/s','mt')
      call ncinfo(ncname(19,:),'wqtt','Total moisture flux','kg/kg m/s','mt')
      call ncinfo(ncname(20,:),'wqls','SFS-liquid water flux','kg/kg m/s','mt')
      call ncinfo(ncname(21,:),'wqlr','Resolved liquid water flux','kg/kg m/s','mt')
      call ncinfo(ncname(22,:),'wqlt','Total liquid water flux','kg/kg m/s','mt')
      call ncinfo(ncname(23,:),'uws','SFS-momentum flux (uw)','m^2/s^2','mt')
      call ncinfo(ncname(24,:),'uwr','Resolved momentum flux (uw)','m^2/s^2','mt')
      call ncinfo(ncname(25,:),'uwt','Total momentum flux (vw)','m^2/s^2','mt')
      call ncinfo(ncname(26,:),'vws','SFS-momentum flux (vw)','m^2/s^2','mt')
      call ncinfo(ncname(27,:),'vwr','Resolved momentum flux (vw)','m^2/s^2','mt')
      call ncinfo(ncname(28,:),'vwt','Total momentum flux (vw)','m^2/s^2','mt')
      call ncinfo(ncname(29,:),'w2s','SFS-TKE','m^2/s^2','mt')
      call ncinfo(ncname(30,:),'w2r','Resolved vertical velocity variance','m^2/s^2','mt')
      !call ncinfo(ncname(31,:),'w2t','Total vertical velocity variance','m^2/s^2','mt')
      call ncinfo(ncname(31,:),'skew','vertical velocity skewness','-','mt')
      call ncinfo(ncname(32,:),'u2r','Resolved horizontal velocity variance (u)','m^2/s^2','tt')
      call ncinfo(ncname(33,:),'v2r','Resolved horizontal velocity variance (v)','m^2/s^2','tt')
      call ncinfo(ncname(34,:),'thl2r','Resolved theta_l variance','K^2','tt')
      call ncinfo(ncname(35,:),'thv2r','Resolved buoyancy variance','K^2','tt')
      call ncinfo(ncname(36,:),'th2r','Resolved theta variance','K^2','tt')
      call ncinfo(ncname(37,:),'qt2r','Resolved total water variance','(kg/kg)^2','tt')
      call ncinfo(ncname(38,:),'ql2r','Resolved liquid water variance','(kg/kg)^2','tt')
      call ncinfo(ncname(39,:),'cs','Smagorinsky constant','-','tt')
      call ncinfo(ncname(40,:),'cfrac','Cloud fraction','-','tt')
      do n=1,nsv
        write (csvname(1:3),'(i3.3)') n
        call ncinfo(ncname(40+7*(n-1)+1,:),'sv'//csvname,'Scalar '//csvname//' specific mixing ratio','(kg/kg)','tt')
        call ncinfo(ncname(40+7*(n-1)+2,:),'svp'//csvname,'Scalar '//csvname//' tendency','(kg/kg/s)','tt')
        call ncinfo(ncname(40+7*(n-1)+3,:),'svpt'//csvname,'Scalar '//csvname//' turbulence tendency','(kg/kg/s)','tt')
        call ncinfo(ncname(40+7*(n-1)+4,:),'sv'//csvname//'2r','Resolved scalar '//csvname//' variance','(kg/kg)^2','tt')
        call ncinfo(ncname(40+7*(n-1)+5,:),'wsv'//csvname//'s','SFS scalar '//csvname//' flux','kg/kg m/s','mt')
        call ncinfo(ncname(40+7*(n-1)+6,:),'wsv'//csvname//'r','Resolved scalar '//csvname//' flux','kg/kg m/s','mt')
        call ncinfo(ncname(40+7*(n-1)+7,:),'wsv'//csvname//'t','Total scalar '//csvname//' flux','kg/kg m/s','mt')
      end do

      if (isurf==1) then
        call open_nc(fname,  ncid,nrec,n3=kmax,ns=ksoilmax)
      else
        call open_nc(fname,  ncid,nrec,n3=kmax)
      endif
      if (nrec == 0) then
        call define_nc( ncid, 1, tncname)
        call writestat_dims_nc(ncid)
      end if
      call define_nc( ncid, NVar, ncname)
    end if


    end if
    !$acc enter data copyin(umn, vmn, thlmn, thvmn, qtmn, qlmn, qlhmn, cfracmn, wthlsmn, wthlrmn, wthltmn, &
    !$acc&                  wthvsmn, wthvrmn, wthvtmn, wqtsmn, wqtrmn, wqttmn, wqlsmn, wqlrmn, wqltmn, &
    !$acc&                  uwtmn, vwtmn, uwrmn, vwrmn, uwsmn, vwsmn, u2mn, v2mn, w2mn, w2submn, skewmn, &
    !$acc&                  qt2mn, thl2mn, thv2mn, th2mn, ql2mn, svmmn, svpmn, svpav, svptmn, svptav, &
    !$acc&                  sv2mn, wsvsmn, wsvrmn, wsvtmn, cszav, cszmn, qlmnlast, wthvtmnlast, qlhav, &
    !$acc&                  wqlsub, wqlres, wthlsub, wthlres, wqtsub, wqtres, wthvsub, wthvres, wqttot, & 
    !$acc&                  wqltot, wthltot, wthvtot, wsvsub, wsvres, wsvtot, uwres, vwres, uwsub, vwsub, &
    !$acc&                  uwtot, vwtot, umav, vmav, thvmav, thlmav, qtmav, qlmav, cfracav, u2av, v2av, &
    !$acc&                  w2av, w2subav, qt2av, thl2av, thv2av, th2av, svmav, svpav, svptav, sv2av, w3av, &
    !$acc&                  ql2av, thvmav, thmav, qlptav, thv0)

  end subroutine initgenstat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine genstat

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
      call do_genstat
    end if
    if (timee>=tnextwrite) then
      tnextwrite = tnextwrite+itimeav
      call writestat
    end if
    dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))
  end subroutine genstat

  subroutine do_genstat

    use modfields, only : u0,v0,w0,um,vm,wm,qtm,thlm,thl0,qt0,qt0h, &
                          ql0,ql0h,thl0h,thv0h,sv0, svm, e12m,exnf,exnh
    use modsurfdata,only: thls,qts,svs,ustar,thlflux,qtflux,svflux
    use modsubgriddata,only : ekm, ekh, csz
    use modglobal, only : i1,ih,j1,jh,k1,kmax,nsv,dzf,dzh,rlv,rv,rd,cp, &
                          ijtot,cu,cv,iadv_sv,iadv_kappa,eps1,dxi,dyi
    use modmpi,    only : comm3d,mpi_sum,mpierr,slabsum,D_MPI_ALLREDUCE,myid
    use advec_kappa, only : halflev_kappa
    use modtimer, only: timer_tic, timer_toc
    implicit none

    real cthl,cqt,den

    integer i, j, k, n, km
    real    tsurf, qsat, c1, c2
    real    qs0h, t0h, ekhalf, euhalf, evhalf
    real    wthls, wthlr, wqts, wqtr, wqls, wqlr, wthvs, wthvr
    real    uws,vws,uwr,vwr
    real    upcu, vpcv
    real    qls

    real :: qlhav_s, wthlsub_s, wqtsub_s, wqlsub_s, wthvsub_s
    real :: uwsub_s, vwsub_s, uwres_s, vwres_s
    real :: wqlres_s, wthlres_s, wthvres_s, wqtres_s
    real :: wsvsub_s, wsvres_s

  !-----------------------------------------------------------------------
  !     1.    INITIALISE LOCAL CONSTANTS
  !     --    --------------------------

  !     --------------------------------------------------------
  !     3.0    RESET ARRAYS FOR SLAB AVERAGES
  !     ---    ------------------------------
  !     --------------------------------------------------------
    !$acc kernels default(present)
    qlhav = 0.0
    u2av = 0.0
    v2av = 0.0
    w2av = 0.0
    w3av = 0.0
    w2subav = 0.0
    qt2av = 0.0
    thl2av = 0.0
    thv2av = 0.0
    th2av = 0.0
    ql2av = 0.0
    thvmav = 0.0
    wqtsub = 0.0
    wqtres = 0.0
    wqlsub = 0.0
    wqlres = 0.0
    umav = 0.0
    vmav = 0.0
    thlmav = 0.0
    thmav = 0.0
    qtmav = 0.0
    qlmav = 0.0
    wthlsub = 0.0
    wthlres = 0.0
    wthvsub = 0.0
    wthvres = 0.0
    sv2av = 0.0
    uwres = 0.0
    vwres = 0.0
    uwsub = 0.0
    vwsub = 0.0
    svmav = 0.0
    qlptav = 0.0
    wqltot = 0.0
    wqttot = 0.0
    wthvtot = 0.0
    wthltot = 0.0
    uwtot = 0.0
    vwtot = 0.0
    wsvsub = 0.0
    wsvres = 0.0
    sv2av = 0.0
    cfracav= 0.0
    cszav = 0.0
    !$acc end kernels

    qlhav_s = 0.0
    wthlsub_s = 0.0
    wqtsub_s = 0.0
    wqlsub_s = 0.0
    wthvsub_s = 0.0
    uwsub_s = 0.0
    vwsub_s = 0.0
    uwres_s = 0.0
    vwres_s = 0.0
    wqlres_s = 0.0
    wthlres_s = 0.0
    wthvres_s = 0.0
    wqtres_s = 0.0
    wsvsub_s = 0.0
    wsvres_s = 0.0

    !$acc parallel loop collapse(3) default(present)
    do  k=1,k1
      do  j=2,j1
        do  i=2,i1
          thv0(i,j,k) = (thl0(i,j,k)+rlv*ql0(i,j,k)/(cp*exnf(k))) &
                        *(1+(rv/rd-1)*qt0(i,j,k)-rv/rd*ql0(i,j,k))
        enddo
      enddo
    enddo

    !$acc parallel loop default(present)
    do k=1,k1
      cfracav(k)    = cfracav(k)+count(ql0(2:i1,2:j1,k)>0)
    end do

    !$acc host_data use_device(cfracav)
    call D_MPI_ALLREDUCE(cfracav,k1,MPI_SUM,comm3d,mpierr)
    !$acc end host_data


    !$acc host_data use_device(umav, um, vmav, vm, thlmav, thlm, qtmav, qtm, &
    !$acc&                     qlmav, ql0, thvmav, thv0)
    call slabsum(umav  ,1,k1,um  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    call slabsum(vmav  ,1,k1,vm  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    call slabsum(thlmav,1,k1,thlm,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    call slabsum(qtmav ,1,k1,qtm ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    call slabsum(qlmav ,1,k1,ql0 ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    call slabsum(thvmav,1,k1,thv0,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    !$acc end host_data

    if (nsv > 0) then
      do n=1,nsv
        !$acc host_data use_device(svmav, svm)
        call slabsum(svmav(1:1,n),1,k1,svm(:,:,:,n),2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
        !$acc end host_data
      enddo
    end if
        
    !$acc kernels default(present)
    umav  = umav  /ijtot + cu
    vmav  = vmav  /ijtot + cv
    thlmav = thlmav/ijtot
    qtmav = qtmav /ijtot
    qlmav = qlmav /ijtot
    thmav  = thlmav + (rlv/cp)*qlmav/exnf
    thvmav = thvmav/ijtot
    svmav = svmav/ijtot
    cfracav = cfracav / ijtot
    !$acc end kernels

    cszav  = csz
  !------------------------------------------------------------------
  !     4     CALCULATE SLAB AVERAGED OF FLUXES AND SEVERAL MOMENTS
  !     -------------------------------------------------------------
  !        4.1 special treatment for lowest level
  !     -------------------------------------------------

    !$acc update self(exnh(1))

    qls   = 0.0 ! hj: no liquid water at the surface
    tsurf = thls*exnh(1)+(rlv/cp)*qls
    qsat  = qts - qls

    if (qls< eps1) then  ! TH: Should always be true
      c1 = 1.+(rv/rd-1)*qts
      c2 = (rv/rd-1)
    else
      c1 = (1.-qts+rv/rd*qsat*(1.+rlv/(rv*tsurf))) &
           / (1.+rlv/(rv*tsurf)*rlv/(cp*tsurf)*qsat)
      c2 = c1*rlv/(tsurf*cp)-1.
    end if

    den = 1. + (rlv**2)*qsat/(rv*cp*(tsurf**2))
    cthl = (exnh(1)*cp/rlv)*((1-den)/den)
    cqt = 1./den
        
    !$acc parallel loop collapse(2) default(present) private(upcu, vpcv) &
    !$acc& reduction(+: qlhav_s, wthlsub_s, wqtsub_s, wthvsub_s, uwsub_s, vwsub_s)
    do j = 2, j1
      do i = 2, i1
        qlhav_s = qlhav_s + ql0h(i,j,1)
        wthlsub_s = wthlsub_s + thlflux(i,j)
        wqtsub_s = wqtsub_s + qtflux (i,j)
        wthvsub_s = wthvsub_s + ( c1*thlflux(i,j)+c2*thls*qtflux(i,j) ) !hj: thv0 replaced by thls

        !Momentum flux
        upcu = um(i, j, 1) + cu
        upcu = sign(1., upcu) * max(abs(upcu), eps1)

        uwsub_s = uwsub_s - (0.5 * (ustar(i,j) + ustar(i-1,j)))**2 &
                    * upcu / sqrt(upcu**2 + ((vm(i,j,1) + vm(i-1,j,1) + vm(i,j+1,1) + vm(i-1,j+1,1)) / 4. + cv)**2)

        vpcv = vm(i, j, 1) + cv
        vpcv = sign(1., vpcv) * max(abs(vpcv), eps1)

        vwsub_s = vwsub_s - (0.5 * (ustar(i,j) + ustar(i,j-1)))**2 &
                    * vpcv / sqrt(vpcv**2 + ((um(i,j,1) + um(i+1,j,1) + um(i,j-1,1) + um(i+1,j-1,1)) / 4. + cu)**2)
      end do
    end do

    !$acc kernels default(present)
    qlhav(1) = qlhav_s
    wthlsub(1) = wthlsub_s
    wqtsub(1) = wqtsub_s
    wqlsub(1) = 0.0
    wthvsub(1) = wthvsub_s
    uwsub(1) = uwsub_s
    vwsub(1) = vwsub_s
    !$acc end kernels

  !      --------------------------
  !      4.2 higher levels
  !      --------------------------
    !$acc parallel loop gang default(present) &
    !$acc& private(qlhav_s, wqlsub_s, wqlres_s, wthlsub_s, wthlres_s, wthvsub_s, wthvres_s, &
    !$acc&         wqtsub_s, wqtres_s, uwres_s, vwres_s, uwsub_s, vwsub_s)
    do k = 2, kmax
      qlhav_s = 0.0
      wthlsub_s = 0.0
      wqtsub_s = 0.0
      wqlsub_s = 0.0
      wthvsub_s = 0.0
      uwsub_s = 0.0
      vwsub_s = 0.0
      uwres_s = 0.0
      vwres_s = 0.0
      wqlres_s = 0.0
      wthlres_s = 0.0
      wthvres_s = 0.0
      wqtres_s = 0.0
      !$acc loop collapse(2) &
      !$acc& private(qs0h, t0h, den, cthl, cqt, c1, c2, ekhalf, euhalf, evhalf, wthls, wthlr, wqts, &
      !$acc&         wqtr, wqls, wqlr, wthvs, wthvr, uwr, vwr, uws, vws) &
      !$acc& reduction(+:qlhav_s, wqlsub_s, wqlres_s, wthlsub_s, wthlres_s, wthvsub_s, wthvres_s, &
      !$acc&             wqtsub_s, wqtres_s, uwres_s, vwres_s, uwsub_s, vwsub_s)
      do j = 2, j1
        do i = 2, i1

  !       --------------------------------------------------------
  !        Calculate half level fields for thl and qt consistent
  !        with the used advection scheme ( kappa or cent. diff.)
  !       ---------------

  !       --------------------------------------------------------
    !     ------------------------------------------------------
    !     calculate ql and thv at time t0 at full and half level
    !      ----------------------------------------------------
          qlhav_s = qlhav_s  + ql0h(i,j,k)
    !       -----------------------------------------------------------
    !       calculate prefactors for subgrid wthv and wql fluxes
    !        at half levels
    !       -----------------------------------------------------------
          qs0h  =  (qt0h(i,j,k) - ql0h(i,j,k))
          t0h   =  exnh(k)*thl0h(i,j,k) + (rlv/cp)*ql0h(i,j,k)

          den   = 1. + (rlv**2)*qs0h/(rv*cp*(t0h**2))
          cthl  = (exnh(k)*cp/rlv)*((1-den)/den)
          cqt   =  (1./den)

          ! TODO: fix the branching here
          if (ql0h(i,j,k)>0) then
            c1    = (1.-qt0h(i,j,k)+rv/rd*qs0h &
                    * (1.+rd/rv*rlv/(rd*t0h)))/den
            c2    =  c1*rlv/(t0h*cp)-1.
          else
            c1 = 1. + (rv/rd-1)*qt0h(i,j,k)
            c2 = (rv/rd-1)
          end if

    !       -----------------------------------------------------------
    !       calculate resolved and subgrid fluxes at half levels
    !       -----------------------------------------------------------

          ekhalf  = (ekh(i,j,k)*dzf(k-1)+ekh(i,j,k-1)*dzf(k))/(2*dzh(k))
          euhalf = ( dzf(k-1) * ( ekm(i,j,k)  + ekm(i-1,j,k)  )  + &
                          dzf(k)  * ( ekm(i,j,k-1) + ekm(i-1,j,k-1) ) ) / &
                      ( 4.   * dzh(k) )
          evhalf = ( dzf(k-1) * ( ekm(i,j,k)  + ekm(i,j-1,k)  )  + &
                          dzf(k)  * ( ekm(i,j,k-1) + ekm(i,j-1,k-1) ) ) / &
                      ( 4.   * dzh(k) )

          wthls    = -ekhalf*(thl0(i,j,k)-thl0(i,j,k-1))/dzh(k)
          wthlr    = w0(i,j,k)*thl0h(i,j,k)

          wqts    = -ekhalf*(qt0(i,j,k)-qt0(i,j,k-1))/dzh(k)
          wqtr    = w0(i,j,k)*qt0h(i,j,k)

          wqls    = cthl*wthls+ cqt*wqts
          wqlr    = w0(i,j,k)*ql0h(i,j,k)

          wthvs    = c1*wthls + c2*thl0h(i,j,k)*wqts
          wthvr    = w0(i,j,k)*thv0h(i,j,k)

          uwr     = (w0(i,j,k)+w0(i-1,j,k)) &
                    *((u0(i,j,k-1)+cu)*dzf(k)+(u0(i,j,k)+cu)*dzf(k-1))/(4*dzh(k))
          vwr     = (w0(i,j,k)+w0(i,j-1,k)) &
                    *((v0(i,j,k-1)+cv)*dzf(k)+(v0(i,j,k)+cv)*dzf(k-1))/(4*dzh(k))
          uws     = -euhalf &
                    *((u0(i,j,k)-u0(i,j,k-1))/dzh(k)+(w0(i,j,k)-w0(i-1,j,k))*dxi)
          vws     = -evhalf &
                    *((v0(i,j,k)-v0(i,j,k-1))/dzh(k)+(w0(i,j,k)-w0(i,j-1,k))*dyi)


          if (ql0h(i,j,k)>0) then
            wqlsub_s = wqlsub_s + wqls
          end if

          wqlres_s = wqlres_s + wqlr

          wthlsub_s = wthlsub_s + wthls
          wthlres_s = wthlres_s + wthlr

          wthvsub_s = wthvsub_s + wthvs
          wthvres_s = wthvres_s + wthvr

          wqtsub_s = wqtsub_s + wqts
          wqtres_s = wqtres_s + wqtr

          uwres_s = uwres_s + uwr
          vwres_s = vwres_s + vwr
          uwsub_s = uwsub_s + uws
          vwsub_s = vwsub_s + vws
        end do
      end do
      wqlsub(k) = wqlsub_s
      wqlres(k) = wqlres_s
      wthlsub(k) = wthlsub_s
      wthlres(k) = wthlres_s
      wqtsub(k) = wqtsub_s
      wqtres(k) = wqtres_s
      uwres(k) = uwres_s
      vwres(k) = vwres_s
      uwsub(k) = uwsub_s
      vwsub(k) = vwsub_s
    end do
  !     -------------------
    call calc_moment(u2av, um, 2, 1, kmax, 2, i1, 2, j1, umav, cu)
    call calc_moment(v2av, vm, 2, 1, kmax, 2, i1, 2, j1, vmav, cv)
    call calc_moment(w2av, wm, 2, 1, kmax, 2, i1, 2, j1)
    call calc_moment(w3av, wm, 3, 1, kmax, 2, i1, 2, j1)
    call calc_moment(w2subav, e12m, 2, 1, kmax, 2, i1, 2, j1)
    call calc_moment(qt2av, qtm, 2, 1, kmax, 2, i1, 2, j1, qtmav)
    call calc_moment(thl2av, thlm, 2, 1, kmax, 2, i1, 2, j1, thlmav)
    call calc_moment(thv2av, thv0, 2, 1, kmax, 2, i1, 2, j1, thvmav)
    call calc_moment(th2av, thlm, 2, 1, kmax, 2, i1, 2, j1, thmav)
    call calc_moment(ql2av, ql0, 2, 1, kmax, 2, i1, 2, j1, qlmav)

    if (nsv > 0) then
      do n = 1, nsv
        call calc_moment(sv2av(:, n), svm(:, :, :, n), 2, 1, kmax, 2, i1, 2, j1, svmav(:, n))
      end do

      do n = 1, nsv
        if (iadv_sv(n)==iadv_kappa) then
           call halflev_kappa(sv0(2-ih:i1+ih,2-jh:j1+jh,1:k1,n),sv0h)
        else
          !$acc parallel loop collapse(3) default(present)
          do k = 2, k1
            do j = 2, j1
              do i = 2, i1
                sv0h(i,j,k) = (sv0(i,j,k,n)*dzf(k-1)+sv0(i,j,k-1,n)*dzf(k))/(2*dzh(k))
              enddo
            enddo
          enddo
          !$acc kernels default(present)
          sv0h(2:i1,2:j1,1) = svs(n)
          !$acc end kernels
        end if

        !$acc parallel loop default(present) private(wsvres_s)
        do k = 2, kmax
          wsvres_s = 0.0
          !$acc loop collapse(2) reduction(+: wsvres)
          do j = 2, j1
            do i = 2, i1
              wsvres_s = wsvres_s + w0(i,j,k)*sv0h(i,j,k)
            end do
          end do
          wsvres(k,n) = wsvres_s
        end do

        wsvsub_s = 0.0
        !$acc kernels default(present)
        do j = 2, j1
          do i = 2, i1
            wsvsub_s = wsvsub_s + svflux(i,j,n)
          end do
        end do
        wsvsub(1,n) = wsvsub_s
        !$acc end kernels

        !$acc parallel loop private(wsvsub_s)
        do k = 2, kmax
          wsvsub_s = 0.0
          !$acc loop collapse(2) private(ekhalf) reduction(+: wsvsub_s)
          do j = 2, j1
            do i= 2, i1
              ekhalf = (ekh(i,j,k)*dzf(k-1)+ekh(i,j,k-1)*dzf(k))/(2*dzh(k))
              wsvsub_s= wsvsub_s-ekhalf*(sv0(i,j,k,n)-sv0(i,j,k-1,n)) / dzh(k)
            end do
          end do
          wsvsub(k,n) = wsvsub_s
        end do
      end do
    end if
  !     -------------------------------
  !     5   CALCULATE MOMENTUM FLUXES
  !     -------------------------------

  !     5.1 special treatment for lowest level
  !     -------------------------------------------------
  !      DEPRECATED


  !     5.2 higher levels by vert. integr. of the mom. tendencies
  !     ---------------------------------------------------------
  !         DEPRECATED

  ! MPI communication
    !$acc host_data use_device(qlhav, wqlsub, wqlres, wthlsub, wthlres, wthvsub, &
    !$acc&                     wthvres, uwsub, vwsub, uwres, vwres, u2av, v2av, &
    !$acc&                     w2av, w3av, w2subav, qt2av, thl2av, thv2av, th2av, &
    !$acc&                     ql2av, qlptav, sv2av, wsvsub, wsvres)
    call D_MPI_ALLREDUCE(qlhav, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(wqlsub, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(wqlres, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(wthlsub, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(wthlres, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(wqtsub, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(wqtres, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(wthvsub, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(wthvres, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(uwsub, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(vwsub, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(uwres, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(vwres, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(u2av, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(v2av, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(w2av, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(w3av, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(w2subav, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(qt2av, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(thl2av, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(thv2av, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(th2av, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(ql2av, k1, MPI_SUM, comm3d,mpierr)
!     call D_MPI_ALLREDUCE(qs2avl, qs2av, k1,     &
!                       MPI_SUM, comm3d,mpierr)
!     call D_MPI_ALLREDUCE(qsavl, qsav, k1,     &
!                       MPI_SUM, comm3d,mpierr)
!     call D_MPI_ALLREDUCE(rhavl, rhav, k1,     &
!                       MPI_SUM, comm3d,mpierr)
!     call D_MPI_ALLREDUCE(ravl, rav, k1,     &
!                       MPI_SUM, comm3d,mpierr)
!     call D_MPI_ALLREDUCE(r2avl, r2av, k1,     &
!                       MPI_SUM, comm3d,mpierr)
!     call D_MPI_ALLREDUCE(r3avl, r3av, k1,     &
!                       MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(qlptav, k1, MPI_SUM, comm3d,mpierr)

    if (nsv > 0) then
      do n=1,nsv
        call D_MPI_ALLREDUCE(sv2av(:,n),k1, MPI_SUM, comm3d,mpierr)
        call D_MPI_ALLREDUCE(wsvsub(:,n), k1, MPI_SUM, comm3d,mpierr)
        call D_MPI_ALLREDUCE(wsvres(:,n), k1, MPI_SUM, comm3d,mpierr)
      end do
    end if
    !$acc end host_data

  !     -----------------------------------------------
  !     6   NORMALIZATION OF THE FIELDS AND FLUXES
  !     -----------------------------------------------
      !$acc kernels default(present)
      qlhav   = qlhav  /ijtot

      wqlsub  = wqlsub /ijtot
      wqlres  = wqlres /ijtot

      wthlsub  = wthlsub /ijtot
      wthlres  = wthlres /ijtot

      wqtsub  = wqtsub /ijtot
      wqtres  = wqtres /ijtot

      wthvsub  = wthvsub /ijtot
      wthvres  = wthvres /ijtot

      wqttot  = wqtres + wqtsub
      wqltot  = wqlres + wqlsub
      wthltot  = wthlres + wthlsub
      wthvtot  = wthvres + wthvsub

      if (nsv > 0) then
        wsvsub = wsvsub /ijtot
        wsvres = wsvres /ijtot
        wsvtot = wsvsub + wsvres
      end if

      uwres    = uwres    /ijtot
      vwres    = vwres    /ijtot
      uwsub    = uwsub    /ijtot
      vwsub    = vwsub    /ijtot
      uwtot    = uwres + uwsub
      vwtot    = vwres + vwsub
      !$acc end kernels

!       qs2av    = qs2av    /ijtot
!       qsav     = qsav     /ijtot
!       qs2av    = qs2av    - qsav**2
!       rhav     = rhav     /ijtot
!       rav      = rav      /ijtot
!       r2av     = r2av     /ijtot
!       r3av     = r3av     /ijtot

  !********************************************************************

  !     4.0   ADD SLAB AVERAGES TO TIME MEAN
  !           ------------------------------
      !$acc kernels default(present)
      umn    = umn   + umav
      vmn    = vmn   + vmav
      thvmn  = thvmn + thvmav
      thlmn  = thlmn + thlmav
      qtmn   = qtmn  + qtmav
      qlmn   = qlmn  + qlmav
      cfracmn= cfracmn+cfracav
      qlhmn  = qlhmn + qlhav

      wthlsmn = wthlsmn + wthlsub
      wthlrmn = wthlrmn + wthlres
      wthltmn = wthltmn + wthltot
      wthvsmn = wthvsmn + wthvsub
      wthvrmn = wthvrmn + wthvres
      wthvtmn = wthvtmn + wthvtot
      wqtsmn = wqtsmn + wqtsub
      wqtrmn = wqtrmn + wqtres
      wqttmn = wqttmn + wqttot
      wqlsmn = wqlsmn + wqlsub
      wqlrmn = wqlrmn + wqlres
      wqltmn = wqltmn + wqltot
      uwtmn  = uwtmn + uwtot
      vwtmn  = vwtmn + vwtot
      uwrmn  = uwrmn + uwres
      vwrmn  = vwrmn + vwres
      uwsmn  = uwsmn + uwsub
      vwsmn  = vwsmn + vwsub
      u2mn     = u2mn     + u2av
      v2mn     = v2mn     + v2av
      w2mn     = w2mn     + w2av
      w2submn  = w2submn  + w2subav
      qt2mn    = qt2mn    + qt2av
      thl2mn   = thl2mn   + thl2av
      thv2mn   = thv2mn   + thv2av
      th2mn    = th2mn    + th2av
!       ql2mn    = ql2mn    + ql2av
!       qs2mn    = qs2mn    + qs2av
!       qsmn     = qsmn     + qsav
!       rhmn     = rhmn     + rhav
!       rmn      = rmn      + rav
!       r2mn     = r2mn     + r2av
!       r3mn     = r3mn     + r3av
      if (nsv > 0) then
        svmmn   = svmmn  + svmav
        svpmn   = svpmn  + svpav
!         svplsmn = svplsmn+ svplsav
        svptmn  = svptmn + svptav

        sv2mn  = sv2mn + sv2av

        wsvsmn = wsvsmn + wsvsub
        wsvrmn = wsvrmn + wsvres
        wsvtmn = wsvtmn + wsvtot
      end if
      skewmn   = skewmn   + w3av/max(w2av**1.5,epsilon(w2av(1)))

      cszmn = cszmn + cszav
      !$acc end kernels

  end subroutine do_genstat

  subroutine calc_moment(prof, var, n, kb, ke, ib, ie, jb, je, mean, c_in)
    use modglobal, only: ijtot

    implicit none

    integer, intent(in) :: kb, ke, ib, ie, jb, je
    integer, intent(in) :: n
    real, intent(out) :: prof(kb:ke)
    real(field_r), intent(in) :: var(:, :, :)
    real(field_r), optional, intent(in) :: mean(kb:ke)
    real(field_r), optional, intent(in) :: c_in !< Translational velocity
    real(field_r) :: c 
    real(field_r) :: prof_s
    integer :: i, j, k

    if (.not.present(c_in)) then
      c = 0.
    else
      c = c_in
    end if

    if (.not.present(mean)) then
      !$acc parallel loop default(present) private(prof_s)
      do k = kb, ke
        !$acc loop collapse(2) reduction(+: prof_s)
        do j = jb, je
          do i = ib, ie
            prof_s = prof_s + (var(i, j, k) + c)**n
          end do
        end do 
        prof(k) = prof_s / ijtot
      end do
    else
      !$acc parallel loop default(present) private(prof_s)
      do k = kb, ke
        !$acc loop collapse(2) reduction(+: prof_s)
        do j = jb, je
          do i = ib, ie
            prof_s = prof_s + (var(i, j, k) + c - mean(k))**n
          end do
        end do 
        prof(k) = prof_s / ijtot
      end do
    end if

  end subroutine calc_moment

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine writestat
      use modglobal, only : kmax,k1,nsv, zh,zf,rtimee,rlv,cp,cexpnr,ifoutput
      use modfields, only : presf,presh,exnf,exnh,rhof,rhobf,rhobh
      use modsubgriddata, only : csz
      use modmpi,    only : myid
      use modstat_nc, only: lnetcdf, writestat_nc
      implicit none


      real,dimension(k1,nvar) :: vars
      real,allocatable, dimension(:) :: tmn, thmn
      integer nsecs, nhrs, nminut,k,n
      real convt, convq
      character(40) :: name

      nsecs   = nint(rtimee)
      nhrs    = int(nsecs/3600)
      nminut  = int(nsecs/60)-nhrs*60
      nsecs   = mod(nsecs,60)
      convt   = 86400.
      convq   = 86400*1000.
      allocate(tmn   (k1), thmn  (k1))

      !$acc kernels default(present)
      umn    = umn    /nsamples
      vmn    = vmn    /nsamples
      thvmn  = thvmn  /nsamples
      thlmn  = thlmn  /nsamples
      qtmn   = qtmn   /nsamples
      qlmn   = qlmn   /nsamples
      cfracmn= cfracmn/nsamples
      qlhmn  = qlhmn  /nsamples


      wthlsmn = wthlsmn/nsamples
      wthlrmn = wthlrmn/nsamples
      wthltmn = wthltmn/nsamples

      wqtsmn = wqtsmn/nsamples
      wqtrmn = wqtrmn/nsamples
      wqttmn = wqttmn/nsamples

      wqlsmn = wqlsmn/nsamples
      wqlrmn = wqlrmn/nsamples
      wqltmn = wqltmn/nsamples

      wthvsmn = wthvsmn/nsamples
      wthvrmn = wthvrmn/nsamples
      wthvtmn = wthvtmn/nsamples

      uwtmn  = uwtmn /nsamples
      vwtmn  = vwtmn /nsamples
      uwrmn  = uwrmn /nsamples
      vwrmn  = vwrmn /nsamples
      uwsmn  = uwsmn /nsamples
      vwsmn  = vwsmn /nsamples

!       rmn      = rmn    /nsamples
!       r2mn     = r2mn   /nsamples
!       rhmn     = rhmn   /nsamples
      w2mn     = w2mn   /nsamples
      skewmn   = skewmn /nsamples
      w2submn  = w2submn/nsamples
      qt2mn    = qt2mn  /nsamples
      v2mn     = v2mn   /nsamples
      u2mn     = u2mn   /nsamples
      thl2mn   = thl2mn /nsamples
      thv2mn   = thv2mn /nsamples
      th2mn    = th2mn  /nsamples
      ql2mn    = ql2mn  /nsamples
!       qs2mn    = qs2mn  /nsamples
!       qsmn     = qsmn   /nsamples

      if (nsv > 0) then
        svmmn   = svmmn  /nsamples
        svpmn   = svpmn  /nsamples
!         svplsmn = svplsmn/nsamples
        svptmn  = svptmn /nsamples

        sv2mn = sv2mn/nsamples

        wsvsmn = wsvsmn/nsamples
        wsvrmn = wsvrmn/nsamples
        wsvtmn = wsvtmn/nsamples
      end if

      cszmn = cszmn / nsamples
      !$acc end kernels


  !     ------------------------------------------
  !     2.0  Construct other time averaged fields
  !     ------------------------------------------

      !$acc kernels default(present) copy(thmn, tmn)
      thmn = thlmn + (rlv/cp)*qlmn/exnf
      tmn  = thmn*exnf
      !$acc end kernels

      !$acc update self(umn, vmn, thvmn, thlmn, qtmn, qlmn, cfracmn, qlhmn, &
      !$acc&            wthlsmn, wthlrmn, wthltmn, wqtsmn, wqtrmn, wqttmn, &
      !$acc&            wqlsmn, wqlrmn, wqltmn, wthvsmn, wthvrmn, wthvtmn, &
      !$acc&            uwtmn, vwtmn, uwrmn, vwrmn, uwsmn, vwsmn, w2mn, skewmn, &
      !$acc&            w2submn, qt2mn, v2mn, u2mn, thl2mn, thv2mn, th2mn, ql2mn, &
      !$acc&            cszmn, cfracmn)

      !$acc update self(svmmn, svpmn, svptmn, sv2mn, wsvsmn, wsvrmn, wsvtmn) if(nsv > 0)
        
  !     ----------------------
  !     2.0  write the fields
  !           ----------------

    if(myid==0)then
      open (ifoutput,file='field.'//cexpnr,position='append')
      write(ifoutput,'(//A,/A,F5.0,A,I4,A,I2,A,I2,A)') &
      '#--------------------------------------------------------'      &
      ,'#',(timeav),'--- AVERAGING TIMESTEP --- '      &
      ,nhrs,':',nminut,':',nsecs      &
      ,'   HRS:MIN:SEC AFTER INITIALIZATION '
      write (ifoutput,'(A/2A/2A)') &
          '#--------------------------------------------------------' &
          ,'#LEV  HGHT    PRES    TEMP       TH_L     THETA      TH_V     ' &
          ,'  QT_AV      QL_AV      U       V   CLOUD FRACTION  CS' &
          ,'#      (M)    (MB)   (----------- (KELVIN) ---------------)    ' &
          ,'(----(G/KG)------)  (--- (M/S ---)   (-----------)  (---)'
      do k=1,kmax
        write(ifoutput,'(I3,F10.2,F7.1,5F10.4,F12.5,3F11.4,F11.5)') &
            k, &
            zf    (k), &
            presf (k)/100., &
            tmn   (k), &
            thlmn (k), &
            thmn  (k), &
            thvmn (k), &
            qtmn  (k)*1000., &
            qlmn  (k)*1000., &
            umn   (k), &
            vmn   (k), &
            cfracmn(k), &
            cszmn(k)
      end do
      close (ifoutput)

  !     -------------------------
  !     6.5   write the fluxes
  !     -----------------------------------------------------------------------

      open (ifoutput,file='flux1.'//cexpnr,position='append')

      write(ifoutput,'(//2A,/A,F5.0,A,I4,A,I2,A,I2,A)') &
            '#-------------------------------------------------------------' &
          ,'---------------------------------)' &
            ,'#',(timeav),'--- AVERAGING TIMESTEP --- ' &
            ,nhrs,':',nminut,':',nsecs &
            ,'   HRS:MIN:SEC AFTER INITIALIZATION '
      write (ifoutput,'(2A/A/2A/2A/2A/2A)') &
            '#---------------------------------------------------------------' &
          ,'---------------------------------)' &
          ,'#                                                               ' &
          ,'#                  |                                  TURBULENT ' &
          ,'FLUXES                           |                             ' &
          ,'#LEV HEIGHT  PRES  |   WTHL_SUB    WTHL_RES    WTHL_TOT        WQ' &
          ,'T_SUB      WQT_RES     WQT_TOT   ' &
          ,'#     (M)    (MB) | (----------   (K M/S)   --------------)     ' &
          ,'(---------- (M/S)   -------------)                             ' &
          ,'#---------------------------------------------------------------' &
          ,'---------------------------------)'

      write(ifoutput,'(I3,F10.2,F7.1,6E13.5)') &
            (k, &
            zh       (k), &
            presh     (k)/100., &
            wthlsmn   (k)                              , &
            wthlrmn   (k)                              , &
            wthltmn   (k)                              , &
            wqtsmn   (k)                              , &
            wqtrmn   (k)                              , &
            wqttmn   (k),&
             k=1,kmax)



      close(ifoutput)

      open (ifoutput,file='flux2.'//cexpnr,position='append')

      write(ifoutput,'(//A,/A,F5.0,A,I4,A,I2,A,I2,A)') &
            '#--------------------------------------------------------' &
            ,'#',(timeav),'--- AVERAGING TIMESTEP --- ' &
            ,nhrs,':',nminut,':',nsecs &
            ,'   HRS:MIN:SEC AFTER INITIALIZATION '
      write (ifoutput,'(A/A/A/3A/4A/2A)') &
            '#(------------------------------------------------------------)' &
          ,'#                                                             ' &
          ,'#                                     TURBULENT FLUXES        ' &
          ,'#LEV HEIGHT  PRES  |    UW_TOT      VW_TOT       UW_SGS   ' &
          ,'   VW_SGS       UW_RES       VW_RES ' &
          ,'      WTH_TOT       WQ_L       WTHV_SUB     WTHV_RES     WTHV_TOT' &
          ,'#     (M)    (MB) |   ' &
          ,'(---------------------------- (M/S)^2  ---------------------------------)' &
          ,'   (-(K M/S)-)   (-(' &
          ,'M/S)-)     (----------  (K M/S)  -----------)' &
          ,'#(------------------------------------------------------------' &
          ,'-----------------------------------------)'

      write(ifoutput,'(I3,F10.2,F7.1,11E13.5)') &
          (k, &
            zh       (k), &
            presh     (k)/100., &
            uwtmn    (k)                              , &
            vwtmn    (k)                              , &
            uwsmn    (k)                              , &
            vwsmn    (k)                              , &
            uwrmn    (k)                              , &
            vwrmn    (k)                              , &
            wthltmn   (k) + wqltmn(k)*(rlv/cp)/exnh(k) , &
            wqltmn   (k)                              , &
            wthvsmn   (k)                              , &
            wthvrmn   (k)                              , &
            wthvtmn   (k)                              , &
            k=1,kmax)


      close(ifoutput)

      open (ifoutput,file='moments.'//cexpnr,position='append')

      write(ifoutput,'(//A,/A,F5.0,A,I4,A,I2,A,I2,A)') &
        '#--------------------------------------------------------'      &
        ,'#',(timeav),'--- AVERAGING TIMESTEP --- '      &
        ,nhrs,':',nminut,':',nsecs      &
        ,'   HRS:MIN:SEC AFTER INITIALIZATION '

      write (ifoutput,'(3A/3A/3A)') &
            '#----------------------------------------------------' &
            ,'---------------------------------------------------' &
            ,'------------------------------' &
            ,'#  LEV   HGHT   PRES     THL**2       THV**2         ' &
            ,'TH**2       QT**2       U*U    ' &
            ,'  V*V     HGHT     W*W     SKEWW     SFS-TKE' &
            ,'#        (M)   (MB)     (--------------(K*K)--------' &
            ,'------)     (G/KG)^2     (-(M/S)' &
            ,'^2)-)      (M)     (M/S)^2  ()        (M/S)^2'

      write(ifoutput,'(I5,F10.2,F7.1,4E13.5,2F9.4,F10.2,3F9.4)') &
            (k, &
            zf    (k), &
            presf (k)/100., &
            thl2mn(k), &
            thv2mn(k), &
            th2mn (k), &
            qt2mn (k)*1e6, &
!             qs2mn (k)*1e6, &
            u2mn  (k), &
            v2mn  (k), &
            zh    (k), &
            w2mn  (k), &
            skewmn(k), &
            w2submn(k), &
            k=1,kmax)

      close(ifoutput)


  !----   Write information about scalar field and its tendencies ------

      do n=1,nsv
        name = 'svnnnfld.'//cexpnr
        write (name(3:5),'(i3.3)') n
        open (ifoutput,file=name,position='append')

        write(ifoutput,'(//2A,/A,F5.0,A,I4,A,I2,A,I2,A)') &
            '#-------------------------------------------------------------' &
            ,'---------------------------------)' &
            ,'#',(timeav),'--- AVERAGING TIMESTEP --- ' &
            ,nhrs,':',nminut,':',nsecs &
            ,'   HRS:MIN:SEC AFTER INITIALIZATION '

        write (ifoutput,'(2A/A/A/A,I2.2,A,/A/A/2A)') &
            '#--------------------------------------------------------- ' &
            ,'--------------' &
            ,'#               --FIELD &  T E N D E N C I E S  --------    ' &
            ,'#                                                           ' &
            ,'#                  |        SCALAR(' &
            ,n &
            ,')                  |' &
            ,'# LEV HEIGHT  PRES      SV(1)      TURB    TOTAL  |' &
            ,'#      (M)   (MB)  |       -----  (KG/KG/DAY) ----- |' &
            ,'#----------------------------------------------------------' &
            ,'-------------'
        write(ifoutput,'(I4,2F10.2,E14.5e3,2F10.2,E14.5E3)') &
            (k, &
              zf       (k), &
              presf    (k)/100., &
              svmmn    (k,n), &
              svptmn   (k,n)  *convt, &
              svpmn    (k,n)  *convt, &
              sv2mn    (k,n), &
              k=1,kmax)

        close(ifoutput)



  !        -----------------------
  !        Write the scalar fluxes
  !        -----------------------
        name = 'svnnnflx.'//cexpnr
        write (name(3:5),'(i3.3)') n
        open (ifoutput,file=name,position='append')

        write(ifoutput,'(//2A,/A,F5.0,A,I4,A,I2,A,I2,A)') &
            '#-------------------------------------------------------------' &
            ,'---------------------------------)' &
            ,'#',(timeav),'--- AVERAGING TIMESTEP --- ' &
            ,nhrs,':',nminut,':',nsecs &
            ,'   HRS:MIN:SEC AFTER INITIALIZATION '
        write (ifoutput,'(2A/A/A,I2.2,A,/A/A/2A)') &
            '#---------------------------------------------------------------' &
            ,'---------------------------------)' &
            ,'#                                                               ' &
            ,'#                 |         TURBULENT FLUXES (SV=',n,')       |' &
            ,'#LEV HEIGHT  PRES       WSV_SUB     WSV_RES     WSV_TOT' &
            ,'#     (M)    (MB) | (----------   (KG/KG M/S)   --------------) ' &
            ,'#---------------------------------------------------------------' &
            ,'---------------------------------)'

        write(ifoutput,'(I3,2F10.2,3E14.5E3)') &
            (k, &
              zh       (k), &
              presh    (k)/100., &
              wsvsmn   (k,n)                            , &
              wsvrmn   (k,n)                            , &
              wsvtmn   (k,n)                            , &
              k=1,kmax)

         close(ifoutput)

      end do
      if (lnetcdf) then
        vars(:, 1)=rhof
        vars(:, 2)=rhobf
        vars(:, 3)=rhobh
        vars(:, 4)=presh
        vars(:, 5)=umn
        vars(:, 6)=vmn
        vars(:, 7)=thlmn
        vars(:, 8)=thvmn
        vars(:, 9)=qtmn
        vars(:,10)=qlmn
        vars(:,11)=wthlsmn
        vars(:,12)=wthlrmn
        vars(:,13)=wthltmn
        vars(:,14)=wthvsmn
        vars(:,15)=wthvrmn
        vars(:,16)=wthvtmn
        vars(:,17)=wqtsmn
        vars(:,18)=wqtrmn
        vars(:,19)=wqttmn
        vars(:,20)=wqlsmn
        vars(:,21)=wqlrmn
        vars(:,22)=wqltmn
        vars(:,23)=uwsmn
        vars(:,24)=uwrmn
        vars(:,25)=uwtmn
        vars(:,26)=vwsmn
        vars(:,27)=vwrmn
        vars(:,28)=vwtmn
        vars(:,29)=w2submn
        vars(:,30)=w2mn
        !vars(:,31)=w2submn+w2mn
        vars(:,31)=skewmn
        vars(:,32)=u2mn
        vars(:,33)=v2mn
        vars(:,34)=thl2mn
        vars(:,35)=thv2mn
        vars(:,36)=th2mn
        vars(:,37)=qt2mn
        vars(:,38)=ql2mn
        vars(:,39)=csz
        vars(:,40)=cfracmn
        do n=1,nsv
          vars(:,40+7*(n-1)+1)=svmmn(:,n)
          vars(:,40+7*(n-1)+2)=svpmn(:,n)
          vars(:,40+7*(n-1)+3)=svptmn(:,n)
          vars(:,40+7*(n-1)+4)=sv2mn(:,n)
          vars(:,40+7*(n-1)+5)=wsvsmn(:,n)
          vars(:,40+7*(n-1)+6)=wsvrmn(:,n)
          vars(:,40+7*(n-1)+7)=wsvtmn(:,n)
        end do
        call writestat_nc(ncid,1,tncname,(/rtimee/),nrec,.true.)
        call writestat_nc(ncid,nvar,ncname,vars(1:kmax,:),nrec,kmax)
      end if

    end if ! end if(myid==0)
      
      !$acc kernels default(present)
      qlmnlast=qlmn
      wthvtmnlast=wthvtmn

      umn      = 0.
      vmn      = 0.
      thlmn    = 0.
      thvmn    = 0.
      qtmn     = 0.
      qlmn     = 0.
      qlhmn    = 0.
      cfracmn  = 0.

      wthlsmn =  0.
      wthlrmn =  0.
      wthltmn =  0.

      wthvsmn =  0.
      wthvrmn =  0.
      wthvtmn =  0.

      wqtsmn =  0.
      wqtrmn =  0.
      wqttmn =  0.

      wqlsmn =  0.
      wqlrmn =  0.
      wqltmn =  0.

      uwtmn  =  0.
      vwtmn  =  0.
      uwrmn  =  0.
      vwrmn  =  0.
      uwsmn  =  0.
      vwsmn  =  0.


      u2mn     = 0.
      v2mn     = 0.
      w2mn     = 0.
      w2submn  = 0.
      skewmn   = 0.
      qt2mn    = 0.
      thl2mn   = 0.
      thv2mn   = 0.
      th2mn    = 0.
      ql2mn    = 0.
!       qs2mn    = 0.
!       qsmn     = 0.
!       rhmn     = 0.
!       rmn      = 0.
!       r2mn     = 0.
!       r3mn     = 0.

      svmmn   = 0.
      svpmn   = 0.
!       svplsmn = 0.
      svptmn  = 0.

      sv2mn = 0.

      wsvsmn = 0.
      wsvrmn = 0.
      wsvtmn = 0.

      cszmn  = 0.
      !$acc end kernels

      deallocate(tmn, thmn)


  end subroutine writestat
  subroutine exitgenstat
    use modmpi, only : myid
    use modstat_nc, only : exitstat_nc,lnetcdf
    implicit none

    if(.not.(lstat)) return

    if(lnetcdf .and. myid==0) call exitstat_nc(ncid)

    deallocate(umn       ,vmn   )
    deallocate(thlmn        ,thvmn )
    deallocate(qtmn         ,qlmn  ,  qlhmn, cfracmn)
    deallocate(wthlsmn ,wthlrmn ,wthltmn)
    deallocate(wthvsmn ,wthvrmn ,wthvtmn)
    deallocate(wqlsmn ,wqlrmn ,wqltmn)
    deallocate(wqtsmn ,wqtrmn ,wqttmn)
    deallocate(wsvsmn ,wsvrmn,wsvtmn)
    deallocate(uwtmn,vwtmn,uwrmn,vwrmn,uwsmn,vwsmn )
!     deallocate(rmn, r2mn   , r3mn , rhmn )
    deallocate(w2mn       , skewmn )
    deallocate(w2submn    )
    deallocate(u2mn       , v2mn  ,     qt2mn)
    deallocate(thl2mn     , thv2mn,     th2mn,     ql2mn)
!     deallocate(qs2mn      , qsmn  )
    deallocate(svmmn,svptmn,svplsmn,svpmn)
    deallocate(sv2mn)

    deallocate(umav )
    deallocate(vmav )
    deallocate(thlmav )
    deallocate(thmav )
    deallocate(qtmav )
    deallocate(qlmav )
    deallocate(cfracav )
    deallocate(svmav )
    deallocate(uptav)
    deallocate(vptav)

    deallocate(thptav)
    deallocate(qlptav)
    deallocate(uwtot )
    deallocate(vwtot )
    deallocate(uwres )
    deallocate(vwres )
    deallocate(uwsub )
    deallocate(vwsub )
    deallocate(wthlsub)
    deallocate(wthlres)
    deallocate(wthltot)
    deallocate(wqtsub)
    deallocate(wqtres)
    deallocate(wqttot)
    deallocate(wqlsub)
    deallocate(wqlres)
    deallocate(wqltot)
    deallocate(wthvtot)
    deallocate(wthvres)
    deallocate(wthvsub)
    deallocate(wsvsub)
    deallocate(wsvres)
    deallocate(wsvtot)
    deallocate(thvhav)
    deallocate(th0av)
    deallocate(svpav)
    deallocate(svptav)

    deallocate(cszmn)
    deallocate(cszav)

    deallocate(qlmnlast)
    deallocate(wthvtmnlast)

    deallocate( qlhav)
    deallocate( u2av    , &
              v2av    , &
              w2av    , &
              w3av    , &
              w2subav , &
              qt2av   , &
              thl2av  , &
              thv2av  , &
              th2av   , &
              ql2av   , &
              sv2av   )
    deallocate(thv0)
    deallocate(thvmav)
    deallocate(sv0h)

  end subroutine exitgenstat


end module modgenstat

