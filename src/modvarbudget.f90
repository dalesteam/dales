! -----------------------------------------------------
!   compute variance budgets
!   TODO: Add support for doing covariance budget calculations
!   TODO: Add all modelled scalars
! -----------------------------------------------------

module modvarbudget

  use modglobal, only : longint
  implicit none
  PRIVATE
  PUBLIC :: initvarbudget, varbudget, exitvarbudget
  save
  
  ! netCDF variables 
  integer,parameter :: nvar = 14
  character(80),dimension(nvar,4) :: ncname

  real    :: dtav, timeav
  integer(kind=longint) :: idtav, itimeav,tnext,tnextwrite
  integer :: nsamples
  logical :: lvarbudget= .false. ! switch for variance budgets

  real, allocatable  :: thl2fav    (:)                 ! variance at full level &
  real, allocatable  :: thl2Prfav  (:), thl2Prfmn(:)   ! resolved production of variance at full level &
  real, allocatable  :: thl2Psfav  (:), thl2Psfmn(:)   ! subgrid production of variance at full level &
  real, allocatable  :: thl2Trfav  (:), thl2Trfmn(:)   ! resolved transport of variance at full level &
  real, allocatable  :: thl2Disfav (:), thl2Disfmn(:)  ! dissipation of thl variance at full level &
  real, allocatable  :: thl2Sfav   (:), thl2Sfmn(:)    ! source of variance at full level &
  real, allocatable  :: thl2tendf  (:)                 ! tendency of thl variance at full level &
  real, allocatable  :: thl2residf (:)                 ! residual of thl variance at full level &
  real, allocatable  :: thl2bf     (:)                 ! variance at beginning of averaging period &

  real, allocatable  :: qt2fav    (:)                 ! variance at full level &
  real, allocatable  :: qt2Prfav  (:), qt2Prfmn(:)    ! resolved production of variance at half level &
  real, allocatable  :: qt2Psfav  (:), qt2Psfmn(:)    ! subgrid production of variance at full level &
  real, allocatable  :: qt2Trfav  (:), qt2Trfmn(:)    ! resolved transport of variance at full level &
  real, allocatable  :: qt2Disfav (:), qt2Disfmn(:)   ! dissipation of qt variance at full level &
  real, allocatable  :: qt2Sfav   (:), qt2Sfmn(:)     ! source of variance at full level &
  real, allocatable  :: qt2tendf  (:)                 ! tendency of qt variance at full level &
  real, allocatable  :: qt2residf (:)                 ! tendency of qt variance at full level &
  real, allocatable  :: qt2bf     (:)                 ! variance at beginning of averaging period &

contains

  subroutine initvarbudget
    use mpi
    use modmpi,     only : myid,mpierr, comm3d,my_real, mpi_logical
    use modglobal,  only : k1,ih,i1,jh,j1,ifnamopt,fname_options, ifoutput,&
                           cexpnr,dtav_glob,timeav_glob,dt_lim,btime,tres,&
                           lwarmstart,checknamelisterror,ladaptive,dtmax
    use modstat_nc, only : lnetcdf,define_nc,ncinfo,nctiminfo,writestat_dims_nc
    use modgenstat, only : idtav_prof=>idtav, itimeav_prof=>itimeav,ncid_prof=>ncid

    implicit none

    integer ierr
    namelist/NAMVARBUDGET/ &
    dtav,timeav,lvarbudget

    dtav=dtav_glob;timeav=timeav_glob

    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMVARBUDGET,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMVARBUDGET')
      write(6 ,NAMVARBUDGET)
      close(ifnamopt)
    end if

    call MPI_BCAST(timeav     ,1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(dtav       ,1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(lvarbudget ,1,MPI_LOGICAL,0,comm3d,mpierr)
    idtav = dtav/tres
    itimeav = timeav/tres

    tnext      = idtav   +btime
    tnextwrite = itimeav +btime
    nsamples = itimeav/idtav
    if(.not.(lvarbudget)) return
    dt_lim = min(dt_lim,tnext)

    if (abs(timeav/dtav-nsamples)>1e-4) then
      stop 'timeav must be a integer multiple of dtav'
    end if
    if (.not. ladaptive .and.abs( dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'dtav should be a integer multiple of dtmax'
    end if

    allocate(thl2fav   (k1))
    allocate(thl2Prfav (k1),thl2Prfmn(k1),thl2Psfav (k1),thl2Psfmn (k1))
    allocate(thl2Trfav (k1),thl2Trfmn(k1),thl2Disfav(k1),thl2Disfmn(k1))
    allocate(thl2Sfav  (k1),thl2Sfmn (k1))
    allocate(thl2tendf (k1))
    allocate(thl2residf(k1))
    allocate(thl2bf    (k1))

    allocate(qt2fav   (k1))
    allocate(qt2Prfav (k1),qt2Prfmn(k1),qt2Psfav (k1),qt2Psfmn (k1))
    allocate(qt2Trfav (k1),qt2Trfmn(k1),qt2Disfav(k1),qt2Disfmn(k1))
    allocate(qt2Sfav  (k1),qt2Sfmn (k1))
    allocate(qt2tendf (k1))
    allocate(qt2residf(k1))
    allocate(qt2bf    (k1))

    thl2fav    = 0.
    thl2Prfav  = 0.; thl2Prfmn = 0.; thl2Psfav  = 0.; thl2Psfmn  = 0.
    thl2Trfav  = 0.; thl2Trfmn = 0.; thl2Disfav = 0.; thl2Disfmn = 0.
    thl2Sfav   = 0.; thl2Sfmn  = 0.
    thl2tendf  = 0.
    thl2residf = 0.
    thl2bf     = 0.

    qt2fav    = 0.
    qt2Prfav  = 0.; qt2Prfmn  = 0.; qt2Psfav  = 0.; qt2Psfmn  = 0.
    qt2Trfav  = 0.; qt2Trfmn  = 0.; qt2Disfav = 0.; qt2Disfmn = 0. 
    qt2Sfav   = 0.; qt2Sfmn   = 0.
    qt2tendf  = 0.
    qt2residf = 0.
    qt2bf     = 0.

    if(myid==0 .and. .not. lwarmstart) then
       open (ifoutput,file='varbudget.'//cexpnr,status='replace')
       close (ifoutput)
    end if
    if (lnetcdf) then
      idtav      = idtav_prof
      itimeav    = itimeav_prof
      tnext      = idtav+btime
      tnextwrite = itimeav+btime
      nsamples   = itimeav/idtav
     if (myid==0) then
        call ncinfo(ncname( 1,:),'thl2tendf','Tendency of thl variance','K^2/s','tt')
        call ncinfo(ncname( 2,:),'thl2Pr','Resolved production of thl variance','K^2/s','tt')
        call ncinfo(ncname( 3,:),'thl2Ps','SFS production of thl variance','K^2/s','tt')
        call ncinfo(ncname( 4,:),'thl2Tr','Resolved transport of thl variance','K^2/s','tt')
        call ncinfo(ncname( 5,:),'thl2D','Dissipation of thl variance','K^2/s','tt')
        call ncinfo(ncname( 6,:),'thl2S','Source of thl variance','K^2/s','tt')
        call ncinfo(ncname( 7,:),'thl2Res','Residual of thl budget','K^2/s','tt')
        call ncinfo(ncname( 8,:),'qt2tendf','Tendency of qt variance','kg^2/kg^2/s','tt')
        call ncinfo(ncname( 9,:),'qt2Pr','Resolved production of qt variance','kg^2/kg^2/s','tt')
        call ncinfo(ncname(10,:),'qt2Ps','SFS production of qt variance','kg^2/kg^2/s','tt')
        call ncinfo(ncname(11,:),'qt2Tr','Resolved transport of qt variance','kg^2/kg^2/s','tt')
        call ncinfo(ncname(12,:),'qt2D','Dissipation of qt variance','kg^2/kg^2/s','tt')
        call ncinfo(ncname(13,:),'qt2S','Source of qt variance','kg^2/kg^2/s','tt')
        call ncinfo(ncname(14,:),'qt2Res','Residual of qt budget','kg^2/kg^2/s','tt')
        call define_nc( ncid_prof, NVar, ncname)
     end if

   end if

  end subroutine initvarbudget

  subroutine varbudget

    use modglobal, only : rk3step,timee,dt_lim
    implicit none
    if (.not. lvarbudget) return
    if (rk3step/=3) return

    if(timee<tnext .and. timee<tnextwrite) then
      dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))
      return
    end if
    if (timee>=tnext) then
      tnext = tnext+idtav
      call do_varbudget
    end if
    if (timee>=tnextwrite) then
      tnextwrite = tnextwrite+itimeav
      call writevarbudget 
    end if
    dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))

  end subroutine varbudget


  subroutine do_varbudget

    use modfields     , only : qt0,thl0
    use modsubgriddata, only : ekh
    use modsurfdata   , only : thlflux,qtflux
    use modmicrodata  , only : qtpmcr
    use modraddata    , only : thlprad
    use modglobal     , only : i1,j1,ih,jh,k1,cp,ijtot, &
                               iadv_thl,iadv_qt
    use modmpi        , only : slabsum
  
    implicit none

  ! Set av values to zero
    thl2fav   = 0.
    thl2Prfav = 0.; thl2Psfav  = 0.
    thl2Trfav = 0.; thl2Disfav = 0.
    thl2Sfav  = 0.;

    qt2fav   = 0.
    qt2Prfav = 0.; qt2Psfav  = 0.
    qt2Trfav = 0.; qt2Disfav = 0.
    qt2Sfav  = 0.

! ---------------------------------------------------------
!  1. COMPUTE VARIANCE BUDGETS
! ---------------------------------------------------------

    call varproduction (thl0,thlflux,thlprad,iadv_thl, &
            thl2fav, &
            thl2Prfav,thl2Psfav,thl2Trfav,thl2Disfav,thl2Sfav)

    call varproduction (qt0,qtflux,qtpmcr,iadv_qt, &
           qt2fav, &
           qt2Prfav,qt2Psfav,qt2Trfav,qt2Disfav,qt2Sfav)

! ---------------------------------------------------------
!  2. ADD SLAB AVERAGES TO TIME MEANS
! ---------------------------------------------------------
    
    thl2Prfmn  = thl2Prfmn  + thl2Prfav
    thl2Psfmn  = thl2Psfmn  + thl2Psfav
    thl2Trfmn  = thl2Trfmn  + thl2Trfav
    thl2Disfmn = thl2Disfmn + thl2Disfav
    thl2Sfmn   = thl2Sfmn   + thl2Sfav

    qt2Prfmn  = qt2Prfmn  + qt2Prfav
    qt2Psfmn  = qt2Psfmn  + qt2Psfav
    qt2Trfmn  = qt2Trfmn  + qt2Trfav
    qt2Disfmn = qt2Disfmn + qt2Disfav
    qt2Sfmn   = qt2Sfmn   + qt2Sfav

  end subroutine do_varbudget

  subroutine varproduction (varxf,varxflux,src,iadv_var, &
                varx2fav, &
                resprodf,subprodf,restranf,disf,srcf)

! -----------------------------------------------------
!   compute variance budget
!   as input a (full level) variable is required
! -----------------------------------------------------
! TODO:  Enable covariance calculation
! TODO:  Speed up by not calculating zero derivatives
! FIXME: SG production term is still wrong

    use modglobal,      only : i1,i2,ih,j1,j2,jh,k1,kmax,     &
                               dzf,dzh,ijtot,dx2i,dy2i,cu,cv, &
                               iadv_cd2,iadv_5th,iadv_52,     &
                               iadv_cd6,iadv_62,iadv_kappa
    use modsubgriddata, only : ekh
    use modsubgrid,     only : diffc
    use modfields,      only : u0,v0,w0,u0av,v0av
    use modmpi,         only : comm3d,my_real,mpi_sum,mpierr, &
                               slabsum
    use mpi,             only : mpi_allreduce
    use advec_2nd,      only : advecc_2nd
    use advec_52,       only : advecc_52
    use advec_5th,      only : advecc_5th
    use advec_62,       only : advecc_62
    use advec_6th,      only : advecc_6th
    use advec_hybrid,   only : advecc_hybrid
    use advec_hybrid_f, only : advecc_hybrid_f
    use advec_kappa,    only : advecc_kappa
    use advec_upw,      only : advecc_upw

    implicit none

    integer i,j,k,im,ip,jm,jp,km,kp       !counter variables
    integer iadv_var

!    ----------- input variables
    real &
        varxf      (2-ih:i1+ih,2-jh:j1+jh,k1), &    !input variable x at full level &
        varxflux   (i2,j2),                    &    !surface flux of varx
        src        (2-ih:i1+ih,2-jh:j1+jh,k1), &    !source of variable x at full level &
        resprodf   (k1),                       &    !resolved production of variance, full level &
        subprodf   (k1),                       &    !'gradient' production, at full level &
        restranf   (k1),                       &    !resolved transport of (co-)variance at full level &
        disf       (k1),                       &    !'real' dissipation , at full level
        srcf       (k1)                             !source term interaction

!    ----------- function variables
    real &
        varxfdev   (2-ih:i1+ih,2-jh:j1+jh,k1), &    !fluctuation of varxf &
        varxfmn    (2-ih:i1+ih,2-jh:j1+jh,k1), &    !3D field of slab-averaged varxf &
        u0_dev     (2-ih:i1+ih,2-jh:j1+jh,k1), &    !fluctuation of u &
        v0_dev     (2-ih:i1+ih,2-jh:j1+jh,k1), &    !fluctuation of v &
        w0_dev     (2-ih:i1+ih,2-jh:j1+jh,k1), &    !fluctuation of w &
        u0_stor    (2-ih:i1+ih,2-jh:j1+jh,k1), &    !container for u &
        v0_stor    (2-ih:i1+ih,2-jh:j1+jh,k1), &    !container for v &
        w0_stor    (2-ih:i1+ih,2-jh:j1+jh,k1), &    !container for w &
        ekh_stor   (2-ih:i1+ih,2-jh:j1+jh,k1), &    !container for ekh &
        term       (2-ih:i1+ih,2-jh:j1+jh,k1), &    !output for adv/diff routines &
        dumfield   (2-ih:i1+ih,2-jh:j1+jh,k1), &    !dummy field to construct terms &
        varxfluxmn (i2,j2),                    &    !2D field of averaged varx surface flux &
        w0av       (k1),                       &    !mean vertical velocity
        ekhav      (k1),                       &    !mean ekh
        term_av    (k1),                       &    !mean term
        varxfav    (k1),                       &    !mean of varxf at full level &
        varx2fav   (k1)                             !variance of varxf at full level &

!    ----------- local (processor) variables
    real varx2favl   (k1)  !variance of varxf at full level &

    real varxfluxavl,varxfluxav

    ! Set to zero at start
    varxfdev    = 0.; varxfmn = 0.
    varxfluxmn  = 0.; varxfluxav = 0.
    w0av        = 0.
    ekhav       = 0.
    term_av     = 0.
    varxfav     = 0.
    varx2fav    = 0.
    resprodf    = 0.
    subprodf    = 0.
    restranf    = 0.
    disf        = 0.
    srcf        = 0.

    varx2favl   = 0.
    varxfluxavl = 0.

 !-------------------------------------------------------------
 ! Store the velocity fields 
 !-------------------------------------------------------------
    u0_stor = u0
    v0_stor = v0
    w0_stor = w0

 !-------------------------------------------------------------
 !      compute slab-averaged and deviation values of varx
 !-------------------------------------------------------------
 
    call slabsum(varxfav ,1,k1,varxf ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    varxfav = varxfav/ijtot

    ! Necessary to have 3D slab-averaged field to calculate production terms
    do k=1,k1
       varxfmn(:,:,k) = varxfav(k)
    enddo
 
    do k=1,k1
       varxfdev(:,:,k) = varxf(:,:,k) - varxfav(k)
    enddo

 !-------------------------------------------------------------
 !      compute mean and deviation values of wind fields
 !-------------------------------------------------------------

    call slabsum(w0av  ,1,k1,w0  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    w0av = w0av / ijtot

    do k=1,k1
       u0_dev(:,:,k) = u0(:,:,k) - (u0av(k) - cu)
       v0_dev(:,:,k) = v0(:,:,k) - (v0av(k) - cv)
       w0_dev(:,:,k) = w0(:,:,k) - w0av(k)
    enddo

 !-------------------------------------------------------------
 !      compute variance at full level (used for storing tendency)
 !-------------------------------------------------------------
    
    do j=2,j1
    do i=2,i1
    do k=1,k1
      varx2favl (k) = varx2favl (k) + varxfdev(i,j,k)**2
    end do
    end do
    end do

    call MPI_ALLREDUCE(varx2favl, varx2fav, k1, MY_REAL, &
                       MPI_SUM, comm3d,mpierr)
    varx2fav = varx2fav / ijtot

  !-------------------------------------------------------------
  !      compute variance production term resprodf
  !      resprodf = 2<varxf' * d/dxj(uj' <varxf>)
  !-------------------------------------------------------------

    ! Need to do this to advect with deviation fields
    u0 = u0_dev
    v0 = v0_dev
    w0 = w0_dev

    term = 0.

    select case(iadv_var)
      case(iadv_cd2)
        call advecc_2nd(varxfmn,term)
      case(iadv_5th)
        call advecc_5th(varxfmn,term)
      case(iadv_52)
        call advecc_52(varxfmn,term)
      case(iadv_cd6)
        call advecc_6th(varxfmn,term)
      case(iadv_62)
        call advecc_62(varxfmn,term)
      case(iadv_kappa)
        call advecc_kappa(varxfmn,term)
      case default
          stop "Unknown advection scheme "
    end select

    ! Reset fields
    u0 = u0_stor
    v0 = v0_stor
    w0 = w0_stor

    ! FIXME modstress applies cyclic BCs here

    dumfield = 2.*varxfdev*term
    call slabsum(resprodf ,1,k1,dumfield ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    resprodf = resprodf/ijtot

  !-------------------------------------------------------------
  !      compute variance transport term restranf:
  !      restranf = 2*varxfdev*d/dxj(uj*varxfdev)
  !      FIXME : Currently advects with u0,v0,w0 -> could be u0dev,v0dev,w0dev?
  !-------------------------------------------------------------

    term = 0.

    select case(iadv_var)
      case(iadv_cd2)
        call advecc_2nd(varxfdev,term)
      case(iadv_5th)
        call advecc_5th(varxfdev,term)
      case(iadv_52)
        call advecc_52(varxfdev,term)
      case(iadv_cd6)
        call advecc_6th(varxfdev,term)
      case(iadv_62)
        call advecc_62(varxfdev,term)
      case(iadv_kappa)
        call advecc_kappa(varxfdev,term)
      case default
          stop "Unknown advection scheme "
    end select

    ! FIXME modstress applies cyclic BCs here

    dumfield = 2.*varxfdev*term

    call slabsum(restranf ,1,k1,dumfield ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    restranf = restranf/ijtot

   !-------------------------------------------------------------
   !      Compute 'real' dissipation disf:
   !      disf = 2 <a' d/dxj (K d/dxj a')>
   !      where
   !      d/dxj (K d/dxj a') = d/dxj (K d/dxj a) - <d/dxj (K d/dxj a)>
   !-------------------------------------------------------------
    
    term = 0.

    call diffc(varxf,term,varxflux)

    call slabsum(term_av ,1,k1,term ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    term_av = term_av / ijtot

    do k=1,k1
       term(:,:,k) = term(:,:,k) - term_av(k)
    enddo

    dumfield = 2.*varxfdev*term

    call slabsum(disf ,1,k1,dumfield ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    disf = disf/ijtot
    
   !-------------------------------------------------------------
   !      Compute subgrid gradient production subprodf:
   !      subprodf = 2 <a' d/dxj (K' d/dxj <a>)>
   !-------------------------------------------------------------

   ! --- Calculate helpers
   ! Mean surface flux
    do i=2,i1
    do j=2,j1
      varxfluxavl = varxfluxavl + varxflux(i,j)
    end do
    end do
    call MPI_ALLREDUCE(varxfluxavl, varxfluxav, 1,    MY_REAL, &
                           MPI_SUM, comm3d,mpierr)
    varxfluxav = varxfluxav/ijtot

    varxfluxmn(:,:) = varxfluxav

    ! ekh'
    ekh_stor = ekh

    call slabsum(ekhav ,1,k1,ekh ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    ekhav = ekhav/ijtot

    do k=1,k1
       ekh(:,:,k) = ekh(:,:,k) - ekhav(k)
    enddo

    ! --- Calculate term
    term = 0.

    call diffc(varxfmn,term,varxfluxmn)

    ekh = ekh_stor

    ! FIXME modstress applies cyclic BCs here

    dumfield = 2.*varxfdev*term

    call slabsum(subprodf ,1,k1,dumfield ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    subprodf = subprodf/ijtot

   !-------------------------------------------------------------
   !      Compute source term srcf:
   !      srcf = 2 <a' S'>
   !-------------------------------------------------------------

    term_av = 0.
    call slabsum(term_av ,1,k1,src ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    term_av = term_av/ijtot

    term = 0.
    do k=1,k1
       term(:,:,k) = src(:,:,k) - term_av(k)
    enddo

    dumfield = 2.*varxfdev*term

    call slabsum(srcf ,1,k1,dumfield ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    srcf = srcf/ijtot

  return

  end subroutine varproduction

  subroutine writevarbudget

    use modglobal,  only : kmax,k1,zh,zf,rtimee,cexpnr,ifoutput
    use modmpi,     only : myid
    use modstat_nc, only : lnetcdf, writestat_nc
    use modgenstat, only : ncid_prof=>ncid,nrec_prof=>nrec
    implicit none

    real,dimension(k1,nvar) :: vars
    integer nsecs, nhrs, nminut,k
    nsecs   = nint(rtimee)
    nhrs    = int(nsecs/3600)
    nminut  = int(nsecs/60)-nhrs*60
    nsecs   = mod(nsecs,60)

! ---------------------------------------------------------
!    compute time-averaged fields
! ---------------------------------------------------------

    thl2tendf  = (thl2fav - thl2bf) / timeav
    thl2Prfmn  = thl2Prfmn          / nsamples
    thl2Psfmn  = thl2Psfmn          / nsamples
    thl2Trfmn  = thl2Trfmn          / nsamples
    thl2Disfmn = thl2Disfmn         / nsamples
    thl2Sfmn   = thl2Sfmn           / nsamples

    qt2tendf   = (qt2fav - qt2bf)   / timeav
    qt2Prfmn   = qt2Prfmn           / nsamples
    qt2Psfmn   = qt2Psfmn           / nsamples
    qt2Trfmn   = qt2Trfmn           / nsamples
    qt2Disfmn  = qt2Disfmn          / nsamples
    qt2Sfmn    = qt2Sfmn            / nsamples

    thl2bf = thl2fav
    qt2bf  = qt2fav

! ---------------------------------------------------------
!     compute residual
! ---------------------------------------------------------

    do k=1,k1
      thl2residf(k) = thl2tendf(k) - thl2Prfmn (k) - thl2Psfmn(k) &
                    - thl2Trfmn(k) - thl2Disfmn(k) - thl2Sfmn (k)

      qt2residf(k) = qt2tendf(k) - qt2Prfmn (k) - qt2Psfmn(k) &
                   - qt2Trfmn(k) - qt2Disfmn(k) - qt2Sfmn (k) 
    end do

 !-------------------------------------------------------------
 !     write output
 !-------------------------------------------------------------

    if(myid==0)then

      open (ifoutput,file='varbudget.'//cexpnr,position='append')

      call writebudget ( &
                nsecs,nhrs,nminut,timeav,'thl', &
                thl2tendf,thl2Prfmn ,thl2Psfmn, &
                thl2Trfmn,thl2Disfmn,thl2Sfmn , &
                thl2residf                    , &
                k1,zh,zf,ifoutput)

      call writebudget ( &
                nsecs,nhrs,nminut,timeav,'qt_', &
                qt2tendf,qt2Prfmn ,qt2Psfmn   , &
                qt2Trfmn,qt2Disfmn,qt2Sfmn    , &
                qt2residf                     , &
                k1,zh,zf,ifoutput)

      if (lnetcdf) then
        vars(:, 1)=thl2tendf
        vars(:, 2)=thl2Prfmn
        vars(:, 3)=thl2Psfmn
        vars(:, 4)=thl2Trfmn
        vars(:, 5)=thl2Disfmn
        vars(:, 6)=thl2Sfmn
        vars(:, 7)=thl2residf
        vars(:, 8)=qt2tendf
        vars(:, 9)=qt2Prfmn
        vars(:,10)=qt2Psfmn
        vars(:,11)=qt2Trfmn
        vars(:,12)=qt2Disfmn
        vars(:,13)=qt2Sfmn
        vars(:,14)=qt2residf
        call writestat_nc(ncid_prof,nvar,ncname,vars(1:kmax,:),nrec_prof,kmax)
      endif

    endif   !endif loop write myid=0

!        reset variables

    thl2Prfmn  = 0.
    thl2Psfmn  = 0.
    thl2Trfmn  = 0.
    thl2Disfmn = 0.
    thl2Sfmn   = 0.
    thl2residf = 0.

    qt2Prfmn  = 0.
    qt2Psfmn  = 0.
    qt2Trfmn  = 0.
    qt2Disfmn = 0.
    qt2Sfmn   = 0.
    qt2residf = 0.

  end subroutine writevarbudget

  subroutine exitvarbudget

    implicit none

    if(.not.(lvarbudget)) return

    deallocate(thl2fav   ,thl2Prfav ,thl2Prfmn,thl2Psfav,thl2Psfmn)
    deallocate(thl2Trfav ,thl2Trfmn ,thl2Sfav ,thl2Sfmn)
    deallocate(thl2Disfav,thl2Disfmn,thl2tendf,thl2residf,thl2bf)

    deallocate(qt2fav   ,qt2Prfav ,qt2Prfmn,qt2Psfav,qt2Psfmn)
    deallocate(qt2Trfav ,qt2Trfmn ,qt2Sfav ,qt2Sfmn)
    deallocate(qt2Disfav,qt2Disfmn,qt2tendf,qt2residf,qt2bf)

  end subroutine exitvarbudget


!---------------------------------------------------------------------------

  subroutine writebudget ( &
       sec,hr,mn,timeav,var, &
       x2tendf,x2Prfmn,x2Psfmn, &
       x2Trfmn,x2Disfmn,x2Sfmn, &
       x2residfmn, &
       ke,zh,zf,fnr)

  implicit none

  integer sec,hr,mn
  integer ke
  integer fnr
  real timeav
  character(3) var
  character(6) varheader

  real &
       x2tendf(ke), x2Prfmn (ke), x2Psfmn(ke), &
       x2Trfmn(ke), x2Disfmn(ke), x2Sfmn (ke), &
       x2residfmn(ke)

  real zh(ke),zf(ke)

  integer k

  if (var=='thl') then
     varheader = 'THETAL'
  else if (var=='qt_') then
     varheader = '- QTOT'
  else if (var(1:2)=='sv') then
     varheader = '-- SV '
     write (varheader(6:6),'(a)') var(3:3)
  endif

  write(fnr,'(A,/F5.0,A,I4,A,I2,A,I2,A)') &
     '--------------------------------------------------------'      &
     ,(timeav),'--- AVERAGING TIMESTEP --- '      &
     ,hr,':',mn,':',sec      &
     ,'   HRS:MIN:SEC AFTER INITIALIZATION '

  write (fnr,'(a/3a/a/2a,17a)') &
       '--------------------------------------------------------', &
       ' ----------------- ' , &
        varheader, ' VARIANCE BUDGET -----------------', &
       '--------------------------------------------------------', &
       ' lev heighth heightf tendencyf Prod_resf Prod_subf' , &
       ' Tran_resf Dissipationf Sourcef Residualf'

  do k=1,ke-1
    write (fnr,'(i3,2f8.2,14e15.7)') k, &
             zh(k), zf(k), &
             x2tendf(k),  &
             x2Prfmn(k),  &
             x2Psfmn(k),  &
             x2Trfmn(k),  &
             x2Disfmn(k), &
             x2Sfmn(k),   &
             x2residfmn(k)
  enddo

  return
  end subroutine writebudget


end module modvarbudget


