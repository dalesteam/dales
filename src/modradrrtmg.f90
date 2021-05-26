module modradrrtmg
  use modraddata
  implicit none

  private
  public :: radrrtmg

contains

  subroutine radrrtmg
    use modglobal,     only : cp,rlv,dzf,&
                              imax,jmax,kmax,i1,j1,k1,&
                              kind_rb,SHR_KIND_R4,boltz
    use modmpi,        only : myid
    use modfields,     only : initial_presh,initial_presf,rhof,exnf,thl0
    use modsurfdata ,  only : tskin
    use rrtmg_lw_init, only : rrtmg_lw_ini
    use rrtmg_lw_rad,  only : rrtmg_lw
    use shr_orb_mod,   only : shr_orb_params
    use rrtmg_sw_init, only : rrtmg_sw_ini
    use rrtmg_sw_rad,  only : rrtmg_sw
    implicit none

    integer                :: npatch    ! Sounding levels above domain
    integer                :: i,j,k,ierr(3)
    logical                :: sunUp
    real(SHR_KIND_R4),save ::  eccen, & ! Earth's eccentricity factor (unitless) (typically 0 to 0.1)
                               obliq, & ! Earth's obliquity angle (deg) (-90 to +90) (typically 22-26)
                               mvelp, & ! Earth's moving vernal equinox at perhelion (deg)(0 to 360.0)
                               !
                               ! Orbital information after processed by orbit_params
                               !
                               obliqr, &  ! Earth's obliquity in radians
                               lambm0, &  ! Mean longitude of perihelion at the vernal equinox (radians)
                               mvelpp     ! Earth's moving vernal equinox longitude
                                          ! of perihelion plus pi (radians)

    real                   :: thlpld,thlplu,thlpsd,thlpsu
    real(KIND=kind_rb)     :: cpdair
!    real(KIND=kind_rb),allocatable,dimension(:,:) :: cloudFrac, &
!                                                     liquidRe,  &
                                                     !iceRe,     &
                                                     !LWP_slice, &
                                                     !IWP_slice



    if(.not.isReadSounding) then
      call readSounding(initial_presh(k1)/100.,npatch_start,npatch_end)

      if(npatch_end.ne.npatch_start) then
        npatch = npatch_end - npatch_start + 1
      else
        if(myid==0) write(*,*) 'No sounding levels above the LES domain, check sounding input file'
        stop 'ERROR: No valid radiation sounding found (modradrrtmg.f90)'
      end if

      nzrad  = kmax + npatch    !old notation
      kradmax = nzrad   !a la kmax, k1
      krad1   = nzrad + 1
      krad2   = nzrad + 2

      isReadSounding = .true.
    end if

    ! Allocate the sizes of the slices once nzrad is known
    !
    !  nzrad = kmax + npatch, nzrad+1
    !
    if(.not.isAllocated_RadInputsOutputs) then
      allocate(layerP      (imax,krad1),       &
               layerT      (imax,krad1),       &
               h2ovmr      (imax,krad1),       &
               o3vmr       (imax,krad1),       &
               co2vmr      (imax,krad1),       &
               ch4vmr      (imax,krad1),       &
               n2ovmr      (imax,krad1),       &
               o2vmr       (imax,krad1),       &
               cfc11vmr    (imax,krad1),       &
               cfc12vmr    (imax,krad1),       &
               cfc22vmr    (imax,krad1),       &
               ccl4vmr     (imax,krad1),       &
               cloudFrac   (imax,krad1),       &
               liquidRe    (imax,krad1),       &
               iceRe       (imax,krad1),       &
!
               LWP_slice   (imax,krad1),       &
               IWP_slice   (imax,krad1),       &
               presh_input      (krad1),       &
                 STAT=ierr(1))

      allocate(tabs_slice  (imax,kradmax),         &
               qv_slice    (imax,kradmax),         &
               qcl_slice   (imax,kradmax),         &
               qci_slice   (imax,kradmax),         &
               o3_slice    (imax,kradmax),         &
               rho_slice   (imax,kradmax),       &
               tg_slice    (imax),               &
               presf_input      (kradmax),         &
!
               interfaceP     (imax,krad2),    &
               interfaceT     (imax,krad2),    &
!
               lwUp_slice     (imax,krad2),    &
               lwDown_slice   (imax,krad2),    &
               lwUpCS_slice   (imax,krad2),    &
               lwDownCS_slice (imax,krad2),    &
               swUp_slice     (imax,krad2),    &
               swDown_slice   (imax,krad2),    &
               swDownDir_slice(imax,krad2),    &
               swDownDif_slice(imax,krad2),    &
               swUpCS_slice   (imax,krad2),    &
               swDownCS_slice (imax,krad2),    &
               lwHR_slice     (imax,krad2),    &
               lwHRCS_slice   (imax,krad2),    &
               swHR_slice     (imax,krad2),    &
               swHRCS_slice   (imax,krad2),    &
!
                 STAT=ierr(2))
      allocate(solarZenithAngleCos(imax), asdir(imax), asdif(imax), aldir(imax),       &
                 aldif(imax),                                                          &
                 STAT=ierr(3))

      if(any(ierr(:)/=0)) then
        if(myid==0) write(*,*) 'Could not allocate input/output arrays in modradrrtmg'
        stop 'ERROR: Radiation variables could not be allocated in modradrrtmg.f90'
      else
        isAllocated_RadInputsOutputs = .true.
      end if
    end if

    if(.not.isReadTraceProfiles) then
      ! Patch sounding profile pressures above domain pressures (convert to hPa!)
      presf_input(1:kmax)   = initial_presf(1:kmax)  /100.
      presh_input(1:k1)     = initial_presh(1:k1)/100.

      if(npatch>0) then
        presf_input(k1  :kradmax) = psnd(npatch_start:npatch_end)
        presh_input(k1+1:kradmax) = 0.5*( psnd(npatch_start:npatch_end-1) &
                                      + psnd(npatch_start+1:npatch_end) )
        presh_input(krad1)      = max( 0.5*psnd(npatch_end),            &
                                      1.5*psnd(npatch_end) - 0.5*psnd(npatch_end-1) )
      end if
      call readTraceProfs

      if(myid==0) write(*,*) 'Trace gas profile have been read'
      isReadTraceProfiles = .true.
    end if

    cpdair = cp
    if(.not.isInitializedRrtmg) then
      if (rad_longw) then
        call rrtmg_lw_ini(cpdair)
      end if
      if (rad_shortw) then
        call shr_orb_params(iyear,eccen,obliq,mvelp,obliqr,lambm0,mvelpp,.false.)
        if (myid==0) write(*,*) 'orb_params = ',eccen,obliq,mvelp,obliqr,lambm0,mvelpp
        call rrtmg_sw_ini(cpdair)
      end if
      isInitializedRrtmg = .true.
    end if

! +=+=+=+=+=+=+=+= End of reading and initialization stage =+=+=++=+=+=+=+=+=+=+=++ !

    ! zero the RRTMG output arrays here, in case rrtmg_lw or rrtmg_sw is not executed
    swUp_slice = 0
    swDown_slice = 0
    swDownDir_slice = 0
    swDownDif_slice = 0
    swUpCS_slice = 0
    swDownCS_slice = 0
    lwUp_slice = 0
    lwDown_slice = 0
    lwUpCS_slice = 0
    lwDownCS_slice = 0

   ! Loop over the slices in the model, in the y direction
    do j=2,j1
      call setupSlicesFromProfiles &
           ( j, npatch_start, &                                           !input
           LWP_slice, IWP_slice, cloudFrac, liquidRe, iceRe )             !output

      if (rad_longw) then
        call rrtmg_lw &
             ( tg_slice, cloudFrac, IWP_slice, LWP_slice, iceRe, liquidRe )!input
        !if(myid==0) write(*,*) 'after call to rrtmg_lw'
      end if
      if (rad_shortw) then
         call setupSW(sunUp)
         if (sunUp) then
           call rrtmg_sw &
                ( tg_slice, cloudFrac, IWP_slice, LWP_slice, iceRe, liquidRe )
         end if
      end if

      lwu(2:i1,j,1:k1) =  lwUp_slice  (1:imax,1:k1)
      lwd(2:i1,j,1:k1) = -lwDown_slice(1:imax,1:k1)
      if (.not. rad_longw) then !we get LW at surface identically to how it is done in sunray subroutine 
        do i=2,i1
          lwd(i,j,1) =  -0.8 * boltz * thl0(i,j,1) ** 4.
          lwu(i,j,1) =  1.0 * boltz * tskin(i,j) ** 4.
        end do
      end if

      swu(2:i1,j,1:k1) =  swUp_slice  (1:imax,1:k1)
      swd(2:i1,j,1:k1) = -swDown_slice(1:imax,1:k1)

      swdir(2:i1,j,1:k1) = -swDownDir_slice(1:imax,1:k1)
      swdif(2:i1,j,1:k1) = -swDownDif_slice(1:imax,1:k1)
      lwc  (2:i1,j,1:k1) =  LWP_slice      (1:imax,1:k1)
 
      lwuca(2:i1,j,1:k1) =  lwUpCS_slice  (1:imax,1:k1)
      lwdca(2:i1,j,1:k1) = -lwDownCS_slice(1:imax,1:k1)
      swuca(2:i1,j,1:k1) =  swUpCS_slice  (1:imax,1:k1)
      swdca(2:i1,j,1:k1) = -swDownCS_slice(1:imax,1:k1)

      SW_up_TOA (2:i1,j) =  swUp_slice  (1:imax,krad2)
      SW_dn_TOA (2:i1,j) = -swDown_slice(1:imax,krad2)
      LW_up_TOA (2:i1,j) =  lwUp_slice  (1:imax,krad2)
      LW_dn_TOA (2:i1,j) = -lwDown_slice(1:imax,krad2)

      SW_up_ca_TOA (2:i1,j) =  swUpCS_slice  (1:imax,krad2)
      SW_dn_ca_TOA (2:i1,j) = -swDownCS_slice(1:imax,krad2)
      LW_up_ca_TOA (2:i1,j) =  lwUpCS_slice  (1:imax,krad2)
      LW_dn_ca_TOA (2:i1,j) = -lwDownCS_slice(1:imax,krad2)

    end do ! Large loop over j=2,j1

    do k=1,kmax
      do j=2,j1
        do i=2,i1
          thlpld          = -(lwd(i,j,k+1)-lwd(i,j,k))
          thlplu          = -(lwu(i,j,k+1)-lwu(i,j,k))
          thlpsd          = -(swd(i,j,k+1)-swd(i,j,k))
          thlpsu          = -(swu(i,j,k+1)-swu(i,j,k))

          !thlprad(i,j,k)  = thlprad(i,j,k) + (thlpld+thlplu+thlpsu+thlpsd)/(rhof(k)*cp*exnf(k)*dzf(k))
          thlprad(i,j,k)  = thlprad(i,j,k)-(lwd(i,j,k+1)-lwd(i,j,k)+lwu(i,j,k+1)-lwu(i,j,k)+swd(i,j,k+1)-swd(i,j,k)+swu(i,j,k+1)-swu(i,j,k)) &
                              /(rhof(k)*cp*exnf(k)*dzf(k))
        end do
      end do
    end do

    !if(myid==0) write(*,*) 'RadiationDone'
!    stop 'FINISHED radrrtmg!!'
  end subroutine radrrtmg

! ==============================================================================;
! ==============================================================================;

  subroutine readSounding(ptop_model,npatch_start,npatch_end)
    use modglobal, only     : cexpnr
    use modmpi, only        : myid
    use netcdf
    implicit none

    real,intent(in)        :: ptop_model
    integer,intent(out)    :: npatch_start,npatch_end        ! the level#s of the sounding above the model

    real,allocatable,dimension(:,:) :: psnd_in,qsnd_in,tsnd_in,o3snd_in
    integer                :: ncid,dimIDp, dimIDt, varID,sts ! netcdf id; netcdf function output
    integer                :: ntime,nlev                     ! dimensions in the netcdf file

    character (len=19)     :: SoundingFileName
    character (len=nf90_max_name) :: tmpName
    integer                :: k,ierr

    SoundingFileName = 'backrad.inp.'//cexpnr//'.nc'

    sts        = nf90_open(trim(SoundingFileName),nf90_nowrite,ncid)
    if (sts.ne.nf90_noerr) stop 'ERROR: Sounding file not found!'

    ! get number of pressure levels
    sts        = nf90_inq_dimid(ncid,"lev",dimIDp)
    sts        = nf90_inquire_dimension(ncid, dimIDp, tmpName, nlev)

    ! get number of time levels
    sts        = nf90_inq_dimid(ncid,"time",dimIDt)
    if(sts.eq.nf90_NoErr) then
      sts        = nf90_inquire_dimension(ncid, dimIDt, tmpName, ntime)
    end if
    if(sts.ne.nf90_NoErr) then
      ntime      = 1
    end if

    if(nlev.gt.nzsnd) then
      write(*,*) 'ERROR in modradrrtmg.f90'
      write(*,*) '***** number of levels in ', TRIM(SoundingFileName)
      write(*,*) '***** exceeds nzsnd.  Reset nzsnd in modraddata.f90'
      stop 'ERROR in readSounding'
    end if

    allocate(psnd_in(nlev,ntime), tsnd_in(nlev,ntime), qsnd_in(nlev,ntime), &
             o3snd_in(nlev,ntime), &
             psnd(nlev), tsnd(nlev), qsnd(nlev), o3snd(nlev), &
             STAT=ierr)
    if(ierr.ne.0) then
      write(*,*) 'ERROR in modradrrtmg.f90'
      write(*,*) '***** Could not allocate arrays to read in sounding'
      stop 'ERROR in readSounding'
    end if

    ! get pressure levels (in Pascal)
    sts          = nf90_inq_varid(ncid,"lev",varID)
    sts          = nf90_get_var(ncid, varID, psnd_in)

    ! get temperature and moisture (K and kg/kg, respectively)
    sts          = nf90_inq_varid(ncid,"T",varID)
    sts          = nf90_get_var(ncid, varID, tsnd_in)
    sts          = nf90_inq_varid(ncid,"q",varID)
    sts          = nf90_get_var(ncid, varID, qsnd_in)
    if (usero3) then
      sts          = nf90_inq_varid(ncid,"o3",varID)
      if (sts/=nf90_NoErr) then
        usero3 = .false.
        if (myid==0) then
          write(*,*) 'WARNING: could not obtain ozon profile from file.'
          write(*,*) 'Using reference profile instead!'
          write(*,*) 'usero3 has now been set to .false.'
        end if
      else
        sts          = nf90_get_var(ncid, varID, o3snd_in)
      end if
    end if

    psnd(1:nlev) = 0. ; tsnd(1:nlev) = 0. ; qsnd(1:nlev) = 0. ; o3snd(1:nlev) = 0.
    ! reverse order from top-down to bottom-up.
    ! psnd(1:nlev) = psnd_in(nlev:1:-1,1)/100. ! convert from Pa to hPa
    ! tsnd(1:nlev) = tsnd_in(nlev:1:-1,1)
    ! qsnd(1:nlev) = qsnd_in(nlev:1:-1,1)
    ! Variables are already bottom-up in the Netcdf file for some reason .. JvdD

    ! No time dependent reference profiles are allowed (yet), first timestep taken
    do k=1,nlev
      psnd(k) = psnd_in(k,1) / 100
      tsnd(k) = tsnd_in(k,1)
      qsnd(k) = qsnd_in(k,1)
      o3snd(k) = o3snd_in(k,1)
    enddo

    ! find whether we need this sounding to be patched on top of model.
    npatch_start = nlev+1
    npatch_end = nlev+1

    do k = 1,nlev
      !if(psnd(k).lt.ptop_model - 10.) then
      if(psnd(k).lt.ptop_model) then
        ! start patch ~10hPa above model top.
        npatch_start = k
        EXIT
      end if
    end do

    if(npatch_start.le.nlev) then
      npatch_end = nlev
    end if

    if(myid == 0) then
      write(*,*)
      write(*,*) 'Background sounding, p (mb), T (K), q (kg/kg)'
      do k = 1,nlev
        if(k.eq.npatch_start) write(*,*) '**** Start Patch Here *****'
        write(*,998) k,psnd(k), tsnd(k), qsnd(k)
998     format(i4,f8.3,f8.3,e12.4)
        if(k.eq.npatch_end) write(*,*) '**** End Patch Here *****'
      end do
      if(npatch_start.gt.nlev) write(*,*) '**** No patching required -- model top is deeper than sounding ****'
    end if

    deallocate(psnd_in, tsnd_in, qsnd_in, o3snd_in, STAT=ierr)
    if(ierr.ne.0) then
      write(*,*) 'ERROR in modradrrtmg.f90'
      write(*,*) '***** Could not allocate arrays to read in sounding'
      stop 'ERROR in readSounding'
    end if

  end subroutine readSounding

! ==============================================================================;
! ==============================================================================;

  subroutine readTraceProfs        ! original tracesini subroutine in rad_driver
    use modglobal, only : kind_rb, kind_im, &
                          grav
    use modmpi, only    : myid
    use rrlw_ncpar, only: getAbsorberIndex
    use netcdf
    implicit none

    integer, parameter :: nTraceGases = 9
!    integer :: nz_tracegases
    integer :: sts(9),ncid,dimIDp,dimIDab,varID,np,nab
    integer :: ab ! Absorber Index
    integer :: m,ks,k,ierr

    real    :: plow, pupp, pmid, wgtlow, wgtupp
!    real    :: godp
    real(kind=kind_rb),allocatable,dimension(:)    :: pMLS  ! Sounding pressure
    real(kind=kind_rb),allocatable,dimension(:,:)  :: trace, trace_in
    real(kind=kind_rb) :: tmppresf(krad1),           &
                          tmppresh(krad2),           &
                          tmpTrace(krad1),           &    ! Temporary trace gas profile
                          trpath(krad2,nTraceGases), &
                          godp(krad1)

    character(len=nf90_max_name) :: tmpName
    character(len=5),dimension(nTraceGases),parameter :: traceGasNameOrder = (/ &
       'O3   ','CO2  ','CH4  ','N2O  ','O2   ','CFC11','CFC12','CFC22','CCL4 '/)

    if(.not.isAllocated_TraceGases) then
      ! allocate trace gas arrays.  These have an extra level for the
      !   mean mass-weighted trace gas concentration in the overlying atmosphere.
      !
!      nz_tracegases = krad1  ! add one level to compute trace gas levels to TOA
      allocate(o3(krad1), co2(krad1), ch4(krad1), &
           n2o(krad1), o2(krad1), cfc11(krad1), &
           cfc12(krad1), cfc22(krad1), ccl4(krad1), &
           STAT=ierr)
      if(ierr.ne.0) then
        write(*,*) 'ERROR: could not allocate trace gas arrays in tracesini'
        stop 'ERROR in readTraceProfs'
      else
        isAllocated_TraceGases=.true.
      end if
    end if

    ! Read profiles from rrtmg data file.
    sts(:)  = nf90_noerr
    sts(1)  = nf90_open('rrtmg_lw.nc',nf90_nowrite,ncid)

    sts(2)  = nf90_inq_dimid(ncid,"Pressure",dimIDp)
    sts(3)  = nf90_inquire_dimension(ncid, dimIDp, tmpName, np)

    sts(4)  = nf90_inq_dimid(ncid,"Absorber",dimIDab)
    sts(5)  = nf90_inquire_dimension(ncid, dimIDab,tmpName, nab)

    if (any(sts.ne.nf90_noerr)) then
      write(*,*) 'ERROR: input file either not found or incorrectly formatted'
      stop 'rrtmg_lw.nc input file either not found or incorrectly formatted'
    end if

    ! allocate local variables and set their value to zero
    allocate(pMLS(np), trace(nTraceGases,np), trace_in(nab,np), STAT=ierr)
    if(ierr.ne.0) then
      write(*,*) 'ERROR: could not declare arrays in tracesini'
      stop 'ERROR: Could not allocate (in readTraceProfs, modradrrtmg.f90)'
    end if
    pMLS=0. ; trace=0. ; trace_in=0.

    sts(6)  = nf90_inq_varid(ncid,"Pressure",varID)
    sts(7)  = nf90_get_var(ncid, varID, pMLS)              !Get sounding pressure levels

    sts(8)  = nf90_inq_varid(ncid,"AbsorberAmountMLS",varID)
    sts(9)  = nf90_get_var(ncid, varID, trace_in)          !Get amounts of absorber gases

    do m = 1,nTraceGases
      call getAbsorberIndex(TRIM(traceGasNameOrder(m)),ab) !trace gas order (in), ab (out)
      trace(m,1:np) = trace_in(ab,1:np)                    !Make sure the profiles are stored
      where (trace(m,:)>2.)                                !the right order
        trace(m,:) = 0.
      end where
    end do

    if(maxval(abs(sts(:))).ne.nf90_noerr) then
      write(*,*) 'Error in reading trace gas sounding from RRTMG_data/rrtmg_lw.nc'
      stop 'ERROR: tracegas data could not be read from rrtmg_lw.nc'
    end if

    ! An extra layer is added to the existing pressure profiles (for better results?)
    tmppresf(1:kradmax)   = presf_input(1:kradmax)
    tmppresh(1:krad1) = presh_input(1:krad1)

    tmppresf(krad1)   = 0.5*presh_input(krad1)
    tmppresh(krad2)   = min(1.e-4_kind_rb,0.25*tmppresf(krad1))

    ! trace gas paths at surface are zero.
    trpath(1,:) = 0.

    ! Loop over 2nd lowest level to top of atmosphere (including added level)
    do k = 2,krad2
      ! start with trace path at interface below.
      trpath(k,:) = trpath(k-1,:)

      ! if pressure greater than trace gas sounding, assume concentration at bottom.
      ! >extrapolate
      if (tmppresh(k-1).gt.pMLS(1)) then
        trpath(k,:) = trpath(k,:) &
             + (tmppresh(k-1) - max(tmppresh(k),pMLS(1)))/grav & ! dp/g
             *trace(:,1)                                 ! *tr
      end if

      ! loop over sounding levels
      do ks = 2,np
        ! limit pMLS(ks:ks-1) so that they are within the model level
        !  tmppresi(k-1:k).
        plow = min(tmppresh(k-1),max(tmppresh(k),pMLS(ks-1)))
        pupp = min(tmppresh(k-1),max(tmppresh(k),pMLS(ks))  )

        ! Interpolation of trace gas concentration
        if(plow.gt.pupp) then
          pmid = 0.5*(plow+pupp)

          wgtlow = (pmid-pMLS(ks))/(pMLS(ks-1)-pMLS(ks))
          wgtupp = (pMLS(ks-1)-pmid)/(pMLS(ks-1)-pMLS(ks))
          ! include this level of the sounding in the trace gas path
          trpath(k,:) = trpath(k,:) &
               + (plow - pupp)/grav*(wgtlow*trace(:,ks-1)+wgtupp*trace(:,ks)) ! dp/g*tr
        end if
      end do ! loop over sounding levels, ks

      ! if pressure is less than trace gas sounding, assume concentration at top
      ! >extrapolate
      if (tmppresh(k).lt.pMLS(np)) then
        trpath(k,:) = trpath(k,:) &
             + (min(tmppresh(k-1),pMLS(np)) - tmppresh(k))/grav &
             *trace(:,np)
      end if

    end do ! loop over levels, k

    do m = 1,nTraceGases
      ! Faster than original loop ??
      godp(:)  = grav / (tmppresh(1:krad1) - tmppresh(2:krad2))
      tmpTrace = ( trpath(2:krad2,m) - trpath(1:krad1,m) ) * godp(:)

!      do k = 1,nzm+1
!        godp = ggr/(tmppresi(k) - tmppresi(k+1))
!        tmpTrace(k) = (trpath(k+1,m) - trpath(k,m))*godp
!      end do

      if    (TRIM(traceGasNameOrder(m))=='O3')    then
        o3(:) = tmpTrace(:)
      elseif(TRIM(traceGasNameOrder(m))=='CO2')   then
        co2(:) = REAL(co2factor,KIND=kind_rb)*tmpTrace(:)
      elseif(TRIM(traceGasNameOrder(m))=='CH4')   then
        ch4(:) = tmpTrace(:)  !cstep ch4(:) = ch4factor for RCE
      elseif(TRIM(traceGasNameOrder(m))=='N2O')   then
        n2o(:) = tmpTrace(:)  !cstep n2o(:) = n2ofactor for RCE
      elseif(TRIM(traceGasNameOrder(m))=='O2')    then
        o2(:) = tmpTrace(:)
      elseif(TRIM(traceGasNameOrder(m))=='CFC11') then
        cfc11(:) = tmpTrace(:)
      elseif(TRIM(traceGasNameOrder(m))=='CFC12') then
        cfc12(:) = tmpTrace(:)
      elseif(TRIM(traceGasNameOrder(m))=='CFC22') then
        cfc22(:) = tmpTrace(:)
      elseif(TRIM(traceGasNameOrder(m))=='CCL4')  then
        ccl4(:) = tmpTrace(:)
      end if
    end do !loop over trace gases, m

    if(myid==0)then
      write(*,*) 'RRTMG rrtmg_lw.nc trace gas profile: number of levels=',np
      write(*,*) 'gas traces vertical profiles (ppmv *10^-6):'
      write(*,*) 'p, hPa', ('       ',traceGasNameOrder(m),m=1,nTraceGases)
      do k=1,krad1
        write(*,*) tmppresf(k),o3(k),co2(k),ch4(k),n2o(k),o2(k), &
             cfc11(k),cfc12(k), cfc22(k),ccl4(k)
      end do
    end if

    deallocate(pMLS, trace, trace_in, STAT=ierr)
    if(ierr.ne.0) then
      write(*,*) 'ERROR: could not deallocate arrays in tracesini'
      stop 'ERROR: Could not deallocate (in readTraceProfs, modradrrtmg.f90)'
    end if
  end subroutine readTraceProfs

! ==============================================================================;
! ==============================================================================;

  subroutine setupSlicesFromProfiles(j,npatch_start, &
           LWP_slice,IWP_slice,cloudFrac,liquidRe,iceRe)
  !=============================================================================!
  ! This subroutine sets up 2D (xz) slices of different variables:              !
  ! tabs,qv,qcl,qci(=0),tg,layerP,interfaceP,layerT,interfaceT,LWP,IWP(=0),     !
  ! cloudFrac,liquidRe,iceRe(=0)                                                !
  ! for the full height of the atmosphere (up to nzrad+2, where the upper level !
  ! is an extrapolation, for more accuracy)                                     !
  !                                                                             !
  ! The variables are stored in modraddata.f90 and allocated in the subroutine  !
  ! from which setupSlicesFromProfiles is called (radrrtmg)                     !
  ! JvdDussen, 24-6-2010                                                        !
  ! ============================================================================!

      use modglobal, only: imax,jmax,kmax,i1,k1,grav,kind_rb,rlv,cp,Rd,pref0
      use modfields, only: thl0,ql0,qt0,exnf
      use modsurfdata, only: tskin,ps
      use modmicrodata, only : Nc_0,sig_g
      use modmpi, only: myid

      implicit none

      integer,intent(in) :: j,npatch_start
      real(KIND=kind_rb),intent(out) ::    LWP_slice(imax,krad1), &
                                           IWP_slice(imax,krad1), &
                                           cloudFrac(imax,krad1), &
                                           liquidRe (imax,krad1), &
                                           iceRe    (imax,krad1)
      integer :: i,k,ksounding,im
      real (KIND=kind_rb) :: exners
      real(KIND=kind_rb) :: layerMass(imax,krad1)
      !real(KIND=kind_rb),dimension(imax,kmax)     :: tabs         ! Absolute temperature
      real(KIND=kind_rb),dimension(imax,jmax)     :: sstxy        ! sea surface temperature
      real   (SHR_KIND_R4), parameter :: pi = 3.14159265358979
      real , parameter :: rho_liq = 1000.

      real :: reff_factor
      real :: ilratio
      real :: tempC  !temperature in celsius
      real :: IWC0 ,B_function !cstep needed for ice effective radius following Eqs. (14) and (35) from Wyser 1998

      IWC0 = 50e-3  !kg/m3, Wyser 1998 Eq. 14 (he gives 50 g/m3)
      reff_factor = 1e6*(3. /(4.*pi*Nc_0*rho_liq) )**(1./3.) * exp(log(sig_g)**2 )

      ! Compute absolute temperature and water contents (without border points)

     ! tabs(:,:) = 0.;
      sstxy(:,:) = 0.
      tabs_slice(:,:) = 0.; qv_slice(:,:) = 0.; qcl_slice(:,:) = 0.; qci_slice(:,:) = 0.;
      rho_slice(:,:) = 0.

      exners = (ps/pref0) ** (rd/cp)

      do i=2,i1
      do k=1,kmax
          im = i-1
          tabs_slice(im,k) = thl0(i,j,k) * exnf(k) &
                          + (rlv / cp) * ql0(i,j,k)
      enddo
      enddo

      !tabs(:,:)     = thl0(2:i1,j,2:k1) * spread(exnf(2:k1),dim=1,ncopies=imax) &
      !                  + (rlv / cp) * ql0(2:i1,j,2:k1)

      do i=2,i1
        im=i-1

        !tg_slice  (im)   = sst
        tg_slice  (im)   = tskin(i,j) * exners  ! Note: tskin = thlskin...

        do k=1,kmax
           qv_slice  (im,k) = max(qt0(i,j,k) - ql0(i,j,k),1e-18) !avoid RRTMG reading negative initial values 
           qcl_slice (im,k) = ql0(i,j,k)
           qci_slice (im,k) = 0.
           o3_slice  (im,k) = o3snd(npatch_start) ! o3 constant below domain top (if usero3!)

           h2ovmr    (im,k) = mwdry/mwh2o * qv_slice(im, k)
!           h2ovmr    (im,k) = mwdry/mwh2o * (qv_slice(im,k)/(1-qv_slice(im,k)))
           layerT    (im,k) = tabs_slice(im,k)
           layerP    (im,k) = presf_input(k)
        enddo
      enddo

     ! Patch sounding on top (no qcl or qci above domain; hard coded)
      do i=1,imax
      ksounding=npatch_start
      do k=kmax+1,kradmax
         tabs_slice(i,k) =  tsnd(ksounding)
         qv_slice  (i,k) =  qsnd(ksounding)
         qcl_slice (i,k) = 0.
         qci_slice (i,k) = 0.
         rho_slice (i,k) = 100*presf_input(k)/(Rd*tabs_slice(i,k)) !cstep factor 100 because pressure in hPa
         ksounding=ksounding+1
      enddo
      enddo

      ! o3 profile provided by user (user03=true) or reference prof from RRTMG
      if (usero3) then
        do i=1,imax
        ksounding=npatch_start
        do k=kmax+1,kradmax
           o3_slice(i,k) = o3snd(ksounding)
           ksounding=ksounding+1
        enddo
        do k=1,kradmax
          o3vmr  (i, k)   = mwdry/mwo3 * o3_slice(i, k)
        enddo
        o3vmr  (i, krad1)   = o3vmr(i,kradmax)
        enddo
      else
        do i=1,imax
        do k=1,krad1
            o3vmr   (i, k) = o3(k)
        enddo
        enddo
      end if

      do i=1,imax
        do k=kmax+1,kradmax

           !h2ovmr  (i, k)    = mwdry/mwh2o * (qv_slice(i,k)/(1-qv_slice(i,k)))
           h2ovmr  (i, k)    = mwdry/mwh2o * qv_slice(i,k)
           layerP(i,k)       = presf_input (k)
           layerT(i,k)       = tabs_slice(i,k)
        enddo
        ! Properly set boundary conditions
!        h2ovmr  (i, krad1)   = mwdry/mwh2o * qv_slice(i,kradmax)
        h2ovmr  (i, krad1)   = h2ovmr(i,kradmax)
        layerP  (i, krad1)   = 0.5*presh_input(krad1)
        layerT  (i, krad1)   = 2.*tabs_slice(i, kradmax) - tabs_slice(i, kradmax-1)
      enddo

      do i=1,imax
        do k=1,krad1
          co2vmr  (i, k) = co2(k)
          ch4vmr  (i, k) = ch4(k)
          n2ovmr  (i, k) = n2o(k)
          o2vmr   (i, k) = o2(k)
          cfc11vmr(i, k) = cfc11(k)
          cfc12vmr(i, k) = cfc12(k)
          cfc22vmr(i, k) = cfc22(k)
          ccl4vmr (i, k) = ccl4(k)

          interfaceP(i,k ) =   presh_input(k)
        enddo

        interfaceP(i, krad2)  = min( 1.e-4_kind_rb , 0.25*layerP(1,krad1) )
        do k=2,krad1
           interfaceT(i, k) = (layerT(i,k-1) + layerT(i, k)) / 2.
        enddo
        interfaceT(i, krad2) = 2.*layerT(i, krad1) - interfaceT(i, krad1)
        interfaceT(i, 1)  = tg_slice(i)
      enddo

      do i=1,imax
        do k=1,kradmax
          layerMass(i,k) = 100.*( interfaceP(i,k) - interfaceP(i,k+1) ) / grav  !of full level
          LWP_slice(i,k) = qcl_slice(i,k)*layerMass(i,k)*1e3
          IWP_slice(i,k) = qci_slice(i,k)*layerMass(i,k)*1e3
          qci_slice(i,k) = qci_slice(i,k)*rho_slice(i,k)   !cstep, qci is now in kg/m3 needed for ice effective radius
        enddo
        layerMass(i,krad1) = 100.*( interfaceP(i,krad1) - interfaceP(i,krad2) ) / grav
        LWP_slice(i,krad1) = 0.
        IWP_slice(i,krad1) = 0.
      enddo

      cloudFrac(:,:) = 0.
      liquidRe (:,:) = 0.
      iceRe    (:,:) = 0.

      do i=1,imax
        do k=1,kradmax
          cloudFrac(i,k) = 0.
          liquidRe (i,k) = 0.
          iceRe    (i,k) = 0.
          if (LWP_slice(i,k).gt.0.) then
!            liquidRe(i,k) = 14.   ! cstep temporary solution, ocean value computeRe_Liquid(real(layerT),merge(0.,1.,ocean))
            cloudFrac(i,k) = 1.

            !cstep liquidRe(i, k) = 1.e6*( 3.*( 1.e-3*LWP_slice(i,k)/layerMass(i,k) ) &
            !cstep                  /(4.*pi*Nc_0*rho_liq) )**(1./3.) * exp(log(sig_g)**2 )
            !cstep: equation above contains function of many constants, are now absorbed in reff_factor

            liquidRe(i, k) = reff_factor  * qcl_slice(i,k)**(1./3.)

!  cstep, 1e6*Nc_0 assumes Nc_0 in cm^-3             /(4.*pi*1.e6*Nc_0*rho_liq) )**(1./3.) * exp(log(sigmag)**2 )
!  cstep              write (6,*) 'reff = ',i,k,qcl_slice(i,k),liquidRe(i,k)

            if (liquidRe(i,k).lt.2.5) then
                liquidRe(i,k) = 2.5
            endif

            if (liquidRe(i,k).gt.60.) then
                liquidRe(i,k) = 60.
            endif
          endif

          if (IWP_slice(i,k).gt.0) then
             cloudFrac(i,k) = 1.
             !cstep Ou Liou: tempC = layerT(i,k)--tmelt
             !cstep Ou Liou  iceRe(i,k) = 326.3 + 12.42 * tempC + 0.197 * tempC**2 + 0.0012 * tempC**3  !cstep : Ou Liou 1995
             B_function =  -2 + 0.001 *(273.-layerT(i,k))**1.5 * log10(qci_slice(i,k)/IWC0) !Eq. 14 Wyser 1998
             iceRe (i,k) = 377.4 + 203.3 * B_function + 37.91 * B_function**2 + 2.3696 * B_function**3 !micrometer, Wyser 1998, Eq. 35
             cloudFrac(i,k) = 1.
             B_function =  -2 + 0.001 *(273.-layerT(i,k))**1.5 * log10(qci_slice(i,k)/IWC0)
             iceRe (i,k) = 377.4 + 203.3 * B_function + 37.91 * B_function**2 + 2.3696 * B_function**3 !micrometer, Wyser 1998

             if (iceRe(i,k).lt.5.) then
                iceRe(i,k) = 5.
             endif

             if (iceRe(i,k).gt.140.) then
                iceRe(i,k) = 140.
             endif
          endif
        enddo
      enddo

  end subroutine setupSlicesFromProfiles

! ==============================================================================;
! ==============================================================================;

  subroutine setupSW(sunUp)

    use modglobal,   only : xday,xlat,xlon,imax,xtime,rtimee
    use shr_orb_mod, only : shr_orb_decl
    use modmpi,      only : myid
    use modsurfdata, only : albedoav

    implicit none

    logical,intent(out) :: sunUp
    real                :: dayForSW

    if(doperpetual) then
        !====================!
        !   To be inserted   !
        !====================!
    else

      if(doseasons) then
        ! The diurnal cycle of insolation will vary
        ! according to time of year of the current day.
        dayForSW = xday + (xtime + rtimee/3600) / 24
      else
        ! The diurnal cycle of insolation from the calendar
        ! day on which the simulation starts (day0) will be
        ! repeated throughout the simulation.
        ! dayForSW = float(floor(day0)) + xday - float(floor(xday))
      end if

      call shr_orb_decl( dayForSW )                      ! Saves some orbital values to modraddata
      solarZenithAngleCos(:) =  &
           zenith(xtime*3600 + rtimee, xday, xlat, xlon) ! Used function in modraddata
!      solarZenithAngleCos(:) =  0.707106781               ! cos 45gr
!      solarZenithAngleCos(:) = 0.087155742747658           ! cos 85gr
       !solarZenithAngleCos(:) = 0.615661475  !cos 52 gr
       !solarZenithAngleCos(:) = 0.356   !cos 69.144 gr
      !solarZenithAngleCos(:) = 1.
       ! solarZenithAngleCos(:) = mu0_cgils

    end if

    sunUp = .false.
    ! if all values in solarZenithAngleCos are >= its smallest positive, non-zero element
    if (all(solarZenithAngleCos(:) >= tiny(solarZenithAngleCos))) then
      sunUp = .true.
      if (lCnstAlbedo) then
        aldir = albedoav
        asdir = albedoav
        aldif = albedoav        ! Specification of the diffuse albedo is also important for the
        asdif = albedoav        ! total surface albedo
      else
        call albedo             ! calculate albedo for the solarZenithAngleCos
      end if

    end if

  end subroutine setupSW

! ==============================================================================;
! ==============================================================================;

  elemental real function computeRe_Liquid(temperature, landfrac, icefrac, snowh) &
     result(rel)
    real,           intent(in) :: temperature, landfrac
    real, optional, intent(in) :: icefrac, snowh  ! Snow depth over land, water equivalent (m)

    real, parameter ::  rliqland  =  8.0, & ! liquid drop size if over land
                        rliqocean = 14.0, & ! liquid drop size if over ocean
                        rliqice   = 14.0    ! liquid drop size if over sea ice

    ! jrm Reworked effective radius algorithm
    ! Start with temperature-dependent value appropriate for continental air
    rel = rliqland + (rliqocean - rliqland) * min(1.0, max(0.0, (tmelt - temperature) * 0.05))

    if(present(snowh)) & ! Modify for snow depth over land
      rel = rel + (rliqocean - rel) * min(1.0, max(0.0, snowh*10.))

    ! Ramp between polluted value over land to clean value over ocean.
    rel = rel + (rliqocean-rel) * min(1.0, max(0.0, 1.0 - landfrac))

    if(present(icefrac)) & ! Ramp between the resultant value and a sea ice value in the presence of ice.
      rel = rel + (rliqice-rel) * min(1.0, max(0.0, icefrac))

  end function computeRe_Liquid

! ==============================================================================;
! ==============================================================================;
  subroutine albedo
    !-----------------------------------------------------------------------
    ! Computes surface albedos over ocean
    ! and the surface (added by Marat Khairoutdinov)

    ! Two spectral surface albedos for direct (dir) and diffuse (dif)
    ! incident radiation are calculated. The spectral intervals are:
    !   s (shortwave)  = 0.2-0.7 micro-meters
    !   l (longwave)   = 0.7-5.0 micro-meters
    !
    ! Uses knowledge of surface type to specify albedo, as follows:
    !
    ! Ocean           Uses solar zenith angle to compute albedo for direct
    !                 radiation; diffuse radiation values constant; albedo
    !                 independent of spectral interval and other physical
    !                 factors such as ocean surface wind speed.
    !
    ! For more details , see Briegleb, Bruce P., 1992: Delta-Eddington
    ! Approximation for Solar Radiation in the NCAR Community Climate Model,
    ! Journal of Geophysical Research, Vol 97, D7, pp7603-7612).

    real, parameter :: adif = 0.06
    !-----------------------------------------------------------------------
    if (ocean) then
      !
      !
      ! Ice-free ocean albedos function of solar zenith angle only, and
      ! independent of spectral interval:
      !
      where(solarZenithAngleCos(:) <= 0)
        aldir(:) = 0; asdir(:) = 0; aldif(:) = 0; asdif(:) = 0
      elsewhere
        aldir(:) = ( .026 / (solarZenithAngleCos(:)**1.7 + .065)) + &
             (.15*(solarZenithAngleCos(:) - 0.10) * (solarZenithAngleCos(:) - 0.50) * &
             (solarZenithAngleCos(:) - 1.00) )
        asdir(:) = aldir(:)
        aldif(:) = adif
        asdif(:) = adif
      end where
    else ! land
      where(solarZenithAngleCos(:) <= 0)
        aldir(:) = 0; asdir(:) = 0; aldif(:) = 0; asdif(:) = 0
      elsewhere
        ! Albedos for land type I (Briegleb)
        asdir(:) = 1.4 * 0.06 / ( 1. + 0.8 * solarZenithAngleCos(:) )
        asdif(:) = 1.2 * 0.06
        aldir(:) = 1.4 * 0.24 / ( 1. + 0.8 * solarZenithAngleCos(:) )
        aldif(:) = 1.2 * 0.24
      end where
    endif

  end subroutine albedo

end module modradrrtmg
