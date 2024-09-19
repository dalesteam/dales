!> \file modfielddump.f90
!>
!!  Dumps 3D fields of several variables, lime fielddump but configurable per MPI process.
!!  In the namelist a region istart,iend,jstart,jend can be specified.
!!  The MPI processes that overlap with this region will output their 3D fields, like a fielddump.
!!  Note: more data than requested may be saved.
!!
!!  Namelist sample:
!!  &NAMLOCALFIELDDUMP
!!  llocalfielddump = .true.
!!  dtav            = 3600
!!  istart          = 100
!!  iend            = 120
!!  jstart          =  30
!!  jend            =  50
!!  /


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
module modlocalfielddump
  use modprecision, only: field_r
  use modglobal, only : longint, nsv

implicit none
private
PUBLIC :: initlocalfielddump, localfielddump,exitlocalfielddump
save
!NetCDF variables
  integer :: nvar = 7
  integer :: ncid,nrec = 0
  character(80) :: fname = 'localdump.xxx.xxx.xxx.nc'
  character(80),dimension(:,:), allocatable :: ncname
  character(80),dimension(1,4) :: tncname

  real    :: dtav, tmin, tmax
  integer(kind=longint) :: idtav,tnext,itmax,itmin
  integer :: klow,khigh,ncoarse=1
  logical :: llocalfielddump= .false. !< switch to enable the fielddump (on/off)
  logical :: lu = .true.         !< switch for saving the u field
  logical :: lv = .true.         !< switch for saving the v field
  logical :: lw = .true.         !< switch for saving the w field
  logical :: lqt = .true.        !< switch for saving the qt field
  logical :: lql = .true.        !< switch for saving the ql field
  logical :: lthl = .true.       !< switch for saving the thl field
  logical :: lbuoy = .true.      !< switch for saving the buoy field
  logical :: lcli = .false.       !< switch for saving the cli field
  logical :: lclw = .false.       !< switch for saving the clw field
  logical :: lta = .false.        !< switch for saving the ta field
  logical :: lplw = .false.       !< switch for saving the plw field
  logical :: lpli = .false.       !< switch for saving the pli field
  logical :: lhus = .false.       !< switch for saving the hus field
  logical :: lhur = .false.       !< switch for saving the hur field
  logical :: ltntr = .false.      !< switch for saving the tntr field
  logical :: ltntrs = .false.     !< switch for saving the tntrs field
  logical :: ltntrl = .false.     !< switch for saving the tntrl field
  logical :: lsv(100) = .true.   !< switches for saving the sv fields

  ! indices for the variables in the netCDF vars array
  integer :: ind, ind_u=-1, ind_v=-1, ind_w=-1, ind_qt=-1, ind_ql=-1, ind_thl=-1, ind_buoy=-1, ind_sv(100)=-1
  integer :: ind_cli=-1, ind_clw=-1, ind_ta=-1, ind_plw=-1, ind_pli=-1, ind_hus=-1, ind_hur=-1, ind_tntr=-1, ind_tntrs=-1, ind_tntrl=-1
  integer :: istart=-1, iend=-1, jstart=-1, jend=-1
  integer :: myistart, myiend, myjstart, myjend
  logical :: localfielddumpactive

contains
!> Initializing fielddump. Read out the namelist, initializing the variables
  subroutine initlocalfielddump
    use modmpi,   only :myid,comm3d,myidx,myidy &
                       , D_MPI_BCAST
    use modglobal,only :imax,jmax,kmax,i1,j1,cexpnr,ifnamopt,fname_options,dtmax,dtav_glob,kmax, ladaptive,dt_lim,btime,tres,&
         checknamelisterror, output_prefix
    use modstat_nc,only : lnetcdf,open_nc, define_nc,ncinfo,nctiminfo,writestat_dims_nc
    use modtracers, only : tracer_prop
    use modmicrodata, only : imicro, imicro_sice, imicro_sice2
    implicit none
    integer :: ierr, n
    character(3) :: csvname

    namelist/NAMLOCALFIELDDUMP/ &
         dtav,llocalfielddump,klow,khigh,ncoarse, tmin, tmax,&
         lu, lv, lw, lqt, lql, lthl, lbuoy, lcli, lclw, lta, lplw, lpli, lhus, lhur, ltntr, ltntrs, ltntrl, lsv,&
         istart, iend, jstart, jend

    dtav=dtav_glob
    klow=1
    khigh=kmax
    tmin = 0.
    tmax = 1e8
    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMLOCALFIELDDUMP,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMLOCALFIELDDUMP')
      write(6 ,NAMLOCALFIELDDUMP)
      close(ifnamopt)

      if ((lcli .or. lclw) .and. .not. (imicro ==  imicro_sice .or. imicro == imicro_sice2)) then
         write (*,*) "LOCALFIELDDUMP: cli and clw output works only with simpleice microphysics. Turning off."
         lcli = .false.
         lclw = .false.
      end if
    end if
    call D_MPI_BCAST(ncoarse         ,1,0,comm3d,ierr)
    call D_MPI_BCAST(klow            ,1,0,comm3d,ierr)
    call D_MPI_BCAST(khigh           ,1,0,comm3d,ierr)
    call D_MPI_BCAST(dtav            ,1,0,comm3d,ierr)
    call D_MPI_BCAST(tmin            ,1,0,comm3d,ierr)
    call D_MPI_BCAST(tmax            ,1,0,comm3d,ierr)
    call D_MPI_BCAST(llocalfielddump ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lu              ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lv              ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lw              ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lqt             ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lql             ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lthl            ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lbuoy           ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lcli            ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lclw            ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lta             ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lplw            ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lpli            ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lhus            ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lhur            ,1,0,comm3d,ierr)
    call D_MPI_BCAST(ltntr           ,1,0,comm3d,ierr)
    call D_MPI_BCAST(ltntrs          ,1,0,comm3d,ierr)
    call D_MPI_BCAST(ltntrl          ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lsv           ,100,0,comm3d,ierr)
    call D_MPI_BCAST(istart          ,1,0,comm3d,ierr)
    call D_MPI_BCAST(iend            ,1,0,comm3d,ierr)
    call D_MPI_BCAST(jstart          ,1,0,comm3d,ierr)
    call D_MPI_BCAST(jend            ,1,0,comm3d,ierr)

    idtav = dtav/tres
    itmin = tmin/tres
    itmax = tmax/tres

    tnext      = idtav   +btime
    if(.not.(llocalfielddump)) return
    dt_lim = min(dt_lim,tnext)

    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'dtav should be a integer multiple of dtmax'
    end if

    ! find the global start and end i,j indices of this domain
    myistart =    2 + myidx*imax
    myiend   =   i1 + myidx*imax
    myjstart =    2 + myidy*jmax
    myjend   =   j1 + myidy*jmax

    ! does this domain overlap with the local fielddump domain?
    localfielddumpactive = .true.
    if (istart > myiend .or. iend < myistart) localfielddumpactive = .false.
    if (jstart > myjend .or. jend < myjstart) localfielddumpactive = .false.

    if (.not. localfielddumpactive) return

    if (lnetcdf) then
      write(fname,'(A,i3.3,A,i3.3,A)') 'localdump.', myidx, '.', myidy, '.xxx.nc'
      fname(19:21) = cexpnr
      nvar = 17+nsv ! maximum number of variables
      allocate(ncname(nvar,4))

      call nctiminfo(tncname(1,:))
      ind = 1
      if (lu) then
         ind_u = ind
         ind = ind + 1
         call ncinfo(ncname(ind_u,:),'u','West-East velocity','m/s','mttt')
      end if
      if (lv) then
         ind_v = ind
         ind = ind + 1
         call ncinfo(ncname(ind_v,:),'v','South-North velocity','m/s','tmtt')
      end if
      if (lw) then
         ind_w = ind
         ind = ind + 1
         call ncinfo(ncname(ind_w,:),'w','Vertical velocity','m/s','ttmt')
      end if
      if (lqt) then
         ind_qt = ind
         ind = ind + 1
         call ncinfo(ncname(ind_qt,:),'qt','Total water specific humidity','kg/kg','tttt')
      end if
      if (lql) then
         ind_ql = ind
         ind = ind + 1
         call ncinfo(ncname(ind_ql,:),'ql','Liquid water specific humidity','kg/kg','tttt')
      end if
      if (lthl) then
         ind_thl = ind
         ind = ind + 1
         call ncinfo(ncname(ind_thl,:),'thl','Liquid water potential temperature','K','tttt')
      end if
      if (lbuoy) then
         ind_buoy = ind
         ind = ind + 1
         call ncinfo(ncname(ind_buoy,:),'buoy','Buoyancy','K','tttt')
      end if
      if (lcli) then
         ind_cli = ind
         ind = ind + 1
         call ncinfo(ncname(ind_cli,   :),'cli','mass fraction of cloud ice','kg/kg','tttt') ! new
      end if
      if (lclw) then
         ind_clw = ind
         ind = ind + 1
         call ncinfo(ncname(ind_clw,   :),'clw','mass fraction of cloud liquid water','kg/kg','tttt') ! new
      end if
      if (lta) then
         ind_ta = ind
         ind = ind + 1
         call ncinfo(ncname(ind_ta,    :),'ta','air temperature ','K','tttt') ! new
      end if
      if (lplw) then
         ind_plw = ind
         ind = ind + 1
         call ncinfo(ncname(ind_plw,   :),'plw','mass fraction of precipitating liquid water','kg/kg','tttt') !new
      end if
      if (lpli) then
         ind_pli = ind
         ind = ind + 1
         call ncinfo(ncname(ind_pli,   :),'pli','mass fraction of precipitating ice','kg/kg','tttt')   !new
      end if
      if (lhus) then
         ind_hus = ind
         ind = ind + 1
         call ncinfo(ncname(ind_hus,   :),'hus','specific humidity','kg/kg','tttt')   !new
      end if
      if (lhur) then
         ind_hur = ind
         ind = ind + 1
         call ncinfo(ncname(ind_hur,   :),'hur','relative humidity','%','tttt')      !new
      end if
      if (ltntr) then
         ind_tntr = ind
         ind = ind + 1
         call ncinfo(ncname(ind_tntr,  :),'tntr','tendency of air temperature due to radiative heating','K/s','tttt')   !new
      end if
      if (ltntrs) then
         ind_tntrs = ind
         ind = ind + 1
         call ncinfo(ncname(ind_tntrs, :),'tntrs','tendency of air temperature due to shortwave radiative heating','K/s','tttt')  !new
      end if
      if (ltntrl) then
         ind_tntrl = ind
         ind = ind + 1
         call ncinfo(ncname(ind_tntrl, :),'tntrl','tendency of air temperature due to longwave radiative heating','K/s','tttt')   !new
      end if

      do n=1,nsv
        if (lsv(n)) then
           ind_sv(n) = ind
           ind = ind + 1
           write (csvname(1:3),'(i3.3)') n
           call ncinfo(ncname(ind_sv(n),:), tracer_prop(n)%tracname, tracer_prop(n)%traclong, tracer_prop(n)%unit, 'tttt')
        end if
     end do

      nvar = ind - 1 ! total number of fields actually in use

      call open_nc(trim(output_prefix)//fname,  ncid,nrec,n1=ceiling(1.0*imax/ncoarse),n2=ceiling(1.0*jmax/ncoarse),n3=khigh-klow+1)

      if (nrec==0) then
        call define_nc( ncid, 1, tncname)
        call writestat_dims_nc(ncid, ncoarse, klow)
     end if
     call define_nc( ncid, NVar, ncname(1:NVar,:))
     ! must slice ncname here because define_nc expects first dimension to be NVar
     ! and NVar may have decreased if some fields are not saved
    end if

  end subroutine initlocalfielddump
!> Do local fielddump. Saves the 3D fields if this MPI domain overlaps with the domain given in NAMOPTIONS.
!> if lnetcdf, write to netCDF (as float32) (the only option)
  subroutine localfielddump
    use modfields, only : u0,v0,w0,thl0,qt0,ql0,sv0,thv0h,thvh,tmp0,rhof,exnf,presf
    use modglobal, only : imax,i1,jmax,j1,rk3step,dzf, &
                          timee,dt_lim,rtimee,cp,tdn,tup
    use modstat_nc, only : lnetcdf, writestat_nc
    use modmicrodata, only : iqr, tuprsg, tdnrsg
    use modraddata, only   :lwu,lwd,swu,swd
    use modthermodynamics, only: qsat_tab
#if defined(_OPENACC)
    use modgpu, only: update_host
#endif
    implicit none

    real, allocatable :: vars(:,:,:,:)
    integer i,j,k,n,ii,jj
    integer :: writecounter = 1
    integer :: reclength


    if (.not. llocalfielddump) return
    if (.not. localfielddumpactive) return
    if (rk3step/=3) return

    if(timee<tnext) then
      dt_lim = min(dt_lim,tnext-timee)
      return
    end if

    tnext = tnext+idtav
    dt_lim = minval((/dt_lim,tnext-timee/))

    ! Only write fields if time is in the range (tmin, tmax)
    if (timee < itmin .or. timee > itmax) return

    !$acc update self(u0) if(lu) async
    !$acc update self(v0) if(lv) async
    !$acc update self(w0) if(lw) async
    !$acc update self(qt0) if(lqt) async
    !$acc update self(ql0) if(lql) async
    !$acc update self(thl0) if(lthl) async
    !$acc update self(sv0) if(any(lsv)) async
    !$acc update self(thv0h, thvh) if(lbuoy) async
    !$acc wait

    if (lnetcdf) allocate(vars(ceiling(1.0*imax/ncoarse),ceiling(1.0*jmax/ncoarse),khigh-klow+1,nvar))

    reclength = ceiling(1.0*imax/ncoarse)*ceiling(1.0*jmax/ncoarse)*(khigh-klow+1)*2


    if (lnetcdf  .and. lu) vars(:,:,:,ind_u) = u0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)


    if (lnetcdf .and. lv) vars(:,:,:,ind_v) = v0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)

    if (lnetcdf .and. lw) vars(:,:,:,ind_w) = w0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)

    if (lnetcdf .and. lqt) vars(:,:,:,ind_qt) = qt0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)

    if (lnetcdf .and. lql) vars(:,:,:,ind_ql) = ql0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)

    if (lnetcdf .and. lthl) vars(:,:,:,ind_thl) = thl0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)

    ! buoyancy
    if (lnetcdf .and. lbuoy) then
      vars(:,:,:,ind_buoy) = thv0h(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
      do k=klow,khigh
        vars(:,:,k-klow+1,ind_buoy) = vars(:,:,k-klow+1,ind_buoy) - thvh(k)
      end do
    end if

    if (lnetcdf .and. lcli) vars(:,:,:,ind_cli) = ql0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh) * &
         (1 - max(0._field_r,min(1._field_r,(tmp0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)-tdn)/(tup-tdn))))

    if (lnetcdf .and. lclw) vars(:,:,:,ind_clw) = ql0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh) * &
         max(0._field_r,min(1._field_r,(tmp0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)-tdn)/(tup-tdn)))

    if (lnetcdf .and. lta) vars(:,:,:,ind_ta) = tmp0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)

    ! liquid and ice precip
    ! assuming simpleice is used
    !rsgratio_= max(0._field_r,min(1._field_r,(tmp0(i,j,k)-tdnrsg)/(tuprsg-tdnrsg))) ! rain vs snow/graupel partitioning   rsg = 1 if t > tuprsg
    if (lnetcdf .and. lplw) vars(:,:,:,ind_plw) = sv0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh,iqr) * &
         max(0._field_r,min(1._field_r,(tmp0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)-tdnrsg)/(tuprsg-tdnrsg)))
    if (lnetcdf .and. lpli) vars(:,:,:,ind_pli) = sv0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh,iqr)  * &
         (1 - max(0._field_r,min(1._field_r,(tmp0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)-tdnrsg)/(tuprsg-tdnrsg))))

    if (lnetcdf .and. lhus) vars(:,:,:,ind_hus) = qt0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh) - ql0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)

    if (lnetcdf .and. lhur) then
       do k = klow, khigh
          do j = 1, ceiling(1.0*jmax/ncoarse)
             jj = 2 + (j-1) * ncoarse
             do i = 1, ceiling(1.0*imax/ncoarse)
                ii = 2 + (i-1) * ncoarse
                vars(i,j,k-klow+1,ind_hur) = 100 * (qt0(ii,jj,k) - ql0(ii,jj,k)) / &
                     qsat_tab(tmp0(ii,jj,k), presf(k))
             end do
          end do
       end do
    end if

    if (lnetcdf .and. ltntr) then
       do k = klow,khigh
          vars(:,:,k-klow+1,ind_tntr) = ( &    !! need to loop over k here
               -swd(2:i1:ncoarse,2:j1:ncoarse,k+1) &
               -swu(2:i1:ncoarse,2:j1:ncoarse,k+1) &
               +swd(2:i1:ncoarse,2:j1:ncoarse,k) &
               +swu(2:i1:ncoarse,2:j1:ncoarse,k) &
               -lwd(2:i1:ncoarse,2:j1:ncoarse,k+1) &
               -lwu(2:i1:ncoarse,2:j1:ncoarse,k+1) &
               +lwd(2:i1:ncoarse,2:j1:ncoarse,k) &
               +lwu(2:i1:ncoarse,2:j1:ncoarse,k) &
               )  / (rhof(k)*exnf(k)*cp*dzf(k))
       end do
    end if

    if (lnetcdf .and. ltntrs) then
       do k = klow,khigh
          vars(:,:,k-klow+1,ind_tntrs) = ( & !! need to loop over k here
               -swd(2:i1:ncoarse,2:j1:ncoarse,k+1) &
               -swu(2:i1:ncoarse,2:j1:ncoarse,k+1) &
               +swd(2:i1:ncoarse,2:j1:ncoarse,k) &
               +swu(2:i1:ncoarse,2:j1:ncoarse,k) &
               )  / (rhof(k)*exnf(k)*cp*dzf(k))
       end do
    end if

    if (lnetcdf .and. ltntrl) then
       do k = klow,khigh
          vars(:,:,k-klow+1,ind_tntrl) = ( & !! need to loop over k here
               -lwd(2:i1:ncoarse,2:j1:ncoarse,k+1) &
               -lwu(2:i1:ncoarse,2:j1:ncoarse,k+1) &
               +lwd(2:i1:ncoarse,2:j1:ncoarse,k) &
               +lwu(2:i1:ncoarse,2:j1:ncoarse,k) &
               )  / (rhof(k)*exnf(k)*cp*dzf(k))
       end do
    end if

    ! scalar variables
    if (lnetcdf) then
       do n=1,nsv
          if (lsv(n)) then
             vars(:,:,:,ind_sv(n)) = sv0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh,n)
          end if
       end do
    end if

    if(lnetcdf) then
      call writestat_nc(ncid,1,tncname,(/rtimee/),nrec,.true.)
      call writestat_nc(ncid,nvar,ncname,vars,nrec,ceiling(1.0*imax/ncoarse),ceiling(1.0*jmax/ncoarse),khigh-klow+1)
    end if


    writecounter=writecounter+1

    if (lnetcdf) deallocate(vars)

  end subroutine localfielddump
!> Clean up when leaving the run
  subroutine exitlocalfielddump
    use modstat_nc, only : exitstat_nc,lnetcdf
    implicit none

    if(llocalfielddump .and. localfielddumpactive .and. lnetcdf) call exitstat_nc(ncid)
  end subroutine exitlocalfielddump

end module modlocalfielddump
