!> \file modfielddump.f90
!!  Dumps 3D fields of several variables

!>
!!  Dumps 3D fields of several variables
!>
!!  Dumps 3D fields of several variables Written to wb*.myidx.myidy.expnr
!! If netcdf is true, this module leads the fielddump.myidx.myidy..expnr.nc output
!!  \author Thijs Heus,MPI-M
!!  \par Revision list
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
module modfielddump
  use fortran_support, only: &
#if FIELD_PRECISION==64
                             t_ptr_3d => t_ptr_3d_dp
#else
                             t_ptr_3d => t_ptr_3d_sp
#endif
  use modfields, only: u0, v0, w0, qt0, ql0, e120, thl0, tmp0, sv0, rhof, &
                       exnf, thv0h, thvh, presf
  use modglobal, only: j1, i1, dzf, cp, tdn, tup, nsv
  use modraddata, only: lwu, lwd, swu, swd
  use modsubgriddata, only: ekh0 => ekh, ekm0 => ekm
  use modtracers, only: get_tracer_index, tracer_prop
  use modthermodynamics, only: calc_qsat
  use modprecision, only: field_r, longint
  use modlogging, only: finish
  use modnetcdf_file_t, only: field_dump_file_t
  use modstat_nc_files, only: add_output_file, is_sampling_timestep

implicit none
character(len=*), parameter :: modname = 'modfielddump'
private

PUBLIC :: initfielddump, fielddump

save
!NetCDF variables

  type(field_dump_file_t) :: ofile    !< Output file object.
  integer                 :: ofile_id !< File ID in the list of output files.

  real    :: dtav, tmin, tmax
  integer(kind=longint) :: idtav,tnext,itmax,itmin
  integer :: klow,khigh,ncoarse=1
  logical :: lfielddump= .false. !< switch to enable the fielddump (on/off)
  logical :: ldiracc   = .false. !< switch for doing direct access writing (on/off)
  logical :: lbinary   = .false. !< switch for doing direct access writing (on/off)
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
  logical :: le12 = .false.       !< switch for saving the e12 field
  logical :: lekh = .false.       !< switch for saving the ekh field
  logical :: lekm = .false.       !< switch for saving the ekm field
  logical :: lsv(100) = .true.   !< switches for saving the sv fields

contains
!> Initializing fielddump. Read out the namelist, initializing the variables
  subroutine initfielddump
    use modmpi,   only :myid,comm3d,myidx,myidy,D_MPI_BCAST
    use modglobal,only :imax,jmax,kmax,cexpnr,ifnamopt,fname_options,dtmax,dtav_glob,kmax, ladaptive,dt_lim,btime,tres,&
         checknamelisterror, output_prefix
    use modstat_nc,only : lnetcdf,open_nc, define_nc,ncinfo,nctiminfo,writestat_dims_nc
    use modtracers, only : tracer_prop, get_tracer_index
    use modmicrodata, only : imicro, imicro_sice, imicro_sice2
    use fortran_support, only: nnml_output
    implicit none

    character(len=*), parameter :: routine = modname//'/initfielddump'

    integer :: ierr, n, iqr
    character(3) :: csvname

    namelist/NAMFIELDDUMP/ &
         dtav,lfielddump,ldiracc,lbinary,klow,khigh,ncoarse, tmin, tmax,&
         lu, lv, lw, lqt, lql, lthl, lbuoy, lcli, lclw, lta, lplw, lpli, lhus, lhur, ltntr, ltntrs, ltntrl, le12, lekh, lekm,  lsv

    dtav=dtav_glob
    klow=1
    khigh=kmax
    tmin = 0.
    tmax = 1e8
    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMFIELDDUMP,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMFIELDDUMP')
      write(nnml_output ,NAMFIELDDUMP)
      close(ifnamopt)

      if ((lcli .or. lclw) .and. .not. (imicro ==  imicro_sice .or. imicro == imicro_sice2)) then
         write (*,*) "FIELDDUMP: cli and clw output works only with simpleice microphysics. Turning off."
         lcli = .false.
         lclw = .false.
      end if
      if (lplw.or.lpli) then
        iqr = get_tracer_index("qr")
      endif
      if ((iqr == 0).and.((lplw.or.lpli))) then
        print *, "lplw or lpli are true but there is no qr tracer. Turning plw and pli output off."
        lplw = .false.
        lpli = .false. 
      endif
    end if
    call D_MPI_BCAST(ncoarse     ,1,0,comm3d,ierr)
    call D_MPI_BCAST(klow        ,1,0,comm3d,ierr)
    call D_MPI_BCAST(khigh       ,1,0,comm3d,ierr)
    call D_MPI_BCAST(dtav        ,1,0,comm3d,ierr)
    call D_MPI_BCAST(tmin        ,1,0,comm3d,ierr)
    call D_MPI_BCAST(tmax        ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lfielddump  ,1,0,comm3d,ierr)
    call D_MPI_BCAST(ldiracc     ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lbinary     ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lu          ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lv          ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lw          ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lqt         ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lql         ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lthl        ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lbuoy       ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lcli        ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lclw        ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lta         ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lplw        ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lpli        ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lhus        ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lhur        ,1,0,comm3d,ierr)
    call D_MPI_BCAST(ltntr       ,1,0,comm3d,ierr)
    call D_MPI_BCAST(ltntrs      ,1,0,comm3d,ierr)
    call D_MPI_BCAST(ltntrl      ,1,0,comm3d,ierr)
    call D_MPI_BCAST(le12        ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lekh        ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lekm        ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lsv       ,100,0,comm3d,ierr)

    idtav = int(dtav / tres, kind=kind(idtav))
    itmin = int(tmin / tres, kind=kind(itmin))
    itmax = int(tmax / tres, kind=kind(itmax))

    tnext      = idtav   +btime
    if(.not.(lfielddump)) return
    dt_lim = min(dt_lim,tnext)

    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      call finish(routine, 'dtav should be a integer multiple of dtmax')
    end if
      
    ofile = field_dump_file_t('fielddump2.nc', nz=kmax, ncoarse=ncoarse, &
                              klo=klow, khi=khigh)

    call add_output_file(ofile, dtav, ofile_id)

    if (lu) call ofile%add_var('u', 'West-East velocity', 'm/s', 'mttt')
    if (lv) call ofile%add_var('v', 'South-North velocity', 'm/s', 'tmtt')
    if (lw) call ofile%add_var('w', 'Vertical velocity', 'm/s', 'ttmt')
    if (lqt) call ofile%add_var('qt', 'Total water specific humidity', 'kg/kg', 'tttt')
    if (lql) call ofile%add_var('ql', 'Liquid water specific humidity', 'kg/kg', 'tttt')
    if (lthl) call ofile%add_var('thl', 'Liquid water potential temperature', 'K', 'tttt')
    if (lbuoy) call ofile%add_var('buoy', 'Buoyancy', 'K', 'tttt')
    if (lcli) call ofile%add_var('cli', 'mass fraction of cloud ice', 'kg/kg', 'tttt')
    if (lclw) call ofile%add_var('clw', 'mass fraction of cloud liquid water', 'kg/kg', 'tttt')
    if (lta) call ofile%add_var('ta', 'air temperature', 'K', 'tttt')
    if (lplw) call ofile%add_var('plw', 'mass fraction of precipitating liquid water', 'kg/kg', 'tttt')
    if (lpli) call ofile%add_var('pli', 'mass fraction of precipitating ice', 'kg/kg', 'tttt')
    if (lhus) call ofile%add_var('hus', 'specific humidity', 'kg/kg', 'tttt')
    if (lhur) call ofile%add_var('hur', 'relative humidity', 'kg/kg', 'tttt')
    if (ltntr) call ofile%add_var('tntr', 'tendency of air temperature due to radiative heating', 'K/s', 'tttt')
    if (ltntrs) call ofile%add_var('tntrs', 'tendency of air temperature due to shortwave radiative heating', 'K/s', 'tttt')
    if (ltntrl) call ofile%add_var('tntrl', 'tendency of air temperature due to longwave radiative heating', 'K/s', 'tttt')
    if (le12) call ofile%add_var('e12', 'square root of turbulent kinetic energy', 'm/s', 'tttt')
    if (lekh) call ofile%add_var('ekh', 'diffusion coefficient for heat and moisture', 'm2/s', 'tttt')
    if (lekm) call ofile%add_var('e12', 'diffusion coefficient for momentum', 'm2/s', 'tttt')

    do n = 1, nsv
      if (lsv(n)) then
        call ofile%add_var(tracer_prop(n)%tracname, tracer_prop(n)%traclong, tracer_prop(n)%unit, 'tttt')
      end if
    end do

  end subroutine initfielddump

  subroutine fielddump

    integer :: i, j, k, n, ii, jj, kk
    integer :: iqr

    real(field_r), pointer :: u(:,:,:)
    real(field_r), pointer :: v(:,:,:)
    real(field_r), pointer :: w(:,:,:)
    real(field_r), pointer :: qt(:,:,:)
    real(field_r), pointer :: ql(:,:,:)
    real(field_r), pointer :: e12(:,:,:)
    real(field_r), pointer :: ekm(:,:,:)
    real(field_r), pointer :: ekh(:,:,:)
    real(field_r), pointer :: thl(:,:,:)
    real(field_r), pointer :: ta(:,:,:)
    real(field_r), pointer :: buoy(:,:,:)
    real(field_r), pointer :: cli(:,:,:)
    real(field_r), pointer :: clw(:,:,:)
    real(field_r), pointer :: pli(:,:,:)
    real(field_r), pointer :: plw(:,:,:)
    real(field_r), pointer :: hus(:,:,:)
    real(field_r), pointer :: hur(:,:,:)
    real(field_r), pointer :: tntr(:,:,:)
    real(field_r), pointer :: tntrs(:,:,:)
    real(field_r), pointer :: tntrl(:,:,:)
    
    type(t_ptr_3d) :: sv(100)

    if (lfielddump .and. is_sampling_timestep(ofile_id)) then
      
      if (lu) then
        call ofile%get_pointer('u', u) 
        !$acc kernels default(present) async
        u(:,:,:) = u0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
        !$acc end kernels
      end if

      if (lv) then
        call ofile%get_pointer('v', v) 
        !$acc kernels default(present) async
        v(:,:,:) = v0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
        !$acc end kernels
      end if

      if (lw) then
        call ofile%get_pointer('w', w)
        !$acc kernels default(present) async
        w(:,:,:) = w0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
        !$acc end kernels
      end if

      if (lqt) then
        call ofile%get_pointer('qt', qt)
        !$acc kernels default(present) async
        qt(:,:,:) = qt0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
        !$acc end kernels
      end if

      if (lql) then
        call ofile%get_pointer('ql', ql)
        !$acc kernels default(present) async
        ql(:,:,:) = ql0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
        !$acc end kernels
      end if

      if (lthl) then
        call ofile%get_pointer('thl', thl)
        !$acc kernels default(present) async
        thl(:,:,:) = thl0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
        !$acc end kernels
      end if

      if (le12) then
        call ofile%get_pointer('e12', e12)
        !$acc kernels default(present) async
        e12(:,:,:) = e120(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
        !$acc end kernels
      end if

      if (lekh) then
        call ofile%get_pointer('ekh', ekh)
        !$acc kernels default(present) async
        ekh(:,:,:) = ekh0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
        !$acc end kernels
      end if

      if (lekm) then
        call ofile%get_pointer('ekm', ekm)
        !$acc kernels default(present) async
        ekm(:,:,:) = ekm0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
        !$acc end kernels
      end if

      if (lta) then
        call ofile%get_pointer('ta', ta)
        !$acc kernels default(present) async
        ta(:,:,:) = tmp0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
        !$acc end kernels
      end if

      do n = 1, nsv
        if (lsv(n)) then
          call ofile%get_pointer(tracer_prop(n)%tracname, sv(n)%p)
          !$acc kernels default(present) async
          sv(n)%p(:,:,:) = sv0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh,n)
          !$acc end kernels
        end if
      end do

      if (lbuoy) then
        call ofile%get_pointer('buoy', buoy)
        !$acc parallel loop collapse(3) default(present) async
        do k = klow, khigh
          do j = 2, j1, ncoarse
            do i = 2, i1, ncoarse
              kk = k - klow + 1
              jj = (j - 2) / ncoarse + 2
              ii = (i - 2) / ncoarse + 2
              buoy(ii,jj,kk) = thv0h(i,j,k) - thvh(k)
            end do
          end do
        end do
      end if

      if (lcli) then
        call ofile%get_pointer('cli', cli)
        !$acc parallel loop collapse(3) default(present) async 
        do k = klow, khigh
          do j = 2, j1, ncoarse 
            do i = 2, i1, ncoarse
              kk = k - klow + 1
              jj = (j - 2) / ncoarse + 2
              ii = (i - 2) / ncoarse + 2
              cli(ii,jj,k) = ql0(i,j,k) * (1 &
                             - max(0.0_field_r, min(1.0_field_r, &
                               (tmp0(i,j,k) - tdn) / (tup - tdn))))
            end do
          end do
        end do
      end if

      if (lclw) then
        call ofile%get_pointer('clw', clw)
        !$acc parallel loop collapse(3) default(present) async 
        do k = klow, khigh
          do j = 2, j1, ncoarse 
            do i = 2, i1, ncoarse
              kk = k - klow + 1
              jj = (j - 2) / ncoarse + 2
              ii = (i - 2) / ncoarse + 2
              clw(ii,jj,kk) = ql0(i,j,k) * max(0.0_field_r, min(1.0_field_r, &
                               (tmp0(i,j,k) - tdn) / (tup - tdn)))
            end do
          end do
        end do
      end if

      if (lpli) then
        call ofile%get_pointer('pli', pli)
        iqr = get_tracer_index('qr')
        !$acc parallell loop collapse(3) default(present) async
        do k = klow, khigh
          do j = 2, j1, ncoarse
            do i = 2, i1, ncoarse
              kk = k - klow + 1
              jj = (j - 2) / ncoarse + 2
              ii = (i - 2) / ncoarse + 2
              pli(ii,jj,kk) = sv0(i,j,k,iqr) * (1 &
                              - max(0.0_field_r, min(1.0_field_r, &
                               (tmp0(i,j,k) - tdn) / (tup - tdn))))
            end do
          end do
        end do
      end if

      if (lplw) then
        call ofile%get_pointer('plw', plw)
        iqr = get_tracer_index('qr')
        !$acc parallel loop collapse(3) default(present) async
        do k = klow, khigh
          do j = 2, j1, ncoarse
            do i = 2, i1, ncoarse
              kk = k - klow + 1
              jj = (j - 2) / ncoarse + 2
              ii = (i - 2) / ncoarse + 2
              pli(ii,jj,kk) = sv0(i,j,k,iqr) - max(0.0_field_r, min(1.0_field_r, &
                               (tmp0(i,j,k) - tdn) / (tup - tdn)))
            end do
          end do
        end do
      end if
 
      if (lhus) then
        call ofile%get_pointer('hus', hus)
        !$acc kernels default(present) async
        hus(:,:,:) = qt0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh) &
                     - ql0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh) 
        !$acc end kernels
      end if

      if (lhur) then 
        call ofile%get_pointer('hur', hur)
        do k = klow, khigh
          do j = 2, j1
            do i = 2, i1
              kk = k - klow + 1
              jj = (j - 2) / ncoarse + 2
              ii = (i - 2) / ncoarse + 2
              hur(ii,jj,kk) = 100 * (qt0(i,j,k) - ql0(i,j,k)) &
                              / calc_qsat(tmp0(i,j,k), presf(k)) 
            end do
          end do
        end do
      end if

      if (ltntr) then 
        call ofile%get_pointer('tntr', tntr)
        !$acc parallel loop collapse(3) default(present)
        do k = klow, khigh
          do j = 2, j1
            do i = 2, i1
              kk = k - klow + 1
              jj = (j - 2) / ncoarse + 2
              ii = (i - 2) / ncoarse + 2
              tntr(ii,jj,kk) = (- swd(i,j,k+1) - swu(i,j,k+1) &
                                + swd(i,j,k)   + swu(i,j,k) &
                                - lwd(i,j,k+1) - lwu(i,j,k+1) &
                                + lwd(i,j,k)   + lwu(i,j,k)) &
                               / (rhof(k) * exnf(k) * cp * dzf(k))
            end do
          end do
        end do
      end if

      if (ltntrs) then 
        call ofile%get_pointer('tntrs', tntrs)
        !$acc parallel loop collapse(3) default(present)
        do k = klow, khigh
          do j = 2, j1
            do i = 2, i1
              kk = k - klow + 1
              jj = (j - 2) / ncoarse + 2
              ii = (i - 2) / ncoarse + 2
              tntrs(ii,jj,kk) = (- swd(i,j,k+1) - swu(i,j,k+1) &
                                 + swd(i,j,k)   + swu(i,j,k)) &
                                / (rhof(k) * exnf(k) * cp * dzf(k))
            end do
          end do
        end do
      end if

      if (ltntrl) then 
        call ofile%get_pointer('tntrl', tntrl)
        !$acc parallel loop collapse(3) default(present)
        do k = klow, khigh
          do j = 2, j1
            do i = 2, i1
              kk = k - klow + 1
              jj = (j - 2) / ncoarse + 2
              ii = (i - 2) / ncoarse + 2
              tntrl(ii,jj,kk) = (- lwd(i,j,k+1) - lwu(i,j,k+1) &
                                 + lwd(i,j,k)   + lwu(i,j,k)) &
                                / (rhof(k) * exnf(k) * cp * dzf(k))
            end do
          end do
        end do
      end if

    end if

  end subroutine fielddump

end module modfielddump
