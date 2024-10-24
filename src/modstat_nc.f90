!> \file modstat_nc.f90
!!  Background routines to write NetCDF output

!>
!!  Background routines to write NetCDF output.
!>
!! All calls to the netcdf library should be directed through here.
!! Inspired on the UCLA-LES routine by Bjorn Stevens.
!!  \author Thijs Heus,MPI-M
!!  \par Revision list
!!  \todo documentation
!!   \todo restartfiles in NetCDF?
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
module modstat_nc
    use netcdf
    use modprecision, only: field_r
    use modmpi,       only: myid
    implicit none
    logical :: lnetcdf = .true.
    logical :: lsync   = .false.     ! Sync NetCDF file after each writestat_*_nc
    logical :: lclassic = .false.    ! Create netCDF in CLASSIC format (less RAM usage, compression not supported)
    integer :: deflate = 2           ! Deflate level for netCDF files (only for NETCDF4 format)

    integer, save :: timeID=0, ztID=0, zmID=0, xtID=0, xmID=0, ytID=0, ymID=0,ztsID=0, zqID=0
    real(kind=4) :: nc_fillvalue = -999.
!> The only interface necessary to write data to netcdf, regardless of the dimensions.
    interface writestat_nc
      module procedure writestat_time_nc
      module procedure writestat_1D_nc
      module procedure writestat_2D_nc
      module procedure writestat_3D_nc
      module procedure writestat_3D_short_nc
    end interface writestat_nc

    !> Read a field from a netCDF file by its name
    interface read_nc_field
        module procedure read_nc_field_1D
    end interface read_nc_field

    interface read_nc_attribute
      module procedure read_nc_attribute_char
      module procedure read_nc_attribute_r4
      module procedure read_nc_attribute_r8
      module procedure read_nc_attribute_logical
    end interface read_nc_attribute

    private :: read_nc_field_1D

    private :: read_nc_attribute_char
    private :: read_nc_attribute_r4
    private :: read_nc_attribute_r8
    private :: read_nc_attribute_logical

contains


  subroutine initstat_nc
    use modglobal, only : ifnamopt,fname_options,checknamelisterror
    use modmpi,    only : mpierr,mpi_logical,comm3d,myid,D_MPI_BCAST
    implicit none

    integer             :: ierr

    namelist/NAMNETCDFSTATS/ &
    lnetcdf, lsync, lclassic, deflate

    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMNETCDFSTATS,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMNETCDFSTATS')
      write(6, NAMNETCDFSTATS)
      close(ifnamopt)
    end if

    call D_MPI_BCAST(lnetcdf    ,1, 0,comm3d,mpierr)
    call D_MPI_BCAST(lsync      ,1, 0,comm3d,mpierr)
    call D_MPI_BCAST(lclassic   ,1, 0,comm3d,mpierr)
    call D_MPI_BCAST(deflate    ,1, 0,comm3d,mpierr)

  end subroutine initstat_nc
!
! ----------------------------------------------------------------------
!> Subroutine Open_NC: Opens a NetCDF File and identifies starting record
!
  subroutine open_nc (fname, ncid,nrec,n1, n2, n3, ns,nq)
    use modglobal, only : author,version,rtimee
    use modversion, only : git_version
    implicit none
    integer, intent (out) :: ncid,nrec
    integer, optional, intent (in) :: n1, n2, n3, ns, nq
    character (len=40), intent (in) :: fname

    character (len=12):: date='',time=''
    integer :: iret,varid,ncall,RecordDimID
    real, allocatable :: xtimes(:)
    logical :: exans

    inquire(file=trim(fname),exist=exans)

    ncall = 0
    if (.not.exans) then
      call date_and_time(date,time)
      if (lclassic) then
         call nchandle_error(nf90_create(fname,NF90_CLASSIC_MODEL,ncid))
      else
         call nchandle_error(nf90_create(fname,NF90_NETCDF4,ncid))
      end if
      call nchandle_error(nf90_put_att(ncid,NF90_GLOBAL,'title',fname))
      call nchandle_error(nf90_put_att(ncid,NF90_GLOBAL,'history','Created on '//trim(date)//' at '//trim(time)))
      call nchandle_error(nf90_put_att(ncid, NF90_GLOBAL, 'Source',trim(version)//' git: '//trim(git_version)))
      call nchandle_error(nf90_put_att(ncid, NF90_GLOBAL, 'Author',trim(author)))
      call nchandle_error(nf90_def_dim(ncID, 'time', NF90_UNLIMITED, timeID))
      if (present(n1)) then
         call nchandle_error(nf90_def_dim(ncID, 'xt', n1, xtID))
         call nchandle_error(nf90_def_dim(ncID, 'xm', n1, xmID))
         call nchandle_error(nf90_def_var(ncID,'xt',NF90_FLOAT,xtID ,VarID))
         call nchandle_error(nf90_put_att(ncID,VarID,'long_name','West-East displacement of cell centers'))
         call nchandle_error(nf90_put_att(ncID,VarID,'units','m'))
         call nchandle_error(nf90_def_var(ncID,'xm',NF90_FLOAT,xmID,VarID))
         call nchandle_error(nf90_put_att(ncID,VarID,'long_name','West-East displacement of cell edges'))
         call nchandle_error(nf90_put_att(ncID,VarID,'units','m'))
      end if
      if (present(n2)) then
         call nchandle_error(nf90_def_dim(ncID, 'yt', n2, ytID))
         call nchandle_error(nf90_def_dim(ncID, 'ym', n2, ymID))
         call nchandle_error(nf90_def_var(ncID,'yt',NF90_FLOAT,ytID ,VarID))
         call nchandle_error(nf90_put_att(ncID,VarID,'long_name','South-North displacement of cell centers'))
         call nchandle_error(nf90_put_att(ncID,VarID,'units','m'))
         call nchandle_error(nf90_def_var(ncID,'ym',NF90_FLOAT,ymID,VarID))
         call nchandle_error(nf90_put_att(ncID,VarID,'long_name','South-North displacement of cell edges'))
         call nchandle_error(nf90_put_att(ncID,VarID,'units','m'))
      end if
      if (present(n3)) then
         call nchandle_error(nf90_def_dim(ncID, 'zt', n3, ztID))
         call nchandle_error(nf90_def_dim(ncID, 'zm', n3, zmID))
         call nchandle_error(nf90_def_var(ncID,'zt',NF90_FLOAT,(/ztID/) ,VarID))
         call nchandle_error(nf90_put_att(ncID,VarID,'long_name','Vertical displacement of cell centers'))
         call nchandle_error(nf90_put_att(ncID,VarID,'units','m'))
         call nchandle_error(nf90_def_var(ncID,'zm',NF90_FLOAT,(/zmID/),VarID))
         call nchandle_error(nf90_put_att(ncID,VarID,'long_name','Vertical displacement of cell edges'))
         call nchandle_error(nf90_put_att(ncID,VarID,'units','m'))
      end if
      if (present(ns)) then
         call nchandle_error(nf90_def_dim(ncID, 'zts', ns, ztsID))
         call nchandle_error(nf90_def_var(ncID,'zts',NF90_FLOAT,(/ztsID/) ,VarID))
         call nchandle_error(nf90_put_att(ncID,VarID,'long_name','Soil level depth of cell centers'))
         call nchandle_error(nf90_put_att(ncID,VarID,'units','m'))
      end if
      if (present(nq)) then
         call nchandle_error(nf90_def_dim(ncID, 'zq', nq, zqID))
         call nchandle_error(nf90_def_var(ncID,'zq',NF90_FLOAT,(/zqID/) ,VarID))
         call nchandle_error(nf90_put_att(ncID,VarID,'long_name','Heights of interface levels'))
         call nchandle_error(nf90_put_att(ncID,VarID,'units','m'))
      end if

    else
       nrec = 0
       ncall= 0
       call nchandle_error(nf90_open (trim(fname), NF90_WRITE, ncid))
       call nchandle_error(nf90_inquire(ncid, unlimitedDimId = RecordDimID))
       call nchandle_error(nf90_inquire_dimension(ncid, RecordDimID, len=nrec))
       if (nrec>0) then
        call nchandle_error(nf90_inq_varid(ncid,'time',timeID))
        allocate (xtimes(nrec))
        call nchandle_error(nf90_get_var(ncid, timeId, xtimes(1:nrec)))

        ! Find the index where writing should continue.
        ! The next record to be written is ncall+1
        ! A warm start run does not write statistics immediately, thus
        ! fields up to and including the current time should be preserved.
        do while(xtimes(ncall+1)  <= rtimee)
            ncall=ncall+1
            if (ncall >= nrec) exit
        end do
        deallocate(xtimes)
       end if
       if (present(n1)) then
         iret = nf90_inq_dimid(ncid,'xt',xtId)
         iret = nf90_inq_dimid(ncid,'xm',xmId)
       end if
       if (present(n2)) then
         iret = nf90_inq_dimid(ncid,'yt',ytId)
         iret = nf90_inq_dimid(ncid,'ym',ymId)
       end if
       if (present(n3)) then
         iret = nf90_inq_dimid(ncid,'zt',ztId)
         iret = nf90_inq_dimid(ncid,'zm',zmId)
       end if
       if (present(ns)) then
         iret = nf90_inq_dimid(ncid,'zts',ztsId)
       end if
       if (present(nq)) then
         iret = nf90_inq_dimid(ncid,'zq',zqId)
       end if
    end if
    nrec = ncall

    iret = nf90_enddef(ncID)
    ! Fails with  NetCDF: Operation not allowed in data mode
    ! when the netCDF files already exist

  end subroutine open_nc

  !
  ! ----------------------------------------------------------------------
  !> Subroutine Define_NC: Defines the structure of the nc file (if not
  !! already open)
  !
  subroutine define_nc(ncID, nVar, sx)
    implicit none
    integer, intent (in) :: nVar, ncID
    character (*), intent (in) :: sx(nVar,4)

    integer, save ::  dim_mttt(4) = 0, dim_tmtt(4) = 0, dim_ttmt(4) = 0, dim_tttt(4) = 0, &
                      dim_tt(2)= 0, dim_mt(2)= 0,dim_t0tt(3)=0,dim_m0tt(3)=0,dim_t0mt(3)=0,dim_tt0t(3)=0, &
                      dim_mt0t(3)=0,dim_tm0t(3)=0,dim_0ttt(3)=0,dim_0mtt(3)=0,dim_0tmt(3)=0,&
                      dim_tts(2)=0,dim_t0tts(3)=0,dim_0ttts(3)=0,dim_tttts(4)=0,dim_qt(2)=0

    integer :: iret, n, VarID
    ! These calls are allowed to fail - not all files have all dimensions
    iret = nf90_inq_dimid(ncid,'time',timeId)
    iret = nf90_inq_dimid(ncid,'xt',xtId)
    iret = nf90_inq_dimid(ncid,'xm',xmId)
    iret = nf90_inq_dimid(ncid,'yt',ytId)
    iret = nf90_inq_dimid(ncid,'ym',ymId)
    iret = nf90_inq_dimid(ncid,'zt',ztId)
    iret = nf90_inq_dimid(ncid,'zm',zmId)
    iret = nf90_inq_dimid(ncid,'zts',ztsId)
    iret = nf90_inq_dimid(ncid,'zq',zqId)
    iret = nf90_redef(ncid)
    dim_tt = (/ztId,timeId/)
    dim_mt = (/zmId,timeId/)

    dim_t0tt= (/xtID,ztID,timeId/)! thermo point
    dim_t0mt= (/xtID,zmID,timeId/)! zpoint
    dim_m0tt= (/xmID,ztID,timeId/)! upoint
    dim_tt0t= (/xtID,ytID,timeId/)! thermo point
    dim_tm0t= (/xtID,ymID,timeId/)! vpoint
    dim_mt0t= (/xmID,ytID,timeId/)! upoint
    dim_0ttt= (/ytID,ztID,timeId/)! thermo point
    dim_0tmt= (/ytID,zmID,timeId/)! wpoint
    dim_0mtt= (/ymID,ztID,timeId/)! vpoint

    dim_tttt= (/xtID,ytID,ztID,timeId/)! thermo point
    dim_ttmt= (/xtID,ytID,zmID,timeId/)! zpoint
    dim_mttt= (/xmID,ytID,ztID,timeId/)! upoint
    dim_tmtt= (/xtID,ymID,ztId,timeId/)! ypoint

    dim_tts = (/ztsId,timeId/)
    dim_t0tts= (/xtID,ztsID,timeId/)! thermo soil point
    dim_0ttts= (/ytID,ztsID,timeId/)! thermo point
    dim_tttts= (/xtID,ytID,ztsID,timeId/)! thermo point

    dim_qt = (/zqId,timeId/)

    do n=1,nVar
      iret = nf90_inq_varid(ncid, trim(sx(n,1)), VarID)
      if (iret == 0) cycle
      select case(trim(sx(n,4)))
        case ('time')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,(/timeID/) ,VarID)
        case ('tt')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_tt ,VarID)
        case ('mt')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_mt,VarID)
  !2D Fields
        case ('t0tt')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_t0tt,VarID)
        case ('t0mt')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_t0mt,VarID)
        case ('m0tt')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_m0tt,VarID)
        case ('tt0t')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_tt0t,VarID)
        case ('tm0t')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_tm0t,VarID)
        case ('mt0t')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_mt0t,VarID)
        case ('0ttt')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_0ttt,VarID)
        case ('0tmt')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_0tmt,VarID)
        case ('0mtt')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_0mtt,VarID)
  !3D Fields
        case ('tttt')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_tttt,VarID)
        case ('mttt')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_mttt,VarID)
        case ('tmtt')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_tmtt,VarID)
        case ('ttmt')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_ttmt,VarID)
!Soil fields
        case ('tts')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_tts ,VarID)
        case ('t0tts')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_t0tts,VarID)
        case ('0ttts')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_0ttts,VarID)
        case ('tttts')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_tttts,VarID)

!Quadrant analysis fields
        case('qt')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_qt ,VarID)
        case default
        print *, 'ABORTING: Bad dimensional information ',sx(n,:)
        stop
        ! call appl_abort(0)

      end select
      if (iret/=0) then
        write (*,*) 'nvar', nvar, sx(n,:)
        call nchandle_error(iret)
      end if

      if (deflate > 0 .and. .not. lclassic) then
         call nchandle_error(nf90_def_var_deflate(ncid,varID, 0, 1, deflate_level = deflate))
         ! NETCDF4 only
      end if
      call nchandle_error(nf90_put_att(ncID,VarID,'long_name',sx(n,2)))
      call nchandle_error(nf90_put_att(ncID,VarID,'units',sx(n,3)))
      call nchandle_error(nf90_put_att(ncid, VarID, '_FillValue',nc_fillvalue))
      !fails "NetCDF: Not a valid data type or _FillValue type mismatch"
      !on Fugaku with netCDF-Fortran 4.5.2, netCDF 4.7.3

    end do
    iret= nf90_enddef(ncID)
  end subroutine define_nc

  subroutine redefine_nc(ncid)
  implicit none
    integer, intent(in) :: ncid
    call nchandle_error(nf90_redef(ncid))
  end subroutine redefine_nc

  subroutine exitstat_nc(ncid)

   implicit none
   integer, intent(in) :: ncid
   integer status

   status = nf90_close(ncid)
   call nchandle_error(status)
 end subroutine exitstat_nc

 subroutine writestat_dims_nc(ncid, ncoarse, klow, proc)
    ! optional arguments ncoarse (coarsegraining in the horizontal directions)
    !                    klow    (lower bound for z. Upper bound is taken from the size of the dimension)
    !                    proc    (if present and length on horizontal cooridinates is 1 include processor starting edges and center)
    use modglobal, only : dx,dy,zf,zh,jmax,imax,kmax
    use modsurfdata, only : zsoilc,isurf
    use modlsm, only : z_soil
    use modmpi, only : myidx,myidy
    implicit none
    integer, intent(in) :: ncid
    integer, optional, intent(in) :: ncoarse, klow
    logical, optional, intent(in) :: proc
    integer             :: i=0,iret,length,varid, nc
    integer             :: kl
    logical             :: lproc
    lproc = .false.
    if (present(ncoarse)) then
      nc = ncoarse
    else
      nc = 1
    end if
    if (present(klow)) then
      kl = klow
    else
      kl = 1
    end if
    if (present(proc)) then
      lproc = .true.
    end if

    if (.not. lproc) then
      iret = nf90_inq_varid(ncid, 'xt', VarID)
      if (iret==0) iret=nf90_inquire_dimension(ncid, xtID, len=length)
      if (iret==0) iret = nf90_put_var(ncid, varID, (/(dx*(0.5+nc*i)+myidx*imax*dx,i=0,length-1)/),(/1/))
      iret = nf90_inq_varid(ncid, 'xm', VarID)
      if (iret==0) iret=nf90_inquire_dimension(ncid, xmID, len=length)
      if (iret==0) iret = nf90_put_var(ncid, varID, (/(dx*nc*i+myidx*imax*dx,i=0,length-1)/),(/1/))

      iret = nf90_inq_varid(ncid, 'yt', VarID)
      if (iret==0) iret=nf90_inquire_dimension(ncid, ytID, len=length)
      if (iret==0) iret = nf90_put_var(ncid, varID, (/(dy*(0.5+nc*i)+myidy*jmax*dy,i=0,length-1)/),(/1/))
      iret = nf90_inq_varid(ncid, 'ym', VarID)
      if (iret==0) iret=nf90_inquire_dimension(ncid, ymID, len=length)
      if (iret==0) iret = nf90_put_var(ncid, varID, (/(dy*nc*i+myidy*jmax*dy,i=0,length-1)/),(/1/))
    else
      iret = nf90_inq_varid(ncid, 'xt', VarID)
      if (iret==0) iret = nf90_put_var(ncid, varID, (/(0.5*dx*imax+myidx*imax*dx)/),(/1/))
      iret = nf90_inq_varid(ncid, 'xm', VarID)
      if (iret==0) iret = nf90_put_var(ncid, varID, (/(myidx*imax*dx)/),(/1/))

      iret = nf90_inq_varid(ncid, 'yt', VarID)
      if (iret==0) iret = nf90_put_var(ncid, varID, (/(0.5*dy*jmax+myidy*jmax*dy)/),(/1/))
      iret = nf90_inq_varid(ncid, 'ym', VarID)
      if (iret==0) iret = nf90_put_var(ncid, varID, (/(myidy*jmax*dy)/),(/1/))
    end if

    iret = nf90_inq_varid(ncid, 'zt', VarID)
    if (iret==0) iret=nf90_inquire_dimension(ncid,ztID, len=length)
    if (iret==0) iret = nf90_put_var(ncid, varID, zf(kl:kl+length-1),(/1/))
    iret = nf90_inq_varid(ncid, 'zm', VarID)
    if (iret==0) iret=nf90_inquire_dimension(ncid, zmID, len=length)
    if (iret==0) iret = nf90_put_var(ncid, varID, zh(kl:kl+length-1),(/1/))
    if (isurf==1) then
      iret = nf90_inq_varid(ncid, 'zts', VarID)
      if (iret==0) iret = nf90_inquire_dimension(ncid, ztsID, len=length)
      if (iret==0) iret = nf90_put_var(ncid, varID, zsoilc(1:length),(/1/))
    else if (isurf==11) then
      iret = nf90_inq_varid(ncid, 'zts', VarID)
      if (iret==0) iret = nf90_inquire_dimension(ncid, ztsID, len=length)
      if (iret==0) iret = nf90_put_var(ncid, varID, z_soil(1:length),(/1/))
    end if

  end subroutine writestat_dims_nc

  subroutine writestat_dims_q_nc(ncid,k1,k2)
    use modglobal, only : dx,dy,zf,zh,jmax,imax
    use modsurfdata, only : zsoilc,isurf
    use modlsm, only : z_soil
    use modmpi, only : myidx,myidy
    implicit none
    integer, intent(in) :: ncid,k1,k2
    integer             :: i=0,iret,length,varid
    iret = nf90_inq_varid(ncid, 'xt', VarID)
    if (iret==0) iret=nf90_inquire_dimension(ncid, xtID, len=length)
    if (iret==0) iret = nf90_put_var(ncid, varID, (/(dx*(0.5+i)+myidx*imax*dx,i=0,length-1)/),(/1/))
    iret = nf90_inq_varid(ncid, 'xm', VarID)
    if (iret==0) iret=nf90_inquire_dimension(ncid, xmID, len=length)
    if (iret==0) iret = nf90_put_var(ncid, varID, (/(dx*i+myidx*imax*dx,i=0,length-1)/),(/1/))

    iret = nf90_inq_varid(ncid, 'yt', VarID)
    if (iret==0) iret=nf90_inquire_dimension(ncid, ytID, len=length)
    if (iret==0) iret = nf90_put_var(ncid, varID, (/(dy*(0.5+i)+myidy*jmax*dy,i=0,length-1)/),(/1/))
    iret = nf90_inq_varid(ncid, 'ym', VarID)
    if (iret==0) iret=nf90_inquire_dimension(ncid, ymID, len=length)
    if (iret==0) iret = nf90_put_var(ncid, varID, (/(dy*i+myidy*jmax*dy,i=0,length-1)/),(/1/))

    iret = nf90_inq_varid(ncid, 'zt', VarID)
    if (iret==0) iret=nf90_inquire_dimension(ncid,ztID, len=length)
    if (iret==0) iret = nf90_put_var(ncid, varID, zf(1:length),(/1/))
    iret = nf90_inq_varid(ncid, 'zm', VarID)
    if (iret==0) iret=nf90_inquire_dimension(ncid, zmID, len=length)
    if (iret==0) iret = nf90_put_var(ncid, varID, zh(1:length),(/1/))
    if (isurf==1) then
      iret = nf90_inq_varid(ncid, 'zts', VarID)
      if (iret==0) iret = nf90_inquire_dimension(ncid, ztsID, len=length)
      if (iret==0) iret = nf90_put_var(ncid, varID, zsoilc(1:length),(/1/))
  else if (isurf==11) then
      iret = nf90_inq_varid(ncid, 'zts', VarID)
      if (iret==0) iret = nf90_inquire_dimension(ncid, ztsID, len=length)
      if (iret==0) iret = nf90_put_var(ncid, varID, z_soil(1:length),(/1/))
    end if

    iret = nf90_inq_varid(ncid, 'zq', VarID)
    if (iret==0) iret=nf90_inquire_dimension(ncid,zqID, len=length)
    if (length .ne. (1+k2-k1)) then
      print *,"k1     = ",k1
      print *,"k2     = ",k2
      print *,"length = ",length
      stop "Problem in writestat_dims_q_nc: not matching lengths"
    endif
    if (iret==0) iret = nf90_put_var(ncid, varID, zh(k1:k2),(/1/))

  end subroutine writestat_dims_q_nc

  subroutine writestat_time_nc(ncid,nvar,ncname,vars,nrec,lraise)
    implicit none
    integer, intent(in)                      :: ncid,nvar
    integer, intent(inout)                   :: nrec
    real,dimension(nvar),intent(in)          :: vars
    character(*), dimension(:,:),intent(in)  :: ncname
    logical, intent(in)                      :: lraise

    integer :: n,varid
    if(lraise) then
      nrec = nrec+1
    end if
    do n=1,nvar
       call nchandle_error(nf90_inq_varid(ncid, ncname(n,1), VarID))
       call nchandle_error(nf90_put_var(ncid, VarID, vars(n), start=(/nrec/)))
    end do
    if (lsync) call sync_nc(ncid)
  end subroutine writestat_time_nc

  subroutine writestat_1D_nc(ncid,nvar,ncname,vars,nrec,dim1)
    implicit none
    integer, intent(in)                      :: ncid,nvar,dim1
    integer, intent(in)                      :: nrec
    real,dimension(dim1,nvar),intent(in)     :: vars
    character(*), dimension(:,:),intent(in)  :: ncname

    integer :: n,varid
    do n=1,nvar
       call nchandle_error(nf90_inq_varid(ncid, ncname(n,1), VarID))
       call nchandle_error(nf90_put_var(ncid, VarID, vars(1:dim1,n),(/1,nrec/),(/dim1,1/)))
    end do
    if (lsync) call sync_nc(ncid)
  end subroutine writestat_1D_nc

  subroutine writestat_2D_nc(ncid,nvar,ncname,vars,nrec,dim1,dim2)
    implicit none
    integer, intent(in)                      :: ncid,nvar,dim1,dim2
    integer, intent(in)                      :: nrec
    real,dimension(:,:,:),intent(in)         :: vars
    character(*), dimension(:,:),intent(in)  :: ncname

    integer :: n,varid
    do n=1,nvar
       call nchandle_error(nf90_inq_varid(ncid, ncname(n,1), VarID))
       call nchandle_error(nf90_put_var(ncid, VarID, vars(1:dim1,1:dim2,n),(/1,1,nrec/),(/dim1,dim2,1/)))
    end do
    if (lsync) call sync_nc(ncid)
  end subroutine writestat_2D_nc
  subroutine writestat_3D_nc(ncid,nvar,ncname,vars,nrec,dim1,dim2,dim3)
    implicit none
    integer, intent(in)                      :: ncid,nvar,dim1,dim2,dim3
    integer, intent(in)                      :: nrec
    real,dimension(dim1,dim2,dim3,nvar),intent(in)       :: vars
    character(*), dimension(:,:),intent(in)  :: ncname

    integer :: n,varid
    do n=1,nvar
       call nchandle_error(nf90_inq_varid(ncid, ncname(n,1), VarID))
       call nchandle_error(nf90_put_var(ncid, VarID, vars(1:dim1,1:dim2,1:dim3,n),(/1,1,1,nrec/),(/dim1,dim2,dim3,1/)))
    end do
    if (lsync) call sync_nc(ncid)
  end subroutine writestat_3D_nc
  subroutine writestat_3D_short_nc(ncid,nvar,ncname,vars,nrec,dim1,dim2,dim3)
    implicit none
    integer, intent(in)                      :: ncid,nvar,dim1,dim2,dim3
    integer, intent(in)                      :: nrec
    integer(KIND=selected_int_kind(4)),dimension(dim1,dim2,dim3,nvar),intent(in)       :: vars
    character(*), dimension(:,:),intent(in)  :: ncname

    integer :: n,varid
    do n=1,nvar
       call nchandle_error(nf90_inq_varid(ncid, ncname(n,1), VarID))
       call nchandle_error(nf90_put_var(ncid, VarID, vars(1:dim1,1:dim2,1:dim3,n),(/1,1,1,nrec/),(/dim1,dim2,dim3,1/)))
    end do
    if (lsync) call sync_nc(ncid)
  end subroutine writestat_3D_short_nc

  !> Read a 1D field from a netCDF by its name.
  !!
  !! Optionally, a fill value can be provided, which the is used in case the
  !! requested variable is not found in the netCDF file.
  !! @param ncid ID of the opened netCDF file.
  !! @param varname Name of the variable to read.
  !! @param array Array to fill.
  !! @param fillvalue Default value to fill array with in case the variable
  !! is not found.
  subroutine read_nc_field_1D(ncid, varname, array, start, count, fillvalue)
    integer,       intent(in)           :: ncid
    character(*),  intent(in)           :: varname
    real(field_r), intent(out)          :: array(:)
    integer,       intent(in), optional :: start
    integer,       intent(in), optional :: count
    real(field_r), intent(in), optional :: fillvalue

    integer :: varid
    integer :: ierr
    integer :: start_(1), count_(1)

    ierr = nf90_inq_varid(ncid, varname, varid)

    if (present(start)) then 
      start_(1) = start
    else
      start_(1) = 1
    end if

    if (present(count)) then
      count_(1) = count
    else
      count_(1) = size(array)
    end if

    select case (ierr)
      case (NF90_ENOTVAR)
        if (present(fillvalue)) then
          array(:) = fillvalue
        else
          call nchandle_error(ierr)
        end if
      case (NF90_NOERR)
        ierr = nf90_get_var(ncid, varid, array, start=start_, count=count_)
      case default
        call nchandle_error(ierr)
    end select

  end subroutine read_nc_field_1D

  subroutine read_nc_attribute_char(ncid, varid, attname, value, default)
    integer,      intent(in)           :: ncid
    integer,      intent(in)           :: varid
    character(*), intent(in)           :: attname
    character(*), intent(out)          :: value
    character(*), intent(in), optional :: default

    integer :: ierr

    ierr = nf90_get_att(ncid, varid, attname, value)

    if (ierr == NF90_ENOTATT) then
      if (present(default)) value = default
    else
      call nchandle_error(ierr)
    end if

  end subroutine read_nc_attribute_char

  subroutine read_nc_attribute_r4(ncid, varid, attname, value, default)
    integer,      intent(in)           :: ncid
    integer,      intent(in)           :: varid
    character(*), intent(in)           :: attname
    real(4),      intent(out)          :: value
    real(4),      intent(in), optional :: default

    integer :: ierr

    ierr = nf90_get_att(ncid, varid, attname, value)

    if (ierr == NF90_ENOTATT) then
      if (present(default)) value = default
    else
      call nchandle_error(ierr)
    end if

  end subroutine read_nc_attribute_r4

  subroutine read_nc_attribute_r8(ncid, varid, attname, value, default)
    integer,      intent(in)           :: ncid
    integer,      intent(in)           :: varid
    character(*), intent(in)           :: attname
    real(8),      intent(out)          :: value
    real(8),      intent(in), optional :: default

    integer :: ierr

    ierr = nf90_get_att(ncid, varid, attname, value)

    if (ierr == NF90_ENOTATT) then
      if (present(default)) value = default
    else
      call nchandle_error(ierr)
    end if

  end subroutine read_nc_attribute_r8

  subroutine read_nc_attribute_logical(ncid, varid, attname, value, default)
    integer,      intent(in)           :: ncid
    integer,      intent(in)           :: varid
    character(*), intent(in)           :: attname
    logical,      intent(out)          :: value
    logical,      intent(in), optional :: default

    integer :: ierr, value_
    
    ierr = nf90_get_att(ncid, varid, attname, value_)

    select case (ierr)
      case (NF90_ENOTATT)
        if (present(default)) then
          value = default
        else
          call nchandle_error(ierr)
        end if
      case (NF90_NOERR)
        select case (value_)
          case (1)
            value = .true.
          case (0)
            value = .false.
          case default
            if (myid == 0) write(*,*) "Invalid value provided for ", attname
            stop
        end select
      case default
        call nchandle_error(ierr)
    end select
    
  end subroutine read_nc_attribute_logical

  subroutine ncinfo(out,in1,in2,in3,in4)

    implicit none
    character(*), dimension(4),intent(out) ::out
    character(*), intent(in) ::in1,in2,in3,in4
    out(1) = in1
    out(2) = in2
    out(3) = in3
    out(4) = in4
  end subroutine ncinfo

  subroutine nchandle_error(status)
    use netcdf
    implicit none

    integer, intent(in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      call abort
    end if

  end subroutine nchandle_error

  subroutine sync_nc(ncid)
    implicit none
    integer, intent(in) :: ncid
    integer :: iret

    iret = nf90_sync(ncid)
    call nchandle_error(iret)
  end subroutine sync_nc

  subroutine nctiminfo(info)
    use modglobal, only: xyear, xday, xtime
    implicit none
    character(*), dimension(4), intent(out) :: info
    character(len=33)                       :: unitstr
    integer                                 :: dt(6)

    if (xyear > 0 .and. xday > 0) then
      call get_date_time(xyear,xday,xtime,dt)
      write(unitstr,'(A14,I4.4,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2)') 'seconds since ',dt(1),'-',dt(2),'-',dt(3),'T',dt(4),':',dt(5),':',dt(6)
      call ncinfo(info(:),'time','Time',unitstr,'time')
    else
      call ncinfo(info(:),'time','Time','s','time')
    end if
  end subroutine nctiminfo

  ! Utility converting year + day to full date
  ! day = 1 is first day of the year
  ! hours > 24 advances the day, and day > 365 or 366 advances the year
  subroutine get_date_time(year, day, hours, datetime)
    implicit none
    integer, intent(in)     :: year
    real, intent(in)        :: day
    real, intent(in)        :: hours
    integer, intent(out)    :: datetime(6)
    real                    :: rh
    integer                 :: dd,hh,mm,ss,date(3)

    rh = (day - floor(day)) * 24 + hours
    hh = floor(rh)
    mm = floor((rh - hh) * 60)
    ss = floor((rh - hh - mm/60.) * 3600)
    dd = hh / 24
    hh = hh - 24 * dd
    call get_date(year,floor(day) + dd, date)
    datetime(1:3) = date(1:3)
    datetime(4) = hh
    datetime(5) = mm
    datetime(6) = ss
  end subroutine get_date_time

  ! Utility converting year + day to full date
  ! day = 1 is first day of the year
  subroutine get_date(year, day, date)
    implicit none
    integer, intent(in)     :: year
    integer, intent(in)     :: day
    integer, intent(out)    :: date(3)
    integer                 :: i,yy,mm,dd,dsum,ndays
    integer                 :: mdays(12)

    yy = year
    dsum = 0
    do
      ndays = 365
      if (leap_year(yy)) then
        ndays = 366
      endif
      if (dsum + ndays >= day) then
        exit
      endif
      dsum = dsum + ndays
      yy = yy + 1
    enddo
    dd = day - dsum
    mdays = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
    if (leap_year(yy)) then
      mdays(2) = 29
    endif
    dsum = 0
    mm = 0
    do i=1,12
      mm = i
      if (dsum + mdays(i) >= dd) then
        exit
      endif
      dd = dd - mdays(i)
    enddo
    date(1) = yy
    date(2) = mm
    date(3) = dd
  end subroutine get_date

  ! Utility checking whether the input year is a leap year
  function leap_year(year) result(isleap)
    implicit none
    integer, intent(in) :: year
    logical             :: isleap

    isleap = ((mod(year, 4) == 0) .and. (mod(year, 100) /= 0)) .or. (mod(year, 400) == 0)
  end function leap_year

end module modstat_nc
