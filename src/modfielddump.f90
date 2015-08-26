!> \file modfielddump.f90
!!  Dumps 3D fields of several variables

!>
!!  Dumps 3D fields of several variables
!>
!!  Dumps 3D fields of several variables Written to wb*.myid.expnr
!! If netcdf is true, this module leads the fielddump.myid.expnr.nc output
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
  use modncrestart, only : nMaxFds,iWhichVars
  use modglobal,    only : longint

implicit none
private
PUBLIC :: initfielddump,fielddump,exitfielddump
save
!NetCDF variables
  integer           :: ncid                     !< Number of variables (if bulkmicro->nvar+2)
  character(80)     :: fname = 'fielddump_hhhmmss.xxx.nc'

  logical          :: lfielddump= .false.       !< switch to enable the fielddump (on/off)
  real             :: beghr(20)=-1., &          !< array of begin and end hours inbetween which
                      endhr(20)=-1.             !< fields are dumped at 'interval' times
  real             :: interval=300.             !< time interval between dumps in sec
  integer          :: kfieldtop=-1              !< top level included in the fielddumps
  integer          :: iInterval                 !< integer interval time
  integer          :: tNextDump                 !< time for the next field dump
  integer          :: nhours                    !< number of beg/end hours
  integer          :: iNextHr=1                 !< index of begin hour >= current time

contains
!> Initializing fielddump. Read out the namelist, initializing the variables
  subroutine initfielddump
    use modmpi,       only : myid,my_real,comm3d,mpi_logical,mpi_integer,cmyid,nprocs
    use modglobal,    only : cexpnr,ifnamopt,fname_options,dtmax,btime,runtime,timee,tres,dt_lim, &
                             imax,jmax,kmax,ladaptive,zf,zh
    use modNcRestart, only : ncDims,ncVars,nDims,ncFd,nFds,&
                             xnr,xhnr,ynr,yhnr,znr,zhnr,zfdnr,zfdhnr,nMaxVars
    use modboundary,  only : ksp
    use modmicrodata, only : imicro,imicro_bulk
    use modstat_nc,   only : lnetcdf,open_nc,define_nc,redefine_nc,ncinfo,writestat_dims_nc
    implicit none
    integer :: ierr,i,iD,iCur

    namelist/NAMFIELDDUMP/ &
    lfielddump,beghr,endhr,interval,kfieldtop,iWhichVars

    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMFIELDDUMP,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions NAMFIELDDUMP'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMFIELDDUMP'
      endif
      write(6 ,NAMFIELDDUMP)
      close(ifnamopt)
    end if
    call MPI_BCAST(lfielddump  , 1,MPI_LOGICAL,0,comm3d,ierr)
    call MPI_BCAST(beghr       ,20,MY_REAL    ,0,comm3d,ierr)
    call MPI_BCAST(endhr       ,20,MY_REAL    ,0,comm3d,ierr)
    call MPI_BCAST(interval    , 1,MY_REAL    ,0,comm3d,ierr)
    call MPI_BCAST(kfieldtop   , 1,MPI_INTEGER,0,comm3d,ierr)
    call MPI_BCAST(iWhichVars  ,nMaxFds,MPI_INTEGER,0,comm3d,ierr)

    if(.not.(lfielddump)) return

    nhours = count(.not. beghr < 0.)
    if (nhours /= count(.not. endhr < 0.) ) then
      if (myid == 0) write(*,*) 'begin hours are ',beghr(1:nhours)
      if (myid == 0) write(*,*) 'end hours are '  ,endhr(1:nhours)
      stop 'ERROR: # of begin hours /= # of end hours; modfielddump'
    end if

    if (any(endhr(1:nhours) - beghr(1:nhours) < 0.)) then
      stop 'ERROR: some end hour(s) smaller than begin hours.'
    end if

    if (kfieldtop<0 .or. kfieldtop>kmax+1) then
      if (myid == 0) then 
        write(*,*) 'WARNING: no or invalid kfieldtop given.'
        write(*,*) 'WARNING: kfieldtop set to just below sponge layer: ',ksp-1
      end if
      kfieldtop = ksp-1
    end if

    ! Ignore all values that are not allowed
    if (imicro==imicro_bulk) then
      if (myid==0) write(*,*) 'WARNING: ignored ',count(iWhichVars>9),' field(s)'
      where(iWhichVars>9) iWhichVars=-1
    else
      if (myid==0) write(*,*) 'WARNING: ignored ',count(iWhichVars>7),' field(s)'
      where(iWhichVars>7) iWhichVars=-1
    end if

    if (all(iWhichVars<0)) then
      if (myid==0) then
        write(*,*) 'WARNING: no variables selected for modfielddump.'
        write(*,*) 'WARNING: doing all fields.'
      end if
      nFds = 9
      if (imicro==imicro_bulk) then
        nFds = 11
      end if
      iWhichVars(1:nFds) = (/ (i,i=1,nFds) /)
    else
      ! Make sure that zfull, zhalf, pfull and phalf are always written
      nFds = count(iWhichVars>0)
      iWhichVars(3:nFds+2) = iWhichVars(1:nFds)+2
      iWhichVars(1:2)      = (/1,2/)
    end if
    write(*,*) 'iWhichVars2 = ',iWhichVars
    ! Count the total number of fields to be written
    nFds = count(iWhichVars>0)

    ! Increase the number of dimensions by 2, to account for the fielddumps that do not go all the way to the top
    nDims = nDims+2
    ncDims(zfdnr)%name = 'zfd';   ncDims(zfdnr)%size = kfieldtop;
    ncDims(zfdnr)%values(1:ncDims(zfdnr)%size)=zf(1:kfieldtop)
    ncDims(zfdnr)%longname = 'Vertical displacement (full levels, fdumps)'; ncDims(zfdnr)%units = 'm'
    ncDims(zfdhnr)%name = 'zfdh';   ncDims(zfdhnr)%size = kfieldtop;
    ncDims(zfdhnr)%values(1:ncDims(zfdhnr)%size)=zh(1:kfieldtop)
    ncDims(zfdhnr)%longname = 'Vertical displacement (half levels, fdumps)'; ncDims(zfdhnr)%units = 'm'

    ! Write the hours in hundreds of seconds
    beghr = beghr * 3600/tres
    endhr = endhr * 3600/tres

    tNextDump = beghr(iNextHr)

    ! Do loop, in case of a restart run. Finds the appropriate next fielddump time
    do while (tNextDump<timee)
      iNextHr   = iNextHr + 1
      tNextDump = beghr(iNextHr)
      if (iNextHr >= nhours) tNextDump = tres*(btime + runtime) + 1.
    end do

    iInterval = interval / tres
    dt_lim = min(dt_lim,tNextDump)

    if (.not. ladaptive .and. abs(interval/dtmax-nint(interval/dtmax))>1e-4) then
      stop 'interval time should be an integer multiple of dtmax; modfielddump'
    end if

    if (lnetcdf) then

      fname(19:21) = cexpnr

      i = 1
      ncFd(i)%name='presf'; ncFd(i)%longname='Pressure at full levels'
        ncFd(i)%units='Pa'; ncFd(i)%dNr(1)=znr;   i=i+1
      ncFd(i)%name='presh'; ncFd(i)%longname='Pressure at half levels'
        ncFd(i)%units='Pa'; ncFd(i)%dNr(1)=zhnr;   i=i+1
      ncFd(i)%name='u'; ncFd(i)%longname='West-East velocity'                       !nr: 1
        !ncFd(i)%units='m/s';     ncFd(i)%dNr(1:3)=(/xhnr,yhnr,zfdnr/); i=i+1
        ncFd(i)%units='m/s';     ncFd(i)%dNr(1:3)=(/xnr,ynr,zfdnr/); i=i+1
      ncFd(i)%name='v'; ncFd(i)%longname='South-North velocity'                     !nr: 2
        !ncFd(i)%units='m/s';     ncFd(i)%dNr(1:3)=(/xhnr,yhnr,zfdnr/); i=i+1
        ncFd(i)%units='m/s';     ncFd(i)%dNr(1:3)=(/xnr,ynr,zfdnr/); i=i+1
      ncFd(i)%name='w'; ncFd(i)%longname='Vertical velocity'                        !nr: 3
        !ncFd(i)%units='m/s';     ncFd(i)%dNr(1:3)=(/xnr,ynr,zfdhnr/);  i=i+1
        ncFd(i)%units='m/s';     ncFd(i)%dNr(1:3)=(/xnr,ynr,zfdnr/);  i=i+1
      ncFd(i)%name='thl'; ncFd(i)%longname='Liquid water potential temperature'     !nr: 4
        ncFd(i)%units='K';       ncFd(i)%dNr(1:3)=(/xnr,ynr,zfdnr/);   i=i+1
      ncFd(i)%name='qt'; ncFd(i)%longname='Total specific humidity'                 !nr: 5
        ncFd(i)%units='kg/kg';   ncFd(i)%dNr(1:3)=(/xnr,ynr,zfdnr/);   i=i+1
      ncFd(i)%name='ql'; ncFd(i)%longname='Liquid water specific humidity'          !nr: 6
        ncFd(i)%units='kg/kg';   ncFd(i)%dNr(1:3)=(/xnr,ynr,zfdnr/);   i=i+1
      ncFd(i)%name='e'; ncFd(i)%longname='Subgrid scale turbulent kinetic energy'   !nr: 7
        ncFd(i)%units='m^2/s^2'; ncFd(i)%dNr(1:3)=(/xnr,ynr,zfdnr/);   i=i+1
      if (imicro==imicro_bulk) then
        ncFd(i)%name='qr'; ncFd(i)%longname='Rain water specific humidity'          !nr: 8
          ncFd(i)%units='kg/kg';   ncFd(i)%dNr(1:3)=(/xnr,ynr,zfdnr/);   i=i+1
        ncFd(i)%name='Nr'; ncFd(i)%longname='Rain droplet number density'           !nr: 9
          ncFd(i)%units='cm^-3';   ncFd(i)%dNr(1:3)=(/xnr,ynr,zfdnr/);   i=i+1
      end if

      if (myid==0) then
        write(*,*) 'MODFIELDDUMP: will write fields: ',ncFd(iWhichVars(1:nFds))%name
      end if

      ! Set all starting points and counts (number of values to write in each direction)
      ! required for the parallel writing
      do i=1,nFds
        iCur=iWhichVars(i)
        do iD=1,4
          if (ncFd(iCur)%dNr(iD)>0) then
            ncFd(iCur)%count(iD) = ncDims(ncFd(iCur)%dNr(iD))%size
            if (iD==2) ncFd(iCur)%count(iD) = ceiling(ncDims(ncFd(iCur)%dNr(iD))%size /real(nprocs))
          end if
        end do
        ncFd(iCur)%start = 1
        ncFd(iCur)%start(2) = jmax*myid+1
      end do

    else ! if lnetcdf==.false.

      if (myid == 0) then
        write(*,*) 'WARNING: lfielddump = .true. while lnetcdf = .false.'
        write(*,*) 'WARNING: No fielddumps will be made. Setting lfielddump = .false.'
        lfielddump = .false.
      end if

    end if ! if lnetcdf
    
  end subroutine initfielddump

!> Do fielddump. 
  subroutine fielddump
    use modglobal,  only : timee,rtimee,dt_lim,tres,timeleft,i1,j1,k1
    use modmpi,     only : myid
    use modstat_nc, only : lnetcdf
    use modncrestart, only : genNcFile,writeNcVar,iWhichVars,nFds,ncFd,closeNcFile
    use modfields,  only : presf,presh,u0,v0,w0,thl0,qt0,ql0,e120,sv0
    use modmicrodata, only : iqr,iNr
    implicit none

    integer :: isec,imin,ihour
    integer :: iV,iCur
    integer            :: t1,t2,count_rate,count_max

    if(.not.(lfielddump)) return
    
    if (timee<tNextDump) then
      dt_lim = min(dt_lim,tNextDump-timee)
    else ! if (timee>=tNextDump)
      ! Start the timing
      call system_clock(t1, count_rate, count_max)
      ! It's time for a fielddump.
      ! Make a filename first,...
      ihour = floor(rtimee/3600.)
      imin  = floor(mod(rtimee,3600.)/60.)
      isec  = floor(mod(rtimee,60.))
      write (fname(11:13),'(i3.3)') ihour
      write (fname(14:15),'(i2.2)') imin
      write (fname(16:17),'(i2.2)') isec
      if (myid==0) write(*,*) 'fname = ',fname
      
!      ! ...then call the subroutine that makes the actual file and puts data in.
      call genNcFile(fname)
      do iV=1,nFds
        iCur=iWhichVars(iV)
        select case(iCur)
          case(1); call writeNcVar(ncFd(iCur),presf(1:k1))
          case(2); call writeNcVar(ncFd(iCur),presh(1:k1))
          !case(3); call writeNcVar(ncFd(iCur),u0(1:i1,1:j1,1:k1))
          case(3); call writeNcVar(ncFd(iCur), &
            (u0(1:i1-1,1:j1-1,1:k1)+u0(2:i1,1:j1-1,1:k1)+u0(2:i1,1:j1-1,1:k1)+u0(2:i1,2:j1,1:k1))/4.)
          !case(4); call writeNcVar(ncFd(iCur),v0(1:i1,1:j1,1:k1))
          case(4); call writeNcVar(ncFd(iCur), &
            (v0(1:i1-1,1:j1-1,1:k1)+v0(2:i1,1:j1-1,1:k1)+v0(2:i1,1:j1-1,1:k1)+v0(2:i1,2:j1,1:k1))/4.)
          !case(5); call writeNcVar(ncFd(iCur),w0(1:i1,1:j1,1:k1))
          case(5); call writeNcVar(ncFd(iCur),(w0(2:i1,2:j1,1:k1-1)+w0(2:i1,2:j1,2:k1))/2.)
          case(6); call writeNcVar(ncFd(iCur),thl0(2:i1,2:j1,1:k1))
          case(7); call writeNcVar(ncFd(iCur),qt0(2:i1,2:j1,1:k1))
          case(8); call writeNcVar(ncFd(iCur),ql0(2:i1,2:j1,1:k1))
          case(9); call writeNcVar(ncFd(iCur),e120(2:i1,2:j1,1:k1))
          case(10); call writeNcVar(ncFd(iCur),sv0(2:i1,2:j1,1:k1,iqr))
          case(11); call writeNcVar(ncFd(iCur),sv0(2:i1,2:j1,1:k1,iNr))
        end select
      end do

      call closeNcFile
      
      call system_clock(t2, count_rate, count_max)

      if (myid==0) write(*,*) 'ELAPSED TIME TOT = ',real(t2-t1)/real(count_rate)

      ! Determine the time for the next fielddump
      tNextDump = tNextDump + iInterval
      if (tNextDump > endhr(iNextHr) ) then
        ! In this case, there is no dump at endhr(iNextHr)
        iNextHr = iNextHr + 1
        if (iNextHr <= nhours) then
          ! Not yet reached the last time interval from namoptions yet
          tNextDump = beghr(iNextHr)
        else
          ! No more times present in beghr array.
          tNextDump = timeleft
          lfielddump = .false.
        end if
      end if
      dt_lim = min(dt_lim,tNextDump-timee)

    end if
 
  end subroutine fielddump
!> Clean up when leaving the run
  subroutine exitfielddump
    use modstat_nc, only : exitstat_nc,lnetcdf
    implicit none

    if(.not.(lfielddump)) return

!    if(lfielddump .and. lnetcdf) call exitstat_nc(ncid)
  end subroutine exitfielddump

end module modfielddump
