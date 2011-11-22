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

  use modglobal, only : longint

implicit none
private
PUBLIC :: initfielddump,fielddump,exitfielddump
save
!NetCDF variables
  integer           :: ncid,nVar = 12           !< Number of variables (if bulkmicro->nVar+2)
  character(80)     :: fname = 'fielddump_hhmmss.xxx.nc'
  character(80),allocatable,dimension(:,:) :: ncfieldinfo
  character(80),dimension(1,4)             :: tncfieldinfo

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
    use modmpi,       only : myid,my_real,comm3d,mpi_logical,mpi_integer,cmyid
    use modglobal,    only : cexpnr,ifnamopt,fname_options,dtmax,btime,runtime,timee,tres,dt_lim, &
                             imax,jmax,kmax,ladaptive
    use modboundary,  only : ksp
    use modmicrodata, only : imicro,imicro_bulk
    use modstat_nc,   only : lnetcdf,open_nc,define_nc,redefine_nc,ncinfo,writestat_dims_nc
    implicit none
    integer :: ierr

    namelist/NAMFIELDDUMP/ &
    lfielddump,beghr,endhr,interval,kfieldtop

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

    if (kfieldtop<0 .or. kfieldtop>kmax) then
      if (myid == 0) then 
        write(*,*) 'WARNING: no or invalid kfieldtop given.'
        write(*,*) 'WARNING: kfieldtop set to just below sponge layer: ',ksp-1
      end if
      kfieldtop = ksp-1
    end if

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

      fname(18:20) = cexpnr

      if (imicro==imicro_bulk) nVar = nVar + 2
      allocate(ncfieldinfo(nVar,4))

      call ncinfo(tncfieldinfo(1,:),'time','Time','s','time')
      call ncinfo(ncfieldinfo( 1,:),'zfull','Height at full levels','m','z')
      call ncinfo(ncfieldinfo( 2,:),'zhalf','Height at half levels','m','zh')
      call ncinfo(ncfieldinfo( 3,:),'pfull','Pressure at full levels','Pa','z')
      call ncinfo(ncfieldinfo( 4,:),'phalf','Pressure at half levels','Pa','zh')
      call ncinfo(ncfieldinfo( 5,:),'u','West-East velocity','m/s','xhyhz')
      call ncinfo(ncfieldinfo( 6,:),'v','South-North velocity','m/s','xhyhz')
      call ncinfo(ncfieldinfo( 7,:),'w','Vertical velocity','m/s','xyzh')
      call ncinfo(ncfieldinfo( 8,:),'thl','Liquid water potential temperature','K','xyz')
      call ncinfo(ncfieldinfo( 9,:),'qt','Total water content','g/kg','xyz')
      call ncinfo(ncfieldinfo(10,:),'ql','Liquid water content','g/kg','xyz')
      call ncinfo(ncfieldinfo(11,:),'e', 'Subgrid scale turbulent kinetic energy','m^2/s^2','xyz')
      call ncinfo(ncfieldinfo(12,:),'p', 'Pressure fluctuation','Pa','xyz')
      
      if (imicro==imicro_bulk) then
        call ncinfo(ncfieldinfo(13,:),'qr','Rain water content','g/kg','xyz')
        call ncinfo(ncfieldinfo(14,:),'Nr','Rain droplet number density','cm^-3','xyz')
      end if

    else

      if (myid == 0) then
        write(*,*) 'WARNING: lfielddump = .true. while lnetcdf = .false.'
        write(*,*) 'WARNING: No fielddumps will be made. Setting lfielddump = .false.'
        lfielddump = .false.
      end if

    end if ! if lnetcdf
    
  end subroutine initfielddump

!> Do fielddump
  subroutine fielddump
    use modglobal,  only : timee,rtimee,dt_lim,tres,timeleft
    use modmpi,     only : myid
    use modstat_nc, only : lnetcdf,ncfielddump
    implicit none

    integer :: imin,ihour,isec

    if(.not.(lfielddump)) return
    
    if (timee<tNextDump) then
      dt_lim = min(dt_lim,tNextDump-timee)
    else ! if (timee>=tNextDump)
      ! It's time for a fielddump.
      ! Make a filename first,...
      ihour = floor(rtimee/3600)
      imin  = floor((rtimee-ihour*3600) /60.)
      isec  = floor((rtimee-ihour*3600-imin*60))
      write (fname(11:12),'(i2.2)') ihour
      write (fname(13:14),'(i2.2)') imin
      write (fname(15:16),'(i2.2)') isec
      if (myid==0) write(*,*) 'fname = ',fname
      
      ! ...then call the subroutine that makes the actual file and puts data in.
      call ncfielddump(fname,nVar,ncfieldinfo,kfieldtop)

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

    deallocate(ncfieldinfo)
!    if(lfielddump .and. lnetcdf) call exitstat_nc(ncid)
  end subroutine exitfielddump

end module modfielddump
