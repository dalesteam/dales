!> \file modtowerdump.f90
!!  Dumps 1D fields of several variables

!>
!!  Dumps 1D fields of several variables
!>
!!  Dumps 1D fields of several variables Written to wb*.myidx.myidy.expnr
!! If netcdf is true, this module leads the towerdump.myidx.myidy..expnr.nc output
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
module modtowerdump

  use modglobal, only : longint, nsv

implicit none
private
PUBLIC :: inittowerdump, towerdump,exittowerdump
save
!NetCDF variables
  integer :: nvar = 11
  integer :: ncid,nrec = 0
  character(80) :: fname = 'towerdump.xxx.nc'
  character(80),dimension(:,:), allocatable :: ncname
  character(80),dimension(1,4) :: tncname

  real    :: dtav, tmin, tmax
  integer(kind=longint) :: idtav,tnext,itmax,itmin
  integer :: klow,khigh
  logical :: ltowerdump= .false. !< switch to enable the towerdump (on/off)

contains
!> Initializing towerdump. Read out the namelist, initializing the variables
  subroutine inittowerdump
    use modmpi,   only :myid,my_real,comm3d,mpi_logical,mpi_integer,myidx,myidy
    use modglobal,only :imax,jmax,kmax,cexpnr,ifnamopt,fname_options,dtmax,dtav_glob,kmax, ladaptive,dt_lim,btime,tres
    use modstat_nc,only : lnetcdf,open_nc, define_nc,ncinfo,writestat_dims_nc
    implicit none
    integer :: ierr, n
    character(3) :: csvname


    namelist/NAMTOWERDUMP/ &
    dtav,ltowerdump,klow,khigh,tmin, tmax

    dtav=dtav_glob
    klow=1
    khigh=kmax
    tmin = 0. 
    tmax = 1e8
    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMTOWERDUMP,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions NAMTOWERDUMP'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMTOWERDUMP'
      endif
      write(6 ,NAMTOWERDUMP)
      close(ifnamopt)
    end if
    call MPI_BCAST(klow        ,1,MPI_INTEGER,0,comm3d,ierr)
    call MPI_BCAST(khigh       ,1,MPI_INTEGER,0,comm3d,ierr)
    call MPI_BCAST(dtav        ,1,MY_REAL   ,0,comm3d,ierr)
    call MPI_BCAST(tmin        ,1,MY_REAL   ,0,comm3d,ierr)
    call MPI_BCAST(tmax        ,1,MY_REAL   ,0,comm3d,ierr)
    call MPI_BCAST(ltowerdump  ,1,MPI_LOGICAL,0,comm3d,ierr)
    idtav = dtav/tres
    itmin = tmin/tres
    itmax = tmax/tres

    tnext      = idtav   +btime
    if(.not.(ltowerdump)) return
    dt_lim = min(dt_lim,tnext)

    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'dtav should be a integer multiple of dtmax'
    end if

    nvar = nvar + nsv
    if (myid==0) then
      if (lnetcdf) then
        write(fname,'(A)') 'towerdump.xxx.nc'
        fname(11:13) = cexpnr
        allocate(ncname(nvar,4))
        call ncinfo(tncname(1,:),'time',                              'Time',    's','time')
        call ncinfo(ncname( 1,:),  'um',         'West-East velocity upwind',  'm/s',  'tt')
        call ncinfo(ncname( 2,:),  'up',       'West-East velocity downwind',  'm/s',  'tt')
        call ncinfo(ncname( 3,:),  'vm',       'South-North velocity upwind',  'm/s',  'tt')
        call ncinfo(ncname( 4,:),  'vp',     'South-North velocity downwind',  'm/s',  'tt')
        call ncinfo(ncname( 5,:),   'w',                 'Vertical velocity',  'm/s',  'mt')
        call ncinfo(ncname( 6,:),  'qt',     'Total water specific humidity','kg/kg',  'tt')
        call ncinfo(ncname( 7,:),  'ql',    'Liquid water specific humidity','kg/kg',  'tt')
        call ncinfo(ncname( 8,:), 'thl','Liquid water potential temperature',    'K',  'tt')
        call ncinfo(ncname( 9,:),'thvh',    'Virt. pot. temp. at half level',    'K',  'mt')
        call ncinfo(ncname(10,:),   'T',              'Absolute temperature',    'K',  'tt')
        call ncinfo(ncname(11,:),   'P',                          'Pressure',   'Pa',  'tt')
        do n=1,nsv
          write (csvname(1:3),'(i3.3)') n
          call ncinfo(ncname(11+n,:),'sv'//csvname,'Scalar '//csvname//' mixing ratio','ppb','tt')
        end do
        call open_nc(fname,  ncid,nrec,n3=khigh-klow+1)
        if (nrec==0) then
          call define_nc( ncid, 1, tncname)
          call writestat_dims_nc(ncid)
        end if
       call define_nc( ncid, NVar, ncname)
      end if
    end if

  end subroutine inittowerdump

!> Do towerdump. Collect data to truncated (2 byte) integers, and write them to file
  subroutine towerdump
    use modfields, only : u0,v0,w0,thl0,qt0,ql0,svm,thv0h,exnf, presf
    use modsurfdata,only : thls,qts,thvs
    use modglobal, only : imax,i1,ih,jmax,j1,jh,k1,rk3step,&
                          timee,dt_lim,cexpnr,ifoutput,rtimee,cp,rlv
    use modmpi,    only : myid,cmyidx, cmyidy
    use modstat_nc, only : lnetcdf, writestat_nc
    use modmicrodata, only : iqr, imicro, imicro_none
    implicit none

    real, allocatable :: vars(:,:)
    integer :: writecounter = 1


    if (.not. ltowerdump) return
    if (rk3step/=3) return

    if(timee<tnext) then
      dt_lim = min(dt_lim,tnext-timee)
      return
    end if

    tnext = tnext+idtav
    dt_lim = minval((/dt_lim,tnext-timee/))

    if (myid==0) then

      allocate(vars(khigh-klow+1,nvar))

      if (lnetcdf) then
         vars(:,      1) =    u0(2,2,klow:khigh)
         vars(:,      2) =    u0(3,2,klow:khigh)
         vars(:,      3) =    v0(2,2,klow:khigh)
         vars(:,      4) =    v0(2,3,klow:khigh)
         vars(:,      5) =    w0(2,2,klow:khigh)
         vars(:,      6) =   qt0(2,2,klow:khigh)
         vars(:,      7) =   ql0(2,2,klow:khigh)
         vars(:,      8) =  thl0(2,2,klow:khigh)
         vars(:,      9) = thv0h(2,2,klow:khigh)
         vars(:,     10) =  thl0(2,2,klow:khigh) * exnf(klow:khigh) + (rlv/cp) * ql0(2,2,klow:khigh) 
         vars(:,     11) = presf(    klow:khigh)
         vars(:,12:nvar) =   svm(2,2,klow:khigh,:)
      endif

      if(lnetcdf) then
        call writestat_nc(ncid,1,tncname,(/rtimee/),nrec,.true.)
        call writestat_nc(ncid,nvar,ncname,vars,nrec,khigh-klow+1)
      end if

      writecounter=writecounter+1

      deallocate(vars)

    end if

  end subroutine towerdump
!> Clean up when leaving the run
  subroutine exittowerdump
    use modstat_nc, only : exitstat_nc,lnetcdf
    use modmpi,    only : myid
    implicit none

    if(ltowerdump .and. (lnetcdf .and. (myid==0))) call exitstat_nc(ncid)
  end subroutine exittowerdump

end module modtowerdump
