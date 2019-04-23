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

  use modglobal, only : longint

implicit none
private
PUBLIC :: initfielddump, fielddump,exitfielddump
save
!NetCDF variables
  integer,parameter :: nvar = 13
  integer :: ncid,nrec = 0
  character(80) :: fname = 'fielddump.xxx.xxx.xxx.nc'
  character(80),dimension(nvar,4) :: ncname
  character(80),dimension(1,4) :: tncname

  real    :: dtav
  integer(kind=longint) :: idtav,tnext
  integer :: klow,khigh
  logical :: lfielddump= .false. !< switch to enable the fielddump (on/off)

contains
!> Initializing fielddump. Read out the namelist, initializing the variables
  subroutine initfielddump
    use modmpi,   only :myid,my_real,comm3d,mpi_logical,mpi_integer,myidx,myidy
    use modglobal,only :imax,jmax,kmax,cexpnr,ifnamopt,fname_options,dtmax,dtav_glob, ladaptive,dt_lim,btime,tres
    use modstat_nc,only : lnetcdf,open_nc, define_nc,ncinfo,writestat_dims_nc
    implicit none
    integer :: ierr


    namelist/NAMFIELDDUMP/ &
    dtav,lfielddump,klow,khigh

    dtav=dtav_glob
    klow=1
    khigh=kmax
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
    call MPI_BCAST(klow        ,1,MPI_INTEGER,0,comm3d,ierr)
    call MPI_BCAST(khigh       ,1,MPI_INTEGER,0,comm3d,ierr)
    call MPI_BCAST(dtav        ,1,MY_REAL   ,0,comm3d,ierr)
    call MPI_BCAST(lfielddump  ,1,MPI_LOGICAL,0,comm3d,ierr)
    idtav = dtav/tres

    tnext      = idtav   +btime
    if(.not.(lfielddump)) return
    dt_lim = min(dt_lim,tnext)

    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'dtav should be a integer multiple of dtmax'
    end if
    if (lnetcdf) then
      write(fname,'(A,i3.3,A,i3.3,A)') 'fielddump.', myidx, '.', myidy, '.xxx.nc'
      fname(19:21) = cexpnr
      call ncinfo(tncname(1,:),'time','Time','s','time')
      call ncinfo(ncname( 1,:),'ua','West-East velocity','m/s','mttt')
      call ncinfo(ncname( 2,:),'va','South-North velocity','m/s','tmtt')
      call ncinfo(ncname( 3,:),'wa','Vertical velocity','m/s','ttmt')
      call ncinfo(ncname( 4,:),'cli','mass fraction of cloud ice','kg/kg','tttt')  
      call ncinfo(ncname( 5,:),'clw','mass fraction of cloud liquid water','kg/kg','tttt')
      call ncinfo(ncname( 6,:),'ta','air temperature ','K','tttt') 
      call ncinfo(ncname( 7,:),'plw','mass fraction of precipitating liquid water','kg/kg','tttt')
      call ncinfo(ncname( 8,:),'pli','mass fraction of precipitating ice','kg/kg','tttt')   
      call ncinfo(ncname( 9,:),'hus','specific humidity','kg/kg','tttt')   
      call ncinfo(ncname( 10,:),'hur','relative humidity','%','tttt')   
      call ncinfo(ncname( 11,:),'tntr','tendency of air temperature due to radiative heating','K/s','tttt')  
      call ncinfo(ncname( 12,:),'tntrs','tendency of air temperature due to shortwave radiative heating','K/s','tttt')  
      call ncinfo(ncname( 13,:),'tntrl','tendency of air temperature due to longwave radiative heating','K/s','tttt')  


      call open_nc(fname,  ncid,nrec,n1=imax,n2=jmax,n3=khigh-klow+1)
      if (nrec==0) then
        call define_nc( ncid, 1, tncname)
        call writestat_dims_nc(ncid)
      end if
     call define_nc( ncid, NVar, ncname)
    end if

  end subroutine initfielddump

!> Do fielddump. Collect data to truncated (2 byte) integers, and write them to file
  subroutine fielddump
    use modfields, only : u0,v0,w0,thl0,qt0,ql0,sv0,thv0h,thvh,&
                          tmp0,qri ,qsat0,rhof,exnf
    use modsurfdata,only : thls,qts,thvs
    use modglobal, only : imax,i1,ih,jmax,j1,jh,kmax,k1,rk3step,&
                          timee,dt_lim,cexpnr,ifoutput,rtimee,&
                           tup,tdn,cp,dzf
    use modmpi,    only : myid,cmyidx, cmyidy
    use modstat_nc, only : lnetcdf, writestat_nc
    use modmicrodata, only : iqr, imicro, imicro_none
    use modraddata, only   :lwu,lwd,swu,swd
    !rce use modmsebudg, only   :msebudg
    implicit none

    real, allocatable :: field(:,:,:),vars(:,:,:,:)
    real ::  il_ratio_
    integer i,j,k
    integer :: writecounter = 1
    integer :: reclength


    if (.not. lfielddump) return
    if (rk3step/=3) return

    if(timee<tnext) then
      dt_lim = min(dt_lim,tnext-timee)
      return
    end if

    tnext = tnext+idtav
    dt_lim = minval((/dt_lim,tnext-timee/))

    allocate(field(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(vars(imax,jmax,khigh-klow+1,nvar))

    reclength = imax*jmax*(khigh-klow+1)*2

    write (6,*) 'go to ncdf'
    if (lnetcdf) then
      field = u0
      vars(:,:,:,1) = field(2:i1,2:j1,klow:khigh)
      write (6,*) 'var 1'

      field = v0
      vars(:,:,:,2) = field(2:i1,2:j1,klow:khigh)
      write (6,*) 'var 2'

      field = w0
      vars(:,:,:,3) = field(2:i1,2:j1,klow:khigh)
      write (6,*) 'var 3'

      do i=2,i1 
      do j=2,j1
      do k=1,k1
       il_ratio_=max(0.,min(1.,(tmp0(i,j,k)-tdn)/(tup-tdn)))! cloud water vs cloud ice partitioning  
       field(i,j,k) = ql0(i,j,k) * il_ratio_
      enddo
      enddo
      enddo
      vars(:,:,:,5) = field(2:i1,2:j1,klow:khigh)
      
      do i=2,i1
      do j=2,j1
      do k=1,k1
       il_ratio_=max(0.,min(1.,(tmp0(i,j,k)-tdn)/(tup-tdn)))! cloud water vs cloud ice partitioning  
       field(i,j,k) = ql0(i,j,k) * (1-il_ratio_)
      enddo
      enddo
      enddo
      vars(:,:,:,4) = field(2:i1,2:j1,klow:khigh)

      write (6,*) 'var 4,5'

      field = tmp0
      vars(:,:,:,6) = field(2:i1,2:j1,klow:khigh)

      field = qri
      vars(:,:,:,8) = field(2:i1,2:j1,klow:khigh)
      write (6,*) 'var 6,8'

      if(imicro/=imicro_none) then
        do i=2,i1
        do j=2,j1
        do k=1,k1
          field(i,j,k) = sv0(i,j,k,iqr)
        enddo
        enddo
        enddo
      else
        field = 0.
      endif
      vars(:,:,:,7) = field(2:i1,2:j1,klow:khigh)
      write (6,*) 'var 7'

      field = qt0 - ql0 
      vars(:,:,:,9)  = field(2:i1,2:j1,klow:khigh) 
 
      do i=2,i1
      do j=2,j1
      do k=1,k1
       field(i,j,k) = (qt0(i,j,k)--ql0(i,j,k))/qsat0(i,j,k)
      enddo
      enddo
      enddo

      vars(:,:,:,10) = field (2:i1,2:j1,klow:khigh)
      write (6,*) 'var 9,10'

      do i=2,i1
      do j=2,j1
      do k=1,k1
        field(i,j,k) = (-lwd(i,j,k+1) - lwu(i,j,k+1) + lwd(i,j,k) + lwu(i,j,k))/(rhof(k)*exnf(k)*cp*dzf(k))
      end do
      end do
      end do
      vars(:,:,:,13)  = field(2:i1,2:j1,klow:khigh)

      do i=2,i1
      do j=2,j1
      do k=1,k1
        field(i,j,k) = (-swd(i,j,k+1) - swu(i,j,k+1) + swd(i,j,k) + swu(i,j,k))/(rhof(k)*exnf(k)*cp*dzf(k))
      end do
      end do
      end do
      vars(:,:,:,12)  = field(2:i1,2:j1,klow:khigh)

      vars(:,:,:,11) =  vars(:,:,:,13) +  vars(:,:,:,12)

      write (6,*) 'var 11,12,13'

      call writestat_nc(ncid,1,tncname,(/rtimee/),nrec,.true.)
      write (6,*) 'time writing done'
      call writestat_nc(ncid,nvar,ncname,vars,nrec,imax,jmax,khigh-klow+1)
      write (6,*) 'writing done'
    end if

    !rce call msebudg


    writecounter=writecounter+1

    deallocate(field,vars)

  end subroutine fielddump
!> Clean up when leaving the run
  subroutine exitfielddump
    use modstat_nc, only : exitstat_nc,lnetcdf
    implicit none

    if(lfielddump .and. lnetcdf) call exitstat_nc(ncid)
  end subroutine exitfielddump

end module modfielddump
