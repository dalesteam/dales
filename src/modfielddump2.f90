!> \file modfielddump2.f90
!!  Dumps 2D fields of several variables

!>
!!  Dumps 2D fields of several variables
!>
!!  Dumps 2D fields of several variables Written to wb*.myidx.myidy.expnr
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
module modfielddump2

  use modglobal, only : longint, &
                        Nvars_2D

implicit none
private
PUBLIC :: initfielddump2, fielddump2,exitfielddump2
save
!NetCDF variables
  integer,parameter :: nvar = 23  != Nvars_2D
  integer,parameter :: nmsevar = 5 
  integer :: ncid2,ncid3,nrec2 = 0, nrec3 = 0
  character(80) :: fname = 'fielddump2.xxx.xxx.xxx.nc'
  character(80),dimension(nvar,4) :: ncname
  character(80),dimension(1,4) :: tncname
  character(80),dimension(nmsevar,4) :: ncmsename
  character(80),dimension(1,4) :: tmsencname

  real    :: dtav
  integer(kind=longint) :: idtav,tnext
  logical :: lfielddump= .false. !< switch to enable the fielddump (on/off)

contains
!> Initializing fielddump. Read out the namelist, initializing the variables
  subroutine initfielddump2
    use modmpi,   only :myid,my_real,comm3d,mpi_logical,mpi_integer,myidx,myidy
    use modglobal,only :imax,jmax,kmax,cexpnr,ifnamopt,fname_options,dtmax,dtav_glob,kmax, ladaptive,dt_lim,btime,tres
    use modstat_nc,only : lnetcdf,open_nc, define_nc,ncinfo,writestat_dims_nc
    implicit none
    integer :: ierr

    namelist/NAMFIELDDUMP2/ &
    dtav,lfielddump

    dtav=dtav_glob
    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMFIELDDUMP2,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions NAMFIELDDUMP2'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMFIELDDUMP2'
      endif
      write(6 ,NAMFIELDDUMP2)
      close(ifnamopt)
    end if
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
      write(fname,'(A,i3.3,A,i3.3,A)') 'fielddump2.', myidx, '.', myidy, '.xxx.nc'    !rce table 5 2D hourly averaged variables
      fname(20:22) = cexpnr
      call ncinfo(tncname(1,:),'time','Time','s','time')
      call ncinfo(ncname( 1,:),'hfls','surface upward latent heat flux','W/m2','ttmt')
      call ncinfo(ncname( 2,:),'hfss','surface upward sensible heat flux','W/m2','ttmt')
      call ncinfo(ncname( 3,:),'rlds','surface downwellling longwave flux','W/m2','ttmt')
      call ncinfo(ncname( 4,:),'rlus','surface upwelling longwave flux','W/m2','ttmt')
      call ncinfo(ncname( 5,:),'rsds','surface downwellling shortwave flux','W/m2','ttmt')
      call ncinfo(ncname( 6,:),'rsus','surface upwellling shortwave flux','W/m2','ttmt')
      call ncinfo(ncname( 7,:),'rsdt','TOA incoming shortwave flux','units','ttmt')
      call ncinfo(ncname( 8,:),'rsut','TOA outgoing shortwave flux','units','ttmt')
      call ncinfo(ncname( 9,:),'rlut','TOA outgoing longwave flux','units','ttmt')
      call ncinfo(ncname(10,:),'rsdscs','surface downwelling shortwave flux - clear sky','W/m2','ttmt')
      call ncinfo(ncname(11,:),'rsuscs','surface upwelling shortwave flux - clear sky','W/m2','ttmt')
      call ncinfo(ncname(12,:),'rldscs','surface downwelling longwave flux - clear sky','W/m2','ttmt')
      call ncinfo(ncname(13,:),'rluscs','surface upwelling longwave flux - clear sky','W/m2','ttmt')
      call ncinfo(ncname(14,:),'rsutcs','TOA outgoing shortwave flux - clear sky','W/m2','ttmt')
      call ncinfo(ncname(15,:),'rlutcs','TOA outgoing longwave flux - clear sky','W/m2','ttmt')
      call ncinfo(ncname(16,:),'prw','water vapor path','kg/m2','ttmt')
      call ncinfo(ncname(17,:),'clwvi','condensed water path','kg/m2','ttmt')
      call ncinfo(ncname(18,:),'clivi','ice water path','kg/m2','ttmt')
      call ncinfo(ncname(19,:),'spwr','saturated water vapor path','kg/m2','ttmt')
      call ncinfo(ncname(20,:),'uabot','eastward wind at lowest model level','m/s','mttt')
      call ncinfo(ncname(21,:),'vabot','northward wind at lowest model level','m/s','tmtt')
      call ncinfo(ncname(22,:),'tabot','air temperature at lowest model level','K','tmtt')
      call ncinfo(ncname(23,:),'pr', 'surface precipitation rate','kg/m2s','ttmt')

      call open_nc(fname,  ncid2,nrec2,n1=imax,n2=jmax,n3=1)
      if (nrec2==0) then
        call define_nc( ncid2, 1, tncname)
        call writestat_dims_nc(ncid2)
      end if
     call define_nc( ncid2, NVar, ncname)
     write (6,*) 'finish 2D nc  initialization' !rce

     !----------               
     write(fname,'(A,i3.3,A,i3.3,A)') 'mse_budget.', myidx, '.', myidy, '.xxx.nc'    !rce table 6 2D hourly averaged variables
      fname(20:22) = cexpnr
      call ncinfo(tmsencname(1,:),'time','Time','s','time')
      call ncinfo(ncmsename( 1,:),'fmse'     ,'mass-weighted vertical integral of frozen moist static energy','J/m2','ttmt')
      call ncinfo(ncmsename( 2,:),'hadvmse'  ,'mass-weighted vertical integral of horizontal advective tendency of frozen moist static energy','J/m2/','ttmt')
      call ncinfo(ncmsename( 3,:),'vadvmse'  ,'mass-weighted vertical integral of vertical advective tendency of frozen moist static energy','J/m2/s','ttmt')
      call ncinfo(ncmsename( 4,:),'tnfmse'   ,'total tendency of mass-weighted vertical integral of frozen moist static energy','J/m2/s','ttmt')
      call ncinfo(ncmsename( 5,:),'tnfmsevar','total tendency of spatial variance of mass-weighted vertical integral of frozen moist static energy','J2/m4/s','ttmt')
      !rce write (6,*) 'done ncinfo'

      call open_nc(fname,  ncid3,nrec3,n1=imax,n2=jmax,n3=1)
      !rce write (6,*) 'done open_nc'
      if (nrec3==0) then
        call define_nc( ncid3, 1, tmsencname)
        !rce write (6,*) 'done define_nc' 
        call writestat_dims_nc(ncid3)
        !rece write (6,*) 'done writestat'
      end if
     call define_nc( ncid3, nmsevar, ncmsename)
     write (6,*) 'finish nc budget initialization' !rce

    end if
  end subroutine initfielddump2

!> Do fielddump. Collect data to truncated (2 byte) integers, and write them to file
  subroutine fielddump2
    use modfields, only : field_2D_mn,field_mse_2D
    use modglobal, only : imax,i1,ih,jmax,j1,jh,k1,rk3step,&
                          timee,dt_lim,cexpnr,ifoutput,rtimee,dt
    use modmpi,    only : myid,cmyidx, cmyidy
    use modstat_nc, only : lnetcdf, writestat_nc
    use modmicrodata, only : iqr, imicro, imicro_none
    use modmsebudg, only: msebudg
    implicit none

    real, allocatable :: vars(:,:,:,:)
    integer i,j,k,ivar
    integer :: writecounter = 1
    integer :: reclength


    if (.not. lfielddump) return
    if (rk3step/=3) return
    write (6,*) 'times ',dt,rtimee

    if(timee<tnext) then
      dt_lim = min(dt_lim,tnext-timee)
      return
    end if

    tnext = tnext+idtav
    dt_lim = minval((/dt_lim,tnext-timee/))

    allocate(vars(imax,jmax,1,nvar))

    reclength = imax*jmax*2

    vars(:,:,:,:) = -999.99
    do ivar=1,nvar
       vars(1:imax,1:jmax,1,ivar) = field_2D_mn (2:i1,2:j1,ivar)
    enddo

    if(lnetcdf) then
      write (6,*) 'ncid2 ',ncid2,rtimee
      call writestat_nc(ncid2,1,tncname,(/rtimee/),nrec2,.true.)
       write (6,*) 'do writestat ncid2 rtimee',rtimee
      call writestat_nc(ncid2,nvar,ncname,vars,nrec2,imax,jmax,1)
    end if

    

    write(6,*) 'do msebudg'
    if (lnetcdf) then
       call msebudg
       do ivar=1,nmsevar
         vars(1:imax,1:jmax,1,ivar) = field_mse_2D (2:i1,2:j1,ivar)
       enddo
       !write (6,*) 'ncid3 ',ncid3,rtimee
       write (6,*) 'ja ja' ,vars(:,:,1,4)
       call writestat_nc(ncid3,1,tmsencname,(/rtimee/),nrec3,.true.)
       !call writestat_nc(ncid3,nmsevar,ncmsename,field_mse_2D,nrec,imax,jmax,1)
       call writestat_nc(ncid3,nmsevar,ncmsename,vars,nrec3,imax,jmax,1)
    end if
    
    write (6,*) 'done msebudg'

    writecounter=writecounter+1

    field_2D_mn = 0.

    deallocate(vars)

  end subroutine fielddump2
!> Clean up when leaving the run
  subroutine exitfielddump2
    use modstat_nc, only : exitstat_nc,lnetcdf
    implicit none

    if(lfielddump .and. lnetcdf) call exitstat_nc(ncid2)
  end subroutine exitfielddump2

end module modfielddump2
