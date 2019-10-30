!> \file modthreedheating.f90
!!  Dumps 3D fields of several variables

!>
!!  Dumps 3D fields of several variables
!>
!!  Dumps 3D fields of several variables Written to wb*.myidx.myidy.expnr
!! If netcdf is true, this module leads the threedheating.myidx.myidy..expnr.nc output
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
module modthreeDheating
#ifndef HAVE_3DHEATING
  contains
  subroutine initthreedheating
    stop
  end subroutine initthreedheating
#else
  use modglobal, only : longint

implicit none
private
PUBLIC :: initthreedheating, threedheating,exitthreedheating
save
!NetCDF variables
  integer,parameter :: nvar = 13
  integer :: ncid3,nrec = 0
  character(80) :: fname = 'threedheating.xxx.xxx.xxx.nc'
  character(80),dimension(nvar,4) :: ncname
  character(80),dimension(1,4) :: tncname

  real    :: dtav
  integer(kind=longint) :: idtav,tnext
  integer :: klow,khigh
  logical :: lthreedheating= .false. !< switch to enable the threedheating (on/off)
  logical :: ldiracc   = .false. !< switch for doing direct access writing (on/off)
  logical :: lbinary   = .false. !< switch for doing direct access writing (on/off)

contains
!> Initializing threedheating. Read out the namelist, initializing the variables
  subroutine initthreedheating
    use modmpi,   only :myid,my_real,comm3d,mpi_logical,mpi_integer,myidx,myidy
    use modglobal,only :imax,jmax,kmax,cexpnr,ifnamopt,fname_options,dtmax,dtav_glob,kmax, ladaptive,dt_lim,btime,tres
    use modstat_nc,only : lnetcdf,open_nc, define_nc,ncinfo,writestat_dims_nc
    implicit none
    integer :: ierr


    namelist/NAMthreedheating/ &
    dtav,lthreedheating,ldiracc,lbinary,klow,khigh

    dtav=dtav_glob
    klow=1
    khigh=kmax
    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMthreedheating,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions NAMthreedheating'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMthreedheating'
      endif
      write(6 ,NAMthreedheating)
      close(ifnamopt)
    end if
    call MPI_BCAST(klow        ,1,MPI_INTEGER,0,comm3d,ierr)
    call MPI_BCAST(khigh       ,1,MPI_INTEGER,0,comm3d,ierr)
    call MPI_BCAST(dtav        ,1,MY_REAL   ,0,comm3d,ierr)
    call MPI_BCAST(lthreedheating  ,1,MPI_LOGICAL,0,comm3d,ierr)
    call MPI_BCAST(ldiracc     ,1,MPI_LOGICAL,0,comm3d,ierr)
    call MPI_BCAST(lbinary     ,1,MPI_LOGICAL,0,comm3d,ierr)
    idtav = dtav/tres
    
    tnext      = idtav   +btime
    if(.not.(lthreedheating)) return
    dt_lim = min(dt_lim,tnext)

    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'dtav should be a integer multiple of dtmax'
    end if
    if (lnetcdf) then
      write(fname,'(A,i3.3,A,i3.3,A)') 'threedheating.', myidx, '.', myidy, '.xxx.nc'
      fname(23:25) = cexpnr
      call ncinfo(tncname(1,:),'time','Time','s','time')
      call ncinfo(ncname( 1,:),'dtrad','theta_l tendency due to total radiaton','K/d','tttt')
      call ncinfo(ncname( 2,:),'dtSW','theta_l tendency due to short wave radiaton','K/d','tttt')
      call ncinfo(ncname( 3,:),'dtLW','theta_l tendency due to long wave radiaton','K/d','tttt')
      call ncinfo(ncname( 4,:),'ql','Liquid water specific humidity','1e-5kg/kg','tttt')
      call ncinfo(ncname( 5,:),'qt','total specific humidity','1e-5kg/kg','tttt')
      call ncinfo(ncname( 6,:),'thl','Liquid water potential temperature','K','tttt')      
      call ncinfo(ncname( 7,:),'swd','shortwave radiatio downward','W/m2','tttt')  
      call ncinfo(ncname( 8,:),'swu','shortwave radiatio downward','W/m2','tttt')
      call ncinfo(ncname( 9,:),'swdir','direct shortwave radiatio downward','W/m2','tttt')
      call ncinfo(ncname( 10,:),'swdif','difuse shortwave radiatio downward','W/m2','tttt')
      call ncinfo(ncname( 11,:),'lwd','shortwave radiatio downward','W/m2','tttt')
      call ncinfo(ncname( 12,:),'lwu','shortwave radiatio downward','W/m2','tttt')
      call ncinfo(ncname( 13,:),'qlrad','liquid water specific humidity','kg/kg','tttt')
    

call open_nc(fname,  ncid3,nrec,n1=imax,n2=jmax,n3=khigh-klow+1)
      if (nrec==0) then
        call define_nc( ncid3, 1, tncname)
        call writestat_dims_nc(ncid3)
      end if
     call define_nc( ncid3, NVar, ncname)
    end if

  end subroutine initthreedheating

!> Do threedheating. Collect data to truncated (2 byte) integers, and write them to file
  subroutine threedheating
    use modraddata,only : thlprad,thlprLW,thlprSW,swd,swu,lwd,lwu,swdir,swdif,iradiation,irad_tenstr
    use modsurfdata,only : thls,qts,thvs
    use modfields,only : ql0,qt0,w0,thl0,qlrad
    use modglobal, only : imax,i1,ih,jmax,j1,jh,k1,rk3step,&
                          timee,dt_lim,cexpnr,ifoutput,rtimee
    use modmpi,    only : myid,cmyidx, cmyidy
    use modstat_nc, only : lnetcdf, writestat_nc
    use modmicrodata, only : iqr, imicro, imicro_none
    use modsampling, only: wqtthrad,wthlthrad
    implicit none

    !integer(KIND=selected_int_kind(4)), allocatable :: field(:,:,:),vars(:,:,:,:)
    real,allocatable :: field(:,:,:),vars(:,:,:,:)
    integer i,j,k
    integer :: writecounter = 1
    integer :: reclength
    
    
    if (.not. lthreedheating) return

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

    field = thlprad*86400.
    if (lnetcdf) vars(:,:,:,1) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbthlprad.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbthlprad.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif
    
    field = thlprSW*86400
    if (lnetcdf) vars(:,:,:,2) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='thlprSW.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='thlprSW.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
  endif
    field = thlprLW*86400
    if (lnetcdf) vars(:,:,:,3) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='thlprLW.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='thlprLW.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif
    field = ql0*1.0e5
    if (lnetcdf) vars(:,:,:,4) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqlql.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqlql.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif
    
    field = qt0*1.0e5
    if (lnetcdf) vars(:,:,:,5) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wb9swdswd.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wb9swdswd.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif
    field = thl0
    if (lnetcdf) vars(:,:,:,6) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wb8swdswd.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wb8swdswd.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif
    
    field = swd
    if (lnetcdf) vars(:,:,:,7) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wb7swdswd.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wb7swdswd.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif
    
 field = swu
    if (lnetcdf) vars(:,:,:,8) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wb1swdswd.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wb1swdswd.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif
  
   if (iradiation == irad_tenstr) then
    field = -1*swdir
   else
    field = swdir
   endif
    if (lnetcdf) vars(:,:,:,9) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wb7swdswdir.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wb7swdswdir.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif
   field = swdif
    if (lnetcdf) vars(:,:,:,10) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wb7swdswdif.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wb7swdswdif.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

 field = lwd
    if (lnetcdf) vars(:,:,:,11) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wb2swdswd.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wb2swdswd.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

 field = lwu
    if (lnetcdf) vars(:,:,:,12) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wb3swdswd.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wb3swdswd.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif
    
 field = qlrad
 if (lnetcdf) vars(:,:,:,13) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbrad3swdswd.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbrad3swdswd.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

   if(lnetcdf) then
      
      call writestat_nc(ncid3,1,tncname,(/rtimee/),nrec,.true.)
      call writestat_nc(ncid3,nvar,ncname,vars,nrec,imax,jmax,khigh-klow+1)
    end if
    writecounter=writecounter+1
    if(lbinary) then
      if (myid==0) then
        open(ifoutput, file='wbthls.'//cexpnr,form='formatted',position='append')
        write(ifoutput,'(F12.1 3F12.5)') timee,thls, qts,thvs
        close(ifoutput)
      end if
    endif
    deallocate(field,vars)

  end subroutine threedheating
!> Clean up when leaving the run
  subroutine exitthreedheating
    use modstat_nc, only : exitstat_nc,lnetcdf
    
    implicit none
   
    if(lthreedheating .and. lnetcdf) call exitstat_nc(ncid3)
  end subroutine exitthreedheating
#endif
end module modthreedheating
