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

  use modglobal, only : longint, nsv

implicit none
private
PUBLIC :: initfielddump, fielddump,exitfielddump
save
!NetCDF variables
  integer :: nvar = 7
  integer :: ncid,nrec = 0
  character(80) :: fname = 'fielddump.xxx.xxx.xxx.nc'
  character(80),dimension(:,:), allocatable :: ncname
  character(80),dimension(1,4) :: tncname

  real    :: dtav, tmin, tmax
  integer(kind=longint) :: idtav,tnext,itmax,itmin
  integer :: klow,khigh,ncoarse=-1
  logical :: lfielddump= .false. !< switch to enable the fielddump (on/off)
  logical :: ldiracc   = .false. !< switch for doing direct access writing (on/off)
  logical :: lbinary   = .false. !< switch for doing direct access writing (on/off)

contains
!> Initializing fielddump. Read out the namelist, initializing the variables
  subroutine initfielddump
    use modmpi,   only :myid,comm3d,mpi_logical,mpi_integer,myidx,myidy &
                       , D_MPI_BCAST
    use modglobal,only :imax,jmax,kmax,cexpnr,ifnamopt,fname_options,dtmax,dtav_glob,kmax, ladaptive,dt_lim,btime,tres,checknamelisterror
    use modstat_nc,only : lnetcdf,open_nc, define_nc,ncinfo,writestat_dims_nc
    implicit none
    integer :: ierr, n
    character(3) :: csvname


    namelist/NAMFIELDDUMP/ &
    dtav,lfielddump,ldiracc,lbinary,klow,khigh,ncoarse, tmin, tmax

    dtav=dtav_glob
    klow=1
    khigh=kmax
    tmin = 0. 
    tmax = 1e8
    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMFIELDDUMP,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMFIELDDUMP')
      write(6 ,NAMFIELDDUMP)
      close(ifnamopt)
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
    if (ncoarse==-1) then
      ncoarse = 1
    end if
    idtav = dtav/tres
    itmin = tmin/tres
    itmax = tmax/tres

    tnext      = idtav   +btime
    if(.not.(lfielddump)) return
    dt_lim = min(dt_lim,tnext)

    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'dtav should be a integer multiple of dtmax'
    end if

    nvar = nvar + nsv
    if (lnetcdf) then
      write(fname,'(A,i3.3,A,i3.3,A)') 'fielddump.', myidx, '.', myidy, '.xxx.nc'
      fname(19:21) = cexpnr
      allocate(ncname(nvar,4))
      call ncinfo(tncname(1,:),'time','Time','s','time')
      call ncinfo(ncname( 1,:),'u','West-East velocity','m/s','mttt')
      call ncinfo(ncname( 2,:),'v','South-North velocity','m/s','tmtt')
      call ncinfo(ncname( 3,:),'w','Vertical velocity','m/s','ttmt')
      call ncinfo(ncname( 4,:),'qt','Total water specific humidity','1e-5kg/kg','tttt')
      call ncinfo(ncname( 5,:),'ql','Liquid water specific humidity','1e-5kg/kg','tttt')
      call ncinfo(ncname( 6,:),'thl','Liquid water potential temperature above 300K','K','tttt')
!       call ncinfo(ncname( 7,:),'qr','Rain water mixing ratio','1e-5kg/kg','tttt')
      call ncinfo(ncname( 7,:),'buoy','Buoyancy','K','tttt')
      do n=1,nsv
        write (csvname(1:3),'(i3.3)') n
        call ncinfo(ncname(7+n,:),'sv'//csvname,'Scalar '//csvname//' specific concentration','(kg/kg)','tttt')
      end do
      call open_nc(fname,  ncid,nrec,n1=ceiling(1.0*imax/ncoarse),n2=ceiling(1.0*jmax/ncoarse),n3=khigh-klow+1)
      if (nrec==0) then
        call define_nc( ncid, 1, tncname)
        call writestat_dims_nc(ncid, ncoarse)
      end if
     call define_nc( ncid, NVar, ncname)
    end if

  end subroutine initfielddump

!> Do fielddump. Collect data to truncated (2 byte) integers, and write them to file
  subroutine fielddump
    use modfields, only : u0,v0,w0,thl0,qt0,ql0,sv0,thv0h,thvh
    use modsurfdata,only : thls,qts,thvs
    use modglobal, only : imax,i1,ih,jmax,j1,jh,k1,rk3step,&
                          timee,dt_lim,cexpnr,ifoutput,rtimee
    use modmpi,    only : myid,cmyidx, cmyidy
    use modstat_nc, only : lnetcdf, writestat_nc
    use modmicrodata, only : iqr, imicro, imicro_none
    implicit none

    integer(KIND=selected_int_kind(4)), allocatable :: field(:,:,:)
    real, allocatable :: vars(:,:,:,:)
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
    allocate(vars(ceiling(1.0*imax/ncoarse),ceiling(1.0*jmax/ncoarse),khigh-klow+1,nvar))

    reclength = ceiling(1.0*imax/ncoarse)*ceiling(1.0*jmax/ncoarse)*(khigh-klow+1)*2

    field = NINT(1.0E3*u0,2)
    if (lnetcdf) vars(:,:,:,1) = u0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbuu.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
      else
        open  (ifoutput,file='wbuu.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1, ncoarse),j=2,j1, ncoarse),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    field = NINT(1.0E3*v0,2)
    if (lnetcdf) vars(:,:,:,2) = v0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbvv.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
      else
        open  (ifoutput,file='wbvv.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1, ncoarse),j=2,j1, ncoarse),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    field = NINT(1.0E3*w0,2)
    if (lnetcdf) vars(:,:,:,3) = w0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbww.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
      else
        open  (ifoutput,file='wbww.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1, ncoarse),j=2,j1, ncoarse),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    field = NINT(1.0E5*qt0,2)
    if (lnetcdf) vars(:,:,:,4) = qt0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqt.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
      else
        open  (ifoutput,file='wbqt.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1, ncoarse),j=2,j1, ncoarse),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    field = NINT(1.0E5*ql0,2)
    if (lnetcdf) vars(:,:,:,5) = ql0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbql.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
      else
        open  (ifoutput,file='wbql.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1, ncoarse),j=2,j1, ncoarse),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    field = NINT(1.0E2*(thl0-300),2)
    if (lnetcdf) vars(:,:,:,6) = thl0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbthl.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
      else
        open  (ifoutput,file='wbthl.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1, ncoarse),j=2,j1, ncoarse),k=klow,khigh)
      end if
      close (ifoutput)
    end if

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1.0E5*sv0(i,j,k,iqr),2)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1, ncoarse),j=2,j1, ncoarse),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    field=0.
    do i=2-ih,i1+ih
    do j=2-jh,j1+jh
    do k=2,k1
      field(i,j,k) = NINT(1.0E2*(thv0h(i,j,k)-thvh(k)),2)
    enddo
    enddo
    enddo
    
    if (lnetcdf) then 
      vars(:,:,:,7) = thv0h(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
      do k=klow,khigh
        vars(:,:,k,7) = vars(:,:,k,7) - thvh(k)
      end do
    end if
    do i=2-ih,i1+ih, ncoarse
    do j=2-jh,j1+jh, ncoarse
    do k=2,k1
      field(i,j,k) = NINT(1.0E2*(thv0h(i,j,k)-thvh(k)),2)
    enddo
    enddo
    enddo

    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbthv.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
      else
        open  (ifoutput,file='wbthv.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1, ncoarse),j=2,j1, ncoarse),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    if (lnetcdf) vars(:,:,:,8:nvar) = sv0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh,:)

    if(lnetcdf) then
      call writestat_nc(ncid,1,tncname,(/rtimee/),nrec,.true.)
      call writestat_nc(ncid,nvar,ncname,vars,nrec,ceiling(1.0*imax/ncoarse),ceiling(1.0*jmax/ncoarse),khigh-klow+1)
    end if

    if(lbinary) then
      if (myid==0) then
        open(ifoutput, file='wbthls.'//cexpnr,form='formatted',position='append')
        write(ifoutput,'(F12.1 3F12.5)') timee,thls, qts,thvs
        close(ifoutput)
      end if
    endif

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
