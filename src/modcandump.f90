!> \file modcandump.f90
!!  Dumps 3D fields of canopy variables

!>
!!  Dumps 3D fields of canopy variables
!>
!!  Dumps 3D fields of several variables Written to wb*.myidx.myidy.expnr
!! If netcdf is true, this module leads the candump.myidx.myidy..expnr.nc output
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
module modcandump

  use modglobal, only : longint, nsv
  use modcanopy, only : ncanopy

implicit none
private
PUBLIC :: initcandump, candump,exitcandump
save
!NetCDF variables
  integer :: nvar = 11
  integer :: ncid,nrec = 0
  character(80) :: fname = 'candump.xxx.xxx.xxx.nc'
  character(80),dimension(:,:), allocatable :: ncname
  character(80),dimension(1,4) :: tncname

  real    :: dtav, tmin, tmax
  integer(kind=longint) :: idtav,tnext,itmax,itmin
  integer :: ncoarse=-1
  integer :: khigh= 10 
  integer :: klow = 1
  logical :: lcandump= .false. !< switch to enable the candump (on/off)
  logical :: ldiracc   = .false. !< switch for doing direct access writing (on/off)
  logical :: lbinary   = .false. !< switch for doing direct access writing (on/off)

contains
!> Initializing candump. Read out the namelist, initializing the variables
  subroutine initcandump
    use modmpi,   only :myid,my_real,comm3d,mpi_logical,mpi_integer,myidx,myidy
    use modglobal,only :imax,jmax,kmax,cexpnr,ifnamopt,fname_options,dtmax,dtav_glob,kmax, ladaptive,dt_lim,btime,tres
    use modstat_nc,only : lnetcdf,open_nc, define_nc,ncinfo,writestat_dims_nc
    use modcanopy, only : ncanopy
    implicit none
    integer :: ierr, n


    namelist/NAMCANDUMP/ &
    dtav,lcandump,ldiracc,lbinary,khigh,klow,ncoarse,tmin,tmax

    dtav  = dtav_glob
    klow  = 1
    khigh = ncanopy
    tmin  = 0. 
    tmax  = 1e8
    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMCANDUMP,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions NAMCANDUMP'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMCANDUMP'
      endif
      write(6 ,NAMCANDUMP)
      close(ifnamopt)
    end if
    call MPI_BCAST(ncoarse     ,1,MPI_INTEGER,0,comm3d,ierr)
    call MPI_BCAST(klow        ,1,MPI_INTEGER,0,comm3d,ierr)
    call MPI_BCAST(khigh       ,1,MPI_INTEGER,0,comm3d,ierr)
    call MPI_BCAST(dtav        ,1,MY_REAL   ,0,comm3d,ierr)
    call MPI_BCAST(tmin        ,1,MY_REAL   ,0,comm3d,ierr)
    call MPI_BCAST(tmax        ,1,MY_REAL   ,0,comm3d,ierr)
    call MPI_BCAST(lcandump    ,1,MPI_LOGICAL,0,comm3d,ierr)
    call MPI_BCAST(ldiracc     ,1,MPI_LOGICAL,0,comm3d,ierr)
    call MPI_BCAST(lbinary     ,1,MPI_LOGICAL,0,comm3d,ierr)
    if (ncoarse==-1) then
      ncoarse = 1
    end if
    idtav = dtav/tres
    itmin = tmin/tres
    itmax = tmax/tres

    tnext      = idtav   +btime
    if(.not.(lcandump)) return
    dt_lim = min(dt_lim,tnext)

    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'dtav should be a integer multiple of dtmax'
    end if

    nvar = nvar 
    if (lnetcdf) then
      write(fname,'(A,i3.3,A,i3.3,A)') 'candump.', myidx, '.', myidy, '.xxx.nc'
      fname(17:19) = cexpnr
      allocate(ncname(nvar,4))
      call ncinfo(tncname(1,:), 'time','Time','s','time')
      call ncinfo(ncname( 1,:), 'sh_can','Sensible heat source','W/m3','tttt')
      call ncinfo(ncname( 2,:), 'le_can','Latent heat source','W/3','tttt')
      call ncinfo(ncname( 3,:), 'Fco2_can','CO2 source','mg CO2 m-3 s-1','tttt')
      call ncinfo(ncname( 4,:), 'qt','Total water specific humidity','kg/kg','tttt')
      call ncinfo(ncname( 5,:), 'thl','Liquid water potential temperature above 300K','K','tttt')
      call ncinfo(ncname( 6,:), 'buoy','Buoyancy','K','tttt')
      call ncinfo(ncname( 7,:), 'CO2','CO2 concentration','ppb','tttt')
      call ncinfo(ncname( 8,:), 'gcc_leafsun', 'Carbon stomatal conductance on sunlit leaves','m s-1','tttt')
      call ncinfo(ncname( 9,:), 'gcc_leafshad','Carbon stomatal conductance on shaded leaves','m s-1','tttt')
      call ncinfo(ncname( 10,:),'ci_leafsun',  'Carbon internal concentration on sunlit leaves','mg m-3','tttt')
      call ncinfo(ncname( 11,:),'ci_leafshad', 'Carbon internal concentration on shaded leaves','mg m-3','tttt')
      !radiation vars still missing
      call open_nc(fname,  ncid,nrec,n1=ceiling(1.0*imax/ncoarse),n2=ceiling(1.0*jmax/ncoarse),n3=khigh-klow+1)
      if (nrec==0) then
        call define_nc( ncid, 1, tncname)
        call writestat_dims_nc(ncid, ncoarse)
      end if
     call define_nc( ncid, NVar, ncname)
    end if

  end subroutine initcandump

!> Do candump. Collect data to truncated (2 byte) integers, and write them to file
  subroutine candump
    use modfields, only : thl0,qt0,thv0h,thvh,sv0
    use modsurfdata,only : indCO2
    use modglobal, only : imax,i1,ih,jmax,j1,jh,k1,rk3step,&
                          timee,dt_lim,cexpnr,ifoutput,rtimee
    use modmpi,    only : myid,cmyidx, cmyidy
    use modstat_nc, only : lnetcdf, writestat_nc
    use modcanopy, only : sh_can,le_can,Fco2_can,gcc_leafsun,gcc_leafshad,ci_leafsun,ci_leafshad
    implicit none

    integer(KIND=selected_int_kind(4)), allocatable :: field(:,:,:)
    real, allocatable :: vars(:,:,:,:)
    integer i,j,k
    integer :: writecounter = 1
    integer :: reclength


    if (.not. lcandump) return
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

    field = NINT(1.0E3*sh_can,2)
    if (lnetcdf) vars(:,:,:,1) = sh_can(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbshcan.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
      else
        open  (ifoutput,file='wbshcan.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1, ncoarse),j=2,j1, ncoarse),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    field = NINT(1.0E3*le_can,2)
    if (lnetcdf) vars(:,:,:,2) = le_can(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wblecan.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
      else
        open  (ifoutput,file='wblecan.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1, ncoarse),j=2,j1, ncoarse),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    field = NINT(1.0E5*Fco2_can,2)
    if (lnetcdf) vars(:,:,:,3) = Fco2_can(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbfco2can.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
      else
        open  (ifoutput,file='wbfco2can.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
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

    field = NINT(1.0E2*(thl0-300),2)
    if (lnetcdf) vars(:,:,:,5) = thl0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
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


    field=0.
    do i=2-ih,i1+ih
    do j=2-jh,j1+jh
    do k=2,k1
      field(i,j,k) = NINT(1.0E2*(thv0h(i,j,k)-thvh(k)),2)
    enddo
    enddo
    enddo
    
    if (lnetcdf) then 
      vars(:,:,:,6) = thv0h(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
      do k=klow,khigh
        vars(:,:,k,6) = vars(:,:,k,6) - thvh(k)
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
    
    field = NINT(sv0(:,:,:,indCO2),2)
    if (lnetcdf) vars(:,:,:,7) = sv0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh,indCO2)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbco2.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
      else
        open  (ifoutput,file='wbco2.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1, ncoarse),j=2,j1, ncoarse),k=klow,khigh)
      end if
      close (ifoutput)
    endif
    
    field = NINT(1.0e5*gcc_leafsun(:,:,:),2)
    if (lnetcdf) vars(:,:,:,8) = gcc_leafsun(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbgcc_leafsun.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
      else
        open  (ifoutput,file='wbgcc_leafsun.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1, ncoarse),j=2,j1, ncoarse),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    field = NINT(1.0e5*gcc_leafshad(:,:,:),2)
    if (lnetcdf) vars(:,:,:,9) = gcc_leafshad(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbgcc_leafshad.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
      else
        open  (ifoutput,file='wbgcc_leafshad.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1, ncoarse),j=2,j1, ncoarse),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    field = NINT(1.0e2*ci_leafsun(:,:,:),2)
    if (lnetcdf) vars(:,:,:,10) = ci_leafsun(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbci_leafsun.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
      else
        open  (ifoutput,file='wbci_leafsun.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1, ncoarse),j=2,j1, ncoarse),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    field = NINT(1.0e2*ci_leafshad(:,:,:),2)
    if (lnetcdf) vars(:,:,:,11) = ci_leafshad(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbci_leafshad.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
      else
        open  (ifoutput,file='wbci_leafshad.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1, ncoarse),j=2,j1, ncoarse),k=klow,khigh)
      end if
      close (ifoutput)
    endif


    if(lnetcdf) then
      call writestat_nc(ncid,1,tncname,(/rtimee/),nrec,.true.)
      call writestat_nc(ncid,nvar,ncname,vars,nrec,ceiling(1.0*imax/ncoarse),ceiling(1.0*jmax/ncoarse),khigh-klow+1)
    end if

    writecounter=writecounter+1

    deallocate(field,vars)

  end subroutine candump
!> Clean up when leaving the run
  subroutine exitcandump
    use modstat_nc, only : exitstat_nc,lnetcdf
    implicit none

    if(lcandump .and. lnetcdf) call exitstat_nc(ncid)
  end subroutine exitcandump

end module modcandump
