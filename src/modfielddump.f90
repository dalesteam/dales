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
  integer,parameter :: nvar = 45
  integer :: ncid,nrec = 0
  character(80) :: fname = 'fielddump.xxx.xxx.xxx.nc'
  character(80),dimension(nvar,4) :: ncname
  character(80),dimension(1,4) :: tncname

  real    :: dtav
  integer(kind=longint) :: idtav,tnext
  integer :: klow,khigh
  logical :: lfielddump= .false. !< switch to enable the fielddump (on/off)
  logical :: ldiracc   = .false. !< switch for doing direct access writing (on/off)
  logical :: lbinary   = .false. !< switch for doing direct access writing (on/off)

contains
!> Initializing fielddump. Read out the namelist, initializing the variables
  subroutine initfielddump
    use modmpi,   only :myid,my_real,comm3d,mpi_logical,mpi_integer,myidx,myidy
    use modglobal,only :imax,jmax,kmax,cexpnr,ifnamopt,fname_options,dtmax,dtav_glob,kmax, ladaptive,dt_lim,btime,tres
    use modstat_nc,only : lnetcdf,open_nc, define_nc,ncinfo,writestat_dims_nc
    implicit none
    integer :: ierr


    namelist/NAMFIELDDUMP/ &
    dtav,lfielddump,ldiracc,lbinary,klow,khigh

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
    call MPI_BCAST(ldiracc     ,1,MPI_LOGICAL,0,comm3d,ierr)
    call MPI_BCAST(lbinary     ,1,MPI_LOGICAL,0,comm3d,ierr)
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

      call ncinfo(ncname( 1,:),'u', 'West-East velocity',  'm/s', 'mttt')
      call ncinfo(ncname( 2,:),'v', 'South-North velocity','m/s', 'tmtt')
      call ncinfo(ncname( 3,:),'w', 'Vertical velocity',   'm/s', 'ttmt')

      call ncinfo(ncname( 4,:),'qt',  'Total water specific humidity',                 '1e-5kg/kg', 'tttt')
      call ncinfo(ncname( 5,:),'ql',  'Liquid water specific humidity',                '1e-6kg/kg', 'tttt')
      call ncinfo(ncname( 6,:),'thl', 'Liquid water potential temperature above 300K', 'K',         'tttt')
      call ncinfo(ncname( 7,:),'qr',  'Rain water specific humidity',                  '1e-6kg/kg', 'tttt')
      call ncinfo(ncname( 8,:),'buoy','Buoyancy',                                      'K',         'tttt')

      call ncinfo(ncname( 9,:),'nusn','Free aerosol number concentration Nucleation soluble',    '1e-4 #/kg', 'tttt')
      call ncinfo(ncname(10,:),'aisn','Free aerosol number concentration Aitken soluble',        '1e-4 #/kg', 'tttt')
      call ncinfo(ncname(11,:),'acsn','Free aerosol number concentration Accumulation soluble',  '1e-4 #/kg', 'tttt')
      call ncinfo(ncname(12,:),'cosn','Free aerosol number concentration Coarse soluble',        '1e-4 #/kg', 'tttt')
      call ncinfo(ncname(13,:),'aiin','Free aerosol number concentration Aitken insoluble',      '1e-4 #/kg', 'tttt')
      call ncinfo(ncname(14,:),'acin','Free aerosol number concentration Accumulation insoluble','1e-4 #/kg', 'tttt')
      call ncinfo(ncname(15,:),'coin','Free aerosol number concentration Coarse insoluble',      '1e-4 #/kg', 'tttt')
      call ncinfo(ncname(16,:),'nc','Cloud droplet number concentration',                        '1e-4 #/kg', 'tttt')
      call ncinfo(ncname(17,:),'nr','Rain droplet number concentration',                         '1e-4 #/kg', 'tttt')

      call ncinfo(ncname(18,:),'so4nus','Free aerosol sulphate mass concentration NUS', '1e-4 ng/kg','tttt')
      call ncinfo(ncname(19,:),'so4ais','Free aerosol sulphate mass concentration AIS', '1e-2 ng/kg','tttt')
      call ncinfo(ncname(20,:),'so4acs','Free aerosol sulphate mass concentration ACS', 'ng/kg','tttt')
      call ncinfo(ncname(21,:),'so4cos','Free aerosol sulphate mass concentration COS', 'ng/kg','tttt')
      call ncinfo(ncname(22,:),'so4cld','In-cloud aerosol sulphate mass concentration', '1e-2 ng/kg','tttt')
      call ncinfo(ncname(23,:),'so4rai','In-rain aerosol sulphate mass concentration',  '1e-2 ng/kg','tttt')

      call ncinfo(ncname(24,:),'bcais', 'Free aerosol black carbon mass concentration AIS', '1e-2 ng/kg','tttt')
      call ncinfo(ncname(25,:),'bcacs', 'Free aerosol black carbon mass concentration ACS', '1e-2 ng/kg','tttt')
      call ncinfo(ncname(26,:),'bccos', 'Free aerosol black carbon mass concentration COS', '1e-2 ng/kg','tttt')
      call ncinfo(ncname(27,:),'bcaii', 'Free aerosol black carbon mass concentration AII', '1e-2 ng/kg','tttt')
      call ncinfo(ncname(28,:),'bccld', 'In-cloud aerosol black carbon mass concentration', '1e-2 ng/kg','tttt')
      call ncinfo(ncname(29,:),'bcrai', 'In-rain aerosol black carbon mass concentration',  '1e-2 ng/kg','tttt')

      call ncinfo(ncname(30,:),'pomais','Free aerosol organic matter mass concentration AIS', '1e-2 ng/kg','tttt')
      call ncinfo(ncname(31,:),'pomacs','Free aerosol organic matter mass concentration ACS', '1e-2 ng/kg','tttt')
      call ncinfo(ncname(32,:),'pomcos','Free aerosol organic matter mass concentration COS', '1e-2 ng/kg','tttt')
      call ncinfo(ncname(33,:),'pomaii','Free aerosol organic matter mass concentration AII', '1e-2 ng/kg','tttt')
      call ncinfo(ncname(34,:),'pomcld','In-cloud aerosol organic matter mass concentration', '1e-2 ng/kg','tttt')
      call ncinfo(ncname(35,:),'pomrai','In-rain aerosol organic matter mass concentration',  '1e-2 ng/kg','tttt')
 
      call ncinfo(ncname(36,:),'ssacs', 'Free aerosol sea salt mass concentration ACS', 'ng/kg','tttt')
      call ncinfo(ncname(37,:),'sscos', 'Free aerosol sea salt mass concentration COS', 'ng/kg','tttt')
      call ncinfo(ncname(38,:),'sscld', 'In-cloud aerosol sea salt mass concentration', 'ng/kg','tttt')
      call ncinfo(ncname(39,:),'ssrai', 'In-rain aerosol sea salt mass concentration',  'ng/kg','tttt')

      call ncinfo(ncname(40,:),'duacs', 'Free aerosol mineral dust mass concentration ACS', '1e-2 ng/kg','tttt')
      call ncinfo(ncname(41,:),'ducos', 'Free aerosol mineral dust mass concentration COS', '1e-2 ng/kg','tttt')
      call ncinfo(ncname(42,:),'duaci', 'Free aerosol mineral dust mass concentration ACI', '1e-2 ng/kg','tttt')
      call ncinfo(ncname(43,:),'ducoi', 'Free aerosol mineral dust mass concentration COI', '1e-2 ng/kg','tttt')
      call ncinfo(ncname(44,:),'ducld', 'In-cloud aerosol mineral dust mass concentration', '1e-2 ng/kg','tttt')
      call ncinfo(ncname(45,:),'durai', 'In-rain aerosol mineral dust mass concentration',  '1e-2 ng/kg','tttt')
 
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
    use modfields,    only : um,vm,wm,thlm,qtm,ql0,svm,thv0h,thvh
    use modsurfdata,  only : thls,qts,thvs
    use modglobal,    only : imax,i1,ih,jmax,j1,jh,k1,rk3step,&
                             timee,dt_lim,cexpnr,ifoutput,rtimee
    use modmpi,       only : myid,cmyidx, cmyidy
    use modstat_nc,   only : lnetcdf, writestat_nc
    use modmicrodata, only : iqr, imicro, imicro_none, &
                             iaer_offset,              &   
                             inus_n,iais_n,iacs_n,icos_n,iaii_n,iaci_n,icoi_n,inc,inr, &   
                             iso4nus,iso4ais,iso4acs,iso4cos,iso4cld,iso4rai, &
                             ibcais, ibcacs, ibccos, ibcaii, ibccld, ibcrai,  &
                             ipomais,ipomacs,ipomcos,ipomaii,ipomcld,ipomrai, &
                             issacs, isscos,                 isscld, issrai,  &
                             iduacs, iducos, iduaci, iducoi, iducld, idurai   
         
    implicit none

    integer(KIND=selected_int_kind(4)), allocatable :: field(:,:,:),vars(:,:,:,:)
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

    field = NINT(1.0E3*um,2)
    if (lnetcdf) vars(:,:,:,1) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbuu.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbuu.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    field = NINT(1.0E3*vm,2)
    if (lnetcdf) vars(:,:,:,2) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbvv.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbvv.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    field = NINT(1.0E3*wm,2)
    if (lnetcdf) vars(:,:,:,3) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbww.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbww.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    field = NINT(1.0E5*qtm,2)
    if (lnetcdf) vars(:,:,:,4) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqt.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqt.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    field = NINT(1.0E6*ql0,2)
    if (lnetcdf) vars(:,:,:,5) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbql.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbql.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    field = NINT(1.0E2*(thlm-300),2)
    if (lnetcdf) vars(:,:,:,6) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbthl.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbthl.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    end if

    ! ===================================================================================================================    

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1.0E6*svm(i,j,k,iqr),2)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif

    if (lnetcdf) vars(:,:,:,7) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    field=0.
    do i=2-ih,i1+ih
    do j=2-jh,j1+jh
    do k=2,k1
      field(i,j,k) = NINT(1.0E2*(thv0h(i,j,k)-thvh(k)),2)
    enddo
    enddo
    enddo

    if (lnetcdf) vars(:,:,:,8) = field(2:i1,2:j1,klow:khigh)

    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbthv.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbthv.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1E-4*svm(i,j,k,inus_n+iaer_offset),4)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif

    if (lnetcdf) vars(:,:,:,9) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1E-4*svm(i,j,k,iais_n+iaer_offset),4)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif

    if (lnetcdf) vars(:,:,:,10) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1E-4*svm(i,j,k,iacs_n+iaer_offset),4)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif

    if (lnetcdf) vars(:,:,:,11) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1E-4*svm(i,j,k,icos_n+iaer_offset),4)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif

    if (lnetcdf) vars(:,:,:,12) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1E-4*svm(i,j,k,iaii_n+iaer_offset),4)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif

    if (lnetcdf) vars(:,:,:,13) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1E-4*svm(i,j,k,iaci_n+iaer_offset),4)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif

    if (lnetcdf) vars(:,:,:,14) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1E-4*svm(i,j,k,icoi_n+iaer_offset),4)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif

    if (lnetcdf) vars(:,:,:,15) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! = CLOUD DROPLETS ===================================================================================================    

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1e-4*svm(i,j,k,inc+iaer_offset),4)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif

    if (lnetcdf) vars(:,:,:,16) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1e-4*svm(i,j,k,inr+iaer_offset),4)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif

    if (lnetcdf) vars(:,:,:,17) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1.0E7*svm(i,j,k,iso4nus+iaer_offset),4)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif

    if (lnetcdf) vars(:,:,:,18) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1.0E5*svm(i,j,k,iso4ais+iaer_offset),4)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif

    if (lnetcdf) vars(:,:,:,19) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1.0E3*svm(i,j,k,iso4acs+iaer_offset),4)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif

    if (lnetcdf) vars(:,:,:,20) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1.0E3*svm(i,j,k,iso4cos+iaer_offset),4)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif

    if (lnetcdf) vars(:,:,:,21) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1.0E5*svm(i,j,k,iso4cld+iaer_offset),4)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif

    if (lnetcdf) vars(:,:,:,22) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1.0E5*svm(i,j,k,iso4rai+iaer_offset),4)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif

    if (lnetcdf) vars(:,:,:,23) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1.0E5*svm(i,j,k,ibcais+iaer_offset),4)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif

    if (lnetcdf) vars(:,:,:,24) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1.0E5*svm(i,j,k,ibcacs+iaer_offset),4)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif

    if (lnetcdf) vars(:,:,:,25) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1.0E5*svm(i,j,k,ibccos+iaer_offset),4)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif

    if (lnetcdf) vars(:,:,:,26) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1.0E5*svm(i,j,k,ibcaii+iaer_offset),4)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif

    if (lnetcdf) vars(:,:,:,27) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1.0E5*svm(i,j,k,ibccld+iaer_offset),4)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif

    if (lnetcdf) vars(:,:,:,28) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1.0E5*svm(i,j,k,ibcrai+iaer_offset),4)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif

    if (lnetcdf) vars(:,:,:,29) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1.0E5*svm(i,j,k,ipomais+iaer_offset),4)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif

    if (lnetcdf) vars(:,:,:,30) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1.0E5*svm(i,j,k,ipomacs+iaer_offset),4)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif

    if (lnetcdf) vars(:,:,:,31) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1.0E5*svm(i,j,k,ipomcos+iaer_offset),4)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif

    if (lnetcdf) vars(:,:,:,32) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1.0E5*svm(i,j,k,ipomaii+iaer_offset),4)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif

    if (lnetcdf) vars(:,:,:,33) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1.0E5*svm(i,j,k,ipomcld+iaer_offset),4)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif

    if (lnetcdf) vars(:,:,:,34) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1.0E5*svm(i,j,k,ipomrai+iaer_offset),4)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif

    if (lnetcdf) vars(:,:,:,35) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1.0E3*svm(i,j,k,issacs+iaer_offset),4)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif

    if (lnetcdf) vars(:,:,:,36) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1.0E3*svm(i,j,k,isscos+iaer_offset),4)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif

    if (lnetcdf) vars(:,:,:,37) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1.0E3*svm(i,j,k,isscld+iaer_offset),4)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif

    if (lnetcdf) vars(:,:,:,38) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1.0E3*svm(i,j,k,issrai+iaer_offset),4)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif

    if (lnetcdf) vars(:,:,:,39) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1.0E5*svm(i,j,k,iduacs+iaer_offset),4)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif

    if (lnetcdf) vars(:,:,:,40) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1.0E5*svm(i,j,k,iducos+iaer_offset),4)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif

    if (lnetcdf) vars(:,:,:,41) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1.0E5*svm(i,j,k,iduaci+iaer_offset),4)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif

    if (lnetcdf) vars(:,:,:,42) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1.0E5*svm(i,j,k,iducoi+iaer_offset),4)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif

    if (lnetcdf) vars(:,:,:,43) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1.0E5*svm(i,j,k,iducld+iaer_offset),4)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif

    if (lnetcdf) vars(:,:,:,44) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    if(imicro/=imicro_none) then
      do i=2-ih,i1+ih
      do j=2-jh,j1+jh
      do k=1,k1
        field(i,j,k) = NINT(1.0E5*svm(i,j,k,idurai+iaer_offset),4)
      enddo
      enddo
      enddo
    else
      field = 0.
    endif

    if (lnetcdf) vars(:,:,:,45) = field(2:i1,2:j1,klow:khigh)
    if (lbinary) then
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! ===================================================================================================================    

    if(lnetcdf) then
      call writestat_nc(ncid,1,tncname,(/rtimee/),nrec,.true.)
      call writestat_nc(ncid,nvar,ncname,vars,nrec,imax,jmax,khigh-klow+1)
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
