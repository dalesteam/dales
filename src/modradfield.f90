!> \file modradfield.f90
!!  Dumps 2D fields of several variables
!>
!!  \author Stephan de Roode,TU Delft
!!  \author Fredrik Jansson,TU Delft
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
module modradfield

  use modglobal, only : longint
  use modprecision, only : field_r

implicit none
private
PUBLIC :: initradfield, radfield, exitradfield
save
!NetCDF variables
  integer,parameter :: nvar = 22
  integer :: ncid2,nrec2 = 0
  integer :: nsamples
  character(80) :: fname
  character(80),dimension(nvar,4) :: ncname
  character(80),dimension(1,4) :: tncname

  real    :: dtav,timeav
  integer(kind=longint) :: idtav,itimeav,tnext,tnextwrite
  logical :: lradfield= .false. !< switch to enable the fielddump (on/off)

  real, allocatable :: field_2D_mn (:,:,:)

contains
!> Initializing fielddump. Read out the namelist, initializing the variables
  subroutine initradfield
    use modmpi,   only :myid,comm3d,myidx,myidy,D_MPI_BCAST
    use modglobal,only :imax,jmax,i1,ih,j1,jh,cexpnr,ifnamopt,fname_options,dtmax,dtav_glob,&
         timeav_glob,ladaptive,dt_lim,btime,tres,checknamelisterror,output_prefix
    use modstat_nc,only : open_nc, define_nc,ncinfo,writestat_dims_nc
    implicit none
    integer :: ierr

    namelist/NAMRADFIELD/ &
         dtav,timeav,lradfield

    dtav=dtav_glob
    timeav=timeav_glob

    if(myid==0)then
       open(ifnamopt,file=fname_options,status='old',iostat=ierr)
       read(ifnamopt,NAMRADFIELD,iostat=ierr)
       call checknamelisterror(ierr, ifnamopt, 'NAMRADFIELD')
       write(6, NAMRADFIELD)
       close(ifnamopt)
    end if
    call D_MPI_BCAST(dtav       ,1,0,comm3d,ierr)
    call D_MPI_BCAST(timeav     ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lradfield  ,1,0,comm3d,ierr)
    idtav = dtav/tres
    itimeav = timeav/tres

    tnext      = idtav   +btime ! sample time
    tnextwrite = itimeav +btime ! write time
    nsamples = itimeav/idtav

    if(.not.(lradfield)) return
    dt_lim = min(dt_lim,tnext)


    if (abs(timeav/dtav-nsamples)>1e-4) then
       stop 'radfield timeav must be a integer multiple of dtav'
    end if

    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
       stop 'radfield dtav should be a integer multiple of dtmax'
    end if

    allocate(field_2D_mn(2-ih:i1+ih,2-jh:j1+jh,nvar))
    field_2D_mn = 0

    write(fname,'(A,i3.3,A,i3.3,A)') 'radfield.', myidx, '.', myidy, '.xxx.nc'    !rce table 5 2D hourly averaged variables
    fname(18:20) = cexpnr
    call ncinfo(tncname(1,:),'time','Time','s','time')
    call ncinfo(ncname( 1,:),'hfls','surface upward latent heat flux','W/m2','tt0t')
    call ncinfo(ncname( 2,:),'hfss','surface upward sensible heat flux','W/m2','tt0t')
    call ncinfo(ncname( 3,:),'rlds','surface downwellling longwave flux','W/m2','tt0t')
    call ncinfo(ncname( 4,:),'rlus','surface upwelling longwave flux','W/m2','tt0t')
    call ncinfo(ncname( 5,:),'rsds','surface downwellling shortwave flux','W/m2','tt0t')
    call ncinfo(ncname( 6,:),'rsus','surface upwellling shortwave flux','W/m2','tt0t')
    call ncinfo(ncname( 7,:),'rsdt','TOM incoming shortwave flux','units','tt0t')
    call ncinfo(ncname( 8,:),'rsut','TOM outgoing shortwave flux','units','tt0t')
    call ncinfo(ncname( 9,:),'rlut','TOM outgoing longwave flux','units','tt0t')
    call ncinfo(ncname(10,:),'rsdscs','surface downwelling shortwave flux - clear sky','W/m2','tt0t')
    call ncinfo(ncname(11,:),'rsuscs','surface upwelling shortwave flux - clear sky','W/m2','tt0t')
    call ncinfo(ncname(12,:),'rldscs','surface downwelling longwave flux - clear sky','W/m2','tt0t')
    call ncinfo(ncname(13,:),'rluscs','surface upwelling longwave flux - clear sky','W/m2','tt0t')
    call ncinfo(ncname(14,:),'rsutcs','TOM outgoing shortwave flux - clear sky','W/m2','tt0t')
    call ncinfo(ncname(15,:),'rlutcs','TOM outgoing longwave flux - clear sky','W/m2','tt0t')
    call ncinfo(ncname(16,:),'rsds_dir','surface downwellling shortwave direct flux','W/m2','tt0t')
    call ncinfo(ncname(17,:),'rsds_dif','surface downwellling shortwave diffuse flux','W/m2','tt0t')

    call ncinfo(ncname(18,:),'prw','water vapor path','kg/m2','tt0t')
    call ncinfo(ncname(19,:),'clwvi','condensed water path','kg/m2','tt0t') ! both liquid and ice
    !call ncinfo(ncname(20,:),'clivi','ice water path','kg/m2','tt0t')
    !call ncinfo(ncname(21,:),'spwr','saturated water vapor path','kg/m2','tt0t')
    call ncinfo(ncname(20,:),'uabot','eastward wind at lowest model level','m/s','mt0t')
    call ncinfo(ncname(21,:),'vabot','northward wind at lowest model level','m/s','tm0t')
    call ncinfo(ncname(22,:),'tabot','air temperature at lowest model level','K','tt0t')

    call open_nc(trim(output_prefix)//fname,  ncid2,nrec2,n1=imax,n2=jmax,n3=1)
    if (nrec2==0) then
       call define_nc( ncid2, 1, tncname)
       call writestat_dims_nc(ncid2)
    end if
    call define_nc( ncid2, nvar, ncname)
  end subroutine initradfield


  subroutine radfield
    use modglobal, only : rk3step,timee,dt_lim
    implicit none

    if (.not. lradfield) return
    if (rk3step/=3) return

    if(timee<tnext .and. timee<tnextwrite) then
      dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))
      return
    end if
    if (timee>=tnext) then
      tnext = tnext+idtav
      call sample_radfield
    end if
    if (timee>=tnextwrite) then
      tnextwrite = tnextwrite+itimeav
      call writestat_radfield
    end if
    dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))
  end subroutine radfield


  subroutine sample_radfield
    use modfields, only : rhof,qt0,ql0,u0,v0,exnf,thl0,qt0
    use modsurfdata, only: qtflux,thlflux
    use modglobal, only: dzf,tup,tdn,i1,j1,kmax,rlv,cp
    use modraddata, only : lwd,lwu,swd,swu,lwdca,lwuca,swdca,swuca,swdir,swdif
    implicit none
    integer :: i,j,k
    real :: tmp !,ilratio,qsat

    ! rcemip , neglect density difference zf and zh (consistent with surface flux parameterization)
    field_2D_mn (2:i1,2:j1,1) = field_2D_mn (2:i1,2:j1,1) + rhof(1) * rlv * qtflux (2:i1,2:j1)
    field_2D_mn (2:i1,2:j1,2) = field_2D_mn (2:i1,2:j1,2) + rhof(1) * cp * thlflux(2:i1,2:j1)

    ! radiation fields are allocated and initialized to 0 even if radiation is off
    field_2D_mn (2:i1,2:j1,3) = field_2D_mn (2:i1,2:j1,3) + lwd(2:i1,2:j1,1)
    field_2D_mn (2:i1,2:j1,4) = field_2D_mn (2:i1,2:j1,4) + lwu(2:i1,2:j1,1)
    field_2D_mn (2:i1,2:j1,5) = field_2D_mn (2:i1,2:j1,5) + swd(2:i1,2:j1,1)
    field_2D_mn (2:i1,2:j1,6) = field_2D_mn (2:i1,2:j1,6) + swu(2:i1,2:j1,1)
    field_2D_mn (2:i1,2:j1,7) = field_2D_mn (2:i1,2:j1,7) + swd(2:i1,2:j1,kmax)
    field_2D_mn (2:i1,2:j1,8) = field_2D_mn (2:i1,2:j1,8) + swu(2:i1,2:j1,kmax)
    field_2D_mn (2:i1,2:j1,9) = field_2D_mn (2:i1,2:j1,9) + lwu(2:i1,2:j1,kmax)
    field_2D_mn (2:i1,2:j1,10) = field_2D_mn (2:i1,2:j1,10) + swdca(2:i1,2:j1,1)
    field_2D_mn (2:i1,2:j1,11) = field_2D_mn (2:i1,2:j1,11) + swuca(2:i1,2:j1,1)
    field_2D_mn (2:i1,2:j1,12) = field_2D_mn (2:i1,2:j1,12) + lwdca(2:i1,2:j1,1)
    field_2D_mn (2:i1,2:j1,13) = field_2D_mn (2:i1,2:j1,13) + lwuca(2:i1,2:j1,1)
    field_2D_mn (2:i1,2:j1,14) = field_2D_mn (2:i1,2:j1,14) + swuca(2:i1,2:j1,kmax)
    field_2D_mn (2:i1,2:j1,15) = field_2D_mn (2:i1,2:j1,15) + lwuca(2:i1,2:j1,kmax)
    field_2D_mn (2:i1,2:j1,16) = field_2D_mn (2:i1,2:j1,16) + swdir(2:i1,2:j1,kmax)
    field_2D_mn (2:i1,2:j1,17) = field_2D_mn (2:i1,2:j1,17) + swdif(2:i1,2:j1,kmax)

    do k=1,kmax
       do j=2,j1
          do i=2,i1
             field_2D_mn (i,j,18) = field_2D_mn (i,j,18) + rhof(k) * (qt0(i,j,k) - ql0(i,j,k)) * dzf(k)
             field_2D_mn (i,j,19) = field_2D_mn (i,j,19) + rhof(k) *  ql0(i,j,k) * dzf(k)

             ! ice water path and qsat-path disabled in Fugaku branch
             ! due to no ice and new thermodynamics where tmp, qsat are not saved as 3d fields

             !tmp = exnf(k)*thl0(i,j,k)  + (rlv/cp) * ql0(i,j,k)
             !ilratio=max(0.,min(1.,(tmp-tdn)/(tup-tdn)))! cloud water vs cloud ice partitioning
             !qsat=
             !field_2D_mn (i,j,20) = field_2D_mn (i,j,20) + rhof(k) *  ql0(i,j,k) * dzf(k) *(1-ilratio)
             !field_2D_mn (i,j,20) = field_2D_mn (i,j,20) + rhof(k) *  ql0(i,j,k) * dzf(k)
             !field_2D_mn (i,j,21) = field_2D_mn (i,j,21) + rhof(k) *  qsat * dzf(k)
          end do
       end do
    end do


    field_2D_mn (2:i1,2:j1,20) = field_2D_mn (2:i1,2:j1,20) + u0(2:i1,2:j1,1)
    field_2D_mn (2:i1,2:j1,21) = field_2D_mn (2:i1,2:j1,21) + v0(2:i1,2:j1,1)

    do j=2,j1
       do i=2,i1
          tmp = exnf(k)*thl0(i,j,1)  + (rlv/cp) * ql0(i,j,1)
          field_2D_mn (2:i1,2:j1,22) = field_2D_mn (2:i1,2:j1,22) + tmp
       end do
    end do

  end subroutine sample_radfield


  subroutine writestat_radfield
    use modglobal, only : imax,jmax,i1,j1,rtimee
    use modstat_nc, only : writestat_nc

    implicit none

    field_2D_mn (2:i1,2:j1,:) = field_2D_mn (2:i1,2:j1,:) / nsamples

    call writestat_nc(ncid2,1,tncname,(/rtimee/),nrec2,.true.)
    call writestat_nc(ncid2,nvar,ncname,field_2D_mn(2:i1,2:j1,:),nrec2,imax,jmax)

    field_2D_mn = 0.
  end subroutine writestat_radfield

!> Clean up when leaving the run
  subroutine exitradfield
    use modstat_nc, only : exitstat_nc,lnetcdf
    implicit none

    if(lradfield) call exitstat_nc(ncid2)
    if(lradfield) deallocate(field_2D_mn)
  end subroutine exitradfield

end module modradfield
