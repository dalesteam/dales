!> \file modcrosssection.f90
!!   Dumps an instantenous crosssection of the field

!>
!! Dumps an instantenous crosssection of the field.
!>
!! Crosssections in the yz-plane and in the xy-plane            |
    !        of u,v,w,thl,thv,qt,ql. Written to movv_*.expnr and movh_*.expnr
!! If netcdf is true, this module leads the cross.myid.expnr.nc output

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
module modcrosssection


  use modglobal, only : longint,kmax

implicit none
private
PUBLIC :: initcrosssection, crosssection,exitcrosssection
save
!NetCDF variables
	
  !____________________
  ! 	START 	Ruben Schulte, 23-04-2021
  ! Removing cross-section outputs and keep w, thl, qt, e120 and scalars
  integer :: nvar = 4
  ! 			END
  !____________________
  integer :: ncid1 = 0
  integer,allocatable :: ncid2(:)
  integer :: ncid3 = 1
  integer :: nrec1 = 0
  integer,allocatable :: nrec2(:)
  integer :: nrec3 = 0
  integer :: crossheight(100)
  integer :: nxy = 0
  integer :: cross
  character(4) :: cheight
  character(80) :: fname1 = 'crossxz.xxxxyxxx.xxx.nc'
  character(80) :: fname2 = 'crossxy.xxxx.xxxxyxxx.xxx.nc'
  character(80) :: fname3 = 'crossyz.xxxxyxxx.xxx.nc'
  ! 30-07-2020 Ruben Schulte start: changed "dimension(nvar,4) :: ncname1" into 
  !		"dimension(:,:), allocatable :: ncname"  
  character(80),dimension(:,:), allocatable :: ncname1
  character(80),dimension(1,4) :: tncname1
  character(80),dimension(:,:), allocatable :: ncname2
  character(80),dimension(1,4) :: tncname2
  character(80),dimension(:,:), allocatable :: ncname3
  character(80),dimension(1,4) :: tncname3
  ! 30-07-2020 Ruben Schulte end: changed "dimension(nvar,4) :: ncname1" into 
  !		"dimension(:,:), allocatable :: ncname"  

  real    :: dtav
  integer(kind=longint) :: idtav,tnext
  logical :: lcross = .false. !< switch for doing the crosssection (on/off)
  logical :: lbinary = .false. !< switch for doing the crosssection (on/off)
  integer :: crossplane = 2 !< Location of the xz crosssection
  integer :: crossortho = 2 !< Location of the yz crosssection

contains
!> Initializing Crosssection. Read out the namelist, initializing the variables
  subroutine initcrosssection
    use modmpi,   only :myid,my_real,mpierr,comm3d,mpi_logical,mpi_integer,cmyid,myidx,myidy
	! 30-07-2020 Ruben Schulte start: added nsv
    use modglobal,only :imax,jmax,ifnamopt,fname_options,dtmax,dtav_glob,ladaptive,j1,kmax,i1,dt_lim,cexpnr,tres,btime, nsv
	! 30-07-2020 Ruben Schulte end: added nsv
    use modstat_nc,only : lnetcdf,open_nc, define_nc,ncinfo,writestat_dims_nc
   implicit none

    integer :: ierr,k
	
	! 30-07-2020 Ruben Schulte start: Added as part of the do loops adding scalar data to netcdf
    character(3) :: csvname		
	integer n
    ! 30-07-2020 Ruben Schulte end: do loop adding scalar data to netcdf	

    namelist/NAMCROSSSECTION/ &
    lcross, lbinary, dtav, crossheight, crossplane, crossortho

    allocate(ncid2(kmax),nrec2(kmax))
    crossheight(1)=2
    crossheight(2:100)=-999
    ncid2(1)=2
    ncid2(2:kmax)=0
    nrec2(1:kmax)=0

    dtav = dtav_glob
    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMCROSSSECTION,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions NAMCROSSSECTION'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMCROSSSECTION'
      endif
      write(6 ,NAMCROSSSECTION)
      close(ifnamopt)
    end if

    call MPI_BCAST(dtav       ,1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(lcross     ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lbinary    ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(crossheight(1:100),100,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(crossplane ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(crossortho ,1,MPI_INTEGER,0,comm3d,mpierr)

    nxy=0
    k=1
    do while (crossheight(k) > 0)
    nxy=nxy+1
    ncid2(k)=k+1
    nrec2(k)=0
    k=k+1
    end do

    idtav = dtav/tres
    tnext   = idtav+btime
    if(.not.(lcross)) return
    dt_lim = min(dt_lim,tnext)

    if(any((crossheight(1:100).gt.kmax)) .or. crossplane>j1 .or. crossortho> i1 ) then
      stop 'CROSSSECTION: crosssection out of range'
    end if
    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'CROSSSECTION: dtav should be a integer multiple of dtmax'
    end if
    if (lnetcdf) then
	! 30-07-2020 Ruben Schulte start: Added in order to add scalar data to netcdf
	nvar = nvar + 1*nsv
	! 30-07-2020 Ruben Schulte end: Added in order to add scalar data to netcdf
    if (myidy==0) then
      fname1(9:16) = cmyid
      fname1(18:20) = cexpnr
	! 30-07-2020 Ruben Schulte start: Added in order to add scalar data to netcdf
      allocate(ncname1(nvar,4))
	! 30-07-2020 Ruben Schulte end: Added in order to add scalar data to netcdf
	
  !____________________
  ! 	START 	Ruben Schulte, 23-04-2021
  ! Removing cross-section outputs and keep w, thl, qt, e120 and scalars
      call ncinfo(tncname1(1,:),'time','Time','s','time')
      call ncinfo(ncname1( 1,:),'wxz', 'xz crosssection of the Vertical velocity','m/s','t0mt')
      call ncinfo(ncname1( 2,:),'thlxz','xz crosssection of the Liquid water potential temperature','K','t0tt')
      call ncinfo(ncname1( 3,:),'qtxz','xz crosssection of the Total water specific humidity','kg/kg','t0tt')
      call ncinfo(ncname1( 4,:),'e120xz','xz crosssection of sqrt(turbulent kinetic energy)','m^2/s^2','t0tt')
	  do n=1,nsv			
          write (csvname(1:3),'(i3.3)') n
		  call ncinfo(ncname1(4+n,:),'sv'//csvname//'xz','xz crosssection of scalar '//csvname//' specific mixing ratio','(kg/kg)','t0tt')
	  end do
  ! 			END
  !____________________
  
	  call open_nc(fname1,  ncid1,nrec1,n1=imax,n3=kmax)
      if (nrec1 == 0) then
        call define_nc( ncid1, 1, tncname1)
        call writestat_dims_nc(ncid1)
      end if
      call define_nc( ncid1, NVar, ncname1)
    end if
    do cross=1,nxy
      write(cheight,'(i4.4)') crossheight(cross)
      fname2(9:12) = cheight
      fname2(14:21) = cmyid
      fname2(23:25) = cexpnr
	! 30-07-2020 Ruben Schulte start: Added in order to add scalar data to netcdf
      allocate(ncname2(nvar,4))
	! 30-07-2020 Ruben Schulte end: Added in order to add scalar data to netcdf
	
  !____________________
  ! 	START 	Ruben Schulte, 23-04-2021
  ! Removing cross-section outputs and keep w, thl, qt, e120 and scalars
      call ncinfo(tncname2(1,:),'time','Time','s','time')
      call ncinfo(ncname2( 1,:),'wxy','xy crosssections of the Vertical velocity','m/s','tt0t')
      call ncinfo(ncname2( 2,:),'thlxy','xy crosssections of the Liquid water potential temperature','K','tt0t')
      call ncinfo(ncname2( 3,:),'qtxy','xy crosssections of the Total water specific humidity','kg/kg','tt0t')
      call ncinfo(ncname2( 4,:),'e120xy','xy crosssection of sqrt(turbulent kinetic energy)','m^2/s^2','tt0t')
	  do n=1,nsv			
          write (csvname(1:3),'(i3.3)') n
		  call ncinfo(ncname2(4+n,:),'sv'//csvname//'xy','xy crosssection of scalar '//csvname//' specific mixing ratio','(kg/kg)','tt0t')
	  end do
  ! 			END
  !____________________
  
      call open_nc(fname2,  ncid2(cross),nrec2(cross),n1=imax,n2=jmax)
      if (nrec2(cross)==0) then
        call define_nc( ncid2(cross), 1, tncname2)
        call writestat_dims_nc(ncid2(cross))
      end if
      call define_nc( ncid2(cross), NVar, ncname2)
   end do
   if (myidx==0) then
    fname3(9:16) = cmyid
    fname3(18:20) = cexpnr
	! 30-07-2020 Ruben Schulte start: Added in order to add scalar data to netcdf
      allocate(ncname3(nvar,4))
	! 30-07-2020 Ruben Schulte end: Added in order to add scalar data to netcdf
	
  !____________________
  ! 	START 	Ruben Schulte, 23-04-2021
  ! Removing cross-section outputs and keep w, thl, qt, e120 and scalars
    call ncinfo(tncname3(1,:),'time','Time','s','time')
    call ncinfo(ncname3( 1,:),'wyz','yz crosssection of the Vertical velocity','m/s','0tmt')
    call ncinfo(ncname3( 2,:),'thlyz','yz crosssection of the Liquid water potential temperature','K','0ttt')
    call ncinfo(ncname3( 3,:),'qtyz','yz crosssection of the Total water specific humidity','kg/kg','0ttt')
    call ncinfo(ncname3( 4,:),'e120yz','yz crosssection of sqrt(turbulent kinetic energy)','m^2/s^2','0ttt')
	  do n=1,nsv			
          write (csvname(1:3),'(i3.3)') n
		  call ncinfo(ncname3(4+n,:),'sv'//csvname//'yz','yz crosssection of scalar '//csvname//' specific mixing ratio','(kg/kg)','0ttt')
	  end do
  ! 			END
  !____________________
  
    call open_nc(fname3,  ncid3,nrec3,n2=jmax,n3=kmax)
    if (nrec3==0) then
      call define_nc( ncid3, 1, tncname3)
      call writestat_dims_nc(ncid3)
    end if
    call define_nc( ncid3, NVar, ncname3)
    end if
 end if


  end subroutine initcrosssection
!>Run crosssection. Mainly timekeeping
  subroutine crosssection
    use modglobal, only : rk3step,timee,dt_lim
    use modstat_nc, only : writestat_nc
    implicit none


    if (.not. lcross) return
    if (rk3step/=3) return
    if(timee<tnext) then
      dt_lim = min(dt_lim,tnext-timee)
      return
    end if
    tnext = tnext+idtav
    dt_lim = minval((/dt_lim,tnext-timee/))

    call wrtvert
    call wrthorz
    call wrtorth

  end subroutine crosssection


!> Do the xz crosssections and dump them to file
  subroutine wrtvert
  use modglobal, only : imax,i1,kmax,nsv,rlv,cp,rv,rd,cu,cv,cexpnr,ifoutput,rtimee
  use modfields, only : um,vm,wm,thlm,qtm,svm,thl0,qt0,ql0,e120,exnf,thvf,cloudnr
  use modmpi,    only : myidy
  use modstat_nc, only : lnetcdf, writestat_nc
  implicit none

  integer i,k,n
  character(20) :: name

  real, allocatable :: thv0(:,:),vars(:,:,:),buoy(:,:)

  if( myidy /= 0 ) return 

  allocate(thv0(2:i1,1:kmax),buoy(2:i1,1:kmax))


    do  i=2,i1
    do  k=1,kmax
      thv0(i,k) = (thl0(i,crossplane,k)+rlv*ql0(i,crossplane,k)/(cp*exnf(k))) &
                    *(1+(rv/rd-1)*qt0(i,crossplane,k)-rv/rd*ql0(i,crossplane,k))
      buoy(i,k) = thv0(i,k)-thvf(k)
    enddo
    enddo

    if(lbinary) then
      open(ifoutput,file='movv_u.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((um(i,crossplane,k)+cu,i=2,i1),k=1,kmax)
      close(ifoutput)

      open(ifoutput,file='movv_v.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((vm(i,crossplane,k)+cv,i=2,i1),k=1,kmax)
      close(ifoutput)

      open(ifoutput,file='movv_w.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((wm(i,crossplane,k),i=2,i1),k=1,kmax)
      close(ifoutput)

      open(ifoutput,file='movv_thl.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((thlm(i,crossplane,k),i=2,i1),k=1,kmax)
      close(ifoutput)

      open(ifoutput,file='movv_thv.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((thv0(i,k),i=2,i1),k=1,kmax)
      close(ifoutput)

      open(ifoutput,file='movv_buoy.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((buoy(i,k),i=2,i1),k=1,kmax)
      close(ifoutput)

      open(ifoutput,file='movv_qt.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((1.e3*qtm(i,crossplane,k),i=2,i1),k=1,kmax)
      close(ifoutput)

      open(ifoutput,file='movv_ql.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((1.e3*ql0(i,crossplane,k),i=2,i1),k=1,kmax)
      close(ifoutput)

      do n = 1,nsv
        name = 'movh_tnn.'//cexpnr
        write(name(7:8),'(i2.2)') n
        open(ifoutput,file=name,position='append',action='write')
        write(ifoutput,'(es12.5)') ((svm(i,crossplane,k,n),i=2,i1),k=1,kmax)
        close(ifoutput)
      end do
    end if

    if (lnetcdf) then
      allocate(vars(1:imax,1:kmax,nvar))
  !____________________
  ! 	START 	Ruben Schulte, 23-04-2021
  ! Removing cross-section outputs and keep w, thl, qt, e120 and scalars
      vars(:,:,1) = wm(2:i1,crossplane,1:kmax)
      vars(:,:,2) = thlm(2:i1,crossplane,1:kmax)
      vars(:,:,3) = qtm(2:i1,crossplane,1:kmax)
      vars(:,:,4) = e120(2:i1,crossplane,1:kmax)
      ! 30-07-2020 Ruben Schulte start: do loop adding scalar data to netcdf
      do n=1,nsv
          vars(:,:,4+n) = svm(2:i1,crossplane,1:kmax,n)
      end do
	  ! 30-07-2020 Ruben Schulte end: do loop adding scalar data to netcdf
  ! 			END
  !____________________
  
      call writestat_nc(ncid1,1,tncname1,(/rtimee/),nrec1,.true.)
      call writestat_nc(ncid1,nvar,ncname1(1:nvar,:),vars,nrec1,imax,kmax)
      deallocate(vars)
    end if
    deallocate(thv0,buoy)

  end subroutine wrtvert

!> Do the xy crosssections and dump them to file
  subroutine wrthorz
    use modglobal, only : imax,jmax,i1,j1,nsv,rlv,cp,rv,rd,cu,cv,cexpnr,ifoutput,rtimee
    use modfields, only : um,vm,wm,thlm,qtm,svm,thl0,qt0,ql0,e120,exnf,thvf,cloudnr
    use modmpi,    only : cmyid
    use modstat_nc, only : lnetcdf, writestat_nc
    use modmicrodata, only : iqr,inr
    implicit none


    ! LOCAL
    integer i,j,n
    character(40) :: name
    real, allocatable :: thv0(:,:,:),vars(:,:,:),buoy(:,:,:)

    allocate(thv0(2:i1,2:j1,nxy),buoy(2:i1,2:j1,nxy))

    do  cross=1,nxy
    do  j=2,j1
    do  i=2,i1
      thv0(i,j,cross) =&
       (thl0(i,j,crossheight(cross))+&
       rlv*ql0(i,j,crossheight(cross))/&
       (cp*exnf(crossheight(cross)))) &
                    *(1+(rv/rd-1)*qt0(i,j,crossheight(cross))&
                    -rv/rd*ql0(i,j,crossheight(cross)))
      buoy(i,j,cross) =thv0(i,j,cross)-thvf(crossheight(cross))
    enddo
    enddo
    enddo

    if(lbinary) then
      do  cross=1,nxy
      write(cheight,'(i4.4)') crossheight(cross)
      open(ifoutput,file='movh_u.'//cheight//'.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((um(i,j,crossheight(cross))+cu,i=2,i1),j=2,j1)
      close(ifoutput)

      open(ifoutput,file='movh_v.'//cheight//'.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((vm(i,j,crossheight(cross))+cv,i=2,i1),j=2,j1)
      close(ifoutput)

      open(ifoutput,file='movh_w.'//cheight//'.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((wm(i,j,crossheight(cross)),i=2,i1),j=2,j1)
      close(ifoutput)

      open(ifoutput,file='movh_thl.'//cheight//'.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((thlm(i,j,crossheight(cross)),i=2,i1),j=2,j1)
      close(ifoutput)

      open(ifoutput,file='movh_thv.'//cheight//'.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((thv0(i,j,cross),i=2,i1),j=2,j1)
      close(ifoutput)

      open(ifoutput,file='movh_buoy.'//cheight//'.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((buoy(i,j,cross),i=2,i1),j=2,j1)
      close(ifoutput)

      open(ifoutput,file='movh_qt.'//cheight//'.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((1.e3*qtm(i,j,crossheight(cross)),i=2,i1),j=2,j1)
      close(ifoutput)

      open(ifoutput,file='movh_ql.'//cheight//'.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((1.e3*ql0(i,j,crossheight(cross)),i=2,i1),j=2,j1)
      close(ifoutput)

      do n = 1,nsv
        name = 'movh_snn.'//trim(cheight)//'.'//cmyid//'.'//cexpnr
        write(name(7:8),'(i2.2)') n
        open(ifoutput,file=name,position='append',action='write')
        write(ifoutput,'(es12.5)') ((svm(i,j,crossheight(cross),n),i=2,i1),j=2,j1)
        close(ifoutput)
      end do
      end do
    endif

    if (lnetcdf) then
      do cross=1,nxy
      allocate(vars(1:imax,1:jmax,nvar))
      vars=0.
  !____________________
  ! 	START 	Ruben Schulte, 23-04-2021
  ! Removing cross-section outputs and keep w, thl, qt, e120 and scalars
      vars(:,:,1) = wm(2:i1,2:j1,crossheight(cross))
      vars(:,:,2) = thlm(2:i1,2:j1,crossheight(cross))
      vars(:,:,3) = qtm(2:i1,2:j1,crossheight(cross))
      vars(:,:,4) = e120(2:i1,2:j1,crossheight(cross))
      ! 30-07-2020 Ruben Schulte start: do loop adding scalar data to netcdf
      do n=1,nsv
          vars(:,:,4+n) = svm(2:i1,2:j1,crossheight(cross),n)
      end do
	  ! 30-07-2020 Ruben Schulte end: do loop adding scalar data to netcdf
  ! 			END
  !____________________
      call writestat_nc(ncid2(cross),1,tncname2,(/rtimee/),nrec2(cross),.true.)
      call writestat_nc(ncid2(cross),nvar,ncname2(1:nvar,:),vars,nrec2(cross),imax,jmax)
      deallocate(vars)
      end do
    end if

    deallocate(thv0,buoy)

  end subroutine wrthorz

  ! yz cross section
  subroutine wrtorth
    use modglobal, only : jmax,kmax,j1,nsv,rlv,cp,rv,rd,cu,cv,cexpnr,ifoutput,rtimee
    use modfields, only : um,vm,wm,thlm,qtm,svm,thl0,qt0,ql0,e120,exnf,thvf,cloudnr
    use modmpi,    only : cmyid, myidx
    use modstat_nc, only : lnetcdf, writestat_nc
    implicit none


    ! LOCAL
    integer j,k,n
    character(21) :: name

    real, allocatable :: thv0(:,:),vars(:,:,:),buoy(:,:)

    if( myidx /= 0 ) return 


    allocate(thv0(1:j1,1:kmax),buoy(1:j1,1:kmax))

    do  j=1,j1
    do  k=1,kmax
      thv0(j,k) =&
       (thl0(crossortho,j,k)+&
       rlv*ql0(crossortho,j,k)/&
       (cp*exnf(k))) &
                    *(1+(rv/rd-1)*qt0(crossortho,j,k)&
                    -rv/rd*ql0(crossortho,j,k))
       buoy(j,k) =thv0(j,k)-thvf(k)
    enddo
    enddo

    if(lbinary) then
      open(ifoutput,file='movo_u.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((um(crossortho,j,k)+cu,j=2,j1),k=1,kmax)
      close(ifoutput)

      open(ifoutput,file='movo_v.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((vm(crossortho,j,k)+cv,j=2,j1),k=1,kmax)
      close(ifoutput)

      open(ifoutput,file='movo_w.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((wm(crossortho,j,k),j=2,j1),k=1,kmax)
      close(ifoutput)

      open(ifoutput,file='movo_thl.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((thlm(crossortho,j,k),j=2,j1),k=1,kmax)
      close(ifoutput)

      open(ifoutput,file='movo_thv.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((thv0(j,k),j=2,j1),k=1,kmax)
      close(ifoutput)

      open(ifoutput,file='movo_buoy.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((buoy(j,k),j=2,j1),k=1,kmax)
      close(ifoutput)

      open(ifoutput,file='movo_qt.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((1.e3*qtm(crossortho,j,k),j=2,j1),k=1,kmax)
      close(ifoutput)

      open(ifoutput,file='movo_ql.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((1.e3*ql0(crossortho,j,k),j=2,j1),k=1,kmax)
      close(ifoutput)

      do n = 1,nsv
        name = 'movh_tnn.'//cmyid//'.'//cexpnr
        write(name(7:8),'(i2.2)') n
        open(ifoutput,file=name,position='append',action='write')
        write(ifoutput,'(es12.5)') ((svm(crossortho,j,k,n),j=2,j1),k=1,kmax)
        close(ifoutput)
      end do
    end if

    if (lnetcdf) then
      allocate(vars(1:jmax,1:kmax,nvar))
  !____________________
  ! 	START 	Ruben Schulte, 23-04-2021
  ! Removing cross-section outputs and keep w, thl, qt, e120 and scalars
      vars(:,:,1) = wm(crossortho,2:j1,1:kmax)
      vars(:,:,2) = thlm(crossortho,2:j1,1:kmax)
      vars(:,:,3) = qtm(crossortho,2:j1,1:kmax)
      vars(:,:,4) = e120(crossortho,2:j1,1:kmax)
      ! 30-07-2020 Ruben Schulte start: do loop adding scalar data to netcdf
      do n=1,nsv
          vars(:,:,4+n) = svm(crossortho,2:j1,1:kmax,n)
      end do
	  ! 30-07-2020 Ruben Schulte end: do loop adding scalar data to netcdf
  ! 			END
  !____________________
      call writestat_nc(ncid3,1,tncname3,(/rtimee/),nrec3,.true.)
      call writestat_nc(ncid3,nvar,ncname3(1:nvar,:),vars,nrec3,jmax,kmax)
      deallocate(vars)
    end if

    deallocate(thv0,buoy)

  end subroutine wrtorth

!> Clean up when leaving the run
  subroutine exitcrosssection
    use modstat_nc, only : exitstat_nc,lnetcdf
    use modmpi, only : myid
    implicit none

    if(lcross .and. lnetcdf) then
    if (myid==0) then
    call exitstat_nc(ncid1)
    end if
    do cross=1,nxy
    call exitstat_nc(ncid2(cross))
    end do
    call exitstat_nc(ncid3)
    end if

  end subroutine exitcrosssection

end module modcrosssection
