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
  integer,parameter :: nvar = 11
  integer :: ncid1(100)
  integer :: ncid2(100)
  integer :: ncid3(100)
  integer :: nrec1(100) = 0
  integer :: nrec2(100) = 0
  integer :: nrec3(100) = 0
  integer :: crossheight(100)
  integer :: nxy = 0, nxz = 0, nyz = 0
  integer :: cross
  character(4) :: cheight
  character(80) :: fname1 = 'crossxz.yyyy.x000.exp.nc'
  character(80) :: fname2 = 'crossxy.zzzz.x000y000.exp.nc'
  character(80) :: fname3 = 'crossyz.xxxx.y000.exp.nc'
  character(80),dimension(nvar,4) :: ncname1
  character(80),dimension(1,4) :: tncname1
  character(80),dimension(nvar,4) :: ncname2
  character(80),dimension(1,4) :: tncname2
  character(80),dimension(nvar,4) :: ncname3
  character(80),dimension(1,4) :: tncname3

  real    :: dtav
  integer(kind=longint) :: idtav,tnext
  logical :: lcross = .false. !< switch for doing the crosssection (on/off)
  logical :: lbinary = .false. !< switch for doing the crosssection (on/off)
  integer :: crossplane(100) !< Location of the xz crosssection
  integer :: crossortho(100) !< Location of the yz crosssection

  logical :: lxy = .true.   !< switch for doing xy crosssections
  logical :: lxz = .true.   !< switch for doing xz crosssections
  logical :: lyz = .true.   !< switch for doing yz crosssections


contains
!> Initializing Crosssection. Read out the namelist, initializing the variables
  subroutine initcrosssection
    use modmpi,   only :myid,mpierr,comm3d,mpi_logical,mpi_integer,cmyid,cmyidx,cmyidy,myidx,myidy &
                       , D_MPI_BCAST
    use modglobal,only :imax,jmax,itot,jtot,ifnamopt,fname_options,dtmax,dtav_glob,ladaptive,j1,kmax,i1,dt_lim,cexpnr,&
                        tres,btime,checknamelisterror,output_prefix
    use modstat_nc,only : lnetcdf,open_nc, define_nc,ncinfo,nctiminfo,writestat_dims_nc

   implicit none

    integer :: ierr,k

    namelist/NAMCROSSSECTION/ &
    lcross, lbinary, dtav, crossheight, crossplane, crossortho, lxy, lxz, lyz

    crossheight(1)=2
    crossheight(2:100)=-999
    crossplane(1)=2
    crossplane(2:100)=-999
    crossortho(1)=2
    crossortho(2:100)=-999

    dtav = dtav_glob
    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMCROSSSECTION,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMCROSSSECTION')
      write(6 ,NAMCROSSSECTION)
      close(ifnamopt)
    end if

    call D_MPI_BCAST(dtav       ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(lcross     ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(lbinary    ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(crossheight(1:100),100,0,comm3d,mpierr)
    call D_MPI_BCAST(crossplane(1:100),100,0,comm3d,mpierr)
    call D_MPI_BCAST(crossortho(1:100),100,0,comm3d,mpierr)
    call D_MPI_BCAST(lxy        ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(lxz        ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(lyz        ,1,0,comm3d,mpierr)

    if(any((crossheight(1:100).gt.kmax)) .or. any(crossplane > jtot+1) .or. any(crossortho > itot+1) ) then
      stop 'CROSSSECTION: crosssection out of range'
    end if
    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'CROSSSECTION: dtav should be a integer multiple of dtmax'
    end if

    k=1
    do while (crossheight(k) > 0)
       nxy=nxy+1
       k=k+1
    end do

    k=1
    do while (crossplane(k) > 0)
       crossplane(k) = crossplane(k) - myidy*jmax  ! convert to local grid index
       nxz=nxz+1
       k=k+1
    end do

    k=1
    do while (crossortho(k) > 0)
       crossortho(k) = crossortho(k) - myidx*imax  ! convert to local grid index
       nyz=nyz+1
       k=k+1
    end do

    ! cross sections with
    ! 2  <= crossplane, <= j1
    ! 2  <= crossortho <= i1
    ! belong to this processor

    idtav = dtav/tres
    tnext   = idtav+btime
    if(.not.(lcross)) return
    dt_lim = min(dt_lim,tnext)

    if (lnetcdf) then
       if (lxz) then
          do cross=1,nxz
             if (crossplane(cross) >= 2 .and. crossplane(cross) <= j1) then
                !      fname1(9:16) = cmyid
                !      fname1(18:20) = cexpnr
                write(cheight,'(i4.4)') crossplane(cross) + myidy*jmax
                fname1(9:12) = cheight
                fname1(15:17) = cmyidx
                fname1(19:21) = cexpnr
                call nctiminfo(tncname1(1,:))
                call ncinfo(ncname1( 1,:),'uxz', 'xz crosssection of the West-East velocity','m/s','m0tt')
                call ncinfo(ncname1( 2,:),'vxz', 'xz crosssection of the South-North velocity','m/s','t0tt')
                call ncinfo(ncname1( 3,:),'wxz', 'xz crosssection of the Vertical velocity','m/s','t0mt')
                call ncinfo(ncname1( 4,:),'thlxz','xz crosssection of the Liquid water potential temperature','K','t0tt')
                call ncinfo(ncname1( 5,:),'thvxz','xz crosssection of the Virtual potential temperature','K','t0tt')
                call ncinfo(ncname1( 6,:),'qtxz','xz crosssection of the Total water specific humidity','kg/kg','t0tt')
                call ncinfo(ncname1( 7,:),'qlxz','xz crosssection of the Liquid water specific humidity','kg/kg','t0tt')
                call ncinfo(ncname1( 8,:),'buoyxz','xz crosssection of the Buoyancy','K','t0tt')
                call ncinfo(ncname1( 9,:),'qrxz','xz crosssection of the Rain water specific humidity','kg/kg','t0tt')
                call ncinfo(ncname1( 10,:),'nrxz','xz crosssection of the Number concentration','-','t0tt')
                call ncinfo(ncname1( 11,:),'e120xz','xz crosssection of sqrt(turbulent kinetic energy)','m^2/s^2','t0tt')
                call open_nc(trim(output_prefix)//fname1,ncid1(cross),nrec1(cross),n1=imax,n3=kmax)
                if (nrec1(cross) == 0) then
                   call define_nc(ncid1(cross), 1, tncname1)
                   call writestat_dims_nc(ncid1(cross))
                end if
                call define_nc(ncid1(cross), NVar, ncname1)
             end if
          end do
    end if
    if (lxy) then
       do cross=1,nxy
          write(cheight,'(i4.4)') crossheight(cross)
          fname2(9:12) = cheight
          fname2(14:21) = cmyid
          fname2(23:25) = cexpnr
          call nctiminfo(tncname2(1,:))
          call ncinfo(ncname2( 1,:),'uxy','xy crosssections of the West-East velocity','m/s','mt0t')
          call ncinfo(ncname2( 2,:),'vxy','xy crosssections of the South-North velocity','m/s','tm0t')
          call ncinfo(ncname2( 3,:),'wxy','xy crosssections of the Vertical velocity','m/s','tt0t')
          call ncinfo(ncname2( 4,:),'thlxy','xy crosssections of the Liquid water potential temperature','K','tt0t')
          call ncinfo(ncname2( 5,:),'thvxy','xy crosssections of the Virtual potential temperature','K','tt0t')
          call ncinfo(ncname2( 6,:),'qtxy','xy crosssections of the Total water specific humidity','kg/kg','tt0t')
          call ncinfo(ncname2( 7,:),'qlxy','xy crosssections of the Liquid water specific humidity','kg/kg','tt0t')
          call ncinfo(ncname2( 8,:),'buoyxy','xy crosssection of the Buoyancy','K','tt0t')
          call ncinfo(ncname2( 9,:),'qrxy','xy crosssection of the Rain water specific humidity','kg/kg','tt0t')
          call ncinfo(ncname2(10,:),'nrxy','xy crosssection of the rain droplet number concentration','-','tt0t')
          call ncinfo(ncname2(11,:),'e120xy','xy crosssection of sqrt(turbulent kinetic energy)','m^2/s^2','tt0t')
          call open_nc(trim(output_prefix)//fname2,ncid2(cross),nrec2(cross),n1=imax,n2=jmax)
          if (nrec2(cross)==0) then
             call define_nc(ncid2(cross), 1, tncname2)
             call writestat_dims_nc(ncid2(cross))
          end if
          call define_nc(ncid2(cross), NVar, ncname2)
       end do
    end if
    if (lyz) then  ! .and. myidx == 0
       do cross=1,nyz
          if (crossortho(cross) >= 2 .and. crossortho(cross) <= i1) then
             !fname3(9:16) = cmyid
             !fname3(18:20) = cexpnr
             write(cheight,'(i4.4)') crossortho(cross) + myidx*imax
             fname3(9:12) = cheight
             fname3(15:17) = cmyidy
             fname3(19:21) = cexpnr
             call nctiminfo(tncname3(1,:))
             call ncinfo(ncname3( 1,:),'uyz','yz crosssection of the West-East velocity','m/s','0ttt')
             call ncinfo(ncname3( 2,:),'vyz','yz crosssection of the South-North velocity','m/s','0mtt')
             call ncinfo(ncname3( 3,:),'wyz','yz crosssection of the Vertical velocity','m/s','0tmt')
             call ncinfo(ncname3( 4,:),'thlyz','yz crosssection of the Liquid water potential temperature','K','0ttt')
             call ncinfo(ncname3( 5,:),'thvyz','yz crosssection of the Virtual potential temperature','K','0ttt')
             call ncinfo(ncname3( 6,:),'qtyz','yz crosssection of the Total water specific humidity','kg/kg','0ttt')
             call ncinfo(ncname3( 7,:),'qlyz','yz crosssection of the Liquid water specific humidity','kg/kg','0ttt')
             call ncinfo(ncname3( 8,:),'buoyyz','yz crosssection of the Buoyancy','K','0ttt')
             call ncinfo(ncname3( 9,:),'qryz','yz crosssection of the Rain water specific humidity','kg/kg','0ttt')
             call ncinfo(ncname3(10,:),'nryz','yz crosssection of the Number concentration','-','0ttt')
             call ncinfo(ncname3(11,:),'e120yz','yz crosssection of sqrt(turbulent kinetic energy)','m^2/s^2','0ttt')
             call open_nc(trim(output_prefix)//fname3,  ncid3(cross),nrec3(cross),n2=jmax,n3=kmax)
             if (nrec3(cross)==0) then
                call define_nc(ncid3(cross), 1, tncname3)
                call writestat_dims_nc(ncid3(cross))
             end if
             call define_nc(ncid3(cross), NVar, ncname3)
          end if
       end do
    end if
 end if


  end subroutine initcrosssection
!>Run crosssection. Mainly timekeeping
  subroutine crosssection
    use modglobal, only : rk3step,timee,dt_lim
    use modstat_nc, only : writestat_nc
#if defined(_OPENACC)
    use modgpu, only: update_host
#endif
    implicit none


    if (.not. lcross) return
    if (rk3step/=3) return
    if(timee<tnext) then
      dt_lim = min(dt_lim,tnext-timee)
      return
    end if
#if defined(_OPENACC)
    call update_host
#endif
    tnext = tnext+idtav
    dt_lim = minval((/dt_lim,tnext-timee/))

    if (lxz) call wrtvert
    if (lxy) call wrthorz
    if (lyz) call wrtorth

  end subroutine crosssection


!> Do the xz crosssections and dump them to file
  subroutine wrtvert
  use modglobal, only : imax,i1,j1,kmax,nsv,rlv,cp,rv,rd,cu,cv,cexpnr,ifoutput,rtimee
  use modfields, only : um,vm,wm,thlm,qtm,svm,thl0,qt0,ql0,e120,exnf,thvf
  use modmpi,    only : myidy
  use modstat_nc, only : lnetcdf, writestat_nc
  implicit none

  integer i,k,n
  character(20) :: name

  real, allocatable :: thv0(:,:),vars(:,:,:),buoy(:,:)

  allocate(thv0(2:i1,1:kmax),buoy(2:i1,1:kmax))

  if(lbinary .and. myidy == 0) then
     ! the old binary format only writes the first cross plane, in the first MPI row.

     do  i=2,i1
        do  k=1,kmax
           thv0(i,k) = (thl0(i,crossplane(1),k)+rlv*ql0(i,crossplane(1),k)/(cp*exnf(k))) &
                *(1+(rv/rd-1)*qt0(i,crossplane(1),k)-rv/rd*ql0(i,crossplane(1),k))
           buoy(i,k) = thv0(i,k)-thvf(k)
        enddo
     enddo

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
      do cross=1,nxz
         if (crossplane(cross) >= 2 .and. crossplane(cross) <= j1) then
            do  i=2,i1
               do  k=1,kmax
                  thv0(i,k) = (thl0(i,crossplane(cross),k)+rlv*ql0(i,crossplane(cross),k)/(cp*exnf(k))) &
                       *(1+(rv/rd-1)*qt0(i,crossplane(cross),k)-rv/rd*ql0(i,crossplane(cross),k))
                  buoy(i,k) = thv0(i,k)-thvf(k)
               enddo
            enddo
            vars(:,:,1) = um(2:i1,crossplane(cross),1:kmax)+cu
            vars(:,:,2) = vm(2:i1,crossplane(cross),1:kmax)+cv
            vars(:,:,3) = wm(2:i1,crossplane(cross),1:kmax)
            vars(:,:,4) = thlm(2:i1,crossplane(cross),1:kmax)
            vars(:,:,5) = thv0(2:i1,1:kmax)
            vars(:,:,6) = qtm(2:i1,crossplane(cross),1:kmax)
            vars(:,:,7) = ql0(2:i1,crossplane(cross),1:kmax)
            vars(:,:,8) = buoy(2:i1,1:kmax)
            if(nsv>1) then
               vars(:,:,9) = svm(2:i1,crossplane(cross),1:kmax,2)
               vars(:,:,10) = svm(2:i1,crossplane(cross),1:kmax,1)
            else
               vars(:,:,9) = 0.
               vars(:,:,10) = 0.
            end if
            vars(:,:,11) = e120(2:i1,crossplane(cross),1:kmax)
            call writestat_nc(ncid1(cross),1,tncname1,(/rtimee/),nrec1(cross),.true.)
            call writestat_nc(ncid1(cross),nvar,ncname1(1:nvar,:),vars,nrec1(cross),imax,kmax)
         end if
      end do
      deallocate(vars)
    end if
    deallocate(thv0,buoy)

  end subroutine wrtvert

!> Do the xy crosssections and dump them to file
  subroutine wrthorz
    use modglobal, only : imax,jmax,i1,j1,nsv,rlv,cp,rv,rd,cu,cv,cexpnr,ifoutput,rtimee
    use modfields, only : um,vm,wm,thlm,qtm,svm,thl0,qt0,ql0,e120,exnf,thvf
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
       allocate(vars(1:imax,1:jmax,nvar))
       do cross=1,nxy
          vars=0.
          vars(:,:,1) = um(2:i1,2:j1,crossheight(cross))+cu
          vars(:,:,2) = vm(2:i1,2:j1,crossheight(cross))+cv
          vars(:,:,3) = wm(2:i1,2:j1,crossheight(cross))
          vars(:,:,4) = thlm(2:i1,2:j1,crossheight(cross))
          vars(:,:,5) = thv0(2:i1,2:j1,cross)
          vars(:,:,6) = qtm(2:i1,2:j1,crossheight(cross))
          vars(:,:,7) = ql0(2:i1,2:j1,crossheight(cross))
          vars(:,:,8) = buoy(2:i1,2:j1,cross)
          if(nsv>1) then
             vars(:,:,9) = svm(2:i1,2:j1,crossheight(cross),iqr)
             vars(:,:,10) = svm(2:i1,2:j1,crossheight(cross),inr)
          else
             vars(:,:,9) = 0.
             vars(:,:,10) = 0.
          end if
          vars(:,:,11) = e120(2:i1,2:j1,crossheight(cross))
          call writestat_nc(ncid2(cross),1,tncname2,(/rtimee/),nrec2(cross),.true.)
          call writestat_nc(ncid2(cross),nvar,ncname2(1:nvar,:),vars,nrec2(cross),imax,jmax)
       end do
       deallocate(vars)
    end if

    deallocate(thv0,buoy)

  end subroutine wrthorz

  ! yz cross section
  subroutine wrtorth
    use modglobal, only : jmax,kmax,i1,j1,nsv,rlv,cp,rv,rd,cu,cv,cexpnr,ifoutput,rtimee
    use modfields, only : um,vm,wm,thlm,qtm,svm,thl0,qt0,ql0,e120,exnf,thvf
    use modmpi,    only : cmyid, myidx
    use modstat_nc, only : lnetcdf, writestat_nc
    implicit none


    ! LOCAL
    integer j,k,n
    character(21) :: name

    real, allocatable :: thv0(:,:),vars(:,:,:),buoy(:,:)

    allocate(thv0(1:j1,1:kmax),buoy(1:j1,1:kmax))

    if(lbinary .and. myidx == 0) then
       ! the old binary format only writes the first cross plane, in the first MPI column
       do  j=1,j1
          do  k=1,kmax
             thv0(j,k) =&
                  (thl0(crossortho(1),j,k)+&
                  rlv*ql0(crossortho(1),j,k)/&
                  (cp*exnf(k))) &
                  *(1+(rv/rd-1)*qt0(crossortho(1),j,k)&
                  -rv/rd*ql0(crossortho(1),j,k))
             buoy(j,k) =thv0(j,k)-thvf(k)
          enddo
       enddo

      open(ifoutput,file='movo_u.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((um(crossortho(1),j,k)+cu,j=2,j1),k=1,kmax)
      close(ifoutput)

      open(ifoutput,file='movo_v.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((vm(crossortho(1),j,k)+cv,j=2,j1),k=1,kmax)
      close(ifoutput)

      open(ifoutput,file='movo_w.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((wm(crossortho(1),j,k),j=2,j1),k=1,kmax)
      close(ifoutput)

      open(ifoutput,file='movo_thl.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((thlm(crossortho(1),j,k),j=2,j1),k=1,kmax)
      close(ifoutput)

      open(ifoutput,file='movo_thv.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((thv0(j,k),j=2,j1),k=1,kmax)
      close(ifoutput)

      open(ifoutput,file='movo_buoy.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((buoy(j,k),j=2,j1),k=1,kmax)
      close(ifoutput)

      open(ifoutput,file='movo_qt.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((1.e3*qtm(crossortho(1),j,k),j=2,j1),k=1,kmax)
      close(ifoutput)

      open(ifoutput,file='movo_ql.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((1.e3*ql0(crossortho(1),j,k),j=2,j1),k=1,kmax)
      close(ifoutput)

      do n = 1,nsv
        name = 'movh_tnn.'//cmyid//'.'//cexpnr
        write(name(7:8),'(i2.2)') n
        open(ifoutput,file=name,position='append',action='write')
        write(ifoutput,'(es12.5)') ((svm(crossortho(1),j,k,n),j=2,j1),k=1,kmax)
        close(ifoutput)
      end do
    end if

    if (lnetcdf) then
       allocate(vars(1:jmax,1:kmax,nvar))
       do cross=1,nyz
          if (crossortho(cross) >= 2 .and. crossortho(cross) <= i1) then
             do  j=1,j1
                do  k=1,kmax
                   thv0(j,k) =&
                        (thl0(crossortho(cross),j,k)+&
                        rlv*ql0(crossortho(cross),j,k)/&
                        (cp*exnf(k))) &
                        *(1+(rv/rd-1)*qt0(crossortho(cross),j,k)&
                        -rv/rd*ql0(crossortho(cross),j,k))
                   buoy(j,k) =thv0(j,k)-thvf(k)
                enddo
             enddo
             vars(:,:,1) = um(crossortho(cross),2:j1,1:kmax)+cu
             vars(:,:,2) = vm(crossortho(cross),2:j1,1:kmax)+cv
             vars(:,:,3) = wm(crossortho(cross),2:j1,1:kmax)
             vars(:,:,4) = thlm(crossortho(cross),2:j1,1:kmax)
             vars(:,:,5) = thv0(2:j1,1:kmax)
             vars(:,:,6) = qtm(crossortho(cross),2:j1,1:kmax)
             vars(:,:,7) = ql0(crossortho(cross),2:j1,1:kmax)
             vars(:,:,8) = buoy(2:j1,1:kmax)
             if(nsv>1) then
                vars(:,:,9) = svm(crossortho(cross),2:j1,1:kmax,2)
                vars(:,:,10) = svm(crossortho(cross),2:j1,1:kmax,1)
             else
                vars(:,:,9) = 0.
                vars(:,:,10) = 0.
             end if
             vars(:,:,11) = e120(crossortho(cross),2:j1,1:kmax)
             call writestat_nc(ncid3(cross),1,tncname3,(/rtimee/),nrec3(cross),.true.)
             call writestat_nc(ncid3(cross),nvar,ncname3(1:nvar,:),vars,nrec3(cross),jmax,kmax)
          end if
       end do
       deallocate(vars)
    end if

    deallocate(thv0,buoy)

  end subroutine wrtorth

!> Clean up when leaving the run
  subroutine exitcrosssection
    use modstat_nc, only : exitstat_nc,lnetcdf
    use modmpi, only : myidx, myidy
    use modglobal, only : i1,j1
    implicit none

    if(lcross .and. lnetcdf) then
       if (lxz) then
          do cross=1,nxz
             if (crossplane(cross) >= 2 .and. crossplane(cross) <= j1) then
                call exitstat_nc(ncid1(cross))
             end if
          end do
       end if
       if (lxy) then
          do cross=1,nxy
             call exitstat_nc(ncid2(cross))
          end do
       end if
       if (lyz) then
          do cross=1,nyz
             if (crossortho(cross) >= 2 .and. crossortho(cross) <= i1) then
                call exitstat_nc(ncid3(cross))
             end if
          end do
       end if
    end if

  end subroutine exitcrosssection

end module modcrosssection
