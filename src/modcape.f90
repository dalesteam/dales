!> \file modcape.f90
!!   Dumps cross sections of CAPEreal, CAPEmax, CINlower,CINupper, CINmax, ALE, W2max, QTcloudbase, THLcloudbase, Wcloudbase, THVcloudbase, QLcloudbase, LWP, RWP, cloud top
!
!>
!! Crosssections in the xy-plane
!! If netcdf is true, this module leads the capereal.myid.expnr.nc output

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
module modcape

  use modglobal, only : longint,kmax

implicit none
private
PUBLIC :: initcape,docape,exitcape
save
!NetCDF variables
  integer,parameter :: nvar = 15
  integer :: ncid = 0
  character(80) :: fname = 'cape.xxx.xxx.nc'
  character(80),dimension(nvar,4) :: ncname
  character(80),dimension(1,4) :: tncname
  integer :: nrec = 0
  real    :: dtav
  integer(kind=longint) :: idtav,tnext
  logical :: lcape = .false. !< switch for doing the crosssection (on/off)

contains

!> Initializing cape crossections. Read out the namelist, initializing the variables
  subroutine initcape
    use modmpi,   only :myid,my_real,mpierr,comm3d,mpi_logical,mpi_integer,cmyid
    use modglobal,only :imax,jmax,ifnamopt,fname_options,dtmax,rk3step, dtav_glob,ladaptive,j1,kmax,i1,dt_lim,cexpnr,tres,btime
    use modstat_nc,only : lnetcdf,open_nc, define_nc, redefine_nc,ncinfo,writestat_dims_nc
   implicit none

    integer :: ierr,k

    namelist/NAMCAPE/ &
    lcape, dtav

    dtav = dtav_glob
    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMCAPE,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions NAMCAPE'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMCAPE'
      endif
      write(6 ,NAMCAPE)
      close(ifnamopt)
    end if

    call MPI_BCAST(dtav       ,1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(lcape     ,1,MPI_LOGICAL,0,comm3d,mpierr)

    idtav = dtav/tres
    tnext   = idtav+btime
    if(.not.(lcape)) return
    dt_lim = min(dt_lim,tnext)

    !CAPEreal, CAPEmax, CINreal, CINmax, ALE, W2max, QTcloudbase, THLcloudbase, Wcloudbase, THVcloudbase, QLcloudbase, LWP, RWP, cloud top

    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'cape: dtav should be a integer multiple of dtmax'
    end if
    if (lnetcdf) then
    fname(6:8) = cmyid
    fname(10:12) = cexpnr
    call ncinfo(tncname(1,:),'time','Time','s','time')
    call ncinfo(ncname( 1,:),'capereal','xy crosssections of actual capereal','J/m^2','tt0t')
    call ncinfo(ncname( 2,:),'cinlower','xy crosssections of actual CIN in CINmax-region','J/m^2','tt0t')
    call ncinfo(ncname( 3,:),'cinupper','xy crosssections of actual CIN above','J/m^2','tt0t')
    call ncinfo(ncname( 4,:),'capemax','xy crosssections of CAPEmax','J/m^2','tt0t')
    call ncinfo(ncname( 5,:),'cinmax','xy crosssections of CIN as in CAPEmax','J/m^2','tt0t')
    call ncinfo(ncname( 6,:),'hw2cb','xy crosssections of 1/2 W^2 at the top of the subcloud layer','m^2/s^2','tt0t')
    call ncinfo(ncname( 7,:),'hw2max','xy crosssections of highest 1/2 W^2','m^2/s^2','tt0t')
    call ncinfo(ncname( 8,:),'qtcb','xy crosssections of qt at cloudbase','kg/kg','tt0t')
    call ncinfo(ncname( 9,:),'thlcb','xy crosssections of thl at cloudbase','K','tt0t')
    call ncinfo(ncname( 10,:),'wcb','xy crosssections of w at cloudbase','m/s','tt0t')
    call ncinfo(ncname( 11,:),'thvcb','xy crosssections thv at cloudbase','K','tt0t')
    call ncinfo(ncname( 12,:),'qlcb','xy crosssections ql at cloudbase','kg/kg','tt0t')
    call ncinfo(ncname( 13,:),'lwp','xy crosssections liquid water path','kg/m^2','tt0t')
    call ncinfo(ncname( 14,:),'rwp','xy crosssections rain water path','kg/m^2','tt0t')
    call ncinfo(ncname( 15,:),'cldtop','xy crosssections cloud top height','m','tt0t')
    call open_nc(fname,  ncid,n1=imax,n2=jmax)
    call define_nc( ncid, 1, tncname)
    call writestat_dims_nc(ncid)
    call redefine_nc(ncid)
    call define_nc(ncid, NVar, ncname)
    end if

  end subroutine initcape

!>Run crosssection.
  subroutine docape
    use modglobal, only : imax,jmax,i1,j1,k1,kmax,nsv,rlv,cp,rv,rd,cu,cv,cexpnr,ifoutput,rk3step,timee,rtimee,dt_lim,grav,eps1,nsv,ttab,esatltab,esatitab,zf,dzf,tup,tdn
    use modfields, only : thl0,qt0,ql0,w0,sv0,exnf,thvf,exnf,presf
    use modmpi,    only : cmyid
    use modstat_nc, only : lnetcdf, writestat_nc
    use modgenstat, only : qlmnlast,wtvtmnlast
    use modmicrodata, only : iqr
    implicit none

    real, allocatable :: capereal(:,:),cinlower(:,:),cinupper(:,:),capemax(:,:),cinmax(:,:),hw2cb(:,:),hw2max(:,:),qtcb(:,:),thlcb(:,:),wcb(:,:),thvcb(:,:),qlcb(:,:),lwp(:,:),rwp(:,:),cldtop(:,:)
    real, allocatable :: thvfull(:,:,:),thvma(:,:,:),qlma(:,:,:),vars(:,:,:)
    integer, allocatable :: capetop(:,:),matop(:,:)
    logical,allocatable :: capemask(:,:,:)

    ! LOCAL VARIABLES
    integer :: i,j,k,kcb,ktest,tlonr,thinr,niter,nitert
    character(40) :: name
    real :: Tnr,Tnr_old,ilratio,tlo,thi,esl1,esi1,qsatur,thlguess,thlguessmin,ttry,qvsl1,qvsi1

    if (.not. lcape) return
    if (rk3step/=3) return
    if(timee<tnext) then
      dt_lim = min(dt_lim,tnext-timee)
      return
    end if
    tnext = tnext+idtav
    dt_lim = minval((/dt_lim,tnext-timee/))

    allocate(capereal(2:i1,2:j1),cinlower(2:i1,2:j1),cinupper(2:i1,2:j1))
    allocate(capemax(2:i1,2:j1),cinmax(2:i1,2:j1),hw2cb(2:i1,2:j1))
    allocate(hw2max(2:i1,2:j1),qtcb(2:i1,2:j1),thlcb(2:i1,2:j1),wcb(2:i1,2:j1))
    allocate(thvcb(2:i1,2:j1),qlcb(2:i1,2:j1),lwp(2:i1,2:j1),rwp(2:i1,2:j1),cldtop(2:i1,2:j1))
    allocate(thvfull(2:i1,2:j1,1:k1),thvma(2:i1,2:j1,1:k1),qlma(2:i1,2:j1,1:k1),capemask(2:i1,2:j1,1:k1),capetop(2:i1,2:j1),matop(2:i1,2:j1))

    ! DETERMINE CLOUD BASE, UNFORTUNATELY HAVE TO USE STATS HERE: END UP JUST BELOW
    kcb=1
    ! find robust minimum of buoyancy flux (determined at half-level!)
    if(any(abs(wtvtmnlast)>1e-10)) then
      do k=kmax-2,1,-1
        if ((wtvtmnlast(k)<wtvtmnlast(k-1)).and.(wtvtmnlast(k)<wtvtmnlast(k+1)).and.(wtvtmnlast(k)<wtvtmnlast(k+2))) then
          ktest = k
        end if
      end do
    end if

    ! Take highest half-level below which it is non-cloudy
    if (ktest>0) then
      do k=2,ktest
        if(qlmnlast(k-1)<0.01) then
        kcb=k
        end if
      end do
    end if

    ! loops over i,j,k
    lwp=0.
    rwp=0.
    cldtop=0.
    hw2max=0.
    do k=1,k1
    do j=2,j1
    do i=2,i1
      capemask(i,j,k)=.false. ! reset
      qlma(i,j,k)=0. ! reset
      thvfull(i,j,k)=(thl0(i,j,k)+rlv*ql0(i,j,k)/(cp*exnf(k))) &
                      *(1+(rv/rd-1)*qt0(i,j,k)-rv/rd*ql0(i,j,k))
      lwp(i,j)=lwp(i,j)+ql0(i,j,k)*dzf(k)
      if (ql0(i,j,k) > eps1) then
        cldtop(i,j)=zf(k)
      endif
      if (w0(i,j,k)**2 > hw2max(i,j)) then
        hw2max(i,j)=0.5*w0(i,j,k)*abs(w0(i,j,k))
      endif
    enddo
    enddo
    enddo

    if(nsv>1) then
      do k=1,k1
      do j=2,j1
      do i=2,i1
      rwp(i,j)=rwp(i,j)+sv0(i,j,k,iqr)*dzf(k)
      enddo
      enddo
      enddo
    endif

    ! Cloud base level quantities and reset tops
    do  j=2,j1
    do  i=2,i1
    thlcb(i,j)=thl0(i,j,kcb)
    thvcb(i,j)=thvfull(i,j,kcb)
    wcb(i,j)=(w0(i,j,kcb)+w0(i,j,kcb+1))/2.
    hw2cb(i,j)=0.5*wcb(i,j)*abs(wcb(i,j))
    qtcb(i,j)=qt0(i,j,kcb)
    qlcb(i,j)=ql0(i,j,kcb)
    capetop(i,j)=0
    matop(i,j)=0
    enddo
    enddo

    !calculate moist adiabat from cloud base: let pressure adjust to slab mean

    do  k=1,k1 
    do  j=2,j1
    do  i=2,i1
    ! full level
      if(matop(i,j)==0) then
        Tnr=exnf(k)*thlcb(i,j)+(rlv/cp)*qlma(i,j,k-1) ! First guess for full level, use ql from below
        Tnr_old=0.
          do while (abs(Tnr-Tnr_old) > 0.002) ! Find T at first level
            niter = niter+1
            Tnr_old=Tnr
            ilratio = max(0.,min(1.,(Tnr-tdn)/(tup-tdn)))
            tlonr=int((Tnr-150.)*5.)
            thinr=tlonr+1
            tlo=ttab(tlonr)
            thi=ttab(thinr)
            esl1=(thi-Tnr)*5.*esatltab(tlonr)+(Tnr-tlo)*5.*esatltab(thinr)
            esi1=(thi-Tnr)*5.*esatitab(tlonr)+(Tnr-tlo)*5.*esatitab(thinr)
            qsatur = ilratio*(rd/rv)*esl1/(presf(k)-(1.-rd/rv)*esl1)+(1.-ilratio)*(rd/rv)*esi1/(presf(k)-(1.-rd/rv)*esi1)
            thlguess = Tnr/exnf(k)-(rlv/(cp*exnf(k)))*max(qtcb(i,j)-qsatur,0.)
    
            ttry=Tnr-0.002
            ilratio = max(0.,min(1.,(ttry-tdn)/(tup-tdn)))
            tlonr=int((Tnr-150.)*5.)
            thinr=tlonr+1
            tlo=ttab(tlonr)
            thi=ttab(thinr)
            esl1=(thi-ttry)*5.*esatltab(tlonr)+(ttry-tlo)*5.*esatltab(thinr)
            esi1=(thi-ttry)*5.*esatitab(tlonr)+(ttry-tlo)*5.*esatitab(thinr)
            qsatur = ilratio*(rd/rv)*esl1/(presf(k)-(1.-rd/rv)*esl1)+(1.-ilratio)*(rd/rv)*esi1/(presf(k)-(1.-rd/rv)*esi1)
            thlguessmin = ttry/exnf(k)-(rlv/(cp*exnf(k)))*max(qtcb(i,j)-qsatur,0.)
    
            Tnr = Tnr - (thlguess-thlcb(i,j))/((thlguess-thlguessmin)*500.)
          enddo
        nitert =max(nitert,niter)
        niter = 0
        ilratio = max(0.,min(1.,(Tnr-tdn)/(tup-tdn)))
        tlonr=int((Tnr-150.)*5.)
        thinr=tlonr+1
        tlo=ttab(tlonr)
        thi=ttab(thinr)
        esl1=(thi-Tnr)*5.*esatltab(tlonr)+(Tnr-tlo)*5.*esatltab(thinr)
        esi1=(thi-Tnr)*5.*esatitab(tlonr)+(Tnr-tlo)*5.*esatitab(thinr)
        qvsl1=rd/rv*esl1/(presf(k)-(1.-rd/rv)*esl1)
        qvsi1=rd/rv*esi1/(presf(k)-(1.-rd/rv)*esi1)
        qsatur = ilratio*qvsl1+(1.-ilratio)*qvsi1
        qlma(i,j,k) = max(qtcb(i,j)-qsatur,0.)
        thvma(i,j,k)=(thlcb(i,j)+(rlv*qlma(i,j,k))/(cp*exnf(k)))*(1+(rv/rd-1)*qtcb(i,j)-rv/rd*qlma(i,j,k)) ! calculate thv, assuming thl conserved
        if(thvma(i,j,k)<thvf(k)-10) then
          matop(i,j)=k
        endif
      endif
    enddo
    enddo
    enddo

    ! calculate top of moist adiabat for capereal calculations
    do k=kcb,k1 
    do j=2,j1
    do i=2,i1
      if(k<matop(i,j)) then
        if(thvma(i,j,k)>thvf(k)) then
        capetop(i,j)=k
        capemask(i,j,k)=.true.
        endif
      endif
    enddo
    enddo
    enddo

    ! capereal and CIN of Moist Adiabat
    ! no nice interpolation yet
    do j=2,j1
    do i=2,i1
    capemax(i,j)=0.
    cinmax(i,j)=0.
    capereal(i,j)=0.
    cinlower(i,j)=0.
    cinupper(i,j)=0.
    enddo
    enddo

    do k=kcb,k1
    do j=2,j1
    do i=2,i1
      if(k<capetop(i,j)) then
        cinmax(i,j)=cinmax(i,j)+max(grav*dzf(k)*(thvf(k)-thvma(i,j,k))/thvf(k),0.)
        capemax(i,j)=capemax(i,j)+max(grav*dzf(k)*(thvma(i,j,k)-thvf(k))/thvf(k),0.)
        if(capemask(i,j,k).eqv..true.) then
          capereal(i,j)=capereal(i,j)+max(grav*dzf(k)*(thvfull(i,j,k)-thvf(k))/thvf(k),0.)
          cinupper(i,j)=cinupper(i,j)+max(grav*dzf(k)*(thvf(k)-thvfull(i,j,k))/thvf(k),0.)
        else
          capereal(i,j)=capereal(i,j)+max(grav*dzf(k)*(thvfull(i,j,k)-thvf(k))/thvf(k),0.)
          cinlower(i,j)=cinlower(i,j)+max(grav*dzf(k)*(thvf(k)-thvfull(i,j,k))/thvf(k),0.)
        endif
      endif
    enddo
    enddo
    enddo

!     open(ifoutput,file='XXX'//cmyid//'.'//cexpnr,position='append',action='write')
!     write(ifoutput,'(es12.5)') ((XXX(i,j),i=2,i1),j=2,j1)
!     close(ifoutput)
!     open(ifoutput,file='XXX'//cmyid//'.'//cexpnr,position='append',action='write')
!     write(ifoutput,'(es12.5)') ((XXX(i,j),i=2,i1),j=2,j1)
!     close(ifoutput)
!     open(ifoutput,file='XXX'//cmyid//'.'//cexpnr,position='append',action='write')
!     write(ifoutput,'(es12.5)') ((XXX(i,j),i=2,i1),j=2,j1)
!     close(ifoutput)
!     open(ifoutput,file='XXX'//cmyid//'.'//cexpnr,position='append',action='write')
!     write(ifoutput,'(es12.5)') ((XXX(i,j),i=2,i1),j=2,j1)
!     close(ifoutput)
!     open(ifoutput,file='XXX'//cmyid//'.'//cexpnr,position='append',action='write')
!     write(ifoutput,'(es12.5)') ((XXX(i,j),i=2,i1),j=2,j1)
!     close(ifoutput)
!     open(ifoutput,file='XXX'//cmyid//'.'//cexpnr,position='append',action='write')
!     write(ifoutput,'(es12.5)') ((XXX(i,j),i=2,i1),j=2,j1)
!     close(ifoutput)
!     open(ifoutput,file='XXX'//cmyid//'.'//cexpnr,position='append',action='write')
!     write(ifoutput,'(es12.5)') ((XXX(i,j),i=2,i1),j=2,j1)
!     close(ifoutput)
!     open(ifoutput,file='XXX'//cmyid//'.'//cexpnr,position='append',action='write')
!     write(ifoutput,'(es12.5)') ((XXX(i,j),i=2,i1),j=2,j1)
!     close(ifoutput)
!     open(ifoutput,file='XXX'//cmyid//'.'//cexpnr,position='append',action='write')
!     write(ifoutput,'(es12.5)') ((XXX(i,j),i=2,i1),j=2,j1)
!     close(ifoutput)
!     open(ifoutput,file='XXX'//cmyid//'.'//cexpnr,position='append',action='write')
!     write(ifoutput,'(es12.5)') ((XXX(i,j),i=2,i1),j=2,j1)
!     close(ifoutput)
!     open(ifoutput,file='XXX'//cmyid//'.'//cexpnr,position='append',action='write')
!     write(ifoutput,'(es12.5)') ((XXX(i,j),i=2,i1),j=2,j1)
!     close(ifoutput)
!     open(ifoutput,file='XXX'//cmyid//'.'//cexpnr,position='append',action='write')
!     write(ifoutput,'(es12.5)') ((XXX(i,j),i=2,i1),j=2,j1)
!     close(ifoutput)

    if (lnetcdf) then
      allocate(vars(1:imax,1:jmax,15))
      vars(:,:,1) = capereal(2:i1,2:j1)
      vars(:,:,2) = cinlower(2:i1,2:j1)
      vars(:,:,3) = cinupper(2:i1,2:j1)
      vars(:,:,4) = capemax(2:i1,2:j1)
      vars(:,:,5) = cinmax(2:i1,2:j1)
      vars(:,:,6) = hw2cb(2:i1,2:j1)
      vars(:,:,7) = hw2max(2:i1,2:j1)
      vars(:,:,8) = qtcb(2:i1,2:j1)
      vars(:,:,9) = thlcb(2:i1,2:j1)
      vars(:,:,10) = wcb(2:i1,2:j1)
      vars(:,:,11) = thvcb(2:i1,2:j1)
      vars(:,:,12) = qlcb(2:i1,2:j1)
      vars(:,:,13) = lwp(2:i1,2:j1)
      vars(:,:,14) = rwp(2:i1,2:j1)
      vars(:,:,15) = cldtop(2:i1,2:j1)
      call writestat_nc(ncid,1,tncname,(/rtimee/),nrec,.true.)
      call writestat_nc(ncid,15,ncname(1:15,:),vars,nrec,imax,jmax)
      deallocate(vars)
    end if

    deallocate(capereal,cinlower,cinupper,capemax,cinmax,hw2cb,hw2max,qtcb,thlcb,wcb,thvcb,qlcb,lwp,rwp,cldtop,thvfull,thvma,qlma,capemask,capetop,matop)

  end subroutine docape

!> Clean up when leaving the run
  subroutine exitcape
    use modstat_nc, only : exitstat_nc,lnetcdf
    use modmpi, only : myid
    implicit none

    if(lcape .and. lnetcdf) then
    call exitstat_nc(ncid)
    end if

  end subroutine exitcape

end module modcape
