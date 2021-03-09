!> \file modcape.f90
!!   Dumps cross sections of dcape, CAPEmax, CINlower,CINupper, CINmax, ALE, W2max,
!!   QTcloudbase, THLcloudbase, Wcloudbase, THVcloudbase, QLcloudbase, LWP, RWP, cloud top
!
!>
!! Crosssections in the xy-plane
!! If netcdf is true, this module leads the dcape.myid.expnr.nc output

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
  integer,parameter :: nvar = 20
  integer :: ncid4 = 0
  integer :: nrec = 0
  character(80) :: fname = 'cape.xxxxyxxx.xxx.nc'
  character(80),dimension(nvar,4) :: ncname
  character(80),dimension(1,4) :: tncname
  real    :: dtav
  integer(kind=longint) :: idtav,tnext
  logical :: lcape = .false. !< switch for doing the crosssection (on/off)

contains

!> Initializing cape crossections. Read out the namelist, initializing the variables
  subroutine initcape
    use modmpi,   only :myid,mpierr,comm3d,mpi_logical,cmyid,D_MPI_BCAST
    use modglobal,only :imax,jmax,ifnamopt,fname_options,dtmax,dtav_glob,ladaptive,dt_lim,cexpnr,tres,btime,checknamelisterror
    use modstat_nc,only : lnetcdf,open_nc, define_nc, redefine_nc,ncinfo,nctiminfo,writestat_dims_nc
   implicit none

    integer :: ierr

    namelist/NAMCAPE/ &
    lcape, dtav

    dtav = dtav_glob
    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMCAPE,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMCAPE')
      write(6 ,NAMCAPE)
      close(ifnamopt)
    end if

    call D_MPI_BCAST(dtav    ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(lcape   ,1,0,comm3d,mpierr)

    idtav = dtav/tres
    tnext   = idtav+btime
    if(.not.(lcape)) return
    dt_lim = min(dt_lim,tnext)

    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'cape: dtav should be a integer multiple of dtmax'
    end if
    if (lnetcdf) then
    fname(6:13) = cmyid
    fname(15:17) = cexpnr
    call nctiminfo(tncname(1,:))
    call ncinfo(ncname( 1,:),'dcape','xy crosssections of actual dcape','J/m^2','tt0t')
    call ncinfo(ncname( 2,:),'dscape','xy crosssections of actual dscape','J/m^2','tt0t')
    call ncinfo(ncname( 3,:),'dcin','xy crosssections of actual CIN between zcb and 1.2 zcb','J/m^2','tt0t')
    call ncinfo(ncname( 4,:),'dscin','xy crosssections of actual CIN up to zcb','J/m^2','tt0t')
    call ncinfo(ncname( 5,:),'dcintot','xy crosssections of total CIN up to dcape level','J/m^2','tt0t')
    call ncinfo(ncname( 6,:),'capemax','xy crosssections of CAPEmax','J/m^2','tt0t')
    call ncinfo(ncname( 7,:),'cinmax','xy crosssections of CIN as in CAPEmax','J/m^2','tt0t')
    call ncinfo(ncname( 8,:),'hw2cb','xy crosssections of 1/2 W^2 at the top of the subcloud layer','m^2/s^2','tt0t')
    call ncinfo(ncname( 9,:),'hw2max','xy crosssections of highest 1/2 W^2','m^2/s^2','tt0t')
    call ncinfo(ncname( 10,:),'qtcb','xy crosssections of qt at cloudbase','kg/kg','tt0t')
    call ncinfo(ncname( 11,:),'thlcb','xy crosssections of thl at cloudbase','K','tt0t')
    call ncinfo(ncname( 12,:),'wcb','xy crosssections of w at cloudbase','m/s','tt0t')
    call ncinfo(ncname( 13,:),'buoycb','xy crosssections buoyancy at cloudbase','K','tt0t')
    call ncinfo(ncname( 14,:),'buoymax','xy crosssections maximum buoyancy','K','tt0t')
    call ncinfo(ncname( 15,:),'qlcb','xy crosssections ql at cloudbase','kg/kg','tt0t')
    call ncinfo(ncname( 16,:),'lwp','xy crosssections liquid water path','kg/m^2','tt0t')
    call ncinfo(ncname( 17,:),'rwp','xy crosssections rain water path','kg/m^2','tt0t')
    call ncinfo(ncname( 18,:),'twp','total water path','kg/m^2','tt0t')
    call ncinfo(ncname( 19,:),'cldtop','xy crosssections cloud top height','m','tt0t')
    call ncinfo(ncname( 20,:),'surfprec','surface precipitation','-','tt0t')
    call open_nc(fname,  ncid4,nrec,n1=imax,n2=jmax)
    if (nrec==0) then
      call define_nc( ncid4, 1, tncname)
      call writestat_dims_nc(ncid4)
    end if
    call define_nc( ncid4, NVar, ncname)
    end if

  end subroutine initcape

!>Run crosssection.
  subroutine docape
    use modglobal, only : imax,jmax,i1,j1,k1,kmax,nsv,rlv,cp,rv,rd,rk3step,timee,rtimee,dt_lim,grav,eps1,&
    nsv,ttab,esatltab,esatitab,zf,dzf,tup,tdn,zh,kcb
    use modfields, only : thl0,qt0,ql0,w0,sv0,exnf,thvf,exnf,presf,rhobf
    use modstat_nc, only : lnetcdf, writestat_nc
    use modgenstat, only : qlmnlast,wthvtmnlast
    use modmicrodata, only : iqr, precep, imicro
    use modmpi
    implicit none

    real, allocatable :: dcape(:,:),dscape(:,:),dcin(:,:),dscin(:,:),dcintot(:,:),capemax(:,:),&
    cinmax(:,:),hw2cb(:,:),hw2max(:,:),qtcb(:,:),&
    thlcb(:,:),wcb(:,:),buoycb(:,:),buoymax(:,:),qlcb(:,:),lwp(:,:),twp(:,:),rwp(:,:),&
    cldtop(:,:),thl200400(:,:),qt200400(:,:),sprec(:,:)
    real, allocatable :: thvfull(:,:,:),thvma(:,:,:),qlma(:,:,:),vars(:,:,:)
    integer, allocatable :: capetop(:,:),matop(:,:)
    logical,allocatable :: capemask(:,:,:)

    ! LOCAL VARIABLES
    integer :: i,j,k,ktest,tlonr,thinr,niter,nitert,kdmax,kdmaxl
    real :: Tnr,Tnr_old,ilratio,tlo,thi,esl1,esi1,qsatur,thlguess,thlguessmin,ttry,qvsl1,qvsi1

    if (.not. lcape) return
    if (rk3step/=3) return
    if(timee<tnext) then
      dt_lim = min(dt_lim,tnext-timee)
      return
    end if
    tnext = tnext+idtav
    dt_lim = minval((/dt_lim,tnext-timee/))

    allocate(dcape(2:i1,2:j1),dscape(2:i1,2:j1),dcin(2:i1,2:j1),dscin(2:i1,2:j1),dcintot(2:i1,2:j1))
    allocate(capemax(2:i1,2:j1),cinmax(2:i1,2:j1),hw2cb(2:i1,2:j1))
    allocate(thl200400(2:i1,2:j1),qt200400(2:i1,2:j1))
    allocate(hw2max(2:i1,2:j1),qtcb(2:i1,2:j1),thlcb(2:i1,2:j1),wcb(2:i1,2:j1))
    allocate(buoycb(2:i1,2:j1),buoymax(2:i1,2:j1),qlcb(2:i1,2:j1),&
             lwp(2:i1,2:j1),rwp(2:i1,2:j1),cldtop(2:i1,2:j1),twp(2:i1,2:j1))
    allocate(thvfull(2:i1,2:j1,1:k1),thvma(2:i1,2:j1,1:k1),qlma(2:i1,2:j1,1:k1),&
             capemask(2:i1,2:j1,1:k1),capetop(2:i1,2:j1),matop(2:i1,2:j1),sprec(2:i1,2:j1))

    ! DETERMINE CLOUD BASE, UNFORTUNATELY HAVE TO USE STATS HERE: END UP JUST BELOW
    kcb=1
    ktest=0
    ! find robust minimum of buoyancy flux (determined at half-level!)
    if(any(abs(wthvtmnlast)>1e-10)) then
      do k=kmax-7,2,-1
        if ((wthvtmnlast(k)<wthvtmnlast(k-1)).and.(wthvtmnlast(k)<wthvtmnlast(k+1)).and.&
           (wthvtmnlast(k)<wthvtmnlast(k+2)).and.(wthvtmnlast(k)<wthvtmnlast(k+3)).and.&
           (wthvtmnlast(k)<wthvtmnlast(k+4))) then
          ktest = k
        end if
      end do
    end if

    ! Take highest half-level below which it is non-cloudy
    if (ktest>0) then
      do k=2,ktest
        if(qlmnlast(k-1)<0.001) then
        kcb=k
        end if
      end do
    end if

    kdmaxl=kcb

    ! loops over i,j,k
    lwp=0.
    twp=0.
    rwp=0.
    sprec=0.
    cldtop=0.
    hw2max=0.
    buoymax=0.
    thl200400=0.
    qt200400=0.
    do k=1,k1
    do j=2,j1
    do i=2,i1
      capemask(i,j,k)=.false. ! reset
      qlma(i,j,k)=0. ! reset
      thvma(i,j,k)=0.
      thvfull(i,j,k)=(thl0(i,j,k)+rlv*ql0(i,j,k)/(cp*exnf(k))) &
                      *(1+(rv/rd-1)*qt0(i,j,k)-rv/rd*ql0(i,j,k))
      lwp(i,j)=lwp(i,j)+rhobf(k)*ql0(i,j,k)*dzf(k)
      twp(i,j)=twp(i,j)+rhobf(k)*qt0(i,j,k)*dzf(k)
      if (ql0(i,j,k) > eps1) then
        cldtop(i,j)=zf(k)
      endif
      if ((ql0(i,j,k)>eps1) .and. (thvfull(i,j,k)-thvf(k)>0)) then
        kdmaxl=max(kdmaxl,k) ! maximum height of buoyant cloud cores
      endif
      if (w0(i,j,k)**2 > hw2max(i,j)) then
        hw2max(i,j)=0.5*w0(i,j,k)*abs(w0(i,j,k))
      endif
      if ((thvfull(i,j,k)-thvf(k))>buoymax(i,j)) then
        buoymax(i,j)=thvfull(i,j,k)-thvf(k)
      endif
      if ((k>1).and.(k<k1)) then
      if ((zh(k+1)>200.).and.(zh(k)<400.)) then
        thl200400(i,j)=thl200400(i,j)+thl0(i,j,k)*max(min(zh(k+1),400.)-max(zh(k),200.),0.)
        qt200400(i,j)=qt200400(i,j)+qt0(i,j,k)*max(min(zh(k+1),400.)-max(zh(k),200.),0.)
      endif
      endif
    enddo
    enddo
    enddo

    call D_MPI_ALLREDUCE(kdmaxl,kdmax,1,MPI_MAX,comm3d,mpierr)

    do j=2,j1
    do i=2,i1
        thl200400(i,j)=thl200400(i,j)/200.
        qt200400(i,j)=qt200400(i,j)/200.
    enddo
    enddo

    if(nsv>1) then
    if(imicro>0) then
      do k=1,k1
      do j=2,j1
      do i=2,i1
      rwp(i,j)=rwp(i,j)+rhobf(k)*sv0(i,j,k,iqr)*dzf(k)
      enddo
      enddo
      enddo
      do j=2,j1
      do i=2,i1
      sprec(i,j)=precep(i,j,1)*rhobf(1) ! correct for density to find total rain mass-flux
      enddo
      enddo
    endif
    endif

    ! Cloud base level quantities and reset tops
    do  j=2,j1
    do  i=2,i1
    thlcb(i,j)=thl0(i,j,kcb)
    buoycb(i,j)=thvfull(i,j,kcb)-thvf(kcb)
    wcb(i,j)=(w0(i,j,kcb)+w0(i,j,kcb+1))/2.
    hw2cb(i,j)=0.5*wcb(i,j)*abs(wcb(i,j))
    qtcb(i,j)=qt0(i,j,kcb)
    qlcb(i,j)=ql0(i,j,kcb)
    capetop(i,j)=1
    matop(i,j)=0
    enddo
    enddo
    !calculate moist adiabat from surface, rather than cloud base: let pressure adjust to slab mean

    nitert=0

    k=1
    do j=2,j1
    do i=2,i1
    ! full level
        Tnr=exnf(k)*thl200400(i,j) ! First guess for full level, use no ql from below
        Tnr_old=0.
        niter = 0
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
            thlguess = Tnr/exnf(k)-(rlv/(cp*exnf(k)))*max(qt200400(i,j)-qsatur,0.)

            ttry=Tnr-0.002
            ilratio = max(0.,min(1.,(ttry-tdn)/(tup-tdn)))
            tlonr=int((Tnr-150.)*5.)
            thinr=tlonr+1
            tlo=ttab(tlonr)
            thi=ttab(thinr)
            esl1=(thi-ttry)*5.*esatltab(tlonr)+(ttry-tlo)*5.*esatltab(thinr)
            esi1=(thi-ttry)*5.*esatitab(tlonr)+(ttry-tlo)*5.*esatitab(thinr)
            qsatur = ilratio*(rd/rv)*esl1/(presf(k)-(1.-rd/rv)*esl1)+(1.-ilratio)*(rd/rv)*esi1/(presf(k)-(1.-rd/rv)*esi1)
            thlguessmin = ttry/exnf(k)-(rlv/(cp*exnf(k)))*max(qt200400(i,j)-qsatur,0.)

            Tnr = Tnr - (thlguess-thl200400(i,j))/((thlguess-thlguessmin)*500.)
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
        qlma(i,j,k) = max(qt200400(i,j)-qsatur,0.)
        thvma(i,j,k)=(thl200400(i,j)+(rlv*qlma(i,j,k))/(cp*exnf(k)))&
                     *(1+(rv/rd-1)*qt200400(i,j)-rv/rd*qlma(i,j,k)) ! calculate thv
    enddo
    enddo
    do k=2,k1
    do j=2,j1
    do i=2,i1
    ! full level
      if(matop(i,j)==0) then
        Tnr=exnf(k)*thl200400(i,j)+(rlv/cp)*qlma(i,j,k-1) ! Guess for full level
        Tnr_old=0.
        niter=0
          do while (abs(Tnr-Tnr_old) > 0.002) ! Find T at level
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
            thlguess = Tnr/exnf(k)-(rlv/(cp*exnf(k)))*max(qt200400(i,j)-qsatur,0.)

            ttry=Tnr-0.002
            ilratio = max(0.,min(1.,(ttry-tdn)/(tup-tdn)))
            tlonr=int((Tnr-150.)*5.)
            thinr=tlonr+1
            tlo=ttab(tlonr)
            thi=ttab(thinr)
            esl1=(thi-ttry)*5.*esatltab(tlonr)+(ttry-tlo)*5.*esatltab(thinr)
            esi1=(thi-ttry)*5.*esatitab(tlonr)+(ttry-tlo)*5.*esatitab(thinr)
            qsatur = ilratio*(rd/rv)*esl1/(presf(k)-(1.-rd/rv)*esl1)+(1.-ilratio)*(rd/rv)*esi1/(presf(k)-(1.-rd/rv)*esi1)
            thlguessmin = ttry/exnf(k)-(rlv/(cp*exnf(k)))*max(qt200400(i,j)-qsatur,0.)

            Tnr = Tnr - (thlguess-thl200400(i,j))/((thlguess-thlguessmin)*500.)
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
        qlma(i,j,k) = max(qt200400(i,j)-qsatur,0.)
        thvma(i,j,k)=(thl200400(i,j)+(rlv*qlma(i,j,k))/(cp*exnf(k)))&
                     *(1+(rv/rd-1)*qt200400(i,j)-rv/rd*qlma(i,j,k)) ! calculate thv
        if(thvma(i,j,k)<thvf(k)-10) then
          matop(i,j)=k
        elseif(k==k1) then
          matop(i,j)=k
        endif
      endif
    enddo
    enddo
    enddo

    ! calculate top of moist adiabat for dcape calculations
    do k=1,k1
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

    ! dcape and CIN of Moist Adiabat
    ! no nice interpolation yet
    do j=2,j1
    do i=2,i1
    capemax(i,j)=0.
    cinmax(i,j)=0.
    dcape(i,j)=0.
    dcin(i,j)=0.
    dscin(i,j)=0.
    dcintot(i,j)=0.
    dscape(i,j)=0.
    enddo
    enddo

    do k=1,k1
    do j=2,j1
    do i=2,i1
      if(k<capetop(i,j)) then
        cinmax(i,j)=cinmax(i,j)+max(grav*dzf(k)*(thvf(k)-thvma(i,j,k))/thvf(k),0.)
        capemax(i,j)=capemax(i,j)+max(grav*dzf(k)*(thvma(i,j,k)-thvf(k))/thvf(k),0.)
      endif
     if(k<kdmax) then
        if(k<kcb) then
            dscape(i,j)=dscape(i,j)+max(grav*dzf(k)*(thvfull(i,j,k)-thvf(k))/thvf(k),0.)
            dscin(i,j)=dscin(i,j)+max(grav*dzf(k)*(thvf(k)-thvfull(i,j,k))/thvf(k),0.)
        else
            dcape(i,j)=dcape(i,j)+max(grav*dzf(k)*(thvfull(i,j,k)-thvf(k))/thvf(k),0.)
        endif
            dcintot(i,j)=dcintot(i,j)+max(grav*dzf(k)*(thvf(k)-thvfull(i,j,k))/thvf(k),0.)
        if((zf(k)<1.2*zf(kcb)).and.(zf(k)>zf(kcb))) then
            dcin(i,j)=dcin(i,j)+max(grav*dzf(k)*(thvf(k)-thvfull(i,j,k))/thvf(k),0.)
        endif
      endif
    enddo
    enddo
    enddo

    if (lnetcdf) then
      allocate(vars(1:imax,1:jmax,nvar))
      vars(:,:,1) = dcape(2:i1,2:j1)
      vars(:,:,2) = dscape(2:i1,2:j1)
      vars(:,:,3) = dcin(2:i1,2:j1)
      vars(:,:,4) = dscin(2:i1,2:j1)
      vars(:,:,5) = dcintot(2:i1,2:j1)
      vars(:,:,6) = capemax(2:i1,2:j1)
      vars(:,:,7) = cinmax(2:i1,2:j1)
      vars(:,:,8) = hw2cb(2:i1,2:j1)
      vars(:,:,9) = hw2max(2:i1,2:j1)
      vars(:,:,10) = qtcb(2:i1,2:j1)
      vars(:,:,11)= thlcb(2:i1,2:j1)
      vars(:,:,12) = wcb(2:i1,2:j1)
      vars(:,:,13) = buoycb(2:i1,2:j1)
      vars(:,:,14) = buoymax(2:i1,2:j1)
      vars(:,:,15) = qlcb(2:i1,2:j1)
      vars(:,:,16) = lwp(2:i1,2:j1)
      vars(:,:,17) = rwp(2:i1,2:j1)
      vars(:,:,18) = twp(2:i1,2:j1)
      vars(:,:,19) = cldtop(2:i1,2:j1)
      vars(:,:,20)= sprec(2:i1,2:j1)
      call writestat_nc(ncid4,1,tncname,(/rtimee/),nrec,.true.)
      call writestat_nc(ncid4,nvar,ncname(1:nvar,:),vars,nrec,imax,jmax)
      deallocate(vars)
    end if

    deallocate(dcape,dscape,dcin,dscin,dcintot,capemax,cinmax,hw2cb,hw2max,qtcb,thlcb,wcb,&
    buoycb,buoymax,qlcb,lwp,twp,rwp,cldtop,thvfull,thvma,qlma,capemask,capetop,matop,thl200400,qt200400,sprec)

  end subroutine docape

!> Clean up when leaving the run
  subroutine exitcape
    use modstat_nc, only : exitstat_nc,lnetcdf
    implicit none

    if(lcape .and. lnetcdf) then
    call exitstat_nc(ncid4)
    end if

  end subroutine exitcape

end module modcape
