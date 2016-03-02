!> \file modsamplingsurf.f90
!!  Calculates surface statistics under conditional criteria

!>
!!  Calculates statistics under conditional criteria
!! Currently implemented criteria for sampling are:
!! - Cloud (ql0>0)
!!
!!  \author Xabier Pedruzo, WUR
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
module modsamplingsurf

use modglobal, only : longint
use modsampsurfdata

implicit none
private
PUBLIC :: initsamplingsurf, samplingsurf, exitsamplingsurf
save
!NetCDF variables
  integer,parameter :: nvar = 32
  integer(kind=longint) :: idtav,itimeav,tnext,tnextwrite
  integer :: nsamples,isamp,isamptot
  character(20),dimension(10) :: samplname,longsamplname
 !tmlsm conditional variables---Xabi
 real, allocatable, dimension(:) :: Qnetavl,Havl,LEavl,G0avl,tendskinavl,rsavl, &
                                    raavl,tskinavl,cliqavl,wlavl,rssoilavl,rsvegavl,Respavl, &
                                    wco2avl,Anavl,gcco2avl,ciavl,co2surfavl,nrsampAGSl,tauavl,swdiravl,swdifavl

contains
!> Initialization routine, reads namelists and inits variables
  subroutine initsamplingsurf

    use modmpi,    only : comm3d, my_real,mpierr,myid,mpi_logical
    use modglobal, only : ladaptive, dtmax,k1,ifnamopt,fname_options,kmax,   &
                           dtav_glob,timeav_glob,btime,tres,cexpnr,ifoutput
    use modgenstat, only : idtav_prof=>idtav, itimeav_prof=>itimeav
    use modsurfdata, only: isurf,lrsAgs
    implicit none

    integer :: ierr

    namelist/NAMSAMPLINGsurf/ &
    dtav,timeav,lsampclsurf,lsampallsurf,lsampbuupsurf,lsampupsurf,lsampclearsurf,& 
    lsampclO10surf,lsampcl5_10surf,lsampcl2_5surf,lsampcl0_2surf
    

    dtav=dtav_glob;timeav=timeav_glob

    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMSAMPLINGsurf,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions NAMSAMPLINGsurf'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMSAMPLINGsurf'
      endif
      write(6 ,NAMSAMPLINGsurf)
      close(ifnamopt)
    end if


    call MPI_BCAST(timeav          ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(dtav            ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(lsampallsurf    ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampupsurf     ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampbuupsurf   ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampclsurf     ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampclearsurf  ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampcl0_2surf  ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampcl2_5surf  ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampcl5_10surf ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampclO10surf  ,1,MPI_LOGICAL,0,comm3d,mpierr)

    isamptot = 0
    if (lsampallsurf) then
      isamptot = isamptot + 1
      samplname (isamptot) = 'all'
      longsamplname(isamptot) = 'All '
    endif
    if (lsampupsurf) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'upd'
      longsamplname(isamptot) = 'Updraft '
    end if
    if (lsampbuupsurf) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'buup'
      longsamplname(isamptot) = 'Buoyant Updraft '
    end if
    if (lsampclsurf) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'cld'
      longsamplname(isamptot) = 'Cloud '
    end if
    if (lsampclearsurf) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'clear'
      longsamplname(isamptot) = 'Clear Sky'
    end if

    if (lsampcl0_2surf) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'cloud0_2'
      longsamplname(isamptot) = 'Cloud with 0<cctau<=2'
    end if
    
    if (lsampcl2_5surf) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'cloud2_5'
      longsamplname(isamptot) = 'Cloud with 2<cctau<=5'
    end if
    
    if (lsampcl5_10surf) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'cloud5_10'
      longsamplname(isamptot) = 'Cloud with 5<cctau<=10'
    end if

    if (lsampclO10surf) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'cloudO10'
      longsamplname(isamptot) = 'Cloud with cctau>10'
    end if

    if(isamptot < 2) return
    idtav = dtav/tres
    itimeav = timeav/tres

    tnext      = idtav   +btime
    tnextwrite = itimeav +btime

    if (abs(timeav/dtav-nint(timeav/dtav))>1e-4) then
      stop 'timeav must be a integer multiple of dtav'
    end if
    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'dtav should be a integer multiple of dtmax'
    end if
    if (isurf/=1) return
    if (lrsAgs) then
      allocate(nrsampAGSl (isamptot),Qnetavl    (isamptot),Havl      (isamptot),LEavl    (isamptot), &
               G0avl      (isamptot),tendskinavl(isamptot),rsavl     (isamptot),raavl    (isamptot), &
               tskinavl   (isamptot),cliqavl    (isamptot),wlavl     (isamptot),rssoilavl(isamptot), &
               rsvegavl   (isamptot),Respavl    (isamptot),wco2avl   (isamptot),Anavl    (isamptot), &
               gcco2avl   (isamptot),ciavl      (isamptot),co2surfavl(isamptot),tauavl   (isamptot), &
               swdiravl   (isamptot),swdifavl   (isamptot))  
  
 ! initialize variables
      
      nrsampAGSl  = 0.0
      Qnetavl     = 0.0
      Havl        = 0.0
      LEavl       = 0.0
      G0avl       = 0.0
      tendskinavl = 0.0
      rsavl       = 0.0
      raavl       = 0.0
      tskinavl    = 0.0
      cliqavl     = 0.0
      wlavl       = 0.0
      rssoilavl   = 0.0
      rsvegavl    = 0.0
      Respavl     = 0.0
      wco2avl     = 0.0
      Anavl       = 0.0
      gcco2avl    = 0.0
      ciavl       = 0.0
      co2surfavl  = 0.0
      tauavl      = 0.0
      swdiravl    = 0.0
      swdifavl    = 0.0
    

      if(myid==0)then
        do isamp = 1,isamptot
          open (ifoutput,file=trim(samplname(isamp))//'tmlsm.'//cexpnr,status='replace')
          write(ifoutput,'(//3A,/A,F5.0,A)') &
             '#-----------------------------',trim(samplname(isamp)),'tmlsm -------------------------------'      &
             ,'#',timeav,'--- AVERAGING TIMESTEP (s) --- '     

          write(ifoutput,'(4a)') &
             '#     time     nr_samples Qnet        H          LE        G0     ', &
             '   tendskin      rs            ra           tskin     cliq   ', &
             '    Wl           rssoil     rsveg       Resp      wco2         An     gc_CO2     ci      CO2_surf',& 
             '    cctau    swdir    swdif'
          write(ifoutput,'(4a)') &
             '#      [s]     []         [W/m2]     [W/m2]     [W/m2]    [W/m2]   ', &
             '    [W/m2]       [s/m]      [s/m]      [K]         [-]   ', &
             '     [m]          [s/m]      [s/m]                                   [mm/s?]               [ppm]',&   
             '    [-]      [w/m2]    [W/m2]'
          close (ifoutput)
        enddo
      endif
    endif !lrsAgs

  end subroutine initsamplingsurf
!> Cleans up after the run
  subroutine exitsamplingsurf
    use modmpi,     only : myid
    use modsurfdata, only: isurf,lrsAGS

    implicit none

    if (isamptot < 2) return
    if (isurf /= 1)   return
    if (lrsAgs) then

      deallocate( nrsampAGSl,Qnetavl,Havl,LEavl,G0avl,tendskinavl,rsavl,raavl, & 
                  tskinavl,cliqavl,wlavl,rssoilavl,rsvegavl,Respavl,wco2avl,Anavl, &
                  gcco2avl,ciavl,co2surfavl,tauavl,swdiravl,swdifavl)  
    end if

  end subroutine exitsamplingsurf

!> General routine, does the timekeeping
  subroutine samplingsurf
    use modglobal, only : rk3step,timee,dt_lim
    implicit none
   if (isamptot<2) return
    if (rk3step/=3) return
    if(timee<tnext .and. timee<tnextwrite) then
      dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))
      return
    end if
    if (timee>=tnext) then
      tnext = tnext+idtav
      do isamp = 1,isamptot
        call dosamplingsurf
      end do
    end if
    if (timee>=tnextwrite) then
      tnextwrite = tnextwrite+itimeav
      call writesamplingsurf
    end if
    dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))

  return
  end subroutine samplingsurf
!> Performs the actual sampling
  subroutine dosamplingsurf
    use modmpi,     only : myid
    use modglobal, only : i1,i2,j1,j2,kmax,k1,ih,jh,&
                          dx,dy,dzh,dzf,cp,rv,rlv,rd,ijtot, &
                          grav,om22,cu,nsv,zh
    use modsurfdata, only:Qnet,H,LE,G0,tendskin,rs,ra,tskin,cliq,wl,rssoil,rsveg, &
                          AnField,RespField,gcco2Field,ciField,indCO2,wco2Field,cctau
    use modfields,   only: w0,thl0,qt0,ql0,thv0h,exnf,svm
    use modraddata,  only:swdir,swdif
    use modmpi,    only : slabsum,my_real,mpi_integer,comm3d,mpierr,mpi_sum
    use modpois,   only : p
    use modsurfdata, only: isurf,lrsAgs

    implicit none

    logical, allocatable, dimension(:,:) :: maskAGS
    real, allocatable, dimension(:,:,:) :: w0f
    real, allocatable, dimension(:,:,:) :: thv0
    real, allocatable, dimension(:) :: thvav

    real, allocatable, dimension(:) :: thvhav

    integer :: i,j,k,km,kp,iih,iif

    if (isurf/=1) return
    if (lrsAgs) then
      allocate(thvav(k1))
      allocate(maskAGS (2-ih:i1+ih,2-jh:j1+jh))   
      allocate(thv0(2-ih:i1+ih,2-jh:j1+jh,k1),&
               w0f  (2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(thvhav(k1))
 
      !next thv0 and w0f calculations are already done in modsampling.AGS and
      !variables could be taken from there to speed up the code
      do k=1,k1
         thv0(2:i1,2:j1,k) = (thl0(2:i1,2:j1,k)+rlv*ql0(2:i1,2:j1,k)/(cp*exnf(k))) &
                           *(1+(rv/rd-1)*qt0(2:i1,2:j1,k)-rv/rd*ql0(2:i1,2:j1,k))
      enddo
      do k=1,kmax
        w0f (2:i1,2:j1,k) = 0.5*(w0 (2:i1,2:j1,k) + w0  (2:i1,2:j1,k+1))
      end do
 
      maskAGS  = .false.
 
      thvav = 0.0
      call slabsum(thvav,1,k1,thv0,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      thvav = thvav/ijtot
 
      thvhav = 0.0
      call slabsum(thvhav,1,k1,thv0h,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      thvhav = thvhav/ijtot
 
      select case (samplname(isamp))
      case ('upd')
        do i=2,i1
        do j=2,j1
        !we consider three positive first levels to be an updraft from surface
        iih=0
        iif=0
        do k=1,3
           if (w0f(i,j,k).gt.0.) then
               iif = iif+1
           elseif (w0(i,j,k).gt.0.) then
               iih = iih+1
          endif
          if (iih==3 .or. iif==3) then
            maskAGS(i,j) = .true.
          endif
        enddo  
        enddo
        enddo
     
      case ('buup')
        do i=2,i1
        do j=2,j1
        iih=0
        iif=0
        do k=1,3
          if (w0f(i,j,k).gt.0.0.and.thv0 (i,j,k) > thvav(k)) then
              iif = iif+1
          elseif (w0(i,j,k).gt.0.0.and.thv0h(i,j,k) > thvhav(k)) then
              iih = iih+1
        endif
          if (iih==3 .or. iif==3) then
            maskAGS(i,j) = .true.
          endif
        enddo
        enddo
        enddo
      
      case ('cld')
        do i=2,i1
        do j=2,j1
        do k=1,kmax
         if (ql0(i,j,k)>epsilon(1.0)) then
             maskAGS(i,j)= .true.
             exit
         end if    
        enddo
        enddo
        enddo
 
      case ('clear')
        maskAGS  = .true.
        do i=2,i1
        do j=2,j1
        do k=1,kmax
         if (ql0(i,j,k)>epsilon(1.0)) then
             maskAGS(i,j)= .false.
             exit
         end if    
        enddo
        enddo
        enddo
 
      case ('cloud0_2')
        do i=2,i1
        do j=2,j1
        if (cctau(i,j)>epsilon(1.0).and.cctau(i,j).le.2.0) then
             maskAGS(i,j)= .true.
        end if    
        enddo
        enddo
 
      case ('cloud2_5')
        do i=2,i1
        do j=2,j1
        if (cctau(i,j).gt.2.0.and.cctau(i,j).lt.5.0) then
             maskAGS(i,j)= .true.
        end if    
        enddo
        enddo
 
      case ('cloud5_10')
        do i=2,i1
        do j=2,j1
        if (cctau(i,j).ge.5.0.and.cctau(i,j).lt.10.0) then
             maskAGS(i,j)= .true.
        end if    
        enddo
        enddo
 
      case ('cloudO10')
        do i=2,i1
        do j=2,j1
        if (cctau(i,j).ge.10.0) then
             maskAGS(i,j)= .true.
        end if    
        enddo
        enddo
 
      case ('all')
          maskAGS = .true.
      end select



    !AGS-tmlsm values
      nrsampAGSl(isamp)         = nrsampAGSl (isamp)+count(maskAGS    (2:i1,2:j1))
 
      Qnetavl    (isamp) = Qnetavl    (isamp)+sum  (Qnet       (2:i1,2:j1),maskAGS(2:i1,2:j1))
      Havl       (isamp) = Havl       (isamp)+sum  (H          (2:i1,2:j1),maskAGS(2:i1,2:j1))
      LEavl      (isamp) = LEavl      (isamp)+sum  (LE         (2:i1,2:j1),maskAGS(2:i1,2:j1))
      G0avl      (isamp) = G0avl      (isamp)+sum  (G0         (2:i1,2:j1),maskAGS(2:i1,2:j1))
      tendskinavl(isamp) = tendskinavl(isamp)+sum  (tendskin   (2:i1,2:j1),maskAGS(2:i1,2:j1))
      rsavl      (isamp) = rsavl      (isamp)+sum  (rs         (2:i1,2:j1),maskAGS(2:i1,2:j1))
      raavl      (isamp) = raavl      (isamp)+sum  (ra         (2:i1,2:j1),maskAGS(2:i1,2:j1))
      tskinavl   (isamp) = tskinavl   (isamp)+sum  (tskin      (2:i1,2:j1),maskAGS(2:i1,2:j1))
      cliqavl    (isamp) = cliqavl    (isamp)+sum  (cliq       (2:i1,2:j1),maskAGS(2:i1,2:j1))
      wlavl      (isamp) = wlavl      (isamp)+sum  (wl         (2:i1,2:j1),maskAGS(2:i1,2:j1))
      rssoilavl  (isamp) = rssoilavl  (isamp)+sum  (rssoil     (2:i1,2:j1),maskAGS(2:i1,2:j1))
      rsvegavl   (isamp) = rsvegavl   (isamp)+sum  (rsveg      (2:i1,2:j1),maskAGS(2:i1,2:j1))
      Respavl    (isamp) = Respavl    (isamp)+sum  (RespField  (2:i1,2:j1),maskAGS(2:i1,2:j1))
      wco2avl    (isamp) = wco2avl    (isamp)+sum  (wco2Field  (2:i1,2:j1),maskAGS(2:i1,2:j1))
      Anavl      (isamp) = Anavl      (isamp)+sum  (AnField    (2:i1,2:j1),maskAGS(2:i1,2:j1))
      gcco2avl   (isamp) = gcco2avl   (isamp)+sum  (gcco2Field (2:i1,2:j1),maskAGS(2:i1,2:j1))
      ciavl      (isamp) = ciavl      (isamp)+sum  (ciField    (2:i1,2:j1),maskAGS(2:i1,2:j1))
      co2surfavl (isamp) = co2surfavl (isamp)+sum  (svm        (2:i1,2:j1,1,indCO2),maskAGS(2:i1,2:j1))
      tauavl     (isamp) = tauavl     (isamp)+sum  (cctau      (2:i1,2:j1),maskAGS(2:i1,2:j1))
      swdiravl   (isamp) = swdiravl   (isamp)+sum  (swdir      (2:i1,2:j1,1),maskAGS(2:i1,2:j1))
      swdifavl (isamp)   = swdifavl   (isamp)+sum  (swdif      (2:i1,2:j1,1),maskAGS(2:i1,2:j1))
 
      deallocate(maskAGS)
    end if !lrsAGs
  end subroutine dosamplingsurf
!> Write the statistics to file
  subroutine writesamplingsurf

    use modglobal, only : rtimee,k1,kmax,zf,zh,cexpnr,ifoutput,ijtot
    use modfields, only : presf,presh
    use modmpi,    only : myid,my_real,comm3d,mpierr,mpi_sum
    use modsurfdata, only:isurf,lrsAgs
    implicit none
    real,dimension(k1,nvar) :: vars

    !variables for AGS-tmlsm
    real                             :: nrsampAGSmn,Qnetmn,Hmn,LEmn,G0mn,tendskinmn,rsmn,ramn, &
                                        tskinmn,cliqmn,wlmn,rssoilmn,rsvegmn,Respmn,wco2mn,&
                                        Anmn,gcco2mn,cimn,co2surfmn,taumn,swdirmean,swdifmean
    real, allocatable, dimension(:)  :: nrsampAGSav,Qnetav,Hav,LEav,G0av,tendskinav,rsav,raav, &
                                        tskinav,cliqav,wlav,rssoilav,rsvegav,Respav,wco2av,&
                                        Anav,gcco2av,ciav,co2surfav,tauav,swdiraver,swdifaver
    integer :: nsecs, nhrs, nminut, k
    integer :: inorm
    if (lrsAgs) then
      allocate( nrsampAGSav(isamptot),Qnetav(isamptot),Hav(isamptot),LEav(isamptot),G0av(isamptot),&
                tendskinav(isamptot),rsav(isamptot),raav(isamptot),tskinav(isamptot), &
                cliqav(isamptot),wlav(isamptot),rssoilav(isamptot),rsvegav(isamptot), &
                Respav(isamptot),wco2av(isamptot),Anav(isamptot),gcco2av(isamptot), &
                ciav(isamptot),co2surfav(isamptot),tauav(isamptot),swdiraver(isamptot), &
                swdifaver(isamptot))
 
      nsecs   = nint(rtimee)
      nhrs    = int(nsecs/3600)
      nminut  = int(nsecs/60)-nhrs*60
      nsecs   = mod(nsecs,60)
      inorm   = nint(ijtot*timeav/dtav)
 
      call MPI_ALLREDUCE(nrsampAGSl  ,nrsampAGSav,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(Qnetavl     ,Qnetav     ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(Havl        ,Hav        ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(LEavl       ,LEav       ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(G0avl       ,G0av       ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(tendskinavl ,tendskinav ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(rsavl       ,rsav       ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(raavl       ,raav       ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(tskinavl    ,tskinav    ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(cliqavl     ,cliqav     ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(wlavl       ,wlav       ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(rssoilavl   ,rssoilav   ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(rsvegavl    ,rsvegav    ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(Respavl     ,Respav     ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(wco2avl     ,wco2av     ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(Anavl       ,Anav       ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(gcco2avl    ,gcco2av    ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(ciavl       ,ciav       ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(co2surfavl  ,co2surfav  ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(tauavl      ,tauav      ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(swdiravl    ,swdiraver  ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(swdifavl    ,swdifaver  ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      
    !reset variables
      nrsampAGSl  = 0.0
      Qnetavl     = 0.0
      Havl        = 0.0
      LEavl       = 0.0
      G0avl       = 0.0
      tendskinavl = 0.0
      rsavl       = 0.0
      raavl       = 0.0
      tskinavl    = 0.0
      cliqavl     = 0.0
      wlavl       = 0.0
      rssoilavl   = 0.0
      rsvegavl    = 0.0
      Respavl     = 0.0
      wco2avl     = 0.0
      Anavl       = 0.0
      gcco2avl    = 0.0
      ciavl       = 0.0
      co2surfavl  = 0.0
      tauavl      = 0.0
      swdiravl    = 0.0
      swdifavl    = 0.0
        
        
      if (myid==0) then
        do isamp = 1,isamptot
 
           nrsampAGSmn  = 0.0
           Qnetmn     = 0.0
           Hmn        = 0.0
           LEmn       = 0.0
           G0mn       = 0.0
           tendskinmn = 0.0
           rsmn       = 0.0
           ramn       = 0.0
           tskinmn    = 0.0
           cliqmn     = 0.0
           wlmn       = 0.0
           rssoilmn   = 0.0
           rsvegmn    = 0.0
           Respmn     = 0.0
           wco2mn     = 0.0
           Anmn       = 0.0
           gcco2mn    = 0.0
           cimn       = 0.0
           co2surfmn  = 0.0
           taumn      = 0.0
           swdirmean  = 0.0
           swdifmean  = 0.0
 
    !normalize variables
 
 
          if (nrsampAGSav(isamp).gt.0) then
 
            Qnetmn     = Qnetav    (isamp)/nrsampAGSav(isamp)
            Hmn        = Hav       (isamp)/nrsampAGSav(isamp)
            LEmn       = LEav      (isamp)/nrsampAGSav(isamp)
            G0mn       = G0av      (isamp)/nrsampAGSav(isamp)
            tendskinmn = tendskinav(isamp)/nrsampAGSav(isamp)
            rsmn       = rsav      (isamp)/nrsampAGSav(isamp)
            ramn       = raav      (isamp)/nrsampAGSav(isamp)
            tskinmn    = tskinav   (isamp)/nrsampAGSav(isamp)
            cliqmn     = cliqav    (isamp)/nrsampAGSav(isamp)
            wlmn       = wlav      (isamp)/nrsampAGSav(isamp)
            rssoilmn   = rssoilav  (isamp)/nrsampAGSav(isamp)
            rsvegmn    = rsvegav   (isamp)/nrsampAGSav(isamp)
            Respmn     = Respav    (isamp)/nrsampAGSav(isamp)
            wco2mn     = wco2av    (isamp)/nrsampAGSav(isamp)
            Anmn       = Anav      (isamp)/nrsampAGSav(isamp)
            gcco2mn    = gcco2av   (isamp)/nrsampAGSav(isamp)
            cimn       = ciav      (isamp)/nrsampAGSav(isamp)
            co2surfmn  = co2surfav (isamp)/nrsampAGSav(isamp)
            taumn      = tauav     (isamp)/nrsampAGSav(isamp)
            swdirmean  = swdiraver (isamp)/nrsampAGSav(isamp)
            swdifmean  = swdifaver (isamp)/nrsampAGSav(isamp)
          endif
          nrsampAGSmn= nrsampAGSav   (isamp)/inorm
 
 
    !write files
         write(ifoutput,'(f10.0, F10.4,6f11.3,f17.3,2f11.3,e13.3, 5f11.3,e13.3,2f9.2,f8.4,2f9.3)') &
          rtimee      , &
          nrsampAGSmn , &
          Qnetmn      , &
          Hmn         , &
          LEmn        , &
          G0mn        , &
          tendskinmn  , &
          rsmn        , &
          ramn        , &
          tskinmn     , &
          cliqmn      , &
          wlmn        , & 
          rssoilmn    , &  
          rsvegmn     , &
          Respmn      , &
          wco2mn      , &
          Anmn        , &
          gcco2mn     , &
          cimn        , &
          co2surfmn/1000   , &
          taumn       , &  
          swdirmean   , &
          swdifmean
         close(ifoutput)
                         
 
        end do
      end if
      deallocate( nrsampAGSav,Qnetav,Hav,LEav,G0av,tendskinav,rsav,raav,tskinav, &
                  cliqav,wlav,rssoilav,rsvegav,Respav,wco2av,Anav,gcco2av,ciav, &
                  co2surfav,tauav, swdiraver,swdifaver)
    end if !lrsAgs
  end subroutine writesamplingsurf

end module modsamplingsurf
