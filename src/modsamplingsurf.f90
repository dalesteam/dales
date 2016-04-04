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
  character(20),dimension(25) :: samplname,longsamplname
 !tmlsm conditional variables---Xabi
 real, allocatable, dimension(:) :: Qnetavl,Havl,LEavl,G0avl,tendskinavl,rsavl, &
                                    raavl,tskinavl,cliqavl,wlavl,rssoilavl,rsvegavl,Respavl, &
                                    wco2avl,Anavl,gcco2avl,ciavl,co2surfavl,nrsampAGSl,tauavl,swdiravl,swdifavl,&
                                    An_stdevl,fAn_stdevl,LE_stdevl,fLE_stdevl

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
    lsampcl0o0_0o5surf,lsampcl0o5_1o0surf,lsampcl1o0_1o5surf,lsampcl1o5_2o0surf, &
    lsampcl2o0_2o5surf,lsampcl2o5_3o0surf,lsampcl3o0_3o5surf,lsampcl3o5_4o0surf, &
    lsampcl4o0_4o5surf,lsampcl4o5_5o0surf,lsampcl5_6surf,lsampcl6_7surf, & 
    lsampcl7_8surf,lsampcl8_9surf,lsampcl9_10surf,lsampcl10_15surf, &
    lsampcl15_20surf,lsampclO20surf 
    
    

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
    call MPI_BCAST(lsampcl0o0_0o5surf  ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampcl0o5_1o0surf  ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampcl1o0_1o5surf  ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampcl1o5_2o0surf  ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampcl2o0_2o5surf  ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampcl2o5_3o0surf  ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampcl3o0_3o5surf  ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampcl3o5_4o0surf  ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampcl4o0_4o5surf  ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampcl4o5_5o0surf  ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampcl5_6surf  ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampcl6_7surf  ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampcl7_8surf  ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampcl8_9surf  ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampcl9_10surf  ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampcl10_15surf  ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampcl15_20surf ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampclO20surf  ,1,MPI_LOGICAL,0,comm3d,mpierr)

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

    if (lsampcl0o0_0o5surf) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'cloud0.0_0.5'
      longsamplname(isamptot) = 'Cloud with 0.0<cctau<=0.5'
    end if

    if (lsampcl0o5_1o0surf) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'cloud0.5_1.0'
      longsamplname(isamptot) = 'Cloud with 0.5<cctau<11.0'
    end if

    if (lsampcl1o0_1o5surf) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'cloud1.0_1.5'
      longsamplname(isamptot) = 'Cloud with 1.0<cctau<1.5'
    end if

    if (lsampcl1o5_2o0surf) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'cloud1.5_2.0'
      longsamplname(isamptot) = 'Cloud with 1.5<cctau<2.0'
    end if

    if (lsampcl2o0_2o5surf) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'cloud2.0_2.5'
      longsamplname(isamptot) = 'Cloud with 2.0<cctau<2.5'
    end if

    if (lsampcl2o5_3o0surf) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'cloud2.5_3.0'
      longsamplname(isamptot) = 'Cloud with 2.5<cctau<3.0'
    end if

    if (lsampcl3o0_3o5surf) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'cloud3.0_3.5'
      longsamplname(isamptot) = 'Cloud with 3.0<cctau<3.5'
    end if

    if (lsampcl3o5_4o0surf) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'cloud3.5_4.0'
      longsamplname(isamptot) = 'Cloud with 3.5<cctau<4.0'
    end if

    if (lsampcl4o0_4o5surf) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'cloud4.0_4.5'
      longsamplname(isamptot) = 'Cloud with 4.0<cctau<4.5'
    end if
    
    if (lsampcl4o5_5o0surf) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'cloud4.5_5.0'
      longsamplname(isamptot) = 'Cloud with 4.5<cctau<=5.0'
    end if

    if (lsampcl5_6surf) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'cloud5_6'
      longsamplname(isamptot) = 'Cloud with 5<cctau<6'
    end if

    if (lsampcl6_7surf) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'cloud6_7'
      longsamplname(isamptot) = 'Cloud with 6<cctau<7'
    end if

    if (lsampcl7_8surf) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'cloud7_8'
      longsamplname(isamptot) = 'Cloud with 7<cctau<8'
    end if

    if (lsampcl8_9surf) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'cloud8_9'
      longsamplname(isamptot) = 'Cloud with 8<cctau<9'
    end if

    if (lsampcl9_10surf) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'cloud9_10'
      longsamplname(isamptot) = 'Cloud with 9<cctau<10'
    end if

    if (lsampcl10_15surf) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'cloud10_15'
      longsamplname(isamptot) = 'Cloud with 10<cctau<15'
    end if

    if (lsampcl15_20surf) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'cloud15_20'
      longsamplname(isamptot) = 'Cloud with 15<cctau<20'
    end if
    
    if (lsampclO20surf) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'cloudO20'
      longsamplname(isamptot) = 'Cloud with cctau>20'
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
    if (isurf/=1)return 
    if (lrsAgs) then
      allocate(nrsampAGSl (isamptot),Qnetavl    (isamptot),Havl      (isamptot),LEavl    (isamptot), &
               G0avl      (isamptot),tendskinavl(isamptot),rsavl     (isamptot),raavl    (isamptot), &
               tskinavl   (isamptot),cliqavl    (isamptot),wlavl     (isamptot),rssoilavl(isamptot), &
               rsvegavl   (isamptot),Respavl    (isamptot),wco2avl   (isamptot),Anavl    (isamptot), &
               gcco2avl   (isamptot),ciavl      (isamptot),co2surfavl(isamptot),tauavl   (isamptot), &
               swdiravl   (isamptot),swdifavl   (isamptot),An_stdevl (isamptot),fAn_stdevl(isamptot), &
               LE_stdevl  (isamptot),fLE_stdevl  (isamptot))  
  
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
      An_stdevl   = 0.0 
      fAn_stdevl  = 0.0
      LE_stdevl   = 0.0
      fLE_stdevl  = 0.0

      if(myid==0)then
        do isamp = 1,isamptot
          open (ifoutput,file=trim(samplname(isamp))//'tmlsm.'//cexpnr,status='replace')
          write(ifoutput,'(//3A,/A,F5.0,A)') &
             '#-----------------------------',trim(samplname(isamp)),'tmlsm -------------------------------'      &
             ,'#',timeav,'--- AVERAGING TIMESTEP (s) --- '     

          write(ifoutput,'(4a)') &
             '#     time     nr_samples Qnet        H          LE        G0     ', &
             '   tendskin      rs            ra           tskin     cliq   ', &
             '    Wl           rssoil     rsveg       Resp      wco2         An     gc_CO2        ci     CO2_surf',& 
             '  cctau    swdir    swdif    stdev_An   Fstdev_An       stdev_LE      Fstdev_LE '
          write(ifoutput,'(4a)') &
             '#      [s]     []         [W/m2]     [W/m2]     [W/m2]    [W/m2]   ', &
             '    [W/m2]       [s/m]      [s/m]      [K]         [-]   ', &
             '     [m]          [s/m]      [s/m]            [ppm m s-1]             [mm/s?]                 [ppm]',&   
             '    [-]      [w/m2]    [W/m2]  [accu mean]  [static mean(f)]    [accu mean,W/m2]  [static mean(f),W/m2] '
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
    if (isurf /= 1) return
    if (lrsAgs) then

      deallocate( nrsampAGSl,Qnetavl,Havl,LEavl,G0avl,tendskinavl,rsavl,raavl, & 
                  tskinavl,cliqavl,wlavl,rssoilavl,rsvegavl,Respavl,wco2avl,Anavl, &
                  gcco2avl,ciavl,co2surfavl,tauavl,swdiravl,swdifavl,An_stdevl, &
                  fAn_stdevl,LE_stdevl,fLE_stdevl)  
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
                          AnField,RespField,gcco2Field,ciField,indCO2,wco2Field,cctau,&
                          fAn_sqdiffl,An_sqdiffl,fLE_sqdiffl,LE_sqdiffl
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
    real    :: fAnavl,fLEavl, fnrsampAGSl,nrsampAGSaver
    real    :: fAnav,fLEav,fnrsampAGSaver,fAnaver,Anaver,fLEaver,LEaver

    if (isurf/=1) return
    if (lrsAgs) then
      allocate(thvav(k1))
      allocate(maskAGS (2-ih:i1+ih,2-jh:j1+jh))   
      allocate(thv0(2-ih:i1+ih,2-jh:j1+jh,k1),&
               w0f  (2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(thvhav(k1))
      allocate (fAn_sqdiffl(2:i1,2:j1),An_sqdiffl(2:i1,2:j1)) 
      allocate (fLE_sqdiffl(2:i1,2:j1),LE_sqdiffl(2:i1,2:j1)) 
      
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
      An_sqdiffl = 0.0
      fAn_sqdiffl = 0.0
      LE_sqdiffl = 0.0
      fLE_sqdiffl = 0.0
 
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
 
      case ('cloud0.0_0.5')
        do i=2,i1
        do j=2,j1
        if (cctau(i,j)>epsilon(1.0).and.cctau(i,j).le.0.5) then
             maskAGS(i,j)= .true.
        end if    
        enddo
        enddo
 
      case ('cloud0.5_1.0')
        do i=2,i1
        do j=2,j1
        if (cctau(i,j).gt.0.5 .and.cctau(i,j).le.1.0) then
             maskAGS(i,j)= .true.
        end if    
        enddo
        enddo
 
      case ('cloud1.0_1.5')
        do i=2,i1
        do j=2,j1
        if (cctau(i,j).gt.1.0 .and.cctau(i,j).le.1.5) then
             maskAGS(i,j)= .true.
        end if    
        enddo
        enddo
 
      case ('cloud1.5_2.0')
        do i=2,i1
        do j=2,j1
        if (cctau(i,j).gt.1.5 .and.cctau(i,j).le.2.0) then
             maskAGS(i,j)= .true.
        end if    
        enddo
        enddo
 
      case ('cloud2.0_2.5')
        do i=2,i1
        do j=2,j1
        if (cctau(i,j).gt.2.0 .and.cctau(i,j).le.2.5) then
             maskAGS(i,j)= .true.
        end if    
        enddo
        enddo
 
      case ('cloud2.5_3.0')
        do i=2,i1
        do j=2,j1
        if (cctau(i,j).gt.2.5 .and.cctau(i,j).le.3.0) then
             maskAGS(i,j)= .true.
        end if    
        enddo
        enddo
 
      case ('cloud3.0_3.5')
        do i=2,i1
        do j=2,j1
        if (cctau(i,j).gt.3.0 .and.cctau(i,j).le.3.5) then
             maskAGS(i,j)= .true.
        end if    
        enddo
        enddo
 
      case ('cloud3.5_4.0')
        do i=2,i1
        do j=2,j1
        if (cctau(i,j).gt.3.5 .and.cctau(i,j).le.4.0) then
             maskAGS(i,j)= .true.
        end if    
        enddo
        enddo
 
      case ('cloud4.0_4.5')
        do i=2,i1
        do j=2,j1
        if (cctau(i,j).gt.4.0 .and.cctau(i,j).le.4.5) then
             maskAGS(i,j)= .true.
        end if    
        enddo
        enddo
 
      case ('cloud4.5_5.0')
        do i=2,i1
        do j=2,j1
        if (cctau(i,j).gt.4.5 .and.cctau(i,j).le.5.0) then
             maskAGS(i,j)= .true.
        end if    
        enddo
        enddo
 
      case ('cloud5_6')
        do i=2,i1
        do j=2,j1
        if (cctau(i,j).gt.5.0.and.cctau(i,j).le.6.0) then
             maskAGS(i,j)= .true.
        end if    
        enddo
        enddo
 
      case ('cloud6_7')
        do i=2,i1
        do j=2,j1
        if (cctau(i,j).gt.6.0.and.cctau(i,j).le.7.0) then
             maskAGS(i,j)= .true.
        end if    
        enddo
        enddo
 
      case ('cloud7_8')
        do i=2,i1
        do j=2,j1
        if (cctau(i,j).gt.7.0.and.cctau(i,j).le.8.0) then
             maskAGS(i,j)= .true.
        end if    
        enddo
        enddo
 
      case ('cloud8_9')
        do i=2,i1
        do j=2,j1
        if (cctau(i,j).gt.8.0.and.cctau(i,j).le.9.0) then
             maskAGS(i,j)= .true.
        end if    
        enddo
        enddo
 
      case ('cloud9_10')
        do i=2,i1
        do j=2,j1
        if (cctau(i,j).gt.9.0.and.cctau(i,j).le.10.0) then
             maskAGS(i,j)= .true.
        end if    
        enddo
        enddo
 
      case ('cloud10_15')
        do i=2,i1
        do j=2,j1
        if (cctau(i,j).gt.10.0.and.cctau(i,j).le.15.0) then
             maskAGS(i,j)= .true.
        end if    
        enddo
        enddo
 
      case ('cloud15_20')
        do i=2,i1
        do j=2,j1
        if (cctau(i,j).gt.15.0.and.cctau(i,j).le.20.0) then
             maskAGS(i,j)= .true.
        end if    
        enddo
        enddo
 
      case ('cloudO20')
        do i=2,i1
        do j=2,j1
        if (cctau(i,j).gt.20.0) then
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
      swdifavl   (isamp) = swdifavl   (isamp)+sum  (swdif      (2:i1,2:j1,1),maskAGS(2:i1,2:j1))

!for standard deviation we first calculate the mean
!!!NOTE: calculations of stdev are perfectly reliable if dtav = timeav 
!else, results are not that statistically strong because the mean used for calculation
! of stdev is dynamic and used previous loops to  update the average for every dosampling loop
      !!!!!
      !!!option 1
      !fixed averaging: we calculate the stdev with respect to the mean  at that
      !timestep
      !!!!!
      !vars to be defined:fAnavl,fnrsampAGSl real     nrsampAGSaver 
      !!                  :fAnav,fnrsampAGSaver,fAnaver,   real  Anaver,nrsampAGSaver
      
      !                   :fAn_sqdiffl, (like Anfield)          An_sqdiffl
      !                   :fAn_stdevl, real isamp dimensions     An_stdevl
      
      fAnavl       = sum  (AnField    (2:i1,2:j1),maskAGS(2:i1,2:j1))
      fLEavl       = sum  (LE         (2:i1,2:j1),maskAGS(2:i1,2:j1))
      fnrsampAGSl  = count(maskAGS    (2:i1,2:j1))
      call MPI_ALLREDUCE(fAnavl          ,fAnav      ,1,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(fLEavl          ,fLEav      ,1,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(fnrsampAGSl     ,fnrsampAGSaver   ,1,MY_REAL,MPI_SUM,comm3d,mpierr)
      
      if (fnrsampAGSaver .gt. 0.0) then
        fAnaver = fAnav/fnrsampAGSaver
        fLEaver = fLEav/fnrsampAGSaver
        fAn_sqdiffl(2:i1,2:j1) = (AnField (2:i1,2:j1) - fAnaver) * (AnField(2:i1,2:j1) - fAnaver)
        fLE_sqdiffl(2:i1,2:j1) = (LE      (2:i1,2:j1) - fLEaver) * (LE     (2:i1,2:j1) - fLEaver)
        fAn_stdevl (isamp) = fAn_stdevl(isamp) + sum(fAn_sqdiffl(2:i1,2:j1),maskAGS(2:i1,2:j1)) 
        fLE_stdevl (isamp) = fLE_stdevl(isamp) + sum(fLE_sqdiffl(2:i1,2:j1),maskAGS(2:i1,2:j1)) 
      end if
!!!!!!!!!!
!!!!!!!!!Option2:if we do dosampling more times than write sampling, every time
!we do dosampling and not writesampling we keep the average values and use them
!for next dosampling stdev calculatuion. --this means there will also be a time
!average, then
!!!!!!!!!!
      !progressive averaging
      call MPI_ALLREDUCE(Anavl       ,Anaver         ,1,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(LEavl       ,LEaver         ,1,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(nrsampAGSl  ,nrsampAGSaver  ,1,MY_REAL,MPI_SUM,comm3d,mpierr)

!the global domain average is:
      if (nrsampAGSaver .gt. 0.0) then
        Anaver= Anaver / nrsampAGSaver
        LEaver= LEaver / nrsampAGSaver
      
      !now we calculate the local square differneces:
      
        An_sqdiffl (2:i1,2:j1)= (AnField (2:i1,2:j1) - Anaver) * (AnField (2:i1,2:j1) - Anaver)
        LE_sqdiffl (2:i1,2:j1)= (LE      (2:i1,2:j1) - LEaver) * (LE (2:i1,2:j1) - LEaver)
      
      !now that we have a field of squared differences, we calculate the local sum of it
      
        An_stdevl (isamp) = An_stdevl(isamp) + sum(An_sqdiffl (2:i1,2:j1),maskAGS(2:i1,2:j1))
        LE_stdevl (isamp) = LE_stdevl(isamp) + sum(LE_sqdiffl (2:i1,2:j1),maskAGS(2:i1,2:j1))
      end if
      

      deallocate(maskAGS)
      deallocate (fAn_sqdiffl,An_sqdiffl)
      deallocate (fLE_sqdiffl,LE_sqdiffl)
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
                                        Anmn,gcco2mn,cimn,co2surfmn,taumn,swdirmean,swdifmean,&
                                        An_stdevmn,fAn_stdevmn,LE_stdevmn,fLE_stdevmn
    real, allocatable, dimension(:)  :: nrsampAGSav,Qnetav,Hav,LEav,G0av,tendskinav,rsav,raav, &
                                        tskinav,cliqav,wlav,rssoilav,rsvegav,Respav,wco2av,&
                                        Anav,gcco2av,ciav,co2surfav,tauav,swdiraver,swdifaver,&
                                        An_stdev,fAn_stdev,LE_stdev,fLE_stdev
    integer :: nsecs, nhrs, nminut, k
    integer :: inorm
    if (lrsAgs) then
      allocate( nrsampAGSav(isamptot),Qnetav(isamptot),Hav(isamptot),LEav(isamptot),G0av(isamptot),&
                tendskinav(isamptot),rsav(isamptot),raav(isamptot),tskinav(isamptot), &
                cliqav(isamptot),wlav(isamptot),rssoilav(isamptot),rsvegav(isamptot), &
                Respav(isamptot),wco2av(isamptot),Anav(isamptot),gcco2av(isamptot), &
                ciav(isamptot),co2surfav(isamptot),tauav(isamptot),swdiraver(isamptot), &
                swdifaver(isamptot),An_stdev(isamptot),fAn_stdev(isamptot),LE_stdev(isamptot), &
                fLE_stdev(isamptot))
 
      nsecs   = nint(rtimee)
      nhrs    = int(nsecs/3600)
      nminut  = int(nsecs/60)-nhrs*60
      nsecs   = mod(nsecs,60)
      inorm   = nint(ijtot*timeav/dtav)
 
      call MPI_ALLREDUCE(nrsampAGSl   ,nrsampAGSav,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(Qnetavl      ,Qnetav     ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(Havl         ,Hav        ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(LEavl        ,LEav       ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(G0avl        ,G0av       ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(tendskinavl  ,tendskinav ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(rsavl        ,rsav       ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(raavl        ,raav       ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(tskinavl     ,tskinav    ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(cliqavl      ,cliqav     ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(wlavl        ,wlav       ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(rssoilavl    ,rssoilav   ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(rsvegavl     ,rsvegav    ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(Respavl      ,Respav     ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(wco2avl      ,wco2av     ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(Anavl        ,Anav       ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(gcco2avl     ,gcco2av    ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(ciavl        ,ciav       ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(co2surfavl   ,co2surfav  ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(tauavl       ,tauav      ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(swdiravl     ,swdiraver  ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(swdifavl     ,swdifaver  ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(An_stdevl    ,An_stdev   ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(fAn_stdevl   ,fAn_stdev  ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(LE_stdevl    ,LE_stdev   ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(fLE_stdevl   ,fLE_stdev  ,isamptot,MY_REAL,MPI_SUM,comm3d,mpierr)
      
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
      An_stdevl   = 0.0
      fAn_stdevl  = 0.0
      LE_stdevl   = 0.0
      fLE_stdevl  = 0.0
        
        
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
           An_stdevmn  = 0.0
           fAn_stdevmn = 0.0
           LE_stdevmn  = 0.0
           fLE_stdevmn = 0.0
  
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
             An_stdevmn  = sqrt(An_stdev (isamp)/nrsampAGSav(isamp))
             fAn_stdevmn = sqrt(fAn_stdev (isamp)/nrsampAGSav(isamp))
             LE_stdevmn  = sqrt(LE_stdev (isamp)/nrsampAGSav(isamp))
             fLE_stdevmn = sqrt(fLE_stdev (isamp)/nrsampAGSav(isamp))
           endif
           nrsampAGSmn= nrsampAGSav   (isamp)/inorm
  
  
    !write files
 !     if (myid==0) then
           open(ifoutput,file=trim(samplname(isamp))//'tmlsm.'//cexpnr,position='append')
           write(ifoutput,'(f10.0, F11.8,6f11.3,f17.3,2f11.3,e13.3, 5f11.3,e13.3,2f9.2,f8.4,2f9.3,2f11.5,2f8.3)') &
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
           swdifmean   , &
           An_stdevmn  , &
           fAn_stdevmn , &
           LE_stdevmn  , &
           fLE_stdevmn
           close(ifoutput)
                         
 
        !end if
        end do
      end if
      deallocate( nrsampAGSav,Qnetav,Hav,LEav,G0av,tendskinav,rsav,raav,tskinav, &
                  cliqav,wlav,rssoilav,rsvegav,Respav,wco2av,Anav,gcco2av,ciav, &
                  co2surfav,tauav,swdiraver,swdifaver,An_stdev,fAn_stdev, &
                  LE_stdev,fLE_stdev)
    end if !lrsAgs
  end subroutine writesamplingsurf

end module modsamplingsurf
