!> \file modstartup.f90
!!  Initializes the run

!>
!! Initializes the run.
!>
!! Modstartup reads the namelists and initial data, sets the fields and calls
!! the inits of the other routines where necessary. Reading and writing of the
!! restart files also live in this module.
!!  \author Chiel van Heerwaarden, Wageningen U.R.
!!  \author Thijs Heus,MPI-M
!!  \todo documentation
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

module modstartup


implicit none
 private
 public :: readnamoptions,checksettings,readinitfiles
save

  integer (KIND=selected_int_kind(6)) :: irandom= 0     !    * number to seed the randomnizer with
  integer :: krand = huge(0)
  real :: randthl= 0.1,randqt=1e-5                 !    * thl and qt amplitude of randomnization

contains
  subroutine readnamoptions

      !-----------------------------------------------------------------|
      !                                                                 |
      !     Reads all general options from namoptions                   |
      !                                                                 |
      !      Chiel van Heerwaarden        15/06/2007                    |
      !      Thijs Heus                   15/06/2007                    |
      !-----------------------------------------------------------------|

    use modglobal,         only : initglobal,iexpnr,runtime, dtmax,dtav_glob,timeav_glob,&
                                  lwarmstart,startfile,lncstartfile,trestart,lncrestart,itrestart,&
                                  nsv,imax,jtot,kmax,xsize,ysize,xlat,xlon,xday,xtime,&
                                  lmoist,lcoriol,igrw_damp,geodamptime,lmomsubs,cu, cv,ifnamopt,fname_options,llsadv,&
                                  ibas_prf,lambda_crit,iadv_mom,iadv_tke,iadv_thl,iadv_qt,iadv_sv,&
                                  courant,peclet,ladaptive,iTimeInt,author
    use modforces,         only : lforce_user
    use modsurfdata,       only : z0,ustin,wtsurf,wqsurf,wsvsurf,ps,thls,isurf,ksoilmax,mpatch, &
                                  tsoilav,tsoildeepav,phiwav,rootfav, &
                                  lmostlocal,lsmoothflux,lneutral,z0mav,z0hav, &
                                  rsisurf2,Cskinav,lambdaskinav,albedoav,Qnetav,cvegav,Wlav,&
                                  rsminav,rssoilminav,LAIav,gDav,lidealised, &
                                  lhetero,xpatches,ypatches,land_use,loldtable
    use modraddata,        only : irad,iradiation,&
                                  rad_ls,rad_longw,rad_shortw,rad_smoke,useMcICA,&
                                  timerad,itimerad,rka,dlwtop,dlwbot,sw0,gc,reff,isvsmoke,&
                                  lCnstZenith,cnstZenith,usero3,co2factor,ocean,lCnstAlbedo
    use modtimedep,        only : ltimedep
    use modboundary,       only : ksp
    use modthermodynamics, only : lqlnr, chi_half
    use modsubgriddata,    only : ldelta, cf,cn,Rigc,Prandtl,lmason,lsmagorinsky
    use modmpi,            only : comm3d,myid, mpi_integer,mpi_logical,my_real,mpierr, mpi_character

    implicit none
    integer :: ierr

    !declare namelists
    namelist/RUN/ &
        iexpnr,lwarmstart,startfile,lncstartfile,runtime,dtmax,dtav_glob,timeav_glob,&
        trestart,lncrestart,irandom,randthl,randqt,krand,nsv,courant,peclet,ladaptive,iTimeInt,author
    namelist/DOMAIN/ &
        imax,jtot,kmax,xsize,ysize,xlat,xlon,xday,xtime,ksp
    namelist/PHYSICS/ &
        !cstep z0,ustin,wtsurf,wqsurf,wsvsurf,ps,thls,chi_half,lmoist,isurf,lneutraldrag,&
        z0,ustin,wtsurf,wqsurf,wsvsurf,ps,thls,lmoist,isurf,chi_half,&
        lcoriol,igrw_damp,geodamptime,lmomsubs,ltimedep,irad,timerad,iradiation,rad_ls,rad_longw,rad_shortw,rad_smoke,useMcICA,&
        rka,dlwtop,dlwbot,sw0,gc,reff,isvsmoke,lCnstZenith,cnstZenith,lCnstAlbedo,usero3,co2factor,ocean,lforce_user
    namelist/DYNAMICS/ &
        llsadv, lqlnr, lambda_crit, cu, cv, ibas_prf, iadv_mom, iadv_tke, iadv_thl, iadv_qt, iadv_sv
    namelist/NAMSURFACE/ & !< Soil related variables
      isurf,tsoilav, tsoildeepav, phiwav, rootfav, &
      ! Land surface related variables
      lmostlocal, lsmoothflux, lneutral, z0mav, z0hav, rsisurf2, Cskinav, lambdaskinav, albedoav, Qnetav, cvegav, Wlav, &
      ! Jarvis-Steward related variables
      rsminav, rssoilminav, LAIav, gDav, &
      ! Prescribed values for isurf 2, 3, 4
      z0, thls, ps, ustin, wtsurf, wqsurf, wsvsurf,lidealised, &
      ! Heterogeneous variables
      lhetero, xpatches, ypatches, land_use, loldtable

    !read namelists
    if(myid==0)then
      if (command_argument_count() >=1) then
        call get_command_argument(1,fname_options)
      end if
      write (*,*) fname_options

      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      if (ierr /= 0) then
        stop 'ERROR:Namoptions does not exist'
      end if
      read (ifnamopt,RUN,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions RUN'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions RUN'
      endif
      write(6 ,RUN)
      rewind(ifnamopt)
      read (ifnamopt,DOMAIN,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions DOMAIN'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions DOMAIN'
      endif
      write(6 ,DOMAIN)
      rewind(ifnamopt)
      read (ifnamopt,PHYSICS,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions PHYSICS'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions PHYSICS'
      endif
      write(6 ,PHYSICS)
      rewind(ifnamopt)
      read (ifnamopt,DYNAMICS,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions DYNAMICS'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions DYNAMICS'
      endif
      write(6 ,DYNAMICS)
      close(ifnamopt)
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMSURFACE,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions NAMSURFACE'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMSURFACE'
      endif
      write(6 ,NAMSURFACE)
      close(ifnamopt)
    end if

    
  !broadcast namelists
    call MPI_BCAST(iexpnr     ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(lwarmstart ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(startfile  ,50,MPI_CHARACTER,0,comm3d,mpierr)
    call MPI_BCAST(lncstartfile,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(author     ,80,MPI_CHARACTER,0,comm3d,mpierr)
    call MPI_BCAST(runtime    ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(trestart   ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(lncrestart ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(dtmax      ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(dtav_glob  ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(timeav_glob,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(nsv        ,1,MPI_INTEGER,0,comm3d,mpierr)

    call MPI_BCAST(imax       ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(jtot       ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(kmax       ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(xsize      ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(ysize      ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(xlat       ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(xlon       ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(xday       ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(xtime      ,1,MY_REAL   ,0,comm3d,mpierr)

    call MPI_BCAST(z0         ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(ustin      ,1,MY_REAL   ,0,comm3d,mpierr)
    !call MPI_BCAST(lneutraldrag ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(wtsurf     ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(wqsurf     ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(wsvsurf(1:nsv),nsv,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(ps         ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(thls       ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(chi_half   ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(lmoist     ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lcoriol    ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(igrw_damp  ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(geodamptime,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(lforce_user,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lmomsubs   ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(ltimedep   ,1,MPI_LOGICAL,0,comm3d,mpierr)

    call MPI_BCAST(irad       ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(timerad    ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(iradiation ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(rad_ls     ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(rad_longw  ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(rad_shortw ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(rad_smoke  ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(useMcIca   ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lCnstZenith,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lCnstAlbedo,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(usero3     ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(co2factor  ,1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(ocean      ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(rka        ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(dlwtop     ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(dlwbot     ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(sw0        ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(cnstZenith ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(gc         ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(reff       ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(isvsmoke   ,1,MPI_INTEGER,0,comm3d,mpierr)

    call MPI_BCAST(llsadv     ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lqlnr      ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lambda_crit,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(cu         ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(cv         ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(ksp        ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(irandom    ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(krand      ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(randthl    ,1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(randqt     ,1,MY_REAL   ,0,comm3d,mpierr)

    call MPI_BCAST(ladaptive  ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(iTimeInt   ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(courant,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(peclet,1,MY_REAL   ,0,comm3d,mpierr)

    call MPI_BCAST(isurf   ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(ibas_prf,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(iadv_mom,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(iadv_tke,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(iadv_thl,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(iadv_qt ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(iadv_sv(1:nsv) ,nsv,MPI_INTEGER,0,comm3d,mpierr)
    ! NAMSURFACE
    call MPI_BCAST(isurf        , 1       , MPI_INTEGER, 0, comm3d, mpierr)
    call MPI_BCAST(tsoilav      , ksoilmax, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(tsoildeepav  , 1       , MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(phiwav       , ksoilmax, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(rootfav      , ksoilmax, MY_REAL, 0, comm3d, mpierr)

    call MPI_BCAST(lmostlocal   , 1, MPI_LOGICAL, 0, comm3d, mpierr)
    call MPI_BCAST(lsmoothflux  , 1, MPI_LOGICAL, 0, comm3d, mpierr)
    call MPI_BCAST(lneutral     , 1, MPI_LOGICAL, 0, comm3d, mpierr)
    call MPI_BCAST(z0mav        , 1, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(z0hav        , 1, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(rsisurf2     , 1, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(Cskinav      , 1, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(lambdaskinav , 1, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(albedoav     , 1, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(Qnetav       , 1, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(lidealised   , 1, MY_REAL, 0, comm3d, mpierr)

    call MPI_BCAST(rsminav      , 1, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(rssoilminav  , 1, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(cvegav       , 1, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(Wlav         , 1, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(LAIav        , 1, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(gDav         , 1, MY_REAL, 0, comm3d, mpierr)

    call MPI_BCAST(z0         ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(ustin      ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(wtsurf     ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(wqsurf     ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(wsvsurf(1:nsv),nsv,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(ps         ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(thls       ,1,MY_REAL   ,0,comm3d,mpierr)

    call MPI_BCAST(lhetero                    ,            1, MPI_LOGICAL, 0, comm3d, mpierr)
    call MPI_BCAST(loldtable                  ,            1, MPI_LOGICAL, 0, comm3d, mpierr)
    call MPI_BCAST(xpatches                   ,            1, MPI_INTEGER, 0, comm3d, mpierr)
    call MPI_BCAST(ypatches                   ,            1, MPI_INTEGER, 0, comm3d, mpierr)
    call MPI_BCAST(land_use(1:mpatch,1:mpatch),mpatch*mpatch, MPI_INTEGER, 0, comm3d, mpierr)


  end subroutine readnamoptions 


  subroutine checksettings
  !-----------------------------------------------------------------|
  !                                                                 |
  !      Thijs Heus   TU Delft  9/2/2006                            |
  !                                                                 |
  !     purpose.                                                    |
  !     --------                                                    |
  !                                                                 |
  !      checks whether crucial parameters are set correctly        |
  !                                                                 |
  !     interface.                                                  |
  !     ----------                                                  |
  !     called from program.f90                                     |
  !                                                                 |
  !-----------------------------------------------------------------|

    use modsurfdata,only : wtsurf,wqsurf,ustin,thls,z0,isurf,ps,lhetero
    use modglobal, only : imax,jtot, ysize,xsize,dtmax,runtime, startfile,lwarmstart,eps1,llsadv
    use modmpi,    only : myid, nprocs,mpierr
    use modtimedep, only : ltimedep

      if(mod(jtot,nprocs) /= 0) then
        if(myid==0)then
          write(6,*)'STOP ERROR IN NUMBER OF PROCESSORS'
          write(6,*)'nprocs must divide jtot!!! '
          write(6,*)'nprocs and jtot are: ',nprocs, jtot
        end if
        call MPI_FINALIZE(mpierr)
        stop
      end if

      if(mod(imax,nprocs)/=0)then
        if(myid==0)then
          write(6,*)'STOP ERROR IN NUMBER OF PROCESSORS'
          write(6,*)'nprocs must divide imax!!! '
          write(6,*)'nprocs and imax are: ',nprocs,imax
        end if
        call MPI_FINALIZE(mpierr)
        stop
      end if

  !Check Namroptions

    if (runtime < 0)stop 'runtime out of range/not set'
    if (dtmax < 0)  stop 'dtmax out of range/not set '
    if (ps < eps1)     stop 'psout of range/not set'
    if (thls < eps1)   stop 'thls out of range/not set'
    if (xsize < 0)  stop 'xsize out of range/not set'
    if (ysize < 0)  stop 'ysize out of range/not set '

    if (lwarmstart) then
      if (startfile == '') stop 'no restartfile set'
    end if
  !isurf
    if (myid == 0) then
      select case (isurf)
      case(1)
      case(2,10)
      case(3:4)
        if (wtsurf == -1)  stop 'wtsurf not set'
        if (wqsurf == -1)  stop 'wqsurf not set'
      case default
        stop 'isurf out of range/not set'
      end select
      if (isurf ==3) then
        if (ustin < 0)  stop 'ustin out of range/not set'
      end if
    end if

    if (ltimedep .and. lhetero) then
      if (myid == 0) write(6,*)&
      'WARNING: You selected to use time dependent (ltimedep) and heterogeneous surface conditions (lhetero) at the same time' 
      if (myid == 0) write(0,*)&
      'WARNING: You selected to use time dependent (ltimedep) and heterogeneous surface conditions (lhetero) at the same time' 
    endif

    if (llsadv) then
      if (myid==0) stop 'llsadv should not be used anymore. Large scale gradients were calculated in a non physical way (and lmomsubs had to be set to true to retain conservation of mass).'
    end if

  end subroutine checksettings

  subroutine readinitfiles
    use modfields,         only : u0,v0,w0,thl0,thl0h,qt0,qt0h,ql0,sv0,e120,ekm,ekh,&
                                  um,vm,wm,qtm,thlm,svm,e12m,&
                                  rhobf,rhobh,drhobdzf,drhobdzh,&
                                  dqtdxls,dqtdyls,dqtdtls,dpdxl,dpdyl,thlpcar,thvh,thvf,&
                                  wfls,whls,ug,vg,uprof,vprof,thlprof,qtprof,e12prof,svprof,&
                                  v0av,u0av,qt0av,ql0av,thl0av,sv0av,exnf,exnh,presf,presh,rhof,&
                                  dthldz,thvs_patch
    use modglobal,         only : i1,ih,j1,jh,kmax,k1,dtmax,idtmax,dt,rdt,runtime,timeleft,tres,&
                                  rtimee,timee,ntimee,ntrun,btime,dt_lim,subDt,nsv,&
                                  zf,zh,dzf,dzh,rv,rd,grav,cp,rlv,pref0,om23_gs,&
                                  rslabs,cu,cv,e12min,dzh,cexpnr,ifinput,lwarmstart,itrestart,&
                                  trestart,iTimeInt,iTimeWicker,iTimeTVD,ladaptive,tnextrestart,ibas_prf
    use modsurfdata,       only : isurf,thls,tskinm,tsoilm,tsoil,phiwm,phiw,&
                                  Wlm,Wl,qts,oblav,thvs,svs,lhetero
    use modboundary,       only : boundary
    use modmpi,            only : slabsum,myid,my_real
    use modthermodynamics, only : avgProfs,thermodynamics
    use modwarmstart,      only : readBinRestart
    use modtstep,          only : rkb
    use moduser,           only : initsurf_user
    use modsurface,        only : qtsurf,surface

    integer i,j,k,n

    ! Read profiles from ASCII input files
    call readInputProfs
     
    if (.not. lwarmstart) then

      ! Set initial times
      rdt = dtmax / 100.
      dt  = floor(rdt/tres)
      timee = 0

      ! Fill up the 3d fields with the input profiles 
      do j=2,j1; do i=2,i1
          thl0(i,j,1:kmax) = thlprof(1:kmax)
          qt0 (i,j,1:kmax) = qtprof (1:kmax)
          u0  (i,j,1:kmax) = uprof  (1:kmax) - cu
          v0  (i,j,1:kmax) = vprof  (1:kmax) - cv
          w0  (i,j,1:kmax) = 0.0
          e120(i,j,1:kmax) = e12prof(1:kmax)
          ekh (i,j,1:kmax) = 0.0
      end do; end do

      if (nsv>0) then
        do j=2,j1; do i=2,i1
          sv0(i,j,1:kmax,1:nsv) = svprof(1:kmax,1:nsv)
        end do; end do
      end if

    !---------------------------------------------------------------
    !  1.2 randomnize fields
    !---------------------------------------------------------------
      krand  = min(krand,kmax)
      do k = 1,krand
        call randomnize(qt0 ,k,randqt ,irandom,ih,jh)
        call randomnize(thl0,k,randthl,irandom,ih,jh)
      end do
      ! apply cyclic and top boundary conditions
      call boundary
      ! After setting the top boundary condition, the averaged profiles have to be calculated again,
      ! otherwise thl0av(k1)==0, which makes thermodynamics crash
      ekm  = ekh

      ! I moved the calculation of the profile averages to before modboundary.
      ! This is because that is more consistent with the determination of the top boundary condition
      ! However, during startup, this causes the k1 levels to be zero, which wreaks havoc in modthermodynamics
    else !if lwarmstart

      call readBinRestart

      if(ladaptive .eqv. .false.) rdt=dtmax
    end if

    ! Initialize m fields (indentical to what happens at a restart)
    ! Previously, the qtm and thlm fields were separately randomnized and therefore
    ! not identical to the qt0 and thl0 fields
    ! The m fields are not allocated for iTimeInt==2
    if (iTimeInt==iTimeWicker .or. iTimeInt==iTimeTVD) then
      thlm = thl0
      qtm  = qt0
      svm  = sv0
      e12m = e120
      um   = u0
      vm   = v0
      wm   = w0
    end if

    ! Calculate the horizontal averages (used in a.o. modthermodynamics)
    call avgProfs 
    ! Make the vertical profile of rho
    call baseprofs

    ! Set the timing variables
    idtmax = floor(dtmax/tres)
    btime   = timee
    timeleft=ceiling(runtime/tres)
    dt_lim = timeleft
    rdt = real(dt)*tres
    subDt = rkb(1)*rdt
    ntrun   = 0
    rtimee  = real(timee)*tres
    ntimee  = nint(timee/dtmax)
    itrestart = floor(trestart/tres)
    tnextrestart = btime + itrestart

  end subroutine readinitfiles

  !=====================================================================!
  ! This routine reads the basic input profiles from the ASCII files:   !
  ! prof.inp., lscale.inp. and scalar.inp (in case nsv>0).              !
  ! JvdD
  !=====================================================================!
  subroutine readInputProfs
    use modglobal, only : ifinput,cexpnr,e12min,kmax,k1,nsv
    use modmpi,    only : myid,my_real,comm3d,mpierr
    use modfields, only : thlprof,qtprof,uprof,vprof,e12prof,svprof,&
                          ug,vg,wfls,dqtdxls,dqtdyls,dqtdtls,thlpcar
    implicit none
    character(80) :: chmess
    real          :: height(k1)
    integer       :: k,n

    if (myid==0) then
      !==== Read prof.inp. ======!
      open (ifinput,file='prof.inp.'//cexpnr)
      read (ifinput,'(a80)') chmess
      write(   *   ,'(a80)') chmess
      read (ifinput,'(a80)') chmess

      ! First read the input profiles...
      do k=1,kmax
        read (ifinput,*) &
          height(k), thlprof(k), qtprof(k), uprof(k), vprof(k), e12prof(k)
      end do
      close(ifinput)
      ! ... then write them as a check to the output file
      write(*,*) 'height    thl      qt         u      v     e12'
      do k=kmax,1,-1
        write (*,'(f7.1,f8.1,e12.4,3f7.1)') &
          height(k), thlprof(k), qtprof(k), uprof(k), vprof(k), e12prof(k)
      end do

      !==== Read lscale.inp. ======!
      open (ifinput,file='lscale.inp.'//cexpnr)
      read (ifinput,'(a80)') chmess ! Header lines
      read (ifinput,'(a80)') chmess ! 

      ! First read the input profiles...
      do  k=1,kmax
        read (ifinput,*) &
          height(k),ug(k),vg(k),wfls(k),dqtdxls(k),dqtdyls(k),dqtdtls(k),thlpcar(k)
      end do
      close(ifinput)
      ! ... then write them as a check to the output file
      write(6,*) ' height u_geo   v_geo    subs     ' &
                    ,'   dqtdx      dqtdy        dqtdtls     thl_rad '
      do k=kmax,1,-1
        write (6,'(3f7.1,5e12.4)') &
          height(k),ug(k),vg(k),wfls(k),dqtdxls(k),dqtdyls(k),dqtdtls(k),thlpcar(k)
      end do

      !==== Read scalar.inp. ======!
      if (nsv>0) then
        open (ifinput,file='scalar.inp.'//cexpnr)
        read (ifinput,'(a80)') chmess ! Header lines
        read (ifinput,'(a80)') chmess !
        ! First read the input profiles...
        do k=1,kmax
          read (ifinput,*) &
            height(k), (svprof(k,n),n=1,nsv)
        end do
        ! ... then write them as a check to the output file
        open (ifinput,file='scalar.inp.'//cexpnr)
        write (6,*) 'height   sv(1) --------- sv(nsv) '
        do k=kmax,1,-1
          write (6,*) &
            height(k), (svprof(k,n),n=1,nsv)
        end do
      end if

      ! Check for unphysical input
      if (minval(e12prof(1:kmax)) < e12min) then
        write(*,*)  'e12 value is zero (or less) in prof.inp'
        write(*,*)  'Replaced by e12min.'
        do k=1,kmax
          e12prof(k) = max(e12prof(k),e12min)
        end do
      end if
    end if ! end if myid==0

    ! Broadcast the profiles to all processors
    !  prof.inp.
    call MPI_BCAST(thlprof,kmax  ,MY_REAL,0,comm3d,mpierr)
    call MPI_BCAST(qtprof ,kmax  ,MY_REAL,0,comm3d,mpierr)
    call MPI_BCAST(uprof  ,kmax  ,MY_REAL,0,comm3d,mpierr)
    call MPI_BCAST(vprof  ,kmax  ,MY_REAL,0,comm3d,mpierr)
    call MPI_BCAST(e12prof,kmax  ,MY_REAL,0,comm3d,mpierr)
    !  lscale.inp.
    call MPI_BCAST(ug     ,kmax  ,MY_REAL,0,comm3d,mpierr)
    call MPI_BCAST(vg     ,kmax  ,MY_REAL,0,comm3d,mpierr)
    call MPI_BCAST(wfls   ,kmax  ,MY_REAL,0,comm3d,mpierr)
    call MPI_BCAST(dqtdxls,kmax  ,MY_REAL,0,comm3d,mpierr)
    call MPI_BCAST(dqtdyls,kmax  ,MY_REAL,0,comm3d,mpierr)
    call MPI_BCAST(dqtdtls,kmax  ,MY_REAL,0,comm3d,mpierr)
    call MPI_BCAST(thlpcar,kmax  ,MY_REAL,0,comm3d,mpierr)
    !  scalar.inp
    call MPI_BCAST(svprof ,k1*nsv,MY_REAL,0,comm3d,mpierr)

  end subroutine readInputProfs

!----------------------------------------------------------------
  subroutine randomnize(field,klev,ampl,ir,ihl,jhl)

    use modmpi,    only :  myid,nprocs
    use modglobal, only : i1,imax,jmax,j1,k1
    integer (KIND=selected_int_kind(6)):: imm, ia, ic,ir
    integer ihl, jhl
    integer i,j,klev
    integer m,mfac
    real ran,ampl
    real field(2-ihl:i1+ihl,2-jhl:j1+jhl,k1)
    parameter (imm = 134456, ia = 8121, ic = 28411)

    if (myid>0) then
      mfac = myid * jmax * imax
      do m =1,mfac
        ir=mod((ir)*ia+ic,imm)

      end do
    end if
    do j=2,j1
    do i=2,i1
      ir=mod((ir)*ia+ic,imm)
      ran=real(ir)/real(imm)
      field(i,j,klev) = field(i,j,klev) + (ran-0.5)*2.0*ampl
    end do
    end do

    if (nprocs-1-myid > 0) then
      mfac = (nprocs-1-myid) * imax * jmax
      do m=1,mfac
        ir=mod((ir)*ia+ic,imm)
      end do
    end if

    return
  end subroutine randomnize

  subroutine baseprofs
    ! Calculates the profiles corresponding to the base state
    ! In the currecnt implementation, neither the base pressure, nor the base virtual temperature plays a role in the dynamics
    ! They are nevertheless calculated and printed to the stdin/baseprof files for user convenience
    use modfields,         only : rhobf,rhobh,drhobdzf,drhobdzh,rhodzi,thlprof,qtprof
    use modglobal,         only : k1,kmax,zf,zh,dzf,dzh,rv,rd,grav,cp,pref0,lwarmstart,ibas_prf,cexpnr,ifinput
    use modsurfdata,       only : ps
    use modmpi,            only : myid,comm3d,mpierr,my_real
    implicit none

    real :: thvb,prsb ! for calculating moist adiabat
    integer :: j,k
    real :: pb(k1),tb(k1)
    character(80) chmess
    real :: zsurf=0.
    real :: tsurf
    real :: psurf=101325.
    real,dimension(4) :: zmat=(/11000.,20000.,32000.,47000./)
    real,dimension(4) :: lapserate=(/-6.5/1000.,0.,1./1000,2.8/1000/) 
    real,dimension(4) :: pmat
    real,dimension(4) :: tmat

    ! For a warmstart, do not calculate a base profile as it is included in the restart file already. You do however need the halflevels and the gradients
    if (lwarmstart) ibas_prf=5

    ! NOTE: I replaced all references to thls and qts by thlprof(1) and qtprof(1)
    ! First the density profile was (slightly) dependent on the choice of the surface scheme
    ! Differences should be small.

    ! theta_v constant and hydrostatic balance
    if(ibas_prf==1) then
      thvb=thlprof(1)*(1+(rv/rd-1)*qtprof(1)) ! using thlprof(1), q_l assumed to be 0 during first time step
      do k=1,k1
        prsb=(ps**(rd/cp)-(grav*zf(k)*pref0**(rd/cp))/(cp*thvb))**(cp/rd) !As in thermodynamics
        rhobf(k)=prsb/(rd*thvb*((prsb/pref0)**(rd/cp)))
      enddo
    ! Quasi-Boussinesq (Similar to Dales 3, except for buoyancy term now depending on slab mean state)
    elseif(ibas_prf==2) then
      thvb=thlprof(1)*(1+(rv/rd-1)*qtprof(1))
      rhobh(1)=ps/(rd*thvb*(ps/pref0)**(rd/cp))
      do k=1,k1
        rhobf(k)=rhobh(1)
      enddo
    ! Use standard atmospheric lapse rate with surface temperature offset
    elseif(ibas_prf==3) then
      tsurf=thlprof(1)*(ps/pref0)**(rd/cp)
      pmat(1)=exp((log(ps)*lapserate(1)*rd+log(tsurf+zsurf*lapserate(1))*grav-&
        log(tsurf+zmat(1)*lapserate(1))*grav)/(lapserate(1)*rd))
      tmat(1)=tsurf+lapserate(1)*(zmat(1)-zsurf);
      
      do j=2,4
        if(abs(lapserate(j))<1e-10) then
          pmat(j)=exp((log(pmat(j-1))*tmat(j-1)*rd+zmat(j-1)*grav-zmat(j)*grav)/(tmat(j-1)*rd))
        else
          pmat(j)=exp((log(pmat(j-1))*lapserate(j)*rd+log(tmat(j-1)+zmat(j-1)*lapserate(j))*grav-&
            log(tmat(j-1)+zmat(j)*lapserate(j))*grav)/(lapserate(j)*rd))
        endif
        tmat(j)=tmat(j-1)+lapserate(j)*(zmat(j)-zmat(j-1));
      enddo
      
      do k=1,k1
        if(zf(k)<zmat(1)) then
          pb(k)=exp((log(ps)*lapserate(1)*rd+log(tsurf+zsurf*lapserate(1))*grav-&
            log(tsurf+zf(k)*lapserate(1))*grav)/(lapserate(1)*rd))
          tb(k)=tsurf+lapserate(1)*(zf(k)-zsurf)
        else
          j=1
          do while(zf(k)>zmat(j))
            j=j+1
          end do
          tb(k)=tmat(j-1)+lapserate(j)*(zf(k)-zmat(j-1))
          if(abs(lapserate(j))<1e-99) then
            pb(k)=exp((log(pmat(j-1))*tmat(j-1)*rd+zmat(j-1)*grav-zf(k)*grav)/(tmat(j-1)*rd))
          else
            pb(k)=exp((log(pmat(j-1))*lapserate(j)*rd+log(tmat(j-1)+zmat(j-1)*lapserate(j))*grav-&
              log(tmat(j-1)+zf(k)*lapserate(j))*grav)/(lapserate(j)*rd))
          endif
        endif
        rhobf(k)=pb(k)/(rd*tb(k)) ! dry estimate
      enddo
    ! Use standard atmospheric lapse rate without surface temperature offset
    elseif(ibas_prf==4) then
      tsurf=288.16
      pmat(1)=exp((log(ps)*lapserate(1)*rd+log(tsurf+zsurf*lapserate(1))*grav-&
        log(tsurf+zmat(1)*lapserate(1))*grav)/(lapserate(1)*rd))
      tmat(1)=tsurf+lapserate(1)*(zmat(1)-zsurf);
      
      do j=2,4
        if(abs(lapserate(j))<1e-10) then
          pmat(j)=exp((log(pmat(j-1))*tmat(j-1)*rd+zmat(j-1)*grav-zmat(j)*grav)/(tmat(j-1)*rd))
        else
          pmat(j)=exp((log(pmat(j-1))*lapserate(j)*rd+log(tmat(j-1)+zmat(j-1)*lapserate(j))*grav-&
            log(tmat(j-1)+zmat(j)*lapserate(j))*grav)/(lapserate(j)*rd))
        endif
        tmat(j)=tmat(j-1)+lapserate(j)*(zmat(j)-zmat(j-1));
      enddo
      
      do k=1,k1
        if(zf(k)<zmat(1)) then
          pb(k)=exp((log(ps)*lapserate(1)*rd+log(tsurf+zsurf*lapserate(1))*grav-&
            log(tsurf+zf(k)*lapserate(1))*grav)/(lapserate(1)*rd))
          tb(k)=tsurf+lapserate(1)*(zf(k)-zsurf)
        else
          j=1
          do while(zf(k)>zmat(j))
            j=j+1
          end do
          tb(k)=tmat(j-1)+lapserate(j)*(zf(k)-zmat(j-1))
          if(abs(lapserate(j))<1e-99) then
            pb(k)=exp((log(pmat(j-1))*tmat(j-1)*rd+zmat(j-1)*grav-zf(k)*grav)/(tmat(j-1)*rd))
          else
            pb(k)=exp((log(pmat(j-1))*lapserate(j)*rd+log(tmat(j-1)+zmat(j-1)*lapserate(j))*grav-&
              log(tmat(j-1)+zf(k)*lapserate(j))*grav)/(lapserate(j)*rd))
          endif
        endif
        rhobf(k)=pb(k)/(rd*tb(k)) ! dry estimate
      enddo
    end if

    ! Extrapolate the density profile to the top level
    rhobf(k1)=rhobf(kmax)+(zf(k1)-zf(kmax))/(zf(kmax)-zf(kmax-1))*(rhobf(kmax)-rhobf(kmax-1))

    ! Interpolate to find half level densities
    do k=2,k1
      rhobh(k) = (rhobf(k)*dzf(k-1)+rhobf(k-1)*dzf(k))/(dzf(k)+dzf(k-1))
    end do
    rhobh(1) = rhobf(1)-(rhobf(2)-rhobf(1))*(zf(1)-zh(1))/(zf(2)-zf(1))

    ! calculate derivatives
    do  k=1,kmax
      drhobdzf(k) = (rhobh(k+1) - rhobh(k))/dzf(k)
    end do
    drhobdzf(k1) = drhobdzf(kmax)
  
    drhobdzh(1) = 2*(rhobf(1)-rhobh(1))/dzh(1)
    do k=2,k1
      drhobdzh(k) = (rhobf(k)-rhobf(k-1))/dzh(k)
    end do

    ! write profiles and derivatives to standard output
    if (myid==0) then
      write (6,*) ' height   rhobf       rhobh'
      do k=k1,1,-1
          write (6,'(1f7.1,2E25.17)') &
                zf (k), &
                rhobf (k), &
                rhobh (k)
      end do
      write (6,*) ' height   drhobdzf    drhobdzh'
      do k=k1,1,-1
          write (6,'(1f7.1,2E25.17)') &
                zf (k), &
                drhobdzf (k), &
                drhobdzh (k)
      end do
    end if

    ! Useful for advection
    rhodzi = 1./(rhobf*dzf)

    ! MPI broadcast variables
    !call MPI_BCAST(rhobf       ,k1,MY_REAL   ,0,comm3d,mpierr)
    !call MPI_BCAST(rhobh       ,k1,MY_REAL   ,0,comm3d,mpierr)
    !call MPI_BCAST(drhobdzf    ,k1,MY_REAL   ,0,comm3d,mpierr)
    !call MPI_BCAST(drhobdzh    ,k1,MY_REAL   ,0,comm3d,mpierr)

  end subroutine baseprofs

end module modstartup
