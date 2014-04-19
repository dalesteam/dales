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
! private
! public :: startup, writerestartfiles,trestart
save

  integer (KIND=selected_int_kind(6)) :: irandom= 0     !    * number to seed the randomnizer with
  integer :: krand = huge(0), krandumin=1,krandumax=0
  real :: randthl= 0.1,randqt=1e-5                 !    * thl and qt amplitude of randomnization
  real :: randu = 0.5

contains
  subroutine startup

      !-----------------------------------------------------------------|
      !                                                                 |
      !     Reads all general options from namoptions                   |
      !                                                                 |
      !      Chiel van Heerwaarden        15/06/2007                    |
      !      Thijs Heus                   15/06/2007                    |
      !-----------------------------------------------------------------|

    use modglobal,         only : initglobal,iexpnr,runtime, dtmax,dtav_glob,timeav_glob,&
                                  lwarmstart,startfile,trestart,itrestart,&
                                  nsv,imax,jtot,kmax,xsize,ysize,xlat,xlon,xday,xtime,&
                                  lmoist,lcoriol,igrw_damp,geodamptime,lmomsubs,cu, cv,ifnamopt,fname_options,llsadv,&
                                  ibas_prf,lambda_crit,iadv_mom,iadv_tke,iadv_thl,iadv_qt,iadv_sv,courant,peclet,ladaptive,author, lrigidlid, unudge
    use modforces,         only : lforce_user
    use modsurfdata,       only : z0,ustin,wtsurf,wqsurf,wsvsurf,ps,thls,isurf
    use modsurface,        only : initsurface
    use modfields,         only : initfields
    use modpois,           only : initpois
    use modradiation,      only : initradiation
    use modraddata,        only : irad,iradiation,&
                                  rad_ls,rad_longw,rad_shortw,rad_smoke,useMcICA,&
                                  timerad,itimerad,rka,dlwtop,dlwbot,sw0,gc,reff,isvsmoke
    use modtimedep,        only : inittimedep,ltimedep
    use modboundary,       only : initboundary,ksp
    use modthermodynamics, only : initthermodynamics,lqlnr, chi_half
    use modmicrophysics,   only : initmicrophysics
    use modsubgrid,        only : initsubgrid
    use modsubgriddata,    only : ldelta, cf,cn,Rigc,Prandtl,lmason,lsmagorinsky
    use modmpi,            only : comm3d,myid, mpi_integer,mpi_logical,my_real,mpierr, mpi_character

    implicit none
    integer :: ierr

    !declare namelists
    namelist/RUN/ &
        iexpnr,lwarmstart,startfile,runtime,dtmax,dtav_glob,timeav_glob,&
        trestart,irandom,randthl,randqt,krand,nsv,courant,peclet,ladaptive,author,&
        krandumin, krandumax, randu
    namelist/DOMAIN/ &
        imax,jtot,kmax,&
        xsize,ysize,&
        xlat,xlon,xday,xtime,ksp
    namelist/PHYSICS/ &
        !cstep z0,ustin,wtsurf,wqsurf,wsvsurf,ps,thls,chi_half,lmoist,isurf,lneutraldrag,&
        z0,ustin,wtsurf,wqsurf,wsvsurf,ps,thls,lmoist,isurf,chi_half,&
        lcoriol,igrw_damp,geodamptime,lmomsubs,ltimedep,irad,timerad,iradiation,rad_ls,rad_longw,rad_shortw,rad_smoke,useMcICA,&
        rka,dlwtop,dlwbot,sw0,gc,reff,isvsmoke,lforce_user,lrigidlid, unudge
    namelist/DYNAMICS/ &
        llsadv, lqlnr, lambda_crit, cu, cv, ibas_prf, iadv_mom, iadv_tke, iadv_thl, iadv_qt, iadv_sv

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
    end if
    
  !broadcast namelists
    call MPI_BCAST(iexpnr     ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(lwarmstart ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(startfile  ,50,MPI_CHARACTER,0,comm3d,mpierr)
    call MPI_BCAST(author     ,80,MPI_CHARACTER,0,comm3d,mpierr)
    call MPI_BCAST(runtime    ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(trestart   ,1,MY_REAL   ,0,comm3d,mpierr)
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
    call MPI_BCAST(lrigidlid  ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(unudge     ,1,MY_REAL   ,0,comm3d,mpierr)

    call MPI_BCAST(irad       ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(timerad    ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(iradiation ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(rad_ls     ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(rad_longw  ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(rad_shortw ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(rad_smoke  ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(useMcIca   ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(rka        ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(dlwtop     ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(dlwbot     ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(sw0        ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(gc         ,1,MY_REAL   ,0,comm3d,mpierr)
    ! CvH call MPI_BCAST(sfc_albedo ,1,MY_REAL   ,0,comm3d,mpierr)
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
    call MPI_BCAST(krandumin  ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(krandumax  ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(randthl    ,1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(randqt     ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(randu      ,1,MY_REAL   ,0,comm3d,mpierr)

    call MPI_BCAST(ladaptive  ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(courant,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(peclet,1,MY_REAL   ,0,comm3d,mpierr)

    call MPI_BCAST(isurf   ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(ibas_prf,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(iadv_mom,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(iadv_tke,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(iadv_thl,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(iadv_qt ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(iadv_sv(1:nsv) ,nsv,MPI_INTEGER,0,comm3d,mpierr)

    ! Allocate and initialize core modules
    call initglobal
    call initfields

    call initboundary
    call initthermodynamics
    call initradiation
    call initsurface
    call initsubgrid
    call initpois
    call initmicrophysics
    call inittimedep !depends on modglobal,modfields, modmpi, modsurf, modradiation

    call checkinitvalues
    call readinitfiles


  end subroutine startup


  subroutine checkinitvalues
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
  !                                                                 |
  !     *checkinitvalues* is called from *program*.                 |
  !                                                                 |
  !-----------------------------------------------------------------|

    use modsurfdata,only : wtsurf,wqsurf,ustin,thls,z0,isurf,ps,lhetero
    use modglobal, only : imax,jtot, ysize,xsize,dtmax,runtime, startfile,lwarmstart,eps1
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

  end subroutine checkinitvalues

  subroutine readinitfiles
    use modfields,         only : u0,v0,w0,um,vm,wm,thlm,thl0,thl0h,qtm,qt0,qt0h,&
                                  ql0,ql0h,thv0h,sv0,svm,e12m,e120,&
                                  rhobf,rhobh,drhobdzf,drhobdzh,&
                                  dudxls,dudyls,dvdxls,dvdyls,dthldxls,dthldyls,&
                                  dqtdxls,dqtdyls,dqtdtls,dpdxl,dpdyl,&
                                  wfls,whls,ug,vg,uprof,vprof,thlprof, qtprof,e12prof, svprof,&
                                  v0av,u0av,qt0av,ql0av,thl0av,sv0av,exnf,exnh,presf,presh,rhof,&
                                  thlpcar,thvh,thvf
    use modglobal,         only : i1,i2,ih,j1,j2,jh,kmax,k1,dtmax,idtmax,dt,rdt,runtime,timeleft,tres,&
                                  rtimee,timee,ntimee,ntrun,btime,dt_lim,nsv,&
                                  zf,zh,dzf,dzh,rv,rd,grav,cp,rlv,pref0,om23_gs,&
                                  rslabs,cu,cv,e12min,dzh,dtheta,dqt,dsv,cexpnr,ifinput,lwarmstart,itrestart,&
                                  trestart, ladaptive,llsadv,tnextrestart,ibas_prf
    use modsubgrid,        only : ekm,ekh
    use modsurfdata,       only : wtsurf,wqsurf,wsvsurf, &
                                  thls,tskin,tskinm,tsoil,tsoilm,phiw,phiwm,Wl,Wlm,thvs,ustin,ps,qts,isurf,svs,obl,oblav,&
                                  thvs_patch,lhetero,qskin 
    use modsurface,        only : surface,qtsurf,dthldz
    use modboundary,       only : boundary
    use modmpi,            only : slabsum,myid,comm3d,mpierr,my_real
    use modthermodynamics, only : thermodynamics,calc_halflev
    use moduser,           only : initsurf_user

    integer i,j,k,n

    real, allocatable :: height(:), th0av(:)
    real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1) :: thv0
    real tv


    character(80) chmess

    allocate (height(k1))
    allocate (th0av(k1))



    if (.not. lwarmstart) then

    !********************************************************************

    !    1.0 prepare initial fields from files 'prog.inp' and 'scalar.inp'
    !    ----------------------------------------------------------------

    !--------------------------------------------------------------------
    !    1.1 read fields
    !-----------------------------------------------------------------

      rdt = dtmax / 100.
      dt  = floor(rdt/tres)
      timee = 0
      if (myid==0) then
        open (ifinput,file='prof.inp.'//cexpnr)
        read (ifinput,'(a80)') chmess
        write(*,     '(a80)') chmess
        read (ifinput,'(a80)') chmess

        do k=1,kmax
          read (ifinput,*) &
                height (k), &
                thlprof(k), &
                qtprof (k), &
                uprof  (k), &
                vprof  (k), &
                e12prof(k)
        end do

        close(ifinput)
        write(*,*) 'height    thl      qt         u      v     e12'
        do k=kmax,1,-1
          write (*,'(f7.1,f8.1,e12.4,3f7.1)') &
                height (k), &
                thlprof(k), &
                qtprof (k), &
                uprof  (k), &
                vprof  (k), &
                e12prof(k)

        end do

        if (minval(e12prof(1:kmax)) < e12min) then
          write(*,*)  'e12 value is zero (or less) in prof.inp'
          do k=1,kmax
            e12prof(k) = max(e12prof(k),e12min)
          end do
        end if

      end if ! end if myid==0
    ! MPI broadcast numbers reading
      call MPI_BCAST(thlprof,kmax,MY_REAL   ,0,comm3d,mpierr)
      call MPI_BCAST(qtprof ,kmax,MY_REAL   ,0,comm3d,mpierr)
      call MPI_BCAST(uprof  ,kmax,MY_REAL   ,0,comm3d,mpierr)
      call MPI_BCAST(vprof  ,kmax,MY_REAL   ,0,comm3d,mpierr)
      call MPI_BCAST(e12prof,kmax,MY_REAL   ,0,comm3d,mpierr)
      do k=1,kmax
      do j=1,j2
      do i=1,i2
        thl0(i,j,k) = thlprof(k)
        thlm(i,j,k) = thlprof(k)
        qt0 (i,j,k) = qtprof (k)
        qtm (i,j,k) = qtprof (k)
        u0  (i,j,k) = uprof  (k) - cu
        um  (i,j,k) = uprof  (k) - cu
        v0  (i,j,k) = vprof  (k) - cv
        vm  (i,j,k) = vprof  (k) - cv
        w0  (i,j,k) = 0.0
        wm  (i,j,k) = 0.0
        e120(i,j,k) = e12prof(k)
        e12m(i,j,k) = e12prof(k)
        ekm (i,j,k) = 0.0
        ekh (i,j,k) = 0.0
      end do
      end do
      end do
    !---------------------------------------------------------------
    !  1.2 randomnize fields
    !---------------------------------------------------------------

      krand  = min(krand,kmax)
      do k = 1,krand
        call randomnize(qtm ,k,randqt ,irandom,ih,jh)
        call randomnize(qt0 ,k,randqt ,irandom,ih,jh)
        call randomnize(thlm,k,randthl,irandom,ih,jh)
        call randomnize(thl0,k,randthl,irandom,ih,jh)
      end do

      do k=krandumin,krandumax
        call randomnize(um  ,k,randu  ,irandom,ih,jh)
        call randomnize(u0  ,k,randu  ,irandom,ih,jh)
        call randomnize(vm  ,k,randu  ,irandom,ih,jh)
        call randomnize(v0  ,k,randu  ,irandom,ih,jh)
        call randomnize(wm  ,k,randu  ,irandom,ih,jh)
        call randomnize(w0  ,k,randu  ,irandom,ih,jh)
      end do

      svprof = 0.
      if(myid==0)then
        if (nsv>0) then
          open (ifinput,file='scalar.inp.'//cexpnr)
          read (ifinput,'(a80)') chmess
          read (ifinput,'(a80)') chmess
          do k=1,kmax
            read (ifinput,*) &
                  height (k), &
                  (svprof (k,n),n=1,nsv)
          end do
          open (ifinput,file='scalar.inp.'//cexpnr)
          write (6,*) 'height   sv(1) --------- sv(nsv) '
          do k=kmax,1,-1
            write (6,*) &
                  height (k), &
                (svprof (k,n),n=1,nsv)
          end do

        end if
      end if ! end if myid==0

      call MPI_BCAST(wsvsurf,nsv   ,MY_REAL   ,0,comm3d,mpierr)

      call MPI_BCAST(svprof ,k1*nsv,MY_REAL   ,0,comm3d,mpierr)
      do k=1,kmax
        do j=1,j2
          do i=1,i2
            do n=1,nsv
              sv0(i,j,k,n) = svprof(k,n)
              svm(i,j,k,n) = svprof(k,n)
            end do
          end do
        end do
      end do

!-----------------------------------------------------------------
!    2.2 Initialize surface layer and base profiles
!-----------------------------------------------------------------

      select case(isurf)
      case(1)
        tskin  = thls
        tskinm = tskin
        tsoilm = tsoil
        phiwm  = phiw
        Wlm    = Wl
      case(2)
        tskin  = thls
      case(3,4)
        thls = thlprof(1)
        qts  = qtprof(1)
        tskin  = thls
        qskin  = qts
      case(10)
        call initsurf_user
      end select

      ! Set initial Obukhov length to -0.1 for iteration
      obl   = -0.1
      oblav = -0.1

      call qtsurf
      
      dthldz = (thlprof(1) - thls)/zf(1)
      thvs = thls * (1. + (rv/rd - 1.) * qts)
      if(lhetero) thvs_patch = thvs  !Needed for initialization: thls_patch and qt_patch not yet calculated

      u0av(1)   = uprof(1)
      thl0av(1) = thlprof(1)
      svs = svprof(1,:)
     
      call baseprofs ! call baseprofs before thermodynamics
      call boundary
      call thermodynamics
      call surface

       ! Gradients at the top are now calculated in modboundary, every timestep
!      dtheta = (thlprof(kmax)-thlprof(kmax-1)) / dzh(kmax)
!      dqt    = (qtprof (kmax)-qtprof (kmax-1)) / dzh(kmax)
!      do n=1,nsv
!        dsv(n) = (svprof(kmax,n)-svprof(kmax-1,n)) / dzh(kmax)
!      end do

      call boundary
      call thermodynamics

    else !if lwarmstart

      call readrestartfiles
      um   = u0
      vm   = v0
      wm   = w0
      thlm = thl0
      qtm  = qt0
      svm  = sv0
      e12m = e120
      call calc_halflev
      exnf = (presf/pref0)**(rd/cp)
      exnh = (presh/pref0)**(rd/cp)

      do  j=2,j1
      do  i=2,i1
      do  k=2,k1
        thv0h(i,j,k) = (thl0h(i,j,k)+rlv*ql0h(i,j,k)/(cp*exnh(k))) &
                      *(1+(rv/rd-1)*qt0h(i,j,k)-rv/rd*ql0h(i,j,k))
      end do
      end do
      end do

      do  j=2,j1
      do  i=2,i1
      do  k=1,k1
        thv0(i,j,k) = (thl0(i,j,k)+rlv*ql0(i,j,k)/(cp*exnf(k))) &
                      *(1+(rv/rd-1)*qt0(i,j,k)-rv/rd*ql0(i,j,k))
      end do
      end do
      end do

      thvh=0.
      call slabsum(thvh,1,k1,thv0h,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1) ! redefine halflevel thv using calculated thv
      thvh = thvh/rslabs

      thvf = 0.0
      call slabsum(thvf,1,k1,thv0,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      thvf = thvf/rslabs

      u0av = 0.0
      v0av = 0.0
      thl0av = 0.0
      th0av  = 0.0
      qt0av  = 0.0
      ql0av  = 0.0
      sv0av = 0.

      call slabsum(u0av  ,1,k1,u0  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(v0av  ,1,k1,v0  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(thl0av,1,k1,thl0,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(qt0av ,1,k1,qt0 ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      call slabsum(ql0av ,1,k1,ql0 ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      do n=1,nsv
        call slabsum(sv0av(1,n),1,k1,sv0(1,1,1,n),2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      end do

      u0av  = u0av  /rslabs + cu
      v0av  = v0av  /rslabs + cv
      thl0av = thl0av/rslabs
      qt0av = qt0av /rslabs
      ql0av = ql0av /rslabs
      sv0av = sv0av /rslabs
      th0av  = thl0av + (rlv/cp)*ql0av/exnf
      thvh(1) = th0av(1)*(1+(rv/rd-1)*qt0av(1)-rv/rd*ql0av(1)) ! override first level
      
      do k=1,k1
        rhof(k) = presf(k)/(rd*thvf(k)*exnf(k))
      end do

      ! CvH - only do this for fixed timestepping. In adaptive dt comes from restartfile
      if(ladaptive .eqv. .false.) rdt=dtmax
      call baseprofs !call baseprofs

    end if

!-----------------------------------------------------------------
!    2.1 read and initialise fields
!-----------------------------------------------------------------


    if(myid==0)then
      open (ifinput,file='lscale.inp.'//cexpnr)
      read (ifinput,'(a80)') chmess
      read (ifinput,'(a80)') chmess
      write(6,*) ' height u_geo   v_geo    subs     ' &
                    ,'   dqtdx      dqtdy        dqtdtls     thl_rad '
      do  k=1,kmax
        read (ifinput,*) &
              height (k), &
              ug     (k), &
              vg     (k), &
              wfls   (k), &
              dqtdxls(k), &
              dqtdyls(k), &
              dqtdtls(k), &
              thlpcar(k)
      end do
      close(ifinput)

      do k=kmax,1,-1
        write (6,'(3f7.1,5e12.4)') &
              height (k), &
              ug     (k), &
              vg     (k), &
              wfls   (k), &
              dqtdxls(k), &
              dqtdyls(k), &
              dqtdtls(k), &
              thlpcar(k)
      end do


    end if ! end myid==0

! MPI broadcast variables read in

    call MPI_BCAST(ug       ,kmax,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(vg       ,kmax,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(wfls     ,kmax,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(dqtdxls  ,kmax,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(dqtdyls  ,kmax,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(dqtdtls  ,kmax,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(thlpcar  ,kmax,MY_REAL   ,0,comm3d,mpierr)

!-----------------------------------------------------------------
!    2.3 make large-scale horizontal pressure gradient
!-----------------------------------------------------------------

!******include rho if rho = rho(z) /= 1.0 ***********

    do k=1,kmax
      dpdxl(k) =  om23_gs*vg(k)
      dpdyl(k) = -om23_gs*ug(k)
    end do

  !-----------------------------------------------------------------
  !    2.5 make large-scale horizontal gradients
  !-----------------------------------------------------------------

    whls(1)  = 0.0
    do k=2,kmax
      whls(k) = ( wfls(k)*dzf(k-1) +  wfls(k-1)*dzf(k) )/(2*dzh(k))
    end do
    whls(k1) = (wfls(kmax)+0.5*dzf(kmax)*(wfls(kmax)-wfls(kmax-1)) &
                                                  /dzh(kmax))

  !******include rho if rho = rho(z) /= 1.0 ***********

    if (llsadv) then
      if (myid==0) stop 'llsadv should not be used anymore. Large scale gradients were calculated in a non physical way (and lmomsubs had to be set to true to retain conservation of mass)'
    end if
    dudxls   = 0.0
    dudyls   = 0.0
    dvdxls   = 0.0
    dvdyls   = 0.0
    dthldxls = 0.0
    dthldyls = 0.0

    idtmax = floor(dtmax/tres)
    btime   = timee
    timeleft=ceiling(runtime/tres)
    dt_lim = timeleft
    rdt = real(dt)*tres
    ntrun   = 0
    rtimee  = real(timee)*tres
    ntimee  = nint(timee/dtmax)
    itrestart = floor(trestart/tres)
    tnextrestart = btime + itrestart
    deallocate (height,th0av)


  end subroutine readinitfiles

  subroutine readrestartfiles

    use modsurfdata, only : ustar,thlflux,qtflux,svflux,dudz,dvdz,dthldz,dqtdz,ps,thls,qts,thvs,oblav,&
                           tsoil,phiw,tskin,Wl,isurf,ksoilmax,Qnet,swdavn,swuavn,lwdavn,lwuavn,nradtime,&
                           obl,xpatches,ypatches,ps_patch,thls_patch,qts_patch,thvs_patch,oblpatch,lhetero,qskin
    use modraddata, only: iradiation, useMcICA
    use modfields,  only : u0,v0,w0,thl0,qt0,ql0,ql0h,e120,dthvdz,presf,presh,sv0,tmp0,esl,qvsl,qvsi
    use modglobal,  only : i1,i2,ih,j1,j2,jh,k1,dtheta,dqt,dsv,startfile,timee,&
                          iexpnr,ntimee,tres,rk3step,ifinput,nsv,runtime,dt,rdt,cu,cv
    use modmpi,     only : cmyid, myid
    use modsubgriddata, only : ekm


    character(50) :: name
    integer i,j,k,n
    !********************************************************************

  !    1.0 Read initfiles
  !-----------------------------------------------------------------
    name = startfile
    name(5:5) = 'd'
    name(12:14)=cmyid
    write(6,*) 'loading ',name
    open(unit=ifinput,file=name,form='unformatted', status='old')

      read(ifinput)  (((u0    (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
!       u0 = u0-cu
      read(ifinput)  (((v0    (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
!       v0 = v0-cv
      read(ifinput)  (((w0    (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      read(ifinput)  (((thl0  (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      read(ifinput)  (((qt0   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      read(ifinput)  (((ql0   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      read(ifinput)  (((ql0h  (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      read(ifinput)  (((e120  (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      read(ifinput)  (((dthvdz(i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      read(ifinput)  (((ekm   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      read(ifinput)  (((tmp0   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      read(ifinput)  (((esl   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      read(ifinput)  (((qvsl   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      read(ifinput)  (((qvsi   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      read(ifinput)   ((ustar (i,j  ),i=1,i2      ),j=1,j2      )
      read(ifinput)   ((thlflux (i,j  ),i=1,i2      ),j=1,j2      )
      read(ifinput)   ((qtflux  (i,j  ),i=1,i2      ),j=1,j2      )
      read(ifinput)   ((dthldz(i,j  ),i=1,i2      ),j=1,j2      )
      read(ifinput)   ((dqtdz (i,j  ),i=1,i2      ),j=1,j2      )
      read(ifinput)  (  presf (    k)                            ,k=1,k1)
      read(ifinput)  (  presh (    k)                            ,k=1,k1)
      read(ifinput)  ps,thls,qts,thvs,oblav
      read(ifinput)  dtheta,dqt,timee,dt,tres
      read(ifinput)   ((obl (i,j  ),i=1,i2      ),j=1,j2      )
      read(ifinput)   ((tskin(i,j ),i=1,i2      ),j=1,j2      )
      read(ifinput)   ((qskin(i,j ),i=1,i2      ),j=1,j2      )
      if(lhetero) then
        read(ifinput)   ((ps_patch  (i,j),i=1,xpatches),j=1,ypatches)
        read(ifinput)   ((thls_patch(i,j),i=1,xpatches),j=1,ypatches)
        read(ifinput)   ((qts_patch (i,j),i=1,xpatches),j=1,ypatches)
        read(ifinput)   ((thvs_patch(i,j),i=1,xpatches),j=1,ypatches)
        read(ifinput)   ((oblpatch  (i,j),i=1,xpatches),j=1,ypatches)
      endif
    close(ifinput)

    if (nsv>0) then
      name(5:5) = 's'
      write(6,*) 'loading ',name
      open(unit=ifinput,file=name,form='unformatted')
      read(ifinput) ((((sv0(i,j,k,n),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1),n=1,nsv)
      read(ifinput) (((svflux(i,j,n),i=1,i2),j=1,j2),n=1,nsv)
      read(ifinput) (dsv(n),n=1,nsv)
      read(ifinput)  timee
      close(ifinput)
    end if

    if (isurf == 1) then
      name(5:5) = 'l'
      write(6,*) 'loading ',name
      open(unit=ifinput,file=name,form='unformatted')
      read(ifinput) (((tsoil(i,j,k),i=1,i2),j=1,j2),k=1,ksoilmax)
      read(ifinput) (((phiw(i,j,k),i=1,i2),j=1,j2),k=1,ksoilmax)
      read(ifinput) ((tskin(i,j),i=1,i2),j=1,j2)
      read(ifinput) ((Wl(i,j),i=1,i2),j=1,j2)
      read(ifinput) ((Qnet(i,j),i=1,i2),j=1,j2)
      if(iradiation == 1 .and. useMcICA) then
        read(ifinput) (((swdavn(i,j,n),i=1,i2),j=1,j2),n=1,nradtime)
        read(ifinput) (((swuavn(i,j,n),i=1,i2),j=1,j2),n=1,nradtime)
        read(ifinput) (((lwdavn(i,j,n),i=1,i2),j=1,j2),n=1,nradtime)
        read(ifinput) (((lwuavn(i,j,n),i=1,i2),j=1,j2),n=1,nradtime)
      end if
      read(ifinput)  timee
      close(ifinput)
    end if

  end subroutine readrestartfiles
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine writerestartfiles
    use modsurfdata,only: ustar,thlflux,qtflux,svflux,dudz,dvdz,dthldz,dqtdz,ps,thls,qts,thvs,oblav,&
                          tsoil,phiw,tskin,Wl,ksoilmax,isurf,ksoilmax,Qnet,swdavn,swuavn,lwdavn,lwuavn,nradtime,&
                          obl,xpatches,ypatches,ps_patch,thls_patch,qts_patch,thvs_patch,oblpatch,lhetero,qskin
    use modraddata, only: iradiation, useMcICA
    use modfields, only : u0,v0,w0,thl0,qt0,ql0,ql0h,e120,dthvdz,presf,presh,sv0,tmp0,esl,qvsl,qvsi
    use modglobal, only : i1,i2,ih,j1,j2,jh,k1,dsv,itrestart,tnextrestart,dt_lim,rtimee,timee,tres,cexpnr,&
                          ntimee,rtimee,rk3step,ifoutput,nsv,timeleft,dtheta,dqt,dt,cu,cv
    use modmpi,    only : cmyid,myid
    use modsubgriddata, only : ekm

    implicit none
    logical :: lexitnow = .false.
    integer imin,ihour
    integer i,j,k,n
    character(20) name,linkname

    if (timee == 0) return
    if (rk3step /=3) return
    name = 'exit_now.'//cexpnr
    inquire(file=trim(name), EXIST=lexitnow)

    if (timee<tnextrestart) dt_lim = min(dt_lim,tnextrestart-timee)
    if (timee>=tnextrestart .or. lexitnow) then
      tnextrestart = tnextrestart+itrestart
      ihour = floor(rtimee/3600)
      imin  = floor((rtimee-ihour * 3600) /3600. * 60.)
      name = 'initd  h  m   .'
      write (name(6:7)  ,'(i2.2)') ihour
      write (name(9:10) ,'(i2.2)') imin
      name(12:14)= cmyid
      name(16:18)= cexpnr
      open  (ifoutput,file=name,form='unformatted',status='replace')

      write(ifoutput)  (((u0 (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      write(ifoutput)  (((v0 (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      write(ifoutput)  (((w0    (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      write(ifoutput)  (((thl0  (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      write(ifoutput)  (((qt0   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      write(ifoutput)  (((ql0   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      write(ifoutput)  (((ql0h  (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      write(ifoutput)  (((e120  (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      write(ifoutput)  (((dthvdz(i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      write(ifoutput)  (((ekm   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      write(ifoutput)  (((tmp0   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      write(ifoutput)  (((esl   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      write(ifoutput)  (((qvsl   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      write(ifoutput)  (((qvsi   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      write(ifoutput)   ((ustar (i,j  ),i=1,i2      ),j=1,j2      )
      write(ifoutput)   ((thlflux (i,j  ),i=1,i2      ),j=1,j2      )
      write(ifoutput)   ((qtflux  (i,j  ),i=1,i2      ),j=1,j2      )
      write(ifoutput)   ((dthldz(i,j  ),i=1,i2      ),j=1,j2      )
      write(ifoutput)   ((dqtdz (i,j  ),i=1,i2      ),j=1,j2      )
      write(ifoutput)  (  presf (    k)                            ,k=1,k1)
      write(ifoutput)  (  presh (    k)                            ,k=1,k1)
      write(ifoutput)  ps,thls,qts,thvs,oblav
      write(ifoutput)  dtheta,dqt,timee,  dt,tres
      write(ifoutput)   ((obl (i,j  ),i=1,i2      ),j=1,j2      )
      write(ifoutput)   ((tskin(i,j ),i=1,i2      ),j=1,j2      )
      write(ifoutput)   ((qskin(i,j ),i=1,i2      ),j=1,j2      )
      if(lhetero) then
        write(ifoutput)  ((ps_patch  (i,j),i=1,xpatches),j=1,ypatches)
        write(ifoutput)  ((thls_patch(i,j),i=1,xpatches),j=1,ypatches)
        write(ifoutput)  ((qts_patch (i,j),i=1,xpatches),j=1,ypatches)
        write(ifoutput)  ((thvs_patch(i,j),i=1,xpatches),j=1,ypatches)
        write(ifoutput)  ((oblpatch  (i,j),i=1,xpatches),j=1,ypatches)
      endif

      close (ifoutput)
      linkname = name
      linkname(6:11) = "latest"
      call system("ln -sf "//name //" "//linkname)

      if (nsv>0) then
        name  = 'inits  h  m   .'
        write (name(6:7)  ,'(i2.2)') ihour
        write (name(9:10) ,'(i2.2)') imin
        name(12:14) = cmyid
        name(16:18) = cexpnr
        open  (ifoutput,file=name,form='unformatted')
        write(ifoutput) ((((sv0(i,j,k,n),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1),n=1,nsv)
        write(ifoutput) (((svflux(i,j,n),i=1,i2),j=1,j2),n=1,nsv)
        write(ifoutput) (dsv(n),n=1,nsv)
        write(ifoutput)  timee

        close (ifoutput)
        linkname = name
        linkname(6:11) = "latest"
        call system("ln -sf "//name //" "//linkname)

      end if

      if (isurf == 1) then
        name  = 'initl  h  m   .'
        write (name(6:7)  ,'(i2.2)') ihour
        write (name(9:10) ,'(i2.2)') imin
        name(12:14) = cmyid
        name(16:18) = cexpnr
        open  (ifoutput,file=name,form='unformatted')
        write(ifoutput) (((tsoil(i,j,k),i=1,i2),j=1,j2),k=1,ksoilmax)
        write(ifoutput) (((phiw(i,j,k),i=1,i2),j=1,j2),k=1,ksoilmax)
        write(ifoutput) ((tskin(i,j),i=1,i2),j=1,j2)
        write(ifoutput) ((Wl(i,j),i=1,i2),j=1,j2)
        write(ifoutput) ((Qnet(i,j),i=1,i2),j=1,j2)
        if(iradiation == 1 .and. useMcICA) then
          write(ifoutput) (((swdavn(i,j,n),i=1,i2),j=1,j2),n=1,nradtime)
          write(ifoutput) (((swuavn(i,j,n),i=1,i2),j=1,j2),n=1,nradtime)
          write(ifoutput) (((lwdavn(i,j,n),i=1,i2),j=1,j2),n=1,nradtime)
          write(ifoutput) (((lwuavn(i,j,n),i=1,i2),j=1,j2),n=1,nradtime)
        end if
        write(ifoutput)  timee

        close (ifoutput)
        linkname = name
        linkname(6:11) = "latest"
        call system("ln -sf "//name //" "//linkname)
      end if
      if (lexitnow) then
        timeleft = 0  !jump out of the time loop
      end if
      if (lexitnow .and. myid == 0 ) then
        open(1, file=trim(name), status='old')
        close(1,status='delete')
        write(*,*) 'Stopped at t=',rtimee
      end if

      if (myid==0) then
        write(*,'(A,F15.7,A,I4)') 'dump at time = ',rtimee,' unit = ',ifoutput
      end if

    end if


  end subroutine writerestartfiles

  subroutine exitmodules
    use modfields,         only : exitfields
    use modglobal,         only : exitglobal
    use modmpi,            only : exitmpi
    use modboundary,       only : exitboundary
    use modmicrophysics,   only : exitmicrophysics
    use modpois,           only : exitpois
    use modtimedep,        only : exittimedep
    use modradiation,      only : exitradiation
    use modsubgrid,        only : exitsubgrid
    use modsurface,        only : exitsurface
    use modthermodynamics, only : exitthermodynamics

    call exittimedep
    call exitthermodynamics
    call exitsurface
    call exitsubgrid
    call exitradiation
    call exitpois
    call exitmicrophysics
    call exitboundary
    call exitfields
    call exitglobal
    call exitmpi

 end subroutine exitmodules
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
    use modfields,         only : rhobf,rhobh,drhobdzf,drhobdzh
    use modglobal,         only : k1,kmax,zf,zh,dzf,dzh,rv,rd,grav,cp,pref0,lwarmstart,ibas_prf,cexpnr,ifinput,ifoutput
    use modsurfdata,       only : thls,ps,qts
    use modmpi,            only : myid,comm3d,mpierr,my_real
    implicit none

    real :: thvb,prsb ! for calculating moist adiabat
    integer :: j,k
    real, allocatable :: height(:),pb(:),tb(:)
    character(80) chmess
    real :: zsurf=0.
    real :: tsurf
    real :: psurf=101325.
    real,dimension(4) :: zmat=(/11000.,20000.,32000.,47000./)
    real,dimension(4) :: lapserate=(/-6.5/1000.,0.,1./1000,2.8/1000/) 
    real,dimension(4) :: pmat
    real,dimension(4) :: tmat

    allocate (height(k1),pb(k1),tb(k1))

    if(myid==0)then

      if((lwarmstart.eqv..true.).and.(ibas_prf /= 5)) then
        ibas_prf = 5
        print *, 'WARNING: warm start requires input files for density. ibas_prf defaulted to 5'
      endif

      if(ibas_prf==1) then !thv constant and hydrostatic balance
        thvb=thls*(1+(rv/rd-1)*qts) ! using thls, q_l assumed to be 0 during first time step
        do k=1,k1
          prsb=(ps**(rd/cp)-(grav*zf(k)*pref0**(rd/cp))/(cp*thvb))**(cp/rd) !As in thermodynamics
          rhobf(k)=prsb/(rd*thvb*((prsb/pref0)**(rd/cp)))
        enddo
        open (ifoutput,file='baseprof.inp.'//cexpnr)
        write(ifoutput,*) '#baseprofiles'
        write(ifoutput,*) '#height rhobf'
        do k=1,kmax
          write (ifoutput,'(1f7.1,E25.17)') &
                zf (k), &
                rhobf (k)
        enddo
        close(ifoutput)
      elseif(ibas_prf==2) then ! Quasi-Boussinesq (Similar to Dales 3, except for buoyancy term now depending on slab mean state)
        thvb=thls*(1+(rv/rd-1)*qts)
        rhobh(1)=ps/(rd*thvb*(ps/pref0)**(rd/cp))
        do k=1,k1
          rhobf(k)=rhobh(1)
        enddo
        open (ifoutput,file='baseprof.inp.'//cexpnr)
        write(ifoutput,*) '#baseprofiles'
        write(ifoutput,*) '#height rhobf'
        do k=1,kmax
          write (ifoutput,'(1f7.1,E25.17)') &
                zf (k), &
                rhobf (k)
        enddo
        close(ifoutput)
      elseif(ibas_prf==3) then! use standard atmospheric lapse rate with surface temperature offset
        tsurf=thls*(ps/pref0)**(rd/cp)
        pmat(1)=exp((log(ps)*lapserate(1)*rd+log(tsurf+zsurf*lapserate(1))*grav-&
          log(tsurf+zmat(1)*lapserate(1))*grav)/(lapserate(1)*rd))
        tmat(1)=tsurf+lapserate(1)*(zmat(1)-zsurf);
        ! write(*,*)(*,*) 'make profiles'
        
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
        open (ifoutput,file='baseprof.inp.'//cexpnr)
        write(ifoutput,*) '#baseprofiles'
        write(ifoutput,*) '#height rhobf'
        do k=1,kmax
          write (ifoutput,'(1f7.1,E25.17)') &
                zf (k), &
                rhobf (k)
        enddo
        close(ifoutput)
      elseif(ibas_prf==4) then! use standard atmospheric lapse rate without surface temperature offset
        tsurf=288.16
        pmat(1)=exp((log(ps)*lapserate(1)*rd+log(tsurf+zsurf*lapserate(1))*grav-&
          log(tsurf+zmat(1)*lapserate(1))*grav)/(lapserate(1)*rd))
        tmat(1)=tsurf+lapserate(1)*(zmat(1)-zsurf);
        ! write(*,*)(*,*) 'make profiles'
        
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
        open (ifoutput,file='baseprof.inp.'//cexpnr)
        write(ifoutput,*) '#baseprofiles'
        write(ifoutput,*) '#height rhobf'
        do k=1,kmax
          write (ifoutput,'(1f7.1,E25.17)') &
                zf (k), &
                rhobf (k)
        enddo
        close(ifoutput)
      end if

      ! Reset all background profiles
      rhobf=0.
      rhobh=0.

      ! Read background profiles in all cases
      open (ifinput,file='baseprof.inp.'//cexpnr)
      read (ifinput,'(a80)') chmess
      read (ifinput,'(a80)') chmess
      
      do k=1,kmax
        read (ifinput,*) &
                height(k), &
                rhobf (k)
      end do
      close(ifinput)

   
      rhobf(k1)=rhobf(kmax)+(zf(k1)-zf(kmax))/(zf(kmax)-zf(kmax-1))*(rhobf(kmax)-rhobf(kmax-1))

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
    write (6,*) ' height   rhobf       rhobh'
    do k=k1,1,-1
        write (6,'(1f7.1,2E25.17)') &
              height (k), &
              rhobf (k), &
              rhobh (k)
    end do
    write (6,*) ' height   drhobdzf    drhobdzh'
    do k=k1,1,-1
        write (6,'(1f7.1,2E25.17)') &
              height (k), &
              drhobdzf (k), &
              drhobdzh (k)
    end do

    end if ! ENDIF MYID=0

    ! MPI broadcast variables
    call MPI_BCAST(rhobf       ,k1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(rhobh       ,k1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(drhobdzf    ,k1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(drhobdzh    ,k1,MY_REAL   ,0,comm3d,mpierr)

    deallocate(pb,tb)

  end subroutine baseprofs

end module modstartup
