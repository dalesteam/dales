!> \file modstartup.f90
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
!  Copyright 1993-2026 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!

module modstartup
use iso_c_binding
use fortran_support, only: int2string, split_string
use modprecision,      only : field_r
use modtimer,        only: timer_init, timer_tic, timer_toc
use modstat_nc
use modchecksim, only: check_array
use modlogging, only: warning, finish, message

implicit none
! private
! public :: startup, writerestartfiles,trestart
  character(len=*), parameter :: modname = "modstartup"
save

  integer (KIND=selected_int_kind(6)) :: irandom= 0     !    * number to seed the randomnizer with
  integer :: krand = huge(0), krandumin=1,krandumax=0
  real :: randthl= 0.1,randqt=1e-5                 !    * thl and qt amplitude of randomnization
  real :: randu = 0.5
  real :: wctime=8640000.   !<     * The maximum wall clock time of a simulation (set to 100 days by default)

interface ! interface to use UNIX C mkdir function. Otherwise different compilers have different incompatible variants
   function mkdir(path,mode) bind(c,name="mkdir")
     use iso_c_binding
     integer(c_int) :: mkdir
     character(kind=c_char,len=1) :: path(*)
     integer(c_int16_t), value :: mode
   end function mkdir
end interface

contains
  subroutine startup(path)

      !-----------------------------------------------------------------|
      !                                                                 |
      !     Reads all general options from namoptions                   |
      !                                                                 |
      !      Chiel van Heerwaarden        15/06/2007                    |
      !      Thijs Heus                   15/06/2007                    |
      !-----------------------------------------------------------------|

    use modglobal,         only : version,initglobal,iexpnr, ltotruntime, runtime, dtmax, dtav_glob,timeav_glob,&
                                  lwarmstart,startfile,trestart,&
                                  nsv,itot,jtot,kmax,xsize,ysize,xlat,xlon,xyear,xday,xtime,&
                                  lcoriol,lpressgrad,igrw_damp,geodamptime,uvdamprate,lmomsubs,cu,cv,&
                                  ifnamopt,fname_options,llsadv, &
                                  ibas_prf,lambda_crit,iadv_mom,iadv_tke,iadv_thl,iadv_qt,iadv_sv,courant,peclet,ladaptive,author,&
                                  lrigidlid,unudge,ntimedep,&
                                  checknamelisterror, &
                                  loutdirs, output_prefix, &
                                  lopenbc,linithetero,lperiodic,dxint,dyint,dzint,dxturb,dyturb,taum,tauh,pbc,&
                                  lsynturb,nmodes,tau,lambda,lambdas,lambdas_x,lambdas_y,lambdas_z,iturb, &
                                  rdt,rk3step,i1,j1,k1,ih,jh,lboundary,iinput,dzf
    use modforces,         only : lforce_user
    use modsurface,        only : initsurface
    use moddatetime,       only : initdatetime
    use modemission,       only : initemission
    use modlsm,            only : initlsm, kmax_soil
    use modslurb,          only : initslurb, preprocess_slurb
    use moddrydeposition,  only : initdrydep
    use modfields,         only : initfields,um,vm,wm,u0,v0,w0,up,vp,wp,rhobf
    use modtracers,        only : inittracers, allocate_tracers, add_tracer
    use modpois,           only : initpois,poisson
    use modradiation,      only : initradiation
    use modraddata,        only : irad,iradiation,&
                                  rad_ls,rad_longw,rad_shortw,rad_smoke,useMcICA,&
                                  timerad,rka,dlwtop,dlwbot,sw0,gc,reff,isvsmoke,lcloudshading
    use modtimedep,        only : inittimedep,ltimedep,ltimedepuv
    use modtimedepsv,      only : inittimedepsv,ltimedepsv
    use modtestbed,        only : inittestbed
    use modboundary,       only : initboundary,ksp
    use modthermodynamics, only : initthermodynamics
    use modmicrophysics,   only : initmicrophysics
    use modsubgrid,        only : initsubgrid
    use modmpi,            only : initmpi,commwrld,myid,myidx,myidy,cmyidy,nprocx,nprocy,mpierr,periods &
                                , D_MPI_BCAST
    use tstep,             only : inittstep
    use modchem,           only : initchem
    use modversion,        only : git_version
    use modopenboundary,   only : initopenboundary,openboundary_divcorr,openboundary_excjs,lbuoytop,&
                                  rhointi, openboundary_phasevelocity
    use modibm,            only : initibm

    use modchecksim,       only : chkdiv
    use modnamelist,       only : read_namelists
    use modspraying,       only : initspraying
    use fortran_support,   only: nnml_output

    implicit none

    character(len=*), parameter :: routine = modname//'/startup'

    integer :: ierr
    logical,dimension(2) :: lper = .false.
    character(256), optional, intent(in) :: path
    real rk3coef


    !declare namelists
    namelist/RUN/ &
        iexpnr,lwarmstart,startfile,ltotruntime, runtime,dtmax,wctime,dtav_glob,timeav_glob,&
        trestart,irandom,randthl,randqt,krand,nsv,courant,peclet,ladaptive,author,&
        krandumin, krandumax, randu,&
        nprocx,nprocy,loutdirs, iinput
    namelist/DOMAIN/ &
        itot,jtot,kmax,kmax_soil,&
        xsize,ysize,&
        xlat,xlon,xyear,xday,xtime,ksp
    namelist/PHYSICS/ &
        !cstep z0,ustin,wtsurf,wqsurf,wsvsurf,ps,thls,chi_half,lmoist,isurf,lneutraldrag,&
        lcoriol,lpressgrad,igrw_damp,geodamptime,uvdamprate,lmomsubs,ltimedep,ltimedepuv,ltimedepsv,ntimedep,&
        irad,timerad,iradiation,rad_ls,rad_longw,rad_shortw,rad_smoke,useMcICA,&
        rka,dlwtop,dlwbot,sw0,gc,reff,isvsmoke,lforce_user,lcloudshading,lrigidlid,unudge
    namelist/DYNAMICS/ &
        llsadv, lambda_crit, cu, cv, ibas_prf, iadv_mom, iadv_tke, iadv_thl, iadv_qt, iadv_sv
    namelist/OPENBC/ &
        lopenbc,linithetero,lper,lbuoytop,dxint,dyint,dzint,dxturb,dyturb,taum,tauh,pbc,lsynturb,iturb,tau,lambda,nmodes,lambdas,lambdas_x,lambdas_y,lambdas_z,lbuoytop


    ! get myid
    ! call MPI_INIT(mpierr)
    ! call MPI_COMM_RANK( MPI_COMM_WORLD, myid, mpierr )

    !read namelists
    if(myid==0)then
      write (*, *) trim(version)//' git: '//trim(git_version)
      if(present(path)) then
          fname_options=path
      else
         if (command_argument_count() >=1) then
            call get_command_argument(1,fname_options)
         end if
      end if
      write (*,*) fname_options

      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      if (ierr /= 0) then
        call finish(routine, 'ERROR:Namoptions does not exist')
      end if
      read (ifnamopt,RUN,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'RUN')
      write(nnml_output ,RUN)
      rewind(ifnamopt)
      read (ifnamopt,DOMAIN,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'DOMAIN')
      write(nnml_output ,DOMAIN)
      rewind(ifnamopt)
      read (ifnamopt,PHYSICS,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'PHYSICS')
      write(nnml_output ,PHYSICS)
      rewind(ifnamopt)
      read (ifnamopt,DYNAMICS,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'DYNAMICS')
      write(nnml_output ,DYNAMICS)
      rewind(ifnamopt)
      read (ifnamopt,OPENBC,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'OPENBC')
      write(nnml_output ,OPENBC)
      close(ifnamopt)
      if(lopenbc) then
        ! Check if grid needs to be periodic
        periods = (/lper(1),lper(2)/)
        lperiodic(1:2) = lper(1)
        lperiodic(3:4) = lper(2)
      endif
      close(ifnamopt)
    end if


    ! these must be shared before initmpi sets up the cartesian grid
    ! commwrld is already set up
    call D_MPI_BCAST(periods,2,0,commwrld,mpierr)
    call D_MPI_BCAST(nprocx ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(nprocy ,1,0,commwrld,mpierr)

    ! Initialize MPI
    call initmpi

    !$acc update device (myidx,myidy)
!$omp target update to(myidx,myidy)

    ! Ignore user-provided nsv, we take care of it ourselves
    nsv = 0

  !broadcast namelists
    call D_MPI_BCAST(iexpnr     ,1,0,commwrld,mpierr) ! RUN
    call D_MPI_BCAST(lwarmstart ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(startfile  ,50,0,commwrld,mpierr)
    call D_MPI_BCAST(author     ,80,0,commwrld,mpierr)
    call D_MPI_BCAST(runtime    ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(trestart   ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(dtmax      ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(dtav_glob  ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(ltotruntime,1,0,commwrld,mpierr)
    call D_MPI_BCAST(wctime     ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(timeav_glob,1,0,commwrld,mpierr)
    call D_MPI_BCAST(nsv        ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(loutdirs   ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(iinput, 1, 0, commwrld, mpierr)

    call D_MPI_BCAST(itot       ,1,0,commwrld,mpierr) ! DOMAIN
    call D_MPI_BCAST(jtot       ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(kmax       ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(kmax_soil  ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(xsize      ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(ysize      ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(xlat       ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(xlon       ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(xyear      ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(xday       ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(xtime      ,1,0,commwrld,mpierr)

    !call D_MPI_BCAST(lneutraldrag ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(lcoriol     ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(lpressgrad  ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(igrw_damp   ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(geodamptime ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(uvdamprate  ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(lforce_user ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(lmomsubs    ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(ntimedep    ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(ltimedep    ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(ltimedepuv  ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(ltimedepsv  ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(lrigidlid   ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(unudge      ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(irad       ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(timerad    ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(iradiation ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(rad_ls     ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(rad_longw  ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(rad_shortw ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(rad_smoke  ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(useMcIca   ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(rka        ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(dlwtop     ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(dlwbot     ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(sw0        ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(gc         ,1,0,commwrld,mpierr)
    ! CvH call D_MPI_BCAST(sfc_albedo ,1,MY_REAL   ,0,commwrld,mpierr)
    call D_MPI_BCAST(reff       ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(isvsmoke   ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(lcloudshading,1,0,commwrld,mpierr)

    call D_MPI_BCAST(llsadv     ,1,0,commwrld,mpierr) ! DYNAMICS
    call D_MPI_BCAST(lambda_crit,1,0,commwrld,mpierr)
    call D_MPI_BCAST(cu         ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(cv         ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(ksp        ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(irandom    ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(krand      ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(krandumin  ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(krandumax  ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(randthl    ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(randqt     ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(randu      ,1,0,commwrld,mpierr)

    call D_MPI_BCAST(ladaptive  ,1,0,commwrld,mpierr) ! RUN
    call D_MPI_BCAST(courant,1,0,commwrld,mpierr)
    call D_MPI_BCAST(peclet,1,0,commwrld,mpierr)

    call D_MPI_BCAST(ibas_prf,1,0,commwrld,mpierr)
    call D_MPI_BCAST(iadv_mom,1,0,commwrld,mpierr)
    call D_MPI_BCAST(iadv_tke,1,0,commwrld,mpierr)
    call D_MPI_BCAST(iadv_thl,1,0,commwrld,mpierr)
    call D_MPI_BCAST(iadv_qt ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(iadv_sv ,1,0,commwrld,mpierr)

    ! Broadcast openboundaries Variables
    call D_MPI_BCAST(lopenbc,    1, 0,commwrld,mpierr)
    call D_MPI_BCAST(linithetero,1, 0,commwrld,mpierr)
    call D_MPI_BCAST(lperiodic,  5, 0,commwrld,mpierr)
    call D_MPI_BCAST(lbuoytop,   1, 0,commwrld,mpierr)
    call D_MPI_BCAST(dxint,      1, 0,commwrld,mpierr)
    call D_MPI_BCAST(dyint,      1, 0,commwrld,mpierr)
    call D_MPI_BCAST(dzint,      1, 0,commwrld,mpierr)
    call D_MPI_BCAST(dxturb,     1, 0,commwrld,mpierr)
    call D_MPI_BCAST(dyturb,     1, 0,commwrld,mpierr)
    call D_MPI_BCAST(taum,       1, 0,commwrld,mpierr)
    call D_MPI_BCAST(tauh,       1, 0,commwrld,mpierr)
    call D_MPI_BCAST(pbc,        1, 0,commwrld,mpierr)
    call D_MPI_BCAST(lsynturb,   1, 0,commwrld,mpierr)
    call D_MPI_BCAST(iturb,      1, 0,commwrld,mpierr)
    call D_MPI_BCAST(lambda,     1, 0,commwrld,mpierr)
    call D_MPI_BCAST(tau,        1, 0,commwrld,mpierr)
    call D_MPI_BCAST(nmodes,     1, 0,commwrld,mpierr)
    call D_MPI_BCAST(lambdas,    1, 0,commwrld,mpierr)
    call D_MPI_BCAST(lambdas_x,  1, 0,commwrld,mpierr)
    call D_MPI_BCAST(lambdas_y,  1, 0,commwrld,mpierr)
    call D_MPI_BCAST(lambdas_z,  1, 0,commwrld,mpierr)

    ! Read all namelists
    call read_namelists(fname_options)

    call testwctime
    ! Allocate and initialize core modules
    call initglobal
    call timer_init
    call timer_tic('modstartup/startup', 0)
    call initfields
    call inittracers
    call initmicrophysics
    call initspraying
    call allocate_tracers ! At this point, all tracers have to be defined
    call inittestbed    !reads initial profiles from scm_in.nc, to be used in readinitfiles
    call inittstep

    call initibm ! keep here as it may overwrite ibas_prf

    if(.not.lopenbc) then
      call initboundary
    else
      call initopenboundary
    endif
    call initthermodynamics
    call initradiation
    call initchem
    call initsurface
    call initdatetime
    call initemission
    call initlsm
    call initdrydep
    call initsubgrid
    call initslurb

    if (loutdirs) then
       output_prefix(1:3) = cmyidy
       output_prefix(4:4) = '/'
       if (myidx == 0) then
            ierr = mkdir(cmyidy//c_null_char, int(o'772',c_int16_t))
       end if
    end if

    call readinitfiles ! moved to obtain the correct btime for the timedependent forcings in case of a warmstart
    call inittimedep !depends on modglobal,modfields, modmpi, modsurf, modradiation, and on modtracers
    call initpois ! hypre solver needs grid and baseprofiles
    if(lopenbc) then  ! Correct boundaries and initial field for divergence
      ! Create 1/int(rho) - must be after rhobf has been initialized
      allocate(rhointi(k1))
      rhointi = 1./(rhobf*dzf)

      call openboundary_phasevelocity() ! needed for initialization, called late in the time loop

      call chkdiv
      call openboundary_divcorr ! Remove divergence from large scale input
      ! Use poisson solver to get rid of divergence in initial field, needs to
      ! be here to avoid cross dependencies between modopenbondaries and modpois
      if(myid==0) print *, 'Start divergence correction initial field'
      call chkdiv
      up = 0.; vp = 0.; wp = 0. ! Set tendencies to zero
      call poisson
      rk3coef = rdt / (4. - dble(rk3step))
      um = um + rk3coef * up
      vm = vm + rk3coef * vp
      wm = wm + rk3coef * wp
      call openboundary_excjs(um   , 2,i1,2,j1,1,k1,ih,jh,.not.lboundary(1:4).or.lperiodic(1:4))
      call openboundary_excjs(vm   , 2,i1,2,j1,1,k1,ih,jh,.not.lboundary(1:4).or.lperiodic(1:4))
      call openboundary_excjs(wm   , 2,i1,2,j1,1,k1,ih,jh,.not.lboundary(1:4).or.lperiodic(1:4))
      u0 = um; v0 = vm; w0 = wm
      call chkdiv
      ! Reset tendencies
      up = 0.; vp = 0.; wp = 0.
      if(myid==0) print *, 'Finished divergence correction initial field'
    endif

    call inittstep
    call preprocess_slurb(lwarmstart)

    call checkinitvalues

    call check_initial_state

    call timer_toc('modstartup/startup')

  end subroutine startup


  !> Checks whether crucial parameters are set correctly
  subroutine checkinitvalues
    use modsurfdata, only: wtsurf, wqsurf, ustin, thls, isurf, ps
    use modglobal,   only: itot, jtot, ysize, xsize, dtmax, runtime, &
                           startfile, lwarmstart, eps1, imax, jmax, ih, jh, &
                           lcoriol, lpressgrad
    use modmpi,      only: myid, nprocx, nprocy, mpierr, MPI_FINALIZE
    use modtimedep,  only: ltimedep

    character(len=*), parameter :: routine = modname//'/checkinitvalues'

    ! Check MPI configuration
    if (mod(jtot, nprocy) /= 0) then
      call finish(routine, 'ERROR: jtot (', jtot, ') is not &
        &divisible by nprocy (', nprocy, '). Please change your MPI &
        &configuration.')
    else
      if(myid == 0) then
        call message(routine, 'jmax = jtot / nprocy = ', jmax)
      end if
    end if

    if (mod(itot, nprocx) /= 0) then
        call finish(routine, 'ERROR: jtot (', itot, ') is not &
          &divisible by nprocy (', nprocx, '). Please change your MPI &
          &configuration.')
    else
      if (myid == 0) then
        call message(routine, 'imax = itot / nprocx = ', imax)
      end if
    end if

    ! Check if we have overlapping ghost cells
    if (ih > imax) then
      call finish(routine, 'ERROR: imax (', imax, ') is smaller &
        &than the required number of ghost cells (', ih, '). Please change &
        &your MPI configuration.')
    end if

    if (jh > jmax) then
      call finish(routine, 'ERROR: jmax (', jmax, ') is smaller &
        &than the required number of ghost cells (', jh, '). Please change &
        &your MPI configuration.')
    end if

    ! Check namoptions
    if (runtime < 0) call finish(routine, 'runtime out of range/not set')
    if (dtmax < 0) call finish(routine, 'dtmax out of range/not set')
    if (ps < eps1) call finish(routine, 'psout of range/not set')
    if (thls < eps1) call finish(routine, 'thls out of range/not set')
    if (xsize < 0) call finish(routine, 'xsize out of range/not set')
    if (ysize < 0) call finish(routine, 'ysize out of range/not set')

    if (lwarmstart) then
      if (startfile == '') call finish(routine, 'no restartfile set')
    end if

    ! Surface
    if (myid == 0) then
      select case (isurf)
        case (1)
        case (2,10)
        case (3:4)
          if (wtsurf < -1E10) call finish(routine, 'wtsurf not set')
          if (wqsurf < -1E10) call finish(routine, 'wqsurf not set')
        case (11)
        case default
          call finish(routine, 'isurf out of range/not set')
      end select

      if (isurf == 3) then
        if (ustin < 0) call finish(routine, 'ustin out of range/not set')
      end if
    end if

    if (lcoriol .and. lpressgrad) then
      if (myid==0) call finish(routine, "Coriolis force (lcoriol) and channel-like pressure gradient (lpressgrad) are mutually exclusive. To use Coriolis force with NO pressure gradient, set geowinds to zero.")
   end if
  end subroutine checkinitvalues

  subroutine readinitfiles
    use modfields,         only : u0,v0,w0,um,vm,wm,thlm,thl0,thl0h,qtm,qt0,qt0h,&
                                  ql0,ql0h,thv0h,sv0,svm,e12m,e120,&
                                  dudxls,dudyls,dvdxls,dvdyls,dthldxls,dthldyls,&
                                  dqtdxls,dqtdyls,dqtdtls,dpdxl,dpdyl,&
                                  wfls,whls,ug,vg,uprof,vprof,thlprof, qtprof,e12prof, svprof,&
                                  v0av,u0av,qt0av,ql0av,thl0av,sv0av,exnf,exnh,presf,presh,initial_presf,initial_presh,rhof,&
                                  thlpcar,thvh,thvf
    use modglobal,         only : i1,i2,ih,j1,j2,jh,kmax,k1,dtmax,idtmax,dt,rdt,runtime,timeleft,tres,&
                                  rtimee,timee,ntrun,btime,dt_lim,nsv,&
                                  zf,dzf,dzh,rv,rd,cp,rlv,pref0,om23_gs,&
                                  ijtot,cu,cv,e12min,dzh,cexpnr,ifinput,lwarmstart,ltotruntime,itrestart,&
                                  trestart, ladaptive,llsadv,tnextrestart,longint,lopenbc,linithetero, &
                                  iinput, input_netcdf, input_ascii, lcoriol, &
                                  dzhi, iadv_thl, iadv_qt, iadv_kappa, eps1
    use modthermodynamics, only : lconstexner,lbaseexner
    use modsubgrid,        only : ekm,ekh
    use modsurfdata,       only : wsvsurf, &
                                  thls,tskin,tskinm,tsoil,tsoilm,phiw,phiwm,Wl,Wlm,thvs,qts,isurf,svs,obl,oblav,&
                                  qskin
    use modsurface,        only : surface,qtsurf,dthldz,ps
    use modlsm,            only : init_lsm_tiles
    use modboundary,       only : boundary
    use modmpi,            only : slabsum,myid,comm3d,mpierr,D_MPI_BCAST, print_info_stderr
    use modthermodynamics, only : thermodynamics,calc_halflev, lmoist
    use moduser,           only : initsurf_user
    use modibmdata,        only : thlibm, qtibm, lapply_ibm, fluid_mask
    use modtestbed,        only : ltestbed,tb_ps,tb_thl,tb_qt,tb_u,tb_v,tb_w,tb_ug,tb_vg,&
                                  tb_dqtdxls,tb_dqtdyls,tb_qtadv,tb_thladv
    use modopenboundary,   only : openboundary_ghost,openboundary_readboundary,openboundary_initfields
    use modtracers,        only : tracer_prop, tracer_profs_from_netcdf, nsv_user
    use utils,             only : to_lower
    use modslabaverage,    only : slabavg
    use modlogging,        only : profile_output

#if defined(DALES_GPU)
    use modgpu, only: update_gpu, update_host, host_is_updated, update_gpu_surface
#endif

    character(len=*), parameter :: routine = modname//"/readinitfiles"

    integer i,j,k,n,ierr
    integer isv, isv_u
    logical negval !switch to allow or not negative values in randomnization

    real(field_r), allocatable :: height(:), th0av(:)
    real(field_r), allocatable :: thv0(:,:,:)

    character(len=512) :: chmess
    character(len=512) :: header_line
    character(len=16)  :: header
    integer, parameter :: maxcol = 50
    integer            :: header_pos(maxcol)
    integer            :: header_len(maxcol)
    logical            :: found
    real               :: vals_at_lev(maxcol)
    integer            :: nheader
    integer            :: start, end

    allocate (height(k1))
    allocate (th0av(k1))
    allocate (thv0(2-ih:i1+ih,2-jh:j1+jh,k1))

    if (.not. lwarmstart) then

    !********************************************************************

    !    1.0 prepare initial fields from files 'prog.inp' and 'scalar.inp'
    !    ----------------------------------------------------------------

    !--------------------------------------------------------------------
    !    1.1 read fields
    !-----------------------------------------------------------------

      rdt = dtmax / 100.
      dt  = max(floor(rdt/tres), 1)
      timee = 0
      if (myid==0) then

        if (ltestbed) then

          write(*,*) 'readinitfiles: testbed mode: profiles for initialization obtained from scm_in.nc'

          do k=1,kmax
            height (k) = zf(k)
            thlprof(k) = tb_thl(1,k)
            qtprof (k) = tb_qt(1,k)
            uprof  (k) = tb_u(1,k)
            vprof  (k) = tb_v(1,k)
            e12prof(k) = e12min
          end do

          ps         = tb_ps(1)

        else if (iinput == input_netcdf) then
          call init_from_netcdf('init.'//cexpnr//'.nc', height, uprof, vprof, &
                                thlprof, qtprof, e12prof, ug, vg, dpdxl, dpdyl, wfls, &
                                dqtdxls, dqtdyls, dqtdtls, thlpcar, kmax)
          if (nsv_user > 0) then
            call tracer_profs_from_netcdf('tracers.'//cexpnr//'.nc', &
                                          tracer_prop, svprof(1:kmax,:))
          end if
        else
          open (ifinput,file='prof.inp.'//cexpnr,status='old',iostat=ierr)
          if (ierr /= 0) then
             call finish(routine, 'Cannot open the file ', 'prof.inp.'//cexpnr)
          end if
          read (ifinput,'(a512)') chmess
          write(*,     '(a512)') chmess
          read (ifinput,'(a512)') chmess

          do k = 1, kmax
            read (ifinput,*) &
                height (k), &
                thlprof(k), &
                qtprof (k), &
                uprof  (k), &
                vprof  (k), &
                e12prof(k)
          end do

          close(ifinput)
        end if   !ltestbed

        write(profile_output,*) 'height    thl      qt         u      v     e12'
        do k = kmax, 1, -1
          write (profile_output,'(f7.1,f8.1,e12.4,3f7.1)') &
                zf     (k), &
                thlprof(k), &
                qtprof (k), &
                uprof  (k), &
                vprof  (k), &
                e12prof(k)
        end do

        if (minval(e12prof(1:kmax)) < e12min) then
          write(*,*)  'e12 value is zero (or less) in prof.inp'
          do k = 1, kmax
            e12prof(k) = max(e12prof(k),e12min)
          end do
        end if

      end if ! end if myid==0

      ! MPI broadcast numbers reading
      call D_MPI_BCAST(thlprof,kmax,0,comm3d,mpierr)
      call D_MPI_BCAST(qtprof ,kmax,0,comm3d,mpierr)
      call D_MPI_BCAST(uprof  ,kmax,0,comm3d,mpierr)
      call D_MPI_BCAST(vprof  ,kmax,0,comm3d,mpierr)
      call D_MPI_BCAST(e12prof,kmax,0,comm3d,mpierr)

      if(myid==0)then
        if (nsv_user>0 .and. iinput == input_ascii) then
          open (ifinput,file='scalar.inp.'//cexpnr,status='old',iostat=ierr)
          if (ierr /= 0) then
             call finish(routine,'Cannot open the file ', 'scalar.inp.'//cexpnr)
          end if

          ! reading header (2 lines)
          read (ifinput,'(a512)') chmess
          read (ifinput,'(a512)') chmess

          header_line(:) = chmess(:)

          call split_string(header_line, nheader, header_pos, header_len)

          ! Try to find profiles
          do isv = 1, nsv
            found = .false.
            do isv_u = 1, nheader
              start = header_pos(isv_u)
              end = header_pos(isv_u) + header_len(isv_u) - 1
              header = header_line(start:end)
              if (trim(tracer_prop(isv)%tracname) == trim(header)) then
                do k = 1, kmax
                  read(ifinput, *, iostat=ierr) vals_at_lev(1:nheader)
                  svprof(k,isv) = vals_at_lev(isv_u)
                end do
                found = .true.
                ! Go back to the start of the file
                rewind(ifinput)
                read(ifinput,'(a512)') chmess
                read(ifinput,'(a512)') chmess
              end if
            end do

            if (.not. found) then
              call warning(routine, "no initial profile found for &
                & "//tracer_prop(isv)%tracname)
            end if
          end do

          close(ifinput)
        end if
      end if

      if (myid == 0) then
        ! Print tracer profiles to stderr
        write(profile_output, '(a9)', advance='no') 'height   '
        do isv = 1, nsv
          write(profile_output, '(a12)', advance='no') tracer_prop(isv)%tracname
        end do

        write(profile_output, *)

        do k = kmax, 1, -1
          write(profile_output,'(f7.1,2x)', advance='no') height(k)
          do isv = 1, nsv
            write(profile_output, '(e10.4,2x)', advance='no') svprof(k,isv)
          end do
          write(profile_output, *)
        end do
      end if

      call D_MPI_BCAST(wsvsurf,        nsv,    0, comm3d, mpierr)
      call D_MPI_BCAST(svprof,         k1*nsv, 0, comm3d, mpierr)

      ! Initialize fields
      if(lopenbc .and. linithetero) then! Openboundaries with heterogeneous initialisation
        call openboundary_initfields(tracer_prop)
        do j = 1,j2
          do i = 1,i2
            wm(i,j,1) = 0.
            w0(i,j,1) = 0.
          end do
        end do
        do k = 1,kmax
          do j = 1,j2
            do i = 1,i2
               ekm(i,j,k) = 0.0
               ekh(i,j,k) = 0.0
            end do
          end do
        end do
      else
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
        if(nsv>0) then
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
        endif
      endif
    !---------------------------------------------------------------
    !  1.2 randomnize fields
    !---------------------------------------------------------------

      krand  = min(krand,kmax)
      negval = .False. ! No negative perturbations for qt (negative moisture is non physical)

      if (irandom < 0) then
         call randomize_new(irandom)
      else
         do k = 1,krand
            call randomnize(qtm ,k,randqt ,irandom,ih,jh,negval)
            call randomnize(qt0 ,k,randqt ,irandom,ih,jh,negval)
         end do
         negval = .True. ! negative perturbations allowed
         do k = 1,krand
            call randomnize(thlm,k,randthl,irandom,ih,jh,negval)
            call randomnize(thl0,k,randthl,irandom,ih,jh,negval)
         end do

         do k=krandumin,krandumax
            call randomnize(um  ,k,randu  ,irandom,ih,jh,negval)
            call randomnize(u0  ,k,randu  ,irandom,ih,jh,negval)
            call randomnize(vm  ,k,randu  ,irandom,ih,jh,negval)
            call randomnize(v0  ,k,randu  ,irandom,ih,jh,negval)
            call randomnize(wm  ,k,randu  ,irandom,ih,jh,negval)
            call randomnize(w0  ,k,randu  ,irandom,ih,jh,negval)
         end do
      end if

      ! when using ibm, overwrite randomnization inside obtacles (velocities, thl, qt)
      if (lapply_ibm) then
        do k=1,kmax
          do j=2,j1
            do i=2,i1
              if (.not.(fluid_mask(i,j,k))) then
                thlm(i,j,k)     = thlibm !set to thlibm value
                thl0(i,j,k)     = thlibm
                qtm(i,j,k)      = qtibm  !set to qtibm value
                qt0(i,j,k)      = qtibm
                um (i:i+1,j,k)  = 0.
                u0 (i:i+1,j,k)  = 0.
                vm (i,j:j+1,k)  = 0.
                v0 (i,j:j+1,k)  = 0.
                wm (i,j,k:k+1)  = 0.
                w0 (i,j,k:k+1)  = 0.
                e12m(i,j,k)     = e12min
                e120(i,j,k)     = e12min
                if (nsv > 0) then
                  do n=1,nsv
                    sv0(i,j,k,n) = 0.
                    svm(i,j,k,n) = 0.
                  end do
                end if
              end if
            end do
          end do
        end do
      end if

      !--------------------------------------------------------------------------
      !    2.2 Check surface settings, initialize surface layer and base profiles
      !--------------------------------------------------------------------------

      ! We call thermodynamics to calculate qtsurf, for which we need to know ps is set correctly.
      if (ps < eps1) call finish(routine, 'ps out of range/not set')

      select case(isurf)
      case(1)
        tskin(:,:)  = thls
        tskinm(:,:) = tskin(:,:)
        tsoilm(:,:,:) = tsoil(:,:,:)
        phiwm(:,:,:)  = phiw(:,:,:)
        Wlm(:,:)    = Wl(:,:)
      case(2)
        tskin(:,:)  = thls
      case(3,4)
        thls   = thlprof(1)
        qts    = qtprof(1)
        tskin(:,:)  = thls
        qskin(:,:)  = qts
      case(11)
        thls   = thlprof(1)
        qts    = qtprof(1)
        tskin(:,:)  = thls
        qskin(:,:)  = qts
        call init_lsm_tiles
      case(10)
        call initsurf_user
      end select
      if(lopenbc) then
        call openboundary_readboundary(tracer_prop)
        call openboundary_ghost
      endif

      ! Set initial Obukhov length to -0.1 for iteration
      obl(:,:) = -0.1
      oblav    = -0.1

      ! qtsurf act on device data.
#if defined(_OPENACC)
      call update_gpu_surface
#endif
      if (lmoist) call qtsurf

      dthldz(:,:) = (thlprof(1) - thls) / zf(1)
      thvs = thls * (1. + (rv/rd - 1.) * qts)

      thl0av(:) = thlprof(:)   ! these are used for the top boundary in modboundary
      qt0av(:)  = qtprof(:)    ! but have not been initialized yet (?)
      sv0av(:,:) = svprof(:,:) !
                               ! TODO: OpenACC??

      u0av(:)   = uprof(:)
      v0av(:)   = vprof(:)

      svs = svprof(1,:)

      call baseprofs ! call baseprofs before thermodynamics

#if defined(DALES_GPU)
      call update_gpu
#endif
      if ( lopenbc ) then
        call openboundary_ghost()
      else
        call boundary(on_gpu=.true.)
      end if

      call thermodynamics
      call surface

      if ( lopenbc ) then
        call openboundary_ghost()
      else
        call boundary(on_gpu=.true.)
      end if

      call thermodynamics

#if defined(DALES_GPU)
      call update_host
      host_is_updated = .false.
#endif

      ! save initial pressure profiles
      ! used for initialising radiation scheme at restart, to reproduce the same state
      initial_presf(:) = presf(:)
      initial_presh(:) = presh(:)

    else !if lwarmstart

      call readrestartfiles
      call baseprofs

      um(:,:,:) = u0(:,:,:)
      vm(:,:,:) = v0(:,:,:)
      wm(:,:,:) = w0(:,:,:)
      thlm(:,:,:) = thl0(:,:,:)
      qtm(:,:,:)  = qt0(:,:,:)
      svm(:,:,:,:)  = sv0(:,:,:,:)
      e12m(:,:,:) = e120(:,:,:)

#if defined(DALES_GPU)
      call update_gpu
#endif

      call calc_halflev(thl0, dzf, dzhi, thls, iadv_thl == iadv_kappa, thl0h)
      call calc_halflev(qt0, dzf, dzhi, qts, iadv_qt == iadv_kappa, qt0h)

#if defined(DALES_GPU)
      call update_host
      host_is_updated = .false.
#endif

      if (.not. lbaseexner) then
         if (lconstexner) then
            exnf(:) = (initial_presf(:)/pref0)**(rd/cp)
            exnh(:) = (initial_presh(:)/pref0)**(rd/cp)
         else
            exnf(:) = (presf(:)/pref0)**(rd/cp)
            exnh(:) = (presh(:)/pref0)**(rd/cp)
         endif
      endif

      do k = 2, k1
        do j = 2, j1
          do i = 2, i1
            thv0h(i,j,k) = (thl0h(i,j,k)+rlv*ql0h(i,j,k)/(cp*exnh(k))) &
                          *(1+(rv/rd-1)*qt0h(i,j,k)-rv/rd*ql0h(i,j,k))
          end do
        end do
      end do

      do k = 1, k1
        do j = 2, j1
          do i = 2 ,i1
            thv0(i,j,k) = (thl0(i,j,k)+rlv*ql0(i,j,k)/(cp*exnf(k))) &
                          *(1+(rv/rd-1)*qt0(i,j,k)-rv/rd*ql0(i,j,k))
          end do
        end do
      end do

      thvh(:) = 0.0
      thvf(:) = 0.0
      u0av(:) = 0.0
      v0av(:) = 0.0
      thl0av(:) = 0.0
      th0av(:) = 0.0
      qt0av(:) = 0.0
      ql0av(:) = 0.0
      sv0av(:,:) = 0.0

      if (.not.(lapply_ibm)) then
        call slabsum(thvh,1,k1,thv0h,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1) ! redefine halflevel thv using calculated thv
        call slabsum(thvf,1,k1,thv0,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)

        call slabsum(u0av  ,1,k1,u0  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
        call slabsum(v0av  ,1,k1,v0  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
        call slabsum(thl0av,1,k1,thl0,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
        call slabsum(qt0av ,1,k1,qt0 ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
        call slabsum(ql0av ,1,k1,ql0 ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
        do n = 1, nsv
          call slabsum(sv0av(1:1,n),1,k1,sv0(:,:,:,n),2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
        end do

        thvh(:) = thvh(:) / ijtot
        thvf(:) = thvf(:) / ijtot
        u0av(:) = u0av(:) / ijtot + cu
        v0av(:) = v0av(:) / ijtot + cv
        thl0av(:) = thl0av(:) / ijtot
        qt0av(:) = qt0av(:) / ijtot
        ql0av(:) = ql0av(:) / ijtot
        sv0av(:,:) = sv0av(:,:) / ijtot
      else
        call slabavg(thv0h,fluid_mask,ih,thvh)
        call slabavg(thv0,fluid_mask,ih,thvf)
        call slabavg(u0,fluid_mask,ih,u0av)
        call slabavg(v0,fluid_mask,ih,v0av)
        call slabavg(thl0,fluid_mask,ih,thl0av)
        call slabavg(qt0,fluid_mask,ih,qt0av)
        call slabavg(ql0,fluid_mask,ih,ql0av)
        do n=1,nsv
          call slabavg(sv0(:,:,:,n),fluid_mask,ih,sv0av(:,n))
        end do
      end if

      th0av(:) = thl0av(:) + (rlv/cp) * ql0av(:) / exnf(:)
      thvh(1) = th0av(1)*(1+(rv/rd-1)*qt0av(1)-rv/rd*ql0av(1)) ! override first level
      rhof(:) = presf(:)/(rd*thvf(:)*exnf(:))

      ! CvH - only do this for fixed timestepping. In adaptive dt comes from restartfile
      if(ladaptive .eqv. .false.) rdt=dtmax

      call baseprofs !call baseprofs
      if(lopenbc) then
        call openboundary_readboundary(tracer_prop)
      endif

#if defined(DALES_GPU)
      call update_gpu
#endif

    end if  ! end if (.not. warmstart)

!-----------------------------------------------------------------
!    2.1 read and initialise fields
!-----------------------------------------------------------------


    if(myid==0)then

      if (ltestbed) then

          write(*,*) 'readinitfiles: testbed mode: profiles for ls forcing obtained from scm_in.nc'

          do k=1,kmax
            height (k) = zf(k)
            ug     (k) = tb_ug(1,k)
            vg     (k) = tb_vg(1,k)
            wfls   (k) = tb_w(1,k)
            dqtdxls(k) = tb_dqtdxls(1,k)
            dqtdyls(k) = tb_dqtdyls(1,k)
            dqtdtls(k) = tb_qtadv(1,k)
            thlpcar(k) = tb_thladv(1,k)
          end do

      else

        if (iinput == input_netcdf) then
          continue ! Profiles have been read by init_from_netcdf
        else
          open (ifinput,file='lscale.inp.'//cexpnr, status='old',iostat=ierr)
          if (ierr /= 0) then
             call finish(routine,'Cannot open the file ', 'lscale.inp.'//cexpnr)
          end if
          read (ifinput,'(a80)') chmess
          read (ifinput,'(a80)') chmess

          ! if coriolis force, read in 2nd and 3rd columns as ug and vg
          if (lcoriol) then
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
          else ! otherwhise read in same columns as pressure gradient
            do  k=1,kmax
              read (ifinput,*) &
                  height (k), &
                  dpdxl  (k), &
                  dpdyl  (k), &
                  wfls   (k), &
                  dqtdxls(k), &
                  dqtdyls(k), &
                  dqtdtls(k), &
                  thlpcar(k)
            end do
          end if
          close(ifinput)
        end if

      end if

      if (lcoriol) then
        write(profile_output,*) ' height u_geo   v_geo    subs     ' &
                  ,'   dqtdx      dqtdy        dqtdtls     thl_rad '
        do k=kmax,1,-1
          write (profile_output,'(3f7.1,5e12.4)') &
                zf     (k), &
                ug     (k), &
                vg     (k), &
                wfls   (k), &
                dqtdxls(k), &
                dqtdyls(k), &
                dqtdtls(k), &
                thlpcar(k)
        end do
      else
        write(profile_output,*) ' height u_geo   v_geo    subs     ' &
        ,'   dqtdx      dqtdy        dqtdtls     thl_rad '
        do k=kmax,1,-1
          write (profile_output,'(3f7.1,5e12.4)') &
                zf     (k), &
                dpdxl  (k), &
                dpdyl  (k), &
                wfls   (k), &
                dqtdxls(k), &
                dqtdyls(k), &
                dqtdtls(k), &
                thlpcar(k)
        end do
      end if

    end if ! end myid==0

    ! MPI broadcast variables read in

    call D_MPI_BCAST(ug       ,kmax,0,comm3d,mpierr)
    call D_MPI_BCAST(vg       ,kmax,0,comm3d,mpierr)
    call D_MPI_BCAST(dpdxl    ,kmax,0,comm3d,mpierr)
    call D_MPI_BCAST(dpdyl    ,kmax,0,comm3d,mpierr)
    call D_MPI_BCAST(wfls     ,kmax,0,comm3d,mpierr)
    call D_MPI_BCAST(dqtdxls  ,kmax,0,comm3d,mpierr)
    call D_MPI_BCAST(dqtdyls  ,kmax,0,comm3d,mpierr)
    call D_MPI_BCAST(dqtdtls  ,kmax,0,comm3d,mpierr)
    call D_MPI_BCAST(thlpcar  ,kmax,0,comm3d,mpierr)

    !-----------------------------------------------------------------
    !    2.3 make large-scale horizontal pressure gradient
    !-----------------------------------------------------------------

    !******include rho if rho = rho(z) /= 1.0 ***********
    if (lcoriol) then  ! only when Coriolis is enabled, calculte pressure gradients via geostrophic wind speeds (otherwise assume they have been set already)
       do k = 1, kmax
          dpdxl(k) =  om23_gs*vg(k)
          dpdyl(k) = -om23_gs*ug(k)
       end do
    end if
    !-----------------------------------------------------------------
    !    2.5 make large-scale horizontal gradients
    !-----------------------------------------------------------------

    whls(1)  = 0.0
    do k = 2, kmax
      whls(k) = ( wfls(k)*dzf(k-1) +  wfls(k-1)*dzf(k) )/(2*dzh(k))
    end do
    whls(k1) = (wfls(kmax)+0.5*dzf(kmax)*(wfls(kmax)-wfls(kmax-1)) &
                                                  /dzh(kmax))

    !******include rho if rho = rho(z) /= 1.0 ***********

    if (llsadv) then
      if (myid==0) call finish(routine, 'llsadv should not be used anymore. Large scale gradients were calculated in a non physical way (and lmomsubs had to be set to true to retain conservation of mass)')
    end if
    dudxls   = 0.0
    dudyls   = 0.0
    dvdxls   = 0.0
    dvdyls   = 0.0
    dthldxls = 0.0
    dthldyls = 0.0

    idtmax = floor(dtmax/tres)
    btime   = timee
    if (.not.(ltotruntime)) then
      runtime = runtime + btime*tres
    end if
    timeleft=ceiling((runtime)/tres-btime,longint)

    dt_lim = timeleft
    rdt = real(dt)*tres
    ntrun   = 0
    rtimee  = real(timee)*tres
    itrestart = floor(trestart/tres,longint)
    tnextrestart = btime + itrestart

    deallocate (height,th0av,thv0)

  end subroutine readinitfiles

  subroutine readrestartfiles

    use modsurfdata, only : ustar,thlflux,qtflux,svflux,dthldz,dqtdz,ps,thls,qts,thvs,oblav,&
                           tsoil,phiw,tskin,Wl,isurf,ksoilmax,Qnet,swdavn,swuavn,lwdavn,lwuavn,nradtime,&
                           obl,qskin
    use modraddata, only: iradiation,useMcICA, tnext_radiation => tnext, &
                          thlprad,swd,swu,lwd,lwu,swdca,swuca,lwdca,lwuca,swdir,swdif,lwc,&
                          SW_up_TOA,SW_dn_TOA,LW_up_TOA,LW_dn_TOA,&
                          SW_up_ca_TOA,SW_dn_ca_TOA,LW_up_ca_TOA,LW_dn_ca_TOA
    use modfields,  only : u0,v0,w0,thl0,qt0,ql0,ql0h,e120,dthvdz,presf,presh,initial_presf,initial_presh,sv0,tmp0,esl,qvsl,qvsi
    use modglobal,  only : i1,i2,ih,j1,j2,jh,k1,startfile,timee,&
                           tres,ifinput,nsv,dt,output_prefix
    use modboundary, only: dqt, dtheta, dsv
    use modmpi,     only : myid, cmyid
    use modsubgriddata, only : ekm,ekh
    use modlsm, only : kmax_soil, tile, nlu
    use modslurb, only : enable_slurb, slurb_tile
    use modslurbdata, only : moist_physics
    use modlogging, only : finish


    character(len=*), parameter :: routine = modname//'/readrestartfiles'
    character(50) :: name
    integer i,j,k,n, ilu
    !********************************************************************

  !    1.0 Read initfiles
  !-----------------------------------------------------------------
    name = startfile
    name(5:5) = 'd'
    name(14:21)=cmyid
    if (myid == 0) write(6,*) 'loading ',name
    open(unit=ifinput,file=trim(output_prefix)//name,form='unformatted', status='old')

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
      read(ifinput)  (((ekh   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
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
      read(ifinput)  (  initial_presf (    k)                            ,k=1,k1)
      read(ifinput)  (  initial_presh (    k)                            ,k=1,k1)
      read(ifinput)  ps,thls,qts,thvs,oblav
      read(ifinput)  dtheta,dqt,timee,dt,tres
      read(ifinput)   ((obl (i,j  ),i=1,i2      ),j=1,j2      )
      read(ifinput)   ((tskin(i,j ),i=1,i2      ),j=1,j2      )
      read(ifinput)   ((qskin(i,j ),i=1,i2      ),j=1,j2      )

!!!!! radiation quantities
      read(ifinput)  tnext_radiation
      read(ifinput)  (((thlprad (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      read(ifinput)  (((swd     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      read(ifinput)  (((swu     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      read(ifinput)  (((lwd     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      read(ifinput)  (((lwu     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      read(ifinput)  (((swdca   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      read(ifinput)  (((swuca   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      read(ifinput)  (((lwdca   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      read(ifinput)  (((lwuca   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      read(ifinput)  (((swdir   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      read(ifinput)  (((swdif   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      read(ifinput)  (((lwc     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)

      read(ifinput)  ((SW_up_TOA    (i,j ),i=1,i2),j=1,j2)
      read(ifinput)  ((SW_dn_TOA    (i,j ),i=1,i2),j=1,j2)
      read(ifinput)  ((LW_up_TOA    (i,j ),i=1,i2),j=1,j2)
      read(ifinput)  ((LW_dn_TOA    (i,j ),i=1,i2),j=1,j2)
      read(ifinput)  ((SW_up_ca_TOA (i,j ),i=1,i2),j=1,j2)
      read(ifinput)  ((SW_dn_ca_TOA (i,j ),i=1,i2),j=1,j2)
      read(ifinput)  ((LW_up_ca_TOA (i,j ),i=1,i2),j=1,j2)
      read(ifinput)  ((LW_dn_ca_TOA (i,j ),i=1,i2),j=1,j2)
!!!!! end of radiation quantities

    close(ifinput)

    if (nsv>0) then
      name(5:5) = 's'
      if (myid == 0) write(6,*) 'loading ',name
      open(unit=ifinput,file=trim(output_prefix)//name,form='unformatted')
      read(ifinput) ((((sv0(i,j,k,n),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1),n=1,nsv)
      read(ifinput) (((svflux(i,j,n),i=1,i2),j=1,j2),n=1,nsv)
      read(ifinput) (dsv(n),n=1,nsv)
      read(ifinput)  timee
      close(ifinput)
    end if

    if (isurf == 1) then
      name(5:5) = 'l'
      if (myid == 0) write(6,*) 'loading ',name
      open(unit=ifinput,file=trim(output_prefix)//name,form='unformatted')
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

    else if (isurf == 11) then
      name(5:5) = 'l'
      if (myid == 0) write(6,*) 'loading ',name
      open(unit=ifinput,file=trim(output_prefix)//name,form='unformatted')
      read(ifinput) (((tsoil(i,j,k), i=1,i2), j=1,j2), k=1,kmax_soil)
      read(ifinput) (((phiw (i,j,k), i=1,i2), j=1,j2), k=1,kmax_soil)
      read(ifinput) ((tskin (i,j),   i=1,i2), j=1,j2)
      read(ifinput) ((Wl    (i,j),   i=1,i2), j=1,j2)

      do ilu=1,nlu
        read(ifinput) ((tile(ilu)%thlskin(i,j), i=1,i2), j=1,j2)
        read(ifinput) ((tile(ilu)%qtskin(i,j), i=1,i2), j=1,j2)
        read(ifinput) ((tile(ilu)%obuk(i,j), i=1,i2), j=1,j2)
      end do

      if (enable_slurb) then
        read(ifinput) ((slurb_tile%t_can_0(i,j), i=1,i2), j=1,j2)
        read(ifinput) ((slurb_tile%t_can_m(i,j), i=1,i2), j=1,j2)
        read(ifinput) (((slurb_tile%t_wall_a_0(k,i,j), i=1,i2), j=1,j2), k=1,size(slurb_tile%t_wall_a_0,1))
        read(ifinput) (((slurb_tile%t_wall_a_m(k,i,j), i=1,i2), j=1,j2), k=1,size(slurb_tile%t_wall_a_m,1))
        read(ifinput) (((slurb_tile%t_wall_b_0(k,i,j), i=1,i2), j=1,j2), k=1,size(slurb_tile%t_wall_b_0,1))
        read(ifinput) (((slurb_tile%t_wall_b_m(k,i,j), i=1,i2), j=1,j2), k=1,size(slurb_tile%t_wall_b_m,1))
        read(ifinput) (((slurb_tile%t_win_a_0(k,i,j), i=1,i2), j=1,j2), k=1,size(slurb_tile%t_win_a_0,1))
        read(ifinput) (((slurb_tile%t_win_a_m(k,i,j), i=1,i2), j=1,j2), k=1,size(slurb_tile%t_win_a_m,1))
        read(ifinput) (((slurb_tile%t_win_b_0(k,i,j), i=1,i2), j=1,j2), k=1,size(slurb_tile%t_win_b_0,1))
        read(ifinput) (((slurb_tile%t_win_b_m(k,i,j), i=1,i2), j=1,j2), k=1,size(slurb_tile%t_win_b_m,1))
        read(ifinput) (((slurb_tile%t_roof_0(k,i,j), i=1,i2), j=1,j2), k=1,size(slurb_tile%t_roof_0,1))
        read(ifinput) (((slurb_tile%t_roof_m(k,i,j), i=1,i2), j=1,j2), k=1,size(slurb_tile%t_roof_m,1))
        read(ifinput) (((slurb_tile%t_road_0(k,i,j), i=1,i2), j=1,j2), k=1,size(slurb_tile%t_road_0,1))
        read(ifinput) (((slurb_tile%t_road_m(k,i,j), i=1,i2), j=1,j2), k=1,size(slurb_tile%t_road_m,1))
        read(ifinput) ((slurb_tile%q_can_0(i,j), i=1,i2), j=1,j2)
        read(ifinput) ((slurb_tile%q_can_m(i,j), i=1,i2), j=1,j2)
        read(ifinput) ((slurb_tile%m_liq_roof_0(i,j), i=1,i2), j=1,j2)
        read(ifinput) ((slurb_tile%m_liq_roof_m(i,j), i=1,i2), j=1,j2)
        read(ifinput) ((slurb_tile%m_liq_road_0(i,j), i=1,i2), j=1,j2)
        read(ifinput) ((slurb_tile%m_liq_road_m(i,j), i=1,i2), j=1,j2)
      end if

      read(ifinput) timee

      close(ifinput)
    end if

  end subroutine readrestartfiles
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ! this function is called from time stepping,
  ! determines when to write a restart file, then calls do_writerestartfiles to do the work
  !  if trestart = 0, no periodic restart files will be written.
  subroutine writerestartfiles
    use modglobal, only : trestart,itrestart,tnextrestart,dt_lim,timee,timeleft,rk3step
#if defined(DALES_GPU)
    use modgpu, only: update_host
#endif
    implicit none

    if (timee == 0) return
    if (rk3Step/=3) return

    if (timee<tnextrestart) dt_lim = min(dt_lim,tnextrestart-timee)

    ! if trestart > 0, write a restartfile every trestart seconds and at the end
    ! if trestart = 0, write restart files only at the end of the simulation
    ! if trestart < 0, don't write any restart files
    if ((timee>=tnextrestart .and. trestart > 0) .or. (timeleft==0 .and. trestart >= 0)) then
      tnextrestart = tnextrestart+itrestart
#if defined(DALES_GPU)
      call update_host
#endif
      call do_writerestartfiles
    end if
  end subroutine writerestartfiles

  ! this function writes a restart file
  ! separated from writerestartfiles to be callable from the library interface
  subroutine do_writerestartfiles
    use modsurfdata,only: ustar,thlflux,qtflux,svflux,dthldz,dqtdz,ps,thls,qts,thvs,oblav,&
                          tsoil,phiw,tskin,Wl,ksoilmax,isurf,ksoilmax,Qnet,swdavn,swuavn,lwdavn,lwuavn,nradtime,&
                          obl,qskin
    use modraddata, only: iradiation,useMcICA, tnext_radiation => tnext, &
                          thlprad,swd,swu,lwd,lwu,swdca,swuca,lwdca,lwuca,swdir,swdif,lwc,&
                          SW_up_TOA,SW_dn_TOA,LW_up_TOA,LW_dn_TOA,&
                          SW_up_ca_TOA,SW_dn_ca_TOA,LW_up_ca_TOA,LW_dn_ca_TOA

    use modfields, only : u0,v0,w0,thl0,qt0,ql0,ql0h,e120,dthvdz,presf,presh,initial_presf,initial_presh,sv0,tmp0,esl,qvsl,qvsi
    use modglobal, only : i1,i2,ih,j1,j2,jh,k1,cexpnr,ifoutput,timee,rtimee,tres,nsv,dt,output_prefix
    use modboundary, only: dqt, dtheta, dsv
    use modmpi,    only : cmyid,myid
    use modsubgriddata, only : ekm,ekh
    use modlsm,    only : kmax_soil, tile, nlu
    use modslurb,  only : enable_slurb, slurb_tile
    use modslurbdata, only : moist_physics

    implicit none
    integer imin,ihour
    integer i,j,k,n, ilu
    character(50) name,linkname

      ihour = floor(rtimee/3600)
      imin  = floor((rtimee-ihour * 3600) /3600. * 60.)
      name = 'initdXXXXhXXmXXXXXXXX.XXX'
      write (name(6:9)  ,'(i4.4)') ihour
      write (name(11:12),'(i2.2)') imin
      name(14:21)= cmyid
      name(23:25)= cexpnr
      open  (ifoutput,file=trim(output_prefix)//name,form='unformatted',status='replace')

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
      write(ifoutput)  (((ekh   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
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
      write(ifoutput)  (  initial_presf (    k)                            ,k=1,k1)
      write(ifoutput)  (  initial_presh (    k)                            ,k=1,k1)
      write(ifoutput)  ps,thls,qts,thvs,oblav
      write(ifoutput)  dtheta,dqt,timee,  dt,tres
      write(ifoutput)   ((obl (i,j  ),i=1,i2      ),j=1,j2      )
      write(ifoutput)   ((tskin(i,j ),i=1,i2      ),j=1,j2      )
      write(ifoutput)   ((qskin(i,j ),i=1,i2      ),j=1,j2      )

!!!!! radiation quantities
      write(ifoutput)  tnext_radiation
      write(ifoutput)  (((thlprad (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      write(ifoutput)  (((swd     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      write(ifoutput)  (((swu     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      write(ifoutput)  (((lwd     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      write(ifoutput)  (((lwu     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      write(ifoutput)  (((swdca   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      write(ifoutput)  (((swuca   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      write(ifoutput)  (((lwdca   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      write(ifoutput)  (((lwuca   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      write(ifoutput)  (((swdir   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      write(ifoutput)  (((swdif   (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
      write(ifoutput)  (((lwc     (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)

      write(ifoutput)  ((SW_up_TOA    (i,j ),i=1,i2),j=1,j2)
      write(ifoutput)  ((SW_dn_TOA    (i,j ),i=1,i2),j=1,j2)
      write(ifoutput)  ((LW_up_TOA    (i,j ),i=1,i2),j=1,j2)
      write(ifoutput)  ((LW_dn_TOA    (i,j ),i=1,i2),j=1,j2)
      write(ifoutput)  ((SW_up_ca_TOA (i,j ),i=1,i2),j=1,j2)
      write(ifoutput)  ((SW_dn_ca_TOA (i,j ),i=1,i2),j=1,j2)
      write(ifoutput)  ((LW_up_ca_TOA (i,j ),i=1,i2),j=1,j2)
      write(ifoutput)  ((LW_dn_ca_TOA (i,j ),i=1,i2),j=1,j2)
!!!!! end of radiation quantities

      close (ifoutput)
      linkname = name
      linkname(6:13) = "_latest_"
      call system("ln -s -f "//name //" "//trim(output_prefix)//linkname)

      if (nsv>0) then
        name(5:5)='s'
        open  (ifoutput,file=trim(output_prefix)//name,form='unformatted')
        write(ifoutput) ((((sv0(i,j,k,n),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1),n=1,nsv)
        write(ifoutput) (((svflux(i,j,n),i=1,i2),j=1,j2),n=1,nsv)
        write(ifoutput) (dsv(n),n=1,nsv)
        write(ifoutput)  timee

        close (ifoutput)
        linkname = name
        linkname(6:13) = "_latest_"
        call system("ln -s -f "//name //" "//trim(output_prefix)//linkname)

      end if

      if (isurf == 1) then
        name(5:5)='l'
        open  (ifoutput,file=trim(output_prefix)//name,form='unformatted')
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
        linkname(6:13) = "_latest_"
        call system("ln -s -f "//name //" "//trim(output_prefix)//linkname)
      else if (isurf == 11) then
        name(5:5)='l'
        open  (ifoutput,file=trim(output_prefix)//name,form='unformatted')
        write(ifoutput) (((tsoil(i,j,k), i=1,i2), j=1,j2), k=1,kmax_soil)
        write(ifoutput) (((phiw (i,j,k), i=1,i2), j=1,j2), k=1,kmax_soil)
        write(ifoutput) ((tskin (i,j),   i=1,i2), j=1,j2)
        write(ifoutput) ((Wl    (i,j),   i=1,i2), j=1,j2)

        ! Sub-grid tiles
        do ilu=1,nlu
          write(ifoutput) ((tile(ilu)%thlskin(i,j), i=1,i2), j=1,j2)
          write(ifoutput) ((tile(ilu)%qtskin(i,j), i=1,i2), j=1,j2)
          write(ifoutput) ((tile(ilu)%obuk(i,j), i=1,i2), j=1,j2)
        end do

        write(ifoutput)  timee

        if (enable_slurb) then
          write(ifoutput) ((slurb_tile%t_can_0(i,j), i=1,i2), j=1,j2)
          write(ifoutput) ((slurb_tile%t_can_m(i,j), i=1,i2), j=1,j2)
          write(ifoutput) (((slurb_tile%t_wall_a_0(k,i,j), i=1,i2), j=1,j2), k=1,size(slurb_tile%t_wall_a_0,1))
          write(ifoutput) (((slurb_tile%t_wall_a_m(k,i,j), i=1,i2), j=1,j2), k=1,size(slurb_tile%t_wall_a_m,1))
          write(ifoutput) (((slurb_tile%t_wall_b_0(k,i,j), i=1,i2), j=1,j2), k=1,size(slurb_tile%t_wall_b_0,1))
          write(ifoutput) (((slurb_tile%t_wall_b_m(k,i,j), i=1,i2), j=1,j2), k=1,size(slurb_tile%t_wall_b_m,1))
          write(ifoutput) (((slurb_tile%t_win_a_0(k,i,j), i=1,i2), j=1,j2), k=1,size(slurb_tile%t_win_a_0,1))
          write(ifoutput) (((slurb_tile%t_win_a_m(k,i,j), i=1,i2), j=1,j2), k=1,size(slurb_tile%t_win_a_m,1))
          write(ifoutput) (((slurb_tile%t_win_b_0(k,i,j), i=1,i2), j=1,j2), k=1,size(slurb_tile%t_win_b_0,1))
          write(ifoutput) (((slurb_tile%t_win_b_m(k,i,j), i=1,i2), j=1,j2), k=1,size(slurb_tile%t_win_b_m,1))
          write(ifoutput) (((slurb_tile%t_roof_0(k,i,j), i=1,i2), j=1,j2), k=1,size(slurb_tile%t_roof_0,1))
          write(ifoutput) (((slurb_tile%t_roof_m(k,i,j), i=1,i2), j=1,j2), k=1,size(slurb_tile%t_roof_m,1))
          write(ifoutput) (((slurb_tile%t_road_0(k,i,j), i=1,i2), j=1,j2), k=1,size(slurb_tile%t_road_0,1))
          write(ifoutput) (((slurb_tile%t_road_m(k,i,j), i=1,i2), j=1,j2), k=1,size(slurb_tile%t_road_m,1))

          if (moist_physics) then
            write(ifoutput) ((slurb_tile%q_can_0(i,j), i=1,i2), j=1,j2)
            write(ifoutput) ((slurb_tile%q_can_m(i,j), i=1,i2), j=1,j2)
            write(ifoutput) ((slurb_tile%m_liq_roof_0(i,j), i=1,i2), j=1,j2)
            write(ifoutput) ((slurb_tile%m_liq_roof_m(i,j), i=1,i2), j=1,j2)
            write(ifoutput) ((slurb_tile%m_liq_road_0(i,j), i=1,i2), j=1,j2)
            write(ifoutput) ((slurb_tile%m_liq_road_m(i,j), i=1,i2), j=1,j2)
          end if
        end if

        close (ifoutput)
        linkname = name
        linkname(6:13) = "_latest_"
        call system("ln -s -f "//name //" "//trim(output_prefix)//linkname)
      end if

      if (myid==0) then
        write(*,'(A,F15.7,A,I4)') 'dump at time = ',rtimee,' unit = ',ifoutput
      end if

  end subroutine do_writerestartfiles

  subroutine testwctime
    use iso_fortran_env, only : real32
    use modmpi,    only : mpi_get_time
    use modglobal, only : timeleft
    implicit none
    real(real32), save :: tstart = -1., tend = -1. !MPI_get_time needs 4-byte real

    if (tstart < 0) then
      call mpi_get_time(tstart)
    else
      call mpi_get_time(tend)
      if (tend-tstart>=wctime) then
        write (*,*) wctime, "NO WALL CLOCK TIME LEFT"
        timeleft=0
      end if
    end if


  end subroutine testwctime

  subroutine exitmodules
    use modfields,         only : exitfields
    use modglobal,         only : exitglobal,lopenbc
    use modmpi,            only : exitmpi
    use modboundary,       only : exitboundary
    use modmicrophysics,   only : exitmicrophysics
    use modpois,           only : exitpois
    use modtimedep,        only : exittimedep
    use modradiation,      only : exitradiation
    use modsubgrid,        only : exitsubgrid
    use modtracers,        only : exittracers
    use modsurface,        only : exitsurface
    use modlsm,            only : exitlsm
    use modslurb,          only : exitslurb
    use moddrydeposition,  only : exitdrydep
    use modthermodynamics, only : exitthermodynamics
    use modemission,       only : exitemission
    use modopenboundary,   only : exitopenboundary
    use modibm,            only : exitibm
    use modchecksim,       only : exitchecksim
    use tstep,             only : exittstep

    call exittimedep
    call exitthermodynamics
    call exittracers
    call exittstep
    call exitchecksim
    call exitsurface
    call exitlsm
    call exitslurb
    call exitdrydep
    call exitsubgrid
    call exitradiation
    call exitpois
    call exitmicrophysics
    call exitemission
    if(lopenbc) then
      call exitopenboundary
    else
      call exitboundary
    endif
    call exitibm
    call exitfields
    call exitglobal
    call exitmpi

 end subroutine exitmodules
!----------------------------------------------------------------
  subroutine randomnize(field,klev,ampl,ir,ihl,jhl,negval)
    ! Adds (pseudo) random noise with given amplitude to the field at level k
    ! Use our own pseudo random function so results are reproducibly the same,
    ! independent of parallelization.

    use modmpi,    only : myidx, myidy
    use modglobal, only : itot,jtot,imax,jmax,i1,j1,k1
    integer (KIND=selected_int_kind(6)):: imm, ia, ic,ir
    integer ihl, jhl
    integer i,j,klev
    integer is,ie,js,je
    real ran,ampl
    real(field_r) field(2-ihl:i1+ihl,2-jhl:j1+jhl,k1)
    parameter (imm = 134456, ia = 8121, ic = 28411)
    logical negval

    is = myidx * imax + 1
    ie = is + imax - 1

    js = myidy * jmax + 1
    je = js + jmax - 1

    do j=1,jtot
    do i=1,itot
        ir=mod((ir)*ia+ic,imm)
        ran=real(ir)/real(imm)
        if (i >= is .and. i <= ie .and. &
            j >= js .and. j <= je) then
            if (.not. negval) then ! Avoid non-physical negative values
              field(i-is+2,j-js+2,klev) = field(i-is+2,j-js+2,klev) + (ran-0.5)*2.0*min(real(ampl, field_r),field(i-is+2,j-js+2,klev))
            else
              field(i-is+2,j-js+2,klev) = field(i-is+2,j-js+2,klev) + (ran-0.5)*2.0*ampl
            endif

        endif
    enddo
    enddo

    return
  end subroutine randomnize

  !> randomize fields using the system random number generator
  !  the random perturbation is reproducible with the same seed only with the
  !  same parallelization.
  subroutine randomize_new (irandom)
    use modglobal,  only : i1,j1,k1,ih,jh
    use modfields,  only : u0,v0,w0,um,vm,wm,thlm,thl0,qtm,qt0
    use modmpi,     only : myid

    integer, intent(in) :: irandom
    integer :: n
    integer, allocatable :: seed(:)
    real(field_r), allocatable :: noise(:,:,:)
    character(len=*), parameter :: routine = modname//'/randomize_new'

    call random_seed(size = n) ! query the seed size
    if (n < 2) then
       call finish(routine, "The random seed size on this system is too small", n)
    end if
    allocate(seed(n))
    ! initialize the seed array with the random seed from the namelist and MPI rank
    ! -> reproducible random numbers if parallelization is not changed
    seed = 0
    seed(1) = irandom
    seed(2) = myid
    call random_seed(put=seed)
    deallocate(seed)

    allocate(noise(2-ih:i1+ih,2-jh:j1+jh,krand))

    call random_number(noise)
    qtm(2:i1,2:j1,1:krand)  = qtm(2:i1,2:j1,1:krand)  + &
         min(randqt,qtm(2:i1,2:j1,1:krand)) * 2 * (noise(2:i1,2:j1,1:krand)-0.5_field_r)  ! avoid negative q
    call random_number(noise)                           
    qt0(2:i1,2:j1,1:krand)  = qt0(2:i1,2:j1,1:krand)  + &
         min(randqt,qt0(2:i1,2:j1,1:krand)) * 2 * (noise(2:i1,2:j1,1:krand)-0.5_field_r)
    call random_number(noise)
    thlm(2:i1,2:j1,1:krand) = thlm(2:i1,2:j1,1:krand) + 2*randthl * (noise(2:i1,2:j1,1:krand)-0.5_field_r)
    call random_number(noise)
    thl0(2:i1,2:j1,1:krand) = thl0(2:i1,2:j1,1:krand) + 2*randthl * (noise(2:i1,2:j1,1:krand)-0.5_field_r)

    ! momentum is by default not ranomized with the old system
    ! disabling for now because it seems to lead to bad divergence
    !call random_number(noise)
    !um(2:i1,2:j1,1:krand)   = um(2:i1,2:j1,1:krand)   + 2*randu   * (noise(2:i1,2:j1,1:krand)-0.5_field_r)
    !call random_number(noise)
    !u0(2:i1,2:j1,1:krand)   = u0(2:i1,2:j1,1:krand)   + 2*randu   * (noise(2:i1,2:j1,1:krand)-0.5_field_r)
    !call random_number(noise)
    !vm(2:i1,2:j1,1:krand)   = vm(2:i1,2:j1,1:krand)   + 2*randu   * (noise(2:i1,2:j1,1:krand)-0.5_field_r)
    !call random_number(noise)
    !v0(2:i1,2:j1,1:krand)   = v0(2:i1,2:j1,1:krand)   + 2*randu   * (noise(2:i1,2:j1,1:krand)-0.5_field_r)
    !call random_number(noise)
    !wm(2:i1,2:j1,1:krand)   = wm(2:i1,2:j1,1:krand)   + 2*randu   * (noise(2:i1,2:j1,1:krand)-0.5_field_r)
    !call random_number(noise)
    !w0(2:i1,2:j1,1:krand)   = w0(2:i1,2:j1,1:krand)   + 2*randu   * (noise(2:i1,2:j1,1:krand)-0.5_field_r)

    deallocate(noise)
  end subroutine randomize_new

  subroutine baseprofs
    ! Calculates the profiles corresponding to the base state
    ! In the current implementation, neither the base pressure, nor the base virtual temperature plays a role in the dynamics
    ! They are nevertheless calculated and printed to the stdin/baseprof files for user convenience
    use modfields,         only : rhobf,rhobh,exnf,exnh
    use modglobal,         only : k1,kmax,zf,zh,dzf,dzh,rv,rd,grav,cp,pref0,lwarmstart,ibas_prf,cexpnr,ifinput,ifoutput, ibas_usr
    use modthermodynamics, only : lbaseexner
    use modsurfdata,       only : thls,ps,qts
    use modmpi,            only : myid,comm3d,mpierr,D_MPI_BCAST
    use modlogging,        only : profile_output
    implicit none

    character(len=*), parameter :: routine = modname//'/baseprofs'

    real :: thvb,prsb ! for calculating moist adiabat
    integer :: j,k
    real(field_r), allocatable :: height(:),pb(:),tb(:),pbh(:)
    character(80) chmess
    real :: zsurf=0.
    real :: tsurf
    real(field_r),dimension(4) :: zmat=(/11000.,20000.,32000.,47000./)
    real(field_r),dimension(4) :: lapserate=(/-6.5/1000.,0.,1./1000,2.8/1000/)
    real(field_r),dimension(4) :: pmat
    real(field_r),dimension(4) :: tmat

    allocate (height(k1),pb(k1),tb(k1),pbh(k1))

    if(myid==0)then

      if( (.not. lwarmstart) .and. (ibas_prf /= ibas_usr) ) then

        if(ibas_prf <= 3 .and. thls < 0) then
          call finish(routine, 'thls has not been initialized but is needed for setting up the base profiles.')
        end if

        if(ibas_prf==1) then !thv constant and hydrostatic balance
          thvb=thls*(1+(rv/rd-1)*qts) ! using thls, q_l assumed to be 0 during first time step
          do k=1,k1
            prsb=(ps**(rd/cp)-(grav*zf(k)*pref0**(rd/cp))/(cp*thvb))**(cp/rd) !As in thermodynamics
            rhobf(k)=prsb/(rd*thvb*((prsb/pref0)**(rd/cp)))
          end do
        else if(ibas_prf==2) then ! Quasi-Boussinesq (Similar to Dales 3, except for buoyancy term now depending on slab mean state)
          thvb=thls*(1+(rv/rd-1)*qts)
          rhobh(1)=ps/(rd*thvb*(ps/pref0)**(rd/cp))
          do k=1,k1
            rhobf(k)=rhobh(1)
          end do
        else if(ibas_prf==3) then! use standard atmospheric lapse rate with surface temperature offset
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
              do while(zf(k)>=zmat(j))
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
          end do
        else if(ibas_prf==4) then! use standard atmospheric lapse rate without surface temperature offset
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
            end if
            tmat(j)=tmat(j-1)+lapserate(j)*(zmat(j)-zmat(j-1));
          end do

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
              end if
            end if
            rhobf(k)=pb(k)/(rd*tb(k)) ! dry estimate
          end do
        end if

        ! Write background profiles in all cases
        open (ifoutput,file='baseprof.inp.'//cexpnr)
        write(ifoutput,*) '#baseprofiles'
        write(ifoutput,*) '#height rhobf'
        do k=1,kmax
          write (ifoutput,'(1f7.1,E25.17)') &
                zf (k), &
                rhobf (k)
        end do
        close(ifoutput)

      else ! lwarmstart or user-specified baseprof

        if (lwarmstart) then
          ibas_prf = ibas_usr
          print *, 'WARNING: warm start requires input files for density. ibas_prf defaulted to 5'
        end if

        ! Read background profiles in case of warmstart or user-provided input
        open (ifinput,file='baseprof.inp.'//cexpnr)
        read (ifinput,'(a80)') chmess
        read (ifinput,'(a80)') chmess

        do k = 1, kmax
          read (ifinput,*) &
                  height(k), &
                  rhobf (k)
        end do
        close(ifinput)

      end if ! end if .not. lwarmstart .and. .not.user-specified baseprof

      ! Set height at k1 equal to kmax for the sake of printing to screen
      height(k1) = height(kmax)

      rhobf(k1)=rhobf(kmax)+(zf(k1)-zf(kmax))/(zf(kmax)-zf(kmax-1))*(rhobf(kmax)-rhobf(kmax-1))
      do k = 2, k1
        rhobh(k) = (rhobf(k)*dzf(k-1)+rhobf(k-1)*dzf(k))/(dzf(k)+dzf(k-1))
      end do
      rhobh(1) = rhobf(1)-(rhobf(2)-rhobf(1))*(zf(1)-zh(1))/(zf(2)-zf(1))

      ! pb is only available for ibas_prf >= 3
      if (ibas_prf >= 3) then
        pbh(1)   = ps
        do k = 2, k1
          pbh(k)   = (   pb(k)*dzf(k-1)+   pb(k-1)*dzf(k))/(dzf(k)+dzf(k-1)) ! interpolate base half-level pressure like half-level base rho
        end do 
      end if

      ! write profiles and derivatives to standard output
      write (profile_output,*) ' height   rhobf       rhobh'
      do k=k1,1,-1
          write (profile_output,'(1f7.1,2E25.17)') &
                height (k), &
                rhobf (k), &
                rhobh (k)
      end do

      ! exner function from base profiles
      ! TODO: pb is not available here on warm start
      if (lbaseexner) then
         exnf = (pb/pref0)**(rd/cp)
         exnh = (pbh/pref0)**(rd/cp)
      end if

    end if ! ENDIF MYID=0

    ! MPI broadcast variables
    call D_MPI_BCAST(rhobf       ,k1,0,comm3d,mpierr)
    call D_MPI_BCAST(rhobh       ,k1,0,comm3d,mpierr)

    if (lbaseexner) then
       call D_MPI_BCAST(exnf        ,k1,0,comm3d,mpierr)
       call D_MPI_BCAST(exnh        ,k1,0,comm3d,mpierr)
    end if

    deallocate(height,pb,tb,pbh)

  end subroutine baseprofs

  !> \brief Read initial profiles from init.XXX.nc
  !!
  !! \param filename Path to the netCDF file to read from.
  !! \param height Vertical levels.
  !! \param uprof Initial eastward velocity profile.
  !! \param vprof Initial northward velocity profile.
  !! \param thlprof Initial liquid water potential temperature profile.
  !! \param qtprof Initial total water mixing ratio profile.
  !! \param e12prof Initial profile of the square root of the turbulence kinetic energy (TKE).
  !! \param ug Geostrophic eastward wind.
  !! \param vg Geostrophic northward wind.
  !! \param wfls Large-scale subsidence.
  !! \param dqtdxls Eastward gradient of the total water mixing ratio due to advection.
  !! \param dqtdyls Northward gradient of the total water mixing ratio due to advection.
  !! \param dqtdtls Tendency of the total water mixing ratio.
  !! \param dthlrad Tendency of the liquid water potential temperature due to radiative heating.
  !! \param kmax Index of highest vertical level.
  !!
  !! \note Tracers are read from tracers.XXX.nc, not here.
  !! \todo Make DEPHY-compatible.
  subroutine init_from_netcdf(filename, height, uprof, vprof, thlprof, qtprof, &
                              e12prof, ug, vg, dpdxl, dpdyl, wfls, dqtdxls, dqtdyls, &
                              dqtdtls, dthlrad, kmax)
    character(*),   intent(in)  :: filename
    real(field_r),  intent(out) :: height(:)
    real(field_r),  intent(out) :: uprof(:)
    real(field_r),  intent(out) :: vprof(:)
    real(field_r),  intent(out) :: thlprof(:)
    real(field_r),  intent(out) :: qtprof(:)
    real(field_r),  intent(out) :: e12prof(:)
    real(field_r),  intent(out) :: ug(:)
    real(field_r),  intent(out) :: vg(:)
    real(field_r),  intent(out) :: dpdxl(:)
    real(field_r),  intent(out) :: dpdyl(:)
    real(field_r),  intent(out) :: wfls(:)
    real(field_r),  intent(out) :: dqtdxls(:)
    real(field_r),  intent(out) :: dqtdyls(:)
    real(field_r),  intent(out) :: dqtdtls(:)
    real(field_r),  intent(out) :: dthlrad(:)
    integer,        intent(in)  :: kmax

    integer :: ncid

    call nchandle_error(nf90_open(filename, NF90_NOWRITE, ncid))

    ! "Regular" prognostic fields
    call read_nc_field(ncid, "ua", uprof, start=1, count=kmax, &
                       fillvalue=0._field_r)
    call read_nc_field(ncid, "va", vprof, start=1, count=kmax, &
                       fillvalue=0._field_r)
    call read_nc_field(ncid, "thetal", thlprof, start=1, count=kmax, &
                       fillvalue=0._field_r)
    call read_nc_field(ncid, "qt", qtprof, start=1, count=kmax, &
                       fillvalue=0._field_r)
    call read_nc_field(ncid, "tke", e12prof, start=1, count=kmax, &
                       fillvalue=0._field_r)
    ! reading with no count in assumes kmax+1 values in the NC file, so we need to set count
    call read_nc_field(ncid, "zh", height, start=1, count=kmax)

    ! Large-scale forcings
    call read_nc_field(ncid, "ug", ug, start=1, count=kmax, &
                       fillvalue=0._field_r)
    call read_nc_field(ncid, "vg", vg, start=1, count=kmax, &
                       fillvalue=0._field_r)
    call read_nc_field(ncid, "dpdx", dpdxl, start=1, count=kmax, &
                       fillvalue=0._field_r)
    call read_nc_field(ncid, "dpdy", dpdyl, start=1, count=kmax, &
                       fillvalue=0._field_r)
    call read_nc_field(ncid, "wa", wfls, start=1, count=kmax, &
                       fillvalue=0._field_r)
    call read_nc_field(ncid, "dqtdxls", dqtdxls, start=1, count=kmax, &
                       fillvalue=0._field_r)
    call read_nc_field(ncid, "dqtdyls", dqtdyls, start=1, count=kmax, &
                       fillvalue=0._field_r)
    call read_nc_field(ncid, "tnqt_adv", dqtdtls, start=1, count=kmax, &
                       fillvalue=0._field_r)
    call read_nc_field(ncid, "tnthetal_rad", dthlrad, start=1, count=kmax, &
                       fillvalue=0._field_r)

    call nchandle_error(nf90_close(ncid))

  end subroutine init_from_netcdf

  !> Check prognostic variables before simulation
  subroutine check_initial_state()
    use modglobal, only: lopenbc, i1, j1, kmax
    use modthermodynamics, only: lmoist
    use modfields, only: u0, v0, w0, thl0, qt0, sv0
    use modchecksim, only: lstop

    ! Weird bug, casting thresholds to _field_r leads to compilation error for
    ! some reason
    integer, parameter :: rkind = kind(u0)

    integer :: s

    call check_array(u0, 'u0', 'startup', &
                     threshold=[real(-100, rkind), real(100, rkind)],stop_if_invalid=lstop, dump_if_invalid=.true.)
    call check_array(v0, 'v0', 'startup', &
                     threshold=[real(-100, rkind), real(100, rkind)],stop_if_invalid=lstop, dump_if_invalid=.true.)
    call check_array(w0, 'w0', 'startup', &
                     threshold=[real(-30, rkind), real(30, rkind)],stop_if_invalid=lstop, dump_if_invalid=.true.)
    if (lopenbc) then
      call check_array(thl0(2:i1,2:j1,1:kmax), 'thl0', 'startup', &
                      threshold=[real(150, rkind), real(2000, rkind)],stop_if_invalid=lstop, dump_if_invalid=.true.)
    else
      call check_array(thl0, 'thl0', 'startup', &
                      threshold=[real(150, rkind), real(2000, rkind)],stop_if_invalid=lstop, dump_if_invalid=.true.)
    end if
    if (lmoist) call check_array(qt0, 'qt0', 'startup', &
                                 threshold=[real(0, rkind), real(1, rkind)],stop_if_invalid=lstop, dump_if_invalid=.true.)

    do s = 1, size(sv0, dim=4)
      call check_array(sv0(:,:,:,s), 'sv0('//int2string(s)//')', 'startup',stop_if_invalid=lstop, dump_if_invalid=.true.)
    end do

  end subroutine check_initial_state

end module modstartup
