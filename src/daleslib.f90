!> \file daleslib.f90
!! Library API for Dales.
!  This file is part of DALESLIB.
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

! Initializes the libary's dales instance. This routine requires a path to the namoptions file to be passed 
! as an argument, and optionally the mpi communicator can be set too (this is particularly useful when 
! nesting into other parallelized codes). If the latter argument is not provided, the communicator will be 
! assumed to be MPI_COMM_WORLD. 

module daleslib
    !use modglobal, only: ifmessages
    
    implicit none

    integer, parameter :: FIELDID_U=1
    integer, parameter :: FIELDID_V=2
    integer, parameter :: FIELDID_W=3
    integer, parameter :: FIELDID_THL=4
    integer, parameter :: FIELDID_QT=5

    integer, parameter :: QT_FORCING_GLOBAL=0
    integer, parameter :: QT_FORCING_LOCAL=1
    integer, parameter :: QT_FORCING_VARIANCE=2
    integer, parameter :: QT_FORCING_STRONG=3
   
    integer :: my_task, master_task

    real, allocatable :: u_tend(:)
    real, allocatable :: v_tend(:)
    real, allocatable :: thl_tend(:)
    real, allocatable :: qt_tend(:)
    real, allocatable :: ql_tend(:)  ! Only used whith local nudging
    real, allocatable :: ql_ref(:)   ! QL profile from global model - not a tendency
    real, allocatable :: qt_alpha(:)

    real, allocatable :: u_nudge(:)
    real, allocatable :: v_nudge(:)
    real, allocatable :: thl_nudge(:)
    real, allocatable :: qt_nudge(:)
    real :: u_nudge_time = 0
    real :: v_nudge_time = 0
    real :: thl_nudge_time = 0
    real :: qt_nudge_time = 0

    
    real :: ps_tend

    integer :: qt_forcing_type = QT_FORCING_GLOBAL

    
    contains

        subroutine initialize(path,mpi_communicator,date,time)

            !!----------------------------------------------------------------
            !!     0.0    USE STATEMENTS FOR CORE MODULES
            !!----------------------------------------------------------------
            use modmpi,             only : initmpicomm,myid,MPI_COMM
            use modstartup,         only : startup

            !----------------------------------------------------------------
            !     0.1     USE STATEMENTS FOR ADDONS STATISTICAL ROUTINES
            !----------------------------------------------------------------
            use modcape,            only : initcape
            use modchecksim,        only : initchecksim
            use modstat_nc,         only : initstat_nc
            !use modspectra2,       only : initspectra2
            use modtimestat,        only : inittimestat
            use modgenstat,         only : initgenstat
            use modradstat,         only : initradstat
            use modlsmstat,         only : initlsmstat
            use modsampling,        only : initsampling
            use modquadrant,        only : initquadrant
            use modcrosssection,    only : initcrosssection 
            use modAGScross,        only : initAGScross
            use modlsmcrosssection, only : initlsmcrosssection
            use modcloudfield,      only : initcloudfield
            use modfielddump,       only : initfielddump
            use modsamptend,        only : initsamptend

            use modbulkmicrostat,   only : initbulkmicrostat
            use modbudget,          only : initbudget
            use modheterostats,     only : initheterostats

            ! modules below are disabled by default to improve compilation time
            !use modstress,         only : initstressbudget

            !use modtilt,           only : inittilt
            !use modparticles,      only : initparticles
            use modnudge,           only : initnudge
            !use modprojection,     only : initprojection
            use modchem,            only : initchem
            use modcanopy,          only : initcanopy
            implicit none

            character(len=256), intent(in)  :: path
            type(MPI_COMM), intent(in), optional :: mpi_communicator
            integer, intent(in), optional   :: date,time
            integer                         :: yy,mo,dd,hh,mm,ss

            !----------------------------------------------------------------
            !     0      INITIALIZE MPI COMMUNICATOR
            !----------------------------------------------------------------

            call initmpicomm(mpi_communicator)

            !----------------------------------------------------------------
            !     1      READ NAMELISTS,INITIALISE GRID, CONSTANTS AND FIELDS
            !----------------------------------------------------------------

            call startup(path)

            ! Initial time overrides

            if (present(date)) then

               if (date /= 0) then
                  yy = date/10000
                  mo = mod(date,10000)/100
                  dd = mod(date,100)
                  hh = 0
                  mm = 0
                  ss = 0
                  if (present(time)) then
                     hh = time/10000
                     mm = mod(time,10000)/100
                     ss = mod(time,100)
                  endif
                  call set_start_time(yy,mo,dd,hh,mm,ss)
               endif
            endif

            !---------------------------------------------------------
            !      2     INITIALIZE STATISTICAL ROUTINES AND ADD-ONS
            !---------------------------------------------------------
            call initchecksim
            call initstat_nc   ! Should be called before stat-routines that might do netCDF
            call inittimestat  ! Timestat must preceed all other timeseries that could write in the same netCDF file (unless stated otherwise
            call initgenstat   ! Genstat must preceed all other statistics that could write in the same netCDF file (unless stated otherwise
            !call inittilt
            call initsampling
            call initquadrant
            call initcrosssection
            call initAGScross
            call initlsmcrosssection
            !call initprojection
            call initcloudfield
            call initfielddump
            call initsamptend
            call initradstat
            call initlsmstat
            !call initparticles
            call initnudge
            call initbulkmicrostat
            call initbudget
            !call initstressbudget
            call initchem
            call initheterostats
            call initcanopy

            !call initspectra2
            call initcape

            !Set additional library information
            my_task=myid
            master_task=0

            call initdaleslib
            
        end subroutine initialize


        ! allocate arrays for tendencies set through the daleslib interface
        ! NOTE don't initialize parameters here - this is called after
        ! parameters have been set
        subroutine initdaleslib
          use modglobal, only: kmax
          ! use modforces, only: lforce_user
          implicit none

          allocate(u_tend(1:kmax))
          allocate(v_tend(1:kmax))
          allocate(thl_tend(1:kmax))
          allocate(qt_tend(1:kmax))
          allocate(ql_tend(1:kmax))
          allocate(ql_ref(1:kmax))
          allocate(qt_alpha(1:kmax))
          u_tend = 0
          v_tend = 0
          thl_tend = 0
          qt_tend = 0
          ql_tend = 0
          ql_ref = 0
          ! lforce_user = .true.
          ps_tend = 0
          qt_alpha = 0

          allocate(u_nudge(1:kmax))
          allocate(v_nudge(1:kmax))
          allocate(thl_nudge(1:kmax))
          allocate(qt_nudge(1:kmax))

          u_nudge = 0
          v_nudge = 0
          thl_nudge = 0
          qt_nudge = 0


        end subroutine initdaleslib

                ! deallocate arrays for tendencies 
        subroutine exitdaleslib
          implicit none

          deallocate(u_tend)
          deallocate(v_tend)
          deallocate(thl_tend)
          deallocate(qt_tend)
          deallocate(ql_tend)
          deallocate(ql_ref)
          deallocate(qt_alpha)
          deallocate(u_nudge)
          deallocate(v_nudge)
          deallocate(thl_nudge)
          deallocate(qt_nudge)
        end subroutine exitdaleslib


        subroutine daleslib_nudge
           use modfields,   only : up,vp,thlp,qtp, u0av,v0av,thl0av,qt0av
           use modglobal,   only : i1,j1,kmax
           integer k
           do k=1,kmax
              if (u_nudge_time /= 0)     up(2:i1,2:j1,k) =   up(2:i1,2:j1,k) + (u_nudge(k) - u0av(k))     / u_nudge_time
              if (v_nudge_time /= 0)     vp(2:i1,2:j1,k) =   vp(2:i1,2:j1,k) + (v_nudge(k) - v0av(k))     / v_nudge_time
              if (thl_nudge_time /= 0) thlp(2:i1,2:j1,k) = thlp(2:i1,2:j1,k) + (thl_nudge(k) - thl0av(k)) / thl_nudge_time
              if (qt_nudge_time /= 0)   qtp(2:i1,2:j1,k) =  qtp(2:i1,2:j1,k) + (qt_nudge(k) - qt0av(k))   / qt_nudge_time
           end do
           !write(*,*) 'nudge(): qt_nudge_time=', qt_nudge_time
           !if (qt_nudge_time /= 0) write(*,*) (qt_nudge - qt0av)   / qt_nudge_time
           
         end subroutine daleslib_nudge
          
        
        subroutine force_tendencies
          use modglobal,   only : i1,j1,imax,jmax,kmax,rdt
          use modfields,   only : up,vp,thlp,qtp,qt0,ql0,ql0av,qt0av         

          implicit none
          integer k
          real alpha, qlt, qvt
          real qtp_local (2:i1, 2:j1), qtp_local_lim (2:i1, 2:j1), qtp_lost, al(1:kmax)

          !write(*,*) "force_tendencies() : qt_forcing_type =", qt_forcing_type
          
          if (qt_forcing_type == QT_FORCING_LOCAL) then
             al = 0
             if (gathersatfrac(al) /= 0) then
                 write(*,*) "MPI error in force_tendencies!"
                 return
             endif
          endif

          do k=1,kmax
             up  (2:i1,2:j1,k) = up  (2:i1,2:j1,k) + u_tend(k) 
             vp  (2:i1,2:j1,k) = vp  (2:i1,2:j1,k) + v_tend(k) 
             thlp(2:i1,2:j1,k) = thlp(2:i1,2:j1,k) + thl_tend(k)


             qtp_lost = 0

             if (qt_forcing_type == QT_FORCING_GLOBAL) then
                qtp_local =  qt_tend(k)
             end if

             if (qt_forcing_type == QT_FORCING_VARIANCE) then
                ! force fluctuations of qt, with a forcing profile given by alpha(k)
                qtp_local = (qt0(2:i1,2:j1,k) - qt0av(k)) * qt_alpha(k) + qt_tend(k)
             endif
             
             if (qt_forcing_type == QT_FORCING_STRONG) then
                ! adjust qt variability to make ql match ql_ref from OpenIFS.
                alpha = 0
                ! not enough clouds on this level.
                ! OpenIFS can have very small but non-zero QL. We ignore these.
                 if (ql_ref(k) > ql0av(k) .and. ql_ref(k) > 1e-6) then
                    alpha = 0.01
                 end if
                 
                 ! we have too much clouds on this level
                 if (ql_ref(k) < ql0av(k) * .95) then
                    alpha = -0.01
                 end if
                 ! write (*,*) "strong nudge", k, ql0av(k), ql_ref(k), alpha
                 qtp_local = (qt0(2:i1,2:j1,k) - qt0av(k)) * alpha + qt_tend(k)
             endif
             
             if (qt_forcing_type == QT_FORCING_LOCAL) then
                qlt = ql_tend(k)

                if (1 - al(k) /= 0) then  ! avoid divide-by-0
                   qvt = (qt_tend(k) + al(k)*qlt)/(1 - al(k))
                end if
                qtp_local = merge(qlt,qvt,ql0(2:i1,2:j1,k) > 0)  ! (ql0>0) ? qlt : qvt
             endif                

             ! limit the tendency, to avoid qt < 0
             qtp_local_lim = max (qtp_local, -qt0(2:i1, 2:j1, k)/(rdt*1.2) )  ! "rdt" here used to be 60. *1.2 for some safety margin.

             qtp_lost = sum(qtp_local - qtp_local_lim)
             ! < 0 if the cut-off was activated
             ! NOTE !  the correction is per thread for simplicity
             
             qtp(2:i1,2:j1,k) = qtp(2:i1,2:j1,k)  +  qtp_local_lim  + ( qtp_lost / (imax * jmax) ) 



          enddo
        end subroutine force_tendencies


!             if (qt_correction_type == QT_CORRECTION_REGULAR) then
!                ! additive correction of qt
!                qtp(2:i1,2:j1,k) = qtp(2:i1,2:j1,k) + ( qt_tend(k) +  qtp_lost / (imax * jmax) )
!             endif
        
             ! multiplicative correcion of qt
!             if (qt_correction_type == QT_CORRECTION_RESCALE) then
!                qtp(2:i1,2:j1,k) = qtp(2:i1,2:j1,k)  +  qt0(2:i1,2:j1,k) / qt0av(k) * ( qt_tend(k) +  qtp_lost / (imax * jmax) )
!             endif

             !    ! do this if we don't have enough clouds on this k-level. OpenIFS can have very low but non-zero QL, 
             !    ! we only care if QL large enough. 1e-4 is too high - never happens.
             !    alpha = 0
             !    if (ql_ref(k) > ql0av(k) .and. ql_ref(k) > 1e-6) then
             !       alpha = 0.01
             !    end if

             !    ! we have too much clouds on this level
             !    if (ql_ref(k) < ql0av(k) * .95) then
             !       alpha = -0.01
             !    end if

             !    if (alpha /= 0) then
             !       ! amplify fluctuations from the mean by factor (1+alpha)
             !       qtp_local = (qt0(2:i1,2:j1,k) - qt0av(k)) * alpha
                   
             !       ! limit the tendency from below, so that it cannot make qt0 negative in a 60s time step
             !       qtp_local_lim = max (qtp_local, -qt0(2:i1, 2:j1, k)/60.0)
             !       qtp_lost = sum(qtp_local - qtp_local_lim) ! < 0 if the cut-off was activated
                   
             !       ! NOTE !  the correction is per thread for simplicity
                   
             !       qtp(2:i1,2:j1,k) = qtp(2:i1,2:j1,k)  +  qtp_local_lim 
                   
             !       write(*,*) k, 'ql_ref', ql_ref(k), ' <-> ', 'ql0av', ql0av(k), 'adjusting qt fluctuations by', alpha, '. qtp_lost:', qtp_lost
             !    endif




        
        ! find min and max value of all tendencies. Check that they are within reasonable boundaries.
        ! if not, print a message including 'msg' and the min/max values of all the tendencies.
        subroutine check_tend(msg)
          use modfields,   only : up,vp,wp,thlp,qtp,e12p
          use modglobal,   only : i1,j1,kmax,rk3step,rtimee
          implicit none
          real thlp_l, qtp_l, e12p_l, up_l, vp_l, wp_l, thlp_h, qtp_h, e12p_h, up_h, vp_h, wp_h
          character (*), intent (in) :: msg
          thlp_l = minval(thlp(2:i1,2:j1,1:kmax))
          thlp_h = maxval(thlp(2:i1,2:j1,1:kmax))
          qtp_l  = minval(qtp(2:i1,2:j1,1:kmax))
          qtp_h  = maxval(qtp(2:i1,2:j1,1:kmax))
          e12p_l = minval(e12p(2:i1,2:j1,1:kmax))
          e12p_h = maxval(e12p(2:i1,2:j1,1:kmax))
          up_l    = minval(up(2:i1,2:j1,1:kmax))
          up_h    = maxval(up(2:i1,2:j1,1:kmax))
          vp_l    = minval(vp(2:i1,2:j1,1:kmax))
          vp_h    = maxval(vp(2:i1,2:j1,1:kmax))
          wp_l    = minval(wp(2:i1,2:j1,1:kmax))
          wp_h    = maxval(wp(2:i1,2:j1,1:kmax))
          if (max(-up_l, up_h, -vp_l, vp_h, -wp_l, wp_h) > 10 .or. max(-qtp_l, qtp_h) > 0.003 .or. max(-thlp_l, thlp_h) > 5 .or. max(-e12p_l, e12p_h) > 1) then
             write(*,*) '--- extreme tendencies found ', msg, ' ---'
             write(*,*) 'thlp', thlp_l, thlp_h
             write(*,*) 'qtp', qtp_l, qtp_h
             write(*,*) 'e12p', e12p_l, e12p_h
             write(*,*) 'up', up_l, up_h
             write(*,*) 'vp', vp_l, vp_h
             write(*,*) 'wp', wp_l, wp_h
             write(*,*) 'rk3step:', rk3step, 'rtimee:', rtimee
             write(*,*) '-----------------'
          endif
          
        end subroutine check_tend

        
        ! Performs a single time step. This is exactly the code inside the time loop of the main program. 
        ! If the adaptive time stepping is enabled, the new time stamp is not guaranteed to be a dt later.
        subroutine step
            !!----------------------------------------------------------------
            !!     0.0    USE STATEMENTS FOR CORE MODULES
            !!----------------------------------------------------------------
            use modstartup,         only : writerestartfiles
            use modtimedep,         only : timedep
            use modboundary,        only : boundary, grwdamp! JvdD ,tqaver
            use modthermodynamics,  only : thermodynamics
            use modmicrophysics,    only : microsources
            use modsurface,         only : surface
            use modsubgrid,         only : subgrid
            use modforces,          only : forces, coriolis, lstend
            use modradiation,       only : radiation
            use modpois,            only : poisson
            !use modedgecold,       only : coldedge
            use tstep,             only : tstep_update,  tstep_integrate

            !----------------------------------------------------------------
            !     0.1     USE STATEMENTS FOR ADDONS STATISTICAL ROUTINES
            !----------------------------------------------------------------
            use modcape,            only : docape
            use modchecksim,        only : checksim
            !use modspectra2,       only : dospecs, tanhfilter
            use modtimestat,        only : timestat
            use modgenstat,         only : genstat
            use modradstat,         only : radstat
            use modlsmstat,         only : lsmstat
            use modsampling,        only : sampling
            use modquadrant,        only : quadrant
            use modcrosssection,    only : crosssection
            use modAGScross,        only : AGScross
            use modlsmcrosssection, only : lsmcrosssection
            use modcloudfield,      only : cloudfield
            use modfielddump,       only : fielddump
            use modsamptend,        only : samptend,tend_start,tend_adv,tend_subg,tend_force,&
                tend_rad,tend_ls,tend_micro,tend_topbound,tend_pois,tend_addon,tend_coriolis,leibniztend

            use modbulkmicrostat,   only : bulkmicrostat
            use modbudget,          only : budgetstat
            use modheterostats,     only : heterostats

            ! modules below are disabled by default to improve compilation time
            !use modstress,         only : stressbudgetstat

            !use modtilt,           only : tiltedgravity, tiltedboundary
            !use modparticles,      only : particles
            use modnudge,           only : nudge
            !use modprojection,     only : projection
            use modchem,            only : twostep
            use modcanopy,          only : canopy
            use modadvection,       only : advection

            implicit none


                      
            call tstep_update                           ! Calculate new timestep
            call timedep
            call update_ps
            
            call samptend(tend_start,firstterm=.true.)

            !-----------------------------------------------------
            !   3.1   RADIATION         
            !-----------------------------------------------------
            call radiation !radiation scheme
            call samptend(tend_rad)

            !-----------------------------------------------------
            !   3.2   THE SURFACE LAYER
            !-----------------------------------------------------
            call surface

            !-----------------------------------------------------
            !   3.3   ADVECTION AND DIFFUSION
            !-----------------------------------------------------
            call advection
            call samptend(tend_adv)
            !call check_tend('after advection')
            call subgrid
            ! call check_tend('after subgrid')
            call canopy
            call samptend(tend_subg)

            !-----------------------------------------------------
            !   3.4   REMAINING TERMS
            !-----------------------------------------------------
            call coriolis !remaining terms of ns equation
            call samptend(tend_coriolis)
            call forces !remaining terms of ns equation
            call force_tendencies         ! NOTE - not standard DALES, these are our own tendencies
            call daleslib_nudge           ! NOTE - not standard DALES, omuse-interface for nudging towards given profiles
            call samptend(tend_force)

            call lstend !large scale forcings
            call samptend(tend_ls)

            !call check_tend('after lstend')
            call microsources !Drizzle etc.
            !call check_tend('after microsources')
            call samptend(tend_micro)

            !------------------------------------------------------
            !   3.4   EXECUTE ADD ONS
            !------------------------------------------------------
            call nudge
            !    call dospecs
            !    call tiltedgravity

            call samptend(tend_addon)

            !-----------------------------------------------------------------------
            !   3.5  PRESSURE FLUCTUATIONS, TIME INTEGRATION AND BOUNDARY CONDITIONS
            !-----------------------------------------------------------------------
            call grwdamp !damping at top of the model
            !JvdD    call tqaver !set thl, qt and sv(n) equal to slab average at level kmax
            call samptend(tend_topbound)
            call poisson
            call check_tend('after poisson')
            call samptend(tend_pois,lastterm=.true.)

            call tstep_integrate                        ! Apply tendencies to all variables
            call boundary
            !call tiltedboundary
            !-----------------------------------------------------
            !   3.6   LIQUID WATER CONTENT AND DIAGNOSTIC FIELDS
            !-----------------------------------------------------
            call thermodynamics
            call leibniztend
            !-----------------------------------------------------
            !   3.7  WRITE RESTARTFILES AND DO STATISTICS
            !------------------------------------------------------
            call twostep
            !call coldedge
            call checksim
            call timestat  !Timestat must preceed all other timeseries that could write in the same netCDF file (unless stated otherwise
            call genstat  !Genstat must preceed all other statistics that could write in the same netCDF file (unless stated otherwise
            call radstat
            call lsmstat
            call sampling
            call quadrant
            call crosssection
            call AGScross
            call lsmcrosssection
            !call tanhfilter
            call docape
            !call projection
            call cloudfield
            call fielddump
            !call particles

            call bulkmicrostat
            call budgetstat
            !call stressbudgetstat
            call heterostats

            call writerestartfiles


        end subroutine step

        ! update the surface pressure according to ps_tend
        ! perform the update every full timestep to avoid introducing ps0, psm
        subroutine update_ps
          use modglobal,          only: rdt, rk3step, lmoist
          use modsurfdata,        only : ps, qts
          use modsurface,         only : qtsurf


          if(rk3step == 1) then
             ps = ps + ps_tend * rdt
             
             ! modtimedep calls qtsurf to update surface values after changing ps, so we'll do the same
             if (lmoist .and. ps_tend /= 0) then
                call qtsurf
             else
                qts = 0.
             endif
          endif
        end subroutine update_ps
        

        ! Performs repeatedly time stepping until the tstop (measured in seconds after the cold start) is reached. 
        ! Note that the resulting model time may be a bit later than tstop.
        subroutine run_to(tstop)

            use modglobal,          only: timee,longint,rk3step

            implicit none

            integer(kind=longint),intent(in) :: tstop

            do while (timee < tstop .or. rk3step < 3)
                call step
            end do

        end subroutine run_to

        ! Returns the model time, measured in s after the cold start.
        function get_model_time() result(ret)

            use modglobal,          only: timee,longint

            implicit none

            integer(kind=longint) :: ret

            ret=timee

        end function get_model_time

        ! Cleans up; closes file handles ans deallocates arrays. This is a copy of the code 
        ! in the main program after the time loop.
        subroutine finalize

            !!----------------------------------------------------------------
            !!     0.0    USE STATEMENTS FOR CORE MODULES
            !!----------------------------------------------------------------
            use modstartup,         only : exitmodules

            !----------------------------------------------------------------
            !     0.1     USE STATEMENTS FOR ADDONS STATISTICAL ROUTINES
            !----------------------------------------------------------------
            use modcape,            only : exitcape
            use modgenstat,         only : exitgenstat
            use modradstat,         only : exitradstat
            use modlsmstat,         only : exitlsmstat
            use modsampling,        only : exitsampling
            use modquadrant,        only : exitquadrant
            use modcrosssection,    only : exitcrosssection  
            use modAGScross,        only : exitAGScross
            use modlsmcrosssection, only : exitlsmcrosssection
            use modcloudfield,      only : cloudfield
            use modfielddump,       only : exitfielddump
            use modsamptend,        only : exitsamptend

            use modbulkmicrostat,   only : exitbulkmicrostat
            use modbudget,          only : exitbudget
            use modheterostats,     only : exitheterostats

            ! modules below are disabled by default to improve compilation time
            !use modstress,         only : exitstressbudget

            !use modtilt,           only : exittilt
            !use modparticles,      only : exitparticles
            use modnudge,           only : exitnudge
            use modcanopy,          only : exitcanopy

            implicit none

            !--------------------------------------------------------
            !    4    FINALIZE ADD ONS AND THE MAIN PROGRAM
            !-------------------------------------------------------
            call exitgenstat
            call exitradstat
            call exitlsmstat
            !call exitparticles
            call exitnudge
            call exitsampling
            call exitquadrant
            call exitsamptend
            call exitbulkmicrostat
            call exitbudget
            !call exitstressbudget
            call exitcrosssection
            call exitAGScross
            call exitlsmcrosssection
            call exitcape
            call exitfielddump
            call exitheterostats
            call exitcanopy
            call exitmodules
            call exitdaleslib

        end subroutine finalize

        function grid_shape() result(s)

            use modglobal, only: itot,jtot,kmax

            implicit none

            integer:: s(3)

            s=(/itot,jtot,kmax/)

        end function

        function allocate_z_axis(a) result(ierr)

            use modglobal, only: kmax

            implicit none

            real, allocatable :: a(:)
            integer           :: ierr

            ierr=0
            
            if(allocated(a)) then
                deallocate(a,stat=ierr)
            endif
            
            if(ierr/=0) then
                return
            endif

            allocate(a(kmax),stat=ierr)

        end function

        function allocate_2d(a) result(ierr)

            use modglobal, only: itot,jtot

            implicit none

            real, allocatable :: a(:,:)
            integer           :: ierr

            ierr=0
            
            if(allocated(a)) then
                deallocate(a,stat=ierr)
            endif
            
            if(ierr/=0) then
                return
            endif

            allocate(a(itot,jtot),stat=ierr)

        end function
        
        function allocate_3d(a) result(ierr)

            use modglobal, only: itot,jtot,kmax

            implicit none

            real, allocatable :: a(:,:,:)
            integer           :: ierr

            ierr=0
            
            if(allocated(a)) then
                deallocate(a,stat=ierr)
            endif
            
            if(ierr/=0) then
                return
            endif

            allocate(a(itot,jtot,kmax),stat=ierr)

        end function


    ! Retrieves the z-layer average of the given array
    ! Al : input field, 3d.  
    ! Ag : output vector, 1d, in Z direction. Assumed to be as high as Al
    ! NOTE averages the FULL input array - if less is wanted, pass a slice !
    function gatherlayeravg(Al,Ag) result(ret)
      use modmpi, only: comm3d, myid, nprocs, D_MPI_REDUCE, MPI_SUM
      real, intent(in)      :: Al(:,:,:)
      real, intent(out)     :: Ag(:)
      integer               :: k, nk, ret

      ret = 0
      nk = size(Al, 3)
      Ag = (/ (sum(Al(:,:,k)), k=1,nk) /)     ! sum layers of Al
      
      !in-place reduction
      if (myid == 0) then
         CALL D_mpi_reduce( Ag, nk, MPI_SUM, 0, comm3d, ret)
      else
         CALL D_mpi_reduce(          Ag, Ag, nk, MPI_SUM, 0, comm3d, ret)
      endif
      
      if (myid == 0) then
         Ag = Ag / (size(Al,1) * size(Al,2) * nprocs)
      endif
    end function gatherlayeravg


    ! Counts the profile of saturated grid cell fraction and scatters the result to all processes
    function gathersatfrac(Ag) result(ret)
      use modmpi, only: comm3d, nprocs, D_MPI_ALLREDUCE, MPI_SUM
      use modglobal,   only : i1, j1
      use modfields, only: ql0
      real,    intent(out)    :: Ag(:)
      integer                 :: k, nk, ret

      ret = 0
      nk = size(ql0, 3)
      Ag = (/ (sum(merge(1., 0., ql0(2:i1,2:j1,k) > 0.)), k=1,nk) /)     ! sum layers of ql

      CALL D_mpi_allreduce( Ag, nk, MPI_SUM, comm3d, ret)

      Ag = Ag / (size(ql0,1) * size(ql0,2) * nprocs)
    
    end function gathersatfrac

    ! Retrieves the z-layer average of the ice fraction of the condensed water ql
    ! ilratio is calculated from the temperature, as in simpleice and icethermo routines.
    ! note: assumes the qi array is large enough (kmax elements)
    function gatherlayericeavg(qi) result(ret)
      use modmpi, only: comm3d, myid, nprocs, D_MPI_REDUCE, MPI_SUM
      use modglobal, only: imax, jmax, kmax, i1, j1, tup, tdn
      use modfields, only: ql0, tmp0
      
      real, intent(out)     :: qi(:)
      integer               :: i,j,k, ret
      real                  :: ilratio

      ret = 0
      do k=1,kmax
         qi(k) = 0
         do j = 2,j1
            do i = 2,i1
               ilratio = max(0.,min(1., (tmp0(i,j,k)-tdn) / (tup-tdn))) ! ice liquid ratio . 0 forice, 1 for liquid
               qi(k) = qi(k) + (1.0 - ilratio) * ql0(i,j,k)              ! amount of ice 
            enddo
         enddo
      enddo
      
      !in-place reduction
      if (myid == 0) then
         CALL D_mpi_reduce( qi, kmax, MPI_SUM, 0, comm3d, ret)
      else
         CALL D_mpi_reduce(          qi, qi, kmax, MPI_SUM, 0, comm3d, ret)
      endif
      
      if (myid == 0) then
         qi = qi / (imax * jmax * nprocs)
      endif
    end function gatherlayericeavg

    ! Retrieves the ice fraction of the condensed water ql
    ! ilratio is calculated from the temperature, as in simpleice and icethermo routines.
    ! note: assumes the qi array is large enough (kmax elements)
    function geticecontent(qi) result(ret)
      use modglobal, only: kmax, i1, j1, tup, tdn
      use modfields, only: ql0, tmp0

      real, intent(out)     :: qi(2:i1,2:j1,kmax)
      integer               :: i,j,k,ret
      real                  :: ilratio

      ret = 0
      do k=1,kmax
         do j = 2,j1
            do i = 2,i1
               ilratio = max(0., min(1., (tmp0(i,j,k) - tdn) / (tup - tdn))) ! ice liquid ratio . 0 forice, 1 for liquid
               qi(i,j,k) = (1.0 - ilratio) * ql0(i,j,k)              ! amount of ice
            enddo
         enddo
      enddo
    end function geticecontent


    ! map global indices to local
    ! gi = 1...itot  <-- one-based global indices
    ! gj = 1...jtot
    ! gk = 1...ktot
    ! i = 2...i1  <- two-based local indices, good for accessing the fields
    ! j = 2...j1
    function localindex(gi, gj, gk, i, j, k) result(ret)
      use modglobal, only: imax, jmax, kmax, i1, j1
      use modmpi, only: myidx, myidy
      integer, intent(in) :: gi, gj, gk
      integer, intent(out) :: i, j, k
      integer :: ret

      i = gi - myidx * imax + 1
      j = gj - myidy * jmax + 1
      k = gk
      
      if (i >= 2 .and. i <= i1 .and. j >= 2 .and. j <= j1 .and. k >= 1 .and. k <= kmax) then
         ret = 1  ! this point belongs to us
      else
         ret = 0
      endif
      
    end function localindex
    
    
    ! getter function for 3D field data
    ! g_i, g_j, g_k are index arrays
    ! data from the specified field is extracted and returned in a
    ! a(r) = field(g_(r), g_j(r), g_k(r))
    !
    ! Assumptions: the indexing arrays are one-based
    ! the array should be sliced to contain ONLY the physical cells before calling,
    ! so the data arrays are 1-based
    function gathervol(g_i,g_j,g_k,a,n,field) result(ret)
      use modglobal, only: imax, jmax
      use modmpi, only: myidx, myidy, comm3d, myid, D_MPI_REDUCE, MPI_SUM
      
      integer, intent(in)                 :: n
      integer, dimension(n), intent(in)   :: g_i,g_j,g_k
      real,    dimension(n), intent(out)  :: a
      real,    intent(in)                 :: field (:,:,:)
      integer                             :: r, i, j, k, is, js, ks, ret

      ret = 0
      is = size(field, 1)
      js = size(field, 2)
      ks = size(field, 3)

      ! store the data that belongs to "my" chunk or 0
      ! this is reduced with SUM, yielding the full data
      do r = 1,n
         i = g_i(r) - myidx * imax
         j = g_j(r) - myidy * jmax
         k = g_k(r)
         if (i >= 1 .and. i <= is .and. j >= 1 .and. j <= js .and. k >= 1 .and. k <= ks) then
            a(r) = field(i,j,k)
            !write(*,*)  '(', g_i(r), g_j(r), g_k(r), ') ', i, j, k, '(', myidx, myidy, ')'
         else
            a(r) = 0
         endif
      enddo
      
      !in-place reduction
      if (myid == 0) then
         CALL D_mpi_reduce( a, n, MPI_SUM, 0, comm3d, ret)
      else
         CALL D_mpi_reduce(           a, a, n, MPI_SUM, 0, comm3d, ret)
      endif
      
    end function gathervol


    ! getter function for 2D field data
    ! g_i, g_j, g_k are index arrays
    ! data from the specified field is extracted and returned in a
    ! a(r) = field(g_(r), g_j(r))
    !
    ! Assumptions: the indexing arrays are one-based
    ! the array should be sliced to contain ONLY the physical cells before calling,
    ! so the data arrays are 1-based
    function gatherlayer(g_i,g_j,a,n,field) result(ret)
      use modglobal, only: imax, jmax
      use modmpi, only: myidx, myidy, comm3d, myid, D_MPI_REDUCE, MPI_SUM

      integer, intent(in)                 :: n
      integer, dimension(n), intent(in)   :: g_i,g_j
      real,    dimension(n), intent(out)  :: a
      real,    intent(in)                 :: field (:,:)
      integer                             :: r, i, j, is, js, ret

      ret = 0
      is = size(field, 1)
      js = size(field, 2)

      ! store the data that belongs to "my" chunk or 0
      ! this is reduced with SUM, yielding the full data
      do r = 1,n
         i = g_i(r) - myidx * imax
         j = g_j(r) - myidy * jmax
         if (i >= 1 .and. i <= is .and. j >= 1 .and. j <= js) then
            a(r) = field(i,j)
         else
            a(r) = 0
         endif
      enddo

      !in-place reduction
      if (myid == 0) then
         CALL D_mpi_reduce( a, n, MPI_SUM, 0, comm3d, ret)
      else
         CALL D_mpi_reduce(           a, a, n, MPI_SUM, 0, comm3d, ret)
      endif
    end function gatherlayer

    ! gather the Liquid Water Path
    ! vertical sum of QL * rho * dz = integral dz * rho * q_l
    ! unit: kg / m^2
    !
    ! g_i, g_j are index arrays
    ! data from the specified field is extracted and returned in a
    ! a(r) = field(g_(r), g_j(r))
    
    ! Assumptions: the indexing arrays are one-based
    !              the data is stored in field(1:imax,1:jmax,1:kmax)
    ! the array should be sliced to contain ONLY the physical cells before calling,
    ! so the data arrays are 1-based
    function gatherlwp(g_i,g_j,a,n,field) result(ret)
      use modglobal, only: imax, jmax, kmax, dzf
      use modfields, only: rhobf
      use modmpi, only: myidx, myidy, comm3d, myid, D_MPI_REDUCE, MPI_SUM
      
      integer, intent(in)                 :: n
      integer, dimension(n), intent(in)   :: g_i,g_j
      real,    dimension(n), intent(out)  :: a
      real,    intent(in)                 :: field (:,:,:)
      integer                             :: r, i, j, is, js, ks, ret

      ret = 0
      is = size(field, 1)
      js = size(field, 2)
      ks = size(field, 3)

      ! store the data that belongs to "my" chunk or 0
      ! this is reduced with SUM, yielding the full data
      do r = 1,n
         i = g_i(r) - myidx * imax
         j = g_j(r) - myidy * jmax
         if (i >= 1 .and. i <= is .and. j >= 1 .and. j <= js ) then
            a(r) = sum(field(i,j,1:kmax) * dzf(1:kmax) * rhobf(1:kmax)) ! vertical sum of the field, weighted by cell height and density
         else
            a(r) = 0
         endif
      enddo
      
      !in-place reduction
      if (myid == 0) then
         CALL D_mpi_reduce( a, n, MPI_SUM, 0, comm3d, ret)
      else
         CALL D_mpi_reduce(           a, a, n, MPI_SUM, 0, comm3d, ret)
      endif
      
    end function gatherlwp

    ! Compute cloud fraction for a slabs of ql
    ! for every slab, count how many columns contain non-zero ql
    ! sum this count over MPI threads, then divide by the number of columns
    
    ! ql : input field, 3d.  
    ! I  : index array defining how layers are grouped.
    ! A  : output vector, 1d, in Z direction. Assumed to be as high as I
    ! A(i) = result for layers I(i-1) ... I(i)   I(0) is assumed to be 1
    

    ! important: the order of summing and counting matters
    ! we assume that any non-zero ql makes the column opaque
    ! then this cloud fraction is the fraction of opaque columns.
    
    ! NOTE averages the FULL input array - if less is wanted, pass a slice !
    function gathercloudfrac(ql,I,A) result(ret)
      use modmpi, only: comm3d, myid, nprocs, D_MPI_REDUCE, MPI_SUM
      real,    intent(in)     :: ql(:,:,:)
      integer, intent(in)     :: I(:)
      real,    intent(out)    :: A(:)
      integer                 :: ii, k1, k2, ret

      ret = 0
      k1 = 1                    ! start of k-range
      do ii = 1, size(I)
         k2 = I(ii)             ! end of k-range
         A(ii) = count (sum (ql(:,:,k1:k2-1), dim=3) > 0)  ! count how many columns in the current slab contains non-zero ql
         k1 = I(ii)
      enddo
      
      write(*,*) 'gatherCloudFrac partial (', myid, '): ', A
      
      !in-place reduction
      if (myid == 0) then
         CALL D_mpi_reduce( A, size(I), MPI_SUM, 0, comm3d, ret)
      else
         CALL D_mpi_reduce(          A,  A, size(I), MPI_SUM, 0, comm3d, ret)
      endif
      
      if (myid == 0) then
         A = A / (size(ql,1) * size(ql,2) * nprocs)
         write(*,*) 'gatherCloudFrac final:', A
      endif
    end function gathercloudfrac

    subroutine set_start_time(year,month,day,hour,minute,second)
        use modstat_nc, only : get_date, leap_year
        use modglobal, only  : xyear, xday, xtime
        implicit none
        integer, intent(in) :: year,month,day,hour,minute,second
        integer             :: mdays(12), mm, dd, i
        
        xyear = year
        if (month < 1 .or. month > 12) then
            write(*,*) 'warning: invalid month, applying modulo to get', mod(month,13)
        endif
        mm = mod(month,13)
        mdays = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
        if (leap_year(xyear)) then
            mdays(2) = 29
        endif
        if (day < 1 .or. day > mdays(mm)) then
            write(*,*) 'warning: invalid day, applying modulo to get', mod(day,mdays(mm) + 1)
        endif
        dd = mod(day,mdays(mm) + 1)
        xday = 0
        i = 1
        do i=1,mm - 1
            xday = xday + mdays(i)
        enddo
        xday = xday + dd
        xtime = float(hour) + minute/60. + second/3600. 
    end subroutine set_start_time
    
end module daleslib
