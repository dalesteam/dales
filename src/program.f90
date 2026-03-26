!> \file program.f90
!! Main program

!>
!! \mainpage
!! Dutch Atmospheric Large Eddy Simulation
!! \section DALES Dutch Atmospheric Large Eddy Simulation
!!
!! @version 5.0.0-beta
!!
!! @author
!! Steef Boing
!! (TU Delft)
!! \author
!! Huug Ouwersloot
!! (Wageningen University)
!! \author
!! Johan van der Dussen
!! (TU Delft)
!! \author
!! Steef B\"oing
!! (TU Delft)
!>
!! \section Log Change log
!! \par New Features
!! \par Main Changes
!! \todo

!! Notes
!! This subversion
!! Huug:
!! - Included heterosurf routine
!! - Statistics for heterosurf routine
!! Steef:
!! - Important note; adapted by Huug: ekm and ekh is again set to just Kh for right calculation of subgrid fluxes
!!   mosts statistic have been adjusted accordingly, however, budgets still need full update
!! - Anelastic baseprofile maker
!! - Anelastic advection
!! - Anelastic poisson solver
!! - Anelastic diffusion
!! - Resolved buoyancy (based on theta_l,q_l -> theta_v), using mean theta_v in divisor
!!   Subtracting mean state theta_v before Poisson solver
!! - Rainwater loading included in buoyancy (modforces)
!! - Simple ice microphysics scheme (Grabowski 98, with switches for autoconversion and graupel)
!! - Updated microstat for bulk and ice scheme
!! - Diagnostic temperature and saturation fields included, used to speed up micro (adjusted restart files accordingly)
!! - Speeded up gamma functions in bulkmicro and ice-micro using tabulation
!! - Reviewed saturation pressure with table lookup formula (Murphy and Koop, unified water/ice)
!! - Analytical functions for surface forcing (currently hard-coded)
!! - Larger fielddump range for temperatures
!! - Fixed statistics for heights above 10000 m
!! - Combined sampling/tendency routine (experimental)
!! - CAPE/CIN etc routine (experimental)
!! - CFL criterion based on pythagorean CFL
!! - Sampling written to separate netcdf files
!! - Modsampling update
!! - Radiation and bulkmicro tendencies exner function correction
!! - Consistent notation of theta_v in output
!! - Radiation negative qt crash
!! - Integrate WENO advection (Johan)
!! - Removed tqaver
!! - Subsidence with local values
!! - top boundary conditions (thl,qt-gradients) time-dependent
!! \par todo (this release)
!! - Scalasca CMake and Marmot options (Johan)
!! - Consistent modbudget and modgenstat with anelastic dynamics (Steef)
!! \par todo (future)
!! - General code cleanup
!! - Unified and simpler diagnostics
!! - Fielddump timing (Johan)
!! - Input header detection (Steef)
!! - Cleanup namoptions, remove dtav and timeav from some of the namoptions
!! - Check warm startup for interactive radiation cases
!! - 2D Parallelization
!! - Use more complicated theta_l formulation, include latent heat of freezing
!! - Adjust buoyancy and subgrid accordingly
!! - Integrate precipitation loading in theta_v
!! - Add 2-moment scheme? (Thijs working on complicated scheme, use Grabowski/Morrison?)
!!
!! \section License License
!!  This file is part of DALES.
!!
!!  DALES is free software; you can redistribute it and/or modify it under the
!! terms of the GNU General Public License as published by the Free Software
!! Foundation; either version 3 of the License, or (at your option) any later
!! version.
!!
!!  DALES is distributed in the hope that it will be useful, but WITHOUT ANY
!! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!! PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License along with
!! this program.  If not, see <http://www.gnu.org/licenses/>.
!!
!!  Copyright 1993-2009 Delft University of Technology, Wageningen University,
!! Utrecht University, KNMI
!!
program DALES

!!----------------------------------------------------------------
!!     0.0    USE STATEMENTS FOR CORE MODULES
!!----------------------------------------------------------------
  use modglobal,         only : rk3step,timeleft,lopenbc
  use modmpi,            only : initmpicomm
  use modstartup,        only : startup, writerestartfiles,testwctime,exitmodules
  use modtimedep,        only : timedep
  use modboundary,       only : boundary, grwdamp! JvdD ,tqaver
  use modthermodynamics, only : thermodynamics
  use modmicrophysics,   only : microphysics
  use modsurface,        only : surface
  use modlsm,            only : lsm
  use moddrydeposition,  only : drydep
  use modsubgrid,        only : subgrid
  use modforces,         only : forces, coriolis, lstend
  use modradiation,      only : radiation
  use modpois,           only : poisson
  use tstep,             only : tstep_update,  tstep_integrate, reset_tendencies
  use modlogging,        only : initlogging, exitlogging
  !use modedgecold,       only : coldedge

!----------------------------------------------------------------
!     0.1     USE STATEMENTS FOR ADDONS STATISTICAL ROUTINES
!----------------------------------------------------------------
  use modscalarpulse,  only : initscalarpulse, scalarpulse
  use modcape,         only : initcape,docape
  use modchecksim,     only : initchecksim, checksim
  use modstat_nc,      only : initstat_nc
  !use modspectra2,     only : dospecs,initspectra2,tanhfilter
  use modtimestat,     only : inittimestat, timestat, exittimestat
  use modgenstat,      only : initgenstat, genstat, exitgenstat
  use modradstat,      only : initradstat ,radstat, exitradstat
  use modlsmstat,      only : initlsmstat ,lsmstat, exitlsmstat
  !use moddepstat,      only : initdepstat ,depstat, exitdepstat
  use modsampling,     only : initsampling, sampling,exitsampling
  use modquadrant,     only : initquadrant, quadrant,exitquadrant
  use modcrosssection, only : initcrosssection, crosssection
  use modAGScross,     only : initAGScross, AGScross,exitAGScross
  use modlsmcrosssection, only : initlsmcrosssection, lsmcrosssection,exitlsmcrosssection
  use moddepcrosssection, only : initdepcrosssection, depcrosssection,exitdepcrosssection
  use modcloudfield,   only : initcloudfield, cloudfield
  use modfielddump,    only : initfielddump, fielddump
  use modradfield,     only : initradfield, radfield
  use modsamptend,     only : initsamptend, samptend,exitsamptend, tend_start,tend_subg,tend_force,&
                              tend_rad,tend_ls,tend_micro, tend_topbound,tend_pois,tend_addon, tend_coriolis,&
                              leibniztend, writesamptend

  use modbudget,       only : initbudget, budgetstat, exitbudget
  use modheterostats,  only : initheterostats, heterostats, exitheterostats
  use modvarbudget,    only : initvarbudget, varbudget, exitvarbudget
  use modmsebudg,      only : initmsebudg, msebudg1, msebudg2, exitmsebudg
  ! modules below are disabled by default to improve compilation time
  !use modstress,       only : initstressbudget, stressbudgetstat, exitstressbudget

  !use modtilt,         only : inittilt, tiltedgravity, tiltedboundary, exittilt
  !use modparticles,    only : initparticles, particles, exitparticles
  use modnudge,        only : initnudge, nudge, exitnudge
  use modnudgeboundary, only : initnudgeboundary, nudgeboundary, exitnudgeboundary
  use modtestbed,      only : testbednudge, exittestbed
  !use modprojection,   only : initprojection, projection
  use modchem,         only : initchem,twostep
  use modcanopy,       only : initcanopy, canopy, exitcanopy
  use modadvection,    only : advection
  use moddatetime,     only : datetime
  use modemission,     only : emission
  use modopenboundary, only : openboundary_ghost,openboundary_tend,openboundary_phasevelocity,openboundary_turb
  use modstat_profiles, only: init_profiles, sample_profiles, write_profiles, exit_profiles
  use modibm,          only : applyibm, zerowallvelocity
  use modibmdata,      only : lpoislast
  use modlatsponge,    only : lateral_sponge
  use modspraying,     only : spraying
  use modprecursor,    only : init_precursor, precursor_nudge_boundary, &
                              swap_fields, exit_precursor, &
                              lprecursor, Nsim, statid, turid, refid
  use modcloudstat,    only: init_cloudstat, do_cloudstat
  use modstat_nc_files, only: stats_limit_timestep, init_output_files, write_output_files
!----------------------------------------------------------------
!     0.2     USE STATEMENTS FOR TIMER MODULE
!----------------------------------------------------------------

  use modtimer,       only : timer_tic, timer_toc, timer_print, timer_write, timer_cleanup

!----------------------------------------------------------------
!     0.3     USE STATEMENTS FOR GPU UTILITIES
!----------------------------------------------------------------

#if defined(_OPENACC)
  use modgpu, only: update_gpu, host_is_updated
#endif

  implicit none

  integer :: istep
  integer :: simid !< Simulation ID, used for precursor simulations

  ! Select CPU for execution of startup routines
!----------------------------------------------------------------
!     1      READ NAMELISTS,INITIALISE GRID, CONSTANTS AND FIELDS
!----------------------------------------------------------------
  ! call initmpi initmpi depends on options in the namelist, call moved to startup
  call initmpicomm
  call initlogging
  call startup

!---------------------------------------------------------
!      2     INITIALIZE STATISTICAL ROUTINES AND ADD-ONS
!---------------------------------------------------------
  call initchecksim
  call initstat_nc   ! Should be called before stat-routines that might do netCDF
  call inittimestat  ! Timestat must preceed all other timeseries that could write in the same netCDF file (unless stated otherwise
  call initgenstat   ! Genstat must preceed all other statistics that could write in the same netCDF file (unless stated otherwise
  !call inittilt
  call initquadrant
  call initcrosssection
  call initAGScross
  call initlsmcrosssection
  call initdepcrosssection
  !call initprojection
  call initcloudfield
  call initradstat
  call initradfield
  call initlsmstat
  !call initdepstat
  !call initparticles
  call initnudge
  call initnudgeboundary
  call initbudget
  call initvarbudget
  call initmsebudg
  !call initstressbudget
  ! call initchem
  call initsampling
  call initfielddump
  call initsamptend
  call initheterostats
  call initcanopy
  !call initspectra2
  call initscalarpulse
  call initcape
  call init_cloudstat

  call init_profiles
  call init_precursor

#if defined(_OPENACC)
  call update_gpu
#endif

  ! Initialize IO
  call init_output_files

!------------------------------------------------------
!   3.0   MAIN TIME LOOP
!------------------------------------------------------
  call testwctime
  istep = 1
  do while (timeleft > 0)
    do simid = 1, Nsim

      if (simid == refid) call tstep_update

      do rk3step = 1, 3
        call timer_tic('program/timestep', istep)
    
    
        ! Calculate new timestep, and reset tendencies to 0.
        call timedep
        call scalarpulse
        call samptend(tend_start,firstterm=.true.)
        call datetime
    
        ! Check if we have to sample profiles this time step
        call sample_profiles
    
        call datetime
    
    !-----------------------------------------------------
    !   3.1   Openboundaries
    !-----------------------------------------------------
        if(lopenbc) then
          call openboundary_turb
          call openboundary_ghost
          call openboundary_tend
        endif
    
    !-----------------------------------------------------
    !   3.2   RADIATION
    !-----------------------------------------------------
        call radiation !radiation scheme
        call samptend(tend_rad)
    
    !-----------------------------------------------------
    !   3.3   THE SURFACE LAYER / LAND-SURFACE
    !-----------------------------------------------------
        call lsm
        call drydep
        call surface
    
    !-----------------------------------------------------
    !   3.4   ADVECTION AND DIFFUSION
    !-----------------------------------------------------
        call advection
        call subgrid
        call canopy
        call samptend(tend_subg)
    
    !-----------------------------------------------------
    !   3.5   REMAINING TERMS
    !-----------------------------------------------------
        call coriolis !remaining terms of ns equation
        call samptend(tend_coriolis)
        call forces !remaining terms of ns equation
        call samptend(tend_force)
    
        call lstend !large scale forcings
        call samptend(tend_ls)
        call microphysics
        call samptend(tend_micro)
        call emission
    
    !------------------------------------------------------
    !   3.6   EXECUTE ADD ONS
    !------------------------------------------------------
        call nudge
        call nudgeboundary
        call testbednudge
        if (simid == turid) call spraying
    !    call dospecs
    !    call tiltedgravity
    
        call samptend(tend_addon)
    
        if (lprecursor .and. simid == turid) call precursor_nudge_boundary
    
    !-----------------------------------------------------------------------
    !   3.7  PRESSURE FLUCTUATIONS, TIME INTEGRATION AND BOUNDARY CONDITIONS
    !-----------------------------------------------------------------------
        call grwdamp !damping at top of the model
    !JvdD    call tqaver !set thl, qt and sv(n) equal to slab average at level kmax
        call samptend(tend_topbound)
    
        ! either apply ibm before or after poisson solver
        if (lpoislast .eqv.  .true.) call applyibm
        if (lpoislast .eqv. .false.) call zerowallvelocity ! put wall velocities to zero before Poisson
        call poisson
    
        if (lpoislast .eqv. .false.) call applyibm ! then only apply IBM after Poisson
    
        call samptend(tend_pois,lastterm=.true.)
        if(lopenbc) call openboundary_phasevelocity()
    
        call lateral_sponge
    
        call tstep_integrate                        ! Apply tendencies to all variables
    
        call msebudg1
        ! NOTE: the tendencies are not zeroed yet, but kept for analysis and statistcis
        !       Do not change them below this point.
        if(lopenbc) then
          call openboundary_ghost
        else
          call boundary
        endif
    
    
        !call tiltedboundary
    !-----------------------------------------------------
    !   3.8   LIQUID WATER CONTENT AND DIAGNOSTIC FIELDS
    !-----------------------------------------------------
        call thermodynamics
        call leibniztend
        call writesamptend
    !-----------------------------------------------------
    !   3.9  WRITE RESTARTFILES AND DO STATISTICS
    !------------------------------------------------------
        if (simid == statid) then
          call stats_limit_timestep
          call twostep
          !call coldedge
          call checksim
          call timestat  !Timestat must preceed all other timeseries that could write in the same netCDF file (unless stated otherwise
          call genstat  !Genstat must preceed all other statistics that could write in the same netCDF file (unless stated otherwise
          call write_profiles
          call radstat
          call lsmstat
          !call depstat
          call sampling
          call quadrant
          call crosssection
          call AGScross
          call lsmcrosssection
          call depcrosssection
          !call tanhfilter
          call docape
          !call projection
          call cloudfield
          call fielddump
          call radfield
          !call particles

          call do_cloudstat
    
          call budgetstat
          call varbudget
          call msebudg2
          !call stressbudgetstat
          call heterostats
    
          call testwctime
          call writerestartfiles

          call write_output_files
        end if

        call reset_tendencies

#if defined(_OPENACC)
        host_is_updated = .false.
#endif
        call timer_toc('program/timestep')
      end do ! rk3step

      if (lprecursor) call swap_fields
    end do ! simid
    istep = istep + 1
  end do ! time loop

!-------------------------------------------------------
!             END OF TIME LOOP
!-------------------------------------------------------

  call timer_print
  call timer_write
  call timer_cleanup

!--------------------------------------------------------
!    4    FINALIZE ADD ONS AND THE MAIN PROGRAM
!-------------------------------------------------------
  call exitgenstat
  call exitradstat
  call exitlsmstat
  !call exitdepstat
  !call exitparticles
  call exitnudge
  call exitnudgeboundary
  call exittestbed
  call exitsampling
  call exitquadrant
  call exitsamptend
  call exitbudget
  call exitvarbudget
  call exitmsebudg
  !call exitstressbudget
  call exitAGScross
  call exitlsmcrosssection
  call exitdepcrosssection
  call exitheterostats
  call exitcanopy
  call exittimestat
  call exitnudgeboundary  !cstep
  call exitmodules
  call exit_profiles
  call exitlogging


end program DALES
