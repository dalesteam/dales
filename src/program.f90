!> \file program.f90
!! Main program

!>
!! \mainpage
!! Dutch Atmospheric Large Eddy Simulation
!! \section DALES Dutch Atmospheric Large Eddy Simulation
!!
!! @version 4.4
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
  use modglobal,         only : rk3step,timeleft,is_starting
  use modmpi,            only : initmpicomm!, myid
  use modstartup,        only : startup, writerestartfiles,testwctime,exitmodules
  use modtimedep,        only : timedep
  use modboundary,       only : boundary, grwdamp! JvdD ,tqaver
  use modthermodynamics, only : thermodynamics
  use modmicrophysics,   only : microsources
  use modsurface,        only : surface
  use modsubgrid,        only : subgrid
  use modforces,         only : forces, coriolis, lstend
  use modradiation,      only : radiation
  use modpois,           only : poisson
  use tstep,             only : tstep_update,  tstep_integrate
  !use modedgecold,       only : coldedge

!----------------------------------------------------------------
!     0.1     USE STATEMENTS FOR ADDONS STATISTICAL ROUTINES
!----------------------------------------------------------------
  use modcape,         only : initcape,exitcape,docape
  use modchecksim,     only : initchecksim, checksim
  use modstat_nc,      only : initstat_nc
  !use modspectra2,     only : dospecs,initspectra2,tanhfilter
  use modtimestat,     only : inittimestat, timestat, exittimestat
  use modgenstat,      only : initgenstat, genstat, exitgenstat
  use modradstat,      only : initradstat ,radstat, exitradstat
  use modlsmstat,      only : initlsmstat ,lsmstat, exitlsmstat
  use modsampling,     only : initsampling, sampling,exitsampling
  use modquadrant,     only : initquadrant, quadrant,exitquadrant
  use modcrosssection, only : initcrosssection, crosssection,exitcrosssection
  use modAGScross,     only : initAGScross, AGScross,exitAGScross
  use modlsmcrosssection, only : initlsmcrosssection, lsmcrosssection,exitlsmcrosssection
  use modcloudfield,   only : initcloudfield, cloudfield
  use modfielddump,    only : initfielddump, fielddump,exitfielddump
  use modradfield,     only : initradfield, radfield, exitradfield
  use modsamptend,     only : initsamptend, samptend,exitsamptend, tend_start,tend_adv,tend_subg,tend_force,&
                              tend_rad,tend_ls,tend_micro, tend_topbound,tend_pois,tend_addon, tend_coriolis,leibniztend

  use modbulkmicrostat,only : initbulkmicrostat, bulkmicrostat,exitbulkmicrostat
  use modbudget,       only : initbudget, budgetstat, exitbudget
  use modheterostats,  only : initheterostats, heterostats, exitheterostats
  use modvarbudget,    only : initvarbudget, varbudget, exitvarbudget
  ! modules below are disabled by default to improve compilation time
  !use modstress,       only : initstressbudget, stressbudgetstat, exitstressbudget

  !use modtilt,         only : inittilt, tiltedgravity, tiltedboundary, exittilt
  !use modparticles,    only : initparticles, particles, exitparticles
  use modnudge,        only : initnudge, nudge, exitnudge
  use modtestbed,      only : testbednudge, exittestbed
  !use modprojection,   only : initprojection, projection
  use modchem,         only : initchem,twostep
  use modcanopy,       only : initcanopy, canopy, exitcanopy
  use modadvection,    only : advection


!----------------------------------------------------------------
!     0.2     USE STATEMENTS FOR NVTX TIMER MODULE
!----------------------------------------------------------------

  use modtimer,       only : timer_tic, timer_toc, timer_print

!----------------------------------------------------------------
!     0.3     USE STATEMENTS FOR GPU UTILITIES
!----------------------------------------------------------------

#if defined(_OPENACC)
  use modgpu, only: update_gpu, host_is_updated
#endif

  implicit none

  integer :: istep

  ! Select CPU for execution of startup routines
!----------------------------------------------------------------
!     1      READ NAMELISTS,INITIALISE GRID, CONSTANTS AND FIELDS
!----------------------------------------------------------------
  ! call initmpi initmpi depends on options in the namelist, call moved to startup
  call initmpicomm
  call startup

!---------------------------------------------------------
!      2     INITIALIZE STATISTICAL ROUTINES AND ADD-ONS
!---------------------------------------------------------
  call timer_tic('Initialization', 0)
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
  call initradfield
  call initlsmstat
  !call initparticles
  call initnudge
  call initbulkmicrostat
  call initbudget
  call initvarbudget
  !call initstressbudget
! call initchem
  call initheterostats
  call initcanopy

  !call initspectra2
  call initcape

#if defined(_OPENACC)
  call update_gpu
#endif

  ! Startup is done, set flag to false
  is_starting = .false.

  call timer_toc('Initialization')

!------------------------------------------------------
!   3.0   MAIN TIME LOOP
!------------------------------------------------------
  call testwctime

  istep = 1
  do while (timeleft>0 .or. rk3step < 3)
    call timer_tic('Time step', istep)
    ! Calculate new timestep, and reset tendencies to 0.
    call timer_tic('Time step update', 0)
    call tstep_update
    call timer_toc('Time step update')
    call timer_tic('Time dep', 0)
    call timedep
    call timer_toc('Time dep')
    call samptend(tend_start,firstterm=.true.)

!-----------------------------------------------------
!   3.1   RADIATION
!-----------------------------------------------------
    call timer_tic('Radiation', 0)
    call radiation !radiation scheme
    call timer_toc('Radiation')
    call samptend(tend_rad)

!-----------------------------------------------------
!   3.2   THE SURFACE LAYER
!-----------------------------------------------------
    call timer_tic('Surface', 0)
    call surface
    call timer_toc('Surface')

!-----------------------------------------------------
!   3.3   ADVECTION AND DIFFUSION
!-----------------------------------------------------
    call timer_tic('Advection', nvtx_color='g')
    call advection
    call timer_toc('Advection')
    call samptend(tend_adv)
    call timer_tic('Subgrid', nvtx_color='b')
    call subgrid
    call timer_toc('Subgrid')
    call timer_tic('Canopy', 0)
    call canopy
    call timer_toc('Canopy')
    call samptend(tend_subg)

!-----------------------------------------------------
!   3.4   REMAINING TERMS
!-----------------------------------------------------
    call timer_tic('Coriolis', 0)
    call coriolis !remaining terms of ns equation
    call timer_toc('Coriolis')
    call samptend(tend_coriolis)
    call timer_tic('Forces', 0)
    call forces !remaining terms of ns equation
    call timer_toc('Forces')
    call samptend(tend_force)

    call timer_tic('LStend', 0)
    call lstend !large scale forcings
    call timer_toc('LStend')
    call samptend(tend_ls)
    call timer_tic('Microsources', 0)
    call microsources !Drizzle etc.
    call timer_toc('Microsources')
    call samptend(tend_micro)

!------------------------------------------------------
!   3.4   EXECUTE ADD ONS
!------------------------------------------------------
    call nudge
    call testbednudge
!    call dospecs
!    call tiltedgravity

    call samptend(tend_addon)

!-----------------------------------------------------------------------
!   3.5  PRESSURE FLUCTUATIONS, TIME INTEGRATION AND BOUNDARY CONDITIONS
!-----------------------------------------------------------------------
    call timer_tic('GRWdamp', 0)
    call grwdamp !damping at top of the model
    call timer_toc('GRWdamp')
!JvdD    call tqaver !set thl, qt and sv(n) equal to slab average at level kmax
    call samptend(tend_topbound)
    call timer_tic('Poisson', nvtx_color='y')
    call poisson
    call timer_toc('Poisson')
    call samptend(tend_pois,lastterm=.true.)

    ! Apply tendencies to all variables
    call timer_tic('Tstep integrate', nvtx_color='c')
    call tstep_integrate
    call timer_toc('Tstep integrate')
    ! NOTE: the tendencies are not zeroed yet, but kept for analysis and statistcis
    !       Do not change them below this point.
    call timer_tic('Boundary', 0)
    call boundary
    call timer_toc('Boundary')
    !call tiltedboundary
!-----------------------------------------------------
!   3.6   LIQUID WATER CONTENT AND DIAGNOSTIC FIELDS
!-----------------------------------------------------
    call timer_tic('Thermodynamics', 0)
    call thermodynamics
    call timer_toc('Thermodynamics')
    call leibniztend
!-----------------------------------------------------
!   3.7  WRITE RESTARTFILES AND DO STATISTICS
!------------------------------------------------------
    call timer_tic('Statistics', 0)
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
    call radfield
    !call particles

    call bulkmicrostat
    call budgetstat
    call varbudget
    !call stressbudgetstat
    call heterostats
    call timer_toc('Statistics')
    
    call timer_tic('Restartfiles', 0)
    call testwctime
    call writerestartfiles
    call timer_toc('Restartfiles')
    call timer_toc('Time step') 
    istep = istep + 1
#if defined(_OPENACC)
    host_is_updated = .false.
#endif
  end do

!-------------------------------------------------------
!             END OF TIME LOOP
!-------------------------------------------------------

  call timer_print

!--------------------------------------------------------
!    4    FINALIZE ADD ONS AND THE MAIN PROGRAM
!-------------------------------------------------------
  call exitgenstat
  call exitradstat
  call exitlsmstat
  !call exitparticles
  call exitnudge
  call exittestbed
  call exitsampling
  call exitquadrant
  call exitsamptend
  call exitbulkmicrostat
  call exitbudget
  call exitvarbudget
  !call exitstressbudget
  call exitcrosssection
  call exitAGScross
  call exitlsmcrosssection
  call exitcape
  call exitfielddump
  call exitradfield
  call exitheterostats
  call exitcanopy
  call exittimestat
  call exitmodules


end program DALES
