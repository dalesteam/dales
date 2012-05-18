!> \file program.f90
!! Main program

!>
!! \mainpage
!! Dutch Atmospheric Large Eddy Simulation
!! \section DALES Dutch Atmospheric Large Eddy Simulation
!!
!! @version 4.0.1alpha
!!
!! @author
!! Steef Boing
!! (TU Delft)
!! \author
!! Chiel van Heerwaarden
!! (Wageningen University)
!! \author
!! Thijs Heus
!! (Max Planck Institute Hamburg)
!! \author
!! Steef B\"oing
!! (TU Delft)
!>
!! \section Log Change log
!! \par New Features
!! \par Main Changes
!! \todo

!! This subversion
!! Huug:
!! - Included heterosurf routine
!! - Statistics for heterosurf routine
!! Steef:
!! - Important note; ekm and ekh now denote rhobf*Kh for computational efficiency
!!    mosts statistic have been adjusted accordingly, however, budgets still need full update
!! - Anelastic baseprofile maker
!! - Anelastic advection
!! - Anelastic poisson solver changes
!! - Anelastic diffusion
!! - Resolved buoyancy (based on theta_l,q_l -> theta_v), using mean theta_v in divisor
!!    Also subtracting mean state theta_v before Poisson solver
!! - Rainwater loading included in buoyancy
!! - Simple ice microphysics scheme (Grabowski 98, with switches for autoconversion and graupel)
!! - Integrated microstat for bulk and ice scheme
!! - Diagnostic temperature field included
!! - Speeded up gamma functions in bulkmicro and ice-micro using tabulation
!! - Reviewed saturation pressure with table lookup formula (Murphy and Koop, unified water/ice)
!! - Analytical functions for surface forcing (may be removed later)
!! - Larger fielddump range for temperatures
!! - Fixed statistics for heights above 10000 m
!! - Combined sampling/tendency routine (experimental)
!! - CAPE/CIN etc routine (experimental)
!! - CFL criterion based on pythagorean CFL
!! - Sampling written to separate netcdf files
!! - Modsampling update (Stephan)
!! \par todo (this version)
!! - Check/fix warm startup (talk to Thijs)
!! - Check tqaver necessity (Johan coordinates?)
!! - Input header detection (Steef?)
!! - Redundant output by all processors (Johan)
!! - Add moddeprecated (Steef), cleaner namoptions
!! - Remove dtav and timeav from some of the namoptions (Steef)
!! - Look into radiative tendencies exner function correction (Steef, Johan)
!! - Radiation negative qt crash (Johan)
!! - Integrate WEENO advection (Johan)
!! - Unified and simpler diagnostics (Johan)
!! - Fielddump timing (Johan)
!! - Scalasca CMake and Marmot options (Johan)
!! - Consistent notation of theta_v in output (Steef)
!! - Consistent modbudget and modgenstat with anelastic dynamics (Steef)
!! - Separate microphysics from scalars (Steef)
!! - 2D Parallelization
!! \par todo (future)
!! Steef:
!! - Use more complicated theta_l formulation, include latent heat of freezing
!! - Adjust buoyancy and subgrid accordingly
!! - Integrate precipitation loading in theta_v
!! - Look at possible extra term in TKE equation (consistent filtering approach)
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
program DALES      !Version 4.0.0alpha

!!----------------------------------------------------------------
!!     0.0    USE STATEMENTS FOR CORE MODULES
!!----------------------------------------------------------------
  use modmpi,            only : myid, initmpi
  use modglobal,         only : rk3step,timee,btime,runtime,timeleft
  use modfields,         only : thl0
  use modstartup,        only : startup, writerestartfiles,exitmodules
  use modtimedep,        only : timedep
  use modboundary,       only : boundary, grwdamp,tqaver
  use modthermodynamics, only : thermodynamics
  use modmicrophysics,   only : microsources
  use modsurface,        only : surface
  use modsubgrid,        only : subgrid
  use modforces,         only : forces, coriolis, lstend
  use modradiation,      only : radiation
  use modpois,           only : poisson
  !use modedgecold,       only : coldedge

!----------------------------------------------------------------
!     0.1     USE STATEMENTS FOR ADDONS STATISTICAL ROUTINES
!----------------------------------------------------------------
  use modcape,         only : initcape,exitcape,docape
  use modchecksim,     only : initchecksim, checksim
  use modstat_nc,      only : initstat_nc
  !use modspectra2,     only : dospecs,initspectra2,tanhfilter
  use modtimestat,     only : inittimestat, timestat
  use modgenstat,      only : initgenstat, genstat, exitgenstat
  use modradstat,      only : initradstat ,radstat, exitradstat
  use modlsmstat,      only : initlsmstat ,lsmstat, exitlsmstat
  use modsampling,     only : initsampling, sampling,exitsampling
  use modcrosssection, only : initcrosssection, crosssection,exitcrosssection  
  use modlsmcrosssection, only : initlsmcrosssection, lsmcrosssection,exitlsmcrosssection
  use modcloudfield,   only : initcloudfield, cloudfield
  use modfielddump,    only : initfielddump, fielddump,exitfielddump
  use modsamptend,     only : initsamptend, samptend,exitsamptend, tend_start,tend_adv,tend_subg,tend_force,&
                              tend_rad,tend_ls,tend_micro, tend_topbound,tend_pois,tend_addon, tend_coriolis,leibniztend

  use modbulkmicrostat,only : initbulkmicrostat, bulkmicrostat,exitbulkmicrostat
  use modbudget,       only : initbudget, budgetstat, exitbudget
  use modheterostats,  only : initheterostats, heterostats, exitheterostats

  ! modules below are disabled by default to improve compilation time
  !use modstress,       only : initstressbudget, stressbudgetstat, exitstressbudget

  !use modtilt,         only : inittilt, tiltedgravity, tiltedboundary, exittilt
  !use modparticles,    only : initparticles, particles, exitparticles
  use modnudge,        only : initnudge, nudge, exitnudge
  !use modprojection,   only : initprojection, projection
  use modchem,         only : initchem,twostep


  implicit none

!----------------------------------------------------------------
!     1      READ NAMELISTS,INITIALISE GRID, CONSTANTS AND FIELDS
!----------------------------------------------------------------
  call initmpi

  call startup

!---------------------------------------------------------
!      2     INITIALIZE STATISTICAL ROUTINES AND ADD-ONS
!---------------------------------------------------------
  call initchecksim
  call initstat_nc   ! Should be called before stat-routines that might do netCDF
  call inittimestat  ! Timestat must preceed all other timeseries that could write in the same netCDF file (unless stated otherwise
  call initgenstat   ! Genstat must preceed all other statistics that could write in the same netCDF file (unless stated otherwise
  !call inittilt
  call initsampling
  call initcrosssection
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

  !call initspectra2
  call initcape


!------------------------------------------------------
!   3.0   MAIN TIME LOOP
!------------------------------------------------------
  write(*,*)'START myid ', myid
  do while (timeleft>0 .or. rk3step < 3)
    call tstep_update
    call timedep
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
    call subgrid
    call samptend(tend_subg)

!-----------------------------------------------------
!   3.4   REMAINING TERMS
!-----------------------------------------------------
    call coriolis !remaining terms of ns equation
    call samptend(tend_coriolis)
    call forces !remaining terms of ns equation
    call samptend(tend_force)

    call lstend !large scale forcings
    call samptend(tend_ls)
    call microsources !Drizzle etc.
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
    call tqaver !set thl, qt and sv(n) equal to slab average at level kmax
    call samptend(tend_topbound)
    call poisson
    call samptend(tend_pois,lastterm=.true.)

    call tstep_integrate
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
    call crosssection
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
  end do

!-------------------------------------------------------
!             END OF TIME LOOP
!-------------------------------------------------------


!--------------------------------------------------------
!    4    FINALIZE ADD ONS AND THE MAIN PROGRAM
!-------------------------------------------------------
  call exitgenstat
  call exitradstat
  call exitlsmstat
  !call exitparticles
  call exitnudge
  call exitsampling
  call exitsamptend
  call exitbulkmicrostat
  call exitbudget
  !call exitstressbudget
  call exitcrosssection
  call exitlsmcrosssection
  call exitcape
  call exitfielddump
  call exitheterostats
  call exitmodules

end program DALES
