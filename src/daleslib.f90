subroutine initialize()

    !!----------------------------------------------------------------
    !!     0.0    USE STATEMENTS FOR CORE MODULES
    !!----------------------------------------------------------------
    use modmpi,             only : initmpi
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
    use modcrosssection,    only : initcrosssection 
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

end subroutine initialize

subroutine step()

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
    use modcrosssection,    only : crosssection
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

    implicit none

    call tstep_update                           ! Calculate new timestep
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
    !JvdD    call tqaver !set thl, qt and sv(n) equal to slab average at level kmax
    call samptend(tend_topbound)
    call poisson
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

end subroutine step

subroutine finalize()

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
    use modcrosssection,    only : exitcrosssection  
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

end subroutine finalize
