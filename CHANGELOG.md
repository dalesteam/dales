Changes in DALES
================

Version 4.3 - 2020
------------------

This version introduces a library interface to DALES, defined in `daleslib.f90.` With this interface,
DALES can be used with the [OMUSE](https://omuse.readthedocs.io/en/latest/) framework, which provides
a Python interface to DALES. Another new feature is the support of iterative Poisson solvers, see
the [Wiki](https://github.com/dalesteam/dales/wiki/Iterative-Poisson-solver)


### Improvements

* Add a namelist option lsync in NAMNETCDFSTATS to synchronize netCDF output files after writing
* Iterative Poisson solver support with the HYPRE library (by J. Attema).
* Don't overwrite ascii output files on warm start. Commit 25fae048
* Warn when moduser routines are invoked from a case, when a custom moduser.f90 has not been compiled in. Commit ae148a3bd3
* Allow advection schemes 5th and 6th with non-uniform grid again, add experimental kappa scheme (77) with better non-uniform grid support. Commit 4ab9a7571b
* Add missing rho in divergence diagnostic, show reason for dt limit. Commit da3a32f951b5
* Add qsat field for statistics, and  statistics for accumulated surface rain flux in simpleice2.  Commit 6036c1a93b
* Add a library interface to DALES, to enable interfacing with OMUSE. Commits 6e0e39f2, 9f91f2d, 4b97aca, 640977e


### Optimization

* modpois: don't exchange the pwp field, ghost cells in it are not needed. Commit c7f72b1b
* If lcoriol is .false., skip coriolis calculations. Commit c343772
* advec_5th: move i,j loops inside if k. Advection calls ~3 times faster. Commit 0badd3a124
* Asynchronous MPI halo exchanges. E,W in parallel, then N,S. By J. Attema. Commits a54b9f8e, 6951323, 0c08a262


### Bugs fixed

* radiation and restarts : https://github.com/dalesteam/dales/issues/40
* modsimpleice2 and initialization order, https://github.com/dalesteam/dales/issues/49 . Commit be3f9711bb6
* timedependent forcings, finding the correct interval for interpolation with time dependent surface variables. https://github.com/dalesteam/dales/issues/48 . Commit d6f2225
* modbudget: remove extra factor rhobf in calculation of sbtke average. From S. de Roode.
* modbudget: add if(myid==0) around netcdf profile writing. Commit b7afc418c5
* Fix ibas_prf=3 initialization when zf(k) is exactly 11000 (1st value in the lapse rate table).
  [Issue #41](https://github.com/dalesteam/dales/issues/41), Commit ee6230bc00
* modradstat: handle sign conventions of different radiation schemes (S. de Roode). Commit f43012
* radiation: initialize {sw,lw}{Up,Down}_slice arrays before rrtmg calls, since the calls are not always performed (S. de Roode). Commit d4e979de7d
* netCDF staistics: fix for determining the index where writing is resumed at warm start.
  [Issue #42](https://github.com/dalesteam/dales/issues/42), Commit 8ff2cced3
* Fix closing of crosssection netcdf files. Commit e72701e08


Version 4.2 - 2019-06-05
------------------------

This is a summary of changes in DALES 4.2. More detials of the chanes follow below.

* Optimization of advection and subgrid schemes
* Improvement surface scheme (phi functions) for stable conditions
* Further testing RRTMG (introducing  aerosol)
* Extension from one-big leaf to two big-leaves (sunlit and shaded)
* Arbitrary distribution of sources and sinks scalar
* Heterogeneity land-surface to account for canopy dynamic effects
* Option for more agressive optimization flags (try export SYST=gnu-fast)



### Bugs fixed

* Typo in advec_5th and advec_52, in the advecv routines.
The 5th scheme has the typo only in the code for the boundaries, while 52 has it also for the bulk.
[Issue #37](https://github.com/dalesteam/dales/issues/37)
Commit c1ac043a038f7

* Obukhov length calculation: catch case Rib = 0 - set L=1e6. Can happen in the first step if the surface flux is 0. The Obukhov length is capped at 1e6 anyway, so we may as well cap it also when Rib == 0.
[Issue #24](https://github.com/dalesteam/dales/issues/24)
Commit a3c5d1b

* adaptive time step: limit dt by ekh, not only ekm - ekh may be larger.
[Issue #29](https://github.com/dalesteam/dales/issues/29)
Commit a67cf8f, 9f8304d83fe4

* Test cases arm_brown, arm_unstable, bomex : correcting namelist entry SUBGRID to NAMSUBGRID. These three test cases had an entry &SUBGRID in the name list. This entry is ignored, since the proper name is &NAMSUBGRID. All three cases also defined subgrid parameters in the name list which are not expected there. These were removed, as the values were set to the defaults.
Commit 12d910a

* Add kind=longint to floor and ceiling functions of time to prevent overflow in timeleft.
Commit 39f75a78e


### Improvements


* Add alternative simpleice - imicro=6. Contains fixes and optimizations but less statistics.  Note that the original simpleice is still preserved and can be used with imicro=5 as before.
Commit 52ca4d5.
    * snow ventilation parameters differ from references [Issue #26](https://github.com/dalesteam/dales/issues/26)
    * many loops were fused -> need fewer arrays, better cache efficiency
    * rsgratio was sometimes used uninitialized
    * avoid some numerical fixes: rsgratio(i,j,k)+1.e-6, qrr/(qrr+1.e-9)
    * performance gain: ~10% faster simulation in total, case-dependent.

* Add support for non-equidistant grid to _52 and _62 advection
schemes. For these two schemes, support can be added with copy-paste
from the 2nd order scheme, since also the 52 and 62 schemes are 2nd
order in the vertical direction. Commits 4534fb8c638b, bc9606e86b5d

* Stop with error message if the chosen advection is incompatible with
non-equidistant grid.  Most advection routines do not support a
non-equidistant grid. Issue: The list of schemes compatible with a non-equidistant grid was overly restrictive. Propose to enable some of them again in the next version.
Commit b8ce2d0df24

* Surface scheme: define phim and phih as functions. Add cap on phi for zeta > 1 to prevent crashes in stable conditions.
[Issue #28](https://github.com/dalesteam/dales/issues/28)
[Issue #30](https://github.com/dalesteam/dales/issues/30)
Commit 8469307

* Make sgs_surface_fix = .false. by default, and remove it from the namelist options.
This fix is not needed after the changes to the phi and psi functions in very stable conditions. It causes non-physical "flares" of high tke at stable conditions.

* Add tke field to the netCDF cross sections.
Having the TKE field available for inspection is good for testing the SGS and surface schemes.
Commit c9b43e4

* Allocate large arrays instead of storing on the stack.
There were only a few of these. Avoids problems with the stack size being too small -> crash.
Commit 1ab3ea1, 1950e32

* Store vertical netCDF cross sections only in the first tile in x and in y.
Commit 22f5643

* Expose more subgrid scheme parameters in namelist: ch1,ch2,cm,ce1,ce2. Add forgotten MPI_BCAST of courantp. We have been discussing whether the Prandtl number could be increased from ⅓ to 1 (Smagorinsky scheme), and correspondingly the ch2 coefficient could be decreased from 2 to 0 (SGS-TKE scheme). Adding ch2 and some other coefficients to the name list for easy testing of this.
Commit 13c2464

* Better restart logic: if trestart=0, write restart file only at the end. If < 0, don't write a restart file. Commit 8067376f9c

* Add macros to custom doxygen latex header.tex, to work with doxygen 1.8.2
Commit 1104f3b

* Add ECMWF build option. Cray fortran fix: avoid -1**real.
Commit 9fd1524

* Remove increment of uninitialized and unused ii variable, in the case of only one process. Tripped ifort's -check all.
Commit 2669b1d


### Optimizations

* Add more aggressive compiler options for gfortran, selectable with SYST=gnu-fast.
Enables -Ofast -march=native. Test on 1h of bomex case: gain 30% speed, results are close.
Also adding another system, SYST=ECMWF-intel, for using ifortran on the ECMWF Cray.
Commit 3434d3372c

* Add a null advection scheme. A way to use simpleice microphysics without advecting the unused scalar. Saves some percents of simulation time for free.
Commit b3d87987e4 or 826c7b5

* Optimization of advec_52, advec_62 - moving if(k==1) out of the k-loop, for a performance gain of about 8% (full simulation).
Commit 3b298b0d03, e18eb56163

* Remove rhoputin from advec_2nd, advec_52, advec_62, advec_5th.
rhoputin was a temporary array used to re-use the result of a multiplication.
~10% whole-simulation performance gain when removing this, due to fewer memory accesses.
Commit afd62e8, 6b9a3e41b7c

* Optimized advec_hybrid scheme. advec_hybrid.f90 is the original implementation (scheme 55), a new file advec_hybrid_f.f90 (scheme 555) has been added with a different implementation of the same advection scheme. The advection routine itself is about 4 times faster, a full simulation using it is about 30% faster. Commit e5e029bc4627a, 8d6fd70

* Modradfull optimization - calculate exp only when needed in coefft0.
Commit 5c712f3

* Optimization of the subgrid scheme save some divisions and intermediate arrays.
Commit 42897d8

* tstep_integrate: vector math instead of do loops faster, slightly shorter.
Commit ddcb29c


Known issues in DALES 4.2
--------------------------

* Runs with radiation do not restart properly.
[Issue #40](https://github.com/dalesteam/dales/issues/40)

* Radiation tendencies are not initialized to 0.

* Tendencies in the radiation statistics (output only) are sometimes calculated with the wrong sign.




Wishlist for future improvements
--------------------------------

* Make the MPI halo exchange less synchronous, to overlap communication and calculations.

* Permit longer time steps.
The CFL criterion is formulated for an equidistant grid, and unnecessarily restrictive.
adaptive time step: longer steps possible when dx,dy,dz non-equal, as in microHH.
see discussion for [issue #29](https://github.com/dalesteam/dales/issues/29).
When dx,dy,dz are non-equal, the diffusion criterion is stronger than strictly necessary. Relaxing it offsets the performance impact of commit a67cf8f. This has not been rigorously tested! As this permits longer time steps, it makes the simulation less stable. This may need compensating by specifying a slightly smaller peclet number e.g. in the namoptions file.
According to Chiel, microHH does the same.
(Commit 1151312)

* Don't do halo exchange of the m-fields
The m-fields ghost cells are seldom needed.  Saves maybe 1 % for a
single process run, and a large amount of MPI communication in a multi-process run.
This was tried but reverted since it interfered with chemistry, which *does* use the m-field halos.
(Commit 0c2cf59)

* Surface statistics 2D fields LWP, RWP, TWP and accumulated rain,
  output as netCDF in the cross-section routine.
  The accumulated rain works only for imicro=6.
  TODO:
    - make this optional, selectable by a flag in the namelist.
    - output the rain field only if imicro==6
    - check if duplicated in modlsmstat or modAGScross.

* Implementation of Mean-state acceleration [Jones et al, JAMES 2015,
https://doi.org/10.1002/2015MS000488]
A simple change in the time step routine, accelerating the rate of
change of the horizontal averages of prognostic quantities. Useful
for speeding up simulations *if* there is a large difference in time
scale for eddies and for horizontal averages.  The acceleration factor
can be set in namoptions. Default is 1, meaning no acceleration.
Huug proposes: add an option to enable acceleration for some time, for spinup.
Handle chemicals sv fields. Add namelist options, maybe a vector for all sv - accelerate or not.


* Merged netCDF writing. It is tedious to stitch together the separate netCDF files produced - one file per MPI task. There are scripts in the Dales repository, but they are from the time of slab parallelization and only stitch in one direction. Consider the parallel netCDF API where many MPI tasks can write together into a single file.

* Synchronize netCDF files periodically. When switching to compressed netCDF v4, the files can be unreadable while DALES is running or after DALES crashes or is stopped in the middle of the run.
