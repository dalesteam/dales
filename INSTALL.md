Installing DALES
================

These are installation instruuctions for DALES. See also the [installation notes](https://github.com/dalesteam/dales/wiki/Installation-notes) in the DALES wiki for further notes on specific systems.

DEPENDENCIES
------------

Dependencies of DALES are
1) MPI 
2) NetCDF (optional)
3) CMake
4) FFTW version 3.x (optional) 
5) HYPRE (optional) for the iterative Poisson solver.
6) Doxygen, DOT and LaTeX for building the documentation

Compiling
---------

Packages are available at <https://github.com/dalesteam/dales>
If git is installed, the best way to obtain DALES is with:

`git clone https://github.com/dalesteam/dales.git`
For updates you can then just use the command 
```
git pull origin master
```

Then follow the following instructions:
```
cd dales 
# git checkout <branch>  optionally switch to another branch
mkdir build 
cd build 
cmake ..
make
# make install # optional
```
This should provide the dales4 binary in the build/src directory. New builds can be made by `make` in the build directory. The documentation is built by typing `make docs` in the build directory.

Compiling with HYPRE or FFTW
----------------------------

To compile DALES with support for HYPRE and/or FFTW, add the following flags on the cmake command line.
For HYPRE:
```
cmake .. -DUSE_HYPRE=True -DHYPRE_LIB=path-to-libHYPRE.so
```
where -DHYPRE_LIB is needed if the library is not automatically found.
For FFTW:
```
cmake .. -DUSE_FFTW=True -DFFTW_LIB=path-to-libfftw3.so
```
where the path again is needed only if the library is not automatically found.

Even when compiled with these libraries, DALES will use the traditional FFT-based Poisson solver unless
one of the alternatives are specified in the namelist, e.g. to use FFTW:

```
&SOLVER
solver_id = 100
/
```
See the [Wiki](https://github.com/dalesteam/dales/wiki/Alternative-Poisson-solvers) for the options for the iterative solver.

Optimization and tuning
-----------------------

Before `cmake`, one can execute
```
export SYST=gnu-fast
```
which triggers compilation with more aggressive optimization options, in particular `-Ofast` and `-march=native`. Other pre-defined sets of compilation options are included in the file `CMakeLists.txt`.

Choosing debug or release builds
--------------------------------

One can change between release and debug options by invoking commandline options of cmake:
```
cmake -DCMAKE_BUILD_TYPE="DEBUG" ..
```
or `"RELEASE"`, if that is the build-type desired.  `"RELEASE"` is the default.
Various case specific routines can be found in de dales/case/<casename> and are built by invoking
```
cmake -DCASE="<casename>"
```
where `<casename>` could be something like "rico". Regardless of the case built, the standard DALES runs are always possible (and default).





