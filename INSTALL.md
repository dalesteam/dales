Installing DALES
================

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
mkdir build 
cd build 
cmake ..
make
# make install # optional
```
This should provide the dales4 binary in the build/src directory. New builds can be made by `make` in the build directory. The documentation is built by typing `make docs` in the build directory.

Optimization and tuning
-----------------------

Before `cmake`, one can execute
```
export SYST=gnu-fast
```
which triggers compilation with more aggressive optimization options, in particular `-Ofast` and `-march=native`. Other pre-defined sets of compilation options are included in the file `CMakeLists.txt`.

ALTERNATIVE BUILDS
==================

One can change between release and debug options by invoking commandline options of cmake:
```
cmake -DCMAKE_BUILD_TYPE="DEBUG" ..
```
or `"RELEASE"`, if that is the build-type desired.
Various case specific routines can be found in de dales/case/<casename> and are built by invoking
```
cmake -DCASE="<casename>"
```
where `<casename>` could be something like "rico". Regardless of the case built, the standard DALES runs are always possible (and default).


DEPENDENCIES
============

Dependencies of DALES are
1) MPI (the only obligatory one)
2) NetCDF (optional)
3) CMake (alternatively, the old makefile icw makedepf90 still works)
4) Doxygen, DOT and LaTeX for building the documentation

See also the [installation notes](https://github.com/dalesteam/dales/wiki/Installation-notes) in the DALES wiki for further notes on specific systems.