### Choose CMAKE Type
if(NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif()

### Set compiler flags
if("$ENV{SYST}" STREQUAL "HUYGENS")
  set(CMAKE_Fortran_COMPILER "mpfort")
  set(CMAKE_Fortran_FLAGS "-qfree=F90 -qrealsize=8  -qwarn64 -qflttrap=en:ov:zero:inv:imp -qflag=w:e" CACHE STRING "")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O4 -qnoipa -qstrict=none:exceptions" CACHE STRING "")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O2 -g -qfullpath -C -qflttrp=enable:nanq:overflow:zerodivide -qsigtrap -qinitauto=ff" CACHE STRING "")
elseif("$ENV{SYST}" STREQUAL "CARTESIUS")
  set(CMAKE_Fortran_COMPILER "mpiifort")
  set(CMAKE_Fortran_FLAGS "-r8 -ftz -extend_source" CACHE STRING "")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3" CACHE STRING "")
  set(CMAKE_Fortran_FLAGS_DEBUG "-traceback -fpe1 -O0 -g -check all" CACHE STRING "")
elseif("$ENV{SYST}" STREQUAL "localpc_ifort")
  set(CMAKE_Fortran_COMPILER "mpif90")
  set(CMAKE_Fortran_FLAGS "-r8 -ftz -extend_source" CACHE STRING "")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3" CACHE STRING "")
  set(CMAKE_Fortran_FLAGS_DEBUG "-traceback -fpe1 -O0 -g -check all" CACHE STRING "")
elseif("$ENV{SYST}" STREQUAL "HYDRA")
  set(CMAKE_Fortran_COMPILER "mpiifort")
  set(CMAKE_Fortran_FLAGS "-r8 -ftz -extend_source" CACHE STRING "")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3" CACHE STRING "")
  set(CMAKE_Fortran_FLAGS_DEBUG "-traceback -fpe1 -O0 -g -check all" CACHE STRING "")
elseif("$ENV{SYST}" STREQUAL "FEDORA")
  set(CMAKE_Fortran_COMPILER "mpif90")
  set(CMAKE_Fortran_FLAGS "-fdefault-real-8 -fdefault-double-8 -march=native -malign-double -ffree-line-length-none -I /usr/lib64/gfortran/modules/mpich/ -std=legacy" CACHE STRING "")
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3" CACHE STRING "")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-finit-real=nan -fbounds-check -fbacktrace -fno-f2c -O0 -g -ffpe-trap=invalid,zero,overflow" CACHE STRING "")
  set (CMAKE_Fortran_FLAGS_DEV     "-fbounds-check -fbacktrace -fno-f2c -O3 -g -ffpe-trap=invalid,zero,overflow" CACHE STRING "")
elseif("$ENV{SYST}" STREQUAL "ECMWF")
  set(CMAKE_Fortran_COMPILER "ftn")
  set(CMAKE_Fortran_FLAGS "-s real64" CACHE STRING "")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3" CACHE STRING "")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -R b -K trap=fp" CACHE STRING "")
elseif("$ENV{SYST}" STREQUAL "ECMWF-intel")
  set(CMAKE_Fortran_COMPILER "ftn")
  set(CMAKE_Fortran_FLAGS "-r8 -ftz -extend_source" CACHE STRING "")
  set(CMAKE_Fortran_FLAGS_RELEASE "-g -traceback -O3 -xHost" CACHE STRING "")
  set(CMAKE_Fortran_FLAGS_DEBUG "-traceback -fpe1 -O0 -g -check all" CACHE STRING "")
elseif("$ENV{SYST}" STREQUAL "ECMWF-gnu")
  set(CMAKE_Fortran_COMPILER "ftn")
  set(CMAKE_Fortran_FLAGS "-W -Wall -fdefault-real-8 -fdefault-double-8 -march=native -malign-double -ffree-line-length-none" CACHE STRING "")
  set(CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -Ofast -g -fbacktrace" CACHE STRING "")
  set(CMAKE_Fortran_FLAGS_DEBUG   "-finit-real=nan -fbounds-check -fbacktrace -fno-f2c -O0 -g -ffpe-trap=invalid,zero,overflow" CACHE STRING "")
elseif("$ENV{SYST}" STREQUAL "gnu-fast")
  set(CMAKE_Fortran_COMPILER "mpif90")
  set(CMAKE_Fortran_FLAGS "-W -Wall -fdefault-real-8 -fdefault-double-8 -march=native -malign-double -ffree-line-length-none -std=legacy" CACHE STRING "")
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -Ofast -g -fbacktrace" CACHE STRING "")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-finit-real=nan -fbounds-check -fbacktrace -fno-f2c -O0 -g -ffpe-trap=invalid,zero,overflow" CACHE STRING "")
elseif("$ENV{SYST}" STREQUAL "lisa-intel")
   set(CMAKE_Fortran_COMPILER "mpiifort")
   set(CMAKE_Fortran_FLAGS "-r8 -ftz -extend_source" CACHE STRING "")
   set(CMAKE_Fortran_FLAGS_RELEASE "-g -traceback -O3 -xHost" CACHE STRING "")
   set(CMAKE_Fortran_FLAGS_DEBUG "-traceback -fpe1 -O0 -g -check all" CACHE STRING "")
elseif("$ENV{SYST}" STREQUAL "FX-Fujitsu")
   set(CMAKE_Fortran_COMPILER "mpifrtpx")
   set(CMAKE_Fortran_FLAGS "-Ad -g" CACHE STRING "")
   set (CMAKE_Fortran_FLAGS_RELEASE "-Kfast" CACHE STRING "")
   set (CMAKE_Fortran_FLAGS_DEBUG   "-Haefosux -O2" CACHE STRING "")
   # Fujitsu fortran compiler
   # -Ad default real is 64 bit
   # -Kfast optimization and SVE vectorization
   # -Haefosux checking: s-subscript range, u-undefinied values
   #                     a-argument mismatches, e-shape conformance
   #                     f-simultaneousOPEN and I/O recursive calls
   #                     o-overlapping dummy arguments and extend undefined
   #                     x-module, submodule, common, and pointer undefined
elseif("$ENV{SYST}" STREQUAL "NO_OVERRIDES")
   # don't set any compilation flags here, allows setting them outside CMake
else()
  set(CMAKE_Fortran_COMPILER "mpif90")
  set(CMAKE_Fortran_FLAGS "-fdefault-real-8 -fdefault-double-8 -march=native -malign-double -ffree-line-length-none -std=legacy" CACHE STRING "")
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3" CACHE STRING "")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-finit-real=nan -fbounds-check -fbacktrace -fno-f2c -O0 -g -ffpe-trap=invalid,zero,overflow" CACHE STRING "")
endif()

# On compiler options
# June 2020: added -std=legacy to gfortran for compatibility with gfortran 10, which introduced stricter type checking
# the old FFT code from Netlib/FFTPACK does not compile with gfortran 10 without this flag.

# For GCC based compilers:
#  * -march=native   : compile for the current CPU, allows use of modern AVX etc.
#  * -malign-double  : align stuff on 64-bit boundaries.
#  * -finit-real-nan : makes the code slower, moved them to debug


## Project parameters
PROJECT(DALES Fortran)
cmake_minimum_required(VERSION 3.2)
set(VERSION_MAJOR "4")
set(VERSION_MINOR "2")
set(VERSION_PATCH "1")

# use nf-config to find paths to netcdf. Broken on Fedora due to module path missing.
execute_process(COMMAND nf-config --prefix OUTPUT_VARIABLE NFCONFIG_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND nf-config --includedir OUTPUT_VARIABLE NFCONFIG_INCLUDE_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
# execute_process(COMMAND nf-config --flibs OUTPUT_VARIABLE NFCONFIG_FORTRAN_LIB OUTPUT_STRIP_TRAILING_WHITESPACE)
# not used at the moment - not just a path but several compiler flags

### If necessary, resort to BASH-methods to find netcdf-directory
EXEC_PROGRAM(${CMAKE_CURRENT_SOURCE_DIR}/findnetcdf OUTPUT_VARIABLE ADDMODULEPATH)

### Find NetCDF files
FIND_PATH(NETCDF_INCLUDE_DIR netcdf.mod NETCDF.mod
  PATHS
  ${NFCONFIG_INCLUDE_DIR}
  $ENV{SARA_NETCDF_INCLUDE}
  $ENV{SURFSARA_NETCDF_INCLUDE}
  $ENV{NETCDF_INCLUDE}
  ${ADDMODULEPATH}/include
  /usr/include
  $ENV{HOME}/include
  /usr/lib64/gfortran/modules
  DOC "NetCDF include directory (must contain netcdf.mod)"
)

FIND_LIBRARY(NETCDF_C_LIB netcdf
  PATHS
  $ENV{SARA_NETCDF_LIB}
  $ENV{SURFSARA_NETCDF_LIB}
  $ENV{NETCDF_LIB}
  ${ADDMODULEPATH}/lib
  ${ADDMODULEPATH}/lib64
  /usr/lib
  /usr/lib64
  $ENV{HOME}/lib
  $ENV{HOME}/lib64
  DOC "NetCDF C library"
)

FIND_LIBRARY(NETCDF_FORTRAN_LIB netcdff
  PATHS
  ${NFCONFIG_PREFIX}/lib
  $ENV{SARA_NETCDF_LIB}
  $ENV{SURFSARA_NETCDF_LIB}
  $ENV{NETCDF_LIB}
  ${ADDMODULEPATH}/lib
  ${ADDMODULEPATH}/lib64
  /usr/lib
  /usr/lib64
  $ENV{HOME}/lib
  $ENV{HOME}/lib64
  DOC "NetCDF Fortran library"
)

# libHYPRE.a for static library
FIND_LIBRARY(HYPRE_LIB HYPRE
  PATHS
  $ENV{HYPRE_LIB}
  ${ADDMODULEPATH}/lib
  ${ADDMODULEPATH}/lib64
  /usr/lib
  /usr/lib64
  $ENV{HOME}/lib
  $ENV{HOME}/lib64
  DOC "Iterative solver library"
)

FIND_PATH(FFTW_INCLUDE_DIR fftw3.f03
  PATHS
  $ENV{FFTW_INCLUDE_DIR}
  /usr/include
  DOC "FFTW include files (fftw3.f03)"
)

#libfftw3.a for static library
FIND_LIBRARY(FFTW_LIB fftw3
  PATHS
  $ENV{FFTW_LIB}
  /usr/lib
  /usr/lib64
  DOC "FFTW static library"
)

if(NETCDF_INCLUDE_DIR)
  include_directories(${NETCDF_INCLUDE_DIR})
else(NETCDF_INCLUDE_DIR)
  MESSAGE(STATUS "WARNING: No NETCDF bindings are found.")
endif(NETCDF_INCLUDE_DIR)

if(FFTW_INCLUDE_DIR)
  include_directories(${FFTW_INCLUDE_DIR})
else(FFTW_INCLUDE_DIR)
  MESSAGE(STATUS "WARNING: No FFTW3 bindings are found (fftw3.f03).")
endif(FFTW_INCLUDE_DIR)

if(NETCDF_C_LIB)
  set(NETCDF_LIBS ${NETCDF_C_LIB})
else(NETCDF_C_LIB)
  MESSAGE(STATUS "WARNING: No NETCDF C bindings are found - may be OK.")
endif(NETCDF_C_LIB)

if(NETCDF_FORTRAN_LIB)
  set(NETCDF_LIBS ${NETCDF_LIBS} ${NETCDF_FORTRAN_LIB})
else(NETCDF_FORTRAN_LIB)
  MESSAGE(STATUS "WARNING: No Fortran NETCDF bindings are found.")
endif(NETCDF_FORTRAN_LIB)

### Iterative solver (HYPRE)
OPTION(USE_HYPRE "Also build iterative solver (optional, needs HYPRE)" OFF)

### FFTW based poisson solver (FFTW)
OPTION(USE_FFTW "Also build FFTW based poisson solver (optional, needs FFTW3)" OFF)

### Documentation
INCLUDE(FindDoxygen)
if(DOXYGEN)
    ADD_SUBDIRECTORY(utils/doc)
else()
    MESSAGE(STATUS "WARNING: Doxygen not found - Reference manual will not be created")
endif()

### Set case
if(NOT CASE)
  set (CASE standard CACHE STRING
      "Set the case."
      FORCE)
endif()

### Add case specific file
FILE(GLOB usrfile "${CMAKE_CURRENT_SOURCE_DIR}/cases/${CASE}/moduser.f90")
if(usrfile STREQUAL "")
  set(usrfile "${CMAKE_CURRENT_SOURCE_DIR}/cases/standard/moduser.f90")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${usrfile} ${CMAKE_CURRENT_SOURCE_DIR}/src/moduser.f90)
MESSAGE(STATUS "Case " ${CASE} " uses " ${usrfile})

ADD_SUBDIRECTORY(src)
