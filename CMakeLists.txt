INCLUDE(CMakeForceCompiler) #Necessary to change between MPI/parallel/profiler compilers, without having to do a make clean

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

if($ENV{SYST} STREQUAL "")
  set($ENV{SYST} STREQUAL "localpc_gfortran")
endif($ENV{SYST} STREQUAL "")

if($ENV{SYST} STREQUAL "localpc_gfortran")
  set(CMAKE_Fortran_COMPILER "mpif90")
  set(CMAKE_Fortran_FLAGS "-finit-real=nan -W -Wall -fdefault-real-8 -ffree-line-length-none -I/include -I/usr/include -I/usr/local/include" CACHE STRING "")
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3" CACHE STRING "")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fbounds-check -fbacktrace -fno-f2c -O0 -g -ffpe-trap=invalid,zero,overflow" CACHE STRING "")
elseif($ENV{SYST}  STREQUAL "HUYGENS")
  set(CMAKE_Fortran_COMPILER "mpfort")
  set(CMAKE_Fortran_FLAGS " -qfree=F90 -qrealsize=8  -qwarn64 -qflttrap=en:ov:zero:inv:imp -qflag=w:e -I$ENV{HOME}/include " CACHE STRING "")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O4 -qnoipa -qstrict=none:exceptions  " CACHE STRING "")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O2 -g -qfullpath -C -qflttrp=enable:nanq:overflow:zerodivide -qsigtrap -qinitauto=ff" CACHE STRING "")
elseif($ENV{SYST}  STREQUAL "localpc_ifort")
  set(CMAKE_Fortran_COMPILER "mpif90")
  set(CMAKE_Fortran_FLAGS "-r8 -ftz -fpe0 -extend_source" CACHE STRING "")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -pg" CACHE STRING "")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -ffpe-trap=invalid,zero,overflow" CACHE STRING "")
endif($ENV{SYST} STREQUAL "localpc_gfortran")
PROJECT(DALES Fortran)
cmake_minimum_required(VERSION 2.6)

# Set Profiler !JvdD from Mesohh
set (PROFILER CACHE STRING "Enable a profiler. Options: SCALASCA or MARMOT")

if(DEFINED PROFILER)
	if(PROFILER STREQUAL "SCALASCA")
		set(CMAKE_Fortran_COMPILER "/home/vddussen/bin/scalasca")
	elseif(PROFILER STREQUAL "MARMOT")
		set(PROFILE_COMPILER "marmotf90")
	endif(PROFILER STREQUAL "SCALASCA")
#	CMAKE_FORCE_Fortran_COMPILER(${PROFILE_COMPILER} mpfort
endif(DEFINED PROFILER)

MESSAGE(STATUS "Compiler " ${CMAKE_Fortran_COMPILER})

set(VERSION_MAJOR "3")
set(VERSION_MINOR "1")
set(VERSION_PATCH "90")
find_path(NETCDF_INCLUDE_DIR netcdf.mod NETCDF.mod
  DOC "NetCDF include directory (must contain netcdf.mod)"
  PATHS
  /usr/include
  $ENV{HOME}/include
  $ENV{SARA_NETCDF_INCLUDE}
PATH_SUFFIXES netcdf-3)
find_library(NETCDF_C_LIB netcdf
  DOC "NetCDF C library"
  PATHS
  /usr/lib
  /usr/lib64
  $ENV{HOME}/lib
  $ENV{HOME}/lib64
  $ENV{SARA_NETCDF_LIB}
)
find_library(NETCDF_FORTRAN_LIB netcdff
  DOC "NetCDF Fortran library"
  PATHS
  /usr/lib
  /usr/lib64
  $ENV{HOME}/lib
  $ENV{HOME}/lib64
  $ENV{SARA_NETCDF_LIB}
)

if(NETCDF_INCLUDE_DIR)
else(NETCDF_INCLUDE_DIR)
  MESSAGE(STATUS "WARNING: No NETCDF bindings are found.")
endif(NETCDF_INCLUDE_DIR)

if(NETCDF_C_LIB)
  set(NETCDF_LIBS ${NETCDF_C_LIB})
else(NETCDF_C_LIB)
  MESSAGE(STATUS "WARNING: No NETCDF bindings are found.")
endif(NETCDF_C_LIB)

if(NETCDF_FORTRAN_LIB)
  set(NETCDF_LIBS ${NETCDF_LIBS} ${NETCDF_FORTRAN_LIB})
else(NETCDF_FORTRAN_LIB)
  MESSAGE(STATUS "WARNING: No Fortran NETCDF bindings are found.")
endif(NETCDF_FORTRAN_LIB)

### Find HDF5 directory
#find_path(HDF5_INCLUDE_DIR hdf5.mod HDF5.mod
#  DOC "HDF5 include directory (must contain hdf5.mod)"
#  PATHS
#  /usr/include
#  $ENV{HOME}/include
#  $ENV{SARA_HDF5_INCLUDE}
#PATH_SUFFIXES hdf5)

#find_library(HDF5_C_LIB hdf5_hl
#  DOC "HDF5 C library"
#  PATHS
#  /usr/lib
#  /usr/lib64
#  $ENV{HOME}/lib
#  $ENV{HOME}/lib64
#  $ENV{SARA_HDF5_LIB}
#)

#find_library(HDF5_FORTRAN_LIB hdf5
#  DOC "HDF5 Fortran library"
#  PATHS
#  /usr/lib
#  /usr/lib64
#  $ENV{HOME}/lib
#  $ENV{HOME}/lib64
#  $ENV{SARA_HDF5_LIB}
#)
#find_library(HDF5_Z_LIB z
#  DOC "HDF5 z library"
#  PATHS
#  /usr/lib
#  /usr/lib64
#  $ENV{HOME}/lib
#  $ENV{HOME}/lib64
#  $ENV{SARA_HDF5_LIB}
#)
#find_library(HDF5_GPFS_LIB gpfs
#  DOC "HDF5 gpfs library"
#  PATHS
#  /usr/lib
#  /usr/lib64
#  $ENV{HOME}/lib
#  $ENV{HOME}/lib64
#  $ENV{SARA_HDF5_LIB}
#)
### Includes
#if(HDF5_INCLUDE_DIR)
#else(HDF5_INCLUDE_DIR)
#  MESSAGE(STATUS "WARNING: No HDF5 bindings are found.")
#endif(HDF5_INCLUDE_DIR)
#
#if(HDF5_C_LIB)
#  set(HDF5_LIBS ${HDF5_C_LIB})
#else(HDF5_C_LIB)
#  MESSAGE(STATUS "WARNING: No C HDF5 bindings are found.")
#endif(HDF5_C_LIB)
#
#if(HDF5_FORTRAN_LIB)
#  set(HDF5_LIBS ${HDF5_LIBS} ${HDF5_FORTRAN_LIB})
#else(HDF5_FORTRAN_LIB)
#  MESSAGE(STATUS "WARNING: No Fortran HDF5 bindings are found.")
#endif(HDF5_FORTRAN_LIB)
#
#if(HDF5_Z_LIB)
#  set(HDF5_LIBS ${HDF5_LIBS} ${HDF5_Z_LIB})
#else(HDF5_Z_LIB)
#  MESSAGE(STATUS "WARNING: No z HDF5 bindings are found.")
#endif(HDF5_Z_LIB)
#
#if(HDF5_GPFS_LIB)
#  set(HDF5_LIBS ${HDF5_LIBS} ${HDF5_GPFS_LIB})
#else(HDF5_GPFS_LIB)
#  MESSAGE(STATUS "WARNING: No Gpfs HDF5 bindings are found.")
#endif(HDF5_GPFS_LIB)

###documentation
INCLUDE(FindDoxygen)
IF (DOXYGEN)
    ADD_SUBDIRECTORY(utils/doc)
ELSE (DOXYGEN)
    MESSAGE(STATUS "WARNING: Doxygen not found - Reference manual will
not be created")
ENDIF (DOXYGEN)

if (CASE)
    FILE(GLOB usrfile "${CMAKE_CURRENT_SOURCE_DIR}/cases/${CASE}/moduser.f90")
    FILE(REMOVE ${CMAKE_CURRENT_SOURCE_DIR}/src/moduser.f90)
    if (usrfile STREQUAL "")
      set(usrfile "${CMAKE_CURRENT_SOURCE_DIR}/cases/standard/moduser.f90")
    endif (usrfile STREQUAL "")
      CONFIGURE_FILE(${usrfile} ${CMAKE_CURRENT_SOURCE_DIR}/src/moduser.f90 COPYONLY)
      MESSAGE(STATUS "Case " ${CASE} " uses " ${usrfile})
else(CASE)
   FILE(GLOB usrfile "${CMAKE_CURRENT_SOURCE_DIR}/src/moduser.f90")
   if (usrfile STREQUAL "")
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/cases/standard/moduser.f90 ${CMAKE_CURRENT_SOURCE_DIR}/src/moduser.f90 COPYONLY)
    endif (usrfile STREQUAL "")
endif(CASE)
#unset(CASE CACHE)


ADD_SUBDIRECTORY( src)
