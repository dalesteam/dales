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
  set(CMAKE_Fortran_FLAGS "-finit-real=nan -fdefault-real-8 -ffree-line-length-none -I /usr/lib64/gfortran/modules/mpich/" CACHE STRING "")
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3" CACHE STRING "")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fbounds-check -fbacktrace -fno-f2c -O0 -g -ffpe-trap=invalid,zero,overflow" CACHE STRING "")
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
elseif("$ENV{SYST}" STREQUAL "gnu-fast")

elseif("$ENV{SYST}" STREQUAL "LMU_MIM")
  set(CMAKE_Fortran_COMPILER "mpiifort")
  set(CMAKE_Fortran_FLAGS "-cpp -r8 -ftz -extend_source -g -sox -no-wrap-margin" CACHE STRING "")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3" CACHE STRING "")
  set(CMAKE_Fortran_FLAGS_DEBUG "-traceback -fpe1 -O0 -g -check all" CACHE STRING "")

else()

  set(CMAKE_Fortran_COMPILER "mpif90")
  set(CMAKE_Fortran_FLAGS "-finit-real=nan -W -Wall -fdefault-real-8 -ffree-line-length-none" CACHE STRING "")
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -Ofast -march=native -g -fbacktrace" CACHE STRING "")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fbounds-check -fbacktrace -fno-f2c -O0 -g -ffpe-trap=invalid,zero,overflow" CACHE STRING "")


endif()

## Project parameters
PROJECT(DALES Fortran)
cmake_minimum_required(VERSION 2.6)
set(VERSION_MAJOR "4")
set(VERSION_MINOR "2")
set(VERSION_PATCH "0")

### If necessary, resort to BASH-methods to find netcdf-directory
EXEC_PROGRAM(${CMAKE_CURRENT_SOURCE_DIR}/findnetcdf OUTPUT_VARIABLE ADDMODULEPATH)

### Find NetCDF files
FIND_PATH(NETCDF_INCLUDE_DIR netcdf.mod NETCDF.mod
  PATHS
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

if(NETCDF_INCLUDE_DIR)
  include_directories(${NETCDF_INCLUDE_DIR})
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

option(WITH_TENSTREAM "Include the tenstream solver for radiation" OFF)
if(WITH_TENSTREAM)
  if(  (NOT TENSTREAM_DIR)
    OR (NOT TENSTREAM_BUILD_DIR) )
    message( STATUS "You asked to include the TENSTREAM solver, yet you did not specify a path to the library -- please set:")
    message( STATUS " TENSTREAM_DIR         (root dir of tenstream solver) ")
    message( STATUS " TENSTREAM_BUILD_DIR   (directory where you issued cmake&&make)" )
    set(TENSTREAM_DIR       "" CACHE STRING "Tenstream base dir"            )
    set(TENSTREAM_BUILD_DIR "" CACHE STRING "Tenstream build or install dir")
    message( FATAL_ERROR " Tenstream Directory not set!" )
  endif((NOT TENSTREAM_DIR) OR (NOT TENSTREAM_BUILD_DIR) )

  set(CMAKE_C_FLAGS               ${USER_C_FLAGS} )
  set(CMAKE_Fortran_FLAGS         ${USER_Fortran_FLAGS} )
  set(CMAKE_Fortran_FLAGS_RELEASE ${USER_Fortran_FLAGS_RELEASE} )
  set(CMAKE_Fortran_FLAGS_DEBUG   ${USER_Fortran_FLAGS_DEBUG})

  # Add Tenstream
  include_directories("${TENSTREAM_BUILD_DIR}/include/")

  set (TENSTREAM_LIBS "pprts_rrtmg" "pprts_rrtm_lw" "pprts_rrtm_sw" "pprts" "pprts_base")
  foreach (LIB ${TENSTREAM_LIBS})
    FIND_LIBRARY(TENSTREAM_LIBRARY_${LIB} "${LIB}" "${TENSTREAM_BUILD_DIR}/lib/")
    list(APPEND LIBS "${TENSTREAM_LIBRARY_${LIB}}")
  endforeach(LIB)
  message( STATUS "Adding TenStream libs: ${LIBS}")

  # Add Petsc
  enable_language(C)
  list(APPEND CMAKE_MODULE_PATH "${TENSTREAM_DIR}/cmake/petsc_cmake_modules/")
  MESSAGE(STATUS "Using PETSC $ENV{PETSC_DIR} / $ENV{PETSC_ARCH} ")
  find_package(PETSc)
  message(STATUS "Found Petsc:    ${PETSC_FOUND}" )
  message(STATUS "Petsc Includes: ${PETSC_INCLUDES}" )
  message(STATUS "Petsc Library : ${PETSC_LIBRARIES}")
  message(STATUS "Petsc compiler: ${PETSC_COMPILER}" )
  message(STATUS "Petsc VERSION : ${PETSC_VERSION}" )
  message(STATUS "Petsc definiti: ${PETSC_DEFINITIONS}")
  message(STATUS "Petsc variables ${petsc_conf_variables}")
  message(STATUS "Petsc variables ${petsc_conf_variables_arch}")

  include_directories(${PETSC_INCLUDES})
  list(APPEND LIBS "${PETSC_LIBRARIES}")

  add_definitions(-DHAVE_TENSTREAM)

  file(STRINGS ${petsc_conf_variables_arch} PETSC_VARIABLES_STRINGS)
  string(REGEX REPLACE ".*PETSC_PRECISION = ([a-z]+).*" "\\1"
    PETSC_PRECISION ${PETSC_VARIABLES_STRINGS})
  message(STATUS " Precision is :: ${PETSC_PRECISION} ")

  string(COMPARE EQUAL ${PETSC_PRECISION} "single" _is_single)

  if(${_is_single})
    add_definitions(-DTENSTREAM_SINGLE)
  endif(${_is_single})
endif(WITH_TENSTREAM)


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
ADD_DEFINITIONS(-DHAVE_3DHEATING)
ADD_SUBDIRECTORY(src)
