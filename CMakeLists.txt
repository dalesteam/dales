############################
# DALES CMAKE FILE
# Supports (thusfar) GNU, Intel and IBM compilers
# To set the compiler to a specific one, set the FC environment variable
# NETCDF can be set in the config/default.cmake, or through command line options
# NETCDF_INCLUDE, NETCDF_LIB, FFTW_INCLUDE and FFTW_LIB, respectively.
# Command line options are fed to CMake using the -D switch
# Options include:
# * MPI (TRUE/FALSE) to toggle the use of MPI. True by default.
# * PROFILER (SCALASCA/MARMOT) to enable a specific profiler. Default is none.
# * CMAKE_BUILD_TYPE (RELEASE/DEBUG) to build in optimized or debug mode. Default is RELEASE
# All options are persistent in the sense that once set, they will last until changed by the user.
# For daily use "cmake .. && make" should suffice to build the code, but a command line could look like:
# export FC="gfortran" && cmake -D SYST=tornado -D MPI=TRUE -D PROFILER=SCALASCA -D CMAKE_BUILD_TYPE=DEBUG .. && make
############################
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/config)

INCLUDE(CMakeForceCompiler) #Necessary to change between MPI/parallel/profiler compilers, without having to do a make clean
cmake_minimum_required (VERSION 2.8.1)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})


#  make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release." FORCE )
else()
  string(TOUPPER ${CMAKE_BUILD_TYPE} TEMP)
  set(CMAKE_BUILD_TYPE ${TEMP} CACHE STRING
      "Choose the type of build, options are: None Debug Release."  FORCE)
endif ()

#Load system specific settings
if (NOT SYST)
  set (SYST default CACHE STRING
      "Choose the location: mpipc thunder fedora archlinux supermuc jugene" FORCE )
endif ()
include (${SYST} OPTIONAL)

#Start the project only after all the variables are set
project (dales Fortran)

get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER_ID} NAME)

if (NOT ${Fortran_COMPILER_WRAPPER}  STREQUAL "")
  CMAKE_FORCE_Fortran_COMPILER(${Fortran_COMPILER_WRAPPER} ${Fortran_COMPILER_NAME})
endif ()
MESSAGE(STATUS "Fortran Compiler " ${CMAKE_Fortran_COMPILER})

set(CMAKE_Fortran_FLAGS ${USER_Fortran_FLAGS} )
set(CMAKE_Fortran_FLAGS_RELEASE ${USER_Fortran_FLAGS_RELEASE} )
set(CMAKE_Fortran_FLAGS_DEBUG ${USER_Fortran_FLAGS_DEBUG})
MESSAGE(STATUS "Build Type " ${CMAKE_BUILD_TYPE})
if (CMAKE_BUILD_TYPE STREQUAL "RELEASE")
  MESSAGE(STATUS "Compiler Flags " ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE})
else()
  MESSAGE(STATUS "Compiler Flags " ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_DEBUG})
endif()


set(INCLUDE_DIRS ${NETCDF_INCLUDE_DIR})

add_subdirectory(src)


#####################
# DOCUMENTATION
#####################
add_custom_target(todo ALL)
ADD_CUSTOM_COMMAND(TARGET todo POST_BUILD
                  COMMAND echo "UCLALES TODO LIST" > TODO
                  COMMAND date  >> TODO
                  COMMAND grep -Rin \\todo  src | sed 's/!.*TODO//I' >>  TODO
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                  )

