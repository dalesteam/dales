cmake_minimum_required( VERSION 3.17 FATAL_ERROR )
list( APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/external/ecbuild/cmake" )
find_package( ecbuild REQUIRED )

### Set default CMAKE build type: Release
if(NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif()

project( DALES VERSION 5.0.0 LANGUAGES Fortran )

# Disable Intel's MKL FFTW
set( FFTW_ENABLE_MKL OFF )

# Options

ecbuild_add_option( FEATURE FFTW
                    DESCRIPTION "Build with FFTW"
                    DEFAULT ON
                    REQUIRED_PACKAGES "NAME FFTW COMPONENTS single double" )

ecbuild_add_option( FEATURE HYPRE
                    DESCRIPTION "Build with HYPRE"
                    DEFAULT OFF )

ecbuild_add_option( FEATURE FP32_FIELDS
                    DESCRIPTION "Use single precision for fields"
                    DEFAULT OFF )

ecbuild_add_option( FEATURE FP32_POIS
                    DESCRIPTION "Use single precision for Poisson solver"
                    DEFAULT OFF )

ecbuild_add_option( FEATURE FP32_RAD
                    DESCRIPTION "Use single precision for RTE-RRTMGP (use with caution!)"
                    DEFAULT OFF )

ecbuild_add_option( FEATURE DALESLIB
                    DESCRIPTION "Build DALES libraries for use with OMUSE"
                    DEFAULT OFF )

ecbuild_add_option( FEATURE ACC
                    DESCRIPTION "Build with GPU support through OpenACC"
                    DEFAULT OFF )

ecbuild_add_option( FEATURE RRTMG
                    DESCRIPTION "Build with legacy RRTMG libraries"
                    ADVANCED DEFAULT ON )

# Check for required dependencies

ecbuild_find_package( NetCDF REQUIRED COMPONENTS Fortran C )

ecbuild_find_package( MPI REQUIRED COMPONENTS Fortran )

# Set precision

if( ENABLE_FP32_FIELDS )
  set ( FIELD_PRECISION 32 )
else()
  set ( FIELD_PRECISION 64 )
endif()

add_compile_definitions( FIELD_PRECISION=${FIELD_PRECISION} )

if( ENABLE_FP32_POIS )
  set ( POIS_PRECISION 32 )
else()
  set ( POIS_PRECISION 64 )
endif()

add_compile_definitions( POIS_PRECISION=${POIS_PRECISION} )

if( ENABLE_FP32_RAD )
  add_compile_definitions( RTE_USE_SP )
  set( ENABLE_RRTMG OFF )
  ecbuild_info( "RRTMG is disabled because single precision is selected for RTE-RRTMGP" )
endif()

# Additional options

if( ENABLE_FFTW )
  add_compile_definitions( USE_FFTW )
endif()

if( ENABLE_HYPRE )
  # For some reason, ecbuild_find_package doesn't find HYPRE (on Mac, at least)
  find_library( HYPRE_LIB HYPRE REQUIRED )
  ecbuild_info( "Found HYPRE: ${HYPRE_LIB}" )
  add_compile_definitions( USE_HYPRE )
endif()

if( ENABLE_RRTMG )
  add_compile_definitions( USE_RRTMG )
endif()

# GPU stuff
include( dales_find_cuda )
include( dales_find_hip )

if( ENABLE_ACC )
  find_package( OpenACC COMPONENTS Fortran QUIET )
  if( OpenACC_Fortran_FOUND )
    set( OpenACC ON )
    ecbuild_info( "Found OpenACC" )
  endif()

  # Try to find CUDA, if it's not available try HIP instead
  if( ENABLE_ACC )
    dales_find_cuda()
    if( NOT HAVE_CUDA )
      dales_find_hip()
    endif()
  endif()

  # For debug builds, enable accelerator info
  if( ${CMAKE_Fortran_COMPILER_ID} MATCHES PGI|NVHPC )
    ecbuild_add_fortran_flags( "-Minfo=accel" BUILD DEBUG ) 
  endif()

  add_compile_definitions( DALES_GPU )

  # Disable automatic asynchronious kernel execution on Cray
  if( ${CMAKE_Fortran_COMPILER_ID} MATCHES Cray )
    ecbuild_add_fortran_flags( "-h acc_model=auto_async_none" )
  endif()
endif()

# Find the fortran-support library (included).
include ( dales_add_fortran_support )

# Set the compiler flags for the selected compiler
include( dales_compile_options )

# Add RRTMGP
include( dales_add_rrtmgp )

add_subdirectory( src )
