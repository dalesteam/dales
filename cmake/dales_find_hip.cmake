macro( dales_find_hip )
  if( NOT DEFINED HIPFORT_INCLUDE_DIRS )
    find_path( HIPFORT_INCLUDE_DIRS
               NAMES hipfort.mod HIPFORT.mod
               HINTS ENV ROCM_PATH
                     ENV ROCM_ROOT
                     ENV ROCM_HOME
                     ENV EBROOTHIPFORT
               PATH_SUFFIXES include/hipfort/amdgcn
               REQUIRED )
  endif()

  if( NOT DEFINED HIPFORT_LIB )
    find_library( HIPFORT_LIB libhipfort-amdgcn.a
                  HINTS ENV ROCM_PATH
                        ENV ROCM_ROOT
                        ENV ROCM_HOME
                        ENV EBROOTHIPFORT
                  PATH_SUFFIXES lib/
                  REQUIRED )
  endif()

  if( NOT DEFINED HIPFFT_LIB )
    find_library( HIPFFT_LIB libhipfft.so
                  HINTS ENV ROCM_PATH
                        ENV ROCM_ROOT
                        ENV ROCM_HOME
                        ENV EBROOTHIPFORT
                  PATH_SUFFIXES lib/
                  REQUIRED )
  endif()

  if( HIPFORT_INCLUDE_DIRS AND HIPFORT_LIB )
    set( HAVE_HIP ON )
    add_compile_definitions( USE_HIP )
    ecbuild_info( "Found hipfort: [${HIPFORT_LIB}]" )
  else()
    ecbuild_info( "Could not find hipfort" )
  endif()
endmacro()
