macro( dales_find_hip )
  if( NOT DEFINED HIPFORT_INCLUDE_DIRS )
    find_path( HIPFORT_INCLUDE_DIRS HIPFORT.mod
               HINTS ENV EBROOTHIPFORT
               PATH_SUFFIXES include/hipfort/amdgcn
               REQUIRED )
  endif()
  
  if( NOT DEFINED HIPFORT_LIB )
    find_library( HIPFORT_LIB libhipfort-amdgcn.a
                  HINTS ENV EBROOTHIPFORT
        	        PATH_SUFFIXES lib/
        	        REQUIRED )
  endif()

  if( HIPFORT_INCLUDE_DIRS AND HIPFORT_LIB )
    set( HAVE_HIP ON )
    ecbuild_info( "Found hipfort: [${HIPFORT_LIB}]" )
  endif()
endmacro()
