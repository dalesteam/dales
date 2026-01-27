macro( dales_find_cuda )
  # Check if a target architecture has been set. If not, default to 80 (A100 & RTX3090)
  # For H100, use 90
  if( NOT DEFINED CMAKE_CUDA_ARCHITECTURES )
    set( CMAKE_CUDA_ARCHITECTURES 80 )
  endif()
  
  if( ${CMAKE_Fortran_COMPILER_ID} MATCHES PGI|NVHPC )
    set( HAVE_CUDA ON )
    ecbuild_add_fortran_flags( "-cudalib=cufft,nvtx3" )
  else()
    ecbuild_info( "Could not find CUDA" )
  endif()

endmacro()

