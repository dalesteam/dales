macro( dales_find_cuda )
  # Check if a target architecture has been set. If not, default to 80 (A100 & RTX3090)
  # For H100, use 90
  if( NOT DEFINED CMAKE_CUDA_ARCHITECTURES )
    set( CMAKE_CUDA_ARCHITECTURES 80 )
  endif()
  
  # Look for the CUDA Toolkit, which has cuFFT and NVTX
  set( HAVE_CUDA ON )
  find_package( CUDAToolkit )

  if( NOT TARGET CUDA::cufft AND ENABLE_ACC )
    set( HAVE_CUDA OFF )
  endif()

  if( NOT TARGET CUDA::nvToolsExt AND ENABLE_ACC )
    set( HAVE_CUDA OFF )
  endif() 

  if( HAVE_CUDA )
    if( ENABLE_ACC )
      ecbuild_info( "CUDA architecture: [${CMAKE_CUDA_ARCHITECTURE}]" )
      ecbuild_info( "Found cuFFT:       [${CUDA_cufft_LIBRARY}]" )
      ecbuild_info( "Found NVTX:        [${CUDA_nvToolsExt_LIBRARY}]" )
      add_compile_definitions( USE_NVTX )
    endif()
  endif()
endmacro()

