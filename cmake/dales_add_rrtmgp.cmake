include( FetchContent )

set( BUILD_C_HEADERS OFF )
set( RTE_ENABLE_SP ${ENABLE_FP32_RAD} )
set( KERNEL_MODE "default" CACHE STRING "Set kernel mode for RTE" )

if( ${ENABLE_ACC} OR ${ENABLE_OMP_GPU} )
    set( KERNEL_MODE "accel" )
endif()

ecbuild_info( "Fetching RTE-RRTMGP" )

FetchContent_Declare(
    rrtmgp
    GIT_REPOSITORY "https://github.com/dindon-sournois/rte-rrtmgp.git"
    GIT_TAG 6b5d8f46aae1058d4bf15de87253b4e3eaf53ac5 # v1.9.2 with OpenMP and AMDFlang fix
)

FetchContent_MakeAvailable( rrtmgp )

# For some reason, RRTMGP doesn't pick up the OpenACC compiler flag,
# so we set it here manually
if( ${ENABLE_ACC} )
  target_compile_options( rrtmgp PUBLIC ${OpenACC_Fortran_FLAGS} )
  target_compile_options( rrtmgpkernels PUBLIC ${OpenACC_Fortran_FLAGS} )
  target_compile_options( rte PUBLIC ${OpenACC_Fortran_FLAGS} )
  target_compile_options( rtekernels PUBLIC ${OpenACC_Fortran_FLAGS} )
endif()

# remove SIMD from pure routine
if( ${CMAKE_Fortran_COMPILER_ID} MATCHES Cray )
  set(FILE_TO_PATCH ${rrtmgp_SOURCE_DIR}/rte-kernels/mo_rte_solver_kernels.F90)
  file(READ ${FILE_TO_PATCH} FILE_CONTENTS)
  string(REPLACE "!$OMP SIMD" "" FILE_CONTENTS "${FILE_CONTENTS}")
  file(WRITE ${FILE_TO_PATCH} "${FILE_CONTENTS}")
endif()
