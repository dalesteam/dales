include( FetchContent )

set( BUILD_C_HEADERS OFF )
set( RTE_ENABLE_SP ${ENABLE_FP32_RAD} )
set( KERNEL_MODE "default" CACHE STRING "Set kernel mode for RTE" )

if( ${ENABLE_ACC} )
    set( KERNEL_MODE "accel" )
endif()

ecbuild_info( "Fetching RTE-RRTMGP" )

FetchContent_Declare(
    rrtmgp
    GIT_REPOSITORY "https://github.com/earth-system-radiation/rte-rrtmgp.git"
    GIT_TAG 77ff83ccf645e5bc404c138ca4e7a6e3abf5d963 # v1.9.2
)

FetchContent_MakeAvailable( rrtmgp )

# For some reason, RRTMGP doesn't pick up the OpenACC compiler flag,
# so we set it here manually
target_compile_options( rrtmgp PUBLIC ${OpenACC_Fortran_FLAGS} )
target_compile_options( rrtmgpkernels PUBLIC ${OpenACC_Fortran_FLAGS} )
target_compile_options( rte PUBLIC ${OpenACC_Fortran_FLAGS} )
target_compile_options( rtekernels PUBLIC ${OpenACC_Fortran_FLAGS} )
