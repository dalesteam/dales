include(FetchContent)

# Set CMake options for fortran-support
set( BUILD_TESTING OFF )
set( BUILD_SHARED_LIBS OFF )
set( FS_ENABLE_OPENACC ${ENABLE_ACC} )

ecbuild_info( "Fetching fortran-support" )

FetchContent_Declare(
    fortran-support
    GIT_REPOSITORY https://gitlab.dkrz.de/icon-libraries/libfortran-support.git
    GIT_TAG 2.2.0
    OVERRIDE_FIND_PACKAGE
)