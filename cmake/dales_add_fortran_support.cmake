include(FetchContent)

# Set CMake options for fortran-support
set( BUILD_TESTING OFF )
set( BUILD_SHARED_LIBS OFF )
set( FS_ENABLE_OPENACC ${ENABLE_ACC} )

ecbuild_info( "Fetching fortran-support" )

FetchContent_Declare(
    fortran-support
    GIT_REPOSITORY "https://gitlab.dkrz.de/icon-libraries/libfortran-support.git"
    GIT_TAG d0d20147dfe96b41b2b8d2a9892e19aaeb3d6bbf # git tag 2.2.0
)

FetchContent_MakeAvailable(fortran-support)
