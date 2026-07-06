include(FetchContent)

# Set CMake options for fortran-support
set( BUILD_TESTING OFF )
set( BUILD_SHARED_LIBS OFF )
set( FS_ENABLE_OPENACC ${ENABLE_ACC} )

ecbuild_info( "Fetching fortran-support" )

FetchContent_Declare(
    fortran-support
    GIT_REPOSITORY "https://gitlab.dkrz.de/icon-libraries/libfortran-support.git"
    GIT_TAG 34246610b17db29f214fb2f95ca1c9087f09b89c # git tag 2.2.1
    PATCH_COMMAND git apply --reject --whitespace=fix ${PROJECT_SOURCE_DIR}/patches/libfortran-support.patch
    UPDATE_DISCONNECTED TRUE
)

FetchContent_MakeAvailable(fortran-support)
