
set(BUILD_TESTING OFF)
set(BUILD_SHARED_LIBS OFF)
set(FS_ENABLE_BACKTRACE_TEST ON)
set(FS_ENABLE_OMP OFF)
set(FS_ENABLE_OPENACC  ENABLE_ACC)
set(FS_ENABLE_MIXED_PRECISION  OFF)
set(FS_ENABLE_SINGLE_PRECISION  FP32_FIELDS)
set(fortran-support_ROOT "external/libfortran-support")


# Where the external project will be installed
set(EXTERNAL_INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/libfortan_install)

# Build the external library in a separate CMake instance
ExternalProject_Add(
    libfortran-support_ep
    GIT_REPOSITORY https://gitlab.dkrz.de/icon-libraries/libfortran-support.git
    GIT_TAG 2.2.0
    CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_DIR}
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DBUILD_SHARED_LIBS=OFF
    TEST_COMMAND ""               # skip tests
)

ecbuild_add_library(TARGET libfortran-support_external TYPE INTERFACE DEPENDS libfortran-support_ep)

target_include_directories(libfortran-support_external INTERFACE
${EXTERNAL_INSTALL_DIR}/include
)
target_link_libraries(libfortran-support_external INTERFACE
${EXTERNAL_INSTALL_DIR}/lib/libfortran-support.a
)
add_library(LibFortranSupport::libfortran-support ALIAS libfortran-support_external)
