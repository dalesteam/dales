# Get version from git, at build time
# Inspiration from CMake Cookbook by Bast and di Remigio,
# https://github.com/dev-cafe/cmake-cookbook/blob/master/chapter-06/recipe-07/cxx-example/CMakeLists.txt
# and
# https://crascit.com/2017/04/18/generated-sources-in-cmake-builds/

set (GIT_VERSION "unknown")
find_package(Git QUIET)
if(GIT_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --abbrev=6 --dirty --always --tags
    OUTPUT_VARIABLE GIT_VERSION
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ) 
endif()

message(STATUS "GIT_VERSION = ${GIT_VERSION} ${TARGET_DIR}/modversion.f90")

configure_file(
  modversion.f90.in
  ${TARGET_DIR}/modversion.f90
  )

message(STATUS "done")
