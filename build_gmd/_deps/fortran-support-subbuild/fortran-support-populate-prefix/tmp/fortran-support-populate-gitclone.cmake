# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

if(EXISTS "/home/annelot/dales-dev/build_gmd/_deps/fortran-support-subbuild/fortran-support-populate-prefix/src/fortran-support-populate-stamp/fortran-support-populate-gitclone-lastrun.txt" AND EXISTS "/home/annelot/dales-dev/build_gmd/_deps/fortran-support-subbuild/fortran-support-populate-prefix/src/fortran-support-populate-stamp/fortran-support-populate-gitinfo.txt" AND
  "/home/annelot/dales-dev/build_gmd/_deps/fortran-support-subbuild/fortran-support-populate-prefix/src/fortran-support-populate-stamp/fortran-support-populate-gitclone-lastrun.txt" IS_NEWER_THAN "/home/annelot/dales-dev/build_gmd/_deps/fortran-support-subbuild/fortran-support-populate-prefix/src/fortran-support-populate-stamp/fortran-support-populate-gitinfo.txt")
  message(STATUS
    "Avoiding repeated git clone, stamp file is up to date: "
    "'/home/annelot/dales-dev/build_gmd/_deps/fortran-support-subbuild/fortran-support-populate-prefix/src/fortran-support-populate-stamp/fortran-support-populate-gitclone-lastrun.txt'"
  )
  return()
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND} -E rm -rf "/home/annelot/dales-dev/build_gmd/_deps/fortran-support-src"
  RESULT_VARIABLE error_code
)
if(error_code)
  message(FATAL_ERROR "Failed to remove directory: '/home/annelot/dales-dev/build_gmd/_deps/fortran-support-src'")
endif()

# try the clone 3 times in case there is an odd git clone issue
set(error_code 1)
set(number_of_tries 0)
while(error_code AND number_of_tries LESS 3)
  execute_process(
    COMMAND "/usr/bin/git"
            clone --no-checkout --config "advice.detachedHead=false" "https://gitlab.dkrz.de/icon-libraries/libfortran-support.git" "fortran-support-src"
    WORKING_DIRECTORY "/home/annelot/dales-dev/build_gmd/_deps"
    RESULT_VARIABLE error_code
  )
  math(EXPR number_of_tries "${number_of_tries} + 1")
endwhile()
if(number_of_tries GREATER 1)
  message(STATUS "Had to git clone more than once: ${number_of_tries} times.")
endif()
if(error_code)
  message(FATAL_ERROR "Failed to clone repository: 'https://gitlab.dkrz.de/icon-libraries/libfortran-support.git'")
endif()

execute_process(
  COMMAND "/usr/bin/git"
          checkout "34246610b17db29f214fb2f95ca1c9087f09b89c" --
  WORKING_DIRECTORY "/home/annelot/dales-dev/build_gmd/_deps/fortran-support-src"
  RESULT_VARIABLE error_code
)
if(error_code)
  message(FATAL_ERROR "Failed to checkout tag: '34246610b17db29f214fb2f95ca1c9087f09b89c'")
endif()

set(init_submodules TRUE)
if(init_submodules)
  execute_process(
    COMMAND "/usr/bin/git" 
            submodule update --recursive --init 
    WORKING_DIRECTORY "/home/annelot/dales-dev/build_gmd/_deps/fortran-support-src"
    RESULT_VARIABLE error_code
  )
endif()
if(error_code)
  message(FATAL_ERROR "Failed to update submodules in: '/home/annelot/dales-dev/build_gmd/_deps/fortran-support-src'")
endif()

# Complete success, update the script-last-run stamp file:
#
execute_process(
  COMMAND ${CMAKE_COMMAND} -E copy "/home/annelot/dales-dev/build_gmd/_deps/fortran-support-subbuild/fortran-support-populate-prefix/src/fortran-support-populate-stamp/fortran-support-populate-gitinfo.txt" "/home/annelot/dales-dev/build_gmd/_deps/fortran-support-subbuild/fortran-support-populate-prefix/src/fortran-support-populate-stamp/fortran-support-populate-gitclone-lastrun.txt"
  RESULT_VARIABLE error_code
)
if(error_code)
  message(FATAL_ERROR "Failed to copy script-last-run stamp file: '/home/annelot/dales-dev/build_gmd/_deps/fortran-support-subbuild/fortran-support-populate-prefix/src/fortran-support-populate-stamp/fortran-support-populate-gitclone-lastrun.txt'")
endif()
