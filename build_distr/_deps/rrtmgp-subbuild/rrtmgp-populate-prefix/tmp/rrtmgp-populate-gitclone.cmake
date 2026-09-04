# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

if(EXISTS "/home/annelot/dales-dev/build_distr/_deps/rrtmgp-subbuild/rrtmgp-populate-prefix/src/rrtmgp-populate-stamp/rrtmgp-populate-gitclone-lastrun.txt" AND EXISTS "/home/annelot/dales-dev/build_distr/_deps/rrtmgp-subbuild/rrtmgp-populate-prefix/src/rrtmgp-populate-stamp/rrtmgp-populate-gitinfo.txt" AND
  "/home/annelot/dales-dev/build_distr/_deps/rrtmgp-subbuild/rrtmgp-populate-prefix/src/rrtmgp-populate-stamp/rrtmgp-populate-gitclone-lastrun.txt" IS_NEWER_THAN "/home/annelot/dales-dev/build_distr/_deps/rrtmgp-subbuild/rrtmgp-populate-prefix/src/rrtmgp-populate-stamp/rrtmgp-populate-gitinfo.txt")
  message(STATUS
    "Avoiding repeated git clone, stamp file is up to date: "
    "'/home/annelot/dales-dev/build_distr/_deps/rrtmgp-subbuild/rrtmgp-populate-prefix/src/rrtmgp-populate-stamp/rrtmgp-populate-gitclone-lastrun.txt'"
  )
  return()
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND} -E rm -rf "/home/annelot/dales-dev/build_distr/_deps/rrtmgp-src"
  RESULT_VARIABLE error_code
)
if(error_code)
  message(FATAL_ERROR "Failed to remove directory: '/home/annelot/dales-dev/build_distr/_deps/rrtmgp-src'")
endif()

# try the clone 3 times in case there is an odd git clone issue
set(error_code 1)
set(number_of_tries 0)
while(error_code AND number_of_tries LESS 3)
  execute_process(
    COMMAND "/usr/bin/git"
            clone --no-checkout --config "advice.detachedHead=false" "https://github.com/earth-system-radiation/rte-rrtmgp.git" "rrtmgp-src"
    WORKING_DIRECTORY "/home/annelot/dales-dev/build_distr/_deps"
    RESULT_VARIABLE error_code
  )
  math(EXPR number_of_tries "${number_of_tries} + 1")
endwhile()
if(number_of_tries GREATER 1)
  message(STATUS "Had to git clone more than once: ${number_of_tries} times.")
endif()
if(error_code)
  message(FATAL_ERROR "Failed to clone repository: 'https://github.com/earth-system-radiation/rte-rrtmgp.git'")
endif()

execute_process(
  COMMAND "/usr/bin/git"
          checkout "77ff83ccf645e5bc404c138ca4e7a6e3abf5d963" --
  WORKING_DIRECTORY "/home/annelot/dales-dev/build_distr/_deps/rrtmgp-src"
  RESULT_VARIABLE error_code
)
if(error_code)
  message(FATAL_ERROR "Failed to checkout tag: '77ff83ccf645e5bc404c138ca4e7a6e3abf5d963'")
endif()

set(init_submodules TRUE)
if(init_submodules)
  execute_process(
    COMMAND "/usr/bin/git" 
            submodule update --recursive --init 
    WORKING_DIRECTORY "/home/annelot/dales-dev/build_distr/_deps/rrtmgp-src"
    RESULT_VARIABLE error_code
  )
endif()
if(error_code)
  message(FATAL_ERROR "Failed to update submodules in: '/home/annelot/dales-dev/build_distr/_deps/rrtmgp-src'")
endif()

# Complete success, update the script-last-run stamp file:
#
execute_process(
  COMMAND ${CMAKE_COMMAND} -E copy "/home/annelot/dales-dev/build_distr/_deps/rrtmgp-subbuild/rrtmgp-populate-prefix/src/rrtmgp-populate-stamp/rrtmgp-populate-gitinfo.txt" "/home/annelot/dales-dev/build_distr/_deps/rrtmgp-subbuild/rrtmgp-populate-prefix/src/rrtmgp-populate-stamp/rrtmgp-populate-gitclone-lastrun.txt"
  RESULT_VARIABLE error_code
)
if(error_code)
  message(FATAL_ERROR "Failed to copy script-last-run stamp file: '/home/annelot/dales-dev/build_distr/_deps/rrtmgp-subbuild/rrtmgp-populate-prefix/src/rrtmgp-populate-stamp/rrtmgp-populate-gitclone-lastrun.txt'")
endif()
