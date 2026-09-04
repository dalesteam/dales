# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/home/annelot/dales-dev/build_gmd/_deps/rrtmgp-src"
  "/home/annelot/dales-dev/build_gmd/_deps/rrtmgp-build"
  "/home/annelot/dales-dev/build_gmd/_deps/rrtmgp-subbuild/rrtmgp-populate-prefix"
  "/home/annelot/dales-dev/build_gmd/_deps/rrtmgp-subbuild/rrtmgp-populate-prefix/tmp"
  "/home/annelot/dales-dev/build_gmd/_deps/rrtmgp-subbuild/rrtmgp-populate-prefix/src/rrtmgp-populate-stamp"
  "/home/annelot/dales-dev/build_gmd/_deps/rrtmgp-subbuild/rrtmgp-populate-prefix/src"
  "/home/annelot/dales-dev/build_gmd/_deps/rrtmgp-subbuild/rrtmgp-populate-prefix/src/rrtmgp-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/annelot/dales-dev/build_gmd/_deps/rrtmgp-subbuild/rrtmgp-populate-prefix/src/rrtmgp-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/annelot/dales-dev/build_gmd/_deps/rrtmgp-subbuild/rrtmgp-populate-prefix/src/rrtmgp-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
