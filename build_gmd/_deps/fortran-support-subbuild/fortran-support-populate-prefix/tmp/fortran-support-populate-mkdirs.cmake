# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/home/annelot/dales-dev/build_gmd/_deps/fortran-support-src"
  "/home/annelot/dales-dev/build_gmd/_deps/fortran-support-build"
  "/home/annelot/dales-dev/build_gmd/_deps/fortran-support-subbuild/fortran-support-populate-prefix"
  "/home/annelot/dales-dev/build_gmd/_deps/fortran-support-subbuild/fortran-support-populate-prefix/tmp"
  "/home/annelot/dales-dev/build_gmd/_deps/fortran-support-subbuild/fortran-support-populate-prefix/src/fortran-support-populate-stamp"
  "/home/annelot/dales-dev/build_gmd/_deps/fortran-support-subbuild/fortran-support-populate-prefix/src"
  "/home/annelot/dales-dev/build_gmd/_deps/fortran-support-subbuild/fortran-support-populate-prefix/src/fortran-support-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/annelot/dales-dev/build_gmd/_deps/fortran-support-subbuild/fortran-support-populate-prefix/src/fortran-support-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/annelot/dales-dev/build_gmd/_deps/fortran-support-subbuild/fortran-support-populate-prefix/src/fortran-support-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
