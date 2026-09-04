#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "fortran-support::fortran-support" for configuration "Release"
set_property(TARGET fortran-support::fortran-support APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(fortran-support::fortran-support PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C;Fortran"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libfortran-support.a"
  )

list(APPEND _cmake_import_check_targets fortran-support::fortran-support )
list(APPEND _cmake_import_check_files_for_fortran-support::fortran-support "${_IMPORT_PREFIX}/lib/libfortran-support.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
