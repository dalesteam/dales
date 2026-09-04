#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "rte-rrtmgp::rte" for configuration "Release"
set_property(TARGET rte-rrtmgp::rte APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(rte-rrtmgp::rte PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "Fortran"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/librte.a"
  )

list(APPEND _cmake_import_check_targets rte-rrtmgp::rte )
list(APPEND _cmake_import_check_files_for_rte-rrtmgp::rte "${_IMPORT_PREFIX}/lib/librte.a" )

# Import target "rte-rrtmgp::rrtmgp" for configuration "Release"
set_property(TARGET rte-rrtmgp::rrtmgp APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(rte-rrtmgp::rrtmgp PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "Fortran"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/librrtmgp.a"
  )

list(APPEND _cmake_import_check_targets rte-rrtmgp::rrtmgp )
list(APPEND _cmake_import_check_files_for_rte-rrtmgp::rrtmgp "${_IMPORT_PREFIX}/lib/librrtmgp.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
