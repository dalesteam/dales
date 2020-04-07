set(target_name dales4)
set(target_lib dales)

# Add all Fortran90 files from this directory
# This will also add some unwated files, which we will remove again below
FILE(GLOB sourcefiles "*.f90")

# Use git-version.cmake to create modversion.f90, containing a version string from git, e.g. "4.2-34-g62b85a-dirty"
add_custom_target(tag_git_version ALL
  COMMAND ${CMAKE_COMMAND} -D TARGET_DIR=${CMAKE_CURRENT_BINARY_DIR} -P ${CMAKE_CURRENT_SOURCE_DIR}/git-version.cmake
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/modversion.f90
  )

if(USE_HYPRE)
  # Build the iterative solver
  list(REMOVE_ITEM sourcefiles ${CMAKE_CURRENT_SOURCE_DIR}/modnohypre.f90)

  # Stand-alone DALES program
  add_executable(${target_name} ${sourcefiles} ${CMAKE_CURRENT_BINARY_DIR}/modversion.f90) 
  target_link_libraries(${target_name} ${NETCDF_LIBS} ${HYPRE_LIB})

else(USE_HYPRE)
  # Do not build the iterative solver, but use (empty) fallback functions
  list(REMOVE_ITEM sourcefiles ${CMAKE_CURRENT_SOURCE_DIR}/modhypre.f90)

  # Stand-alone DALES program
  add_executable(${target_name} ${sourcefiles} ${CMAKE_CURRENT_BINARY_DIR}/modversion.f90) 
  target_link_libraries(${target_name} ${NETCDF_LIBS})

endif(USE_HYPRE)

# DALES library, e.g. for use with OMUSE
add_library(${target_lib} ${sourcefiles} ${CMAKE_CURRENT_BINARY_DIR}/modversion.f90)


# set separate module directories for the program and library targets, to avoid problems during parallel builds (make -j 8)
set_target_properties(${target_name} PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/program_modules")
set_target_properties(${target_lib} PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/library_modules")

target_link_libraries(${target_name} ${NETCDF_LIBS})

install(TARGETS dales4 dales DESTINATION ${CMAKE_BINARY_DIR})


