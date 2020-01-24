set(target_name dales4)
set(target_lib dales)

# Add all Fortran90 files from this directory
# This will also add some unwated files, which we will remove again below
FILE(GLOB sourcefiles "*.f90")

# Use git-version.cmake to create modversion.f90, containing a version string from git, e.g. "4.2-34-g62b85a-dirty"
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/modversion.f90
  COMMAND ${CMAKE_COMMAND} -D TARGET_DIR=${CMAKE_CURRENT_BINARY_DIR} -P ${CMAKE_CURRENT_SOURCE_DIR}/git-version.cmake
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
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


install(TARGETS dales4 dales DESTINATION ${CMAKE_BINARY_DIR})


