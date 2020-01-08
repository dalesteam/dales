set(target_name dales4)
set(target_lib dales)

FILE(GLOB sourcefiles "*.f90")

# Stand-alone DALES program
add_executable(${target_name} ${sourcefiles}) 

# DALES library, e.g. for use with OMUSE
add_library(${target_lib} ${sourcefiles})

# Use git-version.cmake to create modversion.f90, containing a version string from git, e.g. "4.2-34-g62b85a-dirty"
add_custom_command(
  OUTPUT
     ${CMAKE_CURRENT_BINARY_DIR}/modversion.f90
    ALL
  COMMAND
    ${CMAKE_COMMAND} -D TARGET_DIR=${CMAKE_CURRENT_BINARY_DIR} -P ${CMAKE_CURRENT_SOURCE_DIR}/git-version.cmake
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_SOURCE_DIR}
  )

# Rebuild modversion.f90 every time
add_custom_target(
  get_git_version
  ALL
  DEPENDS
    ${CMAKE_CURRENT_BINARY_DIR}/modversion.f90
  )

# modversion.f90 has to be generated before we start building dales
add_dependencies(${target_name} get_git_version)
add_dependencies(${target_lib} get_git_version)


target_link_libraries(${target_name} ${NETCDF_LIBS})
install(TARGETS dales4 dales DESTINATION ${CMAKE_BINARY_DIR})
