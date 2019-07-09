set(target_name dales4)
set(target_lib dales)

FILE(GLOB sourcefiles "*.f90")

# stand-alone dales program
add_executable(${target_name} ${sourcefiles})

# dales library, e.g. for use with OMUSE
add_library(${target_lib} ${sourcefiles})

target_link_libraries(${target_name} ${NETCDF_LIBS})
install(TARGETS dales4 dales DESTINATION ${CMAKE_BINARY_DIR})
