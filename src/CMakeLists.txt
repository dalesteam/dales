set(target_name dales4.4)
set(target_lib dales)

# Add all Fortran90 files from this directory
# This will also add some unwated files, which we will remove again below
FILE(GLOB sourcefiles "*.f90" 
	"RRTMG/RRTMG_LW/modules/*.f90"
	"RRTMG/RRTMG_LW/src/*.f90"
	"RRTMG/RRTMG_SW/modules/*.f90"
	"RRTMG/RRTMG_SW/src/*.f90"
	"RTE-RRTMGP/rte/*.F90"
	"RTE-RRTMGP/rte/kernels/*.F90"
	"RTE-RRTMGP/rrtmgp/*.F90"
	"RTE-RRTMGP/rrtmgp/kernels/*.F90"
	"RTE-RRTMGP/examples/mo_simple_netcdf.F90"
	"RTE-RRTMGP/examples/mo_load_coefficients.F90")

# TODO: make more concise
list(REMOVE_ITEM sourcefiles ${CMAKE_CURRENT_SOURCE_DIR}/test_transposes.f90)
list(REMOVE_ITEM sourcefiles ${CMAKE_CURRENT_SOURCE_DIR}/RRTMG/RRTMG_LW/src/mcica_random_numbers.f90)
list(REMOVE_ITEM sourcefiles ${CMAKE_CURRENT_SOURCE_DIR}/RRTMG/RRTMG_LW/src/mcica_subcol_gen_lw.f90)
list(REMOVE_ITEM sourcefiles ${CMAKE_CURRENT_SOURCE_DIR}/RRTMG/RRTMG_LW/src/mcica_subcol_gen_lw.1col.f90)
list(REMOVE_ITEM sourcefiles ${CMAKE_CURRENT_SOURCE_DIR}/RRTMG/RRTMG_LW/src/rrtmg_lw.1col.f90)
list(REMOVE_ITEM sourcefiles ${CMAKE_CURRENT_SOURCE_DIR}/RRTMG/RRTMG_LW/src/rrtmg_lw_cldprmc.f90)
list(REMOVE_ITEM sourcefiles ${CMAKE_CURRENT_SOURCE_DIR}/RRTMG/RRTMG_LW/src/rrtmg_lw_k_g.f90)
list(REMOVE_ITEM sourcefiles ${CMAKE_CURRENT_SOURCE_DIR}/RRTMG/RRTMG_LW/src/rrtmg_lw_rad.f90)
list(REMOVE_ITEM sourcefiles ${CMAKE_CURRENT_SOURCE_DIR}/RRTMG/RRTMG_LW/src/rrtmg_lw_rrtmc.f90)
list(REMOVE_ITEM sourcefiles ${CMAKE_CURRENT_SOURCE_DIR}/RRTMG/RRTMG_SW/src/mcica_random_numbers.f90)
list(REMOVE_ITEM sourcefiles ${CMAKE_CURRENT_SOURCE_DIR}/RRTMG/RRTMG_SW/src/mcica_subcol_gen_sw.f90)
list(REMOVE_ITEM sourcefiles ${CMAKE_CURRENT_SOURCE_DIR}/RRTMG/RRTMG_SW/src/mcica_subcol_gen_sw.1col.f90)
list(REMOVE_ITEM sourcefiles ${CMAKE_CURRENT_SOURCE_DIR}/RRTMG/RRTMG_SW/src/rrtmg_sw.1col.f90)
list(REMOVE_ITEM sourcefiles ${CMAKE_CURRENT_SOURCE_DIR}/RRTMG/RRTMG_SW/src/rrtmg_sw_cldprmc.f90)
list(REMOVE_ITEM sourcefiles ${CMAKE_CURRENT_SOURCE_DIR}/RRTMG/RRTMG_SW/src/rrtmg_sw_k_g.f90)
list(REMOVE_ITEM sourcefiles ${CMAKE_CURRENT_SOURCE_DIR}/RRTMG/RRTMG_SW/src/rrtmg_sw_rad.f90)
list(REMOVE_ITEM sourcefiles ${CMAKE_CURRENT_SOURCE_DIR}/RRTMG/RRTMG_SW/src/rrtmg_sw_spcvmc.f90)

# Still a hack and only works for GCC
set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/RRTMG/RRTMG_LW/src/rrtmg_lw_init.f90 PROPERTIES COMPILE_FLAGS -Wno-error=implicit-interface)
set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/RRTMG/RRTMG_SW/src/rrtmg_sw_init.f90 PROPERTIES COMPILE_FLAGS -Wno-error=implicit-interface)

# Use git-version.cmake to create modversion.f90, containing a version string from git, e.g. "4.2-34-g62b85a-dirty"
add_custom_target(tag_git_version ALL
  COMMAND ${CMAKE_COMMAND} -D TARGET_DIR=${CMAKE_CURRENT_BINARY_DIR} -P ${CMAKE_CURRENT_SOURCE_DIR}/git-version.cmake
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/modversion.f90
  )

# Stand-alone DALES program
add_executable(${target_name} ${sourcefiles} ${CMAKE_CURRENT_BINARY_DIR}/modversion.f90)
target_link_libraries(${target_name} ${NETCDF_LIBS} ${OPTIONAL_LIBS})
add_dependencies(${target_name} tag_git_version)

# DALES library, e.g. for use with OMUSE
add_library(${target_lib} ${sourcefiles} ${CMAKE_CURRENT_BINARY_DIR}/modversion.f90)
add_dependencies(${target_lib} tag_git_version)

# set separate module directories for the program and library targets, to avoid problems during parallel builds (make -j 8)
set_target_properties(${target_name} PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/program_modules")
set_target_properties(${target_lib} PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/library_modules")


install(TARGETS ${target_name} ${target_lib} DESTINATION ${CMAKE_BINARY_DIR})
