
####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was config.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

####################################################################################

include(${CMAKE_CURRENT_LIST_DIR}/rte-rrtmgp-targets.cmake)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(rte-rrtmgp REQUIRED_VARS rte-rrtmgp_DIR)
