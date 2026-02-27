# Check if the selected NetCDF library has parallel support
macro(dales_query_parallel_netcdf)

  find_program( NC_CONFIG_EXECUTABLE
    NAMES nc-config
    DOC "nc-config executable"
  )
  
  if( NC_CONFIG_EXECUTABLE )
    execute_process( COMMAND ${NC_CONFIG_EXECUTABLE} --has-parallel
      RESULT_VARIABLE nc_config_result
      OUTPUT_VARIABLE nc_config_has_parallel
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if( nc_config_result EQUAL 0)
      string( COMPARE EQUAL "${nc_config_has_parallel}" "yes" has_parallel )
      if( has_parallel )
        ecbuild_info( "NetCDF parallel functionality enabled" )
        add_compile_definitions( NC_HAS_PARALLEL )
      else()
        ecbuild_warn( "NetCDF library does not have parallel functionality" )
      endif()
    endif()
  else()
    ecbuild_warn( "Could not locate nc-config executable - disabling NetCDF parallel functionality" )
  endif()

endmacro()
