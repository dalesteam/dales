# Check the NetCDF library for Zstandard compression support.
macro( dales_check_zstandard )

  ecbuild_info( "Testing NetCDF for Zstandard support" )

  set( TMPDIR "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp" )

  file( WRITE "${TMPDIR}/test_zstd.f90"
        "
        program test_zstd
        use netcdf
        implicit none
        integer :: ncid, dimid, varid, istat
        istat = nf90_create('${TMPDIR}/test.nc', NF90_CLOBBER, ncid)
        istat = nf90_def_dim(ncid, 'x', 1, dimid)
        istat = nf90_def_var(ncid, 'y', NF90_INT, [dimid], varid)
        istat = nf90_def_var_zstandard(ncid, varid, 4)
        stop istat
        end program
        "
  )

  ecbuild_try_run( EXITCODE COMPILED
                   ${TMPDIR} "${TMPDIR}/test_zstd.f90"
                   CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${NetCDF_Fortran_INCLUDE_DIR}"
                   LINK_LIBRARIES ${NetCDF_Fortran_LIBRARY}
  )

  if( ${COMPILED} )
    if( ${EXITCODE} EQUAL 0 )
      ecbuild_info( "Testing NetCDF for Zstandard support - Success" )
      add_compile_definitions( NC_ZSTANDARD )
    else()
      ecbuild_info( "Testing NetCDF for Zstandard support - Failed" )
      ecbuild_info( "   - test_zstd.f90 did not run succesfully. NetCDF-C is likely not compiled with Zstandard support." )
    endif()
  else()
    ecbuild_info( "Testing NetCDF for Zstandard support - Failed" )
    ecbuild_info( "   - test_zstd.f90 did not compile succesfully. You are likely using an old version of NetCDF-Fortran." )
  endif()

endmacro()
