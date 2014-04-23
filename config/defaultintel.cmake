# Thunder
set(CMAKE_Fortran_COMPILER "ifort")
set(Fortran_COMPILER_WRAPPER mpif90)

set(USER_Fortran_FLAGS "-traceback -r8 -ftz -extend_source")
set(USER_Fortran_FLAGS_RELEASE "-O3 -no-prec-div -xHOST -fp-model source")
set(USER_Fortran_FLAGS_DEBUG "-fpe0 -O0 -g -check all -check nopointers -check noarg_temp_created")

set(NETCDF_INCLUDE_DIR "/sw/squeeze-x64/netcdf_fortran-latest-static-intel13/include")
set(NETCDF_LIB_1       "/sw/squeeze-x64/netcdf_fortran-latest-static-intel13/lib/libnetcdff.a")
set(NETCDF_LIB_2       "/sw/squeeze-x64/netcdf-4.2-static/lib/libnetcdf.a")
set(HDF5_LIB_1         "/sw/squeeze-x64/hdf5-1.8.7-static/lib/libhdf5_hl.a")
set(HDF5_LIB_2         "/sw/squeeze-x64/hdf5-1.8.7-static/lib/libhdf5.a")
set(SZIP_LIB           "")
set(LIBS ${NETCDF_LIB_1} ${NETCDF_LIB_2} ${HDF5_LIB_1} ${HDF5_LIB_2} ${SZIP_LIB} m z curl)
