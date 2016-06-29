# Thunder
set(CMAKE_Fortran_COMPILER "ifort")
set(Fortran_COMPILER_WRAPPER mpfort)

set(USER_Fortran_FLAGS "-traceback -r8 -ftz -extend_source")
set(USER_Fortran_FLAGS_RELEASE "-O3 -no-prec-div -fp-model source -xAVX")
set(USER_Fortran_FLAGS_DEBUG "-fpe0 -O0 -g -check all -check nopointers -check noarg_temp_created")

set(NETCDF_INCLUDE_DIR "/lrz/sys/libraries/netcdf/4.2.1.1/include")
set(NETCDF_LIB_1       "/lrz/sys/libraries/netcdf/4.2.1.1/lib/libnetcdff.a")
set(NETCDF_LIB_2       "/lrz/sys/libraries/netcdf/4.2.1.1/lib/libnetcdf.a")
set(HDF5_LIB_1         "/lrz/sys/libraries/netcdf/hdf5_1.8.9/lib/libhdf5_hl.a")
set(HDF5_LIB_2         "/lrz/sys/libraries/netcdf/hdf5_1.8.9/lib/libhdf5.a")
set(SZIP_LIB           "/lrz/sys/libraries/hdf5/szip_2.1_u1/lib/libsz.a")
set(LIBS ${NETCDF_LIB_1} ${NETCDF_LIB_2} ${HDF5_LIB_1} ${HDF5_LIB_2} ${SZIP_LIB} m z curl)
