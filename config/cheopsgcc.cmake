# CHEOPS GCC
set(CMAKE_Fortran_COMPILER "gfortran")
set(Fortran_COMPILER_WRAPPER mpif90)

set(USER_Fortran_FLAGS "-traceback -r8 -ftz -extend_source")
set(USER_Fortran_FLAGS_RELEASE "-O3 -no-prec-div -xHOST -fp-model source")
set(USER_Fortran_FLAGS_DEBUG "-fpe0 -O0 -g -check all -check nopointers -check noarg_temp_created")

set(NETCDF_INCLUDE_DIR "/opt/rrzk/lib/netcdf/4.1.3/include")
set(NETCDF_LIB_1       "/opt/rrzk/lib/netcdf/4.1.3/lib/libnetcdff.a")
set(NETCDF_LIB_2       "/opt/rrzk/lib/netcdf/4.1.3/lib/libnetcdf.a")
set(HDF5_LIB_1         "/opt/rrzk/lib/hdf5/1.8.11/lib/libhdf5_hl.a")
set(HDF5_LIB_2         "/opt/rrzk/lib/hdf5/1.8.11/lib/libhdf5.a")
set(SZIP_LIB           "/opt/rrzk/lib/szip/szip-2.1/lib/libsz.a")
set(LIBS ${NETCDF_LIB_1} ${NETCDF_LIB_2} ${HDF5_LIB_1} ${HDF5_LIB_2} ${SZIP_LIB} m z curl)
