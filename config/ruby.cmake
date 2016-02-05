# Thunder
set(CMAKE_Fortran_COMPILER "ifort")
set(Fortran_COMPILER_WRAPPER mpif90)

set(USER_Fortran_FLAGS "-traceback -r8 -ftz -extend_source ")
set(USER_Fortran_FLAGS_RELEASE " -O3 -no-prec-div -xHost")
set(USER_Fortran_FLAGS_DEBUG "-fpe0 -O0 -g -check all -check nopointers -check noarg_temp_created ")

set(NETCDF_INCLUDE_DIR "/usr/local/netcdf/intel/15/mvapich2/2.1/4.3.3.1/include")
set(NETCDF_LIB_1       "/usr/local/netcdf/intel/15/mvapich2/2.1/4.3.3.1/lib/libnetcdff.a")
set(NETCDF_LIB_2       "/usr/local/netcdf/intel/15/mvapich2/2.1/4.3.3.1/lib/libnetcdf.a")
set(HDF5_LIB_1         "/usr/local/hdf5/intel/15/mvapich2/2.1/1.8.15/lib/libhdf5_hl.a")
set(HDF5_LIB_2         "/usr/local/hdf5/intel/15/mvapich2/2.1/1.8.15/lib/libhdf5.a")
set(SZIP_LIB           "")
set(LIBS ${NETCDF_LIB_1} ${NETCDF_LIB_2} ${HDF5_LIB_1} ${HDF5_LIB_2} ${SZIP_LIB} m z curl)
