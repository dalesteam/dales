# Default GCC
set(CMAKE_Fortran_COMPILER "gfortran")
set(Fortran_COMPILER_WRAPPER mpif90)

set(USER_Fortran_FLAGS "-fbacktrace -finit-real=nan -fdefault-real-8  -fno-f2c -ffree-line-length-none")
set(USER_Fortran_FLAGS_RELEASE "-funroll-all-loops -O3 -march=native -mtune=native")
set(USER_Fortran_FLAGS_DEBUG "-W -Wall -Wuninitialized -fcheck=all -fbacktrace -O0 -g -ffpe-trap=invalid,zero,overflow")

set(NETCDF_INCLUDE_DIR "/sw/squeeze-x64/netcdf-latest-static-gcc47/include")
set(NETCDF_LIB_1       "/sw/squeeze-x64/netcdf-latest-static-gcc47/lib/libnetcdff.a")
set(NETCDF_LIB_2       "/sw/squeeze-x64/netcdf-latest-static-gcc47/lib/libnetcdf.a")
set(HDF5_LIB_1         "/sw/squeeze-x64/hdf5-1.8.7-static/lib/libhdf5_hl.a")
set(HDF5_LIB_2         "/sw/squeeze-x64/hdf5-1.8.7-static/lib/libhdf5.a")
set(SZIP_LIB           "/sw/squeeze-x64/szip-2.1-static/lib/libsz.a")
set(LIBS ${NETCDF_LIB_1} ${NETCDF_LIB_2} ${HDF5_LIB_1} ${HDF5_LIB_2} ${SZIP_LIB} m z curl)
