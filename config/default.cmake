# ARCH Linux
set(CMAKE_Fortran_COMPILER "gfortran")
set(Fortran_COMPILER_WRAPPER mpif90)

set(USER_Fortran_FLAGS "-fbacktrace -finit-real=nan -fdefault-real-8  -fno-f2c -ffree-line-length-none")
set(USER_Fortran_FLAGS_RELEASE "-funroll-all-loops -O3 -march=native -mtune=native")
set(USER_Fortran_FLAGS_DEBUG "-W -Wall -Wuninitialized -fcheck=all -fbacktrace -O0 -g -ffpe-trap=invalid,zero,overflow")

set(NETCDF_INCLUDE_DIR "/usr/include")
set(NETCDF_LIB_1       "/usr/lib64/libnetcdff.so")
set(NETCDF_LIB_2       "/usr/lib64/libnetcdf.so")
set(HDF5_LIB_1         "/usr/lib64/libhdf5_hl.so")
set(HDF5_LIB_2         "/usr/lib64/libhdf5.so")
# set(HDF4_LIB_1         "/usr/lib64/hdf/libdf.a")
# set(HDF4_LIB_2         "/usr/lib64/hdf/libmfhdf.a")
set(HDF4_LIB_1         "")
set(HDF4_LIB_2         "")
set(SZIP_LIB           "")
set(LIBS ${NETCDF_LIB_1} ${NETCDF_LIB_2} ${HDF5_LIB_1} ${HDF5_LIB_2} ${HDF4_LIB_2} ${HDF4_LIB_1} ${SZIP_LIB} dl m z curl)
