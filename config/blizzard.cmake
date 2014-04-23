# ARCH Linux
set(CMAKE_Fortran_COMPILER "xlf")
set(Fortran_COMPILER_WRAPPER mpxlf90)

set(USER_Fortran_FLAGS "-qfree=F90 -qrealsize=8  -qwarn64 -qnosave -qinitauto=FFF00000 -qflttrap=en:ov:zero:inv:imp -qflag=w:e")
set(USER_Fortran_FLAGS_RELEASE "-O4 -qnoipa -qstrict=none:exceptions -qinitauto=ff -qsigtrap")
set(USER_Fortran_FLAGS_DEBUG "-O0 -qfullpath -C -g -qflttrp=enable:inexact:invalid:nanq:overflow:zerodivide -qsigtrap -qinitauto")

set(NETCDF_INCLUDE_DIR "/sw/aix61/netcdf-4.1.2-hdf5-threadsafe/include")
set(NETCDF_LIB_1       "/sw/aix61/netcdf-4.1.2-hdf5-threadsafe/lib/libnetcdff.a")
set(NETCDF_LIB_2       "/sw/aix61/netcdf-4.1.2-hdf5-threadsafe/lib/libnetcdf.a")
set(HDF5_LIB_1         "/sw/aix61/hdf5-1.8.6-threadsafe/lib/libhdf5_hl.a")
set(HDF5_LIB_2         "/sw/aix61/hdf5-1.8.6-threadsafe/lib/libhdf5.a")
set(SZIP_LIB           "")
set(LIBS ${NETCDF_LIB_1} ${NETCDF_LIB_2} ${HDF5_LIB_1} ${HDF5_LIB_2} ${SZIP_LIB} m z)
