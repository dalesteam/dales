list(
  APPEND dales_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/hicfft.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/modnvtx.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/modtimer.f90
)

set( dales_srcs "${dales_srcs}" PARENT_SCOPE )
