list(
  APPEND dales_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/modlsm.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/modlsmcrosssection.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/modlsmdata.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/modlsmstat.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/modslurb.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/modslurb_energybalance.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/modslurb_radiationmodel.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/modslurbdata.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/modslurbhelpers.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/modslurbresistance_stability.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/modsurface.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/modsurfdata.f90
)

set( dales_srcs "${dales_srcs}" PARENT_SCOPE )