include(ExternalProject)

# Tixi library, from ICON. Adds support for XML files
ExternalProject_Add(
  libtixi
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tixi
  CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tixi/configure --prefix=${PROJECT_BINARY_DIR}/tixi
  BUILD_COMMAND make
  INSTALL_COMMAND make install
)
