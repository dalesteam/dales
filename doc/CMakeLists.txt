if( ENABLE_DOCS )

  include(FetchContent)

  ecbuild_find_package( Doxygen )

  if( Doxygen_FOUND )

    doxygen_add_docs( doc ${PROJECT_SOURCE_DIR}/src )

  endif()
endif()