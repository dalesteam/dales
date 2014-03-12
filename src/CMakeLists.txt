FILE(GLOB sourcefiles "*.?90")
add_executable(${PROJECT_NAME} ${sourcefiles})
include_directories(${INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${LIBS})
