cmake_minimum_required(VERSION 3.4)

include_directories(
  ${PROJECT_SOURCE_DIR}/src
  ${BOOST_INCLUDE_DIRS}
)

add_library(armafemIOUtils
  io_utils.cpp 
)
target_link_libraries(armafemIOUtils ${BOOST_LIBRARIES})

add_subdirectory(vtk)
add_subdirectory(gnuplot)
