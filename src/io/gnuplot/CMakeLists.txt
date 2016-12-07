cmake_minimum_required(VERSION 3.4)

include_directories(
  ${PROJECT_SOURCE_DIR}/src
  ${ARMADILLO_INCLUDE_DIRS}
  ${BOOST_INCLUDE_DIRS}
)

add_library(armafemIOGnuplot
  gnuplotExporter.cpp
)

target_link_libraries(armafemIOGnuplot
  armafemIOUtils
  ${ARMADILLO_LIBRARIES}
  ${BOOST_LIBRARIES}
)
