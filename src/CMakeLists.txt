cmake_minimum_required(VERSION 3.4)

include_directories(${ARMADILLO_INCLUDE_DIRS})
add_library(armafem 
  assem.cpp
  beam2e.cpp
  solveq.cpp
)
target_link_libraries(armafem ${ARMADILLO_LIBRARIES})
