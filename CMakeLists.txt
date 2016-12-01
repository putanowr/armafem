cmake_minimum_required(VERSION 3.4)

project(armafem)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

string(REPLACE ":" ";" LIB_PATHS $ENV{LD_LIBRARY_PATH})
string(REPLACE ":" ";" INC_PATHS $ENV{CPLUS_INCLUDE_PATH})

set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${LIB_PATHS})
set(CMAKE_INCLUDE_PATH ${CMAKE_LIBRARY_PATH} ${INC_PATHS})

enable_testing(true)  # Enables unit-testing.

find_package(Armadillo REQUIRED)

add_subdirectory(src)
add_subdirectory(tests)

#include_directories(${ARMADILLO_INCLUDE_DIRS})
#add_executable(armasample01 armasample01.cpp)
#target_link_libraries(armasample01 ${ARMADILLO_LIBRARIES})
