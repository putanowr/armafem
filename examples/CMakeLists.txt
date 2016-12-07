cmake_minimum_required(VERSION 3.4)

include_directories(${PROJECT_SOURCE_DIR}/src)
add_executable (simply_supported_beam simply_supported_beam.cpp)
add_executable (elliptic_1D elliptic_1D.cpp)
target_link_libraries(simply_supported_beam  armafem)
target_link_libraries(elliptic_1D 
                      ${Boost_LIBRARIES}
                      armafem armafemIOVtk armafemIOGnuplot)
