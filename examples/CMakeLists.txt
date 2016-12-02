cmake_minimum_required(VERSION 3.4)

include_directories(${PROJECT_SOURCE_DIR}/src)
add_executable (simply_supported_beam
                simply_supported_beam.cpp
               )

target_link_libraries(simply_supported_beam  armafem)
