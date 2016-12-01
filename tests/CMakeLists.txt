cmake_minimum_required(VERSION 3.4)

include(ExternalProject)
find_package(Git REQUIRED)

ExternalProject_Add(
    catch
    PREFIX ${CMAKE_BINARY_DIR}/catch
    GIT_REPOSITORY https://github.com/philsquared/Catch.git
    TIMEOUT 10
    UPDATE_COMMAND ${GIT_EXECUTABLE} pull
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
   )
# Expose required variable (CATCH_INCLUDE_DIR) to parent scope
ExternalProject_Get_Property(catch source_dir)
set(CATCH_INCLUDE_DIR ${source_dir}/include CACHE INTERNAL "Path to include folder for Catch")

add_library(CatchLib INTERFACE)
target_include_directories(CatchLib INTERFACE ${CATCH_INCLUDE_DIR})
add_dependencies(CatchLib catch)

include_directories(${CATCH_INCLUDE_DIR})
add_executable (armafemtest testmain.cpp test_smoke.cpp)

target_link_libraries(armafemtest CatchLib)

add_test (NAME ArmaFemTest COMMAND armafemtest)
