cmake_minimum_required(VERSION 3.28)

###############################################################################
# Project setup
###############################################################################

project(vincenty
        LANGUAGES C CXX
)

set(CMAKE_CXX_STANDARD 17 CACHE STRING "C++ version selection")
set(CMAKE_CXX_STANDARD_REQUIRED ON)

###############################################################################
# External dependencies setup
###############################################################################

# This needs python development headers and pybind installed
find_package(Python 3.11 COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

find_package(Catch2 3 REQUIRED)

###############################################################################
# C Library and executable
###############################################################################

add_library(vinc SHARED src/vinc.c src/ext/vinc.h)
target_include_directories(vinc PUBLIC src/ext)
target_link_libraries(vinc PUBLIC m)
set_property(TARGET vinc PROPERTY POSITION_INDEPENDENT_CODE ON)

add_library(vinc_cpp SHARED src/vinc.cpp src/ext/vinc.hpp)
target_include_directories(vinc_cpp PUBLIC src/ext)
set_property(TARGET vinc PROPERTY POSITION_INDEPENDEN_CODE ON)

add_executable(main src/main.c)
target_link_libraries(main PRIVATE vinc)

###############################################################################
# C++ Library and Python extension module
###############################################################################

pybind11_add_module(vinc_cpp_pybind11 MODULE src/vinc.cpp src/ext/vinc.hpp)
target_compile_definitions(vinc_cpp_pybind11 PRIVATE VINC_PYBIND11=1)
target_include_directories(vinc_cpp_pybind11 PUBLIC src/ext)

###############################################################################
# Unit Tests
###############################################################################

add_executable(tests_cpp src/tests/test_vinc.cpp)
target_link_libraries(tests_cpp
        PRIVATE
        Catch2::Catch2WithMain
        vinc_cpp
)

add_executable(tests_c src/tests/test_vinc.cpp)
target_link_libraries(tests_c
        PRIVATE
        Catch2::Catch2WithMain
        vinc
)
target_compile_definitions(tests_c PRIVATE VINC_C)

###############################################################################
# Installation
###############################################################################

install(TARGETS main DESTINATION bin)
install(TARGETS vinc vinc_cpp DESTINATION lib)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/src/ext/ DESTINATION include)
