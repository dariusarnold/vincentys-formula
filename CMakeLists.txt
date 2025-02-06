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

###############################################################################
# C Library and executable
###############################################################################

add_library(vinc SHARED src/vinc.c src/ext/vinc.h)
target_include_directories(vinc PUBLIC src/ext)
target_link_libraries(vinc PUBLIC m)
set_property(TARGET vinc PROPERTY POSITION_INDEPENDENT_CODE ON)

add_library(vinc_cpp SHARED src/vinc.cpp src/ext/vinc.hpp)
target_include_directories(vinc PUBLIC src/ext)
set_property(TARGET vinc PROPERTY POSITION_INDEPENDEN_CODE ON)

add_executable(main src/main.c)
target_link_libraries(main PRIVATE vinc)

###############################################################################
# C++ Library and Python extension module
###############################################################################

pybind11_add_module(vinc_cpp_pybind11 MODULE src/vinc.cpp)
target_compile_definitions(vinc_cpp_pybind11 PRIVATE VINC_PYBIND11=1)

###############################################################################
# Installation
###############################################################################

install(TARGETS main DESTINATION bin)
install(TARGETS vinc vinc_cpp DESTINATION lib)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/src/ext/ DESTINATION include)
