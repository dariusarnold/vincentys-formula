cmake_minimum_required(VERSION 3.28)
project(vincenty)

add_library(vinc SHARED src/vinc.c src/ext/vinc.h)
target_include_directories(vinc PUBLIC src/ext)
target_link_libraries(vinc PUBLIC m)
set_property(TARGET vinc PROPERTY POSITION_INDEPENDENT_CODE ON)

add_executable(main src/main.c)
target_link_libraries(main PRIVATE vinc)
