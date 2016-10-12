cmake_minimum_required(VERSION 3.6)

project(proj_5)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(SOURCE_FILES_1 proj_5.c)
add_executable(proj_5 ${SOURCE_FILES_1})

#project(tute9)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#set(SOURCE_FILES_2 tute9.c)
#add_executable(tute9 ${SOURCE_FILES_2})

find_library(M_LIB m X11)
find_library(CPGPLOT_LIB cpgplot pgplot png z)
target_link_libraries(proj_5 ${M_LIB})
target_link_libraries(proj_5 ${CPGPLOT_LIB})