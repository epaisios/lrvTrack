#SET(CMAKE_VERBOSE_MAKEFILE on)
set( CMAKE_EXPORT_COMPILE_COMMANDS 1 )
cmake_minimum_required (VERSION 2.6)
project(lrvTrack)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeScripts")
set(CMAKE_CXX_FLAGS "-Wc++11-extensions -march=native")
set(CMAKE_CXX_FLAGS_DEBUG "-ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
if(UNIX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
endif()


# The version number.
set (LRVTRACK_VERSION_MAJOR 0)
set (LRVTRACK_VERSION_MINOR 9)


add_subdirectory(src lrvTrack)
