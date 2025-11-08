# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered user_drivers subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(demo_drivers C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../install")
endif()
include(CTest)

set(SUBDIRS joe_cool) # <--- Add your own user driver directory here

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving user_drivers subdirectory")
# ------------------------------------------------------------------------------
