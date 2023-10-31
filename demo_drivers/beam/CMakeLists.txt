# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered beam subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(beam C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()
include(CTest)

set(SUBDIRS tensioned_string steady_ring steady_third_ring unsteady_ring)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving beam subdirectory")
# ------------------------------------------------------------------------------
