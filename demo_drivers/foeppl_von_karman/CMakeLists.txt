# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered foeppl_von_karman subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(foeppl_von_karman C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()
include(CTest)

set(SUBDIRS circular_disk)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving foeppl_von_karman subdirectory")
# ------------------------------------------------------------------------------
