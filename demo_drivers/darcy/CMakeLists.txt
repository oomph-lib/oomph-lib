# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered darcy subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(darcy C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()
include(CTest)

set(SUBDIRS unstructured_two_d_circle)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving darcy subdirectory")
# ------------------------------------------------------------------------------
