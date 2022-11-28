# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered spherical_navier_stokes subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(spherical_navier_stokes C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()
include(CTest)

set(SUBDIRS spin_up spherical_couette steady_rot)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving spherical_navier_stokes subdirectory")
# ------------------------------------------------------------------------------
