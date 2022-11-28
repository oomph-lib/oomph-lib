# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered linearised_axisym_navier_stokes subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(linearised_axisym_navier_stokes C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()
include(CTest)

set(SUBDIRS counter_rotating_disks time_periodic_taylor_couette)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving linearised_axisym_navier_stokes subdirectory")
# ------------------------------------------------------------------------------
