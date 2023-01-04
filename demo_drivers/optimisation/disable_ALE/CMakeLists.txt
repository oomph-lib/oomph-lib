# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered disable_ALE subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(disable_ALE C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

set(SUBDIRS unsteady_heat navier_stokes)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving disable_ALE subdirectory")
# ------------------------------------------------------------------------------
