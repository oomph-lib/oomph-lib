# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE
        "Entered generalised_newtonian_axisym_navier_stokes subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(generalised_newtonian_axisym_navier_stokes C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()
include(CTest)

set(SUBDIRS axisym_vibrating_shell)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE
        "Leaving generalised_newtonian_axisym_navier_stokes subdirectory")
# ------------------------------------------------------------------------------
