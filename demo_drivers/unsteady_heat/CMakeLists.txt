# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered unsteady_heat subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(unsteady_heat C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()
include(CTest)

set(SUBDIRS
    space_time_two_d_unsteady_heat
    two_d_unsteady_heat
    two_d_unsteady_heat_2adapt
    two_d_unsteady_heat_ALE
    two_d_unsteady_heat_adapt
    two_d_unsteady_heat_t_adapt)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving unsteady_heat subdirectory")
# ------------------------------------------------------------------------------
