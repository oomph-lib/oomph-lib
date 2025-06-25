# ------------------------------------------------------------------------------
# OOMPH auto-inserted: warn once if user configures inside an existing build.
if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/CMakeCache.txt"
   AND CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  message(
    WARNING
      "You're re-configuring an existing build tree. For a clean configuration, "
      "delete the whole build folder and start with a fresh directory.")
endif()
# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered unsteady_heat subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
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
