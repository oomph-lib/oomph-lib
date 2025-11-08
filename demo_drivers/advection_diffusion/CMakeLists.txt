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
message(VERBOSE "Entered advection_diffusion subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(advection_diffusion C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()

include(CTest)

set(SUBDIRS two_d_adv_diff_flux_bc two_d_adv_diff_adapt two_d_adv_diff_SUPG
    two_d_variable_diff)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving advection_diffusion subdirectory")
# ------------------------------------------------------------------------------
