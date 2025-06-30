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
message(VERBOSE "Entered self_test subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(self_test C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()

include(CTest)

set(SUBDIRS
    generic
    navier_stokes
    poisson
    matrix_matrix_multiply
    line_visualiser
    deliberately_broken_code_for_self_test_test
    locate_zeta
    faceted_surface_vtu_output)

if(OOMPH_HAS_MPI)
  list(APPEND SUBDIRS mpi)
endif()

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving self_test subdirectory")
# ------------------------------------------------------------------------------
