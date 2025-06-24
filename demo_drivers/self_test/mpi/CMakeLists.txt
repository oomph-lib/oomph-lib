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
message(VERBOSE "Entered mpi subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(mpi C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

set(SUBDIRS
    generic_mpi
    line_visualiser
    three_d_mesh_dist
    two_d_mesh_dist
    crdoublematrix_copy_constructor
    distribution_concatenation
    matrix_addition
    matrix_concatenation_without_communication
    matrix_concatenation
    vector_concatenation
    vector_concatenation_without_communication
    vector_split
    vector_split_without_communication
    metis_based_partition_tester)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving mpi subdirectory")
# ------------------------------------------------------------------------------
