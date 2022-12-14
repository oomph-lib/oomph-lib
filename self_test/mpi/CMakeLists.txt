# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered mpi subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(mpi C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()
include(CTest)

set(SUBDIRS two_d_mesh_dist three_d_mesh_dist line_visualiser)

if(OOMPH_WANT_MPI_VARIABLENP)
  list(
    APPEND
    SUBDIRS
    distribution_concatenation
    vector_split_without_communication
    vector_concatenation_without_communication
    matrix_concatenation_without_communication
    vector_split
    vector_concatenation
    matrix_concatenation
    crdoublematrix_copy_constructor
    matrix_addition
    generic_mpi)
endif()

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving mpi subdirectory")
# ------------------------------------------------------------------------------
