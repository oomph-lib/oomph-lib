# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered generic subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(generic C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()
include(CTest)

set(SUBDIRS vector_matrix_test block_selector_test complex_matrices_test
    eigen_solver_test problem_test)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving generic subdirectory")
# ------------------------------------------------------------------------------
