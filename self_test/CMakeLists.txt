# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered self_test subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(self_test C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../install")
endif()
include(CTest)

set(SUBDIRS generic)

if(OOMPH_WANT_EXTENDED_SELF_TESTS)
  list(
    APPEND
    SUBDIRS
    navier_stokes
    poisson
    matrix_matrix_multiply
    line_visualiser
    deliberately_broken_code_for_self_test_test
    locate_zeta)

  if(OOMPH_HAS_MPI)
    list(APPEND SUBDIRS mpi)
  endif()
endif()

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving self_test subdirectory")
# ------------------------------------------------------------------------------
