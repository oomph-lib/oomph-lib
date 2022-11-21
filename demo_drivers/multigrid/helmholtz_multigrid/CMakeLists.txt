# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered helmholtz_multigrid subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(helmholtz_multigrid C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME structured_cubic_point_source
  SOURCES structured_cubic_point_source.cc
  LIBRARIES oomph::pml_helmholtz
  CXX_STANDARD 17)

oomph_add_executable(
  NAME unstructured_two_d_helmholtz
  SOURCES unstructured_two_d_helmholtz.cc
  LIBRARIES oomph::pml_helmholtz
  CXX_STANDARD 17)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME multigrid.helmholtz_multigrid
  TARGET_DEPENDENCIES structured_cubic_point_source unstructured_two_d_helmholtz
  EXTRA_REQUIRES
  LABELS multigrid helmholtz_multigrid)

message(VERBOSE "Leaving helmholtz_multigrid subdirectory")
# ------------------------------------------------------------------------------
