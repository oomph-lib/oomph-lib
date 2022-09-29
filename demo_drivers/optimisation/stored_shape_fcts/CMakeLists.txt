# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered stored_shape_fcts subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(stored_shape_fcts C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME two_d_poisson_stored_shape_fcts
  SOURCES two_d_poisson_stored_shape_fcts.cc
  LIBRARIES oomph::poisson
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME optimisation.stored_shape_fcts
  TARGET_DEPENDENCIES two_d_poisson_stored_shape_fcts
  EXTRA_REQUIRES
  LABELS "optimisation;stored_shape_fcts")

message(VERBOSE "Leaving stored_shape_fcts subdirectory")
# ------------------------------------------------------------------------------
