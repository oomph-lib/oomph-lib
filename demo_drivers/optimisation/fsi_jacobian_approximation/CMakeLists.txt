# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered fsi_jacobian_approximation subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(fsi_jacobian_approximation C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME fsi_jacobian_approximation
  SOURCES fsi_jacobian_approximation.cc
  LIBRARIES oomph::beam oomph::navier_stokes
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME optimisation.fsi_jacobian_approximation
  TARGET_DEPENDENCIES fsi_jacobian_approximation
  EXTRA_REQUIRES
  LABELS optimisation fsi_jacobian_approximation)

message(VERBOSE "Leaving fsi_jacobian_approximation subdirectory")
# ------------------------------------------------------------------------------
