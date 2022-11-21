# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered two_d_poisson_flux_bc subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(two_d_poisson_flux_bc C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME two_d_poisson_flux_bc
  SOURCES two_d_poisson_flux_bc.cc
  LIBRARIES oomph::poisson
  CXX_STANDARD 17)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME poisson.two_d_poisson_flux_bc
  TARGET_DEPENDENCIES two_d_poisson_flux_bc
  EXTRA_REQUIRES
  LABELS poisson two_d_poisson_flux_bc)

message(VERBOSE "Leaving two_d_poisson_flux_bc subdirectory")
# ------------------------------------------------------------------------------
