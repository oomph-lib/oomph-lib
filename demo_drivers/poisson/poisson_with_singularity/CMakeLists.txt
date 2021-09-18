# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered poisson_with_singularity subdirectory")

cmake_minimum_required(VERSION 3.16..3.20 FATAL_ERROR)
project(poisson_with_singularity C CXX Fortran)
if(NOT BUILDING_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME two_d_poisson
  SOURCES two_d_poisson.cc poisson_elements_with_singularity.h
  LIBRARIES oomph::poisson
  CXX_STANDARD 14)

oomph_add_executable(
  NAME two_d_poisson_without_singularity_subtracted_off
  SOURCES two_d_poisson.cc poisson_elements_with_singularity.h
  LIBRARIES oomph::poisson
  CXX_DEFINITIONS DONT_USE_SINGULARITY
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME poisson.poisson_with_singularity
  TARGET_DEPENDENCIES two_d_poisson
                      two_d_poisson_without_singularity_subtracted_off
  EXTRA_REQUIRES
  LABELS "poisson;poisson_with_singularity")

message(VERBOSE "Leaving poisson_with_singularity subdirectory")
# ------------------------------------------------------------------------------
