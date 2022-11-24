# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered unsteady_ring subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(unsteady_ring C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME lin_unsteady_ring
  SOURCES lin_unsteady_ring.cc
  LIBRARIES oomph::beam)

oomph_add_executable(
  NAME unsteady_ring
  SOURCES unsteady_ring.cc
  LIBRARIES oomph::beam)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME beam.unsteady_ring
  TARGET_DEPENDENCIES lin_unsteady_ring unsteady_ring
  EXTRA_REQUIRES
  LABELS beam unsteady_ring)

message(VERBOSE "Leaving unsteady_ring subdirectory")
# ------------------------------------------------------------------------------
