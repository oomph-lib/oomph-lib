# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered steady_ring subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(steady_ring C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME steady_ring
  SOURCES steady_ring.cc
  LIBRARIES oomph::beam
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME beam.steady_ring
  TARGET_DEPENDENCIES steady_ring
  EXTRA_REQUIRES
  LABELS "beam;steady_ring")

message(VERBOSE "Leaving steady_ring subdirectory")
# ------------------------------------------------------------------------------
