# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered osc_ring subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(osc_ring C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME osc_ring_alg
  SOURCES osc_ring_alg.cc
  LIBRARIES oomph::navier_stokes)

oomph_add_executable(
  NAME osc_ring_macro
  SOURCES osc_ring_macro.cc
  LIBRARIES oomph::navier_stokes)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.osc_ring
  TARGET_DEPENDENCIES osc_ring_alg osc_ring_macro
  EXTRA_REQUIRES
  LABELS navier_stokes osc_ring)

message(VERBOSE "Leaving osc_ring subdirectory")
# ------------------------------------------------------------------------------
