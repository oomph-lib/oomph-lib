# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered fsi_osc_ring subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(fsi_osc_ring C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME fsi_osc_ring
  SOURCES fsi_osc_ring.cc
  LIBRARIES oomph::beam oomph::navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME fsi_osc_ring_compare_jacs
  SOURCES fsi_osc_ring_compare_jacs.cc
  LIBRARIES oomph::beam oomph::navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME osc_ring_alg
  SOURCES osc_ring_alg.cc osc_ring_sarah_asymptotics.h
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME osc_ring_macro
  SOURCES osc_ring_macro.cc osc_ring_sarah_asymptotics.h
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME interaction.fsi_osc_ring
  TARGET_DEPENDENCIES fsi_osc_ring fsi_osc_ring_compare_jacs osc_ring_alg
                      osc_ring_macro
  EXTRA_REQUIRES
  LABELS interaction fsi_osc_ring)

message(VERBOSE "Leaving fsi_osc_ring subdirectory")
# ------------------------------------------------------------------------------
