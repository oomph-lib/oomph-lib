# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered rayleigh_instability subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(rayleigh_instability C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME rayleigh_instability
  SOURCES rayleigh_instability.cc
  LIBRARIES oomph::axisym_navier_stokes oomph::fluid_interface oomph::meshes
            oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME axisym_navier_stokes.rayleigh_instability
  TARGET_DEPENDENCIES rayleigh_instability
  EXTRA_REQUIRES
  LABELS axisym_navier_stokes rayleigh_instability)

message(VERBOSE "Leaving rayleigh_instability subdirectory")
# ------------------------------------------------------------------------------
