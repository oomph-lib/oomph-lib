# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered insoluble_surfactant subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(insoluble_surfactant C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME surfactant
  SOURCES surfactant.cc
  LIBRARIES oomph::advection_diffusion oomph::navier_stokes oomph::multi_physics
            oomph::fluid_interface
  CXX_STANDARD 17)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME multi_physics.insoluble_surfactant
  TARGET_DEPENDENCIES surfactant
  EXTRA_REQUIRES
  LABELS multi_physics insoluble_surfactant)

message(VERBOSE "Leaving insoluble_surfactant subdirectory")
# ------------------------------------------------------------------------------
