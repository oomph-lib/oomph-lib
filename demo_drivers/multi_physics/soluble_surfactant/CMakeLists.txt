# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered soluble_surfactant subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(soluble_surfactant C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME soluble_surfactant
  SOURCES double_buoyant_navier_stokes_elements.h soluble_surfactant.cc
  LIBRARIES oomph::advection_diffusion_reaction oomph::navier_stokes
            oomph::fluid_interface)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME multi_physics.soluble_surfactant
  TARGET_DEPENDENCIES soluble_surfactant
  EXTRA_REQUIRES
  LABELS multi_physics soluble_surfactant)

message(VERBOSE "Leaving soluble_surfactant subdirectory")
# ------------------------------------------------------------------------------
