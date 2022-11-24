# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered rayleigh_instability_surfactant subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(rayleigh_instability_surfactant C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME 3d_rayleigh_instability_surfactant
  SOURCES 3d_rayleigh_instability_surfactant.cc
  LIBRARIES oomph::navier_stokes oomph::fluid_interface)

oomph_add_executable(
  NAME rayleigh_instability_insoluble_surfactant
  SOURCES rayleigh_instability_insoluble_surfactant.cc
          axisymmetric_advection_navier_stokes_elements.h
  LIBRARIES oomph::axisym_navier_stokes oomph::fluid_interface)

oomph_add_executable(
  NAME rayleigh_instability_soluble_surfactant
  SOURCES rayleigh_instability_soluble_surfactant.cc
  LIBRARIES oomph::axisym_navier_stokes oomph::fluid_interface oomph::solid
            oomph::constitutive oomph::axisym_advection_diffusion)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME multi_physics.rayleigh_instability_surfactant
  TARGET_DEPENDENCIES
    3d_rayleigh_instability_surfactant rayleigh_instability_insoluble_surfactant
    rayleigh_instability_soluble_surfactant
  EXTRA_REQUIRES
  LABELS multi_physics rayleigh_instability_surfactant)

message(VERBOSE "Leaving rayleigh_instability_surfactant subdirectory")
# ------------------------------------------------------------------------------
