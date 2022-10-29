# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered two_layer_soluble_surfactant subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(two_layer_soluble_surfactant C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME refineable_two_layer_soluble_surfactant
  SOURCES double_buoyant_navier_stokes_elements.h
          refineable_two_layer_soluble_surfactant.cc
  LIBRARIES oomph::advection_diffusion_reaction oomph::navier_stokes
            oomph::fluid_interface oomph::solid oomph::constitutive
  CXX_STANDARD 14)

oomph_add_executable(
  NAME two_layer_soluble_surfactant
  SOURCES soluble_surfactant_transport_equations.h
          double_buoyant_navier_stokes_elements.h
          two_layer_soluble_surfactant.cc
          soluble_surfactant_transport_equations.h
  LIBRARIES oomph::advection_diffusion_reaction oomph::navier_stokes
            oomph::fluid_interface
  CXX_STANDARD 14)

oomph_add_executable(
  NAME refineable_unstructured_two_layer_ss
  SOURCES double_buoyant_navier_stokes_elements.h
          refineable_unstructured_two_layer_ss.cc
  LIBRARIES oomph::advection_diffusion_reaction oomph::navier_stokes
            oomph::fluid_interface oomph::solid oomph::constitutive
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME multi_physics.two_layer_soluble_surfactant
  TARGET_DEPENDENCIES
    refineable_two_layer_soluble_surfactant two_layer_soluble_surfactant
    refineable_unstructured_two_layer_ss
  EXTRA_REQUIRES
  LABELS multi_physics two_layer_soluble_surfactant)

message(VERBOSE "Leaving two_layer_soluble_surfactant subdirectory")
# ------------------------------------------------------------------------------
