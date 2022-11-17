# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered thermal_fibre subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(thermal_fibre C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME thermal_fibre
  SOURCES non_isothermal_axisym_nst_elements.h thermal_fibre.cc
  LIBRARIES oomph::axisym_navier_stokes oomph::fluid_interface
            oomph::steady_axisym_advection_diffusion
  CXX_STANDARD 17)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME multi_physics.thermal_fibre
  TARGET_DEPENDENCIES thermal_fibre
  EXTRA_REQUIRES
  LABELS multi_physics thermal_fibre)

message(VERBOSE "Leaving thermal_fibre subdirectory")
# ------------------------------------------------------------------------------
