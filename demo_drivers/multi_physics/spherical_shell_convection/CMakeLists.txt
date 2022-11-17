# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered spherical_shell_convection subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(spherical_shell_convection C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME spherical_shell_convection
  SOURCES spherical_shell_convection.cc spherical_buoyant_navier_stokes.h
  LIBRARIES oomph::spherical_advection_diffusion oomph::spherical_navier_stokes
  CXX_STANDARD 17)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME multi_physics.spherical_shell_convection
  TARGET_DEPENDENCIES spherical_shell_convection
  EXTRA_REQUIRES
  LABELS multi_physics spherical_shell_convection)

message(VERBOSE "Leaving spherical_shell_convection subdirectory")
# ------------------------------------------------------------------------------
