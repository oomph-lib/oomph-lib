# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered single_layer_free_surface subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(single_layer_free_surface C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME elastic_single_layer
  SOURCES elastic_single_layer.cc
  LIBRARIES oomph::navier_stokes oomph::fluid_interface oomph::solid
            oomph::constitutive
  CXX_STANDARD 17)

oomph_add_executable(
  NAME spine_single_layer
  SOURCES spine_single_layer.cc
  LIBRARIES oomph::navier_stokes oomph::fluid_interface
  CXX_STANDARD 17)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.single_layer_free_surface
  TARGET_DEPENDENCIES elastic_single_layer spine_single_layer
  EXTRA_REQUIRES
  LABELS navier_stokes single_layer_free_surface)

message(VERBOSE "Leaving single_layer_free_surface subdirectory")
# ------------------------------------------------------------------------------
