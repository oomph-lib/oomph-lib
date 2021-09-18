# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered two_layer_interface_axisym subdirectory")

cmake_minimum_required(VERSION 3.16..3.20 FATAL_ERROR)
project(two_layer_interface_axisym C CXX Fortran)
if(NOT BUILDING_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME elastic_two_layer_interface_axisym
  SOURCES elastic_two_layer_interface_axisym.cc
  LIBRARIES oomph::axisym_navier_stokes oomph::navier_stokes
            oomph::fluid_interface oomph::solid oomph::constitutive
  CXX_STANDARD 14)

oomph_add_executable(
  NAME spine_two_layer_interface_axisym
  SOURCES spine_two_layer_interface_axisym.cc
  LIBRARIES oomph::axisym_navier_stokes oomph::navier_stokes
            oomph::fluid_interface
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME axisym_navier_stokes.two_layer_interface_axisym
  TARGET_DEPENDENCIES elastic_two_layer_interface_axisym
                      spine_two_layer_interface_axisym
  EXTRA_REQUIRES
  LABELS "axisym_navier_stokes;two_layer_interface_axisym")

message(VERBOSE "Leaving two_layer_interface_axisym subdirectory")
# ------------------------------------------------------------------------------
