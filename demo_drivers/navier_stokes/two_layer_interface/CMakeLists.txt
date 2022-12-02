# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered two_layer_interface subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(two_layer_interface C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME elastic_two_layer_interface
  SOURCES elastic_two_layer_interface.cc
  LIBRARIES oomph::navier_stokes oomph::fluid_interface oomph::solid
            oomph::constitutive oomph::meshes oomph::generic)

oomph_add_executable(
  NAME spine_two_layer_interface
  SOURCES spine_two_layer_interface.cc
  LIBRARIES oomph::navier_stokes oomph::fluid_interface oomph::meshes
            oomph::generic)

oomph_add_executable(
  NAME refineable_two_layer_interface
  SOURCES refineable_two_layer_interface.cc
  LIBRARIES oomph::navier_stokes oomph::fluid_interface oomph::solid
            oomph::constitutive oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.two_layer_interface
  TARGET_DEPENDENCIES elastic_two_layer_interface spine_two_layer_interface
                      refineable_two_layer_interface
  EXTRA_REQUIRES
  LABELS navier_stokes two_layer_interface)

message(VERBOSE "Leaving two_layer_interface subdirectory")
# ------------------------------------------------------------------------------
